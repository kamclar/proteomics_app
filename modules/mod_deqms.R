# Differential abundance with limma and DEqMS



mod_deqms_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Step 4. Differential Expression (DEqMS / limma)"),
    uiOutput(ns("guard")),
    uiOutput(ns("main_ui"))
  )
}

mod_deqms_server <- function(id, app_state) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    deqms_results <- reactiveVal(list())   # named list: pair -> result df

    output$guard <- renderUI({
      if (!isTRUE(app_state$imputation_done)) {
        tags$div(class = "alert alert-warning",
                 " Please complete Step 2 (Imputation) first.")
      }
    })
    # Comparisons are derived from the cleaned sample names after import.
    available_pairs <- reactive({
      req(isTRUE(app_state$imputation_done))
      snames <- colnames(app_state$imputed_data$intensity)
      all_pairs(infer_groups(snames))
    })

    ttest_pairs <- reactive({
      if (!isTRUE(app_state$upload_done)) return(character(0))
      tt <- app_state$parsed_data$ttest
      if (is.null(tt) || !length(tt)) return(character(0))
      names(tt)
    })

    reverse_pair_name <- function(pair) {
      parts <- strsplit(pair, "_vs_", fixed = TRUE)[[1]]
      if (length(parts) != 2) return(NA_character_)
      paste(parts[2], parts[1], sep = "_vs_")
    }

    ttest_matched_pairs <- reactive({
      pairs <- available_pairs()
      tt <- ttest_pairs()
      if (!length(pairs) || !length(tt)) return(character(0))

      rev_pairs <- vapply(pairs, reverse_pair_name, character(1))
      pairs[pairs %in% tt | rev_pairs %in% tt]
    })

    default_selected_pairs <- reactive({
      matched <- ttest_matched_pairs()
      if (length(ttest_pairs())) matched else available_pairs()
    })

    output$ttest_pair_note <- renderUI({
      pairs <- available_pairs()
      tt <- ttest_pairs()
      matched <- ttest_matched_pairs()

      if (!length(tt)) {
        return(tags$div(
          class = "alert alert-info",
          "No uploaded t-test comparisons were detected, so all DEqMS comparisons are selected by default."
        ))
      }

      unmatched_ttest <- setdiff(tt, c(matched, vapply(matched, reverse_pair_name, character(1))))

      tags$div(
        class = "alert alert-info",
        tags$b("Preselected from uploaded t-test results: "),
        length(matched), " of ", length(pairs), " available DEqMS comparison(s). ",
        "A comparison is checked when the same pair, or its reverse direction, was found in the uploaded t-test columns.",
        if (length(unmatched_ttest)) {
          tagList(
            tags$br(),
            tags$small(
              "T-test comparison(s) not available from current sample groups: ",
              paste(unmatched_ttest, collapse = ", ")
            )
          )
        }
      )
    })

    output$main_ui <- renderUI({
      req(isTRUE(app_state$imputation_done))
      pairs <- available_pairs()
      tagList(
        fluidRow(
          column(4,
            wellPanel(
              h4("Select Comparisons"),
              helpText("Select pairs to compute. Large numbers of pairs may take time."),
              uiOutput(ns("ttest_pair_note")),
              br(),
              div(style = "display:flex; gap:8px; margin-bottom:6px;",
                actionButton(ns("select_all"),   "Select All",   class = "btn-xs btn-default"),
                actionButton(ns("deselect_all"), "Deselect All", class = "btn-xs btn-default")
              ),
              checkboxGroupInput(ns("pairs"), label = NULL,
                                 choices  = pairs,
                                 selected = default_selected_pairs()),
              radioButtons(
                ns("run_mode"), "Run mode",
                choices = c(
                  "Replace all existing DEqMS results" = "replace",
                  "Keep existing results and update selected pairs" = "update"
                ),
                selected = "replace"
              ),
              hr(),
              h4("Thresholds"),
              numericInput(ns("adj_pval"), "Adj. p-value cutoff",
                           value = 0.05, min = 0.001, max = 0.5, step = 0.005),
              numericInput(ns("log2fc"), "log2FC threshold",
                           value = 0.5, min = 0, max = 5, step = 0.1),
              hr(),
              actionButton(ns("btn_run"), "Run DEqMS",
                           class = "btn-primary btn-block"),
              actionButton(ns("btn_clear"), "Clear DEqMS Results",
                           class = "btn-default btn-block"),
              br(),
              uiOutput(ns("run_status"))
            )
          ),
          column(8,
            tabsetPanel(
              tabPanel("Results Table",
                br(),
                uiOutput(ns("result_selector")),
                br(),
                DT::dataTableOutput(ns("result_table"))
              ),
              tabPanel("Volcano Preview",
                br(),
                uiOutput(ns("volcano_selector_ui")),
                br(),
                plotOutput(ns("volcano_preview"), height = "450px")
              ),
              tabPanel("TKO Variance Plot",
                br(),
                uiOutput(ns("tko_selector_ui")),
                br(),
                plotOutput(ns("tko_plot"), height = "420px")
              ),
              tabPanel("Summary",
                br(),
                DT::dataTableOutput(ns("summary_table"))
              )
            )
          )
        ),
        fluidRow(
          column(12,
            div(class = "app-action-row",
              downloadButton(ns("dl_all"), "Download All Results (ZIP)",
                             class = "btn-default"),
              div(class = "app-proceed-wrap",
                uiOutput(ns("proceed_btn"))
              )
            )
          )
        )
      )
    })

    # Select/deselect all
    observeEvent(input$select_all,   {
      updateCheckboxGroupInput(session, "pairs", selected = available_pairs())
    })
    observeEvent(input$deselect_all, {
      updateCheckboxGroupInput(session, "pairs", selected = character(0))
    })

    clear_downstream_deqms_state <- function() {
      app_state$deqms_results <- list()
      app_state$tko_fits <- list()
      app_state$deqms_done <- FALSE
      app_state$enrich_results <- list()
      app_state$enrichment_done <- FALSE
    }

    observeEvent(input$btn_clear, {
      deqms_results(list())
      clear_downstream_deqms_state()
      showNotification("DEqMS results cleared.", type = "message", duration = 4)
    })
    observeEvent(input$btn_run, {
      req(isTRUE(app_state$imputation_done))
      req(length(input$pairs) > 0)

      d          <- app_state$imputed_data
      int        <- as.matrix(d$intensity)
      meta       <- d$meta
      snames     <- colnames(int)
      grps_vec   <- infer_groups(snames)
      groups_map <- split(snames, grps_vec)


      # Peptide counts for DEqMS (use meta$Peptides if available)
      if (!"Peptides" %in% names(meta)) {
        showNotification("DEqMS requires a valid peptide count column. No peptide count column was found in metadata.",
                         type = "error", duration = 12)
        return()
      }
      counts <- suppressWarnings(as.numeric(meta$Peptides))
      bad_counts <- is.na(counts) | counts < 1
      if (any(bad_counts)) {
        showNotification(
          sprintf("DEqMS requires valid peptide counts for every protein. %d protein(s) have missing or invalid peptide counts.",
                  sum(bad_counts)),
          type = "error", duration = 12
        )
        return()
      }
      peptide_counts <- setNames(counts, rownames(meta))

      pairs_to_run <- input$pairs
      n_pairs      <- length(pairs_to_run)
      results_new  <- list()
      tko_new      <- list()

      withProgress(message = "Running DEqMS...", value = 0, {
        for (i in seq_along(pairs_to_run)) {
          pair    <- pairs_to_run[[i]]
          incProgress(1 / n_pairs,
                      detail = paste0("Pair ", i, "/", n_pairs, ": ", pair))

          parts <- strsplit(pair, "_vs_")[[1]]
          g1 <- parts[1]; g2 <- parts[2]

          cols1 <- groups_map[[g1]]
          cols2 <- groups_map[[g2]]
          if (is.null(cols1) || is.null(cols2)) {
            message("Skipping ", pair, " - missing columns")
            next
          }

          res_fwd <- tryCatch(
            run_deqms_pair(int, cols1, cols2, g1, g2, peptide_counts, meta),
            error = function(e) {
              showNotification(paste("DEqMS error on", pair, ":", conditionMessage(e)),
                               type = "warning", duration = 8)
              NULL
            }
          )

          if (!is.null(res_fwd)) {
            # Store TKO fit object for variance plot (forward pair only)
            fit4 <- attr(res_fwd, "fit4")
            if (!is.null(fit4)) {
              tko_new[[pair]] <- fit4
            }

            # Annotate
            res_fwd <- annotate_deqms(res_fwd, g1, g2, meta,
                                      input$adj_pval, input$log2fc)
            results_new[[pair]] <- res_fwd

            # Mirror (flip sign, no re-fit)
            pair_rev  <- paste0(g2, "_vs_", g1)
            res_rev   <- mirror_deqms(res_fwd, g1, g2,
                                       input$adj_pval, input$log2fc)
            results_new[[pair_rev]] <- res_rev
          }
        }
      })

      if (identical(input$run_mode, "update")) {
        existing <- deqms_results()
        existing[names(results_new)] <- results_new
        deqms_results(existing)

        existing_tko <- app_state$tko_fits
        existing_tko[names(tko_new)] <- tko_new
        app_state$tko_fits <- existing_tko
      } else {
        deqms_results(results_new)
        app_state$tko_fits <- tko_new
      }

      app_state$deqms_results  <- deqms_results()
      app_state$deqms_done     <- length(deqms_results()) > 0
      app_state$enrich_results <- list()
      app_state$enrichment_done <- FALSE

      showNotification(
        sprintf("DEqMS run complete: %d result tables available.", length(deqms_results())),
        type = "message",
        duration = 5
      )
    })

    output$run_status <- renderUI({
      res <- deqms_results()
      if (!length(res)) return(NULL)
      tags$div(class = "alert alert-success",
               sprintf(" %d comparisons computed", length(res)))
    })
    output$result_selector <- renderUI({
      res <- deqms_results()
      if (!length(res)) return(tags$p("No results yet."))
      selectInput(ns("selected_result"), "View comparison:",
                  choices = names(res), selected = names(res)[1])
    })

    output$volcano_selector_ui <- renderUI({
      res <- deqms_results()
      if (!length(res)) return(NULL)
      selectInput(ns("volcano_pair"), "Comparison:",
                  choices = names(res), selected = names(res)[1])
    })

    output$tko_selector_ui <- renderUI({
      res <- deqms_results()
      if (!length(res)) return(NULL)
      # TKO stored per forward pair
      fwd_pairs <- names(res)[!grepl("^.*_vs_.*_vs_", names(res))]
      selectInput(ns("tko_pair"), "Comparison:",
                  choices = names(res), selected = names(res)[1])
    })
    output$result_table <- DT::renderDataTable({
      req(length(deqms_results()) > 0, input$selected_result)
      df <- deqms_results()[[input$selected_result]]
      req(!is.null(df))
      show_cols <- intersect(
        c("Genes", "Protein.Names", "logFC", "sca.adj.pval",
          "negLog10_sca_adjPval", "Diff_Abund", "AveExpr", "t", "P.Value"),
        names(df)
      )
      DT::datatable(df[, show_cols, drop = FALSE], rownames = TRUE,
                    options = list(pageLength = 15, scrollX = TRUE, dom = "ftip")) |>
        DT::formatRound(intersect(c("logFC", "AveExpr", "t", "P.Value",
                                     "sca.adj.pval", "negLog10_sca_adjPval"), names(df)), 4) |>
        DT::formatStyle("Diff_Abund",
          backgroundColor = DT::styleEqual(
            c("Non-significant"),
            c("rgba(200,200,200,0.2)")
          )
        )
    })
    # Small static volcano preview; the interactive plotter has its own screen.
    output$volcano_preview <- renderPlot({
      req(length(deqms_results()) > 0, input$volcano_pair)
      df     <- deqms_results()[[input$volcano_pair]]
      req(!is.null(df), "logFC" %in% names(df), "negLog10_sca_adjPval" %in% names(df))

      fc_thr <- input$log2fc
      p_thr  <- -log10(input$adj_pval)

      col_vec <- ifelse(df$negLog10_sca_adjPval >= p_thr & df$logFC >= fc_thr, "#F76075",
                 ifelse(df$negLog10_sca_adjPval >= p_thr & df$logFC <= -fc_thr, "#40BFC1",
                        "#BDBFBE"))

      par(mar = c(5, 5, 3, 2))
      plot(df$logFC, df$negLog10_sca_adjPval,
           col = col_vec, pch = 19, cex = 0.8,
           xlab = "log2 Fold Change",
           ylab = "-log10(adj. p-value)",
           main = paste("DEqMS:", input$volcano_pair))
      abline(h = p_thr, lty = 2, col = "gray40")
      abline(v = c(-fc_thr, fc_thr), lty = 2, col = "gray40")

      # Label top hits
      sig   <- df[df$negLog10_sca_adjPval >= p_thr & abs(df$logFC) >= fc_thr, ]
      top_n <- head(sig[order(sig$negLog10_sca_adjPval, decreasing = TRUE), ], 10)
      if (nrow(top_n) > 0 && "Genes" %in% names(top_n)) {
        text(top_n$logFC, top_n$negLog10_sca_adjPval,
             labels = top_n$Genes, cex = 0.65, pos = 3, col = "black")
      }

      # Legend
      n_up   <- sum(col_vec == "#F76075",  na.rm = TRUE)
      n_down <- sum(col_vec == "#40BFC1",  na.rm = TRUE)
      n_ns   <- sum(col_vec == "#BDBFBE",  na.rm = TRUE)
      legend("topright",
             legend = c(sprintf("Up (%d)", n_up),
                        sprintf("Down (%d)", n_down),
                        sprintf("NS (%d)", n_ns)),
             col = c("#F76075", "#40BFC1", "#BDBFBE"),
             pch = 19, bty = "n", cex = 0.85)
    })
    # TKO variance plot
    output$tko_plot <- renderPlot({
      req(length(deqms_results()) > 0, input$tko_pair)
      # TKO plot is stored as attribute in the fit - re-compute here just the plot portion
      # We store tko_fits in app_state
      fit4 <- app_state$tko_fits[[input$tko_pair]]
      if (is.null(fit4)) {
        plot.new()
        text(0.5, 0.5, "TKO plot not available for this comparison.\n(May be a mirror pair)", cex = 1)
        return()
      }
      DEqMS::VarianceBoxplot(fit4, n = 30,
                              main = paste("TKO Variance:", input$tko_pair),
                              xlab = "Peptide count")
    })
    output$summary_table <- DT::renderDataTable({
      res <- deqms_results()
      if (!length(res)) return(data.frame())
      rows <- lapply(names(res), function(nm) {
        df   <- res[[nm]]
        fc_t <- input$log2fc %||% 0.5
        p_t  <- input$adj_pval %||% 0.05
        up   <- sum(df$logFC >= fc_t & df$sca.adj.pval <= p_t, na.rm = TRUE)
        down <- sum(df$logFC <= -fc_t & df$sca.adj.pval <= p_t, na.rm = TRUE)
        data.frame(Comparison = nm, Up = up, Down = down,
                   Total_Significant = up + down, stringsAsFactors = FALSE)
      })
      df_sum <- do.call(rbind, rows)
      DT::datatable(df_sum, rownames = FALSE,
                    options = list(pageLength = 20, dom = "ftip"))
    })
    output$dl_all <- downloadHandler(
      filename = function() download_filename(app_state, "DEqMS_results", "zip"),
      content  = function(file) {
        res <- deqms_results()
        tmp <- tempdir()
        fnames <- vapply(names(res), function(nm) {
          f <- file.path(tmp, paste0("res_DEqMS_", nm, ".csv"))
          utils::write.csv(res[[nm]], f, row.names = TRUE)
          f
        }, character(1))
        zip::zip(zipfile = file, files = fnames, mode = "cherry-pick")
      }
    )
    output$proceed_btn <- renderUI({
      if (!isTRUE(app_state$deqms_done)) return(NULL)
      actionButton(ns("btn_proceed"), "Proceed to Volcano Plotter",
                   class = "btn-success btn-lg app-proceed-btn")
    })

    observeEvent(input$btn_proceed, {
      app_state$active_tab <- "volcano"
    })
  })
}
# DEqMS helpers
run_deqms_pair <- function(int_mat, cols1, cols2, g1, g2, peptide_counts, meta) {
  dfD   <- int_mat[, c(cols1, cols2), drop = FALSE]
  
  cond  <- factor(c(rep(g1, length(cols1)), rep(g2, length(cols2))),
                  levels = c(g1, g2))
  
  design   <- model.matrix(~0 + cond)
  colnames(design) <- gsub("^cond", "", colnames(design))
  
  contrast <- limma::makeContrasts(
    contrasts = paste0(g1, "-", g2), levels = design
  )
  
  fit1 <- limma::lmFit(dfD, design)
  
  fit2 <- limma::contrasts.fit(fit1, contrasts = contrast)
  
  fit3 <- limma::eBayes(fit2)
  
  rn <- rownames(fit3$coefficients)
  fit3$count <- peptide_counts[rn]
  if (any(is.na(fit3$count) | fit3$count < 1)) {
    stop("Missing or invalid peptide counts for DEqMS fit.")
  }
  
  fit4 <- DEqMS::spectraCounteBayes(fit3)
  
  res <- DEqMS::outputResult(fit4, coef_col = 1)
  
  attr(res, "fit4") <- fit4
  res
}

annotate_deqms <- function(res, g1, g2, meta, adj_p_cut = 0.05, fc_cut = 0.5) {
  pair_name <- paste0(g1, "_vs_", g2)

  res$Diff_Abund <- ifelse(
    res$sca.adj.pval <= adj_p_cut & res$logFC >=  fc_cut, paste("Up in", g1),
    ifelse(res$sca.adj.pval <= adj_p_cut & res$logFC <= -fc_cut, paste("Up in", g2),
           "Non-significant")
  )
  res$negLog10_sca_adjPval <- -log10(res$sca.adj.pval)

  for (col in c("Genes", "First.Protein.Descriptions", "Protein.Ids")) {
    if (col %in% names(meta)) {
      res[[col]] <- meta[rownames(res), col]
    }
  }
  res$Protein.Names <- rownames(res)
  res[order(res$sca.adj.pval), , drop = FALSE]
}

mirror_deqms <- function(res_fwd, g1, g2, adj_p_cut = 0.05, fc_cut = 0.5) {
  res_rev        <- res_fwd
  res_rev$logFC  <- -res_fwd$logFC
  if ("t" %in% names(res_rev)) res_rev$t <- -res_fwd$t

  res_rev$Diff_Abund <- ifelse(
    res_rev$sca.adj.pval <= adj_p_cut & res_rev$logFC >=  fc_cut, paste("Up in", g1),
    ifelse(res_rev$sca.adj.pval <= adj_p_cut & res_rev$logFC <= -fc_cut, paste("Up in", g2),
           "Non-significant")
  )
  res_rev
}

`%||%` <- function(a, b) if (!is.null(a)) a else b
