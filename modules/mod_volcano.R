# Interactive volcano plotter 


mod_volcano_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    h3("Step 5. Volcano Plotter"),
    fluidRow(
      column(
        3,
        wellPanel(
          h4("Data Source"),
          radioButtons(
            ns("source"), NULL,
            choices = c(
              "DEqMS results" = "deqms",
              "T-test (from uploaded file)" = "ttest"
            ),
            selected = "deqms"
          ),
          
          uiOutput(ns("comparison_ui")),
          hr(),
          
          h4("Thresholds"),
          numericInput(
            ns("fc_thr"), "log2 FC threshold:",
            value = 0.5, min = 0, max = 5, step = 0.1
          ),
          numericInput(
            ns("p_thr"), "-log10(p) threshold:",
            value = 1.3, min = 0, max = 10, step = 0.1
          ),
          hr(),
          
          h4("Gene Search"),
          textInput(
            ns("gene_search"),
            "Find gene:",
            placeholder = "Enter gene name..."
          ),
            actionButton(
              ns("btn_search"),
              "Label Gene",
              class = "btn-sm btn-primary"
            ),
          hr(),
          
          h4("Display Options"),
          checkboxInput(
            ns("show_ns"),
            "Show non-significant points",
            value = TRUE
          ),
          sliderInput(
            ns("ns_sample"),
            "NS point sampling %:",
            min = 5, max = 100, value = 20, step = 5
          ),
          actionButton(
            ns("btn_flip"),
            "Flip Sides"
          ),
          actionButton(
            ns("btn_clear"),
            "Clear Labels"
          ),
          hr(),

          h4("Data Flags"),
          checkboxInput(
            ns("show_low_abundant"),
            "Low-intensity partial missing: grey outline",
            value = FALSE
          ),
          checkboxInput(
            ns("show_on_off"),
            "On/off: triangle",
            value = FALSE
          ),
          tags$p(
            class = "help-block",
            "Low-intensity partial missing = lowest 10% raw mean intensity in the selected comparison, with some but not all values missing.",
            tags$br(),
            "On/off = measured in only one group of the selected comparison."
          ),
          uiOutput(ns("flag_counts")),
          actionButton(ns("btn_reactome_cache"), "Build Reactome cache", class = "btn-sm btn-default btn-block"),
          helpText("Adds currently unseen proteins to the local Reactome cache for faster reuse in Volcano and Enrichment."),
          uiOutput(ns("reactome_cache_status")),
          hr(),
          
          h4("Axis Ranges"),
          fluidRow(
            column(6, numericInput(ns("x_min"), "X min:", value = NA)),
            column(6, numericInput(ns("x_max"), "X max:", value = NA))
          ),
          fluidRow(
            column(6, numericInput(ns("y_min"), "Y min:", value = NA)),
            column(6, numericInput(ns("y_max"), "Y max:", value = NA))
          ),
          hr(),
          
          downloadButton(ns("dl_plot"), "Export SVG", class = "btn-sm btn-block"),
          downloadButton(ns("dl_csv"), "Export CSV", class = "btn-sm btn-block"),
          br(),
          uiOutput(ns("proceed_btn"))
        )
      ),
      column(
        9,
        fluidRow(
          column(
            12,
            tags$div(
              class = "alert alert-info",
              style = "margin-bottom: 12px;",
              "Tip: click points in the volcano plot to highlight matching rows in the table, or select rows in the table to label and highlight the same proteins in the plot."
            )
          ),
          column(
            8,
            plotly::plotlyOutput(ns("volcano_plot"), height = "700px")
          ),
          column(
            4,
            h4("Selected Comparison Table"),
            DT::dataTableOutput(ns("volcano_table"))
          ),
          column(
            12,
            conditionalPanel(
              condition = paste0("input['", ns("show_low_abundant"), "'] || input['", ns("show_on_off"), "']"),
              hr(),
              h4("Low-intensity partial missing / on-off"),
              DT::dataTableOutput(ns("flagged_table"))
            )
          )
        )
      )
    )
  )
}

mod_volcano_server <- function(id, app_state) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    flip_state <- reactiveVal(FALSE)
    clicked_points <- reactiveVal(integer(0))   # stores point_id
    volcano_data <- reactiveVal(NULL)
    reactome_cache <- reactiveVal(read_reactome_cache())

    comparison_groups <- function(comp) {
      if (is.null(comp) || !nzchar(comp)) return(NULL)
      parts <- strsplit(comp, "_vs_", fixed = TRUE)[[1]]
      if (length(parts) != 2) return(NULL)
      parts
    }

    volcano_missing_flags <- function(df, comp) {
      if (!isTRUE(app_state$upload_done)) return(df)
      int <- app_state$parsed_data$intensity
      sn  <- app_state$parsed_data$sample_names
      if (is.null(int) || is.null(sn) || !nrow(int)) return(df)

      parts <- comparison_groups(comp)
      if (is.null(parts)) return(df)

      grps <- infer_groups(sn)
      cols1 <- intersect(sn[grps == parts[1]], colnames(int))
      cols2 <- intersect(sn[grps == parts[2]], colnames(int))
      if (!length(cols1) || !length(cols2)) return(df)

      int_mat <- as.matrix(int)
      comp_cols <- unique(c(cols1, cols2))
      comp_mat <- int[, comp_cols, drop = FALSE]
      raw_mean <- rowMeans(comp_mat, na.rm = TRUE)
      n_valid_comp <- rowSums(!is.na(comp_mat))
      n_missing_comp <- rowSums(is.na(comp_mat))
      raw_mean[n_valid_comp == 0] <- NA_real_

      n1 <- rowSums(!is.na(int[, cols1, drop = FALSE]))
      n2 <- rowSums(!is.na(int[, cols2, drop = FALSE]))
      mean1 <- rowMeans(int[, cols1, drop = FALSE], na.rm = TRUE)
      mean2 <- rowMeans(int[, cols2, drop = FALSE], na.rm = TRUE)
      mean1[n1 == 0] <- NA_real_
      mean2[n2 == 0] <- NA_real_
      left_floor <- stats::quantile(as.numeric(int_mat), probs = 0.01, na.rm = TRUE, names = FALSE)
      if (!is.finite(left_floor)) left_floor <- min(as.numeric(int_mat), na.rm = TRUE)
      on_off <- (n1 > 0 & n2 == 0) | (n1 == 0 & n2 > 0)
      on_off_direction <- ifelse(n1 > 0 & n2 == 0, paste0(parts[1], " only"),
                          ifelse(n1 == 0 & n2 > 0, paste0(parts[2], " only"), ""))
      on_off_display_x <- ifelse(n1 > 0 & n2 == 0, mean1 - left_floor,
                          ifelse(n1 == 0 & n2 > 0, left_floor - mean2, NA_real_))
      partial_missing <- n_missing_comp > 0 & n_valid_comp > 0 & !on_off
      low_cut <- stats::quantile(raw_mean[partial_missing], probs = 0.10, na.rm = TRUE, names = FALSE)
      low_abundant <- !is.na(raw_mean) & partial_missing & is.finite(low_cut) & raw_mean <= low_cut

      raw_flags <- data.frame(
        row_key = rownames(int),
        raw_mean = raw_mean,
        low_abundant = low_abundant,
        on_off = on_off,
        on_off_direction = on_off_direction,
        on_off_display_x = on_off_display_x,
        stringsAsFactors = FALSE
      )

      match_idx <- match(rownames(df), raw_flags$row_key)
      if (all(is.na(match_idx)) && nrow(df) == nrow(raw_flags)) {
        match_idx <- seq_len(nrow(df))
      }

      df$raw_mean <- raw_flags$raw_mean[match_idx]
      df$low_abundant <- raw_flags$low_abundant[match_idx]
      df$on_off <- raw_flags$on_off[match_idx]
      df$on_off_direction <- raw_flags$on_off_direction[match_idx]
      df$on_off_display_x <- raw_flags$on_off_display_x[match_idx]

      df$low_abundant[is.na(df$low_abundant)] <- FALSE
      df$on_off[is.na(df$on_off)] <- FALSE
      df$on_off_direction[is.na(df$on_off_direction)] <- ""
      df
    }

    build_volcano_df <- function() {
      src  <- input$source
      comp <- input$comparison
      flip <- flip_state()
      
      if (is.null(src) || is.null(comp)) return(NULL)
      # Use the same plotting controls for both upstream t-tests and DEqMS.
      if (src == "ttest") {
        tt <- app_state$parsed_data$ttest
        if (is.null(tt) || !comp %in% names(tt)) return(NULL)
        
        df <- tt[[comp]]
        if (nrow(df) == 0) return(NULL)
        
        df$x <- df$logFC
        df$y <- df$negLog10p
        df$source_label <- "T-test"
        
      } else {
        res <- app_state$deqms_results
        if (is.null(res) || !comp %in% names(res)) return(NULL)
        
        df <- res[[comp]]
        if (nrow(df) == 0) return(NULL)
        
        df$x <- df$logFC
        df$y <- df$negLog10_sca_adjPval
        df$source_label <- "DEqMS"
      }
      if ("Genes" %in% names(df) && any(!is.na(df$Genes) & trimws(df$Genes) != "")) {
        df$label <- as.character(df$Genes)
      } else if ("Protein.Names" %in% names(df) && any(!is.na(df$Protein.Names) & trimws(df$Protein.Names) != "")) {
        df$label <- as.character(df$Protein.Names)
      } else if (!is.null(rownames(df))) {
        df$label <- rownames(df)
      } else {
        df$label <- rep(NA_character_, nrow(df))
      }
      
      bad_label <- is.na(df$label) | trimws(df$label) == ""
      if (any(bad_label)) {
        df$label[bad_label] <- paste0("row_", which(bad_label))
      }
      
      # Stable row ids keep Plotly clicks and DT selections synchronized.
      df$point_id <- seq_len(nrow(df))
      # Classify with the current thresholds.
      fc_t <- input$fc_thr
      p_t  <- input$p_thr
      
      finite_stat <- is.finite(df$y) & is.finite(df$logFC)
      df$colour <- ifelse(
        finite_stat & df$y >= p_t & df$logFC >= fc_t, "Up",
        ifelse(finite_stat & df$y >= p_t & df$logFC <= -fc_t, "Down", "NS")
      )
      
      # Flip only display x
      df$x <- if (flip) -df$logFC else df$logFC
      df <- volcano_missing_flags(df, comp)
      
      df
    }

    build_volcano_table_df <- function() {
      df <- build_volcano_df()
      if (is.null(df) || !nrow(df)) return(NULL)

      protein_col <- if ("Protein.Names" %in% names(df)) {
        "Protein.Names"
      } else if ("ProtDesc" %in% names(df)) {
        "ProtDesc"
      } else if ("Protein.Accessions" %in% names(df)) {
        "Protein.Accessions"
      } else {
        NA_character_
      }

      out <- data.frame(
        point_id = df$point_id,
        Gene = if ("Genes" %in% names(df)) as.character(df$Genes) else rep(NA_character_, nrow(df)),
        Protein = if (!is.na(protein_col)) as.character(df[[protein_col]]) else as.character(df$label),
        logFC = ifelse(is.finite(df$logFC), round(df$logFC, 4), NA_real_),
        negLog10P = ifelse(is.finite(df$y), round(df$y, 4), NA_real_),
        Low_abundant = if ("low_abundant" %in% names(df)) df$low_abundant else FALSE,
        On_off = if ("on_off" %in% names(df)) df$on_off else FALSE,
        On_off_direction = if ("on_off_direction" %in% names(df)) df$on_off_direction else "",
        stringsAsFactors = FALSE
      )

      out$Gene[is.na(out$Gene) | trimws(out$Gene) == ""] <- NA_character_
      out$Protein[is.na(out$Protein) | trimws(out$Protein) == ""] <- df$label[is.na(out$Protein) | trimws(out$Protein) == ""]
      out <- out[order(out$negLog10P, decreasing = TRUE), , drop = FALSE]
      rownames(out) <- NULL
      out
    }

    build_flagged_table_df <- function() {
      df <- build_volcano_table_df()
      if (is.null(df) || !nrow(df)) return(NULL)

      show_low <- isTRUE(input$show_low_abundant)
      show_on <- isTRUE(input$show_on_off)
      flagged <- df[(show_low & df$Low_abundant) | (show_on & df$On_off), , drop = FALSE]

      if (!nrow(flagged)) {
        selected <- paste(c(if (show_low) "low-intensity partial-missing", if (show_on) "on/off"), collapse = " or ")
        return(data.frame(Message = paste0("No ", selected, " points in this comparison.")))
      }
      flagged$Flag <- mapply(function(low, on) {
        paste(c(if (show_low && low) "Low-intensity partial-missing", if (show_on && on) "On/off"), collapse = "; ")
      }, flagged$Low_abundant, flagged$On_off, USE.NAMES = FALSE)
      flagged$Reactome <- vapply(flagged$Gene, function(gene) {
        if (is.na(gene) || !nzchar(trimws(gene))) return("")
        reactome_membership_hints(gene, cache = reactome_cache())
      }, character(1))
      flagged[, c("Gene", "Protein", "Flag", "On_off_direction", "Reactome", "logFC", "negLog10P"), drop = FALSE]
    }
    # Comparison selector
    output$comparison_ui <- renderUI({
      src <- input$source
      if (is.null(src)) return(NULL)
      
      if (src == "ttest") {
        tt <- if (isTRUE(app_state$upload_done)) app_state$parsed_data$ttest else list()
        
        if (!length(tt)) {
          return(tags$div(
            class = "alert alert-info",
            "No t-test columns detected in upload."
          ))
        }
        
        selectInput(ns("comparison"), "Comparison:", choices = names(tt))
        
      } else {
        res <- if (isTRUE(app_state$deqms_done)) app_state$deqms_results else list()
        
        if (!length(res)) {
          return(tags$div(
            class = "alert alert-info",
            "No DEqMS results yet. Run Step 4."
          ))
        }
        
        selectInput(ns("comparison"), "Comparison:", choices = names(res))
      }
    })
    # Rebuild data when the selected result or thresholds change.
    observeEvent(list(input$comparison, input$source, input$fc_thr, input$p_thr), {
      clicked_points(integer(0))
      volcano_data(build_volcano_df())
    }, ignoreNULL = FALSE)
    
    observeEvent(input$btn_flip, {
      flip_state(!flip_state())
      clicked_points(integer(0))
      volcano_data(build_volcano_df())
    })
    
    observeEvent(input$btn_clear, {
      clicked_points(integer(0))
    })

    observeEvent(input$btn_reactome_cache, {
      req(isTRUE(app_state$upload_done))
      symbols <- dataset_gene_symbols(app_state$parsed_data)
      if (!length(symbols)) {
        showNotification("No gene symbols found for Reactome cache.", type = "warning", duration = 5)
        return()
      }
      withProgress(message = "Building local Reactome cache...", value = 0.5, {
        res <- tryCatch(build_reactome_cache_for_symbols(symbols), error = function(e) {
          showNotification(paste("Reactome cache failed:", conditionMessage(e)), type = "error", duration = 10)
          NULL
        })
      })
      if (!is.null(res)) {
        reactome_cache(res$cache)
        st <- reactome_cache_status(res$cache)
        showNotification(
          sprintf("Reactome cache ready: %d mapped gene(s), %d pathway(s).", st$mapped_genes, st$pathways),
          type = "message",
          duration = 6
        )
      }
    })

    output$reactome_cache_status <- renderUI({
      st <- reactome_cache_status(reactome_cache())
      tags$small(
        sprintf("Cache: %d mapped gene(s), %d pathway(s).", st$mapped_genes, st$pathways)
      )
    })

    output$flag_counts <- renderUI({
      df <- build_volcano_df()
      if (is.null(df) || !nrow(df)) return(NULL)
      tags$small(
        sprintf("Current comparison: %d low-intensity partial-missing, %d on/off.",
                sum(df$low_abundant, na.rm = TRUE),
                sum(df$on_off, na.rm = TRUE))
      )
    })

    output$volcano_table <- DT::renderDataTable({
      df <- build_volcano_table_df()
      req(!is.null(df))

      DT::datatable(
        df[, c("Gene", "Protein", "logFC", "negLog10P"), drop = FALSE],
        rownames = FALSE,
        selection = list(mode = "multiple", target = "row"),
        options = list(
          pageLength = 15,
          scrollY = 640,
          scrollX = TRUE,
          dom = "ftip"
        )
      ) |>
        DT::formatRound(c("logFC", "negLog10P"), 4)
    }, server = FALSE)

    output$flagged_table <- DT::renderDataTable({
      df <- withProgress(message = "Preparing flagged protein table...", value = 0, {
        incProgress(0.35, detail = "Filtering selected flags")
        out <- build_flagged_table_df()
        incProgress(0.45, detail = "Adding Reactome hints")
        incProgress(0.20, detail = "Rendering table")
        out
      })
      req(!is.null(df))

      DT::datatable(
        df,
        rownames = FALSE,
        options = list(
          pageLength = 10,
          scrollY = 280,
          scrollX = TRUE,
          dom = "ftip"
        )
      )
    }, server = FALSE)

    observe({
      df <- build_volcano_table_df()
      req(!is.null(df))

      selected_rows <- which(df$point_id %in% clicked_points())
      DT::selectRows(DT::dataTableProxy(ns("volcano_table")), selected_rows)
    })

    observeEvent(input$volcano_table_rows_selected, {
      df <- build_volcano_table_df()
      rows <- input$volcano_table_rows_selected
      if (is.null(df)) return()

      if (!length(rows)) {
        clicked_points(integer(0))
        return()
      }

      clicked_points(df$point_id[rows])
    }, ignoreNULL = FALSE)
    # Search gene
    observeEvent(input$btn_search, {
      req(input$gene_search)
      
      search_term <- trimws(input$gene_search)
      if (nchar(search_term) == 0) return()
      
      df <- volcano_data()
      req(!is.null(df))
      
      labels_clean <- trimws(as.character(df$label))
      labels_lower <- tolower(labels_clean)
      term_lower   <- tolower(search_term)
      
      matches_idx <- integer(0)
      
      # 1-character query: exact match only
      if (nchar(search_term) == 1) {
        matches_idx <- which(
          !is.na(labels_clean) &
            labels_lower == term_lower
        )
        
      } else {
        # 2+ characters: exact -> starts-with -> contains
        matches_idx <- which(
          !is.na(labels_clean) &
            labels_lower == term_lower
        )
        
        if (length(matches_idx) == 0) {
          matches_idx <- which(
            !is.na(labels_clean) &
              startsWith(labels_lower, term_lower)
          )
        }
        
        if (length(matches_idx) == 0) {
          matches_idx <- which(
            !is.na(labels_clean) &
              grepl(search_term, labels_clean, ignore.case = TRUE, fixed = TRUE)
          )
        }
      }
      
      if (length(matches_idx) == 0) {
        showNotification(
          paste0("Gene '", search_term, "' not found"),
          type = "warning",
          duration = 3
        )
        return()
      }
      
      if (length(matches_idx) > 20) {
        showNotification(
          paste0(
            "Too many matches (", length(matches_idx),
            "). Please type a more specific gene name."
          ),
          type = "warning",
          duration = 4
        )
        return()
      }
      
      matched_ids <- df$point_id[matches_idx]
      current <- clicked_points()
      new_ids <- setdiff(matched_ids, current)
      
      clicked_points(c(current, new_ids))
      
      showNotification(
        paste0("Labeled ", length(new_ids), " match(es) for '", search_term, "'"),
        type = "message",
        duration = 3
      )
    })
    observeEvent(
      suppressWarnings(plotly::event_data("plotly_click", source = "volcano_click", priority = "event")),
      {
        click <- suppressWarnings(plotly::event_data("plotly_click", source = "volcano_click"))
        
        if (!is.null(click) && !is.null(click$customdata) && length(click$customdata) > 0) {
          point_id <- suppressWarnings(as.integer(click$customdata[[1]]))
          
          if (!is.na(point_id)) {
            current <- clicked_points()
            
            if (point_id %in% current) {
              clicked_points(setdiff(current, point_id))
            } else {
              clicked_points(c(current, point_id))
            }
          }
        }
      },
      ignoreInit = TRUE,
      ignoreNULL = TRUE
    )
    output$volcano_plot <- plotly::renderPlotly({
      df <- build_volcano_df()
      req(!is.null(df))
      
      volcano_data(df)
      
      comp <- input$comparison
      flip <- flip_state()
      fc_t <- input$fc_thr
      p_t  <- input$p_thr
      clicked <- clicked_points()
      
      # counts from full data
      n_up   <- sum(df$colour == "Up", na.rm = TRUE)
      n_down <- sum(df$colour == "Down", na.rm = TRUE)
      n_ns   <- sum(df$colour == "NS", na.rm = TRUE)

      on_off_baseline_y <- 0.05

      df$plot_x <- df$x
      df$plot_y <- df$y
      df$plot_y_is_display_only <- FALSE

      if ("on_off_display_x" %in% names(df)) {
        display_x <- df$on_off_display_x
        if (isTRUE(flip)) display_x <- -display_x
        needs_x <- df$on_off & !is.finite(df$plot_x) & is.finite(display_x)
        df$plot_x[needs_x] <- display_x[needs_x]
      }

      needs_y <- df$on_off & !is.finite(df$plot_y)
      df$plot_y[needs_y] <- on_off_baseline_y
      df$plot_y_is_display_only[needs_y] <- TRUE
      
      # Remove rows with missing coordinates
      df_plot <- df[is.finite(df$plot_x) & is.finite(df$plot_y), , drop = FALSE]
      req(nrow(df_plot) > 0)
      
      # Optionally subsample NS points
      if (isTRUE(input$show_ns)) {
        sig_points <- df_plot[df_plot$colour != "NS", , drop = FALSE]
        ns_points  <- df_plot[df_plot$colour == "NS", , drop = FALSE]
        
        sample_pct <- input$ns_sample / 100
        n_sample <- max(100, floor(nrow(ns_points) * sample_pct))
        
        if (nrow(ns_points) > n_sample) {
          set.seed(42)
          ns_points <- ns_points[sample(nrow(ns_points), n_sample), , drop = FALSE]
        }
        
        df_plot <- rbind(sig_points, ns_points)
      } else {
        df_plot <- df_plot[df_plot$colour != "NS", , drop = FALSE]
      }

      flagged_keep <- integer(0)
      if (isTRUE(input$show_low_abundant)) {
        flagged_keep <- c(flagged_keep, df$point_id[!is.na(df$low_abundant) & df$low_abundant])
      }
      if (isTRUE(input$show_on_off)) {
        flagged_keep <- c(flagged_keep, df$point_id[!is.na(df$on_off) & df$on_off])
      }

      # Always keep explicitly selected and enabled flagged points visible,
      # even if they were filtered out by NS sampling.
      keep_ids <- unique(c(clicked, flagged_keep))
      if (length(keep_ids) > 0) {
        keep_rows <- df[df$point_id %in% keep_ids, , drop = FALSE]
        if (nrow(keep_rows)) {
          keep_rows <- keep_rows[!keep_rows$point_id %in% df_plot$point_id, , drop = FALSE]
          if (nrow(keep_rows)) {
            df_plot <- rbind(df_plot, keep_rows)
          }
        }
      }
      
      req(nrow(df_plot) > 0)
      df_plot$x <- df_plot$plot_x
      df_plot$y <- df_plot$plot_y
      
      # Point colors
      color_vec <- character(nrow(df_plot))
      color_vec[df_plot$colour == "Up"]   <- "#F76075"
      color_vec[df_plot$colour == "Down"] <- "#40BFC1"
      color_vec[df_plot$colour == "NS"]   <- "#BDBFBE"
      
      # Hover text
      df_plot$hover <- paste0(
        "<b>", df_plot$label, "</b><br>",
        "logFC: ", round(df_plot$x, 3), "<br>",
        "-log10(p): ", ifelse(df_plot$plot_y_is_display_only, "not available", round(df_plot$y, 3)), "<br>",
        "Raw mean: ", ifelse(is.finite(df_plot$raw_mean), round(df_plot$raw_mean, 3), "not available"), "<br>",
        "Class: ", df_plot$colour,
        ifelse(df_plot$low_abundant, "<br>Flag: low-intensity partial-missing", ""),
        ifelse(df_plot$on_off, paste0("<br>Flag: on/off (", df_plot$on_off_direction, ")"), ""),
        ifelse(df_plot$plot_y_is_display_only, "<br>Display note: no finite p-value; plotted at baseline", "")
      )
      
      x_span <- diff(range(df_plot$x, na.rm = TRUE))
      y_span <- diff(range(df_plot$y, na.rm = TRUE))
      
      x_margin <- if (is.finite(x_span) && x_span > 0) x_span * 0.1 else 1
      y_margin <- if (is.finite(y_span) && y_span > 0) y_span * 0.1 else 1
      
      x_range <- c(
        if (!is.na(input$x_min)) input$x_min else min(df_plot$x, na.rm = TRUE) - x_margin,
        if (!is.na(input$x_max)) input$x_max else max(df_plot$x, na.rm = TRUE) + x_margin
      )
      
      y_range <- c(
        if (!is.na(input$y_min)) input$y_min else min(0, min(df_plot$y, na.rm = TRUE) - y_margin * 0.5),
        if (!is.na(input$y_max)) input$y_max else max(df_plot$y, na.rm = TRUE) + y_margin
      )
      
      title_str <- paste0(
        if (flip) paste0("<- ", comp) else comp,
        if (df_plot$source_label[1] == "T-test") " [T-test]" else " [DEqMS]"
      )
      
      p <- plotly::plot_ly(
        data       = df_plot,
        x          = ~x,
        y          = ~y,
        type       = "scatter",
        mode       = "markers",
        customdata = ~point_id,
        hovertext  = ~hover,
        hoverinfo  = "text",
        marker     = list(
          color   = color_vec,
          size    = 7,
          opacity = 0.75,
          line    = list(width = 0.5, color = "white")
        ),
        showlegend = FALSE,
        source     = "volcano_click"
      ) |>
        plotly::layout(
          title = list(text = title_str, font = list(size = 16)),
          xaxis = list(
            title = "log2 Fold Change",
            range = x_range,
            zeroline = TRUE,
            zerolinecolor = "#888888"
          ),
          yaxis = list(
            title = "-log10(p-value)",
            range = y_range
          ),
          margin = list(l = 60, r = 40, t = 60, b = 60, pad = 10),
          shapes = list(
            list(
              type = "line",
              x0 = x_range[1], x1 = x_range[2],
              y0 = p_t, y1 = p_t,
              line = list(color = "red", width = 1, dash = "dash")
            ),
            list(
              type = "line",
              x0 = fc_t, x1 = fc_t,
              y0 = 0, y1 = y_range[2],
              line = list(color = "blue", width = 1, dash = "dash")
            ),
            list(
              type = "line",
              x0 = -fc_t, x1 = -fc_t,
              y0 = 0, y1 = y_range[2],
              line = list(color = "blue", width = 1, dash = "dash")
            )
          ),
          hovermode = "closest",
          showlegend = FALSE
        )

      if (isTRUE(input$show_low_abundant)) {
        low_rows <- df_plot[df_plot$low_abundant, , drop = FALSE]

        if (nrow(low_rows)) {
          p <- plotly::add_trace(
            p,
            data = low_rows,
            x = ~x,
            y = ~y,
            type = "scatter",
            mode = "markers",
            customdata = ~point_id,
            hovertext = ~hover,
            hoverinfo = "text",
            marker = list(
              symbol = "circle",
              size = 9,
              color = "rgba(0,0,0,0)",
              opacity = 1,
              line = list(width = 1.4, color = "#6B7280")
            ),
            name = "Low-intensity partial-missing",
            showlegend = FALSE,
            inherit = FALSE
          )
        }
      }

      if (isTRUE(input$show_on_off)) {
        onoff_rows <- df_plot[df_plot$on_off, , drop = FALSE]

        if (nrow(onoff_rows)) {
          p <- plotly::add_trace(
            p,
            data = onoff_rows,
            x = ~x,
            y = ~y,
            type = "scatter",
            mode = "markers",
            customdata = ~point_id,
            hovertext = ~hover,
            hoverinfo = "text",
            marker = list(
              symbol = "triangle-up-open",
              size = 17,
              color = "#374151",
              opacity = 1,
              line = list(width = 2.2, color = "#374151")
            ),
            name = "On/off",
            showlegend = FALSE,
            inherit = FALSE
          )
        }
      }
      
      # Keep counts visible even when non-significant points are sampled.
      count_annotations <- list(
        list(
          x = x_range[1],
          y = y_range[2],
          xref = "x",
          yref = "y",
          text = paste0("Down: ", n_down),
          showarrow = FALSE,
          xanchor = "left",
          yanchor = "top",
          font = list(size = 14, color = "#40BFC1"),
          bgcolor = "rgba(255,255,255,0.75)",
          bordercolor = "#40BFC1",
          borderwidth = 1
        ),
        list(
          x = x_range[2],
          y = y_range[2],
          xref = "x",
          yref = "y",
          text = paste0("Up: ", n_up),
          showarrow = FALSE,
          xanchor = "right",
          yanchor = "top",
          font = list(size = 14, color = "#F76075"),
          bgcolor = "rgba(255,255,255,0.75)",
          bordercolor = "#F76075",
          borderwidth = 1
        ),
        list(
          x = mean(x_range),
          y = y_range[2],
          xref = "x",
          yref = "y",
          text = paste0("NS: ", n_ns),
          showarrow = FALSE,
          xanchor = "center",
          yanchor = "top",
          font = list(size = 13, color = "#666666"),
          bgcolor = "rgba(255,255,255,0.70)",
          bordercolor = "#BDBFBE",
          borderwidth = 1
        )
      )
      
      clicked_annotations <- list()
      
      if (length(clicked) > 0) {
        annotation_rows <- df_plot[df_plot$point_id %in% clicked, , drop = FALSE]
        
        if (nrow(annotation_rows) > 0) {
          clicked_annotations <- lapply(seq_len(nrow(annotation_rows)), function(i) {
            row <- annotation_rows[i, ]
            
            list(
              x = row$x,
              y = row$y,
              text = row$label,
              showarrow = TRUE,
              arrowhead = 2,
              arrowsize = 1,
              arrowwidth = 1,
              arrowcolor = "black",
              ax = 20,
              ay = -30,
              font = list(size = 11, color = "black"),
              bgcolor = "rgba(255,255,255,0.85)",
              bordercolor = "black",
              borderwidth = 1
            )
          })
        }
      }
      
      p <- plotly::layout(
        p,
        annotations = c(count_annotations, clicked_annotations)
      )

      p <- plotly::event_register(p, "plotly_click")
      
      plotly::config(
        p,
        displayModeBar = TRUE,
        displaylogo = FALSE,
        toImageButtonOptions = list(
          format = "svg",
          filename = paste0("volcano_", comp)
        ),
        modeBarButtonsToRemove = c("select2d", "lasso2d", "autoScale2d")
      )
    })
    output$dl_plot <- downloadHandler(
      filename = function() {
        comp <- if (is.null(input$comparison) || !nzchar(input$comparison)) "comparison" else input$comparison
        suffix <- paste0("volcano_", gsub("[^A-Za-z0-9_]+", "_", comp))
        download_filename(app_state, suffix, "svg")
      },
      content = function(file) {
        showNotification(
          "Use the camera icon in the plot toolbar to export as SVG",
          type = "message",
          duration = 3
        )
      }
    )
    
    output$dl_csv <- downloadHandler(
      filename = function() {
        comp <- if (is.null(input$comparison) || !nzchar(input$comparison)) "comparison" else input$comparison
        suffix <- paste0("volcano_data_", gsub("[^A-Za-z0-9_]+", "_", comp))
        download_filename(app_state, suffix, "csv")
      },
      content = function(file) {
        df <- volcano_data()
        req(!is.null(df))
        
        export_cols <- intersect(
          c(
            "point_id", "label", "x", "y", "logFC", "colour", "Genes",
            "ProtDesc", "Protein.Names", "negLog10p",
            "negLog10_sca_adjPval", "sca.adj.pval",
            "raw_mean", "low_abundant", "on_off", "on_off_direction"
          ),
          names(df)
        )
        
        utils::write.csv(df[, export_cols, drop = FALSE], file, row.names = FALSE)
      }
    )

    output$proceed_btn <- renderUI({
      if (!isTRUE(app_state$deqms_done)) return(NULL)
      actionButton(ns("btn_proceed"), "Proceed to Enrichment",
                   class = "btn-success btn-block app-proceed-btn")
    })

    observeEvent(input$btn_proceed, {
      app_state$active_tab <- "enrichment"
    })
  })
}
