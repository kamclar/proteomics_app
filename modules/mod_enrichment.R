# GO, KEGG, and Reactome enrichment for DEqMS results


mod_enrichment_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Step 6. Enrichment Analysis"),
    uiOutput(ns("guard")),
    uiOutput(ns("main_ui"))
  )
}

mod_enrichment_server <- function(id, app_state) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    enrich_results <- reactiveVal(list())
    kegg_cache     <- reactiveVal(NULL)
    reactome_cache <- reactiveVal(read_reactome_cache())
    installing     <- reactiveVal(FALSE)

    selected_dbs <- reactive({
      c(
        if (isTRUE(input$run_go_bp)) "GO_BP",
        if (isTRUE(input$run_go_mf)) "GO_MF",
        if (isTRUE(input$run_go_cc)) "GO_CC",
        if (isTRUE(input$run_kegg))  "KEGG",
        if (isTRUE(input$run_react)) "Reactome"
      )
    })

    output$guard <- renderUI({
      if (!isTRUE(app_state$deqms_done)) {
        tags$div(class = "alert alert-warning",
                 " Please complete Step 4 (DEqMS) first.")
      }
    })

    output$main_ui <- renderUI({
      req(isTRUE(app_state$deqms_done))
      tagList(
        fluidRow(
          column(4,
            wellPanel(
              h4("Settings"),

              h5("Comparisons"),
              div(style = "display:flex; gap:8px; margin-bottom:4px;",
                actionButton(ns("sel_all"),   "All",  class = "btn-xs btn-default"),
                actionButton(ns("desel_all"), "None", class = "btn-xs btn-default")
              ),
              uiOutput(ns("pair_checkboxes")),

              hr(),
              h5("Thresholds"),
              numericInput(ns("fc_cut"),  "log2FC cutoff",    value = 0.5, min = 0, step = 0.1),
              numericInput(ns("p_cut"),   "sca.adj.pval cutoff", value = 0.05, min = 0.001, step = 0.005),

              hr(),
              h5("Databases"),
              checkboxInput(ns("run_go_bp"), "GO Biological Process", value = TRUE),
              checkboxInput(ns("run_go_mf"), "GO Molecular Function",  value = TRUE),
              checkboxInput(ns("run_go_cc"), "GO Cellular Component",  value = FALSE),
              checkboxInput(ns("run_kegg"),  "KEGG Pathways",          value = TRUE),
              checkboxInput(ns("run_react"), "Reactome Pathways",      value = TRUE),

              uiOutput(ns("dependency_status")),

              hr(),
              h5("Performance"),
              numericInput(ns("workers"),
                           "Parallel workers (1 = sequential)",
                           value = max(1, parallel::detectCores() - 1),
                           min = 1, max = parallel::detectCores(), step = 1),
              helpText("More workers = faster for many comparisons."),

              hr(),
              h5("KEGG mode"),
              radioButtons(ns("kegg_mode"), label = NULL,
                choices = list(
                  "Prefetch once (fastest for multiple runs)" = "prefetch",
                  "Online (standard enrichKEGG)"              = "online"
                ),
                selected = "prefetch"
              ),
              tags$div(class = "alert alert-warning", style = "font-size:12px; padding:6px;",
                       " KEGG prefetch requires an internet connection and may be slow on first use."),

              hr(),
              h5("Reactome cache"),
              actionButton(ns("btn_reactome_cache"), "Build Reactome cache",
                           class = "btn-default btn-block"),
              uiOutput(ns("reactome_cache_status")),

              hr(),
              actionButton(ns("btn_install_selected"), "Install Selected Enrichment Packages",
                           class = "btn-default btn-block"),
              br(),
              actionButton(ns("btn_run"), "Run Enrichment",
                           class = "btn-primary btn-block"),
              br(),
              uiOutput(ns("run_status")),
              hr(),
              actionButton(ns("btn_check_cores"), "Check CPU cores",
                           class = "btn-default btn-xs btn-block")
            )
          ),
          column(8,
            tabsetPanel(
              tabPanel("Results Table",
                br(),
                uiOutput(ns("result_nav")),
                br(),
                DT::dataTableOutput(ns("result_tbl"))
              ),
              tabPanel("Dotplot",
                br(),
                uiOutput(ns("dotplot_nav")),
                br(),
                plotOutput(ns("dotplot"), height = "500px"),
                br(),
                downloadButton(ns("dl_dotplot"), "Download Plot (SVG)",
                               class = "btn-default btn-sm")
              ),
              tabPanel("Summary",
                br(),
                DT::dataTableOutput(ns("summary_tbl"))
              )
            )
          )
        ),
        fluidRow(
          column(12,
            div(class = "app-action-row",
              downloadButton(ns("dl_all"), "Download All Enrichment CSVs (ZIP)",
                             class = "btn-default"),
              div(class = "app-proceed-wrap",
                uiOutput(ns("proceed_btn"))
              )
            )
          )
        )
      )
    })

    output$dependency_status <- renderUI({
      dbs <- selected_dbs()
      if (!length(dbs)) return(NULL)

      status <- enrichment_dependency_status(dbs)
      if (!length(status$missing)) {
        return(tags$div(class = "alert alert-success", style = "font-size:12px; padding:6px;",
                        " Selected enrichment packages are installed."))
      }

      note <- if ("Reactome" %in% dbs) {
        "Reactome is a large optional download because it includes reactome.db."
      } else if (any(grepl("^GO_", dbs))) {
        "GO is a large optional download because it includes human and GO annotation databases."
      } else {
        "KEGG support is downloaded only when enrichment is used."
      }

      tags$div(class = "alert alert-warning", style = "font-size:12px; padding:6px;",
               tags$b(" Missing optional packages: "),
               paste(status$missing, collapse = ", "),
               tags$br(),
               note)
    })

    observeEvent(input$btn_install_selected, {
      dbs <- selected_dbs()
      if (!length(dbs)) {
        showNotification("Select at least one enrichment database first.", type = "warning")
        return()
      }

      status <- enrichment_dependency_status(dbs)
      if (!length(status$missing)) {
        showNotification("Selected enrichment packages are already installed.", type = "message")
        return()
      }

      installing(TRUE)
      on.exit(installing(FALSE), add = TRUE)

      withProgress(message = "Installing optional enrichment packages...", value = 0.2, {
        ok <- install_enrichment_packages(status$missing)
        incProgress(0.8)
      })

      if (isTRUE(ok)) {
        showNotification("Enrichment packages installed. You can run enrichment now.", type = "message", duration = 8)
      } else {
        showNotification("Package installation failed. Check internet access and Bioconductor availability.", type = "error", duration = 12)
      }
    })
    # Let the user choose exactly which DEqMS result tables to send downstream.
    output$pair_checkboxes <- renderUI({
      req(isTRUE(app_state$deqms_done))
      nms <- names(app_state$deqms_results)
      checkboxGroupInput(ns("pairs"), label = NULL,
                         choices = nms, selected = nms)
    })

    observeEvent(input$sel_all, {
      req(isTRUE(app_state$deqms_done))
      updateCheckboxGroupInput(session, "pairs",
        selected = names(app_state$deqms_results))
    })
    observeEvent(input$desel_all, {
      updateCheckboxGroupInput(session, "pairs", selected = character(0))
    })

    observeEvent(input$btn_check_cores, {
      cores <- parallel::detectCores(logical = TRUE)
      cores_text <- if (is.na(cores)) "Could not detect CPU core count." else {
        paste0("Detected ", cores, " logical CPU core(s).")
      }
      showNotification(
        paste(cores_text, "Current enrichment worker setting:", input$workers),
        type = "message",
        duration = 6
      )
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
      tags$div(class = "alert alert-info", style = "font-size:12px; padding:6px;",
        sprintf("Local cache: %d mapped gene(s), %d pathway(s). Background is always the current dataset.",
                st$mapped_genes, st$pathways))
    })

    observeEvent(input$btn_run, {
      req(isTRUE(app_state$deqms_done), length(input$pairs) > 0)

      dbs_to_run <- selected_dbs()
      if (!length(dbs_to_run)) {
        showNotification("Select at least one database.", type = "warning")
        return()
      }

      status <- enrichment_dependency_status(dbs_to_run)
      if (length(status$missing)) {
        showNotification(
          paste("Install optional enrichment packages first:", paste(status$missing, collapse = ", ")),
          type = "warning",
          duration = 12
        )
        return()
      }

      kegg_data <- kegg_cache()
      if (input$run_kegg && input$kegg_mode == "prefetch" && is.null(kegg_data)) {
        withProgress(message = "Prefetching KEGG data (one-time, may be slow)...", value = 0.5, {
          kegg_data <- tryCatch(prefetch_kegg_hsa_local(), error = function(e) {
            showNotification(paste("KEGG prefetch failed:", conditionMessage(e),
                                   "Switch KEGG mode to Online or retry prefetch."),
                             type = "error", duration = 12)
            NULL
          })
        })
        if (is.null(kegg_data)) return()
        kegg_cache(kegg_data)
      }

      pairs_sel  <- input$pairs
      deqms_res  <- app_state$deqms_results
      fc_cut     <- input$fc_cut
      p_cut      <- input$p_cut
      workers    <- max(1L, as.integer(input$workers))
      kegg_mode  <- input$kegg_mode
      kegg_d     <- kegg_data
      universe_symbols <- dataset_gene_symbols(app_state$parsed_data)
      if (!length(universe_symbols)) {
        showNotification("Enrichment requires a Genes column to define the current dataset background. No gene symbols were found.",
                         type = "error", duration = 12)
        return()
      }
      if ("Reactome" %in% dbs_to_run) {
        cache_res <- tryCatch(build_reactome_cache_for_symbols(universe_symbols), error = function(e) {
          showNotification(paste("Reactome cache failed:", conditionMessage(e)), type = "error", duration = 12)
          NULL
        })
        if (is.null(cache_res)) return()
        if (!is.null(cache_res)) reactome_cache(cache_res$cache)
      }
      universe_ids <- if (requires_orgdb_mapping(dbs_to_run)) {
        symbols_to_entrez_local(universe_symbols)
      } else {
        symbols_to_entrez_kegg_online(universe_symbols)
      }
      if (!length(universe_ids)) {
        showNotification("Enrichment could not map the current dataset background to gene IDs. Analysis was not run.",
                         type = "error", duration = 12)
        return()
      }

      withProgress(message = "Running enrichment analysis...", value = 0, {

        use_future <- workers > 1L && requireNamespace("future", quietly = TRUE)
        if (use_future) {
          future::plan(future::multisession, workers = workers)
          on.exit(future::plan(future::sequential), add = TRUE)
        }

        results_new <- list()

        for (i in seq_along(pairs_sel)) {
          pair <- pairs_sel[[i]]
          incProgress(1 / length(pairs_sel),
                      detail = paste0(i, "/", length(pairs_sel), " - ", pair))

          df <- deqms_res[[pair]]
          if (is.null(df)) next

          # Gene symbols
          up_genes   <- get_sig_genes(df, fc_cut,  p_cut, direction = "up")
          down_genes <- get_sig_genes(df, fc_cut,  p_cut, direction = "down")

          if (requires_orgdb_mapping(dbs_to_run)) {
            up_ids   <- symbols_to_entrez_local(up_genes)
            down_ids <- symbols_to_entrez_local(down_genes)
          } else {
            up_ids   <- symbols_to_entrez_kegg_online(up_genes)
            down_ids <- symbols_to_entrez_kegg_online(down_genes)
          }

          up_enrich   <- run_enrichment_bundle(up_ids,   dbs_to_run, kegg_mode, kegg_d,
                                                universe_ids, reactome_cache())
          down_enrich <- run_enrichment_bundle(down_ids, dbs_to_run, kegg_mode, kegg_d,
                                                universe_ids, reactome_cache())

          results_new[[paste0(pair, "__Up")]]   <- up_enrich
          results_new[[paste0(pair, "__Down")]] <- down_enrich
        }

        existing <- enrich_results()
        enrich_results(c(existing, results_new))
        app_state$enrich_results <- enrich_results()
        app_state$enrichment_done <- TRUE
      })
    })

    output$run_status <- renderUI({
      res <- enrich_results()
      if (!length(res)) return(NULL)
      tags$div(class = "alert alert-success",
               sprintf(" %d comparison-direction bundles computed.", length(res)))
    })
    enrich_choices <- reactive({
      res <- enrich_results()
      if (!length(res)) return(character(0))
      # Flatten the nested result list for a single selector.
      flat <- lapply(names(res), function(pd) {
        bundle <- res[[pd]]
        if (!length(bundle)) return(NULL)
        paste0(pd, "__", names(bundle))
      })
      unlist(Filter(Negate(is.null), flat))
    })

    output$result_nav <- renderUI({
      ch <- enrich_choices()
      if (!length(ch)) return(tags$p("No results yet. Run enrichment first."))
      selectInput(ns("sel_result"), "View result:", choices = ch, selected = ch[1])
    })

    output$dotplot_nav <- renderUI({
      ch <- enrich_choices()
      if (!length(ch)) return(NULL)
      selectInput(ns("sel_dotplot"), "View dotplot:", choices = ch, selected = ch[1])
    })
    get_selected_result <- function(key) {
      if (is.null(key) || !nchar(key)) return(NULL)
      # Keys are built as pair__direction__database.
      parts <- strsplit(key, "__")[[1]]
      if (length(parts) < 3) return(NULL)
      pd  <- paste(parts[1], parts[2], sep = "__")
      db  <- paste(parts[3:length(parts)], collapse = "__")
      enrich_results()[[pd]][[db]]
    }

    output$result_tbl <- DT::renderDataTable({
      req(length(enrich_choices()) > 0, input$sel_result)
      obj <- get_selected_result(input$sel_result)
      if (is.null(obj)) return(data.frame(Message = "No significant terms found."))
      df <- as.data.frame(obj)
      if (!nrow(df)) return(data.frame(Message = "No significant terms found."))
      DT::datatable(df, rownames = FALSE,
                    options = list(pageLength = 15, scrollX = TRUE, dom = "ftip"))
    })
    # Dotplot
    output$dotplot <- renderPlot({
      req(length(enrich_choices()) > 0, input$sel_dotplot)
      obj <- get_selected_result(input$sel_dotplot)
      if (is.null(obj)) {
        plot.new(); text(0.5, 0.5, "No significant terms for this selection.", cex = 1.2)
        return()
      }
      df <- tryCatch(as.data.frame(obj), error = function(e) NULL)
      if (is.null(df) || !nrow(df)) {
        plot.new(); text(0.5, 0.5, "No significant terms for this selection.", cex = 1.2)
        return()
      }
      tryCatch(
        print(enrichplot::dotplot(obj) +
                ggplot2::ggtitle(input$sel_dotplot) +
                ggplot2::theme(text = ggplot2::element_text(size = 9),
                               axis.text.y = ggplot2::element_text(size = 8))),
        error = function(e) {
          plot.new(); text(0.5, 0.5, paste("Plot error:", conditionMessage(e)), cex = 0.9)
        }
      )
    })

    output$dl_dotplot <- downloadHandler(
      filename = function() download_filename(app_state, "dotplot", "svg"),
      content  = function(file) {
        req(input$sel_dotplot)
        obj <- get_selected_result(input$sel_dotplot)
        req(!is.null(obj))
        p <- enrichplot::dotplot(obj) + ggplot2::ggtitle(input$sel_dotplot)
        ggplot2::ggsave(file, plot = p, device = "svg", width = 7, height = 5)
      }
    )
    output$summary_tbl <- DT::renderDataTable({
      res <- enrich_results()
      if (!length(res)) return(data.frame())
      rows <- lapply(names(res), function(pd) {
        bundle <- res[[pd]]
        lapply(names(bundle), function(db) {
          obj <- bundle[[db]]
          n   <- if (!is.null(obj)) nrow(as.data.frame(obj)) else 0L
          data.frame(Comparison_Direction = pd, Database = db,
                     Significant_Terms = n, stringsAsFactors = FALSE)
        })
      })
      df <- do.call(rbind, unlist(rows, recursive = FALSE))
      DT::datatable(df, rownames = FALSE,
                    options = list(pageLength = 20, dom = "ftip"))
    })
    output$dl_all <- downloadHandler(
      filename = function() download_filename(app_state, "enrichment_results", "zip"),
      content  = function(file) {
        res <- enrich_results()
        tmp <- tempdir()
        fnames <- c()
        for (pd in names(res)) {
          bundle <- res[[pd]]
          for (db in names(bundle)) {
            obj <- bundle[[db]]
            if (is.null(obj)) next
            df  <- as.data.frame(obj)
            if (!nrow(df)) next
            fname <- file.path(tmp, paste0(pd, "__", db, ".csv"))
            utils::write.csv(df, fname, row.names = FALSE)
            fnames <- c(fnames, fname)
          }
        }
        if (length(fnames)) zip::zip(zipfile = file, files = fnames, mode = "cherry-pick")
      }
    )
    output$proceed_btn <- renderUI({
      if (!isTRUE(app_state$enrichment_done)) return(NULL)
      actionButton(ns("btn_proceed"), "Proceed to Stable Proteins",
                   class = "btn-success btn-lg app-proceed-btn")
    })

    observeEvent(input$btn_proceed, {
      app_state$active_tab <- "background"
    })
  })
}
# Enrichment helpers
ENRICHMENT_PACKAGE_GROUPS <- list(
  GO = c("clusterProfiler", "org.Hs.eg.db", "AnnotationDbi", "enrichplot"),
  KEGG = c("clusterProfiler", "KEGGREST", "enrichplot"),
  Reactome = c("ReactomePA", "reactome.db", "org.Hs.eg.db", "AnnotationDbi", "enrichplot")
)

enrichment_groups_for_dbs <- function(dbs) {
  groups <- character(0)
  if (any(grepl("^GO_", dbs))) groups <- c(groups, "GO")
  if ("KEGG" %in% dbs) groups <- c(groups, "KEGG")
  if ("Reactome" %in% dbs) groups <- c(groups, "Reactome")
  unique(groups)
}

enrichment_packages_for_dbs <- function(dbs) {
  groups <- enrichment_groups_for_dbs(dbs)
  unique(unlist(ENRICHMENT_PACKAGE_GROUPS[groups], use.names = FALSE))
}

enrichment_dependency_status <- function(dbs) {
  packages <- enrichment_packages_for_dbs(dbs)
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  list(required = packages, missing = missing)
}

install_enrichment_packages <- function(packages) {
  packages <- unique(packages)
  if (!length(packages)) return(TRUE)

  tryCatch({
    repos <- getOption("repos")
    if (is.null(repos) || identical(unname(repos["CRAN"]), "@CRAN@")) {
      options(repos = c(CRAN = "https://cloud.r-project.org"))
    }
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install(packages, ask = FALSE, update = FALSE)
    remaining <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
    !length(remaining)
  }, error = function(e) {
    message("Optional enrichment package installation failed: ", conditionMessage(e))
    FALSE
  })
}

requires_orgdb_mapping <- function(dbs) {
  any(grepl("^GO_", dbs)) || "Reactome" %in% dbs
}

# Symbol caches are populated once per session; annotation lookups are slow.
.VALID_SYMBOLS    <- NULL
.SYMBOL2ENTREZ    <- new.env(parent = emptyenv())
.KEGG_SYMBOL2ENTREZ <- new.env(parent = emptyenv())

get_valid_symbols <- function() {
  if (is.null(.VALID_SYMBOLS)) {
    .VALID_SYMBOLS <<- AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "SYMBOL")
  }
  .VALID_SYMBOLS
}

extract_first_symbol <- function(x) {
  x <- x[!is.na(x)]
  x <- trimws(x)
  x <- x[x != "" & x != "NA"]
  if (!length(x)) return(character(0))
  x <- sub("[;|,].*$", "", x)
  unique(trimws(x)[trimws(x) != ""])
}

symbols_to_entrez_local <- function(symbols) {
  symbols <- unique(symbols)
  valid   <- get_valid_symbols()
  symbols <- symbols[symbols %in% valid]
  if (!length(symbols)) return(character(0))

  have    <- vapply(symbols, exists, logical(1), envir = .SYMBOL2ENTREZ, inherits = FALSE)
  missing <- symbols[!have]

  if (length(missing)) {
    mapped <- suppressMessages(AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db,
      keys     = missing,
      columns  = "ENTREZID",
      keytype  = "SYMBOL"
    ))
    mapped <- mapped[!is.na(mapped$ENTREZID) & mapped$ENTREZID != "", , drop = FALSE]
    sp <- split(mapped$ENTREZID, mapped$SYMBOL)
    for (s in names(sp)) assign(s, unique(as.character(sp[[s]])), envir = .SYMBOL2ENTREZ)
    nohit <- setdiff(missing, names(sp))
    for (s in nohit) assign(s, character(0), envir = .SYMBOL2ENTREZ)
  }

  ids <- unique(unlist(mget(symbols, envir = .SYMBOL2ENTREZ, ifnotfound = list(character(0)))))
  ids <- as.character(ids)
  unique(ids[!is.na(ids) & ids != ""])
}

build_kegg_symbol_map <- function() {
  if (length(ls(envir = .KEGG_SYMBOL2ENTREZ, all.names = TRUE))) return(invisible(TRUE))

  genes <- KEGGREST::keggList("hsa")
  ids <- sub("^hsa:", "", names(genes))
  aliases <- sub(";.*$", "", unname(genes))
  aliases <- strsplit(aliases, ",\\s*")

  for (i in seq_along(ids)) {
    vals <- unique(trimws(aliases[[i]]))
    vals <- vals[nzchar(vals)]
    for (symbol in vals) {
      key <- toupper(symbol)
      current <- get0(key, envir = .KEGG_SYMBOL2ENTREZ, inherits = FALSE, ifnotfound = character(0))
      assign(key, unique(c(current, ids[[i]])), envir = .KEGG_SYMBOL2ENTREZ)
    }
  }

  invisible(TRUE)
}

symbols_to_entrez_kegg_online <- function(symbols) {
  symbols <- unique(symbols)
  symbols <- symbols[!is.na(symbols) & symbols != ""]
  if (!length(symbols)) return(character(0))

  build_kegg_symbol_map()
  keys <- toupper(symbols)
  ids <- unique(unlist(mget(keys, envir = .KEGG_SYMBOL2ENTREZ, ifnotfound = list(character(0)))))
  ids <- as.character(ids)
  unique(ids[!is.na(ids) & ids != ""])
}

get_sig_genes <- function(df, fc_cut, p_cut, direction = "up") {
  req_cols <- c("logFC", "sca.adj.pval", "Genes")
  if (!all(req_cols %in% names(df))) return(character(0))
  if (direction == "up") {
    sub <- df[df$logFC >= fc_cut & df$sca.adj.pval <= p_cut & !is.na(df$Genes), ]
  } else {
    sub <- df[df$logFC <= -fc_cut & df$sca.adj.pval <= p_cut & !is.na(df$Genes), ]
  }
  extract_first_symbol(sub$Genes)
}

prefetch_kegg_hsa_local <- function() {
  link <- KEGGREST::keggLink("pathway", "hsa")
  t2g  <- data.frame(
    term = sub("^path:", "", unname(link)),
    gene = sub("^hsa:", "",  names(link)),
    stringsAsFactors = FALSE
  )
  nm  <- KEGGREST::keggList("pathway", "hsa")
  t2n <- data.frame(
    term = sub("^path:", "", names(nm)),
    name = unname(nm),
    stringsAsFactors = FALSE
  )
  list(t2g = t2g, t2n = t2n)
}

safe_enrich <- function(expr) {
  tryCatch(expr, error = function(e) {
    message("Enrichment error: ", conditionMessage(e))
    NULL
  })
}

MIN_GENES_ENRICH <- 5L

run_enrichment_bundle <- function(gene_ids, dbs, kegg_mode = "prefetch", kegg_data = NULL,
                                  universe_gene_ids = character(0),
                                  reactome_cache_data = read_reactome_cache()) {
  gene_ids <- unique(as.character(gene_ids))
  gene_ids <- gene_ids[!is.na(gene_ids) & gene_ids != ""]
  universe_gene_ids <- unique(as.character(universe_gene_ids))
  universe_gene_ids <- universe_gene_ids[!is.na(universe_gene_ids) & universe_gene_ids != ""]
  if (length(gene_ids) < MIN_GENES_ENRICH) return(list())
  if (!length(universe_gene_ids)) {
    stop("Enrichment requires a non-empty current-dataset universe.")
  }

  res <- list()

  if ("GO_BP" %in% dbs) {
    res$GO_BP <- safe_enrich(clusterProfiler::enrichGO(
      gene = gene_ids, OrgDb = org.Hs.eg.db::org.Hs.eg.db,
      universe = universe_gene_ids,
      keyType = "ENTREZID", ont = "BP",
      pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.20
    ))
  }
  if ("GO_MF" %in% dbs) {
    res$GO_MF <- safe_enrich(clusterProfiler::enrichGO(
      gene = gene_ids, OrgDb = org.Hs.eg.db::org.Hs.eg.db,
      universe = universe_gene_ids,
      keyType = "ENTREZID", ont = "MF",
      pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.20
    ))
  }
  if ("GO_CC" %in% dbs) {
    res$GO_CC <- safe_enrich(clusterProfiler::enrichGO(
      gene = gene_ids, OrgDb = org.Hs.eg.db::org.Hs.eg.db,
      universe = universe_gene_ids,
      keyType = "ENTREZID", ont = "CC",
      pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.20
    ))
  }

  if ("KEGG" %in% dbs) {
    if (kegg_mode == "prefetch" && !is.null(kegg_data)) {
      plain <- sub("^hsa:", "", gene_ids)
      universe_plain <- sub("^hsa:", "", universe_gene_ids)
      res$KEGG <- safe_enrich(clusterProfiler::enricher(
        gene = plain, TERM2GENE = kegg_data$t2g, TERM2NAME = kegg_data$t2n,
        universe = universe_plain,
        pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.20
      ))
    } else if (kegg_mode == "prefetch" && is.null(kegg_data)) {
      stop("KEGG prefetch mode requires prefetched KEGG data.")
    } else {
      args <- list(gene = gene_ids, organism = "hsa", keyType = "ncbi-geneid",
                   universe = universe_gene_ids,
                   pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.20)
      res$KEGG <- safe_enrich(do.call(clusterProfiler::enrichKEGG, args))
    }
  }

  if ("Reactome" %in% dbs) {
    reactome_terms <- reactome_entrez_term_tables_from_cache(reactome_cache_data)
    if (nrow(reactome_terms$t2g)) {
      res$Reactome <- safe_enrich(clusterProfiler::enricher(
        gene = gene_ids,
        universe = universe_gene_ids,
        TERM2GENE = reactome_terms$t2g,
        TERM2NAME = reactome_terms$t2n,
        pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.20
      ))
    } else {
      res$Reactome <- safe_enrich(ReactomePA::enrichPathway(
        gene = gene_ids, organism = "human",
        universe = universe_gene_ids,
        pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.20,
        readable = FALSE
      ))
    }
  }

  Filter(Negate(is.null), res)
}
