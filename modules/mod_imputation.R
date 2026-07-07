# Imputation tab: Perseus-style or imputeLCMD MNAR methods, kNN, or mixed mode.
mod_imputation_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Step 2. Imputation"),
    uiOutput(ns("guard")),
    uiOutput(ns("main_ui"))
  )
}

mod_imputation_server <- function(id, app_state) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    imputed_data <- reactiveVal(NULL)
    output$guard <- renderUI({
      if (!isTRUE(app_state$upload_done)) {
        tags$div(class = "alert alert-warning",
                 " Please upload and parse a file in Step 1 first.")
      }
    })

    output$main_ui <- renderUI({
      req(isTRUE(app_state$upload_done))
      tagList(
        fluidRow(
          column(4,
            wellPanel(
              h4("Imputation Settings"),

              checkboxInput(ns("skip_imputation"), "Skip - data is already imputed",
                            value = FALSE),

              conditionalPanel(
                condition = paste0("!input['", ns("skip_imputation"), "']"),

                hr(),
                radioButtons(ns("method"), "Method",
                  choices = list(
                    "MNAR - left-censored"                                        = "mnar_min",
                    "MNAR - left-censored per group"                              = "mnar_group_min",
                    "MAR - k-Nearest Neighbours (kNN)"                             = "mar_knn",
                    "Mixed - auto-detect MNAR/MAR per protein"                     = "mixed"
                  ),
                  selected = "mnar_min"
                ),

                conditionalPanel(
                  condition = paste0("input['", ns("method"), "'] == 'mnar_min' || ",
                                     "input['", ns("method"), "'] == 'mnar_group_min' || ",
                                     "input['", ns("method"), "'] == 'mixed'"),
                  hr(),
                  h5("MNAR settings"),
                  radioButtons(ns("mnar_backend"), "MNAR method",
                    choices = list("MNAR (Perseus default: width = 0.3, downshift = 1.8)" = "Perseus",
                                   "QRILC (imputeLCMD)"                = "QRILC",
                                   "MinProb (imputeLCMD)"              = "MinProb"),
                    selected = "Perseus"),
                  conditionalPanel(
                    condition = paste0("input['", ns("mnar_backend"), "'] == 'Perseus'"),
                    helpText("Default Perseus setting: width 0.3, downshift 1.8. Missing values are drawn from a low-intensity normal distribution."),
                    tags$details(
                      style = "border: 1px solid #d9e2ec; border-radius: 4px; padding: 8px 10px; background: #f8fafc; margin-top: 6px;",
                      tags$summary(
                        style = "cursor: pointer; font-weight: 600; color: #102a43;",
                        "Advanced Perseus settings (click to expand)"
                      ),
                      sliderInput(ns("mnar_perseus_width"), "Width",
                                  min = 0.05, max = 1.0, value = 0.3, step = 0.05),
                      sliderInput(ns("mnar_perseus_downshift"), "Downshift",
                                  min = 0.1, max = 4.0, value = 1.8, step = 0.1)
                    )
                  ),
                  conditionalPanel(
                    condition = paste0("input['", ns("mnar_backend"), "'] != 'Perseus'"),
                    sliderInput(ns("mnar_tune_sigma"), "Distribution width (tune.sigma)",
                                min = 0.1, max = 2.0, value = 1.0, step = 0.1)
                  ),
                  conditionalPanel(
                    condition = paste0("input['", ns("mnar_backend"), "'] == 'MinProb'"),
                    sliderInput(ns("mnar_q"), "MinProb quantile (q)",
                                min = 0.001, max = 0.05, value = 0.01, step = 0.001)
                  )
                ),

                conditionalPanel(
                  condition = paste0("input['", ns("method"), "'] == 'mar_knn' || input['", ns("method"), "'] == 'mixed'"),
                  hr(),
                  h5("MAR / Mixed settings"),
                  numericInput(ns("knn_k"), "kNN - number of neighbours", value = 5, min = 2, max = 20),
                  conditionalPanel(
                    condition = paste0("input['", ns("method"), "'] == 'mixed'"),
                    sliderInput(ns("mnar_threshold"),
                                "Mixed: missing fraction -> MNAR if above",
                                min = 0.1, max = 0.9, value = 0.5, step = 0.05),
                    uiOutput(ns("mnar_threshold_hint"))
                  )
                ),

                conditionalPanel(
                  condition = paste0("input['", ns("method"), "'] == 'mar_knn' || input['", ns("method"), "'] == 'mixed'"),

                  checkboxInput(ns("per_group_knn"), "Compute per group", value = FALSE),

                  conditionalPanel(
                    condition = paste0("input['", ns("method"), "'] == 'mar_knn' && input['", ns("per_group_knn"), "']"),
                    helpText("When enabled, kNN is computed separately for each group, using only that group's own replicates. A protein with 50% or more missing values within a group is not sent to kNN (unreliable there) - it gets an MNAR downshift instead, and the app shows a notification. Disabled = kNN is computed across all groups pooled together.")
                  ),
                  conditionalPanel(
                    condition = paste0("input['", ns("method"), "'] == 'mixed' && input['", ns("per_group_knn"), "']"),
                    radioButtons(ns("mixed_group_strategy"), "Group strategy",
                      choices = list("Smart / hybrid (recommended)" = "pooled",
                                     "Fully per group"               = "per_group"),
                      selected = "pooled"),
                    helpText("Smart/hybrid: classify missingness separately within each group, then run kNN once on the full matrix for MAR values."),
                    helpText("Fully per group: classify and run kNN separately within each group. This is stricter and can switch many sparse group blocks to MNAR.")
                  )
                ),

                hr(),
                h5("Row pre-filter"),
                sliderInput(ns("min_valid"),
                            "Min. valid values per protein row (across all samples)",
                            min = 1, max = 10, value = 1, step = 1),
                uiOutput(ns("prefilter_note")),
                helpText("Protein rows below this threshold are removed before imputation and from all downstream steps. Default 1 removes rows with no measured values."),
                hr(),
                h5("Reproducibility"),
                numericInput(ns("random_seed"), "Random seed", value = 12345, min = 1, step = 1),
                helpText("The seed makes random MNAR draws reproducible. Reusing it with the same settings gives the same imputed values."),
                hr(),
                tags$details(
                  style = "border: 1px solid #d9e2ec; border-radius: 4px; padding: 8px 10px; background: #f8fafc;",
                  tags$summary(
                    style = "cursor: pointer; font-weight: 600; color: #102a43;",
                    "Method help (click to expand)"
                  ),
                  tags$dl(
                    tags$dt("MNAR (Perseus default: width = 0.3, downshift = 1.8)"),
                    tags$dd("MNAR method using a low-intensity normal distribution with the common Perseus settings. This is the default for method comparison."),
                    tags$dt("QRILC"),
                    tags$dd("MNAR method for left-censored proteomics data. It imputes low missing values from the estimated left tail of the intensity distribution."),
                    tags$dt("MinProb"),
                    tags$dd("MNAR method that imputes low values around a low quantile. It is simpler than QRILC and often stable."),
                    tags$dt("kNN"),
                    tags$dd("MAR method. Missing values are estimated from proteins with similar abundance profiles. Sparse rows are routed to the selected MNAR method."),
                    tags$dt("Mixed"),
                    tags$dd("Classifies missingness first. MNAR cases use the selected MNAR method; MAR cases use kNN."),
                    tags$dt("Per group"),
                    tags$dd("Runs classification or imputation within each condition group. Useful when missingness is condition-specific.")
                  )
                )
              ),

              hr(),
              actionButton(ns("btn_impute"), "Run Imputation",
                           class = "btn-primary btn-block"),
              br(),
              uiOutput(ns("impute_status"))
            )
          ),
          column(8,
            tabsetPanel(
              tabPanel("Before / After Density",
                br(),
                plotOutput(ns("density_plot"), height = "400px")
              ),
              tabPanel("Missing Value Summary",
                br(),
                plotOutput(ns("mv_barplot"), height = "400px"),
                br(),
                DT::dataTableOutput(ns("mv_table"))
              ),
              tabPanel("Group Statistics",
                br(),
                tags$div(class = "alert alert-info",
                  "Per-protein, per-group summary statistics from the raw data (before imputation). ",
                  "Use this to spot low-abundance proteins with scattered missing values, a single ",
                  "unusually high replicate, or proteins missing from a whole group."),
                tags$div(class = "alert alert-warning",
                  tags$b("SD/CV% from few values are unreliable: "),
                  "red = computed from only 2 valid values (SD is systematically biased low, roughly ±75% uncertainty), ",
                  "orange = from 3 values (±50% uncertainty). No highlight = 4 or more values."),
                DT::dataTableOutput(ns("group_stats_table")),
                br(),
                downloadButton(ns("dl_group_stats"), "Download Group Statistics CSV",
                               class = "btn-default btn-sm")
              ),
              tabPanel("Threshold Sensitivity",
                br(),
                tags$div(class = "alert alert-info",
                  "Shows how many partially observed protein-by-group cases switch between kNN and MNAR as the mixed threshold changes."),
                plotOutput(ns("threshold_sensitivity_plot"), height = "350px"),
                uiOutput(ns("threshold_sensitivity_counts")),
                br(),
                DT::dataTableOutput(ns("threshold_sensitivity_table")),
                br(),
                downloadButton(ns("dl_threshold_sensitivity"), "Download Threshold Sensitivity CSV",
                               class = "btn-default btn-sm")
              ),
              tabPanel("Imputed Data Preview",
                br(),
                DT::dataTableOutput(ns("imputed_preview"))
              )
            )
          )
        ),
        fluidRow(
          column(12,
            div(class = "app-action-row",
              downloadButton(ns("dl_imputed"), "Download Imputed CSV",
                             class = "btn-default"),
              div(class = "app-proceed-wrap",
                uiOutput(ns("proceed_btn"))
              )
            )
          )
        )
      )
    })
    observeEvent(input$btn_impute, {
      req(isTRUE(app_state$upload_done))

      int_raw  <- app_state$parsed_data$intensity
      meta     <- app_state$parsed_data$meta

      withProgress(message = "Running imputation...", value = 0, {

        incProgress(0.1, detail = "Pre-filtering proteins")
        # Drop proteins with fewer than min_valid valid values. This removes
        # them from the whole pipeline, not just this step - always report
        # exactly how many and why, never drop rows silently.
        n_valid  <- rowSums(!is.na(int_raw))
        min_valid <- max(1, as.integer(input$min_valid))
        n_dropped <- sum(n_valid < min_valid)
        if (n_dropped > 0) {
          showNotification(
            sprintf("%d protein(s) had fewer than %d valid value(s) across all samples and were excluded before imputation (removed from the whole pipeline).",
                    n_dropped, min_valid),
            type = "warning", duration = 10)
        }
        int_filt <- int_raw[n_valid >= min_valid, , drop = FALSE]
        meta_filt <- meta[rownames(int_filt), , drop = FALSE]

        incProgress(0.2, detail = paste0("Applying method: ", input$method))
        seed <- as.integer(input$random_seed)
        if (is.na(seed)) seed <- 12345L
        set.seed(seed)

        int_imp <- tryCatch({
          switch(input$method,
            mnar_min       = impute_mnar_lcmd(int_filt,
                                              method = input$mnar_backend,
                                              q = input$mnar_q,
                                              tune.sigma = input$mnar_tune_sigma,
                                              perseus_width = input$mnar_perseus_width,
                                              perseus_downshift = input$mnar_perseus_downshift),
            mnar_group_min = impute_mnar_group_min(int_filt,
                                                   app_state$parsed_data$sample_names,
                                                   method = input$mnar_backend,
                                                   q = input$mnar_q,
                                                   tune.sigma = input$mnar_tune_sigma,
                                                   perseus_width = input$mnar_perseus_width,
                                                   perseus_downshift = input$mnar_perseus_downshift),
            mar_knn        = if (isTRUE(input$per_group_knn)) {
                                impute_mar_knn_by_group(int_filt,
                                                        app_state$parsed_data$sample_names,
                                                        k = input$knn_k,
                                                        mnar_method = input$mnar_backend,
                                                        mnar_q = input$mnar_q,
                                                        mnar_tune_sigma = input$mnar_tune_sigma,
                                                        perseus_width = input$mnar_perseus_width,
                                                        perseus_downshift = input$mnar_perseus_downshift)
                              } else {
                                safe_knn_impute(int_filt,
                                                k = input$knn_k,
                                                mnar_method = input$mnar_backend,
                                                mnar_q = input$mnar_q,
                                                mnar_tune_sigma = input$mnar_tune_sigma,
                                                perseus_width = input$mnar_perseus_width,
                                                perseus_downshift = input$mnar_perseus_downshift)
                              },
            mixed          = if (isTRUE(input$per_group_knn)) {
                                impute_mixed_by_group(int_filt,
                                                      app_state$parsed_data$sample_names,
                                                      k = input$knn_k,
                                                      mnar_method = input$mnar_backend,
                                                      mnar_q = input$mnar_q,
                                                      mnar_tune_sigma = input$mnar_tune_sigma,
                                                      perseus_width = input$mnar_perseus_width,
                                                      perseus_downshift = input$mnar_perseus_downshift,
                                                      mnar_thr = input$mnar_threshold,
                                                      knn_scope = input$mixed_group_strategy)
                              } else {
                                impute_mixed(int_filt,
                                             k = input$knn_k,
                                             mnar_method = input$mnar_backend,
                                             mnar_q = input$mnar_q,
                                             mnar_tune_sigma = input$mnar_tune_sigma,
                                             perseus_width = input$mnar_perseus_width,
                                             perseus_downshift = input$mnar_perseus_downshift,
                                             mnar_thr = input$mnar_threshold)
                              }
          )
        }, error = function(e) {
          showNotification(paste("Imputation error:", conditionMessage(e)),
                           type = "error", duration = 12)
          NULL
        })

        incProgress(0.6, detail = "Storing results")
        if (!is.null(int_imp)) {
          result <- list(
            intensity     = int_imp,
            intensity_raw = int_filt,
            meta          = meta_filt,
            method        = input$method,
            mnar_backend  = if (input$method %in% c("mnar_min", "mnar_group_min", "mixed")) input$mnar_backend else NA_character_,
            perseus_width = if (identical(input$mnar_backend, "Perseus")) input$mnar_perseus_width else NA_real_,
            perseus_downshift = if (identical(input$mnar_backend, "Perseus")) input$mnar_perseus_downshift else NA_real_,
            seed          = seed,
            min_valid     = min_valid,
            rows_dropped  = n_dropped
          )
          imputed_data(result)
          app_state$imputed_data  <- result
          app_state$imputation_done <- TRUE
        }
        incProgress(0.1, detail = "Done")
      })
    })

    observeEvent(input$skip_imputation, {
      if (isTRUE(input$skip_imputation) && isTRUE(app_state$upload_done)) {
        p <- app_state$parsed_data
        result <- list(
          intensity     = p$intensity,
          intensity_raw = p$intensity,
          meta          = p$meta,
          method        = "none",
          mnar_backend  = NA_character_,
          seed          = NA_integer_,
          min_valid     = NA_integer_,
          rows_dropped  = 0L
        )
        imputed_data(result)
        app_state$imputed_data    <- result
        app_state$imputation_done <- TRUE
        showNotification("Imputation skipped - using data as-is.", type = "message")
      }
    })

    output$prefilter_note <- renderUI({
      req(isTRUE(app_state$upload_done))
      int <- app_state$parsed_data$intensity
      min_valid <- max(1, as.integer(if (is.null(input$min_valid)) 1 else input$min_valid))
      n_drop <- sum(rowSums(!is.na(int)) < min_valid)
      if (n_drop > 0) {
        tags$div(class = "alert alert-warning", style = "padding: 6px 10px; margin-bottom: 8px;",
          sprintf("%d protein row(s) will be removed at the current pre-filter setting.", n_drop))
      } else {
        tags$div(class = "alert alert-info", style = "padding: 6px 10px; margin-bottom: 8px;",
          "No protein rows will be removed at the current pre-filter setting.")
      }
    })

    output$impute_status <- renderUI({
      d <- imputed_data()
      if (is.null(d)) return(NULL)
      seed_txt <- if (!is.null(d$seed) && !is.na(d$seed)) {
        sprintf(" Random seed: %d.", d$seed)
      } else {
        ""
      }
      drop_txt <- if (!is.null(d$rows_dropped) && d$rows_dropped > 0) {
        sprintf(" Removed %d protein row(s) by pre-filter.", d$rows_dropped)
      } else {
        " No protein rows removed by pre-filter."
      }
      tags$div(class = "alert alert-success",
               sprintf("Done. %d proteins x %d samples.%s%s",
                       nrow(d$intensity), ncol(d$intensity), drop_txt, seed_txt))
    })
    # Missing_frac only takes a handful of discrete values per group (e.g. with
    # 4 replicates: 0.25/0.5/0.75/1.0), so moving the threshold slider within a
    # gap between two of those values changes nothing. This spells out, for
    # the actual group sizes in the current data, exactly how many replicates
    # need to be missing before MNAR kicks in - computed from the current
    # upload, not a hardcoded replicate count.
    output$mnar_threshold_hint <- renderUI({
      req(isTRUE(app_state$upload_done))
      thr   <- if (is.null(input$mnar_threshold)) 0.5 else input$mnar_threshold
      grps  <- infer_groups(app_state$parsed_data$sample_names)
      sizes <- sort(unique(table(grps)))
      lines <- vapply(sizes, function(n) {
        need <- ceiling(thr * n - 1e-9)
        sprintf("groups of %d: MNAR triggers at %d+ missing (%d/%d)", n, need, need, n)
      }, character(1))
      helpText(paste(lines, collapse = "; "))
    })
    # Density plot before / after
    output$density_plot <- renderPlot({
      req(isTRUE(app_state$upload_done))
      d_imp <- imputed_data()

      int_before <- app_state$parsed_data$intensity
      vals_before <- as.vector(as.matrix(int_before))
      vals_before <- vals_before[!is.na(vals_before)]

      par(mar = c(4, 4, 3, 1))
      d_b <- density(vals_before, na.rm = TRUE)
      plot(d_b, col = "steelblue", lwd = 2,
           main = "Intensity Distribution Before/After Imputation",
           xlab = "log2 Intensity", ylab = "Density", ylim = NULL)

      if (!is.null(d_imp)) {
        vals_after <- as.vector(as.matrix(d_imp$intensity))
        vals_after <- vals_after[!is.na(vals_after)]
        d_a <- density(vals_after, na.rm = TRUE)
        lines(d_a, col = "tomato", lwd = 2, lty = 2)
        legend("topright", legend = c("Before imputation", "After imputation"),
               col = c("steelblue", "tomato"), lwd = 2, lty = c(1, 2), bty = "n")
      } else {
        legend("topright", legend = "Before imputation (no imputation run yet)",
               col = "steelblue", lwd = 2, bty = "n")
      }
    })
    # Missing value barplot
    output$mv_barplot <- renderPlot({
      req(isTRUE(app_state$upload_done))
      int <- app_state$parsed_data$intensity
      grps <- infer_groups(colnames(int))
      uniq_grps <- sort(unique(grps))
      mv_pct <- vapply(uniq_grps, function(g) {
        cols <- colnames(int)[grps == g]
        round(100 * mean(is.na(int[, cols, drop = FALSE])), 1)
      }, numeric(1))
      par(mar = c(6, 4, 3, 1))
      bp <- barplot(mv_pct, names.arg = uniq_grps, las = 2, col = "steelblue",
                    main = "Missing Value % per Group", ylab = "Missing (%)",
                    ylim = c(0, max(mv_pct) * 1.2))
      text(bp, mv_pct + 0.5, paste0(mv_pct, "%"), cex = 0.85, font = 2)
    })

    output$mv_table <- DT::renderDataTable({
      req(isTRUE(app_state$upload_done))
      int  <- app_state$parsed_data$intensity
      grps <- infer_groups(colnames(int))
      uniq_grps <- sort(unique(grps))
      rows <- lapply(uniq_grps, function(g) {
        cols <- colnames(int)[grps == g]
        sub  <- int[, cols, drop = FALSE]
        data.frame(
          Group           = g,
          Samples         = length(cols),
          Total_Values    = nrow(sub) * ncol(sub),
          Missing         = sum(is.na(sub)),
          Missing_Pct     = round(100 * mean(is.na(sub)), 2),
          stringsAsFactors = FALSE
        )
      })
      df <- do.call(rbind, rows)
      DT::datatable(df, rownames = FALSE, options = list(dom = "t"))
    })
    # Group statistics - per-protein, per-group descriptive stats on raw data
    group_stats_display_df <- function() {
      req(isTRUE(app_state$upload_done))
      wide <- compute_group_stats_table(app_state$parsed_data$intensity,
                                        app_state$parsed_data$sample_names)
      meta <- app_state$parsed_data$meta
      cbind(meta[rownames(wide), , drop = FALSE], wide)
    }

    output$group_stats_table <- DT::renderDataTable({
      df <- group_stats_display_df()
      dt <- DT::datatable(df, rownames = TRUE,
                          options = list(pageLength = 15, scrollX = TRUE, dom = "ftip"))
      dt <- DT::formatRound(dt, grep("_(Mean|Median|SD|Min|Max)$", names(df), value = TRUE), 3)
      dt <- DT::formatRound(dt, grep("_(Missing_Pct|CV_pct)$", names(df), value = TRUE), 2)

      # Flag SD/CV% cells by how many valid replicates they came from - a
      # sample SD from very few points isn't just imprecise, it's
      # systematically biased low (verified empirically: ~20% low at n=2,
      # ~11% low at n=3, vs <1% at n=32).
      uniq_grps <- unique(sub("_N_valid$", "", grep("_N_valid$", names(df), value = TRUE)))
      for (g in uniq_grps) {
        nvalid_col   <- paste0(g, "_N_valid")
        flagged_cols <- intersect(paste0(g, c("_N_valid", "_SD", "_CV_pct")), names(df))
        dt <- DT::formatStyle(
          dt, flagged_cols, valueColumns = nvalid_col,
          backgroundColor = DT::styleInterval(c(2, 3), c("#f2dede", "#fcf8e3", ""))
        )
      }
      dt
    })

    output$dl_group_stats <- downloadHandler(
      filename = function() download_filename(app_state, "group_statistics", "csv"),
      content  = function(file) {
        utils::write.csv(group_stats_display_df(), file, row.names = TRUE)
      }
    )
    threshold_sensitivity_df <- reactive({
      req(isTRUE(app_state$upload_done))
      compute_threshold_sensitivity(app_state$parsed_data$intensity,
                                    app_state$parsed_data$sample_names)
    })

    output$threshold_sensitivity_plot <- renderPlot({
      df <- threshold_sensitivity_df()
      req(nrow(df) > 0)
      thr <- if (is.null(input$mnar_threshold)) 0.5 else input$mnar_threshold
      par(mar = c(4.5, 4.5, 3, 1))
      plot(df$Threshold, df$MNAR_cases, type = "b", pch = 19, lwd = 2,
           col = "steelblue", ylim = c(0, max(df$Partial_cases, na.rm = TRUE)),
           xlab = "MNAR threshold", ylab = "Protein-by-group cases",
           main = "Mixed threshold sensitivity")
      lines(df$Threshold, df$KNN_cases, type = "b", pch = 19, lwd = 2, col = "tomato")
      abline(v = thr, lty = 2, lwd = 2, col = "gray40")
      legend("right", bty = "n", lwd = 2, lty = c(1, 1, 2), pch = c(19, 19, NA),
             col = c("steelblue", "tomato", "gray40"),
             legend = c("MNAR downshift", "kNN", "Current threshold"))
    })

    output$threshold_sensitivity_counts <- renderUI({
      df <- threshold_sensitivity_df()
      req(nrow(df) > 0)
      thr <- if (is.null(input$mnar_threshold)) 0.5 else input$mnar_threshold
      idx <- which.min(abs(df$Threshold - thr))
      tags$div(style = "text-align:center; font-weight:600; margin: 4px 0 10px 0;",
        sprintf("Current threshold %.2f: %d MNAR, %d kNN, %d partial cases.",
                thr, df$MNAR_cases[idx], df$KNN_cases[idx], df$Partial_cases[idx]))
    })

    output$threshold_sensitivity_table <- DT::renderDataTable({
      df <- threshold_sensitivity_df()
      DT::datatable(df, rownames = FALSE,
                    options = list(pageLength = 10, scrollX = TRUE, dom = "tip"))
    })

    output$dl_threshold_sensitivity <- downloadHandler(
      filename = function() download_filename(app_state, "threshold_sensitivity", "csv"),
      content  = function(file) {
        utils::write.csv(threshold_sensitivity_df(), file, row.names = FALSE)
      }
    )
    # Imputed data preview
    output$imputed_preview <- DT::renderDataTable({
      req(imputed_data())
      d   <- imputed_data()
      df  <- cbind(d$meta, round(d$intensity, 3))
      DT::datatable(df, rownames = FALSE,
                    options = list(pageLength = 10, scrollX = TRUE, dom = "ftip"))
    })
    output$dl_imputed <- downloadHandler(
      filename = function() download_filename(app_state, "imputed_data", "csv"),
      content  = function(file) {
        req(imputed_data())
        d  <- imputed_data()
        df <- cbind(d$meta, d$intensity)
        utils::write.csv(df, file, row.names = TRUE)
      }
    )
    output$proceed_btn <- renderUI({
      if (!isTRUE(app_state$imputation_done)) return(NULL)
      actionButton(ns("btn_proceed"), "Proceed to Exploration",
                   class = "btn-success btn-lg app-proceed-btn")
    })

    observeEvent(input$btn_proceed, {
      app_state$active_tab <- "exploration"
    })

  })
}
# Imputation routines
impute_mnar_lcmd <- function(mat, method = "Perseus", q = 0.01, tune.sigma = 1,
                             perseus_width = 0.3, perseus_downshift = 1.8) {
  mat <- as.data.frame(mat)
  if (!anyNA(mat)) return(mat)

  method <- match.arg(method, c("Perseus", "QRILC", "MinProb"))
  q <- if (is.null(q) || is.na(q)) 0.01 else as.numeric(q)
  tune.sigma <- if (is.null(tune.sigma) || is.na(tune.sigma)) 1 else as.numeric(tune.sigma)
  mat_in <- as.matrix(mat)

  if (method == "Perseus") {
    return(impute_mnar_perseus(mat_in, width = perseus_width, downshift = perseus_downshift))
  }

  if (!requireNamespace("imputeLCMD", quietly = TRUE)) {
    stop(
      paste(
        "Package 'imputeLCMD' is required for QRILC or MinProb MNAR imputation.",
        "Run install_packages.bat from the app folder, then restart the app.",
        paste0("Current working directory: ", getwd()),
        paste0("R library paths: ", paste(.libPaths(), collapse = " | ")),
        sep = "\n"
      )
    )
  }

  res <- NULL
  utils::capture.output({
    res <- if (method == "QRILC") {
      qr <- suppressMessages(imputeLCMD::impute.QRILC(mat_in, tune.sigma = tune.sigma))
      if (is.list(qr)) qr[[1]] else qr
    } else {
      suppressMessages(imputeLCMD::impute.MinProb(mat_in, q = q, tune.sigma = tune.sigma))
    }
  })

  res <- as.data.frame(res)
  names(res) <- names(mat)
  rownames(res) <- rownames(mat)
  res
}

impute_mnar_perseus <- function(mat, width = 0.3, downshift = 1.8) {
  width <- if (is.null(width) || is.na(width)) 0.3 else as.numeric(width)
  downshift <- if (is.null(downshift) || is.na(downshift)) 1.8 else as.numeric(downshift)
  mat <- as.matrix(mat)
  out <- mat
  for (j in seq_len(ncol(out))) {
    observed <- out[, j][!is.na(out[, j])]
    if (!length(observed)) next
    mu <- mean(observed, na.rm = TRUE) - downshift * stats::sd(observed, na.rm = TRUE)
    sigma <- width * stats::sd(observed, na.rm = TRUE)
    if (!is.finite(sigma) || sigma <= 0) sigma <- 0.01
    miss <- which(is.na(out[, j]))
    if (length(miss)) out[miss, j] <- stats::rnorm(length(miss), mean = mu, sd = sigma)
  }
  as.data.frame(out)
}

impute_mnar_group_min <- function(mat, sample_names, method = "Perseus", q = 0.01, tune.sigma = 1,
                                  perseus_width = 0.3, perseus_downshift = 1.8) {
  grps <- infer_groups(sample_names)
  mat  <- as.data.frame(mat)
  for (g in unique(grps)) {
    cols <- sample_names[grps == g]
    cols <- intersect(cols, colnames(mat))
    if (length(cols) && anyNA(mat[, cols, drop = FALSE])) {
      filled <- impute_mnar_lcmd(mat[, cols, drop = FALSE],
                                 method = method, q = q, tune.sigma = tune.sigma,
                                 perseus_width = perseus_width,
                                 perseus_downshift = perseus_downshift)
      block <- mat[, cols, drop = FALSE]
      na_mask <- is.na(block)
      block[na_mask] <- filled[na_mask]
      mat[, cols] <- block
    }
  }
  as.data.frame(mat)
}

impute_mar_knn <- function(mat, k = 5) {
  if (!requireNamespace("impute", quietly = TRUE)) {
    stop("Package 'impute' is required for kNN imputation. Install it before running MAR or mixed kNN methods.")
  }
  res <- NULL
  utils::capture.output({
    res <- impute::impute.knn(as.matrix(mat), k = k)
  })
  as.data.frame(res$data)
}

# Standard entry point for kNN imputation everywhere in this file. Rows that
# would trigger impute.knn's high-missingness fallback are pulled out first and
# filled with the selected MNAR method instead.
safe_knn_impute <- function(mat, k = 5, mnar_method = "Perseus", mnar_q = 0.01,
                            mnar_tune_sigma = 1, rowmax = 0.5,
                            perseus_width = 0.3, perseus_downshift = 1.8) {
  mat       <- as.data.frame(mat)
  miss_frac <- rowMeans(is.na(mat))
  risky     <- rownames(mat)[miss_frac >= rowmax]
  safe      <- setdiff(rownames(mat), risky)

  if (length(risky)) {
    showNotification(
      sprintf("%d protein(s) had %.0f%%+ missing values among the columns being imputed; used %s instead of kNN for these.",
              length(risky), rowmax * 100, mnar_method),
      type = "warning", duration = 10)
    mnar_filled  <- impute_mnar_lcmd(mat, method = mnar_method,
                                     q = mnar_q, tune.sigma = mnar_tune_sigma,
                                     perseus_width = perseus_width,
                                     perseus_downshift = perseus_downshift)
    mat[risky, ] <- mnar_filled[risky, , drop = FALSE]
  }
  if (length(safe)) {
    sub <- mat[safe, , drop = FALSE]
    if (anyNA(sub)) mat[safe, ] <- impute_mar_knn(sub, k = k)
  }
  as.data.frame(mat)
}

# Runs kNN separately within each group's own columns instead of pooled.
impute_mar_knn_by_group <- function(mat, sample_names, k = 5, mnar_method = "Perseus",
                                    mnar_q = 0.01, mnar_tune_sigma = 1,
                                    perseus_width = 0.3, perseus_downshift = 1.8) {
  mat  <- as.data.frame(mat)
  grps <- infer_groups(sample_names)
  for (g in unique(grps)) {
    cols <- intersect(sample_names[grps == g], colnames(mat))
    if (length(cols) && anyNA(mat[, cols, drop = FALSE])) {
      mat[, cols] <- safe_knn_impute(mat[, cols, drop = FALSE], k = k,
                                     mnar_method = mnar_method, mnar_q = mnar_q,
                                     mnar_tune_sigma = mnar_tune_sigma,
                                     perseus_width = perseus_width,
                                     perseus_downshift = perseus_downshift)
    }
  }
  as.data.frame(mat)
}

impute_mixed <- function(mat, k = 5, mnar_method = "Perseus", mnar_q = 0.01,
                         mnar_tune_sigma = 1, perseus_width = 0.3,
                         perseus_downshift = 1.8, mnar_thr = 0.5) {
  mat_orig  <- as.data.frame(mat)
  miss_frac <- rowMeans(is.na(mat_orig))
  mnar_rows <- rownames(mat_orig)[miss_frac >= mnar_thr]
  mar_rows  <- rownames(mat_orig)[miss_frac > 0 & miss_frac < mnar_thr]
  out       <- mat_orig

  if (length(mnar_rows)) {
    mnar_filled <- impute_mnar_lcmd(mat_orig, method = mnar_method,
                                    q = mnar_q, tune.sigma = mnar_tune_sigma,
                                    perseus_width = perseus_width,
                                    perseus_downshift = perseus_downshift)
    block <- out[mnar_rows, , drop = FALSE]
    na_mask <- is.na(block)
    block[na_mask] <- mnar_filled[mnar_rows, , drop = FALSE][na_mask]
    out[mnar_rows, ] <- block
  }

  if (length(mar_rows)) {
    mar_filled <- safe_knn_impute(mat_orig, k = k, mnar_method = mnar_method,
                                  mnar_q = mnar_q, mnar_tune_sigma = mnar_tune_sigma,
                                  perseus_width = perseus_width,
                                  perseus_downshift = perseus_downshift)
    block <- out[mar_rows, , drop = FALSE]
    na_mask <- is.na(block)
    block[na_mask] <- mar_filled[mar_rows, , drop = FALSE][na_mask]
    out[mar_rows, ] <- block
  }

  as.data.frame(out)
}

# Per-protein, per-group MNAR/MAR/complete classification (vs. impute_mixed's
# classification, which pools missingness across all groups at once).
classify_missingness_by_group <- function(mat, sample_names, mnar_thr = 0.5) {
  mat  <- as.data.frame(mat)
  grps <- infer_groups(sample_names)
  uniq_grps <- sort(unique(grps))
  cls <- matrix(NA_character_, nrow(mat), length(uniq_grps),
                dimnames = list(rownames(mat), uniq_grps))
  for (g in uniq_grps) {
    cols <- intersect(sample_names[grps == g], colnames(mat))
    if (!length(cols)) next
    miss_frac <- rowMeans(is.na(mat[, cols, drop = FALSE]))
    cls[, g]  <- ifelse(miss_frac <= 0, "complete",
                  ifelse(miss_frac >= mnar_thr, "mnar", "mar"))
  }
  cls
}

# Both "Mixed, per group" strategies live in one function - they only differ
# in knn_scope (whether kNN neighbour search pools all groups or is restricted
# to each group's own columns). MNAR handling is identical for both.
impute_mixed_by_group <- function(mat, sample_names, k = 5, mnar_method = "Perseus",
                                   mnar_q = 0.01, mnar_tune_sigma = 1,
                                   perseus_width = 0.3, perseus_downshift = 1.8, mnar_thr = 0.5,
                                   knn_scope = c("pooled", "per_group")) {
  knn_scope <- match.arg(knn_scope)
  mat_orig  <- as.data.frame(mat)
  grps      <- infer_groups(sample_names)
  uniq_grps <- sort(unique(grps))
  cls <- classify_missingness_by_group(mat_orig, sample_names, mnar_thr = mnar_thr)

  # Computed once from the untouched matrix, so kNN's neighbour search is
  # never contaminated by MNAR-fabricated values sitting in the same row.
  mar_filled <- if (knn_scope == "pooled") {
    safe_knn_impute(mat_orig, k = k, mnar_method = mnar_method,
                    mnar_q = mnar_q, mnar_tune_sigma = mnar_tune_sigma,
                    perseus_width = perseus_width,
                    perseus_downshift = perseus_downshift)
  } else {
    impute_mar_knn_by_group(mat_orig, sample_names, k = k,
                            mnar_method = mnar_method, mnar_q = mnar_q,
                            mnar_tune_sigma = mnar_tune_sigma,
                            perseus_width = perseus_width,
                            perseus_downshift = perseus_downshift)
  }

  mnar_filled_by_group <- NULL

  out <- mat_orig
  for (g in uniq_grps) {
    cols <- intersect(sample_names[grps == g], colnames(mat_orig))
    if (!length(cols)) next
    mnar_rows <- rownames(cls)[cls[, g] == "mnar"]
    mar_rows  <- rownames(cls)[cls[, g] == "mar"]

    if (length(mnar_rows)) {
      if (is.null(mnar_filled_by_group)) {
        mnar_filled_by_group <- impute_mnar_group_min(mat_orig, sample_names,
                                                      method = mnar_method,
                                                      q = mnar_q,
                                                      tune.sigma = mnar_tune_sigma,
                                                      perseus_width = perseus_width,
                                                      perseus_downshift = perseus_downshift)
      }
      block   <- out[mnar_rows, cols, drop = FALSE]
      na_mask <- is.na(block)
      if (any(na_mask)) {
        filled <- mnar_filled_by_group[mnar_rows, cols, drop = FALSE]
        block[na_mask] <- filled[na_mask]
        out[mnar_rows, cols] <- block
      }
    }
    if (length(mar_rows)) {
      block   <- out[mar_rows, cols, drop = FALSE]
      na_mask <- is.na(block)
      if (any(na_mask)) {
        filled <- mar_filled[mar_rows, cols, drop = FALSE]
        block[na_mask] <- filled[na_mask]
        out[mar_rows, cols] <- block
      }
    }
  }
  as.data.frame(out)
}

# One row per protein, columns repeated per group as "<Group>_<Stat>".
# Computed on RAW (pre-imputation) data - same source as mv_table/mv_barplot.
compute_group_stats_table <- function(mat, sample_names) {
  mat  <- as.data.frame(mat)
  grps <- infer_groups(sample_names)
  uniq_grps <- sort(unique(grps))

  blocks <- lapply(uniq_grps, function(g) {
    cols <- intersect(sample_names[grps == g], colnames(mat))
    sub  <- as.matrix(mat[, cols, drop = FALSE])

    n_total <- ncol(sub)
    n_valid <- rowSums(!is.na(sub))
    all_na  <- n_valid == 0

    mean_   <- rowMeans(sub, na.rm = TRUE); mean_[all_na] <- NA_real_
    median_ <- apply(sub, 1, function(x) if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE))
    sd_     <- apply(sub, 1, function(x) if (sum(!is.na(x)) < 2) NA_real_ else sd(x, na.rm = TRUE))
    min_    <- apply(sub, 1, function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE))
    max_    <- apply(sub, 1, function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))
    cv_pct  <- ifelse(!is.na(sd_) & !is.na(mean_) & mean_ != 0, 100 * sd_ / mean_, NA_real_)

    out <- data.frame(
      N_valid = n_valid, N_total = n_total,
      Missing_Pct = 100 * (1 - n_valid / n_total),
      Mean = mean_, Median = median_, SD = sd_, CV_pct = cv_pct,
      Min = min_, Max = max_,
      row.names = rownames(mat)
    )
    names(out) <- paste0(g, "_", names(out))
    out
  })
  do.call(cbind, blocks)
}

compute_threshold_sensitivity <- function(mat, sample_names, thresholds = seq(0.1, 0.9, by = 0.05)) {
  mat <- as.data.frame(mat)
  grps <- infer_groups(sample_names)
  uniq_grps <- sort(unique(grps))

  partial_missing <- numeric(0)
  groups_seen <- character(0)
  for (g in uniq_grps) {
    cols <- intersect(sample_names[grps == g], colnames(mat))
    if (!length(cols)) next
    sub <- mat[, cols, drop = FALSE]
    n_total <- ncol(sub)
    n_valid <- rowSums(!is.na(sub))
    partial <- n_valid > 0 & n_valid < n_total
    partial_missing <- c(partial_missing, 1 - n_valid[partial] / n_total)
    groups_seen <- c(groups_seen, rep(g, sum(partial)))
  }

  if (!length(partial_missing)) {
    return(data.frame(
      Threshold = thresholds,
      Partial_cases = 0L,
      MNAR_cases = 0L,
      KNN_cases = 0L,
      MNAR_pct = 0,
      KNN_pct = 0,
      Groups_with_partial_cases = 0L
    ))
  }

  rows <- lapply(thresholds, function(thr) {
    n_mnar <- sum(partial_missing >= thr)
    n_total <- length(partial_missing)
    data.frame(
      Threshold = thr,
      Partial_cases = n_total,
      MNAR_cases = n_mnar,
      KNN_cases = n_total - n_mnar,
      MNAR_pct = round(100 * n_mnar / n_total, 1),
      KNN_pct = round(100 * (n_total - n_mnar) / n_total, 1),
      Groups_with_partial_cases = length(unique(groups_seen)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
