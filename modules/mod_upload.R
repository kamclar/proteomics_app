# Upload/import 


mod_upload_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Step 1. Upload & Configure Data"),
    fluidRow(
      column(4,
        wellPanel(
          h4("Upload"),
          fileInput(ns("file"), "Select Spectronaut XLSX file",
                    accept = c(".xlsx", ".xls")),
          uiOutput(ns("sheet_picker")),
          numericInput(ns("header_row"), "Column names are in row", value = 1,
                       min = 1, max = 2, step = 1),
          helpText("Use row 2 when row 1 contains merged group headers."),
          hr(),
          h4("Options"),
          checkboxInput(ns("already_log2"), "Intensities already log2-transformed", value = TRUE),
          helpText("If unchecked, log2 will be applied during import."),
          hr(),
          h4("Columns"),
          uiOutput(ns("column_picker")),
          hr(),
          actionButton(ns("btn_parse"), "Parse File", class = "btn-success btn-lg btn-block app-proceed-btn"),
          br(),
          uiOutput(ns("parse_status"))
        )
      ),
      column(8,
        tabsetPanel(
          tabPanel("Sheet Preview",
            br(),
            uiOutput(ns("preview_summary")),
            br(),
            DT::dataTableOutput(ns("sheet_preview"))
          ),
          tabPanel("Sample Overview",
            br(),
            uiOutput(ns("sample_summary")),
            br(),
            DT::dataTableOutput(ns("sample_table"))
          ),
          tabPanel("Protein Preview",
            br(),
            uiOutput(ns("protein_summary")),
            br(),
            DT::dataTableOutput(ns("meta_preview"))
          ),
          tabPanel("T-test Comparisons Detected",
            br(),
            uiOutput(ns("ttest_summary")),
            br(),
            DT::dataTableOutput(ns("ttest_preview"))
          ),
          tabPanel("Missing Value Map",
            br(),
            plotOutput(ns("mv_heatmap"), height = "420px"),
            br(),
            uiOutput(ns("mv_stats"))
          )
        )
      )
    ),
    fluidRow(
      column(12,
        div(class = "text-right",
          uiOutput(ns("proceed_btn"))
        )
      )
    )
  )
}

mod_upload_server <- function(id, app_state) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    parsed <- reactiveVal(NULL)
    empty_chr <- function(x) if (is.null(x)) character(0) else x
    output$sheet_picker <- renderUI({
      req(input$file)
      sheets <- openxlsx::getSheetNames(input$file$datapath)
      selectInput(ns("sheet"), "Sheet", choices = sheets, selected = sheets[1])
    })

    observeEvent({
      input$file
      input$sheet
    }, {
      req(input$file, input$sheet)
      sheets <- openxlsx::getSheetNames(input$file$datapath)
      sheet_idx <- match(input$sheet, sheets)
      if (is.na(sheet_idx)) sheet_idx <- 1L

      guessed <- guess_header_row(input$file$datapath, sheet = sheet_idx)
      updateNumericInput(session, "header_row", value = guessed)
    }, ignoreInit = FALSE)

    selected_sheet_index <- reactive({
      req(input$file)
      sheets <- openxlsx::getSheetNames(input$file$datapath)
      sheet_name <- input$sheet
      if (is.null(sheet_name) || !(sheet_name %in% sheets)) return(1L)
      match(sheet_name, sheets)
    })

    raw_for_selection <- reactive({
      req(input$file, input$header_row)
      withProgress(message = "Loading workbook preview...", value = 0, style = "notification", {
        incProgress(0.35, detail = "Reading selected sheet")
        raw <- read_proteomics_sheet(
          path = input$file$datapath,
          sheet = selected_sheet_index(),
          header_row = input$header_row
        )
        incProgress(0.65, detail = "Preparing column detection")
        raw
      })
    })

    output$column_picker <- renderUI({
      req(raw_for_selection())
      raw <- raw_for_selection()
      detected <- detect_upload_columns(raw)
      choices <- detected$all_columns

      tagList(
        selectizeInput(
          ns("metadata_cols"), "Metadata / annotation columns",
          choices = choices, selected = detected$metadata, multiple = TRUE,
          options = list(plugins = list("remove_button"))
        ),
        selectizeInput(
          ns("intensity_cols"), "Intensity columns",
          choices = choices, selected = detected$intensity, multiple = TRUE,
          options = list(plugins = list("remove_button"))
        ),
        selectizeInput(
          ns("ttest_cols"), "T-test result columns (optional)",
          choices = choices, selected = detected$ttest, multiple = TRUE,
          options = list(plugins = list("remove_button"))
        )
      )
    })

    output$preview_summary <- renderUI({
      req(raw_for_selection())
      raw <- raw_for_selection()
      detected <- detect_upload_columns(raw)
      tagList(
        tags$p(strong("Columns found: "), ncol(raw)),
        tags$p(
          strong("Auto-detected: "),
          length(detected$metadata), " metadata, ",
          length(detected$intensity), " intensity, ",
          length(detected$ttest), " t-test columns"
        )
      )
    })

    output$sheet_preview <- DT::renderDataTable({
      req(input$file)
      preview <- read_sheet_preview(
        path = input$file$datapath,
        sheet = selected_sheet_index(),
        n_rows = 25
      )
      DT::datatable(
        preview,
        rownames = TRUE,
        options = list(pageLength = 10, scrollX = TRUE, dom = "tip")
      )
    })
    observeEvent(input$btn_parse, {
      req(input$file, input$header_row, input$intensity_cols)

      withProgress(message = "Parsing XLSX...", value = 0, style = "notification", {
        incProgress(0.2, detail = "Reading file")
        tryCatch({
          result <- parse_proteomics_xlsx(
            path         = input$file$datapath,
            sheet        = selected_sheet_index(),
            already_log2 = input$already_log2,
            header_row   = input$header_row,
            source_name  = input$file$name,
            selected_columns = list(
              metadata = empty_chr(input$metadata_cols),
              intensity = empty_chr(input$intensity_cols),
              ttest = empty_chr(input$ttest_cols)
            )
          )
          incProgress(0.6, detail = "Processing columns")
          parsed(result)
          app_state$parsed_data  <- result
          app_state$upload_done  <- TRUE
          incProgress(0.2, detail = "Done")
        }, error = function(e) {
          parsed(NULL)
          app_state$upload_done <- FALSE
          showNotification(paste("Error parsing file:", conditionMessage(e)),
                           type = "error", duration = 10)
        })
      })
    })
    output$parse_status <- renderUI({
      if (is.null(parsed())) return(NULL)
      p <- parsed()
      rename_warning <- NULL
      if (!is.null(p$name_changes) && nrow(p$name_changes) > 0) {
        shown <- head(p$name_changes, 8)
        rename_warning <- tags$div(
          class = "alert alert-warning",
          " Some sample or group names were renamed to valid R names for DEqMS/limma contrasts.",
          tags$ul(lapply(seq_len(nrow(shown)), function(i) {
            tags$li(
              shown$Type[i], ": ",
              tags$code(shown$Original[i]), " -> ",
              tags$code(shown$Renamed[i])
            )
          })),
          if (nrow(p$name_changes) > nrow(shown)) {
            tags$p(tags$small(sprintf(
              "...and %d more rename(s). See the Sample Overview table for active names.",
              nrow(p$name_changes) - nrow(shown)
            )))
          }
        )
      }

      tagList(
        tags$div(class = "alert alert-success",
                 " File parsed successfully"),
        rename_warning
      )
    })
    output$sample_summary <- renderUI({
      req(parsed())
      p  <- parsed()
      ns_names <- p$sample_names
      grps <- infer_groups(ns_names)
      tagList(
        tags$p(strong("Experiment: "), p$experiment_name),
        tags$p(strong("Samples detected: "), length(ns_names)),
        tags$p(strong("Groups detected: "),
               paste(sort(unique(grps)), collapse = ", "))
      )
    })

    output$sample_table <- DT::renderDataTable({
      req(parsed())
      p    <- parsed()
      grps <- infer_groups(p$sample_names)
      df   <- data.frame(
        Sample = p$sample_names,
        Group  = grps,
        stringsAsFactors = FALSE
      )
      DT::datatable(df, rownames = FALSE, options = list(pageLength = 20, dom = "tp"))
    })
    output$protein_summary <- renderUI({
      req(parsed())
      p <- parsed()
      n_complete <- sum(complete.cases(p$intensity))
      n_total    <- nrow(p$intensity)
      pct_mv     <- round(100 * mean(is.na(p$intensity)), 1)
      tagList(
        tags$p(strong("Proteins: "), n_total),
        tags$p(strong("Proteins with no missing values: "), n_complete,
               sprintf(" (%.1f%%)", 100 * n_complete / n_total)),
        tags$p(strong("Overall missing value rate: "), paste0(pct_mv, "%"))
      )
    })

    output$meta_preview <- DT::renderDataTable({
      req(parsed())
      p   <- parsed()
      int <- p$intensity
      mv  <- rowSums(is.na(int))
      df  <- cbind(p$meta,
                   Missing_Values = mv,
                   Valid_Values   = ncol(int) - mv)
      DT::datatable(df, rownames = FALSE,
                    options = list(pageLength = 15, scrollX = TRUE, dom = "ftip"))
    })
    output$ttest_summary <- renderUI({
      req(parsed())
      tt <- parsed()$ttest
      if (!length(tt)) {
        return(tags$div(class = "alert alert-warning",
                        " No t-test columns detected in the first sheet."))
      }
      tags$p(strong("T-test comparisons detected: "), length(tt), ". ",
             paste(names(tt), collapse = ", "))
    })

    output$ttest_preview <- DT::renderDataTable({
      req(parsed())
      tt <- parsed()$ttest
      if (!length(tt)) return(data.frame())
      # Preview the strongest entries from the first detected comparison.
      df <- tt[[1]]
      df <- df[!is.na(df$logFC) & !is.na(df$negLog10p), ]
      df <- head(df[order(df$negLog10p, decreasing = TRUE), ], 200)
      DT::datatable(df, rownames = FALSE,
                    options = list(pageLength = 10, scrollX = TRUE, dom = "ftip"))
    })
    output$mv_heatmap <- renderPlot({
      req(parsed())
      int <- parsed()$intensity
      grps <- infer_groups(colnames(int))

      # Missingness is easier to inspect by biological/sample group than by run.
      uniq_grps <- sort(unique(grps))
      mv_mat <- vapply(uniq_grps, function(g) {
        cols <- colnames(int)[grps == g]
        rowMeans(is.na(int[, cols, drop = FALSE]))
      }, numeric(nrow(int)))

      # Keep the heatmap readable on large exports.
      top_idx <- order(apply(mv_mat, 1, var), decreasing = TRUE)[seq_len(min(100, nrow(mv_mat)))]
      plot_mat <- mv_mat[top_idx, , drop = FALSE]

      # Base graphics keeps this screen independent of optional plotting packages.
      par(mar = c(5, 2, 3, 1))
      image(t(plot_mat), col = colorRampPalette(c("steelblue", "white", "tomato"))(50),
            xaxt = "n", yaxt = "n",
            main = "Missing Value Rate by Group\n(top 100 most variable proteins)",
            xlab = "Group", ylab = "")
      axis(1, at = seq(0, 1, length.out = length(uniq_grps)),
           labels = uniq_grps, las = 2, cex.axis = 0.85)
    })

    output$mv_stats <- renderUI({
      req(parsed())
      int  <- parsed()$intensity
      grps <- infer_groups(colnames(int))
      rows_stats <- lapply(sort(unique(grps)), function(g) {
        cols <- colnames(int)[grps == g]
        pct  <- round(100 * mean(is.na(int[, cols, drop = FALSE])), 1)
        tags$li(strong(g), ": ", paste0(pct, "% missing"))
      })
      tagList(
        tags$p(strong("Missing value rate per group:")),
        tags$ul(rows_stats)
      )
    })
    output$proceed_btn <- renderUI({
      req(parsed())
      actionButton(ns("btn_proceed"), "Proceed to Imputation",
                   class = "btn-success btn-lg app-proceed-btn")
    })

    observeEvent(input$btn_proceed, {
      app_state$active_tab <- "imputation"
    })

  })
}
