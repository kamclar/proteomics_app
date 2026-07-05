# ProteomicsApp
#
# Main Shiny shell. The analysis steps live in modules/ so the individual
# screens stay manageable.

local_lib <- file.path(getwd(), "packages")
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
  message("Using local package library: ", local_lib)
} else if (file.exists("library_config.R")) {
  source("library_config.R")
}

optional_lib <- Sys.getenv("PROTEOMICS_APP_R_LIB", unset = "")
if (!nzchar(optional_lib)) {
  r_minor <- sub("\\..*$", "", R.version$minor)
  optional_lib <- file.path(
    tools::R_user_dir("ProteomicsApp", "data"),
    "R",
    paste0("win-library-", R.version$major, ".", r_minor)
  )
}
dir.create(optional_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(unique(normalizePath(c(
  if (dir.exists(local_lib)) local_lib,
  optional_lib,
  .libPaths()
), winslash = "/", mustWork = FALSE)))
message("Using R package libraries: ", paste(.libPaths(), collapse = " | "))

library(shiny)
library(shinydashboard)
library(DT)
library(plotly)
library(openxlsx)
library(limma)
library(DEqMS)
library(ggplot2)
library(dplyr)

options(shiny.maxRequestSize = 100*1024^2) 
options(shiny.notification.position = "top-center")

big_dot <- function() {
  tags$span(class = "ui-big-dot", HTML("&bull;"))
}

dot_label <- function(label) {
  tagList(big_dot(), label)
}

source("modules/utils_data.R", local = TRUE)
source("modules/mod_upload.R", local = TRUE)
source("modules/mod_imputation.R", local = TRUE)
source("modules/mod_exploration.R", local = TRUE)
source("modules/mod_deqms.R", local = TRUE)
source("modules/mod_volcano.R", local = TRUE)
source("modules/mod_enrichment.R", local = TRUE)
source("modules/mod_background.R", local = TRUE)
source("modules/mod_comparison.R", local = TRUE)

ui <- shinydashboard::dashboardPage(
  skin = "blue",

  header = shinydashboard::dashboardHeader(
    title = tags$span(
      tags$img(src = "logo.png", height = "30px", style = "margin-right:8px;"),
      "ProteomicsApp"
    ),
    titleWidth = 260
  ),

  sidebar = shinydashboard::dashboardSidebar(
    width = 260,
    shinydashboard::sidebarMenu(
      id = "sidebar_tabs",
      shinydashboard::menuItem(dot_label("1 - Upload"),          tabName = "upload"),
      shinydashboard::menuItem(dot_label("2 - Imputation"),       tabName = "imputation"),
      shinydashboard::menuItem(dot_label("3 - Exploration"),      tabName = "exploration"),
      shinydashboard::menuItem(dot_label("4 - DEqMS"),            tabName = "deqms"),
      shinydashboard::menuItem(dot_label("5 - Volcano Plotter"),  tabName = "volcano"),
      shinydashboard::menuItem(dot_label("6 - Enrichment"),       tabName = "enrichment"),
      shinydashboard::menuItem(dot_label("7 - Stable Proteins"),  tabName = "background"),
      shinydashboard::menuItem(dot_label("8 - Comparison"),       tabName = "comparison"),
      shinydashboard::menuItem(dot_label("Help / About"),         tabName = "help")
      
    ),
    tags$div(
      style = "padding: 10px 16px; color: #aaa; font-size: 11px; position: absolute; bottom: 0;",
      "Spectronaut -> DEqMS pipeline",
      tags$br(),
      tags$a(href = "#", style = "color: #aaa;", "v0.4.5")
    )
  ),

  body = shinydashboard::dashboardBody(
    tags$head(
      tags$style(HTML("
        body {
          color: #1f2933;
          background: #eef2f6;
        }
        .skin-blue .main-header .logo,
        .skin-blue .main-header .navbar {
          background: #f7fafc;
          color: #102a43;
          border-bottom: 1px solid #d9e2ec;
          box-shadow: 0 10px 30px rgba(15, 23, 42, 0.05);
        }
        .skin-blue .main-header .logo {
          font-weight: 700;
        }
        .skin-blue .main-header .logo:hover,
        .skin-blue .main-header .navbar:hover {
          background: #f7fafc;
        }
        .skin-blue .main-header .sidebar-toggle {
          color: #486581;
        }
        .skin-blue .main-header .sidebar-toggle:hover {
          background: #e6eef5;
          color: #102a43;
        }
        .skin-blue .main-sidebar,
        .skin-blue .left-side {
          background: linear-gradient(180deg, #65137d 0%, #6d2e80 100%);
          box-shadow: inset -1px 0 0 rgba(255,255,255,0.04);
        }
        .skin-blue .sidebar a {
          color: #d9e2ec;
        }
        .skin-blue .sidebar-menu > li.header {
          color: #9fb3c8;
          background: transparent;
        }
        .skin-blue .sidebar-menu > li > a {
          margin: 4px 10px;
          border-radius: 6px;
          padding: 12px 14px;
          font-size: 13px;
          font-weight: 600;
          transition: background 0.18s ease, color 0.18s ease, transform 0.18s ease;
        }
        .skin-blue .sidebar-menu > li > a:hover {
          background: rgba(255,255,255,0.08);
          color: #f0f4f8;
        }
        .skin-blue .sidebar-menu > li.active > a {
          background: linear-gradient(135deg, #c69adb 0%, #6f8db5 100%);
          color: white;
          box-shadow: 0 10px 24px rgba(72, 101, 129, 0.28);
        }
        .content-wrapper,
        .right-side {
          background: linear-gradient(180deg, #f4f7fb 0%, #edf2f7 100%);
        }
        .content {
          padding: 20px 22px 28px 22px;
        }
        .box,
        .well {
          border: 1px solid #d9e2ec;
          border-radius: 6px;
          background: rgba(255,255,255,0.96);
          box-shadow: 0 10px 26px rgba(15, 23, 42, 0.06);
        }
        .box {
          border-top: 1px solid #d9e2ec;
        }
        .well {
          padding: 18px 18px 16px 18px;
        }
        .btn-block { width: 100%; }
        .btn {
          border-radius: 6px;
          border-width: 1px;
          font-weight: 600;
          letter-spacing: 0;
          transition: all 0.18s ease;
        }
        .btn:focus,
        .btn:active:focus {
          outline: none;
          box-shadow: 0 0 0 3px rgba(20, 184, 166, 0.14);
        }
        .btn.btn-primary,
        .action-button.btn-primary {
          background: linear-gradient(135deg, #14b8a6 0%, #0f9f8f 100%);
          border-color: #0f9f8f;
          color: #fff;
          box-shadow: 0 10px 22px rgba(15, 159, 143, 0.22);
        }
        .btn.btn-primary:hover,
        .btn.btn-primary:active,
        .btn.btn-primary:focus,
        .action-button.btn-primary:hover,
        .action-button.btn-primary:active,
        .action-button.btn-primary:focus {
          background: linear-gradient(135deg, #11998a 0%, #0d8b7c 100%);
          border-color: #0d8b7c;
          color: #fff;
          transform: translateY(-1px);
        }
        .btn.btn-success,
        .btn-default.btn-success,
        .action-button.btn-success {
          background: linear-gradient(135deg, #2f855a 0%, #276749 100%);
          border-color: #276749;
          color: #fff;
          box-shadow: 0 12px 26px rgba(47, 133, 90, 0.24);
        }
        .btn.btn-success:hover,
        .btn.btn-success:focus,
        .btn-default.btn-success:hover,
        .btn-default.btn-success:focus,
        .action-button.btn-success:hover,
        .action-button.btn-success:focus {
          background: linear-gradient(135deg, #2b7a53 0%, #21543b 100%);
          border-color: #21543b;
          color: #fff;
          transform: translateY(-1px);
        }
        .btn-default {
          background: #f8fbff;
          border-color: #cbd5e1;
          color: #243b53;
        }
        .btn-default:hover,
        .btn-default:focus {
          background: #eef5fb;
          border-color: #9fb3c8;
          color: #102a43;
        }
        a.btn-default {
          background: #fff8eb;
          border-color: #f2c879;
          color: #8a5a00 !important;
          box-shadow: 0 6px 16px rgba(217, 119, 6, 0.12);
        }
        a.btn-default:hover,
        a.btn-default:focus {
          background: #fff1d6;
          border-color: #e0ae4f;
          color: #6f4700 !important;
          transform: translateY(-1px);
        }
        .app-proceed-btn {
          min-width: 240px;
          font-size: 16px;
          font-weight: 700;
          box-shadow: 0 16px 34px rgba(47, 133, 90, 0.28);
        }
        .app-proceed-btn:hover,
        .app-proceed-btn:focus {
          box-shadow: 0 18px 38px rgba(39, 103, 73, 0.32);
        }
        .app-action-row {
          display: flex;
          justify-content: space-between;
          align-items: center;
          gap: 12px;
          margin-top: 12px;
          flex-wrap: wrap;
        }
        .app-action-row .btn,
        .app-action-row .shiny-download-link {
          white-space: normal;
        }
        .app-action-row .app-proceed-wrap {
          margin-left: auto;
        }
        .btn-xs {
          padding: 4px 10px;
          font-size: 12px;
        }
        .alert {
          border-radius: 6px;
          padding: 12px 14px;
          border: 1px solid #d9e2ec;
        }
        .alert-info {
          background: #eef8fb;
          color: #114b5f;
          border-color: #b8dfe8;
        }
        .alert-warning {
          background: #fff6e5;
          color: #8d5d00;
          border-color: #f3d39a;
        }
        .alert-success {
          background: #f2faf8;
          color: #91bdb3;
          border-color: #91bdb3;
        }
        .progress {
          border-radius: 999px;
          background: #d9e2ec;
          height: 10px;
        }
        .progress-bar {
          transition: width 0.2s ease;
          background: linear-gradient(90deg, #14b8a6 0%, #0f9f8f 100%);
          border-radius: 999px;
        }
        .shiny-file-input-progress,
        .fileinput-upload {
          display: block;
          width: 100%;
          height: 28px;
          margin: 8px 0 10px;
          border-radius: 999px;
          overflow: hidden;
          background: #d9e2ec;
          visibility: visible;
        }
        .shiny-file-input-progress .progress-bar,
        .fileinput-upload .progress-bar {
          min-width: 8em;
          height: 28px;
          line-height: 28px;
          padding: 0 12px;
          font-size: 12px;
          font-weight: 700;
          text-align: center;
          color: #fff;
          white-space: nowrap;
          overflow: hidden;
          text-overflow: ellipsis;
        }
        .shiny-file-input-progress .progress-text,
        .fileinput-upload .progress-text {
          font-size: 12px;
          color: #486581;
          margin-bottom: 6px;
        }
        .nav-tabs {
          border-bottom: 1px solid #d9e2ec;
          margin-bottom: 14px;
        }
        .nav-tabs > li > a {
          border-radius: 6px 6px 0 0;
          color: #486581;
          font-weight: 600;
          border: 1px solid transparent;
        }
        .nav-tabs > li > a:hover {
          background: #f0f7fa;
          color: #102a43;
        }
        .nav-tabs > li.active > a,
        .nav-tabs > li.active > a:hover,
        .nav-tabs > li.active > a:focus {
          border-top: 2px solid #14b8a6;
          border-left: 1px solid #d9e2ec;
          border-right: 1px solid #d9e2ec;
          border-bottom: 1px solid white;
          background: #fff;
          color: #102a43;
          font-weight: 700;
        }
        h3 {
          margin-top: 0;
          color: #102a43;
          border-bottom: 1px solid #d9e2ec;
          padding-bottom: 10px;
          margin-bottom: 18px;
          font-weight: 700;
        }
        h4 {
          color: #243b53;
          font-weight: 700;
        }
        h5 {
          color: #334e68;
          font-weight: 700;
        }
        .dataTables_wrapper .dataTables_filter input,
        .dataTables_wrapper .dataTables_length select,
        .form-control {
          border-radius: 6px;
          border: 1px solid #cbd5e1;
          box-shadow: none;
        }
        .form-control:focus {
          border-color: #14b8a6;
          box-shadow: 0 0 0 3px rgba(20, 184, 166, 0.12);
        }
        .dataTables_wrapper .dataTables_filter input {
          background: white;
        }
        table.dataTable thead th {
          background: #f8fafc;
          color: #243b53;
          border-bottom: 1px solid #d9e2ec !important;
          font-weight: 700;
        }
        table.dataTable tbody td {
          border-color: #e5edf5 !important;
        }
        .dataTables_wrapper .dataTables_paginate .paginate_button {
          border-radius: 6px !important;
        }
        .ui-big-dot {
          display: inline-block;
          margin-right: 7px;
          color: rgba(255,255,255,0.78);
          font-size: 24px;
          line-height: 0.75;
          vertical-align: -0.08em;
          font-weight: 700;
        }
        .btn .ui-big-dot {
          font-size: 22px;
          margin-right: 6px;
        }
        #shiny-notification-panel {
          position: fixed !important;
          top: 50% !important;
          left: 50% !important;
          right: auto !important;
          bottom: auto !important;
          transform: translate(-50%, -50%) !important;
          width: min(560px, calc(100vw - 32px));
          z-index: 3000 !important;
        }
        #shiny-notification-panel .shiny-notification {
          width: 100%;
          margin: 0 0 12px 0;
          border-radius: 8px;
          box-shadow: 0 18px 42px rgba(15, 23, 42, 0.2);
        }
        #shiny-notification-panel .progress {
          margin-bottom: 0;
        }
      "))
    ),

    shinydashboard::tabItems(
      shinydashboard::tabItem(tabName = "upload",
        mod_upload_ui("upload")
      ),
      shinydashboard::tabItem(tabName = "imputation",
        mod_imputation_ui("imputation")
      ),
      shinydashboard::tabItem(tabName = "exploration",
        mod_exploration_ui("exploration")
      ),
      shinydashboard::tabItem(tabName = "deqms",
        mod_deqms_ui("deqms")
      ),
      shinydashboard::tabItem(tabName = "volcano",
        mod_volcano_ui("volcano")
      ),
      shinydashboard::tabItem(tabName = "enrichment",
        mod_enrichment_ui("enrichment")
      ),
      shinydashboard::tabItem(tabName = "background",
        mod_background_ui("background")
      ),
      shinydashboard::tabItem(tabName = "comparison",
              mod_comparison_ui("comparison")
      ),
      shinydashboard::tabItem(tabName = "help",
        fluidRow(
          column(8, offset = 2,
            wellPanel(
              h3("Help & About"),
              h4("Workflow overview"),
              tags$ol(
                tags$li(tags$b("Upload:"), " Load a Spectronaut XLSX file, preview the sheet, choose the header row, and review the auto-detected intensity, metadata, and t-test columns."),
                tags$li(tags$b("Imputation:"), " Handle missing values. Choose MNAR (Perseus default, QRILC, or MinProb), MAR (kNN), or mixed auto-detection. Download imputed CSV."),
                tags$li(tags$b("Exploration:"), " QC with PCA, MDS, correlation heatmap, hierarchical clustering, and optional UMAP. Outlier detection based on within-group correlation MAD."),
                tags$li(tags$b("DEqMS:"), " Run limma + DEqMS for selected pairwise comparisons. Mirror direction is computed instantly without re-fitting. Download all results as ZIP."),
                tags$li(tags$b("Volcano Plotter:"), " Interactive Plotly volcano. Switch between t-test (from upload) and DEqMS results. Flip sides, adjust thresholds, export SVG."),
                tags$li(tags$b("Enrichment:"), " GO (BP/MF/CC), KEGG, and Reactome enrichment on up- and down-regulated proteins per comparison. Enrichment packages are installed on demand to keep the main app distribution small."),
                tags$li(tags$b("Stable Proteins:"), " Review proteins that are non-significant across selected DEqMS comparisons, have small fold changes, sufficient data coverage, and no on/off pattern. Compare them with changed proteins and Reactome pathways."),
                tags$li(tags$b("Comparison:"), " Compare uploaded t-test results against DEqMS results for matching contrasts, review agreement and discordant genes, and export discordant hits.")
              ),
              hr(),
              h4("Column name conventions"),
              tags$p("The app expects Spectronaut output with:"),
              tags$ul(
                tags$li("Metadata columns prefixed with ", tags$code("PG."), " (stripped automatically)"),
                tags$li("Intensity columns matching pattern: ", tags$code("N_ExperimentName.raw.PG.Quantity")),
                tags$li("T-test columns containing: ", tags$code("Student's T-test Difference X"), " and ", tags$code("-Log Student's T-test p-value X"))
              ),
              hr(),
              tags$p(tags$small("Created by Kamila Clarova, IOCB Prague."))
            )
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  app_state <- reactiveValues(
    upload_done      = FALSE,
    imputation_done  = FALSE,
    exploration_done = FALSE,
    deqms_done       = FALSE,
    enrichment_done  = FALSE,

    parsed_data    = NULL,
    imputed_data   = NULL,
    deqms_results  = list(),
    tko_fits       = list(),
    enrich_results = list(),
    active_tab = "upload"
  )



  mod_upload_server(    "upload",      app_state)
  mod_imputation_server("imputation",  app_state)
  mod_exploration_server("exploration", app_state)
  mod_deqms_server("deqms",       app_state)
  mod_volcano_server(   "volcano",     app_state)
  mod_enrichment_server("enrichment",  app_state)
  mod_background_server("background",  app_state)
  mod_comparison_server("comparison", app_state)

  observeEvent(app_state$active_tab, {
    shinydashboard::updateTabItems(session, "sidebar_tabs", app_state$active_tab)
  }, ignoreInit = TRUE)
}

shinyApp(ui = ui, server = server)
 
