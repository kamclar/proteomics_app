# Compare uploaded Spectronaut t-test columns against the DEqMS results.
mod_comparison_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Step 8. T-test vs DEqMS Comparison"),
    uiOutput(ns("main_ui"))
  )
}

mod_comparison_server <- function(id, app_state) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    comparison_data <- reactiveVal(NULL)
    common_comparisons <- reactive({
      has_ttest <- isTRUE(app_state$upload_done) &&
        !is.null(app_state$parsed_data$ttest) &&
        length(app_state$parsed_data$ttest) > 0

      has_deqms <- isTRUE(app_state$deqms_done) &&
        !is.null(app_state$deqms_results) &&
        length(app_state$deqms_results) > 0

      if (!has_ttest || !has_deqms) return(character(0))
      intersect(names(app_state$parsed_data$ttest), names(app_state$deqms_results))
    })
    
    output$main_ui <- renderUI({
      has_ttest <- isTRUE(app_state$upload_done) &&
                   !is.null(app_state$parsed_data$ttest) &&
                   length(app_state$parsed_data$ttest) > 0
      
      has_deqms <- isTRUE(app_state$deqms_done) &&
                   !is.null(app_state$deqms_results) &&
                   length(app_state$deqms_results) > 0
      
      if (!has_ttest || !has_deqms) {
        return(tags$div(
          class = "alert alert-warning",
          h4("Both T-test and DEqMS results required"),
          p("Please complete:"),
          tags$ul(
            if (!has_ttest) tags$li("Step 1: Upload data with t-test results"),
            if (!has_deqms) tags$li("Step 4: Run DEqMS analysis")
          )
        ))
      }
      
      common_comps <- common_comparisons()
      
      if (length(common_comps) == 0) {
        return(tags$div(
          class = "alert alert-info",
          "No matching comparisons found between t-test and DEqMS results."
        ))
      }
      
      tagList(
        fluidRow(
          column(3,
            wellPanel(
              h4("Comparison Settings"),
              selectInput(ns("comparison"), "Select Comparison:",
                          choices = common_comps,
                          selected = common_comps[1]),
              hr(),
              numericInput(ns("fc_threshold"), "log2FC threshold:",
                           value = 0.5, min = 0, max = 5, step = 0.1),
              numericInput(ns("p_threshold"), "-log10(p) threshold:",
                           value = 1.3, min = 0, max = 10, step = 0.1),
              helpText("Thresholds applied to both t-test and DEqMS results."),
              hr(),
              actionButton(ns("btn_compare"), "Run Comparison",
                           class = "btn-primary btn-block"),
              br(),
              downloadButton(ns("dl_report"), "Download Comparison PDF",
                             class = "btn-default btn-block btn-sm"),
              downloadButton(ns("dl_discordant"), "Download Discordant Genes",
                             class = "btn-default btn-block btn-sm")
            )
          ),
          column(9,
            uiOutput(ns("summary_stats")),
            
            # Overlap summaries
            fluidRow(
              column(6,
                h5("Up-regulated Genes"),
                plotOutput(ns("venn_up"), height = "300px")
              ),
              column(6,
                h5("Down-regulated Genes"),
                plotOutput(ns("venn_down"), height = "300px")
              )
            ),
            
            hr(),
            
            # Agreement plots
            fluidRow(
              column(6,
                h5("log2 Fold Change Comparison"),
                plotOutput(ns("scatter_fc"), height = "400px")
              ),
              column(6,
                h5("Significance Comparison"),
                plotOutput(ns("scatter_pval"), height = "400px")
              )
            ),
            
            hr(),
            
            # Discordant genes
            h5("Discordant Genes"),
            DT::dataTableOutput(ns("discordant_table"))
          )
        )
      )
    })
    build_comparison_result <- function(comp, fc_t, p_t) {
      ttest_df <- app_state$parsed_data$ttest[[comp]]
      deqms_df <- app_state$deqms_results[[comp]]

      if (is.null(ttest_df) || is.null(deqms_df)) {
        stop("Comparison data not found")
      }

      gene_col_ttest <- if ("Genes" %in% names(ttest_df)) "Genes" else NULL
      gene_col_deqms <- if ("Genes" %in% names(deqms_df)) "Genes" else NULL

      if (is.null(gene_col_ttest) || is.null(gene_col_deqms)) {
        stop("Comparison requires a Genes column in both t-test and DEqMS results. Row-name matching is not used for biological comparison.")
      }
      ttest_cols <- intersect(c("Genes", "logFC", "negLog10p"), names(ttest_df))
      deqms_cols <- intersect(c("Genes", "logFC", "negLog10_sca_adjPval", "sca.adj.pval"), names(deqms_df))

      merged <- merge(
        ttest_df[, ttest_cols, drop = FALSE],
        deqms_df[, deqms_cols, drop = FALSE],
        by = "Genes",
        suffixes = c("_ttest", "_deqms")
      )

      if ("negLog10p" %in% names(merged)) {
        names(merged)[names(merged) == "negLog10p"] <- "negLog10p_ttest"
      }
      if ("negLog10_sca_adjPval" %in% names(merged)) {
        names(merged)[names(merged) == "negLog10_sca_adjPval"] <- "negLog10p_deqms"
      }

      merged$ttest_sig_up <- merged$negLog10p_ttest >= p_t & merged$logFC_ttest >= fc_t
      merged$ttest_sig_down <- merged$negLog10p_ttest >= p_t & merged$logFC_ttest <= -fc_t
      merged$deqms_sig_up <- merged$negLog10p_deqms >= p_t & merged$logFC_deqms >= fc_t
      merged$deqms_sig_down <- merged$negLog10p_deqms >= p_t & merged$logFC_deqms <= -fc_t

      merged$ttest_sig <- merged$ttest_sig_up | merged$ttest_sig_down
      merged$deqms_sig <- merged$deqms_sig_up | merged$deqms_sig_down

      merged$category <- ifelse(
        merged$ttest_sig & merged$deqms_sig, "Both Significant",
        ifelse(merged$ttest_sig & !merged$deqms_sig, "T-test Only",
               ifelse(!merged$ttest_sig & merged$deqms_sig, "DEqMS Only",
                      "Both NS"))
      )

      merged$same_direction <- sign(merged$logFC_ttest) == sign(merged$logFC_deqms)

      stats <- list(
        total_genes = nrow(merged),
        ttest_sig_count = sum(merged$ttest_sig, na.rm = TRUE),
        deqms_sig_count = sum(merged$deqms_sig, na.rm = TRUE),
        both_sig_count = sum(merged$ttest_sig & merged$deqms_sig, na.rm = TRUE),
        ttest_only_count = sum(merged$ttest_sig & !merged$deqms_sig, na.rm = TRUE),
        deqms_only_count = sum(!merged$ttest_sig & merged$deqms_sig, na.rm = TRUE),
        up_ttest = sum(merged$ttest_sig_up, na.rm = TRUE),
        up_deqms = sum(merged$deqms_sig_up, na.rm = TRUE),
        up_both = sum(merged$ttest_sig_up & merged$deqms_sig_up, na.rm = TRUE),
        down_ttest = sum(merged$ttest_sig_down, na.rm = TRUE),
        down_deqms = sum(merged$deqms_sig_down, na.rm = TRUE),
        down_both = sum(merged$ttest_sig_down & merged$deqms_sig_down, na.rm = TRUE),
        direction_agreement = sum(merged$same_direction & merged$ttest_sig & merged$deqms_sig, na.rm = TRUE)
      )

      discordant <- merged[
        (merged$ttest_sig & !merged$deqms_sig) |
        (!merged$ttest_sig & merged$deqms_sig) |
        (merged$ttest_sig & merged$deqms_sig & !merged$same_direction),
      ]

      list(
        merged = merged,
        stats = stats,
        discordant = discordant,
        fc_threshold = fc_t,
        p_threshold = p_t,
        comparison = comp
      )
    }


    observeEvent(input$btn_compare, {
      req(input$comparison)
      
      comp <- input$comparison
      fc_t <- input$fc_threshold
      p_t <- input$p_threshold
      
      withProgress(message = "Comparing results...", value = 0, {
        incProgress(0.20, detail = "Loading selected comparison")
        result <- tryCatch(
          build_comparison_result(comp, fc_t, p_t),
          error = function(e) {
            showNotification(paste("Comparison error:", conditionMessage(e)), type = "error")
            NULL
          }
        )
        incProgress(0.55, detail = "Calculating agreement and discordant proteins")
        if (!is.null(result)) comparison_data(result)
        incProgress(0.25, detail = "Preparing plots and tables")
      })
    })
    output$summary_stats <- renderUI({
      data <- comparison_data()
      if (is.null(data)) {
        return(tags$div(
          class = "alert alert-info",
          "Select a comparison and click Run Comparison. Results, plots, and tables will appear here."
        ))
      }
      
      stats <- data$stats
      
      pct_ttest_confirmed <- round(100 * stats$both_sig_count / stats$ttest_sig_count, 1)
      pct_deqms_confirmed <- round(100 * stats$both_sig_count / stats$deqms_sig_count, 1)
      pct_direction_agree <- round(100 * stats$direction_agreement / stats$both_sig_count, 1)
      
      tags$div(
        class = "well",
        h4("Agreement Summary"),
        fluidRow(
          column(3,
            tags$div(style = "text-align: center; padding: 10px; background: #e8f4f8; border-radius: 5px;",
              h3(stats$total_genes, style = "margin: 5px;"),
              p("Total Genes", style = "margin: 0;")
            )
          ),
          column(3,
            tags$div(style = "text-align: center; padding: 10px; background: #d4edda; border-radius: 5px;",
              h3(stats$both_sig_count, style = "margin: 5px;"),
              p("Both Significant", style = "margin: 0;")
            )
          ),
          column(3,
            tags$div(style = "text-align: center; padding: 10px; background: #fff3cd; border-radius: 5px;",
              h3(stats$ttest_only_count, style = "margin: 5px;"),
              p("T-test Only", style = "margin: 0;")
            )
          ),
          column(3,
            tags$div(style = "text-align: center; padding: 10px; background: #f8d7da; border-radius: 5px;",
              h3(stats$deqms_only_count, style = "margin: 5px;"),
              p("DEqMS Only", style = "margin: 0;")
            )
          )
        ),
        hr(),
        fluidRow(
          column(4,
            p(strong("T-test -> DEqMS confirmation:"), paste0(pct_ttest_confirmed, "%")),
            p(paste0(stats$both_sig_count, " of ", stats$ttest_sig_count, " t-test hits confirmed by DEqMS"))
          ),
          column(4,
            p(strong("DEqMS -> T-test confirmation:"), paste0(pct_deqms_confirmed, "%")),
            p(paste0(stats$both_sig_count, " of ", stats$deqms_sig_count, " DEqMS hits found in t-test"))
          ),
          column(4,
            p(strong("Direction agreement:"), paste0(pct_direction_agree, "%")),
            p(paste0(stats$direction_agreement, " of ", stats$both_sig_count, " shared hits have same direction"))
          )
        )
      )
    })
    # Overlap summaries
    output$venn_up <- renderPlot({
      data <- comparison_data()
      if (is.null(data)) {
        plot.new()
        text(0.5, 0.5, "Run comparison to show up-regulated overlap.", cex = 1.1)
        return()
      }
      
      stats <- data$stats
      
      par(mar = c(1, 1, 3, 1))
      plot.new()
      plot.window(xlim = c(0, 10), ylim = c(0, 10))
      
      symbols(x = c(3.5, 6.5), y = c(5, 5), circles = c(2, 2),
              inches = FALSE, add = TRUE, fg = c("#F76075", "#40BFC1"),
              lwd = 3)
      
      text(2, 5, stats$up_ttest - stats$up_both, cex = 2, font = 2)
      text(8, 5, stats$up_deqms - stats$up_both, cex = 2, font = 2)
      text(5, 5, stats$up_both, cex = 2, font = 2)
      
      text(2, 8.5, "T-test", col = "#F76075", cex = 1.2, font = 2)
      text(8, 8.5, "DEqMS", col = "#40BFC1", cex = 1.2, font = 2)
      
      title(paste0("Up-regulated (n=", stats$up_ttest + stats$up_deqms - stats$up_both, ")"),
            cex.main = 1.3)
    })
    
    output$venn_down <- renderPlot({
      data <- comparison_data()
      if (is.null(data)) {
        plot.new()
        text(0.5, 0.5, "Run comparison to show down-regulated overlap.", cex = 1.1)
        return()
      }
      
      stats <- data$stats
      
      par(mar = c(1, 1, 3, 1))
      plot.new()
      plot.window(xlim = c(0, 10), ylim = c(0, 10))
      
      symbols(x = c(3.5, 6.5), y = c(5, 5), circles = c(2, 2),
              inches = FALSE, add = TRUE, fg = c("#F76075", "#40BFC1"),
              lwd = 3)
      
      text(2, 5, stats$down_ttest - stats$down_both, cex = 2, font = 2)
      text(8, 5, stats$down_deqms - stats$down_both, cex = 2, font = 2)
      text(5, 5, stats$down_both, cex = 2, font = 2)
      
      text(2, 8.5, "T-test", col = "#F76075", cex = 1.2, font = 2)
      text(8, 8.5, "DEqMS", col = "#40BFC1", cex = 1.2, font = 2)
      
      title(paste0("Down-regulated (n=", stats$down_ttest + stats$down_deqms - stats$down_both, ")"),
            cex.main = 1.3)
    })
    # Agreement plots
    output$scatter_fc <- renderPlot({
      data <- comparison_data()
      if (is.null(data)) {
        plot.new()
        text(0.5, 0.5, "Run comparison to show fold-change agreement.", cex = 1.1)
        return()
      }
      
      df <- data$merged
      
      colors <- c("Both Significant" = "#2ecc71", "T-test Only" = "#e74c3c",
                  "DEqMS Only" = "#3498db", "Both NS" = "#bdc3c7")
      
      plot(df$logFC_ttest, df$logFC_deqms,
           col = colors[df$category],
           pch = 19, cex = 0.8,
           xlab = "T-test log2FC",
           ylab = "DEqMS log2FC",
           main = "Fold Change Correlation")
      
      abline(a = 0, b = 1, col = "black", lty = 2, lwd = 2)
      abline(h = 0, v = 0, col = "gray", lty = 3)
      
      cor_val <- cor(df$logFC_ttest, df$logFC_deqms, use = "complete.obs")
      legend("topleft", 
             legend = c(paste0("r = ", round(cor_val, 3)),
                       names(colors)),
             col = c("black", colors),
             pch = c(NA, rep(19, 4)),
             lty = c(2, rep(NA, 4)),
             bty = "n")
    })
    
    output$scatter_pval <- renderPlot({
      data <- comparison_data()
      if (is.null(data)) {
        plot.new()
        text(0.5, 0.5, "Run comparison to show significance agreement.", cex = 1.1)
        return()
      }
      
      df <- data$merged
      
      colors <- c("Both Significant" = "#2ecc71", "T-test Only" = "#e74c3c",
                  "DEqMS Only" = "#3498db", "Both NS" = "#bdc3c7")
      
      plot(df$negLog10p_ttest, df$negLog10p_deqms,
           col = colors[df$category],
           pch = 19, cex = 0.8,
           xlab = "T-test -log10(p)",
           ylab = "DEqMS -log10(adj.p)",
           main = "Significance Correlation")
      
      abline(a = 0, b = 1, col = "black", lty = 2, lwd = 2)
      abline(h = data$p_threshold, v = data$p_threshold, 
             col = "red", lty = 3, lwd = 1.5)
      
      cor_val <- cor(df$negLog10p_ttest, df$negLog10p_deqms, use = "complete.obs")
      legend("topleft",
             legend = c(paste0("r = ", round(cor_val, 3)),
                       names(colors)),
             col = c("black", colors),
             pch = c(NA, rep(19, 4)),
             lty = c(2, rep(NA, 4)),
             bty = "n")
    })
    # Discordant genes
    output$discordant_table <- DT::renderDataTable({
      data <- comparison_data()
      if (is.null(data)) {
        return(data.frame(Message = "Run comparison to show discordant proteins."))
      }
      
      df <- data$discordant
      
      if (nrow(df) == 0) {
        return(data.frame(Message = "No discordant genes found!"))
      }
      
      show_cols <- c("Genes", "logFC_ttest", "negLog10p_ttest", 
                     "logFC_deqms", "negLog10p_deqms", "category")
      show_cols <- intersect(show_cols, names(df))
      
      df_show <- df[, show_cols, drop = FALSE]
      
      DT::datatable(df_show, rownames = FALSE,
                    options = list(pageLength = 15, scrollX = TRUE, dom = "ftip")) |>
        DT::formatRound(c("logFC_ttest", "logFC_deqms", 
                         "negLog10p_ttest", "negLog10p_deqms"), 3)
    })
    output$dl_discordant <- downloadHandler(
      filename = function() {
        data <- comparison_data()
        comp <- if (is.null(data$comparison) || !nzchar(data$comparison)) "comparison" else data$comparison
        suffix <- paste0("discordant_", gsub("[^A-Za-z0-9_]+", "_", comp))
        download_filename(app_state, suffix, "csv")
      },
      content = function(file) {
        data <- comparison_data()
        req(!is.null(data))
        utils::write.csv(data$discordant, file, row.names = FALSE)
      }
    )

    output$dl_report <- downloadHandler(
      filename = function() download_filename(app_state, "comparison_report", "pdf"),
      content = function(file) {
        comps <- common_comparisons()
        req(length(comps) > 0)

        report <- vector("list", length(comps))
        withProgress(message = "Building comparison PDF...", value = 0, style = "notification", {
          for (i in seq_along(comps)) {
            incProgress(1 / length(comps), detail = comps[[i]])
            report[[i]] <- build_comparison_result(comps[[i]], input$fc_threshold, input$p_threshold)
          }
        })

        colors <- c(
          "Both Significant" = "#2ecc71",
          "T-test Only" = "#e74c3c",
          "DEqMS Only" = "#3498db",
          "Both NS" = "#bdc3c7"
        )
        total_pages <- 1 + 2 * length(report)
        page_num <- 1L
        experiment_name <- download_experiment_prefix(app_state)
        imp <- app_state$imputed_data
        imputation_lines <- character(0)
        if (!is.null(imp)) {
          imputation_lines <- c(
            paste("Imputation method:", if (!is.null(imp$method)) imp$method else "unknown"),
            if (!is.null(imp$mnar_backend) && !is.na(imp$mnar_backend)) {
              paste("MNAR backend:", imp$mnar_backend, "(imputeLCMD)")
            },
            if (!is.null(imp$seed) && !is.na(imp$seed)) {
              paste("Random seed:", imp$seed, "(reproduces random MNAR draws)")
            },
            if (!is.null(imp$rows_dropped)) {
              paste("Rows removed by pre-filter:", imp$rows_dropped)
            }
          )
        }

        add_page_footer <- function(page_num, total_pages) {
          graphics::mtext(
            paste("Page", page_num, "of", total_pages),
            side = 1, line = 1.5, adj = 1, cex = 0.85, col = "gray35"
          )
        }

        grDevices::pdf(file, width = 11, height = 8.5, onefile = TRUE)
        on.exit(grDevices::dev.off(), add = TRUE)

        graphics::par(mfrow = c(1, 1), mar = c(2, 2, 2, 2))
        graphics::plot.new()
        graphics::rect(0.08, 0.72, 0.92, 0.78, col = "#2C7FB8", border = NA)
        graphics::text(0.5, 0.84, labels = experiment_name, cex = 2.4, font = 2, col = "#0B3C5D")
        graphics::text(0.5, 0.72, labels = "T-test vs DEqMS Comparison Report", cex = 1.6, font = 3)
        graphics::text(
          0.5, 0.56,
          labels = paste("Matching comparisons:", length(report)),
          cex = 1.2
        )
        graphics::text(
          0.5, 0.49,
          labels = paste("Thresholds: |log2FC| >=", input$fc_threshold, ", -log10(p) >=", input$p_threshold),
          cex = 1.1
        )
        graphics::text(
          0.5, 0.42,
          labels = paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M")),
          cex = 1.1
        )
        if (length(imputation_lines)) {
          graphics::text(
            0.5, 0.34,
            labels = paste(imputation_lines, collapse = "\n"),
            cex = 0.95
          )
        }
        graphics::text(
          0.5, 0.24,
          labels = "Each comparison includes overlap summaries, fold-change and significance plots, and discordant-gene review.",
          cex = 0.95
        )
        add_page_footer(page_num, total_pages)
        page_num <- page_num + 1L

        for (data in report) {
          stats <- data$stats
          df <- data$merged

          graphics::par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

          # Page 1, panel 1: summary text
          graphics::plot.new()
          graphics::text(0.5, 0.94, labels = data$comparison, cex = 2, font = 2)
          graphics::text(0.5, 0.86, labels = "T-test vs DEqMS Comparison Report", cex = 1.1, font = 3)
          summary_lines <- c(
            paste("Thresholds: |log2FC| >=", data$fc_threshold, ", -log10(p) >=", data$p_threshold),
            paste("Total genes:", stats$total_genes),
            paste("Both significant:", stats$both_sig_count),
            paste("T-test only:", stats$ttest_only_count),
            paste("DEqMS only:", stats$deqms_only_count),
            paste("Direction agreement:", stats$direction_agreement)
          )
          graphics::text(0.05, seq(0.72, 0.28, length.out = length(summary_lines)),
                         labels = summary_lines, adj = c(0, 0.5), cex = 1)

          # Page 1, panel 2: up venn
          graphics::plot.new()
          graphics::plot.window(xlim = c(0, 10), ylim = c(0, 10))
          graphics::symbols(x = c(3.5, 6.5), y = c(5, 5), circles = c(2, 2),
                            inches = FALSE, add = TRUE, fg = c("#F76075", "#40BFC1"), lwd = 3)
          graphics::text(2, 5, stats$up_ttest - stats$up_both, cex = 1.8, font = 2)
          graphics::text(8, 5, stats$up_deqms - stats$up_both, cex = 1.8, font = 2)
          graphics::text(5, 5, stats$up_both, cex = 1.8, font = 2)
          graphics::text(2, 8.5, "T-test", col = "#F76075", cex = 1.1, font = 2)
          graphics::text(8, 8.5, "DEqMS", col = "#40BFC1", cex = 1.1, font = 2)
          graphics::title("Up-regulated Genes")

          # Page 1, panel 3: down venn
          graphics::plot.new()
          graphics::plot.window(xlim = c(0, 10), ylim = c(0, 10))
          graphics::symbols(x = c(3.5, 6.5), y = c(5, 5), circles = c(2, 2),
                            inches = FALSE, add = TRUE, fg = c("#F76075", "#40BFC1"), lwd = 3)
          graphics::text(2, 5, stats$down_ttest - stats$down_both, cex = 1.8, font = 2)
          graphics::text(8, 5, stats$down_deqms - stats$down_both, cex = 1.8, font = 2)
          graphics::text(5, 5, stats$down_both, cex = 1.8, font = 2)
          graphics::text(2, 8.5, "T-test", col = "#F76075", cex = 1.1, font = 2)
          graphics::text(8, 8.5, "DEqMS", col = "#40BFC1", cex = 1.1, font = 2)
          graphics::title("Down-regulated Genes")

          # Page 1, panel 4: FC scatter
          graphics::plot(df$logFC_ttest, df$logFC_deqms,
                         col = colors[df$category], pch = 19, cex = 0.6,
                         xlab = "T-test log2FC", ylab = "DEqMS log2FC",
                         main = "Fold Change Correlation")
          graphics::abline(a = 0, b = 1, col = "black", lty = 2, lwd = 2)
          graphics::abline(h = 0, v = 0, col = "gray", lty = 3)
          add_page_footer(page_num, total_pages)
          page_num <- page_num + 1L

          # Page 2: significance scatter + discordant table
          graphics::par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
          graphics::plot(df$negLog10p_ttest, df$negLog10p_deqms,
                         col = colors[df$category], pch = 19, cex = 0.6,
                         xlab = "T-test -log10(p)", ylab = "DEqMS -log10(adj.p)",
                         main = paste("Significance Correlation -", data$comparison))
          graphics::abline(a = 0, b = 1, col = "black", lty = 2, lwd = 2)
          graphics::abline(h = data$p_threshold, v = data$p_threshold, col = "red", lty = 3)

          disc <- data$discordant
          graphics::plot.new()
          graphics::text(0.5, 0.94, labels = data$comparison, cex = 2, font = 2)
          graphics::text(0.5, 0.86, labels = "Discordant Genes", cex = 1.2, font = 3)
          if (!nrow(disc)) {
            graphics::text(0.5, 0.5, "No discordant genes found.", cex = 1.2)
          } else {
            show_cols <- intersect(c("Genes", "logFC_ttest", "negLog10p_ttest",
                                     "logFC_deqms", "negLog10p_deqms", "category"), names(disc))
            disc_show <- head(disc[, show_cols, drop = FALSE], 20)
            lines <- c(
              paste(names(disc_show), collapse = " | "),
              apply(format(disc_show, digits = 3), 1, paste, collapse = " | ")
            )
            graphics::text(0.02, seq(0.74, 0.08, length.out = length(lines)),
                           labels = lines,
                           adj = c(0, 0.5), cex = 0.7)
          }
          add_page_footer(page_num, total_pages)
          page_num <- page_num + 1L
        }
      }
    )
  })
}
