# SQuality control: PCA, MDS, correlations, clustering, optional UMAP.


mod_exploration_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Step 3. Exploratory Analysis"),
    uiOutput(ns("guard")),
    uiOutput(ns("main_ui"))
  )
}

mod_exploration_server <- function(id, app_state) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    explore_results <- reactiveVal(NULL)

    output$guard <- renderUI({
      if (!isTRUE(app_state$imputation_done)) {
        tags$div(class = "alert alert-warning",
                 " Please complete Step 2 (Imputation) first.")
      }
    })

    output$main_ui <- renderUI({
      req(isTRUE(app_state$imputation_done))
      tagList(
        fluidRow(
          column(3,
            wellPanel(
              h4("Settings"),
              sliderInput(ns("top_n"), "Top N variable proteins for MDS/QC",
                          min = 100, max = 5000, value = 2000, step = 100),
              checkboxInput(ns("z_score"), "Z-score proteins before PCA", value = TRUE),
              checkboxInput(ns("show_labels"), "Show sample labels on plots", value = TRUE),
              checkboxInput(ns("run_umap"), "Run UMAP (requires 'uwot')", value = FALSE),
              hr(),
              actionButton(ns("btn_run"), "Run Analysis",
                           class = "btn-primary btn-block"),
              br(),
              uiOutput(ns("run_status")),
              hr(),
              downloadButton(ns("dl_stats"), "Download QC Stats CSV",
                             class = "btn-default btn-sm btn-block"),
              downloadButton(ns("dl_all_exploration"), "Download All Exploration ZIP",
                             class = "btn-default btn-sm btn-block")
            )
          ),
          column(9,
            tabsetPanel(
              tabPanel("PCA",
                br(),
                fluidRow(
                  column(8, plotOutput(ns("pca_plot"), height = "420px")),
                  column(4, plotOutput(ns("scree_plot"), height = "420px"))
                ),
                downloadButton(ns("dl_pca_svg"), "Download PCA SVG",
                               class = "btn-default btn-sm")
              ),
              tabPanel("MDS",
                br(),
                plotOutput(ns("mds_plot"), height = "420px"),
                downloadButton(ns("dl_mds_svg"), "Download MDS SVG",
                               class = "btn-default btn-sm")
              ),
              tabPanel("Correlation Heatmap",
                br(),
                plotOutput(ns("cor_heatmap"), height = "480px"),
                downloadButton(ns("dl_cor_svg"), "Download Heatmap SVG",
                               class = "btn-default btn-sm")
              ),
              tabPanel("Sample Dendrogram",
                br(),
                plotOutput(ns("dendrogram"), height = "420px"),
                downloadButton(ns("dl_dendrogram_svg"), "Download Dendrogram SVG",
                               class = "btn-default btn-sm")
              ),
              tabPanel("UMAP",
                br(),
                uiOutput(ns("umap_msg")),
                plotOutput(ns("umap_plot"), height = "420px"),
                downloadButton(ns("dl_umap_svg"), "Download UMAP SVG",
                               class = "btn-default btn-sm")
              ),
              tabPanel("QC / Outlier Table",
                br(),
                uiOutput(ns("outlier_info")),
                br(),
                DT::dataTableOutput(ns("qc_table"))
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
    })
    observeEvent(input$btn_run, {
      req(isTRUE(app_state$imputation_done))
      d <- app_state$imputed_data

      withProgress(message = "Running exploratory analysis...", value = 0, {

        incProgress(0.05, detail = "Preparing matrix")
        int       <- as.matrix(d$intensity)
        int[!is.finite(int)] <- NA

        # PCA uses the full matrix; the other QC views use variable proteins so
        # sample structure is not drowned out by near-constant measurements.
        X_all <- int

        topN      <- min(input$top_n, nrow(int))
        prot_var  <- apply(int, 1, var, na.rm = TRUE)
        top_idx   <- order(prot_var, decreasing = TRUE)[seq_len(topN)]
        X_top     <- int[top_idx, , drop = FALSE]

        # Z-score before PCA so high-abundance proteins do not dominate.
        if (input$z_score) {
          Xz_all <- t(scale(t(X_all), center = TRUE, scale = TRUE))
          Xz_all[!is.finite(Xz_all)] <- 0
        } else {
          Xz_all <- X_all
          Xz_all[is.na(Xz_all)] <- 0
        }

        # Use the same scaling choice for the reduced QC matrix.
        if (input$z_score) {
          Xz <- t(scale(t(X_top), center = TRUE, scale = TRUE))
          Xz[!is.finite(Xz)] <- 0
        } else {
          Xz <- X_top
          Xz[is.na(Xz)] <- 0
        }

        sample_df <- data.frame(
          sample = colnames(Xz_all),
          group  = infer_groups(colnames(Xz_all)),
          stringsAsFactors = FALSE
        )
        incProgress(0.2, detail = "PCA from all proteins")
        pca    <- prcomp(t(Xz_all), center = FALSE, scale. = FALSE)
        pve    <- (pca$sdev^2) / sum(pca$sdev^2)
        pca_df <- data.frame(
          sample = sample_df$sample,
          group  = sample_df$group,
          PC1    = pca$x[, 1],
          PC2    = pca$x[, 2],
          stringsAsFactors = FALSE
        )
        incProgress(0.15, detail = "MDS")
        mds_res <- tryCatch(
          limma::plotMDS(Xz, plot = FALSE),
          error = function(e) NULL
        )
        mds_df <- if (!is.null(mds_res)) {
          data.frame(sample = sample_df$sample, group = sample_df$group,
                     Dim1 = mds_res$x, Dim2 = mds_res$y, stringsAsFactors = FALSE)
        } else NULL
        incProgress(0.15, detail = "Correlation matrix")
        cor_mat  <- cor(Xz, use = "pairwise.complete.obs", method = "pearson")
        dist_mat <- as.dist(1 - cor_mat)
        hc       <- hclust(dist_mat, method = "average")
        incProgress(0.1, detail = "QC stats")
        grp_map <- setNames(sample_df$group, sample_df$sample)
        samples <- colnames(cor_mat)
        stats_rows <- lapply(samples, function(s) {
          g       <- grp_map[[s]]
          within  <- cor_mat[s, names(grp_map)[grp_map == g]]
          within  <- within[names(within) != s]
          between <- cor_mat[s, names(grp_map)[grp_map != g]]
          data.frame(
            Sample         = s,
            Group          = g,
            Within_Median  = round(median(within,  na.rm = TRUE), 4),
            Within_Mean    = round(mean(within,    na.rm = TRUE), 4),
            Within_Min     = round(min(within,     na.rm = TRUE), 4),
            Between_Max    = round(max(between,    na.rm = TRUE), 4),
            Sep_Score      = round(median(within,  na.rm = TRUE) -
                                     max(between,  na.rm = TRUE), 4),
            stringsAsFactors = FALSE
          )
        })
        sample_stats <- do.call(rbind, stats_rows)
        # Flag samples that have unusually weak correlation to their own group.
        sample_stats <- do.call(rbind, lapply(split(sample_stats, sample_stats$Group), function(sub) {
          med <- median(sub$Within_Median, na.rm = TRUE)
          m   <- mad(sub$Within_Median, constant = 1, na.rm = TRUE)
          sub$Robust_Z   <- ifelse(m > 0, (sub$Within_Median - med) / m, NA)
          sub$Outlier    <- !is.na(sub$Robust_Z) & sub$Robust_Z < -3
          sub
        }))
        incProgress(0.1, detail = "UMAP")
        umap_df <- NULL
        if (isTRUE(input$run_umap) && requireNamespace("uwot", quietly = TRUE)) {
          set.seed(42)
          um <- tryCatch(
            uwot::umap(t(Xz), n_neighbors = min(10, ncol(Xz) - 1),
                       min_dist = 0.3, metric = "cosine"),
            error = function(e) NULL
          )
          if (!is.null(um)) {
            umap_df <- data.frame(sample = sample_df$sample, group = sample_df$group,
                                  UMAP1 = um[, 1], UMAP2 = um[, 2], stringsAsFactors = FALSE)
          }
        }

        incProgress(0.1, detail = "Storing results")
        explore_results(list(
          pca_df       = pca_df,
          pve          = pve,
          mds_df       = mds_df,
          cor_mat      = cor_mat,
          hc           = hc,
          sample_stats = sample_stats,
          umap_df      = umap_df,
          show_labels  = input$show_labels
        ))
        app_state$exploration_done <- TRUE
      })
    })

    output$run_status <- renderUI({
      if (is.null(explore_results())) return(NULL)
      tags$div(class = "alert alert-success",
               " Analysis complete")
    })
    # PCA
    pastel_group_palette <- function(groups) {
      n <- length(groups)
      if (n <= 0) return(character(0))

      base <- c(
        RColorBrewer::brewer.pal(8, "Pastel2"),
        RColorBrewer::brewer.pal(9, "Pastel1"),
        RColorBrewer::brewer.pal(12, "Set3")
      )

      if (n <= length(base)) {
        base[seq_len(n)]
      } else {
        grDevices::hcl.colors(n, palette = "Pastel 1")
      }
    }

    draw_pca_plot <- function(res, equal_axes = FALSE) {
      df  <- res$pca_df
      pve <- res$pve
      
      groups <- sort(unique(df$group))
      pal <- pastel_group_palette(groups)
      cols <- pal[match(df$group, groups)]
      
      # Leave room for group labels when many conditions are present.
      par(mar = c(5, 5, 3, 8), xpd = TRUE)
      
      plot_args <- list(
        x = df$PC1,
        y = df$PC2,
        bg = cols,
        col = "#6B7280",
        pch = 21,
        cex = 2,
        lwd = 1.2,
        xlab = sprintf("PC1 (%.1f%%)", 100 * pve[1]),
        ylab = sprintf("PC2 (%.1f%%)", 100 * pve[2]),
        main = "PCA - Samples (all proteins)"
      )
      if (isTRUE(equal_axes)) plot_args$asp <- 1
      do.call(plot, plot_args)
      
      # Keep the legend off the points.
      legend(x = par("usr")[2], y = par("usr")[4],
             legend = groups, pt.bg = pal, col = "#6B7280", pch = 21,
             bty = "n", cex = 0.9, xjust = 0)
    }

    output$pca_plot <- renderPlot({
      req(explore_results())
      draw_pca_plot(explore_results())
    })

    draw_scree_plot <- function(res) {
      pve <- res$pve
      n_pc <- min(10, length(pve))
      bar_cols <- grDevices::hcl.colors(n_pc, palette = "Pastel 1")

      par(mar = c(5, 5, 3, 1))
      barplot(
        100 * pve[seq_len(n_pc)],
        names.arg = paste0("PC", seq_len(n_pc)),
        col = bar_cols,
        border = "white",
        ylab = "Variance explained (%)",
        main = "PCA Scree Plot"
      )
    }

    output$scree_plot <- renderPlot({
      req(explore_results())
      draw_scree_plot(explore_results())
    })
    # MDS
    draw_mds_plot <- function(res) {
      df  <- res$mds_df
      if (is.null(df)) {
        plot.new(); title("MDS not available (limma error)")
        return()
      }
      
      groups <- sort(unique(df$group))
      pal <- pastel_group_palette(groups)
      cols <- pal[match(df$group, groups)]
      
      # Leave room for group labels when many conditions are present.
      par(mar = c(5, 5, 3, 8), xpd = TRUE)
      
      plot(df$Dim1, df$Dim2,
           bg = cols, col = "#6B7280", pch = 21, cex = 2, lwd = 1.2,
           xlab = "Dimension 1", ylab = "Dimension 2",
           main = "MDS (limma)")
      
      # Keep the legend off the points.
      legend(x = par("usr")[2], y = par("usr")[4],
             legend = groups, pt.bg = pal, col = "#6B7280", pch = 21,
             bty = "n", cex = 0.9, xjust = 0)
    }

    output$mds_plot <- renderPlot({
      req(explore_results())
      draw_mds_plot(explore_results())
    })
    # Correlation heatmap
    output$cor_heatmap <- renderPlot({
      req(explore_results())
      cor_mat <- explore_results()$cor_mat
      n       <- nrow(cor_mat)
      par(mar = c(8, 8, 3, 2))
      image(cor_mat[nrow(cor_mat):1, ],
            col = colorRampPalette(c("navy", "white", "firebrick"))(100),
            xaxt = "n", yaxt = "n",
            main = "Sample-Sample Pearson Correlation")
      axis(1, at = seq(0, 1, length.out = n), labels = colnames(cor_mat), las = 2, cex.axis = 0.75)
      axis(2, at = seq(0, 1, length.out = n), labels = rev(rownames(cor_mat)), las = 2, cex.axis = 0.75)
    })
    # Dendrogram
    output$dendrogram <- renderPlot({
      req(explore_results())
      hc <- explore_results()$hc
      par(mar = c(5, 4, 3, 1))
      plot(hc, main = "Hierarchical Clustering of Samples\n(1 - Pearson r)",
           xlab = "", sub = "", cex = 0.85)
    })
    # UMAP
    output$umap_msg <- renderUI({
      req(explore_results())
      df <- explore_results()$umap_df
      if (is.null(df)) {
        if (!requireNamespace("uwot", quietly = TRUE)) {
          tags$div(class = "alert alert-info",
                   "Package 'uwot' not installed. Run ",
                   tags$code("install.packages('uwot')"),
                   " and enable UMAP in settings.")
        } else {
          tags$div(class = "alert alert-info",
                   "Enable 'Run UMAP' in the settings panel and re-run analysis.")
        }
      }
    })

    draw_umap_plot <- function(res) {
      df <- res$umap_df
      req(!is.null(df))
      groups <- sort(unique(df$group))
      pal    <- pastel_group_palette(groups)
      cols   <- pal[match(df$group, groups)]
      par(mar = c(5, 5, 3, 2))
      plot(df$UMAP1, df$UMAP2,
           bg = cols, col = "#6B7280", pch = 21, cex = 1.8, lwd = 1.2,
           xlab = "UMAP 1", ylab = "UMAP 2",
           main = "UMAP")
      if (res$show_labels)
        text(df$UMAP1, df$UMAP2, labels = df$sample, pos = 3, cex = 0.7)
      legend("topright", legend = groups, pt.bg = pal, col = "#6B7280",
             pch = 21, bty = "n", cex = 0.85)
    }

    output$umap_plot <- renderPlot({
      req(explore_results())
      draw_umap_plot(explore_results())
    })

    save_svg_plot <- function(file, draw_fun, width = 8, height = 6) {
      grDevices::svg(file, width = width, height = height)
      on.exit(grDevices::dev.off(), add = TRUE)
      draw_fun()
    }

    output$dl_pca_svg <- downloadHandler(
      filename = function() download_filename(app_state, "exploration_pca", "svg"),
      content = function(file) {
        req(explore_results())
        save_svg_plot(file, function() {
          par(mfrow = c(1, 2))
          draw_pca_plot(explore_results(), equal_axes = TRUE)
          draw_scree_plot(explore_results())
        }, width = 12, height = 6)
      }
    )

    output$dl_mds_svg <- downloadHandler(
      filename = function() download_filename(app_state, "exploration_mds", "svg"),
      content = function(file) {
        req(explore_results())
        save_svg_plot(file, function() draw_mds_plot(explore_results()))
      }
    )

    output$dl_cor_svg <- downloadHandler(
      filename = function() download_filename(app_state, "exploration_correlation_heatmap", "svg"),
      content = function(file) {
        req(explore_results())
        save_svg_plot(file, function() {
          cor_mat <- explore_results()$cor_mat
          n <- nrow(cor_mat)
          par(mar = c(8, 8, 3, 2))
          image(cor_mat[nrow(cor_mat):1, ],
                col = colorRampPalette(c("navy", "white", "firebrick"))(100),
                xaxt = "n", yaxt = "n",
                main = "Sample-Sample Pearson Correlation")
          axis(1, at = seq(0, 1, length.out = n), labels = colnames(cor_mat), las = 2, cex.axis = 0.75)
          axis(2, at = seq(0, 1, length.out = n), labels = rev(rownames(cor_mat)), las = 2, cex.axis = 0.75)
        }, width = 8, height = 8)
      }
    )

    output$dl_dendrogram_svg <- downloadHandler(
      filename = function() download_filename(app_state, "exploration_dendrogram", "svg"),
      content = function(file) {
        req(explore_results())
        save_svg_plot(file, function() {
          hc <- explore_results()$hc
          par(mar = c(5, 4, 3, 1))
          plot(hc, main = "Hierarchical Clustering of Samples\n(1 - Pearson r)",
               xlab = "", sub = "", cex = 0.85)
        })
      }
    )

    output$dl_umap_svg <- downloadHandler(
      filename = function() download_filename(app_state, "exploration_umap", "svg"),
      content = function(file) {
        req(explore_results())
        save_svg_plot(file, function() draw_umap_plot(explore_results()))
      }
    )

    output$dl_all_exploration <- downloadHandler(
      filename = function() download_filename(app_state, "exploration_analysis", "zip"),
      content = function(file) {
        req(explore_results())
        res <- explore_results()
        tmp <- tempfile("exploration_export_")
        dir.create(tmp, recursive = TRUE, showWarnings = FALSE)

        files <- character(0)

        qc_file <- file.path(tmp, "QC_sample_stats.csv")
        utils::write.csv(res$sample_stats, qc_file, row.names = FALSE)
        files <- c(files, qc_file)

        pca_file <- file.path(tmp, "PCA_and_scree.svg")
        save_svg_plot(pca_file, function() {
          par(mfrow = c(1, 2))
          draw_pca_plot(res, equal_axes = TRUE)
          draw_scree_plot(res)
        }, width = 12, height = 6)
        files <- c(files, pca_file)

        mds_file <- file.path(tmp, "MDS.svg")
        save_svg_plot(mds_file, function() draw_mds_plot(res))
        files <- c(files, mds_file)

        cor_file <- file.path(tmp, "Correlation_heatmap.svg")
        save_svg_plot(cor_file, function() {
          cor_mat <- res$cor_mat
          n <- nrow(cor_mat)
          par(mar = c(8, 8, 3, 2))
          image(cor_mat[nrow(cor_mat):1, ],
                col = colorRampPalette(c("navy", "white", "firebrick"))(100),
                xaxt = "n", yaxt = "n",
                main = "Sample-Sample Pearson Correlation")
          axis(1, at = seq(0, 1, length.out = n), labels = colnames(cor_mat), las = 2, cex.axis = 0.75)
          axis(2, at = seq(0, 1, length.out = n), labels = rev(rownames(cor_mat)), las = 2, cex.axis = 0.75)
        }, width = 8, height = 8)
        files <- c(files, cor_file)

        dendrogram_file <- file.path(tmp, "Sample_dendrogram.svg")
        save_svg_plot(dendrogram_file, function() {
          hc <- res$hc
          par(mar = c(5, 4, 3, 1))
          plot(hc, main = "Hierarchical Clustering of Samples\n(1 - Pearson r)",
               xlab = "", sub = "", cex = 0.85)
        })
        files <- c(files, dendrogram_file)

        if (!is.null(res$umap_df)) {
          umap_file <- file.path(tmp, "UMAP.svg")
          save_svg_plot(umap_file, function() draw_umap_plot(res))
          files <- c(files, umap_file)
        }

        zip::zip(zipfile = file, files = files, mode = "cherry-pick")
      }
    )
    # QC table
    output$qc_table <- DT::renderDataTable({
      req(explore_results())
      df <- explore_results()$sample_stats
      
      # The table is meant as a triage view; full stats are available by download.
      df_outliers <- df[df$Outlier == TRUE, ]
      
      DT::datatable(df_outliers, rownames = FALSE,
                    options = list(pageLength = 20, dom = "ftip")) |>
        DT::formatStyle("Outlier",
                        backgroundColor = DT::styleEqual(TRUE, "rgba(255,100,100,0.3)"))
    })
    output$dl_stats <- downloadHandler(
      filename = function() download_filename(app_state, "QC_sample_stats", "csv"),
      content  = function(file) {
        req(explore_results())
        utils::write.csv(explore_results()$sample_stats, file, row.names = FALSE)
      }
    )
    output$proceed_btn <- renderUI({
      if (!isTRUE(app_state$exploration_done)) return(NULL)
      actionButton(ns("btn_proceed"), "Proceed to DEqMS",
                   class = "btn-success btn-lg app-proceed-btn")
    })

    observeEvent(input$btn_proceed, {
      app_state$active_tab <- "deqms"
    })
  })
}
