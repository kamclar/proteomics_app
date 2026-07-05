# Stable/background protein overview across DEqMS comparisons

mod_background_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Step 7. Stable Proteins"),
    uiOutput(ns("guard")),
    uiOutput(ns("main_ui"))
  )
}

mod_background_server <- function(id, app_state) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    reactome_cache <- reactiveVal(read_reactome_cache())

    output$guard <- renderUI({
      if (!isTRUE(app_state$deqms_done)) {
        tags$div(class = "alert alert-warning",
                 " Please complete Step 4 (DEqMS) first.")
      }
    })

    output$main_ui <- renderUI({
      req(isTRUE(app_state$deqms_done))
      nms <- names(app_state$deqms_results)
      tagList(
        fluidRow(
          column(
            3,
            wellPanel(
              h4("Comparisons"),
              div(style = "display:flex; gap:8px; margin-bottom:6px;",
                actionButton(ns("sel_all"), "All", class = "btn-xs btn-default"),
                actionButton(ns("desel_all"), "None", class = "btn-xs btn-default")
              ),
              checkboxGroupInput(ns("pairs"), NULL, choices = nms, selected = nms),
              hr(),
              h4("Definition"),
              numericInput(ns("changed_fc"), "Changed log2FC cutoff", value = 0.5, min = 0, step = 0.1),
              numericInput(ns("changed_p"), "Changed adj. p-value cutoff", value = 0.05, min = 0.001, max = 0.5, step = 0.005),
              numericInput(ns("stable_fc"), "Stable max |log2FC|", value = 0.5, min = 0, step = 0.1),
              numericInput(ns("min_valid"), "Stable min valid samples", value = 0.75, min = 0, max = 1, step = 0.05),
              helpText("Stable candidates are non-significant in all selected comparisons, have small absolute fold change, enough measured values, and no on/off pattern."),
              hr(),
              h4("Reactome"),
              actionButton(ns("btn_reactome_cache"), "Build Reactome cache", class = "btn-sm btn-default btn-block"),
              helpText("Adds currently unseen proteins to the local Reactome cache, then compares pathway membership for stable and changed proteins."),
              uiOutput(ns("reactome_cache_status"))
            )
          ),
          column(
            9,
            fluidRow(
              column(12, uiOutput(ns("definition_note"))),
              column(12, uiOutput(ns("summary_boxes"))),
              column(12, plotOutput(ns("class_plot"), height = "320px")),
              column(12, plotOutput(ns("classification_map"), height = "430px")),
              column(12, plotOutput(ns("reactome_plot"), height = "430px")),
              column(
                12,
                tabsetPanel(
                  tabPanel("Stable candidates", br(), DT::dataTableOutput(ns("stable_table"))),
                  tabPanel("Changed proteins", br(), DT::dataTableOutput(ns("changed_table"))),
                  tabPanel("Other tested", br(), DT::dataTableOutput(ns("other_table"))),
                  tabPanel("Reactome comparison", br(), DT::dataTableOutput(ns("reactome_table")))
                )
              )
            )
          )
        )
      )
    })

    output$definition_note <- renderUI({
      tags$div(
        class = "alert alert-info",
        tags$b("Classification: "),
        "Changed = significant in at least one selected comparison. ",
        "Stable candidate = non-significant in all selected comparisons, small maximum absolute log2FC, sufficient measured values, and no on/off pattern. ",
        "Other tested = measured/tested proteins that are not confidently changed or stable."
      )
    })

    observeEvent(input$sel_all, {
      req(isTRUE(app_state$deqms_done))
      updateCheckboxGroupInput(session, "pairs", selected = names(app_state$deqms_results))
    })

    observeEvent(input$desel_all, {
      updateCheckboxGroupInput(session, "pairs", selected = character(0))
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
        sprintf("Local cache: %d mapped gene(s), %d pathway(s).", st$mapped_genes, st$pathways))
    })

    background_data <- reactive({
      req(isTRUE(app_state$deqms_done), length(input$pairs) > 0)
      withProgress(message = "Preparing Stable Proteins data...", value = 0, {
        incProgress(0.25, detail = "Reading selected DEqMS comparisons")
        dat <- build_stable_background(
          app_state = app_state,
          pairs = input$pairs,
          changed_fc = input$changed_fc %||% 0.5,
          changed_p = input$changed_p %||% 0.05,
          stable_fc = input$stable_fc %||% 0.5,
          min_valid = input$min_valid %||% 0.75
        )
        incProgress(0.55, detail = "Classifying stable, changed, and other proteins")
        incProgress(0.20, detail = "Preparing plots and tables")
        dat
      })
    })

    output$summary_boxes <- renderUI({
      dat <- background_data()
      cls <- dat$summary
      box <- function(label, value, color) {
        tags$div(
          style = paste0("display:inline-block; min-width:150px; margin:0 10px 12px 0; padding:14px 16px; ",
                         "border-radius:6px; border:1px solid #d9e2ec; background:#fff; ",
                         "border-left:5px solid ", color, ";"),
          tags$div(style = "font-size:12px; color:#486581; font-weight:700;", label),
          tags$div(style = "font-size:26px; font-weight:700; color:#102a43;", value)
        )
      }
      tagList(
        box("Stable candidates", cls$Stable, "#2f855a"),
        box("Changed proteins", cls$Changed, "#c53030"),
        box("Other tested", cls$Other, "#718096"),
        box("Selected comparisons", length(input$pairs), "#3182ce")
      )
    })

    output$class_plot <- renderPlot({
      dat <- background_data()
      df <- dat$protein_table
      plot_df <- data.frame(Class = factor(df$Class, levels = c("Stable", "Changed", "Other")))
      ggplot2::ggplot(plot_df, ggplot2::aes(x = Class, fill = Class)) +
        ggplot2::geom_bar(width = 0.62) +
        ggplot2::scale_fill_manual(values = c(Stable = "#2f855a", Changed = "#c53030", Other = "#a0aec0")) +
        ggplot2::labs(x = NULL, y = "Proteins", title = "Stable candidates vs changed proteins") +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(legend.position = "none", plot.title = ggplot2::element_text(face = "bold"))
    })

    output$classification_map <- renderPlot({
      dat <- background_data()
      df <- dat$protein_table
      if (!nrow(df)) {
        plot.new()
        text(0.5, 0.5, "No proteins available for classification.", cex = 1.1)
        return()
      }
      df$Class <- factor(df$Class, levels = c("Stable", "Changed", "Other"))
      df$Plot_negLog10P <- ifelse(is.finite(df$Min_adj_pvalue) & df$Min_adj_pvalue > 0,
                                  -log10(df$Min_adj_pvalue), NA_real_)
      y_cap <- stats::quantile(df$Plot_negLog10P, 0.99, na.rm = TRUE, names = FALSE)
      if (!is.finite(y_cap)) y_cap <- 10
      df$Plot_negLog10P <- pmin(df$Plot_negLog10P, y_cap)
      changed_p <- input$changed_p %||% 0.05
      stable_fc <- input$stable_fc %||% 0.5
      changed_fc <- input$changed_fc %||% 0.5
      p_cut_line <- if (is.finite(changed_p) && changed_p > 0) -log10(changed_p) else NA_real_
      ggplot2::ggplot(
        df,
        ggplot2::aes(x = Max_abs_logFC, y = Plot_negLog10P, color = Class)
      ) +
        ggplot2::geom_point(ggplot2::aes(alpha = Valid_fraction), size = 1.8, na.rm = TRUE) +
        ggplot2::geom_vline(xintercept = stable_fc, linetype = "dashed", color = "#4a5568") +
        ggplot2::geom_vline(xintercept = changed_fc, linetype = "dotted", color = "#4a5568") +
        {
          if (is.finite(p_cut_line)) ggplot2::geom_hline(yintercept = p_cut_line, linetype = "dotted", color = "#4a5568")
        } +
        ggplot2::scale_color_manual(values = c(Stable = "#2f855a", Changed = "#c53030", Other = "#718096")) +
        ggplot2::scale_alpha_continuous(range = c(0.25, 0.9), limits = c(0, 1), na.value = 0.25) +
        ggplot2::labs(
          x = "Maximum absolute log2FC across selected comparisons",
          y = "-log10 minimum adjusted p-value",
          alpha = "Valid fraction",
          title = "Classification map"
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"), legend.title = ggplot2::element_text(face = "bold"))
    })

    output$reactome_plot <- renderPlot({
      dat <- background_data()
      rt <- build_reactome_set_comparison(dat$protein_table, reactome_cache())
      if (!nrow(rt)) {
        plot.new()
        text(0.5, 0.5, "No Reactome cache matches yet. Build the Reactome cache first.", cex = 1.1)
        return()
      }
      top <- rt[order(rt$p_adj, -rt$Set_Count), , drop = FALSE]
      top <- head(top, 16)
      top$Pathway <- factor(top$Pathway, levels = rev(unique(top$Pathway)))
      ggplot2::ggplot(top, ggplot2::aes(x = pmin(-log10(p_adj), 12), y = Pathway, fill = Class)) +
        ggplot2::geom_col(position = "dodge", width = 0.72) +
        ggplot2::scale_fill_manual(values = c(Stable = "#2f855a", Changed = "#c53030")) +
        ggplot2::labs(x = "-log10 adjusted Fisher p-value", y = NULL, title = "Reactome pathways: stable vs changed") +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"), legend.title = ggplot2::element_blank())
    })

    output$stable_table <- DT::renderDataTable({
      dat <- background_data()
      df <- dat$protein_table[dat$protein_table$Class == "Stable", , drop = FALSE]
      df <- df[order(df$Max_abs_logFC, -df$Valid_fraction), , drop = FALSE]
      show <- intersect(c("Gene", "Protein", "Description", "Max_abs_logFC", "Min_adj_pvalue",
                          "Valid_fraction", "Valid_samples", "Total_samples"), names(df))
      DT::datatable(df[, show, drop = FALSE], rownames = FALSE,
                    options = list(pageLength = 15, scrollX = TRUE, dom = "ftip")) |>
        DT::formatRound(intersect(c("Max_abs_logFC", "Min_adj_pvalue", "Valid_fraction"), show), 4)
    })

    output$changed_table <- DT::renderDataTable({
      dat <- background_data()
      df <- dat$protein_table[dat$protein_table$Class == "Changed", , drop = FALSE]
      df <- df[order(df$Min_adj_pvalue, -df$Max_abs_logFC), , drop = FALSE]
      show <- intersect(c("Gene", "Protein", "Description", "Changed_in", "Max_abs_logFC",
                          "Min_adj_pvalue", "Valid_fraction", "On_off_any"), names(df))
      DT::datatable(df[, show, drop = FALSE], rownames = FALSE,
                    options = list(pageLength = 15, scrollX = TRUE, dom = "ftip")) |>
        DT::formatRound(intersect(c("Max_abs_logFC", "Min_adj_pvalue", "Valid_fraction"), show), 4)
    })

    output$other_table <- DT::renderDataTable({
      dat <- background_data()
      df <- dat$protein_table[dat$protein_table$Class == "Other", , drop = FALSE]
      df <- df[order(df$Other_reason, df$Min_adj_pvalue, -df$Max_abs_logFC), , drop = FALSE]
      show <- intersect(c("Gene", "Protein", "Description", "Other_reason", "Max_abs_logFC",
                          "Min_adj_pvalue", "Valid_fraction", "Valid_samples",
                          "Total_samples", "On_off_any"), names(df))
      DT::datatable(df[, show, drop = FALSE], rownames = FALSE,
                    options = list(pageLength = 15, scrollX = TRUE, dom = "ftip")) |>
        DT::formatRound(intersect(c("Max_abs_logFC", "Min_adj_pvalue", "Valid_fraction"), show), 4)
    })

    output$reactome_table <- DT::renderDataTable({
      dat <- background_data()
      rt <- build_reactome_set_comparison(dat$protein_table, reactome_cache())
      if (!nrow(rt)) return(data.frame(Message = "No Reactome cache matches yet. Build the Reactome cache first."))
      rt <- rt[order(rt$p_adj, rt$Class), , drop = FALSE]
      DT::datatable(rt, rownames = FALSE,
                    options = list(pageLength = 15, scrollX = TRUE, dom = "ftip")) |>
        DT::formatRound(c("p_value", "p_adj"), 4)
    })
  })
}

build_stable_background <- function(app_state, pairs, changed_fc = 0.5, changed_p = 0.05,
                                    stable_fc = 0.5, min_valid = 0.75) {
  res <- app_state$deqms_results[pairs]
  res <- res[!vapply(res, is.null, logical(1))]
  if (!length(res)) return(list(protein_table = data.frame(), summary = list(Stable = 0, Changed = 0, Other = 0)))

  common <- Reduce(intersect, lapply(res, rownames))
  common <- common[!is.na(common) & nzchar(common)]
  if (!length(common)) return(list(protein_table = data.frame(), summary = list(Stable = 0, Changed = 0, Other = 0)))

  logfc <- sapply(res, function(df) as.numeric(df[common, "logFC"]))
  pval <- sapply(res, function(df) {
    if ("sca.adj.pval" %in% names(df)) as.numeric(df[common, "sca.adj.pval"]) else rep(NA_real_, length(common))
  })
  if (is.null(dim(logfc))) logfc <- matrix(logfc, ncol = 1)
  if (is.null(dim(pval))) pval <- matrix(pval, ncol = 1)
  colnames(logfc) <- names(res)
  colnames(pval) <- names(res)

  first <- res[[1]][common, , drop = FALSE]
  gene <- if ("Genes" %in% names(first)) as.character(first$Genes) else rep(NA_character_, length(common))
  desc <- if ("First.Protein.Descriptions" %in% names(first)) {
    as.character(first$First.Protein.Descriptions)
  } else {
    rep("", length(common))
  }

  valid_fraction <- rep(NA_real_, length(common))
  valid_samples <- rep(NA_integer_, length(common))
  total_samples <- rep(NA_integer_, length(common))
  on_off_any <- rep(FALSE, length(common))
  if (isTRUE(app_state$upload_done) && !is.null(app_state$parsed_data$intensity)) {
    int <- as.matrix(app_state$parsed_data$intensity)
    matched <- match(common, rownames(int))
    ok <- !is.na(matched)
    valid_samples[ok] <- rowSums(!is.na(int[matched[ok], , drop = FALSE]))
    total_samples[ok] <- ncol(int)
    valid_fraction[ok] <- valid_samples[ok] / pmax(total_samples[ok], 1L)
    on_off_any[ok] <- stable_on_off_any(int[matched[ok], , drop = FALSE],
                                        app_state$parsed_data$sample_names, names(res))
  }

  sig <- is.finite(logfc) & is.finite(pval) & abs(logfc) >= changed_fc & pval <= changed_p
  changed_any <- rowSums(sig, na.rm = TRUE) > 0
  max_abs_logfc <- apply(abs(logfc), 1, function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) NA_real_ else max(x)
  })
  min_adj_p <- apply(pval, 1, function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) NA_real_ else min(x)
  })
  stable_ok <- rowSums(is.finite(logfc) & is.finite(pval), na.rm = TRUE) == ncol(logfc) &
    max_abs_logfc <= stable_fc &
    min_adj_p > changed_p &
    !on_off_any &
    !is.na(valid_fraction) & valid_fraction >= min_valid

  class <- ifelse(stable_ok, "Stable", ifelse(changed_any, "Changed", "Other"))
  changed_in <- apply(sig, 1, function(x) paste(names(res)[which(x)], collapse = "; "))
  finite_stats_count <- rowSums(is.finite(logfc) & is.finite(pval), na.rm = TRUE)
  other_reason <- stable_other_reason(
    class = class,
    finite_stats_count = finite_stats_count,
    total_stats_count = ncol(logfc),
    max_abs_logfc = max_abs_logfc,
    min_adj_p = min_adj_p,
    stable_fc = stable_fc,
    changed_p = changed_p,
    valid_fraction = valid_fraction,
    min_valid = min_valid,
    on_off_any = on_off_any
  )

  protein_table <- data.frame(
    Protein = common,
    Gene = gene,
    Description = desc,
    Class = class,
    Changed_in = changed_in,
    Max_abs_logFC = max_abs_logfc,
    Min_adj_pvalue = min_adj_p,
    Valid_fraction = valid_fraction,
    Valid_samples = valid_samples,
    Total_samples = total_samples,
    On_off_any = on_off_any,
    Other_reason = other_reason,
    stringsAsFactors = FALSE
  )
  protein_table$Changed_in[protein_table$Changed_in == ""] <- NA_character_
  summary <- as.list(setNames(as.integer(table(factor(protein_table$Class, levels = c("Stable", "Changed", "Other")))),
                              c("Stable", "Changed", "Other")))
  list(protein_table = protein_table, summary = summary)
}

stable_other_reason <- function(class, finite_stats_count, total_stats_count, max_abs_logfc,
                                min_adj_p, stable_fc, changed_p, valid_fraction,
                                min_valid, on_off_any) {
  reason <- rep(NA_character_, length(class))
  other <- class == "Other"
  reason[other & finite_stats_count < total_stats_count] <- "Missing DEqMS statistics in one or more selected comparisons"
  reason[other & is.na(reason) & on_off_any] <- "On/off pattern in at least one selected comparison"
  reason[other & is.na(reason) & (is.na(valid_fraction) | valid_fraction < min_valid)] <- "Low data coverage"
  reason[other & is.na(reason) & is.finite(max_abs_logfc) & max_abs_logfc > stable_fc &
           (!is.finite(min_adj_p) | min_adj_p > changed_p)] <- "Large fold change but not significant"
  reason[other & is.na(reason) & is.finite(min_adj_p) & min_adj_p <= changed_p &
           (!is.finite(max_abs_logfc) | max_abs_logfc < stable_fc)] <- "Significant p-value but below fold-change cutoff"
  reason[other & is.na(reason)] <- "Borderline or mixed evidence"
  reason
}

stable_on_off_any <- function(int_mat, sample_names, pairs) {
  grps <- infer_groups(sample_names)
  out <- rep(FALSE, nrow(int_mat))
  for (pair in pairs) {
    parts <- strsplit(pair, "_vs_", fixed = TRUE)[[1]]
    if (length(parts) != 2) next
    c1 <- intersect(sample_names[grps == parts[1]], colnames(int_mat))
    c2 <- intersect(sample_names[grps == parts[2]], colnames(int_mat))
    if (!length(c1) || !length(c2)) next
    n1 <- rowSums(!is.na(int_mat[, c1, drop = FALSE]))
    n2 <- rowSums(!is.na(int_mat[, c2, drop = FALSE]))
    out <- out | ((n1 > 0 & n2 == 0) | (n1 == 0 & n2 > 0))
  }
  out
}

build_reactome_set_comparison <- function(protein_table, cache = read_reactome_cache()) {
  if (!nrow(protein_table) || !nrow(cache)) return(data.frame())
  mapped <- cache[!is.na(cache$pathname) & nzchar(cache$pathname), , drop = FALSE]
  if (!nrow(mapped)) return(data.frame())

  universe <- unique(split_gene_symbols(protein_table$Gene))
  mapped <- mapped[mapped$gene_symbol %in% universe, , drop = FALSE]
  if (!nrow(mapped)) return(data.frame())
  universe_mapped <- unique(mapped$gene_symbol)

  one_set <- function(class_name) {
    genes <- unique(split_gene_symbols(protein_table$Gene[protein_table$Class == class_name]))
    genes <- intersect(genes, universe_mapped)
    if (length(genes) < 3) return(data.frame())
    pathways <- sort(unique(mapped$pathname))
    rows <- lapply(pathways, function(pathway) {
      path_genes <- unique(mapped$gene_symbol[mapped$pathname == pathway])
      a <- length(intersect(genes, path_genes))
      if (a < 2) return(NULL)
      b <- length(setdiff(genes, path_genes))
      c <- length(setdiff(path_genes, genes))
      d <- length(setdiff(universe_mapped, union(genes, path_genes)))
      p <- stats::fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")$p.value
      data.frame(
        Class = class_name,
        Pathway = pathway,
        Set_Count = a,
        Set_Size = length(genes),
        Background_Count = length(path_genes),
        Background_Size = length(universe_mapped),
        p_value = p,
        stringsAsFactors = FALSE
      )
    })
    df <- do.call(rbind, Filter(Negate(is.null), rows))
    if (is.null(df) || !nrow(df)) return(data.frame())
    df$p_adj <- stats::p.adjust(df$p_value, method = "BH")
    df
  }

  out <- rbind(one_set("Stable"), one_set("Changed"))
  if (is.null(out)) data.frame() else out
}
