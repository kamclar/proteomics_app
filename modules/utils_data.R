# utils_data.R - shared data-wrangling helpers
META_COLS_BASE <- c(
  # Standard Spectronaut columns
  "Peptides", "Genes", "Protein.Ids",
  "First.Protein.Descriptions", "Protein.Accessions", "Protein.Groups",
  # Variations found in some exports
  "ProteinAccessions", "ProteinDescriptions", "ProteinNames", 
  "NrOfStrippedSequencesIdentified", "Qvalue", "Cscore"
)
# Strip "PG." prefix from metadata column names
clean_meta_colnames <- function(df) {
  names(df) <- make_unique_colnames(sub("^PG\\.", "", names(df)))
  df
}
# Keep imported column names stable for UI selection and downstream indexing.
# openxlsx can return blanks or duplicate names from unusual header rows.
make_unique_colnames <- function(nms) {
  nms <- trimws(as.character(nms))
  blank <- is.na(nms) | nms == ""
  nms[blank] <- paste0("Column_", which(blank))
  make.unique(nms, sep = "_")
}
# Guess whether the real column names start on row 1 or row 2.
# Spectronaut files with merged group headers often have row 1 mostly empty.
guess_header_row <- function(path, sheet = 1) {
  raw_test <- openxlsx::read.xlsx(path, sheet = sheet, rows = 1:3,
                                  check.names = FALSE, colNames = FALSE)
  if (!nrow(raw_test)) return(1L)

  row1_blank_pct <- sum(is.na(raw_test[1, ]) | trimws(as.character(raw_test[1, ])) == "") /
    ncol(raw_test)
  if (row1_blank_pct > 0.5) 2L else 1L
}
# Read a compact sheet preview for the upload screen.
read_sheet_preview <- function(path, sheet = 1, n_rows = 25) {
  preview <- openxlsx::read.xlsx(path, sheet = sheet, rows = seq_len(n_rows),
                                 check.names = FALSE, colNames = FALSE)
  names(preview) <- paste0("Column_", seq_len(ncol(preview)))
  preview
}
# Read the selected sheet using a chosen header row.
read_proteomics_sheet <- function(path, sheet = 1, header_row = NULL) {
  if (is.null(header_row) || is.na(header_row)) {
    header_row <- guess_header_row(path, sheet)
  }

  raw <- openxlsx::read.xlsx(path, sheet = sheet, startRow = header_row,
                             check.names = FALSE)
  names(raw) <- make_unique_colnames(names(raw))
  raw
}
# Detect intensity columns: have numeric prefix and .raw.PG.Quantity suffix
# Patterns:
#   "26_20251118_BT_2.raw.PG.Quantity" -> "BT_2"
#   "[1] 20251118_B1_1.raw.PG.Quantity" -> "B1_1"
clean_intensity_colname_raw <- function(x) {
  x <- sub("^X\\.", "", x)                           # Remove R's X. prefix
  x <- sub("^\\[[0-9]+\\]\\.?", "", x)               # Remove [1], [2], etc.
  x <- sub("^[0-9]+[._]", "", x)                     # Remove any remaining numeric prefix
  x <- sub("^[0-9]{6,8}[._]", "", x)                 # Remove date prefix (6-8 digits)
  x <- sub("\\.raw.*$", "", x, ignore.case = TRUE)   # Remove .raw.PG.Quantity
  trimws(x)
}

clean_intensity_colname <- function(x) {
  make_valid_sample_name(clean_intensity_colname_raw(x))
}

build_name_change_log <- function(original_samples, safe_samples) {
  sample_changes <- data.frame(
    Type = "Sample",
    Original = original_samples,
    Renamed = safe_samples,
    stringsAsFactors = FALSE
  )
  sample_changes <- sample_changes[sample_changes$Original != sample_changes$Renamed, , drop = FALSE]

  original_groups <- unique(sub("_[0-9]+$", "", original_samples))
  safe_groups <- make_valid_group_name(original_groups)
  group_changes <- data.frame(
    Type = "Group",
    Original = original_groups,
    Renamed = safe_groups,
    stringsAsFactors = FALSE
  )
  group_changes <- group_changes[group_changes$Original != group_changes$Renamed, , drop = FALSE]

  rbind(sample_changes, group_changes)
}

is_intensity_col <- function(colname) {
  grepl("\\.raw\\.PG\\.Quantity", colname, ignore.case = TRUE)
}
# Group names are later used in limma design/contrast matrices, so they must be
# valid R names. Do this at import time so t-test and DEqMS keys stay aligned.
make_valid_group_name <- function(x) {
  make.names(x, unique = FALSE)
}

make_valid_sample_name <- function(x) {
  make.names(x, unique = FALSE)
}
# Extract experiment ID from Spectronaut-style raw/run names.
# Expected pattern: MS followed by digits and an underscore, e.g. MS12345_.
extract_experiment_name <- function(x = character(0), path = NULL) {
  candidates <- c(as.character(x), if (!is.null(path)) basename(path))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]

  if (!length(candidates)) return("experiment")

  patterns <- c("MS[0-9]+_", "MS[0-9]+")
  for (pattern in patterns) {
    hits <- regmatches(candidates, regexpr(pattern, candidates, ignore.case = TRUE))
    hits <- hits[nzchar(hits)]
    if (length(hits)) {
      exp_name <- sub("_$", "", hits[[1]])
      exp_name <- toupper(exp_name)
      return(make.names(exp_name, unique = FALSE))
    }
  }

  "experiment"
}

download_experiment_prefix <- function(app_state, fallback = "experiment") {
  exp_name <- app_state$parsed_data$experiment_name
  if (is.null(exp_name) || is.na(exp_name) || !nzchar(exp_name)) {
    exp_name <- fallback
  }
  exp_name
}

download_filename <- function(app_state, suffix, ext, fallback = "experiment") {
  paste0(download_experiment_prefix(app_state, fallback), "_", suffix, "_", Sys.Date(), ".", ext)
}
# Detect t-test result columns from the first sheet
# Returns list(diff=..., neglogp=..., fc=..., mean=..., median=...)
detect_ttest_cols <- function(df) {
  nms <- names(df)
  list(
    diff    = grep("Student.*T-test.*Difference",   nms, ignore.case = TRUE, value = TRUE),
    neglogp = grep("-Log.*T-test.*p-value",          nms, ignore.case = TRUE, value = TRUE),
    fc      = grep("fold change",                    nms, ignore.case = TRUE, value = TRUE),
    mean    = grep("^Mean\\s",                       nms, ignore.case = TRUE, value = TRUE),
    median  = grep("^Median\\s",                     nms, ignore.case = TRUE, value = TRUE)
  )
}
# Suggest column roles for the upload UI using the same rules as parsing.
detect_upload_columns <- function(df) {
  clean_names <- names(clean_meta_colnames(df))

  intensity <- clean_names[vapply(clean_names, is_intensity_col, logical(1))]
  metadata <- clean_names[
    clean_names %in% META_COLS_BASE |
      grepl("NrOfStrippedSequences|Peptides|PeptideCount", clean_names, ignore.case = TRUE)
  ]

  ttest_pattern <- "Student|fold change|^Mean |^Median |-Log"
  ttest <- grep(ttest_pattern, clean_names, ignore.case = TRUE, value = TRUE)

  list(
    metadata = metadata,
    intensity = intensity,
    ttest = ttest,
    all_columns = clean_names
  )
}
# Parse uploaded XLSX: returns a list with
#   $meta        - metadata data.frame (proteins x meta cols)
#   $intensity   - numeric matrix/df (proteins x samples), log2 or raw
#   $ttest       - list of t-test column vectors, keyed by comparison name
#   $sample_names - cleaned sample names
#   $raw_df      - the full cleaned data.frame (used downstream)
parse_proteomics_xlsx <- function(path, sheet = 1, already_log2 = TRUE,
                                  header_row = NULL, selected_columns = NULL,
                                  source_name = NULL) {

  raw <- read_proteomics_sheet(path, sheet = sheet, header_row = header_row)

  # Clean PG. prefix
  raw <- clean_meta_colnames(raw)

  col_nms <- names(raw)
  suggested <- detect_upload_columns(raw)

  selected_intensity <- suggested$intensity
  selected_metadata <- suggested$metadata
  selected_ttest <- suggested$ttest

  if (!is.null(selected_columns)) {
    if (!is.null(selected_columns$intensity)) {
      selected_intensity <- intersect(selected_columns$intensity, col_nms)
    }
    if (!is.null(selected_columns$metadata)) {
      selected_metadata <- intersect(selected_columns$metadata, col_nms)
    }
    if (!is.null(selected_columns$ttest)) {
      selected_ttest <- intersect(selected_columns$ttest, col_nms)
    }
  }

  int_idx <- match(selected_intensity, col_nms)
  meta_idx <- match(selected_metadata, col_nms)
  ttest_idx <- match(selected_ttest, col_nms)

  if (!length(int_idx)) {
    stop("No intensity columns selected. Select columns ending in .raw.PG.Quantity or choose them manually.")
  }

  int_cols_raw    <- col_nms[int_idx]
  experiment_name <- extract_experiment_name(
    int_cols_raw,
    path = if (!is.null(source_name) && nzchar(source_name)) source_name else path
  )
  int_cols_original_clean <- vapply(int_cols_raw, clean_intensity_colname_raw, character(1))
  int_cols_clean  <- make_valid_sample_name(int_cols_original_clean)
  name_changes <- build_name_change_log(int_cols_original_clean, int_cols_clean)

  df_clean <- raw
  names(df_clean)[int_idx] <- int_cols_clean

  meta_df <- df_clean[, col_nms[meta_idx], drop = FALSE]

  peptide_col_candidates <- grep(
    "NrOfStrippedSequences|Peptides|NrOfPeptides|PeptideCount",
    names(meta_df),
    ignore.case = TRUE,
    value = TRUE
  )
  
  if (length(peptide_col_candidates) > 0) {
    peptide_col <- peptide_col_candidates[1]
    if (peptide_col != "Peptides") {
      meta_df$Peptides <- as.numeric(meta_df[[peptide_col]])
    }
  } else {
    message("WARNING: No peptide count column found. DEqMS will require a valid peptide count column before it can run.")
    meta_df$Peptides <- NA_real_
  }

  
  # Use best available column as rownames
  id_col <- NULL
  if ("Protein.Ids" %in% names(meta_df)) {
    id_col <- "Protein.Ids"
  } else if ("ProteinAccessions" %in% names(meta_df)) {
    id_col <- "ProteinAccessions"
  } else if ("Genes" %in% names(meta_df)) {
    id_col <- "Genes"
  }
  
  if (!is.null(id_col)) {
    rownames(meta_df) <- make.unique(as.character(meta_df[[id_col]]))
  } else {
    rownames(meta_df) <- paste0("Protein_", seq_len(nrow(meta_df)))
  }

  # Intensity matrix (convert to numeric)
  int_df <- df_clean[, int_cols_clean, drop = FALSE]
  int_df[] <- lapply(int_df, function(x) suppressWarnings(as.numeric(as.character(x))))
  rownames(int_df) <- rownames(meta_df)

  if (!already_log2) {
    int_df[int_df <= 0 | is.na(int_df)] <- NA
    int_df[] <- lapply(int_df, function(x) log2(x))
  }

  ttest_cols_raw <- names(df_clean)[ttest_idx]
  groups_detected <- unique(infer_groups(int_cols_clean))

  ttest_list     <- parse_ttest_cols(raw, ttest_cols_raw, groups_detected)

  list(
    meta         = meta_df,
    intensity    = int_df,
    ttest        = ttest_list,
    sample_names = int_cols_clean,
    experiment_name = experiment_name,
    name_changes = name_changes,
    raw_df       = df_clean,
    was_log2     = already_log2
  )
}
# Extract t-test comparison names and build per-comparison data.frames
# Handles patterns like:
#   "Student's T-test Difference B1_BD"
#   "B1_BD fold change"
#   "-Log Student's T-test p-value B1_BD"
parse_ttest_cols <- function(df, ttest_col_names, known_groups = NULL) {
  if (!length(ttest_col_names)) return(list())
  
  diff_cols    <- grep("T[.-]test.*Difference", ttest_col_names, ignore.case = TRUE, value = TRUE)
  neglogp_cols <- grep("Log.*p[.-]value", ttest_col_names, ignore.case = TRUE, value = TRUE)
  fc_cols      <- grep("fold change",          ttest_col_names, ignore.case = TRUE, value = TRUE)

  comparisons <- sub(".*Difference\\.", "", diff_cols)
  comparisons <- trimws(comparisons)

  comparisons_original <- comparisons
  
  if (!is.null(known_groups)) {
    comparisons <- sapply(comparisons, function(comp) {
      parts <- strsplit(comp, "_")[[1]]
      
      for (i in 1:(length(parts)-1)) {
        g1_candidate <- paste(parts[1:i], collapse="_")
        g2_candidate <- paste(parts[(i+1):length(parts)], collapse="_")
        g1_safe <- make_valid_group_name(g1_candidate)
        g2_safe <- make_valid_group_name(g2_candidate)
        
        # Check if both are valid groups
        if (g1_safe %in% known_groups && g2_safe %in% known_groups) {
          return(paste0(g1_safe, "_vs_", g2_safe))
        }
      }
      
      message("Warning: Could not split '", comp, "' using known sample groups; skipping this t-test comparison.")
      NA_character_
    }, USE.NAMES = FALSE)
    keep <- !is.na(comparisons) & nzchar(comparisons)
    comparisons <- comparisons[keep]
    comparisons_original <- comparisons_original[keep]
  } else {
    comparisons <- gsub("^([^_]+)_(.+)$", "\\1_vs_\\2", comparisons)
    comparisons <- vapply(strsplit(comparisons, "_vs_"), function(parts) {
      paste(make_valid_group_name(parts), collapse = "_vs_")
    }, character(1))
  }

  result <- list()
  for (i in seq_along(comparisons)) {
    comp_vs <- comparisons[i]
    comp_orig <- comparisons_original[i]
    
    comp_esc <- gsub("([.])", "\\\\\\1", comp_orig)
    
    d_col  <- grep(paste0("Difference\\.", comp_esc, "$"),  diff_cols,    value = TRUE)[1]
    p_col  <- grep(paste0("p[.-]value.*", comp_esc, "$"),  neglogp_cols, value = TRUE)[1]
    fc_col <- grep(paste0("^", comp_esc, ".*fold.*change"),     fc_cols,      value = TRUE)[1]
    
    if (is.na(d_col) || is.na(p_col)) {
      message("Warning: Could not find columns for comparison '", comp_orig, "'")
      next
    }
    
    sub_df <- data.frame(
      logFC       = suppressWarnings(as.numeric(df[[d_col]])),
      negLog10p   = suppressWarnings(as.numeric(df[[p_col]])),
      stringsAsFactors = FALSE
    )
    if (!is.na(fc_col)) sub_df$foldChange <- suppressWarnings(as.numeric(df[[fc_col]]))
    
    # Carry protein IDs for tooltip - handle both naming conventions
    if ("Genes" %in% names(df)) {
      sub_df$Genes <- df$Genes
    }
    if ("First.Protein.Descriptions" %in% names(df)) {
      sub_df$ProtDesc <- df$First.Protein.Descriptions
    } else if ("ProteinDescriptions" %in% names(df)) {
      sub_df$ProtDesc <- df$ProteinDescriptions
    }
    if ("Protein.Accessions" %in% names(df)) {
      sub_df$Protein.Accessions <- df$Protein.Accessions
    } else if ("ProteinAccessions" %in% names(df)) {
      sub_df$Protein.Accessions <- df$ProteinAccessions
    }
    
    result[[comp_vs]] <- sub_df
  }
  
  result
}
# Infer groups from sample names: strip trailing _N (replicate number)
# Examples:
#   "B1_1" -> "B1"
#   "BD_2" -> "BD"
#   "20251118_B1_1" -> "B1" (if date prefix wasn't already removed)
infer_groups <- function(sample_names) {
  cleaned <- sub("^[0-9]{8}_", "", sample_names)
  groups <- sub("_[0-9]+$", "", cleaned)
  groups <- make_valid_group_name(groups)
  
  groups
}
# All unique pairwise comparisons as "A_vs_B" strings
all_pairs <- function(groups) {
  g <- sort(unique(groups))
  if (length(g) < 2) return(character(0))
  pairs <- combn(g, 2, simplify = FALSE)
  vapply(pairs, function(p) paste(p[1], "vs", p[2], sep = "_"), character(1))
}

# Shared local Reactome cache. The cache stores pathway membership only; formal
# enrichment still uses the current dataset as the statistical background.
reactome_cache_path <- function() {
  file.path("data", "reactome_cache.rds")
}

empty_reactome_cache <- function() {
  data.frame(
    gene_symbol = character(0),
    entrezid = character(0),
    pathid = character(0),
    reactomeid = character(0),
    pathname = character(0),
    species = character(0),
    last_updated = character(0),
    stringsAsFactors = FALSE
  )
}

read_reactome_cache <- function(path = reactome_cache_path()) {
  if (!file.exists(path)) return(empty_reactome_cache())
  cache <- tryCatch(readRDS(path), error = function(e) empty_reactome_cache())
  needed <- names(empty_reactome_cache())
  for (nm in setdiff(needed, names(cache))) cache[[nm]] <- NA_character_
  cache[, needed, drop = FALSE]
}

write_reactome_cache <- function(cache, path = reactome_cache_path()) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  cache <- unique(cache[, names(empty_reactome_cache()), drop = FALSE])
  saveRDS(cache, path)
  invisible(cache)
}

split_gene_symbols <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x)]
  x <- unlist(strsplit(x, "[;|,]"))
  x <- trimws(x)
  x <- x[nzchar(x) & x != "NA"]
  unique(x)
}

dataset_gene_symbols <- function(parsed_data_or_meta) {
  meta <- if (is.list(parsed_data_or_meta) && "meta" %in% names(parsed_data_or_meta)) {
    parsed_data_or_meta$meta
  } else {
    parsed_data_or_meta
  }
  if (is.null(meta) || !"Genes" %in% names(meta)) return(character(0))
  split_gene_symbols(meta$Genes)
}

reactome_cache_status <- function(cache = read_reactome_cache()) {
  mapped <- cache[!is.na(cache$pathid) & nzchar(cache$pathid), , drop = FALSE]
  list(
    rows = nrow(cache),
    genes = length(unique(cache$gene_symbol)),
    mapped_genes = length(unique(mapped$gene_symbol)),
    pathways = length(unique(mapped$pathid))
  )
}

build_reactome_cache_for_symbols <- function(symbols, cache_path = reactome_cache_path()) {
  symbols <- unique(split_gene_symbols(symbols))
  cache <- read_reactome_cache(cache_path)
  cached <- unique(cache$gene_symbol)
  missing <- setdiff(symbols, cached)

  if (!length(missing)) {
    return(list(cache = cache, added = 0L, queried = 0L))
  }

  if (!requireNamespace("AnnotationDbi", quietly = TRUE) ||
      !requireNamespace("org.Hs.eg.db", quietly = TRUE) ||
      !requireNamespace("reactome.db", quietly = TRUE)) {
    stop("Reactome cache requires AnnotationDbi, org.Hs.eg.db, and reactome.db.")
  }

  mapped <- suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = missing,
    columns = "ENTREZID",
    keytype = "SYMBOL"
  ))
  mapped <- mapped[!is.na(mapped$ENTREZID) & nzchar(mapped$ENTREZID), , drop = FALSE]

  react <- empty_reactome_cache()
  if (nrow(mapped)) {
    react_raw <- suppressMessages(AnnotationDbi::select(
      reactome.db::reactome.db,
      keys = unique(as.character(mapped$ENTREZID)),
      columns = c("PATHID", "PATHNAME"),
      keytype = "ENTREZID"
    ))
    react_raw <- react_raw[!is.na(react_raw$PATHID) & nzchar(react_raw$PATHID), , drop = FALSE]
    if (nrow(react_raw)) {
      react <- merge(mapped, react_raw, by = "ENTREZID", all.x = FALSE, all.y = FALSE)
      react <- data.frame(
        gene_symbol = react$SYMBOL,
        entrezid = as.character(react$ENTREZID),
        pathid = as.character(react$PATHID),
        reactomeid = as.character(react$PATHID),
        pathname = sub("^Homo sapiens:\\s*", "", as.character(react$PATHNAME)),
        species = "Homo sapiens",
        last_updated = as.character(Sys.Date()),
        stringsAsFactors = FALSE
      )
    }
  }

  nohit <- setdiff(missing, unique(react$gene_symbol))
  if (length(nohit)) {
    react <- rbind(
      react,
      data.frame(
        gene_symbol = nohit,
        entrezid = NA_character_,
        pathid = NA_character_,
        reactomeid = NA_character_,
        pathname = NA_character_,
        species = NA_character_,
        last_updated = as.character(Sys.Date()),
        stringsAsFactors = FALSE
      )
    )
  }

  cache <- unique(rbind(cache, react))
  write_reactome_cache(cache, cache_path)
  list(cache = cache, added = nrow(react), queried = length(missing))
}

reactome_membership_hints <- function(symbols, max_terms = 3L, cache = read_reactome_cache()) {
  symbols <- split_gene_symbols(symbols)
  if (!length(symbols) || !nrow(cache)) return("")
  sub <- cache[cache$gene_symbol %in% symbols & !is.na(cache$pathname) & nzchar(cache$pathname), , drop = FALSE]
  if (!nrow(sub)) return("")
  counts <- sort(table(sub$pathname), decreasing = TRUE)
  paste(names(counts)[seq_len(min(max_terms, length(counts)))], collapse = "; ")
}

reactome_term_tables_from_cache <- function(symbols, cache = read_reactome_cache()) {
  symbols <- unique(split_gene_symbols(symbols))
  sub <- cache[cache$gene_symbol %in% symbols & !is.na(cache$pathid) & nzchar(cache$pathid), , drop = FALSE]
  if (!nrow(sub)) {
    return(list(t2g = data.frame(), t2n = data.frame()))
  }
  list(
    t2g = unique(data.frame(term = sub$pathid, gene = sub$gene_symbol, stringsAsFactors = FALSE)),
    t2n = unique(data.frame(term = sub$pathid, name = sub$pathname, stringsAsFactors = FALSE))
  )
}

reactome_entrez_term_tables_from_cache <- function(cache = read_reactome_cache()) {
  sub <- cache[!is.na(cache$entrezid) & nzchar(cache$entrezid) &
                 !is.na(cache$pathid) & nzchar(cache$pathid), , drop = FALSE]
  if (!nrow(sub)) {
    return(list(t2g = data.frame(), t2n = data.frame()))
  }
  list(
    t2g = unique(data.frame(term = sub$pathid, gene = sub$entrezid, stringsAsFactors = FALSE)),
    t2n = unique(data.frame(term = sub$pathid, name = sub$pathname, stringsAsFactors = FALSE))
  )
}
