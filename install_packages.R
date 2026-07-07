# install_packages.R
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
app_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}
setwd(app_dir)

local_lib <- file.path(app_dir, "packages")
dir.create(local_lib, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(local_lib, .libPaths()))

cran_repo <- Sys.getenv("PROTEOMICS_APP_CRAN_REPO", "https://cloud.r-project.org")
options(repos = c(CRAN = cran_repo))

installed_local <- function() {
  rownames(installed.packages(lib.loc = local_lib))
}

loadable <- function(packages) {
  vapply(packages, requireNamespace, logical(1), quietly = TRUE)
}

install_cran_missing <- function(packages) {
  missing <- setdiff(packages, installed_local())
  if (!length(missing)) return(invisible(NULL))
  cat("Instaluji CRAN balicky do packages:\n")
  cat(paste(missing, collapse = ", "), "\n")
  install.packages(missing, lib = local_lib, repos = cran_repo)
}

install_bioc_missing <- function(packages) {
  missing <- setdiff(packages, installed_local())
  if (!length(missing)) return(invisible(NULL))
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", lib = local_lib, repos = cran_repo)
  }
  cat("Instaluji Bioconductor balicky do packages:\n")
  cat(paste(missing, collapse = ", "), "\n")
  BiocManager::install(missing, lib = local_lib, ask = FALSE, update = FALSE)
}

reinstall_cran <- function(packages) {
  if (!length(packages)) return(invisible(NULL))
  cat("Preinstalovavam CRAN balicky, ktere R neumi nacist:\n")
  cat(paste(packages, collapse = ", "), "\n")
  install.packages(packages, lib = local_lib, repos = cran_repo)
}

reinstall_bioc <- function(packages) {
  if (!length(packages)) return(invisible(NULL))
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", lib = local_lib, repos = cran_repo)
  }
  cat("Preinstalovavam Bioconductor balicky, ktere R neumi nacist:\n")
  cat(paste(packages, collapse = ", "), "\n")
  BiocManager::install(packages, lib = local_lib, ask = FALSE, update = FALSE, force = TRUE)
}

cran_packages <- c(
  "shiny",
  "shinydashboard",
  "DT",
  "plotly",
  "openxlsx",
  "ggplot2",
  "dplyr",
  "zip",
  "RColorBrewer",
  "imputeLCMD"
)

bioc_packages <- c(
  "limma",
  "DEqMS",
  "impute"
)

install_cran_missing(cran_packages)
install_bioc_missing(bioc_packages)

load_status <- loadable(c(cran_packages, bioc_packages))
not_loadable <- names(load_status)[!load_status]
if (length(not_loadable)) {
  reinstall_cran(intersect(not_loadable, cran_packages))
  reinstall_bioc(intersect(not_loadable, bioc_packages))
}

required <- c(cran_packages, bioc_packages)
missing_after <- setdiff(required, installed_local())
load_after <- loadable(required)
not_loadable_after <- names(load_after)[!load_after]

cat("\n=== HOTOVO ===\n")
cat("Knihovna:", local_lib, "\n")
if (length(missing_after)) {
  cat("Stale chybi:\n")
  print(missing_after)
  quit(status = 1)
}

if (length(not_loadable_after)) {
  cat("Tyto balicky jsou v knihovne, ale R je neumi nacist:\n")
  print(not_loadable_after)
  cat("R library paths:\n")
  print(.libPaths())
  quit(status = 1)
}

cat("Zakladni balicky aplikace jsou dostupne a nacitatelne z packages/.\n")
