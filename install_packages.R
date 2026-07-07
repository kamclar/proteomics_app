# install_packages.R
local_lib <- normalizePath("packages", winslash = "/", mustWork = FALSE)
dir.create(local_lib, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(local_lib, .libPaths()))

cran_repo <- Sys.getenv("PROTEOMICS_APP_CRAN_REPO", "https://cloud.r-project.org")
options(repos = c(CRAN = cran_repo))

installed_local <- function() {
  rownames(installed.packages(lib.loc = local_lib))
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

required <- c(cran_packages, bioc_packages)
missing_after <- setdiff(required, installed_local())

cat("\n=== HOTOVO ===\n")
cat("Knihovna:", local_lib, "\n")
if (length(missing_after)) {
  cat("Stale chybi:\n")
  print(missing_after)
  quit(status = 1)
}

cat("Zakladni balicky aplikace jsou dostupne v packages/.\n")
