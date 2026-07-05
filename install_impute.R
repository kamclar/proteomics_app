# install_impute.R
local_lib <- normalizePath("packages", winslash = "/", mustWork = FALSE)
dir.create(local_lib, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(local_lib, .libPaths()))

cat("Instaluji BiocManager...\n")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = local_lib, repos = "https://cloud.r-project.org")
}

cat("Instaluji impute...\n")
BiocManager::install("impute", lib = local_lib, ask = FALSE, update = FALSE)

cat("Instaluji imputeLCMD...\n")
install.packages("imputeLCMD", lib = local_lib, repos = "https://cloud.r-project.org")

cat("\n=== HOTOVO ===\n")
cat("Balicky v packages:\n")
print(list.files(local_lib))
