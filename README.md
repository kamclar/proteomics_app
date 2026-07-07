# Proteomics App

Shiny application for exploratory proteomics data analysis. The app supports data upload, missing value inspection, imputation, DEqMS and t-test based comparison, volcano plots, enrichment analysis, and a stable protein overview.

## Repository content

- `app.R`: main Shiny application.
- `modules/`: app modules for upload, imputation, exploration, DEqMS, volcano plots, enrichment, stable proteins, and comparison.
- `docs/`: short user documentation for selected app features.
  See `docs/imputation_methods.md` for the current imputation methods.
- `data/reactome_cache.rds`: small Reactome mapping cache used to speed up Reactome annotations.
- `www/`: static app assets.
- `run_app.bat`: Windows launcher for the app.
- `install_r_runtime.bat`: Windows helper script for installing R if no R runtime is available.
- `install_packages.bat` and `install_packages.R`: helper scripts for installing local R packages into the `packages/` folder.

The repository does not include local R package libraries, local R runtime files.

## Running the app

On Windows, run:

```bat
run_app.bat
```

The launcher first looks for a bundled R runtime in `R-runtime/`. If it is not found, it tries to use an installed `Rscript.exe` from the system.

You can also run the app from R:

```r
shiny::runApp(".", launch.browser = TRUE, port = 3838)
```

## Installing R

If the computer does not have R, run:

```bat
install_r_runtime.bat
```

This script checks for a bundled `R-runtime/` first. If no runtime or system R
is found, it installs R for Windows using `winget` package `RProject.R`.

After R is installed, run:

```bat
install_packages.bat
run_app.bat
```

## Required packages

The app expects the main R packages to be available either in the local `packages/` folder or in the active R library.

Core packages include:

- `shiny`
- `shinydashboard`
- `DT`
- `plotly`
- `openxlsx`
- `limma`
- `DEqMS`
- `ggplot2`
- `dplyr`
- `zip`
- `impute`
- `imputeLCMD`
- `RColorBrewer`

## Installing local packages

To install local packages into the `packages/` folder, run:

```bat
install_packages.bat
```

or:

```bat
R-runtime\bin\Rscript.exe install_packages.R
```

The current script installs the core app packages locally, including Shiny, Plotly, DEqMS, limma, `impute`, and `imputeLCMD`. Optional enrichment packages are not installed by default.

## Optional enrichment packages

GO, KEGG, and Reactome enrichment require additional Bioconductor packages. These are intentionally not bundled in the repository because they can be large.

Install them only when enrichment features are needed. The app can also offer installation from the Enrichment tab, depending on the local R setup.

Common optional packages:

- `clusterProfiler`
- `org.Hs.eg.db`
- `ReactomePA`
- `enrichplot`
- `AnnotationDbi`
- `KEGGREST`
- `reactome.db`

## Portable use

For a portable folder, keep these folders next to the app:

- `R-runtime/` with a compatible Windows R runtime
- `packages/` with the required local packages

These folders are ignored by Git and should be distributed separately when needed.
