# Proteomics App

Shiny application for exploratory proteomics data analysis. The app supports data upload, missing value inspection, imputation, DEqMS and t-test based comparison, volcano plots, enrichment analysis, and a stable protein overview.

## Repository content

- `app.R`: main Shiny application.
- `modules/`: app modules for upload, imputation, exploration, DEqMS, volcano plots, enrichment, stable proteins, and comparison.
- `docs/`: short user documentation for selected app features.
- `data/reactome_cache.rds`: small Reactome mapping cache used to speed up Reactome annotations.
- `www/`: static app assets.
- `run_app.bat`: Windows launcher for the app.
- `push_readme.bat`: helper script for committing and pushing README-only changes.
- `push_changes.bat`: helper script for committing and pushing all non-ignored local changes.
- `install_impute.bat` and `install_impute.R`: helper scripts for installing imputation packages into the local `packages/` folder.
- `packages.R`: optional package preparation entry point.

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

## Pushing README changes

To commit and push only `README.md` changes to GitHub, run:

```bat
push_readme.bat
```

The script asks for a commit message and pushes to the `main` branch.

## Pushing app changes

To commit and push all non-ignored local changes to GitHub, run:

```bat
push_changes.bat
```

The script shows the changed files first, asks for confirmation, asks for a commit message, and pushes to the `main` branch.

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

## Installing imputation packages

To install imputation packages into the local `packages/` folder, run:

```bat
install_impute.bat
```

or:

```bat
R-runtime\bin\Rscript.exe install_impute.R
```

This installs `impute` and `imputeLCMD` locally, so the imputation methods can run without using packages from a user library.

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
