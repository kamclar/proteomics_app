# Imputation methods

This note describes the imputation options currently available in the app.

## When to impute

Imputation replaces missing intensity values. It should be treated as a
modeling decision. Missing values can reflect low abundance and detection
limits, but they can also come from acquisition, matching, preprocessing, or
other technical effects.

The app keeps the random seed with the imputation result. Reusing the same seed,
data, and settings gives the same random MNAR values.

## Row pre-filter

Before imputation, the app can remove protein rows with fewer than the selected
number of valid values across all samples. The default is 1, which removes only
rows with no measured values. Removed rows are excluded from downstream steps,
and the app reports how many rows were removed.

## MNAR left-censored

Use this when missing values are expected to represent low-abundance proteins
below the detection limit. The app offers three MNAR backends.

### MNAR (Perseus default: width = 0.3, downshift = 1.8)

This is the default MNAR method. For each sample column, missing values are
drawn from a low-intensity normal distribution:

- mean = observed column mean minus `downshift * observed column SD`
- SD = `width * observed column SD`

Default settings are `width = 0.3` and `downshift = 1.8`, matching the label
shown in the app. Advanced settings allow both parameters to be adjusted.

Use it for simple method comparison and for low-intensity missingness where a
Perseus-compatible workflow is desired.

### QRILC

QRILC is provided by `imputeLCMD::impute.QRILC`. It is designed for
left-censored missing values and uses quantile regression to estimate the left
tail of the intensity distribution.

The app exposes `tune.sigma`, which controls the spread of the imputed MNAR
distribution.

Use it when low-intensity censoring is expected and a package-backed
proteomics-specific method is preferred.

### MinProb

MinProb is provided by `imputeLCMD::impute.MinProb`. It imputes missing values
by random draws from a Gaussian distribution centered near a low sample-specific
quantile.

The app exposes:

- `q`, the low quantile used to estimate the low-intensity center
- `tune.sigma`, which controls the spread of the imputed MNAR distribution

Use it as a simpler MNAR alternative to QRILC.

## MNAR left-censored per group

This uses the selected MNAR backend separately within each condition group. It
is useful when missingness is condition-specific and should not be estimated
from all samples pooled together.

The available backends are the same as for MNAR left-censored:

- MNAR (Perseus default: width = 0.3, downshift = 1.8)
- QRILC
- MinProb

## MAR kNN

kNN imputation is provided by `impute::impute.knn`. It estimates missing values
from proteins with similar abundance profiles.

The app exposes `k`, the number of neighbours. The default is 5.

The app adds a guard around sparse rows. If a protein has 50% or more missing
values among the columns being imputed, it is not sent to kNN. Instead, it is
filled with the selected MNAR backend and the app shows a warning. This avoids
using the package fallback where highly missing rows can be filled by column
means.

If `Compute per group` is enabled, kNN is run separately inside each condition
group. Otherwise, kNN is run once on the pooled matrix.

Use kNN when missingness is expected to be closer to MAR, where similar proteins
carry useful information about the missing values.

## Mixed MNAR/MAR

Mixed mode is an app-level decision layer. The app classifies missing values
first, then uses package-backed or Perseus-style methods for the actual filling.

Without per-group mode:

- proteins with missing fraction at or above the selected threshold use MNAR
- proteins with partial missingness below the threshold use kNN
- complete proteins are left unchanged

With per-group mode:

- missingness is classified separately within each condition group
- `Smart / hybrid` classifies by group, then runs kNN once on the full matrix
  for MAR values
- `Fully per group` classifies and imputes separately inside each group

The threshold sensitivity panel shows how many protein-by-group cases switch
between MNAR and kNN as the mixed threshold changes.

Use mixed mode when the dataset likely contains both low-intensity censoring and
more random missingness.

## Skip imputation

Use `Skip - data is already imputed` when the uploaded data are already complete
or were imputed outside the app. The app then passes the uploaded intensity
matrix forward unchanged.

## Sources

- `imputeLCMD` package manual. Left-censoring is described as an MNAR mechanism
  relevant to proteomics, and the manual documents QRILC, MinProb, and mixed
  MAR/MNAR approaches:
  https://cran.r-project.org/web/packages/imputeLCMD/imputeLCMD.pdf
- `imputeLCMD::impute.QRILC` documentation:
  https://rdrr.io/cran/imputeLCMD/man/impute.QRILC.html
- `imputeLCMD::impute.MinProb` documentation:
  https://rdrr.io/cran/imputeLCMD/man/impute.MinProb.html
- Bioconductor `impute` package:
  https://bioconductor.org/packages/release/bioc/html/impute.html
- `impute::impute.knn` documentation, including nearest-neighbour averaging,
  sparse-row behavior, and reproducibility seed:
  https://rdrr.io/bioc/impute/man/impute.knn.html
- Troyanskaya et al. 2001, missing value estimation by k-nearest neighbours for
  expression matrices:
  https://doi.org/10.1093/bioinformatics/17.6.520
- Missing values in proteomics perspective:
  https://arxiv.org/abs/2304.06654
