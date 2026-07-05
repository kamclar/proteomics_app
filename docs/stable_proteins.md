# Stable Proteins tab

## Purpose

The `Stable Proteins` tab separates selected DEqMS results into three practical
groups:

- `Changed`: proteins that meet the selected adjusted p-value and log2FC
  thresholds in at least one selected comparison.
- `Stable candidate`: proteins that are non-significant in all selected
  comparisons, have small maximum absolute log2FC, sufficient measured values,
  and no on/off pattern.
- `Other tested`: proteins that were tested but are not confidently changed or
  stable. This includes low coverage proteins, on/off proteins, proteins with
  missing DEqMS statistics, and borderline proteins.

`Stable candidate` should not be read as proven unchanged. It is a conservative
QC and interpretation category.

## Added views

- Summary counts for stable, changed, and other tested proteins.
- A class-size bar plot.
- A classification map using maximum absolute log2FC and minimum adjusted
  p-value across the selected comparisons.
- An `Other tested` table with a short reason for each protein.
- Reactome pathway comparison for stable and changed proteins using the local
  Reactome cache.

## Interpretation

Changed proteins are the primary biological signal. Stable candidates provide a
background-like context and can be useful for QC, but they are not formal proof
of no biological change. Other tested proteins are intentionally kept visible
because they often explain ambiguity: low data coverage, missing statistics,
on/off behavior, or fold changes that do not pass the selected statistical
thresholds.

Reactome results for changed proteins can suggest affected pathways. Reactome
results for stable candidates should be interpreted as pathway context among
well-covered proteins, not as evidence that a pathway is actively stable.

## Fallback policy

The app avoids biological fallbacks in this tab. If gene symbols are unavailable,
Reactome mapping is skipped rather than guessed from protein IDs or row names. If
DEqMS statistics or raw coverage are unavailable, proteins are not classified as
stable candidates. They remain visible as `Other tested` with a reason.

## Sources

- DEqMS: differential protein expression analysis for proteomics with peptide
  or PSM count-aware variance modeling.
  https://bioconductor.org/packages/DEqMS/
- limma `decideTests`: notes on multiple testing, adjusted p-values, and
  fold-change thresholds.
  https://rdrr.io/bioc/limma/man/decideTests.html
- limma `eBayes` / `treat`: moderated statistics and testing relative to a
  scientifically meaningful fold-change threshold.
  https://rdrr.io/bioc/limma/man/ebayes.html
- ProtRank paper: missing values can strongly affect differential proteomics and
  should be handled carefully.
  https://arxiv.org/abs/1909.13667
- Missing values in proteomics review/perspective.
  https://arxiv.org/abs/2304.06654
- clusterProfiler `enricher`: enrichment should use an appropriate `universe`
  background when available.
  https://rdrr.io/bioc/clusterProfiler/man/enricher.html
