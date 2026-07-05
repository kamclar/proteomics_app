# Volcano data flags

This note documents the Volcano tab data flags, why they exist, and what they
should not be interpreted as.


## Low-intensity partial-missing flag

In proteomics, low abundance is not usually a fixed binary category. It depends on instrument sensitivity,
sample complexity, preprocessing, and the current comparison.

The Volcano tab therefore uses the more cautious label:

`Low-intensity partial missing`

Current operational definition:

1. Use only samples from the selected comparison.
2. Compute the raw mean intensity across those samples.
3. Keep only proteins with partial missingness in the selected comparison:
   at least one measured value and at least one missing value.
4. Exclude on/off proteins, because those are a separate missingness pattern.
5. Flag proteins in the lowest 10% of comparison-specific raw mean intensity
   among this partial-missing set.

Interpretation:

- This is a QC and missingness-context flag.
- It is not a biological statement that the protein is globally low abundance.
- It should not override DEqMS or t-test evidence.
- A strong DEqMS hit can be low intensity, but should only be flagged here if
  it is also partially missing in the selected comparison.

## On/off flag

The on/off flag marks proteins that are measured in one group of the selected
comparison and absent in the other group.

Interpretation:

- This can be biologically meaningful, especially for proteins expected to be
  condition-specific.
- It can also reflect detection limits or stochastic missingness.
- If a finite p-value cannot be computed, the Volcano tab plots the point near
  the baseline and reports the p-value as not available in the hover text. This
  avoids visually presenting display-only points as extremely significant.

## Proteomics background

Missing values in mass-spectrometry proteomics are often linked to detection
limits and low signal. In missing-data terminology, this is often discussed as
left-censoring or MNAR: the probability of a value being missing can depend on
the unobserved low abundance value. However, missingness can also come from
technical acquisition, matching, or processing effects.

For this reason, the app keeps these as visual QC flags rather than statistical
calls.

Useful references:

- Missing data mechanisms, including MNAR:
  https://en.wikipedia.org/wiki/Missing_data
- ProtRank discussion of missing values and infinite fold-change issues in
  proteomics:
  https://arxiv.org/abs/1909.13667
- Vanderaa and Gatto, "Revisiting the thorny issue of missing values in
  single-cell proteomics":
  https://arxiv.org/abs/2304.06654
- General proteomics context on detection range and low-abundance peptides:
  https://en.wikipedia.org/wiki/Quantitative_proteomics

## Future decisions

- Whether the 10% cutoff should become user-adjustable.
- Whether the flag should use raw mean intensity, median intensity, or a
  detection-limit model.
- Whether the flagged table should include group-level valid counts by default.
