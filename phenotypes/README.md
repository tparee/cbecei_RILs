# Phenotypes Directory

This directory contains phenotypic data for:

- Growth rate
- Sex-specific body size

---

# Growth Rate

## `growthrates.csv`

Contains raw growth rate measurements.

### Columns

- `strain`  
  RIL identifier (`QGXXXX`)

- `hours_to_starve`  
  Number of hours until resource exhaustion.

  Measurements were obtained by placing:

  - two adult females
  - two adult males

  on a plate seeded with *E. coli* and monitoring the population until it peaked and subsequently crashed due to starvation.

  Plates were imaged hourly to estimate the starvation time  
  (see https://doi.org/10.1093/genetics/iyaf073).

- `block`  
  Independent experimental block.

---

## `logGrowthrates_emmeans.csv`

Contains BLUEs (Best Linear Unbiased Estimates) estimated using the `emmeans` R package from the model:

```r
lm(log_hours_to_starve ~ block + strain, data = growthrates)
```

These estimates were used only for plotting Figure 3A.

Raw data were used for:

- GWAS analyses
- Heritability estimates

---

# Sex-Specific Size

## `size_data_raw.csv.gz`

Contains raw worm length measurements.

### Experimental Design

- Worms were imaged 48 hours after the L1 stage.
- Imaging was performed in 96-well plates using an ImageXpress instrument.
- Worm objects were detected and measured using the Andersen lab pipeline  
  (see https://doi.org/10.1371/journal.pone.0252000).

Objects were then classified by sex using an XGBoost model trained to distinguish males and females (see `suppl/ImageXpress_sex_xgboost`).

---

## `size_summarystat.csv`

Contains summary statistics for each assayed plate (i.e., one strain within one experimental block).

### Included Statistics

- Mean worm length (`µm`) for:
  - females
  - males

- Unscaled values:
  - `mean_worm_length_um_f`
  - `mean_worm_length_um_m`

- Scaled values (variance = 1 and mean = 0 within each sex):
  - `scaled_mean_worm_length_um_f`
  - `scaled_mean_worm_length_um_m`

- Standard errors:
  - `(scaled_)SE_f`
  - `(scaled_)SE_m`

### Sexually Divergent and Convergent Axes

The table also includes:

- `div`
- `conv`

These correspond to the coordinates of each strain after a 45° rotation of the male-versus-female size plane.

See:

- Methods section
- https://doi.org/10.1371/journal.pbio.3000244

for an example and methodological details.

---

## `size_emmeans.csv`

Contains BLUEs estimated using the `emmeans` R package from the model:

```r
lm(worm_length_um ~ sex + block + strain:sex + sex:block, data = size)
```

These estimates were used only for plotting Figure 3B.

Summary statistics were used for:

- GWAS analyses
- Heritability estimates
