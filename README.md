# Congruent life stage hotspots identify essential fish habitats to fill biodiversity conservation gaps 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18632166.svg)](https://doi.org/10.5281/zenodo.18632166) [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [![Data: CC BY 4.0](https://img.shields.io/badge/Data%20License-CC%20BY%204.0-brightgreen.svg)](https://creativecommons.org/licenses/by/4.0/)

Open-science code and data to model, map, and analyze biodiversity-based Essential Fish Habitats (EFHs) across the Western Mediterranean Sea using hierarchical joint species distribution models (HMSC) and hurdle models, followed by spatio-temporal prediction, hotspot analysis and MPA coverage assessment.

All figures can be regenerated directly from the provided outputs (once the hurdle stack is downloaded from Zenodo), without rerunning scripts 1–5.

## Purpose ![target](https://img.shields.io/badge/-purpose-37474f?style=flat-square)
- Quantify species distributions and abundances with HMSC (presence/absence and abundance components).
- Generate spatio-temporal predictions and hurdle-model products for multiple life stages.
- Identify emerging hotspots/coldspots and evaluate marine protected area (MPA) coverage of EFH patterns.

## Methods and Workflow ![workflow](https://img.shields.io/badge/-workflow-00695c?style=flat-square)
1. **Model specification**: Assemble inputs and define HMSC models ([1-Define HMSC models.r](1-Define%20HMSC%20models.r)).
2. **Model fitting**: Fit presence/absence and abundance models (supports HPC workflows) ([2_1-Fit presence absence model.r](2_1-Fit%20presence%20absence%20model.r), [2_2-Fit abundance model.r](2_2-Fit%20abundance%20model.r)); Slurm job scripts for cluster runs are available in [Models/](Models).
3. **Diagnostics**: Assess convergence, effective sample sizes, and PSRF ([3-Models convergence and fit.r](3-Models%20convergence%20and%20fit.r)).
4. **Spatial predictions**: Generate posterior samples and means over the spatial grid ([4-Spatial predictions with fitted HMSC.r](4-Spatial%20predictions%20with%20fitted%20HMSC.r)).
5. **Extrapolation assessment**: Compute Multivariate Environmental Similarity Surfaces (MESS) to identify where spatial predictions extrapolate beyond the environmental/anthropogenic conditions used to fit both HMSC models ([Multivariate Environmental Similarity Surfaces.R](Multivariate%20Environmental%20Similarity%20Surfaces.R)).
6. **Hurdle outputs**: Combine presence/absence and abundance to produce hurdle predictions and derived rasters ([5-Hurdle model predictions.r](5-Hurdle%20model%20predictions.r)).
7. **Hotspot analysis**: Emerging Hotspot Analysis (EHSA) and visualization scripts ([6-Emerging Hotspot Analysis.r](6-Emerging%20Hotspot%20Analysis.r), [Figure_3_4.r](Figure_3_4.r), [Figure_5_6.r](Figure_5_6.r)).
8. **MPA coverage**: Quantify overlap of hotspot classes with marine protected areas ([Figure_7.r](Figure_7.r)).
9. **Trend/hotspot cross-tabulation**: Cross-tabulate EHSA classes (hotspot/coldspot/no pattern) against Sen's slope trend significance of hurdle-model abundance, by life stage ([Figure_S_EHSA_Sen_slope_crosstab.r](Figure_S_EHSA_Sen_slope_crosstab.r)).

## Inputs and Data ![data](https://img.shields.io/badge/-data-283593?style=flat-square)
- **Primary inputs**: RData objects and spatial layers in [Inputs_HMSC/](Inputs_HMSC).
  - Species/community matrices, traits, and covariates (e.g., `allData.RData`, `Traits.RData`).
  - Prediction grids and protection layers (shapefiles in `Grid/` and `protectionPLactivities_niv1&2&3&4/`).
- **Model checkpoints**: Fitted models and MCMC outputs in [Models/](Models).
- **Hurdle prediction stack**: Spatio-temporal raster stack available on Zenodo [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18632166.svg)](https://doi.org/10.5281/zenodo.18632166). Use [download_zenodo_spatiotemporal.R](download_zenodo_spatiotemporal.R) to place `r_stack_Hurdle_0.05_0_1000_19992021.Rdata` into [Outputs/Spatio_temporal_prediction/](Outputs/Spatio_temporal_prediction); download this file before running scripts that depend on the hurdle stack (e.g., [6-Emerging Hotspot Analysis.r](6-Emerging%20Hotspot%20Analysis.r), [Figure_2.r](Figure_2.r), [Figure_3_4.r](Figure_3_4.r), [Figure_5_6.r](Figure_5_6.r)).

## Outputs ![outputs](https://img.shields.io/badge/-outputs-5d4037?style=flat-square)
- Posterior samples and means for presence/absence and abundance: [Outputs/Presence_Absence_posteriors/](Outputs/Presence_Absence_posteriors), [Outputs/Abundance_posteriors/](Outputs/Abundance_posteriors), [Outputs/Presence_Absence_mean_posterior/](Outputs/Presence_Absence_mean_posterior), [Outputs/Abundance_mean_posterior/](Outputs/Abundance_mean_posterior). Full posterior outputs are several hundred GB and not hosted on Zenodo; they can be requested from [AUTHOR_NAME].
- Hurdle predictions (combined): [Outputs/Hurdle_prediction/](Outputs/Hurdle_prediction). These combined hurdle outputs are large and not mirrored on Zenodo; contact [AUTHOR_NAME] if you need access.
- Spatio-temporal hurdle stack: [Outputs/Spatio_temporal_prediction/](Outputs/Spatio_temporal_prediction). Source data are hosted on Zenodo (see the hurdle prediction stack entry above).
- EHSA results: [Outputs/EHSA/](Outputs/EHSA).
- Figures and supplementary plots: [Figures/](Figures).

## Installation and Environment ![setup](https://img.shields.io/badge/-setup-546e7a?style=flat-square)
- **R version**: Tested with R (>= 4.2.x).
- **Core packages**: Listed in [0-Load libraries.r](0-Load%20libraries.r); major dependencies include `Hmsc`, `raster`, `sf`, `ggplot2`, `dplyr`, `zenodor`, and others.
- **Install dependencies**: The loader script installs missing CRAN packages; `zenodor` is installed from GitHub. From the repo root:

```r
source("0-Load libraries.r")
```

If `zenodor` is absent, it is pulled via `remotes::install_github("FRBCesab/zenodor")` inside [download_zenodo_spatiotemporal.R](download_zenodo_spatiotemporal.R).

## Basic Usage ![run](https://img.shields.io/badge/-run-1b5e20?style=flat-square)
1. Load packages and data paths:
   ```r
   source("0-Load libraries.r")
   ```
2. Define and prepare models:
   ```r
   source("1-Define HMSC models.r")
   ```
3. Fit models (presence/absence and abundance) — adjust settings for local vs HPC:
   ```r
   source("2_1-Fit presence absence model.r")
   source("2_2-Fit abundance model.r")
   ```
4. Run diagnostics and convergence checks:
   ```r
   source("3-Models convergence and fit.r")
   ```
5. Produce spatial predictions and hurdle outputs:
   ```r
   source("4-Spatial predictions with fitted HMSC.r")
   source("5-Hurdle model predictions.r")
   ```
6. Retrieve the published hurdle stack (optional if not already present):
   ```r
   source("download_zenodo_spatiotemporal.R")
   ```
7. Run Emerging Hotspot Analysis (EHSA):
   ```r
   source("6-Emerging Hotspot Analysis.r")
   ```
8. Generate figures as needed (e.g., `Figure_1.r`, `Figure_2.r`, `Figure_3_4.r`, `Figure_5_6.r`, `Figure_7.r` or `Figure_S*.r`).

## Repository Structure ![folders](https://img.shields.io/badge/-folders-3949ab?style=flat-square)
- [0-Load libraries.r](0-Load%20libraries.r): Dependency management and loading.
- [1-Define HMSC models.r](1-Define%20HMSC%20models.r): Model setup and data preparation.
- [2_1-Fit presence absence model.r](2_1-Fit%20presence%20absence%20model.r), [2_2-Fit abundance model.r](2_2-Fit%20abundance%20model.r): Model fitting workflows.
- [3-Models convergence and fit.r](3-Models%20convergence%20and%20fit.r): Diagnostics.
- [4-Spatial predictions with fitted HMSC.r](4-Spatial%20predictions%20with%20fitted%20HMSC.r): Spatial posterior predictions.
- [Multivariate Environmental Similarity Surfaces.R](Multivariate%20Environmental%20Similarity%20Surfaces.R): MESS extrapolation diagnostics.
- [5-Hurdle model predictions.r](5-Hurdle%20model%20predictions.r): Combine components into hurdle outputs.
- [6-Emerging Hotspot Analysis.r](6-Emerging%20Hotspot%20Analysis.r): Emerging Hotspot Analysis workflow.
- Figure-generation scripts: all [Figure_*.r](Figure_1.r) scripts in the repository root (including main and supplementary figures; e.g., [Figure_7.r](Figure_7.r), [Figure_S8_S9.r](Figure_S8_S9.r), [Figure_S_EHSA_Sen_slope_crosstab.r](Figure_S_EHSA_Sen_slope_crosstab.r)).
- [:file\_folder: **Inputs_HMSC**](Inputs_HMSC/): Input data and spatial layers.
- [:file\_folder: **Models**](Models/): Fitted model objects and initialization files.
- [:file\_folder: **Outputs**](Outputs/): Generated predictions, EHSA results, and hurdle outputs (git-ignored intermediates where appropriate).

## Citation ![cite](https://img.shields.io/badge/-cite-4e342e?style=flat-square)
If you use this code or data, please cite:
- **Dataset**: [AUTHOR_NAME] (2026). *Community-based Essential Fish Habitats in Mediterranean Sea* (Version V1.0.0). Zenodo. [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18632166.svg)](https://doi.org/10.5281/zenodo.18632166)
- **Code**: This GitHub repository (XXXXX/Biodiversity_based_EFHs, main branch).
- **Associated publication**: XXXXXXXX.

## License ![license](https://img.shields.io/badge/-license-263238?style=flat-square)
This project is released under the GNU General Public License v3.0. See [LICENSE](LICENSE) for details.
