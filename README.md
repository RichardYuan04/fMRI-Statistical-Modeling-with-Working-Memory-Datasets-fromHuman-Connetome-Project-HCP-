# fMRI Statistical Modeling with Working Memory Datasets from Human Connectome Project (HCP)

This project investigates how cortical activation across 360 brain regions responds to working memory (WM) task performance under different visual stimuli (faces and places/scenes), using high-quality fMRI data from the Human Connectome Project (HCP).

## 🧠 Project Motivation

Working memory performance involves complex interactions between sensory representation regions (e.g., FFA for faces, PPA for places) and executive function regions (e.g., DLPFC). While previous studies often focused on memory load effects, this project targets **individual differences in performance** and seeks to **characterize both linear and nonlinear activation patterns** in different brain regions, as performance increases.

## 📊 Methodology

- **Data Source:** HCP task-fMRI data (MSM-All registered, HCP-MMP 1.0 parcellation)
- **Behavioral Measure:** BIS (Balanced Integration Score), combining speed (reaction time) and accuracy
- **Modeling:** Region-wise Generalized Additive Models (GAMs) with BIS as smooth term
- **Covariates:** Age, Gender
- **Correction:** ANOVA-based p-values with False Discovery Rate (FDR) correction
- **Feature Extraction:** First and second derivatives of activation-performance curves to classify region types

## 📁 Project Structure

```
├── BrainMap/
│   └── region_type.R        # Scripts for brain region visualization using ggseg and ggplot2
│
├── GAM/
│   ├── fdr_correction.R     # ANOVA-based testing and FDR correction for GAM models
│   ├── gam_functions.R      # Helper functions for GAM fitting, derivative calculation, etc.
│   └── plot_gam.R           # Visualization of fitted GAM curves and scatter plots
```

## 📌 Dependencies

- R (≥ 4.0)
- `mgcv`
- `ggplot2`
- `ggseg`
- `ggpubr`
- `dplyr`, `purrr`, `tidyr`
- Optional: `fdrtool`, `boot` for significance testing

## 📈 Example Output

- Cortical maps showing region-wise effect size， derivative-based region type classifications
- GAM fit results
- GAM-fitted activation curves by BISw
  

## 🙏 Acknowledgments

- Human Connectome Project (HCP) for data access
- HCP-MMP 1.0 parcellation atlas (Glasser et al.)
- R packages: `mgcv`, `ggseg`, `ggplot2`, etc.

## 📬 Contact

Yuqi Yuan  
Email: u3597424@connect.hku.hk
GitHub: [RichardYuan04](https://github.com/RichardYuan04)

Bohan Zhang
Email:
Github:
---

This repository is part of an ongoing academic study. Please cite appropriately if you use any part of the code or methods.
