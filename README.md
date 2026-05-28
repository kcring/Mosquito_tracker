# Mosquito Flight Kinematics Demo Package

This repository is a public-safe subset of the local `mosquito_project` workspace. It is meant to showcase the core code, the R reconstruction layer, and the Shiny review app without publishing internal notes or sensitive handoff materials.

## Included in this repo

- Original MATLAB source used for the Butail-style tracking workflow in `data/butail_et_al_2012_code/` (source folders only)
- R reconstruction modules and pipeline entry points in `r_reconstruction/`
- Interactive Shiny review app in `r_reconstruction/shiny_review/`
- Small sample trajectory data in `data/real_trajectories/`
- Demo figures and example outputs in `figures/`

## Intentionally excluded

- Internal project notes and donor-facing materials in `docs/`, `audit/`, and `deliverables/`
- Local agent instructions in `AGENTS.md`
- Large raw media, audio, and Trackit exports under `synthetic_data/`
- Bundled documentation archives under `data/butail_et_al_2012_code/doc/`

## Quick start

From the repository root in R:

```r
source("r_reconstruction/run_pipeline.R")
source("r_reconstruction/data_examples/synthetic_image_demo.R")

# Synthetic pipeline demo
result <- run_synthetic_demo(n_targets = 5, n_frames = 100)

# Shiny review app
shiny::runApp("r_reconstruction/shiny_review", launch.browser = TRUE)
```

See `r_reconstruction/README.md` for pipeline details, module status, and demo workflows.
