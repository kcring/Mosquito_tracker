# Mosquito Flight Tracker

MATLAB source, an R reconstruction layer, and a Shiny review app for mosquito 3D tracking demos.

## Structure

- `data/butail_et_al_2012_code/` original MATLAB code
- `r_reconstruction/` R pipeline and helpers
- `r_reconstruction/shiny_review/` Shiny review app
- `data/real_trajectories/` sample trajectory data
- `figures/` example outputs

## Quick start

```r
source("r_reconstruction/run_pipeline.R")
source("r_reconstruction/data_examples/synthetic_image_demo.R")

result <- run_synthetic_demo(n_targets = 5, n_frames = 100)
shiny::runApp("r_reconstruction/shiny_review", launch.browser = TRUE)
```

See `r_reconstruction/README.md` for pipeline details.
