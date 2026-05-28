# Mosquito Flight Tracker

An R reconstruction of the Butail et al. (2012) mosquito 3D flight-tracking workflow, plus a
browser-based review app. It takes synchronized stereo camera images, reconstructs each mosquito's
3D flight path, and produces kinematic summaries and figures.

The original system was written in MATLAB (included under `data/`). This repository rebuilds the core
of that workflow in R so it runs without a MATLAB license and is easier to extend.

## Scope (please read)

This is a **Phase 1 reconstruction**, not a finished or validated product.

- The geometry and analysis steps (calibration, stereo matching, triangulation, smoothing, kinematics)
  are faithfully reconstructed.
- The most complex tracking step is a **simplified** version; the original multi-hypothesis and
  occlusion handling are documented but not rebuilt.
- It has **not** been validated against the original MATLAB output or against real camera footage.
  The image-to-3D step has only been run on synthetic images generated for testing.
- The real-trajectory figures show the code runs on real numbers; they are not a biological finding.

## Repository layout

```
README.md                         this file
data/butail_et_al_2012_code/      original MATLAB source (requires MATLAB to run)
data/real_trajectories/           small real Anopheles trajectory samples (CSV)
figures/                          example outputs (PNG/PDF) + the demo data bundle
r_reconstruction/
  R/                              the R pipeline, one file per stage
  run_pipeline.R                  main entry point + synthetic demo
  data_examples/                  runnable demo scripts
  shiny_review/app.R              interactive review app
  tests/                          structural tests
  README.md                       module-level details
```

## Requirements

- R 4.0 or newer.
- R packages:

```r
install.packages(c("tiff", "png", "scatterplot3d", "shiny", "plotly", "DT"))
```

## Setup

```bash
git clone https://github.com/kcring/Mosquito_tracker.git
cd Mosquito_tracker
```

Then start R **from the repository root**. All scripts use paths relative to the root, so set the
working directory there first:

```r
setwd("/path/to/Mosquito_tracker")   # the folder containing this README
```

## How to run

### 1. End-to-end synthetic demo (no images needed)

Simulates 3D tracks, runs the tracking and smoothing pipeline, and writes figures. Good first check
that everything is installed.

```r
source("r_reconstruction/run_pipeline.R")
result <- run_synthetic_demo(n_targets = 5, n_frames = 100)
```

Output: `r_reconstruction/data_examples/demo_output/`.

### 2. Image-to-3D demo (the front half of the pipeline)

Renders synthetic stereo camera images, then runs the real image-processing pipeline on them:
background subtraction, blob detection, stereo matching, triangulation, and tracking. Compares the
recovered 3D positions against the known truth.

```r
source("r_reconstruction/data_examples/synthetic_image_demo.R")
demo_synthetic_images()
```

Output: `figures/synthetic_image_demo/` (detection, triangulation accuracy, and recovered-vs-true
figures, plus the data bundle used by the review app). The input images are synthetic; the processing
code is the real pipeline.

### 3. Real-trajectory demo (the back half of the pipeline)

Runs smoothing, swarm statistics, and plotting on a real *Anopheles coluzzii* trajectory sample that
ships with the repo.

```r
source("r_reconstruction/data_examples/load_real_trajectories.R")
demo_real_data()
```

Output: `figures/real_data_demo/`.

### 4. Interactive review app

A browser-based approximation of the original MATLAB review tool. It loads a prepared data bundle
(calibration + frames + reconstructed tracks) and lets you review tracks against the camera images.

```r
shiny::runApp("r_reconstruction/shiny_review", launch.browser = TRUE)
```

It auto-loads the synthetic demo bundle from `figures/synthetic_image_demo/`. Tabs:

- **Camera overlay**: step through the stereo frames with reconstructed tracks drawn on each image
  (numbered, with trails). Use the frame slider and Prev/Next buttons.
- **3D view**: interactive 3D plot of all tracks.
- **Track table + decisions**: select a track, then Accept / Reject / Join / Swap. Export the
  decisions to a CSV that feeds back into the pipeline.

If you regenerate the bundle (`demo_synthetic_images()` above), the app picks up the new version.

### 5. Structure tests

```r
source("r_reconstruction/tests/test_pipeline_structure.R")
```

### Note: the male-vs-female comparison

`r_reconstruction/data_examples/compare_sex.R` needs the sex-labelled Feugère et al. (2020) recordings,
which are **not** included in this repository (large files, separate license). To run it, download that
dataset (see references below) into a `synthetic_data/` folder at the repo root, then:

```r
source("r_reconstruction/data_examples/compare_sex.R")
compare_sex()
```

## Running on your own real data

To run the full pipeline on real footage you need, per experiment: synchronized stereo image sequences
(left/right), a stereo camera calibration, and the frame rate. See `r_reconstruction/README.md` for the
exact directory layout and the `run_pipeline()` arguments.

## Running the original MATLAB code

The original code in `data/butail_et_al_2012_code/` requires MATLAB and the bundled calibration toolbox.
It is included for reference and validation. The R reconstruction does not depend on it.

## More detail

See `r_reconstruction/README.md` for the module-by-module breakdown and the status of each component.

## References and data sources

- Butail, S., Manoukis, N., Diallo, M., Ribeiro, J.M., Lehmann, T. & Paley, D.A. (2012).
  Reconstructing the flight kinematics of swarming and mating in wild mosquitoes.
  *Journal of the Royal Society Interface*, 9(75), 2624-2638. doi:10.1098/rsif.2012.0150
- Vielma, S. et al. (2025). Sex ratios influence spatial occupancy and kinematic stability of
  *Anopheles coluzzii* mosquito swarms. *Parasites & Vectors*. Dataset: https://osf.io/6nkyq/
- Feugère, L., Gibson, G. & Roux, O. (2020). Audio and 3D flight-track recordings of mosquito
  responses to opposite-sex sound-stimuli. *Dryad*. https://doi.org/10.5061/dryad.9cnp5hqhj
