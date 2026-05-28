# R Reconstruction - Mosquito 3D Flight-Tracking
## Phase 1 Package

> **Note**: This is a Phase 1 reconstruction, not a production system.
> See the status table below for what is implemented, partial, or deferred.
> The back half of the pipeline (smoothing, metrics, plots) is validated on **real**
> *Anopheles coluzzii* trajectory data; the front half (image → 3D) is validated on
> **synthetic stereo images** because no real raw image sequences are public.

---

## What is this?

A modular R reconstruction of the MATLAB mosquito 3D tracking pipeline described in:

> Butail, S., Manoukis, N., Diallo, M., Ribeiro, J.M., Lehmann, T. & Paley, D.A. (2012).
> Reconstructing the flight kinematics of swarming and mating in wild mosquitoes.
> *Journal of the Royal Society Interface*, 9(75), 2624–2638.
> doi:10.1098/rsif.2012.0150

The original MATLAB code is at `data/butail_et_al_2012_code/`.

---

## Option A or B?

**Option A (aggressive hybrid)**. This public repo preserves the original MATLAB source,
reconstructs the deterministic geometry and reporting pipeline in R, and uses a lightweight
Shiny/CSV review layer instead of claiming a full end-to-end rewrite.

The deterministic geometry pipeline (calibration I/O, foreground extraction, streak measurements, epipolar validation, triangulation, gating, Kalman smoothing, trajectory postprocessing) is faithfully reconstructed. The MHT and particle filter are simplified and clearly labeled. The manual review GUI is approximated by a CSV workflow. See the status table below.

---

## Module overview

| File | Role | MATLAB equivalent | Status |
|---|---|---|---|
| `R/io.R` | Configuration, calibration I/O, Xh matrix I/O, climate data | `config.m`, `readOffCamCalib.m`, `t1_save_tracks.m` | translated |
| `R/preprocess.R` | Frame loading, sliding-window background, foreground extraction, blob detection | `getZ.m`, `init_bgparams.m` | translated |
| `R/measurements.R` | Streak model, endpoint extraction, missing-measurement search | `setEndPointVelocities.m`, `adaptive_thresholding.m` | translated |
| `R/stereo_match.R` | Fundamental matrix, epipolar check, DLT triangulation, projection, gating | `get_F_for_stereo.m`, `lsTriangulate.m`, `w2cam_nd.m`, `vclz.m` | translated |
| `R/tracking.R` | CV motion model, particle filter, simplified GNN clustering, per-frame update, MHT/occlusion stubs | `mhtpftracker.m`, `predict1.m`, `p_lfn.m`, `pf_update.m`, `clustr.m` | partial (MHT/occlusion deferred) |
| `R/manual_review.R` | CSV review workflow, join/swap/reject operations | `trackone.m`, `glueTracks.m` | partial (CSV workflow) |
| `shiny_review/app.R` | Interactive review app: camera-image overlay + 3D view + table + accept/reject/join/swap → CSV | `trackone.m` GUI | translated (synthetic images) |
| `R/postprocess.R` | Kalman smoother (RTS), trajectory tidying, swarm metrics, re-triangulation | `filter_and_smooth.m`, `nnf_kf.m`, `compareAutoWithManual.m` | translated |
| `R/plotting.R` | 3D trajectories, speed plots, QC panel, mating event, climate data | `plot_data.m`, `plotMatingData.m`, `compareAutoWithManual.m` | translated |
| `run_pipeline.R` | Orchestrates everything; `run_synthetic_demo()` for testing | `trackemall.m` | partial |

---

## How to run

### Quick start: synthetic demo (no images required)

```r
setwd("path/to/mosquito_project")
source("r_reconstruction/R/io.R")
source("r_reconstruction/R/preprocess.R")
source("r_reconstruction/R/measurements.R")
source("r_reconstruction/R/stereo_match.R")
source("r_reconstruction/R/tracking.R")
source("r_reconstruction/R/manual_review.R")
source("r_reconstruction/R/postprocess.R")
source("r_reconstruction/R/plotting.R")
source("r_reconstruction/run_pipeline.R")

result <- run_synthetic_demo(n_targets = 5, n_frames = 100)
```

Output goes to `r_reconstruction/data_examples/demo_output/`.

### Image → 3D demo (renders synthetic stereo images, runs the full front half)

```r
source("r_reconstruction/data_examples/synthetic_image_demo.R")
demo_synthetic_images()   # → figures/synthetic_image_demo/*.pdf
```

Renders motion-blurred mosquito streaks onto two noisy camera images (PNG on disk),
then runs background subtraction → blob detection → epipolar matching → triangulation →
tracking. Recovers ~98% of mosquitoes at ~1 mm median 3D accuracy. The **input images
are synthetic**, but every processing step is the real pipeline code.

### Real-data demos (run on real *Anopheles coluzzii* trajectories)

```r
# Swarm trajectories + smoothing + plots
source("r_reconstruction/data_examples/load_real_trajectories.R"); demo_real_data()

# Male vs Female kinematics comparison
source("r_reconstruction/data_examples/compare_sex.R"); compare_sex()
```

### Interactive review app (Shiny) - `trackone` replica

```r
# 1) generate the image-review bundle (once)
source("r_reconstruction/data_examples/synthetic_image_demo.R"); demo_synthetic_images()
# 2) launch the app (auto-loads the bundle)
shiny::runApp("r_reconstruction/shiny_review", launch.browser = TRUE)
```

This now reproduces the core of the MATLAB `trackone.m` review GUI:

- **Camera overlay tab**: step through the synchronised stereo frames with the
  reconstructed 3D tracks **projected back onto each camera image** (id-labelled,
  with track trails) plus ground-truth markers - the verify/combine review.
- Interactive 3D track view, sortable track table, and accept/reject/join/swap that
  write a review CSV.

The overlay images are synthetic (rendered by `synthetic_image_demo.R`); point the
bundle at real calibrated footage and the same overlay works unchanged.

### Run the structure test

```r
source("r_reconstruction/tests/test_pipeline_structure.R")
```

### Run on real data (when available)

```r
source("r_reconstruction/run_pipeline.R")

result <- run_pipeline(
  exp_dir   = "/path/to/experiment/",
  frame_dir = "/path/to/experiment/frames/",
  params    = default_params(),
  run_review = TRUE
)
```

**Prerequisites for real data**:
1. A directory with TIFF (or PNG) stereo image sequences, named `L*.tif` and `R*.tif`
2. A `calib/` subdirectory containing a Bouguet-calibration `.m` file (e.g. `calib_20100829.m`)
3. A `calib/expfile.txt` with two fields: `expname` and `image_id` (see `io.R::read_expfile`)
4. Install `tiff` and optionally `scatterplot3d` packages:
   ```r
   install.packages(c("tiff", "scatterplot3d"))
   ```

### Manual review workflow

After automated tracking produces a review CSV:

```r
# The pipeline writes: data_dir/review_<expname>.csv
# Fill in the 'decision' column in a spreadsheet, then:

decisions <- load_review_decisions("path/to/review_exp.csv")
Xh_reviewed <- apply_review_decisions(Xh_auto, decisions)
```

Supported decisions: `accept`, `reject`, `join` (requires `join_with` + `join_from_frame`), `swap` (requires `swap_with`).

---

## What remains dependent on MATLAB / manual review

| Component | Status | Notes |
|---|---|---|
| Original MATLAB tracking | MATLAB required | All source at `data/butail_et_al_2012_code/` |
| `trackone.m` GUI review | Approximated | Interactive Shiny app in `shiny_review/app.R`: camera-image overlay (on synthetic frames) + 3D view + CSV workflow |
| Bouguet calibration toolbox | MATLAB required | `tracker/TOOLBOX_calib/` included in source |
| Full MHT (Murty k-best) | Deferred to Phase 2 | `tracker/murtykbest.m` exists but not ported |
| Occlusion EM splitting | Deferred to Phase 2 | `tracker/emgm.m` exists; requires pixel arrays |
| Real image sequences | Not available | No TIFF data public; `synthetic_image_demo.R` renders stand-in images and runs the full front half |
| Real calibration files | Not available | `get_cam_calib_mar252016.m` is the only example |

---

## Folder structure

```
r_reconstruction/
├── README.md                 ← this file
├── run_pipeline.R            ← main entry point + synthetic demo
├── R/
│   ├── io.R                  ← calibration, I/O, config
│   ├── preprocess.R          ← foreground extraction
│   ├── measurements.R        ← streak measurements
│   ├── stereo_match.R        ← epipolar, triangulation, gating
│   ├── tracking.R            ← motion model, PF, clustering, MHT stubs
│   ├── manual_review.R       ← CSV review workflow, Shiny spec
│   ├── postprocess.R         ← Kalman smoother, metrics
│   └── plotting.R            ← visualisation
├── data_examples/
│   ├── README.md             ← data provenance documentation
│   └── demo_output/          ← generated by run_synthetic_demo()
└── tests/
    └── test_pipeline_structure.R  ← structural validation test
```

---

## Key data format

The central data structure is the `Xh` matrix, stored as RDS and CSV:

- **Xh matrix** (MATLAB-compatible format): `[6 * Nmax_targets × Nframes]`
  - Rows `6(t-1)+1 : 6t` = target `t`: `[x, y, z, vx, vy, vz]` (mm, mm/s)
  - Value `0` = untracked frame

- **Tidy data frame** (R-native): one row per tracked (target, frame) observation
  - Columns: `target_id, frame, time_s, x, y, z, vx, vy, vz`

Functions for converting between formats: `xh_to_tidy()`, `get_target_state()`, `set_target_state()`.

---

## Coordinate system

Follows the original MATLAB system (from paper section 4.3):
- Origin: ground level under the camera rig
- x: East–West (mm)
- y: South–North (mm)
- z: Vertical (mm, positive up)
- Velocities in the same frame (mm/s)

---

## References

- Butail et al. (2012). J. R. Soc. Interface. doi:10.1098/rsif.2012.0150
- Hartley & Zisserman (2004). Multiple View Geometry in Computer Vision. Cambridge UP.
- Reid (1979). IEEE Trans. Autom. Control. (MHT algorithm)
- Gordon, Salmond & Smith (1993). IEEE Proc. Radar Signal Process. (particle filter)
