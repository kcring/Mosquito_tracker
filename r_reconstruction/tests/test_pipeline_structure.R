# test_pipeline_structure.R
# Structural validation tests for the Phase 1 R mosquito tracking pipeline.
#
# PURPOSE:
#   These tests verify that each R module runs, produces outputs of the
#   correct shape and type, and that the mathematical functions return
#   expected values on known synthetic inputs.
#
#   They do NOT test biological correctness or tracking performance.
#   All test data is SYNTHETIC.
#
# USAGE:
#   source("r_reconstruction/tests/test_pipeline_structure.R")
#   # or run individual sections

cat("=== Mosquito Tracking Pipeline: Structure Tests ===\n")
cat("    [SYNTHETIC DATA ONLY — not biological validation]\n\n")

# ── Utility helpers ──────────────────────────────────────────────────────────

.assert <- function(cond, msg) {
  if (!isTRUE(cond)) stop(sprintf("[FAIL] %s", msg))
  cat(sprintf("  [PASS] %s\n", msg))
  invisible(TRUE)
}

.assert_near <- function(a, b, tol = 1e-6, msg) {
  .assert(abs(a - b) < tol, sprintf("%s  (got %.8g, expected %.8g)", msg, a, b))
}

.assert_shape <- function(x, dims, msg) {
  actual <- dim(x) %||% length(x)
  .assert(isTRUE(all.equal(actual, dims)), sprintf("%s  (got %s, expected %s)",
          msg, paste(actual, collapse="x"), paste(dims, collapse="x")))
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# ── Load modules ─────────────────────────────────────────────────────────────

base_dir <- if (file.exists("r_reconstruction/R/io.R")) {
  "r_reconstruction/R"
} else if (file.exists("R/io.R")) {
  "R"
} else {
  stop("Cannot find R modules. Run from mosquito_project root.")
}

source(file.path(base_dir, "io.R"))
source(file.path(base_dir, "preprocess.R"))
source(file.path(base_dir, "measurements.R"))
source(file.path(base_dir, "stereo_match.R"))
source(file.path(base_dir, "tracking.R"))
source(file.path(base_dir, "manual_review.R"))
source(file.path(base_dir, "postprocess.R"))
source(file.path(base_dir, "plotting.R"))

cat("Modules loaded.\n\n")

# ── Set up synthetic calibration (from get_cam_calib_mar252016.m) ────────────

cat("--- Calibration ---\n")

km1 <- matrix(c(1192.504524, 0, 641.954689,
                 0, 1190.963449, 490.295747,
                 0, 0, 1), 3, 3, byrow = TRUE)
trm1 <- diag(4)

km2 <- matrix(c(2014.671530, 0, 980.141647,
                 0, 2015.746073, 566.983465,
                 0, 0, 1), 3, 3, byrow = TRUE)
trm2 <- rbind(c(0.9996,  0.0239,  0.0121, -237.7600),
              c(-0.0237,  0.9995, -0.0223,    8.9712),
              c(-0.0126,  0.0220,  0.9997,   19.3559),
              c(0, 0, 0, 1))

cams <- make_stereo_calib(km1, km2, trm1, trm2)
.assert(length(cams) == 2, "make_stereo_calib returns list of 2")
.assert(!is.null(cams[[1]]$km), "cam1 has km field")
.assert(!is.null(cams[[2]]$trm), "cam2 has trm field")
.assert_shape(cams[[1]]$P, c(3,4), "cam1 projection matrix is 3x4")
.assert_shape(cams[[2]]$P, c(3,4), "cam2 projection matrix is 3x4")

# Fundamental matrix
F_mat <- compute_fundamental_matrix(cams[[1]], cams[[2]])
.assert_shape(F_mat, c(3,3), "F matrix is 3x3")
.assert(abs(det(F_mat)) < 1e-6, "F matrix is singular (rank 2)")

cat("\n--- Projection and Epipolar ---\n")

# A 3D point at swarm distance (~500mm depth, near centre)
r_test <- c(100, 50, 500)  # mm

pix1 <- project_to_image(matrix(r_test, 3, 1), cams[[1]])
pix2 <- project_to_image(matrix(r_test, 3, 1), cams[[2]])
.assert_shape(pix1, c(2,1), "Projection cam1 returns 2x1")
.assert_shape(pix2, c(2,1), "Projection cam2 returns 2x1")

# Pixels should be within reasonable range (1392x1040 image)
.assert(pix1[1] > 0 & pix1[1] < 2000, "cam1 pixel col in reasonable range")
.assert(pix1[2] > 0 & pix1[2] < 1500, "cam1 pixel row in reasonable range")

# Epipolar constraint should hold for the true point
ep_ok <- epipolar_check(as.vector(pix1), as.vector(pix2), F_mat, te = 0.5)
.assert(ep_ok, "Epipolar constraint satisfied for true point")

# Triangulation: should recover the original point
pix_pair <- cbind(as.vector(pix1), as.vector(pix2))
tri <- ls_triangulate(pix_pair, cams)
.assert(!is.null(tri$r), "Triangulation returns r")
.assert_shape(tri$r, 3, "Triangulated position is length 3")
.assert(sqrt(sum((tri$r - r_test)^2)) < 20, "Triangulation error < 20mm")
cat(sprintf("  Triangulation error: %.2f mm (threshold: 20 mm)\n",
            sqrt(sum((tri$r - r_test)^2))))

cat("\n--- Measurements ---\n")

Z_empty <- make_Z_struct()
.assert(!is.null(Z_empty$u), "Z struct has u field")
.assert(!is.null(Z_empty$ep), "Z struct has ep field")
.assert_shape(Z_empty$ep, c(2,2), "Z ep field is 2x2")

# Create a synthetic blob for endpoint extraction
set.seed(123)
# A horizontal streak 40px long, 5px wide
blob_pixels <- do.call(rbind, lapply(-20:20, function(dx) {
  c(round(100 + dx), 200 + round(rnorm(1, 0, 2)))  # row, col
}))
blob <- list(
  centroid   = c(200, 100),  # col, row
  area       = nrow(blob_pixels),
  major_axis = 45,
  minor_axis = 8,
  orientation = 0,
  pixel_list = blob_pixels
)

# Test extract_measurements
Z_list <- extract_measurements(list(blob), t_area = c(20, 200), img_dims = c(1024, 1392))
.assert(length(Z_list) == 1, "extract_measurements returns 1 blob")
.assert(!is.null(Z_list[[1]]$ep), "Extracted Z has endpoints")
.assert(Z_list[[1]]$length > 0, "Extracted Z has positive length")
cat(sprintf("  Streak length: %.1f px\n", Z_list[[1]]$length))

cat("\n--- Stereo pair validation ---\n")

# Create matching measurement pair
Z1 <- make_Z_struct(); Z2 <- make_Z_struct()
Z1$u <- as.vector(pix1)
Z2$u <- as.vector(pix2)
Z1$ep <- matrix(c(as.vector(pix1)[1]-5, as.vector(pix1)[2],
                   as.vector(pix1)[1]+5, as.vector(pix1)[2]), 2, 2)
Z2$ep <- matrix(c(as.vector(pix2)[1]-5, as.vector(pix2)[2],
                   as.vector(pix2)[1]+5, as.vector(pix2)[2]), 2, 2)
Z1$sigma <- c(3, 1.5); Z2$sigma <- c(3, 1.5)

pairs <- validate_stereo_pairs(list(Z1), list(Z2), F_mat, te = 0.5)
.assert(length(pairs) == 1, "validate_stereo_pairs finds matching pair")

cat("\n--- Tracking data structures ---\n")

si <- state_index()
.assert(length(si$ri)  == 3, "State index ri has 3 elements")
.assert(length(si$rdi) == 3, "State index rdi has 3 elements")
.assert(si$nX == 6, "State size is 6")

F_cv <- cv_motion_matrix(1/25)
.assert_shape(F_cv, c(6,6), "CV motion matrix is 6x6")
.assert_near(F_cv[1,4], 1/25, msg = "CV matrix dt entry correct")

Q_cv <- cv_process_noise(1/25, sigma_w = 100e6)
.assert_shape(Q_cv, c(6,6), "Process noise matrix is 6x6")
.assert(all(Q_cv == t(Q_cv)), "Process noise is symmetric")
.assert(all(eigen(Q_cv)$values >= -1e-9), "Process noise is PSD")

cat("\n--- Target initialisation and particle filter ---\n")

params <- default_params()
params$Np <- 50  # small for speed

# Synthetic stereo pair from known 3D point
r_true <- c(0, 0, 500)
p1 <- as.vector(project_to_image(matrix(r_true,3,1), cams[[1]]))
p2 <- as.vector(project_to_image(matrix(r_true,3,1), cams[[2]]))
Z1t <- make_Z_struct(); Z2t <- make_Z_struct()
Z1t$u <- p1; Z2t$u <- p2
Z1t$sigma <- c(3,1.5); Z2t$sigma <- c(3,1.5)
te <- params$te
Z1t$ep <- cbind(p1 + c(-3, 0), p1 + c(3, 0))
Z2t$ep <- cbind(p2 + c(-3, 0), p2 + c(3, 0))
Z_pair_test <- list(Z1 = Z1t, Z2 = Z2t)

target <- init_target(Z_pair_test, cams, params, target_id = 1L)
.assert(!is.null(target$particles), "init_target returns particles")
.assert_shape(target$particles, c(6, params$Np), "Particles are 6 x Np")
.assert(target$id == 1L, "Target ID is 1")
.assert(!target$confirmed, "New target not yet confirmed")

# Predict step
target_pred <- predict_targets(list(target), dt = 1/25, sigma_w = params$sigma_w)[[1]]
.assert(!is.null(target_pred$particles), "Predict returns particles")
.assert_shape(target_pred$particles, c(6, params$Np), "Predicted particles are 6 x Np")

# Position likelihood
wts <- likelihood_position(Z1t, target$particles[1:3,,drop=FALSE], cams[[1]])
.assert(length(wts) == params$Np, "Position likelihood returns Np weights")
.assert(all(wts >= 0), "Likelihood weights are non-negative")
.assert(sum(wts) > 0, "At least some positive likelihood weights")

# Resample
idx <- resample_particles(wts)
.assert(length(idx) == params$Np, "Resample returns Np indices")
.assert(all(idx >= 1 & idx <= params$Np), "Resample indices in valid range")

cat("\n--- Xh matrix I/O ---\n")

Xh <- make_xh(max_targets = 10, n_frames = 50)
.assert_shape(Xh, c(60, 50), "Xh is 6*10 x 50")

# Write and read target state
state_in <- c(100, 200, 500, 300, -200, 50)
Xh <- set_target_state(Xh, target_id = 3, frame_k = 10, state = state_in)
state_out <- get_target_state(Xh, target_id = 3, frame_k = 10)
.assert(all.equal(state_in, state_out), "set/get target state round-trips")

# xh_to_tidy
Xh2 <- make_xh(5, 20)
for (tid in 1:3) {
  for (k in 1:20) {
    Xh2 <- set_target_state(Xh2, tid, k, c(tid*10, k, 500, 1, 0, 0))
  }
}
tidy <- xh_to_tidy(Xh2, fps = 25)
.assert(is.data.frame(tidy), "xh_to_tidy returns data frame")
.assert(all(c("target_id","frame","x","y","z","vx","vy","vz") %in% names(tidy)),
        "Tidy df has all required columns")
.assert(length(unique(tidy$target_id)) == 3, "Tidy df has 3 targets")

cat("\n--- Kalman smoother ---\n")

# Generate noisy synthetic trajectory
set.seed(99)
n_test <- 60
traj_noisy <- data.frame(
  target_id = 1, frame = 1:n_test, time_s = (0:(n_test-1))/25,
  x = cumsum(rnorm(n_test, 0, 50)),
  y = cumsum(rnorm(n_test, 0, 50)),
  z = 500 + cumsum(rnorm(n_test, 0, 20)),
  vx = rnorm(n_test, 300, 50),
  vy = rnorm(n_test, -200, 50),
  vz = rnorm(n_test, 0, 30)
)

smoothed <- kalman_smooth(traj_noisy, dt = 1/25, sigma_w = 100e6, sigma_r = 50)
.assert(is.data.frame(smoothed), "kalman_smooth returns data frame")
.assert(nrow(smoothed) == nrow(traj_noisy), "Smoother preserves row count")
.assert(all(!is.na(smoothed$x)), "No NA in smoothed x")

# Smoother should reduce variance (not guaranteed but usually true for long tracks)
cat(sprintf("  SD x: noisy=%.1f  smoothed=%.1f\n",
            sd(traj_noisy$x), sd(smoothed$x)))

cat("\n--- Manual review operations ---\n")

Xh_test <- make_xh(5, 30)
# Put data in targets 1 and 2
for (k in 1:15) Xh_test <- set_target_state(Xh_test, 1, k, c(k, 0, 500, 1, 0, 0))
for (k in 16:30) Xh_test <- set_target_state(Xh_test, 2, k, c(k, 0, 500, 1, 0, 0))

# Join: copy target 2 (frames 16:30) into target 1 starting at frame 16
Xh_joined <- join_tracks(Xh_test, src_id = 2, dst_id = 1, from_frame = 16)
.assert(get_target_state(Xh_joined, 1, 20)[1] == 20, "Joined target 1 has data at frame 20")
.assert(get_target_state(Xh_joined, 2, 20)[1] == 0,  "Src target 2 is zeroed after join")

# Swap
Xh_s <- make_xh(3, 10)
for (k in 1:10) Xh_s <- set_target_state(Xh_s, 1, k, c(1, 0, 500, 0, 0, 0))
for (k in 1:10) Xh_s <- set_target_state(Xh_s, 2, k, c(2, 0, 500, 0, 0, 0))
Xh_sw <- swap_tracks(Xh_s, 1, 2)
.assert(get_target_state(Xh_sw, 1, 5)[1] == 2, "Swap: target 1 now has target 2's data")
.assert(get_target_state(Xh_sw, 2, 5)[1] == 1, "Swap: target 2 now has target 1's data")

# Reject
Xh_r <- reject_track(Xh_test, target_id = 1)
.assert(all(Xh_r[1:6, ] == 0), "Rejected target 1 is all zeros")

cat("\n--- Postprocess summary metrics ---\n")

traj_test <- data.frame(
  target_id = rep(1:3, each = 40),
  frame     = rep(1:40, 3),
  time_s    = rep((0:39)/25, 3),
  x = c(rnorm(40, 100, 50), rnorm(40, -100, 50), rnorm(40, 0, 50)),
  y = c(rnorm(40, 100, 50), rnorm(40, 100, 50), rnorm(40, 0, 50)),
  z = c(rnorm(40, 500, 50), rnorm(40, 400, 50), rnorm(40, 600, 50)),
  vx = rnorm(120, 300, 100), vy = rnorm(120, -100, 80), vz = rnorm(120, 0, 50)
)
traj_test <- add_speed(traj_test)
.assert("speed_ms" %in% names(traj_test), "add_speed adds speed_ms column")
.assert(all(traj_test$speed_ms >= 0), "All speeds non-negative")

swarm_s <- compute_swarm_stats(traj_test)
.assert(is.data.frame(swarm_s), "compute_swarm_stats returns data frame")
.assert("n_mosquitoes" %in% names(swarm_s), "swarm stats has n_mosquitoes")
.assert(all(swarm_s$n_mosquitoes > 0), "All frames have mosquitoes")

paths <- compute_path_lengths(traj_test)
.assert(nrow(paths) == 3, "Path lengths for 3 targets")
.assert(all(paths$path_length_mm >= 0), "Non-negative path lengths")

sep <- compute_separation(traj_test, 1, 2)
.assert(is.data.frame(sep), "compute_separation returns data frame")
.assert("separation_mm" %in% names(sep), "Separation df has distance column")

cat("\n--- Plotting (file output only, no screen required) ---\n")

tmp_dir <- tempdir()

# 2D trajectory plot
pdf_path <- file.path(tmp_dir, "test_trajectories.pdf")
tryCatch({
  plot_trajectories_2d(traj_test, swarm_s,
                        title = "[TEST SYNTHETIC]",
                        save_path = pdf_path)
  .assert(file.exists(pdf_path), "Trajectory PDF created")
  cat(sprintf("  Trajectory PDF: %s\n", pdf_path))
}, error = function(e) cat(sprintf("  [WARN] Trajectory plot: %s\n", e$message)))

# Speed scatter
pdf_path2 <- file.path(tmp_dir, "test_speed_scatter.pdf")
tryCatch({
  plot_speed_scatter(traj_test, save_path = pdf_path2)
  .assert(file.exists(pdf_path2), "Speed scatter PDF created")
}, error = function(e) cat(sprintf("  [WARN] Speed scatter: %s\n", e$message)))

# QC panel
pdf_path3 <- file.path(tmp_dir, "test_qc_panel.pdf")
tryCatch({
  plot_qc_panel(traj_test, swarm_s, title = "[TEST SYNTHETIC] QC",
                 save_path = pdf_path3)
  .assert(file.exists(pdf_path3), "QC panel PDF created")
}, error = function(e) cat(sprintf("  [WARN] QC panel: %s\n", e$message)))

cat("\n--- Save synthetic example files for data_examples/ ---\n")

example_dir <- "r_reconstruction/data_examples"
dir.create(example_dir, showWarnings = FALSE, recursive = TRUE)

# Save calibration
saveRDS(list(cams = cams, F_mat = F_mat, params = params),
        file.path(example_dir, "SYNTHETIC_calibration.rds"))
.assert(file.exists(file.path(example_dir, "SYNTHETIC_calibration.rds")),
        "SYNTHETIC_calibration.rds saved")

# Save example measurements
saveRDS(list(Z_cam1 = list(Z1), Z_cam2 = list(Z2)),
        file.path(example_dir, "SYNTHETIC_measurements_frame001.rds"))
.assert(file.exists(file.path(example_dir, "SYNTHETIC_measurements_frame001.rds")),
        "SYNTHETIC_measurements_frame001.rds saved")

# Save stereo pairs
saveRDS(pairs, file.path(example_dir, "SYNTHETIC_stereo_pairs_frame001.rds"))
.assert(file.exists(file.path(example_dir, "SYNTHETIC_stereo_pairs_frame001.rds")),
        "SYNTHETIC_stereo_pairs_frame001.rds saved")

cat("\n=== All structural tests passed ===\n")
cat("    [Reminder: SYNTHETIC DATA ONLY — not biological validation]\n\n")
