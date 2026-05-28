# run_pipeline.R
# Main entry point for the Phase 1 R mosquito tracking pipeline.
#
# This script orchestrates the full pipeline on a single experiment,
# from raw TIFF images (or synthetic data) to smoothed 3D trajectories
# and QC plots.
#
# Mirrors MATLAB source: tracker/trackemall.m (batch runner),
#   tracker/mhtpftracker.m (per-experiment tracking loop),
#   tracker/trackone.m (manual review driver)
#
# USAGE (real data):
#   Rscript run_pipeline.R  -- then edit paths at top of script, or:
#   source("run_pipeline.R")
#
# USAGE (synthetic demo — no images required):
#   source("r_reconstruction/run_pipeline.R")
#   run_synthetic_demo()
#
# NOTE:
#   - Computational geometry (calibration, epipolar, triangulation) is
#     faithfully reproduced from the MATLAB source.
#   - Tracking core (CV-KF + simplified GNN + particle filter weighting)
#     is a Phase 1 simplified version; full MHT and occlusion EM are stubs.
#   - Manual review is CSV-driven (not a MATLAB GUI clone).
#   - No biological/empirical validation; only synthetic data available.

# ============================================================
# LOAD R MODULES
# ============================================================

module_dir <- file.path(dirname(sys.frame(1)$ofile %||% "r_reconstruction"), "R")
if (!dir.exists(module_dir)) module_dir <- "r_reconstruction/R"

source(file.path(module_dir, "io.R"))
source(file.path(module_dir, "preprocess.R"))
source(file.path(module_dir, "measurements.R"))
source(file.path(module_dir, "stereo_match.R"))
source(file.path(module_dir, "tracking.R"))
source(file.path(module_dir, "manual_review.R"))
source(file.path(module_dir, "postprocess.R"))
source(file.path(module_dir, "plotting.R"))

`%||%` <- function(a, b) if (is.null(a)) b else a

# ============================================================
# MAIN PIPELINE FUNCTION
# ============================================================

#' Run the full tracking pipeline on a single experiment
#'
#' @param exp_dir  experiment root directory (contains calib/, frames/)
#' @param frame_dir  directory with TIFF frames; NULL = exp_dir/frames/
#' @param params  parameter list from default_params() or modified version
#' @param run_review  logical: pause for CSV review before smoothing (default TRUE)
#' @param review_csv  path for review CSV; NULL = auto-named in data_dir
#' @return list with traj_df (smoothed), Xh_auto, Xh_reviewed, swarm_stats
run_pipeline <- function(exp_dir,
                          frame_dir  = NULL,
                          params     = default_params(),
                          run_review = TRUE,
                          review_csv = NULL) {

  cat("\n=== Mosquito 3D Tracking Pipeline (Phase 1 R) ===\n")
  cat(sprintf("Experiment: %s\n", exp_dir))

  # ── Stage 1: Initialise ─────────────────────────────────
  cat("\n[1/7] Initialising experiment...\n")
  exp <- init_experiment(exp_dir, frame_dir, params)
  cat(sprintf("  Cameras: %d | Frames: %d per camera\n",
              params$nc, length(exp$frame_lists[[1]])))

  # ── Stage 2: Per-frame tracking ─────────────────────────
  cat("\n[2/7] Running automated tracking...\n")

  n_frames <- min(sapply(exp$frame_lists, length))
  Xh <- make_xh(params$max_targets, n_frames)
  targets <- list()
  next_id <- 1L

  for (k in seq_len(n_frames)) {
    if (k %% 25 == 0) cat(sprintf("  Frame %d/%d...\n", k, n_frames))

    # Load image windows for each camera
    Z_all <- lapply(seq_len(params$nc), function(c) {
      imgarr <- tryCatch(
        load_frame_window(exp$frame_lists[[c]], k, params$br0),
        error = function(e) NULL
      )
      if (is.null(imgarr)) return(list())
      bg_res  <- compute_background(imgarr, params$br0, params$br0, params$fg_is_dark)
      fg      <- compute_foreground(bg_res$bg, bg_res$current,
                                     exp$params$binary_t %||% params$binary_t,
                                     params$fg_is_dark)
      blobs   <- extract_blobs(fg)
      extract_measurements(blobs, params$t_area,
                           img_dims = dim(imgarr)[1:2])
    })

    # Stereo pair validation
    Z_pairs <- if (length(Z_all[[1]]) > 0 && length(Z_all[[2]]) > 0)
      validate_stereo_pairs(Z_all[[1]], Z_all[[2]], exp$F_mat, params$te_epipolar)
    else
      list()

    # Per-frame update
    result  <- update_and_initialize(targets, Z_pairs, exp$cams, params, next_id)
    targets <- result$targets
    next_id <- result$next_id

    # Write to Xh
    Xh <- write_targets_to_xh(Xh, targets, k)
  }

  # Save auto-tracked output
  save_tracklets(Xh, exp$expname, exp$data_dir, exp$cams, params, label = "auto")
  Xh_auto <- Xh

  # ── Stage 3: Manual review ──────────────────────────────
  cat("\n[3/7] Manual review...\n")
  traj_auto <- xh_to_tidy(Xh_auto, params$fps)
  summary_df <- show_track_summary(traj_df = traj_auto)

  Xh_reviewed <- Xh_auto

  if (run_review && nrow(traj_auto) > 0) {
    if (is.null(review_csv)) {
      review_csv <- file.path(exp$data_dir,
                              sprintf("review_%s.csv", exp$expname))
    }
    create_review_table(unique(traj_auto$target_id), review_csv)
    cat(sprintf("\n  [!] Review CSV written to:\n      %s\n", review_csv))
    cat("  Fill in decisions, then call:\n")
    cat("    decisions <- load_review_decisions('", review_csv, "')\n", sep="")
    cat("    Xh_reviewed <- apply_review_decisions(Xh_auto, decisions)\n")
    cat("  Then re-run from Stage 4 (or set run_review=FALSE to skip).\n")
  }

  # ── Stage 4: Kalman smoothing ───────────────────────────
  cat("\n[4/7] Kalman smoothing...\n")
  traj_reviewed <- xh_to_tidy(Xh_reviewed, params$fps)
  traj_smoothed <- tryCatch(
    smooth_all_tracks(traj_reviewed, dt = 1/params$fps,
                      sigma_w = params$sigma_w),
    error = function(e) {
      warning("[run_pipeline.R] Smoother failed: ", e$message)
      traj_reviewed
    }
  )
  traj_smoothed <- add_speed(traj_smoothed)
  save_tracklets(Xh_reviewed, exp$expname, exp$data_dir,
                 exp$cams, params, label = "")

  # ── Stage 5: Summary metrics ────────────────────────────
  cat("\n[5/7] Computing summary metrics...\n")
  swarm_stats <- compute_swarm_stats(traj_smoothed)
  path_lengths <- compute_path_lengths(traj_smoothed)
  cat(sprintf("  Tracks: %d | Position points: %d\n",
              nrow(path_lengths), nrow(traj_smoothed)))

  # ── Stage 6: QC plots ───────────────────────────────────
  cat("\n[6/7] Generating plots...\n")
  fig_dir <- file.path(dirname(exp$data_dir), "figures")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

  qc_pdf <- file.path(fig_dir, sprintf("qc_%s.pdf", exp$expname))
  tryCatch(
    plot_qc_panel(traj_smoothed, swarm_stats,
                   title = sprintf("QC: %s", exp$expname),
                   save_path = qc_pdf),
    error = function(e) warning("[run_pipeline.R] QC plot failed: ", e$message)
  )
  cat(sprintf("  QC plot: %s\n", qc_pdf))

  traj_pdf <- file.path(fig_dir, sprintf("trajectories_%s.pdf", exp$expname))
  tryCatch(
    plot_trajectories_2d(traj_smoothed, swarm_stats,
                          title = sprintf("Trajectories: %s", exp$expname),
                          save_path = traj_pdf),
    error = function(e) warning("[run_pipeline.R] Trajectory plot failed: ", e$message)
  )
  cat(sprintf("  Trajectory plot: %s\n", traj_pdf))

  # ── Stage 7: Save final CSV ─────────────────────────────
  cat("\n[7/7] Saving final outputs...\n")
  final_csv <- file.path(exp$data_dir,
                          sprintf("trajectories_smoothed_%s.csv", exp$expname))
  write.csv(traj_smoothed, final_csv, row.names = FALSE)
  swarm_csv <- file.path(exp$data_dir,
                          sprintf("swarm_stats_%s.csv", exp$expname))
  write.csv(swarm_stats, swarm_csv, row.names = FALSE)
  cat(sprintf("  Trajectories: %s\n", final_csv))
  cat(sprintf("  Swarm stats:  %s\n", swarm_csv))

  cat("\n=== Pipeline complete ===\n")

  invisible(list(
    traj_df     = traj_smoothed,
    swarm_stats = swarm_stats,
    path_lengths = path_lengths,
    Xh_auto     = Xh_auto,
    Xh_reviewed = Xh_reviewed,
    exp         = exp
  ))
}


# ============================================================
# SYNTHETIC DEMO (no images required)
# ============================================================

#' Run the full pipeline on synthetic data to demonstrate structure
#'
#' Generates synthetic stereo tracks, runs tracking, smoothing, and plots.
#' All output is clearly labelled SYNTHETIC.
#'
#' @param n_targets  number of synthetic mosquitoes (default 5)
#' @param n_frames  number of frames (default 100 ~ 4 s at 25 fps)
#' @param out_dir  output directory for results
#' @return same list as run_pipeline()
run_synthetic_demo <- function(n_targets = 5, n_frames = 100,
                                out_dir = "r_reconstruction/data_examples/demo_output") {
  cat("\n=== SYNTHETIC DEMO (Phase 1 R Reconstruction) ===\n")
  cat("  [!] All data below is SYNTHETIC and for structural validation only.\n")
  cat("      These are NOT real mosquito trajectories.\n\n")

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # ── Generate synthetic calibration ────────────────────
  source("r_reconstruction/R/io.R", local = TRUE)
  source("r_reconstruction/R/stereo_match.R", local = TRUE)

  # Realistic calibration from get_cam_calib_mar252016.m
  km1 <- matrix(c(1192.5, 0, 641.9,
                   0, 1191.0, 490.3,
                   0, 0, 1), 3, 3, byrow = TRUE)
  trm1 <- diag(4)
  km2 <- matrix(c(2014.7, 0, 980.1,
                   0, 2015.7, 567.0,
                   0, 0, 1), 3, 3, byrow = TRUE)
  trm2 <- rbind(c(0.9996, 0.0239, 0.0121, -237.76),
                c(-0.0237, 0.9995, -0.0223, 8.9712),
                c(-0.0126, 0.0220, 0.9997, 19.356),
                c(0, 0, 0, 1))
  cams <- make_stereo_calib(km1, km2, trm1, trm2)
  F_mat <- compute_fundamental_matrix(cams[[1]], cams[[2]])

  # ── Generate synthetic ground-truth tracks ────────────
  params <- default_params()
  params$fps <- 25
  set.seed(42)

  dt <- 1/params$fps
  frames <- seq_len(n_frames)
  gt_tracks <- list()

  for (tid in seq_len(n_targets)) {
    # Start near swarm centre at ~1.5m depth, 0.5m height
    r0  <- c(rnorm(1, 0, 200), rnorm(1, 0, 200), rnorm(1, 500, 100))
    v0  <- c(rnorm(1, 0, 600), rnorm(1, 0, 600), rnorm(1, 0, 300))

    traj <- matrix(0, 6, n_frames)
    traj[, 1] <- c(r0, v0)
    for (k in 2:n_frames) {
      w <- rnorm(6, 0, sqrt(c(rep(params$sigma_w * dt^3/3, 3),
                               rep(params$sigma_w * dt, 3))))
      traj[, k] <- c(
        traj[1:3, k-1] + traj[4:6, k-1] * dt + w[1:3],
        traj[4:6, k-1] + w[4:6]
      )
      # Reflect at swarm boundaries
      traj[3, k] <- max(100, min(900, traj[3, k]))
    }

    gt_tracks[[tid]] <- data.frame(
      target_id = tid, frame = frames, time_s = (frames-1)/25,
      x = traj[1,], y = traj[2,], z = traj[3,],
      vx = traj[4,], vy = traj[5,], vz = traj[6,]
    )
  }
  traj_gt <- do.call(rbind, gt_tracks)
  cat(sprintf("  Generated %d synthetic tracks over %d frames.\n", n_targets, n_frames))

  # ── Project to pixel observations (simulate measurements) ─
  Z_pairs_all <- lapply(frames, function(k) {
    pairs <- list()
    for (tid in seq_len(n_targets)) {
      r_world <- c(traj_gt$x[traj_gt$target_id==tid & traj_gt$frame==k],
                   traj_gt$y[traj_gt$target_id==tid & traj_gt$frame==k],
                   traj_gt$z[traj_gt$target_id==tid & traj_gt$frame==k])
      v_world <- c(traj_gt$vx[traj_gt$target_id==tid & traj_gt$frame==k],
                   traj_gt$vy[traj_gt$target_id==tid & traj_gt$frame==k],
                   traj_gt$vz[traj_gt$target_id==tid & traj_gt$frame==k])
      te <- params$te
      u1 <- project_to_image(matrix(r_world, 3, 1), cams[[1]]) + rnorm(2, 0, 2)
      u2 <- project_to_image(matrix(r_world, 3, 1), cams[[2]]) + rnorm(2, 0, 2)

      Z1 <- make_Z_struct(); Z2 <- make_Z_struct()
      Z1$u <- as.vector(u1); Z2$u <- as.vector(u2)
      Z1$sigma <- c(3, 1.5); Z2$sigma <- c(3, 1.5)

      ep1_world <- r_world - v_world * te / 2
      ep2_world <- r_world + v_world * te / 2
      Z1$ep <- cbind(as.vector(project_to_image(matrix(ep1_world,3,1), cams[[1]])),
                     as.vector(project_to_image(matrix(ep2_world,3,1), cams[[1]])))
      Z2$ep <- cbind(as.vector(project_to_image(matrix(ep1_world,3,1), cams[[2]])),
                     as.vector(project_to_image(matrix(ep2_world,3,1), cams[[2]])))
      Z1$length <- 5; Z2$length <- 5

      if (epipolar_check(Z1$u, Z2$u, F_mat, te = params$te_epipolar)) {
        pairs[[length(pairs)+1]] <- list(Z1=Z1, Z2=Z2)
      }
    }
    pairs
  })

  # ── Run simplified tracking ────────────────────────────
  cat("  Running tracking on synthetic measurements...\n")
  Xh <- make_xh(params$max_targets, n_frames)
  targets <- list()
  next_id <- 1L

  for (k in frames) {
    result  <- update_and_initialize(targets, Z_pairs_all[[k]], cams, params, next_id)
    targets <- result$targets
    next_id <- result$next_id
    Xh <- write_targets_to_xh(Xh, targets, k)
  }

  # ── Smooth ────────────────────────────────────────────
  traj_auto     <- xh_to_tidy(Xh, fps = 25)
  traj_smoothed <- if (nrow(traj_auto) > 0)
    smooth_all_tracks(traj_auto) else traj_auto
  traj_smoothed <- add_speed(traj_smoothed)
  swarm_stats   <- compute_swarm_stats(traj_smoothed)

  # ── Save ──────────────────────────────────────────────
  write.csv(traj_gt,       file.path(out_dir, "SYNTHETIC_ground_truth.csv"),       row.names=FALSE)
  write.csv(traj_smoothed, file.path(out_dir, "SYNTHETIC_tracked_smoothed.csv"),   row.names=FALSE)
  write.csv(swarm_stats,   file.path(out_dir, "SYNTHETIC_swarm_stats.csv"),        row.names=FALSE)
  saveRDS(list(cams=cams, F_mat=F_mat, params=params), file.path(out_dir, "SYNTHETIC_calibration.rds"))

  # ── Plot ──────────────────────────────────────────────
  fig_dir <- file.path(out_dir, "figures")
  dir.create(fig_dir, showWarnings = FALSE)

  tryCatch(plot_trajectories_2d(traj_gt,
    title = "[SYNTHETIC] Ground-truth trajectories",
    save_path = file.path(fig_dir, "SYNTHETIC_ground_truth_trajectories.pdf")),
    error = function(e) warning("GT plot: ", e$message))

  tryCatch(plot_trajectories_2d(traj_smoothed, swarm_stats,
    title = "[SYNTHETIC] Tracked + smoothed trajectories",
    save_path = file.path(fig_dir, "SYNTHETIC_tracked_trajectories.pdf")),
    error = function(e) warning("Track plot: ", e$message))

  tryCatch(plot_qc_panel(traj_smoothed, swarm_stats,
    title = "[SYNTHETIC] QC Overview",
    save_path = file.path(fig_dir, "SYNTHETIC_qc_panel.pdf")),
    error = function(e) warning("QC plot: ", e$message))

  cat(sprintf("\n  Outputs in: %s\n", out_dir))
  cat("  REMINDER: All tracks above are SYNTHETIC. Not real mosquito data.\n")
  cat("=== Demo complete ===\n\n")

  invisible(list(
    traj_gt       = traj_gt,
    traj_smoothed = traj_smoothed,
    swarm_stats   = swarm_stats,
    cams          = cams,
    F_mat         = F_mat
  ))
}
