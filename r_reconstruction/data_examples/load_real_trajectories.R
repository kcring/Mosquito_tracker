# load_real_trajectories.R
#
# Ingests publicly available Anopheles coluzzii 3D swarm trajectory data
# and runs it through the mosquito project postprocessing + plotting pipeline.
#
# DATA SOURCE:
#   Vielma et al. (2025/2026), "Sex ratios influence spatial occupancy and
#   kinematic stability of Anopheles coluzzii mosquito swarms."
#   Parasites & Vectors. Dataset: https://osf.io/6nkyq/
#
#   File used: group_f200_02.csv  (f200 = 200 females treatment group)
#   A 2 MB sample was downloaded via HTTP byte-range request:
#     data/real_trajectories/Anopheles_coluzzii_swarm_sample.csv
#
# COORDINATE SYSTEM (as per the Trackit/PFMD documentation):
#   x, y: horizontal plane (cm) — converted to mm
#   z:    height above floor (cm) — converted to mm
#   Time resolution: 10 ms per row (100 Hz)
#
# RELATIONSHIP TO BUTAIL PIPELINE:
#   These are *final reconstructed trajectories* (pipeline OUTPUT format).
#   They are NOT raw image data.  However, they allow full validation of
#   every downstream R module: postprocess.R, plotting.R, manual_review.R.
#
# NOTE:
#   - Species:    Anopheles coluzzii (sister species to An. gambiae s.s.)
#   - Setting:    Laboratory-induced swarms (not wild field swarms)
#   - Tracker:    Trackit 3D Fly (SciTrackS GmbH), not the Butail MATLAB tracker
#   - This is real biological data, NOT synthetic.  Swarm kinematics are
#     closely comparable to the Butail (2012) Mali field data.
#
# USAGE:
#   Rscript r_reconstruction/data_examples/load_real_trajectories.R
#   -- or --
#   source("r_reconstruction/data_examples/load_real_trajectories.R")
#   demo_real_data()


# ============================================================
# SETUP
# ============================================================

.find_r_dir <- function() {
  candidates <- c("r_reconstruction/R", "R", "../R")
  for (d in candidates) if (file.exists(file.path(d, "io.R"))) return(d)
  stop("Cannot find R modules. Run from the mosquito_project root directory.")
}

.source_modules <- function(r_dir) {
  for (m in c("io.R", "postprocess.R", "plotting.R")) {
    src <- file.path(r_dir, m)
    if (!file.exists(src)) stop("Missing module: ", src)
    source(src, local = FALSE)
  }
  message("[load_real] Loaded modules: io.R, postprocess.R, plotting.R")
}


# ============================================================
# INGESTION: OSF/Trackit CSV  ->  tidy trajectory data frame
# ============================================================

#' Load and convert OSF mosquito trajectory CSV to tidy format
#'
#' Converts the Trackit-output CSV (id, datetime, x_cm, y_cm, z_cm, ...)
#' into the pipeline's tidy format (target_id, frame, time_s, x, y, z, vx, vy, vz).
#'
#' @param csv_path  path to the CSV file
#' @param fps       target frame rate after resampling (NULL = keep native 100 Hz)
#' @param min_track_frames  minimum number of frames a track must have to be kept
#' @param max_tracks  cap at this many tracks (NULL = all); useful for quick demos
#' @return tidy data frame
load_osf_trajectories <- function(csv_path,
                                   fps              = 100,
                                   min_track_frames = 10,
                                   max_tracks       = NULL) {

  if (!file.exists(csv_path)) stop("CSV not found: ", csv_path)
  message("[load_real] Reading ", basename(csv_path), " ...")

  raw <- read.csv(csv_path, stringsAsFactors = FALSE)

  required_cols <- c("id", "datetime", "x_cm", "y_cm", "z_cm")
  missing <- setdiff(required_cols, names(raw))
  if (length(missing) > 0) stop("Missing columns: ", paste(missing, collapse = ", "))

  # ---- Parse timestamps -> relative seconds ----
  # datetime format: "2022-12-04 01:29:29.083591"
  raw$posix <- as.POSIXct(raw$datetime, format = "%Y-%m-%d %H:%M:%OS", tz = "UTC")
  t0 <- min(raw$posix, na.rm = TRUE)
  raw$time_s_raw <- as.numeric(difftime(raw$posix, t0, units = "secs"))

  # ---- Convert positions: cm -> mm ----
  raw$x_mm <- raw$x_cm * 10
  raw$y_mm <- raw$y_cm * 10
  raw$z_mm <- raw$z_cm * 10

  # ---- Assign integer target IDs ----
  track_ids     <- unique(raw$id)
  id_map        <- setNames(seq_along(track_ids), track_ids)
  raw$target_id <- id_map[raw$id]

  # ---- Filter short tracks ----
  track_lengths <- tapply(raw$target_id, raw$target_id, length)
  keep_ids <- as.integer(names(track_lengths[track_lengths >= min_track_frames]))
  raw <- raw[raw$target_id %in% keep_ids, ]

  # ---- Optionally cap number of tracks (useful for fast demo plots) ----
  if (!is.null(max_tracks) && length(keep_ids) > max_tracks) {
    # keep the longest tracks for best visual output
    tl <- sort(track_lengths[as.character(keep_ids)], decreasing = TRUE)
    keep_ids <- as.integer(names(tl)[seq_len(max_tracks)])
    raw <- raw[raw$target_id %in% keep_ids, ]
  }

  # Re-assign contiguous target IDs after filtering
  track_ids2    <- sort(unique(raw$target_id))
  id_map2       <- setNames(seq_along(track_ids2), track_ids2)
  raw$target_id <- id_map2[as.character(raw$target_id)]

  # ---- Assign frame indices ----
  # Native resolution is 10 ms = 100 Hz.  We assign frame = round(time_s / dt) + 1
  dt <- 1 / fps
  raw$frame <- as.integer(round(raw$time_s_raw / dt)) + 1L

  # ---- Per-track: deduplicate on frame, compute velocities ----
  tidy_list <- lapply(sort(unique(raw$target_id)), function(tid) {
    sub <- raw[raw$target_id == tid, ]
    sub <- sub[order(sub$frame), ]

    # Deduplicate (keep first occurrence per frame)
    sub <- sub[!duplicated(sub$frame), ]

    n <- nrow(sub)
    if (n < 2) return(NULL)

    # Finite-difference velocities (mm/s)
    dt_vec <- diff(sub$time_s_raw)
    dx     <- diff(sub$x_mm)
    dy     <- diff(sub$y_mm)
    dz     <- diff(sub$z_mm)

    # Forward differences; last point copies neighbour
    vx <- c(dx / dt_vec, dx[n-1] / dt_vec[n-1])
    vy <- c(dy / dt_vec, dy[n-1] / dt_vec[n-1])
    vz <- c(dz / dt_vec, dz[n-1] / dt_vec[n-1])

    data.frame(
      target_id = tid,
      frame     = sub$frame,
      time_s    = sub$time_s_raw,
      x  = sub$x_mm,
      y  = sub$y_mm,
      z  = sub$z_mm,
      vx = vx,
      vy = vy,
      vz = vz,
      size_mm2  = if ("size_mm2"  %in% names(sub)) sub$size_mm2  else NA_real_,
      speed_mps = if ("speed_mps" %in% names(sub)) sub$speed_mps else NA_real_
    )
  })

  tidy <- do.call(rbind, Filter(Negate(is.null), tidy_list))
  rownames(tidy) <- NULL

  n_tracks <- length(unique(tidy$target_id))
  n_rows   <- nrow(tidy)
  dur_s    <- diff(range(tidy$time_s))
  message(sprintf("[load_real] Loaded %d tracks, %d position records, %.1f s duration",
                  n_tracks, n_rows, dur_s))
  tidy
}


# ============================================================
# MAIN DEMO FUNCTION
# ============================================================

#' Run the real-data demo
#'
#' Loads the Anopheles coluzzii trajectory sample, smooths, computes stats,
#' and produces figures saved to figures/real_data_demo/.
#'
#' @param csv_path    path to the OSF trajectory CSV sample
#' @param out_dir     directory to save output figures
#' @param max_tracks  max tracks to plot (default 30 for clean figure)
demo_real_data <- function(
  csv_path  = "data/real_trajectories/Anopheles_coluzzii_swarm_sample.csv",
  out_dir   = "figures/real_data_demo",
  max_tracks = 30
) {
  r_dir <- .find_r_dir()
  .source_modules(r_dir)

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # ---------- 1. LOAD --------------------------------------------------
  message("\n=== Step 1: Load real trajectories ===")
  traj <- load_osf_trajectories(csv_path, max_tracks = max_tracks)

  # ---------- 2. SMOOTH ------------------------------------------------
  message("\n=== Step 2: Kalman smooth ===")
  # Native resolution is 100 Hz; paper used 25 Hz.
  # We smooth at native resolution; sigma_r reflects Trackit's ~1mm accuracy
  traj_smooth <- smooth_all_tracks(traj, dt = 1/100, sigma_w = 50e6, sigma_r = 5)

  # ---------- 3. METRICS -----------------------------------------------
  message("\n=== Step 3: Compute metrics ===")
  traj_smooth <- add_speed(traj_smooth)
  swarm_stats <- compute_swarm_stats(traj_smooth)
  path_lengths <- compute_path_lengths(traj_smooth)

  # Print summary
  message(sprintf("  Tracks: %d", length(unique(traj_smooth$target_id))))
  message(sprintf("  Mean speed: %.2f m/s (+/- %.2f)",
                  mean(traj_smooth$speed_ms, na.rm = TRUE),
                  sd(traj_smooth$speed_ms, na.rm = TRUE)))
  message(sprintf("  Mean swarm height: %.0f mm (%.0f cm)",
                  mean(swarm_stats$mean_z, na.rm = TRUE),
                  mean(swarm_stats$mean_z, na.rm = TRUE) / 10))
  message(sprintf("  Swarm centroid range X: %.0f mm, Y: %.0f mm, Z: %.0f mm",
                  diff(range(swarm_stats$mean_x, na.rm = TRUE)),
                  diff(range(swarm_stats$mean_y, na.rm = TRUE)),
                  diff(range(swarm_stats$mean_z, na.rm = TRUE))))

  # ---------- 4. PLOTS -------------------------------------------------
  message("\n=== Step 4: Generate figures ===")

  # 4a. 3D trajectories
  fig3d <- file.path(out_dir, "trajectories_3d.pdf")
  plot_trajectories_3d(
    traj_smooth,
    title     = paste0("Anopheles coluzzii swarm — ", length(unique(traj_smooth$target_id)),
                       " tracks (real data, Vielma et al.)"),
    save_path = fig3d
  )
  message("  Saved: ", fig3d)

  # 4b. 2D projections (x/y/z vs time)
  fig2d <- file.path(out_dir, "trajectories_2d.pdf")
  plot_trajectories_2d(
    traj_smooth,
    title     = "Anopheles coluzzii — real trajectories (x, y, z vs time)",
    save_path = fig2d
  )
  message("  Saved: ", fig2d)

  # 4c. Speed scatter
  fig_speed <- file.path(out_dir, "speed_scatter.pdf")
  plot_speed_scatter(
    traj_smooth,
    title     = "Flight speed distribution (real data)",
    save_path = fig_speed
  )
  message("  Saved: ", fig_speed)

  # 4d. Mosquito count over time
  fig_count <- file.path(out_dir, "mosquito_count.pdf")
  plot_mosquito_count(
    swarm_stats,
    save_path = fig_count
  )
  message("  Saved: ", fig_count)

  # 4e. QC panel
  fig_qc <- file.path(out_dir, "qc_panel.pdf")
  plot_qc_panel(
    traj_smooth, swarm_stats,
    title     = "QC Panel — Anopheles coluzzii real data",
    save_path = fig_qc
  )
  message("  Saved: ", fig_qc)

  # ---------- 5. SAVE TIDY CSV -----------------------------------------
  csv_out <- file.path(out_dir, "tidy_trajectories_smoothed.csv")
  write.csv(traj_smooth, csv_out, row.names = FALSE)
  message("  Saved tidy CSV: ", csv_out)

  csv_stats <- file.path(out_dir, "swarm_stats.csv")
  write.csv(swarm_stats, csv_stats, row.names = FALSE)
  message("  Saved swarm stats: ", csv_stats)

  message("\n=== Demo complete. Output in: ", out_dir, " ===\n")

  invisible(list(
    traj         = traj_smooth,
    swarm_stats  = swarm_stats,
    path_lengths = path_lengths
  ))
}


# ============================================================
# RUN IF CALLED DIRECTLY
# ============================================================

if (!interactive() && identical(environment(), globalenv())) {
  demo_real_data()
}
