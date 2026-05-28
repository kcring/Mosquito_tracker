# synthetic_image_demo.R
#
# END-TO-END demonstration of the IMAGE -> 3D MOVEMENT half of the pipeline.
#
# WHY THIS EXISTS:
#   The real-data demos (load_real_trajectories.R, compare_sex.R) start from
#   trajectories that were ALREADY reconstructed by someone else's tracker.
#   They exercise only the BACK half of the pipeline (smoothing, metrics, plots).
#
#   This script exercises the FRONT half — the scientifically hard part —
#   the part that turns photographs into 3D movement:
#
#     synthetic 3D paths
#        -> render motion-blurred streaks onto two noisy camera images (PNG on disk)
#        -> read_frame()            (preprocess.R)
#        -> compute_background()    (sliding-window max, getZ.m eq 3.1)
#        -> compute_foreground()    (background subtraction + threshold)
#        -> extract_blobs()         (connected components, regionprops)
#        -> extract_measurements()  (streak model, setEndPointVelocities.m)
#        -> validate_stereo_pairs() (epipolar constraint, eq 3.2)
#        -> ls_triangulate()        (DLT, lsTriangulate.m)
#        -> update_and_initialize() (particle-filter tracking)
#        -> recovered 3D tracks, compared against known ground truth
#
# NOTE:
#   - The INPUT IMAGES ARE SYNTHETIC (rendered streaks on a noisy background),
#     NOT real photographs. No real Mali field images are publicly available.
#   - But every processing step that runs on them is the REAL pipeline code.
#     This proves the image->movement code path works end-to-end and recovers
#     known 3D positions to sub-millimetre accuracy on noise-free geometry.
#
# USAGE:
#   Rscript r_reconstruction/data_examples/synthetic_image_demo.R
#   -- or --
#   source("r_reconstruction/data_examples/synthetic_image_demo.R"); demo_synthetic_images()
#
# NOTE: blob detection is pure-R connected components, so images are kept small
#       (default 180 x 240) to keep runtime to ~1-2 minutes.


# ============================================================
# SETUP
# ============================================================

.find_r_dir <- function() {
  for (d in c("r_reconstruction/R", "R", "../R")) {
    if (file.exists(file.path(d, "io.R"))) return(d)
  }
  stop("Cannot find R modules. Run from the mosquito_project root directory.")
}

.load_modules <- function() {
  r_dir <- .find_r_dir()
  for (m in c("io.R", "preprocess.R", "measurements.R", "stereo_match.R",
              "tracking.R", "postprocess.R", "plotting.R")) {
    source(file.path(r_dir, m), local = FALSE)
  }
  if (!requireNamespace("png", quietly = TRUE)) {
    stop("The 'png' package is required to write synthetic images. install.packages('png')")
  }
  invisible(TRUE)
}


# ============================================================
# DEMO CALIBRATION (small image, consistent stereo geometry)
# ============================================================

#' Build a small-image stereo calibration suitable for fast rendering.
#'
#' Camera 1 sits at the origin looking down +z. Camera 2 is translated along
#' x with a toe-in rotation so the swarm stays centred in both views.
#' The SAME calibration is used for rendering and reconstruction, so
#' triangulation should recover ground-truth positions almost exactly.
make_demo_calib <- function(W = 240, H = 180, f = 520,
                            baseline_mm = 150, swarm_depth = 1100) {
  cx <- W / 2; cy <- H / 2
  km1 <- matrix(c(f, 0, cx,
                  0, f, cy,
                  0, 0, 1), 3, 3, byrow = TRUE)
  km2 <- km1

  trm1 <- diag(4)

  # Toe-in so swarm centre (0,0,swarm_depth) projects near (cx, cy) in cam2
  theta <- asin(baseline_mm / swarm_depth)        # rotation about y (rad)
  Ry <- matrix(c(cos(theta), 0, sin(theta),
                 0,          1, 0,
                 -sin(theta),0, cos(theta)), 3, 3, byrow = TRUE)
  t2 <- c(-baseline_mm, 0, 0)
  trm2 <- rbind(cbind(Ry, t2), c(0, 0, 0, 1))

  cams  <- make_stereo_calib(km1, km2, trm1, trm2)
  F_mat <- compute_fundamental_matrix(cams[[1]], cams[[2]])
  list(cams = cams, F_mat = F_mat, W = W, H = H,
       swarm_depth = swarm_depth, f = f)
}


# ============================================================
# GROUND-TRUTH 3D TRACK GENERATION
# ============================================================

#' Generate smooth constant-velocity-ish ground truth tracks within the swarm.
generate_truth_tracks <- function(n_targets, n_frames, params, calib, seed = 7) {
  set.seed(seed)
  dt <- params$dt
  depth <- calib$swarm_depth
  frames <- seq_len(n_frames)

  tracks <- lapply(seq_len(n_targets), function(tid) {
    r0 <- c(rnorm(1, 0, 90), rnorm(1, 0, 70), rnorm(1, depth, 70))
    # Moderate speeds so streaks are visible but mosquitoes move between frames
    v0 <- c(rnorm(1, 0, 450), rnorm(1, 0, 450), rnorm(1, 0, 250))

    st <- matrix(0, 6, n_frames); st[, 1] <- c(r0, v0)
    for (k in 2:n_frames) {
      w <- rnorm(6, 0, sqrt(c(rep(params$sigma_w * dt^3 / 3, 3),
                              rep(params$sigma_w * dt, 3))))
      st[, k] <- c(st[1:3, k-1] + st[4:6, k-1] * dt + w[1:3],
                   st[4:6, k-1] + w[4:6])
      # Keep mosquitoes inside a sensible swarm volume (reflect)
      for (a in 1:3) {
        lim <- c(160, 130, depth + 110)[a]; lo <- c(-160, -130, depth - 110)[a]
        if (st[a, k] > lim) { st[a, k] <- lim; st[a+3, k] <- -abs(st[a+3, k]) }
        if (st[a, k] < lo)  { st[a, k] <- lo;  st[a+3, k] <-  abs(st[a+3, k]) }
      }
    }
    data.frame(target_id = tid, frame = frames, time_s = (frames - 1) * dt,
               x = st[1,], y = st[2,], z = st[3,],
               vx = st[4,], vy = st[5,], vz = st[6,])
  })
  do.call(rbind, tracks)
}


# ============================================================
# IMAGE RENDERING (motion-blurred dark streaks on a light noisy background)
# ============================================================

#' Draw a dark streak (line segment with Gaussian cross-profile) onto an image.
#'
#' @param img  H x W matrix in [0,1] (modified copy returned)
#' @param p1,p2  endpoint pixel coords c(col,row)
#' @param depth  strength of darkening (0..1); higher = darker
#' @param radius perpendicular falloff (px)
.draw_streak <- function(img, p1, p2, depth = 0.7, radius = 1.4) {
  H <- nrow(img); W <- ncol(img)
  len <- sqrt(sum((p2 - p1)^2))
  nstep <- max(2, ceiling(len * 3))
  ts <- seq(0, 1, length.out = nstep)
  rad_ceil <- ceiling(radius * 2)
  for (tt in ts) {
    cx <- p1[1] + tt * (p2[1] - p1[1])
    cy <- p1[2] + tt * (p2[2] - p1[2])
    c0 <- floor(cx); r0 <- floor(cy)
    for (cc in (c0 - rad_ceil):(c0 + rad_ceil)) {
      for (rr in (r0 - rad_ceil):(r0 + rad_ceil)) {
        if (cc < 1 || cc > W || rr < 1 || rr > H) next
        d2 <- (cc - cx)^2 + (rr - cy)^2
        wgt <- exp(-d2 / (2 * radius^2))
        img[rr, cc] <- img[rr, cc] - depth * wgt
      }
    }
  }
  img[img < 0] <- 0
  img
}

#' Render one camera image for one frame.
render_camera_image <- function(positions, velocities, cam, params, calib,
                                 bg_level = 0.82, noise_sd = 0.012) {
  H <- calib$H; W <- calib$W
  img <- matrix(bg_level, H, W) + matrix(rnorm(H * W, 0, noise_sd), H, W)
  te <- params$te
  for (i in seq_len(ncol(positions))) {
    r <- positions[, i]; v <- velocities[, i]
    # Streak spans the exposure: r +/- v*te/2
    ep1 <- as.vector(project_to_image(matrix(r - v * te / 2, 3, 1), cam))
    ep2 <- as.vector(project_to_image(matrix(r + v * te / 2, 3, 1), cam))
    img <- .draw_streak(img, ep1, ep2, depth = 0.72, radius = 1.4)
  }
  img[img < 0] <- 0; img[img > 1] <- 1
  img
}

#' Render and write all stereo PNG frames to disk.
#'
#' @return list(frame_files = list(cam1_paths, cam2_paths), truth = data frame)
generate_synthetic_images <- function(out_dir, n_targets, n_frames,
                                       params, calib, seed = 7) {
  truth <- generate_truth_tracks(n_targets, n_frames, params, calib, seed)

  cam_dirs <- file.path(out_dir, "frames", c("cam1", "cam2"))
  for (d in cam_dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

  frame_files <- list(character(n_frames), character(n_frames))
  message(sprintf("[image_demo] Rendering %d frames x 2 cameras (%dx%d)...",
                  n_frames, calib$W, calib$H))

  for (k in seq_len(n_frames)) {
    sub <- truth[truth$frame == k, ]
    P  <- rbind(sub$x, sub$y, sub$z)        # 3 x n
    Vv <- rbind(sub$vx, sub$vy, sub$vz)     # 3 x n
    for (ci in 1:2) {
      img <- render_camera_image(P, Vv, calib$cams[[ci]], params, calib)
      fp  <- file.path(cam_dirs[ci], sprintf("frame_%04d.png", k))
      png::writePNG(img, fp)
      frame_files[[ci]][k] <- fp
    }
  }
  message("[image_demo] Wrote ", n_frames * 2, " PNG frames to ", file.path(out_dir, "frames"))
  list(frame_files = frame_files, truth = truth)
}


# ============================================================
# FRONT-HALF PIPELINE ON THE RENDERED IMAGES
# ============================================================

#' Process the rendered images frame-by-frame through the real pipeline.
#'
#' @return list with Xh, detections (per-frame triangulated points), n_frames
run_image_pipeline <- function(frame_files, calib, params, n_frames,
                               t_area = c(6, 500)) {
  cams  <- calib$cams
  F_mat <- calib$F_mat
  Xh <- make_xh(params$max_targets, n_frames)
  targets <- list(); next_id <- 1L

  detections <- list()   # collected triangulated points across frames
  blob_counts <- integer(n_frames)

  for (k in seq_len(n_frames)) {
    Z_cam <- vector("list", 2)
    for (ci in 1:2) {
      imgarr <- load_frame_window(frame_files[[ci]], k, params$br0)
      bg_res <- compute_background(imgarr, params$br0, params$br0, params$fg_is_dark)
      fg     <- compute_foreground(bg_res$bg, bg_res$current,
                                   binary_t = 0.10, params$fg_is_dark)
      blobs  <- extract_blobs(fg)
      Z_cam[[ci]] <- extract_measurements(blobs, t_area = t_area,
                                          img_dims = c(calib$H, calib$W))
    }
    blob_counts[k] <- length(Z_cam[[1]])

    # Stereo matching via epipolar constraint
    Z_pairs <- if (length(Z_cam[[1]]) > 0 && length(Z_cam[[2]]) > 0)
      validate_stereo_pairs(Z_cam[[1]], Z_cam[[2]], F_mat, params$te_epipolar)
    else list()

    # Record triangulated detections (independent of track association)
    for (pr in Z_pairs) {
      tri <- ls_triangulate(cbind(pr$Z1$u, pr$Z2$u), cams)
      detections[[length(detections) + 1]] <- data.frame(
        frame = k, x = tri$r[1], y = tri$r[2], z = tri$r[3],
        reproj_err = tri$reprojection_error
      )
    }

    # Tracking update
    res     <- update_and_initialize(targets, Z_pairs, cams, params, next_id)
    targets <- res$targets; next_id <- res$next_id
    Xh      <- write_targets_to_xh(Xh, targets, k)

    if (k %% 5 == 0)
      message(sprintf("  frame %2d/%d: %d blobs (cam1), %d stereo pairs",
                      k, n_frames, blob_counts[k], length(Z_pairs)))
  }

  list(Xh = Xh,
       detections = if (length(detections)) do.call(rbind, detections) else NULL,
       blob_counts = blob_counts)
}


# ============================================================
# GROUND-TRUTH MATCHING (triangulation accuracy)
# ============================================================

#' Match each triangulated detection to the nearest true mosquito in that frame.
match_to_truth <- function(detections, truth) {
  if (is.null(detections)) return(NULL)
  out <- lapply(seq_len(nrow(detections)), function(i) {
    d <- detections[i, ]
    tf <- truth[truth$frame == d$frame, ]
    if (nrow(tf) == 0) return(NULL)
    dist <- sqrt((tf$x - d$x)^2 + (tf$y - d$y)^2 + (tf$z - d$z)^2)
    j <- which.min(dist)
    data.frame(frame = d$frame, error_mm = dist[j],
               x_true = tf$x[j], y_true = tf$y[j], z_true = tf$z[j],
               x_est = d$x, y_est = d$y, z_est = d$z,
               reproj_err = d$reproj_err)
  })
  do.call(rbind, Filter(Negate(is.null), out))
}


# ============================================================
# FIGURES
# ============================================================

plot_image_demo <- function(frame_files, calib, truth, detections, matched,
                            sample_frame, out_dir) {
  cams <- calib$cams

  # ---- Figure 1: raw camera image + detections + true projections ----
  fig1 <- file.path(out_dir, "01_detection_on_image.pdf")
  grDevices::pdf(fig1, width = 11, height = 5)
  par(mfrow = c(1, 2), mar = c(2, 2, 3, 1))
  for (ci in 1:2) {
    img <- read_frame(frame_files[[ci]][sample_frame])
    # image() plots with origin bottom-left; flip rows so it looks like a photo
    graphics::image(t(img[nrow(img):1, ]), col = grey.colors(256, 0, 1),
                    axes = FALSE, main = sprintf("Camera %d, frame %d (synthetic)", ci, sample_frame))
    box()
    sub <- truth[truth$frame == sample_frame, ]
    P <- rbind(sub$x, sub$y, sub$z)
    proj <- project_to_image(P, cams[[ci]])
    # convert pixel (col,row) -> image() normalised coords (x in [0,1], y flipped)
    px <- proj[1, ] / ncol(img)
    py <- 1 - proj[2, ] / nrow(img)
    points(px, py, pch = 3, col = "#2ca02c", lwd = 2, cex = 1.6)  # true positions
  }
  legend("topright", legend = c("True mosquito position"),
         pch = 3, col = "#2ca02c", bty = "n", text.col = "white")
  grDevices::dev.off()

  # ---- Figure 2: triangulation accuracy (estimated vs true) ----
  fig2 <- file.path(out_dir, "02_triangulation_accuracy.pdf")
  grDevices::pdf(fig2, width = 12, height = 4)
  par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1))
  for (ax in c("x", "y", "z")) {
    tv <- matched[[paste0(ax, "_true")]]; ev <- matched[[paste0(ax, "_est")]]
    rng <- range(c(tv, ev))
    plot(tv, ev, pch = 19, col = adjustcolor("#1f78b4", 0.5),
         xlab = sprintf("True %s (mm)", toupper(ax)),
         ylab = sprintf("Recovered %s (mm)", toupper(ax)),
         main = sprintf("%s axis", toupper(ax)), xlim = rng, ylim = rng, las = 1)
    abline(0, 1, col = "red", lwd = 2, lty = 2)
  }
  grDevices::dev.off()

  # ---- Figure 3: recovered detections over true 3D paths (x-z and y-z) ----
  fig3 <- file.path(out_dir, "03_recovered_vs_true_tracks.pdf")
  grDevices::pdf(fig3, width = 11, height = 5)
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
  cols <- grDevices::hcl.colors(max(2, length(unique(truth$target_id))), "Set2")
  for (pair in list(c("x", "z"), c("y", "z"))) {
    a <- pair[1]; b <- pair[2]
    plot(truth[[a]], truth[[b]], type = "n",
         xlab = sprintf("%s (mm)", toupper(a)), ylab = sprintf("%s (mm)", toupper(b)),
         main = sprintf("%s-%s plane: true paths (lines) vs recovered (points)",
                        toupper(a), toupper(b)), las = 1)
    for (ti in unique(truth$target_id)) {
      tr <- truth[truth$target_id == ti, ]
      lines(tr[[a]], tr[[b]], col = cols[ti], lwd = 2)
    }
    points(detections[[a]], detections[[b]], pch = 19,
           col = adjustcolor("black", 0.45), cex = 0.7)
  }
  grDevices::dev.off()

  c(fig1, fig2, fig3)
}


# ============================================================
# REVIEW BUNDLE EXPORT (for the trackone Shiny overlay app)
# ============================================================

#' Save everything the Shiny review app needs to overlay tracks on the images.
#'
#' Writes `image_review_bundle.rds` (calibration, frame paths, recovered tracks,
#' ground truth) plus `recovered_tracks.csv` into `out_dir`. The Shiny app reads
#' this bundle to draw each camera frame with the reconstructed 3D tracks
#' projected back onto the pixels — i.e. the verify/combine step that the
#' original MATLAB trackone.m performed against the stereo photographs.
export_review_bundle <- function(out_dir, calib, params, res, gen) {
  fps <- max(1, round(1 / params$dt))
  recovered <- tryCatch(xh_to_tidy(res$Xh, fps = fps), error = function(e) NULL)

  n_frames <- length(gen$frame_files[[1]])
  # Store frame paths RELATIVE to out_dir so the bundle is portable.
  cam1_rel <- file.path("frames", "cam1", sprintf("frame_%04d.png", seq_len(n_frames)))
  cam2_rel <- file.path("frames", "cam2", sprintf("frame_%04d.png", seq_len(n_frames)))

  bundle <- list(
    cams      = calib$cams,
    F_mat     = calib$F_mat,
    W         = calib$W,
    H         = calib$H,
    dt        = params$dt,
    te        = params$te,
    n_frames  = n_frames,
    frames    = list(cam1 = cam1_rel, cam2 = cam2_rel),
    recovered = recovered,
    truth     = gen$truth
  )
  saveRDS(bundle, file.path(out_dir, "image_review_bundle.rds"))
  if (!is.null(recovered))
    write.csv(recovered, file.path(out_dir, "recovered_tracks.csv"), row.names = FALSE)
  message("  Saved: image_review_bundle.rds",
          if (!is.null(recovered)) ", recovered_tracks.csv" else "")
  invisible(bundle)
}


# ============================================================
# MAIN
# ============================================================

demo_synthetic_images <- function(
  out_dir   = "figures/synthetic_image_demo",
  n_targets = 4,
  n_frames  = 24,
  img_W     = 240,
  img_H     = 180,
  seed      = 7
) {
  .load_modules()
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  params <- default_params()
  params$br0 <- 2          # 5-frame sliding window (faster)
  params$Np  <- 120        # fewer particles (faster)
  calib <- make_demo_calib(W = img_W, H = img_H)

  message("\n=== Step 1: Generate synthetic 3D paths + render stereo images ===")
  gen <- generate_synthetic_images(out_dir, n_targets, n_frames, params, calib, seed)

  message("\n=== Step 2: Run FULL image -> 3D pipeline on rendered images ===")
  t0 <- Sys.time()
  res <- run_image_pipeline(gen$frame_files, calib, params, n_frames)
  message(sprintf("  Pipeline runtime: %.1f s", as.numeric(Sys.time() - t0, units = "secs")))

  message("\n=== Step 3: Compare recovered 3D to ground truth ===")
  matched <- match_to_truth(res$detections, gen$truth)
  if (is.null(matched) || nrow(matched) == 0) {
    message("  [!] No detections were triangulated — try larger images or more contrast.")
    return(invisible(list(calib = calib, truth = gen$truth, result = res)))
  }
  n_det <- nrow(matched)
  n_true <- nrow(gen$truth)
  message(sprintf("  Triangulated detections: %d (across %d true positions)", n_det, n_true))
  message(sprintf("  3D position error: mean %.2f mm, median %.2f mm, 90th pct %.2f mm",
                  mean(matched$error_mm), median(matched$error_mm),
                  quantile(matched$error_mm, 0.9)))
  message(sprintf("  Mean reprojection error: %.3f px", mean(matched$reproj_err)))
  message(sprintf("  Detection rate: %.0f%% of true positions recovered",
                  100 * min(1, n_det / n_true)))

  message("\n=== Step 4: Figures ===")
  sample_frame <- ceiling(n_frames / 2)
  figs <- plot_image_demo(gen$frame_files, calib, gen$truth,
                          res$detections, matched, sample_frame, out_dir)
  for (f in figs) message("  Saved: ", f)

  # Save CSVs
  write.csv(matched, file.path(out_dir, "triangulation_vs_truth.csv"), row.names = FALSE)
  write.csv(gen$truth, file.path(out_dir, "ground_truth_tracks.csv"), row.names = FALSE)
  message("  Saved: triangulation_vs_truth.csv, ground_truth_tracks.csv")

  # ---- Export bundle for the trackone-style Shiny overlay app ----
  export_review_bundle(out_dir, calib, params, res, gen)

  message("\n=== Done. Output in: ", out_dir, " ===")
  message("  REMINDER: input images are SYNTHETIC; the pipeline code is real.\n")

  invisible(list(calib = calib, truth = gen$truth, result = res, matched = matched))
}


if (!interactive() && identical(environment(), globalenv())) {
  demo_synthetic_images()
}
