# postprocess.R
# Postprocessing: Kalman smoothing, trajectory tidying, summary metrics.
#
# Mirrors MATLAB source files:
#   tracker/filter_and_smooth.m  (Kalman smoother pass)
#   tracker/nnf_kf.m             (nearest-neighbour Kalman filter update)
#   tracker/kalmanPredict.m      (KF prediction)
#   tracker/kalmanUpdate.m       (KF update)
#   tracker/sma.m                (moving average)
#   tracker/mydiff.m             (finite differences)
#   analysis/update_data_with_calibration.m  (re-triangulate)
#   analysis/compareAutoWithManual.m (summary statistics)
#
# STATUS: translated
# See audit/MATLAB_TO_R_CROSSWALK.csv for details.

# ============================================================
# KALMAN FILTER PRIMITIVES (mirrors kalmanPredict.m, kalmanUpdate.m)
# ============================================================

#' Kalman filter predict step
#'
#' @param x   n x 1 state mean
#' @param P   n x n state covariance
#' @param F   n x n state transition matrix
#' @param Q   n x n process noise covariance
#' @return list with x_pred, P_pred
kf_predict <- function(x, P, F, Q) {
  x_pred <- F %*% x
  P_pred <- F %*% P %*% t(F) + Q
  list(x = as.vector(x_pred), P = P_pred)
}

#' Kalman filter update step
#'
#' @param x_pred  n x 1 predicted state
#' @param P_pred  n x n predicted covariance
#' @param z       m x 1 measurement
#' @param H       m x n measurement matrix
#' @param R       m x m measurement noise covariance
#' @return list with x_upd, P_upd, K (gain), inov (innovation)
kf_update_step <- function(x_pred, P_pred, z, H, R) {
  S    <- H %*% P_pred %*% t(H) + R
  K    <- P_pred %*% t(H) %*% solve(S)
  inov <- z - H %*% x_pred
  x_upd <- as.vector(x_pred + K %*% inov)
  P_upd <- (diag(nrow(P_pred)) - K %*% H) %*% P_pred
  list(x = x_upd, P = P_upd, K = K, inov = as.vector(inov))
}


# ============================================================
# KALMAN SMOOTHER (mirrors filter_and_smooth.m, nnf_kf.m)
# ============================================================

#' RTS (Rauch-Tung-Striebel) Kalman smoother for a single track
#'
#' Runs a forward Kalman filter pass followed by a backward smoothing pass.
#' This is the standard smoother used in the paper (section 3.3).
#' Mirrors the intent of filter_and_smooth.m / nnf_kf.m.
#'
#' @param track_data  data frame (subset of tidy Xh) for one target, columns:
#'                    frame, x, y, z, vx, vy, vz  (NA for untracked frames)
#' @param dt  time step (s)
#' @param sigma_w  process noise (mm^2/s^4)
#' @param sigma_r  measurement noise std (mm, default 20 mm)
#' @return data frame with smoothed x, y, z, vx, vy, vz columns
kalman_smooth <- function(track_data, dt = 1/25, sigma_w = 100e6, sigma_r = 20) {
  track <- track_data[order(track_data$frame), ]
  n_frames <- nrow(track)
  if (n_frames < 2) return(track)

  F <- cv_motion_matrix(dt)
  Q <- cv_process_noise(dt, sigma_w)
  H <- cbind(diag(3), matrix(0, 3, 3))   # observe position only
  R <- diag(sigma_r^2, 3)

  # --- Forward pass ---
  xf <- matrix(NA_real_, 6, n_frames)
  Pf <- array(NA_real_,  dim = c(6, 6, n_frames))
  xp <- matrix(NA_real_, 6, n_frames)   # predicted
  Pp <- array(NA_real_,  dim = c(6, 6, n_frames))

  # Initialise with first valid frame
  first_valid <- which(!is.na(track$x))[1]
  if (is.na(first_valid)) return(track)  # no valid data

  x0 <- c(track$x[first_valid], track$y[first_valid], track$z[first_valid],
          if (!is.na(track$vx[first_valid])) c(track$vx[first_valid], track$vy[first_valid], track$vz[first_valid])
          else c(0, 0, 0))
  P0 <- diag(c(50^2, 50^2, 50^2, 500^2, 500^2, 500^2))

  x_cur <- x0; P_cur <- P0

  for (k in seq_len(n_frames)) {
    # Predict
    pred    <- kf_predict(x_cur, P_cur, F, Q)
    xp[, k] <- pred$x
    Pp[,,k] <- pred$P

    # Update if we have a valid measurement
    if (!is.na(track$x[k])) {
      z <- c(track$x[k], track$y[k], track$z[k])
      upd <- kf_update_step(pred$x, pred$P, z, H, R)
      x_cur <- upd$x
      P_cur <- upd$P
    } else {
      x_cur <- pred$x
      P_cur <- pred$P
    }
    xf[, k] <- x_cur
    Pf[,,k] <- P_cur
  }

  # --- Backward RTS pass ---
  xs <- xf  # smoothed estimates
  for (k in (n_frames - 1):1) {
    G   <- Pf[,,k] %*% t(F) %*% solve(Pp[,,k+1])
    xs[, k] <- xf[, k] + G %*% (xs[, k+1] - xp[, k+1])
  }

  # Write back
  track$x  <- xs[1, ]; track$y  <- xs[2, ]; track$z  <- xs[3, ]
  track$vx <- xs[4, ]; track$vy <- xs[5, ]; track$vz <- xs[6, ]
  track
}

#' Apply Kalman smoother to all tracks in a tidy trajectory data frame
#'
#' @param traj_df  tidy data frame with columns target_id, frame, x, y, z, vx, vy, vz
#' @param dt  time step (s)
#' @param sigma_w  process noise
#' @param sigma_r  measurement noise
#' @return smoothed data frame
smooth_all_tracks <- function(traj_df, dt = 1/25, sigma_w = 100e6, sigma_r = 20) {
  target_ids <- unique(traj_df$target_id)
  smoothed <- do.call(rbind, lapply(target_ids, function(tid) {
    sub <- traj_df[traj_df$target_id == tid, ]
    if (nrow(sub) < 2) return(sub)
    tryCatch(
      kalman_smooth(sub, dt = dt, sigma_w = sigma_w, sigma_r = sigma_r),
      error = function(e) {
        warning(sprintf("[postprocess.R] Smoother failed for target %d: %s", tid, e$message))
        sub
      }
    )
  }))
  smoothed
}


# ============================================================
# SUMMARY METRICS (mirrors compareAutoWithManual.m, paper Figures 9-13)
# ============================================================

#' Compute speed for each position in a trajectory
#'
#' @param traj_df  tidy trajectory data frame
#' @return same data frame with speed_mms (mm/s), speed_ms (m/s), horizontal_speed_ms, vertical_speed_ms
add_speed <- function(traj_df) {
  traj_df$speed_mms       <- sqrt(traj_df$vx^2 + traj_df$vy^2 + traj_df$vz^2)
  traj_df$speed_ms        <- traj_df$speed_mms / 1000
  traj_df$horizontal_speed_ms <- sqrt(traj_df$vx^2 + traj_df$vy^2) / 1000
  traj_df$vertical_speed_ms   <- abs(traj_df$vz) / 1000
  traj_df
}

#' Compute swarm centroid and 3-sigma bounds per frame
#'
#' Mirrors the centroid/std calculation in compareAutoWithManual.m and
#' paper Figure 9 (3-sigma bounds for all mosquitoes).
#'
#' @param traj_df  tidy trajectory data frame
#' @return data frame with frame, time_s, mean_x/y/z (mm), std_x/y/z (mm), n_mosquitoes
compute_swarm_stats <- function(traj_df) {
  frames <- sort(unique(traj_df$frame))
  result <- do.call(rbind, lapply(frames, function(f) {
    sub <- traj_df[traj_df$frame == f & !is.na(traj_df$x), ]
    n <- nrow(sub)
    if (n == 0) return(NULL)
    data.frame(
      frame        = f,
      time_s       = sub$time_s[1],
      n_mosquitoes = n,
      mean_x = mean(sub$x), mean_y = mean(sub$y), mean_z = mean(sub$z),
      std_x  = if(n>1) sd(sub$x) else 0,
      std_y  = if(n>1) sd(sub$y) else 0,
      std_z  = if(n>1) sd(sub$z) else 0
    )
  }))
  result
}

#' Compute separation distance between two targets over time
#'
#' Mirrors paper Figure 13 (mating event separation distance).
#'
#' @param traj_df  tidy trajectory data frame
#' @param id1, id2  target IDs to compare
#' @return data frame with frame, time_s, separation_mm, separation_m
compute_separation <- function(traj_df, id1, id2) {
  t1 <- traj_df[traj_df$target_id == id1, c("frame","time_s","x","y","z")]
  t2 <- traj_df[traj_df$target_id == id2, c("frame","time_s","x","y","z")]
  merged <- merge(t1, t2, by = c("frame","time_s"), suffixes = c("_1","_2"))
  merged$separation_mm <- sqrt((merged$x_1 - merged$x_2)^2 +
                                (merged$y_1 - merged$y_2)^2 +
                                (merged$z_1 - merged$z_2)^2)
  merged$separation_m  <- merged$separation_mm / 1000
  merged[, c("frame","time_s","separation_mm","separation_m")]
}

#' Compute total path length for each target
#'
#' @param traj_df  tidy trajectory data frame
#' @return data frame with target_id and path_length_mm
compute_path_lengths <- function(traj_df) {
  target_ids <- unique(traj_df$target_id)
  do.call(rbind, lapply(target_ids, function(tid) {
    sub <- traj_df[traj_df$target_id == tid & !is.na(traj_df$x), ]
    sub <- sub[order(sub$frame), ]
    if (nrow(sub) < 2) return(data.frame(target_id = tid, path_length_mm = 0))
    dx <- diff(sub$x); dy <- diff(sub$y); dz <- diff(sub$z)
    data.frame(target_id = tid,
               path_length_mm = sum(sqrt(dx^2 + dy^2 + dz^2)))
  }))
}

#' Re-triangulate tracks with updated calibration
#'
#' Mirrors analysis/update_data_with_calibration.m:
#' Takes existing Xh (positions in old camera frame), projects to pixels
#' using old cams, then triangulates using new cams.
#'
#' @param Xh_old  state matrix with old calibration positions
#' @param cams_old  old camera structs
#' @param cams_new  new camera structs
#' @return Xh_new with re-triangulated positions
update_with_calibration <- function(Xh_old, cams_old, cams_new) {
  n_targets <- nrow(Xh_old) %/% 6
  n_frames  <- ncol(Xh_old)
  Xh_new    <- Xh_old * 0

  for (t in seq_len(n_targets)) {
    rows <- (t-1)*6 + 1:3
    pos_old <- Xh_old[rows, , drop = FALSE]  # 3 x Nframes

    for (k in seq_len(n_frames)) {
      r <- pos_old[, k]
      if (all(r == 0)) next

      pix1 <- as.vector(project_to_image(matrix(r, 3, 1), cams_old[[1]]))
      pix2 <- as.vector(project_to_image(matrix(r, 3, 1), cams_old[[2]]))

      tri <- tryCatch(
        ls_triangulate(cbind(pix1, pix2), cams_new),
        error = function(e) list(r = r)
      )
      Xh_new[rows, k] <- tri$r
    }
  }

  # Copy velocities unchanged
  for (t in seq_len(n_targets)) {
    vel_rows <- (t-1)*6 + 4:6
    Xh_new[vel_rows, ] <- Xh_old[vel_rows, ]
  }
  Xh_new
}


# ============================================================
# HELPER: motion model functions (used in smoother, re-exported)
# ============================================================

#' Discrete constant-velocity state transition matrix (re-exported from tracking.R)
#' Included here so postprocess.R can run standalone.
cv_motion_matrix <- function(dt) {
  F <- diag(6)
  F[1:3, 4:6] <- diag(dt, 3)
  F
}

cv_process_noise <- function(dt, sigma_w = 100e6) {
  Q <- matrix(0, 6, 6)
  Q[1:3, 1:3] <- diag(sigma_w * dt^3 / 3, 3)
  Q[1:3, 4:6] <- diag(sigma_w * dt^2 / 2, 3)
  Q[4:6, 1:3] <- diag(sigma_w * dt^2 / 2, 3)
  Q[4:6, 4:6] <- diag(sigma_w * dt, 3)
  Q
}
