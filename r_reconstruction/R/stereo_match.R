# stereo_match.R
# Stereo geometry: calibration, epipolar constraint, triangulation, projection, gating.
#
# Mirrors MATLAB source files:
#   tracker/get_F_for_stereo.m   (fundamental matrix)
#   tracker/lsTriangulate.m      (DLT least-squares triangulation)
#   tracker/w2cam_nd.m           (world -> image projection, no distortion)
#   tracker/w2cam.m              (world -> image projection for particles)
#   tracker/cam2world.m          (image -> 3D ray)
#   tracker/tra2b.m              (homogeneous transform)
#   tracker/normv.m              (vector normalisation)
#   tracker/vclz.m               (Mahalanobis gating check)
#   tracker/get_undistorted_points.m  (radial distortion correction)
#
# STATUS: translated
# See audit/MATLAB_TO_R_CROSSWALK.csv for details.

# ============================================================
# HOMOGENEOUS TRANSFORMS (mirrors tra2b.m, trpa2b.m)
# ============================================================

#' Transform 3D points by a 4x4 homogeneous matrix
#'
#' Mirrors tra2b.m: res = T * [pts; 1]  (returns 3 x N)
#'
#' @param pts  3 x N matrix of [x; y; z] world coordinates
#' @param T    4x4 homogeneous transform
#' @return 3 x N transformed points
transform_points <- function(pts, T) {
  if (is.vector(pts)) pts <- matrix(pts, nrow = 3)
  N <- ncol(pts)
  homo <- rbind(pts, rep(1, N))  # 4 x N
  res  <- T %*% homo             # 4 x N
  res[1:3, , drop = FALSE]
}


# ============================================================
# PROJECTION (mirrors w2cam_nd.m, w2cam.m)
# ============================================================

#' Project 3D world coordinates to image pixels
#'
#' Mirrors w2cam_nd.m (no lens distortion applied):
#'   cr = trm * [r; 1]  (camera frame)
#'   xn = cr[1]/cr[3], yn = cr[2]/cr[3]  (normalised)
#'   pixel = km * [xn; yn; 1]
#'
#' @param pts_world  3 x N matrix of world coordinates (mm)
#' @param cam  camera struct with km (3x3) and trm (4x4)
#' @param apply_distortion  logical (default FALSE; kc1/kc2 correction)
#' @return 2 x N matrix of pixel coordinates [col; row] = [u; v]
project_to_image <- function(pts_world, cam, apply_distortion = FALSE) {
  if (is.vector(pts_world)) pts_world <- matrix(pts_world, nrow = 3)

  # Transform to camera frame
  cr <- transform_points(pts_world, cam$trm)  # 3 x N

  # Perspective division
  z <- cr[3, ]
  xn <- cr[1, ] / z
  yn <- cr[2, ] / z

  # Apply radial distortion (partial: kc1, kc2 only)
  if (apply_distortion) {
    rn2 <- xn^2 + yn^2
    factor <- 1 + cam$kc1 * rn2 + cam$kc2 * rn2^2
    xn <- xn * factor
    yn <- yn * factor
  }

  # Apply intrinsic matrix
  pix <- cam$km %*% rbind(xn, yn, rep(1, length(xn)))
  pix[1:2, , drop = FALSE]
}

# Alias used in particle filter code
project_particles <- project_to_image

#' Project image pixel to 3D ray direction (mirrors cam2world.m)
#'
#' @param pixel  2 x 1 pixel coordinates [col; row]
#' @param cam  camera struct
#' @return 3 x 1 unit ray direction in world frame
cam_to_world_ray <- function(pixel, cam) {
  pixel_h <- c(pixel, 1)
  # Normalised image coordinates
  xn <- solve(cam$km) %*% pixel_h
  xn <- xn / xn[3]
  # Rotate from camera to world frame
  R_wc <- t(cam$trm[1:3, 1:3])  # world-from-camera rotation
  ray <- R_wc %*% xn[1:3]
  ray / sqrt(sum(ray^2))
}

#' Undistort image points (mirrors get_undistorted_points.m)
#'
#' Applies inverse radial distortion model (kc1, kc2).
#' STATUS: partial — uses kc1/kc2 only; higher-order Bouguet terms omitted.
#'
#' @param pixels  2 x N distorted pixel coordinates
#' @param cam  camera struct with km, kc1, kc2
#' @param max_iter  Newton iterations
#' @return 2 x N undistorted pixel coordinates
undistort_points <- function(pixels, cam, max_iter = 5) {
  if (is.vector(pixels)) pixels <- matrix(pixels, nrow = 2)
  pp  <- cam$km[1:2, 3]        # principal point [cx; cy]
  f   <- c(cam$km[1,1], cam$km[2,2])  # focal lengths

  xn <- (pixels[1, ] - pp[1]) / f[1]
  yn <- (pixels[2, ] - pp[2]) / f[2]

  # Iterative undistortion
  xu <- xn; yu <- yn
  for (i in seq_len(max_iter)) {
    r2 <- xu^2 + yu^2
    factor <- 1 + cam$kc1 * r2 + cam$kc2 * r2^2
    xu <- xn / factor
    yu <- yn / factor
  }

  rbind(xu * f[1] + pp[1],
        yu * f[2] + pp[2])
}


# ============================================================
# FUNDAMENTAL MATRIX (mirrors get_F_for_stereo.m)
# ============================================================

#' Compute the fundamental matrix for a stereo camera pair
#'
#' Mirrors get_F_for_stereo.m exactly:
#'   t = cam2.trm[1:3, 4]  (translation of cam2 relative to cam1)
#'   tx = skew-symmetric cross-product matrix of t
#'   E = tx * R2  (essential matrix)
#'   F = K2^{-T} * E * K1^{-1}
#'
#' Camera 1 is assumed to be at the origin (identity trm).
#'
#' @param cam1  camera struct for left camera
#' @param cam2  camera struct for right camera
#' @return 3x3 fundamental matrix F
compute_fundamental_matrix <- function(cam1, cam2) {
  t <- cam2$trm[1:3, 4]
  R2 <- cam2$trm[1:3, 1:3]

  tx <- matrix(c(0,    -t[3],  t[2],
                  t[3],  0,    -t[1],
                 -t[2],  t[1],  0), nrow = 3, byrow = TRUE)

  E <- tx %*% R2
  F <- solve(t(cam2$km)) %*% E %*% solve(cam1$km)
  F
}


# ============================================================
# EPIPOLAR VALIDATION (mirrors paper eq 3.2)
# ============================================================

#' Test epipolar constraint for a stereo measurement pair
#'
#' Checks: |ũ2^T * F * ũ1| < te
#' where ũ = [u; 1] is homogeneous pixel coordinate.
#' Ref: paper eq 3.2; Hartley & Zisserman 2004.
#'
#' @param u1  2-element pixel midpoint from camera 1 [col, row]
#' @param u2  2-element pixel midpoint from camera 2
#' @param F   3x3 fundamental matrix
#' @param te  epipolar threshold (default 0.5, paper Table 2)
#' @return logical: TRUE if pair satisfies epipolar constraint
epipolar_check <- function(u1, u2, F, te = 0.5) {
  u1h <- c(u1, 1)
  u2h <- c(u2, 1)
  abs(t(u2h) %*% F %*% u1h) < te
}

#' Validate all stereo measurement pairs using the epipolar constraint
#'
#' For each pair (Z1 from cam1, Z2 from cam2), tests the constraint
#' on their midpoints. Returns indices of valid pairs.
#'
#' @param Z_cam1  list of Z measurement structs from camera 1
#' @param Z_cam2  list of Z measurement structs from camera 2
#' @param F       3x3 fundamental matrix
#' @param te      epipolar threshold
#' @return list of valid pairs, each element: list(i1=, i2=, Z1=, Z2=)
validate_stereo_pairs <- function(Z_cam1, Z_cam2, F, te = 0.5) {
  valid_pairs <- list()
  for (i1 in seq_along(Z_cam1)) {
    for (i2 in seq_along(Z_cam2)) {
      if (epipolar_check(Z_cam1[[i1]]$u, Z_cam2[[i2]]$u, F, te)) {
        valid_pairs[[length(valid_pairs) + 1]] <- list(
          i1 = i1, i2 = i2,
          Z1 = Z_cam1[[i1]], Z2 = Z_cam2[[i2]]
        )
      }
    }
  }
  valid_pairs
}


# ============================================================
# TRIANGULATION (mirrors lsTriangulate.m)
# ============================================================

#' Least-squares DLT triangulation of a stereo measurement pair
#'
#' Mirrors lsTriangulate.m exactly:
#'   For each camera c, the constraint is:
#'     [K_c[1,:] - u_c[1]*K_c[3,:]] * trm_c * r_h = 0
#'     [K_c[2,:] - u_c[2]*K_c[3,:]] * trm_c * r_h = 0
#'   Stacked as A * r_h = 0, solved by (A'A)^{-1} A' b.
#'
#' Ref: Hartley & Sturm (1997), paper ref [38].
#'
#' @param pixels  2 x nc matrix of pixel coordinates, one column per camera
#' @param cams  list of camera structs
#' @return list with r [3x1] 3D position (mm) and reprojection error (scalar)
ls_triangulate <- function(pixels, cams) {
  nc <- length(cams)
  if (ncol(pixels) != nc) stop("[stereo_match.R] pixels columns must match number of cameras")

  A <- matrix(0, nrow = 2 * nc, ncol = 4)
  for (c in seq_len(nc)) {
    km  <- cams[[c]]$km
    trm <- cams[[c]]$trm[1:3, 1:4]  # 3 x 4

    u <- pixels[1, c]
    v <- pixels[2, c]

    A[2*(c-1)+1, ] <- (km[1,] - u * km[3,]) %*% trm
    A[2*(c-1)+2, ] <- (km[2,] - v * km[3,]) %*% trm
  }

  b <- -A[, 4]
  Asub <- A[, 1:3]

  r <- tryCatch(
    solve(t(Asub) %*% Asub) %*% t(Asub) %*% b,
    error = function(e) {
      warning("[stereo_match.R] Triangulation singular: using pseudoinverse")
      MASS_pinv(t(Asub) %*% Asub) %*% t(Asub) %*% b
    }
  )

  # Reprojection error
  r_h <- c(r, 1)
  err <- 0
  for (c in seq_len(nc)) {
    reproj <- cams[[c]]$km %*% cams[[c]]$trm[1:3, 1:4] %*% r_h
    reproj_px <- reproj[1:2] / reproj[3]
    err <- err + sqrt(sum((pixels[, c] - reproj_px)^2))
  }

  list(r = as.vector(r), reprojection_error = err / nc)
}

# Minimal pseudoinverse via SVD (avoids MASS dependency)
MASS_pinv <- function(A, tol = 1e-10) {
  svd_res <- svd(A)
  d <- svd_res$d
  d_inv <- ifelse(d > tol, 1/d, 0)
  svd_res$v %*% diag(d_inv, length(d)) %*% t(svd_res$u)
}


# ============================================================
# GATING (mirrors vclz.m)
# ============================================================

#' Mahalanobis gating check between a measurement and a cluster
#'
#' Mirrors vclz.m: computes chi^2 = (z_hat - z)' S^{-1} (z_hat - z)
#' for each target in the cluster, where z_hat is the projected predicted
#' position. Returns the minimum chi^2 across all targets.
#'
#' @param z_meas  2 x 1 measured pixel midpoint
#' @param particle_cloud  3 x Np particle position matrix for the target
#' @param cam  camera struct
#' @param t_gate  chi^2 gating threshold (default 16, paper Table 2)
#' @return chi^2 value (scalar); check: value < t_gate means inside gate
gate_check_single <- function(z_meas, particle_cloud, cam, t_gate = 16) {
  if (is.vector(particle_cloud)) particle_cloud <- matrix(particle_cloud, nrow = 3)

  z_hat_all <- project_to_image(particle_cloud, cam)  # 2 x Np
  z_hat_mean <- rowMeans(z_hat_all)                   # 2 x 1
  S <- cov(t(z_hat_all))                              # 2 x 2
  if (det(S) < 1e-5) {
    S <- S + diag(1e-3, 2)  # regularise singular covariance
  }
  inov <- z_meas - z_hat_mean
  as.numeric(t(inov) %*% solve(S) %*% inov)
}

#' Gating check using predicted state mean (fast, no particles needed)
#'
#' Used when particles are not yet initialised (KF-only path).
#'
#' @param z_meas    2 x 1 pixel measurement
#' @param r_pred    3 x 1 predicted 3D position
#' @param P_pred    3 x 3 position covariance
#' @param cam       camera struct
#' @param t_gate    chi^2 threshold
#' @return chi^2 value
gate_check_kf <- function(z_meas, r_pred, P_pred, cam, t_gate = 16) {
  # Innovation covariance: S = H P H' where H = d(project)/d(r)
  H <- .compute_H(r_pred, cam)   # 2 x 3 measurement Jacobian
  S <- H %*% P_pred %*% t(H) + diag(1, 2)  # add small measurement noise
  z_hat <- project_to_image(r_pred, cam)
  inov  <- z_meas - as.vector(z_hat)
  as.numeric(t(inov) %*% solve(S) %*% inov)
}

# Internal: compute 2x3 Jacobian of perspective projection at point r
.compute_H <- function(r, cam) {
  cr <- cam$trm[1:3, 1:3] %*% r + cam$trm[1:3, 4]
  X <- cr[1]; Y <- cr[2]; Z_c <- cr[3]
  f <- c(cam$km[1,1], cam$km[2,2])

  # dproj/dcr (2x3)
  dpdc <- matrix(c(f[1]/Z_c, 0, -f[1]*X/Z_c^2,
                   0, f[2]/Z_c, -f[2]*Y/Z_c^2), nrow = 2, byrow = TRUE)
  # chain rule: H = dpdc * R
  dpdc %*% cam$trm[1:3, 1:3]
}
