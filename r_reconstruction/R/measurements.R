# measurements.R
# Streak measurement extraction from foreground blobs.
#
# Mirrors MATLAB source files:
#   tracker/strZ.m         (measurement struct definition)
#   tracker/setEndPointVelocities.m  (endpoint extraction)
#   tracker/adaptive_thresholding.m  (missing measurement search)
#   tracker/lineFit.m      (1D line fitting)
#   tracker/normv.m        (vector normalisation)
#   tracker/tra2b.m        (frame transform)
#   tracker/mydiff.m       (finite differences)
#   tracker/sma.m          (simple moving average)
#
# STATUS: translated
# See audit/MATLAB_TO_R_CROSSWALK.csv for details.

# ============================================================
# MEASUREMENT STRUCT (mirrors strZ.m)
# ============================================================

#' Create an empty measurement struct (mirrors strZ)
#'
#' Fields:
#'   u         [2x1] midpoint pixel coordinates [col, row] = [x, y]
#'   v         [2x1] orientation vector (major axis direction)
#'   sigma     [2x1] std dev along major and minor axes
#'   ep        [2x2] endpoints [e_start | e_end], each column [col,row]
#'   area      scalar blob area in pixels^2
#'   length    scalar streak length in pixels
#'   pixel_list [N x 2] pixel positions [row, col]
#'   id        scalar: linear index of midpoint in image
#'
#' @param n  number of measurements to pre-allocate (default 1)
#' @return single measurement list or list of n empty measurement lists
make_Z_struct <- function(n = 1) {
  empty <- list(
    u          = c(0, 0),
    v          = c(0, 0),
    sigma      = c(1, 1),
    ep         = matrix(0, 2, 2),
    area       = 0,
    length     = 0,
    pixel_list = matrix(0, 0, 2),
    id         = 0
  )
  if (n == 1) return(empty)
  replicate(n, empty, simplify = FALSE)
}


# ============================================================
# STREAK MEASUREMENT EXTRACTION (mirrors getZ.m + setEndPointVelocities.m)
# ============================================================

#' Extract streak measurements from a list of blobs
#'
#' Applies area filter, converts each valid blob to a Z measurement struct
#' with midpoint, endpoints, orientation, and length.
#'
#' Mirrors getZ.m + setEndPointVelocities.m combined.
#'
#' @param blobs  list of blob structs from extract_blobs()
#' @param t_area  c(min_area, max_area) pixel^2 filter
#' @param img_dims  c(H, W) image dimensions (for id computation)
#' @return list of Z measurement structs
extract_measurements <- function(blobs, t_area = c(20, 150), img_dims = c(1024, 1392)) {
  if (length(blobs) == 0) return(list())

  Z_list <- list()
  for (blob in blobs) {
    if (blob$area < t_area[1] || blob$area >= t_area[2]) next

    Z <- make_Z_struct()
    Z$u    <- blob$centroid  # [col, row]
    Z$area <- blob$area
    Z$pixel_list <- blob$pixel_list  # [N x 2] row,col

    # Orientation vector (from major axis angle)
    theta_rad <- blob$orientation * pi / 180
    Z$v <- c(blob$major_axis * cos(theta_rad),
             -blob$major_axis * sin(theta_rad))

    # Sigma from major/minor axis lengths
    Z$sigma <- c(
      sqrt(abs(blob$major_axis * cos(theta_rad) - blob$minor_axis * sin(theta_rad))),
      sqrt(abs(blob$major_axis * sin(theta_rad) + blob$minor_axis * cos(theta_rad)))
    )
    Z$sigma[Z$sigma < 0.5] <- 0.5  # minimum sigma to avoid singular matrices

    # Extract endpoints using streak model
    Z <- set_endpoint_velocities(Z)

    # ID = linear index of centroid pixel in [H x W] image
    H <- img_dims[1]; W <- img_dims[2]
    row_idx <- min(H, max(1, round(Z$u[2])))
    col_idx <- min(W, max(1, round(Z$u[1])))
    Z$id <- (col_idx - 1) * H + row_idx  # col-major linear index

    Z_list[[length(Z_list) + 1]] <- Z
  }
  Z_list
}

#' Fit a line through streak pixels and extract endpoints
#'
#' Mirrors setEndPointVelocities.m:
#'   1. Transform pixels into body frame aligned with major axis
#'   2. Fit line y = p0 + p1*x in body frame
#'   3. Project back to image frame
#'   4. Take min/max x as endpoints
#'   5. Update midpoint as intersection of fitted line with x=0
#'
#' @param Z  measurement struct with pixel_list, u, v fields
#' @return updated Z with ep, length fields set
set_endpoint_velocities <- function(Z) {
  if (nrow(Z$pixel_list) < 2) {
    # not enough pixels; use centroid as both endpoints
    Z$ep <- cbind(Z$u, Z$u)
    Z$length <- 0
    return(Z)
  }

  # pixel coords: col=x, row=y in image convention
  x_pix <- Z$pixel_list[, 2]  # col indices = x
  y_pix <- Z$pixel_list[, 1]  # row indices = y

  # Build body frame: x-axis along streak direction, origin at centroid
  xaxis <- normalize_vec(Z$v)
  yaxis <- c(-xaxis[2], xaxis[1])  # 90-degree rotation

  # World-to-body transform: 3x3 homogeneous [xaxis | yaxis | u]
  wTb <- rbind(c(xaxis, Z$u[1]),
               c(yaxis, Z$u[2]),
               c(0, 0, 1))
  bTw <- solve(wTb)

  # Transform pixel coords to body frame
  pts_w <- rbind(x_pix, y_pix)
  pts_b <- transform_points_2d(pts_w, bTw)
  xB <- pts_b[1, ]
  yB <- pts_b[2, ]

  # Fit line y = p[1] + p[2]*x in body frame
  p <- line_fit(xB, yB)
  yB_fit <- p[1] + p[2] * xB

  # Sort by x-coordinate; find min/max extent
  ord  <- order(xB)
  xB_s <- xB[ord]
  yB_s <- yB_fit[ord]

  e1_idx <- 1
  e2_idx <- length(xB_s)

  # Apply simple moving average to smooth pixel path (mirrors sma.m with n=5)
  if (length(xB_s) >= 5) {
    xB_s <- moving_average(xB_s, 5)
    yB_s <- moving_average(yB_s, 5)
    e1_idx <- 1
    e2_idx <- length(xB_s)
  }

  # Transform endpoints back to world frame
  ep_b <- rbind(c(xB_s[e1_idx], xB_s[e2_idx]),
                c(yB_s[e1_idx], yB_s[e2_idx]))
  ep_w <- transform_points_2d(ep_b, wTb)

  Z$ep     <- ep_w              # [2 x 2]: each column is [col, row]
  Z$length <- sqrt(sum((ep_w[, 2] - ep_w[, 1])^2))

  # Update midpoint: project x=0 in body frame back to world
  y_at_zero <- p[1]  # p[1] + p[2]*0
  mid_b <- c(0, y_at_zero)
  mid_w <- wTb[1:2, 1:2] %*% mid_b + wTb[1:2, 3]
  Z$u <- mid_w

  Z
}

#' Simple 1D line fit y ~ x (mirrors lineFit.m)
#'
#' @param x  numeric vector
#' @param y  numeric vector (same length)
#' @return numeric vector c(intercept, slope)
line_fit <- function(x, y) {
  if (length(unique(x)) < 2) return(c(mean(y), 0))
  fit <- lm.fit(cbind(1, x), y)
  fit$coefficients
}

#' Normalise a 2D or 3D vector (mirrors normv.m)
#'
#' @param v  numeric vector
#' @return unit vector, or zero vector if norm is zero
normalize_vec <- function(v) {
  n <- sqrt(sum(v^2))
  if (n < 1e-10) return(v)
  v / n
}

#' Simple moving average (mirrors sma.m)
#'
#' @param x  numeric vector
#' @param n  window size (odd recommended)
#' @return smoothed vector (same length, edge-padded)
moving_average <- function(x, n = 5) {
  if (length(x) < n) return(x)
  k <- floor(n / 2)
  w <- rep(1 / n, n)
  stats::filter(x, w, sides = 2, circular = FALSE)
  # stats::filter returns NA at edges; replace with original values
  result <- as.numeric(stats::filter(x, w, sides = 2, circular = FALSE))
  result[is.na(result)] <- x[is.na(result)]
  result
}

#' Finite difference along a vector (mirrors mydiff.m)
#'
#' @param x  numeric vector or matrix (columns are differentiated)
#' @param n  order of difference (default 1)
#' @return differenced vector/matrix, same length (last entry repeated)
finite_diff <- function(x, n = 1) {
  if (is.vector(x)) {
    d <- diff(x, differences = n)
    c(d, rep(d[length(d)], n))
  } else {
    apply(x, 2, function(col) finite_diff(col, n))
  }
}

#' Transform 2D homogeneous points by a 3x3 matrix (mirrors tra2b.m)
#'
#' @param pts  2 x N matrix of [x; y] points
#' @param T    3x3 homogeneous transform matrix
#' @return 2 x N transformed points
transform_points_2d <- function(pts, T) {
  n <- ncol(pts)
  homo <- rbind(pts, rep(1, n))   # 3 x N
  res  <- T %*% homo              # 3 x N
  res[1:2, , drop = FALSE]
}


# ============================================================
# ADAPTIVE MISSING-MEASUREMENT SEARCH (mirrors adaptive_thresholding.m)
# ============================================================

#' Search for missing measurements in predicted gating region
#'
#' For each predicted target with no measurement within the gating volume,
#' lowers the threshold to 0.75 * binary_t and re-extracts within a
#' bounding box around the predicted pixel position.
#'
#' Mirrors adaptive_thresholding.m / the inner loop in getZ.m.
#'
#' STATUS: partial — logic translated; requires image array to be in memory.
#'
#' @param predicted_pixels  list of predicted [col, row] pixel positions (per camera)
#'                          for targets with no current measurement
#' @param bg  background matrix [H x W]
#' @param current  current frame matrix [H x W]
#' @param binary_t  current threshold
#' @param bbox  c(half_w, half_h) search box half-size in pixels (default c(50,50))
#' @param t_area  area filter c(min, max)
#' @param fg_is_dark  logical
#' @return list of additional Z measurement structs found
find_missing_measurements <- function(predicted_pixels, bg, current,
                                       binary_t, bbox = c(50, 50),
                                       t_area = c(20, 150),
                                       fg_is_dark = TRUE) {
  lowered_t <- 0.75 * binary_t
  H <- nrow(bg); W <- ncol(bg)
  found <- list()

  for (pred in predicted_pixels) {
    col0 <- round(pred[1]); row0 <- round(pred[2])

    # bounding box in image coordinates
    c1 <- max(1, col0 - bbox[1]); c2 <- min(W, col0 + bbox[1])
    r1 <- max(1, row0 - bbox[2]); r2 <- min(H, row0 + bbox[2])

    # extract sub-image
    bg_sub      <- bg[r1:r2, c1:c2]
    current_sub <- current[r1:r2, c1:c2]

    fg_sub <- compute_foreground(bg_sub, current_sub, lowered_t, fg_is_dark)
    blobs  <- extract_blobs(fg_sub)

    if (length(blobs) == 0) next

    # offset pixel coords back to full image frame
    for (b in blobs) {
      b$centroid[1] <- b$centroid[1] + c1 - 1
      b$centroid[2] <- b$centroid[2] + r1 - 1
      if (!is.null(b$pixel_list) && nrow(b$pixel_list) > 0) {
        b$pixel_list[, 1] <- b$pixel_list[, 1] + r1 - 1  # row
        b$pixel_list[, 2] <- b$pixel_list[, 2] + c1 - 1  # col
      }
    }

    new_Z <- extract_measurements(blobs, t_area = t_area, img_dims = c(H, W))
    found <- c(found, new_Z)
  }
  found
}
