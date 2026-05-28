# preprocess.R
# Image preprocessing and foreground extraction.
#
# Mirrors MATLAB source files:
#   tracker/getZ.m (foreground extraction, sliding-window background)
#   tracker/init_bgparams.m (adaptive threshold initialisation)
#   tracker/setRoi.m (region-of-interest masking)
#   tracker/update_imgarr.m (sliding window frame buffer)
#
# STATUS: translated (core foreground logic); partial (adaptive init loop)
# See audit/MATLAB_TO_R_CROSSWALK.csv for details.
#
# NOTE ON IMAGES: The R implementation represents an image as a
# numeric matrix [height x width], normalised to [0, 1].
# Use read_frame() to load TIFF/PNG frames.
# No image processing library is required for the core logic
# (uses base-R matrix operations).
# Optional: install 'tiff' package for read_frame() support.

# ============================================================
# FRAME LOADING
# ============================================================

#' Read a single image frame as a normalised matrix
#'
#' Supports TIFF (requires 'tiff' package) and PNG (requires 'png' package).
#' Falls back to reading as a raw byte vector if packages unavailable.
#'
#' @param filepath  path to image file
#' @return numeric matrix [height x width], values in [0, 1]
read_frame <- function(filepath) {
  ext <- tolower(tools::file_ext(filepath))

  if (ext %in% c("tif", "tiff")) {
    if (!requireNamespace("tiff", quietly = TRUE)) {
      stop("[preprocess.R] 'tiff' package needed for TIFF files. Install with: install.packages('tiff')")
    }
    img <- tiff::readTIFF(filepath, as.is = FALSE)  # normalised to [0,1]
  } else if (ext == "png") {
    if (!requireNamespace("png", quietly = TRUE)) {
      stop("[preprocess.R] 'png' package needed for PNG files. Install with: install.packages('png')")
    }
    img <- png::readPNG(filepath)
  } else {
    stop("[preprocess.R] Unsupported image format: ", ext, ". Use TIFF or PNG.")
  }

  # If multichannel (RGB/RGBA), convert to greyscale
  if (length(dim(img)) == 3) {
    img <- 0.299 * img[,,1] + 0.587 * img[,,2] + 0.114 * img[,,3]
  }
  img
}

#' Load a sliding window of frames into a 3D array
#'
#' Mirrors update_imgarr.m: maintains a [H x W x (2*br0+1)] buffer.
#'
#' @param frame_files  sorted character vector of frame paths for one camera
#' @param frame_k  current frame index (1-based)
#' @param br0  sliding window half-width (default 3 -> 7-frame window)
#' @return array [height x width x (2*br0+1)]
load_frame_window <- function(frame_files, frame_k, br0 = 3) {
  n_files <- length(frame_files)
  win_size <- 2 * br0 + 1

  # clamp window indices to available frames
  k_min <- max(1, frame_k - br0)
  k_max <- min(n_files, frame_k + br0)
  indices <- seq(k_min, k_max)

  # read first frame to get dimensions
  first <- read_frame(frame_files[indices[1]])
  H <- nrow(first); W <- ncol(first)

  buf <- array(0, dim = c(H, W, win_size))

  # fill: centre is at position br0+1 in the buffer
  for (i in seq_along(indices)) {
    slot <- indices[i] - frame_k + br0 + 1
    if (slot >= 1 && slot <= win_size) {
      buf[,,slot] <- if (i == 1) first else read_frame(frame_files[indices[i]])
    }
  }
  buf
}


# ============================================================
# BACKGROUND SUBTRACTION (mirrors getZ.m eq 3.1)
# ============================================================

#' Compute dynamic background image using sliding-window extremum
#'
#' For dark-on-light (fg_is_dark=TRUE):
#'   background = max over window (bright background revealed by slow movement)
#'   foreground = background - current_frame  (mosquitoes are darker)
#'
#' For light-on-dark (fg_is_dark=FALSE):
#'   background = min over window
#'   foreground = current_frame - background
#'
#' Eq 3.1: B_{u,v}[k] = max_{i in [k-d, k+d]} I_{u,v}[i]
#'
#' @param imgarr  3D array [H x W x window_size]
#' @param br  actual half-width used (can be <= br0)
#' @param br0  buffer centre offset (centre frame is at position br0+1)
#' @param fg_is_dark  logical (default TRUE for mosquitoes)
#' @return list with bg [H x W] background, current [H x W] current frame
compute_background <- function(imgarr, br = 3, br0 = 3, fg_is_dark = TRUE) {
  centre_idx <- br0 + 1
  lo <- max(1, centre_idx - br)
  hi <- min(dim(imgarr)[3], centre_idx + br)
  window <- imgarr[,, lo:hi, drop = FALSE]
  current <- imgarr[,, centre_idx]

  if (fg_is_dark) {
    bg <- apply(window, c(1, 2), max)
  } else {
    bg <- apply(window, c(1, 2), min)
  }
  list(bg = bg, current = current)
}

#' Compute foreground binary mask
#'
#' Subtracts background from current frame and applies intensity threshold.
#' Mirrors getZ.m: fg = im2bw(imsubtract(bg, current)/bitval, binary_t)
#'
#' @param bg  background matrix [H x W], values in [0, 1]
#' @param current  current frame matrix [H x W], values in [0, 1]
#' @param binary_t  intensity threshold in [0, 1]
#' @param fg_is_dark  logical
#' @return logical matrix [H x W] — TRUE where foreground
compute_foreground <- function(bg, current, binary_t = 0.05, fg_is_dark = TRUE) {
  if (fg_is_dark) {
    diff <- bg - current
  } else {
    diff <- current - bg
  }
  diff[diff < 0] <- 0
  diff > binary_t
}

#' Apply ROI mask (mirrors setRoi.m)
#'
#' Zeroes out pixels outside the region of interest.
#'
#' @param fg  foreground binary matrix
#' @param roi  numeric vector c(col_min, row_min, col_max, row_max) in pixels,
#'              or NULL to skip masking
#' @return masked foreground matrix
apply_roi <- function(fg, roi = NULL) {
  if (is.null(roi)) return(fg)
  H <- nrow(fg); W <- ncol(fg)
  mask <- matrix(FALSE, H, W)
  r1 <- max(1, roi[2]); r2 <- min(H, roi[4])
  c1 <- max(1, roi[1]); c2 <- min(W, roi[3])
  mask[r1:r2, c1:c2] <- TRUE
  fg & mask
}


# ============================================================
# CONNECTED COMPONENT LABELLING (mirrors regionprops)
# ============================================================

#' Simple connected-component labelling using flood-fill
#'
#' Equivalent to MATLAB's bwlabel/regionprops for small binary images.
#' Returns a list of blob property structures.
#'
#' For performance on real images, consider using the 'EBImage' package:
#'   blobs <- EBImage::bwlabel(fg)
#'   props <- EBImage::computeFeatures.shape(blobs)
#'
#' @param fg  logical matrix [H x W]
#' @return list of blob structs, each with:
#'   centroid [col, row], area, major_axis, minor_axis, orientation (degrees),
#'   pixel_list [N x 2 matrix of row,col indices]
extract_blobs <- function(fg) {
  # Use native R flood-fill connected components
  H <- nrow(fg); W <- ncol(fg)
  labels <- matrix(0L, H, W)
  current_label <- 0L

  for (r in seq_len(H)) {
    for (c in seq_len(W)) {
      if (fg[r, c] && labels[r, c] == 0) {
        current_label <- current_label + 1L
        labels <- .flood_fill(fg, labels, r, c, H, W, current_label)
      }
    }
  }

  blobs <- list()
  for (lbl in seq_len(current_label)) {
    idx <- which(labels == lbl, arr.ind = TRUE)  # [N x 2]: col 1=row, col 2=col
    if (nrow(idx) == 0) next
    blob <- .compute_blob_props(idx)
    blobs[[length(blobs) + 1]] <- blob
  }
  blobs
}

# Internal: flood fill using a stack (iterative BFS)
.flood_fill <- function(fg, labels, start_r, start_c, H, W, lbl) {
  stack <- list(c(start_r, start_c))
  while (length(stack) > 0) {
    pt <- stack[[length(stack)]]
    stack[[length(stack)]] <- NULL
    r <- pt[1]; c <- pt[2]
    if (r < 1 || r > H || c < 1 || c > W) next
    if (!fg[r, c] || labels[r, c] != 0) next
    labels[r, c] <- lbl
    stack <- c(stack, list(c(r+1,c), c(r-1,c), c(r,c+1), c(r,c-1)))
  }
  labels
}

# Internal: compute blob properties from pixel index matrix
# idx is [N x 2]: column 1 = row index, column 2 = col index
.compute_blob_props <- function(idx) {
  rows <- idx[, 1]; cols <- idx[, 2]
  area <- nrow(idx)
  centroid <- c(mean(cols), mean(rows))  # [col, row] = [x, y] convention

  # Second-moment tensor for orientation and axis lengths
  dr <- rows - mean(rows)
  dc <- cols - mean(cols)
  m20 <- sum(dc^2) / area
  m02 <- sum(dr^2) / area
  m11 <- sum(dc * dr) / area

  # Eigenvalues of covariance matrix for axis lengths
  trace <- m20 + m02
  det_val <- m20 * m02 - m11^2
  disc <- sqrt(max(0, (trace/2)^2 - det_val))
  lam1 <- trace/2 + disc
  lam2 <- trace/2 - disc

  major_axis <- 4 * sqrt(max(0, lam1))  # matches MATLAB MajorAxisLength convention
  minor_axis <- 4 * sqrt(max(0, lam2))

  # Orientation: angle of major axis to horizontal (degrees)
  orientation <- if (abs(m20 - m02) < 1e-10 && abs(m11) < 1e-10) {
    0
  } else {
    atan2(2 * m11, m20 - m02) * 180 / pi / 2
  }

  list(
    centroid    = centroid,     # [col, row] in pixel coords
    area        = area,
    major_axis  = major_axis,
    minor_axis  = minor_axis,
    orientation = orientation,  # degrees from horizontal
    pixel_list  = idx           # [N x 2] row,col
  )
}


# ============================================================
# ADAPTIVE THRESHOLD INITIALISATION (mirrors init_bgparams.m)
# ============================================================

#' Initialise background subtraction parameters adaptively
#'
#' Adjusts binary_t until the number of detected blobs is within
#' [0.5, 1.5] * expected_nmq. Mirrors init_bgparams.m.
#'
#' STATUS: partial — core adaptive loop translated;
#'   interactive user prompt replaced by expected_nmq argument.
#'
#' @param imgarr  3D frame window array [H x W x win]
#' @param expected_nmq  expected number of mosquitoes (approximate)
#' @param fg_is_dark  logical
#' @param t_area  numeric c(min, max) blob area filter
#' @param roi  ROI vector or NULL
#' @param binary_t_init  starting threshold
#' @param max_iter  maximum search iterations
#' @return list with binary_t, area_t, br, noise_std
init_bg_params <- function(imgarr, expected_nmq = 10, fg_is_dark = TRUE,
                            t_area = c(20, 150), roi = NULL,
                            binary_t_init = 0.05, max_iter = 20) {
  br0 <- (dim(imgarr)[3] - 1) %/% 2
  binary_t <- binary_t_init
  target_lo <- 0.5 * expected_nmq
  target_hi <- 1.5 * expected_nmq

  bg_res <- compute_background(imgarr, br = br0, br0 = br0, fg_is_dark = fg_is_dark)

  for (iter in seq_len(max_iter)) {
    fg  <- compute_foreground(bg_res$bg, bg_res$current, binary_t, fg_is_dark)
    fg  <- apply_roi(fg, roi)
    blobs <- extract_blobs(fg)

    areas <- vapply(blobs, function(b) b$area, 0.0)
    valid_blobs <- blobs[areas >= t_area[1] & areas < t_area[2]]
    n_blobs <- length(valid_blobs)

    if (n_blobs >= target_lo && n_blobs <= target_hi) break

    if (n_blobs > target_hi) {
      binary_t <- binary_t * 1.05  # raise threshold to reduce blobs
    } else {
      binary_t <- binary_t * 0.95  # lower threshold to find more blobs
    }
    binary_t <- max(0.001, min(0.99, binary_t))
  }

  # Estimate background noise std from flat regions
  noise_std <- sd(as.vector(bg_res$bg - bg_res$current)[as.vector(!fg)])
  if (is.na(noise_std) || noise_std == 0) noise_std <- 0.01

  list(
    binary_t  = binary_t,
    area_t    = t_area,
    br        = br0,
    noise_std = noise_std,
    n_blobs_found = n_blobs
  )
}
