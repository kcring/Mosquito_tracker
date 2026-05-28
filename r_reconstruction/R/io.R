# io.R
# Input/output and configuration for the mosquito tracking R reconstruction.
#
# Mirrors MATLAB source files:
#   config.m, initproc.m, scan_expfile.m, setupdir.m,
#   readOffCamCalib.m, get_cam_calib_*.m, strCam.m,
#   t1_save_tracks.m, getind.m, snipx.m, reformat_data.m,
#   get_climate_data.m
#
# STATUS: translated (geometry/IO) | partial (calibration .m parser, .mat loading)
# See audit/MATLAB_TO_R_CROSSWALK.csv for details.

# ============================================================
# CONFIGURATION (mirrors config.m + paper Table 2)
# ============================================================

#' Default tracking parameters
#'
#' Returns a list of parameters matching config.m, mht2.m, and paper Table 2.
#' All distances in mm, time in s, angles in degrees unless noted.
#'
#' @return named list of parameters
default_params <- function() {
  list(
    # --- camera / timing ---
    te          = 1/40,        # exposure time (s) [config.m: optTrack.te]
    dt          = 1/25,        # frame interval (s); paper says 25 fps
    fps         = 25,
    nc          = 2,           # number of cameras

    # --- foreground extraction (mirrors config.m + getZ.m) ---
    br0         = 3,           # sliding window half-width -> window = 2*br0+1 = 7 frames
    fg_is_dark  = TRUE,        # mosquitoes are dark on light background
    t_area      = c(20, 150),  # min/max blob area (pixels^2) [paper Table 2]
    binary_t    = 0.05,        # default intensity threshold (overridden by init_bg_params)
    bbox        = c(50, 50),   # search bounding box half-size for missing measurements

    # --- stereo / calibration ---
    te_epipolar = 0.5,         # epipolar constraint threshold [paper eq 3.2]

    # --- measurement model ---
    sigma_ep    = diag(c(4, 4)), # endpoint covariance [paper Table 2, pixels^2]

    # --- tracking / MHT ---
    t_gate      = 16,          # Mahalanobis gating threshold [paper Table 2]
    Ns          = 1,           # MHT N-scanback [paper Table 2]
    hyp_sort    = 3,           # max hypotheses to keep per cluster
    P_D         = 0.9,         # probability of detection [mht2.m]
    beta_f      = 0.01,        # false alarm density [mht2.m]
    beta_n      = 0.01,        # new target density [mht2.m]

    # --- particle filter ---
    Np          = 200,         # number of particles [paper Table 2]
    sigma_w     = 100,         # disturbance covariance (mm^2 s^-4) [paper Table 2]
    pos_init_std  = 5,         # new target position std (mm) [paper section 3.3]
    vel_init_std  = 500,       # new target velocity std (mm/s)
    speed_min     = 100,       # mm/s (uniform prior lower bound)
    speed_max     = 4000,      # mm/s (uniform prior upper bound)

    # --- swarm geometry ---
    swarm_d0    = 200,         # min swarm depth (mm) [config.m]
    swarm_s0    = 900,         # max swarm depth (mm)

    # --- output ---
    max_targets = 100          # pre-allocation for Xh [trackone.m optTrack.trackone.max_t]
  )
}


# ============================================================
# CAMERA STRUCTURE (mirrors strCam.m)
# ============================================================

#' Create an empty camera structure
#'
#' Fields match strCam.m:
#'   id  - string identifier
#'   km  - 3x3 intrinsic matrix (focal lengths, principal point)
#'   kc1 - first radial distortion coefficient (Bouguet model)
#'   kc2 - second radial distortion coefficient
#'   trm - 4x4 extrinsic transform [R | t; 0 0 0 1]
#'         Camera 1 = identity; camera 2 = relative pose
#'   P   - 3x4 projection matrix (computed from km and trm)
#'
#' @param id string camera identifier
#' @return named list (camera struct)
make_cam_struct <- function(id = "Camera") {
  list(
    id  = id,
    km  = diag(3),
    kc1 = 0,
    kc2 = 0,
    trm = diag(4),
    P   = cbind(diag(3), rep(0, 3))
  )
}

#' Compute the 3x4 projection matrix from intrinsics and extrinsics
#'
#' P = km * [R | t]  where [R | t] is the first 3 rows of trm
#'
#' @param cam camera struct
#' @return camera struct with P field updated
compute_projection_matrix <- function(cam) {
  Rt <- cam$trm[1:3, 1:4]
  cam$P <- cam$km %*% Rt
  cam
}


# ============================================================
# CALIBRATION I/O (mirrors readOffCamCalib.m + get_cam_calib_*.m)
# ============================================================

#' Parse a Bouguet-style calibration .m file (text-based)
#'
#' The MATLAB Bouguet calibration toolbox generates a function like:
#'   function cam = get_cam_calib_YYYYMMDD(camid)
#'   switch camid
#'   case 1
#'   cam.id = 'L_00';
#'   cam.km = [f11 f12 f13; f21 f22 f23; f31 f32 f33];
#'   cam.kc1 = -0.21;
#'   cam.kc2 = 0.45;
#'   cam.trm = [4x4 matrix];
#'   ...
#'
#' This function reads the text and extracts numeric fields for each camera.
#'
#' STATUS: partial - handles the format of get_cam_calib_mar252016.m;
#'   variations in multi-line layout may need adjustment.
#'   FALLBACK: supply calibration directly as an R list (see make_cam_struct).
#'
#' @param filepath path to the .m calibration file
#' @param nc number of cameras (default 2)
#' @return list of cam structs, length nc
read_cam_calib_m_file <- function(filepath, nc = 2) {
  if (!file.exists(filepath)) {
    stop(sprintf("[io.R] Calibration file not found: %s", filepath))
  }
  lines <- readLines(filepath)

  cams <- vector("list", nc)
  for (camid in seq_len(nc)) {
    cam <- make_cam_struct(sprintf("cam%d", camid))

    # find 'case <camid>' block
    case_line <- grep(sprintf("^case\\s+%d", camid), lines)
    if (length(case_line) == 0) {
      warning(sprintf("[io.R] Could not find case %d in %s", camid, filepath))
      cams[[camid]] <- cam
      next
    }

    # find next 'case' or 'end' to delimit the block
    remaining <- lines[case_line:(length(lines))]
    end_line_rel <- grep("^(case|end)", remaining)
    if (length(end_line_rel) > 1) {
      block <- remaining[1:(end_line_rel[2] - 1)]
    } else {
      block <- remaining
    }
    block_text <- paste(block, collapse = "\n")

    # extract id string
    id_match <- regmatches(block_text, regexpr("cam\\.id\\s*=\\s*'([^']*)'", block_text))
    if (length(id_match) > 0) {
      cam$id <- sub("cam\\.id\\s*=\\s*'([^']*)'", "\\1", id_match)
    }

    # extract scalar fields
    cam$kc1 <- .extract_scalar(block_text, "kc1")
    cam$kc2 <- .extract_scalar(block_text, "kc2")

    # extract 3x3 km matrix
    cam$km <- .extract_matrix(block_text, "km", 3, 3)

    # extract 4x4 trm matrix
    cam$trm <- .extract_matrix(block_text, "trm", 4, 4)

    cam <- compute_projection_matrix(cam)
    cams[[camid]] <- cam
  }
  cams
}

# Internal helper: extract a named scalar from text
.extract_scalar <- function(text, field) {
  pat <- sprintf("cam\\.%s\\s*=\\s*([+-]?[0-9.eE+-]+)", field)
  m <- regmatches(text, regexpr(pat, text))
  if (length(m) == 0 || m == "") return(0)
  as.numeric(sub(sprintf("cam\\.%s\\s*=\\s*", field), "", m))
}

# Internal helper: extract a named matrix from text
# Looks for cam.field=[...] with row separators ; or newline
.extract_matrix <- function(text, field, nrow, ncol) {
  # Find the matrix content between [ and ]
  pat <- sprintf("cam\\.%s\\s*=\\s*\\[([^\\]]+)\\]", field)
  m <- regmatches(text, regexpr(pat, text, perl = TRUE))
  if (length(m) == 0 || m == "") return(diag(nrow))

  content <- sub(sprintf("cam\\.%s\\s*=\\s*\\[", field), "", m)
  content <- sub("\\].*", "", content)

  # Split by rows (semicolons or newlines)
  rows <- strsplit(content, "[;\n]")[[1]]
  rows <- rows[nchar(trimws(rows)) > 0]

  mat <- matrix(0, nrow = nrow, ncol = ncol)
  for (i in seq_along(rows)) {
    if (i > nrow) break
    vals <- as.numeric(strsplit(trimws(rows[[i]]), "[,\\s]+")[[1]])
    vals <- vals[!is.na(vals)]
    mat[i, seq_along(vals)] <- vals
  }
  mat
}

#' Read calibration from a directory containing a *calib*.m file
#'
#' Mirrors readOffCamCalib.m: scans the directory for a file matching
#' '*calib*.m' and calls read_cam_calib_m_file on it.
#'
#' @param calib_dir path to calibration directory
#' @param nc number of cameras
#' @return list of cam structs
read_off_cam_calib <- function(calib_dir, nc = 2) {
  files <- list.files(calib_dir, pattern = "calib.*\\.m$", full.names = TRUE)
  if (length(files) == 0) {
    stop(sprintf("[io.R] No *calib*.m file found in %s", calib_dir))
  }
  read_cam_calib_m_file(files[1], nc = nc)
}

#' Create calibration from raw numeric values (no .m file needed)
#'
#' Use this as a fallback when the .m parser fails, or to supply
#' synthetic calibration in tests.
#'
#' @param km1,km2  3x3 intrinsic matrices
#' @param trm1,trm2  4x4 extrinsic transforms (trm1 usually identity)
#' @param kc1_1,kc2_1  distortion for cam 1
#' @param kc1_2,kc2_2  distortion for cam 2
#' @return list of 2 cam structs
make_stereo_calib <- function(km1, km2, trm1 = diag(4), trm2,
                               kc1_1 = 0, kc2_1 = 0,
                               kc1_2 = 0, kc2_2 = 0) {
  cam1 <- make_cam_struct("cam1")
  cam1$km <- km1; cam1$trm <- trm1; cam1$kc1 <- kc1_1; cam1$kc2 <- kc2_1
  cam1 <- compute_projection_matrix(cam1)

  cam2 <- make_cam_struct("cam2")
  cam2$km <- km2; cam2$trm <- trm2; cam2$kc1 <- kc1_2; cam2$kc2 <- kc2_2
  cam2 <- compute_projection_matrix(cam2)

  list(cam1, cam2)
}


# ============================================================
# EXPERIMENT DIRECTORY (mirrors scan_expfile.m, initproc.m)
# ============================================================

#' Read experiment name and image ID from expfile.txt
#'
#' Format: two whitespace-separated fields on the first line:
#'   <expname>  <image_id>
#'
#' @param exp_dir experiment root directory (contains calib/ subdirectory)
#' @return named list with expname and image_id
read_expfile <- function(exp_dir) {
  fpath <- file.path(exp_dir, "calib", "expfile.txt")
  if (!file.exists(fpath)) {
    stop(sprintf("[io.R] expfile.txt not found: %s\nRun setupdir first.", fpath))
  }
  txt <- read.table(fpath, header = FALSE, col.names = c("expname", "image_id"),
                    stringsAsFactors = FALSE)
  list(expname = txt$expname[1], image_id = txt$image_id[1])
}

#' Set up output directory structure under exp_dir/output/data/
#'
#' @param exp_dir experiment root directory
#' @return path to output/data directory (created if missing)
setup_dirs <- function(exp_dir) {
  out_dir <- file.path(exp_dir, "output", "data")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_dir
}

#' List image files for one camera
#'
#' @param frame_dir directory containing TIFF/BMP frames
#' @param image_id  prefix string for this camera (e.g. "L" or "R")
#' @return sorted character vector of full file paths
list_frame_files <- function(frame_dir, image_id = "") {
  pat <- if (nchar(image_id) > 0) paste0("^", image_id) else ""
  all_files <- list.files(frame_dir, full.names = TRUE)
  if (nchar(pat) > 0) all_files <- all_files[grepl(pat, basename(all_files))]
  sort(all_files)
}

#' Initialise experiment: load everything needed to start tracking
#'
#' @param exp_dir  experiment root (contains calib/, frames/)
#' @param frame_dir  directory with TIFF sequences; if NULL uses exp_dir/frames/
#' @param params   result of default_params() or modified list
#' @return named list with cams, F_mat, frame_lists, expname, data_dir, params
init_experiment <- function(exp_dir, frame_dir = NULL, params = default_params()) {
  if (is.null(frame_dir)) frame_dir <- file.path(exp_dir, "frames")
  exp_info <- tryCatch(read_expfile(exp_dir), error = function(e) {
    list(expname = basename(exp_dir), image_id = "")
  })

  calib_dir <- file.path(exp_dir, "calib")
  cams <- tryCatch(
    read_off_cam_calib(calib_dir, nc = params$nc),
    error = function(e) {
      warning("[io.R] Calibration load failed: ", conditionMessage(e),
              "\nUsing identity calibration. Supply real calibration before tracking.")
      lapply(seq_len(params$nc), function(i) make_cam_struct(sprintf("cam%d", i)))
    }
  )

  # Compute fundamental matrix from camera pair
  F_mat <- tryCatch({
    environment(compute_fundamental_matrix) # ensure stereo_match.R is loaded
    compute_fundamental_matrix(cams[[1]], cams[[2]])
  }, error = function(e) {
    warning("[io.R] Could not compute F: ", conditionMessage(e))
    matrix(0, 3, 3)
  })

  # List frame files for each camera
  cam_ids <- c("L", "R")
  frame_lists <- lapply(seq_len(params$nc), function(i) {
    list_frame_files(frame_dir, cam_ids[i])
  })

  data_dir <- setup_dirs(exp_dir)

  list(
    cams       = cams,
    F_mat      = F_mat,
    frame_lists = frame_lists,
    expname    = exp_info$expname,
    image_id   = exp_info$image_id,
    data_dir   = data_dir,
    params     = params
  )
}


# ============================================================
# XH MATRIX I/O (mirrors t1_save_tracks.m, reformat_data.m)
# ============================================================

#' Xh matrix indexing (mirrors getind.m, snipx.m)
#'
#' The Xh matrix stores all targets stacked vertically:
#'   rows 6(t-1)+1 : 6t  for target t
#'   columns  1 : Nframes
#'
#' @param n_states  number of states per target (6: x y z vx vy vz)
#' @param target_id  1-based target index
#' @param state_idx  which states to retrieve (default 1:6)
#' @return integer vector of row indices into Xh
get_xh_row_range <- function(target_id, state_idx = 1:6, n_states = 6) {
  base <- (target_id - 1) * n_states
  base + state_idx
}

#' Create an empty Xh matrix
#'
#' @param max_targets  maximum number of targets (rows = 6 * max_targets)
#' @param n_frames  number of frames (columns)
#' @return numeric matrix initialised to 0
make_xh <- function(max_targets = 100, n_frames = 1) {
  matrix(0, nrow = 6 * max_targets, ncol = n_frames)
}

#' Get state for target t at frame k from Xh
#'
#' @param Xh  state matrix
#' @param target_id  1-based target index
#' @param frame_k  1-based frame index (NULL = all frames)
#' @return numeric vector [x y z vx vy vz] or matrix [6 x n_frames]
get_target_state <- function(Xh, target_id, frame_k = NULL) {
  rows <- get_xh_row_range(target_id)
  if (is.null(frame_k)) Xh[rows, ] else Xh[rows, frame_k]
}

#' Set state for target t at frame k in Xh
#'
#' @param Xh  state matrix (modified in place by reference-like copy)
#' @param target_id  1-based target index
#' @param frame_k  1-based frame index
#' @param state  numeric vector of length 6
#' @return updated Xh
set_target_state <- function(Xh, target_id, frame_k, state) {
  rows <- get_xh_row_range(target_id)
  Xh[rows, frame_k] <- state
  Xh
}

#' Convert Xh matrix to tidy long-format data frame
#'
#' Rows with all-zero state are treated as "not tracked" and set to NA.
#'
#' @param Xh  state matrix [6*Nmax x Nframes]
#' @param fps  frames per second (for time column)
#' @return data frame with columns: target_id, frame, time_s, x, y, z, vx, vy, vz
xh_to_tidy <- function(Xh, fps = 25) {
  n_targets <- nrow(Xh) %/% 6
  n_frames  <- ncol(Xh)
  frames    <- seq_len(n_frames)

  result <- do.call(rbind, lapply(seq_len(n_targets), function(t) {
    rows <- get_xh_row_range(t)
    if (max(rows) > nrow(Xh)) return(NULL)
    states <- Xh[rows, , drop = FALSE]  # 6 x n_frames
    tracked <- apply(states, 2, function(col) any(col != 0))
    if (!any(tracked)) return(NULL)

    df <- data.frame(
      target_id = t,
      frame     = frames,
      time_s    = (frames - 1) / fps,
      x  = states[1, ],
      y  = states[2, ],
      z  = states[3, ],
      vx = states[4, ],
      vy = states[5, ],
      vz = states[6, ]
    )
    df[!tracked, c("x","y","z","vx","vy","vz")] <- NA
    df
  }))

  result[!is.na(result$x), ]  # drop NA rows for compactness
}

#' Save tracklets to RDS and CSV
#'
#' R equivalent of t1_save_tracks.m:
#'   data_mq_auto_<expname>.rds  — Xh matrix + metadata
#'   data_mq_auto_<expname>.csv  — tidy long-format trajectories
#'
#' @param Xh  state matrix
#' @param expname  experiment name string
#' @param data_dir  output directory
#' @param cams  camera structs
#' @param params  parameter list
#' @param label  "auto" for automated output, "" for reviewed
save_tracklets <- function(Xh, expname, data_dir, cams = NULL,
                            params = NULL, label = "auto") {
  prefix <- if (nchar(label) > 0) sprintf("data_mq_%s_%s", label, expname)
            else sprintf("data_mq_%s", expname)
  rds_path <- file.path(data_dir, paste0(prefix, ".rds"))
  csv_path <- file.path(data_dir, paste0(prefix, ".csv"))

  saveRDS(list(Xh = Xh, cams = cams, params = params, expname = expname), rds_path)

  tidy_df <- xh_to_tidy(Xh, fps = if (!is.null(params)) params$fps else 25)
  write.csv(tidy_df, csv_path, row.names = FALSE)

  message(sprintf("[io.R] Saved %s (RDS + CSV)", prefix))
  invisible(list(rds = rds_path, csv = csv_path))
}

#' Load tracklets from RDS
#'
#' @param rds_path  path to .rds file saved by save_tracklets
#' @return list with Xh, cams, params, expname
load_tracklets <- function(rds_path) {
  if (!file.exists(rds_path)) stop("[io.R] File not found: ", rds_path)
  readRDS(rds_path)
}


# ============================================================
# CLIMATE DATA (mirrors get_climate_data.m, disp_climate_data.m)
# ============================================================

#' Read Kestrel weather station CSV
#'
#' The Kestrel 4500 outputs a CSV with columns including
#' DateTime, WindSpeed, Temperature, Humidity.
#'
#' @param csv_path  path to Kestrel CSV file
#' @return data frame with columns: time, wind_ms, temp_c, humidity_pct
read_climate_data <- function(csv_path) {
  if (!file.exists(csv_path)) stop("[io.R] Climate file not found: ", csv_path)
  raw <- read.csv(csv_path, stringsAsFactors = FALSE)

  # Normalise common Kestrel column names
  colmap <- list(
    time        = c("DateTime", "Date", "Time", "Timestamp"),
    wind_ms     = c("WindSpeed", "Wind.Speed..m.s.", "Wind_Speed_ms"),
    temp_c      = c("Temperature", "Temp..C.", "Temp_C"),
    humidity_pct = c("Humidity", "RH....", "Humidity_pct")
  )

  out <- data.frame(row.names = seq_len(nrow(raw)))
  for (dest in names(colmap)) {
    found <- intersect(colmap[[dest]], names(raw))
    out[[dest]] <- if (length(found) > 0) raw[[found[1]]] else NA
  }
  out
}
