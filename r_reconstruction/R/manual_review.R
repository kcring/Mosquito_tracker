# manual_review.R
# Lightweight CSV-driven manual tracklet review workflow.
#
# This module approximates the human-supervised component of the original
# MATLAB pipeline (trackone.m GUIDE GUI + glueTracks.m).
#
# WHAT THE ORIGINAL MATLAB STEP DOES:
#   - Operator loads data_mq_auto_*.mat in trackone.m
#   - Steps through frames using E/Q keyboard shortcuts
#   - Overlays tracks on raw images in two camera views
#   - Accepts tracklets as-is, rejects them, or combines/swaps them
#   - glueTracks.m implements the join/swap operations on the Xh matrix
#   - Final output: data_mq_EXP.mat (reviewed, combined tracks)
#
# PHASE 1 APPROXIMATION (Option A, lightweight):
#   - Read automated tracklets as a tidy data frame
#   - Write a review CSV with one row per (track_id, decision)
#   - Functions to apply the decisions: accept, reject, join, swap
#   - Plotting helpers to visualise tracks before/after review
#   - A documented Shiny spec stub
#
# HONEST LABEL: This is NOT equivalent to the original GUI review.
#   It enables the same operations but without synchronised image overlay.
#   For full visual review, the original MATLAB GUI or a future Shiny app
#   is required. See SHINY_SPEC below.
#
# Mirrors MATLAB source files:
#   tracker/trackone.m     (GUI driver — operations only, not display)
#   tracker/glueTracks.m   (join + swap operations on Xh)
#   tracker/t1_splice_tracks.m (splice tracklets)
#   tracker/t1_get_mqid.m  (retrieve track IDs)
#
# STATUS: translated (core operations), deferred (GUI)
# See audit/MATLAB_TO_R_CROSSWALK.csv and audit/UNRESOLVED_COMPONENTS.md.

# ============================================================
# REVIEW CSV FORMAT
# ============================================================

#' Create a review decision table from a list of track IDs
#'
#' Writes a CSV file for manual editing. The operator fills in:
#'   decision:      "accept" | "reject" | "join" | "swap"
#'   join_with:     target_id to join this track with (if decision=="join")
#'   join_from_frame: frame at which to start the join
#'   swap_with:     target_id to swap identity with (if decision=="swap")
#'   notes:         free text
#'
#' @param track_ids  integer vector of target IDs to review
#' @param output_csv  path to write the review CSV
#' @return data frame of the review table (also written to file)
create_review_table <- function(track_ids, output_csv) {
  df <- data.frame(
    target_id      = as.integer(track_ids),
    decision       = "accept",       # to be filled by reviewer
    join_with      = NA_integer_,
    join_from_frame = NA_integer_,
    swap_with      = NA_integer_,
    notes          = "",
    stringsAsFactors = FALSE
  )
  write.csv(df, output_csv, row.names = FALSE)
  message(sprintf("[manual_review.R] Review table written to: %s", output_csv))
  message("  Fill in 'decision' column: accept / reject / join / swap")
  message("  For 'join': set join_with = target_id and join_from_frame = frame number")
  message("  For 'swap': set swap_with = target_id")
  invisible(df)
}

#' Load a completed review CSV
#'
#' @param csv_path  path to filled review CSV
#' @return data frame with review decisions
load_review_decisions <- function(csv_path) {
  if (!file.exists(csv_path)) {
    stop("[manual_review.R] Review CSV not found: ", csv_path)
  }
  df <- read.csv(csv_path, stringsAsFactors = FALSE)
  required_cols <- c("target_id", "decision")
  missing <- setdiff(required_cols, names(df))
  if (length(missing) > 0) {
    stop("[manual_review.R] Review CSV missing columns: ", paste(missing, collapse=", "))
  }
  df$decision <- trimws(tolower(df$decision))
  valid_decisions <- c("accept", "reject", "join", "swap")
  bad <- setdiff(unique(df$decision), valid_decisions)
  if (length(bad) > 0) {
    warning("[manual_review.R] Unknown decisions (treated as 'accept'): ",
            paste(bad, collapse=", "))
  }
  df
}


# ============================================================
# TRACK OPERATIONS ON XH (mirrors glueTracks.m, t1_splice_tracks.m)
# ============================================================

#' Apply all review decisions to the Xh matrix
#'
#' Applies decisions in order: swaps first, then joins, then rejects.
#' Mirrors the combined effect of trackone.m GUI operations + glueTracks.m.
#'
#' @param Xh  state matrix [6*Nmax x Nframes]
#' @param decisions  data frame from load_review_decisions()
#' @return updated Xh
apply_review_decisions <- function(Xh, decisions) {
  # 1. Swaps
  swap_rows <- decisions[decisions$decision == "swap" & !is.na(decisions$swap_with), ]
  already_swapped <- integer(0)
  for (i in seq_len(nrow(swap_rows))) {
    t1 <- swap_rows$target_id[i]
    t2 <- swap_rows$swap_with[i]
    if (t1 %in% already_swapped || t2 %in% already_swapped) next
    Xh <- swap_tracks(Xh, t1, t2)
    already_swapped <- c(already_swapped, t1, t2)
    message(sprintf("[manual_review.R] Swapped target %d <-> %d", t1, t2))
  }

  # 2. Joins
  join_rows <- decisions[decisions$decision == "join" & !is.na(decisions$join_with), ]
  for (i in seq_len(nrow(join_rows))) {
    t_src  <- join_rows$target_id[i]
    t_dst  <- join_rows$join_with[i]
    t_from <- if (!is.na(join_rows$join_from_frame[i])) join_rows$join_from_frame[i] else 1L
    Xh <- join_tracks(Xh, src_id = t_src, dst_id = t_dst, from_frame = t_from)
    message(sprintf("[manual_review.R] Joined target %d -> %d from frame %d",
                    t_src, t_dst, t_from))
  }

  # 3. Rejects
  reject_rows <- decisions[decisions$decision == "reject", ]
  for (i in seq_len(nrow(reject_rows))) {
    t_id <- reject_rows$target_id[i]
    Xh <- reject_track(Xh, t_id)
    message(sprintf("[manual_review.R] Rejected target %d", t_id))
  }

  Xh
}

#' Get the row range in Xh for a given target ID
#'
#' Thin wrapper around io.R::get_xh_row_range to avoid dependency issues.
.rows <- function(target_id) {
  base <- (target_id - 1) * 6
  base + 1:6
}

#' Join two tracklets: copy rows from src into dst starting at from_frame
#'
#' Mirrors glueTracks.m join / t1_splice_tracks.m splice logic:
#'   Xh[dst_rows, from_frame:end] <- Xh[src_rows, from_frame:end]
#'   Xh[src_rows, from_frame:end] <- 0
#'
#' @param Xh  state matrix
#' @param src_id  source target ID (will be zeroed from from_frame onwards)
#' @param dst_id  destination target ID (receives data from src)
#' @param from_frame  frame at which to start the splice
#' @return updated Xh
join_tracks <- function(Xh, src_id, dst_id, from_frame = 1L) {
  n_frames <- ncol(Xh)
  if (from_frame > n_frames) {
    warning("[manual_review.R] join_tracks: from_frame > n_frames, no action taken")
    return(Xh)
  }
  frame_range <- from_frame:n_frames
  src_rows <- .rows(src_id)
  dst_rows <- .rows(dst_id)

  Xh[dst_rows, frame_range] <- Xh[src_rows, frame_range]
  Xh[src_rows, frame_range] <- 0
  Xh
}

#' Swap track identities between two targets (all frames)
#'
#' Mirrors glueTracks.m swap logic:
#'   tmp = Xh[t2_rows, :]
#'   Xh[t2_rows, :] = Xh[t1_rows, :]
#'   Xh[t1_rows, :] = tmp
#'
#' @param Xh  state matrix
#' @param t1_id, t2_id  target IDs to swap
#' @param from_frame  frame at which to start the swap (default 1 = all frames)
#' @return updated Xh
swap_tracks <- function(Xh, t1_id, t2_id, from_frame = 1L) {
  n_frames <- ncol(Xh)
  frame_range <- from_frame:n_frames
  r1 <- .rows(t1_id)
  r2 <- .rows(t2_id)

  tmp <- Xh[r2, frame_range, drop = FALSE]
  Xh[r2, frame_range] <- Xh[r1, frame_range]
  Xh[r1, frame_range] <- tmp
  Xh
}

#' Reject a tracklet: zero out its rows in Xh
#'
#' @param Xh  state matrix
#' @param target_id  target to reject
#' @param from_frame  start frame (default 1 = all frames)
#' @return updated Xh
reject_track <- function(Xh, target_id, from_frame = 1L) {
  frame_range <- from_frame:ncol(Xh)
  r <- .rows(target_id)
  Xh[r, frame_range] <- 0
  Xh
}

#' Get the ID and tracking span of all non-zero tracks in Xh
#'
#' Mirrors show_tracked_mq.m summary output.
#'
#' @param Xh  state matrix
#' @return data frame with target_id, first_frame, last_frame, n_frames_tracked
get_track_summary <- function(Xh) {
  n_targets <- nrow(Xh) %/% 6
  result <- do.call(rbind, lapply(seq_len(n_targets), function(t) {
    z_col <- Xh[(t-1)*6 + 3, ]   # z position row
    active <- which(z_col != 0)
    if (length(active) == 0) return(NULL)
    data.frame(
      target_id     = t,
      first_frame   = min(active),
      last_frame    = max(active),
      n_frames_tracked = length(active)
    )
  }))
  if (is.null(result)) return(data.frame())
  result
}


# ============================================================
# SHINY SPEC (future work stub)
# ============================================================

#' STUB: Launch minimal Shiny review app
#'
#' STATUS: deferred — Shiny app is not implemented in Phase 1.
#'
#' SPEC: What a Phase 2 Shiny review app would need:
#'
#' Panel layout:
#'   Left panel:   Left camera image with overlaid track projections
#'   Right panel:  Right camera image with overlaid track projections
#'   Bottom panel: 3D trajectory view (rgl or plotly)
#'   Side panel:   Track list table, decision buttons
#'
#' Interactions:
#'   - Slider: frame navigation (replaces E/Q keyboard)
#'   - Click track in image -> highlight in table
#'   - Button: Accept | Reject | Join (enter target_id + frame) | Swap
#'   - Button: Save decisions (write review CSV)
#'   - Speed filter: adjust threshold, rerun foreground for current frame
#'
#' Data requirements:
#'   - Image sequences accessible to the R session (or subsampled previews)
#'   - Xh matrix (auto-tracked)
#'   - Camera calibration structs (for projection overlay)
#'
#' Implementation estimate: ~3-5 days for a functional MVP.
#'
#' @param Xh  state matrix
#' @param cams  camera structs
#' @param frame_files  list of frame file path vectors per camera
launch_review_shiny <- function(Xh, cams, frame_files = NULL) {
  message("===================================================")
  message("[manual_review.R] Shiny review app: NOT YET IMPLEMENTED")
  message("")
  message("Phase 1 uses the CSV-driven workflow instead:")
  message("  1. summary  <- get_track_summary(Xh)")
  message("  2. decisions_csv <- 'review_decisions.csv'")
  message("  3. create_review_table(summary$target_id, decisions_csv)")
  message("  4. # Edit review_decisions.csv in a spreadsheet")
  message("  5. decisions <- load_review_decisions(decisions_csv)")
  message("  6. Xh_reviewed <- apply_review_decisions(Xh, decisions)")
  message("")
  message("For a visual alternative: use plotting.R::plot_trajectories_3d()")
  message("to inspect individual tracks before filling in the CSV.")
  message("===================================================")
  invisible(NULL)
}
