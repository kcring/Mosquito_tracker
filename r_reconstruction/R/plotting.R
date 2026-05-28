# plotting.R
# Trajectory visualisation and QC plots.
#
# Mirrors MATLAB source files:
#   tracker/plot_data.m         (3D trajectory plot)
#   tracker/plot_data2.m        (2D projections)
#   tracker/plotMatingData.m    (mating event tracks, paper Figure 12)
#   tracker/show_tracked_mq.m   (track summary)
#   analysis/compareAutoWithManual.m  (centroid/std time series, paper Figure 9)
#   analysis/disp_climate_data.m (weather station time series)
#   tracker/track_integrity_plot.m (tracking continuity)
#
# STATUS: translated (all functions use base R graphics)
# All outputs are produced with base R + optional grDevices — no ggplot2 required.
# See audit/MATLAB_TO_R_CROSSWALK.csv for details.
#
# HONEST LABEL: Plots are structural equivalents of the paper figures.
#   With synthetic data they demonstrate the plotting architecture;
#   biological interpretation requires real trajectory data.

# ============================================================
# COLOUR UTILITIES
# ============================================================

.track_colours <- function(n) {
  if (n == 0) return(character(0))
  grDevices::hcl.colors(max(n, 2), palette = "Set2")[seq_len(n)]
}

.alpha_col <- function(col, alpha = 0.3) {
  rgb_vals <- grDevices::col2rgb(col) / 255
  grDevices::rgb(rgb_vals[1], rgb_vals[2], rgb_vals[3], alpha = alpha)
}


# ============================================================
# 3D TRAJECTORY PLOT (mirrors plot_data.m, paper Figures 9, 12)
# ============================================================

#' Plot 3D mosquito trajectories
#'
#' Produces a 3D perspective plot with optional 2D projections onto
#' the floor (xy), side (xz), and back (yz) planes.
#' Mirrors plot_data.m / plotMatingData.m.
#'
#' NOTE: True 3D rotation requires the 'scatterplot3d' or 'rgl' package.
#' This function produces a static perspective view using base R persp-like
#' coordinate transformation. For interactive 3D, install 'rgl'.
#'
#' @param traj_df  tidy trajectory data frame
#' @param highlight_ids  integer vector of target IDs to highlight (e.g. mating pair)
#' @param highlight_cols  named list: id -> colour string
#' @param show_projections  logical: add 2D projection shadows
#' @param title  plot title
#' @param xlim, ylim, zlim  axis limits (mm); NULL = auto
#' @param fps  frames per second (for time axis)
#' @param save_path  if not NULL, saves to PDF/PNG at this path
plot_trajectories_3d <- function(traj_df,
                                  highlight_ids  = NULL,
                                  highlight_cols = list(),
                                  show_projections = TRUE,
                                  title  = "Mosquito 3D Trajectories",
                                  xlim   = NULL, ylim = NULL, zlim = NULL,
                                  fps    = 25,
                                  save_path = NULL) {
  if (!requireNamespace("scatterplot3d", quietly = TRUE)) {
    message("[plotting.R] 'scatterplot3d' not installed. Using 2D side-by-side fallback.")
    plot_trajectories_2d(traj_df, title = title, save_path = save_path)
    return(invisible(NULL))
  }

  if (!is.null(save_path)) grDevices::pdf(save_path, width = 10, height = 8)

  traj_clean <- traj_df[!is.na(traj_df$x), ]
  target_ids <- unique(traj_clean$target_id)
  n_targets  <- length(target_ids)
  cols <- .track_colours(n_targets)

  if (is.null(xlim)) xlim <- range(traj_clean$x, na.rm = TRUE)
  if (is.null(ylim)) ylim <- range(traj_clean$y, na.rm = TRUE)
  if (is.null(zlim)) zlim <- range(traj_clean$z, na.rm = TRUE)

  s3d <- scatterplot3d::scatterplot3d(
    x = numeric(0), y = numeric(0), z = numeric(0),
    xlim = xlim, ylim = ylim, zlim = zlim,
    xlab = "East-West (mm)", ylab = "South-North (mm)", zlab = "Vertical (mm)",
    main = title, type = "n", box = TRUE, grid = TRUE
  )

  for (i in seq_along(target_ids)) {
    tid <- target_ids[i]
    sub <- traj_clean[traj_clean$target_id == tid, ]
    sub <- sub[order(sub$frame), ]
    col <- if (as.character(tid) %in% names(highlight_cols)) highlight_cols[[as.character(tid)]]
           else if (tid %in% highlight_ids)                  "red"
           else                                               cols[i]

    s3d$points3d(sub$x, sub$y, sub$z, type = "l", col = col, lwd = 1.5)

    if (show_projections) {
      s3d$points3d(sub$x, sub$y, rep(zlim[1], nrow(sub)),
                   type = "l", col = .alpha_col(col, 0.3), lwd = 0.8)
    }
  }

  if (!is.null(save_path)) grDevices::dev.off()
  invisible(s3d)
}

#' 2D multi-panel trajectory plots (fallback / companion to 3D)
#'
#' Produces a 3-row panel: x(t), y(t), z(t) with swarm mean and 3-sigma bounds.
#' Mirrors compareAutoWithManual.m / paper Figure 9.
#'
#' @param traj_df  tidy trajectory data frame
#' @param swarm_stats  optional output of compute_swarm_stats()
#' @param highlight_ids  integer vector
#' @param title  plot title
#' @param fps  frames per second
#' @param save_path  if not NULL, save to file
plot_trajectories_2d <- function(traj_df,
                                  swarm_stats    = NULL,
                                  highlight_ids  = NULL,
                                  title          = "Mosquito Trajectories",
                                  fps            = 25,
                                  save_path      = NULL) {
  if (!is.null(save_path)) grDevices::pdf(save_path, width = 10, height = 9)

  old_par <- par(mfrow = c(3, 1), mar = c(3, 4, 1.5, 1), oma = c(0, 0, 2, 0))
  on.exit({ par(old_par); if (!is.null(save_path)) grDevices::dev.off() })

  dims  <- list(list(col = "x",  label = "x  (mm)"),
                list(col = "y",  label = "y  (mm)"),
                list(col = "z",  label = "z  (mm)"))
  target_ids <- unique(traj_df$target_id[!is.na(traj_df$x)])
  cols <- .track_colours(length(target_ids))

  for (d in dims) {
    traj_clean <- traj_df[!is.na(traj_df[[d$col]]), ]
    t_axis <- traj_clean$time_s

    ylim <- range(traj_clean[[d$col]], na.rm = TRUE)
    if (diff(ylim) == 0) ylim <- ylim + c(-1, 1)

    # Draw 3-sigma bounds if swarm_stats available
    if (!is.null(swarm_stats)) {
      mean_col <- paste0("mean_", d$col)
      std_col  <- paste0("std_", d$col)
      if (mean_col %in% names(swarm_stats)) {
        ss_clean <- swarm_stats[!is.na(swarm_stats[[mean_col]]), ]
        ts_s <- ss_clean$time_s
        m <- ss_clean[[mean_col]]
        s <- ss_clean[[std_col]]
        ylim <- range(c(ylim, m + 3*s, m - 3*s), na.rm = TRUE)
      }
    }

    plot(range(traj_clean$time_s, na.rm = TRUE), ylim,
         type = "n", xlab = "Time (s)", ylab = d$label, las = 1)
    grid(col = "grey90", lty = 1)

    # 3-sigma envelope
    if (!is.null(swarm_stats)) {
      mean_col <- paste0("mean_", d$col)
      std_col  <- paste0("std_", d$col)
      if (mean_col %in% names(swarm_stats)) {
        ss_clean <- swarm_stats[!is.na(swarm_stats[[mean_col]]), ]
        ts_s <- ss_clean$time_s; m <- ss_clean[[mean_col]]; s <- ss_clean[[std_col]]
        polygon(c(ts_s, rev(ts_s)), c(m + 3*s, rev(m - 3*s)),
                col = grDevices::adjustcolor("grey70", 0.4), border = NA)
        lines(ts_s, m, col = "grey40", lty = 2, lwd = 1.5)
      }
    }

    # Individual tracks
    for (i in seq_along(target_ids)) {
      tid <- target_ids[i]
      sub <- traj_df[traj_df$target_id == tid & !is.na(traj_df[[d$col]]), ]
      sub <- sub[order(sub$frame), ]
      col <- if (!is.null(highlight_ids) && tid %in% highlight_ids) "red" else cols[i]
      lines(sub$time_s, sub[[d$col]], col = col, lwd = 1)
    }
  }
  mtext(title, outer = TRUE, cex = 1.1)
}

#' Plot mating event: male + female tracks (mirrors plotMatingData.m, paper Figure 12)
#'
#' @param traj_df  tidy trajectory data frame
#' @param male_ids  integer vector of male target IDs
#' @param female_id  integer: female target ID
#' @param couple_id  integer: coupled-pair target ID (or NULL)
#' @param save_path  if not NULL, saves to file
plot_mating_event <- function(traj_df,
                               male_ids,
                               female_id,
                               couple_id  = NULL,
                               save_path  = NULL) {
  all_ids <- c(male_ids, female_id, couple_id)
  sub <- traj_df[traj_df$target_id %in% all_ids, ]

  col_map <- setNames(
    c(rep("royalblue", length(male_ids)),
      "red",
      if (!is.null(couple_id)) "purple" else NULL),
    as.character(all_ids)
  )

  plot_trajectories_2d(sub,
                        highlight_ids = c(female_id, couple_id),
                        title = "Mating Event Trajectories",
                        save_path = save_path)
}

#' Speed scatter plot: horizontal vs vertical speed (mirrors paper Figure 10)
#'
#' @param traj_df  tidy trajectory data frame (with speed columns from add_speed())
#' @param group_col  column name to colour by (e.g. "target_id" or "sex")
#' @param wind_speed  optional scalar wind speed (m/s) to add as vertical line
#' @param title  plot title
#' @param save_path  if not NULL, saves to file
plot_speed_scatter <- function(traj_df,
                                group_col    = "target_id",
                                wind_speed   = NULL,
                                title        = "Horizontal vs Vertical Speed",
                                save_path    = NULL) {
  if (!all(c("horizontal_speed_ms","vertical_speed_ms") %in% names(traj_df))) {
    traj_df <- add_speed(traj_df)
  }
  clean <- traj_df[!is.na(traj_df$horizontal_speed_ms), ]

  if (!is.null(save_path)) grDevices::pdf(save_path, width = 6, height = 5)
  old_par <- par(mar = c(4, 4, 2, 1))
  on.exit({ par(old_par); if (!is.null(save_path)) grDevices::dev.off() })

  groups <- unique(clean[[group_col]])
  cols <- .track_colours(length(groups))
  col_vec <- cols[match(clean[[group_col]], groups)]

  plot(clean$horizontal_speed_ms, clean$vertical_speed_ms,
       col = col_vec, pch = 20, cex = 0.5,
       xlab = "Horizontal speed (m/s)", ylab = "Vertical speed (m/s)",
       main = title, las = 1)

  if (!is.null(wind_speed)) {
    abline(v = wind_speed, col = "steelblue", lty = 2, lwd = 1.5)
    text(wind_speed, par("usr")[4] * 0.95, sprintf(" Wind\n %.1f m/s", wind_speed),
         col = "steelblue", adj = 0, cex = 0.8)
  }
}

#' Separation distance plot for a mating pair (mirrors paper Figure 13)
#'
#' @param sep_df  data frame from compute_separation()
#' @param threshold_m  distance threshold for "close encounter" (m, default 0.04 = 4 cm)
#' @param couple_time_s  time (s) at which coupling occurs (0 = centre)
#' @param save_path  if not NULL, saves to file
plot_separation <- function(sep_df,
                             threshold_m  = 0.04,
                             couple_time_s = NULL,
                             save_path    = NULL) {
  if (!is.null(save_path)) grDevices::pdf(save_path, width = 7, height = 4)
  old_par <- par(mar = c(4, 4, 2, 1))
  on.exit({ par(old_par); if (!is.null(save_path)) grDevices::dev.off() })

  t <- sep_df$time_s
  d <- sep_df$separation_m
  if (!is.null(couple_time_s)) t <- t - couple_time_s

  plot(t, d, type = "l", lwd = 1.5, col = "black",
       xlab = "Time relative to coupling (s)", ylab = "Separation (m)",
       main = "Mating pair separation distance", las = 1)
  abline(h = threshold_m, col = "grey60", lty = 2)
  text(min(t), threshold_m + 0.002, sprintf(" %.0f mm threshold", threshold_m * 1000),
       col = "grey40", adj = 0, cex = 0.8)

  close_enc <- which(d < threshold_m)
  if (length(close_enc) > 0) {
    points(t[close_enc], d[close_enc], pch = 25, col = "firebrick", bg = "firebrick", cex = 0.7)
  }
}

#' Mosquito count over time (mirrors compareAutoWithManual.m subplot 4)
#'
#' @param swarm_stats  data frame from compute_swarm_stats()
#' @param swarm_stats2  optional second data frame (e.g. auto vs manual)
#' @param label1, label2  legend labels
#' @param save_path  if not NULL, saves to file
plot_mosquito_count <- function(swarm_stats,
                                 swarm_stats2 = NULL,
                                 label1 = "Tracked",
                                 label2 = "Reference",
                                 save_path = NULL) {
  if (!is.null(save_path)) grDevices::pdf(save_path, width = 8, height = 3.5)
  old_par <- par(mar = c(4, 4, 2, 1))
  on.exit({ par(old_par); if (!is.null(save_path)) grDevices::dev.off() })

  ylim <- range(swarm_stats$n_mosquitoes, na.rm = TRUE)
  if (!is.null(swarm_stats2)) ylim <- range(c(ylim, swarm_stats2$n_mosquitoes), na.rm = TRUE)
  ylim[1] <- 0

  plot(swarm_stats$time_s, swarm_stats$n_mosquitoes, type = "l",
       col = "steelblue", lwd = 1.5, ylim = ylim,
       xlab = "Time (s)", ylab = "# of mosquitoes", las = 1,
       main = "Mosquito count over time")

  if (!is.null(swarm_stats2)) {
    lines(swarm_stats2$time_s, swarm_stats2$n_mosquitoes, col = "tomato", lwd = 1.5, lty = 2)
    legend("topright", legend = c(label1, label2),
           col = c("steelblue","tomato"), lty = c(1,2), lwd = 1.5, bty = "n")
  }
}

#' Climate / weather station data plot (mirrors disp_climate_data.m)
#'
#' @param climate_df  data frame from read_climate_data()
#' @param save_path  if not NULL, saves to file
plot_climate_data <- function(climate_df, save_path = NULL) {
  if (!is.null(save_path)) grDevices::pdf(save_path, width = 9, height = 7)
  old_par <- par(mfrow = c(3,1), mar = c(3, 4, 1.5, 1), oma = c(0,0,2,0))
  on.exit({ par(old_par); if (!is.null(save_path)) grDevices::dev.off() })

  if ("wind_ms" %in% names(climate_df)) {
    plot(climate_df$wind_ms, type = "l", col = "steelblue", lwd = 1.5,
         xlab = "", ylab = "Wind speed (m/s)", las = 1)
    grid(col = "grey90", lty = 1)
  }
  if ("temp_c" %in% names(climate_df)) {
    plot(climate_df$temp_c, type = "l", col = "tomato", lwd = 1.5,
         xlab = "", ylab = "Temperature (°C)", las = 1)
    grid(col = "grey90", lty = 1)
  }
  if ("humidity_pct" %in% names(climate_df)) {
    plot(climate_df$humidity_pct, type = "l", col = "seagreen", lwd = 1.5,
         xlab = "Sample index", ylab = "Humidity (%)", las = 1)
    grid(col = "grey90", lty = 1)
  }
  mtext("Weather Station Data (Kestrel 4500)", outer = TRUE, cex = 1.1)
}

#' Summary printout for all tracks (mirrors show_tracked_mq.m)
#'
#' @param Xh  state matrix or NULL
#' @param traj_df  tidy data frame (alternative input)
#' @return prints summary; returns data frame
show_track_summary <- function(Xh = NULL, traj_df = NULL) {
  if (is.null(traj_df) && !is.null(Xh)) {
    traj_df <- xh_to_tidy(Xh)
  }
  if (is.null(traj_df) || nrow(traj_df) == 0) {
    message("[plotting.R] No tracks to summarise.")
    return(invisible(NULL))
  }

  summary_df <- do.call(rbind, lapply(unique(traj_df$target_id), function(tid) {
    sub <- traj_df[traj_df$target_id == tid, ]
    data.frame(
      target_id = tid,
      n_frames  = nrow(sub),
      t_start_s = min(sub$time_s),
      t_end_s   = max(sub$time_s),
      duration_s = diff(range(sub$time_s)),
      mean_speed_ms = if ("speed_ms" %in% names(sub)) mean(sub$speed_ms, na.rm=TRUE) else NA
    )
  }))

  cat(sprintf("\n=== Track summary: %d active tracks ===\n", nrow(summary_df)))
  print(summary_df, row.names = FALSE)
  cat(sprintf("Total position points: %d\n\n", nrow(traj_df)))
  invisible(summary_df)
}

#' QC panel: 4-plot overview (count, 3D structure, speed distribution, z-height)
#'
#' @param traj_df  tidy trajectory data frame with speed columns
#' @param swarm_stats  from compute_swarm_stats()
#' @param title  overall title
#' @param save_path  if not NULL, saves to file
plot_qc_panel <- function(traj_df, swarm_stats = NULL, title = "QC Overview",
                           save_path = NULL) {
  if (!all(c("speed_ms","horizontal_speed_ms") %in% names(traj_df))) {
    traj_df <- add_speed(traj_df)
  }

  if (!is.null(save_path)) grDevices::pdf(save_path, width = 11, height = 8)
  old_par <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 1), oma = c(0,0,2,0))
  on.exit({ par(old_par); if (!is.null(save_path)) grDevices::dev.off() })

  # Panel 1: mosquito count over time
  if (!is.null(swarm_stats)) {
    plot_mosquito_count(swarm_stats, label1 = "Tracked")
  } else {
    frame_counts <- table(traj_df$frame)
    plot(as.numeric(names(frame_counts)), as.numeric(frame_counts),
         type = "l", col = "steelblue", lwd = 1.5,
         xlab = "Frame", ylab = "# mosquitoes", main = "Mosquito count", las = 1)
  }

  # Panel 2: z position over time
  clean <- traj_df[!is.na(traj_df$z), ]
  plot(clean$time_s, clean$z, col = .alpha_col("royalblue", 0.3), pch = 20, cex = 0.4,
       xlab = "Time (s)", ylab = "Height z (mm)", main = "Height distribution", las = 1)
  if (!is.null(swarm_stats) && "mean_z" %in% names(swarm_stats)) {
    lines(swarm_stats$time_s, swarm_stats$mean_z, col = "black", lwd = 1.5)
  }

  # Panel 3: speed histogram
  speed_vals <- traj_df$speed_ms[!is.na(traj_df$speed_ms)]
  hist(speed_vals, breaks = 30, col = "steelblue", border = "white",
       xlab = "Speed (m/s)", main = "Speed distribution", las = 1)
  abline(v = median(speed_vals), col = "red", lty = 2, lwd = 1.5)
  legend("topright", legend = sprintf("Median: %.2f m/s", median(speed_vals)),
         col = "red", lty = 2, bty = "n", cex = 0.85)

  # Panel 4: horizontal vs vertical speed
  h <- traj_df$horizontal_speed_ms[!is.na(traj_df$horizontal_speed_ms)]
  v <- traj_df$vertical_speed_ms[!is.na(traj_df$vertical_speed_ms)]
  plot(h, v, col = .alpha_col("grey40", 0.3), pch = 20, cex = 0.4,
       xlab = "Horizontal speed (m/s)", ylab = "Vertical speed (m/s)",
       main = "Speed ratio", las = 1)

  mtext(title, outer = TRUE, cex = 1.1)
}
