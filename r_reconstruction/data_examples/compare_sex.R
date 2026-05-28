# compare_sex.R
#
# Male vs Female flight-kinematics comparison for swarming Anopheles coluzzii.
#
# DATA SOURCE:
#   Feugere L, Gibson G, Roux O (2020). "Audio and 3D flight-track recordings of
#   mosquito responses to opposite-sex sound-stimuli." Dryad.
#   https://doi.org/10.5061/dryad.9cnp5hqhj
#
#   Files are sex-labelled by name:
#     *_mal-col_*  = MALE   Anopheles coluzzii recordings
#     *_fem-col_*  = FEMALE Anopheles coluzzii recordings
#
#   3D trajectories come from the Trackit system as the "...Splined.csv" files:
#     object; time; XSplined; YSplined; ZSplined; VXSplined; VYSplined; VZSplined
#     (positions in metres, velocities in m/s, semicolon-delimited)
#
# NOTE:
#   - Species: Anopheles coluzzii (same complex as An. gambiae s.s. in Butail 2012)
#   - Setting: laboratory swarms responding to sound stimuli
#   - Sex assignment is by RECORDING (each session is single-sex), which is the
#     experimental design of the source study. This is REAL, sex-labelled data.
#   - Sessions tagged "notSwarming" in the folder name are EXCLUDED.
#
# USAGE:
#   Rscript r_reconstruction/data_examples/compare_sex.R
#   -- or --
#   source("r_reconstruction/data_examples/compare_sex.R"); compare_sex()


# ============================================================
# SETUP
# ============================================================

.find_r_dir <- function() {
  for (d in c("r_reconstruction/R", "R", "../R")) {
    if (file.exists(file.path(d, "io.R"))) return(d)
  }
  stop("Cannot find R modules. Run from the mosquito_project root directory.")
}


# ============================================================
# DISCOVER SEX-LABELLED TRAJECTORY FILES
# ============================================================

#' Find all swarming Splined trajectory CSVs, labelled by sex
#'
#' @param data_root  directory containing the *_mal-col_* / *_fem-col_* folders
#' @return data frame with columns: path, sex, recording, session
find_sex_files <- function(data_root) {
  if (!dir.exists(data_root)) stop("Data root not found: ", data_root)

  files <- list.files(data_root, pattern = "PostProc_.*_Splined\\.csv$",
                      recursive = TRUE, full.names = TRUE)
  if (length(files) == 0) stop("No Splined CSV files found under ", data_root)

  # Exclude non-swarming sessions (folder name contains 'notswarming' / 'notSwarming')
  files <- files[!grepl("notswarming", files, ignore.case = TRUE)]

  sex <- ifelse(grepl("mal-col", files), "Male",
                ifelse(grepl("fem-col", files), "Female", NA))
  keep <- !is.na(sex)
  files <- files[keep]; sex <- sex[keep]

  # recording = the top-level dated folder; session = trackit timestamp folder
  recording <- sub("^.*/(\\d{6}_[a-z]+-col_d\\d+)/.*$", "\\1", files)
  session   <- basename(dirname(files))

  data.frame(path = files, sex = sex, recording = recording,
             session = session, stringsAsFactors = FALSE)
}


# ============================================================
# LOAD + CONVERT ONE FILE -> TIDY (positions mm, velocities mm/s)
# ============================================================

#' Load a single Trackit Splined CSV into tidy pipeline format
#'
#' @param path  CSV path
#' @param sex   "Male" or "Female"
#' @param file_idx  integer used to make track IDs globally unique
#' @param min_track_pts  drop tracks shorter than this many points
#' @return tidy data frame (target_id, frame, time_s, x, y, z, vx, vy, vz, sex, recording)
load_splined_file <- function(path, sex, file_idx, min_track_pts = 20) {
  raw <- tryCatch(
    read.csv(path, sep = ";", stringsAsFactors = FALSE, check.names = TRUE),
    error = function(e) { warning("Failed to read ", path, ": ", e$message); NULL }
  )
  if (is.null(raw) || nrow(raw) == 0) return(NULL)

  needed <- c("object", "time", "XSplined", "YSplined", "ZSplined",
              "VXSplined", "VYSplined", "VZSplined")
  if (!all(needed %in% names(raw))) {
    warning("Missing columns in ", basename(path)); return(NULL)
  }

  # Drop short tracks
  tab <- table(raw$object)
  keep_obj <- as.numeric(names(tab[tab >= min_track_pts]))
  raw <- raw[raw$object %in% keep_obj, ]
  if (nrow(raw) == 0) return(NULL)

  # Globally unique integer target IDs: file_idx * 100000 + object
  uid <- file_idx * 100000L + as.integer(raw$object)

  data.frame(
    target_id = uid,
    frame     = NA_integer_,                 # not needed for distribution comparison
    time_s    = raw$time,
    x  = raw$XSplined * 1000,                 # m -> mm
    y  = raw$YSplined * 1000,
    z  = raw$ZSplined * 1000,
    vx = raw$VXSplined * 1000,                # m/s -> mm/s
    vy = raw$VYSplined * 1000,
    vz = raw$VZSplined * 1000,
    sex       = sex,
    recording = sub("^.*/(\\d{6}_[a-z]+-col_d\\d+)/.*$", "\\1", path),
    stringsAsFactors = FALSE
  )
}


#' Load and pool all sex-labelled trajectories
#'
#' @param data_root  directory with the recording folders
#' @param max_files_per_sex  cap files per sex (NULL = all)
#' @return tidy data frame with sex + recording columns, plus speed columns
load_all_sex_data <- function(data_root, max_files_per_sex = NULL) {
  catalog <- find_sex_files(data_root)
  message(sprintf("[compare_sex] Found %d swarming files: %d Male, %d Female",
                  nrow(catalog),
                  sum(catalog$sex == "Male"), sum(catalog$sex == "Female")))

  if (!is.null(max_files_per_sex)) {
    catalog <- do.call(rbind, lapply(split(catalog, catalog$sex), function(d) {
      d[seq_len(min(nrow(d), max_files_per_sex)), ]
    }))
  }

  parts <- lapply(seq_len(nrow(catalog)), function(i) {
    load_splined_file(catalog$path[i], catalog$sex[i], file_idx = i)
  })
  tidy <- do.call(rbind, Filter(Negate(is.null), parts))
  rownames(tidy) <- NULL

  # Speed (m/s)
  tidy$speed_ms   <- sqrt(tidy$vx^2 + tidy$vy^2 + tidy$vz^2) / 1000
  tidy$speed_h_ms <- sqrt(tidy$vx^2 + tidy$vy^2) / 1000
  tidy$speed_v_ms <- abs(tidy$vz) / 1000
  tidy$height_cm  <- tidy$z / 10

  message(sprintf("[compare_sex] Pooled %s position records across %d tracks",
                  format(nrow(tidy), big.mark = ","),
                  length(unique(tidy$target_id))))
  tidy
}


# ============================================================
# PER-TRACK SUMMARY
# ============================================================

#' Summarise each track to one row (mean speed, height, duration, path length)
per_track_summary <- function(tidy) {
  do.call(rbind, lapply(split(tidy, tidy$target_id), function(d) {
    d <- d[order(d$time_s), ]
    n <- nrow(d)
    step_mm <- if (n > 1) sqrt(diff(d$x)^2 + diff(d$y)^2 + diff(d$z)^2) else numeric(0)
    path_mm <- sum(step_mm, na.rm = TRUE)
    data.frame(
      target_id   = d$target_id[1],
      sex         = d$sex[1],
      recording   = d$recording[1],
      n_points    = n,
      duration_s  = diff(range(d$time_s)),
      mean_speed_ms = mean(d$speed_ms, na.rm = TRUE),
      mean_height_cm = mean(d$height_cm, na.rm = TRUE),
      path_length_m  = path_mm / 1000,
      stringsAsFactors = FALSE
    )
  }))
}


# ============================================================
# COMPARISON FIGURES (base R)
# ============================================================

.sex_cols <- c(Male = "#1f78b4", Female = "#e31a1c")

#' Overlaid density comparison for one variable
.density_compare <- function(tidy, var, xlab, title, xlim = NULL) {
  m <- tidy[[var]][tidy$sex == "Male"]
  f <- tidy[[var]][tidy$sex == "Female"]
  m <- m[is.finite(m)]; f <- f[is.finite(f)]
  dm <- density(m); df <- density(f)
  if (is.null(xlim)) xlim <- range(c(dm$x, df$x))
  ylim <- c(0, max(c(dm$y, df$y)))

  plot(dm, col = .sex_cols["Male"], lwd = 2.5, main = title,
       xlab = xlab, xlim = xlim, ylim = ylim, las = 1)
  polygon(dm, col = adjustcolor(.sex_cols["Male"], 0.25), border = NA)
  polygon(df, col = adjustcolor(.sex_cols["Female"], 0.25), border = NA)
  lines(df, col = .sex_cols["Female"], lwd = 2.5)
  abline(v = mean(m), col = .sex_cols["Male"], lty = 2)
  abline(v = mean(f), col = .sex_cols["Female"], lty = 2)
  legend("topright", legend = c(sprintf("Male (mean %.2f)", mean(m)),
                                sprintf("Female (mean %.2f)", mean(f))),
         col = .sex_cols, lwd = 2.5, bty = "n", cex = 0.9)
}

#' Main comparison plot panel
plot_sex_comparison <- function(tidy, track_summary, save_path = NULL) {
  if (!is.null(save_path)) grDevices::pdf(save_path, width = 12, height = 9)
  old <- par(mfrow = c(2, 3), mar = c(4.5, 4.5, 3, 1), oma = c(0, 0, 2.5, 0))
  on.exit({ par(old); if (!is.null(save_path)) grDevices::dev.off() })

  # 1. Flight speed (per-point)
  .density_compare(tidy, "speed_ms", "Flight speed (m/s)",
                   "Instantaneous flight speed",
                   xlim = c(0, quantile(tidy$speed_ms, 0.99, na.rm = TRUE)))

  # 2. Flight height
  .density_compare(tidy, "height_cm", "Height relative to marker (cm)",
                   "Flight height")

  # 3. Horizontal vs vertical mean speed (grouped bars)
  agg_h <- tapply(tidy$speed_h_ms, tidy$sex, mean, na.rm = TRUE)
  agg_v <- tapply(tidy$speed_v_ms, tidy$sex, mean, na.rm = TRUE)
  mat <- rbind(Horizontal = agg_h[c("Male","Female")],
               Vertical   = agg_v[c("Male","Female")])
  barplot(mat, beside = TRUE, col = c("#6baed6", "#fb6a4a"),
          ylab = "Mean speed (m/s)", main = "Horizontal vs vertical speed",
          legend.text = c("Horizontal", "Vertical"),
          args.legend = list(x = "topright", bty = "n"), las = 1)

  # 4. Per-track mean speed boxplot
  boxplot(mean_speed_ms ~ sex, data = track_summary,
          col = adjustcolor(.sex_cols[c("Female","Male")], 0.5),
          ylab = "Per-track mean speed (m/s)", xlab = "",
          main = "Per-track mean speed", las = 1)

  # 5. Per-track mean height boxplot
  boxplot(mean_height_cm ~ sex, data = track_summary,
          col = adjustcolor(.sex_cols[c("Female","Male")], 0.5),
          ylab = "Per-track mean height (cm)", xlab = "",
          main = "Per-track mean height", las = 1)

  # 6. Track duration boxplot
  boxplot(duration_s ~ sex, data = track_summary,
          col = adjustcolor(.sex_cols[c("Female","Male")], 0.5),
          ylab = "Track duration (s)", xlab = "",
          main = "Track duration", las = 1)

  mtext("Anopheles coluzzii: Male vs Female flight kinematics (real data, Feugere et al. 2020)",
        outer = TRUE, cex = 1.1, font = 2)
}


# ============================================================
# SUMMARY STATISTICS TABLE
# ============================================================

sex_summary_table <- function(tidy, track_summary) {
  by_sex <- function(d, fn, col) tapply(d[[col]], d$sex, fn, na.rm = TRUE)
  data.frame(
    metric = c("N recordings", "N tracks", "N position records",
               "Mean flight speed (m/s)", "Median flight speed (m/s)",
               "Mean horizontal speed (m/s)", "Mean vertical speed (m/s)",
               "Mean height (cm)", "Mean track duration (s)",
               "Mean path length (m)"),
    Male = c(
      length(unique(track_summary$recording[track_summary$sex == "Male"])),
      sum(track_summary$sex == "Male"),
      sum(tidy$sex == "Male"),
      round(mean(tidy$speed_ms[tidy$sex == "Male"], na.rm = TRUE), 3),
      round(median(tidy$speed_ms[tidy$sex == "Male"], na.rm = TRUE), 3),
      round(mean(tidy$speed_h_ms[tidy$sex == "Male"], na.rm = TRUE), 3),
      round(mean(tidy$speed_v_ms[tidy$sex == "Male"], na.rm = TRUE), 3),
      round(mean(tidy$height_cm[tidy$sex == "Male"], na.rm = TRUE), 2),
      round(mean(track_summary$duration_s[track_summary$sex == "Male"], na.rm = TRUE), 2),
      round(mean(track_summary$path_length_m[track_summary$sex == "Male"], na.rm = TRUE), 2)
    ),
    Female = c(
      length(unique(track_summary$recording[track_summary$sex == "Female"])),
      sum(track_summary$sex == "Female"),
      sum(tidy$sex == "Female"),
      round(mean(tidy$speed_ms[tidy$sex == "Female"], na.rm = TRUE), 3),
      round(median(tidy$speed_ms[tidy$sex == "Female"], na.rm = TRUE), 3),
      round(mean(tidy$speed_h_ms[tidy$sex == "Female"], na.rm = TRUE), 3),
      round(mean(tidy$speed_v_ms[tidy$sex == "Female"], na.rm = TRUE), 3),
      round(mean(tidy$height_cm[tidy$sex == "Female"], na.rm = TRUE), 2),
      round(mean(track_summary$duration_s[track_summary$sex == "Female"], na.rm = TRUE), 2),
      round(mean(track_summary$path_length_m[track_summary$sex == "Female"], na.rm = TRUE), 2)
    ),
    stringsAsFactors = FALSE
  )
}


# ============================================================
# MAIN
# ============================================================

compare_sex <- function(
  data_root = "synthetic_data",
  out_dir   = "figures/sex_comparison",
  max_files_per_sex = NULL
) {
  r_dir <- .find_r_dir()
  source(file.path(r_dir, "io.R"),          local = FALSE)
  source(file.path(r_dir, "postprocess.R"), local = FALSE)
  source(file.path(r_dir, "plotting.R"),    local = FALSE)

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  message("\n=== Step 1: Load sex-labelled trajectories ===")
  tidy <- load_all_sex_data(data_root, max_files_per_sex = max_files_per_sex)

  message("\n=== Step 2: Per-track summary ===")
  track_summary <- per_track_summary(tidy)

  message("\n=== Step 3: Summary statistics ===")
  stats <- sex_summary_table(tidy, track_summary)
  print(stats, row.names = FALSE)

  message("\n=== Step 4: Statistical test (per-track mean speed) ===")
  tt <- tryCatch(
    t.test(mean_speed_ms ~ sex, data = track_summary),
    error = function(e) NULL
  )
  if (!is.null(tt)) {
    message(sprintf("  Welch t-test on per-track mean speed: t=%.2f, df=%.0f, p=%.4g",
                    tt$statistic, tt$parameter, tt$p.value))
  }

  message("\n=== Step 5: Figures ===")
  fig_main <- file.path(out_dir, "male_vs_female_kinematics.pdf")
  plot_sex_comparison(tidy, track_summary, save_path = fig_main)
  message("  Saved: ", fig_main)

  # Save CSV outputs
  write.csv(stats,         file.path(out_dir, "sex_summary_stats.csv"), row.names = FALSE)
  write.csv(track_summary, file.path(out_dir, "per_track_summary.csv"), row.names = FALSE)
  message("  Saved: sex_summary_stats.csv, per_track_summary.csv")

  message("\n=== Done. Output in: ", out_dir, " ===\n")
  invisible(list(tidy = tidy, track_summary = track_summary, stats = stats, ttest = tt))
}


if (!interactive() && identical(environment(), globalenv())) {
  compare_sex()
}
