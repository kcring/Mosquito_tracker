# app.R
# ============================================================
# Mosquito Tracklet Review — Shiny app (Phase 1, Option A)
# ============================================================
#
# A lightweight, INTERACTIVE replacement for the human-supervised review
# step of the original MATLAB GUI (trackone.m). It approximates the
# verify/combine workflow described in Butail et al. (2012, section 3.3)
# using the operations from manual_review.R (glueTracks.m).
#
# WHAT IT DOES (works on data we have today):
#   - CAMERA OVERLAY (trackone): steps through the synchronised stereo frames
#     and overlays the reconstructed 3D tracks projected back onto each camera
#     image (id-labelled, with track trails) + ground-truth markers. This is the
#     verify/combine review the original MATLAB trackone.m performed against the
#     photographs. Driven by the bundle from synthetic_image_demo.R.
#   - Loads a tidy trajectory CSV (target_id, frame, time_s, x, y, z, ...)
#   - Interactive 3D view of all tracklets (plotly), click/select to inspect
#   - Sortable track table (DT) with per-track stats
#   - Accept / Reject / Join / Swap operations applied live in tidy space
#   - Live preview of the edited tracks (before vs after)
#   - Saves a review-decisions CSV compatible with
#       manual_review.R::load_review_decisions() + apply_review_decisions()
#
# NOTE:
#   - The overlay images are SYNTHETIC (rendered by synthetic_image_demo.R),
#     because no public raw field image sequences exist. Every projection and
#     tracking step that runs on them is the real pipeline code. When real
#     calibrated footage exists, point the bundle at it and the same overlay works.
#
# RUN:
#   # from the mosquito_project root:
#   R -e "shiny::runApp('r_reconstruction/shiny_review', launch.browser = TRUE)"
#   # optionally set a data file first:
#   Sys.setenv(MOSQUITO_TRAJ_CSV = "figures/real_data_demo/tidy_trajectories_smoothed.csv")

suppressPackageStartupMessages({
  library(shiny)
  library(plotly)
  library(DT)
})

# Allow larger uploads than Shiny's 5 MB default (smoothed real-data CSVs can be big).
options(shiny.maxRequestSize = 200 * 1024^2)

# ------------------------------------------------------------
# Bring in the projection math so we can overlay 3D tracks on the
# camera images (the verify/combine step the MATLAB trackone.m did).
# ------------------------------------------------------------
.find_module <- function(fname) {
  cands <- c(file.path("..", "R", fname),                 # shiny sets wd = app dir
             file.path("r_reconstruction", "R", fname),   # run from project root
             file.path("R", fname))
  for (p in cands) if (file.exists(p)) return(p)
  NULL
}
.overlay_available <- FALSE
local({
  sm <- .find_module("stereo_match.R")
  if (!is.null(sm)) { source(sm, local = FALSE); .overlay_available <<- TRUE }
})
.has_png <- requireNamespace("png", quietly = TRUE)

#' Candidate locations of the image-review bundle produced by synthetic_image_demo.R
.default_bundle_candidates <- function() {
  c(Sys.getenv("MOSQUITO_IMAGE_BUNDLE", ""),
    file.path("..", "..", "figures", "synthetic_image_demo", "image_review_bundle.rds"),
    "figures/synthetic_image_demo/image_review_bundle.rds")
}
find_default_bundle <- function() {
  for (p in .default_bundle_candidates())
    if (nzchar(p) && file.exists(p)) return(normalizePath(p))
  ""
}

#' Load an image-review bundle and remember its directory (frames are relative to it).
load_bundle <- function(path) {
  b <- readRDS(path)
  b$dir <- dirname(normalizePath(path))
  b
}

#' Draw both camera frames for one frame index with tracks projected onto pixels.
render_overlay_frame <- function(bundle, frame_k, traj,
                                 highlight = NULL, show_truth = TRUE, trail = TRUE) {
  frame_k <- max(1, min(bundle$n_frames, as.integer(frame_k)))
  all_ids <- sort(unique(traj$target_id))
  pal     <- grDevices::hcl.colors(max(2, length(all_ids)), "Set2")
  idcol   <- function(id) pal[match(id, all_ids)]

  op <- par(mfrow = c(1, 2), mar = c(0.5, 0.5, 2.2, 0.5), bg = "black")
  on.exit(par(op))

  for (ci in 1:2) {
    camname <- c("cam1", "cam2")[ci]
    fpath   <- file.path(bundle$dir, bundle$frames[[camname]][frame_k])
    img <- if (.has_png && file.exists(fpath)) png::readPNG(fpath) else NULL
    if (!is.null(img) && length(dim(img)) == 3) img <- img[, , 1]
    W <- bundle$W; H <- bundle$H
    if (!is.null(img)) { H <- nrow(img); W <- ncol(img) }

    plot(NA, xlim = c(0, W), ylim = c(H, 0), asp = 1, axes = FALSE,
         xlab = "", ylab = "", col.main = "white",
         main = sprintf("Camera %d  —  frame %d / %d", ci, frame_k, bundle$n_frames))
    if (!is.null(img)) rasterImage(as.raster(img), 0, H, W, 0)
    cam <- bundle$cams[[ci]]

    # Ground-truth positions (synthetic demo only)
    if (show_truth && !is.null(bundle$truth)) {
      tf <- bundle$truth[bundle$truth$frame == frame_k, ]
      if (nrow(tf) > 0) {
        pr <- project_to_image(rbind(tf$x, tf$y, tf$z), cam)
        points(pr[1, ], pr[2, ], pch = 3, col = "#39FF14", cex = 1.9, lwd = 1.6)
      }
    }

    # Faint trailing path up to current frame, per track
    if (trail) {
      for (id in all_ids) {
        tr <- traj[traj$target_id == id & traj$frame <= frame_k, ]
        if (nrow(tr) < 2) next
        tr <- tr[order(tr$frame), ]
        pp <- project_to_image(rbind(tr$x, tr$y, tr$z), cam)
        lines(pp[1, ], pp[2, ], col = adjustcolor(idcol(id), 0.55), lwd = 1.4)
      }
    }

    # Current-frame track markers + id labels
    cf <- traj[traj$frame == frame_k, ]
    if (nrow(cf) > 0) {
      pr <- project_to_image(rbind(cf$x, cf$y, cf$z), cam)
      is_hl <- !is.null(highlight) & cf$target_id %in% highlight
      points(pr[1, ], pr[2, ], pch = 21, bg = idcol(cf$target_id),
             col = ifelse(is_hl, "#e31a1c", "white"),
             cex = ifelse(is_hl, 2.6, 1.8), lwd = ifelse(is_hl, 2.8, 1.2))
      text(pr[1, ], pr[2, ] - 6, labels = cf$target_id, col = "white", cex = 0.85)
    }
  }
}

# ------------------------------------------------------------
# Data loading helpers
# ------------------------------------------------------------

# Candidate default datasets (first that exists wins)
.default_csv_candidates <- function() {
  c(
    Sys.getenv("MOSQUITO_TRAJ_CSV", ""),
    "figures/real_data_demo/tidy_trajectories_smoothed.csv",
    "figures/sex_comparison/per_track_summary.csv",
    "r_reconstruction/data_examples/demo_output/SYNTHETIC_tracked_smoothed.csv"
  )
}

find_default_csv <- function() {
  for (p in .default_csv_candidates()) {
    if (nzchar(p) && file.exists(p)) return(normalizePath(p))
  }
  ""
}

#' Read a tidy trajectory CSV and standardise required columns.
read_traj_csv <- function(path) {
  df <- read.csv(path, stringsAsFactors = FALSE)
  # Be flexible about column names
  nm <- names(df)
  pick <- function(cands) { hit <- intersect(cands, nm); if (length(hit)) hit[1] else NA }
  c_id   <- pick(c("target_id", "id", "track_id"))
  c_x    <- pick(c("x", "X", "x_mm"))
  c_y    <- pick(c("y", "Y", "y_mm"))
  c_z    <- pick(c("z", "Z", "z_mm"))
  c_fr   <- pick(c("frame", "Frame"))
  c_t    <- pick(c("time_s", "time", "t"))
  if (any(is.na(c(c_id, c_x, c_y, c_z)))) {
    missing <- c(
      if (is.na(c_id)) "a track id (target_id / id / track_id)",
      if (is.na(c_x))  "x (or x_mm)",
      if (is.na(c_y))  "y (or y_mm)",
      if (is.na(c_z))  "z (or z_mm)"
    )
    stop(
      "This file can't be reviewed: it needs one row per (track, frame) ",
      "with 3D coordinates, but it is missing ", paste(missing, collapse = ", "), ".\n",
      "Columns found: ", paste(nm, collapse = ", "), ".\n",
      "Tip: use a per-frame trajectory file such as ",
      "'tidy_trajectories_smoothed.csv' or 'ground_truth_tracks.csv'. ",
      "Summary files (e.g. swarm_stats, per_track_summary, *_summary_stats) ",
      "do not contain per-frame positions and can't be loaded here.",
      call. = FALSE)
  }

  out <- data.frame(
    target_id = as.integer(factor(df[[c_id]])),       # contiguous ids
    orig_id   = df[[c_id]],
    x = df[[c_x]], y = df[[c_y]], z = df[[c_z]],
    stringsAsFactors = FALSE
  )
  out$frame  <- if (!is.na(c_fr)) df[[c_fr]] else ave(seq_len(nrow(out)), out$target_id, FUN = seq_along)
  out$time_s <- if (!is.na(c_t)) df[[c_t]] else out$frame
  # velocities if present
  for (v in c("vx", "vy", "vz")) if (v %in% nm) out[[v]] <- df[[v]]
  out
}

#' Per-track summary table.
track_table <- function(traj) {
  ids <- sort(unique(traj$target_id))
  do.call(rbind, lapply(ids, function(i) {
    d <- traj[traj$target_id == i, ]
    d <- d[order(d$time_s), ]
    n <- nrow(d)
    spd <- NA_real_
    if (all(c("vx","vy","vz") %in% names(d)))
      spd <- mean(sqrt(d$vx^2 + d$vy^2 + d$vz^2), na.rm = TRUE) / 1000
    else if (n > 1) {
      dt <- diff(d$time_s); dt[dt == 0] <- NA
      spd <- mean(sqrt(diff(d$x)^2 + diff(d$y)^2 + diff(d$z)^2) / dt, na.rm = TRUE) / 1000
    }
    data.frame(
      target_id   = i,
      n_points    = n,
      duration_s  = round(diff(range(d$time_s)), 2),
      mean_speed_ms = round(spd, 3),
      mean_x = round(mean(d$x), 1),
      mean_y = round(mean(d$y), 1),
      mean_z = round(mean(d$z), 1)
    )
  }))
}

# ------------------------------------------------------------
# Decision application in TIDY space (mirrors glueTracks.m semantics)
# ------------------------------------------------------------

apply_decisions_tidy <- function(traj, decisions) {
  out <- traj
  # 1. swaps (relabel ids for frames >= from)
  sw <- decisions[decisions$decision == "swap" & !is.na(decisions$swap_with), ]
  done <- integer(0)
  for (i in seq_len(nrow(sw))) {
    a <- sw$target_id[i]; b <- sw$swap_with[i]
    if (a %in% done || b %in% done) next
    fr <- if (!is.na(sw$join_from_frame[i])) sw$join_from_frame[i] else -Inf
    ia <- out$target_id == a & out$frame >= fr
    ib <- out$target_id == b & out$frame >= fr
    out$target_id[ia] <- b
    out$target_id[ib] <- a
    done <- c(done, a, b)
  }
  # 2. joins (relabel src -> dst for frames >= from)
  jn <- decisions[decisions$decision == "join" & !is.na(decisions$join_with), ]
  for (i in seq_len(nrow(jn))) {
    src <- jn$target_id[i]; dst <- jn$join_with[i]
    fr  <- if (!is.na(jn$join_from_frame[i])) jn$join_from_frame[i] else -Inf
    out$target_id[out$target_id == src & out$frame >= fr] <- dst
  }
  # 3. rejects (drop rows)
  rj <- decisions[decisions$decision == "reject", ]
  if (nrow(rj) > 0) out <- out[!(out$target_id %in% rj$target_id), ]
  out
}

empty_decisions <- function() {
  data.frame(target_id = integer(0), decision = character(0),
             join_with = integer(0), join_from_frame = integer(0),
             swap_with = integer(0), notes = character(0),
             stringsAsFactors = FALSE)
}

# ------------------------------------------------------------
# 3D plot
# ------------------------------------------------------------

plot3d_tracks <- function(traj, highlight = NULL, title = "Tracklets") {
  ids <- sort(unique(traj$target_id))
  pal <- grDevices::hcl.colors(max(2, length(ids)), "Set2")
  p <- plot_ly()
  for (k in seq_along(ids)) {
    d <- traj[traj$target_id == ids[k], ]
    d <- d[order(d$time_s), ]
    is_hl <- !is.null(highlight) && ids[k] %in% highlight
    p <- add_trace(p, x = d$x, y = d$y, z = d$z,
                   type = "scatter3d", mode = "lines+markers",
                   line = list(width = if (is_hl) 8 else 3,
                               color = if (is_hl) "#e31a1c" else pal[k]),
                   marker = list(size = if (is_hl) 4 else 2,
                                 color = if (is_hl) "#e31a1c" else pal[k]),
                   name = paste0("Track ", ids[k]),
                   hoverinfo = "name")
  }
  layout(p, title = title,
         scene = list(xaxis = list(title = "X (mm)"),
                      yaxis = list(title = "Y (mm)"),
                      zaxis = list(title = "Z (mm)"),
                      aspectmode = "data"))
}

# ============================================================
# UI
# ============================================================

ui <- fluidPage(
  titlePanel("Mosquito Tracklet Review (Phase 1 — Option A)"),
  tags$p(style = "color:#555;margin-top:-8px;",
         "Interactive approximation of the MATLAB trackone.m review step. ",
         tags$b("No raw image overlay"), " (public images unavailable) — the 3D view is the substitute."),
  sidebarLayout(
    sidebarPanel(
      width = 4,
      fileInput("csv", "Load trajectory CSV (optional)", accept = ".csv"),
      tags$details(
        tags$summary("Image overlay bundle (trackone)"),
        helpText("Auto-loads figures/synthetic_image_demo/image_review_bundle.rds if present.",
                 " Regenerate it with synthetic_image_demo()."),
        textInput("bundle_path", "Bundle .rds path (blank = default)", value = ""),
        actionButton("load_bundle", "Load image bundle")
      ),
      verbatimTextOutput("datainfo"),
      tags$hr(),
      h4("Track table"),
      DTOutput("tracks"),
      tags$hr(),
      h4("Decision for selected track"),
      radioButtons("decision", NULL,
                   choices = c("accept", "reject", "join", "swap"),
                   selected = "accept", inline = TRUE),
      conditionalPanel(
        "input.decision == 'join'",
        numericInput("join_with", "Join INTO target id", value = NA),
        numericInput("join_from_frame", "From frame (blank = all)", value = NA)
      ),
      conditionalPanel(
        "input.decision == 'swap'",
        numericInput("swap_with", "Swap identity with target id", value = NA),
        numericInput("swap_from_frame", "From frame (blank = all)", value = NA)
      ),
      textInput("notes", "Notes", ""),
      actionButton("add", "Add / update decision", class = "btn-primary"),
      actionButton("clear", "Clear all decisions"),
      tags$hr(),
      downloadButton("save", "Save decisions CSV", class = "btn-success")
    ),
    mainPanel(
      width = 8,
      tabsetPanel(
        tabPanel("Camera overlay (trackone)",
                 br(),
                 fluidRow(
                   column(6, sliderInput("frame", "Frame", min = 1, max = 1, value = 1,
                                         step = 1, width = "100%",
                                         animate = animationOptions(interval = 500))),
                   column(2, br(), actionButton("prev_frame", "\u25C0 Prev")),
                   column(2, br(), actionButton("next_frame", "Next \u25B6")),
                   column(2, br(), checkboxInput("show_truth", "Ground truth", TRUE))
                 ),
                 plotOutput("overlay", height = "430px"),
                 tags$p(style = "color:#777;font-size:12px;",
                        tags$b("Green +"), " = true position (synthetic demo). ",
                        tags$b("Filled circles"), " = reconstructed tracks (id labelled); ",
                        "red outline = selected in the table. Lines = track path so far. ",
                        "Pick a track in the table, then Join broken fragments or Reject false tracks."),
                 verbatimTextOutput("overlay_info")),
        tabPanel("3D — Original", plotlyOutput("p3d_orig", height = "520px")),
        tabPanel("3D — Preview (after decisions)", plotlyOutput("p3d_prev", height = "520px")),
        tabPanel("2D projections", plotOutput("p2d", height = "520px")),
        tabPanel("Decisions",
                 br(), DTOutput("dectable"),
                 br(), verbatimTextOutput("summary"))
      )
    )
  )
)

# ============================================================
# SERVER
# ============================================================

server <- function(input, output, session) {

  rv <- reactiveValues(traj = NULL, decisions = empty_decisions(), bundle = NULL)

  # Initial load: prefer an image bundle (so the overlay works out of the box),
  # otherwise fall back to a default tidy CSV.
  observe({
    if (is.null(rv$traj) && is.null(rv$bundle)) {
      bp <- find_default_bundle()
      if (.overlay_available && nzchar(bp)) {
        b <- tryCatch(load_bundle(bp), error = function(e) NULL)
        if (!is.null(b) && !is.null(b$recovered)) {
          rv$bundle <- b
          rv$traj   <- b$recovered
          rv$source <- basename(bp)
          return()
        }
      }
      p <- find_default_csv()
      if (nzchar(p)) {
        rv$traj <- tryCatch(read_traj_csv(p), error = function(e) NULL)
        rv$source <- p
      }
    }
  })

  # Keep the frame slider in sync with the loaded bundle
  observeEvent(rv$bundle, {
    req(rv$bundle)
    updateSliderInput(session, "frame", min = 1, max = rv$bundle$n_frames, value = 1)
  })

  # Load an image bundle on demand
  observeEvent(input$load_bundle, {
    if (!.overlay_available) {
      showNotification("Projection module (stereo_match.R) not found; overlay unavailable.",
                       type = "error"); return()
    }
    p <- input$bundle_path
    if (!nzchar(p)) p <- find_default_bundle()
    if (!nzchar(p) || !file.exists(p)) {
      showNotification("Bundle .rds not found. Run synthetic_image_demo() first.",
                       type = "error", duration = NULL); return()
    }
    b <- tryCatch(load_bundle(p),
                  error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL })
    if (is.null(b)) return()
    rv$bundle    <- b
    rv$traj      <- b$recovered
    rv$decisions <- empty_decisions()
    rv$source    <- basename(p)
    showNotification(
      sprintf("Loaded image bundle: %d frames, %d reconstructed tracks.",
              b$n_frames, length(unique(b$recovered$target_id))), type = "message")
  })

  # Frame navigation buttons
  observeEvent(input$prev_frame, {
    updateSliderInput(session, "frame", value = max(1, input$frame - 1))
  })
  observeEvent(input$next_frame, {
    mx <- if (!is.null(rv$bundle)) rv$bundle$n_frames else input$frame + 1
    updateSliderInput(session, "frame", value = min(mx, input$frame + 1))
  })

  # File upload
  observeEvent(input$csv, {
    req(input$csv)
    parsed <- tryCatch(
      read_traj_csv(input$csv$datapath),
      error = function(e) {
        showNotification(conditionMessage(e), type = "error", duration = NULL)
        NULL
      }
    )
    if (is.null(parsed)) return(invisible())   # keep previous data, app stays alive
    rv$traj <- parsed
    rv$decisions <- empty_decisions()
    rv$source <- input$csv$name
    showNotification(
      sprintf("Loaded %s: %d tracks, %d points.",
              input$csv$name, length(unique(parsed$target_id)), nrow(parsed)),
      type = "message")
  })

  output$datainfo <- renderText({
    if (is.null(rv$traj)) return("No data loaded. Upload a tidy trajectory CSV.")
    sprintf("Source: %s\nTracks: %d | Points: %d",
            basename(rv$source %||% "?"),
            length(unique(rv$traj$target_id)), nrow(rv$traj))
  })

  tt <- reactive({ req(rv$traj); track_table(rv$traj) })

  output$tracks <- renderDT({
    datatable(tt(), selection = "multiple", rownames = FALSE,
              options = list(pageLength = 8, scrollX = TRUE))
  })

  selected_ids <- reactive({
    sel <- input$tracks_rows_selected
    if (is.null(sel) || length(sel) == 0) return(NULL)
    tt()$target_id[sel]
  })

  # Add / update a decision for selected track(s)
  observeEvent(input$add, {
    ids <- selected_ids()
    if (is.null(ids)) {
      showNotification("Select one or more tracks in the table first.", type = "warning")
      return()
    }
    for (id in ids) {
      newrow <- data.frame(
        target_id       = as.integer(id),
        decision        = input$decision,
        join_with       = if (input$decision == "join") as.integer(input$join_with) else NA_integer_,
        join_from_frame = if (input$decision == "join") as.integer(input$join_from_frame)
                          else if (input$decision == "swap") as.integer(input$swap_from_frame)
                          else NA_integer_,
        swap_with       = if (input$decision == "swap") as.integer(input$swap_with) else NA_integer_,
        notes           = input$notes,
        stringsAsFactors = FALSE
      )
      d <- rv$decisions
      d <- d[d$target_id != id, ]          # replace any existing decision
      rv$decisions <- rbind(d, newrow)
    }
    showNotification(sprintf("Recorded '%s' for track(s): %s",
                             input$decision, paste(ids, collapse = ", ")),
                     type = "message")
  })

  observeEvent(input$clear, { rv$decisions <- empty_decisions() })

  output$dectable <- renderDT({
    datatable(rv$decisions, rownames = FALSE, options = list(pageLength = 10))
  })

  preview <- reactive({
    req(rv$traj)
    if (nrow(rv$decisions) == 0) return(rv$traj)
    apply_decisions_tidy(rv$traj, rv$decisions)
  })

  output$overlay <- renderPlot({
    if (is.null(rv$bundle)) {
      op <- par(bg = "white", mar = c(0,0,0,0)); on.exit(par(op))
      plot.new()
      text(0.5, 0.5, paste0(
        "No image bundle loaded.\n\n",
        "Run  source('r_reconstruction/data_examples/synthetic_image_demo.R'); demo_synthetic_images()\n",
        "then click 'Load image bundle' (or restart the app to auto-load)."),
        cex = 1.05, col = "#444")
      return()
    }
    render_overlay_frame(rv$bundle, input$frame, preview(),
                         highlight = selected_ids(),
                         show_truth = isTRUE(input$show_truth))
  })

  output$overlay_info <- renderPrint({
    if (is.null(rv$bundle)) { cat("No image bundle loaded.\n"); return() }
    cf <- preview()[preview()$frame == input$frame, ]
    cat(sprintf("Frame %d/%d | tracks visible this frame: %d | total tracks: %d\n",
                input$frame, rv$bundle$n_frames, nrow(cf),
                length(unique(preview()$target_id))))
    cat("Reminder: input images are SYNTHETIC; the projection + tracking code is the real pipeline.\n")
  })

  output$p3d_orig <- renderPlotly({
    req(rv$traj); plot3d_tracks(rv$traj, highlight = selected_ids(),
                                title = "Original tracklets")
  })

  output$p3d_prev <- renderPlotly({
    req(rv$traj); plot3d_tracks(preview(),
                                title = "After applying decisions")
  })

  output$p2d <- renderPlot({
    req(rv$traj)
    d <- rv$traj; ids <- sort(unique(d$target_id))
    pal <- grDevices::hcl.colors(max(2, length(ids)), "Set2")
    hl <- selected_ids()
    op <- par(mfrow = c(3, 1), mar = c(4, 4, 1.5, 1)); on.exit(par(op))
    for (ax in c("x", "y", "z")) {
      plot(range(d$time_s), range(d[[ax]]), type = "n",
           xlab = "time (s)", ylab = sprintf("%s (mm)", toupper(ax)))
      for (k in seq_along(ids)) {
        dd <- d[d$target_id == ids[k], ]; dd <- dd[order(dd$time_s), ]
        is_hl <- !is.null(hl) && ids[k] %in% hl
        lines(dd$time_s, dd[[ax]],
              col = if (is_hl) "#e31a1c" else pal[k],
              lwd = if (is_hl) 3 else 1)
      }
    }
  })

  output$summary <- renderPrint({
    cat("Pending decisions:\n")
    if (nrow(rv$decisions) == 0) { cat("  (none)\n"); return() }
    print(table(rv$decisions$decision))
    cat("\nPreview: tracks before =", length(unique(rv$traj$target_id)),
        "| after =", length(unique(preview()$target_id)), "\n")
  })

  # Save decisions CSV (compatible with manual_review.R)
  output$save <- downloadHandler(
    filename = function() sprintf("review_decisions_%s.csv", Sys.Date()),
    content = function(file) {
      d <- rv$decisions
      if (nrow(d) == 0) d <- data.frame(target_id = sort(unique(rv$traj$target_id)),
                                        decision = "accept", join_with = NA_integer_,
                                        join_from_frame = NA_integer_, swap_with = NA_integer_,
                                        notes = "")
      write.csv(d, file, row.names = FALSE)
    }
  )
}

`%||%` <- function(a, b) if (is.null(a)) b else a

shinyApp(ui, server)
