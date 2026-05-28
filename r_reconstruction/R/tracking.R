# tracking.R
# Tracking: constant-velocity model, particle filter, simplified GNN data association,
# clustering, MHT skeleton (simplified), and supporting functions.
#
# Mirrors MATLAB source files:
#   tracker/strXi.m, tracker/strTrk.m      (data structures)
#   tracker/predict1.m, discrete_cv.m      (motion model)
#   tracker/p_lfn.m, p_pos.m, p_mq_velocity.m  (likelihood functions)
#   tracker/pf_update.m, resample.m        (particle filter)
#   tracker/initTarget.m                   (new target init)
#   tracker/clustr.m, mhtcluster.m         (clustering)
#   tracker/munkres.m                      (Hungarian assignment)
#   tracker/hypotheses_reduction.m, nscanback.m  (MHT reduction, simplified)
#   tracker/gen_new_hyp.m, hypothesis_prob.m    (MHT -- stubs)
#   tracker/murtykbest.m                   (Murty k-best -- stub)
#   tracker/emgm.m, mqOccResolve.m         (occlusion EM -- stub)
#   tracker/update_and_initialize.m        (per-frame update)
#   tracker/terminate_track.m              (track termination)
#   tracker/mhtpftracker.m                 (main tracking loop)
#
# STATUS:
#   translated: CV motion model, particle filter (p_lfn, pf_update, resample),
#               new target init, Hungarian assignment, simplified GNN clustering,
#               per-frame update, track termination
#   partial:    N-scanback (Ns=1), cluster combine/split
#   deferred:   full MHT hypothesis enumeration (Murty), occlusion EM splitting
#
# See audit/MATLAB_TO_R_CROSSWALK.csv for full status.

# ============================================================
# STATE INDEX STRUCT (mirrors strXi.m)
# ============================================================

#' State vector indexing (mirrors strXi)
#'
#' The state vector is [x, y, z, vx, vy, vz] (6-dimensional).
#'
#' @return list with ri, rdi, nX (index vectors)
state_index <- function() {
  list(ri = 1:3, rdi = 4:6, nX = 6)
}


# ============================================================
# CONSTANT-VELOCITY MOTION MODEL (mirrors predict1.m, discrete_cv.m)
# ============================================================

#' Build discrete constant-velocity state transition matrix F
#'
#' State: [r; rdot]  6x1
#' F = [[I, dt*I], [0, I]]
#' Mirrors discrete_cv.m.
#'
#' @param dt  time step (s)
#' @return 6x6 state transition matrix
cv_motion_matrix <- function(dt) {
  F <- diag(6)
  F[1:3, 4:6] <- diag(dt, 3)
  F
}

#' Build process noise covariance Q for CV model
#'
#' Singer model: Q_rr = sw*(dt^3/3)*I, Q_rv = sw*(dt^2/2)*I, Q_vv = sw*dt*I
#' Matches the paper: w ~ N(0, sigma_w) with sigma_w = 100 m^2/s^4
#'
#' @param dt  time step (s)
#' @param sigma_w  disturbance covariance (mm^2/s^4, default 100e6 in mm units)
#' @return 6x6 process noise covariance
cv_process_noise <- function(dt, sigma_w = 100e6) {
  Q <- matrix(0, 6, 6)
  Q[1:3, 1:3] <- diag(sigma_w * dt^3 / 3, 3)
  Q[1:3, 4:6] <- diag(sigma_w * dt^2 / 2, 3)
  Q[4:6, 1:3] <- diag(sigma_w * dt^2 / 2, 3)
  Q[4:6, 4:6] <- diag(sigma_w * dt, 3)
  Q
}

#' Predict a single target state (deterministic mean, mirrors predict1.m)
#'
#' @param state  6-element state vector [x y z vx vy vz]
#' @param dt  time step
#' @return predicted state
predict_state <- function(state, dt) {
  F <- cv_motion_matrix(dt)
  as.vector(F %*% state)
}

#' Predict all particle positions forward by one time step
#'
#' Adds Gaussian disturbance to each particle.
#'
#' @param particles  6 x Np particle matrix
#' @param dt  time step (s)
#' @param sigma_w  disturbance variance (mm^2/s^4)
#' @return 6 x Np predicted particles
predict_particles <- function(particles, dt, sigma_w = 100e6) {
  F <- cv_motion_matrix(dt)
  Q <- cv_process_noise(dt, sigma_w)
  L <- tryCatch(t(chol(Q)), error = function(e) diag(sqrt(pmax(diag(Q), 0))))

  Np <- ncol(particles)
  noise <- L %*% matrix(rnorm(6 * Np), 6, Np)
  F %*% particles + noise
}

#' Predict all targets one step forward
#'
#' @param targets  list of target objects
#' @param dt  time step
#' @param sigma_w  disturbance variance
#' @return list of updated target objects
predict_targets <- function(targets, dt = 1/25, sigma_w = 100e6) {
  lapply(targets, function(t) {
    t$state_mean   <- predict_state(t$state_mean, dt)
    t$particles    <- predict_particles(t$particles, dt, sigma_w)
    t$P            <- cv_motion_matrix(dt) %*% t$P %*% t(cv_motion_matrix(dt)) +
                      cv_process_noise(dt, sigma_w)
    t
  })
}


# ============================================================
# LIKELIHOOD FUNCTIONS (mirrors p_pos.m, p_mq_velocity.m, p_lfn.m)
# ============================================================

#' Position likelihood for a particle cloud (mirrors p_pos.m, eq 3.4)
#'
#' P_mp^c(u^c | r) = N(u^c; f^c(r), Sigma_mp)
#' where f^c(r) is the perspective projection and Sigma_mp = diag(sigma_x, sigma_y)
#' comes from the blob ellipse axes.
#'
#' @param Z_meas  Z measurement struct (must have u, sigma fields)
#' @param particles  3 x Np position matrix
#' @param cam  camera struct
#' @return 1 x Np weight vector
likelihood_position <- function(Z_meas, particles, cam) {
  pix <- project_to_image(particles, cam)  # 2 x Np
  sx <- max(Z_meas$sigma[1], 0.5)
  sy <- max(Z_meas$sigma[2], 0.5)
  dnorm(pix[1,], Z_meas$u[1], sx) * dnorm(pix[2,], Z_meas$u[2], sy)
}

#' Velocity (endpoint) likelihood for a particle cloud (mirrors p_mq_velocity.m, eq 3.5)
#'
#' Bimodal because the start/end of the streak are ambiguous:
#'   P_ep^c = N(e1; f(r-), S_ep) * N(e2; f(r+), S_ep) +
#'            N(e1; f(r+), S_ep) * N(e2; f(r-), S_ep)
#' where r+/- = r +/- rdot * te / 2
#'
#' @param Z_meas  Z measurement struct (must have ep [2x2], sigma fields)
#' @param particles  6 x Np state matrix (pos rows 1:3, vel rows 4:6)
#' @param cam  camera struct
#' @param te  exposure time (s)
#' @param sigma_ep  2x2 endpoint covariance (pixels^2)
#' @return 1 x Np weight vector
likelihood_velocity <- function(Z_meas, particles, cam, te = 1/40,
                                 sigma_ep = diag(c(4, 4))) {
  r    <- particles[1:3, , drop = FALSE]
  rdot <- particles[4:6, , drop = FALSE]

  r_minus <- r - rdot * te / 2
  r_plus  <- r + rdot * te / 2

  e1_hat <- project_to_image(r_minus, cam)  # 2 x Np
  e2_hat <- project_to_image(r_plus,  cam)

  ep1 <- Z_meas$ep[, 1]  # measured start endpoint
  ep2 <- Z_meas$ep[, 2]  # measured end endpoint
  se  <- sqrt(diag(sigma_ep))  # std per dimension

  # bimodal: sum of (e1->ep1, e2->ep2) and (e1->ep2, e2->ep1)
  w1 <- dnorm(e1_hat[1,], ep1[1], se[1]) * dnorm(e1_hat[2,], ep1[2], se[2]) *
        dnorm(e2_hat[1,], ep2[1], se[1]) * dnorm(e2_hat[2,], ep2[2], se[2])

  w2 <- dnorm(e1_hat[1,], ep2[1], se[1]) * dnorm(e1_hat[2,], ep2[2], se[2]) *
        dnorm(e2_hat[1,], ep1[1], se[1]) * dnorm(e2_hat[2,], ep1[2], se[2])

  w1 + w2
}

#' Combined position + velocity likelihood across both cameras (mirrors p_lfn.m, eq 3.6)
#'
#' P(Z|X) = prod_{c=1,2} P_mp^c(u^c|r) * P_ep^c(e^c|r, rdot)
#' Also includes a uniform prior on speed (rejects stopped/hypersonic particles).
#'
#' @param Z_pair  list with Z1 and Z2 measurement structs
#' @param particles  6 x Np particle matrix
#' @param cams  list of 2 camera structs
#' @param params  parameter list from default_params()
#' @return 1 x Np weight vector
likelihood_combined <- function(Z_pair, particles, cams, params) {
  te       <- params$te
  sigma_ep <- params$sigma_ep
  v_lo     <- params$speed_min
  v_hi     <- params$speed_max

  wts <- rep(1, ncol(particles))

  for (c in 1:2) {
    Z <- if (c == 1) Z_pair$Z1 else Z_pair$Z2
    cam <- cams[[c]]
    wts <- wts *
      likelihood_position(Z, particles[1:3,, drop=FALSE], cam) *
      likelihood_velocity(Z, particles, cam, te, sigma_ep)
  }

  # Uniform speed prior (100-4000 mm/s)
  speed <- sqrt(colSums(particles[4:6,, drop=FALSE]^2))
  wts <- wts * dunif(speed, v_lo, v_hi)

  wts
}


# ============================================================
# PARTICLE FILTER (mirrors pf_update.m, resample.m)
# ============================================================

#' Systematic resampling (mirrors resample.m)
#'
#' @param weights  1 x Np numeric weight vector (need not sum to 1)
#' @return integer vector of selected particle indices
resample_particles <- function(weights) {
  Np <- length(weights)
  w  <- weights / sum(weights)
  cumw <- cumsum(w)
  u0 <- runif(1, 0, 1/Np)
  u  <- u0 + (0:(Np-1)) / Np
  idx <- rep(1L, Np)
  j <- 1L
  for (i in seq_len(Np)) {
    while (j < Np && cumw[j] < u[i]) j <- j + 1L
    idx[i] <- j
  }
  idx
}

#' Particle filter update step (mirrors pf_update.m)
#'
#' Weights particles by the combined likelihood, resamples,
#' and updates the triangulated 3D position.
#'
#' @param target  target list with particles [6 x Np], state_mean [6], P [6x6]
#' @param Z_pair  list with Z1, Z2 measurement structs (or NULL if no measurement)
#' @param cams  list of 2 camera structs
#' @param params  parameter list
#' @return updated target list
pf_update <- function(target, Z_pair, cams, params) {
  Np <- ncol(target$particles)

  if (!is.null(Z_pair)) {
    wts <- likelihood_combined(Z_pair, target$particles, cams, params)
    wts[is.nan(wts) | is.infinite(wts)] <- 0
    if (sum(wts) < 1e-300) {
      wts <- rep(1/Np, Np)
    }
    idx <- resample_particles(wts)
    target$particles <- target$particles[, idx, drop = FALSE]

    # Update 3D position via triangulation of particle-mean projection
    u_mean <- rowMeans(project_to_image(target$particles[1:3,], cams[[1]]))
    v_mean <- rowMeans(project_to_image(target$particles[1:3,], cams[[2]]))
    tri <- ls_triangulate(cbind(u_mean, v_mean), cams)
    target$state_mean[1:3] <- tri$r
    target$state_mean[4:6] <- rowMeans(target$particles[4:6,])
    target$P[1:3, 1:3]    <- cov(t(target$particles[1:3,]))
    target$P[4:6, 4:6]    <- cov(t(target$particles[4:6,]))
  }

  target
}


# ============================================================
# TARGET INITIALISATION (mirrors initTarget.m)
# ============================================================

#' Initialise a new target from a stereo measurement pair
#'
#' Creates Np particles with small position uncertainty (5 mm) and
#' large velocity uncertainty (500 mm/s), centred on the triangulated point.
#' Mirrors initTarget.m and paper section 3.3.
#'
#' @param Z_pair  list with Z1, Z2 measurement structs
#' @param cams  list of 2 camera structs
#' @param params  parameter list
#' @param target_id  integer ID for this target
#' @return new target list
init_target <- function(Z_pair, cams, params, target_id) {
  u1 <- Z_pair$Z1$u; u2 <- Z_pair$Z2$u
  tri <- ls_triangulate(cbind(u1, u2), cams)
  r0  <- tri$r

  Np  <- params$Np
  pos_std <- params$pos_init_std    # 5 mm
  vel_std <- params$vel_init_std    # 500 mm/s

  pos_particles <- matrix(rnorm(3 * Np, 0, pos_std), 3, Np) +
                    matrix(r0, 3, Np)
  vel_particles <- matrix(rnorm(3 * Np, 0, vel_std), 3, Np)
  particles <- rbind(pos_particles, vel_particles)

  P <- diag(c(rep(pos_std^2, 3), rep(vel_std^2, 3)))

  list(
    id           = target_id,
    state_mean   = c(r0, rep(0, 3)),
    particles    = particles,
    P            = P,
    frames_tracked = 1,
    confirmed    = FALSE,   # confirmed after 3 frames (paper section 3.3)
    active       = TRUE
  )
}

#' Confirm new targets that have been tracked for >= 3 frames
#'
#' @param targets  list of target objects
#' @param min_frames  frames required for confirmation (default 3)
#' @return updated targets
confirm_targets <- function(targets, min_frames = 3) {
  lapply(targets, function(t) {
    if (!t$confirmed && t$frames_tracked >= min_frames) t$confirmed <- TRUE
    t
  })
}

#' Mark a target as inactive from frame k onwards (mirrors terminate_track.m)
#'
#' @param target  target object
#' @return target with active=FALSE
terminate_track <- function(target) {
  target$active <- FALSE
  target
}

#' Remove duplicate tracks that are too close in 3D space
#'
#' Mirrors terminate_duplicate_tracks.m.
#'
#' @param targets  list of target objects
#' @param dist_threshold  mm (default 30 mm ~ 3 body lengths)
#' @return targets with duplicates removed
terminate_duplicates <- function(targets, dist_threshold = 30) {
  active <- which(vapply(targets, function(t) t$active, FALSE))
  marked_remove <- rep(FALSE, length(targets))

  for (i in seq_along(active)) {
    for (j in seq_along(active)) {
      if (i >= j) next
      ti <- targets[[active[i]]]; tj <- targets[[active[j]]]
      d <- sqrt(sum((ti$state_mean[1:3] - tj$state_mean[1:3])^2))
      if (d < dist_threshold) {
        # Remove the newer (higher ID) target
        if (ti$id > tj$id) marked_remove[active[i]] <- TRUE
        else                marked_remove[active[j]] <- TRUE
      }
    }
  }
  for (i in which(marked_remove)) targets[[i]]$active <- FALSE
  targets
}


# ============================================================
# CLUSTERING (mirrors clustr.m, gencost.m, munkres.m)
# ============================================================

#' Simplified greedy nearest-neighbour / GNN clustering
#'
#' Groups stereo measurement pairs and active targets into independent sets
#' (clusters). Each cluster contains the measurements and targets that could
#' interact based on the gating volume.
#'
#' Phase 1 simplification: uses a single round of distance-based grouping.
#' Full MHT cluster combine/split logic is deferred.
#'
#' @param Z_pairs  list of validated stereo measurement pairs
#' @param targets  list of active target objects
#' @param cams  list of camera structs
#' @param t_gate  Mahalanobis gating threshold
#' @return list of cluster objects, each with: z_ids, t_ids
cluster_measurements <- function(Z_pairs, targets, cams, t_gate = 16) {
  nz <- length(Z_pairs)
  nt <- length(targets)

  if (nz == 0 && nt == 0) return(list())

  # Build gate matrix: gate_mat[z, t] = 1 if measurement z is within gate of target t
  gate_mat <- matrix(FALSE, max(nz, 1), max(nt, 1))
  if (nz > 0 && nt > 0) {
    for (z in seq_len(nz)) {
      for (t in seq_len(nt)) {
        if (!targets[[t]]$active) next
        chi2 <- gate_check_single(
          Z_pairs[[z]]$Z1$u,
          targets[[t]]$particles[1:3, , drop = FALSE],
          cams[[1]], t_gate
        )
        gate_mat[z, t] <- chi2 < t_gate
      }
    }
  }

  # Union-Find clustering
  parent <- seq_len(nz + nt)
  find <- function(x) {
    while (parent[x] != x) { parent[x] <<- parent[parent[x]]; x <- parent[x] }
    x
  }
  union <- function(x, y) {
    px <- find(x); py <- find(y)
    if (px != py) parent[px] <<- py
  }

  for (z in seq_len(nz)) {
    for (t in seq_len(nt)) {
      if (gate_mat[z, t]) union(z, nz + t)
    }
  }

  # Collect clusters
  cluster_map <- list()
  all_ids <- seq_len(nz + nt)
  roots <- vapply(all_ids, find, 0L)

  for (root in unique(roots)) {
    members <- all_ids[roots == root]
    z_ids <- members[members <= nz]
    t_ids <- members[members > nz] - nz
    t_ids <- t_ids[t_ids >= 1 & t_ids <= nt]
    if (length(z_ids) > 0 || length(t_ids) > 0) {
      cluster_map[[length(cluster_map) + 1]] <- list(
        z_ids = z_ids,
        t_ids = t_ids
      )
    }
  }

  cluster_map
}

#' Hungarian (Munkres) algorithm for optimal assignment
#'
#' Minimises total cost given a cost matrix.
#' Uses base R optim strategy; for large matrices, package 'clue' is faster.
#'
#' @param cost_mat  nz x nt cost matrix (lower = better)
#' @return integer vector of length nt: assignment[t] = z assigned to t, or 0
hungarian_assign <- function(cost_mat) {
  nz <- nrow(cost_mat)
  nt <- ncol(cost_mat)

  if (nz == 0 || nt == 0) return(rep(0L, max(nt, 1)))
  if (nz == 1 && nt == 1) {
    return(if (is.finite(cost_mat[1,1])) 1L else 0L)
  }

  # Use base R for small matrices; fallback to greedy for larger
  if (nz <= 10 && nt <= 10) {
    return(.hungarian_small(cost_mat))
  }

  # Greedy fallback for larger matrices
  assign <- rep(0L, nt)
  used_z <- rep(FALSE, nz)
  for (t in seq_len(nt)) {
    valid <- which(!used_z & is.finite(cost_mat[, t]))
    if (length(valid) == 0) next
    best_z <- valid[which.min(cost_mat[valid, t])]
    assign[t] <- best_z
    used_z[best_z] <- TRUE
  }
  assign
}

# Small-matrix Hungarian via exhaustive permutation (n<=8)
.hungarian_small <- function(cost_mat) {
  nz <- nrow(cost_mat)
  nt <- ncol(cost_mat)
  n  <- min(nz, nt)
  best_cost <- Inf
  best_assign <- rep(0L, nt)

  perms <- .gen_perms(seq_len(nz), n)
  for (p in perms) {
    cost <- sum(vapply(seq_len(n), function(j) cost_mat[p[j], j], 0.0))
    if (cost < best_cost) {
      best_cost <- cost
      best_assign <- c(p[seq_len(n)], rep(0L, nt - n))
    }
  }
  best_assign
}

.gen_perms <- function(v, k) {
  if (k == 0) return(list(integer(0)))
  result <- list()
  for (i in seq_along(v)) {
    rest <- .gen_perms(v[-i], k - 1)
    for (r in rest) result[[length(result) + 1]] <- c(v[i], r)
  }
  result
}


# ============================================================
# SIMPLIFIED MHT STUBS (deferred from full implementation)
# ============================================================

#' STUB: Generate hypotheses for a cluster (mirrors gen_new_hyp.m)
#'
#' PHASE 1 SIMPLIFICATION: returns a single GNN assignment as a degenerate
#' hypothesis set rather than enumerating all feasible assignments.
#'
#' Full implementation would use Murty's k-best algorithm to enumerate
#' the K highest-probability assignments; see murtykbest.m.
#'
#' STATUS: deferred (see audit/UNRESOLVED_COMPONENTS.md #4)
#'
#' @param cluster  cluster object with z_ids, t_ids
#' @param Z_pairs  list of all stereo pairs
#' @param targets  list of all targets
#' @param cams  camera structs
#' @param params  parameter list
#' @return list of hypothesis objects, each with: assignment (length nt int vec), prob
generate_hypotheses <- function(cluster, Z_pairs, targets, cams, params) {
  z_ids <- cluster$z_ids
  t_ids <- cluster$t_ids
  nz <- length(z_ids); nt <- length(t_ids)

  if (nt == 0) {
    return(list(list(assignment = integer(0), prob = 1)))
  }

  # Build Euclidean cost matrix (projected distances)
  if (nz == 0) {
    # All targets missed — single null hypothesis
    return(list(list(assignment = rep(0L, nt), prob = 1 - params$P_D)))
  }

  cost_mat <- matrix(Inf, nz, nt)
  for (zi in seq_len(nz)) {
    for (ti in seq_len(nt)) {
      z_idx <- z_ids[zi]; t_idx <- t_ids[ti]
      if (!targets[[t_idx]]$active) next
      z_hat <- project_to_image(
        matrix(targets[[t_idx]]$state_mean[1:3], 3, 1), cams[[1]])
      dist <- sqrt(sum((Z_pairs[[z_idx]]$Z1$u - as.vector(z_hat))^2))
      cost_mat[zi, ti] <- dist
    }
  }

  assign <- hungarian_assign(cost_mat)
  prob <- hypothesis_prob(assign, z_ids, t_ids, Z_pairs, targets, cams, params)

  list(list(assignment = assign, prob = prob))
}

#' STUB: Compute hypothesis probability (mirrors hypothesis_prob.m)
#'
#' Simplified probability using position likelihood; full implementation
#' requires innovation covariance, P_D, clutter density per Reid (1979).
#'
#' STATUS: partial — simplified likelihood product; deferred normalisation
#'
#' @param assignment  integer vector: assign[t] = z index assigned to target t (0=missed)
#' @param z_ids, t_ids  indices into Z_pairs and targets
#' @param Z_pairs, targets, cams, params  as usual
#' @return scalar probability (unnormalised)
hypothesis_prob <- function(assignment, z_ids, t_ids, Z_pairs, targets, cams, params) {
  prob <- 1
  for (ti_local in seq_along(assignment)) {
    zi_local <- assignment[ti_local]
    t_idx    <- t_ids[ti_local]
    if (zi_local == 0) {
      prob <- prob * (1 - params$P_D)
    } else {
      z_idx <- z_ids[zi_local]
      wts <- likelihood_position(
        Z_pairs[[z_idx]]$Z1,
        matrix(targets[[t_idx]]$state_mean[1:3], 3, 1),
        cams[[1]]
      )
      prob <- prob * params$P_D * max(wts, 1e-300)
    }
  }
  prob
}

#' STUB: N-scanback hypothesis reduction (mirrors nscanback.m, Ns=1)
#'
#' Commits the highest-probability assignment from the previous step.
#' Full MHT would maintain a tree; Phase 1 commits immediately (Ns=1).
#'
#' @param hypotheses  list of hypothesis objects
#' @return single committed hypothesis
reduce_hypotheses <- function(hypotheses) {
  if (length(hypotheses) == 0) return(list(assignment = integer(0), prob = 1))
  probs <- vapply(hypotheses, function(h) h$prob, 0.0)
  hypotheses[[which.max(probs)]]
}

#' STUB: Occlusion resolution (mirrors mqOccResolve.m / emgm.m)
#'
#' PHASE 1: This function detects potential occlusions but does NOT
#' perform pixel-level EM splitting (requires image arrays).
#' Returns the original measurement list with an occlusion flag.
#'
#' Full implementation would:
#'   1. Detect when two targets share a measurement
#'   2. Use emgm.m to soft-cluster the blob pixels into individual streaks
#'   3. Return updated measurements for each occluded target
#'
#' STATUS: deferred (see audit/UNRESOLVED_COMPONENTS.md #5)
#'
#' @param assignment  committed hypothesis assignment
#' @param Z_pairs  all stereo pairs
#' @param targets  all targets
#' @return list with Z_pairs (unchanged) and occlusions (indices of detected occlusions)
resolve_occlusion <- function(assignment, Z_pairs, targets) {
  # Detect: multiple targets assigned to the same measurement
  assigned_z <- assignment[assignment > 0]
  occlusions <- which(duplicated(assigned_z) | duplicated(assigned_z, fromLast = TRUE))
  if (length(occlusions) > 0) {
    message(sprintf("[tracking.R] %d potential occlusion(s) detected. ",
                    length(unique(assigned_z[occlusions])),
                    "EM blob splitting not implemented in Phase 1 — ",
                    "occluded blobs kept as-is."))
  }
  list(Z_pairs = Z_pairs, occlusions = occlusions)
}


# ============================================================
# PER-FRAME UPDATE (mirrors update_and_initialize.m, mhtpftracker.m)
# ============================================================

#' Update all targets for one frame and initialise new ones
#'
#' Main per-frame tracking step:
#'   1. Build clusters of measurements and targets
#'   2. For each cluster: generate hypotheses, pick best, update targets
#'   3. Initialise new targets from unassigned measurement pairs
#'   4. Predict all targets for next frame
#'   5. Confirm targets tracked for >= 3 frames
#'
#' @param targets  list of current target objects
#' @param Z_pairs  validated stereo measurement pairs for this frame
#' @param cams  camera structs
#' @param params  parameter list
#' @param next_id  integer ID counter for new targets
#' @return list with: targets (updated), next_id (incremented)
update_and_initialize <- function(targets, Z_pairs, cams, params, next_id = 1L) {
  dt <- params$dt

  # Remove inactive targets
  targets <- targets[vapply(targets, function(t) t$active, FALSE)]

  # Cluster
  clusters <- cluster_measurements(Z_pairs, targets, cams, params$t_gate)

  assigned_z <- integer(0)

  for (cl in clusters) {
    if (length(cl$t_ids) == 0) next

    # Generate and reduce hypotheses
    hyps <- generate_hypotheses(cl, Z_pairs, targets, cams, params)
    best <- reduce_hypotheses(hyps)

    # Occlusion check (stub; no pixel splitting)
    occ_result <- resolve_occlusion(best$assignment, Z_pairs, targets)

    # Update each target in this cluster
    for (ti_local in seq_along(cl$t_ids)) {
      t_idx    <- cl$t_ids[ti_local]
      zi_local <- if (ti_local <= length(best$assignment)) best$assignment[ti_local] else 0

      if (zi_local > 0) {
        z_idx <- cl$z_ids[zi_local]
        assigned_z <- c(assigned_z, z_idx)
        targets[[t_idx]] <- pf_update(targets[[t_idx]], Z_pairs[[z_idx]], cams, params)
        targets[[t_idx]]$frames_tracked <- targets[[t_idx]]$frames_tracked + 1
      }
      # If no assignment: target still predicted; don't increment frame count
    }
  }

  # Initialise new targets from unassigned measurements
  unassigned <- setdiff(seq_along(Z_pairs), assigned_z)
  for (z_idx in unassigned) {
    new_t <- init_target(Z_pairs[[z_idx]], cams, params, target_id = next_id)
    targets[[length(targets) + 1]] <- new_t
    next_id <- next_id + 1L
  }

  # Confirm long-lived targets
  targets <- confirm_targets(targets, min_frames = 3)

  # Remove very close duplicates
  targets <- terminate_duplicates(targets)

  # Predict all targets for the next frame
  targets <- predict_targets(targets, dt, params$sigma_w)

  list(targets = targets, next_id = next_id)
}


# ============================================================
# XH MATRIX WRITING (mirrors t1_save_tracks.m index logic)
# ============================================================

#' Write current target states to the Xh matrix for frame k
#'
#' @param Xh  state matrix [6*Nmax x Nframes]
#' @param targets  list of target objects
#' @param frame_k  current frame index (1-based)
#' @return updated Xh
write_targets_to_xh <- function(Xh, targets, frame_k) {
  for (t in targets) {
    if (!t$active || !t$confirmed) next
    t_id <- t$id
    if (t_id > nrow(Xh) %/% 6) next  # target ID exceeds pre-allocated slots
    rows <- get_xh_row_range(t_id)
    Xh[rows, frame_k] <- t$state_mean
  }
  Xh
}
