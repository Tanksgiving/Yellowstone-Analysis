
# state-space model functions for paleots
#
# these functions extend the state-space model code to handle
# covariate-dependent ou models with one or more covariates.
# the optimum is evaluated on the trait time grid and treated
# as piecewise constant over each observed interval.
#
# reported b0 is the intercept on the centered covariate scale.
# the corresponding intercept on the raw covariate scale is stored
# in fit$misc$b0_raw.

# internal helper functions

.symmetrize <- function(M) {
  M <- as.matrix(M)
  0.5 * (M + t(M))
}

.coerce_scalar_matrix <- function(x, nrow = 1, ncol = 1) {
  matrix(x, nrow = nrow, ncol = ncol)
}

.obs_var_from_paleoTS <- function(y) {
  pmax(y$vv / y$nn, 0)
}

# put the covariates into matrix form and center each column
.prepare_covariate_matrix <- function(Z, n_expected) {
  if (is.null(dim(Z))) Z <- matrix(Z, ncol = 1)
  if (!is.matrix(Z)) stop("Z must be a matrix or coercible vector.")
  if (nrow(Z) != n_expected) stop("Number of rows in Z must match length(y$mm).")
  Zc <- scale(Z, center = TRUE, scale = FALSE)
  means <- as.numeric(attr(Zc, "scaled:center"))
  if (is.null(colnames(Zc))) {
    colnames(Zc) <- paste0("x", seq_len(ncol(Zc)))
  }
  list(
    Zc = Zc,
    means = means,
    n_cov = ncol(Zc),
    cov_names = colnames(Zc)
  )
}

.make_cov_par_names <- function(n_cov) {
  c("b0", paste0("b", seq_len(n_cov)))
}

.raw_intercept_from_centered <- function(b0_centered, slopes, means) {
  as.numeric(b0_centered - sum(slopes * means))
}

# try both signs of the hessian and return a usable variance matrix
.safe_hessian_vcov <- function(H) {
  H <- .symmetrize(H)
  out <- tryCatch(solve(H), error = function(e) NULL)
  if (!is.null(out) && all(is.finite(diag(out))) && all(diag(out) >= 0)) {
    return(out)
  }
  out <- tryCatch(solve(-H), error = function(e) NULL)
  if (!is.null(out) && all(is.finite(diag(out))) && all(diag(out) >= 0)) {
    return(out)
  }
  matrix(NA_real_, nrow(H), ncol(H), dimnames = dimnames(H))
}

# regularize the innovation covariance if it is nearly singular
.safe_spd_solve <- function(S, ridge_start = 1e-10, ridge_mult = 10, max_tries = 8) {
  S <- .symmetrize(S)
  q <- nrow(S)
  Iq <- diag(1, q)
  ridge <- 0
  for (i in 0:max_tries) {
    Sc <- if (ridge == 0) S else S + ridge * Iq
    ch <- tryCatch(chol(Sc), error = function(e) NULL)
    if (!is.null(ch)) {
      Sinv <- chol2inv(ch)
      logdet <- 2 * sum(log(diag(ch)))
      return(list(S = Sc, Sinv = Sinv, logdet = logdet, ridge = ridge))
    }
    ridge <- if (ridge == 0) ridge_start else ridge * ridge_mult
  }
  stop("Innovation covariance could not be regularized to SPD form.")
}

.build_kf_object <- function(n, qdim, pdim) {
  list(
    xp    = array(NA_real_, dim = c(pdim, 1, n)),
    Pp    = array(NA_real_, dim = c(pdim, pdim, n)),
    xf    = array(NA_real_, dim = c(pdim, 1, n)),
    Pf    = array(NA_real_, dim = c(pdim, pdim, n)),
    innov = array(NA_real_, dim = c(qdim, 1, n)),
    sig   = array(NA_real_, dim = c(qdim, qdim, n)),
    Kn    = array(NA_real_, dim = c(pdim, qdim, n))
  )
}

# convert the optimization result into the standard paleots fit object
.finalize_simple_fit <- function(fit, y, model_name, logL_fun, K, hess = FALSE, ...) {
  se <- NULL
  if (hess && !is.null(fit$hessian)) {
    # the hessian is on the optimization scale, so map it back to the reported
    # parameterization before returning standard errors
    V <- .safe_hessian_vcov(fit$hessian)
    se <- sqrt(diag(V))
  }
  w <- as.paleoTSfit(
    logL = fit$value,
    parameters = fit$par,
    modelName = model_name,
    method = "SSM",
    K = K,
    n = length(y$mm),
    se = se,
    convergence = fit$convergence,
    logLFunction = logL_fun
  )
  w
}

# functions for the simple state-space models

opt.ssm.GRW <- function(y, pool = TRUE, cl = list(fnscale = -1), hess = FALSE) {
  if (pool) y <- pool.var(y, ret.paleoTS = TRUE)
  if (y$tt[1] != 0) {
    stop("Initial time should be 0. Use as.paleoTS() or read.paleoTS() to correctly process ages.")
  }
  
  # use the observed first value and the mle from the standard paleots fit as starting values
  p0 <- c(y$mm[1], mle.GRW(y))
  if (p0[3] <= 0) p0[3] <- 1e-5
  names(p0) <- c("anc", "mstep", "vstep")
  
  if (is.null(cl$parscale)) cl$parscale <- getAllParscale(y, p0)
  
  fit <- optim(
    p0, fn = logL.ssm.GRW, method = "L-BFGS-B", control = cl,
    lower = c(NA, NA, 1e-6), hessian = hess, y = y
  )
  
  w <- .finalize_simple_fit(
    fit = fit, y = y, model_name = "GRW",
    logL_fun = "logL.ssm.GRW", K = 3, hess = hess
  )
  w$kf <- logL.ssm.GRW(p = w$parameters, y = y, logL.only = FALSE)
  w
}

opt.ssm.URW <- function(y, pool = TRUE, cl = list(fnscale = -1), hess = FALSE) {
  if (pool) y <- pool.var(y, ret.paleoTS = TRUE)
  if (y$tt[1] != 0) {
    stop("Initial time should be 0. Use as.paleoTS() or read.paleoTS() to correctly process ages.")
  }
  
  # use the observed first value and the mle from the standard paleots fit as starting values
  p0 <- c(y$mm[1], mle.URW(y))
  if (p0[2] <= 0) p0[2] <- 1e-5
  names(p0) <- c("anc", "vstep")
  
  if (is.null(cl$parscale)) cl$parscale <- getAllParscale(y, p0)
  
  fit <- optim(
    p0, fn = logL.ssm.URW, method = "L-BFGS-B", control = cl,
    lower = c(NA, 1e-6), hessian = hess, y = y
  )
  
  w <- .finalize_simple_fit(
    fit = fit, y = y, model_name = "URW",
    logL_fun = "logL.ssm.URW", K = 2, hess = hess
  )
  w$kf <- logL.ssm.URW(p = w$parameters, y = y, logL.only = FALSE)
  w
}

opt.ssm.Stasis <- function(y, pool = TRUE, cl = list(fnscale = -1), hess = FALSE) {
  if (pool) y <- pool.var(y, ret.paleoTS = TRUE)
  if (y$tt[1] != 0) {
    stop("Initial time should be 0. Use as.paleoTS() or read.paleoTS() to correctly process ages.")
  }
  
  # take the starting values from the standard stasis fit
  p0 <- mle.Stasis(y)
  if (p0[2] <= 0) p0[2] <- 1e-5
  names(p0) <- c("theta", "omega")
  
  if (is.null(cl$parscale)) cl$parscale <- getAllParscale(y, p0)
  
  fit <- optim(
    p0, fn = logL.ssm.Stasis, method = "L-BFGS-B", control = cl,
    lower = c(NA, 1e-6), hessian = hess, y = y
  )
  
  w <- .finalize_simple_fit(
    fit = fit, y = y, model_name = "Stasis",
    logL_fun = "logL.ssm.Stasis", K = 2, hess = hess
  )
  w$kf <- logL.ssm.Stasis(p = w$parameters, y = y, logL.only = FALSE)
  w
}

opt.ssm.StrictStasis <- function(y, pool = TRUE, cl = list(fnscale = -1), hess = FALSE) {
  if (pool) y <- pool.var(y, ret.paleoTS = TRUE)
  if (y$tt[1] != 0) {
    stop("Initial time should be 0. Use as.paleoTS() or read.paleoTS() to correctly process ages.")
  }
  
  # start strict stasis at the overall mean
  p0 <- mean(y$mm)
  names(p0) <- "theta"
  
  if (is.null(cl$parscale)) cl$parscale <- getAllParscale(y, p0)
  
  fit <- optim(
    p0, fn = logL.ssm.StrictStasis, method = "L-BFGS-B",
    control = cl, lower = c(NA), hessian = hess, y = y
  )
  
  w <- .finalize_simple_fit(
    fit = fit, y = y, model_name = "StrictStasis",
    logL_fun = "logL.ssm.StrictStasis", K = 1, hess = hess
  )
  w$kf <- logL.ssm.StrictStasis(p = w$parameters, y = y, logL.only = FALSE)
  w
}

opt.ssm.OU <- function(y, pool = TRUE, cl = list(fnscale = -1), hess = FALSE) {
  if (pool) y <- pool.var(y, ret.paleoTS = TRUE)
  if (y$tt[1] != 0) {
    stop("Initial time should be 0. Use as.paleoTS() or read.paleoTS() to correctly process ages.")
  }
  
  # start the ou fit from the urw variance, the mean trait value, and a moderate half-life
  p0 <- c(y$mm[1], mle.URW(y), mean(y$mm), log(2) / (0.5 * max(y$tt)))
  if (p0[2] <= 0) p0[2] <- 1e-5
  names(p0) <- c("anc", "vstep", "theta", "alpha")
  
  if (is.null(cl$parscale)) cl$parscale <- getAllParscale(y, p0)
  
  fit <- optim(
    p0, fn = logL.ssm.OU, method = "L-BFGS-B", control = cl,
    lower = c(NA, 1e-6, NA, 1e-6), hessian = hess, y = y
  )
  
  w <- .finalize_simple_fit(
    fit = fit, y = y, model_name = "OU",
    logL_fun = "logL.ssm.OU", K = 4, hess = hess
  )
  w$kf <- logL.ssm.OU(p = w$parameters, y = y, logL.only = FALSE)
  w
}

opt.ssm.ACDC <- function(y, pool = TRUE, cl = list(fnscale = -1), hess = FALSE) {
  if (pool) y <- pool.var(y, ret.paleoTS = TRUE)
  if (y$tt[1] != 0) {
    stop("Initial time should be 0. Use as.paleoTS() or read.paleoTS() to correctly process ages.")
  }
  
  # start acdc at the urw variance with no acceleration or deceleration
  p0 <- c(y$mm[1], mle.URW(y), 0)
  if (p0[2] <= 0) p0[2] <- 1e-3
  names(p0) <- c("anc", "vstep0", "r")
  
  if (is.null(cl$parscale)) cl$parscale <- getAllParscale(y, p0)
  
  fit <- optim(
    p0, fn = logL.ssm.ACDC, method = "L-BFGS-B", control = cl,
    lower = c(NA, 1e-6, NA), hessian = hess, y = y
  )
  
  w <- .finalize_simple_fit(
    fit = fit, y = y, model_name = "ACDC",
    logL_fun = "logL.ssm.ACDC", K = 3, hess = hess
  )
  w$kf <- logL.ssm.ACDC(p = w$parameters, y = y, logL.only = FALSE)
  w
}

# functions for covariate-dependent ou models
#
# these functions extend the usual ou model by letting the optimum move with
# one or more measured covariates. the latent process is still ou, but theta
# is no longer fixed through the whole sequence.

opt.ssm.covOU <- function(y, Z, pool = TRUE, cl = list(fnscale = -1), hess = FALSE) {
  # fit an ou model where the optimum depends on one or more covariates
  # the optimum is evaluated on the observed trait time grid and treated as
  # piecewise constant over each observed interval
  if (pool) y <- pool.var(y, ret.paleoTS = TRUE)
  if (y$tt[1] != 0) {
    stop("Initial time should be 0. Use as.paleoTS() or read.paleoTS().")
  }
  
  # center the covariates here so the slopes describe departures from the
  # average environment and b0 stays interpretable as the optimum at the
  # mean covariate values
  # the covariates are centered in the same way as covOU so the optimum
  # coefficients have the same interpretation in the variance-shift version
  prep <- .prepare_covariate_matrix(Z, n_expected = length(y$mm))
  Zc <- prep$Zc
  raw_means <- prep$means
  n_cov <- prep$n_cov
  
  dt_mean <- mean(diff(y$tt))
  alpha_start <- max(log(2) / (5 * dt_mean), 1e-8)
  vstep_init <- max(1e-6, mle.URW(y))
  
  # use a simple linear model to get reasonable starting values for the
  # optimum function. this does not define the final fit, but it helps the
  # optimizer start in a reasonable part of parameter space
  wr <- stats::lm(y$mm ~ Zc)
  reg_coefs <- unname(coef(wr))
  names(reg_coefs) <- .make_cov_par_names(n_cov)
  
  p0 <- c(
    anc = y$mm[1],
    log_vstep = log(vstep_init),
    log_alpha = log(alpha_start),
    reg_coefs
  )
  
  if (is.null(cl$parscale)) {
    cl$parscale <- tryCatch(
      getAllParscale(y, p0, Zc),
      error = function(e) rep(1, length(p0))
    )
  }
  
  fit <- optim(
    par = p0,
    fn = logL.ssm.covOU.multi,
    method = "L-BFGS-B",
    lower = rep(-Inf, length(p0)),
    control = cl,
    hessian = hess,
    y = y,
    Z = Zc
  )
  
  opt_par <- fit$par
  names(opt_par) <- names(p0)
  
  slopes <- opt_par[paste0("b", seq_len(n_cov))]
  b0_centered <- opt_par["b0"]
  b0_raw <- .raw_intercept_from_centered(b0_centered, slopes, raw_means)
  
  # report vstep and alpha back on their natural scales
  # keep b0 on the centered scale because that is the quantity tied to the
  # average covariate conditions used in the fitted model
  rep_par <- c(
    anc = opt_par["anc"],
    vstep = exp(opt_par["log_vstep"]),
    alpha = exp(opt_par["log_alpha"]),
    b0 = b0_centered,
    slopes
  )
  names(rep_par) <- c("anc", "vstep", "alpha", "b0", paste0("b", seq_len(n_cov)))
  
  w <- as.paleoTSfit(
    logL = fit$value,
    parameters = rep_par,
    modelName = "covOU-multi",
    method = "SSM",
    K = length(p0),
    n = length(y$mm),
    se = NULL,
    convergence = fit$convergence,
    logLFunction = "logL.ssm.covOU.multi"
  )
  
  # save the kalman filter output from the optimized parameter vector so the
  # fitted states, innovations, and residual diagnostics are available later
  w$kf <- logL.ssm.covOU.multi(p = opt_par, y = y, Z = Zc, logL.only = FALSE)
  w$misc <- list(
    Z_means = raw_means,
    Z_centered = TRUE,
    b0_centered = unname(b0_centered),
    b0_raw = unname(b0_raw),
    opt_par = opt_par
  )
  
  if (hess && !is.null(fit$hessian)) {
    # the hessian is on the optimization scale, so map it back to the reported
    # parameterization before returning standard errors
    Vopt <- .safe_hessian_vcov(fit$hessian)
    dimnames(Vopt) <- list(names(p0), names(p0))
    
    rn <- names(rep_par)
    J <- matrix(0, nrow = length(rn), ncol = length(p0), dimnames = list(rn, names(p0)))
    J["anc", "anc"] <- 1
    J["vstep", "log_vstep"] <- rep_par["vstep"]
    J["alpha", "log_alpha"] <- rep_par["alpha"]
    J["b0", "b0"] <- 1
    for (i in seq_len(n_cov)) {
      bi <- paste0("b", i)
      J[bi, bi] <- 1
    }
    
    Vrep <- J %*% Vopt %*% t(J)
    w$se <- sqrt(diag(Vrep))
    names(w$se) <- rn
    w$vcov <- Vrep
  }
  
  w
}

opt.ssm.covOU_vshift <- function(y, Z, gg, pool = TRUE, cl = list(fnscale = -1), hess = FALSE) {
  # fit a covariate-dependent ou model with process variance allowed to shift among segments
  # this keeps the same moving optimum as covOU, but lets stochastic variance
  # differ among predefined parts of the sequence
  if (pool) y <- pool.var(y, ret.paleoTS = TRUE)
  if (y$tt[1] != 0) {
    stop("Initial time should be 0. Use as.paleoTS() or read.paleoTS().")
  }
  if (length(gg) != length(y$mm)) stop("gg must have length equal to length(y$mm).")
  if (any(!is.finite(gg))) stop("gg must be finite.")
  gg <- as.integer(gg)
  nseg <- max(gg)
  if (!all(gg %in% seq_len(nseg))) stop("gg must be coded as integers in 1:max(gg).")
  
  prep <- .prepare_covariate_matrix(Z, n_expected = length(y$mm))
  Zc <- prep$Zc
  raw_means <- prep$means
  n_cov <- prep$n_cov
  
  dt_mean <- mean(diff(y$tt))
  alpha_start <- max(log(2) / (5 * dt_mean), 1e-8)
  vstep_init <- max(1e-6, mle.URW(y))
  
  wr <- stats::lm(y$mm ~ Zc)
  reg_coefs <- unname(coef(wr))
  names(reg_coefs) <- .make_cov_par_names(n_cov)
  
  p0 <- c(
    anc = y$mm[1],
    setNames(rep(log(vstep_init), nseg), paste0("log_vstep", seq_len(nseg))),
    log_alpha = log(alpha_start),
    reg_coefs
  )
  
  if (is.null(cl$parscale)) {
    cl$parscale <- tryCatch(
      getAllParscale(y, p0, Zc),
      error = function(e) rep(1, length(p0))
    )
  }
  
  fit <- optim(
    par = p0,
    fn = logL.ssm.covOU_vshift.multi,
    method = "L-BFGS-B",
    lower = rep(-Inf, length(p0)),
    control = cl,
    hessian = hess,
    y = y,
    Z = Zc,
    gg = gg
  )
  
  opt_par <- fit$par
  names(opt_par) <- names(p0)
  
  slopes <- opt_par[paste0("b", seq_len(n_cov))]
  b0_centered <- opt_par["b0"]
  b0_raw <- .raw_intercept_from_centered(b0_centered, slopes, raw_means)
  vsteps <- exp(opt_par[paste0("log_vstep", seq_len(nseg))])
  
  rep_par <- c(
    anc = opt_par["anc"],
    setNames(vsteps, paste0("vstep", seq_len(nseg))),
    alpha = exp(opt_par["log_alpha"]),
    b0 = b0_centered,
    slopes
  )
  names(rep_par) <- c("anc", paste0("vstep", seq_len(nseg)), "alpha", "b0", paste0("b", seq_len(n_cov)))
  
  w <- as.paleoTSfit(
    logL = fit$value,
    parameters = rep_par,
    modelName = "covOU_vshift-multi",
    method = "SSM",
    K = length(p0),
    n = length(y$mm),
    se = NULL,
    convergence = fit$convergence,
    logLFunction = "logL.ssm.covOU_vshift.multi"
  )
  
  # save the kalman filter output from the optimized fit so the same residual
  # and adequacy diagnostics can be used for the variance-shift model
  w$kf <- logL.ssm.covOU_vshift.multi(p = opt_par, y = y, Z = Zc, gg = gg, logL.only = FALSE)
  w$gg <- gg
  w$misc <- list(
    Z_means = raw_means,
    Z_centered = TRUE,
    b0_centered = unname(b0_centered),
    b0_raw = unname(b0_raw),
    opt_par = opt_par
  )
  
  if (hess && !is.null(fit$hessian)) {
    # the hessian is on the optimization scale, so map it back to the reported
    # parameterization before returning standard errors
    Vopt <- .safe_hessian_vcov(fit$hessian)
    dimnames(Vopt) <- list(names(p0), names(p0))
    
    rn <- names(rep_par)
    J <- matrix(0, nrow = length(rn), ncol = length(p0), dimnames = list(rn, names(p0)))
    J["anc", "anc"] <- 1
    for (i in seq_len(nseg)) {
      J[paste0("vstep", i), paste0("log_vstep", i)] <- rep_par[paste0("vstep", i)]
    }
    J["alpha", "log_alpha"] <- rep_par["alpha"]
    J["b0", "b0"] <- 1
    for (i in seq_len(n_cov)) {
      bi <- paste0("b", i)
      J[bi, bi] <- 1
    }
    
    Vrep <- J %*% Vopt %*% t(J)
    w$se <- sqrt(diag(Vrep))
    names(w$se) <- rn
    w$vcov <- Vrep
  }
  
  w
}

# single-covariate wrappers
#
# these wrappers keep the older calling pattern available while
# sending the work to the revised multi-covariate functions. this lets
# older code keep working without duplicating the newer covariate-ou code.

opt.ssm.covOU.single <- function(y, z, pool = TRUE, cl = list(fnscale = -1), hess = FALSE) {
  opt.ssm.covOU(y = y, Z = matrix(z, ncol = 1), pool = pool, cl = cl, hess = hess)
}

opt.ssm.covOU_vshift.single <- function(y, z, gg, pool = TRUE, cl = list(fnscale = -1), hess = FALSE) {
  opt.ssm.covOU_vshift(y = y, Z = matrix(z, ncol = 1), gg = gg, pool = pool, cl = cl, hess = hess)
}

# log-likelihood functions
#
# each function below builds the time-varying state-space system
# and returns either the log likelihood or the full kalman output.

logL.ssm.URW <- function(p, y, logL.only = TRUE) {
  # random walk with sampling error separated in the observation equation
  anc <- p[1]
  vs <- p[2]
  n <- length(y$mm)
  dt <- c(0, diff(y$tt))
  
  Atv <- array(1, dim = c(1, 1, n))
  Phitv <- array(1, dim = c(1, 1, n))
  mu0 <- anc
  Sigma0 <- .obs_var_from_paleoTS(y)[1]
  Rtv <- array(.obs_var_from_paleoTS(y), dim = c(1, 1, n))
  Qtv <- array(vs * dt, dim = c(1, 1, n))
  
  kf <- Kfiltertv(
    num = n, y = y$mm, Atv = Atv, mu0 = mu0, Sigma0 = Sigma0,
    Phitv = Phitv, Ups = NULL, Gam = NULL, Qtv = Qtv, Rtv = Rtv, input = NULL
  )
  if (logL.only) kf$like else kf
}

logL.ssm.GRW <- function(p, y, logL.only = TRUE) {
  # generalized random walk with a directional step component
  anc <- p[1]
  ms <- p[2]
  vs <- p[3]
  n <- length(y$mm)
  dt <- c(0, diff(y$tt))
  
  Atv <- array(1, dim = c(1, 1, n))
  Phitv <- array(1, dim = c(1, 1, n))
  mu0 <- anc
  Sigma0 <- .obs_var_from_paleoTS(y)[1]
  Rtv <- array(.obs_var_from_paleoTS(y), dim = c(1, 1, n))
  Qtv <- array(vs * dt, dim = c(1, 1, n))
  input <- ms * dt
  
  kf <- Kfiltertv(
    num = n, y = y$mm, Atv = Atv, mu0 = mu0, Sigma0 = Sigma0,
    Phitv = Phitv, Ups = 1, Gam = NULL, Qtv = Qtv, Rtv = Rtv, input = input
  )
  if (logL.only) kf$like else kf
}

logL.ssm.Stasis <- function(p, y, logL.only = TRUE) {
  # stasis represented as white-noise deviations around a fixed mean
  theta <- p[1]
  omega <- p[2]
  n <- length(y$mm)
  
  Atv <- array(1, dim = c(1, 1, n))
  Phitv <- array(0, dim = c(1, 1, n))
  mu0 <- y$mm[1]
  Sigma0 <- .obs_var_from_paleoTS(y)[1]
  Rtv <- array(.obs_var_from_paleoTS(y), dim = c(1, 1, n))
  Qtv <- array(omega, dim = c(1, 1, n))
  input <- rep(1, n)
  
  kf <- Kfiltertv(
    num = n, y = y$mm, Atv = Atv, mu0 = mu0, Sigma0 = Sigma0,
    Phitv = Phitv, Ups = theta, Gam = NULL, Qtv = Qtv, Rtv = Rtv, input = input
  )
  if (logL.only) kf$like else kf
}

logL.ssm.StrictStasis <- function(p, y, logL.only = TRUE) {
  # strict stasis with no process variance
  theta <- p[1]
  n <- length(y$mm)
  
  Atv <- array(1, dim = c(1, 1, n))
  Phitv <- array(0, dim = c(1, 1, n))
  mu0 <- y$mm[1]
  Sigma0 <- .obs_var_from_paleoTS(y)[1]
  Rtv <- array(.obs_var_from_paleoTS(y), dim = c(1, 1, n))
  Qtv <- array(0, dim = c(1, 1, n))
  input <- rep(1, n)
  
  kf <- Kfiltertv(
    num = n, y = y$mm, Atv = Atv, mu0 = mu0, Sigma0 = Sigma0,
    Phitv = Phitv, Ups = theta, Gam = NULL, Qtv = Qtv, Rtv = Rtv, input = input
  )
  if (logL.only) kf$like else kf
}

logL.ssm.OU <- function(p, y, logL.only = TRUE) {
  # ou model with a fixed optimum
  anc <- p[1]
  vs <- p[2]
  theta <- p[3]
  alpha <- p[4]
  
  n <- length(y$mm)
  dt <- c(0, diff(y$tt))
  
  Atv <- array(1, dim = c(1, 1, n))
  e1 <- exp(-alpha * dt)
  e1m1 <- -expm1(-alpha * dt)
  e2m1 <- -expm1(-2 * alpha * dt)
  
  Phitv <- array(e1, c(1, 1, n))
  Qtv <- array((vs / (2 * alpha)) * e2m1, c(1, 1, n))
  mu0 <- anc
  Sigma0 <- .obs_var_from_paleoTS(y)[1]
  Rtv <- array(.obs_var_from_paleoTS(y), c(1, 1, n))
  input <- theta * e1m1
  
  kf <- Kfiltertv(
    num = n, y = y$mm, Atv = Atv, mu0 = mu0, Sigma0 = Sigma0,
    Phitv = Phitv, Ups = 1, Gam = NULL, Qtv = Qtv, Rtv = Rtv, input = input
  )
  if (logL.only) kf$like else kf
}

logL.ssm.ACDC <- function(p, y, logL.only = TRUE) {
  # random walk with step variance changing exponentially through time
  anc <- p[1]
  vs <- p[2]
  r <- p[3]
  n <- length(y$mm)
  dt <- c(0, diff(y$tt))
  si <- y$tt - y$tt[1]
  
  Atv <- array(1, dim = c(1, 1, n))
  Phitv <- array(1, dim = c(1, 1, n))
  mu0 <- anc
  Sigma0 <- .obs_var_from_paleoTS(y)[1]
  Rtv <- array(.obs_var_from_paleoTS(y), dim = c(1, 1, n))
  
  if (abs(r) < 1e-10) {
    Qtv <- array(vs * dt, dim = c(1, 1, n))
  } else {
    V <- vs * (expm1(r * si) / r)
    Qtv <- array(c(0, diff(V)), dim = c(1, 1, n))
  }
  
  kf <- Kfiltertv(
    num = n, y = y$mm, Atv = Atv, mu0 = mu0, Sigma0 = Sigma0,
    Phitv = Phitv, Ups = NULL, Gam = NULL, Qtv = Qtv, Rtv = Rtv, input = NULL
  )
  if (logL.only) kf$like else kf
}

logL.ssm.covOU.multi <- function(p, y, Z, logL.only = TRUE) {
  # ou model with a covariate-defined optimum
  # theta_t is built from the centered covariate matrix and then inserted into
  # the discrete-time ou transition used by the kalman filter
  anc <- p[1]
  vs <- exp(p[2])
  alpha <- exp(p[3])
  
  n_cov <- ncol(Z)
  b_vec <- p[4:(3 + n_cov + 1)]
  b0 <- b_vec[1]
  slopes <- b_vec[-1]
  
  n <- length(y$mm)
  dt <- c(0, diff(y$tt))
  
  Atv <- array(1, dim = c(1, 1, n))
  e1 <- exp(-alpha * dt)
  e1m1 <- -expm1(-alpha * dt)
  e2m1 <- -expm1(-2 * alpha * dt)
  
  Phitv <- array(e1, c(1, 1, n))
  Qtv <- array((vs / (2 * alpha)) * e2m1, c(1, 1, n))
  mu0 <- anc
  Sigma0 <- .obs_var_from_paleoTS(y)[1]
  Rtv <- array(.obs_var_from_paleoTS(y), c(1, 1, n))
  
  # build the optimum at each sampled level from the centered covariates
  # build the optimum at each sampled level from the centered covariates
  theta_t <- b0 + as.vector(Z %*% slopes)
  input <- theta_t * e1m1
  
  kf <- Kfiltertv(
    num = n, y = y$mm, Atv = Atv, mu0 = mu0, Sigma0 = Sigma0,
    Phitv = Phitv, Ups = 1, Gam = NULL, Qtv = Qtv, Rtv = Rtv, input = input
  )
  if (logL.only) kf$like else kf
}

logL.ssm.covOU_vshift.multi <- function(p, y, Z, gg, logL.only = TRUE) {
  # covariate-dependent ou model with variance shifts across predefined segments
  # the optimum is still covariate-defined, but q_t now changes with the
  # segment membership vector gg
  nseg <- max(gg)
  anc <- p[1]
  vsteps <- exp(p[2:(1 + nseg)])
  alpha <- exp(p[2 + nseg])
  
  n_cov <- ncol(Z)
  b_vec <- p[(nseg + 3):(nseg + 3 + n_cov)]
  b0 <- b_vec[1]
  slopes <- b_vec[-1]
  
  n <- length(y$mm)
  dt <- c(0, diff(y$tt))
  vs_samp <- vsteps[gg]
  
  Atv <- array(1, dim = c(1, 1, n))
  e1 <- exp(-alpha * dt)
  e1m1 <- -expm1(-alpha * dt)
  e2m1 <- -expm1(-2 * alpha * dt)
  
  Phitv <- array(e1, c(1, 1, n))
  Qtv <- array((vs_samp / (2 * alpha)) * e2m1, c(1, 1, n))
  mu0 <- anc
  Sigma0 <- .obs_var_from_paleoTS(y)[1]
  Rtv <- array(.obs_var_from_paleoTS(y), c(1, 1, n))
  
  # build the optimum at each sampled level from the centered covariates
  theta_t <- b0 + as.vector(Z %*% slopes)
  input <- theta_t * e1m1
  
  kf <- Kfiltertv(
    num = n, y = y$mm, Atv = Atv, mu0 = mu0, Sigma0 = Sigma0,
    Phitv = Phitv, Ups = 1, Gam = NULL, Qtv = Qtv, Rtv = Rtv, input = input
  )
  if (logL.only) kf$like else kf
}

# wrappers for the older covariate-ou parameterization
# the older natural-scale parameter order expected here is
# covou: c(anc, vstep, b0, b1, alpha)
# covou_vshift: c(anc, vstep1, ..., vstepN, b0, b1, alpha)

logL.ssm.covOU <- function(p, y, z, logL.only = TRUE) {
  Z <- matrix(z, ncol = 1)
  if (length(p) == 5 && !any(grepl("^log_", names(p)))) {
    p_opt <- c(
      anc = p[1],
      log_vstep = log(p[2]),
      log_alpha = log(p[5]),
      b0 = p[3],
      b1 = p[4]
    )
  } else {
    p_opt <- p
  }
  logL.ssm.covOU.multi(p = p_opt, y = y, Z = Z, logL.only = logL.only)
}

logL.ssm.covOU_vshift <- function(p, y, z, gg, logL.only = TRUE) {
  Z <- matrix(z, ncol = 1)
  nseg <- max(gg)
  if (length(p) == nseg + 4 && !any(grepl("^log_", names(p)))) {
    p_opt <- c(
      anc = p[1],
      setNames(log(p[2:(nseg + 1)]), paste0("log_vstep", seq_len(nseg))),
      log_alpha = log(p[nseg + 4]),
      b0 = p[nseg + 2],
      b1 = p[nseg + 3]
    )
  } else {
    p_opt <- p
  }
  logL.ssm.covOU_vshift.multi(p = p_opt, y = y, Z = Z, gg = gg, logL.only = logL.only)
}

logL.ssm.URWshift <- function(p, y, gg, logL.only = TRUE) {
  nseg <- max(gg)
  anc <- p[1]
  vsv <- p[2:(nseg + 1)]
  vss <- vsv[gg]
  
  n <- length(y$mm)
  dt <- c(0, diff(y$tt))
  
  Atv <- array(1, dim = c(1, 1, n))
  Phitv <- array(1, dim = c(1, 1, n))
  mu0 <- anc
  Sigma0 <- .obs_var_from_paleoTS(y)[1]
  Rtv <- array(.obs_var_from_paleoTS(y), dim = c(1, 1, n))
  Qtv <- array(vss * dt, dim = c(1, 1, n))
  
  kf <- Kfiltertv(
    num = n, y = y$mm, Atv = Atv, mu0 = mu0, Sigma0 = Sigma0,
    Phitv = Phitv, Ups = NULL, Gam = NULL, Qtv = Qtv, Rtv = Rtv, input = NULL
  )
  if (logL.only) kf$like else kf
}

# random walk with variance shifts

opt.ssm.URWshift <- function(y, gg, pool = TRUE, cl = list(fnscale = -1), hess = FALSE) {
  # fit a random walk where the step variance can differ among time segments
  if (pool) y <- pool.var(y, ret.paleoTS = TRUE)
  if (y$tt[1] != 0) {
    stop("Initial time should be 0. Use as.paleoTS() or read.paleoTS() to correctly process ages.")
  }
  
  if (length(gg) != length(y$mm)) stop("gg must have length n.")
  nseg <- max(gg)
  if (!all(gg %in% seq_len(nseg))) stop("gg must be integers in 1..max(gg).")
  
  p0 <- c(y$mm[1], rep(mle.URW(y), nseg))
  names(p0) <- c("anc", paste0("vstep", seq_len(nseg)))
  
  if (is.null(cl$parscale)) cl$parscale <- getAllParscale(y, p0)
  
  fit <- optim(
    p0, fn = logL.ssm.URWshift, method = "L-BFGS-B",
    lower = c(NA, rep(1e-6, nseg)), hessian = hess,
    y = y, gg = gg, control = cl
  )
  
  w <- .finalize_simple_fit(
    fit = fit, y = y, model_name = "URW-shift",
    logL_fun = "logL.ssm.URWshift", K = length(p0), hess = hess
  )
  w$kf <- logL.ssm.URWshift(p = w$parameters, y = y, gg = gg, logL.only = FALSE)
  w$gg <- gg
  w
}

# residual diagnostics
#
# by default these use standardized innovations from the kalman filter.

checkSSMresiduals <- function(y, w, show.plot = TRUE,
                              resid.type = c("standardized", "unstandardized"),
                              use.innov = TRUE,
                              ask = interactive()) {
  resid.type <- match.arg(resid.type)
  
  have_innov <- !is.null(w$kf$innov) && !is.null(w$kf$sig)
  
  if (isTRUE(use.innov) && have_innov) {
    rr <- as.vector(w$kf$innov)
    scalev <- sqrt(drop(w$kf$sig))
  } else {
    pred <- drop(w$kf$xp)
    rr <- y$mm - pred
    if (have_innov) {
      scalev <- sqrt(drop(w$kf$sig))
    } else {
      scalev <- sqrt(pmax(drop(w$kf$Pp), 0))
    }
  }
  
  scalev[!is.finite(scalev) | scalev <= 0] <- NA_real_
  if (resid.type == "standardized") rr <- rr / scalev
  
  if (show.plot) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar), add = TRUE)
    par(ask = ask)
    plot(rr, xlab = "Index", ylab = paste(resid.type, "residuals"), main = "SSM Residuals")
    abline(h = 0, lty = 3)
    stats::qqnorm(rr, main = paste("QQ Plot of", resid.type, "residuals"))
    stats::qqline(rr)
  }
  
  invisible(rr)
}

# kalman filter
#
# this version keeps the original structure but uses a more stable
# covariance update and regularizes the innovation covariance when needed.

Kfiltertv <- function(num, y, Atv, mu0, Sigma0, Phitv, Ups, Gam, Qtv, Rtv, input) {
  y <- as.matrix(y)
  qdim <- ncol(y)
  
  Phi1 <- as.matrix(Phitv[, , 1])
  pdim <- nrow(Phi1)
  
  mu0 <- matrix(mu0, nrow = pdim, ncol = 1)
  Sigma0 <- matrix(Sigma0, nrow = pdim, ncol = pdim)
  
  have_input <- !is.null(input)
  if (have_input) {
    input <- as.matrix(input)
    if (nrow(input) != num) {
      input <- matrix(input, nrow = num, byrow = FALSE)
    }
    rdim <- ncol(input)
    if (is.null(Ups)) Ups <- matrix(0, nrow = pdim, ncol = rdim)
    if (is.null(Gam)) Gam <- matrix(0, nrow = qdim, ncol = rdim)
    Ups <- as.matrix(Ups)
    Gam <- as.matrix(Gam)
  }
  
  out <- .build_kf_object(n = num, qdim = qdim, pdim = pdim)
  like_accum <- 0
  
  for (i in seq_len(num)) {
    Phi <- as.matrix(Phitv[, , i])
    Q <- .symmetrize(as.matrix(Qtv[, , i]))
    R <- .symmetrize(as.matrix(Rtv[, , i]))
    B <- matrix(Atv[, , i], nrow = qdim, ncol = pdim)
    
    y_i <- matrix(y[i, ], nrow = qdim, ncol = 1)
    u_i <- if (have_input) matrix(input[i, ], ncol = 1) else NULL
    
    if (i == 1) {
      xp_i <- Phi %*% mu0
      if (have_input) xp_i <- xp_i + Ups %*% u_i
      Pp_i <- .symmetrize(Phi %*% Sigma0 %*% t(Phi) + Q)
    } else {
      xp_i <- Phi %*% out$xf[, , i - 1]
      if (have_input) xp_i <- xp_i + Ups %*% u_i
      Pp_i <- .symmetrize(Phi %*% out$Pf[, , i - 1] %*% t(Phi) + Q)
    }
    
    if (have_input) {
      innov_i <- y_i - B %*% xp_i - Gam %*% u_i
    } else {
      innov_i <- y_i - B %*% xp_i
    }
    
    S_i <- .symmetrize(B %*% Pp_i %*% t(B) + R)
    Sfac <- .safe_spd_solve(S_i)
    Sinv <- Sfac$Sinv
    
    K_i <- Pp_i %*% t(B) %*% Sinv
    xf_i <- xp_i + K_i %*% innov_i
    
    # joseph-stabilized covariance update
    I_p <- diag(1, pdim)
    IKH <- I_p - K_i %*% B
    Pf_i <- .symmetrize(IKH %*% Pp_i %*% t(IKH) + K_i %*% R %*% t(K_i))
    
    quad <- as.numeric(t(innov_i) %*% Sinv %*% innov_i)
    like_accum <- like_accum + Sfac$logdet + quad
    
    out$xp[, , i] <- xp_i
    out$Pp[, , i] <- Pp_i
    out$xf[, , i] <- xf_i
    out$Pf[, , i] <- Pf_i
    out$innov[, , i] <- innov_i
    out$sig[, , i] <- Sfac$S
    out$Kn[, , i] <- K_i
  }
  
  out$like <- -0.5 * like_accum
  out
}