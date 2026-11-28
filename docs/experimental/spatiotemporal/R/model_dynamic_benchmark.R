# Fixed neighbourhoods contain named site IDs and remain fixed across orders.
# For evaluation at new sites their adjacent environmental data may be known,
# without using their community records.
st2_neighbourhood_design <- function(design, env, neighbourhood) {
  if (is.null(neighbourhood)) return(design)
  if (design$k != 1L || design$windows$temporal_order != 1L || !is.list(neighbourhood) ||
      !all(design$sites %in% names(neighbourhood))) stop("Specify fixed named neighbourhoods at order (1,1).")
  dn <- dimnames(env)
  vars <- names(design$basis$center)
  if (!all(design$sites %in% dn[[1L]]) || !identical(dn[[3L]], vars)) {
    stop("Name all neighbourhood environmental site and predictor axes.")
  }
  # Use time labels to align env even when the panel uses a temporal split.
  if (!is.null(dn[[2L]]) && all(as.character(design$times) %in% dn[[2L]])) {
    env <- env[, as.character(design$times), , drop = FALSE]
  }
  if (dim(env)[2L] != length(design$times)) stop("Align neighbourhood environment to the design times.")
  by_site <- lapply(design$sites, function(s) {
    neighbours <- neighbourhood[[s]]
    if (!length(neighbours) || s %in% neighbours || anyDuplicated(neighbours) ||
        !all(neighbours %in% dn[[1L]])) stop("Every neighbourhood needs distinct available neighbours.")
    vapply(seq_along(vars), function(j) {
      gaps <- sweep(matrix(env[neighbours, , j], length(neighbours)), 2L, env[s, , j], "-")
      colMeans(abs(gaps)) / design$basis$scale[j]
    }, numeric(length(design$times)))
  })
  w <- design$rows$window[design$row_id]
  heterogeneity <- function(time) t(vapply(seq_along(w), function(i) {
    by_site[[design$site[i]]][time[w[i], 1L], ]
  }, numeric(length(vars))))
  a <- matrix(heterogeneity(design$windows$from), ncol = length(vars))
  b <- matrix(heterogeneity(design$windows$to), ncol = length(vars))
  design$X <- cbind(design$X, setNames(as.data.frame(a), paste0("neighbourhood.", vars)),
                    setNames(as.data.frame(b - a), paste0("neighbourhood_change.", vars)))
  design$X <- as.matrix(design$X)
  design
}

.st2_frame <- function(X, selected, species) {
  idx <- which(selected, arr.ind = TRUE)
  frame <- as.data.frame(X[idx[, 1L], , drop = FALSE])
  names(frame) <- paste0("x", seq_len(ncol(frame)))
  frame$species <- factor(species[idx[, 2L]], levels = species)
  list(data = frame, index = idx)
}

.st2_formula <- function(frame, hierarchical) {
  terms <- setdiff(names(frame), c("y", "species"))
  if (hierarchical) terms <- c(terms, "(1 | species)")
  stats::as.formula(paste("y ~", if (length(terms)) paste(terms, collapse = " + ") else "1"))
}

.st2_check_rank <- function(frame) {
  X <- cbind(intercept = 1, as.matrix(frame[, !names(frame) %in% c("y", "species"), drop = FALSE]))
  qr(X)$rank == ncol(X)
}

st2_fit_dynamic <- function(x, design, env = NULL, neighbourhood = NULL,
                            interaction = FALSE, nAGQ = 0L) {
  if (design$k != 1L || design$windows$temporal_order != 1L) stop("M2 uses original site and survey states.")
  if (!requireNamespace("lme4", quietly = TRUE)) stop("Install lme4 to fit M2.")
  design <- st2_neighbourhood_design(design, env, neighbourhood)
  states <- st2_states(x, design)
  species <- colnames(states$z)
  terms <- .st2_terms(design, interaction)
  X <- design$X[, terms, drop = FALSE]
  fits <- lapply(c("gain", "loss"), function(event) {
    index <- if (event == "gain") !states$z else states$z
    y <- if (event == "gain") states$y[index] else !states$y[index]
    frame <- .st2_frame(X, index, species)$data
    frame$y <- y
    if (!length(y)) return(list(status = "no_opportunities"))
    if (all(y == y[1L])) return(list(status = "boundary", rate = as.numeric(y[1L])))
    if (!.st2_check_rank(frame)) return(list(status = "aliased"))
    model <- suppressMessages(lme4::glmer(.st2_formula(frame, TRUE), data = frame,
                                          family = stats::binomial(), nAGQ = nAGQ))
    list(status = if (length(model@optinfo$conv$lme4$messages)) "convergence_message" else "estimated",
         model = model, singular = lme4::isSingular(model))
  })
  names(fits) <- c("gain", "loss")
  list(fit = fits, design = design, terms = terms, species = species, neighbourhood = neighbourhood)
}

st2_predict_dynamic <- function(object, design, env = NULL) {
  if (design$k != 1L || design$windows$temporal_order != 1L ||
      !identical(design$basis, object$design$basis)) stop("M2 predictions require the fitted basis at order (1,1).")
  design <- st2_neighbourhood_design(design, env, object$neighbourhood)
  X <- design$X[, object$terms, drop = FALSE]
  use <- matrix(TRUE, nrow(X), length(object$species))
  frame <- .st2_frame(X, use, object$species)$data
  lapply(object$fit, function(f) {
    if (f$status == "boundary") return(matrix(f$rate, nrow(X), ncol(use)))
    if (is.null(f$model)) return(matrix(NA_real_, nrow(X), ncol(use)))
    matrix(stats::predict(f$model, frame, type = "response", allow.new.levels = TRUE), nrow(X))
  })
}

.st2_positive_features <- function(design, states, terms, persistent, chosen) {
  X <- design$X[, terms, drop = FALSE]
  rows <- .st2_frame(X, chosen, colnames(states$z))
  if (design$k > 1L) rows$data$other_occupancy <- states$other[rows$index] / (design$k - 1L)
  if (persistent) rows$data$log_initial <- log(states$before[rows$index])
  rows
}

st2_fit_positive <- function(x, design, interaction = FALSE, hierarchical = FALSE,
                             family = c("lognormal", "truncated_nbinom")) {
  family <- match.arg(family)
  if (design$windows$temporal_order != 1L) stop("Positive abundance uses original survey observations.")
  if (hierarchical && !requireNamespace("lme4", quietly = TRUE)) stop("Install lme4 for hierarchical log abundance.")
  if (family == "truncated_nbinom" && !requireNamespace("glmmTMB", quietly = TRUE)) stop("Install glmmTMB for count sensitivity.")
  states <- st2_states(x, design)
  if (family == "truncated_nbinom" && any(states$after != floor(states$after))) stop("Count fits require integer abundances.")
  terms <- .st2_terms(design, interaction)
  fits <- lapply(c(FALSE, TRUE), function(persistent) {
    chosen <- states$y & states$z == persistent
    if (!any(chosen)) return(list(status = "no_opportunities"))
    records <- .st2_positive_features(design, states, terms, persistent, chosen)
    frame <- records$data
    frame$y <- states$after[records$index]
    if (!.st2_check_rank(frame)) return(list(status = "aliased"))
    formula <- .st2_formula(frame, hierarchical)
    frame$obs_weight <- design$weight[records$index[, 1L]]
    if (family == "truncated_nbinom") {
      model <- glmmTMB::glmmTMB(formula, frame, weights = obs_weight, family = glmmTMB::truncated_nbinom2())
      status <- if (isTRUE(model$sdr$pdHess) && model$fit$convergence == 0L) "estimated" else "nonconverged"
    } else {
      frame$y <- log(frame$y)
      if (hierarchical) {
        model <- suppressMessages(lme4::lmer(formula, frame, weights = obs_weight, REML = FALSE))
        status <- if (length(model@optinfo$conv$lme4$messages)) "convergence_message" else "estimated"
      } else {
        model <- stats::lm(formula, frame, weights = obs_weight)
        status <- "estimated"
      }
    }
    list(status = status, model = model, scale = stats::sigma(model))
  })
  names(fits) <- c("new", "persistent")
  list(fit = fits, terms = terms, design = design, species = colnames(states$z), family = family)
}

st2_predict_positive_parameters <- function(object, design, states) {
  if (design$k != object$design$k || design$windows$temporal_order != 1L ||
      !identical(design$basis, object$design$basis) ||
      !identical(colnames(states$z), object$species)) stop("Positive abundance design differs from training.")
  location <- scale <- matrix(NA_real_, nrow(states$z), ncol(states$z))
  for (persistent in c(FALSE, TRUE)) {
    f <- object$fit[[if (persistent) "persistent" else "new"]]
    chosen <- states$z == persistent
    if (is.null(f$model) || !any(chosen)) next
    records <- .st2_positive_features(design, states, object$terms, persistent, chosen)
    if (object$family == "truncated_nbinom") {
      records$data$obs_weight <- 1
      location[records$index] <- stats::predict(f$model, records$data, type = "link", allow.new.levels = TRUE)
    } else {
      location[records$index] <- if (inherits(f$model, "lm")) stats::predict(f$model, records$data) else
        stats::predict(f$model, records$data, allow.new.levels = TRUE)
    }
    scale[records$index] <- f$scale
  }
  list(location = location, scale = scale, family = object$family)
}
