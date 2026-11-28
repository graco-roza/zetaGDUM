st5_time <- function(x) {
  if (inherits(x,"POSIXt")) return(as.numeric(x)/86400)
  if (inherits(x,"Date")) return(as.numeric(x))
  x
}

# Community rows and environmental rows align with the records, including unavailable records.
st5_input <- function(records, incidence = NULL, abundance = NULL, env = NULL) {
  records <- as.data.frame(records)
  if (!all(c("site", "time") %in% names(records))) stop("Supply site and time columns.")
  records$site <- as.character(records$site)
  records$time <- st5_time(records$time)
  if (!is.numeric(records$time) || any(!is.finite(records$time)) ||
      anyNA(records$site) || any(records$site == "") ||
      anyDuplicated(records[c("site","time")])) stop("Site-time keys must be unique and valid.")
  if (is.null(records$available)) records$available <- TRUE
  if (anyNA(records$available) || !all(records$available %in% c(0,1)))
    stop("Availability must be binary.")
  observed <- which(records$available == 1)
  if (!length(observed)) stop("No observed communities.")
  check <- function(a, binary) {
    a <- as.matrix(a)
    if (nrow(a) != nrow(records) || !ncol(a) ||
        anyDuplicated(colnames(a)) || !(is.numeric(a) || is.logical(a)))
      stop("Invalid aligned community matrix.")
    a <- a[observed,,drop = FALSE]
    if (anyNA(a) || any(!is.finite(a) | a < 0) || (binary && any(!a %in% 0:1)))
      stop("Invalid values in an observed community.")
    a
  }
  if (is.null(incidence) && is.null(abundance)) stop("Supply incidence or abundance.")
  a <- if (!is.null(abundance)) check(abundance, FALSE)
  x <- if (is.null(incidence)) a > 0 else check(incidence, TRUE) > 0
  if (!is.null(a) && (!identical(dim(x),dim(a)) ||
      !identical(colnames(x),colnames(a)) || any(x != (a > 0))))
    stop("Incidence and abundance must identify the same positive occurrences.")
  e <- if (!is.null(env)) {
    e <- as.matrix(env)
    if (!is.numeric(e) || nrow(e) != nrow(records) || is.null(colnames(e)) ||
        anyDuplicated(colnames(e)) || any(colnames(e) == ""))
      stop("Supply named numeric environmental columns aligned to records.")
    e[observed,,drop = FALSE]
  } else matrix(numeric(), length(observed), 0L)
  sites <- sort(unique(records$site))
  times <- sort(unique(records$time))
  index <- matrix(NA_integer_, length(sites), length(times), dimnames = list(sites, times))
  index[cbind(match(records$site[observed],sites),match(records$time[observed],times))] <- seq_along(observed)
  list(records = records[observed,,drop = FALSE], x = x, a = a, env = e,
       sites = sites, times = times, index = index)
}

st5_sets <- function(n, k, budget) {
  if (n < k) return(matrix(integer(),0L,k))
  if (choose(n,k) <= budget) return(t(utils::combn(n,k)))
  if (!is.finite(budget)) stop("Combination space exceeds representable size.")
  # IID uniform unordered combinations, with replacement. Duplicates preserve the draw weights.
  t(vapply(seq_len(budget), function(i) sort(sample.int(n,k)), integer(k)))
}

st5_support <- function(p, type = c("spatial","local","rectangle","change","shared"),
                        k = 2L, m = 1L, pairs = c("consecutive","all"), lag = NULL,
                        windows = NULL, budget = 128L, seed = 15L) {
  type <- match.arg(type)
  pairs <- match.arg(pairs)
  if (any(!is.finite(c(k,m))) || any(c(k,m) < 1 | c(k,m) != floor(c(k,m))) ||
      length(budget) != 1L || is.na(budget) || budget < 1 || budget != floor(budget))
    stop("Orders and budget must be positive integers. The budget may be infinite.")
  if (!is.null(lag) && (any(!is.finite(lag) | lag <= 0) || !length(lag))) stop("Invalid lag.")
  if ((type == "spatial" && (m != 1L || k < 2L)) ||
      (type == "rectangle" && (m < 2L || k < 2L))) stop("Invalid static response orders.")
  if (type == "local") {k <- 1L; m <- 1L}
  old_seed <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit(if (is.null(old_seed)) rm(".Random.seed", envir = .GlobalEnv)
          else assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
  set.seed(seed)
  matched <- function(w) which(rowSums(is.na(p$index[,w,drop = FALSE])) == 0)
  ts <- seq_along(p$times)
  selected <- list()
  make <- function(a,b = integer(), focal = NULL) {
    region <- matched(c(a,b))
    keep <- if (is.null(focal)) region else intersect(region,focal)
    if (length(keep) < k) return(NULL)
    list(a = a, b = b, region = region, focal = keep, sets = matrix(keep[st5_sets(length(keep),k,budget)],ncol = k),
         eligible = choose(length(keep),k))
  }
  validate_window <- function(w) {
    q <- match(st5_time(w),p$times)
    if (length(q) != m || anyNA(q) || is.unsorted(q, strictly = TRUE))
      stop("Windows must use ordered, distinct recorded times and the requested order.")
    q
  }
  consecutive <- function(i,m) {
    if (length(i) < m) return(list())
    lapply(seq_len(length(i)-m+1L),function(j) i[j+seq_len(m)-1L])
  }
  time_pairs <- function(i) {
    if (length(i) < 2L) return(matrix(integer(),0,2))
    z <- if (pairs == "all" || !is.null(lag)) t(combn(i,2)) else cbind(head(i,-1),tail(i,-1))
    if (!is.null(lag)) {
      elapsed <- p$times[z[,2]] - p$times[z[,1]]
      z <- z[vapply(elapsed,function(d) any(abs(d-lag) <= 1e-8 * pmax(1,lag)),logical(1)),,drop=FALSE]
    }
    z
  }
  if (type == "local") {
    keys <- lapply(seq_along(p$sites), function(i) {
      z <- time_pairs(which(!is.na(p$index[i,])))
      data.frame(a=z[,1], b=z[,2],site=rep(i,nrow(z)))
    })
    keys <- do.call(rbind,keys)
    selected <- lapply(split(keys,interaction(keys$a,keys$b,drop=TRUE)), function(z)
      make(z$a[1],z$b[1],z$site))
  } else if (type %in% c("spatial","rectangle")) {
    ws <- if (!is.null(windows)) lapply(windows,validate_window) else consecutive(ts,m)
    selected <- lapply(ws,make)
  } else if (!is.null(windows)) {
    selected <- lapply(windows,function(w) {
      a <- validate_window(w[[1]]); b <- validate_window(w[[2]])
      if (max(a) >= min(b)) stop("Paired windows must be non-overlapping and ordered.")
      make(a,b)
    })
  } else if (m == 1L) {
    z <- time_pairs(ts)
    selected <- lapply(seq_len(nrow(z)), function(j) make(z[j,1],z[j,2]))
  } else {
    w <- consecutive(ts,m)
    selected <- unlist(lapply(seq_along(w),function(i) {
      ii <- if (pairs == "all" || !is.null(lag)) which(seq_along(w) >= i+m) else intersect(i+m,seq_along(w))
      if (!is.null(lag)) ii <- ii[vapply(ii,function(j)
        any(abs(mean(p$times[w[[j]]] - p$times[w[[i]]])-lag) <= 1e-8*pmax(1,lag)),logical(1))]
      lapply(ii,function(j) make(w[[i]],w[[j]]))
    }),recursive=FALSE)
  }
  structure(Filter(Negate(is.null),selected), type=type, k=k,m=m)
}

st5_observation_ids <- function(p,sites,times) {
  if (!length(times)) return(matrix(integer(),nrow(sites),0))
  do.call(cbind,lapply(times,function(t) matrix(p$index[,t][sites],nrow(sites),ncol(sites))))
}

st5_eligibility <- function(p, support) {
  n <- vapply(support,function(g) nrow(g$sets),integer(1))
  eligible <- vapply(support,`[[`,numeric(1),"eligible")
  ids <- unique(unlist(lapply(support,function(g) st5_observation_ids(p,g$sets,c(g$a,g$b)))))
  eligible_ids <- unique(unlist(lapply(support,function(g) {
    sites <- g$focal
    p$index[sites,c(g$a,g$b),drop=FALSE]
  })))
  data.frame(type=attr(support,"type"), k=attr(support,"k"),m=attr(support,"m"),
             periods=length(n), combination_rows=sum(n), eligible_rows=sum(eligible),
             sampled_fraction=if (sum(eligible)) sum(n)/sum(eligible) else NA_real_,
             observation_fraction=length(ids)/nrow(p$x),
             eligible_observation_fraction=length(eligible_ids)/nrow(p$x))
}
