st3_processes <- c("expansion", "contraction", "redistribution", "numerical")

st3_plan <- function(x, before, after) {
  p <- .st_panel(x)
  n <- dim(p$x)[1]
  S <- dim(p$x)[3]
  before <- .st_orders(before, length(p$times), "before")
  after <- .st_orders(after, length(p$times), "after")
  if (length(before) != length(after) || any(diff(before) <= 0) || any(diff(after) <= 0) ||
      max(before) >= min(after)) stop("Use equally sized, ordered, non-overlapping survey windows.")
  m <- length(before)
  species <- dimnames(p$x)[[3L]]
  a <- matrix(p$x[, before, , drop = FALSE], n * m, S, dimnames = list(NULL, species))
  b <- matrix(p$x[, after, , drop = FALSE], n * m, S, dimnames = list(NULL, species))
  mask <- matrix(0L, n * m, S)
  mask[a > 0 & b > 0 & a != b] <- 8L
  ambiguous <- list()
  rows <- vector("list", m)
  for (j in seq_len(m)) {
    cell <- seq_len(n) + n * (j - 1L)
    x0 <- a[cell, , drop = FALSE] > 0
    x1 <- b[cell, , drop = FALSE] > 0
    gain <- !x0 & x1
    loss <- x0 & !x1
    g <- colSums(gain)
    l <- colSums(loss)
    if (any(g + l > 0)) for (s in which(g + l > 0)) {
      G <- cell[gain[, s]] + (s - 1L) * nrow(a)
      L <- cell[loss[, s]] + (s - 1L) * nrow(a)
      if (!length(L)) {
        mask[G] <- 1L
      } else if (!length(G)) {
        mask[L] <- 2L
      } else {
        mask[c(G,L)] <- 4L
        if (g[s] != l[s]) {
          expansion <- g[s] > l[s]
          ambiguous[[length(ambiguous) + 1L]] <- list(
            cells = if (expansion) G else L,
            count = abs(g[s] - l[s]), code = if (expansion) 1L else 2L)
        }
      }
    }
    rows[[j]] <- data.frame(position = j, species, gains = g, losses = l,
                           expansion = pmax(g-l, 0), contraction = pmax(l-g, 0),
                           redistribution = pmin(g,l),
                           regional_appearance = colSums(x0) == 0 & colSums(x1) > 0,
                           regional_disappearance = colSums(x0) > 0 & colSums(x1) == 0,
                           row.names = NULL)
  }
  log_count <- sum(vapply(ambiguous, function(z) lchoose(length(z$cells), z$count), numeric(1L)))
  list(before = a, after = b, mask = mask, ambiguous = ambiguous, log_allocation_count = log_count,
       ledger = do.call(rbind, rows), sites = dimnames(p$x)[[1]],
       alignment = data.frame(initial_index = before, final_index = after,
                              initial = p$times[before], final = p$times[after],
                              lag = p$times[after] - p$times[before]))
}

# Uniform over majority-site subsets. No pairing of individual gains to losses.
st3_allocations <- function(plan, exact_limit = 32L, draws = 32L, seed = 1L) {
  if (length(exact_limit) != 1L || length(draws) != 1L ||
      any(!is.finite(c(exact_limit, draws)) | c(exact_limit, draws) < c(1, 2)) ||
      any(c(exact_limit, draws) != floor(c(exact_limit, draws)))) {
    stop("Use a positive integer exact limit and at least two draws.")
  }
  if (!length(plan$ambiguous)) return(list(masks = list(plan$mask), exact = TRUE))
  choices <- plan$ambiguous
  exact <- plan$log_allocation_count <= log(exact_limit) + 1e-12
  assign_blocks <- function(selected) {
    mask <- plan$mask
    for (i in seq_along(choices)) mask[selected[[i]]] <- choices[[i]]$code
    mask
  }
  if (exact) {
    selections <- lapply(choices, function(z) utils::combn(z$cells, z$count, simplify = FALSE))
    grid <- expand.grid(lapply(selections, seq_along), KEEP.OUT.ATTRS = FALSE)
    masks <- lapply(seq_len(nrow(grid)), function(i) {
      assign_blocks(lapply(seq_along(selections), function(j) selections[[j]][[grid[i,j]]]))
    })
  } else {
    old <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(if (is.null(old)) rm(list = ".Random.seed", envir = .GlobalEnv)
            else assign(".Random.seed", old, envir = .GlobalEnv), add = TRUE)
    set.seed(seed)
    masks <- lapply(seq_len(draws), function(i) {
      assign_blocks(lapply(choices, function(z) z$cells[sample.int(length(z$cells), z$count)]))
    })
  }
  list(masks = masks, exact = exact)
}
