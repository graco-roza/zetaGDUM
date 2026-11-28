st5_basis <- function(p, z, predictors = colnames(p$env)) {
  ids <- sort(unique(unlist(c(z$rows$a,z$rows$b))))
  E <- p$env[,predictors,drop=FALSE]
  supported <- vapply(seq_along(predictors),function(j)
    length(unique(E[ids,j][is.finite(E[ids,j])])) >= 4L,logical(1))
  vars <- predictors[supported]
  # A shared basis sample uses training observation keys with all selected predictors.
  repeat {
    rows <- ids[rowSums(!is.finite(E[ids,vars,drop=FALSE]))==0]
    enough <- vapply(vars,function(v) length(unique(E[rows,v]))>=4L,logical(1))
    if(all(enough)) break
    worst <- which.max(colSums(!is.finite(E[ids,vars,drop=FALSE])))
    vars <- vars[-worst]
  }
  basis <- if (length(vars)) zetaGDUM::ugdm_design(E[rows,vars,drop=FALSE],matrix(1:2,1))$environment_basis
  list(native=basis,vars=vars,unsupported=setdiff(predictors,vars),training_ids=rows)
}

st5_design <- function(p,z,basis,lag_ref=1) {
  if (!is.finite(lag_ref) || lag_ref<=0) stop("Use a positive lag reference.")
  r <- z$rows
  n <- nrow(r)
  paired <- z$type %in% c("change","local","shared")
  k <- unique(r$k)
  if (length(k)!=1L || length(unique(r$m))!=1L) stop("Fit at a single supported order.")
  A <- do.call(rbind,r$a)
  B <- if(paired) do.call(rbind,r$b) else A
  values <- columns <- list()
  extrapolation <- missing <- rep(FALSE,n)
  add <- function(x,name,predictor="",group="nuisance") {
    x <- as.matrix(x)
    labels <- if(ncol(x)>1L) paste0(name,seq_len(ncol(x))) else name
    colnames(x) <- labels
    values[[length(values)+1L]] <<- x
    columns[[length(columns)+1L]] <<- data.frame(name=labels,predictor,group)
  }
  if (z$type=="change") add(z$baseline,"baseline")
  if (paired) add(log(r$lag/lag_ref),"log_lag")
  if (ncol(A)>k) {
    add(log1p(r$span/lag_ref),"initial_duration")
    if (paired) add(log1p(r$final_span/lag_ref)-log1p(r$span/lag_ref),"change_duration")
  }
  for (j in seq_along(basis$vars)) {
    v <- basis$vars[j]
    ref <- basis$native
    env <- (p$env[,v]-ref$center[v])/ref$scale[v]
    knot <- attr(ref$splines[[j]],"Boundary.knots")
    shp <- function(x) {
      ok <- is.finite(x)
      out <- matrix(NA_real_,length(x),3)
      if(any(ok)) out[ok,] <- stats::predict(ref$splines[[j]],pmin(pmax(x[ok],knot[1]),knot[2]))
      out
    }
    summarize <- function(ids,temporal=FALSE) {
      raw <- matrix(env[ids],n,ncol(ids))
      means <- Reduce(`+`,lapply(seq_len(ncol(ids)/k),function(i)
        raw[,seq_len(k)+(i-1)*k,drop=FALSE])) / (ncol(ids)/k)
      spline <- shp(as.vector(means))
      sp <- if(k>1L) {
        pairs <- combn(k,2)
        Reduce(`+`,lapply(seq_len(ncol(pairs)),function(i)
          abs(spline[seq_len(n)+(pairs[1,i]-1)*n,,drop=FALSE] -
              spline[seq_len(n)+(pairs[2,i]-1)*n,,drop=FALSE]))) / ncol(pairs)
      }
      tm <- if(temporal) {
        at <- shp(as.vector(raw))
        tp <- combn(ncol(ids)/k,2)
        Reduce(`+`,lapply(seq_len(ncol(tp)),function(i)
          Reduce(`+`,lapply(seq_len(k),function(s)
            abs(at[seq_len(n)+((tp[1,i]-1)*k+s-1)*n,,drop=FALSE] -
                at[seq_len(n)+((tp[2,i]-1)*k+s-1)*n,,drop=FALSE]))) / k)) / ncol(tp)
      }
      list(position=rowSums(means),spread=sp,temporal=tm,raw=raw)
    }
    a <- summarize(A,z$type=="rectangle")
    b <- if(paired) summarize(B) else a
    raw <- cbind(a$raw,b$raw)
    extrapolation <- extrapolation | rowSums(raw<knot[1] | raw>knot[2],na.rm=TRUE)>0
    missing <- missing | rowSums(!is.finite(raw))>0
    add(a$position,paste0(v,":P"),v,"position")
    if(k>1L) add(a$spread,paste0(v,":S"),v,"spatial_spread")
    if(paired) {
      add(b$position-a$position,paste0(v,":dP"),v,"directional_change")
      if(k>1L) add(b$spread-a$spread,paste0(v,":dS"),v,"change_spatial_spread")
    }
    if(z$type=="local") add(abs(shp(b$raw)-shp(a$raw)),paste0(v,":T"),v,"temporal_spread")
    if(z$type=="rectangle") add(a$temporal,paste0(v,":T"),v,"temporal_spread")
  }
  list(x=if(length(values)) do.call(cbind,values) else matrix(numeric(),n,0),
       columns=if(length(columns)) do.call(rbind,columns) else data.frame(name=character(),predictor=character(),group=character()),
       extrapolation=extrapolation,missing=missing,unsupported=basis$unsupported)
}
