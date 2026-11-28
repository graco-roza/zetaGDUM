source("docs/experimental/spatiotemporal/phase4_swedish_lakes/scripts/common.R")
a <- archive()
sets <- list(A_original=c(1992:1996,2018:2022),A=c(1995:1999,2018:2022),
             B_original=1992:2022,B=1995:2022)
membership <- design <- coverage <- list()
for (group in assemblages) {
  z <- a[[group]]$samples
  ids <- sort(unique(z$id))
  wide <- matrix(FALSE,length(ids),length(years),dimnames=list(ids,years))
  available <- z$status=="observed_nonempty"
  wide[cbind(match(z$id[available],ids),match(z$year[available],years))] <- TRUE
  rownames(wide) <- ids
  write.csv(data.frame(id=ids,wide,check.names=FALSE),processed(paste0(group,"_coverage.csv")),row.names=FALSE)
  coverage[[group]] <- wide
  for (s in names(sets)) {
    chosen <- sets[[s]]
    selected <- wide[,as.character(chosen),drop=FALSE]
    include <- rowSums(selected)==length(chosen)
    reasons <- apply(!selected,1,function(x) paste(chosen[x],collapse=" "))
    membership[[length(membership)+1L]] <- data.frame(assemblage=group,design=s,id=ids,included=include,
                                        missing_or_undefined_years=reasons)
    panel <- if (any(include)) st4_panel(a[[group]],ids[include],chosen) else NULL
    design[[length(design)+1L]] <- data.frame(assemblage=group,design=s,lakes=sum(include),
          taxa=if(is.null(panel))0L else sum(colSums(panel$x>0,dims=2)>0),
          start=min(chosen),end=max(chosen),surveys=length(chosen))
  }
}
membership <- do.call(rbind,membership)
save_csv(membership,"panel_membership")
save_csv(do.call(rbind,design),"panel_summary")
saveRDS(list(membership=membership,years=sets,coverage=coverage),processed("panels.rds"))
