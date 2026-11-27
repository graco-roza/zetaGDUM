# Run this script from the repository root with the development package installed.
source("docs/final_gdum/design.R")
a <- generate(168L, "spatial", 71104L)
coordinates <- cbind(km = seq_len(nrow(a$e)), transect = 0)
block <- cut(coordinates[,1L], 4L, labels = FALSE)
orders <- 2:8
profile <- zeta_order_information(a$x,a$e,orders=orders,response="all",seed=18L)

# Declare the spatial partition before any model comparison.
k <- 4L
cv <- ugdm_validate_components(a$x,a$e,coordinates,order=k,fold=block,
  combinations=1200L,seed=971L)
print(aggregate(D2 ~ architecture + model + component,cv,median,na.rm=TRUE))

# Refit an illustrative model on all sites. Compare prediction surfaces only
# in bins supported by actual site configurations from both training pools.
idx <- sample_sets(nrow(a$x),k,1800L,seed=98L)
eval <- sample_sets(nrow(a$x),k,4000L,seed=242L)
y <- ugdm_component_response(a$x,idx)
v <- ugdm_component_response(a$x,eval)
d <- ugdm_component_design(a$e,idx,"interaction")
z <- ugdm_component_design(a$e,eval,reference=d)
fit <- ugdm_fit_components(d,y,architecture="coherent")
importance <- ugdm_component_importance(fit,d,y,evaluation=z,evaluation_response=v)
print(importance$hierarchy)
print(ugdm_component_diagnostics(fit,d,y)$stages)
surface <- ugdm_component_surface(fit,z,"phosphorus",bins=8,min_sets=8,min_sites=12)
print(table(surface$supported))

# Repeat numerical integration with the same sites and evaluation set.
pools <- lapply(c(47,58,61,78,93),function(seed) {
  j <- sample_sets(nrow(a$x),k,1800L,seed=seed)
  g <- ugdm_component_design(a$e,j,reference=d)
  r <- ugdm_component_response(a$x,j)
  f <- ugdm_fit_components(g,r)
  ugdm_component_importance(f,g,r,evaluation=z,evaluation_response=v,
    resolution="predictor")$predictors
})
confirmation <- lapply(c("total","occupancy","heterogeneity"),function(component) {
  integration <- t(vapply(pools,function(x) {
    a <- x[x$response==component,]
    setNames(a$absolute,a$predictor)
  },setNames(numeric(ncol(a$e)),names(a$e))))
  spatial <- subset(cv,architecture=="coherent" & model=="interaction")
  spatial <- spatial[spatial$component==component,]
  independent_D2 <- if(all(spatial$converged) && all(is.finite(spatial$D2))) median(spatial$D2) else NA_real_
  ugdm_within_order_diagnostic(fit,d,y[,component],integration,draws=NULL,
    heldout_D2=independent_D2,sharing=y[,"total"]<1,component=component,
    min_response_sd=0,min_heldout_D2=.025,max_mc_se=.02,max_width=.15,
    max_condition=1e8,calibrated=FALSE)
})
print(setNames(vapply(confirmation,"[[",character(1),"status"),
               c("total","occupancy","heterogeneity")))

# Show conditional predictions at low and high represented positions.
# Only spread bins supported at both positions are used for comparison.
positions <- sort(unique(surface$position[surface$supported]))
selected <- positions[unique(round(quantile(seq_along(positions),c(.25,.75))))]
a1 <- subset(surface,position==selected[1L] & supported)
a2 <- subset(surface,position==selected[2L] & supported)
common <- sort(intersect(a1$spread,a2$spread))
folder <- Sys.getenv("MS_GDUM_WORKFLOW_OUT","ms_gdum_workflow_output")
dir.create(folder,recursive=TRUE,showWarnings=FALSE)
write.csv(surface,file.path(folder,"phosphorus_surface.csv"),row.names=FALSE)
if (length(common) >= 2L) {
  pdf(file.path(folder,"phosphorus_spread_sections.pdf"),width=9,height=3.6,useDingbats=FALSE)
  par(mfrow=c(1,3),mar=c(4,4,2,1))
  for (component in c("total","occupancy","heterogeneity")) {
    values <- cbind(a1[[component]][match(common,a1$spread)],
                    a2[[component]][match(common,a2$spread)])
    matplot(common,values,type="b",pch=c(16,1),lty=1:2,col=c("#216A8A","#A64B30"),
      main=component,xlab="Mean pairwise phosphorus difference",
      ylab="Predicted dissimilarity")
    legend("bottomright",sprintf("Position %.1f",selected),lty=1:2,
      col=c("#216A8A","#A64B30"),bty="n",cex=.8)
  }
  invisible(dev.off())
}
message("Workflow outputs are in ", folder)
