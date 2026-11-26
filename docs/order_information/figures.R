source("docs/order_information/design.R")
prof <- do.call(rbind,lapply(list.files(output,"^profile_.*csv.gz$",full.names=TRUE),
                            function(f) read.csv(gzfile(f))))
keep <- c("narrow","broad","dominant")
cols <- c(narrow="#1B6586", broad="#A44A23", dominant="#4F794B")
labels <- c(narrow="Narrow niches", broad="Broad niches", dominant="Ubiquitous dominants")
p <- subset(prof, name %in% keep & split=="test" & n==150 & response=="relative")
p <- aggregate(cbind(zeta_remaining,retention_rate,relative_information,position_variance,
    nonlinear_position_variance,spread_variance,spline_spread_variance)~name+order,p,median)
m <- read.csv(file.path(output,"reliability.csv"))
r <- subset(m,name %in% keep & split=="test" & n==150 & response=="relative")
s <- read.csv(file.path(output,"selection_predictions.csv"))
padded_range <- function(z) {
  z <- z[is.finite(z)]
  v <- range(z,0)
  if(diff(v)==0) v <- v+c(0,1)
  v
}
axis_panel <- function(v,ylab,tag,title,xlim=c(2,40),ylim=NULL) {
  if(is.null(ylim)) ylim <- padded_range(v)
  plot(NA_real_,xlim=xlim,ylim=ylim, xlab="Zeta order",ylab=ylab, las=1,
       bty="l",xaxs="i",yaxs="i",ann=TRUE)
  graphics::title(main=paste(tag,title),adj=0,line=.75,cex.main=1)
}
curves <- function(d,vars,which_type) {
  for (name in keep) {
    a <- d[d$name==name,]
    a <- a[order(a$order),]
    for (j in seq_along(vars)) lines(a$order,which_type(a[[vars[j]]]),col=cols[name],lty=j,lwd=1.65)
  }
}
figure <- function() {
  par(mfrow=c(2,2),mar=c(3.4,4.25,2.55,.7),oma=c(0,0,3,0),mgp=c(2.5,.7,0),
      tcl=-.27,cex.axis=.78,cex.lab=.86)
  vars <- c("zeta_remaining","retention_rate","relative_information")
  axis_panel(unlist(p[vars]),"Fraction or ratio to order 2","A","Ecological information")
  curves(p,vars,identity)
  legend(x=24,y=.72,c("Zeta / pairwise zeta","Retention","Response variance / pairwise"),
    lty=1:3,lwd=1.7,col="#454A50",cex=.71,bty="n")
  vars <- c("spread_variance","nonlinear_position_variance")
  axis_panel(NULL,"Log10 design variance","B","Environmental information",
    ylim=range(log10(unlist(p[c(vars,"position_variance")])))*c(1.02,.92))
  curves(p,vars,log10)
  a <- subset(p,name=="dominant")
  a <- a[order(a$order),]
  lines(a$order,log10(a$position_variance),lwd=1.4,lty=3,col="#323437")
  legend("topright",c("Pairwise spread","Nonlinear position","Linear position"),
    lty=1:3,lwd=1.7,cex=.71,bty="n")
  vars <- c("heldout_median","difference_width95","leading")
  axis_panel(unlist(r[vars]),"Explained deviance or proportion","C","Modelling reliability",
    ylim=c(min(0, r$heldout_median, na.rm=TRUE),1.05))
  abline(h=0,col="#A8ADB2")
  curves(r,vars,identity)
  legend(x=24,y=.72,c("Held-out explained deviance","95% contrast interval width","Correct leader"),
         lty=1:3,lwd=1.7,cex=.70,bty="n")

  a <- subset(s,target=="inference" & rule=="combined" & response=="relative" &
    evaluation %in% c("new_communities","unseen_scenarios"))
  a$group <- match(a$name,c(keep,"correlated","spatial","null"))
  regime_cols <- c(inferentially_supported="#257954",prediction_only="#DBA140",unsupported="#BAC2CA")
  par(mar=c(3.4,6.2,2.55,.7))
  plot(NA_real_,xlim=c(2,40),ylim=c(.65,7.8),yaxt="n",xlab="Evaluated zeta order",ylab="",
    xaxs="i",las=1,bty="l")
  title(main="D  Order-selection accuracy",adj=0,line=.75,cex.main=1)
  axis(2,at=6:1,labels=c("Narrow, 150","Broad, 70","Dominant, 250","Correlated, 150",
                           "Spatial, 250","Null, 70"),las=1,tick=FALSE,cex.axis=.74)
  abline(h=1:6,col="#EAEEF0")
  points(a$order,7-a$group,pch=15,col=regime_cols[a$regime],cex=.95)
  points(a$order[a$predicted],7-a$group[a$predicted],pch=0,cex=1.24)
  legend(x=2,y=7.8,c("Inference","Prediction only","Unsupported","Rule recommends"),
         pch=c(15,15,15,0),col=c(unname(regime_cols),"black"),ncol=2,cex=.68,bty="n")
  for (j in seq_along(keep)) {
    mtext(labels[keep[j]], side=3, line=1.1, at=c(.22,.50,.78)[j],
          outer=TRUE, col=cols[keep[j]], font=2, cex=.87)
  }
}
supplement <- function() {
  par(mfrow=c(2,2),mar=c(3.4,4.25,2.55,.7),mgp=c(2.5,.7,0),las=1,tcl=-.27,cex.axis=.78,cex.lab=.86)
  a <- subset(m,split=="principal" & response=="relative")
  vars <- c("gradient_coverage95","secondary_coverage95","difference_coverage95")
  axis_panel(unlist(a[vars]),"Observed coverage","A","Coverage by estimand",ylim=c(0,1.05))
  abline(h=.95,col="#7E848B",lty=4)
  curves(a,vars,identity)
  legend("bottomright",c("Gradient","Secondary","Difference"),lty=1:3,lwd=1.7,cex=.72,bty="n")
  a <- subset(m,response=="relative" & split %in% c("principal","budget") & order<=12)
  axis_panel(a$difference_mc_rmse,"Numerical RMSE","B","Integration budget",xlim=c(2,12))
  for (name in keep) for (split in c("principal","budget")) {
    d <- a[a$name==name & a$split==split,]
    d <- d[order(d$order),]
    lines(d$order,d$difference_mc_rmse,col=cols[name],lty=if(split=="principal") 1 else 2,lwd=1.65)
  }
  legend("topright",c("900 sets","3,600 sets"),lty=1:2,lwd=1.7,bty="n",cex=.72)
  a <- subset(prof,name %in% keep & split=="test" & n==150 & response=="relative")
  a <- aggregate(cbind(response_variance,variance_se)~name+order,a,median)
  axis_panel(NULL,"Log10 dissimilarity variance","C","Absolute compositional variation",
             ylim=range(log10(a$response_variance)))
  curves(a,"response_variance",log10)
  axis_panel(NULL,"Log10 variance estimate SE","D","Diagnostic numerical precision",
             ylim=range(log10(a$variance_se[a$variance_se>0])))
  curves(a,"variance_se",log10)
}
make <- function(filename,drawing,width=9.5,height=7.1) {
  grDevices::pdf(file.path(output,paste0(filename,".pdf")),width,height,useDingbats=FALSE)
  drawing()
  invisible(dev.off())
  grDevices::svg(file.path(output,paste0(filename,".svg")),width,height)
  drawing()
  invisible(dev.off())
}
make("figure3_order_selection",figure)
make("figureS1_order_sensitivity",supplement)
