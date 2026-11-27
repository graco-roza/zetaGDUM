source("docs/final_gdum/design.R")
r <- read.csv(file.path(outdir,"model_summary.csv"))
u <- read.csv(file.path(outdir,"interval_summary.csv"))
emp <- read.csv(file.path(outdir,"empirical_summary.csv"))
surf <- read.csv(gzfile(file.path(outdir,"reference_surfaces.csv.gz")))
colors <- c(total="#244F69",occupancy="#2A765E",heterogeneity="#B15939")
architecture_cols <- c(coherent="#176B82",independent="#A35138")
label <- c(total="Total",occupancy="Occupancy",heterogeneity="Heterogeneity")
tag <- function(letter,title) title(main=paste(letter,title),adj=0,line=.75,cex.main=1)
figure <- function() {
 par(mfrow=c(2,2),mar=c(4.1,4.2,2.1,.45),mgp=c(2.55,.65,0),tcl=-.27,
  cex.axis=.79,cex.lab=.89)
 a <- subset(r,phase=="principal" & predictor=="phosphorus" & model=="additive")
 p <- reshape(a[c("scenario","order","response","architecture","oracle_RMSE")],
   v.names="oracle_RMSE",idvar=c("scenario","order","response"),timevar="architecture",direction="wide")
 p$error <- p$oracle_RMSE.coherent-p$oracle_RMSE.independent
 y <- range(p$error,0)
 plot(NA,xlim=c(.5,3.5),ylim=y+c(-1,1)*.08*diff(y),xaxt="n",xlab="Zeta order",
      ylab="Coherent minus independent RMSE",bty="l",las=1)
 tag("A","Prediction error by architecture")
 axis(1,1:3,c(2,4,8))
 abline(h=0,col="#98A1A7",lwd=1)
 for (j in seq_along(responses)) {
  b <- p[p$response==responses[j],]
  x <- match(b$order,c(2,4,8))+.19*(j-2L)
  points(x,b$error,pch=21,col=colors[j],bg=adjustcolor(colors[j],.55),cex=.7)
  for (k in unique(x)) segments(k-.06,median(b$error[x==k]),k+.06,median(b$error[x==k]),col=colors[j],lwd=2)
 }
 legend("bottomleft",unname(label),pch=21,pt.bg=colors,col=colors,cex=.73,bty="n")

 b <- subset(surf,scenario=="contrasting" & order==2 & supported)
 q <- quantile(b$position,c(.34,.66))
 b$context <- ifelse(b$position<q[1],"Lower phosphorus",
                     ifelse(b$position>q[2],"Higher phosphorus",NA_character_))
 b$bin <- cut(b$spread,breaks=seq(0,max(b$spread)+1e-7,length.out=10),labels=FALSE,include.lowest=TRUE)
 d <- aggregate(cbind(spread,heterogeneity,oracle_heterogeneity)~context+bin+model,b,mean)
 d <- d[order(d$bin),]
 lim <- range(d[c("heterogeneity","oracle_heterogeneity")])
 plot(NA,xlim=range(d$spread),ylim=lim+c(-1,1)*.09*diff(lim),xlab="Mean pairwise phosphorus difference",
      ylab="Heterogeneity dissimilarity",bty="l",las=1)
 tag("B","Context-dependent response")
 envcols <- c("Lower phosphorus"="#216A8A","Higher phosphorus"="#AD5B38")
 for (nm in names(envcols)) {
  z <- d[d$context==nm,]
  y <- z[z$model=="additive",]
  v <- z[z$model=="interaction",]
  lines(y$spread,y$heterogeneity,col=envcols[nm],lwd=1.5,lty=2)
  lines(v$spread,v$heterogeneity,col=envcols[nm],lwd=1.5,lty=1)
  points(v$spread,v$oracle_heterogeneity,pch=16,col=envcols[nm],cex=.61)
 }
 legend("bottomright",names(envcols),col=envcols,lwd=2,cex=.70,bty="n")
 legend("topleft",c("Interaction","Additive","Reference mean"),lty=c(1,2,0),
        pch=c(NA,NA,16),col="#43474A",cex=.70,bty="n")

 q <- subset(u,scenario=="occurrence" & mode=="baseline" & assessment=="surface" &
    estimand %in% paste0("secondary_",responses))
 q$response <- sub("^secondary_","",q$estimand)
 x <- match(q$response,responses)+ifelse(q$order==2,-.12,.12)
 plot(NA,xlim=c(.55,3.45),ylim=c(0,1.05),xaxt="n",xlab="Secondary predictor contribution",
      ylab="95% interval coverage",bty="l",las=1)
 tag("C","Low-order uncertainty")
 axis(1,1:3,unname(label))
 abline(h=.95,lty=2,col="#8C969C")
 segments(x,q$percentile95,x,q$basic95,col=colors[q$response],lwd=1.3)
 points(x,q$percentile95,pch=16,col=colors[q$response],cex=1)
 points(x,q$basic95,pch=1,col=colors[q$response],cex=1)
 if (nrow(q)) text(x, pmax(.03,pmin(.98,pmax(q$basic95,q$percentile95)+.055)),
                     labels=q$order,col="#494D52",cex=.66)
 legend("bottomright",c("Percentile","Basic","Numbers denote order"),pch=c(16,1,NA),cex=.74,bty="n")

 e <- subset(emp,model=="additive" & architecture=="coherent")
 names_order <- c("Finnish birds","Forest Hemiptera","Pond diatoms","Arctic plants","Atlantic marine fish")
 e$j <- match(e$community,names_order)
 e$y <- 11-2*e$j+ifelse(e$order==2,.5,-.5)
 e$x <- match(e$component,responses)
 par(mar=c(4.1,6.55,2.1,.45))
 plot(NA,xlim=c(.49,3.51),ylim=c(-.10,10.55),xaxt="n",yaxt="n",xlab="Dissimilarity component",
      ylab="",bty="n",las=1)
 tag("D","Spatial transfer on manuscript datasets")
 axis(1,1:3,unname(label),tick=FALSE)
 axis(2,seq(9,1,by=-2),names_order,tick=FALSE,las=1,cex.axis=.74)
 z <- ifelse(e$pooled_D2>0,ifelse(e$pooled_D2>.2,"#ACD1CB","#D6E8E6"),
   ifelse(e$pooled_D2< -.2,"#ECC0AE","#F7E4DC"))
 rect(e$x-.49,e$y-.48,e$x+.49,e$y+.48,col=z,border="white")
 text(e$x,e$y,sprintf("%.2f",e$pooled_D2),cex=.72,col="#172E38")
 text(3.5,10.32,"Order 2 above each later order",adj=1,cex=.66,col="#485059")
}
write_figure <- function(ext) {
 path <- file.path(outdir,paste0("figure3_component_validation.",ext))
 if(ext=="pdf") pdf(path,width=10,height=7.2,useDingbats=FALSE) else svg(path,width=10,height=7.2)
 figure()
 invisible(dev.off())
}
write_figure("pdf")
write_figure("svg")
