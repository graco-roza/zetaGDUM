root <- "docs/experimental/spatiotemporal"
source(file.path(root,"phase3_load.R"))
st3_load(root)
out <- Sys.getenv("ST3_OUTPUT",file.path(root,"phase3/results"))
figs <- Sys.getenv("ST3_FIGURES",file.path(root,"phase3/figures"))
dir.create(figs,recursive=TRUE,showWarnings=FALSE)
read <- function(x) read.csv(file.path(out,paste0(x,".csv")))
ink <- c(expansion="#0072B2",contraction="#B94818",redistribution="#00816A",numerical="#75589B",total="#172531")
process <- st3_processes
pretty_process <- c(expansion="Expansion",contraction="Contraction",redistribution="Redistribution",numerical="Numerical",total="Net")
title_color <- "#172D3C"
main <- read("main_panels")
ord <- read("order_processes")
supp <- read("order_support")
short_cases <- c(expansion="Expansion",contraction="Contraction",redistribution="Redistribution",
  numerical="Numerical only",expansion_numerical="Expansion + numerical",contraction_numerical="Contraction + numerical",
  opposing="Opposing processes",filtering="Spatial filtering",widespread="Widespread species",
  restricted="Restricted species",persistence="Persistence change",null="Null",shared_shock="Shared shocks")

draw <- function(number,name,width,height,fn) {
  svg(file.path(figs,sprintf("%02d-%s.svg",number,name)),width,height,family="sans",pointsize=12)
  on.exit(dev.off())
  fn()
}
layout_base <- function(rooms=c(1,1),mar=c(4.2,4.5,3.5,1.0),outer=c(.2,0,2.1,0)) {
  par(mfrow=rooms,mar=mar,oma=outer,mgp=c(2.6,.67,0),bty="l",las=1,
      col.axis="#3B4C57",col.lab=title_color,col.main=title_color)
}
heading <- function(x) mtext(x,3,outer=TRUE,line=.6,font=2,cex=1.12,col=title_color)
heat <- function(tab,field,title,label,lim=NULL,digits=2L) {
  plot_st_surface(tab,field,label,digits,limits=lim)
  mtext(title,3,line=1.56,font=2,cex=.92,col=title_color)
}
legend_process <- function(pos="topright",show_net=FALSE) {
  labs <- pretty_process[if (show_net) c(process,"total") else process]
  cols <- ink[names(labs)]
  legend(pos,unname(labs),bty="n",col=unname(cols),pch=if(show_net)c(15,15,15,15,18) else 15,
         cex=.78)
}
draw(1,"framework",12.3,5.7,function() {
  layout_base(mar=c(.1,.15,1.1,.15))
  plot(0:1,0:1,type="n",xlim=c(0,1),ylim=c(0,1),axes=FALSE,xlab="",ylab="",xaxs="i",yaxs="i")
  box <- function(x,y,w,h,fill,label,secondary=NULL,textcolor=title_color) {
    rect(x-w/2,y-h/2,x+w/2,y+h/2,col=fill,border="#B5C4CE",lwd=1.15)
    text(x,y+if(is.null(secondary))0 else .021,label,col=textcolor,font=2,cex=.95)
    if (!is.null(secondary)) text(x,y-.028,secondary,col=textcolor,cex=.71)
  }
  arrow <- function(x,y,x1,y1) arrows(x,y,x1,y1,length=.09,col="#526975",lwd=1.7)
  box(.115,.525,.188,.16,"#EDF3F7","Matched communities","Absolute abundances at all sites")
  box(.515,.525,.205,.19,"#EDF3F7","Four-block Shapley","Relative abundance recalculated",
      title_color)
  y <- c(.82,.624,.428,.232)
  tint <- c("#D7EAF5","#F8E3D7","#D9EFE9","#EBE1F6")
  for (i in seq_along(process)) {
    box(.318,y[i],.16,.126,tint[i],pretty_process[process[i]])
    arrow(.21,.525,.236,y[i])
    arrow(.401,y[i],.409,.525)
  }
  box(.811,.73,.303,.16,"#EDF3F7","Change in differentiation","Signed change in spatial or joint diversity")
  box(.811,.49,.303,.16,"#EDF3F7","Temporal reorganisation","Shared assemblages and matched sites")
  box(.811,.25,.303,.16,"#EDF3F7","Crossed abundance change","Four processes × two diversity components")
  for (h in c(.73,.49,.25)) arrow(.622,.525,.657,h)
  text(.512,.065,"Environmental models use one combination × contrast row and retain signed diversity units.",
       col="#314650",cex=.91)
  heading("Community-level accounting of compositional change")
})
draw(2,"abundance-components",11.5,7.8,function() {
  layout_base(c(2,2))
  z <- aggregate(cbind(before,after)~spatial_order+temporal_order+response,
     subset(ord,scenario=="filtering"&response %in% c("occupancy","numerical")&process=="expansion"),mean)
  lim <- range(c(z$before,z$after))
  for (status in c("before","after")) for (response in c("occupancy","numerical")) {
    a <- z[z$response==response,]
    heat(a,status,paste(if(status=="before")"Initial" else "Final",response),
         "Abundance dissimilarity units",lim,3)
  }
  heading("Incomplete occupancy and positive numerical heterogeneity")
})
draw(3,"signed-processes",12,7.4,function() {
  layout_base(c(1,2),c(4.35,10.0,3,1.05))
  z <- aggregate(observed_mean~scenario+process,main,mean)
  a <- xtabs(observed_mean~scenario+process,z)
  lim <- range(c(rowSums(pmax(a[,process],0)),rowSums(pmin(a[,process],0))))*1.14
  draw_cases <- function(cases,lab,legend=FALSE) {
    y <- rev(seq_along(cases))
    plot(NA,xlim=lim,ylim=c(.25,length(y)+.95),yaxt="n",xlab="Change in abundance dissimilarity",ylab="",main=lab)
    axis(2,y,short_cases[cases],tick=FALSE,las=1,cex.axis=.84)
    abline(v=0,col="#93A3AF",lwd=1)
    for (i in seq_along(cases)) {
      val <- a[cases[i],process]
      left <- right <- 0
      for (j in seq_along(process)) {
        from <- if(val[j]>=0) right else left
        rect(min(from,from+val[j]),y[i]-.23,max(from,from+val[j]),y[i]+.23,col=ink[process[j]],border=NA)
        if(val[j]>=0) right <- right+val[j] else left <- left+val[j]
      }
      points(a[cases[i],"total"],y[i],pch=18,cex=1.15,col=ink["total"])
    }
    if(legend) legend_process("bottomright",TRUE)
  }
  draw_cases(names(short_cases)[1:7],"A   Isolated and opposing changes")
  draw_cases(names(short_cases)[8:13],"B   Environmental and null cases",TRUE)
  heading("Signed process contributions at spatial order two")
})
draw(4,"crossed-processes",11.6,7.7,function() {
  layout_base(c(2,2),c(3.9,9.4,3,1))
  z <- aggregate(mean~spatial_order+temporal_order+response+process,
     subset(ord,scenario=="filtering"&response %in% c("occupancy","numerical")),mean)
  lim <- max(abs(z$mean))
  palette <- colorRampPalette(c("#215D96","#F8F9FA","#AE3B13"))(101)
  for (j in 1:4) {
    k <- c(2L,6L,2L,6L)[j]
    m <- c(1L,1L,4L,4L)[j]
    a <- z[z$spatial_order==k&z$temporal_order==m,]
    amount <- xtabs(mean~process+response,a)[process,c("occupancy","numerical")]
    plot(NA,xlim=c(.5,2.5),ylim=c(.5,4.5),axes=FALSE,xlab="",ylab="",
         main=sprintf("Spatial order %s  |  Temporal order %s",k,m),xaxs="i",yaxs="i")
    axis(1,1:2,c("Occupancy dissimilarity","Numerical dissimilarity"),tick=FALSE,cex.axis=.83)
    axis(2,4:1,pretty_process[process],tick=FALSE,cex.axis=.8)
    for (c in 1:2) for (r in 1:4) {
      v <- amount[r,c]
      jcol <- round(50*(1+v/lim))+1
      rect(c-.5,4.5-r,c+.5,5.5-r,col=palette[jcol],border="white")
      text(c,5-r,sprintf("%+.3f",v),cex=.91,col=if(abs(v)>.62*lim)"white" else title_color)
    }
  }
  heading("Crossed allocation of abundance dissimilarity change")
})
draw(5,"opposing-processes",12,5.5,function() {
  layout_base(c(1,2))
  q <- subset(main,scenario=="opposing")
  a <- reshape(q[c("replicate","process","observed_mean")],idvar="replicate",timevar="process",direction="wide")
  x <- a$observed_mean.expansion
  y <- a$observed_mean.numerical
  limits <- range(c(-x,y))
  plot(-x,y,asp=1,xlim=limits,ylim=limits,pch=16,cex=.64,col=adjustcolor(ink["numerical"],.52),
       xlab="Expansion homogenisation magnitude",ylab="Numerical differentiation magnitude",
       main="A   Equal and opposing contributions")
  abline(0,1,lty=2,col=title_color)
  means <- c(mean(x),mean(y),mean(a$observed_mean.total))
  lim <- max(abs(c(x,y)))*1.1
  plot(NA,xlim=c(.6,3.4),ylim=c(-lim,lim),xaxt="n",xlab="",ylab="Change in abundance dissimilarity",
       main="B   Across 200 independent panels")
  axis(1,1:3,c("Expansion","Numerical","Net change"),tick=FALSE)
  abline(h=0,col="#A5B2BA",lty=2)
  parts <- list(x,y,a$observed_mean.total)
  for(j in 1:3) {
    ci <- quantile(parts[[j]],c(.025,.975))
    segments(j,ci[1],j,ci[2],lwd=2,col=ink[c("expansion","numerical","total")][j])
    points(j,means[j],pch=c(16,16,18)[j],cex=1.35,col=ink[c("expansion","numerical","total")][j])
  }
  heading("Compositional change beneath a zero net trend")
})
draw(6,"environmental-curves",11.8,7.8,function() {
  layout_base(c(2,2),outer=c(1.35,0,2.1,0))
  z <- read("environmental_curves")
  y <- aggregate(estimate~process+delta,z,function(x)c(mean=mean(x),lo=quantile(x,.1,names=FALSE),hi=quantile(x,.9,names=FALSE)))
  coverage <- aggregate(extrapolation~process+delta,z,mean)
  ylim <- range(y$estimate[,c("lo","hi")])
  for(p in process) {
    a <- y[y$process==p,]
    b <- coverage[coverage$process==p,]
    x <- a$delta
    v <- a$estimate
    plot(x,v[,"mean"],type="n",ylim=ylim,xlab="Mean site climate change (training SD)",
         ylab="Predicted contribution",main=pretty_process[p])
    polygon(c(x,rev(x)),c(v[,"lo"],rev(v[,"hi"])),col=adjustcolor(ink[p],.18),border=NA)
    abline(h=0,col="#9AAAB5",lty=2)
    lines(x,v[,"mean"],col=ink[p],lwd=2)
    if(any(b$extrapolation>.5)) points(x[b$extrapolation>.5],v[b$extrapolation>.5,"mean"],pch=1,cex=.65,col=ink[p])
  }
  heading("Signed environmental response functions and panel variability")
  mtext("Bands show the 10th–90th percentiles across panels. Open circles mark extrapolation in most fits.",
        1,outer=TRUE,line=.25,cex=.78,col="#3B4C57")
})
draw(7,"validation",12.5,8.0,function() {
  layout_base(c(2,2),c(4.5,4.7,3.15,.9))
  t <- read("independent_truth")
  x <- aggregate(estimate~scenario+process,main,mean)
  z <- merge(x,t,by=c("scenario","process"))
  keep <- z$process %in% process & abs(z$truth)>1e-6
  lim <- range(c(z$truth[keep],z$estimate[keep]))
  plot(z$truth[keep],z$estimate[keep],type="n",xlim=lim,ylim=lim,asp=1,xlab="Independent Monte Carlo truth",
       ylab="Mean fitted contrast",main="A   Environmental contrast recovery")
  abline(0,1,lty=2,col="#8B9BA9")
  for (j in seq_along(process)) {
    q <- keep & z$process==process[j]
    points(z$truth[q],z$estimate[q],pch=j+14L,col=ink[process[j]])
  }
  legend_process("topleft")
  h <- read("holdout_panels")
  h$improvement <- h$baseline_rmse-h$rmse
  a <- aggregate(improvement~split+process,subset(h,spatial_order==2&temporal_order==1),mean)
  m <- xtabs(improvement~split+process,a)
  m <- m[c("sites","forward"),c(process,"total")]
  mids <- barplot(m,beside=TRUE,col=c("#7AAAC9","#5A747F"),ylim=range(c(-.003,.011,m)),names.arg=rep("",5),
     main="B   Held-out gain over intercept only",ylab="Decrease in mean RMSE")
  axis(1,colMeans(mids),pretty_process[c(process,"total")],tick=FALSE,cex.axis=.71)
  legend("topleft",c("Complete site histories","Forward surveys"),fill=c("#7AAAC9","#5A747F"),bty="n",cex=.76)
  abline(h=0,col="#243746")
  a <- read("interval_summary")
  for (scenario in c("null","shared_shock")) {
    q <- a[a$scenario==scenario & a$process=="total",]
    method <- c("spatial_block","temporal_block","crossed")
    q <- q[match(method,q$method),]
    labs <- c("Spatial block","Temporal block","Crossed")
    se <- q$coverage_mc_se
    plot(1:3,q$coverage,type="n",ylim=c(.65,1.02),xaxt="n",xlab="",ylab="Coverage of nominal 95% intervals",
         main=if(scenario=="null")"C   Clean null" else "D   Shared shocks")
    axis(1,1:3,labs,tick=FALSE,cex.axis=.85)
    segments(1:3,q$coverage-1.96*se,1:3,pmin(1,q$coverage+1.96*se),col="#576E7D",lwd=2)
    points(1:3,q$coverage,pch=16,cex=1.22,col=ink[if(scenario=="null")"expansion" else "contraction"])
    abline(h=.95,lty=2,col="#81929F")
  }
  heading("Recovery and calibration in independent validation panels")
})
draw(8,"order-selection",11.8,7.7,function() {
  layout_base(c(2,2))
  z <- aggregate(cbind(effective_df,rank,zero_share)~spatial_order+temporal_order,subset(supp,scenario=="filtering"),mean)
  stress <- aggregate(zero_share~spatial_order+temporal_order,subset(supp,scenario=="support_stress"),mean)
  a <- aggregate(common_support~spatial_order+temporal_order,
                 subset(ord,scenario=="support_stress"&response=="temporal_presence"&process=="expansion"),mean)
  b <- aggregate(crossed_contrast_se~spatial_order+temporal_order,subset(supp,scenario=="filtering"),mean)
  heat(z,"effective_df","A   Effective design degrees of freedom","Penalized ridge geometry",digits=1)
  heat(a,"common_support","B   Temporal incidence support","All subsets defined in the site-specific anchor scenario",c(0,1),2)
  heat(stress,"zero_share","C   Initial zero sharing","Initial shared vectors empty in the site-specific anchor scenario",c(0,1),2)
  z <- merge(z,b,all.x=TRUE)
  heat(z,"crossed_contrast_se","D   Crossed-block uncertainty","Abundance contrast SE where at least four windows exist",digits=3)
  heading("Information and support across spatial and temporal orders")
})
