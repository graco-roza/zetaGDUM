source("docs/experimental/spatiotemporal/phase4_swedish_lakes/scripts/common.R")
fig <- Sys.getenv("ST4_FIGURES",file.path(root,"figures"))
dir.create(fig,recursive=TRUE,showWarnings=FALSE)
a <- fread(result("process_summary.csv"))
order_info <- fread(result("informative_orders.csv"))
zeta <- fread(result("zeta_surfaces.csv"))
site <- fread(result("process_combinations.csv.gz"))
trend <- fread(result("B_trajectory.csv"))
bench <- fread(result("pairwise_tatsumi.csv"))
geo <- fread(processed("lake_geography.csv"))
def <- panels()
ink <- c(expansion="#006DA7",contraction="#D25B22",redistribution="#00836D",numerical="#765298")
communities <- c(phytoplankton="Phytoplankton",littoral="Littoral invertebrates",profundal="Profundal invertebrates")
metrics <- c(presence="Incidence",abundance="Relative abundance",dominance="Local dominance")
processes <- st3_processes
layout_ <- function(n=c(1,1),mar=c(3.9,4.2,3.1,.9),outer=c(.7,.4,3.1,.4)) {
  par(mfrow=n,mar=mar,oma=outer,mgp=c(2.5,.6,0),bty="l",las=1,
      col.main="#182D3F",col.lab="#263A49",col.axis="#263A49")
}
title_ <- function(title,note=NULL) {
  mtext(title,side=3,outer=TRUE,line=1.35,font=2,cex=1.18,col="#132D42")
  if(!is.null(note))mtext(note,side=3,outer=TRUE,line=.15,cex=.77,col="#365264")
}
draw <- function(name,code) {
  png(file.path(fig,paste0(name,".png")),width=12,height=8,units="in",res=180,type="cairo")
  code()
  dev.off()
  code()
}
legend_ <- function(where="topright", net=FALSE,cex=.7) {
  legend(where,c(proc_labels,if(net)"Net"),pch=c(rep(15,4),if(net)18),col=c(ink,if(net)"#132D42"),bty="n",cex=cex)
}
heat <- function(x,field,labs,lim=NULL,diverge=FALSE,digits=2L,marks=NULL) {
  ks <- sort(unique(x$spatial_order)); ms <- sort(unique(x$temporal_order))
  v <- matrix(NA_real_,length(ks),length(ms))
  at <- cbind(match(x$spatial_order,ks),match(x$temporal_order,ms))
  v[at] <- x[[field]]
  pal <- if(diverge) colorRampPalette(c("#175E93","#F6F7F7","#A23D29"))(101) else
    colorRampPalette(c("#EFF5FA","#86B6D0","#07527A"))(101)
  if(is.null(lim))lim<-range(v,na.rm=TRUE)
  if(!all(is.finite(lim))||diff(lim)==0)lim<-c(0,1)
  plot(NA,xlim=c(.5,length(ks)+.5),ylim=c(.5,length(ms)+.5),axes=FALSE,
       xlab="Spatial order",ylab="Temporal order",main=labs,xaxs="i",yaxs="i",cex.main=.91)
  axis(1,seq_along(ks),ks,tick=FALSE)
  axis(2,seq_along(ms),ms,tick=FALSE)
  for(i in seq_along(ks))for(j in seq_along(ms)){
    q<-v[i,j]
    r<-if(is.finite(q)) max(0,min(1,(q-lim[1])/diff(lim))) else NA_real_
    col<-if(is.na(r))"#ECEBE6" else pal[1+floor(100*r)]
    rect(i-.5,j-.5,i+.5,j+.5,col=col,border="white")
    show<-if(is.na(q))"NA" else if(digits<0)formatC(q,format="g",digits=2) else formatC(q,format="f",digits=digits)
    if(!is.null(marks)) {
      k<-which(at[,1]==i & at[,2]==j)
      if(length(k)&&isTRUE(marks[k]))show<-paste0(show,"*")
    }
    text(i,j,show,cex=.70,col=if(!is.na(r)&& (r>.76 || diverge&&r<.18))"white" else "#162B3D")
  }
}
cairo_pdf(file.path(fig,"swedish-lakes-atlas.pdf"),width=12,height=8,onefile=TRUE,pointsize=12)
pdf_device <- dev.cur()
draw("01-monitoring",function(){
  layout_(c(2,2))
  plot(NA,xlim=c(10,25),ylim=c(55,69),asp=1.85,xlab="Longitude (°E)",ylab="Latitude (°N)")
  maps::map("world","Sweden",fill=TRUE,col="#E7EFF2",lwd=.5,border="#8397A4",add=TRUE)
  count<-table(def$membership$id[def$membership$design=="A"&def$membership$included])
  c0<-as.integer(count[geo$EU_ID]);c0[is.na(c0)]<-0L
  pal<-c("#A8B3BB","#74AFC4","#397F9E","#084F77")
  points(geo$longitude,geo$latitude,pch=21,bg=pal[c0+1],col="white",cex=.89,lwd=.5)
  legend("topleft",legend=as.character(0:3),title="Matched assemblages",pch=21,pt.bg=pal,col="#536C7C",bty="n",cex=.72)
  title("A   Published lake coordinates")
  for(g in assemblages){
    m<-def$coverage[[g]]
    i<-order(geo$latitude[match(rownames(m),geo$EU_ID)])
    m<-m[i,,drop=FALSE]
    image(years,seq_len(nrow(m)),t(m*1),col=c("#EBE9E2","#2278A4"),
          xlab="Year",ylab="Lakes ordered south to north",main=communities[g],useRaster=TRUE)
    abline(v=c(1994.5,1999.5,2017.5),lty=2,col="#273746")
  }
  title_("Swedish monitoring lakes, 1992–2022","Blue cells identify observed nonempty assemblages. Coverage defines the matched panels.")
})
for(field in c("zeta","spatial_retention","temporal_retention")) {
  number<-match(field,c("zeta","spatial_retention","temporal_retention"))
  draw(c("02a-zeta-decline","02b-spatial-retention","02c-temporal-retention")[number],local({
    f<-field; no<-number
    function() {
      layout_(c(3,3),mar=c(3.4,3.75,2.15,.85))
      for(g in assemblages)for(type in names(metrics)){
        x<-zeta[zeta$assemblage==g &zeta$taxonomy=="harmonised"&zeta$period=="late"&zeta$window=="consecutive"&
               zeta$response==type&(is.na(zeta$threshold)|zeta$threshold==.8),]
        heat(x,f,paste(communities[g],metrics[type],sep=" / "),if(no==1)NULL else c(0,1),
             digits=if(no==1&&type=="abundance") -1 else if(no==1)1 else 2)
      }
      title_(c("Zeta decline across two orders","Retention when adding a lake","Retention when adding a survey")[no],
             "Late window, 2018–2022, consecutive survey subsets. Raw zeta retains its metric-specific units.")
    }
  }))
}
draw("03-trajectories",function(){
  layout_(c(3,1),mar=c(3,4.4,1.9,1.1))
  cols<-c(abundance="#172D3F",occupancy=unname(ink[1]),numerical=unname(ink[4]))
  ks<-c(2L,4L,6L)
  for(g in assemblages){
    d<-trend[trend$assemblage==g&trend$k%in%ks&trend$metric%in%names(cols),]
    plot(range(d$year),c(0,1),type="n",xlab="Year",ylab="Dissimilarity",main=communities[g])
    for(f in names(cols))for(k in ks){
      s<-d[d$metric==f&d$k==k,];s<-s[order(s$year),]
      lines(s$year,s$value,col=cols[f],lty=match(k,ks),lwd=1.8)
    }
    if(g==assemblages[1]){
      legend("topright",c("Total","Occupancy-related","Numerical"),col=cols,lwd=2,bty="n",horiz=TRUE,cex=.72)
      legend("bottomright",paste("Spatial order",ks),lty=1:3,lwd=1.5,bty="n",horiz=TRUE,cex=.7)
    }
  }
  title_("Long-term spatial dissimilarity on fixed lake panels","Total relative dissimilarity and its two additive components, 1995–2022.")
})
signed <- function(d,maxval) {
  labels<-c("Incidence","Relative abundance","Dominance")
  plot(NA,xlim=c(-maxval,maxval),ylim=c(.4,3.8),yaxt="n",xlab="Diversity units",ylab="")
  axis(2,3:1,labels,tick=FALSE,cex.axis=.69)
  abline(v=0,col="#8E9AA3")
  for(j in seq_along(labels)){
    x<-d[d$metric==c("incidence","abundance","dominance_dissimilarity")[j],]
    x<-x[match(processes,x$process),]
    left<-right<-0
    for(i in 1:4){
      rect(left-x$negative[i],4-j-.25,left,4-j+.25,col=ink[i],border=NA)
      rect(right,4-j-.25,right+x$positive[i],4-j+.25,col=ink[i],border=NA)
      left<-left-x$negative[i];right<-right+x$positive[i]
    }
    points(x$delta[1],4-j,pch=18,cex=1.2,col="#132D42")
  }
}
draw("04-processes",function(){
  layout_(c(3,3),mar=c(3.1,7.25,2.0,.9),outer=c(1.15,.4,3.25,.4))
  d<-a[a$temporal_order==1&a$spatial_order%in%c(2,4,6)&a$metric%in%c("incidence","abundance","dominance_dissimilarity"),]
  s <- aggregate(cbind(negative,positive)~assemblage+spatial_order+metric,d,sum)
  margin<-max(c(s$positive,s$negative))*1.1
  for(g in assemblages)for(k in c(2,4,6)){
    signed(d[d$assemblage==g&d$spatial_order==k,],margin)
    title(paste(communities[g]," / k =",k),cex.main=.86)
  }
  title_("Homogenisation and differentiation by process","Order one in time. Negative stacks represent homogenisation, and the diamond is net change.")
  mtext("Blue  Expansion     Orange  Contraction     Green  Redistribution     Purple  Relative shares",
         1,outer=TRUE,line=.1,cex=.77,col="#334E60")
})
draw("05-percentages",function(){
  layout_(c(2,3),mar=c(3.4,3.7,2.0,.8))
  d<-a[a$temporal_order==1&a$spatial_order>=2&a$metric=="abundance",]
  for(side in c("negative_share","positive_share"))for(g in assemblages){
    x<-d[d$assemblage==g,]
    v<-xtabs(reformulate(c("process","spatial_order"),side),x)[processes,,drop=FALSE]
    barplot(v,col=ink,border=NA,ylim=c(0,100),main=communities[g],xlab="Spatial order",ylab="Percent within sign")
    mtext(if(side=="negative_share")"Homogenisation" else "Differentiation",3,line=-.5,cex=.74)
  }
  title_("Relative process contributions to abundance dissimilarity","Each sign uses its own gross magnitude as the denominator. Colours match Figure 4.")
})
draw("06-crossed",function(){
  layout_(c(2,3),mar=c(4.25,6.0,2.2,.75))
  d<-a[a$spatial_order==2&a$temporal_order%in%c(1,5)&a$metric%in%c("occupancy","numerical"),]
  maxv<-max(abs(d$signed_mean),na.rm=TRUE)
  pal<-colorRampPalette(c("#1B6192","#FAF9F5","#A1402B"))(101)
  for(m in c(1,5))for(g in assemblages){
    p<-d[d$temporal_order==m&d$assemblage==g,]
    v<-xtabs(signed_mean~process+metric,p)[processes,c("occupancy","numerical")]
    plot(NA,xlim=c(.5,2.5),ylim=c(.5,4.5),axes=FALSE,xlab="",ylab="",main=paste(communities[g]," / m =",m))
    axis(1,1:2,c("Occupancy","Numerical"),tick=FALSE)
    axis(2,4:1,proc_labels,tick=FALSE,cex.axis=.74)
    for(i in 1:4)for(j in 1:2){
      col<-pal[1+round(50*(1+v[i,j]/maxv))]
      rect(j-.5,4.5-i,j+.5,5.5-i,col=col,border="white")
      text(j,5-i,sprintf("%+.3f",v[i,j]),col=if(abs(v[i,j])>.62*maxv)"white" else "#142C3C")
    }
  }
  title_("Crossed process-by-diversity decomposition","Four processes × two components. All entries are signed dissimilarity units at spatial order two.")
})
draw("07-hidden-reorganisation",function(){
  layout_(c(1,3),mar=c(4.1,4.8,2.1,.8))
  d<-site[site$k==2&site$m==1,]
  d<-aggregate(cbind(before,after)~assemblage+offset+metric,d,mean)
  for(g in assemblages){
    p<-d[d$assemblage==g,]
    x<-p$after[p$metric=="abundance"]-p$before[p$metric=="abundance"]
    plot(NA,xlim=c(-.2,.2),ylim=c(0,1),xlab="Net spatial relative dissimilarity change",
         ylab="Temporal distance",main=communities[g],cex.lab=.88)
    abline(v=0,lty=2,col="#8EA1AC")
    for(i in 1:3){
      fld<-c("site_temporal_abundance","temporal_abundance","site_temporal_dominance")[i]
      rr<-p[p$metric==fld,];rr<-rr[order(rr$offset),]
      points(x,rr$after,pch=c(16,17,1)[i],col=c("#0075AB","#79529A","#AD4F2B")[i],cex=1.2)
    }
    if(g==assemblages[1]) legend("bottomleft",c("Local relative Bray–Curtis","Shared-minima Bray–Curtis","Local dominance Sørensen"),
       bty="n",pch=c(16,17,1),col=c("#0075AB","#79529A","#AD4F2B"),cex=.72)
  }
  title_("Temporal reorganisation can remain high near zero spatial change","Each point represents one matched survey pair on the long-term lake panel.")
})
draw("08-spatiotemporal-orders",function(){
  layout_(c(3,4),mar=c(3.35,3.5,2.15,.7))
  d<-a[a$metric=="abundance"&a$spatial_order>=2,]
  lim<-max(abs(d$signed_mean),na.rm=TRUE)
  for(g in assemblages)for(p in processes){
    x<-d[d$assemblage==g&d$process==p,]
    f<-order_info[order_info$assemblage==g&order_info$response=="presence",]
    take<-match(paste(x$spatial_order,x$temporal_order),paste(f$spatial_order,f$temporal_order))
    x$signed_mean[x$valid_fraction<.8]<-NA
    heat(x,"signed_mean",paste(communities[g],proc_labels[match(p,processes)],sep=" / "),c(-lim,lim),TRUE,2,
         marks=!f$informative_shared[take])
  }
  title_("Signed abundance contributions across both orders","NA marks insufficient decomposition support. * flags weak shared temporal support or taxon information.")
})
draw("09-pairwise-benchmark",function(){
  layout_(c(2,3),mar=c(4.0,4.55,2.1,.8))
  for(type in c("incidence","abundance"))for(g in assemblages){
    s<-bench[bench$metric==type&bench$assemblage==g&bench$tatsumi_defined,]
    n0<-s$tatsumi_loss_homogenisation+s$tatsumi_gain_homogenisation
    p0<-s$tatsumi_loss_differentiation+s$tatsumi_gain_differentiation
    lim<-range(0,n0,p0,s$shapley_homogenisation,s$shapley_differentiation,na.rm=TRUE)
    plot(n0,s$shapley_homogenisation,xlim=lim,ylim=lim,asp=1,pch=16,cex=.48,col=adjustcolor("#00699E",.46),
         xlab="Tatsumi gross magnitude",ylab="Shapley gross magnitude",main=paste(communities[g],type,sep=" / "))
    points(p0,s$shapley_differentiation,pch=1,col=adjustcolor("#B5542D",.62),cex=.55)
    abline(0,1,lty=2,col="#899CA8")
    if(g==assemblages[1])legend("topleft",c("Homogenising","Differentiating"),pch=c(16,1),
      col=c("#00699E","#B5542D"),bty="n",cex=.67)
  }
  title_("Pairwise accounting agrees on net change but can allocate different magnitudes",
         "Same observed lake pairs and relative compositions. Only defined published components are shown.")
})
dev.off(pdf_device)
