source("docs/experimental/spatiotemporal/phase5_model/scripts/common.R")
tests <- list()
verify <- function(label,condition) {
  tests[[label]] <<- isTRUE(condition)
  if(!isTRUE(condition)) stop(label)
}
close <- function(a,b) max(c(0,abs(a-b)),na.rm=TRUE)<1e-10
sim <- st5_simulate(5,12,35,12,"irregular")
r <- st5_missing(sim,.25,"random",33)
p <- st5_input(r,abundance=sim$abundance,env=sim$env)
binary <- st5_input(r,incidence=sim$abundance>0,env=sim$env)
verify("unavailable is absent from index rather than an assemblage",sum(is.na(p$index))==sum(!r$available))
loc <- st5_support(p,"local",budget=Inf)
lr <- st5_responses(p,loc,TRUE)
verify("abundance endpoints do not depend on requesting attribution",
       close(lr$after,st5_responses(p,loc,FALSE)$after))
verify("consecutive observed survey coverage",nrow(lr$rows)==nrow(p$x)-length(unique(p$records$site)))
verify("actual unequal local lags",all(lr$rows$lag==lr$rows$final-lr$rows$initial)&&length(unique(lr$rows$lag))>1)
verify("local directional components",close(lr$rows$gain_component+lr$rows$loss_component,lr$after[,"site_temporal_presence"]))
ch <- st5_support(p,"change",budget=Inf)
a <- st5_responses(p,ch,TRUE)
b <- st5_responses(binary,ch,TRUE)
verify("incidence three-game agrees with four-game dummy",close(a$phi[,"incidence",1:3],b$phi[,"incidence",])&&
          close(a$phi[,"incidence","numerical"],0))
verify("observed common-support efficiency",max(a$efficiency_error,b$efficiency_error,lr$efficiency_error)<1e-10)
verify("spatial additive incidence",close(a$before[,"incidence"],rowSums(a$before[,c("turnover","nestedness")])))
verify("abundance identity",close(a$after[,"abundance"],rowSums(a$after[,c("occupancy","numerical")])))
z <- st5_extract(b,"incidence")
ref <- st5_basis(binary,z)
d <- st5_design(binary,z,ref)
native_rows <- do.call(rbind,z$rows$a)
native <- zetaGDUM::ugdm_design(binary$env,native_rows,ref=list(environment_basis=ref$native))
verify("native position sums",close(d$x[,grep(":P$",colnames(d$x))],native$Xpos))
verify("native spread basis",close(d$x[,grep(":S[123]$",colnames(d$x))],native$Xdir))
pred <- st5_predict(st5_fit(z,d),z,d)
verify("coherent sum",close(rowSums(pred$phi),pred$change))
verify("coherent bounded endpoints",all(pred$endpoint>=0 & pred$endpoint<=1))
crossed <- st5_crossed_data(a)
cd <- st5_design(p,crossed,st5_basis(p,crossed))
cp <- st5_predict(st5_fit(crossed,cd,mode="crossed"),crossed,cd)
verify("crossed nonnegative simplex",close(rowSums(cp$composition),1)&&min(cp$composition)>=0)
verify("crossed rows equal process totals",close(apply(cp$crossed,c(1,2),sum),cp$phi))
verify("crossed columns equal component changes",close(apply(cp$crossed,c(1,3),sum),
      cp$composition[,-1]-crossed$composition[,1:2]))
verify("crossed totals equal endpoint change",close(rowSums(cp$phi),cp$change))
complete <- st5_input(sim$records,incidence=sim$abundance>0,env=sim$env)
chron <- st5_responses(complete,st5_support(complete,"change",budget=Inf),TRUE)
endpoint <- st5_responses(complete,st5_support(complete,"change",windows=list(list(min(sim$times),max(sim$times))),budget=Inf),TRUE)
verify("fixed-support telescoping",close(tapply(chron$after[,"incidence"]-chron$before[,"incidence"],
      rep(seq_len(nrow(endpoint$rows)),length(sim$times)-1),sum),endpoint$after[,"incidence"]-endpoint$before[,"incidence"]))
folds <- st5_folds(z,"forward",blocks=4)
verify("forward folds dates separated",length(folds)>0&&all(vapply(folds,function(f)
  max(z$rows$final[f$train])<min(z$rows$initial[f$test]),logical(1))))
blocks <- spatial_blocks(complete$sites,sim$xy[as.integer(complete$sites),1],sim$xy[as.integer(complete$sites),2])
folds <- st5_folds(z,"spatial",site_blocks=blocks)
verify("spatial fold leakage absent",length(folds)==3&&all(vapply(folds,function(f)
  !length(intersect(unlist(z$rows$sites[f$train]),unlist(z$rows$sites[f$test]))),logical(1))))
v <- st5_validate(binary,z,"spatial",blocks,grid=c(.001,.1))
i <- st5_importance(v)
verify("predictor Shapley additivity",close(i[i$fold=="pooled"&i$response=="joint",]$total_improvement[1],
      sum(i[i$fold=="pooled"&i$response=="joint",]$importance)))
u <- st5_input(data.frame(site=rep(c("a","b"),2),time=rep(1:2,each=2),available=c(1,0,1,1)),
   incidence=rbind(c(0,0),c(NA,NA),c(1,1),c(0,0)))
verify("observed empties allowed",sum(rowSums(u$x)==0)==2)
uu <- st5_responses(u,st5_support(u,"local",budget=Inf),TRUE)
verify("empty baseline common coalition undefined",is.na(uu$before[1,1])&&all(is.na(uu$phi)))
ss <- st5_responses(u,st5_support(u,"spatial",budget=Inf))
cc <- st5_responses(u,st5_support(u,"change",budget=Inf),TRUE)
verify("unsupported matching yields empty response",nrow(cc$rows)==0)
pair_sample <- st5_input(data.frame(site=c("a","b","c","a","b","c"),time=rep(1:2,each=3)),
   incidence=rbind(c(1,1,0),c(1,0,1),c(1,0,0),c(1,0,0),c(1,1,1),c(1,1,0)))
regional <- st5_extract(st5_responses(pair_sample,st5_support(pair_sample,"local",lag=1,budget=Inf),TRUE),
                       "site_temporal_presence")
verify("local decomposition uses the full regional network",
       regional$phi[regional$rows$sites=="a","expansion"]==0 &&
       regional$phi[regional$rows$sites=="a","redistribution"]>0)
verify("one-site temporal response does not use spatial spread",
       !any(st5_design(p,st5_extract(lr,"site_temporal_presence"),ref)$columns$group=="spatial_spread"))
rect <- st5_responses(binary,st5_support(binary,"rectangle",k=2,m=2,budget=Inf))
verify("each rectangle uses four recorded communities",all(lengths(rect$rows$a)==4)&&
       all(unlist(rect$rows$a)%in%seq_len(nrow(binary$x))))
verify("each spatial response uses its observed sites",
       st5_eligibility(binary,st5_support(binary,"spatial",budget=Inf))$eligible_observation_fraction==1)
mix <- st5_support(binary,"change",k=2,m=2,
                  windows=list(list(sim$times[1:2],sim$times[3:4])),budget=Inf)
verify("paired windows keep chronological alignment",all(vapply(mix,function(g) max(g$a)<min(g$b),logical(1))))
verify("native pairwise Sørensen agrees with binary counts",
       close(rect$before[,"incidence"],vapply(rect$rows$a,function(i)
         1-sum(colSums(binary$x[i,])==length(i))/mean(rowSums(binary$x[i,])),numeric(1))))
verify("defined local distance remains available when attribution is undefined",
       st5_extract(uu,"site_temporal_presence")$y==1)
verify("period and importance weights target the stated means",
       close(st5_weights(c(1,1,2),c(1,3,9)),c(.125,.375,.5)))
write_result(data.frame(test=names(tests),passed=unlist(tests)),"mathematical_validation")
