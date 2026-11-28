st5_simulate <- function(seed,n=60L,species=200L,times=30L,
    scenario=c("multiple","expansion","contraction","redistribution","numerical","opposing",
               "null","stable","reversal","irregular"), effect=1) {
  scenario <- match.arg(scenario)
  effect <- rep_len(effect,2L)
  set.seed(seed)
  xy <- matrix(runif(n*2),n,2)
  time <- cumsum(c(0,if(scenario=="irregular") sample(c(.35,1,2.4),times-1,TRUE) else rep(1,times-1)))
  ar <- function() as.numeric(stats::filter(rnorm(times,0,.45),.7,method="recursive"))
  E <- array(0,c(n,times,4),dimnames=list(NULL,NULL,c("E1","E2","noise1","noise2")))
  for(j in 1:4) {
    spatial <- if(j<=2) xy[,j] else rnorm(n)
    E[,,j] <- outer(spatial,ar(),"+")+matrix(rnorm(n*times,0,.35),n)
    E[,,j] <- scale(as.vector(E[,,j]))
  }
  records <- expand.grid(site=as.character(seq_len(n)),time=time,stringsAsFactors=FALSE)
  A <- array(0,c(n,times,species))
  groups <- rep(1:5,length.out=species)
  prev <- matrix(runif(n*species)<rep(runif(species,.18,.78),each=n),n)
  absolute <- matrix(rlnorm(n*species,1,.6),n)
  A[,1,] <- prev*absolute
  active <- switch(scenario,expansion=1,contraction=2,redistribution=3,numerical=4,
     opposing=c(1,2,4),null=1:4,stable=3,1:4)
  for(t in seq_len(times-1)) {
    # Draw every random variable independently of the active mechanism and current states.
    event <- matrix(runif(n*species),n)
    choice <- matrix(runif(n*species),n)
    shock <- matrix(rnorm(n*species,0,.18),n)
    replacement <- matrix(rlnorm(n*species,1,.6),n)
    dx <- E[,t+1,1]-E[,t,1]
    dy <- E[,t+1,2]-E[,t,2]
    driver <- if(scenario%in%c("null","stable")) c(0,0) else effect
    if(scenario=="reversal" && t>times/2) driver <- -driver
    next_state <- prev
    if(1%in%active) {
      i <- which(groups==1)
      next_state[,i] <- prev[,i] | event[,i] < plogis(-4.65+driver[1]*(.9*E[,t,1]+1.1*dx))
    }
    if(2%in%active) {
      i <- which(groups==2)
      next_state[,i] <- prev[,i] & event[,i] >= plogis(-4.5+driver[2]*(.85*E[,t,2]+dy))
    }
    if(3%in%active) for(j in which(groups==3)) {
      from <- which(prev[,j]); to <- which(!prev[,j])
      r <- min(length(from),length(to))
      n_swap <- if(r) sum(event[seq_len(r),j] <
        plogis(if(scenario=="stable") -1.2 else -2.4 + driver[2]*sd(E[,t,2])*.65)) else 0L
      if(n_swap) {
        next_state[from[order(choice[from,j])[seq_len(n_swap)]],j] <- FALSE
        next_state[to[order(choice[to,j]-driver[2]*.12*dy[to])[seq_len(n_swap)]],j] <- TRUE
      }
    }
    absolute[!prev & next_state] <- replacement[!prev & next_state]
    if(4%in%active) {
      i <- which(groups==4)
      absolute[,i] <- absolute[,i]*exp(pmin(pmax(driver[1]*.45*dx+shock[,i],-5),5))
    }
    if(scenario%in%c("multiple","opposing","null","reversal","irregular")) {
      i <- which(groups==5)
      next_state[,i] <- xor(prev[,i],event[,i]<.012)
      new <- !prev[,i] & next_state[,i]
      a <- absolute[,i]; a[new] <- replacement[,i][new]; absolute[,i] <- a
    }
    prev <- next_state
    A[,t+1,] <- prev*absolute
  }
  a <- matrix(A,n*times,species)
  a <- a/pmax(rowSums(a),.Machine$double.xmin)
  list(records=records,abundance=a,env=matrix(E,n*times,4,dimnames=list(NULL,dimnames(E)[[3]])),
       xy=xy,times=time)
}

st5_missing <- function(sim,rate,scheme=c("random","network","period","environment","cluster"),seed) {
  scheme <- match.arg(scheme)
  if(!is.finite(rate) || rate<0 || rate>=1) stop("Invalid missing fraction.")
  set.seed(seed)
  r <- sim$records
  t <- match(r$time,sim$times)
  site <- as.integer(r$site)
  score <- switch(scheme,random=runif(nrow(r)),
    network={sign <- sample(c(-1,1),nrow(sim$xy),TRUE); sign[site]*(t-mean(t))+runif(nrow(r),0,2)},
    period=runif(length(sim$times))[t],
    environment=sim$env[,1]+rnorm(nrow(r),0,.4),
    cluster=as.numeric(sim$xy[site,1]+.15*sim$xy[site,2])+runif(nrow(r),0,.1))
  if(scheme=="period") {
    remove <- head(order(tapply(score,t,mean),decreasing=TRUE),round(rate*length(sim$times)))
    r$available <- !t%in%remove
  } else {
    remove <- head(order(score,decreasing=TRUE),round(rate*nrow(r)))
    r$available <- !seq_len(nrow(r))%in%remove
  }
  r
}
