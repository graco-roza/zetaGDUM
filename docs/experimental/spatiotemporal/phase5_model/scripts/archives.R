st5_lake_archive <- function(lake_root) {
  lake_raw <- file.path(lake_root,"data/raw")
  read_tsv <- function(name) read.delim(file.path(lake_raw,name),check.names=FALSE,na.strings=c("NA",""))
  chem <- read_tsv("01_Physochemical.tsv")
  climate <- read_tsv("02_ClimateData.tsv")
  chem_key <- paste(chem[[1]],chem$Year); climate_key <- paste(climate$EU_ID,climate$Year)
  stopifnot(!anyDuplicated(chem_key),!anyDuplicated(climate_key))
  assemblages <- readRDS(file.path(lake_root,"data/processed/assemblages.rds"))
  geography <- read.csv(file.path(lake_root,"data/processed/lake_geography.csv"))
  sites <- unique(geography$EU_ID)
  coord <- as.matrix(geography[match(sites,geography$EU_ID),c("latitude","longitude")])
  # Contiguous latitude bands, with year blocks used separately.
  spatial_lake <- spatial_blocks(sites,coord[,1],coord[,2],3L)
  lake_inputs <- lapply(assemblages,function(a) {
    s <- a$samples
    keep <- is.finite(rowSums(a$harmonised)) & grepl("^observed",s$status)
    key <- paste(s$id[keep],s$year[keep])
    C <- chem[match(key,chem_key),-c(1,2)]
    L <- climate[match(key,climate_key),-c(1:3)]
    audit <- data.frame(predictor=c(names(C),names(L)),available=colSums(is.finite(as.matrix(cbind(C,L)))),
                        total=length(key),chemistry_keys=sum(key%in%chem_key),climate_keys=sum(key%in%climate_key))
    # Predictor selection is based on eligibility and ecological interpretation before checking skill.
    E <- cbind(pH=C[["pH"]], conductivity=log1p(C[["Conductivity_mS/m"]]),
               air_temperature=L[[1]], precipitation=log1p(L[[2]]))
    p <- st5_input(data.frame(site=s$id[keep],time=s$year[keep]),
                   abundance=a$harmonised[keep,,drop=FALSE],env=E)
    list(p=p,audit=audit)
  })
  list(inputs=lake_inputs,blocks=spatial_lake)
}

st5_moth_archive <- function(moth_root) {
  moth <- data.table::fread(file.path(moth_root,"data/raw/biotime-604.csv"))
  moth[,site:=substring(SAMPLE_DESC,6)]
  obs <- unique(moth[,.(site,time=YEAR)])
  keys <- paste(obs$site,obs$time)
  taxon <- sort(unique(moth$valid_name))
  x <- matrix(FALSE,nrow(obs),length(taxon),dimnames=list(NULL,taxon))
  x[cbind(match(paste(moth$site,moth$YEAR),keys),match(moth$valid_name,taxon))] <- moth$ABUNDANCE>0
  # Duplicate records use any positive detection.
  i <- which(moth$ABUNDANCE>0)
  x[cbind(match(paste(moth$site[i],moth$YEAR[i]),keys),match(moth$valid_name[i],taxon))] <- TRUE
  stopifnot(nrow(obs)==282L)
  st5_input(obs,incidence=x)
}
