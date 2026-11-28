source("docs/experimental/spatiotemporal/phase4_swedish_lakes/scripts/common.R")
taxa_files <- c("03_PhytoBioVtaxa.tsv","05_LittoralTaxa.tsv","07_ProfundalTaxa.tsv")
metric_files <- c("04_PhytoMetrics.tsv","06_LittoralMetrics.tsv","08_ProfundalMetrics.tsv")
lake_name <- function(x) tolower(gsub("[[:space:].-]","",stringi::stri_trans_general(x,"Latin-ASCII")))
taxon_label <- function(x) {
  z <- trimws(stringi::stri_trans_nfc(x))
  z <- gsub("[,.]+_","_",z)
  z <- gsub("(_gr\\.?|-group)$","_group",z,ignore.case=TRUE)
  z
}
genus_label <- function(x) {
  token <- sub("_.*","",x)
  # Preserve higher ranks, non-binomial labels, and complexes spanning genera.
  cross <- grepl("/[A-Z][a-z]+",x) | grepl("^[A-Za-z]+-[A-Z]",x) |
    grepl("_[A-Z][a-z]+_",x)
  higher <- grepl("(idae|inae|aceae|ales|phyceae|ini)$",token) |
    token %in% c("Oligochaeta","Nematoda","Hirudinea","Bivalvia","Gastropoda","Diptera","Trichoptera",
                 "Coleoptera","Monader","Turbellaria","Nemertea","Hydracarina")
  ifelse(cross | higher,x,token)
}
doc <- xml2::read_xml(unz(raw("GCB-31-e70274-s001.docx"),"word/document.xml"))
ns <- xml2::xml_ns(doc)
tab <- xml2::xml_find_first(doc,".//w:tbl",ns)
geo <- do.call(rbind,lapply(xml2::xml_find_all(tab,"./w:tr",ns)[-(1:2)],function(row) {
  cells <- xml2::xml_find_all(row,"./w:tc",ns)
  vapply(cells[1:3],function(c) paste(xml2::xml_text(xml2::xml_find_all(c,".//w:t",ns)),collapse=""),"")
}))
geo <- data.frame(publication_name=trimws(geo[,1]),latitude=as.numeric(geo[,2]),longitude=as.numeric(geo[,3]))
catalogue_lakes <- unique(fread(raw(taxa_files[1]))[,c("EU_ID","Lake_name")])
# Lake names are unique after normalizing punctuation and accents. Explicit
# remaining aliases, if any, are recorded in the join table before analysis.
g <- match(lake_name(catalogue_lakes$Lake_name),lake_name(geo$publication_name))
aliased <- is.na(g)
aliases <- c("V. Rännöbodsjön"="Rännöbodsjön")
g[aliased] <- match(unname(aliases[catalogue_lakes$Lake_name[aliased]]),geo$publication_name)
# Use the RT90 location code only to choose among homonymous Table S1 rows.
normal <- lake_name(catalogue_lakes$Lake_name)
names_geo <- lake_name(geo$publication_name)
homonyms <- which(duplicated(names_geo)[match(normal,names_geo)] |
                  duplicated(names_geo,fromLast=TRUE)[match(normal,names_geo)])
homonyms <- homonyms[is.finite(homonyms)]
codes <- do.call(rbind,strsplit(sub("^SE","",catalogue_lakes$EU_ID),"-",fixed=TRUE))
prefix <- data.frame(x=as.numeric(codes[,2])*10,y=as.numeric(codes[,1])*10)
xy <- sf::st_coordinates(sf::st_transform(sf::st_as_sf(prefix,coords=c("x","y"),crs=3021),4326))
for(i in homonyms) {
  candidates <- which(names_geo==normal[i])
  distance <- (geo$longitude[candidates]-xy[i,1])^2+(geo$latitude[candidates]-xy[i,2])^2
  if (min(distance)>.03^2) stop("Could not disambiguate a homonymous site.")
  g[i] <- candidates[which.min(distance)]
}
stopifnot(!anyNA(g), !anyDuplicated(g))
geography <- cbind(catalogue_lakes,geo[g,,drop=FALSE],row.names=NULL)
geography$match_rule <- ifelse(aliased,"documented Table S1 name alias","normalized name from Table S1")
geography$match_rule[homonyms] <- "homonym resolved using EU_ID RT90 location code"
geography$ecoregion <- NA_character_
# Coordinates supply contiguous latitude blocks independent of outcome.
geography$spatial_block <- NA_integer_
take <- which(is.finite(geography$latitude))
geography$spatial_block[take[order(geography$latitude[take])]] <- ceiling(seq_along(take)*6/length(take))
write.csv(geography,processed("lake_geography.csv"),row.names=FALSE,na="")
datasets <- list()
audits <- mapping <- list()
for(i in seq_along(assemblages)) {
  nm <- assemblages[i]
  t <- fread(raw(taxa_files[i]))
  m <- fread(raw(metric_files[i]))
  key <- paste(t$EU_ID,t$Year)
  metric_key <- paste(m$EU_ID,m$Year)
  if (anyDuplicated(key) || anyDuplicated(metric_key) || anyDuplicated(names(t)) ||
      anyNA(t$EU_ID) || anyNA(t$Year)) stop("Nonunique or invalid lake-year/taxon keys in ",nm)
  response <- as.matrix(t[,-(1:3)])
  if (!is.numeric(response) || any(!is.finite(response) & !is.na(response)) ||
      any(response<0,na.rm=TRUE)) stop("Invalid Hellinger values in ",nm)
  square <- response^2
  mass <- rowSums(square)
  j <- match(key,metric_key)
  richness <- m[[4]][j]
  amount <- m[[5]][j]
  status <- ifelse(is.na(mass),"missing_taxon_values",ifelse(mass>0,"observed_nonempty",
            ifelse(is.finite(amount) & amount>0,"zero_matrix_positive_metrics",
            ifelse(is.finite(amount) & amount==0 & richness==0,"verified_empty","zero_matrix_missing_metrics"))))
  keep <- which(status=="observed_nonempty")
  if (max(abs(mass[keep]-1))>1e-5) stop("Non-Hellinger row sums in ",nm)
  q <- square
  q[keep,] <- close(q[keep,,drop=FALSE])
  q[-keep,] <- NA_real_
  observed <- colSums(q[keep,,drop=FALSE]>0)>0
  names0 <- colnames(q)
  harmonised <- taxon_label(names0)
  genus <- genus_label(harmonised)
  sum_groups <- function(identifiers) {
    a <- matrix(NA_real_,nrow(q),length(unique(identifiers)),
                dimnames=list(NULL,unique(identifiers)))
    a[keep,] <- t(rowsum(t(q[keep,,drop=FALSE]),identifiers,reorder=FALSE))
    stopifnot(max(abs(rowSums(a[keep,,drop=FALSE])-1))<1e-12)
    a
  }
  a <- sum_groups(harmonised)
  g <- sum_groups(genus)
  mapping[[i]] <- data.frame(assemblage=nm,archived_taxon=names0,harmonised,genus_sensitivity=genus,
                   mapping_rule=ifelse(names0==harmonised,"retained archive name","equivalent punctuation/group abbreviation"),
                   observed_anywhere=observed)
  samples <- data.frame(assemblage=nm,id=t$EU_ID,lake=t$Lake_name,year=t$Year,status,
                        metric_present=!is.na(j),metric_richness=richness,hellinger_squared_sum=mass)
  datasets[[nm]] <- list(samples=samples,harmonised=a,genus=g)
  audits[[i]] <- data.frame(assemblage=nm,archived_rows=nrow(t),archived_taxon_columns=ncol(q),
                          nonempty=length(keep),observed_raw_taxa=sum(observed),
                          harmonised_taxa=sum(colSums(a[keep,,drop=FALSE]>0)>0),
                          genus_labels=sum(colSums(g[keep,,drop=FALSE]>0)>0),
                          max_mass_error=max(abs(mass[keep]-1)),duplicate_taxon_keys=0,
                          duplicate_metric_keys=0,missing_metrics=sum(is.na(j)),
                          richness_conflicts=sum(rowSums(response[keep,,drop=FALSE]>0)!=richness[keep],na.rm=TRUE))
}
saveRDS(datasets,processed("assemblages.rds"))
write.csv(do.call(rbind,mapping),processed("taxonomy_mapping.csv"),row.names=FALSE)
save_csv(do.call(rbind,audits),"archive_audit")
save_csv(do.call(rbind,lapply(datasets,`[[`,"samples")),"lake_year_status")
