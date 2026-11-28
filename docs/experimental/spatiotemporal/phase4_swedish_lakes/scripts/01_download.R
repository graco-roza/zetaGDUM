root <- "docs/experimental/spatiotemporal/phase4_swedish_lakes"
files <- c("data/01_Physochemical.tsv","data/02_ClimateData.tsv","data/03_PhytoBioVtaxa.tsv",
           "data/04_PhytoMetrics.tsv","data/05_LittoralTaxa.tsv","data/06_LittoralMetrics.tsv",
           "data/07_ProfundalTaxa.tsv","data/08_ProfundalMetrics.tsv","documentation/Metadata.txt")
base <- "https://api.researchdata.se/dataset/2025-136/1/file/"
dir.create(file.path(root,"data/raw"),recursive=TRUE,showWarnings=FALSE)
verify <- function(path, type) {
  if (!file.exists(path) || file.info(path)$size < 250) stop("Invalid download size: ", path)
  preview <- readLines(path,n=5,warn=FALSE,encoding="UTF-8")
  if (!length(preview) || any(grepl("<!DOCTYPE|<html|^\\s*[{]",preview,ignore.case=TRUE,useBytes=TRUE))) {
    stop("Non-data payload: ",path)
  }
  if (type=="tsv") {
    z <- data.table::fread(path,sep="\t",data.table=FALSE,check.names=FALSE)
    if (!nrow(z) || ncol(z)<3 || !any(grepl("year",names(z),ignore.case=TRUE))) stop("Invalid TSV: ",path)
    return(c(rows=nrow(z),columns=ncol(z)))
  }
  if (!any(grepl("data|lake|Johnson", preview, ignore.case=TRUE,useBytes=TRUE))) stop("Unexpected metadata")
  c(rows=length(readLines(path,warn=FALSE)),columns=NA_integer_)
}
fetch <- function(url,dest) {
  system2("curl",c("-L","--fail","--silent","--show-error","--retry","3",
                    "--connect-timeout","15","--max-time","120",shQuote(url),"-o",shQuote(dest)))
}
zip_url <- paste0(base,"zip")
archive <- Sys.getenv("ST4_CATALOGUE_ZIP")
if (!nzchar(archive)) {
  archive <- tempfile(fileext=".zip")
  if (fetch(zip_url,archive)!=0) stop("Catalogue ZIP could not be obtained for checksum verification.")
}
members <- utils::unzip(archive,list=TRUE)$Name
member <- paste0("2025-136-1/",files)
if (!all(c(member,"sha256.txt") %in% members)) stop("Incorrect catalogue ZIP contents.")
sums <- read.table(unz(archive,"sha256.txt"),col.names=c("sha256","member"),colClasses="character")
expected <- sums$sha256[match(member,sums$member)]
if (anyNA(expected) || any(nchar(expected)!=64L)) stop("Invalid catalogue checksums.")
download <- function(i) {
  rel <- files[i]
  dest <- file.path(root,"data/raw",basename(rel))
  type <- tools::file_ext(dest)
  hash <- function(path) digest::digest(path,algo="sha256",file=TRUE)
  source <- "cached, matched catalogue SHA-256"
  url <- paste0(base,rel)
  if (!file.exists(dest) || Sys.getenv("ST4_REFRESH")=="1" || hash(dest)!=expected[i]) {
    tmp <- tempfile(tmpdir=dirname(dest))
    on.exit(unlink(tmp))
    status <- if (Sys.getenv("ST4_USE_CATALOGUE_ZIP")=="1") 1L else fetch(url,tmp)
    if (status==0 && file.exists(tmp) && hash(tmp)==expected[i]) {
      source <- "individual API path"
    } else {
      extract <- tempfile()
      dir.create(extract)
      on.exit(unlink(extract,recursive=TRUE),add=TRUE)
      extracted <- utils::unzip(archive,files=member[i],exdir=extract)
      if (!file.copy(extracted,tmp,overwrite=TRUE)) stop("Could not extract ",rel)
      source <- "catalogue API ZIP fallback"
    }
    if (hash(tmp)!=expected[i]) stop("Checksum mismatch: ",rel)
    verify(tmp,type)
    if (!file.rename(tmp,dest)) stop("Could not record archive file.")
  }
  dims <- verify(dest,type)
  header <- if(type=="tsv") readLines(dest,n=1,warn=FALSE) else "Archive metadata"
  data.frame(path=rel,url,source,bytes=file.info(dest)$size,sha256=hash(dest),checksum_match=hash(dest)==expected[i],
             rows=unname(dims[1]),columns=unname(dims[2]),
             first_fields=paste(head(strsplit(header,"\t",fixed=TRUE)[[1]],5),collapse=" | "))
}
result <- parallel::mclapply(seq_along(files),download,mc.cores=if(.Platform$OS.type=="windows")1L else 4L)
bad <- vapply(result,inherits,logical(1),"try-error")
if (any(bad)) stop(paste(unlist(result[bad]),collapse="\n"))
write.csv(do.call(rbind,result),file.path(root,"data/raw/manifest.csv"),row.names=FALSE,na="")

supp <- file.path(root,"data/raw/GCB-31-e70274-s001.docx")
if (!file.exists(supp)) {
  zip <- tempfile(fileext=".zip")
  url <- "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12146949/supplementaryFiles"
  if (fetch(url,zip)!=0) stop("Article supplement was unavailable.")
  utils::unzip(zip,files=basename(supp),exdir=dirname(supp))
}
if (!"word/document.xml" %in% utils::unzip(supp,list=TRUE)$Name) stop("Invalid article supplement.")
write.csv(data.frame(file=basename(supp),bytes=file.info(supp)$size,
                     sha256=digest::digest(supp,algo="sha256",file=TRUE),
                     url="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12146949/supplementaryFiles"),
          file.path(root,"data/raw/supplement_manifest.csv"),row.names=FALSE)
