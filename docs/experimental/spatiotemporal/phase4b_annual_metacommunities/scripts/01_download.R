source("docs/experimental/spatiotemporal/phase4b_annual_metacommunities/scripts/common.R")

download <- function(url, name, bytes = NA_real_, sha = NA_character_, header = NULL) {
  path <- file.path(raw, name)
  valid <- function(p) {
    if (!file.exists(p) || file.size(p) == 0) return(FALSE)
    if (is.finite(bytes) && file.size(p) != bytes) return(FALSE)
    if (!is.na(sha) && digest::digest(p, algo = "sha256", file = TRUE) != sha) return(FALSE)
    if (!is.null(header) && !startsWith(readLines(p, n = 1L, warn = FALSE), header)) return(FALSE)
    TRUE
  }
  if (!valid(path)) {
    partial <- tempfile(tmpdir = raw)
    on.exit(unlink(partial), add = TRUE)
    status <- system2("curl", c("-fLsS", "--max-time", "120", shQuote(url), "-o", shQuote(partial)))
    if (status == 0L && valid(partial)) {
      stopifnot(file.rename(partial, path))
    } else {
      return(data.frame(name, url, bytes, sha256 = sha, status = paste("unavailable or unverified (curl", status, ")")))
    }
  }
  data.frame(name, url, bytes = file.size(path),
             sha256 = digest::digest(path, algo = "sha256", file = TRUE), status = "verified")
}

base <- "https://biotime.st-andrews.ac.uk/dl_request.php"
manifest <- list(
  download(paste0(base,"?dl=rawdata&study=604"), "biotime-604.csv", header = "ABUNDANCE,BIOMAS,valid_name,"),
  download(paste0(base,"?dl=metadata&study=604"), "biotime-604-metadata.csv", header = "STUDY_ID,REALM,CLIMATE,"),
  download("https://b2share.eudat.eu/records/ccynf-ak936/files/S048-S060%20S072.xlsx?download=1",
           "biotime-604-original.xlsx", bytes = 2056400,
           sha = "44361064cba1696d9f78c4ea80770938c4f5072f74304280d44c976775fd7c50"))
dryad <- "https://datadryad.org"
groups <- list(
  "220216" = c("observations_coarse.csv", "README_model_validation.txt"),
  "118759" = c("Data_Trends_Moths.txt", "readme_Hällfors_et_al._2021.docx",
               "Data_FIBMS_count.csv", "Data_FIBMS_transect.csv", "Data_FIBMS_visit.csv"))
for (version in names(groups)) {
  index <- download(paste0(dryad, "/api/v2/versions/", version, "/files"), paste0("dryad-", version, ".json"), header = "{")
  manifest <- c(manifest, list(index))
  if (index$status != "verified") next
  records <- jsonlite::fromJSON(file.path(raw, index$name), simplifyVector = FALSE)$`_embedded`$`stash:files`
  for (item in records) {
    if (!item$path %in% groups[[version]]) next
    id <- basename(item$`_links`$self$href)
    url <- paste0(dryad, "/downloads/file_stream/", id)
    record <- download(url, item$path, item$size, item$digest)
    if (record$status != "verified") {
      api <- paste0(dryad, item$`_links`$`stash:download`$href)
      record <- download(api, item$path, item$size, item$digest)
      record$url <- paste(url, api)
    }
    manifest <- c(manifest, list(record))
  }
}
write_result(data.table::rbindlist(manifest), "download_manifest")
stopifnot(manifest[[1]]$status == "verified", manifest[[2]]$status == "verified")
