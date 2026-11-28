source("docs/experimental/spatiotemporal/phase4b_annual_metacommunities/scripts/common.R")
library(data.table)

x <- fread(file.path(raw, "biotime-604.csv"))
meta_lines <- readLines(file.path(raw, "biotime-604-metadata.csv"), n = 2L)
meta_prefix <- sub(",GRAIN_SIZE_TEXT,.*$", "", meta_lines[1L])
meta_row <- sub(",,[0-9.]+,1,1,60,.*$", "", meta_lines[2L])
meta <- read.csv(text = paste(meta_prefix, sub(",$", "", meta_row), sep = "\n"), check.names = FALSE)
stopifnot(nrow(x) == meta$TOTAL,
          all(grepl("^[0-9]{4}_[0-9]{4}$", x$SAMPLE_DESC)),
          all(as.integer(substr(x$SAMPLE_DESC,1,4)) == x$YEAR),
          all(is.finite(x$ABUNDANCE) & x$ABUNDANCE >= 0))
x[, site := substring(SAMPLE_DESC,6)]
coords <- unique(x[, .(site, latitude = LATITUDE, longitude = LONGITUDE)])
stopifnot(!anyDuplicated(coords$site))
ids <- sort(coords$site)
years <- seq.int(min(x$YEAR), max(x$YEAR))
duplicate <- x[, .N, by = .(site,YEAR,valid_name)][N > 1]
write_result(duplicate, "duplicate_site_year_taxa")

original <- file.path(raw,"biotime-604-original.xlsx")
info <- readxl::read_excel(original,"Additional_information", col_types = c(rep("text",4), rep("numeric",3), "text"))
workbook <- as.data.table(readxl::read_excel(original,"Biotic_data", skip = 2,
  col_types = c("text","text","numeric","numeric","skip","skip"), .name_repair = "minimal"))
setnames(workbook, c("site","taxon","year","count"))
workbook <- workbook[!is.na(site) & !is.na(year)]
stopifnot(nrow(workbook) > 0, setequal(ids, info[["trap code"]]), all(is.finite(workbook$count) & workbook$count >= 0))

intersect_coverage <- function(text) {
  segments <- strsplit(gsub(" ", "", text), ",", fixed = TRUE)[[1]]
  lapply(segments, function(segment) {
    limits <- as.integer(strsplit(segment, "-", fixed = TRUE)[[1]])
    if (length(limits) > 2L) return(integer())
    seq.int(limits[1L], tail(limits, 1L))
  })
}
declared <- rbindlist(lapply(seq_len(nrow(info)), function(i) {
  y <- unlist(intersect_coverage(info[["Years with trapping coverage"]][i]))
  data.table(site = info[["trap code"]][i], year = y, workbook_declared = TRUE)
}))
obs <- x[, .(archive_records = .N, positive_records = sum(ABUNDANCE > 0),
             individuals = sum(ABUNDANCE), trap_nights = NA_integer_,
             first_trap_date = NA_character_, last_trap_date = NA_character_), by = .(site, year = YEAR)]
original_years <- workbook[, .(original_records = .N), by = .(site, year)]
effort <- Reduce(function(a,b) merge(a,b,by=c("site","year"),all.x=TRUE),
                 list(CJ(site=ids,year=years), obs, declared, original_years))
effort[, status := fifelse(is.na(archive_records), "no BioTIME observation",
  fifelse(positive_records > 0, "observed community", "archived nonpositive records"))]
effort[, dataset := "BioTIME 604"]
write_result(effort, "sampling_effort")
write_result(merge(coords, as.data.table(info)[, .(site = `trap code`, declared_years = `Years with trapping coverage`)], by="site"),
             "site_identifiers")
write_result(data.table(taxon = sort(unique(x$valid_name)), source = "BioTIME accepted names",
                        operation = "retained verbatim"), "taxonomic_audit")
write_result(data.table(taxon = sort(unique(workbook$taxon)), source = "B2SHARE workbook"), "original_taxon_inventory")

choose_panels <- function(coverage, label) {
  years <- sort(unique(coverage$year))
  intervals <- CJ(start = years, end = years)[end >= start]
  intervals[, duration := end - start + 1L]
  candidate <- lapply(seq_len(nrow(intervals)), function(i) {
    selection <- coverage[year %between% c(intervals$start[i],intervals$end[i]) & eligible,
                          .N, by=site][N == intervals$duration[i], site]
    data.table(sites = length(selection), members = paste(sort(selection), collapse = ";"))
  })
  candidate <- cbind(intervals, rbindlist(candidate))
  candidate[, dataset := label]
  candidate[, tier := fifelse(sites >= 30 & duration >= 15, "target",
                              fifelse(sites >= 20 & duration >= 10, "exploratory", "below minimum"))]
  candidate[]
}

coverage <- effort[, .(site,year,eligible = status == "observed community")]
panels <- choose_panels(coverage,"BioTIME 604")
write_result(panels, "balanced_panel_candidates")
# Archive illustration only. Use coverage and a ten-year minimum to fit both window designs.
illustration <- panels[duration >= 10][order(-sites,-duration,start)][1]
keep <- strsplit(illustration$members,";",fixed = TRUE)[[1]]
yy <- seq.int(illustration$start,illustration$end)
stopifnot(all(coverage[site %in% keep & year %in% yy, eligible]))
write_result(data.table(dataset = "BioTIME 604", site = ids, retained = ids %in% keep,
                        start = min(yy), end = max(yy), use = "archive illustration only"), "balanced_panel_membership")
records <- x[site %in% keep & YEAR %in% yy & ABUNDANCE > 0]
taxa <- sort(unique(records$valid_name))
counts <- array(0, c(length(keep),length(yy),length(taxa)), dimnames = list(keep, yy, taxa))
grouped <- records[, .(count = sum(ABUNDANCE)), by = .(site,YEAR,valid_name)]
counts[cbind(match(grouped$site, keep),match(grouped$YEAR,yy),match(grouped$valid_name,taxa))] <- grouped$count
stopifnot(all(apply(counts > 0, c(1,2), sum) > 0))
saveRDS(list(panel = st_panel(counts > 0, yy), counts = counts,
             stratum = setNames(ceiling(rank(coords[match(keep,site),latitude]) / length(keep) * 4),keep)),
        file.path(processed,"biotime-604.rds"))

manifest <- read_result("download_manifest")
source_status <- function(names) all(names %in% manifest[status == "verified", name])
write_result(data.table(dataset = c("BioTIME 604","Dryad Nocturna","Dryad breeding birds","Dryad FIBMS"),
  files_accessible = c(TRUE, source_status("Data_Trends_Moths.txt"), source_status("observations_coarse.csv"),
                       source_status(c("Data_FIBMS_count.csv","Data_FIBMS_transect.csv","Data_FIBMS_visit.csv"))),
  status = c("below minimum; archived-count illustration",
             "coverage and effort unaudited", "stable routes and protocol unaudited", "visit effort and sites unaudited")),
  "dataset_eligibility")
write_result(data.table(records=nrow(x), advertised_taxa=meta$NUMBER_OF_SPECIES,
  parsed_taxa=uniqueN(x$valid_name), retained_taxa=length(taxa), archive_sites=length(ids),
  positive_site_years=nrow(obs[positive_records>0]), retained_sites=length(keep),
  retained_years=length(yy), duplicate_site_year_taxon_keys=nrow(duplicate),
  dates_available=sum(!is.na(x$DAY) & !is.na(x$MONTH)),
  workbook_only_site_years=effort[!is.na(original_records)&is.na(archive_records),.N],
  archived_duplicates=sum(duplicate$N-1),
  archive_only_names=length(setdiff(unique(x$valid_name),unique(workbook$taxon)))),
  "audit_counts")
