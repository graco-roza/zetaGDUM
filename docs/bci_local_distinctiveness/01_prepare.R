# 01_prepare.R: build quadrat x species matrices for the 1985 (census 2) and 2015 (census 8) BCI censuses at the resolution of Legendre & Condit (2019):
# 1250 quadrats of 20 x 20 m; individual trees with status "A" (alive) and dbh >= 10 mm (verified against their published counts, see out/verification_counts.csv).
# Inputs (not redistributed): raw/dryad/{tree,stem}/bci.{tree,stem}N.rdata, raw/dryad/bci.spptable.rdata (Condit et al. 2019, doi:10.15146/5xcp-0d46),
# raw/habitat.txt (forestgeo/bciex data-raw/from_suzanne_lao/bcifiles/habitat.txt; Harms et al. 2001 habitat map at 20 m).
suppressMessages(library(dplyr)); dir.create("out", FALSE)
h <- read.delim("raw/habitat.txt"); h$quad <- sprintf("x%03d_y%03d", h$x / 20 + 1, h$y / 20 + 1); stopifnot(nrow(h) == 1250, !anyDuplicated(h$quad))
ld <- function(tp, cs) { e <- new.env(); load(sprintf("raw/dryad/%s/bci.%s%d.rdata", tp, tp, cs), e); get(ls(e)[1], e) }
mat <- function(tp, cs, rule = "tree_A_dbh10") { d <- ld(tp, cs); ok <- d$status == "A" & !is.na(d$status) & !is.na(d$gx) & !is.na(d$gy); if (rule != "stem_all") ok <- ok & !is.na(d$dbh) & d$dbh >= 10; d <- d[ok, ]
  qx <- pmin(floor(d$gx / 20), 49); qy <- pmin(floor(d$gy / 20), 24); d$quad <- sprintf("x%03d_y%03d", qx + 1, qy + 1)
  x <- as.matrix(unclass(table(factor(d$quad, h$quad), d$sp))); storage.mode(x) <- "double"; x }
out <- list(); ver <- list()
for (cs in c(2, 8)) { yr <- c(`2` = 1985, `8` = 2015)[[as.character(cs)]]
  for (tp in c("tree", "stem")) { x <- mat(tp, cs); out[[paste(yr, tp, sep = "_")]] <- x
    hb <- h$habitat[match(rownames(x), h$quad)]; ver[[length(ver) + 1]] <- data.frame(year = yr, table = tp, individuals = sum(x), swamp_individuals = sum(x[hb == "swamp", ]), nonswamp_individuals = sum(x[hb != "swamp", ]), species = sum(colSums(x) > 0), empty_quadrats = sum(rowSums(x) == 0)) } }
V <- bind_rows(ver); u <- length(union(colnames(out$`1985_tree`)[colSums(out$`1985_tree`) > 0], colnames(out$`2015_tree`)[colSums(out$`2015_tree`) > 0]))
write.csv(V, "out/verification_counts.csv", row.names = FALSE); print(V); cat("species union 1985/2015 (tree table):", u, "(L&C: 328)\n")
cat("swamp stems 1985 -> 2015 (tree):", sum(out$`1985_tree`[h$habitat == "swamp", ]), "->", sum(out$`2015_tree`[h$habitat == "swamp", ]), " gain", sum(out$`2015_tree`[h$habitat == "swamp", ]) - sum(out$`1985_tree`[h$habitat == "swamp", ]), "(L&C: 1354)\n")
cat("non-swamp loss:", sum(out$`1985_tree`[h$habitat != "swamp", ]) - sum(out$`2015_tree`[h$habitat != "swamp", ]), "(L&C: 35,697)\n")
sp <- local({ e <- new.env(); load("raw/dryad/bci.spptable.rdata", e); get(ls(e)[1], e) }); saveRDS(list(mat = out, habitat = h, spp = sp[, intersect(c("sp", "Latin", "Genus", "Species", "Family"), names(sp))]), "out/bci_matrices.rds")
