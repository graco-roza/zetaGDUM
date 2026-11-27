# 03_lcbd_tests.R: permutation tests of Hellinger LCBD with enough permutations for Holm correction over 1250 quadrats (min p must be < 0.05/1250):
# reproduces the significance pattern of Legendre & Condit (2019, Fig. 2). Tree records, 1985 and 2015.
.libPaths(c(Sys.getenv("ZG_LIB6"), .libPaths())); suppressMessages({library(adespatial); library(dplyr)}); B <- readRDS("out/bci_matrices.rds"); h <- B$habitat; NP <- as.integer(Sys.getenv("NPERM", "99999")); out <- list()
for (nm in c("1985_tree", "2015_tree")) { x <- B$mat[[nm]]; x <- x[, colSums(x) > 0]; set.seed(2019); t0 <- Sys.time(); bd <- beta.div(x, method = "hellinger", nperm = NP, adj = TRUE, save.D = FALSE); cat(nm, "elapsed", format(Sys.time() - t0), "\n")
  d <- data.frame(set = nm, quad = rownames(x), habitat = h$habitat[match(rownames(x), h$quad)], LCBD = bd$LCBD, p = bd$p.LCBD, p_holm = bd$p.adj); out[[nm]] <- d
  cat("significant p<=0.05:", sum(d$p <= .05), " Holm:", sum(d$p_holm <= .05), " swamp sig p:", sum(d$p[d$habitat == "swamp"] <= .05), " swamp Holm:", sum(d$p_holm[d$habitat == "swamp"] <= .05), "\n"); print(table(d$habitat[d$p_holm <= .05])) }
write.csv(bind_rows(out), "out/lcbd_tests.csv", row.names = FALSE)
