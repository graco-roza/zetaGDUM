# lib6.R: shared set-up. Package: zetaGDUM at commit f7e115f installed in ZG_LIB6.
.libPaths(c(Sys.getenv("ZG_LIB6"), .libPaths())); suppressMessages({library(zetaGDUM); library(dplyr); library(tidyr)})
D <- readRDS("out/prepared.rds")
# analysis environment: log10 TP, TN, Chla (log scale; spreads and positions then represent proportional differences)
get_env <- function(d, chla = TRUE) { e <- data.frame(TP = log10(d$env$TP), TN = log10(d$env$TN)); if (chla) e$Chla <- log10(d$env$Chla); e }
get_xy <- function(d) as.matrix(d$env[, c("lon", "lat")])
