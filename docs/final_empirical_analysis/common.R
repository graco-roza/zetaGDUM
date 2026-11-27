# common.R: shared set-up for the final empirical analyses (Figures 3-6). Communities: BCI 1985 (tree records, 1250 quadrats), EPA NLA 2017 macroinvertebrates and zooplankton.
# Package: zetaGDUM at commit f7e115f, unmodified (library path in ZG_LIB6). All stochastic steps use fixed seeds declared in the scripts.
.libPaths(c(Sys.getenv("ZG_LIB6"), .libPaths())); suppressMessages({library(zetaGDUM); library(dplyr); library(tidyr); library(ggplot2); library(patchwork)})
pal <- c(petrol = "#1F5F6B", orange = "#C8641E", sage = "#8AA37B", grey = "#7A7F85")
COMM <- c(BCI = "BCI trees (1985)", BM = "EPA macroinvertebrates", ZP = "EPA zooplankton"); COMM_COL <- c(`BCI trees (1985)` = pal[["sage"]], `EPA macroinvertebrates` = pal[["petrol"]], `EPA zooplankton` = pal[["orange"]])
th <- function(b = 7) theme_classic(base_size = b) + theme(panel.grid = element_blank(), strip.background = element_blank(), strip.text = element_text(face = "bold", size = b), legend.key.size = unit(3, "mm"), legend.margin = margin(0, 0, 0, 0), plot.tag = element_text(face = "bold", size = 9), axis.line = element_line(linewidth = .3), axis.ticks = element_line(linewidth = .3))
sv <- function(p, name, w, h, dir = "figures") { dir.create(dir, FALSE, TRUE); for (e in c("pdf", "svg", "png")) ggsave(file.path(dir, sprintf("%s.%s", name, e)), p, width = w, height = h, dpi = 220, device = if (e == "pdf") cairo_pdf else NULL) }
POOLS <- readRDS("out/pools.rds")       # list(shared = list(BM, ZP), full = list(BM, ZP), bci = list(x, habitat))
get_comm <- function(code, pool = "shared") if (code == "BCI") POOLS$bci$x else POOLS[[pool]][[code]]$x
