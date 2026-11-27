# fig_common.R: shared set-up for the BCI figures
suppressMessages({library(ggplot2); library(patchwork); library(dplyr); library(tidyr)})
pal <- c(petrol = "#1F5F6B", orange = "#C8641E", sage = "#8AA37B", grey = "#7A7F85")
hab_col <- c(swamp = "#C8641E", stream = "#1F5F6B", slope = "#8AA37B", low_plateau = "#7A7F85", hi_plateau = "#B9C4C8", mixed = "#D8C7A0", young = "#6F4E7C")
hab_lab <- c(swamp = "Swamp", stream = "Stream", slope = "Slope", low_plateau = "Low plateau", hi_plateau = "High plateau", mixed = "Mixed", young = "Young forest")
th <- function(b = 7) theme_classic(base_size = b) + theme(panel.grid = element_blank(), strip.background = element_blank(), strip.text = element_text(face = "bold", size = b), legend.key.size = unit(3, "mm"), legend.margin = margin(0, 0, 0, 0), plot.tag = element_text(face = "bold", size = 9), axis.line = element_line(linewidth = .3), axis.ticks = element_line(linewidth = .3))
thm <- function(b = 7) theme_void(base_size = b) + theme(plot.tag = element_text(face = "bold", size = 9), legend.key.size = unit(3, "mm"), legend.margin = margin(0, 0, 0, 0), legend.position = "bottom", legend.key.width = unit(10, "mm"), legend.title = element_text(size = b), plot.title = element_text(size = b, face = "bold"))
S <- readRDS("out/site_metrics.rds"); S$habitat <- factor(S$habitat, names(hab_col)); K <- c(2, 3, 4, 5, 8, 12, 20, 40, 80, 160, 320, 640)
sv <- function(p, name, w, h) for (e in c("pdf", "svg", "png")) ggsave(sprintf("figures/%s.%s", name, e), p, width = w, height = h, dpi = 220, device = if (e == "pdf") cairo_pdf else NULL)
# significance (Holm) from the large permutation test when available
SIG <- if (file.exists("out/lcbd_tests.csv")) read.csv("out/lcbd_tests.csv") else NULL
files_o <- list.files("out", "^lcbd_tests_occ_.*csv$", full.names = TRUE); SIGO <- if (length(files_o)) bind_rows(lapply(files_o, read.csv)) else NULL
add_sig <- function(d, occ = FALSE) { T0 <- if (occ) SIGO else SIG; if (is.null(T0)) { d$holm <- NA; return(d) }; m <- T0[match(paste(d$set, d$quad), paste(T0$set, T0$quad)), ]; d$holm <- m$p_holm <= .05; d }
map_layer <- function(d, v, ttl, sw = TRUE, sig = FALSE, option = "viridis", dir = 1, lab = NULL, trans = "identity", occ = FALSE) {
  d <- add_sig(d, occ); p <- ggplot(d, aes(x_m + 10, y_m + 10)) + geom_tile(aes(fill = .data[[v]]), width = 20, height = 20)
  if (sw) p <- p + geom_tile(data = d[d$habitat == "swamp", ], fill = NA, colour = pal[["orange"]], width = 20, height = 20, linewidth = .25)
  if (sig && !all(is.na(d$holm))) p <- p + geom_point(data = d[which(d$holm), ], shape = 21, size = .8, stroke = .25, fill = "white", colour = "black")
  p + coord_equal(expand = FALSE) + scale_fill_viridis_c(option = option, direction = dir, name = lab %||% v, trans = trans, guide = guide_colourbar(title.position = "top", title.hjust = .5, barwidth = unit(28, "mm"), barheight = unit(2.5, "mm"))) + labs(title = ttl) + thm() }
`%||%` <- function(a, b) if (is.null(a)) b else a
