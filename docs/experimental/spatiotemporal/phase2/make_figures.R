root <- "docs/experimental/spatiotemporal"
source(file.path(root, "phase2_load.R"))
st2_load(root)
out <- Sys.getenv("ST2_OUTPUT", file.path(root, "phase2/results"))
target <- Sys.getenv("ST2_FIGURES", file.path(root, "phase2/figures"))
dir.create(target, recursive = TRUE, showWarnings = FALSE)
read <- function(name) utils::read.csv(file.path(out, paste0(name, ".csv")))
ink <- .st_colours
types <- c(presence = "Incidence", abundance = "Minimum relative abundance", dominance = "Dominance")
titles <- c("Joint zeta under two temporal designs", "Occurrence and conditional abundance",
            "Transitions by initial sampled occupancy", "Matched environmental contrasts",
            "Factorial response recovery", "Prediction and coherence", "Information across orders",
            "Calibrating uncertainty")
draw <- function(number, width, height, fn) {
  svg_path <- file.path(target, sprintf("phase2-%02d.svg", number))
  grDevices::svg(svg_path, width, height, family = "sans", pointsize = 12)
  fn()
  grDevices::dev.off()
  fn()
}
style <- function(layout, mar = c(4, 4.3, 3.7, .9), outer = c(.2, 0, 2.3, 0)) {
  par(mfrow = layout, mar = mar, oma = outer, las = 1, bty = "l",
      mgp = c(2.65, .65, 0), col.lab = "#233440", col.axis = "#364858", col.main = "#142936")
}
heading <- function(s) mtext(s, 3, outer = TRUE, font = 2, line = .6, cex = 1.13)
heatmap <- function(..., main) {
  plot_st_surface(..., main = NULL)
  mtext(main, side = 3, line = 1.43, font = 2, cex = 0.96)
}
plot_line <- function(x, values, names, title, ylab, xlab = "Climate change (SD)", ylim = c(0, 1), integer_ticks = FALSE) {
  plot(range(x), ylim, type = "n", xlab = xlab, ylab = ylab, main = title,
       xaxt = if (integer_ticks) "n" else "s")
  if (integer_ticks) axis(1, at = x)
  palette <- unname(ink[c("blue", "orange", "green", "grey")])
  for (i in seq_len(ncol(values))) lines(x, values[, i], col = palette[i], lty = i, lwd = 2.1)
  legend("topright", names, bty = "n", col = palette, lty = seq_along(names), lwd = 2, cex = .77)
}
grDevices::cairo_pdf(file.path(target, "phase2-diagnostics.pdf"), width = 12, height = 8,
                      family = "sans", onefile = TRUE, pointsize = 12)
pdf <- grDevices::dev.cur()

draw(1, 12, 7.8, function() {
  style(c(2,3))
  z <- read("zeta_surfaces")
  for (mode in c("all", "consecutive")) for (type in names(types)) {
    a <- z[z$response == type & z$window == mode, ]
    lim <- range(z$zeta[z$response == type])
    heatmap(a, digits = if (type == "abundance") 3 else 1, limits = lim,
                    main = types[[type]], legend = if (mode == "all") "All temporal combinations" else "Consecutive survey windows")
  }
  heading(titles[1])
})

draw(2, 12, 4.4, function() {
  style(c(1,3))
  z <- read("decomposition_surface")
  heatmap(z, "occurrence_mean", "Mean across species", digits = 3, main = "A   Probability of sharing")
  heatmap(z, "conditional_minimum", "Sharing-weighted mean", digits = 3, main = "B   Conditional relative minimum")
  heatmap(z, "zeta_abundance", "Expected sum of relative minima", digits = 3, main = "C   Abundance zeta")
  heading(titles[2])
})

draw(3, 11, 5, function() {
  style(c(1,2))
  a <- read("environmental_response")
  a <- a[a$spatial_order == 3 & a$delta == 0, ]
  plot_line(a$r0, a[c("gain_probability","loss_probability")], c("Gain given absence", "Loss given presence"),
            "A   Site risks", "Probability", "Initially occupied sites of three", c(0, .55), integer_ticks = TRUE)
  plot_line(a$r0, a[c("expansion","contraction","pure_redistribution")], c("Expansion", "Contraction", "Pure redistribution"),
            "B   Combination outcomes", "Probability", "Initially occupied sites of three", c(0, .67), integer_ticks = TRUE)
  heading(titles[3])
})

draw(4, 12, 7.8, function() {
  style(c(2,2))
  a <- read("environmental_response")
  for (r0 in 1:2) {
    b <- a[a$spatial_order == 3 & a$r0 == r0, ]
    plot_line(b$delta, b[c("gain_probability","loss_probability")], c("Gain", "Loss"),
              sprintf("%s   Site risks, %d of 3 occupied", LETTERS[2*r0-1], r0),
              "Conditional event probability", ylim = c(0, .75))
    plot_line(b$delta, b[c("expansion","contraction","pure_redistribution","redistribution")],
              c("Expansion", "Contraction", "Pure redistribution", "Expected balanced replacements"),
              sprintf("%s   Combination outcomes, %d of 3 occupied", LETTERS[2*r0], r0),
              "Probability or expected count", ylim = c(0, .85))
  }
  heading(titles[4])
})

draw(5, 14, 7.8, function() {
  style(c(2,4), mar = c(3.8, 3.8, 3.3, .8))
  d <- read("factorial_contrasts")
  a <- subset(d, target > 0 & group == target & architecture == "M2" & estimand == "trajectory_surface" &
                ((scenario %in% c("gain","loss") & response == "presence") |
                   (scenario == "abundance" & response == "abundance")))
  a$scale <- ave(abs(a$truth), a$scenario, FUN = function(t) sqrt(mean(t^2)))
  a$square <- (a$error / a$scale)^2
  z <- aggregate(square ~ target + spatial_order + temporal_order, a, mean)
  z$standardized_rmse <- sqrt(z$square)
  for (i in 1:8) {
    flags <- expand.grid(c("Low","High"), c("Low","High"), c("Low","High"), stringsAsFactors = FALSE)[i, ]
    title <- paste(i, " ", substr(flags[1,1],1,1), substr(flags[1,2],1,1), substr(flags[1,3],1,1))
    heatmap(z[z$target == i, ], "standardized_rmse", "Scaled target-group RMSE",
                    limits = c(0, max(z$standardized_rmse)), main = title)
  }
  heading(paste(titles[5], "| L/H indicate occupancy, autocorrelation and abundance"))
})

draw(6, 12, 8, function() {
  style(c(2,2), mar = c(4.4, 4.7, 3, .8))
  d <- read("prediction_scores")
  a <- subset(d, spatial_order == 2 & temporal_order == 1 & model %in% c("M0","M1","M2") &
                  response %in% c("shared_gain","shared_loss"))
  means <- aggregate(log_loss ~ split + model, a, mean)
  model <- c("M0","M1","M2")
  splits <- c("sites","regions","forward","joint")
  mat <- xtabs(log_loss ~ model + split, means)[model, splits]
  barplot(mat, beside = TRUE, col = ink[c("grey","blue","green")], ylim = c(0,max(mat)*1.25),
          main = "A   Shared-event prediction (k=2, m=1)", ylab = "Mean held-out log loss", xlab = "Holdout")
  legend("topleft", model, fill = ink[c("grey","blue","green")], bty = "n", cex = .79)
  a <- subset(d, spatial_order == 2 & temporal_order == 1 & model %in% c("M1","M1_interaction") &
                  response %in% c("shared_gain","shared_loss") & scenario %in% c("null","interaction"))
  b <- reshape(a[c("replicate","scenario","split","response","model","log_loss")],
               idvar = c("replicate","scenario","split","response"), timevar = "model", direction = "wide")
  b$improvement <- b$log_loss.M1 - b$log_loss.M1_interaction
  boxplot(improvement ~ scenario, b, names = c("Planted interaction","Clean null"),
          ylab = "Decrease in shared-event log loss", main = "B   Added position × change",
          col = c("#B7D7E8","#D8DDE2"))
  abline(h = 0, lty = 2, col = ink["grey"])
  c <- read("prediction_coherence")
  z <- aggregate(cbind(tested,violations) ~ model + comparison, subset(c,tested>0), sum)
  z$rate <- z$violations / z$tested
  m <- xtabs(rate ~ model + comparison, z)[model, c("overlap","space","time")]
  m["M0", "overlap"] <- NA_real_
  mids <- barplot(m, beside = TRUE, col = ink[c("grey","blue","green")], ylim = c(0,1.02),
          names.arg = c("Same site","Nested space","Nested time"), ylab = "Fraction incompatible",
          main = "C   Overlap and nesting checks")
  text(mids[1L,1L], .035, "NA", cex = .68)
  legend("topleft", model, fill = ink[c("grey","blue","green")], bty = "n", cex = .79)
  c <- read("prediction_calibration")
  b <- aggregate(cbind(predicted, observed, n) ~ model + bin,
                 c[c$spatial_order == 2 & c$temporal_order == 1 & c$event == "gain" & c$model %in% model, ], sum)
  plot(NA, xlim=c(0,1), ylim=c(0,1), xlab="Mean predicted gain probability",
       ylab="Observed gain fraction", main="D   Shared-gain calibration (k=2, m=1)")
  abline(0,1,col="#ABB5BD",lty=2)
  for (j in seq_along(model)) {
    q <- subset(b, model == c("M0","M1","M2")[j])
    points(q$predicted/q$n,q$observed/q$n,type="b",pch=c(1,16,17)[j],col=ink[c("grey","blue","green")[j]],lty=j)
  }
  legend("topleft",model,pch=c(1,16,17),col=ink[c("grey","blue","green")],bty="n",cex=.79)
  heading(titles[6])
})

draw(7, 12, 7.9, function() {
  style(c(2,3))
  d <- subset(read("order_diagnostic"), event == "gain")
  for (type in names(types)) {
    a <- d[d$response == type, ]
    heatmap(a, "zero_fraction", "Zero fraction", limits = c(0,1), main = types[[type]])
  }
  for (type in names(types)) {
    a <- d[d$response == type, ]
    initial <- a$profiled_position[a$spatial_order==1 & a$temporal_order==1]
    a$position_log_ratio <- log10(a$profiled_position / initial)
    heatmap(a, "position_log_ratio", "log10 position information ratio",
                    main = types[[type]], digits = 2)
  }
  heading(titles[7])
})

draw(8, 12, 7.8, function() {
  style(c(2,2), mar = c(3.9, 8.7, 3, 1))
  a <- read("null_summary")
  methods <- c("naive","delete_site","delete_spatial_block","moving_temporal","crossed")
  labels <- c("Naive","Delete site","Delete region","Moving time block","Crossed")
  for (scenario in c("null","shock")) for (event in c("gain","loss")) {
    # Use explicit indexing so plot-level scenario/event values are unambiguous.
    z <- a[a$scenario == scenario & a$event == event & a$spatial_order == 2L, ]
    z <- z[match(methods,z$method),]
    y <- rev(seq_along(methods))
    plot(z$coverage, y, type="n", xlim=c(0,1), ylim=c(.6,5.4), yaxt="n",
         xlab="Coverage of nominal 95% intervals", ylab="",
         main=paste(if(scenario=="null") "Clean null" else "Shared shocks", event, sep=" / "))
    axis(2,y,labels,tick=FALSE,las=1,cex.axis=.87)
    abline(v=.95,col="#758893",lty=2)
    segments(pmax(0,z$coverage-1.96*z$mc_se),y,pmin(1,z$coverage+1.96*z$mc_se),y,
             col=ink["grey"],lwd=2)
    points(z$coverage,y,col=ink[if(scenario=="null") "blue" else "orange"],pch=c(1,16,16,16,16))
  }
  heading(paste(titles[8],"| 200 panels per regime, spatial order two"))
})
grDevices::dev.off(pdf)
