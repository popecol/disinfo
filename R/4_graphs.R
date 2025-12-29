
# Code for the paper:
# Zalewska K., Skoracka A., Bonte D., Puchalska E., Lewandowski M., Kuczynski L.
# Is passive dispersal informed? - Experimental evidence for decision-making in phytophagous arthropods

# Charts from the paper using traditional R graphics.


# Data and model ----------------------------------------------------------

library(glmmTMB)
library(emmeans)

load("data/data.RData")
m <- glmmTMB(D/N ~ cue * env * host_spec + (1|line), weights = N, data, family = betabinomial)
car::Anova(m)


# Figure S# ----------------------------------------------------------------

# Main effects of kairomones (a), current environment (b) and host specialisation (c) on dispersal rates. 

all_terms <- attr(terms(formula(m)), "term.labels")
v <- all_terms[!grepl("[:|]", all_terms)]
lab <- data.frame(v = v, lab = c("Kairomone", "Current environment", "Host specialisation"))

# Marginal effects ----
op <- par(mfrow = c(1, 3), mar = c(11, 5, 2, 2))

for (i in seq(length(v))) {
  cf <- emmeans(m, v[i], type = "response", level = 0.95)
  cf50 <- emmeans(m, v[i], type = "response", level = 0.5)
  dat <- print(cf)
  dat50 <- print(cf50)
  
  k <- nrow(dat)
  x <- seq(k)
  xlim <- extendrange(x, f = ifelse(k == 3, 0.3, 0.7))
  y <- dat[, 2]
  lci <- dat[, 5]; uci <- dat[, 6]
  # ylim <- range(lci, uci)
  # ylim <- extendrange(ylim, f = 0.1)
  ylim <- c(0.03, 0.15)
  lci50 <- dat50[, 5]; uci50 <- dat50[, 6]
  plot(x, y, xlim = xlim, ylim = ylim, xaxt = "n", xlab = "", ylab = "Dispersal rate", type = "n", cex.lab = 2)
  axis(1, at = x, labels = dat[, 1], cex.axis = 1.5, las = 2)
  mtext(lab$lab[i], 1, 8, cex = 1.3)
  segments(x, lci50, x, uci50, lwd = 3)
  segments(x, lci, x, uci, lwd = 1)
  points(x, y, cex = 3, col = "white", pch = 20)
  points(x, y, cex = 1.3)
  mtext(paste0("(", letters[i], ")"), side = 2, line = 2.5, padj = -5.4, cex = 1.5, las = 1)
}
par(op)


# Figure 1 ----------------------------------------------------------------

plot_2_way <- function(model, formula, col_vec = c("black", "red", "blue"), 
                       show_sig_top = TRUE, show_sig_bottom = TRUE, 
                       add_legend = TRUE, legend_pos = "topright",
                       legend_buffer = 0.45, 
                       bottom_buffer = NULL,
                       ylab_dist = 3.5, 
                       xlab_dist = 3.0, 
                       title_cex = 1.3, # Size of X and Y titles
                       axis_cex = 1.2,  # Size of tick labels
                       pt_cex = 1.5,    # Size of data points
                       pt_lwd = 2.5,    # Thickness of point borders
                       leg_cex = 1.3,   # Size of legend text
                       verbose = TRUE) {
  
  # Plots 2-way interactions from a model with EMMeans and significance brackets.
  #
  # Arguments:
  #           model: A model object supported by the 'emmeans' package.
  #         formula: A formula specifying the interaction (e.g. '~ x_axis | grouping').
  #         col_vec: Vector of colors for the grouping levels (default black/red/blue).
  #    show_sig_top: Logical. Show significance brackets between groups (top).
  # show_sig_bottom: Logical. Show significance brackets between x-levels (bottom).
  #      add_legend: Logical. Whether to draw the legend.
  #      legend_pos: Position keyword for the legend (default "topright").
  #   legend_buffer: Numeric. Fraction of y-range to reserve at top for legend.
  #   bottom_buffer: Numeric. Fraction of y-range to reserve at bottom for brackets.
  #       ylab_dist: Numeric. Margin line distance for Y-axis title.
  #       xlab_dist: Numeric. Margin line distance for X-axis title.
  #       title_cex: Numeric. Text size multiplier for axis titles.
  #        axis_cex: Numeric. Text size multiplier for axis tick labels.
  #          pt_cex: Numeric. Size of the data points (pch).
  #          pt_lwd: Numeric. Line width of the data point borders.
  #         leg_cex: Numeric. Text size multiplier for the legend.
  #         verbose: Logical. Print summary of significant contrasts to console.
  #
  # Returns:
  #   No return value; produces a plot on the current graphics device.
  #
  # Details:
  #   The function calculates estimated marginal means (95% CI) and performs 
  #   pairwise contrasts. It expands the plotting canvas vertically to accommodate 
  #   significance brackets. Top brackets compare grouping levels within an X-level; 
  #   bottom brackets compare X-levels within a group. The Y-axis limits are strictly 
  #   enforced at 0 for logical consistency, but brackets can be drawn below 0.
  
  # --- 1. Calculate Estimated Marginal Means (95% only) ---
  f <- formula(formula)
  cf <- emmeans(model, f, type = "response", level = 0.95)
  dat <- data.frame(cf)
  
  x_var_name <- names(dat)[1]; by_var_name <- names(dat)[2]
  x_dat <- dat[, 1]; by_dat <- dat[, 2]
  
  levels_x <- levels(x_dat); levels_by <- levels(by_dat)
  n_groups <- length(levels_by)
  
  # --- 2. Calculate Contrasts ---
  sig_data_top <- NULL
  if (show_sig_top) {
    prs_top <- contrast(cf, method = "pairwise", by = x_var_name, adjust = "tukey")
    s_df_top <- as.data.frame(prs_top)
    sig_data_top <- s_df_top[s_df_top$p.value < 0.05, ]
  }
  
  sig_data_bot <- NULL
  if (show_sig_bottom) {
    prs_bot <- contrast(cf, method = "pairwise", by = by_var_name, adjust = "tukey")
    s_df_bot <- as.data.frame(prs_bot)
    sig_data_bot <- s_df_bot[s_df_bot$p.value < 0.05, ]
  }
  
  if (verbose) {
    cat(sprintf("[%s] Significant contrasts: Top=%d, Bottom=%d\n", 
                paste(format(formula), collapse=""),
                if(is.null(sig_data_top)) 0 else nrow(sig_data_top),
                if(is.null(sig_data_bot)) 0 else nrow(sig_data_bot)))
  }
  
  # --- 3. Setup Plot Params ---
  offsets <- seq(-0.15, 0.15, length.out = n_groups)
  cols <- rep(col_vec, length.out = n_groups)
  pch_vec <- c(21, 22, 24, 23, 25) 
  pchs <- rep(pch_vec, length.out = n_groups)
  
  # --- 4. Handle Canvas Expansion ---
  data_max <- max(dat[, 7])
  logical_ylim <- c(0, data_max)
  visual_range <- diff(logical_ylim)
  
  # Top Padding
  top_buffer <- 0.1 
  if(show_sig_top && !is.null(sig_data_top) && nrow(sig_data_top) > 0) {
    top_buffer <- top_buffer + 0.4
  }
  if(add_legend) {
    top_buffer <- top_buffer + legend_buffer 
  }
  pad_top <- top_buffer * visual_range
  
  # Bottom Padding
  if (!is.null(bottom_buffer)) {
    pad_bot <- bottom_buffer * visual_range
  } else {
    # Default: 0.1 base + 0.1 per group
    pad_bot <- if(show_sig_bottom) (0.1 + (0.1 * n_groups)) * visual_range else 0.05 * visual_range
  }
  
  canvas_ylim <- c(logical_ylim[1] - pad_bot, logical_ylim[2] + pad_top)
  
  x_lab <- if(exists("lab") && x_var_name %in% lab$v) lab[lab$v == x_var_name, 2] else x_var_name
  by_lab <- if(exists("lab") && by_var_name %in% lab$v) lab[lab$v == by_var_name, 2] else by_var_name
  
  # --- 5. Draw Plot ---
  plot(1, type = "n", 
       xlim = c(0.5, nlevels(x_dat) + 0.5), 
       ylim = canvas_ylim, 
       xaxt = "n", yaxt = "n", 
       xlab = "", ylab = "", 
       las = 1)
  
  # Titles (Same size)
  mtext("Dispersal rate", side = 2, line = ylab_dist, cex = title_cex)
  mtext(x_lab, side = 1, line = xlab_dist, cex = title_cex)
  
  # Axes
  y_ticks <- pretty(logical_ylim)
  y_ticks <- y_ticks[y_ticks >= 0] 
  axis(2, at = y_ticks, las = 1, cex.axis = axis_cex)
  axis(1, at = 1:nlevels(x_dat), labels = levels_x, cex.axis = axis_cex)
  
  # --- 6. Plot Points ---
  for (i in seq_along(levels_by)) {
    idx <- dat[, 2] == levels_by[i]
    x_pos <- as.numeric(x_dat[idx]) + offsets[i]
    y <- dat[idx, 3]
    lci <- dat[idx, 6]; uci <- dat[idx, 7]
    
    # Draw 95% CI
    segments(x_pos, lci, x_pos, uci, lwd = 2, col = cols[i])
    
    # Draw Points
    points(x_pos, y, pch = pchs[i], bg = "white", col = cols[i], cex = pt_cex, lwd = pt_lwd)
  }
  
  step_h <- visual_range * 0.05 
  
  # --- 7. Draw TOP Brackets ---
  if (show_sig_top && !is.null(sig_data_top) && nrow(sig_data_top) > 0) {
    for (xi in seq_along(levels_x)) {
      curr_level <- levels_x[xi]
      curr_sigs <- sig_data_top[as.character(sig_data_top[[x_var_name]]) == as.character(curr_level), ]
      
      if (nrow(curr_sigs) > 0) {
        local_idx <- as.character(dat[[x_var_name]]) == as.character(curr_level)
        base_y <- max(dat[local_idx, 7]) + step_h
        
        for (r in 1:nrow(curr_sigs)) {
          parts <- unlist(strsplit(as.character(curr_sigs$contrast[r]), "\\s+[-/]\\s+"))
          if (length(parts) >= 2) {
            idx1 <- which(levels_by == trimws(parts[1]))
            idx2 <- which(levels_by == trimws(parts[2]))
            if (length(idx1) > 0 && length(idx2) > 0) {
              x1 <- xi + offsets[idx1]; x2 <- xi + offsets[idx2]
              lines(c(x1, x1, x2, x2), c(base_y, base_y + step_h/3, base_y + step_h/3, base_y), lwd = 1.5)
              p_val <- curr_sigs$p.value[r]
              sym <- ifelse(p_val < 0.001, "***", ifelse(p_val < 0.01, "**", "*"))
              text(mean(c(x1, x2)), base_y + step_h/2, labels = sym, cex = 1.2, pos = 3, offset = 0.1)
              base_y <- base_y + step_h * 1.5
            }
          }
        }
      }
    }
  }
  
  # --- 8. Draw BOTTOM Brackets ---
  if (show_sig_bottom && !is.null(sig_data_bot) && nrow(sig_data_bot) > 0) {
    # Reverted to simple start at 0
    start_y <- 0 
    
    for (gi in seq_along(levels_by)) {
      curr_group <- levels_by[gi]
      curr_sigs <- sig_data_bot[as.character(sig_data_bot[[by_var_name]]) == as.character(curr_group), ]
      
      if (nrow(curr_sigs) > 0) {
        # Standard spacing logic
        lane_y <- start_y - (step_h * 0.5) - ((gi - 1) * (step_h * 1.8))
        
        for (r in 1:nrow(curr_sigs)) {
          parts <- unlist(strsplit(as.character(curr_sigs$contrast[r]), "\\s+[-/]\\s+"))
          if (length(parts) >= 2) {
            idx1 <- which(levels_x == trimws(parts[1]))
            idx2 <- which(levels_x == trimws(parts[2]))
            if (length(idx1) > 0 && length(idx2) > 0) {
              x1 <- idx1 + offsets[gi]; x2 <- idx2 + offsets[gi]
              
              y_top <- lane_y
              y_bottom <- lane_y - (step_h * 0.4)
              lines(c(x1, x1, x2, x2), c(y_top, y_bottom, y_bottom, y_top), lwd = 1.5, col = cols[gi])
              p_val <- curr_sigs$p.value[r]
              sym <- ifelse(p_val < 0.001, "***", ifelse(p_val < 0.01, "**", "*"))
              text(mean(c(x1, x2)), y_bottom, labels = sym, cex = 1.2, col = cols[gi], pos = 1, offset = 0.5) 
            }
          }
        }
      }
    }
  }
  
  # --- 9. Legend ---
  if (add_legend) {
    legend(legend_pos, legend = levels_by, col = cols, pt.bg = "white", pch = pchs, lwd = 2, title = by_lab, bty = "n", cex = leg_cex)
  }
}


op <- par(mfrow = c(1, 3), mar = c(5.5, 6, 3, 2))

# Panel A: Manual bottom buffer
plot_2_way(m, '~ env | cue', legend_buffer = 0.6, bottom_buffer = 0.3)
mtext("(a)", side = 3, line = 1, adj = -0.1, cex = 1.2)

# Panel B: Manual bottom buffer
plot_2_way(m, '~ host_spec | cue', legend_buffer = 0.3, bottom_buffer = 0.25)
mtext("(b)", side = 3, line = 1, adj = -0.1, cex = 1.2)

# Panel C: Manual bottom buffer
plot_2_way(m, '~ env | host_spec', legend_buffer = 0.1, bottom_buffer = 0.2)
mtext("(c)", side = 3, line = 1, adj = -0.1, cex = 1.2)

par(op)


# Figure 4 ----------------------------------------------------------------

library(RColorBrewer)
display.brewer.all(type = "qual")
pal <- brewer.pal(4, "Set1")
col <- adjustcolor(pal, alpha = 0.5)
colt <- adjustcolor(pal, alpha = 0.1)

load("data/data.RData")
data <- subset(data, n_kairo > 0)  # Zero is equivalent to control in the variable ‘cue’.
data <- subset(data, plant == "W") # No mixes for env==‘unknown’.
data <- transform(data, plant = NULL, env = NULL, group = interaction(cue, host_spec, sep = "_"))
data <- droplevels(data)

m <- glmmTMB(D/N ~ 0 + group + group:n_kairo + (1|line), weights = N, data, family = betabinomial)
summary(m)

group <- levels(data$group)
label <- sub("_", " ", group)

n_kairo <- seq(0.9, 3.1, length.out = 100)
nd <- expand.grid(group = group, n_kairo = n_kairo, line = NA, N = 10)

# Excluding unobserved range
idx <- (nd$group == "Unfamiliar_Generalist") & (nd$n_kairo > 2.1)
nd <- nd[!idx, ]

pred <- predict(m, nd, se = TRUE)
fit <- pred$fit; se <- pred$se.fit
q <- -qt(0.025, df.residual(m))
lw <- fit - q * se; up <- fit + q * se
p <- plogis(cbind(fit, lw, up))
nd <- cbind(nd, p)

op <- par(mar = c(5, 5, 3, 2))
plot(fit ~ n_kairo, nd, type = "n", xlim = extendrange(1:3), ylim = range(lw, up), xlab = "No. of kairomones", ylab = "Dispersal rate", cex.lab = 1.5, xaxt = "n")
axis(1, at = 1:3)
for (i in seq(length(group))) {
  g <- group[i]
  dat <- subset(nd, group == g)
  polygon(x = c(dat$n_kairo, rev(dat$n_kairo)), y = c(dat$lw, rev(dat$up)), border = NA, col = colt[i])
  lines(dat$n_kairo, dat$fit, col = col[i], lwd = 3)
}
legend(2, 0.25, legend = label, text.col = "white", bty = "n", col = colt, lwd = 15)
legend(2, 0.25, legend = label, bty = "n", col = pal, lwd = 2)
par(op)


