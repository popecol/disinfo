
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


# Figure S3 ----------------------------------------------------------------

# Main effects of kairomones (a), current environment (b) and host specialisation (c) on dispersal rates. 

plot_1_way <- function(model, formula, col_vec = "black", 
                       show_sig = TRUE, 
                       top_buffer = NULL,
                       ylab_dist = 4, 
                       xlab_dist = 8.0, 
                       title_cex = 1.3, 
                       axis_cex = 1.4,  
                       pt_cex = 1.5,    
                       pt_lwd = 2.5,    
                       star_dist = 0.1,
                       bracket_gap = 0.1, 
                       ylim = NULL, 
                       verbose = TRUE) {
  
  # Plots main effects (1-way) from a model with EMMeans and significance brackets.
  #
  # Arguments:
  #           model: A model object supported by the 'emmeans' package.
  #         formula: A formula specifying the main effect (e.g. '~ x_axis').
  #         col_vec: Color (or vector of colors) for the points.
  #        show_sig: Logical. Show pairwise significance brackets on top.
  #      top_buffer: Numeric. Manual fraction of y-range to reserve at top for brackets.
  #       ylab_dist: Numeric. Margin line distance for Y-axis title.
  #       xlab_dist: Numeric. Margin line distance for X-axis title.
  #       title_cex: Numeric. Text size multiplier for axis titles.
  #        axis_cex: Numeric. Text size multiplier for axis tick labels.
  #          pt_cex: Numeric. Size of the data points.
  #          pt_lwd: Numeric. Line width of the data point borders.
  #       star_dist: Numeric. Offset distance for significance stars above brackets.
  #     bracket_gap: Numeric. Vertical space (fraction of y-range) between data and first bracket.
  #         verbose: Logical. Print summary of significant contrasts to console.
  
  # --- 1. Calculate Estimated Marginal Means (95% only) ---
  f <- formula(formula)
  cf <- emmeans(model, f, type = "response", level = 0.95)
  dat <- data.frame(cf)
  
  # Robust Column Identification
  y_col_idx <- grep("prob|rate|response|emmean", names(dat), ignore.case = TRUE)[1]
  lci_col_idx <- grep("lower|LCL|2.5", names(dat), ignore.case = TRUE)[1]
  uci_col_idx <- grep("upper|UCL|97.5", names(dat), ignore.case = TRUE)[1]
  
  if (is.na(lci_col_idx) || is.na(uci_col_idx)) {
    stop("Could not automatically identify Confidence Interval columns.")
  }
  
  x_var_name <- names(dat)[1]
  x_dat <- dat[, 1]
  levels_x <- levels(x_dat)
  
  # --- 2. Calculate Contrasts ---
  sig_data <- NULL
  if (show_sig) {
    prs <- contrast(cf, method = "pairwise", adjust = "tukey")
    s_df <- as.data.frame(prs)
    sig_data <- s_df[s_df$p.value < 0.05, ]
  }
  
  if (verbose) {
    cat(sprintf("[%s] Significant contrasts: %d\n", 
                paste(format(formula), collapse=""),
                if(is.null(sig_data)) 0 else nrow(sig_data)))
  }
  
  # --- 3. Setup Plot Params ---
  if (length(col_vec) == 1) {
    cols <- rep(col_vec, length(levels_x))
  } else {
    cols <- rep(col_vec, length.out = length(levels_x))
  }
  
  # --- 4. Handle Canvas Expansion ---
  
  # Determine logical limits (What the axis shows)
  if (!is.null(ylim)) {
    logical_ylim <- ylim
    data_max <- ylim[2] # Use the manual max for bracket calculations
  } else {
    data_max <- max(dat[, uci_col_idx])
    logical_ylim <- c(0, data_max)
  }
  
  visual_range <- diff(logical_ylim)
  
  # Top Padding Calculation
  if (!is.null(top_buffer)) {
    pad_top <- top_buffer * visual_range
  } else {
    base_buff <- 0.1
    if (show_sig && !is.null(sig_data) && nrow(sig_data) > 0) {
      base_buff <- base_buff + 0.3
    }
    pad_top <- base_buff * visual_range
  }
  
  pad_bot <- 0.05 * visual_range
  canvas_ylim <- c(logical_ylim[1] - pad_bot, logical_ylim[2] + pad_top)
  
  x_lab <- if(exists("lab") && x_var_name %in% lab$v) lab[lab$v == x_var_name, 2] else x_var_name
  
  # --- 5. Draw Plot ---
  plot(1, type = "n", 
       xlim = c(0.5, length(levels_x) + 0.5), 
       ylim = canvas_ylim, 
       xaxt = "n", yaxt = "n", 
       xlab = "", ylab = "")
  
  # Titles
  mtext("Dispersal rate", side = 2, line = ylab_dist, cex = title_cex)
  mtext(x_lab, side = 1, line = xlab_dist, cex = title_cex)
  
  # Axes
  y_ticks <- pretty(logical_ylim)
  # Ensure we don't draw ticks outside the requested ylim if manual
  if(!is.null(ylim)) y_ticks <- y_ticks[y_ticks <= ylim[2]]
  y_ticks <- y_ticks[y_ticks >= 0] 
  
  axis(2, at = y_ticks, las = 1, cex.axis = axis_cex)
  axis(1, at = 1:length(levels_x), labels = levels_x, cex.axis = axis_cex, las = 2)
  
  # --- 6. Plot Points ---
  for (i in seq_along(levels_x)) {
    x_pos <- i
    y <- dat[i, y_col_idx]
    lci <- dat[i, lci_col_idx]
    uci <- dat[i, uci_col_idx]
    
    segments(x_pos, lci, x_pos, uci, lwd = 2, col = cols[i])
    points(x_pos, y, pch = 21, bg = "white", col = cols[i], cex = pt_cex, lwd = pt_lwd)
  }
  
  # --- 7. Draw Significance Brackets ---
  if (show_sig && !is.null(sig_data) && nrow(sig_data) > 0) {
    step_h <- visual_range * 0.05 
    base_y <- data_max + (visual_range * bracket_gap)
    
    for (r in 1:nrow(sig_data)) {
      parts <- unlist(strsplit(as.character(sig_data$contrast[r]), "\\s+[-/]\\s+"))
      if (length(parts) >= 2) {
        idx1 <- which(levels_x == trimws(parts[1]))
        idx2 <- which(levels_x == trimws(parts[2]))
        
        if (length(idx1) > 0 && length(idx2) > 0) {
          x1 <- idx1; x2 <- idx2
          
          lines(c(x1, x1, x2, x2), 
                c(base_y, base_y + step_h/3, base_y + step_h/3, base_y), 
                lwd = 1.5)
          
          p_val <- sig_data$p.value[r]
          sym <- ifelse(p_val < 0.001, "***", ifelse(p_val < 0.01, "**", "*"))
          
          text(mean(c(x1, x2)), base_y + step_h/2, labels = sym, cex = 1.2, pos = 3, offset = star_dist)
          
          base_y <- base_y + step_h * 1.5
        }
      }
    }
  }
}


mult <- 1.4

# Dimensions in millimetres:
w <- 173; h <- 80

# Dimensions in inches * multiplier:
wi <- round(mult * w / 25.4, 1)
hi <- round(mult * h / 25.4, 1)

cairo_pdf("Figures/Fig_S3.pdf", width = wi, height = hi, symbolfamily = "OpenSymbol")

op <- par(mfrow = c(1, 3), mar = c(11, 6, 3, 2))

# Panel A
plot_1_way(m, '~ cue', ylim = c(0, 0.15))
mtext("(a)", side = 3, line = 1, adj = -0.1, cex = 1.2)

# Panel B
plot_1_way(m, '~ env', top_buffer = 0.2, ylim = c(0, 0.15))
mtext("(b)", side = 3, line = 1, adj = -0.1, cex = 1.2)

# Panel C
plot_1_way(m, '~ host_spec', top_buffer = 0.2, ylim = c(0, 0.15))
mtext("(c)", side = 3, line = 1, adj = -0.1, cex = 1.2)

par(op)

dev.off()



# Figure 1 ----------------------------------------------------------------

# library(RColorBrewer)
# display.brewer.all(8, type = "qual")
# pal <- brewer.pal(8, "Dark2")

plot_2_way <- function(model, formula, 
                       col_vec = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666"), 
                       show_sig_top = TRUE, show_sig_bottom = TRUE, 
                       add_legend = TRUE, legend_pos = "topright",
                       legend_buffer = 0.45, 
                       bottom_buffer = NULL,
                       ylab_dist = 3.7, 
                       xlab_dist = 3.0, 
                       title_cex = 1.3, 
                       axis_cex = 1.4,  
                       pt_cex = 1.5,    
                       pt_lwd = 2.5,    
                       leg_cex = 1.3,   
                       star_dist_top = -0.1,
                       star_dist_bot = 0.4,
                       offset_spread = 0.2, 
                       ylim = NULL,
                       x_tick_labs = NULL, 
                       xlab_text = NULL,   
                       leg_title = NULL,
                       top_bracket_spacing = 1.5,
                       verbose = TRUE) {
  
  # Plots 2-way interactions from a model with EMMeans and significance brackets.
  #
  # Arguments:
  #           model: A model object supported by the 'emmeans' package.
  #         formula: A formula specifying the interaction (e.g. '~ x_axis | grouping').
  #         col_vec: Vector of colors for the grouping levels (default 'Dark2').
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
  #   star_dist_top: Numeric. Offset distance for significance stars above top brackets.
  #   star_dist_bot: Numeric. Offset distance for significance stars below bottom brackets.
  #   offset_spread: Numeric. Width of the spread of points around the x-tick.
  #            ylim: Vector c(min, max). Manual limits for the Y-axis.
  #     x_tick_labs: Character vector. Custom labels for X-axis ticks.
  #       xlab_text: Character string. Custom title for the X-axis.
  #       leg_title: Character string. Custom title for the Legend.
  # top_bracket_spacing: Numeric. Multiplier for vertical gap between top brackets (default 1.5).
  #         verbose: Logical. Print summary of significant contrasts to console.
  
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
  offsets <- seq(-offset_spread, offset_spread, length.out = n_groups)
  cols <- rep(col_vec, length.out = n_groups)
  full_pch_vec <- c(21, 22, 24, 23, 25, 21, 22, 24, 23, 25) 
  pchs <- rep(full_pch_vec, length.out = n_groups)
  
  # --- 4. Handle Canvas Expansion ---
  real_data_max <- max(dat[, 7])
  user_max <- if(!is.null(ylim)) ylim[2] else real_data_max
  calc_max <- max(real_data_max, user_max)
  
  logical_ylim <- c(0, calc_max)
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
    n_bottom_sigs <- if(!is.null(sig_data_bot)) nrow(sig_data_bot) else 0
    pad_bot <- (0.1 + (0.05 * n_bottom_sigs) + (0.05 * n_groups)) * visual_range
  }
  
  canvas_ylim <- c(0 - pad_bot, calc_max + pad_top)
  
  # Logic for X-Axis TITLE
  if (!is.null(xlab_text)) {
    final_x_title <- xlab_text
  } else {
    final_x_title <- if(exists("lab") && x_var_name %in% lab$v) lab[lab$v == x_var_name, 2] else x_var_name
  }
  
  # Logic for Legend TITLE
  if (!is.null(leg_title)) {
    by_lab <- leg_title
  } else {
    by_lab <- if(exists("lab") && by_var_name %in% lab$v) lab[lab$v == by_var_name, 2] else by_var_name
  }
  
  # --- 5. Draw Plot ---
  plot(1, type = "n", 
       xlim = c(0.5, nlevels(x_dat) + 0.5), 
       ylim = canvas_ylim, 
       xaxt = "n", yaxt = "n", 
       xlab = "", ylab = "", 
       las = 1)
  
  # Titles
  mtext("Dispersal rate", side = 2, line = ylab_dist, cex = title_cex)
  mtext(final_x_title, side = 1, line = xlab_dist, cex = title_cex)
  
  # Axes
  y_ticks <- pretty(c(0, user_max))
  if(!is.null(ylim)) y_ticks <- y_ticks[y_ticks <= ylim[2]]
  y_ticks <- y_ticks[y_ticks >= 0] 
  
  axis(2, at = y_ticks, las = 1, cex.axis = axis_cex)
  
  # Logic for X-Axis TICKS
  final_tick_labels <- if(!is.null(x_tick_labs)) x_tick_labs else levels_x
  
  if(length(final_tick_labels) != nlevels(x_dat)) {
    warning(paste("Length of x_tick_labs (", length(final_tick_labels), 
                  ") does not match levels of X variable (", nlevels(x_dat), 
                  "). Using default levels."))
    final_tick_labels <- levels_x
  }
  
  axis(1, at = 1:nlevels(x_dat), labels = final_tick_labels, cex.axis = axis_cex)
  
  # --- 6. Plot Points ---
  for (i in seq_along(levels_by)) {
    idx <- dat[, 2] == levels_by[i]
    x_pos <- as.numeric(x_dat[idx]) + offsets[i]
    y <- dat[idx, 3]
    lci <- dat[idx, 6]; uci <- dat[idx, 7]
    
    segments(x_pos, lci, x_pos, uci, lwd = 2, col = cols[i])
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
              text(mean(c(x1, x2)), base_y + step_h/2, labels = sym, cex = 1.2, pos = 3, offset = star_dist_top)
              
              # Increment with spacing control
              base_y <- base_y + (step_h * top_bracket_spacing)
            }
          }
        }
      }
    }
  }
  
  # --- 8. Draw BOTTOM Brackets ---
  if (show_sig_bottom && !is.null(sig_data_bot) && nrow(sig_data_bot) > 0) {
    current_y <- 0 
    for (gi in seq_along(levels_by)) {
      curr_group <- levels_by[gi]
      curr_sigs <- sig_data_bot[as.character(sig_data_bot[[by_var_name]]) == as.character(curr_group), ]
      
      if (nrow(curr_sigs) > 0) {
        current_y <- current_y - (step_h * 0.5)
        for (r in 1:nrow(curr_sigs)) {
          parts <- unlist(strsplit(as.character(curr_sigs$contrast[r]), "\\s+[-/]\\s+"))
          if (length(parts) >= 2) {
            idx1 <- which(levels_x == trimws(parts[1]))
            idx2 <- which(levels_x == trimws(parts[2]))
            if (length(idx1) > 0 && length(idx2) > 0) {
              x1 <- idx1 + offsets[gi]; x2 <- idx2 + offsets[gi]
              y_top <- current_y
              y_bottom <- current_y - (step_h * 0.4)
              lines(c(x1, x1, x2, x2), c(y_top, y_bottom, y_bottom, y_top), lwd = 1.5, col = cols[gi])
              p_val <- curr_sigs$p.value[r]
              sym <- ifelse(p_val < 0.001, "***", ifelse(p_val < 0.01, "**", "*"))
              text(mean(c(x1, x2)), y_bottom, labels = sym, cex = 1.2, col = cols[gi], pos = 1, offset = star_dist_bot) 
              current_y <- y_bottom - (step_h * 1.2) 
            }
          }
        }
      }
    }
  }
  
  # --- 9. Legend ---
  if (add_legend) {
    leg_ncol <- if (n_groups > 3) 2 else 1
    legend(legend_pos, legend = levels_by, col = cols, pt.bg = "white", pch = pchs, lwd = 2, title = by_lab, bty = "n", cex = leg_cex, ncol = leg_ncol, text.width = max(strwidth(levels_by, cex = leg_cex)) * 1.5)
  }
}


mult <- 1.4

# Dimensions in millimetres:
w <- 173; h <- 70

# Dimensions in inches * multiplier:
wi <- round(mult * w / 25.4, 1)
hi <- round(mult * h / 25.4, 1)

cairo_pdf("Figures/Fig_1.pdf", width = wi, height = hi, symbolfamily = "OpenSymbol")

op <- par(mfrow = c(1, 3), mar = c(5.5, 6, 3, 2))

# Panel A
plot_2_way(m, '~ env | cue', legend_buffer = 0.6, bottom_buffer = 0.3, leg_title = "Kairomone", xlab_text = "Current environment")
mtext("(a)", side = 3, line = 1, adj = -0.1, cex = 1.2)

# Panel B
plot_2_way(m, '~ host_spec | cue', legend_buffer = 0.3, bottom_buffer = 0.25, leg_title = "Kairomone", xlab_text = "Host specialisation")
mtext("(b)", side = 3, line = 1, adj = -0.1, cex = 1.2)

# Panel C
plot_2_way(m, '~ env | host_spec', legend_buffer = 0.1, bottom_buffer = 0.2, leg_title = "Host specialisation", xlab_text = "Current environment")
mtext("(c)", side = 3, line = 1, adj = -0.1, cex = 1.2)

par(op)

dev.off()


# Figure 2 ----------------------------------------------------------------

library(RColorBrewer)
display.brewer.all(type = "qual")
pal <- brewer.pal(4, "Set1")
col <- adjustcolor(pal, alpha = 0.5)
colt <- adjustcolor(pal, alpha = 0.1)


# Data ----

load("data/data.RData")
data <- subset(data, n_kairo > 0)  # Zero is equivalent to control in the variable ‘cue’.
data <- subset(data, plant == "W") # No mixes for env==‘unknown’.
data <- transform(data, group = interaction(cue, host_spec, sep = "_"), q = D / N, SN = n_familiar / n_kairo, regime = NULL, plant = NULL, env = NULL)
data <- droplevels(data)
summary(data)


# Models ----

# No. of kairomones
mnk <- glmmTMB(D/N ~ 0 + group + group:n_kairo + (1|line), weights = N, data, family = betabinomial)

# SN ratio
msn <- glmmTMB(D/N ~ 0 + host_spec + host_spec:SN + (1|line), weights = N, data, family = betabinomial)


# 1. Prepare Data for the Graph ----

# (a) Number of kairomones (Model mnk)
group_levels <- levels(data$group)
n_kairo_seq <- seq(0.9, 3.1, length.out = 100)

# Create prediction grid for (a)
nd_a <- expand.grid(group = group_levels, n_kairo = n_kairo_seq, line = NA, N = 10)

# Exclude unobserved range for Unfamiliar Generalists
idx_a <- (nd_a$group == "Unfamiliar_Generalist") & (nd_a$n_kairo > 2.1)
nd_a <- nd_a[!idx_a, ]

# Predictions and Confidence Intervals for (a)
pred_a <- predict(mnk, nd_a, se = TRUE)
fit_a <- pred_a$fit
se_a <- pred_a$se.fit
q_a <- -qt(0.025, df.residual(mnk))

lw_a <- fit_a - q_a * se_a
up_a <- fit_a + q_a * se_a

# Transform to probability scale
p_a <- plogis(cbind(fit = fit_a, lw = lw_a, up = up_a))
nd_a <- cbind(nd_a, p_a)

# (b) Signal-to-noise ratio (Model msn)
host_spec_levels <- levels(data$host_spec)
SN_seq <- seq(0, 1, length.out = 100)


# Create prediction grid for (b)
nd_b <- expand.grid(host_spec = host_spec_levels, SN = SN_seq, line = NA, N = 10)

# Predictions and Confidence Intervals for (b)
pred_b <- predict(msn, nd_b, se = TRUE)
fit_b <- pred_b$fit
se_b <- pred_b$se.fit
q_b <- -qt(0.025, df.residual(msn))

lw_b <- fit_b - q_b * se_b
up_b <- fit_b + q_b * se_b

# Transform to probability scale
p_b <- plogis(cbind(fit = fit_b, lw = lw_b, up = up_b))
nd_b <- cbind(nd_b, p_b)


# 2. Plotting ----

cex.lab <- 1.3
ylim <- c(0, 0.25)

mult <- 1.4

# Dimensions in millimetres:
w <- 110; h <- 70

# Dimensions in inches * multiplier:
wi <- round(mult * w / 25.4, 1)
hi <- round(mult * h / 25.4, 1)

cairo_pdf("Figures/Fig_2.pdf", width = wi, height = hi, symbolfamily = "OpenSymbol")


op <- par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2.5, 1))

# Panel (a) ----
plot(fit ~ n_kairo, nd_a, type = "n", 
     xlim = extendrange(1:3), 
     ylim = ylim, 
     xlab = "No. of kairomones", 
     ylab = "Dispersal rate", 
     cex.lab = cex.lab,
     cex.axis = 0.8,
     xaxt = "n")
axis(1, at = 1:3)

for (i in seq_along(group_levels)) {
  g <- group_levels[i]
  dat <- subset(nd_a, group == g)
  polygon(x = c(dat$n_kairo, rev(dat$n_kairo)), 
          y = c(dat$lw, rev(dat$up)), 
          border = NA, col = colt[i])
  lines(dat$n_kairo, dat$fit, col = col[i], lwd = 3)
}

# Legend for (a)
# Cleaning names for the legend (removing underscores)
legend_labels <- gsub("_", " ", group_levels)
legend("topright", legend = legend_labels, 
       col = col[1:length(group_levels)], 
       lwd = 3, bty = "n", cex = 0.7)

mtext("(a)", side = 3, line = 1, adj = -0.1, cex = 1.2)


# --- Panel (b) ---
plot(fit ~ SN, nd_b, type = "n", 
     ylim = ylim, 
     xlab = "Signal-to-Noise Ratio", 
     ylab = "", 
     cex.axis = 0.8,
     cex.lab = cex.lab)

for (i in seq_along(host_spec_levels)) {
  h <- host_spec_levels[i]
  dat <- subset(nd_b, host_spec == h)
  polygon(x = c(dat$SN, rev(dat$SN)), 
          y = c(dat$lw, rev(dat$up)), 
          border = NA, col = colt[i]) # Re-using colors 1 and 2
  lines(dat$SN, dat$fit, col = col[i], lwd = 3)
}

# Legend for (b)
legend("topleft", legend = host_spec_levels, 
       col = col[1:length(host_spec_levels)], 
       lwd = 3, bty = "n", cex = 0.7)

mtext("(b)", side = 3, line = 1, adj = -0.1, cex = 1.2)

par(op)

dev.off()


# Figure S4 ---------------------------------------------------------------

# Data ----
load("data/data.RData")
data <- subset(data, n_kairo < 2)
data <- transform(data, W = NULL, B = NULL, O = NULL, S = NULL)
data <- droplevels(data)
summary(data)

# Model ----
m <- glmmTMB(D/N ~ host_spec * plant * kairo + (1|line), weights = N, data, family = betabinomial)
summary(m)
car::Anova(m)


# Figure ----

mult <- 1.5

# Dimensions in millimetres:
w <- 173; h <- 80

# Dimensions in inches * multiplier:
wi <- round(mult * w / 25.4, 1)
hi <- round(mult * h / 25.4, 1)

cairo_pdf("Figures/Fig_S4.pdf", width = wi, height = hi, symbolfamily = "OpenSymbol")
op <- par(mfrow = c(1, 3), mar = c(5.5, 6, 3, 2))

# Panel A
plot_2_way(m, '~ plant | kairo', legend_buffer = 0.4, leg_title = "Kairomone", xlab_text = "Current environment")
mtext("(a)", side = 3, line = 1, adj = -0.1, cex = 1.2)

# Panel B
plot_2_way(m, '~ host_spec | kairo', leg_title = "Kairomone", xlab_text = "Host specialisation", bottom_buffer = 0.4, top_bracket_spacing = 2.3, legend_buffer = 0.5)
mtext("(b)", side = 3, line = 1, adj = -0.1, cex = 1.2)

# Panel C
plot_2_way(m, '~ host_spec | plant', leg_title = "Current environment", xlab_text = "Host specialisation", bottom_buffer = 0.3, legend_buffer = 0.1)
mtext("(c)", side = 3, line = 1, adj = -0.1, cex = 1.2)

par(op)

dev.off()

