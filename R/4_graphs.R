
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


# Figure 1 ----------------------------------------------------------------

# Main effects of kairomones (a), current environment (b) and host specialisation (c) on dispersal rates. 

# Main effects
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
  mtext(paste0("(", letters[i], ")"), side = 2, line = 2.5, padj = -5.5, cex = 1.5, las = 1)
}
par(op)


# Figure 2 ----------------------------------------------------------------

# 2-way interactions ----

plot_2_way <- function(model, formula) {
  f <- formula(formula)
  cf <- emmeans(model, f, type = "response", level = 0.95)
  cf50 <- emmeans(model, f, type = "response", level = 0.5)
  dat <- data.frame(print(cf))
  dat50 <- data.frame(print(cf50))
  xnam <- names(dat)[1]
  
  x_dat <- dat[, 1]
  x <- seq(nlevels(x_dat))
  by <- levels(dat[, 2])
  
  ylim <- range(dat[, 6:7])
  ylim <- extendrange(ylim, f = 0.1)
  # ylim <- c(0.03, 0.15)
  
  op <- par(mfrow = c(1, 3), mar = c(10, 5, 3, 2))
  for (i in seq(length(by))) {
    xlim <- extendrange(x, f = ifelse(k == 3, 0.3, 0.7))
    idx <- dat[, 2] == by[i]
    y <- dat[idx, 3]
    lci <- dat[idx, 6]; uci <- dat[idx, 7]
    lci50 <- dat50[idx, 6]; uci50 <- dat50[idx, 7]
    plot(x, y, xlim = xlim, ylim = ylim, xaxt = "n", xlab = "", ylab = "Dispersal rate", type = "n", cex.lab = 2, main = by[i], cex.main = 1.7)
    axis(1, at = x, labels = levels(x_dat), cex.axis = 1.5, las = 2)
    segments(x, lci50, x, uci50, lwd = 3)
    segments(x, lci, x, uci, lwd = 1)
    points(x, y, cex = 3, col = "white", pch = 20)
    points(x, y, cex = 1.3)
    mtext(paste0("(", letters[i], ")"), side = 2, line = 2.5, padj = -6.5, cex = 1.5, las = 1)
  }
  adj <- ifelse(length(by) == 2, 0.3, 0.5)
  mtext(lab[lab$v == xnam, 2], 1, -2, cex = 1.5, outer = TRUE, adj = adj)
  par(op)
}

plot_2_way(m, "~ env | cue")
plot_2_way(m, "~ host_spec | cue")
plot_2_way(m, "~ env | host_spec")


# Figure 5 ----------------------------------------------------------------

library(RColorBrewer)
display.brewer.all(type = "qual")
pal <- brewer.pal(4, "Set1")
col <- adjustcolor(pal, alpha = 0.5)
colt <- adjustcolor(pal, alpha = 0.1)


hist(data$N); summary(data$N)

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


