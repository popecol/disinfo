
# Code for the paper:
# Zalewska, K., A.Skoracka, D.Bonte, E.Puchalska, M.Lewandowski, and L.Kuczyński. 2026. 
# “Niche Breadth and Olfactory Context Shape Informed Passive Dispersal.” 
# Ecology Letters 29, no. 4: e70373. https://doi.org/10.1111/ele.70373.

# Experiment 1 
# Dispersal in response to a single kairomones in the context of:
# 1. Kairomone (Control, Familiar or Unfamiliar),
# 2. Current environment (Familiar or Unfamiliar),
# 3. Host specialisation (Generalist or Specialist).


# Setup -------------------------------------------------------------------

# devtools::install_github("glmmTMB/glmmTMB/glmmTMB")

library(glmmTMB)
library(DHARMa)
library(emmeans)

em <- function(obj, formula, size = 18, ...) {
  
  # Emmeans Plotting Helper.
  require(ggplot2)
  
  f <- update(formula, pairwise ~ .)
  
  # 1. Calculate Contrasts (Odds Ratios)
  # type = "response": back-transforms to Odds Ratios
  contr <- emmeans(obj, f, type = "response")$contrasts
  
  # 2. Print Summary with BOTH Intervals and P-values
  # infer = c(TRUE, TRUE): Requests both Confidence Intervals and Hypothesis Tests
  print(summary(contr, infer = c(TRUE, TRUE)))
  
  # 3. Generate Plot
  p <- emmip(obj, formula, CIs = TRUE, type = "response", ylab = "Dispersal rate", ...)
  p + theme(text = element_text(size = size))
}


# Data --------------------------------------------------------------------

load("data/data.RData")
summary(data)

#     line: Evolutionary lineage identifier.
#     plant: Host plant (W=wheat, B=barley, S=smooth brome).
#     kairo: Kairomones delivered (control=no kairomones, W=wheat, B=barley, 
#            S=smooth brome, O=oats, and their combinations)
# host_spec: Host specialisation (Generalist, Specialist).
#         W: Was the wheat kairomone delivered? (0=no, 1=yes)
#         B: Was the barley kairomone delivered? (0=no, 1=yes)
#         O: Was the oats kairomone delivered? (0=no, 1=yes)
#         S: Was the smooth brome kairomone delivered? (0=no, 1=yes)
#   n_kairo: The overall number of kairomones delivered.
#         N: The number of individuals at the beginning of the experiment.
#         R: The number of individuals that remained on the plant (i.e. did not disperse).
#         D: The number of dispersers (D=N-R).
#       cue: Have any of the kairomones ever been encountered before? (Control=no kairomone,
#            Familiar=yes, Unfamiliar=no)
#       env: Have the current environment ever been encountered before? (Familiar=yes, Unfamiliar=no)
#         q: Dispersal rate (q=D/N)


# Beta-binomial GLMM ------------------------------------------------------

m1 <- glmmTMB(D/N ~ cue * env * host_spec + (1|line), weights = N, data, family = betabinomial)
summary(m1)
sr <- simulateResiduals(m1); plot(sr)
car::Anova(m1)

# Marginal effects ----
p1 <- em(m1, ~ cue, 16)
p2 <- em(m1, ~ env, 16)
p3 <- em(m1, ~ host_spec, 16)
cowplot::plot_grid(p1, p2, p3, ncol = 3, labels = "AUTO")


# 2-way interactions ----
em(m1, ~ env | cue)
# em(m1, ~ cue | env)

em(m1, ~ host_spec | cue)
# em(m1, ~ cue | host_spec)

em(m1, ~ env | host_spec)
# em(m1, ~ host_spec | env)


# 3-way interaction ----
em(m1, ~ env | cue | host_spec)
# em(m1, ~ cue | env | host_spec)
