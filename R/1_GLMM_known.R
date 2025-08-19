
# Code for the paper:
# Zalewska K., Skoracka A., Bonte D., Puchalska E., Lewandowski M., Kuczynski L.
# Is passive dispersal informed? - Experimental evidence for decision-making in phytophagous arthropods  

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
  # Emmeans Ploting Helper.
  
  require(ggplot2)
  f <- update(formula, pairwise ~ .)
  print(emmeans(obj, f, type = "response")$contrasts)
  p <- emmip(obj, formula, CIs = TRUE, type = "response", ylab = "Dispersal rate", ...)
  p + theme(text = element_text(size = size))
}


# Data --------------------------------------------------------------------

load("data/data.RData")
summary(data)

#      line: Evolutionary lineage identifier.
#     plant: Host plant (W=wheat, B=barley, S=smooth brome).
#     kairo: Kairomones delivered (control=no kairomones, W=wheat, B=barley, S=smooth brome, O=oats, and their combinations)
# host_spec: Host specialisation (Generalist, Specialist).
#         W: Was the wheat kairomone delivered? (0=no, 1=yes)
#         B: Was the barley kairomone delivered? (0=no, 1=yes)
#         O: Was the oats kairomone delivered? (0=no, 1=yes)
#         S: Was the smooth brome kairomone delivered? (0=no, 1=yes)
#   n_kairo: The overall number of kairomones delivered.
#         N: The number of individuals at the beginning of the experiment.
#         R: The number of individuals that remained on the plant (i.e. did not disperse).
#         D: The number of dispersers (D=N-R).
#       cue: Have any of the kairomones ever been encountered before? (Control=no kairomone, Familiar=yes, Unfamiliar=no)
#       env: Have the current environment ever been encountered before? (Familiar=yes, Unfamiliar=no)
#         q: Dispersal rate (q=D/N)


# Beta-binomial GLMM ------------------------------------------------------

m <- glmmTMB(D/N ~ cue * env * host_spec + (1|line), weights = N, data, family = betabinomial)
summary(m)
sr <- simulateResiduals(m); plot(sr)
car::Anova(m)

# Marginal effects ----
p1 <- em(m, ~ cue, 16)
p2 <- em(m, ~ env, 16)
p3 <- em(m, ~ host_spec, 16)
cowplot::plot_grid(p1, p2, p3, ncol = 3, labels = "AUTO")


# 2-way interactions ----
em(m, ~ env | cue)
# em(m, ~ cue | env)

em(m, ~ host_spec | cue)
# em(m, ~ cue | host_spec)

em(m, ~ env | host_spec)
# em(m, ~ host_spec | env)


# 3-way interaction ----
em(m, ~ env | cue | host_spec)
# em(m, ~ cue | env | host_spec)
