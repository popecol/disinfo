
# Code for the paper:
# Zalewska K., Skoracka A., Bonte D., Puchalska E., Lewandowski M., Kuczynski L.
# Is passive dispersal informed? - Experimental evidence for decision-making in phytophagous arthropods

# Experiment 2
# Dispersal in response to the signal noise through a mixture of kairomones from the target environment.

# The analysis is the same as in the '1_GLMM_known.R' script, but in this case the model has been fitted to species-level kairomones.


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

# Removal of mixtures
data <- subset(data, n_kairo < 2)
data <- transform(data, W = NULL, B = NULL, O = NULL, S = NULL)
data <- droplevels(data)
summary(data)

#      line: Evolutionary lineage identifier.
#     plant: Host plant (W=wheat, B=barley, S=smooth brome).
#     kairo: Kairomones delivered (control=no kairomones, W=wheat, B=barley, S=smooth brome, O=oats, and their combinations)
# host_spec: Host specialisation (Generalist, Specialist).
#   n_kairo: The overall number of kairomones delivered.
#         N: The number of individuals at the beginning of the experiment.
#         R: The number of individuals that remained on the plant (i.e. did not disperse).
#         D: The number of dispersers (D=N-R).
#       cue: Have any of the kairomones ever been encountered before? (Control=no kairomone, Familiar=yes, Unfamiliar=no)
#       env: Have the current environment ever been encountered before? (Familiar=yes, Unfamiliar=no)
#         q: Dispersal rate (q=D/N)

table(data$plant, data$kairo)


# Beta-binomial GLMM ------------------------------------------------------

m <- glmmTMB(D/N ~ host_spec * plant * kairo + (1|line), weights = N, data, family = betabinomial)
summary(m)
car::Anova(m)
sr <- simulateResiduals(m); plot(sr)

# Marginal effects
p1 <- em(m, ~ kairo, 16) + xlab("Kairomone")
p2 <- em(m, ~ plant, 16) + xlab("Current environment")
p3 <- em(m, ~ host_spec, 16) + xlab("Specialisation")
cowplot::plot_grid(p1, p2, p3, ncol = 3, labels = "AUTO")

# 2-way interactions
em(m, ~ plant | kairo) + xlab("Current environment")
em(m, ~ host_spec | kairo) + xlab("Specialisation")
em(m, ~ host_spec | plant) + xlab("Specialisation")

# 3-way interaction ----
em(m, ~ host_spec | plant | kairo) + xlab("Specialisation")
