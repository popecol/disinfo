
# Code for the paper:
# Zalewska K., Skoracka A., Bonte D., Puchalska E., Lewandowski M., Kuczynski L.
# Is passive dispersal informed? - Experimental evidence for decision-making in phytophagous arthropods

# Experiment 2
# Dispersal in response to the signal noise through a mixture of kairomones from the target environment 


# Setup -------------------------------------------------------------------

# devtools::install_github("glmmTMB/glmmTMB/glmmTMB")

library(glmmTMB)
library(DHARMa)
library(emmeans)
library(ggplot2)


# Data --------------------------------------------------------------------

load("data/data.RData")

data <- subset(data, n_kairo > 0)  # Zero is equivalent to control in the variable ‘cue’.
data <- subset(data, plant == "W") # No mixes for env==‘unknown’.
data <- transform(data, plant = NULL, env = NULL, group = interaction(cue, host_spec, sep = "_"))
data <- droplevels(data)

summary(data)

#      line: Evolutionary lineage identifier.
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
#         q: Dispersal rate (q=D/N)
#     group: Combination of 'cue' and 'host_spec'.

summary(data$kairo)
table(data$cue, data$n_kairo)
table(data$host_spec, data$n_kairo)
table(data$group, data$n_kairo)


# Beta-binomial GLMM ------------------------------------------------------

# ANCOVA, the 'means' parameterisation:
m <- glmmTMB(D/N ~ 0 + group + group:n_kairo + (1|line), weights = N, data, family = betabinomial)
summary(m)
sr <- simulateResiduals(m); plot(sr)
car::Anova(m)

# Tests for slopes
at <- list(n_kairo = 1:3)
emtrends(m, pairwise ~ group, var = "n_kairo")
p <- emmip(m, group ~ n_kairo, CIs = TRUE, type = "response", at = at, xlab = "No. of kairomones", ylab = "Dispersal rate")
# p <- emmip(m, group ~ n_kairo, type = "response", cov.reduce = range, xlab = "No. of kairomones", ylab = "Dispersal rate")
p + theme(text = element_text(size = 18))
