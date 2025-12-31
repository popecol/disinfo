
# Code for the paper:
# Zalewska K., Skoracka A., Bonte D., Puchalska E., Lewandowski M., Kuczynski L.
# Is passive dispersal informed? - Experimental evidence for decision-making in phytophagous arthropods

# Experiment 2
# Dispersal in response to the signal noise through a mixture of kairomones from the target environment 


# Setup -------------------------------------------------------------------

# install.packages("glmmTMB", type = "source")
# remotes::install_github("glmmTMB/glmmTMB/glmmTMB")
library(glmmTMB)
library(emmeans)
library(ggplot2)


# Data --------------------------------------------------------------------

load("data/data.RData")
data <- subset(data, n_kairo > 0)  # Zero is equivalent to control in the variable ‘cue’.
data <- subset(data, plant == "W") # No mixes for env==‘unknown’.
data <- transform(data, group = interaction(cue, host_spec, sep = "_"), q = D / N, SN = n_familiar / n_kairo, regime = NULL, plant = NULL, env = NULL)
data <- droplevels(data)
summary(data)

#       line: Evolutionary lineage identifier.
#      kairo: Kairomones delivered (control=no kairomones, W=wheat, B=barley, 
#             S=smooth brome, O=oats, and their combinations)
#          W: Was the wheat kairomone delivered? (0=no, 1=yes)
#          B: Was the barley kairomone delivered? (0=no, 1=yes)
#          O: Was the oats kairomone delivered? (0=no, 1=yes)
#          S: Was the smooth brome kairomone delivered? (0=no, 1=yes)
#    n_kairo: The overall number of kairomones delivered.
#          N: The number of individuals at the beginning of the experiment.
#          R: The number of individuals that remained on the plant (i.e. did not disperse).
#  host_spec: Host specialisation (Generalist, Specialist).
#          D: The number of dispersers (D=N-R).
#        cue: Have any of the kairomones ever been encountered before? (Control=no kairomone, 
#             Familiar=yes, Unfamiliar=no)
# n_familiar: The number of familiar kairomones delivered.
#      group: Combination of 'cue' and 'host_spec'.
#          q: Dispersal rate (q=D/N)
#         SN: Signal-to-noise ratio (SN=n_familiar/n_kairo)

summary(data$kairo)
table(data$cue, data$n_kairo)
table(data$host_spec, data$n_kairo)
table(data$group, data$n_kairo)


# Beta-binomial GLMM ------------------------------------------------------
# ANCOVA, the 'means' parameterisation.

# The number of cairomones ----

mnk <- glmmTMB(D/N ~ 0 + group + group:n_kairo + (1|line), weights = N, data, family = betabinomial)

summary(mnk)
car::Anova(mnk)

# Tests for slopes
emtrends(mnk, pairwise ~ group, var = "n_kairo")
p <- emmip(mnk, group ~ n_kairo, CIs = TRUE, type = "response", at = list(n_kairo = sort(unique(data$n_kairo))), xlab = "SN ratio", ylab = "Dispersal rate")
p + theme(text = element_text(size = 18))


# Signal-to-noise ratio ----
msn <- glmmTMB(D/N ~ 0 + host_spec + host_spec:SN + (1|line), weights = N, data, family = betabinomial)

summary(msn)
car::Anova(msn)

# Tests for slopes
emtrends(msn, pairwise ~ host_spec, var = "SN")
p <- emmip(msn, host_spec ~ SN, CIs = TRUE, type = "response", at = list(SN = sort(unique(data$SN))), xlab = "SN ratio", ylab = "Dispersal rate")
p + theme(text = element_text(size = 18))

