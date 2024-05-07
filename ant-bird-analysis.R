library(tidyverse)
library(MASS)
library(pscl)
library(multcomp)
library(emmeans)
library(broom)

# Data cleaning and transformation steps
# For example:
# bird_ind <- bird %>%
#   mutate(
#        rsa = (Rinter > 0),
#        roa = (Rintra > 0),
#        isa = (Iinter > 0),
#        ioa = (Iintra > 0),
#     )
# write_csv(bird_ind, "bird_ind.csv")

# apply factor
hierarchy_fct <- c("Phleg", "M. fort", "Rheg", "D. mer", "G. sal") # most dominant to least dominant
bird_ind <- read_csv("data/bird_ind.csv") %>% mutate(species = factor(species, hierarchy_fct))

# load in dataset
head(bird_ind)
levels(bird_ind$species)

# Best model selection
# nbinom, poisson, zero-infl and quasi forms of these models, as candidates 

# Fit Poisson model
pois_full <- glm(Successes ~ species + rsa + roa + isa + ioa + species:rsa + species:roa + species:isa + species:ioa, family = "poisson", data = bird_ind)
summary(pois_full)

# Fit Quasi Poisson
quasi_full <- glm(Successes ~ species + rsa + roa + isa + ioa + species*rsa + species*roa + species*isa + species*ioa, family = "quasipoisson", data = bird_ind)
summary(quasi_full)
#   Dispersion parameter is not much larger than 1, meaning there isn't much overdispersion to account for. 
#   If we multiply each of the standard errors of the original poisson model to the sqrt(dispersion parameter),
#   there is not much of a change, as can be seen in the model when compared to the poisson model. 

# Fit Negative Binomial model
nbinom_full <- glm.nb(Successes ~ species + rsa + roa + isa + ioa + species:rsa + species:roa + species:isa + species:ioa, data = bird_ind)
summary(nbinom_full)

# Ruled out Zero-Inflated Negative Binomial and Poisson models due to species*ioa and species*isa having NA values

# Negative Binomial or Poisson?
# Upon looking at the shape of data:
bird %>%
  group_by(species) %>%
  summarize(
    meanRinter = mean(Rinter),
    meanRintra = mean(Rintra),
    meanIinter = mean(Iinter),
    meanIntra = mean(Iintra),
    sdRinter = sd(Rinter),
    sdRintra = sd(Rintra),
    sdIinter = sd(Iinter),
    sdIintra = sd(Iintra),
    n = n()
  )

#   It is revealed that Standard deviation counts are higher than means suggesting
#   our mean and variances are not equal....Poisson might not be the best model. 
#   Negative binomial does account for overdispersion where the variance is allowed
#   to be greater than the mean. Quasi poisson does not strike me as any different 
#   from my poisson model.


# Fitting Best Model: Testing Predictor Significance
nbinom_full <- glm.nb(Successes ~ species + rsa + roa + isa + ioa + species*rsa + species*roa + species*isa + species*ioa, data = bird_ind)
summary(nbinom_full)
plot(nbinom_full)

# use stepwise function to determine model with lowest AIC
step.model <- step(nbinom_full, direction="both")

# best model candidates w/ lowest AICs:
bestmod1 <- glm.nb(Successes ~ species + rsa + ioa + species*rsa, data = bird_ind)
bestmod2 <- glm.nb(Successes ~ species + rsa + species*rsa, data = bird_ind)
summary(bestmod1)
summary(bestmod2)

# since there is not a large difference in AIC between the two candidates, we choose the simpler option: bestmod2

# assess goodness of fit for model
gof.pvalue = 1 - pchisq(bestmod2$deviance, bestmod2$df.residual)
gof.pvalue
# model does not fit the data very well...


# Now we continue to Best Model Interpretation:
# Rename
bestmod <- bestmod2

# evaluate series of pairwise comparisons across interaction term rsa; not really sure if this is important
# pairs(emmeans(bestmod, ~ rsa | species))

# visualize the log mean number of successes (linear prediction) for species receiving same species aggression (rsa)
emmip_successes <- emmip(bestmod, species ~ rsa | rsa) + labs(y = "Log Mean Number of Successes", x = "Levels of Rsa")

ggsave("data/log_mean_rsa.jpg", emmip_successes, width = 8, height = 6, units = "in")

# calculate estimated marginal means (EMMs) for pairwise comparisons between species w respect to rsa (grouping var)
emmeans(bestmod, pairwise ~ species | rsa)

# calculate marginal means for the interaction between species and rsa
bird.emm <- emmeans(bestmod, ~ species * rsa)
# create pairwise contrasts between rsa levels, adjust for multiple testing using the multivariate t method
bird.cont <- contrast(bird.emm, "consec", simple = "each", combine = TRUE, adjust = "mvt")

# Estimate marginal means for the interaction
emm_int <- emmeans(bestmod, specs = ~ species * rsa)

# Pairwise comparisons with adjustment for multiple testing
comp <- contrast(emm_int, method = "pairwise", adjust = "bonferroni")
summary(comp)

confint_successes <- 
  plot(comp) + theme_bw() + geom_vline(xintercept = 0, color = "red") + 
  labs(y = "Estimated Log-Count Difference", x = "95% Confint Estimate")
confint_successes
ggsave("data/confint_successes.jpg", plot = confint_successes, width = 7, height = 10, dpi = 300)

# visualization of confidence intervals; any interval that contains 0 is not significant 
plot(comp)



# now let's do it with attempts

# Fit Poisson model
pois_full <- glm(Successes ~ species + rsa + roa + isa + ioa + species:rsa + species:roa + species:isa + species:ioa, family = "poisson", data = bird_ind)
summary(pois_full)

# Fit Quasi Poisson
quasi_full <- glm(Successes ~ species + rsa + roa + isa + ioa + species*rsa + species*roa + species*isa + species*ioa, family = "quasipoisson", data = bird_ind)
summary(quasi_full)
#   Dispersion parameter is not much larger than 1, meaning there isn't much overdispersion to account for. 
#   If we multiply each of the standard errors of the original poisson model to the sqrt(dispersion parameter),
#   there is not much of a change, as can be seen in the model when compared to the poisson model. 

# Fit Negative Binomial model
nbinom_full <- glm.nb(Successes ~ species + rsa + roa + isa + ioa + species:rsa + species:roa + species:isa + species:ioa, data = bird_ind)
summary(nbinom_full)

# Ruled out Zero-Inflated Negative Binomial and Poisson models due to species*ioa and species*isa having NA values

# Negative Binomial or Poisson?
# Upon looking at the shape of data:
bird %>%
  group_by(species) %>%
  summarize(
    meanRinter = mean(Rinter),
    meanRintra = mean(Rintra),
    meanIinter = mean(Iinter),
    meanIntra = mean(Iintra),
    sdRinter = sd(Rinter),
    sdRintra = sd(Rintra),
    sdIinter = sd(Iinter),
    sdIintra = sd(Iintra),
    n = n()
  )

#   It is revealed that Standard deviation counts are higher than means suggesting
#   our mean and variances are not equal....Poisson might not be the best model. 
#   Negative binomial does account for overdispersion where the variance is allowed
#   to be greater than the mean. Quasi poisson does not strike me as any different 
#   from my poisson model.


# Fitting Best Model: Testing Predictor Significance
nbinom_full_a <- glm.nb(attempts ~ species + rsa + roa + isa + ioa + species*rsa + species*roa + species*isa + species*ioa, data = bird_ind)
summary(nbinom_full_a)
plot(nbinom_full_a)

# use stepwise function to determine model with lowest AIC
step.model <- step(nbinom_full_a, direction="both")

# best model candidates w/ lowest AICs:
bestmod1_a <- glm.nb(attempts ~ species + roa + ioa + species*roa, data = bird_ind)
bestmod2_a <- glm.nb(attempts ~ species + rsa + roa + isa + ioa + species:roa, data = bird_ind)
summary(bestmod1_a)
summary(bestmod2_a)


anova(bestmod1_a, bestmod2_a, test = "Chisq")

# Since there's no siginificant difference between the two models' AIcs, and 
# there is no evidence to reject that the models are any different, we're going to chose the first model.

bestmod_a <- bestmod1_a

# assess goodness of fit for model
gof.pvalue = 1 - pchisq(bestmod_a$deviance, bestmod_a$df.residual)
gof.pvalue
# low p-value; model does not fit the data very well...


# Now we continue to Best Model Interpretation:

# evaluate series of pairwise comparisons across interaction term roa; not really sure if this is important
# pairs(emmeans(bestmod, ~ rsa | species))

# visualize the log mean number of successes (linear prediction) for species receiving other species aggression (roa)
emmip(bestmod_a, species ~ roa | roa)

# calculate estimated marginal means (EMMs) for pairwise comparisons between species w respect to roa (grouping var)
emmeans(bestmod_a, pairwise ~ species | roa)

# calculate marginal means for the interaction between species and roa
bird.emm <- emmeans(bestmod_a, ~ species * roa)
# create pairwise contrasts between roa levels, adjust for multiple testing using the multivariate t method
bird.cont <- contrast(bird.emm, "consec", simple = "each", combine = TRUE, adjust = "mvt")

# Estimate marginal means for the interaction
emm_int <- emmeans(bestmod_a, specs = ~ species * roa)

# Pairwise comparisons with adjustment for multiple testing
comp <- contrast(emm_int, method = "pairwise", adjust = "bonferroni")
summary(comp)

# visualization of confidence intervals; any interval that contains 0 is not significant 
confint_attempts <- 
  plot(comp) + theme_bw() + geom_vline(xintercept = 0, color = "red") + 
  labs(y = "Estimated Log-Count Difference", x = "95% Confint Estimate")
confint_attempts
ggsave("data/confint_attempts.jpg", plot = confint_attempts, width = 7, height = 10, dpi = 300)

# let's do the same with our other term, roa or whether a bird recieved outer species aggression;
# visualize the log mean number of successes (linear prediction) for roa
emmip_attempts <- emmip(bestmod_a, species ~ roa | roa) + labs(y = "Log Mean Number of Attempts", x = "Levels of Roa" )

ggsave("data/log_mean_roa.jpg", emmip_attempts, width = 8, height = 6, units = "in")

# calculate estimated marginal means (EMMs) for pairwise comparisons between species w respect to ioa (grouping var)
emmeans(bestmod_a, pairwise ~ species | roa)



# Tidy up the summary output
tidy_output <- tidy(bestmod)

tidy_output_a <- tidy(bestmod_a)

# View the tidied coefficient estimates and standard errors


tidy_output

