# ------------------------------------------------------------------------------
# Model convergence diagnostics and fit evaluation
# ------------------------------------------------------------------------------

# Load libraries and reset workspace
source("Scripts_R/0-Load libraries.R")
rm(list = ls())

### Check MCMC convergence diagnostics
# convert models to coda object
mPA <- get(load("../Models/model_mpa_thin_4000_samples_250_chains_4.RData"))
mABU <- get(load("../Models/model_abu_thin_4000_samples_250_chains_4.RData"))

mpost_PA <- convertToCodaObject(mPA)
mpost_ABU <- convertToCodaObject(mABU)

### Convergence diagnostics (ESS and PSRF)
summary(effectiveSize(mpost_PA$Beta))
summary(effectiveSize(mpost_PA$Gamma))
round(mean(gelman.diag(mpost_PA$Beta, multivariate = FALSE)$psrf, na.rm = TRUE), 2)
round(sd(gelman.diag(mpost_PA$Beta, multivariate = FALSE)$psrf, na.rm = TRUE), 2)
round(mean(gelman.diag(mpost_PA$Gamma, multivariate = FALSE)$psrf, na.rm = TRUE), 2)
round(sd(gelman.diag(mpost_PA$Gamma, multivariate = FALSE)$psrf, na.rm = TRUE), 2)

summary(effectiveSize(mpost_ABU$Beta))
summary(effectiveSize(mpost_ABU$Gamma))
round(mean(gelman.diag(mpost_ABU$Beta, multivariate = FALSE)$psrf, na.rm = TRUE), 2)
round(sd(gelman.diag(mpost_ABU$Beta, multivariate = FALSE)$psrf, na.rm = TRUE), 2)
round(mean(gelman.diag(mpost_ABU$Gamma, multivariate = FALSE)$psrf, na.rm = TRUE), 2)
round(sd(gelman.diag(mpost_ABU$Gamma, multivariate = FALSE)$psrf, na.rm = TRUE), 2)

# ------------------------------------------------------------------------------
# Figure S4: Beta diagnostics (ESS + PSRF)
# ------------------------------------------------------------------------------

png(filename = file.path("Figures", "Figure S4.png"), width = 10, height = 10, units = "in", res = 400)
par(mfrow = c(2, 2))
hist(effectiveSize(mpost_PA$Beta), main = "ess(beta) - PA model", xlab = "Effective sample size")
vioplot(gelman.diag(mpost_PA$Beta, multivariate = FALSE)$psrf, main = "psrf(beta) - PA model")
hist(effectiveSize(mpost_ABU$Beta), main = "ess(beta) - ABU model", xlab = "Effective sample size")
vioplot(gelman.diag(mpost_ABU$Beta, multivariate = FALSE)$psrf, main = "psrf(beta) - ABU model")
dev.off()

# ------------------------------------------------------------------------------
# Figure S5: Gamma diagnostics (ESS + PSRF)
# ------------------------------------------------------------------------------

png(filename = file.path("Figures", "Figure S5.png"), width = 10, height = 10, units = "in", res = 400)
par(mfrow = c(2, 2))
hist(effectiveSize(mpost_PA$Gamma), main = "ess(gamma) - PA model", xlab = "Effective sample size")
vioplot(gelman.diag(mpost_PA$Gamma, multivariate = FALSE)$psrf, main = "psrf(gamma) - PA model")
hist(effectiveSize(mpost_ABU$Gamma), main = "ess(gamma) - ABU model", xlab = "Effective sample size")
vioplot(gelman.diag(mpost_ABU$Gamma, multivariate = FALSE)$psrf, main = "psrf(gamma) - ABU model")
dev.off()

# ------------------------------------------------------------------------------
# Model fit evaluation
# ------------------------------------------------------------------------------

pred.mpa <- computePredictedValues(mPA, expected = TRUE)
MF.mpa <- evaluateModelFit(mPA, predY = pred.mpa)

pred.abu <- computePredictedValues(mABU, expected = TRUE)
MF.abu <- evaluateModelFit(mABU, predY = pred.abu)

df_mpa <- data.frame(Species = mPA$spNames, onto = substr(mPA$spNames, start = 9, stop = nchar(mPA$spNames)), RMSE = MF.mpa$RMSE, AUC = MF.mpa$AUC, TjurR2 = MF.mpa$TjurR2, Prevalence = colSums(mPA$Y)/mPA$ny, Model = "Presence-Absence")
df_mpa <- melt(df_mpa)
df_abu <- data.frame(Species = mABU$spNames, onto = substr(mABU$spNames, start = 9, stop = nchar(mABU$spNames)), RMSE = MF.abu$RMSE, R2 = MF.abu$R2, Prevalence = colSums(((mABU$Y>0)*1)/mABU$ny, na.rm = TRUE), Model = "Conditional Abundance")
df_abu <- melt(df_abu)

tab <- rbind(df_mpa, df_abu)

# ------------------------------------------------------------------------------
# Summaries and tests
# ------------------------------------------------------------------------------

tab_summary <- tab %>%
  dplyr::group_by(onto, Model, variable) %>%
  dplyr::summarise(
    mean_value = mean(value, na.rm = TRUE),
    sd_value   = sd(value, na.rm = TRUE),
    n          = dplyr::n(),
    .groups    = "drop"
  )

tab_summary_no_onto <- tab %>%
  dplyr::group_by(Model, variable) %>%
  dplyr::summarise(
    mean_value = mean(value, na.rm = TRUE),
    sd_value = sd(value, na.rm = TRUE),
    n = dplyr::n()
  ) %>%
  ungroup()

tab_r2 <- tab %>% filter(variable == "R2")
wilcox_result <- wilcox.test(value ~ onto, data = tab_r2)

tab_tjur <- tab %>% filter(variable == "TjurR2")
wilcox_result <- wilcox.test(value ~ onto, data = tab_tjur)

tab_auc <- tab %>% filter(variable == "AUC")
wilcox_result <- wilcox.test(value ~ onto, data = tab_auc)

# ------------------------------------------------------------------------------
# Figure S6: Model fit indicators
# ------------------------------------------------------------------------------

# indic <- c("AUC", "TjurR2", "R2")
# tab.fit <- tab[tab$variable %in% indic, ]

# gg.expl_power <- ggplot(tab.fit, aes(x = variable, y = value, fill = onto)) + 
#  geom_boxplot(alpha = .7, outlier.shape = NA) +
#  geom_jitter(position = position_jitterdodge(jitter.width = .15), alpha = .2, shape = 21) +
#  scale_fill_manual(values = c("#E69F00", "#56B4E9"), name = "Life stage", labels = c("Adult", "Juvenile")) +
#  labs(y = "Value", x = NULL) +
#  theme_bw()

# Save Figure S6
#ggplot2::ggsave(gg.expl_power, filename = file.path("Figures", "Figure S6.png"),
#                width = 30, height = 15, units = "cm", dpi = 400)