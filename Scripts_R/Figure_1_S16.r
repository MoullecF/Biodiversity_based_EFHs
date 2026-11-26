###############################################################################
# Title: Variance Partitioning for HMSC Models (Figure 1)
###############################################################################

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------

rm(list = ls())
source("Scripts_R/0-Load libraries.R")
source("Scripts_R/ggplot_theme.R")

# Load fitted models
Mpa <- get(load("Models/model_mpa_thin_4000_samples_250_chains_4.RData"))
Mabu <- get(load("Models/model_abu_thin_4000_samples_250_chains_4.RData"))

# Variable groups used in variance partitioning
groupnames <- c("Depth", "Temperature", "Salinity", "ChlorophyllA", "Human pressure")
group <- c(
  1, 1, 1,
  2, 2, 2, 2,
  3, 3,
  4, 4,
  5, 5
)

variable_labels <- c(
  "Depth" = "Depth",
  "Temperature" = "Temperature",
  "Salinity" = "Salinity",
  "ChlorophyllA" = "Chlorophyll a",
  "Human.pressure" = "Human pressures",
  "Random..year" = "Year",
  "Random..grid.cell" = "Grid cell"
)

cols <- c("#E69F00", "#56B4E9")

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
compute_vp <- function(model, group, groupnames, variable_labels) {
  # Compute variance partitioning and return a long table with labels applied.
  vp <- computeVariancePartitioning(model, group = group, groupnames = groupnames)
  vp <- data.frame(t(vp$vals))
  vp$Species <- rownames(vp)
  vp <- reshape2::melt(vp)
  vp$Onto <- substr(vp$Species, start = 9, stop = nchar(vp$Species))
  vp <- vp %>% mutate(variable = recode(variable, !!!variable_labels))
  vp
}

summarize_vp <- function(vp) {
  # Summaries used for plot annotations and optional inspection.
  list(
    ordered = vp %>% arrange(variable, desc(value)),
    by_onto = vp %>%
      group_by(variable, Onto) %>%
      summarise(Mean = round(mean(value), 4), SD = round(sd(value), 4)),
    overall = vp %>%
      group_by(variable) %>%
      summarise(Mean = round(mean(value), 2), SD = round(sd(value), 4))
  )
}

plot_vp <- function(vp, vp_summary, title, y_label) {
  ggplot(vp, aes(x = variable, y = value, fill = Onto)) +
    geom_flat_violin(
      position = position_nudge(x = 0.2),
      scale = "width",
      trim = TRUE,
      alpha = 0.6,
      show.legend = TRUE
    ) +
    geom_point(
      aes(color = Onto),
      position = position_jitter(w = 0.1),
      size = 1.25,
      alpha = 0.35,
      show.legend = FALSE
    ) +
    stat_summary(
      fun = mean,
      geom = "point",
      aes(fill = Onto),
      color = "black",
      shape = 23,
      size = 3,
      show.legend = FALSE
    ) +
    geom_vline(xintercept = 5.7, linetype = "dashed") +
    geom_text(
      data = vp_summary$by_onto[vp_summary$by_onto$Onto == "adult", ],
      aes(label = paste0(round(Mean, 2), "±", round(SD, 2)), y = -0.05, x = c(1:7)),
      color = "#E69F00",
      size = 2.75,
      fontface = "bold"
    ) +
    geom_text(
      data = vp_summary$by_onto[vp_summary$by_onto$Onto == "juvenile", ],
      aes(label = paste0(round(Mean, 2), "±", round(SD, 2)), y = -0.08, x = c(1:7)),
      color = "#56B4E9",
      size = 2.75,
      fontface = "bold"
    ) +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols, name = "Life stage", labels = c("Adult", "Juvenile")) +
    labs(y = y_label, x = "") +
    guides(color = "none") +
    theme_bw() +
    scale_x_discrete(expand = expansion(mult = c(0.06, 0.12))) +
    ggtitle(title) +
    theme(
      axis.title.y = element_text(size = 12),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      plot.title = element_text(size = 12)
    )
}

# -----------------------------------------------------------------------------
# Variance partitioning
# -----------------------------------------------------------------------------
VP.Mpa <- compute_vp(Mpa, group, groupnames, variable_labels)
VP.Mabu <- compute_vp(Mabu, group, groupnames, variable_labels)

VP.all <- rbind(
  cbind(VP.Mpa,  Model = "Presence-Absence"),
  cbind(VP.Mabu, Model = "Abundance")
)

# write.csv(VP.all, "Outputs/VP_all.csv", row.names = FALSE)

vpSummary.Mpa <- summarize_vp(VP.Mpa)
vpSummary.Mabu <- summarize_vp(VP.Mabu)

# -----------------------------------------------------------------------------
# Plots
# -----------------------------------------------------------------------------
gg.VP.Mpa <- plot_vp(
  VP.Mpa,
  vpSummary.Mpa,
  title = "Presence-Absence (PA)",
  y_label = "Proportion of explained variation\n"
)

gg.VP.Mabu <- plot_vp(
  VP.Mabu,
  vpSummary.Mabu,
  title = "Abundance conditional on presence (ABU)",
  y_label = ""
) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
  )

# -----------------------------------------------------------------------------
# Combine plots
# -----------------------------------------------------------------------------
gg.VP <- gg.VP.Mpa + gg.VP.Mabu + plot_layout(guides = "collect", widths = c(1, 1))

# -----------------------------------------------------------------------------
# Save plot
# -----------------------------------------------------------------------------
ggplot2::ggsave(gg.VP, filename = "Figures/Figure1.png", width = 30, height = 15, units = "cm", dpi = 400)

# -----------------------------------------------------------------------------
# Compare juvenile and adult β estimates for species with both life stages
# -----------------------------------------------------------------------------
# Figure S16 in supplementary materials
library(tibble)

### Extract β posterior estimates and keep species with both stages
# Extract posterior estimates
postBeta.mpa <- getPostEstimate(Mpa, parName = "Beta")
postBeta.mabu <- getPostEstimate(Mabu, parName = "Beta")

# Posterior mean
beta.mpa <- postBeta.mpa$mean
beta.mabu <- postBeta.mabu$mean

# Add names
rownames(beta.mpa) <- colnames(Mpa$X)
colnames(beta.mpa) <- colnames(postBeta.mpa$mean)
rownames(beta.mabu) <- colnames(Mabu$X)
colnames(beta.mabu) <- colnames(postBeta.mabu$mean)

# Keep only significant beta estimates
support_level <- 0.95
beta_sig_mpa <- beta.mpa
beta_sig_mabu <- beta.mabu
beta_sig_mpa[postBeta.mpa$support < support_level & postBeta.mpa$supportNeg < support_level] <- NA
beta_sig_mabu[postBeta.mabu$support < support_level & postBeta.mabu$supportNeg < support_level] <- NA

# Species and life stage extraction
species <- substr(colnames(beta_sig_mpa), 1, 7)
stage <- sub(".*_", "", colnames(beta_sig_mpa))

# Species with both juvenile and adult
species_both <- names(which(tapply(stage, species, function(x) all(c("juvenile","adult") %in% x))))

# Keep only species with both stages
beta_both_mpa <- beta_sig_mpa[, species %in% species_both]
beta_both_mabu <- beta_sig_mabu[, species %in% species_both]

### Convert β matrices into long format
beta_long_mpa <- beta_both_mpa %>%
  as.data.frame() %>%
  rownames_to_column("environment") %>%
  pivot_longer(cols = -environment,
    names_to = "species_stage",
    values_to = "beta") %>%
  mutate(species = substr(species_stage,1,7),
    stage = sub(".*_", "", species_stage))

beta_long_mabu <- beta_both_mabu %>%
  as.data.frame() %>%
  rownames_to_column("environment") %>%
  pivot_longer(cols = -environment,
    names_to = "species_stage",
    values_to = "beta") %>%
  mutate(species = substr(species_stage,1,7),
    stage = sub(".*_", "", species_stage))

### Pair juvenile and adult β values
beta_pair_mpa <- beta_long_mpa %>%
  dplyr::select(environment, species, stage, beta) %>%
  pivot_wider(names_from = stage,
    values_from = beta) %>%
  filter(!is.na(juvenile),
    !is.na(adult)) %>%
  mutate(beta_difference = adult - juvenile,
    same_direction = sign(adult)==sign(juvenile))

beta_pair_mabu <- beta_long_mabu %>%
  dplyr::select(environment, species, stage, beta) %>%
  pivot_wider(names_from = stage,
    values_from = beta) %>%
  filter(!is.na(juvenile),
    !is.na(adult)) %>%
  mutate(beta_difference = adult - juvenile,
    same_direction = sign(adult)==sign(juvenile))

# Proportion of species with opposite responses
View(beta_pair_mpa %>%
  group_by(environment) %>%
  summarise(n_species = n(),
    proportion_same_sign = mean(same_direction) * 100,
    proportion_opposite_sign = mean(!same_direction) * 100))

View(beta_pair_mabu %>%
  group_by(environment) %>%
  summarise(n_species = n(),
    proportion_same_sign = mean(same_direction) * 100,
    proportion_opposite_sign = mean(!same_direction) * 100))

### Correlation between juvenile and adult β parameters
beta_correlation_mpa <- beta_pair_mpa %>%
  group_by(environment) %>%
  summarise(n_species = n(),
    Pearson_r = round(cor(juvenile, adult, method="pearson"), 2))

beta_correlation_mabu <- beta_pair_mabu %>%
  group_by(environment) %>%
  summarise(n_species = n(),
    Pearson_r = round(cor(juvenile, adult, method="pearson"), 2))

### Rank predictors by magnitude of ontogenic difference
beta_shift_mpa <- beta_pair_mpa %>%
  group_by(environment) %>%
  summarise(mean_absolute_difference = round(mean(abs(beta_difference)), 2), median_absolute_difference = round(median(abs(beta_difference)), 2)) %>%
  arrange(desc(mean_absolute_difference))

beta_shift_mabu <- beta_pair_mabu %>%
  group_by(environment) %>%
  summarise(mean_absolute_difference = mean(abs(beta_difference)), median_absolute_difference = round(median(abs(beta_difference)), 2)) %>%
  arrange(desc(mean_absolute_difference))

### Plot juvenile vs adult β estimates
# remove intercept from the plot
beta_pair_mpa <- beta_pair_mpa %>% filter(environment != "(Intercept)")
beta_pair_mabu <- beta_pair_mabu %>% filter(environment != "(Intercept)")

beta_plot_mpa <- ggplot(beta_pair_mpa, aes(x=juvenile, y=adult)) +
  geom_point(alpha=0.6) +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  facet_wrap(~environment, scales="free") +
  theme_bw()+
  labs(x="Juvenile β estimate", y="Adult β estimate", tag = "A") +
  theme(plot.tag = element_text(size = 12, face = "bold"))

beta_plot_mabu <- ggplot(beta_pair_mabu, aes(x=juvenile, y=adult)) +
  geom_point(alpha=0.6) +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  facet_wrap(~environment, scales="free") +
  theme_bw()+
  labs(x="Juvenile β estimate", y="Adult β estimate", tag = "B") +
  theme(plot.tag = element_text(size = 12, face = "bold"))

### Save plots
figure_beta_comp <- beta_plot_mpa / beta_plot_mabu
ggplot2::ggsave(figure_beta_comp, filename = "Figures/Figure_S16.png", width = 30, height = 30, units = "cm", dpi = 400)
