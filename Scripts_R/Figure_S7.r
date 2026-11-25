###############################################################################
# Title: Trait Contributions to Variance Partitioning (Figure S7)
###############################################################################

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------

source("Scripts_R/0-Load libraries.r")

# -----------------------------------------------------------------------------
# Load fitted models
# -----------------------------------------------------------------------------

mpa_model <- get(load("Models/model_mpa_thin_4000_samples_250_chains_4.RData"))
abu_model <- get(load("Models/model_abu_thin_4000_samples_250_chains_4.RData"))

# -----------------------------------------------------------------------------
# Variance partitioning setup
# -----------------------------------------------------------------------------

group_names <- c(
  "Depth",
  "Temperature bot",
  "Temperature surf",
  "Salinity",
  "ChlorophyllA",
  "Fishing",
  "Gravity"
)

group_index <- c(
  1, 1, 1,  # Depth
  2, 2,     # SBT
  3, 3,     # SST
  4, 4,     # SSS
  5, 5,     # Chla
  6,        # Fishing
  7         # Gravity
)

vp_mpa <- computeVariancePartitioning(mpa_model, group = group_index, groupnames = group_names)
vp_abu <- computeVariancePartitioning(abu_model, group = group_index, groupnames = group_names)

round(vp_mpa$R2T$Y, 2)
round(vp_mpa$R2T$Beta, 2)
round(vp_abu$R2T$Y, 2)
round(vp_abu$R2T$Beta, 2)

# -----------------------------------------------------------------------------
# Trait contributions (% of variance explained by traits per covariate)
# -----------------------------------------------------------------------------

covariates <- c(
  "Intercept", "Depth 1", "Depth 2",
  "SBT 1", "SBT 2",
  "SST 1", "SST 2",
  "SSS 1", "SSS 2",
  "Chla 1", "Chla 2",
  "Fishing", "Gravity"
)

traits_beta_mpa <- data.frame(
  Model     = rep("Presence-Absence", length(covariates)),
  Variable  = covariates,
  value_pct = vp_mpa$R2T$Beta * 100
)

traits_beta_abu <- data.frame(
  Model     = rep("Conditional Abundance", length(covariates)),
  Variable  = covariates,
  value_pct = vp_abu$R2T$Beta * 100
)

# Drop intercept and order by contribution size
traits_beta_mpa <- traits_beta_mpa %>%
  filter(Variable != "Intercept") %>%
  mutate(Variable = factor(Variable)) %>%
  arrange(desc(value_pct)) %>%
  mutate(Variable = factor(Variable, levels = Variable))

traits_beta_abu <- traits_beta_abu %>%
  filter(Variable != "Intercept") %>%
  mutate(Variable = factor(Variable)) %>%
  arrange(desc(value_pct)) %>%
  mutate(Variable = factor(Variable, levels = Variable))

# -----------------------------------------------------------------------------
# Plots
# -----------------------------------------------------------------------------

plot_traits_mpa <- ggplot(traits_beta_mpa, aes(x = Variable, y = value_pct, fill = value_pct)) +
  geom_col(position = position_dodge()) +
  scale_fill_viridis_c() +
  labs(y = "Proportion %", x = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12)) +
  guides(fill = "none") +
  ggtitle("Presence-Absence")

plot_traits_abu <- ggplot(traits_beta_abu, aes(x = Variable, y = value_pct, fill = value_pct)) +
  geom_col(position = position_dodge()) +
  scale_fill_viridis_c() +
  labs(y = NULL, x = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12)) +
  guides(fill = "none") +
  ggtitle("Conditional Abundance")

# Combine plots
combined_traits <- plot_traits_mpa + plot_traits_abu +
  plot_layout(guides = "collect", widths = c(1, 1))

# ggplot2::ggsave(combined_traits, filename = "Figures/Figure_S7.png", width = 30, height = 30, units = "cm", dpi = 400)

# -----------------------------------------------------------------------------
# Gamma parameter heatmaps
# -----------------------------------------------------------------------------

# Define support level
supportLevel <- 0.95

# Gamma heatmap for presence/absence model
postGamma <-  getPostEstimate(mpa_model, parName = "Gamma")

temp.gamma <- as.data.frame(postGamma$mean)
names(temp.gamma) <- mpa_model$trNames
temp.gamma$param <- mpa_model$covNames
temp.gamma <- temp.gamma %>% tidyr::pivot_longer(.,-param, names_to = 'traits', values_to ='value')

sup.gamma <- as.data.frame(postGamma$support) 
names(sup.gamma) <- mpa_model$trNames
sup.gamma$param <- mpa_model$covNames
sup.gamma <- sup.gamma %>% tidyr::pivot_longer(.,-param, names_to = 'traits', values_to ='value')

temp.gamma$support <- sup.gamma$value

temp.gamma$pos.neg <- 0
temp.gamma$pos.neg[temp.gamma$support > supportLevel] <- 1
temp.gamma$pos.neg[temp.gamma$support < (1-supportLevel)] <- -1

temp.gamma$pos.neg <- factor(temp.gamma$pos.neg, levels = c(1,0,-1)) 

# Rename variables
temp.gamma$Environment <- factor(temp.gamma$param)
levels(temp.gamma$Environment)
temp.gamma$Environment <- dplyr::recode(temp.gamma$Environment, 
                                        "Fishing" = "Fishing",
                                        "Gravity" = "Gravity",
                                        "poly(Chla, degree = 2, raw = TRUE)1" = "Chla1",
                                        "poly(Chla, degree = 2, raw = TRUE)2" = "Chla2",
                                        "poly(Depth, degree = 2, raw = TRUE)1" = "Depth1",
                                        "poly(Depth, degree = 2, raw = TRUE)2" = "Depth2",
                                        "poly(SBT, degree = 2, raw = TRUE)1" = "SBT1",
                                        "poly(SBT, degree = 2, raw = TRUE)2" = "SBT2",
                                        "poly(SSS, degree = 2, raw = TRUE)1" = "SSS1",
                                        "poly(SSS, degree = 2, raw = TRUE)2" = "SSS2",
                                        "poly(SST, degree = 2, raw = TRUE)1" = "SST1",
                                        "poly(SST, degree = 2, raw = TRUE)2" = "SST2")

temp.gamma$Environment <-  factor(temp.gamma$Environment,
                                  levels = c("(Intercept)", "Depth1", "Depth2", "SST1", "SST2",
                                             "SBT1", "SBT2", "SSS1","SSS2", "Chla1", "Chla2", 
                                             "Fishing", "Gravity"))

# Rename traits:
temp.gamma$Traits <- factor(temp.gamma$traits)
temp.gamma$Traits <- dplyr::recode(temp.gamma$Traits, 
                                   "Body_shapeeel-like" = "Eel-like (body shape)",
                                   "Body_shapeelongated" = "Elongated (body shape)",
                                   "Body_shapeflat" = "Flat  (body shape)",
                                   "Body_shapefusiform" = "Fusiform  (body shape)",
                                   "Body_shapeshort and/or deep" = "Short/deep (body shape)",
                                   "Trophic_guildgeneralist" =  "Generalist (Trophic guild)",
                                   "Trophic_guildpiscivorous" = "Piscivorous (Trophic guild)",
                                   "Trophic_guildplanktivorous" = "Planktivorous (Trophic guild)",
                                   "Water_column_positionbathypelagic" = "Bathypelagic (habitat)",
                                   "Water_column_positionbenthopelagic" = "Benthopelagic (habitat)",
                                   "Water_column_positiondemersal" = "Demersal (habitat)",
                                   "Water_column_positionpelagic" = "Pelagic (habitat)",
                                   "Water_column_positionreef-associated" = "Reef-associated (habitat)",
                                   "Growth_coefficient" = "Growth (K)",
                                   "Max_size" = "Maximum length",
                                   "Min_size" = "Minimum length",
                                   "Trophic_level" = "Trophic level")

temp.gamma$Traits <-  factor(temp.gamma$Traits,
                             levels = rev(c("Minimum length", "Maximum length", "Growth (K)", "Trophic level",
                                            "Generalist (Trophic guild)", "Piscivorous (Trophic guild)", "Planktivorous (Trophic guild)",
                                            "Eel-like (body shape)", "Elongated (body shape)", "Flat  (body shape)", "Fusiform  (body shape)", "Short/deep (body shape)", 
                                            "Bathypelagic (habitat)", "Benthopelagic (habitat)", "Demersal (habitat)", "Pelagic (habitat)", "Reef-associated (habitat)",
                                            "(Intercept)")))

gamma.heatmap.PA <- ggplot(temp.gamma, aes(x = Environment, y = Traits)) +
  geom_tile(aes(fill=pos.neg, color = pos.neg),color='grey') +
  scale_fill_manual(values=c('#F21A00','white','#3B9AB2'),labels=c('+','','-'), name='', guide=guide_legend(keyheight =4, keywidth = 1)) + xlab('') + ylab('')+
  theme_bw()+
  theme(axis.text.y =element_text(size=10),
        strip.text.y = element_blank(),
        axis.text.x  = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size = 12, face = "bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.grid = element_blank(),
        panel.background = element_blank(), strip.background = element_blank(),
        plot.margin = margin(t = 5,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0), # Left margin
        panel.spacing.y = unit(0.1, "lines"))+
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  ggtitle("Presence-Absence")


# Gamma heatmap for abundance model
postGamma <-  getPostEstimate(abu_model, parName = "Gamma")

temp.gamma <- as.data.frame(postGamma$mean) # "temp" stands for temporary
names(temp.gamma) <- abu_model$trNames
temp.gamma$param <- abu_model$covNames
temp.gamma <- temp.gamma %>% tidyr::pivot_longer(.,-param, names_to = 'traits', values_to ='value')

sup.gamma <- as.data.frame(postGamma$support) 
names(sup.gamma) <- abu_model$trNames
sup.gamma$param <- abu_model$covNames
sup.gamma <- sup.gamma %>% tidyr::pivot_longer(.,-param, names_to = 'traits', values_to ='value')

temp.gamma$support <- sup.gamma$value

temp.gamma$pos.neg <- 0
temp.gamma$pos.neg[temp.gamma$support > supportLevel] <- 1
temp.gamma$pos.neg[temp.gamma$support < (1-supportLevel)] <- -1
temp.gamma$pos.neg <- factor(temp.gamma$pos.neg, levels = c(1,0,-1)) 

# Rename variables
temp.gamma$Environment <- factor(temp.gamma$param)
levels(temp.gamma$Environment)
temp.gamma$Environment <- dplyr::recode(temp.gamma$Environment, 
                                        "Fishing" = "Fishing",
                                        "Gravity" = "Gravity",
                                        "poly(Chla, degree = 2, raw = TRUE)1" = "Chla1",
                                        "poly(Chla, degree = 2, raw = TRUE)2" = "Chla2",
                                        "poly(Depth, degree = 2, raw = TRUE)1" = "Depth1",
                                        "poly(Depth, degree = 2, raw = TRUE)2" = "Depth2",
                                        "poly(SBT, degree = 2, raw = TRUE)1" = "SBT1",
                                        "poly(SBT, degree = 2, raw = TRUE)2" = "SBT2",
                                        "poly(SSS, degree = 2, raw = TRUE)1" = "SSS1",
                                        "poly(SSS, degree = 2, raw = TRUE)2" = "SSS2",
                                        "poly(SST, degree = 2, raw = TRUE)1" = "SST1",
                                        "poly(SST, degree = 2, raw = TRUE)2" = "SST2")

temp.gamma$Environment <-  factor(temp.gamma$Environment,
                                  levels = c("(Intercept)", "Depth1", "Depth2", "SST1", "SST2",
                                             "SBT1", "SBT2", "SSS1","SSS2", "Chla1", "Chla2", 
                                             "Fishing", "Gravity"))

# Rename traits:
temp.gamma$Traits <- factor(temp.gamma$traits)
temp.gamma$Traits <- dplyr::recode(temp.gamma$Traits, 
                                   "Body_shapeeel-like" = "Eel-like (body shape)",
                                   "Body_shapeelongated" = "Elongated (body shape)",
                                   "Body_shapeflat" = "Flat  (body shape)",
                                   "Body_shapefusiform" = "Fusiform  (body shape)",
                                   "Body_shapeshort and/or deep" = "Short/deep (body shape)",
                                   "Trophic_guildgeneralist" =  "Generalist (Trophic guild)",
                                   "Trophic_guildpiscivorous" = "Piscivorous (Trophic guild)",
                                   "Trophic_guildplanktivorous" = "Planktivorous (Trophic guild)",
                                   "Water_column_positionbathypelagic" = "Bathypelagic (habitat)",
                                   "Water_column_positionbenthopelagic" = "Benthopelagic (habitat)",
                                   "Water_column_positiondemersal" = "Demersal (habitat)",
                                   "Water_column_positionpelagic" = "Pelagic (habitat)",
                                   "Water_column_positionreef-associated" = "Reef-associated (habitat)",
                                   "Growth_coefficient" = "Growth (K)",
                                   "Max_size" = "Maximum length",
                                   "Min_size" = "Minimum length",
                                   "Trophic_level" = "Trophic level")

temp.gamma$Traits <-  factor(temp.gamma$Traits,
                             levels = rev(c("Minimum length", "Maximum length", "Growth (K)", "Trophic level",
                                            "Generalist (Trophic guild)", "Piscivorous (Trophic guild)", "Planktivorous (Trophic guild)",
                                            "Eel-like (body shape)", "Elongated (body shape)", "Flat  (body shape)", "Fusiform  (body shape)", "Short/deep (body shape)", 
                                            "Bathypelagic (habitat)", "Benthopelagic (habitat)", "Demersal (habitat)", "Pelagic (habitat)", "Reef-associated (habitat)",
                                            "(Intercept)")))

gamma.heatmap.AB <- ggplot(temp.gamma, aes(x = Environment, y = Traits)) +
  geom_tile(aes(fill=pos.neg, color = pos.neg),color='grey') +
  scale_fill_manual(values=c('#F21A00','white','#3B9AB2'),labels=c('+','','-'), name='',guide=guide_legend(keyheight =4, keywidth = 1)) + xlab('') + ylab('')+
  theme_bw()+
  theme(axis.text.y =element_text(size=6),
        strip.text = element_blank(),
        axis.text.x  = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size = 12, face = "bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.grid = element_blank(),
        panel.background = element_blank(), strip.background = element_blank(),
        plot.margin = margin(t = 5,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0)) + # Left margin
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme(axis.text.y = element_blank(),  # Hide Y-axis labels for the right plot
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.spacing.y = unit(0.1, "lines"))+
  ggtitle("Conditional Abundance")

# -----------------------------------------------------------------------------
# Combine plot for P/A and abundance model
# -----------------------------------------------------------------------------
combined_gamma_plot <- gamma.heatmap.PA + gamma.heatmap.AB + plot_layout(guides = "collect", widths = c(1, 1))

# ggplot2::ggsave(combined_gamma_plot, filename = "Figures/Figure_S7_2.png", width = 30, height = 30, units = "cm", dpi = 400)
