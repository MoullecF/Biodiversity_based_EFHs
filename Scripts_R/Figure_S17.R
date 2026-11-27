###############################################################################
# Title: Figure S17 - Species Associations Heatmap
###############################################################################

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------

source("Scripts_R/0-Load libraries.r")

support_level <- 0.95
eps <- 1e-9

# -----------------------------------------------------------------------------
# Load fitted model
# -----------------------------------------------------------------------------

mpa_model <- get(load("Models/model_mpa_thin_4000_samples_250_chains_4.RData"))

# -----------------------------------------------------------------------------
# Residual associations at the grid-cell level
# -----------------------------------------------------------------------------

omega_cor <- computeAssociations(mpa_model)
names(omega_cor) <- names(mpa_model$ranLevels)
omega <- omega_cor$grid.cell

# Posterior mean correlation where supported (>= 0.95 or <= 0.05), else 0
to_plot <- ((omega$support >= support_level - eps) +
  (omega$support <= (1 - support_level) + eps) > 0) * omega$mean
plot_order <- corrMatOrder(omega$mean, order = "AOE")
corr <- as.matrix(to_plot[plot_order, plot_order])
corr[corr > 0] <- 1
corr[corr < 0] <- -1

n_comb <- nrow(corr)
stopifnot(n_comb == ncol(mpa_model$Y)) # 116 species-life stage combinations
corr[lower.tri(corr, diag = TRUE)] <- NA
melted_cormat <- reshape2::melt(corr, na.rm = TRUE)
melted_cormat$value <- factor(melted_cormat$value, levels = c("1", "0", "-1"))

# -----------------------------------------------------------------------------
# Heatmap
# -----------------------------------------------------------------------------

omega_plot <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value)) +
  scale_x_discrete(limits = rev(levels(melted_cormat$Var2))) +
  geom_tile(color = "grey") +
  labs(x = "", y = "") +
  scale_fill_manual(values = c("#F21A00", "white", "#3B9AB2"), labels = c("+", "", "-"), name = "", guide = guide_legend(keyheight = 4, keywidth = 1)) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 6),
    legend.text = element_text(size = 10),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.25),
    panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
    panel.background = element_blank(), strip.background = element_blank(),
    legend.position.inside = c(0.96, 0.6),
    plot.margin = unit(c(11, 11, 5.5, 5.5), "pt")  )

# ggplot2::ggsave(omega_plot, filename = "Figures/Figure_S17.png", width = 30, height = 30, units = "cm", dpi = 400)

# -----------------------------------------------------------------------------
# Summary counts: pairwise associations between life stages
# -----------------------------------------------------------------------------

melted_cormat <- melted_cormat %>%
  mutate(
    stage1 = ifelse(grepl("juvenile", Var1), "juvenile", "adult"),
    stage2 = ifelse(grepl("juvenile", Var2), "juvenile", "adult"),
    pair_type = case_when(
      stage1 == "juvenile" & stage2 == "juvenile" ~ "juvenile only",
      stage1 == "adult" & stage2 == "adult" ~ "adult only",
      TRUE ~ "juvenile-adult"),
    sign = case_when(value == "1" ~ "positive", value == "-1" ~ "negative", TRUE ~ "none"))

n_possible <- nrow(melted_cormat)
n_pos <- sum(melted_cormat$sign == "positive")
n_neg <- sum(melted_cormat$sign == "negative")
n_sig <- n_pos + n_neg

cat(sprintf("%d significant associations out of %d possible pairs (%.1f%%): %d positive (%.1f%%), %d negative (%.1f%%)\n", n_sig, n_possible, 100 * n_sig / n_possible, n_pos, 100 * n_pos / n_sig, n_neg, 100 * n_neg / n_sig))

summary_tbl <- melted_cormat %>%
  filter(sign != "none") %>%
  dplyr::count(sign, pair_type) %>%
  group_by(sign) %>%
  mutate(percent = round(100 * n / sum(n), 1)) %>%
  ungroup()

print(summary_tbl)
