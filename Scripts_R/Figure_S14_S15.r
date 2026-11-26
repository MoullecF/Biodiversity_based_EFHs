###############################################################################
# Title: Covariate Effects and Life-Stage Responses (Figures S14-S15)
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
# Helpers
# -----------------------------------------------------------------------------

support_level <- 0.95

env_recode <- c(
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
  "poly(SST, degree = 2, raw = TRUE)2" = "SST2"
)

env_levels <- c(
  "(Intercept)", "Depth1", "Depth2", "SST1", "SST2",
  "SBT1", "SBT2", "SSS1", "SSS2", "Chla1", "Chla2",
  "Fishing", "Gravity"
)

build_beta_tables <- function(model) {
  post_beta <- getPostEstimate(model, "Beta")

  beta_mean <- as.data.frame(post_beta$mean) %>%
    mutate(param = model$covNames) %>%
    pivot_longer(-param, names_to = "species", values_to = "value")

  beta_support <- as.data.frame(post_beta$support) %>%
    mutate(param = model$covNames) %>%
    pivot_longer(-param, names_to = "species", values_to = "support")

  beta_mean %>%
    left_join(beta_support, by = c("param", "species")) %>%
    mutate(
      pos_neg = case_when(
        support > support_level ~ 1,
        support < (1 - support_level) ~ -1,
        TRUE ~ 0
      ),
      pos_neg = factor(pos_neg, levels = c(1, 0, -1)),
      Environment = factor(recode(param, !!!env_recode), levels = env_levels),
      Onto = substr(species, 9, nchar(species))
    )
}

ordered_species <- function(model) {
  colnames(model$Y)[order(grepl("adult", colnames(model$Y)))]
}

build_heatmap_plot <- function(df, title, show_y_axis = TRUE) {
  ggplot(df, aes(x = Environment, y = Species_sort)) +
    geom_tile(aes(fill = pos_neg, color = pos_neg), color = "grey") +
    scale_fill_manual(
      values = c("#F21A00", "white", "#3B9AB2"),
      labels = c("+", "", "-"),
      name = "",
      guide = guide_legend(keyheight = 4, keywidth = 1)
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = "", y = "", title = title) +
    theme_bw() +
    theme(
      strip.text.y = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.text = element_text(size = 12, face = "bold"),
      panel.grid = element_blank(),
      panel.background = element_blank(),
      strip.background = element_blank(),
      plot.margin = margin(t = 5, r = 0, b = 0, l = 0),
      panel.spacing.y = unit(0.1, "lines"),
      axis.text.y = if (show_y_axis) element_text(size = 6) else element_blank(),
      axis.title.y = if (show_y_axis) element_text() else element_blank(),
      axis.ticks.y = if (show_y_axis) element_line() else element_blank()
    ) +
    facet_grid(Onto ~ ., scales = "free_y", space = "free_y")
}

build_prop_plot <- function(df, title, show_y_axis = TRUE) {
  df %>%
    count(param, pos_neg, Onto) %>%
    group_by(param) %>%
    mutate(percent = round(100 * n / sum(n), 1)) %>%
    ungroup() %>%
    mutate(
      Onto = factor(Onto, levels = c("juvenile", "adult")),
      pos_neg = factor(pos_neg, levels = c(0, 1, -1)),
      effect_onto = factor(
        paste0(pos_neg, "_", Onto),
        levels = c("0_adult", "0_juvenile", "1_adult", "1_juvenile", "-1_adult", "-1_juvenile")
      ),
      param = recode(param, !!!env_recode),
      param = factor(param, levels = env_levels)
    ) %>%
    filter(param != "(Intercept)") %>%
    ggplot(aes(x = param, y = percent, fill = effect_onto)) +
    geom_col(position = "stack") +
    scale_fill_manual(
      values = c(
        "-1_juvenile" = "#6BAED6",
        "-1_adult" = "#08519C",
        "1_juvenile" = "#F8766D",
        "1_adult" = "#C00000",
        "0_juvenile" = "grey90",
        "0_adult" = "grey85"
      ),
      labels = c(
        "-1_juvenile" = "- (Juveniles)",
        "-1_adult" = "- (Adults)",
        "1_juvenile" = "+ (Juveniles)",
        "1_adult" = "+ (Adults)",
        "0_juvenile" = "NS (Juveniles)",
        "0_adult" = "NS (Adults)"
      ),
      name = "Responses",
      guide = guide_legend(keyheight = 1.5, keywidth = 1)
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = "", y = if (show_y_axis) "Proportion of species-life stages (%)" else NULL, title = title) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = if (show_y_axis) element_text() else element_blank(),
      axis.ticks.y = if (show_y_axis) element_line() else element_blank(),
      plot.margin = margin(t = 5, r = 0, b = 0, l = 0),
      panel.spacing.y = unit(0.1, "lines")
    )
}

# -----------------------------------------------------------------------------
# Presence–absence model
# -----------------------------------------------------------------------------

beta_pa <- build_beta_tables(mpa_model) %>%
  mutate(Species_sort = factor(species, levels = rev(ordered_species(mpa_model))))

heatmap_pa <- build_heatmap_plot(beta_pa, "Presence-Absence", show_y_axis = TRUE)

prop_pa <- build_prop_plot(beta_pa, "Presence–absence", show_y_axis = TRUE)

# -----------------------------------------------------------------------------
# Conditional abundance model
# -----------------------------------------------------------------------------

beta_ab <- build_beta_tables(abu_model) %>%
  mutate(Species_sort = factor(species, levels = rev(ordered_species(abu_model))))

heatmap_ab <- build_heatmap_plot(beta_ab, "Conditional Abundance", show_y_axis = FALSE)

prop_ab <- build_prop_plot(beta_ab, "Abundance conditional on presence", show_y_axis = FALSE)

# -----------------------------------------------------------------------------
# Combine plots and save
# -----------------------------------------------------------------------------

combined_beta_plot <- heatmap_pa + heatmap_ab + plot_layout(guides = "collect", widths = c(1, 1))

# ggplot2::ggsave(combined_beta_plot, filename = "Figures/Figure_S14.png", width = 30, height = 30, units = "cm", dpi = 400)

combined_prop <- prop_pa + prop_ab + plot_layout(guides = "collect", widths = c(1, 1))

# ggplot2::ggsave(combined_prop, filename = "Figures/Figure_S15.png", width = 30, height = 30, units = "cm", dpi = 400)

# Proportion of species with positive, negative, or non-significant responses to a given covariate
summarise_response <- function(df, model_label, variable) {
  df %>%
    filter(Environment == variable) %>%
    count(pos_neg) %>%
    mutate(
      response = recode(
        as.character(pos_neg),
        "1" = "Positive",
        "0" = "Non-significant",
        "-1" = "Negative"
      ),
      percent = round(100 * n / sum(n), 1),
      model = model_label
    ) %>%
    dplyr::select(model, response, n, percent)
}

response_summary <- bind_rows(summarise_response(beta_pa, "Presence-absence", "Gravity"),
  summarise_response(beta_ab, "Abundance", "Gravity"))

