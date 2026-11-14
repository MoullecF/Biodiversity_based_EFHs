###############################################################################
# Title: Spatiotemporal Trends in Hurdle-Model Abundance (Figure 2)
###############################################################################

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------

rm(list = ls())
source("./0-Load libraries.R")
source("./ggplot_theme.R")

# -----------------------------------------------------------------------------
# Data preparation
# -----------------------------------------------------------------------------

# Load raster stack of hurdle-model outputs (species, stage, year)
hurdle_stack <- get(load("./Outputs/Spatio_temporal_prediction/r_stack_Hurdle_0.05_0_1000_19992021.Rdata"))

# Load depth mask
depth_mask <- readRDS("./Inputs_HMSC/depth_mask.RDS")

# Apply depth mask to hurdle stack (to remove potential extrapolation beyond depth range)
hurdle_stack <- mask(hurdle_stack, depth_mask)

# Parse life stage and year out of each layer name
layer_meta <- str_match(names(hurdle_stack), "(.*)_(juvenile|adult)_(\\d{4})")
stage <- layer_meta[, 3]
year <- layer_meta[, 4]

# Sum abundance across species within each year x life-stage combination
group_keys <- paste(year, stage, sep = "_")
summed_stack <- stackApply(hurdle_stack, indices = group_keys, fun = sum, na.rm = TRUE)
names(summed_stack) <- unique(group_keys)

cap_outliers_mad <- function(raster_layer) {
  raster_layer[raster_layer == 0] <- NA
  outliers <- outliers_mad(values(raster_layer))
  raster_layer[raster_layer > outliers$limits[2]] <- outliers$limits[2]
  raster_layer
}

stage_year_stack <- stack(lapply(seq_len(nlayers(summed_stack)), function(i) cap_outliers_mad(summed_stack[[i]])))

juvenile_stack <- subset(stage_year_stack, subset = paste0("X", 1999:2021, "_juvenile"))
adult_stack <- subset(stage_year_stack, subset = paste0("X", 1999:2021, "_adult"))

# -----------------------------------------------------------------------------
# Trend estimation (Sen's slope and p-values)
# -----------------------------------------------------------------------------

calc_sen_stat <- function(raster_stack, stat = c("slope", "pvalue")) {
  stat <- match.arg(stat)
  calc(raster_stack, fun = function(x) {
    if (sum(!is.na(x)) < 2) return(NA)
    sen <- sens.slope(na.omit(x))
    if (stat == "slope") sen$estimates else sen$p.value
  })
}

sen_slope_juvenile <- calc_sen_stat(juvenile_stack, stat = "slope")
sen_pvalue_juvenile <- calc_sen_stat(juvenile_stack, stat = "pvalue")
sen_slope_adult <- calc_sen_stat(adult_stack, stat = "slope")
sen_pvalue_adult <- calc_sen_stat(adult_stack, stat = "pvalue")

create_non_significant_mask <- function(pvalue_raster, threshold = 0.05) {
  non_sig_raster <- pvalue_raster > threshold
  non_sig_raster[non_sig_raster == 0] <- NA

  if (all(is.na(values(non_sig_raster)))) {
    return(NULL)
  }

  non_sig_poly <- rasterToPolygons(
    non_sig_raster,
    fun = function(x) !is.na(x),
    dissolve = TRUE
  )

  sf::st_as_sf(non_sig_poly)
}

non_sig_mask_juvenile <- create_non_significant_mask(sen_pvalue_juvenile)
non_sig_mask_adult <- create_non_significant_mask(sen_pvalue_adult)

# Combine slope, p-value and significance flag into one long-format dataframe per life stage
build_slope_df <- function(sen_slope, sen_pvalue) {
  slope_df <- as.data.frame(sen_slope, xy = TRUE, na.rm = FALSE)
  names(slope_df)[3] <- "slope"

  pvalue_df <- as.data.frame(sen_pvalue, xy = TRUE, na.rm = FALSE)
  names(pvalue_df)[3] <- "pvalue"

  slope_df$pvalue <- pvalue_df$pvalue
  slope_df$significant <- slope_df$pvalue < 0.05
  slope_df
}

slope_df_juvenile <- build_slope_df(sen_slope_juvenile, sen_pvalue_juvenile)
slope_df_adult <- build_slope_df(sen_slope_adult, sen_pvalue_adult)

# Proportion of significant and non significant trends per life stage
prop_significant <- function(slope_df) {
  prop_sig <- round(sum(slope_df$significant, na.rm = TRUE) / sum(!is.na(slope_df$significant)), 2)
  prop_non_sig <- round(1 - prop_sig, 2)
  c(significant = prop_sig, non_significant = prop_non_sig)
}

prop_significant_juvenile <- prop_significant(slope_df_juvenile)
prop_significant_adult <- prop_significant(slope_df_adult)

bin_sen_slope <- function(slope_df, breaks, labels) {
  slope_df$value_binned <- cut(
    slope_df$slope,
    breaks = breaks,
    labels = labels,
    include.lowest = TRUE
  )
  slope_df
}

# Bin thresholds differ between stages to reflect their distinct slope ranges
slope_df_juvenile <- bin_sen_slope(slope_df_juvenile,
  breaks = c(-Inf, -16, -10, -7, -4, -1, 0, 1, 2, 4, Inf),
  labels = c("< -16", "[-16;-10[", "[-10;-7[", "[-7;-4[", "[-4;-1[",
    "[-1;0[", "[0;+1[", "[+1;+2[", "[+2;+4[", "> +4"))

slope_df_adult <- bin_sen_slope(slope_df_adult,
  breaks = c(-Inf, -142, -112, -59, -37, -21, -11, -5, 0, 7, Inf),
  labels = c("< -142", "[-142;-112[", "[-112;-59[", "[-59;-37[", "[-37;-21[",
    "[-21;-11[", "[-11;-5[", "[-5;0[", "[0;+7[", "> +7"))

# -----------------------------------------------------------------------------
# Colour palettes
# -----------------------------------------------------------------------------

col_neg_juvenile <- c("#CBDEF0", "#ABCFE5", "#81BADA", "#58A1CE", "#3787C0", "#1B69AF")
col_neg_adult <- c("#E0ECF7", "#CBDEF0", "#ABCFE5", "#81BADA", "#58A1CE", "#3787C0", "#1B69AF", "#084D96")
col_pos_juvenile <- c("#FFFFCC", "#FEE692", "#FEBF5A", "#FD8D3C")
col_pos_adult <- c("#FFFFCC", "#FEE692")

col_vec_juvenile <- c(rev(col_neg_juvenile), col_pos_juvenile)
col_vec_adult <- c(rev(col_neg_adult), col_pos_adult)

# -----------------------------------------------------------------------------
# Trend map plotting
# -----------------------------------------------------------------------------

plot_trend_map <- function(slope_df, color_values, stage_label, non_sig_mask = NULL) {
  trend_plot <- ggplot() +
    geom_tile(data = slope_df, aes(x = x, y = y, fill = value_binned, alpha = significant))

  # Overlay a hatch pattern on cells with a non-significant trend
  if (!is.null(non_sig_mask) && nrow(non_sig_mask) > 0) {
    trend_plot <- trend_plot +
      ggpattern::geom_sf_pattern(
        data = non_sig_mask,
        inherit.aes = FALSE,
        fill = NA,
        colour = NA,
        pattern = "stripe",
        pattern_angle = 45,
        pattern_spacing = 0.005,
        pattern_density = 0.25,
        pattern_fill = "black",
        pattern_colour = "black",
        pattern_size = 0.1,
        pattern_alpha = 0.5,
        show.legend = FALSE)
  }

  trend_plot +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.25), guide = "none") +
    geom_sf(data = world, fill = "grey90", color = "grey20") +
    annotate("text", x = -3.4, y = 44, label = stage_label, fontface = "bold") +
    coord_sf(xlim = c(min(slope_df$x) - 0.1, max(slope_df$x) + 0.1), ylim = c(min(slope_df$y) - 0.1, max(slope_df$y) + 0.1), expand = FALSE) +
    labs(fill = "Sen's slope (ind/km²/year)", x = "Longitude", y = "Latitude") +
    scale_fill_manual(
      values = color_values,
      na.value = NA,
      na.translate = FALSE,
      limits = levels(slope_df$value_binned),
      guide = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 1, label.position = "bottom", keywidth = unit(1.2, "cm"), keyheight = unit(0.4, "cm"))) +
    theme(
      panel.grid.major = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      panel.background = element_rect(fill = "white"),
      axis.text = element_text(size = 8),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      legend.key = element_rect(color = "white"))
}

# -----------------------------------------------------------------------------
# Generate trend maps for juvenile and adult stages
# -----------------------------------------------------------------------------

juvenile_plot <- plot_trend_map(slope_df_juvenile, col_vec_juvenile, "Juveniles", non_sig_mask_juvenile)
adult_plot <- plot_trend_map(slope_df_adult, col_vec_adult, "Adults", non_sig_mask_adult)

combined_plot <- (juvenile_plot + adult_plot) +
  plot_layout(ncol = 1) +
  plot_annotation(tag_levels = "A")

# -----------------------------------------------------------------------------
# Save figure
# -----------------------------------------------------------------------------

output_file <- "./Figures/Figure_2.png"
# ggplot2::ggsave(filename = output_file, plot = combined_plot, width = 21, height = 30, units = "cm", dpi = 600)
