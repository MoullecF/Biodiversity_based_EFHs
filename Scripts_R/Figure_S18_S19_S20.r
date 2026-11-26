###############################################################################
# Title: Cross-Tabulation of EHSA Classes and Sen's Slope Trends
###############################################################################

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------

rm(list = ls())
source("Scripts_R/0-Load libraries.R")
source("Scripts_R/ggplot_theme.R")

# -----------------------------------------------------------------------------
# Load and prepare hurdle-model raster stack
# -----------------------------------------------------------------------------

hurdle_stack <- get(load("Outputs/Spatio_temporal_prediction/r_stack_Hurdle_0.05_0_1000_19992021.Rdata"))

# Load depth mask
depth_mask <- readRDS("Inputs_HMSC/depth_mask.RDS")

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

# -----------------------------------------------------------------------------
# Load Emerging Hotspot Analysis (EHSA) classifications
# -----------------------------------------------------------------------------

# Projected CRS used for area calculations (km²), matching Figure_5_6.r.
med_albers_crs <- "+proj=aea +lat_1=30 +lat_2=45 +lat_0=37.5 +lon_0=15 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

path_ehsa_juvenile <- "Outputs/EHSA/EHSA_juvenile.rds"
path_ehsa_adult <- "Outputs/EHSA/EHSA_adult.rds"

ehsa_juvenile_90 <- readRDS(path_ehsa_juvenile)
ehsa_adult_90 <- readRDS(path_ehsa_adult)

hotspot_classes_lc <- c("new hotspot", "consecutive hotspot", "intensifying hotspot",
  "persistent hotspot", "diminishing hotspot", "sporadic hotspot")

coldspot_classes_lc <- c("new coldspot", "consecutive coldspot", "intensifying coldspot",
  "persistent coldspot", "diminishing coldspot", "sporadic coldspot")

# -----------------------------------------------------------------------------
# Join EHSA classes with Sen's slope significance classes
# -----------------------------------------------------------------------------

build_ehsa_grid <- function(ehsa_sf, hotspot_classes, coldspot_classes) {
  coords <- st_coordinates(st_centroid(ehsa_sf))

  # Cell area (km²) computed on an equal-area projection
  area_km2 <- ehsa_sf %>%
    st_transform(crs = med_albers_crs) %>%
    st_area() %>%
    units::set_units("km^2") %>%
    as.numeric()

  ehsa_sf %>%
    st_drop_geometry() %>%
    mutate(
      x = round(coords[, 1], 2),
      y = round(coords[, 2], 2),
      area_km2 = area_km2,
      # Fine-grained EHSA sub-category (e.g. "new hotspot", "persistent coldspot")
      ehsa_subclass = classification,
      ehsa_class = case_when(
        classification %in% hotspot_classes ~ "Hotspot",
        classification %in% coldspot_classes ~ "Coldspot",
        TRUE ~ "No pattern")) %>%
    dplyr::select(x, y, ehsa_class, ehsa_subclass, area_km2)
}

classify_sen_slope <- function(slope_df) {
  slope_df %>%
    mutate(
      x = round(x, 2),
      y = round(y, 2),
      sen_class = case_when(
        is.na(slope) | is.na(pvalue) ~ NA_character_,
        pvalue < 0.05 & slope > 0 ~ "Significant positive",
        pvalue < 0.05 & slope < 0 ~ "Significant negative",
        TRUE ~ "Non-significant")) %>%
    dplyr::select(x, y, sen_class)
}

ehsa_grid_juvenile <- build_ehsa_grid(ehsa_juvenile_90, hotspot_classes_lc, coldspot_classes_lc)
ehsa_grid_adult <- build_ehsa_grid(ehsa_adult_90, hotspot_classes_lc, coldspot_classes_lc)

sen_class_juvenile <- classify_sen_slope(slope_df_juvenile)
sen_class_adult <- classify_sen_slope(slope_df_adult)

# Join on grid coordinates
ehsa_sen_juvenile <- inner_join(ehsa_grid_juvenile, sen_class_juvenile, by = c("x", "y")) %>%
  filter(!is.na(sen_class)) %>%
  mutate(stage = "Juveniles")

ehsa_sen_adult <- inner_join(ehsa_grid_adult, sen_class_adult, by = c("x", "y")) %>%
  filter(!is.na(sen_class)) %>%
  mutate(stage = "Adults")

# -----------------------------------------------------------------------------
# Build cross-tabulation
# -----------------------------------------------------------------------------

sen_levels <- c("Significant positive", "Non-significant", "Significant negative")
ehsa_levels <- c("Hotspot", "Coldspot", "No pattern")

# Proportion (%) of area within each EHSA class and stage falling into each Sen's slope significance class
ehsa_sen_crosstab <- bind_rows(ehsa_sen_juvenile, ehsa_sen_adult) %>%
  dplyr::mutate(
    ehsa_class = factor(ehsa_class, levels = ehsa_levels),
    sen_class = factor(sen_class, levels = sen_levels),
    stage = factor(stage, levels = c("Juveniles", "Adults"))) %>%
  dplyr::group_by(stage, ehsa_class, sen_class) %>%
  dplyr::summarise(area_km2 = sum(area_km2), n_pixels = n(), .groups = "drop_last") %>%
  dplyr::mutate(pct_area = round(100 * area_km2 / sum(area_km2), 1)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(stage, ehsa_class, sen_class)

ehsa_sen_crosstab_wide <- ehsa_sen_crosstab %>%
  dplyr::select(stage, ehsa_class, sen_class, pct_area) %>%
  tidyr::pivot_wider(names_from = sen_class, values_from = pct_area, values_fill = 0) %>%
  left_join(
    ehsa_sen_crosstab %>%
      dplyr::group_by(stage, ehsa_class) %>%
      dplyr::summarise(area_km2_total = round(sum(area_km2), 1), .groups = "drop"),
    by = c("stage", "ehsa_class")) %>%
  arrange(stage, ehsa_class)

print(ehsa_sen_crosstab_wide)

# -----------------------------------------------------------------------------
# Build cross-tabulation by fine-grained hotspot/coldspot sub-category
# -----------------------------------------------------------------------------

subclass_levels <- str_to_title(c(hotspot_classes_lc, coldspot_classes_lc))

# Proportion (%) of area within each EHSA sub-category and stage falling into each Sen's slope significance class
ehsa_sen_subclass_crosstab <- bind_rows(ehsa_sen_juvenile, ehsa_sen_adult) %>%
  filter(ehsa_class %in% c("Hotspot", "Coldspot")) %>%
  dplyr::mutate(
    ehsa_subclass = str_to_title(ehsa_subclass),
    ehsa_subclass = factor(ehsa_subclass, levels = subclass_levels),
    sen_class = factor(sen_class, levels = sen_levels),
    stage = factor(stage, levels = c("Juveniles", "Adults"))) %>%
  dplyr::group_by(stage, ehsa_subclass, sen_class) %>%
  dplyr::summarise(area_km2 = sum(area_km2), n_pixels = n(), .groups = "drop_last") %>%
  dplyr::mutate(pct_area = round(100 * area_km2 / sum(area_km2), 1)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(stage, ehsa_subclass, sen_class)

ehsa_sen_subclass_crosstab_wide <- ehsa_sen_subclass_crosstab %>%
  dplyr::select(stage, ehsa_subclass, sen_class, pct_area) %>%
  tidyr::pivot_wider(names_from = sen_class, values_from = pct_area, values_fill = 0) %>%
  left_join(
    ehsa_sen_subclass_crosstab %>%
      dplyr::group_by(stage, ehsa_subclass) %>%
      dplyr::summarise(area_km2_total = round(sum(area_km2), 1), .groups = "drop"),
    by = c("stage", "ehsa_subclass")) %>%
  arrange(stage, ehsa_subclass)

print(ehsa_sen_subclass_crosstab_wide)

# -----------------------------------------------------------------------------
# Plot cross-tabulation as stacked bar plot
# -----------------------------------------------------------------------------

sen_fill_colors <- c("Significant positive" = "#C00000", "Non-significant" = "grey80", "Significant negative" = "#08519C")

ehsa_sen_barplot <- ggplot(ehsa_sen_crosstab, aes(x = ehsa_class, y = pct_area, fill = sen_class)) +
  geom_col(position = "stack") +
  facet_wrap(~ stage) +
  scale_fill_manual(values = sen_fill_colors, name = "Sen's slope") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Proportion of area (%)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"))

# Sub-category version
ehsa_sen_subclass_barplot <- ggplot(ehsa_sen_subclass_crosstab, aes(x = ehsa_subclass, y = pct_area, fill = sen_class)) +
  geom_col(position = "stack") +
  facet_wrap(~ stage) +
  scale_fill_manual(values = sen_fill_colors, name = "Sen's slope") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Proportion of area (%)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"))

# Combined plot of broad EHSA classes and fine-grained sub-categories for each stage.
combined_crosstab_plot <- (ehsa_sen_barplot + ehsa_sen_subclass_barplot) +
  plot_layout(ncol = 1, heights = c(1, 1.2), guides = "collect") &
  theme(legend.position = "bottom")
combined_crosstab_plot <- combined_crosstab_plot +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16, face = "bold"))

# ggplot2::ggsave(combined_crosstab_plot, filename = "Figures/Figure_S18.png", width = 26, height = 20, units = "cm", dpi = 400)

# -----------------------------------------------------------------------------
# Sen's slope magnitude bins, by hotspot/coldspot sub-category and stage
# -----------------------------------------------------------------------------

adult_slope_breaks <- c(-Inf, -142, -112, -59, -37, -21, -11, -5, 0, 7, Inf)
adult_slope_labels <- c("< -142", "[-142;-112[", "[-112;-59[", "[-59;-37[", "[-37;-21[",
  "[-21;-11[", "[-11;-5[", "[-5;0[", "[0;+7[", "> +7")

juvenile_slope_breaks <- c(-Inf, -16, -10, -7, -4, -1, 0, 1, 2, 4, Inf)
juvenile_slope_labels <- c("< -16", "[-16;-10[", "[-10;-7[", "[-7;-4[", "[-4;-1[",
  "[-1;0[", "[0;+1[", "[+1;+2[", "[+2;+4[", "> +4")

bin_slope_values <- function(slope_df, breaks, labels) {
  slope_df %>%
    mutate(x = round(x, 2), y = round(y, 2), slope_bin = cut(slope, breaks = breaks, labels = labels, include.lowest = TRUE)) %>%
    dplyr::select(x, y, slope_bin)
}

slope_bin_adult <- bin_slope_values(slope_df_adult, adult_slope_breaks, adult_slope_labels)
slope_bin_juvenile <- bin_slope_values(slope_df_juvenile, juvenile_slope_breaks, juvenile_slope_labels)

# Proportion (%) of area within each EHSA sub-category for a given stage falling into each Sen's slope magnitude bin
summarise_slope_bins_by_subclass <- function(ehsa_grid, slope_bin_df, slope_labels,
                                              ehsa_class_label, subclass_levels) {
  ehsa_grid %>%
    filter(ehsa_class == ehsa_class_label) %>%
    inner_join(slope_bin_df, by = c("x", "y")) %>%
    filter(!is.na(slope_bin)) %>%
    mutate(ehsa_subclass = str_to_title(ehsa_subclass),
      ehsa_subclass = factor(ehsa_subclass, levels = str_to_title(subclass_levels)),
      slope_bin = factor(slope_bin, levels = slope_labels)) %>%
    group_by(ehsa_subclass, slope_bin, .drop = FALSE) %>%
    summarise(area_km2 = sum(area_km2), n_pixels = n(), .groups = "drop") %>%
    group_by(ehsa_subclass) %>%
    mutate(area_km2_total = sum(area_km2),
      pct_area = if_else(area_km2_total > 0, round(100 * area_km2 / area_km2_total, 1), NA_real_)) %>%
    ungroup() %>%
    arrange(ehsa_subclass, slope_bin)
}

pivot_slope_bin_wide <- function(slope_bin_summary) {
  slope_bin_summary %>%
    dplyr::select(ehsa_subclass, slope_bin, pct_area, area_km2_total) %>%
    tidyr::pivot_wider(names_from = slope_bin, values_from = pct_area, values_fill = 0) %>%
    arrange(ehsa_subclass)
}

# Adults
hotspot_slope_bin_subclass_summary_adult <- summarise_slope_bins_by_subclass(
  ehsa_grid_adult, slope_bin_adult, adult_slope_labels, "Hotspot", hotspot_classes_lc)

coldspot_slope_bin_subclass_summary_adult <- summarise_slope_bins_by_subclass(
  ehsa_grid_adult, slope_bin_adult, adult_slope_labels, "Coldspot", coldspot_classes_lc)

hotspot_slope_bin_subclass_wide_adult <- pivot_slope_bin_wide(hotspot_slope_bin_subclass_summary_adult)
coldspot_slope_bin_subclass_wide_adult <- pivot_slope_bin_wide(coldspot_slope_bin_subclass_summary_adult)

print(hotspot_slope_bin_subclass_wide_adult)
print(coldspot_slope_bin_subclass_wide_adult)

# Juveniles
hotspot_slope_bin_subclass_summary_juvenile <- summarise_slope_bins_by_subclass(
  ehsa_grid_juvenile, slope_bin_juvenile, juvenile_slope_labels, "Hotspot", hotspot_classes_lc)

coldspot_slope_bin_subclass_summary_juvenile <- summarise_slope_bins_by_subclass(
  ehsa_grid_juvenile, slope_bin_juvenile, juvenile_slope_labels, "Coldspot", coldspot_classes_lc)

hotspot_slope_bin_subclass_wide_juvenile <- pivot_slope_bin_wide(hotspot_slope_bin_subclass_summary_juvenile)
coldspot_slope_bin_subclass_wide_juvenile <- pivot_slope_bin_wide(coldspot_slope_bin_subclass_summary_juvenile)

print(hotspot_slope_bin_subclass_wide_juvenile)
print(coldspot_slope_bin_subclass_wide_juvenile)

# -----------------------------------------------------------------------------
# Plot slope-magnitude bins as stacked bar charts, by sub-category and stage
# -----------------------------------------------------------------------------

slope_bin_fill_colors <- function(slope_labels) {
  n_bins <- length(slope_labels)
  setNames(colorRampPalette(c("#08519C", "grey90", "#C00000"))(n_bins), slope_labels)
}

plot_slope_bin_barplot <- function(slope_bin_summary, slope_labels, y_label) {
  ggplot(slope_bin_summary, aes(x = ehsa_subclass, y = pct_area, fill = slope_bin)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = slope_bin_fill_colors(slope_labels), name = "Sen's slope bin", drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    labs(x = NULL, y = y_label) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"))
}

hotspot_slope_bin_barplot_adult <- plot_slope_bin_barplot(hotspot_slope_bin_subclass_summary_adult, adult_slope_labels, "Proportion of adult hotspot area (%)")
coldspot_slope_bin_barplot_adult <- plot_slope_bin_barplot(coldspot_slope_bin_subclass_summary_adult, adult_slope_labels, "Proportion of adult coldspot area (%)")

hotspot_slope_bin_barplot_juvenile <- plot_slope_bin_barplot(hotspot_slope_bin_subclass_summary_juvenile, juvenile_slope_labels, "Proportion of juvenile hotspot area (%)")
coldspot_slope_bin_barplot_juvenile <- plot_slope_bin_barplot(coldspot_slope_bin_subclass_summary_juvenile, juvenile_slope_labels, "Proportion of juvenile coldspot area (%)")

# Combined plot of hotspot and coldspot slope-magnitude bins for each stage
combined_slope_bin_plot_adult <- (hotspot_slope_bin_barplot_adult + coldspot_slope_bin_barplot_adult) +
  plot_layout(ncol = 1, heights = c(1, 1.2), guides = "collect") &
  theme(legend.position = "bottom") &
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 1, label.position = "bottom", keywidth = unit(1.2, "cm"), keyheight = unit(0.4, "cm"))) &
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16, face = "bold"))

# Save
# ggplot2::ggsave(combined_slope_bin_plot_adult, filename = "Figures/Figure_S19.png", width = 26, height = 25, units = "cm", dpi = 400)

combined_slope_bin_plot_juvenile <- (hotspot_slope_bin_barplot_juvenile + coldspot_slope_bin_barplot_juvenile) +
  plot_layout(ncol = 1, heights = c(1, 1.2), guides = "collect") &
  theme(legend.position = "bottom") &
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 1, label.position = "bottom", keywidth = unit(1.2, "cm"), keyheight = unit(0.4, "cm"))) &
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16, face = "bold"))

# Save
# ggplot2::ggsave(combined_slope_bin_plot_juvenile, filename = "Figures/Figure_S20.png", width = 26, height = 25, units = "cm", dpi = 400)

