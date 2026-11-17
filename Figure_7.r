###############################################################################
# Title: MPA Coverage of EHSA Hotspots (Figure 7) and of juvenile and adult
#        hotspots (test of hypothesis 2, "early life stages")
###############################################################################

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------

rm(list = ls())

source("./0-Load libraries.r")
source("./ggplot_theme.R")

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------

# Juvenile and adult EHSA outputs
file_ehsa_juvenile <- "./Outputs/EHSA/EHSA_juvenile.RDS"
file_ehsa_adult    <- "./Outputs/EHSA/EHSA_adult.RDS"

class_col <- "classification"
protection_levels <- c("Unprotected", "Minimally", "Lightly", "Highly", "Fully")

# -----------------------------------------------------------------------------
# Load spatial data
# -----------------------------------------------------------------------------
 
# Marine protected areas with protection levels minimally --> Fully
mpa_polygons <- st_read("./Inputs_HMSC/protectionPLactivities_niv1&2&3&4/protectionPLactivities_niv1&2&3&4.shp")
mpa_polygons <- st_make_valid(mpa_polygons)
 
# Recode numeric protection codes
mpa_polygons <- mpa_polygons %>%
  mutate(protctn = factor(protctn, levels = c(1, 2, 3, 4), labels = c("Minimally", "Lightly", "Highly", "Fully")))
 
# Define color palette for protection levels
protection_palette <- c(
  "Minimally"   = "#fde725",
  "Lightly"     = "#5ec962",
  "Highly"      = "#21918c",
  "Fully"       = "#440154")
 
# Plot MPAs with protection levels
map_MPA <- ggplot() +
  geom_sf(data = st_transform(mpa_polygons, st_crs(world)), aes(fill = protctn), color = NA) +
  geom_sf(data = world, fill = "grey90", color = "grey20") +
  coord_sf(xlim = c(-6, 16.5), ylim = c(35, 45.5), expand = FALSE) +
  scale_fill_manual(values = protection_palette, drop = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  guides(fill = guide_legend(title = "Protection\nlevel")) +
  theme(
    panel.grid.major = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.background = element_rect(fill = "white"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    legend.key = element_rect(color = "white"),
    plot.margin = unit(c(0, 0, 0, 0), "pt"))
 
# ggsave(map_MPA, filename = "./Figures/Map_MPA_protection_levels.png", width = 20, height = 20, units = "cm", dpi = 400)
 
# Combined hotspot and coldspot EHSA outputs
df_spots <- readRDS("./Outputs/EHSA/Combined_HS_CS.RDS")
 
# Align CRS for spatial operations
df_spots <- st_transform(df_spots, st_crs(mpa_polygons))
 
# -----------------------------------------------------------------------------
# Function: area of each feature
# -----------------------------------------------------------------------------
 
# Intersections can leave lines or points where a cell edge touches an MPA boundary
# These have no surface so only the polygonal features are measured and the others are given an area of 0
area_m2 <- function(x) {
  g <- st_geometry(x)
  a <- numeric(length(g))
  d <- st_dimension(g)
  ok <- !is.na(d) & d == 2
  if (any(ok)) a[ok] <- as.numeric(st_area(st_make_valid(g[ok])))
  a
}
 
# -----------------------------------------------------------------------------
# Function: surface of a set of cells in each protection level
# -----------------------------------------------------------------------------
 
# Returns one row per group x protection level, with the area (m2) and the percentage of the group surface
mpa_overlap_summary <- function(cells, group_col, mpa) {
  cells <- st_transform(cells, st_crs(mpa))
  cells <- cells[, group_col]
  names(cells)[1] <- "group"
  cells$group <- as.character(cells$group)
 
  protected <- suppressWarnings(st_intersection(cells, mpa[, "protctn"]))
  protected$area <- area_m2(protected)
  protected <- protected %>%
    st_drop_geometry() %>%
    mutate(protctn = as.character(protctn)) %>%
    group_by(group, protctn) %>%
    summarise(total_area = sum(area), .groups = "drop")
 
  unprotected <- suppressWarnings(st_difference(cells, st_union(st_geometry(mpa))))
  unprotected$area <- area_m2(unprotected)
  unprotected <- unprotected %>%
    st_drop_geometry() %>%
    group_by(group) %>%
    summarise(total_area = sum(area), .groups = "drop") %>%
    mutate(protctn = "Unprotected")
 
  bind_rows(protected, unprotected) %>%
    group_by(group) %>%
    mutate(prop = round(100 * total_area / sum(total_area), 2)) %>%
    ungroup() %>%
    mutate(protctn = factor(protctn, levels = protection_levels)) %>%
    arrange(group, protctn)
}
 
# -----------------------------------------------------------------------------
# Aggregate protection level summaries (congruent hotspots and coldspots)
# -----------------------------------------------------------------------------
 
summary_all <- mpa_overlap_summary(df_spots, "Spot", mpa_polygons) %>%
  rename(Spot = group)
 
print(summary_all, n = 100)
 
# -----------------------------------------------------------------------------
# Juvenile and adult hotspots --> surface in each protection level
# -----------------------------------------------------------------------------
 
ehsa_juv <- readRDS(file_ehsa_juvenile)
ehsa_adu <- readRDS(file_ehsa_adult)
 
# Both layers must hold the same cells in the same order
eq <- st_equals(ehsa_juv, ehsa_adu)
stopifnot(nrow(ehsa_juv) == nrow(ehsa_adu), all(lengths(eq) == 1), all(unlist(eq) == seq_len(nrow(ehsa_juv))))
 
# Broad class from the EHSA class label
broad_class <- function(x) {
  dplyr::case_when(
    grepl("hotspot",  x, ignore.case = TRUE) ~ "Hotspot",
    grepl("coldspot", x, ignore.case = TRUE) ~ "Coldspot",
    TRUE ~ "No pattern")
}
 
cells_stage <- ehsa_juv[, class_col]
cells_stage$juv <- broad_class(as.character(ehsa_juv[[class_col]]))
cells_stage$adu <- broad_class(as.character(ehsa_adu[[class_col]]))
 
cells_stage <- cells_stage %>%
  mutate(
    g_juv = paste("Juvenile", tolower(juv)),
    g_adu = paste("Adult", tolower(adu)),
    g_combo = dplyr::case_when(
      juv == "Hotspot" & adu == "Hotspot" ~ "Congruent hotspot",
      juv == "Hotspot" & adu != "Hotspot" ~ "Juvenile-only hotspot",
      juv != "Hotspot" & adu == "Hotspot" ~ "Adult-only hotspot",
      TRUE ~ "Not a hotspot for either stage"),
    g_all = "Study area")
 
# Number of cells per group
print(table(cells_stage$g_juv))
print(table(cells_stage$g_adu))
print(table(cells_stage$g_combo))
 
mpa_by_group <- bind_rows(
  mpa_overlap_summary(cells_stage, "g_juv",   mpa_polygons) %>% mutate(set = "Juveniles"),
  mpa_overlap_summary(cells_stage, "g_adu",   mpa_polygons) %>% mutate(set = "Adults"),
  mpa_overlap_summary(cells_stage, "g_combo", mpa_polygons) %>% mutate(set = "Juveniles x adults"),
  mpa_overlap_summary(cells_stage, "g_all",   mpa_polygons) %>% mutate(set = "Study area"))
 
# Wide table
mpa_by_group_wide <- mpa_by_group %>%
  group_by(set, group) %>%
  mutate(area_km2 = round(sum(total_area) / 1e6)) %>%
  ungroup() %>%
  dplyr::select(set, group, area_km2, protctn, prop) %>%
  tidyr::pivot_wider(names_from = protctn, values_from = prop, values_fill = 0, names_expand = TRUE) %>%
  mutate(Any_protection = round(100 - Unprotected, 2), Highly_or_fully = round(Highly + Fully, 2))
 
print(as.data.frame(mpa_by_group_wide))
 
# saveRDS(mpa_by_group, "./Outputs/MPA_overlap_by_life_stage.RDS")
# write.csv(mpa_by_group_wide, "./Outputs/MPA_overlap_by_life_stage.csv", row.names = FALSE)
 
# -----------------------------------------------------------------------------
# Per-cell protection coverage
# -----------------------------------------------------------------------------
 
# Compute proportional MPA coverage within each hotspot/coldspot cell
# Assign unique IDs and compute original cell areas
df_spots <- df_spots %>%
  mutate(id = row_number(), spot_area = area_m2(geometry))
 
# Intersect each spot cell with MPAs and compute proportional coverage
cell_intersections <- suppressWarnings(st_intersection(df_spots, mpa_polygons)) %>%
  mutate(
    intersection_area = area_m2(geometry), # Area of overlap
    prop_cover = intersection_area / spot_area ) # Fraction covered
 
# Pivot to wide format
# If multiple MPA polygons overlap a cell --> sum their contributions
coverage_wide <- cell_intersections %>%
  st_set_geometry(NULL) %>%
  group_by(id, protctn) %>%
  summarise(proportion = sum(prop_cover) * 100, .groups = "drop") %>%
  pivot_wider(names_from = protctn, values_from = proportion, values_fill = 0)
 
# Merge coverage back into spatial dataframe
df_spots_with_cover <- df_spots %>%
  left_join(coverage_wide, by = "id")
 
# -----------------------------------------------------------------------------
# Prepare hotspot cells for spatial visualization
# -----------------------------------------------------------------------------
 
# Filter hotspots and compute unprotected proportion
hotspots <- df_spots_with_cover %>%
  filter(Spot == "Hotspots")
 
# Identify new protection columns added by the coverage join
protctn_cols <- setdiff(names(hotspots), names(df_spots))
 
# Remainder after summing all protected proportions
hotspots <- hotspots %>%
  mutate(unprotected = 100 - rowSums(across(all_of(protctn_cols)), na.rm = TRUE))
 
# Reshape to long format for mapping by protection category
hotspots_long <- hotspots %>%
  dplyr::select(id, geometry, all_of(protctn_cols), unprotected) %>%
  pivot_longer(
    cols = c(all_of(protctn_cols), unprotected),
    names_to = "category",
    values_to = "proportion") %>%
  mutate(category = dplyr::recode_values(category, "unprotected" ~ "Unprotected", default = category),
    proportion = pmax(proportion, 0, na.rm = TRUE)) %>% # Clamp negatives to 0
  st_transform(4326)  # WGS84 for display
 
# -----------------------------------------------------------------------------
# Dominant protection map
# -----------------------------------------------------------------------------
 
# Classify each hotspot cell by the protection level with the highest
# proportional coverage, including Unprotected as the lowest "protection level"
dominant_cols <- c(protctn_cols, "unprotected")
# Preference order
dominant_priority <- intersect(c("Fully", "Highly", "Lightly", "Minimally", "unprotected"), dominant_cols)
 
hotspots_dominant <- hotspots %>%
  mutate(across(all_of(dominant_cols), ~ coalesce(.x, 0))) %>%
  mutate(
    max_cover = pmax(!!!syms(dominant_cols), na.rm = TRUE),
    # Pick the first column that attains the row max
    dominant = dominant_priority[max.col(as.matrix(st_drop_geometry(pick(all_of(dominant_priority)))) == max_cover, ties.method = "first")],
    dominant = dplyr::recode_values(dominant, "unprotected" ~ "Unprotected", default = dominant),
    dominant = factor(dominant, levels = protection_levels)) %>%
  st_transform(4326)
 
dominant_count <- hotspots_dominant %>%
  st_drop_geometry() %>%
  count(dominant, name = "n_cells") %>%
  mutate(prop_cells = round((n_cells / sum(n_cells)) * 100, 2))
 
print(dominant_count)
 
# -----------------------------------------------------------------------------
# Palette
# -----------------------------------------------------------------------------
 
protection_palette <- c(
  "Unprotected" = "grey60",
  "Minimally"   = "#fde725",
  "Lightly"     = "#5ec962",
  "Highly"      = "#21918c",
  "Fully"       = "#440154")
 
dominant_levels <- levels(hotspots_dominant$dominant)
dominant_palette_ordered <- protection_palette[dominant_levels]
 
# -----------------------------------------------------------------------------
# Dummy geometry to force legend entry
# -----------------------------------------------------------------------------
legend_dummy_dominant <- hotspots_dominant %>%
  slice(1) %>%
  mutate(dominant = factor("Fully", levels = dominant_levels))
 
dominant_map <- ggplot() +
  geom_sf(data = hotspots_dominant, aes(fill = dominant), color = NA, linewidth = 0.05) +
  geom_sf(data = legend_dummy_dominant, aes(fill = dominant), alpha = 0, show.legend = TRUE) +
  geom_sf(data = world, fill = "grey90", color = "grey20") +
  coord_sf(xlim = c(st_bbox(hotspots_long)[1] - 0.1, st_bbox(hotspots_long)[3] + 0.15), ylim = c(st_bbox(hotspots_long)[2] - 0.1, st_bbox(hotspots_long)[4] + 0.4), expand = FALSE) +
  scale_fill_manual(values = dominant_palette_ordered, breaks = dominant_levels, limits = dominant_levels, drop = FALSE, name = "Dominant protection (highest % of grid cell)", guide = guide_legend(override.aes = list(alpha = 1))) +
  labs(x = "Longitude", y = "Latitude") +
  theme(
    panel.grid.major = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.background = element_rect(fill = "white"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    legend.key = element_rect(color = "white"),
    plot.margin = unit(c(0, 0, 0, 0), "pt"))
 
# Horizontal barplot showing hotspot surface by protection level
bar_data <- summary_all[summary_all$Spot == "Hotspots", ] %>%
  tidyr::complete(protctn, fill = list(prop = 0, total_area = 0, Spot = "Hotspots"))
 
prop_of <- function(level) {
  v <- bar_data$prop[bar_data$protctn == level]
  if (length(v) == 0) 0 else v
}
bar_label <- function(level, text) {
  paste0(text, " (~", format(round(prop_of(level), 1), nsmall = 1), "%)")
}
 
bar_breaks <- round(sapply(c("Fully", "Highly", "Lightly", "Minimally", "Unprotected"), prop_of))
bar_breaks <- sort(unique(c(0, bar_breaks[bar_breaks >= 3])))
bar_break_labels <- ifelse(bar_breaks == 0, "0", paste0(bar_breaks, "%"))
bar_xmax <- max(70, ceiling(prop_of("Unprotected")) + 2)
 
bar_hotspots <- ggplot(bar_data, aes(x = prop, y = protctn, fill = protctn)) +
  geom_bar(stat = "identity", colour = "white", linewidth = 0.4, show.legend = FALSE) +
  scale_fill_manual(values = protection_palette, drop = FALSE) +
  scale_x_continuous(breaks = bar_breaks, labels = bar_break_labels, limits = c(0, bar_xmax), expand = expansion(mult = c(0, 0.02))) +
  scale_y_discrete(drop = FALSE) +
  labs(x = "Proportion of hotspot surface", y = element_blank()) +
  annotate("text", x = 14, y = 5, label = bar_label("Fully", "Fully protected"), size = 13 / .pt) +
  annotate("text", x = 16, y = 4, label = bar_label("Highly", "Highly protected"), size = 13 / .pt) +
  annotate("text", x = 25, y = 3, label = bar_label("Lightly", "Lightly protected"), size = 13 / .pt) +
  annotate("text", x = 39, y = 2, label = bar_label("Minimally", "Minimally protected"), size = 13 / .pt) +
  annotate("text", x = 50, y = 1, label = bar_label("Unprotected", "Unprotected"), size = 13 / .pt) +
  HMSC.theme +
  theme(
    axis.title.x = element_text(face = "plain", size = 14),
    axis.text.x = element_text(face = "plain", size = 12),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "pt"),
    panel.spacing = unit(0, "pt"),
    aspect.ratio = 1.06)
 
dominant_with_bar <- (dominant_map | bar_hotspots) +
  plot_layout(guides = "collect", widths = c(2, 1)) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title.position = "top",
    legend.margin = margin(),
    legend.box.spacing = unit(2, "pt"),
    plot.margin = unit(c(0, 0, 0, 0), "pt"))

# Save
# ggplot2::ggsave(dominant_with_bar, filename = "./Figures/Figure_7.png", width = 45, height = 18, units = "cm", dpi = 400)