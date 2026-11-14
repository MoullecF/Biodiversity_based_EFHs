###############################################################################
# Title: MPA Coverage of EHSA Hotspots (Figure 7)
###############################################################################

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------

rm(list = ls())

source("./0-Load libraries.r")
source("./ggplot_theme.R")

# -----------------------------------------------------------------------------
# Load spatial data
# -----------------------------------------------------------------------------

# Marine protected areas with protection levels minimally --> Fully.
mpa_polygons <- st_read("./Inputs_HMSC/protectionPLactivities_niv1&2&3&4/protectionPLactivities_niv1&2&3&4.shp")
mpa_polygons <- st_make_valid(mpa_polygons)

# Recode numeric protection codes to labelled factor matching the palette names
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
  coord_sf(xlim = c(-6, 16.5), ylim = c(35, 45.5), expand = FALSE ) +
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
# Aggregate protection level summaries
# -----------------------------------------------------------------------------

# Calculate total hotspot/coldspot surface in each protection category
# Intersect spots with MPAs and compute total areas by protection level
protected_intersections <- st_intersection(
  df_spots %>% dplyr::select(Spot),  # Keep only spot type column
  mpa_polygons %>% dplyr::select(protctn)) %>%  # Keep only protection level
  mutate(area = as.numeric(st_area(geometry)))

# Aggregate by spot type and protection level
protected_areas <- protected_intersections %>%
  st_drop_geometry() %>%
  group_by(Spot, protctn) %>%
  summarise(total_area = sum(area), .groups = "drop") %>%
  mutate(protctn = as.character(protctn))

# Compute areas outside any MPA (unprotected)
unprotected_areas <- st_difference(df_spots, st_union(mpa_polygons)) %>%
  mutate(area = as.numeric(st_area(geometry)), protctn = "Unprotected") %>%
  st_drop_geometry() %>%
  dplyr::select(Spot, protctn, total_area = area)

# Combine protected and unprotected summaries and compute proportions
summary_all <- bind_rows(protected_areas, unprotected_areas) %>%
  group_by(Spot) %>%
  mutate(prop = (total_area / sum(total_area, na.rm = TRUE)) * 100) %>%
  group_by(Spot, protctn, .add = TRUE) %>%
  summarise(
    total_area = sum(total_area, na.rm = TRUE),
    prop = round(sum(prop, na.rm = TRUE), 2),
    .groups = "drop") %>%
  mutate(protctn = factor(protctn, levels = c("Unprotected", "Minimally", "Lightly", "Highly", "Fully")))

# -----------------------------------------------------------------------------
# Per-cell protection coverage
# -----------------------------------------------------------------------------

# Compute proportional MPA coverage within each hotspot/coldspot cell
# Assign unique IDs and compute original cell areas
df_spots <- df_spots %>%
  mutate(id = row_number(), spot_area = st_area(geometry))

# Intersect each spot cell with MPAs and compute proportional coverage
cell_intersections <- st_intersection(df_spots, mpa_polygons) %>%
  mutate(
    intersection_area = st_area(geometry),  # Area of overlap
    prop_cover = as.numeric(intersection_area / spot_area)  # Fraction covered
  )

# Pivot to wide format: one column per protection level
# If multiple MPA polygons overlap a cell, sum their contributions
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
# One row per cell × category combination
hotspots_long <- hotspots %>%
  dplyr::select(id, geometry, all_of(protctn_cols), unprotected) %>%
  pivot_longer(
    cols = c(all_of(protctn_cols), unprotected),
    names_to = "category",
    values_to = "proportion") %>%
  mutate(category = recode(category, "unprotected" = "Unprotected"),
    proportion = pmax(proportion, 0, na.rm = TRUE)) %>%  # Clamp negatives to 0
  st_transform(4326)  # WGS84 for display

# -----------------------------------------------------------------------------
# Dominant protection map (highest coverage)
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
    dominant = recode(dominant, "unprotected" = "Unprotected"),
    dominant = factor(dominant, levels = c("Unprotected", "Minimally", "Lightly", "Highly", "Fully"))) %>%
  st_transform(4326)

dominant_count <- hotspots_dominant %>%
  st_drop_geometry() %>%
  count(dominant, name = "n_cells") %>%
  mutate(prop_cells = round((n_cells / sum(n_cells)) * 100, 2))

# -----------------------------------------------------------------------------
# Palette (ALL levels included)
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
  coord_sf(xlim = c(st_bbox(hotspots_long)[1] - 0.1,
             st_bbox(hotspots_long)[3] + 0.15),
    ylim = c(st_bbox(hotspots_long)[2] - 0.1,
             st_bbox(hotspots_long)[4] + 0.15), expand = FALSE) +
  scale_fill_manual(values = dominant_palette_ordered,
    breaks = dominant_levels, limits = dominant_levels,
    drop = FALSE, name = "Dominant protection (highest % of grid cell)",
    guide = guide_legend(override.aes = list(alpha = 1))) +
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

# Horizontal bar chart showing hotspot surface by protection level.
bar_hotspots <- ggplot(summary_all[summary_all$Spot == "Hotspots", ], aes(x = prop, y = protctn, fill = protctn)) +
  geom_bar(stat = "identity", colour = "white", linewidth = 0.4, show.legend = FALSE) +
  scale_fill_manual(values = protection_palette, drop = FALSE) +
  scale_x_continuous(breaks = c(0, 1, 9, 21, 68), labels = c("0", "1%", "9%", "21%", "68%"),
    limits = c(0, 70), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Proportion of hotspot surface", y = element_blank()) +
  annotate("text", x = 14, y = 5, label = "Fully protected (~0.1%)", size = 13 / .pt) +
  annotate("text", x = 16, y = 4, label = "Highly protected (~1.3%)", size = 13 / .pt) +
  annotate("text", x = 25, y = 3, label = "Lightly protected (~9.3%)", size = 13 / .pt) +
  annotate("text", x = 39, y = 2, label = "Minimally protected (~21.2%)", size = 13 / .pt) +
  annotate("text", x = 50, y = 1, label = "Unprotected (~68.1%)", size = 13 / .pt) +
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
  plot_layout(guides = "collect", widths = c(2,1)) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title.position = "top",
    legend.margin = margin(),
    legend.box.spacing = unit(2, "pt"),
    plot.margin = unit(c(0, 0, 0, 0), "pt"))
dominant_with_bar

# Save
ggplot2::ggsave(dominant_with_bar, filename = "./Figures/Figure_7.png", width = 45, height = 18, units = "cm", dpi = 400)
