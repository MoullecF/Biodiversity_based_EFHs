###############################################################################
# Title: EHSA Hotspot and Coldspot Overlap Maps (Figures 5-6 and Table in appendix)
###############################################################################

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------

rm(list = ls())

source("./0-Load libraries.R")
source("./ggplot_theme.R")

# -----------------------------------------------------------------------------
# Global settings and inputs
# -----------------------------------------------------------------------------

# Projected CRS used for area calculations (km²)
med_albers_crs <- "+proj=aea +lat_1=30 +lat_2=45 +lat_0=37.5 +lon_0=15 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

path_ehsa_juvenile <- "./Outputs/EHSA/EHSA_juvenile.rds"
path_ehsa_adult <- "./Outputs/EHSA/EHSA_adult.rds"

ehsa_juvenile_90 <- readRDS(path_ehsa_juvenile)
ehsa_adult_90 <- readRDS(path_ehsa_adult)

classification_order <- c("new hotspot", "consecutive hotspot", "intensifying hotspot",
  "persistent hotspot", "diminishing hotspot", "sporadic hotspot",
  "no pattern detected", "sporadic coldspot", "diminishing colspot",
  "persistent coldspot", "intensifying coldspot", "consecutive coldspot", "new coldspot")

# Target spot classes for each group (excludes "No Pattern Detected")
hotspot_classes <- c("New Hotspot", "Consecutive Hotspot", "Intensifying Hotspot", "Persistent Hotspot",
  "Diminishing Hotspot", "Sporadic Hotspot")

coldspot_classes <- c("New Coldspot", "Consecutive Coldspot", "Intensifying Coldspot", "Persistent Coldspot",
  "Diminishing Coldspot", "Sporadic Coldspot")

# Exclude combined classes tied to "no pattern detected" 
npd_remove <- c("J.New - A.Npd", "J.Int - A.Npd", "J.Per - A.Npd", "J.dim - A.Npd", "J.Spo - A.Npd",
  "J.Npd - A.Spo", "J.Npd - A.Per", "J.Npd - A.Int", "J.Npd - A.New", "J.Npd - A.Con",
  "J.Con - A.Npd")

# Factor level order for hotspot combinations
hs_levels <- c("J.Int - A.Int", "J.Int - A.Per", "J.Int - A.Spo",
  "J.Per - A.Int", "J.Per - A.Per", "J.Per - A.Spo",
  "J.Dim - A.Int", "J.Dim - A.Per", "J.Dim - A.Spo",
  "J.Spo - A.Int", "J.Spo - A.Per", "J.Spo - A.Spo",
  "J.Con - A.Per", "J.New - A.Int", "J.New - A.Per", "J.New - A.Spo",
  "J.Npd - A.Npd")

# Factor level order for coldspot combinations
cs_levels <- c("J.Int - A.Int", "J.Int - A.Per", "J.Int - A.Spo", "J.Int - A.Con", "J.Int - A.New",
  "J.Per - A.Int", "J.Per - A.Per", "J.Per - A.Spo", "J.Per - A.Con", "J.Per - A.New",
  "J.Spo - A.Int", "J.Spo - A.Per", "J.Spo - A.Spo", "J.Spo - A.Con", "J.Spo - A.New",
  "J.Con - A.Int", "J.Con - A.Per", "J.Con - A.Spo", "J.Con - A.Con", "J.Con - A.New",
  "J.New - A.Int", "J.New - A.Per", "J.New - A.Spo", "J.New - A.Con", "J.New - A.New",
  "J.Npd - A.Npd")

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

format_ehsa <- function(ehsa_df, order_levels) {
  ehsa_df %>%
    as.data.frame() %>%
    mutate(classification = str_to_title(classification),
      classification = factor(classification, levels = str_to_title(order_levels))) %>%
    arrange(classification) %>%
    st_as_sf()
}

shorten_classification <- function(classification, prefix, spot_type) {
    replacements <- c(
    "No Pattern Detected" = paste0(prefix, ".Npd"),
    setNames(paste0(prefix, ".New"), paste0("New ", spot_type)),
    setNames(paste0(prefix, ".Con"), paste0("Consecutive ", spot_type)),
    setNames(paste0(prefix, ".Int"), paste0("Intensifying ", spot_type)),
    setNames(paste0(prefix, ".Per"), paste0("Persistent ", spot_type)),
    setNames(paste0(prefix, ".Dim"), paste0("Diminishing ", spot_type)),
    setNames(paste0(prefix, ".Spo"), paste0("Sporadic ", spot_type)))
  str_replace_all(classification, replacements)
}

combine_stage_classes <- function(juvenile_sf, adult_sf, drop_patterns = NULL) {
    out <- st_join(juvenile_sf, adult_sf, largest = TRUE) %>%
    mutate(Classification = paste(classification.x, classification.y, sep = " - ")) %>%
    filter(!str_detect(Classification, "- NA"))
  if (!is.null(drop_patterns)) {
    out <- out %>% filter(!Classification %in% drop_patterns)
  }
  out %>% mutate(Classification = as.factor(Classification))
}

join_combo <- function(juv_sf, adu_sf, combination_label, pattern_status_label) {
    combine_stage_classes(juv_sf, adu_sf) %>%
    mutate(Combination    = combination_label,
      Pattern_status = pattern_status_label)
}

summarize_class_area <- function(sf_data, spot_label, total_surface_area_km2,
                                  group_vars = c("Spot", "Classification")) {
   sf_data %>%
    dplyr::mutate(area_km2 = st_area(geometry) %>% units::set_units("km^2") %>% as.numeric(),
      Spot = spot_label) %>%
    st_drop_geometry() %>%
    dplyr::group_by(across(all_of(group_vars))) %>%
    dplyr::summarise(total_area_km2 = sum(area_km2), n_polygons = n(), .groups = "drop") %>%
    dplyr::mutate(pct_total = (total_area_km2 / total_surface_area_km2) * 100,
      total_area_km2 = round(total_area_km2, 3),
      pct_total = round(pct_total, 3)) %>%
    dplyr::arrange(desc(total_area_km2))
}

collapse_small_classes <- function(area_table) {
  area_table %>%
    dplyr::mutate(classification_group = ifelse(pct_total < 1, "Others (< 1%)", as.character(Classification))) %>%
    dplyr::group_by(classification_group, Spot) %>%
    dplyr::summarize(total_area_km2 = sum(total_area_km2),
      n_polygons = sum(n_polygons),
      pct_total = sum(pct_total),
      .groups = "drop") %>%
    dplyr::rename(Classification = classification_group) %>%
    dplyr::mutate(Classification = as.factor(Classification))
}

make_spot_map <- function(spot_df, palette, label_text) {
  ggplot(data = world) +
    geom_sf(data = spot_df, aes(fill = Classification), color = "white", lwd = 0.05, inherit.aes = TRUE) +
    scale_fill_manual(values = palette) +
    geom_sf(fill = "grey90", color = "grey20") +
    coord_sf(xlim = c(st_bbox(spot_df)[1] - 0.1, st_bbox(spot_df)[3] + 0.1),
      ylim = c(st_bbox(spot_df)[2] - 0.1, st_bbox(spot_df)[4] + 0.1),
      expand = FALSE) +
    labs(fill = paste0(label_text, "\nClassification"), x = "Longitude", y = "Latitude") +
    guides(shape = guide_legend(override.aes = list(size = 0.2)), fill = guide_legend(ncol = 2)) +
    theme(panel.grid.major = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      panel.background = element_rect(fill = "white"),
      axis.text = element_text(size = 8),
      legend.title = element_text(size = 10),
      legend.position = "right",
      legend.key = element_rect(color = "white")) +
    geom_text(data = data.frame(label = label_text, x = -3.4, y = 44),
      aes(label = label, x = x, y = y), size = 5, fontface = "bold")
}

# -----------------------------------------------------------------------------
# Data preparation
# -----------------------------------------------------------------------------

# Harmonize labels and factor order across stages
ehsa_juvenile_90 <- format_ehsa(ehsa_juvenile_90, classification_order)
ehsa_adult_90    <- format_ehsa(ehsa_adult_90, classification_order)

ehsa_juvenile_90_proj <- st_transform(ehsa_juvenile_90, crs = med_albers_crs)

# Use juvenile grid to define the total spatial surface (same grid as adult)
total_surface_area_km2 <- ehsa_juvenile_90_proj %>%
  distinct(geometry, .keep_all = TRUE) %>%
  st_area() %>%
  sum() %>%
  units::set_units("km^2") %>%
  as.numeric()

# Strict hotspot / coldspot subsets (excludes "No Pattern Detected")
hs_juv <- ehsa_juvenile_90 %>% filter(classification %in% hotspot_classes)
hs_adu <- ehsa_adult_90 %>% filter(classification %in% hotspot_classes)
cs_juv <- ehsa_juvenile_90 %>% filter(classification %in% coldspot_classes)
cs_adu <- ehsa_adult_90 %>% filter(classification %in% coldspot_classes)

# "No Pattern Detected" subsets
npd_juv <- ehsa_juvenile_90 %>% filter(classification == "No Pattern Detected")
npd_adu <- ehsa_adult_90 %>% filter(classification == "No Pattern Detected")

# Prefix life stage for compact combined labels
hs_juv$classification <- shorten_classification(hs_juv$classification, "J", "Hotspot")
cs_juv$classification <- shorten_classification(cs_juv$classification, "J", "Coldspot")
hs_adu$classification <- shorten_classification(hs_adu$classification, "A", "Hotspot")
cs_adu$classification <- shorten_classification(cs_adu$classification, "A", "Coldspot")
npd_juv$classification <- shorten_classification(npd_juv$classification, "J", "Hotspot") 
npd_adu$classification <- shorten_classification(npd_adu$classification, "A", "Hotspot")

# -----------------------------------------------------------------------------
# Hotspot-only / Coldspot-only combined maps (Figures 5-6)
# -----------------------------------------------------------------------------

# Join juvenile/adult polygons into combined classes (hotspot x hotspot, coldspot x coldspot only)
# Warning: this step can be memory-intensive and take several minutes to complete
hs_df <- combine_stage_classes(hs_juv, hs_adu, npd_remove)
gc() # Clean up memory
cs_df <- combine_stage_classes(cs_juv, cs_adu, npd_remove)
gc()

# Build smooth palettes for each combined legend
hs_palette_base <- colorRampPalette(brewer.pal(9, "YlGn"))(13)
cs_palette_base <- colorRampPalette(brewer.pal(9, "BuPu"))(22)

hs_legend_col <- c(
  "J.Int - A.Int" = hs_palette_base[12],
  "J.Int - A.Per" = hs_palette_base[11],
  "J.Int - A.Spo" = hs_palette_base[10],
  "J.Per - A.Int" = hs_palette_base[9],
  "J.Per - A.Per" = hs_palette_base[8],
  "J.Per - A.Spo" = hs_palette_base[7],
  "J.Dim - A.Per" = hs_palette_base[6],
  "J.Dim - A.Spo" = hs_palette_base[5],
  "J.Spo - A.Int" = hs_palette_base[4],
  "J.Spo - A.Per" = hs_palette_base[3],
  "J.Spo - A.Spo" = hs_palette_base[2],
  #"J.Con - A.Per" = hs_palette_base[3],
  "J.New - A.Int" = hs_palette_base[1])
  #"J.New - A.Spo" = hs_palette_base[1])

cs_legend_col <- c(
  "J.Int - A.Int" = cs_palette_base[22],
  "J.Int - A.Per" = cs_palette_base[21],
  "J.Int - A.Spo" = cs_palette_base[20],
  "J.Int - A.Con" = cs_palette_base[19],
  "J.Int - A.New" = cs_palette_base[18],
  "J.Per - A.Int" = cs_palette_base[17],
  "J.Per - A.Per" = cs_palette_base[16],
  "J.Per - A.Spo" = cs_palette_base[15],
  "J.Per - A.Con" = cs_palette_base[14],
  "J.Per - A.New" = cs_palette_base[13],
  "J.Spo - A.Int" = cs_palette_base[12],
  "J.Spo - A.Per" = cs_palette_base[11],
  "J.Spo - A.Spo" = cs_palette_base[10],
  "J.Spo - A.Con" = cs_palette_base[9],
  "J.Spo - A.New" = cs_palette_base[8],
  "J.Con - A.Int" = cs_palette_base[7],
  "J.Con - A.Per" = cs_palette_base[6],
  "J.Con - A.Spo" = cs_palette_base[5],
  "J.Con - A.Con" = cs_palette_base[4],
  "J.New - A.Int" = cs_palette_base[3],
  "J.New - A.Per" = cs_palette_base[2],
  "J.New - A.Spo" = cs_palette_base[1])

# Apply factor levels to enforce plotting order
hs_df$Classification <- factor(hs_df$Classification, levels = hs_levels)
cs_df$Classification <- factor(cs_df$Classification, levels = cs_levels)

# Save combined hotspot and coldspot spatial datasets as a single sf object
combined_hs_cs <- bind_rows(
  hs_df %>% mutate(Spot = "Hotspots"),
  cs_df %>% mutate(Spot = "Coldspots"))
# saveRDS(combined_hs_cs, "./Outputs/EHSA/Combined_HS_CS.RDS")

## Figure 5: Hotspot and coldspot maps

gg.HS <- make_spot_map(hs_df, hs_legend_col, "Hotspots")
gg.CS <- make_spot_map(cs_df, cs_legend_col, "Coldspots")

gg.combined.spot <- ggplot(data = world) +
  geom_sf(data = hs_df, aes(fill = Classification), color = "white", lwd = 0.05, inherit.aes = TRUE) +
  scale_fill_manual("Hotspots Pattern", values = hs_legend_col) +
  guides(fill = guide_legend(ncol = 2)) +
  new_scale_fill() +
  geom_sf(data = cs_df, aes(fill = Classification), color = "white", lwd = 0.05, inherit.aes = TRUE) +
  scale_fill_manual("Coldspots Pattern", values = cs_legend_col) +
  geom_sf(fill = "grey90", color = "grey20") +
  coord_sf(xlim = c(st_bbox(hs_df)[1] - 0.5, st_bbox(hs_df)[3] + 0.1), ylim = c(st_bbox(hs_df)[2] - 0.1, st_bbox(hs_df)[4] + 0.4), expand = FALSE) +
  labs(fill = "Classification", x = "Longitude", y = "Latitude") +
  guides(fill = guide_legend(ncol = 2)) +
  theme(
    panel.grid.major = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.background = element_rect(fill = "white"),
    axis.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.position = "right",
    legend.key = element_rect(color = "white"),
    legend.key.size = unit(0.5, "cm"),
    plot.margin = unit(c(0, 0, 0, 0), "cm"))

# ggplot2::ggsave(gg.combined.spot, filename = "./Figures/Figure_5.png", width = 20, height = 15, units = "cm", dpi = 400)

## Figure 6: barplot of hotspot and coldspot areas by combined classification

hs_df_proj <- st_transform(hs_df, crs = med_albers_crs)
cs_df_proj <- st_transform(cs_df, crs = med_albers_crs)

areas_by_classification_hs <- summarize_class_area(hs_df_proj, "Hotspots", total_surface_area_km2)
areas_by_classification_cs <- summarize_class_area(cs_df_proj, "Coldspots", total_surface_area_km2)

tab_hs <- collapse_small_classes(areas_by_classification_hs)
tab_cs <- collapse_small_classes(areas_by_classification_cs)

tab_surf <- bind_rows(tab_hs, tab_cs) %>%
  mutate(Classification = as.factor(Classification))

cold_palette <- c(
  "J.Int - A.Int" = "#4D004B", "J.Int - A.Per" = "#60055D", "J.Int - A.Spo" = "#740B70",
  "J.Per - A.Int" = "#873C99", "J.Per - A.Per" = "#894DA2", "J.Per - A.Spo" = "#8A5DAA",
  "J.Spo - A.Int" = "#8C8DC2", "J.Spo - A.Per" = "#8F9DC9", "J.Spo - A.Spo" = "#96ABD1",
  "Others (< 1%)" = "grey")

hot_palette <- c(
  "J.Int - A.Int" = "#005C32", "J.Int - A.Per" = "#0C713B", "J.Int - A.Spo" = "#238443",
  "J.Per - A.Int" = "#379E54", "J.Per - A.Per" = "#53B466", "J.Per - A.Spo" = "#78C679",
  "J.Spo - A.Per" = "#EDF8B2", "J.Spo - A.Spo" = "#FAFDC8",
  "Others (< 1%)" = "grey")

gg.surf <- ggplot(tab_surf) +
  geom_col(data = filter(tab_surf, Spot == "Coldspots"), aes(x = Spot, y = pct_total, fill = fct_relevel(
        Classification,
        "Others (< 1%)", "J.Int - A.Int", "J.Int - A.Per",
        "J.Per - A.Int", "J.Per - A.Per", "J.Per - A.Spo",
        "J.Spo - A.Int", "J.Spo - A.Per", "J.Spo - A.Spo")), position = position_stack(reverse = FALSE)) +
  scale_fill_manual(values = cold_palette, name = "Combined Coldspots Pattern",
    breaks = c("J.Int - A.Int", "J.Int - A.Per",
      "J.Per - A.Int", "J.Per - A.Per", "J.Per - A.Spo",
      "J.Spo - A.Int", "J.Spo - A.Per", "J.Spo - A.Spo",
      "Others (< 1%)"), drop = FALSE) +
  guides(fill = guide_legend(ncol = 2)) +
  new_scale_fill() +
  geom_col(data = filter(tab_surf, Spot == "Hotspots"), aes(x = Spot, y = pct_total, fill = fct_relevel(
        Classification,
        "Others (< 1%)", "J.Int - A.Int", "J.Int - A.Per", "J.Int - A.Spo",
        "J.Per - A.Int", "J.Per - A.Per", "J.Per - A.Spo",
        "J.Spo - A.Per", "J.Spo - A.Spo")), position = position_stack(reverse = FALSE)) +
  scale_fill_manual(values = hot_palette, name = "Combined Hotspots Pattern",
    breaks = c("J.Int - A.Int", "J.Int - A.Per", "J.Int - A.Spo",
      "J.Per - A.Int", "J.Per - A.Per", "J.Per - A.Spo",
      "J.Spo - A.Per", "J.Spo - A.Spo",
      "Others (< 1%)"), drop = FALSE) +
  scale_x_discrete(limits = c("Coldspots", "Hotspots"), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 27),
    breaks = c(0, 10, 20, 25), labels = c("0%", "10%", "20%", "25%"), expand = c(0, 0)) +
  labs(x = "", y = "Proportion of total surface") +
  HMSC.theme +
  guides(fill = guide_legend(ncol = 2)) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

# ggplot2::ggsave(gg.surf, filename = "./Figures/Figure_6.png", width = 20, height = 15, units = "cm", dpi = 400)

# -----------------------------------------------------------------------------
# Full 9-way juvenile x adult combination analysis
# -----------------------------------------------------------------------------

# Warning: this step can be memory-intensive and take several minutes to complete !!!
# Need at least 35GB of RAM to run this step
combos <- list(
  hs_df %>% mutate(Combination = "Congruent hotspot (J.HS x A.HS)",
                    Pattern_status = "Both stages: pattern detected"),
  cs_df %>% mutate(Combination = "Congruent coldspot (J.CS x A.CS)",
                    Pattern_status = "Both stages: pattern detected"),
  join_combo(hs_juv, cs_adu, "Mismatch (J.HS x A.CS)", "Both stages: pattern detected"),
  join_combo(cs_juv, hs_adu, "Mismatch (J.CS x A.HS)", "Both stages: pattern detected"),
  join_combo(hs_juv, npd_adu, "One-sided no pattern (J.HS x A.Npd)", "One stage: no pattern"),
  join_combo(cs_juv, npd_adu, "One-sided no pattern (J.CS x A.Npd)", "One stage: no pattern"),
  join_combo(npd_juv, hs_adu, "One-sided no pattern (J.Npd x A.HS)", "One stage: no pattern"),
  join_combo(npd_juv, cs_adu, "One-sided no pattern (J.Npd x A.CS)", "One stage: no pattern"),
  join_combo(npd_juv, npd_adu, "Both stages: no pattern (J.Npd x A.Npd)", "Both stages: no pattern"))

gc() # Clean up memory

# Keep only the columns needed then stack all 9 combinations
all_combinations_df <- do.call(rbind, lapply(combos, function(df) df[, c("Classification", "Combination", "Pattern_status")]))

all_combinations_df$Classification <- as.factor(all_combinations_df$Classification)
all_combinations_df$Combination <- as.factor(all_combinations_df$Combination)
all_combinations_df$Pattern_status <- as.factor(all_combinations_df$Pattern_status)

# Project and compute polygon areas
all_combinations_proj <- st_transform(all_combinations_df, crs = med_albers_crs)

# Area per fine-grained classification, within each combination type
areas_by_classification <- all_combinations_proj %>%
  dplyr::mutate(area_km2 = st_area(geometry) %>% units::set_units("km^2") %>% as.numeric()) %>%
  st_drop_geometry() %>%
  dplyr::group_by(Combination, Pattern_status, Classification) %>%
  dplyr::summarise(total_area_km2 = sum(area_km2), n_polygons = n(), .groups = "drop") %>%
  dplyr::mutate(pct_total = round((total_area_km2 / total_surface_area_km2) * 100, 3),
    total_area_km2 = round(total_area_km2, 3)) %>%
  dplyr::arrange(Combination, desc(total_area_km2))

# Area collapsed to the 9 combination categories
summary_by_combination <- areas_by_classification %>%
  dplyr::group_by(Combination) %>%
  dplyr::summarise(total_area_km2 = sum(total_area_km2),
    pct_total = round(sum(pct_total), 3),
    n_polygons = sum(n_polygons),
    .groups = "drop")

# Same (split by pattern status)
summary_by_combination_pattern <- areas_by_classification %>%
  dplyr::group_by(Combination, Pattern_status) %>%
  dplyr::summarise(total_area_km2 = sum(total_area_km2),
    pct_total = round(sum(pct_total), 3),
    n_polygons = sum(n_polygons),
    .groups = "drop")

print(summary_by_combination)
print(summary_by_combination_pattern)

# saveRDS(areas_by_classification, "./Outputs/EHSA/Areas_by_classification_full.RDS")