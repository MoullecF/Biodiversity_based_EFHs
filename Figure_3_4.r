###############################################################################
# Title: EHSA Hotspot Maps and Area Summaries (Figures 3-4)
###############################################################################

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------

rm(list = ls())

source("./0-Load libraries.R")
source("./ggplot_theme.R")

# -----------------------------------------------------------------------------
# Load and order EHSA outputs
# -----------------------------------------------------------------------------

ehsa_juvenile_90 <- readRDS("./Outputs/EHSA/EHSA_juvenile.rds")
ehsa_adult_90 <- readRDS("./Outputs/EHSA/EHSA_adult.rds")

classification_order <- c(
  "intensifying hotspot",
  "persistent hotspot",
  "diminishing hotspot",
  "sporadic hotspot",
  "consecutive hotspot",
  "new hotspot",
  "no pattern detected",
  "new coldspot",
  "consecutive coldspot",
  "sporadic coldspot",
  "diminishing colspot",
  "persistent coldspot",
  "intensifying coldspot"
)

prepare_ehsa <- function(ehsa_df, order_levels) {
  # Standardize labels and enforce a consistent plotting order.
  ehsa_df %>%
    as.data.frame() %>%
    mutate(
      classification = str_to_title(classification),
      classification = factor(classification, levels = str_to_title(order_levels))
    ) %>%
    arrange(classification) %>%
    st_as_sf()
}

ehsa_juvenile_90 <- prepare_ehsa(ehsa_juvenile_90, classification_order)
ehsa_adult_90 <- prepare_ehsa(ehsa_adult_90, classification_order)

palette_ehsa <- c(
  "Intensifying Hotspot" = "#003c30",
  "Persistent Hotspot" = "#005F56",
  "Diminishing Hotspot" = "#23867E",
  "Sporadic Hotspot" = "#5AB2A8",
  "Consecutive Hotspot" = "#97D6CD",
  "New Hotspot" = "#CEEBE7",
  "No Pattern Detected" = "#F5F5F5",
  "New Coldspot" = "#F5EACB",
  "Consecutive Coldspot" = "#E6CE94",
  "Sporadic Coldspot" = "#CFA154",
  "Diminishing Coldspot" = "#AE7121",
  "Persistent Coldspot" = "#824B09",
  "Intensifying Coldspot" = "#543005"
)

# -----------------------------------------------------------------------------
# Figure 3: EHSA maps
# -----------------------------------------------------------------------------

plot_ehsa_map <- function(ehsa_sf, show_legend, stage_label) {
  # Render classification polygons in the original spatial extent
  ggplot(data = world) +
    geom_sf(
      data = ehsa_sf,
      aes(fill = classification),
      color = "white",
      lwd = 0.05,
      show.legend = show_legend
    ) +
    scale_fill_manual(values = palette_ehsa) +
    geom_sf(fill = "grey90", color = "grey20") +
    annotate("text", x = -3.4, y = 44, label = stage_label, fontface = "bold") +
    coord_sf(
      xlim = c(st_bbox(ehsa_sf)[1] - 0.1, st_bbox(ehsa_sf)[3] + 0.1),
      ylim = c(st_bbox(ehsa_sf)[2] - 0.1, st_bbox(ehsa_sf)[4] + 0.1),
      expand = FALSE
    ) +
    labs(fill = "Pattern", x = "Longitude", y = "Latitude") +
    guides(shape = guide_legend(override.aes = list(size = 0.2))) +
    theme(
      panel.grid.major = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      panel.background = element_rect(fill = "white"),
      axis.text = element_text(size = 8),
      legend.title = element_text(size = 10),
      legend.position = "right",
      legend.key = element_rect(color = "white")
    )
}

ehsa_juv_90 <- plot_ehsa_map(ehsa_juvenile_90, show_legend = TRUE, stage_label = "Juveniles")

ehsa_adu_90 <- plot_ehsa_map(ehsa_adult_90, show_legend = FALSE, stage_label = "Adults")

map_ehsa <- ehsa_juv_90 / ehsa_adu_90 +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "right")

# Save
# ggplot2::ggsave(map_ehsa, filename = "./Figures/Figure_3.png", width = 20, height = 15, units = "cm", dpi = 400)

# -----------------------------------------------------------------------------
# Figure 4: Surface areas by classification
# -----------------------------------------------------------------------------

med_albers <- "+proj=aea +lat_1=30 +lat_2=45 +lat_0=37.5 +lon_0=15 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
ehsa_juvenile_90_proj <- st_transform(ehsa_juvenile_90, crs = med_albers)
ehsa_adult_90_proj <- st_transform(ehsa_adult_90, crs = med_albers)

surface_area_km2 <- ehsa_juvenile_90_proj %>%
  st_area() %>%
  sum() %>%
  units::set_units("km^2") %>%
  as.numeric()

summarize_area_by_classification <- function(ehsa_proj, onto_label, surface_area_km2) {
  # Collapse to total area per classification and aggregate small classes
  ehsa_proj %>%
    mutate(
      area_km2 = st_area(geometry) %>%
        units::set_units("km^2") %>%
        as.numeric(),
      onto = onto_label
    ) %>%
    st_drop_geometry() %>%
    group_by(classification, onto) %>%
    summarise(total_area_km2 = sum(area_km2), n_polygons = n(), .groups = "drop") %>%
    mutate(
      pct_total = (total_area_km2 / surface_area_km2) * 100,
      total_area_km2 = round(total_area_km2, 3),
      pct_total = round(pct_total, 1)
    ) %>%
    mutate(
      classification_group = ifelse(pct_total < 1, "Others (< 1%)", as.character(classification))
    ) %>%
    group_by(classification_group, onto) %>%
    summarize(
      total_area_km2 = sum(total_area_km2),
      n_polygons = sum(n_polygons),
      pct_total = sum(pct_total),
      .groups = "drop"
    ) %>%
    rename(classification = classification_group) %>%
    mutate(classification = as.factor(classification)) %>%
    arrange(desc(pct_total))
}

areas_by_classification_juv <- summarize_area_by_classification(
  ehsa_juvenile_90_proj,
  onto_label = "Juveniles",
  surface_area_km2 = surface_area_km2
)

areas_by_classification_adu <- summarize_area_by_classification(
  ehsa_adult_90_proj,
  onto_label = "Adults",
  surface_area_km2 = surface_area_km2
)

tab_surf <- bind_rows(areas_by_classification_juv, areas_by_classification_adu) %>%
  mutate(onto = factor(onto, levels = c("Juveniles", "Adults"))) %>%
  mutate(
    classification = factor(
      classification,
      levels = c(
        "Intensifying Hotspot",
        "Persistent Hotspot",
        "Diminishing Hotspot",
        "Sporadic Hotspot",
        "Consecutive Hotspot",
        "New Hotspot",
        "No Pattern Detected",
        "New Coldspot",
        "Consecutive Coldspot",
        "Sporadic Coldspot",
        "Diminishing Coldspot",
        "Persistent Coldspot",
        "Intensifying Coldspot",
        "Others (< 1%)"
      )
    )
  )

hotspot_summary <- tab_surf %>%
  filter(
    classification %in% c(
      "Persistent Hotspot",
      "Sporadic Hotspot",
      "Intensifying Hotspot",
      "Diminishing Hotspot",
      "New Hotspot",
      "Consecutive Hotspot"
    )
  ) %>%
  group_by(onto) %>%
  summarise(
    total_hotspot_pct = sum(pct_total),
    total_area_km2 = sum(total_area_km2),
    total_polygons = sum(n_polygons)
  ) %>%
  mutate(total_hotspot_pct = round(total_hotspot_pct, 1))
print(hotspot_summary)

coldspot_summary <- tab_surf %>%
  filter(
    classification %in% c(
      "Persistent Coldspot",
      "Sporadic Coldspot",
      "Intensifying Coldspot",
      "Diminishing Coldspot",
      "New Coldspot",
      "Consecutive Coldspot"
    )
  ) %>%
  group_by(onto) %>%
  summarise(
    total_Coldspot_pct = sum(pct_total),
    total_area_km2 = sum(total_area_km2),
    total_polygons = sum(n_polygons)
  ) %>%
  mutate(total_Coldspot_pct = round(total_Coldspot_pct, 1))
print(coldspot_summary)

custom_ylabels <- c(
  "<span style='font-size:12pt'>0%</span>",
  "<span style='font-size:12pt'>25%</span>",
  "<span style='font-size:10pt'>46%</span>",
  "<span style='font-size:10pt'>35.2%</span>",
  "<span style='font-size:12pt'>50%</span>",
  "<span style='font-size:12pt'>75%</span>",
  "<span style='font-size:12pt'>100%</span>"
)

palette_ehsa_with_other <- c(palette_ehsa, "Others (< 1%)" = "grey")

gg_surf_spot <- ggplot(tab_surf, aes(x = onto, y = pct_total, fill = classification)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = palette_ehsa_with_other, guide = guide_legend(reverse = TRUE)) +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 25, 46, 35.2, 50, 75, 100), labels = custom_ylabels, expand = c(0, 0)) +
  scale_x_discrete(limits = c("Juveniles", "Adults"), expand = c(0, 0)) +
  labs(x = "", y = "Proportion of total surface", fill = "Pattern") +
  geom_segment(aes(x = 0.55, xend = 1.45, y = 46, yend = 46),
    linetype = "dashed", color = "black", linewidth = 0.7) +
  geom_segment(aes(x = 0.55, xend = 2.45, y = 35.2, yend = 35.2),
    linetype = "dashed", color = "black", linewidth = 0.7) +
  HMSC.theme +
  theme(plot.margin = unit(c(1, 1, 1, 1), units = , "cm"), axis.text.y = element_markdown())

# Save
# ggplot2::ggsave(gg_surf_spot, filename = "./Figures/Figure_4.png", width = 20, height = 15, units = "cm", dpi = 400)
