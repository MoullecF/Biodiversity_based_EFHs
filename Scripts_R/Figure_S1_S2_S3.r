###############################################################################
# Title: Survey Coverage and Effort (Figures S1–S3)
###############################################################################

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------

source("Scripts_R/0-Load libraries.r")
source("Scripts_R/ggplot_theme.r")

load(file = "Inputs_HMSC/allData.RData")

# -----------------------------------------------------------------------------
# Data preparation
# -----------------------------------------------------------------------------

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

grid <- sf::st_read("Inputs_HMSC/grid/grid_0.05_rectangles.shp")
grid <- grid[!grid$bottom >43 | !grid$right<0,]

S <- S %>%
  mutate(gsa = factor(sub("^([0-9]+).*", "\\1", id),
                      levels = c("1", "2", "5", "6", "7", "8", "9", "10", "11")))

gsa_labels <- c(
  "1" = "(1) Northern\nAlboran Sea",
  "2" = "(2) Alboran\nIsland",
  "5" = "(5) Balearic\nIslands",
  "6" = "(6) Northern\nSpain",
  "7" = "(7) Gulf of\nLion",
  "8" = "(8) Corsica",
  "9" = "(9) Ligurian Sea and\nNorthern Tyrrhenian Sea",
  "10" = "(10) Southern and\ncentral Tyrrhenian Sea",
  "11" = "(11.1 & 11.2) Western and\nEastern Sardinia"
)

gsa_palette <- setNames(brewer.pal(length(gsa_labels), "Set1"), names(gsa_labels))

# -----------------------------------------------------------------------------
# Figure S1: Western Mediterranean hauls (1999–2021)
# -----------------------------------------------------------------------------

gg_map_medits_hauls <- ggplot(world) +
  geom_sf(data = grid, fill = "white", color = "gray90") +
  geom_point(
    data = S,
    aes(x = Long, y = Lat, colour = gsa),
    position = position_jitter(width = 0.02, height = 0.02),
    size = 0.15
  ) +
  scale_colour_manual(values = gsa_palette, breaks = names(gsa_labels), labels = gsa_labels) +
  geom_sf(fill = "grey80", color = "grey20") +
  geom_text(label = "Western\nMediterranean Sea", x = 5.7, y = 39, size = 3) +
  coord_sf(
    xlim = range(S$Long) + c(-0.5, 1.5),
    ylim = range(S$Lat) + c(-1, 1.5),
    expand = FALSE
  ) +
  labs(x = "Longitude", y = "Latitude") +
  guides(colour = guide_legend(title = "Geographical\nSub-Areas", override.aes = list(size = 4))) +
  theme(
    panel.grid.major = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.background = element_rect(fill = "white"),
    plot.margin = unit(c(11, 11, 5.5, 5.5), "pt"),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 12),
    legend.position = "bottom"
  )

# ggplot2::ggsave(gg_map_medits_hauls, filename = "Figures/Figure_S1.png", width = 20, height = 15, units = "cm", dpi = 400)

# -----------------------------------------------------------------------------
# Figure S2: Number of hauls per year
# -----------------------------------------------------------------------------

haul_counts_year <- S %>%
  distinct(Year, id) %>%
  count(Year, name = "haul_count")

mean_haul_count <- mean(haul_counts_year$haul_count, na.rm = TRUE)

gg_number_hauls_year <- ggplot(haul_counts_year, aes(x = Year, y = haul_count)) +
  geom_col(fill = "grey", color = "grey35", alpha = 0.5) +
  geom_hline(yintercept = mean_haul_count, linetype = "dashed") +
  geom_text(
    size = 3,
    x = 2006,
    y = mean_haul_count + 55,
    label = paste0("Average number per year: ~", round(mean_haul_count), " hauls")
  ) +
  theme_bw() +
  labs(x = "Year", y = "Number of hauls")

# ggplot2::ggsave(gg_number_hauls_year, filename = "Figures/Figure_S2.png", width = 20, height = 15, units = "cm", dpi = 400)

# -----------------------------------------------------------------------------
# Figure S3: Haul counts per GSA and year
# -----------------------------------------------------------------------------

haul_counts_gsa <- S %>%
  distinct(Year, gsa, id) %>%
  count(Year, gsa, name = "haul_count")

y_labels <- c(
  "1" = "Northern\nAlboran Sea",
  "2" = "Alboran\nIsland",
  "5" = "Balearis\nIsland",
  "6" = "Northern\nSpain",
  "7" = "Gulf of\nLion",
  "8" = "Corsica",
  "9" = "Ligurian Sea and\nNorthern Tyrrhenian Sea",
  "10" = "Southern and\ncentral Tyrrhenian Sea",
  "11" = "Western and\nEastern Sardinia"
)

gg_number_hauls_gsa <- ggplot(haul_counts_gsa, aes(Year, gsa, fill = haul_count)) +
  geom_tile(size = 0.1, color = "grey80") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(
    breaks = names(y_labels),
    labels = y_labels,
    expand = c(0, 0)
  ) +
  scale_fill_viridis(discrete = FALSE, na.value = "white") +
  labs(x = "Years", y = "Geographical Sub-Areas (GSAs)", fill = "Number\nof hauls") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# ggplot2::ggsave(gg_number_hauls_gsa, filename = "Figures/Figure_S3.png", width = 20, height = 15, units = "cm", dpi = 400)

