###############################################################################
#
# Extrapolation in multivariate environmental space with dsmextra
# Based on: https://densitymodelling.github.io/dsmextra/articles/dsmextra.html
#
###############################################################################

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------

rm(list = ls())

source("Scripts_R/0-Load libraries.R")
# remotes::install_github("densitymodelling/dsmextra")
library(dsmextra)


# -----------------------------------------------------------------------------
# 0. Settings
# -----------------------------------------------------------------------------

model_file <- "Models/model_mpa_thin_4000_samples_250_chains_4.RData"
grid_file  <- "Inputs_HMSC/Hindcast_Env_grid_rect0.05_WMED_0_1000.RData"
out_dir    <- "Outputs/Extrapolation/"

years <- 1999:2021

# Covariates, in the names used in model$XData
covariates <- c("Depth", "SST", "SBT", "SSS", "Chla", "Fishing", "Gravity")

# Same covariates in grid_env (names = XData names, values = grid_env columns)
grid_names <- c(Depth = "Depth", SST = "SST", SBT = "botTemp", SSS = "so",
                Chla = "chl", Fishing = "FPI_tot", Gravity = "Gravity")

# XData stores depth with the opposite sign to grid_env
flip_xdata_depth <- TRUE

my_crs <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs")
nearby_radius <- 1 # default (Mannocci et al. 2018; Virgili et al. 2017)
snap_if_irregular <- TRUE # FALSE = stop instead
snap_digits <- 4
cell_size <- 0.05

# Depth mask. Only prediction cells inside the mask are kept as target data
# depth_mask <- FALSE to run on the full grid
apply_depth_mask <- TRUE
mask_file <- "Inputs_HMSC/depth_mask.RDS"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "Figures"), recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. Reference and target data
# -----------------------------------------------------------------------------

model <- get(load(model_file))
stopifnot(all(covariates %in% names(model$XData)))
segs <- as.data.frame(model$XData)[, covariates]
if (flip_xdata_depth) segs$Depth <- -segs$Depth
rm(model)

load(grid_file)
stopifnot(all(c("X", "Y", "YEAR", grid_names) %in% names(grid_env)))

# Drop potential duplicate rows per (X, Y, YEAR) to avoid redundant computation
n_before <- nrow(grid_env)
grid_env$kx <- round(grid_env$X, 6)
grid_env$ky <- round(grid_env$Y, 6)
n_same_cov <- nrow(dplyr::distinct(cbind(grid_env[, c("kx", "ky", "YEAR")], round(grid_env[, unname(grid_names)], 8))))
grid_env <- dplyr::distinct(grid_env, kx, ky, YEAR, .keep_all = TRUE)
grid_env$kx <- NULL; grid_env$ky <- NULL
message("grid_env: ", n_before, " rows -> ", nrow(grid_env), " after dropping duplicates per (X, Y, YEAR)")
if (n_same_cov != nrow(grid_env)) {
  warning(n_same_cov - nrow(grid_env), " dropped rows differed in >= 1 covariate from the row kept. ",
          "Check what the repeated rows are.")
}
stopifnot(all(c("X", "Y", "YEAR", grid_names) %in% names(grid_env)))
predgrid_all <- data.frame(x = grid_env$X, y = grid_env$Y, YEAR = grid_env$YEAR)
for (v in covariates) predgrid_all[[v]] <- grid_env[[grid_names[[v]]]]
predgrid_all <- predgrid_all[predgrid_all$YEAR %in% years, ]
rm(grid_env)

# -----------------------------------------------------------------------------
# 1b. Depth mask on the target data
# -----------------------------------------------------------------------------
# Keep only the prediction cells that fall in a non-NA cell of depth_mask

if (apply_depth_mask) {
  depth_mask <- readRDS(mask_file)
  mask_values <- raster::getValues(depth_mask)
  n_mask_cells <- sum(!is.na(mask_values))

  # Mask cell containing each prediction point (NA if outside the mask extent)
  mask_cell <- raster::cellFromXY(depth_mask, cbind(predgrid_all$x, predgrid_all$y))

  # Alignment check: points must sit near the centres of the mask cells. A point
  # on a cell edge (offset = 0.5 cell) could be assigned to the wrong cell.
  centre <- raster::xyFromCell(depth_mask, mask_cell)
  offset <- pmax(abs(predgrid_all$x - centre[, 1]) / raster::res(depth_mask)[1],
                 abs(predgrid_all$y - centre[, 2]) / raster::res(depth_mask)[2])
  if (max(offset, na.rm = TRUE) > 0.1) {
    stop("Prediction points are not centred on the depth mask cells (max offset = ",
         round(max(offset, na.rm = TRUE), 3), " cell): the two grids are not aligned.")
  }

  in_mask <- !is.na(mask_cell) & !is.na(mask_values[mask_cell])
  n_cells_before <- length(unique(mask_cell))
  n_rows_before  <- nrow(predgrid_all)
  predgrid_all <- predgrid_all[in_mask, ]
  n_cells_after <- length(unique(mask_cell[in_mask]))

  message("Depth mask: ", n_cells_before, " cells -> ", n_cells_after, " kept (",
          n_cells_before - n_cells_after, " removed, ",
          round(100 * (n_cells_before - n_cells_after) / n_cells_before, 1), " %); ",
          n_rows_before, " rows -> ", nrow(predgrid_all))
  if (n_cells_after != n_mask_cells) {
    warning("The depth mask has ", n_mask_cells, " cells but ", n_cells_after,
            " are present in the prediction grid.")
  }
  rm(depth_mask, mask_values, mask_cell, centre, offset, in_mask)
}

message("Reference: ", nrow(segs), " hauls (rows with NA are dropped by dsmextra: ",
        sum(!complete.cases(segs)), ")")
message("Target: ", nrow(predgrid_all), " rows, ",
        nrow(distinct(predgrid_all[, c("x", "y")])), " distinct x,y")
message("Reference vs grid ranges (check units and Depth sign):")
print(data.frame(covariate = covariates,
                 ref_min  = sapply(segs[covariates], min, na.rm = TRUE),
                 ref_max  = sapply(segs[covariates], max, na.rm = TRUE),
                 grid_min = sapply(predgrid_all[covariates], min, na.rm = TRUE),
                 grid_max = sapply(predgrid_all[covariates], max, na.rm = TRUE),
                 row.names = NULL), digits = 4)

# dsmextra needs a regular grid
grid_snapped <- FALSE
call_dsmextra <- function(fun, pg, ...) {
  res <- tryCatch(fun(prediction.grid = pg, ..., verbose = FALSE), error = function(e) e)
  if (inherits(res, "error")) {
    if (!grepl("regularly spaced", conditionMessage(res)) || !snap_if_irregular) stop(res)
    message("   raw grid not recognised as regular: x, y snapped to ", snap_digits,
            " decimals inside dsmextra only")
    grid_snapped <<- TRUE
    pg$x <- round(pg$x, snap_digits)
    pg$y <- round(pg$y, snap_digits)
    res <- tryCatch(fun(prediction.grid = pg, ..., verbose = FALSE),
                    error = function(e)
                      stop("The grid is not regular even after snapping: ", conditionMessage(e)))
  }
  res
}

# -----------------------------------------------------------------------------
# 2. Run year by year
# -----------------------------------------------------------------------------

cells_list  <- list() # ExDet, class and MIC per row
nearby_list <- list() # % nearby per cell
by_year     <- list()
summary_log <- character()

for (yr in years) {
  message("\n=== ", yr, " ===")
  predgrid <- predgrid_all[predgrid_all$YEAR == yr, c("x", "y", covariates)]
  predgrid <- predgrid[complete.cases(predgrid), ] # keeps rows aligned with dsmextra output
  grid_snapped <- FALSE

  # ExDet: univariate and combinatorial extrapolation
  ex <- call_dsmextra(compute_extrapolation, pg = predgrid, samples = segs, covariate.names = covariates, coordinate.system = my_crs)
  flag_ex <- grid_snapped

  # % of cells by extrapolation type and MIC tables
  summ <- tryCatch(capture.output(summary(ex)), error = function(e) paste("summary() failed:", conditionMessage(e)))
  cat(summ, sep = "\n")
  summary_log <- c(summary_log, paste("=====", yr, "====="), summ, "")

  d <- ex$data$all
  stopifnot(nrow(d) == nrow(predgrid))
  d$x <- predgrid$x; d$y <- predgrid$y # raw coordinates (same row order)
  d$class <- ifelse(d$ExDet < 0, "univariate",
             ifelse(d$ExDet > 1, "combinatorial", "analogue"))
  # mic_* are positions in covariates
  d$mic_label <- ifelse(d$class == "univariate",    covariates[d$mic_univariate],
                 ifelse(d$class == "combinatorial", covariates[d$mic_combinatorial],
                        NA_character_))
  d$YEAR <- yr
  cells_list[[as.character(yr)]] <- d[, c("YEAR", "x", "y", "ExDet", "class", "mic_label")]

  by_year[[as.character(yr)]] <- data.frame(
    YEAR = yr,
    n_rows = nrow(d),
    n_cells_xy = nrow(dplyr::distinct(d[, c("x", "y")])),
    grid_snapped = flag_ex,
    n_univariate    = sum(d$class == "univariate"),
    n_combinatorial = sum(d$class == "combinatorial"),
    n_analogue      = sum(d$class == "analogue"))

  # % of reference data nearby
  nb <- call_dsmextra(compute_nearby, pg = predgrid,
                      samples = segs, covariate.names = covariates,
                      coordinate.system = my_crs, nearby = nearby_radius)
  # compute_nearby returns a raster (one value per cell)
  nb_df <- raster::as.data.frame(nb$raster, xy = TRUE, na.rm = TRUE)
  names(nb_df) <- c("xr", "yr", "perc_nearby")
  nb_df$kx <- round(nb_df$xr, 4); nb_df$ky <- round(nb_df$yr, 4)
  raw_xy <- predgrid %>% mutate(kx = round(x, 4), ky = round(y, 4)) %>%
    distinct(kx, ky, .keep_all = TRUE) %>% dplyr::select(kx, ky, x, y)
  nb_df <- inner_join(nb_df, raw_xy, by = c("kx", "ky")) %>% dplyr::select(x, y, perc_nearby)
  nearby_list[[as.character(yr)]] <- cbind(YEAR = yr, nb_df)

  rm(ex, nb, d, nb_df)
}

cells  <- bind_rows(cells_list)
nearby <- bind_rows(nearby_list)

# -----------------------------------------------------------------------------
# 3. Tables
# -----------------------------------------------------------------------------

# % of cells by year: univariate, combinatorial, analogue
tab_year <- bind_rows(by_year) %>%
  mutate(pct_univariate    = 100 * n_univariate / n_rows,
         pct_combinatorial = 100 * n_combinatorial / n_rows,
         pct_analogue      = 100 * n_analogue / n_rows)

write.csv(tab_year, file.path(out_dir, "extrapolation_by_year.csv"), row.names = FALSE)
print(as.data.frame(tab_year), digits = 4)

# Most influential covariate (MIC) by year and extrapolation type
tab_mic <- cells %>%
  filter(class != "analogue") %>%
  dplyr::count(YEAR, type = class, covariate = mic_label, name = "n") %>%
  group_by(YEAR, type) %>%
  mutate(pct_of_type = 100 * n / sum(n)) %>%
  ungroup()

write.csv(tab_mic, file.path(out_dir, "mic_by_year.csv"), row.names = FALSE)

tab_mic_all <- cells %>%
  filter(class != "analogue") %>%
  dplyr::count(type = class, covariate = mic_label, name = "n") %>%
  group_by(type) %>%
  mutate(pct_of_type = 100 * n / sum(n)) %>%
  ungroup() %>%
  arrange(type, desc(n))

write.csv(tab_mic_all, file.path(out_dir, "mic_all_years.csv"), row.names = FALSE)
print(as.data.frame(tab_mic_all), digits = 4)

writeLines(summary_log, file.path(out_dir, "summary_dsmextra_by_year.txt"))
saveRDS(cells,  file.path(out_dir, "extrapolation_cells_by_year.rds"))
saveRDS(nearby, file.path(out_dir, "nearby_cells_by_year.rds"))

# Extrapolation risk across years, one value per cell
# ExDet is averaged across all years regardless of each year's class, then
# collapsed to a single-sided risk score: 0 for analogue cells, increasing
# with the severity of univariate extrapolation. No combinatorial extrapolation
# occurs in this dataset (checked above), so the score only needs the
# univariate side; if combinatorial cells appear with a different dataset,
# add back a max(0, mean_ExDet - 1) term.
stopifnot(all(cells$class != "combinatorial"))
cells_risk <- cells %>%
  dplyr::group_by(x, y) %>%
  dplyr::summarise(mean_ExDet = mean(ExDet), n_years_total = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(risk = pmax(0, -mean_ExDet))

cells_risk_sf <- sf::st_as_sf(cells_risk, coords = c("x", "y"), crs = my_crs, remove = FALSE)

saveRDS(cells_risk_sf, file.path(out_dir, "extrapolation_risk.rds"))

# -----------------------------------------------------------------------------
# 4. Maps: one multipanel figure per map type, one panel per year
# -----------------------------------------------------------------------------

xlim <- range(predgrid_all$x) + c(-0.3, 0.3)
ylim <- range(predgrid_all$y) + c(-0.3, 0.3)

# Land
world <- tryCatch(rnaturalearth::ne_countries(scale = "large", returnclass = "sf"), error = function(e) NULL)
if (!is.null(world)) {
  world <- sf::st_sf(geometry = sf::st_geometry(world))
  old_s2 <- sf::sf_use_s2()
  suppressMessages(sf::sf_use_s2(FALSE))
  world <- tryCatch(suppressWarnings(sf::st_crop(world, xmin = xlim[1], xmax = xlim[2],
                                                 ymin = ylim[1], ymax = ylim[2])), error = function(e) world)
  suppressMessages(sf::sf_use_s2(old_s2))
}

n_col <- 5

map_layers <- function() {
  list(if (!is.null(world)) geom_sf(data = world, fill = "grey80", colour = "grey30", linewidth = 0.1),
       coord_sf(xlim = xlim, ylim = ylim, expand = FALSE),
       facet_wrap(~YEAR, ncol = n_col),
       labs(x = NULL, y = NULL),
       theme(panel.grid.major = element_blank(),
             panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4),
             panel.background = element_rect(fill = "white"),
             axis.text = element_text(size = 6),
             strip.text = element_text(size = 9),
             legend.title = element_text(size = 9),
             legend.text = element_text(size = 8)))
}

# One tile per cell and year and same colour limits in every panel
cells_d  <- distinct(cells, YEAR, x, y, .keep_all = TRUE)
lim_uni  <- c(min(cells_d$ExDet[cells_d$class == "univariate"], -1), 0)
lim_comb <- c(1, max(cells_d$ExDet[cells_d$class == "combinatorial"], 2))
mic_cols <- c(Depth = "#0072B2", SST = "#D55E00", SBT = "#E69F00", SSS = "#009E73",
              Chla = "#56B4E9", Fishing = "#CC79A7", Gravity = "#000000")

# Analogue, univariate and combinatorial extrapolation on the same map
p1 <- ggplot() +
  geom_tile(data = filter(cells_d, class == "analogue"), aes(x, y, fill = ExDet), width = cell_size, height = cell_size) +
  scale_fill_gradient(name = "Analogue\n(ExDet 0-1)", low = "#e5f5e0", high = "#00441b", limits = c(0, 1)) +
  new_scale_fill() +
  geom_tile(data = filter(cells_d, class == "univariate"), aes(x, y, fill = ExDet), width = cell_size, height = cell_size) +
  scale_fill_gradient(name = "Univariate\n(ExDet < 0)", low = "#08306b", high = "#c6dbef", limits = lim_uni) +
  new_scale_fill() +
  geom_tile(data = filter(cells_d, class == "combinatorial"), aes(x, y, fill = ExDet), width = cell_size, height = cell_size) +
  scale_fill_gradient(name = "Combinatorial\n(ExDet > 1)", low = "#fcbba1", high = "#99000d", limits = lim_comb) +
  map_layers()

ggsave(file.path(out_dir, "Figures", "extrapolation_all_years.png"), p1, width = 38, height = 24, units = "cm", dpi = 300)

# Most influential covariate (univariate and combinatorial cells together and analogue cells in grey)
p2 <- ggplot() +
  geom_tile(data = filter(cells_d, class == "analogue"), aes(x, y), fill = "grey90", width = cell_size, height = cell_size) +
  geom_tile(data = filter(cells_d, class != "analogue"), aes(x, y, fill = mic_label), width = cell_size, height = cell_size) +
  scale_fill_manual(name = "Most influential\ncovariate", values = mic_cols, limits = covariates, drop = FALSE) +
  map_layers()

ggsave(file.path(out_dir, "Figures", "mic_all_years.png"), p2, width = 38, height = 24, units = "cm", dpi = 300)

# 4.3 Percentage of reference data nearby
p3 <- ggplot() +
  geom_tile(data = nearby, aes(x, y, fill = perc_nearby), width = cell_size, height = cell_size) +
  scale_fill_viridis_c(name = "% reference\ndata nearby", limits = c(0, max(nearby$perc_nearby))) +
  map_layers()

ggsave(file.path(out_dir, "Figures", "nearby_all_years.png"), p3, width = 38, height = 24, units = "cm", dpi = 300)

# 4.4 Extrapolation risk across years
p4 <- ggplot() +
  geom_tile(data = cells_risk, aes(x, y, fill = risk),
            width = cell_size, height = cell_size) +
  scale_fill_viridis_c(name = "Extrapolation\nrisk\n(mean ExDet, 0 = analogue)") +
  { if (!is.null(world)) geom_sf(data = world, fill = "grey80", colour = "grey30",
                                 linewidth = 0.1) } +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  labs(x = NULL, y = NULL, title = "Extrapolation risk across years") +
  theme(panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4),
        panel.background = element_rect(fill = "white"))

ggsave(file.path(out_dir, "Figures", "extrapolation_risk_all_years.png"), p4, width = 20, height = 16, units = "cm", dpi = 300)

# -----------------------------------------------------------------------------
# 5. Reproducibility
# -----------------------------------------------------------------------------

writeLines(c(paste("dsmextra", as.character(packageVersion("dsmextra"))), capture.output(sessionInfo())), file.path(out_dir, "sessionInfo.txt"))