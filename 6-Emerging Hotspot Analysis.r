
###############################################################################
# Title: Emerging Hotspot Analysis for HMSC Predictions
###############################################################################

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------

rm(list = ls())

source("./0-Load libraries.R")

# Quick fix for sfdep shuffle_nbs performance issue (to speed up!)
fast_shuffle <- function(i, n, card) {
  s <- sample.int(n - 1L, card, useHash = TRUE)
  s + (s >= i)
}
assignInNamespace("shuffle_nbs", fast_shuffle, ns = "sfdep")

# -----------------------------------------------------------------------------
# Inputs
# -----------------------------------------------------------------------------

# Number of permutations for the simulated p-value of Gi*
nsim <- 999

# Seed for reproducible simulated p-values
seed <- 123

# Load depth mask
depth_mask <- readRDS("./Inputs_HMSC/depth_mask.RDS")

# Load raster stack of hurdle-model outputs
raster_stack <- get(load("./Outputs/Spatio_temporal_prediction/r_stack_Hurdle_0.05_0_1000_19992021.Rdata"))

# Set to c("adult") for a single run, or c("juvenile", "adult") for both
stages_to_run <- c("juvenile", "adult")

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

step_msg <- function(...) message(format(Sys.time(), "%H:%M:%S"), " | ", ...)

build_yearly_sf <- function(stage_stack, depth_mask) {
  # Sum abundance across all layers of each year (one layer per year)
  years <- sub(".*_(\\d{4})$", "\\1", names(stage_stack))
  unique_years <- unique(years)
  yearly <- stackApply(stage_stack,
                       indices = as.integer(factor(years, levels = unique_years)),
                       fun = sum, na.rm = FALSE)
  names(yearly) <- paste0("tot_", unique_years)

  # Apply depth mask (to remove potential extrapolation beyond depth range)
  yearly <- mask(yearly, depth_mask)

  # Convert to polygons: one row per pixel, one column per year
  sf_raster <- st_as_sf(rasterToPolygons(yearly, dissolve = FALSE))

  # Pixels with no neighbours (to be removed, to avoid issues in EHSA)
  # The geometry is the same every year so this is computed only once
  list_nb <- suppressWarnings(poly2nb(sf_raster, queen = TRUE))
  empty_nb <- which(card(list_nb) == 0)
  step_msg(nrow(sf_raster), " pixels, ", length(empty_nb), " without neighbours removed")

  # For each year: zeros to NA, cap upper outliers (MAD), drop isolated pixels
  summed_list <- lapply(unique_years, function(year) {
    summed_ab <- sf_raster[, paste0("tot_", year)]
    names(summed_ab)[1] <- "tot_ab"

    summed_ab$tot_ab[summed_ab$tot_ab == 0] <- NA
    outliers_raster <- outliers_mad(summed_ab$tot_ab)
    upper <- outliers_raster$limits[2]
    summed_ab$tot_ab[which(summed_ab$tot_ab > upper)] <- upper

    summed_ab$Year <- as.integer(year)

    if (length(empty_nb) > 0) {
      summed_ab <- summed_ab[-empty_nb, ]
    }

    summed_ab$pixel_id <- seq_len(nrow(summed_ab))
    summed_ab
  })

  do.call(rbind, summed_list)
}

add_geometry <- function(ehsa, pixels_sf) {
  # Attach the pixel polygons to the EHSA table (location = pixel_id)
  idx <- match(ehsa$location, pixels_sf$pixel_id)
  stopifnot(!anyNA(idx))
  out <- st_sf(ehsa, geometry = st_geometry(pixels_sf)[idx])
  # Keep the Gi* table (include_gi = TRUE), stored as an attribute
  attr(out, "gi_star") <- attr(ehsa, "gi_star")
  out
}

run_ehsa <- function(stage_name, raster_stack, depth_mask, nsim = 199, seed = NULL) {
  # Run EHSA for a single life stage and return results
  step_msg("[", stage_name, "] building yearly sums")
  stage_stack <- raster::subset(raster_stack, grep(stage_name, names(raster_stack), value = TRUE))
  summed_sf <- build_yearly_sf(stage_stack, depth_mask)

  # emerging_hotspot_analysis() does not handle missing values
  n_na <- sum(is.na(summed_sf$tot_ab))
  if (n_na > 0) {
    stop("[", stage_name, "] ", n_na, " missing values in tot_ab (zeros set to NA): ",
         "fill or remove them before running the EHSA")
  }

  # Convert sf object to spacetime object
  step_msg("[", stage_name, "] building spacetime cube")
  spt <- as_spacetime(summed_sf, .loc_col = "pixel_id", .time_col = "Year")
  stopifnot(is_spacetime_cube(spt))

  step_msg("[", stage_name, "] running EHSA (nsim = ", nsim, ")")
  if (!is.null(seed)) set.seed(seed)
  ehsa <- emerging_hotspot_analysis(x = spt, .var = "tot_ab", threshold = 0.1, include_gi = TRUE, nsim = nsim)

  # Add the geometry: one polygon per pixel (taken from the first year)
  pixels_sf <- summed_sf[summed_sf$Year == min(summed_sf$Year), "pixel_id"]
  add_geometry(ehsa, pixels_sf)
}

# -----------------------------------------------------------------------------
# Run EHSA per life stage (~10min per stage with nsim = 999)
# -----------------------------------------------------------------------------

for (stage in stages_to_run) {
  t0 <- Sys.time()
  ehsa <- run_ehsa(stage, raster_stack, depth_mask, nsim = nsim, seed = seed)
  output_path <- paste0("./Outputs/EHSA/ehsa_", stage, "_depth_masked.rds")
  saveRDS(ehsa, file = output_path)
  step_msg("Saved EHSA output to: ", output_path, " (", round(as.numeric(Sys.time() - t0, units = "mins"), 1), " min)")
}

