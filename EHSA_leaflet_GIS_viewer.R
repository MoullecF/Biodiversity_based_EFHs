###############################################################################
# Title: Interactive Leaflet GIS viewer for EHSA / hotspot-coldspot outputs
#
# Modular Leaflet map with:
#   - switchable basemaps (Esri World Topographic [default], Esri World
#     Imagery, Esri World Ocean Base)
#   - an overlay layer control supporting multiple independently togglable
#     polygon layers, each with its own colour palette, legend, popups and
#     opacity
#   - navigation tools (zoom, fullscreen, scale bar, mouse coordinates,
#     home/reset view, zoom-to-layer)
#
# The script is organised in the sections requested: configuration, data
# loading, CRS handling, map initialisation, basemaps, overlays, styling,
# navigation, legends/popups, export.
###############################################################################

# ==============================================================================
# 1. Package loading and configuration
# ==============================================================================

required_pkgs <- c("leaflet", "leaflet.extras", "leafem", "sf", "htmlwidgets", "htmltools", "RColorBrewer")
invisible(lapply(required_pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}))
library(leaflet)
library(leaflet.extras)
library(leafem)   # provides addMouseCoordinates()
library(sf)
library(htmlwidgets)
library(htmltools)

# Small null-coalescing helper used for optional spec fields (e.g. legend_columns)
`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- User configuration -----------------------------------------------------
# Edit this section only. Nothing below needs to change to add/remove layers.

cfg <- list(
  # Input spatial files (.rds/.RDS holding an sf object). Replace/extend as needed.
  # Each entry is a layer spec consumed by build_overlay() in section 6.
  # Leave `path = NULL` to register a placeholder layer with no data yet.
  input_files = list(
    combined_hs_cs = "./Outputs/EHSA/Combined_HS_CS.RDS",
    ehsa_adult     = "./Outputs/EHSA/EHSA_adult_90.rds",
    ehsa_juvenile  = "./Outputs/EHSA/EHSA_juvenile_90.rds",
    extrapolation_risk = "./Outputs/Extrapolation/extrapolation_risk.rds"
  ),

  # Output HTML path
  output_html = "./Outputs/EHSA/EHSA_GIS_viewer.html",
  html_title  = "EHSA Hotspot/Coldspot Interactive GIS Viewer - Western & Central Mediterranean",

  # Target CRS for all layers (leaflet requires EPSG:4326 / longlat)
  target_crs = 4326,

  # Initial map extent: western & central Mediterranean
  init_bbox = list(lng1 = -6, lat1 = 35, lng2 = 20, lat2 = 45),

  # Geometry simplification for large polygon layers (sf::st_simplify).
  # Set to NULL to disable; units are in degrees (target_crs is longlat).
  simplify_tolerance = NULL   # e.g. 0.001 to lightly simplify very dense grids
)

# ==============================================================================
# 2. Input spatial data
# ==============================================================================

# Generic, safe loader: returns NULL (with a message) rather than failing hard
# if a file is missing, so placeholder layers can be wired up before data exists.
load_sf_safely <- function(path, label) {
  if (is.null(path)) {
    message("[", label, "] No input file configured yet (placeholder layer).")
    return(NULL)
  }
  if (!file.exists(path)) {
    message("[", label, "] File not found: ", path, " -- registering as placeholder.")
    return(NULL)
  }
  obj <- tryCatch(readRDS(path), error = function(e) {
    message("[", label, "] Failed to read '", path, "': ", conditionMessage(e))
    NULL
  })
  if (is.null(obj)) return(NULL)
  if (!inherits(obj, "sf")) {
    message("[", label, "] Object at '", path, "' is not an sf object -- skipped.")
    return(NULL)
  }
  obj
}

# ==============================================================================
# 3. CRS validation and transformation
# ==============================================================================

# Validates geometries, drops empty/invalid ones (reporting counts, never silently
# dropping attributes), and reprojects to the target CRS used by Leaflet.
prepare_sf <- function(x, label, target_crs, simplify_tolerance = NULL) {
  if (is.null(x)) return(NULL)

  n0 <- nrow(x)

  if (is.na(sf::st_crs(x))) {
    stop("[", label, "] Input has no CRS defined -- cannot safely reproject. ",
         "Please assign the correct CRS with sf::st_set_crs() before use.")
  }

  x <- sf::st_make_valid(x)

  empty_idx <- sf::st_is_empty(x)
  if (any(empty_idx)) {
    message("[", label, "] Dropping ", sum(empty_idx), " empty geometr",
            ifelse(sum(empty_idx) == 1, "y", "ies"), " out of ", n0, ".")
    x <- x[!empty_idx, ]
  }

  na_geom_idx <- is.na(sf::st_dimension(x))
  if (any(na_geom_idx)) {
    message("[", label, "] Dropping ", sum(na_geom_idx), " row(s) with invalid/NA geometry.")
    x <- x[!na_geom_idx, ]
  }

  x <- sf::st_transform(x, crs = target_crs)

  if (!is.null(simplify_tolerance)) {
    x <- sf::st_simplify(x, dTolerance = simplify_tolerance, preserveTopology = TRUE)
  }

  message("[", label, "] Ready: ", nrow(x), "/", n0, " features retained.")
  x
}

# Converts a POINT sf grid into square POLYGON cells of the given cell size
# (degrees, since data is in longlat), centred on each point -- used for grids
# supplied as point centroids (e.g. extrapolation_risk) so they can be drawn
# with addPolygons() like the other overlay layers.
points_to_grid_cells <- function(x, cell_size, label) {
  if (is.null(x)) return(NULL)
  coords <- sf::st_coordinates(x)
  half <- cell_size / 2
  cells <- lapply(seq_len(nrow(x)), function(i) {
    cx <- coords[i, 1]; cy <- coords[i, 2]
    sf::st_polygon(list(matrix(
      c(cx - half, cy - half,
        cx + half, cy - half,
        cx + half, cy + half,
        cx - half, cy + half,
        cx - half, cy - half),
      ncol = 2, byrow = TRUE
    )))
  })
  out <- sf::st_sf(sf::st_drop_geometry(x), geometry = sf::st_sfc(cells, crs = sf::st_crs(x)))
  message("[", label, "] Converted ", nrow(out), " point(s) to ", cell_size, "-degree square cells.")
  out
}

# ==============================================================================
# 4. Map initialisation
# ==============================================================================

init_map <- function(bbox) {
  leaflet(options = leafletOptions(zoomControl = TRUE, minZoom = 3)) %>%
    fitBounds(lng1 = bbox$lng1, lat1 = bbox$lat1, lng2 = bbox$lng2, lat2 = bbox$lat2)
}

# ==============================================================================
# 5. Background map definitions
# ==============================================================================

# Named list so overlay/basemap controls and JS can refer to consistent labels.
# Only Esri basemaps are kept -- CartoDB/OSM tiles failed to switch reliably
# in testing, so they were dropped to keep the control simple and working.
basemap_names <- c(
  "Esri World Topographic"  = "Esri World Topographic",
  "Esri World Imagery"      = "Esri World Imagery",
  "Esri World Ocean Base"   = "Esri World Ocean Base"
)

add_basemaps <- function(map) {
  map %>%
    # Default basemap added first so it renders initially.
    # unname() is required: indexing a named vector keeps the name attribute,
    # which htmlwidgets serializes as a JSON object (e.g. {"CartoDB Positron":
    # "CartoDB Positron"}) instead of a plain string -- breaking group matching
    # in addLayersControl() and silently disabling basemap switching.
    addProviderTiles(providers$Esri.WorldTopoMap, group = unname(basemap_names["Esri World Topographic"])) %>%
    addProviderTiles(providers$Esri.WorldImagery, group = unname(basemap_names["Esri World Imagery"])) %>%
    # Esri Ocean Basemap: served as a standard Esri tile service (REST), valid as of
    # writing; if Esri retires/renames this service, swap the URL below.
    addTiles(
      urlTemplate = "https://services.arcgisonline.com/ArcGIS/rest/services/Ocean/World_Ocean_Base/MapServer/tile/{z}/{y}/{x}",
      attribution = "Esri, GEBCO, NOAA, National Geographic, DeLorme, HERE, Geonames.org, and other contributors",
      group = unname(basemap_names["Esri World Ocean Base"])
    )
}

# ==============================================================================
# 6. Overlay layer definitions (modular)
# ==============================================================================

# One overlay "spec" fully describes how to build and style a layer. New layers
# are added by appending a list() to `overlay_specs` below (see instructions at
# the bottom of this script) -- no other code needs to change.
#
# Fields:
#   id           - unique internal id (character)
#   label        - name shown in the layer control / legend
#   data         - sf object (POLYGON/MULTIPOLYGON), or NULL for a placeholder
#   fill_attr    - column name used to colour polygons (NULL = single colour)
#   palette      - named character vector (value -> hex colour), or a single
#                  hex colour string when fill_attr is NULL
#   popup_attrs  - character vector of column names to show in the popup
#                  (NULL = show all non-geometry attributes)
#   legend_title - title shown above the legend for this layer
#   default_on   - logical, whether the layer is visible when the map opens
#   border_color - polygon outline colour
#   opacity      - initial fill opacity (0-1); independently adjustable per layer
#                  via the opacity sliders added in section 7

build_popups <- function(df, popup_attrs = NULL) {
  attrs <- sf::st_drop_geometry(df)
  if (!is.null(popup_attrs)) {
    missing_cols <- setdiff(popup_attrs, names(attrs))
    if (length(missing_cols) > 0) {
      warning("Requested popup attribute(s) not found and will be skipped: ",
              paste(missing_cols, collapse = ", "))
    }
    attrs <- attrs[, intersect(popup_attrs, names(attrs)), drop = FALSE]
  }
  if (ncol(attrs) == 0) return(rep("", nrow(df)))
  # unname(): apply() returns a vector named by row names, which htmlwidgets
  # would serialise as a JSON object (jsonlite named-vector warning).
  unname(apply(attrs, 1, function(row) {
    rows_html <- paste0(
      "<tr><td style='font-weight:600;padding-right:8px;'>", htmltools::htmlEscape(names(row)), "</td>",
      "<td>", htmltools::htmlEscape(row), "</td></tr>",
      collapse = ""
    )
    paste0("<table style='font-size:12px;border-collapse:collapse;'>", rows_html, "</table>")
  }))
}

# Placeholder polygon layers: replace `data = NULL` with a real sf object, set
# `fill_attr`/`palette`/`popup_attrs` once the real attributes are known.
overlay_specs <- list()

# ==============================================================================
# Real data layers (Combined_HS_CS, EHSA_adult_90, EHSA_juvenile_90)
# Palettes reproduced from Figure_5_6.r (combined hotspot/coldspot gradients)
# and Figure_3_4.R (adult/juvenile EHSA classification palette).
# ==============================================================================

# --- Figure 5/6 palette (Combined_HS_CS) ------------------------------------
hs_palette_base <- colorRampPalette(RColorBrewer::brewer.pal(9, "YlGn"))(13)
cs_palette_base <- colorRampPalette(RColorBrewer::brewer.pal(9, "BuPu"))(22)

# Expands the abbreviated hotspot/coldspot category labels (e.g. "J.Int",
# "A.Spo") to readable text (e.g. "Juv.Intensifying", "Adu.Sporadic").
expand_hs_cs_labels <- function(x) {
  x <- gsub("J\\.", "Juv.", x)
  x <- gsub("A\\.", "Adu.", x)
  abbr_map <- c(
    "Int" = "Intensifying", "Per" = "Persistent", "Spo" = "Sporadic",
    "Con" = "Consecutive",  "dim" = "Diminishing"
  )
  for (abbr in names(abbr_map)) {
    x <- gsub(paste0("(Juv\\.|Adu\\.)", abbr, "\\b"),
              paste0("\\1", abbr_map[[abbr]]), x)
  }
  x
}

hs_legend_col <- c(
  "J.Int - A.Int" = hs_palette_base[13], "J.Int - A.Per" = hs_palette_base[12],
  "J.Int - A.Spo" = hs_palette_base[11], "J.Per - A.Int" = hs_palette_base[10],
  "J.Per - A.Per" = hs_palette_base[9],  "J.Per - A.Spo" = hs_palette_base[8],
  "J.dim - A.Spo" = hs_palette_base[7],  "J.Spo - A.Int" = hs_palette_base[6],
  "J.Spo - A.Per" = hs_palette_base[5],  "J.Spo - A.Spo" = hs_palette_base[4],
  "J.Con - A.Per" = hs_palette_base[3],  "J.New - A.Int" = hs_palette_base[2],
  "J.New - A.Spo" = hs_palette_base[1]
)

cs_legend_col <- c(
  "J.Int - A.Int" = cs_palette_base[22], "J.Int - A.Per" = cs_palette_base[21],
  "J.Int - A.Spo" = cs_palette_base[20], "J.Int - A.Con" = cs_palette_base[19],
  "J.Int - A.New" = cs_palette_base[18], "J.Per - A.Int" = cs_palette_base[17],
  "J.Per - A.Per" = cs_palette_base[16], "J.Per - A.Spo" = cs_palette_base[15],
  "J.Per - A.Con" = cs_palette_base[14], "J.Per - A.New" = cs_palette_base[13],
  "J.Spo - A.Int" = cs_palette_base[12], "J.Spo - A.Per" = cs_palette_base[11],
  "J.Spo - A.Spo" = cs_palette_base[10], "J.Spo - A.Con" = cs_palette_base[9],
  "J.Spo - A.New" = cs_palette_base[8],  "J.Con - A.Int" = cs_palette_base[7],
  "J.Con - A.Per" = cs_palette_base[6],  "J.Con - A.Spo" = cs_palette_base[5],
  "J.Con - A.Con" = cs_palette_base[4],  "J.New - A.Int" = cs_palette_base[3],
  "J.New - A.Per" = cs_palette_base[2],  "J.New - A.Spo" = cs_palette_base[1]
)

combined_legend_col <- c(hs_legend_col, cs_legend_col)
names(hs_legend_col) <- expand_hs_cs_labels(names(hs_legend_col))
names(cs_legend_col) <- expand_hs_cs_labels(names(cs_legend_col))
combined_legend_col <- c(hs_legend_col, cs_legend_col)

# --- Figure 3/4 palette (EHSA_adult_90 / EHSA_juvenile_90) ------------------
palette_ehsa <- c(
  "intensifying hotspot" = "#003c30", "persistent hotspot" = "#005F56",
  "diminishing hotspot" = "#23867E",  "sporadic hotspot" = "#5AB2A8",
  "consecutive hotspot" = "#97D6CD",  "new hotspot" = "#CEEBE7",
  "no pattern detected" = "#F5F5F5",  "new coldspot" = "#F5EACB",
  "consecutive coldspot" = "#E6CE94", "sporadic coldspot" = "#CFA154",
  "diminishing coldspot" = "#AE7121", "persistent coldspot" = "#824B09",
  "intensifying coldspot" = "#543005"
)

data_combined <- load_sf_safely(cfg$input_files$combined_hs_cs, "combined_hs_cs")
data_adult    <- load_sf_safely(cfg$input_files$ehsa_adult, "ehsa_adult")
data_juvenile <- load_sf_safely(cfg$input_files$ehsa_juvenile, "ehsa_juvenile")
data_extrapolation_risk <- load_sf_safely(cfg$input_files$extrapolation_risk, "extrapolation_risk")

data_combined <- prepare_sf(data_combined, "combined_hs_cs", cfg$target_crs, cfg$simplify_tolerance)
# Keep only the attributes needed for the combined layer, with clearer names.
if (!is.null(data_combined)) {
  data_combined <- data_combined[, c("classification.x", "classification.y", "Classification", "Spot")]
  names(data_combined)[names(data_combined) == "classification.x"] <- "Classification juveniles"
  names(data_combined)[names(data_combined) == "classification.y"] <- "Classification adults"
  names(data_combined)[names(data_combined) == "Classification"]   <- "Combined classification"
}

data_adult    <- prepare_sf(data_adult, "ehsa_adult", cfg$target_crs, cfg$simplify_tolerance)
data_juvenile <- prepare_sf(data_juvenile, "ehsa_juvenile", cfg$target_crs, cfg$simplify_tolerance)
# Keep only the classification attribute, with a clearer name.
if (!is.null(data_adult))    data_adult    <- data_adult[, "classification"]
if (!is.null(data_juvenile)) data_juvenile <- data_juvenile[, "classification"]
if (!is.null(data_adult))    names(data_adult)[names(data_adult) == "classification"]       <- "Classification"
if (!is.null(data_juvenile)) names(data_juvenile)[names(data_juvenile) == "classification"]  <- "Classification"
data_extrapolation_risk <- prepare_sf(data_extrapolation_risk, "extrapolation_risk", cfg$target_crs, cfg$simplify_tolerance)
# Source data is a point grid (0.05-degree spacing) -- convert to square
# polygon cells so it renders with the same addPolygons() pipeline as the
# other overlay layers.
data_extrapolation_risk <- points_to_grid_cells(data_extrapolation_risk, 0.05, "extrapolation_risk")
# Zero means "no extrapolation risk assessed" rather than a real risk value --
# treat as missing so it is excluded from the colour scale and rendered as NA.
if (!is.null(data_extrapolation_risk)) {
  data_extrapolation_risk$risk[data_extrapolation_risk$risk == 0] <- NA
}

# Classes present in the data but absent from the Figure 5/6 palette fall back to grey.
# Classification is a factor with unused levels (29 total, fewer observed) -- coerce to
# character so colorFactor() aligns strictly to the values actually present, not to all
# factor levels (which previously caused a palette/colour misalignment).
#
# Note: hs_legend_col and cs_legend_col share several identical classification
# strings (e.g. "J.Int - A.Int") with DIFFERENT colours -- they are only
# disambiguated by the `Spot` column (Hotspots vs Coldspots), not by the
# classification label itself. The HS-/CS- prefix must therefore be driven by
# `Spot`, not by checking which palette's names the label matches.
if (!is.null(data_combined)) {
  data_combined$`Combined classification` <- as.character(data_combined$`Combined classification`)
  data_combined$`Combined classification` <- expand_hs_cs_labels(data_combined$`Combined classification`)
  data_combined$`Classification juveniles` <- expand_hs_cs_labels(as.character(data_combined$`Classification juveniles`))
  data_combined$`Classification adults`   <- expand_hs_cs_labels(as.character(data_combined$`Classification adults`))

  spot_prefix <- ifelse(data_combined$Spot == "Hotspots", "HS - ",
                         ifelse(data_combined$Spot == "Coldspots", "CS - ", ""))
  data_combined$`Combined classification` <- paste0(spot_prefix, data_combined$`Combined classification`)

  hs_legend_col_tagged <- setNames(hs_legend_col, paste0("HS - ", names(hs_legend_col)))
  cs_legend_col_tagged <- setNames(cs_legend_col, paste0("CS - ", names(cs_legend_col)))
  combined_legend_col <- c(hs_legend_col_tagged, cs_legend_col_tagged)

  missing_classes <- setdiff(unique(data_combined$`Combined classification`), names(combined_legend_col))
  if (length(missing_classes) > 0) {
    combined_legend_col <- c(combined_legend_col, setNames(rep("grey70", length(missing_classes)), missing_classes))
  }
}
if (!is.null(data_adult))    data_adult$Classification    <- as.character(data_adult$Classification)
if (!is.null(data_juvenile)) data_juvenile$Classification <- as.character(data_juvenile$Classification)

# --- Marine protected areas (protection levels, minimally -> fully) ---------
mpa_polygons <- sf::st_read(
  "./Inputs_HMSC/protectionPLactivities_niv1&2&3&4/protectionPLactivities_niv1&2&3&4.shp",
  quiet = TRUE
)
# Fix any invalid geometries (self-intersections, etc.).
mpa_polygons <- sf::st_make_valid(mpa_polygons)

# Recode numeric protection codes to a labelled factor matching the palette names.
mpa_polygons$protctn <- factor(
  mpa_polygons$protctn,
  levels = c(1, 2, 3, 4),
  labels = c("Minimally", "Lightly", "Highly", "Fully")
)
# Only the protection level is plotted -- drop other attributes before building popups/fill.
mpa_polygons <- mpa_polygons["protctn"]
names(mpa_polygons)[names(mpa_polygons) == "protctn"] <- "MPA - Protection level"
mpa_polygons <- prepare_sf(mpa_polygons, "mpa_polygons", cfg$target_crs, cfg$simplify_tolerance)
if (!is.null(mpa_polygons)) mpa_polygons$`MPA - Protection level` <- as.character(mpa_polygons$`MPA - Protection level`)

protection_palette <- c(
  "Minimally"   = "#fde725",
  "Lightly"     = "#5ec962",
  "Highly"      = "#21918c",
  "Fully"       = "#440154"
)

real_specs <- list(
  list(
    id           = "combined_hs_cs",
    label        = "Combined hotspots/coldspots (adult x juvenile)",
    data         = data_combined,
    fill_attr    = "Combined classification",
    palette      = combined_legend_col,
    popup_attrs  = c("Classification juveniles", "Classification adults", "Combined classification", "Spot"),
    legend_title = "Combined classification",
    default_on   = TRUE,
    border_color = "#333333",
    opacity      = 0.75
  ),
  list(
    id           = "ehsa_adult",
    label        = "Adult EHSA (90% CI)",
    data         = data_adult,
    fill_attr    = "Classification",
    palette      = palette_ehsa,
    popup_attrs  = NULL,
    legend_title = "EHSA classification",
    default_on   = FALSE,
    border_color = "#333333",
    opacity      = 0.75
  ),
  list(
    id           = "ehsa_juvenile",
    label        = "Juvenile EHSA (90% CI)",
    data         = data_juvenile,
    fill_attr    = "Classification",
    palette      = palette_ehsa,
    popup_attrs  = NULL,
    legend_title = "EHSA classification",
    default_on   = FALSE,
    border_color = "#333333",
    opacity      = 0.75
  ),
  list(
    id           = "extrapolation_risk",
    label        = "Extrapolation risk",
    data         = data_extrapolation_risk,
    fill_attr    = "risk",
    palette      = c("white", "red"),
    popup_attrs  = NULL,
    legend_title = "Extrapolation risk",
    default_on   = FALSE,
    border_color = "transparent",
    opacity      = 0.75
  ),
  list(
    id           = "mpa_polygons",
    label        = "Protection level",
    data         = mpa_polygons,
    fill_attr    = "MPA - Protection level",
    palette      = protection_palette,
    popup_attrs  = "MPA - Protection level",
    legend_title = "Protection level",
    default_on   = FALSE,
    border_color = "#333333",
    opacity      = 0.75
  )
)

overlay_specs <- c(real_specs, overlay_specs)

# ==============================================================================
# 7. Layer styling and opacity controls
# ==============================================================================

# Adds a single overlay layer (polygons + legend) to the map, wiring up a
# colour palette, popups and group name used by the layer control.
# Populated as a side effect of add_overlay_layer() for every categorical
# (non-numeric) layer, so configure_legends_and_opacity() can build a
# checkbox-based legend (one checkbox per category + select/deselect all)
# instead of a plain colour-swatch addLegend().
category_legend_data <- list()

add_overlay_layer <- function(map, spec) {
  if (is.null(spec$data) || nrow(spec$data) == 0) {
    message("[", spec$id, "] No data available -- layer registered as empty placeholder ",
            "(will appear in the layer control but show nothing on the map).")
    return(map)
  }

  df <- spec$data
  df$.popup <- build_popups(df, spec$popup_attrs)

  is_numeric_fill <- !is.null(spec$fill_attr) && is.numeric(df[[spec$fill_attr]])

  fill_opacity <- spec$opacity

  if (is_numeric_fill) {
    values <- df[[spec$fill_attr]]
    pal_fun <- colorNumeric(palette = spec$palette, domain = values, na.color = "transparent")
    fill_colors <- pal_fun(values)
    # Zero (no extrapolation risk) is fully transparent; other cells use the
    # layer's configured opacity.
    fill_opacity <- ifelse(is.na(values) | values == 0, 0, spec$opacity)
  } else if (!is.null(spec$fill_attr)) {
    if (!spec$fill_attr %in% names(df)) {
      stop("[", spec$id, "] fill_attr '", spec$fill_attr, "' not found in data columns: ",
           paste(names(df), collapse = ", "))
    }
    # Coerce to character so colorFactor() aligns strictly to the values present
    # (a factor's unused levels would otherwise shift the colour assignment)
    values <- as.character(df[[spec$fill_attr]])
    if (is.character(spec$palette) && !is.null(names(spec$palette))) {
      pal_fun <- colorFactor(palette = unname(spec$palette[sort(unique(values))]),
                              domain = values, na.color = "grey80")
    } else {
      pal_fun <- colorFactor(palette = spec$palette, domain = values, na.color = "grey80")
    }
    fill_colors <- pal_fun(values)
  } else {
    fill_colors <- spec$palette   # single flat colour
  }

  # Categorical (non-numeric) layers get a checkbox-based legend (one checkbox
  # per category, plus select/deselect-all) instead of a plain colour-swatch
  # addLegend(). To let JS toggle visibility per category, each category is
  # added as its own addPolygons() call, tagged with a `category` option, while
  # all of them still share `group = spec$label` so the main layer control
  # continues to show/hide the whole layer.
  if (!is.null(spec$fill_attr) && !is_numeric_fill) {
    if (is.character(spec$palette) && !is.null(names(spec$palette))) {
      cat_order <- names(spec$palette)[names(spec$palette) %in% values]
    } else {
      cat_order <- sort(unique(values))
    }

    for (cat in cat_order) {
      idx <- values == cat
      if (!any(idx)) next
      map <- map %>%
        addPolygons(
          data = df[idx, ], group = spec$label,
          fillColor = fill_colors[idx], color = spec$border_color, weight = 0.4,
          fillOpacity = spec$opacity, popup = df$.popup[idx],
          label = values[idx],
          options = c(leaflet::pathOptions(), list(category = cat)),
          highlightOptions = highlightOptions(weight = 2, color = "black", bringToFront = TRUE)
        )
    }

    category_legend_data[[spec$id]] <<- list(
      id           = spec$id,
      label        = spec$label,
      legend_title = spec$legend_title,
      categories   = cat_order,
      colors       = unname(pal_fun(cat_order))
    )
  } else {
    map <- map %>%
      addPolygons(
        data = df, group = spec$label,
        fillColor = fill_colors, color = spec$border_color, weight = 0.4,
        fillOpacity = fill_opacity, popup = ~.popup,
        label = if (!is.null(spec$fill_attr)) values else NULL,
        highlightOptions = highlightOptions(weight = 2, color = "black", bringToFront = TRUE)
      )
  }

  if (is_numeric_fill) {
    map <- map %>%
      addLegend("bottomright", pal = pal_fun, values = values,
                title = spec$legend_title, group = spec$label,
                className = "info legend")
  }

  map
}

add_all_overlays <- function(map, specs) {
  for (spec in specs) map <- add_overlay_layer(map, spec)
  map
}

# ==============================================================================
# 8. Navigation controls
# ==============================================================================

add_navigation_controls <- function(map, overlay_labels) {
  map %>%
    addScaleBar(position = "bottomleft", options = scaleBarOptions(imperial = FALSE)) %>%
    addFullscreenControl(position = "topleft") %>%
    addMouseCoordinates() %>%
    addEasyButton(
      easyButton(
        icon = "fa-home", title = "Reset view (full study area)",
        onClick = JS("function(btn, map){ map.fitBounds(map._initialBounds); }")
      )
    ) %>%
    addLayersControl(
      baseGroups = unname(basemap_names),
      overlayGroups = overlay_labels,
      options = layersControlOptions(collapsed = TRUE),
      position = "topright"
    )
}

# ==============================================================================
# 9. Legend and popup configuration
# ==============================================================================

# Only show legends for overlays that are visible by default; keep others
# hidden until their layer is toggled on, and provide per-layer opacity sliders.
# Uses leaflet::hideGroup() (not manual map.removeLayer()) so the layer
# control's internal state stays in sync and toggling keeps working.
# Hand-rolled JSON encoder for category_legend_data (avoids a jsonlite
# dependency for this simple, fully-known structure: a list of
# id/label/legend_title/categories/colors entries).
encode_category_legend_json <- function(category_legend_data) {
  esc <- function(x) gsub('"', '\\\\"', x)
  entries <- vapply(category_legend_data, function(d) {
    cats_js <- paste0('"', esc(d$categories), '"', collapse = ", ")
    cols_js <- paste0('"', esc(d$colors), '"', collapse = ", ")
    sprintf(
      '{"id": "%s", "label": "%s", "legendTitle": "%s", "categories": [%s], "colors": [%s]}',
      esc(d$id), esc(d$label), esc(d$legend_title), cats_js, cols_js
    )
  }, character(1))
  paste0("[", paste(entries, collapse = ", "), "]")
}

configure_legends_and_opacity <- function(map, specs, bbox, category_legend_data) {
  # Force `map` immediately: in the %>% pipe, `map` is a lazy promise wrapping
  # the whole add_all_overlays() chain. add_overlay_layer() populates
  # category_legend_data as a *side effect* (via <<-) while building that
  # chain, so if `map` isn't forced before category_legend_data is read below,
  # the side effect hasn't happened yet and category_legend_data is captured
  # as still empty -- which silently dropped all checkbox legends.
  force(map)
  spec_labels <- vapply(specs, function(s) s$label, character(1))
  overlay_labels_js <- paste0('"', spec_labels, '"', collapse = ", ")
  default_on_js <- paste0(vapply(specs, function(s) tolower(as.character(isTRUE(s$default_on))), character(1)), collapse = ", ")
  category_legend_json <- encode_category_legend_json(category_legend_data)

  # Labels of the layers that got a plain addLegend() (numeric fill with data),
  # in the order the legends were added to the map. Used in JS to tag each
  # addLegend() legend with ITS OWN layer (previously they were matched to
  # overlayLabels by position, which attached e.g. the 'Extrapolation risk'
  # legend to the first overlay).
  has_numeric_legend <- vapply(specs, function(s) {
    !is.null(s$data) && nrow(s$data) > 0 && !is.null(s$fill_attr) &&
      s$fill_attr %in% names(s$data) && is.numeric(s$data[[s$fill_attr]])
  }, logical(1))
  numeric_legend_labels_js <- if (any(has_numeric_legend)) {
    paste0('"', spec_labels[has_numeric_legend], '"', collapse = ", ")
  } else {
    ""
  }

  # Hide non-default overlays properly via the leaflet R API, before any
  # custom JS runs -- this keeps the layer control checkboxes correct.
  hide_labels <- spec_labels[!vapply(specs, function(s) isTRUE(s$default_on), logical(1))]
  if (length(hide_labels) > 0) {
    map <- map %>% hideGroup(hide_labels)
  }

  init_bounds_js <- sprintf(
    "map._initialBounds = L.latLngBounds([%f, %f], [%f, %f]);",
    bbox$lat1, bbox$lng1, bbox$lat2, bbox$lng2
  )

  map %>%
    # NB: the JS is assembled with paste0(), not sprintf(): sprintf() refuses
    # format strings longer than 8192 characters, which this block exceeds.
    htmlwidgets::onRender(paste0("
      function(el, x) {
        var map = this;
        ", init_bounds_js, "

        var overlayLabels = [", overlay_labels_js, "];
        var defaultOn = [", default_on_js, "];
        var numericLegendLabels = [", numeric_legend_labels_js, "];

        // Returns the leaflet FeatureGroup holding ALL layers of an overlay
        // group, whether or not the group is currently shown on the map.
        // (map.eachLayer() only sees layers currently displayed, so it misses
        // every overlay hidden at load time via hideGroup().)
        function getGroupContainer(name) {
          return (map.layerManager && map.layerManager.getLayerGroup)
            ? map.layerManager.getLayerGroup(name, false) : null;
        }

        // Tag each plain addLegend() legend with its own group and set its
        // initial visibility to match the layer's default on/off state.
        // (Category legends do not exist yet at this point, so this selector
        // only matches the addLegend() ones, in the order they were added.)
        var legends = document.querySelectorAll('.leaflet-control.info.legend');
        legends.forEach(function(lg, idx) {
          var name = numericLegendLabels[idx];
          if (name) {
            lg.setAttribute('data-group', name);
            lg.style.display = defaultOn[overlayLabels.indexOf(name)] ? 'block' : 'none';
          }
        });

        // Show/hide the matching legend whenever its overlay is toggled
        map.on('overlayadd', function(e) {
          document.querySelectorAll('.leaflet-control.info.legend').forEach(function(lg) {
            if (lg.getAttribute('data-group') === e.name) lg.style.display = 'block';
          });
        });
        map.on('overlayremove', function(e) {
          document.querySelectorAll('.leaflet-control.info.legend').forEach(function(lg) {
            if (lg.getAttribute('data-group') === e.name) lg.style.display = 'none';
          });
        });

        // Per-layer opacity sliders, one per overlay with data
        var sliderPanel = L.control({position: 'bottomleft'});
        sliderPanel.onAdd = function(map) {
          var div = L.DomUtil.create('div', 'leaflet-bar');
          div.style.background = 'white';
          div.style.padding = '6px 8px';
          div.style.maxHeight = '160px';
          div.style.overflowY = 'auto';
          var html = '<div style=\"font-size:12px;font-weight:600;margin-bottom:4px;\">Layer opacity</div>';
          overlayLabels.forEach(function(name, i) {
            html += '<label style=\"font-size:11px;display:block;margin-top:4px;\">' + name + '</label>' +
              '<input type=\"range\" min=\"0\" max=\"1\" step=\"0.05\" value=\"0.75\" ' +
              'style=\"width:150px;\" data-group=\"' + name + '\" class=\"opacity-slider\">';
          });
          div.innerHTML = html;
          L.DomEvent.disableClickPropagation(div);
          return div;
        };
        sliderPanel.addTo(map);

        // Cache, once, every styleable layer of each overlay group -- whether
        // the group is visible or not -- so the opacity sliders also apply to
        // layers that are currently switched off and to unchecked categories.
        // Features drawn fully transparent on purpose (fillOpacity 0, i.e.
        // zero/NA extrapolation risk) are left out so they stay transparent.
        var allGroupLayers = {};
        overlayLabels.forEach(function(name) {
          var grp = getGroupContainer(name);
          allGroupLayers[name] = grp ? grp.getLayers().filter(function(layer) {
            return layer.setStyle && !(layer.options && layer.options.fillOpacity === 0);
          }) : [];
        });

        document.querySelectorAll('.opacity-slider').forEach(function(slider) {
          slider.addEventListener('input', function(e) {
            var group = e.target.getAttribute('data-group');
            var val = parseFloat(e.target.value);
            (allGroupLayers[group] || []).forEach(function(layer) {
              layer.setStyle({fillOpacity: val});
            });
          });
        });

        // ------------------------------------------------------------------
        // Checkbox-based category legends (one per categorical layer): a
        // 'Select all' / 'Deselect all' pair plus one checkbox per category,
        // in a scrollable box. All categories are checked by default.
        // Unchecking a category removes that category's polygons from the
        // layer's group container (and re-checking adds them back), so it
        // works whether or not the layer is currently displayed, and the
        // selection survives switching the whole layer off and on again in
        // the main layer control.
        var categoryLegendData = ", category_legend_json, ";

        categoryLegendData.forEach(function(entry) {
          // Collect the per-category leaflet layers added for this spec,
          // keyed by category value (tagged via options.category in R).
          // Read from the group container rather than map.eachLayer(), so
          // layers hidden at load time (default_on = FALSE) are found too.
          var groupContainer = getGroupContainer(entry.label);
          var categoryLayers = {};
          entry.categories.forEach(function(cat) { categoryLayers[cat] = []; });
          if (groupContainer) {
            groupContainer.eachLayer(function(layer) {
              if (layer.options && layer.options.category !== undefined) {
                if (!categoryLayers[layer.options.category]) categoryLayers[layer.options.category] = [];
                categoryLayers[layer.options.category].push(layer);
              }
            });
          }

          var ctrl = L.control({position: 'bottomright'});
          ctrl.onAdd = function(map) {
            var div = L.DomUtil.create('div', 'info legend category-legend');
            div.setAttribute('data-group', entry.label);
            div.style.background = 'white';
            div.style.padding = '6px 8px';
            div.style.fontSize = '11px';
            div.style.lineHeight = '1.3';

            var html = '<div style=\"font-weight:600;margin-bottom:4px;\">' + entry.legendTitle + '</div>';
            html += '<div style=\"margin-bottom:4px;\">' +
              '<a href=\"#\" class=\"select-all\" style=\"margin-right:8px;\">Select all</a>' +
              '<a href=\"#\" class=\"deselect-all\">Deselect all</a></div>';
            html += '<div class=\"category-scroll\" style=\"max-height:180px;overflow-y:auto;\">';
            entry.categories.forEach(function(cat, i) {
              var swatch = 'display:inline-block;width:10px;height:10px;margin-right:4px;' +
                'background:' + entry.colors[i] + ';border:1px solid #888;';
              html += '<label style=\"display:block;white-space:nowrap;\">' +
                '<input type=\"checkbox\" class=\"category-checkbox\" data-category=\"' + cat + '\" checked> ' +
                '<span style=\"' + swatch + '\"></span>' + cat + '</label>';
            });
            html += '</div>';
            div.innerHTML = html;
            L.DomEvent.disableClickPropagation(div);
            L.DomEvent.disableScrollPropagation(div);
            return div;
          };
          ctrl.addTo(map);

          // Use the control's own container (robust even if two layers were
          // ever given the same label).
          var legendEl = ctrl.getContainer();
          legendEl.style.display = defaultOn[overlayLabels.indexOf(entry.label)] ? 'block' : 'none';

          function setCategoryVisible(cat, visible) {
            if (!groupContainer) return;
            (categoryLayers[cat] || []).forEach(function(layer) {
              if (visible) {
                if (!groupContainer.hasLayer(layer)) groupContainer.addLayer(layer);
              } else {
                if (groupContainer.hasLayer(layer)) groupContainer.removeLayer(layer);
              }
            });
          }

          legendEl.querySelectorAll('.category-checkbox').forEach(function(cb) {
            cb.addEventListener('change', function(e) {
              setCategoryVisible(e.target.getAttribute('data-category'), e.target.checked);
            });
          });
          legendEl.querySelector('.select-all').addEventListener('click', function(e) {
            e.preventDefault();
            legendEl.querySelectorAll('.category-checkbox').forEach(function(cb) {
              cb.checked = true;
              setCategoryVisible(cb.getAttribute('data-category'), true);
            });
          });
          legendEl.querySelector('.deselect-all').addEventListener('click', function(e) {
            e.preventDefault();
            legendEl.querySelectorAll('.category-checkbox').forEach(function(cb) {
              cb.checked = false;
              setCategoryVisible(cb.getAttribute('data-category'), false);
            });
          });
        });

        // Show/hide each category legend alongside its layer's main toggle,
        // matching the behaviour already wired up for the plain addLegend()
        // (gradient) legends above.
        map.on('overlayadd', function(e) {
          document.querySelectorAll('.category-legend').forEach(function(lg) {
            if (lg.getAttribute('data-group') === e.name) lg.style.display = 'block';
          });
        });
        map.on('overlayremove', function(e) {
          document.querySelectorAll('.category-legend').forEach(function(lg) {
            if (lg.getAttribute('data-group') === e.name) lg.style.display = 'none';
          });
        });
      }
    "))
}

# ==============================================================================
# Assemble the map
# ==============================================================================

overlay_labels <- vapply(overlay_specs, function(s) s$label, character(1))

m <- init_map(cfg$init_bbox) %>%
  add_basemaps() %>%
  add_all_overlays(overlay_specs) %>%
  add_navigation_controls(overlay_labels) %>%
  configure_legends_and_opacity(overlay_specs, cfg$init_bbox, category_legend_data) %>%
  htmlwidgets::appendContent(htmltools::tags$title(cfg$html_title)) %>%
  # Shrink legend text so it doesn't overlap with the layer-control toggle.
  htmlwidgets::appendContent(htmltools::tags$style(HTML("
    .info.legend {
      font-size: 11px;
      line-height: 1;
      max-height: 220px;
      overflow-y: auto;
    }
    .info.legend .legend-title {
      font-size: 12px;
    }
  ")))

# ==============================================================================
# 10. HTML export
# ==============================================================================

dir.create(dirname(cfg$output_html), recursive = TRUE, showWarnings = FALSE)

# Note on self-containment: basemap tiles (OSM/CartoDB/Esri) are always fetched
# live from their providers -- they cannot be embedded, so an internet
# connection is required to see background imagery when the HTML is opened.
# All overlay polygon data, styling, legends, and JS/CSS for controls ARE
# embedded in the file via selfcontained = TRUE, so collaborators only need
# the single HTML file (plus internet access for basemap tiles).
htmlwidgets::saveWidget(m, file = normalizePath(cfg$output_html, mustWork = FALSE),
                         selfcontained = TRUE, title = cfg$html_title)

message("GIS viewer saved to: ", cfg$output_html)

###############################################################################
# How to add a new overlay layer later
#
# 1. Load your spatial object as an sf polygon layer, e.g.:
#      my_layer <- readRDS("./Inputs/my_new_layer.rds")
#
# 2. Prepare it (CRS check/transform, validity, optional simplification):
#      my_layer <- prepare_sf(my_layer, "my_layer", cfg$target_crs)
#
# 3. Append a new spec to overlay_specs (or real_specs, before the script
#    builds the map), for example:
#
#      list(
#        id           = "my_layer",
#        label        = "My new layer",
#        data         = my_layer,
#        fill_attr    = "my_attribute_column",   # or NULL for a flat colour
#        palette      = c("Class A" = "#1b9e77", "Class B" = "#d95f02"),
#        popup_attrs  = c("my_attribute_column", "some_other_column"),
#        legend_title = "My classification",
#        default_on   = FALSE,
#        border_color = "#333333",
#        opacity      = 0.75
#      )
#
# 4. Re-run the script. The new layer automatically gets a toggle in the
#    layer control, a legend (if fill_attr is set), popups, and an opacity
#    slider -- no other code changes are required.
###############################################################################
