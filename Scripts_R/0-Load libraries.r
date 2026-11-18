# Load all libraries
required_packages <- c(
  "tidyr", "plyr", "dplyr", "ape", "psych", "raster", "geosphere", "sf", "maps",
  "mapdata", "purrr", "phytools", "rprojroot", "Hmsc", "Matrix", "vioplot",
  "ggpubr", "wesanderson", "jsonify", "abind", "corrplot", "stringr",
  "reshape2", "ggplot2", "gridExtra", "grid", "cowplot", "viridis", "ggsci",
  "ggnewscale", "ggthemes", "patchwork", "gtable", "PupillometryR",
  "Routliers", "trend", "RColorBrewer", "ggtext", "sfdep", "spdep", "rnaturalearth",
  "rnaturalearthdata", "forcats", "ggpattern", "zenodor"
)

to_install <- setdiff(required_packages, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

invisible(lapply(required_packages, library, character.only = TRUE))