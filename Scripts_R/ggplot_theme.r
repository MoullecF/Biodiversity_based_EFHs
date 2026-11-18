# ------------------------------------------------------------------------------
# ggplot theme and legend placement helpers
# ------------------------------------------------------------------------------

## World map data (sf)
# ------------------------------------------------------------------------------

world = rnaturalearth::ne_countries(scale = "large", returnclass = "sf")

## THEME ggplot ----
HMSC.theme <- ggplot2::theme(
  legend.position = "right",
  axis.title.x = element_text(face = "bold", size = 14),
  axis.text.x  = element_text(vjust = 0.5, size = 12, colour = "black"),
  axis.title.y = element_text(face = "bold", size = 14),
  axis.text.y  = element_text(vjust = 0.5, size = 12, colour = "black"),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.minor.y = element_blank(),
  panel.grid.major.y = element_blank(),
  panel.background = element_blank(),
  axis.line = element_line(colour = "black"),
  plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5),
  # Use explicit units to avoid parse errors when sourcing.
  plot.margin = unit(c(0, 0, 0, 0), "cm")
)

## Utility: shift legend into empty facet panels
# ------------------------------------------------------------------------------

shift_legend <- function(p) {

  # Check if input is a ggplot or gtable
  if (!"gtable" %in% class(p)) {
    if ("ggplot" %in% class(p)) {
      gp <- ggplotGrob(p)
    } else {
      message("This is neither a ggplot object nor a grob generated from ggplotGrob. Returning original plot.")
      return(p)
    }
  } else {
    gp <- p
  }

  # Find empty facet panels
  facet.panels <- grep("^panel", gp[["layout"]][["name"]])
  empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
  empty.facet.panels <- facet.panels[empty.facet.panels]
  if (length(empty.facet.panels) == 0) {
    message("There are no unfilled facet panels to shift legend into. Returning original plot.")
    return(p)
  }

  # Determine target panel extent
  empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
  empty.facet.panels <- list(
    min(empty.facet.panels[["t"]]), min(empty.facet.panels[["l"]]),
    max(empty.facet.panels[["b"]]), max(empty.facet.panels[["r"]])
  )
  names(empty.facet.panels) <- c("t", "l", "b", "r")

  # Move legend into the empty panel area
  guide.grob <- which(gp[["layout"]][["name"]] == "guide-box-right")
  if (length(guide.grob) == 0) {
    message("There is no legend present. Returning original plot.")
    return(p)
  }
  gp <- gtable_add_grob(
    x = gp,
    grobs = gp[["grobs"]][[guide.grob]],
    t = empty.facet.panels[["t"]],
    l = empty.facet.panels[["l"]],
    b = empty.facet.panels[["b"]],
    r = empty.facet.panels[["r"]],
    name = "new-guide-box"
  )

  # Remove original guide box and squash its space
  guide.grob <- gp[["layout"]][guide.grob, ]
  if (guide.grob[["l"]] == guide.grob[["r"]]) {
    gp <- gtable_squash_cols(gp, cols = guide.grob[["l"]])
  }
  if (guide.grob[["t"]] == guide.grob[["b"]]) {
    gp <- gtable_squash_rows(gp, rows = guide.grob[["t"]])
  }
  gp <- gtable_remove_grobs(gp, "guide-box-right")

  return(gp)
}
