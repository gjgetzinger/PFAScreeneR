library(gtable)
library(cowplot)

shift_legend <- function(p){

  # check if p is a valid object
  if(!"gtable" %in% class(p)){
    if("ggplot" %in% class(p)){
      gp <- ggplotGrob(p) # convert to grob
    } else {
      message("This is neither a ggplot object nor a grob generated from ggplotGrob. Returning original plot.")
      return(p)
    }
  } else {
    gp <- p
  }

  # check for unfilled facet panels
  facet.panels <- grep("^panel", gp[["layout"]][["name"]])
  empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
  empty.facet.panels <- facet.panels[empty.facet.panels]
  if(length(empty.facet.panels) == 0){
    message("There are no unfilled facet panels to shift legend into. Returning original plot.")
    return(p)
  }

  # establish extent of unfilled facet panels (including any axis cells in between)
  empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
  empty.facet.panels <- list(min(empty.facet.panels[["t"]]), min(empty.facet.panels[["l"]]),
                             max(empty.facet.panels[["b"]]), max(empty.facet.panels[["r"]]))
  names(empty.facet.panels) <- c("t", "l", "b", "r")

  # extract legend & copy over to location of unfilled facet panels
  guide.grob <- which(gp[["layout"]][["name"]] == "guide-box")
  if(length(guide.grob) == 0){
    message("There is no legend present. Returning original plot.")
    return(p)
  }
  gp <- gtable_add_grob(x = gp,
                        grobs = gp[["grobs"]][[guide.grob]],
                        t = empty.facet.panels[["t"]],
                        l = empty.facet.panels[["l"]],
                        b = empty.facet.panels[["b"]],
                        r = empty.facet.panels[["r"]],
                        name = "new-guide-box")

  # squash the original guide box's row / column (whichever applicable)
  # & empty its cell
  guide.grob <- gp[["layout"]][guide.grob, ]
  if(guide.grob[["l"]] == guide.grob[["r"]]){
    gp <- gtable_squash_cols(gp, cols = guide.grob[["l"]])
  }
  if(guide.grob[["t"]] == guide.grob[["b"]]){
    gp <- gtable_squash_rows(gp, rows = guide.grob[["t"]])
  }
  gp <- gtable_remove_grobs(gp, "guide-box")

  return(gp)
}



# make discrete gradient scales
# https://github.com/tidyverse/ggplot2/issues/2673#issuecomment-402878574
discrete_gradient_pal <- function(colours, bins = 5) {
  ramp <- scales::colour_ramp(colours)
  function(x) {
    if (length(x) == 0) return(character())

    i <- floor(x * bins)
    i <- ifelse(i > bins-1, bins-1, i)
    ramp(i/(bins-1))
  }
}

scale_colour_discrete_gradient <- function(..., colours, bins = 5, na.value = "grey50", guide = "colourbar", aesthetics = "colour", colors)  {
  colours <- if (missing(colours))
    colors
  else colours
  continuous_scale(
    aesthetics,
    "discrete_gradient",
    discrete_gradient_pal(colours, bins),
    na.value = na.value,
    guide = guide,
    ...
  )
}


plot_graph <- function(graph, x, y, x_lab, y_lab, trans){
  x <- pull(as_tibble(graph), !!x)
  y <- pull(as_tibble(graph), !!y)
  ggraph(
    graph,
    layout = 'manual',
    x = x,
    y = y
  ) +
    geom_node_point(aes(size = `Max area`, color = `Fluorine count`), alpha = 0.75) +
    scale_colour_discrete_gradient(
      colours = viridisLite::viridis(5, begin = 0.1, end = 0.99),
      limits = c(0,25),
      breaks =seq(0, 25, by = 5),
      guide = guide_colourbar(nbin = 200, raster = FALSE, frame.colour = "black", ticks.colour = NA, title.position = "top")
    ) +
    scale_size(trans = trans, guide = guide_legend(title.position = 'top')) +
    geom_edge_link(aes(col = factor(rxn_type)), arrow = arrow(type = "closed", length = unit(0.05, "inches"))) +
    scale_edge_color_brewer(palette = "Set1", guide_edge_colorbar(title = 'Trans. type', title.position = 'top')) +
    geom_node_text(aes(label = label), check_overlap = T, repel = T, color = 'gray66', size = 2) +
    theme_minimal() +
    labs(x = x_lab, y = y_lab) +
    theme(legend.position = 'bottom')
}
