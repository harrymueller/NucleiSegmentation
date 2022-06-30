library(ggplot2)

accurate_plot <- function (data, # dataframe with y,x,value
                           legend_name = "Legend", 
                           filename = "temp.png", 
                           legend_space = 2, # inches to given legend
                           dpi=300, 
                           minres=1500, # minimum resolution in px
                           crop = FALSE, # whether to crop the image using ImageMagick
                           adjust = 1, # whether to limit to a given quantile (1 = no)
                           custom_colours = c(), # vector of colours
                           left_plot = c(), # ggplot obj to plot the left of this
                           black_background = FALSE, # whether to plot over a black background
                           spot_mappings = NULL # pass a spot mapping if required
                          ) {
  if (!is.null(spot_mappings)) {
    names(data) = c("imagerow", "imagecol", "values") # ensures correct names
    data = merge(data, spot_mappings, by.x = c("imagerow", "imagecol"), by.y = c("cy", "cx"))
    data = data[c("y", "x", "values")]
  } else {
    names(data) = c("y", "x", "values")
  }
  
  x = max(data$x)+1
  y = max(data$y)
  
  # diameter (in pixels) of 1 cell
  px = ceiling(minres/min(c(x,y)))
  height = px * y / dpi
  width = px * x / dpi
  
  legend_offset = (legend_space / 2) / width + 1
  
  if (adjust != 1) {
    q = quantile(data[[3]], probs = adjust)
    data[[3]][data[[3]] > q] = q
  }

  background = element_rect(fill=ifelse(black_background, "black", "white"))
  #background = element_rect(fill="grey")
  p = ggplot(data=data, mapping=aes(x=x, y = y)) + 
      geom_raster(aes(fill=values), hjust=0, vjust=0) + 
      coord_equal() +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) + 
      theme(axis.line=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=background,
            plot.margin=grid::unit(c(0,0,0,0), "mm"),
            legend.position = c(legend_offset, 0.5),
            legend.background = element_rect(fill="white"))
  
  p$labels$fill <- legend_name
  
  if (length(custom_colours) > 0)
    p = p + scale_fill_manual(values = custom_colours) #scale_fill_gradient2(low = "yellow", high="red")
  if (length(left_plot) == 0) {
    ggsave(filename, p, height=height, width=width + legend_space*2, dpi=dpi, limitsize = FALSE)
    if (crop) system(sprintf("convert %s -trim +repage %s", filename, filename))
  } else { 
    temp_filename = stringr::str_replace(filename, ".png", ".TEMP.png")

    # save as images -> ensures correct aspect ratio for spatial plots
    ggsave(temp_filename, p, height = height, width = width + legend_space * 2, dpi = dpi, limitsize = FALSE)
    ggsave(filename, left_plot, height = height / 0.8, width = (width + legend_space * 1.5) / 0.8, dpi = dpi*0.8, limitsize = FALSE)

    # crop
    system(sprintf("convert -crop %sx%s+%s+%s +repage %s %s", (width + legend_space * 1.2) * dpi, (height) * dpi, legend_space*0.8*dpi, 0, temp_filename, temp_filename))

    # merge then rm temp file
    system(sprintf("convert %s %s +append %s", filename, temp_filename, filename))
    system(paste0("rm ", temp_filename))

    #system(sprintf("convert %s -trim +repage %s", filename, filename))
  }
}
