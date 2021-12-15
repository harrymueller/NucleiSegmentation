accurate_plot <- function (data, # dataframe with y,x,value
                           legend_name = "Legend", 
                           filename = "temp.png", 
                           legend_space = 2, # inches to given legend
                           dpi=300, 
                           minres=1500, # minimum resolution in px
                           crop = FALSE, # whether to crop the image using ImageMagick
                           adjust = 1, # whether to limit to a given quantile (1 = no)
                           custom_colours = c(), # vector of colours
                           left_plot = FALSE # ggplot obj to plot the left of this
                           ) {
  names(data) = c("y", "x", "values")
  
  x = max(data$x)+1
  y = max(data$y)
  
  # diameter (in pixels) of 1 cell
  px = ceiling(minres/min(c(x,y)))
  height = y / px / dpi * 100
  width = x / px / dpi * 100
  
  legend_offset = (legend_space / 2) / width + 1
  
  if (adjust != 1) {
    q = quantile(data[[3]], probs = adjust)
    data[[3]][data[[3]] > q] = q
  }
  
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
            plot.background=element_blank(),
            legend.position = c(legend_offset, 0.5),
            plot.margin=grid::unit(c(0,0,0,0), "mm"))
  
  p$labels$fill <- legend_name
  
  if (length(custom_colours) > 0)
    p = p + scale_fill_gradientn(colours=custom_colours)
  
  if (!left_plot) {
    ggsave(filename, p, height=height, width=width + legend_space*2, dpi=dpi)
    
    if (crop) system(sprintf("convert %s -trim +repage %s", filename, filename))
  } else {
    ggsave(paste0(filename, ".TEMP"), p, height = height, width = width + legend_space * 2, dpi = dpi)
    ggsave(filename, left_plot, height = height, width = width + legend_space * 2, dpi = dpi)

    # crop accurate plot
    #system(sprintf("convert %s -trim +repage %s", filename, filename))
    system(sprintf("convert %s -trim +repage %s", paste0(filename, ".TEMP"), paste0(filename, ".TEMP")))

    # merge then rm temp file
    system(sprintf("convert %s %s +append %s", filename, paste0(filename, ".TEMP"), filename))
    system(paste0("rm ", filename, ".TEMP"))
  }
}