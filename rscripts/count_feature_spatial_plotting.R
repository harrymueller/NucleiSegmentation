# LIBRARIES
library(Seurat)
library(ggplot2)
library(patchwork)
library(raster)
library(RColorBrewer)

# params
library(argparser)
args <- arg_parser("Plotting nCount, nFeature, Malat1 and Neat1 violin and spatial plots.")
args <- add_argument(args, "--binsize", help = "Bin Size")
args <- add_argument(args, "--id", help = "Tongue ID")
argv <- parse_args(args)

BIN_SIZE = argv$binsize
TONGUE_ID = argv$id

# other consts
INPUT_DIR = "/mnt/data/tongue_STOmics/discovery/gemRDS"
OUTPUT_DIR = sprintf("/mnt/data/tongue_STOmics/discovery/count_feature_plots/%s_bin%s", TONGUE_ID, BIN_SIZE)

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR)

# takes SpatialPlot data from seurat and saves it to a file
accurate_plot <- function (plot_data, tongue_id, filename, calibration = FALSE, binsize = 10) {
  data = plot_data$data
  
  names(data) = c("y", "x", "layer")
  
  if (calibration) {
    xmx = dim(data)[2]
    ymx = dim(data)[1]
    
    # calibration squares -> 10px in and 10px sides
    val = max(data[[3]])
    for (i in 12:22) 
      for (j in 12:22)
        data = rbind(data, c(i, j, val))
    for (i in (ymx-21):(ymx-11)) 
      for (j in (xmx-21):(xmx-11))
        data = rbind(data, c(i, j, val))
  }
  
  # any counts greater than 99.999% quantile are reduced to that quantile
  q = quantile(data[[3]], probs = 0.99999)
  data[[3]][data[[3]] > q] = q
  
  p = ggplot(data=data, mapping=aes(x=x, y = y)) + 
    geom_raster(aes(fill=layer), hjust=0, vjust=0) +
    coord_equal() + 
    scale_fill_gradientn(colours=c("white",rev(brewer.pal(n = 11, name = "Spectral")))) +
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
          plot.background=element_rect(fill="black"))
  gc()
  if (binsize == 1) dpimult = 10
  else dpimult = 1
  # tongue4 15:9:125
  if (tongue_id == "tongue-4") ggsave(filename, height=15, width=9, dpi=125*dpimult)
  # tongue5 15:11:104
  if (tongue_id == "tongue-5") ggsave(filename, height=15, width=11, dpi=104*dpimult)
}

# read in RDS
INPUT = sprintf("%s/%s_bin%s_spatialObj.rds", INPUT_DIR, TONGUE_ID, BIN_SIZE)
obj = readRDS(INPUT)
gc()

# make violin plots
for (t in c("nCount_Spatial", "nFeature_Spatial", "Malat1", "Neat1")) {
  p = VlnPlot(obj, features = t, pt.size = 0.1) + NoLegend()
  ggsave(sprintf("%s/vln_%s.png", OUTPUT_DIR, t), p)
}
gc()

# plot spatial plots
accurate_plot(SpatialPlot(obj, features = "nCount_Spatial"), TONGUE_ID, paste0(OUTPUT_DIR, "/nCount_Spatial.png"), binsize=BIN_SIZE)
gc()
accurate_plot(SpatialPlot(obj, features = "nFeature_Spatial"), TONGUE_ID, paste0(OUTPUT_DIR, "/nFeature_Spatial.png"), binsize=BIN_SIZE)
gc()
accurate_plot(SpatialPlot(obj, features = "Malat1"), TONGUE_ID, paste0(OUTPUT_DIR, "/Malat1_dist.png"), binsize=BIN_SIZE, scale = F)
gc()
accurate_plot(SpatialPlot(obj, features = "Neat1"), TONGUE_ID, paste0(OUTPUT_DIR, "/Neat1_dist.png"), binsize=BIN_SIZE, scale = F)
gc()