# LIBRARIES
library(Seurat)
library(ggplot2)
library(patchwork)
library(raster)
library(RColorBrewer)
library(viridis)

# params
library(argparser)
args <- arg_parser("Plotting nCount, nFeature, Malat1 and Neat1 violin and spatial plots.")
args <- add_argument(args, "--binsize", help = "Bin Size")
args <- add_argument(args, "--id", help = "Tongue ID")
args <- add_argument(args, "--diameter", help="if subsetting, supply the diameter", default=NULL)
argv <- parse_args(args)

BIN_SIZE  = strtoi(argv$binsize)
TONGUE_ID = argv$id
DIAMETER  = NULL # argv$diameter

GENES = c("Krt36", "Krt84", "Krt35", "Krt6b", "Krt6a", "Dmd", "Neb", "Ttn", "Rbfox1", "C1qtnf7", "Dclk1", "Mmp16", "Mecom", "Plcb1", "Flt1", "Fetub", "Ecm1", "Fmo2", "Pappa", "Hs3st5", "St3gal4", "Mki67", "Cenpf", "Cenpe", "Frem2", "Sema3c", "Pcdh7", "Unc5c", "Cpm", "Hhip", "Hephl1") 

# other consts
OUTPUT_DIR = sprintf("/mnt/data/R_analysis/count_feature_plots/%s_bin%s", TONGUE_ID, BIN_SIZE)
if (is.null(DIAMETER)) {
  INPUT_DIR = "/mnt/data/R_analysis/gemRDS"
} else {
  INPUT_DIR = "/mnt/data/R_analysis/subsets"
  OUTPUT_DIR = sprintf("%s_subset%s", OUTPUT_DIR, DIAMETER)
}
source("/mnt/data/scripts/analysis/rscripts/functions/accurate_plot.R")
  
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR)

# read in RDS
if (is.null(DIAMETER)) {
  INPUT = sprintf("%s/%s_bin%s_spatialObj.rds", INPUT_DIR, TONGUE_ID, BIN_SIZE)
} else {
  INPUT = sprintf("%s/%s_bin%s_subset%s.rds", INPUT_DIR, TONGUE_ID, BIN_SIZE, DIAMETER)
}
obj = readRDS(INPUT)

# plots
#colours = c("white",rev(brewer.pal(n = 11, name = "Spectral")))
colours = rev(brewer.pal(n = 11, name = "Spectral"))
i = 0
#for (t in c("nCount_Spatial", "nFeature_Spatial", "Malat1", "Neat1")) {
for (t in GENES) {
  if (BIN_SIZE > 11) {
    #left_plot = VlnPlot(obj, features = t, pt.size = 0.05) + NoLegend()
    isGene = t %in% GENES
    accurate_plot(SpatialPlot(obj, features = t)$data, 
                  filename = paste0(OUTPUT_DIR, "/", t, ".png"), 
                  legend_name = paste0(t, ifelse(!isGene && BIN_SIZE <= 10, "\n(q0.99999)", "")),
                  adjust = ifelse(!isGene && BIN_SIZE <= 10, 0.9999, 1),
                  custom_colours = colours, #left_plot = left_plot,
                  dpi = 750, # increasing DPI increases size of scales etc.
                  minres = 1500)
  } else {
    # plot vln and spatial plot separately for bin <= 2
    #left_plot = VlnPlot(obj, features = t, pt.size = 0.05) + NoLegend()
    #ggsave(paste0(OUTPUT_DIR, "/", t, "_vln.png"), 
    #       left_plot,
    #       width = 4)
    isGene = t %in% GENES
    accurate_plot(SpatialPlot(obj, features = t)$data, 
                  filename = paste0(OUTPUT_DIR, "/", t, "_spatial.png"), 
                  legend_name = paste0(t, ifelse(!isGene && BIN_SIZE <= 10, "\n(q0.99999)", "")),
                  adjust = ifelse(!isGene && BIN_SIZE <= 10, 0.9999, 1),
                  custom_colours = colours,
                  dpi = 750, #2000,
                  minres = 1500, #1,
                  crop = TRUE,
                  black_background = TRUE)
  }
  gc()
} 


