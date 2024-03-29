# LIBRARIES
library(Seurat)
library(ggplot2)
library(patchwork)
library(raster)
library(RColorBrewer)

# params
library(argparser)
args <- arg_parser("Cluster the seurat obj then plot a umap and spatial image")
args <- add_argument(args, "--inputfile", help="input file")
args <- add_argument(args, "--outputdir", help="output dir")
args <- add_argument(args, "--name", help = "name")
args <- add_argument(args, "--method", help="method of normalisation, SCT || LN")
args <- add_argument(args, "--resolution", help="seurat FindCluster resolution", default=0.5)
argv <- parse_args(args)

INPUT       = argv$inputfile
OUTDIR      = argv$outputdir
NAME        = argv$name
RESOLUTION  = as.double(argv$resolution)

METHOD_NAME = ifelse(METHOD == "SCT", "SCTransform", "NormalizeData")

# other consts
source("/mnt/data/scripts/analysis/rscripts/functions/accurate_plot.R")

PLOTS_DIR = paste0(OUTPUT_DIR, "/plots")
RDS_DIR = paste0(OUTPUT_DIR, "/rds")
if (!dir.exists(PLOTS_DIR)) dir.create(PLOTS_DIR)
if (!dir.exists(RDS_DIR)) dir.create(RDS_DIR)

# read in RDS
obj = readRDS(INPUT)

# find clusters and plot
obj <- FindClusters(obj, resolution = RESOLUTION)
p1 <- DimPlot(obj, reduction = "umap", label = TRUE, pt.size = 1)
p2 <- SpatialDimPlot(obj)

accurate_plot(
    p2$data,
    filename = sprintf("%s/%s.png", PLOTS_DIR, NAME),
    legend_name = stringr::str_replace_all(NAME, "_", "\n"),
    left_plot = p1,
    dpi = 500,
    minres = 1000
)

saveRDS(obj, sprintf("%s/%s.rds", RDS_DIR, NAME))
