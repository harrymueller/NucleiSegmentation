# LIBRARIES
library(Seurat)
library(ggplot2)
library(patchwork)
library(raster)
library(RColorBrewer)

# params
library(argparser)
args <- arg_parser("Plotting clusters on UMAP & spatial image")
args <- add_argument(args, "--id", help = "Tongue ID")
args <- add_argument(args, "--diameter", help="if subsetting, supply the diameter", default=NULL)
argv <- parse_args(args)

TONGUE_ID = argv$id
DIAMETER = argv$diameter

# other consts

INPUT_DIR = "/mnt/data"
OUTPUT_DIR = "/mnt/data/umap_clusters"
source("/mnt/data/scripts/rscripts/functions/accurate_plot.R")

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR)

# name, folder, binsizes, resolution
options = list(
    list("NormalizeData", "dimReducedRDS", c(100, 50, 20), c(0.5,0.5,0.5)),
    list("SCTransform", "scDimReducedRDS", c(100, 50), c(0.5, 0.5))
)
if (T) {
    options = list(
        list("NormalizeData", "dimReducedRDS", c(30), c(0.2))#,
        #list("SCTransform", "scDimReducedRDS", c(50), c(0.3))
    )
}
for (reduction in options) {
    bins = reduction[[3]]
    res = reduction[[4]]
    output = paste0(OUTPUT_DIR, "/", reduction[[1]])
    if (!dir.exists(output)) dir.create(output)

    for (i in seq(length(reduction[[3]]))) {
        # read in RDS
        if (is.null(DIAMETER)) {
            INPUT = sprintf("%s/%s/%s_bin%s_red.Rds", INPUT_DIR, reduction[[2]], TONGUE_ID, bins[i])
        } else {
            INPUT = sprintf("%s/%s/%s_bin%s_subset%s_red.Rds", INPUT_DIR, reduction[[2]], TONGUE_ID, bins[i], DIAMETER)
        }
        obj = readRDS(INPUT)
        
        obj <- FindClusters(obj, resolution = res[i])
        p1 <- DimPlot(obj, reduction = "umap", label = TRUE, pt.size = 1)
        p2 <- SpatialDimPlot(obj)

        if (is.null(DIAMETER)) {
            accurate_plot(
                p2$data,
                filename = sprintf("%s/%s_bin%s.png", output, TONGUE_ID, bins[i]),
                legend_name = sprintf("%s\n%s\nbin%s", reduction[[1]], TONGUE_ID, bins[i]),
                left_plot = p1,
                dpi = 500,
                minres = 1000
            )
        } else {
            accurate_plot(
                p2$data,
                filename = sprintf("%s/%s_bin%s_subset%s.png", output, TONGUE_ID, bins[i], DIAMETER),
                legend_name = sprintf("%s\n%s\nbin%s_subset%s", reduction[[1]], TONGUE_ID, bins[i], DIAMETER),
                left_plot = p1,
                dpi = 500,
                minres = 1000
            )
        }
    }
}

