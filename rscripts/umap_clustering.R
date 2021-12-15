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
argv <- parse_args(args)

TONGUE_ID = argv$id

# other consts
if (F) {
  INPUT_DIR = "/data/tongue"
  OUTPUT_DIR = "/data/tongue/umap_clusters"
  source("/data/tongue/scripts/rscripts/accurate_plot.R")
} else {
  INPUT_DIR = "/mnt/data/tongue_STOmics/discovery"
  OUTPUT_DIR = "/mnt/data/tongue_STOmics/discovery/umap_clusters"
  source("/mnt/local/scripts/rscripts/accurate_plot.R")
}
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR)

# name, folder, binsizes, resolution
options = list(
    list("NormalizeData", "dimReducedRDS", c(20, 50, 100), c(0.5,0.5,0.5)),
    list("SCTransform", "scDimReducedRDS", c(50, 100), c(0.5, 0.5))
)
if (F) {
    options = list(
        list("NormalizeData", "dimReducedRDS", c(100), c(0.5)),
        list("SCTransform", "scDimReducedRDS", c(100), c(0.5))
    )
}
for (reduction in options) {
    bins = reduction[[3]]
    res = reduction[[4]]
    output = paste0(OUTPUT_DIR, "/", reduction[[1]])
    if (!dir.exists(output)) dir.create(output)

    for (i in seq(length(reduction[[3]]))) {
        # read in RDS
        INPUT = sprintf("%s/%s/%s_bin%s_red.Rds", INPUT_DIR, reduction[[2]], TONGUE_ID, bins[i])
        obj = readRDS(INPUT)
        
        obj <- FindClusters(obj, resolution = res[i])
        p1 <- DimPlot(obj, reduction = "umap", label = TRUE, pt.size = 1)
        p2 <- SpatialDimPlot(obj)

        accurate_plot(
            p2$data,
            filename = sprintf("%s/%s_bin%s.png", output, TONGUE_ID, bins[i]),
            legend_name = sprintf("%s\n%s\nbin%s", reduction[[1]], TONGUE_ID, bins[i]),
            left_plot = p1,
            dpi = 500,
            minres = 1000
        )
    }
}

