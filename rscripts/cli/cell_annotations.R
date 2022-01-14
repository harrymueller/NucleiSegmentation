##############################
# LIBS + PARAMS
##############################
library(Seurat)
library(SingleR)
library(ggplot2)

RefData = celldex::MouseRNAseqData()

# params
DIR = "/mnt/data"

# arg parsing
library(argparser)
args <- arg_parser("Cluster the seurat obj then plot a umap and spatial image")
args <- add_argument(args, "--id", help="TONGUE_ID")
args <- add_argument(args, "--binsize", help = "binsize")
args <- add_argument(args, "--method", help="method of normalisation, SCT || LN")
args <- add_argument(args, "--diameter", help="supply the diameter, 0 if not subsetting")
args <- add_argument(args, "--resolution", help="seurat FindCluster resolution", default=0.5)
argv <- parse_args(args)

TONGUE_ID   = argv$id
METHOD      = argv$method
BINSIZE     = as.integer(argv$binsize)
DIAMETER    = as.integer(argv$diameter)
RESOLUTION  = as.double(argv$resolution)

# other required variables
METHOD_NAME = ifelse(METHOD == "SCT", "SCTransform", "NormalizeData")
ASSAY_TO_USE = ifelse(METHOD == "SCT", "SCT", "Spatial")

FOLDER_NAME = sprintf("%s_bin%s_subset%s_res%s", TONGUE_ID, BINSIZE, DIAMETER, RESOLUTION)
FILENAME = paste0(FOLDER_NAME, ".rds")
INPUT = paste0(DIR, "/umap_clusters/", METHOD_NAME, "/RDS/", FILENAME)

# output dir
OUTPUT_DIR = sprintf("%s/cell_annotations/%s/%s", DIR, METHOD_NAME, FOLDER_NAME)
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR)

##############################
# FUNCTIONS
##############################
singleR_annotations <- function (bUseMainLabels, bByClusters) {
    if (bIsMain) {
        labels = RefData$label.main
    } else {
        labels = RefData$label.fine
    }

    Idents(obj) = "seurat_clusters"

    # ensure log-normalised counts
    if (bByClusters) {
    pred = SingleR(test = obj@assays$Spatial@data, 
                    ref = RefData, 
                    labels = labels, 
                    clusters = Idents(obj))
    } else {
    pred = SingleR(test = obj@assays$Spatial@data, 
                    ref = RefData, 
                    labels = labels)
    }

    ggsave(paste0(OUTPUT_DIR, "/singleR_heatmap_", ifelse(bIsMain, "main", "fine"), "_", ifelse(bByClusters, "clusters", "cells"), ".png"), p, width = 10, height = 7)

    annotation_label = paste("singleR", 
                            ifelse(bIsMain, "main", "fine"), 
                            ifelse(bByClusters, "clusters", "cells"), 
                            sep="_")
    if (bByClusters) {
    obj@meta.data[[annotation_label]] = factor(unfactor(obj@meta.data$seurat_clusters), 
                                                labels = pred$first.labels)
    } else {
    obj@meta.data[[annotation_label]] = pred$first.labels
    }
    Idents(obj) = annotation_label
    
    p = DimPlot(obj)
    ggsave(paste0(OUTPUT_DIR, "/singleR_umap_", ifelse(bIsMain, "main", "fine"), "_", ifelse(bByClusters, "clusters", "cells"), ".png"), p, width = 8, height = 8)

    return(obj)
}

##############################
# START
##############################
# read file
obj = readRDS(INPUT)

# using mainlabels, annotate clusters and cells
obj = singleR_annotations(T, F)
obj = singleR_annotations(T, T)
#obj = singleR_annotations(F, F)
#obj = singleR_annotations(F, T)

saveRDS(obj, paste0(OUTPUT_DIR, "/singleR.rds"))