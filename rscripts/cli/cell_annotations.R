##############################
# LIBS + PARAMS
##############################
suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(SingleR))
if (!requireNamespace("scibetR", quietly = TRUE))
  devtools::install_github("zwj-tina/scibetR")
suppressMessages(library(scibetR))
suppressMessages(library(dplyr))

if (!requireNamespace("celldex", quietly = TRUE))
  BiocManager::install("celldex")
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
    # refdata labels to use
    if (bUseMainLabels) {
        labels = RefData$label.main
    } else {
        labels = RefData$label.fine
    }

    # ensure ident is seurat clusters
    Idents(obj) = "seurat_clusters"

    # ensure log-normalised counts
    if (bByClusters) {
        pred = SingleR(test = obj@assays[[ASSAY_TO_USE]]@data, 
                        ref = RefData, 
                        labels = labels, 
                        clusters = Idents(obj))
    } else {
        pred = SingleR(test = obj@assays[[ASSAY_TO_USE]]@data, 
                        ref = RefData, 
                        labels = labels)
    }

    # heatmap
    p = plotScoreHeatmap(pred)
    ggsave(paste0(OUTPUT_DIR, "/singleR_heatmap_", ifelse(bUseMainLabels, "main", "fine"), "_", ifelse(bByClusters, "clusters", "cells"), ".png"), p, width = 10, height = 7)

    # labels
    annotation_label = paste("singleR", 
                            ifelse(bUseMainLabels, "main", "fine"), 
                            ifelse(bByClusters, "clusters", "cells"), 
                            sep="_")
    if (bByClusters) {
        obj@meta.data[[annotation_label]] = factor(unfactor(obj@meta.data$seurat_clusters), 
                                                   labels = pred$first.labels)
    } else {
        obj@meta.data[[annotation_label]] = pred$first.labels
    }

    # umap
    Idents(obj) = annotation_label
    p = DimPlot(obj)
    ggsave(paste0(OUTPUT_DIR, "/singleR_umap_", ifelse(bUseMainLabels, "main", "fine"), "_", ifelse(bByClusters, "clusters", "cells"), ".png"), p, width = 8, height = 8)

    return(obj)
}

scibet_annotations <- function (bUseMainLabels) {
    train = as.data.frame(RefData@assays@data)
    train[[1]] = NULL
    train[[1]] = NULL

    if (bUseMainLabels) {
        refLabels = RefData$label.main
        labelName = "scibet_main"
    } else {
        refLabels = RefData$label.fine
        labelName = "scibet_fine"
    }

    # add labels
    train = rbind(train, as.integer(as.factor(refLabels)))
    rownames(train)[length(rownames(train))] = "label"

    # train and test sets
    train = as.data.frame(t(train))
    test = as.data.frame(t(as.matrix(obj@assays[[ASSAY_TO_USE]]@data)))

    # annotate cells
    pred = SciBet_R(train, test)

    # add labels to seurat obj
    obj@meta.data$scibet = pred
    labels = levels(factor(RefData$label.main))
    for (i in seq(length(labels))) {
        obj@meta.data$scibet[obj@meta.data$scibet == i] = labels[i]
    }

    # umap
    Idents(obj) = "scibet"
    p = DimPlot(obj)
    ggsave(paste0(OUTPUT_DIR, "/scibet_umap_", ifelse(bIsMain, "main", "fine"), ".png"), p, width = 8, height = 8)

    # dotpot
    g = SelectGene_R(test, k = 30)
    p = DotPlot(obj, features = g, assay = ASSAY_TO_USE, dot.scale = 3, cols="RdYlBu") + 
            ggtitle(sprintf("Dot plot of 3 markers genes from each cluster\n%s_bin%s%s", 
                            TONGUE_ID, BINSIZE, 
                            ifelse(DIAMETER == 0, "", paste0("_subset",DIAMETER)))) + 
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 14),
                text = element_text(size = 17),
                plot.background = element_rect(fill = "white"))
    ggsave(paste0(OUTPUT_DIR, "/scibet_dotplot_", ifelse(bIsMain, "main", "fine"), ".png"), p, width = 9, height = 7)
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

obj = scibet_annotations(T)

saveRDS(obj, paste0(OUTPUT_DIR, "/annotated_seurat.rds"))