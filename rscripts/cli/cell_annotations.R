##############################
# LIBS
##############################
suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(SingleR))
suppressMessages(library(dplyr))
suppressMessages(library(tibble))
suppressMessages(library(xlsx))

# temp until rebuilt docker
if (!requireNamespace("scibetR", quietly = TRUE))
  devtools::install_github("zwj-tina/scibetR")
suppressMessages(library(scibetR))

# celldex for reference dataset
if (!requireNamespace("celldex", quietly = TRUE))
  BiocManager::install("celldex")

# FN for accurate spatial plotting
source("/mnt/data/scripts/rscripts/functions/accurate_plot.R")

##############################
# PARAMS
##############################
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

print("############################################################")
print(paste0("Starting ", FOLDER_NAME, "..."))
print("############################################################")

# ref data
RefData = celldex::MouseRNAseqData()

##############################
# FUNCTIONS
##############################
# plot umaps and spatial dimensions of both cell annotations and clusters in a 2x2 grid
dim_plots <- function(obj, ident, filename) {
    temp_filename = stringr::str_replace_all(filename, ".png", ".TEMP.png")

    Idents(obj) = ident
    p1 <- DimPlot(obj, reduction = "umap", label = F, pt.size = 0.5) + theme(legend.position = "none")
    p2 <- SpatialDimPlot(obj)

    accurate_plot(
        p2$data,
        filename = filename,
        legend_name = paste0("Cell Annotations\n", FOLDER_NAME), # stringr::str_replace_all(FOLDER_NAME, "_", "\n")), 
        left_plot = p1,
        dpi = 400,
        minres = 1000,
        legend_space = 4
    )

    Idents(obj) = "seurat_clusters"
    p1 <- DimPlot(obj, reduction = "umap", label = F, pt.size = 0.5) + theme(legend.position = "none")
    p2 <- SpatialDimPlot(obj)

    accurate_plot(
        p2$data,
        filename = temp_filename,
        legend_name = "Seurat Clusters", 
        left_plot = p1,
        dpi = 400,
        minres = 1000,
        legend_space = 4
    )

    system(sprintf("convert %s %s -append %s", filename, temp_filename, filename))
    system(paste("rm", temp_filename))
}

# singleR annotations
singleR_annotations <- function (bUseMainLabels, bByClusters) {
    print(sprintf("Annotating using singleR (%s %s)...", bUseMainLabels, bByClusters))
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
    ggsave(paste0(OUTPUT_DIR, "/singleR_heatmap_", ifelse(bUseMainLabels, "main", "fine"), "_", ifelse(bByClusters, "clusters", "cells"), ".png"), p, width = 12, height = 7, dpi = 200)

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
    dim_plots(obj, annotation_label, paste0(OUTPUT_DIR, "/singleR_umap_", ifelse(bUseMainLabels, "main", "fine"), "_", ifelse(bByClusters, "clusters", "cells"), ".png"))

    return(obj)
}

# sciebet annotations
scibet_annotations <- function (bUseMainLabels) {
    print(sprintf("Annotating using SciBet (%s)...", bUseMainLabels))

    # create the training set for the celldex ref
    train = as.data.frame(RefData@assays@data)
    train[[1]] = NULL
    train[[1]] = NULL

    # which ref labels to use, and what to call the metadata
    if (bUseMainLabels) {
        refLabels = RefData$label.main
        labelName = "scibet_main"
    } else {
        refLabels = RefData$label.fine
        labelName = "scibet_fine"
    }

    # add labels to training set
    train = rbind(train, as.integer(as.factor(refLabels)))
    rownames(train)[length(rownames(train))] = "label"

    # train and test sets
    train = as.data.frame(t(train))
    test = as.data.frame(t(as.matrix(obj@assays[[ASSAY_TO_USE]]@data)))

    # annotate cells
    pred = SciBet_R(train, test)

    # add labels to seurat obj
    obj@meta.data[[labelName]] = pred
    labels = levels(factor(refLabels))
    for (i in seq(length(labels))) {
        obj@meta.data[[labelName]][obj@meta.data[[labelName]] == i] = labels[i]
    }

    # umap
    dim_plots(obj, labelName, paste0(OUTPUT_DIR, "/scibet_umap_", ifelse(bUseMainLabels, "main", "fine"), ".png"))

    # add labels to test dataset
    test = cbind(test, obj@meta.data[[labelName]])
    colnames(test)[length(colnames(test))] = "label"

    # find top 50 genes (as determined by SciBet fn), then do a dotplot
    genes = SelectGene_R(train, k = 50)
    
    p = scibet::Marker_heatmap(test, genes) + 
                         ggtitle(sprintf("Dot plot of 50 informative training genes\nas determined by the SciBet Entropy-Test\n%s_bin%s%s", 
                                            TONGUE_ID, BINSIZE, 
                                            ifelse(DIAMETER == 0, "", paste0("_subset",DIAMETER)))) + 
                        theme(text = element_text(size = 26), 
                              axis.text.y = element_text(size = 24),
                              axis.text.x = element_text(vjust = 0.5, size = 20),
                              legend.title = element_text(size = 24),
                              legend.text = element_text(size = 20))
    ggsave(paste0(OUTPUT_DIR, "/scibet_dotplot_", ifelse(bUseMainLabels, "main", "fine"), ".png"), p, width = 20, height = 12, dpi = 200)
    return(obj)
}

# produce a df of the counts of cts
ct_xtabs_df <- function(obj, col) {
  df = data.frame(xtabs(~obj@meta.data[[col]]))
  names(df) = c("Cell Type", col)
  #rownames(df) = df[["Cell Type"]]
  #df[["Cell Type"]] = NULL
  df[["Cell Type"]] = unfactor(df[["Cell Type"]])
  return(df)
}

# save the xtabs of all annotations to a xlsx file
save_annotations <- function(obj, f_name) {
  print("Saving annotations...")
  main = ct_xtabs_df(obj, "singleR_main_cells")
  main = merge.data.frame(main, ct_xtabs_df(obj, "singleR_main_clusters"), all = T)
  main = merge.data.frame(main, ct_xtabs_df(obj, "scibet_main"), all = T)
  rownames(main) = main[["Cell Type"]]
  main[["Cell Type"]] = NULL
  
  fine = ct_xtabs_df(obj, "singleR_fine_cells")
  fine = merge.data.frame(fine, ct_xtabs_df(obj, "singleR_fine_clusters"), all = T)
  fine = merge.data.frame(fine, ct_xtabs_df(obj, "scibet_fine"), all = T)
  rownames(fine) = fine[["Cell Type"]]
  fine[["Cell Type"]] = NULL
  
  # save to file
  wb <- createWorkbook()
  s = createSheet(wb, "Cell Annotations")
  
  # main labels
  addDataFrame(data.frame(paste0("Cell annotations of ", FOLDER_NAME, ", using celldex::MouseRNAseqData() 'main' labels")), 
               s, col.names = F, row.names = F)
  addDataFrame(t(main), s, col.names = T, row.names = T, startRow = 2)
  
  # fine labels
  addDataFrame(data.frame(paste0("Cell annotations of ", FOLDER_NAME, ", using celldex::MouseRNAseqData() 'fine' labels")), 
               s, col.names = F, row.names = F, startRow = 7)
  addDataFrame(t(fine), s, col.names = T, row.names = T, startRow = 8)
  
  saveWorkbook(wb, f_name)
}

##############################
# START
##############################
# read file
obj = readRDS(INPUT)

# annotate
obj = singleR_annotations(T, F)
obj = singleR_annotations(T, T)
obj = singleR_annotations(F, F)
obj = singleR_annotations(F, T)

obj = scibet_annotations(T)
obj = scibet_annotations(F)

save_annotations(obj, paste0(OUTPUT_DIR, "/cell_annotations.xlsx"))
saveRDS(obj, paste0(OUTPUT_DIR, "/annotated_seurat.rds"))

# keeps getting created
system("rm /mnt/data/scripts/Rplots.pdf")

print("############################################################")
print(paste0("Finished ", FOLDER_NAME, "."))
print("############################################################")