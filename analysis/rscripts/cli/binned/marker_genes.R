##############################
# LIBS + PARAMS
##############################
# libraries
library(Seurat)
library(ggplot2)
library(xlsx)
library(dplyr)

library(patchwork)
library(raster)
library(RColorBrewer)

source("/mnt/data/scripts/analysis/rscripts/functions/accurate_plot.R")

# params
DIR = "/mnt/data/R_analysis"

# arg parsing
library(argparser)
args <- arg_parser("Cluster the seurat obj then plot a umap and spatial image")
args <- add_argument(args, "--input_file", help="input file")
args <- add_argument(args, "--output_dir", help="output dir")
args <- add_argument(args, "--method", help="method of normalisation, SCT || LN")
argv <- parse_args(args)

INPUT       = argv$input_file
OUTPUT_DIR  = argv$output_dir
METHOD      = argv$method

# other required variables
ASSAY_TO_USE = ifelse(METHOD == "SCT", "SCT", "Spatial") # KEEP

# output dir
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR)

##############################
# FUNCTIONS
##############################
find_markers <- function(obj) {
  clus.idents = levels(obj)
  markers = list("all_markers" = c("placeholder"), 
                 "top3_markers" = c("placeholder"), 
                 "empty" = c("Clusters with no Markers"))

  for (id in clus.idents) {
    m = FindMarkers(obj, ident.1 = id,
                    logfc.threshold = 1, min.pct = 0.5, test.use = "wilcox")
    markers[[id]] = m
    markers[[id]]$genes = rownames(m)
    
    # adjust and sort
    if (dim(markers[[id]])[1] == 0) {
  		markers[[id]] <- NULL
  		markers$empty = c(markers$empty, id)
  		print(paste0("Found no markers for cluster ", id))
    } else {
  	  markers[[id]] = markers[[id]][markers[[id]]$p_val_adj < 0.05,]
  	  markers[[id]] = markers[[id]][order(markers[[id]]$avg_log2FC, decreasing = T),]
  	  print(paste0("Found markers for cluster ", id))
  	}
  }
  
  return(markers)
}

# concat markers to the format shown in the report
get_all_markers <- function (markers) {
  markers$all_markers = data.frame()
  markers$top3_markers = data.frame()
  
  for (i in seq(4, length(names(markers)))) {
    temp = markers[[i]]
    temp$cluster = names(markers)[[i]]
    
    markers$all_markers = rbind(markers$all_markers, temp)
    
    # sort by abs log2FC val
    temp = temp[order(abs(temp$avg_log2FC), decreasing = T),][1:3,]
    markers$top3_markers = rbind(markers$top3_markers, temp)
  }
  
  markers$all_markers$cluster = factor(markers$all_markers$cluster)
  markers$top3_markers$cluster = factor(markers$top3_markers$cluster)
  
  markers$top3_markers = markers$top3_markers[!is.na(markers$top3_markers$cluster),]
  
  return(markers)
}

# saves marker gene info to an excel file
save_markers <- function(markers, f_name) {
  wb <- createWorkbook()
  
  for (id in names(markers)) {
    s = createSheet(wb, id)
    
    if (id == "empty") {
      addDataFrame(markers[[id]], s, col.names = F, row.names = F)
    } else {
      addDataFrame(markers[[id]], s)
    } 
    
  }
  
  saveWorkbook(wb, f_name)
}

plot_heatmap <- function (obj, markers) {
  p = DoHeatmap(subset(obj, downsample = 1000), 
                features = markers$top3_markers$genes, assay = ASSAY_TO_USE, label=F) +
          ggtitle(sprintf("Heatmap of 3 markers genes from each cluster\n%s_bin%s%s downsampled to 1000", 
                            TONGUE_ID, BINSIZE, 
                            ifelse(DIAMETER == 0, "", paste0("_subset",DIAMETER)))) + 
          theme(text = element_text(size = 14))
  # The following features were omitted as they were not found in the scale.data slot for the SCT assay: Acta1, Tnni2, Myh4
  ggsave(paste0(OUTPUT_DIR, "/heatmap.png"), p, width = 8, height = 6)
}

plot_vln <- function (obj, markers) {
  p = VlnPlot(obj, features=markers$top3_markers$genes, assay = ASSAY_TO_USE, stack = T, pt.size = 0, flip=T) +
        NoLegend() + 
        ggtitle(sprintf("Vln plot of 3 markers genes from each cluster\n%s_bin%s%s", 
                        TONGUE_ID, BINSIZE, 
                        ifelse(DIAMETER == 0, "", paste0("_subset",DIAMETER)))) + 
        theme(text = element_text(size = 28),
              plot.background = element_rect(fill = "white"))
  ggsave(paste0(OUTPUT_DIR, "/vln_plot.png"), p, width = 15, height = 15)
}

dot_plot <- function (obj, markers) {
  p = DotPlot(obj, features = unique(markers$top3_markers$genes), assay = ASSAY_TO_USE, dot.scale = 3, cols="RdYlBu") + 
        ggtitle(sprintf("Dot plot of 3 markers genes from each cluster\n%s_bin%s%s", 
                        TONGUE_ID, BINSIZE, 
                        ifelse(DIAMETER == 0, "", paste0("_subset",DIAMETER)))) + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 14),
              text = element_text(size = 17),
              plot.background = element_rect(fill = "white"))
        
  ggsave(paste0(OUTPUT_DIR, "/dot_plot_of_markers.png"), p, width = 9, height = 7)
}

plot_marker_gene_spatial <- function (obj, markers, colours) {
  # store genes used
  genes = c()
  
  # create and save plots
  for (clus in levels(Idents(obj))) {
    g = markers$top3_markers$genes[markers$top3_markers$cluster == clus]
    if (length(g) != 0) { # some markers
      g = g[1]
      genes = c(genes, g)
      
      accurate_plot(SpatialFeaturePlot(obj, features = g)$data, 
                    filename = paste0(OUTPUT_DIR, "/", g, "_spatial.png"), 
                    legend_name = g,
                    custom_colours = colours,
                    dpi = 500,
                    minres = 500,
                    crop = F,
                    black_background = F,
                    legend_space = 1)
    }
  }
  
  genes = unique(genes)
  
  # combine plots into rows of 4
  n = length(genes)
  for (i in seq(1, n)) { # 4 wide
    base_image_index = (ceiling(i / 4) - 1) * 4 + 1
    if (base_image_index != i) {
      base_file = paste0(OUTPUT_DIR, "/", genes[base_image_index], "_spatial.png")
      new_file = paste0(OUTPUT_DIR, "/", genes[i], "_spatial.png")
      
      system(sprintf("convert %s %s +append %s", base_file, new_file, base_file))
      system(paste0("rm ", new_file))
    }
  }
  
  # combine plots into 1
  if (n > 4) {
    for (i in seq(5, (ceiling(n / 4) - 1) * 4 + 1, by = 4)) {
      base_file = paste0(OUTPUT_DIR, "/", genes[1], "_spatial.png")
      new_file = paste0(OUTPUT_DIR, "/", genes[i], "_spatial.png")
      
      system(sprintf("convert %s %s -append %s", base_file, new_file, base_file))
      system(paste0("rm ", new_file))
    }
  }
  
  system(sprintf("mv %s %s", 
                  paste0(OUTPUT_DIR, "/", genes[1], "_spatial.png"), 
                  paste0(OUTPUT_DIR, "/", "marker_genes_spatial_plot.png")))

}

##############################
# START
##############################
# read file
obj = readRDS(INPUT)

# markers
markers = find_markers(obj)
markers = get_all_markers(markers)

save_markers(markers, paste0(OUTPUT_DIR, "/markers.xlsx"))

# plots
plot_heatmap(obj, markers)
plot_vln(obj, markers)
dot_plot(obj, markers)

colours = rev(brewer.pal(n = 11, name = "Spectral"))
plot_marker_gene_spatial(obj, markers, colours)