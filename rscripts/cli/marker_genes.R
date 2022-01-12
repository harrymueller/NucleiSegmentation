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

source("/mnt/data/scripts/rscripts/functions/accurate_plot.R")

# params
DIR = "/mnt/data"

# to do multiple runs in 1 go
DATA = data.frame(
  "id"          = c("tongue-5", "tongue-5", "tongue-5", "tongue-5", "tongue-4", "tongue-4", "tongue-4", "tongue-5", "tongue-5", "tongue-4", "tongue-4"),
  "binsize"     = c(50, 100, 20, 30, 50, 100, 30, 50, 100, 50, 100),
  "norm_method" = c("SCTransform", "SCTransform", "SCTransform", "SCTransform", "SCTransform", "SCTransform", "SCTransform", "NormalizeData", "NormalizeData", "NormalizeData", "NormalizeData"),
  "diameter"    = c(0, 0, 40, 25, 0, 0, 25, 0, 0, 0, 0)
)
#  "res"         = c(0.4),

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
                features = markers$top3_markers$genes, assay = "SCT", label=F) +
          ggtitle(sprintf("Heatmap of 3 markers genes from each cluster\n%s_bin%s%s downsampled to 1000", 
                            current$id, current$binsize, 
                            ifelse(current$diameter == 0, "", paste0("_subset",current$diameter)))) + 
          theme(text = element_text(size = 14))
  # The following features were omitted as they were not found in the scale.data slot for the SCT assay: Acta1, Tnni2, Myh4
  ggsave(paste0(output_dir, "/heatmap.png"), p, width = 8, height = 6)
}

plot_vln <- function (obj, markers) {
  p = VlnPlot(obj, features=markers$top3_markers$genes, assay = "SCT", stack = T, pt.size = 0, flip=T) +
        NoLegend() + 
        ggtitle(sprintf("Vln plot of 3 markers genes from each cluster\n%s_bin%s%s", 
                        current$id, current$binsize, 
                        ifelse(current$diameter == 0, "", paste0("_subset",current$diameter)))) + 
        theme(text = element_text(size = 28),
              plot.background = element_rect(fill = "white"))
  ggsave(paste0(output_dir, "/vln_plot.png"), p, width = 15, height = 15)
}

dot_plot <- function (obj, markers) {
  p = DotPlot(obj, features = unique(markers$top3_markers$genes), assay = "SCT", dot.scale = 3, cols="RdYlBu") + 
        ggtitle(sprintf("Dot plot of 3 markers genes from each cluster\n%s_bin%s%s", 
                        current$id, current$binsize, 
                        ifelse(current$diameter == 0, "", paste0("_subset",current$diameter)))) + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 14),
              text = element_text(size = 17),
              plot.background = element_rect(fill = "white"))
        
  ggsave(paste0(output_dir, "/dot_plot_of_markers.png"), p, width = 9, height = 7)
}

plot_marker_gene_spatial <- function (obj, markers) {
  # store genes used
  genes = c()
  
  # create and save plots
  for (clus in levels(Idents(obj))) {
    g = markers$top3_markers$genes[markers$top3_markers$cluster == clus]
    if (length(g) != 0) { # some markers
      g = g[1]
      genes = c(genes, g)
      
      accurate_plot(SpatialFeaturePlot(obj, features = g)$data, 
                    filename = paste0(output_dir, "/", g, "_spatial.png"), 
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
      base_file = paste0(output_dir, "/", genes[base_image_index], "_spatial.png")
      new_file = paste0(output_dir, "/", genes[i], "_spatial.png")
      
      system(sprintf("convert %s %s +append %s", base_file, new_file, base_file))
      system(paste0("rm ", new_file))
    }
  }
  
  # combine plots into 1
  if (n > 4) {
    for (i in seq(5, (ceiling(n / 4) - 1) * 4 + 1, by = 4)) {
      base_file = paste0(output_dir, "/", genes[1], "_spatial.png")
      new_file = paste0(output_dir, "/", genes[i], "_spatial.png")
      
      system(sprintf("convert %s %s -append %s", base_file, new_file, base_file))
      system(paste0("rm ", new_file))
    }
    
    system(sprintf("mv %s %s", 
                   paste0(output_dir, "/", genes[1], "_spatial.png"), 
                   paste0(output_dir, "/", "marker_genes_spatial_plot.png")))
  }
}

##############################
# START
##############################
for (it in seq(1, length(rownames(DATA)))) {
    current = DATA[it,]
    print("##############################")
    print(paste0(
        "Starting ", current$id, "_bin", current$binsize, 
        ifelse(current$diameter == 0, "", paste0("_subset",current$diameter)), " ..."
    ))
    print("##############################")

    # create filenames
    if (current$diameter == 0) { # not a subset
        filename = sprintf("%s_bin%s_red.Rds", current$id, current$binsize)
        folder_name = sprintf("%s_bin%s", current$id, current$binsize)
    } else {
        filename = sprintf("%s_bin%s_subset%s_red.Rds", current$id, current$binsize, current$diameter)
        folder_name = sprintf("%s_bin%s_subset%s", current$id, current$binsize, current$diameter)
    }

    # output dir
    output_dir = sprintf("%s/markers/%s", DIR, folder_name)
    if (!dir.exists(output_dir)) dir.create(output_dir)

    # read data
    input = sprintf(
        "%s/%s/%s", DIR, ifelse(current$norm_method == "SCTransform", "scDimReducedRDS", "dimReducedRDS"), filename
    )
    obj = readRDS(input)

    # markers
    markers = find_markers(obj)
    markers = get_all_markers(markers)

    save_markers(markers, paste0(output_dir, "/markers.xlsx"))

    # plots
    plot_heatmap(obj, markers)
    plot_vln(obj, markers)
    dot_plot(obj, markers)

    plot_marker_gene_spatial(obj, markers)

    # done
    print("##############################")
    print(paste0(
        "Finished ", current$id, "_bin", current$binsize, 
        ifelse(current$diameter == 0, "", paste0("_subset",current$diameter), "")
    ))
    print("##############################")
    gc()
}