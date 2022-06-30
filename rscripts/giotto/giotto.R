
<<<<<<< HEAD
OUTPUT_DIR = "/mnt/data/R_analysis/page/tongue-4_corp1"
INPUT_DIR = OUTPUT_DIR
PART = "tongue-4-1.corp1"
SAMPLE_NAME = "tongue-5_bin30"
MARKERS = "/mnt/data/R_analysis/page/giotto_markers.Rds"
=======
OUTPUT_DIR = "/mnt/data/R_analysis_original/page"
INPUT_DIR = OUTPUT_DIR
PART = "Spl5Part1"
SAMPLE_NAME = "tongue-5_bin50"
>>>>>>> 3a3c725c5632f2e36462823a8163b72e27ccf499

### installing packages
if (!requireNamespace("Giotto", quietly=T)) remotes::install_github("RubD/Giotto") 
if (!requireNamespace("rio", quietly=T)) install.packages("rio")

library(Giotto)
library(Seurat)
library(Matrix)
library(ggplot2)
source("/mnt/data/scripts/rscripts/functions/accurate_plot.R")
library(patchwork)
library(raster)
library(RColorBrewer)
library(viridis)

installGiottoEnvironment() 

# DEG from SN Data
giotto_markers = readRDS(MARKERS)

## Convert to PAGE matrix
matrix = makeSignMatrixPAGE(sign_names = names(giotto_markers),
                            sign_list = giotto_markers)

# PAGE Enrichment
## Load spatial data
if (T) {
  counts = t(as(Matrix::readMM(file.path(INPUT_DIR, PART, paste0(PART, ".raw.count.mtx"))), "dgCMatrix"))
  gene_ids = read.csv(file.path(INPUT_DIR, PART, paste0(PART, ".raw.geneid.csv")), row.names = NULL)$row.names

  spatial_locations = read.csv(file.path(INPUT_DIR, PART, paste0(PART, ".raw.binid.csv")))
  # norm
  norm_counts = as(Matrix::readMM(file.path(INPUT_DIR, PART, paste0(PART, ".raw.corrected.count.mtx"))), "dgCMatrix")
  norm_gene_ids = read.csv(file.path(INPUT_DIR, PART, paste0(PART, ".raw.sct.geneid.csv")), row.names = NULL)$X

  norm_spatial_locations = read.csv(file.path(INPUT_DIR, PART, paste0(PART, ".raw.sct.binid.csv")))
  #subsetting
  mask = gene_ids %in% norm_gene_ids
  new_gene_ids = gene_ids[mask]
  new_counts = counts[mask,]
  # new
  giotto = createGiottoObject(new_counts, spatial_locs = spatial_locations[c("x","y")], gene_metadata = gene_ids, norm_expr = norm_counts)
} else {
  seurat = readRDS(file.path("/mnt/data/R_analysis_original/gemRDS", "tongue-5_bin50_spatialObj.rds"))

  counts = seurat@assays$Spatial@counts
  spatial_locations = seurat@images$slice1@coordinates[c("row", "col")]
  colnames(spatial_locations) = c("x", "y")

  giotto = createGiottoObject(counts, spatial_locs = spatial_locations)
}

# normalise then run giotto
giotto <- normalizeGiotto(gobject = giotto)
giotto = runPAGEEnrich(gobject = giotto,
                       sign_matrix = matrix,
                       min_overlap_genes = 2)

saveRDS(giotto, file.path(OUTPUT_DIR, "giotto_complete.Rds"))

# plotting
cell_types = names(giotto@spatial_enrichment$PAGE)[2:16]


page_enrichment = giotto@spatial_enrichment$PAGE
coords = giotto@spatial_locs
page_enrichment = merge(page_enrichment, coords, by = "cell_ID")

saveRDS(page_enrichment, file.path(OUTPUT_DIR, "page_enrichment.Rds"))

# plot enrichment scores
print("plot enrichment scores")
for (ct in cell_types) {
  print(ct)
  dat = cbind(page_enrichment[,c("sdimy", "sdimx")], page_enrichment[[ct]])
  accurate_plot(dat, filename = file.path(OUTPUT_DIR, paste0(SAMPLE_NAME, "_", ct, ".png")), custom_colours = viridis(11))
}

page = page_enrichment
# ranking
ranked_page = page
for (i in seq(2, 16)) {
  ranked_page[[i]] = order(order(page[[i]], decreasing = TRUE))
}

# minimise - ie find highest rank
page$ranks = apply(ranked_page[,2:16],1,which.min)
page$ranks = factor(page$ranks, labels = names(page[,2:16]))

# max
page$max = factor(apply(page[,2:16],1,which.max), labels = names(page[,2:16]))

# normalised
norm_page = page
for (i in seq(2, 16)) {
  norm_page[[i]] = (norm_page[[i]] - mean(norm_page[[i]])) / sd(norm_page[[i]])
}
page$norm = factor(apply(norm_page[,2:16],1,which.max), labels = names(page[,2:16]))

# plotting
custom_colours = c("#F56867",    "#FEB915",    "#C798EE",    "#59BE86",    "#7495D3",
"#D1D1D1",    "#6D1A9C",    "#15821E",    "#3A84E6",    "#997273",
"#787878",    "#DB4C6C",    "#9E7A7A",    "#554236",    "#AF5F3C",
"#93796C",    "#F9BD3F",    "#DAB370",    "#877F6C",    "#268785")
for (n in c("ranks", "max", "norm")) {
  accurate_plot(
    cbind(page[,c("sdimy", "sdimx")], page[[n]]),
    filename = file.path(OUTPUT_DIR, "plots", paste0(PART,"_",n,"_labels.png")),
    legend_name = paste0("PAGE Cell Labels - ", n), 
    dpi = 400,
    minres = 1000,
    legend_space = 4,
    custom_colours = custom_colours
  )
}

vdat = data.frame()
for (i in seq(2, 16)) {
  tempdf = data.frame("cell_type" = names(page)[i], "scores" = page[[i]])
  if (i == 2) {
    vdat = tempdf
  } else {
    vdat = rbind(vdat, tempdf)
  }
}
rm(tempdf)
p = ggplot(vdat, aes(x = cell_type, y = scores, fill = cell_type)) + 
  geom_violin() +
  scale_fill_manual(values = custom_colours)+ 
  theme_minimal() +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks=element_blank())
ggsave(file.path(OUTPUT_DIR, "plots", paste0(PART,"_violins.png")), p, width = 14, dpi = 500)