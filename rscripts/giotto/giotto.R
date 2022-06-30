
OUTPUT_DIR = "/mnt/data/R_analysis_original/page"
INPUT_DIR = OUTPUT_DIR
PART = "Spl5Part1"
SAMPLE_NAME = "tongue-5_bin50"

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
giotto_markers = readRDS(file.path(OUTPUT_DIR, "giotto_markers.Rds"))

## Convert to PAGE matrix
matrix = makeSignMatrixPAGE(sign_names = names(giotto_markers),
                            sign_list = giotto_markers)

# PAGE Enrichment
## Load spatial data
if (F) {
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

for (ct in cell_types) {
  print(ct)
  dat = cbind(page_enrichment[,c("sdimy", "sdimx")], page_enrichment[[ct]])
  accurate_plot(dat, filename = file.path(OUTPUT_DIR, paste0(SAMPLE_NAME, "_", ct, ".png")), custom_colours = viridis(11))
}




