
OUTPUT_DIR = "/mnt/data/R_analysis/page"

PART = "Spl5Part1"

### installing packages
if (!requireNamespace("Giotto", quietly=T)) remotes::install_github("RubD/Giotto") 
if (!requireNamespace("rio", quietly=T)) install.packages("rio")

library(Giotto)
library(Seurat)

# DEG from SN Data
giotto_markers = readRDS(file.path(OUTPUT_DIR, "giotto_markers.Rds"))

## Convert to PAGE matrix
matrix = makeSignMatrixPAGE(sign_names = names(giotto_markers),
                            sign_list = giotto_markers)

# PAGE Enrichment
## Load spatial data
counts = Matrix::readMM(file.path(INPUT_DIR, PART, paste0(PART, ".raw.count.mtx")))
gene_ids = read.csv(file.path(INPUT_DIR, PART, paste0(PART, ".raw.geneid.csv")), row.names = NULL)

spatial_locations = read.csv(file.path(INPUT_DIR, PART, paste0(PART, ".raw.binid.csv")))

giotto = createGiottoObject(t(counts), spatial_locs = spatial_locations[c("x","y")], gene_metadata = gene_ids)

# normalise then run giotto
giotto <- normalizeGiotto(gobject = giotto)
giotto = runPAGEEnrich(gobject = giotto,
                       sign_matrix = matrix,
                       min_overlap_genes = 2)

saveRDS(giotto, file.path(OUTPUT_DIR, "giotto.Rds"))




