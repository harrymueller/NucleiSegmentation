# Create summarised experiment RDS from seurat obj
suppressMessages(library(Seurat))
suppressMessages(library(S4Vectors))
suppressMessages(library(SummarizedExperiment))

INPUT = "/mnt/data/R_analysis/cell_reference/1.0_GlobalWithEpithelialRes0.6.rds"
OUTPUT_PATH = "/mnt/data/R_analysis/cell_reference"

rds = readRDS(INPUT)
Idents(rds) = "ClusterInclEpith06" 
DefaultAssay(rds) = "RNA"

# column data
colData = data.frame(
  "label.main" = unfactor(unique(Idents(rds))),
  "label.fine" = unfactor(unique(Idents(rds)))
)

for (n in c("Zero", "One", "Two", "Three", "Four", "Five", "Six", "Seven", "Eight", "Nine", "Ten", "Eleven", "Twelve", "Thirteen", "Fourteen")) {
  colData$label.main[colData$label.main == n] = "Epith"
}

# SCT
rds = SCTransform(rds)

SCTData = AverageExpression(rds, assays = "SCT")
SCTData = as.matrix(SCTData$SCT)

sct = SummarizedExperiment(assays = list("logcounts" = SCTData), rowData = rownames(SCTData), colData = colData)
saveRDS(sct, paste0(OUTPUT_PATH, "/cell_ref_SCT.rds"))

# LN
DefaultAssay(rds) = "RNA"
rds = NormalizeData(rds, normalization.method = "LogNormalize")

LNData = AverageExpression(rds, assays = "RNA")
LNData = as.matrix(LNData$RNA)

ln = SummarizedExperiment(assays = list("logcounts" = LNData), rowData = rownames(LNData), colData = colData)
saveRDS(ln, paste0(OUTPUT_PATH, "/cell_ref_LN.rds"))
