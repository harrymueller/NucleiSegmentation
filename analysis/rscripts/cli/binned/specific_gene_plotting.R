# LIBRARIES
library(Seurat)
library(ggplot2)
library(patchwork)
library(raster)
library(RColorBrewer)
library(viridis)

source("/mnt/data/scripts/analysis/rscripts/functions/accurate_plot.R")

# params
library(argparser)
args <- arg_parser("Plotting nCount, nFeature, Malat1 and Neat1 violin and spatial plots.")

args <- add_argument(args, "--input_rds", help = "input rds file")
args <- add_argument(args, "--output", help = "output dir")
args <- add_argument(args, "--colour", help = " ")
args <- add_argument(args, "--extra", help = "binary | log | none", default = "none")
args <- add_argument(args, "--name", help = " ")
args <- add_argument(args, "--save_to_csv", help = " ", default = NULL)
argv <- parse_args(args)

INPUT = argv$input_rds
OUTPUT_DIR = argv$output

SAVE_TO_CSV = !is.null(argv$save_to_csv)

# Read in genes from excel file
library(data.table)
library(readxl)
GENES = as.data.frame(read_excel("genes.xlsx", col_names = F, sheet = 1))[,1]
GENES = GENES[2:length(GENES)]
#print(GENES)

# output directory
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR)

colours = switch(argv$colour,
                 "black_magma" = magma(11),
                 "red" = c("#000000", "#FF0000"),
                 "green" = c("#000000", "#00FF00"),
                 "blue" = c("#000000", "#0000FF"),
                 "magma" = rev(magma(11)),
                 viridis(11))
  
# read in RDS
obj = readRDS(INPUT)

# plots
i = 1

df = NULL

mask = (GENES %in% rownames(obj))
this_genes = GENES[mask]

for (t in this_genes) {
    print(sprintf("%d/%d %s", i, length(this_genes), t))
    i = i + 1;

    data = SpatialPlot(obj, features = t)$data
    if (argv$extra == "log") {
      data[,3] = log(data[,3] + 1)
      legend = sprintf("log (%s + 1)", t)
    } else if (argv$extra == "binary") {
      data[,3] = ifelse(data[,3] > 0, 1, 0)
      legend = sprintf("binary %s", t)
    } else {
      legend = t
    }
    
    if (!SAVE_TO_CSV) {
      accurate_plot(data, 
                      filename = paste0(OUTPUT_DIR, "/", t, ".png"), 
                      legend_name = legend,
                      adjust = 1,
                      custom_colours = colours,
                      dpi = 750,
                      minres = 1500,
                      crop = FALSE,
                      black_background = FALSE,
                      use_discrete_colours = TRUE)
    } else {
      if (is.null(df)) {
        df = data
      } else {
        df = cbind(df, data[,3])
        colnames(df)[length(colnames(df))] = t
      }
    }
}
if (SAVE_TO_CSV) {
  write.csv(df, sprintf("%s/expr.csv", OUTPUT_DIR), row.names = FALSE)
}