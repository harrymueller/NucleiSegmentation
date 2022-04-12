# LIBRARIES
library(Seurat)
library(ggplot2)
library(patchwork)
library(raster)
library(RColorBrewer)
library(viridis)

source("/mnt/data/scripts/rscripts/functions/accurate_plot.R")

# params
library(argparser)
args <- arg_parser("Plotting nCount, nFeature, Malat1 and Neat1 violin and spatial plots.")

args <- add_argument(args, "--colour", help = " ")
args <- add_argument(args, "--extra", help = "binary | log | none", default = "none")
args <- add_argument(args, "--name", help = " ")
args <- add_argument(args, "--save_to_csv", help = " ", default = NULL)
argv <- parse_args(args)


BIN_SIZE  = c(1) #, 1)
TONGUE_ID = c("tongue-4", "tongue-5")
DIAMETER  = NULL # argv$diameter
SAVE_TO_CSV = T#!is.null(argv$save_to_csv)

#GENES = c("Krt36", "Krt84", "Krt35", "Krt6b", "Krt6a", "Dmd", "Neb", "Ttn", "Rbfox1", "C1qtnf7", "Dclk1", "Mmp16", "Mecom", "Plcb1", "Flt1", "Fetub", "Ecm1", "Fmo2", "Pappa", "Hs3st5", "St3gal4", "Mki67", "Cenpf", "Cenpe", "Frem2", "Sema3c", "Pcdh7", "Unc5c", "Cpm", "Hhip", "Hephl1") 
GENES = c("Krt4", "Krt6b", "Tchh", "Mki67")
#GENES = c("Tchh", "Krt6b", "Krt4")

# output directory
INPUT_DIR = "/mnt/data/R_analysis/gemRDS"
OUTPUT_DIR = sprintf("/mnt/data/R_analysis/count_feature_plots/gene_expression_plots/%s", argv$name)
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR)

colours = switch(argv$colour,
                 "black_magma" = magma(11),
                 "red" = c("#000000", "#FF0000"),
                 "green" = c("#000000", "#00FF00"),
                 "blue" = c("#000000", "#0000FF"),
                 "magma" = rev(magma(11)),
                 viridis(11))
  
for (id in TONGUE_ID) {
    for (bin_size in BIN_SIZE) {
        # directory
        output_dir = sprintf("%s/%s_bin%s", OUTPUT_DIR, id, bin_size)
        if (!dir.exists(output_dir)) dir.create(output_dir)

        # read in RDS
        obj = readRDS(sprintf("%s/%s_bin%s_spatialObj.rds", INPUT_DIR, id, bin_size))

        # plots
        i = 1
        
        df = NULL
        
        for (t in GENES) {
            print(sprintf("%d/%d %s", i, length(GENES), t))
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
                              filename = paste0(output_dir, "/", t, ".png"), 
                              legend_name = legend,
                              adjust = 1,
                              custom_colours = colours,
                              dpi = 750,
                              minres = 1500,
                              crop = FALSE,
                              black_background = FALSE)
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
          write.csv(df, sprintf("%s/expr.csv", output_dir), row.names = FALSE)
        }
    }
}