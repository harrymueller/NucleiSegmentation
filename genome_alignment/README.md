GitHub: https://github.com/BGIResearch/SAW

Used `saw_5.sh` to run. Change params at top. Uses SAW 4.1.0, can update for latest SAW.

## Reference genome
Reference genome already indexed - available in colab server at /data/Rui_transfer/ref_genome. Gene list available at /data/Rui_transfer/

Or you can index it again as described in the GitHub. Use Gencode reference genome.

## Directory structure
.
├── images -> contains image files after qc = *.json and *.tar.gz
├── mask   -> *.h5 mask file
├── output -> outputs created here
└── reads  -> read files *_read_*.fz.gz

