
"%&%" <- function(a,b) paste0(a,b)

library("data.table")
library("tidyverse")

local.dir <- "/Users/jtorres/Google Drive/Science/Projects/metaxcan_t2d_v2/"
serv.dir <- "/group/im-lab/nas40t2/jason/projects/metaxcan_t2d_v2/"
gwas.dir <- serv.dir %&% "t2d/runs/gwas_conversion/results/gwas/"
gwas.file <- gwas.dir %&% "DIAGRAM_T2D_SCOTT.ma"
write.dir <- "/Users/jtorres/FUSE4/projects/metaxcan_t2d_v2/data_frames/"


gwas.df <- fread(gwas.file)

