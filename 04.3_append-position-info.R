library("data.table")
library("dplyr")

"%&%" <- function(a,b) paste0(a,b)

rescomp.dir <- "/well/mccarthy/users/jason/projects/metaxcan_t2d_v2/" 
#"/Users/jtorres/FUSE5/projects/metaxcan_t2d_v2/"
df.dir <- rescomp.dir %&% "data_frames/"
plot.dir <- rescomp.dir %&% "plots/"
ref.dir <- rescomp.dir %&% "reference_files/"
mod.dir <- rescomp.dir %&% "model_snps/GTEx-V6p-HapMap-2016-09-08/"




mod.df <- fread("cat " %&% mod.dir %&% "model-snps-v6p.txt")
snp.df <- fread("cat " %&% mod.dir %&% "snp_keyfile.txt.gz" %&% " | zmore")
full.df <- inner_join(mod.df,snp.df,by="RSID")
write.table(x=full.df,file=df.dir%&%"model-snps-v6p_full.txt",
            sep="\t",quote=FALSE,row.names=FALSE) 
