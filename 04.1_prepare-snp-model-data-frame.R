

library("data.table")
library("dplyr")

"%&%" <- function(a,b) paste0(a,b)

rescomp.dir <- "/well/mccarthy/users/jason/projects/metaxcan_t2d_v2/" 
#"/Users/jtorres/FUSE5/projects/metaxcan_t2d_v2/"
df.dir <- rescomp.dir %&% "data_frames/"
plot.dir <- rescomp.dir %&% "plots/"
ref.dir <- rescomp.dir %&% "reference_files/"
mod.dir <- rescomp.dir %&% "model_snps/GTEx-V6p-HapMap-2016-09-08/"




# This code is run on Rescomp 



"%&%" <- function(a,b) paste0(a,b)
cur.dir <- "/gpfs2/well/mccarthy/users/jason/projects/" %&% 
  "metaxcan_t2d_v2/model_snps/GTEx-V6p-HapMap-2016-09-08/"
files <- list.files("./")
files <- files[grepl(".txt.gz",files)]
out.df <- c()
pb <- txtProgressBar(min=0,max=length(files),style=3)
for (i in 1:length(files)){
  setTxtProgressBar(pb,i)
  f <- files[i]
  tiss <- (strsplit(f,split=".txt.gz")[[1]] %>% strsplit(.,split="TW_"))[[1]][2]
  sub.df <- fread("cat " %&% cur.dir%&%"TW_"%&%tiss%&%".txt.gz" %&% " | zmore")
  sub.df <- dplyr::select(sub.df,one_of("GENE","RSID2"))
  names(sub.df)[2] <- "RSID"
  sub.df <- sub.df[!duplicated(sub.df),]
  sub.df$Tissue <- tiss
  out.df <- rbind(out.df,sub.df)
}
write.table(x=out.df,file=cur.dir%&%"model-snps-v6p.txt",sep="\t",quote=FALSE,row.names=FALSE)




# Run this code in the model SNP directory 

#module load python/2.7.11
#python JTbuildSNPkey.py 




# Append Position information (Run in Rescomp)




mod.df <- fread("cat " %&% mod.dir %&% "model-snps-v6p.txt")
snp.df <- fread("cat " %&% mod.dir %&% "snp_keyfile.txt.gz" %&% " | zmore")
full.df <- inner_join(mod.df,snp.df,by="RSID")
write.table(x=full.df,file=df.dir%&%"model-snps-v6p_full.txt",
            sep="\t",quote=FALSE,row.names=FALSE) 


