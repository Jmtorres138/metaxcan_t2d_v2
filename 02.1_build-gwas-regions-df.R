
"%&%" <- function(a,b) paste0(a,b)

library("data.table")
library("dplyr")
library("purrr")
library("GenomicRanges")


fuse.dir <- "/home/jason/science/servers/"
rescomp.dir <- fuse.dir %&% "FUSE5/projects/metaxcan_t2d_v2/"
data.dir <- rescomp.dir %&% "data_frames/"
gwas.file <- data.dir %&% "Scott_T2D_DIAGRAM-full.txt.gz"


gwas.df <- fread("cat " %&% gwas.file %&% " | zmore")

thresh <- 5 * 10^(-8)

sig.df <- filter(gwas.df,P<=thresh)

win.gr <- GRanges(seqnames=sig.df$CHR,IRanges(start=sig.df$POS,end=sig.df$POS))

start(win.gr) <- start(win.gr) - 1e6 
end(win.gr) <- end(win.gr) + 1e6 

reduce.gr <- IRanges::reduce(win.gr)

locus.df <- data.frame(chrom=seqnames(reduce.gr),
                       loc.start=start(reduce.gr),loc.end=end(reduce.gr),
                       stringsAsFactors=FALSE)

locus.df$chrom <- gsub("ch","",locus.df$chrom) %>% as.integer(.)
locus.df$locus.id <- 1:dim(locus.df)[1]
write.table(locus.df,data.dir%&%"gwas_windows.txt",quote=FALSE,sep="\t",row.names=F)
