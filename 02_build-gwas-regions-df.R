
"%&%" <- function(a,b) paste0(a,b)

library("data.table")
library("dplyr")
library("purrr")
source("https://bioconductor.org/biocLite.R")
biocLite("SNPlocs.Hsapiens.dbSNP144.GRCh37")
library("SNPlocs.Hsapiens.dbSNP144.GRCh37")

local.dir <- "/Users/jtorres/Google Drive/Science/Projects/metaxcan_t2d_v2/"
serv.dir <- "/group/im-lab/nas40t2/jason/projects/metaxcan_t2d_v2/"
gwas.dir <- serv.dir %&% "t2d/runs/gwas_conversion/results/gwas/"
gwas.file <- gwas.dir %&% "DIAGRAM_T2D_SCOTT.ma"
write.dir <- serv.dir %&% "data_frames/"


gwas.df <- fread(gwas.file)

thresh <- 5 * 10^(-8)

sig.df <- filter(gwas.df,p<=thresh)
rsids <- sig.df$SNP
snps.gr <- rsidsToGRanges(rsids, caching=TRUE)

win.gr <- snps.gr
start(win.gr) <- start(win.gr) - 500000
end(win.gr) <- end(win.gr) + 500000

reduce.gr <- IRanges::reduce(win.gr)

locus.df <- data.frame(chrom=seqnames(reduce.gr),
                       loc.start=start(reduce.gr),loc.end=end(reduce.gr),
                       stringsAsFactors=FALSE)

locus.df$chrom <- gsub("ch","",locus.df$chrom) %>% as.integer(.)
locus.df$locus.id <- 1:dim(locus.df)[1]
write.table(locus.df,write.dir%&%"gwas_windows.txt",quote=FALSE,sep="\t",row.names=F)
