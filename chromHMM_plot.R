

"%&%" <- function(a,b) paste0(a,b)

library("data.table")
library("tidyverse")
library("ggbio")

got2d.dir <- "/Users/jtorres/FUSE/" 
states.dir <- got2d.dir %&% "reference/chromatin_segmentation/varshney_2016/chromatin_states/"
rescomp.dir <- "/Users/jtorres/FUSE5/projects/metaxcan_t2d_v2/"
df.dir <- rescomp.dir %&% "data_frames/"


save_t2d_tissue_states <- function(){
    isl.file <- states.dir %&% "Islets.chromatinStates.bed.gz"
    isl.df <- fread("cat " %&% isl.file %&% " | zmore",sep="\t")
    isl.df$V6 <- "Islets"
    
    adi.file <- states.dir %&% "Adipose.chromatinStates.bed.gz"
    adi.df <- fread("cat " %&% adi.file %&% " | zmore",sep="\t")
    adi.df$V6 <- "Adipose"
    
    liv.file <- states.dir %&% "Liver.chromatinStates.bed.gz"
    liv.df <- fread("cat " %&% liv.file %&% " | zmore",sep="\t")
    liv.df$V6 <- "Liver"
    
    mus.file <- states.dir %&% "SkeletalMuscle.chromatinStates.bed.gz"
    mus.df <- fread("cat " %&% mus.file %&% " | zmore",sep="\t")
    mus.df$V6 <- "Muscle"
    
    state.df <- rbind(isl.df,adi.df,liv.df,mus.df)
    
    state.df$V2 <- state.df$V2 + 1 ; state.df$V3 <- state.df$V3 + 1 
    names(state.df) <- c("CHR","START","END","STATE","RGB","Tissue")
    pb <- txtProgressBar(min=0,max=dim(state.df)[1],style=3)
    state.df$COL <- map(1:length(state.df$RGB),function(i){
      setTxtProgressBar(pb,i)
      rgb <- state.df$RGB[i]
      vec <- strsplit(rgb,split=",")[[1]]
      r <- vec[1]; g <- vec[2]; b <- vec[3]
      col <- rgb(r,g,b,maxColorValue=255)
    }) %>% as.character(.)  
    write.table(state.df,file=df.dir%&%"chromHMM_t2d-tisues.txt",sep="\t",row.names=FALSE,quote=FALSE)
}


chromHMM_plot <- function(state.df,chrom,loc.start,loc.end){
  sub.df <- filter(state.df,CHR==chrom,START<=loc.end,END>=loc.start)
  sub.df$START <- sub.df$START %>% as.integer(.)
  sub.df$END <- sub.df$END %>% as.integer(.)
  sub.df$Y <- ifelse(sub.df$Tissue=="Islets",4,
                     ifelse(sub.df$Tissue=="Liver",3,
                            ifelse(sub.df$Tissue=="Muscle",2,
                                   ifelse(sub.df$Tissue=="Adipose",1,NA))))
  
  
  plt <- ggplot(data=sub.df) + 
    geom_segment(data=sub.df,aes(x=START,xend=END,y=Y,yend=Y),color=sub.df$COL,size=10) + 
    coord_cartesian(ylim=c(0,4.5),expand = FALSE) + 
    theme_clear() + 
    theme(axis.title = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y=element_blank()) 
  
  plt <- ggplot(data=sub.df) + 
    geom_rect(data=sub.df,aes(xmin=START,xmax=END,ymin=Y,ymax=Y+1),color=sub.df$COL,fill=sub.df$COL) + 
    #coord_cartesian(ylim=c(-1,5),expand = FALSE) + 
    theme_clear() + 
    theme(axis.title = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y=element_blank()) 
  return(plt)
}



  
