args  =  commandArgs(TRUE);
argmat  =  sapply(strsplit(args, "="), identity)
print("hi")
for (i in seq.int(length=ncol(argmat))) {
  assign(argmat[1, i], argmat[2, i])
}

# available variables
print(ls())

###################################################################################################

library(tidyverse)
library(cowplot)
library(ggforce)
library(plotly)

#---------------------------------------------------------------------------------------------------------




if( TYPE == "CHIPseq" || TYPE == "DNAseq" ){
  INFILE=paste(INPUT, "_all.man",STRINGENCYext, sep="")
  ALL=read_tsv(INFILE, col_names = TRUE, cols(TE=col_character())) 
  
  nLIBs=length(ALL)-2

  ALL=ALL %>%
    gather(key = "LIB", value = "VALUE", c(3:length(ALL)))

  nROW=3
  nCOL=1
  nPAGES=ceiling(nLIBs/nROW/nCOL)
    
  library(doParallel)
  registerDoParallel(cores = 2)

  foreach(i = unique(ALL$TE), .packages = c('ggplot2', 'dplyr', 'plotly', 'htmlwidgets', 'cowplot', 'ggforce')) %dopar% {
    print(i)
    
    FILENAME=paste(OPENdir, i, ".pdf", sep="")
    
    pdf(FILENAME, width=12, height=6)
    for (PAGE in seq(1,nPAGES)){
        
      p=ggplot( ALL, mapping = aes(x=POS, y=VALUE)) + 
        geom_line( data=filter(ALL,TE==i))+
              facet_wrap( ~LIB, ncol=1)+
              xlab(paste("position on ", i, sep=""))+
              ylab("normalized read-counts")+
              facet_wrap_paginate(~ LIB, ncol = nCOL, nrow = nROW, page = PAGE)+
              theme_cowplot(11)+
              theme(axis.text.x = element_text(size=8),
                    axis.text.y = element_text(size=8),
                    text = element_text(size=8),
                    strip.text = element_text(size = 10))
                  
      print(p)  
    }
    dev.off()  
    
    # p=ggplot( ALL, mapping = aes(x=POS, y=VALUE,color=LIB)) + 
    #   geom_line( data=filter(ALL,TE==i))+
    #         xlab(paste("position on ", i, sep=""))+
    #         ylab("normalized read-counts")+
    #         labs(title = paste(i,"\n samples can be hidden by clicking/double-clicking in the legend",sep="") )+
    #         theme_cowplot(11)+
    #         theme(axis.text.x = element_text(size=8),
    #               axis.text.y = element_text(size=8),
    #               text = element_text(size=8),
    #               strip.text = element_text(size = 10))
   
    # pp = ggplotly(p)
    
    # FILENAME=paste(OPENdir, i, ".html", sep="")
    # htmlwidgets::saveWidget(pp, FILENAME, selfcontained = FALSE, libdir = paste(OPENdir,"/html-dependencies/",i, sep=""))

    gc()

  }

  stopImplicitCluster()
    
  
  
}else{
  INFILE=paste(INPUT, "_sense.man",STRINGENCYext, sep="")
  SENSE=read_tsv(INFILE, col_names = TRUE, cols(TE=col_character())) 

  nLIBs=length(SENSE)-2

  SENSE=SENSE %>%
    gather(key = "LIB", value = "VALUE", c(3:length(SENSE)))
  
  if( TYPE == "sRNAseq" || TYPE == "sRNAseqIP" ){
    INFILE=paste(INPUT, "_sense.man.21mer",STRINGENCYext, sep="")
    SENSE_21mer=read_tsv(INFILE, col_names = TRUE, cols(TE=col_character())) 
    SENSE_21mer = SENSE_21mer %>%
      gather(key = "LIB", value = "VALUE", c(3:length(SENSE_21mer)))
  }


  INFILE=paste(INPUT, "_antisense.man",STRINGENCYext, sep="")
  ANTISENSE=read_tsv(INFILE, col_names = TRUE, cols(TE=col_character()))
  ANTISENSE=ANTISENSE %>%
    gather(key = "LIB", value = "VALUE", c(3:length(ANTISENSE)))

  if( TYPE == "sRNAseq" || TYPE == "sRNAseqIP" ){
    INFILE=paste(INPUT, "_antisense.man.21mer",STRINGENCYext, sep="")
    ANTISENSE_21mer=read_tsv(INFILE, col_names = TRUE, cols(TE=col_character()))
    
    ANTISENSE_21mer = ANTISENSE_21mer %>%
      gather(key = "LIB", value = "VALUE", c(3:length(ANTISENSE_21mer)))
  }

  
  nROW=3
  nCOL=1
  nPAGES=ceiling(nLIBs/nROW/nCOL)

  library(doParallel)
  registerDoParallel(cores = 2)

  foreach(i = unique(SENSE$TE), .packages = c('ggplot2', 'dplyr', 'plotly', 'htmlwidgets', 'cowplot', 'ggforce')) %dopar% {
    print(i)
    FILENAME=paste(OPENdir, i, ".pdf", sep="")
    
    pdf(FILENAME, width=12, height=6)
    for (PAGE in seq(1,nPAGES)){
        nam <- paste("p", PAGE, sep = "")
        
        p=ggplot(SENSE, mapping = aes(x=POS, y=VALUE)) + 
            geom_line(data=filter(SENSE, TE==i)) +
            geom_line(data=filter(ANTISENSE, TE==i)) +
            geom_hline(yintercept=0) +
            facet_wrap(~LIB, ncol=1) +
            xlab(paste("position on ", i, sep="")) +
            ylab("normalized read-counts") +
            facet_wrap_paginate(~ LIB, ncol = nCOL, nrow = nROW, page = PAGE) +  
            theme_cowplot(11)+
            theme(axis.text.x = element_text(size=8),
                  axis.text.y = element_text(size=8),
                  text = element_text(size=8),
                  strip.text = element_text(size = 10))

        if (TYPE == "sRNAseq" || TYPE == "sRNAseqIP") {
          p = p + 
            geom_line(data=filter(SENSE_21mer, TE==i), color="orange") +
            geom_line(data=filter(ANTISENSE_21mer, TE==i), color="orange")
        }
        
        print(p)  
    }
    dev.off()

    # p = ggplot(SENSE, mapping = aes(x=POS, y=VALUE, color=LIB)) + 
    #   geom_line(data=filter(SENSE, TE==i)) +
    #   geom_line(data=filter(ANTISENSE, TE==i)) +
    #   geom_hline(yintercept=0) +
    #   xlab(paste("position on ", i, sep="")) +
    #   ylab("normalized read-counts") +
    #   labs(title = paste(i,"\n samples can be hidden by clicking/double-clicking in the legend",sep="") )+
    #   theme_cowplot(11)+
    #   theme(axis.text.x = element_text(size=8),
    #         axis.text.y = element_text(size=8),
    #         text = element_text(size=8),
    #         strip.text = element_text(size = 10))

    #   if (TYPE == "sRNAseq" || TYPE == "sRNAseqIP") {
    #     p = p + 
    #       geom_line(data=filter(SENSE_21mer, TE==i), aes(linetype = "21mer") )+
    #       geom_line(data=filter(ANTISENSE_21mer, TE==i), aes(linetype = "21mer"))+
    #       scale_linetype_manual(values = c("21mer" = "dotted")) 
    #       # guides(
    #       #   color = guide_legend(order = 1),
    #       #   linetype = guide_legend(order = 2)
    #       # )
      
    #   }

    # pp = ggplotly(p)
    
    # FILENAME = paste(OPENdir, i, ".html", sep="")
    # htmlwidgets::saveWidget(pp, FILENAME, selfcontained = FALSE, libdir = paste(OPENdir,"/html-dependencies/",i, sep=""))

    gc()
  }

  stopImplicitCluster()
  
}