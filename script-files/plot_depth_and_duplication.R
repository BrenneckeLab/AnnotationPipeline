library(ggplot2)
library(reshape2)
library(plyr)
library(RColorBrewer)
library(grid)

###################################################################################################
args  =  commandArgs(TRUE);
argmat  =  sapply(strsplit(args, "="), identity)

for (i in seq.int(length=ncol(argmat))) {
  assign(argmat[1, i], argmat[2, i])
}

# available variables
print(ls())

###################################################################################################
fte_theme <- function() {
  
  # Generate the colors for the chart procedurally with RColorBrewer
  palette <- brewer.pal("Greys", n=9)
  color.background = palette[3]
  color.grid.major = palette[4]
  color.axis.text = palette[7]
  color.axis.title = palette[8]
  color.title = palette[9]
  
  # Begin construction of chart
  theme_bw(base_size=9) +
    
    # Set the entire chart region to a light gray color
    theme(panel.background=element_rect(fill=color.background, color=color.background)) +
    theme(plot.background=element_rect(fill=color.background, color=color.background)) +
    theme(panel.border=element_rect(color=color.background)) +
    
    # Format the grid
    theme(panel.grid.major=element_line(color=color.grid.major,size=.25)) +
    theme(panel.grid.minor=element_blank()) +
    theme(axis.ticks=element_blank()) +
    
    # Format the legend, but hide by default
    #theme(legend.title=element_blank()) +
    theme(legend.background = element_rect(fill=color.background)) +
    theme(legend.text = element_text(size=7,color=color.axis.title)) +
    
    # Set title and axis labels, and format these and tick marks
    theme(plot.title=element_text(color=color.title, size=15, vjust=1.25)) +
    theme(axis.text.x=element_text(size=10,color=color.axis.text, angle=90)) +
    theme(axis.text.y=element_text(size=12,color=color.axis.text)) +
    theme(axis.title.x=element_text(size=14,color=color.axis.title, vjust=0)) +
    theme(axis.title.y=element_text(size=14,color=color.axis.title, vjust=0.5)) +
    
    # Plot margins
    theme(plot.margin = unit(c(0.35, 0.2, 0.3, 0.35), "cm"))
}

roundUp <- function(x) 10^ceiling(log10(x))

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
PAL = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

#PAL= c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black", "gold1", "skyblue2", "palegreen2", "#FDBF6F", "gray70", "maroon", "orchid1", "darkturquoise", "darkorange4", "brown")
###################################################################################################
setwd( TMP)
TABLE = read.table(FILE, header=TRUE)

TABLEmelt = melt(TABLE, idvar = c("annotation"))

MAX = roundUp(max(TABLEmelt$value))

###################################################################################################

p=ggplot ( TABLEmelt, aes(x=annotation, y=value, colour=variable))+
    fte_theme()+
    #theme_bw()+
    theme(panel.grid.major = element_line(colour = "black", linetype = "dotted"), axis.text.x= element_text(angle=90, hjust = 1, vjust=0.5))+
    geom_point(  stat = "identity", alpha=0.7 )+
    scale_y_log10(limits=c(1,MAX),  
                  breaks = c(1,10,100,1000),
                  labels = c(1,10,100,1000)
                  )+
    scale_color_manual(values = PAL)+
    annotation_logticks(sides = "l")+
    labs(y="sequencing-depth", colour="library")


if (CATEGORY == "sequencing-depth" ){
  NAME=paste("4_",CATEGORY,".pdf",sep="") 
}else{
  NAME=paste("5_",CATEGORY,".pdf",sep="")
}

ggsave(NAME, p, height=10, width=10 )

