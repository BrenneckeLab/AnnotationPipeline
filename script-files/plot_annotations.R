library(ggplot2)
library(cowplot)
library(reshape2)
library(plyr)
library(RColorBrewer)
library(plotly)


###################################################################################################
args  =  commandArgs(TRUE);
argmat  =  sapply(strsplit(args, "="), identity)

for (i in seq.int(length=ncol(argmat))) {
  assign(argmat[1, i], argmat[2, i])
}

# available variables
print(ls())

###################################################################################################

setwd( TMP)

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

#define colour blind safe palette
cbPalette <- rev(c("#000000", "#004949", "#DBD100", "#009292", "#FF6DB6", "#FFB6DB", "#490092", "#006DDB", "#B66DFF", "#6DB6FF", "#B6DBFF", "#920000", "#24FF24", "#924900", "#FFFF6D"))

###################################################################################################

TABLE = read.table(FILE, header=TRUE)
nCOLS=length(TABLE)

TABLE$sort = as.numeric(rownames(TABLE))
TABLE$annotation <- reorder(TABLE$annotation, TABLE$sort)

TABLE$annotation <- factor(TABLE$annotation, 
                          levels=TABLE[order(TABLE$sort),]$annotation)


TABLE = subset(TABLE, select=-c(sort))

#--------------------------------------------------------------------------------------------------
#plot filter stats
FILTERcategories = c("adaptor_dimer", "artifact_filtered", "length_filtered_short", "length_filtered_long", "N_filtered", "reads_after_filtering")

TABLEsub = subset(TABLE, annotation %in% FILTERcategories)
TABLEmelt = melt(TABLEsub, idvar = c("annotation"))

#subset colors to get nice plot
PAL = cbPalette[c(9,13,6,14,2,4)]

x <- ggplot(TABLEmelt, aes(x=rev(variable), y=rev(value), text=paste("Category=", rev(annotation), "<br>", "Value=", rev(value) , "<br>", "Sample=", rev(variable)), fill=rev(annotation)))+
  geom_bar(position= "fill",  stat="identity", width=.3)+
  scale_y_continuous(labels = scales::percent)+
  scale_fill_manual(values = PAL)+
  labs( fill="Annotation", x="Library", y="relative fraction")+
  fte_theme()

FILENAME=paste(PLOTdir, "1_annotation-filtered_reads.pdf", sep="")
ggsave(FILENAME, x, height=10, width=2+0.5*nCOLS )
pp=ggplotly(x, tooltip = c("text"))
FILENAME=paste(PLOTdir, "1_annotation-filtered_reads.html", sep="")
htmlwidgets::saveWidget(pp, FILENAME)

#--------------------------------------------------------------------------------------------------
#plot mapping stats
#determine read number after filters were applied
TOTALreads=TABLE[TABLE$annotation == "reads_after_filtering",2:nCOLS]

#filter for all reads not in the filter categories and calculate
##number of reads mapping ot the genome
FILTERcategories = c("adaptor_dimer", "artifact_filtered", "length_filtered_short","length_filtered_long", "N_filtered", "reads_after_filtering")
  TABLEsub = subset(TABLE, !( annotation %in% FILTERcategories))


#extract crap categories and calculate sum 
FILTERcategories = c("reads_not_mapped","rRNA", "mito", "rRNA_AS", "mito_AS", "tRNA", "tRNA_AS")
TABLEsub = subset(TABLE, annotation %in% FILTERcategories)

if (nCOLS == 2){
  CRAPreads=sum(TABLEsub[,2:nCOLS])
}else {
  CRAPreads = colSums(TABLEsub[,2:nCOLS])
}

#calculate number of good-mapping reads
GOODreads = TABLE[2,]
GOODreads[,1] = "good_reads"
GOODreads[,2:nCOLS] = TOTALreads - CRAPreads

#combine data frames and melt for plotting
TABLEsub = rbind(TABLEsub, GOODreads)
TABLEmelt = melt(TABLEsub, idvar = c("annotation"))

PAL = cbPalette[c(1,10,15,12,9,13,6,14,2)]

x <- ggplot(TABLEmelt, aes(x=variable, y=value, text=paste("Category=", annotation, "<br>", "Value=", value , "<br>", "Sample=", variable), fill=annotation))+
  geom_bar(position= "fill",  stat="identity", width=.3)+
  scale_y_continuous(labels = scales::percent)+
  scale_fill_manual(values = PAL)+
  labs( fill="Annotation", x="Library", y="relative fraction")+
  fte_theme()

FILENAME=paste(PLOTdir, "2_annotation-crap-reads.pdf", sep="")
ggsave(FILENAME, x, height=10, width=2+0.5*nCOLS )
pp=ggplotly(x, tooltip = c("text"))
FILENAME=paste(PLOTdir, "2_annotation-crap-reads.html", sep="")
htmlwidgets::saveWidget(pp, FILENAME)


#--------------------------------------------------------------------------------------------------
#plot annotation categories

FILTERcategories = c("adaptor_dimer", "artifact_filtered", "length_filtered_short", "length_filtered_long", "N_filtered", "reads_after_filtering","reads_not_mapped","rRNA","mito", 
                     "rRNA_AS", "mito_AS", "tRNA", "tRNA_AS")

TABLEsub = subset(TABLE, !( annotation %in% FILTERcategories))
TABLEmelt = rev(melt(TABLEsub, id.vars = c("annotation")))

PAL = cbPalette[c(2,11,3,5,4,7,8,1,10,15,12,9,13,6,14)]

x <- ggplot(TABLEmelt, aes(x=rev(variable), y=rev(value), text=paste("Category=", rev(annotation), "<br>", "Value=", rev(value) , "<br>", "Sample=", rev(variable)), fill=rev(annotation)))+
  geom_bar(position= "fill",  stat="identity", width=.3)+
  scale_y_continuous(labels = scales::percent)+
  scale_fill_manual(values = PAL)+
  labs( fill="Annotation", x="Library", y="relative fraction")+
  fte_theme()

FILENAME=paste(PLOTdir, "3_annotation-split-up.pdf", sep="")
ggsave(FILENAME, x, height=10, width=2+0.5*nCOLS )
pp=ggplotly(x, tooltip = c("text"))
FILENAME=paste(PLOTdir, "3_annotation-split-up.html", sep="")
htmlwidgets::saveWidget(pp, FILENAME)
