library(ggplot2)
library(reshape2)
library(plyr)
library(scales)
library(cowplot)
library(RColorBrewer)
library(grid)
library(plotly)

###################################################################################################
args  =  commandArgs(TRUE);
argmat  =  sapply(strsplit(args, "="), identity)

for (i in seq.int(length=ncol(argmat))) {
  assign(argmat[1, i], argmat[2, i])
}

# available variables
print(ls())


fte_theme <- function() {
  
  # Generate the colors for the chart procedurally with RColorBrewer
  palette <- brewer.pal("Greys", n=9)
  color.background = palette[2]
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
    theme(axis.text.x=element_text(size=12,color=color.axis.text, angle=45)) +
    theme(axis.text.y=element_text(size=12,color=color.axis.text)) +
    theme(axis.title.x=element_text(size=14,color=color.axis.title, vjust=0)) +
    theme(axis.title.y=element_text(size=14,color=color.axis.title, vjust=0.25)) +
    
    # Plot margins
    theme(plot.margin = unit(c(0.35, 0.2, 0.3, 0.35), "cm"))
}

#define colour blind safe palette
#define colour blind safe palette
cbPalette <- c("#375B78", "#6EF317", "#F88118", "#F731DB", "#C7F3B5", "#5D1411", "#EAAAEA", "#18702A", "#E1CE01", "#5A1463", 
               "#FD2F48", "#11AAA3", "#E8866C", "#48F89B", "#6EA8FD", "#AD1563", "#0B4FA0", "#E3EDEA", "#B562F5", "#3B3010", 
               "#1C1539", "#865116", "#D1F65E", "#FDDDAD", "#FADC68", "#6DF7E5", "#555D15", "#EB64AB", "#B7821A", "#361A23", 
               "#FEB778", "#F1697A", "#84243E", "#98F571", "#DE9CB8", "#2B8F62", "#2FD3AA", "#3088AA", "#930073", "#2977B9", 
               "#9694F9", "#C1D0ED", "#F869E5", "#771702", "#323275", "#B1E119", "#181D2D", "#368C99", "#6DF874", "#192254")
###################################################################################################
FILE=paste(FOLDER, "size_profiles-", CATEGORY, ".txt", sep="")
TABLE = read.table(FILE, header=TRUE)
nCOLS=length(TABLE)

TABLE[,2:nCOLS] = TABLE[,2:nCOLS] / 1000000

FILLEdsizes=TABLE[apply(TABLE[,-1], 1, function(x) !all(x==0)),]
MIN=min(FILLEdsizes$length)
MIN=16
MAX=max(FILLEdsizes$length)
###################################################################################################
#plot full count size profile

TABLEmelt = melt(TABLE, id.vars = c("length"))

x <- ggplot(TABLEmelt, aes(x=length, y=value, colour=variable))+
       geom_step(stat="identity")+
       labs(x="Length", y="count [*10^6]", title="size-profile-all", colour="Library")+
       scale_x_continuous(limits = c(MIN,MAX))+
       scale_colour_manual(values = cbPalette)+
       fte_theme()

outFILE=paste(PLOTdir, "compound-size_profiles-", CATEGORY, ".pdf", sep="")
ggsave(outFILE, x)
pp=ggplotly(x)
FILENAME=paste(PLOTdir, "compound-size_profiles-", CATEGORY, ".html", sep="")
htmlwidgets::saveWidget(pp, FILENAME)

###################################################################################################
#plot normalized count size profile



TABLEmelt = melt(TABLE, id.vars = c("length"))


head(TABLE)
head(TABLEmelt)

x=ggplot(TABLEmelt, aes(x=length, y=value))+
    geom_bar(stat="identity", )+
    scale_y_continuous()+
    scale_x_continuous(limits=c(MIN,MAX))+
    labs(x="Length", y="count [*10^6]", title="size-profile-all-normalized", colour="Library")+
    facet_wrap( vars(variable) )


outFILE=paste(PLOTdir, "size_profiles-", CATEGORY, ".pdf", sep="")
ggsave(outFILE, x)
pp=ggplotly(x)
FILENAME=paste(PLOTdir, "size_profiles-", CATEGORY, ".html", sep="")
htmlwidgets::saveWidget(pp, FILENAME)

