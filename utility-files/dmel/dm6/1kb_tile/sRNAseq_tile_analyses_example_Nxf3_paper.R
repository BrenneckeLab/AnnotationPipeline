################################################################################
# sRNAseq_tile_analyses_example_Nxf3_paper.R
# Author: Peter Refsing Andersen; Contact: pra@mbg.au.dk

# This code reproduces the analyses and plots of Figure 2B-C and 4F in
# PMID: 31398345; DOI: 10.1016/j.cell.2019.07.007

################################################################################
# SCRIPT OVERVIEW:
# 1. Read in and organize data
# 2. Ratio calculation and master data frame generation
# 3. Make box plots of clusters and Other RD-SL
# 4. Make whole-genome scatter plots with cluster highlights

################################################################################
### Set up working environment:

# Read in required packages
library(tidyverse)

# Set up folder paths
work.dir <- paste(".../sRNAseq_example/", sep="") # <-- Set local path!
data.dir <- paste(work.dir, "data/", sep = "")

# Define the directory for output files
output.dir <- paste(work.dir, "Output_from_R_", Sys.Date(), "/", sep="") 
# Create the output directory in work.dir. Will ignore if it already exists
dir.create(file.path(output.dir), showWarnings = FALSE) 

dataset.shortname <- "Nxf3_sRNAseq_KO_final"

# Set path of the window annotation file
tile.anno.file <- paste(work.dir, "kb_tile_files/",
                        "kbWindows_mainchr_genes_tss_dm6.txt", sep="") 

# Find data files:
# Make a list of the txt files found in the input folder
data.files <- list.files(data.dir, "WindowCounts_*") 
# Make a list of the txt file strand names
strand.names <- sub("\\..*", "", sub(".*_", "", data.files)) 

# Set up color blind-friendly color palettes:
# Source: http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
cbgPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                "#0072B2", "#D55E00", "#CC79A7")
# pie(rep(1,8), col = cbgPalette) # prints a pie chart to see the colors


################################################################################
### 1. Read in and organize data:

# Read in window annotation data:
tile.df <- read.table(tile.anno.file)
# rhino.tile.df <- read.csv(rhino.tile.file)

# Set column names for the window data frame
names(tile.df) <- c("chr", "start", "end", "uniq_plus", "uniq_minus", 
                    "all_plus", "all_minus", "window_name", "transcripts", 
                    "genes", "transcript_TSSs", "gene_TSSs") 

tile.df <- tile.df %>%
  mutate(av.uniq.mapppability = (uniq_plus+uniq_minus)/2)

# Create/clear needed variables:
tile.reads <- list()
strand.id <- vector()


# Read in the data from each tile file:
i=1
for(i in 1:length(data.files)) { # Loop for every txt file
  data.file.temp <- paste(data.dir, data.files[i], sep="") # Path of txt file
  # read data from each file:
  data.head <- read_lines(data.file.temp, n_max = 3) # Get headers input file
  if (i == 1) {
    samples <- unlist(strsplit(data.head[1],"\t")) # Split sample name line
    norm.reads <- as.numeric(strsplit(data.head[2],"\t")[[1]]) # Norm. reads
  }
  # Split the total window read number line into substrings
  tile.reads[[strand.names[i]]] <- as.numeric(strsplit(data.head[3],"\t")[[1]]) 
  ifelse(i == 1,
         tile.reads[["Total"]] <- tile.reads[strand.names[i]][[1]],
         tile.reads[["Total"]] <- (tile.reads[["Total"]] + 
                                     tile.reads[strand.names[i]][[1]]))
  # Read and store data from specific strand
  assign(paste("Data", strand.names[i], "df", sep="."), 
         read.table(data.file.temp, sep="\t", skip=3)) 
}

Data.Both.df <- Data.Minus.df + Data.Plus.df # Sum data from both strands
names(Data.Both.df) <- samples

################################################################################
### 2. Ratio calculation and master data frame generation:

# Normalize to library depth (miRNA counts from pipeline)
norm.values <- norm.reads # Set normalization number to miRNA counts
# Normalize data to library depth
norm.data.df <- as.data.frame(t(t(Data.Both.df) / (norm.values))) 

# Note: in the following analyses we focus on pooled duplicate samples since
# analyses of the individual replicates looked almost identical. The single 
# experiments are, however, still in the data frame to look at.

# Create log2-based fold change data frame relative to control
ref.samp <- grep("w1118_dupl", names(norm.data.df)) # Find samples by name
rhi.samp <- grep("rhino_TH1_dupl", names(norm.data.df))
pseudo.count <- 1 # Set pseudo-count value
norm.log2fc.df <- log2((norm.data.df[, -ref.samp] + pseudo.count) / 
                         (norm.data.df[, ref.samp] + pseudo.count))

# Add "L2FC" to the sample names
names(norm.log2fc.df) <- paste("L2FC", names(norm.log2fc.df), sep="_")

### Assemble calculated data and add mappability average values:
# calculate the average mappability of the + and - strand values for each tile
MapNormVector <- (tile.df$uniq_plus+tile.df$uniq_minus)/2

# Normalize the data to the average mappability index
map.norm.data.df <- norm.data.df / (MapNormVector/100) 

# Add MN (Mappability Normalized) to the sample names
names(map.norm.data.df) <- paste("MN_", names(norm.data.df), sep="") 
 
# Make master table with all normalized data:
AllData <- cbind(norm.data.df, map.norm.data.df, norm.log2fc.df)
SummedMaster <- cbind(tile.df, AllData)

# Annotate main germ-line clusters in a category column:
min.map.filter <- 10 # Set minimum mappability percentage
SummedMaster$Category <- "Other"
SummedMaster$Category[SummedMaster$transcripts != "."] <- "Gene body"
SummedMaster$Category[SummedMaster$transcript_TSSs != "."] <- "Gene TSS"
SummedMaster$Category[grepl("42AB", SummedMaster$transcripts)] <- "Cluster42AB"
SummedMaster$Category[grepl("80F", SummedMaster$transcripts)] <- "Cluster80F"
SummedMaster$Category[grepl("20A", SummedMaster$transcripts)] <- "Cluster20A"
SummedMaster$Category[grepl("38C1", SummedMaster$transcripts)] <- "Cluster38C1"
SummedMaster$Category[grepl("38C2", SummedMaster$transcripts)] <- "Cluster38C2"
SummedMaster$Category[grepl("Flam", SummedMaster$transcripts)] <- "ClusterFlamenco"
FlamStart <- min(SummedMaster$start[SummedMaster$Category == "ClusterFlamenco"])
SummedMaster$Category[SummedMaster$Category == "ClusterFlamenco" & SummedMaster$end > FlamStart + 50000] <- "filtered"
SummedMaster$Category[SummedMaster$av.uniq.mapppability < min.map.filter] <- "filtered"

# inspection_set <- SummedMaster[SummedMaster$Category == "ClusterFlamenco", c("chr", "start", "end", "transcripts", "av.uniq.mapppability")]

# Reformat the Category column to 'factor' type and set the levels
cat.levels <- c("Other","Gene body", "Gene TSS", "Cluster42AB", "Cluster80F", 
                "Cluster38C1", "Cluster38C2", "ClusterFlamenco", "Cluster20A",
                "filtered")
SummedMaster$Category <- factor(SummedMaster$Category, levels=cat.levels) 

# Annotate Rhino-dependent source loci (SL):
MIN_COUNT <- 250 # Set threshold for calling of piRNA source loci in general
MIN_FC <- 10 # Set threshold for calling piRNA source loci 'Rhino-dependent'


# Categorize tiles based on piRNA source locus character:
SummedMaster <- SummedMaster %>%
  # Calculate fold-change values for Rhino KO relative to w1118
  mutate(rhiKO.foldchange = (sRNAseq_w1118_dupl_81504_81505 + 1) / (sRNAseq_rhino_TH1_dupl_81510_81511 + 1)) %>%
  # Define source loci categories (SL_CATS) based on rhino-dependent piRNA production:
  mutate(SL_CATS = ifelse(av.uniq.mapppability < min.map.filter, "filtered",
                          ifelse(MN_sRNAseq_w1118_dupl_81504_81505 > MIN_COUNT & rhiKO.foldchange > MIN_FC, "RD-SL",
                                 ifelse(MN_sRNAseq_w1118_dupl_81504_81505 > MIN_COUNT & rhiKO.foldchange < 2, "RI-SL", 
                                        ifelse(MN_sRNAseq_w1118_dupl_81504_81505 > MIN_COUNT & rhiKO.foldchange <= MIN_FC, "RsemiD-SL", "Non-SL")
                                        )))) %>%
  # Set factor levels for the SL_CATS variable:
  mutate(SL_CATS = factor(SL_CATS, levels=c("Non-SL", "RI-SL", "RsemiD-SL", "RD-SL", "filtered")))


### Aggregate data for summary statistics:
Cat_SanityCheck <- SummedMaster %>% 
  count(Category, SL_CATS) # Sanity check of cluster categories

# Save the sRNAseq table for use in other NGS analyses
sRNAseq.df.file <- paste(output.dir, Sys.Date(), "_sRNAseq_Master_Table.csv", sep="")
# write.csv(SummedMaster, file = sRNAseq.df.file) # To save the processed data

# Aggregate to get total category counts for each sample:
aggr.counts <- aggregate(select_if(SummedMaster, grepl("^sRNAseq", names(SummedMaster))),
                         by=SummedMaster["Category"], FUN = sum)

AGC.dupl <- aggr.counts %>%
  select(one_of(c("Category", "sRNAseq_w1118_dupl_81504_81505", 
                  "sRNAseq_rhino_TH1_dupl_81510_81511",
                  "sRNAseq_nxf3_TH2_dupl_81508_81509")))

AGC.dupl.percent <- cbind(AGC.dupl[,"Category"], AGC.dupl[, -1] / 
                            AGC.dupl[, "sRNAseq_w1118_dupl_81504_81505"] * 100)


################################################################################
### 3. Make box plots of clusters and Other RD-SL

# Select categories to apply cluster annotations:
category.sel <- c("ClusterFlamenco", "Cluster20A", "Cluster38C1",
                  "Cluster38C2", "Cluster80F", "Cluster42AB", "Other RD-SL") 

# Set order of cluster categories + add "Other RD-SL" factor:
category.order <- c("ClusterFlamenco", "Cluster20A", "Cluster80F",
                     "Cluster42AB", "Cluster38C1/2", "Other RD-SL") 

# Set column selection for reduced data frame for gathering:
plot.cols <- c("SL.groups", "L2FC_sRNAseq_nxf3_TH2_dupl_81508_81509", "L2FC_sRNAseq_rhino_TH1_dupl_81510_81511")

# Set genotype order for correct plot output:
genotype.order <- c("nxf3_TH2_dupl", "rhino_TH1_dupl")

# Organize data for plotting:
plotL2FC.df <- SummedMaster %>%
  mutate(Category = factor(Category, levels = c(levels(SummedMaster$Category), "Other RD-SL"))) %>%
  mutate(SL.groups = replace(Category, SL_CATS == "RD-SL" & !Category %in% category.sel, "Other RD-SL")) %>%
  filter(SL.groups %in% category.sel) %>%
  filter(MN_sRNAseq_w1118_dupl_81504_81505 > MIN_COUNT) %>%
  select(one_of(plot.cols)) %>%
  gather(key = "Genotype", value = value, -SL.groups) %>%
  mutate(SL.groups = recode(SL.groups,  "Cluster38C1" = "Cluster38C1/2", "Cluster38C2" = "Cluster38C1/2")) %>%
  mutate(SL.groups = factor(SL.groups, levels = category.order)) %>%
  mutate(Genotype = recode(Genotype,
                           "L2FC_sRNAseq_nxf3_TH2_dupl_81508_81509" = "nxf3_TH2_dupl",
                           "L2FC_sRNAseq_rhino_TH1_dupl_81510_81511" = "rhino_TH1_dupl")) %>%
  mutate(Genotype = factor(Genotype, levels = genotype.order))

  
# Make data frame with plot annotations of median values and category size:
anno.df <- aggregate(plotL2FC.df$value, by = plotL2FC.df[,c("SL.groups", "Genotype")], FUN=median)
anno.df$length <- aggregate(plotL2FC.df$value, by = plotL2FC.df[,c("SL.groups", "Genotype")], FUN=length)[,"x"]
anno.df <- anno.df %>% 
  rename(medians = x) %>%
  mutate(medians = signif(medians, 3)) %>%
  mutate(percentage = paste(round(2^medians, 3)*100, "%")) %>%
  mutate(fold_change = round(2^abs(medians), 3))

# Make box plot:
my.box <- ggplot(plotL2FC.df, aes(x = SL.groups, y = value, fill = SL.groups)) +
  facet_wrap(~ Genotype, scales = "free_y") +
  geom_hline(yintercept = 0, lwd = 0.5, color="grey") +
  geom_boxplot(position = position_dodge(width=0.9), outlier.shape = NA) +
  scale_y_continuous(breaks=seq(-10, 5, 1)) +
  # xlab("") + ylab("") +
  geom_text(data=anno.df, aes(x=SL.groups, y=0.9*medians, label=length), color="#000000", size=2) +
  geom_text(data=anno.df, aes(x=SL.groups, y=1.1*medians, label=medians), color="#000000", size=2) +
  ggtitle(paste("1kb tile sRNAseq normalized log2(fold-changes)\nLT_sRNAseq_Dec18, RPKM >", MIN_COUNT)) +
  scale_fill_manual(values = c(rep(cbgPalette[1], 2), rep(cbgPalette[6], 5)) , name="") +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text=element_text(family="Helvetica"),
        strip.text.x = element_text(size = 10),
        axis.line = element_line(color = 'black'),
        axis.title = element_text(size=14, color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1, color = 'black', size=10),
        axis.text.y = element_text(color = 'black', size=10),
        legend.position = "") +
  guides(fill = guide_legend(ncol = 2))

# Save produced plots:
filename <- paste("sRNAseq_TH_Cluster_horiBoxPlot_free_y.pdf", sep="") # compile a output file name for the plot
ggsave(filename=paste(output.dir, filename, sep=""),
       plot=my.box, width=8, height=6, units="in", useDingbats=FALSE)
  

################################################################################
# 4. Make whole-genome scatter plots with cluster highlights:

# Define category sets for reshaping data for plotting:
plot.cats <- c("Category", "L2FC_sRNAseq_nxf3_TH2_dupl_81508_81509", 
               "L2FC_sRNAseq_rhino_TH1_dupl_81510_81511",
               "L2FC_sRNAseq_CG13741_TH1_dupl_81506_81507")
other.cats <- c("ClusterFlamenco", "filtered", "Gene body", "Gene TSS", "Other", "Other Rhi-bound")
cat.levels <- c("Other", "ClusterFlamenco", "Cluster20A", "Cluster80F", "Cluster42AB", "Cluster38C1/2", "Other RD-SL")

# Reshape data for plotting:
plot.df <- SummedMaster %>%
  mutate(Category = replace(Category, -grep("Cluster", SummedMaster$Category), "Other")) %>%
  mutate(Category = recode(Category,  "Cluster38C1" = "Cluster38C1/2", "Cluster38C2" = "Cluster38C1/2")) %>%
  mutate(Category = factor(Category, level = cat.levels)) %>%
  filter(SL_CATS == "RD-SL" | Category == "Cluster20A" | Category == "ClusterFlamenco") %>%
  select(one_of(plot.cats))
  

UPCOUNT <- dim(plot.df[plot.df$L2FC_sRNAseq_nxf3_TH2_dupl_81508_81509>0,])[1]
DOWNCOUNT <- dim(plot.df[plot.df$L2FC_sRNAseq_nxf3_TH2_dupl_81508_81509<0,])[1]

scat.plot <- ggplot(plot.df, aes(x = L2FC_sRNAseq_nxf3_TH2_dupl_81508_81509,
                                 y = L2FC_sRNAseq_rhino_TH1_dupl_81510_81511)) +
  geom_point(data=subset(plot.df, Category=="Other"), color="#999999", alpha=0.2, size=2) +
  geom_point(data=subset(plot.df, Category!="Other"), aes(color=Category), alpha=1, size=2) +
  ggtitle("sRNA-seq tile analysis of RD-SL\n(log2 fold-change relative to w1118)") +
  scale_color_manual(values = c(cbgPalette[8], cbgPalette[6], cbgPalette[3], cbgPalette[2], cbgPalette[4], cbgPalette[7])) +
  scale_y_continuous(limits = c(-10.2, 2.5), breaks=seq(-10, 10, 2.5)) +
  scale_x_continuous(limits = c(-10.2, 5),breaks=seq(-10, 10, 2.5)) +
  # coord_cartesian(xlim = c(-10, 5), ylim = c(-10, -2), expand=F) +
  # theme(text=element_text(family="Helvetica", size = 14, color="black")) +
  # theme(strip.text = element_text(size = 14)) + 
  # theme(axis.text.x = element_text(size=14, color="black")) + 
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = 'black'),
        text=element_text(family="Helvetica"),
        axis.title = element_text(size=14, color = 'black'),
        axis.text = element_text(color = 'black', size=14),
        strip.text = element_text(size = 14),
        legend.position="right", legend.text = element_text(size=14), 
        legend.title = element_text(size=14)) + 
  guides(fill = guide_legend(ncol = 1)) + coord_fixed() +
  geom_vline(xintercept = 0, color="#999999", linetype = "dashed") +
  annotate(geom="text", y=0, x=-8, label=DOWNCOUNT) + annotate(geom="text", y=0, x=4, label=UPCOUNT)

# Save produced plots:
filename <- paste("sRNAseq_L2FC_scatterplots_w20A_flamAll_Nxf3KO_vs_RhiKO_Min", MIN_COUNT, ".pdf", sep="") # compile a output file name for the plot
ggsave(filename=paste(output.dir, filename, sep=""), plot=scat.plot, width=7, height=5, units="in", useDingbats=FALSE)

  
# Plot for CG13741 KO:

UPCOUNT <- dim(plot.df[plot.df$L2FC_sRNAseq_CG13741_TH1_dupl_81506_81507>0,])[1]
DOWNCOUNT <- dim(plot.df[plot.df$L2FC_sRNAseq_CG13741_TH1_dupl_81506_81507<0,])[1]

scat.plot <- ggplot(plot.df, aes(x = L2FC_sRNAseq_CG13741_TH1_dupl_81506_81507,
                                 y = L2FC_sRNAseq_rhino_TH1_dupl_81510_81511)) +
  geom_point(data=subset(plot.df, Category=="Other"), color="#999999", alpha=0.2, size=2) +
  geom_point(data=subset(plot.df, Category!="Other"), aes(color=Category), alpha=1, size=2) +
  ggtitle("sRNA-seq tile analysis of RD-SL\n(log2 fold-change relative to w1118)") +
  scale_color_manual(values = c(cbgPalette[8], cbgPalette[6], cbgPalette[3], cbgPalette[2], cbgPalette[4], cbgPalette[7])) +
  scale_y_continuous(limits = c(-10.2, 2.5), breaks=seq(-10, 10, 2.5)) +
  scale_x_continuous(limits = c(-10.2, 5),breaks=seq(-10, 10, 2.5)) +
  # coord_cartesian(xlim = c(-10, 5), ylim = c(-10, -2), expand=F) +
  # theme(text=element_text(family="Helvetica", size = 14, color="black")) +
  # theme(strip.text = element_text(size = 14)) + 
  # theme(axis.text.x = element_text(size=14, color="black")) + 
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = 'black'),
        text=element_text(family="Helvetica"),
        axis.title = element_text(size=14, color = 'black'),
        axis.text = element_text(color = 'black', size=14),
        strip.text = element_text(size = 14),
        legend.position="right", legend.text = element_text(size=14), 
        legend.title = element_text(size=14)) + 
  guides(fill = guide_legend(ncol = 1)) + coord_fixed() +
  geom_vline(xintercept = 0, color="#999999", linetype = "dashed") +
  annotate(geom="text", y=0, x=-8, label=DOWNCOUNT) + annotate(geom="text", y=0, x=4, label=UPCOUNT)

# Save produced plots:
filename <- paste("sRNAseq_L2FC_scatterplots_w20A_flam_CG13KO_vs_RhiKO_Min", MIN_COUNT, ".pdf", sep="") # compile a output file name for the plot
ggsave(filename=paste(output.dir, filename, sep=""), plot=scat.plot, width=7, height=5, units="in", useDingbats=FALSE)


# END OF SCRIPT
################################################################################
################################################################################

