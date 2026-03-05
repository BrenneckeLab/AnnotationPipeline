args  =  commandArgs(TRUE);
argmat  =  sapply(strsplit(args, "="), identity)

for (i in seq.int(length=ncol(argmat))) {
  assign(argmat[1, i], argmat[2, i])
}

# available variables
print(ls())

###################################################################################################

library(gridExtra)
library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(plotly)
library(wasabi)
library(Cairo)
library(stringi)

###################################################################################################

setwd( TMP)

###################################################################################################
#read in info
INFO=paste(TMP, "DGE_info_all.txt", sep="")
s2c_full = read.table(INFO, header = TRUE, as.is=TRUE )

T2G=paste(TMP, "transcript_to_gene.txt", sep="")
t2g=read.table(T2G, stringsAsFactors = FALSE)
colnames(t2g)=c("target_id","ens_gene")


#---------------------------------------------------------------------------------------------------------
if (ITERATION == 1 ){
  #prepare results for sleuth using wasabi
  FILES=s2c_full[,3]
  prepare_fish_for_sleuth(FILES)
  
  q()
}
#---------------------------------------------------------------------------------------------------------
#do stuff to get PCA for all genotypes

#analysis using all libraries
s2c=s2c_full

library(sleuth)


so <- sleuth_prep(s2c, ~condition, target_mapping = t2g, aggregation_column = 'ens_gene',
                  extra_bootstrap_summary = TRUE,
                  transformation_function_counts = function(x) log2(x + 0.5), #change to log2 based for b(effect size)
                  read_bootstrap_tpm = TRUE,  ## required if sleuth_fit uses which_var = "obs_tpm"
                  gene_mode = TRUE) ##use old method and not p-value aggregation

FILENAME=paste(FOLDER, "gene_TPM.txt", sep="")
as.tibble(kallisto_table(so)) %>%
  select(-one_of(c("condition","scaled_reads_per_base"))) %>%
  filter(! str_detect(target_id, "^TE:")) %>%
  spread( sample,tpm )%>%
  write_tsv(FILENAME)


FILENAME=paste(FOLDER, "TE_TPM_sense.txt", sep="")
as.tibble(kallisto_table(so)) %>%
  select(-one_of(c("condition","scaled_reads_per_base"))) %>%
  filter( str_detect(target_id, "^TE:"), ! str_detect(target_id, "_AS")) %>%
  spread( sample,tpm )%>%
  mutate(target_id = str_replace(target_id, "TE:", "")) %>%
  arrange(target_id) %>%
  write_tsv(FILENAME)

FILENAME=paste(FOLDER, "TE_TPM_antisense.txt", sep="")
as.tibble(kallisto_table(so)) %>%
  select(-one_of(c("condition","scaled_reads_per_base"))) %>%
  filter( str_detect(target_id, "^TE:") , str_detect(target_id, "_AS")) %>%
  spread( sample,tpm )%>%
mutate(target_id = str_replace(target_id, "TE:", "")) %>%
  arrange(target_id) %>%
  write_tsv(FILENAME)

  
FILENAME=paste(TMP, "sleuth_object.so", sep="")
sleuth_save(so, FILENAME)

FILENAME=paste(FOLDER, "/plots//all-samples_heatmap.pdf", sep="")
p = plot_sample_heatmap(so)
ggsave(FILENAME,p, width=30 , height=20)


rm(so)

