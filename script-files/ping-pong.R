library(tidyverse)
library(cowplot)
library(ggforce)
library(ggpubr)
theme_set(theme_bw())

pdf.options(useDingbats = FALSE)


##################################################################################################
args  =  commandArgs(TRUE)

argmat  =  sapply(strsplit(args, "="), identity)

for (i in seq.int(length = ncol(argmat))) {
  assign(argmat[1, i], argmat[2, i])
}

# available variables
print(ls())

##################################################################################################

#ping-pont analysis
outDIR=paste(baseDIR,"/ping-pong-analysis/",sep="")
setwd(outDIR)
# INFILE="TE_pingpong.txt"

### get the table file name for Brennecke linkage
RAW=read_tsv(INFILE,col_types=list(TE=col_character()))
RAW

### 5'plus-5'minus pingpong --- START ---
#determine number of rows

START=0
END=19
POSITION=9
START.position_B=1-START
LINK.position=POSITION-START+1
SPAN=END-START+1
THRESHOLD=-1

TEs = RAW %>%
  distinct(TE) %>%
  pull(TE)

for(currTE in TEs){
  print(currTE)

  TABLE = RAW %>%
    filter(TE==currTE)

  LENGTH = nrow(TABLE)
  END.position_B = LENGTH-END

  table_a = TABLE$`5END_plus`
  table_b = TABLE$`5END_minus`
  
  #not sure why this is required
  if (POSITION > 0){
    start.position_R = 1
    end.position_R = LENGTH-POSITION
  } else {
    start.position_R = 1-POSITION
    end.position_R = LENGTH
  }
  
  ### define a data.frame and a vector
  offset <- data.frame()
  sum_offset = NULL
  
  ### core of Brennecke linkage calculation
  offset <- sapply(1:SPAN, function(i) {
    sapply(START.position_B:END.position_B, function(j) {
      # Only calculate offset if the initiating 5' end is above threshold
      if(table_a[j] > THRESHOLD) {
        table_a[j]*table_b[j+START+i-1]/sum(table_b[(j+START):(j+END)])
      } else {
        0
      }
    })
  })

  offset[is.na(offset)] <- 0
  
  total_offset_sum = sum(offset)
  
  # Initialize sum_offset3 with zeros
  sum_offset3 = rep(0, SPAN + 1)
  
  # Only calculate z-score and offsets if we have valid pairs
  if(total_offset_sum > 0) {
    for (i in 1:SPAN) {
      sum_offset[i] = sum(offset[,i])/total_offset_sum
    }

    plotNAME = paste("pingpong-offset - ", currTE, " - plus_5end vs minus_5end",sep="")
    fileNAME = paste(outDIR,currTE, ".plus.noFilter.pdf",sep="")
    
    pdf(file=fileNAME,width=7,height=6)
    plot(sum_offset,
         main=plotNAME,
         xlab="offset [nt]",
         ylim=c(0,0.6))
    dev.off()
    
    sum_offset2 = sum_offset[c(1:(LINK.position-1),(LINK.position+1):SPAN)]
    m = mean(sum_offset2)
    s = sd(sum_offset2)
    z = (sum_offset[LINK.position]-m)/s
    sum_offset3 = c(z,sum_offset)
  }

  # Write table regardless of whether we found valid pairs
  tableNAME = paste(outDIR,currTE,".plus.noFilter.txt",sep="")
  write.table(sum_offset3, file=tableNAME, quote=F, row.names=F)
}

### 5'minus->5'plus pingpong --- START ---
#determine number of rows

START=-19
END=0
POSITION=-9
START.position_B=1-START
LINK.position=POSITION-START+1
SPAN=END-START+1
THRESHOLD=-1

TEs = RAW %>%
  distinct(TE) %>%
  pull(TE)

for(currTE in TEs){
  print(currTE)

  TABLE = RAW %>%
    filter(TE==currTE)

  LENGTH = nrow(TABLE)
  END.position_B = LENGTH-END

  table_a = TABLE$`5END_minus`
  table_b = TABLE$`5END_plus`
  
  #not sure why this is required
  if (POSITION > 0){
    start.position_R = 1
    end.position_R = LENGTH-POSITION
  } else {
    start.position_R = 1-POSITION
    end.position_R = LENGTH
  }
  
  ### define a data.frame and a vector
  offset <- data.frame()
  sum_offset = NULL
  
  ### core of Brennecke linkage calculation
  offset <- sapply(1:SPAN, function(i) {
    sapply(START.position_B:END.position_B, function(j) {
      # Only calculate offset if the initiating 5' end is above threshold
      if(table_a[j] > THRESHOLD  ) {
        table_a[j]*table_b[j+START+i-1]/sum(table_b[(j+START):(j+END)])
      } else {
        0
      }
    })
  })

  offset[is.na(offset)] <- 0
  
  total_offset_sum = sum(offset)
  
  # Initialize sum_offset3 with zeros
  sum_offset3 = rep(0, SPAN + 1)
  
  # Only calculate z-score and offsets if we have valid pairs
  if(total_offset_sum > 0) {
    for (i in 1:SPAN) {
      sum_offset[i] = sum(offset[,i])/total_offset_sum
    }

    plotNAME = paste("pingpong-offset - ", currTE, " - plus_5end vs minus_5end",sep="")
    fileNAME = paste(outDIR,currTE, ".minus.noFilter.pdf",sep="")
    
    pdf(file=fileNAME,width=7,height=6)
    STARTx=START-1
    ENDx=END-1
    plot(STARTx:ENDx, sum_offset,
         main=plotNAME,
         xlab="offset [nt]",
         ylim=c(0,0.6))
    dev.off()
    
    sum_offset2 = sum_offset[c(1:(LINK.position-1),(LINK.position+1):SPAN)]
    m = mean(sum_offset2)
    s = sd(sum_offset2)
    z = (sum_offset[LINK.position]-m)/s
    sum_offset3 = c(z,sum_offset)
  }

  # Write table regardless of whether we found valid pairs
  tableNAME = paste(outDIR,currTE,".minus.noFilter.txt",sep="")
  write.table(sum_offset3, file=tableNAME, quote=F, row.names=F)
}


##################################################################################################
#phasing

outDIR=paste(baseDIR,"/phasing-analysis/",sep="")
setwd(outDIR)

#my-approach sequencing count corrected
#phasing with +/-10 window around 3' end
WINDOW_SIZE = 10  # Window size on each side of the 3' end
TARGET_POSITION = 1  # The +1 position we want to calculate z-score for
window_size_total = 2*WINDOW_SIZE + 1  # Total window size (-10 to +10)
NORMfactor=as.numeric(NORMfactor)

for(currTE in TEs){
  print(currTE)
  
  TABLE = RAW %>%
    filter(TE==currTE)
  
  LENGTH=nrow(TABLE)
  
  # Create a function to analyze both strands
  analyze_strand = function(strand_name, table_3end, table_5end, target_pos) {
    # Create a matrix to store counts for each position in the window
    position_counts = matrix(0, nrow=LENGTH, ncol=window_size_total)

    # For each position with a 3' end
    for(pos in 1:LENGTH) {
      if(table_3end[pos]/NORMfactor > 10) {  # If there's a 3' end at this position
        # Look at window around this position
        NORMcount=0
        for(offset in -WINDOW_SIZE:WINDOW_SIZE) {
          window_pos = pos + offset
          
          # Check if window position is within bounds
          if(window_pos >= 1 && window_pos <= LENGTH) {
            # Record 5' end count at this offset, as the fraction of the 3' end count
            NORMcount= NORMcount + table_5end[window_pos]
            position_counts[pos, offset+WINDOW_SIZE+1] = table_5end[window_pos] * table_3end[pos]
          }
        }
        if(NORMcount > 0) {
          position_counts[pos, ] = position_counts[pos, ] / NORMcount
        }else{
          position_counts[pos, ] = rep(0, window_size_total)
        }
      }
    }
    
    # Sum counts for each offset position across all 3' ends
      normalized_counts = colSums(position_counts, na.rm=TRUE)
        
    # Calculate z-score for the target position (+1 for plus, -1 for minus)
    target_index = WINDOW_SIZE + 1 + target_pos  # Index for target position
    
    # Get all values except the target position
    other_positions = normalized_counts[-target_index]
    
    # Calculate mean and standard deviation of other positions
    m = mean(other_positions)
    s = sd(other_positions)
    
    # Calculate z-score for target position
    if(s > 0) {
      z_score = (normalized_counts[target_index] - m) / s
    } else {
      z_score = 0  # Handle case where standard deviation is 0
    }
    
    return(list(
      normalized_counts = normalized_counts,
      z_score = z_score,
      target_index = target_index,
      target_pos = target_pos
    ))
  }
  
  # Analyze sense strand (plus) - evaluate +1 position
  plus_results = analyze_strand(
    "plus", 
    TABLE$`3END_plus`,  # 3' ends
    TABLE$`5END_plus`,  # 5' ends
    TARGET_POSITION     # +1 for plus strand
  )
  
  # Analyze antisense strand (minus) - evaluate -1 position
  minus_results = analyze_strand(
    "minus", 
    TABLE$`3END_minus`,  # 3' ends
    TABLE$`5END_minus`,  # 5' ends
    -1                   # -1 for minus strand
  )
  
  # Prepare data for plotting
  offsets = -WINDOW_SIZE:WINDOW_SIZE
  
  # Create plot with 2 facets
  plotNAME = paste("3prime-5prime linkage - ", currTE, " - window around 3' end", sep="")
  fileNAME = paste(outDIR, currTE, "_3prime_window_both_strands.pdf", sep="")
  
  pdf(file=fileNAME, width=10, height=6)
  par(mfrow=c(1,2), mar=c(5,4,4,2))
  
  # Plot for sense strand (plus)
  plot(offsets, plus_results$normalized_counts,
       main=paste(currTE, "- Sense"),
       xlab="Offset from 3' end [nt]",
       ylab="Normalized 5' end frequency",
       ylim=c(0, max(c(plus_results$normalized_counts, minus_results$normalized_counts))*1.1))
  
  # Highlight the +1 position
  points(plus_results$target_pos, plus_results$normalized_counts[plus_results$target_index], 
         col="red", pch=19, cex=1.5)
  
  # Add z-score to the plot
  text(plus_results$target_pos, plus_results$normalized_counts[plus_results$target_index]*1.05, 
       paste("z =", round(plus_results$z_score, 2)), pos=3, col="red")
  
  # Add a vertical line at position 0 (the 3' end position)
  abline(v=0, lty=2, col="gray")
  
  # Plot for antisense strand (minus)
  plot(offsets, minus_results$normalized_counts,
       main=paste(currTE, "- Antisense"),
       xlab="Offset from 3' end [nt]",
       ylab="Normalized 5' end frequency",
       ylim=c(0, max(c(plus_results$normalized_counts, minus_results$normalized_counts))*1.1))
  
  # Highlight the -1 position
  points(minus_results$target_pos, minus_results$normalized_counts[minus_results$target_index], 
         col="blue", pch=19, cex=1.5)
  
  # Add z-score to the plot
  text(minus_results$target_pos, minus_results$normalized_counts[minus_results$target_index]*1.05, 
       paste("z =", round(minus_results$z_score, 2)), pos=3, col="blue")
  
  # Add a vertical line at position 0 (the 3' end position)
  abline(v=0, lty=2, col="gray")
  
  dev.off()
  
  # Save results to a text file for both strands
  result_data = data.frame(
    Offset = offsets,
    NormalizedCount_Sense = plus_results$normalized_counts,
    ZScore_Sense = c(rep(NA, plus_results$target_index-1), plus_results$z_score, 
                     rep(NA, window_size_total-plus_results$target_index)),
    NormalizedCount_Antisense = minus_results$normalized_counts,
    ZScore_Antisense = c(rep(NA, minus_results$target_index-1), minus_results$z_score, 
                         rep(NA, window_size_total-minus_results$target_index))
  )
  
  tableNAME = paste(outDIR, currTE, "_3prime_window_both_strands.txt", sep="")
  write.table(result_data, file=tableNAME, quote=F, row.names=F, sep="\t")
  
  # Also save just the z-scores for easy reference
  z_score_data = data.frame(
    TE = currTE, 
    ZScore_Sense_Plus1 = plus_results$z_score,
    ZScore_Antisense_Minus1 = minus_results$z_score
  )
  z_score_tableNAME = paste(outDIR, currTE, "_zscores_both_strands.txt", sep="")
  write.table(z_score_data, file=z_score_tableNAME, quote=F, row.names=F, sep="\t")
  
}
