#!/usr/bin/env bash
#$ -S /bin/bash

#$ -clear
#$ -e $JOB_NAME.e$JOB_ID.txt
#$ -o $JOB_NAME.o$JOB_ID.txt
#$ -cwd
#$ -q public.q
#$ -pe smp 10

#SBATCH --cpus-per-task=3
#SBATCH --mem=30g
#SBATCH -e "%x.e.%j.txt"
#SBATCH -o "%x.o.%j.txt"
#XX-SBATCH --qos=short-XX
#SBATCH --time=8:00:00
@SBATCH@

hostname

set -ux

###################################################################################################
#extract variables
VARI=$(echo "$1" | sed 's/,/\t/g;s/"//g')
eval "$VARI"
echo $1 | sed 's/,/\n/g'

TIME=$(date "+%s")
TIMEx=$TIME

#define proper CORES setting
if [[ $GRIDsystem == SLURM ]]; then CORES=$(( SLURM_CPUS_PER_TASK * @HYPER@ )); elif [[ -z ${NSLOTS+x} ]]; then CORES=4; else CORES=$NSLOTS; fi
echo $CORES

#presetting of environment to preload modules and to prepare all files required
#load all functions
source ${SCRIPT_DIR}tools

###################################################################################################
#create local tmp
locTMP="${TMPdir}TMP/"
mkdir -p $locTMP
chmod 777 $locTMP

TMPDIR=$locTMP

cd $locTMP

echo $FOLDER
mkdir ${FOLDER}plots

#prepare input files
cat $FILE_CONTAINING_LIBRARIES | tr ' ' '\t' | cut -f 2 >${locTMP}libNAMES.txt

###################################################################################################

#collect filtered read counts
rm -rf ${FOLDER}annotation_counts.txt
NAMElist=$(cat $FILE_CONTAINING_LIBRARIES | awk '{ print $2 }' | tr '\n' '\t')
NAMElistANN=$(printf "annotation\t$NAMElist")
NAMElistSIZE=$(printf "length\t$NAMElist")

echo $NAMElistANN | tr ' ' '\t' >${FOLDER}annotation_counts.txt

for CATEGORY in adaptor_dimer artifact_filtered length_filtered_short length_filtered_long N_filtered reads_after_filtering reads_not_mapped; do
  SWITCH="N"
  while read LINE; do
    NAME=$(echo $LINE | tr ' ' '\t' | cut -f 2)
    libFOLDER=${FOLDER}/individual-libraries/${NAME}/

    if [[ $SWITCH == N ]]; then
      X=$(cat ${libFOLDER}log.txt | grep $CATEGORY | awk '{ print $NF}')
      Xlist="$CATEGORY $X"
      SWITCH="Y"
    else
      X=$(cat ${libFOLDER}log.txt | grep $CATEGORY | awk '{ print $NF}')
      Xlist="$Xlist $X"
    fi
  done <$FILE_CONTAINING_LIBRARIES

  echo $Xlist | tr ' ' '\t' >>${FOLDER}annotation_counts.txt
done

#---------------------------------------------------------------------------------------------------------
#collect annotation count results and plot values

for CATEGORY in annotation_counts splitup_mRNA splitup_TE; do
  NAMES=""

  if [[ $CATEGORY != annotation_counts ]]; then
    rm -rf ${FOLDER}${EXT}${CATEGORY}.txt
  fi

  awk -v OFS="\t" '{
    X[$1]=X[$1]"\t"$2
    Y[$1]+=NR
  }
  END {
    for (i in X) {
      print Y[i],i,X[i]
    }
  }' <(while read NAME; do
    libFOLDER=${FOLDER}/individual-libraries/${NAME}/
    NAMES="${NAMES} ${NAME}"
    cat ${libFOLDER}${NAME}_${CATEGORY}.txt
  done <${locTMP}libNAMES.txt) |
    sort -k1,1n | cut -f 2- >${locTMP}TMP.txt

  if [[ $CATEGORY != annotation_counts ]]; then
    EXT=annotation_
  else
    EXT=""
  fi

  cat ${locTMP}TMP.txt | tr ' ' '\t' | tr -s "\t" >>${FOLDER}${EXT}${CATEGORY}.txt

  if [[ $TYPE == CHIPseq || $TYPE == DNAseq || $noSTRANDED == Y ]]; then
    mv ${FOLDER}${EXT}${CATEGORY}.txt ${locTMP}${CATEGORY}.txt
    awk -v OFS="\t" '
    BEGIN{
      N=1
    }
    {
      if(NR==1){
        print 1,$0
        nLIBS=NF
      }else{
        for(i=2; i<=NF; i++){
          sub("_AS","", $1)
          X[$1][i]+=$i

          if($1 in Z){
            n=1
          }else{
            Z[$1]=1
            N+=1
            Y[N]=$1
          }
        }
      }
    }
    END{
      for(j=2; j<=N; j++){
        print ">"j,Y[j]
        for(i=2; i<=nLIBS; i++){
          print X[Y[j]][i]
        }
      }
    }' ${locTMP}${CATEGORY}.txt | tr '\n' '\t' | tr '>' '\n' | cut -f 2- >${FOLDER}${EXT}${CATEGORY}.txt
  fi
done

Rscript ${SCRIPT_DIR}plot_annotations.R TMP=$locTMP OUTdir=${FOLDER}plots/ FILE=${FOLDER}annotation_counts.txt PLOTdir=${FOLDER}plots/

#---------------------------------------------------------------------------------------------------------
#collect SLAMseq SNP table

if [[ $SLAM == Y ]]; then
for CATEGORY in all miRNA piRNA; do
  echo "SNP" > ${locTMP}SNPsummary.normalized.${CATEGORY}.txt

  awk -v OFS="\t" '{
    X[$1]=X[$1]"\t"$2
    Y[$1]+=NR
  }
  END {
    for (i in X) {
      print Y[i],i,X[i]
    }
  }' <(while read NAME; do
    libFOLDER=${FOLDER}/individual-libraries/${NAME}/
    echo $NAME >> ${locTMP}SNPsummary.normalized.${CATEGORY}.txt
    cat ${libFOLDER}SNPsummary.normalized.${CATEGORY}.txt
  done <${locTMP}libNAMES.txt) |
    sort -k1,1n | cut -f 2- >${locTMP}TMP.txt

  cat ${locTMP}SNPsummary.normalized.${CATEGORY}.txt | tr '\n' '\t'> ${FOLDER}SNPsummary.normalized.${CATEGORY}.txt
  cat ${locTMP}TMP.txt | tr ' ' '\t' | tr -s "\t" >>${FOLDER}SNPsummary.normalized.${CATEGORY}.txt
done
fi
#---------------------------------------------------------------------------------------------------------
#collect 1kb tile table

if [[ $GENOME_VERSION == dm6 ]]; then
  mkdir ${FOLDER}1kb_tiles
  for STRAND in Plus Minus; do
    awk -v OFS="\t" '
    {
      if($2==-1){$2=0}

      X[$1]=X[$1]"\t"$2
      Y[$1]+=NR
    }
    END {
      for (i in X) {
        print Y[i],i,X[i]
      }
    }' <(while read NAME; do
      libFOLDER=${FOLDER}/individual-libraries/${NAME}/
      NAMES="${NAMES} ${NAME}"
      cat ${libFOLDER}WindowCounts_${STRAND}.txt
    done <${locTMP}libNAMES.txt) | sort -k1,1n | cut -f 2- >${locTMP}TMP.txt

    echo POS $NAMElist | tr ' ' '\t' | tr -s "\t" >${FOLDER}1kb_tiles/WindowCounts_${STRAND}.txt
    cat ${locTMP}TMP.txt | tr ' ' '\t' | tr -s "\t" >>${FOLDER}1kb_tiles/WindowCounts_${STRAND}.txt
  done
  
  cp ${UTILITY_DIR}/dmel/${GENOME_VERSION}/1kb_tile/kbWindows_mainchr_genes_tss_dm6.txt ${FOLDER}1kb_tiles/

fi

#---------------------------------------------------------------------------------------------------------
#collect size profiles and plot results
rawFOLDER=${FOLDER}plots/raw-data/
mkdir -p ${rawFOLDER}

if [[ $TYPE == sRNAseq || $TYPE == CLIPseq || $TYPE == sRNAseqIP || $TYPE == SLAMseq || $TYPE == GROseq ]]; then
  rm -rf ${FOLDER}plots/size-profiles/*
  mkdir ${FOLDER}plots/size-profiles/


  for CATEGORY in all filtered TE rRNA miRNA; do
    SWITCH="N"
    while read NAME; do
      libFOLDER=${FOLDER}/individual-libraries/${NAME}/size-profiles/

      #normalize size profile-tables
      NORM=$(cat ${FOLDER}/individual-libraries/$NAME/normalization.txt | tail -n 1 | tr ' ' '\t' | cut -f 1)

      if [[ $SWITCH == N ]]; then
        awk -v OFS="\t" -v NORM=$NORM '{
          print $1,$2/NORM
        }' ${libFOLDER}size-profile_${NAME}-${CATEGORY}.txt >${locTMP}TMP.txt

        SWITCH="Y"
      else
        awk -v OFS="\t" -v NORM=$NORM '{
          print $1,$2/NORM
        }' ${libFOLDER}size-profile_${NAME}-${CATEGORY}.txt >${locTMP}TMPnorm.txt

        join ${locTMP}TMP.txt ${locTMP}TMPnorm.txt >${locTMP}TMP2.txt
        mv ${locTMP}TMP2.txt ${locTMP}TMP.txt

      fi

    done <${locTMP}libNAMES.txt

    echo $NAMElistSIZE | tr ' ' '\t' >${rawFOLDER}size_profiles-${CATEGORY}.txt
    sort -k1,1n ${locTMP}TMP.txt | tr ' ' '\t' >>${rawFOLDER}size_profiles-${CATEGORY}.txt

    Rscript ${SCRIPT_DIR}plot_compound-size-profile.R TMP=$locTMP CATEGORY=$CATEGORY FOLDER=${rawFOLDER} PLOTdir=${FOLDER}plots/size-profiles/
  done
fi

#---------------------------------------------------------------------------------------------------------
#collect TPM values and plot results

if [[ $TYPE == sRNAseq || $TYPE == sRNAseqIP || $TYPE == CLIPseq ]]; then
  CATEGORIES="TE_RPKM_sense-salmon TE_RPKM_antisense-salmon TE_RPKM_sense-bowtie TE_RPKM_antisense-bowtie TE_counts_sense-bowtie TE_counts_antisense-bowtie"
  if [[ $SLAM == Y ]]; then
    tcCATEGORIES=$(echo $CATEGORIES | sed 's/bowtie/bowtie\.TC/g')
    CATEGORIES="$CATEGORIES $tcCATEGORIES"
  fi
elif [[ $TYPE = RNAseq || $TYPE == RIPseq ]]; then
  CATEGORIES="TE_GeTMM  GeTMM_gene"
elif [[ $FORCEquant_unstranded == Y ]]; then
  CATEGORIES="TE_GeTMM_unistrand GeTMM_gene"
fi

if [[ $TYPE == sRNAseq || $TYPE == sRNAseqIP || $TYPE == CLIPseq || $TYPE = RNAseq || $TYPE == RIPseq || $FORCEquant_unstranded == Y ]]; then
  for CATEGORY in $CATEGORIES; do
    if [[ $CATEGORY == "TE_GeTMM" ]]; then
      head -n 1 ${FOLDER}gene-expression/GeTMM_normalized.txt | tr ' ' '\t' > ${FOLDER}${CATEGORY}_sense.txt
      grep '^TE:' ${FOLDER}gene-expression/GeTMM_normalized.txt | grep -v _AS | tr ' ' '\t' >> ${FOLDER}${CATEGORY}_sense.txt
      head -n 1 ${FOLDER}gene-expression/GeTMM_normalized.txt | tr ' ' '\t' > ${FOLDER}${CATEGORY}_antisense.txt
      grep '^TE:' ${FOLDER}gene-expression/GeTMM_normalized.txt | grep _AS | tr ' ' '\t' >> ${FOLDER}${CATEGORY}_antisense.txt

    elif [[ $CATEGORY == GeTMM_gene ]]; then
      head -n 1 ${FOLDER}gene-expression/GeTMM_normalized.txt | tr ' ' '\t' > ${FOLDER}${CATEGORY}.txt
      grep -v '^TE:' ${FOLDER}gene-expression/GeTMM_normalized.txt | tr ' ' '\t'  >> ${FOLDER}${CATEGORY}.txt
    elif [[ $CATEGORY == TE_GeTMM_unistrand ]]; then
      head -n 1 ${FOLDER}gene-expression/GeTMM_normalized.txt | tr ' ' '\t' > ${FOLDER}${CATEGORY}.txt
      grep '^TE:' ${FOLDER}gene-expression/GeTMM_normalized.txt | tr ' ' '\t' >> ${FOLDER}${CATEGORY}.txt
    else
      NAMES="TE"

      awk -v OFS="\t" '{
        X[$1]=X[$1]"\t"$2
      }
      END {
        for (i in X) {
        print i,X[i]
        }
      }' <(
        while read NAME; do
          libFOLDER=${FOLDER}/individual-libraries/${NAME}/
          NAMES="${NAMES} $NAME"
          cat ${libFOLDER}${CATEGORY}.txt
        done <${locTMP}libNAMES.txt
        echo $NAMES | tr ' ' '\t' >${FOLDER}${CATEGORY}.txt
      ) |
        awk -v OFS="\t"  '{
          X=$1
          gsub("_AS","",X)
          print X,$0
        }' |
        sort -k1,1 | cut --complement -f 1 >${locTMP}TMP.txt

      cat ${locTMP}TMP.txt | sed 's/::/\t/g' | tr ' ' '\t' | tr -s "\t" >>${FOLDER}${CATEGORY}.txt
    fi
  done
fi

#---------------------------------------------------------------------------------------------------------
#collect duplication-rates and sequencing-depth values and plot results
#!not functional at the moment

#!if [[ $TYPE == sRNAseq ]]; then
#!  CATEGORIES="read-duplication sequencing-depth"
#!else
#!  rawFOLDER=${FOLDER}plots/raw-data/
#!  mkdir -p ${rawFOLDER}
#!
#!  CATEGORIES="sequencing-depth"
#!fi
#!
#!NAMES="annotation"
#!
#!for CATEGORY in $CATEGORIES; do
#!  echo $CATEGORY
#!  awk -v OFS="\t" '{
#!    X[$1]=X[$1]"\t"$2
#!  }
#!  END {
#!    for (i in X) {
#!    print i,X[i]
#!    }
#!  }' <(
#!    while read NAME; do
#!      libFOLDER=${FOLDER}/individual-libraries/${NAME}/
#!      NAMES="${NAMES} $NAME"
#!      cat ${libFOLDER}${NAME}_${CATEGORY}.txt
#!    done <${locTMP}libNAMES.txt
#!    echo $NAMES | tr ' ' '\t' >${rawFOLDER}${CATEGORY}.txt
#!  ) |
#!    sort -k1,1n >${locTMP}TMP.txt
#!
#!  cat ${locTMP}TMP.txt | sed 's/::/\t/g' | tr ' ' '\t' | tr -s "\t" >>${rawFOLDER}${CATEGORY}.txt
#!
#!  Rscript ${SCRIPT_DIR}plot_depth_and_duplication.R TMP=$locTMP FILE=${rawFOLDER}${CATEGORY}.txt CATEGORY=$CATEGORY
#!  mv ${locTMP}*${CATEGORY}.pdf ${FOLDER}plots/
#!done
#!

#---------------------------------------------------------------------------------------------------------
#collect spike-in data 

if [[ $spikeINnorm == Y ]]; then
    NAMES=""

    awk -v OFS="\t"  '{
        
        NAME=$1
        if(NAME in X) {
          X[NAME]=X[NAME]"\t"$2
        }else{
          X[NAME]=$2
        }
      }
      END {
        for (i in X) {
          print i,X[i]
        }
      }' <(
      while read NAME; do
        libFOLDER=${FOLDER}/individual-libraries/${NAME}/
        NAMES="${NAMES} ${NAME}"
        cat ${libFOLDER}spike-in.raw.txt
      done <${locTMP}libNAMES.txt
      echo spike-in $NAMES | tr '~' '\t' | tr ' ' '\t' >${FOLDER}plots/raw-data/spike-in.raw.txt
    ) |
      tr '~' '\t' | sort -k2,2rn >>${FOLDER}plots/raw-data/spike-in.raw.txt

    awk -v OFS="\t" -v NORM=$NORM '{
        
        NAME=$1
        if(NAME in X) {
          X[NAME]=X[NAME]"\t"$2/NORM
        }else{
          X[NAME]=$2/NORM
        }
      }
      END {
        for (i in X) {
          print i,X[i]
        }
      }' <(
      while read NAME; do
        libFOLDER=${FOLDER}/individual-libraries/${NAME}/
        NAMES="${NAMES} ${NAME}"
        NORM=$(tail -n 1 ${libFOLDER}normalization.txt | tr ' ' '\t' | cut -f 1)
        cat ${libFOLDER}spike-in.raw.txt
      done <${locTMP}libNAMES.txt
      echo spike-in $NAMES | tr '~' '\t' | tr ' ' '\t' >${FOLDER}plots/raw-data/spike-in.norm.txt
    ) |
    tr '~' '\t' | sort -k2,2rn >>${FOLDER}plots/raw-data/spike-in.norm.txt

fi


#---------------------------------------------------------------------------------------------------------
#collect ping-pong and phasing data

if [[ $PingPong == Y ]] && [[ $TYPE == sRNAseq ]]; then
  #ping-pong 

  #z-score
  # Define the base directory
  BASE_DIR="${FOLDER}individual-libraries/"

  # Get all unique xy directories (extract only folder names, no full path)
  xy_dirs=($(find "$BASE_DIR" -mindepth 2 -maxdepth 2 -type d -name "ping-pong-analysis" | xargs -I{} dirname {} | xargs -n1 basename))

  for STRAND in plus.noFilter minus.noFilter; do
    if [[ $STRAND == "plus"* ]]; then  LINE=12; else LINE=13; fi

    # Get all unique ab.txt file names (only filenames, no full path)
    ab_files=($(find "$BASE_DIR"/*/ping-pong-analysis/  -mindepth 1 -maxdepth 1 -type f -name "*.${STRAND}.txt" -exec basename {} \; | sort -u))

    # Output file
    OUTPUT_FILE="${FOLDER}plots/raw-data/aggregated_ping-pong-data.z-score.${STRAND}.txt"

    # Create a temporary file
    TEMP_FILE=${TMPdir}tmp.txt

    # Print the transposed header (first column = file names)
    echo -e "File\t${xy_dirs[*]}" | tr ' ' '\t' > "$TEMP_FILE"

    # Loop through each ab file (rows)
    for ab in "${ab_files[@]}"; do
        row="$ab"  # Start row with file namee
        # Loop through each xy directory (columns)
        for xy in "${xy_dirs[@]}"; do
            # Construct full path to the expected file
            file_path="$BASE_DIR/$xy/ping-pong-analysis/$ab"

            # Debug: Check if file exists
            if [[ -f "$file_path" ]]; then
                value=$(sed -n '2p' "$file_path")
            else
                echo "WARNING: File not found -> $file_path" >&2
                value="NA"
            fi

            row="$row\t$value"
        done

        # Append row to temp file
        echo -e "$row" >> "$TEMP_FILE"
    done

    # Move the temporary file to the final output file
    mv "$TEMP_FILE" "$OUTPUT_FILE"

    #ping-pong score
    # Output file
    OUTPUT_FILE="${FOLDER}plots/raw-data/aggregated_ping-pong-data.ping-pong-score.${STRAND}.txt"

    # Create a temporary file
    TEMP_FILE=${TMPdir}tmp.txt

    # Print the transposed header (first column = file names)
    echo -e "File\t${xy_dirs[*]}" | tr ' ' '\t' > "$TEMP_FILE"

    # Loop through each ab file (rows)
    for ab in "${ab_files[@]}"; do
        row="$ab"  # Start row with file namee
        # Loop through each xy directory (columns)
        for xy in "${xy_dirs[@]}"; do
            # Construct full path to the expected file
            file_path="$BASE_DIR/$xy/ping-pong-analysis/$ab"

            # Debug: Check if file exists
            if [[ -f "$file_path" ]]; then
                value=$(sed -n "${LINE}p" "$file_path")
            else
                echo "WARNING: File not found -> $file_path" >&2
                value="NA"
            fi

            row="$row\t$value"
        done

        # Append row to temp file
        echo -e "$row" >> "$TEMP_FILE"
    done

    # Move the temporary file to the final output file
    mv "$TEMP_FILE" "$OUTPUT_FILE"
  done
  #---------------------------------------------------------------------------------------------------------
  #phasing
  # Define the base directory
  BASE_DIR="${FOLDER}individual-libraries/"

  # Get all unique xy directories (extract only folder names, no full path)
  xy_dirs=($(find "$BASE_DIR" -mindepth 2 -maxdepth 2 -type d -name "ping-pong-analysis" | xargs -I{} dirname {} | xargs -n1 basename))

  # Get all unique ab.txt file names (only filenames, no full path)
  ab_files=($(find "$BASE_DIR"/*/ping-pong-analysis/  -mindepth 1 -maxdepth 1 -type f -name "*.txt" -exec basename {} \; | sort -u))

  # Output file
  OUTPUT_FILE="${FOLDER}plots/raw-data/aggregated_ping-pong-data.txt"

  # Get all unique xy directories (extract only folder names, no full path)
  xy_dirs=($(find "$BASE_DIR" -mindepth 2 -maxdepth 2 -type d -name "ping-pong-analysis" | xargs -I{} dirname {} | xargs -n1 basename))
 
  # Get all unique ab.txt file names (only filenames, no full path)
  ab_files=($(find "$BASE_DIR"/*/phasing-analysis/  -mindepth 1 -maxdepth 1 -type f -name "*zscores_both_strands.txt" -exec basename {} \; | sort -u))


  for currSTRAND in sense antisense; do
    # Create a temporary file
    TEMP_FILE=${TMPdir}tmp.txt

    # Output file
    OUTPUT_FILE="${FOLDER}plots/raw-data/aggregated_phasing-data_${currSTRAND}.txt"

   if [[ $currSTRAND == "sense" ]]; then INDEX=2; else INDEX=3; fi

    # Print the transposed header (first column = file names)
    echo -e "File\t${xy_dirs[*]}" | tr ' ' '\t' > "$TEMP_FILE"

    # Loop through each ab file (rows)
    for ab in "${ab_files[@]}"; do
        row=$(echo "$ab" | sed 's/_zscores_both_strands.txt//')  # Start row with file name

        # Loop through each xy directory (columns)
        for xy in "${xy_dirs[@]}"; do
            # Construct full path to the expected file
            file_path="$BASE_DIR/$xy/phasing-analysis/$ab"

            # Debug: Check if file exists
            if [[ -f "$file_path" ]]; then
                value=$( sed -n 2p "$file_path" | tr ' ' '\t' | tr -s '\t' | cut -f $INDEX)
            else
                echo "WARNING: File not found -> $file_path" >&2
                value="NA"
            fi

            row="$row\t$value"
        done

        # Append row to temp file
        echo -e "$row" >> "$TEMP_FILE"
    done
  done

  # Move the temporary file to the final output file
  mv "$TEMP_FILE" "$OUTPUT_FILE"
fi

#---------------------------------------------------------------------------------------------------------
#generate TE histograms
if [[ $TYPE == RNAseq ]] || [[ $TYPE == sRNAseq ]] || [[ $TYPE == sRNAseqIP ]] || [[ $TYPE == CHIPseq ]] || [[ $TYPE == CLIPseq ]] || [[ $TYPE == GROseq ]]; then
  if [[ $TYPE == CHIPseq || $TYPE == DNAseq || $noSTRANDED == Y ]]; then
    CATEGORIES="all"
  else
    CATEGORIES="sense antisense"
  fi

  if [[ $TYPE == sRNAseq || $TYPE == sRNAseqIP ]]; then
    CLASSES="man man.21mer"
  else
    CLASSES="man"
  fi

  if [[ ${TYPE} == sRNAseq || ${TYPE} == sRNAseqIP ]]; then
    STRINGENCYvec="normal all"
  else
    STRINGENCYvec="normal"
  fi

  for STRINGENCY in $STRINGENCYvec; do
    if [[ $STRINGENCY == all ]]; then
      STRINGENCYext=".all"
      FOLDERext="all-alignments"
    else
      STRINGENCYext=""
      FOLDERext="stringent-no-interTE-alignments"
    fi

    for CATEGORY in $CATEGORIES; do
      for currCLASS in $CLASSES; do
        NAMES=""

        awk -v OFS="\t" '{
            
            NAME=$1"~"$2
            if(NAME in X) {
              X[NAME]=X[NAME]"\t"$3
            }else{
              X[NAME]=$3
            }
          }
          END {
            for (i in X) {
              print i,X[i]
            }
          }' <(
          while read NAME; do
            libFOLDER=${FOLDER}/individual-libraries/${NAME}/
            NAMES="${NAMES} ${NAME}"
            cat ${libFOLDER}${NAME}_TE-${CATEGORY}_${currCLASS}${STRINGENCYext}.bg
          done <${locTMP}libNAMES.txt
          echo TE POS $NAMES | tr '~' '\t' | tr ' ' '\t' >${FOLDER}plots/raw-data/TE_hist_${CATEGORY}.${currCLASS}${STRINGENCYext}.txt
        ) |
          sort -k2,2n >${locTMP}TMP_TEhist_${currCLASS}${STRINGENCYext}.txt

        echo $locTMP

        cat ${locTMP}TMP_TEhist_${currCLASS}${STRINGENCYext}.txt | tr '~' '\t' | sort -k2,2n | tr ' ' '\t' >>${FOLDER}plots/raw-data/TE_hist_${CATEGORY}.${currCLASS}${STRINGENCYext}.txt
      done
    done
    mkdir -p ${FOLDER}plots/TEhist_${FOLDERext}/html-dependencies/
    Rscript ${SCRIPT_DIR}plot_TEhist.R OPENdir=${FOLDER}plots/TEhist_${FOLDERext}/ INPUT=${FOLDER}plots/raw-data/TE_hist TYPE=$TYPE STRINGENCYext=${STRINGENCYext}.txt
  done


fi

#---------------------------------------------------------------------------------------------------------
#move plots to open directory

rm -rf ${locTMP}Rscripts.pdf

#---------------------------------------------------------------------------------------------------------
#cleanup
if [[ $DEBUG != Y ]]; then
  rm -rf $locTMP
fi
PROCESSED_TIME=$(echo -e $(date "+%s") $TIME | awk '{ print ($1-$2)/60 }')
echo "collect-numbers - processing_time=" ${PROCESSED_TIME} >>${libFOLDER}time-log.txt

exit
