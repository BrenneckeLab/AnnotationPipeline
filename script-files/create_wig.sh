#!/usr/bin/env bash
#$ -S /bin/bash

#$ -clear
#$ -e $JOB_NAME.e$JOB_ID.txt
#$ -o $JOB_NAME.o$JOB_ID.txt
#$ -cwd
#$ -q public.q
#$ -pe smp 5

#SBATCH --cpus-per-task=5
#SBATCH --mem=25g
#SBATCH -e "%x.e.%j.txt"
#SBATCH -o "%x.o.%j.txt"
#XX-SBATCH --qos=short-XX
#SBATCH --time=8:00:00
@SBATCH@

hostname

set -ux

###################################################################################################
#extract variables
VARI=$(echo $1 | sed 's/,/\t/g;s/"//g')
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
locTMP=${TMPdir}TMP${TRACK_TYPE_NR}/

mkdir -p "$locTMP"
chmod 777 "$locTMP"

export TMPDIR=$locTMP

###################################################################################################
#preperation of variables

#define length filtering variables
if [[ $TRACK_TYPE_NR == 1 ]] || [[ $TRACK_TYPE_NR == 5 ]]; then
  EXTENSION=""
elif [[ $TRACK_TYPE_NR == 2 ]]; then
  EXTENSION="_siRNA"
elif [[ $TRACK_TYPE_NR == 3 ]]; then
  EXTENSION="_piRNA"
elif [[ $TRACK_TYPE_NR == 4 ]]; then
  EXTENSION="_5ends"
elif [[ $TRACK_TYPE_NR == 6 ]]; then
  EXTENSION="_CLIPtags"
fi

#-------------------------------------------------------------------------------------------------------
#prepare hub folders
HUB_FOLDER=${TMPdir}wig-files/
mkdir -p ${HUB_FOLDER}/

#create extension for UCSC
if [[ $GENOME_VERSION != ASM ]]; then
  CHRext="chr"
else
  CHRext=""
fi

#-------------------------------------------------------------------------------------------------------
#extract uniq mappers and filter for mito and rRNA

if [[ $TRACK_TYPE_NR -le 4 ]]; then
  gunzip -c ${TMPdir}${NAME}_annotated.bed.gz |
    mawk -v locTMP=$locTMP -v NAME=$NAME -v OFS="\t" -v TRACK_TYPE_NR=$TRACK_TYPE_NR -v CHRext=$CHRext -v only5end=${only5end} '
    { 
      if ( only5end == "Y" ){
        if( $6 == "+" ) $3=$2+1
        if( $6 == "-" ) $2=$3-1
      }
      if ($1 == "dmel_mitochondrion_genome" || $1 == "M" || $1 == "rRNA_precursor" ) ; else {
        if( TRACK_TYPE_NR == 1){  
          $1=CHRext$1
          if( $4 ~ "mapping=u" ) {
            print $0 > locTMP NAME "_uniq.bed"
          }
          print $0 > locTMP NAME "_all.bed"
        } else {
          if( TRACK_TYPE_NR == 2){
            split ($4,SPLIT,/:/)
            if (length(SPLIT[2]) == 21 ) {
              $1=CHRext$1
              if( $4 ~ "mapping=u" ) {
                print $0 > locTMP NAME "_uniq.bed"
              }
              print $0 > locTMP NAME "_all.bed"
            }
          } else {
            if(TRACK_TYPE_NR == 3){
              split ($4,SPLIT,/:/)
              if(length(SPLIT[2]) >= 23 ) {
                $1=CHRext$1
                if( $4 ~ "mapping=u" ) {
                  print $0 > locTMP NAME "_uniq.bed"
                }
                print $0 > locTMP NAME "_all.bed"
              }
            } else {
              if(TRACK_TYPE_NR == 4){
                $1=CHRext$1
                if($6=="+") {$3=$2+1}
                if($6=="-") {$2=$3-1}
                
                if( $4 ~ "mapping=u" ) {
                  print $0 > locTMP NAME "_uniq.bed"
                }
                print $0 > locTMP NAME "_all.bed"
              }
            }  
          }  
        } 
      }
    }'

  if [[ $SLAM == Y ]]; then
    grep TCslam ${locTMP}${NAME}_uniq.bed > ${locTMP}${NAME}_uniq-TC.bed
    grep TCslam ${locTMP}${NAME}_all.bed > ${locTMP}${NAME}_all-TC.bed
  fi
elif [[ $TRACK_TYPE_NR -eq 5 || $TRACK_TYPE_NR -eq 6 ]]; then
  #!did not yet add 5end mode for RNAsq
  gunzip -c ${TMPdir}${NAME}_annotated.bed.gz |
    #filter CLIPtags for clip track
    awk -v OFS="\t" -v TRACK_TYPE_NR=$TRACK_TYPE_NR '{
      if(TRACK_TYPE_NR==6){
        if($4~"CLIPtag=TC" || $4~"CLIPtag=Tdel"){
          print
        }
      }else{
        print
      }
    }' |
    awk -v locTMP=$locTMP -v NAME=$NAME -v OFS="\t" -v CHRext=$CHRext '
    { 
      if ($1 == "dmel_mitochondrion_genome" || $1 == "M" || $1 == "rRNA_precursor" ) ; else {
        fullLINE=gensub(/\t/, "::", "G")
        origNAME=$4":SP:"fullLINE

        n=split($11, splitENDS, /,/)
        split($12, splitSTARTS, /,/)
        coordSTART=$2
        coordEND=$3

        $1=CHRext$1

        if( n>1 ) {
          if($4~"mapping=u") {
            print $0 > locTMP"spliced_reads_uniq.bed12"
          }
          print $0 > locTMP"spliced_reads_all.bed12"

          #split up blocks for HomerTools
          for( i=1; i<=n; i++) {
            $2=coordSTART+splitSTARTS[i]
            $3=$2+splitENDS[i]
            $4=origNAME":PC:"i":PX:"n
            $5=n
            $7=$2
            $8=$3
            $10=1
            $11=splitENDS[i]
            $12=0
            if( $4 ~ "mapping=u" ) {
              print $0 > locTMP NAME "_uniq.bed"
            }
            print $0 > locTMP NAME "_all.bed"

          }
        } else {
          if( $4 ~ "mapping=u" ) {
            print $0 > locTMP NAME "_uniq.bed"
          }
          print $0 > locTMP NAME "_all.bed"

        }
      }
    }'

  for MULTI in uniq all; do
    cat ${UTILITY_LOCATION}chrom.sizes | tr ' ' '\t' >${locTMP}chrom.sizes
    bedtools bedtobam -bed12 -i ${locTMP}spliced_reads_${MULTI}.bed12 -g ${locTMP}chrom.sizes >${locTMP}spliced_reads_sort_${MULTI}.bam
    samtools sort -O bam -o ${HUB_FOLDER}/${NAME}_spliced_${MULTI}.bam ${locTMP}spliced_reads_sort_${MULTI}.bam
    samtools index ${HUB_FOLDER}/${NAME}_spliced_${MULTI}.bam 

    bedtools bedtobam -bed12 -i ${TMPdir}splice-junctions_${MULTI}.bed12 -g ${locTMP}chrom.sizes >${TMPdir}splice-junctions_${MULTI}.bam
    samtools sort -O bam -o ${HUB_FOLDER}/${NAME}_junctions_${MULTI}.bam ${TMPdir}splice-junctions_${MULTI}.bam
    samtools index ${HUB_FOLDER}/${NAME}_junctions_${MULTI}.bam
  done
  
  if [[ $SLAM == Y ]]; then
    grep TCslam ${locTMP}${NAME}_uniq.bed > ${locTMP}${NAME}_uniq-TC.bed
    grep TCslam ${locTMP}${NAME}_all.bed > ${locTMP}${NAME}_all-TC.bed
  fi

fi
#-------------------------------------------------------------------------------------------------------
if [[ $SLAM == Y ]]; then
  TRACKmulti="uniq all uniq-TC all-TC"
else
  TRACKmulti="uniq all"
fi

#-------------------------------------------------------------------------------------------------------
#create TagDirectories
for ID in $TRACKmulti; do
  makeTagDirectory ${locTMP}tagdir_${ID} ${locTMP}${NAME}_${ID}.bed -format bed -keepAll -single -force5th -fragLength given &
done
wait

#-------------------------------------------------------------------------------------------------------
#create wig tracks
if [[ ($TYPE == RNAseq* || $TYPE == sRNAseq || $TYPE == sRNAseqIP || $TYPE == GROseq || $TYPE == CapSeq || $TYPE == CLIPseq || $TYPE == RIPseq ) && ${noSTRANDED} != Y ]]; then
  for ID in $TRACKmulti; do
    makeUCSCfile ${locTMP}tagdir_${ID} -fsize 1e20 -strand + -fragLength given -noadj -normLength 0 |
      mawk -v norm=$NORMfactor -v OFS="\t" '{ print $1,$2,$3,$4/norm }' >${locTMP}${NAME}_${ID}_sense.bedgraph
    bedGraphToBigWig ${locTMP}${NAME}_${ID}_sense.bedgraph ${UTILITY_LOCATION}chrom.sizes ${HUB_FOLDER}/${NAME}_${ID}${EXTENSION}_sense.bw

    makeUCSCfile ${locTMP}tagdir_${ID} -fsize 1e20 -strand - -fragLength given -noadj -neg -normLength 0 |
      mawk -v norm=$NORMfactor -v OFS="\t" '{ print $1,$2,$3,$4/norm }' >${locTMP}${NAME}_${ID}_antisense.bedgraph
    bedGraphToBigWig ${locTMP}${NAME}_${ID}_antisense.bedgraph ${UTILITY_LOCATION}chrom.sizes ${HUB_FOLDER}/${NAME}_${ID}${EXTENSION}_antisense.bw
  done 
elif [[ $TYPE == CHIPseq ]] || [[ $TYPE == DNAseq ]] || [[ $noSTRANDED == Y ]]; then
  for ID in $TRACKmulti; do
    makeUCSCfile ${locTMP}tagdir_${ID} -fsize 1e20 -fragLength given -noadj -normLength 0 |
      mawk -v norm=$NORMfactor -v OFS="\t" '{ print $1,$2,$3,$4/norm }' >${locTMP}${NAME}_${ID}.bedgraph
    bedGraphToBigWig ${locTMP}${NAME}_${ID}.bedgraph ${UTILITY_LOCATION}chrom.sizes ${HUB_FOLDER}/${NAME}_${ID}${EXTENSION}.bw
    if [[ $RATIOtracks == Y ]]; then 
      mkdir -p ${TMPdir}bedgraph/
      mv ${locTMP}${NAME}_${ID}.bedgraph ${TMPdir}bedgraph/
    fi
  done
fi

#---------------------------------------------------------------------------------------------------------

if [[ $DEBUG != Y ]]; then
  rm -rf $locTMP
fi

PROCESSED_TIME=$(echo -e "$(date "+%s")" "$TIME" | mawk '{ print ($1-$2)/60 }')
printf "create wig (TrackNr=${TRACK_TYPE_NR}) - processing_time= ${PROCESSED_TIME} \n" >>"${libFOLDER}time-log.txt"

exit
