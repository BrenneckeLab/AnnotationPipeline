#!/usr/bin/env bash
#$ -S /bin/bash

#$ -clear
#$ -e $JOB_NAME.e$JOB_ID.txt
#$ -o $JOB_NAME.o$JOB_ID.txt
#$ -cwd
#$ -q public.q
#$ -pe smp 6

#SBATCH --cpus-per-task=3
#SBATCH --mem=10g
#SBATCH -e "%x.e.%j.txt"
#SBATCH -o "%x.o.%j.txt"
#XX-SBATCH --qos=short-XX
#SBATCH --time=1:00:00
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
if [[ $GRIDsystem == SLURM ]]; then CORES=$SLURM_CPUS_PER_TASK; elif [[ -z ${NSLOTS+x} ]]; then CORES=4; else CORES=$NSLOTS; fi
echo $CORES

#presetting of environment to preload modules and to prepare all files required
#load all functions
source ${SCRIPT_DIR}tools

###################################################################################################
#create local tmp
locTMP="${TMPdir}TMP_1kb/"

mkdir -p $locTMP
chmod 777 $locTMP

###################################################################################################
#determine normalization factor for libry
NORMfactor=$(tail -n 1 <${libFOLDER}normalization.txt | tr ' ' '\t' | cut -f 1)

#process both strands
for STRAND in Plus Minus; do
  #create strand filtering symbol
  if [[ $STRAND == Plus ]]; then STRANDsym="+"; else STRANDsym="-"; fi

  #quantify tiles
  gunzip -c ${TMPdir}${NAME}_annotated.bed.gz |
    #pre-filter reads
    mawk -v OFS="\t" -v STRANDsym=$STRANDsym '
    {
      split($4, splitNAME, /:|=/)
      if($4~"mapping=u" && $4!~"ann=miRNA|ann=mito|ann=rRNA|ann=snoRNA|ann=snRNA|ann=tRNA" && length(splitNAME[2]) > 22 && $6 == STRANDsym) {
        print "chr"$0
      }
    }' |
    #intersect reads with the tiles
    bedtools intersect -a ${UTILITY_DIR}dmel/dm6/1kb_tile/dm6_1kb_windows_DNAseq_mappability_mainchr_sorted.txt -b stdin -nobuf -sorted -wao -F 0.51 |
    #calculate sum per tile
    awk -v OFS="\t" -v locTMP=$locTMP '
    BEGIN{
      COUNT=0
      totalCOUNT=0
    }
    {
      if(TILE==$8){
        COUNT+=$13
      }else{
        if(NR>1) print TILE,COUNT
        totalCOUNT+=COUNT
        TILE=$8
        COUNT=$13
      }
    }
    END{
      print TILE,COUNT
      print "totalCOUNT",totalCOUNT >locTMP "totalCount.tmp" 
    }
    ' >${locTMP}tileCounts_$STRAND.tmp

  echo NORMfactor $NORMfactor | tr ' ' '\t' > ${libFOLDER}WindowCounts_${STRAND}.txt
  cat ${locTMP}totalCount.tmp | tr ' ' '\t'  >> ${libFOLDER}WindowCounts_${STRAND}.txt
  cat ${locTMP}tileCounts_$STRAND.tmp | tr ' ' '\t'  >> ${libFOLDER}WindowCounts_${STRAND}.txt

done

#---------------------------------------------------------------------------------------------------------

if [[ $DEBUG != Y ]]; then
  rm -rf $locTMP
fi

PROCESSED_TIME=$(echo -e "$(date "+%s")" "$TIME" | mawk '{ print ($1-$2)/60 }')
printf "1kb tile analysis - processing_time= ${PROCESSED_TIME} \n" >>"${libFOLDER}time-log.txt"

exit
