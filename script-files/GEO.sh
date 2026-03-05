#!/usr/bin/env bash
#$ -S /bin/bash

#$ -clear
#$ -e $JOB_NAME.e$JOB_ID.txt
#$ -o $JOB_NAME.o$JOB_ID.txt
#$ -cwd
#$ -q public.q
#$ -pe smp 2

#SBATCH --cpus-per-task=2
#SBATCH --mem=30g
#SBATCH -e "%x.e.%j.txt"
#SBATCH -o "%x.o.%j.txt"
#XX-SBATCH --qos=short-XX
#SBATCH --time=1:00:00
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
locTMP="${TMPdir}TMP_download+demux/"

mkdir -p "$locTMP"
chmod 777 "$locTMP"

export TMPDIR=$locTMP

###################################################################################################
#process individual libraries
rm -rf ${FOLDER}GEO/

splitFILES=$(echo $FILES | tr '~' '\t')

mkdir -p ${FOLDER}GEO/RAW
mkdir -p ${FOLDER}GEO/BW


#execute function in parallel
parallel -j $CORES ". ${SCRIPT_DIR}functions; collectDATA {1} $TMPdir $FOLDER $TYPE $locTMP " ::: $splitFILES

#combine sequence counts for sRNAseq data
if [[ $TYPE == sRNAseq ]] || [[ $TYPE == sRNAseqIP ]]; then
  echo Sequence $splitFILES | tr ' ' '\t' | gzip >${FOLDER}GEO/sRNA_counts.txt.gz

  awk -v OFS="\t" -v FILES=$FILES '
  BEGIN{
    #read in the library names
    nFILES=split(FILES, splitFILES, /~/)
  }
  {
    #put sequences and their counts into an array indexed by sequence then library
    ARRAY[$1][$2]=$3
  }
  END{
    #for all sequences in array
    for(SEQ in ARRAY) {
      # > required for conversion into wide format later on
      print ">"SEQ
      #traverse through libraries by numbers to ensure same order all the times
      for(i=1; i<=nFILES; i++) {
        #if there is a count for this sequence in the current library print it, otherwise print 0
        if(splitFILES[i] in ARRAY[SEQ]) {
          print ARRAY[SEQ][splitFILES[i]]
        }else{
          print 0
        }
      }
    }
  }' ${locTMP}*seqCOUNTs.txt |
    #convert long format into wide format
    mawk -v RS=">" -v OFS="\t" '
    {
      if(NR>1) {
        gsub("\n", "\t")
        print 
      }
    }' | sort -k1,1 | sed 's/\t$//' | gzip >>${FOLDER}GEO/sRNA_counts.txt.gz

  SUM=$(md5sum ${FOLDER}GEO/sRNA_counts.txt.gz)
  echo sRNA_counts.txt.gz $SUM | tr ' ' '\t' >>${FOLDER}GEO/md5_Sum_COUNTS.tmp

elif [[ $TYPE == RNAseq ]]; then
  cat ${FOLDER}GeTMM_gene.txt ${FOLDER}TE_GeTMM_*.txt | gzip >${GEO}RNAseq_GeTMM.txt.gz
  SUM=$(md5sum ${GEO}RNAseq_GeTMM.txt.gz)
  echo RNAseq_GeTMM.txt.gz $SUM | tr ' ' '\t' >>${FOLDER}GEO/md5_Sum_COUNTS.tmp
fi

if [[ -f ${FOLDER}GEO/md5_Sum_RAW.tmp ]]; then
  sort -k1,1 ${FOLDER}GEO/md5_Sum_RAW.tmp | cut -f 1-2 >${FOLDER}GEO/md5_Sum_RAW.txt
  rm -rf ${FOLDER}GEO/md5_Sum_RAW.tmp
fi
if [[ -f ${FOLDER}GEO/md5_Sum_BW.tmp ]]; then
  sort -k1,1 ${FOLDER}GEO/md5_Sum_BW.tmp | cut -f 1-2 >${FOLDER}GEO/md5_Sum_BW.txt
  rm -rf ${FOLDER}GEO/md5_Sum_BW.tmp
fi
if [[ -f ${FOLDER}GEO/md5_Sum_COUNTS.tmp ]]; then
  sort -k1,1 ${FOLDER}GEO/md5_Sum_COUNTS.tmp | cut -f 1-2 >${FOLDER}GEO/md5_Sum_COUNTS.txt
  rm -rf ${FOLDER}GEO/md5_Sum_COUNTS.tmp
fi

###################################################################################################
#clean up
if [[ $DEBUG != Y ]]; then
  rm -rf "$locTMP"
fi
PROCESSED_TIME=$(echo -e "$(date "+%s")" "$TIME" | mawk '{ print ($1-$2)/60 }')
printf "create GEO output - processing_time=	${PROCESSED_TIME} \n" >>"${FOLDER}LOGs/time-log.txt"
