#!/usr/bin/env bash
#$ -S /bin/bash

#$ -clear
#$ -e $JOB_NAME.e$JOB_ID.txt
#$ -o $JOB_NAME.o$JOB_ID.txt
#$ -cwd
#$ -q public.q
#$ -pe smp 6

#SBATCH --cpus-per-task=8
#SBATCH --mem=5g
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
locTMP=${TMPdir}TMP_intersect/

mkdir -p $locTMP
chmod 777 $locTMP

export TMPDIR=$locTMP

###################################################################################################
#intersect mappings with the annotation and split into sense and antisense

mkdir ${TMPdir}split_intersect/
nSPLITS=$(ls ${TMPdir}split_mapping | grep bed | wc -w )
RANGE=$(seq 1 $nSPLITS)

#run intersection function
redCORES=$(( CORES / 3 ))
parallel --tmpdir $TMPDIR -j $redCORES ". ${SCRIPT_DIR}functions; intersect_mappings $TMPdir $NAME ${UTILITY_LOCATION}${VERSION}.bed $UTILITY_DIR {}" ::: $RANGE 
    
#determine # of mappings that get collapsed because they share start and end coordinate
Noriginal=$( cat ${TMPdir}split_mapping/${NAME}_mapped*.bed | wc -l )
Ncollapsed=$( cat ${TMPdir}split_intersect/${NAME}_annotated_*.txt | wc -l )
N=$( echo $Noriginal $Ncollapsed | awk '{print $1-$2}' )

printf "

# of mappings that got collapsed because they originate from the very same read
that got mapped into the same location but had different indel assignment= $N
" >> ${libFOLDER}log.txt

#---------------------------------------------------------------------------------------------------------
#cleanup
if [[ $DEBUG != Y ]]; then
  rm -rf $locTMP
fi

PROCESSED_TIME=$(echo -e "$(date "+%s")" $TIME | mawk '{ print ($1-$2)/60 }' )
printf "intersect_with_annotations and apply ruleset - processing_time= ${PROCESSED_TIME} \n">> ${libFOLDER}time-log.txt

exit

