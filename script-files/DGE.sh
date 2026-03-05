#!/usr/bin/env bash
#$ -S /bin/bash

#$ -clear
#$ -e $JOB_NAME.e$JOB_ID.txt
#$ -o $JOB_NAME.o$JOB_ID.txt
#$ -cwd
#$ -q public.q
#$ -pe smp 1

#SBATCH --cpus-per-task=2
#SBATCH --mem=10g
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

Rscript ${SCRIPT_DIR}DGE.R TMP=$TMPdir refGENO=$refGENO GENO=${GENO} FOLDER=${FOLDER} DGE=$DGE Nexec=$Nexec

###################################################################################################
PROCESSED_TIME=$(echo -e $(date "+%s") "$TIME" | mawk '{ print ($1-$2)/60 }')
echo "DGE - processing_time=" "${PROCESSED_TIME}" >"${LOGs}time-log.txt"
echo $PROCESSED_TIME

exit
