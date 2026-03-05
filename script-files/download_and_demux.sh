#!/usr/bin/env bash
#$ -S /bin/bash

#$ -clear
#$ -e $JOB_NAME.e$JOB_ID.txt
#$ -o $JOB_NAME.o$JOB_ID.txt
#$ -cwd
#$ -q public.q
#$ -pe smp 2

#SBATCH --cpus-per-task=6
#SBATCH --mem=15g
#SBATCH -e "%x.e.%A-%a.txt"
#SBATCH -o "%x.o.%A-%a.txt"
#XX-SBATCH --qos=short-XX
#SBATCH --time=24:00:00

hostname

set -ux

###################################################################################################
#extract variables
VARI=$(echo "$1" | sed 's/,/\t/g;s/"//g')
eval "$VARI"
echo $1 | sed 's/,/\n/g'

TIME=$(date "+%s")
TIMEx=$TIME

if [[ $GRIDsystem == SLURM ]]; then CORES=$(( SLURM_CPUS_PER_TASK * @HYPER@ )); elif [[ -z ${NSLOTS+x} ]]; then CORES=4; else CORES=$NSLOTS; fi
echo $CORES
TIME=$(date "+%s")

####################################################################################################
#presetting of environment to preload modules and to prepare all files required
#load all functions
source ${SCRIPT_DIR}tools

locTMP=${TMPdir}TMP_download+demux/
locTMP=${locTMP}${SLURM_ARRAY_TASK_ID}/
  
mkdir -p $locTMP
rm -rf ${locTMP}/*/error.txt

export TMPDIR=$locTMP

####################################################################################################
#generate variables and local directory

if [[ -f ${TMPdir}URLs.txt ]]; then
  #if download and demux from NGS
  multiLIBRARY=$(cat ${TMPdir}URLs.txt | sed -n ${SLURM_ARRAY_TASK_ID}p | tr ' ' '\t' | cut -f 1)
  sampleIDs=$(cat ${TMPdir}URLs.txt | sed -n ${SLURM_ARRAY_TASK_ID}p | tr ' ' '\t' | cut -f 2)
  runIDs=$(cat ${TMPdir}URLs.txt | sed -n ${SLURM_ARRAY_TASK_ID}p | tr ' ' '\t' | cut -f 3)
  BCtype=$(cat ${TMPdir}URLs.txt | sed -n ${SLURM_ARRAY_TASK_ID}p | tr ' ' '\t' | cut -f 4)
  
  #create local tmp
  echo $locTMP $multiLIBRARY
  
  #create barcode -file
  echo "i7 i5 sRBC library_name" | tr ' ' '\t' >${locTMP}BC_file.txt
  
  cat ${TMPdir}URLs.txt | sed -n ${SLURM_ARRAY_TASK_ID}p |
  awk -v locTMP=${locTMP} -v fwADAPTOR=$fwADAPTOR '
    {
      n=split($3, splitLIBs, /;/)
      split($2, splitNGSID, /;/)
      for(i=1; i<=n; i++) {

        split(splitLIBs[i], splitID, /~/)
        
        #split Barcode tag into i7, i5 and sRBC
        split(splitID[1], splitBC,/:|i/)

        print "BC:Z:"splitBC[1],"B2:Z:"splitBC[2],"sRBC="splitBC[3],splitNGSID[i]"~"splitLIBs[i]
      }
    }' | tr ' ' '\t' | sort | uniq >>${locTMP}BC_file.txt
  
  NcurrLIBRARY=1
  libFOLDER=$LOGs
else
  #if bam-files are supplied for demultiplexing
  locTMP=${TMPdir}TMP_download+demux/
  mkdir -p $locTMP
  multiLIBRARY=$(echo $FILE | tr '~' '\t')
  BCtype="NA"
  
  #determine n of libraries to merge
  NcurrLIBRARY=$(echo $multiLIBRARY | wc -w)
    
  #create sRBC setup
  #create barcode -file
  echo "i7 i5 sRBC library_name" | tr ' ' '\t' >${locTMP}BC_file.txt
  awk -v i7=$BCi7 -v i5=$BCi5 -v sRBC=$sRBC -v NAME=$NAME -v OFS="\t" '
  BEGIN{
    n=split(i7,spliti7,/~/)
    if(n == 0){
      i7OUT="B2:Z:"
    }else{
      for(i=1; i<=n; i++){
        if(i==1){ 
          i7OUT="BC:Z:"spliti7[i]
        }else{
          i7OUT=i7OUT"~BC:Z:"spliti7[i]
        }
      }
    }
    n=split(i5,spliti5,/~/)
    if(n == 0){
      i5OUT="B2:Z:"
    }else{
      for(i=1; i<=n; i++){
        if(i==1){ 
          i5OUT="B2:Z:"spliti5[i]
        }else{
          i5OUT=i5OUT"~BC:Z:"spliti5[i]
        }
      }
    }
    n=split(sRBC,splitsRBC,/~/)
    if(n == 0){
      sRBCOUT="sRBC="
    }else{
      for(i=1; i<=n; i++){
        if(i==1){ 
          sRBCOUT="sRBC="splitsRBC[i]
        }else{
          sRBCOUT=sRBCOUT"~sRBC="splitsRBC[i]
        }
      }
    }
    print i7OUT, i5OUT, sRBCOUT, NAME
  }
  {
    x="a"
  }' ${locTMP}BC_file.txt >>${locTMP}BC_file.txt    
fi
mkdir -p $locTMP


#reverse complement sRBC BC due to different read-direction
if [[ $BCtype == *"Stark_3prime_BC_PROseq" ]]; then
  mawk -v OFS="\t" '
  BEGIN{
    c["A"] = "T"; c["C"] = "G"; c["G"] = "C"; c["T"] = "A" 
  }
  function revcomp(x,  i, o) {
    o = ""
    for(i = length; i > 0; i--){
      o = o c[substr(x, i, 1)] }
    return(o)

  }
  {
    if(NR == 1){
      print $0
    }else{
      split($3, splitBC, /=/)
      REV = revcomp(splitBC[2])

      #!turn of reverse complement to see if this is correct
      #@ $3=splitBC[1]"="REV
      #@ sub(splitBC[2], REV, $4)

      #remove i5 sequence as it is useless for this demultiplexing
      $2="B2:Z:"
      
      print $1,$2,$3,$4
    }
  }' ${locTMP}BC_file.txt > ${locTMP}BC_file.txt.tmp
  cp ${locTMP}BC_file.txt.tmp ${locTMP}BC_file.txt
fi

###################################################################################################
#download, and demultiplex to sam files

nCHUNKS=6
redCORES=$(($CORES - 1))
#wipe old output-bam-files
rm -rf ${locTMP}*~*.bam

#@ if [[ $BCtype == *"Stark_3prime_BC_PROseq" ]]; then
#@   parallel --tmpdir $TMPDIR -j 2 ". ${SCRIPT_DIR}functions; downFunct_Stark_3prime_BC_PROseq {1} {#} $locTMP $redCORES $TMPdir ${SLURM_ARRAY_TASK_ID} $FOLDER $downDIR $LOGs " ::: $multiLIBRARY
#@ else
  parallel --tmpdir $TMPDIR -j 2 ". ${SCRIPT_DIR}functions; downFunct {1} {#} $locTMP $redCORES $TMPdir ${SLURM_ARRAY_TASK_ID} $FOLDER $downDIR $LOGs $BCtype" ::: $multiLIBRARY
#@ fi

#sleep a bit to allow storage to catch up
sleep 15s

#@ #only run code if sRBC dmux required
#@ if [[ -f ${locTMP}allsRBC.txt ]]; then
#@   ###################################################################################################
#@   #demultiplex sRBC
#@   allBC=$(cat ${locTMP}BC_file.txt | grep sRBC | tr ' ' '\t' | cut -f 1 | tr '\n' '\t')
  
#@   redCORES=$(( $CORES / 2 ))
#@   parallel -j $redCORES ". ${SCRIPT_DIR}functions; demuxFunctSRBC {2}  $locTMP {1} {3} " ::: $allBC ::: $(seq 1 $NcurrLIBRARY)  ::: $BCtype
  
#@   #---------------------------------------------------------------------------------------------------------
#@   #split bam by sRBC for next steps
#@   allBC=$(cat ${locTMP}allsRBC.txt | tr ' ' '\t' | cut -f 1 | tr '\n' '\t')
  
#@   parallel -j $redCORES ". ${SCRIPT_DIR}functions; extractFunctSRBC {2}  $locTMP {1} {3} " ::: $allBC ::: $(seq 1 $NcurrLIBRARY)  ::: $BCtype

#@   #test if non-sRBC libraries have been included in the import and convert these to bam
#@   convertBC=$(cat ${locTMP}BC_file.txt | grep -v sRBC | tail -n +2 | tr ' ' '\t' | cut -f 1 | tr '\n' '\t')
#@   echo in
#@   if [[ ! -z $convertBC ]]; then
#@     time parallel -j $CORES "inNAME=${locTMP}*/{}.sam; outNAME=\$(echo \$inNAME | sed 's/.sam/.bam/'); echo \$outNAME ; samtools view -b \$inNAME > \$outNAME" ::: $convertBC
#@   fi
#@ else
  #@ #convert demultiplexed sam to bam if no sRBC basecalling is required
  #@ time parallel --tmpdir $TMPDIR -j $redCORES "outNAME=\$(echo {} | sed 's/.sam$/.bam/'); echo \$outNAME ; samtools view -b -@ 2 {} > \$outNAME" ::: ${locTMP}*/*.sam
#@ fi

if [[ -s ${locTMP}/*/error.txt ]]; then
  cat ${locTMP}/*/error.txt >> ${FOLDER}error_download.txt
  exit
fi

PROCESSED_TIME=$(echo -e $(date "+%s") $TIME | awk '{ print ($1-$2)/60 }')
echo  ${PROCESSED_TIME}

###################################################################################################
#merge demultiplexed chunks
if [[ -f ${TMPdir}URLs.txt ]]; then
#!could be stped up by parallel
   tail -n +2 ${locTMP}BC_file.txt > ${locTMP}BC_file.noHEAD.txt
  while read LINE; do
    FILENAME=$(echo $LINE | tr ' ' '\t' | cut -f 4)
    
    samtools cat -o ${LIB_STORAGE_FOLDER}${FILENAME}.bam.${SLURM_JOB_ID}.tmp ${locTMP}*/${FILENAME}.out.bam
    
    newSUM=$(md5sum ${LIB_STORAGE_FOLDER}${FILENAME}.bam.${SLURM_JOB_ID}.tmp | tr ' ' '\t' | cut -f 1)
    echo ${LIB_STORAGE_FOLDER}${FILENAME}.bam $newSUM >${LIB_STORAGE_FOLDER}${FILENAME}.bam.md5.${SLURM_JOB_ID}.tmp
    
    if [[ (! -f ${LIB_STORAGE_FOLDER}${FILENAME}.bam || $FORCEimport == Y) && -f ${LIB_STORAGE_FOLDER}${FILENAME}.bam.${SLURM_JOB_ID}.tmp ]]; then
      mv ${LIB_STORAGE_FOLDER}${FILENAME}.bam.${SLURM_JOB_ID}.tmp ${LIB_STORAGE_FOLDER}${FILENAME}.bam
      mv ${LIB_STORAGE_FOLDER}${FILENAME}.bam.md5.${SLURM_JOB_ID}.tmp ${LIB_STORAGE_FOLDER}${FILENAME}.bam.md5
    else
      rm -rf ${LIB_STORAGE_FOLDER}${FILENAME}.bam.${SLURM_JOB_ID}.tmp
      rm -rf ${LIB_STORAGE_FOLDER}${FILENAME}.bam.md5.${SLURM_JOB_ID}.tmp
    fi
  done < ${locTMP}BC_file.noHEAD.txt
else
  samtools cat -o ${TMPdir}${NAME}_demuxed.bam ${locTMP}*/${NAME}.out.bam
fi

#---------------------------------------------------------------------------------------------------------
if [[ $DEBUG != Y ]]; then
  rm -rf ${locTMP}
fi

PROCESSED_TIME=$(echo -e $(date "+%s") $TIME | awk '{ print ($1-$2)/60 }')
echo "download and demux - processing_time=" "${PROCESSED_TIME}" >>"${libFOLDER}time-log.txt"
echo $PROCESSED_TIME

exit
