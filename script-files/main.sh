#!/usr/bin/env bash
set -ux

###################################################################################################
###################################################################################################
#extract variables
VARI=$(echo $1 | sed 's/,/\t/g;s/"//g')
eval $VARI
VARI=$1
echo $VARI | sed 's/,/\n/g'

###################################################################################################
##generate storage location varaibles and create folders

#prepare storage and COMPUTING variables
libFOLDER=${FOLDER}/individual-libraries/${NAME}/
TMPdir=${rawTMP}/${NAME}/
VARI=${VARI},libFOLDER=${libFOLDER},TMPdir=${TMPdir}

#---------------------------------------------------------------------------------------------------------
#setting of time variables
TIME=$(date "+%s")
TIMEall=$(date "+%s")

#---------------------------------------------------------------------------------------------------------
#folder preperation

#create folders
mkdir -p ${libFOLDER}/
mkdir -p ${TMPdir}
mkdir -p $FOLDER

rm -rf ${libFOLDER}log.txt
rm -rf ${libFOLDER}time-log.txt
rm -rf ${TMPdir}jobIDs.txt

{
  printf "$NAME\n\n"
  printf "type= ${TYPE} $IPstatus $OXstatus\n\n"
  echo output folder = $FOLDER
  echo tmp folder = $TMPdir
  printf "\n"
} >>${libFOLDER}log.txt

rm -rf ${libFOLDER}LOGs/
mkdir -p ${libFOLDER}LOGs/
cd ${libFOLDER}LOGs/
LOGs=${libFOLDER}LOGs/
VARI=${VARI},LOGs=${LOGs}

###################################################################################################
#submission of the individual subjobs

#---------------------------------------------------------------------------------------------------------
#test if all files are present; if not wait and test again --> if available store md5 sum for comparison later
if [[ $NGS == Y ]]; then
  RELEASE=N
  rm -rf ${TMPdir}old.md5
  while [[ $RELEASE == N ]]; do
    RELEASE=Y
    for currFILE in $(echo $FILE | sed 's/:!:/\t/g'); do
      if [[ ! -f $currFILE ]]; then
        RELEASE=N
      fi
    done
    if [[ $RELEASE != Y ]]; then
      sleep 3s
    fi
  done
  
  #add sleep to ensure complete file transfer
  #@sleep 60s
  
  for currFILE in $(echo $FILE | sed 's/:!:/\t/g'); do
    cat ${currFILE}.md5 >>${TMPdir}old.md5
  done
  
fi
if [[ -f ${FOLDER}error_download.txt ]] ; then
  printf "previous download error detected please check demultiplexing log\n" >>${FOLDER}error.txt
  exit
fi

#---------------------------------------------------------------------------------------------------------
#demultiplex local multiplexed file
if [[ ( $BCi7 != "" || $sRBC != "" ) && $NGS != Y && $demuxFASTA != Y ]]; then
  COMMAND=${SCRIPT_DIR}download_and_demux.sh
  
  if [[ $COMPUTING == C ]]; then
    if [[ $GRIDsystem == SLURM ]]; then sbatch --wait $COMMAND ${VARI},SLURM_ARRAY_TASK_ID=1; else qsub -sync y $COMMAND ${VARI}; fi
  else
    $COMMAND ${VARI}
  fi
  
  FILE=${TMPdir}${NAME}_demuxed.bam
  BAM=Y
  VARI="$VARI,FILE=${FILE},BAM=${BAM}"
  
  #added sleep to prevent problems if storage is too slow
  sleep 30s
  
  PROCESSED_TIME=$(echo -e "$(date "+%s")" "$TIME" | awk '{ print ($1-$2)/60 }')
  printf "download and demux - total_time=	${PROCESSED_TIME} \n\n" >>"${libFOLDER}time-log.txt"
  TIME="$(date "+%s")"
fi

#---------------------------------------------------------------------------------------------------------
#bam --> fasta
if [[ $BAM == Y ]] || [[ $FORCE == Y ]] || [[ $SRR == Y ]] || [[ $NGS == Y ]] || [[ $demuxFASTA == Y ]]; then

  COMMAND="${SCRIPT_DIR}bam_to_fasta.sh"

  # set the maximum number of allowed repeats
  MAX_REPEATS=6

  # initialize the counter variable
  REPEATS=0

  while [[ ! -f ${TMPdir}${NAME}_raw.fq ]]; do
    # increment the counter variable
    REPEATS=$((REPEATS + 1))

    if [[ $COMPUTING == C ]]; then
      if [[ $GRIDsystem == SLURM ]]; then sbatch --wait $COMMAND ${VARI}; else qsub -sync y $COMMAND ${VARI}; fi
    else
      $COMMAND ${VARI}
    fi

    # check if the maximum number of allowed repeats has been reached
    if [[ $REPEATS -gt $MAX_REPEATS ]]; then
      printf "bam-->fasta failed to produce output after $MAX_REPEATS attempts for library ${NAME}\n" >>${FOLDER}error.txt
      exit
    fi

    if [[ $COMPUTING == C ]]; then
      if [[ $GRIDsystem == SLURM ]]; then sbatch --wait $COMMAND ${VARI}; else qsub -sync y $COMMAND ${VARI}; fi
    else
      $COMMAND ${VARI}
    fi

    #test if cutadapt caused problems 
    TEST=$(grep -l ERRO  ${LOG}bam_to_fasta.sh.e.*.txt )
    if [[ -n $TEST ]]; then
      printf "cutadapt failed and output truncated in ${NAME}\n" >>${FOLDER}error.txt
      exit
    fi

    #look up md5 sum in the storage after the bam was converted to fasta
    if [[ $NGS == Y ]]; then
      rm -rf ${TMPdir}new.md5
      for currFILE in $(echo $FILE | sed 's/:!:/\t/g'); do
        echo ${currFILE}.md5
        cat ${currFILE}.md5 >>${TMPdir}new.md5
      done
      
      #added sleep to prevent problems if storage is too slow
      sleep 60s
      
      #test if md5 changed and if so terminate
      DIFF=$(diff ${TMPdir}new.md5 ${TMPdir}old.md5)
      if [[ -n $DIFF ]]; then
        printf "md5 sum changed while bam_to_fasta was running ${NAME}\n" >>${FOLDER}error.txt
        exit
      fi
    else
      #added sleep to prevent problems if storage is too slow
      sleep 60s
    fi
    
    #change file variable to the fasta file
    FILE=${TMPdir}${NAME}_raw.fq
    VARI="$VARI,FILE=${FILE}"
    
    #sanity check output file
    if [[ ! -f $FILE ]]; then
      printf "bam-->fasta did not work properly for library ${NAME}\n" >>${FOLDER}error.txt
      exit
    fi
  done

  PROCESSED_TIME=$(echo -e "$(date "+%s")" "$TIME" | awk '{ print ($1-$2)/60 }')
  printf "bam to fasta - total_time=	${PROCESSED_TIME} \n\n" >>"${libFOLDER}time-log.txt"
  TIME="$(date "+%s")"
fi

#if only raw fasta requested transfer file and exit
if [[ $RAW == Y ]]; then
  mkdir ${FOLDER}fasta-files/
  sbatch --wait --time=2:00:00 --job-name="transfer-file" --mem=1g --wrap="gzip -c -v -f $FILE > ${FOLDER}fasta-files/${NAME}.fa.gz"
  #remove complete TMPdir
  if [[ $keepTMP == N ]]; then
    rm -rf $TMPdir
  fi
  exit
fi

#---------------------------------------------------------------------------------------------------------
#peperation of reads
COMMAND="${SCRIPT_DIR}prepare_reads.sh"

# set the maximum number of allowed repeats
MAX_REPEATS=6

# initialize the counter variable
REPEATS=0
# initialize the TEST variable
TEST=""

while [[ -z $TEST ]]; do
  # increment the counter variable
  REPEATS=$((REPEATS + 1))

  # check if the maximum number of allowed repeats has been reached
  if [[ $REPEATS -gt $MAX_REPEATS ]]; then
    printf "prepare_reads failed to produce output after $MAX_REPEATS attempts for library ${NAME}\n" >>${FOLDER}error.txt
    exit
  fi

  if [[ $COMPUTING == C ]]; then
    if [[ $GRIDsystem == SLURM ]]; then sbatch --wait $COMMAND ${VARI}; else qsub -sync y $COMMAND ${VARI}; fi
  else
    $COMMAND ${VARI}
  fi

  #added sleep to prevent problems if storage is too slow
  sleep 30s

  #sanity check output file
  TEST=$(gunzip -c ${TMPdir}${NAME}.fa.gz | head -n 1)
  #if [[ -z $TEST ]]; then
  #  printf "prepare_reads did not work properly for library ${NAME}\n" >>${FOLDER}error.txt
  #  exit
  #fi
done



#added sleep to prevent problems if storage is too slow
sleep 30s

#sanity check output file
TEST=$(gunzip -c ${TMPdir}${NAME}.fa.gz | head -n 1)
if [[ -z $TEST ]]; then
  printf "prepare_reads did not work properly for library ${NAME}\n" >>${FOLDER}error.txt
  exit
fi

PROCESSED_TIME=$(echo -e "$(date "+%s")" "$TIME" | awk '{ print ($1-$2)/60 }')
printf "prepare_reads - total_time=	${PROCESSED_TIME} \n\n" >>"${libFOLDER}time-log.txt"
TIME="$(date "+%s")"

#----------------------------------------------------------------------------------------------------
#map reads against the genome
if [[ $TYPE == "RNAseq" ]] || [[ $TYPE == RIPseq ]] || [[ $TYPE == GROseq ]]; then
  COMMAND="${SCRIPT_DIR}map_genome_RNAseq.sh"
  elif [[ $TYPE == "CLIPseq" ]]; then
  COMMAND="${SCRIPT_DIR}map_genome_CLIPseq.sh"
else
  COMMAND="${SCRIPT_DIR}map_genome.sh"
fi

if [[ $COMPUTING == C ]]; then
  if [[ $GRIDsystem == SLURM ]]; then sbatch --wait $COMMAND ${VARI}; else qsub -sync y $COMMAND ${VARI}; fi
else
  $COMMAND ${VARI}
fi

PROCESSED_TIME=$(echo -e "$(date "+%s")" "$TIME" | awk '{ print ($1-$2)/60 }')
printf "map reads - total_time=	${PROCESSED_TIME} \n\n" >>"${libFOLDER}time-log.txt"
TIME="$(date "+%s")"

#-----------------------------------------------------------------------------------------------------
#intersect mapping with annotations and apply ruleset

COMMAND="${SCRIPT_DIR}intersect_mapping_with_annotations-apply_rules.sh "

if [[ $COMPUTING == C ]]; then
  if [[ $GRIDsystem == SLURM ]]; then sbatch --wait $COMMAND ${VARI}; else qsub -sync y $COMMAND ${VARI}; fi
else
  $COMMAND ${VARI}
fi

PROCESSED_TIME=$(echo -e "$(date "+%s")" "$TIME" | awk '{ print ($1-$2)/60 }')
printf "intersect_with_annotations and apply ruleset - total_time=	${PROCESSED_TIME} \n\n" >>"${libFOLDER}time-log.txt"
TIME="$(date "+%s")"

#---------------------------------------------------------------------------------------------------------
#collapse multimappers and extract numbers for the charts

COMMAND="${SCRIPT_DIR}finish_counts.sh"

if [[ $COMPUTING == C ]]; then
  if [[ $GRIDsystem == SLURM ]]; then sbatch --wait $COMMAND ${VARI}; else qsub -sync y $COMMAND ${VARI}; fi
else
  $COMMAND ${VARI}
fi

PROCESSED_TIME=$(echo -e "$(date "+%s")" "$TIME" | awk '{ print ($1-$2)/60 }')
printf "collapse annotations for multi-mappers and finish counts - total_time=	${PROCESSED_TIME} \n\n" >>"${libFOLDER}time-log.txt"
TIME="$(date "+%s")"

###################################################################################################
#run various modules in parallel if required

#generate proper normalization factor for norm.txt
if [[ $noNORM == Y ]]; then
  NORMfactor=1
else
  if [[ $TYPE == sRNAseq ]] && [[ $IPstatus != IP ]] && [[ $OXstatus != OX ]] && [[ $WIG_FASTA != Y ]] && [[ $WIG != Y ]]; then
    #for 1st track put miRNA normalization count into file
    LASTline=$(tail -n 1 <"${libFOLDER}normalization.txt")
    if [[ ! "$LASTline" =~ miRNA ]]; then
      miRNAfactor=$(awk '{ if ($1 ~ "miRNA") X+=$2 } END{ print X/1000000 }' <"${libFOLDER}${NAME}_annotation_counts.txt")
      printf "${miRNAfactor} - normalized to 1000000 miRNAs \n" >>"${libFOLDER}normalization.txt"
    fi
  elif [[ "$IPstatus" == IP ]] || [[ "$OXstatus" == OX ]] || [[ "$WIG_FASTA" == Y ]]; then
    LASTline=$(tail -n 1 <"${libFOLDER}normalization.txt")
    if [[ ! "$LASTline" =~ fasta ]]; then
      X=$(awk '{if($0 ~ "normalized to 10000000 reads in fasta after filtering") print }' ${libFOLDER}normalization.txt)
      echo "$X" >>"${libFOLDER}normalization.txt"
    fi
  else
    LASTline=$(tail -n 1 <"${libFOLDER}normalization.txt")
    
    if [[ ! "$LASTline" =~ mapping ]]; then
      X=$(awk '{if($0 ~ "normalized to 10000000 uniquely mapping reads") print }' "${libFOLDER}normalization.txt")
      echo "$X" >>"${libFOLDER}normalization.txt"
    fi
  fi
  NORMfactor=$(tail -n 1 <"${libFOLDER}normalization.txt" | tr ' ' '\t' | cut -f 1)

fi

if [[ $spikeINnorm == Y ]]; then
  NORMfactor=$(cat ${TMPdir}spike-in.norm )
  printf "${NORMfactor} - spike-in based normalization factor \n" >>"${libFOLDER}normalization.txt"
fi

#if external normalisation factor is supplied add it to the normalization.txt file and reset NORMfactor
if [[ ! -z $extNORM ]]; then
  NORMfactor=$extNORM
  printf "${extNORM} - externally supplied normalization factor during submission \n" >>"${libFOLDER}normalization.txt"
fi

#add normalization factor to variables
VARI=${VARI},NORMfactor=${NORMfactor}

#-------------------------------------------------------------------------------------------------------
#if RNAseq run mapping free quant
if [[ $TYPE == RNAseq ]] || [[ $TYPE == sRNAseq ]] || [[ $TYPE == RIPseq ]] || [[ $TYPE == GROseq ]] || [[ $TYPE == CLIPseq ]]||[[ $FORCEquant_unstranded == Y ]]; then
  rm -rf ${TMPdir}done.txt
  
  COMMAND="${SCRIPT_DIR}mapping_free_quant.sh"
  if [[ $COMPUTING == C ]]; then
    if [[ $GRIDsystem == SLURM ]]; then { sbatch --wait $COMMAND ${VARI} &}; else { qsub -sync y $COMMAND ${VARI} &}; fi
  else
    $COMMAND ${VARI}
  fi
  
else
  touch ${TMPdir}done.txt
fi

#-------------------------------------------------------------------------------------------------------
#create normalized SNP-table for SLAMseq
if [[ $SLAM == Y ]]; then
  awk -v NORMfactor=$NORMfactor -v libFOLDER=${libFOLDER} '{
    ALL[$1]+=($4*NORMfactor)
    
    if($2~"final_ANN=miRNA"){
      miRNA[$1]+=($4*NORMfactor)
    }
    if($3>23 && $3<=35){
      piRNA[$1]+=($4*NORMfactor)
    }
  }
  END{
    X="A-C-G-T"
    split(X,Y,/-/)

    for(i in Y){
      for(j in Y){
        Z=Y[i]">"Y[j]
        print Z
        if(Z in ALL){
          print Z,ALL[Z] > libFOLDER "SNPsummary.normalized.all.txt"
        }else{ print Z,0> libFOLDER "SNPsummary.normalized.miRNA.txt"}

        if(Z in miRNA){
          print Z,miRNA[Z] > libFOLDER "SNPsummary.normalized.miRNA.txt"
        }else{ print Z,0 > libFOLDER "SNPsummary.normalized.miRNA.txt"}

        if(Z in piRNA){
          print Z,piRNA[Z] > libFOLDER "SNPsummary.normalized.piRNA.txt"
        }else{ print Z,0 > libFOLDER "SNPsummary.normalized.piRNA.txt"}
      }
    }
  }' ${TMPdir}SNPlist.txt | LC_COLLATE=C sort -k1,1 > ${libFOLDER}SNPsummary.normalized.txt
fi

#-------------------------------------------------------------------------------------------------------
#run wig generateion - multiple jobs if pi/siRNA tracks get generated
rm -rf "${FOLDER}/hubs/standard_WIG/"
n_track_types="1"
if [[ "$TYPE" == sRNAseq ]] && [[ "$IPstatus" != IP ]]; then n_track_types="1 2 3"; fi
if [[ "$TYPE" == CapSeq ]]; then n_track_types="1 4"; fi
if [[ "$TYPE" == RNAseq ]] || [[ $TYPE == RIPseq ]] || [[ $TYPE == GROseq ]] ; then n_track_types="5"; fi
if [[ $TYPE == CLIPseq ]]; then n_track_types="5 6"; fi

#submit WIG-track generation
COMMAND="${SCRIPT_DIR}create_wig.sh"
for TRACK_TYPE in $(echo "${n_track_types}" | tr ' ' '\t'); do
  if [[ $COMPUTING == C ]]; then
    if [[ $GRIDsystem == SLURM ]]; then { sbatch --wait $COMMAND ${VARI},TRACK_TYPE_NR=${TRACK_TYPE} &}; else { qsub -sync y $COMMAND ${VARI},TRACK_TYPE_NR=${TRACK_TYPE} &}; fi
  else
    $COMMAND ${VARI},TRACK_TYPE_NR=${TRACK_TYPE}
  fi
done

#-------------------------------------------------------------------------------------------------------

if [[ $TYPE == RNAseq ]] || [[ $TYPE == sRNAseq ]] || [[ $TYPE == sRNAseqIP ]] || [[ $TYPE == CHIPseq ]] || [[ $TYPE == RIPseq ]] || [[ $TYPE == CLIPseq ]]|| [[ $TYPE == DNAseq ]] || [[ $TYPE == GROseq ]] ; then
  COMMAND="${SCRIPT_DIR}generate_TEbg.sh"
  
  if [[ $COMPUTING == C ]]; then
    if [[ $GRIDsystem == SLURM ]]; then { sbatch --wait $COMMAND ${VARI} &}; else { qsub -sync y $COMMAND ${VARI} &}; fi
  else
    $COMMAND ${VARI}
  fi
fi

#-------------------------------------------------------------------------------------------------------
#compute 1kb tile table
if [[ $TYPE == RNAseq ]] || [[ $TYPE == sRNAseq ]] || [[ $TYPE == sRNAseqIP ]] || [[ $TYPE == CHIPseq ]] || [[ $TYPE == RIPseq ]] || [[ $TYPE == CLIPseq ]]; then
  COMMAND="${SCRIPT_DIR}1kb-tiles.sh"
  
  if [[ $COMPUTING == C ]]; then
    if [[ $GRIDsystem == SLURM ]]; then { sbatch --wait $COMMAND ${VARI} &}; else { qsub -sync y $COMMAND ${VARI} &}; fi
  else
    $COMMAND ${VARI}
  fi
fi

#-------------------------------------------------------------------------------------------------------
#pause scritp until all sub-modules are finished
wait

###################################################################################################
#end of script

#test if SLURM gave an out of memory error
if [[ $GRIDsystem == SLURM ]]; then
  for FILE in $(ls ${libFOLDER}LOGs/*.sh.e.*); do
    X=$(grep "Exceeded job memory limit at some point" $FILE)
    Y=$(grep "Killed" $FILE)
    if [[ ! -z ${X} ]] && [[ ! -z ${Y} ]]; then
      jobID=$(echo $FILE | awk '
      {
        n=split($1,X,/\./)
        print X[n-1], $1
      } ')
      STATE=$(sacct -j $jobID | awk '{if(NR>2) { if($NF != "0:0") X="problem" }} END{print X}')
      #    if [[ ! -z $STATE ]]
      #    then
      for i in 1 2 3 4 5 v w x y z; do
        echo $FILE >>"${FOLDER}${i}_!!!!!memory_exeeded in $NAME - ask Dominik for help!!!!!"
        sacct -j $jobID >>"${FOLDER}${i}_!!!!!memory_exeeded in $NAME - ask Dominik for help!!!!!"
      done
      #    fi
    fi
  done
fi

#test if there was a network interruption
if [[ $GRIDsystem == SLURM ]]; then
  TEST=$(cat ${FOLDER}/LOGs/${NAME}.log | grep "Unable to establish control machine address")
  if [[ ! -z $TEST ]]; then
    printf "notwork interruption - be careful and check results ${NAME}\n" >>${FOLDER}error_network.txt
  fi
fi

#rm -fR -- ${TMPdir}*/

#report total time for analysis
PROCESSED_TIME=$(echo -e "$(date "+%s")" "$TIMEall" | awk '{ print ($1-$2)/60 }')
printf "\ntotal time for analysis - total_time=	${PROCESSED_TIME} \n\n" >>"${libFOLDER}time-log.txt"
