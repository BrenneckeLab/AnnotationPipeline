#!/usr/bin/env bash
#$ -S /bin/bash

#$ -clear
#$ -e $JOB_NAME.e$JOB_ID.txt
#$ -o $JOB_NAME.o$JOB_ID.txt
#$ -cwd
#$ -q public.q
#$ -pe smp 1

#SBATCH --cpus-per-task=1
#SBATCH --mem=2g
#SBATCH --time=2-00:00:00
#SBATCH -e "%x.e.%j.txt"
#SBATCH -o "%x.o.%j.txt"
#XX-SBATCH --qos=short-XX
@SBATCH@

hostname

set -ux

###################################################################################################
#extract variables
VARI=$(echo "$1" | sed 's/,/\t/g;s/"//g')
eval "$VARI"
VARI=$1
echo "$1" | sed 's/,/\n/g;s/"//g'

TIME=$(date "+%s")
TIMEx=$TIME

#presetting of environment to preload modules and to prepare all files required
#load all functions
source ${SCRIPT_DIR}tools

###################################################################################################

COLOR=$(echo "$subCOLOR" | tr '~' ',')

rm -rf ${FOLDER}*XXXXXXXXXXXXX*
rm -rf ${FOLDER}error.txt

###################################################################################################
#prepare ref required

if [[ ${prepareREF} == yes ]]; then
  if [[ $GENOME_VERSION == dm* ]]; then
    COMMAND="${SCRIPT_DIR}prepare_ref.sh"
  elif [[ $GENOME_VERSION == ASM ]]; then
    COMMAND=${SCRIPT_DIR}prepare_ref.asm.sh
  elif [[ $GENOME_VERSION == Cel* ]]; then
    COMMAND=${SCRIPT_DIR}prepare_ref.Cel.sh
  fi
  if [[ $GRIDsystem == SLURM && $COMPUTING == C ]]; then sbatch --wait $COMMAND ${VARI}; else $COMMAND ${VARI}; fi >${FOLDER}LOGs/prepare_ref.log.txt
fi

#---------------------------------------------------------------------------------------------------------
#copy chromosome size-file to logs-folder for access by other pipelines
cp ${UTILITY_LOCATION}chrom.sizes ${LOGs}

#copy sequencedSamples backup
cp ${TMPdir}sequencedSamples.txt ${LOGs}sequencedSamples.backup.txt

#---------------------------------------------------------------------------------------------------------
#prepare TE indexes if required
if [[ ! -z $extraSEQ ]]; then
  seqkit grep -v -p stellate -p Su-Ste --line-width 0 ${UTILITY_LOCATION}TE_annot_new_WO_DUST.fa >${TMPdir}TE_noSplice.fa
  cat ${extraSEQ} >>${TMPdir}TE_noSplice.fa

  bowtieBuild --threads 1 ${TMPdir}TE_noSplice.fa ${TMPdir}TE_index
  seqkit fx2tab --name --length ${TMPdir}TE_noSplice.fa >${TMPdir}TE.chrom.sizes
 fi

###################################################################################################
#create links for ASM version if they do not exist
if [[ $GENOME_VERSION == ASM && ! -s ${GENOMEdir}hub.txt ]]; then
  cp -rl ${ASMdir}/* $GENOMEdir
  unlink ${GENOMEdir}hub.txt
  cat ${ASMdir}/hub.txt |
    grep -v "longLabel "  > ${GENOMEdir}hub.txt
  echo longLabel AP-collection: ${COLLECTIONname} >> ${GENOMEdir}hub.txt

  unlink ${GENOMEdir}genomes.txt
  cat ${ASMdir}/genomes.txt |
    grep -v "description " |
    grep -v "genome "  > ${GENOMEdir}genomes.txt
  
  echo genome ${VERSION}_AP-${COLLECTIONname} >> ${GENOMEdir}genomes.txt
  echo description ${VERSION} - AP-collection: ${COLLECTIONname} - >> ${GENOMEdir}genomes.txt
  unlink ${GENOMEdir}/${VERSION}/trackDb.txt
  cp ${ASMdir}/${VERSION}/trackDb.txt ${GENOMEdir}/${VERSION}/trackDb.txt
fi

###################################################################################################
#download NGS files if required
COMMAND="${SCRIPT_DIR}download_and_demux.sh"

if [[ -f ${TMPdir}URLs.txt ]]; then
  nDOWNLOADS=$(cat ${TMPdir}URLs.txt | wc -l)
  mkdir -p ${TMPdir}TMP_download+demux/
  if [[ $GRIDsystem == SLURM && $COMPUTING == C ]]; then sbatch --array=1-${nDOWNLOADS}%12 $COMMAND "${VARI}"; else $COMMAND "${VARI},SLURM_ARRAY_TASK_ID=1"; fi &>${FOLDER}LOGs/download-files.log.txt
fi

if [[ ${DEMUXonly} == Y ]]; then
  touch ${FOLDER}___---DEMULTIPLEXING-ONLY---___
  exit
fi

###################################################################################################
#run main script to start analysis
TIMEall=$(date "+%s")
INVERTold=$INVERT


for i in $(seq 1 "${nFILES}"); do
  FILE=$(sed -n "${i}"p "${TMPdir}files.txt" | awk '{ print $1 }' | tr ';' '~' )
  NAME=$(sed -n "${i}"p "${TMPdir}files.txt" | awk '{ print $2 }')
  #replace tab with space for proper BC lookup
  LINE=$(sed -n "${i}"p "${TMPdir}files.txt" | tr '\t' ' ' )
  echo $i $LINE
  
  if [[ ${FILE} == NGS* || ${FILE} == TGZ*  ]]; then NGS=Y; else NGS=N; fi
  if [[ ${FILE} == *.bam ]]; then BAM=Y; else BAM=N; fi
  if [[ ${FILE} =~ ^(SRR|ERR|DRR|SRX) ]]; then
    SRR=Y
    SRAlink=$(echo $LINE | tr ' ' '\t' | tr '\t' '\n' | grep SRAlink)
    echo $SRAlink
  else
    SRR=N
    SRAlink=""
  fi
  
  IPstatus=$(echo $LINE | awk -v line=$i '{ for(i=3; i<=NF; i++) {if($i == "IP") { print "IP" }}}')
  OXstatus=$(echo $LINE | awk -v line=$i '{ for(i=3; i<=NF; i++) {if($i == "OX") { print "OX" }}}')
  
  if [[ $LINE == *" INVERT"* ]]; then INVERT=Y; else INVERT=$INVERTold; fi
  #extract sub-sampling from sample submission
  if [[ $LINE == *" SUB="* ]]; then
    currSUB=$(echo $LINE |
    awk '{ for(i=1; i<=NF; i++){ if($i ~ "SUB=") {split($i,splitSUB,/=/); print splitSUB[2];exit } } }')
  else
    currSUB=$SUBSAMPLE
  fi
  #enable external supplied normalization value
  if [[ $LINE == *" NORM="* ]]; then
    extNORM=$(echo $LINE |
    awk '{ for(i=1; i<=NF; i++){ if($i ~ "NORM=") {split($i,splitNORM,/=/); print splitNORM[2];exit } } }')
  else
    extNORM=""
  fi
  #enable external supplied normalization value
  if [[ $LINE == *" spikeIN"* ]]; then
    spikeINnorm=Y
  fi
  #detect barcodes
  if [[ $LINE == *" BCi7="*  && $LINE != *" BCi7=; "* ]]; then
    BCi7=$(echo $LINE |
    awk '{ for(i=1; i<=NF; i++){ if($i ~ "BCi7=") {split($i,splitBC,/=/); print splitBC[2];exit } } }' | tr ';' '~')
  else
    BCi7=""
  fi
  if [[ $LINE == *" BC="*  && $LINE != *" BC=; "* ]]; then
    BCi7=$(echo $LINE |
    awk '{ for(i=1; i<=NF; i++){ if($i ~ "BC=") {split($i,splitBC,/=/); print splitBC[2];exit } } }' | tr ';' '~')
  fi
  if [[ $LINE == *" BCi5="*  && $LINE != *" BCi5=; "* ]]; then
    BCi5=$(echo $LINE |
    awk '{ for(i=1; i<=NF; i++){ if($i ~ "BCi5=") {split($i,splitBC,/=/); print splitBC[2];exit } } }' | tr ';' '~')
  else
    BCi5=""
  fi
  if [[ $LINE == *" sRBC="*  && $LINE != *" sRBC=; "* ]]; then
    sRBC=$(echo $LINE |
    awk '{ for(i=1; i<=NF; i++){ if($i ~ "sRBC=") {split($i,splitRBC,/=/); print splitRBC[2];exit } } }' | tr ';' '~')
  else
    sRBC=""
  fi
  if [[ $LINE == *"UMIi7="*  ]]; then
    UMIi7=$(echo $LINE |
    awk '{ for(i=1; i<=NF; i++){ if($i ~ "UMIi7=") {split($i,splitUMI,/=/); print splitUMI[2];exit } } }' | tr ';' '~')
  else
    UMIi7=0
  fi
  if [[ $LINE == *"UMIi5="*  ]]; then
    UMIi5=$(echo $LINE |
    awk '{ for(i=1; i<=NF; i++){ if($i ~ "UMIi5=") {split($i,splitUMI,/=/); print splitUMI[2];exit } } }' | tr ';' '~')
  else
    UMIi5=0
  fi

  
  #create NGS filename in storage
  if [[ $NGS == Y ]]; then
    NGS=Y
    BAM=Y
    FILE=$(echo $LINE |
      awk -v LIB_STORAGE_FOLDER=$LIB_STORAGE_FOLDER -v BCi7="$BCi7" -v BCi5="$BCi5" -v sRBC="$sRBC" '{
        split(BCi7,splitBCi7,/~/)
        split(BCi5,splitBCi5,/~/)
        split(sRBC,splitsRBC,/~/)
        gsub("NGS","",$1)
        gsub("TGZ","",$1)
        nRUNs=split($1,splitNGS,/;/)
        for(i=3; i<=NF; i++) {
          if($i ~ "selRUNs=") {
            sub("selRUNs=","",$i)
            split($i,splitRUNs,/;/)
          }
        }
        for( n=1; n<=nRUNs; n++){
          if(n>1){
            newFILE=newFILE ":!:"LIB_STORAGE_FOLDER splitNGS[n]"~"splitBCi7[n]":"splitBCi5[n]"i"splitsRBC[n]"~"splitRUNs[n] ".bam"
          }else{
            newFILE=LIB_STORAGE_FOLDER splitNGS[n]"~"splitBCi7[n]":"splitBCi5[n]"i"splitsRBC[n]"~"splitRUNs[n] ".bam"
          }
        }
        print newFILE 
      }'
    )
  fi
  
  ${SCRIPT_DIR}main.sh "FILE=${FILE},NAME=${NAME},NGS=${NGS},IPstatus=${IPstatus},OXstatus=${OXstatus},TYPE=${TYPE},SLAM=${SLAM},FOLDER=${FOLDER},FOLDER_NAME=${FOLDER_NAME},rawTMP=${TMPdir},SINGULARITYdir=${SINGULARITYdir},downDIR=${downDIR},DEBUG=${DEBUG},VERSION=${VERSION},GENOME_VERSION=${GENOME_VERSION},BAM=${BAM},FORCE=${FORCE},SRR=${SRR},SRAlink=${SRAlink},fwADAPTOR=${fwADAPTOR},rvADAPTOR=${rvADAPTOR},N_TRIMM=${N_TRIMM},rawPAIRED=${rawPAIRED},onlyPAIRED=${onlyPAIRED},FASTQout=${FASTQout},FASTQoutRAW=${FASTQoutRAW},currSUB=${currSUB},MIN_LENGTH=${MIN_LENGTH},MAX_LENGTH=${MAX_LENGTH},RAW=${RAW},TRIMM=${TRIMM},FIRST=${FIRST},LAST=${LAST},INVERT=${INVERT},Ychrom=${Ychrom},RANDOMmulti=${RANDOMmulti},MM=${MM},FILTERING_INPUT=${FILTERING_INPUT},WIG=${WIG},WIG_FASTA=${WIG_FASTA},noNORM=${noNORM},only5end=${only5end},EXTEND=${EXTEND},COMPUTING=${COMPUTING},GRIDsystem=${GRIDsystem},SCRIPT_DIR=${SCRIPT_DIR},UTILITY_LOCATION=${UTILITY_LOCATION},UTILITY_DIR=${UTILITY_DIR},SYSTEM=${SYSTEM},nSPLITS=${nSPLITS},SE=${SE},BCi7=${BCi7},BCi5=${BCi5},sRBC=${sRBC},demuxFASTA=${demuxFASTA},extNORM=${extNORM},SE2nd=${SE2nd},DGE=${DGE},GEO=${GEO},noSTRANDED=${noSTRANDED},exportBAM=${exportBAM},exportBAMuncollapsed=${exportBAMuncollapsed},exportSalmon=${exportSalmon},extraSEQ=${extraSEQ},RATIOtracks=${RATIOtracks},maxCOUNT=${maxCOUNT},spikeINnorm=${spikeINnorm},PingPong=${PingPong},UMIi7=${UMIi7},UMIi5=${UMIi5},FORCEquant_unstranded=${FORCEquant_unstranded}"  &>>"${LOGs}/${NAME}.log.txt" &
done
wait

if [[ -f ${FOLDER}error_download.txt ]]; then
    for i in A B C x y z; do
      touch ${FOLDER}${i}XXXXXXXXXXXXX---pipeline_crashed_due to failed download and demux---XXXXXXXXXXXXX
    done
    exit
fi
  
if [[ -f ${FOLDER}error.txt ]]; then
  for i in A B C x y z; do
    touch ${FOLDER}${i}XXXXXXXXXXXXX---pipeline_crashed_ask_Dominik_for_help---XXXXXXXXXXXXX
  done
  exit
fi
if [[ -f ${FOLDER}error_network.txt ]]; then
  for i in A B C x y z; do
    touch ${FOLDER}${i}XXXXXXXXXXXXX---pipeline_might_have_crashed_partially_ask_Dominik_for_help---XXXXXXXXXXXXX
  done
  exit
fi

#exit if only raw-fasta generation requested
if [[ $RAW == Y ]]; then
  exit
fi

#-------------------------------------------------------------------------------------------------------
#perform DGE analysis if requested
if [[ $TYPE == RNAseq || $TYPE == RIPseq ]] || [[ $TYPE == CLIPseq ]] || [[ $TYPE == GROseq ]] || [[ $FORCEquant_unstranded == Y ]]; then
  mkdir ${FOLDER}gene-expression/ 

  if [[ $DGE == Y ]]; then
    if [[ ! -f ${FOLDER}DGE/DGE_info.txt ]]; then
      echo sample condition path | tr ' ' '\t' >${FOLDER}DGE/DGE_info.txt
      cat ${TMPdir}DGE_info.txt >>${FOLDER}DGE/DGE_info.txt
    else
      tail -n +2 ${FOLDER}DGE/DGE_info.txt >${TMPdir}DGE_info.txt
    fi
        
    cp ${UTILITY_LOCATION}transcript_to_gene_${VERSION}.txt ${TMPdir}transcript_to_gene.txt
    
    refGENO=$(head -n 1 ${TMPdir}DGE_info.txt | tr ' ' '\t' | cut -f 2)
    GENO=$(cat ${TMPdir}DGE_info.txt | tr ' ' '\t' | cut -f 2 | uniq | tail -n +2 | sort | sed '/^$/d' | tr "\n" " " | sed 's/[ \t]*$//')
    origGENO=$(echo $refGENO $GENO | tr ' ' '\t')
    GENO=$origGENO
    
    rm -rf ${FOLDER}DGE/*_vs_*
    
    for refGENO in $GENO; do
      GENO=$(echo $GENO | awk -v refGENO=$refGENO '{
        for(i=1; i<=NF; i++){
          if($i==refGENO) {$i=""}
        }
        print
      }')
      for currGENO in $GENO; do
        mkdir ${FOLDER}DGE/${refGENO}_vs_${currGENO}
      done
    done
    
    GENO=$(echo $origGENO | tr ' ' '~')
    echo $refGENO $GENO

    Nexec=1
  else
    #count normalization using edgeR only
    refGENO="none"
    GENO="none"   

    rm -rf ${TMPdir}DGE_info_all.txt
    echo sample condition path | tr ' ' '\t' >${TMPdir}DGE_info_all.txt
    cp ${UTILITY_LOCATION}transcript_to_gene_${VERSION}.txt ${TMPdir}transcript_to_gene.txt

    for i in $(seq 1 "${nFILES}"); do
      NAME=$(sed -n "${i}"p ${TMPdir}files.txt | awk '{ print $2 }')
      if [[ $TYPE == GROseq ]]; then
        Nexec=2
        echo $NAME $NAME ${TMPdir}${NAME}/salmon.quant.CDS/quant.sf >>${TMPdir}DGE_info_all.CDS.txt
        mkdir ${FOLDER}gene-expression-CDS/ 
      else
        Nexec=1
      fi
      echo $NAME $NAME ${TMPdir}${NAME}/salmon.quant/quant.sf >>${TMPdir}DGE_info_all.txt
    done
  fi
  
  
  cd $LOGs
  for i in $(seq 1 ${Nexec}); do
    COMMAND="${SCRIPT_DIR}DGE.sh"
    DGE_VARI="$VARI,refGENO=${refGENO},GENO=${GENO},Nexec=${Nexec}"

    if [[ $COMPUTING == C ]]; then
      if [[ $GRIDsystem == SLURM ]]; then sbatch --wait $COMMAND ${DGE_VARI}; else qsub -sync y $COMMAND ${VARI}; fi
    else
      $COMMAND ${DGE_VARI}
    fi
  done
fi  

#---------------------------------------------------------------------------------------------------------
#collect the results from all libraries and combine to single files
COMMAND=${SCRIPT_DIR}collect_numbers.sh

if [[ $COMPUTING == C ]]; then
  SLURMid=$(sbatch --parsable $COMMAND ${VARI})
else
  $COMMAND ${VARI}
fi

###################################################################################################
#create track hub containing tracks for all the libraries analyzed in this run

if [[ $GENOME_VERSION != ASM ]]; then
  ${SCRIPT_DIR}create_hub.sh $VARI
elif [[ $GENOME_VERSION == ASM ]]; then
  ${SCRIPT_DIR}create_hub.asm.sh $VARI
fi


###################################################################################################
#generate GEO output

if [[ $GEO == Y ]]; then
  
  FILES=$(cat "${TMPdir}files.txt" | cut -f 2 | tr '\n' '~')
  
  COMMAND=${SCRIPT_DIR}GEO.sh
  
  if [[ $COMPUTING == C ]]; then
    if [[ $GRIDsystem == SLURM ]]; then sbatch $COMMAND ${VARI},FILES=$FILES ; else { qsub -sync y $COMMAND ${VARI},FILES=$FILES &}; fi
  else
    $COMMAND ${VARI},FILES=$FILES
  fi
  
fi
###################################################################################################
#wait for collect numbers to finish

sbatch --dependency=afterany:${SLURMid} --wait --wrap="sleep 60s"
wait

###################################################################################################
#change permissins for all the files and FOLDERs to be accessible to everybody
chmod -R 777 $FOLDER
chmod -R 777 $TMPdir

rm -rf ~/.${FOLDER_NAME}.txt

if [[ $keepTMP == N ]]; then
  rm -rf $TMPdir
fi

PROCESSED_TIME=$(echo -e "$(date "+%s")" "$TIMEall" | awk '{ print ($1-$2)/60 }')
echo "total time =  ${PROCESSED_TIME}" min >>"${FOLDER}log.txt"
