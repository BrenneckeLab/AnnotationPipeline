#!/usr/bin/env bash
#$ -S /bin/bash

#$ -clear
#$ -e $JOB_NAME.e$JOB_ID.txt
#$ -o $JOB_NAME.o$JOB_ID.txt
#$ -cwd
#$ -q public.q
#$ -pe smp 10

#SBATCH --cpus-per-task=6
#SBATCH --mem=20g
#SBATCH -e "%x.e.%j.txt"
#SBATCH -o "%x.o.%j.txt"
#XX-SBATCH --qos=short-XX
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
locTMP="${TMPdir}TMP_plot_TEhist/"

mkdir -p $locTMP
chmod 777 $locTMP

export TMPDIR=$locTMP
###################################################################################################

redCORES=$(($CORES - 2))

if [[ $Ychrom == Y ]]; then INDEXext="incl-Y"; else INDEXext="excl-Y"; fi

#! add switch to allow changing the MM setting away from the default 2

if [[ $TYPE == sRNAseq ]] || [[ $TYPE == sRNAseqIP ]] || [[ $TYPE == CHIPseq ]] || [[ $TYPE == DNAseq ]] || [[ $TYPE == RIPseq ]]; then
  if [[ $TYPE == sRNAseq ]] || [[ $TYPE == sRNAseqIP ]]; then
    if [[ $GENOME_VERSION == Cel ]]; then    
      sizeCUTOFF=20
    else
      sizeCUTOFF=23
    fi
  else
    sizeCUTOFF=18
  fi

  if [[ ! -z $extraSEQ ]]; then
    TEindex=${rawTMP}TE_index
  else
    TEindex=${UTILITY_LOCATION}bowtie-TE_${INDEXext}
  fi

  gunzip -c ${TMPdir}${NAME}_annotated.fa.gz |
    mawk -v sizeCUTOFF=$sizeCUTOFF '{
      if($1~">" && $1!~"final_ann=rRNA|final_ann=tRNA|final_ann=snoRNA|final_ann=snRNA"){
        #?@ #removed genome-filtering to allow non-genomic TEs to be analyzed
        #?@ if($1~"mapping=u" || $1~"mapping=m") {
          split($1,splitNAME,/:|=/)
          #filter for piRNAs by length
          if(length(splitNAME[2])>sizeCUTOFF || length(splitNAME[2]) == 21 ) {
            print $1"\n"splitNAME[2]
          }
        #?@ }
      }
    }' |
    bowtie --best --strata -f -v 2 -a --sam -p $redCORES -x ${TEindex} - |
    awk -v OFS="\t" -v TMP=$locTMP '
    BEGIN{
      srand(rseed)
    }
    {
      print $0 > TMP "all.sam"
      if($1 ~ "^@") {
        print $0
      }else{
        if($0~"_randomDistribution"){
          SATELLITE[$1][$3][$4]=$0
        }else{
          if($0~ "XM:i:1") {
            print $0
            print $0
          }else{
            if($0~ "XM:i:2"){
              ARRAY[$1][$3][$4]=$0
            }
          }
        }
      }
    }
    END{
      for(READ in ARRAY) {
        if(length(ARRAY[READ]) >0 ) {
          for(TE in ARRAY[READ]){
            for(POS in ARRAY[READ][TE]) {
              print ARRAY[READ][TE][POS]
            }
          }
        }
      }
      for(READ in SATELLITE) {
        if(length(SATELLITE[READ]) == 1 ) {
          for(TE in SATELLITE[READ]){
            n=0
            for(POS in SATELLITE[READ][TE]) {
                n+=1
                nSATELLITE[n]=POS
                
                print SATELLITE[READ][TE][POS]
                
            }
            split(READ,splitREAD,/:|=/)
            for( i=1; i<=splitREAD[4];i++){
              x=int(rand()*n+1)
              LOOKUP="count="splitREAD[4]
              gsub(LOOKUP,"count=1",SATELLITE[READ][TE][nSATELLITE[x]])
              print SATELLITE[READ][TE][nSATELLITE[x]]
            }
          }
        }
      }
    }' |
    samtools view -bS - |
    bedtools bamtobed -i - | tee ${locTMP}TEmappings_man.fullLength.bed | 
    awk -v OFS="\t" -v only5end=${only5end} '
    {
      if ( only5end == "Y" ){
        if( $6 == "+" ) $3=$2+1
        if( $6 == "-" ) $2=$3-1
      }
      print
    }' >${locTMP}TEmappings_man.bed

    samtools view -bS ${locTMP}all.sam |
    bedtools bamtobed -i - | tee ${locTMP}TEmappings_man.fullLength.bed | 
    awk -v OFS="\t" -v only5end=${only5end} '
    {
      if ( only5end == "Y" ){
        if( $6 == "+" ) $3=$2+1
        if( $6 == "-" ) $2=$3-1
      }
      print
    }' >${locTMP}TEmappings_man.all.bed

else
  gunzip -c ${TMPdir}${NAME}_annotated.fa.gz |
    mawk '{
      if($1~">" && $1!~"final_ann=rRNA|final_ann=tRNA|final_ann=snoRNA|final_ann=snRNA"){
        #?@ #removed genome-filtering to allow non-genomic TEs to be analyzed
        #?@ if($1~"mapping=u" || $1~"mapping=m") {
          split($1,splitNAME,/:|=/)
          #filter for piRNAs by length
          print $1"\n"splitNAME[2]
        #?@ }
      }
    }' >${locTMP}reads.fa

  STAR --runThreadN $redCORES \
    --genomeDir ${UTILITY_LOCATION}star_TE_${INDEXext}/ \
    --readFilesIn ${locTMP}reads.fa \
    --outFileNamePrefix ${locTMP} --outFilterMatchNmin $MIN_LENGTH \
    --outSAMmode NoQS --readFilesCommand cat --alignEndsType EndToEnd --twopassMode Basic \
    --outReadsUnmapped Fastx --outMultimapperOrder Random --outSAMtype SAM --outFilterMultimapNmax 1000 --winAnchorMultimapNmax 2000 \
    --outFilterMismatchNmax ${MM} --seedSearchStartLmax 30 \
    --outFilterType Normal --alignSJoverhangMin 10 --alignSJDBoverhangMin 1 --outStd SAM |
    awk -v OFS="\t" -v TMP=$locTMP '
      BEGIN{
      srand(rseed)
    }
    {
      print $0 > TMP "STAR_all.sam"
      if($1 ~ "^@") {
        print $0
      }else{
        if($0~"_randomDistribution"){
          SATELLITE[$1][$3][$4]=$0
        }else{
          if($1 ~ "^@") {
            print $0
          }else{
            if($0~ "NH:i:1") {
              print $0
              print $0
            }else{
              if($0~ "NH:i:2"){
                ARRAY[$1][$3][$4]=$0
              }
            }
          }
        }
      }
    }
    END{
      for(READ in ARRAY) {
        if(length(ARRAY[READ]) == 1 ) {
          for(TE in ARRAY[READ]){
            for(POS in ARRAY[READ][TE]) {
              print ARRAY[READ][TE][POS]
            }
          }
        }
      }
      for(READ in SATELLITE) {
        if(length(SATELLITE[READ]) == 1 ) {
          for(TE in SATELLITE[READ]){
            n=0
            for(POS in SATELLITE[READ][TE]) {
                n+=1
                nSATELLITE[n]=POS
                
                print SATELLITE[READ][TE][POS]
                
            }
            split(READ,splitREAD,/:|=/)
            for( i=1; i<=splitREAD[4];i++){
              x=int(rand()*n+1)
              LOOKUP="count="splitREAD[4]
              gsub(LOOKUP,"count=1",SATELLITE[READ][TE][nSATELLITE[x]])
              print SATELLITE[READ][TE][nSATELLITE[x]]
            }
          }
        }
      }
    }' |
    tee ${locTMP}STAR_TE_processed.sam |
    samtools view -bS - |
    bedtools bamtobed -split -i - |
    awk -v OFS="\t" -v only5end=${only5end} '
    {
      if ( only5end == "Y" ){
        if( $6 == "+" ) $3=$2+1
        if( $6 == "-" ) $2=$3-1 
      }
      print
    }' >${locTMP}TEmappings_man.bed
fi

if [[ $SLAM == Y ]]; then
  N=2
  grep TCslam ${locTMP}TEmappings_man.bed > ${locTMP}TEmappings_man.TC.bed
else
  N=1
  EXT=""
fi

#---------------------------------------------------------------------------------------------------------
if [[ ! -z $extraSEQ ]]; then
  TElengthfile=${rawTMP}TE.chrom.sizes
else
  TElengthfile=${UTILITY_LOCATION}TE_${INDEXext}.sizes
fi

if [[ ${TYPE} == sRNAseq || ${TYPE} == sRNAseqIP ]]; then
  STRINGENCYvec="normal all"
else
  STRINGENCYvec="normal"
fi

  
for i in $(seq 1 $N); do
  if [[ $i -eq 1 ]]; then
    EXT=""
  elif [[ $i -eq 2 ]]; then
    EXT=".TC"
  fi
  for STRINGENCY in $STRINGENCYvec ; do
    if [[ $STRINGENCY == all ]]; then
      STRINGENCYext=".all"
    else
      STRINGENCYext=""
    fi

    #sort and uncollapse bed-file
    sort --parallel=$CORES -k1,1 -k2,2n ${locTMP}TEmappings_man${EXT}${STRINGENCYext}.bed |
      mawk -v OFS="\t" -v TMP=$locTMP -v TYPE=$TYPE -v EXT=${EXT}${STRINGENCYext} '
      {
        split($4,splitNAME,/:|=/)

        for(i=1; i<=splitNAME[4]; i++) {  
          if((TYPE=="sRNAseq" || TYPE=="sRNAseqIP") && length(splitNAME[2])==21) {
            print  > TMP "TEmappings_man_sort" EXT ".21mer.bed"
          }else{
            print 
          }
        }
      }' >${locTMP}TEmappings_man_sort${EXT}${STRINGENCYext}.bed

    if [[ $TYPE == CHIPseq || $TYPE == DNAseq ]]; then
      bedtools genomecov -d -g ${TElengthfile} -i ${locTMP}TEmappings_man_sort${EXT}${STRINGENCYext}.bed |
        mawk -v OFS="\t" -v NORMfactor=${NORMfactor} ' 
        {
          $3=$3/NORMfactor
          print 
        }' >${libFOLDER}${NAME}_TE-all_man${EXT}${STRINGENCYext}.bg

    else
      #generate bedgraphs
      bedtools genomecov -d -strand + -g ${TElengthfile} -i ${locTMP}TEmappings_man_sort${EXT}${STRINGENCYext}.bed |
        mawk -v OFS="\t" -v NORMfactor=${NORMfactor} ' 
        {
          $3=$3/NORMfactor
          print 
        }' >${libFOLDER}${NAME}_TE-sense_man${EXT}${STRINGENCYext}.bg
      bedtools genomecov -d -strand - -g ${TElengthfile} -i ${locTMP}TEmappings_man_sort${EXT}${STRINGENCYext}.bed |
        mawk -v OFS="\t" -v NORMfactor=${NORMfactor} ' 
        {
          $3=-$3/NORMfactor
          print 
        }' >${libFOLDER}${NAME}_TE-antisense_man${EXT}${STRINGENCYext}.bg

      #generate bedgraphs for 21mers
      if [[ $TYPE == sRNAseq || $TYPE == sRNAseqIP ]]; then
        bedtools genomecov -d -strand + -g ${TElengthfile} -i ${locTMP}TEmappings_man_sort${EXT}${STRINGENCYext}.21mer.bed |
          mawk -v OFS="\t" -v NORMfactor=${NORMfactor} ' 
          {
            $3=$3/NORMfactor
            print 
          }' >${libFOLDER}${NAME}_TE-sense_man${EXT}.21mer${STRINGENCYext}.bg
        bedtools genomecov -d -strand - -g ${TElengthfile} -i ${locTMP}TEmappings_man_sort${EXT}${STRINGENCYext}.21mer.bed |
          mawk -v OFS="\t" -v NORMfactor=${NORMfactor} ' 
          {
            $3=-$3/NORMfactor
            print 
          }' >${libFOLDER}${NAME}_TE-antisense_man${EXT}.21mer${STRINGENCYext}.bg
      fi

      #generate simple-style bowtie mapped TE-count table
      if [[ $TYPE == sRNAseq || $TYPE == sRNAseqIP ]] || [[ $TYPE == CLIPseq ]]; then

        #generate count table
        awk -v OFS="\t" -v NORMfactor=$NORMfactor -v PATH=${libFOLDER} -v TEfile=${TElengthfile} -v EXT=${EXT}${STRINGENCYext} '
        BEGIN{  
          #generate commands for direct output
          FILENAME=PATH "TE_RPKM_sense-bowtie"EXT".txt"
          OUT_sense_RPKM="sort -k1,1 > "FILENAME
          FILENAME=PATH "TE_RPKM_antisense-bowtie"EXT".txt"
          OUT_antisense_RPKM="sort -k1,1 > "FILENAME
          FILENAME=PATH "TE_counts_sense-bowtie"EXT".txt"
          OUT_sense_COUNT="sort -k1,1 > "FILENAME
          FILENAME=PATH "TE_counts_antisense-bowtie"EXT".txt"
          OUT_antisense_COUNT="sort -k1,1 > "FILENAME

          #read in all TE identifiers and add a 0 count
          while((getline TEs < TEfile) > 0) {
            split(TEs, splitTEs, /\t/)
            COUNTte[splitTEs[1]]["+"]=0
            COUNTte[splitTEs[1]]["-"]=0
            LENGTHte[splitTEs[1]]=splitTEs[2]
          }
        }
        {
          #get counts per TE into array -- use read counts
          split($4, splitNAME, /:|=/)
          COUNTte[$1][$6]+=splitNAME[4]
        }
        END{
          #print out final table
          for(TE in COUNTte){
            print TE,COUNTte[TE]["+"]/NORMfactor | OUT_sense_COUNT
            print TE"_AS", COUNTte[TE]["-"]/NORMfactor | OUT_antisense_COUNT
            print TE,COUNTte[TE]["+"]/NORMfactor*1000/LENGTHte[TE] | OUT_sense_RPKM
            print TE"_AS", COUNTte[TE]["-"]/NORMfactor*1000/LENGTHte[TE] | OUT_antisense_RPKM
          }
        }' ${locTMP}TEmappings_man_sort${EXT}${STRINGENCYext}.bed

      fi
    fi


    #store bed-file permanently on open directory
    sort --parallel=$CORES -k1,1 -k2,2n ${locTMP}TEmappings_man_sort${EXT}${STRINGENCYext}.bed | gzip > ${libFOLDER}TEmappings${EXT}${STRINGENCYext}.bed.gz
    #store siRNA bed-file permanently on open directory
    if [[ $TYPE == sRNAseq || $TYPE == sRNAseqIP ]]; then
      sort --parallel=$CORES -k1,1 -k2,2n ${locTMP}TEmappings_man_sort${EXT}${STRINGENCYext}.21mer.bed | gzip > ${libFOLDER}TEmappings${EXT}${STRINGENCYext}.21mer.bed.gz
    fi
  done


  #!remove EXT
  EXT=""

  if [[ ! -z $extraSEQ ]]; then
    TEfile=${rawTMP}TE_noSplice.fa
  else
    TEfile=${UTILITY_LOCATION}TE_annot_new_WO_DUST.fa
  fi

  if [[ $PingPong == Y ]] && [[ $TYPE == sRNAseq ]]; then
      #generate pingpong table

      awk -v OFS="\t" -v NORMfactor=$NORMfactor -v PATH=${locTMP} -v TEfile=${TEfile} -v EXT= '
      BEGIN{  
        #generate commands for direct output
        FILENAME=PATH "TE_pingpong"EXT".txt"
        print "TE","POS","5END_plus","3END_plus","5END_minus","3END_minus" > FILENAME

        #read in all TE identifiers and add a 1 pseudocount at each position of the TE for each strand and each end
        while((getline LINE < TEfile) > 0) {
          if(LINE~">"){
            split(LINE, splitTEs, /\t|>/)
            TEs=splitTEs[2]
          }else{
            for(i=0; i<=length(LINE); i++) {
              COUNTte[TEs]["+"]["5END"][i]=1
              COUNTte[TEs]["+"]["3END"][i]=1
              COUNTte[TEs]["-"]["5END"][i]=1
              COUNTte[TEs]["-"]["3END"][i]=1
            }
          }
        }
      }
      {
        #get counts per TE/strand/end/position into array -- use read counts
        split($4, splitNAME, /:|=/)
        if(length(splitNAME[2])>23 ){
          if($6=="+") {
            COUNTte[$1]["+"]["5END"][$2]+=splitNAME[4]
            COUNTte[$1]["+"]["3END"][$3-1]+=splitNAME[4]
          }else{
            COUNTte[$1]["-"]["5END"][$3-1]+=splitNAME[4]
            COUNTte[$1]["-"]["3END"][$2]+=splitNAME[4]
          }
        }
      }
      END{
        #print out final table
        for(TE in COUNTte){
          for(POS in COUNTte[TE]["+"]["5END"]){
            print TE,POS,COUNTte[TE]["+"]["5END"][POS]-1,COUNTte[TE]["+"]["3END"][POS]-1,COUNTte[TE]["-"]["5END"][POS]-1,COUNTte[TE]["-"]["3END"][POS]-1
          }
        }
      }' ${locTMP}TEmappings_man.fullLength.bed | sort -k1,1 -k2,2n >>${locTMP}TE_pingpong.txt

      mkdir -p ${libFOLDER}ping-pong-analysis
      rm -rf ${libFOLDER}ping-pong-analysis/*
      mkdir -p ${libFOLDER}phasing-analysis
      rm -rf ${libFOLDER}phasing-analysis/*
      Rscript ${SCRIPT_DIR}ping-pong.R baseDIR=${libFOLDER}/ INFILE=${locTMP}TE_pingpong.txt NORMfactor=${NORMfactor}
    fi  

done

###################################################################################################
if [[ $DEBUG != Y ]]; then
  rm -rf $locTMP
fi

PROCESSED_TIME=$(echo -e "$(date "+%s")" "$TIME" | mawk '{ print ($1-$2)/60 }')
printf "generate TE bedgrpahs  - processing_time= ${PROCESSED_TIME} \n" >>"${libFOLDER}time-log.txt"

exit

exit
