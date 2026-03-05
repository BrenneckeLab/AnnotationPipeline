#!/usr/bin/env bash
#$ -S /bin/sh

#$ -clear
#$ -e $JOB_NAME.e$JOB_ID.txt
#$ -o $JOB_NAME.o$JOB_ID.txt
#$ -cwd
#$ -q public.q
#$ -pe smp 5

#SBATCH --cpus-per-task=5
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
if [[ $GRIDsystem == SLURM ]]; then 
  CORES=$(( SLURM_CPUS_PER_TASK * @HYPER@ ))
  MEM=$(scontrol show job $SLURM_JOBID | grep TRES | awk '{ split($NF, X, /,|=|G/);{print X[5]-5}}' | head -n 1)
elif [[ -z ${NSLOTS+x} ]];then 
  CORES=4
  MEM=40
else 
  CORES=$NSLOTS
fi
echo $CORES

#presetting of environment to preload modules and to prepare all files required
#load all functions
source ${SCRIPT_DIR}tools

###################################################################################################
#create local tmp
locTMP="${TMPdir}TMP_prepare-reads/"

mkdir -p "$locTMP"
chmod 777 "$locTMP"

export TMPDIR=$locTMP

###################################################################################################

#create empty files to have an empty file if there would not be a single read filtered
touch ${locTMP}filtered.txt
touch ${locTMP}length_filtered.txt
touch ${locTMP}N_filtered.txt

#---------------------------------------------------------------------------------------------------------
#evaluate if library is fasta or fastq

FASTQ_EVALUATION=$( seqkit head -n 1 $FILE | sed -n 1p)
echo $FASTQ_EVALUATION asd
if [[ $FASTQ_EVALUATION == '@'* ]]; then
  echo library is fastq >>${libFOLDER}log.txt
elif [[ $FASTQ_EVALUATION == '>'* ]]; then
  echo library is fasta >>${libFOLDER}log.txt
else
  echo file $FILE unknown filetype >>${libFOLDER}log.txt
  exit
fi

#---------------------------------------------------------------------------------------------------------
printf "\n" >>${libFOLDER}log.txt

REV="cat"
if [[ $TYPE == RNAseq || $TYPE == RNAseq || $TYPE == RNAseq ]]; then
  if [[ $INVERT != Y ]]; then
    REV="seqkit seq --reverse --complement"
  fi
else
  if [[ $INVERT == Y ]]; then
    REV="seqkit seq --reverse --complement"
  fi
fi


###################################################################################################

#---------------------------------------------------------------------------------------------------------
#run read preperation
if [[ $MIN_LENGTH -lt 18 ]]; then
  WINDOW=$MIN_LENGTH
else
  WINDOW=18
fi

redCORES=$(( CORES - 1 ))
seqkit seq --line-width 0 --only-id $FILE |
seqkit fq2fa | 
  $REV |
  seqkit fx2tab | 
  mawk -v TRIMM=$TRIMM -v FIRST=$FIRST -v LAST="$LAST" '
  {
      if (TRIMM=="Y") {
        $NF=substr($NF,FIRST,LAST-FIRST)
        if($NF=="") $NF="S"
      }
      print
  }' |  
  #remove spaces/tabs and other stuff from header and length-filter reads
  mawk -v MIN=$MIN_LENGTH -v MAX=$MAX_LENGTH -v TMP=$locTMP '
  BEGIN{
    ADAPTOR_DIMER=0
    LENGTH_FILTERED_SHORT=0
    LENGTH_FILTERED_LONG=0
  }
  {
      if($2 == "X"){
        ADAPTOR_DIMER+=1
      }else{
        if (length($2) < MIN  ) {
          LENGTH_FILTERED_SHORT+=1
          print > TMP "filtered.short.txt"
       }else{
          if ( length($2) > MAX){
            LENGTH_FILTERED_LONG+=1
            print > TMP "filtered.long.txt"
          }else{
            print
          }
        }
      }
  }
  END{
    print "adaptor_dimer= "ADAPTOR_DIMER  > TMP "stats.txt"
    print "length_filtered_short= "LENGTH_FILTERED_SHORT  > TMP "stats.txt"
    print "length_filtered_long= "LENGTH_FILTERED_LONG  > TMP "stats.txt"
  }' |
  mawk -v TMP=$locTMP -v TYPE=$TYPE -v INVERT=${INVERT} '
  BEGIN{
    N_FILTERED=0
    N=1
  }
  { 
      #trimm reads if requested
            #filter N containing reads
            if($NF ~ "N" ) {
              N_FILTERED+=1
            } else {
              print ">"$1"\n"$NF
            }
  }
  END{
    print "N_filtered= " N_FILTERED > TMP "stats2.txt" 
  }' |
  bbduk in=stdin.fa out=stdout.fa threads=2 entropy=0.35 entropywindow=$WINDOW entropyk=4 -Xmx4g 2>${locTMP}filtered.txt | 
  seqkit fx2tab | 
  #add UMIs before collapsing
  mawk -v OFS="\t" '
  {
    if($1~"UMI"){
      split($1,splitNAME,/::/)
      print $1,$NF,splitNAME[2]
    }else{
      print
    }
  }' | 
  LC_ALL=C sort --parallel=$redCORES -S${MEM}G -k2,2 -k3,3 |  
  mawk -v OFS="\t" '{
    if($3==""){
      print
    }else{
      if(SEQ == $2 && UMI == $3){
        DUPLI+=1
      }else{
        if(NR>1){
          #print
          print NAME":!:"DUPLI, SEQ, UMI
        }
        NAME=$1
        SEQ=$2
        UMI=$3
        DUPLI=1
      }
    }
  }END{
    print NAME":!:"DUPLI, SEQ, UMI

  }' | 
  awk -v TMP=$locTMP -v maxCOUNT=$maxCOUNT '
  {
      TOTALcount+=1
      
      if(NR==1){
        SEQ=$2
        COUNT=1
        N=1
        if($1~"::") {
          split($1,splitNAME,/::/)
          TAGcount=1
          DUPLI=splitNAME[2]
        }else {
          TAGcount=1
        }
      }else{
        if(SEQ == $2){
          COUNT+=1
          if(NF==3){
            split($1,splitTAG,/:!:/)
            UMIcount+=splitTAG[3]
          }
          #@ if($1~"::") {
          #@   split($1,splitNAME,/::/)
          #@   if(DUPLI!~splitNAME[2]) {
          #@     TAGcount+=1
          #@     DUPLI=DUPLI"-"splitNAME[2]
          #@   }else {
          #@     TAGcount+=1
          #@   }
          #@ }   
        }else{
          if(COUNT > maxCOUNT){
            COUNT=maxCOUNT
          }

          if(NR>1){
            print ">NR_" N "_" COUNT/TAGcount "_UMIcount_" UMIcount ":" SEQ ":count=" COUNT "\n" SEQ
          }
          N=N+1
          SEQ=$2
          COUNT=1
          if(NF==3){
            split($1,splitTAG,/:!:/)
            UMIcount=splitTAG[3]
          }
          #@ if($1~"::") {
          #@   split($1,splitNAME,/::/)
          #@   TAGcount=1
          #@   DUPLI=splitNAME[2]
          #@ }else {
          #@   TAGcount=1
          #@ }
        }
      }
  }
  #report filtered read counts 
  END { 
    if(COUNT > maxCOUNT){
      COUNT=maxCOUNT
    }
    print ">NR_" N "_" COUNT/TAGcount UMItag "_" UMIcount ":" SEQ ":count=" COUNT "\n" SEQ
    N=N+1

    print "NR_of_reads_after_filtering_(total_reads_for_%_calculations)= "TOTALcount > TMP "stats3.txt"
    print "collapsed_reads= "N "\n" > TMP "stats3.txt"
    print "duplication_rate= " TOTALcount/N > TMP "stats3.txt"
  } ' | gzip -2 >${TMPdir}${NAME}.fa.gz

cat ${locTMP}stats2.txt >>${locTMP}stats.txt
cat ${locTMP}stats3.txt >>${locTMP}stats.txt

#export artifact filtered read number to log
ARTIFACT_FILTERED=$(grep -F "Total Removed" ${locTMP}filtered.txt | cut -f 2 | tr ' ' '\t' | cut -f 1)
cp ${locTMP}filtered.txt ${libFOLDER}LOGs/BBDUK-stats.txt
echo artifact_filtered= $ARTIFACT_FILTERED >>${libFOLDER}log.txt
cat ${locTMP}stats.txt >>${libFOLDER}log.txt

#calculate and report normalization factor to 10M reads in fasta and initiate normalization.txt file
filteredREADS=$(grep -F NR_of_reads_after_filtering_ ${locTMP}stats.txt | awk '{ print $NF/10000000}')

printf "#last line in this file represents the normalization factor used for the normalization
##of the tracks generated by the annotation pipeline
###divide by the factor annotated to normalize reads\n 
NORMfactor - description\n" >${libFOLDER}normalization.txt
printf "${filteredREADS} - normalized to 10000000 reads in fasta after filtering\n" >>${libFOLDER}normalization.txt
#---------------------------------------------------------------------------------------------------------

if [[ $DEBUG != Y ]]; then
  rm -rf $locTMP
fi

PROCESSED_TIME=$(echo -e $(date "+%s") $TIME | awk '{ print ($1-$2)/60 }')
echo "prepare_reads - processing_time=" ${PROCESSED_TIME} >>${libFOLDER}time-log.txt

exit
