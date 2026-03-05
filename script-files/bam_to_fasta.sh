#!/usr/bin/env bash
#$ -S /bin/bash

#$ -clear
#$ -e $JOB_NAME.e$JOB_ID.txt
#$ -o $JOB_NAME.o$JOB_ID.txt
#$ -cwd
#$ -q public.q
#$ -pe smp 6

#SBATCH --cpus-per-task=6
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

#define proper CORES
if [[ $GRIDsystem == SLURM ]]; then CORES=$(( SLURM_CPUS_PER_TASK * @HYPER@ )); elif [[ -z ${NSLOTS+x} ]]; then CORES=4; else CORES=$NSLOTS; fi
echo $CORES

#presetting of environment to preload modules and to prepare all files required
#load all functions
source ${SCRIPT_DIR}tools

locTMP=${TMPdir}TMP_bam_to_fasta/

mkdir -p $locTMP
chmod 777 $locTMP
  
export TMPDIR=$locTMP

###################################################################################################
#create local tmp

#redirect tempdir
export TMPDIR=${TMPdir}/${RANDOM}/
export TMP=${TMPdir}/${RANDOM}/
mkdir -p $TMPDIR

printf "library is bam and gets converted to zipped fasta
adaptor-clipping statistics can be found in the adapor-clipping-report.X.log file\n\n" >>${libFOLDER}log.txt

####################################################################################################

#calculate the minimum length required after adaptor clipping - to not run into empty lines in 4N trimming
minLENGTH=$(echo $((N_TRIMM * 2)))

rm -rf ${locTMP}${NAME}_*.tab

####################################################################################################
#pre-processing of SRR
if [[ $SRR == Y ]]; then
  #set nMERGE to 1
  nMERGE=1

  #downlading and processing of SRA files
  printf "downloading of remote file: ${FILE}\n" >>${libFOLDER}log.txt

  # set the maximum number of allowed repeats
  MAX_REPEATS=5
  # initialize the counter variable
  REPEATS=1
  #initialize TEST variable
  TEST="" 
  #initialize BCtype
  BCtype=""
  
  #download file from SRA
  LINKs=$(echo $SRAlink | tr '=~' '\t\t' | cut -f 2-)
  echo $LINKs

  while [[ -z $TEST ]]; do
    # check if the maximum number of allowed repeats has been reached
    if [[ $REPEATS -gt $MAX_REPEATS ]]; then
      printf "SRA download failed to produce output after $MAX_REPEATS attempts for library ${NAME}\n" >>${FOLDER}error.txt
      exit
    fi

    parallel --tmpdir $TMPDIR -j 2 ". ${SCRIPT_DIR}functions; sraFunct {1} $locTMP ${NAME} {#} " ::: $LINKs

    if [[ -f ${locTMP}${NAME}.error ]]; then
      #remove error file
      rm -rf ${locTMP}${NAME}.error
      # increment the counter variable
      REPEATS=$((REPEATS + 1))
      randTIME=$( echo $((  RANDOM % 3000 + 1200 )) )
      sleep  ${randTIME}
    else
      TEST="OK"
    fi
  done
  
  #generate variables
  splitFILES=$(echo $NAME | sed 's/:!:/\t/g')
  splitSRBC=$(echo $sRBC | sed 's/~/\t/g')
  nMERGE=1
else
  #determine number of libraries to process
  splitFILES=$(echo $FILE | sed 's/:!:/\t/g')
  splitSRBC=$(echo $sRBC | sed 's/~/\t/g')
  nMERGE=$(echo $splitFILES | wc -w )
  echo $nMERGE $sRBC

  if [[ $nMERGE -gt 1 ]]; then
    printf "merging of files:\n${FILE}\n" | sed 's/:!:/\n/g' >>"${libFOLDER}log.txt"
  fi
fi

PROCESSED_TIME=$(echo -e $(date "+%s") "$TIMEx" | awk '{ print ($1-$2)/60 }')
echo "bam to fasta - pre-step - processing_time=" "${PROCESSED_TIME}" >"${libFOLDER}time-log2.txt"
TIMEx=$(date "+%s")

####################################################################################################
#loop through requested libraries
origNAME=$NAME
fwADAPTORorig=$fwADAPTOR

for i in $(seq 1 $nMERGE); do
  FILE=$(echo $splitFILES | awk -v i=$i '{print $i}')
  sRBC=$(echo $splitSRBC | awk -v i=$i '{print $i}')
  NAME=${origNAME}~${i}
  
  if [[ $SRR != Y  ]]; then
    #determine correct opening procedure for input-file
    samtools view $FILE 2>${locTMP}test.txt | head >/dev/null
    cat ${locTMP}test.txt
    TEST=$(grep main_samview ${locTMP}test.txt )
    
    if [[ $TEST != *"main_samview"* ]]; then
      echo library $FILE is bam >>"${libFOLDER}log.txt"
      FORMAT="bam"
      READ_COMMAND="samtools view $FILE"
    else    
      FASTQ_EVALUATION=$( seqkit seq --line-width 0 ${FILE} | head -10 | sed -n 1p)
      if [[ "$FASTQ_EVALUATION" == '@'* ]]; then
        FORMAT="fastq"
        echo library $FILE is fastq >>"${libFOLDER}log.txt"
        READ_COMMAND="seqkit fx2tab --only-id $FILE "
        elif [[ "$FASTQ_EVALUATION" == '>'* ]]; then
        FORMAT="fasta"
        echo library $FILE is fasta >>"${libFOLDER}log.txt"
        READ_COMMAND=" seqkit seq --line-width 0 $FILE | fastaToFastq in1=stdin.fa qfake=40 fastawrap=10000000 out1=stdout.fq | seqkit fx2tab "
      else
        echo file "$FILE" unknown filetype >>"${libFOLDER}log.txt"
        exit
      fi
    fi


    #---------------------------------------------------------------------------------------------------------
    #sRBC demultiplexing on fasta-files
    if [[ $demuxFASTA == Y && ! -z $sRBC && ($FORMAT == "fasta" || $FORMAT == "fastq" ) ]]; then
      #assemble sRBD adaptor squene
      eval $READ_COMMAND  |
        mawk -v OFS="\t" -v sRBC=$sRBC '
        BEGIN{
          sRBC=sRBC "AGATCGGAA"

        }
        {
          if($2 ~ sRBC){
            print 
          }
        }' | 
        seqkit tab2fx > ${locTMP}demux.fq

      FILE=${locTMP}demux.fq
      FORMAT="fastq"
      READ_COMMAND="seqkit fx2tab --only-id $FILE"      
    fi


    ####################################################################################################
    #extraction of reads from bam
    ##forward reads (1st mates) are getting adaptor-clipped and trimmed
    ##reverse reads (2nd mates) are stored in a file on TMP location
    
    $READ_COMMAND |  
    mawk -v LOG=${libFOLDER}flat-error.log -v TMP=$locTMP -v NAME=$NAME -v FORMAT=$FORMAT -v GEO=$GEO -v UMIi7=$UMIi7 -v UMIi5=$UMIi5 '
      {
        if (FORMAT=="bam") {
          if( $2 == 77 || $2 == 4) {
            #evaluate UMI
            if (UMIi7 > 0) {
              if(colBC==""){
                for(i=10; i<=NF; i++){
                  if($i ~ "BC:Z:"){
                    colBC=i
                    break
                  }
                }
              }
              split($colBC,UMI,/BC:Z:/)
              seqUMIi7=substr(UMI[2],9,UMIi7)
              seqUMIi7="_UMIi7="seqUMIi7
            }
            if (UMIi5 > 0) {
              if(colBC==""){
                for(i=10; i<=NF; i++){
                  if($i ~ "BC:Z:"){
                    colBC=i
                    break
                  }
                }
              }
              split(colBC,UMI,"BC:Z:")
              seqUMIi5=substr(9,UMIi5,UMI[2])
              seqUMIi5="_UMIi5="seqUMIi5
            }

            print NR"_"$10 seqUMIi7 seqUMIi5 "\t" $10 "\t" $11
            print $1 seqUMIi7 seqUMIi5 "/1\t"$10"\t"$11 > TMP NAME "_paired.tab"
            if(GEO=="Y")print $0 > TMP NAME "_1.sam"
          } else {
            if ( $2 == 141 ) {
              print NR"_"$10 seqUMIi7 seqUMIi5 "\t" $10 "\t" $11 > TMP NAME "_2nd.tab"
              print $1 seqUMIi7 seqUMIi5 "/2\t"$10"\t"$11 > TMP NAME "_paired.tab"
              if(GEO=="Y")print $0 > TMP NAME "_2.sam"
            } else {
              print "unexpected sam-flag detected" >> LOG
              exit 1
            }
          }
        } else {
          print
        }
    }' | tr ' ' '\t' >${locTMP}${NAME}_1st.tab
    
  fi
  PROCESSED_TIME=$(echo -e $(date "+%s") "$TIMEx" | awk '{ print ($1-$2)/60 }')
  echo "bam to fasta - first -step - processing_time=" "${PROCESSED_TIME}" >>"${libFOLDER}time-log2.txt"
  TIMEx=$(date "+%s")

  ####################################################################################################
  #export raw fastq before adaptor trimming
  if [[ $FASTQoutRAW == Y ]]; then
    if [[ ! -s ${locTMP}${NAME}_2nd.tab ]]; then
      seqkit tab2fx ${locTMP}${NAME}_1st.tab | gzip > ${libFOLDER}${NAME}_raw.fq.gz
    else
      seqkit tab2fx ${locTMP}${NAME}_paired.tab | gzip > ${libFOLDER}${NAME}_raw.fq.gz
    fi
  fi

  ####################################################################################################
  #subsample input files

  if [[ ! -z $currSUB ]]; then
    RAND=$RANDOM
    awk -v OFS="\t" -v SEED=$RAND '
      BEGIN{
        srand(SEED)
      }
      {
        print rand(),$0
    }' ${locTMP}${NAME}_1st.tab >${locTMP}${NAME}_1st.tmp
    LC_COLLATE=C sort --parallel=$CORES -k1,1n ${locTMP}${NAME}_1st.tmp | head -n $currSUB | cut -f 2- >${locTMP}${NAME}_1st.tab
    
    if [[ -f ${locTMP}${NAME}_2nd.tab ]]; then
      awk -v OFS="\t" -v SEED=$RAND '
      BEGIN{
        srand(SEED)
      }
      {
        print rand(),$0
      }' ${locTMP}${NAME}_2nd.tab >${locTMP}${NAME}_2nd.tmp
      LC_COLLATE=C sort --parallel=$CORES -k1,1n ${locTMP}${NAME}_2nd.tmp | head -n $currSUB | cut -f 2- >${locTMP}${NAME}_2nd.tab
    fi
  fi

  ####################################################################################################
  #process the 1st read except it is not requested
  if [[ $SE2nd != Y && $onlyPAIRED != Y ]]; then
    
    printf "############################################################################\ntrimming of 1st mate\n\n" >${libFOLDER}adaptor-clipping-report_1st.${i}.log

     if [[ $BCtype == *"Stark_3prime_BC_PROseq" ]]; then
       fwADAPTOR="TGGAATTCTCGGGTGCCAAGG"
       printf "due to presence of Stark 3prime barcoded adaptors, adaptor for clipping was changed to $fwADAPTOR \n\n" >>"${libFOLDER}log.txt"
     elif [[ ! -z $sRBC ]]; then
      fwADAPTOR="NN${sRBC}${fwADAPTORorig}"
      printf "due to presence of sRBC, adaptor for clipping was changed to $fwADAPTOR \n\n" >>"${libFOLDER}log.txt"
    fi
    
    #adaptors are trimmed for the 1st mate and n nucleotides get trimmed
    redCORES=$(( $CORES / 3 ))
    #@ parallel --pipepart --tmpdir $TMPDIR --line-buffer --rrs --round-robin -j $redCORES -a ${locTMP}${NAME}_1st.tab "seqkit tab2fx --threads 1 | cutadapt -a ${fwADAPTOR} --max-n=200 - 2>> ${libFOLDER}adaptor-clipping-report_1st.log " |
    seqkit tab2fx --threads 1 ${locTMP}${NAME}_1st.tab | 
      CUTADAPT -a ${fwADAPTOR} --cores=3 --max-n=200 - 2>> ${libFOLDER}adaptor-clipping-report_1st.${i}.log  |
      seqkit fx2tab |  
      awk -v N_TRIMM=${N_TRIMM} -v TMP=$TMPdir -v minLENGTH=$minLENGTH '
        {
            if(length($2)>minLENGTH) {
              xSTART=N_TRIMM+1
              xEND=length($2)-2*N_TRIMM
              #print N-trimmed reads

              if($1 ~ "UMI" ){
                n=split($1,splitHEAD,"_")
                if( n > 2 ){
                  UMI1=":!:FW_"splitHEAD[3]
                }
                if( n > 3 ){
                  UMI2=":!:FW_"splitHEAD[4]
                }
              }
              print NR"::"substr($2,1,N_TRIMM)substr($2,N_TRIMM+xEND+1,N_TRIMM)"X" UMI1 UMI2 "\t" substr($2,xSTART,xEND) "\t" substr($3,xSTART,xEND)
              
            }else {
              print  NR"::X\tX\t!"
            }
        }
      ' |  seqkit tab2fx > ${TMPdir}${NAME}_raw.fq
  fi

  PROCESSED_TIME=$(echo -e $(date "+%s") "$TIMEx" | awk '{ print ($1-$2)/60 }')
  echo "bam to fasta - second -step - processing_time=" "${PROCESSED_TIME}" >>"${libFOLDER}time-log2.txt"
  TIMEx=$(date "+%s")

  ####################################################################################################
  # if input file was paired-end process the 2nd mate and if required the raw paired file

  if [[ $SE == Y ]]; then rm -rf ${locTMP}${NAME}_2nd.tab; fi

  if [[ -s ${locTMP}${NAME}_2nd.tab ]] && [[ $TYPE != CapSeq ]] && [[ $TYPE != CHIPseq ]] && [[ $TYPE != sRNAseq ]]; then
    #---------------------------------------------------------------------------------------------------------
    #processing of 2nd mate - gets clipped and then reverse complemented and added to the 1st mate file
    
    if [[ $onlyPAIRED != Y ]]; then
      printf "\n############################################################################\ntrimming of 2nd mate\n\n" >>${libFOLDER}adaptor-clipping-report_2nd.${i}.log
      
      redCORES=$(( $CORES / 3 ))
      #@ parallel --tmpdir $TMPDIR --pipepart --line-buffer --rrs --round-robin -j $redCORES -a ${locTMP}${NAME}_2nd.tab " seqkit tab2fx --threads 1 | cutadapt -a ${rvADAPTOR} --max-n=200 - 2>> ${libFOLDER}adaptor-clipping-report_2nd.log" |
      seqkit tab2fx --threads 1 ${locTMP}${NAME}_2nd.tab | 
        CUTADAPT -a ${rvADAPTOR} --cores=3 --max-n=200 - 2>> ${libFOLDER}adaptor-clipping-report_2nd.${i}.log |
        seqkit seq --reverse --complement | seqkit fx2tab | 
        awk -v N_TRIMM=${N_TRIMM} -v TMP=$TMPdir -v minLENGTH=$minLENGTH '
            {
                if(length($2)>minLENGTH) {
                  xSTART=N_TRIMM+1
                  xEND=length($2)-2*N_TRIMM
                  #print N-trimmed reads

                  if($1 ~ "UMI" ){
                    n=split($1,splitHEAD,"_")
                    if( n > 2 ){
                      UMI1=":!:RV_"splitHEAD[3]
                    }
                    if( n > 3 ){
                      UMI2=":!:RV_"splitHEAD[4]
                    }
                  }

                  print NR"::"substr($2,1,N_TRIMM)substr($2,N_TRIMM+xEND+1,N_TRIMM)"X" UMI1 UMI2 "\t"substr($2,xSTART,xEND)"\t"substr($3,xSTART,xEND)

                }else{
                  print  NR"::X\tX\t!"
                }
            }
        ' | 
        #! add reverse complementin using seqkit
        seqkit tab2fx >>${TMPdir}${NAME}_raw.fq
    fi
    
    if [[ $rawPAIRED == Y ]]; then
      #---------------------------------------------------------------------------------------------------------
      #processing of the raw paired file - interleaved paired file gets adaptor-clipped and 4N trimmed
      
      printf "\n############################################################################\ntrimming of read-pairs\n\n" >>${libFOLDER}adaptor-cd clipping-report.${i}.log
      ###replace minLENGTH by the minimum length defined by AnnotationPipeline
      seqkit tab2fx ${locTMP}${NAME}_paired.tab | 
        CUTADAPT -a ${fwADAPTOR} -A ${rvADAPTOR} --cores=3 --max-n=200 -O 4 --interleaved - 2>>${libFOLDER}adaptor-clipping-report.${i}.log | seqkit fx2tab |
        awk -v OFS="\t" -v N_TRIMM=${N_TRIMM} '
          BEGIN{
            xSTART=N_TRIMM+1
          }
          {
              xEND=length($NF)-2*N_TRIMM
              print $1,substr($2,xSTART,xEND),substr($3,xSTART,xEND)

          }
        ' | seqkit tab2fx | gzip >>${libFOLDER}${origNAME}_raw-paired.fq.gz
    fi
    
    elif [[ -s ${locTMP}${NAME}_2nd.fa ]] && [[ $TYPE == CapSeq ]]; then
    printf "due to TYPE = CapSeq 2nd mate removed from analysis\n\n" >>"${libFOLDER}log.txt"
  fi
done

#---------------------------------------------------------------------------------------------------------
#reset name
NAME=$origNAME

#---------------------------------------------------------------------------------------------------------
#merge files
if [[ $nMERGE -gt 1 ]]; then
  cat ${TMPdir}${NAME}~*_raw.fq >${TMPdir}${NAME}_raw.fq
else
  mv ${TMPdir}${NAME}~1_raw.fq ${TMPdir}${NAME}_raw.fq
fi

#---------------------------------------------------------------------------------------------------------
#transfer fastq if requested
if [[ $FASTQout == Y ]]; then
  gzip -c ${TMPdir}${NAME}_raw.fq > ${libFOLDER}${NAME}.pre-processed.fq.gz
fi

#---------------------------------------------------------------------------------------------------------
#prepare sequence file for GEO submission
if [[ $GEO == Y ]]; then
  if [[ -s ${locTMP}${NAME}_2nd.tab ]]; then
    #split bam file into read-mates
    samtools view -bS ${locTMP}${NAME}*_1.sam >${TMPdir}${NAME}_demuxed_1.bam
    SUM=$(md5sum ${TMPdir}${NAME}_demuxed_1.bam)
    echo ${NAME}_raw_1 $SUM | tr ' ' '\t' >${TMPdir}raw_1.md5
    
    samtools view -bS ${locTMP}${NAME}*_2.sam >${TMPdir}${NAME}_demuxed_2.bam
    SUM=$(md5sum ${TMPdir}${NAME}_demuxed_2.bam)
    echo ${NAME}_raw_2 $SUM | tr ' ' '\t' >${TMPdir}raw_2.md5
  else
    #copy full sequencing file and calculate md5sum
    cp $FILE ${TMPdir}${NAME}_demuxed.bam
    SUM=$(md5sum ${TMPdir}${NAME}_demuxed.bam)
    echo ${NAME}_raw $SUM | tr ' ' '\t' >${TMPdir}raw.md5
  fi
fi
#---------------------------------------------------------------------------------------------------------
if [[ $DEBUG != Y ]]; then
  rm -rf "$locTMP"
fi

PROCESSED_TIME=$(echo -e $(date "+%s") "$TIME" | awk '{ print ($1-$2)/60 }')
echo "bam to fasta - processing_time=" "${PROCESSED_TIME}" >>"${libFOLDER}time-log.txt"
echo $PROCESSED_TIME

exit
