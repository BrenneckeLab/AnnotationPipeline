#!/usr/bin/env bash
#$ -S /bin/bash

#$ -clear
#$ -e $JOB_NAME.e$JOB_ID.txt
#$ -o $JOB_NAME.o$JOB_ID.txt
#$ -cwd
#$ -q public.q
#$ -pe smp 12

#SBATCH --cpus-per-task=14
#SBATCH --mem=20g
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
locTMP=${TMPdir}TMP_map-genome/

mkdir -p $locTMP
chmod 777 $locTMP

export TMPDIR=$locTMP

####################################################################################################
##mapping to the rRNA_precursor - unmapped reads will be used for further steps
#reduce cures for bowtie by 2 to have free CPU power for pipe
redCORES=$(( CORES - 4 ))
redCORES=$(( CORES - 1 ))

gunzip --stdout ${TMPdir}${NAME}.fa.gz | \
  bowtie -f -v $MM -k 1 --un ${locTMP}${NAME}_rRNA_unmapped.fa --sam -p $redCORES -x ${UTILITY_LOCATION}rRNA - | \
  samtools view -bS - | \
  bedtools bamtobed -bed12  -i - | \
  awk -v FS="\t" -v TMP=$locTMP '
    BEGIN {OFS = FS} 
    { split($4,SPLIT_NAME,/:|=/)
      $4=$4":mapping=m"
      X[SPLIT_NAME[2]]=SPLIT_NAME[4]
      Y+=SPLIT_NAME[4]
      print 
    }
    END{
      for (i in X) {
        COUNT+=X[i]
      }
      print "ribosomal_mapping= " COUNT > TMP "stats1.tmp"
    }
    '  |
  mawk -v OFS="\t" '{split($4,SPLIT_NAME,/:|=/); $5=SPLIT_NAME[4]; print }' > ${locTMP}${NAME}_mapped.bed
      
  PROCESSED_TIME=$(echo -e $(date "+%s") "$TIMEx" | awk '{ print ($1-$2)/60 }')
  echo "bam to fasta - rRNA - processing_time=" "${PROCESSED_TIME}" >"${libFOLDER}time-log3.txt"
  TIMEx=$(date "+%s")

#---------------------------------------------------------------------------------------------------------
##mapping to the mitochondrial genome - unmapped reads will be used for further steps
#reduce cures for bowtie by 2 to have free CPU power for pipe
redCORES=$(( CORES - 3 ))

bowtie -f -v $MM -k 1 --un ${locTMP}${NAME}_mitochondrial_genome_unmapped.fa --sam -p $redCORES -x ${UTILITY_LOCATION}mitochondrion_DNA ${locTMP}${NAME}_rRNA_unmapped.fa | \
  samtools view -bS - | \
  bedtools bamtobed -bed12  -i - | \
  awk -v FS="\t" -v TMP="$locTMP" '
    BEGIN {OFS = FS} 
    { split($4,SPLIT_NAME,/:|=/)
      $4=$4":mapping=m"
      X[SPLIT_NAME[2]]=SPLIT_NAME[4]
      Y+=SPLIT_NAME[4]
      print 
    }
    END{
      for (i in X) {
        COUNT+=X[i]
      }
      print "mitochondrial_mapping= " COUNT > TMP "stats2.tmp"
    }
    ' |
  mawk -v OFS="\t" '{split($4,SPLIT_NAME,/:|=/); $5=SPLIT_NAME[4]; print }' >> ${locTMP}${NAME}_mapped.bed

rm -rf ${locTMP}${NAME}_rRNA_unmapped.fa

  PROCESSED_TIME=$(echo -e $(date "+%s") "$TIMEx" | awk '{ print ($1-$2)/60 }')
  echo "bam to fasta - mito - processing_time=" "${PROCESSED_TIME}" >>"${libFOLDER}time-log3.txt"
  TIMEx=$(date "+%s")


#---------------------------------------------------------------------------------------------------------
#mapping against the genome with up to # MM - unique mappers

if [[ $Ychrom == Y ]]; then INDEXext="_no-Mito_incl-Y"; else INDEXext="_no-Mito_excl-Y"; fi
if [[ $GENOME_VERSION == ASM ]]; then INDEXext="_$GENOME_VERSION"; fi

bowtie -f -v $MM -m 1 --best --strata --un ${locTMP}${NAME}_genome_unmapped.fa --sam -p $redCORES -x ${UTILITY_LOCATION}genome${INDEXext} ${locTMP}${NAME}_mitochondrial_genome_unmapped.fa | \
  samtools view -bS - | tee ${TMPdir}uniq.bam |\
  bedtools bamtobed -bed12  -i - | \
  awk -v FS="\t" -v TMP="$locTMP" '
    BEGIN {OFS = FS} 
    { split($4,SPLIT_NAME,/:|=/)
      $4=$4":mapping=u"
      COUNT+=SPLIT_NAME[4]
      print 
    }
    END{
      print "genome_unique_mapping_reads= " COUNT > TMP "stats3.tmp"
    }
    ' |
  mawk -v OFS="\t" '{split($4,SPLIT_NAME,/:|=/); $5=SPLIT_NAME[4]; print }' >> ${locTMP}${NAME}_mapped.bed

  PROCESSED_TIME=$(echo -e $(date "+%s") "$TIMEx" | awk '{ print ($1-$2)/60 }')
  echo "bam to fasta - uniq - processing_time=" "${PROCESSED_TIME}" >>"${libFOLDER}time-log3.txt"
  TIMEx=$(date "+%s")


#---------------------------------------------------------------------------------------------------------
#mapping against the genome with up to # MM - all mappers - only process unmapped reads from unique mapping

if [[ $RANDOMmulti == Y ]]; then
  Xvar="-M 1"
  UNCOLLAPSE() {
    seqkit fx2tab |
    awk -v OFS="\t" '{
      split($1,splitNAME,/:|=/)
      
      split($1,splitFORassembly,/count=/)
      
      $1=splitFORassembly[1] "count=" 1 

      for(i=1; i<=splitNAME[4];i++){ 
        print
      }
    }' | seqkit tab2fx 
  }
else
  Xvar="-a"
  UNCOLLAPSE() {
    cat
  }
fi
export -f UNCOLLAPSE

cat ${locTMP}${NAME}_genome_unmapped.fa | UNCOLLAPSE | 
bowtie -f -v $MM $Xvar --best --strata --un ${locTMP}${NAME}_unmapped.fa --sam -p $redCORES -x ${UTILITY_LOCATION}genome${INDEXext} - | \
  samtools view -bS - | tee ${TMPdir}multi.bam |\
  bedtools bamtobed -bed12  -i - | \
  awk -v FS="\t" -v TMP=$locTMP '
    BEGIN {OFS = FS} 
    { split($4,SPLIT_NAME,/:|=/)
      $4=$4":mapping=m"
      X[$4]=SPLIT_NAME[4]
      Y+=SPLIT_NAME[4]
      print 
    }
    END{
      for (i in X) {
        COUNT+=X[i]
        #print size and count for size_profile
        print length(i),X[i] > TMP "length_profile4.tmp"
      }
      print "genome_multi_mapping_reads= " COUNT > TMP "stats4.tmp"
      print "total_mappings_number_of_multi_mapping_reads= " Y > TMP "stats4.tmp"
      
    }
    ' |
  mawk -v OFS="\t" '{split($4,SPLIT_NAME,/:|=/); $5=SPLIT_NAME[4]; print }' >> ${locTMP}${NAME}_mapped.bed

  PROCESSED_TIME=$(echo -e $(date "+%s") "$TIMEx" | awk '{ print ($1-$2)/60 }')
  echo "bam to fasta - multi - processing_time=" "${PROCESSED_TIME}" >>"${libFOLDER}time-log3.txt"
  TIMEx=$(date "+%s")

#---------------------------------------------------------------------------------------------------------
#calculate normalization factor based on spike-ins
if [[ $spikeINnorm == Y ]]; then
 
  #build bowtie index
  bowtieBuild ${UTILITY_DIR}/spike-in_normalization/Mind.sRNA.fa ${locTMP}spike-in.index 

  #map unmapped reads to the spike-ins
  seqkit subseq -r 5:17 ${locTMP}${NAME}_unmapped.fa |
  bowtie -f --sam -p $redCORES -v 0 -m 1 -x ${locTMP}spike-in.index - |
    samtools view -bS - | tee ${TMPdir}multi.bam |\
    bedtools bamtobed -bed12  -i - | 
    LC_COLLATE=C sort -k1,1 -k2,2n |
    awk -v OFS="\t" '{split($4,splitNAME,/:|=/); X[$1]+=splitNAME[4]}END{for(i in X)   print i,X[i]}' | 
    sort -k2,2rn > ${libFOLDER}spike-in.raw.txt


  mawk -v OFS="\t" '
  {
    if($2>1000){
      split($1,splitNAME,/~/)
      X[splitNAME[1]]=$2/(splitNAME[2]*100)
    }
  }
  END{
    for(i in X){
      n+=1
      SUM+=X[i]
    }
    print SUM/n
  }' ${libFOLDER}spike-in.raw.txt > ${TMPdir}spike-in.norm

  NORM=$(cat ${TMPdir}spike-in.norm)

  mawk -v OFS="\t" '
  {
    print $1,$2/NORM
  }' ${libFOLDER}spike-in.raw.txt > ${libFOLDER}spike-in.norm.txt

fi

#---------------------------------------------------------------------------------------------------------
#export additional data if requested

#determine t-->C conversion in SLAMseq-reads
if [[ $SLAM == Y ]]; then
  rm -rf ${libFOLDER}TC_reads.txt 
  for MAPtype in uniq multi; do
    samtools view ${TMPdir}${MAPtype}.bam |
      mawk -v OFS="\t" -v OUTfile=${libFOLDER}TC_reads.txt -v OUTfile2=${TMPdir}SNPlist.txt  '
        BEGIN{
          #initiate all SNP types in array
          preBASES="A-C-G-T"
          split(preBASES,BASES,/-/)
          for(BASE in BASES){
            for (BASE2 in BASES){
              VAR=BASES[BASE]">"BASES[BASE2]
              X[VAR]=0
            }
          }

          #complementing pairs
          c["A"] = "T"; c["C"] = "G"; c["G"] = "C"; c["T"] = "A"; c["N"] = "N"
        }

        #define complementing function
          function comp(x, i, o) {
            o = ""
            for(i=1; i <= length;  i++){
              o = o c[substr(x, i, 1)]
            }
            return(o)
          }
      {
        #prepare the MD tag for looping
        split($13,MD,/:/)
        gsub(/[A-Z]/," & ",MD[3])
        n=split(MD[3],splitMD,/ /)

        #split the name TAG for count extraction
        split($1,splitNAME,/:|=/)

        #add current READ to the READ-array
        NAMEstore[$1]+=1

        #if read mapping contains SNP and read was not evaluated before enter SNP analysis
        if(n>1 &&  NAMEstore[$1]==1 ){
          #clear variables
          POS=0
          nMOD=0

          #loop through the MD tag
          for(i=2; i<=n; i=i+2){
            #get Reference base from the MD tag, calculate position in the seuqnce and get the Read base
            ORIG=splitMD[i]
            POS+=splitMD[i-1]+1
            REPL=substr($10,POS,1)

            #complement bases if read was mapped in antisense
            if($2==16){
              ORIG=comp(ORIG)
              REPL=comp(REPL)
            }

            #define SNP type
            VAR=ORIG">"REPL

            #count read number for SNP type
            X[VAR]+=splitNAME[4]
            #print SNP-type + count for later normalization
            print VAR,$1,length(splitNAME[2]),splitNAME[4] > OUTfile2

            #if T>C count number of T>C within the read
            if(ORIG=="T" && REPL=="C"){
              nMOD+=1
            }
          }

          #if T>C was present in the read - report it for later read filtering
          if(nMOD>0){
            print $1,nMOD >> OUTfile
          }
        }
      }
      END{
        #print all SNP-type numbers
        for(TYPE in X){
            print TYPE,X[TYPE]

        }
      }
      ' | sort -k1,1 > ${libFOLDER}SNPsummary.${MAPtype}.txt 
  done

  PROCESSED_TIME=$(echo -e $(date "+%s") "$TIMEx" | awk '{ print ($1-$2)/60 }')
  echo "bam to fasta - T>C - processing_time=" "${PROCESSED_TIME}" >>"${libFOLDER}time-log3.txt"
  TIMEx=$(date "+%s")
fi


#export bam files to result-directory
if [[ $exportBAM == Y ]]; then 
  cp ${TMPdir}uniq.bam ${libFOLDER}${NAME}.uniq.collapsed.bam
  cp ${TMPdir}multi.bam ${libFOLDER}${NAME}.multi.collapsed.bam
fi

if [[ $exportBAMuncollapsed == Y ]]; then 
  samtools view -h ${TMPdir}uniq.bam | 
    mawk -v OFS="\t" '{
      if($0~":count="){
        split($1, splitNAME, /:|=/)
        for(i=1;i<=splitNAME[4]; i++){ 
          print
        }
      }else{
        print
      }
    }' | 
  samtools view -bS > ${libFOLDER}${NAME}.uniq.bam
fi


  PROCESSED_TIME=$(echo -e $(date "+%s") "$TIMEx" | awk '{ print ($1-$2)/60 }')
  echo "bam to fasta - export bam - processing_time=" "${PROCESSED_TIME}" >>"${libFOLDER}time-log3.txt"
  TIMEx=$(date "+%s")


##for CHIPseq extend reads if required
if [[ $TYPE == CHIPseq ]]
then
  if [[ $EXTEND -gt 0 ]]
  then
    mv ${locTMP}${NAME}_mapped.bed ${locTMP}${NAME}_mapped.tmp
    awk -v EXTEND=$EXTEND -v OFS="\t" '
      BEGIN{
        while (getline < "'"${UTILITY_LOCATION}chrom.sizes"'")
        {
          CHRlength[$1]=$2
        }
        #close(CHRsizes)
      }
      {
        if($6 == "+") {
          $3=$2 + EXTEND
          if($3 >CHRlength["chr"$1]) {$3=CHRlength["chr"$1]}
        } else {
          if($6 == "-") {
            $2=$3 - EXTEND
            if($2<0) $2=0
          }
        }
        print $0
      
      }' ${locTMP}${NAME}_mapped.tmp > ${locTMP}${NAME}_mapped.bed
  fi
fi          

###sort reads according genome position and move to TMP folder

nLINES=$(wc -l "${locTMP}${NAME}_mapped.bed" | tr ' ' '\t' | cut -f 1)
mkdir ${TMPdir}split_mapping/

  PROCESSED_TIME=$(echo -e $(date "+%s") "$TIMEx" | awk '{ print ($1-$2)/60 }')
  echo "bam to fasta - counmt lines - processing_time=" "${PROCESSED_TIME}" >>"${libFOLDER}time-log3.txt"
  TIMEx=$(date "+%s")


LC_ALL=C sort -S15G --parallel=$CORES -k4,4 -k1,1 -k2,2n ${locTMP}${NAME}_mapped.bed | \
  awk -v OFS="\t" -v nLINES=$nLINES -v TMP=${TMPdir}split_mapping/${NAME}_mapped -v nSPLITS=${nSPLITS} '
  BEGIN{
    
    CHUNKsize=nSPLITS
    currREADinCHUNK=1
    currCHUNK=1
    currFILE=TMP currCHUNK ".bed"
  }
  {
    nTOTAL+=1
    if(currREADinCHUNK < CHUNKsize) {
      print > currFILE
      currREADinCHUNK+=1
    }else{
      if(currREADinCHUNK == CHUNKsize){
        print > currFILE
        lastID=$4
        currREADinCHUNK+=1
      }else{
        if(lastID!=$4 && nLINES-nTOTAL > 250000 ){
          currCHUNK+=1
          currREADinCHUNK=1
          currFILE=TMP currCHUNK  ".bed"
          print > currFILE
        }else{
          print > currFILE
          currREADinCHUNK+=1
        }
      }
    }
  }' 

  PROCESSED_TIME=$(echo -e $(date "+%s") "$TIMEx" | awk '{ print ($1-$2)/60 }')
  echo "bam to fasta - sort and process - processing_time=" "${PROCESSED_TIME}" >>"${libFOLDER}time-log3.txt"
  TIMEx=$(date "+%s")


#---------------------------------------------------------------------------------------------------------
cat ${locTMP}stats*.tmp > ${locTMP}stats.txt

#---------------------------------------------------------------------------------------------------------
#move unmapped reads to tmp-folder, add n tag and count on the fly
awk -v TMP=$locTMP '{ if($0~">") {
    print $0":mapping=n:ann=unmapped:final_ann=unmapped:filtered=unmapped:mapcount=0:fine_ann=unmapped"
    split($1,SPLIT_NAME,/:|=/)
    X+=SPLIT_NAME[4]
  } else print }
END{print "unmapped_reads= " X >> TMP "stats.txt" }' ${locTMP}${NAME}_unmapped.fa |\
  gzip > ${TMPdir}${NAME}_unmapped.fa.gz

    PROCESSED_TIME=$(echo -e $(date "+%s") "$TIMEx" | awk '{ print ($1-$2)/60 }')
  echo "bam to fasta - move files - processing_time=" "${PROCESSED_TIME}" >>"${libFOLDER}time-log3.txt"
  TIMEx=$(date "+%s")


#---------------------------------------------------------------------------------------------------------
#report percent statistics
TOTAL_READS=$( grep -F NR_of_reads_after_filtering_ ${libFOLDER}log.txt | awk '{ print $NF}' )
TOTAL_MAPPED_COUNTS=$( awk '{ if($0 !~ "unmapped|total") X+=$2} END{print X}' < ${locTMP}stats.txt  ) 

rRNA_MAPPED=$( grep -F ribosomal_mapping ${locTMP}stats.txt | awk -v TOTAL_READS=$TOTAL_READS '{ print $NF*100/TOTAL_READS}' )
rRNA_MAPPED_COUNTS=$( grep -F ribosomal_mapping "${locTMP}stats.txt" | awk -v TOTAL_READS=$TOTAL_READS '{ print $NF}' )
MITOCHONDRIAL_MAPPED=$( grep -F mitochondrial_mapping ${locTMP}stats.txt | awk -v TOTAL_READS=$TOTAL_READS '{ print $NF*100/TOTAL_READS}' )
MITOCHONDRIAL_MAPPED_COUNTS=$( grep -F mitochondrial_mapping "${locTMP}stats.txt" | awk -v TOTAL_READS=$TOTAL_READS '{ print $NF}' )
GENOME_UNIQ=$( grep -F genome_unique_mapping_reads ${locTMP}stats.txt | awk -v TOTAL_READS=$TOTAL_READS '{ print $NF*100/TOTAL_READS}' )
GENOME_UNIQ_COUNTS=$( grep -F genome_unique_mapping_reads ${locTMP}stats.txt | awk -v TOTAL_READS=$TOTAL_READS '{ print $NF}' )
GENOME_MULTI=$( grep -F genome_multi_mapping_reads ${locTMP}stats.txt | awk -v TOTAL_READS=$TOTAL_READS '{ print $NF*100/TOTAL_READS}' )
GENOME_MULTI_COUNTS=$( grep -F genome_multi_mapping_reads ${locTMP}stats.txt | awk -v TOTAL_READS=$TOTAL_READS '{ print $NF}' )
UNMAPPED=$( grep -F unmapped_reads ${locTMP}stats.txt | awk -v TOTAL_READS=$TOTAL_READS '{ print $NF*100/TOTAL_READS}' )
UNMAPPED_COUNTS=$( grep -F unmapped_reads ${locTMP}stats.txt | awk -v TOTAL_READS=$TOTAL_READS '{ print $NF}' )
TOTAL_MAPPED=$( bc <<< "$rRNA_MAPPED"+"$MITOCHONDRIAL_MAPPED"+"$GENOME_UNIQ"+"$GENOME_MULTI" )
TOTAL=$( bc <<< "$rRNA_MAPPED"+"$MITOCHONDRIAL_MAPPED"+"$GENOME_UNIQ"+"$GENOME_MULTI"+"$UNMAPPED" )

{
  printf "\nnumber_of_mapped_reads= $TOTAL_MAPPED_COUNTS \n" 
  printf "\n%% of reads mapping and read-counts - numbers represent real sequencing-counts, not read identities: \n" 
  printf "reads mapped to rRNA-precursor= ${rRNA_MAPPED}%%\t${rRNA_MAPPED_COUNTS}\n" 
  printf "reads mapped to mitochondrial-genome= ${MITOCHONDRIAL_MAPPED}%%\t${MITOCHONDRIAL_MAPPED_COUNTS}\n" 
  printf "reads mapped uniquely to the genome mapping= ${GENOME_UNIQ}%%\t${GENOME_UNIQ_COUNTS}\n" 
  printf "reads mapped multiple times to the genome= ${GENOME_MULTI}%%\t${GENOME_MULTI_COUNTS}\n"
  printf "reads_not_mapped= ${UNMAPPED}%%\t${UNMAPPED_COUNTS}\n" 
  printf "reads mapped in total= ${TOTAL_MAPPED}%% \n" 
  printf "total percent= ${TOTAL}%% \n" 
} >> ${libFOLDER}log.txt
#---------------------------------------------------------------------------------------------------------
#extraction of genomic unique mappers and calculation of normalization factor
  
NORMALIZATION=$( echo $GENOME_UNIQ_COUNTS | awk '{ print $1/10000000 }' )
echo ${NORMALIZATION} - normalized to 10000000 uniquely mapping reads  >> ${libFOLDER}normalization.txt

###################################################################################################
#clean up
if [[ $DEBUG != Y ]]; then
  rm -rf $locTMP
fi

PROCESSED_TIME=$(echo -e "$(date "+%s")" "$TIME" | awk '{ print ($1-$2)/60 }' )
printf "map reads - processing_time=	${PROCESSED_TIME} \n">> ${libFOLDER}time-log.txt


exit

