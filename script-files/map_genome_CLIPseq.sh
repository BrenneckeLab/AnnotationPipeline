#!/usr/bin/env bash
#$ -S /bin/bash

#$ -clear
#$ -e $JOB_NAME.e$JOB_ID.txt
#$ -o $JOB_NAME.o$JOB_ID.txt
#$ -cwd
#$ -q public.q
#$ -pe smp 12

#SBATCH --cpus-per-task=14
#SBATCH --mem=40g
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

gunzip --stdout ${TMPdir}${NAME}.fa.gz | \
  bowtie -f -v $MM -k 1 --un ${locTMP}${NAME}_rRNA_unmapped.fa --sam -p $redCORES -x ${UTILITY_LOCATION}rRNA - | \
  samtools view -bS - | \
  bedtools bamtobed -bed12 -i - | \
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
        #print size and count for size_profile
        print length(i),X[i] > TMP "length_profile.txt"
      }
      print "ribosomal_mapping= " COUNT > TMP "stats.txt"
    }
    ' |
  mawk -v OFS="\t" '{split($4,SPLIT_NAME,/:|=/); $5=SPLIT_NAME[4]; print }' > ${locTMP}${NAME}_mapped.bed
      
#---------------------------------------------------------------------------------------------------------
##mapping to the mitochondrial genome - unmapped reads will be used for further steps
#reduce cures for bowtie by 2 to have free CPU power for pipe
redCORES=$(( CORES - 3 ))

bowtie -f -v $MM -k 1 --un ${locTMP}${NAME}_mitochondrial_genome_unmapped.fa --sam -p $redCORES -x ${UTILITY_LOCATION}mitochondrion_DNA ${locTMP}${NAME}_rRNA_unmapped.fa | \
  samtools view -bS - | \
  bedtools bamtobed -bed12 -i - | \
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
        #print size and count for size_profile
        print length(i),X[i] >> TMP "length_profile.txt"
      }
      print "mitochondrial_mapping= " COUNT >> TMP "stats.txt"
    }
    ' |
  mawk -v OFS="\t" '{split($4,SPLIT_NAME,/:|=/); $5=SPLIT_NAME[4]; print }' >> ${locTMP}${NAME}_mapped.bed

rm -rf ${locTMP}${NAME}_rRNA_unmapped.fa

#---------------------------------------------------------------------------------------------------------
#mapping against the genome with up to # MM - unique mappers

if [[ $Ychrom == Y ]]; then INDEXext="incl-Y"; else INDEXext="excl-Y"; fi
if [[ $GENOME_VERSION == ASM ]]; then INDEXext=""; fi

rm -rf ${locTMP}star/
mkdir ${locTMP}star/
redCORES=$((CORES - 1))

STAR --runThreadN $redCORES \
  --genomeDir ${UTILITY_LOCATION}star_${VERSION}/${INDEXext}/ \
  --readFilesIn ${locTMP}${NAME}_mitochondrial_genome_unmapped.fa \
  --outFileNamePrefix ${locTMP}star/ --outFilterMatchNmin $MIN_LENGTH \
  --outSAMmode NoQS --readFilesCommand cat --alignEndsType EndToEnd --twopassMode Basic \
  --outReadsUnmapped Fastx --outMultimapperOrder Random --outSAMtype SAM \
  --outFilterMultimapNmax 1000 --winAnchorMultimapNmax 2000 --outSAMattributes NH HI NM MD AS nM jM jI XS \
  --outFilterMismatchNoverLmax 0.10 --seedSearchStartLmax 30 --limitOutSAMoneReadBytes 400000 \
  --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --alignIntronMax 1 \
  --outFilterType BySJout --alignSJoverhangMin 15 --alignSJDBoverhangMin 1 --outStd SAM |
  tee ${locTMP}star/mapped.sam >${locTMP}tmp.sam

  awk -v OFS="\t" -v TMP=$locTMP '
  {
    if($1~"@") {
      print
      print > TMP "uniq.sam"
      print > TMP "multi.sam"
    } else {
      split($1,SPLIT_NAME,/:|=/) 
      if($5==255){
        $1=$1":mapping=u"
        COUNTuniq+=SPLIT_NAME[4]
        #print size and count for size_profile
        print length(SPLIT_NAME[2]),SPLIT_NAME[4] >> TMP "length_profile.txt"
        print > TMP "uniq.sam"
      } else {
        $1=$1":mapping=m"
        X[SPLIT_NAME[2]]=SPLIT_NAME[4]
        Y+=SPLIT_NAME[4]
      }
      print
      print > TMP "multi.sam"
    }
  }
  END{
    print "genome_unique_mapping_reads= " COUNTuniq >> TMP "stats.txt"
    for (i in X) {
      COUNT+=X[i]
      #print size and count for size_profile
      print length(i),X[i] >> TMP "length_profile.txt"
    }
    print "genome_multi_mapping_reads= " COUNT >> TMP "stats.txt"
    print "total_mappings_number_of_multi_mapping_reads= " Y >> TMP "stats.txt"
  }' ${locTMP}tmp.sam  |
  awk -v OFS="\t" -v TMP=$locTMP '
  BEGIN{
    SWITCH="N"
    TC="N"
  }
  {
    if($1~"@") {
      a=b
      print
    } else {
      if($6 !~ /D|I|X|H/ ){
        matchedBASES=$6
        sub("M","",matchedBASES)
        if($0!~"MD:Z:"matchedBASES){
          #prepare stuff to identify T->C conversions
          MD=$17
          sub("MD:Z:","",MD)
          gsub(/[A-Z]/,"~&~",MD)
          nMM=split(MD,splitMD,/~/)
          for(i=2; i<nMM; i=i+2) {
            if(splitMD[i]=="T"){ 
              READnt=substr($10,splitMD[i-1]+1,1)
              if(READnt=="C"){
                TC="Y"
              }
            }
          }

        }
      }else{
        if($17~"\^T"){
          Tdel="Y"
        }
      }
      if(TC=="Y" || Tdel=="Y"){
        print $1,TC,Tdel > TMP "CLIPtags.txt"
      }
      print

      TC="N"
      Tdel="N"
    }    
  }' |
  samtools view -bS - |
  bedtools bamtobed -bed12 -i -  |
  mawk -v OFS="\t" '{split($4,SPLIT_NAME,/:|=/); $5=SPLIT_NAME[4]; print }' >>${locTMP}${NAME}_mapped.bed

cp -r ${locTMP}star/ ${TMPdir}
cp -r ${locTMP}CLIPtags.txt ${TMPdir}

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
      ' | sort -k1,1 >> ${libFOLDER}SNPsummary.txt 
  done

  PROCESSED_TIME=$(echo -e $(date "+%s") "$TIMEx" | awk '{ print ($1-$2)/60 }')
  echo "bam to fasta - T>C - processing_time=" "${PROCESSED_TIME}" >>"${libFOLDER}time-log3.txt"
  TIMEx=$(date "+%s")
fi


#export bam files to result-directory
if [[ $exportBAM == Y ]]; then 
  #sort sam and store as bam file
  MEM=$(scontrol show job $SLURM_JOBID | grep TRES | awk '{ split($NF, X, /,|=|G/);{print X[5]-5}}' | head -n 1)
  xMEM=$(( $MEM / $CORES * 1000 ))
  samtools sort --output-fmt=BAM -m ${xMEM}M --threads=$CORES -o ${libFOLDER}${NAME}.uniq.collapsed.bam ${locTMP}uniq.sam 
  samtools sort --output-fmt=BAM -m ${xMEM}M --threads=$CORES -o ${libFOLDER}${NAME}.multi.collapsed.bam ${locTMP}multi.sam 
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



#---------------------------------------------------------------------------------------------------------
###sort reads according genome position and move to TMP folder

nLINES=$(wc -l ${locTMP}${NAME}_mapped.bed | tr ' ' '\t' | cut -f 1)
mkdir ${TMPdir}split_mapping/

LC_ALL=C sort -S30G --parallel=$CORES -k4,4 ${locTMP}${NAME}_mapped.bed | \
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

#---------------------------------------------------------------------------------------------------------
#create star high confidence splice-junction bed
  awk -v OFS="\t" -v TMPdir=$TMPdir '{
    if($4==1) STRAND="+"; else if($4==2) STRAND="-"; else STRAND="."
    if($5==0) MOTIF="non-canonical"; else if($5==1) MOTIF="GT/AG"; else if($5==2) MOTIF="CT/AC"; else if($5==3) MOTIF="GC/AG"; else if($5==4) MOTIF="CT/GC"; else if($5==5) MOTIF="AT/AC"; else if($5==6) MOTIF="GT/AT"
    if($6==1) ANNOTATED="annotated"; else if ($6==0) ANNOTATED="unannotated"
    if($7>$8) {
      print "chr"$1,$2-5,$3+5,ANNOTATED"_"MOTIF"-junktion_unique-reads="$7"_multi-mappers="$8"_max-spliced-overhang="$9,".",STRAND,$2-5,$3+5,"255,0,0",2,"5,5",0","$3-$2+5 > TMPdir "splice-junctions_uniq.bed12"
    }
    print "chr"$1,$2-5,$3+5,ANNOTATED"_"MOTIF"-junktion_unique-reads="$7"_multi-mappers="$8"_max-spliced-overhang="$9,".",STRAND,$2-5,$3+5,"255,0,0",2,"5,5",0","$3-$2+5 > TMPdir "splice-junctions_all.bed12"
  } ' ${locTMP}star/SJ.out.tab

#---------------------------------------------------------------------------------------------------------
#move unmapped reads to tmp-folder, add n tag and count on the fly
awk -v TMP=$locTMP '{ if($0~">") {
    print $0":mapping=n:ann=unmapped:final_ann=unmapped:filtered=unmapped:mapcount=0:fine_ann=unmapped"
    split($1,SPLIT_NAME,/:|=/)
    X+=SPLIT_NAME[4]
  } else print }
END{print "unmapped_reads= " X >> TMP "stats.txt" }' ${locTMP}star/Unmapped.out.mate1 |\
  gzip > ${TMPdir}${NAME}_unmapped.fa.gz
  
#---------------------------------------------------------------------------------------------------------
#report percent statistics
TOTAL_READS=$( grep -F NR_of_reads_after_filtering_ ${libFOLDER}log.txt | awk '{ print $NF}' )
TOTAL_MAPPED_COUNTS=$( awk '{ if($0 !~ "unmapped|total") X+=$2} END{print X}' < ${locTMP}stats.txt  ) 

rRNA_MAPPED=$( grep -F ribosomal_mapping ${locTMP}stats.txt | awk -v TOTAL_READS=$TOTAL_READS '{ print $NF*100/TOTAL_READS}' )
rRNA_MAPPED_COUNTS=$( grep -F ribosomal_mapping ${locTMP}stats.txt | awk -v TOTAL_READS=$TOTAL_READS '{ print $NF}' )
MITOCHONDRIAL_MAPPED=$( grep -F mitochondrial_mapping ${locTMP}stats.txt | awk -v TOTAL_READS=$TOTAL_READS '{ print $NF*100/TOTAL_READS}' )
MITOCHONDRIAL_MAPPED_COUNTS=$( grep -F mitochondrial_mapping ${locTMP}stats.txt | awk -v TOTAL_READS=$TOTAL_READS '{ print $NF}' )
GENOME_UNIQ=$( grep -F genome_unique_mapping_reads ${locTMP}stats.txt | awk -v TOTAL_READS=$TOTAL_READS '{ print $NF*100/TOTAL_READS}' )
GENOME_UNIQ_COUNTS=$( grep -F genome_unique_mapping_reads ${locTMP}stats.txt | awk -v TOTAL_READS=$TOTAL_READS '{ print $NF}' )
GENOME_MULTI=$( grep -F genome_multi_mapping_reads ${locTMP}stats.txt | awk -v TOTAL_READS=$TOTAL_READS '{ print $NF*100/TOTAL_READS}' )
GENOME_MULTI_COUNTS=$( grep -F genome_multi_mapping_reads ${locTMP}stats.txt | awk -v TOTAL_READS=$TOTAL_READS '{ print $NF}' )
UNMAPPED=$( grep -F unmapped_reads ${locTMP}stats.txt | awk -v TOTAL_READS=$TOTAL_READS '{ print $NF*100/TOTAL_READS}' )
UNMAPPED_COUNTS=$( grep -F unmapped_reads "${locTMP}stats.txt" | awk -v TOTAL_READS=$TOTAL_READS '{ print $NF}' )
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
#measure occurences of potential wrong splicings

printf "

counts of sequences that map multiple times but mappings have at least start or end coordinate in common
readIDs are reported to TMP ${TMPdir}multi-mappers_with_shared-coordinates.txt
" >> ${libFOLDER}log.txt
rm -rf ${TMPdir}multi-mappers_with_shared-coordinates.txt

  awk -v OFS="\t" -v OUTFILE=${libFOLDER}log.txt -v XFILE=${TMPdir}multi-mappers_with_shared-coordinates.txt '{
    if($10>1 && $4~"mapping=m") {
      X[$4][$2]+=1  
      Y[$4][$3]+=1
    }
  } 
  END{ 
    for(i in X) { 
      for (j in X[i]) {
        if(X[i][j]>1) {
          countX+=1
          TOTcountX+=X[i][j] 
          print i >> XFILE
        }
      }
    }
    for(i in Y) { 
      for (j in Y[i]) {
        if(Y[i][j]>1) {
          countY+=1
          TOTcountY+=Y[i][j] 
          print i >> XFILE
        }
      }
    }
    print "same-start","#sequences="countX,"#mappings="TOTcountX >>OUTFILE
    print "same-end","#sequences="countY,"#mappings="TOTcountY >>OUTFILE
  }' ${locTMP}${NAME}_mapped.bed

###################################################################################################
#clean up
if [[ $DEBUG != Y ]]; then
  rm -rf $locTMP
fi

PROCESSED_TIME=$(echo -e "$(date "+%s")" $TIME | awk '{ print ($1-$2)/60 }' )
printf "map reads - processing_time=	${PROCESSED_TIME} \n">> ${libFOLDER}time-log.txt
