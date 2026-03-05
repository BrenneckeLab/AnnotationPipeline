#!/usr/bin/env bash

#SBATCH --cpus-per-task=4
#SBATCH --mem=10g
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

locTMP="${TMPdir}prepare-ref/"
# rm -rf $locTMP
mkdir -p $locTMP

export TMPDIR=$locTMP

rm -rf $UTILITY_LOCATION
mkdir -p $UTILITY_LOCATION
###################################################################################################
#add new annotation to versions.txt
#test if file is not empty and calculate the next ID number
if [[ -s ${BASE_UTILITY_LOCATION}versions.txt ]]; then
  CHECK=$(cat ${BASE_UTILITY_LOCATION}versions.txt | grep $VERSION)
  if [[ -z ${CHECK} ]]; then
    lastID=$(cat ${BASE_UTILITY_LOCATION}versions.txt | sort -k1,1n | tail -n -1 | tr ' ' '\t' | cut -f 1)
    newID=$(($lastID + 1))
    echo $newID $VERSION ASM | tr ' ' '\t' >>${BASE_UTILITY_LOCATION}versions.txt
  fi
#if file was empty start with ID=1
else
  newID=1
  echo $newID $VERSION ASM | tr ' ' '\t' >>${BASE_UTILITY_LOCATION}versions.txt
fi

###################################################################################################
#create assembly-hub
rm -rf ${ASMdir}
mkdir -p ${ASMdir}

#copy assembly to asmdir
time cp -r ${asmHUBpath}/* ${ASMdir}/

if [[ $VERSION != $ASMname ]]; then
  mv  ${ASMdir}${ASMname}  ${ASMdir}${VERSION}
  #change name of assembly in all relevant files
  for FILE in ${ASMdir}hub.txt ${ASMdir}genomes.txt ${ASMdir}${VERSION}/trackDb.txt; do
    sed -i "s/$ASMname/$VERSION/g" $FILE
  done

  #rename files in first level folder
  ALLtoCHANGE=$(ls ${ASMdir}${VERSION} | grep ${ASMname} | tr '\n' '\t')
  for FILE in $ALLtoCHANGE; do
    newNAME=$(echo $FILE | sed "s/$ASMname/$VERSION/")
    mv ${ASMdir}${VERSION}/$FILE ${ASMdir}${VERSION}/$newNAME
  done
  #rename files in annotations level folder
  ALLtoCHANGE=$(ls ${ASMdir}${VERSION}/annotations/ | grep ${ASMname} | tr '\n' '\t')
  for FILE in $ALLtoCHANGE; do
    newNAME=$(echo $FILE | sed "s/$ASMname/$VERSION/")
    mv ${ASMdir}${VERSION}/annotations/$FILE ${ASMdir}${VERSION}/annotations/$newNAME
  done
fi

while read POS TYPE REST; do
  PRIO=$(( 40 + $POS ))
  printf "
    name AP_$TYPE
    label AP_$TYPE
    priority $PRIO
    defaultIsClosed 0 \n\n" | sed '/./,$!d' >> ${ASMdir}groups.txt

    if [[ $TYPE == sRNAseq ]]; then
      PRIOx=$(( PRIO + 20 ))
      printf "
        name AP_${TYPE}_piRNA
        label AP_${TYPE}_piRNA
        priority $PRIOx
        defaultIsClosed 1 \n\n" | sed '/./,$!d' >> ${ASMdir}groups.txt
      PRIOx=$(( PRIO + 21 ))
      printf "
        name AP_${TYPE}_siRNA
        label AP_${TYPE}_siRNA
        priority $PRIOx
        defaultIsClosed 1 \n\n" | sed '/./,$!d' >> ${ASMdir}groups.txt
    fi
    if [[ $TYPE == RNAseq || $TYPE == RIPseq || $TYPE == CLIPseq ]]; then 
      PRIOx=$(( PRIO + 20 ))
      printf "
        name AP_${TYPE}_spliceJunctions
        label AP_${TYPE}_spliceJunctions
        priority $PRIOx
        defaultIsClosed 1 \n\n" | sed '/./,$!d' >> ${ASMdir}groups.txt

    fi
done <${BASE_UTILITY_LOCATION}available_TYPES.txt

###################################################################################################
#run annotation preparation for selected FlyBase version
cd ${locTMP}
rm -rf ${locTMP}annotations.bed

#convert gene-track to bed
bigBedToBed ${ASMdir}${VERSION}/annotations/HQmergedAnnotations.bb   ${locTMP}HQ-annotations.bigGenePred

#extract transcripts
for TRANSCRIPTclass in mRNA pseudogene ncRNA; do
  mawk -v OFS="\t" -v CLASS=$TRANSCRIPTclass '
  {
    if($17 ~ CLASS ){
      print 
    }
  }' ${locTMP}HQ-annotations.bigGenePred | cut -f 1-12 | 
    #remove flamenco entries
    grep -v "lncRNA:flam" > ${locTMP}${TRANSCRIPTclass}.bed12

  #extract exon annotations
  bedtools bed12tobed6 -i ${locTMP}${TRANSCRIPTclass}.bed12 | 
    mawk -v OFS="\t" -v CLASS=$TRANSCRIPTclass '
    {
      $4=CLASS"::EX:exon"
      print
    }' >> ${locTMP}annotations.bed


  #prepare additional annotations related to transcripts
  for TYPE in 3pUTR 5pUTR introns; do
    if [[ $TYPE == 3pUTR ]]; then ANN=three_prime_UTR ; fi
    if [[ $TYPE == 5pUTR ]]; then ANN=five_prime_UTR ; fi
    if [[ $TYPE == introns ]]; then ANN=intron ; fi

    bedparse $TYPE ${locTMP}${TRANSCRIPTclass}.bed12 > ${locTMP}${TRANSCRIPTclass}.${TYPE}.bed12
      bedtools bed12tobed6 -i  ${locTMP}${TRANSCRIPTclass}.${TYPE}.bed12 |   
      mawk -v OFS="\t" -v ANN=$ANN -v CLASS=$TRANSCRIPTclass '
      {
        $4=CLASS"::EX:"ANN
        print
      }' >> ${locTMP}annotations.bed
  done
done
   
for TRANSCRIPTclass in tRNA snRNA snoRNA rRNA miRNA hpRNA; do
  mawk -v OFS="\t" -v CLASS=$TRANSCRIPTclass '
  {
    if($17~CLASS ){
      print $1,$2,$3,$17,0,$6
    }
  }' ${locTMP}HQ-annotations.bigGenePred >>  ${locTMP}annotations.bed
done

#add RepeatMasker TEs to annotations
for CLASS in LINE LTR Satellite DNA; do
  bigBedToBed  ${ASMdir}${VERSION}/annotations/${CLASS}.bb stdout |
    awk -v OFS="\t" -v CLASS=$CLASS '{
      $4="TE::TE:"CLASS
      print $1,$2,$3,"TE::TE:"CLASS,0,$6
    }' >>${locTMP}annotations.bed
done

###################################################################################################
#if bowtie index does not exist - create new
if ! ls ${UTILITY_LOCATION}/*.ebwt 1>/dev/null 2>&1; then
  #remove mito from genome and create mito stand-alone sequence
  seqkit grep -r -p Mito --invert-match ${ASMdir}${VERSION}/${VERSION}.fa >${UTILITY_LOCATION}${VERSION}.genome.fa

  seqkit grep -r -p Mito  ${ASMdir}${VERSION}/${VERSION}.fa >${UTILITY_LOCATION}${VERSION}.mito.fa
  if [[ ! -s ${ASMdir}${VERSION}/annotations/mitochondrion_DNA.bed ]]; then
    cp ${BASE_UTILITY_LOCATION}/ASM/mitochondrial_DNA.fa ${UTILITY_LOCATION}${VERSION}.mito.fa
  fi

  #build genmome bowtie index
  bowtieBuild --threads $CORES ${UTILITY_LOCATION}${VERSION}.genome.fa ${UTILITY_LOCATION}genome_${GENOME_VERSION}

  #build mito bowtie index
  bowtieBuild --threads $CORES ${UTILITY_LOCATION}${VERSION}.mito.fa ${UTILITY_LOCATION}mitochondrion_DNA

  #prepare index for synthetic rRNA precursor
  bowtieBuild --threads $CORES ${BASE_UTILITY_LOCATION}ASM/rRNA_precursor.fa ${UTILITY_LOCATION}rRNA

  #add mito and rRNA to chrom_sizes file
  seqkit fx2tab --name --length --threads $CORES ${UTILITY_LOCATION}${VERSION}.genome.fa >${UTILITY_LOCATION}chrom.sizes
  seqkit fx2tab --name --length --threads $CORES ${BASE_UTILITY_LOCATION}ASM/rRNA_precursor.fa >>${UTILITY_LOCATION}chrom.sizes
  seqkit fx2tab --name --length --threads $CORES ${UTILITY_LOCATION}${VERSION}.mito.fa >>${UTILITY_LOCATION}chrom.sizes
fi

seqkit fx2tab --name --length --threads $CORES ${UTILITY_LOCATION}${VERSION}.mito.fa |
  awk -v OFS="\t" '{ print $1,0,$2,"dmel_mitochondrion_genome",0,"+"}' >>${locTMP}annotations.bed

#add rRNA precursor to annotation
seqkit fx2tab --name --length --threads $CORES ${BASE_UTILITY_LOCATION}ASM/rRNA_precursor.fa |
  awk -v OFS="\t" '{ if( $1~"rRNA_precursor" )print $1,0,$2,"rRNA_precursor",0,"+" }' >>${locTMP}annotations.bed

###################################################################################################
#create transcriptome indexes

##download transcriptome fasta and pre-process
for refTYPE in transcript CDS; do
   wget -O ${locTMP}dmel-all-${refTYPE}.fasta.gz ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current/fasta/dmel-all-${refTYPE}*.fasta.gz 

  #prepare transcript fasta and transcript-to-gene file for mRNAs
  seqkit fx2tab ${locTMP}dmel-all-${refTYPE}.fasta.gz |
    awk -v locTMP=$locTMP -v refTYPE=$refTYPE '
      { 
      split($0,LINE,/; /)
        for (FIELD in LINE) {
          if (LINE[FIELD] ~ "length=") { SIZE=LINE[FIELD] }
          if (LINE[FIELD] ~ "name=") { NAME=LINE[FIELD] }
          if (LINE[FIELD] ~ "parent=") { n=split(LINE[FIELD],splitPARENT,/,/); PARENT=splitPARENT[1]; if(n==2) trID=splitPARENT[2]; else trID=$1 }
        
        } 
        if ($0 ~ "loc=Y") { Ystat="Y";}  else {Ystat="N"}


        trNAME=trID"::"SIZE"::"NAME"::"PARENT"::Y="Ystat
        if( NAME !~ "mt:") {
          print ">"trNAME"\n"$NF > locTMP "RNAs." refTYPE ".fa"

          split(NAME,splitNAME,/=|-R[A-Z]/)
          print trNAME, splitNAME[2] > locTMP refTYPE "_to_gene.orig.txt"
        }
      }' 

  #masking of TE sequences in the transcript file
  ##generate bowtie index for transcriptome
  bowtieBuild --noref --threads $CORES ${locTMP}RNAs.${refTYPE}.fa ${locTMP}${refTYPE}.index

  #generate read-coordinates from TE-consensus
  awk -v RS=">" -v OFS="\t" '
  {
    sub("\n","\t",$0)
    gsub("\n","",$0)
    TElength=length($NF)
    maxCOORD=TElength-30
    for(i=0; i<=maxCOORD;i++) {
      print $1, i,i+30, $1"_"i"_+",0,"+"
    }
  }' ${BASE_UTILITY_LOCATION}ASM/TE_annot_new_WO_DUST.fa >${locTMP}TE_reads.bed

  #get fasta and map to the transcriptome index
  redCORES=$(( CORES - 2 ))
  bedtools getfasta -name -fi ${BASE_UTILITY_LOCATION}ASM/TE_annot_new_WO_DUST.fa -bed ${locTMP}TE_reads.bed -fo - |
    bowtie -a -p $redCORES -v 1 -f -S ${locTMP}${refTYPE}.index - |
    samtools view -bS - |
    bedtools bamtobed -i - >${locTMP}TE-regions.${refTYPE}.bed

  #mask fasta file
  bedtools maskfasta -fi ${locTMP}RNAs.${refTYPE}.fa -bed ${locTMP}TE-regions.${refTYPE}.bed -fo ${locTMP}maskedRNAs.${refTYPE}.tmp

  #generate Y-depleted RNA file
  awk -v RS=">" -v TMP=$locTMP -v refTYPE=${refTYPE} '
    {
      if(NR>1) {
        gsub("\n","\t")
        if($1 !~ "Y=Y") {
          print ">"$1"\n"$2  > TMP "maskedRNAs_noY."refTYPE".fa"
        }
        print ">"$1"\n"$2 > TMP "maskedRNAs."refTYPE".fa"
      }
    }
  ' ${locTMP}maskedRNAs.${refTYPE}.tmp

done

#remove old salmon index
rm -rf ${UTILITY_LOCATION}salmon_quasi_${VERSION}*

#add TE names to transcript-to-gene file
for TYPE in noY inclY; do
  if [[ $TYPE == noY ]]; then
    EXT="_noY"
    finEXT="_excl-Y"
    #sense TE sequences for salmon
    seqkit grep -v -p stellate -p Su-Ste --line-width 0 ${BASE_UTILITY_LOCATION}ASM/TE_annot_new_WO_DUST.fa >${locTMP}RNAs+TE${finEXT}.fa
    #antisense TE sequences for salmon
    seqkit grep -v -p stellate -p Su-Ste ${BASE_UTILITY_LOCATION}ASM/TE_annot_new_WO_DUST.fa |
      seqkit seq -p -r --line-width 0 |
      awk -v OFS="\t" '{ if($1~">") print $1"_AS"; else print }' >>${locTMP}RNAs+TE${finEXT}.fa
    #TE sequences for STAR and bowtie
    seqkit grep -v -p stellate -p Su-Ste --line-width 0 ${BASE_UTILITY_LOCATION}ASM/TE_annot_new_WO_DUST.fa >${locTMP}TE_noSplice.fa
    #create TE length file
    seqkit grep -v -p stellate -p Su-Ste ${BASE_UTILITY_LOCATION}ASM/TE_annot_new_WO_DUST.fa |
      seqkit fx2tab --name --length | 
      #remove excessive tabs
      tr -s "\t" >${UTILITY_LOCATION}TE${finEXT}.sizes
  else
    EXT=""
    finEXT="_incl-Y"
    #sense TE sequences for salmon
    cat ${BASE_UTILITY_LOCATION}ASM/TE_annot_new_WO_DUST.fa >${locTMP}RNAs+TE${finEXT}.fa
    #antisense TE sequences for salmon
    seqkit seq -p -r --line-width 0 ${BASE_UTILITY_LOCATION}ASM/TE_annot_new_WO_DUST.fa |
      awk -v OFS="\t" '{ if($1~">") print $1"_AS"; else print }' >>${locTMP}RNAs+TE${finEXT}.fa
    #TE sequences for STAR and bowtie
    cat ${BASE_UTILITY_LOCATION}ASM/TE_annot_new_WO_DUST.fa >${locTMP}TE_noSplice.fa
    #create TE length file
    seqkit fx2tab --name --length ${BASE_UTILITY_LOCATION}ASM/TE_annot_new_WO_DUST.fa | 
      #remove excessive tabs
      tr -s "\t" >${UTILITY_LOCATION}TE${finEXT}.sizes
  fi
  
  #add transcript sequences to TE sequences
  cat ${locTMP}maskedRNAs${EXT}.transcript.fa >>${locTMP}RNAs+TE${finEXT}.fa
  cat ${UTILITY_LOCATION}${VERSION}.genome.fa >>${locTMP}RNAs+TE${finEXT}.fa
  cat ${UTILITY_LOCATION}${VERSION}.genome.fa >>${locTMP}maskedRNAs${EXT}.CDS.fa

  seqkit fx2tab ${UTILITY_LOCATION}${VERSION}.genome.fa | cut -f 1 |
    sed -e 's/>//g' > ${locTMP}decoys.txt

  #create salmon indexes
  salmon index -t ${locTMP}RNAs+TE${finEXT}.fa -k 21 -d ${locTMP}decoys.txt --threads $CORES -i ${UTILITY_LOCATION}salmon_quasi_${VERSION}${finEXT}
  salmon index -t ${locTMP}RNAs+TE${finEXT}.fa -k 11 -d ${locTMP}decoys.txt --threads $CORES -i ${UTILITY_LOCATION}salmon_quasi-sRNA_${VERSION}${finEXT}
  #index CDS file
  salmon index -t ${locTMP}maskedRNAs${EXT}.CDS.fa -k 11 -d ${locTMP}decoys.txt --threads $CORES -i ${UTILITY_LOCATION}salmon_quasi-CDS_${VERSION}${finEXT}

  #create indices for TE histograms
  rm -rf ${UTILITY_LOCATION}star_TE${finEXT}/
  mkdir ${UTILITY_LOCATION}star_TE${finEXT}/
  STAR --runThreadN 5 --runMode genomeGenerate --genomeDir ${UTILITY_LOCATION}star_TE${finEXT}/ \
    --genomeFastaFiles ${locTMP}TE_noSplice.fa --genomeSAindexNbases 8 \
    --sjdbGTFtagExonParentTranscript Parent

  bowtieBuild --noref --threads $CORES ${locTMP}TE_noSplice.fa ${UTILITY_LOCATION}bowtie-TE${finEXT}

  cp ${locTMP}RNAs+TE${finEXT}.fa ${UTILITY_LOCATION}RNAs+TE${finEXT}.fa
done

#merge transcript-to-gene files
cat ${BASE_UTILITY_LOCATION}ASM/TE_annot_new_WO_DUST.fa | awk '{if($1~">") { sub(">",""); print $1,"TE:"$1; print $1"_AS","TE:"$1"_AS" }}' >${UTILITY_LOCATION}transcript_to_gene_${VERSION}.txt
cat ${locTMP}transcript_to_gene.orig.txt >>${UTILITY_LOCATION}transcript_to_gene_${VERSION}.txt
cp ${locTMP}CDS_to_gene.orig.txt ${UTILITY_LOCATION}CDS_to_gene_${VERSION}.txt


#@ ###################################################################################################
#@ #create index for STAR
#@ rm -rf ${UTILITY_LOCATION}star_${VERSION}/
#@ mkdir -p ${UTILITY_LOCATION}star_${VERSION}/

#@ #convert GFF chromosome names to UCSC names
#@ if [[ $GENOME_VERSION == dm6 ]]; then
#@   awk -v OFS="\t" -v LOC="${UTILITY_LOCATION}star_${VERSION}/" '
#@       FNR == NR {
#@           assoc[ $1 ] = $2;
#@           next;
#@       }
#@       FNR < NR {
#@         for ( i = 1; i <= NF; i++ ) {
#@           if ( $i in assoc ) {
#@             $i = assoc[ $i ]
#@           }
#@         }
#@         if($0 !~ /dmel_mitochondrion_genome/ ) {
#@           if($2=="FlyBase" && $3=="intron") {
#@             if($1!~"Y") {
#@               print $1,$4,$5,$7 > LOC "junctions_excl-Y.txt"
#@             } 
#@             print $1,$4,$5,$7 > LOC "junctions_incl-Y.txt"
#@           }
#@           print $0
#@         }
#@       }' ${BASE_UTILITY_LOCATION}name-conversion.txt ${locTMP}dmel-all-${VERSION}.gff >${UTILITY_LOCATION}star_${VERSION}/${VERSION}.gff
#@ else
#@   awk -v OFS="\t" -v LOC="${UTILITY_LOCATION}star_${VERSION}/" '
#@   {
#@     if($0 !~ /dmel_mitochondrion_genome/ ) {
#@       if($2=="FlyBase" && $3=="intron") {
#@         if($1!~"Y") {
#@           print $1,$4,$5,$7 > LOC "junctions_excl-Y.txt"
#@         }
#@         print $1,$4,$5,$7 > LOC "junctions_incl-Y.txt"
#@       }
#@       print $0
#@     }
#@   }' <${locTMP}dmel-all-${VERSION}.gff >${UTILITY_LOCATION}star_${VERSION}/${VERSION}.gff
#@ fi

#@ mkdir -p ${UTILITY_LOCATION}star_${VERSION}/incl-Y
#@ STAR --runThreadN 5 --runMode genomeGenerate --genomeDir ${UTILITY_LOCATION}star_${VERSION}/incl-Y/ \
#@   --genomeFastaFiles ${UTILITY_LOCATION}genome_no-Mito_incl-Y.fa \
#@   --sjdbGTFtagExonParentTranscript Parent --genomeSAindexNbases 12 \
#@   --sjdbFileChrStartEnd ${UTILITY_LOCATION}star_${VERSION}/junctions_incl-Y.txt

#@ mkdir -p ${UTILITY_LOCATION}star_${VERSION}/excl-Y
#@ STAR --runThreadN 5 --runMode genomeGenerate --genomeDir ${UTILITY_LOCATION}star_${VERSION}/excl-Y/ \
#@   --genomeFastaFiles ${UTILITY_LOCATION}genome_no-Mito_excl-Y.fa \
#@   --sjdbGTFtagExonParentTranscript Parent --genomeSAindexNbases 12 \
#@   --sjdbFileChrStartEnd ${UTILITY_LOCATION}star_${VERSION}/junctions_excl-Y.txt

###################################################################################################
#create index for STAR
rm -rf ${UTILITY_LOCATION}star_${VERSION}/
mkdir ${UTILITY_LOCATION}star_${VERSION}/

mkdir ${UTILITY_LOCATION}star_${VERSION}
STAR --runThreadN 5 --runMode genomeGenerate --genomeDir ${UTILITY_LOCATION}star_${VERSION}/ \
  --genomeFastaFiles ${UTILITY_LOCATION}${VERSION}.genome.fa \
  --sjdbGTFtagExonParentTranscript Parent --genomeSAindexNbases 12
# --sjdbFileChrStartEnd ${UTILITY_LOCATION}star_${VERSION}/junctions_incl-Y.txt

#@ ###################################################################################################
#@ #add flybase miRNAs to annotations
#@ wget --quiet -O ${locTMP}dmel-all-miRNA.fasta.gz ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current/fasta/dmel-all-miRNA-*.fasta.gz
#@ gunzip -c ${locTMP}dmel-all-miRNA.fasta.gz |
#@   awk -v RS=">" '
#@   {
#@     if(NR>1){
#@       sub("\n", "\t")
#@       gsub("\n", "")
#@       print ">"$1"\n"$NF
#@     }
#@   }
#@   ' >${locTMP}dmel-all-miRNA.fasta

#@ redCORES=$(( CORES - 2 ))
#@ STAR --runThreadN $redCORES --genomeDir ${UTILITY_LOCATION}star_${VERSION}/ \
#@   --readFilesIn ${locTMP}dmel-all-miRNA.fasta --outFilterMismatchNoverLmax 0.05 \
#@   --outFilterMatchNmin 14 --seedSearchStartLmax 20 --winAnchorMultimapNmax 2000 \
#@   --seedMultimapNmax 1000000 --scoreDelOpen 0 --scoreDelBase 0 --scoreInsOpen 0 \
#@   --scoreInsBase 0 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 \
#@   --alignIntronMax 1 --outReadsUnmapped Fastx --outFilterMultimapNmax 200 --outStd SAM |
#@   samtools view -bS - |
#@   bedtools bamtobed -split |
#@   awk -v OFS="\t" '{
#@     $4="miRNA"
#@     print 
#@   }' >>${locTMP}annotations.bed

#sort bedfile
LC_ALL=C sort -k1,1 -k2,2n --parallel=$CORES -S7G ${locTMP}annotations.bed >${UTILITY_LOCATION}${VERSION}.bed

###################################################################################################

if [[ $DEBUG != Y ]]; then
  rm -rf "$locTMP"
fi

echo $locTMP
exit
