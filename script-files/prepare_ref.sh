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

###################################################################################################
#add new annotation to versions.txt
#test if file is not empty and calculate the next ID number
if [[ -s ${BASE_UTILITY_LOCATION}versions.txt ]]; then
  CHECK=$(cat ${BASE_UTILITY_LOCATION}versions.txt | grep $VERSION)
  if [[ -z ${CHECK} ]]; then
    lastID=$(cat ${BASE_UTILITY_LOCATION}versions.txt | sort -k1,1n | tail -n -1 | tr ' ' '\t' | cut -f 1)
    newID=$(($lastID + 1))
    echo $newID $VERSION $RELEASE | tr ' ' '\t' >>${BASE_UTILITY_LOCATION}versions.txt
  fi
#if file was empty start with ID=1
else
  newID=1
  echo $newID $VERSION $RELEASE | tr ' ' '\t' >>${BASE_UTILITY_LOCATION}versions.txt
fi

###################################################################################################
#!!!!!!!!!!!!!   name conversion for dm6 only !!!!!!!!!!!!!!!!!!!!
if [[ $GENOME_VERSION == dm6 ]]; then
  #-------------------------------------------------------------------------------------------------
  #get chromosome name conversion file
   #?removed download as the ftp download is not working for peter 
   #?curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_assembly_report.txt |
   cat ${UTILITY_LOCATION}name-conversion.dm6.txt |
    #reformat file into conversion table
    awk -v OFS="\t" '
     {
       if($1 !~ "#") {
          print $1,$NF
       }
     }' | tr -d '\015' | sed 's/chr//' >${BASE_UTILITY_LOCATION}name-conversion.txt
fi

###################################################################################################
#run annotation preparation for selected FlyBase version
cd ${locTMP}
echo $locTMP
#download & unzip FlyBase gff file
  #@ curl --output ${locTMP}dmel-all-${VERSION}.gff.gz ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_${VERSION}_${RELEASE}/gff/dmel-all-${VERSION}.gff.gz 
  gunzip -c ${prepANNOTATIONgff} |
  #remove entries that belong to flamenco FBgn0267704
  grep -a -E -v "FBtr0347221|FBtr0347222|FBtr0347223|FBtr0347224|FBtr0347225|FBtr0347226|FBtr0347227|FBtr0347467|FBtr0347468|FBtr0347469" |
  awk -v OFS="\t" '
  {
    if($2=="FlyBase") {
      print
    } 
  }' >${locTMP}dmel-all-${VERSION}.gff

GFF=${locTMP}dmel-all-${VERSION}.gff

#extract annoations from FlyBase gff file
cat $GFF |
  awk -v OFS="\t" -v GFF=$GFF '
    #read in GFF file and put FlyBase entries into array X
    BEGIN{
      #read in file
      while((getline I < GFF ) > 0) {
        #split into array by tabs
        split(I,splitLINE,/\t/);
        #extract FlyBase transcripts and put into array
        if(splitLINE[2]=="FlyBase" ) {
          split(splitLINE[9],splitTAG,/;/);
          if(splitTAG[1]~"ID=FBtr") {
            split(splitTAG[1],splitID,/=/)
            X[splitID[2]]=splitLINE[3]
            Y[splitID[2]]=splitLINE[9]
          }
        }
      }
    }
    #filter for usable annotations
    {
      #convert to 0 based coordinates
      $4=$4-1
      #only consider FlyBase annotations not located on mitochondrion
      if($2=="FlyBase" && $1!~"mitochondrion_genome") {
        #if in following category extract corresponding parent class from array X
        if($3 ~ /^exon$|^intron$|^five_prime_UTR$|^three_prime_UTR$/ ) {
          split($9,splitTAG,/;/);
          for(j in splitTAG) {
            if(splitTAG[j]~"Parent"){
              split(splitTAG[j],parentID,/=|,/)
              #print corresponding to parent class
              if(X[parentID[2]] ~ /^mRNA$|^pseudogene$|^ncRNA$/) {
                if(X[parentID[2]]=="ncRNA" && Y[parentID[2]]~"hpRNA:"){
                  print $1,$4,$5,"hpRNA",0,$7
                }else{
                  print $1,$4,$5,X[parentID[2]]"::EX:"$3,0,$7
                }
              } else {
                if(X[parentID[2]] == "tRNA") {
                  if($4<20) print $1,0,$5+20,X[parentID[2]],0,$7
                  else print $1,$4-20,$5+20,X[parentID[2]],0,$7
                } else {
                  print $1,$4,$5,X[parentID[2]],0,$7
                }
              }
            }
          }
        } else {
          #create annotation for miRNA stuff
          if($3 ~ /^miRNA$|^pre_miRNA$/) {
            print $1,$4,$5,$3,0,$7
          }
        }
      }
    } ' | grep -F -v dmel_mitochondrion_genome >${locTMP}annotations-tmp.bed
if [[ $GENOME_VERSION == dm6 ]]; then
  #convert to UCSC chromosome names
  awk -v OFS="\t" '
       FNR == NR {
           assoc[ $1 ] = $2;
           next;
       }
       FNR < NR {
           for ( i = 1; i <= NF; i++ ) {
               if ( $i in assoc ) {
                   $i = assoc[ $i ]
               }
           }
           print $0
       }' ${BASE_UTILITY_LOCATION}name-conversion.txt ${locTMP}annotations-tmp.bed >${locTMP}annotations.bed
else
  mv ${locTMP}annotations-tmp.bed ${locTMP}annotations.bed
fi

###################################################################################################
#if bowtie index does not exist - create new

if ! ls ${UTILITY_LOCATION}/genome*.ebwt 1>/dev/null 2>&1; then
  #download flybase fasta file
  #!@ replaced by prepGENOMEfasta
  #!@  curl --output ${locTMP}full-genome.fa.gz ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_${VERSION}_${RELEASE}/fasta/dmel-all-chromosome-${VERSION}.fasta.gz

  #rename chromosomes for dm6 and filter Uextra for dm3
  if [[ $GENOME_VERSION == dm6 ]]; then
    seqkit replace -p "\s.+" ${prepGENOMEfasta} |
      seqkit replace --line-width 0 -p "(^.+)" -r '{kv}' -k ${BASE_UTILITY_LOCATION}name-conversion.txt >${locTMP}genome-all.fa
  elif [[ $GENOME_VERSION == dm3 ]]; then
    seqkit grep --line-width 0 -v -r -p Uextra ${prepGENOMEfasta} >${locTMP}genome-all.fa
  fi

  #extract everything except mitochondrial
  seqkit grep --line-width 0 -v -r -p ^M ${locTMP}genome-all.fa >${UTILITY_LOCATION}genome_no-Mito_incl-Y.fa
  bowtieBuild --threads $CORES ${UTILITY_LOCATION}genome_no-Mito_incl-Y.fa ${UTILITY_LOCATION}genome_no-Mito_incl-Y

  #extract everything except mitochondrial and eliminate Y-chromosome
  seqkit grep --line-width 0 -v -r -p ^M -p ^Y ${locTMP}genome-all.fa >${UTILITY_LOCATION}genome_no-Mito_excl-Y.fa
  bowtieBuild --threads $CORES ${UTILITY_LOCATION}genome_no-Mito_excl-Y.fa ${UTILITY_LOCATION}genome_no-Mito_excl-Y

  #extract mitochondrial DNA
  seqkit grep --line-width 0 -r -p ^M ${locTMP}genome-all.fa >${UTILITY_LOCATION}mitochondrion_DNA.fa
  bowtieBuild --threads $CORES ${UTILITY_LOCATION}mitochondrion_DNA.fa ${UTILITY_LOCATION}mitochondrion_DNA

  #prepare index for synthetic rRNA precursor
  seqkit seq --line-width 0 ${UTILITY_LOCATION}rRNA_precursor.fa >${UTILITY_LOCATION}rRNA.fa
  bowtieBuild --threads $CORES ${UTILITY_LOCATION}rRNA.fa ${UTILITY_LOCATION}rRNA

  #create chrom_sizes file
  seqkit fx2tab --name --length --threads $CORES ${locTMP}genome-all.fa |
    awk -v OFS="\t" '{ print "chr"$0 }'>${UTILITY_LOCATION}chrom.sizes
  seqkit fx2tab --name --length --threads $CORES ${UTILITY_LOCATION}rRNA.fa |
    awk -v OFS="\t" '{ print "chr"$0 }'>>${UTILITY_LOCATION}chrom.sizes
fi

#add mitochondrion genome into annotation
seqkit fx2tab --name --length --threads $CORES ${UTILITY_LOCATION}mitochondrion_DNA.fa |
  awk -v OFS="\t" '{ print $1,0,$2,"dmel_mitochondrion_genome",0,"+"}' >>${locTMP}annotations.bed

#add rRNA precursor to annotation
seqkit fx2tab --name --length --threads $CORES ${UTILITY_LOCATION}rRNA.fa |
  awk -v OFS="\t" '{ if( $1~"rRNA_precursor" )print $1,0,$2,"rRNA_precursor",0,"+" }' >>${locTMP}annotations.bed

#add RepeatMasker TEs to annotations
awk -v OFS="\t" '{ if( $12=="LINE" || $12=="LTR" || $12=="Satellite" || $12=="DNA" ) print $6,$7,$8,"TE::TE:"$12,0,$10}' ${UTILITY_LOCATION}UCSC-RepeatMasker.full | sed 's/chr//' >>${locTMP}annotations.bed

#add RepeatMasker rRNAs to annotations
awk -v OFS="\t" '{ if( $12~"rRNA" ) print $6,$7,$8,"rRNA",0,$10}' ${UTILITY_LOCATION}UCSC-RepeatMasker.full |
  sed 's/chr//' >>${locTMP}annotations.bed

#?added hpRNAs - therefore siRNA clusters are not required anymore
#?also removed siRNA_cluster from ruleset
#?#add siRNA cluster
#?cat ${UTILITY_LOCATION}siRNA_cluster.bed | awk -v OFS="\t" '{ gsub("chr",""); $4="siRNA_cluster"; print }' >>${locTMP}annotations.bed

#sort bedfile
LC_ALL=C sort -k1,1 -k2,2n ${locTMP}annotations.bed >${UTILITY_LOCATION}${VERSION}.bed

###################################################################################################
#create transcriptome indexes

##download transcriptome fasta and pre-process
for refTYPE in transcript CDS; do
  if [[ $refTYPE == "transcript" ]]; then
    cp ${prepTRANSCRIPTOMEfasta} ${locTMP}dmel-all-${refTYPE}.fasta.gz
  fi
  if [[ $refTYPE == "CDS" ]]; then
    cp ${prepCDSfasta} ${locTMP}dmel-all-${refTYPE}.fasta.gz
  fi

  #@ curl --output ${locTMP}dmel-all-${refTYPE}.fasta.gz ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_${VERSION}_${RELEASE}/fasta/dmel-all-${refTYPE}-${VERSION}.fasta.gz
  
  if [[ $refTYPE == transcript ]]; then
    #@ curl --output ${locTMP}dmel-all-ncRNA.fasta.gz ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_${VERSION}_${RELEASE}/fasta/dmel-all-ncRNA-${VERSION}.fasta.gz
    cat ${prepNCRNAfasta} >> ${locTMP}dmel-all-${refTYPE}.fasta.gz
  fi
 

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
          print ">"trNAME"\n"$NF > locTMP "RNAs." refTYPE ".tmp.fa"

          split(NAME,splitNAME,/=|-R[A-Z]/)
          print trNAME, splitNAME[2] > locTMP refTYPE "_to_gene.orig.txt"
        }
      }' 

  minimap2 -t 8 -x asm5  ${locTMP}RNAs.${refTYPE}.tmp.fa  ${locTMP}RNAs.${refTYPE}.tmp.fa >  ${locTMP}RNAs.${refTYPE}.paf

  mawk -v OFS="\t" -v TMP=$locTMP -v refTYPE=${refTYPE} '
  {
    split($1,split1,/::/)
    split($6,split2,/::/)
    if($3<40 && $4>$2-40 && $8<40 && $9>$7-40 && $11 > $7-80 && split1[4]!=split2[4]){
      if( $1 in ALL && $6 in ALL){
        x=b
      }else{
        if($1 in X){
          X[$1]=X[$1]":~~:"$6
          ALL[$6]+=1
        }else{
          if($6 in X){
            X[$6]=X[$6]":~~:"$1
            ALL[$1]+=1
          }else{
            if(split1[4] in PARENTS){
                X[$1]=$6
                ALL[$1]+=1
                ALL[$6]+=1
                PARENTS[split1[4]]+=1
                a+=1
            }else{
              if(split2[4] in PARENTS){
                X[$6]=$1
                ALL[$1]+=1
                ALL[$6]+=1
                PARENTS[split2[4]]+=1
                b+=1
              }else{
                X[$1]=$6
                ALL[$1]+=1
                ALL[$6]+=1
                PARENTS[split1[4]]+=1
                c+=1
              }
            }
          }
        }
      }
    }
  }
  END{
    #@ print a,b,c
    for (i in X){
      print i,X[i]
      print i > TMP "group-representative.txt"
      split(X[i],splitX,/:~~:/)
      for (j in splitX){
        print splitX[j] > TMP "IDs_to_filter."refTYPE".txt"
      }
    }
  }' ${locTMP}RNAs.${refTYPE}.paf  > ${locTMP}${refTYPE}-filter.association.txt 

  seqkit grep --line-width 0 -r --invert-match -f ${locTMP}IDs_to_filter.${refTYPE}.txt ${locTMP}RNAs.${refTYPE}.tmp.fa > ${locTMP}RNAs.${refTYPE}.fa

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
  }' ${UTILITY_LOCATION}TE_annot_new_WO_DUST.fa >${locTMP}TE_reads.bed

  #get fasta and map to the transcriptome index
  redCORES=$(( CORES - 2 ))

  rm -rf ${UTILITY_LOCATION}TE_annot_new_WO_DUST.fa.fai
  bedtools getfasta -name -fi ${UTILITY_LOCATION}TE_annot_new_WO_DUST.fa -bed ${locTMP}TE_reads.bed -fo - |
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
    seqkit grep -v -p stellate -p Su-Ste --line-width 0 ${UTILITY_LOCATION}TE_annot_new_WO_DUST.fa >${locTMP}RNAs+TE${finEXT}.fa
    seqkit grep -v -p stellate -p Su-Ste --line-width 0 ${UTILITY_LOCATION}TE_annot_new_WO_DUST.fa >${locTMP}RNAs+TE${finEXT}.unistrand.fa
    #antisense TE sequences for salmon
    seqkit grep -v -p stellate -p Su-Ste ${UTILITY_LOCATION}TE_annot_new_WO_DUST.fa |
      seqkit seq -p -r --line-width 0 |
      awk -v OFS="\t" '{ if($1~">") print $1"_AS"; else print }' >>${locTMP}RNAs+TE${finEXT}.fa
    #TE sequences for STAR and bowtie
    seqkit grep -v -p stellate -p Su-Ste --line-width 0 ${UTILITY_LOCATION}TE_annot_new_WO_DUST.fa >${locTMP}TE_noSplice.fa

    #create TE length file
    seqkit grep -v -p stellate -p Su-Ste ${UTILITY_LOCATION}TE_annot_new_WO_DUST.fa |
      seqkit fx2tab --name --length | 
      #remove excessive tabs
      tr -s "\t" >${UTILITY_LOCATION}TE${finEXT}.sizes
  else
    EXT=""
    finEXT="_incl-Y"
    #sense TE sequences for salmon
    cat ${UTILITY_LOCATION}TE_annot_new_WO_DUST.fa >${locTMP}RNAs+TE${finEXT}.fa
    cat ${UTILITY_LOCATION}TE_annot_new_WO_DUST.fa >${locTMP}RNAs+TE${finEXT}.unistrand.fa
    #antisense TE sequences for salmon
    seqkit seq -p -r --line-width 0 ${UTILITY_LOCATION}TE_annot_new_WO_DUST.fa |
      awk -v OFS="\t" '{ if($1~">") print $1"_AS"; else print }' >>${locTMP}RNAs+TE${finEXT}.fa
    #TE sequences for STAR and bowtie
    cat ${UTILITY_LOCATION}TE_annot_new_WO_DUST.fa >${locTMP}TE_noSplice.fa
    #create TE length file
    seqkit fx2tab --name --length ${UTILITY_LOCATION}TE_annot_new_WO_DUST.fa | 
      #remove excessive tabs
      tr -s "\t" >${UTILITY_LOCATION}TE${finEXT}.sizes
  fi
  
  #add transcript sequences to TE sequences
  cat ${locTMP}maskedRNAs${EXT}.transcript.fa >>${locTMP}RNAs+TE${finEXT}.fa
  cat ${UTILITY_LOCATION}genome_no-Mito${finEXT}.fa >>${locTMP}RNAs+TE${finEXT}.fa
  cat ${UTILITY_LOCATION}genome_no-Mito${finEXT}.fa >>${locTMP}maskedRNAs${EXT}.CDS.fa

  cat ${locTMP}maskedRNAs${EXT}.transcript.fa >>${locTMP}RNAs+TE${finEXT}.unistrand.fa
  cat ${UTILITY_LOCATION}genome_no-Mito${finEXT}.fa >>${locTMP}RNAs+TE${finEXT}.unistrand.fa

  seqkit fx2tab ${UTILITY_LOCATION}genome_no-Mito${finEXT}.fa | cut -f 1 |
    sed -e 's/>//g' > ${locTMP}decoys.txt

  #create salmon indexes
  salmon index -t ${locTMP}RNAs+TE${finEXT}.fa -k 21 -d ${locTMP}decoys.txt --threads $CORES -i ${UTILITY_LOCATION}salmon_quasi_${VERSION}${finEXT}
  salmon index -t ${locTMP}RNAs+TE${finEXT}.fa -k 11 -d ${locTMP}decoys.txt --threads $CORES -i ${UTILITY_LOCATION}salmon_quasi-sRNA_${VERSION}${finEXT}

  #create index for unistrand salmon quantification
  salmon index -t ${locTMP}RNAs+TE${finEXT}.unistrand.fa -k 21 -d ${locTMP}decoys.txt --threads $CORES -i ${UTILITY_LOCATION}salmon_quasi_${VERSION}${finEXT}_unistrand

  #index CDS file
  salmon index -t ${locTMP}maskedRNAs${EXT}.CDS.fa -k 11 -d ${locTMP}decoys.txt --threads $CORES -i ${UTILITY_LOCATION}salmon_quasi-CDS_${VERSION}${finEXT}

  #create indices for TE histograms
  rm -rf ${UTILITY_LOCATION}star_TE${finEXT}/
  mkdir ${UTILITY_LOCATION}star_TE${finEXT}/
  STAR --runThreadN 5 --runMode genomeGenerate --genomeDir ${UTILITY_LOCATION}star_TE${finEXT}/ \
    --genomeFastaFiles ${locTMP}TE_noSplice.fa --genomeSAindexNbases 8 \
    --sjdbGTFtagExonParentTranscript Parent --sjdbFileChrStartEnd ${UTILITY_LOCATION}TE_splice-junctions.txt

  bowtieBuild --noref --threads $CORES ${locTMP}TE_noSplice.fa ${UTILITY_LOCATION}bowtie-TE${finEXT}

  cp ${locTMP}RNAs+TE${finEXT}.fa ${UTILITY_LOCATION}RNAs+TE${finEXT}.fa
done

#merge transcript-to-gene files
cat ${UTILITY_LOCATION}TE_annot_new_WO_DUST.fa | awk '{if($1~">") { sub(">",""); print $1,"TE:"$1; print $1"_AS","TE:"$1"_AS" }}' >${UTILITY_LOCATION}transcript_to_gene_${VERSION}.txt
cat ${locTMP}transcript_to_gene.orig.txt >>${UTILITY_LOCATION}transcript_to_gene_${VERSION}.txt
cp ${locTMP}CDS_to_gene.orig.txt ${UTILITY_LOCATION}CDS_to_gene_${VERSION}.txt


###################################################################################################
#create index for STAR
rm -rf ${UTILITY_LOCATION}star_${VERSION}/
mkdir -p ${UTILITY_LOCATION}star_${VERSION}/

#convert GFF chromosome names to UCSC names
if [[ $GENOME_VERSION == dm6 ]]; then
  awk -v OFS="\t" -v LOC="${UTILITY_LOCATION}star_${VERSION}/" '
      FNR == NR {
          assoc[ $1 ] = $2;
          next;
      }
      FNR < NR {
        for ( i = 1; i <= NF; i++ ) {
          if ( $i in assoc ) {
            $i = assoc[ $i ]
          }
        }
        if($0 !~ /dmel_mitochondrion_genome/ ) {
          if($2=="FlyBase" && $3=="intron") {
            if($1!~"Y") {
              print $1,$4,$5,$7 > LOC "junctions_excl-Y.txt"
            } 
            print $1,$4,$5,$7 > LOC "junctions_incl-Y.txt"
          }
          print $0
        }
      }' ${BASE_UTILITY_LOCATION}name-conversion.txt ${locTMP}dmel-all-${VERSION}.gff >${UTILITY_LOCATION}star_${VERSION}/${VERSION}.gff
else
  awk -v OFS="\t" -v LOC="${UTILITY_LOCATION}star_${VERSION}/" '
  {
    if($0 !~ /dmel_mitochondrion_genome/ ) {
      if($2=="FlyBase" && $3=="intron") {
        if($1!~"Y") {
          print $1,$4,$5,$7 > LOC "junctions_excl-Y.txt"
        }
        print $1,$4,$5,$7 > LOC "junctions_incl-Y.txt"
      }
      print $0
    }
  }' <${locTMP}dmel-all-${VERSION}.gff >${UTILITY_LOCATION}star_${VERSION}/${VERSION}.gff
fi

mkdir -p ${UTILITY_LOCATION}star_${VERSION}/incl-Y
STAR --runThreadN 5 --runMode genomeGenerate --genomeDir ${UTILITY_LOCATION}star_${VERSION}/incl-Y/ \
  --genomeFastaFiles ${UTILITY_LOCATION}genome_no-Mito_incl-Y.fa \
  --sjdbGTFtagExonParentTranscript Parent --genomeSAindexNbases 12 \
  --sjdbFileChrStartEnd ${UTILITY_LOCATION}star_${VERSION}/junctions_incl-Y.txt

mkdir -p ${UTILITY_LOCATION}star_${VERSION}/excl-Y
STAR --runThreadN 5 --runMode genomeGenerate --genomeDir ${UTILITY_LOCATION}star_${VERSION}/excl-Y/ \
  --genomeFastaFiles ${UTILITY_LOCATION}genome_no-Mito_excl-Y.fa \
  --sjdbGTFtagExonParentTranscript Parent --genomeSAindexNbases 12 \
  --sjdbFileChrStartEnd ${UTILITY_LOCATION}star_${VERSION}/junctions_excl-Y.txt

###################################################################################################
if [[ $DEBUG != Y ]]; then
  rm -rf "$locTMP"
fi

echo $locTMP
exit
