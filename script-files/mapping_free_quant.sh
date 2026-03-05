#!/usr/bin/env bash
#$ -S /bin/bash

#$ -clear
#$ -e $JOB_NAME.e$JOB_ID.txt
#$ -o $JOB_NAME.o$JOB_ID.txt
#$ -cwd
#$ -q public.q
#$ -pe smp 10

#SBATCH --cpus-per-task=4
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
locTMP="${TMPdir}TMP_mapfree/"

mkdir -p $locTMP
chmod 777 $locTMP

export TMPDIR=$locTMP

###################################################################################################
#preperation of variables

if [[ $Ychrom == Y ]]; then INDEXext="_incl-Y"; else INDEXext="_excl-Y"; fi



redCORES=$(( CORES - 1 ))
if [[ $TYPE == sRNAseq ]]; then
    if [[ $GENOME_VERSION == Cel ]]; then    
      sizeCUTOFF=20
    else
      sizeCUTOFF=23
    fi


  gunzip -c ${TMPdir}${NAME}_annotated.fa.gz |
    mawk -v seed=$RANDOM -v sizeCUTOFF=$sizeCUTOFF '
    BEGIN{
      srand(seed)
    }
    {
      if($1~">" && $1!~"final_ann=rRNA|final_ann=tRNA|final_ann=snoRNA|final_ann=snRNA"){
        #?@ if($1~"mapping=u" || $1~"mapping=m") {
          split($1,splitNAME,/:|=/)
          #filter for piRNAs by length
          if(length(splitNAME[2])>sizeCUTOFF) {
            sub(">","")
            #print sequnce n-times according sequencing count
            for(i=1; i<=splitNAME[4]; i++) {
              print rand(),">" NR "-" i "_" $1,splitNAME[2]
            }
          }
        #?@ }
      }
    }'  |
    sort -k1,1n --parallel $redCORES -S15G |
    tr ' ' '\t' | cut -f 2-3 | tr '\t' '\n' >${locTMP}decollapsed.fa
else
  gunzip -c ${TMPdir}${NAME}_annotated.fa.gz |
    mawk -v seed=$RANDOM '
    BEGIN{
      srand(seed)
    }
    {
      if($1~">" && $1!~"final_ann=rRNA|final_ann=tR NA|final_ann=snoRNA|final_ann=snRNA"){
        #?@ if($1~"mapping=u" || $1~"mapping=m") {
          sub(">","")
          split($1,splitNAME,/:|=/)
          
          #print sequnce n-times according sequencing count
          for(i=1; i<=splitNAME[4]; i++) {
            print rand(),">" NR "-" i "_" $1,splitNAME[2]
          }
        #?@ }
      }
    }' | 
    sort -k1,1n --parallel $redCORES -S15G |
    tr ' ' '\t' | cut -f 2-3 | tr '\t' '\n' >${locTMP}decollapsed.fa
fi

#---------------------------------------------------------------------------------------------------------
#salmon quantification of TEs

NORM=$(tail -n 1 ${libFOLDER}normalization.txt | tr ' ' '\t' | cut -f 1)

if [[ $TYPE == sRNAseq ]] || [[ $TYPE == sRNAseqIP ]] ; then
  #salmon quantification of TEs
  salmon quant -i ${UTILITY_LOCATION}salmon_quasi-sRNA_${VERSION}${INDEXext} --writeMappings=${TMPdir}mappings.sam --seqBias --gcBias --useVBOpt -l SF -r ${locTMP}decollapsed.fa -o ${locTMP}salmon.quant -p $CORES --incompatPrior 0.0 --validateMappings

  awk -v OFS="\t" -v locTMP=$locTMP -v NORM=$NORM '
  {
    if($1!~"FBgn" && $0!~"Name") {
      if($1 ~ "_AS") {
        if($5>1) {
          print $1,($5/NORM)*1000/$3 > locTMP "TE_antisense.txt"
        }else{
          print $1,0 > locTMP "TE_antisense.txt"
        }
      } else {
        gsub("_AS", "")
        if($5>1) {
          print $1,($5/NORM)*1000/$3 > locTMP "TE_sense.txt"
        }else{
          print $1,0 > locTMP "TE_sense.txt"
        }
      }
    }
  }' <${locTMP}salmon.quant/quant.sf

  sort -k1,1 ${locTMP}TE_sense.txt >${libFOLDER}TE_RPKM_sense-salmon.txt
  sort -k1,1 ${locTMP}TE_antisense.txt >${libFOLDER}TE_RPKM_antisense-salmon.txt

else
  if [[ $FORCEquant_unstranded == Y ]]; then
    UNISTRANDext="_unistrand"
    icncompatPrior=""
    LIBTYPE="U"
  else
    UNISTRANDext=""
    icncompatPrior=" --incompatPrior 0.0 "
    LIBtype="SF"
  fi
  salmon quant -i ${UTILITY_LOCATION}salmon_quasi_${VERSION}${INDEXext}${UNISTRANDext} -g ${UTILITY_LOCATION}transcript_to_gene_${VERSION}.txt --dumpEqWeights --writeMappings=${TMPdir}mappings.sam --seqBias --gcBias --useVBOpt --numBootstraps 100 -l SF -r ${locTMP}decollapsed.fa -o ${locTMP}salmon.quant -p $CORES $icncompatPrior --validateMappings

  if [[ $FORCEquant_unstranded == Y ]]; then
    awk -v OFS="\t" -v locTMP=$locTMP -v NORM=$NORM '
    {
      if($1!~"FBgn" && $0!~"Name") {
          print $1,($5/NORM)*1000/$3 > locTMP "TE.txt"
        }
      }
    }' <${locTMP}salmon.quant/quant.sf

  else
    awk -v OFS="\t" -v locTMP=$locTMP -v NORM=$NORM '
    {
      if($1!~"FBgn" && $0!~"Name") {
        if($1 ~ "_AS") {
          print $1,($5/NORM)*1000/$3 > locTMP "TE_antisense.txt"
        } else {
          gsub("_AS", "")
          print $1,($5/NORM)*1000/$3 > locTMP "TE_sense.txt"
        }
      }
    }' <${locTMP}salmon.quant/quant.sf
  fi

  sort -k1,1 ${locTMP}TE_sense.txt >${libFOLDER}TE_RPKM_sense.txt
  sort -k1,1 ${locTMP}TE_antisense.txt >${libFOLDER}TE_RPKM_antisense.txt

  if [[ $TYPE == CLIPseq ]]; then
    #salmon quantification of TEs using CLIPtags only

    #filter input reads for CLIPtags
    awk -v RS=">" '{
      if(NR>1) {  
        sub("\n", "\t")
        gsub("\n", "")
        if($1~"CLIPtag=TC" || $1~"CLIPtag=Tdel"){
          print ">"$1,"\n"$2
        }
      }
    }' ${locTMP}decollapsed.fa > ${locTMP}decollapsed_CLIP.fa

    salmon quant -i ${UTILITY_LOCATION}salmon_quasi_${VERSION}${INDEXext} -g ${UTILITY_LOCATION}transcript_to_gene_${VERSION}.txt --dumpEqWeights --writeMappings=${TMPdir}mappings_CLIPtag.sam --seqBias --gcBias --useVBOpt --numBootstraps 100 -l SF -r ${locTMP}decollapsed_CLIP.fa -o ${locTMP}salmon_CLIPtag.quant -p $CORES --incompatPrior 0.0 --validateMappings

    awk -v OFS="\t" -v locTMP=$locTMP -v NORM=$NORM '
    {
      if($1!~"FBgn" && $0!~"Name") {
        if($1 ~ "_AS") {
          if($5>1) {
            print $1,($5/NORM)*1000/$3 > locTMP "TE_antisense_CLIPtag.txt"
          }else{
            print $1,0 > locTMP "TE_antisense_CLIPtag.txt"
          }
        } else {
          gsub("_AS", "")
          if($5>1) {
            print $1,($5/NORM)*1000/$3 > locTMP "TE_sense_CLIPtag.txt"
          }else{
            print $1,0 > locTMP "TE_sense_CLIPtag.txt"
          }
        }
      }
    }' <${locTMP}salmon.quant/quant.sf

    sort -k1,1 ${locTMP}TE_sense_CLIPtag.txt >${libFOLDER}TE_RPKM_sense-salmon_CLIPtag.txt
    sort -k1,1 ${locTMP}TE_antisense_CLIPtag.txt >${libFOLDER}TE_RPKM_antisense-salmon_CLIPtag.txt

    mv ${locTMP}salmon_CLIPtag.quant ${TMPdir}
  fi

  if [[ $TYPE == GROseq ]]; then
    salmon quant -i ${UTILITY_LOCATION}salmon_quasi-CDS_${VERSION}${INDEXext} -g ${UTILITY_LOCATION}CDS_to_gene_${VERSION}.txt --dumpEqWeights --seqBias --gcBias --useVBOpt --numBootstraps 100 -l SF -r ${locTMP}decollapsed.fa -o ${locTMP}salmon.quant.CDS -p $CORES --incompatPrior 0.0 --validateMappings
  fi

fi

#exit script if not RNAseq
if [[ $TYPE == RNAseq || $TYPE == GROseq ]]; then
  if [[ $TYPE == GROseq ]]; then nITER=2; else nITER=1; fi
  for ITER in $(seq 1 $nITER); do
    if [[ $ITER -eq 2 ]]; then EXT=".CDS"; else EXT=""; fi
  
    awk -v OFS="\t" '
    {
      if($1~"FBgn") {

        split($1,splitNAME,/::|=/)

        n=split(splitNAME[5],splitID,/:CG/)
        if (n == 1) {
          NAME=splitNAME[7]
          SYN=$1
        }else{
          NAME="FBgnXX_"splitID[1]
          SYN="x::x::x::x::"splitID[1]
        }

        X[NAME]+=$4
        Y[NAME]=SYN
        Z[NAME]+=$5  
        fractLENGTH=$3/$2
        if(maxEffLength[NAME] < fractLENGTH) {
          maxEffLength[NAME]=fractLENGTH
        }

      }
    }
    END{
      for(i in X) {
        split(Y[i],splitNAME,/::|-R[A-Z]|=/)

        if(Z[i]>10){
        #if( maxEffLength[i] > 0.5) {
          print i"::"splitNAME[5],X[i]
        }else{
          if( maxEffLength[i] > 0.5) {
            print i"::"splitNAME[5],X[i]
          }else{
            print i"::"splitNAME[5],0  
          }
        }      
      }
    }' <${locTMP}salmon.quant${EXT}/quant.sf | sort -k1,1 >${libFOLDER}TPM-table_salmon${EXT}.txt

    mkdir ${rawTMP}salmon.quant.all${EXT}
    cp -r ${locTMP}salmon.quant${EXT} ${rawTMP}salmon.quant.all${EXT}/${NAME}

    #if DGE is requested copy samples to the permanent directory
    if [[ $DGE == Y ]]; then
      mkdir -p ${FOLDER}DGE${EXT}/raw/
      cp -r ${locTMP}salmon.quant${EXT} ${FOLDER}DGE${EXT}/raw/${NAME}
    elif [[ $exportSalmon == Y ]]; then
      mkdir -p ${FOLDER}salmon-raw/
      cp -r ${locTMP}salmon.quant${EXT} ${FOLDER}salmon-raw/${NAME}
    fi

    #copy result to TMP for all-sample normalization
    mkdir -p ${rawTMP}DGE${EXT}/raw/
    rm -rf ${TMPdir}salmon.quant${EXT}
    mv ${locTMP}salmon.quant${EXT} ${TMPdir}salmon.quant${EXT}
  
  done

else
  rm -rf ${TMPdir}salmon.quant
  mv ${locTMP}salmon.quant ${TMPdir}
fi

#---------------------------------------------------------------------------------------------------------

if [[ $DEBUG != Y ]]; then
  rm -rf $locTMP
fi

PROCESSED_TIME=$(echo -e "$(date "+%s")" "$TIME" | mawk '{ print ($1-$2)/60 }')
printf "mapping free quantification - processing_time= ${PROCESSED_TIME} \n" >>"${libFOLDER}time-log.txt"

exit
