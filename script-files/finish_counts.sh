#!/usr/bin/env bash
#$ -S /bin/bash

#$ -clear
#$ -e $JOB_NAME.e$JOB_ID.txt
#$ -o $JOB_NAME.o$JOB_ID.txt
#$ -cwd
#$ -q public.q
#$ -pe smp 10

#SBATCH --cpus-per-task=7
#SBATCH --mem=15g
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
locTMP="${TMPdir}TMP_finish/"

mkdir -p "$locTMP/splits/"
chmod 777 "$locTMP"

export TMPDIR=$locTMP
###################################################################################################

#collapse multi mappers to single annotation
#convert filtering input into pipe spaced

nSPLITS=$(ls ${TMPdir}split_intersect | grep txt | wc -w )
RANGE=$(seq 1 $nSPLITS)

parallel --tmpdir $TMPDIR -j $CORES ". ${SCRIPT_DIR}functions; finish_counts {} ${TMPdir}split_intersect/${NAME}_annotated ${locTMP}splits/${NAME} $FILTERING_INPUT $MIN_LENGTH $MAX_LENGTH $TYPE " ::: $RANGE

###################################################################################################
#collect splits into single files

if [[ $SLAM == Y ]]; then
  sort -k1,1 -k2,2n --parallel $CORES -S12G  ${locTMP}splits/${NAME}*_final_annotated.bed | 
    awk -v OFS="\t" -v FILENAME=${libFOLDER}TC_reads.txt '
    BEGIN{
      while((getline LINE < FILENAME) > 0) {
        split(LINE,ID,/ |\t/)
        IDs[ID[1]]=0
      }
    }
    {
      split($4,splitNAME,/:mapping/)
      if(splitNAME[1] in IDs){
        $4=$4":TCslam"
      }
      print
    }' | gzip > ${locTMP}${NAME}_final_annotated.bed.gz

    seqkit fx2tab ${locTMP}splits/${NAME}*_final_annotated.fa.gz |
    awk -v OFS="\t" -v FILENAME=${libFOLDER}TC_reads.txt '
    BEGIN{
      while((getline LINE < FILENAME) > 0) {
        split(LINE,ID,/ |\t/)
        IDs[ID[1]]=0
      }
    }
    {
      split($1,splitNAME,/:mapping/)
      if(splitNAME[1] in IDs){
        $1=$1":TCslam"
      }
      print
    }' | seqkit tab2fx | gzip > ${locTMP}${NAME}_final_annotated.fa.gz

else
  sort -k1,1 -k2,2n --parallel $CORES -S12G  ${locTMP}splits/${NAME}*_final_annotated.bed | gzip > ${locTMP}${NAME}_final_annotated.bed.gz
  cat ${locTMP}splits/${NAME}*_final_annotated.fa.gz > ${locTMP}${NAME}_final_annotated.fa.gz
fi


for FILEtype in final_counts splitup_mRNA splitup_TE size-profile_all size-profile_filtered size-profile_miRNA size-profile_rRNA size-profile_TE
do
  mawk -v OFS="\t" '{
    X[$1]+=$2
  }
  END{
    for(i in X) {
      print i,X[i]
    }
  }' <(cat ${locTMP}splits/${NAME}_*_${FILEtype}.txt) | sort -k1,1 > ${locTMP}${NAME}_${FILEtype}.txt
done
wait

###################################################################################################
#finish output for annotation counts
##annotation output gets sorted according the ruleset_sort
sort -k1,1 "${locTMP}${NAME}_final_counts.txt" > "${locTMP}${NAME}_final_counts_sort.txt"
tail -n +2 "${UTILITY_DIR}ruleset.txt" | grep -E -v "pseudo|ncRNA" | sort -k2,2  | cut -f 2,4 > "${locTMP}ruleset_sort.txt"

join -a 1 -o auto -e 0 -1 1 -2 1 "${locTMP}ruleset_sort.txt" "${locTMP}${NAME}_final_counts_sort.txt" | \
  mawk -v OFS="\t" '{X[$1]=$3; Y[$1]=$2} END{for (i in X) print i,X[i],Y[i] }' | \
  sort -k3,3n | cut -f 1,2 > "${libFOLDER}${NAME}_annotation_counts.txt"

#---------------------------------------------------------------------------------------------------------
##duplication output gets sorted according the ruleset_sort
#quantify duplication within unmapped reads
gunzip -c "${TMPdir}${NAME}_unmapped.fa.gz" | \
  mawk -v OFS="\t" -v locTMP=${locTMP}${NAME} '
    {
      if($1~">") {
        split($1,splitNAME,/:|_|=/)
        TAGnumber+=splitNAME[3]
        ANNdepth+=splitNAME[6]
        READcount+=1
      }
    }
    END{
      print "unmapped",TAGnumber/READcount,READcount > locTMP "_duplication.txt"
      print "unmapped",ANNdepth/READcount,READcount > locTMP "_depth.txt"
   }
  '

for FILEtype in duplication depth
do
  #collect splits
  mawk -v OFS="\t" -v TMP=${locTMP}${NAME} -v FILEtype=$FILEtype '{
    ANNdepth[$1]+=$2
    READcount[$1]+=$3
  }
  END{
    for(ANN in ANNdepth) {
      totANNdepth+=ANNdepth[ANN]
      totREADcount+=READcount[ANN]

      print ANN,ANNdepth[ANN]/READcount[ANN],READcount[ANN]
    }
    print "average_mapped",totANNdepth/totREADcount,totREADcount > TMP"_" FILEtype "_avg.txt"
  }' <(cat ${locTMP}splits/${NAME}_*_${FILEtype}.txt) | sort -k1,1 >> "${locTMP}${NAME}_${FILEtype}.txt"

  #finalize duplication counts for annotations
  sort -k1,1 "${locTMP}${NAME}_${FILEtype}.txt" > "${locTMP}${NAME}_${FILEtype}_sort.txt"
  tail -n +2 "${UTILITY_DIR}ruleset.txt" | grep -E -v "pseudo|ncRNA" | sort -k2,2  | cut -f 2,4 > "${locTMP}ruleset_sort.txt"

  #sort duplication counts according ruleset and add to the unmapped counts
  join -a 1 -o auto -e 0 -1 1 -2 1 "${locTMP}ruleset_sort.txt" "${locTMP}${NAME}_${FILEtype}_sort.txt" | \
    mawk -v OFS="\t" '{X[$1]=$3; Y[$1]=$2; Z[$1]=$4} END{for (i in X) print i,X[i],Y[i],Z[i] }' | \
    sort -k3,3n | cut -f 1,2,4 >> "${locTMP}${NAME}_${FILEtype}-final.txt"

  #add average duplication count
  cat "${locTMP}${NAME}_${FILEtype}_avg.txt" >> "${locTMP}${NAME}_${FILEtype}-final.txt"
done

if [[ $TYPE != sRNAseq ]]
then
  mv "${locTMP}${NAME}_depth-final.txt" "${libFOLDER}${NAME}_sequencing-depth.txt"
else
  mv "${locTMP}${NAME}_depth-final.txt" "${libFOLDER}${NAME}_sequencing-depth.txt"
  mv "${locTMP}${NAME}_duplication-final.txt" "${libFOLDER}${NAME}_read-duplication.txt"
fi

#---------------------------------------------------------------------------------------------------------
#finish output for TE splitup
for ANNclass in TE mRNA
do
  mawk -v OFS="\t" '{ split($1,X,/::/); split(X[1],Y,/~/); split(X[2],Z,/:/); if( Y[2] ~ "AS") print Z[2]"_"Y[2],$2; else print Z[2],$2 }' < ${locTMP}${NAME}_splitup_${ANNclass}.txt | \
    sort -k1,1 > ${libFOLDER}${NAME}_splitup_${ANNclass}.txt
done
#---------------------------------------------------------------------------------------------------------
#finish output for exon splitup
#sort -k1,1 "${locTMP}${NAME}_splitup_mRNA.txt" > "${libFOLDER}${NAME}_splitup_mRNA.txt"

#---------------------------------------------------------------------------------------------------------
#for CLIPseq libraries add CLIPtag annotation to the read-name; for all others only transfer to open directory

if [[ $TYPE != CLIPseq ]]; then
  mv ${locTMP}${NAME}_final_annotated.bed.gz ${TMPdir}${NAME}_annotated.bed.gz
  mv ${locTMP}${NAME}_final_annotated.fa.gz ${TMPdir}${NAME}_annotated.fa.gz
  cat "${TMPdir}${NAME}_unmapped.fa.gz"  >> ${TMPdir}${NAME}_annotated.fa.gz
else
  #process mapped bed file
  gunzip -c ${locTMP}${NAME}_final_annotated.bed.gz |
    awk -v OFS="\t" -v CLIPfile=${TMPdir}CLIPtags.txt '
    BEGIN{
      while((getline LINE < CLIPfile) > 0) {
        split(LINE, splitLINE,/\t| /)
        CLIPtags[splitLINE[1]]["TC"]=splitLINE[2]
        CLIPtags[splitLINE[1]]["Tdel"]=splitLINE[3]
      }
    }
    {
      split($4,splitNAME, /:ann=/)
      if(splitNAME[1] in CLIPtags){
        if(CLIPtags[splitNAME[1]]["TC"]=="Y"){
          $4=$4":CLIPtag=TC"
        }else{
          if(CLIPtags[splitNAME[1]]["Tdel"]=="Y"){
            $4=$4":CLIPtag=Tdel"
          }
        }
      }else{
        $4=$4":CLIPtag=none"
      }
      print
    }' | gzip > ${TMPdir}${NAME}_annotated.bed.gz 

  #process unmapped reads
  cat ${TMPdir}${NAME}_unmapped.fa.gz >>${locTMP}${NAME}_final_annotated.fa.gz

  gunzip -c ${locTMP}${NAME}_final_annotated.fa.gz |
    awk -v OFS="\t" -v CLIPfile=${TMPdir}CLIPtags.txt '
    BEGIN{
      while((getline LINE < CLIPfile) > 0) {
        split(LINE, splitLINE,/\t| /)
        CLIPtags[splitLINE[1]]["TC"]=splitLINE[2]
        CLIPtags[splitLINE[1]]["Tdel"]=splitLINE[3]
      }
      RS=">" 
    }
    {
      if(NR>1){
        sub("\n", "\t")
        gsub("\n", "")
        split($1,splitNAME, /:ann=/)
        if(splitNAME[1] in CLIPtags){
          if(CLIPtags[splitNAME[1]]["TC"]=="Y"){
            $1=$1":CLIPtag=TC"
          }else{
            if(CLIPtags[splitNAME[1]]["Tdel"]=="Y"){
              $1=$1":CLIPtag=Tdel"
            }
          }
        }else{
          $1=$1":CLIPtag=none"
        }
        print ">"$1"\n"$2
      }
    }' | gzip > ${TMPdir}${NAME}_annotated.fa.gz 
fi

#copy files to open-directory
cp ${TMPdir}${NAME}_annotated.fa.gz ${libFOLDER}${NAME}_annotated.fa.gz
cp ${TMPdir}${NAME}_annotated.bed.gz ${libFOLDER}${NAME}_annotated.bed.gz

#---------------------------------------------------------------------------------------------------------
#move size-histogram files to library-folder
mkdir ${libFOLDER}size-profiles/
for CATEGORY in all filtered miRNA rRNA TE; do

  cp ${locTMP}${NAME}_size-profile_${CATEGORY}.txt ${libFOLDER}size-profiles/size-profile_${NAME}-${CATEGORY}.txt
done

#---------------------------------------------------------------------------------------------------------
#cleanup
if [[ $DEBUG != Y ]]; then
  rm -rf $locTMP
fi
PROCESSED_TIME=$(echo -e "$(date "+%s")" "$TIME" | mawk '{ print ($1-$2)/60 }')
printf "collapse annotations for multi-mappers and finish counts - processing_time= ${PROCESSED_TIME} \n" >>"${libFOLDER}time-log.txt"

exit

