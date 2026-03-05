#!/usr/bin/env bash

hostname=
set -u

###################################################################################################
#extract variables 
VARI=$(echo "$1" | sed 's/,/\t/g;s/"//g' )
eval "$VARI"

TIME=$(date "+%s")
GREP=$(which grep)

COLOR=$(echo "$subCOLOR" | tr '~' ',')

#presetting of environment to preload modules and to prepare all files required
#load all functions
source ${SCRIPT_DIR}tools


###################################################################################################

HUB_FOLDER=${FOLDER}/hubs/
FOLDER_NAME_for_hub=$(echo "$FOLDER" | tr '/' ' ' | awk '{ print $NF }')
#@ rm -rf "$HUB_FOLDER"

if [[ $GENOME_VERSION == Cel ]]; then
  GENOME_TAG="ce11"
else
  GENOME_TAG=$GENOME_VERSION
fi

#---------------------------------------------------------------------------------------------------------
for MULTI in uniq all; do
  HUB_FOLDER2="$HUB_FOLDER/${MULTI}/"
  mkdir -p "${HUB_FOLDER2}/${GENOME_TAG}/"
  rm -rf ${HUB_FOLDER2}/${GENOME_TAG}/trackDb.txt


  #create hub.txt
  printf "hub ${FOLDER_NAME_for_hub}_${MULTI}-mappers\nshortLabel ${FOLDER_NAME_for_hub}_${MULTI}\nlongLabel BW_tracks_of_${FOLDER_NAME_for_hub}_${MULTI}-mappers\ngenomesFile genomes.txt\nemail dominik.handler@imba.oeaw.ac.at" >"${HUB_FOLDER2}/hub.txt"
  #create genomes.txt
  printf "genome ${GENOME_TAG}\ntrackDb ${GENOME_TAG}/trackDb.txt" >"${HUB_FOLDER2}/genomes.txt"

  #delete trackdb.txt
  #---------------------------------------------------------------------------------------------------------
  if [[ $TYPE == "CHIPseq" ]] || [[ $TYPE == "DNAseq" ]] || [[ $noSTRANDED == Y ]]; then
    if [[ $RATIOtracks == Y ]]; then
      mkdir -p ${TMPdir}ratio-tracks/
      while read line; do
        CHIPlibs=$(echo $line | tr ' ' '\t' | sed 's/:!!:/\n/g;s/\t/\n/g'| awk -v TMPdir=$TMPdir -v MULTI=$MULTI '{print TMPdir $1 "/bedgraph/" $1 "_" MULTI ".bedgraph"}' | tr '\n' '\t' )
        bedtools unionbedg -i $CHIPlibs | 
        mawk -v OFS="\t"  -v LINE="$line" -v MULTI=$MULTI -v TMP=${TMPdir}ratio-tracks/ '
        BEGIN{
          n=split(LINE, splitLINE,/\t| |:!!:/)
        }
        {
          if($4>0){
            if(NR == 1){
              #@ for(i=2; i<=n; i++){ 
                #@ print  > TMP "ENRICHMENT_"splitLINE[i]":"splitLINE[1]".ratio."MULTI".bedgraph"
                #@ print  > TMP "ENRICHMENT_"splitLINE[i]":"splitLINE[1]".ratio-log."MULTI".bedgraph"
                #@ print  > TMP "ENRICHMENT_"splitLINE[i]"-"splitLINE[1]".subtracted."MULTI".bedgraph"
              #@ }
              a=b
            }else{
              for(i=2; i<=n; i++){ 
                currCOL=3+i
                if($currCOL>0){
                  print $1,$2,$3,$currCOL/$4 >  TMP "ENRICHMENT_"splitLINE[i]":"splitLINE[1]".ratio."MULTI".bedgraph"
                  print $1,$2,$3,log($currCOL/$4) >  TMP "ENRICHMENT_"splitLINE[i]":"splitLINE[1]".ratio-log."MULTI".bedgraph"
                  print $1,$2,$3,$currCOL-$4 > TMP "ENRICHMENT_"splitLINE[i]"-"splitLINE[1]".subtracted."MULTI".bedgraph"
                }
              }
            }
          }
        }' 
      done < ${TMPdir}RATIO.final.txt

      printf "track ENRICHMENT_ratio_${FOLDER_NAME_for_hub}_${MULTI}\ncompositeTrack on\nallButtonPair on\ntype bigWig 0 10\nshortLabel ENRICHMENT_ratio_${FOLDER_NAME_for_hub}_${MULTI}\nlongLabel ENRICHMENT_ratio_${FOLDER_NAME_for_hub}_${MULTI}\nautoScale on\nalwaysZero on \nvisibility full\nmaxHeightPixels 100:50:8\n\n" >>"${HUB_FOLDER2}/${GENOME_TAG}/trackDb.txt"
      printf "track ENRICHMENT_ratio-log_${FOLDER_NAME_for_hub}_${MULTI}\ncompositeTrack on\nallButtonPair on\ntype bigWig 0 10\nshortLabel ENRICHMENT_ratio-log_${FOLDER_NAME_for_hub}_${MULTI}\nlongLabel ENRICHMENT_ratio-log_${FOLDER_NAME_for_hub}_${MULTI}\nautoScale on\nalwaysZero on \nvisibility full\nmaxHeightPixels 100:50:8\n\n" >>"${HUB_FOLDER2}/${GENOME_TAG}/trackDb.txt"
      printf "track ENRICHMENT_subtracted_${FOLDER_NAME_for_hub}_${MULTI}\ncompositeTrack on\nallButtonPair on\ntype bigWig 0 10\nshortLabel ENRICHMENT_subtracted_${FOLDER_NAME_for_hub}_${MULTI}\nlongLabel ENRICHMENT_subtracted_${FOLDER_NAME_for_hub}_${MULTI}\nautoScale on\nalwaysZero on \nvisibility full\nmaxHeightPixels 100:50:8\n\n" >>"${HUB_FOLDER2}/${GENOME_TAG}/trackDb.txt"

      for BEDGRAPH in $(find ${TMPdir}ratio-tracks/ -name ENRICHMENT_*${MULTI}* ); do
        BGfull=$(basename $BEDGRAPH)
        BGname=$(echo ${BGfull}| sed 's/bedgraph/bw/' )
        echo $BGname
        bedGraphToBigWig $BEDGRAPH ${UTILITY_LOCATION}chrom.sizes ${HUB_FOLDER}/${MULTI}/dm6/${BGname}

        BGpure=$(echo "${BGname%.bw}")
        BGtype=$(echo $BGpure | awk '{n=split($1,X,/\./); print X[n-1]}')

        echo $BGpure $BGtype

        if [[ $BGpure == *subtracted* ]]; then
          currCOLOR="130,219,153"
        elif [[ $BGpure == *ratio-log* ]]; then
          currCOLOR="216,145,203"
        else
          currCOLOR="153,191,247"
        fi

        printf "track ${BGpure}\ntype bigWig \nbigDataUrl ${BGname}\nshortLabel ${BGpure}\nlongLabel ${BGpure}\nparent ENRICHMENT_${BGtype}_${FOLDER_NAME_for_hub}_${MULTI}\ncolor ${currCOLOR}\nwindowingFunction mean\n\n" >>"${HUB_FOLDER2}/${GENOME_TAG}/trackDb.txt"

      done

    fi
    
    printf "track ${FOLDER_NAME_for_hub}_${MULTI}\ncompositeTrack on\nallButtonPair on\ntype bigWig 0 10\nshortLabel ${FOLDER_NAME_for_hub}_${MULTI}\nlongLabel ${FOLDER_NAME_for_hub}_${MULTI}\nviewLimits 0:10\nalwaysZero on \nvisibility full\nmaxHeightPixels 100:50:8\n\n" >>"${HUB_FOLDER2}/${GENOME_TAG}/trackDb.txt"

    for i in $(seq 1 "${nFILES}"); do
      name=$(sed -n "${i}"p "$FILE_CONTAINING_LIBRARIES" | awk '{ print $2 }')
      LINE=$(sed -n "${i}"p "$FILE_CONTAINING_LIBRARIES")
      if [[ $LINE == *"COLOR="* ]]; then currCOLOR=$(echo $LINE | awk '{
        for(i=1; i<=NF; i++){ if($i ~ "COLOR=") {split($i,splitCOLOR,/=/); print splitCOLOR[2];exit } } }'); else currCOLOR=$COLOR; fi

      mv "${TMPdir}/${name}/wig-files/${name}_${MULTI}.bw" "${HUB_FOLDER2}/${GENOME_TAG}/${name}_${MULTI}.bw"

      printf "track ${name}_${MULTI}_all_reads\ntype bigWig 0 25\nbigDataUrl ${name}_${MULTI}.bw\nshortLabel ${name}_${MULTI}\nlongLabel ${name}_${MULTI}\nparent ${FOLDER_NAME_for_hub}_${MULTI}\ncolor ${currCOLOR}\nwindowingFunction mean\n\n" >>"${HUB_FOLDER2}/${GENOME_TAG}/trackDb.txt"
    done

    #---------------------------------------------------------------------------------------------------------
  else
    #---------------------------------------------------------------------------------------------------------
    for FILE_NUMBER in $(seq 1 "${nFILES}"); do

      name=$(sed -n "${FILE_NUMBER}"p "${TMPdir}files.txt" | awk '{ print $2 }')
      LINE=$(sed -n "${FILE_NUMBER}"p "$FILE_CONTAINING_LIBRARIES")
      IPstatus=$(echo $LINE | mawk -v line=$FILE_NUMBER '{ for(i=3; i<=NF; i++) {if($i == "IP") { print "IP" }}}')
      if [[ $LINE == *"COLOR="* ]]; then currCOLOR=$(echo $LINE | awk '{
        for(i=1; i<=NF; i++){ if($i ~ "COLOR=") {split($i,splitCOLOR,/=/); print splitCOLOR[2];exit } } }'); else currCOLOR=$COLOR; fi

      #---------------------------------------------------------------------------------------------------------
      n_track_types="1"
      if [[ "$TYPE" == sRNAseq ]] && [[ "$IPstatus" != IP ]]; then n_track_types="1 2 3"; fi
      if [[ "$TYPE" == CapSeq ]]; then n_track_types="1 4"; fi
      if [[ "$TYPE" == RNAseq ]] || [[ $TYPE == "RIPseq" ]] ; then n_track_types="1 5"; fi

      for track_TYPE_nr in $(echo "${n_track_types}" | tr ' ' '\t'); do

        if [[ $track_TYPE_nr == 1 ]]; then
          extension=""
          FOLDER_extension="no_size_restriction"
          VISIBILITY="show"
          VIEW_LIMITS=$(if [[ ${autoViewLimits} == Y ]]; then echo "autoScale on"; else echo "viewLimits -50:50"; fi)
        elif [[ $track_TYPE_nr == 2 ]]; then
          extension="_siRNA"
          FOLDER_extension="siRNA"
          VISIBILITY=""
          VIEW_LIMITS=$(if [[ ${autoViewLimits} == Y ]]; then echo "autoScale on"; else echo "viewLimits -50:50"; fi)
        elif [[ $track_TYPE_nr == 3 ]]; then
          FOLDER_extension="piRNA"
          extension="_piRNA"
          VISIBILITY=""
          VIEW_LIMITS=$(if [[ ${autoViewLimits} == Y ]]; then echo "autoScale on"; else echo "viewLimits -50:50"; fi)
        elif [[ $track_TYPE_nr == 4 ]]; then
          FOLDER_extension="5ends"
          extension="_5ends"
          VISIBILITY="show"
          VIEW_LIMITS="autoScale on"
        elif [[ $track_TYPE_nr == 5 ]]; then
          FOLDER_extension="spliced_reads"
          extension=""
          VISIBILITY=""
        fi

        #preset trackDb file
        if [[ $track_TYPE_nr -le 4 ]]; then
          TEST=$(grep -w ${FOLDER_NAME_for_hub}${extension}_${MULTI} ${HUB_FOLDER2}/${GENOME_TAG}/trackDb.txt)
          if [[ -z $TEST ]]; then
            #preset trackDb.txt
            printf "track ${FOLDER_NAME_for_hub}${extension}_${MULTI}\nsuperTrack on ${VISIBILITY}\ngroup regulation\nshortLabel ${FOLDER_NAME_for_hub}${extension}_${MULTI}\nlongLabel BW_tracks_of_${FOLDER_NAME_for_hub}${extension}_${MULTI}\n\n" >>"${HUB_FOLDER2}/${GENOME_TAG}/trackDb.txt"
            if [[ $SLAM == Y ]]; then
              printf "track ${FOLDER_NAME_for_hub}${extension}_${MULTI}-TC\nsuperTrack on ${VISIBILITY}\ngroup regulation\nshortLabel ${FOLDER_NAME_for_hub}${extension}_${MULTI}-TC\nlongLabel BW_tracks_of_${FOLDER_NAME_for_hub}${extension}_${MULTI}-TC\n\n" >>"${HUB_FOLDER2}/${GENOME_TAG}/trackDb.txt"
            fi           
          fi
        elif [[ $track_TYPE_nr == 5 ]]; then
          TEST=$(grep -w ${FOLDER_NAME_for_hub}${extension}_spliced-reads_${MULTI} ${HUB_FOLDER2}/${GENOME_TAG}/trackDb.txt)
          if [[ -z $TEST ]]; then
            printf "track ${FOLDER_NAME_for_hub}${extension}_spliced-reads_${MULTI}\nsuperTrack on ${VISIBILITY}\ngroup regulation\nshortLabel ${FOLDER_NAME_for_hub}${extension}_${MULTI}\nlongLabel spliced-reads_from_${FOLDER_NAME_for_hub}${extension}_${MULTI}\n\n" >>"${HUB_FOLDER2}/${GENOME_TAG}/trackDb.txt"
            printf "track ${FOLDER_NAME_for_hub}${extension}_splicejunctions_${MULTI}\nsuperTrack on ${VISIBILITY}\ngroup regulation\nshortLabel ${FOLDER_NAME_for_hub}${extension}_${MULTI}\nlongLabel splice-junctions_from_${FOLDER_NAME_for_hub}${extension}_${MULTI}\n\n" >>"${HUB_FOLDER2}/${GENOME_TAG}/trackDb.txt"
          fi
        fi

        if [[ $track_TYPE_nr -le 4 ]]; then
          cp "${TMPdir}/${name}/wig-files/${name}_${MULTI}${extension}_sense.bw" "${HUB_FOLDER2}/${GENOME_TAG}/${name}${extension}_${MULTI}_sense.bw"
          cp "${TMPdir}/${name}/wig-files/${name}_${MULTI}${extension}_antisense.bw" "${HUB_FOLDER2}/${GENOME_TAG}/${name}${extension}_${MULTI}_antisense.bw"
          if [[ $SLAM == Y ]]; then
            cp "${TMPdir}/${name}/wig-files/${name}_${MULTI}-TC${extension}_sense.bw" "${HUB_FOLDER2}/${GENOME_TAG}/${name}${extension}_${MULTI}-TC_sense.bw"
            cp "${TMPdir}/${name}/wig-files/${name}_${MULTI}-TC${extension}_antisense.bw" "${HUB_FOLDER2}/${GENOME_TAG}/${name}${extension}_${MULTI}-TC_antisense.bw"
          fi

          printf "track ${name}${extension}_${MULTI}\ncontainer multiWig\naggregate transparentOverlay\nshowSubtrackColorOnUi on\nshortLabel ${name}${extension}_${MULTI}\n${VIEW_LIMITS}\nalwaysZero on\nlongLabel ${name}${extension}_${MULTI}\ntype bigWig\nparent ${FOLDER_NAME_for_hub}${extension}_${MULTI}\nvisibility full\nmaxHeightPixels 100:50:8\n\n" >>"${HUB_FOLDER2}/${GENOME_TAG}/trackDb.txt"
          printf "track ${name}${extension}_${MULTI}_+\ntype bigWig\nbigDataUrl ${name}${extension}_${MULTI}_sense.bw\nshortLabel ${name}${extension}_${MULTI}_+\nlongLabel ${name}${extension}_${MULTI}_+\nparent ${name}${extension}_${MULTI}\ncolor ${currCOLOR}\nwindowingFunction mean\n\n" >>"${HUB_FOLDER2}/${GENOME_TAG}/trackDb.txt"
          printf "track ${name}${extension}_${MULTI}_-\ntype bigWig\nbigDataUrl ${name}${extension}_${MULTI}_antisense.bw\nshortLabel ${name}${extension}_${MULTI}_-\nlongLabel ${name}${extension}_${MULTI}_-\nparent ${name}${extension}_${MULTI}\ncolor ${currCOLOR}\nwindowingFunction mean\n\n" >>"${HUB_FOLDER2}/${GENOME_TAG}/trackDb.txt"
          if [[ $SLAM == Y ]]; then
            printf "track ${name}${extension}_${MULTI}-TC\ncontainer multiWig\naggregate transparentOverlay\nshowSubtrackColorOnUi on\nshortLabel ${name}${extension}_${MULTI}-TC\n${VIEW_LIMITS}\nalwaysZero on\nlongLabel ${name}${extension}_${MULTI}-TC\ntype bigWig\nparent ${FOLDER_NAME_for_hub}${extension}_${MULTI}-TC\nvisibility full\nmaxHeightPixels 100:50:8\n\n" >>"${HUB_FOLDER2}/${GENOME_TAG}/trackDb.txt"
            printf "track ${name}${extension}_${MULTI}-TC_+\ntype bigWig\nbigDataUrl ${name}${extension}_${MULTI}-TC_sense.bw\nshortLabel ${name}${extension}_${MULTI}-TC_+\nlongLabel ${name}${extension}_${MULTI}-TC_ +\nparent ${name}${extension}_${MULTI}-TC\ncolor ${currCOLOR}\nwindowingFunction mean\n\n" >>"${HUB_FOLDER2}/${GENOME_TAG}/trackDb.txt"
            printf "track ${name}${extension}_${MULTI}-TC_-\ntype bigWig\nbigDataUrl ${name}${extension}_${MULTI}-TC_antisense.bw\nshortLabel ${name}${extension}_${MULTI}-TC_-\nlongLabel ${name}${extension}_${MULTI}-TC_-\nparent ${name}${extension}_${MULTI}-TC\ncolor ${currCOLOR}\nwindowingFunction mean\n\n" >>"${HUB_FOLDER2}/${GENOME_TAG}/trackDb.txt"
          fi
        elif [[ $track_TYPE_nr == 5 ]]; then
          mv "${TMPdir}/${name}/wig-files/${name}_spliced_${MULTI}.bam" "${HUB_FOLDER2}/${GENOME_TAG}/${name}_spliced_${MULTI}.bam"
          mv "${TMPdir}/${name}/wig-files/${name}_spliced_${MULTI}.bam.bai" "${HUB_FOLDER2}/${GENOME_TAG}/${name}_spliced_${MULTI}.bam.bai"
          mv "${TMPdir}/${name}/wig-files/${name}_junctions_${MULTI}.bam" "${HUB_FOLDER2}/${GENOME_TAG}/${name}_junctions_${MULTI}.bam"
          mv "${TMPdir}/${name}/wig-files/${name}_junctions_${MULTI}.bam.bai" "${HUB_FOLDER2}/${GENOME_TAG}/${name}_junctions_${MULTI}.bam.bai"

          printf "track ${name}_spliced-reads\ntype bam\nbamColorMode strand\nshowNames on\nmaxWindowToDraw 500000\nvisibility squish\nbigDataUrl ${name}_spliced_${MULTI}.bam\nshortLabel ${name}_spliced-reads_${MULTI}\nlongLabel ${name}_spliced-reads_${MULTI}\nparent ${FOLDER_NAME_for_hub}_spliced-reads_${MULTI}\n\n" >>"${HUB_FOLDER2}/${GENOME_TAG}/trackDb.txt"
          printf "track ${name}_splice-junctions\ntype bam\nbamColorMode strand\nshowNames on\nmaxWindowToDraw 500000\nvisibility squish\nbigDataUrl ${name}_junctions_${MULTI}.bam\nshortLabel ${name}_splice-junctions_${MULTI}\nlongLabel ${name}_splicejunctions_${MULTI}\nparent ${FOLDER_NAME_for_hub}_splicejunctions_${MULTI}\n\n" >>"${HUB_FOLDER2}/${GENOME_TAG}/trackDb.txt"

        fi

      done
    done
  fi
done
