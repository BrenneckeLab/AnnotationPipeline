#!/usr/bin/env bash

hostname
set -ux

###################################################################################################
#extract variables
VARI=$(echo "$1" | sed 's/,/\t/g;s/"//g')
eval "$VARI"

TIME=$(date "+%s")
GREP=$(which grep)

COLOR=$(echo "$subCOLOR" | tr '~' ',')

###################################################################################################
#wait if GENOMEdir is busy

#do not continue until file is unblocked by other process
while [[ -f ${GENOMEdir}wait.txt ]]; do
  sleep 10s
done

#block GENOMEdir from other processes
touch ${GENOMEdir}wait.txt

#create unique name for tracks
TRACKname=${TYPE}_${RUNname}

###################################################################################################
#remove entries from all files

awk -v FS="\n" -v RS="\n\n" -v OFS="\t" -v ORS="\n\n" -v NAME=$TRACKname '
  {
    if( $0 !~ NAME ) print
  }' ${GENOMEdir}/${VERSION}/trackDb.txt >${TMPdir}trackDb.tmp
mv ${TMPdir}trackDb.tmp ${GENOMEdir}/${VERSION}/trackDb.txt

HUB_FOLDER="${GENOMEdir}/${VERSION}/AP-tracks/${TYPE}/${RUNname}/"
rm -rf $HUB_FOLDER
mkdir -p $HUB_FOLDER

#---------------------------------------------------------------------------------------------------------
for MULTI in uniq all; do

  if [[ $MULTI == all ]]; then PRIORITY=50; else PRIORITY=10; fi
  #---------------------------------------------------------------------------------------------------------
  if [[ $TYPE == "CHIPseq" ]] || [[ $TYPE == "DNAseq" ]] || [[ $noSTRANDED == Y ]]; then

    printf "track ${TRACKname}_${MULTI}\ncompositeTrack on\nallButtonPair on\ntype bigWig 0 10\nshortLabel ${RUNname}_${MULTI}\nlongLabel ${RUNname}_${MULTI}\npriority $PRIORITY \nviewLimits 0:10\nalwaysZero on\nvisibility full\nmaxHeightPixels 100:50:8\n\n" >>${GENOMEdir}/${VERSION}/trackDb.txt

    for i in $(seq 1 "${nFILES}"); do
      name=$(sed -n "${i}"p "$FILE_CONTAINING_LIBRARIES" | awk '{ print $2 }')
      LINE=$(sed -n "${i}"p "$FILE_CONTAINING_LIBRARIES")
      if [[ $LINE == *"COLOR="* ]]; then currCOLOR=$(echo $LINE | awk '{
        for(i=1; i<=NF; i++){ if($i ~ "COLOR=") {split($i,splitCOLOR,/=/); print splitCOLOR[2];exit } } }'); else currCOLOR=$COLOR; fi

      cp "${TMPdir}/${name}/wig-files/${name}_${MULTI}.bw" "${HUB_FOLDER}/${name}_${MULTI}.bw"

      printf "track ${TRACKname}_${MULTI}_${name}\ntype bigWig 0 25\nbigDataUrl AP-tracks/${TYPE}/${RUNname}/${name}_${MULTI}.bw\nshortLabel ${name}_${MULTI}\nlongLabel ${name}_${MULTI}\nparent ${TRACKname}_${MULTI}\ncolor ${currCOLOR}\nwindowingFunction mean\n\n" >>${GENOMEdir}/${VERSION}/trackDb.txt
    done

    #---------------------------------------------------------------------------------------------------------
  else
    #create supertrack for track #1
    printf "track ${TRACKname}_${MULTI}\nsuperTrack on show\nshortLabel ${RUNname}_${MULTI}\nlongLabel ${RUNname}_${MULTI}\npriority ${PRIORITY}\ngroup AP_${TYPE}\n\n" >>${GENOMEdir}/${VERSION}/trackDb.txt
    
    #preset variables to determine superTrack status for other types
    n2=N
    n3=N
    n4=N
    n5=N
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
      if [[ "$TYPE" == RNAseq ]] || [[ $TYPE == "RIPseq" ]]; then n_track_types="1 5"; fi

      for track_TYPE_nr in $(echo "${n_track_types}" | tr ' ' '\t'); do

        if [[ $track_TYPE_nr == 1 ]]; then
          extension=""
          FOLDER_extension="no_size_restriction"
          VISIBILITY="full"
          VIEW_LIMITS=$(if [[ ${autoViewLimits} == Y ]]; then echo "autoScale on"; else echo "viewLimits -50:50"; fi)
        elif [[ $track_TYPE_nr == 2 ]]; then
          extension="_siRNA"
          FOLDER_extension="siRNA"
          VISIBILITY="hide"
          VIEW_LIMITS=$(if [[ ${autoViewLimits} == Y ]]; then echo "autoScale on"; else echo "viewLimits -50:50"; fi)
          if [[ $n2 == N ]]; then 
            printf "track ${TRACKname}_${MULTI}${extension}\nsuperTrack on show\nshortLabel ${RUNname}_${MULTI}${extension}\nlongLabel ${RUNname}_${MULTI}${extension}\npriority ${PRIORITY}\ngroup AP_${TYPE}${extension}\n\n" >>${GENOMEdir}/${VERSION}/trackDb.txt
            n2=Y
          fi
        elif [[ $track_TYPE_nr == 3 ]]; then
          FOLDER_extension="piRNA"
          extension="_piRNA"
          VISIBILITY="hide"
          VIEW_LIMITS=$(if [[ ${autoViewLimits} == Y ]]; then echo "autoScale on"; else echo "viewLimits -50:50"; fi)
          if [[ $n3 == N ]]; then 
            printf "track ${TRACKname}_${MULTI}${extension}\nsuperTrack on show\nshortLabel ${RUNname}_${MULTI}${extension}\nlongLabel ${RUNname}_${MULTI}${extension}\npriority ${PRIORITY}\ngroup AP_${TYPE}${extension}\n\n" >>${GENOMEdir}/${VERSION}/trackDb.txt
            n3=Y
          fi
        elif [[ $track_TYPE_nr == 4 ]]; then
          FOLDER_extension="5ends"
          extension="_5ends"
          VISIBILITY="hide"
          VIEW_LIMITS="autoScale on"
          if [[ $n4 == N ]]; then 
            printf "track ${TRACKname}_${MULTI}${extension}\nsuperTrack on show\nshortLabel ${RUNname}_${MULTI}${extension}\nlongLabel ${RUNname}_${MULTI}${extension}\npriority ${PRIORITY}\ngroup AP_${TYPE}${extension}\n\n" >>${GENOMEdir}/${VERSION}/trackDb.txt
            n4=Y
          fi
        elif [[ $track_TYPE_nr == 5 ]]; then
          FOLDER_extension="spliced_reads"
          extension="_spliceJunctions"
          VISIBILITY="hide"
          if [[ $n5 == N ]]; then 
            printf "track ${TRACKname}_${MULTI}${extension}\nsuperTrack on show\nshortLabel ${RUNname}_${MULTI}${extension}\nlongLabel ${RUNname}_${MULTI}${extension}\npriority ${PRIORITY}\ngroup AP_${TYPE}${extension}\n\n" >>${GENOMEdir}/${VERSION}/trackDb.txt
            n5=Y
          fi
        fi

        #hide all multi-mapper tracks by default
        if [[ $MULTI == all ]]; then
          VISIBILITY=hide
        fi

        if [[ $track_TYPE_nr -le 4 ]]; then
          cp "${TMPdir}/${name}/wig-files/${name}_${MULTI}${extension}_sense.bw" "${HUB_FOLDER}/${name}_${MULTI}_sense.bw"
          cp "${TMPdir}/${name}/wig-files/${name}_${MULTI}${extension}_antisense.bw" "${HUB_FOLDER}/${name}_${MULTI}_antisense.bw"
          if [[ $SLAM == Y ]]; then
            cp "${TMPdir}/${name}/wig-files/${name}_${MULTI}-TC${extension}_sense.bw" "${HUB_FOLDER}/${name}_${MULTI}-TC_sense.bw"
            cp "${TMPdir}/${name}/wig-files/${name}_${MULTI}-TC${extension}_antisense.bw" "${HUB_FOLDER}/${name}_${MULTI}-TC_antisense.bw"
          fi

          printf "track ${TRACKname}_${MULTI}_${name}${extension}\ncontainer multiWig\naggregate transparentOverlay\nshowSubtrackColorOnUi on\nshortLabel ${name}${extension}_${MULTI}\n${VIEW_LIMITS}\nalwaysZero on\nlongLabel ${name}${extension}_${MULTI}\ntype bigWig\nparent ${TRACKname}_${MULTI}${extension}\nvisibility ${VISIBILITY}\nmaxHeightPixels 100:50:8\n\n" >>${GENOMEdir}/${VERSION}/trackDb.txt
          printf "track ${TRACKname}_${MULTI}_${name}${extension}_+\ntype bigWig\nbigDataUrl AP-tracks/${TYPE}/${RUNname}/${name}_${MULTI}_sense.bw\nshortLabel ${name}${extension}_${MULTI}_+\nlongLabel ${name}${extension}_${MULTI}_+\nparent ${TRACKname}_${MULTI}_${name}${extension}\ncolor ${currCOLOR}\nwindowingFunction mean\n\n" >>${GENOMEdir}/${VERSION}/trackDb.txt
          printf "track ${TRACKname}_${MULTI}_${name}${extension}_-\ntype bigWig\nbigDataUrl AP-tracks/${TYPE}/${RUNname}/${name}_${MULTI}_antisense.bw\nshortLabel ${name}${extension}_${MULTI}_-\nlongLabel ${name}${extension}_${MULTI}_-\nparent ${TRACKname}_${MULTI}_${name}${extension}\ncolor ${currCOLOR}\nwindowingFunction mean\n\n" >>${GENOMEdir}/${VERSION}/trackDb.txt
          if [[ $SLAM == Y ]]; then
            TCcolor=$(echo $currCOLOR | tr ',' '\t' | awk '{for(i=1;i<=3;i++){X=$i*1.4; if(X>255){Y[i]=255;}else{Y[i]=X}}}END{ print sprintf("%3.0f", Y[1]) ,sprintf("%3.0f", Y[2]),sprintf("%3.0f", Y[3])}')
            printf "track ${TRACKname}_${MULTI}_${name}${extension}-TC\ncontainer multiWig\naggregate transparentOverlay\nshowSubtrackColorOnUi on\nshortLabel ${name}${extension}_${MULTI}-TC\n${VIEW_LIMITS}\nalwaysZero on\nlongLabel ${name}${extension}_${MULTI}-TC\ntype bigWig\nparent ${TRACKname}_${MULTI}\nvisibility ${VISIBILITY}\nmaxHeightPixels 100:50:8\n\n" >>${GENOMEdir}/${VERSION}/trackDb.txt
            printf "track ${TRACKname}_${MULTI}_${name}${extension}-TC_+\ntype bigWig\nbigDataUrl AP-tracks/${TYPE}/${RUNname}/${name}_${MULTI}-TC_sense.bw\nshortLabel ${name}${extension}_${MULTI}-TC_+\nlongLabel ${name}${extension}_${MULTI}-TC_ +\nparent ${TRACKname}_${MULTI}_${name}${extension}-TC\ncolor ${TCcolor}\nwindowingFunction mean\n\n" >>${GENOMEdir}/${VERSION}/trackDb.txt
            printf "track ${TRACKname}_${MULTI}_${name}${extension}-TC_-\ntype bigWig\nbigDataUrl AP-tracks/${TYPE}/${RUNname}/${name}_${MULTI}-TC_antisense.bw\nshortLabel ${name}${extension}_${MULTI}-TC_-\nlongLabel ${name}${extension}_${MULTI}-TC_-\nparent ${TRACKname}_${MULTI}_${name}${extension}-TC\ncolor ${TCcolor}\nwindowingFunction mean\n\n" >>${GENOMEdir}/${VERSION}/trackDb.txt
          fi
        elif [[ $track_TYPE_nr == 5 ]]; then
          mv "${TMPdir}/${name}/wig-files/${name}_spliced_${MULTI}.bam" "${HUB_FOLDER}/${name}_spliced_${MULTI}.bam"
          mv "${TMPdir}/${name}/wig-files/${name}_spliced_${MULTI}.bam.bai" "${HUB_FOLDER}/${name}_spliced_${MULTI}.bam.bai"
          #!@ mv "${TMPdir}/${name}/wig-files/${name}_junctions_${MULTI}.bam" "${HUB_FOLDER}/${name}_junctions_${MULTI}.bam"
          #!@ mv "${TMPdir}/${name}/wig-files/${name}_junctions_${MULTI}.bam.bai" "${HUB_FOLDER}/${name}_junctions_${MULTI}.bam.bai"

          printf "track ${TRACKname}_${MULTI}_${name}_spliced-reads\ntype bam\nbamColorMode strand\nshowNames on\nmaxWindowToDraw 500000\nvisibility squish\nbigDataUrl AP-tracks/${TYPE}/${RUNname}/${name}_spliced_${MULTI}.bam\nshortLabel ${name}_spliced-reads_${MULTI}\nlongLabel ${name}_spliced-reads_${MULTI}\nparent ${TRACKname}_${MULTI}${extension}\n\n" >>${GENOMEdir}/${VERSION}/trackDb.txt
          #!splice-junctions not created for asm
          #!@ printf "track ${TRACKname}_${MULTI}_${name}_splice-junctions\ntype bam\nbamColorMode strand\nshowNames on\nmaxWindowToDraw 500000\nvisibility squish\nbigDataUrl AP-tracks/${TYPE}/${RUNname}/${name}_junctions_${MULTI}.bam\nshortLabel ${name}_splice-junctions_${MULTI}\nlongLabel ${name}_splicejunctions_${MULTI}\nparent ${TRACKname}_${MULTI}${extension}\n\n" >>${GENOMEdir}/${VERSION}/trackDb.txt

        fi

      done
    done
  fi
done

#unblock directory
rm -rf ${GENOMEdir}wait.txt
