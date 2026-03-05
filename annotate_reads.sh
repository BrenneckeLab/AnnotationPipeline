#!/usr/bin/env bash

set -ux
#for coding
#GENOME_VERSION="dm3"
#ANNOTATION="r6.07"

#set IMPIMBA_MACHINE_NAME on old cluster
SYSTEM="EXTERN"
GRIDsystem="SLURM"

#test
#---------------------------------------------------------------------------------------------------------

#hardcoded variables
PIPELINEversion="4.0"

default_VERSION_dm3="XX-dm3-VERSION-XX"
default_VERSION_dm6="XX-dm6-VERSION-XX"
default_VERSION_Cel="XX-Cel-VERSION-XX"
BASE_UTILITY_LOCATION="XX-BASE_UTILITY_LOCATION-XX"
BASE_HTTP="XX-HTTP_LINK-XX"
LIB_STORAGE_FOLDER="XX-BASE_LIBRARY_STORAGE_LOCATION-XX"
mkdir -p $LIB_STORAGE_FOLDER

BASE_FOLDER="XX-BASE_RESULTS_FOLDER-XX"
SERVER="SERVER"

BASE_FOLDER_TMP="XX-BASE_FOLDER_TMP-XX"
TMPdir=$BASE_FOLDER_TMP
downDIR=${BASE_FOLDER_TMP}/NGS_downloads/
mkdir -p $downDIR

usage() {
  cat <<EOF
  usage: $0 options
  
  ###############################################################################
  This scripts analyses deep-sequencing libraries and generates
  several statistics and UCSC compatiple tracks.
  It can also convert raw-bam files located at NGS into adaptor clipped and N trimmed
  fasta files. 
  
  Version= ${PIPELINEversion}
  
  Do not start the script in the background by adding "&" to the start command.

  
  to attach the screen again use the command screen -r
  
  A detailed description of the functionality and usage can be found in:
  "XX-SCRIPT_DIR-XX"README.md
      
  OPTIONS:
      -h, --help                Show this message
    
  OBLIGATORY [but only one of the two options at any time]
      -i, --input FILE          File containing the information for all the libraries
      -u, --update PATH         Update old AnnotationPipeline run to newest version [path to folder]

  OPTIONAL - available in interactive mode
      -F, --folder-name NAME    Name attachment for the output folder
      -t, --type TYPE           TYPE of LIBRARY
      -j, --slam                Libraries contain SLAMseq T->C conversions
      -v, --genome-version VER  Genome version (dm3, dm6, ASM...)
      -V, --version VER         Annotation version
      -D, --dge                 Perform DGE analysis (RNAseq only)

    bam -> fasta 
      --demux-only              Only demultiplex NGS data without processing individual libraries
      -x, --force-preprocess    Force preprocessing even if input not bam
      -N, --n-trimm NUM         Trim N random nucleotides from each end
      -P, --raw-paired          Also output raw paired file
      -p, --only-paired         Only paired end file (implies -P and -Q)
      -Q, --fastq-out           Output trimmed fastq
      -q, --fastq-out-raw       Output raw demultiplexed fastq
      -A, --custom-adaptors     Use custom adaptor sequences
      -r, --raw                 Only generate demultiplexed preprocessed files
      -R, --subsample NUM       Subsample fraction of reads

    read-preprocessing        
      -m, --min-length NUM      Minimal read length after trimming
      -M, --max-length NUM      Maximal read length after trimming
      -T, --trimm               Enable trimming
      -f, --first NUM           First base to keep
      -l, --last NUM            Last base to keep 
      --polya               Remove polyA stretches
      -I, --invert              Invert reads
      -S, --single-end          Force single-end processing
      -2, --second-mate         Only process 2nd mate
      -J, --max-count NUM       Max sequence count to consider
      -y, --demux-fasta         Enable sRBC demultiplexing

    mapping & analysis 
      -Y, --y-chrom             Include Y chromosome
      -s, --mismatches NUM      Allowed mismatches
      -E, --random-multi        Distribute multimappers evenly
      -W, --wig-norm            Normalize wig to 10M uniquely mapped
      -w, --wig-fasta-norm      Normalize wig to 10M input reads
      -L, --spike-in-norm       Spike-in based normalization
      -z, --no-norm             Disable normalization
      --extra-seq FILE          fasta-File containing extra sequences to be added to the TE-histogram analysis
      --force-quant-unstranded  Force quantification to GeTMM for unstranded libraries (e.g. ChIPseq)
      -5, --only-5end           Wig from 5' end only
      -1, --ping-pong           Ping‑pong/phasing analysis (sRNAseq)
      -n, --no-stranded         Unstranded analysis
      -e, --extend NUM          Extend reads to fragment length
      -b, --export-bam          Export collapsed bam
      -B, --export-bam-uncollapsed Export uncollapsed bam
      -3, --export-salmon       Export salmon raw files
      -c, --color RGB           Track color (R,G,B)
      -a, --auto-scale          Set flag to set wig tracks to auto scale in UCSC

    user management & data import
      -G, --geo                 Create GEO export
      -U, --user NAME           Run as different user


            
EOF
}

    #@ misc / internal
    #@   -C, --computing           Run in parallel/cluster mode
    #@   -k, --keep-tmp            Keep temporary files
    #@   -O, --folder PATH         Base output folder
    #@   -o, --tmp PATH            Temporary directory
    #@   -H, --http URL            HTTP base URL
    #@   -K, --update-mode         Update mode
    #@   -0, --custom-outdir PATH  Custom output directory
    #@   -d, --debug               Debug mode

# ---- SHORT and LONG option specs (short unchanged, long just added) ----
SHORT_OPTS="hi:u:F:t:jv:V:DxN:PpQqArR:m:M:Tf:l:IS2J:Ys:EWywLzna51e:bB3c:XGU:CkO:o:H:KdZg0:"

LONG_OPTS="help,input:,update:,folder-name:,type:,slam,genome-version:,version:,dge,force-preprocess,n-trimm:,raw-paired,only-paired,fastq-out,fastq-out-raw,custom-adaptors,raw,subsample:,min-length:,max-length:,trimm,first:,last:,polya,invert,single-end,second-mate,max-count:,demux-fasta,demux-only,y-chrom,mismatches:,random-multi,wig-norm,wig-fasta-norm,spike-in-norm,no-norm,only-5end,ping-pong,no-stranded,extend:,export-bam,export-bam-uncollapsed,extra-seq:,export-salmon,color:,auto-scale,force-import,geo,user:,computing,keep-tmp,folder:,tmp:,http:,update-mode,alt-server,ameres,custom-outdir:,force-quant-unstranded,debug"
# ---- Call GNU getopt ----
OPTS=$(getopt -o "$SHORT_OPTS" --long "$LONG_OPTS" -n 'parse-options' -- "$@")
if [ $? != 0 ]; then 
  echo "Failed parsing options." >&2
  usage
  exit 1
fi

eval set -- "$OPTS"

# ---- Initialize variables (unchanged from your script) ----
FILE_CONTAINING_LIBRARIESraw=
UPDATErun=
FOLDER_NAME=
TYPE=
SLAM=N
GENOME_VERSION=
VERSION=
DGE=N
DEMUXonly=N
rawPAIRED=N
onlyPAIRED=N
FASTQout=N
FASTQoutRAW=N
customADAPTORS=
fwADAPTOR=
rvADAPTOR=
RAW=
FORCE=
N_TRIMM=
SUBSAMPLE=
MIN_LENGTH=
MAX_LENGTH=
TRIMM=
FIRST=1
LAST=1000
modTRIMM=N
POLYA=N
INVERT=N
SE=N
SE2nd=
maxCOUNT=1000000
demuxFASTA=N
Ychrom=N
MMinput=
RANDOMmulti=N
WIG=
WIG_FASTA=
spikeINnorm=N
noNORM=
autoViewLimits=
only5end=
PingPong=
noSTRANDED=
EXTEND=0
extraSEQ=
exportBAM=N
exportBAMuncollapsed=N
exportSalmon=N
FORCEimport=
GEO=
USER=
altSERVER=N
AMERES=N
FORCEquant_unstranded=N


COMPUTING="C"
keepTMP="N"
FOLDER=
TMP=
HTTP=
UPDATE=
nSPLITS=1000000
ngsUSER=
ngsPW=
DEBUG=N
custom_OUTdir=

RERUN_OPTS=""

# ---- Parse loop: short and long mapped together ----
while true; do
  case "$1" in
    # no large/small here originally, but we want it in the rerun command
    -i | --input )
      FILE_CONTAINING_LIBRARIESraw=$2
      # or $2 if you’re not using getopts
      shift 2
      ;;

    -F | --folder-name )
      FOLDER_NAME=$2
      shift 2
      ;;

    -t | --type )
      TYPE=$2
      shift 2
      ;;

    -v | --genome-version )
      GENOME_VERSION=$2
      shift 2
      ;;

    -V | --version )
      VERSION=$2
      shift 2
      ;;

    -N | --n-trimm )
      re='^[0-9]+$'
      if ! [[ $2 =~ $re ]] ; then
        echo "Option $2 requires a numeric argument."; exit 1
      fi
      N_TRIMM=$2
      shift 2
      ;;

    -m | --min-length )
      re='^[0-9]+$'
      if ! [[ $2 =~ $re ]] ; then
        echo "Option $2 requires a numeric argument."; exit 1
      fi
      MIN_LENGTH=$2
      shift 2
      ;;

    -M | --max-length )
      re='^[0-9]+$'
      if ! [[ $2 =~ $re ]] ; then
        echo "Option $2 requires a numeric argument."; exit 1
      fi
      MAX_LENGTH=$2
      shift 2
      ;;

    -f | --first )
      re='^[0-9]+$'
      if ! [[ $2 =~ $re ]] ; then
        echo "Option $2 requires a numeric argument."; exit 1
      fi
      FIRST=$2
      modTRIMM=Y
      shift 2
      ;;

    -l | --last )
      re='^[0-9]+$'
      if ! [[ $2 =~ $re ]] ; then
        echo "Option $2 requires a numeric argument."; exit 1
      fi
      LAST=$2
      modTRIMM=Y
      shift 2
      ;;

    -s | --mismatches )
      re='^[0-9]+$'
      if ! [[ $2 =~ $re ]] ; then
        echo "Option $2 requires a numeric argument."; exit 1
      fi
      MMinput=$2
      shift 2
      ;;

    -e | --extend )
      re='^[0-9]+$'
      if ! [[ $2 =~ $re ]] ; then
        echo "Option $2 requires a numeric argument."; exit 1
      fi
      EXTEND=$2
      shift 2
      ;;

    -j | --slam )
      SLAM=Y
      RERUN_OPTS+=" --slam"
      shift
      ;;

    -D | --dge )
      DGE="Y"
      RERUN_OPTS+=" --dge"
      shift
      ;;

    -x | --force-preprocess )
      FORCE="Y"
      RERUN_OPTS+=" --force-preprocess"
      shift
      ;;

    -P | --raw-paired )
      rawPAIRED="Y"
      RERUN_OPTS+=" --raw-paired"
      shift
      ;;

    -p | --only-paired )
      onlyPAIRED="Y"
      rawPAIRED="Y"
      RAW="Y"
      RERUN_OPTS+=" --only-paired"
      shift
      ;;

    -X | --force-import )
      FORCEimport="Y"
      RERUN_OPTS+=" --force-import"
      shift
      ;;
    -Q | --fastq-out )
      FASTQout="Y"
      RERUN_OPTS+=" --fastq-out"
      shift
      ;;

    -q | --fastq-out-raw )
      FASTQoutRAW="Y"
      RERUN_OPTS+=" --fastq-out-raw"
      shift
      ;;

    -A | --custom-adaptors )
      customADAPTORS="Y"
      RERUN_OPTS+=" --custom-adaptors"
      shift
      ;;

    -r | --raw )
      RAW="Y"
      RERUN_OPTS+=" --raw"
      shift
      ;;

    -T | --trimm )
      TRIMM="Y"
      RERUN_OPTS+=" --trimm"
      shift
      ;;

    -I | --invert )
      INVERT=Y
      RERUN_OPTS+=" --invert"
      shift
      ;;

    -S | --single-end )
      SE=Y
      RERUN_OPTS+=" --single-end"
      shift
      ;;

    -2 | --second-mate )
      SE2nd=Y
      RERUN_OPTS+=" --second-mate"
      shift
      ;;

    -J | --max-count )
      maxCOUNT=$2
      RERUN_OPTS+=" --max-count $2"
      shift 2
      ;;

    -y | --demux-fasta )
      demuxFASTA="Y"
      RERUN_OPTS+=" --demux-fasta"
      shift
      ;;

    -Y | --y-chrom )
      Ychrom="Y"
      RERUN_OPTS+=" --y-chrom"
      shift
      ;;

    -E | --random-multi )
      RANDOMmulti=Y
      RERUN_OPTS+=" --random-multi"
      shift
      ;;

    -W | --wig-norm )
      WIG=Y
      RERUN_OPTS+=" --wig-norm"
      shift
      ;;

    -w | --wig-fasta-norm )
      WIG_FASTA=Y
      RERUN_OPTS+=" --wig-fasta-norm"
      shift
      ;;

    -L | --spike-in-norm )
      spikeINnorm=Y
      RERUN_OPTS+=" --spike-in-norm"
      shift
      ;;

    -z | --no-norm )
      noNORM=Y
      RERUN_OPTS+=" --no-norm"
      shift
      ;;

    --polya )
      POLYA=Y
      RERUN_OPTS+=" --polya"
      shift
      ;;

    -5 | --only-5end )
      only5end=Y
      RERUN_OPTS+=" --only-5end"
      shift
      ;;

    --extra-seq )
      extraSEQ=$2
      RERUN_OPTS+=" --extra-seq $2"
      shift 2
      ;;

    -1 | --ping-pong )
      PingPong=Y
      RERUN_OPTS+=" --ping-pong"
      shift
      ;;

    -n | --no-stranded )
      noSTRANDED=Y
      RERUN_OPTS+=" --no-stranded"
      shift
      ;;
    --demux-only )
      DEMUXonly=Y
      RERUN_OPTS+=" --demux-only"
      shift
      ;;
    -b | --export-bam )
      exportBAM=Y
      RERUN_OPTS+=" --export-bam"
      shift
      ;;

    -B | --export-bam-uncollapsed )
      exportBAMuncollapsed=Y
      RERUN_OPTS+=" --export-bam-uncollapsed"
      shift
      ;;

    -3 | --export-salmon )
      exportSalmon=Y
      RERUN_OPTS+=" --export-salmon"
      shift
      ;;

    -a | --auto-scale )
      autoViewLimits=Y
      RERUN_OPTS+=" --auto-scale"
      shift
      ;;

    --force-quant-unstranded )
      FORCEquant_unstranded=Y
      RERUN_OPTS+=" --force-quant-unstranded"
      shift
      ;;

    -G | --geo )
      GEO=Y
      RERUN_OPTS+=" --geo"
      shift
      ;;

    -g | --ameres )
      AMERES="Y"
      RERUN_OPTS+=" --ameres"
      shift
      ;;

    -U | --user )
      USER=$2
      RERUN_OPTS+=" --user \"$2\""
      shift 2
      ;;

    # you can also add computing/debug/etc if you want them in rerun:
    -C | --computing )
      COMPUTING="P"
      RERUN_OPTS+=" --computing"
      shift
      ;;
    -k | --keep-tmp )
      keepTMP="Y"
      RERUN_OPTS+=" --keep-tmp"
      shift
      ;;
    -O | --folder )
      FOLDER=$2
      RERUN_OPTS+=" --folder \"$2\""
      shift 2
      ;;
    -o | --tmp )
      TMP=$2
      RERUN_OPTS+=" --tmp \"$2\""
      shift 2
      ;;
    -H | --http )
      HTTP=$2
      RERUN_OPTS+=" --http \"$2\""
      shift 2
      ;;
    -K | --update-mode )
      UPDATE="Y"
      RERUN_OPTS+=" --update-mode"
      shift
      ;;
    -Z | --official-ngs-server )
      altSERVER="Y"
      RERUN_OPTS+=" --official-ngs-server"
      shift
      ;;
    -0 | --custom-outdir )
      custom_OUTdir=$2
      RERUN_OPTS+=" --custom-outdir \"$2\""
      shift 2
      ;;
    -d | --debug )
      DEBUG="Y"
      RERUN_OPTS+=" --debug"
      shift
      ;;

    -- )
      shift
      break
      ;;
    * )
      echo "Unknown option: $1" >&2
      usage
      exit 1
      ;;
  esac
done

###################################################################################################
#test if script is in background or in foreground
##exit if in background
case $(ps -o stat= -p $$) in
*+*) ;;
*)
  printf "\n\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
  printf "script was started in the background\n"
  printf "restarted the script without the & \n"
  printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
  exit
  ;;
esac


#--------------------------------------------------------------------------------------------------
#reset BASE_FOLDER if custom output directory is set
if [[ ! -z "${custom_OUTdir}" ]]; then
  BASE_FOLDER="${custom_OUTdir}"
  mkdir -p $BASE_FOLDER
fi  

#--------------------------------------------------------------------------------------------------
#determine the location of the top script
##set script-dir automatically

SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$(cd -P "$(dirname "$SOURCE")" && pwd)"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
SCRIPT_DIRraw="$(cd -P "$(dirname "$SOURCE")" && pwd)"
SCRIPT_DIRraw="${SCRIPT_DIRraw}/"
SCRIPT_DIR="${SCRIPT_DIRraw}script-files/"
UTILITY_DIR="${SCRIPT_DIRraw}utility-files/"
SINGULARITYdir="${BASE_UTILITY_LOCATION}simg/"


#!hack to copy singularity containers to the home directory as there is some problem in peters environment executing from network storage
rsync -av --checksum "${SINGULARITYdir}" ~/AP_singu/
SINGULARITYdir=~/AP_singu/


#source all tools from singularity containers
source ${SCRIPT_DIR}tools

#--------------------------------------------------------------------------------------------------
#test if obligatory variables are set
if [[ -z "${FILE_CONTAINING_LIBRARIESraw}" ]] && [[ -z "${UPDATErun}" ]]; then
  #usage
  printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
  printf "          At least option -i or -u have to be set!\n"
  printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n"
  exit
fi

#test if maxCOUNT is numeric
if [[ ! $maxCOUNT -gt 0 ]]; then
  #usage
  printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
  printf "          maxCOUNT in -J has to be a number > 0 !\n"
  printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n"
  exit
fi

#test if extra-seq is fasta file
if [[ ! -z "${extraSEQ}" ]]; then
  if [[ ! -s "${extraSEQ}" ]]; then
    printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
    printf "          extra-seq file does not exist or is empty: ${extraSEQ} \n"
    printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n"
    exit
  fi
  FORMATcheck=$(head -n 1 ${extraSEQ} | cut -c1)
  if [[ $FORMATcheck != ">" ]]; then
    printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
    printf "          extra-seq file is not in fasta format: ${extraSEQ} \n"
    printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n"
    exit
  fi
fi

###################################################################################################
#prepare everything for restart with run-update

if [[ ! -z "${UPDATErun}" ]]; then

  UPDATErun=$(echo "$UPDATErun" | tr -s /)
  FOLDERforUPDATE="${BASE_FOLDER}/${UPDATErun}/"
  TMPforUPDATE="${BASE_FOLDER_TMP}/${UPDATErun}/"

  if [[ -d "${BASE_FOLDER}/${UPDATErun}_updated-to-${PIPELINEversion}/" ]]; then
    while true; do
      read -r -p "
      Run already got updated to version ${PIPELINEversion}. 
      If you re-run the update process, the old updated-version gets deleted and 
      the run is started from the original run settings.
      
      Do you want to update again? [y or n]" yn
      case $yn in
      [Yy])
        mv "${BASE_FOLDER}/${UPDATErun}_updated-to-$PIPELINEversion" ${FOLDERforUPDATE}
        break
        ;;
      [Nn])
        printf "\n\nStopped run due to User input\n\n"
        exit
        ;;
      *) echo "Please answer yes [y] or no [n]." ;;
      esac
    done
  else
    printf "\n\n  no idea what is missing here - ask Dominik \n\n"
    exit
  fi

  mkdir -p $TMPforUPDATE
  if [[ -d $FOLDERforUPDATE ]]; then

    newSTART=$(cat ${FOLDERforUPDATE}log.txt | grep -A 1 "start command:" | tail -n 1 |
      awk -v FOLDER=$FOLDERforUPDATE -v TMP=$TMPforUPDATE -v FILEfile=${TMPforUPDATE}files_for_update.txt -v HTTP=${BASE_HTTP}${UPDATErun} -v SCRIPT=$0 '
      BEGIN{
        print SCRIPT
      }
      {
        for(i=2; i<=NF; i++) {
          if( $i == "-V") {
            i=i+1
          }else{
            if( $i == "-i"){
              print $i, FILEfile
              i=i+1
            }else{
              print $i
            }
          }
        }
      }
      END{
        print "-O",FOLDER,"-K","-o",TMP,"-H",HTTP
      }' | tr '\n ' ' ')

    ##extract liebraries to analyze and compile input-file
    awk '{
      if($0~"libraries are available") {p=0}
      if(p==1){
      print 
      }
      if($0~"libraries analyzed") {p=1}
    }' ${FOLDERforUPDATE}log.txt >${TMPforUPDATE}files_for_update.txt
    FILE_CONTAINING_LIBRARIESraw=${TMPforUPDATE}files_for_update.txt
    ##extract settings and compile the new submission command

  else
    #complain that the supplied folder does not exist
    usage
    printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
    printf "          The folder specified to update does not exist!\n"
    printf "${BASE_HTTP}${UPDATErun}\n"
    printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n"
    exit
  fi
fi

###################################################################################################
#set type and/or validate provided type class

if [[ -z $TYPE ]]; then
  printf "\navailabel TYPE classes:\n"
  cat "${UTILITY_DIR}available_TYPES.txt"
  while true; do
    read -r -p "Which Type do you want to use? give number from list:  " nr
    case $nr in
    [1-9] | [1-9][0-9])
      TYPE=$(cat "${UTILITY_DIR}available_TYPES.txt" | awk -v nr=$nr '{if($1==nr) print $2 }')
      if [[ -z $TYPE ]]; then
        printf "\nselected TYPE does not exist. supply valid value!!!\n\n "
      else
        break
      fi
      ;;
    *) printf "\nsupply valid value \n\n" ;;
    esac
  done
elif [[ $TYPE == "RNAseq" ]] || [[ $TYPE == "sRNAseq" ]] || [[ $TYPE == "sRNAseqIP" ]] || [[ $TYPE == "GROseq" ]] || [[ $TYPE == "CHIPseq" ]] || [[ $TYPE == "DNAseq" ]] || [[ $TYPE == "CLIPseq" ]] || [[ $TYPE == "RNAseq_QuantSeq" ]] || [[ $TYPE == "CapSeq" ]] || [[ $TYPE == "RIPseq" ]]; then
  sleep 0.1s
else
  printf "\n\n\nthe TYPE you provided does not exist!!!!\n"
  printf "\navailabel TYPE classes:\n"
  cat "${UTILITY_DIR}available_TYPES.txt"
  while true; do
    read -r -p "Which Type do you want to use? give nunber from list:  " nr
    case $nr in
    [1-9] | [1-9][0-9])
      TYPE=$(cat "${UTILITY_DIR}available_TYPES.txt" | awk -v nr=$nr '{if($1==nr) print $2 }')
      if [[ -z $TYPE ]]; then
        printf "\nselected TYPE does not exist. supply valid value!!!\n\n "
      else
        printf "as TYPE #$TYPE = $TYPE used\n"
        break
      fi
      ;;
    *) printf "\nsupply valid value \n\n" ;;
    esac
  done
fi

###################################################################################################
#test user name

##determin user-name or set to provided user-variable
##before generating new user-FOLDER ask for permission to reduce missspelling errors
if [[ -z $USER ]]; then
  USER_NAME=$(whoami)
else
  echo "running pipeline for other user"

  USER=$(echo "$USER" | tr '[:upper:]' '[:lower:]')
  if [[ ! -d "${BASE_FOLDER}/${USER}" ]]; then
    while true; do
      read -r -p "User name is new. Are you sure it is correctly spelled? Generate new user FOLDER? [y or n]" yn
      case $yn in
      [Yy])
        USER_NAME=$USER
        break
        ;;
      [Nn])
        printf "\n\nprovided user was: $USER\n\n"
        exit
        ;;
      *) echo "Please answer yes [y] or no [n]." ;;
      esac
    done
  else
    USER_NAME=$USER
  fi
fi

##get date and genome VERSION
DATE_OF_DAY=$(date +%F)
#DATE_OF_DAY="2018-08-29"

###################################################################################################
#determine the genome VERSION to use
if [[ -z $GENOME_VERSION ]]; then
  while true; do
    read -r -p "Which genome VERSION should be used? [dm3; dm6; Cel; or ASM (assembly)]" yn
    case $yn in
    dm3)
      GENOME_VERSION="dm3"
      break
      ;;
    dm6)
      GENOME_VERSION="dm6"
      break
      ;;
    ASM)
      GENOME_VERSION="ASM"
      break
      ;;
    Cel)
      GENOME_VERSION="Cel"
      break
      ;;
    *) echo "Please answer dm3, dm6, Cel or ASM." ;;
    esac
  done
else
  if [[ $GENOME_VERSION == "dm3" || $GENOME_VERSION == "dm6" || $GENOME_VERSION == "Cel" || $GENOME_VERSION == "ASM" ]]; then
    sleep 1s
  else
    printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
    printf "\n  Genome Version set wrong \n"
    printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n"
    while true; do
      read -r -p "Please select a valid genome version! [dm3, dm6 or ASM]" yn
      case $yn in
      dm3)
        GENOME_VERSION="dm3"
        break
        ;;
      dm6)
        GENOME_VERSION="dm6"
        break
        ;;
      ASM)
        GENOME_VERSION="assembly"
        break
        ;;
      Cel)
        GENOME_VERSION="Cel"
        break
        ;;
      *) echo "Please answer dm3, dm6 or ASM." ;;
      esac
    done
  fi
fi

if [[ $GENOME_VERSION == dm* ]]; then
  UTILITY_LOCATION=${BASE_UTILITY_LOCATION}dmel/${GENOME_VERSION}/
elif [[ $GENOME_VERSION == Cel ]]; then
  UTILITY_LOCATION=${BASE_UTILITY_LOCATION}Cel/
elif [[ $GENOME_VERSION == ASM ]]; then
  UTILITY_LOCATION=${BASE_UTILITY_LOCATION}ASM/
fi


###################################################################################################
#determine which annotation VERSION the user wants to use
## or if a new one should be created

#create empty VERSION file if not existing
if [[ ! -f ${BASE_UTILITY_LOCATION}/versions.txt ]]; then touch ${BASE_UTILITY_LOCATION}/VERSIONs.txt; fi

RELEASE=
asmHUBpath=
ASMdir=
ASMname=
newVERSION="no"
prepareREF="no"

#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------
#dm# genome versions from flybase

#---------------------------------------------------------------------------------------------------------
#test if selected version exists, if not ask to create
if [[ ! -z $VERSION ]]; then
  VERSIONexists=$(grep -w $VERSION ${BASE_UTILITY_LOCATION}/versions.txt)
  if [[ -z $VERSIONexists && $GENOME_VERSION != ASM* ]]; then
    while true; do
      read -r -p "Selected genome-annotation version does not exist locally. Do you want to install it automatically? [y or n]" yn
      case $yn in
      [Nn])
        VERSION=""
        break
        ;;
      [Yy])
        newVERSION="yes"
        break
        ;;
      *) echo "Please answer yes [y] or no [n]." ;;
      esac
    done
  elif [[ -z $VERSIONexists ]]; then
    #wipe Version variable for non dm6 genomes if the version does not exist in versions.txt
    VERSION=""
  fi
fi


#---------------------------------------------------------------------------------------------------------
#if no version supplied - ask if default should be used
prepANNOTATIONgff=""
prepGENOMEfasta=""
prepTRANSCRIPTOMEfasta=""
prepCDSfasta=""
prepNCRNAfasta=""

if [[ $GENOME_VERSION != ASM* ]]; then
  if [[ -z $VERSION ]]; then
    #initiate TMP directory for version determination files
    TMP=${BASE_FOLDER_TMP}/${USER_NAME}/
    mkdir -p $TMP

    while true; do
      if [[ $GENOME_VERSION == dm3 ]]; then
        read -r -p "Run pipeline with the default annotation ${default_VERSION_dm3}? [y or n] " yn
      elif [[ $GENOME_VERSION == dm6 ]]; then
        read -r -p "Run pipeline with the default annotation ${default_VERSION_dm6}? [y or n] " yn
      elif [[ $GENOME_VERSION == Cel ]]; then
        read -r -p "Run pipeline with the default annotation ${default_VERSION_Cel}? [y or n] " yn
      fi

      case $yn in
      [Yy])
        if [[ $GENOME_VERSION == dm3 ]]; then VERSION=${default_VERSION_dm3}; elif [[ $GENOME_VERSION == dm6 ]]; then VERSION=${default_VERSION_dm6}; elif [[ $GENOME_VERSION == Cel ]]; then VERSION=${default_VERSION_Cel} ;fi
        VERSIONexists=$(grep -w $VERSION ${BASE_UTILITY_LOCATION}/versions.txt)
        if [[ -z $VERSIONexists ]]; then
          while true; do
            read -r -p "Selected genome-annotation version does not exist locally. Do you want to install it? [y or n]" yn
            case $yn in
            [Nn])
              VERSION=""
              break
              ;;
            [Yy])
              newVERSION="yes"
              break
              ;;
            *) echo "Please answer yes [y] or no [n]." ;;
            esac
          done
        fi
        break
        ;;
      [Nn])
        if [[ $GENOME_VERSION == dm3 ]]; then
          printf "available annotations for dm3: \n"
          grep r5. ${BASE_UTILITY_LOCATION}/versions.txt > ${TMP}available_versions.txt
          cat ${TMP}available_versions.txt
          printf "\n"
        elif [[ $GENOME_VERSION == dm6 ]]; then
          printf "available annotations for dm6: \n"
          grep r6. ${BASE_UTILITY_LOCATION}/versions.txt > ${TMP}available_versions.txt
          cat ${TMP}available_versions.txt
          printf "\n"
        elif [[ $GENOME_VERSION == Cel ]]; then
          printf "available annotations for Cel: \n"
          grep Cel ${BASE_UTILITY_LOCATION}/versions.txt > ${TMP}available_versions.txt
          cat ${TMP}available_versions.txt
          printf "\n"
        fi

        #determine which annotation to use
        while true; do
          read -r -p "Define which genome-annotation of the listed to use. Set 0 if new annotation has to be generated: " ann
          case $ann in
          [0])
            printf "generate new annotation VERSION\n\n"
            newVERSION="yes"
            break
            ;;
          [1-9] | [1-9][0-9] | [1-9][0-9][0-9])
            VERSION=$(cat ${TMP}available_versions.txt | awk -v ANN=$ann '{if($1== ANN) print $2 }')
            if [[ -z $VERSION ]]; then
              printf "\nselected annotation does not exist. supply valid value!!!\n\n "
            else
              printf "annotation #$ann = $VERSION used\n"
              break
            fi
            ;;
          *) printf "\nsupply valid value \n\n" ;;
          esac
        done
        break
        ;;
      *) echo "Please answer yes [y] or no [n]." ;;
      esac
    done
  fi

  #---------------------------------------------------------------------------------------------------------
  #determine which annotation to create
  if [[ $newVERSION == yes ]]; then

    GREP=$(which grep)
    source ${SCRIPT_DIR}tools

    cd $TMP
    
    #get index of available annotations
    if [[ ! -s ${TMP}version-ligint ]]; then
      if [[ $GENOME_VERSION == dm* ]]; then
        #!@ ncftpls ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/* |  grep -F dmel | grep -F -v current >${TMP}version-listing
        curl -s https://s3ftp.flybase.org/genomes/Drosophila_melanogaster/ | grep -o "<a href='[^']*/index.html'>[^<]*</a>" | sed -E "s|.*'([^/]*)/index.html'.*>|\1|" | grep -F dmel_r | sort -u >${TMP}version-listing
      elif [[ $GENOME_VERSION == Cel ]]; then
        ncftpls ftp://ftp.wormbase.org/pub/wormbase/releases/  | grep -F WS | grep -F -v current >${TMP}version-listing
      fi
    fi

    #process full listing to 
    cat ${TMP}version-listing |
      awk -v OFS="\t" -v GENOME_VERSION=$GENOME_VERSION 'BEGIN {N=1; print "#","version","release_year" } { 
        gsub("/","")
        if( GENOME_VERSION == "dm3" && $NF ~ "r5.") {
          split($NF,X,/_/); split(X[2],Y,/\./)
          if( Y[2] > 40) { print N,X[2],X[3]"_"X[4]; N=N+1 }
        }else{
          if( GENOME_VERSION == "dm6" && $NF ~ "r6.") {
            split($NF,X,/_/); split(X[2],Y,/\./); print N,X[2],X[3]"_"X[4]; N=N+1 
          }else{
            if( GENOME_VERSION == "Cel" && $NF ~ "WS") {
              print N,$0; N=N+1
            }
          }
        }
      }'  >${TMP}available_versions.txt

    # ask user to provide a flybase version
    if [[ $GENOME_VERSION == dm6 ]]; then
        while true; do
            read -p "Please provide the flybase ID (e.g. r6.63_FB2025_02): " NR
            # Check if the input matches the expected format
            if [[ $NR =~ ^r6\.[0-9]+_FB[0-9]{4}_[0-9]{2}$ ]]; then
                echo "FlyBase ID accepted: $NR"
                FULLid="$NR"
                break
            else
                echo "Invalid format. Please use the format: r6.63_FB2025_02"
            fi
        done
    fi

    VERSION=$(echo $FULLid | tr '_' '\t' | cut -f 1 )
    RELEASE=$(echo $FULLid | tr '_' '\t' | cut -f 2-3 | tr '\t' '_' )


    #ask user to provide path to the genome-annotation while providing a download link for thme:
    while true; do
      printf "\n\nDownload the genome-annotation file from FlyBase using the link below and store on cluster-accesible storage:\n"
      echo https://s3ftp.flybase.org/genomes/Drosophila_melanogaster/dmel_${VERSION}_${RELEASE}/gff/dmel-all-${VERSION}.gff.gz
      read -p "Please provide the path to the genome-annotation file (e.g. dmel_r6.63_FB2025_02.gff.gz): " prepANNOTATIONgff
      if [[ -f $prepANNOTATIONgff ]]; then
        printf "Annotation file found: $prepANNOTATIONgff \n\n"
        break
      else
        echo "File not found. Please provide a valid path."
      fi
    done


    #same for dmel-all-${refTYPE}-${VERSION}.fasta.gz
    #ask user to provide path to the genome-annotation while providing a download link for thme:
    while true; do
      printf "\n\nDownload the genome-sequence file from FlyBase using the link below and store on cluster-accesible storage:\n"
      echo https://s3ftp.flybase.org/genomes/Drosophila_melanogaster/dmel_${VERSION}_${RELEASE}/fasta/dmel-all-chromosome-${VERSION}.fasta.gz
      read -p "Please provide the path to the genome-sequence file (e.g. dmel-all-chromosome-r6.63.fasta.gz): " prepGENOMEfasta
      if [[ -f $prepGENOMEfasta ]]; then
        printf "Sequence file found: $prepGENOMEfasta \n\n"
        break
      else
        echo "File not found. Please provide a valid path."
      fi
    done

    #dmel-all-${refTYPE}-${VERSION}.fasta.gz
     #ask user to provide path to the genome-annotation while providing a download link for thme:
    while true; do
      printf "\n\nDownload the transcriptome-sequence file from FlyBase using the link below and store on cluster-accesible storage:\n"
      echo https://s3ftp.flybase.org/genomes/Drosophila_melanogaster/dmel_${VERSION}_${RELEASE}/fasta/dmel-all-transcript-${VERSION}.fasta.gz
      read -p "Please provide the path to the transcript-sequence file (e.g. dmel-all-transcript-r6.63.fasta.gz): " prepTRANSCRIPTOMEfasta
      if [[ -f $prepTRANSCRIPTOMEfasta ]]; then
        printf "Sequence file found: $prepTRANSCRIPTOMEfasta \n\n"
        break
      else
        echo "File not found. Please provide a valid path."
      fi
    done

    #dmel-all-${refTYPE}-${VERSION}.fasta.gz
     #ask user to provide path to the genome-annotation while providing a download link for thme:
    while true; do
      printf "\n\nDownload the CDS-sequence file from FlyBase using the link below and store on cluster-accesible storage:\n"
      echo https://s3ftp.flybase.org/genomes/Drosophila_melanogaster/dmel_${VERSION}_${RELEASE}/fasta/dmel-all-CDS-${VERSION}.fasta.gz
      read -p "Please provide the path to the CDS-sequence file (e.g. dmel-all-CDS-r6.63.fasta.gz): " prepCDSfasta
      if [[ -f $prepCDSfasta ]]; then
        printf "Sequence file found: $prepCDSfasta \n\n"
        break
      else
        echo "File not found. Please provide a valid path."
      fi
    done

    while true; do
      printf "\n\nDownload the ncRNA-sequence file from FlyBase using the link below and store on cluster-accesible storage:\n"
      echo https://s3ftp.flybase.org/genomes/Drosophila_melanogaster/dmel_${VERSION}_${RELEASE}/fasta/dmel-all-ncRNA-${VERSION}.fasta.gz
      read -p "Please provide the path to the ncRNA-sequence file (e.g. dmel-all-ncRNA-r6.63.fasta.gz): " prepNCRNAfasta
      if [[ -f $prepNCRNAfasta ]]; then
        printf "Sequence file found: $prepNCRNAfasta \n\n"
        break
      else
        echo "File not found. Please provide a valid path."
      fi
    done



    RELEASEloop=n
    while [[ $RELEASEloop != Y ]]; do
      if [[ -z $VERSION ]]; then
        #print available flybase versions
        cat ${TMP}available_versions.txt

        #ask user to provide a flybase version
        while true; do
          read -p "Which genome-annotation version to use? [# in list above]" NR
          case $NR in
          [1-9] | [1-9][0-9] | [1-9][0-9][0-9]*)
            cat ${TMP}available_versions.txt | awk -v TMP=$TMP -v xNR=$NR ' BEGIN {xswitch=0}; 
                   { if($1==xNR) { xswitch=1; print "genome-annotation version "$2" selected" ; print $2 >TMP"ANNversion.txt" } }
                        END { print xswitch > TMP"switch.txt"}'
            switch=$(cat ${TMP}switch.txt)
            if [[ $switch == 1 ]]; then break; else
              cat ${TMP}available_versions.txt
              printf "\nplease select existing genome-annotation version \n\n"
            fi
            ;;
          *) echo "Please provide valid value." ;;
          esac
        done

        #create variable for download
        RELEASEloop=Y
        VERSION=$(cat ${TMP}ANNversion.txt)
        RELEASE=$(cat ${TMP}available_versions.txt | $GREP -F $VERSION | cut -f 3 | sed 's/\r$//')
      else
        cat ${TMP}available_versions.txt | awk -v TMP=$TMP -v VERSION=$VERSION ' BEGIN {xswitch=0}; 
        { if($2==VERSION) { xswitch=1; print "genome-annotation version "$2" selected" ; print $2 >TMP"ANNversion.txt" } }
             END { print xswitch > TMP"switch.txt"}'
        switch=$(cat ${TMP}switch.txt)
        if [[ $switch == 1 ]]; then
          RELEASEloop=Y
          RELEASE=$(cat ${TMP}available_versions.txt | $GREP -F $VERSION | cut -f 3 | sed 's/\r$//')
        else
          VERSION=""
          printf "\nsupplied version does not exist on FlyBase/Wormbase. Please select existing version! \n\n"
        fi

      fi
    done

    #report selected version into the versions file
    #test if version of annotation already exists
    VERSION_EXISTS_TEST=$(cat ${BASE_UTILITY_LOCATION}versions.txt | $GREP -F $VERSION)
    if [[ ! -z $VERSION_EXISTS_TEST ]]; then
      while true; do
        read -p "Annotation for selected genome-annotation version exists. Do you want to rerun the annotation perparation? [y or n]" yn
        case $yn in
        [Yy]*)
          XY=$(echo $VERSION_EXISTS_TEST | tr ' ' '\t' | cut -f 2-3 | tr '\t' ' ')
          printf "\n\n rerunning annotation perparation for $XY \n\n"
          break
          ;;
        [Nn]*)
          printf "\n\n  Run aborted due to selection of wrong genome-annotation version during annotation prperation \n\n" exit
          break
          ;;
        *) echo "Please answer Y or N." ;;
        esac
      done
    fi

    prepareREF="yes"
  fi
fi

#---------------------------------------------------------------------------------------------------------
#ASM genome versions from assembly

if [[ $GENOME_VERSION == ASM ]]; then
  TMP=${BASE_FOLDER_TMP}/${USER_NAME}/
  mkdir -p $TMP

  #test if selected assembly-version 
  TESTversion=
  if [[ -n $VERSION ]]; then
    TESTversion=$(grep -w $VERSION ${BASE_UTILITY_LOCATION}/versions.txt) 
  fi

  if [[ -n $TESTversion && -d ${BASE_FOLDER}ASSEMBLIES/${VERSION}/ ]]; then 
    echo using existing assembly version $VERSION
  else
    #warn user that assemly version of the pipeline was selected and present existing assemblies
    ##available in the versions.txt file

    printf "available assemblies: \n"
    grep ASM ${BASE_UTILITY_LOCATION}/versions.txt > ${TMP}available_versions.txt
    printf "\n0\tgenerate new ASM index\n\n" >>${TMP}available_versions.txt
    cat ${TMP}available_versions.txt
    printf "\n"

    if [[ -d ${BASE_FOLDER}ASSEMBLIES/${VERSION}/ ]]; then
      #determine which assembly to use
      while true; do
        read -r -p "Define which assembly of the listed to use. Set 0 if index for new assembly has to be generated: " ann
        case $ann in
        [0])
          printf "generate index for new assemly\n\n"
          newVERSION="yes"
          break
          ;;
        [1-9] | [1-9][0-9] | [1-9][0-9][0-9])
          VERSION=$(cat ${TMP}available_versions.txt | awk -v ANN=$ann '{if($1== ANN) print $2 }')
          if [[ -z $VERSION ]]; then
            printf "\nselected assembly does not exist. supply valid value!!!\n\n "
          else
            printf "assembly #$ann = $VERSION used\n"
            break
          fi
          ;;
        *) printf "\nsupply valid value \n\n" ;;
        esac
      done
    else
      echo requested VERSION not built yet 
      newVERSION=yes
    fi

    #if new assemlby is requested ask the user for the required input
    if [[ $newVERSION == yes ]]; then
      #set function to test URLs
      function validate_url() {
        if [[ -s $1 ]]; then
          return 0
        else
          return 1
        fi
      }
      export -f validate_url

      #ask for link to the assembly-hub
      while true; do
        read -r -p "Please provide the path to the hub.txt file from the assembly hub: " hub
        case $hub in
        *)
          if validate_url $hub; then
            #ask user if this is the correct version after processing all variables
            VALIDATEhub=$(cat $hub | head -n 1 | tr ' ' '\t')
            if [[ $VALIDATEhub == hub* ]]; then
              asmHUBpath=$(dirname $hub)
              asmHUBpath="${asmHUBpath}/"
              ASMname=$(echo $VALIDATEhub | tr ' ' '\t' | cut -f 2)
              while true; do
                read -r -p "Assembly is currently named $ASMname. To change the name type the name to keep type n: " keepname
                case $keepname in
                "") ;;

                [Nn])
                  VERSION=$ASMname
                  break
                  ;;
                *)
                  VERSION=$(echo $keepname | tr ' ' '_')
                  break
                  ;;
                esac
              done
              printf "\n\n  Assembly will be added as $VERSION \n\n"
              ASMdir=${BASE_FOLDER}ASSEMBLIES/${VERSION}/
              break
            else
              printf "\nnot a functional link to a hub.txt file!!!\n\n "
            fi
          else
            #report that supplied link is not in a valid format and why
            printf "\nnot a functional link to a hub.txt file!!!\n\n "
          fi
          ;;
        esac
      done
      prepareREF="yes"
    fi
  fi
  rm -rf ${TMP}available_versions.txt
  UTILITY_LOCATION=${UTILITY_LOCATION}${VERSION}/
  ASMdir=${BASE_FOLDER}ASSEMBLIES/${VERSION}/

fi


#determine folder name and if not set ask back if it should get set
if [[ $UPDATE != Y ]]; then
  if [[ -z $FOLDER_NAME ]]; then
    printf "\n\n"
    while true; do
      read -r -p "No folder-name-addition given.
  Results will be located in ${BASE_FOLDER}/${USER_NAME}/${GENOME_VERSION}/${TYPE}/${DATE_OF_DAY}
  Provide name added to the end of aboves name or type [n] to stay with the one above: " yn
      case $yn in
      [Nn]) break ;;
      *)
        FOLDER_NAME=$yn
        break
        ;;
      esac
    done
  fi

  #set FOLDER variables
  COLLECTIONname=
  newCOLLECTION=N
  if [[ $GENOME_VERSION == ASM ]]; then 
    #print available runs and ask if addint to one of these or running new/basic
    if [[ ! -d ${BASE_FOLDER}/${USER_NAME}/${VERSION}/ ]]; then
      printf "\nfirst run for assembly ${VERSION}\n"
      printf "\ndo you want to add this pipeline-run to the default collection or do you want to create a named collection of pipeline-runs?\n"
      while true; do
        read -r -p "Provide name to create new collection or type [n] to stay with the default collection: " yn
        case $yn in
        [Nn]) COLLECTIONname=default; break ;;
        *)
          COLLECTIONname=$yn
          break
          ;;
        esac
      done
      newCOLLECTION=Y
    else
      #determine available collections
      ls ${BASE_FOLDER}/${USER_NAME}/${VERSION}/ | awk 'BEGIN{n=1}{if($1=="default") {print 0,$1;} else {print n,$1; n+=1}}' | sort -k1,1n> ${TMP}existing-collections.txt
      #determine which assembly to use
      while true; do
        printf "\n\n\n"
        cat ${TMP}existing-collections.txt
        printf "\n"
        read -r -p "Please select an existing collection (NR) or type a new collection-name " ann
        case $ann in
        [0-9] | [1-9][0-9] | [1-9][0-9][0-9])
          COLLECTIONname=$(cat ${TMP}existing-collections.txt | awk -v ANN=$ann '{if($1== ANN) print $2 }')
          break
          ;;
        *) 
          while true; do
            read -p "new collection will be named $ann - is this correct? [y/n]"  yn
            case $yn in
            [Yy]*)
              COLLECTIONname=$ann
              newCOLLECTION=Y
              break
              ;;
            [Nn]*)
              break
              ;;
            *) echo "Please answer Y or N." ;;
            esac
          done
          if [[ -n $COLLECTIONname ]]; then
            break
          fi        
        esac
      done
    fi
    VERSIONvari=$VERSION/$COLLECTIONname/AP_runs
    GENOMEdir=${BASE_FOLDER}/${USER_NAME}/$VERSION/$COLLECTIONname/
  else 
    VERSIONvari=$GENOME_VERSION
    GENOMEdir=${BASE_FOLDER}/${USER_NAME}/$GENOME_VERSION/
  fi

  if [[ -z ${FOLDER} ]]; then
    if [[ -z $FOLDER_NAME ]]; then
      FOLDER=${BASE_FOLDER}/${USER_NAME}/${VERSIONvari}/${TYPE}/${DATE_OF_DAY}/
      TMP=${BASE_FOLDER_TMP}/${USER_NAME}/${VERSIONvari}/${TYPE}/${DATE_OF_DAY}/
      HTTP="${BASE_HTTP}${USER_NAME}/${VERSIONvari}/${TYPE}/${DATE_OF_DAY}/"
      RUNname=${DATE_OF_DAY}
    else
      FOLDER_NAME=$(echo $FOLDER_NAME | tr ' ' '_' | tr '/' '_')
      FOLDER=${BASE_FOLDER}/${USER_NAME}/${VERSIONvari}/${TYPE}/${DATE_OF_DAY}-${FOLDER_NAME}/
      TMP=${BASE_FOLDER_TMP}/${USER_NAME}/${VERSIONvari}/${TYPE}/${DATE_OF_DAY}-${FOLDER_NAME}/
      HTTP="${BASE_HTTP}${USER_NAME}/${VERSIONvari}/${TYPE}/${DATE_OF_DAY}-${FOLDER_NAME}/"
      RUNname=${DATE_OF_DAY}-${FOLDER_NAME}
    fi
  fi

  #--------------------------------------------------------------------------------------------------
  #test if FOLDERs already exist and if so promt for replacement
  SWITCH=0
  while [[ $SWITCH -eq 0 ]]; do
    if [[ -d $FOLDER ]]; then
      while true; do
        read -r -p "FOLDER on open directory exists, do you wish to overwrite? [y or n or new-name]" yn
        case $yn in
        [Yy])
          SWITCH=1
          break
          ;;
        [Nn]) exit ;;
        *)
          FOLDER_NAME=$yn
          FOLDER=${BASE_FOLDER}/${USER_NAME}/${VERSIONvari}/${TYPE}/${DATE_OF_DAY}-${FOLDER_NAME}/
          TMP=${BASE_FOLDER_TMP}/${USER_NAME}/${VERSIONvari}/${TYPE}/${DATE_OF_DAY}-${FOLDER_NAME}/
          break
          ;;
        esac
      done
    else
      SWITCH=1
    fi
  done

  #remove old folder structures
  if [[ $COMPUTING == C && $DEBUG == N ]]; then
    rm -rf $FOLDER
    rm -rf $TMP
  fi
fi

#initiate TMP
mkdir -p $TMP

#---------------------------------------------------------------------------------------------------------
#determine color code for TYPE

if [[ -z ${COLOR+x} ]]; then
  COLOR=$(grep -w ${TYPE} ${UTILITY_DIR}track-colors_for_TYPES.txt | awk '{print $NF}')
  if [[ $SLAM == Y ]]; then
    COLOR="0,200,255"
  fi
else
  STATUS=$(echo $COLOR | awk -v FS="," '{
    if(NF != 3) print "nOK"; else {
      for (i=1; i<=3; i++) { if ($i <0 || $i >256) print "xOK" }
    }
    }')
  if [[ ! -z ${STATUS} ]]; then
    usage
    printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
    printf "  provided color in wrong format \n"
    printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n"
    exit 1
  fi
fi

subCOLOR=$(echo $COLOR | tr ',', '~')

#---------------------------------------------------------------------------------------------------------
#test supplied subsampling value

if [[ ! -z ${SUBSAMPLE} ]]; then
  TESTcond='^[0-9]+$'
  if ! [[ $SUBSAMPLE =~ $TESTcond ]]; then
    # usage
    printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
    printf "  provide numeric value for subsampling \n"
    printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n"
    exit 1
  fi
fi

###################################################################################################
##check if variables are filled or create optional variables

#check obligatory variables filled
if [[ -z $FILE_CONTAINING_LIBRARIESraw ]]; then
  usage
  printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
  printf "\n  please provide input libaries - set -i \n\n"
  printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n"
  exit 1
elif [[ ! -f $FILE_CONTAINING_LIBRARIESraw ]]; then
  printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
  printf "\n  path supplied by -i is not valid \n  no file found at: $FILE_CONTAINING_LIBRARIESraw \n\n"
  printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n"
  exit 1
else
  #create temporary directory for library-availability check
  TMPtmp=${TMP}/TMPtmp/
  mkdir -p $TMPtmp

  #switch server to alternate server if requested
  if [[ $altSERVER == Y ]]; then 
    SERVER="altSERVER"
  fi
  
  #clean up of input file and only keep first library for local computing
  if [[ $COMPUTING == P ]]; then
    sed '/^\s*$/d' "$FILE_CONTAINING_LIBRARIESraw" | sed -e '$a\' | head -n 3 | tail -n 5 | tr ' ' '\t' | sed 's/\r$//' |
      awk -v OFS="\t" -v TMP=$TMP '{ if($2=="") print >TMP "missing-column.txt"; else print }' |
      awk -v OFS="\t" -v SERVER=${SERVER} '{ if($1~"http") {split($1,X,/data/); $1="https://" SERVER "/data" X[2]}; print}' |
      sort -k1,1r >"${TMPtmp}files.txt"
  else
    sed '/^\s*$/d' "$FILE_CONTAINING_LIBRARIESraw" | sed -e '$a\' | tr ' ' '\t' | sed 's/\r$//' |
      awk -v OFS="\t" -v TMP=$TMP '{ if($2=="") print X >TMP "missing-column.txt"; else print }' |
      awk -v OFS="\t" -v SERVER=${SERVER} '{ if($1~"http") {split($1,X,/data/); $1="https://" SERVER "/data" X[2]}; print}' |
      sort -k1,1r >"${TMPtmp}files.txt"
  fi

  #test if every line contains at least 2 columns (path + name)
  if [[ -f "${TMP}missing-column.txt" ]]; then
    printf "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
    printf "\n  file provided in -i does not contain the required 2 columns for following libraries  \n\n"
    cat "${TMP}missing-column.txt"
    printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n"
    exit 1
  fi

  #test if duplicated names exist
  rm -rf ${TMP}duplicated.txt
  awk -v TMP=$TMP '{
    X[$2]+=1  
  }
  END{
    for(i in X){
      if(X[i]>1){
        print i > TMP "duplicated.txt"
      }
    }
  }' ${TMPtmp}files.txt

  if [[ -f "${TMP}duplicated.txt" ]]; then
    printf "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
    printf "\nthe following names you provided are provided more than once in the file:   \n\n"
    cat "${TMP}duplicated.txt"
    printf "\nplease fix the names - only unique names allowed   \n\n"
    printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n"
    exit 1
  fi

  cp ${TMPtmp}files.txt ${TMPtmp}files_orig.txt

  #---------------------------------------------------------------------------------------------------------
  #determine correctnes and availability of NGS user credentials

  #test if NGS-login required and if so - ask user for the login details
  ngsUSER=$(
    awk '
    {
      if($1 ~ /^NGS/ || $1 ~ /^http/ || $1 ~ /^TGZ/ ){
        print "YES"
        exit
      }
    }' "${TMPtmp}files.txt"
  )

  LISTING=$(
    awk '
    {
      if($1 ~ /^NGS/ || $1 ~ /^http/ || $1 ~ /^TGZ/ ){
        print "YES"
        exit
      }
    }' "${TMPtmp}files.txt"
  )

  if [[ ${ngsUSER} == YES ]]; then
    SWITCH="N"

    #test if netrc-file already exists
    if [[ -s ~/.AP_login.txt ]]; then
      loginFILE=Y
    else
      loginFILE=N
    fi
    #ask user for login detals if netrc-file not available or login&PW are wrong
    ##continue until valid credentials are entered
    while [[ $SWITCH == N ]]; do
      #if netrc file non-existent
      if [[ $loginFILE == N ]]; then
        #ask user for input of username
        while true; do
          printf "\n\n"
          read -r -p "Please specify your NGS login [type n for exit]: " yn
          case $yn in
          [Nn]) exit ;;
          *)
            ngsUSER="$yn"
            break
            ;;
          esac
        done

        #ask user for input of password
        while true; do
          read -r -s -p "Please specify your NGS password [type n for exit]: " yn
          case $yn in
          [Nn]) exit ;;
          *)
            ngsPW="$yn"
            printf "\n"
            break
            ;;
          esac
        done

        #add username and password to the netrc file
        echo "machine ${SERVER} login $ngsUSER password $ngsPW" >~/.AP_login.txt
        chmod 600 ~/.AP_login.txt

        loginFILE=Y
      fi

      #test if username and password are correct
      TEST=$(curl --silent --head --netrc-file ~/.AP_login.txt https://${SERVER}/ | grep HTTP)

      if [[ $TEST == *"200"* ]] && [[ $loginFILE == Y ]]; then
        SWITCH=Y
      elif [[ $TEST == *"500 Internal Server Error"* ]]; then
        printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
        printf "\n  NGS server not available at the moment - please retry in some minutes  \n\n"
        printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n"
        exit 1
      elif [[ ! -z $TEST ]]; then
        printf "\n\nwrong username or password!!!\n\n"
        loginFILE=N
      else
        printf "\n"
        SWITCH="Y"
      fi
    done
  fi


  #---------------------------------------------------------------------------------------------------------
  #download NGS listing if NGS libraries are requested
  #remove -f statement once done
  if [[ $LISTING == YES ]]; then
    printf "\nplease wait while NGS-listing is getting downloaded - can take up to 1-2 minutes\n"
    curl --silent --retry 5 --retry-delay 300 --connect-timeout 240 --netrc-file ~/.AP_login.txt https://${SERVER}/ >${TMPtmp}listing.orig
    if [[ $AMERES == Y ]]; then
      curl --silent --retry 5 --retry-delay 300 --connect-timeout 240 --netrc-file ~/.AP_login.txt https://${SERVER}/ | tail -n +2 >>${TMPtmp}listing.orig
      curl --silent --retry 5 --retry-delay 300 --connect-timeout 240 --netrc-file ~/.AP_login.txt https://${SERVER}/ | tail -n +2 >>${TMPtmp}listing.orig
    fi

    #prepare backup of sequencedSamples file
    cp ${TMPtmp}listing.orig ${TMP}sequencedSamples.txt

    #cure sequencedSamples table into a useful format
    cat ${TMPtmp}listing.orig | tr ' ' '_' | sed 's/\t\t/\tempty\t/g;s/\t\t/\tempty\t/g;' |
      #added replacement as the barcode field changed
      sed 's/standard//g;s/auto//g;s/random//g'  | grep -v Failed >${TMPtmp}listing.txt

    if [[ ! -s ${TMPtmp}listing.txt ]]; then
      printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
      printf "\n  download of NGS index failed - please retry  \n\n"
      printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n"
      exit 1
    fi
  fi

  #---------------------------------------------------------------------------------------------------------
  #preset variables and files
  rm -rf ${TMPtmp}files.tmp
  rm -rf ${TMPtmp}URLs.txt
  rm -rf ${TMP}URLs.txt
  #test if all library files exist and report missing ones
  BAM=N
  exit_code1=0
  exit_code2=0
  exit_code3=0

  #---------------------------------------------------------------------------------------------------------
  #loop through lines in the input file and pre-process file-path if required

  for LINE in $(cat ${TMPtmp}files.txt | tr ' ' '@' | tr '\t' '@'); do

    #extract relevant variables from lib-file
    LINE=$(echo $LINE | tr '@' '\t')
    LIBRARY=$(echo $LINE | tr ' ' '\t' | awk '{print $1}')
    NAME=$(echo $LINE | tr ' ' '\t' | awk '{print $2}')
    
    #---------------------------------------------------------------------------------------------------------
    #create a local directory for current library
    currTMP=${TMPtmp}${NAME}/
    rm -rf $currTMP
    mkdir -p $currTMP

    SWITCH=N
    libPATH=
    NGS=N

    #---------------------------------------------------------------------------------------------------------
    #process line depending on type
    if [[ $LIBRARY == "/"* ]]; then
      #local file - no modification required
      #!printf "$NAME is a local lib-file - no conversion required\n"
      echo $LINE >>${TMPtmp}files.tmp
    elif [[ $LIBRARY == "http"* ]]; then
      #defined NGS files - no modification required
      #!printf "$NAME is a NGS-lib directly specified - no conversion required \n"
      echo $LINE >>${TMPtmp}files.tmp
    elif [[ $LIBRARY =~ ^(SRR|ERR|DRR|SRX) ]]; then
      #SRA library - create proper download link
      #?@ SRAdir=$(echo $LIBRARY | awk '{ print substr($1,1,6)}')
      #?@ SRAid=$(echo $LIBRARY | awk '{ print substr($1,4,length($1))}')
      #?@ SRAnum=$(echo $SRAid | 
      #?@   awk '{ 
      #?@     if(length($1)==8) {
      #?@     }else{
      #?@       if(length($1)==7) {
      #?@         NUM=substr($1,length($1),1)
      #?@         print "00" NUM
      #?@       }else{
      #?@         if(length($1)==6){
      #?@           NUM=substr($1,length($1)-1,2)
      #?@           print "0" NUM
      #?@         }else{
      #?@           print "weird SRA id= " $1 
      #?@           # TMP "error.txt
      #?@         }
      #?@       }
      #?@     fi
      #?@   }')
      #?@ SRAlink="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${SRAdir}/${SRAnum}/$LIBRARY/"
      #?@ echo $SRAlink $SRAdir $SRAnum $LIBRARY $SRAid 
      Xlink=$(curl --silent https://www.ebi.ac.uk/ena/browser/api/xml/${LIBRARY} | grep CDATA | head -n 1 | sed 's/ //g' | awk '{ split($1,X,/CDATA\[|\]/); print X[2]}' )
      SRAlink=$(curl --silent "$Xlink" |  tail -n +2 | sed 's/ftp\./http:\/\/ftp./g' | tr ' ' '\t' | cut -f 2 | tr ';' '~' )

      truncLINE=$(echo $LINE | tr ' ' '\t' |  cut --complement -f 1)
      #!printf "$NAME is a SRA-archive library - no conversion required \n"
      echo $LINE SRAlink=$SRAlink >>${TMPtmp}files.tmp
    elif [[ $LIBRARY == "NGS"* ]]; then
      #test if multiple libraries requested for mergine
      splitLIBs=$(echo $LINE | tr ' ' '\t' | cut -f 1 | tr ';' '\t')
      nLIBS=$(echo $splitLIBs | wc -w)

      for LIBRARY in $splitLIBs; do
        #---------------------------------------------------------------------------------------------------------
        #create a local directory for current library
        currTMP=${TMPtmp}${NAME}/
        rm -rf $currTMP
        mkdir -p $currTMP

        SWITCH=N
        libPATH=
        NGS=N

        #NGS library defined by ID - download meta-data and process for final -line

        #set BAM variable to enter adaptor selection
        BAM=Y
        #---------------------------------------------------------------------------------------------------------
        #from listing select the relevant lines for the current library
        LIBRARY=$(echo $LIBRARY | sed 's/NGS//')

        #create truncated line for printing
        truncLINE=$(echo $LINE | tr ' ' '\t' | awk '{ $1=""; $2=""; print }' | tr ' ' '\t' | tr '\t' '~')

        #determine if sRBC was supplied on the submission-line
        if [[ $truncLINE == *sRBC* ]]; then
          sRBC=$(echo $truncLINE |
            awk '{ for(i=1; i<=NF; i++){ if($i ~ "sRBC=") {split($i,splitsRBC,/=/); print splitsRBC[2];exit } } }')
        else
          sRBC=""
        fi
# -v ORS="\n"
        ##if only a single valid run and basecall is present - add to file
        awk -v LIBRARY=${LIBRARY} -v NAME=$NAME -v OFS=" " -v ORS="\t" -v TMP=${currTMP}/ -v SERVER=${SERVER} '
      BEGIN{
        nRUNs=0
        OUTcols[1]="readType"
        OUTcols[2]="flowcellId"
        OUTcols[3]="laneNr"
        OUTcols[4]="date"
        OUTcols[5]="result"
        OUTcols[6]="url"
        OUTcols[7]="barcodeString"
        OUTcols[8]="barcodes"
        OUTcols[9]="basecalls"
        OUTcols[10]="spikein"
        OUTcols[11]="SequencerType"
        OUTcols[12]="barcodeType"
        OUTcols[13]="sequencerType" 
      }
      {
        #read header into array
        if(NR==1) {
          for(i=1; i<=NF;  i++) {
            HEADER[$i]=i
          } 
        }else{
          if( $HEADER["sampleId"] == LIBRARY && $HEADER["result"] ~ "use:ok|use:accepted|use:combine|unchecked|discussing|NA" && ($HEADER["url"] ~ "http" || $HEADER["path"] ~ "tar.gz" ) ) {
            #create download link for tar.gz files
            if ( $HEADER["path"] ~ "tar.gz" && $HEADER["path"] ~ "_fqfull_" ){
              if($HEADER["path"] ~ "_"$HEADER["laneNr"]"_"){
                nPATH=split($HEADER["path"], pathArray, "/" )
                $HEADER["url"] = "https://" SERVER "/api/tarfile/file/" pathArray[nPATH]
              } 
            }

            RUNid=$HEADER["flowcellId"]"~"$HEADER["laneNr"]"~"$HEADER["basecalls"]"~"$HEADER["demux"]"~"$HEADER["md5"]
            print > TMP "all-lines_in_listing.txt"
            RUNs[RUNid]=$0
            nRUNs+=1
          }
                
        }
      }
      END{
        #if only one valid entry is present add path to easy.txt
        if ( nRUNs == 0 ) {
          
          #if no valid run present fill error-file
          print "no valid entry in NGS-listing for library "LIBRARY" - "NAME > TMP "no-entry.txt"
        
        } else {
          #print header to output file
          for(i=1; i<=length(OUTcols); i++) {
            print OUTcols[i]
          }
          print "\n"

          if ( nRUNs > 1 ){
            for (LANE in RUNs ) {
              n=split(RUNs[LANE], splitLANE, /\t| /)
              if ( splitLANE[HEADER["demux"]]=="false" || splitLANE[HEADER["url"]]~"tar.gz" ){
                for(i=1; i<=length(OUTcols); i++) {
                  print splitLANE[HEADER[OUTcols[i]]]
                  print splitLANE[HEADER[OUTcols[i]]] > TMP "multiple-entries.txt"
                }
                print "\n"
              }
            }
          }else{
            #print relevant values into output files for each entry
            for (LANE in RUNs ) {
              split(RUNs[LANE], splitLANE, /\t| /)
              for(i=1; i<=length(OUTcols); i++) {
                print splitLANE[HEADER[OUTcols[i]]]
                print splitLANE[HEADER[OUTcols[i]]] > TMP "single-entry.txt"
              }
              print "\n"
            }
          }
        }
      }
     ' ${TMPtmp}listing.txt | sed "s/^[ \t]*//" >${currTMP}selected_listing.txt
     
        #if no entry meets requirements print error message
        TEST=""
        TEST=$(tail -n +2  ${currTMP}selected_listing.txt )
        if [[ -z $altSERVER ]]; then
          SERVERsuggestion="try the official NGS server by adding the --official-ngs-server option"
        else
          SERVERsuggestion="try the default server by removing the --official-ngs-server option"
        fi

        if [[ -z $TEST ]]; then
          printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
          cat ${currTMP}no-entry.txt
          printf "library was not found - maybe try the $SERVERsuggestion \n"
          printf "at the moment it is not possible to combine libraries located on different servers \n"
          printf "If this did not help - try to re-run the script in some minutes. If it still fails --> ask Dominik for help\n"
          printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
          exit
        fi
        #---------------------------------------------------------------------------------------------------------
        #process selected listing
        nRUNs=$(cat ${currTMP}selected_listing.txt | wc -l)

        #if multiple entries were present for the ID - process through runs and ask user for selection

        if [[ $nRUNs -gt 2 ]]; then

          #ask user if libraries should be selected or if merging for all of them is required
          if [[ -z ${autoMERGE+x} ]]; then 
            while true; do
              printf "\n\nmultiple lane/demultiplex files have been detected for at least one library!\n"
              read -r -p "Do you want to automatically merge lane-files for all libraries [ y or n ]? " yn
              case $yn in
              [Nn])
                autoMERGE=N
                break
                ;;
              [Yy])
                autoMERGE=Y
                break
                ;;
              *) echo "please answer with y/Y or n/N" ;;
              esac
            done
          fi

          #create file for user run-selection
          awk -v LIBRARY=${LIBRARY} -v NAME=$NAME -v OFS="\t" -v TMP=${currTMP}/ '
          {
            #read header into array
            if(NR==1) {
              for(i=1; i<=NF;  i++) {
                HEADER[$i]=i
              } 
            }else{
              RUNid=$HEADER["flowcellId"]"~"$HEADER["laneNr"]"~"$HEADER["basecalls"]
              RUNs[RUNid]=$0
            }
          }
          END{
            for(LANE in RUNs) {
              N+=1
            }
            if(N == 1 ) {
              print LANE > TMP "selected-runs.txt"
            }else{
              print "#","readType","flowcellID", "laneNr","spikein","sequencingDate","sequencingResult", "basecalls" > TMP "run-selection.txt"

              for(LANE in RUNs) {
                split(RUNs[LANE], splitLANE, /\t/)
                print splitLANE[HEADER["readType"]],splitLANE[HEADER["flowcellId"]],splitLANE[HEADER["laneNr"]],splitLANE[HEADER["spikein"]],splitLANE[HEADER["date"]],splitLANE[HEADER["result"]],splitLANE[HEADER["basecalls"]]
              }
            }
          }' ${currTMP}selected_listing.txt | sort -k4,4 -k3,3n |
            awk -v OFS="\t" '
            {
              print NR,$0
            }' >>"${currTMP}run-selection.txt"

            #move selected-listing and replace by reduced variant
            mv ${currTMP}selected_listing.txt ${currTMP}selected_listing.old.txt

          if [[ $autoMERGE == N ]]; then
            
            #ask user to select one or more run-files for processing
            if [[ -f ${currTMP}run-selection.txt ]]; then
              printf "\nsample $LIBRARY-$NAME sequenced on multiple lanes - please select the ones you would like to analyze. [multiple choices possible as comma separated list]:\n"
              cat ${currTMP}run-selection.txt

              while true; do
                read -r -p "Please select the runs to analyze [ # or #,#,# ]" yn
                case $yn in
                [1-9] | [1-9][0-9])
                  selRUN=$(awk -v SEL=$yn '{if($1==SEL ) print $3"_"$4 }' ${currTMP}run-selection.txt)
                  if [[ -z $selRUN ]]; then
                    echo selected run does not exist!!! select valid entry
                  else
                    break
                  fi
                  ;;
                *)
                  echo $yn |
                    awk -v INFILE=${currTMP}run-selection.txt -v OUTFILE=${currTMP}selected-runs.txt -v ERRORFILE=${currTMP}error.txt '
                      BEGIN{
                        while((getline RAW_LINES < INFILE) > 0) {split(RAW_LINES,splitLINES,/\t| /); availSELECTION[splitLINES[1]]=splitLINES[3]"_"splitLINES[4] }
                        EXIT="Y"
                      } 
                      { 
                        N=split($1, X, /,/)
                        if( N==1) {
                          print  "\n\n\nplease supply runs in a valid comma separated format"
                        }else{
                          for( i in X){
                            if( X[i] in availSELECTION ) {
                              if(OUTSTRING==""){
                                OUTSTRING=availSELECTION[X[i]]
                              }else{
                                OUTSTRING=OUTSTRING" "availSELECTION[X[i]]
                              }
                            }else{
                              print "selected run #"X[i]" not available for selection" > ERRORFILE
                              EXIT=N
                            }
                          }
                          if(EXIT == "Y" ){
                            print OUTSTRING > OUTFILE 
                          }
                        }
                      }'

                  if [[ -f ${currTMP}selected-runs.txt ]]; then
                    selRUN=$(cat ${currTMP}selected-runs.txt)
                    echo $selRUN
                    break
                  elif [[ -f ${currTMP}error.txt ]]; then
                    printf "\n\n\n"
                    cat ${currTMP}error.txt
                    echo please select only existing run IDs
                    cat ${currTMP}run-selection.txt
                    rm -rf ${currTMP}error.txt
                  fi
                  ;;
                esac
              done
            else
              selRUN=$(cat ${currTMP}selected-runs.txt)
            fi
          else
            selRUN=$(
              mawk -v OFS="\t" '
              {
                if(NR>1){
                  if(X==""){
                    X=$3"_"$4
                  }else{
                    X=X" "$3"_"$4
                  }
                }
              }
              END{
                print X
              }' "${currTMP}run-selection.txt"
            )
          fi

          #filter old selected-listing for the new IDs
          selRUN=$(echo $selRUN | tr ' ' '~')

          mawk -v OFS="\t" -v selRUN=$selRUN '
          BEGIN{
            n=split(selRUN,splitRUNs,/~/)
            for(i in splitRUNs){ X[splitRUNs[i]]=1}
          }
          {
            if(NR==1 ){
              print
              for(i=1; i<=NF;  i++) {
                HEADER[$i]=i
              } 
            }else{
              if($HEADER["flowcellId"]"_"$HEADER["laneNr"] in X){
                print
              }
            }
          }' ${currTMP}selected_listing.old.txt > ${currTMP}selected_listing.txt


          #---------------------------------------------------------------------------------------------------------
          #test if only one basecall per run is present in the NGS-listing
          splitRUNid=$(echo $selRUN | tr '~' ' ')
          for RUNid in $splitRUNid; do
            rm -rf ${currTMP}LINE.txt
            rm -rf ${currTMP}selection.txt
            rm -rf ${currTMP}error.txt
            awk -v LIBRARY=${LIBRARY} -v NAME=$NAME -v OFS="\t" -v currTMP=${currTMP}/ -v TMP=$TMPtmp -v RUNid=$RUNid -v LIBRARY=$LIBRARY '
            BEGIN{
              COUNT=0
            }
            {
              #read header into array
              if(NR==1) {
                for(i=1; i<=NF;  i++) {
                  HEADER[$i]=i
                } 
              }else{
                split(RUNid,splitID,/_|~/)
                if ( $HEADER["flowcellId"] == splitID[1] && $HEADER["laneNr"] == splitID[2]){
                  COUNT+=1
                  FILE[COUNT]=$HEADER["flowcellId"]":!:"$HEADER["laneNr"]":!:"$HEADER["basecalls"]":!:"$HEADER["url"]
                  LINE[COUNT]=$0
                  URL=$HEADER["url"]
                  BC=$HEADER["barcodeString"]
                }
              }
            }
            END{
              if(COUNT>1){
                for(i=1; i<=COUNT; i++){
                  gsub(":!:", "\t", FILE[i])
                  print i,FILE[i] > currTMP "selection.txt"
                  print LINE[i] > currTMP "LINE.txt"
                }
              }           
            }' <( head -n 1 ${currTMP}selected_listing.txt && tail -n +2 ${currTMP}selected_listing.txt | sort -k2,2 -k3,3n ) 

            #if multiple basecalls are present throw error and exit
            if [[ -f ${currTMP}error.txt ]]; then
              cat ${currTMP}error.txt
              exit
            fi

            #if multiple basecalls are present ask user for selection
            if [[ -f ${currTMP}selection.txt ]]; then
              printf "\n\nmultiple basecalls are present for the selected run $RUNid - please select the one you would like to analyze.\n"
              cat ${currTMP}selection.txt
              while true; do
                read -r -p "Please select the basecall to analyze [ # ]: " yn
                case $yn in
                [1-9] | [1-9][0-9])
                  selRUN=$(awk -v SEL=$yn '{if($1==SEL ) print $2 }' ${currTMP}selection.txt)
                  if [[ -z $selRUN ]]; then
                    echo selected run does not exist!!! select valid entry
                  else
                    #remove not selected lines from the selected_listing.txt
                    mawk -v OFS="\t" -v LINES=${currTMP}LINE.txt -v selRUN=$selRUN -v SEL=$yn '
                    BEGIN{
                      X=1
                      while((getline RAW_LINES < LINES) > 0) {
                        
                        if(X != SEL){
                          REMOVE[RAW_LINES]=1
                        }
                        X++
                      }
                    }
                    {

                      if($0 in REMOVE){
                        #do not print
                      }else{
                        print
                      }
                    }' ${currTMP}selected_listing.txt > ${currTMP}selected_listing.tmp.txt
                    
                    mv ${currTMP}selected_listing.tmp.txt ${currTMP}selected_listing.txt
                    break
                  fi
                  ;;
                *)
                  echo selected run does not exist!!! select valid entry
                  ;;
                esac
              done
            fi
          done
        fi

        #extract URL and correct BC information from the selected listing
        awk -v OFS="\t" -v LIBRARY=$LIBRARY -v NAME=$NAME -v truncLINE=$truncLINE -v TMP=$TMPtmp -v sRBC=$sRBC '
        BEGIN{
          gsub("~", "\t", truncLINE)

          #test if BC is provided on the submission line  and use this instead of NGS barcode
          n=split(truncLINE,splitLINE,/\t/)
          for(i=1;i<=n;i++){
            if(splitLINE[i]~"^BC="){
              extBCi7=splitLINE[i]
              sub(extBCi7,"",truncLINE)
              sub("BC=","",extBCi7)
            }
            if(splitLINE[i]~"^BCi7="){
              extBCi7=splitLINE[i]
              sub(extBCi7,"",truncLINE)
              sub("BCi7=","",extBCi7)
            }
            if(splitLINE[i]~"^BCi5="){
              extBCi5=splitLINE[i]
              sub(extBCi5,"",truncLINE)
              sub("BCi5=","",extBCi5)
            }
            if(splitLINE[i]~"^sRBC="){
              extsRBC=splitLINE[i]
              sub(extsRBC,"",truncLINE)
              sub("sRBC=","",extsRBC)
            }
          }

          #reverse complementing function 
          c["A"] = "T"; c["C"] = "G"; c["G"] = "C"; c["T"] = "A" 
        }
        function revcomp(x,  i, o) {
          o = ""
          for(i = length; i > 0; i--){
            o = o c[substr(x, i, 1)] }
          return(o)

        }
        {
          #read header into array
          if(NR==1) {
            for(i=1; i<=NF;  i++) {
              HEADER[$i]=i
            } 
          }else{
            nBC=split($HEADER["barcodes"], splitBC,/,/)
            #go through all BCstrings to find the relevant for this library
            for(n=1; n<=nBC;n++) {
              if(splitBC[n]~LIBRARY){

                #split name and extract BC
                split(splitBC[n], fullBC, /:/)
                BCi5=fullBC[2]

                #remove barcodeString (1st BC) from the full BC and set BC variable
                sub($HEADER["barcodeString"], "", BCi5)

                #set i7 barcode
                BCi7=$HEADER["barcodeString"]

                #if nothing left only print  BCi7; otherwise process further
                if(BCi5!=""){
                  #test if sRBC present
                    #split on i to delimit sRBC
                    delete i5split
                      split(BCi5, i5split, "i")
                      BCi5=i5split[1]
                      sRBC=i5split[2]                      


                  #reverse complement 2nd barcode if sequenced on Nova or NextSeq
                  if( BCi5 != "") {
                    if($HEADER["SequencerType"] ~ "NextSeq" || $HEADER["SequencerType"] ~ "NovaSeq" ) {
                      BCi5=revcomp(BCi5)
                    }
                  }
                }

                if(extBCi7 != "" ){
                  BCi7=extBCi7
                }
                if(extBCi5 != "" ){
                  if (extBCi5 == "NA") {
                    BCi5=""
                  }else{
                    BCi5=extBCi5
                  }
                }
                if(extsRBC != "" ){
                  sRBC=extsRBC
                }
              }
            }

            #remove poteintial UMIs from the BC Tag

            #remove poteintial UMIs from the BC Tag
            nN_BCi7=gsub("N","",BCi7)
            nN_BCi5=gsub("N","",BCi5)
            gsub("N","",sRBC)
            
            if(nN_BCi7 > 0 || nN_BCi5 > 0){
              UMI="\tUMIi7="nN_BCi7"\tUMIi5="nN_BCi5
            }else{
              UMI=""
            }

            print "NGS"LIBRARY,NAME, truncLINE, "BCi7="BCi7, "BCi5="BCi5, "sRBC="sRBC, "selRUNs="$HEADER["flowcellId"]"_"$HEADER["laneNr"]"_"$HEADER["basecalls"],$HEADER["barcodeType"] UMI
            print $HEADER["url"], LIBRARY, BCi7":"BCi5"i"sRBC"~"$HEADER["flowcellId"]"_"$HEADER["laneNr"]"_"$HEADER["basecalls"],$HEADER["barcodeType"] >> TMP "URLs.txt"
          }
        }' ${currTMP}selected_listing.txt | tr -s "\t" >>${TMPtmp}files.tmp

      done
      
    elif [[ $LIBRARY == "TGZ"* ]]; then
      #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if [[ $LINE != *"TGZ="* ]]; then
        printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
        printf "\n  no NGS tar.gz file specified for $NAME   \n\n"
        printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n"
        exit 1
      fi
        #set BAM variable to enter adaptor selection
        BAM=Y

      LIBRARY=$(echo $LIBRARY | sed 's/TGZ//')

      truncLINE=$(echo $LINE | tr ' ' '\t' |  cut --complement -f 1-2 | tr '\t' "\n" | grep -v TGZ= | tr '\n' '\t' )

      TGZfull=$(echo $LINE | tr ' ' '\t' | tr '\t' '\n'  | grep "TGZ=" )
      splitTGZ=$(echo $TGZfull | awk -v OFS="\t" '{ 
        split($1,splitTGZ,"=")
        gsub(";","\t",splitTGZ[2])
        print splitTGZ[2]
      }')
      
      BCtype=$(echo $LINE | awk -v OFS="\t"  '{for(i=1;i<=NF;i++){ if($i ~ "BCtype=") {split($i,splitBCtype,"="); print splitBCtype[2]}}}')
      if [[ -z $BCtype ]]; then
        BCtype="standard"
      fi

      for currTGZ in $splitTGZ; do
        FC=$(echo $currTGZ | awk -v OFS="\t"  '{n=split($1,splitFC,"/"); split(splitFC[n],splitFC2,"_"); print splitFC2[1] }')
        LANE=$(echo $currTGZ | awk -v OFS="\t"  '{n=split($1,splitFC,"/"); split(splitFC[n],splitFC2,"_"); print splitFC2[2] }')

        echo $FC $LANE

        #print URL and correct BC information from the information available listing
        awk -v OFS="\t" -v LIBRARY=$LIBRARY -v NAME=$NAME -v truncLINE="$truncLINE"  -v currTGZ="$currTGZ" -v TMP=$TMPtmp   -v FC=$FC -v LANE=$LANE -v BCtype=$BCtype '
        BEGIN{
          gsub("~", "\t", truncLINE)

          #test if BC is provided on the submission line  and use this instead of NGS barcode
          n=split(truncLINE,splitLINE,/\t/)
          for(i=1;i<=n;i++){
            if(splitLINE[i]~"^BC="){
              extBCi7=splitLINE[i]
              sub(extBCi7,"",truncLINE)
              sub("BC=","",extBCi7)
            }
            if(splitLINE[i]~"^BCi7="){
              extBCi7=splitLINE[i]
              sub(extBCi7,"",truncLINE)
              sub("BCi7=","",extBCi7)
            }
            if(splitLINE[i]~"^BCi5="){
              extBCi5=splitLINE[i]
              sub(extBCi5,"",truncLINE)
              sub("BCi5=","",extBCi5)
            }
            if(splitLINE[i]~"^sRBC="){
              extsRBC=splitLINE[i]
              sub(extsRBC,"",truncLINE)
              sub("sRBC=","",extsRBC)
            }
          }

          #reverse complementing function 
          c["A"] = "T"; c["C"] = "G"; c["G"] = "C"; c["T"] = "A" 
        }
        function revcomp(x,  i, o) {
          o = ""
          for(i = length; i > 0; i--){
            o = o c[substr(x, i, 1)] }
          return(o)
        }
        {
          #read header into array

            nBC=split(extBCi7, splitBC,/,/)
            #go through all BCstrings to find the relevant for this library
            for(n=1; n<=nBC;n++) {

                if(extBCi7 != "" ){
                  BCi7=extBCi7
                }
                if(extBCi5 != "" ){
                  BCi5=extBCi5
                }
                if(extsRBC != "" ){
                  sRBC=extsRBC
                }
              
            }

            #remove poteintial UMIs from the BC Tag

            #remove poteintial UMIs from the BC Tag
            nN_BCi7=gsub("N","",BCi7)
            nN_BCi5=gsub("N","",BCi5)
            gsub("N","",sRBC)
            
            if(nN_BCi7 > 0 || nN_BCi5 > 0){
              UMI=",UMIi7="nN_BCi7",UMIi5="nN_BCi5
            }else{
              UMI=""
            }

            print "TGZ"LIBRARY,NAME, truncLINE, "BCi7="BCi7, "BCi5="BCi5, "sRBC="sRBC, "selRUNs="FC"_"LANE"_NA",BCtype UMI
            print currTGZ, LIBRARY, BCi7":"BCi5"i"sRBC"~"FC"_"LANE"_NA",BCtype >> TMP "URLs.txt"
          
        }' <(echo a) | tr -s "\t" >>${TMPtmp}files.tmp
      done
      #test if multiple libraries requested for merging






      #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    else
      echo line= $LINE
      printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
      printf "\n  file specified for library $NAME is not supplied in a supported format   \n\n"
      printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n"
      exit 1
    fi

    #remove potential double lines due to lib merging
    awk -v OFS="\t" '
    {
      Xcount[$2]+=1
      X[$2][Xcount[$2]]=$0
    }
    END{
      for(LIB in X ){
        if (Xcount[LIB]>1){
          for(i=1; i<=Xcount[LIB]; i++){
            n=split(X[LIB][i],splitLINE,/\t| /)
            NGS=NGS splitLINE[1]
            for(j=3; j<=n; j++){
              if(splitLINE[j] ~ "selRUNs="){
                split(splitLINE[j], splitRUNs, /=/)
                selRUNs=selRUNs splitRUNs[2]
              }else{
                if(splitLINE[j] ~ "BCi7="){
                  split(splitLINE[j], splitRUNs, /=/)
                  BCi7s=BCi7s splitRUNs[2]
                }
                if(splitLINE[j] ~ "BCi5="){
                  split(splitLINE[j], splitRUNs, /=/)
                  BCi5s=BCi5s splitRUNs[2]
                }
                if(splitLINE[j] ~ "sRBC="){
                  split(splitLINE[j], splitRUNs, /=/)
                  sRBCs=sRBCs splitRUNs[2]
                }
                if(splitLINE[j] !~"BCi7" && splitLINE[j] !~"BCi5" && splitLINE[j] !~"sRBC" ){
                  if(resLINE !~ splitLINE[j] ){
                    resLINE=resLINE"\t"splitLINE[j]
                  }
                }
              }
            }
            if(i<Xcount[LIB]){
              NGS=NGS ";"
              selRUNs=selRUNs";"
              BCi7s=BCi7s";"
              BCi5s=BCi5s";"
              sRBCs=sRBCs";"
            }
          }

          #remove poteintial UMIs from the BC Tag

            #remove poteintial UMIs from the BC Tag
            nN_BCi7=gsub("N","",BCi7)
            nN_BCi5=gsub("N","",BCi5)
            gsub("N","",sRBC)
            
            if(nN_BCi7 > 0 || nN_BCi5 > 0){
              UMI=",UMIi7="nN_BCi7",UMIi5="nN_BCi5
            }else{
              UMI=""
            }

          print NGS,splitLINE[2],"BCi7="BCi7s,"BCi5="BCi5s,"sRBC="sRBCs,"selRUNs="selRUNs,resLINE UMI
          NGS=""
          selRUNs=""
          BCi7s=""
          BCi5s=""
          sRBCs=""
          resLINE=""
        }else{
          print X[LIB][1]
        }
      }
    }
    ' ${TMPtmp}files.tmp >${TMP}files.txt

    #reset FILE variable to real file
    FILE_CONTAINING_LIBRARIES=${TMP}files.txt

    #---------------------------------------------------------------------------------------------------------
    #test if supplied color is in correct format
    if [[ $LINE == *"COLOR="* ]]; then
      testCOLOR=$(
        echo $LINE |
          awk '
      {
        for(i=1; i<=NF; i++){
          if($i ~ "COLOR=") {
            split($i,splitCOLOR,/=/)
            print splitCOLOR[2]
            exit
          }
        }
      }'
      )

      STATUS=$(echo $testCOLOR | awk -v FS="," '{
      if(NF != 3) print "nOK"; else {
        for (i=1; i<=3; i++) { if ($i <0 || $i >256) print "xOK" }
      }
      }')
      if [[ ! -z ${STATUS} ]]; then
        printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
        printf "\n  provided color for library $NAME in wrong format \n"
        printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n"
        exit 1
      fi

    fi
    #---------------------------------------------------------------------------------------------------------
    #test if supplied subsampling value is numeric
    if [[ $LINE == *"SUB="* ]]; then
      testSUB=$(
        echo $LINE |
          awk '
      {
        for(i=1; i<=NF; i++){
          if($i ~ "SUB=") {
            split($i,splitSUB,/=/)
            print splitSUB[2]
            exit
          }
        }
      }'
      )

      TESTcond='^[0-9]+$'
      if ! [[ $testSUB =~ $TESTcond ]]; then
        usage
        printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
        printf "  provide numeric value for subsampling for library $NAME \n"
        printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n"
        exit 1
      fi
    fi
  done

  #---------------------------------------------------------------------------------------------------------
  #---------------------------------------------------------------------------------------------------------
  #test if all libraries are available

  printf "\n\nplease wait while all libraries are tested for availability\n"
  printf "if SRA libraries are requestet this will take up to 20 seconds per library\n\n"

  #test files in the URL file if NGS libraries present in run
  if [[ -f ${TMPtmp}URLs.txt ]]; then
    #test if library already available --> if so remove from URLs file
    rm -rf ${TMPtmp}filteredURLs.txt
    if [[ $FORCEimport != Y ]]; then
      while read currURL currSAMPLE runID BCtype; do
        if [[ ! -f ${LIB_STORAGE_FOLDER}${currSAMPLE}~${runID}.bam && ! -f ${LIB_STORAGE_FOLDER}${currSAMPLE}~${runID}.bam.tmp ]]; then
          echo $currURL $currSAMPLE $runID $BCtype >>${TMPtmp}filteredURLs.txt
        fi
      done <${TMPtmp}URLs.txt
    else
      cp ${TMPtmp}URLs.txt ${TMPtmp}filteredURLs.txt
      while read URL currSAMPLE runID BCtype; do
        rm -rf ${LIB_STORAGE_FOLDER}${currSAMPLE}~${runID}.bam*
      done <${TMPtmp}filteredURLs.txt
    fi

    #collapse URLs file to single entry per URL
    if [[ -f ${TMPtmp}filteredURLs.txt ]]; then
      awk -v OFS="\t" -v SERVER=$SERVER '
      {
        if($3~":"){
          DUAL="Y"
        }else{
          DUAL="N"
        }
        split($4,splitBCtype,/_cat_/)
        X[$1][DUAL][splitBCtype[1]][$2]=$3
      }
      END{
        for(i in X) {
          for(DUAL in X[i]){
            IDs=""
            RUNs=""
            SWITCH="N"
            for(BCtype in X[i][DUAL]){
              for(j in X[i][DUAL][BCtype]) {
                if(SWITCH == "N"){
                  IDs=j
                  RUNs=X[i][DUAL][BCtype][j]
                  SWITCH="Y"
                }else{
                  IDs=IDs ";" j
                  RUNs=RUNs ";" X[i][DUAL][BCtype][j]
                }
              }
              #change URL to server specified in the header
              ##required as jenking listing contains gecko files
              
              newURL=i

              print newURL,IDs,RUNs,BCtype
            }
          }
        }
      }' ${TMPtmp}filteredURLs.txt >${TMPtmp}collapsedURLs.txt

      #test if all requested files are available @ NGS

      while read currURL currSAMPLES runID BCtype; do
        if [[ $currURL == *https* ]]; then
          echo $currURL
          if curl --netrc-file ~/.AP_login.txt --output /dev/null --silent --head --fail ${currURL}; then
            x=a
          else
            failedSAMPLES=$(echo $currSAMPLES | sed 's/;/\t/g')
            printf "\nfile for libraries $failedSAMPLES is not available on the NGS portal \n"
            exit_code1=1
          fi
        elif [[ $currURL == *tar.gz* ]]; then
          echo $currURL
          if [[ -s $currURL ]]; then
            x=a
          else
            failedSAMPLES=$(echo $currSAMPLES | sed 's/;/\t/g')
            printf "\nfile for libraries $failedSAMPLES is not available on the NGS portal \n"
            exit_code1=1
          fi
        else
          echo non-defined condition
          exit
        fi
      done <${TMPtmp}collapsedURLs.txt

      #move final URL file to the TMP directory
      mv ${TMPtmp}collapsedURLs.txt ${TMP}URLs.txt
    fi
  fi

  cp ${TMP}files.txt ${TMP}files.txt.old
  while read LINE; do
    if [[ $LINE != *selRUNs=* ]]; then #test if all libraries are available - only non-NGS file submissions
      LIBRARY=$(echo $LINE | tr ' ' '\t' | awk '{print $1}')

      #split libraries for merging into individual paths
      multiLIBRARY=$(echo $LIBRARY | tr ';' '\t')
      #determine n of libraries to merge
      NcurrLIBRARY=$(echo $multiLIBRARY | wc -w)

      for currLIBRARY in ${multiLIBRARY}; do
        if [[ $currLIBRARY =~ (SRR|ERR|DRR|SRX) ]] && [[ $LINE == *"sra.ebi"* ]]; then #test for availability of SRA libraries
          #?@ SRAlinkNEW=$(echo $LINE | 
          #?@ awk '{
          #?@   for(i=1; i<=NF; i++){
          #?@     if($i~"sra.ebi" && $i~"SRAlink=") COLUMN=i
          #?@   }
          #?@   sub("SRAlink=","",$COLUMN)
          #?@   print $COLUMN
          #?@ }' | tr '~' '\t')
          #?@ echo $SRAlinkNEW
          #?@ if ncftpls $SRAlink/* > ${TMP}listing.txt; then 
          #?@   echo
          #?@   newLINK=$(awk -v SRAlink=$SRAlink '{ print SRAlink $1}' ${TMP}listing.txt | tr '\n' '~'| sed 's/~$//')
          #?@   echo $newLINK
          #?@   awk -v SRAlink=$SRAlink -v newLINK=$newLINK '{
          #?@     if($0~SRAlink){
          #?@       for(i=1; i<=NF; i++){
          #?@         if($i~SRAlink){
          #?@           $i="SRAlink="newLINK
          #?@         }
          #?@       }
          #?@     }
          #?@     print
          #?@   }' ${TMP}files.txt > ${TMP}files.tmp
          #?@   cat ${TMP}files.tmp
          #?@   mv ${TMP}files.tmp ${TMP}files.txt
          #?@ else
          #?@   printf "\n sample $currLIBRARY does not exist @ ENA\n" 
          #?@   exit_code1=1
          #?@ fi
          BAM=Y

        else #test availability for all other types of libraries (bam/fasta/fastq)
          #set switch if at least one bam file is available
          if [[ ${currLIBRARY} == *.bam ]]; then BAM=Y; fi
          #report missing libraries
          if [ ! -f "$currLIBRARY" ]; then
            printf "\n file $currLIBRARY does not exist"
            exit_code1=1
          else
            #if libraries get merged test if both are bam
            if [[ ${NcurrLIBRARY} -gt 1 ]]; then

              #!removed test as singularity not available on startup host for peter
              #@ samtools view ${currLIBRARY} 2>${TMPtmp}test.txt | head >/dev/null
              if [[ -s ${TMPtmp}test.txt ]]; then
                printf "\n file $currLIBRARY not a bam file"
                exit_code2=2
              fi
            fi
          fi
        fi
      done
    fi
  done <${TMP}files.txt.old

  if [[ $exit_code1 -eq 1 ]] || [[ $exit_code2 -eq 2 ]] || [[ $exit_code3 -eq 3 ]]; then
    printf "\n\n\n"
    printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  fi

  #if error exists report message and exit script
  if [[ $exit_code1 -eq 1 ]]; then
    printf "\n         not all libraries are avilable"
  fi
  if [[ $exit_code2 -eq 2 ]]; then
    printf "\n         only bam files allowed for merging"
  fi
  if [[ $exit_code3 -eq 3 ]]; then
    printf "\n         no barcode supplied for all NGS-files"
  fi

  if [[ $exit_code1 -eq 1 ]] || [[ $exit_code2 -eq 2 ]] || [[ $exit_code3 -eq 3 ]]; then
    printf "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
    printf "\n\n\n"
    exit 1
  fi

  #evaluate the number of libraries
  nFILES=$(cat "${TMPtmp}files.txt" | wc -l)

  #---------------------------------------------------------------------------------------------------------
  #---------------------------------------------------------------------------------------------------------
  #setup DGE file if requested
  #ask if user wants to perform DGE analysis
  if [[ $TYPE == RNAseq ]] && [[ $DGE != Y ]]; then
    printf "\n\n"
    while true; do
      read -r -p "Do you want to perform a differential gene expression (DGE) analysis (replicates required) [ y or n ]? " yn
      case $yn in
      [Nn])
        DGE=N
        break
        ;;
      [Yy])
        DGE=Y
        break
        ;;
      *) echo "please answer with y/Y or n/N" ;;
      esac
    done
  fi

  if [[ $DGE == Y ]]; then
    dgeSWITCH=N
    dgeCHECK=N
    currNR=1

    newDGE=Y
    if [[ -e ${FOLDER}DGE/DGE_info.txt ]]; then
      cat ${FOLDER}DGE/DGE_info.txt | tr ' ' '\t' | cut -f 1-2
      while true; do
        read -r -p "Do you want to keep the DGE settings listed above? [y/n]: " name
        case $name in
        [Y,y])
          newDGE=N
          break
          ;;
        [N,n])
          rm -rf ${FOLDER}DGE/DGE_info.txt
          newDGE=Y
          break
          ;;
        *) echo "Please answer with Y or N!" ;;
        esac
      done
    fi

    if [[ $newDGE == Y ]]; then
      rm -rf ${TMP}DGE_info.txt

      awk '
      {
        print NR,$2
      }' <(sort -k2,2 ${TMPtmp}files.txt) | tr ' ' '\t' >${TMPtmp}DGElibs.txt

      cp ${TMPtmp}DGElibs.txt ${TMPtmp}DGElibs.tmp

      while [[ $dgeSWITCH == N ]]; do

        #determine name for genotype from user
        if [[ $currNR -eq 1 ]]; then TEXT="reference"; else TEXT=$currNR; fi
        while true; do
          cat ${TMPtmp}DGElibs.tmp
          read -r -p "Please provide a name for the $TEXT genotype or type N to continue analysis without adding more genotypes: " name
          case $name in
          [N,n]) if [[ $currNR -le 2 ]]; then
            while true; do
              read -r -p "not enough genotypes defined for DGE! continue without DGE analysis?: " yn
              case $yn in
              [Y,y])
                dgeCHECK="Y"
                break
                ;;
              [N,n]) break ;;
              *) echo "Please answer with Y or N!" ;;
              esac
            done
            if [[ $dgeCHECK == Y ]]; then break; fi
          else
            dgeCHECK="Y"
            break
          fi ;;
          "") echo "No input provided! Please specify either a genotype name or cancel with N!" ;;
          *)
            TEST=$(echo $name | sed 's/[A-Za-z0-9._]//g')
            if [[ -z ${TEST} ]]; then
              GENO=$name
              break
            else
              echo "Please use only letters, numbers, _ and . !"
            fi
            ;;
          esac
        done
        currNR=$((currNR + 1))

        #assign libraries to genotype
        if [[ $dgeCHECK == N ]]; then
          while true; do

            read -r -p "Please select at least 2 libraries for genotype $GENO [#,#]: " selLIBs
            case $selLIBs in
            *)
              TEST=$(echo $selLIBs | sed 's/[0-9,]//g')
              if [[ ! -z $TEST ]]; then
                echo "Please use only numbers!"
              else
                #reset files
                rm -rf ${TMPtmp}tmp_sel.txt
                rm -rf ${TMPtmp}tmp.err
                rm -rf ${TMPtmp}no-break.tmp

                #test entry for validity and add libaries to DGE file if possible
                awk -v selLIB=$selLIBs -v OFS="\t" -v TMP=$TMPtmp -v GENO=$GENO -v FOLDER=${FOLDER}DGE/raw/ '
                    BEGIN{
                      n=split(selLIB, selLIBsplit, /,/)
                      if( n == 1 ) {
                        print "\n\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
                        print "please select more than 1 library - replicates are required for DGE analysis!"
                        print "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n\n"
                        print "xy" > TMP "no-break.tmp"
                        exit_invoked=1
                        exit
                      }
                      for(i in selLIBsplit) {
                        if( selLIBarray[selLIBsplit[i]] == "IN" ){
                          print "\n\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
                          print "duplicated selection "selLIBsplit[i]"! please select valid replicates for DGE analysis"
                          print "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n\n"
                          print "xy" > TMP "no-break.tmp"
                          exit_invoked=1
                          exit
                        }
                        selLIBarray[selLIBsplit[i]]="IN"
                      }
                    }
                    {
                   
                      if($1 in selLIBarray) {
                        print $2, GENO, FOLDER $2"/quant.sf" > TMP "sel_for_DGE.txt"
                        delete selLIBarray[$1]
                      }else{
                        print $0 > TMP "residual_libs.tmp"
                      }
                    }
                    END{
                      if ( length(selLIBarray) > 0 && exit_invoked != 1) {
                        print "\n\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
                        for(i in selLIBarray){
                          print i > TMP "tmp.err"
                          print "selection " i " does not exist in the library index! Please provide a valid entry"
                        }
                        print "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n\n"
                      }
                    }' ${TMPtmp}DGElibs.tmp

              fi
              if [[ -f ${TMPtmp}tmp.err ]]; then
                BADlibs=$(cat ${TMPtmp}tmp.err | tr '\n' ' ')
              elif [[ ! -f ${TMPtmp}no-break.tmp ]]; then
                cat ${TMPtmp}sel_for_DGE.txt >>${TMP}DGE_info.txt
                if [[ ! -f ${TMPtmp}residual_libs.tmp ]]; then nLIBSres=0; else nLIBSres=$(cat ${TMPtmp}residual_libs.tmp | wc -l); fi
                if [[ $nLIBSres -le 1 ]]; then
                  printf "\n\n not enough libraries left for another genotype!\n\n"
                  dgeCHECK="Y"

                else
                  mv ${TMPtmp}residual_libs.tmp ${TMPtmp}DGElibs.tmp
                fi
                break
              fi
              ;;
            esac
          done
        fi

        if [[ $dgeCHECK == Y ]]; then
          #print DGE_info file
          printf "\n\nLIBRARY \t\t\t GENOTYPE\n"
          cat ${TMP}DGE_info.txt | cut -f 1-2 | sed 's/\t/\t\t\t/'

          #ask user if it is fine, otherwise repeat
          while true; do
            read -r -p "\nIs the sample to genotype assignment correct? If no, slection starts again! [y/n]: " yn
            case $yn in
            [Y,y])
              dgeSWITCH="Y"
              break
              ;;
            [N,n])
              dgeCHECK="N"
              currNR=1
              rm -rf ${TMP}DGE_info.txt
              cp ${TMPtmp}DGElibs.txt ${TMPtmp}DGElibs.tmp
              break
              ;;
            *) echo "Please answer with Y or N!" ;;
            esac
          done

        fi

      done
    fi
  fi


  #---------------------------------------------------------------------------------------------------------
  #---------------------------------------------------------------------------------------------------------
  #setup CHIP ratio tracks
  #ask if user wants to perform DGE analysis
  RATIOtracks=N
  if [[ $TYPE == CHIPseq ]]; then
    printf "\n\n"
    while true; do
      read -r -p "Do you want to generate genome browser ratio tracks [ y or n ]? " yn
      case $yn in
      [Nn])
        RATIOtracks=N
        break
        ;;
      [Yy])
        RATIOtracks=Y
        break
        ;;
      *) echo "please answer with y/Y or n/N" ;;
      esac
    done
  fi

  if [[ $RATIOtracks == Y ]]; then
    ratioSWITCH=N
    ratioCHECK=N
    currNR=1

    newRATIO=Y
    if [[ -e ${TMP}RATIO.final.txt  ]]; then
      cat ${TMP}RATIO.final.txt | tr ' ' '\t' | cut -f 1-2
      while true; do
        read -r -p "Do you want to keep the RATIO track settings listed above? [y/n]: " name
        case $name in
        [Y,y])
          newRATIO=N
          break
          ;;
        [N,n])
          rm -rf ${TMP}RATIO.final.txt
          newDGE=Y
          break
          ;;
        *) echo "Please answer with Y or N!" ;;
        esac
      done
    fi

    if [[ $newRATIO == Y ]]; then
      rm -rf ${TMP}RATIO_info.txt
      rm -rf ${TMPtmp}RATIO.final.txt

      awk '
      {
        print NR,$2
      }' <(sort -k2,2 ${TMPtmp}files.txt) | tr ' ' '\t' >${TMPtmp}RATIOlibs.txt

      cp ${TMPtmp}RATIOlibs.txt ${TMPtmp}RATIOlibs.tmp

      while [[ $ratioSWITCH == N ]]; do

        #determine name for genotype from user
        while true; do
          cat ${TMPtmp}RATIOlibs.tmp
          read -r -p "Please select a CHIP input library [provide number in list] or cancel with N: " INPUT
            case $INPUT in
            [N,n]) 
              ratioCHECK="Y"
              break;;
            "") echo "No input provided! Please specify either a input library  or cancel with N!" ;;
            *)
              TEST=$(echo $INPUT | sed 's/[0-9]//g')
              if [[ ! -z $TEST ]]; then
                echo "Please use only numbers!"
              else
                #reset files
                rm -rf ${TMPtmp}tmp.err

                awk -v OFS="\t"  -v INPUT=$INPUT -v TMP=$TMPtmp '{
                  if($1==INPUT){
                    INPUTlib=$2
                  }
                }
                END{
                  if(INPUTlib == "" ){
                    print "\n\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
                      print INPUT > TMP "tmp.err"
                      print "selection " INPUT " does not exist in the library index! Please provide a valid entry"
                    print "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n\n"
                  }else{
                    print INPUTlib > TMP "currINPUT.txt"
                  }
                }' ${TMPtmp}RATIOlibs.tmp

              fi
              if [[ -f ${TMPtmp}tmp.err ]]; then
                BADlibs=$(cat ${TMPtmp}tmp.err | tr '\n' ' ')
              else
                INPUTlib=$(cat ${TMPtmp}currINPUT.txt )
                grep -v $INPUTlib ${TMPtmp}RATIOlibs.tmp > ${TMPtmp}residual_libs.tmp
                mv ${TMPtmp}residual_libs.tmp ${TMPtmp}RATIOlibs.tmp
                break
              fi
            esac
        done

        #assign libraries to genotype
        if [[ $ratioCHECK == N ]]; then
          while true; do
            cat ${TMPtmp}RATIOlibs.tmp
            read -r -p "Please select at libraries for ratio calculation with library $INPUTlib [#,#]: " selLIBs
            case $selLIBs in
            *)
              TEST=$(echo $selLIBs | sed 's/[0-9,]//g')
              if [[ ! -z $TEST ]]; then
                echo "Please use only numbers!"
              else
                #reset files
                rm -rf ${TMPtmp}tmp_sel.txt
                rm -rf ${TMPtmp}tmp.err
                rm -rf ${TMPtmp}no-break.tmp

                #test entry for validity and add libaries to DGE file if possible
                awk -v selLIB=$selLIBs -v OFS="\t" -v TMP=$TMPtmp -v INPUT=$INPUTlib  '
                    BEGIN{
                      n=split(selLIB, selLIBsplit, /,/)
                      for(i in selLIBsplit) {
                        if( selLIBarray[selLIBsplit[i]] == "IN" ){
                          print "\n\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
                          print "duplicated selection "selLIBsplit[i]"! please select libraries onlyonce "
                          print "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n\n"
                          print "xy" > TMP "no-break.tmp"
                          exit_invoked=1
                          exit
                        }
                        selLIBarray[selLIBsplit[i]]="IN"
                      }
                    }
                    {
                   
                      if($1 in selLIBarray) {
                        if(OUTlib==""){
                          OUTlib=$2
                        }else{
                          OUTlib=$2":!!:"OUTlib
                        }
                        delete selLIBarray[$1]
                      }else{
                        print $0 > TMP "residual_libs.tmp"
                      }
                    }
                    END{
                      if ( length(selLIBarray) > 0 && exit_invoked != 1) {
                        print "\n\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
                        for(i in selLIBarray){
                          print i > TMP "tmp.err"
                          print "selection " i " does not exist in the library index! Please provide a valid entry"
                        }
                        print "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n\n"
                      }else{
                        print INPUT, OUTlib >> TMP "RATIO.final.txt"
                      }
                    }' ${TMPtmp}RATIOlibs.tmp

              fi
              if [[ -f ${TMPtmp}tmp.err ]]; then
                BADlibs=$(cat ${TMPtmp}tmp.err | tr '\n' ' ')
              elif [[ ! -f ${TMPtmp}no-break.tmp ]]; then
                if [[ ! -f ${TMPtmp}residual_libs.tmp ]]; then nLIBSres=0; else nLIBSres=$(cat ${TMPtmp}residual_libs.tmp | wc -l); fi
                if [[ $nLIBSres -le 1 ]]; then
                  printf "\n\n not enough libraries left for more ratio-tracks!\n\n"
                  ratioCHECK="Y"

                else
                  mv ${TMPtmp}residual_libs.tmp ${TMPtmp}RATIOlibs.tmp
                fi
                break
              fi
              ;;
            esac
          done
        fi

        if [[ $ratioCHECK == Y ]]; then
          #print DGE_info file
          printf "\n\nINPUT \t\t\t CHIP\n"
          cat ${TMPtmp}RATIO.final.txt | cut -f 1-2 | sed 's/:!!:/,/g'

          #ask user if it is fine, otherwise repeat
          while true; do
            read -r -p "\nIs the input to sample pairing for ratio building correct! [y/n]: " yn
            case $yn in
            [Y,y])
              ratioSWITCH="Y"
              mv ${TMPtmp}RATIO.final.txt ${TMP}
              break
              ;;
            [N,n])
              ratioCHECK="N"
              currNR=1
              rm -rf ${TMPtmp}RATIO.final.txt
              cp ${TMPtmp}RATIOlibs.txt ${TMPtmp}RATIOlibs.tmp
              break
              ;;
            *) echo "Please answer with Y or N!" ;;
            esac
          done

        fi

      done
    fi
  fi

fi

#--------------------------------------------------------f-------------------------------------------------
#if .bam libraries are present reqired variables need to be checked
if [[ $BAM == Y ]] || [[ $FORCE == Y ]] || [[ $demuxFASTA == Y ]]; then
  #ask for custom adaptor sequences if flag -A is set
  if [[ $customADAPTORS == Y ]]; then
    while true; do
      printf "\n\nDo you really want to use custom adaptor sequences istead of the standard Illumina sequences listed below?
        1st-mate adaptor  AGATCGGAAGAGCACACGTCT
        2nd-mate adaptor  AGATCGGAAGAGCGTCGTGTA\n\n"
      read -r -p "Please specify if you want to use custom adaptors [y or n]" yn
      case $yn in
      [Nn])
        fwADAPTOR="AGATCGGAAGAGCACACGTCT"
        rvADAPTOR="AGATCGGAAGAGCGTCGTGTA"
        break
        ;;
      [Yy])
        while true; do
          read -r -p "Please supply the 1st mate adaptor sequence [to keep standard type n]: " ads
          case $ads in
          [n])
            fwADAPTOR="AGATCGGAAGAGCACACGTCT"
            break
            ;;
          *) if [[ $ads =~ ^[ACGTacgt]+$ ]]; then
            fwADAPTOR=$ads
            break
          else
            echo "None DNA base in adaptor string!!!"
          fi ;;
          esac
        done

        while true; do
          read -r -p "Please supply the 2nd mate adaptor sequence (only required for paired-end) [to keep standard type n]: " ads
          case $ads in
          [n])
            rvADAPTOR="AGATCGGAAGAGCGTCGTGTA"
            break
            ;;
          *) if [[ $ads =~ ^[ACGTacgt]+$ ]]; then
            rvADAPTOR=$ads
            break
          else
            echo "None DNA base in adaptor string!!!"
          fi ;;
          esac
        done

        break
        ;;
      esac
    done
  else
    fwADAPTOR="AGATCGGAAGAGCACACGTCT"
    rvADAPTOR="AGATCGGAAGAGCGTCGTGTA"
  fi

  #ask if paired-fasta is really required
  if [[ $rawPAIRED == Y ]]; then
    while true; do
      printf "\nDo you really want to store a raw paired-end fasta file permanently?
This file is only useful for special paired-end analyses. For all other analyses
the fasta file automatically generated by the AnnotationPipeline can be used.\n\n"
      read -r -p "Generate raw paired-fasta file? [y or n]" yn
      case $yn in
      [Nn])
        echo "no raw-paired-fast file gets stored"
        rawPAIRED=N
        break
        ;;
      [Yy])
        printf "The raw-paired-fasta file gets processed
it will be stored together with the AnnotationPipeline fasta-file
in the individual-library folders.\n\n"
        rawPAIRED=Y
        break
        ;;
      *) echo "Please answer yes [y] or no [n]." ;;
      esac
    done
  fi

  #ask user how many Ns should be trimmed
  if [[ -z $N_TRIMM ]]; then
    while true; do
      read -r -p "Do you need to trimm random nucleotides from the read? [y or n]" yn
      case $yn in
      [Nn])
        N_TRIMM=0
        break
        ;;
      [Yy])
        while true; do
          read -r -p "How many randem nucletides have to be trimmed from each side? [0-8]" NTnr
          case $NTnr in
          [0-8])
            N_TRIMM=$NTnr
            break
            ;;
          *) echo "Please provide valid number." ;;
          esac
        done
        break
        ;;
      *) echo "Please answer yes [y] or no [n]." ;;
      esac
    done
  fi
fi

#---------------------------------------------------------------------------------------------------------
#set min max length variable according to default or user input
if [[ $TYPE == "sRNAseq" ]] || [[ $TYPE == "sRNAseqIP" ]]; then
  if [[ -z $MIN_LENGTH ]]; then MIN_LENGTH=18; fi
  if [[ -z $MAX_LENGTH ]]; then MAX_LENGTH=35; fi
elif [[ $TYPE == CapSeq ]]; then
  if [[ -z $MIN_LENGTH ]]; then MIN_LENGTH=35; fi
  if [[ -z $MAX_LENGTH ]]; then MAX_LENGTH=200; fi
else
  if [[ -z $MIN_LENGTH ]]; then MIN_LENGTH=18; fi
  if [[ -z $MAX_LENGTH ]]; then MAX_LENGTH=1000; fi
fi

#---------------------------------------------------------------------------------------------------------
#evaluate TRIMM input and reset  MAX_LENGTH variable to fit
if [[ $TYPE == RNAseq ]] || [[ $TYPE == GROseq ]] || [[ $TYPE == "RIPseq" ]]; then
  TRIMM="Y"
fi

if [[ $TRIMM == "Y" ]]; then
  if [[ $modTRIMM == N ]]; then
    if [[ $TYPE == RNAseq ]] || [[ $TYPE == GROseq ]] || [[ $TYPE == "RIPseq" ]]; then
      FIRST=6
      LAST=200
    else
      usage
      printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
      printf "\n  first/last nucleotide position for TRIMMing not set\n\n"
      printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n"
      exit 1
    fi
  else
    #set MAX_LENGTH to the LAST base in the TRIMMed output
    MAX_LENGTH=$LAST
  fi
fi
#---------------------------------------------------------------------------------------------------------
#test the MM variable for proper settings
if [[ -z ${MMinput} ]] && [[ $TYPE == CLIPseq ]] ; then
  MM="3"
elif [[ -z ${MMinput} && $TYPE == RNAseq ]]; then
  MM="2"
elif [[ -z ${MMinput} ]]; then
  MM="0"
fi

#override MM in case of SLAMseq data
if [[ $SLAM == Y ]]; then
  MM="3"
fi

#override MM in case of manually supplied MM settings
if [[  -n $MMinput ]]; then
  MM=$MMinput
fi

if [[ $MM -gt 3 ]] || [[ $MM -lt 0 ]]; then
  usage
  printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
  printf "\n  Missmatch setting -s is set wrong \n  only values between 0-3 are allowed - input was $MM\n\n"
  printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n"
  exit 1
fi

#---------------------------------------------------------------------------------------------------------
#test if extension variable only used for CHIPseq and if ok
if [[ $TYPE != CHIPseq ]] && [[ $EXTEND -gt 0 ]]; then
  usage
  printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
  printf "\n  extnding the reads is only allowed for CHIPseq\n\n"
  printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n"
  exit 1
elif [[ $EXTEND -lt 0 ]]; then
  usage
  printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
  printf "\n  read extension has to be >0\n\n"
  printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n"
  exit 1
fi

#---------------------------------------------------------------------------------------------------------
#test if only one wig otion is turned on
if [[ $WIG == Y ]] && [[ $WIG_FASTA == Y ]]; then
  usage
  printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
  printf "\n  both -W and -w are set.\n  Only one way of track normalization is allowed\n\n"
  printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n"
  exit 1
fi

#?@! if [[ $noNORM == Y ]]; then
#?@!   WIGcommand="no normalization - raw read-counts are reported"
#?@! else
if [[ $TYPE == sRNAseq ]] && [[ "$WIG_FASTA" != Y ]] && [[ "$WIG" != Y ]]; then
  WIGcommand="total-sRNA track normalization to 1M miRNAs; IP-sRNA track normalization to 10M reads in the input fasta-file"
  sleep 0.1s
elif [[ $WIG_FASTA == Y ]]; then
  #if option for FASTA normalization is set
  WIG_FASTA="Y"
  WIGcommand="track normalization to 10M reads in the input fasta file"
elif [[ $WIG == Y ]]; then
  #if option for uniq-mapper normalization is set
  WIG=Y
  WIGcommand="track normalization to 10M uniquely mapping reads"
elif [[ $TYPE == CHIPseq ]] || [[ $TYPE == DNAseq ]] || [[ $TYPE == sRNAseqIP ]]; then
  #if fasta-normalization is the default setting
  WIG_FASTA="Y"
  WIGcommand="track normalization to 10M reads in the input fasta file"
else
  #for all other TYPEs automatically normalize to mapped reads
  WIG=Y
  WIGcommand="track normalization to 10M uniquely mapping reads"
fi
#?@! fi
#---------------------------------------------------------------------------------------------------------
#create filtering classes for filtered size-profile depending on defined TYPE of LIBRARY
#alwayse sense and antisense annotations filtered

if [[ $TYPE == sRNAseq ]] || [[ $TYPE == sRNAseqIP ]]; then
  FILTERING_INPUT="rRNA,tRNA,mito,miRNA,pre_miRNA,snRNA,snoRNA,ncRNA"
elif [[ -z ${FILTERING_INPUT_RAW+x} ]]; then
  FILTERING_INPUT="rRNA,tRNA,mito"
fi

#convert commas in filtering variable to :
FILTERING_INPUT=$(echo $FILTERING_INPUT | tr ',' ':')

###################################################################################################

#! hack to stop runs where some libraries are lost

nFILESfinal=$(cat ${TMP}files.txt | wc -l )
if [[ ! $nFILES -eq $nFILESfinal ]]; then
  usage
  printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
  printf "\n  number of libraries changed during setup! initial: $nFILES  final: $nFILESfinal\n\n"
  printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n"
  exit 1
fi
###################################################################################################
#review setup part of the pipeline
RANDOMx=$(echo $RANDOM)
mkdir -p ${BASE_FOLDER_TMP}bam-files/

if [[ $TRIMM == Y ]]; then xTRIMM="Yes $FIRST - $LAST"; else xTRIMM="OFF"; fi
if [[ $BAM == Y ]] || [[ $FORCE == Y ]]; then
  if [[ $N_TRIMM -gt 0 ]]; then xN_TRIMM="$N_TRIMM nucleotides trimmed"; else xN_TRIMM="OFF"; fi
  if [[ $rawPAIRED == Y ]]; then xPAIRED="raw file will be generated"; else xPAIRED="OFF"; fi
  if [[ $customADAPTORS == Y ]]; then
    xADAPTOR="custom adaptor sequenses will be used for adaptor-clipping fw=${fwADAPTOR} rv=${rvADAPTOR}"
  else
    xADAPTOR="standard Illumina sequences will be used for adaptor-clipping fw=${fwADAPTOR} rv=${rvADAPTOR}"
  fi

  if [[ $BAM == Y ]]; then
    printf "file convertion from bam --> fasta
  adaptors for clipping= $xADAPTOR
  trimming of random nucleotides= $xN_TRIMM
  generation of paired-fasta file= $xPAIRED
  
" >${BASE_FOLDER_TMP}bam-files/bam_${RANDOMx}.txt

  elif [[ $FORCE == Y ]]; then
    printf "read preprocessing forced with flag -f
  adaptors for clipping= $xADAPTOR
  trimming of random nucleotides= $xN_TRIMM
  generation of paired-fasta file= $xPAIRED
  
" >${BASE_FOLDER_TMP}bam-files/bam_${RANDOMx}.txt

  fi

  #for CapSeq libraries warn user that 2nd mate will get deleted
  if [[ $TYPE == CapSeq ]]; then
    printf "  for paired-end data 2nd mate will be deleted due to TYPE-setting CapSeq
    
" >>${BASE_FOLDER_TMP}bam-files/bam_${RANDOMx}.txt
  fi

else
  touch ${BASE_FOLDER_TMP}bam-files/bam_${RANDOMx}.txt
fi

if [[ $Ychrom == Y ]]; then xYchrom="included in analysis"; else xYchrom="excluded from analysis"; fi

printf "\n\n############################################################################
annotation pipeline will be started with the following parameters:
general parameters
  library file= $FILE_CONTAINING_LIBRARIESraw - all $nFILES libraries are available
"
if [[ -f ${TMP}URLs.txt ]]; then
  printf "  libraries newly imported from NGS: \n"
  cat ${TMP}URLs.txt |
    awk '{ print "    lane_" NR":", $2 }' | sed 's/;/, /g'
  printf "\n"
fi

printf "  type of libraries= $TYPE 
  storage location= $FOLDER
  tmp storage location= $TMP

"

cat ${BASE_FOLDER_TMP}bam-files/bam_${RANDOMx}.txt
printf "processing parameters
  genome assembly version= $GENOME_VERSION
"

if [[ ${prepareREF} == yes ]]; then
  VERSIONext=" -- will be newly created on the fly"
else
  VERSIONext=""
fi

printf "  genome annotation version= $VERSION $VERSIONext
  Y-chromosome= $xYchrom
  length filtering= $MIN_LENGTH - $MAX_LENGTH
  trimming= $xTRIMM
  MM allowed in mapping= $MM
  track normalization= $WIGcommand
  read extension= $EXTEND
  filtered read classes= $FILTERING_INPUT
  color for track= $COLOR


"

while true; do
  read -r -p "Are the settings printed above correct? [y or n]" yn
  case $yn in
  [Yy]*) break ;;
  [Nn]*) exit ;;
  *) echo "Please answer yes or no." ;;
  esac
done

###################################################################################################
#in case of updating  - move old folder
if [[ $UPDATE == Y ]]; then
  rm -rf ${FOLDER}/individual-libraries/
  FOLDERbase=$(echo $FOLDER | tr -s / | sed 's/\/$//')
  mv $FOLDER ${FOLDERbase}_updated-to-$PIPELINEversion
fi

###################################################################################################

#start processing part of the pipeline

#create folders and move script-files to TMP directory
mkdir -p "$FOLDER"
mkdir -p "$TMP"
mkdir -p $UTILITY_LOCATION

rm -rf ${TMP}script-files/
cp -r $SCRIPT_DIR ${TMP}
SCRIPT_DIR=${TMP}script-files/
mv ${SCRIPT_DIR}head.sh ${SCRIPT_DIR}AP_${FOLDER_NAME}.sh
STARTscript=${SCRIPT_DIR}AP_${FOLDER_NAME}.sh
chmod 777 ${SCRIPT_DIR}*

#---------------------------------------------------------------------------------------------------------
#initiate log Folder/Files
LOGs=${FOLDER}LOGs/
rm -rf $LOGs
mkdir -p $LOGs
cd $LOGs

#---------------------------------------------------------------------------------------------------------
#print settings to log
printf "############################################################################
annotation pipeline will be started with the following parameters:
general parameters
  library file= $FILE_CONTAINING_LIBRARIESraw - all $nFILES libraries are available
  type of libraries= $TYPE 
  storage location= $FOLDER
  tmp storage location= $TMP

" >"${FOLDER}log.txt"
cat ${BASE_FOLDER_TMP}bam-files/bam_${RANDOMx}.txt >>"${FOLDER}log.txt"
printf "processing parameters
  genome assembly version= $GENOME_VERSION
  genome annotation version= $VERSION
  Y-chromosome= $xYchrom
  length filtering= $MIN_LENGTH - $MAX_LENGTH
  trimming= $xTRIMM
  MM allowed in mapping= $MM
  track normalization= $WIGcommand
  read extension= $EXTEND
  filtered read classes= $FILTERING_INPUT
  color for track= $COLOR " >>"${FOLDER}log.txt"

printf "\n\n############################################################################\n\n" >>"${FOLDER}log.txt"

#---------------------------------------------------------------------------------------------------------
#assemble full start command

if [[ ! -z "$USER" ]]; then Y="-U "; else Y=""; fi

if [[ ! -z "$FOLDER_NAME" ]]; then Z="-F "; else Z=""; fi

if [[ ! -z "$N_TRIMM" ]]; then N="-N "; else N=""; fi

printf "start command:
${SCRIPT_DIRraw}annotate_reads.sh  -i $FILE_CONTAINING_LIBRARIESraw $RERUN_OPTS ${Z}$FOLDER_NAME -t $TYPE -v $GENOME_VERSION -V $VERSION ${N}${N_TRIMM} -m $MIN_LENGTH -M $MAX_LENGTH -f $FIRST -l $LAST -s $MM -e ${EXTEND} -J ${maxCOUNT} ${Y}$USER" >>"${FOLDER}log.txt"
printf "\n\nVersion= ${PIPELINEversion}" >>"${FOLDER}log.txt"

commitID=$(git --git-dir ${SCRIPT_DIRraw}/.git merge-base HEAD origin/release)
printf "\n\nRelease Branch CommitID= ${commitID}" >>"${FOLDER}log.txt"
printf "\n\n############################################################################\n\n" >>"${FOLDER}log.txt"

printf "libraries analyzed in this run:\n" >>"${FOLDER}log.txt"
cat $FILE_CONTAINING_LIBRARIES >>"${FOLDER}log.txt"
printf "\n\nall $nFILES libraries are available" >>"${FOLDER}log.txt"
printf "\n\n############################################################################\n\n" >>"${FOLDER}log.txt"

printf "libraries as supplied by the user:\n" >>"${FOLDER}log.txt"
cat ${TMPtmp}files_orig.txt >>"${FOLDER}log.txt"
printf "\n\n############################################################################\n\n" >>"${FOLDER}log.txt"

printf "\n\n\nhelp file of used version:\n\n" >>"${FOLDER}log.txt"
usage >>"${FOLDER}log.txt"
printf "\n\n############################################################################\n\n" >>"${FOLDER}log.txt"

#delete temporary folder used for pre-processing
#?@ if [[  $DEBUG != Y ]]; then
#?@   rm -rf ${TMPtmp}
#?@ fi

###################################################################################################
###################################################################################################
#start analysis
#---------------------------------------------------------------------------------------------------------

COMMAND="$STARTscript"
VARI="USER_NAME=${USER_NAME},TYPE=${TYPE},SLAM=${SLAM},BASE_FOLDER=${BASE_FOLDER},FOLDER=${FOLDER},FOLDER_NAME=${FOLDER_NAME},RUNname=${RUNname},TMPdir=${TMP},LIB_STORAGE_FOLDER=${LIB_STORAGE_FOLDER},SINGULARITYdir=${SINGULARITYdir},downDIR=${downDIR},DEBUG=${DEBUG},VERSION=${VERSION},asmHUBpath=${asmHUBpath},ASMdir=${ASMdir},ASMname=${ASMname},DEMUXonly=${DEMUXonly},BAM=${BAM},fwADAPTOR=${fwADAPTOR},rvADAPTOR=${rvADAPTOR},N_TRIMM=${N_TRIMM},rawPAIRED=${rawPAIRED},onlyPAIRED=${onlyPAIRED},FASTQout=${FASTQout},FASTQoutRAW=${FASTQoutRAW},SUBSAMPLE=${SUBSAMPLE},MIN_LENGTH=${MIN_LENGTH},MAX_LENGTH=${MAX_LENGTH},RAW=${RAW},TRIMM=${TRIMM},FIRST=${FIRST},LAST=${LAST},INVERT=${INVERT},Ychrom=${Ychrom},RANDOMmulti=${RANDOMmulti},MM=${MM},FILTERING_INPUT=${FILTERING_INPUT},WIG=${WIG},WIG_FASTA=${WIG_FASTA},spikeINnorm=${spikeINnorm},noNORM=${noNORM},EXTEND=${EXTEND},COMPUTING=${COMPUTING},GRIDsystem=${GRIDsystem},keepTMP=${keepTMP},SCRIPT_DIR=${SCRIPT_DIR},BASE_UTILITY_LOCATION=${BASE_UTILITY_LOCATION},UTILITY_LOCATION=${UTILITY_LOCATION},RELEASE=${RELEASE},VERSION=${VERSION},UTILITY_DIR=${UTILITY_DIR},nFILES=${nFILES},FILE_CONTAINING_LIBRARIES=${FILE_CONTAINING_LIBRARIES},FOLDER_NAME=${FOLDER_NAME},GENOME_VERSION=${GENOME_VERSION},subCOLOR=${subCOLOR},FORCE=${FORCE},SYSTEM=${SYSTEM},LOGs=${LOGs},prepareREF=${prepareREF},nSPLITS=${nSPLITS},SE=${SE},SE2nd=${SE2nd},maxCOUNT=${maxCOUNT},demuxFASTA=${demuxFASTA},autoViewLimits=${autoViewLimits},only5end=${only5end},PingPong=${PingPong},DGE=${DGE},GEO=${GEO},noSTRANDED=${noSTRANDED},FORCEimport=${FORCEimport},exportBAM=${exportBAM},exportBAMuncollapsed=${exportBAMuncollapsed},exportSalmon=${exportSalmon},RATIOtracks=${RATIOtracks},GENOMEdir=${GENOMEdir},newCOLLECTION=${newCOLLECTION},prepANNOTATIONgff=${prepANNOTATIONgff},prepGENOMEfasta=${prepGENOMEfasta},prepTRANSCRIPTOMEfasta=${prepTRANSCRIPTOMEfasta},prepCDSfasta=${prepCDSfasta},prepNCRNAfasta=${prepNCRNAfasta},extraSEQ=${extraSEQ},FORCEquant_unstranded=${FORCEquant_unstranded}"

#---------------------------------------------------------------------------------------------------------
if [[ $COMPUTING == C ]]; then
  if [[ $GRIDsystem == SLURM ]]; then
    sbatch $COMMAND ${VARI}
  elif [[ $GRIDsystem == GRIDENGINE ]]; then
    qsub $COMMAND ${VARI}
  else
    echo wrong submission system for clusters
    exit
  fi
else
  if [[ -z ${SLURM_JOB_ID+x} ]] && [[ $SYSTEM == CLIP ]]; then
    srun --cpus-per-task=10 --mem-per-cpu=5g --qos=short $COMMAND ${VARI}
  else
    $COMMAND ${VARI}
  fi
  wait
fi
exit
###################################################################################################

printf "
pipeline is running safely on the cluster

"

if [[ $GENOME_VERSION == ASM ]]; then
  printf "
to use assembly $VERSION please load the following hub into UCSC
${HTTP}ASSEMBLIES/${VERSION}/hub.txt"
fi

printf "
your resuls will be deposited into:
$HTTP

lykke til og en fin dag!
\n\n"

exit
