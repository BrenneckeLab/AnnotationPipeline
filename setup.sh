#!/usr/bin/env bash

#determine location of the script-directory
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done

SCRIPT_DIRraw="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
SCRIPT_DIR="${SCRIPT_DIRraw}/"

#ask user if settings.txt has been filled
while true; do
  read -r -p "
  Did you modify settings.txt with the correct paths for your environment? [y or n]" yn
  case $yn in
    [Yy])
      break
    ;;
    [Nn])
      printf "\n\nPlease modify settings.txt and rerun setup.sh\n\n"
      exit
    ;;
    *) echo "Please answer yes [y] or no [n]." ;;
  esac
done

#read in settings file
eval "$(<${SCRIPT_DIR}settings.txt)"

#cd $SCRIPT_DIR
SINGULARITY_DIR="${BASE_UTILITY_LOCATION}/simg/"

#replace place-holders with final settings
sed -i `` "s|XX-dm3-VERSION-XX|${default_VERSION_dm3}|;
            s|XX-dm6-VERSION-XX|${default_VERSION_dm6}|;
            s|XX-Cel-VERSION-XX|${default_VERSION_Cel}|;
            s|XX-BASE_UTILITY_LOCATION-XX|${BASE_UTILITY_LOCATION}|;
            s|XX-BASE_LIBRARY_STORAGE_LOCATION-XX|${BASE_LIBRARY_STORAGE_LOCATION}|;
            s|XX-BASE_SINGULARITY_LOCATION-XX|${SINGULARITY_DIR}|;
            s|XX-BASE_RESULTS_FOLDER-XX|${BASE_FOLDER}|;
            s|XX-BASE_FOLDER_TMP-XX|${BASE_FOLDER_TMP}|;
            s|XX-SCRIPT_DIR-XX|${SCRIPT_DIR}|;
s|XX-HTTP_LINK-XX|${HTTP_PATH}|" ${SCRIPT_DIR}annotate_reads.sh
chmod 777 ${SCRIPT_DIR}annotate_reads.sh
chmod 777 ${SCRIPT_DIR}script-files/*

#replace SBATCH options
SBATCHraw=$(awk '/^#@!@#/,0' ${SCRIPT_DIR}settings.txt | tail -n +2 )
SBATCHrepl="$(echo "${SBATCHraw}" | sed ':a;N;$!ba;s/\n/\\n/g' | sed 's/\$/\\$/g')"
HYPER=$(grep HYPER ${SCRIPT_DIR}settings.txt  | tr '=' '\t' | cut -f 2)
if [[ $HYPER == Y ]]; then HYPERval=2; else HYPERval=1; fi

for FILE in ${SCRIPT_DIR}script-files/* ; do
  if [[ ! $FILE == ${SCRIPT_DIR}script-files/update.sh ]]; then
    sed -i "s|@SBATCH@|${SBATCHrepl}|;
          s|@HYPER@|$HYPERval|" $FILE
  fi
done

#create utility location and copy required files there
if [[ $RUN_INSTALLATION == Y ]]
then
  rm -rf ${BASE_UTILITY_LOCATION}
  mkdir -p ${BASE_UTILITY_LOCATION}
  touch ${BASE_UTILITY_LOCATION}/versions.txt
  chmod 777 ${BASE_UTILITY_LOCATION}
  cp -r ${SCRIPT_DIR}utility-files/* ${BASE_UTILITY_LOCATION}/
  mkdir -p ${BASE_LIBRARY_STORAGE_LOCATION}
  chmod 777 ${BASE_LIBRARY_STORAGE_LOCATION}
  
  #download required singularity images
  rm -rf ${SINGULARITY_DIR}
  mkdir -p ${SINGULARITY_DIR}
  cd ${SINGULARITY_DIR}
  echo $SINGULARITY_DIR
  #@ singularity pull --disable-cache --name  APmaster.simg shub://dominik-handler/AP_singu2:ap_master
  #@ singularity pull --disable-cache --name  R.simg shub://dominik-handler/AP_singu2:r@55415f9eb398dc47af6987eab2aa1e89f1e4aee2
  wget https://brenneckelab.imba.oeaw.ac.at/tmp/SINGU/APmaster.simg
  wget https://brenneckelab.imba.oeaw.ac.at/tmp/SINGU/R.simg  
fi

cd $SCRIPT_DIR

rm -rf ${SCRIPT_DIR}setup.sh
mv ${SCRIPT_DIR}script-files/update.sh ${SCRIPT_DIR}update.sh
chmod 777 ${SCRIPT_DIR}update.sh
cp ${SCRIPT_DIR}settings.txt ${SCRIPT_DIR}.settings.txt

printf "\n\ninstallation successful\n\n"
