#!/usr/bin/env bash
set -u

#determine location of the script-directory
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done

SCRIPT_DIRraw="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
SCRIPT_DIR="${SCRIPT_DIRraw}/"

cd $SCRIPT_DIR


testSettings=$(diff ${SCRIPT_DIR}settings.txt ${SCRIPT_DIR}.settings.txt)
#ask user if settings.txt has been changed on purpose
if [[ ! -z $testSettings ]]; then
while true; do
  read -r -p "
  Did you modify settings.txt on purpose ? [y or n]" yn
  case $yn in
  [Yy])
    cp ${SCRIPT_DIR}settings.txt ${SCRIPT_DIR}.settings.txt
    break
    ;;
  [Nn])
    printf "\n\nUsing old settings\n\n"
    break
    ;;
  *) echo "Please answer yes [y] or no [n]." ;;
  esac
done
fi

git reset --hard
git pull --rebase origin release-extern
if [[ $? -gt 0 ]]; then
  echo $? exit due to git login error
  exit 
fi

#replace settings.txt file
rm -rf ${SCRIPT_DIR}settings.txt
cp ${SCRIPT_DIR}.settings.txt ${SCRIPT_DIR}settings.txt

#read in settings file
eval "$(<${SCRIPT_DIR}settings.txt)"
#cd $SCRIPT_DIR
SINGULARITY_DIR="${BASE_UTILITY_LOCATION}/simg/"

#replace place-holders with final settings
 sed -i `` "s|XX-dm3-VERSION-XX|${default_VERSION_dm3}|;
            s|XX-dm6-VERSION-XX|${default_VERSION_dm6}|;
            s|XX-BASE_UTILITY_LOCATION-XX|${BASE_UTILITY_LOCATION}|;
            s|XX-BASE_LIBRARY_STORAGE_LOCATION-XX|${BASE_LIBRARY_STORAGE_LOCATION}|;
            s|XX-BASE_SINGULARITY_LOCATION-XX|${SINGULARITY_DIR}|;
            s|XX-BASE_RESULTS_FOLDER-XX|${BASE_FOLDER}|;
            s|XX-BASE_FOLDER_TMP-XX|${BASE_FOLDER_TMP}|;
            s|XX-HTTP_LINK-XX|${HTTP_PATH}|" `` ${SCRIPT_DIR}annotate_reads.sh
chmod 777 ${SCRIPT_DIR}annotate_reads.sh 
chmod 777 ${SCRIPT_DIR}script-files/*

#replace SBATCH options
SBATCHraw=$(awk '/^#@!@#/,0' ${SCRIPT_DIR}settings.txt | tail -n +2 )
SBATCHrepl="$(echo "${SBATCHraw}" | sed ':a;N;$!ba;s/\n/\\n/g' | sed 's/\$/\\$/g')"
if [[ -z "${SBATCHrepl}" ]]; then
  #@ echo empty SBATCHrepl
  SBATCHrepl="#"
#@ else 
#@   echo SBATCHrepl=$SBATCHrepl
fi
HYPER=$(grep HYPER ${SCRIPT_DIR}settings.txt  | tr '=' '\t' | cut -f 2)
if [[ $HYPER == Y ]]; then HYPERval=2; else HYPERval=1; fi


for FILE in ${SCRIPT_DIR}script-files/* ; do
  if [[ ! $FILE == ${SCRIPT_DIR}script-files/update.sh ]]; then
    sed -i `` "s|@HYPER@|${HYPERval}|;s|@SBATCH@|${SBATCHrepl}|" $FILE
  fi
done

#execute potential extra commands required
chmod 777 ${SCRIPT_DIR}script-files/.update.extra.sh
${SCRIPT_DIR}script-files/.update.extra.sh

#clean up directory again
rm -rf ${SCRIPT_DIR}setup.sh
mv ${SCRIPT_DIR}script-files/update.sh ${SCRIPT_DIR}update.sh
chmod 777 ${SCRIPT_DIR}update.sh

#update singularity packages
cd ${SINGULARITY_DIR}
#singularity pull --disable-cache --force --name APmaster.simg shub://dominik-handler/AP_singu2:ap_master
#singularity pull --disable-cache --force --name R.simg shub://dominik-handler/AP_singu2:r@55415f9eb398dc47af6987eab2aa1e89f1e4aee2
#!@ rm -rf ${SINGULARITY_DIR}/*.simg

#!@ wget https://brenneckelab.imba.oeaw.ac.at/tmp/SINGU/APmaster.simg
#!@ wget https://brenneckelab.imba.oeaw.ac.at/tmp/SINGU/R.simg
