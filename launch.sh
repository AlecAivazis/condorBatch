#!/bin/bash

# build the environment 
export CMS_PATH=/code/osgcode/cmssoft/cms
export SCRAM_ARCH=slc5_amd64_gcc472
source /code/osgcode/cmssoft/cmsset_default.sh
source /code/osgcode/fgolf/5.30-patches/bin/thisroot.sh
export LD_LIBRARY_PATH=${ROOTSYS}/lib:${LD_LIBRARY_PATH}:${_CONDOR_SCRATCH_DIR}
export PATH=${HOME}/bin:${ROOTSYS}/bin:${PATH}:${_CONDOR_SCRATCH_DIR}
export PYTHONPATH=${ROOTSYS}/lib:${PYTHONPATH}

# set the environment variables to extract the file
export COPYDIR=/hadoop/cms/store/user/aaivazis
# the name given to the file by root
export TEMP_FILE_NAME=baby.root
# the name of the file on the target directory
export OUTPUT=ntuple$1.root

# change the working directory to the condor tmp
cd ${_CONDOR_SCRATCH_DIR}

# unarchive the included archives
for f in *.tar.gz
do
  tar -xzvf ${f}
done

# check that everything made it okay
echo "filesystem: "
pwd
ls

echo "running baby maker"
root -b -q doTest.C

# copy the output directory to hadoop
echo "copying files to ${COPYDIR}"
lcg-cp -b -D srmv2 --vo cms --connect-timeout 2400 --verbose file://`pwd`/${TEMP_FILE_NAME} srm://bsrm-3.t2.ucsd.edu:8443/srm/v2/server?SFN=${COPYDIR}/${OUTPUT}
