# This bash script loops over every file in a specific directory and 
# submits individual HTCondor jobs to process the file. 
# author: alec aivazis

# the directory to loop over
#dir=/hadoop/cms/store/user/aaivazis/samples/signal/raw/1000
dir=/hadoop/cms/store/user/aaivazis/samples
# keep track of the batch counter
batch=0
# store the directory that we start in so we can move files back
loc=$(pwd)

echo ${loc}

echo "[ making the project archive... ]"
cd include
tar czvf include.tar.gz *
mv include.tar.gz ..
cd ..

# loop over every file in the target directory
for f in ${dir}/*.root
do

  echo "[ entering target directory... ]"
  cd ${dir}
  echo $(pwd)
  echo "[ compressing target root file... ]"
  echo $(basename ${f})
  tar -zcf source-${batch}.tar.gz $(basename ${f})
  echo "[ retrieving compressed file... ]"
  mv source-${batch}.tar.gz ${loc}

  cd ${loc}

  echo "[ creating condor submission file... ]"
  echo "batch = ${batch}"
  rm -rf batch_condor
  cat <<EOF > batch_condor
universe = vanilla
Executable = exec.sh
Requirements = Memory >= 199 && OpSys == "LINUX"&& (Arch != "DUMMY" )&& Disk > 1000000
Should_Transfer_Files = IF_NEEDED
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = include.tar.gz, source-${batch}.tar.gz
Output = logs/out_${batch}
Error = logs/error_${batch}
Log = logs/log_${batch}
+DESIRED_Sites="UCSD"
notify_user = alec@aivazis.com
x509userproxy = /tmp/x509up_u31136
Arguments = ${batch}
Queue
EOF

  echo "[ submitting to condor... ]"
  #condor_submit batch_condor

  echo "[ cleaning up... ]"
  rm -rf batch_condor

  # increment the batch counter
  batch=$((${batch}+1))
done
