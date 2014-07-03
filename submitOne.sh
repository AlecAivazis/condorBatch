# the directory to loop over
dir=/hadoop/cms/store/user/aaivazis/samples
# keep track of the batch counter
batch=0
for f in ${dir}/*.root
do
  echo "[ building files directory... ]"
  rm -rf source
  mkdir source

  echo "[ pulling appropriate root files... ]"
  cp ${f} source

  echo "[ making the project archive... ]"
  tar czvf input-${batch}.tar.gz ScanChain.* Include.C  goodrun.* doTest.C libMiniFWLite.so source

  echo "[ creating condor submission file... ]"
  echo "batch = ${batch}"
  rm -rf batch_condor
  cat <<EOF > batch_condor
universe = vanilla
Executable = test.sh
Requirements = Memory >= 199 && OpSys == "LINUX"&& (Arch != "DUMMY" )&& Disk > 1000000
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = input-${batch}.tar.gz, CORE.tar.gz 
Output = logs/out
Error = logs/error
Log = logs/log
+DESIRED_Sites="UCSD"
notify_user = alec@aivazis.com
x509userproxy = /tmp/x509up_u31136
Queue 1
Arguments = ${batch}
EOF

  echo "[ submitting to condor... ]"
  condor_submit batch_condor

  echo "[ cleaning up... ]"
  rm -rf source

  # increment the batch counter
  batch=$((${batch}+1))
done
