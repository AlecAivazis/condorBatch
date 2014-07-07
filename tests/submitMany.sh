# the directory to loop over
dir=/hadoop/cms/store/user/aaivazis/samples
# the number of files per batch
count=2
# grab the batch number of the job from the command line argument
batch=$1

echo "[ building files directory... ]"
rm -rf files
mkdir files

echo "[ pulling appropriate root files... ]"
cp $(ls -1 ${dir}/*.root | tail -n +$((${batch}*${count}+1)) | head -n ${count}) files

echo "[ making the project archive... ]"
tar czvf input-${batch}.tar.gz ScanChain.* Include.C  goodrun.* doTest.C libMiniFWLite.so ntuple1.root CORE.tar.gz

echo "[ creating condor submission file... ]"
rm -rf test_condor
cat <<EOF > submit_condor
universe = vanilla
Executable = test.sh
Requirements = Memory >= 199 && OpSys == "LINUX"&& (Arch != "DUMMY" )&& Disk > 1000000
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = input-${batch}.tar.gz
Output = logs/out
Error = logs/error
Log = logs/log
+DESIRED_Sites="UCSD"
notify_user = alec@aivazis.com
x509userproxy = /tmp/x509up_u31136
Queue 1
EOF

echo "[ submitting to condor... ]"
condor_submit submit_condor

echo "[ cleaning up... ]"
rm -rf files
