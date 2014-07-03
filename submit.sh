echo "[ making the project archive... ]"
tar chzvf input.tar.gz ScanChain.* Include.C  goodrun.* doTest.C libMiniFWLite.so source

echo "[ submitting to condor... ]"
condor_submit test_condor
