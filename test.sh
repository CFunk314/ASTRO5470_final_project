#!/bin/bash

start=$(date +%s)
workingdir=$PWD
cd ${workingdir}

make clean
make '-f' Makefile
printf "\nMake complete. Running test 1...\n"

cp input_test1.nml input.nml

./a.out

python tests.py '-t' '1'

printf "\nTest 1 complete. Running test 2...\n"

cp input_test2.nml input.nml

./a.out

python tests.py '-t' '2'

printf "\nTest 2 complete. Running test 3...\n"

cp input_test3.nml input.nml

./a.out

python tests.py '-t' '3'

printf "\nTest 3 complete. Cleaning up...\n"

make clean
rm *.data

end=$(date +%s)
runtime=$((end-start))
printf "\nTests completed in $runtime seconds.\n"


