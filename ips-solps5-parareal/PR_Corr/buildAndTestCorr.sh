#!/bin/bash

CORRECT_NEUTRALS=1

SRC=src_code_PR_corr

rm *.o *.mod
ftn -c $SRC/*.F

if [ "$CORRECT_NEUTRALS" -eq "1" ]; then
echo "Building with CORRECT_NEUTRALS=1"
ftn $SRC/SOLPS-PR_corr.f90 *.o -o corr_solps
else
echo "Building with CORRECT_NEUTRALS=0"
ftn $SRC/SOLPS-PR_corr_No-neutral-corr.f90 *.o -o corr_solps
fi

rm *.o *.mod
cp corr_solps $SRC/
cd $SRC

TEST_OUTPUT=test_output
rm $TEST_OUTPUT 
./corr_solps fine_out_values.0003.0027 coarse_out_values.0004.0027 coarse_out_values.0003.0027 $TEST_OUTPUT 
rm corr_solps
cd ..
