#!/bin/bash

SRC=src_code_PR_corr

rm *.o *.mod
ftn -c $SRC/*.F
ftn $SRC/SOLPS-PR_corr.f90 *.o -o corr_solps
rm *.o *.mod
cp corr_solps $SRC/
cd $SRC

TEST_OUTPUT=test_output
rm $TEST_OUTPUT 
./corr_solps fine_out_values.0003.0027 coarse_out_values.0004.0027 coarse_out_values.0003.0027 $TEST_OUTPUT 
rm corr_solps
cd ..
