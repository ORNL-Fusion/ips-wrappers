#! /bin/tcsh -f
##$ -cwd
##$ -N Task_share 
##$ -m e
##$ -l h_rt=0:30:00
##$ -pe impi_hydra 1
##$ -S /bin/tcsh

# Script to copy and eventually remove input files for PR_correction code. Need this script as fortran 77 is case insensitive.
#module swap pgi intel

#executable forPR correction
#Fine_file 
#New G _file 
#OLD G _file 
#OUTPUT_file
#SIMPATH is self.suffix in python file
# output_file2 & fine_file2 correspond to the b2.sputter.parameters files, so only the fine is kept
setenv fine_file $1
setenv newG_file $2
setenv oldG_file $3
setenv fine_file2 $4
setenv fine_wall $7
setenv output_file1 $8
setenv output_file2 $9
setenv output_wall $10
setenv SIMPATH $11

echo "$SIMPATH"

cp $fine_file fine$SIMPATH
cp $newG_file newg$SIMPATH
cp $oldG_file oldg$SIMPATH

/pfs/work/dsam/SOLPS_bin/PR_Corr/corr_solps fine$SIMPATH newg$SIMPATH oldg$SIMPATH out$SIMPATH 

mv out$SIMPATH $output_file1
cp $fine_file2 $output_file2
cp $fine_wall $output_wall
rm fine$SIMPATH newg$SIMPATH oldg$SIMPATH

#rm Fine$SIMPATH NewG$SIMPATH OldG$SIMPATH out$SIMPATH

#set SIMPATH=$1





