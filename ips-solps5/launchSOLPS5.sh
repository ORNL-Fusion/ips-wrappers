#!/bin/csh
#
# Takes 2 arguments ... baserun folder location & run folder location
# e.g., 
# ./launchSOLPS5.sh $(readlink -f baserun) $(readlink -f thread_00)

echo There are $#argv arguments
echo These are $argv[*]

set baserunlocation="$argv[1]"
set runlocation="$argv[2]"

cd /project/projectdirs/atom/atom-install-edison/solps-5

# This is required to avoid passing the input args to the SOLPS scripts
set argv=

source setup.csh
sbr
echo $PWD

#cd d3d_test/Donly_profdb

echo "Generating random IPS directory name"
echo "------------------------------------"

set rand=`cat /dev/urandom | tr -cd 'a-f0-9' | head -c 8`
set ipsdir=`echo ips-$rand`
mkdir $ipsdir
cd $ipsdir
echo $PWD

echo "Copying baserun to SOLPS location"
echo "---------------------------------"
cp -r $baserunlocation .
cp -r $runlocation run
cd run

echo "Executing ./cleanRun.sh"
echo "-----------------------"

./cleanRun.sh

echo "Executing b2run b2mn < input.dat > run.log"
echo "------------------------------------------"

b2run b2mn < input.dat > run.log # run the code

echo "Executing eirobjx < input.dat > output1"
echo "---------------------------------------"

eirobjx < input.dat > output1 # not sure what this does?

echo "Executing ./lookAtOutput.sh"
echo "---------------------------"

./lookAtOutput.sh

echo "Executing 'idl -e 'solps_plots'"
echo "-------------------------------"

idl -e 'solps_plots' # creates the ips-solps.nc file

echo "Copying results to ips run location"
echo "-----------------------------------"

set solpsrunlocation="$PWD"
echo 'cp -r $solpsrunlocation/* $runlocation/'
cp -r $solpsrunlocation/* $runlocation/

