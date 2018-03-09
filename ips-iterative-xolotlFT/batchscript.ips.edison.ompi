#!/bin/bash -l
#SBATCH --account=atom
#SBATCH --partition regular
#SBATCH --nodes=2
#SBATCH --time=12:00:00
#SBATCH --output=log.slurm.stdOut

cd $SLURM_SUBMIT_DIR   # optional, since this is the default behavior

# Production
#source /project/projectdirs/atom/atom-install-edison/ips-wrappers/env.ips.edison 
# Devel
#source /project/projectdirs/atom/atom-install-edison/ips-wrappers-devel/env.ips.edison 
# Greendl1
#source /project/projectdirs/atom/users/greendl1/code/ips-wrappers/env.ips.edison 
# tyounkin
source /project/projectdirs/atom/users/$USER/ips-wrappers/env.ips.edison

module swap PrgEnv-intel PrgEnv-gnu

#without task pool:
#In megabucky line: ips.py --config=ips.megabucky.config --platform=$IPS_PLATFORM_FILE --log=log.framework 2>>log.stdErr 1>>log.stdOut
#$IPS_PATH/bin/ips.py --config=ips.config --platform=conf.ips.edison --log=log.framework 2>>log.stdErr 1>>log.stdOut

#to run multiple ftridyn runs in parallel as task pool and Xolotol in multiple nodes
export OMPI_ROOT=/project/projectdirs/atom/users/elwasif/ompi/install_4.0
export PATH=$OMPI_ROOT/bin:$PATH
export LD_LIBRARY_PATH=$OMPI_ROOT/lib64:$OMPI_ROOT/lib:$LD_LIBRARY_PATH
export MANPATH=$OMPI_ROOT/share/man:$MANPATH

srun -n $SLURM_NNODES --ntasks-per-node=1 hostname > .node_names.$$
for n in $(cat .node_names.$$); do echo "$n slots=24" >> .hostfile.$$ ; done
orte-dvm --hostfile ./.hostfile.$$ &

##orte-dvm --daemonize --report-uri ./dvm_uri
#unset ORTE_HNP_DVM_URI
#orte-dvm --debug-daemons --report-uri ./dvm_uri.$$ &
sleep 5
#export ORTE_HNP_DVM_URI=$(cat ./dvm_uri.$$)
##export LD_PRELOAD=$OMPI_ROOT/lib64/libmpi.so:$OMPI_ROOT/lib64/libmpi_mpifh.so
##export LD_PRELOAD=$OMPI_ROOT/lib64/libmpi_mpifh.so

$IPS_PATH/bin/ips.py --config=ftx.config --platform=conf.ips.edison.ompi --log=log.framework 2>>log.stdErr 1>>log.stdOut

egrep -i 'error' log.* > log.errors
./setPermissions.sh

