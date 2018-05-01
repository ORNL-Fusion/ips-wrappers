# How to develop on your own copy of ips-wrappers
## Install the IPS
Skip this if you've already installed the IPS. 

1. Create an IPS directory and clone the IPS-framework, wrappers, and examples repos.
```
mkdir IPS
cd IPS
git clone https://github.com/HPC-SimTools/IPS-framework.git ips-framework
git clone https://github.com/ORNL-Fusion/ips-wrappers.git
git clone https://github.com/ORNL-Fusion/ips-examples.git
```
2. Export the `IPS_DIR` environment variable
```
export IPS_DIR=${PWD}
```
3. Add this to your `.bashrc` or otherwise so it's there next time you open a shell (Note: Adapt for `csh` or otherwise).
```
echo 'export IPS_DIR='${PWD} >> ~/.bashrc 
```

## Run the example

1. Source the IPS environemnt
```
cd $IPS_DIR
source ips-wrappers/env.ips
```
2. Run the ABC example
  * Locally
```
cd ips-examples/ABC_example
ips.py --simulation=ABC_simulation.config --platform=platform.conf
```
  * On a batch system (e.g., Edison at NERSC)
```
cd ips-examples/ABC_example
sbatch Edison_run
```
To clean all the run files and start with just the input deck run 
```
./cleanIpsRun.sh
```
