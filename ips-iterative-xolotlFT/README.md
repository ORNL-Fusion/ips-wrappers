Hello_world is (almost) the simplest possible IPS run, exercising two components 
DRIVER  hello_driver.py and WORKER  hello_worker.py.  
"Almost" because it would be possible to have an IPS simulation with only a driver.
Hello world does not use the Plasma State or use any physics components.

There is one external connection, to an environment file that is sourced
in the slurm batch script, run_slurm.  The environment file env.ips.edison is maintained at

/project/projectdirs/atom/atom-install-edison/ips-wrappers/env.ips.edison
  
The components reside in the AToM project IPS_CSWIM_WRAPPERS directory as defined in env.ips.edison,
but they have been copied here into the source directory for ease of viewing by the user.

To run the "simulation" submit the batch script from the command line

sbatch batchscript.ips.edison

To clean all the run files and start with just the input deck run 

./cleanIpsRun.sh


