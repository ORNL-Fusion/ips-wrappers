module swap PrgEnv-${PE_ENV,,} PrgEnv-gnu
export MPICC="cc -shared"
conda env create -f FTX_env.yml --force
