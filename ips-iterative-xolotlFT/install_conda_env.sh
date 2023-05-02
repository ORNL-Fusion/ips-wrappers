module swap PrgEnv-${PE_ENV,,} PrgEnv-gnu
export MPICC="cc -shared"
CONDA_PKGS_DIRS=$(mktemp -d) conda env create -f FTX_env.yml --force -p /global/homes/a/a7l/.conda/envs/ftx
