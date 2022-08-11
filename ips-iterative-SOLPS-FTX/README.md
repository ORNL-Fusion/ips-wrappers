# Installing miniconda
If conda is not already installed on your system you can install miniconda using the commands
```
curl -O https://repo.anaconda.com/miniconda/Miniconda3-VERSION.sh
bash Miniconda3-VERSION.sh
```
where VERSION is a replaced with a specific version available [from conda](https://docs.conda.io/en/latest/miniconda.html)

# Installing the conda environment.
To install the conda environment, run the command

```
conda env create -f SOLPS-FTX_env.yml --force
```
