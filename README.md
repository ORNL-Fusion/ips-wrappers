# How to develop on your own copy of ips-wrappers

**(Assuming you are on Edison and are an ATOM NERSC team member - if not, just ask)**

```
cd /project/projectdirs/atom/users
mkdir $USER
cd $USER
module load git
git clone https://github.com/ORNL-Fusion/ips-atom.git ips-wrappers
cd ips-wrappers
git checkout dlg-devel
git submodule update --init --recursive
```

You also will want to create your own branch, i.e., 

```
git checkout -b my-branch-name
```

then in your [batchscript.ips.edison](ips-atom/template.batchscript.ips.edison), comment out the default, and add the following

```
# Production
# source /project/projectdirs/atom/atom-install-edison/ips-wrappers/env.ips.edison
# Me
source /project/projectdirs/atom/users/$USER/code/ips-wrappers/env.ips.edison
```


