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

then in your [batchscript.ips.edison](https://github.com/ORNL-Fusion/ips-atom/blob/dlg-devel/template.batchscript.ips.edison), comment out the default, and add the following

```
# Production
# source /project/projectdirs/atom/atom-install-edison/ips-wrappers/env.ips.edison
# Me
source /project/projectdirs/atom/users/$USER/code/ips-wrappers/env.ips.edison
```

## After making some changes you want to make available
Check was git has to say about the things you changed
```
git status
```
Add the files you want to commit to the staging area
```
git add filename1 filename2
```
Commit the files (this is only local - it does not push to a remote like SVN)
```
git commit -m 'My informative commit message describing a new awesome feature'
```
Pull any remote changes prior to pushing
```
git pull
```
Now push your changes to the default remote (probably github unless you already knew what you are doing)
```
git push
```


