#!/bin/bash
setfacl -R -m g:atom:rwx `pwd`
setfacl -R -m o::rx `pwd`
find `pwd` -type d | xargs setfacl -R -m d:g:atom:rwx
find `pwd` -type d | xargs setfacl -R -m d:o::rx
