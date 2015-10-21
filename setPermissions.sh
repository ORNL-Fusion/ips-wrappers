#!/bin/bash
setfacl -R -m g:atom:rwX `pwd`
find `pwd` -type d | xargs setfacl -R -m d:g:atom:rwX
