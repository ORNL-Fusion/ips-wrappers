#!/bin/bash
setfacl -R -m g:atom:rwX `pwd`
find `pwd` -type d | xargs setfacl -R -m d:g:atom:rwX
#chmod -R ug+rwX `pwd`
#chgrp -R atom `pwd`
#find `pwd` -type d -exec chmod g+s '{}' \;
