#!/bin/bash
chmod -R ug+rwX `pwd`
chgrp -R atom `pwd`
find `pwd` -type d -exec chmod g+s '{}' \;
