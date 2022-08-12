#!/bin/sh
echo "Running shell script for Fractal Tridyn Input file FTridyn.IN"
echo "command line is: " $FTRIDYN_PATH"/bin/FTridyn_Clean < FTridyn.IN 1 1 1"
$FTRIDYN_PATH/bin/FTridyn_Clean < FTridyn.IN 1 1 1 
exit 0
