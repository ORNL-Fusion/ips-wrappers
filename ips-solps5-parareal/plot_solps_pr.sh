#! /bin/tcsh -f
###Make sure you have b2time.0001.0001.nc in the current directory.
###Then cp b2time.0001.0001.nc b2time.nc in the current directory.
###Make sure all b2time*.*.nc files are one directory higher than current directory.
###Now run this script to copy subsequent PR runs of .nc files and catenate them
### totpe= ' total number of processors', iteration='parareal iteration'
set totpe = 32
set iteration = 4
set pe = 2
while ( $pe <= $totpe )
##	set k = $pe
	set k = $iteration
	while ( $k >= 1 )
##	        if ($pe <= 9 ) then
        	if (-e ../b2time_fine.000"$k".000"$pe".nc) then 
			cp ../b2time_fine.000"$k".000"$pe".nc .
			/afs/ipp-garching.mpg.de/home/d/dpc/bin/nccat b2time.nc b2time_fine.000"$k".000"$pe".nc
			set k = 1
        	else if (-e ../b2time_fine.000"$k".00"$pe".nc) then 
			cp ../b2time_fine.000"$k".00"$pe".nc .
			/afs/ipp-garching.mpg.de/home/d/dpc/bin/nccat b2time.nc b2time_fine.000"$k".00"$pe".nc
			set k = 1
       		else if (-e ../b2time_fine.00"$k".00"$pe".nc) then 
			cp ../b2time_fine.00"$k".00"$pe".nc .
			/afs/ipp-garching.mpg.de/home/d/dpc/bin/nccat b2time.nc	b2time_fine.00"$k".00"$pe".nc
			set k = 1
        	else if (-e ../b2time_fine.00"$k".000"$pe".nc) then 
			cp ../b2time_fine.00"$k".000"$pe".nc .
			/afs/ipp-garching.mpg.de/home/d/dpc/bin/nccat b2time.nc b2time_fine.00"$k".000"$pe".nc
			set k = 1
		endif
##		endif	
		@ k = $k - 1
	end 		
	echo "Welcome $pe times"
	@ pe = $pe + 1
end

