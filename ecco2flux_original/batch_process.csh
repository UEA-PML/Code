#!/bin/csh -f
set picarro_files = `find /home/scratch/data/workspace/ECCO2Flux/Picarro -name \*.txt | sort -u `
set LPMS_files = `find /home/scratch/data/workspace/ECCO2Flux/LPMS -name \*.dat | sort -u`
set licor_files = `find /home/scratch/data/workspace/ECCO2Flux/Licor -name \*.txt | sort -u`
set sysdon_files = `find /home/scratch/data/workspace/ECCO2Flux/SysDon -name \*.txt | sort -u`
set seatex_dir = /home/scratch/data/workspace/ECCO2Flux/seatex
set gyro_file = /home/scratch/data/workspace/ECCO2Flux/Gyro/gyro_20181018.Raw
set underway_file = /home/scratch/data/workspace/ECCO2Flux/oceanlogger/oceanlogger_20181018.Raw

set counter = 1
foreach wfile(/home/scratch/data/workspace/ECCO2Flux/Wind/*COM12.txt)
    echo $counter
    ./process_ECCO2Flux.csh -w $wfile -p $picarro_files[$counter] -s $sysdon_files[$counter] \
    -l $licor_files[$counter] -a $seatex_dir -g $gyro_file -u $underway_file
    @ counter++
end

exit(0)
