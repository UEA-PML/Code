#!/bin/csh -f
# Module: process_ECCO2Flux.csh
# Author: Tim Smyth
# Date: 05/11/2018
# Version: 1.0
#
# Description: a wrapper script to pass correct arguments to process_ECCO2Flux.py.
#              Modification of the files is required prior to using the python
#              script including:
#              1) conversion from dos to UNIX (dos2unix)
#              2) Ensuring there are no spurious pieces of data in the file
#              3) Removing unneccessary parts of the input (brackets, spaces etc.)
#
# USAGE: process_ECCO2Flux.csh <-w wind file> <-m LPMS file | -s sysdon file | -c CR6 file > <-p picarro file> [-l licor file] <-a seatex directory> <-g gyro file> [-ACO] [-match mm/dd/yyy]
#
# All switches and files in < > must be specified for the code to work

set USAGE = "process_ECCO2Flux.csh <-w wind file> <-m LPMS file | -s sysdon file | -c CR6 file > <-p picarro file> [-l licor file] <-a seatex directory> <-g gyro file> [-ACO] [-match mm/dd/yyy] [-u underway file]"

if ($HOST == "jrlc.jcr.nerc-bas.ac.uk") then
   set PYTHON = /users/autoflux/bin/python2.7
   set execdir = /users/autoflux/scripts
else
   set PYTHON = /bin/python
   set execdir = /users/rsg/tjsm/autoflux/scripts
endif

echo $PYTHON

set wind_flag = 0
set LPMS_flag = 0
set picarro_flag = 0
set licor_flag = 0
set sysdon_flag = 0
set seatex_flag = 0
set gyro_flag = 0
set match_flag = 0
set underway_flag = 0
set CR6_flag = 0
set ACO_flag = 0

while ("$1" =~ -*)
   switch ($1)
   case '-w':
      set wind_file = $2
      set wind_flag = 1
      if !(-e $wind_file) then
         echo "Wind file not found"
         echo $USAGE
         exit(-1)
      else 
         echo "Using wind file: $wind_file"
      endif 
      shift
      breaksw
   case '-m':
      set LPMS_file = $2
      set LPMS_flag = 1
      if !(-e $LPMS_file) then
         echo "LPMS file not found"
         echo $USAGE
         exit(-1)
      else 
         echo "Using LPMS file: $LPMS_file"
      endif 
      shift
      breaksw
   case '-c':
      set CR6_file = $2
      set CR6_flag = 1
      if !(-e $CR6_file) then
         echo "CR6 file not found"
         echo $USAGE
         exit(-1)
      else 
         echo "Using CR6 file: $CR6_file"
      endif 
      shift
      breaksw
   case '-p':
      set picarro_file = $2
      set picarro_flag = 1
      if !(-e $picarro_file) then
         echo "picarro file not found"
         echo $USAGE
         exit(-1)
      else 
         echo "Using picarro file: $picarro_file"
      endif 
      shift
      breaksw
   case '-l':
      set licor_file = $2
      set licor_flag = 1
      if !(-e $licor_file) then
         echo "licor file not found"
         echo $USAGE
         exit(-1)
      else 
         echo "Using licor file: $licor_file"
      endif 
      shift
      breaksw
   case '-s':
      set sysdon_file = $2
      set sysdon_flag = 1
      if !(-e $sysdon_file) then
         echo "sysdon file not found"
         echo $USAGE
         exit(-1)
      else 
         echo "Using sysdon file: $sysdon_file"
      endif 
      shift
      breaksw
   case '-a':
      set seatex_dir = $2
      set seatex_flag = 1
      if !(-e $seatex_dir) then 
         echo "seatex directory not found"
         echo $USAGE
         exit(-1)
      else
         # Find all of the necessary files in the directory
         set seatex_files = `ls $seatex_dir/seatex-{gll,gga,vtg}*`
         echo "Using seatex files: $seatex_files"
      endif
      shift
      breaksw
   case '-g':
      set gyro_file = $2
      set gyro_flag = 1
      if !(-e $gyro_file) then
         echo "Gyro file not found"
         echo $USAGE
         exit(-1)
      else
         echo "Using gyro file: $gyro_file"
      endif
      shift
      breaksw
   case '-u':
      set underway_file = $2
      set underway_flag = 1
      if !(-e $underway_file) then
         echo "Underway file not found"
         echo $USAGE
         exit(-1)
      else
         echo "Using underway file: $underway_file"
      endif
      shift
      breaksw
   case '-match':
      set match_date = $2
      echo "Using $match_date for Navigation files"
      set match_flag = 1
      shift
      breaksw
   case '-ACO':
      set ACO_flag = 1
      breaksw
   default:
      echo $USAGE
      exit(-1)
      breaksw
   endsw
   shift
end

if ($wind_flag == 0 || $picarro_flag == 0 || ($LPMS_flag == 0 && $sysdon_flag == 0 && CR6_flag == 0) || $seatex_flag == 0) then
   echo "++++++++++++++++++++++++++++++"
   echo "Need to correctly specify files and switches as follows:"
   echo "" 
   echo $USAGE
   echo "++++++++++++++++++++++++++++++"
   exit(-1)
endif

echo ""
echo " Modifying files for input to process_ECCO2Flux.py"
if ($wind_flag) then 
   cat $wind_file | sed -e 's/\[//g' -e 's/\]//g' -e 's/M\:x\ \=//g' -e 's/H\:x\ \=//g' -e 's/y\ \=//g' -e 's/z\ \=//g' -e 's/t\ \=//g' \
                  | /usr/bin/tr -cd '\r\11\12\15\40-\176' | /usr/bin/tr -d '\r'  | /usr/bin/tr -s " " \
                  | awk 'BEGIN {FS = " "}; {if (NF==9) {print $0}}' \
                  > /tmp/process_ecco2flux_wind.txt
endif

if ($LPMS_flag) then 
   cat $LPMS_file | /usr/bin/tr -cd '\r\11\12\15\40-\176' | /usr/bin/tr -d '\r'  | /usr/bin/tr -s " " \
                  | awk 'BEGIN {FS = ","}; {if (NF==13) {print $0}}' \
                  > /tmp/process_ecco2flux_LPMS.txt
endif

cat $picarro_file | sed -e 's/\[//g' -e 's/\]//g' -e 's/M\:x\ \=//g' -e 's/y\ \=//g' -e 's/z\ \=//g' -e 's/t\ \=//g' \
                  | /usr/bin/tr -cd '\r\11\12\15\40-\176' | /usr/bin/tr -d '\r'  | /usr/bin/tr -s " " \
                  | awk 'BEGIN {FS = " "}; {if (NF==16) {print $0} if (NF==22 && $1!~/DATE/) {print strftime("%a %b %d",$6),$2,strftime("%Y",$6),$6,$9,$10,$17,$19,$20,$21,$22}}' \
                > /tmp/process_ecco2flux_picarro.txt

if ($licor_flag) then
   cat $licor_file | sed -e 's/\[//g' -e 's/\]//g' -e 's/[()]/\ /g' \
                   | /usr/bin/tr -cd '\r\11\12\15\40-\176' | /usr/bin/tr -d '\r'  | /usr/bin/tr -s " " \
                   | awk 'BEGIN {FS = " "}; {if (NF==26) {print $0}}' \
                   > /tmp/process_ecco2flux_licor.txt
endif

if ($sysdon_flag) then
   cat $sysdon_file | /usr/bin/tr -cd '\r\11\12\15\40-\176' | /usr/bin/tr -d '\r'  | /usr/bin/tr -s " " \
                    | awk 'BEGIN {FS = " "}; {if (NF==8) {print $0}}' \
                    > /tmp/process_ecco2flux_SysDon.txt
endif

if ($CR6_flag) then
   cat $CR6_file | /usr/bin/tr -cd '\r\11\12\15\40-\176' | /usr/bin/tr -d '\r'  | /usr/bin/tr -s " " \
                 | /usr/bin/tr -d '"' | tail -n +5 | sed -e 's/^[ \t]*//' \
                 > /tmp/process_ecco2flux_CR6.txt
endif

# Grab the specific date from the Picarro file
set filelength = `wc -l /tmp/process_ecco2flux_picarro.txt | awk '{print int($1/2)}'`
set ofile_stem = `head -$filelength /tmp/process_ecco2flux_picarro.txt | tail -1 | awk 'BEGIN {FS = " "}; {print strftime("%Y%m%d",$6)}'`
set obinaryfile_stem = `head -1 /tmp/process_ecco2flux_picarro.txt | awk 'BEGIN {FS = " "}; {print strftime("%Y%m%d-%H%M%S",$6)}'`

if ($match_flag == 0) then 
   # Determine the length of the file and then get the date from the mid-point
   set match_date = `head -$filelength /tmp/process_ecco2flux_picarro.txt |tail -1 | awk 'BEGIN {FS = " "}; {print strftime("%m/%d/%Y",$6)}'`
endif

echo $match_date

if ($ACO_flag == 0) then
   cat $seatex_dir/seatex-gll* | grep ${match_date} \
                      | /usr/bin/tr -cd '\r\11\12\15\40-\176' | /usr/bin/tr -d '\r'  | /usr/bin/tr -s " " \
                      | awk 'BEGIN {FS = ","}; {if (NF==10) {print $1,$2,$4,$5,$6,$7}}' \
                      > /tmp/process_ecco2flux_seatex-gll.txt

   cat $seatex_dir/seatex-gga* | grep ${match_date} \
                      | /usr/bin/tr -cd '\r\11\12\15\40-\176' | /usr/bin/tr -d '\r'  | /usr/bin/tr -s " " \
                      | awk 'BEGIN {FS = ","}; {if (NF==14) {print $1,$2,$4,$5,$6,$7}}' \
                      > /tmp/process_ecco2flux_seatex-gga.txt

   cat $seatex_dir/seatex-vtg* | grep $match_date \
                   | /usr/bin/tr -cd '\r\11\12\15\40-\176' | /usr/bin/tr -d '\r'  | /usr/bin/tr -s " " \
                   | awk 'BEGIN {FS = ","}; {if (NF==12) {print $1,$2,$4,$8}}' \
                   > /tmp/process_ecco2flux_seatex-vtg.txt

   cat $gyro_file | grep ${match_date} \
                  | /usr/bin/tr -cd '\r\11\12\15\40-\176' | /usr/bin/tr -d '\r'  | /usr/bin/tr -s " " \
                  | awk 'BEGIN {FS = ","}; {if (NF==5) {print $1,$2,$4}}' \
                  > /tmp/process_ecco2flux_gyro.txt
               
   cat $underway_file | grep ${match_date} \
                  | /usr/bin/tr -cd '\r\11\12\15\40-\176' | /usr/bin/tr -d '\r'  | /usr/bin/tr -s " " \
                  | awk 'BEGIN {FS = ","}; {if (NF==25) {print $1,$2,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$17}}' \
                  > /tmp/process_ecco2flux_oceanlogger.txt             
else
   set MM = `echo $match_date | cut -c1-2`
   set DD = `echo $match_date | cut -c4-5`
   set YYYY = `echo $match_date | cut -c7-10`
   set jmatch_date = `julian $DD $MM $YYYY`
   set aco_match = `echo $jmatch_date[2]","$jmatch_date[1]"."`
   echo $aco_match

   cat $seatex_dir/seatex-gga* | grep ${aco_match} \
                      | /usr/bin/tr -cd '\r\11\12\15\40-\176' | /usr/bin/tr -d '\r'  | /usr/bin/tr -s " " \
                      | awk 'BEGIN {FS = ","}; {if (NF==13) {print $1,$2,$6,$7}}' \
                      > /tmp/process_ecco2flux_seatex-gga.txt

   cat $seatex_dir/seatex-vtg* | grep ${aco_match} \
                   | /usr/bin/tr -cd '\r\11\12\15\40-\176' | /usr/bin/tr -d '\r'  | /usr/bin/tr -s " " \
                   | awk 'BEGIN {FS = ","}; {if (NF==8) {print $1,$2,$5,$7}}' \
                   > /tmp/process_ecco2flux_seatex-vtg.txt

   cat $gyro_file | grep ${aco_match} \
                  | /usr/bin/tr -cd '\r\11\12\15\40-\176' | /usr/bin/tr -d '\r'  | /usr/bin/tr -s " " \
                  | awk 'BEGIN {FS = ","}; {if (NF==5) {print $1,$2,$5}}' \
                  > /tmp/process_ecco2flux_gyro.txt
               
   cat $underway_file | grep ${aco_match} \
                  | /usr/bin/tr -cd '\r\11\12\15\40-\176' | /usr/bin/tr -d '\r'  | /usr/bin/tr -s " " \
                  | awk 'BEGIN {FS = ","}; {if (NF==25) {print $1,$2,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$17,$18}}' \
                  > /tmp/process_ecco2flux_oceanlogger.txt             
endif


echo "  Running python script process_ECCO2Flux.py"

if ($ACO_flag) then
   set ship_nav_option = "--ACO"
else
   set ship_nav_option = ""
endif

#### GOT TO HERE - need to read in the CR6 data and change the way the ships nav is read in


if ($LPMS_flag == 1 && $sysdon_flag == 0 && $CR6_flag == 0) then 
   if ($licor_flag == 1) then 
      echo "Processing with LPMS motion and Licor"
      $PYTHON process_ECCO2Flux.py $ship_nav_option $wind_option --underway --LPMS --ofile $execdir/output/${ofile_stem}_fluxdiagnostics.txt  --L0 $execdir/output/L0/${obinaryfile_stem}_L0.pkl
   else
      echo "Processing with LPMS motion and no Licor"
      $PYTHON process_ECCO2Flux.py $ship_nav_option $wind_option --underway --LPMS --no_licor --ofile $execdir/output/${ofile_stem}_fluxdiagnostics.txt --L0 $execdir/output/L0/${obinaryfile_stem}_L0.pkl
   endif   
endif

if ($sysdon_flag == 1 && $LPMS_flag == 0 && $CR6_flag == 0) then
   if ($licor_flag == 1) then 
      echo "Processing with SysDon motion and Licor"
      $PYTHON process_ECCO2Flux.py $ship_nav_option $wind_option --underway --ofile $execdir/output/${ofile_stem}_fluxdiagnostics.txt --L0 $execdir/output/L0/${obinaryfile_stem}_L0.pkl
   else
      echo "Processing with SysDon motion and no Licor"
      $PYTHON process_ECCO2Flux.py $ship_nav_option $wind_option --underway --no_licor --ofile $execdir/output/${ofile_stem}_fluxdiagnostics.txt --L0 $execdir/output/L0/${obinaryfile_stem}_L0.pkl
   endif   
endif

if ($CR6_flag == 1 && $LPMS_flag == 0 && $sysdon_flag == 0) then 
   echo "This is where we will be executing for the CR6 files"
   if ($licor_flag == 1) then 
      echo "Processing with CR6 motion and Licor"
      echo "$PYTHON process_ECCO2Flux.py $ship_nav_option $wind_option --underway --CR6 --ofile $execdir/output/${ofile_stem}_fluxdiagnostics.txt --L0 $execdir/output/L0/${obinaryfile_stem}_L0.pkl"
      $PYTHON process_ECCO2Flux.py $ship_nav_option $wind_option --underway --CR6 --ofile $execdir/output/${ofile_stem}_fluxdiagnostics.txt --L0 $execdir/output/L0/${obinaryfile_stem}_L0.pkl
   else
      echo "Processing with CR6 motion and no Licor"
      $PYTHON process_ECCO2Flux.py $ship_nav_option $wind_option --underway --CR6 --no_licor --ofile $execdir/output/${ofile_stem}_fluxdiagnostics.txt --L0 $execdir/output/L0/${obinaryfile_stem}_L0.pkl
   endif   


endif

exit(0)

