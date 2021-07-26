#!/bin/bash

##Setup environment
source /cvmfs/csee.gsi.de/bin/go4login
export ROOT_INCLUDE_PATH=/lustre/gamma/DESPEC_S460_NEARLINE
echo "DESPEC Kronos Started at `date`"

##Set data location
#dpath=~/lustre/gamma/d004/ts/aida/

##Read in list of files to run. Format names seperated by space,tab,newline
LISTFILE="/lustre/gamma/DESPEC_S460_NEARLINE/Cluster_Submission/files_list_61.txt"



##Count number of files
NFILES=$(cat ${LISTFILE} | wc -l)
echo "Analysing" $NFILES "Files"

##Read names from list file
declare -a array
while IFS= read -r line
do
    array+=($line)
done < "$LISTFILE"

echo "Array is $SLURM_ARRAY_TASK_ID"
part=(  "${array[@]:$SLURM_ARRAY_TASK_ID:1}" ) # :5 number of files to put together -> Has to be the same in the 2 .sh scripts

echo "Running Go4!"



go4analysis -file ${part[*]} -enable-asf 12000 -asf /lustre/gamma/DESPEC_S460_NEARLINE/Cluster_Submission/Nearline_Histograms/NearlineS460_all_PIDgates_f61_new_$SLURM_ARRAY_TASK_ID.root


