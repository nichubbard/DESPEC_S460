#!/bin/bash
##Read in list of files to run
LISTFILE="files_list_61.txt"
declare -a size
while IFS= read -r line
do
    size+=($line)
done < "$LISTFILE"

##Submit job
sbatch -J despec_go4_61_s460 -D /lustre/gamma/DESPEC_S460_NEARLINE/ -o logs/go4_%A_%a.out.log -e logs/go4_%A_%a.err.log \
	--time=8:00:00 --mem-per-cpu=4G --array=0-${#size[@]}:1 -- ./Cluster_Submission/go4_launcher_nearline.sh

