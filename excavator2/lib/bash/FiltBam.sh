#!/bin/bash

# exit when any command fails
set -e

################################################
#
# .bam file filtering 
#
################################################

Bam_File=$1
MAPQ=$2
Output_Folder=$3
Program_Folder=$4
Sample_Name=$5
Target_Folder=$6

mkdir -p "$Output_Folder/.tmp"

mychr="$Target_Folder/*_chromosome.txt"
mychr=$(< $mychr)

for i in $mychr
 do
  samtools view -F 1028 "$Bam_File" "$i" \
  | cut -f2,3,4,5 \
  | cut -f2,3,4 \
  | awk '$3 == 0 || $3 >= "/$MAPQ/"' \
  | perl -lane 'print "$F[0]\t$F[1]"' \
  > "$Output_Folder/.tmp/.FilteredBamtmp.000" # Shouldn't it be a progresive number like PID???
  R --slave --args \
    "$Sample_Name" \
	  "$Output_Folder" \
	  "$i" \
	  "$Program_Folder" \
    "$Target_Folder" \
  < "$Program_Folder/lib/R/MakeReadCount.R"
 done
