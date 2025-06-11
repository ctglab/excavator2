#!/bin/bash

# exit when any command fails
set -e

############################################################################

outputFolder="$7"
allchr=$(cat "$outputFolder"/*_chromosome.txt)

############################################################################  	if output files already exist print warning

if [[ -e "$outputFolder"/MAP/*RData ]]; then 
  echo 'Cleaning old files from  Mappability folder!' >&2
  rm -f "$outputFolder"/MAP/*RData
fi

if [[ -e "$outputFolder"/GCC/*RData ]]; then 
	echo 'Cleaning old files from GC-Content folder!' >&2
	rm -f "$outputFolder"/GCC/*RData
fi

############################################################################  	fasta file indexing if needed

if [[ ! -e "$6.fai" ]]; then 
	samtools faidx "$6"
fi

check=$(echo $allchr | grep "chr" | wc -c)
if (($check == 0)); then
  for i in $allchr
   do
    grep -w "$i" "$2" | awk 'BEGIN {FS="\t"} {OFS="\t"}; {print "chr"$1,$2-1,$3,$4}' > "$outputFolder/.temp.bed"
    bigWigAverageOverBed "$1" "$outputFolder/.temp.bed" "$outputFolder/MAP/Mapout.txt"
  	R --slave --args "$i" "$outputFolder/MAP" < "$3/lib/R/SaveMap.R"
    rm -f "$outputFolder/.temp.bed"
    rm -f "$outputFolder/MAP/Mapout.txt"
    grep -w "$i" "$2" | awk 'BEGIN {FS="\t"} {OFS="\t"}; {print "chr"$1,$2-1,$3,$4}' > "$outputFolder/.temp.bed"
    bedtools nuc -fi "$6" -bed "$outputFolder/.temp.bed" | cut -f6 | sed '1d' > "$outputFolder/GCC/GCC.txt"
    R --slave --args "$i" "$outputFolder/GCC" < "$3/lib/R/SaveGCC.R"
    rm -f "$outputFolder/.temp.bed"
    rm -f "$outputFolder/GCC/GCC.txt"
    grep -w "$i" "$2" | awk 'BEGIN {FS="\t"} {OFS="\t"}; {print "chr"$1,$2,$2+1,$4}' > "$outputFolder/.temp.bed"
    fastaFromBed -fi "$6" -bed "$outputFolder/.temp.bed" -fo "$outputFolder/FRB/FRB.txt"
    R --slave --args "$i" "$outputFolder/FRB" < "$3/lib/R/SaveFRB.R"
    rm -f "$outputFolder/.temp.bed"
    rm -f "$outputFolder/FRB/FRB.txt"
    echo "Processed chromosome $i" >&2
   done
else
  for i in $allchr
   do
    grep -w "$i" "$2" | awk 'BEGIN {FS="\t"} {OFS="\t"}; {print $1,$2-1,$3,$4}' > "$outputFolder/.temp.bed"
    bigWigAverageOverBed "$1" "$outputFolder/.temp.bed" "$outputFolder/MAP/Mapout.txt"
    R --slave --args "$i" "$outputFolder/MAP" < "$3/lib/R/SaveMap.R"
    bedtools nuc -fi "$6" -bed "$outputFolder/.temp.bed" | cut -f6 | sed '1d' > "$outputFolder/GCC/GCC.txt"
    R --slave --args "$i" "$outputFolder/GCC" < "$3/lib/R/SaveGCC.R"
    rm -f "$outputFolder/.temp.bed"
    rm -f "$outputFolder/GCC/GCC.txt"
    rm -f "$outputFolder/MAP/Mapout.txt"
    grep -w "$i" "$2" | awk 'BEGIN {FS="\t"} {OFS="\t"}; {print $1,$2,$2+1,$4}' > "$outputFolder/.temp.bed"
    fastaFromBed -fi "$6" -bed "$outputFolder/.temp.bed" -fo "$outputFolder/FRB/FRB.txt"
    R --slave --args "$i" "$outputFolder/FRB" < "$3/lib/R/SaveFRB.R"
    rm -f "$outputFolder/.temp.bed"
    rm -f "$outputFolder/FRB/FRB.txt"
    echo "Processed chromosome $i" >&2
  done
fi
