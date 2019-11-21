#!/bin/bash -l

#The purpose of this wrapper script is to submit
#the map_and_dedup script using all files in file_names.txt
#Attempt to submit large job as an array. :)


files=$(awk 'END{print NR}' ./file_names.txt)
echo $files

#Submit kraken script as an array
#Submit script 1 to length of patient_array number of times.
qsub -t 1-$files:2 \
	-N NIPV_Dedup \
	./map_and_dedup_Bangledash_NEBNextV1.sh \

echo "Wrapper script done"

