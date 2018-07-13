#!/bin/bash
# organizes a list of paired reads in sample-specific directories
# input a list of sample names

LIST=$1
exec < "$LIST"
while read LINE
do
        echo "$LINE"
	mkdir $LINE
	mv $LINE*.fastq* $LINE
done

