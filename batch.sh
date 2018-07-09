#!/bin/bash
# runs a script supplied by the user on a list of samples
# submitting each job to the job schedule
# this version works on CGRB (Oregon State University)
# usage: batch.sh script list jobname

SCRIPT=$1
LIST=$2
NAME=$3

exec < "$LIST"
while read LINE
do
        echo "$LINE"
	cd $LINE
	rm -rf $NAME.$LINE
	SGE_Batch -c "../$SCRIPT $LINE" -r $NAME.$LINE
	cd ../
done

