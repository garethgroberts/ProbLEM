#!/bin/bash

list=`seq 1 10 | tr '\n' '\t'`

for i in $list; do
	echo "working on run" $i
	mkdir run_$i
	python ProbLEM_force_OU_evo_single_fig.py ; mv *.csv ./run_$i
	done