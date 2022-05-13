#!/bin/bash -l

#bash start.sh cell_type

list=./list.txt


cat ${list}|while read line
do
    echo qsub star.qsub ${line} 
    eval qsub star.qsub ${line} 
done