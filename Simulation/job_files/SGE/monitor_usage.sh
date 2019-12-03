#!/bin/bash
for i in {1..180..1}
do 
qhost -h devel-compute-5-0 >> usage_stats.txt
sleep 30
done

