#!/bin/bash

./runExp.sh -p -i ./workloads -t "1 2 4 8" -x "C E A" -c 1000000\
         -d "2 4 6 8 10 12 14 16 18 20 22 24"\
         -s "bskytree hybrid"
         