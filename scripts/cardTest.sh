#!/bin/bash

./runExp.sh -p -i ./workloads -t "1 2 4 8" -x "C E A" -d 12\
         -c "500000 1000000 2000000 4000000 8000000"\
         -s "bskytree hybrid"
