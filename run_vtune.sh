#!/bin/bash

resultsdir="./vtune_results$1/"
rm -r $resultsdir
./fix_ptrace_scope.sh vtune -collect hotspots -r $resultsdir -- mpirun -np 4 ./src/simple_house2 input_file.txt > fout_vtune

# to see results in GUI: vtune-gui ./vtune_results/vtune.vtune
