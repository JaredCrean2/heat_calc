#!/bin/bash

rm -v ./output?
rm -v ./errfile?
rm -v ./errout?
mpirun --outfile-pattern output%r $@

tail ./output?
