#!/bin/bash

mpirun -np 2 -pmi-port  ~/scripts_local/term --wait -- lldb --source lldb_commands -- ./test/unit_tests --gtest_filter=DiscVectorTester.Sizes
