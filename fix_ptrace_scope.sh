#!/bin/bash

val_orig=`cat /proc/sys/kernel/yama/ptrace_scope`
sudo sysctl -w kernel.yama.ptrace_scope=0
$@
sudo sysctl -w kernel.yama.ptrace_scope=$val_orig
