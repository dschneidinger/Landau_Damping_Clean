#!/bin/bash
osirispath=~/osiris/osiris-1.0.0
cd ${osirispath}

for n in $*;
do
    echo Starting simulation for k = $n
    cp ~/vscode/Landau_damping_clean/vth_0.055_amp_0.002/k_$n/k_$n.1d ${osirispath}/input_file_test.txt
    ./config/docker/osiris mpirun -n 10 bin/osiris-1D.e input_file_test.txt
    mv -f MS TIMINGS HIST run-info ~/vscode/Landau_damping_clean/vth_0.055_amp_0.002/k_$n/
done