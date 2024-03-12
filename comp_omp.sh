#!/usr/bin/bash
# Compiling using G++-13 including OpenMP.
mycmd="/home/linuxbrew/.linuxbrew/bin/g++-13 -fdiagnostics-color=always -std=c++20 -g -fopenmp "
N=$#
echo "OpenMP. Files included:"
for n in $@
do
    echo $n
    mycmd+="$PWD/${n} "
done
mycmd+="-o outfile.exe"
eval $mycmd