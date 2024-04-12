#!/bin/bash
# Compiling using G++-13.
# mycmd="/home/linuxbrew/.linuxbrew/bin/g++-13 -fdiagnostics-color=always -std=c++20 -g "
mycmd="g++ "
N=$#
echo "Files included:"
for n in $@
do
    echo $n
    mycmd+="$PWD/${n} "
    # echo $mycmd$'\n'
done
mycmd+="-o outfile.exe"
# echo $mycmd
eval $mycmd