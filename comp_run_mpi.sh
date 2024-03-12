#!/usr/bin/bash
# Compiling using OpenMPI.
# $1, $2 and so on are the cmd line args in order
# $@ is the array of all cmd line args
# $# is the size of $@

mycmd="mpic++ "
N=$#
echo "Compiling with OpenMPI. Files included:"
for n in ${@:1:$#-1}
do
    echo $n
    mycmd+="$PWD/${n} "
done
mycmd+="-o outfile.exe"
# echo "$mycmd"
eval $mycmd

printf "\nRunning the output file with ${@: -1} processes...\n"
mycmd2="mpiexec -n ${@: -1} outfile.exe"
# echo "$mycmd2"
eval $mycmd2