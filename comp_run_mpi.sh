#!/bin/bash
# Compiling using OpenMPI.
# $1, $2 and so on are the cmd line args in order
# $@ is the array of all cmd line args
# $# is the size of $@

if [ "${@: -2: 1}" != "0" ]; then
mycmd="mpic++ "
N=$#
echo "Compiling with OpenMPI. Files included:"
for n in ${@:1:$#-2}
do
    echo $n
    # mycmd+="$PWD/${n} "
    mycmd+="${n} "
done
mycmd+="-o outfile.exe"
echo "$mycmd"
eval $mycmd
printf "\n"
fi

printf "Running the output file with ${@: -1} processes...\n"
mycmd2="mpiexec -n ${@: -1} outfile.exe"
# echo "$mycmd2"
eval $mycmd2