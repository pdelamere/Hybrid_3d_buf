#!/bin/bash
# Script for running the hybrid code.
# Ensures that no data files are overwritten and that the exact version
# of the hybrid code being used is stored along with its output.

USAGE="Usage: run.sh [-i] <num-proc>\n"
if [[ $# -ne 1 ]]; then
    if [[ $# -ne 2 ]]; then
        printf "Wrong number of command line arguments.\n"
        printf $USAGE
        exit 10
    elif [[ $1 -ne "-i" ]]; then
        printf "Invalid option"
        printf $USAGE
        exit 10
    else # [[ $# -eq 2 && $1 -eq "-i" ]]; then
        IGNORE=true
        NUM_PROC=$2
    fi
else
    IGNORE=false
    NUM_PROC=$1
fi

re='^[0-9]+$'
if ! [[ $NUM_PROC =~ $re ]]; then
    printf "Invalid number of processors"
    printf $USAGE
    exit 10
fi

# Check if working directory is clean
# if not, exit with error message.
if [ "$IGNORE" = false ]; then
    git diff-index --quiet HEAD -- ||\
        { printf "Working directory has uncommitted changes.\nPlease commit and retry.\n"; exit 1; }
fi

# Make sure the program is rebuilt correctly.
printf "Attempting to rebuild the working directory from scratch.\n"
make clean
make || { printf "\nBuild failed, aborting\n"; exit 2; }

# Make a folder to save all the data. Error if it already exists.
mkdir -p $HOME/data/ || { printf "There was a problem making the data folder.\n"; exit 3; }
DATA_FOLDER=$HOME/data/pluto.$(date +"%Y-%m-%d-%T")
mkdir $DATA_FOLDER || { printf "There was a problem making the folder for this run.\n"; exit 4; }

# Copy required files into the new data folder
cp hybrid $DATA_FOLDER/hybrid || { printf "Error while copying executable\n"; exit 5; }
cp inputs.dat $DATA_FOLDER/inputs.dat || { printf "Error while copying inputs.dat\n"; exit 6; }
cp fileShrinker.py $DATA_FOLDER/fileShrinker.py || { printf "Error while copying fileShrinker.py\n"; exit 7; }

# Record the version of the code being used.
if [ "$IGNORE" = false ]; then
    git rev-parse --verify HEAD > $DATA_FOLDER/version.sha1 || { printf "Error getting the git commit sha1.\n"; exit 8; }
fi

# Finally, run the program from the data folder
cd $DATA_FOLDER
mpirun -n $NUM_PROC $DATA_FOLDER/hybrid > $DATA_FOLDER/output 2> $DATA_FOLDER/error &

echo Done
