#!/bin/bash
# Script for running the hybrid code.
# Ensures that no data files are overwritten and that the exact version
# of the hybrid code being used is stored along with its output.
COMMAND_LINE="${0##*/} $@"
usage() { 
    echo "Usage: "
    echo "$0 [-h]"
    echo "$0 [-i][-p <mpi-path>][-n][-d <data-folder>][-f <flags>][-m <message>] <num-proc>"
    echo "-h                :   Print this help message and exit"
    echo "-i                :   Ignore uncommited changes in the working directory"
    echo "-p <mpi-path>     :   Specify the path to the mpi installation directory"
    echo "-n                :   Disable recompiling the program (nocompile)"
    echo "-d <data-folder>  :   Specify the path where data will be stored"
    echo "-f <flags>        :   Override the default compiler flags"
    echo "-m <message>      :   Provide a message to be saved in the data folder"
    exit 10
}
# Default values
IGNORE=false
MPI=""
BUILD=true
MAIN_DATA="$HOME/data"
unset FFLAGS
unset MESSAGE
while getopts ":hip:nd:f:m:" opt; do
    case $opt in
        i)# Ignore the fact that the changes are not commited
            IGNORE=true
            ;;
        p)# Specify the path to mpi
            MPI="${OPTARG}/bin/"
            ;;
        n)# Disable recompiling the program
            BUILD=false
            ;;
        d)# Set folder to hold data
            MAIN_DATA="$OPTARG"
            ;;
        f)# Set flags for the compiler. Overrides all defaults
            FFLAGS="$OPTARG"
            ;;
        m)# Save a message in the data folder
            MESSAGE="$OPTARG"
            ;;
        \?|h)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

NUM_PROC=$1
re='^[0-9]+$'
if ! [[ "$NUM_PROC" =~ $re ]]; then
    printf "Invalid number of processors\n"
    usage
fi

# Check if working directory is clean
# if not, exit with error message.
if [ "$IGNORE" = false ]; then
    git diff-index --quiet HEAD -- ||\
        { printf "Working directory has uncommitted changes.\nPlease commit and retry.\n"; exit 1; }
fi

# Make sure the program is rebuilt correctly.
if [ "$BUILD" = true ]; then
    make clean
    if [ -z ${FFLAGS+x} ]; then # Check if FFLAGS is set
        make FC="${MPI}mpif90" || { printf "\nBuild failed, aborting\n"; exit 2; }
    else
        make FC="${MPI}mpif90" FFLAGS="$FFLAGS" || { printf "\nBuild failed, aborting\n"; exit 2; }
    fi

fi

# Make a folder to save all the data. Error if it already exists.
DATE="$(date +"%Y-%m-%d")"
TIME="$(date +"%T")"
mkdir -p "$HOME/data/$DATE" || { printf "There was a problem making the data folder.\n"; exit 3; }
DATA_FOLDER="$HOME/data/$DATE/pluto.$TIME"
mkdir "$DATA_FOLDER" || { printf "There was a problem making the folder for this run.\n"; exit 4; }

# Copy required files into the new data folder
cp hybrid "$DATA_FOLDER/hybrid" || { printf "Error while copying executable\n"; exit 5; }
cp inputs.dat "$DATA_FOLDER/inputs.dat" || { printf "Error while copying inputs.dat\n"; exit 6; }
cp fileShrinker.py "$DATA_FOLDER/fileShrinker.py" || { printf "Error while copying fileShrinker.py\n"; exit 7; }

printf "$COMMAND_LINE" > "$DATA_FOLDER/invocation"

if ! [ -z ${MESSAGE+x} ]; then
    printf "$MESSAGE" > "$DATA_FOLDER/message"
fi

# Record the version of the code being used.
if [ "$IGNORE" = false ]; then
    git rev-parse --verify HEAD > "$DATA_FOLDER/version.sha1" || { printf "Error getting the git commit sha1.\n"; exit 8; }
fi

# Finally, run the program from the data folder
cd "$DATA_FOLDER"
eval "${MPI}mpirun -n $NUM_PROC $DATA_FOLDER/hybrid > $DATA_FOLDER/output 2> $DATA_FOLDER/error &"

printf "Done\n"
