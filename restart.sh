#!/bin/bash
usage() { 
    echo "Usage: "
    echo "$0 [-h]"
    echo "$0 [-m <message>]"
    echo "-h                :   Print this help message and exit"
    echo "-m <message>      :   Provide a message to be saved in the data folder"
    exit 10
}

# Default values
unset MESSAGE
while getopts ":hm:" opt; do
    case $opt in
        m)# Save a message in the data folder
            MESSAGE="$OPTARG"
            ;;
        \?|h)
            usage
            ;;
    esac
done
shift $((OPTIND-1))
if [[ "$@" -ne "" ]]; then
    usage
fi

if [ -e "mpipath" ]; then
    MPI="$(cat mpipath)"/mpirun
else
    MPI="mpirun"
fi

date +"%Y-%m-%d-%T" >> restarts
if ! [ -z ${MESSAGE+x} ]; then
    echo "\t$MESSAGE" >> restarts
fi

cp --backup=numbered output old-output
cp --backup=numbered error old-error
"$MPI" ./hybrid restart > output 2> error
