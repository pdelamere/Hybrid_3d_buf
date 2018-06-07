#!/usr/local/pkg/python-anaconda/anaconda3/bin/python

shebang = "#!/bin/bash\n"

batch_commands = """
#SBATCH --job-name={jobname}
#SBATCH --partition={partition}
#SBATCH --ntasks={ntasks}
#SBATCH --mail-user=npbarnes@alaska.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
"""

script = """
if [ $# -eq 0 ]; then
    echo "This is not a restart"
    RESTART=false
else
    echo "This is a restart"
    RESTART=true
    PREVIOUS_JOBID=$1
    if [ $PREVIOUS_JOBID = "manual" ]; then
        echo "Previous job restarted manually"
    else
        echo "Checking if the previous job failed or timed out"
        if ! sacct -j $PREVIOUS_JOBID --format=State |grep -q TIMEOUT; then
            echo "Previous job ($PREVIOUS_JOBID) failed"
            echo "Previous job ($PREVIOUS_JOBID) failed" > failed
            exit 1
        else
            echo "Previous job ($PREVIOUS_JOBID) timed out"
        fi
    fi
fi

date +"%Y-%m-%d-%T" >> restarts

if [ "$RESTART" = true ] ; then
    if [ $PREVIOUS_JOBID = "manual" ]; then
        echo -e "\\tManual restart" >> restarts
    else
        echo -e "\\tAutomatic restart" >> restarts
    fi
else
    echo -e "\\tInitial start" >> restarts
fi
    
echo "Submit restart job"
sbatch -d afternotok:$SLURM_JOBID pluto.slurm $SLURM_JOBID &> restart_job_submission
sleep 1
NEXT_JOBID=$(cat restart_job_submission |cut -f 4 -d " ")

echo "Setting ulimits"
ulimit -s unlimited
ulimit -f unlimited
ulimit -l unlimited

echo "Setting up modules"
. /etc/profile.d/modules.sh
module purge
module load slurm
module load compiler/ifort/2016.3.210-GCC-5.4.0-2.26
module load openmpi/intel/3.0.0 
module load lang/Anaconda3/2.5.0

echo "Generate machinefile"
cd $SLURM_SUBMIT_DIR
srun -l /bin/hostname | sort -n | awk '{print $2}' > ./nodes.$SLURM_JOB_ID
echo "Start Hybrid Code"
if [ "$RESTART" = true ]; then
    mpirun -np $SLURM_NTASKS --machinefile ./nodes.$SLURM_JOB_ID hybrid restart > output.$SLURM_JOBID 2> error.$SLURM_JOBID
else
    mpirun -np $SLURM_NTASKS --machinefile ./nodes.$SLURM_JOB_ID hybrid > output.$SLURM_JOBID 2> error.$SLURM_JOBID
fi

RESULT=$?
echo "Removing the machinefile"
rm ./nodes.$SLURM_JOB_ID

echo "Cancel restart"
scancel $NEXT_JOBID

echo "Finish batch script"
exit $RESULT
"""

def gen_batch(filename, jobname, partition, ntasks, tasks_per_node=None, time=None, test=False, cont=False):
    extra_batch_commands = ""
    if tasks_per_node is not None:
        extra_batch_commands += "#SBATCH --tasks-per-node={}\n".format(tasks_per_node)
    if time is not None:
        extra_batch_commands += "#SBATCH --time={}\n".format(time)
    if test:
        extra_batch_commands += "#SBATCH --test-only\n"
    if cont:
        extra_batch_commands += "#SBATCH --contiguous\n"

    with open(filename, mode='w') as f:
        f.write(shebang
                + batch_commands.format(jobname=jobname, partition=partition, ntasks=ntasks)
                + extra_batch_commands
                + script)

if __name__ == '__main__':
    gen_batch('test.slurm', 'pluto', 'debug', 115, restart_filename='restart.slurm')

