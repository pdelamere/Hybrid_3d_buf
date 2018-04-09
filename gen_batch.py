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

restart_preamble = """
./restart.sh -m "Automatic restart"
if sacct -j $1 --format=State |grep -q TIMEOUT
then
"""

preamble = """
echo "Submit restart job"
sbatch -d afternotok:$SLURM_JOBID restart.slurm $SLURM_JOBID >& restart_job_submission
sleep 1
set NEWJOB=$(cat restart_job_submission |cut -f 4 -d " ")

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
"""

job_fmt = 'mpirun -np $SLURM_NTASKS --machinefile ./nodes.$SLURM_JOB_ID hybrid {} > output 2> error'
job = job_fmt.format('')
restart_job = job_fmt.format('restart')

conclusion = """
RESULT=$?
echo "Removing the machinefile"
rm ./nodes.$SLURM_JOB_ID

echo "Cancel restart"
scancel $NEWJOB

echo "Finish batch script"
exit $RESULT
"""

end_restart ="""
else
echo "Previous job failed" > failed
fi
"""


script = (preamble
        + job
        + conclusion)

restart_script = (restart_preamble
        + preamble
        + restart_job
        + conclusion
        + end_restart)

def gen_batch(filename, jobname, partition, ntasks, tasks_per_node=None, time=None, test=False, cont=False, restart_filename=None):
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

    if restart_filename is not None:
        with open(restart_filename, mode='w') as f:
            f.write(shebang
                    + batch_commands.format(jobname=jobname, partition=partition, ntasks=ntasks)
                    + extra_batch_commands
                    + restart_script)



if __name__ == '__main__':
    gen_batch('test.slurm', 'pluto', 'debug', 115, restart_filename='restart.slurm')

