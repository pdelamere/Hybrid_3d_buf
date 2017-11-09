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
echo "Setting ulimits"
ulimit -s unlimited
ulimit -f unlimited
ulimit -l unlimited

echo "Setting up modules"
module purge
module load slurm
module load lang/Anaconda3/2.5.0

echo "Generate machinefile"
srun -l /bin/hostname | sort -n | awk '{print $2}' > ./nodes.$SLURM_JOB_ID
echo "Start Hybrid Code"
mpirun -np $SLURM_NTASKS --machinefile ./nodes.$SLURM_JOB_ID "hybrid" > output 2> error
RESULT=$?
echo "Removing the machinefile"
rm ./nodes.$SLURM_JOB_ID

echo "Finish batch script"
exit $RESULT
"""

def gen_batch(filename, jobname, partition, ntasks, tasks_per_node=None, time=None, test=False):
    extra_batch_commands = ""
    if tasks_per_node is not None:
        extra_batch_commands += "#SBATCH --tasks-per-node={}\n".format(tasks_per_node)
    if time is not None:
        extra_batch_commands += "#SBATCH --time={}\n".format(time)
    if test:
        extra_batch_commands += "#SBATCH --test-only\n"


    with open(filename, mode='w') as f:
        f.write(shebang
                + batch_commands.format(jobname=jobname, partition=partition, ntasks=ntasks)
                + extra_batch_commands
                + script)

if __name__ == '__main__':
    gen_batch('test.slurm', 'pluto', 'debug', 115)

