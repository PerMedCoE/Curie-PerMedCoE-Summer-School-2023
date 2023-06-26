#SBATCH --job-name=physiboss    # Job name
#SBATCH --output=jobname.%j.out # Stdout (%j expands to jobId)
#SBATCH --error=jobname.%j.err # Stderr (%j expands to jobId)
#SBATCH --ntasks=10     # Number of tasks(processes)
#SBATCH --nodes=1     # Number of nodes requested
#SBATCH --ntasks-per-node=10     # Tasks per node
#SBATCH --cpus-per-task=1     # Threads per task
#SBATCH --time=01:00:00   # walltime
#SBATCH --mem=10G   # memory per NODE
#SBATCH --partition=test
#SBATCH --account=project_2007898

# set the number of threads based on --cpus-per-task
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun ./project ./config/PhysiCell_settings.xml

