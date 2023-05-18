#!/bin/bash
#SBATCH -J NNtrain
#SBATCH -p C064M0256G
#SBATCH --nodes=4
##SBATCH --time=72:00:00 
#SBATCH --ntasks-per-node=64
#SBATCH --error=NN_%J.err          ###err log###
#SBATCH --output=NN_%J.out         ###out log###

module load LASP

# gen machinefile
srun hostname -s | sort -n > slurm.hosts

EXEC=lasp 

mpirun -machinefile slurm.hosts $EXEC


scontrol show job $SLURM_JOBID



