#!/bin/bash
#SBATCH -J coor_forset
#SBATCH -p g1_shang 
#SBATCH --nodes=1
##SBATCH --time=72:00:00 
#SBATCH --ntasks-per-node=56
#SBATCH --error=%J.err          ###err log###
#SBATCH --output=%J.out         ###out log###

source /es01/software/profile.d/intel2018.sh

a=/es01/home/caoxiaoming/apps/lasp/INTER_pro3.3.4_intel12_2cw/
if [ -d $I_MPI_ROOT/intel64 ]; then
  b="$I_MPI_ROOT/intel64"
elif [ -d $I_MPI_ROOT/lib ]; then
  b="$I_MPI_ROOT"
fi
if [ ! -f libmpigf.so.4 ]; then ln -s $b/lib/libmpifort.so.12 libmpigf.so.4 ; fi
if [ ! -f libmpi_dbg.so.4 ] ; then ln -s $b/lib/debug_mt/libmpi.so.12 libmpi_dbg.so.4 ; fi
export LD_LIBRARY_PATH=$a:$LD_LIBRARY_PATH


unset I_MPI_PMI_LIBRARY
export I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=0

LIBVASP=${a}/Libraries/1.vasplib
export LD_LIBRARY_PATH=$LIBVASP:$LD_LIBRARY_PATH

EXEC=${a}/Src/lasp 

rm JobDone
echo 'Running' > JobRunning
echo `date ` >> JobRunning

ROOTDIR=$( pwd )
for i in AuPdO*
do

cd ${i}
mpirun -np $SLURM_NPROCS $EXEC
cd $ROOTDIR
done

echo 'Done' >> JobRunning
echo `date` >> JobRunning
mv JobRunning JobDone

echo "LASP" >> ~/finished
echo `date` >> ~/finished
echo $(pwd) >> ~/finished

scontrol show job $SLURM_JOBID



