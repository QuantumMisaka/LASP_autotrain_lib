#!/bin/bash
#PBS -N LASP-autotrain
#PBS -r n
#PBS -j oe
#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=24
#PBS -q xmcao

source /data/intel-env/intel15u6.sh

EXEC=/home/apps/LASP/lasp
export LD_LIBRARY_PATH=/home/apps/LASP/INTER_pro3.2.0_intel12_2cw/Libraries/1.vasplib
EXEC=/home/james/apps/LASP/lasp_NN
export PYTHONPATH=/home/james/apps/anaconda3/bin:$PYTHONPATH
alias python=/home/james/apps/anaconda3/bin/python

POTNAME=H2O
MAXCYCLE=100

# go to workdir
cd $PBS_O_WORKDIR
if [ ! -f ~/.mpd.conf ]; then
/bin/echo "secretword=dfadfs" >> ~/.mpd.conf
/bin/chmod 600 ~/.mpd.conf
fi

NP= `cat $PBS_NODEFILE | wc -l`
echo "Numbers of Processors:;  $NP"
echo "----------------------------"
echo `date`

cat $PBS_NODEFILE | uniq # print the running node

# setup mpi env
# export OMP_NUM_THREADS=1
# export P4_GLOBMEMSIZE=1073741824
export I_MPI_PIN_DOMAIN=auto
export MPD_CON_EXT=$PBS_JOBID

echo "Job Starting at `date`" >> JobRunning

ROOTDIR= $PBS_O_WORKDIR

for((i=1;i<=$MAXCYCLE;i++))
do
# SSW-NN and VASP-DFT running
printf " Start SSW-DFT-NN Cycle %d\n" $i 
echo `date`

if [ $i = "1" ]; then
    python auto_local.py > output
else
    python auto_local.py >> output
fi
# VASP result backup
if [ -d $ROOTDIR/VASP/cycle-0 ]; then
  mv $ROOTDIR/VASP/cycle-0 $ROOTDIR/VASP/para-0
fi
if [ -d $ROOTDIR/VASP/cycle-1 ]; then
  mv $ROOTDIR/VASP/cycle-1 $ROOTDIR/VASP/para-$i
fi
# NN training job
cd NN
sleep 30
echo --- Start NN Training ---
echo `date`
Ntotal=$( grep Energy TrainStr.txt -c )
echo 'Num of Str/For Training data is:'
echo $Ntotal
sed -i '/^NNtrain/d' lasp.in
sed -i '/^Ntrain/d' lasp.in
sed -i '1a \Ntrain '$Ntotal lasp.in

# running NN training
mpirun -r ssh -np $NP $EXEC 2>&1 > output

sleep 20
# input is pot-format tmp file?
# should modify pot-name here
sed -i 's/newrun/'"${POTNAME}"'.input/g' ${POTNAME}.pot
cp ${POTNAME}.pot ${POTNAME}.pot-$i
cp ${POTNAME}.pot ${POTNAME}.input
cp lasp.out lasp.out-$i
cd $PBS_O_WORKDIR

printf " End SSW-DFT-NN Cycle %d\n" $i
echo `date`

done

echo "Job Finished at `date`" >> JobRunning
mv JobRunning JobDone_$PBS_JOBID

echo "#############################################################################" >> $HOME/finished
echo  `date` >> $HOME/finished
echo  `pwd` >> $HOME/finished
echo  "LASP" >> $HOME/finished
echo  "  " >> $HOME/finished