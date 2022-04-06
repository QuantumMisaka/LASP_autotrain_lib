#!/bin/bash

ROOTDIR=$(pwd)
source /home/james/app_envs/lasp_env.sh
PEXEC=/home/james/apps/LASP/lasp_NN
NSLOTS=16
# PYTHONEXEC=/home/james/apps/anaconda3/bin/python
export PYTHONPATH=/home/james/apps/anaconda3/bin:$PYTHONPATH
alias python=/home/james/apps/anaconda3/bin/python
MAXCYCLE=100
POTNAME=H2O

for i in {1..${MAXCYCLE}}
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
#  for jj in $( cat ../hostfile | sort -u )
#  do 
#    ssh $jj 'pkill -9 lasp'
#  dones
sleep 30
echo --- Start NN Training ---
echo `date`
Ntotal=$( grep Energy TrainStr.txt -c )
sed -i '/^NNtrain/d' lasp.in
sed -i '/^Ntrain/d' lasp.in
sed -i '1a \Ntrain '$Ntotal lasp.in
mpirun -rsh=ssh -np $NSLOTS $PEXEC > output

sleep 20
# input is pot-format tmp file?
sed -i 's/newrun/'"${POTNAME}"'.input/g' ${POTNAME}.pot
cp ${POTNAME}.pot ${POTNAME}.pot-$i
cp ${POTNAME}.pot ${POTNAME}.input
cd $ROOTDIR

printf " End SSW-DFT-NN Cycle %d\n" $i
echo `date`

done
