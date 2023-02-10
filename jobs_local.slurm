#!/bin/bash
#SBATCH -J auto-train
#SBATCH -p g1_user
#SBATCH --nodes=1
##SBATCH --time=72:00:00 
#SBATCH --ntasks-per-node=4
#SBATCH --error=%J.err          ###err log###
#SBATCH --output=%J.out         ###out log###

# setup intel mpi
source /es01/software/profile.d/intel2018.sh
# setup python
module load anaconda3

# setting
LASP_DIR=/es01/home/caoxiaoming/apps/lasp/INTER_pro3.3.4_intel12_2cw

if [ -d $I_MPI_ROOT/intel64 ]; then
  b="$I_MPI_ROOT/intel64"
elif [ -d $I_MPI_ROOT/lib ]; then
  b="$I_MPI_ROOT"
fi

if [ ! -f libmpigf.so.4 ]; then ln -s $b/lib/libmpifort.so.12 libmpigf.so.4 ; fi
if [ ! -f libmpi_dbg.so.4 ] ; then ln -s $b/lib/debug_mt/libmpi.so.12 libmpi_dbg.so.4 ; fi
export LD_LIBRARY_PATH=${LASP_DIR}:$LD_LIBRARY_PATH

unset I_MPI_PMI_LIBRARY
export I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=0

LIBVASP=${LASP_DIR}/Libraries/1.vasplib
export LD_LIBRARY_PATH=$LIBVASP:$LD_LIBRARY_PATH

EXEC=${LASP_DIR}/Src/lasp 

POTNAME="AuPdO"
COOR_PATTERNS=1 # for update coor-pattern from train_data every cycle
MAXCYCLE=100
NP=$SLURM_NPROCS
ROOTDIR=$(pwd)

SSWDIR=$ROOTDIR/SSW
NNDIR=$ROOTDIR/NN
VASPDIR=$ROOTDIR/VASP

cp $NNDIR/${POTNAME}.pot $SSWDIR/sourcedir/
cd $SSWDIR/sourcedir
for((i=1;i<=28;i++))
do
    # get 0-28 input file
    cp input.0 input.${i}
done
cd $ROOTDIR

for((i=1;i<=$MAXCYCLE;i++))
do
    # SSW-NN and VASP-DFT running
    printf " Start SSW-DFT-NN Cycle %d\n" $i 
    echo `date`

    if [ $i = "1" ]; then
        python auto.py > output
    else
        python auto.py >> output
    fi
    # VASP result backup
    if [ -d $ROOTDIR/VASP/cycle-0 ]; then
        mv $ROOTDIR/VASP/cycle-0 $ROOTDIR/VASP/para-0
    fi
    if [ -d $ROOTDIR/VASP/cycle-1 ]; then
        mv $ROOTDIR/VASP/cycle-1 $ROOTDIR/VASP/para-$i
    fi
    # NN training job
    cd ${NNDIR} # seems to be very important
    sleep 10
    echo --- Start NN Training ---
    echo `date`
    Ntotal=$(grep Energy TrainStr.txt -c)
    echo 'Num of Str/For Training data is:'
    echo $Ntotal
    sed -i '/^NNtrain/d' lasp.in
    sed -i '/^Ntrain/d' lasp.in
    sed -i '1a \Ntrain '$Ntotal lasp.in

    # running NN training in slurm
    # mpirun -np $NP $EXEC | tee run_print_out
    sbatch lasp.slurm
    # check if NNtrain is finished
    sleep 300
    NN_END=""
    while [[ -z $NN_END ]]
    do
        NN_END=$(grep Finished lasp.out)
        if [[ -z $NN_END ]]
	then
            sleep 20
        else
            echo "NN finished!"
        fi
    done
    sleep 20
    # input is pot-format tmp file?
    # should modify pot-name here
    sed -i 's/newrun/'"${POTNAME}"'.input/g' ${POTNAME}.pot
    cp ${POTNAME}.pot ${POTNAME}.pot-$i-${Ntotal}
    cp ${POTNAME}.pot ${POTNAME}.input
    mv lasp.out lasp.out-$i
    # pot and trainfor/trainstr must be backup
    cp TrainStr.txt TrainStr.txt-${i}-${Ntotal}
    cp TrainFor.txt TrainFor.txt-${i}-${Ntotal}
    # coor_patterns: not necessary
    # if (( COOR_PATTERNS ))
    # then
    #     python $ROOTDIR/get_coordination_patterns.py
    #     mv patterns_db.json ../traindata_patterns.json
    # fi
    cd $ROOTDIR

printf " End SSW-DFT-NN Cycle %d\n" $i
echo `date`

done

echo "#############################################################################" >> $HOME/finished
echo  `date` >> $HOME/finished
echo  `pwd` >> $HOME/finished
echo  "LASP" >> $HOME/finished
echo  "  " >> $HOME/finished
