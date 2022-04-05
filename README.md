# LASP SSW-DFT-NN auto-train python-lib

## Author
Original: ZPLiu's Group (SDHuang, ZPLiu et al)
Modified: James.Misaka.Bourbon.Liu
Last Change: 2022-04-05
Version: V1.0

## How to run it directly

### 1. modify console file

#### 1.1 some following are parameters often need to change

StartfromVASP 0    #  0 start from SSW sampling provided with NN pot 
                   #  1 start from allstr.arc-0 in VASP dir, which is often used to train the first NN pot

Nbad   40    # structures for VASP every cycle

cpupernode 96   # CPU total core (not suggested running between 2 or more CPUs)
SSWcheckcycle  600   # SSW time clock 600 seconds
%block cpuperjob
SSW  24         # cores per SSW job
VASP 24         # cores per VASP job
NN   0          # should be designated in jobs.sh
%endblock cpuperjob

#### 1.2 provide the binary program

%block prog
SSW  /home10/bin/lasp-1.0-release/lasp
VASP  /home10/bin/lasp-1.0-release/lasp
VASPgamma  /home10/bin/lasp-1.0-release/lasp.gamma
NN  /home10/bin/lasp-1.0-release/lasp
%endblock prog

#### 1.3 provide the element name in console file

%block base
O   0.0
H   0.0
%endblock base


### 2. modify jobs.sh

a. make sure the name of NN pot is correct
   e.g. sed -i 's/H2O/PtOH.pot/g' jobs.sh
b. modify the number of cycles, default is 100
...
for i in {1..100}
...

c. modify the cpu/cores required for your computing cluster
(modify it in jobs.sh)

### 3. creat arc file for SSW sampling:  allstr-ini.arc
    SSW/sourcedir/allstr-ini.arc
    you may get allstr-ini.arc from the examples of structure
    which you need to add and train in your pot


### 4. check NN directory
   In rootdir/NN you shouldp prepare:
   lasp.in
   H2O.pot          # not required if start from scratch
   H2O.input        # if start from scratch, use "newrun" for pot
   TrainStr.txt     # if start from scratch, just creat an empty file
   TrainFor.txt     # if start from scratch, just creat an empty file
   adjust_factor

!  lasp already has a lot of Train*.txt files for different systems
!  please first download TrainStr.txt TrainFor.txt from www.lasphub.com

### 5. make sure add python exec path in jobs.sh
alias python=/data/apps/intel/intelpython3/bin/python
or 
export PYTHONPATH=/data/apps/intel/intelpython3/bin:$PYTHONPATH


### 6. qsub jobs.sh


### 7. some tips

1. SSW/sourcedir have input.i for SSW-NN input file, may need to check
2. remember check SSW VASP NN dir before finally running
3. auto.py and SSW-DFT-NN auto still need to be tested

