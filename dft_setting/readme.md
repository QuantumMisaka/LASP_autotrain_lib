# DFT-SETTING DIR FOR VASP-LABEL JOBS

## Which File Need?

INCAR: for single-dft vasp-calc. 

! NOTICE: KPOINTS will automatically get for automash-25 (in Str.genKPOINTS), but if you want to do single-dft to cluster system, please set KPOINTS to gamma

lasp.in: input setting file for lasp-inter vasp

POTCAR.[element]: POTCAR of this element, program will cat these POTCAR automatically

radius.csv: radius dataset for coordination pattern calc. (NOT NECESSARY)

## Function

DFT setting dir for doing vasp single-dft, like rootdir/VASP/sourcedir but can have some difference in INCAR or POTCAR.[element]

Directly used by coor_verlet_sample.py