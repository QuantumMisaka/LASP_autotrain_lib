#! /usr/bin/env  python
# noted by JamesBourbon: run on slurm-server
# many differences from original auto.py
# last change in 20220430, V1.0-testing
# should take care of Pool(process=poolsize)
# start by jobs_local.slurm, test pass
# prepara lasp_ssw.script in SSW, vasp.script in VASP, lasp.slurm in NN
# insert coordination patterns sampling method in nodejob_coor.py and update_patterns.py
# still need test

import os
import sys
import glob
import shutil
from multiprocessing import Pool
import time
from allstr_new import BadStr
from allstr_new import AllStr as AllStr_new # need to be noted
import numpy as np
from update_patterns import get_patterns
# for coor-patterns
SSW_choosing_mode = 1 # 0 for random 1 for coor-pattern
ROOTDIR = os.getcwd()
patterns_DB = f"{ROOTDIR}/traindata_patterns.json"

#  used in SSW-sampling collection
def nodejob_to_collect_allstr(rootdir, workdir, nbadstr):
    '''to parallelly collect allstr.arc result from workdirs, used by pool.apply_async()
    
    collect result will be print-out as outstr.arc in workdirs
    '''
    command = f'python {rootdir}/nodejob.py {workdir} {nbadstr}'
    os.system(command)
    
    
def nodejob_coordination_collect(rootdir, workdir, nbadstr):
    '''to parallelly collect allstr.arc result from workdirs, by coordination pattern choosing
    
    collect result will be print-out as outstr.arc in workdirs
    '''
    command =  f'python {rootdir}/nodejob_coor.py {workdir} {nbadstr}'
    os.system(command)


class RunSSW:
    def __init__(self,SSWdir,cpuperjob,prog,masternode,choose_mode=1):
        self.dir = SSWdir
        self.cpuperjob =int(cpuperjob)
        self.prog = prog
        self.nbadneed = 50
        self.masternode= masternode # masternode status setting
        self.poolsize =int(total_cpu/cpuperjob) # default
        self.choose_mode = choose_mode
        # poolsize parallel now only used for collect data
    
    
    def build_SSW_folder(self, SSWdir:str, njob:int, ncycle:int, allstr:int):
        '''generate new sturcture for SSW sampling in SSW folder
            randomly choose SSW_result and create folder with input.arc and input
        
        Args:
            SSWdir: SSW running dir
            njob(int): SSW para_running jobs number
            ncycle(int): SSW-DFT-NN cycle number
            allstr(0/1 bool): reading state of allstr?
            
        Return:
            workdir(list): chosen results to SSW-sampling
            Natom(list): atom numbers of each Str for SSW-sampling
        
        noted by JamesBourbon in 20220401
        '''
        os.chdir(SSWdir)
        # all thing done in SSWdir
        if os.path.exists('../NN/%s.pot'%(jobname)): # in LASP_Example jobname=H2O
            shutil.copy('../NN/%s.pot'%(jobname),'./sourcedir/') # copy trained pot to SSWfolder
        workdir = [] # SSW_rand_get file in workdir list
        Natom = [] # atom list for each SSW_rand_list work object

        AllStr0 = AllStr_new() # create AllStr object
        AllStr0.arcinit([0,0],'./sourcedir/allstr-ini.arc')
        AllStr0.random_arange(200) # random set allstr-ini.arc for 200 times
        for i in range(28): # 28 is a magical number, equal to ntype
            AllStr0.gen_arc([i],'outstr.arc') # get i_index Str from AllStr
            os.system('mv outstr.arc ./sourcedir/input.arc_%s'%(i)) # new str naming
        inputfile = glob.glob('./sourcedir/input.arc*') # return list of all sampling structure
        print(inputfile)
        inputlist = [] # new structure list by str id
        for item in inputfile:
            inputlist.append(item.split('_')[-1]) # id of new structure
        ntype = len(inputlist) # length of new structure list
        
        # have done: random get new str from allstr-ini.arc
        # something omitted

        for i in range(njob):
            # njob = SSWpara, indicate number for parallel_running jobs
            # typically njob= 2 or 4
            itype = i%ntype 
            dir ="SSW-%d-%d"%(ncycle,i) # mkdir contain SSW-rand result for VASP-DFT
            workdir.append(dir)
            if not glob.glob(dir): os.mkdir(dir) # if not have, mkdir
            # set input.arc (lasp.str)
            shutil.copy('./sourcedir/input.arc_%s'%(inputlist[itype]),'%s/input.arc'%dir)
            # set input file
            # modified by JamesBourbon in 20220401
            inputfile = "./sourcedir/input.%s"%(inputlist[itype])
            # use sed -i to add parameter
            # input file already in SSW/sourcedir
            # sed d-parameter to delete by key words
            os.system('sed -i "/^SSW.printevery/d" %s '%(inputfile))
            os.system('sed -i "/^SSW.printselect/d" %s '%(inputfile))
            os.system('sed -i "/^SSW.printdelay/d" %s '%(inputfile))
            os.system('sed -i "/^SSW.Safe_hardcurv/d" %s '%(inputfile))
            os.system('sed -i "/^supercell/d" %s '%(inputfile))
            if allstr==0:
                os.system('sed -i "1i\SSW.printevery F" %s'%(inputfile))
            if allstr==1:
                # add parameter to input file
                # cannot run this sed on macOS
                os.system('sed -i "1i\SSW.printevery T" %s'%(inputfile))
                os.system('sed -i "1i\SSW.printselect 6" %s'%(inputfile))
                os.system('sed -i "1i\SSW.printdelay  2" %s'%(inputfile))
                os.system('sed -i "1i\SSW.Safe_hardcurv  150" %s'%(inputfile))

            shutil.copy('%s'%(inputfile),'%s/input'%dir)
            shutil.copy('./sourcedir/%s.pot'%(jobname),dir)
            Natom.append(int(os.popen('cat ./%s/input.arc | wc -l'%dir).readline().strip())-7)
            # calc number of atoms in input.arc file, append to Natom list
            if Natom[i]<40: os.system('sed -i "1i\supercell 2 1 1" ./sourcedir/input.%s'%(inputlist[itype]))
            # intentially train some large structure, supercell for enlarge cell

        return workdir,Natom

    
    def run(self,njob:int ,ncycle:int, nbadneed:int, maxtime:float, allstr:int, checkcycle:int, script="lasp_ssw.script"):
        '''main SSW running func on local processing
        
        Args:
            njob: parallel running jobs number
            ncycle: SSW_DFT_NN cycle number
            nbadneed: number of badStr needed?
            maxtime: running maxtime
            allstr: str status?
            checkcycle: SSWcheckcycle in console
            
        Returns:
            nbad: structures for iter-NN-training
        
        noted by JamesBourbon in 20220403
        '''
        # self.poolsize = njob
        workdirs,n_atom_list = self.build_SSW_folder(self.dir,njob,ncycle,allstr)
        # get init allstr.arc folders for SSW-sampling by build_SSW_folder 
        # workdirs: project folders for SSW-sampling
        # n_atom_list: list of number of atoms for each SSW-sampling project 
        os.system('rm -rf killsignal')
        exit = False
        self.nbadneed=nbadneed
        alltime = 0
        # running start
        jobids = []
        for i in range(njob):
            # Ncore=int(math.ceil(n_atom_list[i]/(math.ceil(float(n_atom_list[i])/self.cpuperjob))))
            # adjust Ncore to fit the atom numbers of structure: why?
            shutil.copy(script, workdirs[i])
            os.chdir(workdirs[i])
            command = "sbatch -J %s %s"%(workdirs[i], script)
            # use os.popen to get jobid
            jobids.append(int(os.popen(command).readline().strip().split()[-1]))
            os.chdir(self.dir)

        # SSW terminate judgement
        while not exit:
            time.sleep(checkcycle)
            alltime = alltime + checkcycle
            # alltime to get running time?
            self.collect_data(workdirs,ncycle,allstr)
            # get allstr-arc-i which is already calc-ed
            # time.sleep(10)
            nBad = 0
            for _tmp in range(5):
            #    if glob.glob('allstr.arc-%d'%ncycle):
                try: 
                    nBad = int(os.popen('grep Energy allstr.arc-%d -c'%ncycle).readline().strip())
                    #if (allstr==0 and nBad > nbadneed) or (allstr==1 and nBad > nbadneed*5) : break
                    if (allstr==0 and nBad > nbadneed) or (allstr==1 and nBad > nbadneed) : break
                except: time.sleep(10)
#           nBad = int(os.popen('grep Energy allstr.arc-%d -c'%ncycle).readline().strip())
            #if (allstr==0 and nBad > nbadneed) or (allstr==1 and nBad > nbadneed*5) :
            if (allstr==0 and nBad > nbadneed) or (allstr==1 and nBad > nbadneed) :
                print('terminate: enough nBad')
                os.system('echo SSW terminate: enough nBad')
                self.send_stop_signal(workdirs,ncycle)
                # time.sleep(60)
                exit= True
            if alltime >= maxtime:
                print('teminate: too long time')
                os.system('echo SSW terminate: too long time')
                self.send_stop_signal(workdirs,ncycle)
                # time.sleep(60)
                exit =True
            if glob.glob('killsignal'):
                send_kill_signal(workdirs)
                print('External signal: will be terminated in few seconds')
                os.system('echo SSW terminate: external signal')
                # time.sleep(60)
                exit=True
                
        # terminate SSW
        for jobid in jobids:
            cancel_command = "scancel %d"%(jobid)
            os.system(cancel_command)

        if alltime >= maxtime:
            print('-------suggest change setting---------')
            os.system('echo "suggest change setting"')
        time.sleep(30)
        
        self.collect_data(workdirs,ncycle, allstr)
        # collect data from SSW_samling running result in SSW workdir
        # prepara for VASP-DFT
        time.sleep(10)
        for _tmp in range(5):
            try:
                nBad = int(os.popen('grep Energy allstr.arc-%d -c'%ncycle).readline().strip())
                #if (allstr==0 and nBad > nbadneed) or (allstr==1 and nBad > nbadneed*5) : break
                if (allstr==0 and nBad > nbadneed) or (allstr==1 and nBad > nbadneed) : break
            except: time.sleep(10)
        choose_noted = 'Num of Structure to VASP-DFT: %d '%(nBad) 
        print(choose_noted)
        os.system('echo %s'%(choose_noted))
        return nBad   
 
    def send_stop_signal(self,workdirs, ncycles):
        '''for terminate SSWjobs if time-out or enough nbad is got'''
        os.chdir(self.dir)
        for i in range(len(workdirs)):
            os.system('echo stop > %s/softstop'%(workdirs[i]))
        return 


    def collect_data(self,workdirs:list,ncycle:int,allstr:int, choose_mode = 1):
        '''for collect running result
            merge and move SSW-NN result to VASP workdir
            
        noting by JamesBourbon in 20220402
            
        '''
        os.chdir(self.dir)
        os.system('rm -f *.arc*') # remove last things

        # using process pool to collect data: not necessary
        if not os.path.exists('nodejob.py'):
            os.system ('ln -s ../nodejob.py .')
        pool = Pool(processes=8)
        for wdir in workdirs: 
            if allstr == 0: 
                os.system('cat %s/Badstr.arc >> allstr.arc-%d'%(wdir,ncycle))
            else:
                # switch random and coor-pattern choosing method
                if choose_mode == 1:
                    pool.apply_async(nodejob_coordination_collect, (self.dir, wdir, self.nbadneed))
                else:
                    pool.apply_async(nodejob_to_collect_allstr, (self.dir, wdir, self.nbadneed))
        time.sleep(60)
        
        for i in range(len(workdirs)):
            os.system('cat %s/Badstr.arc >> Badstr.arc-%d'%(workdirs[i],ncycle))
            # collect badstr
            if glob.glob('%s/outstr.arc'%(workdirs[i])): 
                os.system('cat %s/outstr.arc >> allstr.arc-%d'%(workdirs[i],ncycle))
            # get allstr collect result - from nodejob.py
        os.system('rm -f ../VASP/*.arc-%d; ln -s ../SSW/allstr.arc-%d ../VASP/'%(ncycle,ncycle))
        os.system('ln -s ../SSW/Badstr.arc-%d ../VASP/'%(ncycle))
        pool.close()
        pool.join()
        return

class RunVASP:
    def __init__(self,vaspdir,cpuperjob:int,prog, masternode:int, limit=200):
        self.dir =vaspdir
        self.cpuperjob =int(cpuperjob)
        self.prog = prog
#       self.base =base
        self.masternode = masternode
        self.poolsize = int(total_cpu/cpuperjob)
        self.NELM_limit = limit


    def run(self,ncycle:int,nbadneed:int,maxF:float,maxtime:float,script="vasp.script"):
        '''main VASP-run methods'''
        self.pooluse =0
        self.nmax = nbadneed
        self.maxF = maxF

        os.chdir(self.dir)
        self.screen_data('allstr.arc-%d'%ncycle, nmax= nbadneed)
        if not glob.glob('cycle-%d'%ncycle): 
            os.mkdir('cycle-%d'%ncycle)
        if glob.glob('outstr.arc'): 
            os.system('mv outstr.arc ./cycle-%d/allstr.arc-%d'%(ncycle,ncycle))

        if glob.glob('Badstr.arc-%d'%ncycle):
            self.screen_data('Badstr.arc-%d'%ncycle, nmax= nbadneed)
            if glob.glob('outstr.arc'): 
                os.system('cat outstr.arc >> ./cycle-%d/allstr.arc-%d'%(ncycle,ncycle))

        if glob.glob('AddVASPcal.arc'):
            os.system('cat AddVASPcal.arc >> ./cycle-%d/allstr.arc-%d'%(ncycle,ncycle))
            os.system('mv AddVASPcal.arc AddVASPcal.arc-addedcycle%d'%(ncycle))
        # all above is for generate structures for VASP_DFT
        # from allstr.arc, Badstr.arc, AddVASPcal.arc 
        
#       os.system('rm -f allstr.arc-%d'%ncycle)

        workdir = os.path.join(self.dir,'cycle-%d'%ncycle)
        # VASP-DFT work dir: VASP/cycle-i
        
        #shutil.copy('allstr.arc-%d'%ncycle,workdir)
        os.chdir(workdir)
        os.system('rm -rf killsignal')

        # running setting
        AllStr= AllStr_new()
        AllStr.readfile('allstr.arc-%d'%ncycle)

        workdirs=[]
        exit =False
        alltime =0
        job_nums = len(AllStr)
        jobids = []

        for i in range(job_nums):
            # running in cycle-i/para-i dir
            AllStr.print_list([i],'input.arc_%d'%i)
            if not glob.glob('para%d'%(i+1)): 
                os.system('mkdir para%d'%(i+1))
            if glob.glob('para%d/allstr.arc'%(i+1)) : 
                # if VASP-DFT is calc-ed and good, ignore it 
                icontrol= int(os.popen('cat para%d/OSZICAR | wc -l'%(i+1)).readline().strip())
                # icontrol0= int(os.popen('cat para%d/allstr.arc | wc -l'%(i+1)).readline().strip())
                if icontrol < self.NELM_limit: continue
           
            # for VASP-DFT we need POSCAR(not necc in lasp-vasp) and POTCAR
            # should put atom POTCAR in VASP/sourcedir
            AllStr[i].outPOSCAR('POSCAR_%d'%i)
            AllStr[i].genPOTCAR('../sourcedir/','POTCAR_%d'%i)

            # self.prog=para['prog']['VASP']
            # for molecular, automatically use VASPgamma
            # others use pre-KPOINTS
            # get KPOINTS (from prepared sourcedir) 
            if AllStr[i].abc[0]== AllStr[i].abc[1] and \
               AllStr[i].abc[1]==AllStr[i].abc[2] and \
               AllStr[i].abc[1] > 9.99 :
                os.system('\cp  ../sourcedir/KPOINTS_gamma  para%d/KPOINTS'%(i+1))
                self.prog=para['prog']['VASPgamma']
            else:
                AllStr[i].genKPOINTS('KPOINTS')
                os.system('mv KPOINTS para%d/'%(i+1))
            # I write a Allstr.genkpoints -- JamesBourbon
            # but still can directly copy KPOINTS: use mesh auto 25
            # os.system('cp  ../sourcedir/KPOINTS  para%d'%(i+1))
            os.system('rm -f para%d/input; cp  ../sourcedir/input  para%d'%(i+1,i+1))
            os.system('cp  ../sourcedir/INCAR  para%d'%(i+1))
            #os.system('sed -i "/^LDAU /d" para%d/INCAR'%(i+1))
            os.system('mv input.arc_%d para%d/input.arc'%(i,i+1))
            os.system('mv POSCAR_%d para%d/POSCAR'%(i,i+1))
            os.system('mv POTCAR_%d para%d/POTCAR'%(i,i+1))

# something omitted
            rundir = '%s/para%d'%(workdir, i+1)
            workdirs.append(rundir)
            # print('------runjob%d-----'%(i+1))
            os.system("cp ../%s ."%(script))
            # running vasp job in VASP/cycle-1
            shutil.copy(script, rundir)
            os.chdir(rundir)
            # in slurm
            command = "sbatch -J dft-%d %s"%(i+1, script)
            jobids.append(int(os.popen(command).readline().strip().split()[-1]))
            # need to back to workdir
            os.chdir(workdir)
            # running cycle end

        time.sleep(300)    

        # stop or not judgement
        while not exit:
            time.sleep(60)
            alltime= alltime + 60
            os.chdir(workdir)
            if alltime >= maxtime:
                print('teminate: too long time')
                send_kill_signal(workdirs)
                time.sleep(60)
                exit = True
            if glob.glob('killsignal'):
                send_kill_signal(workdirs)
                print('External signal: will be terminqated in few seconds')
                time.sleep(60)
                exit=True
            # judge vasp-calc finished or not
            endjobs = int(os.popen("grep E0 para*/OSZICAR | wc -l").readline().strip())
            print("end jobs: %d"%endjobs)
            if endjobs >= job_nums:
                print("VASP-calc. totally done")
                time.sleep(60)
                exit = True
                
        # just in case
        for jobid in jobids:
            print("cancel %d"%(jobid))
            cancel_command = "scancel %d"%(jobid)
            os.system(cancel_command)


        # collect data
        print('---------Start collect data---------')
        self.collect_data(ncycle,len(AllStr)) # give labeled TrainStr and TrainFor
        #if(Lallstr==0): self.compare(ncycle)
        self.compare(ncycle)
        # self.screen_data('1.arc',forcefile ='2.arc')
        naddstr = int(os.popen("grep Energy TrainStr.txt -c").readline().strip())
        # update Coordination Patterns Database if set it
        # using get_patterns in update_patterns.py to update patterns
        if SSW_choosing_mode == 1:
            print("-------Start Update Coordination Patterns Database--------")
            get_patterns(init_db=patterns_DB) # give patterns_db.json
            os.system(f"cp patterns_db.json {patterns_DB}")  
        # add TrainStr and TrainFor to NNdir
        # the Structure before will also include in NNtrain
        os.system('cat TrainStr.txt >> ../../NN/TrainStr.txt  ')
        os.system('cat TrainFor.txt >> ../../NN/TrainFor.txt  ')
        print ('---VASP run done---')
        time.sleep(10)
        return naddstr


    def compare(self,ncycle):
        '''compare NN_single and VASP-DFT result to get RMSE
            aimed to judge NN_train can finish or not 
        
        noted by JamesBourbon in 20220317
        '''
        # before VASP
        nn = AllStr_new()
        nn.readfile('allstr.arc-%d'%ncycle)
        # After VASP
        vasp = AllStr_new()
        vasp.readfile('allstr.arc')
        # compare energy
        if (len(nn)!= len(vasp)):
            print('some str failed to dft cal')
        else :
            err_all = 0
            for i in range(len(nn)):
                err_all=err_all+np.square(nn[i].energy-vasp[i].energy)
            rmse= np.sqrt(err_all/len(nn))
            print('----------- rmsE %14.8f eV ---------'%rmse)


    def collect_data(self,ncycle:int,nstr:int):
        '''collect VASP-DFT data gen-ed by vasp.6.2
        
        Args:
            ncycle; cycle number
            nstr: structure number
        '''
        os.chdir(self.dir)
        os.chdir('cycle-%d'%ncycle)
        os.system('rm -rf allstr.arc')
        os.system('rm -rf TrainStr.txt')
        os.system('rm -rf TrainFor.txt')
        # os.system('rm -rf 2.arc')
        for i in range(nstr):
            if glob.glob('para%d/OSZICAR'%(i+1)):
                icontrol= int(os.popen('cat para%d/OSZICAR | wc -l'%(i+1)).readline().strip())
                if(icontrol < self.NELM_limit): # if DFT-single end properly
                    os.chdir("para%d"%(i+1))
                    os.system("python ../../../vasp2lasptrain.py")
                    os.system("python ../../../shiftformat.py train2arc")
                    os.chdir("%s/cycle-%d"%(self.dir, ncycle))
                    os.system('cat para%d/allstr.arc >> allstr.arc'%(i+1))
                    os.system('cat para%d/TrainStr.txt >> TrainStr.txt'%(i+1))
                    os.system('cat para%d/TrainFor.txt >> TrainFor.txt'%(i+1))
                    # os.system('cat para%d/allfor.arc >> 2.arc'%(i+1))
        return

    def screen_data(self,strfile,forcefile = False,nmax =999999):
        '''for screen structure in strfile (but for what?) and force if exist
            key function for VASP-run str_file setting
            
        if finished, outstr.arc will be print-out (if force input, outfor.arc be print-out)
            
        noting by JamesBourbon in 20220402
        
        '''        
        AllStr = AllStr_new()
        if forcefile:
            AllStr.arcinit([0,0],strfile,forcefile)
        else :
            AllStr.arcinit([0,0],strfile)
        # Here can set HighE,MaxAngle,MinAngle
        if len(AllStr)==0: return
        # for screen structure/force data
        b=BadStr()
        #   b.HighE=-3.0
        b.MaxFor = self.maxF
        b.MaxLat = 40
        b.MinLat = 2.2
        #
        AllStr = AllStr.filter(b)

        # main filter function use allstr_new.py and structure_new.py
        if(len(AllStr) > nmax+50):
            AllStr.random_arange(200)
            AllStr = AllStr_new(AllStr[:(nmax+50)])

        print('All Str:',len(AllStr))
        #print 'present force',AllStr[0].Lfor
    
        if len(AllStr)==0: return
        if(len(AllStr) > nmax):
            AllStr.random_arange(200)
            AllStr = AllStr_new(AllStr[:(nmax)])
    
        #   AllStr.sort_by_energy()
    
        print('Final Dump Str:',len(AllStr))
        if len(AllStr) >0:
            AllStr.gen_arc(range(len(AllStr)),'outstr.arc',2)
            if AllStr[0].Lfor: 
                AllStr.gen_forarc(range(len(AllStr)),'outfor.arc',2)



class RunNNTraining:
    def __init__(self,NNdir,cpuperjob,prog):
        self.dir = NNdir
        self.cpuperjob = int(cpuperjob)
        self.prog = prog
        #self.poolsize = 2 # bailan coding
        #self.prog='/home7/hsd/program/traincnn-1.0/trainCNN-1.0'
    
    def run(self,ncycle,nadd,NNstd = False, NNepoch = 2000):
        self.set_file(ncycle, nadd, NNepoch)
        os.system('rm -rf killsignal')

#       while not exit:
#           time.sleep(60)
#           alltime= alltime +60
#           #if alltime >= maxtime:
#           #    print 'teminate: too long time'
#           #    self.sendstopsignal(workdirs,ncycle)
#           #    time.sleep(60)
#           #    exit =True
#           #    #pool.terminate()
#           if glob.glob('killsignal'):
#               sendkillsignal(['.'])
#               print 'External signal: will be terminated in few sceonds'
#               time.sleep(60)
#               exit=True
#           process_count=len(pool._cache)
#           if process_count ==0:
#               print 'all work done'
#               exit=True
#           if NNstd:
#               try:
#                   MAXinfo=os.popen('grep MAX hTrainOutput | tail -1').readline().split()
#                   maxE = float(MAXinfo[3])
#                   maxF = float(MAXinfo[5])
#                   RMSinfo=os.popen('grep RMS hTrainOutput | tail -1').readline().split()
#                   rmsE = float(RMSinfo[3])
#                   rmsF = float(RMSinfo[5])
#
#                   if ((maxE < NNstd['maxE']) and (maxF < NNstd['maxF']) and \
#                       (rmsE < NNstd['rmsE']) and (rmsF < NNstd['rmsF'])):
#                       sendkillsignal(['.'])
#                       print 'terminate: has reached NN std'
#                       time.sleep(60)
#                       exit=True
#               except:
#                   continue
#           #if glob.glob('terminatesignal'):
#           #    print 'forced termination of pool'
#           #    exit=True
#           #    pool.terminate()
#           #    time.sleep(60)
#
#
#
#       process_count=len(pool._cache)
#       if (exit) and (process_count!=0):
#           time.sleep(60)
#           process_count=len(pool._cache)
#           if process_count!=0:
#               sendkillsignal(['.'])
#               time.sleep(60)
#               print 'forced termination of pool'
#               pool.terminate()
#
#
#
#       pool.close()
#       pool.join()


    def set_file(self, ncycle, nadd, NNepoch):
        '''setting NNtrain calc in rootdir/NN
        
        NN_training will be running is NNdir: rootdir/NN 
        '''
        os.chdir(self.dir)
        os.mkdir('cycle-%d'%ncycle)
        os.system('cp TrainFor.txt TrainStr.txt cycle-%d'%(ncycle))
        os.system('cp %s.pot %s.input'%(jobname,jobname))
        os.system('cp %s.pot cycle-%d'%(jobname,ncycle))
        os.system('mv lasp.out cycle-%d'%(ncycle))
        os.system('mv SavePot/ cycle-%d'%(ncycle))
        os.system("sed -i 's/NNepochs.*/NNepochs %d/g' lasp.in"%(NNepoch))


        self.nstr = self.nstr +nadd
        #os.system("sed -i 's/Ntrain.*/Ntrain  %d/g' lasp.in"%self.nstr)

    def get_job_info(self):
        '''get job information including jobname and autobase'''
        os.chdir(self.dir)
        #line = os.popen('grep Jobname NNinput').readline().split()
        #self.jobname = line[-1]
        f = open('lasp.in','r')
        line=f.readline()
        while line:
            if len(line.split()) == 0:
                line=f.readline()
                continue
            elif line.split()[0]=='Jobname':
                self.jobname = line.split()[1]
            elif line.split()[0]=='Ntrain':
                self.nstr= int(line.split()[1])
            elif line.split()[0] =='%block':
                blockname = line.split()[-1]
                autobase={}
                blockline=f.readline().split()
                while blockline[0] != '%endblock':
                    autobase[blockline[0]]= 0.0
                    blockline=f.readline().split()
                print(autobase)
            line=f.readline()

        f.close()
        return self.jobname, autobase


def pool_shell(poolsize):
    pool = Pool(processes=poolsize)
    # create precess pool for multiprocessing python
    return pool


def read_block(lines):
    '''read block_lines format
    
    noting by JamesBourbon
    '''
    dict = {}
    for line in lines:
        try:
            # split number and string
            dict[line[0]]= float(line[1])
            if dict[line[0]].is_integer():
                dict[line[0]]=int(dict[line[0]])
        except:
            dict[line[0]]= line[1]
    return dict

def read_para():
    '''compute parameter setting
        read info from console file
        
    noted by JamesBourbon
    '''
    f= open('console','r')
    # variable para definition and its default setting
    para ={'Nbad':200,'maxF':100,'maxSSWtime':99999,'maxVASPtime':99999,\
            'maxtimeperVASP':99999,'maxcycle':99999,'masternode':0}
    # reading info from console
    line=f.readline()
    while line:
    # for line in f:
        # read empty
        if len(line.split()) == 0:
            line=f.readline()
            continue
        # read Nband, maxF, maxSSWtime, cyclecontrol
        # read maxVASPtime, maxtimeperVASP, SSWprog, VASPprog, NNprog
        # read SSWcpu, VASPcpu, (NN not use)
        # read masternode, maxSSWpara, Allstr, cpupernode
        # read NNepoch, StartfromVASP, SSWcheckcycle
        elif line.split()[0]=='Nbad':
            para['Nbad'] = int(line.split()[1])
        elif line.split()[0]=='maxF':
            para['maxF'] = int(line.split()[1])
        elif line.split()[0]=='maxSSWtime':
            para['maxSSWtime']= int(line.split()[1])
        elif line.split()[0] =='cyclecontrol':
            para['maxcycle']= int(line.split()[1])
        elif line.split()[0] == 'maxVASPtime':
            para['maxVASPtime']=int(line.split()[1])
        elif line.split()[0] == 'maxtimeperVASP':
            para['maxtimeperVASP']=int(line.split()[1])
        elif line.split()[0] == 'SSWprog':
            para['SSWprog']= line.split()[1]
        elif line.split()[0] == 'VASPprog':
            para['VASPprog']= line.split()[1]
        elif line.split()[0] == 'NNprog':
            para['NNprog']= line.split()[1]
        elif line.split()[0] == 'SSWcpu':
            para['SSWcpu'] =int(line.split()[1])
        elif line.split()[0] == 'VASPcpu':
            para['VASPcpu'] =int(line.split()[1])
        elif line.split()[0] =='masternode':
            para['masternode'] = int(line.split()[1])
        elif line.split()[0] =='maxSSWpara':
            para['maxSSWpara'] = int(line.split()[1])
        elif line.split()[0] =='Allstr':
            para['Allstr'] = int(line.split()[1])
        elif line.split()[0] =='cpupernode':
            para['cpupernode'] = int(line.split()[1])
        elif line.split()[0] == 'NNepoch' :
            para['NNepoch'] = int(line.split()[1])
        elif line.split()[0] == 'StartfromVASP' :
            para['StartfromVASP'] = int(line.split()[1])
        elif line.split()[0] == 'SSWcheckcycle' :
            para['SSWcheckcycle'] = int(line.split()[1])
        # block read
        elif line.split()[0] =='%block':
            blockname = line.split()[-1]
            lines = []
            blockline=f.readline().split()
            while blockline[0] != '%endblock':
                lines.append(blockline)
                blockline=f.readline().split()
            para[blockname]=read_block(lines)
        line=f.readline()
    f.close()
    return para



def send_kill_signal(workdirs):
    for i in range(len(workdirs)):
        os.system('echo kill > %s/killsignal'%(workdirs[i]))
    return


if __name__ == "__main__":
    newstart = 1

    global para
    para = read_para() # read parameter from console
    
    rootdir = ROOTDIR
    # check SSW VASP NN file
    if not glob.glob('SSW/sourcedir'):
        print('not find SSW file')
        sys.exit()
    if not glob.glob('VASP/sourcedir'):
        print('not find vasp file')
        sys.exit()
    if not glob.glob('NN'):
        print('not find training file')
        sys.exit()


    global total_cpu # use it to get poolsize for Pool
    total_cpu = para['cpupernode']
    # directory definition
    SSWdir = os.path.join(rootdir,'SSW')
    VASPdir = os.path.join(rootdir,'VASP')
    NNdir = os.path.join(rootdir,'NN')

    # clean file 
    # suggest building the new folder for new job
    os.chdir(rootdir)
    if(newstart == 1):
        os.system('rm -rf SSW/SSW*')
        # os.system('rm -rf VASP/cycle-*')
        os.system('rm -rf NN/cycle-*')

    # setting SSW-NN
    SSW = RunSSW(SSWdir,para['cpuperjob']['SSW'],para['prog']['SSW'],para['masternode'], SSW_choosing_mode)
    # totalproc =SSW.set_hostfile()
    # setting VASP-DFT
    VASP = RunVASP(VASPdir,para['cpuperjob']['VASP'],para['prog']['VASP'],para['masternode'])
    # VASP.set_hostfile()
    # setting NNtraining
    NN = RunNNTraining(NNdir,para['cpuperjob']['NN'],para['prog']['NN'])
    # NN.set_hostfile()

    global jobname
    global autobase
    jobname,autobase = NN.get_job_info()

    # global poolcount
    # poolcount =0d



    SSWpara = int(min(para['maxSSWpara'],total_cpu/para['cpuperjob']['SSW']))
    # define njob-SSW
    SSW_note = '==========================start SSW sampling============================='
    VASP_note = '==========================start VASP calculation========================='
    NN_note = '===========================start NN training============================='

# for add some str
# if StartfromVASP
    if (para['StartfromVASP']== 1):
        print('add some str')
        print(VASP_note)
        os.system(f' echo {VASP_note}')
        para['Nbad'] = 30000 # ?
        #nadd= VASP.run(0,para['Nbad'],para['maxF'],para['maxVASPtime'],para['maxtimeperVASP'])
        nadd = VASP.run(0,99999,para['maxF'],para['maxVASPtime'],para['maxtimeperVASP'])
        # poolcount = poolcount + VASP.poolsize
        print('===========================start NN training=============================')
        NN.run(0,nadd, para['NNstd'],para['NNepoch'] )
        # poolcount=poolcount+ NN.poolsize
        #because now NN is running from jobs.sh
        os.chdir(rootdir)
        os.system("sed -i 's/StartfromVASP.*/StartfromVASP 0/g' console")
        sys.exit()

# if StartfromNN: directly do
    startcycle = 1 # question remained

    for i in range(1):
        icycle = startcycle + i

        os.chdir(rootdir)
        # para = read_para() # seemed to be replica
        if icycle > para['maxcycle']: 
            break
        print (' ')
        print ('-------------------------->   Start cycle   <--------------------------')
        print (' ')
        print(para)
        print (' ')
        print(SSW_note)
        os.system(f'echo {SSW_note}')
        nbadcol = SSW.run(SSWpara,icycle,para['Nbad'],para['maxSSWtime'], para['Allstr'],para['SSWcheckcycle'])
        # poolcount = poolcount + SSW.poolsize
        #nbadcol = 1
        if(nbadcol > 0) :
            print(VASP_note)
            os.system(f' echo {VASP_note}')
            nadd= VASP.run(icycle,para['Nbad'],para['maxF'],para['maxVASPtime'])
            # poolcount=poolcount+ VASP.poolsize
            print(NN_note)
            os.system(f'echo {NN_note}')
#           os.system('sed -i "/^NNepochs " NN/NNinput' )
#           os.system('sed -i "2i\NNepochs 2000" NN/NNinput' )
            NN.run(icycle, nadd, para['NNstd'],para['NNepoch'])
            # NN.run only for NNfile preparing
            # poolcount = poolcount + NN.poolsize
        time.sleep(5)
    print('NNfile prepared')
    # finished reading in 20220403
    # print('auto-train end')
