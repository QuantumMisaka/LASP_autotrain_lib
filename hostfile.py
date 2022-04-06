# python exec related to run auto.py in cluster by host
# noted by JamesBourbon in 20220316; not finished
# remained questions
import os
import time
import multiprocessing
import glob
import traceback
import subprocess

class HostFile(object):
    def __init__(self, workdir, cpuperjob, setmasternode = 0):
        self.workdir = workdir
        self.cpuperjob = cpuperjob
        self.setmasternode = setmasternode

    def getProc(self):
        # get hostdict
        # should be specified by different host?
        # very big question: how to get hostfile and environment variable
        hostdict = {}
        
        if self.setmasternode != 0:
            first = True
            print('set masternode')
        else:
            first = False
        
        # get environment variable of host file list in group cluster
        # for SGE in zpliu group
        '''
        hostfile = os.environ["PE_HOSTFILE"]
        print(hostfile)
        fp = open(hostfile, "r")
        for x in fp:
            line = x.split()
            hostdict[line[0]] = int(line[1])
        fp.close()
        '''
        
        #the following part are designed for fudannhpcc system
        # for PBS?
        hostfile = os.environ["PBS_NODEFILE"]
        print(hostfile)
        fp = open(hostfile,"r")
        for x in fp:
            if first:
                masternode = x.split()[0]
                first = False
                continue
            line = x.split()
            hostdict[line[0]] = 12 #int(line[1]) zpliu
        fp.close()
        if self.setmasternode: hostdict.pop(masternode)
        return hostdict
    
    def alloProc(self,hostdict,size=1):
        # construct a list of possible hosts
        totProc = 0
        hostList = []
        for key,val in hostdict.items():
            totProc += val
            hostList.extend([key]*val)
        print("Availiable proc number: ",totProc)
        if size ==0:
            size = totProc
        return self.chunks(hostList,size),totProc

    
    def set_hostfile(self):
        hostDict = self.getProc()

        divHost,totalproc = self.alloProc(hostDict,self.cpuperjob)
        print(divHost)
        # get the pool size
        poolsize= len(divHost)
        print(poolsize)
        # dump the hostfiles
        self.hostInfo = self.dumpHost(divHost)
        return self.hostInfo,poolsize,totalproc

    # chuck the list
    def chunks(self,arr, n):
        all = len(arr)-len(arr)%n
        return [tuple(arr[i:i+n]) for i in range(0, all, n)]


    # dump hosts info into files
    # return a list of file names and corresponding host number

    def dumpHost(self,divHost):
        os.chdir(self.workdir)
        hostInfo = []
        ft = open(".hostfile", "w")
        for i,record in enumerate(divHost):
            fn = ".hostfile_%03d" %i
            fp = open(fn,"w")
            for line in record:
                fp.write(line+"\n")
                ft.write(line+"\n")
            fp.close()
            # why plus one ? in fact, I don't know either
            hostInfo.append((fn, len(set(record))))
        ft.close()
        return hostInfo


def runprog_local(workdir, prog, ncpus):
    '''main function for running LASP program on localhost
        main running process: mpirun LASP
    
    Args: 
        workdir(string): project dir
        prog(string): running program
        ncpus(int): number of cpu to use
    
    noting by JamesBourbon in 20220402
    '''
    try:
        os.chdir(workdir)
        mpiprog = "mpirun -np %d "%(ncpus) + prog
        fout = open('output','w')
        subprocess.call(mpiprog,stdout = fout, stderr=fout,shell=True,executable='/bin/bash')
        fout.close()
        return
    except Exception as err:
        # in py2: except Exception, e
        traceback.print_exc()
        raise err
    except:
        print('error')
        return 


def runprog_cluster(workdir,prog,ncpus,hostInfo,rootdir,env,poolcount = 0):

    try:
        os.chdir(rootdir)
        cwd = os.getcwd()
        exit = False
        # print(multiprocessing.current_process().name)
        nodeInfo = hostInfo[int(multiprocessing.current_process().name.split("-")[-1])-1-poolcount]
        os.chdir(workdir)
        mf = os.path.join(cwd,nodeInfo[0])
        #mpiprog = "/home/software/mpi/intel/impi/4.0.1.007/bin64/mpirun --rsh=ssh -machinefile %s -np %d "%(mf, ncpus) + prog
        mpiprog = "mpirun -r ssh -machinefile %s -np %d "%(mf, ncpus) + prog
        fout1 = open("proginfo","w")
        fout1.write("Current process: "+ multiprocessing.current_process().name +"\n")
        fout1.write("Host file: " + mf + "\n")
        fout = open('output','w')
        #totalrun = 'source ~/.bashrc; '+ mpiprog
        totalrun = mpiprog
        #nodename = os.popen('head -1 %s'%mf).readline().strip()
        #totalrun = 'source ~/.bashrc; ssh '+nodename +"; "+mpiprog
        fout1.write(totalrun+"\n")
        fout1.close()
        #child= subprocess.Popen(mpiprog,stdout = fout, stderr=fout,shell=True,executable='/bin/bash',preexec_fn = os.setpgrp)
        subprocess.call(mpiprog,stdout = fout, stderr=fout,shell=True,executable='/bin/bash')
        fout.close()
        return
    # except Exception,e:
    except Exception as err:
        traceback.print_exc()
        raise err
    except:
        print('error')
        return


def runprog_cluster_manual(workdir,prog,hostInfo,ncpus,rootdir,env,poolcount =0,maxtime= False):
    '''main function for running LASP program on cluster
        main running process: mpirun LASP
    
    Args: 
        workdir(string): project dir
        prog(string): running program
        hostInfo
        ncpus(int): number of cpu to use
        rootdir(string): program running rootdir
        env(string): not used
        poolcount(int): process pool using count?
    
    noting by JamesBourbon in 20220402
    '''
    try:
        os.chdir(rootdir)
        running_root = os.getcwd()
        exit = False # exit status
        print(multiprocessing.current_process().name)
        nodeInfo = hostInfo[int(multiprocessing.current_process().name.split("-")[-1])-1-poolcount]
        os.chdir(workdir)
        node_path = os.path.join(running_root,nodeInfo[0])
        #mpiprog = "/opt/intel/impi/5.0.2.044/intel64/bin/mpirun -machinefile %s -np %d "%(node_path, ncpus) + prog
        #mpiprog = "mpirun -machinefile %s -env I_MPI_DEVICE rdma:OpenIB-cma -np %d "%(node_path, ncpus) + prog
#        mpiprog = "source /home2/shang/.bashrc; mpirun -machinefile %s -np %d "%(node_path, ncpus) + prog
        
        mpiprog = " mpirun  -rsh=ssh -machinefile %s -np %d "%(node_path, ncpus) + prog # in ZPLiu Group
        # main running code including mpirun
        # key parameter: -machinefile, linked to nodeInfo and hostInfo
        # -machinefile for node indicating
        # mpiprog = " mpirun  -rsh=ssh -np %d "%(ncpus) + prog # in CCC-ECUST 138
#       
        fout1 = open("proginfo","w")
        fout1.write("Current process: "+ multiprocessing.current_process().name +"\n")
        fout1.write("Host file: " + node_path + "\n")
        fout = open('output','w')
        
        totalrun = mpiprog
        #nodename = os.popen('head -1 %s'%mf).readline().strip()
        #totalrun = 'source ~/.bashrc; ssh '+nodename +"; "+mpiprog
        fout1.write(totalrun+"\n")
        child= subprocess.Popen(mpiprog,stdout = fout, stderr=fout,shell=True,executable='/bin/bash',preexec_fn = os.setpgrp)
        # create a process link to mpiprog to run 'mpirun LASP'
        print('start run job in %s'%workdir)
        fout1.write('pid   %d\n'%child.pid)
        fout1.close()
        pid = child.pid
        if maxtime:
            alltime = 0
        while not exit:
            time.sleep(30)
            returnCode= child.poll()
            if glob.glob('killsignal'):
                os.kill(-pid,9)
                time.sleep(3)
                fout.write('kill %s\n'%pid)
                a=os.waitpid(pid,0)
                print(a)
                time.sleep(3)
                exit = True
            if isinstance(returnCode,int):
                # judge running successfully or not
                if returnCode == 0:
                    fout.write('successfully done\n')
                else:
                    fout.write('something wrong: returnCode  %d\n'%returnCode)
                exit = True
            if maxtime:
                alltime = alltime+30
                if alltime >maxtime:
                    os.kill(-pid,9)
                    os.system('pkill -9 %s'%prog)
                    time.sleep(3)
                    fout.write('time out :kill %s\n'%pid)
                    a=os.waitpid(pid,0)
                    print(a)
                    time.sleep(3)
                    exit =True

        #fout1.write(str(child.poll())+'\n')
        fout.write('exit\n')
        fout.close()
        return
    except Exception as err:
        traceback.print_exc()
        raise err



if __name__ == "__main__":
    # Host = HostFile(rootdir,ncpu,0)
    # hostInfo,poolsize,totalproc= Host.set_hostfile()
    pass


