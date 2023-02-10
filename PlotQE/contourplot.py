# python3 script
# JamesBourbin in 20220522 finished V1.0
# contour and bubble for Energy-Stein_OP-DOS 2D Energy profile

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import sys
# from scipy import interpolate
# from skimage import transform

# mpl.rcParams['agg.path.chunksize'] = 1000
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams["font.size"] = 16

# inputfile="out_204060.txt"
inputfile = "out_080824b.txt"
# inputfile format is fixed; gen-ed from computeQ.py
class QEdata:
    def __init__(self, data_OP:dict, xname="OP", yname="Energy(eV)", zname="DOS"):
        self.data_Q = data_OP
        self.xname = xname
        self.yname = yname
        self.zname = zname
        self.data_df = pd.DataFrame()

    def set_dataframe(self):
        '''data transfer to Q-E-DOS format, table: OPs,E,dE,DOS 
        
        '''
        # sort Q_pair and output
        self.delta_yname = "Delta-"+self.yname
        data_list = list(self.data_Q.items())
        # sort by energy
        data_list.sort(key=lambda x:x[0][1])
        min_energy = data_list[0][0][1] # get minimum
        all_data = []
        for Q_pair, count in data_list:
            one_data = np.array([Q_pair[0], Q_pair[1], np.round(Q_pair[1]-min_energy,2), count])
            all_data.append(one_data)
        all_data = np.array(all_data)
        self.global_mini_energy = min_energy
        self.data_df = pd.DataFrame(all_data, columns=[self.xname, self.yname, self.delta_yname, self.zname])
        
    def print_dataframe(self):
        filename = self.xname+"-"+self.yname[:6]+"-"+self.zname+".csv"
        self.data_df.to_csv(filename)
        
    def bubble_scatter(self):
        '''plot bubble scatter instead of contour, for fewer points less than 100000'''
        figname="%s-%s-%s_bubble.png"%(self.xname, self.yname[:6], self.zname)
        x = self.data_df[self.xname]
        y = self.data_df[self.delta_yname]
        z = self.data_df[self.zname]
        
        plt.figure(figsize=(8,12))
        plt.scatter(x,y, s=14, c=np.log(z), cmap='jet', marker='o')
        cb = plt.colorbar()
        cb.set_label("lg(%s)"%self.zname)
        plt.xlabel(self.xname)
        plt.ylabel(self.yname)
        plt.savefig(figname, dpi=200)
    
    def coutour_plot(self):
        '''plot imshow-contour, for points more than 10w'''
        figname_20="%s-%s-%s_contour_20.png"%(self.xname, self.yname[:6], self.zname)
        figname_100="%s-%s-%s_contour_100.png"%(self.xname, self.yname[:6], self.zname)
        x = self.data_df[self.xname]
        y = self.data_df[self.delta_yname]
        # z = self.data_df[self.zname]
        x_sorted = sorted(list(set(x)))
        y_sorted = sorted(list(set(y)))
        X, Y = np.meshgrid(x_sorted, y_sorted)
        Z = np.zeros(np.shape(X))
        for i, xi in enumerate(x_sorted):
            for j, yj in enumerate(y_sorted):
                z_point = np.array((self.data_df[
                    (self.data_df[self.xname]==xi) & 
                        (self.data_df[self.delta_yname]==yj) ][self.zname]),dtype=int)
                if len(z_point) == 0:
                    Z[j][i] = 1
                else:
                    Z[j][i] = int(z_point[0])
        # split 20
        plt.figure(figsize=(8,12))
        plt.contourf(X, Y, np.log(Z), 20, alpha=0.75, cmap='jet')
        plt.xlabel(self.xname)
        plt.ylabel(self.yname)
        cb = plt.colorbar()
        cb.set_label("lg(%s)"%self.zname)
        plt.savefig(figname_20, dpi=200)
        # split 100
        plt.figure(figsize=(8,12))
        plt.contourf(X, Y, np.log(Z), 100, alpha=0.75, cmap='jet')
        plt.xlabel(self.xname)
        plt.ylabel(self.yname)
        cb = plt.colorbar()
        cb.set_label("lg(%s)" % self.zname)
        plt.savefig(figname_100, dpi=200)

        

# read data
if __name__ == "__main__":
    if len(sys.argv) >= 2:
        filename = sys.argv[1]
    else:
        filename = inputfile
    data_SteinQ_OP2 = {}
    data_SteinQ_OP4 = {}
    data_SteinQ_OP6 = {}
    data_DistSteinQ_OP2 = {}
    data_DistSteinQ_OP4 = {}
    data_DistSteinQ_OP6 = {}
    # adjust round precision
    energy_round = 1
    Q_round = 2
    print("read data from file %s"%filename)
    with open(filename, 'r') as fo:
        for line in fo:
            # read stein-Q data
            if "of Stein-Q" in line:
                fo.readline()
                data_list = fo.readline().strip().split()
                while len(data_list) == 5:
                    energy = np.round(np.float64(data_list[1]),energy_round)
                    Q2 = np.round(np.float64(data_list[2]),Q_round)
                    Q4 = np.round(np.float64(data_list[3]),Q_round)
                    Q6 = np.round(np.float64(data_list[4]),Q_round)
                    Q2_pair = (Q2, energy)
                    Q4_pair = (Q4, energy)
                    Q6_pair = (Q6, energy)
                    data_SteinQ_OP2[Q2_pair] = data_SteinQ_OP2.get(Q2_pair, 0) + 1
                    data_SteinQ_OP4[Q4_pair] = data_SteinQ_OP4.get(Q4_pair, 0) + 1
                    data_SteinQ_OP6[Q6_pair] = data_SteinQ_OP6.get(Q6_pair, 0) + 1
                    data_list = fo.readline().strip().split()
            # read distance-weighted stein-Q data
            elif "of Distance-weighted Stein-Q" in line:
                fo.readline()
                data_list = fo.readline().strip().split()
                while len(data_list) == 5:
                    energy = np.round(np.float64(data_list[1]),energy_round)
                    Q2 = np.round(np.float64(data_list[2]),Q_round)
                    Q4 = np.round(np.float64(data_list[3]),Q_round)
                    Q6 = np.round(np.float64(data_list[4]),Q_round)
                    Q2_pair = (Q2, energy)
                    Q4_pair = (Q4, energy)
                    Q6_pair = (Q6, energy)
                    data_DistSteinQ_OP2[Q2_pair] = data_DistSteinQ_OP2.get(Q2_pair, 0) + 1
                    data_DistSteinQ_OP4[Q4_pair] = data_DistSteinQ_OP4.get(Q4_pair, 0) + 1
                    data_DistSteinQ_OP6[Q6_pair] = data_DistSteinQ_OP6.get(Q6_pair, 0) + 1
                    data_list = fo.readline().strip().split()
    # get data print
    all = {
        "OP2": data_SteinQ_OP2,
        "OP4": data_SteinQ_OP4,
        "OP6": data_SteinQ_OP6,
        "dist_OP2": data_DistSteinQ_OP2,
        "dist_OP4": data_DistSteinQ_OP4,
        "dist_OP6": data_DistSteinQ_OP6
    }
    for xname, data in all.items():
        print("analysis Energy-%s-DOS data, print csv and plot bubble/contourf"%xname)
        QE = QEdata(data_OP=data, xname=xname)
        QE.set_dataframe()
        QE.print_dataframe()
        QE.bubble_scatter()
        QE.coutour_plot()
    print("Done!")
