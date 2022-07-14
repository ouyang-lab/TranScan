# TranScan software for the paper entitled "Translocation Detection from Hi-C data via Scan Statistics" by Anthony Cheng, Disheng Mao, Yuping Zhang, Joseph Glaz, Zhengqing Ouyang
# The TranScan R package detects translocations in genome-wide proximity ligation data such as Hi-C. The accompanying translocation breakpoint finder is post-processing scripts that are written in Python.
# License: The implementations written for this project is covered by the GNU General Public License, version 3.0 (GPL-3.0)

import matplotlib
matplotlib.use("Agg")

import sys
import argparse
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


def get_sizes(f_sizes):
    sizes = {}
    with open(f_sizes, "r") as f:
        for line in f:
            row = line.strip("\r\n").split("\t")
            sizes[row[0]] = int(row[1])
    return sizes


def load_breakpoints(f_bp):
    dict_bp = {}
    with open(f_bp, "r") as f:
        for line in f:
            row = line.strip("\r\n").split("\t")
            c1, s1, c2, s2 = row[0], int(row[1]), row[2], int(row[3])
            dict_bp[c1] = s1
            dict_bp[c2] = s2
    return dict_bp

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", default=None, help="RR.txt")
    parser.add_argument("-k", default=None, help="KDE.txt")
    parser.add_argument("-sizes", default=None, help="sizes")
    parser.add_argument("-c1", default=None, help="chrom1")
    parser.add_argument("-c2", default=None, help="chrom2")
    parser.add_argument("-o", default=None, help="output")
    parser.add_argument("-bp", default=None, help="Breakpoints File (if any)")
    
    args = parser.parse_args()

    p_chrom1 = args.c1
    p_chrom2 = args.c2
    f_txt = args.i
    f_kde = args.k
    f_sizes = args.sizes
    f_output = args.o
    f_bp = args.bp
    
    p_gap = 2500000
    p_quantitest = 256
    p_half = p_gap/2
    p_bw_win = 1
    p_threshold = 0  # 36
    
    
    if f_bp:
        dict_bp = load_breakpoints(f_bp)
    
    
    data = np.loadtxt(f_txt)
    kde = np.loadtxt(f_kde)
    
    sizes = get_sizes(f_sizes)
    
    len_chrA = sizes[p_chrom1]
    len_chrB = sizes[p_chrom2]

    p_eff_res = (len_chrA + len_chrB + p_gap)/float(p_quantitest)
    prop_chrA = len_chrA/float(len_chrA+len_chrB)
    prop_chrB = len_chrB/float(len_chrA+len_chrB)

    inter_chrA = prop_chrA * (len_chrA + len_chrB + p_gap)
    inter_chrB = prop_chrB * (len_chrB + len_chrB + p_gap)
    

    bins_chrA = int(round(inter_chrA/p_eff_res))
    bins_chrB = int(p_quantitest - bins_chrA)
    inter = data[:bins_chrA,bins_chrA:]
    
    
    rows_ignore = np.where(inter[p_bw_win:-p_bw_win,:].sum(axis=1) == 0)[0]
    cols_ignore = np.where(inter[:,p_bw_win:-p_bw_win].sum(axis=0) == 0)[0]

    interkde = kde[:bins_chrA,bins_chrA:][p_bw_win:-p_bw_win,p_bw_win:-p_bw_win]    
    
    
    X = []
    Y = []
    pos_X = []
    pos_Y = []
    XY = []
    for i in range(p_bw_win, bins_chrA-1):
        state_change_xb = np.count_nonzero(inter[i-p_bw_win,:] != inter[i,:])
        state_change_xf = np.count_nonzero(inter[i,:] != inter[i+p_bw_win,:])
        binth1 = int((i-1)*p_eff_res)
        binth2 = int(i*p_eff_res)
        X.append(state_change_xb+state_change_xf)
        pos_X.append((p_chrom1, binth1, binth2, i))
        
    for i in range(p_bw_win, bins_chrB-1):
        state_change_yb = np.count_nonzero(inter[:,i-p_bw_win] != inter[:,i])
        state_change_yf = np.count_nonzero(inter[:,i] != inter[:,i+p_bw_win])
        binth1 = int((i-1)*p_eff_res)
        binth2 = int(i*p_eff_res)
        Y.append(state_change_yb+state_change_yf)
        pos_Y.append((p_chrom2, binth1, binth2, i))

    
    # calc state change product
    XY = np.outer(X, Y)
    
    # ignore those with contiguous gaps
    XY[rows_ignore,:] = 0
    XY[:,cols_ignore] = 0
    
    # zero out non-rejected region
    XY[np.where(inter[p_bw_win:-p_bw_win,p_bw_win:-p_bw_win] == 0)] = 0

    # look for the maximum point
    index = np.unravel_index(np.argmax(XY, axis=None), XY.shape)
    no_of_pixels = np.max(XY)
    
    # output locus pair state change product crosses threshold
    if no_of_pixels >= p_threshold:
        print "\t".join(map(str, 
                            list(pos_X[index[0]][:3])+list(pos_Y[index[1]][:3])+[no_of_pixels]))
    
    """ Output """
    if f_output:
        plt.imshow(inter)

        if f_bp:  # plot original breakpoint location
            x = dict_bp[p_chrom1]/p_eff_res
            y = dict_bp[p_chrom2]/p_eff_res
            plt.scatter(x, y, s=5, color="black")
        
        plt.axhline(pos_X[index[0]][3], linestyle="--", linewidth=1.0)
        plt.axvline(pos_Y[index[1]][3], linestyle="--", linewidth=1.0)
        plt.savefig(f_output, dpi=300)
    
