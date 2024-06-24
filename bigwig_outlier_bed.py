import os
import pybigtools
import numpy as np
from scipy.signal import find_peaks


class findOut():

    def __init__(self, bwname="test.bw", bedname="test.bed", sd_lo=3, sd_hi=3, bedwin=10):
        self.bedwin = bedwin
        self.bwf = pybigtools.open(bwname)
        self.chrlist = self.bwf.chroms()
        self.bedf = open(bedname, "w")
        self.sd_lo = sd_lo
        self.sd_hi = sd_hi
        self.makeBed()

    def makeBed(self):
        bed = []
        chrs = list(self.chrlist.keys())
        chrs.sort()
        print(chrs)

        for chr in chrs:
            chrlen = self.chrlist[chr]
            print(chr, chrlen)
            nb = int(chrlen/self.bedwin) - 1
            bw = self.bwf.values(chr)
            bwtop = np.quantile(bw,0.99)
            bwbot = np.quantile(bw,0.01)
            bwmean = np.mean(bw)
            bwstd = np.std(bw)
            bwmax = np.max(bw)
            bwhi, bwhiprop = find_peaks(bw, height=bwtop, width=self.bedwin) 
            print("make bwinv")
            bwinv = bw - bwmax
            print("make bwlo")
            bwlo, bwloprop = find_peaks(bwinv, height=bwbot-bwmax, width=self.bedwin) 
            print('n=%d, got %d hi, %d low, topcut=%f, botcut=%f, mean=%f, std=%f' % (np.size(bw), len(bwhi), len(bwlo), bwtop, bwbot, bwmean, bwstd))
            print('bwhi', str(bwhi[:10]))
            print('hiwidths', bwhiprop["widths"][:10])
            print('bwlo', str(bwlo[:10]))
            print('lowidths', bwloprop["widths"][:10])



if __name__ == "__main__":
    fo = findOut(bwname="test.bw", bedname="test.bed", sd_lo=3, sd_hi=3, bedwin=20)
