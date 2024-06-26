import argparse
import numpy as np
import pybigtools
import sys

from scipy.signal import find_peaks
from scipy.signal import find_peaks_cwt

"""
find bigwig regions above and below a chosen percentile point.
0.01 works well for very large sequences with a minimum bed window size of 10 
Tricksy numpy method from http://gregoryzynda.com/python/numpy/contiguous/interval/2019/11/29/contiguous-regions.html
takes about 95 seconds for a 17MB test wiggle
"""

class findOut():

    def __init__(self, bwname="test.bw", qfrac=0.01, bedwin=10, bedout="test.bw.bed"):
        self.bedwin = bedwin
        self.qfrac = qfrac
        self.bwf = pybigtools.open(bwname)
        self.chrlist = self.bwf.chroms()
        self.bedf = open(bedout, "w")
        self.makeBed()
        self.bedf.close()

    def makeBedPeaks(self):
        # use scipy find_peaks - not looking for signal "peaks" but interesting...
        bed = []
        chrs = list(self.chrlist.keys())
        chrs.sort()
        print(chrs)

        for chr in chrs:
            chrlen = self.chrlist[chr]
            print(chr, chrlen)
            nb = int(chrlen/self.bedwin) - 1
            bw = self.bwf.values(chr)
            bw = bw[~np.isnan(bw)] # some have NaN if parts of a contig not covered
            bwtop = np.quantile(bw,0.99)
            bwbot = np.quantile(bw,0.01)
            bwmean = np.mean(bw)
            bwstd = np.std(bw)
            bwmax = np.max(bw)
            bwhi, bwhiprop = find_peaks(bw, height=bwtop, width=self.bedwin, distance=20, prominence=(bwtop-bwbot)/10, wlen=100*self.bedwin) 
            bwlo, bwloprop = find_peaks(-bw, height= -1*bwmax, width=self.bedwin, prominence=(bwtop-bwbot)/10, distance=20, wlen=100*self.bedwin) 
            print('n=%d, got %d hi, %d low, topcut=%f, botcut=%f, mean=%f, std=%f, max=%f' % (np.size(bw), len(bwhi), len(bwlo), bwtop, bwbot, bwmean, bwstd, bwmax))
            for i, peak in enumerate(bwhi):
                bed.append((chr, peak, peak+bwhiprop["widths"][i],'%s_%d' % (chr,i), 1))
            for i, peak in enumerate(bwlo):
                bed.append((chr, peak, peak+bwloprop["widths"][i],'%s_%d' % (chr,i),-1))
        bed.sort()
        beds = ['%s\t%d\t%d\t%s\t%d' % x for x in bed]
        self.bedf.write('\n'.join(beds))
        self.bedf.write('\n')
        
    def makeBed(self):
        bed = []
        chrs = list(self.chrlist.keys())
        chrs.sort()
        restab = ["contig\tn\tmean\tstd\tmin\tmax\t"]
        for chr in chrs:
            chrlen = self.chrlist[chr]
            nb = int(chrlen/self.bedwin) - 1
            bw = self.bwf.values(chr)
            bw = bw[~np.isnan(bw)] # some have NaN if parts of a contig not covered
            bwtop = np.quantile(bw, 1.0-self.qfrac)
            bwbot = np.quantile(bw, self.qfrac)
            # http://gregoryzynda.com/python/numpy/contiguous/interval/2019/11/29/contiguous-regions.html
            bwhiex = np.r_[False, bw >= bwtop, False] # extend with 0s
            bwloex = np.r_[False, bw <= bwbot, False]
            bwhid = np.diff(bwhiex)
            bwlod = np.diff(bwloex)
            bwhii = bwhid.nonzero()[0]
            bwloi = bwlod.nonzero()[0]
            bwhi = np.reshape(bwhii, (-1,2))
            bwlo = np.reshape(bwloi, (-1,2))  
            for i, seg in enumerate(bwhi):
                if seg[1] - seg[0] >= self.bedwin:
                    bed.append((chr, seg[0], seg[1], '%s_%d' % (chr,i), bwtop))
            for i, seg in enumerate(bwlo):
                if seg[1] - seg[0] >= self.bedwin:
                    bed.append((chr, seg[0], seg[1], '%s_%d' % (chr,i), bwbot))
            bwmean = np.mean(bw)
            bwstd = np.std(bw)
            bwmax = np.max(bw)
            nrow = np.size(bw)
            bwmin = np.min(bw)
            restab.append('%s\t%d\t%f\t%f\t%f\t%f' % (chr,nrow,bwmean,bwstd,bwmin,bwmax))
        bed.sort()
        beds = ['%s\t%d\t%d\t%s\t%d' % x for x in bed]
        self.bedf.write('\n'.join(beds))
        self.bedf.write('\n')
        self.bedf.close()
        print('\n'.join(restab), '\n')
        print('Wrote %d bed regions' % len(bed))
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    a = parser.add_argument
    a('--minwin',default=10, type=int)
    a('--qfrac',default=0.01, type=float)
    a('--bigwig',default="test.bw")
    a('--bedout',default="test.bw.bed")
    args = parser.parse_args()
    fo = findOut(bwname=args.bigwig, bedwin=args.minwin, qfrac=args.qfrac, bedout=args.bedout)
