"""
Ross Lazarus June 2024 for VGP 

Bigwigs are great, but hard to reliably "see" small low coverage or small very high coverage regions.
Colouring is not feasible without a new plugin, so this code will find bigwig regions above and below a chosen percentile point.

Multiple bigwigs **with the same reference** can be combined - bed segments will be named appropriately
Combining multiple references works but is silly because only display will rely on one reference so others will not be shown...

0.01 works well for very large sequences with a minimum bed window size of 10 
Tricksy numpy method from http://gregoryzynda.com/python/numpy/contiguous/interval/2019/11/29/contiguous-regions.html
takes about 95 seconds for a 17MB test wiggle

JBrowse2 bed normally displays ignore the score, so probably need options for separate low/high bed file outputs.

Update june 30 2024: wrote a 'no-build' plugin for beds to display red/blue if >0/<0 so those are used for scores
Bed interval naming must be short for JB2 but needs input bigwig name and (lo or hi).

"""


import argparse
import numpy as np
import pybigtools
import sys

class findOut():

    def __init__(self, bwname=["test.bw",], bwlabels=["test",], qlo=0.01, qhi=0.99, bedwin=10, bedout="test.bw.bed"):
        self.bedwin = bedwin
        self.qhi = qhi
        self.qlo = qlo
        self.makeBed(bwname[0], bwlabels[0], bedout)

    def processVals(self, bw, isTop):
        # http://gregoryzynda.com/python/numpy/contiguous/interval/2019/11/29/contiguous-regions.html
        if isTop:
            bwex = np.r_[False, bw >= self.bwtop, False] # extend with 0s
        else:
            bwex = np.r_[False, bw <= self.bwbot, False]
        bwexd = np.diff(bwex)
        bwexdnz = bwexd.nonzero()[0]
        bwregions = np.reshape(bwexdnz, (-1,2))
        return bwregions
        
    def makeBed(self, bwnames, bwlabels, bedout):
        bedf = open(bedout, "w")
        bed = []
        bwlabels = bwlabels.split(',')
        bwnames = bwnames.split(',')
        print('bwnames=', bwnames, "bwlabs=", bwlabels)
        for i, bwname in enumerate(bwnames):
            bwlabel = bwlabels[i].replace(" ",'')
            bwf = pybigtools.open(bwname)
            chrlist = bwf.chroms()
            chrs = list(chrlist.keys())
            chrs.sort()
            restab = ["contig\tn\tmean\tstd\tmin\tmax\tqtop\tqbot"]
            for chr in chrs:
                bw = bwf.values(chr)
                bw = bw[~np.isnan(bw)] # some have NaN if parts of a contig not covered
                self.bwtop = np.quantile(bw, self.qhi)
                self.bwbot = np.quantile(bw, self.qlo)
                if self.qhi is not None:
                    bwhi = self.processVals(bw, isTop=True)
                    for i, seg in enumerate(bwhi):
                        if seg[1] - seg[0] >= self.bedwin:
                            bed.append((chr, seg[0], seg[1], '%s_hi' % (bwlabel), 1))
                if self.qlo is not None:
                    bwlo = self.processVals(bw, isTop=False)            
                    for i, seg in enumerate(bwlo):
                        if seg[1] - seg[0] >= self.bedwin:
                            bed.append((chr, seg[0], seg[1], '%s_lo' % (bwlabel), -1))
                bwmean = np.mean(bw)
                bwstd = np.std(bw)
                bwmax = np.max(bw)
                nrow = np.size(bw)
                bwmin = np.min(bw)
                restab.append('%s\t%d\t%f\t%f\t%f\t%f\t%f\t%f' % (chr,nrow,bwmean,bwstd,bwmin,bwmax,self.bwtop,self.bwbot))
            print('\n'.join(restab), '\n')
        bed.sort()
        beds = ['%s\t%d\t%d\t%s\t%d' % x for x in bed]
        bedf.write('\n'.join(beds))
        bedf.write('\n')
        bedf.close()
        print('Wrote %d bed regions' % len(bed))
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    a = parser.add_argument
    a('-m','--minwin',default=10, type=int)
    a('-l', '--qlo',default=None, type=float)
    a('-i', '--qhi',default=0.99, type=float)
    a('-b', '--bigwig', nargs='+')
    a('-n', '--bigwiglabels', nargs='+')
    a('-o', '--bedout', required=True)
    args = parser.parse_args()
    fo = findOut(bwname=args.bigwig, bwlabels=args.bigwiglabels, bedwin=args.minwin, qlo=args.qlo, qhi=args.qhi, bedout=args.bedout)
