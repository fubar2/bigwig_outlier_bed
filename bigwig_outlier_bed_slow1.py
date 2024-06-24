import os
import pyBigWig


class findOut():

    def __init__(self, bwname="test.bw", bedname="test.bed", sd_lo=3, sd_hi=3, bedwin=10):
        self.bedwin = bedwin
        self.bwf = pyBigWig.open(bwname, "r")
        self.chrlist = self.bwf.chroms()
        self.bedf = open(bedname, "w")
        self.sd_lo = sd_lo
        self.sd_hi = sd_hi
        self.makeBed(self.chrlist)

    def makeBed(self, chrlist):
        bed = []
        for chr in chrlist:
            print(chr)
            gm = self.bwf.stats(chr, 0, type="mean")[0]
            gstd = self.bwf.stats(chr, 0, type="std")[0]
            cutlow = gm - (self.sd_lo * gstd)
            cuthi = gm + (self.sd_hi * gstd)
            chr_len = self.chrlist[chr]
            nb = chr_len // self.bedwin
            means = self.bwf.stats(chr, 0, chr_len, type="mean", nBins=nb)
            inlo = False
            inhi = False
            reg_start = None
            reg_end = None
            reg_means = []
            print('got %d means, gm=%f, lo=%f, hi=%f' % (len(means), gm, cutlow, cuthi))
            for i, m in enumerate(means):
                bend = min(chr_len, (i + 1) * self.bedwin)
                featname = "%s_%d" % (chr,i)
                if m and (m < cutlow or m > cuthi):
                    if inlo:
                        if m < cutlow:  # extend
                            reg_end = bend
                            reg_means.append(m)
                        else:  # high so close
                            rm = sum(reg_means) / len(reg_means)
                            bed.append(
                                "%s\t%d\t%d\t%s\t%.3f\n"
                                % (chr, reg_start, reg_end, featname, rm)
                            )
                            inlo = False
                            reg_means = []
                    elif inhi:
                        if m > cuthi:  # extend
                            reg_end = bend
                            reg_means.append(m)
                        else:
                            rm = sum(reg_means) / len(reg_means)
                            bed.append(
                                "%s\t%d\t%d\t%s\t%.3f\n"
                                % (chr, reg_start, reg_end, featname, rm)
                            )
                            inhi = False
                            reg_means = []
                    elif m < cutlow:  # start new low region
                        inlo = True
                        reg_start = i * self.bedwin
                        reg_end = bend
                        reg_means = [m]
                    elif m > cuthi:
                        inhi = True
                        reg_start = i * self.bedwin
                        reg_end = bend
                        reg_means = [m]
                else:  # not out of range - write current extended bed region
                    if inhi or inlo:
                        inhi = False
                        inlo = False
                        rm = sum(reg_means) / len(reg_means)
                        bed.append(
                            "%s\t%d\t%d\t%s\t%.3f\n"
                            % (chr, reg_start, reg_end, featname, rm)
                        )
                        reg_means = []
                        reg_start = None
                        reg_end = None
        self.bedf.write(''.join(bed))
        self.bedf.close()


if __name__ == "__main__":
    fo = findOut(bwname="test.bw", bedname="test.bed", sd_lo=2, sd_hi=2, bedwin=100)
