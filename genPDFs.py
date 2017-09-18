#!/usr/bin/env python
import os, sys
import numpy as np
import scipy as sp
import lhapdf


def genPDFs(filename, pdfset, flavor, Q, nlines):
    f = open(filename, "w")
    f.write("x\t")
    for i in range(nlines+1):
        f.write("%d\t" %i)
        f.write("\n")
    pdflist = np.arange((nlines+2)*1000, dtype=np.float64).reshape(1000, nlines+2)
    xlist = 10**np.linspace(np.log10(2.035e-7), 0, 1000)
    for i in range(1000):
        pdflist[i,0] = xlist[i]
    for id in range(nlines+1):
        pdf = lhapdf.mkPDF(pdfset, id)
        for i in range(1000):
            pdflist[i, id+1] = pdf.xfxQ(flavor, xlist[i], Q) / xlist[i];
    for i in range(1000):
        for j in range(nlines+2):
            f.write("%.6E\t" % pdflist[i,j])
            f.write("\n")
    f.close()
    return


if __name__=="__main__":
    if len(sys.argv) < 6:
        print "./genPDFs <filename> <pdfset> <flavor> <Q> <nlines>"
        sys.exit()
    else:
        genPDFs(sys.argv[1], sys.argv[2], int(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5]))
        print "Bye!"
        


