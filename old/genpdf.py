#!/usr/bin/env python
import os, sys
import numpy as np
import scipy as sp
import lhapdf

nnpdf = lhapdf.mkPDF("NNPDF30_lo_as_0118", 0)
mmht = lhapdf.mkPDF("MMHT2014lo68cl", 0)
cj15 = lhapdf.mkPDF("CJ15lo", 0)
ct14 = lhapdf.mkPDF("CT14lo", 0)


f = open("pdfgrid.txt", "w")
f.write("x \t NNPDF3 \t MMHT14 \t CJ15 \t CT14\n")

xlist = 10**np.linspace(np.log10(2.035e-7), 0, 2000)
nnpdflist = xlist.copy()
mmhtlist = xlist.copy()
cj15list = xlist.copy()
ct14list = xlist.copy()

for i in range(2000):
    nnpdflist[i] = nnpdf.xfxQ(3, xlist[i], 2.0) - nnpdf.xfxQ(-3, xlist[i], 2.0)
    mmhtlist[i] = mmht.xfxQ(3, xlist[i], 2.0) - mmht.xfxQ(-3, xlist[i], 2.0)
    cj15list[i] = cj15.xfxQ(3, xlist[i], 2.0) - cj15.xfxQ(-3, xlist[i], 2.0)
    ct14list[i] = ct14.xfxQ(3, xlist[i], 2.0) - ct14.xfxQ(-3, xlist[i], 2.0)

for i in range(2000):
    f.write("%.6E \t %.6E \t %.6E \t %.6E \t %.6E\n" %(xlist[i], nnpdflist[i]/xlist[i], mmhtlist[i]/xlist[i], cj15list[i]/xlist[i], ct14list[i]/xlist[i]))

f.close()

print "Bye!"

