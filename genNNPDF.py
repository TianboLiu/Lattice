#!/usr/bin/env python
import os, sys
import numpy as np
import scipy as sp
import lhapdf

nnpdf = lhapdf.mkPDF("NNPDF30_lo_as_0118", 0)

f = open("nnpdfgrid10.txt", "w")
f.write("x")
for i in range(101):
    f.write("\t%d" %i)
f.write("\n")

nnpdflist = np.arange(102*1000, dtype=np.float64).reshape(1000, 102)

xlist = 10**np.linspace(np.log10(0.515e-7), 0, 1000)

for i in range(1000):
    nnpdflist[i,0] = xlist[i]

Q = 10.0
for id in range(101):
    nnpdf = lhapdf.mkPDF("NNPDF30_lo_as_0118", id)
    for i in range(1000):
        nnpdflist[i, id+1] = (nnpdf.xfxQ(3, xlist[i], Q) - nnpdf.xfxQ(-3, xlist[i], Q)) / xlist[i];

for i in range(1000):
    for j in range(102):
        f.write("%.6E\t" % nnpdflist[i,j])
    f.write("\n")

f.close()

print "Bye!"


