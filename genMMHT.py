#!/usr/bin/env python
import os, sys
import numpy as np
import scipy as sp
import lhapdf

nnpdf = lhapdf.mkPDF("MMHT2014lo68cl", 0)

f = open("mmhtgrid10.txt", "w")
f.write("x")
for i in range(51):
    f.write("\t%d" %i)
f.write("\n")

nnpdflist = np.arange(52*1000, dtype=np.float64).reshape(1000, 52)

xlist = 10**np.linspace(np.log10(2.035e-7), 0, 1000)

for i in range(1000):
    nnpdflist[i,0] = xlist[i]

Q = 10.0
for id in range(51):
    nnpdf = lhapdf.mkPDF("MMHT2014lo68cl", id)
    for i in range(1000):
        nnpdflist[i, id+1] = (nnpdf.xfxQ(3, xlist[i], Q) - nnpdf.xfxQ(-3, xlist[i], Q)) / xlist[i];

for i in range(1000):
    for j in range(52):
        f.write("%.6E\t" % nnpdflist[i,j])
    f.write("\n")

f.close()

print "Bye!"


