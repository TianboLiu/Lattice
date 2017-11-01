#!/usr/bin/env python
import os, sys
import numpy as np
import scipy as sp
import lhapdf


def genPDFs(filename, pdfset, Q, Nsets):
    f = open(filename, "w")
    f.write("%s\t %.1f\n" %(pdfset,Q))
    f.write("x\t xf(uv-dv) dxf(uv-dv)\n")
    pdfs={}
    for i in range(Nsets):
        pdfs[i] = lhapdf.mkPDF(pdfset, i)
    X=np.linspace(0.005, 0.995, num=100, endpoint=True)
    for x in X:
        xf = (pdfs[0].xfxQ(2, x, Q) - pdfs[0](-2, x, Q)) - (pdfs[0].xfxQ(1, x, Q) - pdfs[0](-1, x, Q))
        xfall = []
        for i in range(Nsets):
            xfall.append(pdfs[i].xfxQ(2, x, Q) - pdfs[i](-2, x, Q)) - (pdfs[i].xfxQ(1, x, Q) - pdfs[i](-1, x, Q))
        dxf = np.std(xfall)
        f.write("%.3E\t%.3E\t%.3E\n" % (x, xf, dxf))
    f.close()
    return


if __name__=="__main__":
    if len(sys.argv) < 5:
        print "./genPDFs <filename> <pdfset> <Q> <Nsets>"
        sys.exit()
    else:
        genPDFs(sys.argv[1], sys.argv[2], float(sys.argv[3]), int(sys.argv[4]))
        print "Bye!"

