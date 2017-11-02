#!/usr/bin/env python
import os, sys
import numpy as np
import scipy as sp
import lhapdf


def genPDFs(filename, pdfset, Q, Nsets):
    f = open(filename, "w")
    f.write("%s\t %.1f\n" %(pdfset,Q))
    f.write("x\t x*uv\t dx*uv\t x*dv\t dx*dv\t x*(uv-dv)\t dx*(uv-dv)\n")
    pdfs={}
    for i in range(Nsets):
        pdfs[i] = lhapdf.mkPDF(pdfset, i)
    X=np.linspace(0.0005, 0.9995, num=1000, endpoint=True)
    for x in X:
        xf = pdfs[0].xfxQ(2, x, Q) - pdfs[0].xfxQ(-2, x, Q) - pdfs[0].xfxQ(1, x, Q) + pdfs[0].xfxQ(-1, x, Q)
        xfu = pdfs[0].xfxQ(2, x, Q) - pdfs[0].xfxQ(-2, x, Q)
        xfd = pdfs[0].xfxQ(1, x, Q) - pdfs[0].xfxQ(-1, x, Q)
        xfall = []
        xfallu = []
        xfalld = []
        for i in range(Nsets):
            xfall.append((pdfs[i].xfxQ(2, x, Q) - pdfs[i].xfxQ(-2, x, Q)) - (pdfs[i].xfxQ(1, x, Q) - pdfs[i].xfxQ(-1, x, Q)))
            xfallu.append(pdfs[i].xfxQ(2, x, Q) - pdfs[i].xfxQ(-2, x, Q))
            xfalld.append(pdfs[i].xfxQ(1, x, Q) - pdfs[i].xfxQ(-1, x, Q))
        dxf = np.std(np.array(xfall))
        dxfu = np.std(np.array(xfallu))
        dxfd = np.std(np.array(xfalld))
        f.write("%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\t%.3E\n" % (x, xfu, dxfu, xfd, dxfd, xf, dxf))
    f.close()
    return


if __name__=="__main__":
    if len(sys.argv) < 5:
        print "./genPDFs <filename> <pdfset> <Q> <Nsets>"
        sys.exit()
    else:
        genPDFs(sys.argv[1], sys.argv[2], float(sys.argv[3]), int(sys.argv[4]))
        print "Bye!"

