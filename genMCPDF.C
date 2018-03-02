#include <iostream>
#include <fstream>
#include <cmath>

#include "LHAPDF/LHAPDF.h"

using namespace std;

int main(const int argc, const char * argv[]){

  if (argc < 3){
    cout << "./genPDF <pdfset> <NMC> <Q2>" << endl;
    return 0;
  }

  const int NMC = atoi(argv[2]);
  const LHAPDF::PDF * xf0 = LHAPDF::mkPDF(argv[1], 0);
  const LHAPDF::PDF * xf[NMC];
  for (int i = 0; i < NMC; i++)
    xf[i] = LHAPDF::mkPDF(argv[1], i+1);

  double Q2 = atof(argv[3]);
  double X1[500], X2[500];
  for (int i = 0; i < 500; i++){
    X1[i] = pow(10.0, -4.0 + 0.004 * i);
    X2[i] = 0.01 + (1.0 - 0.01) / 499 * i;
  }

  FILE * fs = fopen("xs.dat", "w");
  fprintf(fs, "%s\t%.2f GeV^2\n", argv[1], Q2);
  fprintf(fs, "x  xs-  delta_xs-\n");

  double x, xsm, Exsm;

  for (int i = 0; i < 500; i++){
    x = X1[i];
    xsm = xf0->xfxQ2(3, x, Q2) - xf0->xfxQ2(-3, x, Q2);
    Exsm = 0;
    for (int j = 0; j < NMC; j++){
      Exsm += pow(xf[j]->xfxQ2(3, x, Q2) - xf[j]->xfxQ2(-3, x, Q2) - xsm, 2) / NMC;
    }
    fprintf(fs, "%.6E  %.6E  %.6E\n",
	    x, xsm, sqrt(Exsm));
  }

  for (int i = 0; i < 500; i++){
    x = X2[i];
    xsm = xf0->xfxQ2(3, x, Q2) - xf0->xfxQ2(-3, x, Q2);
    Exsm = 0;
    for (int j = 0; j < NMC; j++){
      Exsm += pow(xf[j]->xfxQ2(3, x, Q2) - xf[j]->xfxQ2(-3, x, Q2) - xsm, 2) / NMC;
    }
    fprintf(fs, "%.6E  %.6E  %.6E\n",
	    x, xsm, sqrt(Exsm));
  }

  fclose(fs);
  
  return 0;
}
