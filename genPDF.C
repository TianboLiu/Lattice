#include <iostream>
#include <cmath>

#include "LHAPDF/LHAPDF.h"

#include "Math/GSLIntegrator.h"


using namespace std;

LHAPDF::PDF * xpdf[101];

double ssbar(const double x, void * par){
  double Q = ((double *) par)[0];
  //if (x < 1e-8) return 0.0;
  double value = (xpdf[0]->xfxQ(3, x, Q) - xpdf[0]->xfxQ(-3, x, Q)) / x;
  return value;
}

double check(double xmin, double Q){
  ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, 0.0, 1e-6);
  ig.SetFunction(&ssbar, &Q);
  double region1 = ig.Integral(xmin, 1e-4);
  double region2 = ig.Integral(1e-4, 1.0);
  return region1 + region2;
}
  
int main(const int argc, const char * argv[]){

  if (argc < 5){
    cout << "./genPDF <opt> <PDFset> <Q2> <xmin>" << endl;
    return 0;
  }

  const int opt = atoi(argv[1]);

  if (opt == 0){
    xpdf[0] = LHAPDF::mkPDF(argv[2], 0);
    double Q = sqrt(atof(argv[3]));
    double xmin = atof(argv[4]);
    cout << Q << "   " << xmin << endl;
    cout << check(xmin, Q) << endl;
    return 0;
  }

  if (opt == 1){
    if (argc < 6){
      cout << "./genPDF <opt> <PDFset> <Q2> <xmin> <Ncol>" << endl;
      return 0;
    }
    double Q = sqrt(atof(argv[3]));
    double xmin = atof(argv[4]);
    int Ncol = atoi(argv[5]);
    for (int i = 0; i <= Ncol; i++){
      xpdf[i] = LHAPDF::mkPDF(argv[2], i);
    }
    FILE * fs = fopen("fs.dat", "w");
    FILE * fsbar = fopen("fsbar.dat", "w");
    fprintf(fs, "x  ");
    fprintf(fsbar, "x  ");
    for (int i = 0; i <= Ncol; i++){
      fprintf(fs, "%d  ", i);
      fprintf(fsbar, "%d  ", i);
    }
    fprintf(fs, "\n");
    fprintf(fsbar, "\n");

    double X1[100], X2[100];
    double steplog = (log(0.1) -  log(xmin)) / 100.0;
    double step = (1.0 - 0.1) / 100.0;
    for (int i = 0; i < 100; i++){
      X1[i] = exp(log(xmin) + steplog * i);
      X2[i] = 0.1 + step * i;
    }
    for (int i = 0; i < 100; i++){
      fprintf(fs, "%.6E  ", X1[i]);
      fprintf(fsbar, "%.6E  ", X1[i]);
      for (int j = 0; j <= Ncol; j++){
	fprintf(fs, "%.6E  ", xpdf[j]->xfxQ(3, X1[i], Q) / X1[i]);
	fprintf(fsbar, "%.6E  ", xpdf[j]->xfxQ(-3, X1[i], Q) / X1[i]);
      }
      fprintf(fs, "\n");
      fprintf(fsbar, "\n");
    }
    for (int i = 0; i < 100; i++){
      fprintf(fs, "%.6E  ", X2[i]);
      fprintf(fsbar, "%.6E  ", X2[i]);
      for (int j = 0; j <= Ncol; j++){
	fprintf(fs, "%.6E  ", xpdf[j]->xfxQ(3, X2[i], Q) / X2[i]);
	fprintf(fsbar, "%.6E  ", xpdf[j]->xfxQ(-3, X2[i], Q) / X2[i]);
      }
      fprintf(fs, "\n");
      fprintf(fsbar, "\n");
    }
    fprintf(fs, "%.6E  ", 1.0);
    fprintf(fsbar, "%.6E  ", 1.0);
    for (int i = 0; i <= Ncol; i++){
      fprintf(fs, "%.6E  ", 0.0);
      fprintf(fsbar, "%.6E  ", 0.0);
    }
    fprintf(fs, "\n");
    fprintf(fsbar, "\n");
    fclose(fs);
    fclose(fsbar);
  }
  
  if (opt == 2){
    if (argc < 6){
      cout << "./genPDF <opt> <PDFset> <Q2> <xmin> <Ncol>" << endl;
      return 0;
    }
    double Q = sqrt(atof(argv[3]));
    double xmin = atof(argv[4]);
    int Ncol = atoi(argv[5]);
    for (int i = 0; i <= Ncol; i++){
      xpdf[i] = LHAPDF::mkPDF(argv[2], i);
    }
    FILE * fsm = fopen("fsm.dat", "w");
    fprintf(fsm, "x  ");
    for (int i = 0; i <= Ncol; i++){
      fprintf(fsm, "%d  ", i);
    }
    fprintf(fsm, "\n");
 
    double X1[100], X2[100];
    double steplog = (log(0.1) -  log(xmin)) / 100.0;
    double step = (1.0 - 0.1) / 100.0;
    for (int i = 0; i < 100; i++){
      X1[i] = exp(log(xmin) + steplog * i);
      X2[i] = 0.1 + step * i;
    }
    for (int i = 0; i < 100; i++){
      fprintf(fsm, "%.6E  ", X1[i]);
      for (int j = 0; j <= Ncol; j++){
	fprintf(fsm, "%.6E  ", (xpdf[j]->xfxQ(3, X1[i], Q) - xpdf[j]->xfxQ(-3, X1[i], Q)) / X1[i]);
      }
      fprintf(fsm, "\n");
    }
    for (int i = 0; i < 100; i++){
      fprintf(fsm, "%.6E  ", X2[i]);
      for (int j = 0; j <= Ncol; j++){
	fprintf(fsm, "%.6E  ", (xpdf[j]->xfxQ(3, X2[i], Q) - xpdf[j]->xfxQ(-3, X2[i], Q)) / X2[i]);
      }
      fprintf(fsm, "\n");
    }
    fprintf(fsm, "%.6E  ", 1.0);
    for (int i = 0; i <= Ncol; i++){
      fprintf(fsm, "%.6E  ", 0.0);
    }
    fprintf(fsm, "\n");
    fclose(fsm);
  } 

  return 0;
}
