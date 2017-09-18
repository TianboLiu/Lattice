#include <iostream>
#include <cmath>

#include "LHAPDF/LHAPDF.h"

#include "Math/GSLIntegrator.h"


using namespace std;

//LHAPDF::PDF * xpdf = LHAPDF::mkPDF("MMHT2014lo68cl", 0); 
LHAPDF::PDF * xpdf = LHAPDF::mkPDF("NNPDF30_lo_as_0118", 0);
double Q = 10.0;

double ssbar(const double x, void * par){
  if (x < 1e-8) return 0.0;
  double value = (xpdf->xfxQ(3, x, Q) - xpdf->xfxQ(-3, x, Q)) / x;
  return value;
}



int main(const int argc, const char * argv[]){

  double par = 0.0;
  ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, 0.0, 1.0e-6);
  ig.SetFunction(&ssbar, &par);
  double result = ig.Integral(atof(argv[1]), atof(argv[2]));

  cout << result << endl;

  return 0;
}
