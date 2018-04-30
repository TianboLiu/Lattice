#include <iostream>
#include <fstream>
#include <cmath>

#include "LHAPDF/LHAPDF.h"

#include "Math/GSLIntegrator.h"
#include "Math/Integrator.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/WrappedTF1.h"
#include "Math/WrappedParamFunction.h"
#include "TF1.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
#include "Math/Minimizer.h"

using namespace std;

LHAPDF::PDF * xf;

double Q2_[10], Fs_[10], dFs_[10];
int LoadData(){
  ifstream infile("f1smod.txt");
  infile >> Q2_[0] >> Fs_[0] >> dFs_[0];
  for (int i = 0; i < 10; i++)
    infile >> Q2_[i] >> Fs_[i] >> dFs_[i];
  infile.close();
  return 0;
}

double s_minus(const double x){
  double mu = 1.0;//scale
  if (x < 1e-7 || x > 1.0) return 0;
  return (xf->xfxQ(3, x, mu) - xf->xfxQ(-3, x, mu)) / x;
}

double model_integrand(const double * var, const double * par){
  double alpha = par[0];
  double As = par[1];
  double x = var[0];
  double Q2 = par[2];
  if (x < 1e-7 || x > 1.0) return 0;
  double f_plus = 2.0 * alpha * (1.0 - x) * log(1.0 / x) + 2.0 * As * x * pow(1.0 - x, 2);
  return s_minus(x) * exp(-0.5 * Q2 * f_plus);
}

double model(const double Q2, const double * par){
  TF1 f("model_integrand", &model_integrand, 0.0, 1.0, 3);
  ROOT::Math::WrappedTF1 wf(f);
  double para[3] = {par[0], par[1], Q2};
  wf.SetParameters(para);
  ROOT::Math::IntegratorOneDim ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, 1e-5, 1e-5, 1000000);
  ig.SetFunction(wf);
  double result = ig.Integral(1e-7, 1.0);
  return result;
}

double chi2(const double * par){
  double sum = 0.0;
  for (int i = 0; i < 10; i++){
    sum += pow(model(Q2_[i], par) - Fs_[i], 2) / dFs_[i];
  }
  return sum;
}

double Minimize(const double * init, double * par){
  ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  min->SetMaxFunctionCalls(100000);
  min->SetTolerance(1e-4);
  min->SetPrintLevel(0);
  ROOT::Math::Functor f(&chi2, 2);
  min->SetFunction(f);
  min->SetLimitedVariable(0, "alpha", init[0], 1e-4, 0.7, 1.0);
  min->SetLimitedVariable(1, "As", init[1], 1e-4, 0.0, 5.0);
  min->Minimize();
  const double * xs = min->X();
  par[0] = xs[0];
  par[1] = xs[1];
  return min->MinValue();
}

double Minimize2(const double * init, double * par){
  ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  min->SetMaxFunctionCalls(100000);
  min->SetTolerance(1e-4);
  min->SetPrintLevel(0);
  ROOT::Math::Functor f(&chi2, 2);
  min->SetFunction(f);
  min->SetFixedVariable(0, "alpha", init[0]);
  min->SetVariable(1, "As", init[1], 1e-4);
  min->Minimize();
  const double * xs = min->X();
  par[0] = xs[0];
  par[1] = xs[1];
  return min->MinValue();
}

int main(){

  LoadData();
  double init[2] = {0.85, 2.0};
  double par[2], chi2min;

  FILE * file = fopen("results/parameters2.dat", "w");
  fprintf(file, "i\tchi2\talpha\tAs\n");
  for (int i = 0; i <= 100; i++){
    //if (i != 19) continue;
    xf = LHAPDF::mkPDF("NNPDF30_nlo_as_0118", i);
    chi2min = Minimize2(init, par);
    fprintf(file, "%d\t%.3E\t%.3E\t%.3E\n", i, chi2min, par[0], par[1]);
  }

  fclose(file);
  
  return 0;
}
  
