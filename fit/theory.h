#ifndef _THEORY_H_
#define _THEORY_H_

#include <iostream>
#include <fstream>
#include <cmath>

#include "Math/GSLIntegrator.h"
#include "Math/Integrator.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
#include "Math/WrappedTF1.h"
#include "Math/WrappedParamFunction.h"
#include "TF1.h"
#include "Math/Minimizer.h"
#include "TString.h"

#include "LHAPDF/LHAPDF.h"

using namespace std;

LHAPDF::PDF * xpdf;

double alpha, As, Asbar, Bs, Bsbar;

double xs(const double x, const double Q2){
  if (x > 1 || x < 1e-8) return 0;
  return xpdf->xfxQ2(3, x, Q2);
}

double xsbar(const double x, const double Q2){
  if (x > 1 || x < 1e-8) return 0;
  return xpdf->xfxQ2(-3, x, Q2);
}

double fs(const double x){
  return (1.0 - x) * alpha * log(1.0 / x) + As * x * pow(1.0 - x, 2) + Bs * (1.0 - x);
}

double fsbar(const double x){
  return (1.0 - x) * alpha * log(1.0 / x) + Asbar * x * pow(1.0 - x, 2) + Bsbar * (1.0 - x);
}

double xHs(const double x, const double t, const double Q2){
  return xs(x, Q2) * exp(t * fs(x));
}

double xHsbar(const double x, const double t, const double Q2){
  return xsbar(x, Q2) * exp(t * fsbar(x));
}

double F_integrand(const double * var, const double * par){
  double x = var[0];
  double t = par[0];
  double Q2 = par[1];
  double result = xHs(x, t, Q2) - xHsbar(x, t, Q2);
  return result / x;
}

double Fs(const double * par){
  //par: t, Q2
  TF1 f("integrand", &F_integrand, 0.0, 1.0, 2);
  ROOT::Math::WrappedTF1 wf(f);
  wf.SetParameters(par);
  ROOT::Math::IntegratorOneDim ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, 1e-6, 1e-5, 500);
  ig.SetFunction(wf);
  double result = ig.Integral(1e-6, 1e-4) + ig.Integral(1e-4, 1e-2) + ig.Integral(1e-2, 1.0);
  return result;
}


#endif
