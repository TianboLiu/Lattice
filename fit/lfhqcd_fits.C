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
#include "TMatrixDEigen.h"
#include "Math/SpecFuncMathCore.h"

using namespace std;

double Variables[10], Values[10], Errors[10];
int LoadData(){
  ifstream file("f1smod.txt");
  file.ignore(300, '\n');
  for (int i = 0; i < 10; i++){
    file >> Variables[i] >> Values[i] >> Errors[i];
    cout << Variables[i] << endl;
  }
  file.close();
  return 0;
}

double Beta(const double x, const double y){
  return ROOT::Math::beta(x, y);
}

double Gamma(const double x){
  return ROOT::Math::tgamma(x);
}

double Ftau(const double tau, const double Q2, const double lambda){
  double N = sqrt(M_PI) * Gamma(tau - 1.0) / Gamma(tau - 0.5);
  return Beta(tau - 1.0, 0.5 + Q2 / (4.0 * lambda)) / N;
}
	    
double F1s(const double Q2, const double * par){
  double A = par[0];
  double lambda = pow(par[1], 2);
  return A * (Ftau(5, Q2, lambda) - Ftau(6, Q2, lambda));
}

double wx(const double x){
  double a = 0.531;
  return pow(x, 1.0 - x) * exp(-a * pow(1.0 - x, 2));
}

double dwx(const double x){
  double a = 0.531;
  return 2.0 * a * (1.0 - x) * pow(x, 1.0 - x) * exp(-a * pow(1.0 - x, 2)) + pow(x, 1.0 - x) * exp(-a * pow(1.0 - x, 2)) * ((1.0 - x) / x - log(x));
}

double qxtau(const double tau, const double x){
  double N = sqrt(M_PI) * Gamma(tau - 1.0) / Gamma(tau - 0.5);
  return pow(1.0 - wx(x), tau - 2.0) * pow(wx(x), -0.5) * dwx(x) / N;
}

double qs(const double x, const double * par){
  double A = par[0];
  return A * (qxtau(5.0, x) - qxtau(6.0, x));
}

double Chi2(const double * par){
  double sum = 0;
  double theory = 0;
  for (int i = 0; i < 10; i++){
    theory = F1s(Variables[i], par);
    sum += pow(theory - Values[i], 2) / pow(Errors[i], 2);
  }
  return sum;
}

TMatrixD Minimize(const double * init, double * cent){
  ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  min->SetMaxFunctionCalls(100000);
  min->SetTolerance(1e-4);
  min->SetPrintLevel(2);
  ROOT::Math::Functor f(&Chi2, 2);
  min->SetFunction(f);
  min->SetVariable(0, "N", init[0], 1e-4);
  min->SetVariable(1, "k", init[1], 1e-4);
  min->Minimize();
  min->PrintResults();
  double cov[4];
  min->GetCovMatrix(cov);
  TMatrixD Cov(2,2);
  Cov(0,0) = cov[0];
  Cov(0,1) = cov[1];
  Cov(1,0) = cov[2];
  Cov(1,1) = cov[3];
  const double * xs = min->X();
  cent[0] = xs[0];
  cent[1] = xs[1];
  return Cov;
}

int Calculate_F1s(TMatrixD Cov, double * cent){
  TMatrixDEigen eigen(Cov);
  TMatrixD val = eigen.GetEigenValues();
  TMatrixD vec = eigen.GetEigenVectors();
  double pars[2] = {cent[0], cent[1]};
  double pars1p[2] = {cent[0] + sqrt(val(0,0)) * vec(0,0),
		      cent[1] + sqrt(val(0,0)) * vec(1,0)};
  double pars1m[2] = {cent[0] - sqrt(val(0,0)) * vec(0,0),
		      cent[1] - sqrt(val(0,0)) * vec(1,0)};
  double pars2p[2] = {cent[0] + sqrt(val(1,1)) * vec(0,1),
		      cent[1] + sqrt(val(1,1)) * vec(1,1)};
  double pars2m[2] = {cent[0] - sqrt(val(1,1)) * vec(0,1),
		      cent[1] - sqrt(val(1,1)) * vec(1,1)};
  double Q2, F, dF;
  FILE * file = fopen("fs.dat", "w");
  fprintf(file, "Q2\tFs\tdFs\n");
  for (int i = 0; i < 1000; i++){
    Q2 = 0.01 * i;
    F = F1s(Q2, pars);
    dF = sqrt(pow(F1s(Q2, pars1p) - F1s(Q2, pars1m), 2) + pow(F1s(Q2, pars2p) - F1s(Q2, pars2m), 2)) / 2.0;
    fprintf(file, "%.4E\t%.4E\t%.4E\n", Q2, F, dF);
  }
  fclose(file);
  return 0;
}

int Single_curve(){
  double pars[2] = {0.0127, 0.5482};
  double Q2, F;
  FILE * file = fopen("fs.dat", "w");
  fprintf(file, "Q2\tFs\n");
  for (int i = 0; i < 1000; i++){
    Q2 = 0.01 * i;
    F = F1s(Q2, pars);
    fprintf(file, "%.4E\t%.4E\n", Q2, F);
  }
  fclose(file);
  return 0;
}

int Calculate_qs(TMatrixD Cov, double * cent){
  TMatrixDEigen eigen(Cov);
  TMatrixD val = eigen.GetEigenValues();
  TMatrixD vec = eigen.GetEigenVectors();
  double pars[2] = {cent[0], cent[1]};
  double pars1p[2] = {cent[0] + sqrt(val(0,0)) * vec(0,0),
		      cent[1] + sqrt(val(0,0)) * vec(1,0)};
  double pars1m[2] = {cent[0] - sqrt(val(0,0)) * vec(0,0),
		      cent[1] - sqrt(val(0,0)) * vec(1,0)};
  double pars2p[2] = {cent[0] + sqrt(val(1,1)) * vec(0,1),
		      cent[1] + sqrt(val(1,1)) * vec(1,1)};
  double pars2m[2] = {cent[0] - sqrt(val(1,1)) * vec(0,1),
		      cent[1] - sqrt(val(1,1)) * vec(1,1)};
  double x, F, dF;
  FILE * file = fopen("fss.dat", "w");
  fprintf(file, "x\ts-sbar\tdelta\n");
  double X[1000];
  for (int i = 0; i < 500; i++){
    X[i] = pow(10.0, -4.0 + 0.004 * i);
    X[i+500] = 0.01 + (1.0 - 1e-10 - 0.01) / 499 * i;
  }
  for (int i = 0; i < 1000; i++){
    x = X[i];
    //cout << x << endl;
    F = qs(x, pars);
    dF = sqrt(pow(qs(x, pars1p) - qs(x, pars1m), 2) + pow(qs(x, pars2p) - qs(x, pars2m), 2)) / 2.0;
    fprintf(file, "%.4E\t%.4E\t%.4E\n", x, F, dF);
  }
  fclose(file);
  return 0;
}


int main(){

  LoadData();

  double init[2] = {0.013, 0.5};
  double cent[2] = {0, 0};

  //cout << Chi2(init) << endl;

  TMatrixD Cov = Minimize(init, cent);

  Cov.Print();

  Calculate_F1s(Cov, cent);
  //Calculate_qs(Cov, cent);

  //Single_curve();

  //cout << Gamma(1) << endl;
  return 0;
}
