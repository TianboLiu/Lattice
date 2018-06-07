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

using namespace std;

const double mL = 1.115683;
const double mK = 0.493677;
const double mq = 0.330;
const double ms = 0.480;
const double mD = 0.600;

double (* FFs)(const double Q2, const double * par);
double (* fssbar)(const double x, const double * par);

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

double Chi2(const double * par){
  double sum = 0;
  double theory = 0;
  for (int i = 0; i < 10; i++){
    theory = FFs(Variables[i], par);
    sum += pow(theory - Values[i], 2) / pow(Errors[i], 2);
  }
  return sum;
}

TMatrixD Minimize(const double * init, double * cent){
  ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  min->SetMaxFunctionCalls(100000);
  min->SetTolerance(1e-4);
  min->SetPrintLevel(0);
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

double Fx_1(const double * var, const double * par){
  //par: Q2, N, k, m1, m2
  double x = var[0];
  double Q2 = par[0];
  double k = par[2];
  double m1 = par[3];
  double m2 = par[4];
  double mu2 = m1 * m1 / x + m2 * m2 / (1.0 - x);
  double result = k * k / (4.0 * M_PI * M_PI) * x * (1.0 - x) * exp(- mu2 / (4.0 * k * k)) * exp(- (1.0 - x) / (16.0 * x * k * k) * Q2);
  return result;
}

double F_1(const double * par){
  //par: Q2, N, k, m1, m2
  TF1 f("Fx_1", &Fx_1, 0.0, 1.0, 5);
  ROOT::Math::WrappedTF1 wf(f);
  wf.SetParameters(par);
  ROOT::Math::IntegratorOneDim ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, 0.0, 1e-5, 1000);
  ig.SetFunction(wf);
  double result = ig.Integral(1e-6, 1.0 - 1e-6);
  return result;
}

double FFs_1(const double Q2, const double * par){
  //par: N, k
  double BM[5] = {Q2, par[0], par[1], mL, mK};
  double BM0[5] = {0.0, par[0], par[1], mL, mK};
  double MB[5] = {Q2, par[0], par[1], mK, mL};
  double MB0[5] = {0.0, par[0], par[1], mK, mL};
  double sD[5] = {Q2, par[0], par[1], ms, mD};
  double sD0[5] = {0.0, par[0], par[1], ms, mD};
  double sbarq[5] = {Q2, par[0], par[1], ms, mq};
  double sbarq0[5] = {0.0, par[0], par[1], ms, mq};
  double Fs = F_1(BM) / F_1(BM0) * F_1(sD) / F_1(sD0);
  double Fsbar = F_1(MB) / F_1(MB0) * F_1(sbarq) / F_1(sbarq0);
  double result = par[0] * (Fs - Fsbar);
  return result;
}

double Fx_2(const double * var, const double * par){
  //par: Q2, N, k, m1, m2
  double x = var[0];
  double Q2 = par[0];
  double k = par[2];
  double m1 = par[3];
  double m2 = par[4];
  double mu2 = m1 * m1 / x + m2 * m2 / (1.0 - x);
  double result = k * k / (16.0 * M_PI * M_PI) * exp(-mu2 / (k * k)) * exp(- (1.0 - x) / (4.0 * x * k * k) * Q2);
  return result;
}

double F_2(const double * par){
  //par: Q2, N, k, m1, m2
  TF1 f("Fx_2", &Fx_2, 0.0, 1.0, 5);
  ROOT::Math::WrappedTF1 wf(f);
  wf.SetParameters(par);
  ROOT::Math::IntegratorOneDim ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, 0.0, 1e-5, 1000);
  ig.SetFunction(wf);
  double result = ig.Integral(1e-6, 1.0 - 1e-6);
  return result;
}

double FFs_2(const double Q2, const double * par){
  //par: N, k
  double BM[5] = {Q2, par[0], par[1], mL, mK};
  double BM0[5] = {0.0, par[0], par[1], mL, mK};
  double MB[5] = {Q2, par[0], par[1], mK, mL};
  double MB0[5] = {0.0, par[0], par[1], mK, mL};
  double sD[5] = {Q2, par[0], par[1], ms, mD};
  double sD0[5] = {0.0, par[0], par[1], ms, mD};
  double sbarq[5] = {Q2, par[0], par[1], ms, mq};
  double sbarq0[5] = {0.0, par[0], par[1], ms, mq};
  double Fs = F_2(BM) / F_2(BM0) * F_2(sD) / F_2(sD0);
  double Fsbar = F_2(MB) / F_2(MB0) * F_2(sbarq) / F_2(sbarq0);
  double result = par[0] * (Fs - Fsbar);
  return result;
}

double Fx_3(const double * var, const double * par){
  //par: Q2, tau, k, m1, m2
  double x = var[0];
  double Q2 = par[0];
  double tau = par[1];
  double k = par[2];
  double m1 = par[3];
  double m2 = par[4];
  double mu2 = m1 * m1 / x + m2 * m2 / (1.0 - x);
  double result;
  if (tau == 2.0){
    result = exp(- x * log(1.0 / x) / (k * k * (1.0 - x)) * mu2) * exp(- log(1.0 / x) / (4.0 * k * k) * Q2);
  }
  else {
    result = pow(1.0 - x, tau - 2) * exp(- x * log(1.0 / x) / (k * k * (1.0 - x)) * mu2) * exp(- log(1.0 / x) / (4.0 * k * k) * Q2);
  }
  return result;
}

double F_3(const double * par){
  //par: Q2, tau, k, m1, m2
  TF1 f("Fx_3", &Fx_3, 0.0, 1.0, 5);
  ROOT::Math::WrappedTF1 wf(f);
  wf.SetParameters(par);
  ROOT::Math::IntegratorOneDim ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, 0.0, 1e-5, 1000);
  ig.SetFunction(wf);
  double result = ig.Integral(1e-6, 1.0 - 1e-6);
  return result;
}

double FFs_3(const double Q2, const double * par){
  //par: N, k
  double BM[5] = {Q2, 5.0, par[1], mL, mK};
  double BM0[5] = {0.0, 5.0, par[1], mL, mK};
  double MB[5] = {Q2, 5.0, par[1], mK, mL};
  double MB0[5] = {0.0, 5.0, par[1], mK, mL};
  double sD[5] = {Q2, 3.0, par[1], ms, mD};
  double sD0[5] = {0.0, 3.0, par[1], ms, mD};
  double sbarq[5] = {Q2, 2.0, par[1], ms, mq};
  double sbarq0[5] = {0.0, 2.0, par[1], ms, mq};
  double Fs = F_3(BM) / F_3(BM0) * F_3(sD) / F_3(sD0);
  double Fsbar = F_3(MB) / F_3(MB0) * F_3(sbarq) / F_3(sbarq0);
  double result = par[0] * (Fs - Fsbar);
  return result;
}

int Calculate_FFs(TMatrixD Cov, double * cent){
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
  for (int i = 0; i < 100; i++){
    Q2 = 0.01 * i;
    F = FFs(Q2, pars);
    dF = sqrt(pow(FFs(Q2, pars1p) - FFs(Q2, pars1m), 2) + pow(FFs(Q2, pars2p) - FFs(Q2, pars2m), 2)) / 2.0;
    fprintf(file, "%.4E\t%.4E\t%.4E\n", Q2, F, dF);
  }
  fclose(file);
  return 0;
}

double fs_1_integrand(const double * y, const double * par){
  //par: N, k
  double BM[5] = {0.0, par[0], par[1], mL, mK};
  double sD[5] = {0.0, par[0], par[1], ms, mD};
  double xq = par[2] / y[0];
  return Fx_1(y, BM) / F_1(BM) * Fx_1(&xq, sD) / F_1(sD) / y[0];
}

double fsbar_1_integrand(const double * y, const double * par){
  //par: N, k
  double MB[5] = {0.0, par[0], par[1], mK, mL};
  double sbarq[5] = {0.0, par[0], par[1], ms, mq};
  double xq = par[2] / y[0];
  return Fx_1(y, MB) / F_1(MB) * Fx_1(&xq, sbarq) / F_1(sbarq) / y[0];
}

double fssbar_1(const double x, const double * pars){
  //par: N, k
  double par[3] = {pars[0], pars[1], x};
  TF1 _fs("fs_1", &fs_1_integrand, 0.0, 1.0, 5);
  ROOT::Math::WrappedTF1 _wfs(_fs);
  _wfs.SetParameters(par);
  TF1 _fsbar("fsbar_1", &fsbar_1_integrand, 0.0, 1.0, 5);
  ROOT::Math::WrappedTF1 _wfsbar(_fsbar);
  _wfsbar.SetParameters(par);
  ROOT::Math::IntegratorOneDim ig(ROOT::Math::IntegrationOneDim::kADAPTIVE,0.0, 1e-5, 1000);
  ig.SetFunction(_wfs);
  double s = ig.Integral(x + 1e-12, 1.0 - 1e-12);
  ig.SetFunction(_wfsbar);
  double sb = ig.Integral(x + 1e-12, 1.0 - 1e-12);
  return (s - sb) * pars[0];
}

double fs_2_integrand(const double * y, const double * par){
  //par: N, k
  double BM[5] = {0.0, par[0], par[1], mL, mK};
  double sD[5] = {0.0, par[0], par[1], ms, mD};
  double xq = par[2] / y[0];
  return Fx_2(y, BM) / F_2(BM) * Fx_2(&xq, sD) / F_2(sD) / y[0];
}

double fsbar_2_integrand(const double * y, const double * par){
  //par: N, k
  double MB[5] = {0.0, par[0], par[1], mK, mL};
  double sbarq[5] = {0.0, par[0], par[1], ms, mq};
  double xq = par[2] / y[0];
  return Fx_2(y, MB) / F_2(MB) * Fx_2(&xq, sbarq) / F_2(sbarq) / y[0];
}

double fssbar_2(const double x, const double * pars){
  //par: N, k
  double par[3] = {pars[0], pars[1], x};
  TF1 _fs("fs_2", &fs_2_integrand, 0.0, 1.0, 5);
  ROOT::Math::WrappedTF1 _wfs(_fs);
  _wfs.SetParameters(par);
  TF1 _fsbar("fsbar_2", &fsbar_2_integrand, 0.0, 1.0, 5);
  ROOT::Math::WrappedTF1 _wfsbar(_fsbar);
  _wfsbar.SetParameters(par);
  ROOT::Math::IntegratorOneDim ig(ROOT::Math::IntegrationOneDim::kADAPTIVE,0.0, 1e-5, 1000);
  ig.SetFunction(_wfs);
  double s = ig.Integral(x + 1e-12, 1.0 - 1e-12);
  ig.SetFunction(_wfsbar);
  double sb = ig.Integral(x + 1e-12, 1.0 - 1e-12);
  return (s - sb) * pars[0];
}

double fs_3_integrand(const double * y, const double * par){
  //par: N, k
  double BM[5] = {0.0, 5.0, par[1], mL, mK};
  double sD[5] = {0.0, 3.0, par[1], ms, mD};
  double xq = par[2] / y[0];
  return Fx_3(y, BM) / F_3(BM) * Fx_3(&xq, sD) / F_3(sD) / y[0];
}

double fsbar_3_integrand(const double * y, const double * par){
  //par: N, k
  double MB[5] = {0.0, 5.0, par[1], mK, mL};
  double sbarq[5] = {0.0, 2.0, par[1], ms, mq};
  double xq = par[2] / y[0];
  return Fx_3(y, MB) / F_3(MB) * Fx_3(&xq, sbarq) / F_3(sbarq) / y[0];
}

double fssbar_3(const double x, const double * pars){
  //par: N, k
  double par[3] = {pars[0], pars[1], x};
  TF1 _fs("fs_3", &fs_3_integrand, 0.0, 1.0, 5);
  ROOT::Math::WrappedTF1 _wfs(_fs);
  _wfs.SetParameters(par);
  TF1 _fsbar("fsbar_3", &fsbar_3_integrand, 0.0, 1.0, 5);
  ROOT::Math::WrappedTF1 _wfsbar(_fsbar);
  _wfsbar.SetParameters(par);
  ROOT::Math::IntegratorOneDim ig(ROOT::Math::IntegrationOneDim::kADAPTIVE,0.0, 1e-5, 1000);
  ig.SetFunction(_wfs);
  double s = ig.Integral(x + 1e-12, 1.0 - 1e-12);
  ig.SetFunction(_wfsbar);
  double sb = ig.Integral(x + 1e-12, 1.0 - 1e-12);
  return (s - sb) * pars[0];
}

int Calculate_fssbar(TMatrixD Cov, double * cent){
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
    F = fssbar(x, pars);
    dF = sqrt(pow(fssbar(x, pars1p) - fssbar(x, pars1m), 2) + pow(fssbar(x, pars2p) - fssbar(x, pars2m), 2)) / 2.0;
    fprintf(file, "%.4E\t%.4E\t%.4E\n", x, F, dF);
  }
  fclose(file);
  return 0;
}

int Calculate_moment(TMatrixD Cov, double * cent){
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
  double x;
  double X[1000];
  double F = 0;
  double F1p = 0;
  double F1m = 0;
  double F2p = 0;
  double F2m = 0;
  for (int i = 0; i < 1000; i++){
    X[i] = 0.5 / 1000 + i / 1000.0;
  }
  for (int i = 0; i < 1000; i++){
    x = X[i];
    //cout << x << endl;
    F += x * fssbar(x, pars) / 1000;
    F1p += x * fssbar(x, pars1p) / 1000;
    F1m += x * fssbar(x, pars1m) / 1000;
    F2p += x * fssbar(x, pars2p) / 1000;
    F2m += x * fssbar(x, pars2m) / 1000;
  }
  cout << "<xS->: " << F << "  +-  " << sqrt(pow(F1p - F1m, 2) + pow(F2p - F2m, 2)) / 2.0 << endl;
  return 0;
}
  

///////

int main(const int argc, const char * argv[]){

  if (argc < 2){
    cout << "./model_fits <opt>" << endl;
    return 0;
  }
  
  LoadData();

  const int opt = atoi(argv[1]);

  if (opt == 1){
    FFs = & FFs_1;
    fssbar = & fssbar_1;
  }
  else if (opt == 2){
    FFs = & FFs_2;
    fssbar = & fssbar_2;
  }
  else if (opt == 3){
    FFs = & FFs_3;
    fssbar = & fssbar_3;
  }

  double init[2] = {0.013, 0.35};
  double cent[2] = {0, 0};

  TMatrixD Cov = Minimize(init, cent);

  Calculate_FFs(Cov, cent);
  //Calculate_fssbar(Cov, cent);
  Calculate_moment(Cov, cent);

  Cov.Print();

  return 0;
}
