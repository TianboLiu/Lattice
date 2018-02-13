#include "theory.h"

double Variables[10], Values[10], Errors[10];

int LoadData(){
  char tmp[200];
  ifstream file("f1smod.txt");
  file.getline(tmp, 200);
  for (int i = 0; i < 10; i++){
    file >> Variables[i] >> Values[i] >> Errors[i];
    cout << Variables[i] << "  " << Values[i] << "  " << Errors[i] << endl;
  }
  file.close();
  return 0;
}

double mu2 = 4.0;
double Chi2(const double * parameters){
  alpha = parameters[0];
  As = parameters[1];
  Asbar = parameters[2];
  Bs = parameters[3];
  Bsbar = parameters[4];
  double par[2];
  par[1] = mu2;
  double sum = 0.0;
  for (int i = 0; i < 10; i++){
    par[0] = - Variables[i];
    sum += pow(Fs(par) - Values[i], 2) / pow(Errors[i], 2);
  }
  return sum;
}

double Minimize0(const double * init){
  ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  min->SetMaxFunctionCalls(100000);
  min->SetTolerance(1e-4);
  min->SetPrintLevel(1);
  ROOT::Math::Functor f(&Chi2, 5);
  min->SetFunction(f);
  min->SetVariable(0, "alpha", init[0], 1e-4);
  //min->SetLimitedVariable(0, "alpha", init[0], 1e-4, 0.85, 0.95);
  min->SetVariable(1, "As", init[1], 1e-4);
  min->SetVariable(2, "Asbar", init[2], 1e-4);
  min->SetVariable(3, "Bs", init[3], 1e-4);
  min->SetVariable(4, "Bsbar", init[4], 1e-4);
  min->Minimize();
  const double * xs = min->X();
  alpha = xs[0];
  As = xs[1];
  Asbar = xs[2];
  Bs = xs[3];
  Bsbar = xs[4];
  return min->MinValue();
}

double Minimize1(const double * init){
  ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  min->SetMaxFunctionCalls(100000);
  min->SetTolerance(1e-3);
  min->SetPrintLevel(0);
  ROOT::Math::Functor f(&Chi2, 5);
  min->SetFunction(f);
  //min->SetVariable(0, "alpha", init[0], 1e-4);
  min->SetLimitedVariable(0, "alpha", init[0], 1e-4, 0.85, 0.95);
  min->SetLowerLimitedVariable(1, "As", init[1], 1e-4, 0.0);
  min->SetLowerLimitedVariable(2, "Asbar", init[2], 1e-4, 0.0);
  //min->SetLimitedVariable(2, "Asbar", init[2], 1e-4, 0.0, 10.0);
  min->SetFixedVariable(3, "Bs", init[3]);
  min->SetFixedVariable(4, "Bsbar", init[4]);
  min->Minimize();
  const double * xs = min->X();
  alpha = xs[0];
  As = xs[1];
  Asbar = xs[2];
  Bs = xs[3];
  Bsbar = xs[4];
  return min->MinValue();
} 

int main(const int argc, const char * argv[]){

  if (argc < 2){
    cout << "./fit <opt>" << endl;
    return 0;
  }

  const int opt = atoi(argv[1]);

  if (opt == 0){//test
    xpdf = LHAPDF::mkPDF("NNPDF30_nlo_as_0118", 0);
    LoadData();
    double init[5] = {0.90, 1.0, 1.0, 0.0, 0.0};
    cout << Minimize1(init) << endl;  
    cout << alpha << " " << As << " " << Asbar << " " << Bs << " " << Bsbar << endl;
  }

  if (opt == 1){//MMHT
    if (argc < 3){
      cout << "./fit <opt=1> <mu2>" << endl;
      return 0;
    }	      
    mu2 = atof(argv[2]);
    FILE * fs = fopen("mmht1.dat", "w");
    fprintf(fs, "MMHT2014 NLO Q2 = %.1f GeV^2\n", mu2);
    fprintf(fs, "id\tchi2\talpha\tAs\tAsbar\tBs\tBsbar\n");
    LoadData();
    double chi2;
    double init[5] = {0.90, 1.0, 1.0, 0.0, 0.0};
    for (int i = 0; i <= 50; i++){
      cout << ">>> Running curve:  " << i << endl;
      xpdf = LHAPDF::mkPDF("MMHT2014nlo68cl", i);
      chi2 = Minimize1(init);
      cout << chi2 << endl;
      fprintf(fs, "%d\t%.2f\t%.2E\t%.2E\t%.2E\t%.2E\t%.2E\n",
	      i, chi2, alpha, As, Asbar, Bs, Bsbar);
    }
    fclose(fs);
  }

  if (opt == 2){//NNPDF
    if (argc < 3){
      cout << "./fit <opt=2> <mu2>" << endl;
      return 0;
    }	      
    mu2 = atof(argv[2]);
    FILE * fs = fopen("nnpdf1.dat", "w");
    fprintf(fs, "NNPDF3.0 NLO Q2 = %.1f GeV^2\n", mu2);
    fprintf(fs, "id\tchi2\talpha\tAs\tAsbar\tBs\tBsbar\n");
    LoadData();
    double chi2;
    double init[5] = {0.90, 1.0, 1.0, 0.0, 0.0};
    for (int i = 0; i <= 100; i++){
      cout << ">>> Running curve:  " << i << endl;
      xpdf = LHAPDF::mkPDF("NNPDF30_nlo_as_0118", i);
      chi2 = Minimize1(init);
      cout << chi2 << endl;
      fprintf(fs, "%d\t%.2f\t%.2E\t%.2E\t%.2E\t%.2E\t%.2E\n",
	      i, chi2, alpha, As, Asbar, Bs, Bsbar);
    }
    fclose(fs);
  }

  


  return 0;
}
  
