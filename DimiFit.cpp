#include "DimiFit.h"
#include "CommonAnaFunctions.h"
#include "EnvVars.h"
#include "DLM_CppTools.h"
#include "DLM_OmpTools.h"
#include "DLM_Random.h"
#include "DLM_RootWrapper.h"
#include "DLM_HistoAnalysis.h"
#include "DLM_Source.h"
#include "DLM_Potentials.h"
#include "DLM_Histo.h"
#include "DLM_MathFunctions.h"

#include <iostream>
#include <unistd.h>
#include <omp.h>

#include "TREPNI.h"
#include "CATS.h"
#include "CATSconstants.h"
#include "CATStools.h"
#include "CECA.h"
#include "DLM_Fit.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TF1.h"
#include "TGenPhaseSpace.h"
#include "TRandom3.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TNtuple.h"


//constant function
void basics_1_model(const std::vector<float>& pars, std::vector<DLM_Histo<float>*>& model){
  for(unsigned uMod=0; uMod<model.size(); uMod++){
    for(unsigned uBin=0; uBin<model.at(uMod)->GetNbins(); uBin++){
      model.at(uMod)->SetBinContent(uBin,pars.at(0)*(1.+uMod));
      model.at(uMod)->SetBinError(uBin,0);
    }
  }
}

void basics_1(){

  DLM_Histo<float> data1;
  data1.SetUp(1);
  data1.SetUp(0,4,0,4);
  data1.Initialize();

  DLM_Histo<float> data2;
  data2.SetUp(1);
  data2.SetUp(0,2,0,2);
  data2.Initialize();

  const float TrueVal = 1.01;
  data1.SetBinContent(unsigned(0),1.1);
  data1.SetBinError(unsigned(0),0.1);
  data1.SetBinContent(unsigned(1),0.9);
  data1.SetBinError(unsigned(1),0.1);
  data1.SetBinContent(unsigned(2),1.0);
  data1.SetBinError(unsigned(2),0.1);
  data1.SetBinContent(unsigned(3),1.05);
  data1.SetBinError(unsigned(3),0.1);

  data2.SetBinContent(unsigned(0),2.12);
  data2.SetBinError(unsigned(0),0.1);
  data2.SetBinContent(unsigned(1),2.02);
  data2.SetBinError(unsigned(1),0.1);

  float TrueChi2 = 0;
  for(unsigned uBin=0; uBin<data1.GetNbins(); uBin++){
    TrueChi2 += pow((data1.GetBinContent(uBin)-TrueVal)/(data1.GetBinError(uBin)),2.);
  }


  DLM_Fit fit1;

  fit1.SetUp(2,1);
  fit1.SetData(0,data1);

  fit1.SetData(1,data2);
  fit1.SetFitFnct(basics_1_model);

  std::vector<float> pars;
  pars.push_back(TrueVal);

  std::vector<DLM_Histo<float>*> result = fit1.Eval(pars);
  for(unsigned uBin=0; uBin<result.at(0)->GetNbins(); uBin++){
    printf("%u -> %.2f\n",uBin,result.at(0)->GetBinContent(uBin));
  }
  printf("chi2 = %.5f vs %.5f / %u\n",fit1.Chi2(),TrueChi2,fit1.Npts());



}

int DIMI_FIT(int narg, char** ARGS){
  basics_1();
  printf("ciao DIMI_FIT\n");
  return 0;
}
