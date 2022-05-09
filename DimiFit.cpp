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
  for(unsigned uBin=0; uBin<data2.GetNbins(); uBin++){
    TrueChi2 += pow((data2.GetBinContent(uBin)-TrueVal*2.)/(data2.GetBinError(uBin)),2.);
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

  fit1.SetParLimits(0,-10,10);
  fit1.SetParOrdMag(0,-2,1);
  fit1.SetNumBestSols(128);

  fit1.DEBUG_PrepareForWalk();

  TFile fOutput(TString::Format("%s/DimiFit/basics_1.root",GetFemtoOutputFolder()),"recreate");
  TGraph gSols;
  gSols.SetName("gSols");
  gSols.SetMarkerSize(1.0);
  gSols.SetMarkerStyle(20);
  gSols.SetLineWidth(0);

  TGraph gSolsExp_Pos;
  gSolsExp_Pos.SetName("gSolsExp_Pos");
  gSolsExp_Pos.SetMarkerSize(1.0);
  gSolsExp_Pos.SetMarkerStyle(20);
  gSolsExp_Pos.SetMarkerColor(kGreen+1);
  gSolsExp_Pos.SetLineWidth(0);

  TGraph gSolsExp_Neg;
  gSolsExp_Neg.SetName("gSolsExp_Neg");
  gSolsExp_Neg.SetMarkerSize(1.0);
  gSolsExp_Neg.SetMarkerStyle(20);
  gSolsExp_Neg.SetMarkerColor(kRed+1);
  gSolsExp_Neg.SetLineWidth(0);

  int PosSol=0;
  int TotSol=0;
  for(unsigned uSol=0; uSol<fit1.DEBUG_GetSolution()->size(); uSol++){
    int sol_exp;
    float sol = fit1.DEBUG_GetSolution()->at(uSol).Par->at(0);
    frexp(sol,&sol_exp);
    gSols.SetPoint(uSol,sol,0);
    if(sol>=0) gSolsExp_Pos.SetPoint(uSol,sol_exp,0);
    else gSolsExp_Neg.SetPoint(uSol,sol_exp,0);
    if(sol>=0) PosSol++;
    TotSol++;
  }
  printf("%.0f%% of (+) sols\n",double(PosSol)/double(TotSol)*100.);
  gSols.Write();
  gSolsExp_Pos.Write();
  gSolsExp_Neg.Write();

  fit1.DEBUG_WanderAround();
  float best_sol;
  float best_chi2 = 1e37;
  for(unsigned uSol=0; uSol<fit1.DEBUG_GetSolution()->size(); uSol++){
    float sol = fit1.DEBUG_GetSolution()->at(uSol).Par->at(0);
    float chi2 = fit1.DEBUG_GetSolution()->at(uSol).Chi2;
    if(chi2<best_chi2){
      best_chi2 = chi2;
      best_sol = sol;
    }
    //printf("%u: %.4f chi2=%.2f\n",uSol,sol,chi2);
  }
  printf("best_sol: %.4f chi2=%.2f\n",best_sol,best_chi2);


  //printf("we should quit\n");
}

//Gauss
void basics_2_model(const std::vector<float>& pars, std::vector<DLM_Histo<float>*>& model){
  const double PI = 3.14159265359;
  double xVal;
  double gauss_val;
//printf("pars size %i\n",pars.size());
  for(unsigned uMod=0; uMod<model.size(); uMod++){
    for(unsigned uBin=0; uBin<model.at(uMod)->GetNbins(); uBin++){
      xVal = model.at(uMod)->GetBinCenter(0,uBin);
      gauss_val = pars.at(0)*1./pars.at(2)/sqrt(2.*PI)*exp(-0.5*pow((xVal-pars.at(1))/pars.at(2),2.));
      model.at(uMod)->SetBinContent(uBin,gauss_val);
      model.at(uMod)->SetBinError(uBin,0);
      //printf("bin %u (%e) : %e\n",uBin,xVal,gauss_val);
    }
  }
}

void basics_2(){
  const unsigned NBins = 256;
  const float MinX = -10;
  const float MaxX = 10;
  const unsigned NumIter = 10*1000;

  DLM_Histo<float> data1;
  data1.SetUp(1);
  data1.SetUp(0,NBins,MinX,MaxX);
  data1.Initialize();

  TRandom3 RanGen(11);

  for(unsigned uIter=0; uIter<NumIter; uIter++){
    data1.Fill(RanGen.Gaus(0,1.5));
  }
  data1.ComputeError();

  TFile fOutput(TString::Format("%s/DimiFit/basics_2.root",GetFemtoOutputFolder()),"recreate");
  TH1F* hdata1 = Convert_DlmHisto_TH1F(&data1,"hdata1");
  TF1* fdata1 = new TF1("fdata1","[0]*TMath::Gaus(x,[1],[2],true)",MinX,MaxX);
  fdata1->SetParameter(0,NumIter);
  fdata1->SetParameter(1,0);
  fdata1->SetParameter(2,1);
  fdata1->SetParLimits(2,0,10);
  hdata1->Fit(fdata1,"S, N, R, M");

  TF1* fdimi1 = new TF1("fdimi1","[0]*TMath::Gaus(x,[1],[2],true)",MinX,MaxX);

  DLM_Fit fit1;
  fit1.SetUp(1,3);
  fit1.SetData(0,data1);
  fit1.SetFitFnct(basics_2_model);

  std::vector<float> pars;
  pars.push_back(fdata1->GetParameter(0));
  pars.push_back(fdata1->GetParameter(1));
  pars.push_back(fdata1->GetParameter(2));
  std::vector<DLM_Histo<float>*> result = fit1.Eval(pars);
  for(unsigned uBin=0; uBin<result.at(0)->GetNbins(); uBin++){
    //printf("%u -> %.2f\n",uBin,result.at(0)->GetBinContent(uBin));
  }
  printf("chi2 = %.5f / %u\n",fit1.Chi2(),fit1.Npts());

//return;
  //fit1.SetParLimits(0,pars.at(0)*0.9,pars.at(0)*1.1);
  fit1.SetParLimits(0,0,NumIter*10.);
  fit1.SetParOrdMag(0,0,37);
  fit1.SetParLimits(1,-1,1);
  fit1.SetParOrdMag(1,-5,0);
  //fit1.SetParLimits(1,0,0);
  //fit1.SetParLimits(2,pars.at(2)*0.9,pars.at(2)*1.1);
  fit1.SetParLimits(2,0.1,10);
  fit1.SetParOrdMag(2,-1,1);
  //fit1.SetNumBestSols(32);
  //fit1.SetNumWildCards(4);
  fit1.DEBUG_PrepareForWalk();
  fit1.DEBUG_WanderAround();
  //fit1.DEBUG_WanderAround();
  float best_sol[3];
  float best_chi2 = 1e37;
  for(unsigned uSol=0; uSol<fit1.DEBUG_GetSolution()->size(); uSol++){
    float chi2 = fit1.DEBUG_GetSolution()->at(uSol).Chi2;
    if(chi2<best_chi2){
      best_chi2 = chi2;
      best_sol[0] = fit1.DEBUG_GetSolution()->at(uSol).Par->at(0);
      best_sol[1] = fit1.DEBUG_GetSolution()->at(uSol).Par->at(1);
      best_sol[2] = fit1.DEBUG_GetSolution()->at(uSol).Par->at(2);
    }
    //printf("%u: %.4f chi2=%.2f\n",uSol,sol,chi2);
  }
  fdimi1->FixParameter(0,best_sol[0]);
  fdimi1->FixParameter(1,best_sol[1]);
  fdimi1->FixParameter(2,best_sol[2]);
  fdimi1->SetLineColor(kGreen+1);
  printf("best_sol: chi2=%.3f\n",best_chi2);
  printf("   %.2e vs %.2e\n",best_sol[0],pars.at(0));
  printf("   %.2e vs %.2e\n",best_sol[1],pars.at(1));
  printf("   %.2e vs %.2e\n",best_sol[2],pars.at(2));

  for(int i=0; i<256; i++){
    if(fit1.DEBUG_PrepareForWalk()) break;
    fit1.DEBUG_WanderAround();
  }

  printf("all done\n");
  std::vector<DLM_FitSolution>* FitSol = fit1.DEBUG_GetSolution();
  //for(int i=0; i<FitSol->size(); i++)
  //  printf("%f vs %f\n",(*FitSol)[i].Chi2,FitSol->at(i).Chi2);

  best_chi2 = 1e37;
  for(unsigned uSol=0; uSol<fit1.DEBUG_GetSolution()->size(); uSol++){
    float chi2 = fit1.DEBUG_GetSolution()->at(uSol).Chi2;
    if(chi2<best_chi2){
      best_chi2 = chi2;
      best_sol[0] = fit1.DEBUG_GetSolution()->at(uSol).Par->at(0);
      best_sol[1] = fit1.DEBUG_GetSolution()->at(uSol).Par->at(1);
      best_sol[2] = fit1.DEBUG_GetSolution()->at(uSol).Par->at(2);
    }
    //printf("%u: %.4f chi2=%.2f\n",uSol,sol,chi2);
  }
  printf("best_sol: chi2=%.3f\n",best_chi2);
  printf("   %.2e vs %.2e\n",best_sol[0],pars.at(0));
  printf("   %.2e vs %.2e\n",best_sol[1],pars.at(1));
  printf("   %.2e vs %.2e\n",best_sol[2],pars.at(2));


  fOutput.cd();
  hdata1->Write();
  fdata1->Write();
  fdimi1->Write();
  delete hdata1;
}



//Gauss (sum of 3)
void basics_3_model(const std::vector<float>& pars, std::vector<DLM_Histo<float>*>& model){
  const double PI = 3.14159265359;
  double xVal;
  double gauss_val;
//printf("pars size %i\n",pars.size());
  for(unsigned uMod=0; uMod<model.size(); uMod++){
    for(unsigned uBin=0; uBin<model.at(uMod)->GetNbins(); uBin++){
      xVal = model.at(uMod)->GetBinCenter(0,uBin);
      gauss_val  = pars.at(0)*1./pars.at(2)/sqrt(2.*PI)*exp(-0.5*pow((xVal-pars.at(1))/pars.at(2),2.));
      gauss_val += pars.at(3)*1./pars.at(5)/sqrt(2.*PI)*exp(-0.5*pow((xVal-pars.at(4))/pars.at(5),2.));
      gauss_val += pars.at(6)*1./pars.at(8)/sqrt(2.*PI)*exp(-0.5*pow((xVal-pars.at(7))/pars.at(8),2.));
      model.at(uMod)->SetBinContent(uBin,gauss_val);
      model.at(uMod)->SetBinError(uBin,0);
      //printf("bin %u (%e) : %e\n",uBin,xVal,gauss_val);
    }
  }
}

void basics_3(){
  const unsigned NBins = 256;
  const float MinX = -10;
  const float MaxX = 10;
  const unsigned NumIter = 10*1000;

  DLM_Histo<float> data1;
  data1.SetUp(1);
  data1.SetUp(0,NBins,MinX,MaxX);
  data1.Initialize();

  TRandom3 RanGen(11);

  for(unsigned uIter=0; uIter<NumIter; uIter++){
    data1.Fill(RanGen.Gaus(0,1.0));
    if(RanGen.Uniform()<0.35){
      data1.Fill(RanGen.Gaus(3.5,1.0));
      data1.Fill(RanGen.Gaus(-3.5,1.0));
    }
  }
  data1.ComputeError();

  TFile fOutput(TString::Format("%s/DimiFit/basics_3.root",GetFemtoOutputFolder()),"recreate");
  TH1F* hdata1 = Convert_DlmHisto_TH1F(&data1,"hdata1");
  TF1* fdata1 = new TF1("fdata1","[0]*TMath::Gaus(x,[1],[2],true)+[3]*TMath::Gaus(x,[4],[5],true)+[6]*TMath::Gaus(x,[7],[8],true)",MinX,MaxX);
  fdata1->SetParameter(0,NumIter);
  fdata1->SetParameter(1,0);
  fdata1->SetParameter(2,1);
  fdata1->SetParLimits(2,0,10);

  fdata1->SetParameter(3,NumIter);
  fdata1->SetParameter(4,0);
  fdata1->SetParameter(5,1);
  fdata1->SetParLimits(5,0,10);

  fdata1->SetParameter(6,NumIter);
  fdata1->SetParameter(7,0);
  fdata1->SetParameter(8,1);
  fdata1->SetParLimits(8,0,10);

  hdata1->Fit(fdata1,"S, N, R, M");


  TF1* fdimi1 = new TF1("fdimi1","[0]*TMath::Gaus(x,[1],[2],true)+[3]*TMath::Gaus(x,[4],[5],true)+[6]*TMath::Gaus(x,[7],[8],true)",MinX,MaxX);


  DLM_Fit fit1;
  fit1.SetUp(1,9);
  fit1.SetData(0,data1);
  fit1.SetFitFnct(basics_3_model);

  std::vector<float> pars;
  pars.push_back(fdata1->GetParameter(0));
  pars.push_back(fdata1->GetParameter(1));
  pars.push_back(fdata1->GetParameter(2));
  pars.push_back(fdata1->GetParameter(3));
  pars.push_back(fdata1->GetParameter(4));
  pars.push_back(fdata1->GetParameter(5));
  pars.push_back(fdata1->GetParameter(6));
  pars.push_back(fdata1->GetParameter(7));
  pars.push_back(fdata1->GetParameter(8));
  std::vector<DLM_Histo<float>*> result = fit1.Eval(pars);
  for(unsigned uBin=0; uBin<result.at(0)->GetNbins(); uBin++){
    //printf("%u -> %.2f\n",uBin,result.at(0)->GetBinContent(uBin));
  }
  printf("chi2 = %.5f / %u\n",fit1.Chi2(),fit1.Npts());

  for(int i=0; i<9; i++){
    fit1.SetParLimits(i,fdata1->GetParameter(i)*0.9,fdata1->GetParameter(i)*1.1);
  }


/*
  fit1.SetParLimits(0,fdata1->GetParameter(0)*0.99,fdata1->GetParameter(0)*1.01);
  fit1.SetParOrdMag(0,0,37);
  fit1.SetParLimits(1,-0.1,0.1);
  fit1.SetParOrdMag(1,-5,0);
  fit1.SetParLimits(2,0.5,1.5);
  fit1.SetParOrdMag(2,-1,1);

  fit1.SetParLimits(3,fdata1->GetParameter(3)*0.99,fdata1->GetParameter(3)*1.01);
  fit1.SetParOrdMag(3,0,37);
  fit1.SetParLimits(4,-4,-3);
  fit1.SetParOrdMag(4,-5,0);
  fit1.SetParLimits(5,0.5,1.5);
  fit1.SetParOrdMag(5,-1,1);

  fit1.SetParLimits(6,fdata1->GetParameter(6)*0.99,fdata1->GetParameter(6)*1.01);
  fit1.SetParOrdMag(6,0,37);
  fit1.SetParLimits(7,3,4);
  fit1.SetParOrdMag(7,-5,0);
  fit1.SetParLimits(8,0.5,1.5);
  fit1.SetParOrdMag(8,-1,1);
*/
  fit1.SetNumBestSols(32);
  fit1.SetNumWildCards(4);
  fit1.DEBUG_PrepareForWalk();
  fit1.DEBUG_WanderAround();
  //fit1.DEBUG_WanderAround();
  float best_sol[3];
  float best_chi2 = 1e37;
  for(unsigned uSol=0; uSol<fit1.DEBUG_GetSolution()->size(); uSol++){
    float chi2 = fit1.DEBUG_GetSolution()->at(uSol).Chi2;
    if(chi2<best_chi2){
      best_chi2 = chi2;
      best_sol[0] = fit1.DEBUG_GetSolution()->at(uSol).Par->at(0);
      best_sol[1] = fit1.DEBUG_GetSolution()->at(uSol).Par->at(1);
      best_sol[2] = fit1.DEBUG_GetSolution()->at(uSol).Par->at(2);
      best_sol[3] = fit1.DEBUG_GetSolution()->at(uSol).Par->at(3);
      best_sol[4] = fit1.DEBUG_GetSolution()->at(uSol).Par->at(4);
      best_sol[5] = fit1.DEBUG_GetSolution()->at(uSol).Par->at(5);
      best_sol[6] = fit1.DEBUG_GetSolution()->at(uSol).Par->at(6);
      best_sol[7] = fit1.DEBUG_GetSolution()->at(uSol).Par->at(7);
      best_sol[8] = fit1.DEBUG_GetSolution()->at(uSol).Par->at(8);
    }
    //printf("%u: %.4f chi2=%.2f\n",uSol,sol,chi2);
  }
  fdimi1->FixParameter(0,best_sol[0]);
  fdimi1->FixParameter(1,best_sol[1]);
  fdimi1->FixParameter(2,best_sol[2]);
  fdimi1->FixParameter(3,best_sol[3]);
  fdimi1->FixParameter(4,best_sol[4]);
  fdimi1->FixParameter(5,best_sol[5]);
  fdimi1->FixParameter(6,best_sol[6]);
  fdimi1->FixParameter(7,best_sol[7]);
  fdimi1->FixParameter(8,best_sol[8]);
  fdimi1->SetLineColor(kGreen+1);
  printf("best_sol: chi2=%.3f\n",best_chi2);
  printf("   %.2e vs %.2e\n",best_sol[0],pars.at(0));
  printf("   %.2e vs %.2e\n",best_sol[1],pars.at(1));
  printf("   %.2e vs %.2e\n",best_sol[2],pars.at(2));

  for(int i=0; i<64; i++){
    if(fit1.DEBUG_PrepareForWalk()) break;
    fit1.DEBUG_WanderAround();
  }

  printf("all done\n");
  std::vector<DLM_FitSolution>* FitSol = fit1.DEBUG_GetSolution();
  //for(int i=0; i<FitSol->size(); i++)
  //  printf("%f vs %f\n",(*FitSol)[i].Chi2,FitSol->at(i).Chi2);

  best_chi2 = 1e37;
  for(unsigned uSol=0; uSol<fit1.DEBUG_GetSolution()->size(); uSol++){
    float chi2 = fit1.DEBUG_GetSolution()->at(uSol).Chi2;
    if(chi2<best_chi2){
      best_chi2 = chi2;
      best_sol[0] = fit1.DEBUG_GetSolution()->at(uSol).Par->at(0);
      best_sol[1] = fit1.DEBUG_GetSolution()->at(uSol).Par->at(1);
      best_sol[2] = fit1.DEBUG_GetSolution()->at(uSol).Par->at(2);
    }
    //printf("%u: %.4f chi2=%.2f\n",uSol,sol,chi2);
  }
  printf("best_sol: chi2=%.3f\n",best_chi2);
  printf("   %.2e vs %.2e\n",best_sol[0],pars.at(0));
  printf("   %.2e vs %.2e\n",best_sol[1],pars.at(1));
  printf("   %.2e vs %.2e\n",best_sol[2],pars.at(2));


  fOutput.cd();
  hdata1->Write();
  fdata1->Write();
  fdimi1->Write();
  delete hdata1;
}


int DIMI_FIT(int narg, char** ARGS){
  //basics_1();
  //basics_2();
  basics_3();
  printf("ciao DIMI_FIT\n");
  return 0;
}
