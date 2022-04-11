
#include "Deuteron.h"

#include "EnvVars.h"
#include "CATS.h"
#include "DLM_Source.h"
#include "DLM_Potentials.h"
#include "CommonAnaFunctions.h"
#include "CATSconstants.h"
#include "DLM_Histo.h"

#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TF1.h"


//s-wave only
void p_pn_cumulant(){

  //estimated, based on pN = 1/2 pd, that the mT of the pN (N is p or n from the deuteron)
  //that for mT(pd) = 1.65 we have mT(pN) of c.a. 1.3, so an effective radius of c.a. 1.25
  const double SourceSize = 1.25;
  const double lambda = 0.8/8.;
  const double NumMomBins = 300;
  const unsigned Reb = 20;
  const double kMin = 0;
  const double kMax = 600;

  CATSparameters cSource(CATSparameters::tSource,2,true);
  cSource.SetParameter(0,SourceSize);
  cSource.SetParameter(1,2);

  CATSparameters cPot_pp_1s0(CATSparameters::tPotential,8,true);
  double pp_pars_1s0[8]={NN_AV18,v18_Coupled3P2,1,1,1,0,0,0};
  cPot_pp_1s0.SetParameters(pp_pars_1s0);

  CATSparameters cPot_pn_1s0(CATSparameters::tPotential,8,true);
  //#,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
  double pn_pars_1s0[8]={NN_AV18,v18_Coupled3P2,1,1,-1,0,0,0};
  cPot_pn_1s0.SetParameters(pn_pars_1s0);

  CATSparameters cPot_pn_3s1(CATSparameters::tPotential,8,true);
  double pn_pars_3s1[8]={NN_AV18,v18_Coupled3P2,1,1,-1,1,0,1};
  cPot_pn_3s1.SetParameters(pn_pars_3s1);

  CATS Cat_pp;
  Cat_pp.SetMomBins(NumMomBins,kMin,kMax);
  Cat_pp.SetThetaDependentSource(false);
  Cat_pp.SetAnaSource(GaussSource, cSource);
  Cat_pp.SetUseAnalyticSource(true);
  Cat_pp.SetMomentumDependentSource(false);
  Cat_pp.SetExcludeFailedBins(false);
  Cat_pp.SetQ1Q2(1);
  Cat_pp.SetQuantumStatistics(true);
  Cat_pp.SetRedMass( 0.5*Mass_p );
  Cat_pp.SetNumChannels(2);
  Cat_pp.SetNumPW(0,1);
  Cat_pp.SetNumPW(1,0);
  Cat_pp.SetSpin(0,0);
  Cat_pp.SetSpin(1,1);
  Cat_pp.SetChannelWeight(0, 1./4.);
  Cat_pp.SetChannelWeight(1, 3./4.);
  Cat_pp.SetShortRangePotential(0,0,fDlmPot,cPot_pp_1s0);
  Cat_pp.SetEpsilonConv(10e-8);
  Cat_pp.SetEpsilonProp(10e-8);
  Cat_pp.KillTheCat();

  CATS Cat_pn;
  Cat_pn.SetMomBins(NumMomBins,kMin,kMax);
  Cat_pn.SetThetaDependentSource(false);
  Cat_pn.SetAnaSource(GaussSource, cSource);
  Cat_pn.SetUseAnalyticSource(true);
  Cat_pn.SetMomentumDependentSource(false);
  Cat_pn.SetExcludeFailedBins(false);
  Cat_pn.SetQ1Q2(0);
  //Cat_pn.SetPdgId(2212, 2112);
  Cat_pn.SetQuantumStatistics(false);
  Cat_pn.SetRedMass( (Mass_p*Mass_n)/(Mass_p+Mass_n) );
  Cat_pn.SetNumChannels(2);
  Cat_pn.SetNumPW(0,1);
  Cat_pn.SetNumPW(1,1);
  Cat_pn.SetSpin(0,0);
  Cat_pn.SetSpin(1,1);
  Cat_pn.SetChannelWeight(0, 1./4.);
  Cat_pn.SetChannelWeight(1, 3./4.);
  Cat_pn.SetShortRangePotential(0,0,fDlmPot,cPot_pn_1s0);
  Cat_pn.SetShortRangePotential(1,0,fDlmPot,cPot_pn_3s1);
  Cat_pn.SetEpsilonConv(10e-8);
  Cat_pn.SetEpsilonProp(10e-8);
  Cat_pn.KillTheCat();


  TFile fInputME(TString::Format("%s/Deuteron/pd_files/UnfoldedME.root",GetFemtoOutputFolder()),"read");
  TDirectoryFile* dirResults=(TDirectoryFile*)(fInputME.FindObjectAny("GraphHistoVar_0"));
  printf("dirResults = %p\n",dirResults);
  TH1F* hMeUnf=NULL;
  dirResults->GetObject("hTotalME _Rebinned_10_Reweighted",hMeUnf);
  printf("hMeUnf = %p\n",hMeUnf);
  hMeUnf->Sumw2();
  hMeUnf->Scale(1./hMeUnf->Integral(),"width");
  TF1* fit_Me = new TF1("fit_Me",GaussSourceTF1,0,600,1);
  fit_Me->SetParameter(0,400);
  hMeUnf->Fit(fit_Me,"S, N, R, M");

  TFile fOutput(TString::Format("%s/Deuteron/p_pn_cumulant.root",GetFemtoOutputFolder()),"recreate");
  hMeUnf->Write();
  fit_Me->Write();

  TGraph graph_pp;
  graph_pp.SetName("graph_pp");
  graph_pp.SetLineWidth(5);
  graph_pp.SetLineColor(kBlue);
  graph_pp.SetLineStyle(7);

  TGraph graph_pn;
  graph_pn.SetName("graph_pn");
  graph_pn.SetLineWidth(5);
  graph_pn.SetLineColor(kGreen+2);
  graph_pn.SetLineStyle(2);

  TGraph graph_pp_x_pn;
  graph_pp_x_pn.SetName("graph_pp_x_pn");
  graph_pp_x_pn.SetLineWidth(5);
  graph_pp_x_pn.SetLineColor(kBlack);

  TH1F* weighted_pp_x_pn = new TH1F("weighted_pp_x_pn","weighted_pp_x_pn",NumMomBins/Reb*2.,kMin,kMax*2.);
  weighted_pp_x_pn->SetLineWidth(5);
  weighted_pp_x_pn->SetLineColor(kAzure+2);

  TH1F* weighted2_pp_x_pn = new TH1F("weighted2_pp_x_pn","weighted2_pp_x_pn",NumMomBins/Reb*2.,kMin,kMax*2.);
  weighted2_pp_x_pn->SetLineWidth(5);
  weighted2_pp_x_pn->SetLineColor(kAzure+2);

  for(unsigned uBin=0; uBin<NumMomBins/Reb*2.; uBin++){
    weighted_pp_x_pn->SetBinContent(uBin+1,1);
    weighted2_pp_x_pn->SetBinContent(uBin+1,1);
  }

  double kstar,c_pp,c_pn;
  int wBin=0;
  double CkInt=0;
  double CkNorm=0;
  for(unsigned uBin=0; uBin<NumMomBins; uBin++){
    kstar = Cat_pp.GetMomentum(uBin);

    c_pp = lambda*Cat_pp.EvalCorrFun(kstar*0.75)+1.-lambda;
    c_pn = lambda*Cat_pn.EvalCorrFun(kstar*0.75)+1.-lambda;
    graph_pp.SetPoint(uBin,kstar,c_pp);
    graph_pn.SetPoint(uBin,kstar,c_pn);
    graph_pp_x_pn.SetPoint(uBin,kstar,c_pp*c_pn);

    CkInt += (c_pp*c_pn)*fit_Me->Eval(kstar);
    CkNorm += fit_Me->Eval(kstar);
    //the last bin
    if((uBin+1)%Reb==0){
      double BinCenter = 20;
      BinCenter += (40.)*double(wBin);
      if(wBin==0) BinCenter = 30;
      if(BinCenter<kMax/0.75){
        weighted_pp_x_pn->SetBinContent(wBin+1,CkInt/CkNorm);
        weighted2_pp_x_pn->SetBinContent(wBin+1,graph_pp_x_pn.Eval(BinCenter));
      }
      else{
        weighted_pp_x_pn->SetBinContent(wBin+1,1);
        weighted2_pp_x_pn->SetBinContent(wBin+1,1);
      }

      CkInt=0;
      CkNorm=0;
      wBin++;
    }
  }

  graph_pp.Write();
  graph_pn.Write();
  graph_pp_x_pn.Write();
  weighted_pp_x_pn->Write();
  weighted2_pp_x_pn->Write();

  delete weighted_pp_x_pn;
}


int DEUTERON_MAIN(int argc, char *argv[]){
  p_pn_cumulant();
  return 0;
}
