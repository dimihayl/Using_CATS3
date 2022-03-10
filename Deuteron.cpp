
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

//s-wave only
void p_pn_cumulant(){

  const double SourceSize = 1.2;
  const double lambda = 0.8;
  const double NumMomBins = 100;
  const double kMin = 0;
  const double kMax = 400;

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
  Cat_pp.SetPdgId(2212, 2212);
  Cat_pp.SetRedMass( 0.5*Mass_p );
  Cat_pp.SetNumChannels(2);
  Cat_pp.SetNumPW(0,1);
  Cat_pp.SetNumPW(1,0);
  Cat_pp.SetSpin(0,0);
  Cat_pp.SetSpin(1,1);
  Cat_pp.SetChannelWeight(0, 1./4.);
  Cat_pp.SetChannelWeight(1, 3./4.);
  Cat_pp.SetShortRangePotential(0,0,fDlmPot,cPot_pp_1s0);
  Cat_pp.KillTheCat();


  CATS Cat_pn;
  Cat_pn.SetMomBins(NumMomBins,kMin,kMax);
  Cat_pn.SetThetaDependentSource(false);
  Cat_pn.SetAnaSource(GaussSource, cSource);
  Cat_pn.SetUseAnalyticSource(true);
  Cat_pn.SetMomentumDependentSource(false);
  Cat_pn.SetExcludeFailedBins(false);
  Cat_pn.SetQ1Q2(0);
  Cat_pn.SetPdgId(2212, 2112);
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
  Cat_pn.KillTheCat();

  TFile fOutput(TString::Format("%s/Deuteron/p_pn_cumulant.root",GetFemtoOutputFolder()),"recreate");

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

  double kstar,c_pp,c_pn;
  for(unsigned uBin=0; uBin<NumMomBins; uBin++){
    kstar = Cat_pp.GetMomentum(uBin);

    c_pp = lambda*Cat_pp.EvalCorrFun(kstar/0.75)+1.-lambda;
    c_pn = lambda*Cat_pn.EvalCorrFun(kstar/0.75)+1.-lambda;
    graph_pp.SetPoint(uBin,kstar,c_pp);
    graph_pn.SetPoint(uBin,kstar,c_pn);
    graph_pp_x_pn.SetPoint(uBin,kstar,c_pp*c_pn);
  }

  graph_pp.Write();
  graph_pn.Write();
  graph_pp_x_pn.Write();

}


int DEUTERON_MAIN(int argc, char *argv[]){
  p_pn_cumulant();
  return 0;
}
