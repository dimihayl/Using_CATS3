
#include <fstream>

#include "Deuteron.h"

#include "EnvVars.h"
#include "CATS.h"
#include "CATStools.h"
#include "DLM_Source.h"
#include "DLM_Potentials.h"
#include "CommonAnaFunctions.h"
#include "CATSconstants.h"
#include "DLM_Histo.h"
#include "DLM_RootWrapper.h"
#include "DLM_MathFunctions.h"
#include "DLM_Ck.h"
#include "DLM_CkDecomposition.h"
#include "DLM_HistoAnalysis.h"

#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"

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



void SetUpKaonDeuteronFSI(CATS& cat, TString FSI_type){

}

void SetUpKaonDeuteronSource(CATS& cat, TString Source_type){

}


//takes the output of SetUp_Kdp_Kd, to produce the correlation functions
//in addition, it creates some basic correlations to compare to, see below (some extra files are set up there !)
// /that is 9% (so 91% the genuine)
// pol1 to fit
void MyOwn_Kd_v1(){

  //1.04 original, 0.96 updated
  const double SourceSize = 1.04;

  //const double lam_gen = 0.91;//real values
  const double lam_gen = 1.0;
  const double avg_mt = 1500;

  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());

  double kMin = 0;
  double kMax = 400;
  unsigned kSteps = 20;

  TGraph* gCk_ER_CECA_Tau = NULL;
/*
  gCk_ER_CECA_Tau = new TGraph [21];
//here we create the graphs for different taus
  DLM_Histo<KdpPars> dlm_kdp_3D;
  dlm_kdp_3D.QuickLoad(TString::Format("%s/dKaon/SetUp_Kdp_Kd/Full_v1/dlm_kdp_3D.dlm",GetCernBoxDimi()).Data());
  //printf("dim = %u\n",dlm_kdp_3D.GetDim());
  //for(int iDim=0; iDim<3; iDim++){
  //  printf(" D%i: %i; %.2f <> %.2f\n",iDim, dlm_kdp_3D.GetNbins(iDim), dlm_kdp_3D.GetLowEdge(iDim), dlm_kdp_3D.GetUpEdge(iDim));
  //}
  //dlm_kdp_3D.GetBinContent(int(0)).Print();
  //return;
  DLM_MtKstar_KdpSource source_kdp(&dlm_kdp_3D);
  CATS Kd_ER_CECA_Tau;
  Kd_ER_CECA_Tau.SetMomBins(kSteps, kMin, kMax);
  AnalysisObject.SetUpCats_Kd(Kd_ER_CECA_Tau,"DG_ER","",0,0);//NLO_Coupled_S
  Kd_ER_CECA_Tau.SetUseAnalyticSource(true);
  Kd_ER_CECA_Tau.SetMomentumDependentSource(true);
  Kd_ER_CECA_Tau.SetAnaSource(CatsSourceForwarder, &source_kdp, 2);
  Kd_ER_CECA_Tau.SetAnaSource(0, avg_mt);
  for(int iTau=0; iTau<21; iTau++){
    double Tau = double(iTau)*0.5;
    //if(Tau==0) Tau = 0.01;
    gCk_ER_CECA_Tau[iTau].SetName(TString::Format("gCk_ER_CECA_Tau_%.2f",Tau));
    gCk_ER_CECA_Tau[iTau].SetLineColor(kOrange+iTau);
    gCk_ER_CECA_Tau[iTau].SetLineWidth(3);

    Kd_ER_CECA_Tau.SetAnaSource(1, Tau);
    Kd_ER_CECA_Tau.KillTheCat();

    for(unsigned uPts=0; uPts<kSteps; uPts++){
      double kstar = Kd_ER_CECA_Tau.GetMomentum(uPts);
      gCk_ER_CECA_Tau[iTau].SetPoint(uPts,kstar,Kd_ER_CECA_Tau.GetCorrFun(uPts)*lam_gen+1-lam_gen);
    }
  }
*/



//here we define from which file we would like to take a CECA kstar integrated source
  //TFile fSrcInput(TString::Format("%s/FunWithCeca/Ceca_Kd_EffFix/old/CECA_Kd_Source.root",GetFemtoOutputFolder()),"read");
  //TFile fSrcInput(TString::Format("%s/dKaon/Source/Full/CECA_Kd_Source_t3.root",GetCernBoxDimi()),"read");
  TFile fSrcInput(TString::Format("%s/FunWithCeca/Ceca_KdReso_FAST/Eta0.8_KdReso_ET1_PR1_DD0.0_EF-1402.root",GetFemtoOutputFolder()),"read");

  TH1F* GhettoFemto_rstar  = (TH1F*)fSrcInput.Get("GhettoFemto_rstar");
  DLM_Histo<float>* dlmGhettoFemto_rstar = Convert_TH1F_DlmHisto(GhettoFemto_rstar);
  DLM_HistoSource dlmCecaSource(dlmGhettoFemto_rstar);

  CATS Kd_ER;
  Kd_ER.SetMomBins(kSteps, kMin, kMax);
  AnalysisObject.SetUpCats_Kd(Kd_ER,"DG_ER","McGauss_ResoTM",0,202);//NLO_Coupled_S
  Kd_ER.SetAnaSource(0, SourceSize);
  Kd_ER.KillTheCat();

  CATS Kd_ER_G;
  Kd_ER_G.SetMomBins(kSteps, kMin, kMax);
  AnalysisObject.SetUpCats_Kd(Kd_ER_G,"DG_ER","DoubleGauss",0,0);//NLO_Coupled_S
  //Kd_ER_G.SetAnaSource(0, 1.1);//old source
  //Kd_ER_G.SetAnaSource(1, 2.14);//old source
  Kd_ER_G.SetAnaSource(0, 1.04);//updated source
  Kd_ER_G.SetAnaSource(1, 2.01);//updated source
  Kd_ER_G.SetAnaSource(2, 0.76);
  Kd_ER_G.KillTheCat();

  CATS Kd_ER_CECA;
  Kd_ER_CECA.SetMomBins(kSteps, kMin, kMax);
  AnalysisObject.SetUpCats_Kd(Kd_ER_CECA,"DG_ER","",0,0);//NLO_Coupled_S
  Kd_ER_CECA.SetUseAnalyticSource(true);
  Kd_ER_CECA.SetAnaSource(CatsSourceForwarder, &dlmCecaSource, 0);
  Kd_ER_CECA.KillTheCat();

  CATS Kd_FCA;
  Kd_FCA.SetMomBins(kSteps, kMin, kMax);
  AnalysisObject.SetUpCats_Kd(Kd_FCA,"DG_FCA","McGauss_ResoTM",0,202);//NLO_Coupled_S
  Kd_FCA.SetAnaSource(0, SourceSize);
  Kd_FCA.KillTheCat();

  CATS Kd_FCA_G;
  Kd_FCA_G.SetMomBins(kSteps, kMin, kMax);
  AnalysisObject.SetUpCats_Kd(Kd_FCA_G,"DG_ER","DoubleGauss",0,0);//NLO_Coupled_S
  //Kd_FCA_G.SetAnaSource(0, 1.1);//old source
  //Kd_FCA_G.SetAnaSource(1, 2.14);//old source
  Kd_FCA_G.SetAnaSource(0, 1.04);//updated source
  Kd_FCA_G.SetAnaSource(1, 2.01);//updated source
  Kd_FCA_G.SetAnaSource(2, 0.76);
  Kd_FCA_G.KillTheCat();

  TGraph gCk_ER;
  gCk_ER.SetName("gCk_ER");
  gCk_ER.SetLineColor(kCyan);
  gCk_ER.SetLineWidth(6);

  TGraph gCk_ER_G;
  gCk_ER_G.SetName("gCk_ER_G");
  gCk_ER_G.SetLineColor(kCyan+1);
  gCk_ER_G.SetLineWidth(5);

  TGraph gCk_ER_CECA;
  gCk_ER_CECA.SetName("gCk_ER_CECA");
  gCk_ER_CECA.SetLineColor(kCyan+2);
  gCk_ER_CECA.SetLineWidth(5);

  TGraph gCk_FCA;
  gCk_FCA.SetName("gCk_FCA");
  gCk_FCA.SetLineColor(kAzure+1);
  gCk_FCA.SetLineWidth(4);

  TGraph gCk_FCA_G;
  gCk_FCA_G.SetName("gCk_FCA_G");
  gCk_FCA_G.SetLineColor(kAzure+2);
  gCk_FCA_G.SetLineWidth(5);

  double rMin = 0;
  double rMax = 64;
  unsigned rSteps = 4096;
  TH1F* hSrc_RSM = new TH1F("hSrc_RSM","hSrc_RSM",rSteps,rMin,rMax);
  TH1F* hSrc_DG = new TH1F("hSrc_DG","hSrc_DG",rSteps,rMin,rMax);

  for(unsigned uPts=0; uPts<kSteps; uPts++){
    double kstar = Kd_ER.GetMomentum(uPts);
    gCk_ER.SetPoint(uPts,kstar,Kd_ER.GetCorrFun(uPts)*lam_gen+1-lam_gen);
    gCk_ER_G.SetPoint(uPts,kstar,Kd_ER_G.GetCorrFun(uPts)*lam_gen+1-lam_gen);
    gCk_ER_CECA.SetPoint(uPts,kstar,Kd_ER_CECA.GetCorrFun(uPts)*lam_gen+1-lam_gen);
    gCk_FCA.SetPoint(uPts,kstar,Kd_FCA.GetCorrFun(uPts)*lam_gen+1-lam_gen);
    gCk_FCA_G.SetPoint(uPts,kstar,Kd_FCA_G.GetCorrFun(uPts)*lam_gen+1-lam_gen);
  }


  for(unsigned uRad=0; uRad<rSteps; uRad++){
    double rstar = hSrc_RSM->GetBinCenter(uRad+1);
    double src_rsm = Kd_ER.EvaluateTheSource(1,rstar,0);
    double src_dg = Kd_ER_G.EvaluateTheSource(1,rstar,0);

    hSrc_RSM->SetBinContent(uRad+1,src_rsm);
    hSrc_DG->SetBinContent(uRad+1,src_dg);
  }

  TFile fOutput(TString::Format("%s/Deuteron/MyOwn_Kd_v1/MyOwn_Kd_v1_Full.root",GetFemtoOutputFolder()),"recreate");
  gCk_ER.Write();
  gCk_ER_G.Write();
  gCk_ER_CECA.Write();
  gCk_FCA.Write();
  gCk_FCA_G.Write();


  for(int iTau=0; iTau<21; iTau++){
    if(gCk_ER_CECA_Tau)
      gCk_ER_CECA_Tau[iTau].Write();
  }

  hSrc_RSM->Write();
  hSrc_DG->Write();

  if(gCk_ER_CECA_Tau)
    delete [] gCk_ER_CECA_Tau;
  delete hSrc_RSM;
  delete hSrc_DG;
  delete dlmGhettoFemto_rstar;
}


//in frascatti, compare to LL, produce the needed output for Oton
void MyOwn_Kd_v2(){

  //1.04 original, 0.96 updated
  const double SourceSize = 1.04;

  //const double lam_gen = 0.91;//real values
  const double lam_gen = 1.0;
  const double avg_mt = 1500;

  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());

  double kMin = 0;
  double kMax = 400;
  unsigned kSteps = 100;

  double ScatLen, EffRan;
  TH1F* hFit;
  TF1* fitSP;

  CATS Kd_ER;
  Kd_ER.SetMomBins(kSteps, kMin, kMax);
  AnalysisObject.SetUpCats_Kd(Kd_ER,"DG_ER","Gauss",0,202);//NLO_Coupled_S
  Kd_ER.SetAnaSource(0, SourceSize);
  Kd_ER.KillTheCat();

  CATS Kd_ER_G;
  Kd_ER_G.SetMomBins(kSteps, kMin, kMax);
  AnalysisObject.SetUpCats_Kd(Kd_ER_G,"DG_ER","DoubleGauss",0,0);//NLO_Coupled_S
  //Kd_ER_G.SetAnaSource(0, 1.1);//old source
  //Kd_ER_G.SetAnaSource(1, 2.14);//old source
  Kd_ER_G.SetAnaSource(0, 1.04);//updated source
  Kd_ER_G.SetAnaSource(1, 2.01);//updated source
  Kd_ER_G.SetAnaSource(2, 0.76);
  Kd_ER_G.SetQ1Q2(0);
  Kd_ER_G.KillTheCat();
  GetScattParameters(Kd_ER_G, ScatLen, EffRan, hFit, fitSP, 3, false, false, 0);
  printf("ER DG Scattering Parameters: f0 = %.2f, d0 = %.2f\n", ScatLen, EffRan);  
  Kd_ER_G.SetQ1Q2(1);
  Kd_ER_G.KillTheCat();

  CATS Kd_ER_SG_1fm;
  Kd_ER_SG_1fm.SetMomBins(kSteps, kMin, kMax);
  AnalysisObject.SetUpCats_Kd(Kd_ER_SG_1fm,"DG_ER","Gauss",0,0);//NLO_Coupled_S
  Kd_ER_SG_1fm.SetAnaSource(0, 1);
  Kd_ER_SG_1fm.KillTheCat();

  CATS Kd_ER_SG_3fm;
  Kd_ER_SG_3fm.SetMomBins(kSteps, kMin, kMax);
  AnalysisObject.SetUpCats_Kd(Kd_ER_SG_3fm,"DG_ER","Gauss",0,0);//NLO_Coupled_S
  Kd_ER_SG_3fm.SetAnaSource(0, 3);
  Kd_ER_SG_3fm.KillTheCat();

  CATS Kd_FCA;
  Kd_FCA.SetMomBins(kSteps, kMin, kMax);
  AnalysisObject.SetUpCats_Kd(Kd_FCA,"DG_FCA","McGauss_ResoTM",0,202);//NLO_Coupled_S
  Kd_FCA.SetAnaSource(0, SourceSize);
  Kd_FCA.KillTheCat();

  CATS Kd_FCA_G;
  Kd_FCA_G.SetMomBins(kSteps, kMin, kMax);
  AnalysisObject.SetUpCats_Kd(Kd_FCA_G,"DG_FCA","DoubleGauss",0,0);//NLO_Coupled_S
  //Kd_FCA_G.SetAnaSource(0, 1.1);//old source
  //Kd_FCA_G.SetAnaSource(1, 2.14);//old source
  Kd_FCA_G.SetAnaSource(0, 1.04);//updated source
  Kd_FCA_G.SetAnaSource(1, 2.01);//updated source
  Kd_FCA_G.SetAnaSource(2, 0.76);
  Kd_FCA_G.KillTheCat();

  CATS Kd_FCA_SG;
  Kd_FCA_SG.SetMomBins(kSteps, kMin, kMax);
  AnalysisObject.SetUpCats_Kd(Kd_FCA_SG,"DG_FCA","Gauss",0,0);//NLO_Coupled_S
  Kd_FCA_SG.SetAnaSource(0, 3);
  Kd_FCA_SG.KillTheCat();

  //square well potential,see mail from Johann @ CERN on 27.11.2024
  CATS Kd_SW_G;
  Kd_SW_G.SetMomBins(kSteps, kMin, kMax);
  AnalysisObject.SetUpCats_Kd(Kd_SW_G,"SW_ER","DoubleGauss",0,0);//NLO_Coupled_S
  //Kd_ER_G.SetAnaSource(0, 1.1);//old source
  //Kd_ER_G.SetAnaSource(1, 2.14);//old source
  Kd_SW_G.SetAnaSource(0, 1.04);//updated source
  Kd_SW_G.SetAnaSource(1, 2.01);//updated source
  Kd_SW_G.SetAnaSource(2, 0.76);
  Kd_SW_G.SetQ1Q2(0);
  Kd_SW_G.KillTheCat();
  GetScattParameters(Kd_SW_G, ScatLen, EffRan, hFit, fitSP, 3, false, false, 0);
  printf("Square Well Scattering Parameters: f0 = %.2f, d0 = %.2f\n", ScatLen, EffRan);
  Kd_SW_G.SetQ1Q2(1);
  Kd_SW_G.KillTheCat();

  CATS Kd_SW_SG_1fm;
  Kd_SW_SG_1fm.SetMomBins(kSteps, kMin, kMax);
  AnalysisObject.SetUpCats_Kd(Kd_SW_SG_1fm,"SW_ER","Gauss",0,0);//NLO_Coupled_S
  Kd_SW_SG_1fm.SetAnaSource(0, 1);
  Kd_SW_SG_1fm.KillTheCat();

  CATS Kd_SW_SG_3fm;
  Kd_SW_SG_3fm.SetMomBins(kSteps, kMin, kMax);
  AnalysisObject.SetUpCats_Kd(Kd_SW_SG_3fm,"SW_ER","Gauss",0,0);//NLO_Coupled_S
  Kd_SW_SG_3fm.SetAnaSource(0, 3);
  Kd_SW_SG_3fm.KillTheCat();


  TGraph gCk_ER;
  gCk_ER.SetName("gCk_ER");
  gCk_ER.SetLineColor(kCyan);
  gCk_ER.SetLineWidth(6);

  TGraph gCk_ER_G;
  gCk_ER_G.SetName("gCk_ER_G");
  gCk_ER_G.SetLineColor(kCyan+1);
  gCk_ER_G.SetLineWidth(5);

  TGraph gCk_ER_SG_1fm;
  gCk_ER_SG_1fm.SetName("gCk_ER_SG_1fm");
  gCk_ER_SG_1fm.SetLineColor(kCyan+3);
  gCk_ER_SG_1fm.SetLineWidth(5);

  TGraph gCk_ER_SG_3fm;
  gCk_ER_SG_3fm.SetName("gCk_ER_SG_3fm");
  gCk_ER_SG_3fm.SetLineColor(kCyan-2);
  gCk_ER_SG_3fm.SetLineWidth(5);

  TGraph gCk_ER_CECA;
  gCk_ER_CECA.SetName("gCk_ER_CECA");
  gCk_ER_CECA.SetLineColor(kCyan+2);
  gCk_ER_CECA.SetLineWidth(5);

  TGraph gCk_SWER_G;
  gCk_SWER_G.SetName("gCk_SWER_G");
  gCk_SWER_G.SetLineColor(kCyan-1);
  gCk_SWER_G.SetLineWidth(5);


  TGraph gCk_SWER_SG_1fm;
  gCk_SWER_SG_1fm.SetName("gCk_SWER_SG_1fm");
  gCk_SWER_SG_1fm.SetLineColor(kCyan-3);
  gCk_SWER_SG_1fm.SetLineWidth(5);

  TGraph gCk_SWER_SG_3fm;
  gCk_SWER_SG_3fm.SetName("gCk_SWER_SG_3fm");
  gCk_SWER_SG_3fm.SetLineColor(kCyan-4);
  gCk_SWER_SG_3fm.SetLineWidth(5);

  TGraph gCk_FCA;
  gCk_FCA.SetName("gCk_FCA");
  gCk_FCA.SetLineColor(kAzure+1);
  gCk_FCA.SetLineWidth(4);

  TGraph gCk_FCA_G;
  gCk_FCA_G.SetName("gCk_FCA_G");
  gCk_FCA_G.SetLineColor(kAzure+2);
  gCk_FCA_G.SetLineWidth(5);

  TGraph gCk_FCA_SG;
  gCk_FCA_SG.SetName("gCk_FCA_SG");
  gCk_FCA_SG.SetLineColor(kAzure+3);
  gCk_FCA_SG.SetLineWidth(5);

  double rMin = 0;
  double rMax = 64;
  unsigned rSteps = 4096;
  TH1F* hSrc_RSM = new TH1F("hSrc_RSM","hSrc_RSM",rSteps,rMin,rMax);
  TH1F* hSrc_DG = new TH1F("hSrc_DG","hSrc_DG",rSteps,rMin,rMax);

  for(unsigned uPts=0; uPts<kSteps; uPts++){
    double kstar = Kd_ER.GetMomentum(uPts);
    gCk_ER.SetPoint(uPts,kstar,Kd_ER.GetCorrFun(uPts)*lam_gen+1-lam_gen);
    gCk_ER_G.SetPoint(uPts,kstar,Kd_ER_G.GetCorrFun(uPts)*lam_gen+1-lam_gen);
    gCk_ER_SG_1fm.SetPoint(uPts,kstar,Kd_ER_SG_1fm.GetCorrFun(uPts)*lam_gen+1-lam_gen);
    gCk_ER_SG_3fm.SetPoint(uPts,kstar,Kd_ER_SG_3fm.GetCorrFun(uPts)*lam_gen+1-lam_gen);
    gCk_FCA.SetPoint(uPts,kstar,Kd_FCA.GetCorrFun(uPts)*lam_gen+1-lam_gen);
    gCk_FCA_G.SetPoint(uPts,kstar,Kd_FCA_G.GetCorrFun(uPts)*lam_gen+1-lam_gen);
    gCk_FCA_SG.SetPoint(uPts,kstar,Kd_FCA_SG.GetCorrFun(uPts)*lam_gen+1-lam_gen);
    gCk_SWER_G.SetPoint(uPts,kstar,Kd_SW_G.GetCorrFun(uPts)*lam_gen+1-lam_gen);
    gCk_SWER_SG_1fm.SetPoint(uPts,kstar,Kd_SW_SG_1fm.GetCorrFun(uPts)*lam_gen+1-lam_gen);
    gCk_SWER_SG_3fm.SetPoint(uPts,kstar,Kd_SW_SG_3fm.GetCorrFun(uPts)*lam_gen+1-lam_gen);

  }


  for(unsigned uRad=0; uRad<rSteps; uRad++){
    double rstar = hSrc_RSM->GetBinCenter(uRad+1);
    double src_rsm = Kd_ER.EvaluateTheSource(1,rstar,0);
    double src_dg = Kd_ER_G.EvaluateTheSource(1,rstar,0);

    hSrc_RSM->SetBinContent(uRad+1,src_rsm);
    hSrc_DG->SetBinContent(uRad+1,src_dg);
  }

  TFile fOutput(TString::Format("%s/Deuteron/MyOwn_Kd_v2/MyOwn_Kd_v2_Full.root",GetFemtoOutputFolder()),"recreate");
  gCk_ER.Write();
  gCk_ER_G.Write();
  gCk_ER_SG_1fm.Write();
  gCk_ER_SG_3fm.Write();
  gCk_SWER_G.Write();
  gCk_SWER_SG_1fm.Write();
  gCk_SWER_SG_3fm.Write();
  gCk_FCA.Write();
  gCk_FCA_G.Write();
  gCk_FCA_SG.Write();

  hSrc_RSM->Write();
  hSrc_DG->Write();

  delete hSrc_RSM;
  delete hSrc_DG;
}

//takes the CECA output, and creates a KDP file with certain settings for the kstar and mT.
//the output will contain the full mt binnig so that we can better estimate the source for a specific mT, 
//but has a single kstar bin, that is integrated up to a cutoff value (100 MeV to mimic CECA paper)
//the input file is NOT the root, but the DLM_Histo
void MyOwn_Kd_v2_CreateKdp(double kstar_cutoff=100){

  TFile fSummary(TString::Format("%s/dKaon/Frascati/Sources/kdp_maps_full_v2.root",GetCernBoxDimi()), "recreate");

  std::vector<double> DD = {0,0.5,1.0,1.5,2.0,2.25,2.5,3.0,3.5,4.0,4.5,4.75,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0};
  //std::vector<double> DD = {0,3.0};
  std::vector<double> DD_bins;
  DD_bins.push_back( DD.at(0) - (DD.at(1)-DD.at(0))*0.5 );
  for(unsigned uEl=0; uEl<DD.size()-1; uEl++){
    DD_bins.push_back( (DD.at(uEl+1)+DD.at(uEl))*0.5 );
  }
  DD_bins.push_back( DD.at(DD.size()-1) + (DD.at(DD.size()-1)-DD.at(DD.size()-2))*0.5 );
  //for(unsigned uEl=0; uEl<DD.size()+1; uEl++){
  //  printf("%.3f ", DD_bins.at(uEl));
  //}
  //printf("\n");
//return;

  std::vector<float> KR = {47.16, 52.40, 57.64};
  std::vector<TString> SPAR = {"FIT","NLO19"};
  std::vector<TString> TYPE = {"Kd","KdReso"};

  for(TString& sTYPE : TYPE){
    printf("sTYPE = %s\n",sTYPE.Data());
    for(TString& sSPAR : SPAR){
      printf(" sSPAR = %s\n",sSPAR.Data());
      for(float& fKR : KR){
        printf("  fKR = %.3f\n",fKR);
        //this will be a 2D kdp map in mt and d_delay
        DLM_Histo<KdpPars>* dlm_kdp_map = NULL;
        for(double& fDD : DD){
          printf("   fDD = %.3f\n",fDD);
          TString InputFileName = TString::Format("%s/dKaon/Frascati/Sources/Raw/%s/Eta0.8_%s_PR35.78_KR%.2f_DD%.2f.Ghetto_kstar_rstar_mT",
            GetCernBoxDimi(), sSPAR.Data(), sTYPE.Data(), fKR, fDD);
          DLM_Histo<float> dlmSrc;
          dlmSrc.QuickLoad(InputFileName.Data());
          dlmSrc.ComputeError();

          DLM_Histo<float> dlmSrc_trimmed;
          dlmSrc_trimmed.SetUp(3);
          dlmSrc_trimmed.SetUp(0,1,0,kstar_cutoff);
          dlmSrc_trimmed.SetUp(1,dlmSrc,1);//rstar
          dlmSrc_trimmed.SetUp(2,dlmSrc,2);//mT
          dlmSrc_trimmed.Initialize();
          for(unsigned uKstar=0; uKstar<dlmSrc.GetNbins(0); uKstar++){
            double kstar = dlmSrc.GetBinCenter(0, uKstar);
            if(kstar>kstar_cutoff) break;
            for(unsigned uRstar=0; uRstar<dlmSrc.GetNbins(1); uRstar++){
              for(unsigned uMt=0; uMt<dlmSrc.GetNbins(2); uMt++){
                double current_value = dlmSrc_trimmed.GetBinContent(0, uRstar, uMt);
                double current_error = dlmSrc_trimmed.GetBinError(0, uRstar, uMt);
                dlmSrc_trimmed.SetBinContent(0, uRstar, uMt, current_value + dlmSrc.GetBinContent(uKstar, uRstar, uMt));
                dlmSrc_trimmed.SetBinError(0, uRstar, uMt, sqrt(current_error*current_error + pow(dlmSrc.GetBinError(uKstar, uRstar, uMt),2.)));
              }
            }
          }

          //takes a 3D histo of Kstar Rstar Mt and returns a 2D kdp histo of Mt Kstar (in that order)
          DLM_Histo<KdpPars>* dlm_kdp_2D = Convert_3Dsource_Kdp(dlmSrc_trimmed,true,3,96,3);
          //DLM_Histo<KdpPars>* dlm_kdp_2D_kstar = Convert_3Dsource_Kdp(dlmSrc,true,3,96,3);
          //we create the 2D kdp histo, including tau, at the first iteration
          if(dlm_kdp_map==NULL){
            dlm_kdp_map = new DLM_Histo<KdpPars>();
            dlm_kdp_map->SetUp(2);
            dlm_kdp_map->SetUp(0, *dlm_kdp_2D, 0);//mT
            double* a_DD_bins = &DD_bins[0];
            double* a_DD = &DD[0];
            dlm_kdp_map->SetUp(1, DD.size(), a_DD_bins, a_DD);
            dlm_kdp_map->Initialize();
          }
          //we fill up the whole big histo
          for(unsigned uMt=0; uMt<dlm_kdp_2D->GetNbins(0); uMt++){
            dlm_kdp_map->SetBinContent(uMt, dlm_kdp_map->GetBin(1, fDD), dlm_kdp_2D->GetBinContent(uMt, 0));
          }

          TString base_name = TString::Format("%s_%s_PR35.78_KR%.2f_DD%.2f", sSPAR.Data(), sTYPE.Data(), fKR, fDD);

          fSummary.cd();
          TH1F* hSource = new TH1F("hSource_"+base_name, "hSource_"+base_name, 4096, 0, 64);
          TGraph* gSource = new TGraph();
          gSource->SetName("gSource_"+base_name);
          double eval_at[2];
          eval_at[0] = 1500;
          eval_at[1] = kstar_cutoff*0.5;
          KdpPars current_kdp = dlm_kdp_2D->Eval(eval_at);
          //KdpPars current_kdp = dlm_kdp_2D->GetBinContent(2, 0);
          //dlm_kdp_2D->GetBinContent(1, 0).Print();
          //dlm_kdp_2D->GetBinContent(2, 0).Print();
          //dlm_kdp_2D->GetBinContent(3, 0).Print();
          //current_kdp.Print();
          //printf("%f\n", dlm_kdp_2D->GetBinCenter(0, 0));
          //printf("%f\n", dlm_kdp_2D->GetBinCenter(0, 1));
          //printf("%f\n", dlm_kdp_2D->GetBinCenter(0, 2));
          //printf("%f\n", dlm_kdp_2D->GetBinCenter(0, 3));
          //printf("%f\n", dlm_kdp_2D->GetBinCenter(0, 4));
          for(unsigned uRad=0; uRad<hSource->GetNbinsX(); uRad++){
            double rstar = hSource->GetBinCenter(uRad+1);
            double src = PoissonSum(rstar, current_kdp);
            hSource->SetBinContent(uRad+1, src);
            gSource->SetPoint(uRad, rstar, src);
          }
          DLM_Histo<float>* dlmSource = Convert_TH1F_DlmHisto(hSource);

          DLM_CommonAnaFunctions AnalysisObject;
          AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());
          double kMin = 0;
          double kMax = 400;
          unsigned kSteps = 40;

          //printf("hello\n");
          CATS Kd_DG_ER_CECA;
          Kd_DG_ER_CECA.SetMomBins(kSteps, kMin, kMax);
          AnalysisObject.SetUpCats_Kd(Kd_DG_ER_CECA,"DG_ER","",0,0);
          DLM_HistoSource dlmCecaSource(dlmSource);
          Kd_DG_ER_CECA.SetAnaSource(CatsSourceForwarder, &dlmCecaSource, 0);
          Kd_DG_ER_CECA.SetUseAnalyticSource(true);
          Kd_DG_ER_CECA.SetAutoNormSource(false);
          Kd_DG_ER_CECA.SetNormalizedSource(true);
          Kd_DG_ER_CECA.SetNotifications(CATS::nWarning);   
          Kd_DG_ER_CECA.KillTheCat();
          TGraph gCk_DG_ER;
          gCk_DG_ER.SetName("gCk_DG_ER_"+base_name);
          for(unsigned uKstar=0; uKstar<kSteps; uKstar++){
            double kstar = Kd_DG_ER_CECA.GetMomentum(uKstar);
            gCk_DG_ER.SetPoint(uKstar, kstar, Kd_DG_ER_CECA.GetCorrFun(uKstar));
          }


          CATS Kd_SW_ER_CECA;
          Kd_SW_ER_CECA.SetMomBins(kSteps, kMin, kMax);
          AnalysisObject.SetUpCats_Kd(Kd_SW_ER_CECA,"SW_ER","",0,0);
          //DLM_HistoSource dlmCecaSource(dlmSource);
          Kd_SW_ER_CECA.SetAnaSource(CatsSourceForwarder, &dlmCecaSource, 0);
          Kd_SW_ER_CECA.SetUseAnalyticSource(true);
          Kd_SW_ER_CECA.SetAutoNormSource(false);
          Kd_SW_ER_CECA.SetNormalizedSource(true);
          Kd_SW_ER_CECA.SetNotifications(CATS::nWarning);   
          Kd_SW_ER_CECA.KillTheCat();
          TGraph gCk_SW_ER;
          gCk_SW_ER.SetName("gCk_SW_ER_"+base_name);
          for(unsigned uKstar=0; uKstar<kSteps; uKstar++){
            double kstar = Kd_SW_ER_CECA.GetMomentum(uKstar);
            gCk_SW_ER.SetPoint(uKstar, kstar, Kd_SW_ER_CECA.GetCorrFun(uKstar));
          }
/*
          CATS Kd_DG_ER_CECA_kstar;
          Kd_DG_ER_CECA_kstar.SetMomBins(kSteps, kMin, kMax);
          AnalysisObject.SetUpCats_Kd(Kd_DG_ER_CECA_kstar,"DG_ER","",0,0);
          DLM_MtKstar_KdpSource source_kdp(dlm_kdp_2D_kstar);
          Kd_DG_ER_CECA_kstar.SetUseAnalyticSource(true);
          Kd_DG_ER_CECA_kstar.SetMomentumDependentSource(true);
          Kd_DG_ER_CECA_kstar.SetAnaSource(CatsSourceForwarder, &source_kdp, 2);
          Kd_DG_ER_CECA_kstar.SetAnaSource(0, 1500);
          Kd_DG_ER_CECA_kstar.SetAutoNormSource(false);
          Kd_DG_ER_CECA_kstar.SetNormalizedSource(true);
          Kd_DG_ER_CECA_kstar.SetNotifications(CATS::nWarning);   
          Kd_DG_ER_CECA_kstar.KillTheCat();
          TGraph gCk_DG_ER_kstar;
          gCk_DG_ER_kstar.SetName("gCk_DG_ER_kstar_"+base_name);
          for(unsigned uKstar=0; uKstar<kSteps; uKstar++){
            double kstar = Kd_DG_ER_CECA_kstar.GetMomentum(uKstar);
            gCk_DG_ER_kstar.SetPoint(uKstar, kstar, Kd_DG_ER_CECA_kstar.GetCorrFun(uKstar));
          }
*/
          CATS Kd_DG_FCA_CECA;
          Kd_DG_FCA_CECA.SetMomBins(kSteps, kMin, kMax);
          AnalysisObject.SetUpCats_Kd(Kd_DG_FCA_CECA,"DG_FCA","",0,0);
          //DLM_HistoSource dlmCecaSource(dlmSource);
          Kd_DG_FCA_CECA.SetAnaSource(CatsSourceForwarder, &dlmCecaSource, 0);
          Kd_DG_FCA_CECA.SetUseAnalyticSource(true);
          Kd_DG_FCA_CECA.SetAutoNormSource(false);
          Kd_DG_FCA_CECA.SetNormalizedSource(true);
          Kd_DG_FCA_CECA.SetNotifications(CATS::nWarning);   
          Kd_DG_FCA_CECA.KillTheCat();
          TGraph gCk_DG_FCA;
          gCk_DG_FCA.SetName("gCk_DG_FCA_"+base_name);
          for(unsigned uKstar=0; uKstar<kSteps; uKstar++){
            double kstar = Kd_DG_FCA_CECA.GetMomentum(uKstar);
            gCk_DG_FCA.SetPoint(uKstar, kstar, Kd_DG_FCA_CECA.GetCorrFun(uKstar));
          }

          //printf("hello\n");
          fSummary.cd();
          gSource->Write();
          hSource->Write();
          //gCk_DG_ER_kstar.Write();
          gCk_DG_ER.Write();
          gCk_SW_ER.Write();
          gCk_DG_FCA.Write();
          //printf("hello2\n");
          delete hSource;
          delete gSource;
          delete dlmSource;
          delete dlm_kdp_2D;
          //printf("hello3\n");
        }//fDD
        //printf("hello4\n");
        dlm_kdp_map->QuickWrite(TString::Format("%s/dKaon/Frascati/Sources/kdp_map_%s_%s_PR35.78_KR%.2f.dlm",
                  GetCernBoxDimi(),sSPAR.Data(), sTYPE.Data(), fKR), true);
        delete dlm_kdp_map;
        


      }
    }
  }


}

//it creates a TGraph, which though will have to be deleted later by you
//example:  TGraph* MyGraph = MyOwn_Kd_v2_InterpolateCk(...);
//          --- do something ---
//          delete MyGraph;
//InputFileName -> full path to the kdp_maps file that you use
//FSI_type -> DG_ER or SW_ER or DG_FCA
//CECA_type -> FIT or NLO19
//d_type -> Kd or KdReso
//KR_val -> 47.16 or 52.40  or 57.64
//binning -> the currectly available delay_step
TGraph* MyOwn_Kd_v2_InterpolateCk(TString InputFileName, TString FSI_type, TString CECA_type, TString d_type, double KR_val, double delay_val, double delay_step = 0.5){
  TFile InputFile(InputFileName, "read");

  TGraph* gCkFinal = NULL;

  //if we are almost demanding a specific existing step
  if( fabs((delay_val/delay_step) - TMath::Nint(delay_val/delay_step))<1e-4 ){
    InputFile.Get(TString::Format("gCk_"));

    return gCkFinal;
  }

}

void Check_Source_Means(){
  TFile fSummary(TString::Format("%s/dKaon/Frascati/Sources/kdp_maps_full_v2.root",GetCernBoxDimi()), "read");

  TH1F* hSource_A_up = (TH1F*)fSummary.Get("hSource_FIT_KdReso_PR35.78_KR47.16_DD1.00");
  TH1F* hSource_A_low = (TH1F*)fSummary.Get("hSource_FIT_KdReso_PR35.78_KR57.64_DD1.00");
  TH1F* hSource_B_up = (TH1F*)fSummary.Get("hSource_FIT_Kd_PR35.78_KR47.16_DD1.00");
  TH1F* hSource_B_low = (TH1F*)fSummary.Get("hSource_FIT_Kd_PR35.78_KR57.64_DD1.00"); 

  double mean_A_up = hSource_A_up->GetMean();


  TF1* fSource = new TF1("fSource",NormDoubleGaussSourceTF1,0,100,4);
  fSource->FixParameter(3,1);
  fSource->FixParameter(0,1.10+0.04);
  fSource->FixParameter(1,2.14+0.03);
  fSource->FixParameter(2,0.76);
  double prx_mean_up = fSource->Mean(0, 100);

  fSource->FixParameter(3,1);
  fSource->FixParameter(0,1.10-0.04);
  fSource->FixParameter(1,2.14-0.07);
  fSource->FixParameter(2,0.76);
  double prx_mean_low = fSource->Mean(0, 100);

  double prx_mean = 0.5*(prx_mean_up+prx_mean_low);
  double prx_err = 0.5*(prx_mean_up-prx_mean_low);

  double SceA_up = hSource_A_up->GetMean();
  double SceA_low = hSource_A_low->GetMean();
  double SceA_mean = 0.5*(SceA_up+SceA_low);
  double SceA_err = 0.5*(SceA_up-SceA_low);

  double SceB_up = hSource_B_up->GetMean();
  double SceB_low = hSource_B_low->GetMean();
  double SceB_mean = 0.5*(SceB_up+SceB_low);
  double SceB_err = 0.5*(SceB_up-SceB_low);

  printf("PRX: %.3f +/- %.3f\n",prx_mean,prx_err);
  printf("SCA: %.3f +/- %.3f\n",SceA_mean,SceA_err);
  printf("SCB: %.3f +/- %.3f\n",SceB_mean,SceB_err);


  //TH1F* hSource_DG_up;
  //TH1F* hSource_DG_low;
}

void MyOwn_Kd_v2_CreateCk(TString InputFile, TString OutputFile){

}

//creates a KDP source, where the dimensions are mt, kstar, tau(delay)
void SetUp_Kdp_Kd(){

  DLM_Histo<KdpPars>* dlm_kdp_3D = NULL;
  //cout << sizeof(KdpPars) << endl;
  
  double Chi2;
  TH1F* hSrc;

  int tau_min = 0;
  int tau_max = 0;
  int tau_absolute_max = 0;

  std::vector<TString> list_of_files;
  for(int iTau=tau_min; iTau<=tau_max; iTau++){
    //printf("iTau = %i\n",iTau);
    list_of_files.push_back(TString::Format("%s/dKaon/Source/Full/Eta0.8_Kd_ET1_PR1_DD%i.0_EF-1402.Ghetto_kstar_rstar_mT",GetCernBoxDimi(),iTau));
  }

  for(int iTau=tau_min; iTau<=tau_max; iTau++){
    DLM_Histo<float> dlmSrc;
    dlmSrc.QuickLoad(list_of_files.at(iTau));
    dlmSrc.ComputeError();

    //takes a 3D histo of Kstar Rstar Mt and returns a 2D kdp histo of Mt Kstar (in that order)
    printf("hi\n");
    DLM_Histo<KdpPars>* dlm_kdp_2D = Convert_3Dsource_Kdp(dlmSrc,true,3,96,3);
    //we create the 3D kdp histo, including tau, at the first iteration
    if(dlm_kdp_3D==NULL){
      printf("A\n");
      dlm_kdp_3D = new DLM_Histo<KdpPars>();
      dlm_kdp_3D->SetUp(3);
      printf("B\n");

      dlm_kdp_3D->SetUp(0, *dlm_kdp_2D, 0);
      printf("C\n");
      dlm_kdp_3D->SetUp(1, *dlm_kdp_2D, 1);
      printf("D\n");
      dlm_kdp_3D->SetUp(2, tau_absolute_max+1, -0.5, double(tau_absolute_max)+0.5);
      printf("E\n");
      dlm_kdp_3D->Initialize();
    }

    //we fill up the whole big histo
    for(unsigned uMt=0; uMt<dlm_kdp_2D->GetNbins(0); uMt++){
      for(unsigned uKstar=0; uKstar<dlm_kdp_2D->GetNbins(1); uKstar++){
        //printf("%u %u %i\n",uMt, uKstar, iTau);
        dlm_kdp_3D->SetBinContent(uMt, uKstar, iTau, dlm_kdp_2D->GetBinContent(uMt, uKstar));
      }
    }
  
    delete dlm_kdp_2D;
  }

  dlm_kdp_3D->QuickWrite(TString::Format("%s/dKaon/SetUp_Kdp_Kd/Full_v1/dlm_kdp_3D.dlm",GetCernBoxDimi()).Data());
  delete dlm_kdp_3D;
}


DLM_Ck* fitter_Kd_Ck=NULL;
double fitter_Kd(double* x, double* pars){
  return fitter_Kd_Ck?DLM_Baseline(x,pars)*(fitter_Kd_Ck->Eval(*x)*pars[5] + 1 - pars[5]):0;
}

//a fit with a kstar integrated CECA source
//we take stat + syst in quadr., and use a castrated pol3 for the BL.
//further, we make full bootstrap, including variations on:
// - fit range
// - lambda parameters (phys. hypothesis)
// - tau parameter (phys. hypothesis)
// - dlmck cutoff (phys. hyptothesis)
// - potential (phys. hyptothesis)
//for each tau (in steps of 0.5 fm) we plot the corresponding best chi2, and evaluate the allowed limits
void Fit_WithCECA_kstarInt(const unsigned NumBootstrap, const unsigned SEED){

  std::vector<TString> src_var = {"ceca_default"};
  std::vector<double> tau = {0.0, 0.5, 1.0};//
  std::vector<double> lam_gen = {0.901, 0.921, 0.941};//
  //std::vector<double> fit_max = {1200, 1500, 1800};
  std::vector<double> fit_max = {800, 1000, 1200, 1400, 1600, 1800, 2000};

  std::vector<double> femto_converge = {500, 700, 900};
  const double cats_min = 0;
  const double cats_max = 400;
  const unsigned cats_steps = 20;
  const double ck_min = 0;
  const double ck_max = 2000;
  const unsigned ck_steps = 100;

  std::vector<TString> pot_var = {"DG_ER", "DG_FCA"};

  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());

  //here we have [variation][tau]
  DLM_HistoSource*** dlmCecaSource = new DLM_HistoSource** [src_var.size()];
  for(unsigned uSrc=0; uSrc<src_var.size(); uSrc++){
    dlmCecaSource[uSrc] = new DLM_HistoSource* [tau.size()];
    for(unsigned uTau=0; uTau<tau.size(); uTau++){
      TFile fSrcInput(TString::Format("%s/FunWithCeca/Ceca_KdReso_FAST/Eta0.8_KdReso_ET1_PR1_DD%.1f_EF-1402.root",GetFemtoOutputFolder(), tau.at(uTau)),"read");
      TH1F* GhettoFemto_rstar  = (TH1F*)fSrcInput.Get("GhettoFemto_rstar");
      DLM_Histo<float>* dlmGhettoFemto_rstar = Convert_TH1F_DlmHisto(GhettoFemto_rstar);
      dlmCecaSource[uSrc][uTau] = new DLM_HistoSource(*dlmGhettoFemto_rstar);
      delete dlmGhettoFemto_rstar;
      delete GhettoFemto_rstar;
      fSrcInput.Close();
    }
  }

  TGraphErrors* gKd_stat;
  TGraphErrors* gKd_syst;
  TGraphErrors* gKd_total;

  TFile fData(TString::Format("%s/dKaon/UnfoldedDataAndSyst.root",GetCernBoxDimi()), "read");
  gKd_stat = (TGraphErrors*)fData.Get("gDATAunfolded");
  gKd_syst = (TGraphErrors*)fData.Get("gSYSTunfolded");
  gKd_total = new TGraphErrors();
  gKd_total->SetName("gKd_total");
  printf("read up stuff\n");
  for(unsigned uPts=0; uPts<gKd_stat->GetN(); uPts++){
    double kstar;
    double ck;
    double stat;
    double syst;
    gKd_stat->GetPoint(uPts,kstar,ck);
    stat = gKd_stat->GetErrorY(uPts);
    syst = gKd_syst->GetErrorY(uPts);
    gKd_total->SetPoint(uPts, kstar, ck);
    gKd_total->SetPointError(uPts, 0, sqrt(stat*stat+syst*syst));
  }
  printf("done\n");

  TRandom3 rangen(SEED);

  TGraph GDATA;
  TGraph GFIT;
  TGraph GBL;
  TGraph GFEMTO;

  Float_t FIT_MAX;
  Float_t Kd_LAM_GEN;
  Float_t CKCUTOFF;
  TString Kd_POT_TYPE;
  TString Kd_SRC_TYPE;
  Float_t CHI2;
  Float_t DELAY;

  TFile fOutput(TString::Format("%s/dKaon/Fit_v1/fit_result.root",GetCernBoxDimi()), "recreate");
  TTree* Kd_Tree = new TTree("Kd_Tree","Kd_Tree");
  printf("1\n");
  Kd_Tree->Branch("gData","TGraph",&GDATA,32000,0);//
  printf("2\n");
  Kd_Tree->Branch("gFit","TGraph",&GFIT,32000,0);//
  printf("3\n");
  Kd_Tree->Branch("gBl","TGraph",&GBL,32000,0);//
  Kd_Tree->Branch("gFemto","TGraph",&GFEMTO,32000,0);//
  Kd_Tree->Branch("FitMax",&FIT_MAX,"FitMax/F");//
  Kd_Tree->Branch("DELAY",&DELAY,"DELAY/F");//
  Kd_Tree->Branch("Kd_LAM_GEN",&Kd_LAM_GEN,"Kd_LAM_GEN/F");//
  Kd_Tree->Branch("CkCutOff",&CKCUTOFF,"CkCutOff/F");//
  Kd_Tree->Branch("Kd_pot_type","TString",&Kd_POT_TYPE,8000,0);//
  Kd_Tree->Branch("Kd_src_type","TString",&Kd_SRC_TYPE,8000,0);//
  Kd_Tree->Branch("Chi2",&CHI2,"Chi2/F");//

  printf("tree done\n");

  //not really needed, just for QA in the output
  TH2F* h2d_data = new TH2F("h2d_data", "h2d_data", ck_steps, ck_min, ck_max, 2048, 0, 2);
  TH2F* h2d_fit = new TH2F("h2d_fit", "h2d_fit", ck_steps, ck_min, ck_max, 2048, 0, 2);

  printf("about to iter\n");
  for(unsigned uSrc=0; uSrc<src_var.size(); uSrc++){
    Kd_SRC_TYPE = src_var.at(uSrc);
    for(unsigned uTau=0; uTau<tau.size(); uTau++){
      DELAY = tau.at(uTau);
      for(unsigned uPot=0; uPot<pot_var.size(); uPot++){
        Kd_POT_TYPE = pot_var.at(uPot);
        CATS Kd_Cat;
        Kd_Cat.SetMomBins(cats_steps, cats_min, cats_max);
        AnalysisObject.SetUpCats_Kd(Kd_Cat,pot_var.at(uPot),"",0,0);//NLO_Coupled_S
        Kd_Cat.SetUseAnalyticSource(true);
        Kd_Cat.SetAnaSource(CatsSourceForwarder, dlmCecaSource[uSrc][uTau], 0);
        Kd_Cat.SetNotifications(CATS::nWarning);
        Kd_Cat.KillTheCat();

        DLM_Ck Kd_Ck(1,9,Kd_Cat,ck_steps,ck_min,ck_max);
        fitter_Kd_Ck = &Kd_Ck;

        for(unsigned uLam=0; uLam<lam_gen.size(); uLam++){
          Kd_LAM_GEN = lam_gen.at(uLam);
          for(unsigned uBoot=0; uBoot<NumBootstrap; uBoot++){
            CKCUTOFF = femto_converge.at(rangen.Integer(femto_converge.size()));
            Kd_Ck.SetCutOff(cats_max, CKCUTOFF);
            Kd_Ck.Update();
            //printf("uSrc=%u; uTau=%u; uPot=%u; uLam=%u; uBoot=%u\n",uSrc,uTau,uPot,uLam,uBoot);
            //p0.p1.p2.p3.p4.lam
            FIT_MAX = fit_max.at(rangen.Integer(fit_max.size()));
            TF1* fit_funtion = new TF1("fit_funtion",fitter_Kd,ck_min,FIT_MAX,6);
            
            fit_funtion->SetParameter(0,1);
            fit_funtion->SetParLimits(0,0.9,1.1);

            fit_funtion->FixParameter(1,0);

            fit_funtion->SetParameter(2,50);
            fit_funtion->SetParLimits(2,1e-5,1e5);

            fit_funtion->SetParameter(3,0);
            fit_funtion->SetParLimits(3,-1e-7,1e-7);

            fit_funtion->FixParameter(4,-1e6);

            fit_funtion->FixParameter(5,Kd_LAM_GEN);

            TGraphErrors* gKd_to_fit = new TGraphErrors();
            gKd_to_fit->SetName("gKd_to_fit");

            for(unsigned uPts=0; uPts<gKd_stat->GetN(); uPts++){
              double kstar;
              double ck;
              double ck_boot = -1;
              double err;
              gKd_total->GetPoint(uPts,kstar,ck);
              if(kstar>FIT_MAX) break;
              err = gKd_stat->GetErrorY(uPts);
              while(ck_boot<=0){
                ck_boot = rangen.Gaus(ck, err);
              }
              gKd_to_fit->SetPoint(uPts, kstar, ck_boot);
              gKd_to_fit->SetPointError(uPts, 0, err);
            }
            gKd_to_fit->Fit(fit_funtion,"Q, S, N, R, M");


            //TGraph* gTree_Kd_data = new TGraph();
            //gTree_Kd_data->SetName("gTree_Kd_data");
            //TGraph* gTree_Kd_fit = new TGraph();
            //gTree_Kd_fit->SetName("gTree_Kd_fit");
            //TGraph* gTree_Kd_bl = new TGraph();
            //gTree_Kd_bl->SetName("gTree_Kd_bl");            
            //TGraph* gTree_Kd_femto = new TGraph();
            //gTree_Kd_femto->SetName("gTree_Kd_femto");

            TF1* bl_function = new TF1("bl_function",fitter_Kd,ck_min,FIT_MAX,6);
            TF1* femto_function = new TF1("femto_function",fitter_Kd,ck_min,FIT_MAX,6);
            for(unsigned uPar=0; uPar<6; uPar++){
              bl_function->FixParameter(uPar, fit_funtion->GetParameter(uPar));
              femto_function->FixParameter(uPar, fit_funtion->GetParameter(uPar));
            }
            bl_function->FixParameter(5, 0);

            femto_function->FixParameter(0, 1);
            femto_function->FixParameter(1, 0);
            femto_function->FixParameter(2, -1e6);
            femto_function->FixParameter(3, -1e6);
            femto_function->FixParameter(4, -1e6);

            for(unsigned uPts=0; uPts<gKd_stat->GetN(); uPts++){
              double kstar;
              double ck;  
              gKd_to_fit->GetPoint(uPts,kstar,ck);
              if(kstar>FIT_MAX) break;
              GDATA.SetPoint(uPts, kstar, ck);
              GFIT.SetPoint(uPts, kstar, fit_funtion->Eval(kstar));
              GBL.SetPoint(uPts, kstar, bl_function->Eval(kstar));
              GFEMTO.SetPoint(uPts, kstar, femto_function->Eval(kstar));

              h2d_data->Fill(kstar, ck);
              h2d_fit->Fill(kstar, fit_funtion->Eval(kstar));
            }


            //GDATA = gTree_Kd_data;
            //GFIT = gTree_Kd_fit;
            //GBL = gTree_Kd_bl;
            //GFEMTO = gTree_Kd_femto;

            CHI2 = fit_funtion->GetChisquare();

            Kd_Tree->Fill();

            delete fit_funtion;
            delete gKd_to_fit;
            //delete gTree_Kd_data;
            //delete gTree_Kd_fit;
            //delete gTree_Kd_bl;
            //delete gTree_Kd_femto;
            delete bl_function;
            delete femto_function;
          }
        }
      }
    }
  }

  fOutput.cd();
  Kd_Tree->Write();
  h2d_data->Write();
  h2d_fit->Write();

  //printf("we are out\n");
  for(unsigned uSrc=0; uSrc<src_var.size(); uSrc++){
    for(unsigned uTau=0; uTau<tau.size(); uTau++){
      delete dlmCecaSource[uSrc][uTau];
    }
    delete [] dlmCecaSource[uSrc];
  }
  delete [] dlmCecaSource;
  delete gKd_total;
  delete Kd_Tree;
  delete h2d_data;
  delete h2d_fit;
}


void pi_d_source(){


  double reff_min=1000; double reff_max=0; double reff_def;


  unsigned NumRadBins = 1024;
  double rMin = 0;
  double rMax = 12;
  TH1F* hSourceUpper = new TH1F("hSourceUpper","hSourceUpper",NumRadBins,rMin,rMax);
  TH1F* hSourceLower = new TH1F("hSourceLower","hSourceLower",NumRadBins,rMin,rMax);
  TH1F* hSource = new TH1F("hSource","hSource",NumRadBins,rMin,rMax);

  for(unsigned uRad=0; uRad<NumRadBins; uRad++){
    hSourceUpper->SetBinContent(uRad+1,0);
    hSourceLower->SetBinContent(uRad+1,1e6);
  }

  std::vector<float> core_radii_mTint = {1.036,1.075,1.112};
  std::vector<float> core_radii_mT0 = {1.163,1.216,1.271};
  std::vector<float> core_radii_mT1 = {1.068,1.109,1.150};
  std::vector<float> core_radii_mT2 = {0.996,1.033,1.068};
  std::vector<float> core_radii_mT3 = {0.932,0.968,1.003};
  std::vector<float> core_radii_mT4 = {0.782,0.831,0.880};
  std::vector<float> vectors[] = {core_radii_mTint, core_radii_mT0, core_radii_mT1, core_radii_mT2, core_radii_mT3, core_radii_mT4};

  for(int uA=0; uA<1; uA++){for(int uB=0; uB<1; uB++){for(int uE=0; uE<3; uE++){
        //#for the sake of simplicity, make it 1.074 +/- 0.39 and min/avg/max = 1.035/1.074/1.113
        //std::vector<float> core_radii = {1.035,1.074,1.113};//mt integrated
        std::vector<float> core_radii = {0.782,0.831,0.880};
        //for (const auto& vec : vectors) {
        //    core_radii.insert(core_radii.end(), vec.begin(), vec.end());
        //}

        DLM_CommonAnaFunctions AnalysisObject;
        AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());

        unsigned NumMomBins = 60;
        double kMin = 0;
        double kMax = 300;
        CATS Kitty;
        Kitty.SetMomBins(NumMomBins, kMin, kMax);
        AnalysisObject.SetUpCats_pi_d(Kitty, "", "McGauss_ResoTM", 0, uA*10000+uB*1000+200+uE);
        DLM_CleverMcLevyResoTM* MagicSource = AnalysisObject.GetCleverMcLevyResoTM_pi_d();
        printf("%i %i %i\n",uA,uB,uE);
        
        for(unsigned uSrc=0; uSrc<core_radii.size(); uSrc++){
          if(uSrc%3==0){
            printf("mT%i-----------------\n",int(uSrc)/3-1);
          }
          double r_eff = GetReff(*MagicSource, core_radii.at(uSrc));
          if(uA==0&&uB==0&&uE==0&&uSrc==1) reff_def = r_eff;
          printf("rcore = %.3f; reff = %.3f\n", core_radii.at(uSrc), r_eff);
          if(reff_min>r_eff) reff_min = r_eff;
          if(reff_max<r_eff) reff_max = r_eff;
        }

        for(unsigned uSrc=0; uSrc<core_radii.size(); uSrc++){
          Kitty.SetAnaSource(0, core_radii.at(uSrc));
          for(unsigned uRad=0; uRad<NumRadBins; uRad++){
            double rstar = hSource->GetBinCenter(uRad+1);
            double sr = Kitty.EvaluateTheSource(10, rstar, 0);
            double sr_lower = hSourceLower->GetBinContent(uRad+1);
            double sr_upper = hSourceLower->GetBinContent(uRad+1);
            if(uA==0&&uB==0&&uE==0&&uSrc==1){
              hSource->SetBinContent(uRad+1, sr);
            }
            if(sr>sr_upper){
              hSourceUpper->SetBinContent(uRad+1,sr);
            }
            if(sr<sr_lower){
              hSourceLower->SetBinContent(uRad+1,sr);
            }
          }
        }

        //Kitty.KillTheCat();
      }
    }
  }

  TF1* fEffSrc = new TF1("fEffSrc",GaussSourceTF1,rMin,rMax,1);
  fEffSrc->SetParameter(0,reff_def);

  TFile fOutput(TString::Format("%s/pi_d/SourceEval/pi_d_source.root",GetCernBoxDimi()),"recreate");
  hSource->Write();
  hSourceLower->Write();
  hSourceUpper->Write();
  fEffSrc->Write();



  //EvaluateTheSource

  printf("reff_min = %.3f; reff_max = %.3f\n", reff_min, reff_max);
  printf("=> %.3f +/- %.3f\n", (reff_max+reff_min)*0.5, (reff_max-reff_min)*0.5);

  delete hSource;
  delete hSourceLower;
  delete hSourceUpper;

}

//a fitter inspired by marcel
//we create the shape of the delta, and can fit it
//in addition, we have a dimi-style baseline
//! CHECK ABOUT THE UNFOLDING !!
void pi_d_DLM_fitter1(){

}

void test_sill_ps(){
  TF1* fSill = new TF1("fSill",SillPhaseSpaceKstar,0,1000,7);
//par[0/1] = masses of the daughters
//par[2] = mass (mother)
//par[3] = width
//par[4] = Resonance pT, if zero no PS
//par[5] = Temperature, if zero no PS
  fSill->SetParameter(0, 1);
  fSill->SetParameter(1, Mass_pic);
  fSill->SetParameter(2, Mass_p);
  fSill->SetParameter(3, 1215);
  fSill->SetParameter(4, 100);
  fSill->SetParameter(5, 1000);
  fSill->SetParameter(6, 30);

//    TFile* fOutput = new TFile(TString::Format("%s/OtherTasks/Georgios_LXi_ResoTest/fOutput_%i_%.2f.root",GetFemtoOutputFolder(),int(SmoothSampling),CoreSize),"recreate");

  double Integral = fSill->Integral(0,1000);
  printf("Integral = %.2e\n",Integral);
  fSill->SetParameter(0,1./Integral);
  Integral = fSill->Integral(0,1000);
  printf("Integral = %.2e\n",Integral);

  TFile fOutput(TString::Format("%s/Deuteron/test_sill_ps.root",GetFemtoOutputFolder()), "recreate");
  fSill->Write();
}


void test_sill_invM_to_kstar(){
  double MASS = 1232;
  double WIDTH = 120;
  double AVG_PT = 1000;
  double EFF_TEMP = 30;
  double MASS_D1 = Mass_pic;
  double MASS_D2 = Mass_p;
  TH1F* hInvMass = new TH1F("hInvMass","hInvMass",1024,1000,3000);
  TH1F* hKstar = new TH1F("hKstar","hKstar",1024,0,2000);

  TH1F* hBoltzmannInvMass = new TH1F("hBoltzmannInvMass","hBoltzmannInvMass",1024,1000,3000);
  TH1F* hBoltzmannKstar = new TH1F("hBoltzmannKstar","hBoltzmannKstar",1024,0,2000);

  TH1F* hSillBoltzmann_IM = new TH1F("hSillBoltzmann_IM","hSillBoltzmann_IM",1024,1000,3000);
  TH1F* hSillBoltzmann_kstar = new TH1F("hSillBoltzmann_kstar","hSillBoltzmann_kstar",1024,0,2000);

  TF1* fSill_IM = new TF1("fSill_IM",Sill_IM,1000,3000,5);
  fSill_IM->SetParameter(0, -1e6);
  fSill_IM->SetParameter(1, MASS_D1);
  fSill_IM->SetParameter(2, MASS_D2);
  fSill_IM->SetParameter(3, MASS);
  fSill_IM->SetParameter(4, WIDTH);  
  fSill_IM->SetNpx(1024);

  TF1* fSill_kstar = new TF1("fSill_kstar",Sill_kstar,0,2000,5);
  fSill_kstar->SetParameter(0, -1e6);
  fSill_kstar->SetParameter(1, MASS_D1);
  fSill_kstar->SetParameter(2, MASS_D2);
  fSill_kstar->SetParameter(3, MASS);
  fSill_kstar->SetParameter(4, WIDTH);  
  fSill_kstar->SetNpx(1024);


  TF1* fBoltzmann_IM = new TF1("fBoltzmann_IM",Boltzmann_IM,1000,3000,3);
  fBoltzmann_IM->SetParameter(0, -1e6);
  fBoltzmann_IM->SetParameter(1, AVG_PT);
  fBoltzmann_IM->SetParameter(2, EFF_TEMP);
  fBoltzmann_IM->SetNpx(1024);

//par[0] = NORM
//par[1/2] = masses of the daughters
//par[3] = mass (mother)
//par[4] = width
//par[5] = avg pt, if zero no PS
//par[6] = Temperature, if zero no PS

  TF1* fSillBoltzmann_IM = new TF1("fSillBoltzmann_IM",SillBoltzmann_IM,1000,3000,7);
  fSillBoltzmann_IM->SetParameter(0, -1e6);
  fSillBoltzmann_IM->SetParameter(1, MASS_D1);
  fSillBoltzmann_IM->SetParameter(2, MASS_D2);
  fSillBoltzmann_IM->SetParameter(3, MASS);
  fSillBoltzmann_IM->SetParameter(4, WIDTH);
  fSillBoltzmann_IM->SetParameter(5, AVG_PT);
  fSillBoltzmann_IM->SetParameter(6, EFF_TEMP);  
  fSillBoltzmann_IM->SetNpx(1024);

  TF1* fSillBoltzmann_kstar = new TF1("fSillBoltzmann_kstar",SillBoltzmann_kstar,0,2000,7);
  fSillBoltzmann_kstar->SetParameter(0, -1e6);
  fSillBoltzmann_kstar->SetParameter(1, MASS_D1);
  fSillBoltzmann_kstar->SetParameter(2, MASS_D2);
  fSillBoltzmann_kstar->SetParameter(3, MASS);
  fSillBoltzmann_kstar->SetParameter(4, WIDTH);
  fSillBoltzmann_kstar->SetParameter(5, AVG_PT);
  fSillBoltzmann_kstar->SetParameter(6, EFF_TEMP);  
  fSillBoltzmann_kstar->SetNpx(1024);

  TRandom3 rangen(23);

  double my_convert = sqrt(pow(0.5*MASS - 0.5*(MASS_D1*MASS_D1+MASS_D2*MASS_D2)/MASS,2.) - pow(MASS_D1*MASS_D2/MASS,2.));
  cout << IM_to_kstar(MASS, MASS_D1, MASS_D2) << " vs " << my_convert << endl;
//return ;
  int NumIter = 100*1000;
  for(int iIter=0; iIter<NumIter; iIter++){
    double ran_im = fSill_IM->GetRandom(1000,3000,&rangen);
    double ran_kstar = IM_to_kstar(ran_im, MASS_D1, MASS_D2);
    hInvMass->Fill(ran_im);
    hKstar->Fill(ran_kstar);


    ran_im = fSillBoltzmann_IM->GetRandom(1000,3000,&rangen);
    ran_kstar = IM_to_kstar(ran_im, MASS_D1, MASS_D2);
    hSillBoltzmann_IM->Fill(ran_im);
    hSillBoltzmann_kstar->Fill(ran_kstar);


    ran_im = fBoltzmann_IM->GetRandom(1000,3000,&rangen);
    ran_kstar = IM_to_kstar(ran_im, MASS_D1, MASS_D2);
    hBoltzmannInvMass->Fill(ran_im);
    hBoltzmannKstar->Fill(ran_kstar);

    
    

  }

  hInvMass->Scale(fSill_IM->Integral(1000,3000)/hInvMass->Integral(),"width");
  hKstar->Scale(fSill_kstar->Integral(0,2000)/hKstar->Integral(),"width");

  hBoltzmannInvMass->Scale(fBoltzmann_IM->Integral(1000,3000)/hBoltzmannInvMass->Integral(),"width");
  hBoltzmannKstar->Scale(1./hBoltzmannKstar->Integral(),"width");

  hSillBoltzmann_IM->Scale(fSillBoltzmann_IM->Integral(1000,3000)/hSillBoltzmann_IM->Integral(),"width");
  hSillBoltzmann_kstar->Scale(fSillBoltzmann_kstar->Integral(0,2000)/hSillBoltzmann_kstar->Integral(),"width");

  TFile fOutput(TString::Format("%s/Deuteron/test_sill_invM_to_kstar.root",GetFemtoOutputFolder()), "recreate");
  hInvMass->Write();
  fSill_IM->Write();
  hKstar->Write();
  fSill_kstar->Write();

  hSillBoltzmann_IM->Write();
  fSillBoltzmann_IM->Write();
  hSillBoltzmann_kstar->Write();
  fSillBoltzmann_kstar->Write();

  hBoltzmannInvMass->Write();
  fBoltzmann_IM->Write();
  hBoltzmannKstar->Write();

}

//par[0] = NORM
//par[1/2] = masses of the daughters
//par[3] = mass (mother)
//par[4] = width
//par[5] = avg pt, if zero no PS
//par[6] = Temperature, if zero no PS


void clean_graph(TGraph* gInput){
  int num_pts = gInput->GetN();
  for(int iPts=num_pts-1; iPts>=0; iPts--){
    gInput->RemovePoint(iPts);
  }
}

//par[0] = source size
//par[1] = alpha par
//par[2] = lambda_genuine
//par[3] = fraction_d_from_Delta = frac_D
//par[4] = delta_amplitude
//par[5] = CkCutOff
//par[6] = CkConv
//par[7] = treat it as a flag for the fit. 0 = my style where we do NOT multiply C_Delta with C_femto, and 1 where we do so
//par[8] = Delta normalization, should be a dummy fixed to -1e6, giving instructions to the code to renorm
//par[9] = mass of daughter 1 (say the pion)
//par[10] = mass of daughter 2 (say the proton)
//par[11] = mass of the delta
//par[12] = width of the delta
//par[13] = avg pT of the daughters
//par[14] = effective temperature
//par[15] = just in case
//par[16] = norm
//par[17] = 0 to castrate the pol3
//par[18] = position of the max of the pol3
//par[19] = p3 parameter of the pol3
//par[20] = -1e6 to switch off the pol4
//NEW in JULY 2025:
//par[21] = amplitude of extra sill reso
//par[22] = mass of daughter 1 (d)
//par[23] = mass of daughter 1 (pi)
//par[24] = mass of an extra sill reso
//par[25] = width of an extra sill reso
//the total fit function is:
//Baseline*(lambda_femto*(1-frac_d)*Ck_femto + norm_delta*Ck_delta*(either 1 or Ck_femto) + lamda_flat)
DLM_CkDecomposition* Bulgaristan_CkDec=NULL;
TGraph* Kstar_Modifier = NULL;
double Bulgaristan_fit(double* kstar, double* par){
  if(!Bulgaristan_CkDec) return 0;
  double& MOM = kstar[0];
  //double lambda_flat = par[3];
  //double lambda_fd = 1.-lambda_flat;
  //double lambda_femto = par[2]*lambda_fd/(1.+par[2]);
  //double lambda_delta = lambda_fd - lambda_femto;
  

  Bulgaristan_CkDec->GetCk()->SetSourcePar(0, par[0]);
  Bulgaristan_CkDec->GetCk()->SetCutOff(par[5],par[6]);
  Bulgaristan_CkDec->Update();

  double Baseline = DLM_Baseline(kstar, &par[16]);
  //pid / pip kstar is around 1.1 => as kstar[0] is the pid, modified is pip, we need to divide here
  double modified_kstar = kstar[0]/1.1;//roughly what happens after the boost @ 120 MeV, see the original calc.
  //if(kstar[0]<=50) modified_kstar = kstar[0]/1.065;
  //else if(kstar[0]>=100) modified_kstar = kstar[0]/1.075;
  //else if(kstar[0]>=200) modified_kstar = kstar[0]/1.11;
  //else if(kstar[0]>=300) modified_kstar = kstar[0]/1.13;
  //else 
  modified_kstar = kstar[0]/Kstar_Modifier->Eval(kstar[0]);
  
  double Delta;
  //kaon shit
  if(par[11]<500){
    Delta = TMath::BreitWigner(MOM, par[11], par[12]);
  }
  else{
    Delta = SillBoltzmann_kstar(&modified_kstar, &par[8]);
  }
  double Femto = Bulgaristan_CkDec->EvalCk(MOM);
  double MadeUpResonance = 0;
  if(TMath::Nint(par[7])>=10){
    //par[0] = NORM
    //par[1/2] = masses of the daughters
    //par[3] = mass (mother)
    //par[4] = width
    MadeUpResonance = Sill_kstar(kstar, &par[21]);
  }
  //printf("MadeUpResonance = %f\n",MadeUpResonance);
  //printf("%f %f\n",MOM, Femto);
  //if(MOM==10){
  //  printf("%.3e * (%.3e*%.3e + %.3e*%.3e + %.3e)\n",Baseline, par[2], Femto, par[3], Delta, 1.-par[2]);
  //}
  //N.B. the Delta goes to zero, hence does not count as a proper lambda par!
  if(TMath::Nint(par[7])==0){
    return Baseline*(par[2]*(1.-par[3])*Femto + par[4]*Delta + (1.-par[2]) + MadeUpResonance);
  }
  else{
    return Baseline*(par[2]*(1.-par[3])*Femto + par[4]*Delta*Femto + (1.-par[2]) + MadeUpResonance);
  }
  
}

//if the |mT_bin|>100 => we deal with pi^- d
//if its >0, we deal with pi^+ d
//-1 = pip mt int
//-101 = pim mt int
void BulgarianIndianGhetto(TString Description, std::vector<int> fit_types, int mt_bin, int NumIter, bool Bootstrap=true, bool DataVar=true, bool FitVar=true, int SEED=0){
printf("WHEN YOU FIT THE NEW DATA USE THE MT RANGES IN SUBFOLDER *2024-10-22\n AND UPDATE THE WIDTH AND TEMP BASED ON THE NEW RANGES!!! USE PythonPlot_Pi_d.py to get them out\n");
  if(mt_bin%100<-1 || mt_bin%100>5){
    printf("Bad mT bin\n");
    abort();
  }
  //these are the first file you got from Bhawani. On 29th Oct 2024 this was updated and now completely
  //converd by the if statements below, these paths here are for bookkeeping only and have no affect right now!
  TString InputDataFileName = 
    TString::Format("%s/pi_d/FitOct2024_Files/Systematics_pp.root",GetCernBoxDimi());
  TString DefaultHistoName = "histDefault";
  TString VarDirName = "Raw";
  TString VarHistoBaseName = "histVar_";
  TString MeListName1 = "PairDist";
  TString MeListName2 = "PairReweighted";
  TString MeHistName = "hTotalMEpairs _Rebinned_5_Reweighted";
  TString SeHistName = "";
  //pip d mt diff, as used for Approval at CF
  /*
  if(mt_bin>=0 && mt_bin<100 && false){
    InputDataFileName = 
    TString::Format("%s/pi_d/FitOct2024_Files/mT_2024-10-11/mT06102024Signal_sum_Mt_%i.root",GetCernBoxDimi(),mt_bin+1);
    DefaultHistoName = "hCk_ReweightedMeV_2";
    if(DataVar==true){
      printf("ERROR: At the moment we do not have mt variations\n");
      abort();
    }
  }
  //pim d mt diff, as used for Approval at CF
  else if(mt_bin>=100 && mt_bin<200 && false){
    InputDataFileName = 
    TString::Format("%s/pi_d/FitOct2024_Files/AntiPionDeuteron/22102024FilesSignal_antisum_Mt_%i.root",GetCernBoxDimi(),mt_bin%100+1);
    DefaultHistoName = "hCk_ReweightedMeV_2";
    MeListName1 = "PairDist";
    MeListName2 = "PairReweighted";
    //hTotalMEapairs _Rebinned_5_Reweighted
    //MeHistName = TString::Format("MEmTMult_%i_Particle1_Particle20 _Rebinned_5_Reweighted", mt_bin%100+1);
    MeHistName = TString::Format("hTotalMEapairs _Rebinned_5_Reweighted");
    if(DataVar==true){
      printf("ERROR: At the moment we do not have mt variations\n");
      abort();
    }
  }
  //pip d mt diff, as used for Approval at PF, includes systematics
  else if(mt_bin>=0 && mt_bin<100 && true){
    InputDataFileName = 
    TString::Format("%s/pi_d/FitNov2024_Files/Files/PionDeuteron/mT_differential/mTBin%i/Systematics_pp.root",GetCernBoxDimi(),abs(mt_bin)%100+1);
  }
  //pim d mt diff, as used for Approval at CF, includes systematics
  else if(mt_bin>=100 && mt_bin<200 && true){
    InputDataFileName = 
    TString::Format("%s/pi_d/FitNov2024_Files/Files/AntiPionDeuteron/mT_differential/mTBin%i/Systematics_pp.root",GetCernBoxDimi(),abs(mt_bin)%100+1);
  }
  //mt integrated pip - d, 2024-10-29
  else if(mt_bin==-1){
  //  TString::Format("%s/pi_d/FitOct2024_Files/2024-10-29/20MeVBins_pip_d/CF_Piondeuteron_Var0.root",GetCernBoxDimi());
  //  DefaultHistoName = "hCk_ReweightedPidVar0MeV_0";
  //  VarDirName = "";
  //  VarHistoBaseName = "histVar_";

  }
  //mt integrated pim - d, 2024-10-29
  else if(mt_bin==-101){

  }
  else{
    printf("Bad mT bin %i\n", mt_bin);
    abort();    
  }

  //TString InputConversionFileName = 
  //  TString::Format("%s/pi_d/Bhawani/ForDimi/ForDimi/Output/ResonanceOutput_2024-10-10.root",GetCernBoxDimi());

  //N.B. PURE vars, the default does not count here
  double NumDataVars = 26;
  if(mt_bin==-1){
    //NumDataVars = 27;
  }
  else if(mt_bin==-101){
    //NumDataVars = 35;
  }
  */
  //N.B. Sometimes variations are missing if they did not pass some criteria
  //what this is is the MAX number of data variations, if we try to open non-existing
  //variation we will resample
  int NumDataVars = 44;


  const double max_kstar = 700;
  const unsigned max_mom_bins = 35;
  //std::vector<float> fit_range = {520, 440, 600};
  std::vector<float> fit_range = {600, 500};
  std::vector<float> eff_source_radii;// = {1.51, 1.39, 1.63};
  std::vector<float> lambda_gen_relvar = {1.0, 0.95, 1.05};
  //the ratio of the lambda pars (yields) of pairs counting as primary
  //or decay products of the delta. If negative, the value is fitted, with 
  //starting value of fabs(1st elements) and limits fabs(2nd and 3rd element)
  std::vector<float> amplitude_delta = {-1e8, -1, -1e14};
  //these are the scenario B and A (B = 2/3 became my default it seems) for how the pT of the p-pi relates to d-pi
  std::vector<float> pT_scale = {2./3., 3./4.};
  std::vector<TString> InputConversionFileName = 
    { TString::Format("%s/pi_d/Bhawani/ForDimi/ForDimi/Output/ResonanceOutput_B.root",GetCernBoxDimi()), 
      TString::Format("%s/pi_d/Bhawani/ForDimi/ForDimi/Output/ResonanceOutput_A.root",GetCernBoxDimi())};
  float CKCUTOFF = 380;
  const double CkConv = 700;

  //0 - do NOT multiply C_Delta
  //1 - do multiply C_Delta by C_genuine
  //std::vector<int> fit_types = {1,0};

  //based on ME of a local test, for pairs below k* 300 MeV
  std::vector<float> piN_pT = {560, 800, 1000, 1200, 1500};
  float piN_pT_avg = 900;
  //these are all per mt bin
  std::vector<std::vector<float>> mT_low;
  std::vector<std::vector<float>> mT_up;
  std::vector<float> mT_low_pip_v1 = {1030,1160,1280,1400,1550};
  std::vector<float> mT_up_pip_v1 = {1160,1280,1400,1550,2400};
  std::vector<float> mT_low_pim_v1 = {1030,1160,1240,1380,1520};
  std::vector<float> mT_up_pim_v1 = {1160,1240,1380,1520,2400};

  mT_low.push_back(mT_low_pip_v1);
  mT_low.push_back(mT_low_pim_v1);

  mT_up.push_back(mT_up_pip_v1);
  mT_up.push_back(mT_up_pim_v1);

  double avg_mt = 1270;
  if(mt_bin>=0 && mt_bin<1000){
    avg_mt = 0.5*(mT_up.at(mt_bin/100).at(mt_bin%100)+mT_low.at(mt_bin/100).at(mt_bin%100));
  }
  double avg_mass_pid = 0.5*(Mass_d+Mass_pic);
  double kay_tee = sqrt(avg_mt*avg_mt-avg_mass_pid*avg_mass_pid);
  
  std::vector<float> lambda_gen_avg = {0.850*0.88, 0.874*0.88, 0.888*0.88, 0.897*0.88, 0.918*0.88};
  double lambda_gen_mTint = 0.88*0.88;//avg_mt of 1.27

  //take 0.75 or 0 weight for the Delta++, based on the system (pip or pim)
  //N.B. for pim we should actually have 1/2 D0 and 1/2 D-, but we dont have access to it
  double Dpp_weight = 0.75;
  if(fabs(mt_bin)>=100 && mt_bin<1000) Dpp_weight = 0;
  std::vector<float> marcel_gamma_D0 = {74.11, 74.06, 74.02, 70.19, 90.29};
  std::vector<float> marcel_temp_D0 = {18.51, 18.27, 18.04, 16.49, 17.41};
  std::vector<float> marcel_gamma_Dpp = {109.63, 106.28, 103.21, 99.56, 86.01};
  std::vector<float> marcel_temp_Dpp = {26.49, 26.71, 26.92, 26.30, 23.25};
  //these values are avaraged over Delta++ and 0, with weight of 0.75 for the ++
  std::vector<float> marcel_gamma = {100.75, 98.23, 95.92, 92.22, 87.53};
  std::vector<float> marcel_temp = {24.49, 24.60, 24.70, 23.84, 21.79};
  float marcel_avg_gamma = 95;
  float marcel_avg_temp = 24.5;

  std::vector<float> reff_avg_val = {1.673, 1.552, 1.464, 1.390, 1.226};
  std::vector<float> reff_avg_err = {0.132, 0.119, 0.115, 0.118, 0.141};
  double reff_mTint_val = 1.513;//avg_mt of 1.27
  double reff_mTint_err = 0.117;
  if(mt_bin>=0 && mt_bin<1000){
    eff_source_radii.push_back(reff_avg_val.at(mt_bin%100));
    eff_source_radii.push_back(reff_avg_val.at(mt_bin%100)-reff_avg_err.at(mt_bin%100));
    eff_source_radii.push_back(reff_avg_val.at(mt_bin%100)+reff_avg_err.at(mt_bin%100));
  }
  else if(mt_bin==1000){
    eff_source_radii.push_back(1.35);
  }
  else{
    eff_source_radii.push_back(reff_mTint_val);
    eff_source_radii.push_back(reff_mTint_val-reff_mTint_err);
    eff_source_radii.push_back(reff_mTint_val+reff_mTint_err);
  }


  TH2F* hResolution_pd = NULL;
  TH2F* hTemp = NULL;
  //printf("...\n");

  

  //PF
/*
  TFile* fReso = new TFile(TString::Format("%s/pi_d/Bhawani/PdMomResoMEinjected.root",GetCernBoxDimi()), "read");
  //printf("%p\n",fReso);
  hTemp = (TH2F*)fReso->Get("MomentumResolutionME_Particle0_Particle2");
  //printf("%p\n",hTemp);
  gROOT->cd();
  hResolution_pd = (TH2F*)hTemp->Clone("hResolution_pd");
  hResolution_pd->GetXaxis()->SetLimits(hTemp->GetXaxis()->GetXmin()*1000.,hTemp->GetXaxis()->GetXmax()*1000.);
  hResolution_pd->GetYaxis()->SetLimits(hTemp->GetYaxis()->GetXmin()*1000.,hTemp->GetYaxis()->GetXmax()*1000.);
*/

   
  //December 2024
  TFile* fReso = new TFile(TString::Format("%s/pi_d/Bhawani/AnalysisResults3996.root",GetCernBoxDimi()), "read");
//fReso->ls();
  TDirectoryFile* fDir_reso = NULL;
  TList* fList1_reso = NULL;
  TList* fList2_reso = NULL;
  TList* fList3_reso = NULL;

  fDir_reso = (TDirectoryFile*)(fReso->FindObjectAny("MBResultsQA0"));
//fDir_reso->ls();
  fDir_reso->GetObject("MBResultsQA0", fList1_reso);
//fList1_reso->Print();
  fList2_reso = (TList*)fList1_reso->FindObject("PairQA");
//fList2_reso->Print();
  fList3_reso = (TList*)fList2_reso->FindObject("QA_Particle0_Particle2");
//fList3_reso->Print();
//printf("-----\n");
  hTemp = (TH2F*)fList3_reso->FindObject("MomentumResolutionSE_Particle0_Particle2");
  //TFile fTMP("wtf.root","recreate");
  gROOT->cd();
  hResolution_pd = (TH2F*)hTemp->Clone("hResolution_pd");
//printf("1-----\n");
  fList3_reso = (TList*)fList2_reso->FindObject("QA_Particle1_Particle3");
  hTemp = (TH2F*)fList3_reso->FindObject("MomentumResolutionSE_Particle1_Particle3");
  gROOT->cd();
  hResolution_pd->Add(hTemp);
//printf("2-----\n");
  fList3_reso = (TList*)fList2_reso->FindObject("QA_Particle0_Particle3");
  hTemp = (TH2F*)fList3_reso->FindObject("MomentumResolutionSE_Particle0_Particle3");
  gROOT->cd();
  hResolution_pd->Add(hTemp);
//printf("3-----\n");
  fList3_reso = (TList*)fList2_reso->FindObject("QA_Particle1_Particle2");
  hTemp = (TH2F*)fList3_reso->FindObject("MomentumResolutionSE_Particle1_Particle2");
  gROOT->cd();
  hResolution_pd->Add(hTemp);
//printf("4-----\n");
  hResolution_pd->GetXaxis()->SetLimits(hTemp->GetXaxis()->GetXmin()*1000.,hTemp->GetXaxis()->GetXmax()*1000.);
  hResolution_pd->GetYaxis()->SetLimits(hTemp->GetYaxis()->GetXmin()*1000.,hTemp->GetYaxis()->GetXmax()*1000.);
//printf("5-----\n");
  //fTMP.cd();
  //hResolution_pd->Write();
  gROOT->cd();
  
  delete fReso;

  TRandom3 rangen(SEED);

/*
  //printf("file %s\n",InputDataFileName.Data());
  TFile fInputData(InputDataFileName, "read");
  TDirectoryFile* fDir = NULL;
  TH1F* hData = NULL;
  TList* fList1 = NULL;
  TList* fList2 = NULL;
  TH1D* hME = NULL;
  fList1 = (TList*)(fInputData.FindObjectAny(MeListName1));
  if(fList1) fList2 = (TList*)fList1->FindObject(MeListName2);
  //fList2->ls();
  if(fList2) hME = (TH1D*)fList2->FindObject(MeHistName);
  if(hME){
    hME->GetXaxis()->SetLimits(hME->GetXaxis()->GetXmin()*1000.,hME->GetXaxis()->GetXmax()*1000.);
  }
*/

  //GET THE PTs
  TFile fInput_pT(TString::Format("%s/pi_d/FitOct2024_Files/AnalysisResults_for_pT.root",GetCernBoxDimi()), "read");
  //printf("hi\n");
  TDirectory* fDir_temp = NULL;
  fDir_temp = (TDirectoryFile*)(fInput_pT.FindObjectAny("HMResults0"));
  //printf("%p\n",fDir_temp);
  TList* fList_temp1;
  //fDir_temp->ls();
  fDir_temp->GetObject("HMResults0", fList_temp1);
  //printf("%p\n",fList_temp1);

  //TList* fList_temp2 = (TList*)fList_temp1->FindObject("HMResults0");
  //    printf("%p\n",fList_temp2);

  TList* fList_P02 = (TList*)fList_temp1->FindObject("Particle0_Particle2");
  TList* fList_P13 = (TList*)fList_temp1->FindObject("Particle1_Particle3");
  TList* fList_P03 = (TList*)fList_temp1->FindObject("Particle0_Particle3");
  TList* fList_P12 = (TList*)fList_temp1->FindObject("Particle1_Particle2");
  int which_mt_bin = abs(mt_bin%100) + 1;
  //if(abs(mt_bin)>=100){which_mt_bin -= 100;}
  if(mt_bin<0) {which_mt_bin = 3;}

  //printf("%p, %p, %p, %p\n", fList_P02, fList_P13, fList_P03, fList_P12);
  printf("which_mt_bin = %i\n",which_mt_bin);
  TH2F* h_P02 = NULL;
  TH2F* h_P13 = NULL;
  TH2F* h_P03 = NULL;
  TH2F* h_P12 = NULL;

  if(mt_bin!=1000){
    h_P02 = (TH2F*)fList_P02->FindObject(TString::Format("MEmT_%i_pT_PionNucleon0_2_vs_kStar", which_mt_bin));
    h_P13 = (TH2F*)fList_P13->FindObject(TString::Format("MEmT_%i_pT_PionNucleon1_3_vs_kStar", which_mt_bin));
    h_P03 = (TH2F*)fList_P03->FindObject(TString::Format("MEmT_%i_pT_PionNucleon0_3_vs_kStar", which_mt_bin));
    h_P12 = (TH2F*)fList_P12->FindObject(TString::Format("MEmT_%i_pT_PionNucleon1_2_vs_kStar", which_mt_bin));
  }

  double mean_pNucl_pT_below300=0;
  TH2F* h_pT_Nucl = NULL;
  //pim_d
  if(mt_bin==1000){

  }
  else if(abs(mt_bin)>=100){
    h_pT_Nucl = (TH2F*)h_P03->Clone("h_pT_Nucl");
    h_pT_Nucl->Add(h_P12);
    //printf("A\n");
  }
  else{
    h_pT_Nucl = (TH2F*)h_P02->Clone("h_pT_Nucl");
    h_pT_Nucl->Add(h_P13);
    //printf("B\n");
  }

  TH1F* hProj = NULL;
  if(mt_bin!=1000){
    hProj = (TH1F*)h_pT_Nucl->ProjectionY(TString::Format("h_pT_Nucl_Proj"),h_pT_Nucl->GetXaxis()->FindBin(0.5*0.001),h_pT_Nucl->GetXaxis()->FindBin(299.5*0.001));
    mean_pNucl_pT_below300 = hProj->GetMean()*1000.;
  }
  //TCanvas c1("can1", "can1", 1);
  //hResolution_pd->Draw("colz");
  ////hProj->Draw();
  //c1.SaveAs("/home/dimihayl/Desktop/hProj.png");
  //hProj->Draw();
  //usleep(10e6);
  //printf("mean_pNucl_pT_below300 = %.3f (%i %i)\n",mean_pNucl_pT_below300,h_pT_Nucl->GetXaxis()->FindBin(0.5*0.001),h_pT_Nucl->GetXaxis()->FindBin(299.5*0.001));
  delete hProj;

  //printf("%p %p %p\n",fList1,fList2,hME);
  //printf("hello\n");

  TGraphErrors  gData(max_mom_bins);
  gData.SetName("gData");
  gData.SetLineColor(kBlack);
  gData.SetLineWidth(4);

  TGraph gFit(max_mom_bins);
  gFit.SetName("gFit");
  gFit.SetLineColor(kPink-9);
  gFit.SetLineWidth(4);

  TGraph gBaseline(max_mom_bins);
  gBaseline.SetName("gBaseline");
  gBaseline.SetLineColor(kGreen-5);
  gBaseline.SetLineWidth(4);

  TGraph gFemto(max_mom_bins);
  gFemto.SetName("gFemto");
  gFemto.SetLineColor(kBlue+2);
  gFemto.SetLineWidth(4);

  TGraph gDelta(max_mom_bins);
  gDelta.SetName("gDelta");
  gDelta.SetLineColor(kOrange-1);
  gDelta.SetLineWidth(4);

  TGraph gMadeUpReso(max_mom_bins);
  gMadeUpReso.SetName("gMadeUpReso");
  gMadeUpReso.SetLineColor(kAzure+10);
  gMadeUpReso.SetLineWidth(4); 

  float MT_LOW;
  float MT_UP;
  float RADIUS;
  int DATA_TYPE;
  int FIT_TYPE;
  int MT_ID = mt_bin;
  float FIT_RANGE;
  float LAM_GEN;
  float FRAC_D;
  float NUM_DELTAS;
  float NUM_MADEUP_RESO;
  float AMP_DELTA;
  float DELTA_MASS;
  float DELTA_WIDTH;
  float PS_PT;
  float PS_TEMP;
  float CHI2_500;
  float PT_SCALE;

  if(mt_bin==1000){
    MT_LOW = 0;
    MT_UP = 10000;  
  }
  else if(mt_bin==-1){
    MT_LOW = mT_low.at(0).at(0);
    MT_UP = mT_up.at(0).at(mT_up.at(0).size()-1);
  }
  else if(mt_bin==-101){
    MT_LOW = mT_low.at(1).at(0);
    MT_UP = mT_up.at(1).at(mT_up.at(0).size()-1);
  }
  else if(mt_bin%100>=0 && mt_bin%100<=4){
    MT_LOW = mT_low.at(mt_bin/100).at(mt_bin%100);
    MT_UP = mT_up.at(mt_bin/100).at(mt_bin%100);    
  }
  else{
    printf("Silly mT bin\n");
    abort();
  }

  //int SEED, int mt_bin, int NumIter, bool Bootstrap=true, bool DataVar=true, bool FitVar=true
  TFile fOutputFile(TString::Format("%s/pi_d/FitMar2025_Files/Results_v2/%s_mT%i_B%i_DV%i_FV%i_S%i.root",
    GetCernBoxDimi(), Description.Data(), mt_bin, Bootstrap, DataVar, FitVar, SEED), "recreate");
  TTree* pi_d_Tree = new TTree("pi_d_Tree","pi_d_Tree");
  //pi_d_Tree->Branch("seed",&SEED,"seed/I")
  pi_d_Tree->Branch("gData","TGraphErrors",&gData,32000,0);//
  pi_d_Tree->Branch("gFit","TGraph",&gFit,32000,0);//
  pi_d_Tree->Branch("gBaseline","TGraph",&gBaseline,32000,0);//
  pi_d_Tree->Branch("gFemto","TGraph",&gFemto,32000,0);//
  pi_d_Tree->Branch("gDelta","TGraph",&gDelta,32000,0);//
  pi_d_Tree->Branch("gMadeUpReso","TGraph",&gMadeUpReso,32000,0);//
  pi_d_Tree->Branch("mT_low",&MT_LOW,"mT_low/F");//
  pi_d_Tree->Branch("mT_up",&MT_UP,"mT_up/F");//
  pi_d_Tree->Branch("mT_id",&MT_ID,"mT_id/I");//
  pi_d_Tree->Branch("rad",&RADIUS,"rad/F");//
  pi_d_Tree->Branch("data_type",&DATA_TYPE,"data_type/I");//ABC, A=0 => no boot, A=1 = bootstrap, BC = variation ID. If negative => default
  pi_d_Tree->Branch("fit_type",&FIT_TYPE,"fit_type/I");
  pi_d_Tree->Branch("fit_range",&FIT_RANGE,"fit_range/F");//
  pi_d_Tree->Branch("lam_gen",&LAM_GEN,"lam_gen/F");//
  pi_d_Tree->Branch("frac_D",&FRAC_D,"frac_D/F");//
  pi_d_Tree->Branch("num_deltas",&NUM_DELTAS,"num_deltas/F");//
  pi_d_Tree->Branch("num_madeup_reso",&NUM_MADEUP_RESO,"num_madeup_reso/F");//
  pi_d_Tree->Branch("amp_delta",&AMP_DELTA,"amp_delta/F");//
  pi_d_Tree->Branch("CkCutOff",&CKCUTOFF,"CkCutOff/F");//
  pi_d_Tree->Branch("delta_mass",&DELTA_MASS,"delta_mass/F");//
  pi_d_Tree->Branch("delta_width",&DELTA_WIDTH,"delta_width/F");//
  pi_d_Tree->Branch("ps_pT",&PS_PT,"ps_pT/F");//
  pi_d_Tree->Branch("ps_Temp",&PS_TEMP,"ps_Temp/F");//
  pi_d_Tree->Branch("pT_scale",&PT_SCALE,"pT_scale/F");//
  pi_d_Tree->Branch("chi2_500",&CHI2_500,"chi2_500/F");//

/*
  ppTree->Branch("gBl","TGraph",&GBL,32000,0);//
  ppTree->Branch("gFemto","TGraph",&GFEMTO,32000,0);//
  ppTree->Branch("BaselineVar",&BASELINEVAR,"BaselineVar/I");//
  ppTree->Branch("FemtoRegion",&FEMTOREGION,"FemtoRegion/F");//
  ppTree->Branch("pS0",&PS0,"pS0/I");//
  ppTree->Branch("SmearStrategy",&SMEARSTRATEGY,"SmearStrategy/I");//
  ppTree->Branch("pp_lam_var",&PP_LAM_VAR,"pp_lam_var/I");//
  ppTree->Branch("MomSmearVar","TString",&MOMSMEARVAR,8000,0);//
  ppTree->Branch("CkCutOff",&CKCUTOFF,"CkCutOff/F");//
  ppTree->Branch("pL_pot_var",&PL_POT_VAR,"pL_pot_var/I");//
  ppTree->Branch("SourceVar","TString",&SOURCEVAR,8000,0);//
  ppTree->Branch("lam_pp",&LAM_PP,"lam_pp/F");//
  ppTree->Branch("lam_ppL",&LAM_PPL,"lam_ppL/F");//
  ppTree->Branch("nsig",&NSIG,"nsig/F");//
  ppTree->Branch("rad",&RADIUS,"rad/F");//
  ppTree->Branch("raderr",&RADERR,"raderr/F");//
*/
  unsigned NumMomBins = 20;
  double kCatMin = 0;
  double kCatMax = 400;
  CATS Kitty;
  Kitty.SetMomBins(NumMomBins, kCatMin, kCatMax);
  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());
  if(mt_bin==1000){
    AnalysisObject.SetUpCats_Kd(Kitty, "", "Gauss", 0, 0);
    Kitty.SetQ1Q2(-1);
  }
  if(fabs(mt_bin)>=100){
    AnalysisObject.SetUpCats_pi_d(Kitty, "Gauss", "Gauss", -1, 0);
    Kitty.SetQ1Q2(-1);
  }
  else{
    AnalysisObject.SetUpCats_pi_d(Kitty, "", "Gauss", 0, 0);
    Kitty.SetQ1Q2(+1);    
  }
  Kitty.SetAnaSource(0, eff_source_radii.at(0));
  if(Kitty.GetNumSourcePars()>1) Kitty.SetAnaSource(1, 2);
  Kitty.KillTheCat();
  Kitty.SetNotifications(CATS::nWarning);

  DLM_Ck CkKitty_Default(Kitty.GetNumSourcePars(),0,Kitty,max_mom_bins,0,max_kstar);
  DLM_Ck CkKitty_Wioleta(max_mom_bins,0,max_kstar);
  //DLM_Ck CkKitty(Kitty.GetNumSourcePars(),0,Kitty,max_mom_bins,0,max_kstar);
  if(mt_bin==1000){
    //  TFile fOutputFile(TString::Format("%s/pi_d/FitNov2024_Files/Results/%s_mT%i_B%i_DV%i_FV%i_S%i.root",

    TFile fWioleta(TString::Format("%s/pi_d/Kd/LL/K-dWioletaAN_1p351.root",GetCernBoxDimi()));
    TGraph* gCgen = (TGraph*)fWioleta.Get("gCgen");
    for(unsigned uMom=0; uMom<max_mom_bins; uMom++){
      double MOM = CkKitty_Wioleta.GetBinCenter(0, uMom);
      CkKitty_Wioleta.SetBinContent(uMom, gCgen->Eval(MOM));
      //printf("%f %f\n",MOM, gCgen->Eval(MOM));
    }
    fWioleta.Close();
  }
  else{//the default stuff
    CkKitty_Default.SetSourcePar(0,Kitty.GetAnaSourcePar(0));
    if(Kitty.GetNumSourcePars()>1) CkKitty_Default.SetSourcePar(1,Kitty.GetAnaSourcePar(1));
  }
  DLM_Ck& CkKitty = mt_bin==1000?CkKitty_Wioleta:CkKitty_Default;
  
  CkKitty.SetCutOff(CKCUTOFF,CkConv);
  CkKitty.Update();

  DLM_CkDecomposition CkDecKitty("pi_d",0,CkKitty,hResolution_pd);
  //DLM_CkDecomposition CkDecKitty("pi_d",0,CkKitty,NULL);
  CkDecKitty.Update(true, true);
  Bulgaristan_CkDec = &CkDecKitty;


  double Progress = 0;
  bool data_saved = false;
  for(int uIter=0; uIter<NumIter; uIter++){
    printf("uIter = %u\n",uIter);
    TFile* fInputData = NULL;
    TH1F* hData = NULL;
    TDirectoryFile* fDir = NULL;
    TList* fList1 = NULL;
    TList* fList2 = NULL;
    TH1F* hME = NULL;
    TH1F* hSE = NULL;

    //-1 is default, -2 is default without bootstrap
    int iDataVar = -1;
    do{
      if(DataVar==true){
        //if we get a -1, it is the default
        iDataVar = rangen.Integer(NumDataVars+1)-1;
      }
      if(!data_saved){
        iDataVar = -2;
        data_saved = true;
      }


      //mt int pip_d used for aprovals at CF and PF
      if(mt_bin == -1){
        if(iDataVar==-1 || iDataVar==-2){
          InputDataFileName = TString::Format("%s/pi_d/FitOct2024_Files/2024-10-29/20MeVBins_pip_d/CF_Piondeuteron_Var0.root",GetCernBoxDimi());
          DefaultHistoName = "hCk_ReweightedPidVar0MeV_0";
          MeListName1 = "PairDist";
          MeListName2 = "PairReweighted";
          MeHistName = "hTotalME _Rebinned_5_Reweighted";
          SeHistName = "hTotalSE _Rebinned_5_Reweighted";
        }
        else{
          InputDataFileName = TString::Format("%s/pi_d/FitOct2024_Files/2024-10-29/20MeVBins_pip_d/CF_Pid_Var%i.root",GetCernBoxDimi(), iDataVar+1);
          DefaultHistoName = TString::Format("hCk_ReweightedReweightedPidVar%iMeV_1MeV_0", iDataVar+1);
          MeListName1 = "PairDist";
          MeListName2 = "PairReweighted";
          MeHistName = TString::Format("hTotalME_PidVar%i_Rebinned_5_Reweighted", iDataVar+1);
          SeHistName = TString::Format("hTotalSE_PidVar%i_Rebinned_5_Reweighted", iDataVar+1);
        }
        std::ifstream file(InputDataFileName.Data());
        if(file.good()) fInputData = new TFile(InputDataFileName, "read");
        if(fInputData) hData = (TH1F*)fInputData->Get(DefaultHistoName);

        if(fInputData) fList1 = (TList*)(fInputData->FindObjectAny(MeListName1));
        if(fList1) fList2 = (TList*)fList1->FindObject(MeListName2);
        if(fList2) hME = (TH1F*)fList2->FindObject(MeHistName);
        if(hME){
          hME->GetXaxis()->SetLimits(hME->GetXaxis()->GetXmin()*1000.,hME->GetXaxis()->GetXmax()*1000.);
        }
        if(fList2) hSE = (TH1F*)fList2->FindObject(SeHistName);
        if(hSE){
          hSE->GetXaxis()->SetLimits(hSE->GetXaxis()->GetXmin()*1000.,hSE->GetXaxis()->GetXmax()*1000.);
        }

      }
      //mt int pim_d used for approvals at CF and PF
      else if(mt_bin == -101){
        if(iDataVar==-1 || iDataVar==-2){
          InputDataFileName = TString::Format("%s/pi_d/FitOct2024_Files/2024-10-29/mTIntegrated_pim_d/CF_AntiPiondeuteron_Var0.root",GetCernBoxDimi());
          DefaultHistoName = "hCk_ReweightedPidVar0MeV_0";
          MeListName1 = "PairDist";
          MeListName2 = "PairReweighted";
          MeHistName = "hTotalME _Rebinned_5_Reweighted";
          SeHistName = "hTotalSE _Rebinned_5_Reweighted";
        }
        else{
          InputDataFileName = TString::Format("%s/pi_d/FitOct2024_Files/2024-10-29/mTIntegrated_pim_d/CF_AntiPid_Var%i.root",GetCernBoxDimi(), iDataVar+1);
          DefaultHistoName = TString::Format("hCk_ReweightedReweightedPidVar%iMeV_1MeV_0", iDataVar+1);
          MeListName1 = "PairDist";
          MeListName2 = "PairReweighted";
          MeHistName = TString::Format("hTotalME_PidVar%i_Rebinned_5_Reweighted", iDataVar+1);    
          SeHistName = TString::Format("hTotalSE_PidVar%i_Rebinned_5_Reweighted", iDataVar+1);      
        }
        std::ifstream file(InputDataFileName.Data());
        if(file.good()) fInputData = new TFile(InputDataFileName, "read");
        if(fInputData) hData = (TH1F*)fInputData->Get(DefaultHistoName);

        if(fInputData) fList1 = (TList*)(fInputData->FindObjectAny(MeListName1));
        if(fList1) fList2 = (TList*)fList1->FindObject(MeListName2);
        if(fList2) hME = (TH1F*)fList2->FindObject(MeHistName);
        if(hME){
          hME->GetXaxis()->SetLimits(hME->GetXaxis()->GetXmin()*1000.,hME->GetXaxis()->GetXmax()*1000.);
        }
        if(fList2) hSE = (TH1F*)fList2->FindObject(SeHistName);
        if(hSE){
          hSE->GetXaxis()->SetLimits(hSE->GetXaxis()->GetXmin()*1000.,hSE->GetXaxis()->GetXmax()*1000.);
        }
      }
      //mt diff pip_d used for approvals at CF
      else if(mt_bin>=0 && mt_bin<100 && false){
        InputDataFileName = TString::Format("%s/pi_d/FitOct2024_Files/mT_2024-10-11/mT06102024Signal_sum_Mt_%i.root",GetCernBoxDimi(),mt_bin+1);
        DefaultHistoName = "hCk_ReweightedMeV_2";
        if(DataVar==true){
          printf("ERROR: At the moment we do not have mt variations\n");
          abort();
        }
        //VarDirName = "Raw";
        MeListName1 = "PairDist";
        MeListName2 = "PairReweighted";
        MeHistName = "hTotalMEpairs _Rebinned_5_Reweighted";
        SeHistName = "hTotalMEpairs _Rebinned_5_Reweighted";//?
        fInputData = new TFile(InputDataFileName, "read");
        if(fInputData) hData = (TH1F*)fInputData->Get(DefaultHistoName);
        //printf("fInputData = %p\n",fInputData);
        //fDir = (TDirectoryFile*)(fInputData->FindObjectAny(VarDirName));
        //printf("fDir = %p\n",fDir);
        //fDir->GetObject(VarHistoBaseName+TString::Format("%i",iDataVar),hData);
        //printf("hData = %p\n",hData);
        if(fInputData) fList1 = (TList*)(fInputData->FindObjectAny(MeListName1));
        //printf("fList1 = %p\n",fList1);
        if(fList1) fList2 = (TList*)fList1->FindObject(MeListName2);
        if(fList2) hME = (TH1F*)fList2->FindObject(MeHistName);
        if(hME){
          hME->GetXaxis()->SetLimits(hME->GetXaxis()->GetXmin()*1000.,hME->GetXaxis()->GetXmax()*1000.);
        }
        if(fList2) hSE = (TH1F*)fList2->FindObject(SeHistName);
        if(hSE){
          hSE->GetXaxis()->SetLimits(hSE->GetXaxis()->GetXmin()*1000.,hSE->GetXaxis()->GetXmax()*1000.);
        }
      }
      //mt diff pim_d used for approvals at CF
      else if(mt_bin>=100 && mt_bin<200 && false){
        InputDataFileName = TString::Format("%s/pi_d/FitOct2024_Files/AntiPionDeuteron/22102024FilesSignal_antisum_Mt_%i.root",GetCernBoxDimi(),mt_bin%100+1);
        DefaultHistoName = "hCk_ReweightedMeV_2";
        if(DataVar==true){
          printf("ERROR: At the moment we do not have mt variations\n");
          abort();
        }
        VarDirName = "Raw";
        MeListName1 = "PairDist";
        MeListName2 = "PairReweighted";
        MeHistName = "hTotalMEapairs _Rebinned_5_Reweighted";
        fInputData = new TFile(InputDataFileName, "read");
        if(fInputData) hData = (TH1F*)fInputData->Get(DefaultHistoName);
        //fDir = (TDirectoryFile*)(fInputData->FindObjectAny(VarDirName));
        //fDir->GetObject(VarHistoBaseName+TString::Format("%i",iDataVar),hData);
        if(fInputData) fList1 = (TList*)(fInputData->FindObjectAny(MeListName1));
        if(fList1) fList2 = (TList*)fList1->FindObject(MeListName2);
        if(fList2) hME = (TH1F*)fList2->FindObject(MeHistName);
        if(hME){
          hME->GetXaxis()->SetLimits(hME->GetXaxis()->GetXmin()*1000.,hME->GetXaxis()->GetXmax()*1000.);
        }
        if(fList2) hSE = (TH1F*)fList2->FindObject(SeHistName);
        if(hSE){
          hSE->GetXaxis()->SetLimits(hSE->GetXaxis()->GetXmin()*1000.,hSE->GetXaxis()->GetXmax()*1000.);
        }
      }
      //mt diff pip_d used for approvals at PF
      else if(mt_bin>=0 && mt_bin<100 && true){
        if(iDataVar==-1 || iDataVar==-2){
          InputDataFileName = TString::Format("%s/pi_d/FitNov2024_Files/Files/PionDeuteron/mT_differential/mTBin%i/CF_Piondeuteron_Var0.root",GetCernBoxDimi(),abs(mt_bin)%100+1);
          DefaultHistoName = "hCk_ReweightedPidVar0MeV_0";
          MeListName1 = "PairDist";
          MeListName2 = "PairReweighted";
          MeHistName = "hTotalME _Rebinned_5_Reweighted";
        }
        else{
          InputDataFileName = TString::Format("%s/pi_d/FitNov2024_Files/Files/PionDeuteron/mT_differential/mTBin%i/CF_Pid_Var%i.root",
            GetCernBoxDimi(),abs(mt_bin)%100+1,iDataVar+1);
          DefaultHistoName = TString::Format("hCk_ReweightedReweightedPidVar%iMeV_1MeV_0", iDataVar+1);
          MeListName1 = "PairDist";
          MeListName2 = "PairReweighted";
          MeHistName = TString::Format("hTotalME_PidVar%i_Rebinned_5_Reweighted", iDataVar+1);
        }
        std::ifstream file(InputDataFileName.Data());
        if(file.good()) fInputData = new TFile(InputDataFileName, "read");
        //else printf("bad file %s\n", InputDataFileName.Data());
        if(fInputData) hData = (TH1F*)fInputData->Get(DefaultHistoName);

        if(fInputData) fList1 = (TList*)(fInputData->FindObjectAny(MeListName1));
        if(fList1) fList2 = (TList*)fList1->FindObject(MeListName2);
        if(fList2) hME = (TH1F*)fList2->FindObject(MeHistName);
        if(hME){
          hME->GetXaxis()->SetLimits(hME->GetXaxis()->GetXmin()*1000.,hME->GetXaxis()->GetXmax()*1000.);
        }
        if(fList2) hSE = (TH1F*)fList2->FindObject(SeHistName);
        if(hSE){
          hSE->GetXaxis()->SetLimits(hSE->GetXaxis()->GetXmin()*1000.,hSE->GetXaxis()->GetXmax()*1000.);
        }
//printf("hi there %p %p %s\n",hSE,hME,SeHistName.Data());
//fList2->Print();

      } 
      //mt diff pim_d used for approvals at PF
      else if(mt_bin>=100 && mt_bin<200 && true){

        if(iDataVar==-1 || iDataVar==-2){
          InputDataFileName = TString::Format("%s/pi_d/FitNov2024_Files/Files/AntiPionDeuteron/mT_differential/mTBin%i/CF_AntiPiondeuteron_Var0.root",GetCernBoxDimi(),abs(mt_bin)%100+1);
          DefaultHistoName = "hCk_ReweightedPidVar0MeV_0";
          MeListName1 = "PairDist";
          MeListName2 = "PairReweighted";
          MeHistName = "hTotalME _Rebinned_5_Reweighted";
        }
        else{
          InputDataFileName = TString::Format("%s/pi_d/FitNov2024_Files/Files/AntiPionDeuteron/mT_differential/mTBin%i/CF_AntiPid_Var%i.root",
            GetCernBoxDimi(),abs(mt_bin)%100+1,iDataVar+1);
          DefaultHistoName = TString::Format("hCk_ReweightedReweightedPidVar%iMeV_1MeV_0", iDataVar+1);
          MeListName1 = "PairDist";
          MeListName2 = "PairReweighted";
          MeHistName = TString::Format("hTotalME_PidVar%i_Rebinned_5_Reweighted", iDataVar+1);
        }
        std::ifstream file(InputDataFileName.Data());
        if(file.good()) fInputData = new TFile(InputDataFileName, "read");
        //else printf("bad file %s\n", InputDataFileName.Data());
        if(fInputData) hData = (TH1F*)fInputData->Get(DefaultHistoName);

        if(fInputData) fList1 = (TList*)(fInputData->FindObjectAny(MeListName1));
        if(fList1) fList2 = (TList*)fList1->FindObject(MeListName2);
        if(fList2) hME = (TH1F*)fList2->FindObject(MeHistName);
        if(hME){
          hME->GetXaxis()->SetLimits(hME->GetXaxis()->GetXmin()*1000.,hME->GetXaxis()->GetXmax()*1000.);
        }
        if(fList2) hSE = (TH1F*)fList2->FindObject(SeHistName);
        if(hSE){
          hSE->GetXaxis()->SetLimits(hSE->GetXaxis()->GetXmin()*1000.,hSE->GetXaxis()->GetXmax()*1000.);
        }

      }
      //this is the ULTRA-GHETTO Kminus-d fit, demanded by Laura on 14th Nov
      else if(mt_bin==1000){
        //printf("Kd\n");
        InputDataFileName = TString::Format("%s/pi_d/Kd/outKD_mydeuteronsOpenPIDDCA_std.root",GetCernBoxDimi());
        std::ifstream file(InputDataFileName.Data());
        if(file.good()) fInputData = new TFile(InputDataFileName, "read");
        //fInputData->ls();
        TGraphErrors* gData_Kd = (TGraphErrors*)fInputData->Get("gCF");
        //printf("%p\n",gData_Kd);
        gROOT->cd();
        if(hData) delete hData;
        hData = new TH1F("hData", "hData", 40, 0, 800);
        for(unsigned uMom=0; uMom<hData->GetNbinsX(); uMom++){
          double xval, yval;
          gData_Kd->GetPoint(uMom, xval, yval);
          //printf("%f %f\n",xval, yval);
          hData->SetBinContent(uMom+1, yval);
          hData->SetBinError(uMom+1, gData_Kd->GetErrorY(uMom));
          //printf("%f %f %f\n",xval,yval,gData_Kd->GetErrorY(uMom));
        }
        hME = (TH1F*)fInputData->Get("ME");
        hSE = (TH1F*)fInputData->Get("SE");
      }

        if(!hME){
          //printf("fuck %i! %p (%s) %p %p\n", iDataVar, fInputData, InputDataFileName.Data(), fList1, fList2);
          //if(fList2){
          //  fList2->ls();
          //  printf(" ME_name = %s\n",MeHistName.Data());
          //}
        }
      //tyj


    }
    //while(hData==NULL && iDataVar!=-1);//I dont know why this was like this
    while(hData==NULL);
    //printf("GOING FORWARD\n");

    //TDirectoryFile* fDir = NULL;
    //TH1F* hData = NULL;
    //TList* fList1 = NULL;
    //TList* fList2 = NULL;
    //TH1D* hME = NULL;
    //fList1 = (TList*)(fInputData.FindObjectAny(MeListName1));
    //if(fList1) fList2 = (TList*)fList1->FindObject(MeListName2);
    ////fList2->ls();
    //if(fList2) hME = (TH1D*)fList2->FindObject(MeHistName);
    //if(hME){
    //  hME->GetXaxis()->SetLimits(hME->GetXaxis()->GetXmin()*1000.,hME->GetXaxis()->GetXmax()*1000.);
    //}

/*
    fInputData.cd();
    //printf("iDataVar=%i\n",iDataVar);
    //fInputData.ls();
    if(iDataVar==-1){
      hData = (TH1F*)fInputData.Get(DefaultHistoName);
    }
    else{
      fDir = (TDirectoryFile*)(fInputData.FindObjectAny(VarDirName));
      fDir->GetObject(VarHistoBaseName+TString::Format("%i",iDataVar),hData);
    }
*/

    if(!hData){
      printf("!hData\n");
      abort();
    }

    gROOT->cd();
    TH1F* hDataToFit = (TH1F*)hData->Clone("hDataToFit");;
    if(Bootstrap==true && iDataVar!=-2){
      double MOM = 0;
      unsigned uMom = 0;
      while(MOM<max_kstar){
        MOM = hData->GetBinCenter(uMom+1);
        double new_value;
        do{
          new_value = rangen.Gaus(hData->GetBinContent(uMom+1), hData->GetBinError(uMom+1));
        }
        while(new_value<0);
        //printf("%.3f +/- %.3f --> %.3f\n",hData->GetBinContent(uMom+1), hData->GetBinError(uMom+1), new_value);
        hDataToFit->SetBinContent(uMom+1, new_value);
        uMom++;
      }
    }

    int pt_scale_int = 0;
    FIT_RANGE = fit_range.at(0);
    RADIUS = eff_source_radii.at(0);
    PT_SCALE = pT_scale.at(pt_scale_int);
    FIT_TYPE = fit_types.at(0);
    if(mt_bin<0) LAM_GEN = lambda_gen_mTint;
    else if(mt_bin==1000) LAM_GEN = 0.9;
    else LAM_GEN = lambda_gen_avg.at(mt_bin%100);
    AMP_DELTA = amplitude_delta.at(0);
    if(FitVar){
      int rndint = rangen.Integer(fit_range.size());
      FIT_RANGE = fit_range.at(rndint);
      //printf("%f %i %i\n", FIT_RANGE, rndint, fit_range.size());
      RADIUS = eff_source_radii.at(rangen.Integer(eff_source_radii.size()));
      if(LAM_GEN>=0){
        if(mt_bin<0) LAM_GEN = lambda_gen_mTint*lambda_gen_relvar.at(rangen.Integer(lambda_gen_relvar.size()));
        else if(mt_bin==1000) LAM_GEN = 0.9;
        else LAM_GEN = lambda_gen_avg.at(mt_bin%100)*lambda_gen_relvar.at(rangen.Integer(lambda_gen_relvar.size()));
      }
      if(AMP_DELTA>=0) AMP_DELTA = amplitude_delta.at(rangen.Integer(amplitude_delta.size()));

      rndint = rangen.Integer(pT_scale.size());
      PT_SCALE = pT_scale.at(rndint);


      rndint = rangen.Integer(fit_types.size());
      FIT_TYPE = fit_types.at(rndint);
    }
    if(mt_bin==1000){
      RADIUS = 1.35;
    }

    TFile fConversion(InputConversionFileName.at(pt_scale_int), "read");
    Kstar_Modifier = (TGraph*)fConversion.Get("gR_ppi_dpi_c");


    Bulgaristan_CkDec->AddPhaseSpace(hME);



    DATA_TYPE = abs(iDataVar);
    if(Bootstrap==true) DATA_TYPE += 100;
    if(iDataVar<0) DATA_TYPE = -DATA_TYPE;
    //overrides the other stuff, this is default data with no boot
    if(iDataVar==-2) DATA_TYPE = -2;


    //par[0] = source size
    //par[1] = alpha par
    //par[2] = lambda_genuine
    //par[3] = lambda_d
    //par[4] = delta_amplitude
    //par[5] = CkCutOff
    //par[6] = CkConv
    //par[7] = just in case
    //par[8] = Delta normalization, should be a dummy fixed to -1e6, giving instructions to the code to renorm
    //par[9] = mass of daughter 1 (say the pion)
    //par[10] = mass of daughter 2 (say the proton)
    //par[11] = mass of the delta
    //par[12] = width of the delta
    //par[13] = avg pT of the daughters
    //par[14] = effective temperature
    //par[15] = just in case
    //par[16] = norm
    //par[17] = 0 to castrate the pol3
    //par[18] = position of the max of the pol3
    //par[19] = p3 parameter of the pol3
    //par[20] = -1e6 to switch off the pol4
    //[21-25] parameters of the extra resonances added during the review at Nature, see Bulgaristan_fit for details 
    const int NumFitPar = 26;
    TF1* fData = new TF1("fData",Bulgaristan_fit,0,FIT_RANGE,NumFitPar);
    TF1* fDataDummy = new TF1("fDataDummy",Bulgaristan_fit,0,FIT_RANGE,NumFitPar);


    //hDataToFit->Fit(fData,"Q, S, N, R, M");

    //minimum is one, if you whant to refine the FRAC_D recursively, increase it
    int RepeatFit = 3;
    int iRepeat = 0;
    FRAC_D = 1.5e-2;
    do{


      fData->FixParameter(0, RADIUS);
      fData->FixParameter(1, 2);
      if(LAM_GEN>=0) fData->FixParameter(2, LAM_GEN);
      else{
        if(mt_bin<0){
          fData->SetParameter(2, fabs(lambda_gen_mTint));
          fData->SetParLimits(2, fabs(lambda_gen_mTint*lambda_gen_relvar.at(0)), fabs(lambda_gen_mTint*lambda_gen_relvar.at(1)));        
        }
        else if(mt_bin==1000){
          fData->FixParameter(2, LAM_GEN);
        }
        else{
          fData->SetParameter(2, fabs(lambda_gen_avg.at(mt_bin%100)));
          fData->SetParLimits(2, fabs(lambda_gen_avg.at(mt_bin%100)*lambda_gen_relvar.at(0)), fabs(lambda_gen_avg.at(mt_bin%100)*lambda_gen_relvar.at(1)));        
        }
      }
      //fData->SetParameter(3,0.5);
      //fData->SetParLimits(3,0,1);
      //3 will be fixed later on (do-while)
      //fData->FixParameter(3,1.5e-02);
      if(mt_bin==1000){
        fData->SetParameter(4, 1);
        fData->SetParLimits(4, 1e-4, 10);        
      }
      else if(AMP_DELTA>=0) fData->FixParameter(4, AMP_DELTA);
      else{
        fData->SetParameter(4, fabs(amplitude_delta.at(0)));
        fData->SetParLimits(4, fabs(amplitude_delta.at(1)), fabs(amplitude_delta.at(2)));
      }
      fData->FixParameter(5, CKCUTOFF);
      fData->FixParameter(6, CkConv);

      fData->FixParameter(7, FIT_TYPE);

      fData->FixParameter(8, -1e6);
      fData->FixParameter(9, Mass_pic);
      fData->FixParameter(10, Mass_p);
      fData->SetParameter(11,1232);
      fData->SetParLimits(11,1180,1260);
  //fData->FixParameter(11,1232); 
  fData->FixParameter(11,1215);
  if(mt_bin==1000){
    fData->SetParameter(11,275);
    fData->SetParLimits(11,250,300);
  }
      fData->SetParameter(12,100);
      fData->SetParLimits(12,40,200);
  //fData->FixParameter(12,90); 
      //this can be perhaps fixed. Marcel does it, but based on the mT
      /*
      mT Bin: 1
      Average mT in GeV Delta: 0.675343
      Average kT in GeV Delta: 0.407003
      Average pT in MeV Delta: 814.007
      mT Bin: 2
      Average mT in GeV Delta: 0.838604
      Average kT in GeV Delta: 0.642511
      Average pT in MeV Delta: 1285.02
      mT Bin: 3
      Average mT in GeV Delta: 1.0609
      Average kT in GeV Delta: 0.913824
      Average pT in MeV Delta: 1827.65
      mT Bin: 4
      Average mT in GeV Delta: 1.33761
      Average kT in GeV Delta: 1.22424
      Average pT in MeV Delta: 2448.49
      mT Bin: 5
      Average mT in GeV Delta: 1.69694
      Average kT in GeV Delta: 1.60909
      Average pT in MeV Delta: 3218.18
      mT Bin: 6
      Average mT in GeV Delta: 2.14387
      Average kT in GeV Delta: 2.07503
      Average pT in MeV Delta: 4150.05
      */

      double kT_Delta;
      //kT_Delta = 2.*kay_tee*PT_SCALE; //estimate from the mT converted
      //fData->SetParameter(13,400);
      //fData->SetParLimits(13,150,1000);
      //printf("kT_Delta = %f\n",kT_Delta);
      //printf("kay_tee = %f\n",kay_tee);
      //printf("PT_SCALE = %f\n",PT_SCALE);
      
      //if(mt_bin<0){
      //  kT_Delta = piN_pT_avg;
      //}
      //else{
      //  kT_Delta = piN_pT.at(mt_bin%100);
      //}
      //fData->FixParameter(13, kT_Delta);

      fData->SetParameter(13,mean_pNucl_pT_below300);
      fData->SetParLimits(13,mean_pNucl_pT_below300*0.8,mean_pNucl_pT_below300*1.2);

      //fData->SetParameter(13,kT_Delta);
      //fData->SetParLimits(13,kT_Delta*0.7,kT_Delta*1.3);
      fData->SetParameter(14,15);
      fData->SetParLimits(14,5,60);

  //fix temp and gamma based on marcel
      if(mt_bin<0){
        fData->FixParameter(12,marcel_avg_gamma);//gamma
        fData->SetParameter(14,marcel_avg_temp);//temp
      }
      else if(mt_bin==1000){
        fData->SetParameter(12,15.6);//gamma
        fData->SetParLimits(12,13.1,20.6);
        fData->FixParameter(14,0);//temp
      }
      else{
        fData->FixParameter(12,marcel_gamma_Dpp.at(mt_bin%100)*Dpp_weight + marcel_gamma_D0.at(mt_bin%100)*(1.-Dpp_weight));//gamma
        fData->SetParameter(14,marcel_temp_Dpp.at(mt_bin%100)*Dpp_weight + marcel_temp_D0.at(mt_bin%100)*(1.-Dpp_weight));//temp
      }

      fData->FixParameter(15, 0);


      if(mt_bin==1000){
        fData->SetParameter(16,1);
        fData->SetParLimits(16,0.9,1.1);
        fData->SetParameter(17, 300);
        fData->SetParLimits(17, 150,450);
        fData->SetParameter(18,2000);
        fData->SetParLimits(18,450,50000);
        fData->SetParameter(19,0);
        fData->SetParLimits(19,-1e-8,1e-8);
        fData->FixParameter(20, 1);
      }
      else{
        fData->SetParameter(16,1);
        fData->SetParLimits(16,0.9,1.1);
        fData->FixParameter(17, 0);
        fData->SetParameter(18,400);
        fData->SetParLimits(18,1,100000);
        fData->SetParameter(19,0);
        fData->SetParLimits(19,-1e-8,1e-8);
        fData->FixParameter(20, -1e6);
      }

      if(FIT_TYPE<10){
        fData->FixParameter(21, 0);
        fData->FixParameter(22, 1);
        fData->FixParameter(23, 1);
        fData->FixParameter(24, 1);
        fData->FixParameter(25, 1);
      }
      else{
        //fData->FixParameter(21, 5*1.07719e-02);
        fData->SetParameter(21, 1e-4);
        fData->SetParLimits(21, 0,10);
        fData->FixParameter(22, Mass_d);
        fData->FixParameter(23, Mass_pic);
        fData->FixParameter(24, 2114);
        fData->FixParameter(25, 20);
      }
//printf("FRAC_D = %f\n",FRAC_D);
      fData->FixParameter(3,FRAC_D);
      hDataToFit->Fit(fData,"Q, S, N, R, M");   

      RADIUS = fData->GetParameter(0);
      LAM_GEN = fData->GetParameter(2);
      //FRAC_D = fData->GetParameter(3);
      AMP_DELTA = fData->GetParameter(4);
      CKCUTOFF = fData->GetParameter(5);
      DELTA_MASS = fData->GetParameter(11);
      DELTA_WIDTH = fData->GetParameter(12);
      PS_PT = fData->GetParameter(13);
      PS_TEMP = fData->GetParameter(14);
      CHI2_500 = 0;
      for(unsigned uBin=0; uBin<hDataToFit->GetNbinsX(); uBin++){
        double MOM = hDataToFit->GetBinCenter(uBin+1);
        if(MOM>500.001) break;
        if(hDataToFit->GetBinError(uBin+1))
          CHI2_500 += pow((hDataToFit->GetBinContent(uBin+1) - fData->Eval(MOM))/hDataToFit->GetBinError(uBin+1),2.);
        else
          CHI2_500 += pow((hDataToFit->GetBinContent(uBin+1) - fData->Eval(MOM))/1.,2.);
      }
      //printf("chi500 = %e\n", CHI2_500);

      clean_graph(&gData);
      clean_graph(&gFit);
      clean_graph(&gFemto);
      clean_graph(&gBaseline);
      clean_graph(&gDelta);
      clean_graph(&gMadeUpReso);

      for(unsigned uBin=0; uBin<hDataToFit->GetNbinsX(); uBin++){
        double MOM = hDataToFit->GetBinCenter(uBin+1);
        if(MOM>FIT_RANGE) break;
        gData.SetPoint(uBin, MOM, hDataToFit->GetBinContent(uBin+1));
        gData.SetPointError(uBin, hDataToFit->GetBinWidth(uBin+1)*0.5, hDataToFit->GetBinError(uBin+1));
        gFit.SetPoint(uBin, MOM, fData->Eval(MOM));
      }


      fOutputFile.cd();
      for(int iPar=0; iPar<NumFitPar; iPar++){
        fDataDummy->FixParameter(iPar, fData->GetParameter(iPar));
      }
      hDataToFit->SetName(hDataToFit->GetName()+TString::Format("_uIter%i",uIter));
      fData->SetName(fData->GetName()+TString::Format("_uIter%i",uIter));
      //hDataToFit->Write();
      //fData->Write();

      //return Baseline*(par[2]*Femto + par[3]*Delta + (1.-par[2]));
      double par_bl0 = fData->GetParameter(16);
      double par_bl1 = fData->GetParameter(17);
      double par_bl2 = fData->GetParameter(18);
      double par_bl3 = fData->GetParameter(19);
      double par_bl4 = fData->GetParameter(20);
      double par_fmt = fData->GetParameter(2);//lambda_gen
      double par_dlt = fData->GetParameter(4);
      double par_frc = fData->GetParameter(3);
      double par_er_n = fData->GetParameter(21);
      double par_er_Mm = fData->GetParameter(22);
      double par_er_Md1 = fData->GetParameter(23);
      double par_er_Md2 = fData->GetParameter(24);
      double par_er_w = fData->GetParameter(25);

      fDataDummy->FixParameter(2, 0);
      fDataDummy->FixParameter(4, 0);
      for(unsigned uBin=0; uBin<hDataToFit->GetNbinsX(); uBin++){
        double MOM = hDataToFit->GetBinCenter(uBin+1);
        if(MOM>FIT_RANGE) break;
        gBaseline.SetPoint(uBin, MOM, fDataDummy->Eval(MOM));
      }

      fDataDummy->FixParameter(16, 1);
      fDataDummy->FixParameter(17, 0);
      fDataDummy->FixParameter(18, 0);
      fDataDummy->FixParameter(19, 0);
      fDataDummy->FixParameter(20, 0);
      fDataDummy->FixParameter(21, 0);
      fDataDummy->FixParameter(2, par_fmt);
      fDataDummy->FixParameter(4, 0);
      fDataDummy->FixParameter(3, 0);
      for(unsigned uBin=0; uBin<hDataToFit->GetNbinsX(); uBin++){
        double MOM = hDataToFit->GetBinCenter(uBin+1);
        if(MOM>FIT_RANGE) break;
        gFemto.SetPoint(uBin, MOM, fDataDummy->Eval(MOM));
      }

      fDataDummy->FixParameter(2, 0);
      fDataDummy->FixParameter(3, par_frc);
      fDataDummy->FixParameter(4, par_dlt);
      for(unsigned uBin=0; uBin<hDataToFit->GetNbinsX(); uBin++){
        double MOM = hDataToFit->GetBinCenter(uBin+1);
        if(MOM>FIT_RANGE) break;
        gDelta.SetPoint(uBin, MOM, fDataDummy->Eval(MOM));
      }

      //fDataDummy->FixParameter(3, 0);
      fDataDummy->FixParameter(4, 0);
      fDataDummy->FixParameter(21, par_er_n);
      for(unsigned uBin=0; uBin<hDataToFit->GetNbinsX(); uBin++){
        double MOM = hDataToFit->GetBinCenter(uBin+1);
        if(MOM>FIT_RANGE) break;
        gMadeUpReso.SetPoint(uBin, MOM, fDataDummy->Eval(MOM));
        //printf("%f %f\n", MOM, fDataDummy->Eval(MOM));
      }


  //if(TMath::Nint(par[7])==0){
  //  return Baseline*(par[2]*(1.-par[3])*Femto + par[4]*Delta + (1.-par[2]) + MadeUpResonance);
  //}
  //else{
  //  return Baseline*(par[2]*(1.-par[3])*Femto + par[4]*Delta*Femto + (1.-par[2]) + MadeUpResonance);
  //}


      //evaluate the integrals of the SE related to d from Delta and d not from Delta
      //the norm of the integral is not relevant, as we only care about the ratio FRAC_D
      double Integral_Primary = 0;
      double Integral_FromDelta = 0;
      double Integral_MadeUpReso = 0;
      FRAC_D = 0;
      if(hME){
        for(unsigned uBin=0; uBin<hME->GetNbinsX(); uBin++){
          double MOM = hME->GetBinCenter(uBin+1);
          if(MOM<FIT_RANGE){
            Integral_Primary += hME->GetBinContent(uBin+1) * (gFemto.Eval(MOM) - (1.-par_fmt));//original
            //Integral_Primary += hME->GetBinContent(uBin+1) * (gFemto.Eval(MOM) - par_fmt);
            Integral_FromDelta += hME->GetBinContent(uBin+1) * (gDelta.Eval(MOM)-1);
            Integral_MadeUpReso += hME->GetBinContent(uBin+1) * (gMadeUpReso.Eval(MOM)-1);
            //printf("at %.0f += %.3e (%.3e) and %.3e -> %.3e vs %.3e\n",MOM,gFemto.Eval(MOM)- (1.-par_fmt),gFemto.Eval(MOM)- par_fmt,gDelta.Eval(MOM)-1,Integral_Primary,Integral_FromDelta);
          }
          else{
            //we assume gFemto is 1, while gDelta is zero
            Integral_Primary += hME->GetBinContent(uBin+1);
          }
        }   
        //N.B. this is NOT used in the pi-d analysis, something for the Kd only
        //it still plays some role for the normalization I feel like, this is our par[3]
        FRAC_D = Integral_FromDelta/(Integral_FromDelta+Integral_Primary);   
      }
      //printf("%i FRAC_D = %.3e at %p\n",iRepeat,FRAC_D, hME);
      iRepeat++;
    }
    while(iRepeat<RepeatFit);
    
    //printf("OUT OF WHILE\n");
    //printf("%p %p\n", hSE, hME);


    //N.B. ACTUALLY, gDelta has kind of the lambda effectively in. I mean, no explicitely, but it it is multiplied
    //by a separate fit factor, where lambda can be factorized out
    double IntME = hME->Integral();
    double IntSE = hSE?hSE->Integral():0;
    NUM_DELTAS = 0;
    NUM_MADEUP_RESO = 0;
    double NUM_DELTAS_v1 = 0;//original, but scale by the lambda par
    double NUM_DELTAS_v2 = 0;
    double NUM_DELTAS_v3 = 0;//aprox fsi is not there
    for(unsigned uBin=0; uBin<hME->GetNbinsX(); uBin++){
      //checked, this is correct
      double MOM = hME->GetBinCenter(uBin+1);
      NUM_DELTAS += hME->GetBinContent(uBin+1) * (gDelta.Eval(MOM)-1);
      NUM_MADEUP_RESO += hME->GetBinContent(uBin+1) * (gMadeUpReso.Eval(MOM)-1);
      NUM_DELTAS_v1 += hME->GetBinContent(uBin+1) * LAM_GEN*(gDelta.Eval(MOM)-1);
    }
    for(unsigned uBin=0; uBin<hSE->GetNbinsX(); uBin++){
      double MOM = hSE->GetBinCenter(uBin+1);
      double CF = gFemto.Eval(MOM);
      double CD = (gDelta.Eval(MOM)-1)/LAM_GEN;
      //printf("%f %f\n",(1.-FRAC_D),LAM_GEN);
      //printf("k=%.0f, SE=%f, CD=%f, CF=%f, x%f\n",MOM,hSE->GetBinContent(uBin+1),CD,CF,(CF*CD)/(CF*CD+CF*(1.-FRAC_D)+1-LAM_GEN));
      NUM_DELTAS_v2 += hSE->GetBinContent(uBin+1)*(CF*CD)/(CF*CD+CF*(1.-FRAC_D)+1-LAM_GEN);
      NUM_DELTAS_v3 += hSE->GetBinContent(uBin+1)*(LAM_GEN*CD)/(LAM_GEN*CD+LAM_GEN*(1.-FRAC_D)+1-LAM_GEN);
    }

    NUM_DELTAS *= IntSE;
    NUM_DELTAS /= IntME;

    NUM_MADEUP_RESO *= IntSE;
    NUM_MADEUP_RESO /= IntME;    

    //printf("NUM_DELTAS = %f vs %f vs %f vs %f\n",NUM_DELTAS,NUM_DELTAS_v1,NUM_DELTAS_v2,NUM_DELTAS_v3);
    //printf("NUM_MADEUP_RESO = %f\n", NUM_MADEUP_RESO);


    fOutputFile.cd();
    pi_d_Tree->Fill();

    //gData.Write();
    //gFit.Write();
    //gFemto.Write();
    //gDelta.Write();
    //gBaseline.Write();

    //Progress = double(uIter)/double(NumIter+1);

    //printf("end of loop\n");
    delete fData;
    delete fDataDummy;
    delete hDataToFit;
    delete fInputData;
    if(hData && mt_bin==1000) {delete hData; hData=NULL;}
  }//for iter

  fOutputFile.cd();
  pi_d_Tree->Write();

  //fInputData.Close();
  fOutputFile.Close();
}

//the SI is set to -0.037 fm as by Wioleta
void pim_d_Coulomb_vs_RealSI(){
  CATS Kitty;
  Kitty.SetMomBins(100,0,400);
  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());
  AnalysisObject.SetUpCats_pi_d(Kitty, "Gauss", "Gauss", -1, 0);
  Kitty.SetQ1Q2(-1);
  Kitty.KillTheCat();

  CATS KittyC;
  KittyC.SetMomBins(100,0,400);
  AnalysisObject.SetUpCats_pi_d(KittyC, "", "Gauss", 0, 0);
  KittyC.SetQ1Q2(-1);
  KittyC.KillTheCat();

  TGraph g_pim_d_SIplusC;
  g_pim_d_SIplusC.SetName("g_pim_d_SIplusC");

  TGraph g_pim_d_C;
  g_pim_d_C.SetName("g_pim_d_C");

  for(unsigned uMom=0; uMom<Kitty.GetNumMomBins(); uMom++){
    g_pim_d_SIplusC.SetPoint(uMom, Kitty.GetMomentum(uMom), Kitty.GetCorrFun(uMom));
    g_pim_d_C.SetPoint(uMom, KittyC.GetMomentum(uMom), KittyC.GetCorrFun(uMom));
  }

  TFile fOutputFile(TString::Format("%s/pi_d/RealPotential/pim_d.root", GetCernBoxDimi()), "recreate");
  g_pim_d_SIplusC.Write();
  g_pim_d_C.Write();
}


//something bhawani gave you to run
void PiDCF()
{
    TString FileName = "PiDCF";
    TString FileNameOutput = TString::Format("%s/Deuteron/%s.root",GetFemtoOutputFolder(),FileName.Data());
    TFile *OutputFile = new TFile(FileNameOutput, "recreate");
    printf("File Created--- For pp AV18\n");
    const double SourceSize = 1.51;
    const double kMin = 0;
    const double kMax = 300;
    const unsigned NumMomBins = 60;
    const double Md = 1875.61294257;
    const double Mpi = 139.57039;
    const double Mred = Md * Mpi / (Md + Mpi);
    CATS Kitty_piDP;
    CATS Kitty_piDM;
    Kitty_piDP.SetMomBins(NumMomBins, kMin, kMax);
    Kitty_piDP.SetExcludeFailedBins(false);
    Kitty_piDP.SetQ1Q2(1);
    Kitty_piDP.SetPdgId(211, 1000010020);
    Kitty_piDP.SetRedMass(Mred);
    Kitty_piDP.SetNumChannels(1);
    Kitty_piDP.SetNumPW(0, 1);
    Kitty_piDP.SetSpin(0, 1);
    Kitty_piDP.SetChannelWeight(0, 1.0);
    CATSparameters cPars(CATSparameters::tSource, 1, true);
    cPars.SetParameter(0, SourceSize);
    Kitty_piDP.SetAnaSource(GaussSource, cPars);
    Kitty_piDP.SetUseAnalyticSource(true);
    Kitty_piDP.SetAutoNormSource(false);
    Kitty_piDP.KillTheCat();
    DLM_Ck Ck_piDP(Kitty_piDP.GetNumSourcePars(), 0, Kitty_piDP);
    Ck_piDP.Update();
    Kitty_piDM.SetMomBins(NumMomBins, kMin, kMax);
    Kitty_piDM.SetExcludeFailedBins(false);
    Kitty_piDM.SetQ1Q2(-1);
    Kitty_piDM.SetPdgId(211, 1000010020);
    Kitty_piDM.SetRedMass(Mred);
    Kitty_piDM.SetNumChannels(1);
    Kitty_piDM.SetNumPW(0, 1);
    Kitty_piDM.SetSpin(0, 1);
    Kitty_piDM.SetChannelWeight(0, 1.0);
    CATSparameters cPars2(CATSparameters::tSource, 1, true);
    cPars2.SetParameter(0, SourceSize);
    Kitty_piDM.SetAnaSource(GaussSource, cPars2);
    Kitty_piDM.SetUseAnalyticSource(true);
    Kitty_piDM.SetAutoNormSource(false);
    Kitty_piDM.KillTheCat();
    DLM_Ck Ck_piDM(Kitty_piDM.GetNumSourcePars(), 0, Kitty_piDM);
    Ck_piDM.Update();
    printf("File for CF  Created\n");
    TGraph gKitty_piDP;
    TGraph gKitty_piDM;
    gKitty_piDP.SetName(TString::Format("gKittypi-dPositive"));
    gKitty_piDM.SetName(TString::Format("gKittypi-dNeg"));
    for (unsigned uBin = 0; uBin < NumMomBins; uBin++)
    {
        double kst = Kitty_piDP.GetMomentum(uBin);
        double CkP = Ck_piDP.Eval(kst);
        double CkM = Ck_piDM.Eval(kst);
        printf("C(%.2f) = %.2f\n", kst, CkP);
        gKitty_piDP.SetPoint(uBin, kst, CkP);
        gKitty_piDM.SetPoint(uBin, kst, CkM);
    }
    gKitty_piDP.Write();
    gKitty_piDM.Write();
    delete OutputFile;
}


void Error_propagation_final_result(){
  double f_Delta_exp = 0.54;
  double err_f_Delta_exp = 0.06;

  double P_Delta_shm = 0.4305;
  double err_P_Delta_shm = 0.0065;

  double P_Reso_shm = 0.645;
  double err_P_Reso_shm = 0.003;

  unsigned NumIter = 1000*1000;

  double eff = 0.7;
  double err_eff = 0.05;

  double f_Delta;
  double f_Reso;
  TRandom3 rangen(23);
  TH1F* h_f_Delta = new TH1F("h_f_Delta", "h_f_Delta", 4096, 0, 110);
  TH1F* h_f_Reso = new TH1F("h_f_Reso", "h_f_Reso", 4096, 0, 110);
  TH1F* h_f_Delta_shm = new TH1F("h_f_Delta_shm", "h_f_Delta_shm", 4096, 0, 110);
  TH1F* h_f_Reso_shm = new TH1F("h_f_Reso_shm", "h_f_Reso_shm", 4096, 0, 110);
  for(unsigned uIter=0; uIter<NumIter; uIter++){
    double fdp = rangen.Gaus(f_Delta_exp, err_f_Delta_exp);
    double eps = rangen.Gaus(eff, err_eff);
    f_Delta = fdp/eps*(1.-fdp*0.25/eps);
    
    double pds = rangen.Gaus(P_Delta_shm, err_P_Delta_shm);
    double prs = rangen.Gaus(P_Reso_shm, err_P_Reso_shm);
    double f_Delta_shm = 2*pds*(1.-pds*0.5);
    double f_Reso_shm = 2*prs*(1.-prs*0.5);


    f_Reso = f_Delta*f_Reso_shm/f_Delta_shm;
    h_f_Delta->Fill(f_Delta*100);
    h_f_Reso->Fill(f_Reso*100);

    h_f_Delta_shm->Fill(f_Delta_shm*100);
    h_f_Reso_shm->Fill(f_Reso_shm*100);    
  }

  double f_Delta_Min;
  double f_Delta_Max;
  GetCentralInterval(*h_f_Delta, 0.68, f_Delta_Min, f_Delta_Max, true);
  f_Delta = (f_Delta_Max+f_Delta_Min)*0.5;
  double err_f_Delta = (f_Delta_Max-f_Delta_Min)*0.5;

  double f_Reso_Min;
  double f_Reso_Max;
  GetCentralInterval(*h_f_Reso, 0.68, f_Reso_Min, f_Reso_Max, true); 
  f_Reso = (f_Reso_Max+f_Reso_Min)*0.5;
  double err_f_Reso = (f_Reso_Max-f_Reso_Min)*0.5;


  double f_Delta_Min_shm;
  double f_Delta_Max_shm;
  GetCentralInterval(*h_f_Delta_shm, 0.68, f_Delta_Min_shm, f_Delta_Max_shm, true);
  double f_Delta_shm = (f_Delta_Max_shm+f_Delta_Min_shm)*0.5;
  double err_f_Delta_shm = (f_Delta_Max_shm-f_Delta_Min_shm)*0.5;

  double f_Reso_Min_shm;
  double f_Reso_Max_shm;
  GetCentralInterval(*h_f_Reso_shm, 0.68, f_Reso_Min_shm, f_Reso_Max_shm, true); 
  double f_Reso_shm = (f_Reso_Max_shm+f_Reso_Min_shm)*0.5;
  double err_f_Reso_shm = (f_Reso_Max_shm-f_Reso_Min_shm)*0.5;

  TFile fOutput(TString::Format("%s/Deuteron/Error_propagation_final_result.root",GetFemtoOutputFolder()), "recreate");
  h_f_Delta->Write();
  h_f_Reso->Write();
  h_f_Delta_shm->Write();
  h_f_Reso_shm->Write();

  printf("f_Delta = %.1f +/- %.1f %% vs %.1f +/- %.1f %%\n", f_Delta, err_f_Delta, f_Delta_shm, err_f_Delta_shm);
  printf("f_Reso = %.1f +/- %.1f %% vs %.1f +/- %.1f %%\n", f_Reso, err_f_Reso, f_Reso_shm, err_f_Reso_shm);

}

void GypsyGhettoForAcceptance(){
  const double eta_min = -0.8;
  const double eta_max = 0.8;

  const double d_pt_min = 800;
  const double d_pt_max = 2400;

  const double pic_pt_min = 140;
  const double pic_pt_max = 4000;

  //const double Mass_Delta = 1232.-100;
  const double Mass_Delta = 1360;
  //const double Mass_Delta = Mass_p+Mass_pic+1;
  const double mom_mean = Mass_Delta*0.6;

  const unsigned NumIter = 100*1000;

  TRandom3 rangen(23);

  unsigned NumAccDeuterons = 0;//these are essentially all our deuterons from Delta
  unsigned NumAccPions = 0;

  unsigned NumAccBoth = 0;//this is the number of reconstructed Deltas
  unsigned NumAccDeuteronOnly = 0;
  unsigned NumAccPionOnly = 0;

  TH1F* h_d_pt = new TH1F("h_d_pt", "h_d_pt", 1024, 0, 4000);
  TH1F* h_d_eta = new TH1F("h_d_eta", "h_d_eta", 1024, -1.5, 1.5);
  TH1F* h_D_pt = new TH1F("h_D_pt", "h_D_pt", 1024, 0, 4000);

  for(unsigned uIter=0; uIter<NumIter; uIter++){
    double rnd_D_eta;
    double rnd_D_phi;
    double rnd_D_px, rnd_D_py, rnd_D_pz;
    double rnd_D_p;

    rnd_D_eta = rangen.Uniform(eta_min*1.5, eta_max*1.5);
    rnd_D_phi = rangen.Uniform(0, 2.*Pi);

    rnd_D_px = rangen.Gaus(0, mom_mean);
    rnd_D_py = rangen.Gaus(0, mom_mean);
    rnd_D_pz = rangen.Gaus(0, mom_mean);
    rnd_D_p = sqrt(rnd_D_px*rnd_D_px+rnd_D_py*rnd_D_py+rnd_D_pz*rnd_D_pz);

    double D_sin_th = 2.*exp(-rnd_D_eta)/(1.+exp(-2.*rnd_D_eta));
    double D_cotg_th = (1.-exp(-2.*rnd_D_eta))/(2.*exp(-rnd_D_eta));
    double D_cos_th = (1.-exp(-2.*rnd_D_eta))/(1.+exp(-2.*rnd_D_eta));

    double D_px = rnd_D_p*cos(rnd_D_phi)*D_sin_th;
    double D_py = rnd_D_p*sin(rnd_D_phi)*D_sin_th;
    double D_pt = sqrt(D_px*D_px+D_py*D_py);
    double D_pz = D_pt*D_cotg_th;


    double daughter_mass[2];
    daughter_mass[0] = Mass_p;
    daughter_mass[1] = Mass_pic;
    TLorentzVector tlv_D;
    //tlv_D.SetPtEtaPhiM(D_pt,rnd_D_eta,rnd_D_phi,Mass_Delta);
    tlv_D.SetXYZM(rnd_D_px,rnd_D_py,rnd_D_pz,Mass_Delta);

    TGenPhaseSpace GenBod;
    GenBod.SetDecay(tlv_D, 2, daughter_mass);
    GenBod.Generate();
    TLorentzVector& tlv_N = *GenBod.GetDecay(0);
    TLorentzVector& tlv_pic = *GenBod.GetDecay(1);
    TLorentzVector tlv_d;
    tlv_d.SetPtEtaPhiM(tlv_N.Pt()*2,tlv_N.Eta(),tlv_N.Phi(),Mass_d);

    TLorentzVector sum_vector = tlv_d + tlv_pic;;
    TVector3 boost_vector = -sum_vector.BoostVector();

    TLorentzVector tlvb_d(tlv_d);
    tlvb_d.Boost(boost_vector);

    TLorentzVector tlvb_pic(tlv_pic);
    tlvb_pic.Boost(boost_vector);
    
    //TLorentzVector tlvb_sum = tlvb_d+tlvb_pic;
    //tlvb_sum.Print();




    double d_eta = tlv_d.Eta();
    double d_pt = tlv_d.Pt();
    double N_eta = tlv_N.Eta();
    double pic_pt = tlv_pic.Pt();
    double pic_eta = tlv_pic.Eta();

    h_D_pt->Fill(D_pt);
    //printf("%f %f\n",d_pt,d_eta);
    bool PionInAcceptance = pic_eta>eta_min && pic_eta<eta_max && pic_pt>pic_pt_min && pic_pt<pic_pt_max;
    bool DeuteronInAcceptance = d_eta>eta_min && d_eta<eta_max && d_pt>d_pt_min && d_pt<d_pt_max;

    /*
    tlv_N.Print();
    tlv_pic.Print();
    tlv_d.Print();
    printf("PionInAcceptance = %i\n",PionInAcceptance);
    printf("DeuteronInAcceptance = %i\n",DeuteronInAcceptance);
    printf("\n");
    */
    double kstar = tlvb_d.P();
    if(uIter==0){
      printf("kstar = %.0f\n",kstar);
    }
      


    //printf("kstar = %.0f\n",kstar);

    if(DeuteronInAcceptance){
      h_d_pt->Fill(d_pt);
      h_d_eta->Fill(d_eta);
      NumAccDeuterons++;
    }

    if(DeuteronInAcceptance && PionInAcceptance){
      NumAccBoth++;
    }
    if(DeuteronInAcceptance && !PionInAcceptance){
      NumAccDeuteronOnly++;
    }
    if(!DeuteronInAcceptance && PionInAcceptance){
      NumAccPionOnly++;
    }

    //delete tlv_N; tlv_N=NULL;
    //delete tlv_pic; tlv_pic=NULL;
  }


  printf("acceptance effect on the Delta: %.4f\n", double(NumAccBoth)/double(NumAccDeuterons));
  printf("NumAccDeuterons = %u\n",NumAccDeuterons);
  printf("NumAccBoth = %u\n",NumAccBoth);
  printf("NumAccDeuteronOnly = %u\n",NumAccDeuteronOnly);
  printf("NumAccPionOnly = %u\n",NumAccPionOnly);

  TFile fOutput(TString::Format("%s/Deuteron/GypsyGhettoForAcceptance.root",GetFemtoOutputFolder()), "recreate");

  h_d_pt->Write();
  h_d_eta->Write();
  h_D_pt->Write();

  delete h_d_pt;
  delete h_d_eta;
  delete h_D_pt;
}






  
//par[5] = Femto
//par[6] = norm
//par[7] = 0 to castrate the pol3
//par[8] = position of the max of the pol3
//par[9] = p3 parameter of the pol3
//par[10] = -1e6 to switch off the pol4
double MC_Dfrac_closure_fit(double* kstar, double* pars){
  double Femto = pars[5];
    //par[0] = NORM
    //par[1/2] = masses of the daughters
    //par[3] = mass (mother)
    //par[4] = width
  double Delta = Sill_kstar(kstar, pars);
  double Baseline = DLM_Baseline(kstar, &pars[6]);
  return (Femto + Delta) * Baseline;
}



void MC_Dfrac_closure(){
  
  TFile fInput(TString::Format("%s/pi_d/Maxi/Merged_MnAM_CC.root",GetCernBoxDimi()), "read");
  TH1F* hSE = (TH1F*)fInput.Get("hSameEvtDPiFullmT");
  TH1F* hME = (TH1F*)fInput.Get("hMixedEvtDPiFull");

  hSE->GetXaxis()->SetLimits(hSE->GetXaxis()->GetXmin()*1000.,hSE->GetXaxis()->GetXmax()*1000.);
  hME->GetXaxis()->SetLimits(hME->GetXaxis()->GetXmin()*1000.,hME->GetXaxis()->GetXmax()*1000.);

  double SE_int = hSE->Integral();
  hSE->Scale(1./hSE->Integral());
  hME->Scale(1./hME->Integral());

  TFile fOutput(TString::Format("%s/pi_d/MC_closure/MC_closure.root",GetCernBoxDimi()), "recreate");

  TH1F* hCk = (TH1F*)hSE->Clone("hCk");
  hCk->Divide(hME);

  TF1* total_fit = new TF1("total_fit",MC_Dfrac_closure_fit,0,1000,11);


  total_fit->SetParameter(0, 1);
  total_fit->FixParameter(1, Mass_p);
  total_fit->FixParameter(2, Mass_pic);
  total_fit->SetParameter(3, 1232);
  total_fit->SetParameter(4, 120);

  total_fit->FixParameter(5, 1);

  total_fit->SetParameter(6, 1);
  total_fit->FixParameter(7, 0);
  total_fit->SetParameter(8, 200);
  total_fit->SetParLimits(8, 5,1300);
  total_fit->SetParameter(9, 0);
  total_fit->SetParLimits(9, -1e-3, 1e-3);
  total_fit->FixParameter(10, -1e6);

  hCk->Fit(total_fit,"S, N, R, M");

  total_fit->FixParameter(5, 0);
  total_fit->FixParameter(6, 1);
  total_fit->FixParameter(7, 0);
  total_fit->FixParameter(8, 0);
  total_fit->FixParameter(9, 0);
  total_fit->FixParameter(10, 0);


  double Integral_FromDelta = 0;
  double Integral_Primary = 0;

  double FRAC_D = 0;

  for(unsigned uBin=0; uBin<hME->GetNbinsX(); uBin++){
    double MOM = hME->GetBinCenter(uBin+1);
    if(MOM<1000){
      Integral_Primary += hME->GetBinContent(uBin+1)*SE_int * (1);//original
      Integral_FromDelta += hME->GetBinContent(uBin+1)*SE_int * (total_fit->Eval(MOM));
    }
    else{
      Integral_Primary += hME->GetBinContent(uBin+1);
    }
  }
  

  FRAC_D = Integral_FromDelta/(Integral_FromDelta+Integral_Primary);   

  printf("FRAC_D = %f\n", FRAC_D);
  printf("Integral_FromDelta = %f\n", Integral_FromDelta);



  hSE->Write();
  hME->Write();
  hCk->Write();
  total_fit->Write();

}

int DEUTERON_MAIN(int argc, char *argv[]){
  //PiDCF();
  //Error_propagation_final_result();

  //p_pn_cumulant();
  //MyOwn_Kd_v2();

  //MyOwn_Kd_v2_CreateKdp();
  //Check_Source_Means();

  //SetUp_Kdp_Kd();
  //Fit_WithCECA_kstarInt(1024,11);

  //pi_d_source();

  //test_sill_invM_to_kstar();
  //test_sill_ps();
  //TString Description, types, int mt_bin, int NumIter, bool Bootstrap=true, bool DataVar=true, bool FitVar=true, int SEED=0
  //BulgarianIndianGhetto("ghetto_output", atoi(argv[1]),atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atoi(argv[5]),atoi(argv[6]));
  //USE THIS:
 // BulgarianIndianGhetto("pid_withMomResoME_2024-12-06", atoi(argv[1]),atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atoi(argv[5]),atoi(argv[6]));
  //0 - do NOT multiply C_Delta
  //1 - do multiply C_Delta by C_genuine
  //11 - as 1, but now also include the extra 2110 reso
 std::vector<int> fit_types = {11};
 //BulgarianIndianGhetto("ExtraReso2", fit_types, atoi(argv[1]),atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atoi(argv[5]),atoi(argv[6]));
 //GypsyGhettoForAcceptance();
 //BulgarianIndianGhetto("TEST", atoi(argv[1]),atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atoi(argv[5]),atoi(argv[6]));
  //BulgarianIndianGhetto(-1,10,true,true,true,23);
//for the mt_bin:
/*
      //mt diff pip_d used for approvals at PF
      else if(mt_bin>=0 && mt_bin<100 && true){
        //mt int pip_d used for aprovals at CF and PF
  if(mt_bin == -1){
  
        //mt diff pim_d used for approvals at PF
  else if(mt_bin>=100 && mt_bin<200 && true){
        //mt int pim_d used for approvals at CF and PF
  else if(mt_bin == -101){
  
       //this is the ULTRA-GHETTO Kminus-d fit, demanded by Laura on 14th Nov
  else if(mt_bin==1000){
*/


  MC_Dfrac_closure();

  //pim_d_Coulomb_vs_RealSI();

  return 0;
}
