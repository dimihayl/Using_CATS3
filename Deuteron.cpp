
#include "Deuteron.h"

#include "EnvVars.h"
#include "CATS.h"
#include "DLM_Source.h"
#include "DLM_Potentials.h"
#include "CommonAnaFunctions.h"
#include "CATSconstants.h"
#include "DLM_Histo.h"
#include "DLM_RootWrapper.h"
#include "DLM_MathFunctions.h"
#include "DLM_Ck.h"

#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TTree.h"

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

  const double lam_gen = 0.91;
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
  Kd_ER_G.SetAnaSource(0, 1.1);
  Kd_ER_G.SetAnaSource(1, 2.14);
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

    //takes a 3D histo of Kstar Rstar Mt and returns a 2D kdp histo of Mt Kstar
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


  for(int uA=0; uA<1; uA++){for(int uB=0; uB<1; uB++){for(int uE=0; uE<3; uE++){
        //#for the sake of simplicity, make it 1.074 +/- 0.39 and min/avg/max = 1.035/1.074/1.113
        std::vector<float> core_radii = {1.035,1.074,1.113};

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

  delete hSource;
  delete hSourceLower;
  delete hSourceUpper;

}



int DEUTERON_MAIN(int argc, char *argv[]){
  //p_pn_cumulant();
  //MyOwn_Kd_v1();
  //SetUp_Kdp_Kd();
  //Fit_WithCECA_kstarInt(1024,11);

  pi_d_source();

  return 0;
}
