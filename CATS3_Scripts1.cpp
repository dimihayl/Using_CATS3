
#include "CATS3_Scripts1.h"
#include "math.h"
#include <stdio.h>
#include "CATS.h"
#include "DLM_Source.h"
#include "DLM_Potentials.h"
#include "DLM_WfModel.h"
#include "DLM_DrawingTools.h"
#include "DLM_Random.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TPaveText.h"

double UsmaniWithExpCore(double* Parameters)
{
  double& r = Parameters[0];
  const int Spin = round(Parameters[2]);
  double& Core_V = Parameters[3];
  double& Core_mu = Parameters[4];

  //Values for the potential
  const double vbar = 6.2;

  const double vsigma = 0.25;

  //const double wc = 2137;

  double x=r*0.7;
  double vc = Parameters[5]/(1+exp((r-Parameters[6])/Parameters[7]));
  double tpi = (1.0+3.0/x+3.0/(x*x)) * (exp(-x)/x) * pow(1.-exp(-2.*r*r),2.);

  double v = 0.;

  if (Spin == 0) v = vc - (vbar + 0.75*vsigma)*tpi*tpi;//Usmani singlet
  else if (Spin == 1)  v = vc - (vbar - 0.25*vsigma)*tpi*tpi;//Usmani triplet
  else printf ("wrong polarization\n");

  v += Core_V*exp(-pow(r/Core_mu,2));

  return v;
}

void UsmaniPlusExpCore(){

    const double kMin_pL=0;
    const double kMax_pL = 300;
    const unsigned NumMomBins_pL = 150;
    const double kStep_pL = (kMax_pL-kMin_pL)/double(NumMomBins_pL);

    //starting value, do not worry about it too much
    const double GaussSourceSize = 0.85;
    const double Mass_p=938.272;
    const double Mass_L=1115.683;

    CATSparameters cPars1S0(CATSparameters::tPotential,6,true);
    cPars1S0.SetParameter(0,0);
    cPars1S0.SetParameter(1,0);
    cPars1S0.SetParameter(2,0.1);
    cPars1S0.SetParameter(3,2137);
    cPars1S0.SetParameter(4,0.5);
    cPars1S0.SetParameter(5,0.2);
    CATSparameters cPars3S1(CATSparameters::tPotential,6,true);
    cPars3S1.SetParameter(0,1);
    cPars3S1.SetParameter(1,0);
    cPars3S1.SetParameter(2,0.5);
    cPars3S1.SetParameter(3,2137);
    cPars3S1.SetParameter(4,0.5);
    cPars3S1.SetParameter(5,0.2);

    CATS AB_pL;
    double Pars_pL[4] = {0,0,0,GaussSourceSize};
    CATSparameters cPars_pL(CATSparameters::tSource,1,true);
    cPars_pL.SetParameters(&Pars_pL[3]);

    AB_pL.SetAnaSource(GaussSource, cPars_pL);
    AB_pL.SetUseAnalyticSource(true);

    AB_pL.SetExcludeFailedBins(false);
    AB_pL.SetMomBins(NumMomBins_pL,kMin_pL,kMax_pL);

    AB_pL.SetNumChannels(1);
    AB_pL.SetNumPW(0,1);
    //AB_pL.SetNumPW(1,1);
    AB_pL.SetSpin(0,0);
    //AB_pL.SetSpin(1,1);
    AB_pL.SetChannelWeight(0, 1.0);
    //AB_pL.SetChannelWeight(1, 0.0);

    AB_pL.SetQ1Q2(0);
    AB_pL.SetPdgId(2212, 3122);
    AB_pL.SetRedMass( (Mass_p*Mass_L)/(Mass_p+Mass_L) );

    AB_pL.SetShortRangePotential(0,0,UsmaniWithExpCore,cPars1S0);
    //AB_pL.SetShortRangePotential(1,0,UsmaniWithExpCore,cPars3S1);

    AB_pL.SetEpsilonConv(1e-9);
    AB_pL.SetEpsilonProp(1e-9);

    AB_pL.KillTheCat();

    double kFitFrom = 4;
    const unsigned FitFromBin = AB_pL.GetMomBin(kFitFrom);
    kFitFrom = AB_pL.GetMomentum(FitFromBin)-0.5*kStep_pL;
    double kFitTo = 48;
    const unsigned FitToBin = AB_pL.GetMomBin(kFitTo);
    kFitTo = AB_pL.GetMomentum(FitToBin)+0.5*kStep_pL;
    const unsigned NumFitBins = FitToBin-FitFromBin+1;

    TFile* fOut = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/UsmaniPlusExpCore/fOut.root","recreate");

    TH1F* hPhaseShift0 = new TH1F("hPhaseShift0", "hPhaseShift0", NumFitBins, kFitFrom, kFitTo);
    for(unsigned uBin=0; uBin<NumFitBins; uBin++){
        //hPhaseShift0->SetBinContent(uBin+1, tan(AB_pL.GetPhaseShift(uBin+FitFromBin,0,0)/AB_pL.GetMomentum(uBin+FitFromBin)));
        hPhaseShift0->SetBinContent(uBin+1, AB_pL.GetMomentum(uBin+FitFromBin)/tan(AB_pL.GetPhaseShift(uBin+FitFromBin,0,0)));
        hPhaseShift0->SetBinError(uBin+1, hPhaseShift0->GetBinContent(uBin+1)*0.0001);
    }

    //this is the eq. for eff. range, but inverted. We use the femto sign convention!
    //TF1* fPhaseShift0 = new TF1("fPhaseShift0", "2.*[0]/197.327/(2.+x*x*[0]*[1]/197.327/197.327)", kFitFrom, kFitTo);
    TF1* fPhaseShift0 = new TF1("fPhaseShift0", "(197.327/[0]+0.5*x*x*[1]/197.327)", kFitFrom, kFitTo);
    fPhaseShift0->SetParameter(0, 0);
    fPhaseShift0->SetParLimits(0, -5,5);
    fPhaseShift0->SetParameter(1, 0);
    fPhaseShift0->SetParLimits(1, -10,10);

    hPhaseShift0->Fit(fPhaseShift0,"S, N, R, M");

    TH1F* hCorrFun0 = new TH1F("hCorrFun0","hCorrFun0",NumMomBins_pL,kMin_pL,kMax_pL);
    for(unsigned uBin=0; uBin<NumMomBins_pL; uBin++){
        hCorrFun0->SetBinContent(uBin+1,AB_pL.GetCorrFun(uBin));
    }

    hPhaseShift0->SetLineColor(kBlue+1);
    fPhaseShift0->SetLineColor(kBlue+1);
    hCorrFun0->SetLineColor(kBlue+1);

    hPhaseShift0->Write();
    fPhaseShift0->Write();
    hCorrFun0->Write();

    cPars1S0.SetParameter(1,10000);
    cPars1S0.SetParameter(2,0.1);
    //cPars1S0.SetParameter(3,1050);
    //cPars1S0.SetParameter(4,0.5);
    //cPars1S0.SetParameter(5,0.2);
    AB_pL.KillTheCat(CATS::kPotentialChanged);

    TH1F* hPhaseShift1 = new TH1F("hPhaseShift1", "hPhaseShift1", NumFitBins, kFitFrom, kFitTo);
    for(unsigned uBin=0; uBin<NumFitBins; uBin++){
        //hPhaseShift1->SetBinContent(uBin+1, tan(AB_pL.GetPhaseShift(uBin+FitFromBin,0,0)/AB_pL.GetMomentum(uBin+FitFromBin)));
        hPhaseShift1->SetBinContent(uBin+1, AB_pL.GetMomentum(uBin+FitFromBin)/tan(AB_pL.GetPhaseShift(uBin+FitFromBin,0,0)));
        hPhaseShift1->SetBinError(uBin+1, hPhaseShift1->GetBinContent(uBin+1)*0.0001);
    }

    //this is the eq. for eff. range, but inverted. We use the femto sign convention!
    //TF1* fPhaseShift0 = new TF1("fPhaseShift0", "2.*[0]/197.327/(2.+x*x*[0]*[1]/197.327/197.327)", kFitFrom, kFitTo);
    TF1* fPhaseShift1 = new TF1("fPhaseShift1", "(197.327/[0]+0.5*x*x*[1]/197.327)", kFitFrom, kFitTo);
    fPhaseShift1->SetParameter(0, 0);
    fPhaseShift1->SetParLimits(0, -5,5);
    fPhaseShift1->SetParameter(1, 0);
    fPhaseShift1->SetParLimits(1, -10,10);

    hPhaseShift1->Fit(fPhaseShift1,"S, N, R, M");

    TH1F* hCorrFun1 = new TH1F("hCorrFun1","hCorrFun1",NumMomBins_pL,kMin_pL,kMax_pL);
    for(unsigned uBin=0; uBin<NumMomBins_pL; uBin++){
        hCorrFun1->SetBinContent(uBin+1,AB_pL.GetCorrFun(uBin));
    }

    hPhaseShift1->SetLineColor(kRed);
    fPhaseShift1->SetLineColor(kRed);
    hCorrFun1->SetLineColor(kRed);

    hPhaseShift1->Write();
    fPhaseShift1->Write();
    hCorrFun1->Write();

    delete hPhaseShift0;
    delete fPhaseShift0;
    delete hCorrFun0;
    delete hPhaseShift1;
    delete fPhaseShift1;
    delete hCorrFun1;
    delete fOut;
}

void Test_Tom_pXi(){
    const double kMax=320;
    complex<double>**** WaveFunctionU; double* RadBins; unsigned NumRadBins=0;

    CATSparameters pSource(CATSparameters::tSource,1,true);
    pSource.SetParameter(0,3.4);

    CATS Kitty;
    Kitty.SetUseAnalyticSource(true);
    Kitty.SetAnaSource(GaussSource, pSource);
    Kitty.SetExcludeFailedBins(false);

    InitESC16_v2("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/Tom/Tom_10122018wfs/",Kitty,&WaveFunctionU,&RadBins,NumRadBins,kMax);

    Kitty.SetChannelWeight(0,1./4.);
    Kitty.SetChannelWeight(1,3./4.);
    Kitty.SetChannelWeight(2,0);
    Kitty.SetChannelWeight(3,0);
    //Kitty.SetGamow(true);

    Kitty.KillTheCat();

    TFile* OutFile = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/Test_Tom_pXi/OutFile.root","recreate");

    TGraph grCk;
    grCk.Set(Kitty.GetNumMomBins());
    grCk.SetName("grCk");
    for(unsigned uBin=0; uBin<Kitty.GetNumMomBins(); uBin++){
        grCk.SetPoint(uBin,Kitty.GetMomentum(uBin),Kitty.GetCorrFun(uBin));
    }
    grCk.Write();

    Kitty.SetChannelWeight(0,1./4.);
    Kitty.SetChannelWeight(1,0);
    Kitty.SetChannelWeight(2,0);
    Kitty.SetChannelWeight(3,0);
    Kitty.KillTheCat();
    TGraph grCk_S0;
    grCk_S0.Set(Kitty.GetNumMomBins());
    grCk_S0.SetName("grCk_S0");
    for(unsigned uBin=0; uBin<Kitty.GetNumMomBins(); uBin++){
        grCk_S0.SetPoint(uBin,Kitty.GetMomentum(uBin),Kitty.GetCorrFun(uBin));
    }
    grCk_S0.Write();

    Kitty.SetChannelWeight(0,0);
    Kitty.SetChannelWeight(1,3./4.);
    Kitty.SetChannelWeight(2,0);
    Kitty.SetChannelWeight(3,0);
    Kitty.KillTheCat();
    TGraph grCk_S1;
    grCk_S1.Set(Kitty.GetNumMomBins());
    grCk_S1.SetName("grCk_S1");
    for(unsigned uBin=0; uBin<Kitty.GetNumMomBins(); uBin++){
        grCk_S1.SetPoint(uBin,Kitty.GetMomentum(uBin),Kitty.GetCorrFun(uBin));
    }
    grCk_S1.Write();

    delete OutFile;
    //CleanHaidenbauer(Kitty,&WaveFunctionU,NULL,&RadBins);
}


void Test_Gamow_pXi(){
    const double kMin = 0;
    const double kMax = 320;
    const unsigned NumMomBins = 80;
    CATS Kitty;
    Kitty.SetUseAnalyticSource(true);
    CATSparameters pSource(CATSparameters::tSource,1,true);
    pSource.SetParameter(0, 1.2);
    Kitty.SetAnaSource(GaussSource, pSource);

    Kitty.SetExcludeFailedBins(false);
    Kitty.SetMomBins(NumMomBins,kMin,kMax);

    Kitty.SetNumChannels(4);
    Kitty.SetNumPW(0,1);
    Kitty.SetNumPW(1,1);
    Kitty.SetNumPW(2,1);
    Kitty.SetNumPW(3,1);
    Kitty.SetSpin(0,0);//I=0; S=0
    Kitty.SetSpin(1,1);//I=0; S=1
    Kitty.SetSpin(2,0);//I=1; S=0
    Kitty.SetSpin(3,1);//I=1; S=1
    Kitty.SetChannelWeight(0, 1./8.);
    Kitty.SetChannelWeight(1, 3./8.);
    Kitty.SetChannelWeight(2, 1./8.);
    Kitty.SetChannelWeight(3, 3./8.);
    Kitty.SetQ1Q2(-1);
    Kitty.SetPdgId(2212, 3312);
    const double Mass_p = 938.272;
    const double Mass_Xim = 1321.7;
    Kitty.SetRedMass( (Mass_p*Mass_Xim)/(Mass_p+Mass_Xim) );
    Kitty.SetMaxRad(64);
    Kitty.SetMaxRho(32);
    Kitty.SetEpsilonConv(1e-8);
    Kitty.SetEpsilonProp(1e-8);

    double pXimPotParsI0S0[11]={0,0,pXim_HALQCD1,12,0,-1,1,0,0,0,0};
    double pXimPotParsI0S1[11]={0,0,pXim_HALQCD1,12,0,-1,1,1,0,1,0};
    double pXimPotParsI1S0[11]={0,0,pXim_HALQCD1,12,1,1,1,0,0,0,0};
    double pXimPotParsI1S1[11]={0,0,pXim_HALQCD1,12,1,1,1,1,0,1,0};

    CATSparameters c_pXimPotParsI0S0(CATSparameters::tPotential,9,true);
    c_pXimPotParsI0S0.SetParameters(&pXimPotParsI0S0[2]);
    CATSparameters c_pXimPotParsI0S1(CATSparameters::tPotential,9,true);
    c_pXimPotParsI0S1.SetParameters(&pXimPotParsI0S1[2]);
    CATSparameters c_pXimPotParsI1S0(CATSparameters::tPotential,9,true);
    c_pXimPotParsI1S0.SetParameters(&pXimPotParsI1S0[2]);
    CATSparameters c_pXimPotParsI1S1(CATSparameters::tPotential,9,true);
    c_pXimPotParsI1S1.SetParameters(&pXimPotParsI1S1[2]);

    Kitty.SetShortRangePotential(0,0,fDlmPot,c_pXimPotParsI0S0);
    Kitty.SetShortRangePotential(1,0,fDlmPot,c_pXimPotParsI0S1);
    Kitty.SetShortRangePotential(2,0,fDlmPot,c_pXimPotParsI1S0);
    Kitty.SetShortRangePotential(3,0,fDlmPot,c_pXimPotParsI1S1);

    Kitty.KillTheCat();

    TFile* OutFile = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/Test_Gamow_pXi/Test_Gamow_pXi.root","recreate");

    TGraph grCk;
    grCk.Set(Kitty.GetNumMomBins());
    grCk.SetName("grCk");
    for(unsigned uBin=0; uBin<Kitty.GetNumMomBins(); uBin++){
        grCk.SetPoint(uBin,Kitty.GetMomentum(uBin),Kitty.GetCorrFun(uBin));
    }
    grCk.Write();

    Kitty.SetGamow(true);
    Kitty.KillTheCat();

    TGraph grCk_Gamow;
    grCk_Gamow.Set(Kitty.GetNumMomBins());
    grCk_Gamow.SetName("grCk_Gamow");
    for(unsigned uBin=0; uBin<Kitty.GetNumMomBins(); uBin++){
        grCk_Gamow.SetPoint(uBin,Kitty.GetMomentum(uBin),Kitty.GetCorrFun(uBin));
    }
    grCk_Gamow.Write();

    delete OutFile;
}

void pXi_PhaseShifts(const int WhichPotential){
    //+11,12,13 => the new ones
    //-1 the old potential, with configuration 12,12,6,6

    const TString OutputFolder = TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pXi_PhaseShifts/090119_%i/",WhichPotential);

    const double kMin = 0;
    const double kMax = 240;
    const unsigned NumMomBins = 120;

    CATS Kitty;
    Kitty.SetUseAnalyticSource(true);
    CATSparameters pSource(CATSparameters::tSource,1,true);
    pSource.SetParameter(0, 1.4);
    Kitty.SetAnaSource(GaussSource, pSource);

    Kitty.SetExcludeFailedBins(false);
    Kitty.SetMomBins(NumMomBins,kMin,kMax);
    const unsigned short NumChannels = 4;
    enum CHANNELS {I0S0,I0S1,I1S0,I1S1};
    const TString ChannelName[NumChannels] = {"I0S0","I0S1","I1S0","I1S1"};
    Kitty.SetNumChannels(NumChannels);
    Kitty.SetNumPW(0,1);
    Kitty.SetNumPW(1,1);
    Kitty.SetNumPW(2,1);
    Kitty.SetNumPW(3,1);
    Kitty.SetSpin(0,0);//I=0; S=0
    Kitty.SetSpin(1,1);//I=0; S=1
    Kitty.SetSpin(2,0);//I=1; S=0
    Kitty.SetSpin(3,1);//I=1; S=1
    Kitty.SetChannelWeight(0, 1./8.);
    Kitty.SetChannelWeight(1, 3./8.);
    Kitty.SetChannelWeight(2, 1./8.);
    Kitty.SetChannelWeight(3, 3./8.);
    Kitty.SetQ1Q2(0);
    Kitty.SetPdgId(2212, 3312);
    const double Mass_p = 938.272;
    const double Mass_Xim = 1321.7;
    Kitty.SetRedMass( (Mass_p*Mass_Xim)/(Mass_p+Mass_Xim) );
    Kitty.SetMaxRad(64);
    Kitty.SetMaxRho(32);
    Kitty.SetEpsilonConv(5e-9);
    Kitty.SetEpsilonProp(5e-9);

    double pXimPotParsI0S0[10]={0,0,pXim_Lattice,12,0,-1,1,0,0,0};
    double pXimPotParsI0S1[10]={0,0,pXim_Lattice,12,0,-1,1,1,0,1};
    double pXimPotParsI1S0[10]={0,0,pXim_Lattice,6,1,1,1,0,0,0};
    double pXimPotParsI1S1[10]={0,0,pXim_Lattice,6,1,1,1,1,0,1};
    CATSparameters c_pXimPotParsI0S0(CATSparameters::tPotential,9,true);
    CATSparameters c_pXimPotParsI0S1(CATSparameters::tPotential,9,true);
    CATSparameters c_pXimPotParsI1S0(CATSparameters::tPotential,9,true);
    CATSparameters c_pXimPotParsI1S1(CATSparameters::tPotential,9,true);

    c_pXimPotParsI0S0.SetParameter(2,0);
    c_pXimPotParsI0S0.SetParameter(3,-1);
    c_pXimPotParsI0S0.SetParameter(4,1);
    c_pXimPotParsI0S0.SetParameter(5,0);
    c_pXimPotParsI0S0.SetParameter(6,0);
    c_pXimPotParsI0S0.SetParameter(7,0);

    c_pXimPotParsI0S1.SetParameter(2,0);
    c_pXimPotParsI0S1.SetParameter(3,-1);
    c_pXimPotParsI0S1.SetParameter(4,1);
    c_pXimPotParsI0S1.SetParameter(5,1);
    c_pXimPotParsI0S1.SetParameter(6,0);
    c_pXimPotParsI0S1.SetParameter(7,1);

    c_pXimPotParsI1S0.SetParameter(2,1);
    c_pXimPotParsI1S0.SetParameter(3,1);
    c_pXimPotParsI1S0.SetParameter(4,1);
    c_pXimPotParsI1S0.SetParameter(5,0);
    c_pXimPotParsI1S0.SetParameter(6,0);
    c_pXimPotParsI1S0.SetParameter(7,0);

    c_pXimPotParsI1S1.SetParameter(2,1);
    c_pXimPotParsI1S1.SetParameter(3,1);
    c_pXimPotParsI1S1.SetParameter(4,1);
    c_pXimPotParsI1S1.SetParameter(5,1);
    c_pXimPotParsI1S1.SetParameter(6,0);
    c_pXimPotParsI1S1.SetParameter(7,1);

    if(WhichPotential==-1){
        c_pXimPotParsI0S0.SetParameter(0,pXim_Lattice);
        c_pXimPotParsI0S0.SetParameter(1,12);

        c_pXimPotParsI0S1.SetParameter(0,pXim_Lattice);
        c_pXimPotParsI0S1.SetParameter(1,12);

        c_pXimPotParsI1S0.SetParameter(0,pXim_Lattice);
        c_pXimPotParsI1S0.SetParameter(1,6);

        c_pXimPotParsI1S1.SetParameter(0,pXim_Lattice);
        c_pXimPotParsI1S1.SetParameter(1,6);
    }
    else if(WhichPotential>0){
        c_pXimPotParsI0S0.SetParameter(0,pXim_HALQCD1);
        c_pXimPotParsI0S0.SetParameter(1,WhichPotential);

        c_pXimPotParsI0S1.SetParameter(0,pXim_HALQCD1);
        c_pXimPotParsI0S1.SetParameter(1,WhichPotential);

        c_pXimPotParsI1S0.SetParameter(0,pXim_HALQCD1);
        c_pXimPotParsI1S0.SetParameter(1,WhichPotential);

        c_pXimPotParsI1S1.SetParameter(0,pXim_HALQCD1);
        c_pXimPotParsI1S1.SetParameter(1,WhichPotential);
    }
    else{
        printf("Един голем!\n");
        return;
    }
    //double pXimPotParsI0S0[11]={0,0,pXim_HALQCD1,FlagQCD,0,-1,1,0,0,0,0};
    //double pXimPotParsI0S1[11]={0,0,pXim_HALQCD1,FlagQCD,0,-1,1,1,0,1,0};
    //double pXimPotParsI1S0[11]={0,0,pXim_HALQCD1,FlagQCD,1,1,1,0,0,0,0};
    //double pXimPotParsI1S1[11]={0,0,pXim_HALQCD1,FlagQCD,1,1,1,1,0,1,0};

    Kitty.SetShortRangePotential(0,0,fDlmPot,c_pXimPotParsI0S0);
    Kitty.SetShortRangePotential(1,0,fDlmPot,c_pXimPotParsI0S1);
    Kitty.SetShortRangePotential(2,0,fDlmPot,c_pXimPotParsI1S0);
    Kitty.SetShortRangePotential(3,0,fDlmPot,c_pXimPotParsI1S1);

    Kitty.SetMaxNumThreads(8);
    Kitty.KillTheCat();

    Kitty.SetNotifications(CATS::nSilent);
    /*
    for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
        Kitty.SetChannelWeight(0, 1./8.);
        Kitty.SetChannelWeight(1, 3./8.);
        Kitty.SetChannelWeight(2, 1./8.);
        Kitty.SetChannelWeight(3, 3./8.);
        Kitty.KillTheCat();
        printf("C(%.0f) = %.3f\n",Kitty.GetMomentum(uMomBin),Kitty.GetCorrFun(uMomBin));
        Kitty.SetChannelWeight(0, 1.);
        Kitty.SetChannelWeight(1, 0.);
        Kitty.SetChannelWeight(2, 0.);
        Kitty.SetChannelWeight(3, 0.);
        Kitty.KillTheCat();
        printf("          %.3f\n",Kitty.GetCorrFun(uMomBin));
        Kitty.SetChannelWeight(0, 0.);
        Kitty.SetChannelWeight(1, 1.);
        Kitty.SetChannelWeight(2, 0.);
        Kitty.SetChannelWeight(3, 0.);
        Kitty.KillTheCat();
        printf("          %.3f\n",Kitty.GetCorrFun(uMomBin));
        Kitty.SetChannelWeight(0, 0.);
        Kitty.SetChannelWeight(1, 0.);
        Kitty.SetChannelWeight(2, 1.);
        Kitty.SetChannelWeight(3, 0.);
        Kitty.KillTheCat();
        printf("          %.3f\n",Kitty.GetCorrFun(uMomBin));
        Kitty.SetChannelWeight(0, 0.);
        Kitty.SetChannelWeight(1, 0.);
        Kitty.SetChannelWeight(2, 0.);
        Kitty.SetChannelWeight(3, 1.);
        Kitty.KillTheCat();
        printf("          %.3f\n",Kitty.GetCorrFun(uMomBin));
    }
    */


    DLM_DtColor DlmCol;

    TF1* ZeroLine = new TF1("ZeroLine", "0", kMin, kMax);
    ZeroLine->SetLineColor(kBlack);
    ZeroLine->SetLineStyle(2);
    ZeroLine->SetLineWidth(4);

    const double kFitFrom = 3;
    const double kFitTo = 39;
    const unsigned NumFitBins = Kitty.GetMomBin(kFitTo)-Kitty.GetMomBin(kFitFrom)+1;
    unsigned FitFromBin = Kitty.GetMomBin(kFitFrom);
    unsigned FitToBin = Kitty.GetMomBin(kFitTo);
    TH1F** hPhaseShift = new TH1F* [NumChannels];
    TF1** fPhaseShift = new TF1* [NumChannels];
    TGraph** gDelta = new TGraph* [NumChannels];
    TCanvas* cDelta = new TCanvas("cDelta", "cDelta", 1);
    cDelta->cd(0);
    cDelta->SetCanvasSize(1280, 720);
    cDelta->SetMargin(0.15,0.05,0.2,0.05);//lrbt
    cDelta->SetGrid();
    cDelta->cd(0);
    TH1F* hAxis = new TH1F("hAxis","hAxis", NumMomBins, kMin, kMax);
    hAxis->SetStats(false);
    hAxis->SetTitle("");
    hAxis->GetXaxis()->SetLabelSize(0.06);
    hAxis->GetXaxis()->SetTitle("k (MeV)");
    hAxis->GetXaxis()->CenterTitle();
    hAxis->GetXaxis()->SetTitleOffset(1.15);
    hAxis->GetXaxis()->SetLabelOffset(0.02);
    hAxis->GetXaxis()->SetTitleSize(0.07);
    hAxis->GetYaxis()->SetLabelSize(0.06);
    hAxis->GetYaxis()->SetTitle("#delta(k) (deg)");
    hAxis->GetYaxis()->CenterTitle();
    hAxis->GetYaxis()->SetTitleOffset(0.75);
    hAxis->GetYaxis()->SetTitleSize(0.07);
    hAxis->GetYaxis()->SetRangeUser(-28, 55);
    hAxis->GetYaxis()->SetLimits(-28, 55);
    hAxis->Draw("");

    TCanvas* cDeltaZoom = new TCanvas("cDeltaZoom", "cDeltaZoom", 1);
    cDeltaZoom->cd(0);
    cDeltaZoom->SetCanvasSize(1280, 720);
    cDeltaZoom->SetMargin(0.15,0.05,0.2,0.05);//lrbt
    cDeltaZoom->SetGrid();
    cDeltaZoom->cd(0);
    TH1F* hAxisZoom = new TH1F("hAxisZoom","hAxisZoom", NumMomBins, kMin, kMax);
    hAxisZoom->SetStats(false);
    hAxisZoom->SetTitle("");
    hAxisZoom->GetXaxis()->SetLabelSize(0.06);
    hAxisZoom->GetXaxis()->SetTitle("k (MeV)");
    hAxisZoom->GetXaxis()->CenterTitle();
    hAxisZoom->GetXaxis()->SetTitleOffset(1.15);
    hAxisZoom->GetXaxis()->SetLabelOffset(0.02);
    hAxisZoom->GetXaxis()->SetTitleSize(0.07);
    hAxisZoom->GetYaxis()->SetLabelSize(0.06);
    hAxisZoom->GetYaxis()->SetTitle("#delta(k) (deg)");
    hAxisZoom->GetYaxis()->CenterTitle();
    hAxisZoom->GetYaxis()->SetTitleOffset(0.75);
    hAxisZoom->GetYaxis()->SetTitleSize(0.07);
    hAxisZoom->GetYaxis()->SetRangeUser(-10, 40);
    hAxisZoom->GetYaxis()->SetLimits(-10, 40);
    hAxisZoom->Draw("");

    for(unsigned short usCh=0; usCh<NumChannels; usCh++){
        hPhaseShift[usCh] = new TH1F("hPhaseShift"+ChannelName[usCh], "hPhaseShift"+ChannelName[usCh], NumFitBins, Kitty.GetMomBinLowEdge(FitFromBin), Kitty.GetMomBinUpEdge(FitToBin));
        for(unsigned uFitBin=0; uFitBin<NumFitBins; uFitBin++){
            hPhaseShift[usCh]->SetBinContent(uFitBin+1, Kitty.GetMomentum(uFitBin+FitFromBin)/tan(Kitty.GetPhaseShift(uFitBin+FitFromBin,usCh,0)));
            hPhaseShift[usCh]->SetBinError(uFitBin+1, hPhaseShift[usCh]->GetBinContent(uFitBin+1)*0.0001);
        }

        fPhaseShift[usCh] = new TF1("fPhaseShift"+ChannelName[usCh], "(197.327/[0]+0.5*x*x*[1]/197.327)", kFitFrom, kFitTo);
        fPhaseShift[usCh]->SetParameter(0, 0);
        fPhaseShift[usCh]->SetParLimits(0, -10,10);
        fPhaseShift[usCh]->SetParameter(1, 2);
        fPhaseShift[usCh]->SetParLimits(1, 0,256);
        hPhaseShift[usCh]->Fit(fPhaseShift[usCh],"S, N, R, M");

        //gDelta[usCh] = new TH1F("hDelta"+ChannelName[usCh], "hDelta"+ChannelName[usCh], NumMomBins, kMin, kMax);
        gDelta[usCh] = new TGraph();
        gDelta[usCh]->SetName("gDelta"+ChannelName[usCh]);
        gDelta[usCh]->Set(NumMomBins);
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            //gDelta[usCh]->SetBinContent(uMomBin+1,Kitty.GetPhaseShift(uMomBin,usCh,0)*180./TMath::Pi());
            gDelta[usCh]->SetPoint(uMomBin,Kitty.GetMomentum(uMomBin),Kitty.GetPhaseShift(uMomBin,usCh,0)*180./TMath::Pi());
        }
        gDelta[usCh]->SetLineWidth(6);
        gDelta[usCh]->SetLineColor(DlmCol.GetColor(usCh+1));

        cDelta->cd(0);
        gDelta[usCh]->Draw("same,C");
        cDeltaZoom->cd(0);
        gDelta[usCh]->Draw("same,C");
    }

    TLegend* myLegend = new TLegend(0.15,0.9,0.95,0.99);//lbrt
    myLegend->SetName("myLegend");
    myLegend->SetTextSize(0.055);
    myLegend->SetNColumns(NumChannels);
    for(unsigned short usCh=0; usCh<NumChannels; usCh++){
        myLegend->AddEntry(gDelta[usCh],ChannelName[usCh]);
    }
    cDelta->cd(0);
    myLegend->Draw("same");
    ZeroLine->Draw("same");
    cDeltaZoom->cd(0);
    myLegend->Draw("same");
    ZeroLine->Draw("same");

    TFile* fOut = new TFile(OutputFolder+"fOut.root","recreate");
    cDelta->Write();
    cDelta->SaveAs(OutputFolder+"cDelta.png");
    cDeltaZoom->Write();
    cDeltaZoom->SaveAs(OutputFolder+"cDeltaZoom.png");

    for(unsigned short usCh=0; usCh<NumChannels; usCh++){
        hPhaseShift[usCh]->Write();
        fPhaseShift[usCh]->Write();
        gDelta[usCh]->Write();
    }

//////////////////////////////////////////////////////////////////////////////

    TGraph** gCorrFun = new TGraph* [NumChannels];
    for(unsigned short usCh=0; usCh<NumChannels; usCh++){
        Kitty.SetChannelWeight(0, 0.);
        Kitty.SetChannelWeight(1, 0.);
        Kitty.SetChannelWeight(2, 0.);
        Kitty.SetChannelWeight(3, 0.);
        Kitty.SetChannelWeight(usCh, 1.);
        Kitty.KillTheCat();

        gCorrFun[usCh] = new TGraph();
        gCorrFun[usCh]->SetName("gCorrFun"+ChannelName[usCh]);
        gCorrFun[usCh]->Set(NumMomBins);
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            gCorrFun[usCh]->SetPoint(uMomBin,Kitty.GetMomentum(uMomBin),Kitty.GetCorrFun(uMomBin));
        }
        gCorrFun[usCh]->SetLineWidth(6);
        gCorrFun[usCh]->SetLineColor(DlmCol.GetColor(usCh+1));
    }

    TF1* UnitLine = new TF1("UnitLine", "1", kMin, kMax);
    UnitLine->SetLineColor(kBlack);
    UnitLine->SetLineStyle(2);
    UnitLine->SetLineWidth(4);

    TPaveText* Info1 = new TPaveText(0.65,0.65,0.95,0.9, "blNDC");//lbrt
    Info1->SetName("Info1");
    Info1->SetBorderSize(1);
    Info1->SetTextSize(0.05);
    Info1->SetFillColor(kWhite);
    Info1->SetTextFont(22);
    for(unsigned short usCh=0; usCh<NumChannels; usCh++){
        Info1->AddText(TString::Format("f_{0}(%s) = %.2f fm",ChannelName[usCh].Data(),fPhaseShift[usCh]->GetParameter(0)));
    }

    TCanvas* c_Ck = new TCanvas("c_Ck", "c_Ck", 1);
    c_Ck->cd(0);
    c_Ck->SetCanvasSize(1280, 720);
    c_Ck->SetMargin(0.15,0.05,0.2,0.05);//lrbt
    c_Ck->SetGrid();
    c_Ck->cd(0);

    TH1F* hAxis_Ck = new TH1F("hAxis_Ck","hAxis_Ck",NumMomBins,kMin,kMax);
    hAxis_Ck->SetStats(false);
    hAxis_Ck->SetTitle("");
    hAxis_Ck->GetXaxis()->SetLabelSize(0.06);
    hAxis_Ck->GetXaxis()->SetTitle("k (MeV)");
    hAxis_Ck->GetXaxis()->CenterTitle();
    hAxis_Ck->GetXaxis()->SetTitleOffset(1.15);
    hAxis_Ck->GetXaxis()->SetLabelOffset(0.02);
    hAxis_Ck->GetXaxis()->SetTitleSize(0.07);
    hAxis_Ck->GetYaxis()->SetLabelSize(0.06);
    hAxis_Ck->GetYaxis()->SetTitle("C(k)");
    hAxis_Ck->GetYaxis()->CenterTitle();
    hAxis_Ck->GetYaxis()->SetTitleOffset(0.75);
    hAxis_Ck->GetYaxis()->SetTitleSize(0.07);
    hAxis_Ck->GetYaxis()->SetRangeUser(0.5, 4);
    hAxis_Ck->GetYaxis()->SetLimits(0.5, 4);
    hAxis_Ck->Draw("");
    for(unsigned short usCh=0; usCh<NumChannels; usCh++){
        gCorrFun[usCh]->Draw("same,C");
    }
    UnitLine->Draw("same");
    myLegend->Draw("same");
    Info1->Draw("same");
    c_Ck->SaveAs(OutputFolder+"c_Ck.png");
    c_Ck->Write();

    TCanvas* cZoom_Ck = new TCanvas("cZoom_Ck", "cZoom_Ck", 1);
    cZoom_Ck->cd(0);
    cZoom_Ck->SetCanvasSize(1280, 720);
    cZoom_Ck->SetMargin(0.15,0.05,0.2,0.05);//lrbt
    cZoom_Ck->SetGrid();
    cZoom_Ck->cd(0);

    TH1F* hAxisZoom_Ck = new TH1F("hAxisZoom_Ck","hAxisZoom_Ck",NumMomBins,kMin,kMax);
    hAxisZoom_Ck->SetStats(false);
    hAxisZoom_Ck->SetTitle("");
    hAxisZoom_Ck->GetXaxis()->SetLabelSize(0.06);
    hAxisZoom_Ck->GetXaxis()->SetTitle("k (MeV)");
    hAxisZoom_Ck->GetXaxis()->CenterTitle();
    hAxisZoom_Ck->GetXaxis()->SetTitleOffset(1.15);
    hAxisZoom_Ck->GetXaxis()->SetLabelOffset(0.02);
    hAxisZoom_Ck->GetXaxis()->SetTitleSize(0.07);
    hAxisZoom_Ck->GetYaxis()->SetLabelSize(0.06);
    hAxisZoom_Ck->GetYaxis()->SetTitle("C(k)");
    hAxisZoom_Ck->GetYaxis()->CenterTitle();
    hAxisZoom_Ck->GetYaxis()->SetTitleOffset(0.75);
    hAxisZoom_Ck->GetYaxis()->SetTitleSize(0.07);
    hAxisZoom_Ck->GetYaxis()->SetRangeUser(0.8, 1.25);
    hAxisZoom_Ck->GetYaxis()->SetLimits(0.8, 1.25);
    hAxisZoom_Ck->Draw("");
    for(unsigned short usCh=0; usCh<NumChannels; usCh++){
        gCorrFun[usCh]->Draw("same,C");
    }
    UnitLine->Draw("same");
    myLegend->Draw("same");
    Info1->Draw("same");
    cZoom_Ck->SaveAs(OutputFolder+"cZoom_Ck.png");
    cZoom_Ck->Write();
//////////////////////////////////////////////////////////////////////////////

    const unsigned NumRadBins=1024;
    const double RadMin=0;
    const double RadMax=12;
    const double RadStepSize=(RadMax-RadMin)/double(NumRadBins);
    double Radius;
    TH1F*** hSwave = new TH1F** [NumChannels];
    TH1F*** hWf2 = new TH1F** [NumChannels];
    for(unsigned short usCh=0; usCh<NumChannels; usCh++){
        hSwave[usCh] = new TH1F* [NumMomBins];
        hWf2[usCh] = new TH1F* [NumMomBins];
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            hSwave[usCh][uMomBin] = new TH1F("hSwave"+ChannelName[usCh]+TString::Format("_%.0f",Kitty.GetMomentum(uMomBin)),
                                             "hSwave"+ChannelName[usCh]+TString::Format("_%.0f",Kitty.GetMomentum(uMomBin)),
                                             NumRadBins,RadMin,RadMax);
            hWf2[usCh][uMomBin] = new TH1F("hWf2"+ChannelName[usCh]+TString::Format("_%.0f",Kitty.GetMomentum(uMomBin)),
                                             "hWf2"+ChannelName[usCh]+TString::Format("_%.0f",Kitty.GetMomentum(uMomBin)),
                                             NumRadBins,RadMin,RadMax);
            for(unsigned uRad=0; uRad<NumRadBins; uRad++){
                Radius = hSwave[usCh][uMomBin]->GetBinCenter(uRad+1);
                hSwave[usCh][uMomBin]->SetBinContent(uRad+1,real(Kitty.EvalRadialWaveFunction(uMomBin,usCh,0,Radius,true)));
                hWf2[usCh][uMomBin]->SetBinContent(uRad+1,abs(Kitty.EvalWaveFun2(uMomBin,Radius,usCh)));
            }
            hSwave[usCh][uMomBin]->Write();
            hWf2[usCh][uMomBin]->Write();
        }
    }

    for(unsigned short usCh=0; usCh<NumChannels; usCh++){
        delete hPhaseShift[usCh];
        delete fPhaseShift[usCh];
        delete gDelta[usCh];
        delete gCorrFun[usCh];
        delete hSwave[usCh];
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            delete hWf2[usCh][uMomBin];
        }
        delete [] hWf2[usCh];
    }
    delete [] hPhaseShift;
    delete [] fPhaseShift;
    delete [] gDelta;
    delete [] gCorrFun;
    delete [] hSwave;
    delete [] hWf2;

    delete hAxis;
    delete hAxisZoom;

    delete hAxis_Ck;
    delete hAxisZoom_Ck;
    delete UnitLine;
    delete ZeroLine;

    delete myLegend;
    delete Info1;

    delete cDelta;
    delete cDeltaZoom;

    delete fOut;
}


//we should have as input:
//the effective tau and masses for p, Lambda and Xi resonances
//the weights of the different mT (or kT or whatever we take) bins
//the radii for the different mT bins
//the output:
//the full S(r) for the pp, pL systems (mT integrated)
//attempt to fit the full S(r) with a Gaussian and see if there is a difference
//is the resulting Gaussian the same in pp and pL ?
void CompareDoubleGaussToGaussExp(){

    enum eSyst { pp, pL, pXim, LL };

    const TString OutputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/CompareDoubleGaussToGaussExp/";
    const TString InputFile = "";
    const TString FitName = "";

    const unsigned NumSystems = 4;
    TString* SystemName = new TString [NumSystems];
    SystemName[0] = "p-p";
    SystemName[1] = "p-#Lambda";
    SystemName[2] = "p-#Xi^{#minus}";
    SystemName[3] = "#Lambda-#Lambda";

    TString* ShortName = new TString [NumSystems];
    ShortName[0] = "pp";
    ShortName[1] = "pL";
    ShortName[2] = "pXi";
    ShortName[3] = "LL";

//Revision in November 2018
    const double ResMass_p = 1361.52;
    const double ResMass_L = 1462.93;
    const double ResMass_Xim = 1823;    // no idea
    const double ResTau_p = 1.65;
    const double ResTau_L = 4.69;
    const double ResTau_Xim = 8.220833333;    // no idea
    const double Prim_p = 0.3578;
    const double Prim_L = 0.3562;
    const double Prim_Xim = 0.9545;   // no idea

    const double Mass_p = 938.272;
    const double Mass_pi = 138; //kind of avg between pi+ and pi0
    const double Mass_L = 1115.683;
    const double Mass_Xim = 1321.7;

    //[system][paricle1/2]
    double** ResMass = new double* [NumSystems];
    double** ResTau = new double* [NumSystems];
    double** PrimFrac = new double* [NumSystems];
    double** Mass0 = new double* [NumSystems];
    double** Mass1 = new double* [NumSystems];
    for(unsigned uSyst=0; uSyst<NumSystems; uSyst++){
        ResMass[uSyst] = new double [2];
        ResTau[uSyst] = new double [2];
        PrimFrac[uSyst] = new double [2];
        Mass0[uSyst] = new double [2];
        Mass1[uSyst] = new double [2];
    }
    ResMass[pp][0] = ResMass_p;
    ResTau[pp][0] = ResTau_p;
    PrimFrac[pp][0] = Prim_p;
    Mass0[pp][0] = Mass_p;
    Mass1[pp][0] = Mass_pi;
    ResMass[pp][1] = ResMass_p;
    ResTau[pp][1] = ResTau_p;
    PrimFrac[pp][1] = Prim_p;
    Mass0[pp][1] = Mass_p;
    Mass1[pp][1] = Mass_pi;

    ResMass[pL][0] = ResMass_p;
    ResTau[pL][0] = ResTau_p;
    PrimFrac[pL][0] = Prim_p;
    Mass0[pL][0] = Mass_p;
    Mass1[pL][0] = Mass_pi;
    ResMass[pL][1] = ResMass_L;
    ResTau[pL][1] = ResTau_L;
    PrimFrac[pL][1] = Prim_L;
    Mass0[pL][1] = Mass_L;
    Mass1[pL][1] = Mass_pi;

    ResMass[pXim][0] = ResMass_p;
    ResTau[pXim][0] = ResTau_p;
    PrimFrac[pXim][0] = Prim_p;
    Mass0[pXim][0] = Mass_p;
    Mass1[pXim][0] = Mass_pi;
    ResMass[pXim][1] = ResMass_Xim;
    ResTau[pXim][1] = ResTau_Xim;
    PrimFrac[pXim][1] = Prim_Xim;
    Mass0[pXim][1] = Mass_Xim;
    Mass1[pXim][1] = Mass_pi;

    ResMass[LL][0] = ResMass_L;
    ResTau[LL][0] = ResTau_L;
    PrimFrac[LL][0] = Prim_L;
    Mass0[LL][0] = Mass_L;
    Mass1[LL][0] = Mass_pi;
    ResMass[LL][1] = ResMass_L;
    ResTau[LL][1] = ResTau_L;
    PrimFrac[LL][1] = Prim_L;
    Mass0[LL][1] = Mass_L;
    Mass1[LL][1] = Mass_pi;

    unsigned* Num_mT = new unsigned [NumSystems];
    Num_mT[pp] = 2;
    Num_mT[pL] = 2;
    Num_mT[pXim] = 2;
    Num_mT[LL] = 2;

    double** Mean_mT = new double* [NumSystems];
    double** Weight_mT = new double* [NumSystems];
    double** Radius_mT = new double* [NumSystems];
    for(unsigned uSyst=0; uSyst<NumSystems; uSyst++){
        Mean_mT[uSyst] = new double [Num_mT[uSyst]];
        Weight_mT[uSyst] = new double [Num_mT[uSyst]];
        Radius_mT[uSyst] = new double [Num_mT[uSyst]];
    }
    //!SET THE Mean_mT and Weight_mT values
    //Values from Andi on 19th September 2018

    Mean_mT[0][0] = 1;   Weight_mT[0][0] = 0.563;
    Mean_mT[0][1] = 2;   Weight_mT[0][1] = 0.437;

    Mean_mT[1][0] = 1;   Weight_mT[1][0] = 0.563;
    Mean_mT[1][1] = 2;   Weight_mT[1][1] = 0.437;

    Mean_mT[2][0] = 1;   Weight_mT[2][0] = 0.563;
    Mean_mT[2][1] = 2;   Weight_mT[2][1] = 0.437;

    Mean_mT[3][0] = 1;   Weight_mT[3][0] = 0.563;
    Mean_mT[3][1] = 2;   Weight_mT[3][1] = 0.437;

    const unsigned NumRadBins = 256;
    const double rMin = 0;
    const double rMax = 8;
    //used for the fake random source
    const double RelError = 0.000;
    const double AbsError = 0.001;

    TFile* InFile = new TFile(InputFile,"read");
    TF1* f_mT_R = InFile?(TF1*)InFile->Get(FitName):NULL;
    bool Del_f_mT_R=false;
    double fun0=0.515*1.00;
    double fun1=1.105*1.00;
    double xVal0=1;
    double xVal1=2;
    double parA = (fun1-fun0)/(xVal1-xVal0);
    double parB = (fun0*xVal1-fun1*xVal0)/(xVal1-xVal0);
    if(!f_mT_R){
        printf("\033[1;31mERROR:\033[0m The file %s could not be opened!\n",InputFile.Data());
        printf(" Initializing f_mT_R with some default values.\n");
        f_mT_R = new TF1("f_mT_R","[0]+x*[1]",rMin,rMax);
        //November revision
        f_mT_R->FixParameter(0,parA);
        f_mT_R->FixParameter(1,parB);
        Del_f_mT_R=true;
    }

    for(unsigned uSyst=0; uSyst<NumSystems; uSyst++){
        for(unsigned umT=0; umT<Num_mT[uSyst]; umT++){
            Radius_mT[uSyst][umT] = f_mT_R->Eval(Mean_mT[uSyst][umT]);
        }
    }

    //[System][kT]
    TH1F*** hkT_Source = new TH1F** [NumSystems];
    //the fit with a Gauss
    TF1*** fkT_GaussSource = new TF1** [NumSystems];;
    //[System]
    TH1F** h_Source = new TH1F* [NumSystems];;
    //the fit with a Gauss
    TF1** f_GaussSource = new TF1* [NumSystems];;
    //the individual Gauss source added up according to the weights
    //TH1F** h_AddedGauss = new TH1F* [NumSystems];

    double Radius;
    double S_Val;
    double* Pars_GETS2 = new double [14];

    for(unsigned uSyst=0; uSyst<NumSystems; uSyst++){
        h_Source[uSyst] = new TH1F(TString::Format("h_Source_%s",ShortName[uSyst].Data()),TString::Format("h_Source_%s",ShortName[uSyst].Data()),NumRadBins, rMin, rMax);
        f_GaussSource[uSyst] = new TF1(TString::Format("f_GaussSource_%s",ShortName[uSyst].Data()),GaussSourceTF1,rMin,rMax,4);
        f_GaussSource[uSyst]->FixParameter(0,0);f_GaussSource[uSyst]->FixParameter(1,0);f_GaussSource[uSyst]->FixParameter(2,0);
        f_GaussSource[uSyst]->SetParameter(3,1.0);
        f_GaussSource[uSyst]->SetParLimits(3,0.2,1.6);
        hkT_Source[uSyst] = new TH1F* [Num_mT[uSyst]];
        fkT_GaussSource[uSyst] = new TF1* [Num_mT[uSyst]];
        for(unsigned umT=0; umT<Num_mT[uSyst]; umT++){
            hkT_Source[uSyst][umT] = new TH1F(TString::Format("hkT_Source_%s_%u",ShortName[uSyst].Data(),umT),TString::Format("h_Source_%s_%u",ShortName[uSyst].Data(),umT),NumRadBins,rMin,rMax);
            fkT_GaussSource[uSyst][umT] = new TF1(TString::Format("fkT_GaussSource_%s_%u",ShortName[uSyst].Data(),umT),GaussSourceTF1,rMin,rMax,4);
            fkT_GaussSource[uSyst][umT]->FixParameter(0,0);fkT_GaussSource[uSyst][umT]->FixParameter(1,0);fkT_GaussSource[uSyst][umT]->FixParameter(2,0);
            fkT_GaussSource[uSyst][umT]->SetParameter(3,1.0);
            fkT_GaussSource[uSyst][umT]->SetParLimits(3,0.2,1.6);
            for(unsigned uRad=0; uRad<NumRadBins; uRad++){
                Radius = hkT_Source[uSyst][umT]->GetBinCenter(uRad+1);
                Pars_GETS2[1] = Radius;
                Pars_GETS2[3] = Radius_mT[uSyst][umT];
                Pars_GETS2[4] = ResTau[uSyst][0];
                Pars_GETS2[5] = PrimFrac[uSyst][0];
                Pars_GETS2[6] = ResMass[uSyst][0];
                Pars_GETS2[7] = Mass0[uSyst][0];
                Pars_GETS2[8] = Mass1[uSyst][0];
                Pars_GETS2[9] = ResTau[uSyst][1];
                Pars_GETS2[10] = PrimFrac[uSyst][1];
                Pars_GETS2[11] = ResMass[uSyst][1];
                Pars_GETS2[12] = Mass0[uSyst][1];
                Pars_GETS2[13] = Mass1[uSyst][1];
                S_Val = GaussExpTotSimple_2body(Pars_GETS2);

                hkT_Source[uSyst][umT]->SetBinContent(uRad+1,S_Val);
                hkT_Source[uSyst][umT]->SetBinError(uRad+1,RelError*S_Val+AbsError);

                h_Source[uSyst]->SetBinContent(uRad+1,h_Source[uSyst]->GetBinContent(uRad+1)+S_Val*Weight_mT[uSyst][umT]);
                h_Source[uSyst]->SetBinError(uRad+1,h_Source[uSyst]->GetBinError(uRad+1)+RelError*S_Val*Weight_mT[uSyst][umT]+AbsError);
            }
            hkT_Source[uSyst][umT]->Fit(fkT_GaussSource[uSyst][umT],"Q, S, N, R, M");
        }
        h_Source[uSyst]->Fit(f_GaussSource[uSyst],"S, N, R, M");
//new TF1("FitKitty",Dlm_AliceFitter4,0,1,18);
    }

    printf("------QA------\n");
    for(unsigned uSyst=0; uSyst<NumSystems; uSyst++){
        printf("In %s we have:\n", ShortName[uSyst].Data());
        for(unsigned umT=0; umT<Num_mT[uSyst]; umT++){
            printf(" mT=%.0f MeV; w=%.2f%%; R=%.2f fm\n",Mean_mT[uSyst][umT],Weight_mT[uSyst][umT]*100.,Radius_mT[uSyst][umT]);
        }
    }
    printf("\n");

    TFile* fOut = new TFile(OutputFolder+"fOut.root","recreate");
    for(unsigned uSyst=0; uSyst<NumSystems; uSyst++){
        for(unsigned umT=0; umT<Num_mT[uSyst]; umT++){
            hkT_Source[uSyst][umT]->Write();
            fkT_GaussSource[uSyst][umT]->Write();
        }
        h_Source[uSyst]->Write();
        f_GaussSource[uSyst]->Write();
    }

    TH1F* hAxisRad = new TH1F("hAxisRad", "hAxisRad", NumRadBins,rMin,rMax);
    hAxisRad->SetStats(false);
    hAxisRad->SetTitle("");
    hAxisRad->GetXaxis()->SetLabelSize(0.065);
    hAxisRad->GetXaxis()->SetTitle("r (fm)");
    hAxisRad->GetXaxis()->CenterTitle();
    hAxisRad->GetXaxis()->SetTitleOffset(1.2);
    hAxisRad->GetXaxis()->SetLabelOffset(0.02);
    hAxisRad->GetXaxis()->SetTitleSize(0.075);
    hAxisRad->GetYaxis()->SetLabelSize(0.065);
    hAxisRad->GetYaxis()->SetTitle("S(r) (fm^{-1})");
    hAxisRad->GetYaxis()->CenterTitle();
    hAxisRad->GetYaxis()->SetTitleOffset(1.00);
    hAxisRad->GetYaxis()->SetTitleSize(0.075);

    hAxisRad->GetYaxis()->SetRangeUser(0, 0.7);
    //hAxisRad->GetXaxis()->SetNdivisions(506);

    int COLORS[4] = {kBlue+2,kRed+1,kCyan+2,kOrange-1};
    TCanvas* cInfo1 = new TCanvas("cInfo1", "cInfo1", 1);
    cInfo1->cd(0); cInfo1->SetCanvasSize(1920, 1080); cInfo1->SetMargin(0.15,0.05,0.2,0.05);//lrbt
    hAxisRad->Draw();

    TPaveText* PT1 = new TPaveText(0.4,0.9,0.975,0.975, "blNDC");//lbrt
    PT1->SetName("PT1");
    PT1->SetBorderSize(1);
    PT1->SetTextSize(0.04);
    PT1->SetFillColor(kWhite);
    PT1->SetTextFont(22);
    PT1->AddText("The solid lines are Gaussian fits to the 'complete' source");

    TLegend* lInfo1 = new TLegend(0.7,0.6,0.975,0.9);//lbrt
    lInfo1->SetName(TString::Format("lInfo1"));
    lInfo1->SetTextSize(0.045);
    for(unsigned uSyst=0; uSyst<NumSystems; uSyst++){
        h_Source[uSyst]->SetLineColor(COLORS[uSyst]);
        h_Source[uSyst]->SetMarkerColor(COLORS[uSyst]);
        h_Source[uSyst]->SetMarkerSize(0);
        h_Source[uSyst]->SetLineWidth(2);
        f_GaussSource[uSyst]->SetLineColor(COLORS[uSyst]);
        f_GaussSource[uSyst]->SetMarkerColor(COLORS[uSyst]);
        f_GaussSource[uSyst]->SetLineWidth(4);
        lInfo1->AddEntry(f_GaussSource[uSyst],TString::Format("%s (R_{G} = %.2f fm)",SystemName[uSyst].Data(),f_GaussSource[uSyst]->GetParameter(3)));
        h_Source[uSyst]->Draw("same");
        f_GaussSource[uSyst]->Draw("same");
    }
    lInfo1->Draw("same");
    PT1->Draw("same");


    cInfo1->Write();
    cInfo1->SaveAs(OutputFolder+"cInfo1.png");

    if(Del_f_mT_R) delete f_mT_R;
    delete hAxisRad;
    delete cInfo1;
    delete fOut;
}

void TestCleverLevy1(){
double PARS_TEMP[5];
PARS_TEMP[1] = 0.007812;
PARS_TEMP[3] = 0.5;
PARS_TEMP[4] = 1.01;
printf("LevySource3D_2particle(PARS_TEMP)=%f\n",LevySource3D_2particle(PARS_TEMP));

//return;
    const TString OutputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/TestCleverLevy1/";

    const unsigned NumStabSteps = 21;
    const double StabMin = 1;
    const double StabMax = 2;
    const double StabStep = NumStabSteps>1?(StabMax-StabMin)/double(NumStabSteps-1):0;

    const unsigned NumScaleSteps = 31;
    const double ScaleMin = NumScaleSteps>1?0.5:1;
    const double ScaleMax = 1.5;
    const double ScaleStep = NumScaleSteps>1?(ScaleMax-ScaleMin)/double(NumScaleSteps-1):0;

    const unsigned NumBins = 256;
    const double rMin = 0;
    const double rMax = 16;
    TFile* fOutput = new TFile(OutputFolder+"fOutput.root","recreate");
    //this is the mse (mean square error) normalized to the mean, i.e. in percent
    TH2F* DeviationMap = new TH2F("DeviationMap","DeviationMap",NumScaleSteps,ScaleMin,ScaleMax,NumStabSteps,StabMin,StabMax);

    double Scale;
    double Stability;
    unsigned CLEV_STAB_BINS = 11;
    double CLEV_STAB_BW = (StabMax-StabMin)/double(CLEV_STAB_BINS-1);
    unsigned CLEV_SCAL_BINS = 11;
    double CLEV_SCAL_BW = (ScaleMax-ScaleMin)/double(CLEV_SCAL_BINS-1);
    unsigned CLEV_RAD_BINS = 64;
    double CLEV_RAD_BW = (rMax-rMin)/double(CLEV_RAD_BINS-1);
    DLM_CleverLevy CleverLevy;
    CleverLevy.InitStability(CLEV_STAB_BINS,StabMin,StabMax);
    CleverLevy.InitScale(CLEV_SCAL_BINS,ScaleMin,ScaleMax);
    CleverLevy.InitRad(CLEV_RAD_BINS,rMin,rMax);
    for(unsigned uScale=0; uScale<NumScaleSteps; uScale++){
printf("uScale=%u\n",uScale);
        Scale = ScaleMin+double(uScale)*ScaleStep;
        for(unsigned uStab=0; uStab<NumStabSteps; uStab++){
printf(" uStab=%u\n",uStab);
            Stability = StabMin+double(uStab)*StabStep;
            double PARS[5];
            TH1F* hAnaSource = new TH1F(TString::Format("hAnaSource_%.2f_%.2f",Scale,Stability), TString::Format("hAnaSource_%.2f_%.2f",Scale,Stability),
                                        NumBins, rMin, rMax);
            TH1F* hCleverSource = new TH1F(TString::Format("hCleverSource_%.2f_%.2f",Scale,Stability),
                                           TString::Format("hCleverSource_%.2f_%.2f",Scale,Stability),
                                        NumBins, rMin, rMax);
            double MEAN_DEV=0;
            double ValueAna;
            double ValueClever;
            for(unsigned uBin=0; uBin<NumBins; uBin++){
//printf("  uBin=%u\n",uBin);
                PARS[1] = hAnaSource->GetBinCenter(uBin+1);
                PARS[3] = Scale;
                PARS[4] = Stability;
                ValueAna = LevySource3D_2particle(PARS);
                ValueClever = CleverLevy.Eval(PARS);
double DEV = fabs(ValueClever-ValueAna)/ValueAna*100.;
//if(DEV>20){
//printf("DEV=%.1f(%.4f vs %.4f); r=%f; σ=%.2f; α=%.2f\n",DEV,ValueAna,ValueClever,PARS[1],Scale,Stability);
//}
                hAnaSource->SetBinContent(uBin+1,ValueAna);
                hCleverSource->SetBinContent(uBin+1,ValueClever);
                MEAN_DEV += fabs(ValueClever-ValueAna)/(0.5*(ValueClever+ValueAna))*100.;
            }
            MEAN_DEV /= double(NumBins);
            DeviationMap->SetBinContent(uScale+1,uStab+1,MEAN_DEV);
            hAnaSource->Write();
            hCleverSource->Write();
            delete hAnaSource;
            delete hCleverSource;
        }
    }


    DeviationMap->Write();

    delete DeviationMap;
    delete fOutput;
}

//test if single or pair distributions are the one corresponding to the standard
//definitions used in femtoscopy
void TestCleverLevy2(){
    const TString OutputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/TestCleverLevy2/";

    const unsigned NumStabSteps = 5;
    const double StabMin = 1;
    const double StabMax = 2;
    const double StabStep = NumStabSteps>1?(StabMax-StabMin)/double(NumStabSteps-1):0;

    const unsigned NumScaleSteps = 3;
    const double ScaleMin = NumScaleSteps>1?1.0:1;
    const double ScaleMax = 2.0;
    const double ScaleStep = NumScaleSteps>1?(ScaleMax-ScaleMin)/double(NumScaleSteps-1):0;

    const double kMin = 0;
    const double kMax = 280;
    const unsigned NumMomBins = 70;

    CATSparameters SourcePars(CATSparameters::tSource,2,true);
    SourcePars.SetParameter(0,0.5*(ScaleMax+ScaleMin));
    SourcePars.SetParameter(1,0.5*(StabMax+StabMin));

    DLM_CleverLevy CleverSingleLevy;
    CleverSingleLevy.InitStability(NumStabSteps,StabMin,StabMax);
    CleverSingleLevy.InitScale(NumScaleSteps,ScaleMin,ScaleMax);
    CleverSingleLevy.InitRad(256,0,64);
    CleverSingleLevy.InitType(0);

    DLM_CleverLevy CleverPairLevy;
    CleverPairLevy.InitStability(20,StabMin,StabMax);
    CleverPairLevy.InitScale(30,ScaleMin,ScaleMax);
    CleverPairLevy.InitRad(256,0,64);
    CleverPairLevy.InitType(1);

    CATS CatSingleLevy;
    CatSingleLevy.SetMaxNumThreads(1);
    CatSingleLevy.SetAnaSource(MemberSourceForwarder, &CleverSingleLevy, 2);
    CatSingleLevy.SetAnaSource(0,0.5*(ScaleMax+ScaleMin),false);
    CatSingleLevy.SetAnaSource(1,0.5*(StabMax+StabMin),false);
    CatSingleLevy.SetUseAnalyticSource(true);
    CatSingleLevy.SetThetaDependentSource(false);
    CatSingleLevy.SetMomentumDependentSource(false);
    CatSingleLevy.SetExcludeFailedBins(false);
    CatSingleLevy.SetMomBins(NumMomBins,kMin,kMax);
    CatSingleLevy.SetNumChannels(1);
    CatSingleLevy.SetSpin(0,0);
    CatSingleLevy.SetChannelWeight(0, 1);
    CatSingleLevy.SetQ1Q2(0);
    CatSingleLevy.SetPdgId(2212, 2212);
    const double Mass_pion = 134.9766;
    CatSingleLevy.SetRedMass( 0.5*Mass_pion );
    CatSingleLevy.KillTheCat();

    CATS CatPairLevy;
    CatPairLevy.SetMaxNumThreads(1);
    CatPairLevy.SetAnaSource(MemberSourceForwarder, &CleverPairLevy, 2);
    CatPairLevy.SetAnaSource(0,0.5*(ScaleMax+ScaleMin),false);
    CatPairLevy.SetAnaSource(1,0.5*(StabMax+StabMin),false);
    CatPairLevy.SetUseAnalyticSource(true);
    CatPairLevy.SetThetaDependentSource(false);
    CatPairLevy.SetMomentumDependentSource(false);
    CatPairLevy.SetExcludeFailedBins(false);
    CatPairLevy.SetMomBins(NumMomBins,kMin,kMax);
    CatPairLevy.SetNumChannels(1);
    CatPairLevy.SetSpin(0,0);
    CatPairLevy.SetChannelWeight(0, 1);
    CatPairLevy.SetQ1Q2(0);
    CatPairLevy.SetPdgId(2212, 2212);
    CatPairLevy.SetRedMass( 0.5*Mass_pion );
    CatPairLevy.KillTheCat();

    TFile* fOutput = new TFile(OutputFolder+"fOutput.root","recreate");

    double Scale;
    double Stability;
    CatSingleLevy.SetNotifications(CATS::nSilent);
    CatPairLevy.SetNotifications(CATS::nSilent);
    for(unsigned uScale=0; uScale<NumScaleSteps; uScale++){
        Scale = ScaleMin+double(uScale)*ScaleStep;
        for(unsigned uStab=0; uStab<NumStabSteps; uStab++){
            Stability = StabMin+double(uStab)*StabStep;
            printf("\r\033[K σ=%.2f; α=%.2f",Scale,Stability);
            cout << flush;

            CatSingleLevy.SetAnaSource(0,Scale,true);
            CatSingleLevy.SetAnaSource(1,Stability,true);
            CatSingleLevy.KillTheCat();
            CatPairLevy.SetAnaSource(0,Scale,true);
            CatPairLevy.SetAnaSource(1,Stability,true);
            CatPairLevy.KillTheCat();
            TGraph gSingleLevy;
            gSingleLevy.Set(NumMomBins);
            gSingleLevy.SetName(TString::Format("gSingleLevy_%.2f_%.2f",Scale,Stability));
            for(unsigned uBin=0; uBin<NumMomBins; uBin++){
                gSingleLevy.SetPoint(uBin,CatSingleLevy.GetMomentum(uBin),CatSingleLevy.GetCorrFun(uBin));
            }
            gSingleLevy.SetLineWidth(5);
            gSingleLevy.SetLineColor(kRed);
            TGraph gPairLevy;
            gPairLevy.Set(NumMomBins);
            gPairLevy.SetName(TString::Format("gPairLevy_%.2f_%.2f",Scale,Stability));
            for(unsigned uBin=0; uBin<NumMomBins; uBin++){
                gPairLevy.SetPoint(uBin,CatPairLevy.GetMomentum(uBin),CatPairLevy.GetCorrFun(uBin));
            }
            gPairLevy.SetLineWidth(4);
            gPairLevy.SetLineColor(kBlue);
            TGraph gAnaLevy;
            gAnaLevy.Set(NumMomBins);
            gAnaLevy.SetName(TString::Format("gAnaLevy_%.2f_%.2f",Scale,Stability));
            for(unsigned uBin=0; uBin<NumMomBins; uBin++){
                double Momentum = CatPairLevy.GetMomentum(uBin);
                gAnaLevy.SetPoint(uBin,Momentum,1.+exp(-pow(2.*Momentum*Scale*CatPairLevy.FmNu(),Stability)));
            }
            gAnaLevy.SetLineWidth(3);
            gAnaLevy.SetLineColor(kGreen+1);
            TGraph gAnaLevy2;
            gAnaLevy2.Set(NumMomBins);
            gAnaLevy2.SetName(TString::Format("gAnaLevy2_%.2f_%.2f",Scale,Stability));
            for(unsigned uBin=0; uBin<NumMomBins; uBin++){
                double Momentum = CatPairLevy.GetMomentum(uBin);
                gAnaLevy2.SetPoint(uBin,Momentum,1.+exp(-pow(2.*Momentum*Scale/sqrt(Stability)*CatPairLevy.FmNu(),Stability)));
            }
            gAnaLevy2.SetLineWidth(2.5);
            gAnaLevy2.SetLineColor(kViolet);
            gSingleLevy.Write();
            gPairLevy.Write();
            gAnaLevy.Write();
            gAnaLevy2.Write();
        }
    }
    delete fOutput;
}


void TestManyLevyParametrizations(const double& Stability, const double& Scale){
    DLM_Random RanGen(11);
    const unsigned NumBins=256;

    TH1F* hStable3D_0 = new TH1F("hStable3D_0","hStable3D_0",NumBins,0,16);
    hStable3D_0->SetLineColor(kRed);
    TH1F* hStable3D_1 = new TH1F("hStable3D_1","hStable3D_1",NumBins,0,16);
    hStable3D_1->SetLineColor(kBlue);
    TH1F* hGauss3D = new TH1F("hGauss3D","hGauss3D",NumBins,0,16);
    TH1F* hCauchy3D = new TH1F("hCauchy3D","hCauchy3D",NumBins,0,16);
    TGraph gAnaLevy_0;
    gAnaLevy_0.Set(NumBins);
    gAnaLevy_0.SetName("gAnaLevy_0");
    TGraph gAnaLevy_1;
    gAnaLevy_1.Set(NumBins);
    gAnaLevy_1.SetName("gAnaLevy_1");
    TGraph gAnaLevy_2;
    gAnaLevy_2.Set(NumBins);
    gAnaLevy_2.SetName("gAnaLevy_2");

    const unsigned NumIter = 400000;
    for(unsigned uIter=0; uIter<NumIter; uIter++){
        hStable3D_0->Fill(RanGen.StableR(3,Stability,0,Scale,0));
        hStable3D_1->Fill(RanGen.StableDiffR(3,Stability,0,Scale,0));
        hGauss3D->Fill(RanGen.GaussR(3,0,Scale));
        hCauchy3D->Fill(RanGen.CauchyR(3,0,Scale));
    }
    for(unsigned uBin=0; uBin<NumBins; uBin++){
        double PARS[5];
        PARS[1] = hStable3D_0->GetBinCenter(uBin+1);
        PARS[3] = Scale;
        PARS[4] = Stability;
        gAnaLevy_0.SetPoint(uBin,PARS[1],LevySource3D_single(PARS));
        gAnaLevy_1.SetPoint(uBin,PARS[1],LevySource3D_2particle(PARS));
        gAnaLevy_2.SetPoint(uBin,PARS[1],LevySource3D(PARS));
    }

    hStable3D_0->Scale(1./double(NumIter),"width");
    hStable3D_1->Scale(1./double(NumIter),"width");
    hGauss3D->Scale(1./double(NumIter),"width");
    hCauchy3D->Scale(1./double(NumIter),"width");

    printf("Integral hStable3D_0 = %f\n",hStable3D_0->Integral(1,hStable3D_0->GetNbinsX()));
    printf("Integral hStable3D_1 = %f\n",hStable3D_1->Integral(1,hStable3D_1->GetNbinsX()));
    printf("Integral hGauss3D = %f\n",hGauss3D->Integral(1,hGauss3D->GetNbinsX()));
    printf("Integral hCauchy3D = %f\n",hCauchy3D->Integral(1,hCauchy3D->GetNbinsX()));

    TFile* OutFile = new TFile(
        TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/TestManyLevyParametrizations/OutFile_%.2f_%.2f.root",Stability,Scale)
                               ,"recreate");
    hStable3D_0->Write();
    hStable3D_1->Write();
    hGauss3D->Write();
    hCauchy3D->Write();
    gAnaLevy_0.Write();
    gAnaLevy_1.Write();
    gAnaLevy_2.Write();

    delete hStable3D_0;
    delete hStable3D_1;
    delete hGauss3D;
    delete hCauchy3D;
    delete OutFile;
}
