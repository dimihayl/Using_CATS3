#include "PionAnalysis.h"

#include "CATS.h"
#include "CATSconstants.h"
#include "DLM_Random.h"
#include "DLM_Potentials.h"
#include "DLM_Source.h"
#include "DLM_CkDecomposition.h"
#include "DLM_Ck.h"
#include "TREPNI.h"
#include "CATStools.h"
#include "CECA.h"
#include "DLM_RootWrapper.h"

#include "CommonAnaFunctions.h"
#include "EnvVars.h"

#include "TString.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TPaveText.h"
#include "TNtuple.h"
#include "TVector3.h"
#include "TAxis.h"
#include "TStyle.h"



void pipi_test(){
    const unsigned NumMomBins = 100;
    const double kMin = 0;
    const double kMax = 400;
    CATS Kitty;
    Kitty.SetMomBins(NumMomBins,kMin,kMax);
    DLM_CommonAnaFunctions AnaObject;
    AnaObject.SetUpCats_pipi(Kitty,"McGauss_ResoTM",202);
    //AnaObject.SetUpCats_pipi(Kitty,"Gauss",0);
    //AnaObject.SetUpCats_pipi(Kitty,"Cauchy",0);
    Kitty.SetAnaSource(0,2.0);
    Kitty.SetAutoNormSource(false);
    Kitty.KillTheCat();

    TGraph g_pipi;
    g_pipi.SetName("g_pipi");
    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        g_pipi.SetPoint(uBin,Kitty.GetMomentum(uBin),Kitty.GetCorrFun(uBin));
    }

    TGraph gs_pipi;
    gs_pipi.SetName("gs_pipi");
    int uRad=0;
    for(float fRad=0; fRad<20; fRad+=0.1){
        uRad++;
        gs_pipi.SetPoint(uRad,fRad,Kitty.EvaluateTheSource(10,fRad,0));
    }

    TFile fOutput("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/Pions/pipi_alarm/pipi_test_noOmega.root","recreate");
    g_pipi.Write();
    gs_pipi.Write();

}

void pipi_alarm(){
    const unsigned NumMomBins = 100;
    const double kMin = 0;
    const double kMax = 400;
    const double SourceSize = 2.0;
    const double prim_fraction = 0.938;
    const double long_fraction = 0.055;
    const double omega_fraction = 0.075/1.0;
    const double purity = 0.99;
    //this is the lambda of primaries, including the omega
    //this is to be used if the omega is in the source
    const double prim_lambda = purity*prim_fraction*(1.-long_fraction)*purity*prim_fraction*(1.-long_fraction);
    //this is the lambda of the primaries, which do not include the lambda
    //this is to be used if the omega is to be part of the feed-down
    const double prim_omega_lambda = purity*prim_fraction*(1.-long_fraction-omega_fraction)*purity*prim_fraction*(1.-long_fraction-omega_fraction);
    const double misid_lambda = (1.-purity)*(1.-purity + purity + purity);
    const double feed_lambda = 1.-prim_lambda-misid_lambda;
    const double feed_omega_lambda = 1.-prim_omega_lambda-misid_lambda;

    printf("prim_lambda = %f\n",prim_lambda);
    printf("prim_omega_lambda = %f\n",prim_omega_lambda);
    printf("feed_lambda = %f\n",feed_lambda);
    printf("feed_omega_lambda = %f\n",feed_omega_lambda);
    printf("lambda + feed = %f vs %f\n",prim_lambda+feed_lambda,prim_omega_lambda+feed_omega_lambda);
    printf("misid_lambda = %f\n",misid_lambda);

    CATS KittyOld;
    KittyOld.SetMomBins(NumMomBins,kMin,kMax);
    DLM_CommonAnaFunctions AnaObjectBroken;
    AnaObjectBroken.SetUpCats_pipi_broken(KittyOld,"McGauss_ResoTM",202);
    //AnaObject.SetUpCats_pipi(Kitty,"Gauss",0);
    //AnaObject.SetUpCats_pipi(Kitty,"Cauchy",0);
    KittyOld.SetAnaSource(0,SourceSize);
    KittyOld.SetAutoNormSource(false);
    KittyOld.KillTheCat();

    TGraph g_pipiOld;
    g_pipiOld.SetName("g_pipiOld");
    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        g_pipiOld.SetPoint(uBin,KittyOld.GetMomentum(uBin),(prim_omega_lambda)*KittyOld.GetCorrFun(uBin)+(1.-prim_omega_lambda));
    }
    TGraph gs_pipiOld;
    gs_pipiOld.SetName("gs_pipiOld");
    int uRad=0;
    for(float fRad=0; fRad<20; fRad+=0.1){
        uRad++;
        gs_pipiOld.SetPoint(uRad,fRad,KittyOld.EvaluateTheSource(10,fRad,0));
    }

    CATS KittyNew;
    KittyNew.SetMomBins(NumMomBins,kMin,kMax);
    DLM_CommonAnaFunctions AnaObject;
    AnaObject.SetUpCats_pipi(KittyNew,"McGauss_ResoTM",204);
    //AnaObject.SetUpCats_pipi(Kitty,"Gauss",0);
    //AnaObject.SetUpCats_pipi(Kitty,"Cauchy",0);
    KittyNew.SetAnaSource(0,SourceSize);
    KittyNew.SetAutoNormSource(false);
    KittyNew.KillTheCat();

    TGraph g_pipiNew;
    g_pipiNew.SetName("g_pipiNew");
    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        g_pipiNew.SetPoint(uBin,KittyNew.GetMomentum(uBin),(prim_omega_lambda)*KittyNew.GetCorrFun(uBin)+(1.-prim_omega_lambda));
    }
    TGraph gs_pipiNew;
    gs_pipiNew.SetName("gs_pipiNew");
    TGraph gs_pipiNew_lambda;
    gs_pipiNew_lambda.SetName("gs_pipiNew_lambda");
    uRad=0;
    for(float fRad=0; fRad<20; fRad+=0.1){
        uRad++;
        gs_pipiNew.SetPoint(uRad,fRad,KittyNew.EvaluateTheSource(10,fRad,0));
        gs_pipiNew_lambda.SetPoint(uRad,fRad,KittyNew.EvaluateTheSource(10,fRad,0));
    }


    CATS KittyNew_omega;
    KittyNew_omega.SetMomBins(NumMomBins,kMin,kMax);
    DLM_CommonAnaFunctions AnaObject_omega;
    AnaObject_omega.SetUpCats_pipi(KittyNew_omega,"McGauss_ResoTM",202);
    //AnaObject.SetUpCats_pipi(Kitty,"Gauss",0);
    //AnaObject.SetUpCats_pipi(Kitty,"Cauchy",0);
    KittyNew_omega.SetAnaSource(0,SourceSize);
    KittyNew_omega.SetAutoNormSource(false);
    KittyNew_omega.KillTheCat();

    TGraph g_pipiNew_omega;
    g_pipiNew_omega.SetName("g_pipiNew_omega");
    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        g_pipiNew_omega.SetPoint(uBin,KittyNew_omega.GetMomentum(uBin),(prim_lambda)*KittyNew_omega.GetCorrFun(uBin)+(1.-prim_lambda));
    }
    TGraph gs_pipiNew_omega;
    gs_pipiNew_omega.SetName("gs_pipiNew_omega");
    uRad=0;
    for(float fRad=0; fRad<20; fRad+=0.1){
        uRad++;
        gs_pipiNew_omega.SetPoint(uRad,fRad,KittyNew_omega.EvaluateTheSource(10,fRad,0));
    }

    g_pipiOld.SetLineColor(kRed);
    g_pipiOld.SetLineWidth(2.5);
    gs_pipiOld.SetLineColor(kRed);
    gs_pipiOld.SetLineWidth(2.5);

    g_pipiNew.SetLineColor(kGreen+1);
    g_pipiNew.SetLineWidth(2.5);
    gs_pipiNew.SetLineColor(kGreen+1);
    gs_pipiNew.SetLineWidth(2.5);

    g_pipiNew_omega.SetLineColor(kAzure+1);
    g_pipiNew_omega.SetLineWidth(2.5);
    gs_pipiNew_omega.SetLineColor(kAzure+1);
    gs_pipiNew_omega.SetLineWidth(2.5);

    TFile fOutput("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/Pions/pipi_alarm/pipi_alarm.root","recreate");
    g_pipiOld.Write();
    gs_pipiOld.Write();

    g_pipiNew.Write();
    gs_pipiNew.Write();

    g_pipiNew_omega.Write();
    gs_pipiNew_omega.Write();

}

//constant between the two pars
double BoxySource(double* Pars){
    double& Radius = Pars[1];
    double& FROM = Pars[3];
    double& UNTIL = Pars[4];
    return 1./(UNTIL-FROM);
}

void ProveWeDoNotNeedOmega(){
    //A) a gauss core without any resonances (so Gauss), but include 50% of lambda
    //B) a gauss core with the resonances from max included in the source with a cTau of 20 fm, with 50% strength, lambda = 1

    const double SourceSize = 1.0;
    const double fraction = 0.1;
    const double lambda = fraction*fraction+2.*(1.-fraction)*fraction;
    printf("lambda = %f\n",lambda);
    //const double lambda = fraction;
    const unsigned NumMomBins = 100;
    const double kMin = 0;
    const double kMax = 200;

    CATS KittyLambda;
    CATSparameters cParsLambda(CATSparameters::tSource,1,true);
    cParsLambda.SetParameter(0,SourceSize);
    KittyLambda.SetMomBins(NumMomBins,kMin,kMax);

    KittyLambda.SetAnaSource(GaussSource, cParsLambda);
    KittyLambda.SetAnaSource(0, SourceSize);
    KittyLambda.SetUseAnalyticSource(true);
    KittyLambda.SetThetaDependentSource(false);
    KittyLambda.SetMomentumDependentSource(false);

    KittyLambda.SetExcludeFailedBins(false);
    KittyLambda.SetQ1Q2(1);
    KittyLambda.SetPdgId(211, 211);
    KittyLambda.SetRedMass( 0.5*Mass_pic );
    KittyLambda.SetNumChannels(1);
    KittyLambda.SetNumPW(0,0);
    KittyLambda.SetSpin(0,0);
    KittyLambda.SetChannelWeight(0, 1.);

    KittyLambda.KillTheCat();

    TGraph g_pipi_lambda;
    g_pipi_lambda.SetName("g_pipi_lambda");
    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        g_pipi_lambda.SetPoint(uBin,KittyLambda.GetMomentum(uBin),(1.-lambda)*KittyLambda.GetCorrFun(uBin)+lambda);
        //g_pipi_lambda.SetPoint(uBin,KittyLambda.GetMomentum(uBin),KittyLambda.GetCorrFun(uBin));
    }

    CATS BoxyKitty;
    CATSparameters cParsBoxy(CATSparameters::tSource,2,true);
    cParsBoxy.SetParameter(0,20);
    cParsBoxy.SetParameter(0,60);
    //BoxyKitty.SetSourceMaxRange(130);
    BoxyKitty.SetMomBins(NumMomBins,kMin,kMax);

    BoxyKitty.SetAnaSource(BoxySource, cParsBoxy);
    BoxyKitty.SetUseAnalyticSource(true);
    BoxyKitty.SetThetaDependentSource(false);
    BoxyKitty.SetMomentumDependentSource(false);

    BoxyKitty.SetExcludeFailedBins(false);
    BoxyKitty.SetQ1Q2(1);
    BoxyKitty.SetPdgId(211, 211);
    BoxyKitty.SetRedMass( 0.5*Mass_pic );
    BoxyKitty.SetNumChannels(1);
    BoxyKitty.SetNumPW(0,0);
    BoxyKitty.SetSpin(0,0);
    BoxyKitty.SetChannelWeight(0, 1.);

    BoxyKitty.KillTheCat();

    TGraph g_pipi_boxy;
    g_pipi_boxy.SetName("g_pipi_boxy");
    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        g_pipi_boxy.SetPoint(uBin,BoxyKitty.GetMomentum(uBin),BoxyKitty.GetCorrFun(uBin));
    }

    CATS KittyReso;
    CATSparameters cParsReso(CATSparameters::tSource,1,true);
    KittyReso.SetMomBins(NumMomBins,kMin,kMax);
    //KittyReso.SetSourceMaxRange(128);
    //KittyReso.SetMaxRad(128);

    DLM_CleverMcLevyResoTM CleverMcLevyResoTM;
    CleverMcLevyResoTM.InitStability(1,2-1e-6,2+1e-6);
    CleverMcLevyResoTM.InitScale(38,0.15,2.0);
    CleverMcLevyResoTM.InitRad(512,0,64);
    CleverMcLevyResoTM.InitType(2);
    CleverMcLevyResoTM.SetUpReso(0,fraction);
    CleverMcLevyResoTM.SetUpReso(1,fraction);
    const double k_CutOff = 200;
    Float_t k_D,fP1,fP2,fM1,fM2,Tau1,Tau2,AngleRcP1,AngleRcP2,AngleP1P2;
    DLM_Random RanGen(11);
    double RanVal1,RanVal2,RanVal3;
    TFile* F_EposDisto_p_pReso = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/Max/ForMax_pi_piReso.root");
    TNtuple* T_EposDisto_p_pReso = (TNtuple*)F_EposDisto_p_pReso->Get("InfoTuple_ClosePairs");
    unsigned N_EposDisto_p_pReso = T_EposDisto_p_pReso->GetEntries();
    T_EposDisto_p_pReso->SetBranchAddress("k_D",&k_D);T_EposDisto_p_pReso->SetBranchAddress("P1",&fP1);T_EposDisto_p_pReso->SetBranchAddress("P2",&fP2);
    T_EposDisto_p_pReso->SetBranchAddress("M1",&fM1);T_EposDisto_p_pReso->SetBranchAddress("M2",&fM2);T_EposDisto_p_pReso->SetBranchAddress("Tau1",&Tau1);
    T_EposDisto_p_pReso->SetBranchAddress("Tau2",&Tau2);T_EposDisto_p_pReso->SetBranchAddress("AngleRcP1",&AngleRcP1);
    T_EposDisto_p_pReso->SetBranchAddress("AngleRcP2",&AngleRcP2);T_EposDisto_p_pReso->SetBranchAddress("AngleP1P2",&AngleP1P2);
    for(unsigned uEntry=0; uEntry<N_EposDisto_p_pReso; uEntry++){
        T_EposDisto_p_pReso->GetEntry(uEntry);
        Tau1 = 0;
        Tau2 = 23;
        fM2 = 1124;
        if(k_D>k_CutOff) continue;
        RanVal1 = RanGen.Exponential(fM2/(fP2*Tau2));
        CleverMcLevyResoTM.AddBGT_PR(RanVal1,-cos(AngleRcP2));
        CleverMcLevyResoTM.AddBGT_RP(RanVal1,cos(AngleRcP2));
    }
    delete F_EposDisto_p_pReso;

    TFile* F_EposDisto_pReso_pReso = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/Max/ForMax_piReso_piReso.root");
    TNtuple* T_EposDisto_pReso_pReso = (TNtuple*)F_EposDisto_pReso_pReso->Get("InfoTuple_ClosePairs");
    unsigned N_EposDisto_pReso_pReso = T_EposDisto_pReso_pReso->GetEntries();
    T_EposDisto_pReso_pReso->SetBranchAddress("k_D",&k_D);T_EposDisto_pReso_pReso->SetBranchAddress("P1",&fP1);T_EposDisto_pReso_pReso->SetBranchAddress("P2",&fP2);
    T_EposDisto_pReso_pReso->SetBranchAddress("M1",&fM1);T_EposDisto_pReso_pReso->SetBranchAddress("M2",&fM2);T_EposDisto_pReso_pReso->SetBranchAddress("Tau1",&Tau1);
    T_EposDisto_pReso_pReso->SetBranchAddress("Tau2",&Tau2);T_EposDisto_pReso_pReso->SetBranchAddress("AngleRcP1",&AngleRcP1);
    T_EposDisto_pReso_pReso->SetBranchAddress("AngleRcP2",&AngleRcP2);T_EposDisto_pReso_pReso->SetBranchAddress("AngleP1P2",&AngleP1P2);
    for(unsigned uEntry=0; uEntry<N_EposDisto_pReso_pReso; uEntry++){
        T_EposDisto_pReso_pReso->GetEntry(uEntry);
        Tau1 = 23;
        Tau2 = 23;
        fM1 = 1124;
        fM2 = 1124;
        if(k_D>k_CutOff) continue;
        RanVal1 = RanGen.Exponential(fM1/(fP1*Tau1));
        RanVal2 = RanGen.Exponential(fM2/(fP2*Tau2));
        CleverMcLevyResoTM.AddBGT_RR(RanVal1,cos(AngleRcP1),RanVal2,cos(AngleRcP2),cos(AngleP1P2));
    }
    delete F_EposDisto_pReso_pReso;

    CleverMcLevyResoTM.InitNumMcIter(1000000);
    KittyReso.SetAnaSource(CatsSourceForwarder, &CleverMcLevyResoTM, 2);
    KittyReso.SetAnaSource(0,SourceSize);
    KittyReso.SetAnaSource(1,2.0);

    KittyReso.SetUseAnalyticSource(true);
    KittyReso.SetThetaDependentSource(false);
    KittyReso.SetMomentumDependentSource(false);
    KittyReso.SetAutoNormSource(false);

    KittyReso.SetExcludeFailedBins(false);
    KittyReso.SetQ1Q2(1);
    KittyReso.SetPdgId(211, 211);
    KittyReso.SetRedMass( 0.5*Mass_pic );
    KittyReso.SetNumChannels(1);
    KittyReso.SetNumPW(0,0);
    KittyReso.SetSpin(0,0);
    KittyReso.SetChannelWeight(0, 1.);


    KittyReso.KillTheCat();

    TH1F* hResoSource = new TH1F("hResoSource","hResoSource",512,0,64);
    TH1F* hResoSourceCats = new TH1F("hResoSourceCats","hResoSourceCats",512,0,64);
    TH1F* hResoSourceElder = new TH1F("hResoSourceElder","hResoSourceElder",512,0,64);
    for(unsigned uBin=0; uBin<hResoSource->GetNbinsX(); uBin++){
        double rad = hResoSource->GetBinCenter(uBin+1);
        double PARS[5];
        PARS[1] = rad;
        PARS[3] = KittyReso.GetAnaSourcePar(0);
        PARS[4] = 2;
        hResoSource->SetBinContent(uBin+1,CleverMcLevyResoTM.Eval(PARS));
        hResoSourceCats->SetBinContent(uBin+1,KittyReso.EvaluateTheSource(10,PARS[1],0));
    }

    double IntegralOfTheSource = hResoSource->Integral()*hResoSource->GetBinWidth(1);
    double IntegralOfTheSourceCats = hResoSourceCats->Integral()*hResoSourceCats->GetBinWidth(1);

    printf("IntegralOfTheSource = %f\n",IntegralOfTheSource);
    printf("IntegralOfTheSourceCats = %f\n",IntegralOfTheSourceCats);

    TGraph g_pipi_reso;
    g_pipi_reso.SetName("g_pipi_reso");
    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        g_pipi_reso.SetPoint(uBin,KittyReso.GetMomentum(uBin),KittyReso.GetCorrFun(uBin));
    }

    TFile fOutput("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/Pions/ProveWeDoNotNeedOmega/fOutput.root","recreate");
    g_pipi_lambda.Write();
    g_pipi_reso.Write();
    g_pipi_boxy.Write();
    hResoSource->Write();
    hResoSourceCats->Write();

    delete hResoSource;
    delete hResoSourceCats;
}


void pipi_omega(){
    const double SourceSize = 1.0;
    const double fraction = 0.1;
    const double lambda = fraction*fraction+2.*(1.-fraction)*fraction;
    //const double lambda = fraction;
    printf("lambda = %f\n",lambda);

    const unsigned NumMomBins = 100;
    const double kMin = 0;
    const double kMax = 200;

    CATS KittyLambda;
    KittyLambda.SetMomBins(NumMomBins,kMin,kMax);
    DLM_CommonAnaFunctions AnaObject;
    AnaObject.SetUpCats_pipi(KittyLambda,"McGauss_ResoTM",204);
    KittyLambda.SetAutoNormSource(false);
    KittyLambda.KillTheCat();

    TGraph g_pipi_lambda;
    g_pipi_lambda.SetName("g_pipi_lambda");
    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        g_pipi_lambda.SetPoint(uBin,KittyLambda.GetMomentum(uBin),(1.-lambda)*KittyLambda.GetCorrFun(uBin)+lambda);
    }

    TH1F* hSourceCats = new TH1F("hSourceCats","hSourceCats",512,0,64);
    for(unsigned uBin=0; uBin<hSourceCats->GetNbinsX(); uBin++){
        double rad = hSourceCats->GetBinCenter(uBin+1);
        double PARS[5];
        PARS[1] = rad;
        PARS[3] = KittyLambda.GetAnaSourcePar(0);
        PARS[4] = 2;
        hSourceCats->SetBinContent(uBin+1,KittyLambda.EvaluateTheSource(10,PARS[1],0));
    }

    double IntegralOfTheSourceCats = hSourceCats->Integral()*hSourceCats->GetBinWidth(1);

    CATS KittyReso;
    KittyReso.SetMomBins(NumMomBins,kMin,kMax);
    AnaObject.SetUpCats_pipi(KittyReso,"McGauss_ResoTM",202);
    KittyReso.SetAutoNormSource(false);
    KittyReso.KillTheCat();

    TGraph g_pipi_reso;
    g_pipi_reso.SetName("g_pipi_reso");
    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        g_pipi_reso.SetPoint(uBin,KittyReso.GetMomentum(uBin),KittyReso.GetCorrFun(uBin));
    }

    TH1F* hResoSourceCats = new TH1F("hResoSourceCats","hResoSourceCats",512,0,64);
    for(unsigned uBin=0; uBin<hResoSourceCats->GetNbinsX(); uBin++){
        double rad = hResoSourceCats->GetBinCenter(uBin+1);
        double PARS[5];
        PARS[1] = rad;
        PARS[3] = KittyReso.GetAnaSourcePar(0);
        PARS[4] = 2;
        hResoSourceCats->SetBinContent(uBin+1,KittyReso.EvaluateTheSource(10,PARS[1],0));
    }

    double IntegralOfTheResoSourceCats = hResoSourceCats->Integral()*hResoSourceCats->GetBinWidth(1);

    //printf("IntegralOfTheSource = %f\n",IntegralOfTheSource);
    printf("IntegralOfTheSourceCats = %f\n",IntegralOfTheSourceCats);
    printf("IntegralOfTheResoSourceCats = %f\n",IntegralOfTheResoSourceCats);

    TFile fOutput("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/Pions/pipi_omega/fOutput.root","recreate");
    g_pipi_lambda.Write();
    g_pipi_reso.Write();
    //hResoSource->Write();
    hSourceCats->Write();
    hResoSourceCats->Write();

    //delete hResoSource;
    delete hSourceCats;
    delete hResoSourceCats;
}

void EffectOnTailFor_pp_pL(){
    const unsigned NumMomBins = 100;
    const double kMin = 0;
    const double kMax = 400;

    CATS Kitty_pp;
    Kitty_pp.SetMomBins(NumMomBins,kMin,kMax);
    DLM_CommonAnaFunctions DCAF;
    DCAF.SetUpCats_pp(Kitty_pp,"AV18","McGauss_ResoTM",0,202);
    Kitty_pp.SetAutoNormSource(true);
    Kitty_pp.KillTheCat();
    TGraph g_pp_on;
    g_pp_on.SetName("g_pp_on");
    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        g_pp_on.SetPoint(uBin,Kitty_pp.GetMomentum(uBin),Kitty_pp.GetCorrFun(uBin));
    }

    CATS Kitty_ppoff;
    Kitty_ppoff.SetMomBins(NumMomBins,kMin,kMax);
    DCAF.SetUpCats_pp(Kitty_ppoff,"AV18","McGauss_ResoTM",0,202);
    Kitty_ppoff.SetAutoNormSource(false);
    Kitty_ppoff.KillTheCat();
    TGraph g_pp_off;
    g_pp_off.SetName("g_pp_off");
    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        g_pp_off.SetPoint(uBin,Kitty_ppoff.GetMomentum(uBin),Kitty_ppoff.GetCorrFun(uBin));
    }

    CATS Kitty_pL;
    Kitty_pL.SetMomBins(NumMomBins,kMin,kMax);
    DCAF.SetUpCats_pL(Kitty_pL,"NLO_Coupled_SPD","McGauss_ResoTM",0,202);
    Kitty_pL.SetAutoNormSource(true);
    Kitty_pL.KillTheCat();
    TGraph g_pL_on;
    g_pL_on.SetName("g_pL_on");
    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        g_pL_on.SetPoint(uBin,Kitty_pL.GetMomentum(uBin),Kitty_pL.GetCorrFun(uBin));
    }

    CATS Kitty_pLoff;
    Kitty_pLoff.SetMomBins(NumMomBins,kMin,kMax);
    DCAF.SetUpCats_pL(Kitty_pLoff,"NLO_Coupled_SPD","McGauss_ResoTM",0,202);
    Kitty_pLoff.SetAutoNormSource(true);
    Kitty_pLoff.KillTheCat();
    Kitty_pLoff.SetAutoNormSource(false);
    Kitty_pLoff.KillTheCat();
    TGraph g_pL_off;
    g_pL_off.SetName("g_pL_off");
    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        g_pL_off.SetPoint(uBin,Kitty_pLoff.GetMomentum(uBin),Kitty_pLoff.GetCorrFun(uBin));
    }

    TFile fOutput("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/Pions/EffectOnTailFor_pp_pL.root","recreate");
    g_pp_on.Write();
    g_pp_off.Write();
    g_pL_on.Write();
    g_pL_off.Write();
}









void SetUpReso_pipi_CUSTOM(CATS& Kitty, DLM_CleverMcLevyResoTM& CleverMcLevyResoTM, const double& cTau, const double& avgmass, const double& fraction, const double& fractionomega, const double& fractionlong, const bool& randomangle=false){

    Kitty.SetThetaDependentSource(false);

    CleverMcLevyResoTM.InitStability(1,2-1e-6,2+1e-6);
    CleverMcLevyResoTM.InitScale(38,0.15,2.0);
    CleverMcLevyResoTM.InitRad(512,0,64);
    CleverMcLevyResoTM.InitType(2);
    CleverMcLevyResoTM.SetUpReso(0,fraction+fractionlong);
    CleverMcLevyResoTM.SetUpReso(1,fraction+fractionlong);

    const double k_CutOff = 200;
    Float_t k_D,fP1,fP2,fM1,fM2,Tau1,Tau2,AngleRcP1,AngleRcP2,AngleP1P2;
    DLM_Random RanGen(11);
    double RanVal1,RanVal2,RanVal3;
    //TFile* F_EposDisto_p_pReso = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/Max/ForMax_pi_piReso_withOmega.root");
    TFile* F_EposDisto_p_pReso = new TFile(TString::Format("%s/CatsFiles/Source/EposAngularDist/ForMax_pi_piReso_withOmega.root",GetCernBoxDimi()));
    TNtuple* T_EposDisto_p_pReso = (TNtuple*)F_EposDisto_p_pReso->Get("InfoTuple_ClosePairs");
    unsigned N_EposDisto_p_pReso = T_EposDisto_p_pReso->GetEntries();
N_EposDisto_p_pReso = 10000;
    T_EposDisto_p_pReso->SetBranchAddress("k_D",&k_D);T_EposDisto_p_pReso->SetBranchAddress("P1",&fP1);T_EposDisto_p_pReso->SetBranchAddress("P2",&fP2);
    T_EposDisto_p_pReso->SetBranchAddress("M1",&fM1);T_EposDisto_p_pReso->SetBranchAddress("M2",&fM2);T_EposDisto_p_pReso->SetBranchAddress("Tau1",&Tau1);
    T_EposDisto_p_pReso->SetBranchAddress("Tau2",&Tau2);T_EposDisto_p_pReso->SetBranchAddress("AngleRcP1",&AngleRcP1);
    T_EposDisto_p_pReso->SetBranchAddress("AngleRcP2",&AngleRcP2);T_EposDisto_p_pReso->SetBranchAddress("AngleP1P2",&AngleP1P2);
    for(unsigned uEntry=0; uEntry<N_EposDisto_p_pReso; uEntry++){
///        if((fM2<782||fM2>783)&&RanGen.Uniform(0,1)<0.33) continue;//fix the omega yield
///        T_EposDisto_p_pReso->GetEntry(uEntry);
        Tau1 = 0;Tau2 = cTau;fM2 = avgmass;
        //put some fraction of longer lived stuff
        double RanNum = RanGen.Uniform(0,fraction+fractionomega+fractionlong);
        if(RanNum<fraction){}
        else if(RanNum<fraction+fractionomega){Tau2 = 23.24;fM2 = 782.6;}
        else {Tau2 = 1000;fM2 = avgmass;}
///        if(k_D>k_CutOff) continue;
        if(randomangle){
            TVector3 vRCORE(0,0,1.);TVector3 vSR2;
            vSR2.SetMagThetaPhi(1,asin(RanGen.Uniform(0,1)),RanGen.Uniform(0,2.*Pi));
            //vSR2.SetXYZ(0,0,-1);
            AngleRcP2 = vRCORE.Angle(vSR2);
            fP2 = 3000;
        }
        if(Tau2*fP2) RanVal1 = RanGen.Exponential(fM2/(fP2*Tau2)); else RanVal1=0;
fM2 = 1000;
//RanVal1 = fP2*Tau2/fM2;
        CleverMcLevyResoTM.AddBGT_PR(RanVal1,-cos(AngleRcP2));
        CleverMcLevyResoTM.AddBGT_RP(RanVal1,cos(AngleRcP2));
    }
    delete F_EposDisto_p_pReso;

    //TFile* F_EposDisto_pReso_pReso = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/Max/ForMax_piReso_piReso_withOmega.root");
    TFile* F_EposDisto_pReso_pReso = new TFile(TString::Format("%s/CatsFiles/Source/EposAngularDist/ForMax_piReso_piReso_withOmega.root",GetCernBoxDimi()));
    TNtuple* T_EposDisto_pReso_pReso = (TNtuple*)F_EposDisto_pReso_pReso->Get("InfoTuple_ClosePairs");
    unsigned N_EposDisto_pReso_pReso = T_EposDisto_pReso_pReso->GetEntries();
N_EposDisto_pReso_pReso = 10000;
    T_EposDisto_pReso_pReso->SetBranchAddress("k_D",&k_D);T_EposDisto_pReso_pReso->SetBranchAddress("P1",&fP1);T_EposDisto_pReso_pReso->SetBranchAddress("P2",&fP2);
    T_EposDisto_pReso_pReso->SetBranchAddress("M1",&fM1);T_EposDisto_pReso_pReso->SetBranchAddress("M2",&fM2);T_EposDisto_pReso_pReso->SetBranchAddress("Tau1",&Tau1);
    T_EposDisto_pReso_pReso->SetBranchAddress("Tau2",&Tau2);T_EposDisto_pReso_pReso->SetBranchAddress("AngleRcP1",&AngleRcP1);
    T_EposDisto_pReso_pReso->SetBranchAddress("AngleRcP2",&AngleRcP2);T_EposDisto_pReso_pReso->SetBranchAddress("AngleP1P2",&AngleP1P2);
    for(unsigned uEntry=0; uEntry<N_EposDisto_pReso_pReso; uEntry++){
///        if((fM1<782||fM1>783)&&RanGen.Uniform(0,1)<0.33) continue;//fix the omega yield
///        if((fM2<782||fM2>783)&&RanGen.Uniform(0,1)<0.33) continue;//fix the omega yield
///        T_EposDisto_pReso_pReso->GetEntry(uEntry);
        Tau1 = cTau; Tau2 = cTau; fM1 = avgmass; fM2 = avgmass;
        double RanNum = RanGen.Uniform(0,fraction+fractionomega+fractionlong);
        if(RanNum<fraction){}
        else if(RanNum<fraction+fractionomega){Tau1 = 23.24;fM1 = 782.6;}
        else {Tau1 = 1000;fM1 = avgmass;}
        RanNum = RanGen.Uniform(0,fraction+fractionomega+fractionlong);
        if(RanNum<fraction){}
        else if(RanNum<fraction+fractionomega){Tau2 = 23.24;fM2 = 782.6;}
        else {Tau2 = 1000;fM2 = avgmass;}
///        if(k_D>k_CutOff) continue;
        if(randomangle){
            TVector3 vRCORE;TVector3 vSR1;TVector3 vSR2;
            vRCORE.SetMagThetaPhi(1,asin(RanGen.Uniform(-1,1))+Pi*0.5,RanGen.Uniform(0,2.*Pi));
            vSR1.SetMagThetaPhi(1,asin(RanGen.Uniform(-1,1)),RanGen.Uniform(0,2.*Pi));
            vSR2.SetMagThetaPhi(1,asin(RanGen.Uniform(-1,1)),RanGen.Uniform(0,2.*Pi));
            //vSR1.SetXYZ(0,0,-1);
            //vSR2.SetXYZ(0,0,1);
            AngleRcP1 = vRCORE.Angle(vSR1);AngleRcP2 = vRCORE.Angle(vSR2);AngleP1P2 = vSR1.Angle(vSR2);
            fP1 = 3000;
            fP2 = 3000;
        }
        if(Tau1*fP1) RanVal1 = RanGen.Exponential(fM1/(fP1*Tau1)); else RanVal1=0;
        if(Tau2*fP2) RanVal2 = RanGen.Exponential(fM2/(fP2*Tau2)); else RanVal2=0;
fM1 = 1000;
fM2 = 1000;
//RanVal1 = fP1*Tau1/fM1;
//RanVal2 = fP2*Tau2/fM2;
        CleverMcLevyResoTM.AddBGT_RR(RanVal1,cos(AngleRcP1),RanVal2,cos(AngleRcP2),cos(AngleP1P2));
    }
    delete F_EposDisto_pReso_pReso;

    CleverMcLevyResoTM.InitNumMcIter(1000000);
    Kitty.SetAnaSource(CatsSourceForwarder, &CleverMcLevyResoTM, 2);
    Kitty.SetAnaSource(0,1.0);
    Kitty.SetAnaSource(1,2.0);
    Kitty.SetUseAnalyticSource(true);

    Kitty.SetMomentumDependentSource(false);
    Kitty.SetExcludeFailedBins(false);

    Kitty.SetQ1Q2(1);
    Kitty.SetPdgId(211, 211);
    Kitty.SetRedMass( 0.5*Mass_pic );

    Kitty.SetNumChannels(1);
    Kitty.SetNumPW(0,0);
    Kitty.SetSpin(0,0);
    Kitty.SetChannelWeight(0, 1.);

}



void pion_core_effect_plots(const TString FolderName, const TString NameAddOn){
printf("pion_core_effect_plots\n");
    TH1F* hAxis = new TH1F("hAxis","hAxis",128,0,160);
    hAxis->SetStats(false);
    hAxis->SetTitle("");
    hAxis->GetXaxis()->SetTitle("#it{k*} (MeV/#it{c})");
    hAxis->GetXaxis()->SetTitleSize(0.06);
    hAxis->GetXaxis()->SetLabelSize(0.06);
    hAxis->GetXaxis()->SetTitleOffset(1.3);
    hAxis->GetXaxis()->SetLabelOffset(0.02);

    hAxis->GetYaxis()->SetTitle("#it{C(k*)}");
    hAxis->GetYaxis()->SetTitleSize(0.06);
    hAxis->GetYaxis()->SetLabelSize(0.06);
    hAxis->GetYaxis()->SetTitleOffset(1.10);
    hAxis->GetYaxis()->SetRangeUser(0.95,1.85);

    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    TFile myfile(FolderName+"fOutput"+NameAddOn+".root","update");
    TGraph* gCore = (TGraph*)myfile.Get("KittyCore");
    TGraph* gShort = (TGraph*)myfile.Get("KittyShort");
    TGraph* gShortRho = (TGraph*)myfile.Get("KittyShortRho");
    TGraph* gNoLong = (TGraph*)myfile.Get("KittyNoLong");
    TGraph* gOmega = (TGraph*)myfile.Get("KittyOmega");
    TGraph* gOmegaLambda = (TGraph*)myfile.Get("KittyOmegaLambda");
    TGraph* gFull = (TGraph*)myfile.Get("KittyFull");
    TGraph* gLongLambda = (TGraph*)myfile.Get("KittyLongLambda");
    TGraph* gCauchy = (TGraph*)myfile.Get("KittyCauchy");

    gCore->SetLineColor(kBlack);
    gCore->SetLineStyle(2);
    gCore->SetLineWidth(5.0);

    gShort->SetLineColor(kViolet);
    gShort->SetLineStyle(3);
    gShort->SetLineWidth(4.0);

    gShortRho->SetLineColor(kBlue);
    gShortRho->SetLineStyle(4);
    gShortRho->SetLineWidth(4.0);

    gNoLong->SetLineColor(kCyan+1);
    gNoLong->SetLineStyle(8);
    gNoLong->SetLineWidth(4.0);

    gOmega->SetLineColor(kGreen+1);
    gOmega->SetLineStyle(9);
    gOmega->SetLineWidth(4.0);

    gOmegaLambda->SetLineColor(kGray+1);
    gOmegaLambda->SetLineStyle(6);
    gOmegaLambda->SetLineWidth(3.0);

    gFull->SetLineColor(kRed+1);
    gFull->SetLineStyle(1);
    gFull->SetLineWidth(4.0);

    gLongLambda->SetLineColor(kGray+1);
    gLongLambda->SetLineStyle(2);
    gLongLambda->SetLineWidth(3.0);

    gCauchy->SetLineColor(kGray+1);
    gCauchy->SetLineStyle(10);
    gCauchy->SetLineWidth(3.0);

    TCanvas* cCk_Core = new TCanvas("cCk_Core", "cCk_Core", 1);
    cCk_Core->cd(0);
    cCk_Core->SetCanvasSize(1280, 720);
    cCk_Core->SetMargin(0.15,0.05,0.2,0.05);//lrbt
    TLegend* l_Core = new TLegend(0.58,0.68,0.92,0.91);//lbrt
    l_Core->SetName(TString::Format("l_Core"));
    l_Core->SetTextSize(0.05);
    l_Core->SetBorderSize(0);
    l_Core->AddEntry(gCore,"r_{core} = 1.5 fm");
    hAxis->Draw("axis");
    gCore->Draw("C,same");
    l_Core->Draw("same");
    //gCore->Write();
    cCk_Core->Write();
    cCk_Core->SaveAs(FolderName+"cCk_Core"+NameAddOn+".pdf");
    cCk_Core->SaveAs(FolderName+"cCk_Core"+NameAddOn+".png");

    TCanvas* cCk_Short = new TCanvas("cCk_Short", "cCk_Short", 1);
    cCk_Short->cd(0);
    cCk_Short->SetCanvasSize(1280, 720);
    cCk_Short->SetMargin(0.15,0.05,0.2,0.05);//lrbt
    TLegend* l_Short = new TLegend(0.58,0.68,0.92,0.91);//lbrt
    l_Short->SetName(TString::Format("l_Short"));
    l_Short->SetTextSize(0.05);
    l_Short->SetBorderSize(0);
    l_Short->AddEntry(gCore,"r_{core} = 1.5 fm");
    l_Short->AddEntry(gShort,"r_{core} #oplus c#tau < 1 fm");
    hAxis->Draw("axis");
    gCore->Draw("C,same");
    gShort->Draw("C,same");
    l_Short->Draw("same");
    //gShort->Write();
    cCk_Short->Write();
    cCk_Short->SaveAs(FolderName+"cCk_Short"+NameAddOn+".pdf");
    cCk_Short->SaveAs(FolderName+"cCk_Short"+NameAddOn+".png");

    TCanvas* cCk_ShortRho = new TCanvas("cCk_ShortRho", "cCk_ShortRho", 1);
    cCk_ShortRho->cd(0);
    cCk_ShortRho->SetCanvasSize(1280, 720);
    cCk_ShortRho->SetMargin(0.15,0.05,0.2,0.05);//lrbt
    TLegend* l_ShortRho = new TLegend(0.58,0.68,0.92,0.91);//lbrt
    l_ShortRho->SetName(TString::Format("l_ShortRho"));
    l_ShortRho->SetTextSize(0.05);
    l_ShortRho->SetBorderSize(0);
    l_ShortRho->AddEntry(gCore,"r_{core} = 1.5 fm");
    l_ShortRho->AddEntry(gShortRho,"r_{core} #oplus c#tau #leq c#tau(#rho)");
    hAxis->Draw("axis");
    gCore->Draw("C,same");
    gShortRho->Draw("C,same");
    l_ShortRho->Draw("same");
    //gShortRho->Write();
    cCk_ShortRho->Write();
    cCk_ShortRho->SaveAs(FolderName+"cCk_ShortRho"+NameAddOn+".pdf");
    cCk_ShortRho->SaveAs(FolderName+"cCk_ShortRho"+NameAddOn+".png");

    TCanvas* cCk_NoLong = new TCanvas("cCk_NoLong", "cCk_NoLong", 1);
    cCk_NoLong->cd(0);
    cCk_NoLong->SetCanvasSize(1280, 720);
    cCk_NoLong->SetMargin(0.15,0.05,0.2,0.05);//lrbt
    TLegend* l_NoLong = new TLegend(0.58,0.68,0.92,0.91);//lbrt
    l_NoLong->SetName(TString::Format("l_NoLong"));
    l_NoLong->SetTextSize(0.05);
    l_NoLong->SetBorderSize(0);
    l_NoLong->AddEntry(gCore,"r_{core} = 1.5 fm");
    l_NoLong->AddEntry(gNoLong,"r_{core} #oplus c#tau < 5 fm");
    hAxis->Draw("axis");
    gCore->Draw("C,same");
    gNoLong->Draw("C,same");
    l_NoLong->Draw("same");
    //gNoLong->Write();
    cCk_NoLong->Write();
    cCk_NoLong->SaveAs(FolderName+"cCk_NoLong"+NameAddOn+".pdf");
    cCk_NoLong->SaveAs(FolderName+"cCk_NoLong"+NameAddOn+".png");

    TCanvas* cCk_Omega = new TCanvas("cCk_Omega", "cCk_Omega", 1);
    cCk_Omega->cd(0);
    cCk_Omega->SetCanvasSize(1280, 720);
    cCk_Omega->SetMargin(0.15,0.05,0.2,0.05);//lrbt
    TLegend* l_Omega = new TLegend(0.58,0.68,0.92,0.91);//lbrt
    l_Omega->SetName(TString::Format("l_Omega"));
    l_Omega->SetTextSize(0.05);
    l_Omega->SetBorderSize(0);
    l_Omega->AddEntry(gCore,"r_{core} = 1.5 fm");
    l_Omega->AddEntry(gOmega,"r_{core} #oplus c#tau #leq c#tau(#omega)");
    hAxis->Draw("axis");
    gCore->Draw("C,same");
    gOmega->Draw("C,same");
    l_Omega->Draw("same");
    //gOmega->Write();
    cCk_Omega->Write();
    cCk_Omega->SaveAs(FolderName+"cCk_Omega"+NameAddOn+".pdf");
    cCk_Omega->SaveAs(FolderName+"cCk_Omega"+NameAddOn+".png");

    TCanvas* cCk_OmegaLambda = new TCanvas("cCk_OmegaLambda", "cCk_OmegaLambda", 1);
    cCk_OmegaLambda->cd(0);
    cCk_OmegaLambda->SetCanvasSize(1280, 720);
    cCk_OmegaLambda->SetMargin(0.15,0.05,0.2,0.05);//lrbt
    TLegend* l_OmegaLambda = new TLegend(0.58,0.68,0.92,0.91);//lbrt
    l_OmegaLambda->SetName(TString::Format("l_OmegaLambda"));
    l_OmegaLambda->SetTextSize(0.05);
    l_OmegaLambda->SetBorderSize(0);
    l_OmegaLambda->AddEntry(gCore,"r_{core} = 1.5 fm");
    l_OmegaLambda->AddEntry(gOmega,"r_{core} #oplus c#tau #leq c#tau(#omega)");
    l_OmegaLambda->AddEntry(gOmegaLambda,"Flat #omega contribution");
    hAxis->Draw("axis");
    gCore->Draw("C,same");
    gOmega->Draw("C,same");
    gOmegaLambda->Draw("C,same");
    l_OmegaLambda->Draw("same");
    //gOmegaLambda->Write();
    cCk_OmegaLambda->Write();
    cCk_OmegaLambda->SaveAs(FolderName+"cCk_OmegaLambda"+NameAddOn+".pdf");
    cCk_OmegaLambda->SaveAs(FolderName+"cCk_OmegaLambda"+NameAddOn+".png");

    TCanvas* cCk_Full = new TCanvas("cCk_Full", "cCk_Full", 1);

    cCk_Full->cd(0);
    cCk_Full->SetCanvasSize(1280, 720);
    cCk_Full->SetMargin(0.15,0.05,0.2,0.05);//lrbt
    TLegend* l_Full = new TLegend(0.58,0.68,0.92,0.91);//lbrt
    l_Full->SetName(TString::Format("l_Full"));
    l_Full->SetTextSize(0.05);
    l_Full->SetBorderSize(0);
    l_Full->AddEntry(gCore,"#it{r}_{core} = 1.5 fm");
    l_Full->AddEntry(gFull,"#it{r}_{core} #oplus all resonances");
    hAxis->Draw("axis");
    gCore->Draw("C,same");
    gFull->Draw("C,same");
    l_Full->Draw("same");

    //gFull->Write();
    cCk_Full->Write();
    cCk_Full->SaveAs(FolderName+"cCk_Full"+NameAddOn+".pdf");
    cCk_Full->SaveAs(FolderName+"cCk_Full"+NameAddOn+".png");

    TCanvas* cCk_LongLambda = new TCanvas("cCk_LongLambda", "cCk_LongLambda", 1);
    cCk_LongLambda->cd(0);
    cCk_LongLambda->SetCanvasSize(1280, 720);
    cCk_LongLambda->SetMargin(0.15,0.05,0.2,0.05);//lrbt
    TLegend* l_LongLambda = new TLegend(0.58,0.68,0.92,0.91);//lbrt
    l_LongLambda->SetName(TString::Format("l_LongLambda"));
    l_LongLambda->SetTextSize(0.05);
    l_LongLambda->SetBorderSize(0);
    l_LongLambda->AddEntry(gCore,"#it{r}_{core} = 1.5 fm");
    l_LongLambda->AddEntry(gFull,"#it{r}_{core} #oplus all resonances");
    l_LongLambda->AddEntry(gLongLambda,"c#tau > c#tau(#omega) assumed flat");
    hAxis->Draw("axis");
    gCore->Draw("C,same");
    gFull->Draw("C,same");
    gLongLambda->Draw("C,same");
    l_LongLambda->Draw("same");
    //gLongLambda->Write();
    cCk_LongLambda->Write();
    cCk_LongLambda->SaveAs(FolderName+"cCk_LongLambda"+NameAddOn+".pdf");
    cCk_LongLambda->SaveAs(FolderName+"cCk_LongLambda"+NameAddOn+".png");

    TCanvas* cCk_Cauchy = new TCanvas("cCk_Cauchy", "cCk_Cauchy", 1);
    cCk_Cauchy->cd(0);
    cCk_Cauchy->SetCanvasSize(1280, 720);
    cCk_Cauchy->SetMargin(0.15,0.05,0.2,0.05);//lrbt
    TLegend* l_Cauchy = new TLegend(0.58,0.68,0.92,0.91);//lbrt
    l_Cauchy->SetName(TString::Format("l_Cauchy"));
    l_Cauchy->SetTextSize(0.05);
    l_Cauchy->SetBorderSize(0);
    l_Cauchy->AddEntry(gCore,"#it{r}_{core} = 1.5 fm");
    l_Cauchy->AddEntry(gFull,"#it{r}_{core} #oplus all resonances");
    l_Cauchy->AddEntry(gCauchy,"#it{r}_{exp} = 2.0 fm; #it{#lambda} = 0.9");
    hAxis->Draw("axis");
    gCore->Draw("C,same");
    gFull->Draw("C,same");
    gCauchy->Draw("C,same");
    l_Cauchy->Draw("same");
    //gCauchy->Write();
    cCk_Cauchy->Write();
    cCk_Cauchy->SaveAs(FolderName+"cCk_Cauchy"+NameAddOn+".pdf");
    cCk_Cauchy->SaveAs(FolderName+"cCk_Cauchy"+NameAddOn+".png");

    //even number is linear axis, odd number is log scale
    //0-1 is 4piS(r)
    //2-3 is S(r)
    for(int iLog=0; iLog<4; iLog++){
        TH1F* hSrAxis = new TH1F("hSrAxis","hSrAxis",128,0,15);
        hSrAxis->SetStats(false);
        hSrAxis->SetTitle("");
        hSrAxis->GetXaxis()->SetTitle("#it{r*} (fm)");
        hSrAxis->GetXaxis()->SetTitleSize(0.06);
        hSrAxis->GetXaxis()->SetLabelSize(0.06);
        hSrAxis->GetXaxis()->SetTitleOffset(1.3);
        hSrAxis->GetXaxis()->SetLabelOffset(0.02);

        if(iLog<=1) hSrAxis->GetYaxis()->SetTitle("4#it{#pir*}^{2}#it{S(r*)} (1/fm)");
        else hSrAxis->GetYaxis()->SetTitle("S(r*) (1/fm^{3})");
        hSrAxis->GetYaxis()->SetTitleSize(0.06);
        hSrAxis->GetYaxis()->SetLabelSize(0.06);
        hSrAxis->GetYaxis()->SetTitleOffset(1.10);
        if(iLog<=1) hSrAxis->GetYaxis()->SetRangeUser(iLog%2?0.005:0,0.3);
        else hSrAxis->GetYaxis()->SetRangeUser(iLog%2?0.000001:0,0.2);

        TGraph* gSrCore = (TGraph*)myfile.Get("SourceCore");
        TGraph* gSrShort = (TGraph*)myfile.Get("SourceShort");
        TGraph* gSrShortRho = (TGraph*)myfile.Get("SourceShortRho");
        TGraph* gSrNoLong = (TGraph*)myfile.Get("SourceNoLong");
        TGraph* gSrOmega = (TGraph*)myfile.Get("SourceOmega");
        TGraph* gSrOmegaLambda = (TGraph*)myfile.Get("SourceOmegaLambda");
        TGraph* gSrFull = (TGraph*)myfile.Get("SourceFull");
        TGraph* gSrLongLambda = (TGraph*)myfile.Get("SourceLongLambda");
        TGraph* gSrCauchy = (TGraph*)myfile.Get("SourceCauchy");
        if(iLog>=2){
            double xVal,yVal;
            unsigned NumPts = gSrCore->GetN();
            for(unsigned uPts=0; uPts<NumPts; uPts++){
                gSrCore->GetPoint(uPts,xVal,yVal);
                gSrCore->SetPoint(uPts,xVal,yVal/(4.*TMath::Pi()*xVal*xVal));

                gSrShort->GetPoint(uPts,xVal,yVal);
                gSrShort->SetPoint(uPts,xVal,yVal/(4.*TMath::Pi()*xVal*xVal));

                gSrShortRho->GetPoint(uPts,xVal,yVal);
                gSrShortRho->SetPoint(uPts,xVal,yVal/(4.*TMath::Pi()*xVal*xVal));

                gSrNoLong->GetPoint(uPts,xVal,yVal);
                gSrNoLong->SetPoint(uPts,xVal,yVal/(4.*TMath::Pi()*xVal*xVal));

                gSrOmega->GetPoint(uPts,xVal,yVal);
                gSrOmega->SetPoint(uPts,xVal,yVal/(4.*TMath::Pi()*xVal*xVal));

                gSrOmegaLambda->GetPoint(uPts,xVal,yVal);
                gSrOmegaLambda->SetPoint(uPts,xVal,yVal/(4.*TMath::Pi()*xVal*xVal));

                gSrFull->GetPoint(uPts,xVal,yVal);
                gSrFull->SetPoint(uPts,xVal,yVal/(4.*TMath::Pi()*xVal*xVal));

                gSrLongLambda->GetPoint(uPts,xVal,yVal);
                gSrLongLambda->SetPoint(uPts,xVal,yVal/(4.*TMath::Pi()*xVal*xVal));

                gSrCauchy->GetPoint(uPts,xVal,yVal);
                gSrCauchy->SetPoint(uPts,xVal,yVal/(4.*TMath::Pi()*xVal*xVal));
            }
        }

        gSrCore->SetLineColor(kBlack);
        gSrCore->SetLineStyle(2);
        gSrCore->SetLineWidth(5.0);

        gSrShort->SetLineColor(kViolet);
        gSrShort->SetLineStyle(3);
        gSrShort->SetLineWidth(4.0);

        gSrShortRho->SetLineColor(kBlue);
        gSrShortRho->SetLineStyle(4);
        gSrShortRho->SetLineWidth(4.0);

        gSrNoLong->SetLineColor(kCyan+1);
        gSrNoLong->SetLineStyle(8);
        gSrNoLong->SetLineWidth(4.0);

        gSrOmega->SetLineColor(kGreen+1);
        gSrOmega->SetLineStyle(9);
        gSrOmega->SetLineWidth(4.0);

        gSrOmegaLambda->SetLineColor(kGray+1);
        gSrOmegaLambda->SetLineStyle(6);
        gSrOmegaLambda->SetLineWidth(3.0);

        gSrFull->SetLineColor(kRed+1);
        gSrFull->SetLineStyle(1);
        gSrFull->SetLineWidth(4.0);

        gSrLongLambda->SetLineColor(kGray+1);
        gSrLongLambda->SetLineStyle(2);
        gSrLongLambda->SetLineWidth(3.0);

        gSrCauchy->SetLineColor(kGray+1);
        gSrCauchy->SetLineStyle(10);
        gSrCauchy->SetLineWidth(3.0);

        TString BaseCanvasName;
        if(iLog<=1)BaseCanvasName = iLog%2?"c4piSrLog_":"c4piSr_";
        else BaseCanvasName = iLog%2?"cSrLog_":"cSr_";
        BaseCanvasName += "Core";
        BaseCanvasName += NameAddOn;

        TString BCN = BaseCanvasName+"Core";
        TCanvas* cSr_Core = new TCanvas(BCN, BCN+"Core", 1);
        cSr_Core->cd(0);
        cSr_Core->SetCanvasSize(1280, 720);
        cSr_Core->SetMargin(0.15,0.05,0.2,0.05);//lrbt
        hSrAxis->Draw("axis");
        gSrCore->Draw("C,same");
        l_Core->Draw("same");
        if(iLog%2) cSr_Core->SetLogy();
        cSr_Core->Write();
        cSr_Core->SaveAs(FolderName+BCN+".pdf");
        cSr_Core->SaveAs(FolderName+BCN+".png");

        BCN = BaseCanvasName+"Short";
        TCanvas* cSr_Short = new TCanvas(BCN,BCN, 1);
        cSr_Short->cd(0);
        cSr_Short->SetCanvasSize(1280, 720);
        cSr_Short->SetMargin(0.15,0.05,0.2,0.05);//lrbt
        hSrAxis->Draw("axis");
        gSrCore->Draw("C,same");
        gSrShort->Draw("C,same");
        l_Short->Draw("same");
        if(iLog%2) cSr_Short->SetLogy();
        cSr_Short->Write();
        cSr_Short->SaveAs(FolderName+BCN+".pdf");
        cSr_Short->SaveAs(FolderName+BCN+".png");

        BCN = BaseCanvasName+"ShortRho";
        TCanvas* cSr_ShortRho = new TCanvas(BCN,BCN, 1);
        cSr_ShortRho->cd(0);
        cSr_ShortRho->SetCanvasSize(1280, 720);
        cSr_ShortRho->SetMargin(0.15,0.05,0.2,0.05);//lrbt
        hSrAxis->Draw("axis");
        gSrCore->Draw("C,same");
        gSrShortRho->Draw("C,same");
        l_ShortRho->Draw("same");
        if(iLog%2) cSr_ShortRho->SetLogy();
        cSr_ShortRho->Write();
        cSr_ShortRho->SaveAs(FolderName+BCN+".pdf");
        cSr_ShortRho->SaveAs(FolderName+BCN+".png");

        BCN = BaseCanvasName+"NoLong";
        TCanvas* cSr_NoLong = new TCanvas(BCN,BCN, 1);
        cSr_NoLong->cd(0);
        cSr_NoLong->SetCanvasSize(1280, 720);
        cSr_NoLong->SetMargin(0.15,0.05,0.2,0.05);//lrbt
        hSrAxis->Draw("axis");
        gSrCore->Draw("C,same");
        gSrNoLong->Draw("C,same");
        l_NoLong->Draw("same");
        if(iLog%2) cSr_NoLong->SetLogy();
        cSr_NoLong->Write();
        cSr_NoLong->SaveAs(FolderName+BCN+".pdf");
        cSr_NoLong->SaveAs(FolderName+BCN+".png");

        BCN = BaseCanvasName+"Omega";
        TCanvas* cSr_Omega = new TCanvas(BCN, BCN, 1);
        cSr_Omega->cd(0);
        cSr_Omega->SetCanvasSize(1280, 720);
        cSr_Omega->SetMargin(0.15,0.05,0.2,0.05);//lrbt
        hSrAxis->Draw("axis");
        gSrCore->Draw("C,same");
        gSrOmega->Draw("C,same");
        l_Omega->Draw("same");
        if(iLog%2) cSr_Omega->SetLogy();
        cSr_Omega->Write();
        cSr_Omega->SaveAs(FolderName+BCN+".pdf");
        cSr_Omega->SaveAs(FolderName+BCN+".png");

        BCN = BaseCanvasName+"OmegaLambda";
        TCanvas* cSr_OmegaLambda = new TCanvas(BCN,BCN, 1);
        cSr_OmegaLambda->cd(0);
        cSr_OmegaLambda->SetCanvasSize(1280, 720);
        cSr_OmegaLambda->SetMargin(0.15,0.05,0.2,0.05);//lrbt
        hSrAxis->Draw("axis");
        gSrCore->Draw("C,same");
        gSrOmega->Draw("C,same");
        gSrOmegaLambda->Draw("C,same");
        l_OmegaLambda->Draw("same");
        if(iLog%2) cSr_OmegaLambda->SetLogy();
        cSr_OmegaLambda->Write();
        cSr_OmegaLambda->SaveAs(FolderName+BCN+".pdf");
        cSr_OmegaLambda->SaveAs(FolderName+BCN+".png");

        BCN = BaseCanvasName+"Full";
        TCanvas* cSr_Full = new TCanvas(BCN,BCN, 1);
        cSr_Full->cd(0);
        cSr_Full->SetCanvasSize(1280, 720);
        cSr_Full->SetMargin(0.15,0.05,0.2,0.05);//lrbt
        hSrAxis->Draw("axis");
        gSrCore->Draw("C,same");
        gSrFull->Draw("C,same");
        l_Full->Draw("same");
        if(iLog%2) cSr_Full->SetLogy();
        cSr_Full->Write();
        cSr_Full->SaveAs(FolderName+BCN+".pdf");
        cSr_Full->SaveAs(FolderName+BCN+".png");

        BCN = BaseCanvasName+"LongLambda";
        TCanvas* cSr_LongLambda = new TCanvas(BCN,BCN, 1);
        cSr_LongLambda->cd(0);
        cSr_LongLambda->SetCanvasSize(1280, 720);
        cSr_LongLambda->SetMargin(0.15,0.05,0.2,0.05);//lrbt
        hSrAxis->Draw("axis");
        gSrCore->Draw("C,same");
        gSrFull->Draw("C,same");
        gSrLongLambda->Draw("C,same");
        l_LongLambda->Draw("same");
        if(iLog%2) cSr_LongLambda->SetLogy();
        cSr_LongLambda->Write();
        cSr_LongLambda->SaveAs(FolderName+BCN+".pdf");
        cSr_LongLambda->SaveAs(FolderName+BCN+".png");

        BCN = BaseCanvasName+"Cauchy";
        TCanvas* cSr_Cauchy = new TCanvas(BCN,BCN,1);
        cSr_Cauchy->cd(0);
        cSr_Cauchy->SetCanvasSize(1280, 720);
        cSr_Cauchy->SetMargin(0.15,0.05,0.2,0.05);//lrbt
        hSrAxis->Draw("axis");
        gSrCore->Draw("C,same");
        gSrFull->Draw("C,same");
        gSrCauchy->Draw("C,same");
        l_Cauchy->Draw("same");
        if(iLog%2) cSr_Cauchy->SetLogy();
        cSr_Cauchy->Write();
        cSr_Cauchy->SaveAs(FolderName+BCN+".pdf");
        cSr_Cauchy->SaveAs(FolderName+BCN+".png");

        delete hSrAxis;
        delete cSr_Core;
        delete cSr_Short;
        delete cSr_ShortRho;
        delete cSr_NoLong;
        delete cSr_Omega;
        delete cSr_OmegaLambda;
        delete cSr_Full;
        delete cSr_LongLambda;
        delete cSr_Cauchy;
    }



}

void pion_core_effect(const bool& randomangle){
    const double SourceSize = 1.5;
    const unsigned NumMomBins = 100;
    const double kMin = 0;
    const double kMax = 200;
    const unsigned NumRadBins = 256;//128
    const double rMin = 0;
    const double rMax = 15;

//0.270 (0.270) prim
//0.055 (0.325) super long
//0.075 (0.400) omega
//0.175 (0.575) rho
//0.175 (0.750) 1 - 2
//0.100 (0.850) 2 - 5
//0.150 (1.000) <1
    const double prim = 0.270;
    const double reso_short = 0.150;
    const double reso_rho = 0.175;
    const double reso_inter = 0.175;
    const double reso_long = 0.100;
    const double reso_omega = 0.075;
    const double reso_infinite = 0.055;

    const double ctau_short = 0.5;
    const double ctau_rho = 1.3;
    const double ctau_ShortRho = (ctau_short*reso_short+ctau_rho*reso_rho)/(reso_short+reso_rho);
    const double ctau_long = 3.35;
    const double ctau_NoLong = (ctau_short*reso_short+ctau_rho*reso_rho+ctau_long*reso_long)/(reso_short+reso_rho+reso_long);

    TString OutputFolderName = TString::Format("%s/PionAnalysis/pion_core_effect/",GetFemtoOutputFolder());
    //TString OutputFolderName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/Pions/pion_core_effect/";
    //TString OutputFolderName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/Pions/JerkingOff/";
    TString NameAddOn = "";
    if(randomangle) NameAddOn = "_RA";
    TFile fInput(OutputFolderName+"fOutput"+NameAddOn+".root","recreate");
    fInput.Close();

    CATS KittyCore;
    DLM_CleverMcLevyResoTM CMLRTM_Core;
    KittyCore.SetMomBins(NumMomBins,kMin,kMax);
    SetUpReso_pipi_CUSTOM(KittyCore, CMLRTM_Core, 0, 1180, 0, 0, 0, randomangle);
    KittyCore.SetAnaSource(0,SourceSize);
    KittyCore.SetAutoNormSource(false);
    KittyCore.KillTheCat();
    DLM_Ck CatCore(1,0,KittyCore); RootFile_DlmCk(OutputFolderName+"fOutput"+NameAddOn+".root", "KittyCore", &CatCore);
    RootFile_DlmSource(OutputFolderName+"fOutput"+NameAddOn+".root","SourceCore",&KittyCore,NumRadBins,rMin,rMax);

    CATS KittyShort;
    DLM_CleverMcLevyResoTM CMLRTM_Short;
    KittyShort.SetMomBins(NumMomBins,kMin,kMax);
    SetUpReso_pipi_CUSTOM(KittyShort, CMLRTM_Short, ctau_short, 1300, reso_short, 0, 0, randomangle);
    printf("Short: cTau = %.2f; f = %.2f\n",ctau_short,reso_short);
    KittyShort.SetAnaSource(0,SourceSize);
    KittyShort.SetAutoNormSource(false);
    KittyShort.KillTheCat();
    DLM_Ck CatShort(1,0,KittyShort); RootFile_DlmCk(OutputFolderName+"fOutput"+NameAddOn+".root", "KittyShort", &CatShort);
    RootFile_DlmSource(OutputFolderName+"fOutput"+NameAddOn+".root","SourceShort",&KittyShort,NumRadBins,rMin,rMax);

    CATS KittyShortRho;
    DLM_CleverMcLevyResoTM CMLRTM_ShortRho;
    KittyShortRho.SetMomBins(NumMomBins,kMin,kMax);
    //SetUpReso_pipi(KittyShortRho, CMLRTM_ShortRho, ctau_ShortRho, 1000, reso_short+reso_rho, 0, 0, randomangle);
    SetUpReso_pipi_CUSTOM(KittyShortRho, CMLRTM_ShortRho, 1.0, 1000, 1.0, 0, 0, randomangle);
    printf("ShortRho: cTau = %.2f; f = %.2f\n",ctau_ShortRho,reso_short+reso_rho);
    KittyShortRho.SetAnaSource(0,SourceSize);
    KittyShortRho.SetAutoNormSource(false);
    KittyShortRho.KillTheCat();
    DLM_Ck CatShortRho(1,0,KittyShortRho); RootFile_DlmCk(OutputFolderName+"fOutput"+NameAddOn+".root", "KittyShortRho", &CatShortRho);
    RootFile_DlmSource(OutputFolderName+"fOutput"+NameAddOn+".root","SourceShortRho",&KittyShortRho,NumRadBins,rMin,rMax);

    CATS KittyNoLong;
    DLM_CleverMcLevyResoTM CMLRTM_NoLong;
    KittyNoLong.SetMomBins(NumMomBins,kMin,kMax);
    SetUpReso_pipi_CUSTOM(KittyNoLong, CMLRTM_NoLong, ctau_NoLong, 1180, reso_short+reso_rho+reso_inter+reso_long, 0, 0, randomangle);
    printf("NoLong: cTau = %.2f; f = %.2f\n",ctau_NoLong,reso_short+reso_rho+reso_inter+reso_long);
    KittyNoLong.SetAnaSource(0,SourceSize);
    KittyNoLong.SetAutoNormSource(false);
    KittyNoLong.KillTheCat();
    DLM_Ck CatNoLong(1,0,KittyNoLong); RootFile_DlmCk(OutputFolderName+"fOutput"+NameAddOn+".root", "KittyNoLong", &CatNoLong);
    RootFile_DlmSource(OutputFolderName+"fOutput"+NameAddOn+".root","SourceNoLong",&KittyNoLong,NumRadBins,rMin,rMax);

    CATS KittyOmega;
    DLM_CleverMcLevyResoTM CMLRTM_Omega;
    KittyOmega.SetMomBins(NumMomBins,kMin,kMax);
    SetUpReso_pipi_CUSTOM(KittyOmega, CMLRTM_Omega, 1.5, 1180, reso_short+reso_rho+reso_inter+reso_long, reso_omega, 0, randomangle);
    printf("Omega: cTau = %.2f; f = %.2f\n",1.5,reso_short+reso_rho+reso_inter+reso_long+reso_omega);
    KittyOmega.SetAnaSource(0,SourceSize);
    KittyOmega.SetAutoNormSource(false);
    KittyOmega.KillTheCat();
    DLM_Ck CatOmega(1,0,KittyOmega); RootFile_DlmCk(OutputFolderName+"fOutput"+NameAddOn+".root", "KittyOmega", &CatOmega);
    RootFile_DlmSource(OutputFolderName+"fOutput"+NameAddOn+".root","SourceOmega",&KittyOmega,NumRadBins,rMin,rMax);

    CATS KittyOmegaLambda;
    DLM_CleverMcLevyResoTM CMLRTM_OmegaLambda;
    KittyOmegaLambda.SetMomBins(NumMomBins,kMin,kMax);
    SetUpReso_pipi_CUSTOM(KittyOmegaLambda, CMLRTM_OmegaLambda, 1.5, 1180, reso_short+reso_rho+reso_inter+reso_long, 0, 0, randomangle);
    printf("OmegaLambda: cTau = %.2f; f = %.2f\n",1.5,reso_short+reso_rho+reso_inter+reso_long);
    KittyOmegaLambda.SetAnaSource(0,SourceSize);
    KittyOmegaLambda.SetAutoNormSource(false);
    KittyOmegaLambda.KillTheCat();
    DLM_Ck CatOmegaLambda(1,0,KittyOmegaLambda); RootFile_DlmCk(OutputFolderName+"fOutput"+NameAddOn+".root", "KittyOmegaLambda", &CatOmegaLambda);
    double lambda_par = (1.-reso_omega)*(1.-reso_omega);
    CatOmegaLambda.Scale(lambda_par);
    CatOmegaLambda.AddToAll(1.-lambda_par);
    RootFile_DlmCk(OutputFolderName+"fOutput"+NameAddOn+".root", "KittyOmegaLambda", &CatOmegaLambda);
    RootFile_DlmSource(OutputFolderName+"fOutput"+NameAddOn+".root","SourceOmegaLambda",&KittyOmegaLambda,NumRadBins,rMin,rMax,lambda_par);

    CATS KittyFull;
    DLM_CleverMcLevyResoTM CMLRTM_Full;
    KittyFull.SetMomBins(NumMomBins,kMin,kMax);
    SetUpReso_pipi_CUSTOM(KittyFull, CMLRTM_Full, 1.5, 1180, reso_short+reso_rho+reso_inter+reso_long, reso_omega, reso_infinite, randomangle);
    printf("Full: cTau = %.2f; f = %.2f\n",1.5,reso_short+reso_rho+reso_inter+reso_long+reso_omega+reso_infinite);
    KittyFull.SetAnaSource(0,SourceSize);
    KittyFull.SetAutoNormSource(false);
    KittyFull.KillTheCat();
    DLM_Ck CatFull(1,0,KittyFull); RootFile_DlmCk(OutputFolderName+"fOutput"+NameAddOn+".root", "KittyFull", &CatFull);
    RootFile_DlmSource(OutputFolderName+"fOutput"+NameAddOn+".root","SourceFull",&KittyFull,NumRadBins,rMin,rMax);

    CATS KittyLongLambda;
    DLM_CleverMcLevyResoTM CMLRTM_LongLambda;
    KittyLongLambda.SetMomBins(NumMomBins,kMin,kMax);
    SetUpReso_pipi_CUSTOM(KittyLongLambda, CMLRTM_LongLambda, 1.5, 1180, reso_short+reso_rho+reso_inter+reso_long, reso_omega, 0, randomangle);
    printf("LongLambda: cTau = %.2f; f = %.2f\n",1.5,reso_short+reso_rho+reso_inter+reso_long);
    KittyLongLambda.SetAnaSource(0,SourceSize);
    KittyLongLambda.SetAutoNormSource(false);
    KittyLongLambda.KillTheCat();
    DLM_Ck CatLongLambda(1,0,KittyLongLambda);
    lambda_par = (1.-reso_infinite)*(1.-reso_infinite);
    CatLongLambda.Scale(lambda_par);
    CatLongLambda.AddToAll(1.-lambda_par);
    RootFile_DlmCk(OutputFolderName+"fOutput"+NameAddOn+".root", "KittyLongLambda", &CatLongLambda);
    RootFile_DlmSource(OutputFolderName+"fOutput"+NameAddOn+".root","SourceLongLambda",&KittyLongLambda,NumRadBins,rMin,rMax,lambda_par);

    CATS KittyCauchy;
    KittyCauchy.SetMomBins(NumMomBins,kMin,kMax);
    DLM_CommonAnaFunctions AnaObject;
    AnaObject.SetUpCats_pipi(KittyCauchy,"Cauchy",0);
    KittyCauchy.SetAnaSource(0,2.0);
    KittyCauchy.SetAutoNormSource(false);
    KittyCauchy.KillTheCat();
    DLM_Ck CatCauchy(1,0,KittyCauchy);
    lambda_par = 0.9;
    CatCauchy.Scale(lambda_par);
    CatCauchy.AddToAll(1.-lambda_par);
    RootFile_DlmCk(OutputFolderName+"fOutput"+NameAddOn+".root", "KittyCauchy", &CatCauchy);
    RootFile_DlmSource(OutputFolderName+"fOutput"+NameAddOn+".root","SourceCauchy",&KittyCauchy,NumRadBins,rMin,rMax,lambda_par);

    pion_core_effect_plots(OutputFolderName,NameAddOn);

}


//0.26666666666666666, 1.7637104743083003
//0.41449275362318844, 1.8945775691699602
//0.6, 1.890367671277997
//0.7710144927536232, 1.5843626482213435
//1.0927536231884059, 1.0363451086956519

void FAST_MT_PLOTS(){
    double xAxis [] = {0.27,0.41,0.60,0.77,1.09};
    TString InputFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/Pions/FAST_MT_PLOTS/FitsFrom0/FitsFrom0.root";
    TString GraphName[3];
    GraphName[0] = "FitRadiusValues_1";
    GraphName[1] = "FitRadiusValues_2";
    GraphName[2] = "FitRadiusValues_3";

    TGraph gmt[3];
    gmt[0].SetName("gmt0");
    gmt[1].SetName("gmt1");
    gmt[2].SetName("gmt2");

    TFile InputFile(InputFileName,"read");

    for(unsigned uMult=0; uMult<3; uMult++){
        TGraph* graph = (TGraph*)InputFile.Get(GraphName[uMult]);
        for(unsigned uKt=0; uKt<5; uKt++){
            double yVal,xVal;
            graph->GetPoint(uKt,xVal,yVal);
            if(uMult!=0||uKt<=1) gmt[uMult].SetPoint(uKt,xAxis[uKt],yVal);
        }
    }

    TH1F* hAxis = new TH1F("hAxis","hAxis",128,0.2,1.25);
    hAxis->SetStats(false);
    hAxis->SetTitle("");
    hAxis->GetXaxis()->SetTitle("m_{T} (GeV)");
    hAxis->GetXaxis()->SetTitleSize(0.06);
    hAxis->GetXaxis()->SetLabelSize(0.06);
    hAxis->GetXaxis()->SetTitleOffset(1.3);
    hAxis->GetXaxis()->SetLabelOffset(0.02);

    hAxis->GetYaxis()->SetTitle("r_{core} (fm)");
    hAxis->GetYaxis()->SetTitleSize(0.06);
    hAxis->GetYaxis()->SetLabelSize(0.06);
    hAxis->GetYaxis()->SetTitleOffset(1.10);
    hAxis->GetYaxis()->SetRangeUser(0.8,2.4);


    gmt[0].SetMarkerColor(kGreen+1);
    gmt[1].SetMarkerColor(kBlue+1);
    gmt[2].SetMarkerColor(kRed+1);

    TCanvas** canvas = new TCanvas* [3];
    for(unsigned uMult=0; uMult<3; uMult++){

        gmt[uMult].SetLineWidth(0);
        gmt[uMult].SetLineColor(kWhite);
        gmt[uMult].SetMarkerSize(3);
        gmt[uMult].SetMarkerStyle(21);

        canvas[uMult] = new TCanvas(TString::Format("canvas%u",uMult), TString::Format("canvas%u",uMult), 1);
        canvas[uMult]->cd(0);
        canvas[uMult]->SetCanvasSize(1280, 720);
        canvas[uMult]->SetMargin(0.15,0.05,0.2,0.05);//lrbt
        hAxis->Draw("axis");
        gmt[uMult].Draw("P,same");

        canvas[uMult]->SaveAs(TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/Pions/FAST_MT_PLOTS/FitsFrom0/canvas%u.pdf",uMult));
    }




}

//evaluate the r_core based on mt
//create RSM prediction for r_eff
//create CECA r_eff based on some fixed disp or Jaime pars.

void pion_proton_FirstLook(){
  const double avg_mt = 1.1;
  const double rsm_core = 0;
  const double hadr_t = 3.55;
  const double time_d = 2.30;
  const double displ = 0.25;
  const double ceca_sp = 1.0;

  const double tau_p = 1.65;
  const double resoM_p = 1362;
  const double resoF_p = 64.22;

  const double tau_pi = 1.50;
  const double resoM_pi = 1180;
  const double resoF_pi = 68.2;

  //PP - easy
  //R - non-comon reso;
  TREPNI Database(0);
  Database.SetSeed(11);
  std::vector<TreParticle*> ParticleList;
  ParticleList.push_back(Database.NewParticle("Proton"));
  ParticleList.push_back(Database.NewParticle("Pion"));
  ParticleList.push_back(Database.NewParticle("Daughter"));//decay daughters (pions)
  ParticleList.push_back(Database.NewParticle("KaonReso"));
  ParticleList.push_back(Database.NewParticle("PionReso"));

}

//To get a fast first look: make few correlations using a Gaussian source
//we assume 100% purity
//double SourceSize, double LambdaGen, double Lambda_P1, double Lambda_P2, double Lambda_P3
//AVG_MT is 0.96
void pion_proton_FastLook(){
  //this seems okay, and perhaps the effective source will have long enough of a tail to put an extra lambda par
  double SourceSize=1.2;
  double SourceNorm=0.9;
  //the 0.876 (0.824) is for the pions, and I think it includes the long lived reso
  //this way L_gen is around 66.7%
  double LambdaGen=0.823*0.824*0.994*0.99*SourceNorm; double Lambda_P1=0.004; double Lambda_P2=0.002; double Lambda_P3=0.001;
  //N.B. we rely on the 4 MeV bin width!!!
  const unsigned NumMomBins = 200;
  const double kMin = 0;
  const double kMax = 800;
  //const double Weight0 = 0.5;
  //const double LambdaGen = 0.6;
  const double LambdaFlat = 1. - LambdaGen;

  //1.06012 / 1.006 = 1.54
  //-0.000143092
  //const double Pol0 = 1.06012 / 1.006;
  //const double Pol1 = -0.000143092 / 1.006;
  const double Pol0 = 1;
  const double Pol1 = 0;
  const double Pol2 = 0;
  const double Pol3 = 0;
  //const double Pol0 = 1.006;
  //const double Pol1 = 0;


  //const double Pol0 = 1.034;
  //const double Pol1 = 0;
  //const double Pol2 = 3.21e-7;
  //const double Pol3 = 3.92e-10;



  CATS Kitty;
  Kitty.SetMomBins(NumMomBins,kMin,kMax);
  Kitty.SetThetaDependentSource(false);
  CATSparameters cPars(CATSparameters::tSource,1,true);
  cPars.SetParameter(0,SourceSize);
  Kitty.SetAnaSource(GaussSource, cPars);
  Kitty.SetAnaSource(0,SourceSize);
  Kitty.SetUseAnalyticSource(true);
  Kitty.SetAutoNormSource(false);
  Kitty.SetMomentumDependentSource(false);

  Kitty.SetExcludeFailedBins(false);
  //same charged
  Kitty.SetQ1Q2(1);
  Kitty.SetQuantumStatistics(false);
  Kitty.SetRedMass( (Mass_p*Mass_pic)/(Mass_p+Mass_pic) );

  Kitty.SetNumChannels(1);
  Kitty.SetNumPW(0,1);
  //Kitty.SetNumPW(1,1);
  Kitty.SetSpin(0,0);
  //Kitty.SetSpin(1,1);
  Kitty.SetChannelWeight(0, 1.0);
  //Kitty.SetChannelWeight(1, 1.-Weight0);

  Kitty.SetEpsilonConv(5e-9);
  Kitty.SetEpsilonProp(5e-9);
  Kitty.SetMaxRad(96);
  Kitty.SetMaxRho(32);


  DLM_Ck Ck_P1(NumMomBins,kMin,kMax);
  DLM_Ck Ck_P2(NumMomBins,kMin,kMax);
  DLM_Ck Ck_P3(NumMomBins,kMin,kMax);

  TFile fIn_ppim(TString::Format("%s/Proton_Pion/Resonances/Proton_AntiPion/proton_pim_GeV.root",GetCernBoxDimi()),"read");
  TH1D* h_ppim_GeV = (TH1D*)fIn_ppim.Get("CF_Total_ProtonAntiPion_Norm");
  h_ppim_GeV->GetXaxis()->SetLimits(h_ppim_GeV->GetXaxis()->GetBinLowEdge(1)*1000,h_ppim_GeV->GetXaxis()->GetBinUpEdge(h_ppim_GeV->GetNbinsX())*1000);

  TFile fIn_ppip(TString::Format("%s/Proton_Pion/Resonances/Proton_AntiPion/proton_pip_GeV.root",GetCernBoxDimi()),"read");
  TH1D* h_ppip_GeV = (TH1D*)fIn_ppip.Get("CF_Total_ProtonPion_Norm");
  h_ppip_GeV->GetXaxis()->SetLimits(h_ppip_GeV->GetXaxis()->GetBinLowEdge(1)*1000,h_ppip_GeV->GetXaxis()->GetBinUpEdge(h_ppip_GeV->GetNbinsX())*1000);

  TH1D* h_ppi_ratio = (TH1D*)h_ppim_GeV->Clone("h_ppi_ratio");
  h_ppi_ratio->Divide(h_ppip_GeV);

  //N.B. we rely on the 4 MeV bin width!!!
  TFile fInput(TString::Format("%s/Proton_Pion/Resonances/Proton_AntiPion/FISTReso_Output_BenjaminBW_2.root",GetCernBoxDimi()),"read");
  TH1D* h_P1 = (TH1D*)fInput.Get("CF_FirstPeak");
  TH1D* h_P2 = (TH1D*)fInput.Get("CF_SecondPeak");
  TH1D* h_P3 = (TH1D*)fInput.Get("CF_ThirdPeak");
  for(unsigned uMom=0; uMom<NumMomBins; uMom++){
if(uMom>=12*0)
    //Ck_P1.SetBinContent(uMom,h_P1->GetBinContent(uMom+1));
    Ck_P1.SetBinContent(uMom-12*0,h_P1->GetBinContent(uMom+1));
    //Ck_P1.SetBinContent(uMom,1);
    Ck_P2.SetBinContent(uMom,h_P2->GetBinContent(uMom+1));
    //Ck_P2.SetBinContent(uMom,1);
    Ck_P3.SetBinContent(uMom,h_P3->GetBinContent(uMom+1));
    //printf("uB%u: %f\n",uMom,h_P1->GetBinContent(uMom+1));
    //usleep(250e3);
  }
  fInput.Close();

  Ck_P1.Update();
  Ck_P2.Update();
  Ck_P3.Update();

  DLM_CkDecomposition CkDec_P1("P1",0,Ck_P1,NULL);
  DLM_CkDecomposition CkDec_P2("P2",0,Ck_P2,NULL);
  DLM_CkDecomposition CkDec_P3("P3",0,Ck_P3,NULL);
  CkDec_P1.Update();
  CkDec_P2.Update();
  CkDec_P3.Update();

//SET UP THE ABOVE DLM_Ck

  DLM_Ck Ck_ppi(Kitty.GetNumSourcePars(),0,Kitty,NumMomBins,kMin,kMax);
  Ck_ppi.Update();
  //the contributions (apart from genuine):
  //0 : FlatFeed

  TFile fOutput(TString::Format("%s/Proton_Pion/Resonances/pion_proton_FastLook_SS=%.2f_GEN-%.0f_P1-%.0f_P2-%.0f_P3-%.0f_PerMil.root",
                        GetCernBoxDimi(),SourceSize,LambdaGen*1000.,Lambda_P1*1000.,Lambda_P2*1000.,Lambda_P3*1000.),"recreate");
  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles/",GetCernBoxDimi()).Data());
  TH2F* hResolution_pp = AnalysisObject.GetResolutionMatrix("pp13TeV_HM_DimiJun20","pp");
//printf("1\n");
  //DLM_Histo<float>* dlm_RESO = Convert_TH2F_DlmHisto(hResolution_pp);

  fOutput.cd();
  h_ppim_GeV->Write();
  h_ppip_GeV->Write();
  h_ppi_ratio->Write();
  double reso_kminX = hResolution_pp->GetXaxis()->GetBinLowEdge(1);
  double reso_kmaxX = hResolution_pp->GetXaxis()->GetBinUpEdge(hResolution_pp->GetNbinsX());
  double reso_kminY = hResolution_pp->GetYaxis()->GetBinLowEdge(1);
  double reso_kmaxY = hResolution_pp->GetYaxis()->GetBinUpEdge(hResolution_pp->GetNbinsY());
  TH2F* h_RESO = new TH2F("h_RESO","h_RESO",hResolution_pp->GetNbinsX(),reso_kminX,reso_kmaxX,hResolution_pp->GetNbinsY(),reso_kminY,reso_kmaxY);
  for(unsigned uBinX=0; uBinX<hResolution_pp->GetNbinsX(); uBinX++){
    for(unsigned uBinY=0; uBinY<hResolution_pp->GetNbinsY(); uBinY++){
      h_RESO->SetBinContent(uBinX+1,uBinY+1,hResolution_pp->GetBinContent(uBinX+1,uBinY+1));
    }
  }

  delete hResolution_pp; hResolution_pp=NULL;
//printf("2\n");
  //fOutput.cd();
  //h_RESO->Write();
  //return;

  //1,2,3 : The 3 peaks from resonances (excluding Lambda decay!)
//!!!!! DLM_CkDecomposition CkDec_ppi("ppi",1,Ck_ppi,hResolution_ppi);
  DLM_CkDecomposition CkDec_ppi("ppi",4,Ck_ppi,h_RESO);
//0.02,0.004,0.002
  CkDec_ppi.QA_LambdaPar(false);
  CkDec_ppi.AddContribution(0,LambdaFlat,DLM_CkDecomposition::cFeedDown);
  printf("1\n");
  CkDec_ppi.AddContribution(1,0.004,DLM_CkDecomposition::cFeedDown,&CkDec_P1);
  //CkDec_ppi.AddContribution(1,Lambda_P1,DLM_CkDecomposition::cFeedDown,&CkDec_P1);
  printf("2\n");
  CkDec_ppi.AddContribution(2,0.002,DLM_CkDecomposition::cFeedDown,&CkDec_P2);
  //CkDec_ppi.AddContribution(2,Lambda_P2,DLM_CkDecomposition::cFeedDown,&CkDec_P2);
  printf("3\n");
  CkDec_ppi.AddContribution(3,0.001,DLM_CkDecomposition::cFeedDown,&CkDec_P3);
  //CkDec_ppi.AddContribution(3,Lambda_P3,DLM_CkDecomposition::cFeedDown,&CkDec_P3);
  printf("4\n");

  CkDec_ppi.SetLambdaMain(LambdaGen);
//return;
  //same charged, coulomb only
  Kitty.KillTheCat();
  Ck_ppi.Update(true);
  CkDec_ppi.Update(true,true);

  TGraph gCk_ppi_SameCh_Coul;
  gCk_ppi_SameCh_Coul.SetName("gCk_ppi_SameCh_Coul");
  gCk_ppi_SameCh_Coul.SetMarkerColor(kRed);
  gCk_ppi_SameCh_Coul.SetMarkerStyle(3);
  gCk_ppi_SameCh_Coul.SetMarkerSize(1);
  gCk_ppi_SameCh_Coul.SetLineColor(kRed);
  gCk_ppi_SameCh_Coul.SetLineWidth(4);

  TGraph dCk_ppi_SameCh_Coul;
  dCk_ppi_SameCh_Coul.SetName("dCk_ppi_SameCh_Coul");
  dCk_ppi_SameCh_Coul.SetMarkerColor(kRed);
  dCk_ppi_SameCh_Coul.SetMarkerStyle(3);
  dCk_ppi_SameCh_Coul.SetMarkerSize(1);
  dCk_ppi_SameCh_Coul.SetLineColor(kRed);
  dCk_ppi_SameCh_Coul.SetLineWidth(4);

  TGraph ThCk_ppi_SameCh_Coul;
  ThCk_ppi_SameCh_Coul.SetName("ThCk_ppi_SameCh_Coul");
  ThCk_ppi_SameCh_Coul.SetMarkerColor(kRed);
  ThCk_ppi_SameCh_Coul.SetMarkerStyle(3);
  ThCk_ppi_SameCh_Coul.SetMarkerSize(1);
  ThCk_ppi_SameCh_Coul.SetLineColor(kRed);
  ThCk_ppi_SameCh_Coul.SetLineWidth(4);
  for(unsigned uMom=0; uMom<NumMomBins; uMom++){
      double MOM = Kitty.GetMomentum(uMom);
      double BL = (Pol0+Pol1*MOM+Pol2*MOM+Pol3*MOM);
      gCk_ppi_SameCh_Coul.SetPoint(uMom,MOM,BL*CkDec_ppi.EvalCk(MOM));
      ThCk_ppi_SameCh_Coul.SetPoint(uMom,MOM,Kitty.GetCorrFun(uMom));
      dCk_ppi_SameCh_Coul.SetPoint(uMom,MOM,h_ppip_GeV->GetBinContent(uMom+1)/BL/CkDec_ppi.EvalCk(MOM));
  }
  fOutput.cd();
  gCk_ppi_SameCh_Coul.Write();
  dCk_ppi_SameCh_Coul.Write();
  ThCk_ppi_SameCh_Coul.Write();

  //opposite charged, coulomb only
  Kitty.SetQ1Q2(-1);
  Kitty.KillTheCat();
  Ck_ppi.Update(true);
  CkDec_ppi.Update(true,true);

  TGraph gCk_ppi_OppCh_Coul;
  gCk_ppi_OppCh_Coul.SetName("gCk_ppi_OppCh_Coul");
  gCk_ppi_OppCh_Coul.SetMarkerColor(kBlue);
  gCk_ppi_OppCh_Coul.SetMarkerStyle(3);
  gCk_ppi_OppCh_Coul.SetMarkerSize(1);
  gCk_ppi_OppCh_Coul.SetLineColor(kBlue);
  gCk_ppi_OppCh_Coul.SetLineWidth(4);

  TGraph dCk_ppi_OppCh_Coul;
  dCk_ppi_OppCh_Coul.SetName("dCk_ppi_OppCh_Coul");
  dCk_ppi_OppCh_Coul.SetMarkerColor(kBlue);
  dCk_ppi_OppCh_Coul.SetMarkerStyle(3);
  dCk_ppi_OppCh_Coul.SetMarkerSize(1);
  dCk_ppi_OppCh_Coul.SetLineColor(kBlue);
  dCk_ppi_OppCh_Coul.SetLineWidth(4);


  TGraph ThCk_ppi_OppCh_Coul;
  ThCk_ppi_OppCh_Coul.SetName("ThCk_ppi_OppCh_Coul");
  ThCk_ppi_OppCh_Coul.SetMarkerColor(kBlue);
  ThCk_ppi_OppCh_Coul.SetMarkerStyle(3);
  ThCk_ppi_OppCh_Coul.SetMarkerSize(1);
  ThCk_ppi_OppCh_Coul.SetLineColor(kBlue);
  ThCk_ppi_OppCh_Coul.SetLineWidth(4);

  for(unsigned uMom=0; uMom<NumMomBins; uMom++){
    double MOM = Kitty.GetMomentum(uMom);
    double BL = (Pol0+Pol1*MOM+Pol2*MOM+Pol3*MOM / 1.006);
      gCk_ppi_OppCh_Coul.SetPoint(uMom,MOM,BL*CkDec_ppi.EvalCk(MOM));
      ThCk_ppi_OppCh_Coul.SetPoint(uMom,MOM,Kitty.GetCorrFun(uMom));
      dCk_ppi_OppCh_Coul.SetPoint(uMom,MOM,h_ppim_GeV->GetBinContent(uMom+1)/BL/CkDec_ppi.EvalCk(MOM));
  }
  fOutput.cd();
  gCk_ppi_OppCh_Coul.Write();
  dCk_ppi_OppCh_Coul.Write();
  ThCk_ppi_OppCh_Coul.Write();

  CATSparameters pPars0(CATSparameters::tPotential,4,true);
  //pPars0.SetParameter(0,1.973303e+02);
  //pPars0.SetParameter(1,6.187704e-01);
  //pPars0.SetParameter(2,-1.100288e+03);
  //pPars0.SetParameter(3,4.147391e-01);
  Kitty.SetShortRangePotential(0,0,DoubleGaussSum,pPars0);
  Kitty.SetShortRangePotential(0,0,0,-3.227447e+01);
  Kitty.SetShortRangePotential(0,0,1,1.077665e+00);
  Kitty.SetShortRangePotential(0,0,2,-2.228376e+02);
  Kitty.SetShortRangePotential(0,0,3,9.697892e-02 );

  //CATSparameters pPars1(CATSparameters::tPotential,4,true);
  //pPars1.SetParameter(0,7.495875e+02);
  //pPars1.SetParameter(1,3.944096e-01);
  //pPars1.SetParameter(2,6.549811e+02);
  //pPars1.SetParameter(3,4.493379e-01);

  //opposite charged, SI + Coulomb
  Kitty.KillTheCat();
  Ck_ppi.Update(true);
  CkDec_ppi.Update(true,true);

  TGraph gCk_ppi_OppCh_Full;
  gCk_ppi_OppCh_Full.SetName("gCk_ppi_OppCh_Full");
  gCk_ppi_OppCh_Full.SetMarkerColor(kBlue+2);
  gCk_ppi_OppCh_Full.SetMarkerStyle(3);
  gCk_ppi_OppCh_Full.SetMarkerSize(1);
  gCk_ppi_OppCh_Full.SetLineColor(kBlue+2);
  gCk_ppi_OppCh_Full.SetLineWidth(4);

  TGraph dCk_ppi_OppCh_Full;
  dCk_ppi_OppCh_Full.SetName("dCk_ppi_OppCh_Full");
  dCk_ppi_OppCh_Full.SetMarkerColor(kBlue+2);
  dCk_ppi_OppCh_Full.SetMarkerStyle(3);
  dCk_ppi_OppCh_Full.SetMarkerSize(1);
  dCk_ppi_OppCh_Full.SetLineColor(kBlue+2);
  dCk_ppi_OppCh_Full.SetLineWidth(4);

  TGraph ThCk_ppi_OppCh_Full;
  ThCk_ppi_OppCh_Full.SetName("ThCk_ppi_OppCh_Full");
  ThCk_ppi_OppCh_Full.SetMarkerColor(kBlue+2);
  ThCk_ppi_OppCh_Full.SetMarkerStyle(3);
  ThCk_ppi_OppCh_Full.SetMarkerSize(1);
  ThCk_ppi_OppCh_Full.SetLineColor(kBlue+2);
  ThCk_ppi_OppCh_Full.SetLineWidth(4);

  for(unsigned uMom=0; uMom<NumMomBins; uMom++){
    double MOM = Kitty.GetMomentum(uMom);
    double BL = (Pol0+Pol1*MOM+Pol2*MOM+Pol3*MOM / 1.006);
      gCk_ppi_OppCh_Full.SetPoint(uMom,MOM,BL*CkDec_ppi.EvalCk(MOM));
      ThCk_ppi_OppCh_Full.SetPoint(uMom,MOM,Kitty.GetCorrFun(uMom));
      dCk_ppi_OppCh_Full.SetPoint(uMom,MOM,h_ppim_GeV->GetBinContent(uMom+1)/BL/CkDec_ppi.EvalCk(MOM));
  }
  fOutput.cd();
  gCk_ppi_OppCh_Full.Write();
  dCk_ppi_OppCh_Full.Write();
  ThCk_ppi_OppCh_Full.Write();


  Kitty.SetQ1Q2(1);
  Kitty.SetShortRangePotential(0,0,0,-1.298001e+01);
  Kitty.SetShortRangePotential(0,0,1,1.453815e+00);
  Kitty.SetShortRangePotential(0,0,2,-2.392445e+01);
  Kitty.SetShortRangePotential(0,0,3,1.114436e-01);
  //same charged, SI + Coulomb
  Kitty.KillTheCat();
  Ck_ppi.Update(true);
  CkDec_ppi.Update(true,true);

  TGraph gCk_ppi_SameCh_Full;
  gCk_ppi_SameCh_Full.SetName("gCk_ppi_SameCh_Full");
  gCk_ppi_SameCh_Full.SetMarkerColor(kRed+2);
  gCk_ppi_SameCh_Full.SetMarkerStyle(3);
  gCk_ppi_SameCh_Full.SetMarkerSize(1);
  gCk_ppi_SameCh_Full.SetLineColor(kRed+2);
  gCk_ppi_SameCh_Full.SetLineWidth(4);

  TGraph dCk_ppi_SameCh_Full;
  dCk_ppi_SameCh_Full.SetName("dCk_ppi_SameCh_Full");
  dCk_ppi_SameCh_Full.SetMarkerColor(kRed+2);
  dCk_ppi_SameCh_Full.SetMarkerStyle(3);
  dCk_ppi_SameCh_Full.SetMarkerSize(1);
  dCk_ppi_SameCh_Full.SetLineColor(kRed+2);
  dCk_ppi_SameCh_Full.SetLineWidth(4);

  TGraph ThCk_ppi_SameCh_Full;
  ThCk_ppi_SameCh_Full.SetName("ThCk_ppi_SameCh_Full");
  ThCk_ppi_SameCh_Full.SetMarkerColor(kRed+2);
  ThCk_ppi_SameCh_Full.SetMarkerStyle(3);
  ThCk_ppi_SameCh_Full.SetMarkerSize(1);
  ThCk_ppi_SameCh_Full.SetLineColor(kRed+2);
  ThCk_ppi_SameCh_Full.SetLineWidth(4);
  for(unsigned uMom=0; uMom<NumMomBins; uMom++){
    double MOM = Kitty.GetMomentum(uMom);
    double BL = (Pol0+Pol1*MOM+Pol2*MOM+Pol3*MOM / 1.006);
      gCk_ppi_SameCh_Full.SetPoint(uMom,MOM,BL*CkDec_ppi.EvalCk(MOM));
      ThCk_ppi_SameCh_Full.SetPoint(uMom,MOM,Kitty.GetCorrFun(uMom));
      dCk_ppi_SameCh_Full.SetPoint(uMom,MOM,h_ppip_GeV->GetBinContent(uMom+1)/BL/CkDec_ppi.EvalCk(MOM));
  }
  fOutput.cd();
  gCk_ppi_SameCh_Full.Write();
  dCk_ppi_SameCh_Full.Write();
  ThCk_ppi_SameCh_Full.Write();
  h_RESO->Write();

  delete h_RESO;
  delete h_ppi_ratio;
}


int PION_ANA(int narg, char** ARGS){
    //FAST_MT_PLOTS();
    //pion_core_effect(false);
    pion_core_effect(true);
    //pipi_alarm();
    //pipi_test();
    //pipi_omega();
    //ProveWeDoNotNeedOmega();
    //EffectOnTailFor_pp_pL();

    //pion_proton_FirstLook();
    //pion_proton_FastLook(1.0,0.6,0.02,0.004,0.002);
    //pion_proton_FastLook();
    //pion_proton_FastLook(1.0,0.6,0.1,0.1,0.1);
    //pion_proton_FastLook(1.2,0.6,0.1,0.1,0.1);

    return 0;
}
