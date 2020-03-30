#include "PionAnalysis.h"

#include "CATS.h"
#include "DLM_Random.h"
#include "DLM_Potentials.h"
#include "DLM_Source.h"

#include "CommonAnaFunctions.h"

#include "TString.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TPaveText.h"
#include "TNtuple.h"

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

    TFile fOutput("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pipi_test/fOutput2_kickedOmega.root","recreate");
    g_pipi.Write();
    gs_pipi.Write();

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

int PION_ANA(int narg, char** ARGS){
    pipi_test();
    //pipi_omega();
    //ProveWeDoNotNeedOmega();
    //EffectOnTailFor_pp_pL();
}
