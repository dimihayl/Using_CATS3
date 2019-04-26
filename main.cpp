
#include "CommonAnaFunctions.h"
#include "GentleDimi.h"
#include "FemtoBoyzScripts.h"


#include <iostream>
#include <stdio.h>
#include <string.h>
#include <omp.h>
#include <complex>

#include "CATS.h"
#include "CATStools.h"
#include "DLM_Potentials.h"
#include "DLM_Source.h"
#include "DLM_CppTools.h"
#include "DLM_CkDecomposition.h"
#include "DLM_Random.h"
#include "DLM_Bessel.h"
#include "DLM_Integration.h"
#include "DLM_WfModel.h"
#include "DLM_Histo.h"
#include "pLambda_1.h"
#include "ForBernie.h"
#include "ForGerhard.h"

#include "TGraph.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TNtuple.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TNtuple.h"


using namespace std;

/*
void TestCATS3_NewExtWf(TString POT_TYPE){
    //const double kMin=0;
    //const double kMax = 320;
    unsigned NumMomBins;
    double* MomentumBins=NULL;
    double* FitRegion=NULL;

    DLM_CommonAnaFunctions AnalysisObject;

    AnalysisObject.SetUpBinning_pL("pp13TeV_MB_Run2paper",NumMomBins,MomentumBins,FitRegion);

    printf("NumMomBins=%u\n",NumMomBins);
    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        printf(" Bin %u: %.2f -> %.2f\n", uBin, MomentumBins[uBin], MomentumBins[uBin+1]);
    }

    CATS Kitty;
    Kitty.SetMomBins(NumMomBins,MomentumBins);
    AnalysisObject.SetUpCats_pL(Kitty,POT_TYPE,"Gauss");
    //printf("KILL\n");
    Kitty.KillTheCat();

    TFile* OutputFile = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/TestCATS3_NewExtWf/"+POT_TYPE+".root","recreate");
    TGraph gCk;
    gCk.SetName(TString::Format("gCk_%s",POT_TYPE.Data()));

    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        gCk.SetPoint(uBin,Kitty.GetMomentum(uBin),Kitty.GetCorrFun(uBin));

        TGraph gS0,gS1;
        gS0.SetName(TString::Format("gOldS0_%s_%.2f",POT_TYPE.Data(),Kitty.GetMomentum(uBin)));
        gS0.Set(100);
        gS1.SetName(TString::Format("gOldS1_%s_%.2f",POT_TYPE.Data(),Kitty.GetMomentum(uBin)));
        gS1.Set(100);

        TGraph gS0wf2,gS1wf2;
        gS0wf2.SetName(TString::Format("gOldS0wf2_%.2f",Kitty.GetMomentum(uBin)));
        gS0wf2.Set(100);
        gS1wf2.SetName(TString::Format("gOldS1wf2_%.2f",Kitty.GetMomentum(uBin)));
        gS1wf2.Set(100);

        //unsigned uRad=0;
        unsigned NumRadBins = Kitty.SavedWaveFunBins[uBin][0][0];
        for(unsigned uRad=0; uRad<NumRadBins; uRad++){
            double RADIUS = 0.5*(Kitty.WaveFunRad[uBin][0][0][uRad]+Kitty.WaveFunRad[uBin][0][0][uRad+1])*Kitty.NuFm();
            gS0.SetPoint(uRad,Kitty.WaveFunRad[uBin][0][0][uRad],abs(Kitty.WaveFunctionU[uBin][0][0][uRad]));
//printf(" %f vs %f\n",abs(Kitty.WaveFunctionU[uBin][0][0][uRad]),abs(Kitty.EvalRadialWaveFunction(uBin,0,0,RADIUS,false)));
            gS0wf2.SetPoint(uRad,Kitty.WaveFunRad[uBin][0][0][uRad],Kitty.EvalWaveFun2(uBin,Kitty.WaveFunRad[uBin][0][0][uRad]*Kitty.NuFm(),0));
        }
        NumRadBins = Kitty.SavedWaveFunBins[uBin][1][0];
        for(unsigned uRad=0; uRad<NumRadBins; uRad++){
            gS1.SetPoint(uRad,Kitty.WaveFunRad[uBin][1][0][uRad],abs(Kitty.WaveFunctionU[uBin][1][0][uRad]));
            gS1wf2.SetPoint(uRad,Kitty.WaveFunRad[uBin][1][0][uRad],Kitty.EvalWaveFun2(uBin,Kitty.WaveFunRad[uBin][1][0][uRad]*Kitty.NuFm(),1));
        }

        gS0.Write();
        gS1.Write();
        gS0wf2.Write();
        gS1wf2.Write();
    }

    gCk.Write();
    delete OutputFile;
}
*/
/*
void TestCommonInit(){
    CATS NormalCat;

    CATSparameters* cPars = NULL;
    DLM_CleverLevy* CleverLevy = NULL;
    CATSparameters* cPotPars1S0 = NULL;
    CATSparameters* cPotPars3P0 = NULL;
    CATSparameters* cPotPars3P1 = NULL;
    CATSparameters* cPotPars3P2 = NULL;

    cPars = new CATSparameters(CATSparameters::tSource,1,true);
    cPars->SetParameter(0,1.2);

    NormalCat.SetAnaSource(GaussSource, *cPars);
    NormalCat.SetUseAnalyticSource(true);

    //#,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
    double PotPars1S0[8]={NN_AV18,v18_Coupled3P2,1,1,1,0,0,0};
    double PotPars3P0[8]={NN_AV18,v18_Coupled3P2,1,1,1,1,1,0};
    double PotPars3P1[8]={NN_AV18,v18_Coupled3P2,1,1,1,1,1,1};
    double PotPars3P2[8]={NN_AV18,v18_Coupled3P2,1,1,1,1,1,2};
    cPotPars1S0 = new CATSparameters(CATSparameters::tPotential,8,true); cPotPars1S0->SetParameters(PotPars1S0);
    cPotPars3P0 = new CATSparameters(CATSparameters::tPotential,8,true); cPotPars3P0->SetParameters(PotPars3P0);
    cPotPars3P1 = new CATSparameters(CATSparameters::tPotential,8,true); cPotPars3P1->SetParameters(PotPars3P1);
    cPotPars3P2 = new CATSparameters(CATSparameters::tPotential,8,true); cPotPars3P2->SetParameters(PotPars3P2);
    NormalCat.SetMomBins(94,0,376);
    NormalCat.SetMomentumDependentSource(false);
    NormalCat.SetThetaDependentSource(false);
    NormalCat.SetExcludeFailedBins(false);

    NormalCat.SetQ1Q2(1);
    NormalCat.SetPdgId(2212, 2212);
    NormalCat.SetRedMass( 0.5*Mass_p );

    NormalCat.SetNumChannels(4);
    NormalCat.SetNumPW(0,2);
    NormalCat.SetNumPW(1,2);
    NormalCat.SetNumPW(2,2);
    NormalCat.SetNumPW(3,2);
    NormalCat.SetSpin(0,0);
    NormalCat.SetSpin(1,1);
    NormalCat.SetSpin(2,1);
    NormalCat.SetSpin(3,1);
    NormalCat.SetChannelWeight(0, 3./12.);
    NormalCat.SetChannelWeight(1, 1./12.);
    NormalCat.SetChannelWeight(2, 3./12.);
    NormalCat.SetChannelWeight(3, 5./12.);

    if(cPotPars1S0) NormalCat.SetShortRangePotential(0,0,fDlmPot,*cPotPars1S0);
    if(cPotPars3P0) NormalCat.SetShortRangePotential(1,1,fDlmPot,*cPotPars3P0);
    if(cPotPars3P1) NormalCat.SetShortRangePotential(2,1,fDlmPot,*cPotPars3P1);
    if(cPotPars3P2) NormalCat.SetShortRangePotential(3,1,fDlmPot,*cPotPars3P2);


    printf("NormalCat is set up\n");

    CATS CommonCat;
    CommonCat.SetMomBins(94,0,376);
    SetUpCats_pp(CommonCat,"AV18","Gauss");
    printf("CommonCat is set up\n");

    NormalCat.KillTheCat();
    printf("NormalCat was killed\n");

    CommonCat.KillTheCat();
    printf("CommonCat was killed\n");

    for(unsigned uBin=0; uBin<94; uBin++){
       printf("C(%.0f==%.0f) = %.2f (%.2f)\n",NormalCat.GetMomentum(uBin),CommonCat.GetMomentum(uBin),
              NormalCat.GetCorrFun(uBin),CommonCat.GetCorrFun(uBin));
    }


    delete cPars;
    delete cPotPars1S0;
    delete cPotPars3P0;
    delete cPotPars3P1;
    delete cPotPars3P2;

}
*/
void plot_pp(){

    const double SourceSize = 1.0;
    const unsigned NumStabilitySteps = 1;
    const double StabMin = 1.6;
    const double StabMax = 1.6;
    const double StabStep = NumStabilitySteps>1?(StabMax-StabMin)/double(NumStabilitySteps-1):0;

    const double kMin=0;
    const double kFine=200;
    const unsigned NumFineBins = 40;
    const double kFineStep = (kFine-kMin)/double(NumFineBins);
    const double kMax=500;
    const unsigned NumCoarseBins = 30;
    const double kCoarseStep = (kMax-kFine)/double(NumCoarseBins);
    const unsigned NumMomBins=NumFineBins+NumCoarseBins;
    double* MomBins = new double[NumMomBins];
    for(unsigned uBin=0; uBin<=NumFineBins; uBin++){
        MomBins[uBin] = kMin+double(uBin)*kFineStep;
    }
    for(unsigned uBin=NumFineBins+1; uBin<=NumMomBins; uBin++){
        MomBins[uBin] = MomBins[uBin-1]+kCoarseStep;;
    }

    TFile* OutputFile = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/plot_pp/OutputFile.root","recreate");

    for(unsigned uStab=0; uStab<NumStabilitySteps; uStab++){

        DLM_CommonAnaFunctions AnalysisObject;

        double Stability = StabMin+double(uStab)*StabStep;
        printf("\nStability = %.2f\n",Stability);
        CATS Kitty;
        Kitty.SetMaxNumThreads(1);
        Kitty.SetMomBins(NumMomBins,MomBins);
        AnalysisObject.SetUpCats_pp(Kitty,"AV18","CleverLevy_Nolan");
        Kitty.SetAnaSource(0,SourceSize);
        Kitty.SetAnaSource(1,Stability);
        Kitty.SetEpsilonProp(1e-8);
        Kitty.KillTheCat();
        TGraph gKitty;
        gKitty.SetName(TString::Format("gKitty_%.2f",Stability));
        gKitty.Set(NumMomBins);
        for(unsigned uBin=0; uBin<NumMomBins; uBin++){
            printf("C(%.2f) = %.2f\n",Kitty.GetMomentum(uBin),Kitty.GetCorrFun(uBin));
            gKitty.SetPoint(uBin,Kitty.GetMomentum(uBin),Kitty.GetCorrFun(uBin));
        }
        gKitty.Write();
    }

    delete [] MomBins;
    delete OutputFile;
}


void plot_pL(){

    const double SourceSize = 1.3;
    const unsigned NumStabilitySteps = 1;
    const double StabMin = 2.0;
    const double StabMax = 2.0;
    const double StabStep = NumStabilitySteps>1?(StabMax-StabMin)/double(NumStabilitySteps-1):0;

    const double kMin=0;
    const double kFine=200;
    const unsigned NumFineBins = 20;
    const double kFineStep = (kFine-kMin)/double(NumFineBins);
    const double kMax=340;
    const unsigned NumCoarseBins = 14;
    const double kCoarseStep = (kMax-kFine)/double(NumCoarseBins);
    const unsigned NumMomBins=NumFineBins+NumCoarseBins;
    double* MomBins = new double[NumMomBins];
    for(unsigned uBin=0; uBin<=NumFineBins; uBin++){
        MomBins[uBin] = kMin+double(uBin)*kFineStep;
    }
    for(unsigned uBin=NumFineBins+1; uBin<=NumMomBins; uBin++){
        MomBins[uBin] = MomBins[uBin-1]+kCoarseStep;;
    }

    TFile* OutputFile = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/plot_pL/OutputFile.root","recreate");

    for(unsigned uStab=0; uStab<NumStabilitySteps; uStab++){

        DLM_CommonAnaFunctions AnalysisObject;

        double Stability = StabMin+double(uStab)*StabStep;
        printf("\nStability = %.2f\n",Stability);
        CATS KittyNLO;
        KittyNLO.SetMaxNumThreads(1);
        KittyNLO.SetMomBins(NumMomBins,MomBins);
        AnalysisObject.SetUpCats_pL(KittyNLO,"NLO_Coupled_S","Gauss");
        KittyNLO.SetAnaSource(0,SourceSize);
        KittyNLO.SetAnaSource(1,Stability);
        KittyNLO.SetEpsilonProp(1e-8);
        KittyNLO.KillTheCat();
        TGraph gKittyNLO;
        gKittyNLO.SetName(TString::Format("gKittyNLO_%.2f",Stability));
        gKittyNLO.Set(NumMomBins);
        for(unsigned uBin=0; uBin<NumMomBins; uBin++){
            printf("C(%.2f) = %.3f\n",KittyNLO.GetMomentum(uBin),KittyNLO.GetCorrFun(uBin));
            gKittyNLO.SetPoint(uBin,KittyNLO.GetMomentum(uBin),KittyNLO.GetCorrFun(uBin));
        }
        gKittyNLO.Write();

        CATS KittyLO;
        KittyLO.SetMaxNumThreads(1);
        KittyLO.SetMomBins(NumMomBins,MomBins);
        AnalysisObject.SetUpCats_pL(KittyLO,"LO","Gauss");
        KittyLO.SetAnaSource(0,SourceSize);
        KittyLO.SetAnaSource(1,Stability);
        KittyLO.SetEpsilonProp(1e-8);
        KittyLO.KillTheCat();
        TGraph gKittyLO;
        gKittyLO.SetName(TString::Format("gKittyLO_%.2f",Stability));
        gKittyLO.Set(NumMomBins);
        for(unsigned uBin=0; uBin<NumMomBins; uBin++){
            printf("C(%.2f) = %.3f\n",KittyLO.GetMomentum(uBin),KittyLO.GetCorrFun(uBin));
            gKittyLO.SetPoint(uBin,KittyLO.GetMomentum(uBin),KittyLO.GetCorrFun(uBin));
        }
        gKittyLO.Write();
    }

    delete [] MomBins;
    delete OutputFile;
}

void plot_pSigma0(){

    const double SourceSize = 1.2;
    const double SourceStability = 1.5;

    const double kMin=0;
    const double kMax=320;
    const unsigned NumMomBins = 64;

    TFile* OutputFile = new TFile(
                TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/plot_pSigma0/OutputFile_%.2f_%.2f.root",
                                SourceStability,SourceSize),"recreate");
    CATS Kitty;
    Kitty.SetMomBins(NumMomBins,kMin,kMax);

    DLM_CleverLevy CleverLevy;
    CleverLevy.InitStability(20,1,2);
    CleverLevy.InitScale(35,0.25,2.0);
    CleverLevy.InitRad(256,0,64);
    CleverLevy.InitType(2);
    Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy, 2);
    Kitty.SetAnaSource(0,SourceSize);
    Kitty.SetAnaSource(1,SourceStability);
    Kitty.SetUseAnalyticSource(true);
    Kitty.SetMomentumDependentSource(false);
    Kitty.SetThetaDependentSource(false);
    Kitty.SetExcludeFailedBins(false);

    DLM_Histo<complex<double>>*** ExternalWF=NULL;
    ExternalWF = Init_pSigma0_Haidenbauer("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/Haidenbauer/pSigma0/",Kitty);
    for(unsigned uCh=0; uCh<Kitty.GetNumChannels(); uCh++){
        Kitty.SetExternalWaveFunction(uCh,0,ExternalWF[0][uCh][0],ExternalWF[1][uCh][0]);
    }

    Kitty.KillTheCat();
    TGraph gKitty;
    gKitty.SetName(TString::Format("gKitty"));
    gKitty.Set(NumMomBins);
    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        printf("C(%.2f) = %.2f\n",Kitty.GetMomentum(uBin),Kitty.GetCorrFun(uBin));
        gKitty.SetPoint(uBin,Kitty.GetMomentum(uBin),Kitty.GetCorrFun(uBin));
    }
    gKitty.Write();

    CleanUpWfHisto(Kitty,ExternalWF);

    delete OutputFile;
}

void McLevyTest1(){

    const double MaxRad=64;
    const unsigned NumRadPts=512;

    const double SourceScale = 1.0;
    const double SourceStability = 2.0;

    DLM_CleverMcLevyReso McSourceReso;
    McSourceReso.InitScale(1,SourceScale-0.01,SourceScale+0.01);
    McSourceReso.InitStability(1,SourceStability-0.01,SourceStability+0.01);
    McSourceReso.InitRad(NumRadPts,0,MaxRad);
    McSourceReso.InitType(2);
    McSourceReso.InitReso(0,1);
    McSourceReso.InitReso(1,1);

    McSourceReso.SetUpReso(0,0,1.-0.3578,1361.52,1.65,938.272,138);
    //McSourceReso.SetUpReso(0,0,1.-0.3562,1462.93,4.69,1115.683,138);
    McSourceReso.SetUpReso(1,0,1.-0.3578,1361.52,1.65,938.272,138);
    //McSourceReso.SetUpReso(1,0,1.-0.3562,1462.93,4.69,1115.683,138);
    McSourceReso.InitNumMcIter(400000);

    DLM_CleverMcLevyReso McSourceLevy;
    McSourceLevy.InitScale(1,SourceScale-0.01,SourceScale+0.01);
    McSourceLevy.InitStability(1,SourceStability-0.01,SourceStability+0.01);
    McSourceLevy.InitRad(NumRadPts,0,MaxRad);
    McSourceLevy.InitType(2);
    McSourceLevy.InitNumMcIter(400000);

    DLM_StableDistribution DlmStable;
    DlmStable.SetScale(SourceScale);
    DlmStable.SetStability(SourceStability);

    MS_GaussExp_mT_Simple MSGEMTS_pp;
    MSGEMTS_pp.SetNum_mT(1);
    MSGEMTS_pp.SetMean_mT(0,1);
    MSGEMTS_pp.SetWeight_mT(0,1);
    MSGEMTS_pp.SetCustomFunction(0,1.0);
    MSGEMTS_pp.SetMass(0,938.272);
    //MSGEMTS_pp.SetMass(1,938.272);
    MSGEMTS_pp.SetMass(1,1115.683);
    MSGEMTS_pp.SetMassR(0,1361.52);
    //MSGEMTS_pp.SetMassR(1,1361.52);
    MSGEMTS_pp.SetMassR(1,1462.93);
    MSGEMTS_pp.SetMassD(0,138);
    MSGEMTS_pp.SetMassD(1,138);
    MSGEMTS_pp.SetTau(0,1.65);
    //MSGEMTS_pp.SetTau(1,1.65);
    MSGEMTS_pp.SetTau(1,4.69);
    MSGEMTS_pp.SetResonanceWeight(0,1.-0.3578);
    //MSGEMTS_pp.SetResonanceWeight(1,1.-0.3578);
    MSGEMTS_pp.SetResonanceWeight(1,1.-0.3562);

    TFile* OutputFile = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/McLevyTest1/OutputFile.root","recreate");
    TH1F* hMcReso = new TH1F("hMcReso","hMcReso",NumRadPts,0,MaxRad);
    TH1F* hMcLevy = new TH1F("hMcLevy","hMcLevy",NumRadPts,0,MaxRad);
    TH1F* hOldCalc = new TH1F("hOldCalc","hOldCalc",NumRadPts,0,MaxRad);
    TH1F* hStabDist = new TH1F("hStabDist","hStabDist",NumRadPts,0,MaxRad);
    TH1F* hLevy = new TH1F("hLevy","hLevy",NumRadPts,0,MaxRad);
    TH1F* hRanGenLevy = new TH1F("hRanGenLevy","hRanGenLevy",NumRadPts,0,MaxRad);

    double PARS[6];
    double& Radius = PARS[1];
    double& Scale = PARS[3];
    double& Stability = PARS[4];
    PARS[5]=3;
    Scale = SourceScale;
    Stability = SourceStability;

    double PARSst[7];


    DLM_Random RanGen(7);

    for(unsigned uRad=0; uRad<NumRadPts; uRad++){
        Radius = hMcReso->GetBinCenter(uRad+1);
        hMcReso->SetBinContent(uRad+1,McSourceReso.Eval(PARS));
        hMcLevy->SetBinContent(uRad+1,McSourceLevy.Eval(PARS));
        hOldCalc->SetBinContent(uRad+1,MSGEMTS_pp.Eval(PARS));
        PARSst[1]=Radius; PARSst[3]=Stability; PARSst[4]=0; PARSst[5]=Scale; PARSst[6]=0;
        hStabDist->SetBinContent(uRad+1,DlmStable.Eval(PARSst));
        hLevy->SetBinContent(uRad+1,LevySource3D(PARS));
    }
    //for(unsigned uIter=0; uIter<400000; uIter++){
    //    hRanGenLevy->Fill(RanGen.StableNolan(3,Stability,0,Scale,0));
    //}
    //hRanGenLevy->Scale(1./400000.,"width");

    hMcReso->Write();
    hMcLevy->Write();
    hOldCalc->Write();
    hStabDist->Write();
    hLevy->Write();
    hRanGenLevy->Write();

    delete hMcReso;
    delete hMcLevy;
    delete hOldCalc;
    delete hLevy;
    delete hRanGenLevy;
    delete hStabDist;

    delete OutputFile;
}


void testCats() {

  //  const double Weight1S0 = 3. / 12.;

  //  const double Weight3P0 = 1. / 12.;

  //  const double Weight3P1 = 3. / 12.;

  //  const double Weight3P2 = 5. / 12.;

  //

  //  double PotPars1S0[8] = { NN_AV18, v18_Coupled3P2, 1, 1, 1, 0, 0, 0 };

  //  double PotPars3P0[8] = { NN_AV18, v18_Coupled3P2, 1, 1, 1, 1, 1, 0 };

  //  double PotPars3P1[8] = { NN_AV18, v18_Coupled3P2, 1, 1, 1, 1, 1, 1 };

  //  double PotPars3P2[8] = { NN_AV18, v18_Coupled3P2, 1, 1, 1, 1, 1, 2 };

  //

  //  CATSparameters *cPotPars1S0 = new CATSparameters(CATSparameters::tPotential,

  //                                                   8, true);

  //  cPotPars1S0->SetParameters(PotPars1S0);

  //  CATSparameters *cPotPars3P0 = new CATSparameters(CATSparameters::tPotential,

  //                                                   8, true);

  //  cPotPars3P0->SetParameters(PotPars3P0);

  //  CATSparameters *cPotPars3P1 = new CATSparameters(CATSparameters::tPotential,

  //                                                   8, true);

  //  cPotPars3P1->SetParameters(PotPars3P1);

  //  CATSparameters *cPotPars3P2 = new CATSparameters(CATSparameters::tPotential,

  //                                                   8, true);

  //  cPotPars3P2->SetParameters(PotPars3P2);



  unsigned int momBins = 50;

  const double kMin = 0;

  const double kMax = 200;



  CATS AB_pp;

  CATSparameters* cPars;

  cPars = new CATSparameters(CATSparameters::tSource, 1, true);

  cPars->SetParameter(0, 1.2);

  AB_pp.SetAnaSource(GaussSource, *cPars);

  AB_pp.SetUseAnalyticSource(true);

  AB_pp.SetMomentumDependentSource(false);

  AB_pp.SetThetaDependentSource(false);

  AB_pp.SetExcludeFailedBins(false);

  AB_pp.SetMomBins(momBins, kMin, kMax);

  AB_pp.SetQ1Q2(1);

  AB_pp.SetPdgId(2212, 2212);

  AB_pp.SetNumChannels(1);

  AB_pp.SetNumPW(0, 1);

  AB_pp.SetSpin(0, 0);

  //AB_pp.SetNumChannels(4);

  //  AB_pp->SetNumPW(0, 2);

  //  AB_pp->SetNumPW(1, 2);

  //  AB_pp->SetNumPW(2, 2);

  //  AB_pp->SetNumPW(3, 2);

  //  AB_pp->SetSpin(0, 0);

  //  AB_pp->SetSpin(1, 1);

  //  AB_pp->SetSpin(2, 1);

  //  AB_pp->SetSpin(3, 1);

  //  AB_pp->SetChannelWeight(0, Weight1S0);

  //  AB_pp->SetChannelWeight(1, Weight3P0);

  //  AB_pp->SetChannelWeight(2, Weight3P1);

  //  AB_pp->SetChannelWeight(3, Weight3P2);

  //  AB_pp->SetShortRangePotential(0, 0, fDlmPot, *cPotPars1S0);

  //  AB_pp->SetShortRangePotential(1, 1, fDlmPot, *cPotPars3P0);

  //  AB_pp->SetShortRangePotential(2, 1, fDlmPot, *cPotPars3P1);

  //  AB_pp->SetShortRangePotential(3, 1, fDlmPot, *cPotPars3P2);

    AB_pp.KillTheCat();

  return;

}

//[0] radius, [1] momentum
//[2] IsoSpin, [3] Spin
double ESC16_pXim_EXAMPLE(double* Parameters){

    if(Parameters[0]>2) return 0;

    const TString FileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/Tom/Potential/ESC16_190319.root";
    const double& Isospin = Parameters[2];
    const double& Spin = Parameters[3];

    TFile fIn(FileName,"read");
    TGraph* gPot;

    if(Isospin==0&&Spin==0){
        gPot = (TGraph*)fIn.Get("EffPotI0S0");
    }
    else if(Isospin==0&&Spin==1){
        gPot = (TGraph*)fIn.Get("EffPotI0S1");
    }
    else if(Isospin==1&&Spin==0){
        gPot = (TGraph*)fIn.Get("EffPotI1S0");
    }
    else if(Isospin==1&&Spin==1){
        gPot = (TGraph*)fIn.Get("EffPotI1S1");
    }
    else{
        printf(" ESC16_pXim_EXAMPLE says ????????????????");
        return 0;
    }

    return gPot->Eval(Parameters[0]);
}


void TestTomPotential(){

    CATSparameters PotParsI0S0(CATSparameters::tSource,2,true);
    PotParsI0S0.SetParameter(0,0);
    PotParsI0S0.SetParameter(1,0);

    CATSparameters PotParsI0S1(CATSparameters::tSource,2,true);
    PotParsI0S1.SetParameter(0,0);
    PotParsI0S1.SetParameter(1,1);

    CATSparameters PotParsI1S0(CATSparameters::tSource,2,true);
    PotParsI1S0.SetParameter(0,1);
    PotParsI1S0.SetParameter(1,0);

    CATSparameters PotParsI1S1(CATSparameters::tSource,2,true);
    PotParsI1S1.SetParameter(0,1);
    PotParsI1S1.SetParameter(1,1);

    CATS Kitty;
    //define the usual stuff
    //Kitty.SetShortRangePotential(0,0,ESC16_pXim_EXAMPLE,&PotParsI0S0);
    //Kitty.SetShortRangePotential(1,0,ESC16_pXim_EXAMPLE,&PotParsI0S1);
    //Kitty.SetShortRangePotential(2,0,ESC16_pXim_EXAMPLE,&PotParsI1S0);
    //Kitty.SetShortRangePotential(3,0,ESC16_pXim_EXAMPLE,&PotParsI1S1);
}

void TestTomWF1(){
    const unsigned NumMomBins = 35;
    const double kMin = 0;
    const double kMax = 350;
    DLM_Histo<complex<double>>*** ExternalWF=NULL;
    CATS Kitty;
    Kitty.SetMomBins(NumMomBins,kMin,kMax);

    CATSparameters cPars(CATSparameters::tSource,1,true);
    cPars.SetParameter(0,0.8);

    Kitty.SetUseAnalyticSource(true);
    Kitty.SetAnaSource(GaussSource,cPars);
    ExternalWF = Init_pXi_ESC16_Iavg_Coulomb("//home/dmihaylov/CernBox/CATS_potentials/Tom/WithCoulomb_Iavg/",Kitty);
    Kitty.SetExternalWaveFunction(0,0,ExternalWF[0][0][0],ExternalWF[1][0][0]);
    Kitty.SetExternalWaveFunction(1,0,ExternalWF[0][1][0],ExternalWF[1][1][0]);

    Kitty.SetMomentumDependentSource(false);
    Kitty.SetThetaDependentSource(false);
    Kitty.SetExcludeFailedBins(false);

    Kitty.KillTheCat();

    TGraph gKitty;
    gKitty.Set(Kitty.GetNumMomBins());
    gKitty.SetName("gKitty");

    for(unsigned uBin=0; uBin<Kitty.GetNumMomBins(); uBin++){
        printf("C(%.2f) = %.3f\n",Kitty.GetMomentum(uBin),Kitty.GetCorrFun(uBin));
        gKitty.SetPoint(uBin,Kitty.GetMomentum(uBin),Kitty.GetCorrFun(uBin));
    }

    TFile* OutputFile = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/TestTomWF1/OutputFile.root","recreate");
    gKitty.Write();

    delete OutputFile;
}


void TestHaide_pL_pWaves(){

    DLM_CommonAnaFunctions AnalysisObject;

    double* MomBins_pL = NULL;
    double* FitRegion_pL = NULL;
    unsigned NumMomBins_pL=30;
    double kMin=0;
    double kMax=300;
    CATS AB_pL;
    AB_pL.SetMomBins(NumMomBins_pL,kMin,kMax);
    AnalysisObject.SetUpCats_pL(AB_pL,"NLO_sp","Gauss");
    AB_pL.KillTheCat();

    for(unsigned uBin=0; uBin<AB_pL.GetNumMomBins(); uBin++){
        printf("C(%.2f) = %.4f\n",AB_pL.GetMomentum(uBin),AB_pL.GetCorrFun(uBin));
        //gKitty.SetPoint(uBin,AB_pL.GetMomentum(uBin),AB_pL.GetCorrFun(uBin));
    }

    delete [] MomBins_pL;
    delete [] FitRegion_pL;
}

void TestSigma0(){
    const unsigned NumMomBins = 35;
    const double kMin = 0;
    const double kMax = 350;
    DLM_Histo<complex<double>>*** ExternalWF=NULL;
    CATS Kitty;
    Kitty.SetMomBins(NumMomBins,kMin,kMax);

    CATSparameters cPars(CATSparameters::tSource,1,true);
    cPars.SetParameter(0,1.0);

    Kitty.SetUseAnalyticSource(true);
    Kitty.SetAnaSource(GaussSource,cPars);
    ExternalWF = Init_pS0_ESC08("//home/dmihaylov/CernBox/CATS_potentials/Tom/pSigma0/030419/",Kitty);
    Kitty.SetExternalWaveFunction(0,0,ExternalWF[0][0][0],ExternalWF[1][0][0]);
    Kitty.SetExternalWaveFunction(1,0,ExternalWF[0][1][0],ExternalWF[1][1][0]);

    Kitty.SetMomentumDependentSource(false);
    Kitty.SetThetaDependentSource(false);
    Kitty.SetExcludeFailedBins(false);
//Kitty.SetChannelWeight(0,0);
//Kitty.SetChannelWeight(1,1);
    Kitty.KillTheCat();

    TFile* OutputFile = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/TestSigma0/OutputFile.root","recreate");

    TGraph gKitty;
    gKitty.Set(Kitty.GetNumMomBins());
    gKitty.SetName("gKitty");

    for(unsigned uBin=0; uBin<Kitty.GetNumMomBins(); uBin++){
        printf("C(%.2f) = %.3f\n",Kitty.GetMomentum(uBin),Kitty.GetCorrFun(uBin));
        gKitty.SetPoint(uBin,Kitty.GetMomentum(uBin),Kitty.GetCorrFun(uBin));

        TGraph gWF1S0;
        gWF1S0.SetName(TString::Format("gWF1S0_%.0f",Kitty.GetMomentum(uBin)));
        TGraph gWF3S1;
        gWF3S1.SetName(TString::Format("gWF3S1_%.0f",Kitty.GetMomentum(uBin)));
        int COUNTER=0;
        for(double RAD=0.05; RAD<15; RAD+=0.05){
            gWF1S0.SetPoint(COUNTER,RAD,real(Kitty.EvalRadialWaveFunction(uBin,0,0,RAD,false)));
            gWF3S1.SetPoint(COUNTER,RAD,real(Kitty.EvalRadialWaveFunction(uBin,1,0,RAD,false)));
            COUNTER++;
        }
        OutputFile->cd();
        gWF1S0.Write();
        gWF3S1.Write();
    }

    OutputFile->cd();
    gKitty.Write();

    delete OutputFile;

}

int main(int argc, char *argv[])
{

    TFile* f1 = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Temp/temp.root","recreate");
    TNtuple* nt1 = new TNtuple("nt1","nt1",
                                "UITERfirst:UITERlast:f0Inv:d0:Chi2fit:NDFfit:Chi2refit:dChi2:dChi2Ndf");

    //omp_set_dynamic(0);     // Explicitly disable dynamic teams
    //omp_set_num_threads(1); // Use 4 threads for all consecutive parallel regions

    DLM_Timer TIMER;
    TCanvas c1;

    char** ARGV = NULL;
    if(argc) ARGV=new char* [argc];
    for(int iARG=1; iARG<argc; iARG++){
        ARGV[iARG] = new char [128];
        strcpy(ARGV[iARG],argv[iARG]);
    }

    //TestHaide_pL_pWaves();

    //TestTomPotential();

    //ExecuteCFmT(argc,ARGV);

    //DimiExecuteCFmT();
    //TestTomWF1();
/*
DLM_CommonAnaFunctions AO;
double lambdapars[5];

AO.SetUpLambdaPars_pp("pp13TeV_HM_March19",0,lambdapars);
printf("%.3f\n",lambdapars[0]*100.);
printf("%.3f\n",lambdapars[1]*100.);
printf("%.3f\n",lambdapars[2]*100.);
printf("%.3f\n\n",lambdapars[3]*100.);

AO.SetUpLambdaPars_pL("pp13TeV_HM_March19",0,0,lambdapars);
printf("%.3f\n",lambdapars[0]*100.);
printf("%.3f\n",lambdapars[1]*100.);
printf("%.3f\n",lambdapars[2]*100.);
printf("%.3f\n",lambdapars[3]*100.);
printf("%.3f\n",lambdapars[4]*100.);
*/

    PLAMBDA_1_MAIN(argc,ARGV);
    //CALL_BERNIE_AND_VALE();
    //FEMTOBOYZ_MAIN(argc,ARGV);
    //GerhardMAIN(argc,ARGV);


    //TestCommonInit();
    //plot_pp();
    //plot_pL();
    //plot_pSigma0();
    //McLevyTest1();
    //testCats();
    //TestSigma0();

    //TestCATS3_NewExtWf("Usmani");
    //TestCATS3_NewExtWf("NLO");
    //TestCATS3_NewExtWf("NLO_Coupled_S");
    for(int iARG=1; iARG<argc; iARG++){
        delete [] ARGV[iARG];
    }
    if(ARGV) delete [] ARGV;

    long long ExeTime = TIMER.Stop()/1000.;
    char* strtime = new char [128];
    ShowTime(ExeTime,strtime,0,true,6);
    printf("The script terminated after: %s\n",strtime);

    delete [] strtime;

    delete nt1;
    delete f1;
    return 0;
}



