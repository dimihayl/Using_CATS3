
#include "CommonAnaFunctions.h"
#include "GentleDimi.h"
#include "FemtoBoyzScripts.h"
#include "pSigma.h"

#include <iostream>
#include <stdio.h>
#include <string.h>
//#include <omp.h>
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
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveText.h"

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
/*
void plot_pSigma0_Tom(){

    const double SourceSize = 1.2;
    const double SourceStability = 2.0;

    const double kMin=0;
    const double kMax=320;
    const unsigned NumMomBins = 64;

    TFile* OutputFile = new TFile(
                TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/plot_pSigma0_Tom/OutputFile_%.2f_%.2f.root",
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
    ExternalWF = Init_pS0_ESC08("/home/dmihaylov/CernBox/CATS_potentials/Tom/pSigma0/260419/",Kitty);
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
*/

void TestQS(){

    CATS Kitty;
    CATSparameters* cPars;

    const unsigned NumMomBins = 256;

    cPars = new CATSparameters(CATSparameters::tSource, 1, true);
    cPars->SetParameter(0, 0.5);
    Kitty.SetAnaSource(GaussSource, *cPars);
    Kitty.SetUseAnalyticSource(true);
    Kitty.SetMomentumDependentSource(false);
    Kitty.SetThetaDependentSource(false);
    Kitty.SetExcludeFailedBins(false);
    Kitty.SetMomBins(NumMomBins, 0, 512);
    Kitty.SetQ1Q2(0);
    Kitty.SetPdgId(2212, 2212);
    Kitty.SetNumChannels(1);
    Kitty.SetNumPW(0, 1);
    Kitty.SetSpin(0, 1);
    Kitty.SetChannelWeight(0,1);
    Kitty.KillTheCat();

    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        printf("C(%.3f) = %.3f\n",Kitty.GetMomentum(uBin),Kitty.GetCorrFun(uBin));
    }

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
    ExternalWF = Init_pS0_ESC08("/home/dmihaylov/CernBox/CATS_potentials/Tom/pSigma0/WithPhase/",Kitty);
    Kitty.SetExternalWaveFunction(0,0,ExternalWF[0][0][0],ExternalWF[1][0][0]);
    Kitty.SetExternalWaveFunction(1,0,ExternalWF[0][1][0],ExternalWF[1][1][0]);

    Kitty.SetMomentumDependentSource(false);
    Kitty.SetThetaDependentSource(false);
    Kitty.SetExcludeFailedBins(false);
Kitty.SetChannelWeight(0,0.25);
Kitty.SetChannelWeight(1,0.75);
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
            gWF1S0.SetPoint(COUNTER,RAD,abs(Kitty.EvalRadialWaveFunction(uBin,0,0,RAD,false)));
            gWF3S1.SetPoint(COUNTER,RAD,abs(Kitty.EvalRadialWaveFunction(uBin,1,0,RAD,false)));
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

void TestRandBW(){

    const unsigned NumIter = 100000;

    TH1F* hBW_root = new TH1F("hBW_root","hBW_root",256,-8,8);
    TH1F* hBW_dimi = new TH1F("hBW_dimi","hBW_dimi",256,-8,8);

    TRandom3 RanRoot(11);
    DLM_Random RanDimi(11);

    double Mean = 0;
    double Gamma = 2;

    for(unsigned uIter=0; uIter<NumIter; uIter++){
        hBW_root->Fill(RanRoot.BreitWigner(Mean,Gamma));
        hBW_dimi->Fill(RanDimi.Cauchy(Mean,Gamma/sqrt(2)));
    }
    hBW_root->SetLineColor(kRed);
    hBW_root->SetLineWidth(3);
    hBW_root->Sumw2();

    hBW_dimi->SetLineColor(kBlue);
    hBW_dimi->SetLineWidth(2);
    hBW_dimi->Sumw2();

    TCanvas* cBW = new TCanvas("cBW", "cBW", 1);
    cBW->cd(0); cBW->SetCanvasSize(960, 540); cBW->SetMargin(0.15,0.05,0.2,0.05);//lrbt

    hBW_root->Draw();
    hBW_dimi->Draw("same");

    cBW->SaveAs("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/TestRandBW/cBW.png");

    delete hBW_root;
    delete hBW_dimi;
    delete cBW;

}

void SourceSmearing(){

    unsigned NumMassSmear = 2;
    bool MassSmear[NumMassSmear] = {false,true};

    //unsigned NumMomSmear = 5;
    //double MomSmear[NumMomSmear] = {0,0.05,0.10,0.15,0.20};

    //unsigned NumMomSmear = 3;
    //double MomSmear[NumMomSmear] = {0,0.05,0.10};

    unsigned NumMomSmear = 2;
    double MomSmear[NumMomSmear] = {0,0.20};

    unsigned NumIter = NumMassSmear*NumMomSmear;

    DLM_CommonAnaFunctions AnaObject;

    const unsigned NumMomBins = 50;
    const double kMin = 0;
    const double kMax = 250;

    double GaussSourceSize_pp = 1.32;//1.35
    double GaussSourceSize_pL = 1.22;//1.55
    double GaussSourceSize_pXim = 1.09;//1.85
    double GaussSourceSize_pOmegam = 1.01;//2.17

    const unsigned NumSourceBins = 128;
    const double rMin = 0.1;
    const double rMax = 8.0;
    const double rStep = (rMax-rMin)/double(NumSourceBins);

    CATS AB_pp;
    AB_pp.SetMomBins(NumMomBins,kMin,kMax);
    AnaObject.SetUpCats_pp(AB_pp,"AV18","McGauss_Reso");
    AB_pp.SetAnaSource(0,0.97);
    AB_pp.SetEpsilonConv(1e-8);
    AB_pp.SetEpsilonProp(1e-8);

    CATS AB_pL;
    AB_pL.SetMomBins(NumMomBins,kMin,kMax);
    AnaObject.SetUpCats_pL(AB_pL,"NLO_Coupled_S","McGauss_Reso");
    AB_pL.SetAnaSource(0,0.92);
    AB_pL.SetEpsilonConv(1e-8);
    AB_pL.SetEpsilonProp(1e-8);

    CATS AB_pXim;
    AB_pXim.SetMomBins(NumMomBins,kMin,kMax);
    AnaObject.SetUpCats_pXim(AB_pXim,"pXim_HALQCD1","McGauss_Reso");
    AB_pXim.SetAnaSource(0,0.79);
    AB_pXim.SetEpsilonConv(1e-8);
    AB_pXim.SetEpsilonProp(1e-8);

    CATS AB_pOmegam;
    AB_pOmegam.SetMomBins(NumMomBins,kMin,kMax);
    AnaObject.SetUpCats_pOmegam(AB_pOmegam,"pOmega_Lattice","McGauss_Reso");
    AB_pOmegam.SetAnaSource(0,0.73);
    AB_pOmegam.SetEpsilonConv(1e-8);
    AB_pOmegam.SetEpsilonProp(1e-8);

    AB_pp.SetNotifications(CATS::nWarning);
    AB_pp.SetMaxNumThreads(4);
    AB_pL.SetNotifications(CATS::nWarning);
    AB_pL.SetMaxNumThreads(4);
    AB_pXim.SetNotifications(CATS::nWarning);
    AB_pXim.SetMaxNumThreads(4);
    AB_pOmegam.SetNotifications(CATS::nWarning);
    AB_pOmegam.SetMaxNumThreads(4);

    TGraph* gCk_pp = new TGraph[NumIter];
    TGraph* gCk_pL = new TGraph[NumIter];
    TGraph* gCk_pXim = new TGraph[NumIter];
    TGraph* gCk_pOmegam = new TGraph[NumIter];

    TGraph* gCkRatio_pp = new TGraph[NumIter];
    TGraph* gCkRatio_pL = new TGraph[NumIter];
    TGraph* gCkRatio_pXim = new TGraph[NumIter];
    TGraph* gCkRatio_pOmegam = new TGraph[NumIter];

    TGraph gpp_low;
    gpp_low.SetName("gpp_low");
    TGraph gpp_up;
    gpp_up.SetName("gpp_up");

    TGraph gpL_low;
    gpL_low.SetName("gpL_low");
    TGraph gpL_up;
    gpL_up.SetName("gpL_up");

    TGraph gpXim_low;
    gpXim_low.SetName("gpXim_low");
    TGraph gpXim_up;
    gpXim_up.SetName("gpXim_up");

    TGraph gpOmegam_low;
    gpOmegam_low.SetName("gpOmegam_low");
    TGraph gpOmegam_up;
    gpOmegam_up.SetName("gpOmegam_up");

    TGraph gSourceGaussReso_pp;
    TGraph gSourceGaussReso_pL;
    TGraph gSourceGaussReso_pXim;
    TGraph gSourceGaussReso_pOmegam;

    TGraph gSourceGauss_pp;
    TGraph gSourceGauss_pL;
    TGraph gSourceGauss_pXim;
    TGraph gSourceGauss_pOmegam;

    gSourceGaussReso_pp.SetName("gSourceGaussReso_pp");
    gSourceGaussReso_pL.SetName("gSourceGaussReso_pL");
    gSourceGaussReso_pXim.SetName("gSourceGaussReso_pXim");
    gSourceGaussReso_pOmegam.SetName("gSourceGaussReso_pOmegam");

    gSourceGauss_pp.SetName("gSourceGauss_pp");
    gSourceGauss_pL.SetName("gSourceGauss_pL");
    gSourceGauss_pXim.SetName("gSourceGauss_pXim");
    gSourceGauss_pOmegam.SetName("gSourceGauss_pOmegam");

    unsigned Counter=0;

    //double MaxLow_pp=1;
    //double MaxUp_pp=1;

    //double MaxLow_pL=1;
    //double MaxUp_pL=1;

    //double MaxLow_pXim=1;
    //double MaxUp_pXim=1;

    //double MaxLow_pOmegam=1;
    //double MaxUp_pOmegam=1;

    double xVal,yVal,xVal2,yVal2;

    for(unsigned uMassSmear=0; uMassSmear<NumMassSmear; uMassSmear++){
        printf("uMassSmear = %u\n",uMassSmear);
        for(unsigned uMomSmear=0; uMomSmear<NumMomSmear; uMomSmear++){
            printf(" uMomSmear = %u\n",uMomSmear);

            gCk_pp[Counter].SetName(TString::Format("gCk_pp_%u_%u",uMassSmear,uMomSmear));
            gCk_pL[Counter].SetName(TString::Format("gCk_pL_%u_%u",uMassSmear,uMomSmear));
            gCk_pXim[Counter].SetName(TString::Format("gCk_pXim_%u_%u",uMassSmear,uMomSmear));
            gCk_pOmegam[Counter].SetName(TString::Format("gCk_pOmegam_%u_%u",uMassSmear,uMomSmear));

            gCkRatio_pp[Counter].SetName(TString::Format("gCkRatio_pp_%u_%u",uMassSmear,uMomSmear));
            gCkRatio_pL[Counter].SetName(TString::Format("gCkRatio_pL_%u_%u",uMassSmear,uMomSmear));
            gCkRatio_pXim[Counter].SetName(TString::Format("gCkRatio_pXim_%u_%u",uMassSmear,uMomSmear));
            gCkRatio_pOmegam[Counter].SetName(TString::Format("gCkRatio_pOmegam_%u_%u",uMassSmear,uMomSmear));

            AnaObject.GetCleverMcLevyReso_pp()->SetUpReso(0,0,1.-0.3578,1361.52,1.65,Mass_p,Mass_pic,MomSmear[uMomSmear],MassSmear[uMassSmear]);
            AnaObject.GetCleverMcLevyReso_pp()->SetUpReso(1,0,1.-0.3578,1361.52,1.65,Mass_p,Mass_pic,MomSmear[uMomSmear],MassSmear[uMassSmear]);

            AnaObject.GetCleverMcLevyReso_pL()->SetUpReso(0,0,1.-0.3578,1361.52,1.65,Mass_p,Mass_pic,MomSmear[uMomSmear],MassSmear[uMassSmear]);
            AnaObject.GetCleverMcLevyReso_pL()->SetUpReso(1,0,1.-0.3562,1462.93,4.69,Mass_L,Mass_pic,MomSmear[uMomSmear],MassSmear[uMassSmear]);

            AnaObject.GetCleverMcLevyReso_pXim()->SetUpReso(0,0,1.-0.3578,1361.52,1.65,Mass_p,Mass_pic,MomSmear[uMomSmear],MassSmear[uMassSmear]);

            AnaObject.GetCleverMcLevyReso_pOmegam()->SetUpReso(0,0,1.-0.3578,1361.52,1.65,Mass_p,Mass_pic,MomSmear[uMomSmear],MassSmear[uMassSmear]);

            AB_pp.KillTheCat(uMassSmear+uMomSmear?CATS::kSourceChanged:CATS::kAllChanged);
            AB_pL.KillTheCat(uMassSmear+uMomSmear?CATS::kSourceChanged:CATS::kAllChanged);
            AB_pXim.KillTheCat(uMassSmear+uMomSmear?CATS::kSourceChanged:CATS::kAllChanged);
            AB_pOmegam.KillTheCat(uMassSmear+uMomSmear?CATS::kSourceChanged:CATS::kAllChanged);

            for(unsigned uBin=0; uBin<AB_pp.GetNumMomBins(); uBin++){
                gCk_pp[Counter].SetPoint(uBin,AB_pp.GetMomentum(uBin),AB_pp.GetCorrFun(uBin));
                if(uMassSmear==0&&uMomSmear==0){
                    gpp_low.SetPoint(uBin,AB_pp.GetMomentum(uBin),AB_pp.GetCorrFun(uBin));
                    gpp_up.SetPoint(uBin,AB_pp.GetMomentum(uBin),AB_pp.GetCorrFun(uBin));
                }
                else{
                    gpp_low.GetPoint(uBin,xVal,yVal);
                    if(AB_pp.GetCorrFun(uBin)<yVal){gpp_low.SetPoint(uBin,AB_pp.GetMomentum(uBin),AB_pp.GetCorrFun(uBin));}
                    gpp_up.GetPoint(uBin,xVal,yVal);
                    if(AB_pp.GetCorrFun(uBin)>yVal){gpp_up.SetPoint(uBin,AB_pp.GetMomentum(uBin),AB_pp.GetCorrFun(uBin));}
                }
            }
            for(unsigned uBin=0; uBin<AB_pL.GetNumMomBins(); uBin++){
                gCk_pL[Counter].SetPoint(uBin,AB_pL.GetMomentum(uBin),AB_pL.GetCorrFun(uBin));
                if(uMassSmear==0&&uMomSmear==0){
                    gpL_low.SetPoint(uBin,AB_pL.GetMomentum(uBin),AB_pL.GetCorrFun(uBin));
                    gpL_up.SetPoint(uBin,AB_pL.GetMomentum(uBin),AB_pL.GetCorrFun(uBin));
                }
                else{
                    gpL_low.GetPoint(uBin,xVal,yVal);
                    if(AB_pL.GetCorrFun(uBin)<yVal){gpL_low.SetPoint(uBin,AB_pL.GetMomentum(uBin),AB_pL.GetCorrFun(uBin));}
                    gpL_up.GetPoint(uBin,xVal,yVal);
                    if(AB_pL.GetCorrFun(uBin)>yVal){gpL_up.SetPoint(uBin,AB_pL.GetMomentum(uBin),AB_pL.GetCorrFun(uBin));}
                }
            }
            for(unsigned uBin=0; uBin<AB_pXim.GetNumMomBins(); uBin++){
                gCk_pXim[Counter].SetPoint(uBin,AB_pXim.GetMomentum(uBin),AB_pXim.GetCorrFun(uBin));
                if(uMassSmear==0&&uMomSmear==0){
                    gpXim_low.SetPoint(uBin,AB_pXim.GetMomentum(uBin),AB_pXim.GetCorrFun(uBin));
                    gpXim_up.SetPoint(uBin,AB_pXim.GetMomentum(uBin),AB_pXim.GetCorrFun(uBin));
                }
                else{
                    gpXim_low.GetPoint(uBin,xVal,yVal);
                    if(AB_pXim.GetCorrFun(uBin)<yVal){gpXim_low.SetPoint(uBin,AB_pXim.GetMomentum(uBin),AB_pXim.GetCorrFun(uBin));}
                    gpXim_up.GetPoint(uBin,xVal,yVal);
                    if(AB_pXim.GetCorrFun(uBin)>yVal){gpXim_up.SetPoint(uBin,AB_pXim.GetMomentum(uBin),AB_pXim.GetCorrFun(uBin));}
                }
            }
            for(unsigned uBin=0; uBin<AB_pOmegam.GetNumMomBins(); uBin++){
                gCk_pOmegam[Counter].SetPoint(uBin,AB_pOmegam.GetMomentum(uBin),AB_pOmegam.GetCorrFun(uBin));
                if(uMassSmear==0&&uMomSmear==0){
                    gpOmegam_low.SetPoint(uBin,AB_pOmegam.GetMomentum(uBin),AB_pOmegam.GetCorrFun(uBin));
                    gpOmegam_up.SetPoint(uBin,AB_pOmegam.GetMomentum(uBin),AB_pOmegam.GetCorrFun(uBin));
                }
                else{
                    gpOmegam_low.GetPoint(uBin,xVal,yVal);
                    if(AB_pOmegam.GetCorrFun(uBin)<yVal){gpOmegam_low.SetPoint(uBin,AB_pOmegam.GetMomentum(uBin),AB_pOmegam.GetCorrFun(uBin));}
                    gpOmegam_up.GetPoint(uBin,xVal,yVal);
                    if(AB_pOmegam.GetCorrFun(uBin)>yVal){gpOmegam_up.SetPoint(uBin,AB_pOmegam.GetMomentum(uBin),AB_pOmegam.GetCorrFun(uBin));}
                }
            }

            for(unsigned uBin=0; uBin<AB_pp.GetNumMomBins(); uBin++){
                gCk_pp[0].GetPoint(uBin,xVal,yVal);
                gCk_pp[Counter].GetPoint(uBin,xVal2,yVal2);
                gCkRatio_pp[Counter].SetPoint(uBin,AB_pp.GetMomentum(uBin),yVal2/yVal);
            }
            for(unsigned uBin=0; uBin<AB_pL.GetNumMomBins(); uBin++){
                gCk_pL[0].GetPoint(uBin,xVal,yVal);
                gCk_pL[Counter].GetPoint(uBin,xVal2,yVal2);
                gCkRatio_pL[Counter].SetPoint(uBin,AB_pL.GetMomentum(uBin),yVal2/yVal);
            }
            for(unsigned uBin=0; uBin<AB_pXim.GetNumMomBins(); uBin++){
                gCk_pXim[0].GetPoint(uBin,xVal,yVal);
                gCk_pXim[Counter].GetPoint(uBin,xVal2,yVal2);
                gCkRatio_pXim[Counter].SetPoint(uBin,AB_pXim.GetMomentum(uBin),yVal2/yVal);
            }
            for(unsigned uBin=0; uBin<AB_pOmegam.GetNumMomBins(); uBin++){
                gCk_pOmegam[0].GetPoint(uBin,xVal,yVal);
                gCk_pOmegam[Counter].GetPoint(uBin,xVal2,yVal2);
                gCkRatio_pOmegam[Counter].SetPoint(uBin,AB_pOmegam.GetMomentum(uBin),yVal2/yVal);
            }

            Counter++;
        }
    }

    for(unsigned uBin=0; uBin<AB_pp.GetNumMomBins(); uBin++){
        gCk_pp[0].GetPoint(uBin,xVal,yVal);
        gpp_low.GetPoint(uBin,xVal2,yVal2);
        gpp_low.SetPoint(uBin,xVal,yVal2/yVal);

        gpp_up.GetPoint(uBin,xVal2,yVal2);
        gpp_up.SetPoint(uBin,xVal,yVal2/yVal);

    }
    for(unsigned uBin=0; uBin<AB_pL.GetNumMomBins(); uBin++){
        gCk_pL[0].GetPoint(uBin,xVal,yVal);
        gpL_low.GetPoint(uBin,xVal2,yVal2);
        gpL_low.SetPoint(uBin,xVal,yVal2/yVal);

        gpL_up.GetPoint(uBin,xVal2,yVal2);
        gpL_up.SetPoint(uBin,xVal,yVal2/yVal);
    }
    for(unsigned uBin=0; uBin<AB_pXim.GetNumMomBins(); uBin++){
        gCk_pXim[0].GetPoint(uBin,xVal,yVal);
        gpXim_low.GetPoint(uBin,xVal2,yVal2);
        gpXim_low.SetPoint(uBin,xVal,yVal2/yVal);

        gpXim_up.GetPoint(uBin,xVal2,yVal2);
        gpXim_up.SetPoint(uBin,xVal,yVal2/yVal);
    }
    for(unsigned uBin=0; uBin<AB_pOmegam.GetNumMomBins(); uBin++){
        gCk_pOmegam[0].GetPoint(uBin,xVal,yVal);
        gpOmegam_low.GetPoint(uBin,xVal2,yVal2);
        gpOmegam_low.SetPoint(uBin,xVal,yVal2/yVal);

        gpOmegam_up.GetPoint(uBin,xVal2,yVal2);
        gpOmegam_up.SetPoint(uBin,xVal,yVal2/yVal);
    }

    TFile* fOutputFile = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/SourceSmearing/fOutputFile.root","recreate");
    for(unsigned uIter=0; uIter<NumIter; uIter++){
        gCk_pp[uIter].Write();
    }
    for(unsigned uIter=0; uIter<NumIter; uIter++){
        gCkRatio_pp[uIter].Write();
    }
    for(unsigned uIter=0; uIter<NumIter; uIter++){
        gCk_pL[uIter].Write();
    }
    for(unsigned uIter=0; uIter<NumIter; uIter++){
        gCkRatio_pL[uIter].Write();
    }
    for(unsigned uIter=0; uIter<NumIter; uIter++){
        gCk_pXim[uIter].Write();
    }
    for(unsigned uIter=0; uIter<NumIter; uIter++){
        gCkRatio_pXim[uIter].Write();
    }
    for(unsigned uIter=0; uIter<NumIter; uIter++){
        gCk_pOmegam[uIter].Write();
    }
    for(unsigned uIter=0; uIter<NumIter; uIter++){
        gCkRatio_pOmegam[uIter].Write();
    }

    for(unsigned uRad=0; uRad<NumSourceBins; uRad++){
        double PARS[4];
        PARS[1] = rMin+double(uRad)*rStep;

        PARS[3] = GaussSourceSize_pp;
        gSourceGauss_pp.SetPoint(uRad,PARS[1],GaussSource(PARS));
        gSourceGaussReso_pp.SetPoint(uRad,PARS[1],AB_pp.EvaluateTheSource(10,PARS[1],0));

        PARS[3] = GaussSourceSize_pL;
        gSourceGauss_pL.SetPoint(uRad,PARS[1],GaussSource(PARS));
        gSourceGaussReso_pL.SetPoint(uRad,PARS[1],AB_pL.EvaluateTheSource(10,PARS[1],0));

        PARS[3] = GaussSourceSize_pXim;
        gSourceGauss_pXim.SetPoint(uRad,PARS[1],GaussSource(PARS));
        gSourceGaussReso_pXim.SetPoint(uRad,PARS[1],AB_pXim.EvaluateTheSource(10,PARS[1],0));

        PARS[3] = GaussSourceSize_pOmegam;
        gSourceGauss_pOmegam.SetPoint(uRad,PARS[1],GaussSource(PARS));
        gSourceGaussReso_pOmegam.SetPoint(uRad,PARS[1],AB_pOmegam.EvaluateTheSource(10,PARS[1],0));
    }


    TH1F* hAxis = new TH1F("hAxis", "hAxis", 256, kMin, kMax);
    hAxis->SetStats(false);
    hAxis->SetTitle("");
    hAxis->GetXaxis()->SetTitle("k* (MeV/c)");
    hAxis->GetXaxis()->SetTitleSize(0.06);
    hAxis->GetXaxis()->SetLabelSize(0.06);
    hAxis->GetXaxis()->CenterTitle();
    hAxis->GetXaxis()->SetTitleOffset(1.3);
    hAxis->GetXaxis()->SetLabelOffset(0.02);
    hAxis->GetYaxis()->SetTitle("C_{Smear}(k*)/C_{Default}(k*)");
    hAxis->GetYaxis()->SetTitleSize(0.06);
    hAxis->GetYaxis()->SetLabelSize(0.06);
    hAxis->GetYaxis()->CenterTitle();
    hAxis->GetYaxis()->SetTitleOffset(1.00);
    hAxis->GetYaxis()->SetRangeUser(0.975,1.025);
    hAxis->GetYaxis()->SetNdivisions(404);

    TCanvas* cRatio_pp = new TCanvas("cRatio_pp", "cRatio_pp", 1);
    cRatio_pp->cd(0); cRatio_pp->SetCanvasSize(1920, 1080); cRatio_pp->SetMargin(0.15,0.05,0.2,0.05);//lrbt
    cRatio_pp->SetGrid(true);
    TLegend* lLegend_pp = new TLegend(0.5,0.75,0.95,0.95);//lbrt
    lLegend_pp->SetName(TString::Format("lLegend"));
    lLegend_pp->SetTextSize(0.045);
    lLegend_pp->AddEntry(&gCkRatio_pp[NumMomSmear],"Effect of resonance width");
    lLegend_pp->AddEntry(&gCkRatio_pp[NumMomSmear-1],"Momentum smearing of 20%");
    TPaveText* PT_pp = new TPaveText(0.2,0.85,0.5,0.95, "blNDC");//lbrt
    PT_pp->SetName("PT_pp");
    PT_pp->SetBorderSize(1);
    PT_pp->SetTextSize(0.045);
    PT_pp->SetFillColor(kWhite);
    PT_pp->SetTextFont(22);
    PT_pp->AddText("p-p correlation");

    hAxis->Draw("axis");
    gCkRatio_pp[NumMomSmear].SetLineWidth(6);
    gCkRatio_pp[NumMomSmear].SetLineColor(kBlue);
    gCkRatio_pp[NumMomSmear-1].SetLineWidth(6);
    gCkRatio_pp[NumMomSmear-1].SetLineColor(kRed);

    gCkRatio_pp[NumMomSmear].Draw("same,c");//only mass smearing
    gCkRatio_pp[NumMomSmear-1].Draw("same,c");//max momentum smear
    lLegend_pp->Draw("same");
    PT_pp->Draw("same");

    gStyle->SetLineWidth(2.5);
    cRatio_pp->SaveAs(TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/SourceSmearing/cRatio_pp.png"));

	gStyle->SetLineWidth(1);
	gCkRatio_pp[NumMomSmear].SetLineWidth(gCkRatio_pp[NumMomSmear].GetLineWidth()/2.5);
    gCkRatio_pp[NumMomSmear-1].SetLineWidth(gCkRatio_pp[NumMomSmear-1].GetLineWidth()/2.5);
    cRatio_pp->SaveAs(TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/SourceSmearing/cRatio_pp.pdf"));
	gStyle->SetLineWidth(2.5);
	gCkRatio_pp[NumMomSmear].SetLineWidth(gCkRatio_pp[NumMomSmear].GetLineWidth()*2.5);
    gCkRatio_pp[NumMomSmear-1].SetLineWidth(gCkRatio_pp[NumMomSmear-1].GetLineWidth()*2.5);


    TCanvas* cRatio_pL = new TCanvas("cRatio_pL", "cRatio_pL", 1);
    cRatio_pL->cd(0); cRatio_pL->SetCanvasSize(1920, 1080); cRatio_pL->SetMargin(0.15,0.05,0.2,0.05);//lrbt
    cRatio_pL->SetGrid(true);
    TLegend* lLegend_pL = new TLegend(0.5,0.75,0.95,0.95);//lbrt
    lLegend_pL->SetName(TString::Format("lLegend"));
    lLegend_pL->SetTextSize(0.045);
    lLegend_pL->AddEntry(&gCkRatio_pL[NumMomSmear],"Effect of resonance width");
    lLegend_pL->AddEntry(&gCkRatio_pL[NumMomSmear-1],"Momentum smearing of 20%");
    TPaveText* PT_pL = new TPaveText(0.2,0.85,0.5,0.95, "blNDC");//lbrt
    PT_pL->SetName("PT_pL");
    PT_pL->SetBorderSize(1);
    PT_pL->SetTextSize(0.045);
    PT_pL->SetFillColor(kWhite);
    PT_pL->SetTextFont(22);
    PT_pL->AddText("p-#Lambda correlation");

    hAxis->Draw("axis");
    gCkRatio_pL[NumMomSmear].SetLineWidth(6);
    gCkRatio_pL[NumMomSmear].SetLineColor(kBlue);
    gCkRatio_pL[NumMomSmear-1].SetLineWidth(6);
    gCkRatio_pL[NumMomSmear-1].SetLineColor(kRed);

    gCkRatio_pL[NumMomSmear].Draw("same,c");//only mass smearing
    gCkRatio_pL[NumMomSmear-1].Draw("same,c");//max momentum smear
    lLegend_pL->Draw("same");
    PT_pL->Draw("same");

    cRatio_pL->SaveAs(TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/SourceSmearing/cRatio_pL.png"));

	gStyle->SetLineWidth(1);
	gCkRatio_pL[NumMomSmear].SetLineWidth(gCkRatio_pL[NumMomSmear].GetLineWidth()/2.5);
    gCkRatio_pL[NumMomSmear-1].SetLineWidth(gCkRatio_pL[NumMomSmear-1].GetLineWidth()/2.5);
    cRatio_pL->SaveAs(TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/SourceSmearing/cRatio_pL.pdf"));
	gStyle->SetLineWidth(2.5);
	gCkRatio_pL[NumMomSmear].SetLineWidth(gCkRatio_pL[NumMomSmear].GetLineWidth()*2.5);
    gCkRatio_pL[NumMomSmear-1].SetLineWidth(gCkRatio_pL[NumMomSmear-1].GetLineWidth()*2.5);


    TCanvas* cRatio_pXim = new TCanvas("cRatio_pXim", "cRatio_pXim", 1);
    cRatio_pXim->cd(0); cRatio_pXim->SetCanvasSize(1920, 1080); cRatio_pXim->SetMargin(0.15,0.05,0.2,0.05);//lrbt
    cRatio_pXim->SetGrid(true);
    TLegend* lLegend_pXim = new TLegend(0.5,0.75,0.95,0.95);//lbrt
    lLegend_pXim->SetName(TString::Format("lLegend"));
    lLegend_pXim->SetTextSize(0.045);
    lLegend_pXim->AddEntry(&gCkRatio_pXim[NumMomSmear],"Effect of resonance width");
    lLegend_pXim->AddEntry(&gCkRatio_pXim[NumMomSmear-1],"Momentum smearing of 20%");
    TPaveText* PT_pXim = new TPaveText(0.2,0.85,0.5,0.95, "blNDC");//lbrt
    PT_pXim->SetName("PT_pXim");
    PT_pXim->SetBorderSize(1);
    PT_pXim->SetTextSize(0.045);
    PT_pXim->SetFillColor(kWhite);
    PT_pXim->SetTextFont(22);
    PT_pXim->AddText("p-#Xi^{#minus} correlation");

    hAxis->Draw("axis");
    gCkRatio_pXim[NumMomSmear].SetLineWidth(6);
    gCkRatio_pXim[NumMomSmear].SetLineColor(kBlue);
    gCkRatio_pXim[NumMomSmear-1].SetLineWidth(6);
    gCkRatio_pXim[NumMomSmear-1].SetLineColor(kRed);

    gCkRatio_pXim[NumMomSmear].Draw("same,c");//only mass smearing
    gCkRatio_pXim[NumMomSmear-1].Draw("same,c");//max momentum smear
    lLegend_pXim->Draw("same");
    PT_pXim->Draw("same");

    cRatio_pXim->SaveAs(TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/SourceSmearing/cRatio_pXim.png"));

	gStyle->SetLineWidth(1);
	gCkRatio_pXim[NumMomSmear].SetLineWidth(gCkRatio_pXim[NumMomSmear].GetLineWidth()/2.5);
    gCkRatio_pXim[NumMomSmear-1].SetLineWidth(gCkRatio_pXim[NumMomSmear-1].GetLineWidth()/2.5);
    cRatio_pXim->SaveAs(TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/SourceSmearing/cRatio_pXim.pdf"));
	gStyle->SetLineWidth(2.5);
	gCkRatio_pXim[NumMomSmear].SetLineWidth(gCkRatio_pXim[NumMomSmear].GetLineWidth()*2.5);
    gCkRatio_pXim[NumMomSmear-1].SetLineWidth(gCkRatio_pXim[NumMomSmear-1].GetLineWidth()*2.5);


    TCanvas* cRatio_pOmegam = new TCanvas("cRatio_pOmegam", "cRatio_pOmegam", 1);
    cRatio_pOmegam->cd(0); cRatio_pOmegam->SetCanvasSize(1920, 1080); cRatio_pOmegam->SetMargin(0.15,0.05,0.2,0.05);//lrbt
    cRatio_pOmegam->SetGrid(true);
    TLegend* lLegend_pOmegam = new TLegend(0.5,0.75,0.95,0.95);//lbrt
    lLegend_pOmegam->SetName(TString::Format("lLegend"));
    lLegend_pOmegam->SetTextSize(0.045);
    lLegend_pOmegam->AddEntry(&gCkRatio_pOmegam[NumMomSmear],"Effect of resonance width");
    lLegend_pOmegam->AddEntry(&gCkRatio_pOmegam[NumMomSmear-1],"Momentum smearing of 20%");
    TPaveText* PT_pOmegam = new TPaveText(0.2,0.85,0.5,0.95, "blNDC");//lbrt
    PT_pOmegam->SetName("PT_pOmegam");
    PT_pOmegam->SetBorderSize(1);
    PT_pOmegam->SetTextSize(0.045);
    PT_pOmegam->SetFillColor(kWhite);
    PT_pOmegam->SetTextFont(22);
    PT_pOmegam->AddText("p-#Omega^{#minus} correlation");

    hAxis->Draw("axis");
    gCkRatio_pOmegam[NumMomSmear].SetLineWidth(6);
    gCkRatio_pOmegam[NumMomSmear].SetLineColor(kBlue);
    gCkRatio_pOmegam[NumMomSmear-1].SetLineWidth(6);
    gCkRatio_pOmegam[NumMomSmear-1].SetLineColor(kRed);

    gCkRatio_pOmegam[NumMomSmear].Draw("same,c");//only mass smearing
    gCkRatio_pOmegam[NumMomSmear-1].Draw("same,c");//max momentum smear
    lLegend_pOmegam->Draw("same");
    PT_pOmegam->Draw("same");

    cRatio_pOmegam->SaveAs(TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/SourceSmearing/cRatio_pOmegam.png"));

	gStyle->SetLineWidth(1);
	gCkRatio_pOmegam[NumMomSmear].SetLineWidth(gCkRatio_pOmegam[NumMomSmear].GetLineWidth()/2.5);
    gCkRatio_pOmegam[NumMomSmear-1].SetLineWidth(gCkRatio_pOmegam[NumMomSmear-1].GetLineWidth()/2.5);
    cRatio_pOmegam->SaveAs(TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/SourceSmearing/cRatio_pOmegam.pdf"));
	gStyle->SetLineWidth(2.5);
	gCkRatio_pOmegam[NumMomSmear].SetLineWidth(gCkRatio_pOmegam[NumMomSmear].GetLineWidth()*2.5);
    gCkRatio_pOmegam[NumMomSmear-1].SetLineWidth(gCkRatio_pOmegam[NumMomSmear-1].GetLineWidth()*2.5);


    TH1F* hAxisSource = new TH1F("hAxisSource", "hAxisSource", 256, 0, 8);
    hAxisSource->SetStats(false);
    hAxisSource->SetTitle("");
    hAxisSource->GetXaxis()->SetTitle("r (fm)");
    hAxisSource->GetXaxis()->SetTitleSize(0.06);
    hAxisSource->GetXaxis()->SetLabelSize(0.06);
    hAxisSource->GetXaxis()->CenterTitle();
    hAxisSource->GetXaxis()->SetTitleOffset(1.3);
    hAxisSource->GetXaxis()->SetLabelOffset(0.02);
    hAxisSource->GetYaxis()->SetTitle("4#pir^{2}S(r) (1/fm)");
    hAxisSource->GetYaxis()->SetTitleSize(0.06);
    hAxisSource->GetYaxis()->SetLabelSize(0.06);
    hAxisSource->GetYaxis()->CenterTitle();
    hAxisSource->GetYaxis()->SetTitleOffset(1.00);
    hAxisSource->GetYaxis()->SetRangeUser(0.0,0.60);
    //hAxisSource->GetYaxis()->SetNdivisions(404);

    TCanvas* cSource_pp = new TCanvas("cSource_pp", "cSource_pp", 1);
    cSource_pp->cd(0); cSource_pp->SetCanvasSize(1920, 1080); cSource_pp->SetMargin(0.15,0.05,0.2,0.05);//lrbt
    cSource_pp->SetGrid(true);
    TLegend* sLegend_pp = new TLegend(0.5,0.75,0.95,0.95);//lbrt
    sLegend_pp->SetName(TString::Format("lLegend"));
    sLegend_pp->SetTextSize(0.045);
    sLegend_pp->AddEntry(&gSourceGauss_pp,"Gaussian source");
    sLegend_pp->AddEntry(&gSourceGaussReso_pp,"Gaussian core + resonances");
    TPaveText* PTs_pp = new TPaveText(0.2,0.85,0.5,0.95, "blNDC");//lbrt
    PTs_pp->SetName("PTs_pp");
    PTs_pp->SetBorderSize(1);
    PTs_pp->SetTextSize(0.045);
    PTs_pp->SetFillColor(kWhite);
    PTs_pp->SetTextFont(22);
    PTs_pp->AddText("p-p source");

    hAxisSource->Draw("axis");
    gSourceGauss_pp.SetLineWidth(6);
    gSourceGauss_pp.SetLineColor(kBlue);
    gSourceGaussReso_pp.SetLineWidth(6);
    gSourceGaussReso_pp.SetLineColor(kRed);

    gSourceGauss_pp.Draw("same,c");//only mass smearing
    gSourceGaussReso_pp.Draw("same,c");//max momentum smear
    sLegend_pp->Draw("same");
    PTs_pp->Draw("same");

    gStyle->SetLineWidth(2.5);
    cSource_pp->SaveAs(TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/SourceSmearing/cSource_pp.png"));

	gStyle->SetLineWidth(1);
	gSourceGauss_pp.SetLineWidth(gSourceGauss_pp.GetLineWidth()/2.5);
    gSourceGaussReso_pp.SetLineWidth(gSourceGaussReso_pp.GetLineWidth()/2.5);
    cSource_pp->SaveAs(TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/SourceSmearing/cSource_pp.pdf"));
	gStyle->SetLineWidth(2.5);
	gSourceGauss_pp.SetLineWidth(gSourceGauss_pp.GetLineWidth()*2.5);
    gSourceGaussReso_pp.SetLineWidth(gSourceGaussReso_pp.GetLineWidth()*2.5);


    TCanvas* cSource_pL = new TCanvas("cSource_pL", "cSource_pL", 1);
    cSource_pL->cd(0); cSource_pL->SetCanvasSize(1920, 1080); cSource_pL->SetMargin(0.15,0.05,0.2,0.05);//lrbt
    cSource_pL->SetGrid(true);
    TLegend* sLegend_pL = new TLegend(0.5,0.75,0.95,0.95);//lbrt
    sLegend_pL->SetName(TString::Format("lLegend"));
    sLegend_pL->SetTextSize(0.045);
    sLegend_pL->AddEntry(&gSourceGauss_pL,"Gaussian source");
    sLegend_pL->AddEntry(&gSourceGaussReso_pL,"Gaussian core + resonances");
    TPaveText* PTs_pL = new TPaveText(0.2,0.85,0.5,0.95, "blNDC");//lbrt
    PTs_pL->SetName("PTs_pL");
    PTs_pL->SetBorderSize(1);
    PTs_pL->SetTextSize(0.045);
    PTs_pL->SetFillColor(kWhite);
    PTs_pL->SetTextFont(22);
    PTs_pL->AddText("p-#Lambda source");

    hAxisSource->Draw("axis");
    gSourceGauss_pL.SetLineWidth(6);
    gSourceGauss_pL.SetLineColor(kBlue);
    gSourceGaussReso_pL.SetLineWidth(6);
    gSourceGaussReso_pL.SetLineColor(kRed);

    gSourceGauss_pL.Draw("same,c");//only mass smearing
    gSourceGaussReso_pL.Draw("same,c");//max momentum smear
    sLegend_pL->Draw("same");
    PTs_pL->Draw("same");

    gStyle->SetLineWidth(2.5);
    cSource_pL->SaveAs(TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/SourceSmearing/cSource_pL.png"));

	gStyle->SetLineWidth(1);
	gSourceGauss_pL.SetLineWidth(gSourceGauss_pL.GetLineWidth()/2.5);
    gSourceGaussReso_pL.SetLineWidth(gSourceGaussReso_pL.GetLineWidth()/2.5);
    cSource_pL->SaveAs(TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/SourceSmearing/cSource_pL.pdf"));
	gStyle->SetLineWidth(2.5);
	gSourceGauss_pL.SetLineWidth(gSourceGauss_pL.GetLineWidth()*2.5);
    gSourceGaussReso_pL.SetLineWidth(gSourceGaussReso_pL.GetLineWidth()*2.5);


    TCanvas* cSource_pXim = new TCanvas("cSource_pXim", "cSource_pXim", 1);
    cSource_pXim->cd(0); cSource_pXim->SetCanvasSize(1920, 1080); cSource_pXim->SetMargin(0.15,0.05,0.2,0.05);//lrbt
    cSource_pXim->SetGrid(true);
    TLegend* sLegend_pXim = new TLegend(0.5,0.75,0.95,0.95);//lbrt
    sLegend_pXim->SetName(TString::Format("lLegend"));
    sLegend_pXim->SetTextSize(0.045);
    sLegend_pXim->AddEntry(&gSourceGauss_pXim,"Gaussian source");
    sLegend_pXim->AddEntry(&gSourceGaussReso_pXim,"Gaussian core + resonances");
    TPaveText* PTs_pXim = new TPaveText(0.2,0.85,0.5,0.95, "blNDC");//lbrt
    PTs_pXim->SetName("PTs_pXim");
    PTs_pXim->SetBorderSize(1);
    PTs_pXim->SetTextSize(0.045);
    PTs_pXim->SetFillColor(kWhite);
    PTs_pXim->SetTextFont(22);
    PTs_pXim->AddText("p-#Xi^{#minus} source");

    hAxisSource->Draw("axis");
    gSourceGauss_pXim.SetLineWidth(6);
    gSourceGauss_pXim.SetLineColor(kBlue);
    gSourceGaussReso_pXim.SetLineWidth(6);
    gSourceGaussReso_pXim.SetLineColor(kRed);

    gSourceGauss_pXim.Draw("same,c");//only mass smearing
    gSourceGaussReso_pXim.Draw("same,c");//max momentum smear
    sLegend_pXim->Draw("same");
    PTs_pXim->Draw("same");

    gStyle->SetLineWidth(2.5);
    cSource_pXim->SaveAs(TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/SourceSmearing/cSource_pXim.png"));

	gStyle->SetLineWidth(1);
	gSourceGauss_pXim.SetLineWidth(gSourceGauss_pXim.GetLineWidth()/2.5);
    gSourceGaussReso_pXim.SetLineWidth(gSourceGaussReso_pXim.GetLineWidth()/2.5);
    cSource_pXim->SaveAs(TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/SourceSmearing/cSource_pXim.pdf"));
	gStyle->SetLineWidth(2.5);
	gSourceGauss_pXim.SetLineWidth(gSourceGauss_pXim.GetLineWidth()*2.5);
    gSourceGaussReso_pXim.SetLineWidth(gSourceGaussReso_pXim.GetLineWidth()*2.5);


    TCanvas* cSource_pOmegam = new TCanvas("cSource_pOmegam", "cSource_pOmegam", 1);
    cSource_pOmegam->cd(0); cSource_pOmegam->SetCanvasSize(1920, 1080); cSource_pOmegam->SetMargin(0.15,0.05,0.2,0.05);//lrbt
    cSource_pOmegam->SetGrid(true);
    TLegend* sLegend_pOmegam = new TLegend(0.5,0.75,0.95,0.95);//lbrt
    sLegend_pOmegam->SetName(TString::Format("lLegend"));
    sLegend_pOmegam->SetTextSize(0.045);
    sLegend_pOmegam->AddEntry(&gSourceGauss_pOmegam,"Gaussian source");
    sLegend_pOmegam->AddEntry(&gSourceGaussReso_pOmegam,"Gaussian core + resonances");
    TPaveText* PTs_pOmegam = new TPaveText(0.2,0.85,0.5,0.95, "blNDC");//lbrt
    PTs_pOmegam->SetName("PTs_pOmegam");
    PTs_pOmegam->SetBorderSize(1);
    PTs_pOmegam->SetTextSize(0.045);
    PTs_pOmegam->SetFillColor(kWhite);
    PTs_pOmegam->SetTextFont(22);
    PTs_pOmegam->AddText("p-#Omega^{#minus} source");

    hAxisSource->Draw("axis");
    gSourceGauss_pOmegam.SetLineWidth(6);
    gSourceGauss_pOmegam.SetLineColor(kBlue);
    gSourceGaussReso_pOmegam.SetLineWidth(6);
    gSourceGaussReso_pOmegam.SetLineColor(kRed);

    gSourceGauss_pOmegam.Draw("same,c");//only mass smearing
    gSourceGaussReso_pOmegam.Draw("same,c");//max momentum smear
    sLegend_pOmegam->Draw("same");
    PTs_pOmegam->Draw("same");

    gStyle->SetLineWidth(2.5);
    cSource_pOmegam->SaveAs(TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/SourceSmearing/cSource_pOmegam.png"));

	gStyle->SetLineWidth(1);
	gSourceGauss_pOmegam.SetLineWidth(gSourceGauss_pOmegam.GetLineWidth()/2.5);
    gSourceGaussReso_pOmegam.SetLineWidth(gSourceGaussReso_pOmegam.GetLineWidth()/2.5);
    cSource_pOmegam->SaveAs(TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/SourceSmearing/cSource_pOmegam.pdf"));
	gStyle->SetLineWidth(2.5);
	gSourceGauss_pOmegam.SetLineWidth(gSourceGauss_pOmegam.GetLineWidth()*2.5);
    gSourceGaussReso_pOmegam.SetLineWidth(gSourceGaussReso_pOmegam.GetLineWidth()*2.5);


    gpp_low.Write();
    gpp_up.Write();

    gpL_low.Write();
    gpL_up.Write();

    gpXim_low.Write();
    gpXim_up.Write();

    gpOmegam_low.Write();
    gpOmegam_up.Write();

    gSourceGauss_pp.Write();
    gSourceGaussReso_pp.Write();

    gSourceGauss_pL.Write();
    gSourceGaussReso_pL.Write();

    gSourceGauss_pXim.Write();
    gSourceGaussReso_pXim.Write();

    gSourceGauss_pOmegam.Write();
    gSourceGaussReso_pOmegam.Write();

    delete [] gCk_pp;
    delete [] gCk_pL;
    delete [] gCk_pXim;
    delete [] gCk_pOmegam;

    delete [] gCkRatio_pp;
    delete [] gCkRatio_pL;
    delete [] gCkRatio_pXim;
    delete [] gCkRatio_pOmegam;

    delete lLegend_pp;
    delete PT_pp;
    delete cRatio_pp;

    delete lLegend_pL;
    delete PT_pL;
    delete cRatio_pL;

    delete lLegend_pXim;
    delete PT_pXim;
    delete cRatio_pXim;

    delete lLegend_pOmegam;
    delete PT_pOmegam;
    delete cRatio_pOmegam;

    delete hAxis;

    delete sLegend_pp;
    delete sLegend_pL;
    delete sLegend_pXim;
    delete sLegend_pOmegam;
    delete PTs_pp;
    delete PTs_pL;
    delete PTs_pXim;
    delete PTs_pOmegam;
    delete hAxisSource;

    delete cSource_pp;
    delete cSource_pL;
    delete cSource_pXim;
    delete cSource_pOmegam;

    delete fOutputFile;
}


void Test_pOmega_Potentials(){

    unsigned NumPots = 6;
    int PotFlag[NumPots] = {11,12,13,121,122,123};

    //unsigned NumPots = 3;
    //int PotFlag[NumPots] = {121,122,123};

    TGraph* gCk = new TGraph [NumPots];
    TGraph* gCkSI = new TGraph [NumPots];
    TGraph* gPot3S1 = new TGraph [NumPots];
    TGraph* gPot5S2 = new TGraph [NumPots];

    const double kMin = 0;
    const double kMax = 280;
    const unsigned NumMomBins = 56;

    const double rMin = 0.01;
    const double rMax = 1.6;
    const unsigned NumRadBins = 128;
    const double rStep = (rMax-rMin)/double(NumRadBins);

    DLM_CommonAnaFunctions AnaObject;

    for(unsigned uPot=0; uPot<NumPots; uPot++){
        gCk[uPot].SetName(TString::Format("gCk_%i",PotFlag[uPot]));
        gCkSI[uPot].SetName(TString::Format("gCkSI_%i",PotFlag[uPot]));
        gPot3S1[uPot].SetName(TString::Format("gPot3S1_%i",PotFlag[uPot]));
        gPot5S2[uPot].SetName(TString::Format("gPot5S2_%i",PotFlag[uPot]));
        CATS AB_pOmegam;
        AB_pOmegam.SetMomBins(NumMomBins,kMin,kMax);
        AnaObject.SetUpCats_pOmegam(AB_pOmegam,TString::Format("pOmega_Lattice_%i",PotFlag[uPot]),"Gauss");
        AB_pOmegam.SetAnaSource(0,2.5);
        AB_pOmegam.SetEpsilonConv(1e-7);
        AB_pOmegam.SetEpsilonProp(1e-7);
        AB_pOmegam.SetShortRangePotential(0,0,0,1e5);
        AB_pOmegam.SetShortRangePotential(0,0,1,2.0);
        AB_pOmegam.SetShortRangePotential(0,0,2,0.000001);
        AB_pOmegam.SetMaxNumThreads(4);
        //AB_pOmegam.SetChannelWeight(0,0);
        //AB_pOmegam.SetChannelWeight(1,1);
        AB_pOmegam.KillTheCat();

        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            gCk[uPot].SetPoint(uMomBin,AB_pOmegam.GetMomentum(uMomBin),AB_pOmegam.GetCorrFun(uMomBin));
        }
        for(unsigned uRadBin=0; uRadBin<NumRadBins; uRadBin++){
            double RADIUS = rMin+double(uRadBin)*rStep;
            gPot3S1[uPot].SetPoint(uRadBin,RADIUS,AB_pOmegam.EvaluateThePotential(0,0,0,RADIUS));
            gPot5S2[uPot].SetPoint(uRadBin,RADIUS,AB_pOmegam.EvaluateThePotential(1,0,0,RADIUS));
        }

        AB_pOmegam.SetQ1Q2(0);
        AB_pOmegam.KillTheCat();
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            gCkSI[uPot].SetPoint(uMomBin,AB_pOmegam.GetMomentum(uMomBin),AB_pOmegam.GetCorrFun(uMomBin));
        }

    }

    TFile* fOutputFile = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/Test_pOmega_Potentials/fOutputFile.root", "recreate");
    for(unsigned uPot=0; uPot<NumPots; uPot++){
        gCk[uPot].Write();
    }
    for(unsigned uPot=0; uPot<NumPots; uPot++){
        gCkSI[uPot].Write();
    }
    for(unsigned uPot=0; uPot<NumPots; uPot++){
        gPot3S1[uPot].Write();
    }
    for(unsigned uPot=0; uPot<NumPots; uPot++){
        gPot5S2[uPot].Write();
    }

    delete [] gCk;
    delete [] gCkSI;
    delete [] gPot3S1;
    delete [] gPot5S2;
    delete fOutputFile;
}

void Test_pp_Dwaves(){

    const unsigned NumDataPts = 11;
    const double DataPoints1S0[NumDataPts] = {32.68,54.74,55.09,48.51,38.78,25.01,15,6.99,0.23,-5.64,-10.86};//1S0
    //const double DataPoints3S1_pn[NumDataPts] = {147.747,118.178,102.611,80.63,62.77,43.23,30.72,21.22,13.39,6.60,0.502};//3S1 for pn
    const double DataPoints3P0[NumDataPts] = {0.134,1.582,3.729,8.575,11.47,9.45,4.74,-0.37,-5.43,-10.39,-15.30};//3P0
    const double DataPoints3P1[NumDataPts] = {-0.081,-0.902,-2.060,-4.932,-8.317,-13.258,-17.434,-21.25,-24.77,-27.99,-30.89};//3P1
    const double DataPoints3P2[NumDataPts] = {0.014,0.214,0.651,2.491,5.855,11.013,13.982,15.63,16.59,17.17,17.54};//3P2
    const double DataPoints1D2[NumDataPts] = {0.001,0.043,0.165,0.696,1.711,3.790,5.606,7.058,8.27,9.42,10.69};//1D2
    //const double DataPoints3F4[NumDataPts] = {0,0,0.001,0.02,0.108,0.478,1.032,1.678,2.325,2.89,3.3};//3F4
    const double DataPtsT[NumDataPts] = {1,5,10,25,50,100,150,200,250,300,350};
    double DataPtsK[NumDataPts];
    for(unsigned uPts=0; uPts<NumDataPts; uPts++){
        DataPtsK[uPts] = tLab_pCm(DataPtsT[uPts],Mass_p,Mass_p);
    }

    const unsigned NumFineBins = 50;
    const double FineStep = 4;
    const unsigned NumCoarseBins = 20;
    const double CoarsStep = 16;
    const unsigned NumMomBins = NumFineBins+NumCoarseBins;
    double* MomBins = new double [NumMomBins+1];
    MomBins[0] = 0;
    for(unsigned uBin=1; uBin<=NumMomBins; uBin++){
        if(uBin<=NumFineBins){
            MomBins[uBin] = MomBins[uBin-1]+FineStep;
        }
        else{
            MomBins[uBin] = MomBins[uBin-1]+CoarsStep;
        }
        //printf("%u = %f\n",uBin,MomBins[uBin]);
    }

    const TString SOURCE = "McGauss_Reso";

    CATSparameters cPars(CATSparameters::tSource,1,true);
    CATSparameters cPotPars1S0(CATSparameters::tPotential,8,true);
    CATSparameters cPotPars3P0(CATSparameters::tPotential,8,true);
    CATSparameters cPotPars3P1(CATSparameters::tPotential,8,true);
    CATSparameters cPotPars3P2(CATSparameters::tPotential,8,true);
    CATSparameters cPotPars1D2(CATSparameters::tPotential,8,true);

    DLM_CleverMcLevyReso CleverMcLevyReso;

    CATS Kitty;

    Kitty.SetMomBins(NumMomBins,MomBins);
    Kitty.SetEpsilonConv(1e-7);
    Kitty.SetEpsilonProp(1e-7);

    if(SOURCE=="Gauss"){
        cPars.SetParameter(0,1.2);
        Kitty.SetAnaSource(GaussSource, cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="McLevyNolan_Reso"){
        CleverMcLevyReso.InitStability(21,1,2);
        CleverMcLevyReso.InitScale(38,0.15,2.0);
        CleverMcLevyReso.InitRad(257,0,64);
        CleverMcLevyReso.InitType(2);
        CleverMcLevyReso.InitReso(0,1);
        CleverMcLevyReso.InitReso(1,1);
        CleverMcLevyReso.SetUpReso(0,0,1.-0.3578,1361.52,1.65,Mass_p,Mass_pic);
        CleverMcLevyReso.SetUpReso(1,0,1.-0.3578,1361.52,1.65,Mass_p,Mass_pic);
        CleverMcLevyReso.InitNumMcIter(1000000);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverMcLevyReso, 2);
        Kitty.SetAnaSource(0,1.0);
        Kitty.SetAnaSource(1,1.6);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="McGauss_Reso"){
        CleverMcLevyReso.InitStability(1,2-1e-6,2+1e-6);
        CleverMcLevyReso.InitScale(38,0.15,2.0);
        CleverMcLevyReso.InitRad(257,0,64);
        CleverMcLevyReso.InitType(2);
        CleverMcLevyReso.InitReso(0,1);
        CleverMcLevyReso.InitReso(1,1);
        CleverMcLevyReso.SetUpReso(0,0,1.-0.3578,1361.52,1.65,Mass_p,Mass_pic);
        CleverMcLevyReso.SetUpReso(1,0,1.-0.3578,1361.52,1.65,Mass_p,Mass_pic);
        CleverMcLevyReso.InitNumMcIter(1000000);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverMcLevyReso, 2);
        Kitty.SetAnaSource(0,0.995);
        Kitty.SetAnaSource(1,2.0);
        Kitty.SetUseAnalyticSource(true);
    }
    else{
        printf("\033[1;31mERROR:\033[0m Non-existing source '%s'\n",SOURCE.Data());
        return;
    }


    //#,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
    double PotPars1S0[8]={NN_AV18,v18_Coupled3P2,1,1,1,0,0,0};
    double PotPars3P0[8]={NN_AV18,v18_Coupled3P2,1,1,1,1,1,0};
    double PotPars3P1[8]={NN_AV18,v18_Coupled3P2,1,1,1,1,1,1};
    double PotPars3P2[8]={NN_AV18,v18_SingleChannelMagic,1,1,1,1,1,2};
    double PotPars1D2[8]={NN_AV18,v18_Coupled3P2,1,1,1,0,2,2};
    cPotPars1S0.SetParameters(PotPars1S0);
    cPotPars3P0.SetParameters(PotPars3P0);
    cPotPars3P1.SetParameters(PotPars3P1);
    cPotPars3P2.SetParameters(PotPars3P2);
    cPotPars1D2.SetParameters(PotPars1D2);

    Kitty.SetMomentumDependentSource(false);
    Kitty.SetThetaDependentSource(false);
    Kitty.SetExcludeFailedBins(false);

    Kitty.SetQ1Q2(1);
    Kitty.SetPdgId(2212, 2212);
    Kitty.SetRedMass( 0.5*Mass_p );

    Kitty.SetNumChannels(4);
    Kitty.SetNumPW(0,3);
    Kitty.SetNumPW(1,2);
    Kitty.SetNumPW(2,2);
    Kitty.SetNumPW(3,2);
    Kitty.SetSpin(0,0);
    Kitty.SetSpin(1,1);
    Kitty.SetSpin(2,1);
    Kitty.SetSpin(3,1);
    Kitty.SetChannelWeight(0, 3./12.);
    Kitty.SetChannelWeight(1, 1./12.);
    Kitty.SetChannelWeight(2, 3./12.);
    Kitty.SetChannelWeight(3, 5./12.);

    Kitty.SetShortRangePotential(0,0,fDlmPot,cPotPars1S0);
    Kitty.SetShortRangePotential(1,1,fDlmPot,cPotPars3P0);
    Kitty.SetShortRangePotential(2,1,fDlmPot,cPotPars3P1);
    Kitty.SetShortRangePotential(3,1,fDlmPot,cPotPars3P2);

    TGraph gSP;
    gSP.SetName("gSP");
    TGraph gSPD;
    gSPD.SetName("gSPD");
    TGraph gSPvsSPD;
    gSPvsSPD.SetName("gSPvsSPD");

    TGraph gPS_1S0;
    TGraph gPS_3P0;
    TGraph gPS_3P1;
    TGraph gPS_3P2;
    TGraph gPS_1D2;
    gPS_1S0.SetName("gPS_1S0");
    gPS_3P0.SetName("gPS_3P0");
    gPS_3P1.SetName("gPS_3P1");
    gPS_3P2.SetName("gPS_3P2");
    gPS_1D2.SetName("gPS_1D2");

    TGraph gPSth_1S0;
    TGraph gPSth_3P0;
    TGraph gPSth_3P1;
    TGraph gPSth_3P2;
    TGraph gPSth_1D2;
    gPSth_1S0.SetName("gPSth_1S0");
    gPSth_3P0.SetName("gPSth_3P0");
    gPSth_3P1.SetName("gPSth_3P1");
    gPSth_3P2.SetName("gPSth_3P2");
    gPSth_1D2.SetName("gPSth_1D2");

    for(unsigned uPts=0; uPts<NumDataPts; uPts++){
        gPSth_1S0.SetPoint(uPts,DataPtsK[uPts],DataPoints1S0[uPts]);
        gPSth_3P0.SetPoint(uPts,DataPtsK[uPts],DataPoints3P0[uPts]);
        gPSth_3P1.SetPoint(uPts,DataPtsK[uPts],DataPoints3P1[uPts]);
        gPSth_3P2.SetPoint(uPts,DataPtsK[uPts],DataPoints3P2[uPts]);
        gPSth_1D2.SetPoint(uPts,DataPtsK[uPts],DataPoints1D2[uPts]);
    }


    Kitty.KillTheCat();
    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        //printf("C(%.2f) = %.2f\n",Kitty.GetMomentum(uBin),Kitty.GetCorrFun(uBin));
        gSP.SetPoint(uBin,Kitty.GetMomentum(uBin),Kitty.GetCorrFun(uBin));
    }

    Kitty.SetShortRangePotential(0,2,fDlmPot,cPotPars1D2);
    Kitty.KillTheCat();
    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        //printf("C(%.2f) = %.2f\n",Kitty.GetMomentum(uBin),Kitty.GetCorrFun(uBin));
        gSPD.SetPoint(uBin,Kitty.GetMomentum(uBin),Kitty.GetCorrFun(uBin));
    }


    TF1* fPhaseShift1S0 = new TF1("fPhaseShift1S0", "[2]*x*x*x*x+0.5*[1]*x*x-1./[0]", 2, 60);
    TF1* fPhaseShift3P0 = new TF1("fPhaseShift3P0", "[2]*x*x*x*x+0.5*[1]*x*x-1./[0]", 2, 60);
    TF1* fPhaseShift3P1 = new TF1("fPhaseShift3P1", "[2]*x*x*x*x+0.5*[1]*x*x-1./[0]", 2, 60);
    TF1* fPhaseShift3P2 = new TF1("fPhaseShift3P2", "[2]*x*x*x*x+0.5*[1]*x*x-1./[0]", 2, 60);
    TF1* fPhaseShift1D2 = new TF1("fPhaseShift1D2", "[2]*x*x*x*x+0.5*[1]*x*x-1./[0]", 2, 60);

    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        gPS_1S0.SetPoint(uBin,Kitty.GetMomentum(uBin),Kitty.GetPhaseShift(uBin,0,0)*180./3.1415);
        gPS_3P0.SetPoint(uBin,Kitty.GetMomentum(uBin),Kitty.GetPhaseShift(uBin,1,1)*180./3.1415);
        gPS_3P1.SetPoint(uBin,Kitty.GetMomentum(uBin),Kitty.GetPhaseShift(uBin,2,1)*180./3.1415);
        gPS_3P2.SetPoint(uBin,Kitty.GetMomentum(uBin),Kitty.GetPhaseShift(uBin,3,1)*180./3.1415);
        gPS_1D2.SetPoint(uBin,Kitty.GetMomentum(uBin),Kitty.GetPhaseShift(uBin,0,2)*180./3.1415);
    }

/*
    for(unsigned iBin=0; iBin<NumFineSteps; iBin++){
        hPhaseShift0->SetBinContent(iBin+1, NalaFine.GetMomentum(iBin)/tan(NalaFine.GetPhaseShift(iBin,0,0)));
        //printf("k = %.2f; delta = %.3f; whatever/tan=%.3f\n",NalaFine.GetMomentum(iBin),NalaFine.GetPhaseShift(iBin,0,0)*180./3.14,
        //       NalaFine.GetMomentum(iBin)/tan(NalaFine.GetPhaseShift(iBin,0,0)));
        hPhaseShift0->SetBinError(iBin+1, hPhaseShift0->GetBinContent(iBin+1)*0.0001);
    }

    TF1* fPhaseShift0 = new TF1("fPhaseShift0", "[2]*x*x*x*x+0.5*[1]*x*x-1./[0]", kMinFine, kMaxFine);
*/




    double xVal_SP, yVal_SP;
    double xVal_SPD, yVal_SPD;

    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        gSP.GetPoint(uBin,xVal_SP,yVal_SP);
        gSPD.GetPoint(uBin,xVal_SPD,yVal_SPD);
        gSPvsSPD.SetPoint(uBin,Kitty.GetMomentum(uBin),yVal_SPD/yVal_SP);
    }

    TFile* fOut = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/Test_pp_Dwaves/fOut.root","recreate");

    gSP.Write();
    gSPD.Write();
    gSPvsSPD.Write();

    gPSth_1S0.Write();
    gPSth_3P0.Write();
    gPSth_3P1.Write();
    gPSth_3P2.Write();
    gPSth_1D2.Write();

    gPS_1S0.Write();
    gPS_3P0.Write();
    gPS_3P1.Write();
    gPS_3P2.Write();
    gPS_1D2.Write();


    delete fOut;

    delete [] MomBins;
}

void Test_pSigma0(){

    const unsigned NumMomBins = 70;
    const double kMin = 0;
    const double kMax = 350;
    DLM_Histo<complex<double>>*** ExternalWF=NULL;
    CATS Kitty;
    Kitty.SetMomBins(NumMomBins,kMin,kMax);

    CATSparameters cPars(CATSparameters::tSource,1,true);
    cPars.SetParameter(0,1.0);

    Kitty.SetUseAnalyticSource(true);
    Kitty.SetAnaSource(GaussSource,cPars);
    ExternalWF = Init_pS0_ESC16("/home/dmihaylov/CernBox/CATS_potentials/Tom/pSigma0/DimiValeNorm170519/",Kitty);
    Kitty.SetExternalWaveFunction(0,0,ExternalWF[0][0][0],ExternalWF[1][0][0]);
    Kitty.SetExternalWaveFunction(1,0,ExternalWF[0][1][0],ExternalWF[1][1][0]);

    Kitty.SetMomentumDependentSource(false);
    Kitty.SetThetaDependentSource(false);
    Kitty.SetExcludeFailedBins(false);
    Kitty.SetChannelWeight(0,0.25);
    Kitty.SetChannelWeight(1,0.75);
    Kitty.KillTheCat();

    TFile* OutputFile = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/Test_pSigma0/OutputFile.root","recreate");

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
        for(double RAD=0.05; RAD<20; RAD+=0.05){
            gWF1S0.SetPoint(COUNTER,RAD,abs(Kitty.EvalRadialWaveFunction(uBin,0,0,RAD,false)));
            gWF3S1.SetPoint(COUNTER,RAD,abs(Kitty.EvalRadialWaveFunction(uBin,1,0,RAD,false)));
            COUNTER++;
        }
        OutputFile->cd();
        gWF1S0.Write();
        gWF3S1.Write();
    }

    TGraph gKitty1S0;
    gKitty1S0.Set(Kitty.GetNumMomBins());
    gKitty1S0.SetName("gKitty1S0");
    Kitty.SetChannelWeight(0,1.0);
    Kitty.SetChannelWeight(1,0.0);
    Kitty.KillTheCat();
    for(unsigned uBin=0; uBin<Kitty.GetNumMomBins(); uBin++){
        gKitty1S0.SetPoint(uBin,Kitty.GetMomentum(uBin),Kitty.GetCorrFun(uBin));
    }

    TGraph gKitty3S1;
    gKitty3S1.Set(Kitty.GetNumMomBins());
    gKitty3S1.SetName("gKitty3S1");
    Kitty.SetChannelWeight(0,0.0);
    Kitty.SetChannelWeight(1,1.0);
    Kitty.KillTheCat();
    for(unsigned uBin=0; uBin<Kitty.GetNumMomBins(); uBin++){
        gKitty3S1.SetPoint(uBin,Kitty.GetMomentum(uBin),Kitty.GetCorrFun(uBin));
    }

    OutputFile->cd();
    gKitty.Write();
    gKitty1S0.Write();
    gKitty3S1.Write();

    delete OutputFile;
}

void TestVALE_pLambda(){

    DLM_CommonAnaFunctions DCAF;
    CATS TestCat;
    TestCat.SetMomBins(60,0,300);
    DCAF.SetUpCats_pL(TestCat,"LO","Gauss",0,0);

//return;
  unsigned NumMomBins = 60;
  double kMin = 0;
  double kMax = 300;
  const double MassProton = 938.272;
  const double MassLambda = 1115.683;
  double redmasspLambda = (MassProton*MassLambda)/(MassProton+MassLambda);

  const TString OutputFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Temp/TVPL.root";

  const TString POT = "LO";

  CATSparameters* cPars = NULL;
  CATSparameters* pPars = NULL;

  CATSparameters* cPotPars1S0 = NULL;
  CATSparameters* cPotPars3S1 = NULL;

  DLM_Histo<complex<double>>*** ExternalWF=NULL;
  unsigned NumChannels=0;


  double* RadBins;
  unsigned NumRadBins;
  const double MaxMomentum=300;

  complex<double> WaveFun[NumChannels];




  double rmin = 0.05;
  double rmax = 15;
  unsigned count = 200;
  double step = (rmax-rmin)/double(count);
  //labeling: WF I,S

   TGraph** grWFpL = new TGraph* [2];
   grWFpL[0] = new TGraph [NumMomBins];
   grWFpL[1] = new TGraph [NumMomBins];

  CATS pLambda;

  double SourceParameters[1]={1.};
  CATSparameters sParameters(CATSparameters::tSource,1,true);
  sParameters.SetParameters(SourceParameters);
  pLambda.SetAnaSource(GaussSource, sParameters);
  pLambda.SetUseAnalyticSource(true);

    pLambda.SetMomBins(NumMomBins,kMin,kMax);

    pLambda.SetMomentumDependentSource(false);
    pLambda.SetThetaDependentSource(false);
    pLambda.SetExcludeFailedBins(false);


  if(POT=="LO"){
    std::cout<<"In LO 1\n"<<std::endl;
        ExternalWF = Init_pL_Haidenbauer("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/Haidenbauer/pLambdaLO_600/",
                                  pLambda, 0, 600);
          NumChannels=2;
    std::cout<<"In LO 2\n"<<std::endl;
      }
      else if(POT=="LO_Coupled_S"){

     ExternalWF =  Init_pL_Haidenbauer("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/Haidenbauer/pLambdaLO_Coupling/",
                                  pLambda, 1, 600);
          NumChannels=4;
      }
      else if(POT=="NLO"){


     ExternalWF =  Init_pL_Haidenbauer("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/Haidenbauer/pLambdaNLO/",
                                  pLambda, 10, 600);
          NumChannels=2;
      }
      //s and p waves
      else if(POT=="NLO_sp"){
     ExternalWF =  Init_pL_Haidenbauer("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/Haidenbauer/pLambdaNLO/",
                                  pLambda, 12, 600);
          NumChannels=4;
      }
      else if(POT=="NLO_Coupled_S"){

     ExternalWF =   Init_pL_Haidenbauer("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/Haidenbauer/pLambdaNLO_Coupling/",
                                  pLambda, 11, 600);
          NumChannels=4;
      }
      else if(POT=="Usmani"){

          //#,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
          double PotPars1S0[8]={pL_UsmaniOli,0,0,0,0,0,0,0};
          double PotPars3S1[8]={pL_UsmaniOli,0,0,0,0,1,0,1};
          cPotPars1S0 = new CATSparameters(CATSparameters::tPotential,8,true); cPotPars1S0->SetParameters(PotPars1S0);
          cPotPars3S1 = new CATSparameters(CATSparameters::tPotential,8,true); cPotPars3S1->SetParameters(PotPars3S1);
          NumChannels=2;
      }

      //p-Omega system
      pLambda.SetQ1Q2(0);
      pLambda.SetPdgId(2212,3122);
      pLambda.SetRedMass(redmasspLambda);
      pLambda.SetNumChannels(NumChannels);
      pLambda.SetChannelWeight(0,0.25);
      pLambda.SetChannelWeight(1,0.75);
      pLambda.SetNumPW(0,1);
      pLambda.SetNumPW(1,1);
      pLambda.SetSpin(0, 0);
      pLambda.SetSpin(1, 1);
    pLambda.SetExternalWaveFunction(0,0,ExternalWF[0][0][0],ExternalWF[1][0][0]);
    pLambda.SetExternalWaveFunction(1,0,ExternalWF[0][1][0],ExternalWF[1][1][0]);
    pLambda.KillTheCat();

          //DCAF.SetUpCats_pL(pLambda,"LO","Gauss",0,0);

      TFile* outFilepLambda = new TFile(OutputFileName,"recreate");
      std::cout<<"In LO 4\n"<<std::endl;
      double Momentum;
      double Radius;


    for(unsigned uCh=0; uCh<NumChannels; uCh++){
    std::cout<<"In LO 5\n"<<std::endl;

      for(unsigned uMom=0; uMom<NumMomBins; uMom++){
        Momentum = pLambda.GetMomentum(uMom);
    std::cout<<"In LO 6\n"<<std::endl;
        grWFpL[uCh][uMom].SetName(TString::Format("WF_%.0f",Momentum));
        grWFpL[uCh][uMom].SetTitle(TString::Format("k^*=%.0f [Mev/c], channel=%c",Momentum,uCh));
        grWFpL[uCh][uMom].GetHistogram()->GetYaxis()->SetTitle("W.F. ");
        grWFpL[uCh][uMom].GetHistogram()->GetXaxis()->SetTitle("r (fm)");
        grWFpL[uCh][uMom].Set(count);
    std::cout<<"In LO 7\n"<<std::endl;
    //     if(!ExternalWF){
    //         pLambda.SetNumPW(uCh,1);
    //         pLambda.SetSpin(uCh, uCh%2==0?0:1);
    //         pLambda.SetChannelWeight(uCh, uCh%2==0?0.25:0.75);
    //
    //     if(cPotPars1S0&&uCh==0)pLambda.SetShortRangePotential(uCh,0,fDlmPot,*cPotPars1S0);
    //     else if(cPotPars3S1&&uCh==1) pLambda.SetShortRangePotential(uCh,0,fDlmPot,*cPotPars3S1);
    //     printf("1) Running pLambda:\n");
    //     pLambda.KillTheCat();
    //
    //     for (unsigned i = 0; i < count; i++) {
    //       Radius = rmin+step*double(i);
    //
    //       WaveFun[uCh] = pLambda.EvalRadialWaveFunction(uMom,uCh,0,Radius,true);
    //       grWFpL[uCh][uMom].SetPoint(i,Radius,real(WaveFun[uCh]));
    //
    //
    //   }
    // }
    //   else
    if(ExternalWF){
            //for(unsigned uMomBin=0; uMomBin<Kitty.GetNumMomBins(); uMomBin++){
                //Kitty.UseExternalWaveFunction(uMomBin,uCh,0,WaveFunctionU[uMomBin][uCh][0], NumRadBins, RadBins, PhaseShifts[uMomBin][uCh][0]);
//printf("Look at that view (%u)!\n",uCh);
    std::cout<<"In LO 8\n"<<std::endl;
//    void SetExternalWaveFunction(const unsigned& usCh, const unsigned& usPW, DLM_Histo<complex<double>>& histWF, DLM_Histo<complex<double>>& histPS);
                //pLambda.SetExternalWaveFunction(uCh,0,ExternalWF[0][uCh][0],ExternalWF[1][uCh][0]);

    std::cout<<"In LO 9\n"<<std::endl;
                printf("1) Running pLambda:\n");

                //pLambda.KillTheCat();

                for (unsigned i = 0; i < count; i++) {
                  Radius = rmin+step*double(i);
    std::cout<<"In LO 10\n"<<std::endl;
                  WaveFun[uCh] = pLambda.EvalRadialWaveFunction(uMom,uCh,0,Radius,true);
                      std::cout<<"In LO 11\n"<<std::endl;
                  grWFpL[uCh][uMom].SetPoint(i,Radius,real(WaveFun[uCh]));
                  std::cout<<"In LO 12\n"<<std::endl;

              }
              CleanUpWfHisto(pLambda,ExternalWF);

            }
    outFilepLambda->cd();

    grWFpL[uCh]->Write();
    }


  }

    delete outFilepLambda;
    for(unsigned uc=0;uc<NumChannels;uc++){
      delete [] grWFpL[uc];
    }
    delete [] grWFpL;

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
    //Main_pSigma();

    //TestCommonInit();
    //plot_pp();
    //plot_pL();
    //plot_pSigma0();
    //McLevyTest1();
    //testCats();
    //TestSigma0();
    //TestRandBW();
    //SourceSmearing();
    //TestQS();
    //Test_pOmega_Potentials();
    //Test_pp_Dwaves();
    //Test_pSigma0();

    //TestCATS3_NewExtWf("Usmani");
    //TestCATS3_NewExtWf("NLO");
    //TestCATS3_NewExtWf("NLO_Coupled_S");

    //TestVALE_pLambda();

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



