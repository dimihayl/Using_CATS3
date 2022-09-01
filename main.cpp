
#include "ThesisPlots.h"
#include "CommonAnaFunctions.h"
#include "GentleDimi.h"
#include "FemtoBoyzScripts.h"
#include "pSigma.h"
#include "GenBod.h"
#include "MixedEvents.h"
#include "SourceStudies.h"
#include "KaonProton.h"
#include "OtherTasks.h"
#include "DimiMcPlayground.h"
#include "PionAnalysis.h"

#include <iostream>
#include <stdio.h>
#include <string.h>
//#include <omp.h>
#include <complex>

#include "CATS.h"
#include "CATSconstants.h"
#include "CATStools.h"
#include "DLM_Potentials.h"
#include "DLM_Source.h"
#include "DLM_Fitters.h"
#include "DLM_CppTools.h"
#include "DLM_CkDecomposition.h"
#include "CommonAnaFunctions.h"
#include "DLM_Random.h"
#include "DLM_Bessel.h"
#include "DLM_Integration.h"
#include "DLM_WfModel.h"
#include "DLM_Histo.h"
#include "pLambda_1.h"
#include "ForBernie.h"
#include "ForGerhard.h"
#include "UnfoldRoot.h"
#include "DLM_CkModels.h"

#include "TGraph.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TNtuple.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TNtuple.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TLorentzVector.h"
#include "TVector3.h"


#include "ppbar.h"

using namespace std;

/*
void TestCATS3_NewExtWf(TString POT_TYPE){
    //const double kMin=0;
    //const double kMax = 320;
    unsigned NumMomBins;
    double* MomentumBins=NULL;
    double* FitRegion=NULL;

    DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder("/home/dmihaylov/CernBox/CatsFiles");

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

        DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder("/home/dmihaylov/CernBox/CatsFiles");

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

void plot_pp_2(){

    const double SourceSize = 1.0;

    const double kMin=0;
    const double kMax=600;
    const unsigned NumMomBins = 150;

    TFile* OutputFile = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/plot_pp/OutputFile2.root","recreate");

    DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder("/home/dmihaylov/CernBox/CatsFiles");
    CATS Kitty;
    Kitty.SetMaxNumThreads(1);
    Kitty.SetMomBins(NumMomBins,kMin,kMax);
    AnalysisObject.SetUpCats_pp(Kitty,"AV18","Gauss");
    Kitty.SetAnaSource(0,SourceSize);
    //Kitty.SetEpsilonProp(1e-8);
    Kitty.KillTheCat();
    TGraph gKitty;
    gKitty.SetName(TString::Format("gKitty"));
    gKitty.Set(NumMomBins);
    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        //printf("C(%.2f) = %.2f\n",Kitty.GetMomentum(uBin),Kitty.GetCorrFun(uBin));
        gKitty.SetPoint(uBin,Kitty.GetMomentum(uBin),Kitty.GetCorrFun(uBin));
    }
    gKitty.Write();

    delete OutputFile;
}

void plot_pL(){

    const double SourceSize = 1.2;
    const unsigned NumStabilitySteps = 1;
    const double StabMin = 2.0;
    const double StabMax = 2.0;
    const double StabStep = NumStabilitySteps>1?(StabMax-StabMin)/double(NumStabilitySteps-1):0;

    const double kMin=0;
    const double kFine=200;
    const unsigned NumFineBins = 40;
    const double kFineStep = (kFine-kMin)/double(NumFineBins);
    const double kMax=340;
    const unsigned NumCoarseBins = 28;
    const double kCoarseStep = (kMax-kFine)/double(NumCoarseBins);
    const unsigned NumMomBins=NumFineBins+NumCoarseBins;
    double* MomBins = new double[NumMomBins];
    for(unsigned uBin=0; uBin<=NumFineBins; uBin++){
        MomBins[uBin] = kMin+double(uBin)*kFineStep;
    }
    for(unsigned uBin=NumFineBins+1; uBin<=NumMomBins; uBin++){
        MomBins[uBin] = MomBins[uBin-1]+kCoarseStep;;
    }

    TFile* OutputFile = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/plot_pL/NLO.root","recreate");

    for(unsigned uStab=0; uStab<NumStabilitySteps; uStab++){

        DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder("/home/dmihaylov/CernBox/CatsFiles");

        double Stability = StabMin+double(uStab)*StabStep;
        printf("\nStability = %.2f\n",Stability);
        CATS KittyNLO;
        KittyNLO.SetMaxNumThreads(1);
        KittyNLO.SetMomBins(NumMomBins,MomBins);
        AnalysisObject.SetUpCats_pL(KittyNLO,"NLO","Gauss");//NLO_Coupled_S
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

    DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder("/home/dmihaylov/CernBox/CatsFiles");

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


void PlotResonanceAngles(const double& CoreRadius, const double& ResoFreePath){
    const double ThetaMin=0;
    const double ThetaMax=360;
    const unsigned NumThetaSteps = 361;
    const double ThetaStepSize = (ThetaMax-ThetaMin)/double(NumThetaSteps);

    double Theta;
    double Alpha;
    double DeltaRadius;
    double Radius;
    double PathModification;

    TGraph g_Theta_PathModif;
    g_Theta_PathModif.SetName("g_Theta_PathModif");
    g_Theta_PathModif.Set(NumThetaSteps);
    TGraph g_Theta_Radius;
    g_Theta_Radius.SetName("g_Theta_Radius");
    g_Theta_Radius.Set(NumThetaSteps);
    TGraph g_Theta_DeltaRadius;
    g_Theta_DeltaRadius.SetName("g_Theta_DeltaRadius");
    g_Theta_DeltaRadius.Set(NumThetaSteps);
    TGraph g_Theta_Alpha;
    g_Theta_Alpha.SetName("g_Theta_Alpha");
    g_Theta_Alpha.Set(NumThetaSteps);

    for(unsigned uTheta=0; uTheta<NumThetaSteps; uTheta++){
        Theta = ThetaMin + double(uTheta)*ThetaStepSize;
        //Alpha = 0.5*atan(ResoFreePath*sin(Theta*DegToRad)/(CoreRadius+ResoFreePath*cos(Theta*DegToRad)))*RadToDeg+90.;
        //DeltaRadius = ResoFreePath*(sin(Theta*DegToRad)/tan(Alpha*DegToRad)-cos(Theta*DegToRad));
        ////DeltaRadius = (sin(Theta*DegToRad-Alpha*DegToRad)/sin(Alpha*DegToRad));
        //Radius = CoreRadius + DeltaRadius;
        //PathModification = DeltaRadius/ResoFreePath;

        //g_Theta_PathModif.SetPoint(uTheta,Theta,PathModification==PathModification?PathModification:0);
        //g_Theta_Radius.SetPoint(uTheta,Theta,Radius==Radius?Radius:0);
        //g_Theta_DeltaRadius.SetPoint(uTheta,Theta,DeltaRadius==DeltaRadius?DeltaRadius:0);
        //g_Theta_Alpha.SetPoint(uTheta,Theta,Alpha==Alpha?Alpha:0);


        Radius = sqrt(CoreRadius*CoreRadius+ResoFreePath*ResoFreePath-2.*CoreRadius*ResoFreePath*cos(Theta*DegToRad));
        PathModification = Radius/(CoreRadius+ResoFreePath);
        g_Theta_Radius.SetPoint(uTheta,Theta,Radius==Radius?Radius:0);
        g_Theta_PathModif.SetPoint(uTheta,Theta,PathModification==PathModification?PathModification:0);
    }

    TFile* fOut = new TFile(TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/TestThetaSource/fAngle_rc%.2f_s%.2f.root",CoreRadius,ResoFreePath),"recreate");

    g_Theta_PathModif.Write();
    g_Theta_Radius.Write();
    //g_Theta_DeltaRadius.Write();
    //g_Theta_Alpha.Write();

    delete fOut;
}

void TestThetaSource(){

    const unsigned NumMomBins = 20;
    const double kMin = 0;
    const double kMax = 200;

    DLM_CommonAnaFunctions AnaObject;

    CATS AB_pp_G;
    AB_pp_G.SetMomBins(NumMomBins,kMin,kMax);//Gauss
    AnaObject.SetUpCats_pp(AB_pp_G,"AV18","McGauss_Reso");
    //AB_pp_G.SetAnaSource(0,1.2);

    CATS AB_pp_GT;
    AB_pp_GT.SetMomBins(NumMomBins,kMin,kMax);//GaussTheta
    AnaObject.SetUpCats_pp(AB_pp_GT,"AV18","GaussTheta");
    //AB_pp_GT.SetAnaSource(0,1.2);

    AB_pp_G.SetMaxNumThreads(1);
    AB_pp_G.KillTheCat();
    AB_pp_GT.SetMaxNumThreads(1);
    AB_pp_GT.KillTheCat();

    TFile* fOut = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/TestThetaSource/fOut.root","recreate");

    TGraph gKitty_G;
    TGraph gKitty_GT;
    gKitty_G.SetName(TString::Format("gKitty_G"));
    gKitty_G.Set(NumMomBins);
    gKitty_GT.SetName(TString::Format("gKitty_GT"));
    gKitty_GT.Set(NumMomBins);
    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        printf("C(%.2f) = %.3f vs %.3f\n",AB_pp_G.GetMomentum(uBin),AB_pp_G.GetCorrFun(uBin),AB_pp_GT.GetCorrFun(uBin));
        gKitty_G.SetPoint(uBin,AB_pp_G.GetMomentum(uBin),AB_pp_G.GetCorrFun(uBin));
        gKitty_GT.SetPoint(uBin,AB_pp_GT.GetMomentum(uBin),AB_pp_GT.GetCorrFun(uBin));

    }
    gKitty_G.Write();
    gKitty_GT.Write();


    delete fOut;
}

void MickeyMouseEffectOnSourceBasedOnTheta(){
    const double Core_pXi_Back = 0.80;
    const double Core_pXi_Epos = 0.88;
    const double Core_pXi_Rand = 1.06;
    const double Core_pXi_TM_Back = 0.80;
    const double Core_pXi_TM_Epos = 0.98;

    const double Core_pOmega_Back = 0.73;
    const double Core_pOmega_Epos = 0.82;
    const double Core_pOmega_Rand = 0.96;
    const double Core_pOmega_TM_Back = 0.73;
    const double Core_pOmega_TM_Epos = 0.921;//0.89

    const double Core_pp_Back = 0.73;
    const double Core_pp_Epos = 0.82;
    const double Core_pp_Rand = 0.96;
    const double Core_pp_TM_Back = 0.73;
    //const double Core_pp_TM_Epos = 0.921;
    //const double Core_pp_TM_Epos = 1.11483;
    // const double Core_pp_TM_Epos = 1.07645;
    //const double Core_pp_TM_Epos = 1.15142;
    //const double Core_pp_TM_Epos = 1.2;
    //for pAp
    //const double Core_pp_TM_Epos = 1.099;
    //const double Core_pp_TM_Epos = 1.063;
    const double Core_pp_TM_Epos = 1.025;

    const double Core_pL_TM_Epos = 1.2;

    //const double Core_pXi_Back = 1.0;
    //const double Core_pXi_Epos = 1.0;
    //const double Core_pXi_Rand = 1.0;

    //const double Core_pOmega_Back = 0.80;
    //const double Core_pOmega_Epos = 0.80;
    //const double Core_pOmega_Rand = 0.80;


    const double LambdaPar_pXi = 0.53;
    const double LambdaPar_pOmega = 0.62;
    const double LambdaPar_pp = 0.72;
    const double LambdaPar_pL = 0.5;

    const double kMin = 0;
    const double kMax = 300;

    const unsigned NumMomBins = 30;

    DLM_CommonAnaFunctions AO_pXi_Back;
    CATS Kitty_pXi_Back;
    Kitty_pXi_Back.SetMaxNumThreads(1);
    Kitty_pXi_Back.SetMomBins(NumMomBins,kMin,kMax);//pXim_HALQCD1
    AO_pXi_Back.SetUpCats_pXim(Kitty_pXi_Back,"pXim_HALQCD1","McGauss_Reso",0,0);//McGauss_Reso
    //AO_pXi_Back.SetUpCats_pp(Kitty_pXi_Back,"AV18","McGauss_Reso",0,0);
    Kitty_pXi_Back.SetAnaSource(0,Core_pXi_Back);
    Kitty_pXi_Back.KillTheCat();

    DLM_CommonAnaFunctions AO_pXi_Rand;
    CATS Kitty_pXi_Rand;
    Kitty_pXi_Rand.SetMaxNumThreads(1);
    Kitty_pXi_Rand.SetMomBins(NumMomBins,kMin,kMax);
    AO_pXi_Rand.SetUpCats_pXim(Kitty_pXi_Rand,"pXim_HALQCD1","McGauss_Reso",0,1);
    //AO_pXi_Rand.SetUpCats_pp(Kitty_pXi_Rand,"AV18","McGauss_Reso",0,2);
    Kitty_pXi_Rand.SetAnaSource(0,Core_pXi_Rand);
    Kitty_pXi_Rand.KillTheCat();

    DLM_CommonAnaFunctions AO_pXi_Epos;
    CATS Kitty_pXi_Epos;
    Kitty_pXi_Epos.SetMaxNumThreads(1);
    Kitty_pXi_Epos.SetMomBins(NumMomBins,kMin,kMax);
    AO_pXi_Epos.SetUpCats_pXim(Kitty_pXi_Epos,"pXim_HALQCD1","McGauss_Reso",0,2);//McGauss_Reso
    //AO_pXi_Epos.SetUpCats_pp(Kitty_pXi_Epos,"AV18","McGauss_Reso",0,0);
    Kitty_pXi_Epos.SetAnaSource(0,Core_pXi_Epos);
    Kitty_pXi_Epos.KillTheCat();

    DLM_CommonAnaFunctions AO_pXi_TM_Back;
    CATS Kitty_pXi_TM_Back;
    Kitty_pXi_TM_Back.SetMaxNumThreads(1);
    Kitty_pXi_TM_Back.SetMomBins(NumMomBins,kMin,kMax);//pXim_HALQCD1
    AO_pXi_TM_Back.SetUpCats_pXim(Kitty_pXi_TM_Back,"pXim_HALQCD1","McGauss_ResoTM",0,1);//McGauss_Reso
    //AO_pXi_TM_Back.SetUpCats_pp(Kitty_pXi_TM_Back,"AV18","McGauss_Reso",0,0);
    Kitty_pXi_TM_Back.SetAnaSource(0,Core_pXi_TM_Back);
    Kitty_pXi_TM_Back.KillTheCat();

    DLM_CommonAnaFunctions AO_pXi_TM_Epos;
    CATS Kitty_pXi_TM_Epos;
    Kitty_pXi_TM_Epos.SetMaxNumThreads(1);
    Kitty_pXi_TM_Epos.SetMomBins(NumMomBins,kMin,kMax);
    AO_pXi_TM_Epos.SetUpCats_pXim(Kitty_pXi_TM_Epos,"pXim_HALQCD1","McGauss_ResoTM",0,202);//McGauss_Reso
    //AO_pXi_TM_Epos.SetUpCats_pp(Kitty_pXi_TM_Epos,"AV18","McGauss_Reso",0,0);
    Kitty_pXi_TM_Epos.SetAnaSource(0,Core_pXi_TM_Epos);
    Kitty_pXi_TM_Epos.KillTheCat();

    DLM_CommonAnaFunctions AO_pOmega_Back;
    CATS Kitty_pOmega_Back;
    Kitty_pOmega_Back.SetMaxNumThreads(1);
    Kitty_pOmega_Back.SetMomBins(NumMomBins,kMin,kMax);
    AO_pOmega_Back.SetUpCats_pOmegam(Kitty_pOmega_Back,"pOmega_Lattice","McGauss_Reso",0,0);
    Kitty_pOmega_Back.SetAnaSource(0,Core_pOmega_Back);
    Kitty_pOmega_Back.KillTheCat();

    DLM_CommonAnaFunctions AO_pOmega_Rand;
    CATS Kitty_pOmega_Rand;
    Kitty_pOmega_Rand.SetMaxNumThreads(1);
    Kitty_pOmega_Rand.SetMomBins(NumMomBins,kMin,kMax);
    AO_pOmega_Rand.SetUpCats_pOmegam(Kitty_pOmega_Rand,"pOmega_Lattice","McGauss_Reso",0,1);
    Kitty_pOmega_Rand.SetAnaSource(0,Core_pOmega_Rand);
    Kitty_pOmega_Rand.KillTheCat();

    DLM_CommonAnaFunctions AO_pOmega_Epos;
    CATS Kitty_pOmega_Epos;
    Kitty_pOmega_Epos.SetMaxNumThreads(1);
    Kitty_pOmega_Epos.SetMomBins(NumMomBins,kMin,kMax);
    AO_pOmega_Epos.SetUpCats_pOmegam(Kitty_pOmega_Epos,"pOmega_Lattice","McGauss_Reso",0,2);
    Kitty_pOmega_Epos.SetAnaSource(0,Core_pOmega_Epos);
    Kitty_pOmega_Epos.KillTheCat();

    DLM_CommonAnaFunctions AO_pOmega_TM_Back;
    CATS Kitty_pOmega_TM_Back;
    Kitty_pOmega_TM_Back.SetMaxNumThreads(1);
    Kitty_pOmega_TM_Back.SetMomBins(NumMomBins,kMin,kMax);
    AO_pOmega_TM_Back.SetUpCats_pOmegam(Kitty_pOmega_TM_Back,"pOmega_Lattice","McGauss_ResoTM",0,1);
    Kitty_pOmega_TM_Back.SetAnaSource(0,Core_pOmega_TM_Back);
    Kitty_pOmega_TM_Back.KillTheCat();

    DLM_CommonAnaFunctions AO_pOmega_TM_Epos;
    CATS Kitty_pOmega_TM_Epos;
    Kitty_pOmega_TM_Epos.SetMaxNumThreads(1);
    Kitty_pOmega_TM_Epos.SetMomBins(NumMomBins,kMin,kMax);
    AO_pOmega_TM_Epos.SetUpCats_pOmegam(Kitty_pOmega_TM_Epos,"pOmega_Lattice","McGauss_ResoTM",0,202);
    Kitty_pOmega_TM_Epos.SetAnaSource(0,Core_pOmega_TM_Epos);
    Kitty_pOmega_TM_Epos.SetAnaSource(1,2);
    Kitty_pOmega_TM_Epos.KillTheCat();

    DLM_CommonAnaFunctions AO_pp_Back;
    CATS Kitty_pp_Back;
    Kitty_pp_Back.SetMaxNumThreads(1);
    Kitty_pp_Back.SetMomBins(NumMomBins,kMin,kMax);
    AO_pp_Back.SetUpCats_pp(Kitty_pp_Back,"AV18","McGauss_Reso",0,0);
    Kitty_pp_Back.SetAnaSource(0,Core_pp_Back);
    Kitty_pp_Back.KillTheCat();

    DLM_CommonAnaFunctions AO_pp_Rand;
    CATS Kitty_pp_Rand;
    Kitty_pp_Rand.SetMaxNumThreads(1);
    Kitty_pp_Rand.SetMomBins(NumMomBins,kMin,kMax);
    AO_pp_Rand.SetUpCats_pp(Kitty_pp_Rand,"AV18","McGauss_Reso",0,1);
    Kitty_pp_Rand.SetAnaSource(0,Core_pp_Rand);
    Kitty_pp_Rand.KillTheCat();

    DLM_CommonAnaFunctions AO_pp_Epos;
    CATS Kitty_pp_Epos;
    Kitty_pp_Epos.SetMaxNumThreads(1);
    Kitty_pp_Epos.SetMomBins(NumMomBins,kMin,kMax);
    AO_pp_Epos.SetUpCats_pp(Kitty_pp_Epos,"AV18","McGauss_Reso",0,2);
    Kitty_pp_Epos.SetAnaSource(0,Core_pp_Epos);
    Kitty_pp_Epos.KillTheCat();

    DLM_CommonAnaFunctions AO_pp_TM_Back;
    CATS Kitty_pp_TM_Back;
    Kitty_pp_TM_Back.SetMaxNumThreads(1);
    Kitty_pp_TM_Back.SetMomBins(NumMomBins,kMin,kMax);
    AO_pp_TM_Back.SetUpCats_pp(Kitty_pp_TM_Back,"AV18","McGauss_ResoTM",0,1);
    Kitty_pp_TM_Back.SetAnaSource(0,Core_pp_TM_Back);
    Kitty_pp_TM_Back.KillTheCat();

    DLM_CommonAnaFunctions AO_pp_TM_Epos;
    CATS Kitty_pp_TM_Epos;
    Kitty_pp_TM_Epos.SetMaxNumThreads(1);
    Kitty_pp_TM_Epos.SetMomBins(NumMomBins,kMin,kMax);
    AO_pp_TM_Epos.SetUpCats_pp(Kitty_pp_TM_Epos,"AV18","McGauss_ResoTM",0,202);
    Kitty_pp_TM_Epos.SetAnaSource(0,Core_pp_TM_Epos);
    Kitty_pp_TM_Epos.SetAnaSource(1,2);
    Kitty_pp_TM_Epos.KillTheCat();

    DLM_CommonAnaFunctions AO_pL_TM_Epos;
    CATS Kitty_pL_TM_Epos;
    Kitty_pL_TM_Epos.SetMaxNumThreads(1);
    Kitty_pL_TM_Epos.SetMomBins(NumMomBins,kMin,kMax);
    AO_pL_TM_Epos.SetUpCats_pL(Kitty_pL_TM_Epos,"Usmani","McGauss_ResoTM",0,202);
    Kitty_pL_TM_Epos.SetAnaSource(0,Core_pL_TM_Epos);
    Kitty_pL_TM_Epos.SetAnaSource(1,2);
    Kitty_pL_TM_Epos.KillTheCat();

    TFile* fOut = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MickeyMouseEffectOnSourceBasedOnTheta/fOut_TM.root","recreate");

    TGraph gKitty_pXi_Back;
    TGraph gKitty_pXi_Epos;
    TGraph gKitty_pXi_Rand;
    TGraph gKitty_pXi_TM_Back;
    TGraph gKitty_pXi_TM_Epos;
    gKitty_pXi_Back.SetName(TString::Format("gKitty_pXi_Back"));
    gKitty_pXi_Epos.SetName(TString::Format("gKitty_pXi_Epos"));
    gKitty_pXi_Rand.SetName(TString::Format("gKitty_pXi_Rand"));
    gKitty_pXi_TM_Back.SetName(TString::Format("gKitty_pXi_TM_Back"));
    gKitty_pXi_TM_Epos.SetName(TString::Format("gKitty_pXi_TM_Epos"));
    gKitty_pXi_Back.Set(NumMomBins);
    gKitty_pXi_Epos.Set(NumMomBins);
    gKitty_pXi_Rand.Set(NumMomBins);
    gKitty_pXi_TM_Back.Set(NumMomBins);
    gKitty_pXi_TM_Epos.Set(NumMomBins);

    TGraph gKitty_pOmega_Back;
    TGraph gKitty_pOmega_Epos;
    TGraph gKitty_pOmega_Rand;
    TGraph gKitty_pOmega_TM_Back;
    TGraph gKitty_pOmega_TM_Epos;
    gKitty_pOmega_Back.SetName(TString::Format("gKitty_pOmega_Back"));
    gKitty_pOmega_Epos.SetName(TString::Format("gKitty_pOmega_Epos"));
    gKitty_pOmega_Rand.SetName(TString::Format("gKitty_pOmega_Rand"));
    gKitty_pOmega_TM_Back.SetName(TString::Format("gKitty_pOmega_TM_Back"));
    gKitty_pOmega_TM_Epos.SetName(TString::Format("gKitty_pOmega_TM_Epos"));
    gKitty_pOmega_Back.Set(NumMomBins);
    gKitty_pOmega_Epos.Set(NumMomBins);
    gKitty_pOmega_Rand.Set(NumMomBins);
    gKitty_pOmega_TM_Back.Set(NumMomBins);
    gKitty_pOmega_TM_Epos.Set(NumMomBins);

    TGraph gKitty_pp_Back;
    TGraph gKitty_pp_Epos;
    TGraph gKitty_pp_Rand;
    TGraph gKitty_pp_TM_Back;
    TGraph gKitty_pp_TM_Epos;
    TGraph gKitty_pL_TM_Epos;
    gKitty_pp_Back.SetName(TString::Format("gKitty_pp_Back"));
    gKitty_pp_Epos.SetName(TString::Format("gKitty_pp_Epos"));
    gKitty_pp_Rand.SetName(TString::Format("gKitty_pp_Rand"));
    gKitty_pp_TM_Back.SetName(TString::Format("gKitty_pp_TM_Back"));
    gKitty_pp_TM_Epos.SetName(TString::Format("gKitty_pp_TM_Epos"));
    gKitty_pL_TM_Epos.SetName(TString::Format("gKitty_pL_TM_Epos"));
    gKitty_pp_Back.Set(NumMomBins);
    gKitty_pp_Epos.Set(NumMomBins);
    gKitty_pp_Rand.Set(NumMomBins);
    gKitty_pp_TM_Back.Set(NumMomBins);
    gKitty_pp_TM_Epos.Set(NumMomBins);
    gKitty_pL_TM_Epos.Set(NumMomBins);

    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        //printf("C(%.2f) = %.3f vs %.3f\n",AB_pp_G.GetMomentum(uBin),AB_pp_G.GetCorrFun(uBin),AB_pp_GT.GetCorrFun(uBin));
        gKitty_pXi_Back.SetPoint(uBin,Kitty_pXi_Back.GetMomentum(uBin),LambdaPar_pXi*Kitty_pXi_Back.GetCorrFun(uBin)+1.-LambdaPar_pXi);
        gKitty_pXi_Epos.SetPoint(uBin,Kitty_pXi_Epos.GetMomentum(uBin),LambdaPar_pXi*Kitty_pXi_Epos.GetCorrFun(uBin)+1.-LambdaPar_pXi);
        gKitty_pXi_Rand.SetPoint(uBin,Kitty_pXi_Rand.GetMomentum(uBin),LambdaPar_pXi*Kitty_pXi_Rand.GetCorrFun(uBin)+1.-LambdaPar_pXi);
        gKitty_pXi_TM_Back.SetPoint(uBin,Kitty_pXi_TM_Back.GetMomentum(uBin),LambdaPar_pXi*Kitty_pXi_TM_Back.GetCorrFun(uBin)+1.-LambdaPar_pXi);
        gKitty_pXi_TM_Epos.SetPoint(uBin,Kitty_pXi_TM_Epos.GetMomentum(uBin),LambdaPar_pXi*Kitty_pXi_TM_Epos.GetCorrFun(uBin)+1.-LambdaPar_pXi);
        gKitty_pOmega_Back.SetPoint(uBin,Kitty_pOmega_Back.GetMomentum(uBin),LambdaPar_pOmega*Kitty_pOmega_Back.GetCorrFun(uBin)+1.-LambdaPar_pOmega);
        gKitty_pOmega_Epos.SetPoint(uBin,Kitty_pOmega_Epos.GetMomentum(uBin),LambdaPar_pOmega*Kitty_pOmega_Epos.GetCorrFun(uBin)+1.-LambdaPar_pOmega);
        gKitty_pOmega_Rand.SetPoint(uBin,Kitty_pOmega_Rand.GetMomentum(uBin),LambdaPar_pOmega*Kitty_pOmega_Rand.GetCorrFun(uBin)+1.-LambdaPar_pOmega);
        gKitty_pOmega_TM_Back.SetPoint(uBin,Kitty_pOmega_TM_Back.GetMomentum(uBin),LambdaPar_pOmega*Kitty_pOmega_TM_Back.GetCorrFun(uBin)+1.-LambdaPar_pOmega);
        gKitty_pOmega_TM_Epos.SetPoint(uBin,Kitty_pOmega_TM_Epos.GetMomentum(uBin),LambdaPar_pOmega*Kitty_pOmega_TM_Epos.GetCorrFun(uBin)+1.-LambdaPar_pOmega);
        gKitty_pp_Back.SetPoint(uBin,Kitty_pp_Back.GetMomentum(uBin),LambdaPar_pp*Kitty_pp_Back.GetCorrFun(uBin)+1.-LambdaPar_pp);
        gKitty_pp_Epos.SetPoint(uBin,Kitty_pp_Epos.GetMomentum(uBin),LambdaPar_pp*Kitty_pp_Epos.GetCorrFun(uBin)+1.-LambdaPar_pp);
        gKitty_pp_Rand.SetPoint(uBin,Kitty_pp_Rand.GetMomentum(uBin),LambdaPar_pp*Kitty_pp_Rand.GetCorrFun(uBin)+1.-LambdaPar_pp);
        gKitty_pp_TM_Back.SetPoint(uBin,Kitty_pp_TM_Back.GetMomentum(uBin),LambdaPar_pp*Kitty_pp_TM_Back.GetCorrFun(uBin)+1.-LambdaPar_pp);
        gKitty_pp_TM_Epos.SetPoint(uBin,Kitty_pp_TM_Epos.GetMomentum(uBin),LambdaPar_pp*Kitty_pp_TM_Epos.GetCorrFun(uBin)+1.-LambdaPar_pp);
        gKitty_pL_TM_Epos.SetPoint(uBin,Kitty_pL_TM_Epos.GetMomentum(uBin),LambdaPar_pL*Kitty_pL_TM_Epos.GetCorrFun(uBin)+1.-LambdaPar_pL);
    }

    TGraph gSource_pXi_Back;
    TGraph gSource_pXi_Epos;
    TGraph gSource_pXi_Rand;
    TGraph gSource_pXi_TM_Back;
    TGraph gSource_pXi_TM_Epos;
    gSource_pXi_Back.SetName(TString::Format("gSource_pXi_Back"));
    gSource_pXi_Epos.SetName(TString::Format("gSource_pXi_Epos"));
    gSource_pXi_Rand.SetName(TString::Format("gSource_pXi_Rand"));
    gSource_pXi_TM_Back.SetName(TString::Format("gSource_pXi_TM_Back"));
    gSource_pXi_TM_Epos.SetName(TString::Format("gSource_pXi_TM_Epos"));

    TGraph gSource_pOmega_Back;
    TGraph gSource_pOmega_Epos;
    TGraph gSource_pOmega_Rand;
    TGraph gSource_pOmega_TM_Back;
    TGraph gSource_pOmega_TM_Epos;
    gSource_pOmega_Back.SetName(TString::Format("gSource_pOmega_Back"));
    gSource_pOmega_Epos.SetName(TString::Format("gSource_pOmega_Epos"));
    gSource_pOmega_Rand.SetName(TString::Format("gSource_pOmega_Rand"));
    gSource_pOmega_TM_Back.SetName(TString::Format("gSource_pOmega_TM_Back"));
    gSource_pOmega_TM_Epos.SetName(TString::Format("gSource_pOmega_TM_Epos"));

    const unsigned NumRadBins = 256;
    const double rMin = 0;
    const double rMax = 6;
    TH1F* hSource_pXi_Back = new TH1F("hSource_pXi_Back","hSource_pXi_Back",NumRadBins,rMin,rMax);
    TH1F* hSource_pXi_Epos = new TH1F("hSource_pXi_Epos","hSource_pXi_Epos",NumRadBins,rMin,rMax);
    TH1F* hSource_pXi_Rand = new TH1F("hSource_pXi_Rand","hSource_pXi_Rand",NumRadBins,rMin,rMax);
    TH1F* hSource_pXi_TM_Back = new TH1F("hSource_pXi_TM_Back","hSource_pXi_TM_Back",NumRadBins,rMin,rMax);
    TH1F* hSource_pXi_TM_Epos = new TH1F("hSource_pXi_TM_Epos","hSource_pXi_TM_Epos",NumRadBins,rMin,rMax);

    TH1F* hSource_pOmega_Back = new TH1F("hSource_pOmega_Back","hSource_pOmega_Back",NumRadBins,rMin,rMax);
    TH1F* hSource_pOmega_Epos = new TH1F("hSource_pOmega_Epos","hSource_pOmega_Epos",NumRadBins,rMin,rMax);
    TH1F* hSource_pOmega_Rand = new TH1F("hSource_pOmega_Rand","hSource_pOmega_Rand",NumRadBins,rMin,rMax);
    TH1F* hSource_pOmega_TM_Back = new TH1F("hSource_pOmega_TM_Back","hSource_pOmega_TM_Back",NumRadBins,rMin,rMax);
    TH1F* hSource_pOmega_TM_Epos = new TH1F("hSource_pOmega_TM_Epos","hSource_pOmega_TM_Epos",NumRadBins,rMin,rMax);

    TH1F* hSource_pp_Back = new TH1F("hSource_pp_Back","hSource_pp_Back",NumRadBins,rMin,rMax);
    TH1F* hSource_pp_Epos = new TH1F("hSource_pp_Epos","hSource_pp_Epos",NumRadBins,rMin,rMax);
    TH1F* hSource_pp_Rand = new TH1F("hSource_pp_Rand","hSource_pp_Rand",NumRadBins,rMin,rMax);
    TH1F* hSource_pp_TM_Back = new TH1F("hSource_pp_TM_Back","hSource_pp_TM_Back",NumRadBins,rMin,rMax);
    TH1F* hSource_pp_TM_Epos = new TH1F("hSource_pp_TM_Epos","hSource_pp_TM_Epos",NumRadBins,rMin,rMax);

    TH1F* hSource_pL_TM_Epos = new TH1F("hSource_pL_TM_Epos","hSource_pL_TM_Epos",NumRadBins,rMin,rMax);

    //const unsigned NumRadBins = 128;
    //const double rMin = 0.1;
    //const double rMax = 10;
    //const double rStep = (rMax-rMin)/double(NumRadBins-1);
    //EvaluateTheSource
    for(unsigned uRad=0; uRad<NumRadBins; uRad++){
        //double RAD = rMin+double(uRad)*rStep;
        double RAD = hSource_pXi_Back->GetBinCenter(uRad+1);
        gSource_pXi_Back.SetPoint(uRad,RAD,Kitty_pXi_Back.EvaluateTheSource(0,RAD,0));
        gSource_pXi_Epos.SetPoint(uRad,RAD,Kitty_pXi_Epos.EvaluateTheSource(0,RAD,0));
        gSource_pXi_Rand.SetPoint(uRad,RAD,Kitty_pXi_Rand.EvaluateTheSource(0,RAD,0));
        gSource_pXi_TM_Back.SetPoint(uRad,RAD,Kitty_pXi_TM_Back.EvaluateTheSource(0,RAD,0));
        gSource_pXi_TM_Epos.SetPoint(uRad,RAD,Kitty_pXi_TM_Epos.EvaluateTheSource(0,RAD,0));
        gSource_pOmega_Back.SetPoint(uRad,RAD,Kitty_pOmega_Back.EvaluateTheSource(0,RAD,0));
        gSource_pOmega_Epos.SetPoint(uRad,RAD,Kitty_pOmega_Epos.EvaluateTheSource(0,RAD,0));
        gSource_pOmega_Rand.SetPoint(uRad,RAD,Kitty_pOmega_Rand.EvaluateTheSource(0,RAD,0));
        gSource_pOmega_TM_Back.SetPoint(uRad,RAD,Kitty_pOmega_TM_Back.EvaluateTheSource(0,RAD,0));
        gSource_pOmega_TM_Epos.SetPoint(uRad,RAD,Kitty_pOmega_TM_Epos.EvaluateTheSource(0,RAD,0));

        hSource_pXi_Back->SetBinContent(uRad+1,Kitty_pXi_Back.EvaluateTheSource(0,RAD,0));
        hSource_pXi_Epos->SetBinContent(uRad+1,Kitty_pXi_Epos.EvaluateTheSource(0,RAD,0));
        hSource_pXi_Rand->SetBinContent(uRad+1,Kitty_pXi_Rand.EvaluateTheSource(0,RAD,0));
        hSource_pXi_TM_Back->SetBinContent(uRad+1,Kitty_pXi_TM_Back.EvaluateTheSource(0,RAD,0));
        hSource_pXi_TM_Epos->SetBinContent(uRad+1,Kitty_pXi_TM_Epos.EvaluateTheSource(0,RAD,0));

        hSource_pOmega_Back->SetBinContent(uRad+1,Kitty_pOmega_Back.EvaluateTheSource(0,RAD,0));
        hSource_pOmega_Epos->SetBinContent(uRad+1,Kitty_pOmega_Epos.EvaluateTheSource(0,RAD,0));
        hSource_pOmega_Rand->SetBinContent(uRad+1,Kitty_pOmega_Rand.EvaluateTheSource(0,RAD,0));
        hSource_pOmega_TM_Back->SetBinContent(uRad+1,Kitty_pOmega_TM_Back.EvaluateTheSource(0,RAD,0));
        hSource_pOmega_TM_Epos->SetBinContent(uRad+1,Kitty_pOmega_TM_Epos.EvaluateTheSource(0,RAD,0));

        hSource_pp_Back->SetBinContent(uRad+1,Kitty_pp_Back.EvaluateTheSource(0,RAD,0));
        hSource_pp_Epos->SetBinContent(uRad+1,Kitty_pp_Epos.EvaluateTheSource(0,RAD,0));
        hSource_pp_Rand->SetBinContent(uRad+1,Kitty_pp_Rand.EvaluateTheSource(0,RAD,0));
        hSource_pp_TM_Back->SetBinContent(uRad+1,Kitty_pp_TM_Back.EvaluateTheSource(0,RAD,0));
        hSource_pp_TM_Epos->SetBinContent(uRad+1,Kitty_pp_TM_Epos.EvaluateTheSource(0,RAD,0));

        hSource_pL_TM_Epos->SetBinContent(uRad+1,Kitty_pL_TM_Epos.EvaluateTheSource(0,RAD,0));

        hSource_pXi_Back->SetBinError(uRad+1,hSource_pXi_Back->GetBinContent(uRad+1)*0.001+0.001);
        hSource_pXi_Epos->SetBinError(uRad+1,hSource_pXi_Epos->GetBinContent(uRad+1)*0.001+0.001);
        hSource_pXi_Rand->SetBinError(uRad+1,hSource_pXi_Rand->GetBinContent(uRad+1)*0.001+0.001);
        hSource_pXi_TM_Back->SetBinError(uRad+1,hSource_pXi_TM_Back->GetBinContent(uRad+1)*0.001+0.001);
        hSource_pXi_TM_Epos->SetBinError(uRad+1,hSource_pXi_TM_Epos->GetBinContent(uRad+1)*0.001+0.001);

        hSource_pOmega_Back->SetBinError(uRad+1,hSource_pOmega_Back->GetBinContent(uRad+1)*0.001+0.001);
        hSource_pOmega_Epos->SetBinError(uRad+1,hSource_pOmega_Epos->GetBinContent(uRad+1)*0.001+0.001);
        hSource_pOmega_Rand->SetBinError(uRad+1,hSource_pOmega_Rand->GetBinContent(uRad+1)*0.001+0.001);
        hSource_pOmega_TM_Back->SetBinError(uRad+1,hSource_pOmega_TM_Back->GetBinContent(uRad+1)*0.001+0.001);
        hSource_pOmega_TM_Epos->SetBinError(uRad+1,hSource_pOmega_TM_Epos->GetBinContent(uRad+1)*0.001+0.001);

        hSource_pp_Back->SetBinError(uRad+1,hSource_pp_Back->GetBinContent(uRad+1)*0.001+0.001);
        hSource_pp_Epos->SetBinError(uRad+1,hSource_pp_Epos->GetBinContent(uRad+1)*0.001+0.001);
        hSource_pp_Rand->SetBinError(uRad+1,hSource_pp_Rand->GetBinContent(uRad+1)*0.001+0.001);
        hSource_pp_TM_Back->SetBinError(uRad+1,hSource_pp_TM_Back->GetBinContent(uRad+1)*0.001+0.001);
        hSource_pp_TM_Epos->SetBinError(uRad+1,hSource_pp_TM_Epos->GetBinContent(uRad+1)*0.001+0.001);

        hSource_pL_TM_Epos->SetBinError(uRad+1,hSource_pL_TM_Epos->GetBinContent(uRad+1)*0.001+0.001);
    }

    TF1* fSource_pXi_Back = new TF1(TString::Format("fSource_pXi_Back"),GaussSourceTF1,rMin,rMax,1);
    fSource_pXi_Back->SetParameter(0,1.0);
    fSource_pXi_Back->SetParLimits(0,0.5,2.0);
    hSource_pXi_Back->Fit(fSource_pXi_Back,"Q, S, N, R, M");

    TF1* fSource_pXi_Epos = new TF1(TString::Format("fSource_pXi_Epos"),GaussSourceTF1,rMin,rMax,1);
    fSource_pXi_Epos->SetParameter(0,1.0);
    fSource_pXi_Epos->SetParLimits(0,0.5,2.0);
    hSource_pXi_Epos->Fit(fSource_pXi_Epos,"Q, S, N, R, M");

    TF1* fSource_pXi_Rand = new TF1(TString::Format("fSource_pXi_Rand"),GaussSourceTF1,rMin,rMax,1);
    fSource_pXi_Rand->SetParameter(0,1.0);
    fSource_pXi_Rand->SetParLimits(0,0.5,2.0);
    hSource_pXi_Rand->Fit(fSource_pXi_Rand,"Q, S, N, R, M");

    TF1* fSource_pXi_TM_Back = new TF1(TString::Format("fSource_pXi_TM_Back"),GaussSourceTF1,rMin,rMax,1);
    fSource_pXi_TM_Back->SetParameter(0,1.0);
    fSource_pXi_TM_Back->SetParLimits(0,0.5,2.0);
    hSource_pXi_TM_Back->Fit(fSource_pXi_TM_Back,"Q, S, N, R, M");

    TF1* fSource_pXi_TM_Epos = new TF1(TString::Format("fSource_pXi_TM_Epos"),GaussSourceTF1,rMin,rMax,1);
    fSource_pXi_TM_Epos->SetParameter(0,1.0);
    fSource_pXi_TM_Epos->SetParLimits(0,0.5,2.0);
    hSource_pXi_TM_Epos->Fit(fSource_pXi_TM_Epos,"Q, S, N, R, M");


    TF1* fSource_pOmega_Back = new TF1(TString::Format("fSource_pOmega_Back"),GaussSourceTF1,rMin,rMax,1);
    fSource_pOmega_Back->SetParameter(0,1.0);
    fSource_pOmega_Back->SetParLimits(0,0.5,2.0);
    hSource_pOmega_Back->Fit(fSource_pOmega_Back,"Q, S, N, R, M");

    TF1* fSource_pOmega_Epos = new TF1(TString::Format("fSource_pOmega_Epos"),GaussSourceTF1,rMin,rMax,1);
    fSource_pOmega_Epos->SetParameter(0,1.0);
    fSource_pOmega_Epos->SetParLimits(0,0.5,2.0);
    hSource_pOmega_Epos->Fit(fSource_pOmega_Epos,"Q, S, N, R, M");

    TF1* fSource_pOmega_Rand = new TF1(TString::Format("fSource_pOmega_Rand"),GaussSourceTF1,rMin,rMax,1);
    fSource_pOmega_Rand->SetParameter(0,1.0);
    fSource_pOmega_Rand->SetParLimits(0,0.5,2.0);
    hSource_pOmega_Rand->Fit(fSource_pOmega_Rand,"Q, S, N, R, M");

    TF1* fSource_pOmega_TM_Back = new TF1(TString::Format("fSource_pOmega_TM_Back"),GaussSourceTF1,rMin,rMax,1);
    fSource_pOmega_TM_Back->SetParameter(0,1.0);
    fSource_pOmega_TM_Back->SetParLimits(0,0.5,2.0);
    hSource_pOmega_TM_Back->Fit(fSource_pOmega_TM_Back,"Q, S, N, R, M");

    TF1* fSource_pOmega_TM_Epos = new TF1(TString::Format("fSource_pOmega_TM_Epos"),GaussSourceTF1,rMin,rMax,1);
    fSource_pOmega_TM_Epos->SetParameter(0,1.0);
    fSource_pOmega_TM_Epos->SetParLimits(0,0.5,2.0);
    hSource_pOmega_TM_Epos->Fit(fSource_pOmega_TM_Epos,"Q, S, N, R, M");


    TF1* fSource_pp_Back = new TF1(TString::Format("fSource_pp_Back"),GaussSourceTF1,rMin,rMax,1);
    fSource_pp_Back->SetParameter(0,1.0);
    fSource_pp_Back->SetParLimits(0,0.5,2.0);
    hSource_pp_Back->Fit(fSource_pp_Back,"Q, S, N, R, M");

    TF1* fSource_pp_Epos = new TF1(TString::Format("fSource_pp_Epos"),GaussSourceTF1,rMin,rMax,1);
    fSource_pp_Epos->SetParameter(0,1.0);
    fSource_pp_Epos->SetParLimits(0,0.5,2.0);
    hSource_pp_Epos->Fit(fSource_pp_Epos,"Q, S, N, R, M");

    TF1* fSource_pp_Rand = new TF1(TString::Format("fSource_pp_Rand"),GaussSourceTF1,rMin,rMax,1);
    fSource_pp_Rand->SetParameter(0,1.0);
    fSource_pp_Rand->SetParLimits(0,0.5,2.0);
    hSource_pp_Rand->Fit(fSource_pp_Rand,"Q, S, N, R, M");

    TF1* fSource_pp_TM_Back = new TF1(TString::Format("fSource_pp_TM_Back"),GaussSourceTF1,rMin,rMax,1);
    fSource_pp_TM_Back->SetParameter(0,1.0);
    fSource_pp_TM_Back->SetParLimits(0,0.5,2.0);
    hSource_pp_TM_Back->Fit(fSource_pp_TM_Back,"Q, S, N, R, M");

    TF1* fSource_pp_TM_Epos = new TF1(TString::Format("fSource_pp_TM_Epos"),GaussSourceTF1,rMin,rMax,1);
    fSource_pp_TM_Epos->SetParameter(0,1.0);
    fSource_pp_TM_Epos->SetParLimits(0,0.5,2.0);
    hSource_pp_TM_Epos->Fit(fSource_pp_TM_Epos,"Q, S, N, R, M");

    TF1* fSource_pL_TM_Epos = new TF1(TString::Format("fSource_pL_TM_Epos"),GaussSourceTF1,rMin,rMax,1);
    fSource_pL_TM_Epos->SetParameter(0,1.0);
    fSource_pL_TM_Epos->SetParLimits(0,0.5,2.0);
    hSource_pL_TM_Epos->Fit(fSource_pL_TM_Epos,"Q, S, N, R, M");

    printf("Core vs effective sources:\n");
    printf(" pp (Back): %.3f --> %.3f fm\n",Core_pp_Back,fSource_pp_Back->GetParameter(0));
    printf(" pp (Epos): %.3f --> %.3f fm\n",Core_pp_TM_Epos,fSource_pp_TM_Epos->GetParameter(0));
    printf(" pΛ (Epos): %.3f --> %.3f fm\n",Core_pL_TM_Epos,fSource_pL_TM_Epos->GetParameter(0));
    printf(" pΞ (Back): %.3f --> %.3f fm\n",Core_pXi_Back,fSource_pXi_Back->GetParameter(0));
    printf(" pΞ (Epos): %.3f --> %.3f fm\n",Core_pXi_TM_Epos,fSource_pXi_TM_Epos->GetParameter(0));
    printf(" pΩ (Back): %.3f --> %.3f fm\n",Core_pOmega_Back,fSource_pOmega_Back->GetParameter(0));
    printf(" pΩ (Epos): %.3f --> %.3f fm\n",Core_pOmega_TM_Epos,fSource_pOmega_TM_Epos->GetParameter(0));


    gKitty_pXi_Back.Write();
    gKitty_pXi_Epos.Write();
    gKitty_pXi_Rand.Write();
    gKitty_pXi_TM_Back.Write();
    gKitty_pXi_TM_Epos.Write();
    gKitty_pOmega_Back.Write();
    gKitty_pOmega_Epos.Write();
    gKitty_pOmega_Rand.Write();
    gKitty_pOmega_TM_Back.Write();
    gKitty_pOmega_TM_Epos.Write();
    gKitty_pp_Back.Write();
    gKitty_pp_Epos.Write();
    gKitty_pp_Rand.Write();
    gKitty_pp_TM_Back.Write();
    gKitty_pp_TM_Epos.Write();
    gKitty_pL_TM_Epos.Write();

    gSource_pXi_Back.Write();
    gSource_pXi_Epos.Write();
    gSource_pXi_Rand.Write();
    gSource_pXi_TM_Back.Write();
    gSource_pXi_TM_Epos.Write();
    gSource_pOmega_Back.Write();
    gSource_pOmega_Epos.Write();
    gSource_pOmega_Rand.Write();
    gSource_pOmega_TM_Back.Write();
    gSource_pOmega_TM_Epos.Write();

    hSource_pXi_Back->Write();
    hSource_pXi_Epos->Write();
    hSource_pXi_Rand->Write();
    hSource_pXi_TM_Back->Write();
    hSource_pXi_TM_Epos->Write();
    hSource_pOmega_Back->Write();
    hSource_pOmega_Epos->Write();
    hSource_pOmega_Rand->Write();
    hSource_pOmega_TM_Back->Write();
    hSource_pOmega_TM_Epos->Write();
    hSource_pp_Back->Write();
    hSource_pp_Epos->Write();
    hSource_pp_Rand->Write();
    hSource_pp_TM_Back->Write();
    hSource_pp_TM_Epos->Write();
    hSource_pL_TM_Epos->Write();

    fSource_pXi_Back->Write();
    fSource_pXi_Epos->Write();
    fSource_pXi_Rand->Write();
    fSource_pXi_TM_Back->Write();
    fSource_pXi_TM_Epos->Write();
    fSource_pOmega_Back->Write();
    fSource_pOmega_Epos->Write();
    fSource_pOmega_Rand->Write();
    fSource_pOmega_TM_Back->Write();
    fSource_pOmega_TM_Epos->Write();
    fSource_pp_Back->Write();
    fSource_pp_Epos->Write();
    fSource_pp_Rand->Write();
    fSource_pp_TM_Back->Write();
    fSource_pp_TM_Epos->Write();
    fSource_pL_TM_Epos->Write();

    delete hSource_pXi_Back;
    delete hSource_pXi_Epos;
    delete hSource_pXi_Rand;
    delete hSource_pXi_TM_Back;
    delete hSource_pXi_TM_Epos;
    delete hSource_pOmega_Back;
    delete hSource_pOmega_Epos;
    delete hSource_pOmega_Rand;
    delete hSource_pOmega_TM_Back;
    delete hSource_pOmega_TM_Epos;
    delete hSource_pp_Back;
    delete hSource_pp_Epos;
    delete hSource_pp_Rand;
    delete hSource_pp_TM_Back;
    delete hSource_pp_TM_Epos;
    delete hSource_pL_TM_Epos;
    delete fSource_pXi_Back;
    delete fSource_pXi_Epos;
    delete fSource_pXi_Rand;
    delete fSource_pXi_TM_Back;
    delete fSource_pXi_TM_Epos;
    delete fSource_pOmega_Back;
    delete fSource_pOmega_Epos;
    delete fSource_pOmega_Rand;
    delete fSource_pOmega_TM_Back;
    delete fSource_pOmega_TM_Epos;
    delete fSource_pp_Back;
    delete fSource_pp_Epos;
    delete fSource_pp_Rand;
    delete fSource_pp_TM_Back;
    delete fSource_pp_TM_Epos;
    delete fSource_pL_TM_Epos;
    delete fOut;
}


void GenerateFakeProtonLambda(){

    const double SourceSize = 1.5;
    const double LambdaPar = 0.6;
    const double Norm = 1.05;

    const unsigned NumMomBins = 80;
    const double kMin = 0;
    const double kMax = 320;

    DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder("/home/dmihaylov/CernBox/CatsFiles");

    CATS Kitty;
    Kitty.SetMomBins(NumMomBins,kMin,kMax);
    AnalysisObject.SetUpCats_pL(Kitty,"Usmani","Gauss");
    Kitty.SetAnaSource(0,SourceSize);
    Kitty.KillTheCat();

    TH1F* hDummyProtonLambda = new TH1F("hDummyProtonLambda","hDummyProtonLambda",NumMomBins,kMin,kMax);

    TRandom3 rangen(8);

    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        double CkVal = Kitty.GetCorrFun(uBin);
        double Momentum = Kitty.GetMomentum(uBin);
        double CkErr = 0.0001*pow(400.-Momentum,1.1);
        hDummyProtonLambda->SetBinContent(uBin+1,Norm*(rangen.Gaus(CkVal,CkErr)*LambdaPar+1.-LambdaPar));
        hDummyProtonLambda->SetBinError(uBin+1,CkErr);
    }

    TFile OutFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CATS_TUTORIAL_2019/Files/DummyProtonLambda2.root","recreate");

    hDummyProtonLambda->Write();

    delete hDummyProtonLambda;
}

void pp_in_txtfile(){

    const TString OutputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pp_in_txtfile/";
    const double RadiusSize = 1.2;
    const double WF_Mom = 100;

    const double kMin = 0;
    const double kMax = 300;
    const unsigned NumMomBins = 150;

    const double rMin = 0;
    const double rMax = 16;
    const unsigned NumRadBins = 320;

    DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder("/home/dmihaylov/CernBox/CatsFiles");

    CATS Kitty_SI;
    Kitty_SI.SetMaxNumThreads(4);
    Kitty_SI.SetMomBins(NumMomBins,kMin,kMax);
    AnalysisObject.SetUpCats_pp(Kitty_SI,"AV18","Gauss");
    Kitty_SI.SetAnaSource(0,1.2);
    Kitty_SI.SetEpsilonProp(2e-8);
    Kitty_SI.SetQ1Q2(0);
    Kitty_SI.KillTheCat();

    CATS Kitty_C;
    Kitty_C.SetMaxNumThreads(4);
    Kitty_C.SetMomBins(NumMomBins,kMin,kMax);
    AnalysisObject.SetUpCats_pp(Kitty_C,"AV18","Gauss");
    Kitty_C.SetAnaSource(0,1.2);
    Kitty_C.SetEpsilonProp(2e-8);
    Kitty_C.SetQ1Q2(1);
    Kitty_C.RemoveShortRangePotential(0,0);
    Kitty_C.RemoveShortRangePotential(0,2);
    Kitty_C.RemoveShortRangePotential(1,0);
    Kitty_C.RemoveShortRangePotential(1,1);
    Kitty_C.RemoveShortRangePotential(1,2);
    Kitty_C.KillTheCat();

    CATS Kitty_Full;
    Kitty_Full.SetMaxNumThreads(4);
    Kitty_Full.SetMomBins(NumMomBins,kMin,kMax);
    AnalysisObject.SetUpCats_pp(Kitty_Full,"AV18","Gauss");
    Kitty_Full.SetAnaSource(0,1.2);
    Kitty_Full.SetEpsilonProp(2e-8);
    Kitty_Full.KillTheCat();

    TH1F* hV_AV18 = new TH1F("hV_AV18","hV_AV18",NumRadBins,rMin,rMax);
    TH1F* hV_C = new TH1F("hV_C","hV_C",NumRadBins,rMin,rMax);
    TH1F* hV_Tot = new TH1F("hV_Tot","hV_Tot",NumRadBins,rMin,rMax);
    TH1F* hWF_AV18_1S0 = new TH1F("hWF_AV18_1S0","hWF_AV18_1S0",NumRadBins,rMin,rMax);
    TH1F* hWF_C_1S0 = new TH1F("hWF_C_1S0","hWF_C_1S0",NumRadBins,rMin,rMax);
    TH1F* hWF_Full_1S0 = new TH1F("hWF_Full_1S0","hWF_Full_1S0",NumRadBins,rMin,rMax);
    TH1F* hWFu_AV18_1S0 = new TH1F("hWFu_AV18_1S0","hWFu_AV18_1S0",NumRadBins,rMin,rMax);
    TH1F* hWFu_C_1S0 = new TH1F("hWFu_C_1S0","hWFu_C_1S0",NumRadBins,rMin,rMax);
    TH1F* hWFu_Full_1S0 = new TH1F("hWFu_Full_1S0","hWFu_Full_1S0",NumRadBins,rMin,rMax);
    TH1F* hWF_AV18 = new TH1F("hWF_AV18","hWF_AV18",NumRadBins,rMin,rMax);
    TH1F* hWF_C = new TH1F("hWF_C","hWF_C",NumRadBins,rMin,rMax);
    TH1F* hWF_Full = new TH1F("hWF_Full","hWF_Full",NumRadBins,rMin,rMax);
    TH1F* hCk_AV18 = new TH1F("hCk_AV18","hCk_AV18",NumMomBins,kMin,kMax);
    TH1F* hCk_C = new TH1F("hCk_C","hCk_C",NumMomBins,kMin,kMax);
    TH1F* hCk_Full = new TH1F("hCk_Full","hCk_Full",NumMomBins,kMin,kMax);

    //AV18 potential
    //Coulomb potential
    //Wave function at 30 MeV
    FILE * ppFilePotentials;
    ppFilePotentials = fopen (TString::Format("%sppFilePotentials.txt",OutputFolder.Data()),"w");
    fprintf (ppFilePotentials, "%16s","r (fm)");
    fprintf (ppFilePotentials, "%16s","V_AV18(r)");
    fprintf (ppFilePotentials, "%16s","V_C(r)");
    fprintf (ppFilePotentials, "%16s","V_Tot(r)");
    fprintf (ppFilePotentials, "%16s","Re_u1S0_AV18(r)");
    fprintf (ppFilePotentials, "%16s","Re_u1S0_C(r)");
    fprintf (ppFilePotentials, "%16s","Re_u1S0_Tot(r)");
    fprintf (ppFilePotentials, "%17s","Re_ψ1S0_AV18(r)");
    fprintf (ppFilePotentials, "%17s","Re_ψ1S0_C(r)");
    fprintf (ppFilePotentials, "%17s","Re_ψ1S0_Tot(r)");
    fprintf (ppFilePotentials, "%17s","|ψ|^2_AV18(r)");
    fprintf (ppFilePotentials, "%17s","|ψ|^2_C(r)");
    fprintf (ppFilePotentials, "%17s","|ψ|^2_Tot(r)");
    fprintf (ppFilePotentials,"\n");
    for(unsigned uRad=0; uRad<NumRadBins; uRad++){
        double RAD = hV_AV18->GetBinCenter(uRad+1);
        double V_AV18 = Kitty_SI.EvaluateThePotential(0,0,WF_Mom,RAD);
        double V_C = Kitty_Full.EvaluateCoulombPotential(RAD);
        double V_Tot = V_AV18+V_C;
        hV_AV18->SetBinContent(uRad+1,V_AV18);
        hV_C->SetBinContent(uRad+1,V_C);
        hV_Tot->SetBinContent(uRad+1,V_Tot);
        unsigned WhichMomBin = Kitty_SI.GetMomBin(WF_Mom);
        double WF_AV18_1S0 = real(Kitty_SI.EvalRadialWaveFunction(WhichMomBin,0,0,RAD,true));
        double WF_C_1S0 = real(Kitty_C.EvalRadialWaveFunction(WhichMomBin,0,0,RAD,true));
        double WF_Full_1S0 = real(Kitty_Full.EvalRadialWaveFunction(WhichMomBin,0,0,RAD,true));
        double WFu_AV18_1S0 = real(Kitty_SI.EvalRadialWaveFunction(WhichMomBin,0,0,RAD,false));
        double WFu_C_1S0 = real(Kitty_C.EvalRadialWaveFunction(WhichMomBin,0,0,RAD,false));
        double WFu_Full_1S0 = real(Kitty_Full.EvalRadialWaveFunction(WhichMomBin,0,0,RAD,false));
        double WF_AV18 = Kitty_SI.EvalWaveFun2(WhichMomBin,RAD,0);
        double WF_C = Kitty_C.EvalWaveFun2(WhichMomBin,RAD,0);
        double WF_Full = Kitty_Full.EvalWaveFun2(WhichMomBin,RAD,0);
        hWF_AV18_1S0->SetBinContent(uRad+1,WF_AV18_1S0);
        hWF_C_1S0->SetBinContent(uRad+1,WF_C_1S0);
        hWF_Full_1S0->SetBinContent(uRad+1,WF_Full_1S0);
        hWFu_AV18_1S0->SetBinContent(uRad+1,WFu_AV18_1S0);
        hWFu_C_1S0->SetBinContent(uRad+1,WFu_C_1S0);
        hWFu_Full_1S0->SetBinContent(uRad+1,WFu_Full_1S0);
        hWF_AV18->SetBinContent(uRad+1,WF_AV18);
        hWF_C->SetBinContent(uRad+1,WF_C);
        hWF_Full->SetBinContent(uRad+1,WF_Full);
        fprintf (ppFilePotentials, "%16.4e",RAD);
        fprintf (ppFilePotentials, "%16.4e",V_AV18);
        fprintf (ppFilePotentials, "%16.4e",V_C);
        fprintf (ppFilePotentials, "%16.4e",V_Tot);
        fprintf (ppFilePotentials, "%16.4e",WFu_AV18_1S0);
        fprintf (ppFilePotentials, "%16.4e",WFu_C_1S0);
        fprintf (ppFilePotentials, "%16.4e",WFu_Full_1S0);
        fprintf (ppFilePotentials, "%16.4e",WF_AV18_1S0);
        fprintf (ppFilePotentials, "%16.4e",WF_C_1S0);
        fprintf (ppFilePotentials, "%16.4e",WF_Full_1S0);
        fprintf (ppFilePotentials, "%16.4e",WF_AV18);
        fprintf (ppFilePotentials, "%16.4e",WF_C);
        fprintf (ppFilePotentials, "%16.4e",WF_Full);
        fprintf (ppFilePotentials,"\n");
    }

    FILE * ppFileCk;
    ppFileCk = fopen (TString::Format("%sppFileCk.txt",OutputFolder.Data()),"w");
    fprintf (ppFileCk, "%16s","k (MeV)");
    fprintf (ppFileCk, "%16s","Ck_AV18");
    fprintf (ppFileCk, "%16s","Ck_C");
    fprintf (ppFileCk, "%16s","Ck_Full");
    fprintf (ppFileCk,"\n");
    //Ck (full)
    for(unsigned uMom=0; uMom<NumMomBins; uMom++){
        double MOM = Kitty_Full.GetMomentum(uMom);
        double Ck_AV18 = Kitty_SI.GetCorrFun(uMom);
        double Ck_C = Kitty_C.GetCorrFun(uMom);
        double Ck_Full = Kitty_Full.GetCorrFun(uMom);
        fprintf (ppFileCk, "%16.4e",MOM);
        fprintf (ppFileCk, "%16.4e",Ck_AV18);
        fprintf (ppFileCk, "%16.4e",Ck_C);
        fprintf (ppFileCk, "%16.4e",Ck_Full);
        fprintf (ppFileCk,"\n");
        hCk_AV18->SetBinContent(uMom+1,Ck_AV18);
        hCk_C->SetBinContent(uMom+1,Ck_C);
        hCk_Full->SetBinContent(uMom+1,Ck_Full);
    }

    TFile* fOutput = new TFile(OutputFolder+"fOutput_pp.root","recreate");
    hV_AV18->Write();
    hV_C->Write();
    hV_Tot->Write();
    hWFu_AV18_1S0->Write();
    hWFu_C_1S0->Write();
    hWFu_Full_1S0->Write();
    hWF_AV18_1S0->Write();
    hWF_C_1S0->Write();
    hWF_Full_1S0->Write();
    hWF_AV18->Write();
    hWF_C->Write();
    hWF_Full->Write();
    hCk_AV18->Write();
    hCk_C->Write();
    hCk_Full->Write();

    delete hV_AV18;
    delete hV_C;
    delete hV_Tot;
    delete hWFu_AV18_1S0;
    delete hWFu_C_1S0;
    delete hWFu_Full_1S0;
    delete hWF_AV18_1S0;
    delete hWF_C_1S0;
    delete hWF_Full_1S0;
    delete hWF_AV18;
    delete hWF_C;
    delete hWF_Full;
    delete hCk_AV18;
    delete hCk_C;
    delete hCk_Full;
    delete fOutput;
}

void pp_pLambda_pXi_Ratios(){
    const TString OutputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pp_pLambda_pXi_Ratios/";
    const TString AnaResultFile = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample10HM/AnalysisResults.root";
    const unsigned NumMomBinsFine = 80;
    const double kMinFine = 0;
    const double kMaxFine = 320;
    const double kStepFine = (kMaxFine-kMinFine)/(NumMomBinsFine);
    const unsigned NumMomBinsCoarse = 8;
    const double kMaxCoarse = 400;
    const double kStepCoarse = (kMaxCoarse-kMaxFine)/(NumMomBinsCoarse);
    const unsigned NumMomBins = NumMomBinsFine + NumMomBinsCoarse;
    const unsigned NumSyst = 3;
    enum Systems { pp, pL, pXim };
    TString SysName[NumSyst];
    SysName[0] = "pp";
    SysName[1] = "pL";
    SysName[2] = "pXim";

    double* MomBins = new double [NumMomBins+1];
    MomBins[0] = kMinFine;
    for(unsigned uBin=1; uBin<=NumMomBins; uBin++){
        if(uBin<=NumMomBinsFine){
            MomBins[uBin] = MomBins[uBin-1]+kStepFine;
        }
        else{
            MomBins[uBin] = MomBins[uBin-1]+kStepCoarse;
        }
        //printf("MomBins[%u]=%.1f\n",uBin,MomBins[uBin]);
    }

    DLM_CommonAnaFunctions AnaObject;

    CATS Kitty[NumSyst];

    Kitty[pp].SetMomBins(NumMomBins,MomBins);
    AnaObject.SetUpCats_pp(Kitty[pp],"AV18","Gauss");
    Kitty[pp].SetAnaSource(0,1.2);
    Kitty[pp].SetMaxNumThreads(4);

    Kitty[pL].SetMomBins(NumMomBins,MomBins);
    AnaObject.SetUpCats_pL(Kitty[pL],"NLO_Coupled_S","Gauss");//NLO_Coupled_S
    Kitty[pL].SetAnaSource(0,1.4);
    Kitty[pL].SetMaxNumThreads(4);

    Kitty[pXim].SetMomBins(NumMomBins,MomBins);
    AnaObject.SetUpCats_pXim(Kitty[pXim],"pXim_HALQCD1","Gauss");
    Kitty[pXim].SetAnaSource(0,1.0);
    Kitty[pXim].SetMaxNumThreads(4);

    Kitty[pp].KillTheCat();
    Kitty[pL].KillTheCat();
    Kitty[pXim].KillTheCat();


    TGraph graph[NumSyst];
    for(unsigned uSys=0; uSys<NumSyst; uSys++){
        graph[uSys].SetName("graph_"+SysName[uSys]);
        graph[uSys].Set(NumMomBins);
        for(unsigned uBin=0; uBin<NumMomBins; uBin++){
            graph[uSys].SetPoint(uBin,Kitty[uSys].GetMomentum(uBin),Kitty[uSys].GetCorrFun(uBin));
        }
    }

    const unsigned NumRatioPlots = NumSyst*(NumSyst-1)/2;
    TGraph graphRatio[NumRatioPlots];
    double CkVal1,CkVal2,MomVal1,MomVal2;
    unsigned uRat=0;
    for(unsigned uSys1=0; uSys1<NumSyst; uSys1++){
        for(unsigned uSys2=uSys1+1; uSys2<NumSyst; uSys2++){
            graphRatio[uRat].SetName("graphRatio_"+SysName[uSys1]+"_"+SysName[uSys2]);
            graphRatio[uRat].Set(NumMomBins);
            for(unsigned uBin=0; uBin<NumMomBins; uBin++){
                graph[uSys1].GetPoint(uBin,MomVal1,CkVal1);
                graph[uSys2].GetPoint(uBin,MomVal2,CkVal2);
                graphRatio[uRat].SetPoint(uBin,MomVal1,CkVal1/CkVal2);
                if(MomVal1!=MomVal2){
                    printf("TROUBLE!!!!!\n");
                }
            }
            uRat++;
        }
    }



    //exp data
    char prefix[16];
    char addon[16];
    strcpy(prefix,"HM");
    strcpy(addon,"");

    TFile fInput(AnaResultFile,"read");
    TDirectoryFile *dirResults=(TDirectoryFile*)(fInput.FindObjectAny(Form("%sResults%s", prefix, addon)));
    TList *Results;
    dirResults->GetObject(Form("%sResults%s", prefix, addon),Results);

    TList* tmpFolder1;
    TList* tmpFolder2;
    const unsigned NumSpecies = 3;
    TString SpeciesDescr[NumSpecies];
    SpeciesDescr[0] = "p";
    SpeciesDescr[1] = "Lam";
    SpeciesDescr[2] = "Xim";
    const unsigned NumCkCombo = NumSpecies*(NumSpecies+1)/2;
    TString ComboDescr[NumCkCombo];
    TH1F** histSE = new TH1F* [NumSpecies];
    TH1F** histME = new TH1F* [NumSpecies];
    TH1F** histCk = new TH1F* [NumSpecies];

    TH1F** histSE_PP = new TH1F* [NumSpecies];
    TH1F** histME_PP = new TH1F* [NumSpecies];
    TH1F** histCk_PP = new TH1F* [NumSpecies];

    TH1F** histSE_APAP = new TH1F* [NumSpecies];
    TH1F** histME_APAP = new TH1F* [NumSpecies];
    TH1F** histCk_APAP = new TH1F* [NumSpecies];

    const unsigned NumCkRatios = NumCkCombo*(NumCkCombo-1)/2;
    TH1F** histCkRatio = new TH1F* [NumCkRatios];
    unsigned uCombo=0;
    TFile fOutput(OutputFolder+"fOutput.root","recreate");
    for(unsigned uSpec1=0; uSpec1<NumSpecies; uSpec1++){
        for(unsigned uSpec2=uSpec1; uSpec2<NumSpecies; uSpec2++){
            fInput.cd();
            //printf("\n%u %u\n",uSpec1,uSpec2);
            tmpFolder1=(TList*)Results->FindObject(TString::Format("Particle%u_Particle%u",uSpec1*2,uSpec2*2));
            //printf(" tmpFolder1=%p\n",tmpFolder1);
            TH1F* tmpSE1 = NULL; tmpSE1 = (TH1F*)tmpFolder1->FindObject(TString::Format("SEDist_Particle%u_Particle%u",uSpec1*2,uSpec2*2));
            //printf(" tmpSE1=%p\n",tmpSE1);
            //printf(" tmpSE1->GetNbinsX()=%i\n",tmpSE1->GetNbinsX());
            fOutput.cd();
            histSE[uCombo] = new TH1F(TString::Format("histSE_Particle%u_Particle%u",uSpec1*2,uSpec2*2),TString::Format("histSE_Particle%u_Particle%u",uSpec1*2,uSpec2*2),
                                      tmpSE1->GetNbinsX(),tmpSE1->GetBinLowEdge(1),tmpSE1->GetXaxis()->GetBinUpEdge(tmpSE1->GetNbinsX()));
//            histSE_PP[uCombo] = new TH1F(TString::Format("histSE_PP_%u_%u",uSpec1*2,uSpec2*2),TString::Format("histSE_PP_%u_%u",uSpec1*2,uSpec2*2),
//                                      tmpSE1->GetNbinsX(),tmpSE1->GetBinLowEdge(1),tmpSE1->GetXaxis()->GetBinUpEdge(tmpSE1->GetNbinsX()));
            fInput.cd();
//            histSE_APAP[uCombo] = new TH1F(TString::Format("histSE_APAP_Particle%u_Particle%u",uSpec1*2,uSpec2*2),TString::Format("histSE_APAP_Particle%u_Particle%u",uSpec1*2,uSpec2*2),
//                                      tmpSE1->GetNbinsX(),tmpSE1->GetBinLowEdge(1),tmpSE1->GetXaxis()->GetBinUpEdge(tmpSE1->GetNbinsX()));
            for(unsigned uBin=0; uBin<=tmpSE1->GetNbinsX()+1; uBin++){
                histSE[uCombo]->SetBinContent(uBin,tmpSE1->GetBinContent(uBin));
                histSE[uCombo]->SetBinError(uBin,tmpSE1->GetBinError(uBin));
//                histSE_PP[uCombo]->SetBinContent(uBin,tmpSE1->GetBinContent(uBin));
//                histSE_PP[uCombo]->SetBinError(uBin,tmpSE1->GetBinError(uBin));
            }
            //histSE[uCombo] = (TH1F*)tmpSE1->Clone(TString::Format("histSE_Particle%u_Particle%u",uSpec1*2,uSpec2*2));
            //printf("Done\n");
            TH1F* tmpME1 = NULL; tmpME1 = (TH1F*)tmpFolder1->FindObject(TString::Format("MEDist_Particle%u_Particle%u",uSpec1*2,uSpec2*2));
            //printf(" tmpME1=%p\n",tmpME1);
            fOutput.cd();
            histME[uCombo] = new TH1F(TString::Format("histME_Particle%u_Particle%u",uSpec1*2,uSpec2*2),TString::Format("histME_Particle%u_Particle%u",uSpec1*2,uSpec2*2),
                                      tmpME1->GetNbinsX(),tmpME1->GetBinLowEdge(1),tmpME1->GetXaxis()->GetBinUpEdge(tmpME1->GetNbinsX()));
//            histME_PP[uCombo] = new TH1F(TString::Format("histME_PP_Particle%u_Particle%u",uSpec1*2,uSpec2*2),TString::Format("histME_PP_Particle%u_Particle%u",uSpec1*2,uSpec2*2),
//                                      tmpME1->GetNbinsX(),tmpME1->GetBinLowEdge(1),tmpME1->GetXaxis()->GetBinUpEdge(tmpME1->GetNbinsX()));
//            histME_APAP[uCombo] = new TH1F(TString::Format("histME_APAP_Particle%u_Particle%u",uSpec1*2,uSpec2*2),TString::Format("histME_APAP_Particle%u_Particle%u",uSpec1*2,uSpec2*2),
//                                      tmpME1->GetNbinsX(),tmpME1->GetBinLowEdge(1),tmpME1->GetXaxis()->GetBinUpEdge(tmpME1->GetNbinsX()));
            fInput.cd();
            for(unsigned uBin=0; uBin<=tmpME1->GetNbinsX()+1; uBin++){
                histME[uCombo]->SetBinContent(uBin,tmpME1->GetBinContent(uBin));
                histME[uCombo]->SetBinError(uBin,tmpME1->GetBinError(uBin));
//                histME_PP[uCombo]->SetBinContent(uBin,tmpME1->GetBinContent(uBin));
//                histME_PP[uCombo]->SetBinError(uBin,tmpME1->GetBinError(uBin));
            }
            //histME[uCombo] = (TH1F*)tmpME1->Clone(TString::Format("histME_Particle%u_Particle%u",uSpec1*2,uSpec2*2));
            tmpFolder2=(TList*)Results->FindObject(TString::Format("Particle%u_Particle%u",uSpec1*2+1,uSpec2*2+1));
            //printf(" tmpFolder2=%p\n",tmpFolder2);
            TH1F* tmpSE2 = NULL; tmpSE2 = (TH1F*)tmpFolder2->FindObject(TString::Format("SEDist_Particle%u_Particle%u",uSpec1*2+1,uSpec2*2+1));
            //printf(" tmpSE2=%p\n",tmpSE2);
            TH1F* tmpME2 = NULL; tmpME2 = (TH1F*)tmpFolder2->FindObject(TString::Format("MEDist_Particle%u_Particle%u",uSpec1*2+1,uSpec2*2+1));
            //printf(" tmpME2=%p\n",tmpME2);
            histSE[uCombo]->Add(tmpSE2);
            histME[uCombo]->Add(tmpME2);
            histME[uCombo]->Scale(histSE[uCombo]->Integral(0,histSE[uCombo]->GetNbinsX()+1)/histME[uCombo]->Integral(0,histME[uCombo]->GetNbinsX()+1));

            for(unsigned uBin=0; uBin<=tmpSE2->GetNbinsX()+1; uBin++){
//                histSE_APAP[uCombo]->SetBinContent(uBin,tmpSE2->GetBinContent(uBin));
//                histSE_APAP[uCombo]->SetBinError(uBin,tmpSE2->GetBinError(uBin));
            }
            for(unsigned uBin=0; uBin<=tmpME2->GetNbinsX()+1; uBin++){
//                histME_APAP[uCombo]->SetBinContent(uBin,tmpME2->GetBinContent(uBin));
//                histME_APAP[uCombo]->SetBinError(uBin,tmpME2->GetBinError(uBin));
            }

            //printf("hSE = %e\n",histSE[uCombo]->Integral(0,histSE[uCombo]->GetNbinsX()+1));
            //printf("hME = %e\n",histME[uCombo]->Integral(0,histME[uCombo]->GetNbinsX()+1));
            fOutput.cd();
            histCk[uCombo] = new TH1F(TString::Format("histCk_Particle%u_Particle%u",uSpec1*2,uSpec2*2),TString::Format("histCk_Particle%u_Particle%u",uSpec1*2,uSpec2*2),
                                      histSE[uCombo]->GetNbinsX(),histSE[uCombo]->GetBinLowEdge(1),histSE[uCombo]->GetXaxis()->GetBinUpEdge(histSE[uCombo]->GetNbinsX()));
//            histCk_PP[uCombo] = new TH1F(TString::Format("histCk_PP_Particle%u_Particle%u",uSpec1*2,uSpec2*2),TString::Format("histCk_PP_Particle%u_Particle%u",uSpec1*2,uSpec2*2),
//                                      histSE[uCombo]->GetNbinsX(),histSE[uCombo]->GetBinLowEdge(1),histSE[uCombo]->GetXaxis()->GetBinUpEdge(histSE[uCombo]->GetNbinsX()));
//            histCk_APAP[uCombo] = new TH1F(TString::Format("histCk_APAP_Particle%u_Particle%u",uSpec1*2,uSpec2*2),TString::Format("histCk_APAP_Particle%u_Particle%u",uSpec1*2,uSpec2*2),
//                                      histSE[uCombo]->GetNbinsX(),histSE[uCombo]->GetBinLowEdge(1),histSE[uCombo]->GetXaxis()->GetBinUpEdge(histSE[uCombo]->GetNbinsX()));
            fInput.cd();
            for(unsigned uBin=0; uBin<=histSE[uCombo]->GetNbinsX()+1; uBin++){
                histCk[uCombo]->SetBinContent(uBin,histSE[uCombo]->GetBinContent(uBin));
                histCk[uCombo]->SetBinError(uBin,histSE[uCombo]->GetBinError(uBin));

//                histCk_PP[uCombo]->SetBinContent(uBin,histSE_PP[uCombo]->GetBinContent(uBin));
//                histCk_PP[uCombo]->SetBinError(uBin,histSE_PP[uCombo]->GetBinError(uBin));

//                histCk_APAP[uCombo]->SetBinContent(uBin,histSE_APAP[uCombo]->GetBinContent(uBin));
//                histCk_APAP[uCombo]->SetBinError(uBin,histSE_APAP[uCombo]->GetBinError(uBin));
            }
            //histCk[uCombo] = (TH1F*)histSE[uCombo]->Clone(TString::Format("histCk_Particle%u_Particle%u",uSpec1*2,uSpec2*2));
            histCk[uCombo]->Divide(histME[uCombo]);
//            histCk_PP[uCombo]->Divide(histME_PP[uCombo]);
//            histCk_APAP[uCombo]->Divide(histME_APAP[uCombo]);
            ComboDescr[uCombo] = SpeciesDescr[uSpec1]+SpeciesDescr[uSpec2];
            uCombo++;
        }
    }
printf("Hello\n");
    unsigned uRatio=0;
    fOutput.cd();
    for(unsigned uCombo1=0; uCombo1<NumCkCombo; uCombo1++){
        for(unsigned uCombo2=uCombo1+1; uCombo2<NumCkCombo; uCombo2++){
            histCkRatio[uRatio] = new TH1F("histRatio_"+ComboDescr[uCombo1]+"_"+ComboDescr[uCombo2],"histRatio_"+ComboDescr[uCombo1]+"_"+ComboDescr[uCombo2],
                                      histCk[uCombo1]->GetNbinsX(),histCk[uCombo1]->GetBinLowEdge(1),histCk[uCombo1]->GetXaxis()->GetBinUpEdge(histCk[uCombo1]->GetNbinsX()));
            for(unsigned uBin=0; uBin<=histCk[uCombo1]->GetNbinsX()+1; uBin++){
                histCkRatio[uRatio]->SetBinContent(uBin,histCk[uCombo1]->GetBinContent(uBin));
                histCkRatio[uRatio]->SetBinError(uBin,histCk[uCombo1]->GetBinError(uBin));
            }
            //histCkRatio[uRatio] = (TH1F*)histCk[uCombo1]->Clone("histRatio_"+ComboDescr[uCombo1]+"_"+ComboDescr[uCombo2]);
            histCkRatio[uRatio]->Divide(histCk[uCombo2]);
            uRatio++;
        }
    }
printf("Hello 0\n");

    for(unsigned uSys=0; uSys<NumSyst; uSys++){
//        graph[uSys].Write();
    }
printf("Hello 1\n");
    for(uRat=0; uRat<NumRatioPlots; uRat++){
//        graphRatio[uRat].Write();
    }
printf("Hello 2\n");
    for(uCombo=0; uCombo<NumCkCombo; uCombo++){
        histSE[uCombo]->Write();
//        histSE_PP[uCombo]->Write();
//        histSE_APAP[uCombo]->Write();
//        histME[uCombo]->Write();
//        histME_PP[uCombo]->Write();
//        histME_APAP[uCombo]->Write();
//        histCk[uCombo]->Write();
//        histCk_PP[uCombo]->Write();
//        histCk_APAP[uCombo]->Write();
    }
printf("Hello 3\n");
    for(uRatio=0; uRatio<NumCkRatios; uRatio++){
//        histCkRatio[uRatio]->Write();
    }
printf("Hello 4\n");

    delete [] MomBins;
    for(uCombo=0; uCombo<NumCkCombo; uCombo++){
//        delete histCk[uCombo];
//        delete histCk_PP[uCombo];
//        delete histCk_APAP[uCombo];
    }
//    delete [] histCk;
//    delete [] histCk_PP;
//    delete [] histCk_APAP;

    for(uRatio=0; uRatio<NumCkRatios; uRatio++){
//        delete histCkRatio[uRatio];
    }
//    delete [] histCkRatio;
}

void TestRandom(){
    TRandom3 Ran3(11);
    DLM_Random RanD(11);
    TH1F* hist3 = new TH1F("hist3","hist3",100,0,10);
    TH1F* histD = new TH1F("histD","histD",100,0,10);

    for(unsigned uIter=0; uIter<100000; uIter++){
        double rv3 = Ran3.Exp(1./3.);
        double rvd = RanD.Exponential(3);
        //printf("rv = %f (%f)\n",rv3,rvd);
        hist3->Fill(rv3);
        histD->Fill(rvd);
    }
    TFile fOut("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Temp/fRandom.root","recreate");
    hist3->Write();
    histD->Write();
}

void SmearTest(){
/*
    CATS Kitty;
    Kitty.SetMaxNumThreads(1);
    Kitty.SetMomBins(30,0,120);
    AnalysisObject.SetUpCats_pp(Kitty,"AV18","Gauss");
    Kitty.SetAnaSource(0,1.0);
    Kitty.SetAnaSource(1,2);
    Kitty.SetEpsilonProp(1e-8);
    Kitty.KillTheCat();

    DLM_Ck* Ck_pp = new DLM_Ck(Kitty.GetNumSourcePars(),0,Kitty);
    Ck_pp->Update();

    TH2F* hResolution_pp = AnalysisObject.GetResolutionMatrix(DataSample,"pp");
    DLM_CkDecomposition CkDec_pp("pp",0,*Ck_pp,hResolution_pp);

    CkDec_pp.Update();
*/
}

void TestCauchy(){

    DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder("/home/dmihaylov/CernBox/CatsFiles");

    CATS AB_Cauchy;
    AB_Cauchy.SetMomBins(25,0,100);
    AnalysisObject.SetUpCats_pp(AB_Cauchy,"AV18","Cauchy",0,0);

    CATS AB_McCauchy;
    AB_McCauchy.SetMomBins(25,0,100);
    AnalysisObject.SetUpCats_pp(AB_McCauchy,"AV18","McLevy_ResoTM",0,0);

    CATS AB_McAlmostCauchy;
    AB_McAlmostCauchy.SetMomBins(25,0,100);
    AnalysisObject.SetUpCats_pp(AB_McAlmostCauchy,"AV18","McLevy_ResoTM",0,0);

    AB_Cauchy.SetAnaSource(0,1.0);

    AB_McCauchy.SetAnaSource(0,1.0);
    AB_McCauchy.SetAnaSource(1,1.0);

    AB_McAlmostCauchy.SetAnaSource(0,1.0);
    AB_McAlmostCauchy.SetAnaSource(1,1.01);

    unsigned NumSourcePts = 512;
    double rMin = 0;
    double rMax = 16;
    double RAD;

    TFile fOut("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/TestCauchy/fOut.root","recreate");

    TH1F* hRad_Cauchy = new TH1F("hRad_Cauchy","hRad_Cauchy",NumSourcePts,rMin,rMax);
    TH1F* hRad_McCauchy = new TH1F("hRad_McCauchy","hRad_McCauchy",NumSourcePts,rMin,rMax);
    TH1F* hRad_McAlmostCauchy = new TH1F("hRad_McAlmostCauchy","hRad_McAlmostCauchy",NumSourcePts,rMin,rMax);

    for(unsigned uRad=0; uRad<NumSourcePts; uRad++){
        RAD = hRad_Cauchy->GetBinCenter(uRad+1);
        hRad_Cauchy->SetBinContent(uRad+1,AB_Cauchy.EvaluateTheSource(10,RAD,0));
        hRad_McCauchy->SetBinContent(uRad+1,AB_McCauchy.EvaluateTheSource(10,RAD,0));
        hRad_McAlmostCauchy->SetBinContent(uRad+1,AB_McAlmostCauchy.EvaluateTheSource(10,RAD,0));
    }
    hRad_Cauchy->Write();
    hRad_McCauchy->Write();
    hRad_McAlmostCauchy->Write();

    delete hRad_Cauchy;
    delete hRad_McCauchy;
    delete hRad_McAlmostCauchy;
}


void Test_New_pLambda(){


    DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder("/home/dmihaylov/CernBox/CatsFiles");
    const double SourceSize = 1.125;
    const double LambdaPar = 0.4713;

    CATS KittyNLO_Old;
    KittyNLO_Old.SetMaxNumThreads(1);
    KittyNLO_Old.SetMomBins(85,0,340);
    AnalysisObject.SetUpCats_pL(KittyNLO_Old,"NLO_Coupled_S","Gauss");
//KittyNLO_Old.SetChannelWeight(2,0.25*2./3.);
//KittyNLO_Old.SetChannelWeight(3,0.75*2./3.);
    KittyNLO_Old.SetAnaSource(0,SourceSize);
    KittyNLO_Old.KillTheCat();

    CATS KittyNLO;
    KittyNLO.SetMaxNumThreads(1);
    KittyNLO.SetMomBins(85,0,340);
    AnalysisObject.SetUpCats_pL(KittyNLO,"NLO_Coupled_SPD","Gauss");
    KittyNLO.SetAnaSource(0,SourceSize);
    for(unsigned short usCh=1; usCh<=6; usCh++){
        //remove the p-waves
        KittyNLO.RemoveExternalWaveFunction(usCh,1);
    }
    for(unsigned short usCh=1; usCh<=3; usCh++){
        //remove the d-waves
        KittyNLO.RemoveExternalWaveFunction(usCh,2);
    }
    for(unsigned short usCh=7; usCh<=15; usCh++){
        //remove all coupling
        //KittyNLO.RemoveExternalWaveFunction(usCh,0);

        //remove all coupling, but the s-wave SN->LN (3S1)
        //if(usCh!=8) KittyNLO.RemoveExternalWaveFunction(usCh,0);

        //remove all coupling, but the s-wave SN->LN (1S0 and 3S1)
        if(usCh!=7&&usCh!=8) KittyNLO.RemoveExternalWaveFunction(usCh,0);

        //remove all coupling, but the s,d-wave SN->LN
        //if(usCh!=8&&usCh!=13) KittyNLO.RemoveExternalWaveFunction(usCh,0);

        //if(usCh!=8&&usCh!=10&&usCh!=13) KittyNLO.RemoveExternalWaveFunction(usCh,0);

        //remove all coupling to p-waves
        //if(usCh==11||usCh==12) KittyNLO.RemoveExternalWaveFunction(usCh,0);
    }


//KittyNLO.SetChannelWeight(8,KittyNLO.GetChannelWeight(8)*2./3.);
//KittyNLO.SetChannelWeight(13,KittyNLO.GetChannelWeight(13)*2./3.);
//KittyNLO.SetChannelWeight(8,KittyNLO.GetChannelWeight(8)*2./3.);
//KittyNLO.SetChannelWeight(10,3./4.);
//KittyNLO.SetChannelWeight(13,3./4.);
    KittyNLO.KillTheCat();

    unsigned short NumChannels = KittyNLO.GetNumChannels();
    double* OriginalWeight = new double [NumChannels];

    TFile fOutput("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/Test_New_pLambda/fOutput.root","recreate");

    TGraph gOld;
    gOld.SetName("gOld");
    gOld.Set(KittyNLO_Old.GetNumMomBins());
    for(unsigned uBin=0; uBin<KittyNLO_Old.GetNumMomBins(); uBin++){
        gOld.SetPoint(uBin,KittyNLO_Old.GetMomentum(uBin),LambdaPar*KittyNLO_Old.GetCorrFun(uBin)+1.-LambdaPar);
    }

    TGraph gFull;
    gFull.SetName("gFull");
    gFull.Set(KittyNLO.GetNumMomBins());
    for(unsigned uBin=0; uBin<KittyNLO.GetNumMomBins(); uBin++){
        gFull.SetPoint(uBin,KittyNLO.GetMomentum(uBin),LambdaPar*KittyNLO.GetCorrFun(uBin)+1.-LambdaPar);
    }

    const unsigned NumMomForRad = 6;
    double MomentaForRad[NumMomForRad];
    MomentaForRad[0] = 50;
    MomentaForRad[1] = 100;
    MomentaForRad[2] = 150;
    MomentaForRad[3] = 200;
    MomentaForRad[4] = 250;
    MomentaForRad[5] = 300;
    const unsigned NumRadPts = 256;
    const double RadMin = 0;
    const double RadMax = 10;
    const double RadBinWidth = (RadMax-RadMin)/double(NumRadPts);

    TGraph* gWF_3S1 = new TGraph[NumMomForRad];
    TGraph* gAS_3S1 = new TGraph[NumMomForRad];

    TGraph* gWF_3P1 = new TGraph[NumMomForRad];
    TGraph* gAS_3P1 = new TGraph[NumMomForRad];

    TGraph* gWF_3D1 = new TGraph[NumMomForRad];
    TGraph* gAS_3D1 = new TGraph[NumMomForRad];
    for(unsigned uMomRad=0; uMomRad<NumMomForRad; uMomRad++){
        gWF_3S1[uMomRad].SetName(TString::Format("gWF_3S1_%.0f",MomentaForRad[uMomRad]));
        gWF_3S1[uMomRad].Set(RadBinWidth);

        gAS_3S1[uMomRad].SetName(TString::Format("gAS_3S1_%.0f",MomentaForRad[uMomRad]));
        gAS_3S1[uMomRad].Set(RadBinWidth);

        gWF_3P1[uMomRad].SetName(TString::Format("gWF_3P1_%.0f",MomentaForRad[uMomRad]));
        gWF_3P1[uMomRad].Set(RadBinWidth);

        gAS_3P1[uMomRad].SetName(TString::Format("gAS_3P1_%.0f",MomentaForRad[uMomRad]));
        gAS_3P1[uMomRad].Set(RadBinWidth);

        gWF_3D1[uMomRad].SetName(TString::Format("gWF_3D1_%.0f",MomentaForRad[uMomRad]));
        gWF_3D1[uMomRad].Set(RadBinWidth);

        gAS_3D1[uMomRad].SetName(TString::Format("gAS_3D1_%.0f",MomentaForRad[uMomRad]));
        gAS_3D1[uMomRad].Set(RadBinWidth);

        for(unsigned uRad=0; uRad<NumRadPts; uRad++){
            double RADIUS = RadBinWidth*0.5+RadBinWidth*double(uRad);

            gWF_3S1[uMomRad].SetPoint(uRad,RADIUS,
                                      std::abs(KittyNLO.EvalRadialWaveFunction(KittyNLO.GetMomBin(MomentaForRad[uMomRad]),1,0,RADIUS,true)));
            gAS_3S1[uMomRad].SetPoint(uRad,RADIUS,
                                      std::abs(KittyNLO.EvalAsymptoticRadialWF(KittyNLO.GetMomBin(MomentaForRad[uMomRad]),1,0,RADIUS,true)));

            gWF_3P1[uMomRad].SetPoint(uRad,RADIUS,
                                      std::abs(KittyNLO.EvalRadialWaveFunction(KittyNLO.GetMomBin(MomentaForRad[uMomRad]),1,1,RADIUS,true)));
            gAS_3P1[uMomRad].SetPoint(uRad,RADIUS,
                                      std::abs(KittyNLO.EvalAsymptoticRadialWF(KittyNLO.GetMomBin(MomentaForRad[uMomRad]),1,1,RADIUS,true)));

            gWF_3D1[uMomRad].SetPoint(uRad,RADIUS,
                                      std::abs(KittyNLO.EvalRadialWaveFunction(KittyNLO.GetMomBin(MomentaForRad[uMomRad]),1,2,RADIUS,true)));
            gAS_3D1[uMomRad].SetPoint(uRad,RADIUS,
                                      std::abs(KittyNLO.EvalAsymptoticRadialWF(KittyNLO.GetMomBin(MomentaForRad[uMomRad]),1,2,RADIUS,true)));
        }
    }


    TGraph* gChannel = new TGraph[NumChannels];
    TGraph* gWeightedChannel = new TGraph[NumChannels];

    for(unsigned short usCh=0; usCh<NumChannels; usCh++){
         OriginalWeight[usCh] = KittyNLO.GetChannelWeight(usCh);
    }



//haidenbauer has:
//channel 0 s-wave only 1/4
//channel 1 s-wave only 3/4
//channel 8 3/4 ?
//channel 10 3/4 ?
/*
for(unsigned short usCh=0; usCh<NumChannels; usCh++){
    KittyNLO.SetChannelWeight(usCh,0);
}
KittyNLO.SetChannelWeight(0,1./4.);
KittyNLO.RemoveExternalWaveFunction(0,1);

KittyNLO.SetChannelWeight(1,3./4.);
KittyNLO.RemoveExternalWaveFunction(1,1);
KittyNLO.RemoveExternalWaveFunction(1,2);

KittyNLO.SetChannelWeight(8,3./4.);
KittyNLO.SetChannelWeight(10,3./4.);

KittyNLO.KillTheCat();
*/
    for(unsigned short usCh=0; usCh<NumChannels; usCh++){
        gChannel[usCh].SetName(TString::Format("gChannel_%u",usCh));
        gChannel[usCh].Set(KittyNLO.GetNumMomBins());

        gWeightedChannel[usCh].SetName(TString::Format("gWeightedChannel_%u",usCh));
        gWeightedChannel[usCh].Set(KittyNLO.GetNumMomBins());

        for(unsigned short usCh2=0; usCh2<NumChannels; usCh2++){
            KittyNLO.SetChannelWeight(usCh2,0);
        }

        double BASELINE = (usCh>=7);

        KittyNLO.SetChannelWeight(usCh,1);
        KittyNLO.KillTheCat();
        for(unsigned uBin=0; uBin<KittyNLO.GetNumMomBins(); uBin++){
            gChannel[usCh].SetPoint(uBin,KittyNLO.GetMomentum(uBin),KittyNLO.GetCorrFun(uBin)+BASELINE);
        }

        KittyNLO.SetChannelWeight(usCh,OriginalWeight[usCh]);
        KittyNLO.KillTheCat();
        for(unsigned uBin=0; uBin<KittyNLO.GetNumMomBins(); uBin++){
            gWeightedChannel[usCh].SetPoint(uBin,KittyNLO.GetMomentum(uBin),KittyNLO.GetCorrFun(uBin)+BASELINE);
        }

        gChannel[usCh].Write();
        gWeightedChannel[usCh].Write();
    }

    gFull.Write();
    gOld.Write();

    for(unsigned uMomRad=0; uMomRad<NumMomForRad; uMomRad++){
        gWF_3S1[uMomRad].Write();
        gAS_3S1[uMomRad].Write();
    }
    for(unsigned uMomRad=0; uMomRad<NumMomForRad; uMomRad++){
        gWF_3P1[uMomRad].Write();
        gAS_3P1[uMomRad].Write();
    }
    for(unsigned uMomRad=0; uMomRad<NumMomForRad; uMomRad++){
        gWF_3D1[uMomRad].Write();
        gAS_3D1[uMomRad].Write();
    }

    delete [] gWF_3S1;
    delete [] gAS_3S1;
    delete [] gWF_3P1;
    delete [] gAS_3P1;
    delete [] gWF_3D1;
    delete [] gAS_3D1;
    delete [] OriginalWeight;
    delete [] gChannel;
    delete [] gWeightedChannel;
}
/*
void pn_B2(){
    CATSparameters cPars(CATSparameters::tSource,1,true);
    cPars.SetParameter(0,1.2);
    CATS Kitty;
    Kitty.SetMomBins(30,0,300);
    Kitty.SetThetaDependentSource(false);
    Kitty.SetAnaSource(GaussSource, cPars);
    Kitty.SetUseAnalyticSource(true);
    Kitty.SetMomentumDependentSource(false);
    //Kitty.SetThetaDependentSource(false);
    Kitty.SetExcludeFailedBins(false);

    Kitty.SetQ1Q2(0);
    Kitty.SetPdgId(2212, 2213);
    Kitty.SetRedMass( 0.5*Mass_p );

    Kitty.SetNumChannels(1);
    Kitty.SetNumPW(0,1);
    Kitty.SetSpin(0,1);
    Kitty.SetChannelWeight(0, 1.);

    CATSparameters pPars(CATSparameters::tPotential,8,true);
    double PotPars[8]={NN_AV18,v18_Coupled3P2,0,-1,1,1,0,1};
    pPars.SetParameters(PotPars);
    Kitty.SetShortRangePotential(0,0,fDlmPot,pPars);
    Kitty.KillTheCat();

    TFile fOutput("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pn_B2/fOutput.root","recreate");

    TGraph gCk;
    gCk.SetName("gCk");
    TGraph gB2_Gauss;
    gB2_Gauss.SetName("gB2_Gauss");
    TGraph gB2_Hulthen;
    gB2_Hulthen.SetName("gB2_Hulthen");

    for(unsigned uBin=0; uBin<Kitty.GetNumMomBins(); uBin++){
        gCk.SetPoint(uBin,Kitty.GetMomentum(uBin),Kitty.GetCorrFun(uBin));
        gB2_Gauss.SetPoint(uBin,Kitty.GetMomentum(uBin),Kitty.GetB2_Gauss(uBin));
        gB2_Hulthen.SetPoint(uBin,Kitty.GetMomentum(uBin),Kitty.GetB2_Hulthen(uBin));
    }

    gCk.Write();
    gB2_Gauss.Write();
    gB2_Hulthen.Write();
}
*/


CATS* ReturnCk_CAT;
double ReturnCk(double* x){
    return (0.91*ReturnCk_CAT->EvalCorrFun(*x)-1.00);
}

//goal: see if Int C(k)-1 = 0
void Integrate_Ck(){

    const double kMin = 0;
    const double kMax = 1500;
    const unsigned NumMomBins = 500;

    DLM_CommonAnaFunctions AO_pL;
    CATS Kitty_pL;
    Kitty_pL.SetMaxNumThreads(4);
    Kitty_pL.SetMomBins(NumMomBins,kMin,kMax);
    AO_pL.SetUpCats_pL(Kitty_pL,"Usmani","Gauss",0,202);
    Kitty_pL.SetAnaSource(0,1.2);
    Kitty_pL.SetAnaSource(1,2);
    Kitty_pL.KillTheCat();

    ReturnCk_CAT = &Kitty_pL;
    double Pars[1];
    DLM_INT_SetFunction(ReturnCk,Pars,0);
    double NORM = DLM_INT_aSimpsonWiki(kMin,kMax);
    printf("NORM = %f\n",NORM);

}

void Silly(){

    int a = 11;
    int b = 13;
    int* pointer;
    printf("a = (%p) %i\n",&a,a);
    printf("b = (%p) %i\n",&b,b);
    pointer = &a;
    printf("pointer = (%p) %i\n",pointer,*pointer);
    *pointer *= 2;
    printf("    a = (%p) %i\n",&a,a);
    printf("    b = (%p) %i\n",&b,b);
    pointer = &b;
    printf("pointer = (%p) %i\n",pointer,*pointer);
    *pointer *= 3;
    printf("    a = (%p) %i\n",&a,a);
    printf("    b = (%p) %i\n",&b,b);

}

void Test_CatsParticle(){

    const unsigned NumIter = 1000000;

    double momx,momy,momz;
    const double MomSpread = 0.8;
    const double Mass = 0.14;
    const double MassProd1 = 0.938;
    const double MassProd2 = 0.938;


    TRandom3 rangen(11);

    TH2F* hTLV_k_E = new TH2F("hTLV_k_E","hTLV_k_E",1024,0,4,128,0,4);
    TH2F* hCats_k_E = new TH2F("hCats_k_E","hCats_k_E",1024,0,4,128,0,4);

    TH2F* hTLV_ks_Ecm = new TH2F("hTLV_ks_Ecm","hTLV_ks_Ecm",1024,0,4,128,0,4);
    TH2F* hCats_ks_Ecm = new TH2F("hCats_ks_Ecm","hCats_ks_Ecm",1024,0,4,128,0,4);

    TH2F* hTLV_ksP_Ecm = new TH2F("hTLV_ksP_Ecm","hTLV_ksP_Ecm",1024,0,4,128,0,4);


    for(unsigned uIter=0; uIter<NumIter; uIter++){
        momx = rangen.Gaus(0,MomSpread);
        momy = rangen.Gaus(0,MomSpread);
        momz = rangen.Gaus(0,MomSpread);
        TLorentzVector TLV_pion1;
        TLV_pion1.SetXYZM(momx,momy,momz,Mass);
        CatsParticle Cats_pion1;
        Cats_pion1.Set(0,0,0,0,sqrt(momx*momx+momy*momy+momz*momz+Mass*Mass),momx,momy,momz);

        momx = rangen.Gaus(0,MomSpread);
        momy = rangen.Gaus(0,MomSpread);
        momz = rangen.Gaus(0,MomSpread);
        TLorentzVector TLV_pion2;
        TLV_pion2.SetXYZM(momx,momy,momz,Mass);
        CatsParticle Cats_pion2;
        Cats_pion2.Set(0,0,0,0,sqrt(momx*momx+momy*momy+momz*momz+Mass*Mass),momx,momy,momz);

        TLorentzVector TLV_Sum = TLV_pion1+TLV_pion2;
        TLorentzVector TLV_Diff = TLV_pion1-TLV_pion2;

        TVector3 BoostVec = -TLV_Sum.BoostVector();

        TLorentzVector TLV_Sum_Boost = TLV_Sum;
        TLorentzVector TLV_Diff_Boost = TLV_Diff;

        TLV_Sum_Boost.Boost(BoostVec);
        TLV_Diff_Boost.Boost(BoostVec);

        hTLV_k_E->Fill(TLV_Diff.P()*0.5,TLV_Sum.M());
        hTLV_ks_Ecm->Fill(TLV_Diff_Boost.P()*0.5,TLV_Sum_Boost.M());

        CatsParticlePair PionPair;
        PionPair.SetPair(Cats_pion1,Cats_pion2,false,false);

        CatsParticlePair PionPairBoost;
        PionPairBoost.SetPair(Cats_pion1,Cats_pion2,false,true);

        hCats_k_E->Fill(PionPair.GetP()*0.5,PionPair.GetSum().Mag());
        hCats_ks_Ecm->Fill(PionPairBoost.GetP()*0.5,PionPairBoost.GetSum().Mag());

        double E_prod = PionPair.GetSum().Mag();
        double Ks_prod = sqrt(pow(E_prod,4.)+pow(MassProd1*MassProd1-MassProd2*MassProd2,2.)-2.*E_prod*E_prod*(MassProd1*MassProd1+MassProd2*MassProd2))/(2.*E_prod);
        hTLV_ksP_Ecm->Fill(Ks_prod,TLV_Sum_Boost.M());
    }

    TFile fOutput("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Temp/Test_CatsParticle.root","recreate");
    hTLV_k_E->Write();
    hCats_k_E->Write();
    hTLV_ks_Ecm->Write();
    hCats_ks_Ecm->Write();
    hTLV_ksP_Ecm->Write();

    delete hTLV_k_E;
    delete hCats_k_E;
    delete hTLV_ks_Ecm;
    delete hCats_ks_Ecm;
    delete hTLV_ksP_Ecm;
}

void Test_OliSigma(){
    const unsigned NumBins =  16;
    const double kMin = 0;
    const double kMax = 160;
    TH1F* hDummy = new TH1F("hDummy","hDummy",NumBins,kMin,kMax);

    double PotPar = 0;
    double SourceSize = 1.5;
    double Momentum;

    for(unsigned uBin=0; uBin<NumBins; uBin++){
        Momentum = hDummy->GetBinCenter(uBin+1);
        //printf(" ROOT = %f; C++ = %f\n",
               //Lednicky_gauss_Sigma0_ROOT(Momentum,&SourceSize,&PotPar),
               //Lednicky_gauss_Sigma0(Momentum,&SourceSize,&PotPar));
    }

    delete hDummy;
}

void Test_pXi_Yukawa(){

    const unsigned NumBins = 200;
    const double kMin = 0;
    const double kMax = 400;

    DLM_CommonAnaFunctions AnaObj;

    CATS KittyOld;
    KittyOld.SetMomBins(NumBins,kMin,kMax);
    AnaObj.SetUpCats_pXim(KittyOld,"pXim_HALQCD1","Gauss",-12,0);
    KittyOld.SetAnaSource(0,1.0);
    KittyOld.KillTheCat();

    CATS KittyNew;
    KittyNew.SetMomBins(NumBins,kMin,kMax);
    AnaObj.SetUpCats_pXim(KittyNew,"pXim_HALQCD1","Gauss",12,0);
    KittyNew.SetAnaSource(0,1.0);
    KittyNew.KillTheCat();

    TFile fOutput("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Temp/Test_pXi_Yukawa.root","recreate");

    TGraph gKittyOld;
    gKittyOld.SetName(TString::Format("gKittyOld"));
    gKittyOld.Set(NumBins);
    for(unsigned uBin=0; uBin<NumBins; uBin++){
        gKittyOld.SetPoint(uBin,KittyOld.GetMomentum(uBin),KittyOld.GetCorrFun(uBin));
    }
    gKittyOld.Write();

    TGraph gKittyNew;
    gKittyNew.SetName(TString::Format("gKittyNew"));
    gKittyNew.Set(NumBins);
    for(unsigned uBin=0; uBin<NumBins; uBin++){
        gKittyNew.SetPoint(uBin,KittyNew.GetMomentum(uBin),KittyNew.GetCorrFun(uBin));
    }
    gKittyNew.Write();

}

void QA_pd(){
    const double SourceSize = 1.2;

    const double kMin=0;
    const double kMax=130;
    const unsigned NumMomBins = 26;

    CATSparameters* cPars;
    cPars = new CATSparameters(CATSparameters::tSource, 1, true);
    cPars->SetParameter(0, SourceSize);
    CATS Kitty;
    Kitty.SetAnaSource(GaussSource, *cPars);
    Kitty.SetMomBins(NumMomBins,kMin,kMax);
    Kitty.SetAnaSource(GaussSource, *cPars);
    Kitty.SetAnaSource(0,SourceSize);
    Kitty.SetUseAnalyticSource(true);
    Kitty.SetMomentumDependentSource(false);
    Kitty.SetThetaDependentSource(false);
    Kitty.SetExcludeFailedBins(false);

printf("Fuck\n");
    DLM_Histo<complex<double>>*** ExternalWF=NULL;
    ExternalWF = Init_pd_Sebastian("/home/dmihaylov/CernBox/CatsFiles/Interaction/pd/Sebastian/",Kitty,1,400);
    for(unsigned uCh=0; uCh<Kitty.GetNumChannels(); uCh++){
        Kitty.SetExternalWaveFunction(uCh,0,ExternalWF[0][uCh][0],ExternalWF[1][uCh][0]);
    }
printf("WF ready\n");
    Kitty.SetChannelWeight(0,0);
    Kitty.SetChannelWeight(1,1);
    Kitty.KillTheCat();

    TFile* OutputFile = new TFile(
                TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/QA_pd/QA_pd_Q.root"),"recreate");
printf("File Created\n");
    TGraph gKitty;
    gKitty.SetName(TString::Format("gKitty"));
    gKitty.Set(NumMomBins);
    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        printf("C(%.2f) = %.2f\n",Kitty.GetMomentum(uBin),Kitty.GetCorrFun(uBin));
        gKitty.SetPoint(uBin,Kitty.GetMomentum(uBin),Kitty.GetCorrFun(uBin));
    }
    gKitty.Write();

    TGraph* gWF_Re_D = new TGraph [NumMomBins];
    TGraph* gWF_Re_Q = new TGraph [NumMomBins];
    TGraph* gWF_Re_R = new TGraph [NumMomBins];
    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        unsigned COUNTER=0;
        gWF_Re_D[uBin].SetName(TString::Format("gWF_Re_D_%.0f",Kitty.GetMomentum(uBin)));
        gWF_Re_Q[uBin].SetName(TString::Format("gWF_Re_Q_%.0f",Kitty.GetMomentum(uBin)));
        gWF_Re_R[uBin].SetName(TString::Format("gWF_Re_R_%.0f",Kitty.GetMomentum(uBin)));
        for(double RAD=0.05; RAD<100; RAD+=0.05){
            gWF_Re_D[uBin].SetPoint(COUNTER,RAD,abs(Kitty.EvalRadialWaveFunction(uBin,0,0,RAD,true)));
            gWF_Re_Q[uBin].SetPoint(COUNTER,RAD,abs(Kitty.EvalRadialWaveFunction(uBin,1,0,RAD,true)));
            gWF_Re_R[uBin].SetPoint(COUNTER,RAD,abs(Kitty.EvalReferenceRadialWF(uBin,0,RAD,true)));

            COUNTER++;
        }
        gWF_Re_D[uBin].Write();
        gWF_Re_Q[uBin].Write();
        gWF_Re_R[uBin].Write();
    }


    CleanUpWfHisto(Kitty,ExternalWF);

    delete [] gWF_Re_D;
    delete [] gWF_Re_Q;
    delete [] gWF_Re_R;
    delete OutputFile;
    delete cPars;
}

void pL_lampar(){

    for(unsigned uVarSL=0; uVarSL<5; uVarSL++){
        for(unsigned uVarXI=0; uVarXI<2; uVarXI++){
            double lampar[5];
            DLM_CommonAnaFunctions AnalysisObject;
            AnalysisObject.SetCatsFilesFolder("/home/dmihaylov/CernBox/CatsFiles");
            AnalysisObject.SetUpLambdaPars_pL("pp13TeV_HM_DimiJun20",0,uVarSL+10*(3*uVarXI),lampar);
            printf("%u %u:\n",uVarSL,uVarXI);
            printf(" genuine: %.4f\n",lampar[0]*100.);
            printf(" pSigma0: %.4f\n",lampar[1]*100.);
            printf("   (S:L): %.4f\n",lampar[1]/lampar[0]);
            printf("    pXim: %.4f\n",lampar[2]*100.);
            printf("    pXi0: %.4f\n",lampar[2]*100.);
            printf("  (XI:LS): %.4f\n",2.*lampar[2]/(lampar[0]+lampar[1]));
            printf("    flat: %.4f\n",lampar[3]*100.-lampar[2]*100.);
            printf("   misid: %.4f\n",lampar[4]*100.);
            printf("     sum: %.4f\n",(lampar[0]+lampar[1]+lampar[2]+lampar[3]+lampar[4])*100.);
        }
    }

}

void TestLedniShit(){
    const unsigned NumMomBins = 10;
    const double kMin = 0;
    const double kMax = 240;
    //object for the parameters to be used by the source function
    CATSparameters SOURCE_PARS(CATSparameters::tSource, 2, true);
    //set the source size to 1.2 fm
    SOURCE_PARS.SetParameter(0, 1.2);
    //here we define the DLM_Ck object for the Lednicky interaction. The arguments are:
    //#source parameters, #interaction parameters, k min, k max, pointer to the function to be used
    DLM_Ck Ck_Lednicky(1, 4, NumMomBins, kMin, kMax, Lednicky_SingletTriplet);
    //set the source and interaction (potential) parameters.
    //for the interaction we have used the scattering lengths (pars 0 and 2) and the effective ranges (1 and 3)
    //for the Usmani potential.
    Ck_Lednicky.SetSourcePar(0, SOURCE_PARS.GetParameter(0));
    Ck_Lednicky.SetPotPar(0, 2.88);
    Ck_Lednicky.SetPotPar(1, 2.92);
    Ck_Lednicky.SetPotPar(2, 1.66);
    Ck_Lednicky.SetPotPar(3, 3.78);
    Ck_Lednicky.Update();


    //object for the parameters to be used by the source function
    CATSparameters SOURCE_PARS_2(CATSparameters::tSource, 2, true);
    //set the source size to 1.2 fm
    SOURCE_PARS_2.SetParameter(0, 1.2);
    //here we define the DLM_Ck object for the Lednicky interaction. The arguments are:
    //#source parameters, #interaction parameters, k min, k max, pointer to the function to be used
    DLM_Ck Ck_Lednicky_2(1, 6, NumMomBins, kMin, kMax, Lednicky_2channel);
    //set the source and interaction (potential) parameters.
    //for the interaction we have used the scattering lengths (pars 0 and 2) and the effective ranges (1 and 3)
    //for the Usmani potential.
    Ck_Lednicky_2.SetSourcePar(0, SOURCE_PARS_2.GetParameter(0));
    Ck_Lednicky_2.SetPotPar(0, 2.88);
    Ck_Lednicky_2.SetPotPar(1, 2.92);
    Ck_Lednicky_2.SetPotPar(2, 1.66);
    Ck_Lednicky_2.SetPotPar(3, 3.78);
    Ck_Lednicky_2.SetPotPar(4, 0.25);
    Ck_Lednicky_2.SetPotPar(5, 0.75);
    Ck_Lednicky_2.Update();
/*
    for(){

    }
*/
}

void CkHisto_test(){
    //10 bins, 0-100 MeV
    DLM_Ck* CkTest = new DLM_Ck(10,0,100);

    for(unsigned uBin=0; uBin<10; uBin++){
        CkTest->SetBinContent(uBin,uBin*uBin);
    }

    printf(" CkTest->Eval(5)=%.2f\n",CkTest->Eval(5));
    printf(" CkTest->Eval(50)=%.2f\n",CkTest->Eval(50));
    printf(" CkTest->Eval(95)=%.2f\n",CkTest->Eval(95));

    delete CkTest;
}

int main(int argc, char *argv[])
{

    //TFile* f1 = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Temp/temp.root","recreate");
    //TNtuple* nt1 = new TNtuple("nt1","nt1",
    //                            "UITERfirst:UITERlast:f0Inv:d0:Chi2fit:NDFfit:Chi2refit:dChi2:dChi2Ndf");

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

    //CkHisto_test();

    //TestLedniShit();
    //cout << file_status("/home/dmihaylov/Dudek_Ubuntu/pNb_cosThCM_0_1.root") << endl;
    //pL_lampar();
    //QA_pd();
    //pipi_test();
    //TestRandom();
    //pp_in_txtfile();
    //pp_pLambda_pXi_Ratios();
    //cout<<tLab_pCm(200,938,938)<<endl;
    //Test_CatsParticle();
    //Test_OliSigma();


    //! FOR THE CATS TUTORIAL 2019
    //GenerateFakeProtonLambda();

//double* lampar = new double[1000];
//DLM_CommonAnaFunctions FUCK;
//FUCK.SetUpLambdaPars_pL("pp13TeV_HM_March19",0,0,lampar);
//printf("primary %f\n",lampar[0]);
//printf("pSigma0->pL %f\n",lampar[1]);
//printf("pXim->pL %f\n",lampar[2]);
//printf("flat %f\n",lampar[3]);
//printf("misid %f\n",lampar[4]);

//CATS CAT;
//Init_pantip_Haidenbauer("/home/dmihaylov/CernBox/CATS_potentials/Haidenbauer/p_antip/WithCoulomb/", CAT, 1);
//return 0;

    //TestHaide_pL_pWaves();
    //Test_pXi_Yukawa();

    //TestTomPotential();

    //ExecuteCFmT(argc,ARGV);
    //TestSetSameParInFit(); return 0;


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

    //PION_ANA(argc,ARGV);
    //THESIS_PLOTS(argc,ARGV);
    //PLAMBDA_1_MAIN(argc,ARGV);
    //CALL_BERNIE_AND_VALE();
    //FEMTOBOYZ_MAIN(argc,ARGV);
    //GENBOD(argc,ARGV);
    //GerhardMAIN(argc,ARGV);
    //Main_pSigma();
    //MIXEDEVENTS(argc,ARGV);
    //SOURCESTUDIES(argc,ARGV);
    //KAONPROTON_MAIN(argc,ARGV);
    //UNFOLD_MAIN(argc,ARGV);
    //OTHERTASKS(argc,ARGV);
    //DimiMcPlayground_MAIN(argc,ARGV);

    //Silly();
/*

void ScatParsFromRandPotential(const TString OutputFolder,
                               const unsigned RandomSeed, const unsigned& NumIter,
                               const double& V0_min, const double& V0_max,
                               const double& mu0_min, const double& mu0_max,
                               const double& V1_min, const double& V1_max,
                               const double& mu1_min, const double& mu1_max,
                               const double& V2_min, const double& V2_max,
                               const double& mu2_min, const double& mu2_max)
*/
//!Looking at Fig,. 1 I realized that also NCS97f yields a nice description of the cups in p-Lambda. The nice thing is that Haidenbauer has the potentials at hand, and thus super easy can compute the relative wave functions in the same format as for chiEFT

    //TestThetaSource();
    //PlotResonanceAngles(1.0,0.5);
    //PlotResonanceAngles(1.0,1.0);
    //PlotResonanceAngles(1.0,2.0);
    //PlotResonanceAngles(1.0,4.0);
    //PlotResonanceAngles(1.0,8.0);

    //MickeyMouseEffectOnSourceBasedOnTheta();

    //pn_B2();

    //TestCommonInit();
    //plot_pp();
    //plot_pp_2();
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

    //TestCauchy();
    //Test_New_pLambda();
    //Integrate_Ck();


    for(int iARG=1; iARG<argc; iARG++){
        delete [] ARGV[iARG];
    }
    if(ARGV) delete [] ARGV;

    long long ExeTime = TIMER.Stop()/1000.;
//for(unsigned uSec=0; uSec<10; uSec++){
//    printf("%f\n",double(ExeTime)/1000.);
//    usleep(1000e3);
//    ExeTime = TIMER.Stop()/1000.;
//}
    char* strtime = new char [128];
    ShowTime(ExeTime,strtime,0,true,6);
    printf("The script terminated after: %s\n",strtime);

    delete [] strtime;

    //delete nt1;
    //delete f1;
    return 0;
}
