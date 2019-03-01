
#include "CommonAnaFunctions.h"

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

#include "TGraph.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TNtuple.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TNtuple.h"


using namespace std;


void TestCATS3_NewExtWf(TString POT_TYPE){
    //const double kMin=0;
    //const double kMax = 320;
    unsigned NumMomBins;
    double* MomentumBins=NULL;
    double* FitRegion=NULL;
    SetUpBinning_pL("pp13TeV_MB_Run2paper",NumMomBins,MomentumBins,FitRegion);

    printf("NumMomBins=%u\n",NumMomBins);
    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        printf(" Bin %u: %.2f -> %.2f\n", uBin, MomentumBins[uBin], MomentumBins[uBin+1]);
    }

    CATS Kitty;
    Kitty.SetMomBins(NumMomBins,MomentumBins);
    SetUpCats_pL(Kitty,POT_TYPE,"Gauss");
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

void plot_pp(){

    const double SourceSize = 1.0;
    const unsigned NumStabilitySteps = 6;
    const double StabMin = 1;
    const double StabMax = 2;
    const double StabStep = (StabMax-StabMin)/double(NumStabilitySteps-1);

    const double kMin=0;
    const double kFine=200;
    const unsigned NumFineBins = 20;
    const double kFineStep = (kFine-kMin)/double(NumFineBins);
    const double kMax=800;
    const unsigned NumCoarseBins = 20;
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
        double Stability = StabMin+double(uStab)*StabStep;
        printf("\nStability = %.2f\n",Stability);
        CATS Kitty;
        Kitty.SetMomBins(NumMomBins,MomBins);
        SetUpCats_pp(Kitty,"AV18","CleverLevy_Nolan");
        Kitty.SetAnaSource(0,SourceSize);
        Kitty.SetAnaSource(1,Stability);
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

    PLAMBDA_1_MAIN(argc,ARGV);
    //CALL_BERNIE_AND_VALE();

    //TestCommonInit();
    //plot_pp();

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
