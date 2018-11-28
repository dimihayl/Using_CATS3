
#include "LamLamMcExclusionMAIN.h"

#include <iostream>
#include <stdio.h>
#include <string.h>
//#include <omp.h>

#include "CATS.h"
#include "CATStools.h"
#include "DLM_Potentials.h"
#include "DLM_Source.h"
#include "DLM_CppTools.h"

#include "TGraph.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TNtuple.h"

using namespace std;


void plot_pp_Ck(){
    const double kMin_pp=0;
    const double kMax_pp = 400;
    const unsigned NumMomBins_pp = 100;

    //starting value, do not worry about it too much
    const double GaussSourceSize = 1.15;
    const double Mass_p=938.272;

//ADVANCED***
    //#,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
    double PotPars1S0[10]={0,0,NN_AV18,v18_Coupled3P2,1,1,1,0,0,0};
    double PotPars3P0[10]={0,0,NN_AV18,v18_Coupled3P2,1,1,1,1,1,0};
    double PotPars3P1[10]={0,0,NN_AV18,v18_Coupled3P2,1,1,1,1,1,1};
    double PotPars3P2[10]={0,0,NN_AV18,v18_Coupled3P2,1,1,1,1,1,2};

    CATSparameters cPotPars1S0(CATSparameters::tPotential,8,true);
    cPotPars1S0.SetParameters(&PotPars1S0[2]);
    CATSparameters cPotPars3P0(CATSparameters::tPotential,8,true);
    cPotPars3P0.SetParameters(&PotPars3P0[2]);
    CATSparameters cPotPars3P1(CATSparameters::tPotential,8,true);
    cPotPars3P1.SetParameters(&PotPars3P1[2]);
    CATSparameters cPotPars3P2(CATSparameters::tPotential,8,true);
    cPotPars3P2.SetParameters(&PotPars3P2[2]);

    const double Weight1S0 = 3./12.;
    const double Weight3P0 = 1./12.;
    const double Weight3P1 = 3./12.;
    const double Weight3P2 = 5./12.;

    CATS AB_pp;
    double Pars_pp[4] = {0,0,0,GaussSourceSize};
    CATSparameters cPars_pp(CATSparameters::tSource,1,true);
    cPars_pp.SetParameters(&Pars_pp[3]);

    AB_pp.SetAnaSource(GaussSource, cPars_pp);
    AB_pp.SetUseAnalyticSource(true);

    AB_pp.SetExcludeFailedBins(false);
    AB_pp.SetMomBins(NumMomBins_pp,kMin_pp,kMax_pp);

    AB_pp.SetNumChannels(4);
    AB_pp.SetNumPW(0,2);
    AB_pp.SetNumPW(1,2);
    AB_pp.SetNumPW(2,2);
    AB_pp.SetNumPW(3,2);
    AB_pp.SetSpin(0,0);
    AB_pp.SetSpin(1,1);
    AB_pp.SetSpin(2,1);
    AB_pp.SetSpin(3,1);
    AB_pp.SetChannelWeight(0, Weight1S0);
    AB_pp.SetChannelWeight(1, Weight3P0);
    AB_pp.SetChannelWeight(2, Weight3P1);
    AB_pp.SetChannelWeight(3, Weight3P2);

    AB_pp.SetQ1Q2(1);
    AB_pp.SetPdgId(2212, 2212);
    AB_pp.SetRedMass( 0.5*Mass_p );

    AB_pp.SetShortRangePotential(0,0,fDlmPot,cPotPars1S0);
    AB_pp.SetShortRangePotential(1,1,fDlmPot,cPotPars3P0);
    AB_pp.SetShortRangePotential(2,1,fDlmPot,cPotPars3P1);
    AB_pp.SetShortRangePotential(3,1,fDlmPot,cPotPars3P2);
//AB_pp.SetGamow(true);
    AB_pp.KillTheCat();


    CATS AB_pp_2;
    AB_pp_2.SetAnaSource(GaussSource, cPars_pp);
    AB_pp_2.SetUseAnalyticSource(true);

    AB_pp_2.SetExcludeFailedBins(false);
    AB_pp_2.SetMomBins(NumMomBins_pp,kMin_pp,kMax_pp);

    AB_pp_2.SetNumChannels(2);
    AB_pp_2.SetNumPW(0,1);
    AB_pp_2.SetNumPW(1,1);
    AB_pp_2.SetSpin(0,0);
    AB_pp_2.SetSpin(1,1);
    AB_pp_2.SetChannelWeight(0, 0.25);
    AB_pp_2.SetChannelWeight(1, 0.75);

    AB_pp_2.SetQ1Q2(1);
    AB_pp_2.SetPdgId(2212, 2212);
    AB_pp_2.SetRedMass( 0.5*Mass_p );

    AB_pp_2.SetShortRangePotential(0,0,fDlmPot,cPotPars1S0);

//AB_pp_2.SetGamow(true);
    AB_pp_2.KillTheCat();

    TGraph gSPwaves;
    gSPwaves.Set(NumMomBins_pp);
    gSPwaves.SetName("gSPwaves");
    TGraph gSwave;
    gSwave.Set(NumMomBins_pp);
    gSwave.SetName("gSwave");

    for(unsigned uBin=0; uBin<NumMomBins_pp; uBin++){
        printf("C(%.2f)=%.3f vs %.3f\n",AB_pp.GetMomentum(uBin),AB_pp.GetCorrFun(uBin),AB_pp_2.GetCorrFun(uBin));
        gSPwaves.SetPoint(uBin,AB_pp.GetMomentum(uBin),AB_pp.GetCorrFun(uBin));
        gSwave.SetPoint(uBin,AB_pp_2.GetMomentum(uBin),AB_pp_2.GetCorrFun(uBin));
    }

    TFile* OutputFile = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Temp/Output.root","recreate");

    gSPwaves.Write();
    gSwave.Write();

    delete OutputFile;
}

void Test_pL_1(){
    const double kMin_pL=0;
    const double kMax_pL = 500;
    const unsigned NumMomBins_pp = 500;

    //starting value, do not worry about it too much
    const double GaussSourceSize = 1.15;
    const double Mass_p=938.272;
    const double Mass_L=1115.683;

    double pLamPotPars1S0[10]={0,0,pL_UsmaniOli,0,0,0,0,0,0,0};
    double pLamPotPars3S1[10]={0,0,pL_UsmaniOli,0,0,0,0,1,0,1};

    CATSparameters cPars1S0(CATSparameters::tPotential,8,true);
    cPars1S0.SetParameters(&pLamPotPars1S0[2]);
    CATSparameters cPars3S1(CATSparameters::tPotential,8,true);
    cPars3S1.SetParameters(&pLamPotPars3S1[2]);

    CATS AB_pL;
    double Pars_pL[4] = {0,0,0,GaussSourceSize};
    CATSparameters cPars_pL(CATSparameters::tSource,1,true);
    cPars_pL.SetParameters(&Pars_pL[3]);

    AB_pL.SetAnaSource(GaussSource, cPars_pL);
    AB_pL.SetUseAnalyticSource(true);

    AB_pL.SetExcludeFailedBins(false);
    AB_pL.SetMomBins(NumMomBins_pp,kMin_pL,kMax_pL);

    AB_pL.SetNumChannels(2);
    AB_pL.SetNumPW(0,1);
    AB_pL.SetNumPW(1,1);
    AB_pL.SetSpin(0,0);
    AB_pL.SetSpin(1,1);
    AB_pL.SetChannelWeight(0, 0.25);
    AB_pL.SetChannelWeight(1, 0.75);

    AB_pL.SetQ1Q2(0);
    AB_pL.SetPdgId(2212, 3122);
    AB_pL.SetRedMass( (Mass_p*Mass_L)/(Mass_p+Mass_L) );

    AB_pL.SetShortRangePotential(0,0,fDlmPot,cPars1S0);
    AB_pL.SetShortRangePotential(1,0,fDlmPot,cPars3S1);

    AB_pL.KillTheCat();

    for(unsigned uBin=0; uBin<NumMomBins_pp; uBin++){
        printf("C(%.2f)=%.3f\n",AB_pL.GetMomentum(uBin),AB_pL.GetCorrFun(uBin));
    }
}

int main(int argc, char *argv[])
{

    //TH1F h1;

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

    LamLamMcExclusionMAIN(argc,ARGV);

    //plot_pp_Ck();
    //Test_pL_1();

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
    return 0;
}
