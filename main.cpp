
#include "LamLamMcExclusionMAIN.h"
#include "CATS3_Scripts1.h"
#include "ForBernie.h"
#include "ForGerhard.h"

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <omp.h>

#include "CATS.h"
#include "CATStools.h"
#include "DLM_Potentials.h"
#include "DLM_Source.h"
#include "DLM_CppTools.h"
#include "DLM_CkDecomposition.h"
#include "DLM_Random.h"
#include "DLM_Bessel.h"
#include "DLM_Integration.h"

#include "TGraph.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TNtuple.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TNtuple.h"


using namespace std;


void plot_pp_Ck(){
    const double kMin_pp=0;
    const double kMax_pp = 200;
    const unsigned NumMomBins_pp = 50;

    //starting value, do not worry about it too much
    const double GaussSourceSize = 1.2;
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

    AB_pp.SetQ1Q2(0);
    AB_pp.SetPdgId(2212, 2212);
    AB_pp.SetRedMass( 0.5*Mass_p );

    //AB_pp.SetShortRangePotential(0,0,fDlmPot,cPotPars1S0);
    //AB_pp.SetShortRangePotential(1,1,fDlmPot,cPotPars3P0);
    //AB_pp.SetShortRangePotential(2,1,fDlmPot,cPotPars3P1);
    //AB_pp.SetShortRangePotential(3,1,fDlmPot,cPotPars3P2);
//AB_pp.SetGamow(true);
    //AB_pp.SetNotifications(CATS::nSilent);
    //AB_pp.SetMaxNumThreads(1);
    AB_pp.KillTheCat();

/*
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
*/
    DLM_Ck Ck_pp(1,0,AB_pp);
    printf("Eval at 100: %f\n", Ck_pp.Eval(100));
    printf("Eval at 200: %f\n", Ck_pp.Eval(200));
    printf("Eval at 400: %f\n", Ck_pp.Eval(400));
    printf("Eval at 401: %f\n", Ck_pp.Eval(401));

    TGraph gSPwaves;
    gSPwaves.Set(NumMomBins_pp);
    gSPwaves.SetName("gSPwaves");

    TGraph gSwave;
    gSwave.Set(NumMomBins_pp);
    gSwave.SetName("gSwave");

    TGraph gTotWF2_98_0;
    gTotWF2_98_0.Set(NumMomBins_pp);
    gTotWF2_98_0.SetName("gTotWF2_98_0");

    TGraph gTotWF2_98_1;
    gTotWF2_98_1.Set(NumMomBins_pp);
    gTotWF2_98_1.SetName("gTotWF2_98_1");

    for(unsigned uBin=0; uBin<NumMomBins_pp; uBin++){
        //printf("C(%.2f)=%.3f vs %.3f\n",AB_pp.GetMomentum(uBin),AB_pp.GetCorrFun(uBin),AB_pp_2.GetCorrFun(uBin));
        printf("C(%.2f)=%.3f\n",AB_pp.GetMomentum(uBin),AB_pp.GetCorrFun(uBin));
        gSPwaves.SetPoint(uBin,AB_pp.GetMomentum(uBin),AB_pp.GetCorrFun(uBin));
        //gSwave.SetPoint(uBin,AB_pp_2.GetMomentum(uBin),AB_pp_2.GetCorrFun(uBin));
    }

    const unsigned WhichBin = AB_pp.GetMomBin(98);
    const unsigned NumRadBins = 1024;
    const double RadMin=0;
    const double RadMax=16;
    DLM_Histo1D<double> dlmRad(NumRadBins,RadMin,RadMax);
    double Radius;
    for(unsigned uRad=0; uRad<NumRadBins; uRad++){
        Radius = dlmRad.GetBinCenter(uRad);
        //gTotWF2_98_0.SetPoint(uRad,Radius,AB_pp.EvalWaveFun2(WhichBin,Radius,0)*pow(Radius*AB_pp.FmNu()*98.,2));
        //gTotWF2_98_1.SetPoint(uRad,Radius,AB_pp.EvalWaveFun2(WhichBin,Radius,1)*pow(Radius*AB_pp.FmNu()*98.,2));
        gTotWF2_98_0.SetPoint(uRad,Radius,AB_pp.EvalWaveFun2(WhichBin,Radius,0));
        gTotWF2_98_1.SetPoint(uRad,Radius,AB_pp.EvalWaveFun2(WhichBin,Radius,1));
    }



    TFile* OutputFile = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Temp/Output.root","recreate");

    gSPwaves.Write();
    gSwave.Write();
    gTotWF2_98_0.Write();
    gTotWF2_98_1.Write();

    delete OutputFile;
}

void plot_pXi_Ck(){
    const double kMin_pp=0;
    const double kMax_pp = 400;
    const unsigned NumMomBins_pp = 100;

    //starting value, do not worry about it too much
    const double GaussSourceSize = 1.15;
    const double Mass_p=938.272;
    const double Mass_Xim = 1321.7;

    double pXimPotParsI0S0[11]={0,0,pXim_Lattice,12,0,-1,1,0,0,0,0};
    double pXimPotParsI0S1[11]={0,0,pXim_Lattice,12,0,-1,1,1,0,1,0};
    double pXimPotParsI1S0[11]={0,0,pXim_Lattice,6,1,1,1,0,0,0,0};
    double pXimPotParsI1S1[11]={0,0,pXim_Lattice,6,1,1,1,1,0,1,0};

    CATSparameters c_pXimPotParsI0S0(CATSparameters::tPotential,9,true);
    c_pXimPotParsI0S0.SetParameters(&pXimPotParsI0S0[2]);
    CATSparameters c_pXimPotParsI0S1(CATSparameters::tPotential,9,true);
    c_pXimPotParsI0S1.SetParameters(&pXimPotParsI0S1[2]);
    CATSparameters c_pXimPotParsI1S0(CATSparameters::tPotential,9,true);
    c_pXimPotParsI1S0.SetParameters(&pXimPotParsI1S0[2]);
    CATSparameters c_pXimPotParsI1S1(CATSparameters::tPotential,9,true);
    c_pXimPotParsI1S1.SetParameters(&pXimPotParsI1S1[2]);

    CATS AB_pXi;
    double Pars_pXi[4] = {0,0,0,GaussSourceSize};
    CATSparameters cPars_pXi(CATSparameters::tSource,1,true);
    cPars_pXi.SetParameters(&Pars_pXi[3]);

    AB_pXi.SetAnaSource(GaussSource, cPars_pXi);
    AB_pXi.SetUseAnalyticSource(true);

    AB_pXi.SetExcludeFailedBins(false);
    AB_pXi.SetMomBins(NumMomBins_pp,kMin_pp,kMax_pp);

    AB_pXi.SetNumChannels(4);
    AB_pXi.SetNumPW(0,1);
    AB_pXi.SetNumPW(1,1);
    AB_pXi.SetNumPW(2,1);
    AB_pXi.SetNumPW(3,1);
    AB_pXi.SetSpin(0,0);
    AB_pXi.SetSpin(1,1);
    AB_pXi.SetSpin(2,0);
    AB_pXi.SetSpin(3,1);
    AB_pXi.SetChannelWeight(0, 1./8.);
    AB_pXi.SetChannelWeight(1, 3./8.);
    AB_pXi.SetChannelWeight(2, 1./8.);
    AB_pXi.SetChannelWeight(3, 3./8.);

    AB_pXi.SetQ1Q2(-1);
    AB_pXi.SetPdgId(2212, 3122);//same as Lambda, in case we want to use EPOS pL source
    AB_pXi.SetRedMass( (Mass_p*Mass_Xim)/(Mass_p+Mass_Xim) );

    AB_pXi.SetShortRangePotential(0,0,fDlmPot,c_pXimPotParsI0S0);
    AB_pXi.SetShortRangePotential(1,0,fDlmPot,c_pXimPotParsI0S1);
    AB_pXi.SetShortRangePotential(2,0,fDlmPot,c_pXimPotParsI1S0);
    AB_pXi.SetShortRangePotential(3,0,fDlmPot,c_pXimPotParsI1S1);
    //AB_pXi.SetNotifications(CATS::nSilent);
    //AB_pXi.SetMaxNumThreads(1);
    AB_pXi.KillTheCat();

    DLM_Ck Ck_pXi(1,0,AB_pXi);
    printf("Eval at 100: %f\n", Ck_pXi.Eval(100));
    printf("Eval at 200: %f\n", Ck_pXi.Eval(200));
    printf("Eval at 400: %f\n", Ck_pXi.Eval(400));
    printf("Eval at 401: %f\n", Ck_pXi.Eval(401));

    TGraph gSPwaves;
    gSPwaves.Set(NumMomBins_pp);
    gSPwaves.SetName("gSPwaves");
    TGraph gSwave;
    gSwave.Set(NumMomBins_pp);
    gSwave.SetName("gSwave");

    for(unsigned uBin=0; uBin<NumMomBins_pp; uBin++){
        //printf("C(%.2f)=%.3f vs %.3f\n",AB_pp.GetMomentum(uBin),AB_pp.GetCorrFun(uBin),AB_pp_2.GetCorrFun(uBin));
        gSPwaves.SetPoint(uBin,AB_pXi.GetMomentum(uBin),AB_pXi.GetCorrFun(uBin));
        //gSwave.SetPoint(uBin,AB_pp_2.GetMomentum(uBin),AB_pp_2.GetCorrFun(uBin));
    }

    TFile* OutputFile = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Temp/Output_pXi.root","recreate");

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

void RandomTest1(const double& Stability, const double& Scale){
    DLM_Random RanGen(11);
    const unsigned NumBins=256;
    //for stab 1.5 => mult. by 1.25
    //for stab 1.2 => mult. by 1.1
    //=> something wrong with the def. of stab, to get to the eq. of the paper,
    //I need to mult. the hStable3D (sigma) by factor of 1+(Stab-1)*0.5 = 0.5*(Stab+1)

    TH1F* hUniform = new TH1F("hUniform","hUniform",NumBins,-16,16);
    TH1F* hStable = new TH1F("hStable","hStable",NumBins,-16,16);
    TH1F* hGauss = new TH1F("hGauss","hGauss",NumBins,-16,16);
    TH1F* hStable3D = new TH1F("hStable3D","hStable3D",NumBins,0,16);
    hStable3D->SetLineColor(kRed);
    TH1F* hStable3D_Diff = new TH1F("hStable3D_Diff","hStable3D_Diff",NumBins,0,16);
    hStable3D_Diff->SetLineColor(kBlue);
    TH1F* hGauss3D_Diff = new TH1F("hGauss3D_Diff","hGauss3D_Diff",NumBins,0,16);
    TH1F* hCauchy3D_Diff = new TH1F("hCauchy3D_Diff","hCauchy3D_Diff",NumBins,0,16);
    TH1F* hStable3D_Ana = new TH1F("hStable3D_Ana","hStable3D_Ana",NumBins,0,16);
    hStable3D_Ana->SetLineColor(kGreen+1);
    TH1F* hStable3D_AnaInt = new TH1F("hStable3D_AnaInt","hStable3D_AnaInt",NumBins,0,16);
    const unsigned NumIter = 400000;
    for(unsigned uIter=0; uIter<NumIter; uIter++){
        //hUniform->Fill(RanGen.Uniform(-0.5,1.5));
        //hStable->Fill(RanGen.Stable(2,0,1,0));
        //hGauss->Fill(RanGen.Gauss(0,1));
        //hStable3D->Fill(RanGen.StableR(3,Stability,0,Scale*(Stability+1)*0.5,0));
        hStable3D->Fill(RanGen.StableR(3,Stability,0,Scale*sqrt(Stability)*2./Stability,0));
        //hStable3D_Diff->Fill(RanGen.StableDiffR(3,Stability,0,Scale/sqrt(2.*(3.-Stability)),0));
        hStable3D_Diff->Fill(RanGen.StableDiffR(3,Stability,0,Scale,0));
        hGauss3D_Diff->Fill(RanGen.GaussDiffR(3,0,Scale));
        hCauchy3D_Diff->Fill(RanGen.CauchyDiffR(3,0,Scale));
    }
    for(unsigned uBin=1; uBin<=NumBins; uBin++){
        double PARS[6];
        PARS[1] = hStable3D_Ana->GetBinCenter(uBin);
        //PARS[3] = Scale*2./Stability;
        //PARS[3] = Scale/sqrt(Stability);
        PARS[3] = Scale;
        PARS[4] = Stability;
        PARS[5] = 3;
        hStable3D_Ana->SetBinContent(uBin,LevySource3D_2particle(PARS));
        //hStable3D_Ana->SetBinContent(uBin,LevySource_A(PARS));
        hStable3D_Ana->SetBinError(uBin,0);
    }
    for(unsigned uBin=1; uBin<=NumBins; uBin++){
        double PARS[6];
        PARS[1] = 16;
        PARS[2] = hStable3D_AnaInt->GetBinCenter(uBin);
        PARS[3] = Scale;
        PARS[4] = Stability;
        PARS[5] = 3;
        //hStable3D_AnaInt->SetBinContent(uBin,LevyIntegral1D(PARS));

        //DLM_INT_SetFunction(LevyIntegral1D,PARS,2);
        //printf("Integral = %f\n",DLM_INT_aSimpsonWiki(0.0,32.));
    }

    hUniform->Scale(1./double(NumIter),"width");
    hStable->Scale(1./double(NumIter),"width");
    hGauss->Scale(1./double(NumIter),"width");
    hStable3D->Scale(1./double(NumIter),"width");
    hStable3D_Diff->Scale(1./double(NumIter),"width");
    hGauss3D_Diff->Scale(1./double(NumIter),"width");
    hCauchy3D_Diff->Scale(1./double(NumIter),"width");
    //hStable3D_Ana->Scale(1./double(hStable3D_Ana->Integral(1,hStable3D_Ana->GetNbinsX())),"width");

    /*
    hUniform->Scale(1./double(hUniform->Integral(1,hUniform->GetNbinsX())),"width");
    hStable->Scale(1./double(hStable->Integral(1,hStable->GetNbinsX())),"width");
    hGauss->Scale(1./double(hGauss->Integral(1,hGauss->GetNbinsX())),"width");
    hStable3D->Scale(1./double(hStable3D->Integral(1,hStable3D->GetNbinsX())),"width");
    hStable3D_Diff->Scale(1./double(hStable3D_Diff->Integral(1,hStable3D_Diff->GetNbinsX())),"width");
    hGauss3D_Diff->Scale(1./double(hGauss3D_Diff->Integral(1,hGauss3D_Diff->GetNbinsX())),"width");
    hCauchy3D_Diff->Scale(1./double(hCauchy3D_Diff->Integral(1,hCauchy3D_Diff->GetNbinsX())),"width");
    hStable3D_Ana->Scale(1./double(hStable3D_Ana->Integral(1,hStable3D_Ana->GetNbinsX())),"width");
    */
/*
    const unsigned NumMomBins_pp = 100;
    const double kMin_pp = 0;
    const double kMax_pp = 400;

    CATS AB_pp;
    AB_pp.SetMaxNumThreads(1);

    DLM_StableDistribution StabDistr;
    //StabDistr.SetStability(Stability);
    //StabDistr.SetLocation(0);
    //StabDistr.SetScale(Scale);
    //StabDistr.SetSkewness(0);
    AB_pp.SetAnaSource(MemberSourceForwarder, &StabDistr, 4);
    AB_pp.SetAnaSource(0,Stability);
    AB_pp.SetAnaSource(1,0);
    AB_pp.SetAnaSource(2,Scale);
    AB_pp.SetAnaSource(3,0);
    AB_pp.SetUseAnalyticSource(true);
    AB_pp.SetThetaDependentSource(false);

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
    AB_pp.SetChannelWeight(0, 3./12.);
    AB_pp.SetChannelWeight(1, 1./12.);
    AB_pp.SetChannelWeight(2, 3./12.);
    AB_pp.SetChannelWeight(3, 5./12.);

    const double Mass_p=938.272;
    AB_pp.SetQ1Q2(1);
    AB_pp.SetPdgId(2212, 2212);
    AB_pp.SetRedMass( 0.5*Mass_p );

    CATSparameters** cPars_pp = new CATSparameters* [4];
    for(unsigned short usCh=0; usCh<4; usCh++){
        cPars_pp[usCh] = new CATSparameters(CATSparameters::tPotential,8,true);
        cPars_pp[usCh]->SetParameter(0,NN_AV18);
        cPars_pp[usCh]->SetParameter(1,v18_Coupled3P2);
        cPars_pp[usCh]->SetParameter(2,1);
        cPars_pp[usCh]->SetParameter(3,1);
        cPars_pp[usCh]->SetParameter(4,1);
        switch(usCh){
        case 0 :
            cPars_pp[usCh]->SetParameter(5,0);
            cPars_pp[usCh]->SetParameter(6,0);
            cPars_pp[usCh]->SetParameter(7,0);
            break;
        case 1 :
            cPars_pp[usCh]->SetParameter(5,1);
            cPars_pp[usCh]->SetParameter(6,1);
            cPars_pp[usCh]->SetParameter(7,0);
            break;
        case 2 :
            cPars_pp[usCh]->SetParameter(5,1);
            cPars_pp[usCh]->SetParameter(6,1);
            cPars_pp[usCh]->SetParameter(7,1);
            break;
        case 3 :
            cPars_pp[usCh]->SetParameter(5,1);
            cPars_pp[usCh]->SetParameter(6,1);
            cPars_pp[usCh]->SetParameter(7,2);
            break;
        default : break;
        }
        AB_pp.SetShortRangePotential(usCh,usCh==0?0:1,fDlmPot,*cPars_pp[usCh]);
    }
    AB_pp.KillTheCat();

    TGraph ppCkStable3D_Diff;
    ppCkStable3D_Diff.SetName("ppCkStable3D_Diff");
    ppCkStable3D_Diff.Set(NumMomBins_pp);
    for(unsigned uBin=0; uBin<NumMomBins_pp; uBin++){
        ppCkStable3D_Diff.SetPoint(uBin,AB_pp.GetMomentum(uBin),AB_pp.GetCorrFun(uBin));
        printf("C(%.0f)=%f\n",AB_pp.GetMomentum(uBin),AB_pp.GetCorrFun(uBin));
    }
*/
    printf("Integral hGauss3D_Diff = %f\n",hGauss3D_Diff->Integral(1,hGauss3D_Diff->GetNbinsX()));
    printf("Integral hStable3D_Diff = %f\n",hStable3D_Diff->Integral(1,hStable3D_Diff->GetNbinsX()));
    printf("Integral hStable3D_Ana = %f\n",hStable3D_Ana->Integral(1,hStable3D_Ana->GetNbinsX()));
/*
    hUniform->Scale(1.,"width");
    hStable->Scale(1.,"width");
    hGauss->Scale(1.,"width");
    hStable3D->Scale(1.,"width");
    hStable3D_Diff->Scale(1.,"width");
    hGauss3D_Diff->Scale(1.,"width");
    hCauchy3D_Diff->Scale(1.,"width");
*/
    TFile* OutFile = new TFile(
        TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/RandomTest1/OutFile_%.2f_%.2f.root",Stability,Scale)
                               ,"recreate");
    hUniform->Write();
    hStable->Write();
    hGauss->Write();
    hStable3D->Write();
    hStable3D_Diff->Write();
    //ppCkStable3D_Diff.Write();
    hGauss3D_Diff->Write();
    hCauchy3D_Diff->Write();
    hStable3D_Ana->Write();
    hStable3D_AnaInt->Write();
    delete hUniform;
    delete hStable;
    delete hGauss;
    delete hStable3D;
    delete hStable3D_Diff;
    delete hGauss3D_Diff;
    delete hCauchy3D_Diff;
    delete hStable3D_Ana;
    delete hStable3D_AnaInt;
    delete OutFile;

    //delete cPars_pp[0];
    //delete cPars_pp[1];
    //delete cPars_pp[2];
    //delete cPars_pp[3];
    //delete [] cPars_pp;
}




double FitWithLevy3D(double* x, double* pars){
    double PARS[6];
    PARS[1] = *x;
    PARS[3] = pars[1];//Scale
    PARS[4] = pars[2];//Stability
    PARS[5] = 3;
    return pars[0]*LevySource3D_2particle(PARS);
}

void RandomTest2(){

    const unsigned NumStabSteps = 11;
    const double StabMin = 1;
    const double StabMax = 2;
    const double StabStep = NumStabSteps>1?(StabMax-StabMin)/double(NumStabSteps-1):0;

    const unsigned NumScaleSteps = 3;
    const double ScaleMin = NumScaleSteps>1?0.5:1;
    const double ScaleMax = 1.5;
    const double ScaleStep = NumScaleSteps>1?(ScaleMax-ScaleMin)/double(NumScaleSteps-1):0;

    DLM_Random RanGen(11);
    const unsigned NumBins=128;

    TFile* fOutput = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/RandomTest2/fOutput.root","recreate");
    TNtuple* ntInfo = new TNtuple("ntInfo", "ntInfo", "Scale:Stability:ModifScale:Chi2Ndf");

    const unsigned NumIter = 100000;
    double Stability;
    double Scale;
    Float_t Buffer[4];
    for(unsigned uScale=0; uScale<NumScaleSteps; uScale++){
        Scale = ScaleMin+double(uScale)*ScaleStep;
        for(unsigned uStab=0; uStab<NumStabSteps; uStab++){
            Stability = StabMin+double(uStab)*StabStep;
//printf("uStab=%u (%f=%f+%f*%f)\n",uStab,Stability,StabMin,double(uStab),StabStep);

            TH1F* hStable3D_Diff = new TH1F(TString::Format("hStable3D_Diff_%.2f_%.2f",Scale,Stability),
                                            TString::Format("hStable3D_Diff_%.2f_%.2f",Scale,Stability),
                                            NumBins,0,16);
            TF1* fitter = new TF1(TString::Format("fitter_%.2f_%.2f",Scale,Stability),FitWithLevy3D,0,16,3);
            fitter->FixParameter(0,1);
            fitter->SetParameter(1,Scale);
            fitter->SetParLimits(1,Scale*0.25,Scale*4.);
            fitter->FixParameter(2,Stability);
            for(unsigned uIter=0; uIter<NumIter; uIter++){
                hStable3D_Diff->Fill(RanGen.StableDiffR(3,Stability,0,Scale,0));
            }
            hStable3D_Diff->Scale(1./double(NumIter),"width");
            hStable3D_Diff->Fit(fitter,"S, N, R, M");
            Buffer[0] = Scale;
            Buffer[1] = Stability;
            Buffer[2] = fitter->GetParameter(1)/Scale;
            Buffer[3] = fitter->GetChisquare()/double(fitter->GetNDF());
            ntInfo->Fill(Buffer);
            hStable3D_Diff->Write();
            fitter->Write();
            delete hStable3D_Diff;
            delete fitter;

        }
    }

    ntInfo->Write();

    delete ntInfo;
    delete fOutput;
}

void PermutationTest(const unsigned short& Dim){
        unsigned NumPermutations = 1;
        for(unsigned short sDim=0; sDim<Dim; sDim++){
            NumPermutations *= 2;
        }
        printf("\nPermutationTest in %uD:\n",Dim);
        char** PER = new char* [NumPermutations];
        char UpdatedPer = 0;
        for(unsigned uPer=0; uPer<NumPermutations; uPer++){
            PER[uPer] = new char [Dim];
            bool Increase = true;
            for(unsigned short sDim=0; sDim<Dim; sDim++){
                if(!uPer) PER[uPer][sDim]=0;
                else{
                    PER[uPer][sDim]=PER[uPer-1][sDim];
                    if(Increase && PER[uPer][sDim]==0){
                        PER[uPer][sDim]=1;
                        Increase=false;
                    }
                    else if(Increase){
                        PER[uPer][sDim]=0;
                        Increase=true;
                    }
                }
            }
        }

        printf(" NumPermutations=%u\n",NumPermutations);
        for(unsigned uPer=0; uPer<NumPermutations; uPer++){
            //printf("  %u =",uPer);
            //for(unsigned short sDim=0; sDim<Dim; sDim++) printf(" %i",PER[uPer][sDim]);
            //printf("\n");
        }

}

void TestDimiHisto_1(){
    DLM_Histo<double> histo1;

    histo1.SetUp(3);
    unsigned NumBinsX = 4;
    double xMin = 0;
    double xMax = 20;
    unsigned NumBinsY = 4;
    double yMin = 0;
    double yMax = 40;
    unsigned NumBinsZ = 4;
    double zMin = -5;
    double zMax = 5;

    histo1.SetUp(0,NumBinsX,xMin,xMax);
    histo1.SetUp(1,NumBinsY,yMin,yMax);
    histo1.SetUp(2,NumBinsZ,zMin,zMax);
    histo1.Initialize();

    for(unsigned uBinX=0; uBinX<NumBinsX; uBinX++){
        printf("histo1.GetBinLowEdge(0,%u)=%f\n",uBinX,histo1.GetBinLowEdge(0,uBinX));
        printf("histo1.GetBinCenter(0,%u)=%f\n",uBinX,histo1.GetBinCenter(0,uBinX));
        printf("histo1.GetBinUpEdge(0,%u)=%f\n",uBinX,histo1.GetBinUpEdge(0,uBinX));
        printf("\n");
    }
    for(unsigned uBinY=0; uBinY<NumBinsY; uBinY++){
        printf("histo1.GetBinLowEdge(1,%u)=%f\n",uBinY,histo1.GetBinLowEdge(1,uBinY));
        printf("histo1.GetBinCenter(1,%u)=%f\n",uBinY,histo1.GetBinCenter(1,uBinY));
        printf("histo1.GetBinUpEdge(1,%u)=%f\n",uBinY,histo1.GetBinUpEdge(1,uBinY));
        printf("\n");
    }
    for(unsigned uBinZ=0; uBinZ<NumBinsZ; uBinZ++){
        printf("histo1.GetBinLowEdge(1,%u)=%f\n",uBinZ,histo1.GetBinLowEdge(2,uBinZ));
        printf("histo1.GetBinCenter(1,%u)=%f\n",uBinZ,histo1.GetBinCenter(2,uBinZ));
        printf("histo1.GetBinUpEdge(1,%u)=%f\n",uBinZ,histo1.GetBinUpEdge(2,uBinZ));
        printf("\n");
    }
    unsigned WhichBin[3];
    unsigned WhichBinTmp[3];
    double BinValues[3];
    double EvalAt[3];

    for(unsigned uBinX=0; uBinX<NumBinsX; uBinX++){
        for(unsigned uBinY=0; uBinY<NumBinsY; uBinY++){
            for(unsigned uBinZ=0; uBinZ<NumBinsZ; uBinZ++){
                WhichBin[0]=uBinX;
                WhichBin[1]=uBinY;
                WhichBin[2]=uBinZ;
                histo1.SetBinContent(WhichBin,2.*histo1.GetBinCenter(0,uBinX)+1.*histo1.GetBinCenter(1,uBinY)-1.5*histo1.GetBinCenter(2,uBinZ));
                unsigned TotBin = histo1.GetTotBin(WhichBin);
                histo1.GetBinCoordinates(TotBin,WhichBinTmp);
                printf("histo1.GetTotBin(%u %u %u)=%u\n",WhichBin[0],WhichBin[1],WhichBin[2],TotBin);
                printf("histo1.GetBinCoordinates(%u)=%u %u %u\n",TotBin,WhichBinTmp[0],WhichBinTmp[1],WhichBinTmp[2]);

                //printf("histo1.GetBinContent(%u %u)=%f\n",WhichBin[0],WhichBin[1],histo1.GetBinContent(WhichBin));
                //;histo1.GetBinContent(WhichBin)
                printf("histo1.GetBinContent(%u,%u,%u=%u)=%f (%f)\n",uBinX,uBinY,uBinZ,histo1.GetTotBin(WhichBin),histo1.GetBinContent(WhichBin),
                       2.*histo1.GetBinCenter(0,uBinX)+1.*histo1.GetBinCenter(1,uBinY)-1.5*histo1.GetBinCenter(2,uBinZ));
            }
        }
    }
    EvalAt[0] = 20;
    EvalAt[1] = 51;
    EvalAt[2] = -7;
    printf(" Eval(%f %f %f) = %f (%f)\n",EvalAt[0],EvalAt[1],EvalAt[2],histo1.Eval(EvalAt),2.*EvalAt[0]+1.*EvalAt[1]-1.5*EvalAt[2]);
}

void TestDimiHisto_2(){
    DLM_Histo<double> histo1;

    histo1.SetUp(2);
    unsigned NumBinsX = 3;
    double BinsX[NumBinsX+1] = {0,40,60,70};
    //double xMin = 0;
    //double xMax = 200;
    unsigned NumBinsY = 4;
    double BinsY[NumBinsY+1] = {0,1,1.5,2.5,5};
    //double yMin = 0;
    //double yMax = 2.5;

    histo1.SetUp(0,NumBinsX,BinsX);
    histo1.SetUp(1,NumBinsY,BinsY);
    histo1.Initialize();

    TRandom3 RanGen(11);
    unsigned BinId[2];
    for(unsigned uBin=0; uBin<histo1.GetNbins(); uBin++){
        histo1.SetBinContent(uBin,RanGen.Integer(1000));
    }
    histo1.ComputeError();
    printf("ComputeError\n");
    for(unsigned uBin=0; uBin<histo1.GetNbins(); uBin++){
        histo1.GetBinCoordinates(uBin,BinId);
        printf(" %u: (%u, %u) = %.3f +/- %.3f\n",uBin,BinId[0],BinId[1],histo1.GetBinContent(uBin),histo1.GetBinError(uBin));
    }
    histo1.Scale(0.5);
    printf("Scale(0.5)\n");
    for(unsigned uBin=0; uBin<histo1.GetNbins(); uBin++){
        histo1.GetBinCoordinates(uBin,BinId);
        printf(" %u: (%u, %u) = %.3f +/- %.3f\n",uBin,BinId[0],BinId[1],histo1.GetBinContent(uBin),histo1.GetBinError(uBin));
    }
    histo1.ScaleToBinSize();
    printf("ScaleToBinSize\n");
    for(unsigned uBin=0; uBin<histo1.GetNbins(); uBin++){
         printf(" histo1.GetBinSize(%u)=%f\n",uBin,histo1.GetBinSize(uBin));
    }
    for(unsigned uBin=0; uBin<histo1.GetNbins(); uBin++){
        histo1.GetBinCoordinates(uBin,BinId);
        printf(" %u: (%u, %u) = %.3f +/- %.3f\n",uBin,BinId[0],BinId[1],histo1.GetBinContent(uBin),histo1.GetBinError(uBin));
    }
    histo1.RescaleAxis(0,0.001,true);
    printf("RescaleAxis(0,0.001)\n");
    for(unsigned uBin=0; uBin<histo1.GetNbins(); uBin++){
        histo1.GetBinCoordinates(uBin,BinId);
        printf(" %u: (%u, %u) = %.3f +/- %.3f\n",uBin,BinId[0],BinId[1],histo1.GetBinContent(uBin),histo1.GetBinError(uBin));
    }
    histo1.RescaleAxis(1,1000.,true);
    printf("RescaleAxis(1,1000.)\n");
    for(unsigned uBin=0; uBin<histo1.GetNbins(); uBin++){
        histo1.GetBinCoordinates(uBin,BinId);
        printf(" %u: (%u, %u) = %.3f +/- %.3f\n",uBin,BinId[0],BinId[1],histo1.GetBinContent(uBin),histo1.GetBinError(uBin));
    }

    DLM_Histo<double> histo2(histo1);
    printf("histo2\n");
    histo2.Scale(2);
    for(unsigned uBin=0; uBin<histo2.GetNbins(); uBin++){
        histo2.GetBinCoordinates(uBin,BinId);
        printf(" %u: (%u, %u) = %.3f +/- %.3f\n",uBin,BinId[0],BinId[1],histo2.GetBinContent(uBin),histo2.GetBinError(uBin));
    }
    DLM_Histo<double> histo3=histo2-histo1;
    printf("histo3=histo2+histo1\n");
    printf("histo3=%p\n",&histo3);
    for(unsigned uBin=0; uBin<histo3.GetNbins(); uBin++){
        histo3.GetBinCoordinates(uBin,BinId);
        printf(" %u: (%u, %u) = %.3f +/- %.3f\n",uBin,BinId[0],BinId[1],histo3.GetBinContent(uBin),histo3.GetBinError(uBin));
    }
    histo3.AddToAll(3);
    histo3/=histo2;
    printf("histo3/=histo2\n");
    for(unsigned uBin=0; uBin<histo3.GetNbins(); uBin++){
        histo3.GetBinCoordinates(uBin,BinId);
        printf(" %u: (%u, %u) = %.3f +/- %.3f\n",uBin,BinId[0],BinId[1],histo3.GetBinContent(uBin),histo3.GetBinError(uBin));
    }
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

    CALL_BERNIE_AND_VALE();
    //LamLamMcExclusionMAIN(argc,ARGV);

    //plot_pp_Ck();
    //Test_pL_1();
    //plot_pXi_Ck();
    //Test_Tom_pXi();
    //Test_Gamow_pXi();
    //pXi_PhaseShifts(11);
    //pXi_PhaseShifts(12);
    //pXi_PhaseShifts(-1);

    //UsmaniPlusExpCore();

    //GerhardMAIN(argc,argv);
    //RandomTest1(1.00,1.00);
    //RandomTest1(1.01,1.00);
    //RandomTest1(1.2,1.0);
    //RandomTest1(1.4,1.0);
    //RandomTest1(1.6,1.0);
    //RandomTest1(1.8,1.0);
    //RandomTest1(2.0-1e-6,1.0);
    //RandomTest2();
    //CompareDoubleGaussToGaussExp();

    //PermutationTest(1);
    //PermutationTest(2);
    //PermutationTest(4);
    //PermutationTest(8);
    //PermutationTest(16);
    //PermutationTest(24);

    //TestDimiHisto_1();
    //TestDimiHisto_2();
    //TestCleverLevy1();

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
