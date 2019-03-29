#include "ForGerhard.h"

#include <stdio.h>

#include "CATS.h"
#include "DLM_Potentials.h"
#include "DLM_Source.h"
#include "DLM_CkDecomposition.h"
#include "DLM_CkModels.h"
#include "DLM_Fitters.h"
#include "DLM_Histo.h"

#include "TString.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TFile.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TCanvas.h"

const bool DEBUG_MODE=true;
const double Mass_p=938.272;
const double Mass_L=1115.683;
const double PI = 3.141592653589793;

void EposFit_pp_1(){

    const TString OutputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/ForGerhard/";

    //there are 2 CATS objects, one for the fitting (4 MeV bin width) and one for the differential look into the source (40 MeV bin width)
    const unsigned MaxNumPairs[2] = {unsigned(256e6),unsigned(16e6)};//maximum number of pairs to iterate over in the input file
    const unsigned MaxNumPairsPerBin[2] = {unsigned(16e3),unsigned(16e3)};//maximum number of pairs to be loaded into each momentum bin
    const short MixingDepth = 8;//!put to one for the actual source. To enhance statistics at the expense of 5-10% uncertainty, put to 8
    const unsigned NumMomBins_pp[2] = {94,10};
    //const unsigned NumMomBins_pp = 30;
    const unsigned kMin_pp[2] = {0,0};
    const unsigned kMax_pp[2] = {NumMomBins_pp[0]*4, NumMomBins_pp[1]*40};
    //-1 is EPOS
    //0 is Gauss
    //1 is Stable distribution
    const int SourceType=1;
    const double GaussSourceSize = 1.15;
    const bool AddSystematics = false;
    const bool InludeBaseline = false;
    const bool AllowSourceRescaling = true;

    //0 => the first analysis
    //1 => as we did it for the LL paper
    const int DataFlag = 1;

    const TString DataPeriod = "pp13TeV";
    TString ppDataFileName;
    TString ppResolutionFile;
    TString ppFeedDownFile;
    // /scratch6/dmihaylov/EPOS_FILE_READER/OutputFiles/f19/Scratch9_OSCAR1997_AllLrzProtonFiles.f19
    TString EPOSfile_pp;
    TString SystErrFileName_pp;

    if(DataPeriod=="pp13TeV"){
        ppDataFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/January2019_PaperDrafts/CFOutput_pp.root";
        ppResolutionFile = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/ResolutionMatrices/Sample6_MeV_compact.root";
        ppFeedDownFile = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/DecayMatrices/run2_decay_matrices_old.root";
        //EPOSfile_pp = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/Scratch9_OSCAR1997_AllLrzProtonFiles.f19";
        //EPOSfile_pp = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/7TeV/pp_nist0_32Ki.f19";
        //EPOSfile_pp = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/EPOS_LBF_pp90/pp_nist0_256Ki.f19";
        EPOSfile_pp = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/EPOS_LBF_pp91/pp_nist0_256Ki.f19";
        //EPOSfile_pp = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/EPOS_LBF_pp91/pp_nist0_seco_256Ki.f19";
        SystErrFileName_pp = "/home/dmihaylov/CernBox/Femto_pp13/C2totalsysPP.root";
    }
    else{
        printf("\033[1;31mERROR:\033[0m %s is an unknown data set\n",DataPeriod.Data());
        return;
    }

    const unsigned short NumChannels_pp = 4;
    CATSparameters*** cPars_pp = new CATSparameters** [2];
    CATSparameters** pp_Source = new CATSparameters* [2];

    CATS* AB_pp = new CATS [2];
    DLM_StableDistribution* StabDistr = new DLM_StableDistribution[2];
    for(short sCat=0; sCat<2; sCat++){
        pp_Source[sCat] = new CATSparameters(CATSparameters::tSource,1,true);
        pp_Source[sCat]->SetParameter(0, GaussSourceSize);
        AB_pp[sCat].SetMaxNumThreads(1);
        if(SourceType==-1){
            AB_pp[sCat].SetUseAnalyticSource(false);
            AB_pp[sCat].SetInputFileName(EPOSfile_pp.Data());
            AB_pp[sCat].SetMaxPairsToRead(MaxNumPairs[sCat]);
            AB_pp[sCat].SetMaxPairsPerBin(MaxNumPairsPerBin[sCat]);
            AB_pp[sCat].SetMixingDepth(MixingDepth);
            AB_pp[sCat].SetThetaDependentSource(false);
            AB_pp[sCat].SetTransportRenorm(1.0);
            AB_pp[sCat].SetTauCorrection(false);
        }
        else if(SourceType==0){
            AB_pp[sCat].SetAnaSource(GaussSource, *pp_Source[sCat]);
            AB_pp[sCat].SetUseAnalyticSource(true);
            AB_pp[sCat].SetThetaDependentSource(false);
        }
        else if(SourceType==1){
            //StabDistr[sCat].SetStability(1.5);
            //StabDistr[sCat].SetLocation(0);
            //StabDistr[sCat].SetScale(1.0);
            //StabDistr[sCat].SetSkewness(0);
            AB_pp[sCat].SetAnaSource(CatsSourceForwarder, &StabDistr[sCat], 4);
            AB_pp[sCat].SetAnaSource(0,1.5,false);
            AB_pp[sCat].SetAnaSource(1,0,false);
            AB_pp[sCat].SetAnaSource(2,1.0,false);
            AB_pp[sCat].SetAnaSource(3,0,false);
            AB_pp[sCat].SetUseAnalyticSource(true);
            AB_pp[sCat].SetThetaDependentSource(false);
        }
        else{
            printf("\033[1;31mERROR:\033[0m %u is an unknown source type\n",SourceType);
        }

        AB_pp[sCat].SetExcludeFailedBins(false);
        AB_pp[sCat].SetMomBins(NumMomBins_pp[sCat],kMin_pp[sCat],kMax_pp[sCat]);
        AB_pp[sCat].SetNumChannels(NumChannels_pp);
        AB_pp[sCat].SetNumPW(0,2);
        AB_pp[sCat].SetNumPW(1,2);
        AB_pp[sCat].SetNumPW(2,2);
        AB_pp[sCat].SetNumPW(3,2);
        AB_pp[sCat].SetSpin(0,0);
        AB_pp[sCat].SetSpin(1,1);
        AB_pp[sCat].SetSpin(2,1);
        AB_pp[sCat].SetSpin(3,1);
        AB_pp[sCat].SetChannelWeight(0, 3./12.);
        AB_pp[sCat].SetChannelWeight(1, 1./12.);
        AB_pp[sCat].SetChannelWeight(2, 3./12.);
        AB_pp[sCat].SetChannelWeight(3, 5./12.);

        AB_pp[sCat].SetQ1Q2(1);
        AB_pp[sCat].SetPdgId(2212, 2212);
        AB_pp[sCat].SetRedMass( 0.5*Mass_p );

        cPars_pp[sCat] = new CATSparameters* [NumChannels_pp];
        for(unsigned short usCh=0; usCh<NumChannels_pp; usCh++){
            cPars_pp[sCat][usCh] = new CATSparameters(CATSparameters::tPotential,8,true);

            cPars_pp[sCat][usCh]->SetParameter(0,NN_AV18);
            cPars_pp[sCat][usCh]->SetParameter(1,v18_Coupled3P2);
            cPars_pp[sCat][usCh]->SetParameter(2,1);
            cPars_pp[sCat][usCh]->SetParameter(3,1);
            cPars_pp[sCat][usCh]->SetParameter(4,1);

            switch(usCh){
            case 0 :
                cPars_pp[sCat][usCh]->SetParameter(5,0);
                cPars_pp[sCat][usCh]->SetParameter(6,0);
                cPars_pp[sCat][usCh]->SetParameter(7,0);
                break;
            case 1 :
                cPars_pp[sCat][usCh]->SetParameter(5,1);
                cPars_pp[sCat][usCh]->SetParameter(6,1);
                cPars_pp[sCat][usCh]->SetParameter(7,0);
                break;
            case 2 :
                cPars_pp[sCat][usCh]->SetParameter(5,1);
                cPars_pp[sCat][usCh]->SetParameter(6,1);
                cPars_pp[sCat][usCh]->SetParameter(7,1);
                break;
            case 3 :
                cPars_pp[sCat][usCh]->SetParameter(5,1);
                cPars_pp[sCat][usCh]->SetParameter(6,1);
                cPars_pp[sCat][usCh]->SetParameter(7,2);
                break;
            default : break;
            }

            AB_pp[sCat].SetShortRangePotential(usCh,usCh==0?0:1,fDlmPot,*cPars_pp[sCat][usCh]);
        }
        AB_pp[sCat].KillTheCat();
    }



    TH1F* hSource=NULL;
    TH1F** hkSource=NULL;
    const double rMin=0;
    const double rMax=16;
    unsigned NumPairs[2];
    NumPairs[0]=0; NumPairs[1]=0;
    if(SourceType==-1){
        NumPairs[0] = AB_pp[0].GetNumPairs();
        NumPairs[1] = AB_pp[1].GetNumPairs();
        unsigned NumRadPts = 1024;
        hSource = new TH1F("hSource","hSource",NumRadPts,rMin,rMax);
        hkSource = new TH1F* [NumMomBins_pp[1]];
        for(unsigned uMomBin=0; uMomBin<NumMomBins_pp[1]; uMomBin++){
            hkSource[uMomBin] = new TH1F(TString::Format("hkSource_%.0f",AB_pp[1].GetMomentum(uMomBin)),
                                        TString::Format("hkSource_%.0f",AB_pp[1].GetMomentum(uMomBin)),
                                        NumRadPts/8,rMin,rMax);
        }
        double RelMom;
        double RelPos;
        double RelCosTh;
        double TotMom;
        unsigned WhichMomBin;
        for(unsigned uPair=0; uPair<NumPairs[0]; uPair++){
            AB_pp[0].GetPairInfo(uPair,RelMom,RelPos,RelCosTh,TotMom);
            hSource->Fill(RelPos);
        }
        for(unsigned uPair=0; uPair<NumPairs[1]; uPair++){
            AB_pp[1].GetPairInfo(uPair,RelMom,RelPos,RelCosTh,TotMom);
            WhichMomBin = AB_pp[1].GetMomBin(RelMom);
            hkSource[WhichMomBin]->Fill(RelPos);
        }
    }
    if(hSource) hSource->Scale(1./double(NumPairs[0]),"width");
    for(unsigned uMomBin=0; uMomBin<NumMomBins_pp[1]; uMomBin++){
        if(!hkSource) break;
        hkSource[uMomBin]->Scale(1./double(hkSource[uMomBin]->Integral(1,hkSource[uMomBin]->GetNbinsX())),"width");
    }
/*
_91
In those files I found 8204200 events!
In those files I found 6.251389e+09 particles!
I have selected 958804 particles with PDGid=2212

_91_Prim
In those files I found 6965100 events!
In those files I found 5.878838e+09 particles!
I have selected 142032 particles with PDGid=2212

_91_Seco
In those files I found 7420800 events!
In those files I found 6.056624e+09 particles!
I have selected 416018 particles with PDGid=2212

btw, this means c.a. 26.7% of primary particles!


_91_Seco_2
In those files I found 8177000 events!
In those files I found 6.231744e+09 particles!
I have selected 458159 particles with PDGid=2212

*/

    if(DEBUG_MODE){
        for(unsigned uMomBin=0; uMomBin<NumMomBins_pp[0]; uMomBin++){
            printf("C(%.0f) = %.3f +/- %.3f\n",AB_pp[0].GetMomentum(uMomBin),AB_pp[0].GetCorrFun(uMomBin),AB_pp[0].GetCorrFunErr(uMomBin));
        }
    }

    double PurityProton;
    double PurityLambda;
    double PurityXi;
    double pp_f0;
    double pp_f1;
    double pL_f0;
    double pL_f1;
    double pL_f2;

    //(single particle quantities)
    if(DataPeriod=="pp13TeV"){
        //as used initially
        if(DataFlag==0){
            PurityProton = 0.991213;//991213
            PurityLambda = 0.965964;
            PurityXi = 0.956;

            pp_f0 = 0.874808;
            pp_f1 = 0.0876342;//fraction of Lambda

            pL_f0 = 0.619493;//fraction of primary Lambdas
            pL_f1 = 0.206498;//fraction of Sigma0
            pL_f2 = 0.0870044;//fractions of Xi0/m
        }
        else{
            PurityProton = 0.989859;
            PurityLambda = 0.96768;
            PurityXi = 0.956;//is this correct?

            pp_f0 = 0.87397;
            pp_f1 = 0.0882211;//fraction of Lambda

            pL_f0 = 0.601008;//fraction of primary Lambdas
            pL_f1 = 0.200336;//fraction of Sigma0
            pL_f2 = 0.0993283;//fractions of Xi0/m
        }
    }
    else if(DataPeriod=="pPb5TeV"){
        if(DataFlag==0){
            PurityProton = 0.984177 ;//pPb 5 TeV
            PurityLambda = 0.93171;
            PurityXi = 0.9;//new cuts

            pp_f0 = 0.862431;
            pp_f1 = 0.0962984;

            pL_f0 = 0.744516;//fraction of primary Lambdas
            pL_f1 = 0.248172;//fraction of Sigma0
            pL_f2 = 0.003656;//fractions of Xi0/m



            PurityProton = 0.984265 ;//pPb 5 TeV
            PurityLambda = 0.937761;
            PurityXi = 0.88;//new cuts

            pp_f0 = 0.862814;
            pp_f1 = 0.09603;

            pL_f0 = 0.521433;
            pL_f1 = 0.173811;
            pL_f2 = 0.152378;
        }
        else{
            //from 20th Aug 2018
            PurityProton = 0.984265 ;//pPb 5 TeV
            PurityLambda = 0.937761;
            PurityXi = 0.88;//new cuts

            pp_f0 = 0.862814;
            pp_f1 = 0.09603;

            pL_f0 = 0.521433;
            pL_f1 = 0.173811;
            pL_f2 = 0.152378;
        }
    }
    else{
        printf("This should never happen!\n");
        return;
    }

//DO NOT CHANGE***
    double ProtonPrim = pp_f0;
    double arrayPercLamProton[3]={pp_f1/(1.-pp_f0)*0.8,pp_f1/(1.-pp_f0),pp_f1/(1.-pp_f0)*1.2};//+/- 20%

    const double SigLambdaPrimDir = pL_f0+pL_f1;
    double arrayPercSigLambda[3]={0.8*pL_f1/pL_f0,pL_f1/pL_f0,1.2*pL_f1/pL_f0};//1/3 +/- 20%
    double arrayPercXiLambda[3]={pL_f2/(1.-pL_f0-pL_f1)*0.8,pL_f2/(1.-pL_f0-pL_f1),pL_f2/(1.-pL_f0-pL_f1)*1.2};//+/- 20%

    //ratio Xi-(1530) to Xi-
    const double Xim1530_to_Xim = 0.32*(1./3.);
    //ratio Xi0(1530) to Xi0 (n=neutral)
    const double Xin1530_to_Xim = 0.32*(2./3.);
    const double Omegam_to_Xim = 0.1;
    const double OmegamXim_BR = 0.086;


    //following my lambda pars with the 3 possible modifications
    //for the proton:
    //0 = primary
    //1 = from Lambda
    //2 = other feeddown (flat)
    //3 = missidentified
    const unsigned NumChannels_p = 4;
    double** Purities_p = new double* [3];
    double** Fraction_p = new double* [3];
    for(unsigned uVar=0; uVar<3; uVar++){
        Purities_p[uVar] = new double [NumChannels_p];
        Fraction_p[uVar] = new double [NumChannels_p];

        Purities_p[uVar][0] = PurityProton;
        Purities_p[uVar][1] = PurityProton;
        Purities_p[uVar][2] = PurityProton;
        Purities_p[uVar][3] = 1.-PurityProton;

        Fraction_p[uVar][0] = ProtonPrim;
        Fraction_p[uVar][1] = (1.-ProtonPrim)*(arrayPercLamProton[uVar]);
        Fraction_p[uVar][2] = (1.-ProtonPrim)*(1.-arrayPercLamProton[uVar]);
        Fraction_p[uVar][3] = 1.;
    }

    //for the Lambda:
    //0 = primary
    //1 = from Sigma0
    //2 = from Xim
    //3 = from Xi0
    //4 = missidentified
    const unsigned NumChannels_L = 5;
    double*** Purities_L = new double** [3];
    double*** Fraction_L = new double** [3];
    for(unsigned uVarSL=0; uVarSL<3; uVarSL++){
        Purities_L[uVarSL] = new double* [3];
        Fraction_L[uVarSL] = new double* [3];
        for(unsigned uVarXi=0; uVarXi<3; uVarXi++){
            Purities_L[uVarSL][uVarXi] = new double [NumChannels_L];
            Fraction_L[uVarSL][uVarXi] = new double [NumChannels_L];

            Purities_L[uVarSL][uVarXi][0] = PurityLambda;
            Purities_L[uVarSL][uVarXi][1] = PurityLambda;
            Purities_L[uVarSL][uVarXi][2] = PurityLambda;
            Purities_L[uVarSL][uVarXi][3] = PurityLambda;
            Purities_L[uVarSL][uVarXi][4] = 1.-PurityLambda;

            //the array is r = S/L, and S+L=1 are the fractions of Sigmas and Lambdas
            double FracOfLambda = 1./(1.+arrayPercSigLambda[uVarSL]);
            Fraction_L[uVarSL][uVarXi][0] = SigLambdaPrimDir*FracOfLambda;
            Fraction_L[uVarSL][uVarXi][1] = SigLambdaPrimDir*(1.-FracOfLambda);
            Fraction_L[uVarSL][uVarXi][2] = (1.-SigLambdaPrimDir)*(arrayPercXiLambda[uVarXi]);
            Fraction_L[uVarSL][uVarXi][3] = (1.-SigLambdaPrimDir)*(1.-arrayPercXiLambda[uVarXi]);
            Fraction_L[uVarSL][uVarXi][4] = 1.;
        }
    }

    //for the Xi:
    //0 = primary
    //1 = from Xi-(1530)
    //2 = from Xi0(1530)
    //3 = from Omega
    //4 = missidentified
    const unsigned NumChannels_Xim = 5;
    double** Purities_Xim = new double* [3];
    double** Fraction_Xim = new double* [3];
    for(unsigned uVar=0; uVar<3; uVar++){
        Purities_Xim[uVar] = new double [NumChannels_Xim];
        Fraction_Xim[uVar] = new double [NumChannels_Xim];

        Purities_Xim[uVar][0] = PurityXi;
        Purities_Xim[uVar][1] = PurityXi;
        Purities_Xim[uVar][2] = PurityXi;
        Purities_Xim[uVar][3] = PurityXi;
        Purities_Xim[uVar][4] = 1.-PurityXi;

        //the ratios that we have for Xis are referred to the total number of Xi particles (which already include all contributions)
        //hence Xi1530_to_Xi indeed is simply the number of Xis that stem from a Xi1530
        Fraction_Xim[uVar][0] = 1.-Xim1530_to_Xim-Xin1530_to_Xim-Omegam_to_Xim*OmegamXim_BR;
        Fraction_Xim[uVar][1] = Xim1530_to_Xim;
        Fraction_Xim[uVar][2] = Xin1530_to_Xim;
        Fraction_Xim[uVar][3] = Omegam_to_Xim*OmegamXim_BR;
        Fraction_Xim[uVar][4] = 1.;
    }
//***

    if(DEBUG_MODE) {printf("Purities computed\n");}

//ADVANCED***
    //!DO NOT TOUCH UNLESS YOU CHANGE THE ABOVE FILES, IN WHICH CASE ASK DIMI FOR HELP
    //1/FractionOfBins th number of original bins of the correction matrix are taken into account
    //originally: 1000 MeV => 1/2 should be good most of the time
    const int Fraction_Res = 2;
    const int Fraction_Sig = 1;
    const double UnitConv_Res = 1;
    const double UnitConv_Sig = 1;

    TH2F* hFeed_pp_pL;
    TH2F* hReso_pp;

    TFile* FileFeed = new TFile(ppFeedDownFile, "read");
    TFile* FileReso = new TFile(ppResolutionFile, "read");
printf("FileFeed=%p; FileReso=%p\n",FileFeed,FileReso);
    //for the Xi we will make the assumption that all residuals are flat since we do not know better
    FileFeed->cd();
    hFeed_pp_pL = (TH2F*)FileFeed->Get("hRes_pp_pL");
printf("hFeed_pp_pL=%p\n",hFeed_pp_pL);
    FileReso->cd();

    hReso_pp = (TH2F*)FileReso->Get("hSigmaMeV_Proton_Proton");
printf("hReso_pp=%p\n",hReso_pp);
    TH2F* hFeed_pp_pL_MeV = new TH2F("hFeed_pp_pL_MeV", "hFeed_pp_pL_MeV",
        hFeed_pp_pL->GetNbinsX()/Fraction_Res, hFeed_pp_pL->GetXaxis()->GetBinLowEdge(1)*UnitConv_Res,hFeed_pp_pL->GetXaxis()->GetBinUpEdge(hFeed_pp_pL->GetNbinsX()/Fraction_Res)*UnitConv_Res,
        hFeed_pp_pL->GetNbinsY()/Fraction_Res, hFeed_pp_pL->GetYaxis()->GetBinLowEdge(1)*UnitConv_Res,hFeed_pp_pL->GetXaxis()->GetBinUpEdge(hFeed_pp_pL->GetNbinsY()/Fraction_Res)*UnitConv_Res);

    TH2F* hReso_pp_MeV = new TH2F("hReso_pp_MeV", "hReso_pp_MeV",
        hReso_pp->GetNbinsX()/Fraction_Sig, hReso_pp->GetXaxis()->GetBinLowEdge(1)*UnitConv_Sig,hReso_pp->GetXaxis()->GetBinUpEdge(hReso_pp->GetNbinsX()/Fraction_Sig)*UnitConv_Sig,
        hReso_pp->GetNbinsY()/Fraction_Sig, hReso_pp->GetYaxis()->GetBinLowEdge(1)*UnitConv_Sig,hReso_pp->GetXaxis()->GetBinUpEdge(hReso_pp->GetNbinsY()/Fraction_Sig)*UnitConv_Sig);

    for(int iBinX=1; iBinX<=hFeed_pp_pL->GetNbinsX()/Fraction_Res; iBinX++){
        for(int iBinY=1; iBinY<=hFeed_pp_pL->GetNbinsY()/Fraction_Res; iBinY++){
            hFeed_pp_pL_MeV->SetBinContent(iBinX, iBinY, hFeed_pp_pL->GetBinContent(iBinX, iBinY));
        }
    }

    for(int iBinX=1; iBinX<=hReso_pp->GetNbinsX()/Fraction_Sig; iBinX++){
        for(int iBinY=1; iBinY<=hReso_pp->GetNbinsY()/Fraction_Sig; iBinY++){
            hReso_pp_MeV->SetBinContent(iBinX, iBinY, hReso_pp->GetBinContent(iBinX, iBinY));
        }
    }

    //! DATA FILE
    TString DataHistoName_pp = "hCk_ReweightedMeV_0";
    TFile* DataFile_pp = ppDataFileName!=""?new TFile(ppDataFileName, "read"):NULL;
printf("DataFile_pp=%p\n",DataFile_pp);
    TH1F* DataHisto_pp = DataFile_pp?(TH1F*)DataFile_pp->Get(DataHistoName_pp):NULL;
printf("DataHisto_pp=%p\n",DataHisto_pp);
    TString SystErrHistoName_pp = "C2totalsysPP";
    TFile* SystErrFile_pp = SystErrFileName_pp!=""?new TFile(SystErrFileName_pp, "read"):NULL;
    TH1F* SystErrHisto_pp = SystErrFile_pp?(TH1F*)SystErrFile_pp->Get(SystErrHistoName_pp):NULL;
    int NumSEB_pp = SystErrHisto_pp==NULL?0:SystErrHisto_pp->GetNbinsX();
    if(NumSEB_pp>DataHisto_pp->GetNbinsX()) NumSEB_pp=DataHisto_pp->GetNbinsX();
    for(int iBin=0; iBin<NumSEB_pp; iBin++){
        if(AddSystematics) continue;
        DataHisto_pp->SetBinError(iBin+1, sqrt(pow(DataHisto_pp->GetBinError(iBin+1),2.)+pow(SystErrHisto_pp->GetBinContent(iBin+1),2.)));
    }

    if(DEBUG_MODE) {printf("Data loaded\n");}

    unsigned NumSourcePars;
    //this is here if it needs to be expanded for further cases
    if(SourceType==-1) NumSourcePars = 1;
    else if(SourceType==0) NumSourcePars = 1;
    else if(SourceType==1) NumSourcePars = 4;
    else NumSourcePars = 1;

    //this way you define a correlation function using a CATS object.
    //needed inputs: num source/pot pars, CATS obj
    DLM_Ck* Ck_pp = new DLM_Ck(NumSourcePars,0,AB_pp[0]);
    //for now pL is only modeled with Lednicky (NLO)
    DLM_Ck* Ck_pL = new DLM_Ck(1,4,NumMomBins_pp[0],kMin_pp[0],kMax_pp[0],Lednicky_SingletTriplet);
//printf("?\n");
    Ck_pL->SetSourcePar(0,GaussSourceSize);
//printf("Answer!\n");
    const int vMod_pL = 1;//1 is NLO, 2 is LO
    if(vMod_pL==1){
        Ck_pL->SetPotPar(0,2.91);
        Ck_pL->SetPotPar(1,2.78);
        Ck_pL->SetPotPar(2,1.54);
        Ck_pL->SetPotPar(3,2.72);
    }
    else if(vMod_pL==2){
        Ck_pL->SetPotPar(0,1.91);
        Ck_pL->SetPotPar(1,1.4);
        Ck_pL->SetPotPar(2,1.23);
        Ck_pL->SetPotPar(3,2.13);
    }
    Ck_pp->Update();
    Ck_pL->Update();

    DLM_CkDecomposition CkDec_pp("pp",3,*Ck_pp,hReso_pp_MeV);
    DLM_CkDecomposition CkDec_pL("pLambda",2,*Ck_pL,NULL);

    double lam_pp =   Purities_p[0][0]*Fraction_p[0][0]*
                            Purities_p[0][0]*Fraction_p[0][0];
    double lam_pp_pL =    Purities_p[0][0]*Fraction_p[0][0]*
                                Purities_p[0][1]*Fraction_p[0][1]*2;
    double lam_pp_fake =  Purities_p[0][3]*Purities_p[0][0]+
                                Purities_p[0][0]*Purities_p[0][3]+
                                Purities_p[0][3]*Purities_p[0][3];

    int vFrac_pp_pL = 1;
    int vFrac_pL_pSigma0 = 1;
    int vFrac_pL_pXim = 1;
    double lam_pL =   Purities_p[vFrac_pp_pL][0]*Fraction_p[vFrac_pp_pL][0]*
                            Purities_L[vFrac_pL_pSigma0][vFrac_pL_pXim][0]*Fraction_L[vFrac_pL_pSigma0][vFrac_pL_pXim][0];
    double lam_pL_pS0 =   Purities_p[vFrac_pp_pL][0]*Fraction_p[vFrac_pp_pL][0]*
                                Purities_L[vFrac_pL_pSigma0][vFrac_pL_pXim][1]*Fraction_L[vFrac_pL_pSigma0][vFrac_pL_pXim][1];
    double lam_pL_pXm =   Purities_p[vFrac_pp_pL][0]*Fraction_p[vFrac_pp_pL][0]*
                                Purities_L[vFrac_pL_pSigma0][vFrac_pL_pXim][2]*Fraction_L[vFrac_pL_pSigma0][vFrac_pL_pXim][2];
    double lam_pL_fake =  Purities_p[vFrac_pp_pL][3]*Purities_L[vFrac_pL_pSigma0][vFrac_pL_pXim][0]+
                                Purities_p[vFrac_pp_pL][0]*Purities_L[vFrac_pL_pSigma0][vFrac_pL_pXim][4]+
                                Purities_p[vFrac_pp_pL][3]*Purities_L[vFrac_pL_pSigma0][vFrac_pL_pXim][4];
    double lam_pL_AllFeed = 1.-lam_pL-lam_pL_fake;

    printf("Cross check the lambda parameters:\n");
    printf(" lam_pp = %.3f\n", lam_pp);
    printf(" lam_pp_pL = %.3f\n", lam_pp_pL);
    printf(" lam_pp_fake = %.3f\n", lam_pp_fake);
    printf("\n");

    printf("lam_pL=%.3f\n",lam_pL);
    //printf("lam_pL_pS0=%.3f\n",lam_pL_pS0);
    //printf("lam_pL_pXm=%.3f\n",lam_pL_pXm);
    printf("lam_pL_AllFeed=%.3f\n",lam_pL_AllFeed);
    printf("lam_pL_fake=%.3f\n",lam_pL_fake);
    printf("\n");

    CkDec_pp.AddContribution(0,lam_pp_pL,DLM_CkDecomposition::cFeedDown,&CkDec_pL,hFeed_pp_pL_MeV);
    CkDec_pp.AddContribution(1,1.-lam_pp-lam_pp_pL-lam_pp_fake,DLM_CkDecomposition::cFeedDown);
    CkDec_pp.AddContribution(2,lam_pp_fake,DLM_CkDecomposition::cFake);

    CkDec_pL.AddContribution(0,lam_pL_AllFeed,DLM_CkDecomposition::cFeedDown);
    CkDec_pL.AddContribution(1,lam_pL_fake,DLM_CkDecomposition::cFake);

    DLM_Fitter1* fitter = new DLM_Fitter1(1);

    fitter->SetSystem(0,*DataHisto_pp,1,CkDec_pp,kMin_pp[0],kMax_pp[0],kMax_pp[0],kMax_pp[0]);
    fitter->SetSeparateBL(0,false);

    fitter->SetParameter("pp",DLM_Fitter1::p_a,1.0,0.7,1.3);
    if(InludeBaseline) fitter->SetParameter("pp",DLM_Fitter1::p_b,1e-4,-2e-3,2e-3);
    else fitter->FixParameter("pp",DLM_Fitter1::p_b,0);
    fitter->FixParameter("pp",DLM_Fitter1::p_c,0);
    fitter->FixParameter("pp",DLM_Fitter1::p_Cl,-1);

    CkDec_pp.Update();
    printf("Fitting in progress...\n");
    const unsigned NumFitIter_r = 17;
    const unsigned NumFitIter_alpha = 17;
    const double RadiusFrom = 0.6;
    const double RadiusTo = 1.4;
    const double AlphaFrom = 0.8;
    const double AlphaTo = 2.0;

    const double StepSizeRadius = NumFitIter_r>1?(RadiusTo-RadiusFrom)/double(NumFitIter_r-1):0.001;
    const double StepSizeAlpha = NumFitIter_alpha>1?(AlphaTo-AlphaFrom)/double(NumFitIter_alpha-1):0.001;
    double CurrentRadius;
    double CurrentAlpha;
    double BestChi2Ndf=1e6;
    double BestChi2;
    const double BestExpectedChi2=115;
    const double MaxAllowedChi2=BestExpectedChi2+2.3;
    double BestPval=0;
    double BestRadius=0;
    double BestAlpha=0;
    TGraph gChi2Ndf;
    gChi2Ndf.SetName("gChi2Ndf");
    gChi2Ndf.Set(NumFitIter_r);
    TFile* OutputFile = new TFile(OutputFolder+"EposFit_pp_1.root","recreate");
    TH2F* hChi2Map = new TH2F("hChi2Map","hChi2Map",NumFitIter_r,RadiusFrom-StepSizeRadius*0.5,RadiusTo+StepSizeRadius*0.5,
                              NumFitIter_alpha,AlphaFrom-StepSizeAlpha*0.5,AlphaTo+StepSizeAlpha*0.5);

    for(unsigned uFit_r=0; uFit_r<NumFitIter_r; uFit_r++){
        //printf("uIter=%u(%u)\n",uFit_r,NumFitIter_r);
        CurrentRadius = RadiusFrom+double(uFit_r)*StepSizeRadius;
        for(unsigned uFit_alpha=0; uFit_alpha<NumFitIter_alpha; uFit_alpha++){
            CurrentAlpha = AlphaFrom+double(uFit_alpha)*StepSizeAlpha;
            printf("r=%f; α=%f\n",CurrentRadius,CurrentAlpha);
            if(SourceType==-1){
                if(AllowSourceRescaling) fitter->SetParameter("pp",DLM_Fitter1::p_sor0,1.,1.,2.);
                else fitter->FixParameter("pp",DLM_Fitter1::p_sor0,CurrentRadius);
                //fitter->FixParameter("pp",DLM_Fitter1::p_sor0,CurrentRadius);
            }
            else if(SourceType==0){
                fitter->SetParameter("pp",DLM_Fitter1::p_sor0,(RadiusTo+RadiusFrom)*0.5,RadiusFrom,RadiusTo);
                //fitter->FixParameter("pp",DLM_Fitter1::p_sor0,CurrentRadius);
            }
            else{
                //fitter->SetParameter("pp",DLM_Fitter1::p_sor0,1.,0.0,2.0);
                fitter->FixParameter("pp",DLM_Fitter1::p_sor0,CurrentAlpha);
                fitter->FixParameter("pp",DLM_Fitter1::p_sor1,0);
                //fitter->SetParameter("pp",DLM_Fitter1::p_sor2,(RadiusTo+RadiusFrom)*0.5,RadiusFrom,RadiusTo);
                fitter->FixParameter("pp",DLM_Fitter1::p_sor2,CurrentRadius);//1.20847
                fitter->FixParameter("pp",DLM_Fitter1::p_sor3,0);
            }
            fitter->GoBabyGo();
            if(fitter->GetChi2Ndf()<BestChi2Ndf){
                BestChi2Ndf=fitter->GetChi2Ndf();
                BestPval=fitter->GetPval();
                BestChi2=fitter->GetChi2();
                if(SourceType==-1)BestRadius=fitter->GetParameter("pp",DLM_Fitter1::p_sor0);
                else if(SourceType==0)BestRadius=fitter->GetParameter("pp",DLM_Fitter1::p_sor0);
                else BestRadius=fitter->GetParameter("pp",DLM_Fitter1::p_sor2);
            }
            hChi2Map->SetBinContent(uFit_r+1,uFit_alpha+1,fitter->GetChi2Ndf());
        }
        //gChi2Ndf.SetPoint(uFit_r,CurrentRadius,fitter->GetChi2Ndf());
    }
    printf("Fitting completed!\n");

    if(SourceType==-1){
        printf("Renorm = %.3f+/-%.3f\n",fitter->GetParameter("pp",DLM_Fitter1::p_sor0),fitter->GetParError("pp",DLM_Fitter1::p_sor0));
    }
    else if(SourceType==0){
        printf("r0 = %.3f+/-%.3f\n",fitter->GetParameter("pp",DLM_Fitter1::p_sor0),fitter->GetParError("pp",DLM_Fitter1::p_sor0));
    }
    else{
        printf("α = %.3f+/-%.3f\n",fitter->GetParameter("pp",DLM_Fitter1::p_sor0),fitter->GetParError("pp",DLM_Fitter1::p_sor0));
        printf("r0 = %.3f+/-%.3f\n",fitter->GetParameter("pp",DLM_Fitter1::p_sor2),fitter->GetParError("pp",DLM_Fitter1::p_sor2));
    }

    printf("BestRadius = %.3f\n",BestRadius);
    printf("BestAlpha = %.3f\n",BestAlpha);
    printf("BestChi2Ndf = %.2f\n",BestChi2Ndf);
    printf("BestChi2 = %.2f\n",BestChi2);
    printf("BestPval = %.2f\n",BestPval);
    printf("BestNsigma = %.3f\n",sqrt(2)*TMath::ErfcInverse(BestPval));

    printf("pp(p_a) = %.3f+/-%.3f\n",fitter->GetParameter("pp",DLM_Fitter1::p_a),fitter->GetParError("pp",DLM_Fitter1::p_a));
    printf("pp(p_b) = %.3f+/-%.3f\n",fitter->GetParameter("pp",DLM_Fitter1::p_b)*1000.,fitter->GetParError("pp",DLM_Fitter1::p_b)*1000.);
    printf("pp(p_sor0) = %.3f+/-%.3f\n",fitter->GetParameter("pp",DLM_Fitter1::p_sor0),fitter->GetParError("pp",DLM_Fitter1::p_sor0));
    if(SourceType==1) printf("pp(p_sor2) = %.3f+/-%.3f\n",fitter->GetParameter("pp",DLM_Fitter1::p_sor2),fitter->GetParError("pp",DLM_Fitter1::p_sor2));
    printf("chi2/ndf = %.3f/%i = %.3f\n",fitter->GetChi2(),fitter->GetNdf(),fitter->GetChi2Ndf());
    printf("pval = %.3f\n",fitter->GetPval());
    printf("nsigma = %.3f\n",sqrt(2)*TMath::ErfcInverse(fitter->GetPval()));



    TGraph FitResult_pp;
    FitResult_pp.SetName(TString::Format("FitResult_pp"));
    fitter->GetFitGraph(0, FitResult_pp);

    DataHisto_pp->Write();
    FitResult_pp.Write();
    gChi2Ndf.Write();
    hChi2Map->Write();
    if(hSource) hSource->Write();
    for(unsigned uMomBin=0; uMomBin<NumMomBins_pp[1]; uMomBin++){
        if(!hkSource) break;
        if(hkSource[uMomBin]) hkSource[uMomBin]->Write();
    }

    if(DEBUG_MODE){printf("Deleting stuff... ");}

    for(short sCat=0; sCat<2; sCat++){
        for(unsigned short usCh=0; usCh<NumChannels_pp; usCh++){
            delete cPars_pp[sCat][usCh];
        }
        delete [] cPars_pp[sCat];
        delete pp_Source[sCat];
    }
    delete [] cPars_pp;
    delete [] pp_Source;
    delete [] AB_pp;

    for(unsigned uVar=0; uVar<3; uVar++){
        delete [] Purities_p[uVar];
        delete [] Fraction_p[uVar];
        delete [] Purities_L[uVar];
        delete [] Fraction_L[uVar];
    }
    delete [] Purities_p;
    delete [] Fraction_p;
    delete [] Purities_L;
    delete [] Fraction_L;
    delete hFeed_pp_pL_MeV;
    delete hReso_pp_MeV;
    delete FileFeed;
    delete FileReso;
    if(DataFile_pp) delete DataFile_pp;
    if(SystErrFile_pp) delete SystErrFile_pp;
    delete Ck_pp;
    delete Ck_pL;

    delete hChi2Map;

    delete OutputFile;
    if(hSource){delete hSource;hSource=NULL;}
    for(unsigned uMomBin=0; uMomBin<NumMomBins_pp[1]; uMomBin++){
        if(!hkSource) break;
        delete hkSource[uMomBin];
    }
    if(hkSource){delete hkSource;hkSource=NULL;}
    delete [] StabDistr;

    if(DEBUG_MODE){printf("Done!\n");}
}


double fitterRad_Gauss(double* rad, double* par){
    double& Radius = rad[0];
    if(Radius<=0) return 0;
    double& Size = par[0];
    if(Size<=0) return 0;
    return 4.*PI*Radius*Radius*pow(4.*PI*Size*Size,-1.5)*exp(-(Radius*Radius)/(4.*Size*Size));
}
double fitterRad_Cauchy(double* rad, double* par){
    double& Radius = rad[0];
    if(Radius<=0) return 0;
    double& Size = par[0];
    return 2.00*2.*Size*Radius*Radius/PI*pow(Radius*Radius+0.25*2.00*2.00*Size*Size,-2.);
}
double fitterRad_GC(double* rad, double* par){
    double& Radius = rad[0];
    if(Radius<=0) return 0;
    double& Size1 = par[0];
    if(Size1<=0) return 0;
    double& Size2 = par[1];
    if(Size2<=0) return 0;
    double& Weight1 = par[2];
    if(Weight1<0 || Weight1>1) return 0;
    return  Weight1 *4.*PI*Radius*Radius*pow(4.*PI*Size1*Size1,-1.5)*exp(-(Radius*Radius)/(4.*Size1*Size1))+
            (1.-Weight1)*2.00*2.*Size2*Radius*Radius/PI*pow(Radius*Radius+0.25*2.00*2.00*Size2*Size2,-2.);
}
double fitterRad_GG(double* rad, double* par){
    double& Radius = rad[0];
    if(Radius<=0) return 0;
    double& Size1 = par[0];
    if(Size1<=0) return 0;
    double& Size2 = par[1];
    if(Size2<=0) return 0;
    double& Weight1 = par[2];
    if(Weight1<0 || Weight1>1) return 0;
    return      Weight1 *4.*PI*Radius*Radius*pow(4.*PI*Size1*Size1,-1.5)*exp(-(Radius*Radius)/(4.*Size1*Size1))+
            (1.-Weight1)*4.*PI*Radius*Radius*pow(4.*PI*Size2*Size2,-1.5)*exp(-(Radius*Radius)/(4.*Size2*Size2));
}
double fitterRad_CC(double* rad, double* par){
    double& Radius = rad[0];
    if(Radius<=0) return 0;
    double& Size1 = par[0];
    if(Size1<=0) return 0;
    double& Size2 = par[1];
    if(Size2<=0) return 0;
    double& Weight1 = par[2];
    if(Weight1<0 || Weight1>1) return 0;
    return      Weight1 *2.00*2.*Size1*Radius*Radius/PI*pow(Radius*Radius+0.25*2.00*2.00*Size1*Size1,-2.)+
            (1.-Weight1)*2.00*2.*Size2*Radius*Radius/PI*pow(Radius*Radius+0.25*2.00*2.00*Size2*Size2,-2.);
}
/*
double fitterRad_CleverLevy(double* rad, double* par){
    double& Radius = rad[0];
    if(Radius<=0) return 0;
    double& Scale = par[0];
    if(Size1<=0) return 0;
    double& Stability = par[1];

}
*/
    //const TString System = "pp";
    //const TString System = "pL";
    //const TString SubSample = "Full";
    //const TString SubSample = "Prim";
void FitEposSource_pp(const TString& System, const TString& SubSample){

    TString OutputFolder;
    if(System=="pp") OutputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/ForGerhard/Test/pp/";
    else if(System=="pL") OutputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/ForGerhard/Test/pL/";

    //there are 2 CATS objects, one for the fitting (4 MeV bin width) and one for the differential look into the source (40 MeV bin width)
    const unsigned MaxNumPairs = 128e6;
    const unsigned MaxNumPairsPerBin = 1e6;
    const short MixingDepth = 32;//!put to one for the actual source. To enhance statistics at the expense of 5-10% uncertainty, put to 8
    const unsigned NumMomBins_pp = 1;
    //const unsigned NumMomBins_pp = 30;
    const unsigned kMin_pp = 0;
    const unsigned kMax_pp = 150;
    const double GaussSourceSize = 1.2;

    const unsigned NumRadBins=128;
    const double RadMin=0;
    const double RadMax=8;

    TString EPOSfile_pp;
    if(System=="pp"){
        //EPOSfile_pp = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/Scratch9_OSCAR1997_AllLrzProtonFiles.f19";
        //EPOSfile_pp = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/7TeV/pp_nist0_32Ki.f19";
        //EPOSfile_pp = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/EPOS_LBF_pp90/pp_nist0_256Ki.f19";
        //EPOSfile_pp = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/EPOS_LBF_pp91/pp_nist0_256Ki.f19";
        //EPOSfile_pp = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/EPOS_LBF_pp91/pp_nist0_prim_256Ki.f19";
        //EPOSfile_pp = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/EPOS_LBF_pp91/pp_nist0_seco_256Ki.f19";
        if(SubSample=="Full") EPOSfile_pp = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/EPOS_LBF_pp91/pp_nist0_030219.f19";
        else if(SubSample=="Prim") EPOSfile_pp = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/EPOS_LBF_pp91/pp_nist0_prim_030219.f19";
    }
    else if(System=="pL"){
        if(SubSample=="Full") EPOSfile_pp = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/EPOS_LBF_pp91/pL_nist0_030219.f19";
        else if(SubSample=="Prim") EPOSfile_pp = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/EPOS_LBF_pp91/pL_nist0_prim_030219.f19";
    }

    TString SourceDiscr;
    if(SubSample=="Prim"){
        SourceDiscr = "prim.";
    }
    else if(SubSample=="Seco"){
        SourceDiscr = "second.";
    }
    else if(SubSample=="Full"){
        SourceDiscr = "full";
    }
    else if(EPOSfile_pp == "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/Scratch9_OSCAR1997_AllLrzProtonFiles.f19"){
        SourceDiscr = "7TeV full";
    }
    else{
        SourceDiscr = "???";
    }

    CATS AB_pp;
    CATSparameters pp_Source(CATSparameters::tSource,1,true);
    pp_Source.SetParameter(0, GaussSourceSize);
    AB_pp.SetMaxNumThreads(1);
    AB_pp.SetUseAnalyticSource(false);
    AB_pp.SetInputFileName(EPOSfile_pp.Data());
    AB_pp.SetMaxPairsToRead(MaxNumPairs);
    AB_pp.SetMaxPairsPerBin(MaxNumPairsPerBin);
    AB_pp.SetMixingDepth(MixingDepth);
    AB_pp.SetThetaDependentSource(false);
    AB_pp.SetTransportRenorm(1.0);
    AB_pp.SetTauCorrection(false);

    const unsigned NumChannels_pp=4;
    AB_pp.SetExcludeFailedBins(false);
    AB_pp.SetMomBins(NumMomBins_pp,kMin_pp,kMax_pp);
    AB_pp.SetNumChannels(NumChannels_pp);
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

    AB_pp.SetQ1Q2(1);
    if(System=="pp") AB_pp.SetPdgId(2212, 2212);
    else if(System=="pL") AB_pp.SetPdgId(2212, 3122);
    AB_pp.SetRedMass( 0.5*Mass_p );

    CATSparameters** cPars_pp = new CATSparameters* [NumChannels_pp];
    for(unsigned short usCh=0; usCh<NumChannels_pp; usCh++){
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

    TFile* fOutput = new TFile(OutputFolder+"fOutput_"+SubSample+".root","recreate");
    TH1F* hRad = new TH1F("hRad","hRad",NumRadBins,RadMin,RadMax);
//AB_pp.SourceHistoTemp.GetNbins();
//AB_pp.SourceHistoTemp.GetBinLowEdge(0,0);
//TH1F* hRadTest = new TH1F("hRadTest","hRadTest",AB_pp.SourceHistoTemp.GetNbins(),AB_pp.SourceHistoTemp.GetBinLowEdge(0,0),
//                          AB_pp.SourceHistoTemp.GetBinUpEdge(0,AB_pp.SourceHistoTemp.GetNbins()-1));
    const unsigned NumPairs = AB_pp.GetNumPairs();
    double RelMom,RelPos,RelCosTh,TotMom;
    for(unsigned uPair=0; uPair<NumPairs; uPair++){
        AB_pp.GetPairInfo(uPair,RelMom,RelPos,RelCosTh,TotMom);
        hRad->Fill(RelPos);
    }
//double NumTestPairs=0;
//for(unsigned uBin=0; uBin<AB_pp.SourceHistoTemp.GetNbins(); uBin++){
//hRadTest->SetBinContent(uBin+1,AB_pp.SourceHistoTemp.GetBinContent(uBin));
//NumTestPairs+=AB_pp.SourceHistoTemp.GetBinContent(uBin);
//}
    hRad->Sumw2();
    hRad->Scale(1./double(NumPairs),"width");
//hRadTest->Sumw2();
//hRadTest->Scale(1./NumTestPairs,"width");

    TF1* fitRad_Gauss = new TF1("fitRad_Gauss",fitterRad_Gauss,RadMin,RadMax,1);
    TF1* fitRad_Cauchy = new TF1("fitRad_Cauchy",fitterRad_Cauchy,RadMin,RadMax,1);
    TF1* fitRad_GC = new TF1("fitRad_GC",fitterRad_GC,RadMin,RadMax,3);
    TF1* fitRad_GG = new TF1("fitRad_GG",fitterRad_GG,RadMin,RadMax,3);
    TF1* fitRad_CC = new TF1("fitRad_CC",fitterRad_CC,RadMin,RadMax,3);

    fitRad_Gauss->SetParameter(0,1);
    fitRad_Gauss->SetParLimits(0,0.1,2);
    fitRad_Gauss->SetLineColor(kRed);
    fitRad_Gauss->SetLineWidth(4);

    fitRad_Cauchy->SetParameter(0,1);
    fitRad_Cauchy->SetParLimits(0,0.1,2);
    fitRad_Cauchy->SetLineColor(kRed);
    fitRad_Cauchy->SetLineWidth(4);

    fitRad_GC->SetParameter(0,1);
    fitRad_GC->SetParLimits(0,0.1,2);
    fitRad_GC->SetParameter(1,1);
    fitRad_GC->SetParLimits(1,0.1,2);
    fitRad_GC->SetParameter(2,0.5);
    fitRad_GC->SetParLimits(2,0,1);
    fitRad_GC->SetLineColor(kRed);
    fitRad_GC->SetLineWidth(4);

    fitRad_GG->SetParameter(0,1);
    fitRad_GG->SetParLimits(0,0.1,2);
    fitRad_GG->SetParameter(1,1);
    fitRad_GG->SetParLimits(1,0.1,2);
    fitRad_GG->SetParameter(2,0.5);
    fitRad_GG->SetParLimits(2,0,1);
    fitRad_GG->SetLineColor(kRed);
    fitRad_GG->SetLineWidth(4);

    fitRad_CC->SetParameter(0,1);
    fitRad_CC->SetParLimits(0,0.1,2);
    fitRad_CC->SetParameter(1,1);
    fitRad_CC->SetParLimits(1,0.1,2);
    fitRad_CC->SetParameter(2,0.5);
    fitRad_CC->SetParLimits(2,0,1);
    fitRad_CC->SetLineColor(kRed);
    fitRad_CC->SetLineWidth(4);

    hRad->Fit(fitRad_Gauss, "S, N, R, M");
    hRad->Fit(fitRad_Cauchy, "S, N, R, M");
    hRad->Fit(fitRad_GC, "S, N, R, M");
    hRad->Fit(fitRad_GG, "S, N, R, M");
    hRad->Fit(fitRad_CC, "S, N, R, M");

    hRad->SetStats(false);
    hRad->SetTitle("");
    hRad->GetXaxis()->SetLabelSize(0.065);
    hRad->GetXaxis()->SetTitle("r (fm)");
    hRad->GetXaxis()->CenterTitle();
    hRad->GetXaxis()->SetTitleOffset(1.25);
    hRad->GetXaxis()->SetLabelOffset(0.02);
    hRad->GetXaxis()->SetTitleSize(0.075);
    hRad->GetYaxis()->SetLabelSize(0.065);
    hRad->GetYaxis()->SetTitle("4#pi r^{2}S(r) (1/fm)");
    hRad->GetYaxis()->CenterTitle();
    hRad->GetYaxis()->SetTitleOffset(0.8);
    hRad->GetYaxis()->SetTitleSize(0.075);
    hRad->GetYaxis()->SetRangeUser(0, 0.65);
    hRad->SetLineColor(kBlack);
    hRad->SetLineWidth(3);

    double Chi2Ndf_Gauss = fitRad_Gauss->GetChisquare()/double(fitRad_Gauss->GetNDF());
    double Chi2Ndf_Cauchy = fitRad_Cauchy->GetChisquare()/double(fitRad_Cauchy->GetNDF());
    double Chi2Ndf_GC = fitRad_GC->GetChisquare()/double(fitRad_GC->GetNDF());
    double Chi2Ndf_GG = fitRad_GG->GetChisquare()/double(fitRad_GG->GetNDF());
    double Chi2Ndf_CC = fitRad_CC->GetChisquare()/double(fitRad_CC->GetNDF());

    printf("Chi2Ndf_Gauss = %.2f\n",Chi2Ndf_Gauss);
    printf("Chi2Ndf_Cauchy = %.2f\n",Chi2Ndf_Cauchy);
    printf("Chi2Ndf_GC = %.2f\n",Chi2Ndf_GC);
    printf("Chi2Ndf_GG = %.2f\n",Chi2Ndf_GG);
    printf("Chi2Ndf_CC = %.2f\n",Chi2Ndf_CC);

    TLegend* myLegend = new TLegend(0.65,0.80,0.95,0.95);//lbrt
    myLegend->SetName("myLegend");
    myLegend->SetTextSize(0.040);
    myLegend->AddEntry(hRad, TString::Format("EPOS source (%s)",SourceDiscr.Data()));
    myLegend->AddEntry(fitRad_Gauss, "Fit");

    TCanvas* cGauss = new TCanvas("cGauss", "cGauss", 1);
    cGauss->cd(0);
    cGauss->SetCanvasSize(1280, 720);
    cGauss->SetMargin(0.15,0.05,0.2,0.05);//lrbt
    hRad->Draw();
    fitRad_Gauss->Draw("same");
    myLegend->Draw("same");

    TPaveText* PT_Gauss = new TPaveText(0.65,0.55,0.95,0.80, "blNDC");//lbrt
    PT_Gauss->SetName("PT_Gauss");
    PT_Gauss->SetBorderSize(1);
    PT_Gauss->SetTextSize(0.045);
    PT_Gauss->SetFillColor(kWhite);
    PT_Gauss->SetTextFont(22);
    PT_Gauss->AddText("Gauss source:");
    PT_Gauss->AddText(TString::Format("#chi^{2}_{ndf}=%.2f\n",Chi2Ndf_Gauss));
    PT_Gauss->AddText(TString::Format("R=%.3f fm\n",fitRad_Gauss->GetParameter(0)));
    PT_Gauss->Draw("same");


    TCanvas* cCauchy = new TCanvas("cCauchy", "cCauchy", 1);
    cCauchy->cd(0);
    cCauchy->SetCanvasSize(1280, 720);
    cCauchy->SetMargin(0.15,0.05,0.2,0.05);//lrbt
    hRad->Draw();
    fitRad_Cauchy->Draw("same");
    myLegend->Draw("same");

    TPaveText* PT_Cauchy = new TPaveText(0.65,0.55,0.95,0.80, "blNDC");//lbrt
    PT_Cauchy->SetName("PT_Cauchy");
    PT_Cauchy->SetBorderSize(1);
    PT_Cauchy->SetTextSize(0.045);
    PT_Cauchy->SetFillColor(kWhite);
    PT_Cauchy->SetTextFont(22);
    PT_Cauchy->AddText("Cauchy source:");
    PT_Cauchy->AddText(TString::Format("#chi^{2}_{ndf}=%.2f\n",Chi2Ndf_Cauchy));
    PT_Cauchy->AddText(TString::Format("R=%.3f fm\n",fitRad_Cauchy->GetParameter(0)));
    PT_Cauchy->Draw("same");

    TCanvas* cGC = new TCanvas("cGC", "cGC", 1);
    cGC->cd(0);
    cGC->SetCanvasSize(1280, 720);
    cGC->SetMargin(0.15,0.05,0.2,0.05);//lrbt
    hRad->Draw();
    fitRad_GC->Draw("same");
    myLegend->Draw("same");

    TPaveText* PT_GC = new TPaveText(0.65,0.50,0.95,0.80, "blNDC");//lbrt
    PT_GC->SetName("PT_GC");
    PT_GC->SetBorderSize(1);
    PT_GC->SetTextSize(0.045);
    PT_GC->SetFillColor(kWhite);
    PT_GC->SetTextFont(22);
    PT_GC->AddText("Gauss+Cauchy source:");
    PT_GC->AddText(TString::Format("#chi^{2}_{ndf}=%.2f\n",Chi2Ndf_GC));
    PT_GC->AddText(TString::Format("RG/RC=%.3f/%.3f fm\n",fitRad_GC->GetParameter(0),fitRad_GC->GetParameter(1)));
    PT_GC->AddText(TString::Format("wG/wC=%.3f/%.3f\n",fitRad_GC->GetParameter(2),1.-fitRad_GC->GetParameter(2)));
    PT_GC->Draw("same");

    TCanvas* cGG = new TCanvas("cGG", "cGG", 1);
    cGG->cd(0);
    cGG->SetCanvasSize(1280, 720);
    cGG->SetMargin(0.15,0.05,0.2,0.05);//lrbt
//cGG->SetLogy();
    hRad->Draw();
    fitRad_GG->Draw("same");
    myLegend->Draw("same");

    TPaveText* PT_GG = new TPaveText(0.65,0.50,0.95,0.80, "blNDC");//lbrt
    PT_GG->SetName("PT_GG");
    PT_GG->SetBorderSize(1);
    PT_GG->SetTextSize(0.045);
    PT_GG->SetFillColor(kWhite);
    PT_GG->SetTextFont(22);
    PT_GG->AddText("Gauss+Gauss source:");
    PT_GG->AddText(TString::Format("#chi^{2}_{ndf}=%.2f\n",Chi2Ndf_GG));
    PT_GG->AddText(TString::Format("RG/RG=%.3f/%.3f fm\n",fitRad_GG->GetParameter(0),fitRad_GG->GetParameter(1)));
    PT_GG->AddText(TString::Format("wG/wG=%.3f/%.3f\n",fitRad_GG->GetParameter(2),1.-fitRad_GG->GetParameter(2)));
    PT_GG->Draw("same");

    TCanvas* cCC = new TCanvas("cCC", "cCC", 1);
    cCC->cd(0);
    cCC->SetCanvasSize(1280, 720);
    cCC->SetMargin(0.15,0.05,0.2,0.05);//lrbt
    hRad->Draw();
    fitRad_CC->Draw("same");
    myLegend->Draw("same");

    TPaveText* PT_CC = new TPaveText(0.65,0.50,0.95,0.80, "blNDC");//lbrt
    PT_CC->SetName("PT_CC");
    PT_CC->SetBorderSize(1);
    PT_CC->SetTextSize(0.045);
    PT_CC->SetFillColor(kWhite);
    PT_CC->SetTextFont(22);
    PT_CC->AddText("Cauchy+Cauchy source:");
    PT_CC->AddText(TString::Format("#chi^{2}_{ndf}=%.2f\n",Chi2Ndf_CC));
    PT_CC->AddText(TString::Format("RC/RC=%.3f/%.3f fm\n",fitRad_CC->GetParameter(0),fitRad_CC->GetParameter(1)));
    PT_CC->AddText(TString::Format("wC/wC=%.3f/%.3f\n",fitRad_CC->GetParameter(2),1.-fitRad_CC->GetParameter(2)));
    PT_CC->Draw("same");

    hRad->Write();
//hRadTest->Write();
    fitRad_Gauss->Write();
    fitRad_Cauchy->Write();
    fitRad_GC->Write();
    fitRad_GG->Write();
    fitRad_CC->Write();
    cGauss->Write();
    cCauchy->Write();
    cGC->Write();
    cGG->Write();
    cCC->Write();

    cGauss->SaveAs(OutputFolder+"cGauss_"+SubSample+".png");
    cCauchy->SaveAs(OutputFolder+"cCauchy_"+SubSample+".png");
    cGC->SaveAs(OutputFolder+"cGC_"+SubSample+".png");
    cGG->SaveAs(OutputFolder+"cGG_"+SubSample+".png");
    cCC->SaveAs(OutputFolder+"cCC_"+SubSample+".png");

    delete hRad;
//delete hRadTest;
    delete fitRad_Gauss;
    delete fitRad_Cauchy;
    delete fitRad_GC;
    delete fitRad_GG;
    delete fitRad_CC;
    delete PT_Gauss;
    delete PT_Cauchy;
    delete PT_GC;
    delete PT_GG;
    delete PT_CC;
    delete myLegend;
    delete cGauss;
    delete cCauchy;
    delete cGC;
    delete cGG;
    delete cCC;

    delete fOutput;
    for(unsigned short usCh=0; usCh<NumChannels_pp; usCh++){delete cPars_pp[usCh];}
    delete [] cPars_pp;

}

void GetSource(const TString& System, const TString& SubSample, const bool& TauCorrection, const bool& BoostToCm,
               const double& kMin, const double& kMax){
    TString OutputFolder;
    if(System=="pp") OutputFolder = TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/ForGerhard/GetSource120319/%.0f_%.0f/pp/",kMin,kMax);
    else if(System=="pL") OutputFolder = TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/ForGerhard/GetSource120319/%.0f_%.0f/pL/",kMin,kMax);
    else if(System=="LL") OutputFolder = TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/ForGerhard/GetSource120319/%.0f_%.0f/LL/",kMin,kMax);
    else if(System=="pXim") OutputFolder = TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/ForGerhard/GetSource120319/%.0f_%.0f/pXim/",kMin,kMax);

    //there are 2 CATS objects, one for the fitting (4 MeV bin width) and one for the differential look into the source (40 MeV bin width)
    const unsigned MaxNumEvents = 250000;
    const short MixingDepth = 32;//!put to one for the actual source. To enhance statistics at the expense of 5-10% uncertainty, put to 8
    //const unsigned NumMomBins_pp = 30;

    const unsigned NumRadBins=128;
    const double RadMin=0;
    const double RadMax=8;

    const unsigned NumMomBins=1024;

    const double Rescale = 1;

    int pdgID[2];

    TString EPOSfile;
    if(System=="pp"){
        if(SubSample=="Full") EPOSfile = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/EPOS_LBF_pp91/pp_nist0_030219.f19";
        else if(SubSample=="Prim") EPOSfile = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/EPOS_LBF_pp91/pp_nist0_prim_030219.f19";
        pdgID[0] = 2212;
        pdgID[1] = 2212;
    }
    else if(System=="pL"){
        if(SubSample=="Full") EPOSfile = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/EPOS_LBF_pp91/pL_nist0_030219.f19";
        else if(SubSample=="Prim") EPOSfile = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/EPOS_LBF_pp91/pL_nist0_prim_030219.f19";
        pdgID[0] = 2212;
        pdgID[1] = 3122;
    }
    else if(System=="LL"){
        if(SubSample=="Full") EPOSfile = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/EPOS_LBF_pp91/.f19";
        else if(SubSample=="Prim") EPOSfile = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/EPOS_LBF_pp91/.f19";
        pdgID[0] = 2212;
        pdgID[1] = 3122;
    }
    else if(System=="pXim"){
        if(SubSample=="Full") EPOSfile = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/EPOS_LBF_pp91/.f19";
        else if(SubSample=="Prim") EPOSfile = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/EPOS_LBF_pp91/.f19";
        pdgID[0] = 2212;
        pdgID[1] = 3122;
    }

    TString SourceDiscr;
    if(SubSample=="Prim"){
        SourceDiscr = "prim.";
    }
    else if(SubSample=="Seco"){
        SourceDiscr = "second.";
    }
    else if(SubSample=="Full"){
        SourceDiscr = "full";
    }
    else if(EPOSfile == "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/Scratch9_OSCAR1997_AllLrzProtonFiles.f19"){
        SourceDiscr = "7TeV full";
    }
    else{
        SourceDiscr = "???";
    }

    DLM_Timer dlmTimer;
    double Time;

    bool ProgressBar=false;
    char* cdummy = new char [512];

    FILE *InFile;
    InFile = fopen(EPOSfile.Data(), "r");
    if(!InFile){
        printf("          \033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", EPOSfile.Data());
        return;
    }

    fseek ( InFile , 0 , SEEK_END );
    long EndPos;
    EndPos = ftell (InFile);
    fseek ( InFile , 0 , SEEK_SET );
    long CurPos;

    const unsigned short NumBlankHeaderLines=3;
    //Read the header lines
    for(unsigned short us=0; us<NumBlankHeaderLines; us++){
        if(!fgets(cdummy, 255, InFile)){
            printf("\033[1;33mWARNING!\033[0m Possible bad input-file, error when reading from %s!\n",EPOSfile.Data());
            continue;
        }
    }

    if(feof(InFile)){
        printf("\033[1;31m          ERROR:\033[0m Trying to read past end of file %s\n", EPOSfile.Data());
        printf("         No particle pairs were loaded :(\n");
        return;
    }

    unsigned NumPairs=0;
    unsigned NumSePairs=0;
    unsigned NumTotalPairs=0;
    unsigned NumTotalSePairs=0;
    unsigned NumSingleParticles1=0;
    unsigned NumSingleParticles2=0;

    int EventNumber;
    int NumPartInEvent;
    double ImpPar;
    double fDummy;

    unsigned TotalNumSePairs=0;
    unsigned NumEvents=0;
    unsigned TotNumEvents = 0;
    //unsigned* NumEvPart = new unsigned [NumIpBins];

    CatsParticle KittyParticle;
    CatsDataBuffer* KittyBuffer = new CatsDataBuffer(MixingDepth,pdgID[0],pdgID[1]);
    unsigned uBuffer=0;
    //CatsDataBuffer KittyBuffer(MixingDepth,pdgID[0],pdgID[1]);

    CatsEvent DummyEvent(pdgID[0],pdgID[1]);
    CatsEvent** KittyEvent = new CatsEvent* [MixingDepth];
    for(unsigned uDepth=0; uDepth<MixingDepth; uDepth++){
        KittyEvent[uDepth] = new CatsEvent(pdgID[0],pdgID[1]);
    }

    //progress
    //percentage of the file read. The reading speed should be more or less constant,
    //so this value can always give an accurate maximum ETA
    float pFile;
    float pEvent;
    float ProgressLoad;

    short pTotal;
    short pTotalOld=0;

    bool NewInterestingEvent;

    unsigned SelectedSePairs;

    unsigned RejectedHighMultEvents = 0;
    const int HighMultLimit = 128;

    TFile* fOutput = new TFile(OutputFolder+"fOutput_"+SubSample+"_"+TauCorrection+"_"+BoostToCm+".root","recreate");
    TH1F** hRad = new TH1F* [2];
    hRad[0] = new TH1F("hRadSe","hRadSe",NumRadBins,RadMin,RadMax);
    hRad[1] = new TH1F("hRadMe","hRadMe",NumRadBins,RadMin,RadMax);

    TH1F** hSingleParticle = new TH1F* [2];
    hSingleParticle[0] = new TH1F("hSingleParticle0","hSingleParticle0",NumRadBins,RadMin,RadMax);
    hSingleParticle[1] = new TH1F("hSingleParticle1","hSingleParticle1",NumRadBins,RadMin,RadMax);

    TH1F** hMomDist = new TH1F* [2];
    hMomDist[0] = new TH1F("hMomDistSe","hMomDistSe",NumMomBins,kMin,kMax);
    hMomDist[1] = new TH1F("hMomDistMe","hMomDistMe",NumMomBins,kMin,kMax);

    //!---Iteration over all events---
    while(!feof(InFile)){

        if(!fscanf(InFile,"%i %i %lf %lf",&EventNumber,&NumPartInEvent,&ImpPar,&fDummy)){
            printf("\033[1;33mWARNING!\033[0m Possible bad input-file, error when reading from %s!\n",EPOSfile.Data());
            continue;
        }
        NewInterestingEvent = false;
        ImpPar = fabs(ImpPar);
        TotNumEvents++;
//printf("TotNumEvents=%u\n",TotNumEvents);
        if(TotNumEvents>=MaxNumEvents)  break;
        if(NumPartInEvent>HighMultLimit) RejectedHighMultEvents++;

        //!---Iteration over all particles in this event---
        for(int iPart=0; iPart<NumPartInEvent; iPart++){
            KittyParticle.ReadFromOscarFile(InFile);
            if(Rescale!=1){
                KittyParticle.RenormSpacialCoordinates(Rescale);
            }

            if(NumPartInEvent>HighMultLimit) continue;

            //sometimes one might go beyond the limit of the file
            if(ftell(InFile)>=EndPos)continue;

            if(KittyParticle.GetE()==0){
                printf("\033[1;33mWARNING!\033[0m Possible bad input-file, there are particles with zero energy!\n");
                continue;
            }

            if(KittyParticle.GetPid()!=pdgID[0] && KittyParticle.GetPid()!=pdgID[1])
                continue; //don't save this particle if it is of the wrong type

            if(!NewInterestingEvent){
                NewInterestingEvent = true;
                NumEvents++;
            }
            KittyEvent[uBuffer]->AddParticle(KittyParticle);
        }//for(int iPart=0; iPart<NumPartInEvent; iPart++)

        for(unsigned uPart=0; uPart<KittyEvent[uBuffer]->GetNumParticles1(); uPart++){
            hSingleParticle[0]->Fill(KittyEvent[uBuffer]->GetParticleType1(uPart).GetR());
//if(NumEvents){
//printf(" %u R=%f\n",uPart,KittyEvent[uBuffer]->GetParticleType1(uPart).GetR());
//}
            NumSingleParticles1++;
        }
        for(unsigned uPart=0; uPart<KittyEvent[uBuffer]->GetNumParticles2(); uPart++){
            hSingleParticle[1]->Fill(KittyEvent[uBuffer]->GetParticleType2(uPart).GetR());
            NumSingleParticles2++;
        }

        KittyEvent[uBuffer]->ComputeParticlePairs(TauCorrection,BoostToCm);
//KittyEvent[uBuffer]->ComputeParticlePairs(false,false);
        KittyBuffer->SetEvent(uBuffer, *KittyEvent[uBuffer]);

        uBuffer++;

        double RelPos;
        double RelCosTh;
        double RelMom;
        double RedMomMeV;

        //if the buffer is full -> empty it!
        //note that if it happens the we leave the while loop before emptying the buffer,
        //uBuffer will be != than zero! use this condition to empty the buffer when exiting the loop!
        if(uBuffer==MixingDepth){
            KittyBuffer->GoBabyGo(TauCorrection,BoostToCm);
            NumTotalPairs+=KittyBuffer->GetNumPairs();
            NumTotalSePairs+=KittyBuffer->GetNumPairsSameEvent();
            uBuffer=0;
            const CatsParticlePair* PairDif;
            for(unsigned uPair=0; uPair<KittyBuffer->GetNumPairs(); uPair++){
                PairDif = KittyBuffer->GetPair(uPair);
                RelPos = PairDif->GetR();
                RelCosTh = (PairDif->GetPx()*PairDif->GetX()+
                               PairDif->GetPy()*PairDif->GetY()+
                               PairDif->GetPz()*PairDif->GetZ())/
                               (PairDif->GetP()*PairDif->GetR());
                RelMom = PairDif->GetP();
                RedMomMeV = 500.*RelMom;
                if(RedMomMeV>=kMin&&RedMomMeV<=kMax){
                    hRad[1]->Fill(RelPos);
                    hMomDist[1]->Fill(RedMomMeV);
                    NumPairs++;
                }
            }
            for(unsigned uPair=0; uPair<KittyBuffer->GetNumPairsSameEvent(); uPair++){
                PairDif = KittyBuffer->GetSePair(uPair);
                RelPos = PairDif->GetR();
                RelCosTh = (PairDif->GetPx()*PairDif->GetX()+
                               PairDif->GetPy()*PairDif->GetY()+
                               PairDif->GetPz()*PairDif->GetZ())/
                               (PairDif->GetP()*PairDif->GetR());
                RelMom = PairDif->GetP();
                RedMomMeV = 500.*RelMom;
                if(RedMomMeV>=kMin&&RedMomMeV<=kMax){
                    hRad[0]->Fill(RelPos);
                    hMomDist[0]->Fill(RedMomMeV);
//if(RelPos>RadMax){
//printf("---> RelPos=%e\n",RelPos);
//}
                    NumSePairs++;
                }
            }
            for(unsigned uDepth=0; uDepth<MixingDepth; uDepth++){
                KittyEvent[uDepth]->Reset();
            }
        }

        CurPos = ftell (InFile);
        pEvent = double(NumEvents)/double(MaxNumEvents);//
//        pMaxPairsToLoad = double(NumPairs)/double(MaxTotPairs);
        pFile = double(CurPos)/double(EndPos);//what fraction of the file has been read
        ProgressLoad = pEvent;
        if(pFile>ProgressLoad) ProgressLoad=pFile;

        pTotal = int(ProgressLoad*100);
        if(pTotal!=pTotalOld){
            Time = double(dlmTimer.Stop())/1e6;
            //EtaPerBin = round((1./pMaxPairsPerBin-1.)*Time);
            //EtaToLoad = round((1./pMaxPairsToLoad-1.)*Time);
            Time = round((1./ProgressLoad-1.)*Time);
            ShowTime((long long)(Time), cdummy, 2, true, 5);
            printf("\r\033[K          Progress %3d%%, ETA %s",pTotal,cdummy);
            cout << flush;
            pTotalOld = pTotal;
        }
    }//while(!feof(InFile))


    printf("TotNumEvents = %u\n",TotNumEvents);
//printf("NumSingleParticles1=%u\n",NumSingleParticles1);
//printf("NumSingleParticles2=%u\n",NumSingleParticles2);


    printf("\r\033[K");

    if(RejectedHighMultEvents){
        printf("\033[1;33m          WARNING:\033[0m Cannot handle very high multiplicity events!\n"
               "                   As a result %u events have been rejected!\n", RejectedHighMultEvents);
        if(double(RejectedHighMultEvents)/double(TotNumEvents)>0.01)
            printf("                   We are sorry for the inconvenience :'(\n"
                   "                   If you get in touch with us we will try to find a fix for you!\n");
    }

    if(!NumPairs){
        printf("\033[1;31m          WARNING:\033[0m There were no pairs loaded! The computation cannot proceed!\n");
    }
    else if(!NumSePairs){
        printf("\033[1;31m          WARNING:\033[0m There were no same-events pairs found! The computation cannot proceed!\n");
    }
    fclose(InFile);

    for(unsigned uDepth=0; uDepth<MixingDepth; uDepth++){
        delete KittyEvent[uDepth];
    }
    delete KittyBuffer;
    delete [] KittyEvent;
    delete [] cdummy;
//return;

    hRad[0]->Sumw2();
    hRad[0]->Scale(1./double(NumSePairs),"width");
    //hRad[0]->Scale(1./hRad[0]->Integral(1,hRad[0]->GetNbinsX()),"width");

    hRad[1]->Sumw2();
    hRad[1]->Scale(1./double(NumPairs),"width");
    //hRad[1]->Scale(1./hRad[1]->Integral(1,hRad[1]->GetNbinsX()),"width");

    hMomDist[0]->Sumw2();
    hMomDist[0]->Scale(1./double(NumSePairs),"width");

    hMomDist[1]->Sumw2();
    hMomDist[1]->Scale(1./double(NumPairs),"width");

    printf("Number of same events pairs: %u\n",NumSePairs);
    printf("Number of mixed events pairs: %u\n",NumPairs);

    //cheap way of getting rid of the Sigma0
    if(System=="pL"&&SubSample=="Full"){
        hRad[0]->Scale(4./3.);
        hRad[1]->Scale(4./3.);
    }

printf("hRad[1]->Integral(1,hRad[1]->GetNbinsX())=%f (%f)\n",hRad[1]->Integral(1,hRad[1]->GetNbinsX()),hRad[1]->Integral(1,hRad[1]->GetNbinsX()+1));

    hSingleParticle[0]->Sumw2();
    hSingleParticle[0]->Scale(1./double(NumSingleParticles1),"width");

    hSingleParticle[1]->Sumw2();
    hSingleParticle[1]->Scale(1./double(NumSingleParticles2),"width");

    TF1** fitRad_Gauss = new TF1* [2];
    TF1** fitRad_Cauchy = new TF1* [2];
    TF1** fitRad_GC = new TF1* [2];
    TF1** fitRad_GG = new TF1* [2];
    TF1** fitRad_CC = new TF1* [2];
    TF1** fitRad_CleverLevy = new TF1* [2];

    double Chi2Ndf_Gauss[2];
    double Chi2Ndf_Cauchy[2];
    double Chi2Ndf_GC[2];
    double Chi2Ndf_GG[2];
    double Chi2Ndf_CC[2];
    double Chi2Ndf_CleverLevy[2];

    TLegend** myLegend = new TLegend*[2];
    TCanvas** cGauss = new TCanvas*[2];
    TPaveText** PT_Gauss = new TPaveText*[2];

    TCanvas** cCauchy = new TCanvas*[2];
    TPaveText** PT_Cauchy = new TPaveText*[2];

    TCanvas** cGC = new TCanvas*[2];
    TPaveText** PT_GC = new TPaveText*[2];

    TCanvas** cGG = new TCanvas*[2];
    TPaveText** PT_GG = new TPaveText*[2];

    TCanvas** cCC = new TCanvas*[2];
    TPaveText** PT_CC = new TPaveText*[2];

    TCanvas** cCleverLevy = new TCanvas*[2];
    TPaveText** PT_CleverLevy = new TPaveText*[2];

    for(int iSeMe=0; iSeMe<2; iSeMe++){
        TString AddOn;
        TString SystemName;

        if(System=="pp"){SystemName="pp";}
        else if(System=="pL"){SystemName="p#Lambda";}
        else if(System=="LL"){SystemName="#Lambda#Lambda";}
        else if(System=="pXim"){SystemName="p#Xi^{#minus}";}

        AddOn = iSeMe==0?"Se":"Me";
        if(TauCorrection) AddOn += "_Tau";
        else AddOn += "_NoTau";

        if(BoostToCm) AddOn += "_Cm";
        else AddOn += "_Lab";

        fitRad_Gauss[iSeMe] = new TF1("fitRad_Gauss_"+AddOn,fitterRad_Gauss,RadMin,RadMax,1);
        fitRad_Cauchy[iSeMe] = new TF1("fitRad_Cauchy_"+AddOn,fitterRad_Cauchy,RadMin,RadMax,1);
        fitRad_GC[iSeMe] = new TF1("fitRad_GC_"+AddOn,fitterRad_GC,RadMin,RadMax,3);
        fitRad_GG[iSeMe] = new TF1("fitRad_GG_"+AddOn,fitterRad_GG,RadMin,RadMax,3);
        fitRad_CC[iSeMe] = new TF1("fitRad_CC_"+AddOn,fitterRad_CC,RadMin,RadMax,3);

        DLM_CleverLevy CleverLevy;
        CleverLevy.InitType(1);
        fitRad_CleverLevy[iSeMe] = new TF1("fitRad_CleverLevy"+AddOn,&CleverLevy,&DLM_CleverLevy::RootEval,RadMin,RadMax,2,"DLM_CleverLevy","RootEval");
        //fitRad_CleverLevy[iSeMe] = new TF1("fitRad_CleverLevy"+AddOn,fitterRad_CleverLevy,RadMin,RadMax,2);

        fitRad_Gauss[iSeMe]->SetParameter(0,1);
        fitRad_Gauss[iSeMe]->SetParLimits(0,0.1,2);
        fitRad_Gauss[iSeMe]->SetLineColor(kRed);
        fitRad_Gauss[iSeMe]->SetLineWidth(4);

        fitRad_Cauchy[iSeMe]->SetParameter(0,1);
        fitRad_Cauchy[iSeMe]->SetParLimits(0,0.1,2);
        fitRad_Cauchy[iSeMe]->SetLineColor(kRed);
        fitRad_Cauchy[iSeMe]->SetLineWidth(4);

        fitRad_GC[iSeMe]->SetParameter(0,1);
        fitRad_GC[iSeMe]->SetParLimits(0,0.1,2);
        fitRad_GC[iSeMe]->SetParameter(1,1);
        fitRad_GC[iSeMe]->SetParLimits(1,0.1,2);
        fitRad_GC[iSeMe]->SetParameter(2,0.5);
        fitRad_GC[iSeMe]->SetParLimits(2,0,1);
        fitRad_GC[iSeMe]->SetLineColor(kRed);
        fitRad_GC[iSeMe]->SetLineWidth(4);

        fitRad_GG[iSeMe]->SetParameter(0,1);
        fitRad_GG[iSeMe]->SetParLimits(0,0.1,2);
        fitRad_GG[iSeMe]->SetParameter(1,1);
        fitRad_GG[iSeMe]->SetParLimits(1,0.1,2);
        fitRad_GG[iSeMe]->SetParameter(2,0.5);
        fitRad_GG[iSeMe]->SetParLimits(2,0,1);
        fitRad_GG[iSeMe]->SetLineColor(kRed);
        fitRad_GG[iSeMe]->SetLineWidth(4);

        fitRad_CC[iSeMe]->SetParameter(0,1);
        fitRad_CC[iSeMe]->SetParLimits(0,0.1,2);
        fitRad_CC[iSeMe]->SetParameter(1,1);
        fitRad_CC[iSeMe]->SetParLimits(1,0.1,2);
        fitRad_CC[iSeMe]->SetParameter(2,0.5);
        fitRad_CC[iSeMe]->SetParLimits(2,0,1);
        fitRad_CC[iSeMe]->SetLineColor(kRed);
        fitRad_CC[iSeMe]->SetLineWidth(4);

        fitRad_CleverLevy[iSeMe]->SetParameter(0,1);
        fitRad_CleverLevy[iSeMe]->SetParLimits(0,0.1,2);
        fitRad_CleverLevy[iSeMe]->SetParameter(1,1.5);
        fitRad_CleverLevy[iSeMe]->SetParLimits(1,1,2);
        fitRad_CleverLevy[iSeMe]->SetLineColor(kRed);
        fitRad_CleverLevy[iSeMe]->SetLineWidth(4);

        hRad[iSeMe]->Fit(fitRad_Gauss[iSeMe], "Q, S, N, R, M");
        hRad[iSeMe]->Fit(fitRad_Cauchy[iSeMe], "Q, S, N, R, M");
        hRad[iSeMe]->Fit(fitRad_GC[iSeMe], "Q, S, N, R, M");
        hRad[iSeMe]->Fit(fitRad_GG[iSeMe], "Q, S, N, R, M");
        hRad[iSeMe]->Fit(fitRad_CC[iSeMe], "Q, S, N, R, M");
        hRad[iSeMe]->Fit(fitRad_CleverLevy[iSeMe], "Q, S, N, R, M");

        hRad[iSeMe]->SetStats(false);
        hRad[iSeMe]->SetTitle("");
        hRad[iSeMe]->GetXaxis()->SetLabelSize(0.065);
        hRad[iSeMe]->GetXaxis()->SetTitle("r (fm)");
        hRad[iSeMe]->GetXaxis()->CenterTitle();
        hRad[iSeMe]->GetXaxis()->SetTitleOffset(1.25);
        hRad[iSeMe]->GetXaxis()->SetLabelOffset(0.02);
        hRad[iSeMe]->GetXaxis()->SetTitleSize(0.075);
        hRad[iSeMe]->GetYaxis()->SetLabelSize(0.065);
        hRad[iSeMe]->GetYaxis()->SetTitle("4#pi r^{2}S(r) (1/fm)");
        hRad[iSeMe]->GetYaxis()->CenterTitle();
        hRad[iSeMe]->GetYaxis()->SetTitleOffset(0.8);
        hRad[iSeMe]->GetYaxis()->SetTitleSize(0.075);
        hRad[iSeMe]->GetYaxis()->SetRangeUser(0, 0.65);
        hRad[iSeMe]->SetLineColor(kBlack);
        hRad[iSeMe]->SetLineWidth(3);

        Chi2Ndf_Gauss[iSeMe] = fitRad_Gauss[iSeMe]->GetChisquare()/double(fitRad_Gauss[iSeMe]->GetNDF());
        Chi2Ndf_Cauchy[iSeMe] = fitRad_Cauchy[iSeMe]->GetChisquare()/double(fitRad_Cauchy[iSeMe]->GetNDF());
        Chi2Ndf_GC[iSeMe] = fitRad_GC[iSeMe]->GetChisquare()/double(fitRad_GC[iSeMe]->GetNDF());
        Chi2Ndf_GG[iSeMe] = fitRad_GG[iSeMe]->GetChisquare()/double(fitRad_GG[iSeMe]->GetNDF());
        Chi2Ndf_CC[iSeMe] = fitRad_CC[iSeMe]->GetChisquare()/double(fitRad_CC[iSeMe]->GetNDF());
        Chi2Ndf_CleverLevy[iSeMe] = fitRad_CleverLevy[iSeMe]->GetChisquare()/double(fitRad_CleverLevy[iSeMe]->GetNDF());

        printf("%s:\n",AddOn.Data());
        printf("Chi2Ndf_Gauss = %.2f\n",Chi2Ndf_Gauss[iSeMe]);
        printf("Chi2Ndf_Cauchy = %.2f\n",Chi2Ndf_Cauchy[iSeMe]);
        printf("Chi2Ndf_GC = %.2f\n",Chi2Ndf_GC[iSeMe]);
        printf("Chi2Ndf_GG = %.2f\n",Chi2Ndf_GG[iSeMe]);
        printf("Chi2Ndf_CC = %.2f\n",Chi2Ndf_CC[iSeMe]);
        printf("Chi2Ndf_CleverLevy = %.2f\n",Chi2Ndf_CleverLevy[iSeMe]);
        printf("\n");

        myLegend[iSeMe] = new TLegend(0.65,0.80,0.95,0.95);//lbrt
        myLegend[iSeMe]->SetName("myLegend"+AddOn);
        myLegend[iSeMe]->SetTextSize(0.040);
        myLegend[iSeMe]->AddEntry(hRad[iSeMe], TString::Format("EPOS source (%s)",SourceDiscr.Data()));
        myLegend[iSeMe]->AddEntry(fitRad_Gauss[iSeMe], "Fit");

        cGauss[iSeMe] = new TCanvas("cGauss_"+AddOn, "cGauss_"+AddOn, 1);
        cGauss[iSeMe]->cd(0);
        cGauss[iSeMe]->SetCanvasSize(1280, 720);
        cGauss[iSeMe]->SetMargin(0.15,0.05,0.2,0.05);//lrbt
        hRad[iSeMe]->Draw();
        fitRad_Gauss[iSeMe]->Draw("same");
        myLegend[iSeMe]->Draw("same");

        PT_Gauss[iSeMe] = new TPaveText(0.65,0.55,0.95,0.80, "blNDC");//lbrt
        PT_Gauss[iSeMe]->SetName("PT_Gauss"+AddOn);
        PT_Gauss[iSeMe]->SetBorderSize(1);
        PT_Gauss[iSeMe]->SetTextSize(0.045);
        PT_Gauss[iSeMe]->SetFillColor(kWhite);
        PT_Gauss[iSeMe]->SetTextFont(22);
        PT_Gauss[iSeMe]->AddText("Gauss source ("+SystemName+"):");
        PT_Gauss[iSeMe]->AddText(TString::Format("#chi^{2}_{ndf}=%.2f\n",Chi2Ndf_Gauss[iSeMe]));
        PT_Gauss[iSeMe]->AddText(TString::Format("R=%.3f fm\n",fitRad_Gauss[iSeMe]->GetParameter(0)));
        PT_Gauss[iSeMe]->Draw("same");


        cCauchy[iSeMe] = new TCanvas("cCauchy_"+AddOn, "cCauchy_"+AddOn, 1);
        cCauchy[iSeMe]->cd(0);
        cCauchy[iSeMe]->SetCanvasSize(1280, 720);
        cCauchy[iSeMe]->SetMargin(0.15,0.05,0.2,0.05);//lrbt
        hRad[iSeMe]->Draw();
        fitRad_Cauchy[iSeMe]->Draw("same");
        myLegend[iSeMe]->Draw("same");

        PT_Cauchy[iSeMe] = new TPaveText(0.65,0.55,0.95,0.80, "blNDC");//lbrt
        PT_Cauchy[iSeMe]->SetName("PT_Cauchy"+AddOn);
        PT_Cauchy[iSeMe]->SetBorderSize(1);
        PT_Cauchy[iSeMe]->SetTextSize(0.045);
        PT_Cauchy[iSeMe]->SetFillColor(kWhite);
        PT_Cauchy[iSeMe]->SetTextFont(22);
        PT_Cauchy[iSeMe]->AddText("Cauchy source ("+SystemName+"):");
        PT_Cauchy[iSeMe]->AddText(TString::Format("#chi^{2}_{ndf}=%.2f\n",Chi2Ndf_Cauchy[iSeMe]));
        PT_Cauchy[iSeMe]->AddText(TString::Format("R=%.3f fm\n",fitRad_Cauchy[iSeMe]->GetParameter(0)));
        PT_Cauchy[iSeMe]->Draw("same");


        cGC[iSeMe] = new TCanvas("cGC_"+AddOn, "cGC_"+AddOn, 1);
        cGC[iSeMe]->cd(0);
        cGC[iSeMe]->SetCanvasSize(1280, 720);
        cGC[iSeMe]->SetMargin(0.15,0.05,0.2,0.05);//lrbt
        hRad[iSeMe]->Draw();
        fitRad_GC[iSeMe]->Draw("same");
        myLegend[iSeMe]->Draw("same");

        PT_GC[iSeMe] = new TPaveText(0.65,0.50,0.95,0.80, "blNDC");//lbrt
        PT_GC[iSeMe]->SetName("PT_GC"+AddOn);
        PT_GC[iSeMe]->SetBorderSize(1);
        PT_GC[iSeMe]->SetTextSize(0.045);
        PT_GC[iSeMe]->SetFillColor(kWhite);
        PT_GC[iSeMe]->SetTextFont(22);
        PT_GC[iSeMe]->AddText("GC source ("+SystemName+"):");
        PT_GC[iSeMe]->AddText(TString::Format("#chi^{2}_{ndf}=%.2f\n",Chi2Ndf_GC[iSeMe]));
        PT_GC[iSeMe]->AddText(TString::Format("RG/RC=%.3f/%.3f fm\n",fitRad_GC[iSeMe]->GetParameter(0),fitRad_GC[iSeMe]->GetParameter(1)));
        PT_GC[iSeMe]->AddText(TString::Format("wG/wC=%.3f/%.3f fm\n",fitRad_GC[iSeMe]->GetParameter(2),1.-fitRad_GC[iSeMe]->GetParameter(2)));
        PT_GC[iSeMe]->Draw("same");

        cGG[iSeMe] = new TCanvas("cGG_"+AddOn, "cGG_"+AddOn, 1);
        cGG[iSeMe]->cd(0);
        cGG[iSeMe]->SetCanvasSize(1280, 720);
        cGG[iSeMe]->SetMargin(0.15,0.05,0.2,0.05);//lrbt
        hRad[iSeMe]->Draw();
        fitRad_GG[iSeMe]->Draw("same");
        myLegend[iSeMe]->Draw("same");

        PT_GG[iSeMe] = new TPaveText(0.65,0.50,0.95,0.80, "blNDC");//lbrt
        PT_GG[iSeMe]->SetName("PT_GG"+AddOn);
        PT_GG[iSeMe]->SetBorderSize(1);
        PT_GG[iSeMe]->SetTextSize(0.045);
        PT_GG[iSeMe]->SetFillColor(kWhite);
        PT_GG[iSeMe]->SetTextFont(22);
        PT_GG[iSeMe]->AddText("GG source ("+SystemName+"):");
        PT_GG[iSeMe]->AddText(TString::Format("#chi^{2}_{ndf}=%.2f\n",Chi2Ndf_GG[iSeMe]));
        PT_GG[iSeMe]->AddText(TString::Format("RG1/RG2=%.3f/%.3f fm\n",fitRad_GG[iSeMe]->GetParameter(0),fitRad_GG[iSeMe]->GetParameter(1)));
        PT_GG[iSeMe]->AddText(TString::Format("wG1/wG2=%.3f/%.3f fm\n",fitRad_GG[iSeMe]->GetParameter(2),1.-fitRad_GG[iSeMe]->GetParameter(2)));
        PT_GG[iSeMe]->Draw("same");

        cCC[iSeMe] = new TCanvas("cCC_"+AddOn, "cCC_"+AddOn, 1);
        cCC[iSeMe]->cd(0);
        cCC[iSeMe]->SetCanvasSize(1280, 720);
        cCC[iSeMe]->SetMargin(0.15,0.05,0.2,0.05);//lrbt
        hRad[iSeMe]->Draw();
        fitRad_CC[iSeMe]->Draw("same");
        myLegend[iSeMe]->Draw("same");

        PT_CC[iSeMe] = new TPaveText(0.65,0.55,0.95,0.80, "blNDC");//lbrt
        PT_CC[iSeMe]->SetName("PT_CC"+AddOn);
        PT_CC[iSeMe]->SetBorderSize(1);
        PT_CC[iSeMe]->SetTextSize(0.045);
        PT_CC[iSeMe]->SetFillColor(kWhite);
        PT_CC[iSeMe]->SetTextFont(22);
        PT_CC[iSeMe]->AddText("CC source ("+SystemName+"):");
        PT_CC[iSeMe]->AddText(TString::Format("#chi^{2}_{ndf}=%.2f\n",Chi2Ndf_CC[iSeMe]));
        PT_CC[iSeMe]->AddText(TString::Format("RC1/RC2=%.3f/%.3f fm\n",fitRad_CC[iSeMe]->GetParameter(0),fitRad_CC[iSeMe]->GetParameter(1)));
        PT_CC[iSeMe]->AddText(TString::Format("wC1/wC2=%.3f/%.3f fm\n",fitRad_CC[iSeMe]->GetParameter(2),1.-fitRad_CC[iSeMe]->GetParameter(2)));
        PT_CC[iSeMe]->Draw("same");


        cCleverLevy[iSeMe] = new TCanvas("cCleverLevy_"+AddOn, "cCleverLevy_"+AddOn, 1);
        cCleverLevy[iSeMe]->cd(0);
        cCleverLevy[iSeMe]->SetCanvasSize(1280, 720);
        cCleverLevy[iSeMe]->SetMargin(0.15,0.05,0.2,0.05);//lrbt
        hRad[iSeMe]->Draw();
        fitRad_CleverLevy[iSeMe]->Draw("same");
        myLegend[iSeMe]->Draw("same");

        PT_CleverLevy[iSeMe] = new TPaveText(0.65,0.55,0.95,0.80, "blNDC");//lbrt
        PT_CleverLevy[iSeMe]->SetName("PT_CleverLevy"+AddOn);
        PT_CleverLevy[iSeMe]->SetBorderSize(1);
        PT_CleverLevy[iSeMe]->SetTextSize(0.045);
        PT_CleverLevy[iSeMe]->SetFillColor(kWhite);
        PT_CleverLevy[iSeMe]->SetTextFont(22);
        PT_CleverLevy[iSeMe]->AddText("Levy source ("+SystemName+"):");
        PT_CleverLevy[iSeMe]->AddText(TString::Format("#chi^{2}_{ndf}=%.2f\n",Chi2Ndf_CC[iSeMe]));
        PT_CleverLevy[iSeMe]->AddText(TString::Format("#sigma=%.3f#pm %.3f fm\n",fitRad_CleverLevy[iSeMe]->GetParameter(0),fitRad_CleverLevy[iSeMe]->GetParError(0)));
        PT_CleverLevy[iSeMe]->AddText(TString::Format("#alpha=%.3f#pm %.3f fm\n",fitRad_CleverLevy[iSeMe]->GetParameter(1),fitRad_CleverLevy[iSeMe]->GetParError(1)));
        PT_CleverLevy[iSeMe]->Draw("same");


        hRad[iSeMe]->Write();
        hMomDist[iSeMe]->Write();
        fitRad_Gauss[iSeMe]->Write();
        fitRad_Cauchy[iSeMe]->Write();
        fitRad_GC[iSeMe]->Write();
        fitRad_GG[iSeMe]->Write();
        fitRad_CC[iSeMe]->Write();
        fitRad_CleverLevy[iSeMe]->Write();
        cGauss[iSeMe]->Write();
        cCauchy[iSeMe]->Write();
        cGC[iSeMe]->Write();
        cGG[iSeMe]->Write();
        cCC[iSeMe]->Write();
        cCleverLevy[iSeMe]->Write();

        cGauss[iSeMe]->SaveAs(OutputFolder+"cGauss_"+SubSample+"_"+AddOn+".png");
        cCauchy[iSeMe]->SaveAs(OutputFolder+"cCauchy_"+SubSample+"_"+AddOn+".png");
        cGC[iSeMe]->SaveAs(OutputFolder+"cGC_"+SubSample+"_"+AddOn+".png");
        cGG[iSeMe]->SaveAs(OutputFolder+"cGG_"+SubSample+"_"+AddOn+".png");
        cCC[iSeMe]->SaveAs(OutputFolder+"cCC_"+SubSample+"_"+AddOn+".png");
        cCleverLevy[iSeMe]->SaveAs(OutputFolder+"cCleverLevy_"+SubSample+"_"+AddOn+".png");
    }

    TF1** fitRad_Gauss_SP = new TF1* [2];
    TF1** fitRad_Cauchy_SP = new TF1* [2];
    TF1** fitRad_GC_SP = new TF1* [2];
    TF1** fitRad_GG_SP = new TF1* [2];
    TF1** fitRad_CC_SP = new TF1* [2];

    double Chi2Ndf_Gauss_SP[2];
    double Chi2Ndf_Cauchy_SP[2];
    double Chi2Ndf_GC_SP[2];
    double Chi2Ndf_GG_SP[2];
    double Chi2Ndf_CC_SP[2];

    TLegend** myLegend_SP = new TLegend*[2];
    TCanvas** cGauss_SP = new TCanvas*[2];
    TPaveText** PT_Gauss_SP = new TPaveText*[2];

    TCanvas** cCauchy_SP = new TCanvas*[2];
    TPaveText** PT_Cauchy_SP = new TPaveText*[2];

    TCanvas** cGC_SP = new TCanvas*[2];
    TPaveText** PT_GC_SP = new TPaveText*[2];

    TCanvas** cGG_SP = new TCanvas*[2];
    TPaveText** PT_GG_SP = new TPaveText*[2];

    TCanvas** cCC_SP = new TCanvas*[2];
    TPaveText** PT_CC_SP = new TPaveText*[2];

    for(int iType=0; iType<2; iType++){
        TString AddOn;
        TString ParticleName;
        if(System=="pp"){
            if(iType==0) {AddOn = "p"; ParticleName = "p";}
            else {AddOn = "p2"; ParticleName = "p";}        }
        else if(System=="pL"){
            if(iType==0) {AddOn = "p"; ParticleName = "p";}
            else {AddOn = "L"; ParticleName = "#Lambda";}
        }
        else if(System=="LL"){
            if(iType==0) {AddOn = "L"; ParticleName = "#Lambda";}
            else {AddOn = "L2"; ParticleName = "#Lambda";}
        }
        else if(System=="pXim"){
            if(iType==0) {AddOn = "p"; ParticleName = "p";}
            else {AddOn = "Xim";  ParticleName = "#Xi^{#minus}";}
        }

        if(TauCorrection) AddOn += "_Tau";
        else AddOn += "_NoTau";

        if(BoostToCm) AddOn += "_Cm";
        else AddOn += "_Lab";

        fitRad_Gauss_SP[iType] = new TF1("fitRad_Gauss_"+AddOn,fitterRad_Gauss,RadMin,RadMax,1);
        fitRad_Cauchy_SP[iType] = new TF1("fitRad_Cauchy_"+AddOn,fitterRad_Cauchy,RadMin,RadMax,1);
        fitRad_GC_SP[iType] = new TF1("fitRad_GC_"+AddOn,fitterRad_GC,RadMin,RadMax,3);
        fitRad_GG_SP[iType] = new TF1("fitRad_GG_"+AddOn,fitterRad_GG,RadMin,RadMax,3);
        fitRad_CC_SP[iType] = new TF1("fitRad_CC_"+AddOn,fitterRad_CC,RadMin,RadMax,3);

        fitRad_Gauss_SP[iType]->SetParameter(0,1);
        fitRad_Gauss_SP[iType]->SetParLimits(0,0.1,2);
        fitRad_Gauss_SP[iType]->SetLineColor(kRed);
        fitRad_Gauss_SP[iType]->SetLineWidth(4);

        fitRad_Cauchy_SP[iType]->SetParameter(0,1);
        fitRad_Cauchy_SP[iType]->SetParLimits(0,0.1,2);
        fitRad_Cauchy_SP[iType]->SetLineColor(kRed);
        fitRad_Cauchy_SP[iType]->SetLineWidth(4);

        fitRad_GC_SP[iType]->SetParameter(0,1);
        fitRad_GC_SP[iType]->SetParLimits(0,0.1,2);
        fitRad_GC_SP[iType]->SetParameter(1,1);
        fitRad_GC_SP[iType]->SetParLimits(1,0.1,2);
        fitRad_GC_SP[iType]->SetParameter(2,0.5);
        fitRad_GC_SP[iType]->SetParLimits(2,0,1);
        fitRad_GC_SP[iType]->SetLineColor(kRed);
        fitRad_GC_SP[iType]->SetLineWidth(4);

        fitRad_GG_SP[iType]->SetParameter(0,1);
        fitRad_GG_SP[iType]->SetParLimits(0,0.1,2);
        fitRad_GG_SP[iType]->SetParameter(1,1);
        fitRad_GG_SP[iType]->SetParLimits(1,0.1,2);
        fitRad_GG_SP[iType]->SetParameter(2,0.5);
        fitRad_GG_SP[iType]->SetParLimits(2,0,1);
        fitRad_GG_SP[iType]->SetLineColor(kRed);
        fitRad_GG_SP[iType]->SetLineWidth(4);

        fitRad_CC_SP[iType]->SetParameter(0,1);
        fitRad_CC_SP[iType]->SetParLimits(0,0.1,2);
        fitRad_CC_SP[iType]->SetParameter(1,1);
        fitRad_CC_SP[iType]->SetParLimits(1,0.1,2);
        fitRad_CC_SP[iType]->SetParameter(2,0.5);
        fitRad_CC_SP[iType]->SetParLimits(2,0,1);
        fitRad_CC_SP[iType]->SetLineColor(kRed);
        fitRad_CC_SP[iType]->SetLineWidth(4);

        hSingleParticle[iType]->Fit(fitRad_Gauss_SP[iType], "Q, S, N, R, M");
        hSingleParticle[iType]->Fit(fitRad_Cauchy_SP[iType], "Q, S, N, R, M");
        hSingleParticle[iType]->Fit(fitRad_GC_SP[iType], "Q, S, N, R, M");
        hSingleParticle[iType]->Fit(fitRad_GG_SP[iType], "Q, S, N, R, M");
        hSingleParticle[iType]->Fit(fitRad_CC_SP[iType], "Q, S, N, R, M");

        hSingleParticle[iType]->SetStats(false);
        hSingleParticle[iType]->SetTitle("");
        hSingleParticle[iType]->GetXaxis()->SetLabelSize(0.065);
        hSingleParticle[iType]->GetXaxis()->SetTitle("r (fm)");
        hSingleParticle[iType]->GetXaxis()->CenterTitle();
        hSingleParticle[iType]->GetXaxis()->SetTitleOffset(1.25);
        hSingleParticle[iType]->GetXaxis()->SetLabelOffset(0.02);
        hSingleParticle[iType]->GetXaxis()->SetTitleSize(0.075);
        hSingleParticle[iType]->GetYaxis()->SetLabelSize(0.065);
        hSingleParticle[iType]->GetYaxis()->SetTitle("4#pi r^{2}S(r) (1/fm)");
        hSingleParticle[iType]->GetYaxis()->CenterTitle();
        hSingleParticle[iType]->GetYaxis()->SetTitleOffset(0.8);
        hSingleParticle[iType]->GetYaxis()->SetTitleSize(0.075);
        hSingleParticle[iType]->GetYaxis()->SetRangeUser(0, 0.65);
        hSingleParticle[iType]->SetLineColor(kBlack);
        hSingleParticle[iType]->SetLineWidth(3);

        Chi2Ndf_Gauss_SP[iType] = fitRad_Gauss_SP[iType]->GetChisquare()/double(fitRad_Gauss_SP[iType]->GetNDF());
        Chi2Ndf_Cauchy_SP[iType] = fitRad_Cauchy_SP[iType]->GetChisquare()/double(fitRad_Cauchy_SP[iType]->GetNDF());
        Chi2Ndf_GC_SP[iType] = fitRad_GC_SP[iType]->GetChisquare()/double(fitRad_GC_SP[iType]->GetNDF());
        Chi2Ndf_GG_SP[iType] = fitRad_GG_SP[iType]->GetChisquare()/double(fitRad_GG_SP[iType]->GetNDF());
        Chi2Ndf_CC_SP[iType] = fitRad_CC_SP[iType]->GetChisquare()/double(fitRad_CC_SP[iType]->GetNDF());
//printf("Gauss Chi2/NDF %f/%u\n",fitRad_Gauss_SP[iType]->GetChisquare(),fitRad_Gauss_SP[iType]->GetNDF());
        printf("%s:\n",AddOn.Data());
        printf("Chi2Ndf_Gauss = %.2f\n",Chi2Ndf_Gauss_SP[iType]);
        printf("Chi2Ndf_Cauchy = %.2f\n",Chi2Ndf_Cauchy_SP[iType]);
        printf("Chi2Ndf_GC = %.2f\n",Chi2Ndf_GC_SP[iType]);
        printf("Chi2Ndf_GG = %.2f\n",Chi2Ndf_GG_SP[iType]);
        printf("Chi2Ndf_CC = %.2f\n",Chi2Ndf_CC_SP[iType]);
        printf("\n");

        myLegend_SP[iType] = new TLegend(0.65,0.80,0.95,0.95);//lbrt
        myLegend_SP[iType]->SetName("myLegend_SB_"+AddOn);
        myLegend_SP[iType]->SetTextSize(0.040);
        myLegend_SP[iType]->AddEntry(hSingleParticle[iType], TString::Format("EPOS source (%s)",SourceDiscr.Data()));
        myLegend_SP[iType]->AddEntry(fitRad_Gauss_SP[iType], "Fit");

        cGauss_SP[iType] = new TCanvas("cGauss_SP_"+AddOn, "cGauss_SP_"+AddOn, 1);
        cGauss_SP[iType]->cd(0);
        cGauss_SP[iType]->SetCanvasSize(1280, 720);
        cGauss_SP[iType]->SetMargin(0.15,0.05,0.2,0.05);//lrbt
        hSingleParticle[iType]->Draw();
        fitRad_Gauss_SP[iType]->Draw("same");
        myLegend_SP[iType]->Draw("same");

        PT_Gauss_SP[iType] = new TPaveText(0.65,0.55,0.95,0.80, "blNDC");//lbrt
        PT_Gauss_SP[iType]->SetName("PT_Gauss_SP"+AddOn);
        PT_Gauss_SP[iType]->SetBorderSize(1);
        PT_Gauss_SP[iType]->SetTextSize(0.045);
        PT_Gauss_SP[iType]->SetFillColor(kWhite);
        PT_Gauss_SP[iType]->SetTextFont(22);
        PT_Gauss_SP[iType]->AddText("Gauss source ("+ParticleName+"):");
        PT_Gauss_SP[iType]->AddText(TString::Format("#chi^{2}_{ndf}=%.2f\n",Chi2Ndf_Gauss_SP[iType]));
        PT_Gauss_SP[iType]->AddText(TString::Format("R=%.3f fm\n",fitRad_Gauss_SP[iType]->GetParameter(0)));
        PT_Gauss_SP[iType]->Draw("same");


        cCauchy_SP[iType] = new TCanvas("cCauchy_SP_"+AddOn, "cCauchy_SP_"+AddOn, 1);
        cCauchy_SP[iType]->cd(0);
        cCauchy_SP[iType]->SetCanvasSize(1280, 720);
        cCauchy_SP[iType]->SetMargin(0.15,0.05,0.2,0.05);//lrbt
        hSingleParticle[iType]->Draw();
        fitRad_Cauchy_SP[iType]->Draw("same");
        myLegend_SP[iType]->Draw("same");

        PT_Cauchy_SP[iType] = new TPaveText(0.65,0.55,0.95,0.80, "blNDC");//lbrt
        PT_Cauchy_SP[iType]->SetName("PT_Cauchy_SP");
        PT_Cauchy_SP[iType]->SetBorderSize(1);
        PT_Cauchy_SP[iType]->SetTextSize(0.045);
        PT_Cauchy_SP[iType]->SetFillColor(kWhite);
        PT_Cauchy_SP[iType]->SetTextFont(22);
        PT_Cauchy_SP[iType]->AddText("Cauchy source ("+ParticleName+"):");
        PT_Cauchy_SP[iType]->AddText(TString::Format("#chi^{2}_{ndf}=%.2f\n",Chi2Ndf_Cauchy_SP[iType]));
        PT_Cauchy_SP[iType]->AddText(TString::Format("R=%.3f fm\n",fitRad_Cauchy_SP[iType]->GetParameter(0)));
        PT_Cauchy_SP[iType]->Draw("same");


        cGC_SP[iType] = new TCanvas("cGC_SP_"+AddOn, "cGC_SP_"+AddOn, 1);
        cGC_SP[iType]->cd(0);
        cGC_SP[iType]->SetCanvasSize(1280, 720);
        cGC_SP[iType]->SetMargin(0.15,0.05,0.2,0.05);//lrbt
        hSingleParticle[iType]->Draw();
        fitRad_GC_SP[iType]->Draw("same");
        myLegend_SP[iType]->Draw("same");

        PT_GC_SP[iType] = new TPaveText(0.65,0.55,0.95,0.80, "blNDC");//lbrt
        PT_GC_SP[iType]->SetName("PT_GC_SP");
        PT_GC_SP[iType]->SetBorderSize(1);
        PT_GC_SP[iType]->SetTextSize(0.045);
        PT_GC_SP[iType]->SetFillColor(kWhite);
        PT_GC_SP[iType]->SetTextFont(22);
        PT_GC_SP[iType]->AddText("GC source ("+ParticleName+"):");
        PT_GC_SP[iType]->AddText(TString::Format("#chi^{2}_{ndf}=%.2f\n",Chi2Ndf_GC_SP[iType]));
        PT_GC_SP[iType]->AddText(TString::Format("R=%.3f fm\n",fitRad_GC_SP[iType]->GetParameter(0)));
        PT_GC_SP[iType]->Draw("same");

        cGG_SP[iType] = new TCanvas("cGG_SP_"+AddOn, "cGG_SP_"+AddOn, 1);
        cGG_SP[iType]->cd(0);
        cGG_SP[iType]->SetCanvasSize(1280, 720);
        cGG_SP[iType]->SetMargin(0.15,0.05,0.2,0.05);//lrbt
        hSingleParticle[iType]->Draw();
        fitRad_GG_SP[iType]->Draw("same");
        myLegend_SP[iType]->Draw("same");

        PT_GG_SP[iType] = new TPaveText(0.65,0.55,0.95,0.80, "blNDC");//lbrt
        PT_GG_SP[iType]->SetName("PT_GG_SP");
        PT_GG_SP[iType]->SetBorderSize(1);
        PT_GG_SP[iType]->SetTextSize(0.045);
        PT_GG_SP[iType]->SetFillColor(kWhite);
        PT_GG_SP[iType]->SetTextFont(22);
        PT_GG_SP[iType]->AddText("GG source ("+ParticleName+"):");
        PT_GG_SP[iType]->AddText(TString::Format("#chi^{2}_{ndf}=%.2f\n",Chi2Ndf_GG_SP[iType]));
        PT_GG_SP[iType]->AddText(TString::Format("R=%.3f fm\n",fitRad_GG_SP[iType]->GetParameter(0)));
        PT_GG_SP[iType]->Draw("same");

        cCC_SP[iType] = new TCanvas("cCC_SP_"+AddOn, "cCC_SP_"+AddOn, 1);
        cCC_SP[iType]->cd(0);
        cCC_SP[iType]->SetCanvasSize(1280, 720);
        cCC_SP[iType]->SetMargin(0.15,0.05,0.2,0.05);//lrbt
        hSingleParticle[iType]->Draw();
        fitRad_CC_SP[iType]->Draw("same");
        myLegend_SP[iType]->Draw("same");

        PT_CC_SP[iType] = new TPaveText(0.65,0.55,0.95,0.80, "blNDC");//lbrt
        PT_CC_SP[iType]->SetName("PT_CC_SP");
        PT_CC_SP[iType]->SetBorderSize(1);
        PT_CC_SP[iType]->SetTextSize(0.045);
        PT_CC_SP[iType]->SetFillColor(kWhite);
        PT_CC_SP[iType]->SetTextFont(22);
        PT_CC_SP[iType]->AddText("CC source ("+ParticleName+"):");
        PT_CC_SP[iType]->AddText(TString::Format("#chi^{2}_{ndf}=%.2f\n",Chi2Ndf_CC_SP[iType]));
        PT_CC_SP[iType]->AddText(TString::Format("R=%.3f fm\n",fitRad_CC_SP[iType]->GetParameter(0)));
        PT_CC_SP[iType]->Draw("same");

        hSingleParticle[iType]->Write();
        //fitRad_Gauss_SP[iType]->Write();
        //fitRad_Cauchy_SP[iType]->Write();
        //fitRad_GC_SP[iType]->Write();
        //fitRad_GG_SP[iType]->Write();
        //fitRad_CC_SP[iType]->Write();
        cGauss_SP[iType]->Write();
        cCauchy_SP[iType]->Write();
        cGC_SP[iType]->Write();
        cGG_SP[iType]->Write();
        cCC_SP[iType]->Write();

        cGauss_SP[iType]->SaveAs(OutputFolder+"cGauss_"+SubSample+"_"+AddOn+".png");
        cCauchy_SP[iType]->SaveAs(OutputFolder+"cCauchy_"+SubSample+"_"+AddOn+".png");
        cGC_SP[iType]->SaveAs(OutputFolder+"cGC_"+SubSample+"_"+AddOn+".png");
        cGG_SP[iType]->SaveAs(OutputFolder+"cGG_"+SubSample+"_"+AddOn+".png");
        cCC_SP[iType]->SaveAs(OutputFolder+"cCC_"+SubSample+"_"+AddOn+".png");
    }

    for(int iSeMe=0; iSeMe<2; iSeMe++){
        delete hRad[iSeMe];
        delete hMomDist[iSeMe];
        delete fitRad_Gauss[iSeMe];
        delete fitRad_Cauchy[iSeMe];
        delete fitRad_GC[iSeMe];
        delete fitRad_GG[iSeMe];
        delete fitRad_CC[iSeMe];
        delete PT_Gauss[iSeMe];
        delete PT_Cauchy[iSeMe];
        delete PT_GC[iSeMe];
        delete PT_GG[iSeMe];
        delete PT_CC[iSeMe];
        delete myLegend[iSeMe];

        delete cGauss[iSeMe];
        delete cCauchy[iSeMe];
        delete cGC[iSeMe];
        delete cGG[iSeMe];
        delete cCC[iSeMe];
    }
    delete [] hRad;
    delete [] hMomDist;
    delete [] fitRad_Gauss;
    delete [] fitRad_Cauchy;
    delete [] fitRad_GC;
    delete [] fitRad_GG;
    delete [] fitRad_CC;
    delete [] PT_Gauss;
    delete [] PT_Cauchy;
    delete [] PT_GC;
    delete [] PT_GG;
    delete [] PT_CC;
    delete [] myLegend;
    delete [] cGauss;
    delete [] cCauchy;
    delete [] cGC;
    delete [] cGG;
    delete [] cCC;

//return;
    for(int iType=0; iType<2; iType++){
        delete hSingleParticle[iType];
        delete fitRad_Gauss_SP[iType];
        delete fitRad_Cauchy_SP[iType];
        delete fitRad_GC_SP[iType];
        delete fitRad_GG_SP[iType];
        delete fitRad_CC_SP[iType];
        delete PT_Gauss_SP[iType];
        delete PT_Cauchy_SP[iType];
        delete PT_GC_SP[iType];
        delete PT_GG_SP[iType];
        delete PT_CC_SP[iType];
        delete myLegend_SP[iType];

        delete cGauss_SP[iType];
        delete cCauchy_SP[iType];
        delete cGC_SP[iType];
        delete cGG_SP[iType];
        delete cCC_SP[iType];
    }

    delete [] hSingleParticle;
    delete [] fitRad_Gauss_SP;
    delete [] fitRad_Cauchy_SP;
    delete [] fitRad_GC_SP;
    delete [] fitRad_GG_SP;
    delete [] fitRad_CC_SP;
    delete [] PT_Gauss_SP;
    delete [] PT_Cauchy_SP;
    delete [] PT_GC_SP;
    delete [] PT_GG_SP;
    delete [] PT_CC_SP;
    delete [] myLegend_SP;
    delete [] cGauss_SP;
    delete [] cCauchy_SP;
    delete [] cGC_SP;
    delete [] cGG_SP;
    delete [] cCC_SP;

    delete fOutput;
}

void GetSource_Compare(){

}

int GerhardMAIN(int narg, char** ARGS){
    printf("Running GerhardMAIN\n");
    //EposFit_pp_1();
    //FitEposSource_pp("pp","Full");
    //FitEposSource_pp("pp","Prim");
    //FitEposSource_pp("pL","Full");
    //FitEposSource_pp("pL","Prim");

    const double KMIN[12] = {   0,  50, 100, 150, 200, 250, 150, 500,1000,2000,   0,   0};
    const double KMAX[12] = {  50, 100, 150, 200, 250, 300, 500,1000,2000,8000, 150,8000};

    //GetSource("pp","Full",false,true,0,8000);
    //GetSource("pp","Prim",false,true,0,8000);

    //GetSource("pp","Full",false,false,0,8000);
    //GetSource("pp","Prim",false,false,0,8000);

    //GetSource("pL","Full",false,true,0,150);
    //GetSource("pL","Prim",false,true,0,150);
//return 0;
    for(int iMom=0; iMom<12; iMom++){
if(iMom!=10) continue;
        GetSource("pp","Full",false,true,KMIN[iMom],KMAX[iMom]);
        GetSource("pp","Prim",false,true,KMIN[iMom],KMAX[iMom]);
        GetSource("pL","Full",false,true,KMIN[iMom],KMAX[iMom]);
        GetSource("pL","Prim",false,true,KMIN[iMom],KMAX[iMom]);

        GetSource("pp","Full",true,true,KMIN[iMom],KMAX[iMom]);
        GetSource("pp","Prim",true,true,KMIN[iMom],KMAX[iMom]);
        GetSource("pL","Full",true,true,KMIN[iMom],KMAX[iMom]);
        GetSource("pL","Prim",true,true,KMIN[iMom],KMAX[iMom]);

        GetSource("pp","Full",false,false,KMIN[iMom],KMAX[iMom]);
        GetSource("pp","Prim",false,false,KMIN[iMom],KMAX[iMom]);
        GetSource("pL","Full",false,false,KMIN[iMom],KMAX[iMom]);
        GetSource("pL","Prim",false,false,KMIN[iMom],KMAX[iMom]);

        GetSource("pp","Full",true,false,KMIN[iMom],KMAX[iMom]);
        GetSource("pp","Prim",true,false,KMIN[iMom],KMAX[iMom]);
        GetSource("pL","Full",true,false,KMIN[iMom],KMAX[iMom]);
        GetSource("pL","Prim",true,false,KMIN[iMom],KMAX[iMom]);

    }

    return 0;
}
