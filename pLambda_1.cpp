
#include "pLambda_1.h"
#include "FemtoBoyzScripts.h"

#include "TString.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TNtuple.h"
#include "TRandom3.h"
#include "TLatex.h"
#include "TStyle.h"

#include "CATS.h"
#include "DLM_CkDecomposition.h"
#include "DLM_CkModels.h"
#include "CommonAnaFunctions.h"
#include "DLM_Fitters.h"

#include <unistd.h>
#include <fcntl.h>
#include <cerrno>

TH1F* pL_CleverRebin(TH1F* hData_pL, const unsigned& NumMomBins_pL, const double* MomBins_pL){
    //the one with var binning
    //unsigned NumMomBinsData_pL = NumMomBins_pL+20;
    unsigned NumMomBinsData_pL = hData_pL->GetNbinsX();
    double* MomBinsData_pL = new double [NumMomBinsData_pL+1];
    double StepSize=0;
    double StepSize0=0;
    for(unsigned uBin=0; uBin<=NumMomBinsData_pL; uBin++){
        if(!uBin){StepSize=MomBins_pL[uBin+1]-MomBins_pL[uBin];}
        if(uBin<=NumMomBins_pL){
            MomBinsData_pL[uBin]=MomBins_pL[uBin];
        }
        else{
            MomBinsData_pL[uBin]=MomBinsData_pL[uBin-1]+StepSize;
        }
    }
    TH1F* hDataClever_pL = new TH1F("hDataClever_pL","hDataClever_pL",NumMomBinsData_pL,MomBinsData_pL);
    for(unsigned uBin=1; uBin<=NumMomBinsData_pL; uBin++){
        StepSize=hDataClever_pL->GetBinWidth(uBin);
        int Which0bin = hData_pL->FindBin(hDataClever_pL->GetBinLowEdge(uBin));
        int Initial0bin = Which0bin;
        StepSize0=hData_pL->GetBinWidth(Which0bin);
        if(fabs(StepSize-StepSize0)<1e-6){
            hDataClever_pL->SetBinContent(uBin,hData_pL->GetBinContent(Which0bin));
            hDataClever_pL->SetBinError(uBin,hData_pL->GetBinError(Which0bin));
        }
        else{
            while(hData_pL->GetBinCenter(Which0bin)<hDataClever_pL->GetXaxis()->GetBinUpEdge(uBin)){
                Which0bin++;
            }
            double BinError=0;
            double BinContent=0;
            for(int uFineBin=Initial0bin; uFineBin<Which0bin; uFineBin++){
                BinContent+=hData_pL->GetBinContent(uFineBin);
                BinError+=hData_pL->GetBinError(uFineBin)*hData_pL->GetBinError(uFineBin);
            }
            hDataClever_pL->SetBinContent(uBin,BinContent/double(Which0bin-Initial0bin));
            hDataClever_pL->SetBinError(uBin,sqrt(BinError)/double(Which0bin-Initial0bin));
        }
    }

//TFile* fTemp = new TFile(OutputFolder+"fTemp.root","recreate");
//hData_pL->Write();
//hDataClever_pL->Write();
//delete fTemp;
    return hDataClever_pL;
}

//a fitter to fit pp and pL using a specific source
//we can fit them together or separately
//For feed-down we will use only pSigma0 and pXi, where the former we model as Oli,
//the latter we model with HALQCD, by fixing the source and radius to a Gauss of size 1.18 (-20% for the systematic variation)
//Fitting mode:
    //"Norm": only norm
    //"Baseline": baseline but only in the femto region
    //"Longbaseline": baseline extending to higher k
    //"Longbaseline_prefit": prefit the baseline assuming Cl=1

void Fit_pp_pL(const TString& DataSample, const TString& SourceType, const TString& FittingMode_pp, const TString& FittingMode_pL){

    DLM_CommonAnaFunctions AnalysisObject;

    const TString OutputFolder = TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pp_pL/Test/");

    const bool SameSourceSize = false;
    const bool SameStability = false;

    //the residuals that are not fitted, are assumed to all have a Gaussian source of some size;
    double ResidualSourceSize=0;
    if(DataSample=="pp13TeV_MB_Run2paper") ResidualSourceSize=1.1;
    else if(DataSample=="pPb5TeV_Run2paper"||DataSample=="pPb5TeV_CPR_Mar19") ResidualSourceSize=1.4;
    else if(DataSample=="pp13TeV_HM_March19") ResidualSourceSize=1.4;
    else printf("\033[1;31mERROR:\033[0m The data sample '%s' does not exist\n",DataSample.Data());

    const unsigned Num_pp_pots = 1;
    const unsigned Num_pL_pots = 6;
    const unsigned Num_pXim_pots = 2;

    const TString pp_Pot[Num_pp_pots] =     {"AV18"};
    const bool Use_pp_Pot[Num_pp_pots] =    {true};
    const TString pL_Pot[Num_pL_pots] =     {"LO",  "LO_Coupled_S", "NLO", "NLO_sp",  "NLO_Coupled_S",    "Usmani"};
    const bool Use_pL_Pot[Num_pL_pots] =    {false,  false,          false,   false,     false,              true};
    const TString pXim_Pot[Num_pXim_pots] = {"pXim_HALQCD1","pXim_Lattice"};
    const bool Use_pXim_Pot[Num_pXim_pots] =  {true,          false};

    double lambda_pp[4];
    AnalysisObject.SetUpLambdaPars_pp(DataSample,0,lambda_pp);

    double lambda_pL[5];
    AnalysisObject.SetUpLambdaPars_pL(DataSample,0,0,lambda_pL);

    double lambda_pXim[5];
    AnalysisObject.SetUpLambdaPars_pXim(DataSample,0,0,lambda_pXim);

    double* MomBins_pp = NULL;
    double* FitRegion_pp = NULL;
    unsigned NumMomBins_pp;

    AnalysisObject.SetUpBinning_pp(DataSample,NumMomBins_pp,MomBins_pp,FitRegion_pp);

    double* MomBins_pL = NULL;
    double* FitRegion_pL = NULL;
    unsigned NumMomBins_pL;
    AnalysisObject.SetUpBinning_pL(DataSample,NumMomBins_pL,MomBins_pL,FitRegion_pL,0,0);

    DLM_Ck* Ck_pSigma0 = new DLM_Ck(1,0,NumMomBins_pL,MomBins_pL,Lednicky_gauss_Sigma0);
    Ck_pSigma0->SetSourcePar(0,ResidualSourceSize);

    TH2F* hResolution_pp = AnalysisObject.GetResolutionMatrix(DataSample,"pp");
    TH2F* hResolution_pL = AnalysisObject.GetResolutionMatrix(DataSample,"pLambda");

    TH2F* hResidual_pp_pL = AnalysisObject.GetResidualMatrix("pp","pLambda");
    TH2F* hResidual_pL_pSigma0 = AnalysisObject.GetResidualMatrix("pLambda","pSigma0");
    TH2F* hResidual_pL_pXim = AnalysisObject.GetResidualMatrix("pLambda","pXim");
//GetAliceExpCorrFun(const TString& DataSample,const TString& System,const TString& CutVar,const int& iReb, const bool& AddSyst=false,const int mTbin=-1)
    TH1F* hData_pp = AnalysisObject.GetAliceExpCorrFun(DataSample,"pp","_0",0);
//! this has to be changed for finer binning
    TH1F* hData_pL = AnalysisObject.GetAliceExpCorrFun(DataSample,"pLambda","_0",2);

    //TH1F* hDataClever_pL = pL_CleverRebin(hData_pL,NumMomBins_pL,MomBins_pL);

    double lam_pp[4];
    double lam_pL[5];
    double lam_pXim[5];

    AnalysisObject.SetUpLambdaPars_pp(DataSample,0,lam_pp);
    AnalysisObject.SetUpLambdaPars_pL(DataSample,0,0,lam_pL);
    AnalysisObject.SetUpLambdaPars_pXim(DataSample,0,0,lam_pXim);

    TFile* OutputFile = new TFile(OutputFolder+TString::Format("Output_%s_%s_%s_%s.root",DataSample.Data(),SourceType.Data(),
                                                               FittingMode_pp.Data(),FittingMode_pL.Data()),"recreate");
    hData_pp->Write();
    hData_pL->Write();
    //hDataClever_pL->Write();

    //iterate over all potentials
    for(unsigned uPot_pp=0; uPot_pp<Num_pp_pots; uPot_pp++){
        if(!Use_pp_Pot[uPot_pp]) continue;
        CATS AB_pp;
        AB_pp.SetMomBins(NumMomBins_pp,MomBins_pp);
        AnalysisObject.SetUpCats_pp(AB_pp,pp_Pot[uPot_pp],SourceType,0,2);
        AB_pp.SetNotifications(CATS::nWarning);
        AB_pp.KillTheCat();
        DLM_Ck* Ck_pp = new DLM_Ck(AB_pp.GetNumSourcePars(),0,AB_pp);
        for(unsigned uPot_pL=0; uPot_pL<Num_pL_pots; uPot_pL++){
            if(!Use_pL_Pot[uPot_pL]) continue;
            CATS AB_pL;
            AB_pL.SetMomBins(NumMomBins_pL,MomBins_pL);
            AnalysisObject.SetUpCats_pL(AB_pL,pL_Pot[uPot_pL],SourceType,0,3);
            AB_pL.SetNotifications(CATS::nWarning);
            AB_pL.KillTheCat();
            DLM_Ck* Ck_pL = new DLM_Ck(AB_pL.GetNumSourcePars(),0,AB_pL);
            for(unsigned uPot_pXim=0; uPot_pXim<Num_pXim_pots; uPot_pXim++){
                if(!Use_pXim_Pot[uPot_pXim]) continue;
                printf("Current iteration:----------------------\n");
                printf("  Data sample: %s\n",DataSample.Data());
                printf("  Source: %s\n",SourceType.Data());
                printf("  Within this iteration I am at:\n");
                printf("    pp: %s\n",pp_Pot[uPot_pp].Data());
                printf("    pΛ: %s\n",pL_Pot[uPot_pL].Data());
                printf("    pΞ: %s\n",pXim_Pot[uPot_pXim].Data());
                CATS AB_pXim;
                //same binning as pL, as we only use pXim as feed-down
                AB_pXim.SetMomBins(NumMomBins_pL,MomBins_pL);
                AnalysisObject.SetUpCats_pXim(AB_pXim,pXim_Pot[uPot_pXim],"Gauss");
                AB_pXim.SetAnaSource(0,ResidualSourceSize);
                AB_pXim.SetNotifications(CATS::nWarning);
                AB_pXim.KillTheCat();
                DLM_Ck* Ck_pXim = new DLM_Ck(AB_pXim.GetNumSourcePars(),0,AB_pXim);

                Ck_pp->Update();
                Ck_pL->Update();
                Ck_pSigma0->Update();
                Ck_pXim->Update();

                DLM_CkDecomposition CkDec_pp("pp",3,*Ck_pp,hResolution_pp);
                DLM_CkDecomposition CkDec_pL("pLambda",4,*Ck_pL,hResolution_pL);
                DLM_CkDecomposition CkDec_pSigma0("pSigma0",0,*Ck_pSigma0,NULL);
                DLM_CkDecomposition CkDec_pXim("pXim",2,*Ck_pXim,NULL);

                CkDec_pp.AddContribution(0,lam_pp[1],DLM_CkDecomposition::cFeedDown,&CkDec_pL,hResidual_pp_pL);
                CkDec_pp.AddContribution(1,lam_pp[2],DLM_CkDecomposition::cFeedDown);
                CkDec_pp.AddContribution(2,lam_pp[3],DLM_CkDecomposition::cFake);

                CkDec_pL.AddContribution(0,lam_pL[1],DLM_CkDecomposition::cFeedDown,&CkDec_pSigma0,hResidual_pL_pSigma0);
                CkDec_pL.AddContribution(1,lam_pL[2],DLM_CkDecomposition::cFeedDown,&CkDec_pXim,hResidual_pL_pXim);
                CkDec_pL.AddContribution(2,lam_pL[3],DLM_CkDecomposition::cFeedDown);
                CkDec_pL.AddContribution(3,lam_pL[4],DLM_CkDecomposition::cFake);//0.03

                //for Xim we simplify a bit and take ALL feed-down as flat
                CkDec_pXim.AddContribution(0,lam_pXim[1]+lam_pXim[2]+lam_pXim[3],DLM_CkDecomposition::cFeedDown);
                CkDec_pXim.AddContribution(1,lam_pXim[4],DLM_CkDecomposition::cFake);

                DLM_Fitter1* fitter = new DLM_Fitter1(2);

                if(FittingMode_pp=="Norm"||FittingMode_pp.Contains("Baseline")){
//printf("FitRegion_pp[0]=%f FitRegion_pp[1]=%f\n",FitRegion_pp[0],FitRegion_pp[1]);
                    fitter->SetSystem(0,*hData_pp,1,CkDec_pp,
                                  FitRegion_pp[0],FitRegion_pp[1],FitRegion_pp[1],FitRegion_pp[1]);
                }
                else if(FittingMode_pp.Contains("Longbaseline")){
                    fitter->SetSystem(0,*hData_pp,1,CkDec_pp,
                                  FitRegion_pp[0],FitRegion_pp[1],FitRegion_pp[2],FitRegion_pp[3]);
                }
                else{
                    printf("\033[1;31mERROR:\033[0m Unknown fitting mode '%s'",FittingMode_pp.Data());
                }

                if(FittingMode_pL=="Norm"||FittingMode_pL.Contains("Baseline")){
                    fitter->SetSystem(1,*hData_pL,1,CkDec_pL,
                                  FitRegion_pL[0],FitRegion_pL[1],FitRegion_pL[1],FitRegion_pL[1]);
                }
                else if(FittingMode_pL.Contains("Longbaseline")){
                    fitter->SetSystem(1,*hData_pL,1,CkDec_pL,
                                  FitRegion_pL[0],FitRegion_pL[1],FitRegion_pL[2],FitRegion_pL[3]);
                }
                else{
                    printf("\033[1;31mERROR:\033[0m Unknown fitting mode '%s'",FittingMode_pL.Data());
                }

                fitter->SetSeparateBL(0,false);
                fitter->SetSeparateBL(1,false);

                //!p_a
                fitter->SetParameter("pp",DLM_Fitter1::p_a,1.0,0.7,1.3);
                fitter->SetParameter("pLambda",DLM_Fitter1::p_a,1.0,0.7,1.3);

                //!p_b
                if(FittingMode_pp.Contains("Baseline")||FittingMode_pp.Contains("Longbaseline"))
                    fitter->SetParameter("pp",DLM_Fitter1::p_b,0,-2e-3,2e-3);
                else
                    fitter->FixParameter("pp",DLM_Fitter1::p_b,0);
                if(FittingMode_pL.Contains("Baseline")||FittingMode_pL.Contains("Longbaseline"))
                    fitter->SetParameter("pLambda",DLM_Fitter1::p_b,0,-2e-3,2e-3);
                else
                    fitter->FixParameter("pLambda",DLM_Fitter1::p_b,0);


                //!p_c
                if(FittingMode_pp.Contains("Baseline2")||FittingMode_pp.Contains("Longbaseline2"))
                    fitter->SetParameter("pp",DLM_Fitter1::p_c,0,-2e-3,2e-3);
                else
                    fitter->FixParameter("pp",DLM_Fitter1::p_c,0);
                if(FittingMode_pL.Contains("Baseline2")||FittingMode_pL.Contains("Longbaseline2"))
                    fitter->SetParameter("pLambda",DLM_Fitter1::p_c,0,-2e-3,2e-3);
                else
                    fitter->FixParameter("pLambda",DLM_Fitter1::p_c,0);

                fitter->FixParameter("pp",DLM_Fitter1::p_Cl,-1);
                fitter->FixParameter("pLambda",DLM_Fitter1::p_Cl,-1);

                //!p_sor0
                if(SourceType=="Gauss"){

                    fitter->SetParameter("pp",DLM_Fitter1::p_sor0,1.1,0.6,1.6);
                    fitter->SetParameter("pLambda",DLM_Fitter1::p_sor0,1.2,0.8,1.8);
                    if(SameSourceSize){
                        fitter->AddSameParameter("pLambda",DLM_Fitter1::p_sor0,"pp",DLM_Fitter1::p_sor0);
                    }

                    //fitter->FixParameter("pp",DLM_Fitter1::p_sor0,1.41);
                    //fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,1.41);
                }
                else if(SourceType=="CleverLevy_Nolan"||SourceType=="CleverLevy_Single"||SourceType=="CleverLevy_Diff"){
                    fitter->SetParameter("pp",DLM_Fitter1::p_sor0,1.0,0.3,1.7);
                    fitter->SetParameter("pp",DLM_Fitter1::p_sor1,1.6,1.0,2.0);
                    fitter->SetParameter("pLambda",DLM_Fitter1::p_sor0,1.0,0.3,1.7);
                    fitter->SetParameter("pLambda",DLM_Fitter1::p_sor1,1.4,1.0,2.0);
                    if(SameSourceSize){
                        fitter->AddSameParameter("pLambda",DLM_Fitter1::p_sor0,"pp",DLM_Fitter1::p_sor0);
                    }
                    if(SameStability){
                        fitter->AddSameParameter("pLambda",DLM_Fitter1::p_sor1,"pp",DLM_Fitter1::p_sor1);
                    }
                }
                else if(SourceType=="McLevyNolan_Reso"){
                    /*
                    fitter->SetParameter("pp",DLM_Fitter1::p_sor0,0.8,0.2,1.4);
                    fitter->SetParameter("pp",DLM_Fitter1::p_sor1,1.6,1.0,2.0);
                    fitter->SetParameter("pLambda",DLM_Fitter1::p_sor0,0.8,0.2,1.4);
                    fitter->SetParameter("pLambda",DLM_Fitter1::p_sor1,1.4,1.0,2.0);
                    if(SameSourceSize){
                        fitter->AddSameParameter("pLambda",DLM_Fitter1::p_sor0,"pp",DLM_Fitter1::p_sor0);
                    }
                    if(SameStability){
                        fitter->AddSameParameter("pLambda",DLM_Fitter1::p_sor1,"pp",DLM_Fitter1::p_sor1);
                    }
                    */

                    if(DataSample=="pp13TeV_MB_Run2paper"){
                        fitter->FixParameter("pp",DLM_Fitter1::p_sor0,0.82);//at mT 1.27
                        fitter->FixParameter("pp",DLM_Fitter1::p_sor1,1.5);
                        fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,0.71);//at mT 1.44
                        fitter->FixParameter("pLambda",DLM_Fitter1::p_sor1,1.5);
                    }
                    else if(DataSample=="pPb5TeV_Run2paper"||DataSample=="pPb5TeV_CPR_Mar19"){
                        fitter->FixParameter("pp",DLM_Fitter1::p_sor0,1.12);//at mT 1.32
                        fitter->FixParameter("pp",DLM_Fitter1::p_sor1,1.45);
                        fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,1.02);//at mT 1.52
                        fitter->FixParameter("pLambda",DLM_Fitter1::p_sor1,1.45);
                    }
                    else if(DataSample=="pp13TeV_HM_March19"){
                        //fitter->FixParameter("pp",DLM_Fitter1::p_sor0,1.03);//at mT 1.35
                        //fitter->FixParameter("pp",DLM_Fitter1::p_sor1,1.5);
                        //fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,0.92);//at mT 1.55
                        //fitter->FixParameter("pLambda",DLM_Fitter1::p_sor1,1.5);
                        fitter->SetParameter("pp",DLM_Fitter1::p_sor0,1.0,0.3,1.7);
                        fitter->SetParameter("pp",DLM_Fitter1::p_sor1,1.6,1.0,2.0);
                        fitter->SetParameter("pLambda",DLM_Fitter1::p_sor0,1.0,0.3,1.7);
                        fitter->SetParameter("pLambda",DLM_Fitter1::p_sor1,1.4,1.0,2.0);
                    }
                    else printf("WHAT HAPPEND?\n");

                }
                else if(SourceType=="McGauss_Reso"){
                    if(DataSample=="pp13TeV_MB_Run2paper"){
                        fitter->FixParameter("pp",DLM_Fitter1::p_sor0,0.79);//at mT 1.27
                        fitter->FixParameter("pp",DLM_Fitter1::p_sor1,2.0);
                        fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,0.71);//at mT 1.44
                        fitter->FixParameter("pLambda",DLM_Fitter1::p_sor1,2.0);
                    }
                    else if(DataSample=="pPb5TeV_Run2paper"||DataSample=="pPb5TeV_CPR_Mar19"){
                        fitter->FixParameter("pp",DLM_Fitter1::p_sor0,1.01);//at mT 1.32
                        fitter->FixParameter("pp",DLM_Fitter1::p_sor1,2.0);
                        fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,0.95);//at mT 1.52
                        fitter->FixParameter("pLambda",DLM_Fitter1::p_sor1,2.0);
                    }
                    else if(DataSample=="pp13TeV_HM_March19"){
                        fitter->FixParameter("pp",DLM_Fitter1::p_sor0,0.96);//at mT 1.35
                        fitter->FixParameter("pp",DLM_Fitter1::p_sor1,2.0);
                        fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,0.87);//at mT 1.55
                        fitter->FixParameter("pLambda",DLM_Fitter1::p_sor1,2.0);
                    }
                    else printf("WHAT HAPPEND?\n");

                    if(SameSourceSize){
                        //fitter->AddSameParameter("pLambda",DLM_Fitter1::p_sor0,"pp",DLM_Fitter1::p_sor0);
                    }
                }
                else{
                    printf("\033[1;31mERROR:\033[0m '%s' does not exist",SourceType.Data());
                }


//fitter->AddSameParameter("pp",DLM_Fitter1::p_sor0,"pLambda",DLM_Fitter1::p_sor0);
//fitter->AddSameParameter("pp",DLM_Fitter1::p_sor1,"pLambda",DLM_Fitter1::p_sor1);
//fitter->AddSameParameter("pp",DLM_Fitter1::p_a,"pLambda",DLM_Fitter1::p_a);
//fitter->AddSameParameter("pp",DLM_Fitter1::p_b,"pLambda",DLM_Fitter1::p_b);
//fitter->AddSameParameter("pp",DLM_Fitter1::p_c,"pLambda",DLM_Fitter1::p_c);

fitter->FixParameter("pp",DLM_Fitter1::p_sor0,1.05);
fitter->FixParameter("pp",DLM_Fitter1::p_sor1,1.6);

fitter->AddSameParameter("pLambda",DLM_Fitter1::p_sor0,"pp",DLM_Fitter1::p_sor0);
fitter->AddSameParameter("pLambda",DLM_Fitter1::p_sor1,"pp",DLM_Fitter1::p_sor1);
fitter->AddSameParameter("pLambda",DLM_Fitter1::p_a,"pp",DLM_Fitter1::p_a);
fitter->AddSameParameter("pLambda",DLM_Fitter1::p_b,"pp",DLM_Fitter1::p_b);
fitter->AddSameParameter("pLambda",DLM_Fitter1::p_c,"pp",DLM_Fitter1::p_c);



                //AB_pp.KillTheCat();
                //AB_pL.KillTheCat();

                CkDec_pp.Update();
                CkDec_pL.Update();
                CkDec_pSigma0.Update();
                CkDec_pXim.Update();

                fitter->GoBabyGo(false);

                TGraph FitResult_pp;
                FitResult_pp.SetName(TString::Format("pp_%s_%s_%s",pp_Pot[uPot_pp].Data(),pL_Pot[uPot_pL].Data(),pXim_Pot[uPot_pXim].Data()));
                fitter->GetFitGraph(0, FitResult_pp);

                TGraph FitResult_pL;
                FitResult_pL.SetName(TString::Format("pL_%s_%s_%s",pp_Pot[uPot_pp].Data(),pL_Pot[uPot_pL].Data(),pXim_Pot[uPot_pXim].Data()));
                fitter->GetFitGraph(1, FitResult_pL);
                //TString Description = TString::Format("%s_%s_%s",pp_Pot[uPot_pp].Data(),pL_Pot[uPot_pL].Data(),pXim_Pot[uPot_pXim].Data());
                OutputFile->cd();
                FitResult_pp.Write();
                FitResult_pL.Write();

                printf("χ2/ndf = %.2f/%i = %.2f\n",fitter->GetChi2(),fitter->GetNdf(),fitter->GetChi2Ndf());
                printf("R(pp) = %.3f +/- %.3f\n",fitter->GetParameter("pp",DLM_Fitter1::p_sor0),fitter->GetParError("pp",DLM_Fitter1::p_sor0));
                printf(" R(pΛ)= %.3f +/- %.3f\n",fitter->GetParameter("pLambda",DLM_Fitter1::p_sor0),fitter->GetParError("pLambda",DLM_Fitter1::p_sor0));
                if(SourceType=="CleverLevy_Nolan"||SourceType=="CleverLevy_Single"||SourceType=="CleverLevy_Diff"||
                   SourceType=="McLevyNolan_Reso"||SourceType=="McGauss_Reso"){
                    printf("α(pp) = %.3f +/- %.3f\n",fitter->GetParameter("pp",DLM_Fitter1::p_sor1),fitter->GetParError("pp",DLM_Fitter1::p_sor1));
                    printf(" α(pΛ)= %.3f +/- %.3f\n",fitter->GetParameter("pLambda",DLM_Fitter1::p_sor1),fitter->GetParError("pLambda",DLM_Fitter1::p_sor1));
                }
                printf("Info on the baseline:\n");
                printf("a(pp) = %.3f +/- %.3f\n",fitter->GetParameter("pp",DLM_Fitter1::p_a),fitter->GetParError("pp",DLM_Fitter1::p_a));
                printf(" a(pΛ)= %.3f +/- %.3f\n",fitter->GetParameter("pLambda",DLM_Fitter1::p_a),fitter->GetParError("pLambda",DLM_Fitter1::p_a));
                printf("a(pp) = %.3fe-3 +/- %.3f\n",fitter->GetParameter("pp",DLM_Fitter1::p_b)*1e3,fitter->GetParError("pp",DLM_Fitter1::p_b)*1e3);
                printf(" a(pΛ)= %.3fe-3 +/- %.3f\n",fitter->GetParameter("pLambda",DLM_Fitter1::p_b)*1e3,fitter->GetParError("pLambda",DLM_Fitter1::p_b)*1e3);
                printf("a(pp) = %.3fe-6 +/- %.3f\n",fitter->GetParameter("pp",DLM_Fitter1::p_c)*1e6,fitter->GetParError("pp",DLM_Fitter1::p_c)*1e6);
                printf(" a(pΛ)= %.3fe-6 +/- %.3f\n",fitter->GetParameter("pLambda",DLM_Fitter1::p_c)*1e6,fitter->GetParError("pLambda",DLM_Fitter1::p_c)*1e6);
                //printf("AB_pp.EvalCorrFun(25)=%.4f\n",AB_pp.EvalCorrFun(25));
                //printf("AB_pp.EvaluateTheSource(25,2,0)=%.4f\n",AB_pp.EvaluateTheSource(25,2,0));
                //printf("pars: %.3f, %.3f\n",AB_pp.AnaSourcePar->GetParameter(0),AB_pp.AnaSourcePar->GetParameter(1));

                delete Ck_pXim;
                delete fitter;
            }
            delete Ck_pL;
        }
        delete Ck_pp;
    }

    delete Ck_pSigma0;
    delete [] MomBins_pp;
    delete [] MomBins_pL;
    delete [] FitRegion_pp;
    delete [] FitRegion_pL;
    //delete hDataClever_pL;
    delete OutputFile;
}

//For feed-down we will use only pSigma0 and pXi, where the former we model as Oli,
//the latter we model with HALQCD, by fixing the source and radius to a Gauss of size 1.18 (-20% for the systematic variation)
//Fitting mode:
    //"Norm": only norm
    //"Baseline": baseline but only in the femto region
    //"Longbaseline": baseline extending to higher k
//if SourceScale||SourceStability==0 => we take the values we have obtained from the pp fit (hard-coded at the moment)
//SourceScale == -1, fit within 10% of the expected value (based on mT or whatever)
//SourceScale == -11, completely free fit
//pL_Pot: "LO",  "LO_Coupled_S", "NLO", "NLO_sp",  "NLO_Coupled_S",    "Usmani"
//in VARIATIONS:
//[0] = femto range
//[1] = fit range
//[2] = fractions for p
//[3] = fractions for L
//[4] = kc parameter
//[5] = which mT bin
//[6] = in case of Reso source: smearing of the mass of the resonance (on or off)//NOT ACTIVE YET
//[7] = in case of Reso source: smearing of the momentum of the resonance (fraction, from 0 to 1)//NOT ACTIVE YET
void Fit_pL(DLM_CommonAnaFunctions& AnalysisObject, const TString& OutputFolder, const TString& OutputFileName,
            const TString& DataSample, const TString& DataVar, const TString& SourceType, const double& SourceScale, const double& SourceStability,
            const TString& pL_Pot,const TString& FittingMode_pL, const int* VARIATIONS, const long& UniqueID=-1, const TString ntFileName="", const TString ntName=""){

    //printf("VARIATIONS[2]=%i\n",VARIATIONS[2]);


    //the residuals that are not fitted, are assumed to all have a Gaussian source of some size;
    double ResidualSourceSize=0;
    if(DataSample=="pp13TeV_MB_Run2paper") ResidualSourceSize=1.1;
    else if(DataSample=="pPb5TeV_Run2paper"||DataSample=="pPb5TeV_CPR_Mar19") ResidualSourceSize=1.4;
    else if(DataSample=="pp13TeV_HM_March19") ResidualSourceSize=1.4;
    else printf("\033[1;31mERROR:\033[0m The data sample '%s' does not exist\n",DataSample.Data());

    double* MomBins_pL = NULL;
    double* FitRegion_pL = NULL;
    unsigned NumMomBins_pL;
    AnalysisObject.SetUpBinning_pL(DataSample,NumMomBins_pL,MomBins_pL,FitRegion_pL,VARIATIONS[0],VARIATIONS[1]);
    //if(pL_Pot=="LO"||pL_Pot=="NLO"||pL_Pot=="Usmani"){
    //    FitRegion_pL[1] = 192;
    //    FitRegion_pL[2] = 312;
    //    FitRegion_pL[3] = 372;
    //}

    DLM_Ck* Ck_pSigma0 = new DLM_Ck(1,0,NumMomBins_pL,MomBins_pL,Lednicky_gauss_Sigma0);
    Ck_pSigma0->SetSourcePar(0,ResidualSourceSize);

    TH2F* hResolution_pL = AnalysisObject.GetResolutionMatrix(DataSample,"pLambda");

    TH2F* hResidual_pL_pSigma0 = AnalysisObject.GetResidualMatrix("pLambda","pSigma0");
    TH2F* hResidual_pL_pXim = AnalysisObject.GetResidualMatrix("pLambda","pXim");

    TH1F* hData_pL = AnalysisObject.GetAliceExpCorrFun(DataSample,"pLambda",DataVar,2,false,VARIATIONS[5]);

    double lam_pL[5];
    double lam_pXim[5];

    AnalysisObject.SetUpLambdaPars_pL(DataSample,VARIATIONS[2],VARIATIONS[3],lam_pL);
    AnalysisObject.SetUpLambdaPars_pXim(DataSample,VARIATIONS[2],0,lam_pXim);

    TString SourceDescription = SourceType;
    if(SourceScale<0) {SourceDescription+="ScaFree";}
    else if(SourceScale==0) {SourceDescription+="ScaDefault";}
    else {SourceDescription+=TString::Format("Sca%.2f",SourceScale);}
    if(SourceStability<0||SourceStability>2) {SourceDescription+="StaFree";}
    else if(SourceStability==0) {SourceDescription+="StaDefault";}
    else {SourceDescription+=TString::Format("Sta%.2f",SourceStability);}
    TString OutFileName = OutputFileName;
    if(OutFileName==""){
        OutFileName=TString::Format("AUTO%li_%s_%s_%s_%s_%i_%i_%i_%i_%i_mT%i.root",UniqueID,
                                                               DataSample.Data(),SourceDescription.Data(),
                                                               pL_Pot.Data(),FittingMode_pL.Data(),
                                                               VARIATIONS[0],VARIATIONS[1],VARIATIONS[2],VARIATIONS[3],VARIATIONS[4],VARIATIONS[5]);
    }
    TFile* OutputFile = new TFile(OutputFolder+OutFileName,"recreate");
    hData_pL->Write();

    CATS AB_pL;
    DLM_Ck* Ck_pL;
    if(pL_Pot.Contains("Lednicky_")){
        Ck_pL = AnalysisObject.SetUpLednicky_pL(NumMomBins_pL,MomBins_pL,pL_Pot);
        Ck_pL->SetSourcePar(0,SourceScale);
    }
    else{
        AB_pL.SetMomBins(NumMomBins_pL,MomBins_pL);
        AnalysisObject.SetUpCats_pL(AB_pL,pL_Pot,SourceType,0,3);
        AB_pL.SetNotifications(CATS::nWarning);
        AB_pL.KillTheCat();
        Ck_pL = new DLM_Ck(AB_pL.GetNumSourcePars(),0,AB_pL);
        Ck_pL->SetSourcePar(0,SourceScale);
        if(AB_pL.GetNumSourcePars()>1){
            Ck_pL->SetSourcePar(0,SourceStability>=1?SourceStability:1.5);
        }
    }

    CATS AB_pXim;
    //same binning as pL, as we only use pXim as feed-down
    AB_pXim.SetMomBins(NumMomBins_pL,MomBins_pL);
    AnalysisObject.SetUpCats_pXim(AB_pXim,"pXim_HALQCD1","Gauss");
    AB_pXim.SetAnaSource(0,ResidualSourceSize);
    AB_pXim.SetNotifications(CATS::nWarning);
    AB_pXim.KillTheCat();
    DLM_Ck* Ck_pXim = new DLM_Ck(AB_pXim.GetNumSourcePars(),0,AB_pXim);
    Ck_pL->Update();
    Ck_pSigma0->Update();
    Ck_pXim->Update();

    DLM_CkDecomposition CkDec_pL("pLambda",4,*Ck_pL,hResolution_pL);
    DLM_CkDecomposition CkDec_pSigma0("pSigma0",0,*Ck_pSigma0,NULL);
    DLM_CkDecomposition CkDec_pXim("pXim",2,*Ck_pXim,NULL);

    CkDec_pL.AddContribution(0,lam_pL[1],DLM_CkDecomposition::cFeedDown,&CkDec_pSigma0,hResidual_pL_pSigma0);
    CkDec_pL.AddContribution(1,lam_pL[2],DLM_CkDecomposition::cFeedDown,&CkDec_pXim,hResidual_pL_pXim);
    CkDec_pL.AddContribution(2,lam_pL[3],DLM_CkDecomposition::cFeedDown);
    CkDec_pL.AddContribution(3,lam_pL[4],DLM_CkDecomposition::cFake);//0.03

    //for Xim we simplify a bit and take ALL feed-down as flat
    CkDec_pXim.AddContribution(0,lam_pXim[1]+lam_pXim[2]+lam_pXim[3],DLM_CkDecomposition::cFeedDown);
    CkDec_pXim.AddContribution(1,lam_pXim[4],DLM_CkDecomposition::cFake);

    DLM_Fitter1* fitter = new DLM_Fitter1(1);
    if(FittingMode_pL.Contains("Norm")||FittingMode_pL.Contains("Baseline")){
        fitter->SetSystem(0,*hData_pL,1,CkDec_pL,
                    FitRegion_pL[0],FitRegion_pL[1],FitRegion_pL[1],FitRegion_pL[1]);
    }
    else if(FittingMode_pL.Contains("Longbaseline")||FittingMode_pL.Contains("Spline")){
        fitter->SetSystem(0,*hData_pL,1,CkDec_pL,
                        FitRegion_pL[0],FitRegion_pL[1],FitRegion_pL[2],FitRegion_pL[3]);
    }
    else{
        printf("\033[1;31mERROR:\033[0m Unknown fitting mode '%s'",FittingMode_pL.Data());
    }

    if(FittingMode_pL.Contains("_prefit")){
        fitter->SetSeparateBL(0,true);
        fitter->FixParameter("pLambda",DLM_Fitter1::p_Cl,-1);
    }
    else{
        fitter->SetSeparateBL(0,false);
        //fitter->FixParameter("pLambda",DLM_Fitter1::p_Cl,-1);
        //fitter->SetParameter("pLambda",DLM_Fitter1::p_Cl,-1,-1,-0.95);

        if(pL_Pot.Contains("Lednicky_")){
            fitter->FixParameter("pLambda",DLM_Fitter1::p_Cl,1);
        }
        else if(VARIATIONS[4]<0){
            fitter->FixParameter("pLambda",DLM_Fitter1::p_Cl,-1);
        }
        else{
            //fitter->SetParameter("pLambda",DLM_Fitter1::p_Cl,1,0.9,1.1);
            if(VARIATIONS[4]>0&&FittingMode_pL!="Norm") fitter->FixParameter("pLambda",DLM_Fitter1::p_kc,VARIATIONS[4]);
            else fitter->SetParameter("pLambda",DLM_Fitter1::p_kc,2000,400,10000);
            //fitter->SetParameter("pLambda",DLM_Fitter1::p_kc,2000,400,10000);
        }

    }


    printf("FittingMode_pL=%s\n",FittingMode_pL.Data());

    //!p_a
    fitter->SetParameter("pLambda",DLM_Fitter1::p_a,1.0,0.7,1.3);
    if(FittingMode_pL.Contains("Spline")){
        fitter->FixParameter("pLambda",DLM_Fitter1::p_a,1.0);
    }

    //!p_b
    if(FittingMode_pL.Contains("Baseline")||FittingMode_pL.Contains("Longbaseline"))
        fitter->SetParameter("pLambda",DLM_Fitter1::p_b,0,-2e-3,2e-3);
    else
        fitter->FixParameter("pLambda",DLM_Fitter1::p_b,0);

    //!p_c
    if(FittingMode_pL.Contains("Baseline2")||FittingMode_pL.Contains("Longbaseline2"))
        fitter->SetParameter("pLambda",DLM_Fitter1::p_c,0,-2e-4,2e-4);
    else
        fitter->FixParameter("pLambda",DLM_Fitter1::p_c,0);

    if(FittingMode_pL.Contains("Spline")){
        int NumKnots=0;
        if(FittingMode_pL.Contains("_2")) NumKnots=2;
        if(FittingMode_pL.Contains("_3")) NumKnots=3;
        if(FittingMode_pL.Contains("_4")) NumKnots=4;
        if(FittingMode_pL.Contains("_5")) NumKnots=5;
        if(FittingMode_pL.Contains("_6")) NumKnots=6;
        if(FittingMode_pL.Contains("_7")) NumKnots=7;
        if(FittingMode_pL.Contains("_8")) NumKnots=8;
        if(FittingMode_pL.Contains("_9")) NumKnots=9;
        if(FittingMode_pL.Contains("_10")) NumKnots=10;
        if(FittingMode_pL.Contains("_11")) NumKnots=11;
        if(FittingMode_pL.Contains("_Tune1")) NumKnots=5;
        fitter->FixParameter("pLambda",DLM_Fitter1::p_spline,NumKnots);
        fitter->FixParameter("pLambda",DLM_Fitter1::p_spline+1,0);
        fitter->SetParameter("pLambda",DLM_Fitter1::p_spline+2,0);

        double* Nodes_x = new double [NumKnots];
        Nodes_x[0] = FitRegion_pL[0];
        Nodes_x[NumKnots-1] = FitRegion_pL[3];
        if(FittingMode_pL.Contains("_Tune1")){
            Nodes_x[0] = 0;
            Nodes_x[1] = 90;
            Nodes_x[2] = 180;
            Nodes_x[3] = 360;
            Nodes_x[4] = 540;
            fitter->FixParameter("pLambda",DLM_Fitter1::p_spline+2,0.00022);
        }
        else{
            double NodeLength = (Nodes_x[NumKnots-1]-Nodes_x[0])/double(NumKnots-1);
            //printf("NodeLength=%f\n",NodeLength);
            for(int iKnot=1; iKnot<NumKnots-1; iKnot++){
                Nodes_x[iKnot] = Nodes_x[iKnot-1]+NodeLength;
                //printf("Nodes_x[%i]=%f\n",iKnot,Nodes_x[iKnot]);
            }
        }

        for(int iKnot=0; iKnot<NumKnots; iKnot++){
            fitter->FixParameter("pLambda",DLM_Fitter1::p_spline+3+iKnot,Nodes_x[iKnot]);//x
            fitter->SetParameter("pLambda",DLM_Fitter1::p_spline+3+NumKnots+iKnot,1,0,5);//y
        }
        delete [] Nodes_x;
    }

    //!p_sor0
    if(SourceType=="Gauss"){
        if(SourceScale==0){fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,1.2);}
        else if(SourceScale<0) {fitter->SetParameter("pLambda",DLM_Fitter1::p_sor0,1.2,0.5,2.0);}
        else {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,SourceScale);}
    }
    else if(SourceType=="CleverLevy_Nolan"||SourceType=="CleverLevy_Single"||SourceType=="CleverLevy_Diff"){
        fitter->SetParameter("pLambda",DLM_Fitter1::p_sor0,1.0,0.3,1.7);
        fitter->SetParameter("pLambda",DLM_Fitter1::p_sor1,1.4,1.0,2.0);
    }
    else if(SourceType=="McLevyNolan_Reso"){
        if(DataSample=="pp13TeV_MB_Run2paper"){
            if(SourceScale==0) {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,0.71);}//at mT 1.44
            else if(SourceScale==-1) {fitter->SetParameter("pLambda",DLM_Fitter1::p_sor0,0.71,0.71/2.,0.71*2.);}
            else if(SourceScale==-11) {fitter->SetParameter("pLambda",DLM_Fitter1::p_sor0,0.71,0.71/2.,0.71*2.);}
            else {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,SourceScale);}

            if(SourceStability==0) {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor1,1.5);}
            else if(SourceStability==-1) {fitter->SetParameter("pLambda",DLM_Fitter1::p_sor1,1.5,1.5*0.90,1.5*1.10);}
            else if(SourceStability==-11) {fitter->SetParameter("pLambda",DLM_Fitter1::p_sor1,1.5,1.0,2.0);}
            else {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor1,SourceStability);}
        }
        else if(DataSample=="pPb5TeV_Run2paper"||DataSample=="pPb5TeV_CPR_Mar19"){
            if(SourceScale==0) {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,1.02);}//at mT 1.44
            else if(SourceScale==-1) {fitter->SetParameter("pLambda",DLM_Fitter1::p_sor0,1.02,1.02*0.90,1.02*1.10);}
            else if(SourceScale==-11) {fitter->SetParameter("pLambda",DLM_Fitter1::p_sor0,1.02,1.02/2.,1.02*2.);}
            else {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,SourceScale);}

            if(SourceStability==0) {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor1,1.45);}
            else if(SourceStability==-1) {fitter->SetParameter("pLambda",DLM_Fitter1::p_sor1,1.45,1.45*0.90,1.45*1.10);}
            else if(SourceStability==-11) {fitter->SetParameter("pLambda",DLM_Fitter1::p_sor1,1.5,1.0,2.0);}
            else {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor1,SourceStability);}
        }
        else if(DataSample=="pp13TeV_HM_March19"){
            if(SourceScale==0) {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,0.92);}//at mT 1.55
            else if(SourceScale==-1) {fitter->SetParameter("pLambda",DLM_Fitter1::p_sor0,0.92,0.92*0.90,0.92*1.10);}
            else if(SourceScale==-11) {fitter->SetParameter("pLambda",DLM_Fitter1::p_sor0,0.92,0.92/2.,0.92*2.);}
            else {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,SourceScale);}

            if(SourceStability==0) {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor1,1.5);}
            else if(SourceStability==-1) {fitter->SetParameter("pLambda",DLM_Fitter1::p_sor1,1.5,1.5*0.90,1.5*1.10);}
            else if(SourceStability==-11) {fitter->SetParameter("pLambda",DLM_Fitter1::p_sor1,1.5,1.0,2.0);}
            else {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor1,SourceStability);}
        }
        else printf("WHAT HAPPEND?\n");
    }
    else if(SourceType=="McGauss_Reso"){
        if(DataSample=="pp13TeV_MB_Run2paper"){
            if(SourceScale==0) {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,0.71);}//at mT 1.44
            else if(SourceScale==-1) {fitter->SetParameter("pLambda",DLM_Fitter1::p_sor0,0.71,0.71*0.90,0.71*1.10);}
            else if(SourceScale==-11) {fitter->SetParameter("pLambda",DLM_Fitter1::p_sor0,0.71,0.71/2.,0.71*2.);}
            else {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,SourceScale);}
            fitter->FixParameter("pLambda",DLM_Fitter1::p_sor1,2.0);
        }
        else if(DataSample=="pPb5TeV_Run2paper"||DataSample=="pPb5TeV_CPR_Mar19"){
            if(SourceScale==0) {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,0.95);}//at mT 1.52
            else if(SourceScale==-1) {fitter->SetParameter("pLambda",DLM_Fitter1::p_sor0,0.95,0.95*0.90,0.95*1.10);}
            else if(SourceScale==-11) {fitter->SetParameter("pLambda",DLM_Fitter1::p_sor0,0.95,0.95/2.,0.95*2.);}
            else {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,SourceScale);}
            fitter->FixParameter("pLambda",DLM_Fitter1::p_sor1,2.0);
        }
        else if(DataSample=="pp13TeV_HM_March19"){
            //if(SourceScale==0) {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,0.87);}//at mT 1.55
            //else if(SourceScale==-1) {fitter->SetParameter("pLambda",DLM_Fitter1::p_sor0,0.87,0.87*0.90,0.87*1.10);}
            //else if(SourceScale==-11) {fitter->SetParameter("pLambda",DLM_Fitter1::p_sor0,0.87,0.87/2.,0.87*2.);}
            //else {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,SourceScale);}

            //these are the value shown at the PF previews
            //if(SourceScale==0) {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,0.918078);}//at mT 1.55
            //else if(SourceScale==-1) {fitter->SetParameter("pLambda",DLM_Fitter1::p_sor0,0.918078,0.90218,0.933344);}
            if(SourceScale==0) {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,0.918);}//at mT 1.55
            else if(SourceScale==-1) {fitter->SetParameter("pLambda",DLM_Fitter1::p_sor0,0.918,0.871,0.965);}//3 sigma
            else if(SourceScale==-11) {fitter->SetParameter("pLambda",DLM_Fitter1::p_sor0,0.87,0.87/2.,0.87*2.);}
            else {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,SourceScale);}

            fitter->FixParameter("pLambda",DLM_Fitter1::p_sor1,2.0);
        }
        else printf("WHAT HAPPEND?\n");
    }
    else{
        printf("\033[1;31mERROR:\033[0m '%s' does not exist",SourceType.Data());
    }

    if(pL_Pot.Contains("Lednicky_")){
        fitter->FixParameter("pLambda",DLM_Fitter1::p_pot0,Ck_pL->GetPotPar(0));
        fitter->FixParameter("pLambda",DLM_Fitter1::p_pot1,Ck_pL->GetPotPar(1));
        fitter->FixParameter("pLambda",DLM_Fitter1::p_pot2,Ck_pL->GetPotPar(2));
        fitter->FixParameter("pLambda",DLM_Fitter1::p_pot3,Ck_pL->GetPotPar(3));
        //fitter->SetFullCkForBaseline(0,true);
        //fitter->FixParameter("pLambda",DLM_Fitter1::p_Cl,-1);
    }

    CkDec_pL.Update();
    CkDec_pSigma0.Update();
    CkDec_pXim.Update();
    fitter->GoBabyGo(false);
    TGraph FitResult_pL;
    FitResult_pL.SetName(TString::Format("FitResult_pL"));
    fitter->GetFitGraph(0, FitResult_pL);
    //TString Description = TString::Format("%s_%s_%s",pp_Pot[uPot_pp].Data(),pL_Pot[uPot_pL].Data(),pXim_Pot[uPot_pXim].Data());
    OutputFile->cd();
    FitResult_pL.Write();

    TGraph FitBaseline_pL;
    FitBaseline_pL.SetName(TString::Format("FitBaseline_pL"));
    fitter->GetMultBaselineGraph(0, FitBaseline_pL);
    //TString Description = TString::Format("%s_%s_%s",pp_Pot[uPot_pL].Data(),pL_Pot[uPot_pL].Data(),pXim_Pot[uPot_pXim].Data());
    OutputFile->cd();
    FitBaseline_pL.Write();

    double CHI2_312=0;
    unsigned NDF_312=0;
    double xVal;
    double yVal;
    unsigned hBin;
    for(unsigned uPoint=0; uPoint<FitResult_pL.GetN(); uPoint++){
        FitResult_pL.GetPoint(uPoint,xVal,yVal);
        if(xVal>312) break;
        hBin = hData_pL->FindBin(xVal);
        CHI2_312 += pow( (hData_pL->GetBinContent(hBin)-yVal)/(hData_pL->GetBinError(hBin)) ,2.);
        NDF_312++;
    }

    printf("χ2/ndf = %.2f/%i = %.2f\n",fitter->GetChi2(),fitter->GetNdf(),fitter->GetChi2Ndf());
    printf("χ2/ndf (312) = %.2f/%i = %.2f\n",CHI2_312,NDF_312,CHI2_312/double(NDF_312));
    printf("R(pΛ) = %.3f +/- %.3f\n",fitter->GetParameter("pLambda",DLM_Fitter1::p_sor0),fitter->GetParError("pLambda",DLM_Fitter1::p_sor0));
    if(AB_pL.GetNumSourcePars()>1){
        printf("α(pΛ) = %.3f +/- %.3f\n",fitter->GetParameter("pLambda",DLM_Fitter1::p_sor1),fitter->GetParError("pLambda",DLM_Fitter1::p_sor1));
    }
    printf("a = %.2e +/- %.2e\n",fitter->GetParameter("pLambda",DLM_Fitter1::p_a),fitter->GetParError("pLambda",DLM_Fitter1::p_a));
    printf("b = %.2e +/- %.2e\n",fitter->GetParameter("pLambda",DLM_Fitter1::p_b),fitter->GetParError("pLambda",DLM_Fitter1::p_b));
    printf("c = %.2e +/- %.2e\n",fitter->GetParameter("pLambda",DLM_Fitter1::p_c),fitter->GetParError("pLambda",DLM_Fitter1::p_c));
    printf("Cl = %.3f +/- %.3f\n",fitter->GetParameter("pLambda",DLM_Fitter1::p_Cl),fitter->GetParError("pLambda",DLM_Fitter1::p_Cl));
    printf("kc = %.3f +/- %.3f\n",fitter->GetParameter("pLambda",DLM_Fitter1::p_kc),fitter->GetParError("pLambda",DLM_Fitter1::p_kc));
    printf("mT = %i\n",VARIATIONS[5]);

    TGraph gChi2Ndf;
    gChi2Ndf.SetName("gChi2Ndf");
    gChi2Ndf.Set(1);
    gChi2Ndf.SetPoint(0,fitter->GetChi2(),fitter->GetNdf());
    gChi2Ndf.SetMarkerStyle(2);
    gChi2Ndf.SetMarkerSize(4);
    gChi2Ndf.Write();

    TF1* fBaseline = new TF1("fBaseline","[0]+[1]*x+[2]*x*x",0,600);
    fBaseline->FixParameter(0,fitter->GetParameter("pLambda",DLM_Fitter1::p_a));
    fBaseline->FixParameter(1,fitter->GetParameter("pLambda",DLM_Fitter1::p_b));
    fBaseline->FixParameter(2,fitter->GetParameter("pLambda",DLM_Fitter1::p_c));
    fBaseline->SetLineColor(kCyan+4);
    fBaseline->SetLineWidth(3);
    fBaseline->SetLineStyle(7);
    fBaseline->Write();

    TGraphErrors gScaleStability;
    gScaleStability.SetName("gScaleStability");
    gScaleStability.Set(1);
    double SCALE = fitter->GetParameter("pLambda",DLM_Fitter1::p_sor0);
    double dSCALE = fitter->GetParError("pLambda",DLM_Fitter1::p_sor0);
    double STABILITY = AB_pL.GetNumSourcePars()>1?fitter->GetParameter("pLambda",DLM_Fitter1::p_sor1):2;
    double dSTABILITY = AB_pL.GetNumSourcePars()>1?fitter->GetParError("pLambda",DLM_Fitter1::p_sor1):0;
    gScaleStability.SetPoint(0,SCALE,STABILITY);
    gScaleStability.SetPointError(0,dSCALE,dSTABILITY);
    gScaleStability.SetMarkerStyle(20);
    gScaleStability.SetMarkerSize(1);
    gScaleStability.SetLineWidth(2);
    gScaleStability.Write();

    //TFile* OutputFileNT = new TFile(OutputFolder+TString::Format("NT_%s_%s_%s_%s_%i_%i_%i_%i_%i.root",
    //                                                           DataSample.Data(),SourceDescription.Data(),
    //                                                           pL_Pot.Data(),FittingMode_pL.Data(),
    //                                                           VARIATIONS[0],VARIATIONS[1],VARIATIONS[2],VARIATIONS[3],VARIATIONS[4]),"recreate");


    TFile* NtFile = NULL;
    TNtuple* ntResult = NULL;
    if(ntFileName!=""&&ntName!=""){
        bool FileIsOpen=false;
        bool FileDoesNotExist=false;
        double WaitedInSec=0;
        do{
            if(FileIsOpen) {usleep(1000); WaitedInSec+=1e-3;}//sleep for a millisecond
            if(WaitedInSec>60){
                printf("\033[1;31mERROR: \033[0mCannot access the file %s\n", ntFileName.Data());
                printf("\033[1;31m ABORTING THE EXECUTABLE!\033[0m\n");
                abort();
            }
            int fd = open(ntFileName.Data(), O_RDONLY);
            if (fd < 0) {
                //perror("open");
                //file does not exist
                FileIsOpen = false;
                FileDoesNotExist = true;
            }
            if (fcntl(fd, F_SETLEASE, F_WRLCK) && EAGAIN == errno) {
                FileIsOpen = true;
            }
            else {
                fcntl(fd, F_SETLEASE, F_UNLCK);
                FileIsOpen = false;
            }
            close(fd);
        }
        while(FileIsOpen);
        if(FileDoesNotExist) NtFile = new TFile(ntFileName,"RECREATE");
        else NtFile = new TFile(ntFileName,"UPDATE");
        ntResult = (TNtuple*)NtFile->Get(ntName);
        if(!ntResult){
            //printf("\033[1;31mERROR:\033[0m For whatever reason the ");
            ntResult = new TNtuple("ntResult", "ntResult",
            "IterM:Iter:Config:Data:SourceType:SourceScale:SourceStability:Potential:Baseline:FemRan:FitRan:pFrac:LamFrac:kcVar:mTbin:FemtoMin:FemtoMax:BlMin:BlMax:"
            "p_a:e_a:p_b:e_b:p_c:e_c:p_Cl:e_Cl:p_kc:e_kc:p_sor0:e_sor0:p_sor1:e_sor1:chi2:ndf:chi2_312:ndf_312");
        }
    }

    if(NtFile&&ntResult){
        Float_t buffer[37];
        long Iter = UniqueID/1000;
        long Config = UniqueID%1000;
        buffer[0] = Iter/long(1000000);
        buffer[1] = Iter%long(1000000);
        buffer[2] = UniqueID%1000;
        TString Temp = DataVar;
        Temp.ReplaceAll("_","");
        buffer[3] = Temp.Atoi();
        if(SourceType=="Gauss") buffer[4] = 0;
        if(SourceType=="McGauss_Reso") buffer[4] = 1;
        if(SourceType=="McLevyNolan_Reso") buffer[4] = 2;

        buffer[5] = SourceScale;
        buffer[6] = SourceStability;

        if(pL_Pot=="LO") buffer[7] = 0;
        if(pL_Pot=="NLO") buffer[7] = 10;
        if(pL_Pot=="LO_Coupled_S") buffer[7] = 1;
        if(pL_Pot=="NLO_Coupled_S") buffer[7] = 11;
        if(pL_Pot=="NLO_sp") buffer[7] = 12;
        if(pL_Pot=="Usmani") buffer[7] = 100;
        if(pL_Pot=="Lednicky_ND") buffer[7] = 1001;
        if(pL_Pot=="Lednicky_NF") buffer[7] = 1002;
        if(pL_Pot=="Lednicky_NSC89") buffer[7] = 1003;
        if(pL_Pot=="Lednicky_NSC97a") buffer[7] = 1004;
        if(pL_Pot=="Lednicky_NSC97b") buffer[7] = 1005;
        if(pL_Pot=="Lednicky_NSC97c") buffer[7] = 1006;
        if(pL_Pot=="Lednicky_NSC97d") buffer[7] = 1007;
        if(pL_Pot=="Lednicky_NSC97e") buffer[7] = 1008;
        if(pL_Pot=="Lednicky_NSC97f") buffer[7] = 1009;
        if(pL_Pot=="Lednicky_ESC08") buffer[7] = 1010;
        if(pL_Pot=="Lednicky_XeftLO") buffer[7] = 1011;
        if(pL_Pot=="Lednicky_XeftNLO") buffer[7] = 1012;
        if(pL_Pot=="Lednicky_JulichA") buffer[7] = 1013;
        if(pL_Pot=="Lednicky_JulichJ04") buffer[7] = 1014;
        if(pL_Pot=="Lednicky_JulichJ04c") buffer[7] = 1015;

        if(FittingMode_pL=="Norm") buffer[8] = 0;
        if(FittingMode_pL=="Baseline") buffer[8] = 1;
        if(FittingMode_pL=="Baseline2") buffer[8] = 2;
        if(FittingMode_pL=="Longbaseline") buffer[8] = 11;
        if(FittingMode_pL=="Longbaseline2") buffer[8] = 12;
        if(FittingMode_pL=="Norm_prefit") buffer[8] = 100;
        if(FittingMode_pL=="Baseline_prefit") buffer[8] = 101;
        if(FittingMode_pL=="Baseline2_prefit") buffer[8] = 102;
        if(FittingMode_pL=="Longbaseline_prefit") buffer[8] = 111;
        if(FittingMode_pL=="Longbaseline2_prefit") buffer[8] = 112;
        if(FittingMode_pL=="Spline3_3") buffer[8] = 403;
        if(FittingMode_pL=="Spline3_4") buffer[8] = 404;
        if(FittingMode_pL=="Spline3_5") buffer[8] = 405;
        if(FittingMode_pL=="Spline3_6") buffer[8] = 406;
        if(FittingMode_pL=="Spline3_7") buffer[8] = 407;
        if(FittingMode_pL=="Spline3_8") buffer[8] = 408;
        if(FittingMode_pL=="Spline3_9") buffer[8] = 409;
        if(FittingMode_pL=="Spline3_10") buffer[8] = 410;
        if(FittingMode_pL=="Spline3_11") buffer[8] = 411;
        if(FittingMode_pL=="Spline3_Tune1") buffer[8] = 451;

        buffer[9] = VARIATIONS[0];
        buffer[10] = VARIATIONS[1];
        buffer[11] = VARIATIONS[2];
        buffer[12] = VARIATIONS[3];
        buffer[13] = VARIATIONS[4];
        buffer[14] = VARIATIONS[5];

        buffer[15] = FitRegion_pL[0];
        buffer[16] = FitRegion_pL[1];
        buffer[17] = FitRegion_pL[2];
        buffer[18] = FitRegion_pL[3];

        buffer[19] = fitter->GetParameter("pLambda",DLM_Fitter1::p_a);
        buffer[20] = fitter->GetParError("pLambda",DLM_Fitter1::p_a);

        buffer[21] = fitter->GetParameter("pLambda",DLM_Fitter1::p_b);
        buffer[22] = fitter->GetParError("pLambda",DLM_Fitter1::p_b);

        buffer[23] = fitter->GetParameter("pLambda",DLM_Fitter1::p_c);
        buffer[24] = fitter->GetParError("pLambda",DLM_Fitter1::p_c);

        buffer[25] = fitter->GetParameter("pLambda",DLM_Fitter1::p_Cl);
        buffer[26] = fitter->GetParError("pLambda",DLM_Fitter1::p_Cl);

        buffer[27] = fitter->GetParameter("pLambda",DLM_Fitter1::p_kc);
        buffer[28] = fitter->GetParError("pLambda",DLM_Fitter1::p_kc);

        buffer[29] = fitter->GetParameter("pLambda",DLM_Fitter1::p_sor0);
        buffer[30] = fitter->GetParError("pLambda",DLM_Fitter1::p_sor0);

        buffer[31] = fitter->GetParameter("pLambda",DLM_Fitter1::p_sor1);
        buffer[32] = fitter->GetParError("pLambda",DLM_Fitter1::p_sor1);

        buffer[33] = fitter->GetChi2();
        buffer[34] = fitter->GetNdf();

        buffer[35] = CHI2_312;
        buffer[36] = NDF_312;

        NtFile->cd();
        ntResult->Fill(buffer);
        ntResult->Write("",TObject::kOverwrite);
    }


    delete Ck_pXim;
    delete fitter;

    delete Ck_pL;

    delete Ck_pSigma0;
    delete [] MomBins_pL;
    delete [] FitRegion_pL;
    ////delete hDataClever_pL;
    delete hResolution_pL;
    delete hResidual_pL_pSigma0;
    delete hResidual_pL_pXim;
    delete hData_pL;
    delete OutputFile;
    if(ntResult){delete ntResult; ntResult=NULL;}
    if(NtFile){delete NtFile; NtFile=NULL;}
}



//For feed-down we will use only pSigma0 and pXi, where the former we model as Oli,
//the latter we model with HALQCD, by fixing the source and radius to a Gauss of size 1.18 (-20% for the systematic variation)
//Fitting mode:
    //"Norm": only norm
    //"Baseline": baseline but only in the femto region
    //"Longbaseline": baseline extending to higher k
//if SourceScale||SourceStability==0 => we take the values we have obtained from the pp fit (hard-coded at the moment)
//SourceScale == -1, fit within 10% of the expected value (based on mT or whatever)
//SourceScale == -11, completely free fit
//pL_Pot: "LO",  "LO_Coupled_S", "NLO", "NLO_sp",  "NLO_Coupled_S",    "Usmani"
//in VARIATIONS:
//[0] = femto range
//[1] = fit range
//[2] = fractions for p
//[3] = fractions for L
//[4] = kc parameter
//[5] = which mT bin
//[6] = in case of Reso source: smearing of the mass of the resonance (on or off)//NOT ACTIVE YET
//[7] = in case of Reso source: smearing of the momentum of the resonance (fraction, from 0 to 1)//NOT ACTIVE YET
void Fit_pp(DLM_CommonAnaFunctions& AnalysisObject, const TString& OutputFolder, const TString& OutputFileName,
            const TString& DataSample, const TString& DataVar, const TString& SourceType, const double& SourceScale, const double& SourceStability,
            const TString& pp_Pot,const TString& FittingMode_pp, const int* VARIATIONS, const int& UniqueID=-1, const TString ntFileName="", const TString ntName=""){

    //printf("VARIATIONS[2]=%i\n",VARIATIONS[2]);


    //the residuals that are not fitted, are assumed to all have a Gaussian source of some size;
    double ResidualSourceSize=0;
    if(DataSample=="pp13TeV_MB_Run2paper") ResidualSourceSize=1.1;
    else if(DataSample=="pPb5TeV_Run2paper"||DataSample=="pPb5TeV_CPR_Mar19") ResidualSourceSize=1.4;
    else if(DataSample=="pp13TeV_HM_March19") ResidualSourceSize=1.4;
    else printf("\033[1;31mERROR:\033[0m The data sample '%s' does not exist\n",DataSample.Data());

    double* MomBins_pp = NULL;
    double* FitRegion_pp = NULL;
    unsigned NumMomBins_pp;
    AnalysisObject.SetUpBinning_pp(DataSample,NumMomBins_pp,MomBins_pp,FitRegion_pp,VARIATIONS[0],VARIATIONS[1]);

    TH2F* hResolution_pp = AnalysisObject.GetResolutionMatrix(DataSample,"pp");
    TH2F* hResidual_pp_pL = AnalysisObject.GetResidualMatrix("pp","pLambda");
    TH1F* hData_pp = AnalysisObject.GetAliceExpCorrFun(DataSample,"pp",DataVar,0,false,VARIATIONS[5]);

    double lam_pp[5];
    double lam_pL[5];

    AnalysisObject.SetUpLambdaPars_pp(DataSample,VARIATIONS[2],lam_pp);
    AnalysisObject.SetUpLambdaPars_pL(DataSample,VARIATIONS[2],VARIATIONS[3],lam_pL);

    TString SourceDescription = SourceType;
    if(SourceScale<0) {SourceDescription+="ScaFree";}
    else if(SourceScale==0) {SourceDescription+="ScaDefault";}
    else {SourceDescription+=TString::Format("Sca%.2f",SourceScale);}
    if(SourceStability<0||SourceStability>2) {SourceDescription+="StaFree";}
    else if(SourceStability==0) {SourceDescription+="StaDefault";}
    else {SourceDescription+=TString::Format("Sta%.2f",SourceStability);}
    TString OutFileName = OutputFileName;
    if(OutFileName==""){
        OutFileName=TString::Format("AUTO%i_%s_%s_%s_%s_%i_%i_%i_%i_%i_mT%i.root",UniqueID,
                                                               DataSample.Data(),SourceDescription.Data(),
                                                               pp_Pot.Data(),FittingMode_pp.Data(),
                                                               VARIATIONS[0],VARIATIONS[1],VARIATIONS[2],VARIATIONS[3],VARIATIONS[4],VARIATIONS[5]);
    }
    TFile* OutputFile = new TFile(OutputFolder+OutFileName,"recreate");
    hData_pp->Write();

    CATS AB_pp;
    DLM_Ck* Ck_pp;
    AB_pp.SetMomBins(NumMomBins_pp,MomBins_pp);
AB_pp.SetNotifications(CATS::nWarning);
    AnalysisObject.SetUpCats_pp(AB_pp,pp_Pot,SourceType,0,2);
    if(SourceScale<0) AB_pp.SetAnaSource(0,1.0);
    else AB_pp.SetAnaSource(0,SourceScale);
    if(AB_pp.GetNumSourcePars()>1){
        AB_pp.SetAnaSource(1,SourceStability);
    }
    //AB_pp.SetNotifications(CATS::nWarning);
    Ck_pp = new DLM_Ck(AB_pp.GetNumSourcePars(),0,AB_pp);
    Ck_pp->SetSourcePar(0,SourceStability);
    if(AB_pp.GetNumSourcePars()>1){
//printf("SourceStability=%f\n",SourceStability);
        Ck_pp->SetSourcePar(1,SourceStability);
    }
    Ck_pp->Update();

    CATS AB_pL;
    DLM_Ck* Ck_pL;
    AB_pL.SetMomBins(NumMomBins_pp,MomBins_pp);
    AnalysisObject.SetUpCats_pL(AB_pL,"NLO_Coupled_S","Gauss");
    AB_pL.SetAnaSource(0,ResidualSourceSize);
    AB_pL.SetNotifications(CATS::nWarning);
    //AB_pL.KillTheCat();
    Ck_pL = new DLM_Ck(AB_pL.GetNumSourcePars(),0,AB_pL);
    Ck_pL->SetSourcePar(0,ResidualSourceSize);
    if(AB_pL.GetNumSourcePars()>1){
        Ck_pL->SetSourcePar(1,SourceStability);
    }
    Ck_pL->Update();

//!
    DLM_CkDecomposition CkDec_pp("pp",3,*Ck_pp,hResolution_pp);
    DLM_CkDecomposition CkDec_pL("pLambda",2,*Ck_pL,NULL);

    CkDec_pp.AddContribution(0,lam_pp[1],DLM_CkDecomposition::cFeedDown,&CkDec_pL,hResidual_pp_pL);
    CkDec_pp.AddContribution(1,lam_pp[2],DLM_CkDecomposition::cFeedDown);
    CkDec_pp.AddContribution(2,lam_pp[3],DLM_CkDecomposition::cFake);

    CkDec_pL.AddContribution(2,lam_pL[1]+lam_pL[2]+lam_pL[3],DLM_CkDecomposition::cFeedDown);
    CkDec_pL.AddContribution(3,lam_pL[4],DLM_CkDecomposition::cFake);//0.03

    DLM_Fitter1* fitter = new DLM_Fitter1(1);
    if(FittingMode_pp.Contains("Norm")||FittingMode_pp.Contains("Baseline")){
        fitter->SetSystem(0,*hData_pp,1,CkDec_pp,
                    FitRegion_pp[0],FitRegion_pp[1],FitRegion_pp[1],FitRegion_pp[1]);
    }
    else if(FittingMode_pp.Contains("Longbaseline")||FittingMode_pp.Contains("Spline")){
        fitter->SetSystem(0,*hData_pp,1,CkDec_pp,
                        FitRegion_pp[0],FitRegion_pp[1],FitRegion_pp[2],FitRegion_pp[3]);
    }
    else{
        printf("\033[1;31mERROR:\033[0m Unknown fitting mode '%s'",FittingMode_pp.Data());
    }

    if(FittingMode_pp.Contains("_prefit")||VARIATIONS[4]<0){
        fitter->SetSeparateBL(0,true);
        fitter->FixParameter("pp",DLM_Fitter1::p_Cl,-1);
    }
    else{
        fitter->SetSeparateBL(0,false);

        if(VARIATIONS[4]>0&&FittingMode_pp!="Norm") fitter->FixParameter("pp",DLM_Fitter1::p_kc,VARIATIONS[4]);
        else fitter->SetParameter("pp",DLM_Fitter1::p_kc,2000,400,10000);

    }

    printf("FittingMode_pp=%s\n",FittingMode_pp.Data());

    //!p_a
    fitter->SetParameter("pp",DLM_Fitter1::p_a,1.0,0.7,1.3);
    if(FittingMode_pp.Contains("Spline"))
        fitter->FixParameter("pp",DLM_Fitter1::p_a,1.0);

    //!p_b
    if(FittingMode_pp.Contains("Baseline")||FittingMode_pp.Contains("Longbaseline"))
        fitter->SetParameter("pp",DLM_Fitter1::p_b,0,-2e-3,2e-3);
    else
        fitter->FixParameter("pp",DLM_Fitter1::p_b,0);

    //!p_c
    if(FittingMode_pp.Contains("Baseline2")||FittingMode_pp.Contains("Longbaseline2"))
        fitter->SetParameter("pp",DLM_Fitter1::p_c,0,-2e-4,2e-4);
    else
        fitter->FixParameter("pp",DLM_Fitter1::p_c,0);

    if(FittingMode_pp.Contains("Spline")){
        int NumKnots=0;
        if(FittingMode_pp.Contains("_2")) NumKnots=2;
        if(FittingMode_pp.Contains("_3")) NumKnots=3;
        if(FittingMode_pp.Contains("_4")) NumKnots=4;
        if(FittingMode_pp.Contains("_5")) NumKnots=5;
        if(FittingMode_pp.Contains("_6")) NumKnots=6;
        if(FittingMode_pp.Contains("_7")) NumKnots=7;
        if(FittingMode_pp.Contains("_8")) NumKnots=8;
        if(FittingMode_pp.Contains("_9")) NumKnots=9;
        if(FittingMode_pp.Contains("_10")) NumKnots=10;
        if(FittingMode_pp.Contains("_11")) NumKnots=11;
        if(FittingMode_pp.Contains("_Tune1")) NumKnots=5;
        fitter->FixParameter("pp",DLM_Fitter1::p_spline,NumKnots);
        fitter->FixParameter("pp",DLM_Fitter1::p_spline+1,0);
        fitter->SetParameter("pp",DLM_Fitter1::p_spline+2,0);

        double* Nodes_x = new double [NumKnots];
        Nodes_x[0] = FitRegion_pp[0];
        Nodes_x[NumKnots-1] = FitRegion_pp[3];
        if(FittingMode_pp.Contains("_Tune1")){
            Nodes_x[0] = 0;
            Nodes_x[1] = 90;
            Nodes_x[2] = 180;
            Nodes_x[3] = 360;
            Nodes_x[4] = 540;
            fitter->FixParameter("pp",DLM_Fitter1::p_spline+2,0.00031);
        }
        else{
            double NodeLength = (Nodes_x[NumKnots-1]-Nodes_x[0])/double(NumKnots-1);
            //printf("NodeLength=%f\n",NodeLength);
            for(int iKnot=1; iKnot<NumKnots-1; iKnot++){
                Nodes_x[iKnot] = Nodes_x[iKnot-1]+NodeLength;
                //printf("Nodes_x[%i]=%f\n",iKnot,Nodes_x[iKnot]);
            }
        }

        for(int iKnot=0; iKnot<NumKnots; iKnot++){
            fitter->FixParameter("pp",DLM_Fitter1::p_spline+3+iKnot,Nodes_x[iKnot]);//x
            fitter->SetParameter("pp",DLM_Fitter1::p_spline+3+NumKnots+iKnot,1,0,5);//y
        }
        delete [] Nodes_x;
    }

    //!p_sor0
    if(SourceType=="Gauss"){
        if(SourceScale==0){fitter->FixParameter("pp",DLM_Fitter1::p_sor0,1.2);}
        else if(SourceScale<0) {fitter->SetParameter("pp",DLM_Fitter1::p_sor0,1.2,0.5,2.0);}
        else {fitter->FixParameter("pp",DLM_Fitter1::p_sor0,SourceScale);}
    }
    else if(SourceType=="CleverLevy_Nolan"||SourceType=="CleverLevy_Single"||SourceType=="CleverLevy_Diff"){
        fitter->SetParameter("pp",DLM_Fitter1::p_sor0,1.0,0.3,1.7);
        fitter->SetParameter("pp",DLM_Fitter1::p_sor1,1.4,1.0,2.0);
    }
    else if(SourceType=="McLevyNolan_Reso"){
        if(DataSample=="pp13TeV_MB_Run2paper"){
            if(SourceScale==0) {fitter->FixParameter("pp",DLM_Fitter1::p_sor0,0.71);}//at mT 1.44
            else if(SourceScale==-1) {fitter->SetParameter("pp",DLM_Fitter1::p_sor0,0.71,0.71/2.,0.71*2.);}
            else if(SourceScale==-11) {fitter->SetParameter("pp",DLM_Fitter1::p_sor0,0.71,0.71/2.,0.71*2.);}
            else {fitter->FixParameter("pp",DLM_Fitter1::p_sor0,SourceScale);}

            if(SourceStability==0) {fitter->FixParameter("pp",DLM_Fitter1::p_sor1,1.5);}
            else if(SourceStability==-1) {fitter->SetParameter("pp",DLM_Fitter1::p_sor1,1.5,1.5*0.90,1.5*1.10);}
            else if(SourceStability==-11) {fitter->SetParameter("pp",DLM_Fitter1::p_sor1,1.5,1.0,2.0);}
            else {fitter->FixParameter("pp",DLM_Fitter1::p_sor1,SourceStability);}
        }
        else if(DataSample=="pPb5TeV_Run2paper"||DataSample=="pPb5TeV_CPR_Mar19"){
            if(SourceScale==0) {fitter->FixParameter("pp",DLM_Fitter1::p_sor0,1.02);}//at mT 1.44
            else if(SourceScale==-1) {fitter->SetParameter("pp",DLM_Fitter1::p_sor0,1.02,1.02*0.90,1.02*1.10);}
            else if(SourceScale==-11) {fitter->SetParameter("pp",DLM_Fitter1::p_sor0,1.02,1.02/2.,1.02*2.);}
            else {fitter->FixParameter("pp",DLM_Fitter1::p_sor0,SourceScale);}

            if(SourceStability==0) {fitter->FixParameter("pp",DLM_Fitter1::p_sor1,1.45);}
            else if(SourceStability==-1) {fitter->SetParameter("pp",DLM_Fitter1::p_sor1,1.45,1.45*0.90,1.45*1.10);}
            else if(SourceStability==-11) {fitter->SetParameter("pp",DLM_Fitter1::p_sor1,1.5,1.0,2.0);}
            else {fitter->FixParameter("pp",DLM_Fitter1::p_sor1,SourceStability);}
        }
        else if(DataSample=="pp13TeV_HM_March19"){
            if(SourceScale==0) {fitter->FixParameter("pp",DLM_Fitter1::p_sor0,0.92);}//at mT 1.55
            else if(SourceScale==-1) {fitter->SetParameter("pp",DLM_Fitter1::p_sor0,0.92,0.92*0.90,0.92*1.10);}
            else if(SourceScale==-11) {fitter->SetParameter("pp",DLM_Fitter1::p_sor0,0.92,0.92/2.,0.92*2.);}
            else {fitter->FixParameter("pp",DLM_Fitter1::p_sor0,SourceScale);}

            if(SourceStability==0) {fitter->FixParameter("pp",DLM_Fitter1::p_sor1,1.5);}
            else if(SourceStability==-1) {fitter->SetParameter("pp",DLM_Fitter1::p_sor1,1.5,1.5*0.90,1.5*1.10);}
            else if(SourceStability==-11) {fitter->SetParameter("pp",DLM_Fitter1::p_sor1,1.5,1.0,2.0);}
            else {fitter->FixParameter("pp",DLM_Fitter1::p_sor1,SourceStability);}
        }
        else printf("WHAT HAPPEND?\n");
    }
    else if(SourceType=="McGauss_Reso"){
        if(DataSample=="pp13TeV_MB_Run2paper"){
            if(SourceScale==0) {fitter->FixParameter("pp",DLM_Fitter1::p_sor0,0.71);}//at mT 1.44
            else if(SourceScale==-1) {fitter->SetParameter("pp",DLM_Fitter1::p_sor0,0.71,0.71*0.90,0.71*1.10);}
            else if(SourceScale==-11) {fitter->SetParameter("pp",DLM_Fitter1::p_sor0,0.71,0.71/2.,0.71*2.);}
            else {fitter->FixParameter("pp",DLM_Fitter1::p_sor0,SourceScale);}
            fitter->FixParameter("pp",DLM_Fitter1::p_sor1,2.0);
        }
        else if(DataSample=="pPb5TeV_Run2paper"||DataSample=="pPb5TeV_CPR_Mar19"){
            if(SourceScale==0) {fitter->FixParameter("pp",DLM_Fitter1::p_sor0,0.95);}//at mT 1.52
            else if(SourceScale==-1) {fitter->SetParameter("pp",DLM_Fitter1::p_sor0,0.95,0.95*0.90,0.95*1.10);}
            else if(SourceScale==-11) {fitter->SetParameter("pp",DLM_Fitter1::p_sor0,0.95,0.95/2.,0.95*2.);}
            else {fitter->FixParameter("pp",DLM_Fitter1::p_sor0,SourceScale);}
            fitter->FixParameter("pp",DLM_Fitter1::p_sor1,2.0);
        }
        else if(DataSample=="pp13TeV_HM_March19"){
            //if(SourceScale==0) {fitter->FixParameter("pp",DLM_Fitter1::p_sor0,0.87);}//at mT 1.55
            //else if(SourceScale==-1) {fitter->SetParameter("pp",DLM_Fitter1::p_sor0,0.87,0.87*0.90,0.87*1.10);}
            //else if(SourceScale==-11) {fitter->SetParameter("pp",DLM_Fitter1::p_sor0,0.87,0.87/2.,0.87*2.);}
            //else {fitter->FixParameter("pp",DLM_Fitter1::p_sor0,SourceScale);}

            //these are the value shown at the PF previews
            //if(SourceScale==0) {fitter->FixParameter("pp",DLM_Fitter1::p_sor0,0.918078);}//at mT 1.55
            //else if(SourceScale==-1) {fitter->SetParameter("pp",DLM_Fitter1::p_sor0,0.918078,0.90218,0.933344);}
            if(SourceScale==0) {fitter->FixParameter("pp",DLM_Fitter1::p_sor0,0.918);}//at mT 1.55
            else if(SourceScale==-1) {fitter->SetParameter("pp",DLM_Fitter1::p_sor0,0.918,0.871,0.965);}//3 sigma
            else if(SourceScale==-11) {fitter->SetParameter("pp",DLM_Fitter1::p_sor0,0.87,0.87/2.,0.87*2.);}
            else {fitter->FixParameter("pp",DLM_Fitter1::p_sor0,SourceScale);}

            fitter->FixParameter("pp",DLM_Fitter1::p_sor1,2.0);
        }
        else printf("WHAT HAPPEND?\n");
    }
    else{
        printf("\033[1;31mERROR:\033[0m '%s' does not exist",SourceType.Data());
    }

    CkDec_pp.Update();
    CkDec_pL.Update();
    fitter->GoBabyGo(false);
    TGraph FitResult_pp;
    FitResult_pp.SetName(TString::Format("FitResult_pp"));
    fitter->GetFitGraph(0, FitResult_pp);
    //TString Description = TString::Format("%s_%s_%s",pp_Pot[uPot_pp].Data(),pL_Pot[uPot_pL].Data(),pXim_Pot[uPot_pXim].Data());
    OutputFile->cd();
    FitResult_pp.Write();

    TGraph FitBaseline_pp;
    FitBaseline_pp.SetName(TString::Format("FitBaseline_pp"));
    fitter->GetMultBaselineGraph(0, FitBaseline_pp);
    //TString Description = TString::Format("%s_%s_%s",pp_Pot[uPot_pp].Data(),pL_Pot[uPot_pL].Data(),pXim_Pot[uPot_pXim].Data());
    OutputFile->cd();
    FitBaseline_pp.Write();

    double CHI2_352=0;
    unsigned NDF_352=0;
    double xVal;
    double yVal;
    unsigned hBin;
    for(unsigned uPoint=0; uPoint<FitResult_pp.GetN(); uPoint++){
        FitResult_pp.GetPoint(uPoint,xVal,yVal);
        if(xVal>352) break;
        hBin = hData_pp->FindBin(xVal);
        CHI2_352 += pow( (hData_pp->GetBinContent(hBin)-yVal)/(hData_pp->GetBinError(hBin)) ,2.);
        NDF_352++;
    }

    printf("χ2/ndf = %.2f/%i = %.2f\n",fitter->GetChi2(),fitter->GetNdf(),fitter->GetChi2Ndf());
    printf("χ2/ndf (352) = %.2f/%i = %.2f\n",CHI2_352,NDF_352,CHI2_352/double(NDF_352));
    printf("R(pΛ) = %.3f +/- %.3f\n",fitter->GetParameter("pp",DLM_Fitter1::p_sor0),fitter->GetParError("pp",DLM_Fitter1::p_sor0));
    if(AB_pp.GetNumSourcePars()>1){
        printf("α(pΛ) = %.3f +/- %.3f\n",fitter->GetParameter("pp",DLM_Fitter1::p_sor1),fitter->GetParError("pp",DLM_Fitter1::p_sor1));
    }
    printf("a = %.2e +/- %.2e\n",fitter->GetParameter("pp",DLM_Fitter1::p_a),fitter->GetParError("pp",DLM_Fitter1::p_a));
    printf("b = %.2e +/- %.2e\n",fitter->GetParameter("pp",DLM_Fitter1::p_b),fitter->GetParError("pp",DLM_Fitter1::p_b));
    printf("c = %.2e +/- %.2e\n",fitter->GetParameter("pp",DLM_Fitter1::p_c),fitter->GetParError("pp",DLM_Fitter1::p_c));
    printf("Cl = %.3f +/- %.3f\n",fitter->GetParameter("pp",DLM_Fitter1::p_Cl),fitter->GetParError("pp",DLM_Fitter1::p_Cl));
    printf("kc = %.3f +/- %.3f\n",fitter->GetParameter("pp",DLM_Fitter1::p_kc),fitter->GetParError("pp",DLM_Fitter1::p_kc));
    printf("mT = %i\n",VARIATIONS[5]);

    TGraph gChi2Ndf;
    gChi2Ndf.SetName("gChi2Ndf");
    gChi2Ndf.Set(1);
    gChi2Ndf.SetPoint(0,fitter->GetChi2(),fitter->GetNdf());
    gChi2Ndf.SetMarkerStyle(2);
    gChi2Ndf.SetMarkerSize(4);
    gChi2Ndf.Write();

    TF1* fBaseline = new TF1("fBaseline","[0]+[1]*x+[2]*x*x",0,600);
    fBaseline->FixParameter(0,fitter->GetParameter("pp",DLM_Fitter1::p_a));
    fBaseline->FixParameter(1,fitter->GetParameter("pp",DLM_Fitter1::p_b));
    fBaseline->FixParameter(2,fitter->GetParameter("pp",DLM_Fitter1::p_c));
    fBaseline->SetLineColor(kCyan+4);
    fBaseline->SetLineWidth(3);
    fBaseline->SetLineStyle(7);
    fBaseline->Write();

    TGraphErrors gScaleStability;
    gScaleStability.SetName("gScaleStability");
    gScaleStability.Set(1);
    double SCALE = fitter->GetParameter("pp",DLM_Fitter1::p_sor0);
    double dSCALE = fitter->GetParError("pp",DLM_Fitter1::p_sor0);
    double STABILITY = AB_pp.GetNumSourcePars()>1?fitter->GetParameter("pp",DLM_Fitter1::p_sor1):2;
    double dSTABILITY = AB_pp.GetNumSourcePars()>1?fitter->GetParError("pp",DLM_Fitter1::p_sor1):0;
    gScaleStability.SetPoint(0,SCALE,STABILITY);
    gScaleStability.SetPointError(0,dSCALE,dSTABILITY);
    gScaleStability.SetMarkerStyle(20);
    gScaleStability.SetMarkerSize(1);
    gScaleStability.SetLineWidth(2);
    gScaleStability.Write();

    //TFile* OutputFileNT = new TFile(OutputFolder+TString::Format("NT_%s_%s_%s_%s_%i_%i_%i_%i_%i.root",
    //                                                           DataSample.Data(),SourceDescription.Data(),
    //                                                           pL_Pot.Data(),FittingMode_pL.Data(),
    //                                                           VARIATIONS[0],VARIATIONS[1],VARIATIONS[2],VARIATIONS[3],VARIATIONS[4]),"recreate");

    TFile* NtFile = NULL;
    TNtuple* ntResult = NULL;
    if(ntFileName!=""&&ntName!=""){
        bool FileIsOpen=false;
        bool FileDoesNotExist=false;
        double WaitedInSec=0;
        do{
            if(FileIsOpen) {usleep(1000); WaitedInSec+=1e-3;}//sleep for a millisecond
            if(WaitedInSec>60){
                printf("\033[1;31mERROR: \033[0mCannot access the file %s\n", ntFileName.Data());
                printf("\033[1;31m ABORTING THE EXECUTABLE!\033[0m\n");
                abort();
            }
            int fd = open(ntFileName.Data(), O_RDONLY);
            if (fd < 0) {
                //perror("open");
                //file does not exist
                FileIsOpen = false;
                FileDoesNotExist = true;
            }
            if (fcntl(fd, F_SETLEASE, F_WRLCK) && EAGAIN == errno) {
                FileIsOpen = true;
            }
            else {
                fcntl(fd, F_SETLEASE, F_UNLCK);
                FileIsOpen = false;
            }
            close(fd);
        }
        while(FileIsOpen);
        if(FileDoesNotExist) NtFile = new TFile(ntFileName,"RECREATE");
        else NtFile = new TFile(ntFileName,"UPDATE");
        ntResult = (TNtuple*)NtFile->Get(ntName);
        if(!ntResult){
            //printf("\033[1;31mERROR:\033[0m For whatever reason the ");
            ntResult = new TNtuple("ntResult", "ntResult",
            "IterM:Iter:Config:Data:SourceType:SourceScale:SourceStability:Potential:Baseline:FemRan:FitRan:pFrac:LamFrac:kcVar:mTbin:FemtoMin:FemtoMax:BlMin:BlMax:"
            "p_a:e_a:p_b:e_b:p_c:e_c:p_Cl:e_Cl:p_kc:e_kc:p_sor0:e_sor0:p_sor1:e_sor1:chi2:ndf:chi2_352:ndf_352");
        }
    }

    if(NtFile&&ntResult){
        Float_t buffer[37];
        long Iter = UniqueID/1000;
        long Config = UniqueID%1000;
        buffer[0] = Iter/long(1000000);
        buffer[1] = Iter%long(1000000);
        buffer[2] = UniqueID%1000;
        TString Temp = DataVar;
        Temp.ReplaceAll("_","");
        buffer[3] = Temp.Atoi();
        if(SourceType=="Gauss") buffer[4] = 0;
        if(SourceType=="McGauss_Reso") buffer[4] = 1;
        if(SourceType=="McLevyNolan_Reso") buffer[4] = 2;

        buffer[5] = SourceScale;
        buffer[6] = SourceStability;

        if(pp_Pot=="AV18") buffer[7] = 0;

        if(FittingMode_pp=="Norm") buffer[8] = 0;
        if(FittingMode_pp=="Baseline") buffer[8] = 1;
        if(FittingMode_pp=="Baseline2") buffer[8] = 2;
        if(FittingMode_pp=="Longbaseline") buffer[8] = 11;
        if(FittingMode_pp=="Longbaseline2") buffer[8] = 12;
        if(FittingMode_pp=="Norm_prefit") buffer[8] = 100;
        if(FittingMode_pp=="Baseline_prefit") buffer[8] = 101;
        if(FittingMode_pp=="Baseline2_prefit") buffer[8] = 102;
        if(FittingMode_pp=="Longbaseline_prefit") buffer[8] = 111;
        if(FittingMode_pp=="Longbaseline2_prefit") buffer[8] = 112;
        if(FittingMode_pp=="Spline3_3") buffer[8] = 403;
        if(FittingMode_pp=="Spline3_4") buffer[8] = 404;
        if(FittingMode_pp=="Spline3_5") buffer[8] = 405;
        if(FittingMode_pp=="Spline3_6") buffer[8] = 406;
        if(FittingMode_pp=="Spline3_7") buffer[8] = 407;
        if(FittingMode_pp=="Spline3_8") buffer[8] = 408;
        if(FittingMode_pp=="Spline3_9") buffer[8] = 409;
        if(FittingMode_pp=="Spline3_10") buffer[8] = 410;
        if(FittingMode_pp=="Spline3_11") buffer[8] = 411;
        if(FittingMode_pp=="Spline3_Tune1") buffer[8] = 451;

        buffer[9] = VARIATIONS[0];
        buffer[10] = VARIATIONS[1];
        buffer[11] = VARIATIONS[2];
        buffer[12] = VARIATIONS[3];
        buffer[13] = VARIATIONS[4];
        buffer[14] = VARIATIONS[5];

        buffer[15] = FitRegion_pp[0];
        buffer[16] = FitRegion_pp[1];
        buffer[17] = FitRegion_pp[2];
        buffer[18] = FitRegion_pp[3];

        buffer[19] = fitter->GetParameter("pp",DLM_Fitter1::p_a);
        buffer[20] = fitter->GetParError("pp",DLM_Fitter1::p_a);

        buffer[21] = fitter->GetParameter("pp",DLM_Fitter1::p_b);
        buffer[22] = fitter->GetParError("pp",DLM_Fitter1::p_b);

        buffer[23] = fitter->GetParameter("pp",DLM_Fitter1::p_c);
        buffer[24] = fitter->GetParError("pp",DLM_Fitter1::p_c);

        buffer[25] = fitter->GetParameter("pp",DLM_Fitter1::p_Cl);
        buffer[26] = fitter->GetParError("pp",DLM_Fitter1::p_Cl);

        buffer[27] = fitter->GetParameter("pp",DLM_Fitter1::p_kc);
        buffer[28] = fitter->GetParError("pp",DLM_Fitter1::p_kc);

        buffer[29] = fitter->GetParameter("pp",DLM_Fitter1::p_sor0);
        buffer[30] = fitter->GetParError("pp",DLM_Fitter1::p_sor0);

        buffer[31] = fitter->GetParameter("pp",DLM_Fitter1::p_sor1);
        buffer[32] = fitter->GetParError("pp",DLM_Fitter1::p_sor1);

        buffer[33] = fitter->GetChi2();
        buffer[34] = fitter->GetNdf();

        buffer[35] = CHI2_352;
        buffer[36] = NDF_352;

        NtFile->cd();
        ntResult->Fill(buffer);
        ntResult->Write("",TObject::kOverwrite);
    }


    delete Ck_pp;
    delete fitter;

    delete Ck_pL;

    delete [] MomBins_pp;
    delete [] FitRegion_pp;
    ////delete hDataClever_pL;
    delete hResolution_pp;
    delete hResidual_pp_pL;
    delete hData_pp;
    delete OutputFile;
    if(ntResult){delete ntResult; ntResult=NULL;}
    if(NtFile){delete NtFile; NtFile=NULL;}

}





/*
void Plot_pp_pL(const TString& DataSample, const TString& SourceType, const TString& FittingMode_pp, const TString& FittingMode_pL){
    TFile* InputFile = new TFile(OutputFolder+TString::Format("Output_%s_%s_%s.root",DataSample.Data(),SourceType.Data(),FittingMode_pp.Data()),"read");




    delete InputFile;
}
*/

//if OnlyFraction>1 => we have 1/OnlyFraction probability to accept an event
//Mode==-1 => the default cut comb only (ignoring all first, last iter etc.)
//Mode==0 => normal
//Mode==1 => we only print the number of iter to be performed
void pL_SystematicsHM(const TString& OutputFolder, const int& WhichConfiguration, const int& FirstIter, const int& LastIter, const unsigned& OnlyFraction,
                      const int& RANDOMSEED, const int& Mode){

    //const bool Perform_mT_Scan = true;

    printf("Output: %s\n",OutputFolder.Data());
    printf("Iter: %i -> %i\n",FirstIter,LastIter);
    printf("Fraction: %i\n",OnlyFraction);
    printf("Seed: %i\n",RANDOMSEED);
    DLM_CommonAnaFunctions AnalysisObject;
    int NumSourceVars;
    int NumSourceScaleVars;
    int NumSourceStabilityVars;
    int NumPotVars;
    int NumBaselineVars;
    int NumFemtoRangeVars;
    int NumFitRangeVars;
    int NumProtonFracVars;
    int NumLambdaFracVars;
    int NumKcVars;
    int NumMtVars;
    int NumDataVars;//these are the cut variations

    int DefSource;
    int DefSourceScale;
    int DefSourceStability;
    int DefPot;
    int DefBaseline;
    int DefFemtoRange;
    int DefFitRange;
    int DefProtonFrac;
    int DefLambdaFrac;
    int DefKc;
    int DefMt;
    int DefData;

    TString* Source;
    double* SourceScale;
    double* SourceStability;
    TString* Potential;
    TString* Baseline;
    int* FemtoRangeVars;
    int* FitRangeVars;
    int* pFracVars;
    int* LamFracVars;
    double* kcVars;
    int* MtVars;
    TString* DataVars;
    int VARIATIONS[6];
    int TotNumIter=0;

    switch(WhichConfiguration){
        //this is the case where we only check for the Gauss+Reso and Gauss
        case 0 :
            NumSourceVars = 1;
            NumSourceScaleVars = 3;
            NumSourceStabilityVars = 1;
            NumPotVars = 2;
            NumBaselineVars = 3;
            NumFemtoRangeVars = 3;
            NumFitRangeVars = 3;
            NumProtonFracVars = 3;
            NumLambdaFracVars = 9;
            NumKcVars = 3;
            NumMtVars = 1;
            NumDataVars = 45;
            break;
        //similar to 0, only that now we include Levy+Reso as an option
        case 1 :
            NumSourceVars = 1;
            NumSourceScaleVars = 4;
            NumSourceStabilityVars = 3;
            NumPotVars = 2;
            NumBaselineVars = 3;
            NumFemtoRangeVars = 3;
            NumFitRangeVars = 3;
            NumProtonFracVars = 3;
            NumLambdaFracVars = 9;
            NumKcVars = 3;
            NumMtVars = 1;
            NumDataVars = 45;
            break;
        case 2 :
            NumSourceVars = 1;
            NumSourceScaleVars = 3;
            NumSourceStabilityVars = 1;
            NumPotVars = 2;
            NumBaselineVars = 3;
            NumFemtoRangeVars = 3;
            NumFitRangeVars = 3;
            NumProtonFracVars = 3;
            NumLambdaFracVars = 9;
            NumKcVars = 3;
            NumMtVars = 1;
            NumDataVars = 45;
            break;
        //without the peak region
        case 10 :
            NumSourceVars = 1;
            NumSourceScaleVars = 4;
            NumSourceStabilityVars = 1;
            NumPotVars = 2;
            NumBaselineVars = 3;
            NumFemtoRangeVars = 3;
            NumFitRangeVars = 3;
            NumProtonFracVars = 3;
            NumLambdaFracVars = 9;
            NumKcVars = 3;
            NumMtVars = 1;
            NumDataVars = 1;
            break;
        //similar to 0, only that now we include Levy+Reso as an option
        case 11 :
            NumSourceVars = 1;
            NumSourceScaleVars = 4;
            NumSourceStabilityVars = 3;
            NumPotVars = 2;
            NumBaselineVars = 3;
            NumFemtoRangeVars = 3;
            NumFitRangeVars = 3;
            NumProtonFracVars = 3;
            NumLambdaFracVars = 9;
            NumKcVars = 3;
            NumMtVars = 1;
            NumDataVars = 1;
            break;
        case 12 :
            NumSourceVars = 1;
            NumSourceScaleVars = 3;
            NumSourceStabilityVars = 1;
            NumPotVars = 2;
            NumBaselineVars = 3;
            NumFemtoRangeVars = 3;
            NumFitRangeVars = 3;
            NumProtonFracVars = 3;
            NumLambdaFracVars = 9;
            NumKcVars = 3;
            NumMtVars = 1;
            NumDataVars = 1;
            break;
        //this is the case where we only check for the Gauss+Reso and Gauss
        //we fit with mT bins
        case 20 :
            NumSourceVars = 1;
            NumSourceScaleVars = 1;
            NumSourceStabilityVars = 1;
            NumPotVars = 2;
            NumBaselineVars = 3;
            NumFemtoRangeVars = 3;
            NumFitRangeVars = 3;
            NumProtonFracVars = 3;
            NumLambdaFracVars = 9;
            NumKcVars = 3;
            NumMtVars = 7;
            NumDataVars = 45;
            break;
        case 21 :
            NumSourceVars = 1;
            NumSourceScaleVars = 1;
            NumSourceStabilityVars = 1;
            NumPotVars = 2;
            NumBaselineVars = 3;
            NumFemtoRangeVars = 3;
            NumFitRangeVars = 3;
            NumProtonFracVars = 3;
            NumLambdaFracVars = 9;
            NumKcVars = 3;
            NumMtVars = 7;
            NumDataVars = 45;
            break;
        case 22 :
            NumSourceVars = 1;
            NumSourceScaleVars = 1;
            NumSourceStabilityVars = 1;
            NumPotVars = 2;
            NumBaselineVars = 3;
            NumFemtoRangeVars = 3;
            NumFitRangeVars = 3;
            NumProtonFracVars = 3;
            NumLambdaFracVars = 9;
            NumKcVars = 3;
            NumMtVars = 7;
            NumDataVars = 45;
            break;
        //this is the case where we only check for the Gauss+Reso and Gauss
        //we fit with mT bins, but fit without the cusp and check LO,NLO with and without coupled channel, as well as Usmani
        case 30 :
            NumSourceVars = 1;
            NumSourceScaleVars = 1;
            NumSourceStabilityVars = 1;
            NumPotVars = 5;
            NumBaselineVars = 3;
            NumFemtoRangeVars = 3;
            NumFitRangeVars = 3;
            NumProtonFracVars = 3;
            NumLambdaFracVars = 9;
            NumKcVars = 1;
            NumMtVars = 7;
            NumDataVars = 1;//45
            break;
        case 31 :
            NumSourceVars = 1;
            NumSourceScaleVars = 1;
            NumSourceStabilityVars = 1;
            NumPotVars = 5;
            NumBaselineVars = 3;
            NumFemtoRangeVars = 3;
            NumFitRangeVars = 3;
            NumProtonFracVars = 3;
            NumLambdaFracVars = 9;
            NumKcVars = 1;
            NumMtVars = 7;
            NumDataVars = 1;//45
            break;
        case 32 :
            NumSourceVars = 1;
            NumSourceScaleVars = 1;
            NumSourceStabilityVars = 1;
            NumPotVars = 5;
            NumBaselineVars = 3;
            NumFemtoRangeVars = 3;
            NumFitRangeVars = 3;
            NumProtonFracVars = 3;
            NumLambdaFracVars = 9;
            NumKcVars = 1;
            NumMtVars = 7;
            NumDataVars = 1;//45
            break;
        //like 30 but for experimenting
        //e.g. now it is 1 mT bin, default data, all variations as bernie
        case 35 :
            NumSourceVars = 1;
            NumSourceScaleVars = 1;
            NumSourceStabilityVars = 1;
            NumPotVars = 3;
            NumBaselineVars = 2;
            NumFemtoRangeVars = 3;
            NumFitRangeVars = 3;
            NumProtonFracVars = 3;
            NumLambdaFracVars = 9;
            NumKcVars = 1;
            NumMtVars = 1;
            NumDataVars = 1;//45
            break;
        //splines
        case 40 :
            NumSourceVars = 1;
            NumSourceScaleVars = 1;
            NumSourceStabilityVars = 1;
            NumPotVars = 2;
            NumBaselineVars = 1;
            NumFemtoRangeVars = 3;
            NumFitRangeVars = 3;
            NumProtonFracVars = 3;
            NumLambdaFracVars = 9;
            NumKcVars = 3;
            NumMtVars = 1;
            NumDataVars = 1;
            break;
        case 41 :
            NumSourceVars = 1;
            NumSourceScaleVars = 1;
            NumSourceStabilityVars = 1;
            NumPotVars = 1;
            NumBaselineVars = 1;
            NumFemtoRangeVars = 3;
            NumFitRangeVars = 3;
            NumProtonFracVars = 3;
            NumLambdaFracVars = 9;
            NumKcVars = 3;
            NumMtVars = 1;
            NumDataVars = 1;
            break;
        case 42 :
            NumSourceVars = 1;
            NumSourceScaleVars = 1;
            NumSourceStabilityVars = 1;
            NumPotVars = 1;
            NumBaselineVars = 1;
            NumFemtoRangeVars = 3;
            NumFitRangeVars = 3;
            NumProtonFracVars = 3;
            NumLambdaFracVars = 9;
            NumKcVars = 3;
            NumMtVars = 1;
            NumDataVars = 1;
            break;
        //the case of meson exchange models
        //without the peak region (Lednicky)
        case 112 :
            NumSourceVars = 1;
            NumSourceScaleVars = 3;
            NumSourceStabilityVars = 1;
            NumPotVars = 15;
            NumBaselineVars = 3;
            NumFemtoRangeVars = 3;
            NumFitRangeVars = 3;
            NumProtonFracVars = 3;
            NumLambdaFracVars = 9;

            NumKcVars = 1;
            NumMtVars = 1;
            NumDataVars = 1;
            break;
        //same as 0,1,2 but we FIT the radius (not fix it!)
        case 200 :
            NumSourceVars = 1;
            NumSourceScaleVars = 1;
            NumSourceStabilityVars = 1;
            NumPotVars = 2;
            NumBaselineVars = 3;
            NumFemtoRangeVars = 3;
            NumFitRangeVars = 3;
            NumProtonFracVars = 3;
            NumLambdaFracVars = 9;
            NumKcVars = 3;
            NumMtVars = 1;
            NumDataVars = 45;
            break;
        case 201 :
            NumSourceVars = 1;
            NumSourceScaleVars = 1;
            NumSourceStabilityVars = 3;
            NumPotVars = 2;
            NumBaselineVars = 3;
            NumFemtoRangeVars = 3;
            NumFitRangeVars = 3;
            NumProtonFracVars = 3;
            NumLambdaFracVars = 9;
            NumKcVars = 3;
            NumMtVars = 1;
            NumDataVars = 45;
            break;
        case 202 :
            NumSourceVars = 1;
            NumSourceScaleVars = 1;
            NumSourceStabilityVars = 1;
            NumPotVars = 2;
            NumBaselineVars = 3;
            NumFemtoRangeVars = 3;
            NumFitRangeVars = 3;
            NumProtonFracVars = 3;
            NumLambdaFracVars = 9;
            NumKcVars = 3;
            NumMtVars = 1;
            NumDataVars = 45;
            break;
        default : return;
    }

    Source = new TString [NumSourceVars];
    SourceScale = new double [NumSourceScaleVars];
    SourceStability = new double [NumSourceStabilityVars];
    Potential = new TString [NumPotVars];
    Baseline = new TString [NumBaselineVars];
    FemtoRangeVars = new int [NumFemtoRangeVars];
    FitRangeVars = new int [NumFitRangeVars];
    pFracVars = new int [NumProtonFracVars];
    LamFracVars = new int [NumLambdaFracVars];
    kcVars = new double [NumLambdaFracVars];
    MtVars = new int [NumMtVars];
    DataVars = new TString [NumDataVars];

    int** DEFAULT_VAR = new int* [NumSourceVars];
    for(int iSource=0; iSource<NumSourceVars; iSource++){
        DEFAULT_VAR[iSource]=new int [NumPotVars];
        for(int iPot=0; iPot<NumPotVars; iPot++)DEFAULT_VAR[iSource][iPot]=-1;
    }

    //TFile* NTfile = NULL;
    //TNtuple* ntResult = NULL;

    if(WhichConfiguration==0){
        Source[0] = "McGauss_Reso";
        DefSource = 0;
        //SourceScale[0] = -1;//free fit, but within the uncertainties
        //initially it was 0.87 +/- 10 %
        //with wrong systematics: 0.871,0.918,0.965
        SourceScale[0] = 0.853;
        SourceScale[1] = 0.878;
        SourceScale[2] = 0.902;
        SourceScale[2] = -11;
        // pLambda:   low 0.90218  mean 0.918078 up 0.933344 //pp + pL
        // pLambda:   low 0.872613  mean 0.882538 up 0.892688 //pp
        // new syst: pLambda(1.55):   low 0.853396 mean 0.877717 up 0.90239
        DefSourceScale = 0;
        SourceStability[0] = 0;//dummy
        DefSourceStability = 0;
        Potential[0] = "LO_Coupled_S";
        Potential[1] = "NLO_Coupled_S";
        DefPot = 1;
        Baseline[0] = "Norm";
        Baseline[1] = "Longbaseline";
        Baseline[2] = "Longbaseline2";
        DefBaseline = 2;

        FemtoRangeVars[0] = 0;//default
        FemtoRangeVars[1] = 1;//down
        FemtoRangeVars[2] = 2;//up
        DefFemtoRange = 0;

        FitRangeVars[0] = 0;//default
        FitRangeVars[1] = 1;//down
        FitRangeVars[2] = 2;//up
        DefFitRange = 0;

        pFracVars[0] = 0;//default
        pFracVars[1] = 1;
        pFracVars[2] = 2;
        DefProtonFrac = 0;

        LamFracVars[0] = 0;//default
        LamFracVars[1] = 1;
        LamFracVars[2] = 2;

        LamFracVars[3] = 10;
        LamFracVars[4] = 11;
        LamFracVars[5] = 12;

        LamFracVars[6] = 20;
        LamFracVars[7] = 21;
        LamFracVars[8] = 22;

        DefLambdaFrac = 0;

        kcVars[0] = 420;
        kcVars[1] = 580;
        kcVars[2] = 740;
        DefKc = 0;

        MtVars[0] = -1;
        DefMt = 0;

        for(int iData=0; iData<NumDataVars; iData++){
            DataVars[iData] = TString::Format("_%i",iData);
        }
        DefData = 0;
    }
    else if(WhichConfiguration==1){
        Source[0] = "McLevyNolan_Reso";
        DefSource = 0;
        SourceScale[0] = -1;//free fit, but within the 15% around 0.92
        SourceScale[1] = 0.92*0.90;
        SourceScale[2] = 0.92;
        SourceScale[3] = 0.92*1.10;
        DefSourceScale = 0;
        SourceStability[0] = -1;
        SourceStability[1] = 1.5*0.90;
        SourceStability[2] = 1.5*1.10;
        DefSourceStability = 0;
        Potential[0] = "LO_Coupled_S";
        Potential[1] = "NLO_Coupled_S";
        DefPot = 1;
        Baseline[0] = "Norm";
        Baseline[1] = "Longbaseline";
        Baseline[2] = "Longbaseline2";
        DefBaseline = 2;

        FemtoRangeVars[0] = 0;//default
        FemtoRangeVars[1] = 1;//down
        FemtoRangeVars[2] = 2;//up
        DefFemtoRange = 0;

        FitRangeVars[0] = 0;//default
        FitRangeVars[1] = 1;//down
        FitRangeVars[2] = 2;//up
        DefFitRange = 0;

        pFracVars[0] = 0;//default
        pFracVars[1] = 1;
        pFracVars[2] = 2;
        DefProtonFrac = 0;

        LamFracVars[0] = 0;//default
        LamFracVars[1] = 1;
        LamFracVars[2] = 2;

        LamFracVars[3] = 10;
        LamFracVars[4] = 11;
        LamFracVars[5] = 12;

        LamFracVars[6] = 20;
        LamFracVars[7] = 21;
        LamFracVars[8] = 22;
        DefLambdaFrac = 0;

        kcVars[0] = 700;
        kcVars[1] = 400;
        kcVars[2] = 1000;
        DefKc = 0;

        MtVars[0] = -1;
        DefMt = 0;
        for(int iData=0; iData<NumDataVars; iData++){
            DataVars[iData] = TString::Format("_%i",iData);
        }
        DefData = 0;
    }
    //pure gauss
    else if(WhichConfiguration==2){
        Source[0] = "Gauss";
        DefSource = 0;
        SourceScale[0] = 1.31;
        SourceScale[1] = 1.31*0.90;
        SourceScale[2] = 1.31*1.10;
        DefSourceScale = 0;
        SourceStability[0] = 0;//dummy
        DefSourceStability = 0;
        Potential[0] = "LO_Coupled_S";
        Potential[1] = "NLO_Coupled_S";
        DefPot = 1;
        Baseline[0] = "Norm";
        Baseline[1] = "Longbaseline";
        Baseline[2] = "Longbaseline2";
        DefBaseline = 2;

        FemtoRangeVars[0] = 0;//default
        FemtoRangeVars[1] = 1;//down
        FemtoRangeVars[2] = 2;//up
        DefFemtoRange = 0;

        FitRangeVars[0] = 0;//default
        FitRangeVars[1] = 1;//down
        FitRangeVars[2] = 2;//up
        DefFitRange = 0;

        pFracVars[0] = 0;//default
        pFracVars[1] = 1;
        pFracVars[2] = 2;
        DefProtonFrac = 0;

        LamFracVars[0] = 0;//default
        LamFracVars[1] = 1;
        LamFracVars[2] = 2;

        LamFracVars[3] = 10;
        LamFracVars[4] = 11;
        LamFracVars[5] = 12;

        LamFracVars[6] = 20;
        LamFracVars[7] = 21;
        LamFracVars[8] = 22;

        DefLambdaFrac = 0;

        kcVars[0] = 700;
        kcVars[1] = 400;
        kcVars[2] = 1000;

        DefKc = 0;

        MtVars[0] = -1;
        DefMt = 0;
        for(int iData=0; iData<NumDataVars; iData++){
            DataVars[iData] = TString::Format("_%i",iData);
        }
        DefData = 0;

    }
    //without the peak
    else if(WhichConfiguration==10){
        Source[0] = "McGauss_Reso";
        DefSource = 0;
        SourceScale[0] = -1;//free fit, but within the 15% around 0.87
        SourceScale[1] = 0.87*0.90;
        SourceScale[2] = 0.87;
        SourceScale[3] = 0.87*1.10;
        DefSourceScale = 0;
        SourceStability[0] = 0;//dummy
        DefSourceStability = 0;
        Potential[0] = "LO";
        Potential[1] = "NLO";
        DefPot = 1;
        Baseline[0] = "Norm";
        Baseline[1] = "Longbaseline_prefit";
        Baseline[2] = "Longbaseline2_prefit";
        DefBaseline = 2;

        FemtoRangeVars[0] = 10;//default
        FemtoRangeVars[1] = 11;//down
        FemtoRangeVars[2] = 12;//up
        DefFemtoRange = 0;

        FitRangeVars[0] = 10;//default
        FitRangeVars[1] = 11;//down
        FitRangeVars[2] = 12;//up
        DefFitRange = 0;

        pFracVars[0] = 0;//default
        pFracVars[1] = 1;
        pFracVars[2] = 2;
        DefProtonFrac = 0;

        LamFracVars[0] = 0;//default
        LamFracVars[1] = 1;
        LamFracVars[2] = 2;

        LamFracVars[3] = 10;
        LamFracVars[4] = 11;
        LamFracVars[5] = 12;

        LamFracVars[6] = 20;
        LamFracVars[7] = 21;
        LamFracVars[8] = 22;

        DefLambdaFrac = 0;

        kcVars[0] = 420;
        kcVars[1] = 580;
        kcVars[2] = 740;

        DefKc = 0;

        MtVars[0] = -1;
        DefMt = 0;
        for(int iData=0; iData<NumDataVars; iData++){
            DataVars[iData] = TString::Format("_%i",iData);
        }
        DefData = 0;

    }
    else if(WhichConfiguration==11){
        Source[0] = "McLevyNolan_Reso";
        DefSource = 0;
        SourceScale[0] = -1;//free fit, but within the 15% around 0.92
        SourceScale[1] = 0.92*0.90;
        SourceScale[2] = 0.92;
        SourceScale[3] = 0.92*1.10;
        DefSourceScale = 0;
        SourceStability[0] = -1;
        SourceStability[1] = 1.5*0.90;
        SourceStability[2] = 1.5*1.10;
        DefSourceStability = 0;
        Potential[0] = "LO";
        Potential[1] = "NLO";
        DefPot = 1;
        Baseline[0] = "Norm";
        Baseline[1] = "Longbaseline_prefit";
        Baseline[2] = "Longbaseline2_prefit";
        DefBaseline = 2;

        FemtoRangeVars[0] = 10;//default
        FemtoRangeVars[1] = 11;//down
        FemtoRangeVars[2] = 12;//up
        DefFemtoRange = 0;

        FitRangeVars[0] = 10;//default
        FitRangeVars[1] = 11;//down
        FitRangeVars[2] = 12;//up
        DefFitRange = 0;

        pFracVars[0] = 0;//default
        pFracVars[1] = 1;
        pFracVars[2] = 2;
        DefProtonFrac = 0;

        LamFracVars[0] = 0;//default
        LamFracVars[1] = 1;
        LamFracVars[2] = 2;

        LamFracVars[3] = 10;
        LamFracVars[4] = 11;
        LamFracVars[5] = 12;

        LamFracVars[6] = 20;
        LamFracVars[7] = 21;
        LamFracVars[8] = 22;
        DefLambdaFrac = 0;

        kcVars[0] = 700;
        kcVars[1] = 400;
        kcVars[2] = 1000;
        DefKc = 0;

        MtVars[0] = -1;
        DefMt = 0;
        for(int iData=0; iData<NumDataVars; iData++){
            DataVars[iData] = TString::Format("_%i",iData);
        }
        DefData = 0;

    }
    //pure gauss
    else if(WhichConfiguration==12){
        Source[0] = "Gauss";
        DefSource = 0;
        SourceScale[0] = 1.31;
        SourceScale[1] = 1.31*0.90;
        SourceScale[2] = 1.31*1.10;
        DefSourceScale = 0;
        SourceStability[0] = 0;//dummy
        DefSourceStability = 0;
        Potential[0] = "LO";
        Potential[1] = "NLO";
        DefPot = 1;
        Baseline[0] = "Norm";
        Baseline[1] = "Longbaseline_prefit";
        Baseline[2] = "Longbaseline2_prefit";
        DefBaseline = 2;

        FemtoRangeVars[0] = 10;//default
        FemtoRangeVars[1] = 11;//down
        FemtoRangeVars[2] = 12;//up
        DefFemtoRange = 0;

        FitRangeVars[0] = 10;//default
        FitRangeVars[1] = 11;//down
        FitRangeVars[2] = 12;//up
        DefFitRange = 0;

        pFracVars[0] = 0;//default
        pFracVars[1] = 1;
        pFracVars[2] = 2;
        DefProtonFrac = 0;

        LamFracVars[0] = 0;//default
        LamFracVars[1] = 1;
        LamFracVars[2] = 2;

        LamFracVars[3] = 10;
        LamFracVars[4] = 11;
        LamFracVars[5] = 12;

        LamFracVars[6] = 20;
        LamFracVars[7] = 21;
        LamFracVars[8] = 22;

        DefLambdaFrac = 0;

        kcVars[0] = 700;
        kcVars[1] = 400;
        kcVars[2] = 1000;
        DefKc = 0;

        MtVars[0] = -1;
        DefMt = 0;
        for(int iData=0; iData<NumDataVars; iData++){
            DataVars[iData] = TString::Format("_%i",iData);
        }
        DefData = 0;
    }
    else if(WhichConfiguration==20||WhichConfiguration==21||WhichConfiguration==22){
        if(WhichConfiguration==20) Source[0] = "McGauss_Reso";
        else if(WhichConfiguration==21) Source[0] = "McLevyNolan_Reso";
        else Source[0] = "Gauss";

        DefSource = 0;
        SourceScale[0] = -11;//free fit
        DefSourceScale = 0;
        SourceStability[0] = WhichConfiguration==21?-11:0;//dummy
        DefSourceStability = 0;
        Potential[0] = "LO_Coupled_S";
        Potential[1] = "NLO_Coupled_S";
        DefPot = 1;
        Baseline[0] = "Norm";
        Baseline[1] = "Longbaseline";
        Baseline[2] = "Longbaseline2";
        DefBaseline = 2;

        FemtoRangeVars[0] = 0;//default
        FemtoRangeVars[1] = 1;//down
        FemtoRangeVars[2] = 2;//up
        DefFemtoRange = 0;

        FitRangeVars[0] = 0;//default
        FitRangeVars[1] = 1;//down
        FitRangeVars[2] = 2;//up
        DefFitRange = 0;

        pFracVars[0] = 0;//default
        pFracVars[1] = 1;
        pFracVars[2] = 2;
        DefProtonFrac = 0;

        LamFracVars[0] = 0;//default
        LamFracVars[1] = 1;
        LamFracVars[2] = 2;

        LamFracVars[3] = 10;
        LamFracVars[4] = 11;
        LamFracVars[5] = 12;

        LamFracVars[6] = 20;
        LamFracVars[7] = 21;
        LamFracVars[8] = 22;

        DefLambdaFrac = 0;

        kcVars[0] = 420;
        kcVars[1] = 580;
        kcVars[2] = 740;
        DefKc = 0;

        MtVars[0] = -1;
        MtVars[1] = 0;
        MtVars[2] = 1;
        MtVars[3] = 2;
        MtVars[4] = 3;
        MtVars[5] = 4;
        MtVars[6] = 5;
        DefMt = 0;
        for(int iData=0; iData<NumDataVars; iData++){
            DataVars[iData] = TString::Format("_%i",iData);
        }
        DefData = 0;
    }
    else if(WhichConfiguration==30||WhichConfiguration==31||WhichConfiguration==32){
        if(WhichConfiguration==30) Source[0] = "McGauss_Reso";
        else if(WhichConfiguration==31) Source[0] = "McLevyNolan_Reso";
        else Source[0] = "Gauss";

        DefSource = 0;
        SourceScale[0] = -11;//free fit
        DefSourceScale = 0;
        SourceStability[0] = WhichConfiguration==31?-11:0;//dummy
        DefSourceStability = 0;
        Potential[0] = "LO_Coupled_S";
        Potential[1] = "NLO_Coupled_S";
        Potential[2] = "LO";
        Potential[3] = "NLO";
        Potential[4] = "Usmani";
        DefPot = 1;
        Baseline[0] = "Norm";
        Baseline[1] = "Longbaseline";
        Baseline[2] = "Longbaseline2";
        DefBaseline = 2;

        FemtoRangeVars[0] = 10;//default
        FemtoRangeVars[1] = 11;//down
        FemtoRangeVars[2] = 12;//up
        DefFemtoRange = 0;

        FitRangeVars[0] = 10;//default
        FitRangeVars[1] = 11;//down
        FitRangeVars[2] = 12;//up
        DefFitRange = 0;

        pFracVars[0] = 0;//default
        pFracVars[1] = 1;
        pFracVars[2] = 2;
        DefProtonFrac = 0;

        LamFracVars[0] = 0;//default
        LamFracVars[1] = 1;
        LamFracVars[2] = 2;

        LamFracVars[3] = 10;
        LamFracVars[4] = 11;
        LamFracVars[5] = 12;

        LamFracVars[6] = 20;
        LamFracVars[7] = 21;
        LamFracVars[8] = 22;

        DefLambdaFrac = 0;

        kcVars[0] = -1;
        DefKc = 0;

        MtVars[0] = -1;
        MtVars[1] = 0;
        MtVars[2] = 1;
        MtVars[3] = 2;
        MtVars[4] = 3;
        MtVars[5] = 4;
        MtVars[6] = 5;
        DefMt = 0;
        for(int iData=0; iData<NumDataVars; iData++){
            DataVars[iData] = TString::Format("_%i",iData);
        }
        DefData = 0;
    }
    else if(WhichConfiguration==35){
        Source[0] = "McGauss_Reso";
        DefSource = 0;
        SourceScale[0] = -11;//free fit
        DefSourceScale = 0;
        SourceStability[0] = 0;//dummy
        DefSourceStability = 0;
        Potential[0] = "LO";
        Potential[1] = "NLO";
        Potential[2] = "Usmani";
        DefPot = 1;
        Baseline[0] = "Norm";
        Baseline[1] = "Longbaseline2";
        DefBaseline = 2;

        FemtoRangeVars[0] = 10;//default
        FemtoRangeVars[1] = 11;//down
        FemtoRangeVars[2] = 12;//up
        DefFemtoRange = 0;

        FitRangeVars[0] = 10;//default
        FitRangeVars[1] = 11;//down
        FitRangeVars[2] = 12;//up
        DefFitRange = 0;

        pFracVars[0] = 0;//default
        pFracVars[1] = 1;
        pFracVars[2] = 2;
        DefProtonFrac = 0;

        LamFracVars[0] = 0;//default
        LamFracVars[1] = 1;
        LamFracVars[2] = 2;

        LamFracVars[3] = 10;
        LamFracVars[4] = 11;
        LamFracVars[5] = 12;

        LamFracVars[6] = 20;
        LamFracVars[7] = 21;
        LamFracVars[8] = 22;

        DefLambdaFrac = 0;

        kcVars[0] = -1;
        DefKc = 0;

        MtVars[0] = 0;
        DefMt = 0;
        for(int iData=0; iData<NumDataVars; iData++){
            DataVars[iData] = TString::Format("_%i",iData);
        }
        DefData = 0;
    }
    else if(WhichConfiguration==40||WhichConfiguration==41||WhichConfiguration==42){
        if(WhichConfiguration==40) Source[0] = "McGauss_Reso";
        else if(WhichConfiguration==41) Source[0] = "McLevyNolan_Reso";
        else Source[0] = "Gauss";

        DefSource = 0;
        SourceScale[0] = -11;//free fit
        DefSourceScale = 0;
        SourceStability[0] = 0;//dummy
        if(WhichConfiguration==41) SourceStability[0] = -11;
        DefSourceStability = 0;
        Potential[0] = "NLO_Coupled_S";
        if(NumPotVars>1) Potential[1] = "LO_Coupled_S";
        DefPot = 0;
        Baseline[0] = "Spline3_Tune1";
        //Baseline[0] = "Spline3_3";
        //Baseline[1] = "Spline3_4";
        //Baseline[2] = "Spline3_5";
        //Baseline[3] = "Spline3_7";
        //Baseline[4] = "Spline3_9";
        //Baseline[5] = "Spline3_11";
        DefBaseline = 0;

        FemtoRangeVars[0] = 0;//default
        FemtoRangeVars[1] = 1;//down
        FemtoRangeVars[2] = 2;//up
        DefFemtoRange = 0;

        FitRangeVars[0] = 0;//default
        FitRangeVars[1] = 1;//down
        FitRangeVars[2] = 2;//up
        DefFitRange = 0;

        pFracVars[0] = 0;//default
        pFracVars[1] = 1;
        pFracVars[2] = 2;
        DefProtonFrac = 0;

        LamFracVars[0] = 0;//default
        LamFracVars[1] = 1;
        LamFracVars[2] = 2;

        LamFracVars[3] = 10;
        LamFracVars[4] = 11;
        LamFracVars[5] = 12;

        LamFracVars[6] = 20;
        LamFracVars[7] = 21;
        LamFracVars[8] = 22;

        DefLambdaFrac = 0;

        kcVars[0] = 420;
        kcVars[1] = 580;
        kcVars[2] = 740;
        DefKc = 0;

        MtVars[0] = -1;
        DefMt = 0;
        DataVars[0] = TString::Format("_%i",0);
        DefData = 0;
    }
    //meson exchange models with Lednicky, without the peak region
    else if(WhichConfiguration==112){
        Source[0] = "Gauss";
        DefSource = 0;
        SourceScale[0] = 1.31;
        SourceScale[1] = 1.31*0.90;
        SourceScale[2] = 1.31*1.10;
        DefSourceScale = 0;
        SourceStability[0] = 0;//dummy
        DefSourceStability = 0;
        Potential[0] = "Lednicky_ND";
        Potential[1] = "Lednicky_NF";
        Potential[2] = "Lednicky_NSC89";
        Potential[3] = "Lednicky_NSC97a";
        Potential[4] = "Lednicky_NSC97b";
        Potential[5] = "Lednicky_NSC97c";
        Potential[6] = "Lednicky_NSC97d";
        Potential[7] = "Lednicky_NSC97e";
        Potential[8] = "Lednicky_NSC97f";
        Potential[9] = "Lednicky_ESC08";
        Potential[10] = "Lednicky_XeftLO";
        Potential[11] = "Lednicky_XeftNLO";
        Potential[12] = "Lednicky_JulichA";
        Potential[13] = "Lednicky_JulichJ04";
        Potential[14] = "Lednicky_JulichJ04c";
        DefPot = 11;

        Baseline[0] = "Norm";
        Baseline[1] = "Longbaseline";
        Baseline[2] = "Longbaseline2";
        DefBaseline = 2;

        FemtoRangeVars[0] = 10;//default
        FemtoRangeVars[1] = 11;//down
        FemtoRangeVars[2] = 12;//up
        DefFemtoRange = 0;

        FitRangeVars[0] = 10;//default
        FitRangeVars[1] = 11;//down
        FitRangeVars[2] = 12;//up
        DefFitRange = 0;

        pFracVars[0] = 0;//default
        pFracVars[1] = 1;
        pFracVars[2] = 2;
        DefProtonFrac = 0;

        LamFracVars[0] = 0;//default
        LamFracVars[1] = 1;
        LamFracVars[2] = 2;

        LamFracVars[3] = 10;
        LamFracVars[4] = 11;
        LamFracVars[5] = 12;

        LamFracVars[6] = 20;
        LamFracVars[7] = 21;
        LamFracVars[8] = 22;

        DefLambdaFrac = 0;

        kcVars[0] = 1e6;
        DefKc = 0;

        MtVars[0] = -1;
        DefMt = 0;
        DataVars[0] = "";
        DefData = 0;
    }
    else if(WhichConfiguration==200){
        Source[0] = "McGauss_Reso";
        DefSource = 0;
        SourceScale[0] = -11;
        DefSourceScale = 0;
        SourceStability[0] = 0;//dummy
        DefSourceStability = 0;
        Potential[0] = "LO_Coupled_S";
        Potential[1] = "NLO_Coupled_S";
        DefPot = 1;
        Baseline[0] = "Norm";
        Baseline[1] = "Longbaseline";
        Baseline[2] = "Longbaseline2";
        DefBaseline = 2;

        FemtoRangeVars[0] = 0;//default
        FemtoRangeVars[1] = 1;//down
        FemtoRangeVars[2] = 2;//up
        DefFemtoRange = 0;

        FitRangeVars[0] = 0;//default
        FitRangeVars[1] = 1;//down
        FitRangeVars[2] = 2;//up
        DefFitRange = 0;

        pFracVars[0] = 0;//default
        pFracVars[1] = 1;
        pFracVars[2] = 2;
        DefProtonFrac = 0;

        LamFracVars[0] = 0;//default
        LamFracVars[1] = 1;
        LamFracVars[2] = 2;

        LamFracVars[3] = 10;
        LamFracVars[4] = 11;
        LamFracVars[5] = 12;

        LamFracVars[6] = 20;
        LamFracVars[7] = 21;
        LamFracVars[8] = 22;

        DefLambdaFrac = 0;

        kcVars[0] = 420;
        kcVars[1] = 580;
        kcVars[2] = 740;
        DefKc = 0;

        MtVars[0] = -1;
        DefMt = 0;

        for(int iData=0; iData<NumDataVars; iData++){
            DataVars[iData] = TString::Format("_%i",iData);
        }
        DefData = 0;
    }
    else if(WhichConfiguration==201){
        Source[0] = "McLevyNolan_Reso";
        DefSource = 0;
        SourceScale[0] = -11;
        DefSourceScale = 0;
        SourceStability[0] = -1;
        SourceStability[1] = 1.5*0.90;
        SourceStability[2] = 1.5*1.10;
        DefSourceStability = 0;
        Potential[0] = "LO_Coupled_S";
        Potential[1] = "NLO_Coupled_S";
        DefPot = 1;
        Baseline[0] = "Norm";
        Baseline[1] = "Longbaseline";
        Baseline[2] = "Longbaseline2";
        DefBaseline = 2;

        FemtoRangeVars[0] = 0;//default
        FemtoRangeVars[1] = 1;//down
        FemtoRangeVars[2] = 2;//up
        DefFemtoRange = 0;

        FitRangeVars[0] = 0;//default
        FitRangeVars[1] = 1;//down
        FitRangeVars[2] = 2;//up
        DefFitRange = 0;

        pFracVars[0] = 0;//default
        pFracVars[1] = 1;
        pFracVars[2] = 2;
        DefProtonFrac = 0;

        LamFracVars[0] = 0;//default
        LamFracVars[1] = 1;
        LamFracVars[2] = 2;

        LamFracVars[3] = 10;
        LamFracVars[4] = 11;
        LamFracVars[5] = 12;

        LamFracVars[6] = 20;
        LamFracVars[7] = 21;
        LamFracVars[8] = 22;
        DefLambdaFrac = 0;

        kcVars[0] = 700;
        kcVars[1] = 400;
        kcVars[2] = 1000;
        DefKc = 0;

        MtVars[0] = -1;
        DefMt = 0;
        for(int iData=0; iData<NumDataVars; iData++){
            DataVars[iData] = TString::Format("_%i",iData);
        }
        DefData = 0;
    }
    //pure gauss
    else if(WhichConfiguration==202){
        Source[0] = "Gauss";
        DefSource = 0;
        SourceScale[0] = -11;
        DefSourceScale = 0;
        SourceStability[0] = 0;//dummy
        DefSourceStability = 0;
        Potential[0] = "LO_Coupled_S";
        Potential[1] = "NLO_Coupled_S";
        DefPot = 1;
        Baseline[0] = "Norm";
        Baseline[1] = "Longbaseline";
        Baseline[2] = "Longbaseline2";
        DefBaseline = 2;

        FemtoRangeVars[0] = 0;//default
        FemtoRangeVars[1] = 1;//down
        FemtoRangeVars[2] = 2;//up
        DefFemtoRange = 0;

        FitRangeVars[0] = 0;//default
        FitRangeVars[1] = 1;//down
        FitRangeVars[2] = 2;//up
        DefFitRange = 0;

        pFracVars[0] = 0;//default
        pFracVars[1] = 1;
        pFracVars[2] = 2;
        DefProtonFrac = 0;

        LamFracVars[0] = 0;//default
        LamFracVars[1] = 1;
        LamFracVars[2] = 2;

        LamFracVars[3] = 10;
        LamFracVars[4] = 11;
        LamFracVars[5] = 12;

        LamFracVars[6] = 20;
        LamFracVars[7] = 21;
        LamFracVars[8] = 22;

        DefLambdaFrac = 0;

        kcVars[0] = 700;
        kcVars[1] = 400;
        kcVars[2] = 1000;

        DefKc = 0;

        MtVars[0] = -1;
        DefMt = 0;
        for(int iData=0; iData<NumDataVars; iData++){
            DataVars[iData] = TString::Format("_%i",iData);
        }
        DefData = 0;

    }
    else{
        printf("????????????\n");
    }

    //if(!JustNumIter){
    //    gROOT->cd();
    //    ntResult = new TNtuple("ntResult", "ntResult",
    //        "SourceType:SourceScale:SourceStability:Potential:Baseline:"
    //        "p_a:e_a:p_b:e_b:p_c:e_c:p_Cl:e_Cl:p_kc:e_kc:p_sor0:e_sor0:p_sor1:e_sor1:chi2:ndf");
    //}

    TRandom3 RanGen(RANDOMSEED);

    TString ntFileName = TString::Format("%sNTfile_%i.root",OutputFolder.Data(),WhichConfiguration);


    int IterationCounter = 0;
    int ComputedIterations = 0;
    for(int iSource=0; iSource<NumSourceVars; iSource++){
    for(int iSourceScal=0; iSourceScal<NumSourceScaleVars; iSourceScal++){
    for(int iSourceStab=0; iSourceStab<NumSourceStabilityVars; iSourceStab++){
    for(int iPot=0; iPot<NumPotVars; iPot++){
    for(int iBl=0; iBl<NumBaselineVars; iBl++){
    for(int iFemRan=0; iFemRan<NumFemtoRangeVars; iFemRan++){
    for(int iFitRan=0; iFitRan<NumFitRangeVars; iFitRan++){
    for(int ipFrac=0; ipFrac<NumProtonFracVars; ipFrac++){
    for(int iLamFrac=0; iLamFrac<NumLambdaFracVars; iLamFrac++){
    for(int ikc=0; ikc<NumKcVars; ikc++){
    for(int iMt=0; iMt<NumMtVars; iMt++){
    for(int iData=0; iData<NumDataVars; iData++){
        if(Baseline[iBl]=="Norm"&&ikc) {continue;}//this variation is only relevant in case of a baseline
        if(Baseline[iBl]=="Norm"&&iFitRan) {continue;}//this variation is only relevant in case of a baseline
        double UNI = RanGen.Uniform();
        //if(IterationCounter>=FirstIter&&IterationCounter<=LastIter){
        //    printf(" UNI=%.3f vs %.3f\n",UNI,1./double(OnlyFraction));
        //    printf(" IterationCounter=%i vs %i\n",IterationCounter,ComputedIterations);
        //}
        bool DefaultVar=true;
        if(iSourceScal!=DefSourceScale||iSourceStab!=DefSourceStability||
                    iFemRan!=DefFemtoRange||iFitRan!=DefFitRange||
                    ipFrac!=DefProtonFrac||iLamFrac!=DefLambdaFrac||ikc!=DefKc||iData!=DefData){DefaultVar=false;}
//if(iPot==4&&DefaultVar==true) printf("HELLO!\n")
        if(DefaultVar) {DEFAULT_VAR[iSource][iPot] = IterationCounter; UNI=-1;}//we always compute the default variations (for all potentials) and for all baselines
        if(Mode==-1){
            if(DefaultVar){IterationCounter++; continue;}
            else{
                DEFAULT_VAR[iSource][iPot] = IterationCounter;
            }
        }
        else if(Mode==1){
            if(DefaultVar){IterationCounter++; continue;}
            else{
                DEFAULT_VAR[iSource][iPot] = IterationCounter;
                IterationCounter++;
                continue;
            }
        }
        else{
            if(Mode!=-1&&(IterationCounter<FirstIter||IterationCounter>LastIter)) {IterationCounter++; continue;}
            unsigned EffectiveFraction = OnlyFraction;
            if(Baseline[iBl]=="Norm"){
                EffectiveFraction /= NumFemtoRangeVars;
                if(!Potential[iPot].Contains("Lednicky_")) EffectiveFraction /= NumFitRangeVars;
                if(EffectiveFraction<0||EffectiveFraction>OnlyFraction) EffectiveFraction=1;
            }
            //this is here in order to make sure we take more Norm events (now they are much less due to some conditions above)
            if(EffectiveFraction>1&&UNI>1./double(EffectiveFraction)) {IterationCounter++; continue;}
        }


        //TString IterationName = TString::Format("%spLam_%s_%s_Sca%.2f_Sta%.2f_%s_%s_",
        //                                        OutputFolder.Data(),"pp13TeV_HM_March19",Source[iSource].Data(),
        //                                        SourceScale[iSourceScal],SourceStability[iSourceStab],Potential[iPot],Baseline[iBl],
        //                                        );

        //[0] = femto range
        //[1] = fit range
        //[2] = fractions for p
        //[3] = fractions for L
        VARIATIONS[0] = FemtoRangeVars[iFemRan];
        VARIATIONS[1] = FitRangeVars[iFitRan];
        VARIATIONS[2] = pFracVars[ipFrac];
        VARIATIONS[3] = LamFracVars[iLamFrac];
        VARIATIONS[4] = kcVars[ikc];
        VARIATIONS[5] = MtVars[iMt];

        TString SourceDescription = Source[iSource];
        if(SourceScale[iSourceScal]<0) {SourceDescription+="ScaFree";}
        else if(SourceScale[iSourceScal]==0) {SourceDescription+="ScaDefault";}
        else {SourceDescription+=TString::Format("Sca%.2f",SourceScale[iSourceScal]);}
        if(SourceStability[iSourceStab]<0||SourceStability[iSourceStab]>2) {SourceDescription+="StaFree";}
        else if(SourceStability[iSourceStab]==0) {SourceDescription+="StaDefault";}
        else {SourceDescription+=TString::Format("Sta%.2f",SourceStability[iSourceStab]);}
        /*
        TString FileName = TString::Format("Iter%i_%s_%s_%s_%s_%i_%i_%i_%i_%i.root",IterationCounter,
                                                               "pp13TeV_HM_March19",SourceDescription.Data(),
                                                               Potential[iPot].Data(),Baseline[iBl].Data(),
                                                               VARIATIONS[0],VARIATIONS[1],VARIATIONS[2],VARIATIONS[3],VARIATIONS[4]);
        */
        TString FileName = TString::Format("Config%i_Iter%i.root",WhichConfiguration,IterationCounter);
        TFile* CheckIfFile = new TFile(OutputFolder+FileName,"READ");
        //! we ignore iterations for which we already have an input file!
        //! => delete the folder in case you want to overwrite
        if(CheckIfFile&&CheckIfFile->IsOpen()){IterationCounter++; continue;}
        delete CheckIfFile;
        CheckIfFile = NULL;

        Fit_pL(AnalysisObject, OutputFolder, FileName,"pp13TeV_HM_March19",DataVars[iData],Source[iSource],SourceScale[iSourceScal],SourceStability[iSourceStab],
               Potential[iPot],Baseline[iBl],VARIATIONS,long(IterationCounter)*long(1000)+long(WhichConfiguration),ntFileName,"ntResult");

        IterationCounter++;
        ComputedIterations++;
    }
    }
    }
    }
    }
    }
    }
    }
    }
    }
    }
    }

    printf("This set up had a total of %i iterations.\n",IterationCounter);

     for(int iSource=0; iSource<NumSourceVars; iSource++) for(int iPot=0; iPot<NumPotVars; iPot++)
        printf(" The default variation (%i) is %i.\n",iPot,DEFAULT_VAR[iSource][iPot]);
    printf(" The function call evaluated %i iterations.\n",ComputedIterations);

    //if(NTfile){
    //    NTfile->cd();
    //    ntResult->Write();
    //}


    //if(ntResult) delete ntResult;
    //if(NTfile) delete NTfile;

    delete [] Source;
    delete [] SourceScale;
    delete [] SourceStability;
    delete [] Potential;
    delete [] Baseline;
    delete [] FemtoRangeVars;
    delete [] FitRangeVars;
    delete [] pFracVars;
    delete [] LamFracVars;
    delete [] kcVars;
    delete [] MtVars;
}


//if OnlyFraction>1 => we have 1/OnlyFraction probability to accept an event
//Mode==-1 => the default cut comb only (ignoring all first, last iter etc.)
//Mode==0 => normal
//Mode==1 => we only print the number of iter to be performed
void pp_SystematicsHM(const TString& OutputFolder, const int& WhichConfiguration, const int& FirstIter, const int& LastIter, const unsigned& OnlyFraction,
                      const int& RANDOMSEED, const int& Mode){

    //const bool Perform_mT_Scan = true;

    printf("Output: %s\n",OutputFolder.Data());
    printf("Iter: %i -> %i\n",FirstIter,LastIter);
    printf("Fraction: %i\n",OnlyFraction);
    printf("Seed: %i\n",RANDOMSEED);
    DLM_CommonAnaFunctions AnalysisObject;
    int NumSourceVars;
    int NumSourceScaleVars;
    int NumSourceStabilityVars;
    int NumPotVars;
    int NumBaselineVars;
    int NumFemtoRangeVars;
    int NumFitRangeVars;
    int NumProtonFracVars;
    int NumLambdaFracVars;
    int NumKcVars;
    int NumMtVars;
    int NumDataVars;//these are the cut variations

    int DefSource;
    int DefSourceScale;
    int DefSourceStability;
    int DefPot;
    int DefBaseline;
    int DefFemtoRange;
    int DefFitRange;
    int DefProtonFrac;
    int DefLambdaFrac;
    int DefKc;
    int DefMt;
    int DefData;

    TString* Source;
    double* SourceScale;
    double* SourceStability;
    TString* Potential;
    TString* Baseline;
    int* FemtoRangeVars;
    int* FitRangeVars;
    int* pFracVars;
    int* LamFracVars;
    double* kcVars;
    int* MtVars;
    TString* DataVars;
    int VARIATIONS[6];
    int TotNumIter=0;

    switch(WhichConfiguration){
        //this is the case where we only check for the Gauss+Reso and Gauss
        case 0 :
            NumSourceVars = 1;
            NumSourceScaleVars = 1;
            NumSourceStabilityVars = 1;
            NumPotVars = 1;
            NumBaselineVars = 5;
            NumFemtoRangeVars = 3;
            NumFitRangeVars = 3;
            NumProtonFracVars = 3;
            NumLambdaFracVars = 1;
            NumKcVars = 3;
            NumMtVars = 1;
            NumDataVars = 45;
            break;
        //similar to 0, only that now we include Levy+Reso as an option
        case 1 :
            NumSourceVars = 1;
            NumSourceScaleVars = 1;
            NumSourceStabilityVars = 1;
            NumPotVars = 1;
            NumBaselineVars = 5;
            NumFemtoRangeVars = 3;
            NumFitRangeVars = 3;
            NumProtonFracVars = 3;
            NumLambdaFracVars = 1;
            NumKcVars = 3;
            NumMtVars = 1;
            NumDataVars = 45;
            break;
        case 2 :
            NumSourceVars = 1;
            NumSourceScaleVars = 1;
            NumSourceStabilityVars = 1;
            NumPotVars = 1;
            NumBaselineVars = 5;
            NumFemtoRangeVars = 3;
            NumFitRangeVars = 3;
            NumProtonFracVars = 3;
            NumLambdaFracVars = 1;
            NumKcVars = 3;
            NumMtVars = 1;
            NumDataVars = 45;
            break;
        //this is the case where we only check for the Gauss+Reso and Gauss
        //we fit with mT bins
        case 20 :
            NumSourceVars = 1;
            NumSourceScaleVars = 1;
            NumSourceStabilityVars = 1;
            NumPotVars = 1;
            NumBaselineVars = 5;
            NumFemtoRangeVars = 3;
            NumFitRangeVars = 3;
            NumProtonFracVars = 3;
            NumLambdaFracVars = 1;
            NumKcVars = 3;
            NumMtVars = 8;
            NumDataVars = 1;//45
            break;
        case 21 :
            NumSourceVars = 1;
            NumSourceScaleVars = 1;
            NumSourceStabilityVars = 1;
            NumPotVars = 1;
            NumBaselineVars = 5;
            NumFemtoRangeVars = 3;
            NumFitRangeVars = 3;
            NumProtonFracVars = 3;
            NumLambdaFracVars = 1;
            NumKcVars = 3;
            NumMtVars = 8;
            NumDataVars = 1;//45
            break;
        case 22 :
            NumSourceVars = 1;
            NumSourceScaleVars = 1;
            NumSourceStabilityVars = 1;
            NumPotVars = 1;
            NumBaselineVars = 5;
            NumFemtoRangeVars = 3;
            NumFitRangeVars = 3;
            NumProtonFracVars = 3;
            NumLambdaFracVars = 1;
            NumKcVars = 3;
            NumMtVars = 8;
            NumDataVars = 1;//45
            break;
        //the test fro splines
        case 42 :
            NumSourceVars = 1;
            NumSourceScaleVars = 1;
            NumSourceStabilityVars = 1;
            NumPotVars = 1;
            NumBaselineVars = 1;
            NumFemtoRangeVars = 3;
            NumFitRangeVars = 3;
            NumProtonFracVars = 3;
            NumLambdaFracVars = 1;
            NumKcVars = 3;
            NumMtVars = 1;
            NumDataVars = 1;
            break;
        //this is the case where we only check for the Gauss+Reso and Gauss, default data
        case 10 :
            NumSourceVars = 1;
            NumSourceScaleVars = 1;
            NumSourceStabilityVars = 1;
            NumPotVars = 1;
            NumBaselineVars = 5;
            NumFemtoRangeVars = 3;
            NumFitRangeVars = 3;
            NumProtonFracVars = 3;
            NumLambdaFracVars = 1;
            NumKcVars = 3;
            NumMtVars = 1;
            NumDataVars = 1;
            break;
        //similar to 0, only that now we include Levy+Reso as an option
        case 11 :
            NumSourceVars = 1;
            NumSourceScaleVars = 1;
            NumSourceStabilityVars = 1;
            NumPotVars = 1;
            NumBaselineVars = 5;
            NumFemtoRangeVars = 3;
            NumFitRangeVars = 3;
            NumProtonFracVars = 3;
            NumLambdaFracVars = 1;
            NumKcVars = 3;
            NumMtVars = 1;
            NumDataVars = 1;
            break;
        case 12:
            NumSourceVars = 1;
            NumSourceScaleVars = 1;
            NumSourceStabilityVars = 1;
            NumPotVars = 1;
            NumBaselineVars = 5;
            NumFemtoRangeVars = 3;
            NumFitRangeVars = 3;
            NumProtonFracVars = 3;
            NumLambdaFracVars = 1;
            NumKcVars = 3;
            NumMtVars = 1;
            NumDataVars = 1;
            break;
        default : return;
    }


    Source = new TString [NumSourceVars];
    SourceScale = new double [NumSourceScaleVars];
    SourceStability = new double [NumSourceStabilityVars];
    Potential = new TString [NumPotVars];
    Baseline = new TString [NumBaselineVars];
    FemtoRangeVars = new int [NumFemtoRangeVars];
    FitRangeVars = new int [NumFitRangeVars];
    pFracVars = new int [NumProtonFracVars];
    LamFracVars = new int [NumLambdaFracVars];
    kcVars = new double [NumLambdaFracVars];
    MtVars = new int [NumMtVars];
    DataVars = new TString [NumDataVars];

    int** DEFAULT_VAR = new int* [NumSourceVars];
    for(int iSource=0; iSource<NumSourceVars; iSource++){
        DEFAULT_VAR[iSource]=new int [NumPotVars];
        for(int iPot=0; iPot<NumPotVars; iPot++)DEFAULT_VAR[iSource][iPot]=-1;
    }

    //TFile* NTfile = NULL;
    //TNtuple* ntResult = NULL;

    if(WhichConfiguration==0){
        Source[0] = "McGauss_Reso";
        DefSource = 0;
        SourceScale[0] = -11;
        DefSourceScale = 0;
        SourceStability[0] = 2.0;//dummy
        DefSourceStability = 0;
        Potential[0] = "AV18";
        DefPot = 0;
        Baseline[0] = "Norm";
        Baseline[1] = "Baseline";
        Baseline[2] = "Baseline2";
        Baseline[3] = "Longbaseline";
        Baseline[4] = "Longbaseline2";
        DefBaseline = 0;

        FemtoRangeVars[0] = 0;//default
        FemtoRangeVars[1] = 1;//down
        FemtoRangeVars[2] = 2;//up
        DefFemtoRange = 0;

        FitRangeVars[0] = 0;//default
        FitRangeVars[1] = 1;//down
        FitRangeVars[2] = 2;//up
        DefFitRange = 0;

        pFracVars[0] = 0;//default
        pFracVars[1] = 1;
        pFracVars[2] = 2;
        DefProtonFrac = 0;

        LamFracVars[0] = 0;//default

        DefLambdaFrac = 0;

        kcVars[0] = 500;
        kcVars[1] = 750;
        kcVars[2] = 1000;
        DefKc = 0;

        MtVars[0] = -1;
        DefMt = 0;

        for(int iData=0; iData<NumDataVars; iData++){
            DataVars[iData] = TString::Format("_%i",iData);
        }
        DefData = 0;
    }
    else if(WhichConfiguration==1){
        Source[0] = "McLevyNolan_Reso";
        DefSource = 0;
        SourceScale[0] = -11;
        DefSourceScale = 0;
        SourceStability[0] = -11;
        DefSourceStability = 0;
        Potential[0] = "AV18";
        DefPot = 0;
        Baseline[0] = "Norm";
        Baseline[1] = "Baseline";
        Baseline[2] = "Baseline2";
        Baseline[3] = "Longbaseline";
        Baseline[4] = "Longbaseline2";
        DefBaseline = 0;

        FemtoRangeVars[0] = 0;//default
        FemtoRangeVars[1] = 1;//down
        FemtoRangeVars[2] = 2;//up
        DefFemtoRange = 0;

        FitRangeVars[0] = 0;//default
        FitRangeVars[1] = 1;//down
        FitRangeVars[2] = 2;//up
        DefFitRange = 0;

        pFracVars[0] = 0;//default
        pFracVars[1] = 1;
        pFracVars[2] = 2;
        DefProtonFrac = 0;

        LamFracVars[0] = 0;//default

        DefLambdaFrac = 0;

        kcVars[0] = 500;
        kcVars[1] = 750;
        kcVars[2] = 1000;
        DefKc = 0;

        MtVars[0] = -1;
        DefMt = 0;

         for(int iData=0; iData<NumDataVars; iData++){
            DataVars[iData] = TString::Format("_%i",iData);
        }
        DefData = 0;
    }
    //pure gauss
    else if(WhichConfiguration==2){
        Source[0] = "Gauss";
        DefSource = 0;
        SourceScale[0] = -11;
        DefSourceScale = 0;
        SourceStability[0] = 2.0;//dummy
        DefSourceStability = 0;
        Potential[0] = "AV18";
        DefPot = 0;
        Baseline[0] = "Norm";
        Baseline[1] = "Baseline";
        Baseline[2] = "Baseline2";
        Baseline[3] = "Longbaseline";
        Baseline[4] = "Longbaseline2";
        DefBaseline = 0;

        FemtoRangeVars[0] = 0;//default
        FemtoRangeVars[1] = 1;//down
        FemtoRangeVars[2] = 2;//up
        DefFemtoRange = 0;

        FitRangeVars[0] = 0;//default
        FitRangeVars[1] = 1;//down
        FitRangeVars[2] = 2;//up
        DefFitRange = 0;

        pFracVars[0] = 0;//default
        pFracVars[1] = 1;
        pFracVars[2] = 2;
        DefProtonFrac = 0;

        LamFracVars[0] = 0;//default

        DefLambdaFrac = 0;

        kcVars[0] = 500;
        kcVars[1] = 750;
        kcVars[2] = 1000;

        DefKc = 0;

        MtVars[0] = -1;
        DefMt = 0;

        for(int iData=0; iData<NumDataVars; iData++){
            DataVars[iData] = TString::Format("_%i",iData);
        }
        DefData = 0;

    }
    else if(WhichConfiguration==20||WhichConfiguration==21||WhichConfiguration==22){
        if(WhichConfiguration==20) Source[0] = "McGauss_Reso";
        else if(WhichConfiguration==21) Source[0] = "McLevyNolan_Reso";
        else Source[0] = "Gauss";

        DefSource = 0;
        SourceScale[0] = -11;//free fit
        DefSourceScale = 0;
        SourceStability[0] = WhichConfiguration==21?-11:2.0;//dummy
        DefSourceStability = 0;
        Potential[0] = "AV18";
        DefPot = 0;
        Baseline[0] = "Norm";
        Baseline[1] = "Baseline";
        Baseline[2] = "Baseline2";
        Baseline[3] = "Longbaseline";
        Baseline[4] = "Longbaseline2";
        DefBaseline = 0;

        FemtoRangeVars[0] = 0;//default
        FemtoRangeVars[1] = 1;//down
        FemtoRangeVars[2] = 2;//up
        DefFemtoRange = 0;

        FitRangeVars[0] = 0;//default
        FitRangeVars[1] = 1;//down
        FitRangeVars[2] = 2;//up
        DefFitRange = 0;

        pFracVars[0] = 0;//default
        pFracVars[1] = 1;
        pFracVars[2] = 2;
        DefProtonFrac = 0;

        LamFracVars[0] = 0;//default

        DefLambdaFrac = 0;

        kcVars[0] = 500;
        kcVars[1] = 750;
        kcVars[2] = 1000;
        DefKc = 0;

        MtVars[0] = -1;
        MtVars[1] = 0;
        MtVars[2] = 1;
        MtVars[3] = 2;
        MtVars[4] = 3;
        MtVars[5] = 4;
        MtVars[6] = 5;
        MtVars[7] = 6;
        DefMt = 0;
        for(int iData=0; iData<NumDataVars; iData++){
            DataVars[iData] = TString::Format("_%i",iData);
        }
        DefData = 0;
    }
    else if(WhichConfiguration==42){
        Source[0] = "Gauss";
        DefSource = 0;
        SourceScale[0] = -11;
        DefSourceScale = 0;
        SourceStability[0] = 2.0;//dummy
        DefSourceStability = 0;
        Potential[0] = "AV18";
        DefPot = 0;
        Baseline[0] = "Spline3_Tune1";
        //Baseline[0] = "Spline3_3";
        //Baseline[1] = "Spline3_4";
        //Baseline[2] = "Spline3_5";
        //Baseline[3] = "Spline3_7";
        //Baseline[4] = "Spline3_9";
        //Baseline[5] = "Spline3_11";
        DefBaseline = 0;

        FemtoRangeVars[0] = 0;//default
        FemtoRangeVars[1] = 1;//down
        FemtoRangeVars[2] = 2;//up
        DefFemtoRange = 0;

        FitRangeVars[0] = 0;//default
        FitRangeVars[1] = 1;//down
        FitRangeVars[2] = 2;//up
        DefFitRange = 0;

        pFracVars[0] = 0;//default
        pFracVars[1] = 1;
        pFracVars[2] = 2;
        DefProtonFrac = 0;

        LamFracVars[0] = 0;//default

        DefLambdaFrac = 0;

        kcVars[0] = 500;
        kcVars[1] = 750;
        kcVars[2] = 1000;

        DefKc = 0;

        MtVars[0] = -1;
        DefMt = 0;

        DataVars[0] = "_0";
        DefData = 0;

    }
    else if(WhichConfiguration==10){
        Source[0] = "McGauss_Reso";
        DefSource = 0;
        SourceScale[0] = -11;
        DefSourceScale = 0;
        SourceStability[0] = 2.0;//dummy
        DefSourceStability = 0;
        Potential[0] = "AV18";
        DefPot = 0;
        Baseline[0] = "Norm";
        Baseline[1] = "Baseline";
        Baseline[2] = "Baseline2";
        Baseline[3] = "Longbaseline";
        Baseline[4] = "Longbaseline2";
        DefBaseline = 0;

        FemtoRangeVars[0] = 0;//default
        FemtoRangeVars[1] = 1;//down
        FemtoRangeVars[2] = 2;//up
        DefFemtoRange = 0;

        FitRangeVars[0] = 0;//default
        FitRangeVars[1] = 1;//down
        FitRangeVars[2] = 2;//up
        DefFitRange = 0;

        pFracVars[0] = 0;//default
        pFracVars[1] = 1;
        pFracVars[2] = 2;
        DefProtonFrac = 0;

        LamFracVars[0] = 0;//default

        DefLambdaFrac = 0;

        kcVars[0] = 500;
        kcVars[1] = 750;
        kcVars[2] = 1000;
        DefKc = 0;

        MtVars[0] = -1;
        DefMt = 0;

        DataVars[0] = TString::Format("_%i",0);
        DefData = 0;
    }
    else if(WhichConfiguration==11){
        Source[0] = "McLevyNolan_Reso";
        DefSource = 0;
        SourceScale[0] = -11;
        DefSourceScale = 0;
        SourceStability[0] = -11;
        DefSourceStability = 0;
        Potential[0] = "AV18";
        DefPot = 0;
        Baseline[0] = "Norm";
        Baseline[1] = "Baseline";
        Baseline[2] = "Baseline2";
        Baseline[3] = "Longbaseline";
        Baseline[4] = "Longbaseline2";
        DefBaseline = 0;

        FemtoRangeVars[0] = 0;//default
        FemtoRangeVars[1] = 1;//down
        FemtoRangeVars[2] = 2;//up
        DefFemtoRange = 0;

        FitRangeVars[0] = 0;//default
        FitRangeVars[1] = 1;//down
        FitRangeVars[2] = 2;//up
        DefFitRange = 0;

        pFracVars[0] = 0;//default
        pFracVars[1] = 1;
        pFracVars[2] = 2;
        DefProtonFrac = 0;

        LamFracVars[0] = 0;//default

        DefLambdaFrac = 0;

        kcVars[0] = 500;
        kcVars[1] = 750;
        kcVars[2] = 1000;
        DefKc = 0;

        MtVars[0] = -1;
        DefMt = 0;
        DataVars[0] = TString::Format("_%i",0);
        DefData = 0;
    }
    //pure gauss
    else if(WhichConfiguration==12){
        Source[0] = "Gauss";
        DefSource = 0;
        SourceScale[0] = -11;
        DefSourceScale = 0;
        SourceStability[0] = 2.0;//dummy
        DefSourceStability = 0;
        Potential[0] = "AV18";
        DefPot = 0;
        Baseline[0] = "Norm";
        Baseline[1] = "Baseline";
        Baseline[2] = "Baseline2";
        Baseline[3] = "Longbaseline";
        Baseline[4] = "Longbaseline2";
        DefBaseline = 0;

        FemtoRangeVars[0] = 0;//default
        FemtoRangeVars[1] = 1;//down
        FemtoRangeVars[2] = 2;//up
        DefFemtoRange = 0;

        FitRangeVars[0] = 0;//default
        FitRangeVars[1] = 1;//down
        FitRangeVars[2] = 2;//up
        DefFitRange = 0;

        pFracVars[0] = 0;//default
        pFracVars[1] = 1;
        pFracVars[2] = 2;
        DefProtonFrac = 0;

        LamFracVars[0] = 0;//default

        DefLambdaFrac = 0;

        kcVars[0] = 500;
        kcVars[1] = 750;
        kcVars[2] = 1000;

        DefKc = 0;

        MtVars[0] = -1;
        DefMt = 0;
        DataVars[0] = TString::Format("_%i",0);
        DefData = 0;
    }
    else{
        printf("????????????\n");
    }

    //if(!JustNumIter){
    //    gROOT->cd();
    //    ntResult = new TNtuple("ntResult", "ntResult",
    //        "SourceType:SourceScale:SourceStability:Potential:Baseline:"
    //        "p_a:e_a:p_b:e_b:p_c:e_c:p_Cl:e_Cl:p_kc:e_kc:p_sor0:e_sor0:p_sor1:e_sor1:chi2:ndf");
    //}

    TRandom3 RanGen(RANDOMSEED);

    TString ntFileName = TString::Format("%sNTfile_%i.root",OutputFolder.Data(),WhichConfiguration);


    int IterationCounter = 0;
    int ComputedIterations = 0;
    for(int iSource=0; iSource<NumSourceVars; iSource++){
    for(int iSourceScal=0; iSourceScal<NumSourceScaleVars; iSourceScal++){
    for(int iSourceStab=0; iSourceStab<NumSourceStabilityVars; iSourceStab++){
    for(int iPot=0; iPot<NumPotVars; iPot++){
    for(int iBl=0; iBl<NumBaselineVars; iBl++){
    for(int iFemRan=0; iFemRan<NumFemtoRangeVars; iFemRan++){
    for(int iFitRan=0; iFitRan<NumFitRangeVars; iFitRan++){
    for(int ipFrac=0; ipFrac<NumProtonFracVars; ipFrac++){
    for(int iLamFrac=0; iLamFrac<NumLambdaFracVars; iLamFrac++){
    for(int ikc=0; ikc<NumKcVars; ikc++){
    for(int iMt=0; iMt<NumMtVars; iMt++){
    for(int iData=0; iData<NumDataVars; iData++){
        if((Baseline[iBl]=="Norm"||Baseline[iBl].Contains("Baseline"))&&ikc) {continue;}//this variation is only relevant in case of a Longbaseline or Splines
        if((Baseline[iBl]=="Norm"||Baseline[iBl].Contains("Baseline"))&&iFitRan) {continue;}//this variation is only relevant in case of a Longbaseline or Splines
        double UNI = RanGen.Uniform();
        //if(IterationCounter>=FirstIter&&IterationCounter<=LastIter){
        //    printf(" UNI=%.3f vs %.3f\n",UNI,1./double(OnlyFraction));
        //    printf(" IterationCounter=%i vs %i\n",IterationCounter,ComputedIterations);
        //}
        bool DefaultVar=true;
        if(iSourceScal!=DefSourceScale||iSourceStab!=DefSourceStability||iBl!=DefBaseline||
                    iFemRan!=DefFemtoRange||iFitRan!=DefFitRange||
                    ipFrac!=DefProtonFrac||iLamFrac!=DefLambdaFrac||ikc!=DefKc||iData!=DefData){DefaultVar=false;}
//printf("VAR %i = %i\n",IterationCounter,int(DefaultVar));
        if(DefaultVar) {DEFAULT_VAR[iSource][iPot] = IterationCounter; UNI=-1;}//we always compute the default variations (for all potentials) and for all baselines
        if(Mode==-1){
            if(DefaultVar){
                IterationCounter++;
                DEFAULT_VAR[iSource][iPot] = IterationCounter;
                continue;
            }
            else{
                continue;
            }
        }
        else if(Mode==1){
            if(DefaultVar){
                DEFAULT_VAR[iSource][iPot] = IterationCounter;
                IterationCounter++;
                continue;
            }
            else{
                IterationCounter++;
                continue;
            }
        }
        else{
            if(Mode!=-1&&(IterationCounter<FirstIter||IterationCounter>LastIter)) {IterationCounter++; continue;}
            unsigned EffectiveFraction = OnlyFraction;
            if((Baseline[iBl]=="Norm"||Baseline[iBl].Contains("Baseline"))){
                EffectiveFraction /= NumFemtoRangeVars;
                if(!Potential[iPot].Contains("Lednicky_")) EffectiveFraction /= NumFitRangeVars;
                if(EffectiveFraction<0||EffectiveFraction>OnlyFraction) EffectiveFraction=1;
            }
            //this is here in order to make sure we take more Norm events (now they are much less due to some conditions above)
            if(EffectiveFraction>1&&UNI>1./double(EffectiveFraction)) {IterationCounter++; continue;}
        }


        //TString IterationName = TString::Format("%spLam_%s_%s_Sca%.2f_Sta%.2f_%s_%s_",
        //                                        OutputFolder.Data(),"pp13TeV_HM_March19",Source[iSource].Data(),
        //                                        SourceScale[iSourceScal],SourceStability[iSourceStab],Potential[iPot],Baseline[iBl],
        //                                        );

        //[0] = femto range
        //[1] = fit range
        //[2] = fractions for p
        //[3] = fractions for L
        VARIATIONS[0] = FemtoRangeVars[iFemRan];
        VARIATIONS[1] = FitRangeVars[iFitRan];
        VARIATIONS[2] = pFracVars[ipFrac];
        VARIATIONS[3] = LamFracVars[iLamFrac];
        VARIATIONS[4] = kcVars[ikc];
        VARIATIONS[5] = MtVars[iMt];

        TString SourceDescription = Source[iSource];
        if(SourceScale[iSourceScal]<0) {SourceDescription+="ScaFree";}
        else if(SourceScale[iSourceScal]==0) {SourceDescription+="ScaDefault";}
        else {SourceDescription+=TString::Format("Sca%.2f",SourceScale[iSourceScal]);}
        if(SourceStability[iSourceStab]<0||SourceStability[iSourceStab]>2) {SourceDescription+="StaFree";}
        else if(SourceStability[iSourceStab]==0) {SourceDescription+="StaDefault";}
        else {SourceDescription+=TString::Format("Sta%.2f",SourceStability[iSourceStab]);}
        /*
        TString FileName = TString::Format("Iter%i_%s_%s_%s_%s_%i_%i_%i_%i_%i.root",IterationCounter,
                                                               "pp13TeV_HM_March19",SourceDescription.Data(),
                                                               Potential[iPot].Data(),Baseline[iBl].Data(),
                                                               VARIATIONS[0],VARIATIONS[1],VARIATIONS[2],VARIATIONS[3],VARIATIONS[4]);
        */
        TString FileName = TString::Format("Config%i_Iter%i.root",WhichConfiguration,IterationCounter);
        TFile* CheckIfFile = new TFile(OutputFolder+FileName,"READ");
        //! we ignore iterations for which we already have an input file!
        //! => delete the folder in case you want to overwrite
        if(CheckIfFile&&CheckIfFile->IsOpen()){IterationCounter++; continue;}
        delete CheckIfFile;
        CheckIfFile = NULL;

        Fit_pp(AnalysisObject, OutputFolder, FileName,"pp13TeV_HM_March19",DataVars[iData],Source[iSource],SourceScale[iSourceScal],SourceStability[iSourceStab],
               Potential[iPot],Baseline[iBl],VARIATIONS,long(IterationCounter)*long(1000)+long(WhichConfiguration),ntFileName,"ntResult");

        IterationCounter++;
        ComputedIterations++;
    }
    }
    }
    }
    }
    }
    }
    }
    }
    }
    }
    }

    printf("This set up had a total of %i iterations.\n",IterationCounter);

     for(int iSource=0; iSource<NumSourceVars; iSource++){
        for(int iPot=0; iPot<NumPotVars; iPot++){
            printf(" The default variation (%i,%i) is %i.\n",iSource,iPot,DEFAULT_VAR[iSource][iPot]);
        }
     }


    printf(" The function call evaluated %i iterations.\n",ComputedIterations);

    //if(NTfile){
    //    NTfile->cd();
    //    ntResult->Write();
    //}


    //if(ntResult) delete ntResult;
    //if(NTfile) delete NTfile;

    delete [] Source;
    delete [] SourceScale;
    delete [] SourceStability;
    delete [] Potential;
    delete [] Baseline;
    delete [] FemtoRangeVars;
    delete [] FitRangeVars;
    delete [] pFracVars;
    delete [] LamFracVars;
    delete [] kcVars;
    delete [] MtVars;
}


void Plot_mT_Scale(const TString OutputFolder, const TString ntFileName, const TString mtFileName, const TString pp_FileName,
                   const TString Descr1, const TString Descr2){

    TFile* mtFile = new TFile(mtFileName,"read");
    TGraphErrors* AveragemT = (TGraphErrors*)mtFile->Get("AveragemT");
    unsigned Num_mT_bins = AveragemT->GetN();
    TH2F* pL_LO_mT_sor0 = new TH2F("pL_LO_mT_sor0","pL_LO_mT_sor0",Num_mT_bins,-0.5,double(Num_mT_bins)-0.5,512,0.2,1.8);
    TH2F* pL_LO_mT_esor0 = new TH2F("pL_LO_mT_esor0","pL_LO_mT_esor0",Num_mT_bins,-0.5,double(Num_mT_bins)-0.5,512,0.0,0.4);

    TH2F* pL_NLO_mT_sor0 = new TH2F("pL_NLO_mT_sor0","pL_NLO_mT_sor0",Num_mT_bins,-0.5,double(Num_mT_bins)-0.5,512,0.2,1.8);
    TH2F* pL_NLO_mT_esor0 = new TH2F("pL_NLO_mT_esor0","pL_NLO_mT_esor0",Num_mT_bins,-0.5,double(Num_mT_bins)-0.5,512,0.0,0.4);

    TH2F* pL_mT_sor0 = new TH2F("pL_mT_sor0","pL_mT_sor0",Num_mT_bins,-0.5,double(Num_mT_bins)-0.5,512,0.2,1.8);
    TH2F* pL_mT_esor0 = new TH2F("pL_mT_esor0","pL_mT_esor0",Num_mT_bins,-0.5,double(Num_mT_bins)-0.5,512,0.0,0.4);


    TFile* ntFile = new TFile(ntFileName,"read");
    TNtuple* ntResult = (TNtuple*)ntFile->Get("ntResult");
    unsigned NumNtEntries = ntResult->GetEntries();

    Float_t IterM;
    Float_t Iter;
    Float_t Config;
    Float_t Data;
    Float_t SourceType;
    Float_t SourceScale;
    Float_t Potential;
    Float_t Baseline;
    Float_t FemRan;
    Float_t FitRan;
    Float_t pFrac;
    Float_t LamFrac;
    Float_t kcVar;
    Float_t mTbin;
    Float_t FemtoMin;
    Float_t FemtoMax;
    Float_t BlMin;
    Float_t BlMax;
    Float_t p_a;
    Float_t e_a;
    Float_t p_b;
    Float_t e_b;
    Float_t p_c;
    Float_t e_c;
    Float_t p_Cl;
    Float_t e_Cl;
    Float_t p_kc;
    Float_t e_kc;
    Float_t p_sor0;
    Float_t e_sor0;
    Float_t p_sor1;
    Float_t e_sor1;
    Float_t chi2;
    Float_t ndf;
    Float_t chi2_312;
    Float_t ndf_312;

    ntResult->SetBranchAddress("IterM",&IterM);
    ntResult->SetBranchAddress("Iter",&Iter);
    ntResult->SetBranchAddress("Config",&Config);
    ntResult->SetBranchAddress("Data",&Data);
    ntResult->SetBranchAddress("SourceType",&SourceType);
    ntResult->SetBranchAddress("SourceScale",&SourceScale);
    ntResult->SetBranchAddress("Potential",&Potential);
    ntResult->SetBranchAddress("Baseline",&Baseline);
    ntResult->SetBranchAddress("FemRan",&FemRan);
    ntResult->SetBranchAddress("FitRan",&FitRan);
    ntResult->SetBranchAddress("pFrac",&pFrac);
    ntResult->SetBranchAddress("LamFrac",&LamFrac);
    ntResult->SetBranchAddress("kcVar",&kcVar);
    ntResult->SetBranchAddress("mTbin",&mTbin);
    ntResult->SetBranchAddress("FemtoMin",&FemtoMin);
    ntResult->SetBranchAddress("FemtoMax",&FemtoMax);
    ntResult->SetBranchAddress("BlMin",&BlMin);
    ntResult->SetBranchAddress("BlMax",&BlMax);
    ntResult->SetBranchAddress("p_a",&p_a);
    ntResult->SetBranchAddress("e_a",&e_a);
    ntResult->SetBranchAddress("p_b",&p_b);
    ntResult->SetBranchAddress("e_b",&e_b);
    ntResult->SetBranchAddress("p_c",&p_c);
    ntResult->SetBranchAddress("e_c",&e_c);
    ntResult->SetBranchAddress("p_Cl",&p_Cl);
    ntResult->SetBranchAddress("e_Cl",&e_Cl);
    ntResult->SetBranchAddress("p_kc",&p_kc);
    ntResult->SetBranchAddress("e_kc",&e_kc);
    ntResult->SetBranchAddress("p_sor0",&p_sor0);
    ntResult->SetBranchAddress("e_sor0",&e_sor0);
    ntResult->SetBranchAddress("p_sor1",&p_sor1);
    ntResult->SetBranchAddress("e_sor1",&e_sor1);
    ntResult->SetBranchAddress("chi2",&chi2);
    ntResult->SetBranchAddress("ndf",&ndf);
    ntResult->SetBranchAddress("chi2_312",&chi2_312);
    ntResult->SetBranchAddress("ndf_312",&ndf_312);
printf("NumNtEntries=%u\n",NumNtEntries);
    for(unsigned uEntry=0; uEntry<NumNtEntries; uEntry++){
        ntFile->cd();
        ntResult->GetEntry(uEntry);
        //int Config;
        //Config = int(Iter)%100;
        //int ITER = int(Iter)/100;
        int Id_LO;
        int Id_NLO;
        if(Config<10){
            Id_LO = 1;
            Id_NLO = 11;
        }
        else if(Config<100){
            Id_LO = 0;
            Id_NLO = 10;
        }
//printf("Iter=%i; Config=%i; LO=%i; NLO=%i\n",ITER,Config,Id_LO,Id_NLO);
        //! Conditions

        //if(Config!=0&&Config!=2) continue; //select Gauss+Reso and Gauss
        //if(Config!=20) continue; //select Gauss+Reso
        //if(Config!=1) continue; //select Levy+Reso
        //if(Config!=2) continue; //select Gauss
        //if(round(Baseline)!=0) continue; //select only norm
        //if(int(Baseline)!=11) continue; //select only pol1
        //if(int(Baseline)!=12) continue; //select only pol2
        if(int(Baseline)==11) continue; // select pol0 and pol2
        //if(SourceScale>0) continue;
        //if(SourceScale<0.8&&SourceScale>0) continue;

        //if(mTbin<0) continue;

        if(Potential==1){
            pL_LO_mT_sor0->Fill(mTbin,p_sor0);
            pL_LO_mT_esor0->Fill(mTbin,e_sor0);
//printf("mTbin,p_sor0 = %f %f\n",mTbin,p_sor0);
        }
        if(Potential==11){
            pL_NLO_mT_sor0->Fill(mTbin,p_sor0);
            pL_NLO_mT_esor0->Fill(mTbin,e_sor0);
        }
        pL_mT_sor0->Fill(mTbin,p_sor0);
        pL_mT_esor0->Fill(mTbin,e_sor0);
    }

    TGraphErrors mTSyst_pL_LO;
    mTSyst_pL_LO.SetName("mTSyst_pL_LO");
    mTSyst_pL_LO.Set(Num_mT_bins);
    mTSyst_pL_LO.SetFillColorAlpha(kGreen+3,0.3);
    mTSyst_pL_LO.SetLineColor(kGreen+3);
    mTSyst_pL_LO.SetLineWidth(3);

    TGraphErrors mTSyst_pL_NLO;
    mTSyst_pL_NLO.SetName("mTSyst_pL_NLO");
    mTSyst_pL_NLO.Set(Num_mT_bins);
    mTSyst_pL_NLO.SetFillColorAlpha(kRed+1,0.3);
    mTSyst_pL_NLO.SetLineColor(kRed+1);
    mTSyst_pL_NLO.SetLineWidth(3);

    TGraphErrors mTSyst_pL;
    mTSyst_pL.SetName("mTSyst_pL");
    mTSyst_pL.Set(Num_mT_bins);
    mTSyst_pL.SetFillColorAlpha(kOrange+4,0.3);
    mTSyst_pL.SetLineColor(kOrange+4);
    mTSyst_pL.SetLineWidth(3);

    TGraphErrors mTStat_pL_LO;
    mTStat_pL_LO.SetName("mTStat_pL_LO");
    mTStat_pL_LO.Set(Num_mT_bins);
    mTStat_pL_LO.SetFillColorAlpha(kGreen+3,0.3);
    mTStat_pL_LO.SetLineColor(kGreen+3);
    mTStat_pL_LO.SetLineWidth(3);

    TGraphErrors mTStat_pL_NLO;
    mTStat_pL_NLO.SetName("mTStat_pL_NLO");
    mTStat_pL_NLO.Set(Num_mT_bins);
    mTStat_pL_NLO.SetFillColorAlpha(kRed+1,0.3);
    mTStat_pL_NLO.SetLineColor(kRed+1);
    mTStat_pL_NLO.SetLineWidth(3);

    TGraphErrors mTStat_pL;
    mTStat_pL.SetName("mTStat_pL");
    mTStat_pL.Set(Num_mT_bins);
    mTStat_pL.SetFillColorAlpha(kOrange+4,0.3);
    mTStat_pL.SetLineColor(kOrange+4);
    mTStat_pL.SetLineWidth(3);

    TGraphErrors* mTSyst_pp=NULL;
    TGraphErrors* mTStat_pp=NULL;

    TFile* ppFile = NULL;
    if(pp_FileName!=""){
        ppFile = new TFile(pp_FileName,"read");
        mTSyst_pp = (TGraphErrors*)ppFile->Get("mTRadiusSyst");
        mTStat_pp = (TGraphErrors*)ppFile->Get("mTRadiusStat");
        //mTSyst_pp = (TGraphErrors*)ppFile->Get("Radius_pp");
        //mTStat_pp = (TGraphErrors*)ppFile->Get("Radius_pp");
    }

    if(mTSyst_pp){
        mTSyst_pp->SetName("mTSyst_pp");
        mTSyst_pp->SetFillColorAlpha(kBlue+2,0.3);
        mTSyst_pp->SetLineColor(kBlue+2);
        mTSyst_pp->SetLineWidth(3);
    }

    if(mTStat_pp){
        mTStat_pp->SetName("mTStat_pp");
        mTStat_pp->SetFillColorAlpha(kBlue+2,0.3);
        mTStat_pp->SetLineColor(kBlue+2);
        mTStat_pp->SetLineWidth(3);
        for(unsigned uPoint=0; uPoint<mTStat_pp->GetN(); uPoint++){
            mTStat_pp->SetPointError(uPoint,0,mTStat_pp->GetErrorY(uPoint));
        }
    }


    for(unsigned umT=0; umT<Num_mT_bins; umT++){
        double AvgmT;
        double AvgmTerr;
        double Dummy;

        TH1D* sorLO = pL_LO_mT_sor0->ProjectionY("sorLO",umT+1,umT+1);
        TH1D* sorNLO = pL_NLO_mT_sor0->ProjectionY("sorNLO",umT+1,umT+1);
        TH1D* sor = pL_mT_sor0->ProjectionY("sor",umT+1,umT+1);
        TH1D* esorLO = pL_LO_mT_esor0->ProjectionY("esorLO",umT+1,umT+1);
        TH1D* esorNLO = pL_NLO_mT_esor0->ProjectionY("esorNLO",umT+1,umT+1);
        TH1D* esor = pL_mT_esor0->ProjectionY("esor",umT+1,umT+1);
        double StatErrLO = esorLO->GetMean();
        double SystErrLO = sorLO->GetStdDev();
        double StatErrNLO = esorNLO->GetMean();
        double SystErrNLO = sorNLO->GetStdDev();
        double StatErr = esor->GetMean();
        double SystErr = sor->GetStdDev();

        AveragemT->GetPoint(umT,Dummy,AvgmT);
        AvgmTerr = AveragemT->GetErrorY(umT);

        mTSyst_pL_LO.SetPoint(umT,AvgmT,sorLO->GetMean());//
        mTSyst_pL_LO.SetPointError(umT,AvgmTerr,SystErrLO);

        mTSyst_pL_NLO.SetPoint(umT,AvgmT,sorNLO->GetMean());//
        mTSyst_pL_NLO.SetPointError(umT,AvgmTerr,SystErrNLO);

        mTSyst_pL.SetPoint(umT,AvgmT,sor->GetMean());//
        mTSyst_pL.SetPointError(umT,AvgmTerr,SystErr);

//printf("sorLO->GetMean()=%f; sorNLO->GetMean()=%f\n",sorLO->GetMean(),sorNLO->GetMean());
        mTStat_pL_LO.SetPoint(umT,AvgmT,sorLO->GetMean());//
        mTStat_pL_LO.SetPointError(umT,0,StatErrLO);
        mTStat_pL_NLO.SetPoint(umT,AvgmT,sorNLO->GetMean());//
        mTStat_pL_NLO.SetPointError(umT,0,StatErrNLO);
        mTStat_pL.SetPoint(umT,AvgmT,sor->GetMean());//
        mTStat_pL.SetPointError(umT,0,StatErr);
    }

    DLM_Histo<double> DimiStat_pp;
    DLM_Histo<double> DimiStat_pL_LO;
    DLM_Histo<double> DimiStat_pL_NLO;
    DLM_Histo<double> DimiStat_pL;

    DLM_Histo<double> DimiStatEff_pL_NLO;
    DLM_Histo<double> DimiStatRatio_pL_NLO;

    DLM_Histo<double> DimiStatEff_pL;
    DLM_Histo<double> DimiStatRatio_pL;

    double xValLow;
    double yValLow;
    double xValUp;
    double yValUp;

    double xVal;
    double yVal;

    DimiStat_pp.SetUp(1);
    mTStat_pp->GetPoint(0,xValLow,yValLow);
    mTStat_pp->GetPoint(0,xValUp,yValUp);
    DimiStat_pp.SetUp(0,mTStat_pp->GetN(),xValLow-0.001,xValUp+0.001);
    DimiStat_pp.Initialize();
printf("mT Radius Stat Syst\n");
printf(" pp AV18:\n");
    for(unsigned uBin=0; uBin<mTStat_pp->GetN(); uBin++){
        mTStat_pp->GetPoint(uBin,xVal,yVal);
        DimiStat_pp.SetBinCenter(0,uBin,xVal);
        DimiStat_pp.SetBinContent(uBin,yVal);
        DimiStat_pp.SetBinError(uBin,mTStat_pp->GetErrorY(uBin));
printf("  %.3f   %.3f   %.3f   %.3f\n",xVal,yVal,mTStat_pp->GetErrorY(uBin),mTSyst_pp->GetErrorY(uBin));
//printf(" DimiStat_pp: cen/con/err = %.3f/%.3f/%.3f\n",DimiStat_pp.GetBinCenter(0,uBin),DimiStat_pp.GetBinContent(uBin),DimiStat_pp.GetBinError(uBin));
    }

    DimiStat_pL_LO.SetUp(1);
    mTStat_pL_LO.GetPoint(0,xValLow,yValLow);
    mTStat_pL_LO.GetPoint(0,xValUp,yValUp);
    DimiStat_pL_LO.SetUp(0,mTStat_pL_LO.GetN(),xValLow-0.001,xValUp+0.001);
    DimiStat_pL_LO.Initialize();
printf("mT Radius Stat Syst\n");
printf(" pL LO:\n");
    for(unsigned uBin=0; uBin<mTStat_pL_LO.GetN(); uBin++){
        mTStat_pL_LO.GetPoint(uBin,xVal,yVal);
        DimiStat_pL_LO.SetBinCenter(0,uBin,xVal);
        DimiStat_pL_LO.SetBinContent(uBin,yVal);
        DimiStat_pL_LO.SetBinError(uBin,mTStat_pL_LO.GetErrorY(uBin));
printf("  %.3f   %.3f   %.3f   %.3f\n",xVal,yVal,mTStat_pL_LO.GetErrorY(uBin),mTSyst_pL_LO.GetErrorY(uBin));
    }

    DimiStat_pL_NLO.SetUp(1);
    mTStat_pL_NLO.GetPoint(0,xValLow,yValLow);
    mTStat_pL_NLO.GetPoint(0,xValUp,yValUp);
    DimiStat_pL_NLO.SetUp(0,mTStat_pL_NLO.GetN(),xValLow-0.001,xValUp+0.001);
    DimiStat_pL_NLO.Initialize();
printf("mT Radius Stat Syst\n");
printf(" pL NLO:\n");
    for(unsigned uBin=0; uBin<mTStat_pL_NLO.GetN(); uBin++){
        mTStat_pL_NLO.GetPoint(uBin,xVal,yVal);
        DimiStat_pL_NLO.SetBinCenter(0,uBin,xVal);
        DimiStat_pL_NLO.SetBinContent(uBin,yVal);
        DimiStat_pL_NLO.SetBinError(uBin,mTStat_pL_NLO.GetErrorY(uBin));
printf("  %.3f   %.3f   %.3f   %.3f\n",xVal,yVal,mTStat_pL_NLO.GetErrorY(uBin),mTSyst_pL_NLO.GetErrorY(uBin));
//printf(" cen/con/err = %.3f/%.3f/%.3f\n",DimiStat_pL_NLO.GetBinCenter(0,uBin),DimiStat_pL_NLO.GetBinContent(uBin),DimiStat_pL_NLO.GetBinError(uBin));
    }


    DimiStatEff_pL_NLO.Copy(DimiStat_pp);
printf("DimiStatEff_pL_NLO.GetNbins()=%u\n",DimiStatEff_pL_NLO.GetNbins());
    for(unsigned uBin=0; uBin<DimiStatEff_pL_NLO.GetNbins(); uBin++){
        double mT_VALUE = DimiStatEff_pL_NLO.GetBinCenter(0,uBin);
        DimiStatEff_pL_NLO.SetBinContent(uBin,DimiStat_pL_NLO.Eval(&mT_VALUE));
        DimiStatEff_pL_NLO.SetBinError(uBin,DimiStat_pL_NLO.EvalError(&mT_VALUE));
//printf(" NLO: mT_VALUE=%.3f; Value=%.3f+/-%.3f\n",mT_VALUE,DimiStat_pL_NLO.Eval(&mT_VALUE),DimiStat_pL_NLO.EvalError(&mT_VALUE));
    }

    DimiStatRatio_pL_NLO.Copy(DimiStatEff_pL_NLO);
    DimiStatRatio_pL_NLO.DivideHisto(DimiStat_pp);

    TGraphErrors gStatRatio_pL_NLO;
    gStatRatio_pL_NLO.SetName("gStatRatio_pL_NLO");
    gStatRatio_pL_NLO.Set(DimiStatRatio_pL_NLO.GetNbins());
    for(unsigned uBin=0; uBin<DimiStatRatio_pL_NLO.GetNbins(); uBin++){
        gStatRatio_pL_NLO.SetPoint(uBin,DimiStatRatio_pL_NLO.GetBinCenter(0,uBin),DimiStatRatio_pL_NLO.GetBinContent(uBin));
        gStatRatio_pL_NLO.SetPointError(uBin,mTSyst_pp->GetErrorX(uBin),DimiStatRatio_pL_NLO.GetBinError(uBin));
    }


    DimiStat_pL.SetUp(1);
    mTStat_pL.GetPoint(0,xValLow,yValLow);
    mTStat_pL.GetPoint(0,xValUp,yValUp);
    DimiStat_pL.SetUp(0,mTStat_pL.GetN(),xValLow-0.001,xValUp+0.001);
    DimiStat_pL.Initialize();
printf("mT Radius Stat Syst\n");
printf(" pL NLO:\n");
    for(unsigned uBin=0; uBin<mTStat_pL.GetN(); uBin++){
        mTStat_pL.GetPoint(uBin,xVal,yVal);
        DimiStat_pL.SetBinCenter(0,uBin,xVal);
        DimiStat_pL.SetBinContent(uBin,yVal);
        DimiStat_pL.SetBinError(uBin,mTStat_pL.GetErrorY(uBin));
printf("  %.3f   %.3f   %.3f   %.3f\n",xVal,yVal,mTStat_pL.GetErrorY(uBin),mTSyst_pL.GetErrorY(uBin));
//printf(" cen/con/err = %.3f/%.3f/%.3f\n",DimiStat_pL.GetBinCenter(0,uBin),DimiStat_pL.GetBinContent(uBin),DimiStat_pL.GetBinError(uBin));
    }


    DimiStatEff_pL.Copy(DimiStat_pp);
printf("DimiStatEff_pL.GetNbins()=%u\n",DimiStatEff_pL.GetNbins());
    for(unsigned uBin=0; uBin<DimiStatEff_pL.GetNbins(); uBin++){
        double mT_VALUE = DimiStatEff_pL.GetBinCenter(0,uBin);
        DimiStatEff_pL.SetBinContent(uBin,DimiStat_pL.Eval(&mT_VALUE));
        DimiStatEff_pL.SetBinError(uBin,DimiStat_pL.EvalError(&mT_VALUE));
//printf(" NLO: mT_VALUE=%.3f; Value=%.3f+/-%.3f\n",mT_VALUE,DimiStat_pL.Eval(&mT_VALUE),DimiStat_pL.EvalError(&mT_VALUE));
    }


    DimiStatRatio_pL.Copy(DimiStatEff_pL);
    DimiStatRatio_pL.DivideHisto(DimiStat_pp);

    TGraphErrors gStatRatio_pL;
    gStatRatio_pL.SetName("gStatRatio_pL");
    gStatRatio_pL.Set(DimiStatRatio_pL.GetNbins());
    for(unsigned uBin=0; uBin<DimiStatRatio_pL.GetNbins(); uBin++){
        gStatRatio_pL.SetPoint(uBin,DimiStatRatio_pL.GetBinCenter(0,uBin),DimiStatRatio_pL.GetBinContent(uBin));
        gStatRatio_pL.SetPointError(uBin,mTSyst_pp->GetErrorX(uBin),DimiStatRatio_pL.GetBinError(uBin));
    }

    TFile* fOut = new TFile(OutputFolder+"fOut.root","recreate");

    TH1F* hAxis = new TH1F("hAxis", "hAxis", 128, 0.9, 2.7);
    hAxis->SetStats(false);
    hAxis->SetTitle("");
    hAxis->GetXaxis()->SetTitle("<#it{m}_{T}> (GeV/#it{c}^{2})");
    hAxis->GetXaxis()->SetTitleSize(0.06);
    hAxis->GetXaxis()->SetLabelSize(0.06);
    //hAxis->GetXaxis()->CenterTitle();
    hAxis->GetXaxis()->SetTitleOffset(1.3);
    hAxis->GetXaxis()->SetLabelOffset(0.02);

    //if(LevySource) hAxis->GetYaxis()->SetTitle("Levy core (fm)");
    //else hAxis->GetYaxis()->SetTitle("Gaussian core (fm)");
    TString FileDescr;
    if(Descr1.Contains("core")){
        hAxis->GetYaxis()->SetTitle("r_{core} (fm)");
        hAxis->GetYaxis()->SetRangeUser(0.55, 1.35);
        FileDescr = "GaussReso";
    }
    else{
        hAxis->GetYaxis()->SetTitle("r_{0} (fm)");
        hAxis->GetYaxis()->SetRangeUser(0.90, 1.80);
        FileDescr = "Gauss";
    }
    //
    hAxis->GetYaxis()->SetTitleSize(0.06);
    hAxis->GetYaxis()->SetLabelSize(0.06);
    //hAxis->GetYaxis()->CenterTitle();
    hAxis->GetYaxis()->SetTitleOffset(0.90);

    //hAxis->GetXaxis()->SetNdivisions(506);


    //

    //TLegend* lLegend = new TLegend(0.45,0.85,0.95,0.95);//lbrt
    //lLegend->SetName(TString::Format("lLegend"));
    //lLegend->SetTextSize(0.04);
    //lLegend->SetNColumns(3);
    TLegend* lLegend = new TLegend(0.67,0.48,0.90,0.71);//lbrt
    lLegend->SetName(TString::Format("lLegend"));
    lLegend->SetTextSize(0.04);
    lLegend->SetBorderSize(0);
    if(mTStat_pp)lLegend->AddEntry(mTStat_pp,"p#minusp (AV18)");
    lLegend->AddEntry(&mTSyst_pL_LO,"p#minus#Lambda (LO)");
    lLegend->AddEntry(&mTSyst_pL_NLO,"p#minus#Lambda (NLO)");
    //if(mT_pXim)lLegend->AddEntry(mT_pXim,"p#minus#Xi^{#minus} (HAL QCD)");

    TCanvas* cmT = new TCanvas("cmT", "cmT", 1);
    cmT->cd(0); cmT->SetCanvasSize(1920/2, 1080/2); cmT->SetMargin(0.15,0.05,0.2,0.05);//lrbt

    hAxis->Draw("axis");
    mTSyst_pL_NLO.Draw("e2,same");
    mTSyst_pL_LO.Draw("e2,same");
    if(mTSyst_pp) mTSyst_pp->Draw("e2,same");
    mTStat_pL_NLO.Draw("e,same");
    mTStat_pL_LO.Draw("e,same");
    if(mTStat_pp) mTStat_pp->Draw("e,same");
    lLegend->Draw("same");

    //TLatex BeamText;
    //BeamText.SetTextSize(gStyle->GetTextSize()*0.90);
    //BeamText.SetTextFont(42);
    //BeamText.SetNDC(kTRUE);
    //BeamText.DrawLatex(0.68, 0.76, "ALICE Preliminary");
    //BeamText.DrawLatex(0.68, 0.71, "high-mult. (0-0.072% INEL) pp #sqrt{#it{s}} = 13 TeV");

    TLatex BeamText;
    BeamText.SetTextSize(0.04);
    BeamText.SetTextFont(42);
    BeamText.SetNDC(kTRUE);
    BeamText.DrawLatex(0.677, 0.90, "ALICE Preliminary");
    TLatex BeamText1;
    BeamText1.SetTextSize(0.04);
    BeamText1.SetTextFont(42);
    BeamText1.SetNDC(kTRUE);
    BeamText1.DrawLatex(0.677, 0.85, "high-mult. (0-0.072% INEL)");
    BeamText1.DrawLatex(0.677, 0.81, "pp #sqrt{#it{s}} = 13 TeV");

    TLatex BeamText2;
    BeamText2.SetTextSize(0.035);
    BeamText2.SetTextFont(42);
    BeamText2.SetNDC(kTRUE);
    if(Descr2!=""){
        BeamText2.DrawLatex(0.677, 0.76, Descr1);
        BeamText2.DrawLatex(0.677, 0.72, Descr2);
    }
    else{
        BeamText2.DrawLatex(0.677, 0.74, Descr1);
    }


    TLegend* lLegend_TOT = new TLegend(0.67,0.48,0.90,0.71);//lbrt
    lLegend_TOT->SetName(TString::Format("lLegend_TOT"));
    lLegend_TOT->SetTextSize(0.04);
    lLegend_TOT->SetBorderSize(0);
    if(mTStat_pp)lLegend_TOT->AddEntry(mTStat_pp,"p#minusp (AV18)");
    lLegend_TOT->AddEntry(&mTSyst_pL,"p#minus#Lambda (LO,NLO,Usmani)");
    //if(mT_pXim)lLegend->AddEntry(mT_pXim,"p#minus#Xi^{#minus} (HAL QCD)");

    TCanvas* cmT_TOT = new TCanvas("cmT_TOT", "cmT_TOT", 1);
    cmT_TOT->cd(0); cmT_TOT->SetCanvasSize(1920/2, 1080/2); cmT_TOT->SetMargin(0.15,0.05,0.2,0.05);//lrbt

    hAxis->Draw("axis");
    mTSyst_pL.Draw("e2,same");
    if(mTSyst_pp) mTSyst_pp->Draw("e2,same");
    mTStat_pL.Draw("e,same");
    if(mTStat_pp) mTStat_pp->Draw("e,same");
    lLegend_TOT->Draw("same");

    //TLatex BeamText;
    //BeamText.SetTextSize(gStyle->GetTextSize()*0.90);
    //BeamText.SetTextFont(42);
    //BeamText.SetNDC(kTRUE);
    //BeamText.DrawLatex(0.68, 0.76, "ALICE Preliminary");
    //BeamText.DrawLatex(0.68, 0.71, "high-mult. (0-0.072% INEL) pp #sqrt{#it{s}} = 13 TeV");

    //TLatex BeamText;
    //BeamText.SetTextSize(0.04);
    //BeamText.SetTextFont(42);
    //BeamText.SetNDC(kTRUE);
    BeamText.DrawLatex(0.677, 0.90, "ALICE Preliminary");
    //TLatex BeamText1;
    //BeamText1.SetTextSize(0.04);
    //BeamText1.SetTextFont(42);
    //BeamText1.SetNDC(kTRUE);
    BeamText1.DrawLatex(0.677, 0.85, "high-mult. (0-0.072% INEL)");
    BeamText1.DrawLatex(0.677, 0.81, "pp #sqrt{#it{s}} = 13 TeV");

    //TLatex BeamText2;
    //BeamText2.SetTextSize(0.035);
    //BeamText2.SetTextFont(42);
    //BeamText2.SetNDC(kTRUE);
    if(Descr2!=""){
        BeamText2.DrawLatex(0.677, 0.76, Descr1);
        BeamText2.DrawLatex(0.677, 0.72, Descr2);
    }
    else{
        BeamText2.DrawLatex(0.677, 0.74, Descr1);
    }




    fOut->cd();
    mTSyst_pL_LO.Write();
    mTSyst_pL_NLO.Write();
    mTSyst_pL.Write();
    mTStat_pL_LO.Write();
    mTStat_pL_NLO.Write();
    mTStat_pL.Write();
    if(mTSyst_pp) mTSyst_pp->Write();
    if(mTStat_pp) mTStat_pp->Write();
    gStatRatio_pL_NLO.Write();
    gStatRatio_pL.Write();
    cmT_TOT->Write();
    cmT_TOT->SaveAs(OutputFolder+"cmT_"+FileDescr+"_TOT.pdf");


    const double AvgMt_pL = 1.55;
    unsigned ClosestBin = 3;
    double ExpectedStat;
    double ExpectedSyst;
    double ExpectedTot;
    double LowestStat=100;
    double HighestStat=0;
    double LowestSyst=100;
    double HighestSyst=0;


/*
    mTSyst_pL_LO.GetPoint(ClosestBin,xVal,yVal);
    yErr = mTSyst_pL_LO.GetErrorY(ClosestBin);
    if(LowestStat>(yVal-yErr)){LowestStat=yVal-yErr;}

    mTSyst_pL_NLO.GetPoint(ClosestBin,xVal,yVal);
    yErr = mTSyst_pL_NLO.GetErrorY(ClosestBin);
    if(LowestStat>(yVal-yErr)){LowestStat=yVal-yErr;}

    mTSyst_pp->GetPoint(ClosestBin,xVal,yVal);
    yErr = mTSyst_pp->GetErrorY(ClosestBin);
    if(LowestStat>(yVal-yErr)){LowestStat=yVal-yErr;}

    mTSyst_pL_LO.GetPoint(ClosestBin,xVal,yVal);
    yErr = mTSyst_pL_LO.GetErrorY(ClosestBin);
    if(HighestStat<(yVal+yErr)){HighestStat=yVal+yErr;}

    mTSyst_pL_NLO.GetPoint(ClosestBin,xVal,yVal);
    yErr = mTSyst_pL_NLO.GetErrorY(ClosestBin);
    if(HighestStat<(yVal+yErr)){HighestStat=yVal+yErr;}

    mTSyst_pp->GetPoint(ClosestBin,xVal,yVal);
    yErr = mTSyst_pp->GetErrorY(ClosestBin);
    if(HighestStat<(yVal+yErr)){HighestStat=yVal+yErr;}

    mTSyst_pL_LO.GetPoint(ClosestBin,xVal,yVal);
    yErr = mTSyst_pL_LO.GetErrorY(ClosestBin);
    if(LowestSyst>(yVal-yErr)){LowestSyst=yVal-yErr;}

    mTSyst_pL_NLO.GetPoint(ClosestBin,xVal,yVal);
    yErr = mTSyst_pL_NLO.GetErrorY(ClosestBin);
    if(LowestSyst>(yVal-yErr)){LowestSyst=yVal-yErr;}

    mTSyst_pp->GetPoint(ClosestBin,xVal,yVal);
    yErr = mTSyst_pp->GetErrorY(ClosestBin);
    if(LowestSyst>(yVal-yErr)){LowestSyst=yVal-yErr;}

    mTSyst_pL_LO.GetPoint(ClosestBin,xVal,yVal);
    yErr = mTSyst_pL_LO.GetErrorY(ClosestBin);
    if(HighestSyst<(yVal+yErr)){HighestSyst=yVal+yErr;}

    mTSyst_pL_NLO.GetPoint(ClosestBin,xVal,yVal);
    yErr = mTSyst_pL_NLO.GetErrorY(ClosestBin);
    if(HighestSyst<(yVal+yErr)){HighestSyst=yVal+yErr;}

    mTSyst_pp->GetPoint(ClosestBin,xVal,yVal);
    yErr = mTSyst_pp->GetErrorY(ClosestBin);
    if(HighestSyst<(yVal+yErr)){HighestSyst=yVal+yErr;}
*/



    delete pL_LO_mT_sor0;
    delete pL_LO_mT_esor0;
    delete pL_NLO_mT_sor0;
    delete pL_NLO_mT_esor0;
    delete pL_mT_sor0;
    delete pL_mT_esor0;
    delete hAxis;
    delete lLegend;
    delete lLegend_TOT;
    delete cmT;
    delete cmT_TOT;
    delete ntFile;
    delete mtFile;
    delete fOut;
    if(ppFile){delete ppFile; ppFile=NULL;}

}

void Plot_mT_Scale_pp_pL(const TString OutputFolder, const TString ntFileName_pp, const TString mtFileName_pp, const TString ntFileName_pL, const TString mtFileName_pL,
                   const TString Descr1, const TString Descr2){

    TFile* mtFile_pp = new TFile(mtFileName_pp,"read");
    TGraphErrors* AveragemT_pp = (TGraphErrors*)mtFile_pp->Get("AveragemT");
    unsigned Num_mT_bins_pp = AveragemT_pp->GetN();
    //TH2F* pp_AV18_mT_sor0 = new TH2F("pp_AV18_mT_sor0","pp_AV18_mT_sor0",Num_mT_bins_pp,-0.5,double(Num_mT_bins_pp)-0.5,512,0.2,1.8);
    //TH2F* pp_AV18_mT_esor0 = new TH2F("pp_AV18_mT_esor0","pp_AV18_mT_esor0",Num_mT_bins_pp,-0.5,double(Num_mT_bins_pp)-0.5,512,0.0,0.4);

    TH2F* pp_mT_sor0 = new TH2F("pp_mT_sor0","pp_mT_sor0",Num_mT_bins_pp,-0.5,double(Num_mT_bins_pp)-0.5,512,0.2,1.8);
    TH2F* pp_mT_esor0 = new TH2F("pp_mT_esor0","pp_mT_esor0",Num_mT_bins_pp,-0.5,double(Num_mT_bins_pp)-0.5,512,0.0,0.4);
    TH2F* pp_mT_sor1 = new TH2F("pp_mT_sor1","pp_mT_sor1",Num_mT_bins_pp,-0.5,double(Num_mT_bins_pp)-0.5,512,1.0,2.0);
    TH2F* pp_mT_esor1 = new TH2F("pp_mT_esor1","pp_mT_esor1",Num_mT_bins_pp,-0.5,double(Num_mT_bins_pp)-0.5,512,0.0,1.0);

    TFile* ntFile_pp = new TFile(ntFileName_pp,"read");
    TNtuple* ntResult_pp = (TNtuple*)ntFile_pp->Get("ntResult");
    unsigned NumNtEntries_pp = ntResult_pp->GetEntries();

    Float_t IterM_pp;
    Float_t Iter_pp;
    Float_t Config_pp;
    Float_t Data_pp;
    Float_t SourceType_pp;
    Float_t SourceScale_pp;
    Float_t Potential_pp;
    Float_t Baseline_pp;
    Float_t FemRan_pp;
    Float_t FitRan_pp;
    Float_t pFrac_pp;
    Float_t LamFrac_pp;
    Float_t kcVar_pp;
    Float_t mTbin_pp;
    Float_t FemtoMin_pp;
    Float_t FemtoMax_pp;
    Float_t BlMin_pp;
    Float_t BlMax_pp;
    Float_t p_a_pp;
    Float_t e_a_pp;
    Float_t p_b_pp;
    Float_t e_b_pp;
    Float_t p_c_pp;
    Float_t e_c_pp;
    Float_t p_Cl_pp;
    Float_t e_Cl_pp;
    Float_t p_kc_pp;
    Float_t e_kc_pp;
    Float_t p_sor0_pp;
    Float_t e_sor0_pp;
    Float_t p_sor1_pp;
    Float_t e_sor1_pp;
    Float_t chi2_pp;
    Float_t ndf_pp;
    Float_t chi2_352_pp;
    Float_t ndf_352_pp;

    ntResult_pp->SetBranchAddress("IterM",&IterM_pp);
    ntResult_pp->SetBranchAddress("Iter",&Iter_pp);
    ntResult_pp->SetBranchAddress("Config",&Config_pp);
    ntResult_pp->SetBranchAddress("Data",&Data_pp);
    ntResult_pp->SetBranchAddress("SourceType",&SourceType_pp);
    ntResult_pp->SetBranchAddress("SourceScale",&SourceScale_pp);
    ntResult_pp->SetBranchAddress("Potential",&Potential_pp);
    ntResult_pp->SetBranchAddress("Baseline",&Baseline_pp);
    ntResult_pp->SetBranchAddress("FemRan",&FemRan_pp);
    ntResult_pp->SetBranchAddress("FitRan",&FitRan_pp);
    ntResult_pp->SetBranchAddress("pFrac",&pFrac_pp);
    ntResult_pp->SetBranchAddress("LamFrac",&LamFrac_pp);
    ntResult_pp->SetBranchAddress("kcVar",&kcVar_pp);
    ntResult_pp->SetBranchAddress("mTbin",&mTbin_pp);
    ntResult_pp->SetBranchAddress("FemtoMin",&FemtoMin_pp);
    ntResult_pp->SetBranchAddress("FemtoMax",&FemtoMax_pp);
    ntResult_pp->SetBranchAddress("BlMin",&BlMin_pp);
    ntResult_pp->SetBranchAddress("BlMax",&BlMax_pp);
    ntResult_pp->SetBranchAddress("p_a",&p_a_pp);
    ntResult_pp->SetBranchAddress("e_a",&e_a_pp);
    ntResult_pp->SetBranchAddress("p_b",&p_b_pp);
    ntResult_pp->SetBranchAddress("e_b",&e_b_pp);
    ntResult_pp->SetBranchAddress("p_c",&p_c_pp);
    ntResult_pp->SetBranchAddress("e_c",&e_c_pp);
    ntResult_pp->SetBranchAddress("p_Cl",&p_Cl_pp);
    ntResult_pp->SetBranchAddress("e_Cl",&e_Cl_pp);
    ntResult_pp->SetBranchAddress("p_kc",&p_kc_pp);
    ntResult_pp->SetBranchAddress("e_kc",&e_kc_pp);
    ntResult_pp->SetBranchAddress("p_sor0",&p_sor0_pp);
    ntResult_pp->SetBranchAddress("e_sor0",&e_sor0_pp);
    ntResult_pp->SetBranchAddress("p_sor1",&p_sor1_pp);
    ntResult_pp->SetBranchAddress("e_sor1",&e_sor1_pp);
    ntResult_pp->SetBranchAddress("chi2",&chi2_pp);
    ntResult_pp->SetBranchAddress("ndf",&ndf_pp);
    ntResult_pp->SetBranchAddress("chi2_352",&chi2_352_pp);
    ntResult_pp->SetBranchAddress("ndf_352",&ndf_352_pp);
printf("NumNtEntries_pp=%u\n",NumNtEntries_pp);
    for(unsigned uEntry=0; uEntry<NumNtEntries_pp; uEntry++){
        ntFile_pp->cd();
        ntResult_pp->GetEntry(uEntry);

        //! Conditions

        //if(Config!=0&&Config!=2) continue; //select Gauss+Reso and Gauss
        //if(Config!=20) continue; //select Gauss+Reso
        //if(Config!=1) continue; //select Levy+Reso
        //if(Config!=2) continue; //select Gauss
        //if(round(Baseline)!=0) continue; //select only norm
        //if(int(Baseline)!=11) continue; //select only pol1
        //if(int(Baseline)!=12) continue; //select only pol2
        //if(int(Baseline_pp)==11) continue; // select pol0 and pol2
        //if(SourceScale>0) continue;
        //if(SourceScale<0.8&&SourceScale>0) continue;

        //if(mTbin<0) continue;

        //pp_AV18_mT_sor0->Fill(mTbin_pp,p_sor0_pp);
        //pp_AV18_mT_esor0->Fill(mTbin_pp,e_sor0_pp);

        pp_mT_sor0->Fill(mTbin_pp,p_sor0_pp);
        pp_mT_esor0->Fill(mTbin_pp,e_sor0_pp);
        pp_mT_sor1->Fill(mTbin_pp,p_sor1_pp);
        pp_mT_esor1->Fill(mTbin_pp,e_sor1_pp);
    }

    /////////////////////////////////

    TFile* mtFile_pL = new TFile(mtFileName_pL,"read");
    TGraphErrors* AveragemT_pL = (TGraphErrors*)mtFile_pL->Get("AveragemT");
    unsigned Num_mT_bins_pL = AveragemT_pL->GetN();
    TH2F* pL_LO_mT_sor0 = new TH2F("pL_LO_mT_sor0","pL_LO_mT_sor0",Num_mT_bins_pL,-0.5,double(Num_mT_bins_pL)-0.5,512,0.2,1.8);
    TH2F* pL_LO_mT_esor0 = new TH2F("pL_LO_mT_esor0","pL_LO_mT_esor0",Num_mT_bins_pL,-0.5,double(Num_mT_bins_pL)-0.5,512,0.0,0.4);

    TH2F* pL_NLO_mT_sor0 = new TH2F("pL_NLO_mT_sor0","pL_NLO_mT_sor0",Num_mT_bins_pL,-0.5,double(Num_mT_bins_pL)-0.5,512,0.2,1.8);
    TH2F* pL_NLO_mT_esor0 = new TH2F("pL_NLO_mT_esor0","pL_NLO_mT_esor0",Num_mT_bins_pL,-0.5,double(Num_mT_bins_pL)-0.5,512,0.0,0.4);

    TH2F* pL_mT_sor0 = new TH2F("pL_mT_sor0","pL_mT_sor0",Num_mT_bins_pL,-0.5,double(Num_mT_bins_pL)-0.5,512,0.2,1.8);
    TH2F* pL_mT_esor0 = new TH2F("pL_mT_esor0","pL_mT_esor0",Num_mT_bins_pL,-0.5,double(Num_mT_bins_pL)-0.5,512,0.0,0.4);

    TH2F* pL_LO_mT_sor1 = new TH2F("pL_LO_mT_sor1","pL_LO_mT_sor1",Num_mT_bins_pL,-0.5,double(Num_mT_bins_pL)-0.5,512,0.2,1.8);
    TH2F* pL_LO_mT_esor1 = new TH2F("pL_LO_mT_esor1","pL_LO_mT_esor1",Num_mT_bins_pL,-0.5,double(Num_mT_bins_pL)-0.5,512,0.0,0.4);

    TH2F* pL_NLO_mT_sor1 = new TH2F("pL_NLO_mT_sor1","pL_NLO_mT_sor1",Num_mT_bins_pL,-0.5,double(Num_mT_bins_pL)-0.5,512,0.2,1.8);
    TH2F* pL_NLO_mT_esor1 = new TH2F("pL_NLO_mT_esor1","pL_NLO_mT_esor1",Num_mT_bins_pL,-0.5,double(Num_mT_bins_pL)-0.5,512,0.0,0.4);

    TH2F* pL_mT_sor1 = new TH2F("pL_mT_sor1","pL_mT_sor1",Num_mT_bins_pL,-0.5,double(Num_mT_bins_pL)-0.5,512,0.2,1.8);
    TH2F* pL_mT_esor1 = new TH2F("pL_mT_esor1","pL_mT_esor1",Num_mT_bins_pL,-0.5,double(Num_mT_bins_pL)-0.5,512,0.0,0.4);

    TFile* ntFile_pL = new TFile(ntFileName_pL,"read");
    TNtuple* ntResult_pL = (TNtuple*)ntFile_pL->Get("ntResult");
    unsigned NumNtEntries_pL = ntResult_pL->GetEntries();

    Float_t IterM_pL;
    Float_t Iter_pL;
    Float_t Config_pL;
    Float_t Data_pL;
    Float_t SourceType_pL;
    Float_t SourceScale_pL;
    Float_t Potential_pL;
    Float_t Baseline_pL;
    Float_t FemRan_pL;
    Float_t FitRan_pL;
    Float_t pFrac_pL;
    Float_t LamFrac_pL;
    Float_t kcVar_pL;
    Float_t mTbin_pL;
    Float_t FemtoMin_pL;
    Float_t FemtoMax_pL;
    Float_t BlMin_pL;
    Float_t BlMax_pL;
    Float_t p_a_pL;
    Float_t e_a_pL;
    Float_t p_b_pL;
    Float_t e_b_pL;
    Float_t p_c_pL;
    Float_t e_c_pL;
    Float_t p_Cl_pL;
    Float_t e_Cl_pL;
    Float_t p_kc_pL;
    Float_t e_kc_pL;
    Float_t p_sor0_pL;
    Float_t e_sor0_pL;
    Float_t p_sor1_pL;
    Float_t e_sor1_pL;
    Float_t chi2_pL;
    Float_t ndf_pL;
    Float_t chi2_312_pL;
    Float_t ndf_312_pL;

    ntResult_pL->SetBranchAddress("IterM",&IterM_pL);
    ntResult_pL->SetBranchAddress("Iter",&Iter_pL);
    ntResult_pL->SetBranchAddress("Config",&Config_pL);
    ntResult_pL->SetBranchAddress("Data",&Data_pL);
    ntResult_pL->SetBranchAddress("SourceType",&SourceType_pL);
    ntResult_pL->SetBranchAddress("SourceScale",&SourceScale_pL);
    ntResult_pL->SetBranchAddress("Potential",&Potential_pL);
    ntResult_pL->SetBranchAddress("Baseline",&Baseline_pL);
    ntResult_pL->SetBranchAddress("FemRan",&FemRan_pL);
    ntResult_pL->SetBranchAddress("FitRan",&FitRan_pL);
    ntResult_pL->SetBranchAddress("pFrac",&pFrac_pL);
    ntResult_pL->SetBranchAddress("LamFrac",&LamFrac_pL);
    ntResult_pL->SetBranchAddress("kcVar",&kcVar_pL);
    ntResult_pL->SetBranchAddress("mTbin",&mTbin_pL);
    ntResult_pL->SetBranchAddress("FemtoMin",&FemtoMin_pL);
    ntResult_pL->SetBranchAddress("FemtoMax",&FemtoMax_pL);
    ntResult_pL->SetBranchAddress("BlMin",&BlMin_pL);
    ntResult_pL->SetBranchAddress("BlMax",&BlMax_pL);
    ntResult_pL->SetBranchAddress("p_a",&p_a_pL);
    ntResult_pL->SetBranchAddress("e_a",&e_a_pL);
    ntResult_pL->SetBranchAddress("p_b",&p_b_pL);
    ntResult_pL->SetBranchAddress("e_b",&e_b_pL);
    ntResult_pL->SetBranchAddress("p_c",&p_c_pL);
    ntResult_pL->SetBranchAddress("e_c",&e_c_pL);
    ntResult_pL->SetBranchAddress("p_Cl",&p_Cl_pL);
    ntResult_pL->SetBranchAddress("e_Cl",&e_Cl_pL);
    ntResult_pL->SetBranchAddress("p_kc",&p_kc_pL);
    ntResult_pL->SetBranchAddress("e_kc",&e_kc_pL);
    ntResult_pL->SetBranchAddress("p_sor0",&p_sor0_pL);
    ntResult_pL->SetBranchAddress("e_sor0",&e_sor0_pL);
    ntResult_pL->SetBranchAddress("p_sor1",&p_sor1_pL);
    ntResult_pL->SetBranchAddress("e_sor1",&e_sor1_pL);
    ntResult_pL->SetBranchAddress("chi2",&chi2_pL);
    ntResult_pL->SetBranchAddress("ndf",&ndf_pL);
    ntResult_pL->SetBranchAddress("chi2_312",&chi2_312_pL);
    ntResult_pL->SetBranchAddress("ndf_312",&ndf_312_pL);
printf("NumNtEntries_pL=%u\n",NumNtEntries_pL);
    for(unsigned uEntry=0; uEntry<NumNtEntries_pL; uEntry++){
        ntFile_pL->cd();
        ntResult_pL->GetEntry(uEntry);
        //int Config;
        //Config = int(Iter)%100;
        //int ITER = int(Iter)/100;
        int Id_LO;
        int Id_NLO;
        if(Config_pL<10){
            Id_LO = 1;
            Id_NLO = 11;
        }
        else if(Config_pL<100){
            Id_LO = 0;
            Id_NLO = 10;
        }
//printf("Iter=%i; Config=%i; LO=%i; NLO=%i\n",ITER,Config,Id_LO,Id_NLO);
        //! Conditions

        //if(Config!=0&&Config!=2) continue; //select Gauss+Reso and Gauss
        //if(Config!=20) continue; //select Gauss+Reso
        //if(Config!=1) continue; //select Levy+Reso
        //if(Config!=2) continue; //select Gauss
        //if(round(Baseline)!=0) continue; //select only norm
        //if(int(Baseline)!=11) continue; //select only pol1
        //if(int(Baseline)!=12) continue; //select only pol2
        if(int(Baseline_pL)==11) continue; // select pol0 and pol2
        //if(SourceScale>0) continue;
        //if(SourceScale<0.8&&SourceScale>0) continue;

        //if(mTbin<0) continue;

        if(Potential_pL==1){
            pL_LO_mT_sor0->Fill(mTbin_pL,p_sor0_pL);
            pL_LO_mT_esor0->Fill(mTbin_pL,e_sor0_pL);
//printf("mTbin,p_sor0 = %f %f\n",mTbin,p_sor0);
        }
        if(Potential_pL==11){
            pL_NLO_mT_sor0->Fill(mTbin_pL,p_sor0_pL);
            pL_NLO_mT_esor0->Fill(mTbin_pL,e_sor0_pL);
        }
        pL_mT_sor0->Fill(mTbin_pL,p_sor0_pL);
        pL_mT_esor0->Fill(mTbin_pL,e_sor0_pL);
        pL_mT_sor1->Fill(mTbin_pL,p_sor1_pL);
        pL_mT_esor1->Fill(mTbin_pL,e_sor1_pL);
    }

    TGraphErrors mTSyst_pL_LO;
    mTSyst_pL_LO.SetName("mTSyst_pL_LO");
    mTSyst_pL_LO.Set(Num_mT_bins_pL);
    mTSyst_pL_LO.SetFillColorAlpha(kGreen+3,0.3);
    mTSyst_pL_LO.SetLineColor(kGreen+3);
    mTSyst_pL_LO.SetLineWidth(3);

    TGraphErrors mTSyst_pL_NLO;
    mTSyst_pL_NLO.SetName("mTSyst_pL_NLO");
    mTSyst_pL_NLO.Set(Num_mT_bins_pL);
    mTSyst_pL_NLO.SetFillColorAlpha(kRed+1,0.3);
    mTSyst_pL_NLO.SetLineColor(kRed+1);
    mTSyst_pL_NLO.SetLineWidth(3);

    TGraphErrors mTSyst_pL;
    mTSyst_pL.SetName("mTSyst_pL");
    mTSyst_pL.Set(Num_mT_bins_pL);
    mTSyst_pL.SetFillColorAlpha(kOrange+4,0.3);
    mTSyst_pL.SetLineColor(kOrange+4);
    mTSyst_pL.SetLineWidth(3);

    TGraphErrors mTSyst_pp;
    mTSyst_pp.SetName("mTSyst_pp");
    mTSyst_pp.Set(Num_mT_bins_pp);
    mTSyst_pp.SetFillColorAlpha(kBlue+2,0.3);
    mTSyst_pp.SetLineColor(kBlue+2);
    mTSyst_pp.SetLineWidth(3);

    TGraphErrors mTStat_pL_LO;
    mTStat_pL_LO.SetName("mTStat_pL_LO");
    mTStat_pL_LO.Set(Num_mT_bins_pL);
    mTStat_pL_LO.SetFillColorAlpha(kGreen+3,0.3);
    mTStat_pL_LO.SetLineColor(kGreen+3);
    mTStat_pL_LO.SetLineWidth(3);

    TGraphErrors mTStat_pL_NLO;
    mTStat_pL_NLO.SetName("mTStat_pL_NLO");
    mTStat_pL_NLO.Set(Num_mT_bins_pL);
    mTStat_pL_NLO.SetFillColorAlpha(kRed+1,0.3);
    mTStat_pL_NLO.SetLineColor(kRed+1);
    mTStat_pL_NLO.SetLineWidth(3);

    TGraphErrors mTStat_pL;
    mTStat_pL.SetName("mTStat_pL");
    mTStat_pL.Set(Num_mT_bins_pL);
    mTStat_pL.SetFillColorAlpha(kOrange+4,0.3);
    mTStat_pL.SetLineColor(kOrange+4);
    mTStat_pL.SetLineWidth(3);

    TGraphErrors mTStat_pp;
    mTStat_pp.SetName("mTStat_pp");
    mTStat_pp.Set(Num_mT_bins_pp);
    mTStat_pp.SetFillColorAlpha(kBlue+2,0.3);
    mTStat_pp.SetLineColor(kBlue+2);
    mTStat_pp.SetLineWidth(3);

    TGraphErrors mTSystSor1_pL_LO;
    mTSystSor1_pL_LO.SetName("mTSystSor1_pL_LO");
    mTSystSor1_pL_LO.Set(Num_mT_bins_pL);
    mTSystSor1_pL_LO.SetFillColorAlpha(kGreen+3,0.3);
    mTSystSor1_pL_LO.SetLineColor(kGreen+3);
    mTSystSor1_pL_LO.SetLineWidth(3);

    TGraphErrors mTSystSor1_pL_NLO;
    mTSystSor1_pL_NLO.SetName("mTSystSor1_pL_NLO");
    mTSystSor1_pL_NLO.Set(Num_mT_bins_pL);
    mTSystSor1_pL_NLO.SetFillColorAlpha(kRed+1,0.3);
    mTSystSor1_pL_NLO.SetLineColor(kRed+1);
    mTSystSor1_pL_NLO.SetLineWidth(3);

    TGraphErrors mTSystSor1_pL;
    mTSystSor1_pL.SetName("mTSystSor1_pL");
    mTSystSor1_pL.Set(Num_mT_bins_pL);
    mTSystSor1_pL.SetFillColorAlpha(kOrange+4,0.3);
    mTSystSor1_pL.SetLineColor(kOrange+4);
    mTSystSor1_pL.SetLineWidth(3);

    TGraphErrors mTSystSor1_pp;
    mTSystSor1_pp.SetName("mTSystSor1_pp");
    mTSystSor1_pp.Set(Num_mT_bins_pp);
    mTSystSor1_pp.SetFillColorAlpha(kBlue+2,0.3);
    mTSystSor1_pp.SetLineColor(kBlue+2);
    mTSystSor1_pp.SetLineWidth(3);

    TGraphErrors mTStatSor1_pL_LO;
    mTStatSor1_pL_LO.SetName("mTStatSor1_pL_LO");
    mTStatSor1_pL_LO.Set(Num_mT_bins_pL);
    mTStatSor1_pL_LO.SetFillColorAlpha(kGreen+3,0.3);
    mTStatSor1_pL_LO.SetLineColor(kGreen+3);
    mTStatSor1_pL_LO.SetLineWidth(3);

    TGraphErrors mTStatSor1_pL_NLO;
    mTStatSor1_pL_NLO.SetName("mTStatSor1_pL_NLO");
    mTStatSor1_pL_NLO.Set(Num_mT_bins_pL);
    mTStatSor1_pL_NLO.SetFillColorAlpha(kRed+1,0.3);
    mTStatSor1_pL_NLO.SetLineColor(kRed+1);
    mTStatSor1_pL_NLO.SetLineWidth(3);

    TGraphErrors mTStatSor1_pL;
    mTStatSor1_pL.SetName("mTStatSor1_pL");
    mTStatSor1_pL.Set(Num_mT_bins_pL);
    mTStatSor1_pL.SetFillColorAlpha(kOrange+4,0.3);
    mTStatSor1_pL.SetLineColor(kOrange+4);
    mTStatSor1_pL.SetLineWidth(3);

    TGraphErrors mTStatSor1_pp;
    mTStatSor1_pp.SetName("mTStatSor1_pp");
    mTStatSor1_pp.Set(Num_mT_bins_pp);
    mTStatSor1_pp.SetFillColorAlpha(kBlue+2,0.3);
    mTStatSor1_pp.SetLineColor(kBlue+2);
    mTStatSor1_pp.SetLineWidth(3);

    for(unsigned umT=0; umT<Num_mT_bins_pp; umT++){
        double AvgmT;
        double AvgmTerr;
        double Dummy;

        TH1D* sor = pp_mT_sor0->ProjectionY("sor",umT+1,umT+1);
        TH1D* esor = pp_mT_esor0->ProjectionY("esor",umT+1,umT+1);
        double StatErr = esor->GetMean();
        double SystErr = sor->GetStdDev();

        AveragemT_pp->GetPoint(umT,Dummy,AvgmT);
        AvgmTerr = AveragemT_pp->GetErrorY(umT);

        mTSyst_pp.SetPoint(umT,AvgmT,sor->GetMean());//
        mTSyst_pp.SetPointError(umT,AvgmTerr,SystErr);

        mTStat_pp.SetPoint(umT,AvgmT,sor->GetMean());//
        mTStat_pp.SetPointError(umT,0,StatErr);

        TH1D* sor1 = pp_mT_sor1->ProjectionY("sor1",umT+1,umT+1);
        TH1D* esor1 = pp_mT_esor1->ProjectionY("esor1",umT+1,umT+1);
        double StatErrSor1 = esor1->GetMean();
        double SystErrSor1 = sor1->GetStdDev();

        //AveragemT_pp->GetPoint(umT,Dummy,AvgmT);
        //AvgmTerr = AveragemT_pp->GetErrorY(umT);

        mTSystSor1_pp.SetPoint(umT,AvgmT,sor1->GetMean());//
        mTSystSor1_pp.SetPointError(umT,AvgmTerr,SystErrSor1);

        mTStatSor1_pp.SetPoint(umT,AvgmT,sor1->GetMean());//
        mTStatSor1_pp.SetPointError(umT,0,StatErrSor1);
    }

    for(unsigned umT=0; umT<Num_mT_bins_pL; umT++){
        double AvgmT;
        double AvgmTerr;
        double Dummy;

        TH1D* sorLO = pL_LO_mT_sor0->ProjectionY("sorLO",umT+1,umT+1);
        TH1D* sorNLO = pL_NLO_mT_sor0->ProjectionY("sorNLO",umT+1,umT+1);
        TH1D* sor = pL_mT_sor0->ProjectionY("sor",umT+1,umT+1);
        TH1D* esorLO = pL_LO_mT_esor0->ProjectionY("esorLO",umT+1,umT+1);
        TH1D* esorNLO = pL_NLO_mT_esor0->ProjectionY("esorNLO",umT+1,umT+1);
        TH1D* esor = pL_mT_esor0->ProjectionY("esor",umT+1,umT+1);
        double StatErrLO = esorLO->GetMean();
        double SystErrLO = sorLO->GetStdDev();
        double StatErrNLO = esorNLO->GetMean();
        double SystErrNLO = sorNLO->GetStdDev();
        double StatErr = esor->GetMean();
        double SystErr = sor->GetStdDev();

        AveragemT_pL->GetPoint(umT,Dummy,AvgmT);
        AvgmTerr = AveragemT_pL->GetErrorY(umT);

        mTSyst_pL_LO.SetPoint(umT,AvgmT,sorLO->GetMean());//
        mTSyst_pL_LO.SetPointError(umT,AvgmTerr,SystErrLO);

        mTSyst_pL_NLO.SetPoint(umT,AvgmT,sorNLO->GetMean());//
        mTSyst_pL_NLO.SetPointError(umT,AvgmTerr,SystErrNLO);

        mTSyst_pL.SetPoint(umT,AvgmT,sor->GetMean());//
        mTSyst_pL.SetPointError(umT,AvgmTerr,SystErr);

//printf("sorLO->GetMean()=%f; sorNLO->GetMean()=%f\n",sorLO->GetMean(),sorNLO->GetMean());
        mTStat_pL_LO.SetPoint(umT,AvgmT,sorLO->GetMean());//
        mTStat_pL_LO.SetPointError(umT,0,StatErrLO);
        mTStat_pL_NLO.SetPoint(umT,AvgmT,sorNLO->GetMean());//
        mTStat_pL_NLO.SetPointError(umT,0,StatErrNLO);
        mTStat_pL.SetPoint(umT,AvgmT,sor->GetMean());//
        mTStat_pL.SetPointError(umT,0,StatErr);

        TH1D* sor1LO = pL_LO_mT_sor1->ProjectionY("sor1LO",umT+1,umT+1);
        TH1D* sor1NLO = pL_NLO_mT_sor1->ProjectionY("sor1NLO",umT+1,umT+1);
        TH1D* sor1 = pL_mT_sor1->ProjectionY("sor1",umT+1,umT+1);
        TH1D* esor1LO = pL_LO_mT_esor1->ProjectionY("esor1LO",umT+1,umT+1);
        TH1D* esor1NLO = pL_NLO_mT_esor1->ProjectionY("esor1NLO",umT+1,umT+1);
        TH1D* esor1 = pL_mT_esor1->ProjectionY("esor1",umT+1,umT+1);
        double StatErrLOSor1 = esor1LO->GetMean();
        double SystErrLOSor1 = sor1LO->GetStdDev();
        double StatErrNLOSor1 = esor1NLO->GetMean();
        double SystErrNLOSor1 = sor1NLO->GetStdDev();
        double StatErrSor1 = esor1->GetMean();
        double SystErrSor1 = sor1->GetStdDev();

        mTSystSor1_pL_LO.SetPoint(umT,AvgmT,sor1LO->GetMean());//
        mTSystSor1_pL_LO.SetPointError(umT,AvgmTerr,SystErrLOSor1);

        mTSystSor1_pL_NLO.SetPoint(umT,AvgmT,sor1NLO->GetMean());//
        mTSystSor1_pL_NLO.SetPointError(umT,AvgmTerr,SystErrNLOSor1);

        mTSystSor1_pL.SetPoint(umT,AvgmT,sor1->GetMean());//
        mTSystSor1_pL.SetPointError(umT,AvgmTerr,SystErrSor1);

printf("umT=%u; sor1->GetMean()=%f\n",umT,sor1->GetMean());
        mTStatSor1_pL_LO.SetPoint(umT,AvgmT,sor1LO->GetMean());//
        mTStatSor1_pL_LO.SetPointError(umT,0,StatErrLOSor1);
        mTStatSor1_pL_NLO.SetPoint(umT,AvgmT,sor1NLO->GetMean());//
        mTStatSor1_pL_NLO.SetPointError(umT,0,StatErrNLOSor1);
        mTStatSor1_pL.SetPoint(umT,AvgmT,sor1->GetMean());//
        mTStatSor1_pL.SetPointError(umT,0,StatErrSor1);
    }

    DLM_Histo<double> DimiStat_pp;
    DLM_Histo<double> DimiStat_pL_LO;
    DLM_Histo<double> DimiStat_pL_NLO;
    DLM_Histo<double> DimiStat_pL;

    DLM_Histo<double> DimiStatEff_pL_NLO;
    DLM_Histo<double> DimiStatRatio_pL_NLO;

    DLM_Histo<double> DimiStatEff_pL;
    DLM_Histo<double> DimiStatRatio_pL;

    double xValLow;
    double yValLow;
    double xValUp;
    double yValUp;

    double xVal;
    double yVal;

    DimiStat_pp.SetUp(1);
    mTStat_pp.GetPoint(0,xValLow,yValLow);
    mTStat_pp.GetPoint(0,xValUp,yValUp);
    DimiStat_pp.SetUp(0,mTStat_pp.GetN(),xValLow-0.001,xValUp+0.001);
    DimiStat_pp.Initialize();
printf("mT Radius Stat Syst\n");
printf(" pp AV18:\n");
    for(unsigned uBin=0; uBin<mTStat_pp.GetN(); uBin++){
        mTStat_pp.GetPoint(uBin,xVal,yVal);
        DimiStat_pp.SetBinCenter(0,uBin,xVal);
        DimiStat_pp.SetBinContent(uBin,yVal);
        DimiStat_pp.SetBinError(uBin,mTStat_pp.GetErrorY(uBin));
printf("  %.3f   %.3f   %.3f   %.3f\n",xVal,yVal,mTStat_pp.GetErrorY(uBin),mTSyst_pp.GetErrorY(uBin));
//printf(" DimiStat_pp: cen/con/err = %.3f/%.3f/%.3f\n",DimiStat_pp.GetBinCenter(0,uBin),DimiStat_pp.GetBinContent(uBin),DimiStat_pp.GetBinError(uBin));
    }

    DimiStat_pL_LO.SetUp(1);
    mTStat_pL_LO.GetPoint(0,xValLow,yValLow);
    mTStat_pL_LO.GetPoint(0,xValUp,yValUp);
    DimiStat_pL_LO.SetUp(0,mTStat_pL_LO.GetN(),xValLow-0.001,xValUp+0.001);
    DimiStat_pL_LO.Initialize();
printf("mT Radius Stat Syst\n");
printf(" pL LO:\n");
    for(unsigned uBin=0; uBin<mTStat_pL_LO.GetN(); uBin++){
        mTStat_pL_LO.GetPoint(uBin,xVal,yVal);
        DimiStat_pL_LO.SetBinCenter(0,uBin,xVal);
        DimiStat_pL_LO.SetBinContent(uBin,yVal);
        DimiStat_pL_LO.SetBinError(uBin,mTStat_pL_LO.GetErrorY(uBin));
printf("  %.3f   %.3f   %.3f   %.3f\n",xVal,yVal,mTStat_pL_LO.GetErrorY(uBin),mTSyst_pL_LO.GetErrorY(uBin));
    }

    DimiStat_pL_NLO.SetUp(1);
    mTStat_pL_NLO.GetPoint(0,xValLow,yValLow);
    mTStat_pL_NLO.GetPoint(0,xValUp,yValUp);
    DimiStat_pL_NLO.SetUp(0,mTStat_pL_NLO.GetN(),xValLow-0.001,xValUp+0.001);
    DimiStat_pL_NLO.Initialize();
printf("mT Radius Stat Syst\n");
printf(" pL NLO:\n");
    for(unsigned uBin=0; uBin<mTStat_pL_NLO.GetN(); uBin++){
        mTStat_pL_NLO.GetPoint(uBin,xVal,yVal);
        DimiStat_pL_NLO.SetBinCenter(0,uBin,xVal);
        DimiStat_pL_NLO.SetBinContent(uBin,yVal);
        DimiStat_pL_NLO.SetBinError(uBin,mTStat_pL_NLO.GetErrorY(uBin));
printf("  %.3f   %.3f   %.3f   %.3f\n",xVal,yVal,mTStat_pL_NLO.GetErrorY(uBin),mTSyst_pL_NLO.GetErrorY(uBin));
//printf(" cen/con/err = %.3f/%.3f/%.3f\n",DimiStat_pL_NLO.GetBinCenter(0,uBin),DimiStat_pL_NLO.GetBinContent(uBin),DimiStat_pL_NLO.GetBinError(uBin));
    }


    DimiStatEff_pL_NLO.Copy(DimiStat_pp);
printf("DimiStatEff_pL_NLO.GetNbins()=%u\n",DimiStatEff_pL_NLO.GetNbins());
    for(unsigned uBin=0; uBin<DimiStatEff_pL_NLO.GetNbins(); uBin++){
        double mT_VALUE = DimiStatEff_pL_NLO.GetBinCenter(0,uBin);
        DimiStatEff_pL_NLO.SetBinContent(uBin,DimiStat_pL_NLO.Eval(&mT_VALUE));
        DimiStatEff_pL_NLO.SetBinError(uBin,DimiStat_pL_NLO.EvalError(&mT_VALUE));
//printf(" NLO: mT_VALUE=%.3f; Value=%.3f+/-%.3f\n",mT_VALUE,DimiStat_pL_NLO.Eval(&mT_VALUE),DimiStat_pL_NLO.EvalError(&mT_VALUE));
    }

    DimiStatRatio_pL_NLO.Copy(DimiStatEff_pL_NLO);
    DimiStatRatio_pL_NLO.DivideHisto(DimiStat_pp);

    TGraphErrors gStatRatio_pL_NLO;
    gStatRatio_pL_NLO.SetName("gStatRatio_pL_NLO");
    gStatRatio_pL_NLO.Set(DimiStatRatio_pL_NLO.GetNbins());
    for(unsigned uBin=0; uBin<DimiStatRatio_pL_NLO.GetNbins(); uBin++){
        gStatRatio_pL_NLO.SetPoint(uBin,DimiStatRatio_pL_NLO.GetBinCenter(0,uBin),DimiStatRatio_pL_NLO.GetBinContent(uBin));
        gStatRatio_pL_NLO.SetPointError(uBin,mTSyst_pp.GetErrorX(uBin),DimiStatRatio_pL_NLO.GetBinError(uBin));
    }


    DimiStat_pL.SetUp(1);
    mTStat_pL.GetPoint(0,xValLow,yValLow);
    mTStat_pL.GetPoint(0,xValUp,yValUp);
    DimiStat_pL.SetUp(0,mTStat_pL.GetN(),xValLow-0.001,xValUp+0.001);
    DimiStat_pL.Initialize();
printf("mT Radius Stat Syst\n");
printf(" pL NLO:\n");
    for(unsigned uBin=0; uBin<mTStat_pL.GetN(); uBin++){
        mTStat_pL.GetPoint(uBin,xVal,yVal);
        DimiStat_pL.SetBinCenter(0,uBin,xVal);
        DimiStat_pL.SetBinContent(uBin,yVal);
        DimiStat_pL.SetBinError(uBin,mTStat_pL.GetErrorY(uBin));
printf("  %.3f   %.3f   %.3f   %.3f\n",xVal,yVal,mTStat_pL.GetErrorY(uBin),mTSyst_pL.GetErrorY(uBin));
//printf(" cen/con/err = %.3f/%.3f/%.3f\n",DimiStat_pL.GetBinCenter(0,uBin),DimiStat_pL.GetBinContent(uBin),DimiStat_pL.GetBinError(uBin));
    }


    DimiStatEff_pL.Copy(DimiStat_pp);
printf("DimiStatEff_pL.GetNbins()=%u\n",DimiStatEff_pL.GetNbins());
    for(unsigned uBin=0; uBin<DimiStatEff_pL.GetNbins(); uBin++){
        double mT_VALUE = DimiStatEff_pL.GetBinCenter(0,uBin);
        DimiStatEff_pL.SetBinContent(uBin,DimiStat_pL.Eval(&mT_VALUE));
        DimiStatEff_pL.SetBinError(uBin,DimiStat_pL.EvalError(&mT_VALUE));
//printf(" NLO: mT_VALUE=%.3f; Value=%.3f+/-%.3f\n",mT_VALUE,DimiStat_pL.Eval(&mT_VALUE),DimiStat_pL.EvalError(&mT_VALUE));
    }


    DimiStatRatio_pL.Copy(DimiStatEff_pL);
    DimiStatRatio_pL.DivideHisto(DimiStat_pp);

    TGraphErrors gStatRatio_pL;
    gStatRatio_pL.SetName("gStatRatio_pL");
    gStatRatio_pL.Set(DimiStatRatio_pL.GetNbins());
    for(unsigned uBin=0; uBin<DimiStatRatio_pL.GetNbins(); uBin++){
        gStatRatio_pL.SetPoint(uBin,DimiStatRatio_pL.GetBinCenter(0,uBin),DimiStatRatio_pL.GetBinContent(uBin));
        gStatRatio_pL.SetPointError(uBin,mTSyst_pp.GetErrorX(uBin),DimiStatRatio_pL.GetBinError(uBin));
    }

    TFile* fOut = new TFile(OutputFolder+"fOut.root","recreate");

    TH1F* hAxis = new TH1F("hAxis", "hAxis", 128, 0.9, 2.7);
    hAxis->SetStats(false);
    hAxis->SetTitle("");
    hAxis->GetXaxis()->SetTitle("<#it{m}_{T}> (GeV/#it{c}^{2})");
    hAxis->GetXaxis()->SetTitleSize(0.06);
    hAxis->GetXaxis()->SetLabelSize(0.06);
    //hAxis->GetXaxis()->CenterTitle();
    hAxis->GetXaxis()->SetTitleOffset(1.3);
    hAxis->GetXaxis()->SetLabelOffset(0.02);

    //if(LevySource) hAxis->GetYaxis()->SetTitle("Levy core (fm)");
    //else hAxis->GetYaxis()->SetTitle("Gaussian core (fm)");
    TString FileDescr;
    if(Descr1.Contains("Levy core")){
        hAxis->GetYaxis()->SetTitle("r_{core} (fm)");
        hAxis->GetYaxis()->SetRangeUser(0.55, 1.65);
        FileDescr = "LevyReso";
    }
    else if(Descr1.Contains("Gaussian core")&&Descr2!="Random angle"){
        hAxis->GetYaxis()->SetTitle("r_{core} (fm)");
        hAxis->GetYaxis()->SetRangeUser(0.55, 1.35);
        FileDescr = "GaussReso";
    }
    else if(Descr2=="Random angle"){
        hAxis->GetYaxis()->SetTitle("r_{core} (fm)");
        hAxis->GetYaxis()->SetRangeUser(0.75, 1.65);
        FileDescr = "RandomAngle";
    }
    else{
        hAxis->GetYaxis()->SetTitle("r_{0} (fm)");
        hAxis->GetYaxis()->SetRangeUser(0.90, 1.80);
        FileDescr = "Gauss";
    }
    //
    hAxis->GetYaxis()->SetTitleSize(0.06);
    hAxis->GetYaxis()->SetLabelSize(0.06);
    //hAxis->GetYaxis()->CenterTitle();
    hAxis->GetYaxis()->SetTitleOffset(0.90);

    //hAxis->GetXaxis()->SetNdivisions(506);

    //

    //TLegend* lLegend = new TLegend(0.45,0.85,0.95,0.95);//lbrt
    //lLegend->SetName(TString::Format("lLegend"));
    //lLegend->SetTextSize(0.04);
    //lLegend->SetNColumns(3);
    TLegend* lLegend = new TLegend(0.67,0.48,0.90,0.71);//lbrt
    lLegend->SetName(TString::Format("lLegend"));
    lLegend->SetTextSize(0.04);
    lLegend->SetBorderSize(0);
    lLegend->AddEntry(&mTStat_pp,"p#minusp (AV18)");
    lLegend->AddEntry(&mTSyst_pL_LO,"p#minus#Lambda (LO)");
    lLegend->AddEntry(&mTSyst_pL_NLO,"p#minus#Lambda (NLO)");
    //if(mT_pXim)lLegend->AddEntry(mT_pXim,"p#minus#Xi^{#minus} (HAL QCD)");

    TCanvas* cmT = new TCanvas("cmT", "cmT", 1);
    cmT->cd(0); cmT->SetCanvasSize(1920/2, 1080/2); cmT->SetMargin(0.15,0.05,0.2,0.05);//lrbt

    hAxis->Draw("axis");
    mTSyst_pL_NLO.Draw("e2,same");
    mTSyst_pL_LO.Draw("e2,same");
    mTSyst_pp.Draw("e2,same");
    mTStat_pL_NLO.Draw("e,same");
    mTStat_pL_LO.Draw("e,same");
    mTStat_pp.Draw("e,same");
    lLegend->Draw("same");

    //TLatex BeamText;
    //BeamText.SetTextSize(gStyle->GetTextSize()*0.90);
    //BeamText.SetTextFont(42);
    //BeamText.SetNDC(kTRUE);
    //BeamText.DrawLatex(0.68, 0.76, "ALICE Preliminary");
    //BeamText.DrawLatex(0.68, 0.71, "high-mult. (0-0.072% INEL) pp #sqrt{#it{s}} = 13 TeV");

    TLatex BeamText;
    BeamText.SetTextSize(0.04);
    BeamText.SetTextFont(42);
    BeamText.SetNDC(kTRUE);
    BeamText.DrawLatex(0.677, 0.90, "ALICE Preliminary");
    TLatex BeamText1;
    BeamText1.SetTextSize(0.04);
    BeamText1.SetTextFont(42);
    BeamText1.SetNDC(kTRUE);
    BeamText1.DrawLatex(0.677, 0.85, "high-mult. (0-0.072% INEL)");
    BeamText1.DrawLatex(0.677, 0.81, "pp #sqrt{#it{s}} = 13 TeV");

    TLatex BeamText2;
    BeamText2.SetTextSize(0.035);
    BeamText2.SetTextFont(42);
    BeamText2.SetNDC(kTRUE);
    if(Descr2!=""){
        BeamText2.DrawLatex(0.677, 0.76, Descr1);
        BeamText2.DrawLatex(0.677, 0.72, Descr2);
    }
    else{
        BeamText2.DrawLatex(0.677, 0.74, Descr1);
    }

    TLegend* lLegend_TOT = new TLegend(0.67,0.48,0.90,0.71);//lbrt
    lLegend_TOT->SetName(TString::Format("lLegend_TOT"));
    lLegend_TOT->SetTextSize(0.04);
    lLegend_TOT->SetBorderSize(0);
    lLegend_TOT->AddEntry(&mTStat_pp,"p#minusp (AV18)");
    lLegend_TOT->AddEntry(&mTSyst_pL,"p#minus#Lambda (LO,NLO,Usmani)");
    //if(mT_pXim)lLegend->AddEntry(mT_pXim,"p#minus#Xi^{#minus} (HAL QCD)");

    TCanvas* cmT_TOT = new TCanvas("cmT_TOT", "cmT_TOT", 1);
    cmT_TOT->cd(0); cmT_TOT->SetCanvasSize(1920/2, 1080/2); cmT_TOT->SetMargin(0.15,0.05,0.2,0.05);//lrbt

    hAxis->Draw("axis");
    mTSyst_pL.Draw("e2,same");
    mTSyst_pp.Draw("e2,same");
    mTStat_pL.Draw("e,same");
    mTStat_pp.Draw("e,same");
    lLegend_TOT->Draw("same");

    //TLatex BeamText;
    //BeamText.SetTextSize(gStyle->GetTextSize()*0.90);
    //BeamText.SetTextFont(42);
    //BeamText.SetNDC(kTRUE);
    //BeamText.DrawLatex(0.68, 0.76, "ALICE Preliminary");
    //BeamText.DrawLatex(0.68, 0.71, "high-mult. (0-0.072% INEL) pp #sqrt{#it{s}} = 13 TeV");

    //TLatex BeamText;
    //BeamText.SetTextSize(0.04);
    //BeamText.SetTextFont(42);
    //BeamText.SetNDC(kTRUE);
    BeamText.DrawLatex(0.677, 0.90, "ALICE Preliminary");
    //TLatex BeamText1;
    //BeamText1.SetTextSize(0.04);
    //BeamText1.SetTextFont(42);
    //BeamText1.SetNDC(kTRUE);
    BeamText1.DrawLatex(0.677, 0.85, "high-mult. (0-0.072% INEL)");
    BeamText1.DrawLatex(0.677, 0.81, "pp #sqrt{#it{s}} = 13 TeV");

    //TLatex BeamText2;
    //BeamText2.SetTextSize(0.035);
    //BeamText2.SetTextFont(42);
    //BeamText2.SetNDC(kTRUE);
    if(Descr2!=""){
        BeamText2.DrawLatex(0.677, 0.76, Descr1);
        BeamText2.DrawLatex(0.677, 0.72, Descr2);
    }
    else{
        BeamText2.DrawLatex(0.677, 0.74, Descr1);
    }

    fOut->cd();
    mTSyst_pL_LO.Write();
    mTSyst_pL_NLO.Write();
    mTSyst_pL.Write();
    mTStat_pL_LO.Write();
    mTStat_pL_NLO.Write();
    mTStat_pL.Write();
    mTSyst_pp.Write();
    mTStat_pp.Write();
    gStatRatio_pL_NLO.Write();
    gStatRatio_pL.Write();
    cmT_TOT->Write();
    cmT_TOT->SaveAs(OutputFolder+"cmT_"+FileDescr+"_TOT.pdf");




////////////////////////
    TH1F* hAxisSor1 = new TH1F("hAxisSor1", "hAxisSor1", 128, 0.9, 2.7);
    hAxisSor1->SetStats(false);
    hAxisSor1->SetTitle("");
    hAxisSor1->GetXaxis()->SetTitle("<#it{m}_{T}> (GeV/#it{c}^{2})");
    hAxisSor1->GetXaxis()->SetTitleSize(0.06);
    hAxisSor1->GetXaxis()->SetLabelSize(0.06);
    //hAxisSor1->GetXaxis()->CenterTitle();
    hAxisSor1->GetXaxis()->SetTitleOffset(1.3);
    hAxisSor1->GetXaxis()->SetLabelOffset(0.02);

    //if(LevySource) hAxis->GetYaxis()->SetTitle("Levy core (fm)");
    //else hAxis->GetYaxis()->SetTitle("Gaussian core (fm)");
    TString FileDescrSor1;
    if(Descr1.Contains("Levy core")){
        hAxisSor1->GetYaxis()->SetTitle("#alpha");
        hAxisSor1->GetYaxis()->SetRangeUser(0.8,2.65);
        FileDescr = "LevyReso";
    }
    else if(Descr1.Contains("Gaussian core")&&Descr2!="Random angle"){
        hAxisSor1->GetYaxis()->SetTitle("#alpha");
        hAxisSor1->GetYaxis()->SetRangeUser(0.8,2.65);
        FileDescrSor1 = "GaussReso";
    }
    else if(Descr2=="Random angle"){
        hAxisSor1->GetYaxis()->SetTitle("#alpha");
        hAxisSor1->GetYaxis()->SetRangeUser(0.85,2.65);
        FileDescrSor1 = "RandomAngle";
    }
    else{
        hAxisSor1->GetYaxis()->SetTitle("#alpha");
        hAxisSor1->GetYaxis()->SetRangeUser(0.85,2.65);
        FileDescrSor1 = "Gauss";
    }
    //
    hAxisSor1->GetYaxis()->SetTitleSize(0.06);
    hAxisSor1->GetYaxis()->SetLabelSize(0.06);
    //hAxisSor1->GetYaxis()->CenterTitle();
    hAxisSor1->GetYaxis()->SetTitleOffset(0.90);

    //hAxisSor1->GetXaxis()->SetNdivisions(506);

    //

    //TLegend* lLegend = new TLegend(0.45,0.85,0.95,0.95);//lbrt
    //lLegend->SetName(TString::Format("lLegend"));
    //lLegend->SetTextSize(0.04);
    //lLegend->SetNColumns(3);
    TLegend* lLegendSor1 = new TLegend(0.67,0.48,0.90,0.71);//lbrt
    lLegendSor1->SetName(TString::Format("lLegendSor1"));
    lLegendSor1->SetTextSize(0.04);
    lLegendSor1->SetBorderSize(0);
    lLegendSor1->AddEntry(&mTStat_pp,"p#minusp (AV18)");
    lLegendSor1->AddEntry(&mTSyst_pL_LO,"p#minus#Lambda (LO)");
    lLegendSor1->AddEntry(&mTSyst_pL_NLO,"p#minus#Lambda (NLO)");
    //if(mT_pXim)lLegendSor1->AddEntry(mT_pXim,"p#minus#Xi^{#minus} (HAL QCD)");

    TCanvas* cmTSor1 = new TCanvas("cmTSor1", "cmTSor1", 1);
    cmTSor1->cd(0); cmTSor1->SetCanvasSize(1920/2, 1080/2); cmTSor1->SetMargin(0.15,0.05,0.2,0.05);//lrbt

    hAxisSor1->Draw("axis");
    mTSystSor1_pL_NLO.Draw("e2,same");
    mTSystSor1_pL_LO.Draw("e2,same");
    mTSystSor1_pp.Draw("e2,same");
    mTStatSor1_pL_NLO.Draw("e,same");
    mTStatSor1_pL_LO.Draw("e,same");
    mTStatSor1_pp.Draw("e,same");
    lLegendSor1->Draw("same");

    //TLatex BeamText;
    //BeamText.SetTextSize(gStyle->GetTextSize()*0.90);
    //BeamText.SetTextFont(42);
    //BeamText.SetNDC(kTRUE);
    //BeamText.DrawLatex(0.68, 0.76, "ALICE Preliminary");
    //BeamText.DrawLatex(0.68, 0.71, "high-mult. (0-0.072% INEL) pp #sqrt{#it{s}} = 13 TeV");

    TLatex BeamTextSor1;
    BeamTextSor1.SetTextSize(0.04);
    BeamTextSor1.SetTextFont(42);
    BeamTextSor1.SetNDC(kTRUE);
    BeamTextSor1.DrawLatex(0.677, 0.90, "ALICE Preliminary");
    TLatex BeamText1Sor1;
    BeamText1Sor1.SetTextSize(0.04);
    BeamText1Sor1.SetTextFont(42);
    BeamText1Sor1.SetNDC(kTRUE);
    BeamText1Sor1.DrawLatex(0.677, 0.85, "high-mult. (0-0.072% INEL)");
    BeamText1Sor1.DrawLatex(0.677, 0.81, "pp #sqrt{#it{s}} = 13 TeV");

    TLatex BeamText2Sor1;
    BeamText2Sor1.SetTextSize(0.035);
    BeamText2Sor1.SetTextFont(42);
    BeamText2Sor1.SetNDC(kTRUE);
    if(Descr2!=""){
        BeamText2Sor1.DrawLatex(0.677, 0.76, Descr1);
        BeamText2Sor1.DrawLatex(0.677, 0.72, Descr2);
    }
    else{
        BeamText2Sor1.DrawLatex(0.677, 0.74, Descr1);
    }

    TLegend* lLegendSor1_TOT = new TLegend(0.67,0.70,0.90,0.93);//lbrt
    lLegendSor1_TOT->SetName(TString::Format("lLegendSor1_TOT"));
    lLegendSor1_TOT->SetTextSize(0.04);
    lLegendSor1_TOT->SetBorderSize(0);
    lLegendSor1_TOT->AddEntry(&mTStatSor1_pp,"p#minusp (AV18)");
    lLegendSor1_TOT->AddEntry(&mTSystSor1_pL,"p#minus#Lambda (LO,NLO,Usmani)");
    //if(mT_pXim)lLegend->AddEntry(mT_pXim,"p#minus#Xi^{#minus} (HAL QCD)");

    TCanvas* cmTSor1_TOT = new TCanvas("cmTSor1_TOT", "cmTSor1_TOT", 1);
    cmTSor1_TOT->cd(0); cmTSor1_TOT->SetCanvasSize(1920/2, 1080/2); cmTSor1_TOT->SetMargin(0.15,0.05,0.2,0.05);//lrbt

    hAxisSor1->Draw("axis");
    mTSystSor1_pL.Draw("e2,same");
    mTSystSor1_pp.Draw("e2,same");
    mTStatSor1_pL.Draw("e,same");
    mTStatSor1_pp.Draw("e,same");
    lLegendSor1_TOT->Draw("same");

    //TLatex BeamText;
    //BeamText.SetTextSize(gStyle->GetTextSize()*0.90);
    //BeamText.SetTextFont(42);
    //BeamText.SetNDC(kTRUE);
    //BeamText.DrawLatex(0.68, 0.76, "ALICE Preliminary");
    //BeamText.DrawLatex(0.68, 0.71, "high-mult. (0-0.072% INEL) pp #sqrt{#it{s}} = 13 TeV");

    //TLatex BeamText;
    //BeamText.SetTextSize(0.04);
    //BeamText.SetTextFont(42);
    //BeamText.SetNDC(kTRUE);
    BeamTextSor1.DrawLatex(0.377, 0.90, "ALICE Preliminary");
    //TLatex BeamText1;
    //BeamText1.SetTextSize(0.04);
    //BeamText1.SetTextFont(42);
    //BeamText1.SetNDC(kTRUE);
    BeamText1Sor1.DrawLatex(0.377, 0.85, "high-mult. (0-0.072% INEL)");
    BeamText1Sor1.DrawLatex(0.377, 0.81, "pp #sqrt{#it{s}} = 13 TeV");

    //TLatex BeamText2;
    //BeamText2.SetTextSize(0.035);
    //BeamText2.SetTextFont(42);
    //BeamText2.SetNDC(kTRUE);
    if(Descr2!=""){
        BeamText2Sor1.DrawLatex(0.377, 0.76, Descr1);
        BeamText2Sor1.DrawLatex(0.377, 0.72, Descr2);
    }
    else{
        BeamText2Sor1.DrawLatex(0.377, 0.74, Descr1);
    }

    fOut->cd();
    mTSystSor1_pL_LO.Write();
    mTSystSor1_pL_NLO.Write();
    mTSystSor1_pL.Write();
    mTStatSor1_pL_LO.Write();
    mTStatSor1_pL_NLO.Write();
    mTStatSor1_pL.Write();
    mTSystSor1_pp.Write();
    mTStatSor1_pp.Write();
    //gStatRatio_pL_NLO.Write();
    //gStatRatio_pL.Write();
    cmTSor1_TOT->Write();
    cmTSor1_TOT->SaveAs(OutputFolder+"cmTSor1_"+FileDescr+"_TOT.pdf");

    const double AvgMt_pL = 1.55;
    unsigned ClosestBin = 3;
    double ExpectedStat;
    double ExpectedSyst;
    double ExpectedTot;
    double LowestStat=100;
    double HighestStat=0;
    double LowestSyst=100;
    double HighestSyst=0;

    delete pL_LO_mT_sor0;
    delete pL_LO_mT_esor0;
    delete pL_NLO_mT_sor0;
    delete pL_NLO_mT_esor0;
    delete pL_mT_sor0;
    delete pL_mT_esor0;
    delete hAxis;
    delete lLegend;
    delete lLegend_TOT;
    delete cmT;
    delete cmT_TOT;

    delete pL_LO_mT_sor1;
    delete pL_LO_mT_esor1;
    delete pL_NLO_mT_sor1;
    delete pL_NLO_mT_esor1;
    delete pL_mT_sor1;
    delete pL_mT_esor1;
    delete hAxisSor1;
    delete lLegendSor1;
    delete lLegendSor1_TOT;
    delete cmTSor1;
    delete cmTSor1_TOT;

    delete ntFile_pp;
    delete ntFile_pL;
    delete mtFile_pp;
    delete mtFile_pL;


    delete fOut;
}

void ComputeDataSystematics(){
    const TString BaseName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample12HM/CFOutput_pL";
    const TString DataHistName = "hCk_ReweightedMeV_2";
    const unsigned NumVars = 45;
    const TString SystBaseName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample12HM/Systematics_pL";
    const TString SystRelHistName = "hRelSyst";
    const TString SystRelFitName = "SystErrRel";
    const unsigned NumCkBins = 1024;
    TString FileName;
    TFile* fData;
    TFile* fData0;
    TH1F* hData;
    TH1F* hData0;

    FileName = BaseName+"_0.root";
    fData0 = new TFile(FileName,"read");
    hData0 = (TH1F*)fData0->Get(DataHistName);
    const unsigned NumMomBins = hData0->GetNbinsX();
    //the absolute dev. from the default
    TH2F* hRelSyst2D = new TH2F("hRelSyst2D","hRelSyst2D",
                             NumMomBins,hData0->GetBinLowEdge(1),hData0->GetXaxis()->GetBinUpEdge(hData0->GetNbinsX()),
                             NumCkBins,0.9,1.1);

    for(unsigned uVar=0; uVar<NumVars; uVar++){
        FileName = TString::Format("%s_%u.root",BaseName.Data(),uVar);
        fData = new TFile(FileName,"read");
        hData = (TH1F*)fData->Get(DataHistName);
        for(unsigned uBin=0; uBin<NumMomBins; uBin++){
            hRelSyst2D->Fill(hData->GetBinCenter(uBin+1),hData->GetBinContent(uBin+1)/hData0->GetBinContent(uBin+1));
        }
        delete fData;
    }

    TH1F* hRelSyst = new TH1F(SystRelHistName,SystRelHistName,NumMomBins,hData0->GetBinLowEdge(1),hData0->GetXaxis()->GetBinUpEdge(hData0->GetNbinsX()));

    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        //TH1D* hProjection = hRelSyst2D->ProjectionY("hProjectionY",uBin+1,uBin+1);
        double MinValue = 1e6;
        double MaxValue = -1e6;
        double Mean=0;
        double Stdev=0;
        unsigned NumEntries=0;
        for(unsigned uCk=0; uCk<NumCkBins; uCk++){
            //if(hProjection->GetBinContent(uCk+1)<MinValue){MinValue=hProjection->GetBinContent(uCk+1);}
            //if(hProjection->GetBinContent(uCk+1)>MaxValue){MaxValue=hProjection->GetBinContent(uCk+1);}
            if(hRelSyst2D->GetBinContent(uBin+1,uCk+1)){
                if(hRelSyst2D->GetYaxis()->GetBinCenter(uCk+1)<MinValue){MinValue=hRelSyst2D->GetYaxis()->GetBinCenter(uCk+1);}
                if(hRelSyst2D->GetYaxis()->GetBinCenter(uCk+1)>MaxValue){MaxValue=hRelSyst2D->GetYaxis()->GetBinCenter(uCk+1);}
                Mean+=hRelSyst2D->GetYaxis()->GetBinCenter(uCk+1);
                Stdev+=pow(hRelSyst2D->GetYaxis()->GetBinCenter(uCk+1),2.);
                NumEntries++;
            }
        }
        Mean/=double(NumEntries);
        Stdev/=double(NumEntries);
        Stdev = sqrt(Stdev-Mean*Mean);

        //printf("uBin=%u; min=%f; max=%f\n",uBin,hRelSyst2D->ProjectionY("hProjectionY",uBin+1,uBin+1)->GetMinimum(),hRelSyst2D->ProjectionY("hProjectionY",uBin+1,uBin+1)->GetMaximum());
        //printf(" uBin=%u; min=%f; max=%f\n",uBin,MinValue,MaxValue);
        double Interval = MaxValue-MinValue;
        double UniformStdv = Interval/sqrt(12.);
        hRelSyst->SetBinContent(uBin+1,UniformStdv);
        hRelSyst->SetBinError(uBin+1,0.00001);
        //hRelSyst->SetBinContent(uBin+1,Stdev);
    }

    hRelSyst->Rebin(4);
    hRelSyst->Scale(0.25);

    TF1* fPol6 = new TF1(SystRelFitName,"[0]+[1]*x+[2]*pow(x,2.)+[3]*pow(x,3.)+[4]*pow(x,4.)+[5]*pow(x,5.)+[6]*pow(x,6.)",hData0->GetBinLowEdge(1),hData0->GetXaxis()->GetBinUpEdge(hData0->GetNbinsX()));
    fPol6->SetParameter(0,0.011);
    fPol6->SetParameter(1,-0.00015);
    fPol6->SetParameter(2,7.8e-7);
    fPol6->SetParameter(3,-2.12e-9);
    fPol6->SetParameter(4,3.1e-12);
    fPol6->SetParameter(5,-2.2e-15);
    fPol6->SetParameter(6,6.4e-19);
    //fPol6->FixParameter(5,0);
    //fPol6->FixParameter(6,0);
    hRelSyst->Fit(fPol6,"S, N, R, M");
    fPol6->SetLineColor(kGreen+2);
    fPol6->SetLineWidth(4);

    hRelSyst->SetStats(false);
    hRelSyst->SetTitle("");
    hRelSyst->GetXaxis()->SetTitle("#it{k*} (MeV/#it{c})");
    hRelSyst->GetXaxis()->SetTitleSize(0.06);
    hRelSyst->GetXaxis()->SetLabelSize(0.06);
    hRelSyst->GetXaxis()->CenterTitle();
    hRelSyst->GetXaxis()->SetTitleOffset(1.3);
    hRelSyst->GetXaxis()->SetLabelOffset(0.02);

    hRelSyst->GetYaxis()->SetTitle("Rel. syst. uncert.");
    hRelSyst->GetYaxis()->SetTitleSize(0.06);
    hRelSyst->GetYaxis()->SetLabelSize(0.06);
    hRelSyst->GetYaxis()->CenterTitle();
    hRelSyst->GetYaxis()->SetTitleOffset(1.20);
    //hRelSyst->GetXaxis()->SetNdivisions(506);
    hRelSyst->GetYaxis()->SetRangeUser(0.0, 0.015);
    //hRelSyst->GetYaxis()->SetRangeUser(0.90, 1.80);
    hRelSyst->SetLineWidth(4);

    TLegend* lLegend = new TLegend(0.55,0.75,0.95,0.95);//lbrt
    lLegend->SetName(TString::Format("lLegend"));
    lLegend->SetTextSize(0.045);
    lLegend->AddEntry(hRelSyst,"Syst. uncert., x4 rebinning");
    lLegend->AddEntry(fPol6,"Fit with pol6");

    TCanvas* cSyst = new TCanvas("cSyst", "cSyst", 1);
    cSyst->cd(0); cSyst->SetCanvasSize(1920/2, 1080/2); cSyst->SetMargin(0.15,0.05,0.2,0.05);//lrbt
    hRelSyst->Draw();
    fPol6->Draw("same");
    lLegend->Draw("same");
    cSyst->SaveAs(BaseName+".png");

    FileName = SystBaseName+".root";
    TFile* SystFile = new TFile(FileName,"recreate");
    hRelSyst2D->Write();
    hRelSyst->Write();
    fPol6->Write();

    delete hRelSyst2D;
    delete SystFile;
    delete hRelSyst;
    delete fPol6;
    delete lLegend;
    delete cSyst;
    delete fData0;

}

void Sigma0_To_Lambda(){

    DLM_CommonAnaFunctions AnalysisObject;
    const TString DataSample = "pp13TeV_HM_March19";
    const TString SourceType = "Gauss";
    int VARIATIONS[10];
    for(int i=0; i<10; i++) VARIATIONS[i]=0;
    const TString DataVar = "_0";
    const double SourceScale = 1.3;
    const double SourceStability = 2.0;
    const TString pL_Pot = "Usmani";

    //the residuals that are not fitted, are assumed to all have a Gaussian source of some size;
    double ResidualSourceSize=0;
    if(DataSample=="pp13TeV_MB_Run2paper") ResidualSourceSize=1.1;
    else if(DataSample=="pPb5TeV_Run2paper"||DataSample=="pPb5TeV_CPR_Mar19") ResidualSourceSize=1.4;
    else if(DataSample=="pp13TeV_HM_March19") ResidualSourceSize=1.4;
    else printf("\033[1;31mERROR:\033[0m The data sample '%s' does not exist\n",DataSample.Data());

    double* MomBins_pL = NULL;
    double* FitRegion_pL = NULL;
    unsigned NumMomBins_pL;
    AnalysisObject.SetUpBinning_pL(DataSample,NumMomBins_pL,MomBins_pL,FitRegion_pL,VARIATIONS[0],VARIATIONS[1]);

    DLM_Ck* Ck_pSigma0 = new DLM_Ck(1,0,NumMomBins_pL,MomBins_pL,Lednicky_gauss_Sigma0);
    Ck_pSigma0->SetSourcePar(0,ResidualSourceSize);

    TH2F* hResolution_pL = AnalysisObject.GetResolutionMatrix(DataSample,"pLambda");

    TH2F* hResidual_pL_pSigma0 = AnalysisObject.GetResidualMatrix("pLambda","pSigma0");
    TH2F* hResidual_pL_pXim = AnalysisObject.GetResidualMatrix("pLambda","pXim");

    TH1F* hData_pL = AnalysisObject.GetAliceExpCorrFun(DataSample,"pLambda",DataVar,2,false,VARIATIONS[5]);

    double lam_pL[5];
    double lam_pXim[5];

    AnalysisObject.SetUpLambdaPars_pL(DataSample,VARIATIONS[2],VARIATIONS[3],lam_pL);
    AnalysisObject.SetUpLambdaPars_pXim(DataSample,VARIATIONS[2],0,lam_pXim);

    TString SourceDescription = SourceType;
    if(SourceScale<0) {SourceDescription+="ScaFree";}
    else if(SourceScale==0) {SourceDescription+="ScaDefault";}
    else {SourceDescription+=TString::Format("Sca%.2f",SourceScale);}
    if(SourceStability<0||SourceStability>2) {SourceDescription+="StaFree";}
    else if(SourceStability==0) {SourceDescription+="StaDefault";}
    else {SourceDescription+=TString::Format("Sta%.2f",SourceStability);}

    CATS AB_pL;
    DLM_Ck* Ck_pL;
    if(pL_Pot.Contains("Lednicky_")){
        Ck_pL = AnalysisObject.SetUpLednicky_pL(NumMomBins_pL,MomBins_pL,pL_Pot);
        Ck_pL->SetSourcePar(0,SourceScale);
    }
    else{
        AB_pL.SetMomBins(NumMomBins_pL,MomBins_pL);
        AnalysisObject.SetUpCats_pL(AB_pL,pL_Pot,SourceType,0,3);
        AB_pL.SetNotifications(CATS::nWarning);
        AB_pL.KillTheCat();
        Ck_pL = new DLM_Ck(AB_pL.GetNumSourcePars(),0,AB_pL);
        Ck_pL->SetSourcePar(0,SourceScale);
        if(AB_pL.GetNumSourcePars()>1){
            Ck_pL->SetSourcePar(0,SourceStability);
        }
    }

    CATS AB_pXim;
    //same binning as pL, as we only use pXim as feed-down
    AB_pXim.SetMomBins(NumMomBins_pL,MomBins_pL);
    AnalysisObject.SetUpCats_pXim(AB_pXim,"pXim_HALQCD1","Gauss",0,1);
    AB_pXim.SetAnaSource(0,ResidualSourceSize);
    AB_pXim.SetNotifications(CATS::nWarning);
    AB_pXim.KillTheCat();
    DLM_Ck* Ck_pXim = new DLM_Ck(AB_pXim.GetNumSourcePars(),0,AB_pXim);
    Ck_pL->Update();
    Ck_pSigma0->Update();
    Ck_pXim->Update();

    DLM_CkDecomposition CkDec_pL("pLambda",4,*Ck_pL,hResolution_pL);
    DLM_CkDecomposition CkDec_pSigma0("pSigma0",0,*Ck_pSigma0,NULL);
    DLM_CkDecomposition CkDec_pXim("pXim",2,*Ck_pXim,NULL);

    CkDec_pL.AddContribution(0,lam_pL[1],DLM_CkDecomposition::cFeedDown,&CkDec_pSigma0,hResidual_pL_pSigma0);
    CkDec_pL.AddContribution(1,lam_pL[2],DLM_CkDecomposition::cFeedDown,&CkDec_pXim,hResidual_pL_pXim);
    CkDec_pL.AddContribution(2,lam_pL[3],DLM_CkDecomposition::cFeedDown);
    CkDec_pL.AddContribution(3,lam_pL[4],DLM_CkDecomposition::cFake);//0.03

    //for Xim we simplify a bit and take ALL feed-down as flat
    CkDec_pXim.AddContribution(0,lam_pXim[1]+lam_pXim[2]+lam_pXim[3],DLM_CkDecomposition::cFeedDown);
    CkDec_pXim.AddContribution(1,lam_pXim[4],DLM_CkDecomposition::cFake);

    CkDec_pL.Update(true);

    TFile* fOut = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/Sigma0_To_Lambda/fOut.root","recreate");

    TGraph TotalSignal;
    TGraph LambdaSignal;
    TGraph Sigma0Signal;
    TGraph XimSignal;
    TGraph FlatSignal;
    TGraph MissSignal;

    TotalSignal.SetName("TotalSignal");
    TotalSignal.Set(AB_pL.GetNumMomBins());

    LambdaSignal.SetName("LambdaSignal");
    LambdaSignal.Set(AB_pL.GetNumMomBins());

    Sigma0Signal.SetName("Sigma0Signal");
    Sigma0Signal.Set(AB_pL.GetNumMomBins());

    XimSignal.SetName("XimSignal");
    XimSignal.Set(AB_pL.GetNumMomBins());

    FlatSignal.SetName("FlatSignal");
    FlatSignal.Set(AB_pL.GetNumMomBins());

    MissSignal.SetName("MissSignal");
    MissSignal.Set(AB_pL.GetNumMomBins());

    for(unsigned uBin=0; uBin<AB_pL.GetNumMomBins(); uBin++){
        double Momentum = AB_pL.GetMomentum(uBin);
        double TotCk = (CkDec_pL.EvalCk(Momentum)-1.);
        double LamCk = (CkDec_pL.EvalSmearedMain(Momentum)-1.)*lam_pL[0];
        double SigCk = (CkDec_pL.EvalSmearedFeed(0,Momentum)-1.)*lam_pL[1];
        double XimCk = (CkDec_pL.EvalSmearedFeed(1,Momentum)-1.)*lam_pL[2];
        double FltCk = 0;
        double MisId = 0;
        TotalSignal.SetPoint(uBin,Momentum,TotCk);
        LambdaSignal.SetPoint(uBin,Momentum,LamCk);
        Sigma0Signal.SetPoint(uBin,Momentum,SigCk);
        XimSignal.SetPoint(uBin,Momentum,XimCk);
        FlatSignal.SetPoint(uBin,Momentum,FltCk);
        MissSignal.SetPoint(uBin,Momentum,MisId);
    }

    TotalSignal.Write();
    LambdaSignal.Write();
    Sigma0Signal.Write();
    XimSignal.Write();
    FlatSignal.Write();
    MissSignal.Write();

    delete Ck_pXim;

    delete Ck_pL;

    delete Ck_pSigma0;
    delete [] MomBins_pL;
    delete [] FitRegion_pL;
    ////delete hDataClever_pL;
    delete hResolution_pL;
    delete hResidual_pL_pSigma0;
    delete hResidual_pL_pXim;
    delete hData_pL;
}

void STUPED_TEST(){
    DLM_CommonAnaFunctions AnalysisObject;

    int VARIATIONS[6];
        VARIATIONS[0] = 10;
        VARIATIONS[1] = 10;
        VARIATIONS[2] = 0;
        VARIATIONS[3] = 0;
        VARIATIONS[4] = -1;//!
        VARIATIONS[5] = 0;

        Fit_pL(AnalysisObject,
               "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/STUPID_TEST/",
               "STUPUD_TEST.root",
               "pp13TeV_HM_March19",
               "_0",//datavar
               "McGauss_Reso",//source
               -11,//scale
               2,//stab
               "Usmani",
               "Longbaseline2",//BL
               VARIATIONS,
               121,
               "ntResult.root",
               "ntResult");
//DLM_CommonAnaFunctions& AnalysisObject, const TString& OutputFolder, const TString& OutputFileName,
//            const TString& DataSample, const TString& DataVar, const TString& SourceType, const double& SourceScale, const double& SourceStability,
//            const TString& pL_Pot,const TString& FittingMode_pL, const int* VARIATIONS, const long& UniqueID=-1, const TString ntFileName="", const TString ntName=""
}

int PLAMBDA_1_MAIN(int argc, char *argv[]){
printf("PLAMBDA_1_MAIN\n");

//STUPED_TEST();
//return 0;

//    Plot_mT_Scale_pp_pL("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/Plot_mT_Scale_pp_pL/",
//                        "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pp/Systematics_180919/NTfile_20.root",
//                        "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample10HM_ver2/CFOutputALL_mT_pp_HM.root",
//                        "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pL/Systematics_180919/NTfile_30.root",
//                        "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample10HM/CFOutputALL_mT_pL_HM.root",
//                        "Gaussian core", "EPOS angle");

    //Sigma0_To_Lambda();
//return 0;
//SystematicsAdd_100419_2
//const TString& OutputFolder, const int& WhichConfiguration, const int& FirstIter, const int& LastIter, const unsigned& OnlyFraction,
                      //const int& RANDOMSEED, const bool& JustNumIter
//run for 20 and 22
    //pp_SystematicsHM("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pp/Test/",
    //                 21,0,0,1,1,0);
    //pp_SystematicsHM("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pp/Systematics_180919/",
    //                 atoi(argv[1]),atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atoi(argv[5]),atoi(argv[6]));
    //pp_SystematicsHM("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pp/SystematicsAdd_300419/",
    //                 atoi(argv[1])+2,atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atoi(argv[5]),atoi(argv[6]));
    //pp_SystematicsHM("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pp/SystematicsAdd_300419/",
    //                 atoi(argv[1])+1,atoi(argv[2])*3,(atoi(argv[3])+1)*3-1,atoi(argv[4])*3,atoi(argv[5]),atoi(argv[6]));
//return 0;
//run for 30 and 32
    //pL_SystematicsHM("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pL/Systematics_020919/",
    //                 30,0,19845,200,1,0);
    pL_SystematicsHM("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pL/Systematics_1mtbin_IMPROVED_RESO/",
                     atoi(argv[1]),atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atoi(argv[5]),atoi(argv[6]));
    //pL_SystematicsHM("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pL/SystematicsAdd_300419/",
    //                 atoi(argv[1])+2,atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atoi(argv[5]),atoi(argv[6]));
    //pL_SystematicsHM("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pL/SystematicsAdd_300419/",
    //                 atoi(argv[1])+1,atoi(argv[2])*3,(atoi(argv[3])+1)*3-1,atoi(argv[4])*3,atoi(argv[5]),atoi(argv[6]));

    //pL_SystematicsHM("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pL/SystematicsAdd_070419/",
    //                 0,0,4,1,1,0);
return 0;
//Systematics_080519 in prelim.
    Plot_mT_Scale(  "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pL/Systematics_190619/PLOT/",
                    "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pL/Systematics_190619/NTfile_30.root",
                    "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample10HM/CFOutputALL_mT_pL_HM.root",
                    "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pL/Systematics_080519/PLOT/PRELIM/mTRad_pp_GaussReso.root",
                    "Gaussian core","with resonances");
    Plot_mT_Scale(  "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pL/Systematics_190619/PLOT/",
                    "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pL/Systematics_190619/NTfile_32.root",
                    "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample10HM/CFOutputALL_mT_pL_HM.root",
                    "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pL/Systematics_080519/PLOT/PRELIM/mTRad_pp_Gauss.root",
                    "Gaussian source","");

    //ComputeDataSystematics();

    return 0;

    int VARIATIONS[6];for(int i=0;i<4;i++)VARIATIONS[i]=0;
    VARIATIONS[4] = 580;//kc
    VARIATIONS[5] = -1;//mT

    DLM_CommonAnaFunctions AnalysisObject;
    const TString OutputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pL/TEMP/";
    //Fit_pL(OutputFolder,"pp13TeV_HM_March19", "Gauss", false, "NLO_Coupled_S","Longbaseline",VARIATIONS);
    //Fit_pL(OutputFolder,"pp13TeV_HM_March19", "McGauss_Reso", false, "NLO_Coupled_S","Longbaseline",VARIATIONS);
    //Fit_pL(AnalysisObject,OutputFolder,"","pp13TeV_HM_March19", "McLevyNolan_Reso", 0, 0, "NLO_Coupled_S","Longbaseline",VARIATIONS);

    //Fit_pL(OutputFolder,"pp13TeV_HM_March19", "Gauss", false, "LO_Coupled_S","Longbaseline",VARIATIONS);
    //Fit_pL(OutputFolder,"pp13TeV_HM_March19", "McGauss_Reso", false, "LO_Coupled_S","Longbaseline",VARIATIONS);
    //Fit_pL(AnalysisObject,OutputFolder,"","pp13TeV_HM_March19", "McLevyNolan_Reso", 0, 0, "LO_Coupled_S","Longbaseline",VARIATIONS);

    //Fit_pL(OutputFolder,"pp13TeV_HM_March19", "Gauss", false, "NLO_Coupled_S","Longbaseline2",VARIATIONS);
    //Fit_pL(OutputFolder,"pp13TeV_HM_March19", "McGauss_Reso", false, "NLO_Coupled_S","Longbaseline2",VARIATIONS);
    //Fit_pL(AnalysisObject,OutputFolder,"","pp13TeV_HM_March19", "McLevyNolan_Reso", 0, 0, "NLO_Coupled_S","Longbaseline2",VARIATIONS);
    //Fit_pL(OutputFolder,"pp13TeV_HM_March19", "Gauss", false, "LO_Coupled_S","Longbaseline2",VARIATIONS);
    //Fit_pL(OutputFolder,"pp13TeV_HM_March19", "McGauss_Reso", false, "LO_Coupled_S","Longbaseline2",VARIATIONS);
    //Fit_pL(AnalysisObject,OutputFolder,"","pp13TeV_HM_March19", "McLevyNolan_Reso", 0, 0, "LO_Coupled_S","Longbaseline2",VARIATIONS);

    //Fit_pL("pp13TeV_HM_March19", "Gauss", false, "NLO",VARIATIONS);
    //Fit_pL("pp13TeV_HM_March19", "McGauss_Reso", false, "NLO",VARIATIONS);
    //Fit_pL("pp13TeV_HM_March19", "McLevyNolan_Reso", false, "NLO",VARIATIONS);

    //Fit_pL("pp13TeV_HM_March19", "Gauss", false, "LO",VARIATIONS);
    //Fit_pL("pp13TeV_HM_March19", "McGauss_Reso", false, "LO",VARIATIONS);
    //Fit_pL("pp13TeV_HM_March19", "McLevyNolan_Reso", false, "LO",VARIATIONS);

    //Fit_pL(AnalysisObject,OutputFolder,"","pp13TeV_HM_March19", "_0", "McGauss_Reso",-11, 2.0,"NLO_Coupled_S","Longbaseline2", VARIATIONS);
    //Fit_pL(AnalysisObject,OutputFolder,"","pp13TeV_HM_March19", "_0", "McLevyNolan_Reso",-11, -11,"NLO_Coupled_S","Longbaseline2", VARIATIONS);

//    return 0;



    //Fit_pp_pL("pp13TeV_HM_March19","Gauss","Norm","Norm");
    /*
    Plot_pL_FAST("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pp_pL/020419/",
                 "pp13TeV_HM_March19","Gauss","Norm","Norm",
                 "LO_Coupled_S", "NLO_Coupled_S",
                 //"#it{r} = whatever", "#alpha = whatever",
                 "Gauss source fixed to p-p", "",
                 1, 0, "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample10HM/Systematics_pL.root"
                 );
    */

/*
    //Fit_pp_pL("pp13TeV_HM_March19","McGauss_Reso","Norm","Norm");
    Plot_pL_FAST("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pp_pL/250319/NLO_Coupled/",
                 "pp13TeV_HM_March19","McGauss_Reso","Norm","Norm",
                 "LO_Coupled_S", "NLO_Coupled_S",
                 "Gauss core fixed to p-p", "",
                 1, 0, "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample10HM/Systematics_pL.root"
                 );
*/
/*
    Fit_pp_pL("pp13TeV_HM_March19","McLevyNolan_Reso","Norm","Norm");
    Plot_pL_FAST("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pp_pL/020419/",
                 "pp13TeV_HM_March19","McLevyNolan_Reso","Norm","Norm",
                 "LO", "NLO",
                 "Levy core fixed to p-p", "",
                 1, 0, "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample10HM/Systematics_pL.root"
                 );
*/


/*
    Fit_pp_pL("pp13TeV_MB_Run2paper","Gauss","Norm","Norm");
    Plot_pL_FAST("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pp_pL/250319/NLO_Coupled/",
                 "pp13TeV_MB_Run2paper","Gauss","Norm","Norm",
                 "LO_Coupled_S", "NLO_Coupled_S",
                 "Gauss source fixed to p-p", "",
                 1, 0, "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample6/C2totalsysPL.root"
                 );

    Fit_pp_pL("pp13TeV_MB_Run2paper","McGauss_Reso","Norm","Norm");
    Plot_pL_FAST("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pp_pL/250319/NLO_Coupled/",
                 "pp13TeV_MB_Run2paper","McGauss_Reso","Norm","Norm",
                 "LO_Coupled_S", "NLO_Coupled_S",
                 "Gauss core fixed to p-p", "",
                 1, 0, "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample6/C2totalsysPL.root"
                 );

    Fit_pp_pL("pp13TeV_MB_Run2paper","McLevyNolan_Reso","Norm","Norm");
    Plot_pL_FAST("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pp_pL/250319/NLO_Coupled/",
                 "pp13TeV_MB_Run2paper","McLevyNolan_Reso","Norm","Norm",
                 "LO_Coupled_S", "NLO_Coupled_S",
                 "Levy core fixed to p-p", "",
                 1, 0, "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample6/C2totalsysPL.root"
                 );
*/


/*
    Fit_pp_pL("pPb5TeV_Run2paper","Gauss","Norm","Norm");
    Plot_pL_FAST("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pp_pL/250319/NLO_Coupled/",
                 "pPb5TeV_Run2paper","Gauss","Norm","Norm",
                 "LO_Coupled_S", "NLO_Coupled_S",
                 "Gauss source fixed to p-p", "",
                 1, 0, "/home/dmihaylov/CernBox/pPb/Systematics/C2totalsysPL.root"
                 );

    Fit_pp_pL("pPb5TeV_Run2paper","McGauss_Reso","Norm","Norm");
    Plot_pL_FAST("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pp_pL/250319/NLO_Coupled/",
                 "pPb5TeV_Run2paper","McGauss_Reso","Norm","Norm",
                 "LO_Coupled_S", "NLO_Coupled_S",
                 "Gauss core fixed to p-p", "",
                 1, 0, "/home/dmihaylov/CernBox/pPb/Systematics/C2totalsysPL.root"
                 );

    Fit_pp_pL("pPb5TeV_Run2paper","McLevyNolan_Reso","Norm","Norm");
    Plot_pL_FAST("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pp_pL/250319/NLO_Coupled/",
                 "pPb5TeV_Run2paper","McLevyNolan_Reso","Norm","Norm",
                 "LO_Coupled_S", "NLO_Coupled_S",
                 "Levy core fixed to p-p", "",
                 1, 0, "/home/dmihaylov/CernBox/pPb/Systematics/C2totalsysPL.root"
                 );
*/
/*
    Fit_pp_pL("pPb5TeV_CPR_Mar19","Gauss","Norm","Norm");
    Plot_pL_FAST("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pp_pL/020419/",
                 "pPb5TeV_CPR_Mar19","Gauss","Norm","Norm",
                 "LO_Coupled_S", "NLO_Coupled_S",
                 "Gauss source fixed to p-p (1.41 fm)", "",
                 1, 0, "/home/dmihaylov/CernBox/SystematicsAndCalib/pPbRun2_MB_StreamFIX/C2totalsysPL.root"
                 );
*/
/*
    Fit_pp_pL("pPb5TeV_CPR_Mar19","McGauss_Reso","Norm","Norm");
    Plot_pL_FAST("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pp_pL/020419/",
                 "pPb5TeV_CPR_Mar19","McGauss_Reso","Norm","Norm",
                 "LO_Coupled_S", "NLO_Coupled_S",
                 "Gauss core fixed to p-p", "",
                 1, 0, "/home/dmihaylov/CernBox/SystematicsAndCalib/pPbRun2_MB_StreamFIX/C2totalsysPL.root"
                 );
*/
/*
    Fit_pp_pL("pPb5TeV_CPR_Mar19","McLevyNolan_Reso","Norm","Norm");
    Plot_pL_FAST("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pp_pL/250319/NLO_Coupled/",
                 "pPb5TeV_CPR_Mar19","McLevyNolan_Reso","Norm","Norm",
                 "LO_Coupled_S", "NLO_Coupled_S",
                 "Levy core fixed to p-p", "",
                 1, 0, "/home/dmihaylov/CernBox/SystematicsAndCalib/pPbRun2_MB_StreamFIX/C2totalsysPL.root"
                 );
*/

    //Fit_pp_pL("pp13TeV_HM_March19","McLevyNolan_Reso","Longbaseline2","Longbaseline2");

    //Fit_pp_pL("pp13TeV_MB_Run2paper","Gauss","Norm","Norm");
    //Fit_pp_pL("pp13TeV_MB_Run2paper","CleverLevy_Nolan","Norm","Norm");
    //Fit_pp_pL("pp13TeV_MB_Run2paper","McLevyNolan_Reso","Norm","Norm");

    //Fit_pp_pL("pPb5TeV_Run2paper","Gauss","Norm","Norm");
    //Fit_pp_pL("pPb5TeV_Run2paper","CleverLevy_Nolan","Norm","Norm");
    //Fit_pp_pL("pPb5TeV_Run2paper","McLevyNolan_Reso","Norm","Norm");

    //Fit_pp_pL("pp13TeV_HM_March19","Gauss","Longbaseline","Longbaseline");
    //Fit_pp_pL("pp13TeV_HM_March19","CleverLevy_Nolan","Longbaseline","Longbaseline");

    //Fit_pp_pL("pp13TeV_MB_Run2paper","Gauss","Longbaseline","Longbaseline");
    //Fit_pp_pL("pp13TeV_MB_Run2paper","CleverLevy_Nolan","Longbaseline","Longbaseline");

    //Fit_pp_pL("pPb5TeV_Run2paper","Gauss","Longbaseline","Longbaseline");
    //Fit_pp_pL("pPb5TeV_Run2paper","CleverLevy_Nolan","Longbaseline","Longbaseline");


    return 0;
}
