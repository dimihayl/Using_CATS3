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

    const TString OutputFolder = TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pp_pL/040419/");

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
    const bool Use_pL_Pot[Num_pL_pots] =    {false,  true,          false,   false,     true,              false};
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

    TH1F* hData_pp = AnalysisObject.GetAliceExpCorrFun(DataSample,"pp",0);
//! this has to be changed for finer binning
    TH1F* hData_pL = AnalysisObject.GetAliceExpCorrFun(DataSample,"pLambda",2);

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
        AnalysisObject.SetUpCats_pp(AB_pp,pp_Pot[uPot_pp],SourceType);
        AB_pp.SetNotifications(CATS::nWarning);
        AB_pp.KillTheCat();
        DLM_Ck* Ck_pp = new DLM_Ck(AB_pp.GetNumSourcePars(),0,AB_pp);
        for(unsigned uPot_pL=0; uPot_pL<Num_pL_pots; uPot_pL++){
            if(!Use_pL_Pot[uPot_pL]) continue;
            CATS AB_pL;
            AB_pL.SetMomBins(NumMomBins_pL,MomBins_pL);
            AnalysisObject.SetUpCats_pL(AB_pL,pL_Pot[uPot_pL],SourceType);
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

                if(FittingMode_pp=="Norm"||FittingMode_pp=="Baseline"){
//printf("FitRegion_pp[0]=%f FitRegion_pp[1]=%f\n",FitRegion_pp[0],FitRegion_pp[1]);
                    fitter->SetSystem(0,*hData_pp,1,CkDec_pp,
                                  FitRegion_pp[0],FitRegion_pp[1],FitRegion_pp[1],FitRegion_pp[1]);
                }
                else if(FittingMode_pp=="Longbaseline"){
                    fitter->SetSystem(0,*hData_pp,1,CkDec_pp,
                                  FitRegion_pp[0],FitRegion_pp[1],FitRegion_pp[2],FitRegion_pp[3]);
                }
                else{
                    printf("\033[1;31mERROR:\033[0m Unknown fitting mode '%s'",FittingMode_pp.Data());
                }

                if(FittingMode_pL=="Norm"||FittingMode_pL=="Baseline"){
                    fitter->SetSystem(1,*hData_pL,1,CkDec_pL,
                                  FitRegion_pL[0],FitRegion_pL[1],FitRegion_pL[1],FitRegion_pL[1]);
                }
                else if(FittingMode_pL=="Longbaseline"){
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
                if(FittingMode_pp=="Baseline"||FittingMode_pp=="Longbaseline")
                    fitter->SetParameter("pp",DLM_Fitter1::p_b,0,-2e-3,2e-3);
                else
                    fitter->FixParameter("pp",DLM_Fitter1::p_b,0);
                if(FittingMode_pL=="Baseline"||FittingMode_pL=="Longbaseline")
                    fitter->SetParameter("pLambda",DLM_Fitter1::p_b,0,-2e-3,2e-3);
                else
                    fitter->FixParameter("pLambda",DLM_Fitter1::p_b,0);

                fitter->FixParameter("pp",DLM_Fitter1::p_c,0);
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
                printf("R(pΛ) = %.3f +/- %.3f\n",fitter->GetParameter("pLambda",DLM_Fitter1::p_sor0),fitter->GetParError("pLambda",DLM_Fitter1::p_sor0));
                if(SourceType=="CleverLevy_Nolan"||SourceType=="CleverLevy_Single"||SourceType=="CleverLevy_Diff"||
                   SourceType=="McLevyNolan_Reso"||SourceType=="McGauss_Reso"){
                    printf("α(pp) = %.3f +/- %.3f\n",fitter->GetParameter("pp",DLM_Fitter1::p_sor1),fitter->GetParError("pp",DLM_Fitter1::p_sor1));
                    printf("α(pΛ) = %.3f +/- %.3f\n",fitter->GetParameter("pLambda",DLM_Fitter1::p_sor1),fitter->GetParError("pLambda",DLM_Fitter1::p_sor1));
                }
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
void Fit_pL(DLM_CommonAnaFunctions& AnalysisObject, const TString& OutputFolder, const TString& OutputFileName,
            const TString& DataSample, const TString& SourceType, const double& SourceScale, const double& SourceStability,
            const TString& pL_Pot,const TString& FittingMode_pL, const int* VARIATIONS, const int& UniqueID=-1, const TString ntFileName="", const TString ntName=""){

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

    TH1F* hData_pL = AnalysisObject.GetAliceExpCorrFun(DataSample,"pLambda",2,true,VARIATIONS[5]);

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
        OutFileName=TString::Format("AUTO%i_%s_%s_%s_%s_%i_%i_%i_%i_%i_mT%i.root",UniqueID,
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
        AnalysisObject.SetUpCats_pL(AB_pL,pL_Pot,SourceType);
        AB_pL.SetNotifications(CATS::nWarning);
        AB_pL.KillTheCat();
        DLM_Ck* Ck_pL = new DLM_Ck(AB_pL.GetNumSourcePars(),0,AB_pL);
        Ck_pL->SetSourcePar(0,SourceScale);
        if(AB_pL.GetNumSourcePars()>1){
            Ck_pL->SetSourcePar(0,SourceStability);
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
    else if(FittingMode_pL.Contains("Longbaseline")){
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



    //!p_sor0
    if(SourceType=="Gauss"){
        if(SourceScale==0){fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,1.2);}
        else if(SourceScale<0) {fitter->SetParameter("pLambda",DLM_Fitter1::p_sor0,1.2,0.8,1.8);}
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
            else if(SourceStability<1||SourceStability>2) {fitter->SetParameter("pLambda",DLM_Fitter1::p_sor1,1.5,1.5*0.90,1.5*1.10);}
            else {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor1,SourceStability);}
        }
        else if(DataSample=="pPb5TeV_Run2paper"||DataSample=="pPb5TeV_CPR_Mar19"){
            if(SourceScale==0) {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,1.02);}//at mT 1.44
            else if(SourceScale==-1) {fitter->SetParameter("pLambda",DLM_Fitter1::p_sor0,1.02,1.02*0.90,1.02*1.10);}
            else if(SourceScale==-11) {fitter->SetParameter("pLambda",DLM_Fitter1::p_sor0,1.02,1.02/2.,1.02*2.);}
            else {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,SourceScale);}

            if(SourceStability==0) {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor1,1.45);}
            else if(SourceStability<1||SourceStability>2) {fitter->SetParameter("pLambda",DLM_Fitter1::p_sor1,1.45,1.45*0.90,1.45*1.10);}
            else {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor1,SourceStability);}
        }
        else if(DataSample=="pp13TeV_HM_March19"){
            if(SourceScale==0) {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,0.92);}//at mT 1.55
            else if(SourceScale==-1) {fitter->SetParameter("pLambda",DLM_Fitter1::p_sor0,0.92,0.92*0.90,0.92*1.10);}
            else if(SourceScale==-11) {fitter->SetParameter("pLambda",DLM_Fitter1::p_sor0,0.92,0.92/2.,0.92*2.);}
            else {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,SourceScale);}

            if(SourceStability==0) {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor1,1.5);}
            else if(SourceStability<1||SourceStability>2) {fitter->SetParameter("pLambda",DLM_Fitter1::p_sor1,1.5,1.5*0.90,1.5*1.10);}
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
            if(SourceScale==0) {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,0.87);}//at mT 1.55
            else if(SourceScale==-1) {fitter->SetParameter("pLambda",DLM_Fitter1::p_sor0,0.87,0.87*0.90,0.87*1.10);}
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

    printf("χ2/ndf = %.2f/%i = %.2f\n",fitter->GetChi2(),fitter->GetNdf(),fitter->GetChi2Ndf());
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
            "Iter:Config:SourceType:SourceScale:SourceStability:Potential:Baseline:FemRan:FitRan:pFrac:LamFrac:kcVar:mTbin:FemtoMin:FemtoMax:BlMin:BlMax:"
            "p_a:e_a:p_b:e_b:p_c:e_c:p_Cl:e_Cl:p_kc:e_kc:p_sor0:e_sor0:p_sor1:e_sor1:chi2:ndf");
        }
    }

    if(NtFile&&ntResult){
        Float_t buffer[33];
        buffer[0] = UniqueID/1000;
        buffer[1] = UniqueID%1000;
        if(SourceType=="Gauss") buffer[2] = 0;
        if(SourceType=="McGauss_Reso") buffer[2] = 1;
        if(SourceType=="McLevyNolan_Reso") buffer[2] = 2;

        buffer[3] = SourceScale;
        buffer[4] = SourceStability;

        if(pL_Pot=="LO") buffer[5] = 0;
        if(pL_Pot=="NLO") buffer[5] = 10;
        if(pL_Pot=="LO_Coupled_S") buffer[5] = 1;
        if(pL_Pot=="NLO_Coupled_S") buffer[5] = 11;
        if(pL_Pot=="NLO_sp") buffer[5] = 12;
        if(pL_Pot=="Usmani") buffer[5] = 100;
        if(pL_Pot=="Lednicky_ND") buffer[5] = 1001;
        if(pL_Pot=="Lednicky_NF") buffer[5] = 1002;
        if(pL_Pot=="Lednicky_NSC89") buffer[5] = 1003;
        if(pL_Pot=="Lednicky_NSC97a") buffer[5] = 1004;
        if(pL_Pot=="Lednicky_NSC97b") buffer[5] = 1005;
        if(pL_Pot=="Lednicky_NSC97c") buffer[5] = 1006;
        if(pL_Pot=="Lednicky_NSC97d") buffer[5] = 1007;
        if(pL_Pot=="Lednicky_NSC97e") buffer[5] = 1008;
        if(pL_Pot=="Lednicky_NSC97f") buffer[5] = 1009;
        if(pL_Pot=="Lednicky_ESC08") buffer[5] = 1010;
        if(pL_Pot=="Lednicky_XeftLO") buffer[5] = 1011;
        if(pL_Pot=="Lednicky_XeftNLO") buffer[5] = 1012;
        if(pL_Pot=="Lednicky_JulichA") buffer[5] = 1013;
        if(pL_Pot=="Lednicky_JulichJ04") buffer[5] = 1014;
        if(pL_Pot=="Lednicky_JulichJ04c") buffer[5] = 1015;

        if(FittingMode_pL=="Norm") buffer[6] = 0;
        if(FittingMode_pL=="Baseline") buffer[6] = 1;
        if(FittingMode_pL=="Baseline2") buffer[6] = 2;
        if(FittingMode_pL=="Longbaseline") buffer[6] = 11;
        if(FittingMode_pL=="Longbaseline2") buffer[6] = 12;
        if(FittingMode_pL=="Norm_prefit") buffer[6] = 100;
        if(FittingMode_pL=="Baseline_prefit") buffer[6] = 101;
        if(FittingMode_pL=="Baseline2_prefit") buffer[6] = 102;
        if(FittingMode_pL=="Longbaseline_prefit") buffer[6] = 111;
        if(FittingMode_pL=="Longbaseline2_prefit") buffer[6] = 112;

        buffer[7] = VARIATIONS[0];
        buffer[8] = VARIATIONS[1];
        buffer[9] = VARIATIONS[2];
        buffer[10] = VARIATIONS[3];
        buffer[11] = VARIATIONS[4];
        buffer[12] = VARIATIONS[5];

        buffer[13] = FitRegion_pL[0];
        buffer[14] = FitRegion_pL[1];
        buffer[15] = FitRegion_pL[2];
        buffer[16] = FitRegion_pL[3];

        buffer[17] = fitter->GetParameter("pLambda",DLM_Fitter1::p_a);
        buffer[18] = fitter->GetParError("pLambda",DLM_Fitter1::p_a);

        buffer[19] = fitter->GetParameter("pLambda",DLM_Fitter1::p_b);
        buffer[20] = fitter->GetParError("pLambda",DLM_Fitter1::p_b);

        buffer[21] = fitter->GetParameter("pLambda",DLM_Fitter1::p_c);
        buffer[22] = fitter->GetParError("pLambda",DLM_Fitter1::p_c);

        buffer[23] = fitter->GetParameter("pLambda",DLM_Fitter1::p_Cl);
        buffer[24] = fitter->GetParError("pLambda",DLM_Fitter1::p_Cl);

        buffer[25] = fitter->GetParameter("pLambda",DLM_Fitter1::p_kc);
        buffer[26] = fitter->GetParError("pLambda",DLM_Fitter1::p_kc);

        buffer[27] = fitter->GetParameter("pLambda",DLM_Fitter1::p_sor0);
        buffer[28] = fitter->GetParError("pLambda",DLM_Fitter1::p_sor0);

        buffer[29] = fitter->GetParameter("pLambda",DLM_Fitter1::p_sor1);
        buffer[30] = fitter->GetParError("pLambda",DLM_Fitter1::p_sor1);

        buffer[31] = fitter->GetChi2();
        buffer[32] = fitter->GetNdf();

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
    int VARIATIONS[6];
    int TotNumIter=0;

    switch(WhichConfiguration){
        //this is the case where we only check for the Gauss+Reso and Gauss
        case 0 :
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
            NumMtVars = 6;
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
        SourceScale[0] = -1;//free fit, but within the 15% around 0.87
        SourceScale[1] = 0.87*0.90;
        SourceScale[2] = 0.87;
        SourceScale[3] = 0.87*1.10;
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

    }
    else if(WhichConfiguration==20){
        Source[0] = "McGauss_Reso";
        DefSource = 0;
        SourceScale[0] = -11;//free fit, but within the 15% around 0.87
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
        MtVars[1] = 0;
        MtVars[2] = 1;
        MtVars[3] = 2;
        MtVars[4] = 3;
        MtVars[5] = 4;
        DefMt = 0;
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
                    ipFrac!=DefProtonFrac||iLamFrac!=DefLambdaFrac||ikc!=DefKc){DefaultVar=false;}
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

        Fit_pL(AnalysisObject, OutputFolder, FileName,"pp13TeV_HM_March19",Source[iSource],SourceScale[iSourceScal],SourceStability[iSourceStab],
               Potential[iPot],Baseline[iBl],VARIATIONS,IterationCounter*1000+WhichConfiguration,ntFileName,"ntResult");

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

void Plot_mT_Scale(const TString OutputFolder, const TString ntFileName, const TString mtFileName, const TString pp_FileName,
                   const TString Title){

    TFile* mtFile = new TFile(mtFileName,"read");
    TGraphErrors* AveragemT = (TGraphErrors*)mtFile->Get("AveragemT");
    unsigned Num_mT_bins = AveragemT->GetN();
    TH2F* pL_LO_mT_sor0 = new TH2F("pL_LO_mT_sor0","pL_LO_mT_sor0",Num_mT_bins,-0.5,double(Num_mT_bins)-0.5,512,0.2,1.8);
    TH2F* pL_LO_mT_esor0 = new TH2F("pL_LO_mT_esor0","pL_LO_mT_esor0",Num_mT_bins,-0.5,double(Num_mT_bins)-0.5,512,0.0,0.4);

    TH2F* pL_NLO_mT_sor0 = new TH2F("pL_NLO_mT_sor0","pL_NLO_mT_sor0",Num_mT_bins,-0.5,double(Num_mT_bins)-0.5,512,0.2,1.8);
    TH2F* pL_NLO_mT_esor0 = new TH2F("pL_NLO_mT_esor0","pL_NLO_mT_esor0",Num_mT_bins,-0.5,double(Num_mT_bins)-0.5,512,0.0,0.4);

    TFile* ntFile = new TFile(ntFileName,"read");
    TNtuple* ntResult = (TNtuple*)ntFile->Get("ntResult");
    unsigned NumNtEntries = ntResult->GetEntries();

    Float_t Iter;
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

    ntResult->SetBranchAddress("Iter",&Iter);
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
printf("NumNtEntries=%u\n",NumNtEntries);
    for(unsigned uEntry=0; uEntry<NumNtEntries; uEntry++){
        ntFile->cd();
        ntResult->GetEntry(uEntry);
        int Config;
        Config = int(Iter)%100;
        int ITER = int(Iter)/100;
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
        if(Config!=20) continue; //select Gauss+Reso
        //if(Config!=1) continue; //select Levy+Reso
        //if(Config!=2) continue; //select Gauss
        //if(round(Baseline)!=0) continue; //select only norm
        //if(int(Baseline)!=11) continue; //select only pol1
        //if(int(Baseline)!=12) continue; //select only pol2
        if(int(Baseline)==11) continue; // select pol0 and pol2
        //if(SourceScale>0) continue;
        //if(SourceScale<0.8&&SourceScale>0) continue;

        if(mTbin<0) continue;

        if(Potential==1){
            pL_LO_mT_sor0->Fill(mTbin,p_sor0);
            pL_LO_mT_esor0->Fill(mTbin,e_sor0);
        }
        if(Potential==11){
            pL_NLO_mT_sor0->Fill(mTbin,p_sor0);
            pL_NLO_mT_esor0->Fill(mTbin,e_sor0);
        }
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

    TGraphErrors* mTSyst_pp=NULL;
    TGraphErrors* mTStat_pp=NULL;

    TFile* ppFile = NULL;
    if(pp_FileName!=""){
        ppFile = new TFile(pp_FileName,"read");
        mTSyst_pp = (TGraphErrors*)ppFile->Get("mTRadiusSyst");
        mTStat_pp = (TGraphErrors*)ppFile->Get("mTRadiusStat");
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
        TH1D* esorLO = pL_LO_mT_esor0->ProjectionY("esorLO",umT+1,umT+1);
        TH1D* esorNLO = pL_NLO_mT_esor0->ProjectionY("esorNLO",umT+1,umT+1);
        double StatErrLO = esorLO->GetMean();
        double SystErrLO = sorLO->GetStdDev();
        double StatErrNLO = esorNLO->GetMean();
        double SystErrNLO = sorNLO->GetStdDev();

        AveragemT->GetPoint(umT,Dummy,AvgmT);
        AvgmTerr = AveragemT->GetErrorY(umT);

        mTSyst_pL_LO.SetPoint(umT,AvgmT,sorLO->GetMean());//
        mTSyst_pL_LO.SetPointError(umT,AvgmTerr,SystErrLO);

        mTSyst_pL_NLO.SetPoint(umT,AvgmT,sorNLO->GetMean());//
        mTSyst_pL_NLO.SetPointError(umT,AvgmTerr,SystErrNLO);
//printf("sorLO->GetMean()=%f; sorNLO->GetMean()=%f\n",sorLO->GetMean(),sorNLO->GetMean());
        mTStat_pL_LO.SetPoint(umT,AvgmT,sorLO->GetMean());//
        mTStat_pL_LO.SetPointError(umT,0,StatErrLO);
        mTStat_pL_NLO.SetPoint(umT,AvgmT,sorNLO->GetMean());//
        mTStat_pL_NLO.SetPointError(umT,0,StatErrNLO);
    }

    DLM_Histo<double> DimiStat_pp;
    DLM_Histo<double> DimiStat_pL_LO;
    DLM_Histo<double> DimiStat_pL_NLO;

    DLM_Histo<double> DimiStatEff_pL_NLO;
    DLM_Histo<double> DimiStatRatio_pL_NLO;

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
        mTStat_pL_LO.GetPoint(0,xVal,yVal);
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

    TFile* fOut = new TFile(OutputFolder+"fOut.root","recreate");

    TH1F* hAxis = new TH1F("hAxis", "hAxis", 128, 0.9, 2.6);
    hAxis->SetStats(false);
    hAxis->SetTitle("");
    hAxis->GetXaxis()->SetTitle("<m_{T}> (GeV/#it{c}^{2})");
    hAxis->GetXaxis()->SetTitleSize(0.06);
    hAxis->GetXaxis()->SetLabelSize(0.06);
    hAxis->GetXaxis()->CenterTitle();
    hAxis->GetXaxis()->SetTitleOffset(1.3);
    hAxis->GetXaxis()->SetLabelOffset(0.02);

    //if(LevySource) hAxis->GetYaxis()->SetTitle("Levy core (fm)");
    //else hAxis->GetYaxis()->SetTitle("Gaussian core (fm)");
    hAxis->GetYaxis()->SetTitle("Gaussian core (fm)");
    hAxis->GetYaxis()->SetTitleSize(0.06);
    hAxis->GetYaxis()->SetLabelSize(0.06);
    hAxis->GetYaxis()->CenterTitle();
    hAxis->GetYaxis()->SetTitleOffset(0.90);

    //hAxis->GetXaxis()->SetNdivisions(506);
    hAxis->GetYaxis()->SetRangeUser(0.55, 1.35);

    TLegend* lLegend = new TLegend(0.65,0.65,0.95,0.95);//lbrt
    lLegend->SetName(TString::Format("lLegend"));
    lLegend->SetTextSize(0.045);
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

    fOut->cd();
    mTSyst_pL_LO.Write();
    mTSyst_pL_NLO.Write();
    mTStat_pL_LO.Write();
    mTStat_pL_NLO.Write();
    if(mTSyst_pp) mTSyst_pp->Write();
    if(mTStat_pp) mTStat_pp->Write();
    gStatRatio_pL_NLO.Write();
    cmT->Write();
    cmT->SaveAs(OutputFolder+"cmT.pdf");


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
    delete hAxis;
    delete lLegend;
    delete cmT;
    delete ntFile;
    delete mtFile;
    delete fOut;
    if(ppFile){delete ppFile; ppFile=NULL;}

}

int PLAMBDA_1_MAIN(int argc, char *argv[]){

//SystematicsAdd_100419_2
//const TString& OutputFolder, const int& WhichConfiguration, const int& FirstIter, const int& LastIter, const unsigned& OnlyFraction,
                      //const int& RANDOMSEED, const bool& JustNumIter
    pL_SystematicsHM("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pL/SystematicsAdd_250419/",
                     atoi(argv[1]),atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atoi(argv[5]),atoi(argv[6]));
    //pL_SystematicsHM("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pL/SystematicsAdd_250419/",
    //                 atoi(argv[1])+2,atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atoi(argv[5]),atoi(argv[6]));
    //pL_SystematicsHM("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pL/SystematicsAdd_250419/",
    //                 atoi(argv[1])+1,atoi(argv[2])*3,(atoi(argv[3])+1)*3-1,atoi(argv[4])*3,atoi(argv[5]),atoi(argv[6]));

    //pL_SystematicsHM("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pL/SystematicsAdd_070419/",
    //                 0,0,4,1,1,0);


    //Plot_mT_Scale(  "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Plot_mT_Scale/",
    //                "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pL/SystematicsAdd_150419/NTfile_20.root",
    //                "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample10HM/CFOutputALL_mT_pL_HM.root",
    //                "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Plot_mT_Scale/mTRad_pp.root",
    //                "mT scaling for p#minusp and p#minus#Lambda");


    return 0;

    int VARIATIONS[4];for(int i=0;i<4;i++)VARIATIONS[i]=0;

    DLM_CommonAnaFunctions AnalysisObject;
    const TString OutputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pL/050419/";
    //Fit_pL(OutputFolder,"pp13TeV_HM_March19", "Gauss", false, "NLO_Coupled_S","Longbaseline",VARIATIONS);
    //Fit_pL(OutputFolder,"pp13TeV_HM_March19", "McGauss_Reso", false, "NLO_Coupled_S","Longbaseline",VARIATIONS);
    Fit_pL(AnalysisObject,OutputFolder,"","pp13TeV_HM_March19", "McLevyNolan_Reso", 0, 0, "NLO_Coupled_S","Longbaseline",VARIATIONS);

    //Fit_pL(OutputFolder,"pp13TeV_HM_March19", "Gauss", false, "LO_Coupled_S","Longbaseline",VARIATIONS);
    //Fit_pL(OutputFolder,"pp13TeV_HM_March19", "McGauss_Reso", false, "LO_Coupled_S","Longbaseline",VARIATIONS);
    Fit_pL(AnalysisObject,OutputFolder,"","pp13TeV_HM_March19", "McLevyNolan_Reso", 0, 0, "LO_Coupled_S","Longbaseline",VARIATIONS);

    //Fit_pL(OutputFolder,"pp13TeV_HM_March19", "Gauss", false, "NLO_Coupled_S","Longbaseline2",VARIATIONS);
    //Fit_pL(OutputFolder,"pp13TeV_HM_March19", "McGauss_Reso", false, "NLO_Coupled_S","Longbaseline2",VARIATIONS);
    Fit_pL(AnalysisObject,OutputFolder,"","pp13TeV_HM_March19", "McLevyNolan_Reso", 0, 0, "NLO_Coupled_S","Longbaseline2",VARIATIONS);
    //Fit_pL(OutputFolder,"pp13TeV_HM_March19", "Gauss", false, "LO_Coupled_S","Longbaseline2",VARIATIONS);
    //Fit_pL(OutputFolder,"pp13TeV_HM_March19", "McGauss_Reso", false, "LO_Coupled_S","Longbaseline2",VARIATIONS);
    Fit_pL(AnalysisObject,OutputFolder,"","pp13TeV_HM_March19", "McLevyNolan_Reso", 0, 0, "LO_Coupled_S","Longbaseline2",VARIATIONS);

    //Fit_pL("pp13TeV_HM_March19", "Gauss", false, "NLO",VARIATIONS);
    //Fit_pL("pp13TeV_HM_March19", "McGauss_Reso", false, "NLO",VARIATIONS);
    //Fit_pL("pp13TeV_HM_March19", "McLevyNolan_Reso", false, "NLO",VARIATIONS);

    //Fit_pL("pp13TeV_HM_March19", "Gauss", false, "LO",VARIATIONS);
    //Fit_pL("pp13TeV_HM_March19", "McGauss_Reso", false, "LO",VARIATIONS);
    //Fit_pL("pp13TeV_HM_March19", "McLevyNolan_Reso", false, "LO",VARIATIONS);

    return 0;



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
