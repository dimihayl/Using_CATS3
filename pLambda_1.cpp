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

#include "CATS.h"
#include "DLM_CkDecomposition.h"
#include "DLM_CkModels.h"
#include "CommonAnaFunctions.h"
#include "DLM_Fitters.h"

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

    TFile* OutputFile = new TFile(OutputFolder+TString::Format("Output_%s_%s_%s_%s.root",DataSample.Data(),SourceType.Data(),FittingMode_pp.Data(),FittingMode_pL.Data()),"recreate");
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
//pL_Pot: "LO",  "LO_Coupled_S", "NLO", "NLO_sp",  "NLO_Coupled_S",    "Usmani"
//in VARIATIONS:
//[0] = femto range
//[1] = fit range
//[2] = fractions for p
//[3] = fractions for L
//[4] = kc parameter
void Fit_pL(DLM_CommonAnaFunctions& AnalysisObject, const TString& OutputFolder,
            const TString& DataSample, const TString& SourceType, const double& SourceScale, const double& SourceStability,
            const TString& pL_Pot,const TString& FittingMode_pL, const int* VARIATIONS, TNtuple* ntResult){

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
    if(pL_Pot=="LO"||pL_Pot=="NLO"||pL_Pot=="Usmani"){
        FitRegion_pL[1] = 192;
        FitRegion_pL[2] = 312;
        FitRegion_pL[3] = 372;
    }

    DLM_Ck* Ck_pSigma0 = new DLM_Ck(1,0,NumMomBins_pL,MomBins_pL,Lednicky_gauss_Sigma0);
    Ck_pSigma0->SetSourcePar(0,ResidualSourceSize);

    TH2F* hResolution_pL = AnalysisObject.GetResolutionMatrix(DataSample,"pLambda");

    TH2F* hResidual_pL_pSigma0 = AnalysisObject.GetResidualMatrix("pLambda","pSigma0");
    TH2F* hResidual_pL_pXim = AnalysisObject.GetResidualMatrix("pLambda","pXim");

    TH1F* hData_pL = AnalysisObject.GetAliceExpCorrFun(DataSample,"pLambda",2);

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
    TFile* OutputFile = new TFile(OutputFolder+TString::Format("Output_%s_%s_%s_%s_%i_%i_%i_%i_%i.root",
                                                               DataSample.Data(),SourceDescription.Data(),
                                                               pL_Pot.Data(),FittingMode_pL.Data(),
                                                               VARIATIONS[0],VARIATIONS[1],VARIATIONS[2],VARIATIONS[3],VARIATIONS[4]),"recreate");
    hData_pL->Write();

    CATS AB_pL;
    AB_pL.SetMomBins(NumMomBins_pL,MomBins_pL);
    AnalysisObject.SetUpCats_pL(AB_pL,pL_Pot,SourceType);
    AB_pL.SetNotifications(CATS::nWarning);
    AB_pL.KillTheCat();

    DLM_Ck* Ck_pL = new DLM_Ck(AB_pL.GetNumSourcePars(),0,AB_pL);

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

    if(FittingMode_pL.Contains("_prefit")) {fitter->SetSeparateBL(0,true);}
    else {fitter->SetSeparateBL(0,false);}


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

    //fitter->FixParameter("pLambda",DLM_Fitter1::p_Cl,-0.995);
    //fitter->SetParameter("pLambda",DLM_Fitter1::p_Cl,-1,-1,-0.95);

    //fitter->SetParameter("pLambda",DLM_Fitter1::p_Cl,1,0.9,1.1);
    if(VARIATIONS[4]>0&&FittingMode_pL!="Norm") fitter->FixParameter("pLambda",DLM_Fitter1::p_kc,VARIATIONS[4]);
    else fitter->SetParameter("pLambda",DLM_Fitter1::p_kc,2000,400,10000);
    //fitter->SetParameter("pLambda",DLM_Fitter1::p_kc,2000,400,10000);

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
            else if(SourceScale<0) {fitter->SetParameter("pLambda",DLM_Fitter1::p_sor0,0.71,0.71*0.85,0.71*1.15);}
            else {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,SourceScale);}

            if(SourceStability==0) {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor1,1.5);}
            else if(SourceStability<1||SourceStability>2) {fitter->SetParameter("pLambda",DLM_Fitter1::p_sor1,1.5,1.5*0.85,1.5*1.15);}
            else {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor1,SourceStability);}
        }
        else if(DataSample=="pPb5TeV_Run2paper"||DataSample=="pPb5TeV_CPR_Mar19"){
            if(SourceScale==0) {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,1.02);}//at mT 1.44
            else if(SourceScale<0) {fitter->SetParameter("pLambda",DLM_Fitter1::p_sor0,1.02,1.02*0.85,1.02*1.15);}
            else {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,SourceScale);}

            if(SourceStability==0) {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor1,1.45);}
            else if(SourceStability<1||SourceStability>2) {fitter->SetParameter("pLambda",DLM_Fitter1::p_sor1,1.45,1.45*0.85,1.45*1.15);}
            else {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor1,SourceStability);}
        }
        else if(DataSample=="pp13TeV_HM_March19"){
            if(SourceScale==0) {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,0.92);}//at mT 1.55
            else if(SourceScale<0) {fitter->SetParameter("pLambda",DLM_Fitter1::p_sor0,0.92,0.92*0.85,0.92*1.15);}
            else {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,SourceScale);}

            if(SourceStability==0) {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor1,1.5);}
            else if(SourceStability<1||SourceStability>2) {fitter->SetParameter("pLambda",DLM_Fitter1::p_sor1,1.5,1.5*0.85,1.5*1.15);}
            else {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor1,SourceStability);}
        }
        else printf("WHAT HAPPEND?\n");
    }
    else if(SourceType=="McGauss_Reso"){
        if(DataSample=="pp13TeV_MB_Run2paper"){
            if(SourceScale==0) {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,0.71);}//at mT 1.44
            else if(SourceScale<0) {fitter->SetParameter("pLambda",DLM_Fitter1::p_sor0,0.71,0.71*0.85,0.71*1.15);}
            else {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,SourceScale);}
            fitter->FixParameter("pLambda",DLM_Fitter1::p_sor1,2.0);
        }
        else if(DataSample=="pPb5TeV_Run2paper"||DataSample=="pPb5TeV_CPR_Mar19"){
            if(SourceScale==0) {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,0.95);}//at mT 1.52
            else if(SourceScale<0) {fitter->SetParameter("pLambda",DLM_Fitter1::p_sor0,0.95,0.95*0.85,0.95*1.15);}
            else {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,SourceScale);}
            fitter->FixParameter("pLambda",DLM_Fitter1::p_sor1,2.0);
        }
        else if(DataSample=="pp13TeV_HM_March19"){
            if(SourceScale==0) {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,0.87);}//at mT 1.55
            else if(SourceScale<0) {fitter->SetParameter("pLambda",DLM_Fitter1::p_sor0,0.87,0.87*0.85,0.87*1.15);}
            else {fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,SourceScale);}
            fitter->FixParameter("pLambda",DLM_Fitter1::p_sor1,2.0);
        }
        else printf("WHAT HAPPEND?\n");
    }
    else{
        printf("\033[1;31mERROR:\033[0m '%s' does not exist",SourceType.Data());
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



    if(ntResult){
        Float_t buffer[21];
        if(SourceType=="Gauss") buffer[0] = 0;
        if(SourceType=="McGauss_Reso") buffer[0] = 1;
        if(SourceType=="McLevyNolan_Reso") buffer[0] = 2;

        buffer[1] = SourceScale;
        buffer[2] = SourceStability;

        if(pL_Pot=="LO") buffer[3] = 0;
        if(pL_Pot=="NLO") buffer[3] = 10;
        if(pL_Pot=="LO_Coupled_S") buffer[3] = 1;
        if(pL_Pot=="NLO_Coupled_S") buffer[3] = 11;
        if(pL_Pot=="NLO_sp") buffer[3] = 12;
        if(pL_Pot=="Usmani") buffer[3] = 100;

        if(FittingMode_pL=="Norm") buffer[4] = 0;
        if(FittingMode_pL=="Baseline") buffer[4] = 1;
        if(FittingMode_pL=="Baseline2") buffer[4] = 2;
        if(FittingMode_pL=="Longbaseline") buffer[4] = 11;
        if(FittingMode_pL=="Longbaseline2") buffer[4] = 12;
        if(FittingMode_pL=="Norm_prefit") buffer[4] = 100;
        if(FittingMode_pL=="Baseline_prefit") buffer[4] = 101;
        if(FittingMode_pL=="Baseline2_prefit") buffer[4] = 102;
        if(FittingMode_pL=="Longbaseline_prefit") buffer[4] = 111;
        if(FittingMode_pL=="Longbaseline2_prefit") buffer[4] = 112;

        buffer[5] = fitter->GetParameter("pLambda",DLM_Fitter1::p_a);
        buffer[6] = fitter->GetParError("pLambda",DLM_Fitter1::p_a);

        buffer[7] = fitter->GetParameter("pLambda",DLM_Fitter1::p_b);
        buffer[8] = fitter->GetParError("pLambda",DLM_Fitter1::p_b);

        buffer[9] = fitter->GetParameter("pLambda",DLM_Fitter1::p_c);
        buffer[10] = fitter->GetParError("pLambda",DLM_Fitter1::p_c);

        buffer[11] = fitter->GetParameter("pLambda",DLM_Fitter1::p_Cl);
        buffer[12] = fitter->GetParError("pLambda",DLM_Fitter1::p_Cl);

        buffer[13] = fitter->GetParameter("pLambda",DLM_Fitter1::p_kc);
        buffer[14] = fitter->GetParError("pLambda",DLM_Fitter1::p_kc);

        buffer[15] = fitter->GetParameter("pLambda",DLM_Fitter1::p_sor0);
        buffer[16] = fitter->GetParError("pLambda",DLM_Fitter1::p_sor0);

        buffer[17] = fitter->GetParameter("pLambda",DLM_Fitter1::p_sor1);
        buffer[18] = fitter->GetParError("pLambda",DLM_Fitter1::p_sor1);

        buffer[19] = fitter->GetChi2();
        buffer[20] = fitter->GetNdf();

        ntResult->Fill(buffer);
    }




    delete Ck_pXim;
    delete fitter;

    delete Ck_pL;

    delete Ck_pSigma0;
    delete [] MomBins_pL;
    delete [] FitRegion_pL;
    //delete hDataClever_pL;
    delete OutputFile;
}


/*
void Plot_pp_pL(const TString& DataSample, const TString& SourceType, const TString& FittingMode_pp, const TString& FittingMode_pL){
    TFile* InputFile = new TFile(OutputFolder+TString::Format("Output_%s_%s_%s.root",DataSample.Data(),SourceType.Data(),FittingMode_pp.Data()),"read");




    delete InputFile;
}
*/

//JustNumIter==true => we only print the number of iter to be performed
void pL_SystematicsHM(const TString& OutputFolder, const int& WhichConfiguration, const int& FirstIter, const int& LastIter, const bool& JustNumIter=false){
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
    int VARIATIONS[5];
    int TotNumIter=0;

    switch(WhichConfiguration){
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

    TFile* NTfile = NULL;
    TNtuple* ntResult = NULL;

    if(WhichConfiguration==0){
        Source[0] = "McGauss_Reso";
        SourceScale[0] = -1;//free fit, but within the 15% around 0.92
        SourceScale[1] = 0.92*0.85;
        SourceScale[2] = 0.92*1.15;
        SourceStability[0] = 0;//dummy
        Potential[0] = "LO_Coupled_S";
        Potential[1] = "NLO_Coupled_S";
        Baseline[0] = "Norm";
        Baseline[1] = "Longbaseline";
        Baseline[2] = "Longbaseline2";

        FemtoRangeVars[0] = 0;//default
        FemtoRangeVars[1] = 1;//down
        FemtoRangeVars[2] = 2;//up

        FitRangeVars[0] = 0;//default
        FitRangeVars[1] = 1;//down
        FitRangeVars[2] = 2;//up

        pFracVars[0] = 0;//default
        pFracVars[1] = 1;
        pFracVars[2] = 2;

        LamFracVars[0] = 0;//default
        LamFracVars[1] = 1;
        LamFracVars[2] = 2;

        LamFracVars[3] = 10;
        LamFracVars[4] = 11;
        LamFracVars[5] = 12;

        LamFracVars[6] = 20;
        LamFracVars[7] = 21;
        LamFracVars[8] = 22;

        kcVars[0] = 700;
        kcVars[1] = 400;
        kcVars[2] = 1000;

        if(!JustNumIter){
            NTfile = new TFile(TString::Format("%sNTfile_%i_%i_%i.root",OutputFolder.Data(),WhichConfiguration,FirstIter,LastIter),"recreate");
            ntResult = new TNtuple("ntResult", "ntResult",
                "SourceType:SourceScale:SourceStability:Potential:Baseline:"
                "p_a:e_a:p_b:e_b:p_c:e_c:p_Cl:e_Cl:p_kc:e_kc:p_sor0:e_sor0:p_sor1:e_sor1:chi2:ndf");
        }
    }
    else{
        printf("????????????\n");
    }

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
        if(Baseline[iBl]=="Norm"&&ikc) {continue;}//this variation is only relevant in case of a baseline
        if(Baseline[iBl]=="Norm"&&iFitRan) {continue;}//this variation is only relevant in case of a baseline
        if(JustNumIter){IterationCounter++; continue;}
        if(IterationCounter<FirstIter||IterationCounter>LastIter) {IterationCounter++; continue;}

        //[0] = femto range
        //[1] = fit range
        //[2] = fractions for p
        //[3] = fractions for L
        VARIATIONS[0] = FemtoRangeVars[iFemRan];
        VARIATIONS[1] = FitRangeVars[iFitRan];
        VARIATIONS[2] = pFracVars[ipFrac];
        VARIATIONS[3] = LamFracVars[iLamFrac];
        VARIATIONS[4] = kcVars[ikc];
        Fit_pL(AnalysisObject, OutputFolder,"pp13TeV_HM_March19",Source[iSource],SourceScale[iSourceScal],SourceStability[iSourceStab],
               Potential[iPot],Baseline[iBl],VARIATIONS,ntResult);

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

    printf("This set up had a total of %i iterations.\n",IterationCounter);
    printf("The function call evaluated %i iterations.\n",ComputedIterations);

    NTfile->cd();
    ntResult->Write();

    if(ntResult) delete ntResult;
    if(NTfile) delete NTfile;
}

int PLAMBDA_1_MAIN(int argc, char *argv[]){

    pL_SystematicsHM("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pL/Systematics_050419/",0,
                     6600,6603,false);
    return 0;

    int VARIATIONS[4];for(int i=0;i<4;i++)VARIATIONS[i]=0;

    DLM_CommonAnaFunctions AnalysisObject;
    const TString OutputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pL/050419/";
    //Fit_pL(OutputFolder,"pp13TeV_HM_March19", "Gauss", false, "NLO_Coupled_S","Longbaseline",VARIATIONS);
    //Fit_pL(OutputFolder,"pp13TeV_HM_March19", "McGauss_Reso", false, "NLO_Coupled_S","Longbaseline",VARIATIONS);
    Fit_pL(AnalysisObject,OutputFolder,"pp13TeV_HM_March19", "McLevyNolan_Reso", 0, 0, "NLO_Coupled_S","Longbaseline",VARIATIONS,NULL);

    //Fit_pL(OutputFolder,"pp13TeV_HM_March19", "Gauss", false, "LO_Coupled_S","Longbaseline",VARIATIONS);
    //Fit_pL(OutputFolder,"pp13TeV_HM_March19", "McGauss_Reso", false, "LO_Coupled_S","Longbaseline",VARIATIONS);
    Fit_pL(AnalysisObject,OutputFolder,"pp13TeV_HM_March19", "McLevyNolan_Reso", 0, 0, "LO_Coupled_S","Longbaseline",VARIATIONS,NULL);

    //Fit_pL(OutputFolder,"pp13TeV_HM_March19", "Gauss", false, "NLO_Coupled_S","Longbaseline2",VARIATIONS);
    //Fit_pL(OutputFolder,"pp13TeV_HM_March19", "McGauss_Reso", false, "NLO_Coupled_S","Longbaseline2",VARIATIONS);
    Fit_pL(AnalysisObject,OutputFolder,"pp13TeV_HM_March19", "McLevyNolan_Reso", 0, 0, "NLO_Coupled_S","Longbaseline2",VARIATIONS,NULL);
    //Fit_pL(OutputFolder,"pp13TeV_HM_March19", "Gauss", false, "LO_Coupled_S","Longbaseline2",VARIATIONS);
    //Fit_pL(OutputFolder,"pp13TeV_HM_March19", "McGauss_Reso", false, "LO_Coupled_S","Longbaseline2",VARIATIONS);
    Fit_pL(AnalysisObject,OutputFolder,"pp13TeV_HM_March19", "McLevyNolan_Reso", 0, 0, "LO_Coupled_S","Longbaseline2",VARIATIONS,NULL);

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
