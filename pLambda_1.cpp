#include "pLambda_1.h"

#include "TString.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TGraph.h"

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
void Fit_pp_pL(const TString& DataSample, const TString& SourceType){

    const TString OutputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pp_pL/HM/LongBaseline/";

    //"Norm": only norm
    //"Baseline": baseline but only in the femto region
    //"LongBaseline": baseline extending to higher k
    const TString FittingMode_pp="LongBaseline";
    const TString FittingMode_pL="LongBaseline";

    const bool SameSourceSize = true;
    const bool SameStability = false;

    //the residuals that are not fitted, are assumed to all have a Gaussian source of some size;
    double ResidualSourceSize=0;
    if(DataSample=="pp13TeV_MB_Run2paper") ResidualSourceSize=1.1;
    else if(DataSample=="pPb5TeV_Run2paper") ResidualSourceSize=1.4;
    else if(DataSample=="pp13TeV_HM_March19") ResidualSourceSize=1.4;
    else printf("\033[1;31mERROR:\033[0m The data sample '%s' does not exist\n",DataSample.Data());

    const unsigned Num_pp_pots = 1;
    const unsigned Num_pL_pots = 4;
    const unsigned Num_pXim_pots = 2;

    const TString pp_Pot[Num_pp_pots] =     {"AV18"};
    const bool Use_pp_Pot[Num_pp_pots] =    {true};
    const TString pL_Pot[Num_pL_pots] =     {"LO",  "NLO",  "NLO_Coupled_S",    "Usmani"};
    const bool Use_pL_Pot[Num_pL_pots] =    {true,  false,   true,               false};
    const TString pXim_Pot[Num_pXim_pots] = {"pXim_HALQCD1","pXim_Lattice"};
    const bool Use_pXim_Pot[Num_pXim_pots] =  {true,          false};

    double lambda_pp[4];
    SetUpLambdaPars_pp(DataSample,0,lambda_pp);

    double lambda_pL[5];
    SetUpLambdaPars_pL(DataSample,0,0,lambda_pL);

    double lambda_pXim[5];
    SetUpLambdaPars_pXim(DataSample,0,0,lambda_pXim);

    double* MomBins_pp = NULL;
    double* FitRegion_pp = NULL;
    unsigned NumMomBins_pp;
    SetUpBinning_pp(DataSample,NumMomBins_pp,MomBins_pp,FitRegion_pp);

    double* MomBins_pL = NULL;
    double* FitRegion_pL = NULL;
    unsigned NumMomBins_pL;
    SetUpBinning_pL(DataSample,NumMomBins_pL,MomBins_pL,FitRegion_pL);

    DLM_Ck* Ck_pSigma0 = new DLM_Ck(1,0,NumMomBins_pL,MomBins_pL,Lednicky_gauss_Sigma0);
    Ck_pSigma0->SetSourcePar(0,ResidualSourceSize);

    TH2F* hResolution_pp = GetResolutionMatrix(DataSample,"pp");
    TH2F* hResolution_pL = GetResolutionMatrix(DataSample,"pLambda");

    TH2F* hResidual_pp_pL = GetResidualMatrix("pp","pLambda");
    TH2F* hResidual_pL_pSigma0 = GetResidualMatrix("pLambda","pSigma0");
    TH2F* hResidual_pL_pXim = GetResidualMatrix("pLambda","pXim");

    TH1F* hData_pp = GetAliceExpCorrFun(DataSample,"pp",0);
    TH1F* hData_pL = GetAliceExpCorrFun(DataSample,"pLambda",1);

    TH1F* hDataClever_pL = pL_CleverRebin(hData_pL,NumMomBins_pL,MomBins_pL);

    double lam_pp[4];
    double lam_pL[5];
    double lam_pXim[5];

    SetUpLambdaPars_pp(DataSample,0,lam_pp);
    SetUpLambdaPars_pL(DataSample,0,0,lam_pL);
    SetUpLambdaPars_pXim(DataSample,0,0,lam_pXim);

    TFile* OutputFile = new TFile(OutputFolder+TString::Format("Output_%s_%s.root",DataSample.Data(),SourceType.Data()),"recreate");
    hData_pp->Write();
    hData_pL->Write();
    hDataClever_pL->Write();

    //iterate over all potentials
    for(unsigned uPot_pp=0; uPot_pp<Num_pp_pots; uPot_pp++){
        if(!Use_pp_Pot[uPot_pp]) continue;
        CATS AB_pp;
        AB_pp.SetMomBins(NumMomBins_pp,MomBins_pp);
        SetUpCats_pp(AB_pp,pp_Pot[uPot_pp],SourceType);
        AB_pp.SetNotifications(CATS::nWarning);
        AB_pp.KillTheCat();
        DLM_Ck* Ck_pp = new DLM_Ck(AB_pp.GetNumSourcePars(),0,AB_pp);
        for(unsigned uPot_pL=0; uPot_pL<Num_pL_pots; uPot_pL++){
            if(!Use_pL_Pot[uPot_pL]) continue;
            CATS AB_pL;
            AB_pL.SetMomBins(NumMomBins_pL,MomBins_pL);
            SetUpCats_pL(AB_pL,pL_Pot[uPot_pL],SourceType);
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
                SetUpCats_pXim(AB_pXim,pXim_Pot[uPot_pXim],"Gauss");
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
                    fitter->SetSystem(0,*hData_pp,1,CkDec_pp,
                                  FitRegion_pp[0],FitRegion_pp[1],FitRegion_pp[1],FitRegion_pp[1]);
                }
                else if(FittingMode_pp=="LongBaseline"){
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
                else if(FittingMode_pL=="LongBaseline"){
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
                if(FittingMode_pp=="Baseline"||FittingMode_pp=="LongBaseline")
                    fitter->SetParameter("pp",DLM_Fitter1::p_b,0,-2e-3,2e-3);
                else
                    fitter->FixParameter("pp",DLM_Fitter1::p_b,0);
                if(FittingMode_pL=="Baseline"||FittingMode_pL=="LongBaseline")
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
                }
                else if(SourceType=="CleverLevy_Nolan"||"CleverLevy_Single"||"CleverLevy_Diff"){
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
                else{
                    printf("\033[1;31mERROR:\033[0m '%s' does not exist",SourceType.Data());
                }
                CkDec_pp.Update();
                CkDec_pL.Update();
                CkDec_pSigma0.Update();
                CkDec_pXim.Update();
                fitter->GoBabyGo(true);

                TGraph FitResult_pp;
                FitResult_pp.SetName(TString::Format("pp_%s_%s_%s",pp_Pot[uPot_pp].Data(),pL_Pot[uPot_pL].Data(),pXim_Pot[uPot_pXim].Data()));
                fitter->GetFitGraph(0, FitResult_pp);

                TGraph FitResult_pL;
                FitResult_pL.SetName(TString::Format("pL_%s_%s_%s",pp_Pot[uPot_pp].Data(),pL_Pot[uPot_pL].Data(),pXim_Pot[uPot_pXim].Data()));
                fitter->GetFitGraph(1, FitResult_pL);

                TString Description = TString::Format("%s_%s_%s",pp_Pot[uPot_pp].Data(),pL_Pot[uPot_pL].Data(),pXim_Pot[uPot_pXim].Data());
                //TCanvas* PLOT = Plot_pL_FitResults(hData_pp,hDataClever_pL,Description);

                OutputFile->cd();
                FitResult_pp.Write();
                FitResult_pL.Write();
                //PLOT->Write();

                delete Ck_pXim;
                delete fitter;
                //delete PLOT;
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
    delete hDataClever_pL;
    delete OutputFile;
}



int PLAMBDA_1_MAIN(int argc, char *argv[]){

    //Fit_pp_pL("pp13TeV_HM_March19","Gauss");
    Fit_pp_pL("pp13TeV_HM_March19","CleverLevy_Nolan");

    return 0;
}






