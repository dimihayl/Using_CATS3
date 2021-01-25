
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
#include "TTree.h"
#include "TRandom3.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TAttFill.h"
#include "TSpline.h"

#include "CATS.h"
#include "DLM_CkDecomposition.h"
#include "DLM_CkModels.h"
#include "CommonAnaFunctions.h"
#include "DLM_Fitters.h"
#include "EnvVars.h"

#include "DLM_HistoAnalysis.h"
#include "DLM_SubPads.h"
#include "DLM_Sort.h"
#include "DLM_Integration.h"
#include "DLM_RootWrapper.h"
#include "DLM_Histo.h"
#include "DLM_DecayMatrix.h"
#include "DLM_Potentials.h"
#include "DLM_Source.h"

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

    DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder("/mnt/Ubuntu_Data/CernBox/Sync/CatsFiles");

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
    else if(DataSample=="pp13TeV_HM_March19"||DataSample=="pp13TeV_HM_Dec19") ResidualSourceSize=1.4;
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
    if(FittingMode_pL.Contains("Norm")||FittingMode_pL.Contains("Baseline")||
       (FittingMode_pL.Contains("S2020")&&!FittingMode_pL.Contains("pol3"))){
        fitter->SetSystem(0,*hData_pL,1,CkDec_pL,
                    FitRegion_pL[0],FitRegion_pL[1],FitRegion_pL[1],FitRegion_pL[1]);
    }
    else if(FittingMode_pL.Contains("Longbaseline")||FittingMode_pL.Contains("Spline")||
            FittingMode_pL.Contains("S2020_pol3")){
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

    fitter->GetUnfoldedCk(0,7,OutputFolder+"Unfold_"+OutFileName);
    OutputFile->cd();

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
//! FOR THE PRELIMINARIES, THE DATA SET "pp13TeV_HM_March19" WAS USED, CHANGE IT WHEN CALLING Fit_pL
void pL_SystematicsHM(const TString& OutputFolder, const int& WhichConfiguration, const int& FirstIter, const int& LastIter, const unsigned& OnlyFraction,
                      const int& RANDOMSEED, const int& Mode){

    //const bool Perform_mT_Scan = true;

    printf("Output: %s\n",OutputFolder.Data());
    printf("Iter: %i -> %i\n",FirstIter,LastIter);
    printf("Fraction: %i\n",OnlyFraction);
    printf("Seed: %i\n",RANDOMSEED);
    DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder("/mnt/Ubuntu_Data/CernBox/Sync/CatsFiles");
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
        //the plots for Schleching, which include a Gauss and pol3 (with ^1 par ==0) BL
        //we do not even consider pol0 anymore
        //the source is the ORIGINAL Gauss+Reso, since that was the approved plot at the time
        //we only look at NLO with the coupling included (d-waves) and reduced peak height to 1/3
        case 50 :
            NumSourceVars = 1;
            NumSourceScaleVars = 3;
            NumSourceStabilityVars = 1;
            NumPotVars = 1;
            NumBaselineVars = 3;
            NumFemtoRangeVars = 3;
            NumFitRangeVars = 2;
            NumProtonFracVars = 3;
            NumLambdaFracVars = 9;
            NumKcVars = 1;
            NumMtVars = 1;
            NumDataVars = 45;
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
    else if(WhichConfiguration==50){
        Source[0] = "McGauss_Reso";
        DefSource = 0;
        //SourceScale[0] = -1;//free fit, but within the uncertainties
        //initially it was 0.87 +/- 10 %
        //with wrong systematics: 0.871,0.918,0.965
        SourceScale[0] = 0.853;
        SourceScale[1] = 0.878;
        SourceScale[2] = 0.902;
        // pLambda:   low 0.90218  mean 0.918078 up 0.933344 //pp + pL
        // pLambda:   low 0.872613  mean 0.882538 up 0.892688 //pp
        // new syst: pLambda(1.55):   low 0.853396 mean 0.877717 up 0.90239
        DefSourceScale = 0;
        SourceStability[0] = 0;//dummy
        DefSourceStability = 0;
        Potential[0] = "Chiral_Coupled_SPD";
        DefPot = 0;
        Baseline[0] = "S2020_pol1";
        Baseline[1] = "S2020_pol2";
        Baseline[2] = "S2020_pol3";
        //Baseline[3] = "S2020_Gauss";//not included in the fitter yet
        DefBaseline = 2;

        FemtoRangeVars[0] = 0;//default
        FemtoRangeVars[1] = 1;//down
        FemtoRangeVars[2] = 2;//up
        DefFemtoRange = 0;

        FitRangeVars[0] = 50;//default (only FemtoReg)
        FitRangeVars[1] = 51;//extended fit up to c.a.450 used for the pol3
        DefFitRange = 50;

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

        kcVars[0] = 600;
        DefKc = 0;

        MtVars[0] = -1;
        DefMt = 0;

        for(int iData=0; iData<NumDataVars; iData++){
            DataVars[iData] = TString::Format("_%i",iData);
        }
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
        if((Baseline[iBl]=="Norm"||(Baseline[iBl].Contains("S2020")&&!Baseline[iBl].Contains("pol3")))&&ikc) {continue;}//this variation is only relevant in case of a baseline
        if((Baseline[iBl]=="Norm"||(Baseline[iBl].Contains("S2020")&&!Baseline[iBl].Contains("pol3")))&&iFitRan) {continue;}//this variation is only relevant in case of a baseline
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
            if((Baseline[iBl]=="Norm"||(Baseline[iBl].Contains("S2020")&&!Baseline[iBl].Contains("pol3")))){
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

//! FOR THE PRELIMINARIES, THE DATA SET "pp13TeV_HM_March19" WAS USED, CHANGE IT WHEN CALLING Fit_pL
        Fit_pL(AnalysisObject, OutputFolder, FileName,"pp13TeV_HM_Dec19",DataVars[iData],Source[iSource],SourceScale[iSourceScal],SourceStability[iSourceStab],
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
    DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder("/mnt/Ubuntu_Data/CernBox/Sync/CatsFiles");
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
        if(WhichConfiguration==20) Source[0] = "McGauss_ResoTM";
        else if(WhichConfiguration==21) Source[0] = "McLevyNolan_ResoTM";
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
    lLegend->AddEntry(&mTSyst_pL_LO,"p#minus #Lambda (LO)");
    lLegend->AddEntry(&mTSyst_pL_NLO,"p#minus #Lambda (NLO)");
    //if(mT_pXim)lLegend->AddEntry(mT_pXim,"p#minus #Xi^{#minus} (HAL QCD)");

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
    lLegend_TOT->AddEntry(&mTSyst_pL,"p#minus #Lambda (LO,NLO,Usmani)");
    //if(mT_pXim)lLegend->AddEntry(mT_pXim,"p#minus #Xi^{#minus} (HAL QCD)");

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
    lLegend->AddEntry(&mTSyst_pL_LO,"p#minus #Lambda (LO)");
    lLegend->AddEntry(&mTSyst_pL_NLO,"p#minus #Lambda (NLO)");
    //if(mT_pXim)lLegend->AddEntry(mT_pXim,"p#minus #Xi^{#minus} (HAL QCD)");

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
    lLegend_TOT->AddEntry(&mTSyst_pL,"p#minus #Lambda (LO,NLO,Usmani)");
    //if(mT_pXim)lLegend->AddEntry(mT_pXim,"p#minus #Xi^{#minus} (HAL QCD)");

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
    lLegendSor1->AddEntry(&mTSyst_pL_LO,"p#minus #Lambda (LO)");
    lLegendSor1->AddEntry(&mTSyst_pL_NLO,"p#minus #Lambda (NLO)");
    //if(mT_pXim)lLegendSor1->AddEntry(mT_pXim,"p#minus #Xi^{#minus} (HAL QCD)");

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
    lLegendSor1_TOT->AddEntry(&mTSystSor1_pL,"p#minus #Lambda (LO,NLO,Usmani)");
    //if(mT_pXim)lLegend->AddEntry(mT_pXim,"p#minus #Xi^{#minus} (HAL QCD)");

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

    DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder("/mnt/Ubuntu_Data/CernBox/Sync/CatsFiles");
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
    DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder("/mnt/Ubuntu_Data/CernBox/Sync/CatsFiles");

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



DLM_CkDecomposition* pLambda_1_Dec;
CATS* pLambda_1_Cat;
CATS* pLambda_1_CatFine;

//parameters:
//[0] = source size
//[1] = source alpha
//[2] = cusp strength
//[3] = Ck convergence point
//baseline of the type: []*(1+[]*exp([mean],[sigma]))*(1+[]*x+[]*x*x+[]*x*x*x)
//the par in-front of x^3 is defined such as to represent the k* value at which the pol3 has an extrema
//[4] = norm
//[5],[6],[7] GAUSS
//[8],[9],[10] POL
double dimi_pL_Schleching_FitterBl(double* x, double* par){
    double& MOM = *x;
    double GAUSS = 1.+par[1]*exp(-pow((MOM-par[2])/par[3],2.));
    double POL = 1.+par[4]*MOM+par[5]*MOM*MOM-(2.*par[5]+par[4])/(3.*par[6])*MOM*MOM*MOM;
    return par[0]*GAUSS*POL;
}
double dimi_pL_Schleching_Fitter(double* x, double* par){
    double& MOM = *x;
    pLambda_1_Cat->SetChannelWeight(7,1./4.*par[2]);//1S0 SN(s) -> LN(s)
    pLambda_1_Cat->SetChannelWeight(8,3./4.*par[2]);//3S1 SN(s) -> LN(s)
    pLambda_1_Cat->SetChannelWeight(10,3./4.*par[2]);//3S1 SN(d) -> LN(s)
    pLambda_1_Cat->SetChannelWeight(13,3./20.*par[2]);//3D1 SN(d) -> LN(d)
    pLambda_1_Cat->SetChannelWeight(15,3./20.*par[2]);//3D1 SN(s) -> LN(d)
    pLambda_1_Cat->SetAnaSource(0,par[0],true);
    if(pLambda_1_Cat->GetNumSourcePars()>1){
        pLambda_1_Cat->SetAnaSource(1,par[1],true);
    }
    pLambda_1_Cat->KillTheCat();
    pLambda_1_Dec->GetCk()->SetCutOff(340,par[3]);
    pLambda_1_Dec->Update(true);

    double FEMTO = pLambda_1_Dec->EvalCk(MOM);
    return dimi_pL_Schleching_FitterBl(x,&par[4])*FEMTO;
}

//if par[0]==-1e6 : use splines
//else
//par[0] is an overall normalization
//than we have a pol4 = p0*(1+p1*k+p2*k^2+p3*k^4+p4*k^4), which has 3 free arguments and the following properties
//par4!=0 (pol4 flat at 0)
//	par1,par2 the two extrema, par3 is the p4, par4 is dummy
//par4==0&&par3!=0 (pol3)
//	par1,par2 the two extrema, par3 is p3
//par4==0&&par3==0&&par2!=0 (pol2)
//	par1 is the extrema, par2 is p2
//par4==0&&par3==0&&par2==0&&par1!=0 (pol1)
//	par1 is p1
//to avoid problems with a starting parameter of zero, to switch the order of the par we use -1e6 as a value
//a Mathematica computation of the equations is in your Femto folder
double dimi_pL_May2020_FitterBl(double* x, double* par){
    double& k = *x;
    double& p0 = par[0];
    //splines
    if(p0==-1e6){
        //double Spline = DLM_FITTER2_FUNCTION_SPLINE3(x,&par[1]);
        //printf("hi %f\n",Spline);
        return DLM_FITTER2_FUNCTION_SPLINE3(x,&par[1]);
    }
    //constrained polynomials
    else{
        double p1;
        double p2;
        double p3;
        double p4;
        if(par[4]!=-1e6){
            p4 = par[3];
            p3 = -4./3.*(par[1]+par[2])*p4;
            p2 = 2.*par[1]*par[2]*p4;
            p1 = 0;
        }
        else if(par[3]!=-1e6){
            p4 = 0;
            p3 = par[3];
            p2 = -1.5*(par[1]+par[2])*p3;
            p1 = 3.*par[1]*par[2]*p3;
        }
        else if(par[2]!=-1e6){
            p4 = 0;
            p3 = 0;
            p2 = par[2];
            p1 = -2.*par[1]*p2;
        }
        else{
            p4 = 0;
            p3 = 0;
            p2 = 0;
            p1 = par[1];
        }
        return p0*(1.+p1*k+p2*pow(k,2)+p3*pow(k,3)+p4*pow(k,4));
    }
}
double dimi_pL_May2020_FitterFemto(double* x, double* par){
    double& MOM = *x;
    const bool FINE = (MOM>264&&MOM<300);

    pLambda_1_Cat->SetChannelWeight(7,1./4.*par[2]);//1S0 SN(s) -> LN(s)
    pLambda_1_Cat->SetChannelWeight(8,3./4.*par[2]);//3S1 SN(s) -> LN(s)
    pLambda_1_Cat->SetChannelWeight(10,3./4.*par[2]);//3S1 SN(d) -> LN(s)
    pLambda_1_Cat->SetChannelWeight(13,3./20.*par[2]);//3D1 SN(d) -> LN(d)
    pLambda_1_Cat->SetChannelWeight(15,3./20.*par[2]);//3D1 SN(s) -> LN(d)
    pLambda_1_Cat->SetAnaSource(0,par[0],true);
    if(pLambda_1_Cat->GetNumSourcePars()>1){
        pLambda_1_Cat->SetAnaSource(1,par[1],true);
    }
    pLambda_1_Cat->KillTheCat();
    pLambda_1_Dec->GetCk()->SetCutOff(par[3],par[4]);

    //just around the cusp, we make the integrated bin value
    //we ASSUME there is no smearing effect present, so this works well only for unfolded data

    if(FINE){
//printf("FINE:\n");
        pLambda_1_CatFine->SetChannelWeight(7,1./4.*par[2]);//1S0 SN(s) -> LN(s)
        pLambda_1_CatFine->SetChannelWeight(8,3./4.*par[2]);//3S1 SN(s) -> LN(s)
        pLambda_1_CatFine->SetChannelWeight(10,3./4.*par[2]);//3S1 SN(d) -> LN(s)
        pLambda_1_CatFine->SetChannelWeight(13,3./20.*par[2]);//3D1 SN(d) -> LN(d)
        pLambda_1_CatFine->SetChannelWeight(15,3./20.*par[2]);//3D1 SN(s) -> LN(d)
        pLambda_1_CatFine->SetAnaSource(0,par[0],true);
        if(pLambda_1_CatFine->GetNumSourcePars()>1){
            pLambda_1_CatFine->SetAnaSource(1,par[1],true);
        }
        pLambda_1_CatFine->KillTheCat();

        const unsigned NumIntSteps = 7;
        const double IntWidth = 6;
        //only for the bin near the cusp, we make the avg. bin value based on multiple evaluation points, else its the bin center
        const unsigned NumAvgSteps = 9;

        double CkGenuine_Fine=0;
        double mom_val;
        //+/- this value, so 1/2 bin width. Should work approx. for different (but similar, say 8-16 MeV) binning as well.
        const double AvgRange = 6;
        const double AvgStep = (2.*AvgRange)/double(NumAvgSteps-1);
        for(unsigned uBin=0; uBin<NumAvgSteps; uBin++){
            mom_val = (MOM-AvgRange)+double(uBin)*AvgStep;
            CkGenuine_Fine += pLambda_1_CatFine->EvalCorrFun(mom_val);
        }
        //up to know, this is the genuine Ck, without any feed-down etc
        CkGenuine_Fine /= double(NumAvgSteps);
        //we will use the result from the full decomposition, to estimate the effect of the feeddown.
        //here the assumption (should be super good) is that apart from the cusp, all other contributions are
        //linear inside the bin in question
        pLambda_1_Dec->Update(true,true);

        //this is the contribution (C-1) from non-genuine stuff, scaled for lambda parameters
        double SignalNonGenuine = pLambda_1_Dec->EvalSignal(MOM)-pLambda_1_Dec->EvalSignalMain(MOM);
        //this is the same observable, but for the Genuine part
        double SignalGenuineFine = (CkGenuine_Fine-1.)*pLambda_1_Dec->GetLambdaMain();
//double SignalGenuine = (pLambda_1_Cat->EvalCorrFun(MOM)-1.)*pLambda_1_Dec->GetLambdaMain();
//double LambdaMain = pLambda_1_Dec->GetLambdaMain();
//DLM_Ck WTF_CK(pLambda_1_Cat->GetNumSourcePars(),0,*pLambda_1_Cat);
//WTF_CK.Update(true);
//printf(" QA Ck:\n");
//printf("  C(%.0f) = %.4f\n",MOM,pLambda_1_Dec->EvalCk(MOM));
//printf("  Ck = %.4f\n",pLambda_1_Dec->GetCk()->Eval(MOM));
//printf("   CL = %.4f -> %.4f\n",pLambda_1_Dec->EvalSignalMain(MOM)+LambdaMain,(pLambda_1_Dec->EvalSignalMain(MOM)+LambdaMain)/LambdaMain);
//printf("  CATS = %.4f -> %.4f\n",pLambda_1_Cat->EvalCorrFun(MOM)*LambdaMain,pLambda_1_Cat->EvalCorrFun(MOM));
//printf("  WTF_CK = %.4f\n",WTF_CK.Eval(MOM));
//printf(" QA (%.0f):\n",MOM);
//printf("  Default value: %.4f\n",pLambda_1_Dec->EvalCk(MOM));
//printf("       QA value: %.4f\n",SignalGenuine+SignalNonGenuine+1.);
//printf("  Refined value: %.4f\n",SignalGenuineFine+SignalNonGenuine+1.);
//usleep(200e3);

        return SignalGenuineFine+SignalNonGenuine+1.;
    }
    else{
        pLambda_1_Dec->Update(true);
        return pLambda_1_Dec->EvalCk(MOM);
    }
}
double dimi_pL_May2020_Fitter(double* x, double* par){
    return dimi_pL_May2020_FitterFemto(x,par)*dimi_pL_May2020_FitterBl(x,&par[5]);
}

//the radius is fixed, and we use the latest Core+Reso
//if the seed is 1, the default variations are computed first
void pL_SchlechingSystematics(const unsigned SEED, const unsigned NumVars_PerBL, const int& POT_VAR){


    //for Schleching, only
    // 0 = NLO13 600 MeV
    //10 = NLO19 600 MeV
    //int POT_VAR = 10;

    const bool DataSyst = false;
    //Dimi is pol3 up to higher fit range
    enum BLTYPE { pol0,pol1,pol2,pol3,gaus,Dimi };
    const unsigned NumBlTypes = 6;
    bool blType[NumBlTypes];
    blType[pol0] = true;
    blType[pol1] = true;
    blType[pol2] = true;
    blType[pol3] = true;
    blType[gaus] = false;
    blType[Dimi] = true;

    //const unsigned NumVars_PerBL = DataSyst?1024:128;

    // new syst: pLambda(1.55):   low 0.853396 mean 0.877717 up 0.90239
    //the values below is for the updated epos version of the source
    const double SourceSize = (1.05676+0.980842)*0.5;
    const double SourceSizeErr = (1.05676-0.980842)*0.5;

    double ResidualSourceSizeSigma0=1.25;
    double ResidualSourceSizeXi=1.0;

    unsigned CompletedVars[NumBlTypes];
    CompletedVars[pol0] = 0;
    CompletedVars[pol1] = 0;
    CompletedVars[pol2] = 0;
    CompletedVars[pol3] = 0;
    CompletedVars[gaus] = 0;
    CompletedVars[Dimi] = 0;


    TRandom3 rangen(SEED);

    TString DataSample = "pp13TeV_HM_Dec19";
    //TString DataSample = "pp13TeV_HM_RotPhiDec19";

    TString OutputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/pL_AllHaide_May2020/";
    TString OutFileName = TString::Format("Output_%s_POT%i_%u.root",DataSample.Data(),POT_VAR,SEED);
    TFile* OutputFile = new TFile(OutputFolder+OutFileName,"recreate");

    double* MomBins_pL = NULL;
    double* FitRegion_pL = NULL;
    const unsigned NumBinsCk = 45; // for dlm_ck
    const double BinWidthCk = 12;
    const double MaxBinValCk = double(NumBinsCk)*BinWidthCk;
    unsigned NumMomBins_pL = 30;// for CATS
    const double MaxBinValCats = double(NumMomBins_pL)*BinWidthCk;
    unsigned uDummy;

    //used for plotting the systematics, keeps track of the mean and stdv of the value of the fit in each bin
    //the last entry is independent of the baseline type
    float** Fit_Mean = new float* [NumBlTypes+1];
    float** Fit_Stdv = new float* [NumBlTypes+1];
    float** Fit_Min = new float* [NumBlTypes+1];
    float** Fit_Max = new float* [NumBlTypes+1];
    unsigned** Fit_Entries = new unsigned* [NumBlTypes+1];
    float** Bl_Mean = new float* [NumBlTypes+1];
    float** Bl_Stdv = new float* [NumBlTypes+1];
    float** Bl_Min = new float* [NumBlTypes+1];
    float** Bl_Max = new float* [NumBlTypes+1];

    float** Ratio_Min = new float* [NumBlTypes+1];
    float** Ratio_Max = new float* [NumBlTypes+1];


    //TH1F* hRadius = new TH1F("hRadius","hRadius",256,0.9,1.1);
    //TH1F* hAlpha = new TH1F("hAlpha","hAlpha",256,0.95,2.05);
    TH1F** hnSigma = new TH1F* [NumBlTypes+1];
    //TH1F** hRatio = new TH1F* [NumBlTypes+1];
    TH1F** hRadius = new TH1F* [NumBlTypes+1];
    TH1F** hCS = new TH1F* [NumBlTypes+1];

//[0] = source size
//[1] = source alpha
//[2] = cusp strength
//[3] = Ck convergence point
//baseline of the type: []*(1+[]*exp([mean],[sigma]))*(1+[]*x+[]*x*x+[]*x*x*x)
//the par in-front of x^3 is defined such as to represent the k* value at which the pol3 has an extrema
//[4] = norm
//[5],[6],[7] GAUSS
//[8],[9],[10] POL
    TNtuple* ntResult = new TNtuple("ntResult", "ntResult",
            "SourceSize:SourceAlpha:CuspStrength:CkConv:Norm:GaussAmpl:GaussMu:GaussSig:p1:p2:p3:"
            "POT_VAR:WhichBl:WhichFemtoRange:WhichFitRange:WhichProtonVar:WhichLambdaVar:WhichData:pval");

    for(unsigned uBl=0; uBl<=NumBlTypes; uBl++){
        Fit_Mean[uBl] = new float [NumBinsCk];
        Fit_Stdv[uBl] = new float [NumBinsCk];
        Fit_Min[uBl] = new float [NumBinsCk];
        Fit_Max[uBl] = new float [NumBinsCk];
        Bl_Mean[uBl] = new float [NumBinsCk];
        Bl_Stdv[uBl] = new float [NumBinsCk];
        Bl_Min[uBl] = new float [NumBinsCk];
        Bl_Max[uBl] = new float [NumBinsCk];
        Ratio_Min[uBl] = new float [NumBinsCk];
        Ratio_Max[uBl] = new float [NumBinsCk];
        Fit_Entries[uBl] = new unsigned [NumBinsCk];
        for(unsigned uBin=0; uBin<NumBinsCk; uBin++){
            Fit_Mean[uBl][uBin] = 0;
            Fit_Stdv[uBl][uBin] = 0;
            Fit_Min[uBl][uBin] = 10000;
            Fit_Max[uBl][uBin] = 0;
            Bl_Mean[uBl][uBin] = 0;
            Bl_Stdv[uBl][uBin] = 0;
            Bl_Min[uBl][uBin] = 10000;
            Bl_Max[uBl][uBin] = 0;
            Ratio_Min[uBl][uBin] = 10000;
            Ratio_Max[uBl][uBin] = 0;
            Fit_Entries[uBl][uBin] = 0;
        }
        TString hnSigmaName = "hnSigma";
        if(uBl<NumBlTypes) hnSigmaName = TString::Format("hnSigma_%u",uBl);
        hnSigma[uBl] = new TH1F (hnSigmaName,hnSigmaName,1024,0,20.0);
        //TString hRatioName = "hRatio";
        //if(uBl<NumBlTypes) hRatioName = TString::Format("hRatio_%u",uBl);
        //hRatio[uBl] = new TH1F (hRatioName,hRatioName,38,0,456);
        TString hRadiusName = "hRadius";
        if(uBl<NumBlTypes) hRadiusName = TString::Format("hRadius_%u",uBl);
        hRadius[uBl] = new TH1F (hRadiusName,hRadiusName,256,0.6,1.4);
        TString hCSName = "hCS";
        if(uBl<NumBlTypes) hCSName = TString::Format("hCS_%u",uBl);
        hCS[uBl] = new TH1F (hCSName,hCSName,256,0.0,1.0);
    }


    DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder("/mnt/Ubuntu_Data/CernBox/Sync/CatsFiles");
    AnalysisObject.SetUpBinning_pL(DataSample,uDummy,MomBins_pL,FitRegion_pL,0,0);

    DLM_Ck* Ck_pSigma0 = new DLM_Ck(1,0,NumBinsCk,0,MaxBinValCk,Lednicky_gauss_Sigma0);
    Ck_pSigma0->SetSourcePar(0,ResidualSourceSizeSigma0);

    TH2F* hResolution_pL = AnalysisObject.GetResolutionMatrix(DataSample,"pLambda");
    TH2F* hResidual_pL_pSigma0 = AnalysisObject.GetResidualMatrix("pLambda","pSigma0");
    TH2F* hResidual_pL_pXim = AnalysisObject.GetResidualMatrix("pLambda","pXim");

    //TString SourceDescription = "Gauss";
    TString SourceDescription = "McGauss_ResoTM";
    //TString SourceDescription = "McLevy_ResoTM";


    //hData_pL->Write();

    CATS AB_pL;
    DLM_Ck* Ck_pL;
    AB_pL.SetMomBins(NumMomBins_pL,0,MaxBinValCats);
    AnalysisObject.SetUpCats_pL(AB_pL,"Chiral_Coupled_SPD",SourceDescription,POT_VAR,202);//NLO_Coupled_S

    double CUSP_WEIGHT = 0.33;//0.54
    AB_pL.SetChannelWeight(7,1./4.*CUSP_WEIGHT);//1S0 SN(s) -> LN(s)
    AB_pL.SetChannelWeight(8,3./4.*CUSP_WEIGHT);//3S1 SN(s) -> LN(s)
    AB_pL.SetChannelWeight(10,3./4.*CUSP_WEIGHT);//3S1 SN(d) -> LN(s)
    AB_pL.SetChannelWeight(13,3./20.*CUSP_WEIGHT);//3D1 SN(d) -> LN(d)
    AB_pL.SetChannelWeight(15,3./20.*CUSP_WEIGHT);//3D1 SN(s) -> LN(d)

    AB_pL.SetAnaSource(0,1.4);
    if(SourceDescription.Contains("Mc")){
        AB_pL.SetAnaSource(0,1.10);//c.a. 10% smaller compared to p-p due to the mT scaling
        AB_pL.SetAnaSource(1,2.0);
    }
    AB_pL.SetNotifications(CATS::nError);
    AB_pL.KillTheCat();
    Ck_pL = new DLM_Ck(AB_pL.GetNumSourcePars(),0,AB_pL,NumBinsCk,0,MaxBinValCk);
    Ck_pL->SetSourcePar(0,AB_pL.GetAnaSourcePar(0));
    if(SourceDescription.Contains("Mc")){
        Ck_pL->SetSourcePar(0,AB_pL.GetAnaSourcePar(0));
        Ck_pL->SetSourcePar(1,AB_pL.GetAnaSourcePar(1));
    }
    Ck_pL->SetCutOff(340,600);

    CATS AB_pXim;
    //same binning as pL, as we only use pXim as feed-down
    AB_pXim.SetMomBins(NumMomBins_pL,0,MaxBinValCats);
    AnalysisObject.SetUpCats_pXim(AB_pXim,"pXim_HALQCD1","Gauss");
    AB_pXim.SetAnaSource(0,ResidualSourceSizeXi);
    AB_pXim.SetNotifications(CATS::nError);
    AB_pXim.KillTheCat();
    DLM_Ck* Ck_pXim = new DLM_Ck(AB_pXim.GetNumSourcePars(),0,AB_pXim,NumBinsCk,0,MaxBinValCk);
    Ck_pL->Update();
    Ck_pSigma0->Update();
    Ck_pXim->Update();

    bool WHILE_CONDITON = false;
    for(unsigned uBl=0; uBl<NumBlTypes; uBl++){
        WHILE_CONDITON += (CompletedVars[uBl]<NumVars_PerBL&&blType[uBl]);
    }
    bool DefaultDataSaved = false;

    unsigned Progress = 0;

    while(WHILE_CONDITON){
        const unsigned WhichBl = rangen.Integer(NumBlTypes);
        if(!blType[WhichBl] || CompletedVars[WhichBl]>=NumVars_PerBL) continue;

        printf("\r\033[K Progress=%u",Progress);
        cout << flush;

        bool DefaultVariation = true;
        if(CompletedVars[WhichBl]||SEED!=1) DefaultVariation = false;

        //printf("WhichBl = %u\n",WhichBl);
        //printf(" DefVar = %i\n",DefaultVariation);

        unsigned WhichFemtoRange = rangen.Integer(3);
        if(DefaultVariation) WhichFemtoRange = 0;
//WhichFemtoRange = 0;
        unsigned WhichFitRange = 50;
        //if(WhichBl==3) WhichFitRange = 50+rangen.Integer(2);
        if(WhichBl==5) WhichFitRange = 51;
        if(DefaultVariation) WhichFitRange = 50;
        AnalysisObject.SetUpBinning_pL(DataSample,uDummy,MomBins_pL,FitRegion_pL,WhichFemtoRange,WhichFitRange);
//WhichFitRange = 50;
        //the cusp strength
        //0 is 33%
        //1 is 27% (20% lower, kind of compatible with experiment)
        //2 is 40% (20% larger)
        //3 is free fit
        unsigned WhichCuspStrength = rangen.Integer(3);
        if(DefaultVariation) WhichCuspStrength = 0;
//WhichCuspStrength=3;
//WhichCuspStrength = 0;
        switch(WhichCuspStrength){
            case 0 : CUSP_WEIGHT = 0.33; break;
            case 1 : CUSP_WEIGHT = 0.27; break;
            case 2 : CUSP_WEIGHT = 0.40; break;
            default : CUSP_WEIGHT = 0.33; break;
        }

//NO FREE
        unsigned WhichSourceRad = rangen.Integer(3);
        double SourceRad;
        if(DefaultVariation) {WhichSourceRad=0;}
//WhichSourceRad=3;
//WhichSourceRad=CompletedVars[WhichBl];

        switch(WhichSourceRad){
            case 0 : SourceRad = SourceSize; break;
            case 1 : SourceRad = SourceSize-SourceSizeErr; break;
            case 2 : SourceRad = SourceSize+SourceSizeErr; break;
            default : SourceRad = SourceSize; break;
        }

        double lam_pL[5];
        double lam_pXim[5];

        int WhichProtonVar = rangen.Integer(3);
        if(DefaultVariation) WhichProtonVar = 0;
        int WhichLambdaVar = rangen.Integer(3);
//WhichProtonVar=0;
//WhichLambdaVar=0;
        if(DefaultVariation) WhichLambdaVar = 0;
        AnalysisObject.SetUpLambdaPars_pL(DataSample,WhichProtonVar,WhichLambdaVar,lam_pL);
        AnalysisObject.SetUpLambdaPars_pXim(DataSample,0,0,lam_pXim);


        DLM_CkDecomposition CkDec_pL("pLambda",4,*Ck_pL,hResolution_pL);
        DLM_CkDecomposition CkDec_pSigma0("pSigma0",0,*Ck_pSigma0,NULL);
        DLM_CkDecomposition CkDec_pXim("pXim",2,*Ck_pXim,NULL);

        //CkDec_pL.AddContribution(0,lam_pL[1],DLM_CkDecomposition::cFeedDown,&CkDec_pSigma0,hResidual_pL_pSigma0);
        //CkDec_pL.AddContribution(1,lam_pL[2],DLM_CkDecomposition::cFeedDown,&CkDec_pXim,hResidual_pL_pXim);
        CkDec_pL.AddContribution(0,lam_pL[1],DLM_CkDecomposition::cFeedDown);
        CkDec_pL.AddContribution(1,lam_pL[2],DLM_CkDecomposition::cFeedDown);
        CkDec_pL.AddContribution(2,lam_pL[3],DLM_CkDecomposition::cFeedDown);
        CkDec_pL.AddContribution(3,lam_pL[4],DLM_CkDecomposition::cFake);//0.03

        //for Xim we simplify a bit and take ALL feed-down as flat
        CkDec_pXim.AddContribution(0,lam_pXim[1]+lam_pXim[2]+lam_pXim[3],DLM_CkDecomposition::cFeedDown);
        CkDec_pXim.AddContribution(1,lam_pXim[4],DLM_CkDecomposition::cFake);

        CkDec_pL.Update();

        unsigned WhichData = rangen.Integer(45);
        if(DefaultVariation||DataSyst==false) WhichData = 0;
//WhichData = 0;
        TString DataVar = TString::Format("_%i",WhichData);
        //TH1F* hData_pL = AnalysisObject.GetAliceExpCorrFun(DataSample,"pLambda","_0",1,false,-1);
        TH1F* hData_pL = AnalysisObject.GetAliceExpCorrFun(DataSample,"pLambda",DataVar,2,false,-1);
        if(!DefaultDataSaved){
            OutputFile->cd();
            TH1F *histoCopy = (TH1F*)hData_pL->Clone("hData");
            histoCopy->Write();
            delete histoCopy;
            DefaultDataSaved = true;
        }

        pLambda_1_Dec = &CkDec_pL;
        pLambda_1_Cat = &AB_pL;

        //parameters:
        //[0] = source size
        //[1] = source alpha
        //[2] = cusp strength
        //[3] = Ck convergence point
        //baseline of the type: []*(1+[]*exp([mean],[sigma]))*(1+[]*x+[]*x*x+[]*x*x*x)
        //the par in-front of x^3 is defined such as to represent the k* value at which the pol3 has an extrema
        //[4] = norm
        //[5],[6],[7] GAUSS
        //[8],[9],[10] POL
        TF1* fit_pL = new TF1("fit_pL",dimi_pL_Schleching_Fitter,FitRegion_pL[0],FitRegion_pL[3],11);

        //source
        if(WhichSourceRad!=3){
            fit_pL->FixParameter(0,SourceRad);
        }
        else{
            fit_pL->SetParameter(0,SourceRad);
            fit_pL->SetParLimits(0,SourceRad*0.8,SourceRad*1.2);
        }
        fit_pL->FixParameter(1,2.0);

        //cusp
        if(WhichCuspStrength==3){
            fit_pL->SetParameter(2,CUSP_WEIGHT);
            fit_pL->SetParLimits(2,0,1.);
        }
        else fit_pL->FixParameter(2,CUSP_WEIGHT);

        //cutoff
        fit_pL->FixParameter(3,600);

        //norm
        fit_pL->SetParameter(4,1);

        //gauss and pol
        if(WhichBl!=4){
            //gauss
            fit_pL->FixParameter(5,0);
            fit_pL->FixParameter(6,0);
            fit_pL->FixParameter(7,300);

            //pol
            fit_pL->SetParameter(8,0);
            fit_pL->SetParLimits(8,-0.01,0.01);
            if(WhichBl==3||WhichBl==5) fit_pL->FixParameter(8,0);
            fit_pL->SetParameter(9,-2.75e-7);
            fit_pL->SetParameter(10,250);
            fit_pL->SetParLimits(10,150,350);
            if(WhichBl<=2) fit_pL->FixParameter(10,1e9);
            if(WhichBl<=1) fit_pL->FixParameter(9,0);
            if(WhichBl<=0) fit_pL->FixParameter(8,0);
        }
        else{
            //gauss
            fit_pL->SetParameter(5,0.05);
            fit_pL->FixParameter(6,0);
            fit_pL->SetParameter(7,300);
            fit_pL->SetParLimits(7,100,500);

            //pol
            fit_pL->FixParameter(8,0);
            fit_pL->FixParameter(9,0);
            fit_pL->FixParameter(10,300);
        }
        hData_pL->Fit(fit_pL,"Q, S, N, R, M");
        fit_pL->SetNpx(1024);

        TF1* fit_pL_Bl = new TF1("fit_pL_Bl",dimi_pL_Schleching_FitterBl,FitRegion_pL[0],FitRegion_pL[3],7);
        for(unsigned uPar=0; uPar<7; uPar++) fit_pL_Bl->FixParameter(uPar,fit_pL->GetParameter(4+uPar));

        for(unsigned uBin=0; uBin<NumBinsCk; uBin++){
            if(Ck_pL->GetBinCenter(0,uBin)>FitRegion_pL[3]) continue;
            //if(Ck_pL->GetBinCenter(0,uBin)>430){
            //    printf("b = %u\n",uBin);
            //    printf("bc = %f\n",Ck_pL->GetBinCenter(0,uBin));
            //    printf("fr = %f\n",FitRegion_pL[3]);
            //}

            Fit_Mean[WhichBl][uBin] += fit_pL->Eval(Ck_pL->GetBinCenter(0,uBin));
            Fit_Stdv[WhichBl][uBin] += fit_pL->Eval(Ck_pL->GetBinCenter(0,uBin))*fit_pL->Eval(Ck_pL->GetBinCenter(0,uBin));
            if(Fit_Min[WhichBl][uBin]>fit_pL->Eval(Ck_pL->GetBinCenter(0,uBin)))
                Fit_Min[WhichBl][uBin] = fit_pL->Eval(Ck_pL->GetBinCenter(0,uBin));
            if(Fit_Max[WhichBl][uBin]<fit_pL->Eval(Ck_pL->GetBinCenter(0,uBin)))
                Fit_Max[WhichBl][uBin] = fit_pL->Eval(Ck_pL->GetBinCenter(0,uBin));
            Fit_Entries[WhichBl][uBin] ++;

            Bl_Mean[WhichBl][uBin] += fit_pL_Bl->Eval(Ck_pL->GetBinCenter(0,uBin));
            Bl_Stdv[WhichBl][uBin] += fit_pL_Bl->Eval(Ck_pL->GetBinCenter(0,uBin))*fit_pL_Bl->Eval(Ck_pL->GetBinCenter(0,uBin));
            if(Bl_Min[WhichBl][uBin]>fit_pL->Eval(Ck_pL->GetBinCenter(0,uBin)))
                Bl_Min[WhichBl][uBin] = fit_pL->Eval(Ck_pL->GetBinCenter(0,uBin));
            if(Bl_Max[WhichBl][uBin]<fit_pL->Eval(Ck_pL->GetBinCenter(0,uBin)))
                Bl_Max[WhichBl][uBin] = fit_pL->Eval(Ck_pL->GetBinCenter(0,uBin));

            if(Ratio_Min[WhichBl][uBin]>hData_pL->GetBinContent(uBin+1)/fit_pL->Eval(Ck_pL->GetBinCenter(0,uBin)))
                Ratio_Min[WhichBl][uBin] = hData_pL->GetBinContent(uBin+1)/fit_pL->Eval(Ck_pL->GetBinCenter(0,uBin));
            if(Ratio_Max[WhichBl][uBin]<hData_pL->GetBinContent(uBin+1)/fit_pL->Eval(Ck_pL->GetBinCenter(0,uBin)))
                Ratio_Max[WhichBl][uBin] = hData_pL->GetBinContent(uBin+1)/fit_pL->Eval(Ck_pL->GetBinCenter(0,uBin));


            Fit_Mean[NumBlTypes][uBin] += fit_pL->Eval(Ck_pL->GetBinCenter(0,uBin));
            Fit_Stdv[NumBlTypes][uBin] += fit_pL->Eval(Ck_pL->GetBinCenter(0,uBin))*fit_pL->Eval(Ck_pL->GetBinCenter(0,uBin));
            if(Fit_Min[NumBlTypes][uBin]>fit_pL->Eval(Ck_pL->GetBinCenter(0,uBin)))
                Fit_Min[NumBlTypes][uBin] = fit_pL->Eval(Ck_pL->GetBinCenter(0,uBin));
            if(Fit_Max[NumBlTypes][uBin]<fit_pL->Eval(Ck_pL->GetBinCenter(0,uBin)))
                Fit_Max[NumBlTypes][uBin] = fit_pL->Eval(Ck_pL->GetBinCenter(0,uBin));
            Fit_Entries[NumBlTypes][uBin] ++;

            Bl_Mean[NumBlTypes][uBin] += fit_pL_Bl->Eval(Ck_pL->GetBinCenter(0,uBin));
            Bl_Stdv[NumBlTypes][uBin] += fit_pL_Bl->Eval(Ck_pL->GetBinCenter(0,uBin))*fit_pL_Bl->Eval(Ck_pL->GetBinCenter(0,uBin));
            if(Bl_Min[NumBlTypes][uBin]>fit_pL->Eval(Ck_pL->GetBinCenter(0,uBin)))
                Bl_Min[NumBlTypes][uBin] = fit_pL->Eval(Ck_pL->GetBinCenter(0,uBin));
            if(Bl_Max[NumBlTypes][uBin]<fit_pL->Eval(Ck_pL->GetBinCenter(0,uBin)))
                Bl_Max[NumBlTypes][uBin] = fit_pL->Eval(Ck_pL->GetBinCenter(0,uBin));
            if(Ratio_Min[NumBlTypes][uBin]>hData_pL->GetBinContent(uBin+1)/fit_pL->Eval(Ck_pL->GetBinCenter(0,uBin)))
                Ratio_Min[NumBlTypes][uBin] = hData_pL->GetBinContent(uBin+1)/fit_pL->Eval(Ck_pL->GetBinCenter(0,uBin));
            if(Ratio_Max[NumBlTypes][uBin]<hData_pL->GetBinContent(uBin+1)/fit_pL->Eval(Ck_pL->GetBinCenter(0,uBin)))
                Ratio_Max[NumBlTypes][uBin] = hData_pL->GetBinContent(uBin+1)/fit_pL->Eval(Ck_pL->GetBinCenter(0,uBin));


        }

        hnSigma[WhichBl]->Fill(sqrt(2)*TMath::ErfcInverse(fit_pL->GetProb()));
        hnSigma[NumBlTypes]->Fill(sqrt(2)*TMath::ErfcInverse(fit_pL->GetProb()));

        hRadius[WhichBl]->Fill(fit_pL->GetParameter(0));
        hRadius[NumBlTypes]->Fill(fit_pL->GetParameter(0));

        hCS[WhichBl]->Fill(fit_pL->GetParameter(2));
        hCS[NumBlTypes]->Fill(fit_pL->GetParameter(2));

        OutputFile->cd();

        //hData_pL->Write();
        //fit_pL->Write();
        //fit_pL_Bl->Write();

        Float_t* buffer = new Float_t [19];
        buffer[0] = fit_pL->GetParameter(0);
        buffer[1] = fit_pL->GetParameter(1);
        buffer[2] = fit_pL->GetParameter(2);
        buffer[3] = fit_pL->GetParameter(3);
        buffer[4] = fit_pL->GetParameter(4);
        buffer[5] = fit_pL->GetParameter(5);
        buffer[6] = fit_pL->GetParameter(6);
        buffer[7] = fit_pL->GetParameter(7);
        buffer[8] = fit_pL->GetParameter(8);
        buffer[9] = fit_pL->GetParameter(9);
        buffer[10] = fit_pL->GetParameter(10);
        buffer[11] = POT_VAR;
        buffer[12] = WhichBl;
        buffer[13] = WhichFemtoRange;
        buffer[14] = WhichFitRange;
        buffer[15] = WhichProtonVar;
        buffer[16] = WhichLambdaVar;
        buffer[17] = WhichData;
        buffer[18] = fit_pL->GetProb();

        ntResult->Fill(buffer);

        delete fit_pL;
        delete fit_pL_Bl;

        CompletedVars[WhichBl]++;
        Progress++;
        WHILE_CONDITON = false;
        for(unsigned uBl=0; uBl<NumBlTypes; uBl++){
            WHILE_CONDITON += (CompletedVars[uBl]<NumVars_PerBL&&blType[uBl]);
        }

        delete [] buffer;
    }
printf("\n");

    for(unsigned uBl=0; uBl<=NumBlTypes; uBl++){
        for(unsigned uBin=0; uBin<NumBinsCk; uBin++){
            if(Fit_Entries[uBl][uBin]){
                Fit_Mean[uBl][uBin] /= double(Fit_Entries[uBl][uBin]);
                Fit_Stdv[uBl][uBin] /= double(Fit_Entries[uBl][uBin]);
                Bl_Mean[uBl][uBin] /= double(Fit_Entries[uBl][uBin]);
                Bl_Stdv[uBl][uBin] /= double(Fit_Entries[uBl][uBin]);
            }
        }
    }

    OutputFile->cd();
    //ntResult->Write("",TObject::kOverwrite);
    ntResult->Write();

    TGraphErrors* graph_Fit = new TGraphErrors[NumBlTypes+1];
    TGraphErrors* graph_Bl = new TGraphErrors[NumBlTypes+1];

    TGraphErrors* graph_FitMAX = new TGraphErrors[NumBlTypes+1];
    TGraphErrors* graph_BlMAX = new TGraphErrors[NumBlTypes+1];
    TGraphErrors* graph_RatioMAX = new TGraphErrors[NumBlTypes+1];


    for(unsigned uBl=0; uBl<=NumBlTypes; uBl++){
        if(uBl<NumBlTypes){
            graph_Fit[uBl].SetName(TString::Format("graph_Fit_%u",uBl));
            graph_Bl[uBl].SetName(TString::Format("graph_Bl_%u",uBl));
            graph_FitMAX[uBl].SetName(TString::Format("graph_FitMAX_%u",uBl));
            graph_BlMAX[uBl].SetName(TString::Format("graph_BlMAX_%u",uBl));
            graph_RatioMAX[uBl].SetName(TString::Format("graph_RatioMAX_%u",uBl));
        }
        else{
            graph_Fit[uBl].SetName(TString::Format("graph_Fit"));
            graph_Bl[uBl].SetName(TString::Format("graph_Bl"));
            graph_FitMAX[uBl].SetName(TString::Format("graph_FitMAX"));
            graph_BlMAX[uBl].SetName(TString::Format("graph_BlMAX"));
            graph_RatioMAX[uBl].SetName(TString::Format("graph_RatioMAX"));
        }
        for(unsigned uBin=0; uBin<NumBinsCk; uBin++){
            graph_Fit[uBl].SetPoint(uBin,Ck_pL->GetBinCenter(0,uBin),Fit_Mean[uBl][uBin]);
            graph_Fit[uBl].SetPointError(uBin,0,sqrt(Fit_Stdv[uBl][uBin]-Fit_Mean[uBl][uBin]*Fit_Mean[uBl][uBin]));

            graph_Bl[uBl].SetPoint(uBin,Ck_pL->GetBinCenter(0,uBin),Bl_Mean[uBl][uBin]);
            graph_Bl[uBl].SetPointError(uBin,0,sqrt(Bl_Stdv[uBl][uBin]-Bl_Mean[uBl][uBin]*Bl_Mean[uBl][uBin]));

            graph_FitMAX[uBl].SetPoint(uBin,Ck_pL->GetBinCenter(0,uBin),(Fit_Max[uBl][uBin]+Fit_Min[uBl][uBin])*0.5);
            graph_FitMAX[uBl].SetPointError(uBin,0,(Fit_Max[uBl][uBin]-Fit_Min[uBl][uBin])*0.5);

            graph_BlMAX[uBl].SetPoint(uBin,Ck_pL->GetBinCenter(0,uBin),(Bl_Max[uBl][uBin]+Bl_Min[uBl][uBin])*0.5);
            graph_BlMAX[uBl].SetPointError(uBin,0,(Bl_Max[uBl][uBin]-Bl_Min[uBl][uBin])*0.5);

            graph_RatioMAX[uBl].SetPoint(uBin,Ck_pL->GetBinCenter(0,uBin),(Ratio_Max[uBl][uBin]+Ratio_Min[uBl][uBin])*0.5);
            graph_RatioMAX[uBl].SetPointError(uBin,0,(Ratio_Max[uBl][uBin]-Ratio_Min[uBl][uBin])*0.5);
        }
        graph_Fit[uBl].Write();
        graph_Bl[uBl].Write();
        graph_FitMAX[uBl].Write();
        graph_BlMAX[uBl].Write();
        graph_RatioMAX[uBl].Write();
        hnSigma[uBl]->Write();
        hRadius[uBl]->Write();
        hCS[uBl]->Write();
    }


    for(unsigned uBl=0; uBl<=NumBlTypes; uBl++){
        delete [] Fit_Mean[uBl];
        delete [] Fit_Stdv[uBl];
        delete [] Bl_Mean[uBl];
        delete [] Bl_Stdv[uBl];
        delete [] Fit_Min[uBl];
        delete [] Fit_Max[uBl];
        delete [] Bl_Min[uBl];
        delete [] Bl_Max[uBl];
        delete [] Ratio_Min[uBl];
        delete [] Ratio_Max[uBl];

        delete hnSigma[uBl];
        delete hRadius[uBl];
        delete hCS[uBl];
    }

    delete [] Fit_Mean;
    delete [] Fit_Stdv;
    delete [] Bl_Mean;
    delete [] Bl_Stdv;
    delete [] Fit_Min;
    delete [] Fit_Max;
    delete [] Bl_Min;
    delete [] Bl_Max;
    delete [] Ratio_Min;
    delete [] Ratio_Max;


    delete [] hnSigma;
    delete [] hRadius;
    delete [] hCS;

    delete Ck_pSigma0;
    delete Ck_pL;
    delete Ck_pXim;
    delete OutputFile;

}
void pL_SchlechingPlots(const TString InputFileName){

    const bool DataOnly = false;

    const unsigned NumBlTypes = 5;
    enum BLTYPE { pol0,pol1,pol2,pol3,gaus };

    //bool blType[NumBlTypes];
    //blType[pol0] = true;
    //blType[pol1] = true;
    //blType[pol2] = true;
    //blType[pol3] = true;
    //blType[gaus] = true;

    double ResidualSourceSizeSigma0=1.25;
    double ResidualSourceSizeXi=1.0;

    //TRandom3 rangen(SEED);

    TString OutputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/pL_AllHaide_May2020/";
    //TString OutFileName = "OutputFile.root";
    //TFile* OutputFile = new TFile(OutputFolder+OutFileName,"recreate");

    double* MomBins_pL = NULL;
    double* FitRegion_pL = NULL;
    const unsigned NumBinsCk = 45; // for dlm_ck
    const double BinWidthCk = 12;
    const double MaxBinValCk = double(NumBinsCk)*BinWidthCk;
    unsigned NumMomBins_pL = 30;// for CATS
    const double MaxBinValCats = double(NumMomBins_pL)*BinWidthCk;
    unsigned uDummy;

    float** Fit_Min = new float* [NumBlTypes+1];
    float** Fit_Max = new float* [NumBlTypes+1];
    unsigned** Fit_Entries = new unsigned* [NumBlTypes+1];
    float** Bl_Min = new float* [NumBlTypes+1];
    float** Bl_Max = new float* [NumBlTypes+1];

    float** nSigma_Min = new float* [NumBlTypes+1];
    float** nSigma_Max = new float* [NumBlTypes+1];

    float* nSigma_Best = new float [NumBlTypes+1];
    float* nSigma_Worst = new float [NumBlTypes+1];

    for(unsigned uBl=0; uBl<=NumBlTypes; uBl++){
        Fit_Min[uBl] = new float [NumBinsCk];
        Fit_Max[uBl] = new float [NumBinsCk];
        Bl_Min[uBl] = new float [NumBinsCk];
        Bl_Max[uBl] = new float [NumBinsCk];
        nSigma_Min[uBl] = new float [NumBinsCk];
        nSigma_Max[uBl] = new float [NumBinsCk];
        Fit_Entries[uBl] = new unsigned [NumBinsCk];
        for(unsigned uBin=0; uBin<NumBinsCk; uBin++){
            Fit_Min[uBl][uBin] = 10000;
            Fit_Max[uBl][uBin] = 0;
            Bl_Min[uBl][uBin] = 10000;
            Bl_Max[uBl][uBin] = 0;
            nSigma_Min[uBl][uBin] = 10000;
            nSigma_Max[uBl][uBin] = 0;
            Fit_Entries[uBl][uBin] = 0;
        }
        nSigma_Best[uBl] = 10000;
        nSigma_Worst[uBl] = 0;
    }

    CATS AB_Dummy;
    AB_Dummy.SetMomBins(NumMomBins_pL,0,MaxBinValCats);
    DLM_Ck* Ck_Dummy;
    Ck_Dummy = new DLM_Ck(0,0,AB_Dummy,NumBinsCk,0,MaxBinValCk);

    TString DataSample = "pp13TeV_HM_Dec19";
    DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder("/mnt/Ubuntu_Data/CernBox/Sync/CatsFiles");
    AnalysisObject.SetUpBinning_pL(DataSample,uDummy,MomBins_pL,FitRegion_pL,0,0);


        DLM_Ck* Ck_pSigma0 = new DLM_Ck(1,0,NumBinsCk,0,MaxBinValCk,Lednicky_gauss_Sigma0);
        Ck_pSigma0->SetSourcePar(0,ResidualSourceSizeSigma0);

        TH2F* hResolution_pL = AnalysisObject.GetResolutionMatrix(DataSample,"pLambda");
        TH2F* hResidual_pL_pSigma0 = AnalysisObject.GetResidualMatrix("pLambda","pSigma0");
        TH2F* hResidual_pL_pXim = AnalysisObject.GetResidualMatrix("pLambda","pXim");

        double lam_pL[5];
        double lam_pXim[5];

        AnalysisObject.SetUpLambdaPars_pL(DataSample,0,0,lam_pL);
        AnalysisObject.SetUpLambdaPars_pXim(DataSample,0,0,lam_pXim);

        TString SourceDescription = "McGauss_ResoTM";

        CATS AB_pL;
        DLM_Ck* Ck_pL;
        AB_pL.SetMomBins(NumMomBins_pL,0,MaxBinValCats);
        AnalysisObject.SetUpCats_pL(AB_pL,"Chiral_Coupled_SPD",SourceDescription,0,202);//NLO_Coupled_S

        double CUSP_WEIGHT = 0.33;//0.54
        AB_pL.SetChannelWeight(7,1./4.*CUSP_WEIGHT);//1S0 SN(s) -> LN(s)
        AB_pL.SetChannelWeight(8,3./4.*CUSP_WEIGHT);//3S1 SN(s) -> LN(s)
        AB_pL.SetChannelWeight(10,3./4.*CUSP_WEIGHT);//3S1 SN(d) -> LN(s)
        AB_pL.SetChannelWeight(13,3./20.*CUSP_WEIGHT);//3D1 SN(d) -> LN(d)
        AB_pL.SetChannelWeight(15,3./20.*CUSP_WEIGHT);//3D1 SN(s) -> LN(d)

        AB_pL.SetAnaSource(0,1.4);
        if(SourceDescription.Contains("Mc")){
            AB_pL.SetAnaSource(0,1.10);//c.a. 10% smaller compared to p-p due to the mT scaling
            AB_pL.SetAnaSource(1,2.0);
        }
        AB_pL.SetNotifications(CATS::nError);
        AB_pL.KillTheCat();
        Ck_pL = new DLM_Ck(AB_pL.GetNumSourcePars(),0,AB_pL,NumBinsCk,0,MaxBinValCk);
        Ck_pL->SetSourcePar(0,AB_pL.GetAnaSourcePar(0));
        if(SourceDescription.Contains("Mc")){
            Ck_pL->SetSourcePar(0,AB_pL.GetAnaSourcePar(0));
            Ck_pL->SetSourcePar(1,AB_pL.GetAnaSourcePar(1));
        }
        Ck_pL->SetCutOff(340,600);

        CATS AB_pXim;
        //same binning as pL, as we only use pXim as feed-down
        AB_pXim.SetMomBins(NumMomBins_pL,0,MaxBinValCats);
        AnalysisObject.SetUpCats_pXim(AB_pXim,"pXim_HALQCD1","Gauss");
        AB_pXim.SetAnaSource(0,ResidualSourceSizeXi);
        AB_pXim.SetNotifications(CATS::nError);
        AB_pXim.KillTheCat();
        DLM_Ck* Ck_pXim = new DLM_Ck(AB_pXim.GetNumSourcePars(),0,AB_pXim,NumBinsCk,0,MaxBinValCk);
        Ck_pL->Update();
        Ck_pSigma0->Update();
        Ck_pXim->Update();


    unsigned Progress = 0;

    Float_t SourceSize,SourceAlpha,CuspStrength,CkConv,Norm,GaussAmpl,GaussMu,GaussSig,p1,p2,p3;
    Float_t fPOT_VAR,fWhichBl,fWhichFemtoRange,fWhichFitRange,fWhichProtonVar,fWhichLambdaVar,fWhichData,pval;
    unsigned WhichBl,WhichFemtoRange,WhichFitRange,WhichData;

    TFile* InputFile = new TFile(InputFileName,"read");
    TNtuple* ntResult = (TNtuple*)InputFile->Get("ntResult");

    ntResult->SetBranchAddress("SourceSize",&SourceSize);
    ntResult->SetBranchAddress("SourceAlpha",&SourceAlpha);
    ntResult->SetBranchAddress("CuspStrength",&CuspStrength);
    ntResult->SetBranchAddress("CkConv",&CkConv);
    ntResult->SetBranchAddress("Norm",&Norm);
    ntResult->SetBranchAddress("GaussAmpl",&GaussAmpl);
    ntResult->SetBranchAddress("GaussMu",&GaussMu);
    ntResult->SetBranchAddress("GaussSig",&GaussSig);
    ntResult->SetBranchAddress("p1",&p1);
    ntResult->SetBranchAddress("p2",&p2);
    ntResult->SetBranchAddress("p3",&p3);
    ntResult->SetBranchAddress("POT_VAR",&fPOT_VAR);
    ntResult->SetBranchAddress("WhichBl",&fWhichBl);
    ntResult->SetBranchAddress("WhichFemtoRange",&fWhichFemtoRange);
    ntResult->SetBranchAddress("WhichFitRange",&fWhichFitRange);
    ntResult->SetBranchAddress("WhichProtonVar",&fWhichProtonVar);
    ntResult->SetBranchAddress("WhichLambdaVar",&fWhichLambdaVar);
    ntResult->SetBranchAddress("WhichData",&fWhichData);
    ntResult->SetBranchAddress("pval",&pval);

    unsigned NumEntries = ntResult->GetEntries();

    TH1F* hData_pL_Stat;


    for(unsigned uEntry=0; uEntry<NumEntries; uEntry++){

        ntResult->GetEntry(uEntry);

        const unsigned WhichBl = TMath::Nint(fWhichBl);
        const int POT_VAR = TMath::Nint(fPOT_VAR);







//if(WhichBl==3) continue;

//printf("WhichBl=%u (%f)\n",WhichBl,fWhichBl);
        printf("\r\033[K Progress=%u",Progress);
        cout << flush;

        unsigned WhichFemtoRange = TMath::Nint(fWhichFemtoRange);
        unsigned WhichFitRange = TMath::Nint(fWhichFitRange);
        unsigned WhichProtonVar = TMath::Nint(fWhichProtonVar);
        unsigned WhichLambdaVar = TMath::Nint(fWhichLambdaVar);

        AnalysisObject.SetUpBinning_pL(DataSample,uDummy,MomBins_pL,FitRegion_pL,WhichFemtoRange,WhichFitRange);
        AnalysisObject.SetUpLambdaPars_pL(DataSample,WhichProtonVar,WhichLambdaVar,lam_pL);

        //the cusp strength
        //0 is 33%
        //1 is 27% (20% lower, kind of compatible with experiment)
        //2 is 40% (20% larger)
        //3 is free fit
        CUSP_WEIGHT = CuspStrength;
        //double SourceRad = SourceSize;

        DLM_CkDecomposition CkDec_pL("pLambda",4,*Ck_pL,hResolution_pL);
        DLM_CkDecomposition CkDec_pSigma0("pSigma0",0,*Ck_pSigma0,NULL);
        DLM_CkDecomposition CkDec_pXim("pXim",2,*Ck_pXim,NULL);

        //CkDec_pL.AddContribution(0,lam_pL[1],DLM_CkDecomposition::cFeedDown,&CkDec_pSigma0,hResidual_pL_pSigma0);
        //CkDec_pL.AddContribution(1,lam_pL[2],DLM_CkDecomposition::cFeedDown,&CkDec_pXim,hResidual_pL_pXim);
        CkDec_pL.AddContribution(0,lam_pL[1],DLM_CkDecomposition::cFeedDown);
        CkDec_pL.AddContribution(1,lam_pL[2],DLM_CkDecomposition::cFeedDown);
        CkDec_pL.AddContribution(2,lam_pL[3],DLM_CkDecomposition::cFeedDown);
        CkDec_pL.AddContribution(3,lam_pL[4],DLM_CkDecomposition::cFake);//0.03

        //for Xim we simplify a bit and take ALL feed-down as flat
        CkDec_pXim.AddContribution(0,lam_pXim[1]+lam_pXim[2]+lam_pXim[3],DLM_CkDecomposition::cFeedDown);
        CkDec_pXim.AddContribution(1,lam_pXim[4],DLM_CkDecomposition::cFake);

        CkDec_pL.Update();

        unsigned WhichData = TMath::Nint(fWhichData);
        TString DataVar = TString::Format("_%i",WhichData);
        //TH1F* hData_pL = AnalysisObject.GetAliceExpCorrFun(DataSample,"pLambda","_0",1,false,-1);
        TH1F* hData_pL = AnalysisObject.GetAliceExpCorrFun(DataSample,"pLambda",DataVar,2,false,-1);

        pLambda_1_Dec = &CkDec_pL;
        pLambda_1_Cat = &AB_pL;

        //parameters:
        //[0] = source size
        //[1] = source alpha
        //[2] = cusp strength
        //[3] = Ck convergence point
        //baseline of the type: []*(1+[]*exp([mean],[sigma]))*(1+[]*x+[]*x*x+[]*x*x*x)
        //the par in-front of x^3 is defined such as to represent the k* value at which the pol3 has an extrema
        //[4] = norm
        //[5],[6],[7] GAUSS
        //[8],[9],[10] POL
        TF1* fit_pL = new TF1("fit_pL",dimi_pL_Schleching_Fitter,FitRegion_pL[0],FitRegion_pL[3],11);

        fit_pL->FixParameter(0,SourceSize);
        fit_pL->FixParameter(1,SourceAlpha);
        fit_pL->FixParameter(2,CuspStrength);
        fit_pL->FixParameter(3,CkConv);
        fit_pL->FixParameter(4,Norm);
        fit_pL->FixParameter(5,GaussAmpl);
        fit_pL->FixParameter(6,GaussMu);
        fit_pL->FixParameter(7,GaussSig);
        fit_pL->FixParameter(8,p1);
        fit_pL->FixParameter(9,p2);
        fit_pL->FixParameter(10,p3);

        hData_pL->Fit(fit_pL,"Q, S, N, R, M");
        fit_pL->SetNpx(1024);

        TF1* fit_pL_Bl = new TF1("fit_pL_Bl",dimi_pL_Schleching_FitterBl,FitRegion_pL[0],FitRegion_pL[3],7);
        for(unsigned uPar=0; uPar<7; uPar++) fit_pL_Bl->FixParameter(uPar,fit_pL->GetParameter(4+uPar));

        if(nSigma_Best[WhichBl]>fabs(sqrt(2)*TMath::ErfcInverse(pval))&&fabs(sqrt(2)*TMath::ErfcInverse(pval))>1e-6)
            nSigma_Best[WhichBl] = fabs(sqrt(2)*TMath::ErfcInverse(pval));
        if(nSigma_Worst[WhichBl]<fabs(sqrt(2)*TMath::ErfcInverse(pval)))
            nSigma_Worst[WhichBl] = fabs(sqrt(2)*TMath::ErfcInverse(pval));

        if(nSigma_Best[NumBlTypes]>fabs(sqrt(2)*TMath::ErfcInverse(pval))&&fabs(sqrt(2)*TMath::ErfcInverse(pval))>1e-6)
            nSigma_Best[NumBlTypes] = fabs(sqrt(2)*TMath::ErfcInverse(pval));
        if(nSigma_Worst[NumBlTypes]<fabs(sqrt(2)*TMath::ErfcInverse(pval)))
            nSigma_Worst[NumBlTypes] = fabs(sqrt(2)*TMath::ErfcInverse(pval));

        for(unsigned uBin=0; uBin<NumBinsCk; uBin++){
            if(Ck_pL->GetBinCenter(0,uBin)>FitRegion_pL[3]) continue;

            if(Fit_Min[WhichBl][uBin]>fit_pL->Eval(Ck_pL->GetBinCenter(0,uBin)))
                Fit_Min[WhichBl][uBin] = fit_pL->Eval(Ck_pL->GetBinCenter(0,uBin));
            if(Fit_Max[WhichBl][uBin]<fit_pL->Eval(Ck_pL->GetBinCenter(0,uBin)))
                Fit_Max[WhichBl][uBin] = fit_pL->Eval(Ck_pL->GetBinCenter(0,uBin));
            Fit_Entries[WhichBl][uBin] ++;

            if(Bl_Min[WhichBl][uBin]>fit_pL_Bl->Eval(Ck_pL->GetBinCenter(0,uBin)))
                Bl_Min[WhichBl][uBin] = fit_pL_Bl->Eval(Ck_pL->GetBinCenter(0,uBin));
            if(Bl_Max[WhichBl][uBin]<fit_pL_Bl->Eval(Ck_pL->GetBinCenter(0,uBin)))
                Bl_Max[WhichBl][uBin] = fit_pL_Bl->Eval(Ck_pL->GetBinCenter(0,uBin));

            if(nSigma_Min[WhichBl][uBin]>(fit_pL->Eval(Ck_pL->GetBinCenter(0,uBin))-hData_pL->GetBinContent(uBin+1))/hData_pL->GetBinError(uBin+1))
                nSigma_Min[WhichBl][uBin] = (fit_pL->Eval(Ck_pL->GetBinCenter(0,uBin))-hData_pL->GetBinContent(uBin+1))/hData_pL->GetBinError(uBin+1);
            if(nSigma_Max[WhichBl][uBin]<(fit_pL->Eval(Ck_pL->GetBinCenter(0,uBin))-hData_pL->GetBinContent(uBin+1))/hData_pL->GetBinError(uBin+1))
                nSigma_Max[WhichBl][uBin] = (fit_pL->Eval(Ck_pL->GetBinCenter(0,uBin))-hData_pL->GetBinContent(uBin+1))/hData_pL->GetBinError(uBin+1);


            if(Fit_Min[NumBlTypes][uBin]>fit_pL->Eval(Ck_pL->GetBinCenter(0,uBin)))
                Fit_Min[NumBlTypes][uBin] = fit_pL->Eval(Ck_pL->GetBinCenter(0,uBin));
            if(Fit_Max[NumBlTypes][uBin]<fit_pL->Eval(Ck_pL->GetBinCenter(0,uBin)))
                Fit_Max[NumBlTypes][uBin] = fit_pL->Eval(Ck_pL->GetBinCenter(0,uBin));
            Fit_Entries[NumBlTypes][uBin] ++;

            if(Bl_Min[NumBlTypes][uBin]>fit_pL_Bl->Eval(Ck_pL->GetBinCenter(0,uBin)))
                Bl_Min[NumBlTypes][uBin] = fit_pL_Bl->Eval(Ck_pL->GetBinCenter(0,uBin));
            if(Bl_Max[NumBlTypes][uBin]<fit_pL_Bl->Eval(Ck_pL->GetBinCenter(0,uBin)))
                Bl_Max[NumBlTypes][uBin] = fit_pL_Bl->Eval(Ck_pL->GetBinCenter(0,uBin));
            if(nSigma_Min[NumBlTypes][uBin]>(fit_pL->Eval(Ck_pL->GetBinCenter(0,uBin))-hData_pL->GetBinContent(uBin+1))/hData_pL->GetBinError(uBin+1))
                nSigma_Min[NumBlTypes][uBin] = (fit_pL->Eval(Ck_pL->GetBinCenter(0,uBin))-hData_pL->GetBinContent(uBin+1))/hData_pL->GetBinError(uBin+1);
            if(nSigma_Max[NumBlTypes][uBin]<(fit_pL->Eval(Ck_pL->GetBinCenter(0,uBin))-hData_pL->GetBinContent(uBin+1))/hData_pL->GetBinError(uBin+1))
                nSigma_Max[NumBlTypes][uBin] = (fit_pL->Eval(Ck_pL->GetBinCenter(0,uBin))-hData_pL->GetBinContent(uBin+1))/hData_pL->GetBinError(uBin+1);

        }

        delete fit_pL;
        delete fit_pL_Bl;

        delete Ck_pSigma0;
        delete Ck_pL;
        delete Ck_pXim;

        Progress++;

    }
printf("\n");

        TGraphErrors* graph_FitMAX = new TGraphErrors[NumBlTypes+1];
        TGraphErrors* graph_BlMAX = new TGraphErrors[NumBlTypes+1];
        TGraphErrors* graph_nSigmaMAX = new TGraphErrors[NumBlTypes+1];

        for(unsigned uBl=0; uBl<=NumBlTypes; uBl++){
            if(uBl<NumBlTypes){
                graph_FitMAX[uBl].SetName(TString::Format("graph_FitMAX_%u",uBl));
                graph_BlMAX[uBl].SetName(TString::Format("graph_BlMAX_%u",uBl));
                graph_nSigmaMAX[uBl].SetName(TString::Format("graph_nSigmaMAX_%u",uBl));
            }
            else{
                graph_FitMAX[uBl].SetName(TString::Format("graph_FitMAX"));
                graph_BlMAX[uBl].SetName(TString::Format("graph_BlMAX"));
                graph_nSigmaMAX[uBl].SetName(TString::Format("graph_nSigmaMAX"));
            }
            for(unsigned uBin=0; uBin<NumBinsCk; uBin++){
                graph_FitMAX[uBl].SetPoint(uBin,Ck_Dummy->GetBinCenter(0,uBin),(Fit_Max[uBl][uBin]+Fit_Min[uBl][uBin])*0.5);
                graph_FitMAX[uBl].SetPointError(uBin,0,(Fit_Max[uBl][uBin]-Fit_Min[uBl][uBin])*0.5);
//printf("%u %f %f\n",uBin,Ck_pL->GetBinCenter(0,uBin),(Fit_Max[uBl][uBin]+Fit_Min[uBl][uBin])*0.5);
                graph_BlMAX[uBl].SetPoint(uBin,Ck_Dummy->GetBinCenter(0,uBin),(Bl_Max[uBl][uBin]+Bl_Min[uBl][uBin])*0.5);
                graph_BlMAX[uBl].SetPointError(uBin,0,(Bl_Max[uBl][uBin]-Bl_Min[uBl][uBin])*0.5);

                graph_nSigmaMAX[uBl].SetPoint(uBin,Ck_Dummy->GetBinCenter(0,uBin),(nSigma_Max[uBl][uBin]+nSigma_Min[uBl][uBin])*0.5);
                graph_nSigmaMAX[uBl].SetPointError(uBin,0,(nSigma_Max[uBl][uBin]-nSigma_Min[uBl][uBin])*0.5);

//printf("%u %f %f\n",uBin,Ck_pL->GetBinCenter(0,uBin),(nSigma_Max[uBl][uBin]+nSigma_Min[uBl][uBin])*0.5);
            }

            //graph_FitMAX[uBl].Write();
            //graph_BlMAX[uBl].Write();
            //graph_nSigmaMAX[uBl].Write();
        }


///////////////////////////////////////////////

    TFile* fPlot = new TFile(OutputFolder+"fPlot_"+".root","recreate");

    hData_pL_Stat = AnalysisObject.GetAliceExpCorrFun(DataSample,"pLambda","_0",2,false,-1);

    gStyle->SetCanvasPreferGL(1);
    SetStyle();
    fPlot->cd();
hData_pL_Stat->Write();
graph_FitMAX[0].Write();
graph_FitMAX[1].Write();
graph_FitMAX[2].Write();
graph_FitMAX[3].Write();
graph_FitMAX[4].Write();
graph_FitMAX[5].Write();
    const float right = 0.025;
    const float top = 0.025;

    graph_FitMAX[NumBlTypes].SetFillColorAlpha(kRed+1,0.67);
    graph_FitMAX[NumBlTypes].SetLineColor(kRed+1);
    graph_FitMAX[NumBlTypes].SetLineWidth(5);

    graph_BlMAX[NumBlTypes].SetFillColorAlpha(kCyan+1,0.5);
    graph_BlMAX[NumBlTypes].SetLineColor(kCyan+1);
    graph_BlMAX[NumBlTypes].SetLineWidth(5);

    DLM_SubPads DlmPad(720,720);
    DlmPad.AddSubPad(0,1,0.33,1);
    DlmPad.AddSubPad(0,1,0,0.33);
    DlmPad.SetMargin(0,0.12,0.02,0.0,0.02);
    DlmPad.SetMargin(1,0.12,0.02,0.09,0.0);
    DlmPad.cd(0);

    hData_pL_Stat->SetTitle("; #it{k*} (MeV/#it{c}); #it{C}(#it{k*})");
    hData_pL_Stat->GetXaxis()->SetRangeUser(0, 336);
    hData_pL_Stat->GetXaxis()->SetNdivisions(505);
    hData_pL_Stat->GetYaxis()->SetRangeUser(0.9, 2.0);
    hData_pL_Stat->SetFillColor(kGray+1);
    SetStyleHisto2(hData_pL_Stat,2,0);
    //hData_pL_Stat->GetYaxis()->SetTitleOffset(1.0);
    hData_pL_Stat->Draw();

    TFile* DataSystFile = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample10HM/Systematics_pL.root");
    TH1F* DataSystHisto = (TH1F*)DataSystFile->Get("SystErrRel");
    TGraphErrors *Tgraph_syserror = DrawSystematicError_FAST(hData_pL_Stat, DataSystHisto, NULL, 3);
    Tgraph_syserror->SetLineColor(kWhite);

    //baselineLL->Draw("same");

    if(!DataOnly){
    graph_BlMAX[NumBlTypes].Draw("3 same");
    graph_FitMAX[NumBlTypes].Draw("3 same");
    }


    hData_pL_Stat->Draw("same");

    Tgraph_syserror->SetFillColorAlpha(kBlack, 0.4);
    Tgraph_syserror->Draw("2 same");
    //hData->Draw("pe same");

    TString LegendSource_line1 = "Gaussian core + resonances";
    TString LegendSource_line2 = TString::Format("n#sigma#in(%.1f,%.1f)",nSigma_Best[NumBlTypes],nSigma_Worst[NumBlTypes]);

    unsigned NumRows=4;
    TLegend *legend = new TLegend(0.39,0.73-0.055*NumRows,0.73,0.73);//lbrt
    legend->SetBorderSize(0);
    legend->SetTextFont(42);
    legend->SetTextSize(gStyle->GetTextSize()*0.90);
    TH1F* hCk_Fake;
    hCk_Fake = (TH1F*)hData_pL_Stat->Clone("hCk_Fake");
    hCk_Fake->SetName("hCk_Fake");
    hCk_Fake->SetLineColor(hCk_Fake->GetFillColor());

    legend->AddEntry(hCk_Fake, "p#minus #Lambda #oplus #bar{p}#minus #bar{#Lambda} pairs", "fpe");
    if(!DataOnly){
    legend->AddEntry(&graph_FitMAX[NumBlTypes],"Femtoscopic fit (NLO)","l");
    legend->AddEntry(&graph_BlMAX[NumBlTypes],"Baseline","l");
    }
    legend->Draw("same");
    TLatex BeamText;
    BeamText.SetTextSize(gStyle->GetTextSize()*0.90);
    BeamText.SetNDC(kTRUE);
    BeamText.DrawLatex(0.40, 0.915, "ALICE work in progress");
    BeamText.DrawLatex(0.40, 0.860, "high-mult. (0#minus 0.17% INEL>0) pp #sqrt{#it{s}} = 13 TeV");

    TLatex BeamTextSource;
    BeamTextSource.SetTextSize(gStyle->GetTextSize()*0.90);
    BeamTextSource.SetNDC(kTRUE);
    BeamTextSource.DrawLatex(0.40, 0.805, LegendSource_line1);
    BeamTextSource.DrawLatex(0.40, 0.750, LegendSource_line2);

//INLET -------------------------------------------------------------------------------------------------------------------

    TH1F* DataHisto_Inlet = (TH1F*)hData_pL_Stat->Clone("DataHisto_Inlet");
    DataHisto_Inlet->SetMarkerSize(hData_pL_Stat->GetMarkerSize()*0.67);
    DataHisto_Inlet->SetLineWidth(hData_pL_Stat->GetLineWidth()*0.67);
    DataHisto_Inlet->GetXaxis()->SetTitleSize(hData_pL_Stat->GetXaxis()->GetTitleSize()*1.75);
    DataHisto_Inlet->GetXaxis()->SetLabelSize(hData_pL_Stat->GetXaxis()->GetLabelSize()*1.75);
    DataHisto_Inlet->GetXaxis()->SetRangeUser(100, 336);
    DataHisto_Inlet->GetXaxis()->SetNdivisions(505);

    DataHisto_Inlet->GetYaxis()->SetTitleSize(hData_pL_Stat->GetYaxis()->GetTitleSize()*1.75);
    DataHisto_Inlet->GetYaxis()->SetLabelSize(hData_pL_Stat->GetYaxis()->GetLabelSize()*1.75);
    DataHisto_Inlet->GetYaxis()->SetTitleOffset(hData_pL_Stat->GetYaxis()->GetTitleOffset()*0.67);
    DataHisto_Inlet->GetYaxis()->SetRangeUser(0.98, 1.05);

    TGraph* grFemto_Inlet = (TGraph*)graph_FitMAX[NumBlTypes].Clone("grFemto_Inlet");
    grFemto_Inlet->SetLineWidth(graph_FitMAX[NumBlTypes].GetLineWidth()*0.67);


    const double fXMinInlet=0.30;
    const double fYMinInlet=0.12;
    const double fXMaxInlet=0.95;
    const double fYMaxInlet=0.50;
    TPad *inset_pad = new TPad("insert", "insertPad", fXMinInlet, fYMinInlet,
                             fXMaxInlet, fYMaxInlet);
    inset_pad->SetTopMargin(0.01);
    inset_pad->SetRightMargin(0.05);
    inset_pad->SetBottomMargin(0.28);
    inset_pad->SetLeftMargin(0.28);
    inset_pad->SetFillStyle(4000);
    if(!DataOnly) inset_pad->Draw();
    inset_pad->cd();
    DataHisto_Inlet->Draw();
    if(!DataOnly){
    graph_BlMAX[NumBlTypes].Draw("3 same");
    grFemto_Inlet->Draw("3 same");
    }
    DataHisto_Inlet->Draw("same");
    Tgraph_syserror->Draw("2 same");

    DlmPad.cd(1);
    TH1F* hAxis = new TH1F("hAxis", "hAxis", 28, 0, 336);
    hAxis->SetStats(false);
    if(!DataOnly) hAxis->SetTitle("; #it{k*} (MeV/#it{c}); #it{n_{#sigma}}");
    else hAxis->SetTitle("; #it{k*} (MeV/#it{c}); #it{C}(#it{k*})");
    hAxis->GetXaxis()->SetRangeUser(0, 336);
    if(!DataOnly) hAxis->GetYaxis()->SetRangeUser(-6, 7);
    else hAxis->GetYaxis()->SetRangeUser(0.98, 1.05);
    //hData->SetTitle("; #it{k*} (MeV/#it{c}); #it{C}(#it{k*})");
    //hData->GetXaxis()->SetRangeUser(0, 320);
    hAxis->GetXaxis()->SetNdivisions(505);
    //hData->GetYaxis()->SetRangeUser(0.85, 2.3);
    //hData->SetFillColor(fFillColors[0]);
    SetStyleHisto2(hAxis,2,0,2);
    //hData->GetYaxis()->SetTitleOffset(1.0);
    hAxis->Draw("");

    graph_nSigmaMAX[NumBlTypes].SetFillColorAlpha(kRed+1,0.67);
    if(!DataOnly) graph_nSigmaMAX[NumBlTypes].Draw("3 same");
    else{
    hData_pL_Stat->Draw("same");
    Tgraph_syserror->Draw("2 same");
    }

    DlmPad.GetCanvas()->SaveAs(OutputFolder+"DlmPad.pdf");
/*
    delete hVarDeviationAV18;
    delete hVarDevRatioAV18;
    delete legend;
    delete hCk_Fake;
    delete Tgraph_syserror;
    delete grFemto_AV18;
    delete grOuterBl_AV18;
    //delete Can_CF_pL;
    delete hChi2NdfAV18;
    delete hNsigmaAV18;
    delete fPlot;
    delete DataFile;
    delete SystFile;
    delete ntFile;
*/













/////////////////////////////


    delete [] Fit_Min;
    delete [] Fit_Max;
    delete [] Bl_Min;
    delete [] Bl_Max;
    delete [] nSigma_Min;
    delete [] nSigma_Max;


    delete Ck_Dummy;
}

//the radius is fixed, and we use the latest Core+Reso
/// if the SEED<100 we take the default iteration first
//i.e. to avoid the default iteration, use SEED>=100
//Data type digits: 321
//1: corrected for SB (0/1 for no/yes)
//2: unfolded (0/1 for no/yes)
//3: corrected for flat feed and misid (0/1 for no/yes) -> done on the fly at the PLOTTING, no actual corrected datasets are available
//   3 is not used by this function at all
void pL_SystematicsMay2020(unsigned SEED, unsigned BASELINE_VAR, int POT_VAR, int Sigma0_Feed, int Data_Type,
                           bool DataSyst, bool FitSyst, bool Bootstrap, unsigned NumIter,
                           const char* CatsFileFolder, const char* OutputFolder){

    //if we go beyond the 1 hour 50 minutes mark, we stop
    //safety for the batch farm
    const double TIME_LIMIT = 110;
    //do extra variations
    const bool ExtendedSyst = true;
    DLM_Timer TIMER_SYST;

    //pol(0/1/2/3)s: pol(0/1/2/3) with a small fit range
    //dpol(2/3/4)s: pol(2/3/4) that is flat at zero, with a small fit range (up to c.a. 300 MeV)
    //dpol(2/3/4)e: pol(2/3/4) that is flat at zero, with an extended fit range (up to c.a. 450 MeV)
    //the latter case has an additional fit parameter, that is the convergence point of the correlation
///The current choice:
/// as cross check: pol0s,pol1s,pol2s,dpol2s,pol3e
/// as systematics: dpol3s,dpol3e,dpol4e
/// for fun: spl1 (splines at 0-48-96-192-288-432)
    enum BLTYPE { pol0s,pol1s,pol2s,pol3s,dpol2s,dpol3s,dpol4s,pol2e,pol3e,dpol2e,dpol3e,dpol4e,spl1 };

    // new syst: pLambda(1.55):   low 0.853396 mean 0.877717 up 0.90239
    //the values below is for the updated epos version of the source
    const float SourceSizeErr = (1.05676-0.980842)*0.5;
    const float SourceSize = (1.05676+0.980842)*0.5+2.*SourceSizeErr;

    float SourceAlpha = 2.0;

    double ResidualSourceSizeSigma0=1.25;
    double ResidualSourceSizeXi=1.0;

    TRandom3 rangen(SEED);

    //TString DataSample = "pp13TeV_HM_Dec19";
    //TString DataSample = "pp13TeV_HM_RotPhiDec19";
    TString DataSample;
    //N.B. for this data we do not have the unfolded version
    if(Data_Type%10==0) DataSample = "pp13TeV_HM_Dec19";
    //corrected for SB, folded
    else if(Data_Type%10==1&&Data_Type/10==0) DataSample = "pp13TeV_HM_DimiJun20";
    else DataSample = "pp13TeV_HM_DimiJul20";

    //TString SourceDescription = "Gauss";
    TString SourceDescription = "McGauss_ResoTM";
    //TString SourceDescription = "McLevy_ResoTM";

    TString OutFileName = TString::Format("Output_%s_POT%i_BL%i_SIG%i_%u.root",DataSample.Data(),POT_VAR,BASELINE_VAR,Sigma0_Feed,SEED);
    TFile* OutputFile = new TFile(TString(OutputFolder)+OutFileName,"recreate");

    double* MomBins_pL = NULL;
    double* FitRegion_pL = NULL;
    const unsigned NumBinsCk = 45; // for dlm_ck
    const double BinWidthCk = 12;
    const double MaxBinValCk = double(NumBinsCk)*BinWidthCk;
    unsigned NumMomBins_pL = 30;// for CATS
    const double MaxBinValCats = double(NumMomBins_pL)*BinWidthCk;
    unsigned uDummy;

    TGraphErrors gData;
    gData.SetName("gData");
    TGraphErrors fitData;
    fitData.SetName("fitData");
    TGraphErrors fitFemto;
    fitFemto.SetName("fitFemto");
    TGraphErrors fitGenuine;
    fitGenuine.SetName("fitGenuine");
    TGraphErrors fitSignal_pL;
    fitSignal_pL.SetName("fitSignal_pL");
    TGraphErrors fitSignal_pL_pS0;
    fitSignal_pL_pS0.SetName("fitSignal_pL_pS0");
    TGraphErrors fitSignal_pL_pXim;
    fitSignal_pL_pXim.SetName("fitSignal_pL_pXim");
    TGraphErrors fitSignal_pL_pXi0;
    fitSignal_pL_pXi0.SetName("fitSignal_pL_pXi0");
    TGraphErrors fitBaseline;
    fitBaseline.SetName("fitBaseline");

    float SourceRad;
    float CuspWeight = 0.33;
    float CkConv = 700;
    float CkCutOff = 360;
    float Norm;
    float pp1;
    float pp2;
    float pp3;
    float pp4;
    float kFemtoMin;
    float kFemtoMax;
    float kLongMin;
    float kLongMax;
    float lam_L_genuine;
    float lam_L_Sig0;
    float lam_L_Xim;
    float lam_L_Xi0;
    float lam_L_Flat;
    float lam_S0_genuine;
    float lam_S0_Flat;
    float lam_Xim_genuine;
    float lam_Xim_Flat;
    float pval;
    unsigned WhichData;
    float SBl;//fraction of left sideband
    float SBpur;//the purity of the sample, used for the correction with sidebands
    bool DefaultVariation;
    //int Sigma0_Feed;
    int Xim_Feed;
    float POT_P1;
    float POT_P2;
    float POT_P3;

/// WHAT ABOUT MISID AND SIDEBANDS
//Lambda purity is around 95%, protons is 99.4%
//for protons: we can prob ignore (pi-Lambda or K-Lambda are not strong enough to pop-out in the correlation)
//for Lambdas: sideband? The sideband should cover most of 4.6% lambda par inside, if we parametrize it we can
//include it as a systematic variation

    TTree* plambdaTree = new TTree("plambdaTree","plambdaTree");
    plambdaTree->Branch("gData","TGraphErrors",&gData,32000,0);//
    plambdaTree->Branch("fitData","TGraphErrors",&fitData,32000,0);//
    plambdaTree->Branch("fitFemto","TGraphErrors",&fitFemto,32000,0);//
    plambdaTree->Branch("fitGenuine","TGraphErrors",&fitGenuine,32000,0);//
    plambdaTree->Branch("fitSignal_pL","TGraphErrors",&fitSignal_pL,32000,0);//
    plambdaTree->Branch("fitSignal_pL_pS0","TGraphErrors",&fitSignal_pL_pS0,32000,0);//
    plambdaTree->Branch("fitSignal_pL_pXim","TGraphErrors",&fitSignal_pL_pXim,32000,0);//
    plambdaTree->Branch("fitSignal_pL_pXi0","TGraphErrors",&fitSignal_pL_pXi0,32000,0);//
    plambdaTree->Branch("fitBaseline","TGraphErrors",&fitBaseline,32000,0);//
    plambdaTree->Branch("DataSample","TString",&DataSample,8000,0);//
    plambdaTree->Branch("SourceDescription","TString",&SourceDescription,8000,0);//
    plambdaTree->Branch("SEED", &SEED, "SEED/i");//
    plambdaTree->Branch("BASELINE_VAR", &BASELINE_VAR, "BASELINE_VAR/i");//
    plambdaTree->Branch("POT_VAR", &POT_VAR, "POT_VAR/I");//
    //plambdaTree->Branch("POT_P1", &POT_P1, "POT_P1/I");//
    //plambdaTree->Branch("POT_P2", &POT_P2, "POT_P2/I");//
    //plambdaTree->Branch("POT_P3", &POT_P3, "POT_P3/I");//
    plambdaTree->Branch("Sigma0_Feed", &Sigma0_Feed, "Sigma0_Feed/I");//
    plambdaTree->Branch("Xim_Feed", &Xim_Feed, "Xim_Feed/I");//
    plambdaTree->Branch("DataSyst", &DataSyst, "DataSyst/O");//
    plambdaTree->Branch("FitSyst", &FitSyst, "FitSyst/O");//
    plambdaTree->Branch("Bootstrap", &Bootstrap, "Bootstrap/O");//
    plambdaTree->Branch("SourceSize", &SourceRad, "SourceSize/F");//
    plambdaTree->Branch("SourceAlpha", &SourceAlpha, "SourceAlpha/F");//
    plambdaTree->Branch("CuspWeight", &CuspWeight, "CuspWeight/F");//
    plambdaTree->Branch("CkConv", &CkConv, "CkConv/F");//
    plambdaTree->Branch("CkCutOff", &CkCutOff, "CkCutOff/F");//
    plambdaTree->Branch("Norm", &Norm, "Norm/F");
    plambdaTree->Branch("pp1", &pp1, "pp1/F");
    plambdaTree->Branch("pp2", &pp2, "pp2/F");
    plambdaTree->Branch("pp3", &pp3, "pp3/F");
    plambdaTree->Branch("pp4", &pp4, "pp4/F");
    plambdaTree->Branch("kFemtoMin", &kFemtoMin, "kFemtoMin/F");
    plambdaTree->Branch("kFemtoMax", &kFemtoMax, "kFemtoMax/F");
    plambdaTree->Branch("kLongMin", &kLongMin, "kLongMin/F");
    plambdaTree->Branch("kLongMax", &kLongMax, "kLongMax/F");
    plambdaTree->Branch("lam_L_genuine", &lam_L_genuine, "lam_L_genuine/F");
    plambdaTree->Branch("lam_L_Sig0", &lam_L_Sig0, "lam_L_Sig0/F");
    plambdaTree->Branch("lam_L_Xim", &lam_L_Xim, "lam_L_Xim/F");
    plambdaTree->Branch("lam_L_Xi0", &lam_L_Xi0, "lam_L_Xi0/F");
    plambdaTree->Branch("lam_L_Flat", &lam_L_Flat, "lam_L_Flat/F");
    plambdaTree->Branch("lam_S0_genuine", &lam_S0_genuine, "lam_S0_genuine/F");
    plambdaTree->Branch("lam_S0_Flat", &lam_S0_Flat, "lam_S0_Flat/F");
    plambdaTree->Branch("lam_Xim_genuine", &lam_Xim_genuine, "lam_Xim_genuine/F");
    plambdaTree->Branch("lam_Xim_Flat", &lam_Xim_Flat, "lam_Xim_Flat/F");
    plambdaTree->Branch("WhichData", &WhichData, "WhichData/i");
    plambdaTree->Branch("SBl", &SBl, "SBl/F");
    plambdaTree->Branch("SBpur", &SBpur, "SBpur/F");
    plambdaTree->Branch("pval", &pval, "pval/F");
    plambdaTree->Branch("DefaultVariation", &DefaultVariation, "DefaultVariation/O");

    DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder(CatsFileFolder);
    AnalysisObject.SetUpBinning_pL(DataSample,uDummy,MomBins_pL,FitRegion_pL,0,0);

    TH2F* hResolution_pL = AnalysisObject.GetResolutionMatrix(DataSample,"pLambda");
    TH2F* hResidual_pL_pSigma0 = AnalysisObject.GetResidualMatrix("pLambda","pSigma0");
    TH2F* hResidual_pL_pXim = AnalysisObject.GetResidualMatrix("pLambda","pXim");
    TH2F* hResidual_pL_pXi0 = AnalysisObject.GetResidualMatrix("pLambda","pXi0");

    CATS AB_pL;
    DLM_Ck* Ck_pL;
    AB_pL.SetMomBins(NumMomBins_pL,0,MaxBinValCats);
    AnalysisObject.SetUpCats_pL(AB_pL,"Chiral_Coupled_SPD",SourceDescription,POT_VAR,202);//NLO_Coupled_S

    CuspWeight = 0.33;//0.54
    if(abs(POT_VAR)>1000){
        AB_pL.SetChannelWeight(7,1./4.*CuspWeight);//1S0 SN(s) -> LN(s)
        AB_pL.SetChannelWeight(8,3./4.*CuspWeight);//3S1 SN(s) -> LN(s)
        AB_pL.SetChannelWeight(10,3./4.*CuspWeight);//3S1 SN(d) -> LN(s)
        AB_pL.SetChannelWeight(13,3./20.*CuspWeight);//3D1 SN(d) -> LN(d)
        AB_pL.SetChannelWeight(15,3./20.*CuspWeight);//3D1 SN(s) -> LN(d)
    }

    AB_pL.SetAnaSource(0,1.4);
    if(SourceDescription.Contains("Mc")){
        AB_pL.SetAnaSource(0,1.10);//c.a. 10% smaller compared to p-p due to the mT scaling
        AB_pL.SetAnaSource(1,SourceAlpha);
    }
    AB_pL.SetNotifications(CATS::nError);
    AB_pL.KillTheCat();

    const unsigned NumFineBins = 25;//it works out like that if we want 1.5 MeV step size, as defined above in the fit function
    const double FineValueMin = 264;
    const double FineValueMax = 300;
    const double FineWidth = (FineValueMax-FineValueMin)/double(NumFineBins-1);
    CATS ABF_pL;
    ABF_pL.SetMomBins(NumFineBins,FineValueMin-FineWidth*0.5,FineValueMax+FineWidth*0.5);
    AnalysisObject.SetUpCats_pL(ABF_pL,"Chiral_Coupled_SPD",SourceDescription,POT_VAR,202);//NLO_Coupled_S

    CuspWeight = 0.33;//0.54
    if(abs(POT_VAR)>1000){
        ABF_pL.SetChannelWeight(7,1./4.*CuspWeight);//1S0 SN(s) -> LN(s)
        ABF_pL.SetChannelWeight(8,3./4.*CuspWeight);//3S1 SN(s) -> LN(s)
        ABF_pL.SetChannelWeight(10,3./4.*CuspWeight);//3S1 SN(d) -> LN(s)
        ABF_pL.SetChannelWeight(13,3./20.*CuspWeight);//3D1 SN(d) -> LN(d)
        ABF_pL.SetChannelWeight(15,3./20.*CuspWeight);//3D1 SN(s) -> LN(d)
    }

    ABF_pL.SetAnaSource(0,1.4);
    if(SourceDescription.Contains("Mc")){
        ABF_pL.SetAnaSource(0,1.10);//c.a. 10% smaller compared to p-p due to the mT scaling
        ABF_pL.SetAnaSource(1,SourceAlpha);
    }
    ABF_pL.SetNotifications(CATS::nError);
    ABF_pL.KillTheCat();

    //usmani fit = 1
    const unsigned PotFitPar = POT_VAR==101?1:0;
    Ck_pL = new DLM_Ck(AB_pL.GetNumSourcePars(),PotFitPar,AB_pL,NumBinsCk,0,MaxBinValCk);
    Ck_pL->SetSourcePar(0,AB_pL.GetAnaSourcePar(0));
    if(SourceDescription.Contains("Mc")){
        Ck_pL->SetSourcePar(0,AB_pL.GetAnaSourcePar(0));
        Ck_pL->SetSourcePar(1,AB_pL.GetAnaSourcePar(1));
    }
    Ck_pL->SetCutOff(CkCutOff,CkConv);

    CATS AB_pXim;
    //same binning as pL, as we only use pXim as feed-down
    AB_pXim.SetMomBins(NumMomBins_pL,0,MaxBinValCats);
    AnalysisObject.SetUpCats_pXim(AB_pXim,"pXim_HALQCDPaper2020","Gauss");
    AB_pXim.SetAnaSource(0,ResidualSourceSizeXi);
    AB_pXim.SetNotifications(CATS::nError);
    AB_pXim.KillTheCat();
    DLM_Ck* Ck_pXim = new DLM_Ck(AB_pXim.GetNumSourcePars(),0,AB_pXim,NumBinsCk,0,MaxBinValCk);
    Ck_pXim->SetCutOff(360,700);

    CATS AB_pXi0;
    //same binning as pL, as we only use pXim as feed-down
    AB_pXi0.SetMomBins(NumMomBins_pL,0,MaxBinValCats);
    AnalysisObject.SetUpCats_pXi0(AB_pXi0,"pXim_HALQCDPaper2020","Gauss");
    AB_pXi0.SetAnaSource(0,ResidualSourceSizeXi);
    AB_pXi0.SetNotifications(CATS::nError);
    AB_pXi0.KillTheCat();
    DLM_Ck* Ck_pXi0 = new DLM_Ck(AB_pXi0.GetNumSourcePars(),0,AB_pXi0,NumBinsCk,0,MaxBinValCk);
    Ck_pXi0->SetCutOff(360,700);

    CATS AB_pS0_Chiral;
    //the minus one is to to avoid going above 350 MeV, since we do not have the WF there
    AB_pS0_Chiral.SetMomBins(NumMomBins_pL-1,0,MaxBinValCats-BinWidthCk);
    AnalysisObject.SetUpCats_pS0(AB_pS0_Chiral,"Chiral","Gauss");
    AB_pS0_Chiral.SetAnaSource(0,ResidualSourceSizeSigma0);
    AB_pS0_Chiral.SetNotifications(CATS::nError);
    AB_pS0_Chiral.KillTheCat();

    CATS AB_pS0_ESC16;
    AB_pS0_ESC16.SetMomBins(NumMomBins_pL-1,0,MaxBinValCats-BinWidthCk);
    AnalysisObject.SetUpCats_pS0(AB_pS0_ESC16,"ESC16","Gauss");
    AB_pS0_ESC16.SetAnaSource(0,ResidualSourceSizeSigma0);
    AB_pS0_ESC16.SetNotifications(CATS::nError);
    AB_pS0_ESC16.KillTheCat();

    DLM_Ck* Ck_pS0_Chiral = new DLM_Ck(AB_pS0_Chiral.GetNumSourcePars(),0,AB_pS0_Chiral,NumBinsCk,0,MaxBinValCk);
    Ck_pS0_Chiral->SetSourcePar(0,ResidualSourceSizeSigma0);
    Ck_pS0_Chiral->SetCutOff(330,700);

    DLM_Ck* Ck_pS0_ESC16 = new DLM_Ck(AB_pS0_ESC16.GetNumSourcePars(),0,AB_pS0_ESC16,NumBinsCk,0,MaxBinValCk);
    Ck_pS0_ESC16->SetSourcePar(0,ResidualSourceSizeSigma0);
    Ck_pS0_ESC16->SetCutOff(330,700);

    Ck_pL->Update();
    Ck_pS0_Chiral->Update();
    Ck_pS0_ESC16->Update();
    Ck_pXim->Update();
    Ck_pXi0->Update();

    for(unsigned uIter=0; uIter<NumIter; uIter++){

        //in ms
        long long ExeTime = TIMER_SYST.Stop()/1000.;
        if(TIME_LIMIT<double(ExeTime)/1000./60.) break;

        printf("\r\033[K Progress=%.0f%%",double(uIter)/double(NumIter)*100.);
        cout << flush;

        DefaultVariation = ( (uIter==0&&SEED<100) );

        unsigned WhichFemtoRange = rangen.Integer(3);
        if(DefaultVariation||FitSyst==false) WhichFemtoRange = 0;
        unsigned WhichFitRange;
        if(BASELINE_VAR<pol2e||BASELINE_VAR>spl1) WhichFitRange=50;
        //for the Dimi baseline, fit to larger range
        else WhichFitRange = 51+WhichFemtoRange;

        AnalysisObject.SetUpBinning_pL(DataSample,uDummy,MomBins_pL,FitRegion_pL,WhichFemtoRange,WhichFitRange);
        kFemtoMin = FitRegion_pL[0];
        kFemtoMax = FitRegion_pL[1];//plays no role
        kLongMin = FitRegion_pL[2];//play no role
        kLongMax = FitRegion_pL[3];

        //the cusp strength
        //0 is 33%
        //1 is 27% (20% lower, kind of compatible with experiment)
        //2 is 40% (20% larger)
        //3 is free fit in the range [17%, 50%]
        unsigned WhichCuspStrength = rangen.Integer(ExtendedSyst?5:3);
        if(DefaultVariation||FitSyst==false) WhichCuspStrength = 0;
//0.2,0.267,0.333,0.4,0.467
        switch(WhichCuspStrength){
            case 0 : CuspWeight = 1./3.; break;
            case 1 : CuspWeight = 1./3.*0.8; break;//0.27
            case 2 : CuspWeight = 1./3.*1.2; break;//0.40
            case 3 : CuspWeight = 1./3.*0.6; break;
            case 4 : CuspWeight = 1./3.*1.4; break;
            default : CuspWeight = 1./3.; break;
        }

        //this parameter only plays a role for the extended fit range
        unsigned WhichCkConv = rangen.Integer(3);
        if(BASELINE_VAR<pol2e||BASELINE_VAR>spl1) WhichCkConv=0;
        if(DefaultVariation||FitSyst==false) WhichCkConv=0;
        switch(WhichCkConv){
            case 0 : CkConv = 700; break;
            case 1 : CkConv = 500; break;
            case 2 : CkConv = 1000; break;
            default : CkConv = 700; break;
        }

        unsigned WhichSourceRad = rangen.Integer(ExtendedSyst?5:3);
        if(DefaultVariation||FitSyst==false) {WhichSourceRad=0;}
        //if(DefaultVariation) SourceRad = SourceSize;
        //else SourceRad = rangen.Gaus(SourceSize,SourceSizeErr);
        switch(WhichSourceRad){
            case 0 : SourceRad = SourceSize; break;
            case 1 : SourceRad = SourceSize-SourceSizeErr; break;
            case 2 : SourceRad = SourceSize+SourceSizeErr; break;
            case 3 : SourceRad = SourceSize+2.*SourceSizeErr; break;
            case 4 : SourceRad = SourceSize-2.*SourceSizeErr; break;
            default : SourceRad = SourceSize; break;
        }

        unsigned WhichSourceAlpha = rangen.Integer(5);
        if(DefaultVariation||SourceDescription.Contains("Gauss")||FitSyst==false) {WhichSourceAlpha=0;}
        switch(WhichSourceAlpha){
            case 0 : SourceAlpha = 2.0; break;
            case 1 : SourceAlpha = 1.9; break;
            case 2 : SourceAlpha = 1.8; break;
            case 3 : SourceAlpha = 1.7; break;
            case 4 : SourceAlpha = 1.6; break;
            default : SourceAlpha = 2.0; break;
        }

        double lam_pL[5];
        double lam_pXim[5];


        //has no effect actually
        //int WhichProtonVar = rangen.Integer(3);
        //if(DefaultVariation||FitSyst==false) WhichProtonVar = 0;
        int WhichProtonVar = 0;
        int WhichLambdaVar = rangen.Integer(ExtendedSyst?5:3);
        if(DefaultVariation||FitSyst==false) WhichLambdaVar = 0;

        //the modification of the amount of Xi secondaries with respect to lambda+sigma0
        if(DefaultVariation||FitSyst==false) WhichLambdaVar += 10*(0+3.*0.);
        else WhichLambdaVar += 10*(0+3.*rangen.Integer(2));

        //has no effect if you run with the data that is SB corrected
        if(DefaultVariation||FitSyst==false||!ExtendedSyst) WhichLambdaVar += 100;//the new purities
        else WhichLambdaVar += 100*(1+rangen.Integer(3));//new purities with variations

        AnalysisObject.SetUpLambdaPars_pL(DataSample,WhichProtonVar,WhichLambdaVar,lam_pL);
        AnalysisObject.SetUpLambdaPars_pXim(DataSample,0,0,lam_pXim);

        lam_L_genuine = 1.-lam_pL[1]-lam_pL[2]-lam_pL[3]-lam_pL[4];
        lam_L_Sig0 = lam_pL[1];
        lam_L_Xim = lam_pL[2];
        lam_L_Xi0 = lam_pL[2];
        lam_L_Flat = lam_pL[3]-lam_L_Xi0;

        //has no effect in fact (not used)
        int WhichSigmaVar = rangen.Integer(3);
        if(DefaultVariation||FitSyst==false) WhichSigmaVar = 0;
        switch(WhichSigmaVar){
        case 1 : lam_S0_Flat = 0.14; break;
        case 2 : lam_S0_Flat = 0.21; break;
        default : lam_S0_Flat = 0.18; break;
        }
        lam_S0_genuine = 1.-lam_S0_Flat;

        lam_Xim_genuine = 1.-lam_pXim[1]-lam_pXim[2]-lam_pXim[3]-lam_pXim[4];
        lam_Xim_Flat = lam_pXim[3];

        //Sigma0_Feed = rangen.Integer(3);
        //if(DefaultVariation||FitSyst==false) Sigma0_Feed=1;
        DLM_Ck* Ck_pS0;
        if(Sigma0_Feed==0) Ck_pS0 = NULL;
        else if(Sigma0_Feed==1) Ck_pS0 = Ck_pS0_Chiral;
        else Ck_pS0 = Ck_pS0_ESC16;

        //Xim_Feed = rangen.Integer(2);
        Xim_Feed = 1;
        if(DefaultVariation||FitSyst==false) Xim_Feed=1;

        DLM_Ck* Ck_pXiminus;
        if(Xim_Feed==0) Ck_pXiminus = NULL;
        else Ck_pXiminus = Ck_pXim;

        DLM_Ck* Ck_pXineutral;
        if(Xim_Feed==0) Ck_pXineutral = NULL;
        else Ck_pXineutral = Ck_pXi0;

        DLM_CkDecomposition* CkDec_pL = new DLM_CkDecomposition("pLambda",5,*Ck_pL,hResolution_pL);
        DLM_CkDecomposition* CkDec_pS0 = NULL;
        if(Ck_pS0) CkDec_pS0 = new DLM_CkDecomposition("pSigma0",1,*Ck_pS0,NULL);
        DLM_CkDecomposition* CkDec_pXim = NULL;
        if(Ck_pXiminus) CkDec_pXim = new DLM_CkDecomposition("pXim",2,*Ck_pXiminus,NULL);
        DLM_CkDecomposition* CkDec_pXi0 = NULL;
        if(Ck_pXineutral) CkDec_pXi0 = new DLM_CkDecomposition("pXi0",2,*Ck_pXineutral,NULL);

        CkDec_pL->AddContribution(0,lam_L_Sig0,DLM_CkDecomposition::cFeedDown,CkDec_pS0?CkDec_pS0:NULL,CkDec_pS0?hResidual_pL_pSigma0:NULL);
        CkDec_pL->AddContribution(1,lam_L_Xim,DLM_CkDecomposition::cFeedDown,CkDec_pXim?CkDec_pXim:NULL,CkDec_pXim?hResidual_pL_pXim:NULL);
        CkDec_pL->AddContribution(2,lam_L_Xi0,DLM_CkDecomposition::cFeedDown,CkDec_pXi0?CkDec_pXi0:NULL,CkDec_pXi0?hResidual_pL_pXi0:NULL);
        //CkDec_pL->AddContribution(0,lam_pL[1],DLM_CkDecomposition::cFeedDown);
        //CkDec_pL->AddContribution(1,lam_pL[2],DLM_CkDecomposition::cFeedDown);
        CkDec_pL->AddContribution(3,lam_L_Flat,DLM_CkDecomposition::cFeedDown);
        CkDec_pL->AddContribution(4,lam_pL[4],DLM_CkDecomposition::cFake);//0.03

        if(CkDec_pS0){
            CkDec_pS0->AddContribution(0,lam_S0_Flat,DLM_CkDecomposition::cFeedDown);
        }
        if(CkDec_pXim){
            //for Xim we simplify a bit and take ALL feed-down as flat
            CkDec_pXim->AddContribution(0,lam_pXim[1]+lam_pXim[2]+lam_pXim[3],DLM_CkDecomposition::cFeedDown);
            CkDec_pXim->AddContribution(1,lam_pXim[4],DLM_CkDecomposition::cFake);
        }
        if(CkDec_pXi0){
            //for Xim we simplify a bit and take ALL feed-down as flat
            CkDec_pXi0->AddContribution(0,lam_pXim[1]+lam_pXim[2]+lam_pXim[3],DLM_CkDecomposition::cFeedDown);
            CkDec_pXi0->AddContribution(1,lam_pXim[4],DLM_CkDecomposition::cFake);
        }

        CkDec_pL->Update(true,true);
//printf("(%u) Eval at 366:\n",uIter);
//printf(" CkDec_pL = %f\n",CkDec_pL->EvalCk(366));
//printf(" CkDec_pSigma0 = %f\n",CkDec_pSigma0->EvalCk(366));
//printf(" CkDec_pXim = %f\n",CkDec_pXim->EvalCk(366));

        TString DataVar;
        //here we have new file naming, which allows for additional variations,
        //namely for the sideband correction, we can change the purity and the left-right SB fraction
        if(DataSample=="pp13TeV_HM_DimiJun20"||DataSample=="pp13TeV_HM_DimiJul20"){
            //data set variation
            WhichData = rangen.Integer(45);
            //purity variation (significant at the level of the stat. uncertainties)
            unsigned WhichSBpur = rangen.Integer(2);
            //fraction L:R variation (found insignificant)
            unsigned WhichSBl = rangen.Integer(2);

            if(DefaultVariation){
                WhichData = 0;
                WhichSBpur = 1;
                WhichSBl = 0;
            }
            if(DataSyst==false){
                WhichData = 0;
            }

            if(WhichSBpur==0){//double Gaussian IM peak
                SBpur = 0.953;
                if(WhichSBl==0){SBl=0.520;}//fit left/right in the LambdaWindow
                else{SBl=0.554;}//fit left/right outside of the LambdaWindow (width==LambdaWindow, so like 4 MeV)
            }
            else{//triple Gaussian IM peak
                SBpur = 0.963;
                if(WhichSBl==0){SBl=0.529;}
                else{SBl=0.574;}
            }
            DataVar = TString::Format("L%.0f_SL4_SR6_P%.0f_%i",SBl*100.,SBpur*100.,WhichData);
        }
        else{
            WhichData = rangen.Integer(45);
            if(DefaultVariation||DataSyst==false) WhichData = 0;
            DataVar = TString::Format("_%i",WhichData);
        }

        //TH1F* hData_pL = AnalysisObject.GetAliceExpCorrFun(DataSample,"pLambda","_0",1,false,-1);
        TH1F* hData = AnalysisObject.GetAliceExpCorrFun(DataSample,"pLambda",DataVar,2,false,-1);
        if(Bootstrap&&!DefaultVariation){
            float NewValue;
            for(unsigned uBin=0; uBin<hData->GetNbinsX(); uBin++){
                NewValue = rangen.Gaus(hData->GetBinContent(uBin+1),hData->GetBinError(uBin+1));
                if(NewValue>0) hData->SetBinContent(uBin+1, NewValue);
            }
        }

        pLambda_1_Dec = CkDec_pL;
        pLambda_1_Cat = &AB_pL;
        pLambda_1_CatFine = &ABF_pL;

        const int NumKnots = 6;
        // = new double [NumKnots];
        const double Nodes_x[] = {0.,48.,96.,192.,288.,432.};
        //1 norm/flag, than whatever number of parameters
        const unsigned NumBlPars = BASELINE_VAR==spl1?1+3+NumKnots*2:5;
        //5 default pars, outside of the BL and norm
        const unsigned NumFemtoPars = 5;
        const unsigned NumFitPars = NumFemtoPars+NumBlPars;

        TF1* fData = new TF1("fData",dimi_pL_May2020_Fitter,FitRegion_pL[0],FitRegion_pL[3],NumFitPars);

        //modify the source rad such as to represent the 2 particle radius (NOT Nolan)
        //in this way the maximum of the S(r) stays more or less at the same location
        //regardless of Alpha

        //source
        //if(WhichSourceRad!=3){
            fData->FixParameter(0,SourceRad*2./SourceAlpha);
        //}
        //else{
        //    fData->SetParameter(0,SourceRad);
        //    fData->SetParLimits(0,SourceRad*0.8,SourceRad*1.2);
        //}
        fData->FixParameter(1,SourceAlpha);

        //cusp
        //if(WhichCuspStrength==3){
        //    fData->SetParameter(2,CuspWeight);
        //    fData->SetParLimits(2,0.2,0.5);
        //}
        //else fData->FixParameter(2,CuspWeight);
        fData->FixParameter(2,CuspWeight);

        //cutoff
        fData->FixParameter(3,CkCutOff);

        fData->FixParameter(4,CkConv);

        //norm
        if(BASELINE_VAR==spl1){
            fData->FixParameter(5,-1e6);
        }
        else{
            fData->SetParameter(5,1);
            fData->SetParLimits(5,0.5,2.0);
        }

        if(BASELINE_VAR==pol0s){
            for(int i=7; i<=9; i++) fData->FixParameter(i,-1e6);
            fData->FixParameter(6,0);
        }
        else if(BASELINE_VAR==pol1s){
            for(int i=7; i<=9; i++) fData->FixParameter(i,-1e6);
            fData->SetParameter(6,0);
            fData->SetParLimits(6,-1e-2,1e-2);
        }
        else if(BASELINE_VAR==pol2s){
            for(int i=8; i<=9; i++) fData->FixParameter(i,-1e6);
            fData->SetParameter(6,200);
            fData->SetParLimits(6,0,400);
            fData->SetParameter(7,0);
            fData->SetParLimits(7,-1e-4,1e-4);
        }
        else if(BASELINE_VAR==pol3s||BASELINE_VAR==pol3e){
            for(int i=9; i<=9; i++) fData->FixParameter(i,-1e6);
            fData->SetParameter(6,0);
            fData->SetParLimits(6,-100000,100);
            fData->SetParameter(7,100);
            fData->SetParLimits(7,0,400);
            fData->SetParameter(8,0);
            fData->SetParLimits(8,-1e-6,1e-6);
        }
        else if(BASELINE_VAR==dpol2s||BASELINE_VAR==dpol2e){
            for(int i=8; i<=9; i++) fData->FixParameter(i,-1e6);
            fData->FixParameter(6,0);
            fData->SetParameter(7,0);
            fData->SetParLimits(7,-1e-4,1e-4);
        }
        else if(BASELINE_VAR==dpol3s||BASELINE_VAR==dpol3e){
            for(int i=9; i<=9; i++) fData->FixParameter(i,-1e6);
            fData->FixParameter(6,0);
            fData->SetParameter(7,200);
            fData->SetParLimits(7,0,400);
            fData->SetParameter(8,0);
            fData->SetParLimits(8,-1e-6,1e-6);
        }
        //
        else if(BASELINE_VAR==dpol4e){
            //fData->SetParameter(6,0);
            //fData->SetParLimits(6,-100000,400);
            //fData->SetParameter(7,275);
            //fData->SetParLimits(7,100,1000);
            //fData->SetParameter(8,0);
            //fData->SetParLimits(8,-1e-7,1e-7);
            fData->SetParameter(6,90);
            fData->SetParLimits(6,0,180);
            fData->SetParameter(7,240);
            fData->SetParLimits(7,180,360);
            fData->SetParameter(8,0);
            fData->SetParLimits(8,-1e-6,1e-6);
            fData->FixParameter(9,1);
        }
        //spl1
        else{
            fData->FixParameter(6,NumKnots);
            //derivative at the firs knot
            fData->FixParameter(7,0);
            fData->SetParameter(8,0);
            fData->SetParLimits(8,-1e-3,1e-3);
            for(unsigned uKnot=0; uKnot<NumKnots; uKnot++){
                double HistVal = hData->GetBinContent(hData->FindBin(Nodes_x[uKnot]));
                fData->FixParameter(6+3+uKnot,Nodes_x[uKnot]);
                fData->SetParameter(6+3+NumKnots+uKnot,HistVal);
                fData->SetParLimits(6+3+NumKnots+uKnot,0,HistVal*2);
            }
        }

        hData->Fit(fData,"Q, S, N, R, M");
        fData->SetNpx(1024);

//usleep(2000e3);
//delete fData;
//delete CkDec_pL;
//delete CkDec_pSigma0;
//delete CkDec_pXim;
//continue;

        CkDec_pL->Update(true,true);

        SourceRad = fData->GetParameter(0)/2.*SourceAlpha;
        CuspWeight = fData->GetParameter(2);
        CkCutOff = fData->GetParameter(3);
        CkConv = fData->GetParameter(4);
        Norm = fData->GetParameter(5);
        pp1 = fData->GetParameter(6);
        pp2 = fData->GetParameter(7);
        pp3 = fData->GetParameter(8);
        pp4 = fData->GetParameter(9);

        pval = fData->GetProb();

        TF1* fBaseline = new TF1("fBaseline",dimi_pL_May2020_FitterBl,FitRegion_pL[0],FitRegion_pL[3],NumBlPars);
        for(unsigned uPar=0; uPar<NumBlPars; uPar++) fBaseline->FixParameter(uPar,fData->GetParameter(NumFemtoPars+uPar));

        TF1* fFemto = new TF1("fFemto",dimi_pL_May2020_FitterFemto,FitRegion_pL[0],FitRegion_pL[3],NumFemtoPars);
        for(unsigned uPar=0; uPar<NumFemtoPars; uPar++) fFemto->FixParameter(uPar,fData->GetParameter(uPar));

        for(unsigned uBin=0; uBin<hData->GetNbinsX(); uBin++){

            double MOM = hData->GetBinCenter(uBin+1);

            if(MOM>500) break;

            gData.SetPoint(uBin,MOM,hData->GetBinContent(uBin+1));
            gData.SetPointError(uBin,0,hData->GetBinError(uBin+1));

            fitData.SetPoint(uBin,MOM,fData->Eval(MOM));
            fitData.SetPointError(uBin,0,0);

            fitFemto.SetPoint(uBin,MOM,fFemto->Eval(MOM));
            fitFemto.SetPointError(uBin,0,0);

            fitGenuine.SetPoint(uBin,MOM,Ck_pL->Eval(MOM));
            fitGenuine.SetPointError(uBin,0,0);

            fitSignal_pL.SetPoint(uBin,MOM,CkDec_pL->EvalSignalSmearedMain(MOM)+1);
            fitSignal_pL.SetPointError(uBin,0,0);

            fitSignal_pL_pS0.SetPoint(uBin,MOM,CkDec_pL->EvalSignalSmearedChild(0,MOM)+1);
            fitSignal_pL_pS0.SetPointError(uBin,0,0);

            fitSignal_pL_pXim.SetPoint(uBin,MOM,CkDec_pL->EvalSignalSmearedChild(1,MOM)+1);
            fitSignal_pL_pXim.SetPointError(uBin,0,0);

            fitSignal_pL_pXi0.SetPoint(uBin,MOM,CkDec_pL->EvalSignalSmearedChild(2,MOM)+1);
            fitSignal_pL_pXi0.SetPointError(uBin,0,0);

            fitBaseline.SetPoint(uBin,MOM,fBaseline->Eval(MOM));
            fitBaseline.SetPointError(uBin,0,0);

        }

        OutputFile->cd();
        plambdaTree->Fill();

        delete hData;
        delete fData;
        delete fBaseline;
        delete fFemto;
        delete CkDec_pL;
        if(CkDec_pS0) delete CkDec_pS0;
        if(CkDec_pXim) delete CkDec_pXim;
        if(CkDec_pXi0) delete CkDec_pXi0;
    }
printf("\n");

    OutputFile->cd();
    //ntResult->Write("",TObject::kOverwrite);
    plambdaTree->Write();

    delete Ck_pL;
    delete Ck_pXim;
    delete Ck_pXi0;
    delete Ck_pS0_Chiral;
    delete Ck_pS0_ESC16;
    if(hResolution_pL){delete hResolution_pL;hResolution_pL=NULL;}
    delete hResidual_pL_pSigma0;
    delete hResidual_pL_pXim;
    delete hResidual_pL_pXi0;
    delete plambdaTree;
    delete OutputFile;
    delete [] MomBins_pL;
    delete [] FitRegion_pL;

}


void GetIterCombo(unsigned* WhichBin,const float& Sigma0_Feed,const float& SourceSize,const float& lam_L_genuine){
    WhichBin[0] = Sigma0_Feed;
    if(SourceSize<1.01) WhichBin[1] = 0;
    else if(SourceSize>1.01&&SourceSize<1.03) WhichBin[1] = 1;
    else if(SourceSize>1.05&&SourceSize<1.07) WhichBin[1] = 2;
    else WhichBin[1] = 3;
    if(lam_L_genuine<0.44) WhichBin[2] = 0;
    else if(lam_L_genuine>0.44&&lam_L_genuine<0.46) WhichBin[2] = 1;
    else WhichBin[2] = 2;
}
void GetIterCombo270520(unsigned* WhichBin,const float& Sigma0_Feed,const float& SourceSize,const float& lam_L_genuine,
                        const float& CuspWeight, const float& SourceAlpha){
    WhichBin[0] = Sigma0_Feed;
    if(SourceSize<1.00) WhichBin[1] = 0;
    else if(SourceSize<1.04) WhichBin[1] = 1;
    else if(SourceSize<1.08) WhichBin[1] = 2;
    else WhichBin[1] = 3;
    if(lam_L_genuine<0.42) WhichBin[2] = 0;
    else if(lam_L_genuine<0.44) WhichBin[2] = 1;
    else if(lam_L_genuine<0.46) WhichBin[2] = 2;
    else if(lam_L_genuine<0.48) WhichBin[2] = 3;
    else WhichBin[2] = 4;
    if(CuspWeight<0.23) WhichBin[3] = 0;
    else if(CuspWeight<0.30) WhichBin[3] = 1;
    else if(CuspWeight<0.37) WhichBin[3] = 2;
    else if(CuspWeight<0.44) WhichBin[3] = 3;
    else WhichBin[3] = 4;
    if(SourceAlpha>1.9) WhichBin[4] = 0;
    else if(SourceAlpha<1.7) WhichBin[4] = 1;
    else WhichBin[4] = 2;
}
void GetIterCombo290520(unsigned* WhichBin,const float& SourceSize,const float& lam_L_genuine,
                        const float& CuspWeight, const float& SourceAlpha){
    if(SourceSize<1.00) WhichBin[1] = 0;
    else if(SourceSize<1.04) WhichBin[1] = 1;
    else if(SourceSize<1.08) WhichBin[1] = 2;
    else WhichBin[1] = 3;
    if(lam_L_genuine<0.42) WhichBin[2] = 0;
    else if(lam_L_genuine<0.44) WhichBin[2] = 1;
    else if(lam_L_genuine<0.46) WhichBin[2] = 2;
    else if(lam_L_genuine<0.48) WhichBin[2] = 3;
    else WhichBin[2] = 4;
    if(CuspWeight<0.23) WhichBin[3] = 0;
    else if(CuspWeight<0.30) WhichBin[3] = 1;
    else if(CuspWeight<0.37) WhichBin[3] = 2;
    else if(CuspWeight<0.44) WhichBin[3] = 3;
    else WhichBin[3] = 4;
    if(SourceAlpha>1.9) WhichBin[4] = 0;
    else if(SourceAlpha<1.7) WhichBin[4] = 1;
    else WhichBin[4] = 2;
}
void GetIterCombo040620(unsigned* WhichBin,const float& SourceSize,const float& lam_L_genuine,
                        const float& CuspWeight, const float& SourceAlpha){
    if(SourceSize<0.96) WhichBin[0] = 0;
    else if(SourceSize<1.00) WhichBin[0] = 1;
    else if(SourceSize<1.04) WhichBin[0] = 2;
    else if(SourceSize<1.08) WhichBin[0] = 3;
    else WhichBin[0] = 4;
    if(lam_L_genuine<0.4) WhichBin[1] = 0;
    else if(lam_L_genuine<0.402) WhichBin[1] = 1;
    else if(lam_L_genuine<0.405) WhichBin[1] = 2;
    else if(lam_L_genuine<0.419) WhichBin[1] = 3;
    else if(lam_L_genuine<0.422) WhichBin[1] = 4;
    else if(lam_L_genuine<0.424) WhichBin[1] = 5;
    else if(lam_L_genuine<0.440) WhichBin[1] = 6;
    else if(lam_L_genuine<0.443) WhichBin[1] = 7;
    else if(lam_L_genuine<0.445) WhichBin[1] = 8;
    else if(lam_L_genuine<0.463) WhichBin[1] = 9;
    else if(lam_L_genuine<0.466) WhichBin[1] = 10;
    else if(lam_L_genuine<0.468) WhichBin[1] = 11;
    else if(lam_L_genuine<0.489) WhichBin[1] = 12;
    else if(lam_L_genuine<0.492) WhichBin[1] = 13;
    else WhichBin[1] = 14;
    if(CuspWeight<0.23) WhichBin[2] = 0;
    else if(CuspWeight<0.30) WhichBin[2] = 1;
    else if(CuspWeight<0.37) WhichBin[2] = 2;
    else if(CuspWeight<0.44) WhichBin[2] = 3;
    else WhichBin[2] = 4;
    if(SourceAlpha>1.65) WhichBin[3] = 0;
    else if(SourceAlpha<1.75) WhichBin[3] = 1;
    else if(SourceAlpha<1.85) WhichBin[3] = 2;
    else if(SourceAlpha<1.95) WhichBin[3] = 3;
    else WhichBin[3] = 4;
}
void GetIterCombo170620(unsigned* WhichBin,const float& SourceSize,const float& lam_L_genuine,
                        const float& CuspWeight, const float& SourceAlpha){
    if(SourceSize<0.96) WhichBin[0] = 0;
    else if(SourceSize<1.00) WhichBin[0] = 1;
    else if(SourceSize<1.04) WhichBin[0] = 2;
    else if(SourceSize<1.08) WhichBin[0] = 3;
    else WhichBin[0] = 4;
    if(lam_L_genuine<0.44) WhichBin[1] = 0;
    else if(lam_L_genuine<0.46) WhichBin[1] = 1;
    else if(lam_L_genuine<0.48) WhichBin[1] = 2;
    else if(lam_L_genuine<0.51) WhichBin[1] = 3;
    else WhichBin[1] = 4;
    if(CuspWeight<0.23) WhichBin[2] = 0;
    else if(CuspWeight<0.30) WhichBin[2] = 1;
    else if(CuspWeight<0.37) WhichBin[2] = 2;
    else if(CuspWeight<0.44) WhichBin[2] = 3;
    else WhichBin[2] = 4;
    if(SourceAlpha>1.65) WhichBin[3] = 0;
    else if(SourceAlpha<1.75) WhichBin[3] = 1;
    else if(SourceAlpha<1.85) WhichBin[3] = 2;
    else if(SourceAlpha<1.95) WhichBin[3] = 3;
    else WhichBin[3] = 4;
}

//sig_lam SIG to LAM ratio
//xi_siglam XI (Xi0+Xim) to (Lam+Sig) ratio
void GetIterCombo040720(unsigned* WhichBin,const float& SourceSize,
                        const float& SBpur,const float& sig_lam,const float& xi_siglam,
                        const float& CuspWeight, const float& SourceAlpha, const float& CkConv){
    if(SourceSize<0.96) WhichBin[0] = 0;
    else if(SourceSize<1.00) WhichBin[0] = 1;
    else if(SourceSize<1.04) WhichBin[0] = 2;
    else if(SourceSize<1.08) WhichBin[0] = 3;
    else WhichBin[0] = 4;
    if(SBpur<0.96) WhichBin[1] = 0;
    else WhichBin[1] = 1;
    if(sig_lam<0.23) WhichBin[2] = 0;
    else if(sig_lam<0.30) WhichBin[2] = 1;
    else if(sig_lam<0.37) WhichBin[2] = 2;
    else if(sig_lam<0.44) WhichBin[2] = 3;
    else WhichBin[2] = 4;
    if(xi_siglam<0.33) WhichBin[3] = 0;
    else WhichBin[3] = 1;
    if(CuspWeight<0.23) WhichBin[4] = 0;
    else if(CuspWeight<0.30) WhichBin[4] = 1;
    else if(CuspWeight<0.37) WhichBin[4] = 2;
    else if(CuspWeight<0.44) WhichBin[4] = 3;
    else WhichBin[4] = 4;
    if(SourceAlpha>1.65) WhichBin[5] = 0;
    else if(SourceAlpha<1.75) WhichBin[5] = 1;
    else if(SourceAlpha<1.85) WhichBin[5] = 2;
    else if(SourceAlpha<1.95) WhichBin[5] = 3;
    else WhichBin[5] = 4;
    if(CkConv<600) WhichBin[6] = 0;
    else if(CkConv<800) WhichBin[6] = 1;
    else WhichBin[6] = 2;
}

/*
Sigma0_Feed = large (3x)
SourceSize = large (4x)
CuspWeight = small (4x, default var is best anyways)
CkConv = small (3x, and its fitted, so it play no role actually)
lam_L_genuine = large (3x)
lam_S0_genuine = small (3x, best solution is the lowest var, but the default is almost identical. Only the largest var sticks a bit out, exclude?)
//by the way, lam_S0_genuine only plays a role for non-flat Sigma0_Feed, and the best solution is the flat guy. So, in fact, we can ignore this as well

Get the best solution out of the variations with largest effect on the final result (Sigma0_Feed, SourceSize, lam_L_genuine)

out of this, build the chi2/ndf distribution. The mean value should be a good estimate on the 'true' chi2 value, and hence sigma
(essentially the rest of the above variations can be neglected, while all other variations should be 'random' in nature)
*/
//WhichSourceAlpha=0 Gauss
//WhichSourceAlpha=1 1.8
//WhichSourceAlpha=2 1.6
void Plot_pL_SystematicsMay2020(const int& SIGMA_FEED, const int& WhichSourceAlpha, TString InputFolder, TString InputFileName){
    //TString InputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/pL_SystematicsMay2020/BatchFarm/200520/";
    //TString InputFileName = "Gauss_4000.root";
    //TString InputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/pL_SystematicsMay2020/BatchFarm/220520_Levy/";
    //TString InputFileName = "Levy_8192.root";
    //TString InputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/pL_SystematicsMay2020/190520/";
    //TString InputFileName = "Levy_1319x128_19x128.root";

    TString OutputFolder;
    if(SIGMA_FEED==0) OutputFolder = InputFolder+"Flat_Sigma/";
    else if(SIGMA_FEED==1) OutputFolder = InputFolder+"NLO_Sigma/";
    else if(SIGMA_FEED==2) OutputFolder = InputFolder+"ESC16_Sigma/";
    else if(SIGMA_FEED==-1) OutputFolder = InputFolder+"Full_Sigma/";
    else OutputFolder = InputFolder;
    /// as cross check: pol0s,pol1s,pol2s,dpol2s,pol3e
    /// as systematics: dpol3s,dpol3e,dpol4e
    /// for fun: spl1 (splines at 0-48-96-192-288-432)
    enum BLTYPE              { pol0s,pol1s,pol2s,pol3s,dpol2s,dpol3s,dpol4s,pol2e,pol3e,dpol2e,dpol3e,dpol4e,spl1 };
    //const bool IncludeBl[] = { true, true, true, false,false, true,  false, false,false,false, true,  false};
    const unsigned NumBaselines = 9;
    int* BlFlag = new int [NumBaselines];
    bool* IncludeBl = new bool [NumBaselines];
    BlFlag[0] = pol0s; IncludeBl[0] = true;
    BlFlag[1] = pol1s; IncludeBl[1] = true;
    BlFlag[2] = pol2s; IncludeBl[2] = true;
    BlFlag[3] = dpol2s; IncludeBl[3] = false;
    BlFlag[4] = pol3e; IncludeBl[4] = false;
    BlFlag[5] = dpol3s; IncludeBl[5] = true;
    BlFlag[6] = dpol3e; IncludeBl[6] = true;
    BlFlag[7] = dpol4e; IncludeBl[7] = false;
    BlFlag[8] = spl1; IncludeBl[8] = false;
    TString* BlName1 = new TString [NumBaselines+1];
    TString* BlName2 = new TString [NumBaselines+1];
    TString* BlDescr = new TString [NumBaselines+1];
    for(unsigned uBl=0; uBl<=NumBaselines; uBl++) BlName1[uBl] = "";
    for(unsigned uBl=0; uBl<=NumBaselines; uBl++) BlName2[uBl] = "";

    BlName1[0] = "Constant baseline";
    BlName2[0] = "Standard fit range";
    BlDescr[0] = "pol0s";

    BlName1[1] = "Linear baseline";
    BlName2[1] = "Standard fit range";
    BlDescr[1] = "pol1s";

    BlName1[2] = "Quadratic baseline";
    BlName2[2] = "Standard fit range";
    BlDescr[2] = "pol2s";

    BlName1[3] = "Quadratic baseline (constrained)";
    BlName2[3] = "Standard fit range";
    BlDescr[3] = "dpol2s";

    BlName1[4] = "Cubic baseline";
    BlName2[4] = "Extended fit range";
    BlDescr[4] = "pol3e";

    BlName1[5] = "Cubic baseline (constrained)";
    BlName2[5] = "Standard fit range";
    BlDescr[5] = "dpol3s";

    BlName1[6] = "Cubic baseline (constrained)";
    BlName2[6] = "Extended fit range";
    BlDescr[6] = "dpol3e";

    BlName1[7] = "Quartic baseline (constrained)";
    BlName2[7] = "Extended fit range";
    BlDescr[7] = "dpol4e";

    BlName1[8] = "Spline (type 1)";
    BlName2[8] = "Extended fit range";
    BlDescr[8] = "spl1";

    BlName1[9] = "All baselines";
    BlName2[9] = "";
    BlDescr[9] = "AllBl";

    const unsigned NumPotVars = 11;
    TString* PotName1 = new TString [NumPotVars+1];
    PotName1[0] = "NLO13 (500)";
    PotName1[1] = "NLO13 (550)";
    PotName1[2] = "NLO13 (600)";
    PotName1[3] = "NLO13 (650)";
    PotName1[4] = "NLO19 (500)";
    PotName1[5] = "NLO19 (550)";
    PotName1[6] = "NLO19 (600)";
    PotName1[7] = "NLO19 (650)";
    PotName1[8] = "LO13 (600)";
    PotName1[9] = "Usmani (original)";
    PotName1[10] = "Usmani (modified)";
    PotName1[11] = "All (broken?)";
    TString* PotName2 = new TString [NumPotVars+1];
    for(unsigned uPot=0; uPot<=NumPotVars; uPot++) PotName2[uPot] = "s,d waves";
    int* PotFlag = new int [NumPotVars];
    bool* IncludePot = new bool [NumPotVars];
    PotFlag[0] = 1500; IncludePot[0] = true;
    PotFlag[1] = 1550; IncludePot[1] = true;
    PotFlag[2] = 1600; IncludePot[2] = true;
    PotFlag[3] = 1650; IncludePot[3] = true;
    PotFlag[4] = 11500; IncludePot[4] = true;
    PotFlag[5] = 11550; IncludePot[5] = true;
    PotFlag[6] = 11600; IncludePot[6] = true;
    PotFlag[7] = 11650; IncludePot[7] = true;
    PotFlag[8] = -11600; IncludePot[8] = true;
    TString* PotDescr = new TString [NumPotVars+1];
    PotDescr[0] = "NLO13-500";
    PotDescr[1] = "NLO13-550";
    PotDescr[2] = "NLO13-600";
    PotDescr[3] = "NLO13-650";
    PotDescr[4] = "NLO19-500";
    PotDescr[5] = "NLO19-550";
    PotDescr[6] = "NLO19-600";
    PotDescr[7] = "NLO19-650";
    PotDescr[8] = "LO13-600";
    PotDescr[9] = "Usmani";
    PotDescr[10] = "UsmaniFit";
    PotDescr[11] = "All";

    const unsigned NumSourceVars = 1;
    TString* SourceName1 = new TString [NumSourceVars];
    if(WhichSourceAlpha==0) SourceName1[0] = "Gaussian core source";
    else if(WhichSourceAlpha==1) SourceName1[0] = "Levy core source (#alpha=1.8)";
    else if(WhichSourceAlpha==2) SourceName1[0] = "Levy core source (#alpha=1.6)";
    else SourceName1[0] = "Levy core source";
    //SourceName1[1] = "Levy core source";

    TFile* InputFile = new TFile(InputFolder+InputFileName,"read");
    TTree* plambdaTree = (TTree*)InputFile->Get("plambdaTree");

	TGraphErrors* gData = NULL;
	TGraphErrors* fitData = NULL;
	TGraphErrors* fitFemto = NULL;
	TGraphErrors* fitGenuine = NULL;
	TGraphErrors* fitSignal_pL = NULL;
	TGraphErrors* fitSignal_pL_pS0 = NULL;
	TGraphErrors* fitSignal_pL_pXim = NULL;
	TGraphErrors* fitBaseline = NULL;
	TString* DataSample = NULL;
	TString* SourceDescription = NULL;
	unsigned SEED;
	unsigned BASELINE_VAR;
	int POT_VAR;
	bool DataSyst;
	bool FitSyst;
	bool Bootstrap;
	float SourceSize;
	float SourceAlpha;
    float CuspWeight;
    float CkConv;
    float CkCutOff;
    float Norm;
    float pp1;
    float pp2;
    float pp3;
    float pp4;
    float kFemtoMin;
    float kFemtoMax;
    float kLongMin;
    float kLongMax;
    float lam_L_genuine;
    float lam_L_Sig0;
    float lam_L_Xim;
    float lam_L_Xi0;
    float lam_L_Flat;
    float lam_S0_genuine;
    float lam_S0_Flat;
    float lam_Xim_genuine;
    float lam_Xim_Flat;
    float pval;
    unsigned WhichData;
    bool DefaultVariation;
    int Sigma0_Feed;
    int Xim_Feed;

	plambdaTree->SetBranchAddress("gData",&gData);
	plambdaTree->SetBranchAddress("fitData",&fitData);
	plambdaTree->SetBranchAddress("fitFemto",&fitFemto);
	plambdaTree->SetBranchAddress("fitGenuine",&fitGenuine);
	plambdaTree->SetBranchAddress("fitSignal_pL",&fitSignal_pL);
	plambdaTree->SetBranchAddress("fitSignal_pL_pS0",&fitSignal_pL_pS0);
	plambdaTree->SetBranchAddress("fitSignal_pL_pXim",&fitSignal_pL_pXim);
	plambdaTree->SetBranchAddress("fitBaseline",&fitBaseline);
	plambdaTree->SetBranchAddress("DataSample",&DataSample);
	plambdaTree->SetBranchAddress("SourceDescription",&SourceDescription);
	plambdaTree->SetBranchAddress("SEED",&SEED);
	plambdaTree->SetBranchAddress("BASELINE_VAR",&BASELINE_VAR);
	plambdaTree->SetBranchAddress("POT_VAR",&POT_VAR);
	plambdaTree->SetBranchAddress("Sigma0_Feed",&Sigma0_Feed);
	plambdaTree->SetBranchAddress("Xim_Feed",&Xim_Feed);
	plambdaTree->SetBranchAddress("DataSyst",&DataSyst);
    plambdaTree->SetBranchAddress("FitSyst", &FitSyst);//
    plambdaTree->SetBranchAddress("Bootstrap", &Bootstrap);//
	plambdaTree->SetBranchAddress("SourceSize",&SourceSize);
	plambdaTree->SetBranchAddress("SourceAlpha",&SourceAlpha);
	plambdaTree->SetBranchAddress("CuspWeight",&CuspWeight);
	plambdaTree->SetBranchAddress("CkCutOff",&CkCutOff);
	plambdaTree->SetBranchAddress("CkConv",&CkConv);
	plambdaTree->SetBranchAddress("Norm",&Norm);
	plambdaTree->SetBranchAddress("pp1",&pp1);
	plambdaTree->SetBranchAddress("pp2",&pp2);
	plambdaTree->SetBranchAddress("pp3",&pp3);
	plambdaTree->SetBranchAddress("pp4",&pp4);
	plambdaTree->SetBranchAddress("kFemtoMin",&kFemtoMin);
	plambdaTree->SetBranchAddress("kFemtoMax",&kFemtoMax);
	plambdaTree->SetBranchAddress("kLongMin",&kLongMin);
	plambdaTree->SetBranchAddress("kLongMax",&kLongMax);
	plambdaTree->SetBranchAddress("lam_L_genuine",&lam_L_genuine);
	plambdaTree->SetBranchAddress("lam_L_Sig0",&lam_L_Sig0);
	plambdaTree->SetBranchAddress("lam_L_Xim",&lam_L_Xim);
	plambdaTree->SetBranchAddress("lam_L_Xi0",&lam_L_Xi0);
	plambdaTree->SetBranchAddress("lam_L_Flat",&lam_L_Flat);
	plambdaTree->SetBranchAddress("lam_S0_genuine",&lam_S0_genuine);
	plambdaTree->SetBranchAddress("lam_S0_Flat",&lam_S0_Flat);
	plambdaTree->SetBranchAddress("lam_Xim_genuine",&lam_Xim_genuine);
	plambdaTree->SetBranchAddress("lam_Xim_Flat",&lam_Xim_Flat);
	plambdaTree->SetBranchAddress("WhichData",&WhichData);
	plambdaTree->SetBranchAddress("pval",&pval);
	plambdaTree->SetBranchAddress("DefaultVariation",&DefaultVariation);

	const unsigned NumEntries = plambdaTree->GetEntries();

    const unsigned MaxNumBins = 50;
    const double xBinWidth = 12;
    TH2F* hdata_toterr = new TH2F("hdata_toterr","hdata_toterr",MaxNumBins,0,MaxNumBins,2048,0,4);
    TGraphErrors* gdata_toterr = new TGraphErrors();


    unsigned NumBestSolutions = 0;
    //in bins of Sigma0_Feed, SourceSize, lam_L_genuine

    const unsigned Num_Sigma0_Feed=3;
    const unsigned Num_SourceSize=4;
    const unsigned Num_lam_L_genuine=5;
    const unsigned Num_CuspWeight=5;
    const unsigned Num_SourceAlpha=1;

    //DLM_Histo<float> NDF_300;
    //NDF_300.SetUp(3);
    //NDF_300.SetUp(0,Num_Sigma0_Feed,-0.5,double(Num_Sigma0_Feed)+0.5);
    //NDF_300.SetUp(1,Num_SourceSize,-0.5,double(Num_Sigma0_Feed)+0.5);
    //NDF_300.SetUp(2,Num_lam_L_genuine,-0.5,double(Num_Sigma0_Feed)+0.5);
    //NDF_300.Initialize();

    for(unsigned uEntry=0; uEntry<NumEntries; uEntry++){
        plambdaTree->GetBranch("gData")->GetEntry(uEntry);
        double xVal,yVal;
        for(unsigned uBin=0; uBin<gData->GetN(); uBin++){
            if(uBin>=MaxNumBins) break;
            gData->GetPoint(uBin,xVal,yVal);
//if(uBin==0){
//printf("DATA: C(%.0f) = %.3f\n",xVal,yVal);
//usleep(0.25e6);
//}
            gdata_toterr->SetPoint(uBin,xVal,yVal);
            hdata_toterr->Fill(double(uBin)+0.5,yVal);
        }
    }
    for(unsigned uBin=0; uBin<gdata_toterr->GetN(); uBin++){
        double xVal,yVal;
        double yMean,yStddev;
        TH1* hProj = hdata_toterr->ProjectionY("hProj",uBin+1,uBin+1);
        yMean = hProj->GetMean();
        yStddev = hProj->GetStdDev();
        gdata_toterr->GetPoint(uBin,xVal,yVal);
        //printf("uBin %u: C(%.0f) = %.3f +/- %.3f\n",uBin,xVal,yMean,yStddev);
        gdata_toterr->SetPoint(uBin,xVal,yMean);
        gdata_toterr->SetPointError(uBin,xBinWidth*0.5,yStddev);
        delete hProj;
    }
//for(unsigned uBl=0; uBl<=NumBaselines; uBl++){
//printf("IBL[%u] = %i\n",uBl,IncludeBl[uBl]);
//}
    for(unsigned uBl=0; uBl<=NumBaselines; uBl++){
        if(!IncludeBl[uBl]&&uBl<NumBaselines) continue;
//if(uBl!=1)continue;
        for(unsigned uPot=0; uPot<=NumPotVars; uPot++){
            if(!IncludePot[uPot]&&uPot<NumPotVars) continue;
//if(uPot!=6)continue;
            TH2F* hnsigma = new TH2F("hnsigma","hnsigma",MaxNumBins,0,MaxNumBins,2048,-32,32);
            TH2F* hnsigma_der = new TH2F("hnsigma_der","hnsigma_der",MaxNumBins,0,MaxNumBins,2048,-32,32);
            TH2F* hfdata = new TH2F("hfdata","hfdata",MaxNumBins,0,MaxNumBins,2048,0,4);
            TH2F* hbl = new TH2F("hbl","hbl",MaxNumBins,0,MaxNumBins,2048,0.5,1.5);

            TH2F* hb_nsigma = new TH2F("hb_nsigma","hb_nsigma",MaxNumBins,0,MaxNumBins,2048,-32,32);
            TH2F* hb_nsigma_der = new TH2F("hb_nsigma_der","hb_nsigma_der",MaxNumBins,0,MaxNumBins,2048,-32,32);
            TH2F* hb_fdata = new TH2F("hb_fdata","hb_fdata",MaxNumBins,0,MaxNumBins,2048,0,4);
            TH2F* hb_bl = new TH2F("hb_bl","hb_bl",MaxNumBins,0,MaxNumBins,2048,0.5,1.5);

            TH1F* hnsigma_0_100 = new TH1F("hnsigma_0_100","hnsigma_0_100",2048,-32,32);
            TH1F* hnsigma_100_200 = new TH1F("hnsigma_100_200","hnsigma_100_200",2048,-32,32);
            TH1F* hnsigma_200_300 = new TH1F("hnsigma_200_300","hnsigma_200_300",2048,-32,32);
            TH1F* hnsigma_0_300 = new TH1F("hnsigma_0_300","hnsigma_0_300",2048,-32,32);
            TH1F* hnsigma_der_0_100 = new TH1F("hnsigma_der_0_100","hnsigma_der_0_100",2048,-32,32);
            TH1F* hnsigma_der_100_200 = new TH1F("hnsigma_der_100_200","hnsigma_der_100_200",2048,-32,32);
            TH1F* hnsigma_der_200_300 = new TH1F("hnsigma_der_200_300","hnsigma_der_200_300",2048,-32,32);
            TH1F* hnsigma_der_0_300 = new TH1F("hnsigma_der_0_300","hnsigma_der_0_300",2048,-32,32);

            TH1F* htot_nsigma_300 = new TH1F("htot_nsigma_300","htot_nsigma_300",2048,-32,32);
            TH1F* htot_nsigma_der_300 = new TH1F("htot_nsigma_der_300","htot_nsigma_der_300",2048,-32,32);
            float nsigma_0_100;
            float nsigma_100_200;
            float nsigma_200_300;
            float nsigma_0_300;
            float nsigma_der_0_100;
            float nsigma_der_100_200;
            float nsigma_der_200_300;
            float nsigma_der_0_300;
            int NDF_0_100;
            int NDF_100_200;
            int NDF_200_300;
            int NDF_0_300;
            for(unsigned uBin=0; uBin<MaxNumBins; uBin++){
                for(unsigned uNSB=0; uNSB<2048; uNSB++){
                    hnsigma->SetBinContent(uBin+1,uNSB+1,0);
                    hnsigma_der->SetBinContent(uBin+1,uNSB+1,0);
                    hb_nsigma->SetBinContent(uBin+1,uNSB+1,0);
                    hb_nsigma_der->SetBinContent(uBin+1,uNSB+1,0);
                }
            }


            //original dimensions: 3 -> Sigma0_Feed, SourceSize,lam_L_genuine
            //now I have added CuspWeight and SourceAlpha
            //switch off some contributions by setting their number to 1
            DLM_Histo<float> Chi2_300;
            Chi2_300.SetUp(5);
            Chi2_300.SetUp(0,Num_Sigma0_Feed,-0.5,double(Num_Sigma0_Feed)+0.5);
            Chi2_300.SetUp(1,Num_SourceSize,-0.5,double(Num_Sigma0_Feed)+0.5);
            Chi2_300.SetUp(2,Num_lam_L_genuine,-0.5,double(Num_Sigma0_Feed)+0.5);
            Chi2_300.SetUp(3,Num_CuspWeight,-0.5,double(Num_CuspWeight)+0.5);
            Chi2_300.SetUp(4,Num_SourceAlpha,-0.5,double(Num_SourceAlpha)+0.5);
            Chi2_300.Initialize();

            DLM_Histo<float> NDF_300;
            NDF_300.SetUp(3);
            NDF_300.SetUp(0,Num_Sigma0_Feed,-0.5,double(Num_Sigma0_Feed)+0.5);
            NDF_300.SetUp(1,Num_SourceSize,-0.5,double(Num_Sigma0_Feed)+0.5);
            NDF_300.SetUp(2,Num_lam_L_genuine,-0.5,double(Num_Sigma0_Feed)+0.5);
            NDF_300.SetUp(3,Num_CuspWeight,-0.5,double(Num_CuspWeight)+0.5);
            NDF_300.SetUp(4,Num_SourceAlpha,-0.5,double(Num_SourceAlpha)+0.5);
            NDF_300.Initialize();

            DLM_Histo<float> NumChi2Entries_300;
            NumChi2Entries_300.SetUp(3);
            NumChi2Entries_300.SetUp(0,Num_Sigma0_Feed,-0.5,double(Num_Sigma0_Feed)+0.5);
            NumChi2Entries_300.SetUp(1,Num_SourceSize,-0.5,double(Num_Sigma0_Feed)+0.5);
            NumChi2Entries_300.SetUp(2,Num_lam_L_genuine,-0.5,double(Num_Sigma0_Feed)+0.5);
            NumChi2Entries_300.SetUp(3,Num_CuspWeight,-0.5,double(Num_CuspWeight)+0.5);
            NumChi2Entries_300.SetUp(4,Num_SourceAlpha,-0.5,double(Num_SourceAlpha)+0.5);
            NumChi2Entries_300.Initialize();

            unsigned UsedNumBins = 0;
            unsigned NumUsefulEntries=0;
            for(unsigned uEntry=0; uEntry<NumEntries; uEntry++){
                //Chi2_300[uEntry] = 1e6;
                //NDF_300[uEntry] = 0;
                plambdaTree->GetEntry(uEntry);
                //if(uPot<NumPotVars&&POT_VAR!=PotFlag[uPot]) continue;
                //if(uBl<NumBaselines&&BASELINE_VAR!=BlFlag[uBl]) continue;
                if(uPot==NumPotVars){
                    bool UnknownPot=true;
                    for(unsigned uPot=0; uPot<NumPotVars; uPot++){
                        if(POT_VAR==PotFlag[uPot]&&IncludePot[uPot]) {UnknownPot=false; break;}
                    }
                    if(UnknownPot) continue;
                }
                else if(POT_VAR!=PotFlag[uPot]) continue;
                //if(uPot<NumPotVars&&POT_VAR!=PotFlag[uPot]) continue;
                if(uBl==NumBaselines){
                    bool UnknownBl=true;
                    for(unsigned uBl=0; uBl<NumBaselines; uBl++){
                        if(BASELINE_VAR==BlFlag[uBl]&&IncludeBl[uBl]) {UnknownBl=false; break;}
                    }
                    if(UnknownBl) continue;
                }
                else if(BASELINE_VAR!=BlFlag[uBl]) continue;
                if(Sigma0_Feed!=SIGMA_FEED&&SIGMA_FEED>=0&&SIGMA_FEED<=2) continue;
if(SourceAlpha<1.99&&WhichSourceAlpha==0) continue;
else if((SourceAlpha<1.79||SourceAlpha>1.81)&&WhichSourceAlpha==1) continue;
else if((SourceAlpha<1.59||SourceAlpha>1.61)&&WhichSourceAlpha==2) continue;
                NumUsefulEntries++;
                unsigned NumBins = gData->GetN();
                if(UsedNumBins<NumBins) UsedNumBins=NumBins;
                double data_val[2];
                double data_err[2];
                double fit_val[2];
                double mom_val[2];
                double nsigma_val;
                double delta_Ck;
                double delta_k;
                double fit_der;
                double data_der;
                double data_der_err;
                double nsigma_der_val;
                nsigma_0_100=0;
                nsigma_100_200=0;
                nsigma_200_300=0;
                nsigma_0_300=0;
                nsigma_der_0_100=0;
                nsigma_der_100_200=0;
                nsigma_der_200_300=0;
                nsigma_der_0_300=0;
                NDF_0_100=0;
                NDF_100_200=0;
                NDF_200_300=0;
                NDF_0_300=0;
                for(unsigned uBin=0; uBin<NumBins; uBin++){
                    gData->GetPoint(uBin,mom_val[uBin%2],data_val[uBin%2]);
                    if(mom_val[uBin%2]>kLongMax) continue;
                    data_err[uBin%2] = gData->GetErrorY(uBin);
                    fitData->GetPoint(uBin,mom_val[uBin%2],fit_val[uBin%2]);
                    nsigma_val = (fit_val[uBin%2]-data_val[uBin%2])/data_err[uBin%2];
                    //first we compute the chi (and chi2) for the first der
                    fit_der = (fit_val[uBin%2]-fit_val[(uBin+1)%2])/(mom_val[uBin%2]-mom_val[(uBin+1)%2]);
                    data_der = (data_val[uBin%2]-data_val[(uBin+1)%2])/(mom_val[uBin%2]-mom_val[(uBin+1)%2]);
                    data_der_err = sqrt(pow(data_err[uBin%2],2.)+pow(data_err[(uBin+1)%2],2.))/(mom_val[uBin%2]-mom_val[(uBin+1)%2]);
                    nsigma_der_val = (fit_der-data_der)/data_der_err;
//if(nsigma_der_val>10){
//printf("nsigma_der_val=%f at uBin=%u\n",nsigma_der_val,uBin);
//}
                    if(mom_val[uBin%2]<=100){
                        nsigma_0_100 += nsigma_val*nsigma_val;
                        nsigma_der_0_100 += nsigma_der_val*nsigma_der_val;
                        NDF_0_100++;
                    }
                    else if(mom_val[uBin%2]<=200){
                        nsigma_100_200 += nsigma_val*nsigma_val;
                        nsigma_der_100_200 += nsigma_der_val*nsigma_der_val;
                        NDF_100_200++;
                    }
                    else if(mom_val[uBin%2]<=300){
                        nsigma_200_300 += nsigma_val*nsigma_val;
                        nsigma_der_200_300 += nsigma_der_val*nsigma_der_val;
                        NDF_200_300++;
                    }
//CAREFUL, BY DEFAULT 300
                    if(mom_val[uBin%2]<=100){
                        nsigma_0_300 += nsigma_val*nsigma_val;
                        nsigma_der_0_300 += nsigma_der_val*nsigma_der_val;
                        NDF_0_300++;
                    }
                    else break;
                }
//printf("nsigma_0_300=%f\n",nsigma_0_300);
//if(nsigma_der_0_300>500||nsigma_der_0_300<=0||nsigma_der_0_300!=nsigma_der_0_300){
//    printf("nsigma_0_300=%f\n",nsigma_0_300);
//    printf(" nsigma_der_0_300=%f\n",nsigma_der_0_300);
//    usleep(3000e3);
//}
                unsigned WhichBin[5];
                //GetIterCombo(WhichBin,Sigma0_Feed,SourceSize,lam_L_genuine);
                GetIterCombo270520(WhichBin,Sigma0_Feed,SourceSize,lam_L_genuine,CuspWeight,SourceAlpha);
                if(Num_Sigma0_Feed==1) WhichBin[0]=0;
                if(Num_SourceSize==1) WhichBin[1]=0;
                if(Num_lam_L_genuine==1) WhichBin[2]=0;
                if(Num_CuspWeight==1) WhichBin[3]=0;
                if(Num_SourceAlpha==1) WhichBin[4]=0;

                Chi2_300.Add(Chi2_300.GetTotBin(WhichBin),nsigma_0_300);
                NDF_300.Add(NDF_300.GetTotBin(WhichBin),NDF_0_300);
                NumChi2Entries_300.Add(NumChi2Entries_300.GetTotBin(WhichBin),1);

                nsigma_0_300 = sqrt(2)*TMath::ErfcInverse(TMath::Prob(nsigma_0_300,NDF_0_300));
                nsigma_der_0_300 = sqrt(2)*TMath::ErfcInverse(TMath::Prob(nsigma_der_0_300,NDF_0_300));
                htot_nsigma_300->Fill(nsigma_0_300);
                htot_nsigma_der_300->Fill(nsigma_der_0_300);
//if(uPot==6&&uBl==1&&WhichBin[0]==0&&WhichBin[1]==2&&WhichBin[2]==2){
//printf(" nsigma_0_300=%.2f\n",nsigma_0_300);
//}
                //NDF_300.Add(NDF_300.GetTotBin(WhichBin),NDF_0_300);
            }

            //turn the Chi2_300 into a mean Chi2_300
            Chi2_300.DivideHisto(NumChi2Entries_300,false);
            NDF_300.DivideHisto(NumChi2Entries_300,false);
            for(unsigned uBin=0; uBin<NDF_300.GetNbins()-1; uBin++){
                if(NDF_300.GetBinContent(uBin)!=NDF_300.GetBinContent(uBin+1)){
                    ///printf("wtf %f %f\n",NDF_300.GetBinContent(uBin),NDF_300.GetBinContent(uBin+1));
                }
            }


            //find the best solution
            float MinChi2;
            unsigned* Bin_MinChi2;//must be deleted
            unsigned Num_MinChi2 = Chi2_300.FindMinima(MinChi2,Bin_MinChi2);
            int MinNdf = TMath::Nint(NDF_300.GetBinContent(Bin_MinChi2[0]));
            float MinNsigma = sqrt(2)*TMath::ErfcInverse(TMath::Prob(MinChi2,MinNdf));
            //printf("We have %u minima for the chi2=%.2f/%i = %.2f (nsigma = %.2f)\n",Num_MinChi2,MinChi2,MinNdf,MinChi2/float(MinNdf),MinNsigma);
            for(unsigned uMin=0; uMin<Num_MinChi2; uMin++){
                unsigned WhichBin[5];
                printf(" Bin_MinChi2=%p\n",Bin_MinChi2);
                printf(" Bin_MinChi2[0]=%u\n",Bin_MinChi2[0]);
                Chi2_300.GetBinCoordinates(Bin_MinChi2[uMin],WhichBin);
                printf(" Solution %u: S0=%u, SZ=%u, lam=%u, CW=%u, SA=%u\n",uMin,WhichBin[0],WhichBin[1],WhichBin[2],WhichBin[3],WhichBin[4]);

                //WhichBin[0]=0; WhichBin[1]=0; WhichBin[2]=0;
                //printf("             S0=%u, SZ=%u, lam=%u = %.2f\n",WhichBin[0],WhichBin[1],WhichBin[2],Chi2_300.GetBinContent(WhichBin));
                //WhichBin[0]=0; WhichBin[1]=1; WhichBin[2]=0;
                //printf("             S0=%u, SZ=%u, lam=%u = %.2f\n",WhichBin[0],WhichBin[1],WhichBin[2],Chi2_300.GetBinContent(WhichBin));
                //WhichBin[0]=0; WhichBin[1]=2; WhichBin[2]=0;
                //printf("             S0=%u, SZ=%u, lam=%u = %.2f\n",WhichBin[0],WhichBin[1],WhichBin[2],Chi2_300.GetBinContent(WhichBin));
                //WhichBin[0]=0; WhichBin[1]=3; WhichBin[2]=0;
                //printf("             S0=%u, SZ=%u, lam=%u = %.2f\n",WhichBin[0],WhichBin[1],WhichBin[2],Chi2_300.GetBinContent(WhichBin));
            }

            //printf("Proceeding...\n");

/*
            DLM_Sort < double, int > SortTool;
            SortTool.SetData(Chi2_300,NumEntries);
            SortTool.MergeSort();
            SortTool.GetSortedData(Chi2_300,Chi2_300);
            int* Temp = new int[NumEntries];
            for(unsigned uEl=0; uEl<NumEntries; uEl++) Temp[uEl] = NDF_300[SortTool.GetKey()[uEl]];
            for(unsigned uEl=0; uEl<NumEntries; uEl++) NDF_300[uEl] = Temp[uEl];
            delete [] Temp;
            Chi2_300[NumEntries] = 0;
            NDF_300[NumEntries] = 0;
            for(unsigned uBest=0; uBest<NumBestSolutions; uBest++){
                Chi2_300[NumEntries] += Chi2_300[uBest];
                NDF_300[NumEntries] += NDF_300[uBest];
                printf("%u: CHI2=%.2f; NDF=%i; NSIGMA=%.1f\n",uBest,Chi2_300[uBest],NDF_300[uBest],sqrt(2)*TMath::ErfcInverse(TMath::Prob(Chi2_300[uBest],NDF_300[uBest])));
            }
            NSIGMA_300 = sqrt(2)*TMath::ErfcInverse(TMath::Prob(Chi2_300[NumEntries],NDF_300[NumEntries]));
            printf(" CHI2=%.2f; NDF=%i; NSIGMA=%.1f\n",Chi2_300[NumEntries],NDF_300[NumEntries],NSIGMA_300);
*/

            if(!NumUsefulEntries){
                delete hnsigma;
                delete hnsigma_der;
                delete hfdata;
                delete hbl;
                delete hb_nsigma;
                delete hb_nsigma_der;
                delete hb_fdata;
                delete hb_bl;
                delete hnsigma_0_100;
                delete hnsigma_100_200;
                delete hnsigma_200_300;
                delete hnsigma_0_300;
                delete hnsigma_der_0_100;
                delete hnsigma_der_100_200;
                delete hnsigma_der_200_300;
                delete hnsigma_der_0_300;
                delete htot_nsigma_300;
                delete htot_nsigma_der_300;
                delete [] Bin_MinChi2;
                continue;
            }

            double nsigma_300_max;
            double nsigma_300_min;
            //I want to keep the best 68%. => I should get rid of the upper most 32%,
            //a dirty trick here is to just select the 36% central interval, leaving out 32% on both sides of the distribution.
            double nsigma_300_med = GetCentralInterval(*htot_nsigma_300,0.32,nsigma_300_min,nsigma_300_max,true);
            double nsigma_der_300_max;
            double nsigma_der_300_min;
            double nsigma_der_300_med = GetCentralInterval(*htot_nsigma_der_300,0.32,nsigma_der_300_min,nsigma_der_300_max,true);
//nsigma_300_max=nsigma_300_min;
            for(unsigned uEntry=0; uEntry<NumEntries; uEntry++){
                plambdaTree->GetEntry(uEntry);
                //if(uPot<NumPotVars&&POT_VAR!=PotFlag[uPot]) continue;
                //if(uBl<NumBaselines&&BASELINE_VAR!=BlFlag[uBl]) continue;
                if(uPot==NumPotVars){
                    bool UnknownPot=true;
                    for(unsigned uPot=0; uPot<NumPotVars; uPot++){
                        if(POT_VAR==PotFlag[uPot]) {UnknownPot=false; break;}
                    }
                    if(UnknownPot) continue;
                }
                else if(POT_VAR!=PotFlag[uPot]) continue;
                //if(uPot<NumPotVars&&POT_VAR!=PotFlag[uPot]) continue;
                if(uBl==NumBaselines){
                    bool UnknownBl=true;
                    for(unsigned uBl=0; uBl<NumBaselines; uBl++){
                        if(BASELINE_VAR==BlFlag[uBl]) {UnknownBl=false; break;}
                    }
                    if(UnknownBl) continue;
                }
                else if(BASELINE_VAR!=BlFlag[uBl]) continue;
                //if(uBl<NumBaselines&&BASELINE_VAR!=BlFlag[uBl]) continue;
                if(Sigma0_Feed!=SIGMA_FEED&&SIGMA_FEED>=0&&SIGMA_FEED<=2) continue;
if(SourceAlpha<1.99&&WhichSourceAlpha==0) continue;
else if((SourceAlpha<1.79||SourceAlpha>1.81)&&WhichSourceAlpha==1) continue;
else if((SourceAlpha<1.59||SourceAlpha>1.61)&&WhichSourceAlpha==2) continue;
                unsigned NumBins = gData->GetN();
                double data_val[2];
                double data_err[2];
                double fit_val[2];
                double mom_val[2];
                double bl_val;
                double nsigma_val;
                double delta_Ck;
                double delta_k;
                double fit_der;
                double data_der;
                double data_der_err;
                double nsigma_der_val;
                nsigma_0_100=0;
                nsigma_100_200=0;
                nsigma_200_300=0;
                nsigma_0_300=0;
                nsigma_der_0_100=0;
                nsigma_der_100_200=0;
                nsigma_der_200_300=0;
                nsigma_der_0_300=0;
                NDF_0_100=0;
                NDF_100_200=0;
                NDF_200_300=0;
                NDF_0_300=0;
                for(unsigned uBin=0; uBin<NumBins; uBin++){
                    gData->GetPoint(uBin,mom_val[uBin%2],data_val[uBin%2]);
                    if(mom_val[uBin%2]>kLongMax) continue;
                    data_err[uBin%2] = gData->GetErrorY(uBin);
                    fitData->GetPoint(uBin,mom_val[uBin%2],fit_val[uBin%2]);
                    nsigma_val = (fit_val[uBin%2]-data_val[uBin%2])/data_err[uBin%2];
                    //first we compute the chi (and chi2) for the first der
                    fit_der = (fit_val[uBin%2]-fit_val[(uBin+1)%2])/(mom_val[uBin%2]-mom_val[(uBin+1)%2]);
                    data_der = (data_val[uBin%2]-data_val[(uBin+1)%2])/(mom_val[uBin%2]-mom_val[(uBin+1)%2]);
                    data_der_err = sqrt(pow(data_err[uBin%2],2.)+pow(data_err[(uBin+1)%2],2.))/(mom_val[uBin%2]-mom_val[(uBin+1)%2]);
                    nsigma_der_val = (fit_der-data_der)/data_der_err;
                    if(mom_val[uBin%2]<=100){
                        nsigma_0_100 += nsigma_val*nsigma_val;
                        nsigma_der_0_100 += nsigma_der_val*nsigma_der_val;
                        NDF_0_100++;
                    }
                    else if(mom_val[uBin%2]<=200){
                        nsigma_100_200 += nsigma_val*nsigma_val;
                        nsigma_der_100_200 += nsigma_der_val*nsigma_der_val;
                        NDF_100_200++;
                    }
                    else if(mom_val[uBin%2]<=300){
                        nsigma_200_300 += nsigma_val*nsigma_val;
                        nsigma_der_200_300 += nsigma_der_val*nsigma_der_val;
                        NDF_200_300++;
                    }
                    if(mom_val[uBin%2]<=300){
                        nsigma_0_300 += nsigma_val*nsigma_val;
                        nsigma_der_0_300 += nsigma_der_val*nsigma_der_val;
                        NDF_0_300++;
                    }
                    else break;
                }
                nsigma_0_300 = sqrt(2)*TMath::ErfcInverse(TMath::Prob(nsigma_0_300,NDF_0_300));
                htot_nsigma_300->Fill(nsigma_0_300);
                nsigma_0_100 = sqrt(2)*TMath::ErfcInverse(TMath::Prob(nsigma_0_100,NDF_0_100));
                nsigma_100_200 = sqrt(2)*TMath::ErfcInverse(TMath::Prob(nsigma_100_200,NDF_100_200));
                nsigma_200_300 = sqrt(2)*TMath::ErfcInverse(TMath::Prob(nsigma_200_300,NDF_200_300));
//if(nsigma_0_300<1.5){
//printf("upot ubl = %u %u\n",uPot,uBl);
//printf(" POT_VAR BASELINE_VAR = %u %u\n",POT_VAR,BASELINE_VAR);
//}

                nsigma_der_0_300 = sqrt(2)*TMath::ErfcInverse(TMath::Prob(nsigma_der_0_300,NDF_0_300));
                htot_nsigma_der_300->Fill(nsigma_der_0_300);
                nsigma_der_0_100 = sqrt(2)*TMath::ErfcInverse(TMath::Prob(nsigma_der_0_100,NDF_0_100));
                nsigma_der_100_200 = sqrt(2)*TMath::ErfcInverse(TMath::Prob(nsigma_der_100_200,NDF_100_200));
                nsigma_der_200_300 = sqrt(2)*TMath::ErfcInverse(TMath::Prob(nsigma_der_200_300,NDF_200_300));

                //if(nsigma_0_300>nsigma_300_max) continue;
//if(!DefaultVariation) continue;
//if(nsigma_200_300>1.0) continue;
//if(SourceSize>1.03||SourceSize<1.01) continue;
//if(lam_L_genuine>0.46||SourceSize<0.45) continue;
//if(Sigma0_Feed==0) continue;
//if(nsigma_0_100>4) continue;
                for(unsigned uBin=0; uBin<NumBins; uBin++){
                    gData->GetPoint(uBin,mom_val[uBin%2],data_val[uBin%2]);
                    if(mom_val[uBin%2]>kLongMax) continue;
                    data_err[uBin%2] = gData->GetErrorY(uBin);
                    fitData->GetPoint(uBin,mom_val[uBin%2],fit_val[uBin%2]);
                    fitBaseline->GetPoint(uBin,mom_val[uBin%2],bl_val);
                    nsigma_val = (fit_val[uBin%2]-data_val[uBin%2])/data_err[uBin%2];
                    hnsigma->SetBinContent(uBin+1,hnsigma->GetYaxis()->FindBin(nsigma_val),hnsigma->GetBinContent(uBin+1)+1);
                    hnsigma_der->SetBinContent(uBin+1,hnsigma_der->GetYaxis()->FindBin(nsigma_val),hnsigma_der->GetBinContent(uBin+1)+1);
                    hfdata->SetBinContent(uBin+1,hfdata->GetYaxis()->FindBin(fit_val[uBin%2]),hfdata->GetBinContent(uBin+1)+1);
                    hbl->SetBinContent(uBin+1,hbl->GetYaxis()->FindBin(bl_val),hbl->GetBinContent(uBin+1)+1);
                    unsigned WhichBinCurrent[3];
                    GetIterCombo(WhichBinCurrent,Sigma0_Feed,SourceSize,lam_L_genuine);
                    bool OneOfTheBest=false;
                    for(unsigned uMin=0; uMin<Num_MinChi2; uMin++){
                        unsigned WhichBinBest[3];
                        Chi2_300.GetBinCoordinates(Bin_MinChi2[uMin],WhichBinBest);
                        if(WhichBinBest[0]==WhichBinCurrent[0]&&
                           WhichBinBest[1]==WhichBinCurrent[1]&&
                           WhichBinBest[2]==WhichBinCurrent[2]) OneOfTheBest=true;
                        //printf(" Solution %u: S0=%u, SZ=%u, lam=%u\n",uMin,WhichBinBest[0],WhichBinBest[1],WhichBinBest[2]);
                    }
                    if(OneOfTheBest){
                        hb_nsigma->SetBinContent(uBin+1,hnsigma->GetYaxis()->FindBin(nsigma_val),hnsigma->GetBinContent(uBin+1)+1);
                        hb_nsigma_der->SetBinContent(uBin+1,hnsigma_der->GetYaxis()->FindBin(nsigma_val),hnsigma_der->GetBinContent(uBin+1)+1);
                        hb_fdata->SetBinContent(uBin+1,hfdata->GetYaxis()->FindBin(fit_val[uBin%2]),hfdata->GetBinContent(uBin+1)+1);
                        hb_bl->SetBinContent(uBin+1,hbl->GetYaxis()->FindBin(bl_val),hbl->GetBinContent(uBin+1)+1);
                    }
                }
                hnsigma_0_100->Fill(nsigma_0_100);
                hnsigma_100_200->Fill(nsigma_100_200);
                hnsigma_200_300->Fill(nsigma_200_300);
                hnsigma_0_300->Fill(nsigma_0_300);

                hnsigma_der_0_100->Fill(nsigma_der_0_100);
                hnsigma_der_100_200->Fill(nsigma_der_100_200);
                hnsigma_der_200_300->Fill(nsigma_der_200_300);
                hnsigma_der_0_300->Fill(nsigma_der_0_300);

            }

            DLM_CommonAnaFunctions AnalysisObject;
            AnalysisObject.SetCatsFilesFolder("/mnt/Ubuntu_Data/CernBox/Sync/CatsFiles");
            TH1F* hData_pL_Stat = AnalysisObject.GetAliceExpCorrFun(DataSample[0],"pLambda","_0",2,false,-1);

            TGraphErrors* ge_Bl = new TGraphErrors();
            ge_Bl->SetName("ge_Bl");
            TGraphErrors* ge_Fit = new TGraphErrors();
            ge_Fit->SetName("ge_Fit");
            TGraphErrors* ge_SigmaSlice = new TGraphErrors();
            ge_SigmaSlice->SetName("ge_SigmaSlice");
            TGraphErrors* ge_SigmaSlice_der = new TGraphErrors();
            ge_SigmaSlice_der->SetName("ge_SigmaSlice_der");

            TGraphErrors* geb_Bl = new TGraphErrors();
            geb_Bl->SetName("geb_Bl");
            TGraphErrors* geb_Fit = new TGraphErrors();
            geb_Fit->SetName("geb_Fit");
            TGraphErrors* geb_SigmaSlice = new TGraphErrors();
            geb_SigmaSlice->SetName("geb_SigmaSlice");
            TGraphErrors* geb_SigmaSlice_der = new TGraphErrors();
            geb_SigmaSlice_der->SetName("geb_SigmaSlice_der");

            TF1* fe_Sigma = new TF1("fe_Sigma","[0]",0,300);

            fe_Sigma->SetName("fe_Sigma");
            for(unsigned uBin=0; uBin<UsedNumBins; uBin++){
                double MOM = hData_pL_Stat->GetBinCenter(uBin+1);
                double MinVal,MaxVal,MedVal;

                TH1* hProj;

                hProj = hfdata->ProjectionY("hProj",uBin+1,uBin+1);
                if(!hProj->Integral()) continue;
                MedVal = GetCentralInterval(*hProj,0.68,MinVal,MaxVal,true);
                ge_Fit->SetPoint(uBin,MOM,(MaxVal+MinVal)*0.5);
                ge_Fit->SetPointError(uBin,0,(MaxVal-MinVal)*0.5);
                delete hProj;

                hProj = hbl->ProjectionY("hProj",uBin+1,uBin+1);
                MedVal = GetCentralInterval(*hProj,0.68,MinVal,MaxVal,true);
                ge_Bl->SetPoint(uBin,MOM,(MaxVal+MinVal)*0.5);
                ge_Bl->SetPointError(uBin,0,(MaxVal-MinVal)*0.5);
                delete hProj;

                hProj = hnsigma->ProjectionY("hProj",uBin+1,uBin+1);
                MedVal = GetCentralInterval(*hProj,0.68,MinVal,MaxVal,true);
                ge_SigmaSlice->SetPoint(uBin,MOM,(MaxVal+MinVal)*0.5);
                ge_SigmaSlice->SetPointError(uBin,0,(MaxVal-MinVal)*0.5);
                delete hProj;

                hProj = hnsigma_der->ProjectionY("hProj",uBin+1,uBin+1);
                MedVal = GetCentralInterval(*hProj,0.68,MinVal,MaxVal,true);
                ge_SigmaSlice_der->SetPoint(uBin,MOM,(MaxVal+MinVal)*0.5);
                ge_SigmaSlice_der->SetPointError(uBin,0,(MaxVal-MinVal)*0.5);
                delete hProj;

                hProj = hb_fdata->ProjectionY("hProj",uBin+1,uBin+1);
                if(!hProj->Integral()) continue;
                MedVal = GetCentralInterval(*hProj,0.68,MinVal,MaxVal,true);
                geb_Fit->SetPoint(uBin,MOM,(MaxVal+MinVal)*0.5);
                geb_Fit->SetPointError(uBin,0,(MaxVal-MinVal)*0.5);
                delete hProj;

                hProj = hb_bl->ProjectionY("hProj",uBin+1,uBin+1);
                MedVal = GetCentralInterval(*hProj,0.68,MinVal,MaxVal,true);
                geb_Bl->SetPoint(uBin,MOM,(MaxVal+MinVal)*0.5);
                geb_Bl->SetPointError(uBin,0,(MaxVal-MinVal)*0.5);
                delete hProj;

                hProj = hb_nsigma->ProjectionY("hProj",uBin+1,uBin+1);
                MedVal = GetCentralInterval(*hProj,0.68,MinVal,MaxVal,true);
                geb_SigmaSlice->SetPoint(uBin,MOM,(MaxVal+MinVal)*0.5);
                geb_SigmaSlice->SetPointError(uBin,0,(MaxVal-MinVal)*0.5);
                delete hProj;

                hProj = hb_nsigma_der->ProjectionY("hProj",uBin+1,uBin+1);
                MedVal = GetCentralInterval(*hProj,0.68,MinVal,MaxVal,true);
                geb_SigmaSlice_der->SetPoint(uBin,MOM,(MaxVal+MinVal)*0.5);
                geb_SigmaSlice_der->SetPointError(uBin,0,(MaxVal-MinVal)*0.5);
                delete hProj;

                fe_Sigma->FixParameter(0,MinNsigma);

/*
                if(MOM<=100){
                    MedVal = GetCentralInterval(*hnsigma_0_100,0.68,MinVal,MaxVal,true);
//MaxVal = MinVal+0.25;
                    ge_SigmaSlice->SetPoint(uBin,MOM,(MaxVal+MinVal)*0.5);
                    ge_SigmaSlice->SetPointError(uBin,0,(MaxVal-MinVal)*0.5);

                    MedVal = GetCentralInterval(*hnsigma_der_0_100,1.-1e-6,MinVal,MaxVal,true);
//MaxVal = MinVal+0.25;
                    ge_SigmaSlice_der->SetPoint(uBin,MOM,(MaxVal+MinVal)*0.5);
                    ge_SigmaSlice_der->SetPointError(uBin,0,(MaxVal-MinVal)*0.5);
                }
                else if(MOM<=200){
                    MedVal = GetCentralInterval(*hnsigma_100_200,1.-1e-6,MinVal,MaxVal,true);
//MaxVal = MinVal+0.25;
                    ge_SigmaSlice->SetPoint(uBin,MOM,(MaxVal+MinVal)*0.5);
                    ge_SigmaSlice->SetPointError(uBin,0,(MaxVal-MinVal)*0.5);

                    MedVal = GetCentralInterval(*hnsigma_der_100_200,1.-1e-6,MinVal,MaxVal,true);
//MaxVal = MinVal+0.25;
                    ge_SigmaSlice_der->SetPoint(uBin,MOM,(MaxVal+MinVal)*0.5);
                    ge_SigmaSlice_der->SetPointError(uBin,0,(MaxVal-MinVal)*0.5);
                }
                else if(MOM<=300){
                    MedVal = GetCentralInterval(*hnsigma_200_300,1.-1e-6,MinVal,MaxVal,true);
//MaxVal = MinVal+0.25;
                    ge_SigmaSlice->SetPoint(uBin,MOM,(MaxVal+MinVal)*0.5);
                    ge_SigmaSlice->SetPointError(uBin,0,(MaxVal-MinVal)*0.5);

                    MedVal = GetCentralInterval(*hnsigma_der_200_300,1.-1e-6,MinVal,MaxVal,true);
//MaxVal = MinVal+0.25;
                    ge_SigmaSlice_der->SetPoint(uBin,MOM,(MaxVal+MinVal)*0.5);
                    ge_SigmaSlice_der->SetPointError(uBin,0,(MaxVal-MinVal)*0.5);
                }
                if(MOM<=300){
                    MedVal = GetCentralInterval(*hnsigma_0_300,1.-1e-6,MinVal,MaxVal,true);
                    fe_Sigma->FixParameter(0,MinVal);
                    //fe_Sigma->FixParameter(0,NSIGMA_300);
//printf("MinVal = %f\n",MinVal);
                }
                delete hProj;
*/
            }


            gStyle->SetCanvasPreferGL(1);
            SetStyle();

            ge_Fit->SetFillColorAlpha(kRed+1,0.40);
            ge_Fit->SetLineColor(kRed+1);
            ge_Fit->SetLineWidth(5);

            ge_Bl->SetFillColorAlpha(kCyan+1,0.40);
            ge_Bl->SetLineColor(kCyan+1);
            ge_Bl->SetLineWidth(5);

            ge_SigmaSlice->SetFillColorAlpha(kRed+1,0.40);
            ge_SigmaSlice_der->SetFillColorAlpha(kBlue+1,0.40);

            geb_Fit->SetFillColorAlpha(kRed+1,0.70);
            geb_Fit->SetLineColor(kRed+1);
            geb_Fit->SetLineWidth(3);

            geb_Bl->SetFillColorAlpha(kCyan+1,0.70);
            geb_Bl->SetLineColor(kCyan+1);
            geb_Bl->SetLineWidth(5);

            geb_SigmaSlice->SetFillColorAlpha(kRed+1,0.70);
            geb_SigmaSlice_der->SetFillColorAlpha(kBlue+1,0.70);

            fe_Sigma->SetLineColor(kBlack);
            fe_Sigma->SetLineWidth(3);
            fe_Sigma->SetLineStyle(2);

            DLM_SubPads* DlmPad = new DLM_SubPads(720,720);
            DlmPad->AddSubPad(0,1,0.33,1);
            DlmPad->AddSubPad(0,1,0,0.33);
            DlmPad->SetMargin(0,0.12,0.02,0.0,0.02);
            DlmPad->SetMargin(1,0.12,0.02,0.09,0.0);
            DlmPad->cd(0);

            hData_pL_Stat->SetTitle("; #it{k*} (MeV/#it{c}); #it{C}(#it{k*})");
            hData_pL_Stat->GetXaxis()->SetRangeUser(0, 456);
            hData_pL_Stat->GetXaxis()->SetNdivisions(505);
            hData_pL_Stat->GetYaxis()->SetRangeUser(0.9, 2.0);
            hData_pL_Stat->SetFillColor(kGray+1);
            SetStyleHisto2(hData_pL_Stat,2,0);
            //hData_pL_Stat->GetYaxis()->SetTitleOffset(1.0);
            hData_pL_Stat->Draw();

            //TFile* DataSystFile = new TFile("/mnt/Ubuntu_Data/CernBox/Sync/CatsFiles/ExpData/ALICE_pp_13TeV_HM/Sample10HM/Systematics_pL.root");
            //TH1F* DataSystHisto = (TH1F*)DataSystFile->Get("SystErrRel");
            //TGraphErrors *Tgraph_syserror = DrawSystematicError_FAST(hData_pL_Stat, DataSystHisto, NULL, 3);
            TGraphErrors *Tgraph_syserror = (TGraphErrors*)gdata_toterr->Clone("Tgraph_syserror");
            for(unsigned uBin=0; uBin<Tgraph_syserror->GetN(); uBin++){
                double tot_err,stat_err,syst_err,k_err,Ck_val;
                double xVal,yVal;
                Ck_val = hData_pL_Stat->GetBinContent(uBin+1);
                k_err = Tgraph_syserror->GetErrorX(uBin);
                tot_err = Tgraph_syserror->GetErrorY(uBin);
                stat_err = hData_pL_Stat->GetBinError(uBin+1);
                syst_err = tot_err>stat_err?sqrt(tot_err*tot_err-stat_err*stat_err):0;
                Tgraph_syserror->GetPoint(uBin,xVal,yVal);
                Tgraph_syserror->SetPoint(uBin,xVal,Ck_val);
                Tgraph_syserror->SetPointError(uBin,k_err,syst_err);
            }
            Tgraph_syserror->SetLineColor(kWhite);
            //baselineLL->Draw("same");
            //if(!DataOnly){
            ge_Bl->Draw("3 same");
            ge_Fit->Draw("3 same");

            geb_Bl->Draw("3 same");
            geb_Fit->Draw("3L same");
            //}
            hData_pL_Stat->Draw("same");

            Tgraph_syserror->SetFillColorAlpha(kBlack, 0.4);
            Tgraph_syserror->Draw("2 same");
            //hData->Draw("pe same");

            TString LegendSource_line1 = SourceName1[0];

            unsigned NumRows=4;
            TGraph DummyLegendEntry;
            DummyLegendEntry.SetName("DummyLegendEntry");
            DummyLegendEntry.SetLineColor(kWhite);
            DummyLegendEntry.SetMarkerColor(kWhite);
            TGraph DummyLegendEntry2;
            DummyLegendEntry2.SetName("DummyLegendEntry2");
            DummyLegendEntry2.SetLineColor(kWhite);
            DummyLegendEntry2.SetMarkerColor(kWhite);
            TLegend *legend = new TLegend(0.39,0.73-0.055*NumRows,0.73,0.78);//lbrt
            legend->SetBorderSize(0);
            legend->SetTextFont(42);
            legend->SetTextSize(gStyle->GetTextSize()*0.90);
            TH1F* hCk_Fake;
            hCk_Fake = (TH1F*)hData_pL_Stat->Clone("hCk_Fake");
            hCk_Fake->SetName("hCk_Fake");
            hCk_Fake->SetLineColor(hCk_Fake->GetFillColor());

            legend->AddEntry(hCk_Fake, "p#minus #Lambda #oplus #bar{p}#minus #bar{#Lambda} pairs", "fpe");
            //if(!DataOnly){
            legend->AddEntry(ge_Fit,"Femtoscopic fit "+PotName1[uPot],"l");
            legend->AddEntry(&DummyLegendEntry,PotName2[uPot],"l");
            legend->AddEntry(ge_Bl,BlName1[uBl],"l");
            legend->AddEntry(&DummyLegendEntry2,BlName2[uBl],"l");
            //}
            legend->Draw("same");
            TLatex BeamText;
            BeamText.SetTextSize(gStyle->GetTextSize()*0.90);
            BeamText.SetNDC(kTRUE);
            BeamText.DrawLatex(0.40, 0.915, "ALICE work in progress");
            BeamText.DrawLatex(0.40, 0.860, "high-mult. (0#minus 0.17% INEL>0) pp #sqrt{#it{s}} = 13 TeV");

            TLatex BeamTextSource;
            BeamTextSource.SetTextSize(gStyle->GetTextSize()*0.90);
            BeamTextSource.SetNDC(kTRUE);
            BeamTextSource.DrawLatex(0.40, 0.805, LegendSource_line1);

        //INLET -------------------------------------------------------------------------------------------------------------------

            TH1F* DataHisto_Inlet = (TH1F*)hData_pL_Stat->Clone("DataHisto_Inlet");
            DataHisto_Inlet->SetMarkerSize(hData_pL_Stat->GetMarkerSize()*0.67);
            DataHisto_Inlet->SetLineWidth(hData_pL_Stat->GetLineWidth()*0.67);
            DataHisto_Inlet->GetXaxis()->SetTitleSize(hData_pL_Stat->GetXaxis()->GetTitleSize()*1.75);
            DataHisto_Inlet->GetXaxis()->SetLabelSize(hData_pL_Stat->GetXaxis()->GetLabelSize()*1.75);
            DataHisto_Inlet->GetXaxis()->SetRangeUser(0, 456);
            DataHisto_Inlet->GetXaxis()->SetNdivisions(505);

            DataHisto_Inlet->GetYaxis()->SetTitleSize(hData_pL_Stat->GetYaxis()->GetTitleSize()*1.75);
            DataHisto_Inlet->GetYaxis()->SetLabelSize(hData_pL_Stat->GetYaxis()->GetLabelSize()*1.75);
            DataHisto_Inlet->GetYaxis()->SetTitleOffset(hData_pL_Stat->GetYaxis()->GetTitleOffset()*0.67);
            DataHisto_Inlet->GetYaxis()->SetRangeUser(0.98, 1.05);

            TGraph* grFemto_Inlet = (TGraph*)ge_Fit->Clone("grFemto_Inlet");
            grFemto_Inlet->SetLineWidth(ge_Fit->GetLineWidth()*0.67);

            TGraph* grbFemto_Inlet = (TGraph*)geb_Fit->Clone("grbFemto_Inlet");
            grbFemto_Inlet->SetLineWidth(geb_Fit->GetLineWidth()*0.67);

            const double fXMinInlet=0.17;
            const double fYMinInlet=0.12;
            const double fXMaxInlet=0.97;
            const double fYMaxInlet=0.50;
            TPad *inset_pad = new TPad("insert", "insertPad", fXMinInlet, fYMinInlet,
                                     fXMaxInlet, fYMaxInlet);
            inset_pad->SetTopMargin(0.01);
            inset_pad->SetRightMargin(0.05);
            inset_pad->SetBottomMargin(0.28);
            inset_pad->SetLeftMargin(0.28);
            inset_pad->SetFillStyle(4000);
            //if(!DataOnly)
            inset_pad->Draw();
            inset_pad->cd();
            DataHisto_Inlet->Draw();
            //if(!DataOnly){
            ge_Bl->Draw("3 same");
            grFemto_Inlet->Draw("3 same");
            geb_Bl->Draw("3 same");
            grbFemto_Inlet->Draw("3 same");
            //}
            DataHisto_Inlet->Draw("same");
            Tgraph_syserror->Draw("2 same");

            DlmPad->cd(1);
            TH1F* hAxis = new TH1F("hAxis", "hAxis", 456/12, 0, 456);
            hAxis->SetStats(false);
            //if(!DataOnly)
            hAxis->SetTitle("; #it{k*} (MeV/#it{c}); #it{n_{#sigma}}");
            //else //hAxis->SetTitle("; #it{k*} (MeV/#it{c}); #it{C}(#it{k*})");
            hAxis->GetXaxis()->SetRangeUser(0, 456);
            //if(!DataOnly)
            hAxis->GetYaxis()->SetRangeUser(-6.2,6.2);
            hAxis->GetXaxis()->SetNdivisions(505);
            //else hAxis->GetYaxis()->SetRangeUser(0.98, 1.05);
            //hData->SetTitle("; #it{k*} (MeV/#it{c}); #it{C}(#it{k*})");
            //hData->GetXaxis()->SetRangeUser(0, 320);
            //hData->GetYaxis()->SetRangeUser(0.85, 2.3);
            //hData->SetFillColor(fFillColors[0]);
            SetStyleHisto2(hAxis,2,0,2);
            hAxis->GetYaxis()->SetNdivisions(504);
            //hData->GetYaxis()->SetTitleOffset(1.0);
            hAxis->Draw("");

            //ge_SigmaSlice->SetFillColorAlpha(kRed+1,0.67);
            //ge_SigmaSlice->GetYaxis()->SetNdivisions(204);
            //if(!DataOnly)
            ge_SigmaSlice->Draw("3 same");
            geb_SigmaSlice->Draw("3 same");
            /// a plot of the deviation of the first derivative
            //ge_SigmaSlice_der->Draw("3 same");
            //fe_Sigma->Draw("same");
            //else{
            //hData_pL_Stat->Draw("same");
            //Tgraph_syserror->Draw("2 same");
            //}

            TLatex BeamText_nsigma;
            BeamText_nsigma.SetTextSize(gStyle->GetTextSize()*1.8);
            BeamText_nsigma.SetNDC(kTRUE);
            //BeamText_nsigma.DrawLatex(0.80, 0.90, TString::Format("n_{#sigma} = %.1f",MinNsigma));
            if(MinChi2>=99.5) BeamText_nsigma.DrawLatex(0.673, 0.90, TString::Format("#chi^{2}/NDF = %.0f/%i = %.1f",MinChi2,MinNdf,MinChi2/double(MinNdf)));
            else BeamText_nsigma.DrawLatex(0.688, 0.90, TString::Format("#chi^{2}/NDF = %.0f/%i = %.1f",MinChi2,MinNdf,MinChi2/double(MinNdf)));
            TLatex BeamText_nsigma_info;
            BeamText_nsigma_info.SetTextSize(gStyle->GetTextSize()*1.2);
            BeamText_nsigma_info.SetNDC(kTRUE);
            //BeamText_nsigma_info.DrawLatex(0.80, 0.83, TString::Format("k*#in[0, 300] MeV/#it{c}"));
            BeamText_nsigma_info.DrawLatex(0.73, 0.83, TString::Format("k*#in[0, 300] MeV/#it{c}, n_{#sigma} = %.1f",MinNsigma));

            DlmPad->GetCanvas()->SaveAs(OutputFolder+TString::Format("DlmPad_%s_%s.pdf",PotDescr[uPot].Data(),BlDescr[uBl].Data()));


//TCanvas* c_hb_bl = new TCanvas("c_hb_bl","c_hb_bl",1)
//c_hb_bl->cd(0); c_hb_bl->SetCanvasSize(1280, 720); c_hb_bl->SetMargin(0.15,0.05,0.2,0.05);//lrbt
//hb_bl->Draw("COLZ");

            delete hnsigma;
            delete hb_nsigma;
            delete hnsigma_0_300;
            delete hnsigma_0_100;
            delete hnsigma_100_200;
            delete hnsigma_200_300;
            delete htot_nsigma_300;
            delete hnsigma_der;
            delete hb_nsigma_der;
            delete hnsigma_der_0_300;
            delete hnsigma_der_0_100;
            delete hnsigma_der_100_200;
            delete hnsigma_der_200_300;
            delete htot_nsigma_der_300;
            delete hfdata;
            delete hbl;
            delete hb_fdata;
            delete hb_bl;

            delete ge_SigmaSlice;
            delete ge_SigmaSlice_der;
            delete geb_SigmaSlice;
            delete geb_SigmaSlice_der;
            delete fe_Sigma;

            delete hData_pL_Stat;
            delete ge_Bl;
            delete ge_Fit;
            delete geb_Bl;
            delete geb_Fit;
            delete hAxis;
            delete hCk_Fake;
            delete legend;
            delete DataHisto_Inlet;
            delete grFemto_Inlet;
            delete grbFemto_Inlet;
            delete inset_pad;
            delete DlmPad;

            delete [] Bin_MinChi2;
        }
    }

    //delete [] Chi2_300;
    //delete [] NDF_300;
    //delete [] NSIGMA_300;
    delete [] IncludePot;
    delete [] IncludeBl;
    delete [] BlFlag;
    delete [] BlName1;
    delete [] BlName2;
    delete [] BlDescr;
    delete [] PotName1;
    delete [] PotName2;
    delete [] PotFlag;
    delete [] PotDescr;
    delete [] SourceName1;
    //delete [] Chi2Ndf;
    delete hdata_toterr;
    delete gdata_toterr;
    delete plambdaTree;


    delete InputFile;
}

//the files are separated for everything BUT the SourceAlpha
void Plot_pL_SystematicsMay2020_2(const int& SIGMA_FEED,
                                  const int& WhichBaseline,
                                  const int& WhichPotential,
                                  const float& ValSourceAlpha,
                                  TString InputFolder, TString InputFileName, TString OutputFolder, const int& WhichDataSet=-1){
//printf("debug\n");
    //if 0, we have much more info on the plots
    //if 1: as intended for the paper (cleaner)
    //if 2: as intended for the thesis
    const int PlotsType = 1;
    int Panel_X = 0;
    if(PlotsType==1){
      if(WhichPotential==11600 && SIGMA_FEED==1)Panel_X=0;
      if(WhichPotential==11600 && SIGMA_FEED==0)Panel_X=1;
      if(WhichPotential==1600 && SIGMA_FEED==1)Panel_X=2;
      if(WhichPotential==1600 && SIGMA_FEED==0)Panel_X=3;
    }
    const bool WorkInProgress = false;
    const bool DataOnly = false;

    enum BLTYPE { pol0s,pol1s,pol2s,pol3s,dpol2s,dpol3s,dpol4s,pol2e,pol3e,dpol2e,dpol3e,dpol4e,spl1 };
    TString BlName1;
    TString BlName2;
    TString BlDescr;
    const float MinRad = PlotsType?0.96:0.96;//0.96
    const float MaxRad = PlotsType?1.08:1.08;
    const float MinOmega = PlotsType?0.35:0.35;//0.25
//const float MinOmega = PlotsType?0.25:0.25;//0.25
    const float MaxOmega = PlotsType?0.45:0.45;
    const float MinCkConv = PlotsType?600:600;
    const float MaxCkConv = PlotsType?800:800;
    //const float Use_SBpur = 0.963;
    const float Use_SBpur = PlotsType?0.963:0.0;
    const int MinGoodCut = PlotsType?-100:-100;
    const int MaxGoodCut = PlotsType?100:100;
    const float FractionOfSolutions = PlotsType?1.:1.;

    const bool Same_omega_siglam = PlotsType?false:false;
    //const bool COMPARE_TO_LO = true;
    bool COMPARE_TO_LO = true;
    if(PlotsType&&WhichPotential==1600) COMPARE_TO_LO = true;
    else if(PlotsType) COMPARE_TO_LO = false;
    if(WhichPotential==-11600) COMPARE_TO_LO = false;
//COMPARE_TO_LO = false;
    //3 digits  #1 = plot the sum of both
    //          #2 = plot xim
    //          #3 = plot xi0
    const int XI_PLOT = PlotsType?100:100;
    const bool PlotLines = PlotsType?false:false;
    const bool DefaultBinValue = PlotsType?true:true;
    int* GoodCutVar = new int [45];
    //for(unsigned uCut=0; uCut<45; uCut++) GoodCutVar[uCut] = 25;
    //exclude all (12) cut vars with an yield change of more than 10%
    GoodCutVar[0] = 0;
    GoodCutVar[1] = -15;
    GoodCutVar[2] = -15;
    GoodCutVar[3] = -10;
    GoodCutVar[4] = -10;
    GoodCutVar[5] = 5;
    GoodCutVar[6] = -10;
    GoodCutVar[7] = -5;
    GoodCutVar[8] = 10;
    GoodCutVar[9] = -10;
    GoodCutVar[10] = 5;
    GoodCutVar[11] = -15;
    GoodCutVar[12] = 5;
    GoodCutVar[13] = -5;
    GoodCutVar[14] = 5;
    GoodCutVar[15] = -25;
    GoodCutVar[16] = 5;
    GoodCutVar[17] = -5;
    GoodCutVar[18] = -1;
    GoodCutVar[19] = -15;
    GoodCutVar[20] = 10;
    GoodCutVar[21] = -5;
    GoodCutVar[22] = 5;
    GoodCutVar[23] = -10;
    GoodCutVar[24] = -5;
    GoodCutVar[25] = 5;
    GoodCutVar[26] = -20;
    GoodCutVar[27] = -10;
    GoodCutVar[28] = -10;
    GoodCutVar[29] = -5;
    GoodCutVar[30] = -5;
    GoodCutVar[31] = 10;
    GoodCutVar[32] = -20;
    GoodCutVar[33] = -15;
    GoodCutVar[34] = -10;
    GoodCutVar[35] = -5;
    GoodCutVar[36] = -15;
    GoodCutVar[37] = -1;
    GoodCutVar[38] = -15;
    GoodCutVar[39] = 1;
    GoodCutVar[40] = -15;
    GoodCutVar[41] = -5;
    GoodCutVar[42] = -25;
    GoodCutVar[43] = -10;
    GoodCutVar[44] = -5;
    if(WhichDataSet>=0){
        for(unsigned uCut=0; uCut<45; uCut++){
            GoodCutVar[uCut] = 1000000;
        }
        GoodCutVar[WhichDataSet] = MinGoodCut;
    }


    switch(WhichBaseline){
        case pol0s :    BlName1 = "Constant baseline";
                        BlName2 = "Standard fit range";
                        BlDescr = "pol0s";
                        break;
        case pol1s :    BlName1 = "Linear baseline";
                        BlName2 = "Standard fit range";
                        BlDescr = "pol1s";
                        break;
        case pol2s :    BlName1 = "Quadratic baseline";
                        BlName2 = "Standard fit range";
                        BlDescr = "pol2s";
                        break;
        case pol3s :    BlName1 = "Cubic baseline";
                        BlName2 = "Standard fit range";
                        BlDescr = "pol3s";
                        break;
        case dpol2s :   BlName1 = "Quadratic baseline";
                        BlName2 = "Standard fit range";
                        BlDescr = "dpol2s";
                        break;
        case dpol3s :   BlName1 = "Cubic baseline";
                        BlName2 = "Standard fit range";
                        BlDescr = "dpol3s";
                        break;
        case dpol4s :   BlName1 = "Quartic baseline";
                        BlName2 = "Standard fit range";
                        BlDescr = "dpol4s";
                        break;
        case pol2e :    BlName1 = "Quadratic baseline";
                        BlName2 = "Extended fit range";
                        BlDescr = "pol2e";
                        break;//
        case pol3e :    BlName1 = "Cubic baseline";
                        BlName2 = "Extended fit range";
                        BlDescr = "pol3e";
                        break;
        case dpol2e :   BlName1 = "Quadratic baseline";
                        BlName2 = "Extended fit range";
                        BlDescr = "dpol2e";
                        break;
        case dpol3e :   BlName1 = "Cubic baseline";
                        BlName2 = "Extended fit range";
                        BlDescr = "dpol3e";
                        break;
        case dpol4e :   BlName1 = "Quartic baseline";
                        BlName2 = "Extended fit range";
                        BlDescr = "dpol4e";
                        break;
        case spl1   :   BlName1 = "Spline baseline (1)";
                        BlName2 = "Extended fit range";
                        BlDescr = "spl1";
                        break;
        default :       BlName1 = "Unknown baseline";
                        BlName2 = "";
                        BlDescr = "Unknown";
                        break;
    }

    TString PotName1;
    TString PotName2;
    TString PotName3;
    TString PotDescr;
    TString SigDescr;
    TString SigName1;
    TString XimName1 = "Residual p#minus #Xi^{#minus}: HAL QCD";
    TString Xi0Name1 = "Residual p#minus #Xi^{0}: HAL QCD";
    TString XiName1 = "Residual p#minus #Xi^{#minus}#oplusp#minus #Xi^{0}: HAL QCD";
    if(PlotsType==1) XiName1 = "Residual p#minus #Xi^{#minus}#oplusp#minus #Xi^{0}";
    int PotFlag;
    switch(WhichPotential){
        case 1500 : PotName1 = "NLO13 (500)";
                    PotName2 = "s,d waves";
                    PotDescr = "NLO13-500";
                    PotFlag = 1500;
                    break;
        case 1550 : PotName1 = "NLO13 (550)";
                    PotName2 = "s,d waves";
                    PotDescr = "NLO13-550";
                    PotFlag = 1550;
                    break;
        case 1600 : PotName1 = "NLO13 (600)";
                    PotName2 = "s,d waves";
                    PotDescr = "NLO13-600";
                    PotFlag = 1600;
                    break;
        case 1650 : PotName1 = "NLO13 (650)";
                    PotName2 = "s,d waves";
                    PotDescr = "NLO13-650";
                    PotFlag = 1650;
                    break;
        case 11500 :PotName1 = "NLO19 (500)";
                    PotName2 = "s,d waves";
                    PotDescr = "NLO19-500";
                    PotFlag = 11500;
                    break;
        case 11550 :PotName1 = "NLO19 (550)";
                    PotName2 = "s,d waves";
                    PotDescr = "NLO19-550";
                    PotFlag = 11550;
                    break;
        case 11600 :PotName1 = "NLO19 (600)";
                    PotName2 = "s,d waves";
                    PotDescr = "NLO19-600";
                    PotFlag = 11600;
                    break;
        case 11650 :PotName1 = "NLO19 (650)";
                    PotName2 = "s,d waves";
                    PotDescr = "NLO19-650";
                    PotFlag = 11650;
                    break;//
        case -11600:PotName1 = "LO13 (600)";
                    PotName2 = "s,d waves";
                    PotDescr = "LO13-600";
                    PotFlag = -11600;
                    break;//
        case 100:   PotName1 = "Usmani";
                    PotName2 = "original";
                    PotDescr = "Usmani";
                    PotFlag = -11600;
                    break;//
        case 101:   PotName1 = "Usmani";
                    PotName2 = "w_{c} = ";
                    PotDescr = "UsmaniFit";
                    PotFlag = -11600;
                    break;//
        default :   PotName1 = "Unknown potential";
                    PotName2 = "";
                    PotDescr = "Unknown";
                    break;
    }
    switch(SIGMA_FEED){
        case 0 :    PotName3 = "Residual p#minus #Sigma^{0}: Flat";
                    SigDescr = "SigmaFlat";
                    SigName1 = PotName3;
                    break;
        case 1 :    PotName3 = "Residual p#minus #Sigma^{0}: #chiEFT";
                    SigDescr = "SigmaChiral";
                    SigName1 = PotName3;
                    break;
        case 2 :    PotName3 = "Residual p#minus #Sigma^{0}: ESC16";
                    SigDescr = "SigmaESC16";
                    SigName1 = PotName3;
                    break;
        default :   PotName3 = "Residual p#minus #Sigma^{0}: unknown";
                    SigDescr = "Unknown";
                    SigName1 = PotName3;
                    break;
    }

    TString SourceName1;
    if(ValSourceAlpha==2) SourceName1 = "Gaussian core source";
    else SourceName1 = TString::Format("Levy core source (#alpha=%.1f)",ValSourceAlpha);

    TFile* InputFile = NULL;
    InputFile = new TFile(InputFolder+InputFileName,"read");
    TTree* plambdaTree = NULL;
    plambdaTree = (TTree*)InputFile->Get("plambdaTree");

    TGraphErrors* fitLoDummy = NULL;

	TGraphErrors* gData = NULL;
	TGraphErrors* fitData = NULL;
	TGraphErrors* fitFemto = NULL;
	TGraphErrors* fitGenuine = NULL;
	TGraphErrors* fitSignal_pL = NULL;
	TGraphErrors* fitSignal_pL_pS0 = NULL;
	TGraphErrors* fitSignal_pL_pXim = NULL;
	TGraphErrors* fitSignal_pL_pXi0 = NULL;
	TGraphErrors* fitBaseline = NULL;
	TString* DataSample = NULL;
	TString* SourceDescription = NULL;
	unsigned SEED;
	unsigned BASELINE_VAR;
	int POT_VAR;
	bool DataSyst;
	bool FitSyst;
	bool Bootstrap;
	float SourceSize;
	float SourceAlpha;
    float CuspWeight;
    float CkConv;
    float CkCutOff;
    float Norm;
    float pp1;
    float pp2;
    float pp3;
    float pp4;
    float kFemtoMin;
    float kFemtoMax;
    float kLongMin;
    float kLongMax;
    float lam_L_genuine;
    float lam_L_Sig0;
    float lam_L_Xim;
    float lam_L_Xi0;
    float lam_L_Flat;
    float lam_S0_genuine;
    float lam_S0_Flat;
    float lam_Xim_genuine;
    float lam_Xim_Flat;
    float pval;
    unsigned WhichData;
    float SBl;
    float SBpur;
    bool DefaultVariation;
    int Sigma0_Feed;
    int Xim_Feed;

	plambdaTree->SetBranchAddress("gData",&gData);
	plambdaTree->SetBranchAddress("fitData",&fitData);
	plambdaTree->SetBranchAddress("fitFemto",&fitFemto);
	plambdaTree->SetBranchAddress("fitGenuine",&fitGenuine);
	plambdaTree->SetBranchAddress("fitSignal_pL",&fitSignal_pL);
	plambdaTree->SetBranchAddress("fitSignal_pL_pS0",&fitSignal_pL_pS0);
	plambdaTree->SetBranchAddress("fitSignal_pL_pXim",&fitSignal_pL_pXim);
	plambdaTree->SetBranchAddress("fitSignal_pL_pXi0",&fitSignal_pL_pXi0);
	plambdaTree->SetBranchAddress("fitBaseline",&fitBaseline);
	plambdaTree->SetBranchAddress("DataSample",&DataSample);
	plambdaTree->SetBranchAddress("SourceDescription",&SourceDescription);
	plambdaTree->SetBranchAddress("SEED",&SEED);
	plambdaTree->SetBranchAddress("BASELINE_VAR",&BASELINE_VAR);
	plambdaTree->SetBranchAddress("POT_VAR",&POT_VAR);
	plambdaTree->SetBranchAddress("Sigma0_Feed",&Sigma0_Feed);
	plambdaTree->SetBranchAddress("Xim_Feed",&Xim_Feed);
	plambdaTree->SetBranchAddress("DataSyst",&DataSyst);
    plambdaTree->SetBranchAddress("FitSyst", &FitSyst);//
    plambdaTree->SetBranchAddress("Bootstrap", &Bootstrap);//
	plambdaTree->SetBranchAddress("SourceSize",&SourceSize);
	plambdaTree->SetBranchAddress("SourceAlpha",&SourceAlpha);
	plambdaTree->SetBranchAddress("CuspWeight",&CuspWeight);
	plambdaTree->SetBranchAddress("CkCutOff",&CkCutOff);
	plambdaTree->SetBranchAddress("CkConv",&CkConv);
	plambdaTree->SetBranchAddress("Norm",&Norm);
	plambdaTree->SetBranchAddress("pp1",&pp1);
	plambdaTree->SetBranchAddress("pp2",&pp2);
	plambdaTree->SetBranchAddress("pp3",&pp3);
	plambdaTree->SetBranchAddress("pp4",&pp4);
	plambdaTree->SetBranchAddress("kFemtoMin",&kFemtoMin);
	plambdaTree->SetBranchAddress("kFemtoMax",&kFemtoMax);
	plambdaTree->SetBranchAddress("kLongMin",&kLongMin);
	plambdaTree->SetBranchAddress("kLongMax",&kLongMax);
	plambdaTree->SetBranchAddress("lam_L_genuine",&lam_L_genuine);
	plambdaTree->SetBranchAddress("lam_L_Sig0",&lam_L_Sig0);
	plambdaTree->SetBranchAddress("lam_L_Xim",&lam_L_Xim);
	plambdaTree->SetBranchAddress("lam_L_Xi0",&lam_L_Xi0);
	plambdaTree->SetBranchAddress("lam_L_Flat",&lam_L_Flat);
	plambdaTree->SetBranchAddress("lam_S0_genuine",&lam_S0_genuine);
	plambdaTree->SetBranchAddress("lam_S0_Flat",&lam_S0_Flat);
	plambdaTree->SetBranchAddress("lam_Xim_genuine",&lam_Xim_genuine);
	plambdaTree->SetBranchAddress("lam_Xim_Flat",&lam_Xim_Flat);
	plambdaTree->SetBranchAddress("WhichData",&WhichData);
    plambdaTree->SetBranchAddress("SBl", &SBl);
    plambdaTree->SetBranchAddress("SBpur", &SBpur);
	plambdaTree->SetBranchAddress("pval",&pval);
	plambdaTree->SetBranchAddress("DefaultVariation",&DefaultVariation);

	const unsigned NumEntries = plambdaTree->GetEntries();

    const unsigned MaxNumBins = 50;
    const double xBinWidth = 12;
    TH2F* hdata_toterr = new TH2F("hdata_toterr","hdata_toterr",MaxNumBins,0,MaxNumBins,2048,0,4);
    TGraphErrors* gdata_toterr = new TGraphErrors();

    unsigned NumBestSolutions = 0;
    const unsigned Num_SourceSize=5;
    //const unsigned Num_lam_L_genuine=InputFileName.Contains("pp13TeV_HM_DimiJun20")?5:5*3;
    const unsigned Num_pur=2;
    const unsigned Num_sig_lam=5;
    const unsigned Num_xi_siglam=2;
    const unsigned Num_CuspWeight=5;
    const unsigned Num_SourceAlpha=5;
    const unsigned Num_CkConv=3;

    DLM_CommonAnaFunctions AnalysisObject;
    AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());
    TH1F* hData_pL_Stat;
    //if(DataSample[0]=="pp13TeV_HM_DimiJun20")
    //hData_pL_Stat = AnalysisObject.GetAliceExpCorrFun("pp13TeV_HM_DimiJun20","pLambda",TString::Format("L%.0f_SL4_SR6_P%.0f_0",0.529*100,0.963*100),2,false,-1);
    hData_pL_Stat = NULL;
    //else hData_pL_Stat = AnalysisObject.GetAliceExpCorrFun(DataSample[0],"pLambda","_0",2,false,-1);

    TH2F* hnsigma = new TH2F("hnsigma","hnsigma",MaxNumBins,0,MaxNumBins,2048,-32,32);
    TH2F* hnsigma_der = new TH2F("hnsigma_der","hnsigma_der",MaxNumBins,0,MaxNumBins,2048,-32,32);
    TH2F* hfdata = new TH2F("hfdata","hfdata",MaxNumBins,0,MaxNumBins,2048,0,4);
    TH2F* hbl = new TH2F("hbl","hbl",MaxNumBins,0,MaxNumBins,2048,0.5,1.5);
    TH2F* hsig = new TH2F("hsig","hsig",MaxNumBins,0,MaxNumBins,2048,0.5,1.5);
    //used for the sum of xim and xi0
    TH2F* hxim = new TH2F("hxim","hxim",MaxNumBins,0,MaxNumBins,2048,0.5,1.5);
    TH2F* hxi0 = new TH2F("hxi0","hxi0",MaxNumBins,0,MaxNumBins,2048,0.5,1.5);
    TH2F* hxi = new TH2F("hxi","hxi",MaxNumBins,0,MaxNumBins,2048,0.5,1.5);

    TH2F* hb_nsigma = new TH2F("hb_nsigma","hb_nsigma",MaxNumBins,0,MaxNumBins,2048,-32,32);
    TH2F* hb_nsigma_der = new TH2F("hb_nsigma_der","hb_nsigma_der",MaxNumBins,0,MaxNumBins,2048,-32,32);
    TH2F* hb_fdata = new TH2F("hb_fdata","hb_fdata",MaxNumBins,0,MaxNumBins,2048,0,4);
    TH2F* hb_bl = new TH2F("hb_bl","hb_bl",MaxNumBins,0,MaxNumBins,2048,0.5,1.5);
    TH2F* hb_sig = new TH2F("hb_sig","hb_sig",MaxNumBins,0,MaxNumBins,2048,0.5,1.5);
    TH2F* hb_xim = new TH2F("hb_xim","hb_xim",MaxNumBins,0,MaxNumBins,2048,0.5,1.5);
    TH2F* hb_xi0 = new TH2F("hb_xi0","hb_xi0",MaxNumBins,0,MaxNumBins,2048,0.5,1.5);
    TH2F* hb_xi = new TH2F("hb_xi","hb_xi",MaxNumBins,0,MaxNumBins,2048,0.5,1.5);
    //added ad-hoc to check Juergens question about fitting with pol3 (non-constrained)
    //this is the resulting best value of the pol1 parameter
    TH2F* hb_pol1 = new TH2F("hb_pol1","hb_pol1",MaxNumBins,0,MaxNumBins,2048,-1,1);

    TH1F* hnsigma_0_100 = new TH1F("hnsigma_0_100","hnsigma_0_100",2048,-32,32);
    TH1F* hnsigma_100_200 = new TH1F("hnsigma_100_200","hnsigma_100_200",2048,-32,32);
    TH1F* hnsigma_200_300 = new TH1F("hnsigma_200_300","hnsigma_200_300",2048,-32,32);
    TH1F* hnsigma_0_300 = new TH1F("hnsigma_0_300","hnsigma_0_300",2048,-32,32);
    TH1F* hnsigma_der_0_100 = new TH1F("hnsigma_der_0_100","hnsigma_der_0_100",2048,-32,32);
    TH1F* hnsigma_der_100_200 = new TH1F("hnsigma_der_100_200","hnsigma_der_100_200",2048,-32,32);
    TH1F* hnsigma_der_200_300 = new TH1F("hnsigma_der_200_300","hnsigma_der_200_300",2048,-32,32);
    TH1F* hnsigma_der_0_300 = new TH1F("hnsigma_der_0_300","hnsigma_der_0_300",2048,-32,32);

    TH1F* htot_nsigma_300 = new TH1F("htot_nsigma_300","htot_nsigma_300",2048,-32,32);
    TH1F* htot_nsigma_der_300 = new TH1F("htot_nsigma_der_300","htot_nsigma_der_300",2048,-32,32);
    float nsigma_0_100;
    float nsigma_100_200;
    float nsigma_200_300;
    float nsigma_0_300;
    float nsigma_der_0_100;
    float nsigma_der_100_200;
    float nsigma_der_200_300;
    float nsigma_der_0_300;
    int NDF_0_100;
    int NDF_100_200;
    int NDF_200_300;
    int NDF_0_300;
    for(unsigned uBin=0; uBin<MaxNumBins; uBin++){
        for(unsigned uNSB=0; uNSB<2048; uNSB++){
            hnsigma->SetBinContent(uBin+1,uNSB+1,0);
            hnsigma_der->SetBinContent(uBin+1,uNSB+1,0);
            hb_nsigma->SetBinContent(uBin+1,uNSB+1,0);
            hb_nsigma_der->SetBinContent(uBin+1,uNSB+1,0);
        }
    }

    DLM_Histo<float> Chi2_300;
    Chi2_300.SetUp(7);
    Chi2_300.SetUp(0,Num_SourceSize,-0.5,double(Num_SourceSize)+0.5);
    Chi2_300.SetUp(1,Num_pur,-0.5,double(Num_pur)+0.5);
    Chi2_300.SetUp(2,Num_sig_lam,-0.5,double(Num_sig_lam)+0.5);
    Chi2_300.SetUp(3,Num_xi_siglam,-0.5,double(Num_xi_siglam)+0.5);
    Chi2_300.SetUp(4,Num_CuspWeight,-0.5,double(Num_CuspWeight)+0.5);
    Chi2_300.SetUp(5,Num_SourceAlpha,-0.5,double(Num_SourceAlpha)+0.5);
    Chi2_300.SetUp(6,Num_CkConv,-0.5,double(Num_CkConv)+0.5);
    Chi2_300.Initialize();

    DLM_Histo<float> NDF_300;
    NDF_300.SetUp(7);
    NDF_300.SetUp(0,Num_SourceSize,-0.5,double(Num_SourceSize)+0.5);
    NDF_300.SetUp(1,Num_pur,-0.5,double(Num_pur)+0.5);
    NDF_300.SetUp(2,Num_sig_lam,-0.5,double(Num_sig_lam)+0.5);
    NDF_300.SetUp(3,Num_xi_siglam,-0.5,double(Num_xi_siglam)+0.5);
    NDF_300.SetUp(4,Num_CuspWeight,-0.5,double(Num_CuspWeight)+0.5);
    NDF_300.SetUp(5,Num_SourceAlpha,-0.5,double(Num_SourceAlpha)+0.5);
    NDF_300.SetUp(6,Num_CkConv,-0.5,double(Num_CkConv)+0.5);
    NDF_300.Initialize();

    DLM_Histo<float> NumChi2Entries_300;
    NumChi2Entries_300.SetUp(7);
    NumChi2Entries_300.SetUp(0,Num_SourceSize,-0.5,double(Num_SourceSize)+0.5);
    NumChi2Entries_300.SetUp(1,Num_pur,-0.5,double(Num_pur)+0.5);
    NumChi2Entries_300.SetUp(2,Num_sig_lam,-0.5,double(Num_sig_lam)+0.5);
    NumChi2Entries_300.SetUp(3,Num_xi_siglam,-0.5,double(Num_xi_siglam)+0.5);
    NumChi2Entries_300.SetUp(4,Num_CuspWeight,-0.5,double(Num_CuspWeight)+0.5);
    NumChi2Entries_300.SetUp(5,Num_SourceAlpha,-0.5,double(Num_SourceAlpha)+0.5);
    NumChi2Entries_300.SetUp(6,Num_CkConv,-0.5,double(Num_CkConv)+0.5);
    NumChi2Entries_300.Initialize();

    unsigned UsedNumBins = 0;
    float SigLamFrac;
    float XiSigLamFrac;
    TRandom3 rangenfrac1(121);
    TRandom3 rangenfrac2(121);
    for(unsigned uEntry=0; uEntry<NumEntries; uEntry++){
        //Chi2_300[uEntry] = 1e6;
        //NDF_300[uEntry] = 0;
        plambdaTree->GetEntry(uEntry);
        if(!hData_pL_Stat){
            if(*DataSample=="pp13TeV_HM_DimiJun20"||*DataSample=="pp13TeV_HM_DimiJul20")
              hData_pL_Stat = AnalysisObject.GetAliceExpCorrFun(DataSample[0],"pLambda",TString::Format("L%.0f_SL4_SR6_P%.0f_0",0.529*100,0.963*100),2,false,-1);
            else
              hData_pL_Stat = AnalysisObject.GetAliceExpCorrFun(DataSample[0],"pLambda",TString::Format("_0"),2,false,-1);
        }
        SigLamFrac=lam_L_Sig0/lam_L_genuine;
        XiSigLamFrac=(lam_L_Xi0+lam_L_Xim)/(lam_L_genuine+lam_L_Sig0);
        unsigned WhichBin[7];
        //if(*DataSample=="pp13TeV_HM_DimiJun20") GetIterCombo170620(WhichBin,SourceSize,lam_L_genuine,CuspWeight,SourceAlpha);
        if(*DataSample=="pp13TeV_HM_DimiJun20"||*DataSample=="pp13TeV_HM_DimiJul20") {}//GetIterCombo040720(WhichBin,SourceSize,SBpur,SigLamFrac,XiSigLamFrac,CuspWeight,SourceAlpha,CkConv);
        else {printf("WARNING trouble, the function MIGHT not work unless you checkout an older version (before 4th July 2020\n)"); //abort();
        }
        GetIterCombo040720(WhichBin,SourceSize,SBpur,SigLamFrac,XiSigLamFrac,CuspWeight,SourceAlpha,CkConv);
        //else GetIterCombo040620(WhichBin,SourceSize,lam_L_genuine,CuspWeight,SourceAlpha);
//printf("0\n");
        if(Sigma0_Feed!=SIGMA_FEED) continue;
//printf("1\n");
        if(BASELINE_VAR!=WhichBaseline) continue;
//printf("2\n");
        if(POT_VAR!=WhichPotential) continue;
//printf("3\n");
        if( fabs(SourceAlpha-ValSourceAlpha)>0.01 ) continue;
//printf("4\n");
        if(SourceSize<MinRad||SourceSize>MaxRad) continue;
        if(CkConv<MinCkConv||CkConv>MaxCkConv) continue;
        if(Use_SBpur&&fabs(Use_SBpur-SBpur)>0.001) continue;
//printf("5\n");
        //if(lam_L_genuine<MinLambda||lam_L_genuine>MaxLambda) continue;
//printf("6\n");
        if(CuspWeight<MinOmega||CuspWeight>MaxOmega) continue;
//printf("7\n");
        if( fabs(SigLamFrac-CuspWeight)>0.001 && Same_omega_siglam ) continue;
        if(GoodCutVar[WhichData]<MinGoodCut||GoodCutVar[WhichData]>MaxGoodCut) continue;
        if(rangenfrac1.Uniform()>FractionOfSolutions) continue;
//if(uEntry<20){
//    printf("uEntry=%u\n",uEntry);
//}
//printf("8\n");
//printf(" PROCEED!\n");
        unsigned NumBins = gData->GetN();
        if(UsedNumBins<NumBins) UsedNumBins=NumBins;
        double data_val[2];
        double data_err[2];
        double fit_val[2];
        double mom_val[2];
        double nsigma_val;
        double delta_Ck;
        double delta_k;
        double fit_der;
        double data_der;
        double data_der_err;
        double nsigma_der_val;
        nsigma_0_100=0;
        nsigma_100_200=0;
        nsigma_200_300=0;
        nsigma_0_300=0;
        nsigma_der_0_100=0;
        nsigma_der_100_200=0;
        nsigma_der_200_300=0;
        nsigma_der_0_300=0;
        NDF_0_100=0;
        NDF_100_200=0;
        NDF_200_300=0;
        NDF_0_300=0;
        for(unsigned uBin=0; uBin<NumBins; uBin++){
            gData->GetPoint(uBin,mom_val[uBin%2],data_val[uBin%2]);
            if(mom_val[uBin%2]>kLongMax) continue;
            data_err[uBin%2] = gData->GetErrorY(uBin);
            fitData->GetPoint(uBin,mom_val[uBin%2],fit_val[uBin%2]);
            nsigma_val = (fit_val[uBin%2]-data_val[uBin%2])/data_err[uBin%2];
            //first we compute the chi (and chi2) for the first der
            fit_der = (fit_val[uBin%2]-fit_val[(uBin+1)%2])/(mom_val[uBin%2]-mom_val[(uBin+1)%2]);
            data_der = (data_val[uBin%2]-data_val[(uBin+1)%2])/(mom_val[uBin%2]-mom_val[(uBin+1)%2]);
            data_der_err = sqrt(pow(data_err[uBin%2],2.)+pow(data_err[(uBin+1)%2],2.))/(mom_val[uBin%2]-mom_val[(uBin+1)%2]);
            nsigma_der_val = (fit_der-data_der)/data_der_err;

            if(mom_val[uBin%2]<=100){
                nsigma_0_100 += nsigma_val*nsigma_val;
                nsigma_der_0_100 += nsigma_der_val*nsigma_der_val;
                NDF_0_100++;
            }
            else if(mom_val[uBin%2]<=200){
                nsigma_100_200 += nsigma_val*nsigma_val;
                nsigma_der_100_200 += nsigma_der_val*nsigma_der_val;
                NDF_100_200++;
            }
            else if(mom_val[uBin%2]<=300){
                nsigma_200_300 += nsigma_val*nsigma_val;
                nsigma_der_200_300 += nsigma_der_val*nsigma_der_val;
                NDF_200_300++;
            }
            if(mom_val[uBin%2]<=300){
                nsigma_0_300 += nsigma_val*nsigma_val;
                nsigma_der_0_300 += nsigma_der_val*nsigma_der_val;
                NDF_0_300++;
            }
            else break;
        }

        //unsigned WhichBin[4];
        //GetIterCombo(WhichBin,Sigma0_Feed,SourceSize,lam_L_genuine);
        //if(*DataSample=="pp13TeV_HM_DimiJun20") GetIterCombo170620(WhichBin,SourceSize,lam_L_genuine,CuspWeight,SourceAlpha);
        //else GetIterCombo040620(WhichBin,SourceSize,lam_L_genuine,CuspWeight,SourceAlpha);
        GetIterCombo040720(WhichBin,SourceSize,SBpur,SigLamFrac,XiSigLamFrac,CuspWeight,SourceAlpha,CkConv);
        //if(Num_SourceSize==1) WhichBin[0]=0;
        //if(Num_lam_L_genuine==1) WhichBin[1]=0;
        //if(Num_CuspWeight==1) WhichBin[2]=0;
        //if(Num_SourceAlpha==1) WhichBin[3]=0;

        Chi2_300.Add(Chi2_300.GetTotBin(WhichBin),nsigma_0_300);
        NDF_300.Add(NDF_300.GetTotBin(WhichBin),NDF_0_300);
        NumChi2Entries_300.Add(NumChi2Entries_300.GetTotBin(WhichBin),1);

        nsigma_0_300 = sqrt(2)*TMath::ErfcInverse(TMath::Prob(nsigma_0_300,NDF_0_300));
        nsigma_der_0_300 = sqrt(2)*TMath::ErfcInverse(TMath::Prob(nsigma_der_0_300,NDF_0_300));
        htot_nsigma_300->Fill(nsigma_0_300);
        htot_nsigma_der_300->Fill(nsigma_der_0_300);
    }

    //turn the Chi2_300 into a mean Chi2_300
    Chi2_300.DivideHisto(NumChi2Entries_300,false);
    NDF_300.DivideHisto(NumChi2Entries_300,false);
    for(unsigned uBin=0; uBin<NDF_300.GetNbins()-1; uBin++){
        if(NDF_300.GetBinContent(uBin)!=NDF_300.GetBinContent(uBin+1)){
            ///printf("wtf %f %f\n",NDF_300.GetBinContent(uBin),NDF_300.GetBinContent(uBin+1));
        }
    }

    //find the best solution
    float MinChi2;
    unsigned* Bin_MinChi2;//must be deleted
    unsigned Num_MinChi2 = Chi2_300.FindMinima(MinChi2,Bin_MinChi2);
    int MinNdf = TMath::Nint(NDF_300.GetBinContent(Bin_MinChi2[0]));
    float MinNsigma = sqrt(2)*TMath::ErfcInverse(TMath::Prob(MinChi2,MinNdf));
    //printf("We have %u minima for the chi2=%.2f/%i = %.2f (nsigma = %.2f)\n",Num_MinChi2,MinChi2,MinNdf,MinChi2/float(MinNdf),MinNsigma);
    for(unsigned uMin=0; uMin<Num_MinChi2; uMin++){
        unsigned WhichBin[7];
        //printf(" Bin_MinChi2=%p\n",Bin_MinChi2);
        //printf(" Bin_MinChi2[0]=%u\n",Bin_MinChi2[0]);
        Chi2_300.GetBinCoordinates(Bin_MinChi2[uMin],WhichBin);
        //printf(" Solution %u: SZ=%u, SBP=%u, S:L=%u, XI:LS=%u, CW=%u, SA=%u, CC=%u\n",
        //       uMin,WhichBin[0],WhichBin[1],WhichBin[2],WhichBin[3],WhichBin[4],WhichBin[5],WhichBin[6]);
    }

    double nsigma_300_max;
    double nsigma_300_min;
    //I want to keep the best 68%. => I should get rid of the upper most 32%,
    //a dirty trick here is to just select the 36% central interval, leaving out 32% on both sides of the distribution.
    double nsigma_300_med = GetCentralInterval(*htot_nsigma_300,0.32,nsigma_300_min,nsigma_300_max,true);
    double nsigma_der_300_max;
    double nsigma_der_300_min;
    double nsigma_der_300_med = GetCentralInterval(*htot_nsigma_der_300,0.32,nsigma_der_300_min,nsigma_der_300_max,true);

    float Best_SourceSize;
    float Best_SourceAlpha;
    float Best_CuspWeight;
    float Best_lam_L_genuine;
    float Best_SigLamFrac;
    float Best_XiSigLamFrac;
    float Best_Purity;
    float Best_CkConv;

    for(unsigned uEntry=0; uEntry<NumEntries; uEntry++){
        plambdaTree->GetEntry(uEntry);
        SigLamFrac=lam_L_Sig0/lam_L_genuine;
        XiSigLamFrac=(lam_L_Xi0+lam_L_Xim)/(lam_L_genuine+lam_L_Sig0);
        unsigned WhichBin[7];
        //if(*DataSample=="pp13TeV_HM_DimiJun20") GetIterCombo170620(WhichBin,SourceSize,lam_L_genuine,CuspWeight,SourceAlpha);
        //else GetIterCombo040620(WhichBin,SourceSize,lam_L_genuine,CuspWeight,SourceAlpha);
        GetIterCombo040720(WhichBin,SourceSize,SBpur,SigLamFrac,XiSigLamFrac,CuspWeight,SourceAlpha,CkConv);

        unsigned NumBins = gData->GetN();
        for(unsigned uBin=0; uBin<NumBins; uBin++){
            if(GoodCutVar[WhichData]<MinGoodCut||GoodCutVar[WhichData]>MaxGoodCut) break;
            if(uBin>=MaxNumBins) break;
            double xVal,yVal,yErr;
            gData->GetPoint(uBin,xVal,yVal);
            yErr = gData->GetErrorY(uBin);
            //the data is only plotted for the default iter (no bootstrap)
            if(DefaultVariation){
                hData_pL_Stat->SetBinContent(uBin+1,yVal);
                hData_pL_Stat->SetBinError(uBin+1,yErr);
            }
            //the error is computed based on all iterations
            //we include the purity variations in the systematics, to avoid a
            //model dependent representation of the data!
            gdata_toterr->SetPoint(uBin,xVal,yVal);
            hdata_toterr->Fill(double(uBin)+0.5,yVal);
        }

        if(Sigma0_Feed!=SIGMA_FEED) continue;
        if(BASELINE_VAR!=WhichBaseline) continue;
        if(POT_VAR!=WhichPotential) continue;
        if( fabs(SourceAlpha-ValSourceAlpha)>0.01 ) continue;
        if(SourceSize<MinRad||SourceSize>MaxRad) continue;
        if(CkConv<MinCkConv||CkConv>MaxCkConv) continue;
        if(Use_SBpur&&fabs(Use_SBpur-SBpur)>0.001) continue;
//printf("%f %f %f\n",MinLambda,lam_L_genuine,MaxLambda);
        //if(lam_L_genuine<MinLambda||lam_L_genuine>MaxLambda) continue;
        if(CuspWeight<MinOmega||CuspWeight>MaxOmega) continue;
        if( fabs(SigLamFrac-CuspWeight)>0.001 && Same_omega_siglam ) continue;
        if(GoodCutVar[WhichData]<MinGoodCut||GoodCutVar[WhichData]>MaxGoodCut) continue;
        if(rangenfrac2.Uniform()>FractionOfSolutions) continue;
//if(uEntry<20){
//    printf("uEntry=%u\n",uEntry);
//}
//printf("GO ON\n");

        //the [2] was used for saving the previous bin for the derivative
        double data_val[2];
        double data_err[2];
        double fit_val[2];
        double mom_val[2];
        double bl_val;
        double sig_val;
        double xim_val;
        double xi0_val;
        double xi_val;
        double nsigma_val;
        double delta_Ck;
        double delta_k;
        double fit_der;
        double data_der;
        double data_der_err;
        double nsigma_der_val;
        nsigma_0_100=0;
        nsigma_100_200=0;
        nsigma_200_300=0;
        nsigma_0_300=0;
        nsigma_der_0_100=0;
        nsigma_der_100_200=0;
        nsigma_der_200_300=0;
        nsigma_der_0_300=0;
        NDF_0_100=0;
        NDF_100_200=0;
        NDF_200_300=0;
        NDF_0_300=0;
        for(unsigned uBin=0; uBin<NumBins; uBin++){
            gData->GetPoint(uBin,mom_val[uBin%2],data_val[uBin%2]);
            if(mom_val[uBin%2]>kLongMax) continue;
            data_err[uBin%2] = gData->GetErrorY(uBin);
            fitData->GetPoint(uBin,mom_val[uBin%2],fit_val[uBin%2]);
            nsigma_val = (fit_val[uBin%2]-data_val[uBin%2])/data_err[uBin%2];
            //first we compute the chi (and chi2) for the first der
            fit_der = (fit_val[uBin%2]-fit_val[(uBin+1)%2])/(mom_val[uBin%2]-mom_val[(uBin+1)%2]);
            data_der = (data_val[uBin%2]-data_val[(uBin+1)%2])/(mom_val[uBin%2]-mom_val[(uBin+1)%2]);
            data_der_err = sqrt(pow(data_err[uBin%2],2.)+pow(data_err[(uBin+1)%2],2.))/(mom_val[uBin%2]-mom_val[(uBin+1)%2]);
            nsigma_der_val = (fit_der-data_der)/data_der_err;
            if(mom_val[uBin%2]<=100){
                nsigma_0_100 += nsigma_val*nsigma_val;
                nsigma_der_0_100 += nsigma_der_val*nsigma_der_val;
                NDF_0_100++;
            }
            else if(mom_val[uBin%2]<=200){
                nsigma_100_200 += nsigma_val*nsigma_val;
                nsigma_der_100_200 += nsigma_der_val*nsigma_der_val;
                NDF_100_200++;
            }
            else if(mom_val[uBin%2]<=300){
                nsigma_200_300 += nsigma_val*nsigma_val;
                nsigma_der_200_300 += nsigma_der_val*nsigma_der_val;
                NDF_200_300++;
            }
            if(mom_val[uBin%2]<=300){
                nsigma_0_300 += nsigma_val*nsigma_val;
                nsigma_der_0_300 += nsigma_der_val*nsigma_der_val;
                NDF_0_300++;
            }
            else break;
        }
        nsigma_0_300 = sqrt(2)*TMath::ErfcInverse(TMath::Prob(nsigma_0_300,NDF_0_300));
        htot_nsigma_300->Fill(nsigma_0_300);
        nsigma_0_100 = sqrt(2)*TMath::ErfcInverse(TMath::Prob(nsigma_0_100,NDF_0_100));
        nsigma_100_200 = sqrt(2)*TMath::ErfcInverse(TMath::Prob(nsigma_100_200,NDF_100_200));
        nsigma_200_300 = sqrt(2)*TMath::ErfcInverse(TMath::Prob(nsigma_200_300,NDF_200_300));

        nsigma_der_0_300 = sqrt(2)*TMath::ErfcInverse(TMath::Prob(nsigma_der_0_300,NDF_0_300));
        htot_nsigma_der_300->Fill(nsigma_der_0_300);
        nsigma_der_0_100 = sqrt(2)*TMath::ErfcInverse(TMath::Prob(nsigma_der_0_100,NDF_0_100));
        nsigma_der_100_200 = sqrt(2)*TMath::ErfcInverse(TMath::Prob(nsigma_der_100_200,NDF_100_200));
        nsigma_der_200_300 = sqrt(2)*TMath::ErfcInverse(TMath::Prob(nsigma_der_200_300,NDF_200_300));

        for(unsigned uBin=0; uBin<NumBins; uBin++){
            gData->GetPoint(uBin,mom_val[uBin%2],data_val[uBin%2]);
            if(mom_val[uBin%2]>kLongMax) continue;
            data_err[uBin%2] = gData->GetErrorY(uBin);
            fitData->GetPoint(uBin,mom_val[uBin%2],fit_val[uBin%2]);
            fitBaseline->GetPoint(uBin,mom_val[uBin%2],bl_val);
//printf("uBin=%u; mom_val=%f; bl_val=%f\n",uBin,mom_val[uBin%2],bl_val);
            fitSignal_pL_pS0->GetPoint(uBin,mom_val[uBin%2],sig_val);
            sig_val *= bl_val;
            fitSignal_pL_pXim->GetPoint(uBin,mom_val[uBin%2],xim_val);
            if(fitSignal_pL_pXi0) fitSignal_pL_pXi0->GetPoint(uBin,mom_val[uBin%2],xi0_val);
            else xi0_val=1;
            xi_val = (xim_val-1)+(xi0_val-1)+1;
            xim_val *= bl_val;
            xi0_val *= bl_val;
            xi_val *= bl_val;

            nsigma_val = (fit_val[uBin%2]-data_val[uBin%2])/data_err[uBin%2];

            //hnsigma->SetBinContent(uBin+1,hnsigma->GetYaxis()->FindBin(nsigma_val),hnsigma->GetBinContent(uBin+1)+1);
            //hnsigma_der->SetBinContent(uBin+1,hnsigma_der->GetYaxis()->FindBin(nsigma_val),hnsigma_der->GetBinContent(uBin+1)+1);
            //hfdata->SetBinContent(uBin+1,hfdata->GetYaxis()->FindBin(fit_val[uBin%2]),hfdata->GetBinContent(uBin+1)+1);
            //hbl->SetBinContent(uBin+1,hbl->GetYaxis()->FindBin(bl_val),hbl->GetBinContent(uBin+1)+1);
            //hsig->SetBinContent(uBin+1,hsig->GetYaxis()->FindBin(sig_val),hsig->GetBinContent(uBin+1)+1);
            ////hxim->SetBinContent(uBin+1,hxim->GetYaxis()->FindBin( (xim_val-1) + (xi0_val-1) + 1 ),hxim->GetBinContent(uBin+1)+1);
            //hxim->SetBinContent(uBin+1,hxim->GetYaxis()->FindBin( xim_val ),hxim->GetBinContent(uBin+1)+1);
            //hxi0->SetBinContent(uBin+1,hxi0->GetYaxis()->FindBin( xi0_val ),hxi0->GetBinContent(uBin+1)+1);
            //hxi->SetBinContent(uBin+1,hxi->GetYaxis()->FindBin( xi_val ),hxi->GetBinContent(uBin+1)+1);

            hnsigma->Fill(double(uBin)+0.5,nsigma_val);
            hnsigma_der->Fill(double(uBin)+0.5,nsigma_val);
            hfdata->Fill(double(uBin)+0.5,fit_val[uBin%2]);
            hbl->Fill(double(uBin)+0.5,bl_val);
            hsig->Fill(double(uBin)+0.5,sig_val);
            hxim->Fill(double(uBin)+0.5,xim_val);
            hxi0->Fill(double(uBin)+0.5,xi0_val);
            hxi->Fill(double(uBin)+0.5,xi_val);

//if(!uEntry)
//printf("%.1f: xim_val=%.4f xi0_val=%.4f xi_val=%.4f\n",mom_val[uBin%2],xim_val,xi0_val,xi_val);
            unsigned WhichBinCurrent[7];
            //if(*DataSample=="pp13TeV_HM_DimiJun20") GetIterCombo170620(WhichBinCurrent,SourceSize,lam_L_genuine,CuspWeight,SourceAlpha);
            //else GetIterCombo040620(WhichBinCurrent,SourceSize,lam_L_genuine,CuspWeight,SourceAlpha);
            GetIterCombo040720(WhichBinCurrent,SourceSize,SBpur,SigLamFrac,XiSigLamFrac,CuspWeight,SourceAlpha,CkConv);
            //if(Num_SourceSize==1) WhichBinCurrent[0]=0;
            //if(Num_lam_L_genuine==1) WhichBinCurrent[1]=0;
            //if(Num_CuspWeight==1) WhichBinCurrent[2]=0;
            //if(Num_SourceAlpha==1) WhichBinCurrent[3]=0;
            bool OneOfTheBest=false;
            for(unsigned uMin=0; uMin<Num_MinChi2; uMin++){
                unsigned WhichBinBest[7];
                Chi2_300.GetBinCoordinates(Bin_MinChi2[uMin],WhichBinBest);
                if(WhichBinBest[0]==WhichBinCurrent[0]&&
                    WhichBinBest[1]==WhichBinCurrent[1]&&
                    WhichBinBest[2]==WhichBinCurrent[2] &&
                    WhichBinBest[3]==WhichBinCurrent[3] &&
                    WhichBinBest[4]==WhichBinCurrent[4] &&
                    WhichBinBest[5]==WhichBinCurrent[5] &&
                    WhichBinBest[6]==WhichBinCurrent[6]
                    ) OneOfTheBest=true;
                //if(!uBin) printf(" Solution %u: SZ=%u, SBP=%u, S:L=%u, XI:LS=%u, CW=%u, SA=%u, CC=%u\n",uMin,
                //WhichBinBest[0],WhichBinBest[1],WhichBinBest[2],WhichBinBest[3],WhichBinBest[4],WhichBinBest[5],WhichBinBest[6]);

                //best purity iteration -> use to plot the data
                //if(WhichBinBest[1]==WhichBinCurrent[1]){
                    //for(unsigned uEntry=0; uEntry<NumEntries; uEntry++){
/// initially the syst err on the data was computed here
                    //}
                //}
            }

            if(OneOfTheBest){
//printf("OneOfTheBest\n");
//printf("WhichData=%u\n",WhichData);

                if(!uBin) NumBestSolutions++;
                //hb_nsigma->SetBinContent(uBin+1,hb_nsigma->GetYaxis()->FindBin(nsigma_val),hb_nsigma->GetBinContent(uBin+1)+1);
                //hb_nsigma_der->SetBinContent(uBin+1,hb_nsigma_der->GetYaxis()->FindBin(nsigma_val),hb_nsigma_der->GetBinContent(uBin+1)+1);
                //hb_fdata->SetBinContent(uBin+1,hb_fdata->GetYaxis()->FindBin(fit_val[uBin%2]),hb_fdata->GetBinContent(uBin+1)+1);
                //hb_bl->SetBinContent(uBin+1,hb_bl->GetYaxis()->FindBin(bl_val),hb_bl->GetBinContent(uBin+1)+1);
                //hb_sig->SetBinContent(uBin+1,hb_sig->GetYaxis()->FindBin(sig_val),hb_sig->GetBinContent(uBin+1)+1);
                //hb_xim->SetBinContent(uBin+1,hb_xim->GetYaxis()->FindBin(xim_val),hb_xim->GetBinContent(uBin+1)+1);
                //hb_xi0->SetBinContent(uBin+1,hb_xi0->GetYaxis()->FindBin(xi0_val),hb_xi0->GetBinContent(uBin+1)+1);
                //hb_xi->SetBinContent(uBin+1,hb_xi->GetYaxis()->FindBin(xi_val),hb_xi->GetBinContent(uBin+1)+1);

                hb_nsigma->Fill(double(uBin)+0.5,nsigma_val);
                hb_nsigma_der->Fill(double(uBin)+0.5,nsigma_val);
                hb_fdata->Fill(double(uBin)+0.5,fit_val[uBin%2]);
                hb_bl->Fill(double(uBin)+0.5,bl_val);
                hb_sig->Fill(double(uBin)+0.5,sig_val);
                hb_xim->Fill(double(uBin)+0.5,xim_val);
                hb_xi0->Fill(double(uBin)+0.5,xi0_val);
                hb_xi->Fill(double(uBin)+0.5,xi_val);
                hb_pol1->Fill(double(uBin)+0.5,1000.*3.*pp1*pp2*pp3);


/*
if(!uBin){
    printf("pp2 = %f\n",pp2);
    printf("r = %f\n",SourceSize);
    printf("p = %f\n",SBpur);
    printf("S:L = %f\n",SigLamFrac);
    printf("XI:LS = %f\n",XiSigLamFrac);
    printf("CW = %f\n",CuspWeight);
    printf("SA = %f\n",SourceAlpha);
    printf("CC = %f\n",CkConv);
}
printf("k=%.0f, bl=%.5f\n",mom_val[uBin%2],bl_val);
*/
                Best_SourceSize = SourceSize;
                Best_SourceAlpha = SourceAlpha;
                Best_CuspWeight = CuspWeight;
                Best_lam_L_genuine = lam_L_genuine;
                Best_SigLamFrac = SigLamFrac;
                Best_XiSigLamFrac = XiSigLamFrac;
                Best_Purity = SBpur;
                Best_CkConv = CkConv;

//printf("We are the best! (%i, %p)\n",COMPARE_TO_LO,fitLoDummy);
                if(COMPARE_TO_LO&&!fitLoDummy){
                    //DUMMY LO. Compute the genuine LO, smeared for momentum
                    CATS AB_pL;
                    DLM_Ck* Ck_pL;
                    AB_pL.SetMomBins(32,0,384);
                    DLM_CommonAnaFunctions AnalysisObject;
                    AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());
                    TH2F* hResolution_pL = AnalysisObject.GetResolutionMatrix(*DataSample,"pLambda");//pp13TeV_HM_Dec19
                    AnalysisObject.SetUpCats_pL(AB_pL,"Chiral_Coupled_SPD","McGauss_ResoTM",-11600,202);//NLO_Coupled_S

                    AB_pL.SetChannelWeight(7,1./4.*Best_CuspWeight);//1S0 SN(s) -> LN(s)
                    AB_pL.SetChannelWeight(8,3./4.*Best_CuspWeight);//3S1 SN(s) -> LN(s)
                    AB_pL.SetChannelWeight(10,3./4.*Best_CuspWeight);//3S1 SN(d) -> LN(s)
                    AB_pL.SetChannelWeight(13,3./20.*Best_CuspWeight);//3D1 SN(d) -> LN(d)
                    AB_pL.SetChannelWeight(15,3./20.*Best_CuspWeight);//3D1 SN(s) -> LN(d)

                    AB_pL.SetAnaSource(0,Best_SourceSize);
                    if(SourceDescription->Contains("Mc")){
                        AB_pL.SetAnaSource(0,Best_SourceSize);
                        AB_pL.SetAnaSource(1,Best_SourceAlpha);
                    }
                    AB_pL.SetNotifications(CATS::nError);
                    AB_pL.KillTheCat();

                    Ck_pL = new DLM_Ck(AB_pL.GetNumSourcePars(),0,AB_pL,45,0,540);
                    Ck_pL->SetSourcePar(0,AB_pL.GetAnaSourcePar(0));
                    if(SourceDescription->Contains("Mc")){
                        Ck_pL->SetSourcePar(0,AB_pL.GetAnaSourcePar(0));
                        Ck_pL->SetSourcePar(1,AB_pL.GetAnaSourcePar(1));
                    }
                    Ck_pL->SetCutOff(CkCutOff,CkConv);
                    Ck_pL->Update();
                    DLM_CkDecomposition* CkDec_pL = new DLM_CkDecomposition("pLambda",0,*Ck_pL,hResolution_pL);
                    CkDec_pL->Update();

                    fitLoDummy = new TGraphErrors();
                    fitLoDummy->SetName("fitLoDummy");
                    double KSTAR_VAL;
                    double LO_VAL;
                    double SIG_VAL;
                    double XIM_VAL;
                    double XI0_VAL;
                    double BL_VAL;
                    for(unsigned uBin=0; uBin<40; uBin++){
                        KSTAR_VAL = 6.+12.*double(uBin);
                        LO_VAL = (CkDec_pL->EvalCk(KSTAR_VAL)-1.)*lam_L_genuine;
                        SIG_VAL = fitSignal_pL_pS0->Eval(KSTAR_VAL)-1.;
                        XIM_VAL = fitSignal_pL_pXim->Eval(KSTAR_VAL)-1.;
                        if(fitSignal_pL_pXi0) XI0_VAL = fitSignal_pL_pXi0->Eval(KSTAR_VAL)-1.;
                        else XI0_VAL=0;
                        BL_VAL = fitBaseline->Eval(KSTAR_VAL);
                        fitLoDummy->SetPoint(uBin,KSTAR_VAL,((LO_VAL+SIG_VAL+XIM_VAL+XI0_VAL)+1.)*BL_VAL);
                        //fitLoDummy->SetPointError(uBin,0,);
                        //printf("C(%.2f) = %.3f\n",KSTAR_VAL,((LO_VAL+SIG_VAL+XI_VAL)+1.)*BL_VAL);
                        //printf(" LO = %.3f\n",LO_VAL);
                        //printf(" SIG_VAL = %.3f\n",SIG_VAL);
                        //printf(" XI_VAL = %.3f\n",XI_VAL);
                        //printf(" lam_L_Flat = %.3f\n",lam_L_Flat);
                        //printf(" BL_VAL = %.3f\n",BL_VAL);
                        //printf(" CATS = %.3f\n",AB_pL.EvalCorrFun(KSTAR_VAL));
                        //usleep(250e3);
                    }
                    delete Ck_pL;
                    delete CkDec_pL;
                }
            }
        }
        hnsigma_0_100->Fill(nsigma_0_100);
        hnsigma_100_200->Fill(nsigma_100_200);
        hnsigma_200_300->Fill(nsigma_200_300);
        hnsigma_0_300->Fill(nsigma_0_300);

        hnsigma_der_0_100->Fill(nsigma_der_0_100);
        hnsigma_der_100_200->Fill(nsigma_der_100_200);
        hnsigma_der_200_300->Fill(nsigma_der_200_300);
        hnsigma_der_0_300->Fill(nsigma_der_0_300);

    }

    printf("NumBestSolutions = %u\n",NumBestSolutions);



    TGraphErrors* ge_Bl = new TGraphErrors();
    ge_Bl->SetName("ge_Bl");
    TGraphErrors* ge_Sig = new TGraphErrors();
    ge_Sig->SetName("ge_Sig");
    TGraphErrors* ge_Xim = new TGraphErrors();
    ge_Xim->SetName("ge_Xim");
    TGraphErrors* ge_Xi0 = new TGraphErrors();
    ge_Xi0->SetName("ge_Xi0");
    TGraphErrors* ge_Xi = new TGraphErrors();
    ge_Xi->SetName("ge_Xi");
    TGraphErrors* ge_Fit = new TGraphErrors();
    ge_Fit->SetName("ge_Fit");
    TGraphErrors* ge_SigmaSlice = new TGraphErrors();
    ge_SigmaSlice->SetName("ge_SigmaSlice");
    TGraphErrors* ge_SigmaSlice_der = new TGraphErrors();
    ge_SigmaSlice_der->SetName("ge_SigmaSlice_der");

    TGraphErrors* geb_Bl = new TGraphErrors();
    geb_Bl->SetName("geb_Bl");
    TGraphErrors* geb_Sig = new TGraphErrors();
    geb_Sig->SetName("geb_Sig");
    TGraphErrors* geb_Xim = new TGraphErrors();
    geb_Xim->SetName("geb_Xim");
    TGraphErrors* geb_Xi0 = new TGraphErrors();
    geb_Xi0->SetName("geb_Xi0");
    TGraphErrors* geb_Xi = new TGraphErrors();
    geb_Xi->SetName("geb_Xi");
    TGraphErrors* geb_Fit = new TGraphErrors();
    geb_Fit->SetName("geb_Fit");
    TGraphErrors* geb_SigmaSlice = new TGraphErrors();
    geb_SigmaSlice->SetName("geb_SigmaSlice");
    TGraphErrors* geb_SigmaSlice_der = new TGraphErrors();
    geb_SigmaSlice_der->SetName("geb_SigmaSlice_der");

    TF1* fe_Sigma = new TF1("fe_Sigma","[0]",0,300);

    fe_Sigma->SetName("fe_Sigma");
    for(unsigned uBin=0; uBin<UsedNumBins; uBin++){
        double MOM = hData_pL_Stat->GetBinCenter(uBin+1);
        double MinVal,MaxVal,MedVal;

        TH1* hProj;

        hProj = hfdata->ProjectionY("hProj",uBin+1,uBin+1);
        if(!hProj->Integral()) continue;
        MedVal = GetCentralInterval(*hProj,0.68,MinVal,MaxVal,true);
        ge_Fit->SetPoint(uBin,MOM,(MaxVal+MinVal)*0.5);
        ge_Fit->SetPointError(uBin,0,(MaxVal-MinVal)*0.5);
        delete hProj;

        hProj = hbl->ProjectionY("hProj",uBin+1,uBin+1);
        MedVal = GetCentralInterval(*hProj,0.68,MinVal,MaxVal,true);
        ge_Bl->SetPoint(uBin,MOM,(MaxVal+MinVal)*0.5);
        ge_Bl->SetPointError(uBin,0,(MaxVal-MinVal)*0.5);
        delete hProj;

        hProj = hsig->ProjectionY("hProj",uBin+1,uBin+1);
        MedVal = GetCentralInterval(*hProj,0.68,MinVal,MaxVal,true);
        ge_Sig->SetPoint(uBin,MOM,(MaxVal+MinVal)*0.5);
        ge_Sig->SetPointError(uBin,0,(MaxVal-MinVal)*0.5);
        delete hProj;

        hProj = hxim->ProjectionY("hProj",uBin+1,uBin+1);
        MedVal = GetCentralInterval(*hProj,0.68,MinVal,MaxVal,true);
        ge_Xim->SetPoint(uBin,MOM,(MaxVal+MinVal)*0.5);
        ge_Xim->SetPointError(uBin,0,(MaxVal-MinVal)*0.5);
        delete hProj;

        hProj = hxi0->ProjectionY("hProj",uBin+1,uBin+1);
        MedVal = GetCentralInterval(*hProj,0.68,MinVal,MaxVal,true);
        ge_Xi0->SetPoint(uBin,MOM,(MaxVal+MinVal)*0.5);
        ge_Xi0->SetPointError(uBin,0,(MaxVal-MinVal)*0.5);
        delete hProj;

        hProj = hxi->ProjectionY("hProj",uBin+1,uBin+1);
        MedVal = GetCentralInterval(*hProj,0.68,MinVal,MaxVal,true);
        ge_Xi->SetPoint(uBin,MOM,(MaxVal+MinVal)*0.5);
        ge_Xi->SetPointError(uBin,0,(MaxVal-MinVal)*0.5);
        delete hProj;

        hProj = hnsigma->ProjectionY("hProj",uBin+1,uBin+1);
        MedVal = GetCentralInterval(*hProj,0.68,MinVal,MaxVal,true);
        ge_SigmaSlice->SetPoint(uBin,MOM,(MaxVal+MinVal)*0.5);
        ge_SigmaSlice->SetPointError(uBin,0,(MaxVal-MinVal)*0.5);
        delete hProj;

        hProj = hnsigma_der->ProjectionY("hProj",uBin+1,uBin+1);
        MedVal = GetCentralInterval(*hProj,0.68,MinVal,MaxVal,true);
        ge_SigmaSlice_der->SetPoint(uBin,MOM,(MaxVal+MinVal)*0.5);
        ge_SigmaSlice_der->SetPointError(uBin,0,(MaxVal-MinVal)*0.5);
        delete hProj;

        hProj = hb_fdata->ProjectionY("hProj",uBin+1,uBin+1);
        if(!hProj->Integral()) continue;
        MedVal = GetCentralInterval(*hProj,0.68,MinVal,MaxVal,true);
        geb_Fit->SetPoint(uBin,MOM,(MaxVal+MinVal)*0.5);
        geb_Fit->SetPointError(uBin,0,(MaxVal-MinVal)*0.5);
        delete hProj;

        hProj = hb_bl->ProjectionY("hProj",uBin+1,uBin+1);
        MedVal = GetCentralInterval(*hProj,0.68,MinVal,MaxVal,true);
//printf("uBin=%u; MOM=%f; MedVal=%f; MinVal=%f; MaxVal=%f\n",
//       uBin,MOM,MedVal,MinVal,MaxVal);
        geb_Bl->SetPoint(uBin,MOM,(MaxVal+MinVal)*0.5);
        geb_Bl->SetPointError(uBin,0,(MaxVal-MinVal)*0.5);
        delete hProj;

        hProj = hb_sig->ProjectionY("hProj",uBin+1,uBin+1);
        MedVal = GetCentralInterval(*hProj,0.68,MinVal,MaxVal,true);
        geb_Sig->SetPoint(uBin,MOM,(MaxVal+MinVal)*0.5);
        geb_Sig->SetPointError(uBin,0,(MaxVal-MinVal)*0.5);
        delete hProj;

        hProj = hb_xim->ProjectionY("hProj",uBin+1,uBin+1);
        MedVal = GetCentralInterval(*hProj,0.68,MinVal,MaxVal,true);
        geb_Xim->SetPoint(uBin,MOM,(MaxVal+MinVal)*0.5);
        geb_Xim->SetPointError(uBin,0,(MaxVal-MinVal)*0.5);
        delete hProj;

        hProj = hb_xi0->ProjectionY("hProj",uBin+1,uBin+1);
        MedVal = GetCentralInterval(*hProj,0.68,MinVal,MaxVal,true);
        geb_Xi0->SetPoint(uBin,MOM,(MaxVal+MinVal)*0.5);
        geb_Xi0->SetPointError(uBin,0,(MaxVal-MinVal)*0.5);
        delete hProj;

        hProj = hb_xi->ProjectionY("hProj",uBin+1,uBin+1);
        MedVal = GetCentralInterval(*hProj,0.68,MinVal,MaxVal,true);
        geb_Xi->SetPoint(uBin,MOM,(MaxVal+MinVal)*0.5);
        geb_Xi->SetPointError(uBin,0,(MaxVal-MinVal)*0.5);
        delete hProj;

        hProj = hb_nsigma->ProjectionY("hProj",uBin+1,uBin+1);
        MedVal = GetCentralInterval(*hProj,0.68,MinVal,MaxVal,true);
        geb_SigmaSlice->SetPoint(uBin,MOM,(MaxVal+MinVal)*0.5);
        geb_SigmaSlice->SetPointError(uBin,0,(MaxVal-MinVal)*0.5);
        delete hProj;

        hProj = hb_nsigma_der->ProjectionY("hProj",uBin+1,uBin+1);
        MedVal = GetCentralInterval(*hProj,0.68,MinVal,MaxVal,true);
        geb_SigmaSlice_der->SetPoint(uBin,MOM,(MaxVal+MinVal)*0.5);
        geb_SigmaSlice_der->SetPointError(uBin,0,(MaxVal-MinVal)*0.5);
        delete hProj;

        fe_Sigma->FixParameter(0,MinNsigma);
    }
/*
    for(unsigned uEntry=0; uEntry<NumEntries; uEntry++){
        plambdaTree->GetBranch("gData")->GetEntry(uEntry);
        double xVal,yVal;
        for(unsigned uBin=0; uBin<gData->GetN(); uBin++){
            if(uBin>=MaxNumBins) break;
            gData->GetPoint(uBin,xVal,yVal);
            gdata_toterr->SetPoint(uBin,xVal,yVal);
            hdata_toterr->Fill(double(uBin)+0.5,yVal);
        }
    }
*/
    for(unsigned uBin=0; uBin<gdata_toterr->GetN(); uBin++){
        double xVal,yVal;
        double yMean,yStddev;
        TH1* hProj = hdata_toterr->ProjectionY("hProj",uBin+1,uBin+1);
        yMean = hProj->GetMean();
        yStddev = hProj->GetStdDev();
        gdata_toterr->GetPoint(uBin,xVal,yVal);
        gdata_toterr->SetPoint(uBin,xVal,yMean);
        gdata_toterr->SetPointError(uBin,xBinWidth*0.5,yStddev);
        delete hProj;
    }



    gStyle->SetCanvasPreferGL(1);
    if(PlotsType==1) SetStyle_pLambda();
    else SetStyle();


    //TColor SigmaColor(6000,132,112,255);
    int ColorNLO13 = kRed;//+1
    int ColorNLO19 = kCyan;
    int ColorLO13 = kGreen;//+1
    int ColorSigmaNLO = kAzure;//kBlue-7 or //kTeal
    int ColorXiLattice = kPink+1;//kOrange+2 or +1
    int ColorBaseline = kGray+1;//kGray+1

    int ColorInteraction =  WhichPotential<0?ColorLO13:
                            WhichPotential<10000?ColorNLO13:ColorNLO19;
    //if(WhichPotential==-11600) ColorInteraction = kGreen;

    if(fitLoDummy){
        fitLoDummy->SetFillColorAlpha(ColorLO13,0.70);
        fitLoDummy->SetLineColor(ColorLO13);
        fitLoDummy->SetLineWidth(4);
    }

    ge_Fit->SetFillColorAlpha(ColorInteraction,0.70);
    ge_Fit->SetLineColor(ColorInteraction);
    ge_Fit->SetLineWidth(5);

    ge_Bl->SetFillColorAlpha(ColorBaseline,0.70);//kCyan-8
    ge_Bl->SetLineColor(ColorBaseline);//kCyan-8
    ge_Bl->SetLineWidth(5);

    ge_Sig->SetFillColorAlpha(ColorSigmaNLO,0.70);//kTeal
    ge_Sig->SetLineColor(ColorSigmaNLO);//kTeal
    ge_Sig->SetLineWidth(5);

    ge_Xim->SetFillColorAlpha(ColorXiLattice,0.70);
    ge_Xim->SetLineColor(ColorXiLattice);
    ge_Xim->SetLineWidth(5);
    ge_Xim->SetLineStyle(2);


    ge_Xi0->SetFillColorAlpha(ColorXiLattice,0.70);
    ge_Xi0->SetLineColor(ColorXiLattice);
    ge_Xi0->SetLineWidth(5);
    ge_Xi0->SetLineStyle(5);

    ge_Xi->SetFillColorAlpha(ColorXiLattice,0.70);
    ge_Xi->SetLineColor(ColorXiLattice);
    ge_Xi->SetLineWidth(5);

    ge_SigmaSlice->SetFillColorAlpha(ColorInteraction,0.70);
    ge_SigmaSlice_der->SetFillColorAlpha(kBlue+1,0.70);

    geb_Fit->SetFillColorAlpha(ColorInteraction,0.70);
    geb_Fit->SetLineColor(ColorInteraction);
    geb_Fit->SetLineWidth(3);

    geb_Bl->SetFillColorAlpha(ColorBaseline,0.70);//kCyan-8
    geb_Bl->SetLineColor(ColorBaseline);//kCyan-8
    geb_Bl->SetLineWidth(5);

    geb_Sig->SetFillColorAlpha(ColorSigmaNLO,0.70);//kTeal
    geb_Sig->SetLineColor(ColorSigmaNLO);//kTeal
    geb_Sig->SetLineWidth(5);

    geb_Xim->SetFillColorAlpha(ColorXiLattice,0.70);
    geb_Xim->SetLineColor(ColorXiLattice);
    geb_Xim->SetLineWidth(5);
    geb_Xim->SetLineStyle(2);

    geb_Xi0->SetFillColorAlpha(ColorXiLattice,0.70);
    geb_Xi0->SetLineColor(ColorXiLattice);
    geb_Xi0->SetLineWidth(5);
    geb_Xi0->SetLineStyle(5);

    geb_Xi->SetFillColorAlpha(ColorXiLattice,0.70);
    geb_Xi->SetLineColor(ColorXiLattice);
    geb_Xi->SetLineWidth(5);

    geb_SigmaSlice->SetFillColorAlpha(ColorInteraction,0.70);
    geb_SigmaSlice->SetLineColor(ColorInteraction);
    geb_SigmaSlice->SetLineWidth(3);
    geb_SigmaSlice_der->SetFillColorAlpha(kBlue+1,0.70);


    fe_Sigma->SetLineColor(kBlack);
    fe_Sigma->SetLineWidth(3);
    fe_Sigma->SetLineStyle(2);

    DLM_SubPads* DlmPad = new DLM_SubPads(720,1080);
    DlmPad->AddSubPad(0,1,0.5,1);
    DlmPad->AddSubPad(0,1,0.25,0.5);
    DlmPad->AddSubPad(0,1,0.,0.25);
    if(PlotsType==1){
      DlmPad->SetMargin(0,0.18,0.02,0.0,0.02);
      DlmPad->SetMargin(1,0.18,0.02,0.0,0.0);
      DlmPad->SetMargin(2,0.18,0.02,0.1,0.0);
    }
    else{
      DlmPad->SetMargin(0,0.14,0.02,0.0,0.02);
      DlmPad->SetMargin(1,0.14,0.02,0.0,0.0);
      DlmPad->SetMargin(2,0.14,0.02,0.07,0.0);
    }

    DlmPad->cd(0);

    hData_pL_Stat->SetTitle("; #it{k*} (MeV/#it{c}); #it{C}(#it{k*})");
    hData_pL_Stat->GetXaxis()->SetRangeUser(0, 456);
    hData_pL_Stat->GetXaxis()->SetNdivisions(505);
    hData_pL_Stat->GetYaxis()->SetRangeUser(0.9, 2.0);
    hData_pL_Stat->SetFillColor(kGray+1);
    if(PlotsType==1)SetStyleHisto_pLambda(hData_pL_Stat,2,0,1,Panel_X==0?false:true);
    else SetStyleHisto2a(hData_pL_Stat,2,0);

    TGraphErrors *Tgraph_syserror = (TGraphErrors*)gdata_toterr->Clone("Tgraph_syserror");
    TGraph gRelStatError;
    gRelStatError.SetName("gRelStatError");
    TGraph gRelSystError;
    gRelSystError.SetName("gRelSystError");
    TGraphErrors DUMMY_FORCE_ROOT_HAVE_SAME_COLORS_PDFPNF;
    for(unsigned uBin=0; uBin<Tgraph_syserror->GetN(); uBin++){
        double tot_err,stat_err,syst_err,k_err,Ck_val;
        double xVal,yVal;
        gdata_toterr->GetPoint(uBin,xVal,yVal);
        if(DefaultBinValue) Ck_val = hData_pL_Stat->GetBinContent(uBin+1);
        else Ck_val = yVal;
        k_err = Tgraph_syserror->GetErrorX(uBin);
        tot_err = Tgraph_syserror->GetErrorY(uBin);
        stat_err = hData_pL_Stat->GetBinError(uBin+1);
        syst_err = tot_err>stat_err?sqrt(tot_err*tot_err-stat_err*stat_err):0;
//syst_err=tot_err;

        hData_pL_Stat->SetBinContent(uBin+1,Ck_val);
        hData_pL_Stat->SetBinError(uBin+1,stat_err);

        Tgraph_syserror->GetPoint(uBin,xVal,yVal);
        Tgraph_syserror->SetPoint(uBin,xVal,Ck_val);
        Tgraph_syserror->SetPointError(uBin,k_err,syst_err);

        gRelStatError.SetPoint(uBin,xVal,stat_err/Ck_val*100.);
        gRelSystError.SetPoint(uBin,xVal,syst_err/Ck_val*100.);
        DUMMY_FORCE_ROOT_HAVE_SAME_COLORS_PDFPNF.SetPoint(uBin,xVal,-1);
    }
    gRelStatError.SetLineColor(kBlack);
    gRelStatError.SetLineWidth(6);
    gRelStatError.SetMarkerColor(kBlack);
    gRelStatError.SetMarkerStyle(0);
    gRelStatError.SetMarkerSize(0);

    gRelSystError.SetLineColor(kGray+1);
    gRelSystError.SetLineWidth(6);
    gRelSystError.SetMarkerColor(kGray+1);
    gRelSystError.SetMarkerStyle(21);
    gRelSystError.SetMarkerSize(1.6);

    DUMMY_FORCE_ROOT_HAVE_SAME_COLORS_PDFPNF.SetFillColorAlpha(ColorInteraction,0.70);
    DUMMY_FORCE_ROOT_HAVE_SAME_COLORS_PDFPNF.SetLineColor(ColorInteraction);
    DUMMY_FORCE_ROOT_HAVE_SAME_COLORS_PDFPNF.SetLineWidth(3);

    Tgraph_syserror->SetLineColor(kWhite);

    hData_pL_Stat->Draw();
    TString DrawOptions = PlotLines?"3L same":"3 same";
    if(!DataOnly){
        if(XI_PLOT%100) geb_Xim->Draw("3L same");//3L to put a line
        if((XI_PLOT/10)%10 && fitSignal_pL_pXi0) geb_Xi0->Draw("3L same");
        if(XI_PLOT/100) geb_Xi->Draw("3L same");
        if(PlotsType!=1 || SIGMA_FEED!=0) geb_Sig->Draw("3L same");

        if(fitLoDummy) fitLoDummy->Draw("3L same");

        geb_Bl->Draw("3L same");
        geb_Fit->Draw("3L same");
    }

    //}

    hData_pL_Stat->Draw("same");

    Tgraph_syserror->SetFillColorAlpha(kBlack, 0.4);
    Tgraph_syserror->Draw("2 same");
    //hData->Draw("pe same");

    TString LegendSource_line1 = SourceName1;

    unsigned NumRows=DataOnly?3:5;
    TGraph DummyLegendEntry;
    DummyLegendEntry.SetName("DummyLegendEntry");
    DummyLegendEntry.SetLineColor(kWhite);
    DummyLegendEntry.SetMarkerColor(kWhite);
    TGraph DummyLegendEntry2;
    DummyLegendEntry2.SetName("DummyLegendEntry2");
    DummyLegendEntry2.SetLineColor(kWhite);
    DummyLegendEntry2.SetMarkerColor(kWhite);
    const float TextLeft = 0.32;
    const float TextTop = 0.90;
    TLegend *legend;
    if(PlotsType==1 && Panel_X==0) legend = new TLegend(TextLeft+0.01+0.03,TextTop-0.58,0.73,TextTop-0.22);//lbrt
    else if(PlotsType==1 && Panel_X==1) legend = new TLegend(TextLeft+0.01+0.03,TextTop-0.22,0.73,TextTop+0.02);//lbrt
    else if(PlotsType==1 && Panel_X==2) legend = new TLegend(TextLeft+0.01+0.03,TextTop-0.22,0.73,TextTop+0.02);//lbrt
    else if(PlotsType==1 && Panel_X==3) legend = new TLegend(TextLeft+0.01+0.03,TextTop-0.22,0.73,TextTop+0.02);//lbrt
    else if(PlotsType) legend = new TLegend(TextLeft-0.01,0.73-0.064*NumRows,0.73,TextTop-0.135);//lbrt
    else legend = new TLegend(TextLeft-0.01,0.73-0.054*NumRows,0.73,TextTop-0.135);//lbrt
    legend->SetBorderSize(0);
    legend->SetTextFont(42);
    if(PlotsType==1) legend->SetTextSize(gStyle->GetTextSize()*1.3);
    else legend->SetTextSize(gStyle->GetTextSize()*0.90);

    TH1F* hCk_Fake;
    hCk_Fake = (TH1F*)hData_pL_Stat->Clone("hCk_Fake");
    hCk_Fake->SetName("hCk_Fake");
    hCk_Fake->SetLineColor(hCk_Fake->GetFillColor());

    if(PlotsType==1){
      if(Panel_X==0){
        legend->AddEntry(hCk_Fake, "p#minus #Lambda #oplus #bar{p}#minus #bar{#Lambda} pairs", "fpe");
        legend->AddEntry(ge_Fit,TString::Format("Fit %s",PotName1.Data()),"l");
        legend->AddEntry(ge_Sig,SigName1,"l");
      }
      if(Panel_X==1){
        legend->AddEntry(ge_Bl,BlName1,"l");
        legend->AddEntry(ge_Xi,XiName1,"l");
      }
      if(Panel_X==2){
        legend->AddEntry(ge_Fit,TString::Format("Fit %s",PotName1.Data()),"l");
        legend->AddEntry(fitLoDummy,"LO13 (600)","l");
      }
    }
    else{
      legend->AddEntry(hCk_Fake, "p#minus #Lambda #oplus #bar{p}#minus #bar{#Lambda} pairs", "fpe");
      if(!DataOnly){
          //legend->AddEntry(ge_Fit,TString::Format("Fit %s   #chi^{2}_{ndf} = %.1f",PotName1.Data(),MinChi2/double(MinNdf)),"l");
          legend->AddEntry(ge_Fit,TString::Format("Fit %s",PotName1.Data()),"l");
          if(COMPARE_TO_LO&&fitLoDummy) legend->AddEntry(fitLoDummy,"LO13 (600)","l");
          else if(!PlotsType) legend->AddEntry(&DummyLegendEntry,PotName2,"l");
          else legend->AddEntry(&DummyLegendEntry," ","l");
          //legend->AddEntry(&DummyLegendEntry,PotName3,"l");
          legend->AddEntry(ge_Bl,BlName1,"l");
          legend->AddEntry(ge_Sig,SigName1,"l");
          if(XI_PLOT/100) legend->AddEntry(ge_Xi,XiName1,"l");
          else if(XI_PLOT==10) legend->AddEntry(ge_Xi0,Xi0Name1,"l");
          else if(XI_PLOT==1) legend->AddEntry(ge_Xim,XimName1,"l");
      }
      else{
          legend->AddEntry(&gRelStatError,TString::Format("Statistical uncertainty"),"l");
          legend->AddEntry(&gRelSystError,TString::Format("Systematic uncertainty"),"p");
      }
    }

    //legend->AddEntry(&DummyLegendEntry2,BlName2,"l");
    //}

    const double Panel_L = 0.25;
    const double Panel_T = 0.85;
    TLatex Panel_A;
    Panel_A.SetTextSize(gStyle->GetTextSize()*1.6);
    Panel_A.SetNDC(kTRUE);
    Panel_A.SetTextFont(22);
    if(PlotsType==1 && Panel_X==0) Panel_A.DrawLatex(Panel_L, Panel_T, "a)");
    if(PlotsType==1 && Panel_X==1) Panel_A.DrawLatex(Panel_L, Panel_T, "b)");
    if(PlotsType==1 && Panel_X==2) Panel_A.DrawLatex(Panel_L, Panel_T, "c)");
    if(PlotsType==1 && Panel_X==3) Panel_A.DrawLatex(Panel_L, Panel_T, "d)");

    /*
    Panel_A.DrawLatex(Panel_L, Panel_T, "a)");
    TLatex Panel_B;
    Panel_B.SetTextSize(gStyle->GetTextSize()*1.20);
    Panel_B.SetNDC(kTRUE);
    Panel_B.SetTextFont(22);
    Panel_B.DrawLatex(Panel_L, Panel_T, "b)");
    TLatex Panel_C;
    Panel_C.SetTextSize(gStyle->GetTextSize()*1.20);
    Panel_C.SetNDC(kTRUE);
    Panel_C.SetTextFont(22);
    Panel_C.DrawLatex(Panel_L, Panel_T, "c)");
    TLatex Panel_D;
    Panel_D.SetTextSize(gStyle->GetTextSize()*1.20);
    Panel_D.SetNDC(kTRUE);
    Panel_D.SetTextFont(22);
    Panel_D.DrawLatex(Panel_L, Panel_T, "d)");
    */

    if(PlotsType==1){
      if(Panel_X!=3) legend->Draw("same");
    }
    else if(PlotsType!=2 || WhichPotential==11600) legend->Draw("same");

    TLatex BeamText;
    BeamText.SetNDC(kTRUE);
    if(PlotsType==1){
      if(Panel_X==0){
        BeamText.SetTextSize(gStyle->GetTextSize()*1.25);
        BeamText.DrawLatex(TextLeft+0.03, TextTop-0.05, "ALICE");
        BeamText.DrawLatex(TextLeft+0.03, TextTop-0.19, "high-mult. (0#minus 0.17% INEL>0)");
        BeamText.DrawLatex(TextLeft+0.03, TextTop-0.12, "pp #sqrt{#it{s}} = 13 TeV");
      }
    }
    else{
      BeamText.SetTextSize(gStyle->GetTextSize()*0.90);
      if(!PlotsType||WorkInProgress) BeamText.DrawLatex(TextLeft, TextTop, "ALICE work in progress");
      else BeamText.DrawLatex(TextLeft, TextTop, "ALICE");
      BeamText.DrawLatex(TextLeft, TextTop-0.055, "high-mult. (0#minus 0.17% INEL>0) pp #sqrt{#it{s}} = 13 TeV");
    }

    TLatex BeamTextSource;
    BeamTextSource.SetTextSize(gStyle->GetTextSize()*0.90);
    BeamTextSource.SetNDC(kTRUE);
    if(!PlotsType) BeamTextSource.DrawLatex(TextLeft, TextTop-0.11, LegendSource_line1);

    TLatex BeamTextHypCaption;
    BeamTextHypCaption.SetTextSize(gStyle->GetTextSize()*0.90);
    BeamTextHypCaption.SetNDC(kTRUE);
    if(!PlotsType) BeamTextHypCaption.DrawLatex(TextLeft+0.05, TextTop-0.49, "Best fit hypothesis:");

    TLatex BeamTextHyp;
    BeamTextHyp.SetTextSize(gStyle->GetTextSize()*0.80);
    BeamTextHyp.SetNDC(kTRUE);
    //BeamTextHyp.DrawLatex(0.45, 0.35, TString::Format("For k*#in[0, 300] MeV/#it{c}:"));
    if(!PlotsType) BeamTextHyp.DrawLatex(TextLeft+0.1, TextTop-0.54, TString::Format("#chi^{2}/NDF = %.0f/%i = %.1f; n_{#sigma} = %.1f",MinChi2,MinNdf,MinChi2/double(MinNdf),MinNsigma));
    if(!PlotsType) BeamTextHyp.DrawLatex(TextLeft+0.1, TextTop-0.59, TString::Format("r_{core} = %.2f fm",Best_SourceSize));
    if(!PlotsType) BeamTextHyp.DrawLatex(TextLeft+0.1, TextTop-0.64, TString::Format("P(#Lambda) = %.1f%% (#lambda=%.1f%%)",Best_Purity*100,Best_lam_L_genuine*100));
    if(Same_omega_siglam){
        if(!PlotsType) BeamTextHyp.DrawLatex(TextLeft+0.1, TextTop-0.69, TString::Format("#Sigma^{0}:#Lambda = #omega = %.2f; #Xi:(#Lambda+#Sigma^{0}) = %.2f",Best_SigLamFrac,Best_XiSigLamFrac));
    }
    else{
        if(!PlotsType) BeamTextHyp.DrawLatex(TextLeft+0.1, TextTop-0.69, TString::Format("#Sigma^{0}:#Lambda = %.2f; #omega = %.2f;  #Xi:(#Lambda+#Sigma^{0}) = %.2f",Best_SigLamFrac,Best_CuspWeight,Best_XiSigLamFrac));
    }
    if(!PlotsType) BeamTextHyp.DrawLatex(TextLeft+0.1, TextTop-0.74, TString::Format("k*_{c} = %.0f MeV",Best_CkConv));

    //INLET -------------------------------------------------------------------------------------------------------------------

    TH1F* DataHisto_Inlet = (TH1F*)hData_pL_Stat->Clone("DataHisto_Inlet");
    DataHisto_Inlet->SetMarkerSize(hData_pL_Stat->GetMarkerSize()*0.67);
    DataHisto_Inlet->SetLineWidth(hData_pL_Stat->GetLineWidth()*0.67);
    //DataHisto_Inlet->GetXaxis()->SetTitleSize(hData_pL_Stat->GetXaxis()->GetTitleSize()*1.75);
    //DataHisto_Inlet->GetXaxis()->SetLabelSize(hData_pL_Stat->GetXaxis()->GetLabelSize()*1.75);
    DataHisto_Inlet->GetXaxis()->SetRangeUser(0, 456);

    if(PlotsType==1)SetStyleHisto_pLambda(DataHisto_Inlet,2,0,2,Panel_X==0?false:true);
    else SetStyleHisto2a(DataHisto_Inlet,2,0,2);
    DataHisto_Inlet->GetXaxis()->SetNdivisions(505);
    DataHisto_Inlet->GetYaxis()->SetNdivisions(205);

    //DataHisto_Inlet->GetYaxis()->SetTitleSize(hData_pL_Stat->GetYaxis()->GetTitleSize()*1.75);
    //DataHisto_Inlet->GetYaxis()->SetLabelSize(hData_pL_Stat->GetYaxis()->GetLabelSize()*1.75);
    //DataHisto_Inlet->GetYaxis()->SetTitleOffset(hData_pL_Stat->GetYaxis()->GetTitleOffset()*0.67);
    if(*DataSample=="pp13TeV_HM_DimiJun20"||*DataSample=="pp13TeV_HM_DimiJul20") DataHisto_Inlet->GetYaxis()->SetRangeUser(0.965, 1.065);
    else DataHisto_Inlet->GetYaxis()->SetRangeUser(0.985, 1.085);

    TGraph* grFemto_Inlet = (TGraph*)ge_Fit->Clone("grFemto_Inlet");
    grFemto_Inlet->SetLineWidth(ge_Fit->GetLineWidth()*0.67);

    TGraph* grbFemto_Inlet = (TGraph*)geb_Fit->Clone("grbFemto_Inlet");
    grbFemto_Inlet->SetLineWidth(geb_Fit->GetLineWidth()*0.67);

    /*
    const double fXMinInlet=0.17;
    const double fYMinInlet=0.12;
    const double fXMaxInlet=0.97;
    const double fYMaxInlet=0.50;

    TPad *inset_pad = new TPad("insert", "insertPad", fXMinInlet, fYMinInlet,
                                fXMaxInlet, fYMaxInlet);
    inset_pad->SetTopMargin(0.01);
    inset_pad->SetRightMargin(0.05);
    inset_pad->SetBottomMargin(0.28);
    inset_pad->SetLeftMargin(0.28);
    inset_pad->SetFillStyle(4000);
    //if(!DataOnly)
    */
/*
    inset_pad->Draw();
    inset_pad->cd();
    DataHisto_Inlet->Draw();
    //if(!DataOnly){
    ge_Bl->Draw("3 same");
    grFemto_Inlet->Draw("3 same");
    geb_Bl->Draw("3 same");
    grbFemto_Inlet->Draw("3 same");
    //}
    DataHisto_Inlet->Draw("same");
    Tgraph_syserror->Draw("2 same");
*/

    DlmPad->cd(1);
    DataHisto_Inlet->Draw();
    if(!DataOnly){
        if(XI_PLOT%100) geb_Xim->Draw(DrawOptions);
        if((XI_PLOT/10)%10) geb_Xi0->Draw(DrawOptions);
        if(XI_PLOT/100) geb_Xi->Draw(DrawOptions);
        if(PlotsType!=1 || SIGMA_FEED!=0) geb_Sig->Draw(DrawOptions);
        if(fitLoDummy) fitLoDummy->Draw("3L same");
        geb_Bl->Draw(DrawOptions);
        grbFemto_Inlet->Draw(DrawOptions);
    }
    DataHisto_Inlet->Draw("same");
    Tgraph_syserror->Draw("2 same");

    DlmPad->cd(2);
    TH1F* hAxis = new TH1F("hAxis", "hAxis", 456/12, 0, 456);
    hAxis->SetStats(false);
    hAxis->GetXaxis()->SetRangeUser(0, 456);
    if(!DataOnly){
        hAxis->SetTitle("; #it{k*} (MeV/#it{c}); #it{n_{#sigma}}");
        //else //hAxis->SetTitle("; #it{k*} (MeV/#it{c}); #it{C}(#it{k*})");
        hAxis->GetYaxis()->SetRangeUser(-6.2,6.2);
    }
    else{
        DlmPad->SetLogy(2);
        DlmPad->SetGridy(2);
        hAxis->SetTitle("; #it{k*} (MeV/#it{c}); #Delta#it{C}/#it{C} (%)");
        hAxis->GetYaxis()->SetRangeUser(0.06,6);
    }
    hAxis->GetXaxis()->SetNdivisions(505);
    if(PlotsType==1)SetStyleHisto_pLambda(hAxis,2,0,2,Panel_X==0?false:true);
    else SetStyleHisto2a(hAxis,2,0,2);
    hAxis->GetYaxis()->SetNdivisions(504);
    hAxis->Draw("");
    if(!DataOnly) geb_SigmaSlice->Draw(DrawOptions);
    else{
        DUMMY_FORCE_ROOT_HAVE_SAME_COLORS_PDFPNF.Draw(DrawOptions);
        gRelStatError.Draw("L same");
        gRelSystError.Draw("P same");
    }

    if(WhichDataSet<0)
        DlmPad->GetCanvas()->SaveAs(OutputFolder+TString::Format("DlmPad%s_%s_%s_%s_%i.pdf",DataOnly?"Data":"",PotDescr.Data(),BlDescr.Data(),SigDescr.Data(),
                                                             TMath::Nint(ValSourceAlpha*10.)));
    else
        DlmPad->GetCanvas()->SaveAs(OutputFolder+TString::Format("DlmPad%s_%s_%s_%s_%i_DS%i.pdf",DataOnly?"Data":"",PotDescr.Data(),BlDescr.Data(),SigDescr.Data(),
                                                             TMath::Nint(ValSourceAlpha*10.),WhichDataSet));

    float MinChi2NDF = MinChi2/double(MinNdf);

    TFile* InfoFile = NULL;
    TTree* InfoTree = NULL;
    char* InfoFileName = new char [512];
    strcpy(InfoFileName,OutputFolder.Data());
    strcat(InfoFileName,"Info.root");
    int InfoFileStatus = file_status(InfoFileName);
    //the file does not exist
    if(InfoFileStatus==0){
        InfoFile = new TFile(InfoFileName,"recreate");
        InfoTree = new TTree("InfoTree","InfoTree");
        InfoTree->Branch("BASELINE_VAR", &BASELINE_VAR, "BASELINE_VAR/i");//
        InfoTree->Branch("POT_VAR", &POT_VAR, "POT_VAR/I");//
        InfoTree->Branch("Sigma0_Feed", &Sigma0_Feed, "Sigma0_Feed/I");//
        InfoTree->Branch("MinChi2NDF", &MinChi2NDF, "MinChi2NDF/F");
        InfoTree->Branch("MinNsigma", &MinNsigma, "MinNsigma/F");
        InfoTree->Branch("Best_SourceSize", &Best_SourceSize, "Best_SourceSize/F");
        InfoTree->Branch("Best_SourceAlpha", &Best_SourceAlpha, "Best_SourceAlpha/F");
        InfoTree->Branch("Best_Purity", &Best_Purity, "Best_Purity/F");
        InfoTree->Branch("Best_lam_L_genuine", &Best_lam_L_genuine, "Best_lam_L_genuine/F");
        InfoTree->Branch("Best_SigLamFrac", &Best_SigLamFrac, "Best_SigLamFrac/F");
        InfoTree->Branch("Best_XiSigLamFrac", &Best_XiSigLamFrac, "Best_XiSigLamFrac/F");
        InfoTree->Branch("Best_CuspWeight", &Best_CuspWeight, "Best_CuspWeight/F");
        InfoTree->Branch("Best_CkConv", &Best_CkConv, "Best_CkConv/F");
    }
    else{
        DLM_Timer FileTimer;
        long long FileWaitTime;//in micros
        //we wait until the file is closed. If this is not the case in 10s => error
        while(InfoFileStatus==-1&&FileWaitTime<10e6){
            //sleep for 10 ms
            usleep(10e3);
            FileWaitTime = FileTimer.Stop();
            InfoFileStatus = file_status(InfoFileName);
        }
        if(InfoFileStatus==-1){
            printf("\033[1;31mERROR:\033[0m Waited more than 10s to close the file %s\n",InfoFileName);
            abort();
        }

        InfoFile = new TFile(InfoFileName,"update");
        InfoTree = (TTree*)InfoFile->Get("InfoTree");
        InfoTree->SetBranchAddress("BASELINE_VAR",&BASELINE_VAR);
        InfoTree->SetBranchAddress("POT_VAR",&POT_VAR);
        InfoTree->SetBranchAddress("Sigma0_Feed",&Sigma0_Feed);
        InfoTree->SetBranchAddress("MinChi2NDF",&MinChi2NDF);
        InfoTree->SetBranchAddress("MinNsigma",&MinNsigma);
        InfoTree->SetBranchAddress("Best_SourceSize",&Best_SourceSize);
        InfoTree->SetBranchAddress("Best_SourceAlpha",&Best_SourceAlpha);
        InfoTree->SetBranchAddress("Best_Purity",&Best_Purity);
        InfoTree->SetBranchAddress("Best_lam_L_genuine",&Best_lam_L_genuine);
        InfoTree->SetBranchAddress("Best_SigLamFrac",&Best_SigLamFrac);
        InfoTree->SetBranchAddress("Best_XiSigLamFrac",&Best_XiSigLamFrac);
        InfoTree->SetBranchAddress("Best_CuspWeight",&Best_CuspWeight);
        InfoTree->SetBranchAddress("Best_CkConv",&Best_CkConv);
    }

    InfoTree->Fill();

    geb_Xim->Write();
    geb_Sig->Write();
    geb_Bl->Write();
    geb_Fit->Write();
    hData_pL_Stat->Write();
    Tgraph_syserror->Write();
    gRelStatError.Write();
    gRelSystError.Write();
    geb_SigmaSlice->Write();
    InfoTree->Write("",TObject::kOverwrite);

    TH2F* HB_BL = (TH2F*)hb_bl->Clone("HB_BL");
    if(WhichDataSet<0)
        HB_BL->SetName(TString::Format("HB_BL%s_%s_%s_%s_%i",DataOnly?"Data":"",PotDescr.Data(),BlDescr.Data(),SigDescr.Data(),TMath::Nint(ValSourceAlpha*10.)));
    else
        HB_BL->SetName(TString::Format("HB_BL%s_%s_%s_%s_%i_DS%i",DataOnly?"Data":"",PotDescr.Data(),BlDescr.Data(),SigDescr.Data(),
                                                             TMath::Nint(ValSourceAlpha*10.),WhichDataSet));
    HB_BL->Write("",TObject::kOverwrite);
     delete HB_BL;

    TH2F* HB_POL1 = (TH2F*)hb_pol1->Clone("HB_POL1");
    if(WhichDataSet<0)
        HB_POL1->SetName(TString::Format("HB_POL1%s_%s_%s_%s_%i",DataOnly?"Data":"",PotDescr.Data(),BlDescr.Data(),SigDescr.Data(),TMath::Nint(ValSourceAlpha*10.)));
    else
        HB_POL1->SetName(TString::Format("HB_POL1%s_%s_%s_%s_%i_DS%i",DataOnly?"Data":"",PotDescr.Data(),BlDescr.Data(),SigDescr.Data(),
                                                             TMath::Nint(ValSourceAlpha*10.),WhichDataSet));
    HB_POL1->Write("",TObject::kOverwrite);
    delete HB_POL1;



printf("Delete\n");
    delete hnsigma;
    delete hb_nsigma;
    delete hnsigma_0_300;
    delete hnsigma_0_100;
    delete hnsigma_100_200;
    delete hnsigma_200_300;
    delete htot_nsigma_300;
    delete hnsigma_der;
    delete hb_nsigma_der;
    delete hnsigma_der_0_300;
    delete hnsigma_der_0_100;
    delete hnsigma_der_100_200;
    delete hnsigma_der_200_300;
    delete htot_nsigma_der_300;
    delete hfdata;
    delete hbl;
    delete hsig;
    delete hxim;
    delete hb_fdata;
    delete hb_bl;
    delete hb_sig;
    delete hb_xim;
    delete hb_pol1;
printf("Delete 1\n");
    delete ge_SigmaSlice;
    delete ge_SigmaSlice_der;
    delete geb_SigmaSlice;
    delete geb_SigmaSlice_der;
    delete fe_Sigma;
printf("Delete 2\n");
    delete hData_pL_Stat;
    delete ge_Bl;
    delete ge_Sig;
    delete ge_Xim;
    delete ge_Fit;
    delete geb_Bl;
    delete geb_Sig;
    delete geb_Xim;
    delete geb_Fit;
    delete hAxis;
    delete hCk_Fake;
    delete legend;
    delete DataHisto_Inlet;
    delete grFemto_Inlet;
    delete grbFemto_Inlet;
    //delete inset_pad;
    delete DlmPad;
//printf("Delete 3\n");
    delete [] Bin_MinChi2;
//printf("Delete 3a\n");
    //delete [] Chi2_300;
    //delete [] NDF_300;
    //delete [] NSIGMA_300;
    //delete [] IncludePot;
    //delete [] IncludeBl;
    //delete [] BlFlag;

    //delete [] BlName1;
//printf("Delete 3b\n");
    //delete [] BlName2;
//printf("Delete 3c\n");
    //delete [] BlDescr;
//printf("Delete 3d\n");
    //delete [] PotName1;
//printf("Delete 3e\n");
    //delete [] PotName2;
    //delete [] PotFlag;
    //delete [] PotDescr;

    //delete [] SourceName1;
    //delete [] Chi2Ndf;
    delete hdata_toterr;
    delete gdata_toterr;
    delete plambdaTree;
//printf("Delete 4\n");
    if(fitLoDummy) delete fitLoDummy;
    delete InputFile;
    delete [] InfoFileName;
    delete [] GoodCutVar;
    if(InfoFileStatus==0){
        delete InfoTree;
    }
    if(InfoFile) {delete InfoFile; InfoFile=NULL;}
}

void Quick_pLambda_plotter(){

	TString OutputFolder = "/mnt/Ubuntu_Data/CernBox/Sync/pLambda/100720_Unfolded/PlotsERC/";

	TGraphErrors* gfit1;
	TGraphErrors* gfit2;
	TGraphErrors* gfit3;

	TGraphErrors* gnsig1;
	TGraphErrors* gnsig2;
	TGraphErrors* gnsig3;

	TH1F* hData_pL_Stat;
	TGraphErrors *Tgraph_syserror;
	TGraphErrors *gRelStatError;
	TGraphErrors *gRelSystError;

    //geb_Xim->Write();
    //geb_Sig->Write();
    //geb_Bl->Write();
    //geb_Fit->Write();
    //hData_pL_Stat->Write();
	//Tgraph_syserror->Write();
	//gRelStatError.Write();
	//gRelSystError.Write();

	TString InputFileName1 = "/mnt/Ubuntu_Data/CernBox/Sync/pLambda/100720_Unfolded/PlotsERC/Info500.root";
	TFile* InputFile1 = new TFile(InputFileName1,"read");
	gfit1 = (TGraphErrors*)InputFile1->Get("geb_Fit");
	gnsig1 = (TGraphErrors*)InputFile1->Get("geb_SigmaSlice");
	hData_pL_Stat = (TH1F*)InputFile1->Get("hCkS_Norm_12MeV");
	Tgraph_syserror = (TGraphErrors*)InputFile1->Get("Tgraph_syserror");
	gRelStatError = (TGraphErrors*)InputFile1->Get("gRelStatError");
	gRelSystError = (TGraphErrors*)InputFile1->Get("gRelSystError");
  printf("%p\n",gfit1);
  printf("%p\n",gnsig1);
  printf("%p\n",hData_pL_Stat);
	//
	TString InputFileName2 = "/mnt/Ubuntu_Data/CernBox/Sync/pLambda/100720_Unfolded/PlotsERC/Info650.root";
	TFile* InputFile2 = new TFile(InputFileName2,"read");
	gfit2 = (TGraphErrors*)InputFile2->Get("geb_Fit");
	gnsig2 = (TGraphErrors*)InputFile2->Get("geb_SigmaSlice");
	//
	TString InputFileName3 = "/mnt/Ubuntu_Data/CernBox/Sync/pLambda/100720_Unfolded/PlotsERC/Info600.root";
	TFile* InputFile3 = new TFile(InputFileName3,"read");
	gfit3 = (TGraphErrors*)InputFile3->Get("geb_Fit");
	gnsig3 = (TGraphErrors*)InputFile3->Get("geb_SigmaSlice");
	//

    DLM_SubPads* DlmPad = new DLM_SubPads(720,1080);
    DlmPad->AddSubPad(0,1,0.5,1);
    DlmPad->AddSubPad(0,1,0.25,0.5);
    DlmPad->AddSubPad(0,1,0.,0.25);
    DlmPad->SetMargin(0,0.14,0.02,0.0,0.02);
    DlmPad->SetMargin(1,0.14,0.02,0.0,0.0);
    DlmPad->SetMargin(2,0.14,0.02,0.07,0.0);
    DlmPad->cd(0);

    hData_pL_Stat->SetStats(false);
/*
    hData_pL_Stat->SetTitle("; #it{k*} (MeV/#it{c}); #it{C}(#it{k*})");
    hData_pL_Stat->GetXaxis()->SetRangeUser(0, 456);
    hData_pL_Stat->GetXaxis()->SetNdivisions(505);
    hData_pL_Stat->GetYaxis()->SetRangeUser(0.9, 2.0);
    hData_pL_Stat->SetFillColor(kGray+1);
    SetStyleHisto2a(hData_pL_Stat,2,0);

    TGraph gRelStatError;
    gRelStatError.SetName("gRelStatError");
    TGraph gRelSystError;
    gRelSystError.SetName("gRelSystError");
    TGraphErrors DUMMY_FORCE_ROOT_HAVE_SAME_COLORS_PDFPNF;
    for(unsigned uBin=0; uBin<Tgraph_syserror->GetN(); uBin++){
        double tot_err,stat_err,syst_err,k_err,Ck_val;
        double xVal,yVal;
        gdata_toterr->GetPoint(uBin,xVal,yVal);
        if(DefaultBinValue) Ck_val = hData_pL_Stat->GetBinContent(uBin+1);
        else Ck_val = yVal;
        k_err = Tgraph_syserror->GetErrorX(uBin);
        tot_err = Tgraph_syserror->GetErrorY(uBin);
        stat_err = hData_pL_Stat->GetBinError(uBin+1);
        syst_err = tot_err>stat_err?sqrt(tot_err*tot_err-stat_err*stat_err):0;

        hData_pL_Stat->SetBinContent(uBin+1,Ck_val);
        hData_pL_Stat->SetBinError(uBin+1,stat_err);

        Tgraph_syserror->GetPoint(uBin,xVal,yVal);
        Tgraph_syserror->SetPoint(uBin,xVal,Ck_val);
        Tgraph_syserror->SetPointError(uBin,k_err,syst_err);

        gRelStatError.SetPoint(uBin,xVal,stat_err/Ck_val*100.);
        gRelSystError.SetPoint(uBin,xVal,syst_err/Ck_val*100.);
        DUMMY_FORCE_ROOT_HAVE_SAME_COLORS_PDFPNF.SetPoint(uBin,xVal,-1);
    }
    gRelStatError.SetLineColor(kBlack);
    gRelStatError.SetLineWidth(6);
    gRelStatError.SetMarkerColor(kBlack);
    gRelStatError.SetMarkerStyle(0);
    gRelStatError.SetMarkerSize(0);

    gRelSystError.SetLineColor(kGray+1);
    gRelSystError.SetLineWidth(6);
    gRelSystError.SetMarkerColor(kGray+1);
    gRelSystError.SetMarkerStyle(21);
    gRelSystError.SetMarkerSize(1.6);

    DUMMY_FORCE_ROOT_HAVE_SAME_COLORS_PDFPNF.SetFillColorAlpha(ColorInteraction,0.70);
    DUMMY_FORCE_ROOT_HAVE_SAME_COLORS_PDFPNF.SetLineColor(ColorInteraction);
    DUMMY_FORCE_ROOT_HAVE_SAME_COLORS_PDFPNF.SetLineWidth(3);

    Tgraph_syserror->SetLineColor(kWhite);
*/
//InputFile1->cd();
    hData_pL_Stat->Draw();
    //TString DrawOptions = PlotLines?"3L same":"3 same";
	TString DrawOptions = "3 same";

  gfit1->SetFillColorAlpha(kBlue+1,0.70);
  gfit1->SetLineColor(kBlue+1);
  gfit1->SetLineWidth(3);

  gfit2->SetFillColorAlpha(kYellow+1,0.70);
  gfit2->SetLineColor(kYellow+1);
  gfit2->SetLineWidth(3);

  gfit3->SetFillColorAlpha(kRed+1,0.70);
  gfit3->SetLineColor(kRed+1);
  gfit3->SetLineWidth(3);

  gnsig1->SetFillColorAlpha(kBlue+1,0.70);
  gnsig1->SetLineColor(kBlue+1);
  gnsig1->SetLineWidth(3);

  gnsig2->SetFillColorAlpha(kYellow+1,0.70);
  gnsig2->SetLineColor(kYellow+1);
  gnsig2->SetLineWidth(3);

  gnsig3->SetFillColorAlpha(kRed+1,0.70);
  gnsig3->SetLineColor(kRed+1);
  gnsig3->SetLineWidth(3);



	gfit1->Draw("3L same");
//return;
	gfit2->Draw("3L same");
	gfit3->Draw("3L same");
    hData_pL_Stat->Draw("same");

    Tgraph_syserror->SetFillColorAlpha(kBlack, 0.4);
    Tgraph_syserror->Draw("2 same");

    TString LegendSource_line1 = "source name";

    unsigned NumRows=5;
    TGraph DummyLegendEntry;
    DummyLegendEntry.SetName("DummyLegendEntry");
    DummyLegendEntry.SetLineColor(kWhite);
    DummyLegendEntry.SetMarkerColor(kWhite);
    TGraph DummyLegendEntry2;
    DummyLegendEntry2.SetName("DummyLegendEntry2");
    DummyLegendEntry2.SetLineColor(kWhite);
    DummyLegendEntry2.SetMarkerColor(kWhite);
    const float TextLeft = 0.32;
    const float TextTop = 0.90;
    TLegend *legend;
    //if(PlotsType) legend = new TLegend(TextLeft-0.01,0.73-0.064*NumRows,0.73,TextTop-0.135);//lbrt
    legend = new TLegend(TextLeft-0.01,0.73-0.054*NumRows,0.73,TextTop-0.135);//lbrt
    legend->SetBorderSize(0);
    legend->SetTextFont(42);
    legend->SetTextSize(gStyle->GetTextSize()*0.90);
    TH1F* hCk_Fake;
    hCk_Fake = (TH1F*)hData_pL_Stat->Clone("hCk_Fake");
    hCk_Fake->SetName("hCk_Fake");
    hCk_Fake->SetLineColor(hCk_Fake->GetFillColor());

    legend->AddEntry(hCk_Fake, "p#minus #Lambda #oplus #bar{p}#minus #bar{#Lambda} pairs", "fpe");
    legend->AddEntry(gfit3, "NLO19-600");
    legend->AddEntry(gfit1, "NLO19-500");
    legend->AddEntry(gfit2, "NLO19-650");
/*
    if(!DataOnly){
        legend->AddEntry(ge_Fit,TString::Format("Fit %s",PotName1.Data()),"l");
        if(COMPARE_TO_LO&&fitLoDummy) legend->AddEntry(fitLoDummy,"LO13 (600)","l");
        else if(!PlotsType) legend->AddEntry(&DummyLegendEntry,PotName2,"l");
        else legend->AddEntry(&DummyLegendEntry," ","l");
        //legend->AddEntry(&DummyLegendEntry,PotName3,"l");
        legend->AddEntry(ge_Bl,BlName1,"l");
        legend->AddEntry(ge_Sig,SigName1,"l");
        if(XI_PLOT/100) legend->AddEntry(ge_Xi,XiName1,"l");
        else if(XI_PLOT==10) legend->AddEntry(ge_Xi0,Xi0Name1,"l");
        else if(XI_PLOT==1) legend->AddEntry(ge_Xim,XimName1,"l");
    }
    else{
        legend->AddEntry(&gRelStatError,TString::Format("Statistical uncertainty"),"l");
        legend->AddEntry(&gRelSystError,TString::Format("Systematic uncertainty"),"p");
    }
*/
    legend->Draw("same");
    TLatex BeamText;
    BeamText.SetTextSize(gStyle->GetTextSize()*0.90);
    BeamText.SetNDC(kTRUE);
    //if(!PlotsType)
	BeamText.DrawLatex(TextLeft, TextTop, "ALICE work in progress");
    //else BeamText.DrawLatex(TextLeft, TextTop, "ALICE");
    BeamText.DrawLatex(TextLeft, TextTop-0.055, "high-mult. (0#minus 0.17% INEL>0) pp #sqrt{#it{s}} = 13 TeV");

    TLatex BeamTextSource;
    BeamTextSource.SetTextSize(gStyle->GetTextSize()*0.90);
    BeamTextSource.SetNDC(kTRUE);
    //if(!PlotsType)
	//BeamTextSource.DrawLatex(TextLeft, TextTop-0.11, LegendSource_line1);

    TLatex BeamTextHypCaption;
    BeamTextHypCaption.SetTextSize(gStyle->GetTextSize()*0.90);
    BeamTextHypCaption.SetNDC(kTRUE);
    //if(!PlotsType)
	//BeamTextHypCaption.DrawLatex(TextLeft+0.05, TextTop-0.49, "Best fit hypothesis:");

    //TLatex BeamTextHyp;
    //BeamTextHyp.SetTextSize(gStyle->GetTextSize()*0.80);
    //BeamTextHyp.SetNDC(kTRUE);
    //if(!PlotsType) BeamTextHyp.DrawLatex(TextLeft+0.1, TextTop-0.54, TString::Format("#chi^{2}/NDF = %.0f/%i = %.1f; n_{#sigma} = %.1f",MinChi2,MinNdf,MinChi2/double(MinNdf),MinNsigma));
    //if(!PlotsType) BeamTextHyp.DrawLatex(TextLeft+0.1, TextTop-0.59, TString::Format("r_{core} = %.2f fm",Best_SourceSize));
    //if(!PlotsType) BeamTextHyp.DrawLatex(TextLeft+0.1, TextTop-0.64, TString::Format("P(#Lambda) = %.1f%% (#lambda=%.1f%%)",Best_Purity*100,Best_lam_L_genuine*100));
    //if(Same_omega_siglam){
    //    if(!PlotsType) BeamTextHyp.DrawLatex(TextLeft+0.1, TextTop-0.69, TString::Format("#Sigma^{0}:#Lambda = #omega = %.2f; #Xi:(#Lambda+#Sigma^{0}) = %.2f",Best_SigLamFrac,Best_XiSigLamFrac));
    //}
    //else{
    //    if(!PlotsType) BeamTextHyp.DrawLatex(TextLeft+0.1, TextTop-0.69, TString::Format("#Sigma^{0}:#Lambda = %.2f; #omega = %.2f;  #Xi:(#Lambda+#Sigma^{0}) = %.2f",Best_SigLamFrac,Best_CuspWeight,Best_XiSigLamFrac));
    //}
    //if(!PlotsType) BeamTextHyp.DrawLatex(TextLeft+0.1, TextTop-0.74, TString::Format("k*_{c} = %.0f MeV",Best_CkConv));

    //INLET -------------------------------------------------------------------------------------------------------------------

    TH1F* DataHisto_Inlet = (TH1F*)hData_pL_Stat->Clone("DataHisto_Inlet");
    DataHisto_Inlet->SetMarkerSize(hData_pL_Stat->GetMarkerSize()*0.67);
    DataHisto_Inlet->SetLineWidth(hData_pL_Stat->GetLineWidth()*0.67);
    DataHisto_Inlet->GetXaxis()->SetRangeUser(0, 456);
    SetStyleHisto2a(DataHisto_Inlet,2,0,2);
    DataHisto_Inlet->GetXaxis()->SetNdivisions(505);
    DataHisto_Inlet->GetYaxis()->SetNdivisions(205);

	DataHisto_Inlet->GetYaxis()->SetRangeUser(0.965, 1.065);

    TGraph* grbFemto_Inlet1 = (TGraph*)gfit1->Clone("grbFemto_Inlet1");
    grbFemto_Inlet1->SetLineWidth(gfit1->GetLineWidth()*0.67);
    TGraph* grbFemto_Inlet2 = (TGraph*)gfit2->Clone("grbFemto_Inlet2");
    grbFemto_Inlet2->SetLineWidth(gfit2->GetLineWidth()*0.67);
    TGraph* grbFemto_Inlet3 = (TGraph*)gfit3->Clone("grbFemto_Inlet3");
    grbFemto_Inlet3->SetLineWidth(gfit3->GetLineWidth()*0.67);

    DlmPad->cd(1);
    DataHisto_Inlet->SetStats(false);
    DataHisto_Inlet->Draw();
	grbFemto_Inlet1->Draw(DrawOptions);
	grbFemto_Inlet2->Draw(DrawOptions);
	grbFemto_Inlet3->Draw(DrawOptions);
    DataHisto_Inlet->Draw("same");
    Tgraph_syserror->Draw("2 same");

    DlmPad->cd(2);
    TH1F* hAxis = new TH1F("hAxis", "hAxis", 456/12, 0, 456);
    hAxis->SetStats(false);
    hAxis->GetXaxis()->SetRangeUser(0, 456);
    //if(!DataOnly){
        hAxis->SetTitle("; #it{k*} (MeV/#it{c}); #it{n_{#sigma}}");
        //else //hAxis->SetTitle("; #it{k*} (MeV/#it{c}); #it{C}(#it{k*})");
        hAxis->GetYaxis()->SetRangeUser(-6.2,6.2);
    //}
    //else{
    //    DlmPad->SetLogy(2);
    //    DlmPad->SetGridy(2);
    //    hAxis->SetTitle("; #it{k*} (MeV/#it{c}); #Delta#it{C}/#it{C} (%)");
    //    hAxis->GetYaxis()->SetRangeUser(0.06,6);
    //}
    hAxis->GetXaxis()->SetNdivisions(505);
    SetStyleHisto2a(hAxis,2,0,2);
    hAxis->GetYaxis()->SetNdivisions(504);
    hAxis->Draw("");
    //if(!DataOnly){
		gnsig1->Draw(DrawOptions);
		gnsig2->Draw(DrawOptions);
		gnsig3->Draw(DrawOptions);
	//}
    //else{
    //    DUMMY_FORCE_ROOT_HAVE_SAME_COLORS_PDFPNF.Draw(DrawOptions);
    //    gRelStatError.Draw("L same");
    //    gRelSystError.Draw("P same");
    //}

    //if(WhichDataSet<0)
        DlmPad->GetCanvas()->SaveAs(OutputFolder+TString::Format("DlmQPad.pdf"));
    //else
    //    DlmPad->GetCanvas()->SaveAs(OutputFolder+TString::Format("DlmQPad%s_%s_%s_%s_%i_DS%i.pdf",DataOnly?"Data":"",PotDescr.Data(),BlDescr.Data(),SigDescr.Data(),
    //                                                         TMath::Nint(ValSourceAlpha*10.),WhichDataSet));

}

void Quick_pLambda_plotter_NLO13_vs_LO13(){

	TString OutputFolder = "/mnt/Ubuntu_Data/CernBox/Sync/pLambda/100720_Unfolded/Laura_Review/";

	TGraphErrors* gfit1;
	TGraphErrors* gfit2;

	TH1F* hData_pL_Stat;
	TGraphErrors *Tgraph_syserror;
	TGraphErrors *gRelStatError;
	TGraphErrors *gRelSystError;

	TString InputFileName1 = "/mnt/Ubuntu_Data/CernBox/Sync/pLambda/100720_Unfolded/Laura_Review/InfoNLO13.root";
	TFile* InputFile1 = new TFile(InputFileName1,"read");
	gfit1 = (TGraphErrors*)InputFile1->Get("geb_Fit");
	hData_pL_Stat = (TH1F*)InputFile1->Get("hCkS_Norm_12MeV");
	Tgraph_syserror = (TGraphErrors*)InputFile1->Get("Tgraph_syserror");
	gRelStatError = (TGraphErrors*)InputFile1->Get("gRelStatError");
	gRelSystError = (TGraphErrors*)InputFile1->Get("gRelSystError");
  printf("%p\n",gfit1);
  printf("%p\n",hData_pL_Stat);
	//
	TString InputFileName2 = "/mnt/Ubuntu_Data/CernBox/Sync/pLambda/100720_Unfolded/Laura_Review/InfoLO13.root";
	TFile* InputFile2 = new TFile(InputFileName2,"read");
	gfit2 = (TGraphErrors*)InputFile2->Get("geb_Fit");
	//

  DLM_SubPads* DlmPad = new DLM_SubPads(1280,960);
  DlmPad->AddSubPad(0,1,0,1);
  DlmPad->SetMargin(0,0.14,0.02,0.13,0.02);//lrbt
  DlmPad->cd(0);

  hData_pL_Stat->SetStats(false);

  hData_pL_Stat->Draw();

	TString DrawOptions = "3 same";

  gfit1->SetFillColorAlpha(kRed+1,0.70);
  gfit1->SetLineColor(kRed+1);
  gfit1->SetLineWidth(3);

  gfit2->SetFillColorAlpha(kGreen+1,0.70);
  gfit2->SetLineColor(kGreen+1);
  gfit2->SetLineWidth(3);

	gfit1->Draw("3L same");
	gfit2->Draw("3L same");
  hData_pL_Stat->Draw("same");

  Tgraph_syserror->SetFillColorAlpha(kBlack, 0.4);
  Tgraph_syserror->Draw("2 same");

  TString LegendSource_line1 = "source name";

  unsigned NumRows=3;
  TGraph DummyLegendEntry;
  DummyLegendEntry.SetName("DummyLegendEntry");
  DummyLegendEntry.SetLineColor(kWhite);
  DummyLegendEntry.SetMarkerColor(kWhite);
  TGraph DummyLegendEntry2;
  DummyLegendEntry2.SetName("DummyLegendEntry2");
  DummyLegendEntry2.SetLineColor(kWhite);
  DummyLegendEntry2.SetMarkerColor(kWhite);
  const float TextLeft = 0.32;
  const float TextTop = 0.93;
  TLegend *legend;
  //if(PlotsType) legend = new TLegend(TextLeft-0.01,0.73-0.064*NumRows,0.73,TextTop-0.135);//lbrt
  legend = new TLegend(TextLeft-0.01,0.74-0.050*NumRows,0.73,TextTop-0.095);//lbrt
  legend->SetBorderSize(0);
  legend->SetTextFont(42);
  legend->SetTextSize(gStyle->GetTextSize()*0.85);
  TH1F* hCk_Fake;
  hCk_Fake = (TH1F*)hData_pL_Stat->Clone("hCk_Fake");
  hCk_Fake->SetName("hCk_Fake");
  hCk_Fake->SetLineColor(hCk_Fake->GetFillColor());

  legend->AddEntry(hCk_Fake, "p#minus #Lambda #oplus #bar{p}#minus #bar{#Lambda} pairs", "fpe");
  legend->AddEntry(gfit1, "NLO13-600");
  legend->AddEntry(gfit2, "LO13-600");

  legend->Draw("same");
  TLatex BeamText;
  BeamText.SetTextSize(gStyle->GetTextSize()*0.90);
  BeamText.SetNDC(kTRUE);
	BeamText.DrawLatex(TextLeft, TextTop, "ALICE work in progress");
  BeamText.DrawLatex(TextLeft, TextTop-0.055, "high-mult. (0#minus 0.17% INEL>0) pp #sqrt{#it{s}} = 13 TeV");

  TLatex BeamTextSource;
  BeamTextSource.SetTextSize(gStyle->GetTextSize()*0.90);
  BeamTextSource.SetNDC(kTRUE);

  TLatex BeamTextHypCaption;
  BeamTextHypCaption.SetTextSize(gStyle->GetTextSize()*0.90);
  BeamTextHypCaption.SetNDC(kTRUE);



///////////////////////////

//INLET -------------------------------------------------------------------------------------------------------------------

  TH1F* DataHisto_Inlet = (TH1F*)hData_pL_Stat->Clone("DataHisto_Inlet");
  DataHisto_Inlet->SetMarkerSize(hData_pL_Stat->GetMarkerSize()*0.67);
  DataHisto_Inlet->SetLineWidth(hData_pL_Stat->GetLineWidth()*0.67);
  DataHisto_Inlet->GetXaxis()->SetTitleSize(hData_pL_Stat->GetXaxis()->GetTitleSize()*1.75);
  DataHisto_Inlet->GetXaxis()->SetLabelSize(hData_pL_Stat->GetXaxis()->GetLabelSize()*1.75);
  DataHisto_Inlet->GetXaxis()->SetRangeUser(100, 336);
  DataHisto_Inlet->GetXaxis()->SetNdivisions(505);

  DataHisto_Inlet->GetYaxis()->SetTitleSize(hData_pL_Stat->GetYaxis()->GetTitleSize()*1.75);
  DataHisto_Inlet->GetYaxis()->SetLabelSize(hData_pL_Stat->GetYaxis()->GetLabelSize()*1.75);
  DataHisto_Inlet->GetYaxis()->SetTitleOffset(hData_pL_Stat->GetYaxis()->GetTitleOffset()*0.67);
  DataHisto_Inlet->GetYaxis()->SetRangeUser(0.97, 1.025);

  //TGraph* grFemto_Inlet = (TGraph*)graph_FitMAX[NumBlTypes].Clone("grFemto_Inlet");
  //grFemto_Inlet->SetLineWidth(graph_FitMAX[NumBlTypes].GetLineWidth()*0.67);

  const double fXMinInlet=0.30;
  const double fYMinInlet=0.21;
  const double fXMaxInlet=0.95;
  const double fYMaxInlet=0.59;
  TPad *inset_pad = new TPad("insert", "insertPad", fXMinInlet, fYMinInlet,
                           fXMaxInlet, fYMaxInlet);
  inset_pad->SetTopMargin(0.01);
  inset_pad->SetRightMargin(0.05);
  inset_pad->SetBottomMargin(0.28);
  inset_pad->SetLeftMargin(0.28);
  inset_pad->SetFillStyle(4000);
  inset_pad->Draw();
  inset_pad->cd();
  DataHisto_Inlet->Draw();

  gfit1->Draw("3 same");
	gfit2->Draw("3 same");
  Tgraph_syserror->Draw("2 same");


///////////////////////////////


    //INLET -------------------------------------------------------------------------------------------------------------------
/*
  TH1F* DataHisto_Inlet = (TH1F*)hData_pL_Stat->Clone("DataHisto_Inlet");
  DataHisto_Inlet->SetMarkerSize(hData_pL_Stat->GetMarkerSize()*0.67);
  DataHisto_Inlet->SetLineWidth(hData_pL_Stat->GetLineWidth()*0.67);
  DataHisto_Inlet->GetXaxis()->SetRangeUser(0, 456);
  SetStyleHisto2a(DataHisto_Inlet,2,0,2);
  DataHisto_Inlet->GetXaxis()->SetNdivisions(505);
  DataHisto_Inlet->GetYaxis()->SetNdivisions(205);

	DataHisto_Inlet->GetYaxis()->SetRangeUser(0.965, 1.065);

    TGraph* grbFemto_Inlet1 = (TGraph*)gfit1->Clone("grbFemto_Inlet1");
    grbFemto_Inlet1->SetLineWidth(gfit1->GetLineWidth()*0.67);
    TGraph* grbFemto_Inlet2 = (TGraph*)gfit2->Clone("grbFemto_Inlet2");
    grbFemto_Inlet2->SetLineWidth(gfit2->GetLineWidth()*0.67);
    TGraph* grbFemto_Inlet3 = (TGraph*)gfit3->Clone("grbFemto_Inlet3");
    grbFemto_Inlet3->SetLineWidth(gfit3->GetLineWidth()*0.67);

    DlmPad->cd(1);
    DataHisto_Inlet->SetStats(false);
    DataHisto_Inlet->Draw();
	grbFemto_Inlet1->Draw(DrawOptions);
	grbFemto_Inlet2->Draw(DrawOptions);
	grbFemto_Inlet3->Draw(DrawOptions);
    DataHisto_Inlet->Draw("same");
    Tgraph_syserror->Draw("2 same");

    DlmPad->cd(2);
    TH1F* hAxis = new TH1F("hAxis", "hAxis", 456/12, 0, 456);
    hAxis->SetStats(false);
    hAxis->GetXaxis()->SetRangeUser(0, 456);
    //if(!DataOnly){
        hAxis->SetTitle("; #it{k*} (MeV/#it{c}); #it{n_{#sigma}}");
        //else //hAxis->SetTitle("; #it{k*} (MeV/#it{c}); #it{C}(#it{k*})");
        hAxis->GetYaxis()->SetRangeUser(-6.2,6.2);
    //}
    //else{
    //    DlmPad->SetLogy(2);
    //    DlmPad->SetGridy(2);
    //    hAxis->SetTitle("; #it{k*} (MeV/#it{c}); #Delta#it{C}/#it{C} (%)");
    //    hAxis->GetYaxis()->SetRangeUser(0.06,6);
    //}
    hAxis->GetXaxis()->SetNdivisions(505);
    SetStyleHisto2a(hAxis,2,0,2);
    hAxis->GetYaxis()->SetNdivisions(504);
    hAxis->Draw("");
    //if(!DataOnly){
		gnsig1->Draw(DrawOptions);
		gnsig2->Draw(DrawOptions);
		gnsig3->Draw(DrawOptions);
	//}
    //else{
    //    DUMMY_FORCE_ROOT_HAVE_SAME_COLORS_PDFPNF.Draw(DrawOptions);
    //    gRelStatError.Draw("L same");
    //    gRelSystError.Draw("P same");
    //}
*/

  DlmPad->GetCanvas()->SaveAs(OutputFolder+TString::Format("DlmGhettoPlot.pdf"));


}



//with the chi2 from the Info.root of the above function
//Compact = no ESC16
void MakeLATEXtable(TString InputFolder, bool Compact=false){
    TString OutputFileName;
    TString OutputFileNameNsig;
    TString OutputFileNameSigDev;
    if(Compact) OutputFileName = InputFolder+"Chi2LaTeX_compact.tex";
    else OutputFileName = InputFolder+"Chi2LaTeX_full.tex";
    if(Compact) OutputFileNameNsig = InputFolder+"nSigLaTeX_compact.tex";
    else OutputFileNameNsig = InputFolder+"nSigLaTeX_full.tex";
    if(Compact) OutputFileNameSigDev = InputFolder+"SigDevLaTeX_compact.tex";
    else OutputFileNameSigDev = InputFolder+"SigDevLaTeX_full.tex";

    char* InfoFileName = new char [512];
    strcpy(InfoFileName,InputFolder.Data());
    strcat(InfoFileName,"Info.root");

    const unsigned NumLamVars = 9;
    TString* PotName_pL = new TString [NumLamVars];
    PotName_pL[0] = "LO13-600";
    PotName_pL[1] = "NLO13-500";
    PotName_pL[2] = "NLO13-550";
    PotName_pL[3] = "NLO13-600";
    PotName_pL[4] = "NLO13-650";
    PotName_pL[5] = "NLO19-500";
    PotName_pL[6] = "NLO19-550";
    PotName_pL[7] = "NLO19-600";
    PotName_pL[8] = "NLO19-650";
    double BestChi2_const = 1e6;
    double BestChi2_cubic = 1e6;

    int* PotFlag_pL = new int [NumLamVars];
    PotFlag_pL[0] = -11600;
    PotFlag_pL[1] = 1500;
    PotFlag_pL[2] = 1550;
    PotFlag_pL[3] = 1600;
    PotFlag_pL[4] = 1650;
    PotFlag_pL[5] = 11500;
    PotFlag_pL[6] = 11550;
    PotFlag_pL[7] = 11600;
    PotFlag_pL[8] = 11650;

    //const unsigned NumSigVars = 3;
    //TString* PotName_pS0 = new TString [NumSigVars];
    //PotName_pS0[0] = "\\pmb{\\chiEFT}";
    //PotName_pS0[1] = "ESC16";
    //PotName_pS0[2] = "Flat";
    //int* PotFlag_pS0 = new int [NumSigVars];
    //PotFlag_pS0[0] = 1;
    //PotFlag_pS0[1] = 2;
    //PotFlag_pS0[2] = 0;

    unsigned BASELINE_VAR;
    int POT_VAR;
    int Sigma0_Feed;
//! hard coded for 12 MeV bins up to 300 MeV
    const float NDF = 25;
    float MinChi2NDF;
    float MinNsigma;
    float Best_SourceSize;
    float Best_SourceAlpha;
    float Best_Purity;
    float Best_lam_L_genuine;
    float Best_SigLamFrac;
    float Best_XiSigLamFrac;
    float Best_CuspWeight;
    float Best_CkConv;

    TFile* InfoFile = new TFile(InfoFileName,"read");
    TTree* InfoTree = (TTree*)InfoFile->Get("InfoTree");
    InfoTree->SetBranchAddress("BASELINE_VAR",&BASELINE_VAR);
    InfoTree->SetBranchAddress("POT_VAR",&POT_VAR);
    InfoTree->SetBranchAddress("Sigma0_Feed",&Sigma0_Feed);
    InfoTree->SetBranchAddress("MinChi2NDF",&MinChi2NDF);
    InfoTree->SetBranchAddress("MinNsigma",&MinNsigma);
    InfoTree->SetBranchAddress("Best_SourceSize",&Best_SourceSize);
    InfoTree->SetBranchAddress("Best_SourceAlpha",&Best_SourceAlpha);
    InfoTree->SetBranchAddress("Best_Purity",&Best_Purity);
    InfoTree->SetBranchAddress("Best_lam_L_genuine",&Best_lam_L_genuine);
    InfoTree->SetBranchAddress("Best_SigLamFrac",&Best_SigLamFrac);
    InfoTree->SetBranchAddress("Best_XiSigLamFrac",&Best_XiSigLamFrac);
    InfoTree->SetBranchAddress("Best_CuspWeight",&Best_CuspWeight);
    InfoTree->SetBranchAddress("Best_CkConv",&Best_CkConv);
    unsigned NumEntries = InfoTree->GetEntries();

    FILE *fptr;
    fptr=fopen(OutputFileName.Data(),"w");
    if(Compact){
        fprintf(fptr,"\\begin{table} \n");
        fprintf(fptr,"\\begin{center} \n");
        fprintf(fptr,"\\begin{tabular}{r|c|c|} \n");
        fprintf(fptr,"\\hline \\hline \n");
        fprintf(fptr,"$\\text{p--}\\Sigma^0$ ($\\rightarrow$)& \\pmb{$\\chi$EFT} & Flat \\\\ $\\text{p--}\\Lambda$ ($\\downarrow$) & & \\\\ \n");
        fprintf(fptr,"\\hline \n");
    }
    else{
        fprintf(fptr,"\\begin{table} \n");
        fprintf(fptr,"\\begin{center} \n");
        fprintf(fptr,"\\begin{tabular}{r|c|c|c|} \n");
        fprintf(fptr,"\\hline \\hline \n");
        fprintf(fptr,"$\\text{p--}\\Sigma^0$ ($\\rightarrow$)& \\pmb{$\\chi$EFT} & ESC16 & Flat \\\\ $\\text{p--}\\Lambda$ ($\\downarrow$) & & & \\\\ \n");
        fprintf(fptr,"\\hline \n");
    }

    //char* buffer = new char [512];
    for(unsigned upL=0; upL<NumLamVars; upL++){
        float val_chiral_cubic=0;
        float val_esc16_cubic=0;
        float val_flat_cubic=0;
        float val_chiral_const=0;
        float val_esc16_const=0;
        float val_flat_const=0;
        for(unsigned uEntry=0; uEntry<NumEntries; uEntry++){
            InfoTree->GetEntry(uEntry);
            if(POT_VAR!=PotFlag_pL[upL]) continue;
            if(Sigma0_Feed==0&&BASELINE_VAR==0){
                val_flat_const = MinChi2NDF;
            }
            else if(Sigma0_Feed==0&&BASELINE_VAR==10){
                val_flat_cubic = MinChi2NDF;
            }
            else if(Sigma0_Feed==1&&BASELINE_VAR==0){
                val_chiral_const = MinChi2NDF;
            }
            else if(Sigma0_Feed==1&&BASELINE_VAR==10){
                val_chiral_cubic = MinChi2NDF;
            }
            else if(Sigma0_Feed==2&&BASELINE_VAR==0){
                val_esc16_const = MinChi2NDF;
            }
            else if(Sigma0_Feed==2&&BASELINE_VAR==10){
                val_esc16_cubic = MinChi2NDF;
            }
        }

        if(PotFlag_pL[upL]==11600) fprintf(fptr,"\\pmb{");
        fprintf(fptr,"%s",PotName_pL[upL].Data());
        if(PotFlag_pL[upL]==11600) fprintf(fptr,"}");
        if(Compact){
            fprintf(fptr," & %.1f (%.1f) & %.1f (%.1f)",
                    val_chiral_cubic,val_chiral_const,
                    val_flat_cubic,val_flat_const);
        }
        else{
            fprintf(fptr," & %.1f (%.1f) & %.1f (%.1f) & %.1f (%.1f)",
                    val_chiral_cubic,val_chiral_const,
                    val_esc16_cubic,val_esc16_const,
                    val_flat_cubic,val_flat_const);
        }

        fprintf(fptr," \\\\ \n");

        if(BestChi2_const>val_chiral_const*NDF) BestChi2_const=val_chiral_const*NDF;
        if(BestChi2_const>val_flat_const*NDF) BestChi2_const=val_flat_const*NDF;

        if(BestChi2_cubic>val_chiral_cubic*NDF) BestChi2_cubic=val_chiral_cubic*NDF;
        if(BestChi2_cubic>val_flat_cubic*NDF) BestChi2_cubic=val_flat_cubic*NDF;
    }
    fprintf(fptr,"\\hline \\hline \n");
    fprintf(fptr,"\\end{tabular} \n");
    fprintf(fptr,"\\caption{Values of $\\chi^2/$NDF for the different interaction hypotheses "
            "of $\\text{p--}\\Lambda$ and $\\text{p--}\\Sigma^0$, evaluated for $k^{*}\\in[0,300]~$MeV. "
            "The default values correspond to the fit with a cubic baseline, the values in brackets represent the results from using a constant baseline. "
            "The default model (in bold) is the $\\chi$EFT NLO calculation, at a cutoff of 600~MeV. } \n");
    fprintf(fptr,"\\label{tab:chi2_2020} \n");
    fprintf(fptr,"\\end{center} \n");
    fprintf(fptr,"\\end{table} \n");

    if(fptr) fclose(fptr);

    printf("BestChi2_const = %.1f (%.1f)\n",BestChi2_const,BestChi2_const/NDF);
    printf("BestChi2_cubic = %.1f (%.1f)\n",BestChi2_cubic,BestChi2_cubic/NDF);


    fptr=fopen(OutputFileNameNsig.Data(),"w");
    if(Compact){
        fprintf(fptr,"\\begin{table} \n");
        fprintf(fptr,"\\begin{center} \n");
        fprintf(fptr,"\\begin{tabular}{r|c|c|} \n");
        fprintf(fptr,"\\hline \\hline \n");
        fprintf(fptr,"$\\text{p--}\\Sigma^0$ ($\\rightarrow$)& \\pmb{$\\chi$EFT} & Flat \\\\ $\\text{p--}\\Lambda$ ($\\downarrow$) & & \\\\ \n");
        fprintf(fptr,"\\hline \n");
    }
    else{
        fprintf(fptr,"\\begin{table} \n");
        fprintf(fptr,"\\begin{center} \n");
        fprintf(fptr,"\\begin{tabular}{r|c|c|c|} \n");
        fprintf(fptr,"\\hline \\hline \n");
        fprintf(fptr,"$\\text{p--}\\Sigma^0$ ($\\rightarrow$)& \\pmb{$\\chi$EFT} & ESC16 & Flat \\\\ $\\text{p--}\\Lambda$ ($\\downarrow$) & & & \\\\ \n");
        fprintf(fptr,"\\hline \n");
    }

    //char* buffer = new char [512];
    for(unsigned upL=0; upL<NumLamVars; upL++){
        float val_chiral_cubic=0;
        float val_esc16_cubic=0;
        float val_flat_cubic=0;
        float val_chiral_const=0;
        float val_esc16_const=0;
        float val_flat_const=0;

        float nsig_chiral_cubic=0;
        float nsig_esc16_cubic=0;
        float nsig_flat_cubic=0;
        float nsig_chiral_const=0;
        float nsig_esc16_const=0;
        float nsig_flat_const=0;

        for(unsigned uEntry=0; uEntry<NumEntries; uEntry++){
            InfoTree->GetEntry(uEntry);
            if(POT_VAR!=PotFlag_pL[upL]) continue;
            if(Sigma0_Feed==0&&BASELINE_VAR==0){
                val_flat_const = MinChi2NDF;
            }
            else if(Sigma0_Feed==0&&BASELINE_VAR==10){
                val_flat_cubic = MinChi2NDF;
            }
            else if(Sigma0_Feed==1&&BASELINE_VAR==0){
                val_chiral_const = MinChi2NDF;
            }
            else if(Sigma0_Feed==1&&BASELINE_VAR==10){
                val_chiral_cubic = MinChi2NDF;
            }
            else if(Sigma0_Feed==2&&BASELINE_VAR==0){
                val_esc16_const = MinChi2NDF;
            }
            else if(Sigma0_Feed==2&&BASELINE_VAR==10){
                val_esc16_cubic = MinChi2NDF;
            }
        }
//double nsigma = sqrt(2)*TMath::ErfcInverse(TMath::Prob(Chi2,TMath::Nint(NDF)));
        double pval;

        pval = TMath::Prob(val_chiral_cubic*NDF-BestChi2_cubic,TMath::Nint(NDF)); if(pval<1e-24)pval=1e-24; if(pval>1)pval=1;
        nsig_chiral_cubic = sqrt(2)*TMath::ErfcInverse(pval);

        pval = TMath::Prob(val_chiral_const*NDF-BestChi2_const,TMath::Nint(NDF)); if(pval<1e-24)pval=1e-24; if(pval>1)pval=1;
        nsig_chiral_const = sqrt(2)*TMath::ErfcInverse(pval);

        pval = TMath::Prob(val_esc16_cubic*NDF-BestChi2_cubic,TMath::Nint(NDF)); if(pval<1e-24)pval=1e-24; if(pval>1)pval=1;
        nsig_esc16_cubic = sqrt(2)*TMath::ErfcInverse(pval);

        pval = TMath::Prob(val_esc16_const*NDF-BestChi2_const,TMath::Nint(NDF)); if(pval<1e-24)pval=1e-24; if(pval>1)pval=1;
        nsig_esc16_const = sqrt(2)*TMath::ErfcInverse(pval);

        pval = TMath::Prob(val_flat_cubic*NDF-BestChi2_cubic,TMath::Nint(NDF)); if(pval<1e-24)pval=1e-24; if(pval>1)pval=1;
        nsig_flat_cubic = sqrt(2)*TMath::ErfcInverse(pval);

        pval = TMath::Prob(val_flat_const*NDF-BestChi2_const,TMath::Nint(NDF)); if(pval<1e-24)pval=1e-24; if(pval>1)pval=1;
        nsig_flat_const = sqrt(2)*TMath::ErfcInverse(pval);

        printf("chiral_cubic:\n");
        printf("  delta = %.2f\n",val_chiral_cubic*NDF-BestChi2_cubic);
        printf("  p_val = %.2e\n",TMath::Prob(val_chiral_cubic*NDF-BestChi2_cubic,TMath::Nint(NDF)));
        printf("  nsig = %.3f\n",nsig_chiral_cubic);
        printf("chiral_const:\n");
        printf("  delta = %.2f\n",val_chiral_const*NDF-BestChi2_const);
        printf("  p_val = %.2e\n",TMath::Prob(val_chiral_const*NDF-BestChi2_const,TMath::Nint(NDF)));
        printf("  nsig = %.3f\n",nsig_chiral_const);
        printf("flat_cubic:\n");
        printf("  delta = %.2f\n",val_flat_cubic*NDF-BestChi2_cubic);
        printf("  p_val = %.2e\n",TMath::Prob(val_flat_cubic*NDF-BestChi2_cubic,TMath::Nint(NDF)));
        printf("  nsig = %.3f\n",nsig_flat_cubic);
        printf("flat_const:\n");
        printf("  delta = %.2f\n",val_flat_const*NDF-BestChi2_const);
        printf("  p_val = %.2e\n",TMath::Prob(val_flat_const*NDF-BestChi2_const,TMath::Nint(NDF)));
        printf("  nsig = %.3f\n",nsig_flat_const);

        if(PotFlag_pL[upL]==11600) fprintf(fptr,"\\pmb{");
        fprintf(fptr,"%s",PotName_pL[upL].Data());
        if(PotFlag_pL[upL]==11600) fprintf(fptr,"}");
        if(Compact){
            fprintf(fptr," & %.1f (%.1f) & %.1f (%.1f)",
                    nsig_chiral_cubic,nsig_chiral_const,
                    nsig_flat_cubic,nsig_flat_const);
        }
        else{
            fprintf(fptr," & %.1f (%.1f) & %.1f (%.1f) & %.1f (%.1f)",
                    nsig_chiral_cubic,nsig_chiral_const,
                    nsig_esc16_cubic,nsig_esc16_const,
                    nsig_flat_cubic,nsig_flat_const);
        }

        fprintf(fptr," \\\\ \n");

    }
    fprintf(fptr,"\\hline \\hline \n");
    fprintf(fptr,"\\end{tabular} \n");
    fprintf(fptr,"\\caption{Exclusion in $n_{\\sigma}$ for the different interaction hypotheses "
            "of $\\text{p--}\\Lambda$ and $\\text{p--}\\Sigma^0$. "
            "The default values correspond to the fit with a cubic baseline, the values in brackets represent the results from using a constant baseline. "
            "The default model (in bold) is the $\\chi$EFT NLO calculation, at a cutoff of 600~MeV. } \n");
    fprintf(fptr,"\\label{tab:nsig_2020} \n");
    fprintf(fptr,"\\end{center} \n");
    fprintf(fptr,"\\end{table} \n");






    fptr=fopen(OutputFileNameSigDev.Data(),"w");
    if(Compact){
        fprintf(fptr,"\\begin{table} \n");
        fprintf(fptr,"\\begin{center} \n");
        fprintf(fptr,"\\begin{tabular}{r|c|c|} \n");
        fprintf(fptr,"\\hline \\hline \n");
        fprintf(fptr,"$\\text{p--}\\Sigma^0$ ($\\rightarrow$)& \\pmb{$\\chi$EFT} & Flat \\\\ $\\text{p--}\\Lambda$ ($\\downarrow$) & & \\\\ \n");
        fprintf(fptr,"\\hline \n");
    }
    else{
        fprintf(fptr,"\\begin{table} \n");
        fprintf(fptr,"\\begin{center} \n");
        fprintf(fptr,"\\begin{tabular}{r|c|c|c|} \n");
        fprintf(fptr,"\\hline \\hline \n");
        fprintf(fptr,"$\\text{p--}\\Sigma^0$ ($\\rightarrow$)& \\pmb{$\\chi$EFT} & ESC16 & Flat \\\\ $\\text{p--}\\Lambda$ ($\\downarrow$) & & & \\\\ \n");
        fprintf(fptr,"\\hline \n");
    }

    //char* buffer = new char [512];
    for(unsigned upL=0; upL<NumLamVars; upL++){
        float val_chiral_cubic=0;
        float val_esc16_cubic=0;
        float val_flat_cubic=0;
        float val_chiral_const=0;
        float val_esc16_const=0;
        float val_flat_const=0;

        float dsig_chiral_cubic=0;
        float dsig_esc16_cubic=0;
        float dsig_flat_cubic=0;
        float dsig_chiral_const=0;
        float dsig_esc16_const=0;
        float dsig_flat_const=0;

        for(unsigned uEntry=0; uEntry<NumEntries; uEntry++){
            InfoTree->GetEntry(uEntry);
            if(POT_VAR!=PotFlag_pL[upL]) continue;
            if(Sigma0_Feed==0&&BASELINE_VAR==0){
                val_flat_const = MinChi2NDF;
            }
            else if(Sigma0_Feed==0&&BASELINE_VAR==10){
                val_flat_cubic = MinChi2NDF;
            }
            else if(Sigma0_Feed==1&&BASELINE_VAR==0){
                val_chiral_const = MinChi2NDF;
            }
            else if(Sigma0_Feed==1&&BASELINE_VAR==10){
                val_chiral_cubic = MinChi2NDF;
            }
            else if(Sigma0_Feed==2&&BASELINE_VAR==0){
                val_esc16_const = MinChi2NDF;
            }
            else if(Sigma0_Feed==2&&BASELINE_VAR==10){
                val_esc16_cubic = MinChi2NDF;
            }
        }
//double nsigma = sqrt(2)*TMath::ErfcInverse(TMath::Prob(Chi2,TMath::Nint(NDF)));
        double pval;

        pval = TMath::Prob(val_chiral_cubic*NDF,TMath::Nint(NDF)); if(pval<1e-24)pval=1e-24; if(pval>1)pval=1;
        dsig_chiral_cubic = sqrt(2)*TMath::ErfcInverse(pval);

        pval = TMath::Prob(val_chiral_const*NDF,TMath::Nint(NDF)); if(pval<1e-24)pval=1e-24; if(pval>1)pval=1;
        dsig_chiral_const = sqrt(2)*TMath::ErfcInverse(pval);

        pval = TMath::Prob(val_esc16_cubic*NDF,TMath::Nint(NDF)); if(pval<1e-24)pval=1e-24; if(pval>1)pval=1;
        dsig_esc16_cubic = sqrt(2)*TMath::ErfcInverse(pval);

        pval = TMath::Prob(val_esc16_const*NDF,TMath::Nint(NDF)); if(pval<1e-24)pval=1e-24; if(pval>1)pval=1;
        dsig_esc16_const = sqrt(2)*TMath::ErfcInverse(pval);

        pval = TMath::Prob(val_flat_cubic*NDF,TMath::Nint(NDF)); if(pval<1e-24)pval=1e-24; if(pval>1)pval=1;
        dsig_flat_cubic = sqrt(2)*TMath::ErfcInverse(pval);

        pval = TMath::Prob(val_flat_const*NDF,TMath::Nint(NDF)); if(pval<1e-24)pval=1e-24; if(pval>1)pval=1;
        dsig_flat_const = sqrt(2)*TMath::ErfcInverse(pval);

        printf("chiral_cubic:\n");
        printf("  chi2 = %.2f\n",val_chiral_cubic);
        printf("  NDF = %.0f\n",NDF);
        printf("  p_val = %.2e\n",TMath::Prob(val_chiral_cubic*NDF,TMath::Nint(NDF)));
        printf("  nsig = %.3f\n",dsig_chiral_cubic);
        printf("chiral_const:\n");

        printf("flat_cubic:\n");

        printf("flat_const:\n");


        if(PotFlag_pL[upL]==11600) fprintf(fptr,"\\pmb{");
        fprintf(fptr,"%s",PotName_pL[upL].Data());
        if(PotFlag_pL[upL]==11600) fprintf(fptr,"}");
        if(Compact){
            fprintf(fptr," & %.1f (%.1f) & %.1f (%.1f)",
                    dsig_chiral_cubic,dsig_chiral_const,
                    dsig_flat_cubic,dsig_flat_const);
        }
        else{
            fprintf(fptr," & %.1f (%.1f) & %.1f (%.1f) & %.1f (%.1f)",
                    dsig_chiral_cubic,dsig_chiral_const,
                    dsig_esc16_cubic,dsig_esc16_const,
                    dsig_flat_cubic,dsig_flat_const);
        }

        fprintf(fptr," \\\\ \n");

    }
    fprintf(fptr,"\\hline \\hline \n");
    fprintf(fptr,"\\end{tabular} \n");
    fprintf(fptr,"\\caption{Deviation in $n_{\\sigma}$ for the different interaction hypotheses "
            "of $\\text{p--}\\Lambda$ and $\\text{p--}\\Sigma^0$. "
            "The default values correspond to the fit with a cubic baseline, the values in brackets represent the results from using a constant baseline. "
            "The default model (in bold) is the $\\chi$EFT NLO calculation, at a cutoff of 600~MeV. } \n");
    fprintf(fptr,"\\label{tab:dsig_2020} \n");
    fprintf(fptr,"\\end{center} \n");
    fprintf(fptr,"\\end{table} \n");




    delete [] PotName_pL;
    //delete [] PotName_pS0;
    delete [] PotFlag_pL;
    //delete [] PotFlag_pS0;
    delete [] InfoFileName;
    delete InfoFile;
}


//compare the NLO results for the computation including s and d waves for:
//NLO(13/19) cut off (500/550/600/650)
//make a ratio plot to compare all results to NLO19 600
void CompareChiralNLO_pLambda(){

    const TString OutputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/CompareChiralNLO_pLambda/";
    const unsigned NumTypes = 3;
    const unsigned NumCutOff = 4;
    int* CutOff = new int [NumCutOff];
    CutOff[0]=500; CutOff[1]=550; CutOff[2]=600; CutOff[3]=650;
    int* LineStyle = new int [NumCutOff];
    LineStyle[0]=2; LineStyle[1]=3; LineStyle[2]=1; LineStyle[3]=4;
    TString* TypeName = new TString [NumTypes];
    TypeName[0] = "NLO13"; TypeName[1] = "NLO19"; TypeName[2] = "LO13";
    int* Type = new int [NumTypes];
    Type[0]=0; Type[1]=1;
    int* LineColor = new int [NumTypes];
    LineColor[0] = kRed+2; LineColor[1] = kAzure+7; LineColor[2] = kGreen+1;
    const int DefType=1;
    const int DefCutOff=2;

    const double kMin=0;
    const double kMax=360;
    const unsigned NumMomBins=30;

    DLM_Histo<float>** dlmChiral = new DLM_Histo<float>* [NumTypes];
    DLM_Histo<float>** dlmChiralRatio = new DLM_Histo<float>* [NumTypes];

    TGraph** gChiral = new TGraph* [NumTypes];
    TGraph** gChiralRatio = new TGraph* [NumTypes];

    for(unsigned uType=0; uType<NumTypes; uType++){
        dlmChiral[uType] = new DLM_Histo<float> [NumCutOff];
        dlmChiralRatio[uType] = new DLM_Histo<float> [NumCutOff];
        gChiral[uType] = new TGraph [NumCutOff];
        gChiralRatio[uType] = new TGraph [NumCutOff];

        for(unsigned uCutOff=0; uCutOff<NumCutOff; uCutOff++){
            printf("\r\033[K uType=%u; uCutOff=%u",uType,uCutOff);
            cout << flush;

            dlmChiral[uType][uCutOff].SetUp(1);
            dlmChiral[uType][uCutOff].SetUp(0,NumMomBins,kMin,kMax);
            dlmChiral[uType][uCutOff].Initialize();

            gChiral[uType][uCutOff].SetName(TString::Format("gChiral_%s_%i",TypeName[uType].Data(),CutOff[uCutOff]));
            gChiral[uType][uCutOff].SetLineWidth(4);
            gChiral[uType][uCutOff].SetLineColor(LineColor[uType]);
            gChiral[uType][uCutOff].SetLineStyle(LineStyle[uCutOff]);

            gChiralRatio[uType][uCutOff].SetName(TString::Format("gChiralRatio_%s_%i",TypeName[uType].Data(),CutOff[uCutOff]));
            gChiralRatio[uType][uCutOff].SetLineWidth(4);
            gChiralRatio[uType][uCutOff].SetLineColor(LineColor[uType]);
            gChiralRatio[uType][uCutOff].SetLineStyle(LineStyle[uCutOff]);

            DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder("/mnt/Ubuntu_Data/CernBox/Sync/CatsFiles");
            CATS Kitty;
            Kitty.SetMomBins(NumMomBins,kMin,kMax);
            if(TypeName[uType].Contains("NLO")) AnalysisObject.SetUpCats_pL(Kitty,"Chiral_Coupled_SPD","Gauss",Type[uType]*10000+1*1000+CutOff[uCutOff],0);
            //else if(CutOff[uCutOff]==600) AnalysisObject.SetUpCats_pL(Kitty,"LO_Coupled_S","Gauss",0,0);
            else if(CutOff[uCutOff]==600) AnalysisObject.SetUpCats_pL(Kitty,"Chiral_Coupled_SPD","Gauss",-11600,0);
            else continue;
            Kitty.SetNotifications(CATS::nWarning);
            Kitty.KillTheCat();

            for(unsigned uBin=0; uBin<NumMomBins; uBin++){
                dlmChiral[uType][uCutOff].SetBinContent(uBin,Kitty.GetCorrFun(uBin));
            }
            dlmChiralRatio[uType][uCutOff] = dlmChiral[uType][uCutOff];
        }
    }
    printf("\r\033[K");

    for(unsigned uType=0; uType<NumTypes; uType++){
        for(unsigned uCutOff=0; uCutOff<NumCutOff; uCutOff++){
            if(CutOff[uCutOff]!=600&&!(TypeName[uType].Contains("NLO"))) continue;
            dlmChiralRatio[uType][uCutOff] /= dlmChiral[DefType][DefCutOff];
            for(unsigned uBin=0; uBin<NumMomBins; uBin++){
                gChiral[uType][uCutOff].SetPoint(uBin,
                                                 dlmChiral[uType][uCutOff].GetBinCenter(0,uBin),
                                                 dlmChiral[uType][uCutOff].GetBinContent(uBin));
                gChiralRatio[uType][uCutOff].SetPoint(  uBin,
                                                        dlmChiralRatio[uType][uCutOff].GetBinCenter(0,uBin),
                                                        dlmChiralRatio[uType][uCutOff].GetBinContent(uBin));
            }
        }
    }

    TFile OutputFile(OutputFolder+"OutputFile.root","recreate");
    for(unsigned uType=0; uType<NumTypes; uType++){
        for(unsigned uCutOff=0; uCutOff<NumCutOff; uCutOff++){
            gChiral[uType][uCutOff].Write();
            gChiralRatio[uType][uCutOff].Write();
        }
    }

    delete [] CutOff;
    delete [] LineStyle;
    delete [] Type;
    delete [] TypeName;
    delete [] LineColor;
    for(unsigned uType=0; uType<NumTypes; uType++){
        delete [] dlmChiral[uType];
    }
    delete [] dlmChiral;

}


double Fit_pLambda_Splines_Cusp(double* xVal, double* pars){
    double& Mom = *xVal;
    double Cusp = 0;
    if(-(Mom-pars[1])>=0) Cusp = pars[0]/sqrt(-(Mom-pars[1]));
    else if((Mom-pars[1])>=0) Cusp = pars[0]/sqrt((Mom-pars[1]));
    double Splines = DLM_FITTER2_FUNCTION_SPLINE3(xVal,&pars[2]);
    return Splines+Cusp;
}

void pLambda_Spline_Fit_Test(){
    DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder("/mnt/Ubuntu_Data/CernBox/Sync/CatsFiles");
    TH1F* hData_pL = AnalysisObject.GetAliceExpCorrFun("pp13TeV_HM_Dec19","pLambda","_0",0,false,-1);
    //0,60,120,180,240,320,400,500,600
    const int NumKnots = 14;
    double* Nodes_x = new double [NumKnots];
    Nodes_x[0] = 0;
    Nodes_x[1] = 50;
    Nodes_x[2] = 100;
    Nodes_x[3] = 150;
    Nodes_x[4] = 200;
    Nodes_x[5] = 230;
    Nodes_x[6] = 260;
    Nodes_x[7] = 280;
    Nodes_x[8] = 290;
    Nodes_x[9] = 300;
    Nodes_x[10] = 320;
    Nodes_x[11] = 400;
    Nodes_x[12] = 500;
    Nodes_x[13] = 600;

    TF1* fit_pL = new TF1("fit_pL",Fit_pLambda_Splines_Cusp,0,600,2+3+NumKnots*2);
    fit_pL->FixParameter(0,0.00);
    fit_pL->FixParameter(1,275);
    fit_pL->FixParameter(2,NumKnots);
    //derivative at the firs knot
    fit_pL->FixParameter(3,0);
    fit_pL->SetParameter(4,0);
    fit_pL->SetParLimits(4,-1e-3,1e-3);
    for(unsigned uKnot=0; uKnot<NumKnots; uKnot++){
        double HistVal = hData_pL->GetBinContent(hData_pL->FindBin(Nodes_x[uKnot]));
        fit_pL->FixParameter(3+2+uKnot,Nodes_x[uKnot]);
        fit_pL->SetParameter(3+2+NumKnots+uKnot,HistVal);
        fit_pL->SetParLimits(3+2+NumKnots+uKnot,0,HistVal*2);
    }
    fit_pL->SetNpx(1024);
    hData_pL->Fit(fit_pL,"S, N, R, M");
    printf("chi2/ndf = %.2f / %i\n",fit_pL->GetChisquare(),fit_pL->GetNDF());
    printf("prob = %.4f\n",fit_pL->GetProb());
    printf("nsigma = %.2f\n",sqrt(2)*TMath::ErfcInverse(fit_pL->GetProb()));

    TFile fOutput("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Unfolding/Test1/pLambda_Spline_Fit_Test.root",
                  "recreate");
    hData_pL->Write();
    fit_pL->Write();

    delete [] Nodes_x;
    delete hData_pL;
    delete fit_pL;
}



DLM_CkDecomposition* pLambda_Spline_CkFit_CKDEC;
double pLambda_Spline_CkFit(double* xVal, double* pars){
    for(unsigned uPar=0; uPar<pLambda_Spline_CkFit_CKDEC->GetCk()->GetNumSourcePar(); uPar++){
        pLambda_Spline_CkFit_CKDEC->GetCk()->SetSourcePar(uPar,pars[uPar]);
    }
    pLambda_Spline_CkFit_CKDEC->Update();
    return pLambda_Spline_CkFit_CKDEC->EvalCk(*xVal);
}
double pLambda_Spline_CkFitRAW(double* xVal, double* pars){
    return DLM_FITTER2_FUNCTION_SPLINE3(xVal,pars);
}


void pLambda_Spline_Fit_Test2(){
    DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder("/mnt/Ubuntu_Data/CernBox/Sync/CatsFiles");
    TH1F* hData_pL = AnalysisObject.GetAliceExpCorrFun("pp13TeV_HM_Dec19","pLambda","_0",0,false,-1);
    //0,60,120,180,240,320,400,500,600
    const int NumKnots = 14;
    double* Nodes_x = new double [NumKnots];
    Nodes_x[0] = 0;
    Nodes_x[1] = 50;
    Nodes_x[2] = 100;
    Nodes_x[3] = 150;
    Nodes_x[4] = 200;
    Nodes_x[5] = 230;
    Nodes_x[6] = 260;
    Nodes_x[7] = 280;
    Nodes_x[8] = 290;
    Nodes_x[9] = 300;
    Nodes_x[10] = 320;
    Nodes_x[11] = 400;
    Nodes_x[12] = 500;
    Nodes_x[13] = 600;

    TF1* fit_pL = new TF1("fit_pL",pLambda_Spline_CkFitRAW,0,600,3+NumKnots*2);
    fit_pL->FixParameter(0,NumKnots);
    //derivative at the firs knot
    fit_pL->FixParameter(1,0);
    fit_pL->SetParameter(2,0);
    fit_pL->SetParLimits(2,-1e-3,1e-3);
    for(unsigned uKnot=0; uKnot<NumKnots; uKnot++){
        double HistVal = hData_pL->GetBinContent(hData_pL->FindBin(Nodes_x[uKnot]));
        fit_pL->FixParameter(3+uKnot,Nodes_x[uKnot]);
        fit_pL->SetParameter(3+NumKnots+uKnot,HistVal);
        fit_pL->SetParLimits(3+NumKnots+uKnot,0,HistVal*2);
    }
    fit_pL->SetNpx(1024);
    hData_pL->Fit(fit_pL,"S, N, R, M");
    printf("chi2/ndf = %.2f / %i\n",fit_pL->GetChisquare(),fit_pL->GetNDF());
    printf("prob = %.4f\n",fit_pL->GetProb());
    printf("nsigma = %.2f\n",sqrt(2)*TMath::ErfcInverse(fit_pL->GetProb()));

    TFile fOutput("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Unfolding/Test1/pLambda_Spline_Fit_Test2.root",
                  "recreate");
    hData_pL->Write();
    fit_pL->Write();

    delete [] Nodes_x;
    delete hData_pL;
    delete fit_pL;
}


double* pLambda_Spline_Ck_ARRAY_SPLINE3_X=NULL;
double* pLambda_Spline_Ck_ARRAY_SPLINE3_Y=NULL;
double pLambda_Spline_Ck(const double& Momentum, const double* SourcePar, const double* PotPar){
    //[0] = NumKnots
    //[1] = der at 0
    //[2] = der at last
    //[3]... posX
    //[...]... poxY
    const int MAX_KNOTS = 20;
    int NumKnots = TMath::Nint(SourcePar[0]);
    if(NumKnots<2) NumKnots=2;
    if(NumKnots>MAX_KNOTS) NumKnots=MAX_KNOTS;
    if(!pLambda_Spline_Ck_ARRAY_SPLINE3_X) pLambda_Spline_Ck_ARRAY_SPLINE3_X = new double [MAX_KNOTS];
    if(!pLambda_Spline_Ck_ARRAY_SPLINE3_Y) pLambda_Spline_Ck_ARRAY_SPLINE3_Y = new double [MAX_KNOTS];
    for(int iKnot=0; iKnot<NumKnots; iKnot++){
        pLambda_Spline_Ck_ARRAY_SPLINE3_X[iKnot] = SourcePar[3+iKnot];
        pLambda_Spline_Ck_ARRAY_SPLINE3_Y[iKnot] = SourcePar[3+NumKnots+iKnot];
        //fix to the previous one of the value is fixed to 1e6
        if(pLambda_Spline_Ck_ARRAY_SPLINE3_Y[iKnot]==1e6&&iKnot) pLambda_Spline_Ck_ARRAY_SPLINE3_Y[iKnot] = SourcePar[3+NumKnots+iKnot-1];
    }
    const double& derStart = SourcePar[1];
    const double& derEnd = SourcePar[2];
    TSpline3 sp3("sp3", pLambda_Spline_Ck_ARRAY_SPLINE3_X, pLambda_Spline_Ck_ARRAY_SPLINE3_Y, NumKnots, "b1e1", derStart, derEnd);
    return sp3.Eval(Momentum);
}
void pLambda_Spline_Fit_Unfold1(){

    const double kMin=0;
    const double kMax=600;
    const unsigned NumMomBins=75;

    //0,60,120,180,240,320,400,500,600
    const int NumKnots = 14;
    double* Nodes_x = new double [NumKnots];
    Nodes_x[0] = kMin;
    Nodes_x[1] = 50;
    Nodes_x[2] = 100;
    Nodes_x[3] = 150;
    Nodes_x[4] = 200;
    Nodes_x[5] = 240;
    Nodes_x[6] = 270;
    Nodes_x[7] = 280;
    Nodes_x[8] = 290;
    Nodes_x[9] = 300;
    Nodes_x[10] = 350;
    Nodes_x[11] = 400;
    Nodes_x[12] = 500;
    Nodes_x[13] = kMax;

    DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder("/mnt/Ubuntu_Data/CernBox/Sync/CatsFiles");
    TH1F* hData_pL = AnalysisObject.GetAliceExpCorrFun("pp13TeV_HM_Dec19","pLambda","_0",1,false,-1);

    TF1* fit_pL_RAW = new TF1("fit_pL_RAW",pLambda_Spline_CkFitRAW,kMin,kMax,3+NumKnots*2);
    fit_pL_RAW->FixParameter(0,NumKnots);
    //derivative at the firs knot
    fit_pL_RAW->FixParameter(1,0);
    fit_pL_RAW->SetParameter(2,0);
    fit_pL_RAW->SetParLimits(2,-1e-3,1e-3);
    for(unsigned uKnot=0; uKnot<NumKnots; uKnot++){
        double HistVal = hData_pL->GetBinContent(hData_pL->FindBin(Nodes_x[uKnot]));
        fit_pL_RAW->FixParameter(3+uKnot,Nodes_x[uKnot]);
        fit_pL_RAW->SetParameter(3+NumKnots+uKnot,HistVal);
        fit_pL_RAW->SetParLimits(3+NumKnots+uKnot,0,HistVal*2);
    }
    fit_pL_RAW->SetNpx(1024);
    hData_pL->Fit(fit_pL_RAW,"S, N, R, M");

    DLM_Ck* Ck_pL = new DLM_Ck(3+NumKnots*2,0,NumMomBins,kMin,kMax,pLambda_Spline_Ck);
    //Ck_pL->SetSourcePar(0,ResidualSourceSize);
    TH2F* hResolution_pL = AnalysisObject.GetResolutionMatrix("pp13TeV_HM_Dec19","pLambda");
    DLM_CkDecomposition* CkDec_pL = new DLM_CkDecomposition("pLambda",0,*Ck_pL,hResolution_pL);
    pLambda_Spline_CkFit_CKDEC = CkDec_pL;

    TF1* fit_pL = new TF1("fit_pL",pLambda_Spline_CkFit,kMin,kMax,3+NumKnots*2);
    fit_pL->FixParameter(0,NumKnots);
    //derivative at the firs knot
    fit_pL->FixParameter(1,0);
    fit_pL->SetParameter(2,0);
    fit_pL->SetParLimits(2,-1e-3,1e-3);
    for(unsigned uKnot=0; uKnot<NumKnots; uKnot++){
        double HistVal = hData_pL->GetBinContent(hData_pL->FindBin(Nodes_x[uKnot]));
        fit_pL->FixParameter(3+uKnot,Nodes_x[uKnot]);
        fit_pL->SetParameter(3+NumKnots+uKnot,HistVal);
        fit_pL->SetParLimits(3+NumKnots+uKnot,0,HistVal*2);
    }
    //fit_pL->SetNpx(1024);

    for(unsigned uPar=0; uPar<CkDec_pL->GetCk()->GetNumSourcePar(); uPar++){
        CkDec_pL->GetCk()->SetSourcePar(uPar,fit_pL->GetParameter(uPar));
    }
    CkDec_pL->Update();

    hData_pL->Fit(fit_pL,"S, N, R, M");
    printf("chi2/ndf = %.2f / %i\n",fit_pL->GetChisquare(),fit_pL->GetNDF());
    printf("prob = %.4f\n",fit_pL->GetProb());
    printf("nsigma = %.2f\n",sqrt(2)*TMath::ErfcInverse(fit_pL->GetProb()));

    TRandom3 rangen(11);
    //TH1F* hData_Unfolded = (TH1F*)hData_pL->Clone("hData_Unfolded");
    DLM_Ck* Ck_Unfolded = new DLM_Ck(0,0,NumMomBins,kMin,kMax,NULL);
    DLM_CkDecomposition* CkDec_Unfolded = new DLM_CkDecomposition("pLambda",0,*Ck_Unfolded,hResolution_pL);
    const unsigned NumIterPerStep = 1024;
    const unsigned MaxStepsWithoutImprovement = 64;
    const unsigned MaxTotalIter = 1024*8;
    unsigned TotalIter = 0;
    unsigned StepsWithoutImprovement=0;
    const unsigned NumSteps = 9;
    const double Chi2Goal = 10;
    DLM_Ck* Ck_Best = new DLM_Ck(0,0,NumMomBins,kMin,kMax,NULL);
    DLM_Ck* Ck_Sample = new DLM_Ck(0,0,NumMomBins,kMin,kMax,NULL);
    double BestChi2 = 1e6;
    unsigned NumGoodSolutions = 0;
    for(unsigned uStep=0; uStep<NumSteps; uStep++){
        double Scale;
        double Momentum;
        //double DataValue;
        double OriginalValue;
        double RandomValue;
        double Error;
        double CkVal;
        double CkValCorrected;
        double Chi2;
        double Chi2_Bin;
        double Chi2_BestFit;
        //i.e. the smeared correlation lies within the exp. data with uncertainties
        //we demand that >68% of the bins lie within 1sigma, and there are no bins lying outside of 2sigma
        bool GoodSolution;
        int NumBins1Sigma=0;
        int NumBins2Sigma=0;
        int NumBinsAnySigma=0;
        int NumBins1Sigma_Th=0;
        int NumBins2Sigma_Th=0;
        int NumBinsAnySigma_Th=0;
        //StepsWithoutImprovement=0;
        if(BestChi2==1e6){
            Ck_Sample->Copy(*Ck_pL);
            //bootstrap the sampling histogram ones, to avoid a bias in the chi2 determination (where the original solution is used as a benchmark)
            for(unsigned uBin=0; uBin<Ck_Sample->GetNbins(); uBin++){
                OriginalValue = Ck_Sample->GetBinContent(uBin);
                Error = hData_pL->GetBinError(uBin+1);
                //RandomValue = rangen.Uniform(OriginalValue-0.5*Error,OriginalValue+0.5*Error);
                RandomValue = rangen.Gaus(OriginalValue,Error);
                Ck_Sample->SetBinContent(uBin,RandomValue);
            }
        }
        else Ck_Sample->Copy(*Ck_Best);
        for(unsigned uIter=0; uIter<NumIterPerStep; uIter++){
            if(!uStep&&!uIter) Scale=0;
            else Scale = pow(0.5,double(uStep));
            bool ChangeNextFewBins=false;
            unsigned ChangeBinsUpTo;
            for(unsigned uBin=0; uBin<Ck_Unfolded->GetNbins(); uBin++){
                Momentum = Ck_Unfolded->GetBinCenter(0,uBin);
                //if(BestChi2>1e5) OriginalValue = fit_pL->Eval(Momentum);
                //else OriginalValue = Ck_Sample->GetBinContent(uBin);
                OriginalValue = Ck_Sample->GetBinContent(uBin);
                Error = hData_pL->GetBinError(uBin+1);
                if(Momentum<kMin) continue;
                if(Momentum>kMax) break;

                //change only some bins (clusters of 5). We have 20% prob to change the next five bins
                if(!ChangeNextFewBins&&rangen.Uniform(1)<0.2){
                    ChangeNextFewBins=true;
                    ChangeBinsUpTo=uBin+5;
                }
                else if(uBin>=ChangeBinsUpTo){
                    ChangeNextFewBins=false;
                }

                //printf("-> ERR=%.4f; SCA=%.3f\n",Error,Scale);
                if(!ChangeNextFewBins) RandomValue=OriginalValue;
                //half of the time we fine tune
                else if(uIter%2) RandomValue = rangen.Uniform(OriginalValue-0.5*Error*Scale,OriginalValue+0.5*Error*Scale);
                //half of the time we allow for a larger spread
                else RandomValue = rangen.Gaus(OriginalValue,Error*Scale);
                //if(uBin==0)
                //printf(" OV = %.4f; RV = %.4f (ERR=%.4f; SCA=%.3f)\n",OriginalValue,RandomValue,Error,Scale);
                Ck_Unfolded->SetBinContent(uBin,RandomValue);
            }
            Ck_Unfolded->Update(true);
            CkDec_Unfolded->Update(true);

//printf("BIN0: %f -> %f\n",Ck_Unfolded->GetBinContent(unsigned(0)),CkDec_Unfolded->EvalCk(Ck_Unfolded->GetBinCenter(0,0)));
            Chi2=0;
            GoodSolution=true;
            for(unsigned uBin=0; uBin<hData_pL->GetNbinsX(); uBin++){
                Momentum = hData_pL->GetBinCenter(uBin+1);
                if(Momentum<kMin) continue;
                if(Momentum>kMax) break;
                CkValCorrected = CkDec_Unfolded->EvalCk(Momentum);
                CkVal = hData_pL->GetBinContent(uBin+1);
                Error = hData_pL->GetBinError(uBin+1);
                //this is the chi2 between the data and the corrected solution (smeared)
                Chi2_Bin = pow((CkVal-CkValCorrected)/(Error),2.);
                Chi2 += Chi2_Bin;
                //the Chi2 (estimated based on the uncertainty of the data) between the theoretical FIT curve
                //and the current solution (theoretical Ck). We want to have good agreement here as well,
                //to ensure a smooth corrected solution.
                Chi2_BestFit = pow((Ck_Unfolded->GetBinContent(uBin)-Ck_pL->GetBinContent(uBin))/(Error),2.);
                Chi2 += Chi2_BestFit;
                NumBins1Sigma+=(Chi2_Bin<=1);
                NumBins2Sigma+=(Chi2_Bin<=4);
                NumBinsAnySigma++;
                NumBins1Sigma_Th+=(Chi2_BestFit<=1);
                NumBins2Sigma_Th+=(Chi2_BestFit<=4);
                NumBinsAnySigma_Th++;
            }
            if(double(NumBins1Sigma)/double(NumBinsAnySigma)<0.68||NumBins2Sigma!=NumBinsAnySigma) GoodSolution=false;
            if(double(NumBins1Sigma_Th)/double(NumBinsAnySigma_Th)<0.68||NumBins2Sigma_Th!=NumBinsAnySigma_Th) GoodSolution=false;
            NumGoodSolutions+=GoodSolution;
            if(Chi2<BestChi2){
                BestChi2=Chi2;
                Ck_Best[0].Copy(Ck_Unfolded[0]);
                Ck_Sample[0].Copy(Ck_Unfolded[0]);
                StepsWithoutImprovement=0;
            }
            else if(Chi2<BestChi2*1.1&&StepsWithoutImprovement>64){
                Ck_Sample[0].Copy(Ck_Unfolded[0]);
                printf("CHANGED SAMPLING\n");
                StepsWithoutImprovement++;
            }
            else if(uStep&&StepsWithoutImprovement>480&&StepsWithoutImprovement<1536&&uIter>480){
                uStep--;
                uIter=0;
                StepsWithoutImprovement++;
            }
            else if(StepsWithoutImprovement==1024){
                Ck_Sample[0].Copy(Ck_Best[0]);
                StepsWithoutImprovement++;
            }
            else{
                StepsWithoutImprovement++;
                /*
                //change randomly 20% of all bins, withing the current uncertainty
                if(StepsWithoutImprovement==MaxStepsWithoutImprovement){
                    for(unsigned uBin=0; uBin<Ck_Unfolded->GetNbins(); uBin++){
                        if(rangen.Uniform(1)>0.2) continue;
                        OriginalValue = Ck_pL->GetBinContent(uBin);
                        RandomValue = rangen.Uniform(OriginalValue-Error*Scale,OriginalValue+Error*Scale);
                        Ck_Sample->SetBinContent(uBin,RandomValue);
                    }
                    StepsWithoutImprovement=0;
                }
                */
            }
            printf(" us=%u; ui=%i; BestChi2=%.6f(%.6f) (TI=%u, GS=%u)\n",uStep,uIter,BestChi2,Chi2,TotalIter,NumGoodSolutions);
            TotalIter++;
            if(TotalIter>=MaxTotalIter) break;
            if(BestChi2<=Chi2Goal) break;
        }
    }

    Ck_Unfolded[0].Copy(Ck_Best[0]);
    Ck_Unfolded->Update(true);
    CkDec_Unfolded->Update(true);


    TFile fOutput("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Unfolding/Test1/pLambda_Spline_Fit_Unfold1.root",
                  "recreate");

    TGraph CkTheory;
    CkTheory.SetName("CkTheory");
    TGraph CkFit;
    CkFit.SetName("CkFit");
    TGraph Ck_Theory_Best;
    Ck_Theory_Best.SetName("Ck_Theory_Best");
    TGraph CkFit_Best;
    CkFit_Best.SetName("CkFit_Best");
    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        double MOM = Ck_pL->GetBinCenter(0,uBin);
        CkTheory.SetPoint(uBin,MOM,Ck_pL->GetBinContent(uBin));
        CkFit.SetPoint(uBin,MOM,CkDec_pL->EvalCk(MOM));
        Ck_Theory_Best.SetPoint(uBin,MOM,Ck_Unfolded->GetBinContent(uBin));
        CkFit_Best.SetPoint(uBin,MOM,CkDec_Unfolded->EvalCk(MOM));
    }
//printf("The hard part is done\n");
//usleep(5e6);
    hData_pL->Write();
    //hData_Unfolded->Write();
    fit_pL_RAW->Write();
    fit_pL->Write();
    CkTheory.Write();
    CkFit.Write();
    Ck_Theory_Best.Write();
    CkFit_Best.Write();

printf("Deleting\n");
    delete [] Nodes_x;
printf(" 1\n");
    delete Ck_pL;
printf(" 2\n");
    delete CkDec_pL;
printf(" 3\n");
    delete Ck_Unfolded;
printf(" 4\n");
    delete CkDec_Unfolded;
printf(" 5\n");
    delete Ck_Best;
    delete Ck_Sample;
printf(" 6\n");
    delete hData_pL;
printf(" 7\n");
    //delete hData_Unfolded;
    delete fit_pL;
printf(" 8\n");
    delete fit_pL_RAW;
printf(" 9\n");
}

//the story so far:
//if I use pp SE matrix, pp ME matrix or pLambda matrix I get essentially the same result
//if I use the original coarse pL (10MeV) SE matrix, I get deviation below 30 MeV (the unfolded Ck sits higher)
//and around the cusp I get some fluctuating bins.
//Conclusion: the coarse binning and low statistics in the old matrix was introducing a bit of bias... use the ME pL matrix instead!

//the idea of this function:
//model the theory Ck with spline. Fit the correlation (including smearing) to find out a somewhat smooth initial theory Ck that fits the problem
//than iterate over all bins, where you take bin N and N+X bins and bootstrap the theory solution randomly, transforming it back to the smeared data.
//if the description is better => accept the solution. A problem was that the solution often wandered off very much away from the original fit,
//leading to weird spikes. To avoid this, and additional constraint is that within each block of fitted bins, the maximum deviation (bin-by-bin)
//between the new smeared function and the original function should be larger than the maximum deviation between the original theory fit and the
//bootstraped solution. This way, the bootstraped solution is forced to remain close to the original fit value.
//if the SEED is zero : NO bootstrap
//else : bootstrap the data
//we use all seeds from SEEDmin to < SEEDmin+NumIter
//for a new computation, make sure the output file is deleted by hand (otherwise it just accumulates statistics)
void pLambda_Spline_Fit_Unfold2(const unsigned& SEEDmin, const unsigned& NumIter, const unsigned& TimeLimit,
                                const double& Perfect_chi2ndf, const double& VeryGood_chi2ndf, const double& Unacceptable_chi2ndf,
                                const double& BinWidth, const TString& DataVariation,
                                const char* CatsFileFolder, const TString& OutputFolder){

    //PUT FALST FOR THE BATCH FARM
    const bool PC_OUTPUT = false;

    //if true, the output is saved in a separate file for each seed
    //the 1D histograms are not filled, however if we run with
    ///NumIter==0 => we will only compute the 1D histos without adding anything new
    //bool SafeFileOpening = true;
    const double kMin=0;
    double kMax=540;
    //const unsigned NumMomBins=45;//12MeV
    const unsigned NumMomBins=TMath::Nint(kMax/BinWidth);//4MeV
    kMax = double(NumMomBins)*BinWidth;
    const bool PrePreFit = false;
    const long long JobTimeLimit = TimeLimit*60;
    bool TerminateASAP = false;
    DLM_Timer JobTime;

    //if you reach such a chi2, ignore all other constraints and terminate
    //const double Perfect_chi2ndf = 0.1;
    //similar as the perfect one, however this condition is only verified at the end of each polishing run
    //const double VeryGood_chi2ndf = 0.2;
    //any solution above that limit will be rejected
    //const double Unacceptable_chi2ndf = 0.6;

    //if it takes too long, we will settle for this one
    //const double Okayish_chi2ndf = 0.4;
    double Avg_chi2ndf=0;
    //the time (in seconds) beyond which we settle for Okayish_chi2ndf;
    long long PatienceForPerfection;
    //whatever we have after that much time, we take
    long long PatienceMax;
    const unsigned NumBackAndForth = 8;//this is a maximum value
    const unsigned MaxStepsWithoutImprovement1 = 24;//32
    const unsigned MaxStepsWithoutImprovement2 = 12;//16
    const unsigned ErrorDepth = 3;//3
    const unsigned BinDepth = TMath::Nint(60./BinWidth);

    char* strtime = new char [128];
    char* ctext1 = new char [128];
    char* ctext2 = new char [128];
    char* ctext3 = new char [128];

    const int NumKnots = 20;
    double* Nodes_x = new double [NumKnots];
    Nodes_x[0] = kMin;
    Nodes_x[1] = 24;
    Nodes_x[2] = 48;
    Nodes_x[3] = 72;
    Nodes_x[4] = 96;
    Nodes_x[5] = 120;
    Nodes_x[6] = 144;
    Nodes_x[7] = 168;
    Nodes_x[8] = 192;
    Nodes_x[9] = 216;
    Nodes_x[10] = 240;
    Nodes_x[11] = 266;
    Nodes_x[12] = 278;
    Nodes_x[13] = 290;
    Nodes_x[14] = 302;
    //Nodes_x[15] = 326;
    Nodes_x[15] = 350;
    Nodes_x[16] = 400;
    Nodes_x[17] = 450;
    Nodes_x[18] = 500;
    Nodes_x[19] = kMax;

    DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder(CatsFileFolder);
    //TH1F* hData_pL = AnalysisObject.GetAliceExpCorrFun("pp13TeV_HM_Dec19","pLambda","_0",TMath::Nint(BinWidth/4)-1,false,-1);
    TH1F* hData_pL = AnalysisObject.GetAliceExpCorrFun("pp13TeV_HM_DimiJun20","pLambda",DataVariation.Data(),TMath::Nint(BinWidth/4)-1,false,-1);
    TString OriginalName = hData_pL->GetName();

    TF1* fit_pL_RAW = new TF1("fit_pL_RAW",pLambda_Spline_CkFitRAW,kMin,kMax,3+NumKnots*2);
    fit_pL_RAW->FixParameter(0,NumKnots);
    //derivative at the firs knot
    fit_pL_RAW->FixParameter(1,0);
    fit_pL_RAW->SetParameter(2,0);
    fit_pL_RAW->SetParLimits(2,-1e-3,1e-3);
    for(unsigned uKnot=0; uKnot<NumKnots; uKnot++){
        double HistVal = hData_pL->GetBinContent(hData_pL->FindBin(Nodes_x[uKnot]));
        fit_pL_RAW->FixParameter(3+uKnot,Nodes_x[uKnot]);
        fit_pL_RAW->SetParameter(3+NumKnots+uKnot,HistVal);
        fit_pL_RAW->SetParLimits(3+NumKnots+uKnot,0,HistVal*2);
    }
    fit_pL_RAW->SetNpx(1024);
    //printf("About to fit hData_pL\n");
    hData_pL->Fit(fit_pL_RAW,"Q, S, N, R, M");
    //printf("hData_pL was fitted!\n");

    double f_chi2 = fit_pL_RAW->GetChisquare();
    int f_ndf = fit_pL_RAW->GetNumberFitPoints();
    double f_chi2ndf = f_chi2/double(f_ndf);
    double f_pval = TMath::Prob(f_chi2,f_ndf);
    double f_ns = sqrt(2)*TMath::ErfcInverse(f_pval+1e-32);
    printf("Pre-fit (folded):\n");
    printf(" chi2/ndf = %.1f/%i\n",f_chi2,f_ndf);
    //printf(" pval = %.4f\n",f_pval);
    //printf(" nsigma = %.2f\n",f_ns);


    DLM_Ck* Ck_pL = new DLM_Ck(3+NumKnots*2,0,NumMomBins,kMin,kMax,pLambda_Spline_Ck);
    //Ck_pL->SetSourcePar(0,ResidualSourceSize);
    //TH2F* hResolution_pL = AnalysisObject.GetResolutionMatrix("pp13TeV_HM_Dec19","pLambda");
    TH2F* hResolution_pL = AnalysisObject.GetResolutionMatrix("pp13TeV_HM_DimiJun20","pLambda");
//TH2F* hResolution_pL = AnalysisObject.GetResolutionMatrix("pp13TeV_HM_Dec19","pp");
//TH2F* hResolution_pL = AnalysisObject.GetResolutionMatrix("pp13TeV_HM_DimiJun20","pp");
    DLM_CkDecomposition* CkDec_pL = new DLM_CkDecomposition("pLambda",0,*Ck_pL,hResolution_pL);
    pLambda_Spline_CkFit_CKDEC = CkDec_pL;

    TF1* fit_pL = new TF1("fit_pL",pLambda_Spline_CkFit,kMin,kMax,3+NumKnots*2);
    fit_pL->FixParameter(0,NumKnots);
    //derivative at the firs knot
    fit_pL->FixParameter(1,0);
    fit_pL->SetParameter(2,0);
    fit_pL->SetParLimits(2,-1e-3,1e-3);
    for(unsigned uKnot=0; uKnot<NumKnots; uKnot++){
        double HistVal = hData_pL->GetBinContent(hData_pL->FindBin(Nodes_x[uKnot]));
        fit_pL->FixParameter(3+uKnot,Nodes_x[uKnot]);
        fit_pL->SetParameter(3+NumKnots+uKnot,HistVal);
        fit_pL->SetParLimits(3+NumKnots+uKnot,0,HistVal*2);
    }
    //fit_pL->SetNpx(1024);

    for(unsigned uPar=0; uPar<CkDec_pL->GetCk()->GetNumSourcePar(); uPar++){
        CkDec_pL->GetCk()->SetSourcePar(uPar,fit_pL->GetParameter(uPar));
    }
    CkDec_pL->Update();

    //printf("About to fit hData_pL\n");
    hData_pL->Fit(fit_pL,"Q, S, N, R, M");
    f_chi2 = fit_pL->GetChisquare();
    f_ndf = fit_pL->GetNumberFitPoints();
    f_chi2ndf = f_chi2/double(f_ndf);
    f_pval = TMath::Prob(f_chi2,f_ndf);
    f_ns = sqrt(2)*TMath::ErfcInverse(f_pval+1e-32);
    printf("Pre-fit (unfolded):\n");
    printf(" chi2/ndf = %.1f/%i\n",f_chi2,f_ndf);
    //printf(" pval = %.4f\n",f_pval);
    //printf(" nsigma = %.2f\n",f_ns);

    //2 on the first attempt, 4 when we retry
    bool AcceptSolution = (f_chi2ndf<3);
    if(!AcceptSolution){
        delete Ck_pL;
        delete CkDec_pL;
        delete fit_pL;
        delete [] Nodes_x;
        printf("REJECTED (pre-fit)\n\n");
        return;
    }

    printf("\nProceeding to the bootstraping and polishing:\n\n");

    double IterTime = 0;

    for(unsigned uSeed=SEEDmin; uSeed<SEEDmin+NumIter; uSeed++){
        bool TerminateIter = false;
        DLM_Timer ComputingTime;
        printf(" uSeed = %u (%u)\n",uSeed,SEEDmin+NumIter);
        TH1F* Chi2Ndf_Unfold = new TH1F("Chi2Ndf_Unfold","Chi2Ndf_Unfold",128,0,Unacceptable_chi2ndf);
        TH1F* CPU_Unfold = new TH1F("CPU_Unfold","CPU_Unfold",128,0,4.*double(JobTimeLimit)/double(NumIter)/60.);
        TH2F* Chi2Ndf_CPU_Unfold = new TH2F("Chi2Ndf_CPU_Unfold","Chi2Ndf_CPU_Unfold",
                                            128,0,Unacceptable_chi2ndf,
                                            128,0,4.*double(JobTimeLimit)/double(NumIter)/60.);
        long long TotalTime = JobTime.Stop()/1000000.;
        double TimeRemaining = double(JobTimeLimit-TotalTime);
        double IterRemaining = double(SEEDmin+NumIter-uSeed);
        double AvgTimePerSeedRemaining = TimeRemaining/IterRemaining;
        double NumIterDone = double(uSeed-SEEDmin);//completed iterations (not counting the current)
        double NumIterRemain = double(SEEDmin+NumIter-uSeed-1);//remaining iteration (not counting the current)
        double AvgIterTime = NumIterDone?IterTime/NumIterDone:0;
        //a bonus that we add to the remaining average time, based on the difference between the
        //remaining avg. and the current average time. If we are bad on time, this bonus becomes negative
        //and serves as a regulator to bring us back on schedule
        double BufferTime = AvgTimePerSeedRemaining-AvgIterTime;
        //the maximum allowed time is: on the first iter, the average
        //afterwards, it represents the remaining time, minus the expected time needed for the remaining iterations,
        //based on the average run time so far. This is in case we do good on time
        double MaxTime = NumIterDone?TimeRemaining-AvgIterTime*NumIterRemain:AvgTimePerSeedRemaining;
        //in case we are doing bad on average, we try to give some extra bonus on the current runs,
        //by increasing by up to 33% the allowed run time
        double BonusFactor = Avg_chi2ndf/VeryGood_chi2ndf;
        if(BonusFactor<1) BonusFactor=1;
        if(BonusFactor>1.33) BonusFactor=1.33;
        //if we do bad on time, the upper definition can lead to very very small max time.
        //if too small, we set the max time to its 'minimal' max time, which is the average time remaining per iter
        if(MaxTime<AvgTimePerSeedRemaining*BonusFactor) MaxTime = AvgTimePerSeedRemaining*BonusFactor;

        //double MaxTime = NumIterDone?TimeRemaining/(NumIterDone+1.):AvgTimePerSeedRemaining;
        double AllowedTime = AvgTimePerSeedRemaining+BufferTime;
        AllowedTime *= BonusFactor;
        if(AllowedTime>MaxTime) AllowedTime=MaxTime;
        PatienceMax = (long long)(AllowedTime);
        PatienceForPerfection = (long long)(AvgTimePerSeedRemaining);
        ShowTime(PatienceMax,strtime,2,true,5);
        printf("  Time limit: %s\n",strtime);

        //ShowTime((long long)(TimeRemaining),strtime,2,true,5);
        //printf("    TR: %s\n",strtime);
        //ShowTime((long long)(AvgTimePerSeedRemaining),strtime,2,true,5);
        //printf("  AITR: %s\n",strtime);
        //ShowTime((long long)(AvgIterTime),strtime,2,true,5);
        //printf("   AIT: %s\n",strtime);
        //ShowTime((long long)(BufferTime),strtime,2,true,5);
        //printf("  BUFF: %s\n",strtime);
        //ShowTime((long long)(MaxTime),strtime,2,true,5);
        //printf("   MAX: %s\n",strtime);

        TRandom3 rangen(uSeed+1);//avoid zero
        double OriginalValue;
        double OriginalValueMin;
        double OriginalValueMax;
        double Error;
        double RandomValue;
        double RandomMean;
        double MOM;

        TH1F *hData_boot = (TH1F*)hData_pL->Clone("hData_boot");
        for(unsigned uBin=0; uBin<NumMomBins; uBin++){
            MOM = hData_pL->GetBinCenter(uBin+1);
            OriginalValueMin = hData_pL->GetBinContent(uBin+1);
            OriginalValueMax = fit_pL_RAW->Eval(MOM);
            if(OriginalValueMax<OriginalValueMin){
                OriginalValue = OriginalValueMin;
                OriginalValueMin = OriginalValueMax;
                OriginalValueMax = OriginalValue;
            }
            Error = hData_pL->GetBinError(uBin+1);
            //if(uSeed) RandomValue = rangen.Gaus(OriginalValue,Error);
            //else RandomValue = OriginalValue;
            //this is not so trivial. If I bootstrap from the data, I get more difficult to fit solutions (kind of normal,
            //as we add fluctuations on top of the already fluctuating data points), if we sample from the best fit things look better,
            //apart from just near the cusp, where some points look to be fluctuating away from the data by too much.
            //THE SOLUTION: pre-bootstrap from a uniform distribution between the data and the fit solution, which is than used as
            //the mean value for the actual bootstrap
            if(uSeed){
                RandomMean = rangen.Uniform(OriginalValueMin,OriginalValueMax);
                RandomValue = rangen.Gaus(RandomMean,Error);
            }
            else RandomValue = hData_pL->GetBinContent(uBin+1);
            hData_boot->SetBinContent(uBin+1,RandomValue);
            hData_boot->SetBinError(uBin+1,Error);
        }

        //TH1F* hData_Unfolded = (TH1F*)hData_boot->Clone("hData_Unfolded");
        DLM_Ck* Ck_Unfolded = new DLM_Ck(0,0,NumMomBins,kMin,kMax,NULL);
        DLM_CkDecomposition* CkDec_Unfolded = new DLM_CkDecomposition("pLambda",0,*Ck_Unfolded,hResolution_pL);

        double BestChi2=1e6;
        double BestChi2_Data=1e6;
        double BestChi2_Fit=1e6;
        double BestChi2_BinData=1e6;
        double BestNDF=0;
        DLM_Ck* Ck_Best = new DLM_Ck(0,0,NumMomBins,kMin,kMax,NULL);
        DLM_Ck* Ck_Sample = new DLM_Ck(0,0,NumMomBins,kMin,kMax,NULL);
        DLM_CkDecomposition* CkDec_Best = new DLM_CkDecomposition("pLambda",0,*Ck_Best,hResolution_pL);
    //64.5, 64.6 for the fast round
    // New best fit: Chi2=85.568 (48.694+36.875) at uED=3 for the pLambda_Spline_Fit_Unfold2_TEST1_Chi2All.root


        unsigned StepsWithoutImprovement=0;



        //bootstrap the sampling histogram ones, to avoid a bias in the chi2 determination (where the original solution is used as a benchmark)
        for(unsigned uBin=0; uBin<Ck_Sample->GetNbins(); uBin++){
            OriginalValue = Ck_pL->GetBinContent(uBin);
            //Error = hData_boot->GetBinError(uBin+1);
            //RandomValue = rangen.Gaus(OriginalValue,Error);
//what the hell is this??
    RandomValue=OriginalValue;
    ///
            Ck_Sample->SetBinContent(uBin,RandomValue);
        }

        for(unsigned uBF=0; uBF<NumBackAndForth; uBF++){
            sprintf(ctext1,"   uBF = %u (%u)",uBF,NumBackAndForth);
            strcpy(ctext2,"");
            strcpy(ctext3,"");
            cout << ctext1 << flush;
            bool REFRESHED_INFO = false;
            for(unsigned uBin=0; uBin<2*Ck_Sample->GetNbins(); uBin++){
                //we iterate twice, from below and from above.
                //The second iteration is for 'polishing',
                //which is performed for fewer iterations, but higher depth
                bool Polishing = uBin>=Ck_Sample->GetNbins();
                //keep track of time, in case we exceed our limit, terminate asap
                TotalTime = JobTime.Stop()/1000000.;
                if(TotalTime>JobTimeLimit){
                    TerminateASAP = true;
                    break;
                }

                unsigned WhichBin = (!Polishing)?uBin:2*Ck_Sample->GetNbins()-uBin-1;
                if(!Polishing&&WhichBin>Ck_Sample->GetNbins()-BinDepth) continue;
                if(Polishing&&WhichBin<BinDepth-1) continue;
                ///printf("WhichBin=%u\n",WhichBin);
        //if(!Polishing) continue;
                unsigned MaxStepsWithoutImprovement = (!Polishing)?MaxStepsWithoutImprovement1:MaxStepsWithoutImprovement2;
                //if we are in a tough spot, increase the MaxSteps
                MaxStepsWithoutImprovement *= (uBF/2+1);
                for(unsigned uED=1*unsigned(Polishing); uED<ErrorDepth+1*unsigned(Polishing); uED++){
                    double Scale = pow(0.5,double(uED));
                    StepsWithoutImprovement=0;
                    unsigned uBin_From =(!Polishing)?WhichBin:WhichBin-BinDepth+1;
                    unsigned uBin_To = (!Polishing)?WhichBin+BinDepth:WhichBin+1;
                    while(StepsWithoutImprovement<MaxStepsWithoutImprovement){
                        Ck_Unfolded[0].Copy(Ck_Sample[0]);
                        //we iterate over the desired bin range and bootstrap a bit
                        for(    unsigned uBin2=uBin_From; uBin2<uBin_To; uBin2++)
                                //unsigned uBin2=WhichBin;
                                //(!Polishing)?uBin2<WhichBin+BinDepth:uBin2>WhichBin-BinDepth;
                                //(!Polishing)?uBin2++:uBin2--)
                                    {
        //printf(" uBin2=%u\n",uBin2);
        //usleep(50e3);
                            OriginalValue = Ck_Sample->GetBinContent(WhichBin);
                            Error = hData_boot->GetBinError(WhichBin+1);
                            if(StepsWithoutImprovement%2) RandomValue = rangen.Uniform(OriginalValue-0.5*Error*Scale,OriginalValue+0.5*Error*Scale);
                            //half of the time we allow for a larger spread
                            else RandomValue = rangen.Gaus(OriginalValue,Error*Scale);
                            Ck_Unfolded->SetBinContent(WhichBin,RandomValue);
                        }
                        Ck_Unfolded->Update(true);
                        CkDec_Unfolded->Update(true);
                        double Chi2=0;
                        double Chi2_Data=0;
                        double Chi2_Fit=0;
                        double Momentum;
                        double CkValCorrected;
                        double CkValData;
                        double Error;
                        double Chi2_BinData;
                        double Chi2_BinFit;
                        double Chi2_BinDataMax=0;
                        double Chi2_BinFitMax=0;
                        double Ndf=0;
                        for(unsigned uBin=0; uBin<hData_boot->GetNbinsX(); uBin++){
                            Momentum = hData_boot->GetBinCenter(uBin+1);
                            if(Momentum<kMin) continue;
                            if(Momentum>kMax) break;
                            CkValCorrected = CkDec_Unfolded->EvalCk(Momentum);
                            CkValData = hData_boot->GetBinContent(uBin+1);
                            Error = hData_boot->GetBinError(uBin+1);
                            //folded - current solution vs data
                            Chi2_BinData = pow((CkValData-CkValCorrected)/(Error),2.);
                            //unfolded - current solution vs fit
                            Chi2_BinFit = pow((Ck_Unfolded->GetBinContent(uBin)-Ck_pL->GetBinContent(uBin))/(Error),2.);
        //printf(" uBin=%u: Chi2_BinData=%.2f; Chi2_BinFit=%.2f for an Error of %.4f\n",uBin,Chi2_BinData,Chi2_BinFit,Error);
        //printf(" Data %.4f; Corrected=%.4f; Fit %.4f; Unfolded %.4f\n",CkValData,CkValCorrected,Ck_pL->GetBinContent(uBin),Ck_Unfolded->GetBinContent(uBin));
                            Chi2 += Chi2_BinData;
                            Chi2 += Chi2_BinFit;
                            Chi2_Data += Chi2_BinData;
                            Chi2_Fit += Chi2_BinFit;
                            if(uBin>=uBin_From&&uBin<uBin_To){
                                if(Chi2_BinData>Chi2_BinDataMax) Chi2_BinDataMax=Chi2_BinData;
                                if(Chi2_BinFit>Chi2_BinFitMax) Chi2_BinFitMax=Chi2_BinFit;
                            }
                            Ndf++;
        ///
                        }

                        //if(Chi2<BestChi2&&Chi2_Data<BestChi2_Data){
                        //if(Chi2_Data<BestChi2_Data){
                        //if(Chi2_Data<BestChi2_Data&&(Chi2_Fit)/double(NumMomBins)<=1.0){
                        //the idea: on the first iteration, make the solution close to the expectation (find the desired global min)
                        //          on the next iterations, loosen up the condition to explore local minima
                        if(Chi2_Data<BestChi2_Data&&Chi2_BinFitMax<=double(uBF+1)*Chi2_BinDataMax){
                        //if(Chi2_BinDataMax<BestChi2_BinData){
                            sprintf(ctext3," --> %.1f",Chi2_Data);
                            if(PC_OUTPUT){
                                printf("\r\033[K%s%s%s",ctext1,ctext2,ctext3);
                                cout << flush;
                            }
                            if(!REFRESHED_INFO){
                                sprintf(ctext2," χ2 = %.1f",Chi2_Data);
                                if(!PC_OUTPUT) cout << ctext2 << flush;
                                REFRESHED_INFO = true;
                            }
                            BestChi2=Chi2;
                            BestChi2_Data=Chi2_Data;
                            BestChi2_BinData=Chi2_BinDataMax;
                            BestChi2_Fit=Chi2_Fit;
                            BestNDF=Ndf;
                            Ck_Best[0].Copy(Ck_Unfolded[0]);
                            Ck_Sample[0].Copy(Ck_Unfolded[0]);
                            StepsWithoutImprovement=0;
                            ///printf(" New best fit: Chi2=%.3f (%.3f+%.3f) at uED=%u\n",Chi2,Chi2_Data,Chi2_Fit,uED);
                            ///printf("  Chi2 Max (%.2f %.2f)\n",Chi2_BinDataMax,Chi2_BinFitMax);

                            if(BestChi2_Data/BestNDF<Perfect_chi2ndf){
                                TerminateIter = true;
                                break;
                            }
                            long long ExeTime = ComputingTime.Stop()/1000000.;
                            if(ExeTime>PatienceMax){
                                TerminateIter = true;
                                break;
                            }
                        }
                        else{
                            //printf(" Current fit: Chi2=%.2f at uED=%u\n",Chi2,uED);
                            StepsWithoutImprovement++;
                        }
                    }//while
                    if(TerminateIter) break;
                }//uED
                if(TerminateIter) break;
            }//uBin
            long long ExeTime = ComputingTime.Stop()/1000000.;
            f_chi2 = BestChi2_Data;
            f_ndf = NumMomBins;
            f_chi2ndf = f_chi2/double(f_ndf);
            f_pval = TMath::Prob(f_chi2,f_ndf);
            f_ns = sqrt(2)*TMath::ErfcInverse(f_pval+1e-32);

            if(!REFRESHED_INFO){
                sprintf(ctext2," χ2 = %.1f",f_chi2);
                cout << ctext2 << flush;
                REFRESHED_INFO = true;
            }
            if(PC_OUTPUT){printf("\r\033[K%s%s",ctext1,ctext2);}
            printf(" --> %.1f\n",f_chi2);
/*
            if(f_chi2ndf<Perfect_chi2ndf){
                ShowTime(ExeTime,strtime,2,true,5);
                printf("  Perfect solution found at %s\n",strtime);
                break;
            }
            if(f_chi2ndf<Okayish_chi2ndf && ExeTime>PatienceForPerfection){
                ShowTime(ExeTime,strtime,2,true,5);
                printf("  Okayish solution found at %s\n",strtime);
                break;
            }
            if(ExeTime>PatienceMax){
                ShowTime(ExeTime,strtime,2,true,5);
                printf("  Over the time limit at %s\n",strtime);
                break;
            }
*/

            //end if we have a better chi2 than the current average,
            //but we are above the average iteration time
            //if(f_chi2ndf<Avg_chi2ndf&&double(ExeTime)>AvgTimePerSeedRemaining){
//printf("%f %f\n",f_chi2ndf,Avg_chi2ndf);

            //terminate ASAP
            if(f_chi2ndf<Perfect_chi2ndf){
                ShowTime(ExeTime,strtime,2,true,5);
                printf("  Perfect solution found at %s\n",strtime);
                break;
            }
            //terminate regardless of the time
            if(f_chi2ndf<VeryGood_chi2ndf){
                ShowTime(ExeTime,strtime,2,true,5);
                printf("  Very good solution found at %s\n",strtime);
                break;
            }
            //terminate if we are slow. If we are fast, try to get to a very good solution
            if(f_chi2ndf<Avg_chi2ndf&&double(ExeTime)>AvgTimePerSeedRemaining&&double(ExeTime)>AvgIterTime){
                ShowTime(ExeTime,strtime,2,true,5);
                printf("  Good solution found at %s\n",strtime);
                break;
            }
            //terminate ASAP if we have exceeded the time limit
            if(ExeTime>PatienceMax){
                ShowTime(ExeTime,strtime,2,true,5);
                printf("  Over the time limit at %s\n",strtime);
                break;
            }


            TotalTime = JobTime.Stop()/1000000.;
            if(TotalTime>JobTimeLimit){
                ShowTime(TotalTime,strtime,2,true,5);
                printf("  Time limit exceeded. Terminating at %s\n",strtime);
                TerminateASAP = true;
                break;
            }
        }//uBF

        long long ExeTime = ComputingTime.Stop()/1000000.;
        Avg_chi2ndf = (Avg_chi2ndf*NumIterDone+f_chi2ndf)/(NumIterDone+1);
        IterTime += ExeTime;
        Chi2Ndf_Unfold->Fill(f_chi2ndf);
        CPU_Unfold->Fill(double(ExeTime)/60.);
        Chi2Ndf_CPU_Unfold->Fill(f_chi2ndf,double(ExeTime)/60.);

        printf("  Final fit quality:\n");
        printf("   chi2/ndf = %.2f/%i (%.2f)\n",f_chi2,f_ndf,Avg_chi2ndf*f_ndf);
        //printf("   pval = %.4f\n",f_pval);
        //printf("   nsigma = %.3f\n",f_ns);

        AcceptSolution = (f_chi2ndf<Unacceptable_chi2ndf);
        if(!AcceptSolution){
            //delete Ck_pL;
            //delete CkDec_pL;
            //delete fit_pL;
            delete Ck_Unfolded;
            delete CkDec_Unfolded;
            delete Ck_Best;
            delete Ck_Sample;
            delete CkDec_Best;
            delete hData_boot;
            printf(" REJECTED\n\n");
            continue;
        }
        else{
            printf(" ACCEPTED\n\n");
        }

        //printf("Ending soon:\n");

        Ck_Best->Update(true);
        CkDec_Best->Update(true);

    //    TFile fOutput(TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Unfolding/Test1/pLambda_Spline_Fit_Unfold2_%.0f.root",BinWidth),
    //                  "recreate");


        hData_pL->SetName(OriginalName+"_Original");
        fit_pL_RAW->SetName(OriginalName+"_fit_pL_RAW");
        fit_pL->SetName(OriginalName+"_fit_pL");


    /*
        TGraph CkTheoryFit;
        CkTheoryFit.SetName(OriginalName+"_TheoryFit");
        TGraph CkSmearedFit;
        CkSmearedFit.SetName(OriginalName+"_SmearedFit");
        TGraph CkTheoryBest;
        //this is the histogram to use
        CkTheoryBest.SetName(OriginalName);
        TGraph CkSmearedBest;
        CkSmearedBest.SetName(OriginalName+"_Smeared");
        for(unsigned uBin=0; uBin<NumMomBins; uBin++){
            double MOM = Ck_pL->GetBinCenter(0,uBin);
            CkTheoryFit.SetPoint(uBin,MOM,Ck_pL->GetBinContent(uBin));
            CkSmearedFit.SetPoint(uBin,MOM,CkDec_pL->EvalCk(MOM));
            CkTheoryBest.SetPoint(uBin,MOM,Ck_Best->GetBinContent(uBin));
            CkSmearedBest.SetPoint(uBin,MOM,CkDec_Best->EvalCk(MOM));
        }
    //printf("The hard part is done\n");
    //usleep(5e6);
        hData_pL->Write();
        fit_pL_RAW->Write();
        fit_pL->Write();
        CkTheoryFit.Write();
        CkSmearedFit.Write();
        CkTheoryBest.Write();
        CkSmearedBest.Write();
    */

    //create ones and rewrite on each iter
        TString OutputFileName;
        //if(SafeFileOpening)
            OutputFileName = OutputFolder+TString::Format("CkSB_pL_%s_Unfolded_%u.root",DataVariation.Data(),uSeed);
        //else OutputFileName = OutputFolder+TString::Format("CkSB_pL_%s_Unfolded.root",DataVariation.Data());

        TH2F* CkBootstrap = NULL;
        TH2F* CkTheoryFit = NULL;
        TH2F* CkSmearedFit = NULL;
        TH2F* CkTheoryBest = NULL;
        TH2F* CkSmearedBest = NULL;

        //TH1F* hCkBootstrap = NULL;
        //TH1F* hCkTheoryFit = NULL;
        //TH1F* hCkSmearedFit = NULL;
        //TH1F* hCkTheoryBest = NULL;
        //TH1F* hCkSmearedBest = NULL;
        //TH1D* hProj = NULL;

        //printf("file_status:\n");
/*
        int InfoFileStatus = file_status(OutputFileName.Data());
        DLM_Timer FileTimer;
        long long FileWaitTime=0;//in micros
        //we wait until the file is closed. If this is not the case in 10s => error
        while(InfoFileStatus==-1&&FileWaitTime<10e6){
            //sleep for 10 ms
            usleep(10e3);
            FileWaitTime = FileTimer.Stop();
            InfoFileStatus = file_status(OutputFileName.Data());
        }

        if(InfoFileStatus==-1){
            printf("\033[1;31mERROR:\033[0m Waited more than 10s to close the file %s\n",OutputFileName.Data());
            abort();
        }
*/
        //printf("InfoFileStatus: %i\n",InfoFileStatus);

        TFile* fOutput = new TFile(OutputFileName,"update");
        //printf("fOutput: %p\n",fOutput);
        //the file does not exist => create the histograms
        //if(InfoFileStatus==0){
            CkBootstrap = new TH2F(OriginalName+"_Bootstrap2D",OriginalName+"_Bootstrap2D",NumMomBins,kMin,kMax,8192,0.5,3);
            CkTheoryFit = new TH2F(OriginalName+"_TheoryFit2D",OriginalName+"_TheoryFit2D",NumMomBins,kMin,kMax,8192,0.5,3);
            CkSmearedFit = new TH2F(OriginalName+"_SmearedFit2D",OriginalName+"_SmearedFit2D",NumMomBins,kMin,kMax,8192,0.5,3);
            CkTheoryBest = new TH2F(OriginalName+"_2D",OriginalName+"_2D",NumMomBins,kMin,kMax,8192,0.5,3);
            CkSmearedBest = new TH2F(OriginalName+"_Smeared2D",OriginalName+"_Smeared2D",NumMomBins,kMin,kMax,8192,0.5,3);

            //hCkBootstrap = (TH1F*)hData_pL->Clone(OriginalName+"_Bootstrap");
            //hCkTheoryFit = (TH1F*)hData_pL->Clone(OriginalName+"_TheoryFit");
            //hCkSmearedFit = (TH1F*)hData_pL->Clone(OriginalName+"_SmearedFit");
            //hCkTheoryBest = (TH1F*)hData_pL->Clone(OriginalName);
            //hCkSmearedBest = (TH1F*)hData_pL->Clone(OriginalName+"_Smeared");
        //}
        //else{
            //CkBootstrap = (TH2F*)fOutput->Get(OriginalName+"_Bootstrap2D");
            //CkTheoryFit = (TH2F*)fOutput->Get(OriginalName+"_TheoryFit2D");
            //CkSmearedFit = (TH2F*)fOutput->Get(OriginalName+"_SmearedFit2D");
            //CkTheoryBest = (TH2F*)fOutput->Get(OriginalName+"_2D");
            //CkSmearedBest = (TH2F*)fOutput->Get(OriginalName+"_Smeared2D");

            //if(!SafeFileOpening){
            //    hCkBootstrap = (TH1F*)fOutput->Get(OriginalName+"_Bootstrap");
            //    hCkTheoryFit = (TH1F*)fOutput->Get(OriginalName+"_TheoryFit");
            //    hCkSmearedFit = (TH1F*)fOutput->Get(OriginalName+"_SmearedFit");
            //    hCkTheoryBest = (TH1F*)fOutput->Get(OriginalName);
            //    hCkSmearedBest = (TH1F*)fOutput->Get(OriginalName+"_Smeared");
            //}
        //}
        //printf("Histos are setup\n");
        for(unsigned uBin=0; uBin<NumMomBins; uBin++){
            double MOM = Ck_pL->GetBinCenter(0,uBin);
    //printf(" MOM=%.0f\n",MOM);
    //printf(" CkBootstrap=%p\n",CkBootstrap);
            CkBootstrap->Fill(MOM,hData_boot->GetBinContent(uBin+1));
    //printf(" CkTheoryFit=%p\n",CkTheoryFit);
            CkTheoryFit->Fill(MOM,Ck_pL->GetBinContent(uBin));
    //printf(" CkSmearedFit=%p\n",CkSmearedFit);
            CkSmearedFit->Fill(MOM,CkDec_pL->EvalCk(MOM));
    //printf(" CkTheoryBest=%p\n",CkTheoryBest);
            CkTheoryBest->Fill(MOM,Ck_Best->GetBinContent(uBin));
    //printf(" CkSmearedBest=%p\n",CkSmearedBest);
            CkSmearedBest->Fill(MOM,CkDec_Best->EvalCk(MOM));
    //printf(" hProj...\n");
            /*
            if(!SafeFileOpening){
                hProj = CkBootstrap->ProjectionY("hProj",uBin+1,uBin+1);
                hCkBootstrap->SetBinContent(uBin+1,hProj->GetMean());
                hCkBootstrap->SetBinError(uBin+1,hProj->GetStdDev());
                delete hProj;

                hProj = CkTheoryFit->ProjectionY("hProj",uBin+1,uBin+1);
                hCkTheoryFit->SetBinContent(uBin+1,hProj->GetMean());
                hCkTheoryFit->SetBinError(uBin+1,hProj->GetStdDev());
                delete hProj;

                hProj = CkSmearedFit->ProjectionY("hProj",uBin+1,uBin+1);
                hCkSmearedFit->SetBinContent(uBin+1,hProj->GetMean());
                hCkSmearedFit->SetBinError(uBin+1,hProj->GetStdDev());
                delete hProj;

                hProj = CkTheoryBest->ProjectionY("hProj",uBin+1,uBin+1);
                hCkTheoryBest->SetBinContent(uBin+1,hProj->GetMean());
                hCkTheoryBest->SetBinError(uBin+1,hProj->GetStdDev());
                delete hProj;

                hProj = CkSmearedBest->ProjectionY("hProj",uBin+1,uBin+1);
                hCkSmearedBest->SetBinContent(uBin+1,hProj->GetMean());
                hCkSmearedBest->SetBinError(uBin+1,hProj->GetStdDev());
                delete hProj;
            }
            */
        }
//printf("About to write:\n");
        //if(uSeed==SEEDmin) hData_pL->Write("",TObject::kOverwrite);
        CkBootstrap->Write("",TObject::kOverwrite);
        CkTheoryFit->Write("",TObject::kOverwrite);
        CkSmearedFit->Write("",TObject::kOverwrite);
        CkTheoryBest->Write("",TObject::kOverwrite);
        CkSmearedBest->Write("",TObject::kOverwrite);

        Chi2Ndf_Unfold->Write("",TObject::kOverwrite);
        CPU_Unfold->Write("",TObject::kOverwrite);
        Chi2Ndf_CPU_Unfold->Write("",TObject::kOverwrite);

        /*
        if(!SafeFileOpening||InfoFileStatus==0){
            hCkBootstrap->Write("",TObject::kOverwrite);
            hCkTheoryFit->Write("",TObject::kOverwrite);
            hCkSmearedFit->Write("",TObject::kOverwrite);
            hCkTheoryBest->Write("",TObject::kOverwrite);
            hCkSmearedBest->Write("",TObject::kOverwrite);
        }
        */

//printf("About to delete:\n");
        //if(InfoFileStatus==0){
        delete CkBootstrap;
        delete CkTheoryFit;
        delete CkSmearedFit;
        delete CkTheoryBest;
        delete CkSmearedBest;
        delete Chi2Ndf_Unfold;
        delete CPU_Unfold;
        delete Chi2Ndf_CPU_Unfold;
            //delete hCkBootstrap;
            //delete hCkTheoryFit;
            //delete hCkSmearedFit;
            //delete hCkTheoryBest;
            //delete hCkSmearedBest;
        //}
//printf("About to delete fOutput:\n");
        delete fOutput;
//printf("About to delete rest:\n");

        delete Ck_Unfolded;
        delete CkDec_Unfolded;
        delete Ck_Best;
        delete CkDec_Best;
        delete Ck_Sample;
        delete hData_boot;

        if(TerminateASAP) break;
//printf("COMPLETE!\n");
    }//uSeed

    delete fit_pL_RAW;
    delete fit_pL;
    delete Ck_pL;
    delete CkDec_pL;
    delete [] Nodes_x;
    delete hData_pL;
    delete [] strtime;
    delete [] ctext1;
    delete [] ctext2;
    delete [] ctext3;
}

void UpdateUnfoldFile(const char* CatsFileFolder, const TString& InputFolderName,
                      const TString& InputFileName, const TString& DataVariation,
                      const int& BinWidth){

    DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder(CatsFileFolder);
    //TH1F* hData_pL = AnalysisObject.GetAliceExpCorrFun("pp13TeV_HM_Dec19","pLambda","_0",TMath::Nint(BinWidth/4)-1,false,-1);
    TH1F* hData_pL = AnalysisObject.GetAliceExpCorrFun("pp13TeV_HM_DimiJun20","pLambda",DataVariation.Data(),TMath::Nint(float(BinWidth)/4.)-1,false,-1);
    TString OriginalName = hData_pL->GetName();
    hData_pL->SetName(OriginalName+"_Original");

    TFile* fOutput = new TFile(InputFolderName+InputFileName,"update");

    //fOutput->ls();
    TH2F* CkBootstrap=NULL;
    TH2F* CkTheoryFit=NULL;
    TH2F* CkSmearedFit=NULL;
    TH2F* CkTheoryBest=NULL;
    TH2F* CkSmearedBest=NULL;

    TH1F* hCkBootstrap=NULL;
    TH1F* hCkTheoryFit=NULL;
    TH1F* hCkSmearedFit=NULL;
    TH1F* hCkTheoryBest=NULL;
    TH1F* hCkSmearedBest=NULL;

    CkBootstrap = (TH2F*)fOutput->Get(OriginalName+"_Bootstrap2D");
    CkTheoryFit = (TH2F*)fOutput->Get(OriginalName+"_TheoryFit2D");
    CkSmearedFit = (TH2F*)fOutput->Get(OriginalName+"_SmearedFit2D");
    CkTheoryBest = (TH2F*)fOutput->Get(OriginalName+"_2D");
    CkSmearedBest = (TH2F*)fOutput->Get(OriginalName+"_Smeared2D");

    hCkBootstrap = (TH1F*)hData_pL->Clone(OriginalName+"_Bootstrap");
    hCkTheoryFit = (TH1F*)hData_pL->Clone(OriginalName+"_TheoryFit");
    hCkSmearedFit = (TH1F*)hData_pL->Clone(OriginalName+"_SmearedFit");
    hCkTheoryBest = (TH1F*)hData_pL->Clone(OriginalName);
    hCkSmearedBest = (TH1F*)hData_pL->Clone(OriginalName+"_Smeared");

    const unsigned NumMomBins = CkTheoryBest->GetXaxis()->GetNbins();
    TH1D* hProj;

    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        hProj = CkBootstrap->ProjectionY("hProj",uBin+1,uBin+1);
        hCkBootstrap->SetBinContent(uBin+1,hProj->GetMean());
        hCkBootstrap->SetBinError(uBin+1,hProj->GetStdDev());
        delete hProj;

        hProj = CkTheoryFit->ProjectionY("hProj",uBin+1,uBin+1);
        hCkTheoryFit->SetBinContent(uBin+1,hProj->GetMean());
        hCkTheoryFit->SetBinError(uBin+1,hProj->GetStdDev());
        delete hProj;

        hProj = CkSmearedFit->ProjectionY("hProj",uBin+1,uBin+1);
        hCkSmearedFit->SetBinContent(uBin+1,hProj->GetMean());
        hCkSmearedFit->SetBinError(uBin+1,hProj->GetStdDev());
        delete hProj;

        hProj = CkTheoryBest->ProjectionY("hProj",uBin+1,uBin+1);
        hCkTheoryBest->SetBinContent(uBin+1,hProj->GetMean());
        hCkTheoryBest->SetBinError(uBin+1,hProj->GetStdDev());
        delete hProj;

        hProj = CkSmearedBest->ProjectionY("hProj",uBin+1,uBin+1);
        hCkSmearedBest->SetBinContent(uBin+1,hProj->GetMean());
        hCkSmearedBest->SetBinError(uBin+1,hProj->GetStdDev());
        delete hProj;
    }

    hData_pL->Write("",TObject::kOverwrite);
    hCkBootstrap->Write("",TObject::kOverwrite);
    hCkTheoryFit->Write("",TObject::kOverwrite);
    hCkSmearedFit->Write("",TObject::kOverwrite);
    hCkTheoryBest->Write("",TObject::kOverwrite);
    hCkSmearedBest->Write("",TObject::kOverwrite);

    delete hData_pL;
    delete fOutput;
}



double FitSmearMatrixNormHelper(double* par){
    double& x = par[0];
    double& Amplitude = par[1];
    double& Mean = par[2];
    double& Sigma = par[3];
    double cpi = TMath::Pi();
    //double Norm = (2.*exp(Mean/Sigma)*sqrt(1./Sigma/Sigma))/
    //(sqrt(cpi)*(exp(Mean/Sigma)+
    //            exp(Mean/Sigma)*TMath::Erf(Mean*sqrt(1./Sigma/Sigma))+
    //            exp(0.25)*TMath::Erfc(0.5*sqrt(1./Sigma/Sigma)*(-2.*Mean+Sigma))));
    return Amplitude*exp(-pow((x-Mean)/Sigma,2.))*(1.-exp(-pow(x/Sigma,2.)));
}
double FitSmearMatrixNorm(double* x, double* par){
    double PAR[4];
    PAR[0] = *x; PAR[1] = par[0]; PAR[2] = par[1]; PAR[3] = par[2];
    double& Amplitude = par[0];
    double& Mean = par[1];
    double& Sigma = par[2];
    const double cpi = TMath::Pi();
    double IntMin = Mean-5*Sigma;
    if(IntMin<=0) IntMin=1e-6;
    double IntMax = Mean+5*Sigma;
    DLM_INT_SetFunction(FitSmearMatrixNormHelper,PAR,0);
    double Norm = DLM_INT_aSimpsonWiki(IntMin,IntMax);
    double Value = FitSmearMatrixNormHelper(PAR);
    return Norm*Value;
}

double FitSmearMatrixRaw(double* x, double* par){
    double& Amplitude = par[0];
    double& Mean = par[1];
    double& Sigma = par[2];
    return Amplitude*exp(-pow((*x-Mean)/Sigma,2.))*(1.-exp(-pow(*x/Sigma,2.)));
    //return Amplitude*exp(-pow((*x-Mean)/Sigma,2.));
}

//
void ParametrizeSmearMatrix1_pL(){

    DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder("/mnt/Ubuntu_Data/CernBox/Sync/CatsFiles");
    TH2F* hReso = AnalysisObject.GetResolutionMatrix("pp13TeV_HM_Dec19","pLambda");
    printf("hReso = %p\n",hReso);
    //hReso->Sumw2();

    TH1D* hProj = hReso->ProjectionY(TString::Format("hReso"),1,1);
    TF1* fProj = new TF1(TString::Format("fReso"),FitSmearMatrixRaw,0,1000,3);
    fProj->SetNpx(1024);
    float Amplitude;
    float Mean;
    float Sigma;
    float Chi2_NDF;
    float Momentum;

    TFile* fOutput = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/ParametrizeSmearMatrix1_pL/fOutput.root","recreate");
    printf("fOutput = %p\n",fOutput);
    TTree* pSmearTree = new TTree("pSmearTree","pSmearTree");
    printf("pSmearTree = %p\n",pSmearTree);
    pSmearTree->Branch("hProj","TH1D",&hProj,32000,0);//
    pSmearTree->Branch("fProj","TF1",&fProj,32000,0);//
    pSmearTree->Branch("Momentum", &Momentum, "Momentum/F");//
    pSmearTree->Branch("Amplitude", &Amplitude, "Amplitude/F");//
    pSmearTree->Branch("Mean", &Mean, "Mean/F");//
    pSmearTree->Branch("Sigma", &Sigma, "Sigma/F");//
    pSmearTree->Branch("Chi2_NDF", &Chi2_NDF, "Chi2_NDF/F");
    printf("TTree ready\n");

    delete hProj;
    delete fProj;

    const unsigned NumBins = hReso->GetXaxis()->GetNbins();
    for(unsigned uBinX=0; uBinX<NumBins; uBinX++){
        printf("uBinX=%u\n",uBinX);
        Momentum = hReso->GetXaxis()->GetBinCenter(uBinX+1);
        hProj = hReso->ProjectionY(TString::Format("hReso_%u",uBinX),uBinX+1,uBinX+1);
        hProj->Sumw2();
        hProj->Scale(1./hProj->Integral(),"width");
        double FitFrom = hProj->GetMean()-3.*hProj->GetStdDev();
        if(FitFrom<0) FitFrom=0;
        double FitTo = hProj->GetMean()+3.*hProj->GetStdDev();
        fProj = new TF1(TString::Format("fReso_%u",uBinX),FitSmearMatrixRaw,FitFrom,FitTo,3);
        fProj->SetParameter(0,0.05);
        fProj->SetParameter(1,Momentum);
        fProj->SetParLimits(1,0.5*Momentum,2.*Momentum);
        fProj->SetParameter(2,5);
        fProj->SetParLimits(2,0,20);
        hProj->Fit(fProj,"Q, S, N, R, M");
        Amplitude = fProj->GetParameter(0);
        Mean = fProj->GetParameter(1);
        Sigma = fProj->GetParameter(2);
        Chi2_NDF = fProj->GetChisquare()/double(fProj->GetNDF());
//printf("Chi2_NDF=%f\n",Chi2_NDF);
        pSmearTree->Fill();
        delete hProj;
        delete fProj;
    }

    pSmearTree->Write();

    delete hReso;
    delete pSmearTree;
    delete fOutput;
}


//for each potential/BL/Sigma variation:
// * find out the chi2 for any of the different r,α,λ,ω variations
// * printout the best chi2 FOR ALL potentials/BL/Sigma, and use it as a baseline to compare all other hypothesis, assuming 4 fit parameters
// * printout the nsigma for each potential,BL,Sigma var, and also for each specific r,α,λ,ω
// * these are many combinations, to reduce the output when presenting, just print:
//      * by how much a specific value for each r,α,λ,ω observables are excluded (for given hypothesis)
void pLambda_Study_Hypotheses(){
    TString FolderName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/pL_SystematicsMay2020/BatchFarm/290520_Levy/";
    TString BaseInputName = "Merged_pp13TeV_HM_Dec19_";

    enum BLTYPE              { pol0s,pol1s,pol2s,pol3s,dpol2s,dpol3s,dpol4s,pol2e,pol3e,dpol2e,dpol3e,dpol4e };
    //const bool IncludeBl[] = { true, true, true, false,false, true,  false, false,false,false, true,  false};
    const unsigned NumBaselines = 8;
    int* BlFlag = new int [NumBaselines];
    bool* IncludeBl = new bool [NumBaselines];
    BlFlag[0] = pol0s; IncludeBl[0] = false;
    BlFlag[1] = pol1s; IncludeBl[1] = false;
    BlFlag[2] = pol2s; IncludeBl[2] = false;
    BlFlag[3] = dpol2s; IncludeBl[3] = false;
    BlFlag[4] = pol3e; IncludeBl[4] = false;
    BlFlag[5] = dpol3s; IncludeBl[5] = false;
    BlFlag[6] = dpol3e; IncludeBl[6] = true;
    BlFlag[7] = dpol4e; IncludeBl[7] = false;
    TString* BlName1 = new TString [NumBaselines+1];
    TString* BlName2 = new TString [NumBaselines+1];
    TString* BlDescr = new TString [NumBaselines+1];
    for(unsigned uBl=0; uBl<=NumBaselines; uBl++) BlName1[uBl] = "";
    for(unsigned uBl=0; uBl<=NumBaselines; uBl++) BlName2[uBl] = "";

    BlName1[0] = "Constant baseline";
    BlName2[0] = "Standard fit range";
    BlDescr[0] = "pol0s";

    BlName1[1] = "Linear baseline";
    BlName2[1] = "Standard fit range";
    BlDescr[1] = "pol1s";

    BlName1[2] = "Quadratic baseline";
    BlName2[2] = "Standard fit range";
    BlDescr[2] = "pol2s";

    BlName1[3] = "Quadratic baseline (constrained)";
    BlName2[3] = "Standard fit range";
    BlDescr[3] = "dpol2s";

    BlName1[4] = "Cubic baseline";
    BlName2[4] = "Extended fit range";
    BlDescr[4] = "pol3e";

    BlName1[5] = "Cubic baseline (constrained)";
    BlName2[5] = "Standard fit range";
    BlDescr[5] = "dpol3s";

    BlName1[6] = "Cubic baseline (constrained)";
    BlName2[6] = "Extended fit range";
    BlDescr[6] = "dpol3e";

    BlName1[7] = "Quartic baseline (constrained)";
    BlName2[7] = "Extended fit range";
    BlDescr[7] = "dpol4e";

    BlName1[8] = "All baselines";
    BlName2[8] = "";
    BlDescr[8] = "AllBl";

    const unsigned NumPotVars = 8;
    TString* PotName1 = new TString [NumPotVars+1];
    PotName1[0] = "NLO13 (500)";
    PotName1[1] = "NLO13 (550)";
    PotName1[2] = "NLO13 (600)";
    PotName1[3] = "NLO13 (650)";
    PotName1[4] = "NLO19 (500)";
    PotName1[5] = "NLO19 (550)";
    PotName1[6] = "NLO19 (600)";
    PotName1[7] = "NLO19 (650)";
    PotName1[8] = "NLO13/19 (500-650)";
    TString* PotName2 = new TString [NumPotVars+1];
    for(unsigned uPot=0; uPot<=NumPotVars; uPot++) PotName2[uPot] = "s,d waves";
    int* PotFlag = new int [NumPotVars];
    bool* IncludePot = new bool [NumPotVars];
    PotFlag[0] = 1500; IncludePot[0] = true;
    PotFlag[1] = 1550; IncludePot[1] = true;
    PotFlag[2] = 1600; IncludePot[2] = true;
    PotFlag[3] = 1650; IncludePot[3] = true;
    PotFlag[4] = 11500; IncludePot[4] = true;
    PotFlag[5] = 11550; IncludePot[5] = true;
    PotFlag[6] = 11600; IncludePot[6] = true;
    PotFlag[7] = 11650; IncludePot[7] = true;
    TString* PotDescr = new TString [NumPotVars+1];
    PotDescr[0] = "NLO13-500";
    PotDescr[1] = "NLO13-550";
    PotDescr[2] = "NLO13-600";
    PotDescr[3] = "NLO13-650";
    PotDescr[4] = "NLO19-500";
    PotDescr[5] = "NLO19-550";
    PotDescr[6] = "NLO19-600";
    PotDescr[7] = "NLO19-650";
    PotDescr[8] = "NLO";

    const unsigned NumSigVars = 3;
    int* SigFlag = new int [NumSigVars];
    bool* IncludeSig = new bool [NumSigVars];
    SigFlag[0] = 0; IncludeSig[0] = true;
    SigFlag[1] = 1; IncludeSig[1] = true;
    SigFlag[2] = 2; IncludeSig[2] = true;
    TString* SigName1 = new TString [NumSigVars+1];
    SigName1[0] = "Flat";
    SigName1[1] = "#chiEFT (NLO)";
    SigName1[2] = "ESC 16";
    SigName1[3] = "All #Sigma models";
    TString* SigDescr = new TString [NumSigVars+1];
    SigDescr[0] = "Flat";
    SigDescr[1] = "chEFT";
    SigDescr[2] = "ESC16";
    SigDescr[3] = "AllSigma";

    const unsigned NumFitPars = 4;
    DLM_Histo<float> dlm_chi2;
    DLM_Histo<float> dlm_ndf;


    dlm_chi2.SetUp(3+NumFitPars);
    dlm_chi2.SetUp(0,8,-0.5,4.5);//pot
    dlm_chi2.SetUp(1,1,-0.5,0.5);//bl
    dlm_chi2.SetUp(2,3,-0.5,2.5);//sig
    dlm_chi2.SetUp(3,5,0.92,1.12);//r
    dlm_chi2.SetUp(4,5,1.55,2.05);//α
    dlm_chi2.SetUp(5,5,0.405,0.505);//λ
    dlm_chi2.SetUp(6,5,0.19,0.47);//ω
    dlm_chi2.Initialize();
    dlm_chi2.SetBinContentAll(1e6);

    dlm_ndf.SetUp(3+NumFitPars);
    dlm_ndf.SetUp(0,8,-0.5,4.5);//pot
    dlm_ndf.SetUp(1,1,-0.5,0.5);//bl
    dlm_ndf.SetUp(2,3,-0.5,2.5);//sig
    dlm_ndf.SetUp(3,5,0.92,1.12);//r
    dlm_ndf.SetUp(4,5,1.55,2.05);//α
    dlm_ndf.SetUp(5,5,0.405,0.505);//λ
    dlm_ndf.SetUp(6,5,0.19,0.47);//ω
    dlm_ndf.Initialize();
    dlm_ndf.SetBinContentAll(1);

    unsigned* WhichBin = new unsigned [3+NumFitPars];
    unsigned WhichTotBin;
    for(unsigned uPot=0; uPot<NumPotVars; uPot++){
        bool BinCentersAreOkay = false;
        if(!IncludePot[uPot]) continue;
        WhichBin[uPot] = uPot;
        for(unsigned uBl=0; uBl<NumBaselines; uBl++){
            if(!IncludeBl[uBl]) continue;
            WhichBin[uBl] = uBl;
            for(unsigned uSig=0; uSig<NumSigVars; uSig++){
                if(!IncludeSig[uSig]) continue;

                DLM_Histo<float> dlmPBS_chi2;
                DLM_Histo<float> dlmPBS_ndf;
                DLM_Histo<float> dlmPBS_num;

                dlmPBS_chi2.SetUp(NumFitPars);
                dlmPBS_chi2.SetUp(0,5,0.92,1.12);//r
                dlmPBS_chi2.SetUp(1,5,1.55,2.05);//α
                dlmPBS_chi2.SetUp(2,5,0.405,0.505);//λ
                dlmPBS_chi2.SetUp(3,5,0.19,0.47);//ω
                dlmPBS_chi2.Initialize();
                dlmPBS_chi2.SetBinContentAll(1e6);

                dlmPBS_ndf.SetUp(NumFitPars);
                dlmPBS_ndf.SetUp(0,5,0.92,1.12);//r
                dlmPBS_ndf.SetUp(1,5,1.55,2.05);//α
                dlmPBS_ndf.SetUp(2,5,0.405,0.505);//λ
                dlmPBS_ndf.SetUp(3,5,0.19,0.47);//ω
                dlmPBS_ndf.Initialize();
                dlmPBS_ndf.SetBinContentAll(1);

                dlmPBS_num.SetUp(NumFitPars);
                dlmPBS_num.SetUp(0,5,0.92,1.12);//r
                dlmPBS_num.SetUp(1,5,1.55,2.05);//α
                dlmPBS_num.SetUp(2,5,0.405,0.505);//λ
                dlmPBS_num.SetUp(3,5,0.19,0.47);//ω
                dlmPBS_num.Initialize();
                dlmPBS_num.SetBinContentAll(1);

                WhichBin[uSig] = uSig;
                TFile* fInput = new TFile(FolderName+BaseInputName+
                                TString::Format("POT%i_BL%i_SIG%i.root",PotFlag[uPot],BlFlag[uBl],SigFlag[uSig]),"read");
                TTree* plambdaTree = (TTree*)fInput->Get("plambdaTree");

                TGraphErrors* gData = NULL;TGraphErrors* fitData = NULL;TGraphErrors* fitFemto = NULL;TGraphErrors* fitGenuine = NULL;
                TGraphErrors* fitSignal_pL = NULL;TGraphErrors* fitSignal_pL_pS0 = NULL;TGraphErrors* fitSignal_pL_pXim = NULL;TGraphErrors* fitBaseline = NULL;
                TString* DataSample = NULL;TString* SourceDescription = NULL;
                unsigned SEED;unsigned BASELINE_VAR;int POT_VAR;
                bool DataSyst;bool FitSyst;bool Bootstrap;
                float SourceSize;float SourceAlpha;float CuspWeight;float CkConv;float CkCutOff;
                float Norm;float pp1;float pp2;float pp3;float pp4;
                float kFemtoMin;float kFemtoMax;float kLongMin;float kLongMax;
                float lam_L_genuine;float lam_L_Sig0;float lam_L_Xim;float lam_L_Xi0;float lam_L_Flat;
                float lam_S0_genuine;float lam_S0_Flat;float lam_Xim_genuine;float lam_Xim_Flat;
                float pval;unsigned WhichData;bool DefaultVariation;int Sigma0_Feed;int Xim_Feed;

                plambdaTree->SetBranchAddress("gData",&gData);
                plambdaTree->SetBranchAddress("fitData",&fitData);
                plambdaTree->SetBranchAddress("fitFemto",&fitFemto);
                plambdaTree->SetBranchAddress("fitGenuine",&fitGenuine);
                plambdaTree->SetBranchAddress("fitSignal_pL",&fitSignal_pL);
                plambdaTree->SetBranchAddress("fitSignal_pL_pS0",&fitSignal_pL_pS0);
                plambdaTree->SetBranchAddress("fitSignal_pL_pXim",&fitSignal_pL_pXim);
                plambdaTree->SetBranchAddress("fitBaseline",&fitBaseline);
                plambdaTree->SetBranchAddress("DataSample",&DataSample);
                plambdaTree->SetBranchAddress("SourceDescription",&SourceDescription);
                plambdaTree->SetBranchAddress("SEED",&SEED);
                plambdaTree->SetBranchAddress("BASELINE_VAR",&BASELINE_VAR);
                plambdaTree->SetBranchAddress("POT_VAR",&POT_VAR);
                plambdaTree->SetBranchAddress("Sigma0_Feed",&Sigma0_Feed);
                plambdaTree->SetBranchAddress("Xim_Feed",&Xim_Feed);
                plambdaTree->SetBranchAddress("DataSyst",&DataSyst);
                plambdaTree->SetBranchAddress("FitSyst", &FitSyst);//
                plambdaTree->SetBranchAddress("Bootstrap", &Bootstrap);//
                plambdaTree->SetBranchAddress("SourceSize",&SourceSize);
                plambdaTree->SetBranchAddress("SourceAlpha",&SourceAlpha);
                plambdaTree->SetBranchAddress("CuspWeight",&CuspWeight);
                plambdaTree->SetBranchAddress("CkCutOff",&CkCutOff);
                plambdaTree->SetBranchAddress("CkConv",&CkConv);
                plambdaTree->SetBranchAddress("Norm",&Norm);
                plambdaTree->SetBranchAddress("pp1",&pp1);
                plambdaTree->SetBranchAddress("pp2",&pp2);
                plambdaTree->SetBranchAddress("pp3",&pp3);
                plambdaTree->SetBranchAddress("pp4",&pp4);
                plambdaTree->SetBranchAddress("kFemtoMin",&kFemtoMin);
                plambdaTree->SetBranchAddress("kFemtoMax",&kFemtoMax);
                plambdaTree->SetBranchAddress("kLongMin",&kLongMin);
                plambdaTree->SetBranchAddress("kLongMax",&kLongMax);
                plambdaTree->SetBranchAddress("lam_L_genuine",&lam_L_genuine);
                plambdaTree->SetBranchAddress("lam_L_Sig0",&lam_L_Sig0);
                plambdaTree->SetBranchAddress("lam_L_Xim",&lam_L_Xim);
                plambdaTree->SetBranchAddress("lam_L_Xi0",&lam_L_Xi0);
                plambdaTree->SetBranchAddress("lam_L_Flat",&lam_L_Flat);
                plambdaTree->SetBranchAddress("lam_S0_genuine",&lam_S0_genuine);
                plambdaTree->SetBranchAddress("lam_S0_Flat",&lam_S0_Flat);
                plambdaTree->SetBranchAddress("lam_Xim_genuine",&lam_Xim_genuine);
                plambdaTree->SetBranchAddress("lam_Xim_Flat",&lam_Xim_Flat);
                plambdaTree->SetBranchAddress("WhichData",&WhichData);
                plambdaTree->SetBranchAddress("pval",&pval);
                plambdaTree->SetBranchAddress("DefaultVariation",&DefaultVariation);

                const unsigned NumEntries = plambdaTree->GetEntries();

                for(unsigned uEntry=0; uEntry<NumEntries; uEntry++){
                    plambdaTree->GetEntry(uEntry);
                    //if(!BinCentersAreOkay){
                        dlm_chi2.SetBinCenter(3,dlm_chi2.GetBin(3,SourceSize),SourceSize);
                        dlm_chi2.SetBinCenter(4,dlm_chi2.GetBin(4,SourceAlpha),SourceAlpha);
                        dlm_chi2.SetBinCenter(5,dlm_chi2.GetBin(5,lam_L_genuine),lam_L_genuine);
                        dlm_chi2.SetBinCenter(6,dlm_chi2.GetBin(6,CuspWeight),CuspWeight);

                        dlmPBS_chi2.SetBinCenter(0,dlmPBS_chi2.GetBin(0,SourceSize),SourceSize);
                        dlmPBS_chi2.SetBinCenter(1,dlmPBS_chi2.GetBin(1,SourceAlpha),SourceAlpha);
                        dlmPBS_chi2.SetBinCenter(2,dlmPBS_chi2.GetBin(2,lam_L_genuine),lam_L_genuine);
                        dlmPBS_chi2.SetBinCenter(3,dlmPBS_chi2.GetBin(3,CuspWeight),CuspWeight);

                        dlmPBS_ndf.SetBinCenter(0,dlmPBS_ndf.GetBin(0,SourceSize),SourceSize);
                        dlmPBS_ndf.SetBinCenter(1,dlmPBS_ndf.GetBin(1,SourceAlpha),SourceAlpha);
                        dlmPBS_ndf.SetBinCenter(2,dlmPBS_ndf.GetBin(2,lam_L_genuine),lam_L_genuine);
                        dlmPBS_ndf.SetBinCenter(3,dlmPBS_ndf.GetBin(3,CuspWeight),CuspWeight);

                        dlmPBS_num.SetBinCenter(0,dlmPBS_num.GetBin(0,SourceSize),SourceSize);
                        dlmPBS_num.SetBinCenter(1,dlmPBS_num.GetBin(1,SourceAlpha),SourceAlpha);
                        dlmPBS_num.SetBinCenter(2,dlmPBS_num.GetBin(2,lam_L_genuine),lam_L_genuine);
                        dlmPBS_num.SetBinCenter(3,dlmPBS_num.GetBin(3,CuspWeight),CuspWeight);
                    //}
                    WhichBin[3] = dlm_chi2.GetBin(3,SourceSize);
                    WhichBin[4] = dlm_chi2.GetBin(4,SourceAlpha);
                    WhichBin[5] = dlm_chi2.GetBin(5,lam_L_genuine);
                    WhichBin[6] = dlm_chi2.GetBin(6,CuspWeight);


if(WhichBin[5]!=WhichBin[6]) continue;
//if(WhichBin[3]!=2&&WhichBin[4]!=4) continue;
if(WhichBin[3]==0||WhichBin[3]==4) continue;
if(WhichBin[4]!=3) continue;
if(WhichBin[5]==0||WhichBin[5]==4) continue;
if(WhichBin[6]==0||WhichBin[6]==4) continue;

                    float Chi2=0;float NDF=0;
                    double Momentum;double DataValue;double DataErr;double FitValue;

                    for(unsigned uBin=0; uBin<fitData->GetN(); uBin++){
                        fitData->GetPoint(uBin,Momentum,FitValue);
                        if(Momentum>300) break;
                        gData->GetPoint(uBin,Momentum,DataValue);
                        DataErr = gData->GetErrorY(uBin);
                        Chi2 += pow((DataValue-FitValue)/DataErr,2.);
                        NDF++;
                    }

                    if(dlmPBS_chi2.GetBinContent(&WhichBin[3])==1e6){
                        dlm_chi2.SetBinContent(WhichBin,Chi2);
                        dlm_ndf.SetBinContent(WhichBin,NDF);

                        dlmPBS_chi2.SetBinContent(&WhichBin[3],Chi2);
                        dlmPBS_ndf.SetBinContent(&WhichBin[3],NDF);
                        dlmPBS_num.SetBinContent(&WhichBin[3],1);
                    }

                    else{
                        WhichTotBin = dlm_chi2.GetTotBin(WhichBin);
                        dlm_chi2.Add(WhichTotBin,Chi2);
                        dlm_ndf.Add(WhichTotBin,NDF);

                        WhichTotBin = dlmPBS_chi2.GetTotBin(&WhichBin[3]);
//printf("WhichTotBin=%u %f %f\n",WhichTotBin,Chi2,NDF);
//printf(" %f ->",dlmPBS_chi2.GetBinContent(WhichTotBin));
                        dlmPBS_chi2.Add(WhichTotBin,Chi2);
//printf(" %f / %f = %f\n",dlmPBS_chi2.GetBinContent(WhichTotBin),dlmPBS_num.GetBinContent(WhichTotBin),dlmPBS_chi2.GetBinContent(WhichTotBin)/dlmPBS_num.GetBinContent(WhichTotBin));
                        dlmPBS_ndf.Add(WhichTotBin,NDF);
                        dlmPBS_num.Add(WhichTotBin,1.0);
//usleep(100e3);
                    }

                }

                dlmPBS_chi2.DivideHistoBinByBin(dlmPBS_num,false);
                dlmPBS_ndf.DivideHistoBinByBin(dlmPBS_num,false);

                float BestChi2;
                float BestNdf;
                unsigned* BestBin = NULL;
                unsigned* BestWhichBin = new unsigned [NumFitPars];
                unsigned NumBesties;
                NumBesties = dlmPBS_chi2.FindMinima(BestChi2,BestBin);
                dlmPBS_chi2.GetBinCoordinates(BestBin[0],BestWhichBin);
                BestNdf = dlmPBS_ndf.GetBinContent(BestBin[0]);
                printf("Mainstream hypothesis:\n  %s p-Λ, %s baseline, %s p-Σ0\n",PotDescr[uPot].Data(),BlDescr[uBl].Data(),SigDescr[uSig].Data());
                printf("  Best solution: r=%.2f(%u); α=%.2f(%u); λ=%.2f(%u); ω=%.2f(%u)\n",
                       dlmPBS_chi2.GetBinCenter(0,BestWhichBin[0]),BestWhichBin[0],
                       dlmPBS_chi2.GetBinCenter(1,BestWhichBin[1]),BestWhichBin[1],
                       dlmPBS_chi2.GetBinCenter(2,BestWhichBin[2]),BestWhichBin[2],
                       dlmPBS_chi2.GetBinCenter(3,BestWhichBin[3]),BestWhichBin[3]);
                printf("    with χ2/NDF = %.1f/%.0f = %.2f\n",BestChi2,BestNdf,BestChi2/BestNdf);

                delete plambdaTree;
                delete fInput;
                delete [] BestBin;
                delete [] BestWhichBin;
                BinCentersAreOkay = true;
            }
        }
    }

    delete [] BlFlag;
    delete [] IncludeBl;
    delete [] BlName1;
    delete [] BlName2;
    delete [] BlDescr;
    delete [] PotDescr;
    delete [] PotName1;
    delete [] PotName2;
    delete [] PotFlag;
    delete [] IncludePot;
    delete [] SigFlag;
    delete [] IncludeSig;
    delete [] SigName1;
    delete [] SigDescr;
    delete [] WhichBin;
}

void CompareTo_pp(){

    TString InputFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/CompareTo_pp/OutFileVarpp.root";
    //TString InputFileName = "/home/dmihaylov/Downloads/OutFileVarpp.root";

    const double kMin=0;
    //Laura said 200
    const double kMax=350;

    double Chi2=0;
    double NDF=0;
    unsigned NumIter=0;

    const unsigned Which_ModPL = 2;
    const unsigned Which_PolBL = 0;
    //0.1915, 0.2014 0.2115 err 0.001
    //0 is integrated
    const float Which_lam_ppL = 0.0;
    //const float Which_lam_ppL = 0.1915;
    //const float Which_lam_ppL = 0.2014;
    //const float Which_lam_ppL = 0.2115;

    unsigned ModPL;
    unsigned PolBL;
    float lam_ppL;

    TFile* fInput = new TFile(InputFileName,"read");
    printf("fInput = %p\n",fInput);
    TTree* AndiTree = (TTree*)fInput->Get("ppTree");
    TGraph* FitResult = NULL;
    TH1F* CorrHist = NULL;
    printf("AndiTree = %p\n",AndiTree);
    //AndiTree->Print();

	AndiTree->SetBranchAddress("CorrHist",&CorrHist);
	AndiTree->SetBranchAddress("FitResult",&FitResult);
	AndiTree->SetBranchAddress("ModPL",&ModPL);
	AndiTree->SetBranchAddress("PolBL",&PolBL);
	AndiTree->SetBranchAddress("lam_ppL",&lam_ppL);

	unsigned NumEntries = AndiTree->GetEntries();
    double CkData;
    double CkFit;
    double xAxis;
    double CkErr;
    double NumChi2Entries=0;

    printf("Iterating over %u events\n",NumEntries);

	for(unsigned uEntry=0; uEntry<NumEntries; uEntry++){
        //printf("uEntry=%u\n",uEntry);
        fInput->cd();
        AndiTree->GetEntry(uEntry);
        //printf(" got it!\n");
        if(ModPL!=Which_ModPL) continue;
        if(PolBL!=Which_PolBL) continue;
        if(Which_lam_ppL&&(lam_ppL<Which_lam_ppL-0.001||lam_ppL>Which_lam_ppL+0.001)) continue;
        //printf(" CorrHist=%p\n",CorrHist);
        for(unsigned uBin=0; uBin<CorrHist->GetNbinsX(); uBin++){
            double Momentum = CorrHist->GetBinCenter(uBin+1);
            if(Momentum>kMax) break;
            CkData = CorrHist->GetBinContent(uBin+1);
            CkErr = CorrHist->GetBinError(uBin+1);
            FitResult->GetPoint(uBin,xAxis,CkFit);
            if(uEntry==0){
                printf("%u: k*=%.2f %.3f+/-%.3f vs %.3f\n",uBin,Momentum,CkData,CkErr,CkFit);
            }
            Chi2 += pow((CkData-CkFit)/CkErr,2.);
            NDF++;
        }
        NumChi2Entries++;
	}

    Chi2 /= NumChi2Entries;
    NDF /= NumChi2Entries;
    double nsigma = sqrt(2)*TMath::ErfcInverse(TMath::Prob(Chi2,TMath::Nint(NDF)));

    printf("<Chi2/NDF> = %.2f/%.0f = %.2f\n",Chi2,NDF,Chi2/NDF);
    printf("nsigma = %.2f\n",nsigma);

    delete fInput;
}


double Fit_IMS_Signal(double* xVal, double* pars){
    double& IM = *xVal;
    double& Ampl1 = pars[0];
    double& Mean1 = pars[1];
    double& Stdv1 = pars[2];
    double& Ampl2 = pars[3];
    double Mean2 = pars[4]; if(Mean2==0) Mean2=Mean1;
    double& Stdv2 = pars[5];
    return Ampl1*TMath::Gaus(IM,Mean1,Stdv1,1)+Ampl2*TMath::Gaus(IM,Mean2,Stdv2,1);
}
double Fit_IMS_Splines(double* xVal, double* pars){
    double Peak = Fit_IMS_Signal(xVal,pars);
    double Splines = DLM_FITTER2_FUNCTION_SPLINE3(xVal,&pars[6]);
    return Splines+Peak;
}
double Fit_IMS_Signal3(double* xVal, double* pars){
    double& IM = *xVal;
    double& Ampl1 = pars[0];
    double& Mean1 = pars[1];
    double& Stdv1 = pars[2];
    double& Ampl2 = pars[3];
    double Mean2 = pars[4]; if(Mean2==0) Mean2=Mean1;
    double& Stdv2 = pars[5];
    double& Ampl3 = pars[6];
    double Mean3 = pars[7]; if(Mean3==0) Mean3=Mean1;
    double& Stdv3 = pars[8];
    return Ampl1*TMath::Gaus(IM,Mean1,Stdv1,1)+Ampl2*TMath::Gaus(IM,Mean2,Stdv2,1)+Ampl3*TMath::Gaus(IM,Mean3,Stdv3,1);
}
double Fit_IMS_Splines3(double* xVal, double* pars){
    double Peak = Fit_IMS_Signal3(xVal,pars);
    double Splines = DLM_FITTER2_FUNCTION_SPLINE3(xVal,&pars[9]);
    return Splines+Peak;
}
double Fit_IMS_Pol2(double* xVal, double* pars){
    double& x = *xVal;
    double Peak = Fit_IMS_Signal(xVal,pars);
    double Pol2 = pars[6]+pars[7]*x+pars[8]*x*x;
    return Pol2+Peak;
}

//0 is se, 1 is me
void Purity_vs_kstar(const int& MEflag, const unsigned& FirstBoot, const unsigned& LastBoot){
//PurityTree->Draw("Purity:Momentum>>h2(20,0,480,64,0.94, 0.98)", "Scenario==2", "COLZ")
    //const double ME_yield_bias = 17;
    const double kMax = 48;
    //const unsigned BootstrapIter = 128;
    const double pdgLambda = 1115.683;
    const double imWindow = 4.0;

    const int RebinKstar = 6;
    const int RebinIMS = 8;

    TString OutputFileName = TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/"
                            "pLambda_1/Purity_vs_kstar/Test/fOutput%s_%u_%u.root",
                            MEflag?"ME":"SE",FirstBoot,LastBoot);


/*
    TString InputFileName = "/mnt/Ubuntu_Data/CernBox/Sync/CatsFiles_Dimi/pLambda/AnalysisResults_Purity01062020.root";
    TString DirectoryName = "HMDimiResultsQA0";
    TString ListName1 = "HMDimiResultsQA0";
    TString ListName2 = "PairQA";
    TString ListNamePP = "QA_Particle0_Particle2";
    TString ListNameAPAP = "QA_Particle1_Particle3";
    TString HistoNamePP = "MassQA_Particle2_2";
    TString HistoNameAPAP = "MassQA_Particle3_2";
*/

    TString InputFileName;
    TString DirectoryName;
    TString ListName1;
    TString ListName2;
    TString ListNamePP;
    TString ListNameAPAP;
    TString HistoNamePP;
    TString HistoNameAPAP;
    //used to get the number of single lambdas
    TString ListNameV0_1;
    TString ListNameV0_2;
    TString HistoName_V0IMS;
    TString ListNameAV0_1;
    TString ListNameAV0_2;
    TString HistoName_AV0IMS;
    InputFileName = "/mnt/Ubuntu_Data/CernBox/Sync/CatsFiles_Dimi/pLambda/AnalysisResults_MEpurity.root";
    DirectoryName = "HMDimiResultsQAS40";
    ListName1 = "HMDimiResultsQAS40";
    ListName2 = "PairQA";
    ListNamePP = "QA_Particle0_Particle2";
    ListNameAPAP = "QA_Particle1_Particle3";
    if(MEflag==0){
        HistoNamePP = "MassQA_Particle2_2";
        HistoNameAPAP = "MassQA_Particle3_2";
    }
    else{
        HistoNamePP = "MEMassQA_Particle2_2";
        HistoNameAPAP = "MEMassQA_Particle3_2";
    }
    ListNameV0_1 = "HMDimiv0CutsS40";
    ListNameV0_2 = "v0Cuts";
    HistoName_V0IMS = "InvMassPt";
    ListNameAV0_1 = "HMDimiAntiv0CutsS40";
    ListNameAV0_2 = "v0Cuts";
    HistoName_AV0IMS = "InvMassPt";

    TFile* fInput = new TFile(InputFileName,"read");
    printf("fInput = %p\n",fInput);
    //fInput->Print();
    TDirectoryFile* dInput = (TDirectoryFile*)(fInput->FindObjectAny(DirectoryName));
    printf("dInput = %p\n",dInput);
    //dInput->ls();
    TList* lInput1 = NULL;
    dInput->GetObject(ListName1,lInput1);
    printf("lInput1 = %p\n",lInput1);
    //lInput1->ls();
    TList* lInput2 = (TList*)lInput1->FindObject(ListName2);
    printf("lInput2 = %p\n",lInput2);
    TList* lInputPP = (TList*)lInput2->FindObject(ListNamePP);
    printf("lInputPP = %p\n",lInputPP);
    lInputPP->ls();
    TList* lInputAPAP = (TList*)lInput2->FindObject(ListNameAPAP);
    printf("lInputAPAP = %p\n",lInputAPAP);
    lInputAPAP->ls();
    TH2F* hPP_IMS_k = (TH2F*)lInputPP->FindObject(HistoNamePP);
    printf("hPP_IMS_k = %p\n",hPP_IMS_k);
    TH2F* hAPAP_IMS_k = (TH2F*)lInputAPAP->FindObject(HistoNameAPAP);
    printf("hAPAP_IMS_k = %p\n",hAPAP_IMS_k);

    dInput = NULL;
    lInput1 = NULL;
    lInput2 = NULL;
    dInput = (TDirectoryFile*)(fInput->FindObjectAny(ListNameV0_1));
    dInput->GetObject(ListNameV0_1,lInput1);
    lInput2 = (TList*)lInput1->FindObject(ListNameV0_2);
    TH1F* hV0IMS = (TH1F*)lInput2->FindObject(HistoName_V0IMS);
    const double NumV0 = hV0IMS->Integral();

    dInput = NULL;
    lInput1 = NULL;
    lInput2 = NULL;
    dInput = (TDirectoryFile*)(fInput->FindObjectAny(ListNameAV0_1));
    dInput->GetObject(ListNameAV0_1,lInput1);
    lInput2 = (TList*)lInput1->FindObject(ListNameAV0_2);
    TH1F* hAV0IMS = (TH1F*)lInput2->FindObject(HistoName_AV0IMS);
    const double NumAV0 = hAV0IMS->Integral();

    const double ME_yield_bias = (hPP_IMS_k->Integral()+hAPAP_IMS_k->Integral())/(NumV0+NumAV0);
    //const double ME_yield_bias = 1;
    printf("ME_yield_bias = %f\n",ME_yield_bias);

    hPP_IMS_k->Sumw2();
    hAPAP_IMS_k->Sumw2();

    hPP_IMS_k->Rebin2D(RebinIMS,RebinKstar);
    hAPAP_IMS_k->Rebin2D(RebinIMS,RebinKstar);

    DLM_Histo<float>* dlm_PP_IMS_k = Convert_TH2F_DlmHisto(hPP_IMS_k);
//return;
    DLM_Histo<float>* dlm_APAP_IMS_k = Convert_TH2F_DlmHisto(hAPAP_IMS_k);
    DLM_Histo<float>* dlm_IMS_k = new DLM_Histo<float>();

    dlm_PP_IMS_k->RescaleAxis(0,1000,false);
    dlm_PP_IMS_k->RescaleAxis(1,1000,false);
    dlm_APAP_IMS_k->RescaleAxis(0,1000,false);
    dlm_APAP_IMS_k->RescaleAxis(1,1000,false);

    dlm_IMS_k->Copy(*dlm_PP_IMS_k);
    *dlm_IMS_k += *dlm_APAP_IMS_k;

    TH2F* hdimi_PP_IMS_k;
    hdimi_PP_IMS_k = Convert_DlmHisto_TH2F(dlm_PP_IMS_k,"hdimi_PP_IMS_k");

    TH2F* hdimi_APAP_IMS_k;
    hdimi_APAP_IMS_k = Convert_DlmHisto_TH2F(dlm_APAP_IMS_k,"hdimi_APAP_IMS_k");

    TH2F* hdimi_IMS_k;
    hdimi_IMS_k = Convert_DlmHisto_TH2F(dlm_IMS_k,"hdimi_IMS_k");



    double DefaultValue;
    double Error;
    double RandomValue;

    //TGraphErrors* gPurity = new TGraphErrors [3];
    TH2F* hInput;
    TH1D* hIMS;

    TFile* fOutput = new TFile(OutputFileName,"recreate");

    float Momentum;
    float Purity;
    float Chi2;
    float NDF;

    float Purity_p2;
    float Chi2_p2;
    float NDF_p2;

    unsigned Scenario;
    TGraph gFit;
    gFit.SetName("gFit");
    TGraph gBgr;
    gBgr.SetName("gBgr");
    TGraphErrors gData;
    gData.SetName("gData");

    TGraph gFit_p2;
    gFit_p2.SetName("gFit_p2");
    TGraph gBgr_p2;
    gBgr_p2.SetName("gBgr_p2");

    TTree* PurityTree = new TTree("PurityTree","PurityTree");
    PurityTree->Branch("Scenario", &Scenario, "Scenario/i");
    PurityTree->Branch("gData","TGraphErrors",&gData,32000,0);//
    PurityTree->Branch("gFit","TGraph",&gFit,32000,0);//
    PurityTree->Branch("gBgr","TGraph",&gBgr,32000,0);//
    PurityTree->Branch("Momentum", &Momentum, "Momentum/F");//
    PurityTree->Branch("Purity", &Purity, "Purity/F");//
    PurityTree->Branch("Chi2", &Chi2, "Chi2/F");//
    PurityTree->Branch("NDF", &NDF, "NDF/F");//

    PurityTree->Branch("gFit_p2","TGraph",&gFit_p2,32000,0);//
    PurityTree->Branch("gBgr_p2","TGraph",&gBgr_p2,32000,0);//
    PurityTree->Branch("Purity_p2", &Purity_p2, "Purity/F");//
    PurityTree->Branch("Chi2_p2", &Chi2_p2, "Chi2/F");//
    PurityTree->Branch("NDF_p2", &NDF_p2, "NDF/F");//

    fInput->cd();

    //3 scenarios: PP, APAP, Full
    for(unsigned uSce=2; uSce<3; uSce++){
printf("uSce=%u\n",uSce);
        Scenario = uSce;
        fInput->cd();
        if(uSce==0) hInput = hdimi_PP_IMS_k;
        else if(uSce==1) hInput = hdimi_APAP_IMS_k;
        else hInput = hdimi_IMS_k;
        for(unsigned uBin=0; uBin<hInput->GetYaxis()->GetNbins(); uBin++){
            fInput->cd();
            Momentum = hInput->GetYaxis()->GetBinCenter(uBin+1);
            if(Momentum>kMax) continue;
printf(" uBin=%u k*=%.0f MeV\n",uBin,Momentum);
            for(unsigned uBoot=FirstBoot; uBoot<LastBoot; uBoot++){
                TRandom3 rangen(uBoot+1);
//printf("  uBoot=%u\n",uBoot);
                fInput->cd();
                hIMS = hInput->ProjectionX("hIMS",uBin+1,uBin+1);
                //adjust the uncertainty based on the discrepancy between the number of single particles and particles in the ME
                for(unsigned uBinIMS=0; uBinIMS<hIMS->GetNbinsX(); uBinIMS++){
                    if(!MEflag) break;
                    hIMS->SetBinError(uBinIMS+1,hIMS->GetBinError(uBinIMS+1)*sqrt(ME_yield_bias));
                    if(hIMS->GetBinCenter(uBinIMS+1)>1106&&hIMS->GetBinCenter(uBinIMS+1)<1126){
                        hIMS->SetBinError(uBinIMS+1,hIMS->GetBinContent(uBinIMS+1)*2);
                    }
                }
                const double dataIntegral = hIMS->Integral(hIMS->FindBin(pdgLambda-2.*imWindow),hIMS->FindBin(pdgLambda+2.*imWindow));
                hIMS->Scale(1.,"width");
                ///used for the PAG presentation
                //double IM_min = 1080;
                //double IM_max = hIMS->GetXaxis()->GetBinUpEdge(hIMS->GetNbinsX());
                double IM_min = 1092;
                double IM_max = 1140;

                for(unsigned uBinIMS=0; uBinIMS<hIMS->GetNbinsX(); uBinIMS++){
                    fInput->cd();
                    if(hIMS->GetBinCenter(uBinIMS+1)<IM_min||hIMS->GetBinCenter(uBinIMS+1)>IM_max){
                        gData.SetPoint(uBinIMS,hIMS->GetBinCenter(uBinIMS+1),0);
                        gData.SetPointError(uBinIMS,0,0);
                        continue;
                    }
                    if(uBoot){
                        DefaultValue = hIMS->GetBinContent(uBinIMS+1);
                        Error = hIMS->GetBinError(uBinIMS+1);
                        RandomValue = 0;
                        while(RandomValue<1){
                            RandomValue = rangen.Gaus(DefaultValue,Error);
                        }
                        hIMS->SetBinContent(uBinIMS+1,RandomValue);
                    }
                    gData.SetPoint(uBinIMS,hIMS->GetBinCenter(uBinIMS+1),hIMS->GetBinContent(uBinIMS+1));
                    gData.SetPointError(uBinIMS,0,hIMS->GetBinError(uBinIMS+1));
                }
/*
                const int NumKnots = 9;
                double* Nodes_x = new double [NumKnots];
                Nodes_x[0] = IM_min;
                Nodes_x[1] = 1090;
                Nodes_x[2] = 1098;
                Nodes_x[3] = 1106;
                Nodes_x[4] = 1124;
                Nodes_x[5] = 1128;
                Nodes_x[6] = 1132;
                Nodes_x[7] = 1142;
                Nodes_x[8] = IM_max;
*/
/*
                const int NumKnots = 6;
                double* Nodes_x = new double [NumKnots];
                Nodes_x[0] = IM_min;
                Nodes_x[1] = 1090;
                Nodes_x[2] = 1100;
                Nodes_x[3] = 1130;
                Nodes_x[4] = 1140;
                Nodes_x[5] = IM_max;
*/
/*
                const int NumKnots = 9;
                double* Nodes_x = new double [NumKnots];
                Nodes_x[0] = IM_min;
                Nodes_x[1] = 1090;
                Nodes_x[2] = 1097;
                Nodes_x[3] = 1104;
                Nodes_x[4] = pdgLambda;
                Nodes_x[5] = 1124;
                Nodes_x[6] = 1130;
                Nodes_x[7] = 1140;
                Nodes_x[8] = IM_max;
*/
/*
                const int NumKnots = 8;
                double* Nodes_x = new double [NumKnots];
                Nodes_x[0] = IM_min;
                Nodes_x[1] = 1090;
                Nodes_x[2] = 1097;
                Nodes_x[3] = 1104;
                Nodes_x[4] = 1124;
                Nodes_x[5] = 1130;
                Nodes_x[6] = 1140;
                Nodes_x[7] = IM_max;
*/
/*
                const int NumKnots = 6;
                double* Nodes_x = new double [NumKnots];
                Nodes_x[0] = IM_min;
                Nodes_x[1] = 1090;
                Nodes_x[2] = 1100;
                Nodes_x[3] = 1132;
                Nodes_x[4] = 1140;
                Nodes_x[5] = IM_max;
*/

/*
///the one used for the PAG presentation
                const int NumKnots = 9;
                double* Nodes_x = new double [NumKnots];
                Nodes_x[0] = IM_min;
                Nodes_x[1] = 1090;
                Nodes_x[2] = 1104;
                Nodes_x[3] = 1115;
                Nodes_x[4] = 1126;
                Nodes_x[5] = 1130;
                Nodes_x[6] = 1135;
                Nodes_x[7] = 1140;
                Nodes_x[8] = IM_max;

                const int NumKnots_s = 6;
                double* Nodes_sx = new double [NumKnots_s];
                Nodes_sx[0] = IM_min;
                Nodes_sx[1] = 1090;
                Nodes_sx[2] = 1100;
                Nodes_sx[3] = 1132;
                Nodes_sx[4] = 1140;
                Nodes_sx[5] = IM_max;
*/

//1084-1092-1098-1104-1110-1116-1122-1128-1134-1140-1146
//1084-1098-1110-1122-1134-1146
//1104-1128
                const int NumKnots = 6;
                double* Nodes_x = new double [NumKnots];
                Nodes_x[0] = IM_min;
                Nodes_x[1] = 1098;
                Nodes_x[2] = 1104;
                Nodes_x[3] = 1128;
                Nodes_x[4] = 1134;
                Nodes_x[5] = IM_max;

                const int NumKnots_s = 4;
                double* Nodes_sx = new double [NumKnots_s];
                Nodes_sx[0] = IM_min;
                Nodes_sx[1] = 1104;
                Nodes_sx[2] = 1128;
                Nodes_sx[3] = IM_max;

/*
                const int NumKnots_s = 8;
                double* Nodes_sx = new double [NumKnots_s];
                Nodes_sx[0] = IM_min;
                Nodes_sx[1] = 1090;
                Nodes_sx[2] = 1097;
                Nodes_sx[3] = 1104;
                Nodes_sx[4] = 1124;
                Nodes_sx[5] = 1130;
                Nodes_sx[6] = 1140;
                Nodes_sx[7] = IM_max;
*/

                TF1* fData;

                TF1* fData_def = new TF1("fData_def",Fit_IMS_Splines,IM_min,IM_max,6+3+NumKnots*2);
                fData_def->SetParameter(0,dataIntegral*0.5);
                fData_def->SetParLimits(0,0,dataIntegral);
                fData_def->SetParameter(1,pdgLambda);
                fData_def->SetParLimits(1,pdgLambda-1.,pdgLambda+1.);
                fData_def->SetParameter(2,1.);
                fData_def->SetParLimits(2,0,6);

                fData_def->SetParameter(3,dataIntegral*0.5);
                fData_def->SetParLimits(3,0,dataIntegral);
                fData_def->FixParameter(4,0.);
                fData_def->SetParameter(5,5.);
                fData_def->SetParLimits(5,0,6);

fData_def->FixParameter(0,0);
fData_def->FixParameter(1,1);
fData_def->FixParameter(2,1);
fData_def->FixParameter(3,0);
fData_def->FixParameter(4,1);
fData_def->FixParameter(5,1);

                fData_def->FixParameter(6,NumKnots);
                //derivative at the firs knot
                fData_def->SetParameter(7,-10.);
                fData_def->SetParLimits(7,-100000,100000);
                fData_def->SetParameter(8,-10.);
                fData_def->SetParLimits(8,-100000,100000);

                for(unsigned uKnot=0; uKnot<NumKnots; uKnot++){
                    fInput->cd();
                    double HistVal = hIMS->GetBinContent(hIMS->FindBin(Nodes_x[uKnot]));
                    double HistErr = hIMS->GetBinError(hIMS->FindBin(Nodes_x[uKnot]));
                    //if(MEflag==1) HistErr *= sqrt(ME_yield_bias);
                    fData_def->FixParameter(9+uKnot,Nodes_x[uKnot]);
/*
                    if(uKnot==3){
                        double LOW = hIMS->GetBinContent(hIMS->FindBin(Nodes_x[uKnot-1]));
                        double UP = hIMS->GetBinContent(hIMS->FindBin(Nodes_x[uKnot+1]));
                        if(LOW>UP){
                            UP = LOW;
                            LOW = hIMS->GetBinContent(hIMS->FindBin(Nodes_x[uKnot+1]));
                        }
                        fData_def->SetParameter(9+NumKnots+uKnot,(UP+LOW)*0.5);
                        fData_def->SetParLimits(9+NumKnots+uKnot,LOW,UP);
                    }
                    else{
                        fData_def->SetParameter(9+NumKnots+uKnot,HistVal);
                        fData_def->SetParLimits(9+NumKnots+uKnot,HistVal-3.*HistErr,HistVal+3.*HistErr);
                    }
*/
/*
                    if(Nodes_x[uKnot]>=1110&&Nodes_x[uKnot]<=1122){
                        double LOW = hIMS->GetBinContent(hIMS->FindBin(1106));
                        double UP = hIMS->GetBinContent(hIMS->FindBin(1126));
                        if(LOW>UP){
                            LOW = hIMS->GetBinContent(hIMS->FindBin(1126));
                            UP = hIMS->GetBinContent(hIMS->FindBin(1106));
                        }
                        fData_def->SetParameter(9+NumKnots+uKnot,(UP+LOW)*0.5);
                        fData_def->SetParLimits(9+NumKnots+uKnot,(UP+LOW)*0.45,(UP+LOW)*0.55);
                    }
*/
/*
                    if(uKnot==5){
                        double LOW = hIMS->GetBinContent(hIMS->FindBin(1106));
                        double UP = hIMS->GetBinContent(hIMS->FindBin(1126));
                          fData_def->FixParameter(9+NumKnots+uKnot,(UP+LOW)*0.5);
                    }
                    else{
                        fData_def->SetParameter(9+NumKnots+uKnot,HistVal);
                        fData_def->SetParLimits(9+NumKnots+uKnot,0,2.*HistVal);
                    }
*/
                    fData_def->SetParameter(9+NumKnots+uKnot,HistVal);
                    fData_def->SetParLimits(9+NumKnots+uKnot,0,HistVal*2);
                }
                fData_def->SetNpx(1024);

                TF1* fData_s = new TF1("fData_s",Fit_IMS_Splines,IM_min,IM_max,6+3+NumKnots_s*2);
                fData_s->SetParameter(0,dataIntegral*0.5);
                fData_s->SetParLimits(0,0,dataIntegral);
                fData_s->SetParameter(1,pdgLambda);
                fData_s->SetParLimits(1,pdgLambda-1.,pdgLambda+1.);
                fData_s->SetParameter(2,1.);
                fData_s->SetParLimits(2,0,6);

                fData_s->SetParameter(3,dataIntegral*0.5);
                fData_s->SetParLimits(3,0,dataIntegral);
                fData_s->FixParameter(4,0.);
                fData_s->SetParameter(5,5.);
                fData_s->SetParLimits(5,0,6);

fData_s->FixParameter(0,0);
fData_s->FixParameter(1,1);
fData_s->FixParameter(2,1);
fData_s->FixParameter(3,0);
fData_s->FixParameter(4,1);
fData_s->FixParameter(5,1);

                fData_s->FixParameter(6,NumKnots_s);
                //derivative at the firs knot
                fData_s->SetParameter(7,-10.);
                fData_s->SetParLimits(7,-100000,100000);
                fData_s->SetParameter(8,-10.);
                fData_s->SetParLimits(8,-100000,100000);

                for(unsigned uKnot=0; uKnot<NumKnots_s; uKnot++){
                    fInput->cd();
                    double HistVal = hIMS->GetBinContent(hIMS->FindBin(Nodes_sx[uKnot]));
                    fData_s->FixParameter(9+uKnot,Nodes_sx[uKnot]);
                    fData_s->SetParameter(9+NumKnots_s+uKnot,HistVal);
                    fData_s->SetParLimits(9+NumKnots_s+uKnot,0,HistVal*2);
                }
                fData_s->SetNpx(1024);

                TF1* fit_Pol2_DG = new TF1("fit_Pol2_DG",Fit_IMS_Pol2,pdgLambda-25,pdgLambda+25,9);
                fit_Pol2_DG->SetParameter(0,dataIntegral*0.5);
                fit_Pol2_DG->SetParLimits(0,0,dataIntegral);
                fit_Pol2_DG->SetParameter(1,pdgLambda);
                fit_Pol2_DG->SetParLimits(1,pdgLambda-1.,pdgLambda+1.);
                fit_Pol2_DG->SetParameter(2,1.);
                fit_Pol2_DG->SetParLimits(2,0,6);

                fit_Pol2_DG->SetParameter(3,dataIntegral*0.5);
                fit_Pol2_DG->SetParLimits(3,0,dataIntegral);
                fit_Pol2_DG->FixParameter(4,0.);
                fit_Pol2_DG->SetParameter(5,5.);
                fit_Pol2_DG->SetParLimits(5,0,6);

                fit_Pol2_DG->SetParameter(6,0.);
                fit_Pol2_DG->SetParameter(7,0.);
                fit_Pol2_DG->SetParameter(8,0.);

                //hIMS->Fit(fData_def,"Q, S, N, R, M");
                hIMS->Fit(fData_s,"Q, S, N, R, M");
                hIMS->Fit(fit_Pol2_DG,"Q, S, N, R, M");
                Chi2 = fData_s->GetChisquare();
                NDF = fData_s->GetNDF();
                Chi2_p2 = fit_Pol2_DG->GetChisquare();
                NDF_p2 = fit_Pol2_DG->GetNDF();
                fData = fData_s;

                if(Chi2/NDF>2.0){
                    hIMS->Fit(fData_def,"S, N, R, M");
                    Chi2 = fData_def->GetChisquare();
                    NDF = fData_def->GetNDF();
                    fData = fData_def;
                    //printf("    refitting\n");
                }
                fData->SetName("fData");

                TF1* fSig = new TF1("fSig",Fit_IMS_Signal,IM_min,IM_max,6);
                for(unsigned uPar=0; uPar<6; uPar++) fSig->FixParameter(uPar,fData->GetParameter(uPar));

                TF1* fSigPol2 = new TF1("fSigPol2",Fit_IMS_Signal,pdgLambda-25,pdgLambda+25,6);
                for(unsigned uPar=0; uPar<6; uPar++) fSigPol2->FixParameter(uPar,fit_Pol2_DG->GetParameter(uPar));

                TF1* fPol2 = new TF1("fPol2","[0]+x*[1]+x*x*[2]",pdgLambda-25,pdgLambda+25);
                for(unsigned uPar=0; uPar<3; uPar++) fPol2->FixParameter(uPar,fit_Pol2_DG->GetParameter(6+uPar));

                for(unsigned uBinIMS=0; uBinIMS<hIMS->GetNbinsX(); uBinIMS++){
                    fInput->cd();
                    if(hIMS->GetBinCenter(uBinIMS+1)<IM_min||hIMS->GetBinCenter(uBinIMS+1)>IM_max){
                        gFit.SetPoint(uBinIMS,hIMS->GetBinCenter(uBinIMS+1),0);
                        gBgr.SetPoint(uBinIMS,hIMS->GetBinCenter(uBinIMS+1),0);
                        gFit_p2.SetPoint(uBinIMS,hIMS->GetBinCenter(uBinIMS+1),0);
                        gBgr_p2.SetPoint(uBinIMS,hIMS->GetBinCenter(uBinIMS+1),0);
                        continue;
                    }
                    gFit.SetPoint(uBinIMS,hIMS->GetBinCenter(uBinIMS+1),fData->Eval(hIMS->GetBinCenter(uBinIMS+1)));
                    gBgr.SetPoint(uBinIMS,hIMS->GetBinCenter(uBinIMS+1),fData->Eval(hIMS->GetBinCenter(uBinIMS+1))-fSig->Eval(hIMS->GetBinCenter(uBinIMS+1)));
                    gFit_p2.SetPoint(uBinIMS,hIMS->GetBinCenter(uBinIMS+1),fit_Pol2_DG->Eval(hIMS->GetBinCenter(uBinIMS+1)));
                    gBgr_p2.SetPoint(uBinIMS,hIMS->GetBinCenter(uBinIMS+1),fit_Pol2_DG->Eval(hIMS->GetBinCenter(uBinIMS+1))-fSigPol2->Eval(hIMS->GetBinCenter(uBinIMS+1)));
                }

                double Signal = fSig->Integral(pdgLambda-imWindow,pdgLambda+imWindow);
                double Total = fData->Integral(pdgLambda-imWindow,pdgLambda+imWindow);

                double SignalPol2 = fSigPol2->Integral(pdgLambda-imWindow,pdgLambda+imWindow);
                double TotalPol2 = fit_Pol2_DG->Integral(pdgLambda-imWindow,pdgLambda+imWindow);


                Purity = Signal/Total;
                Purity_p2 = SignalPol2/TotalPol2;
                fOutput->cd();
                PurityTree->Fill();

//printf("  Chi2/NDF = %.0f/%.0f = %.2f\n",Chi2,NDF,Chi2/NDF);
//printf("  Purity = %.2f\n",Purity*100.);
//printf("   For pol2:\n");
//printf("    Chi2/NDF = %.0f/%i = %.2f\n",fit_Pol2_DG->GetChisquare(),fit_Pol2_DG->GetNDF(),fit_Pol2_DG->GetChisquare()/double(fit_Pol2_DG->GetNDF()));
//printf("    Purity = %.2f\n",SignalPol2/TotalPol2*100.);

                hIMS->Write();
                fData->Write();
                fSig->Write();
                fit_Pol2_DG->Write();
                fSigPol2->Write();
                fPol2->Write();
                gFit.Write();
                gBgr.Write();

                fInput->cd();

                delete fData_def;
                delete fData_s;
                delete fSig;
                delete hIMS;
                delete fit_Pol2_DG;
                delete fPol2;
                delete fSigPol2;
            }
        }
    }

    fOutput->cd();
    PurityTree->Write();


    delete dlm_PP_IMS_k;
    delete dlm_APAP_IMS_k;
    delete dlm_IMS_k;
    delete hdimi_PP_IMS_k;
    delete hdimi_APAP_IMS_k;
    delete hdimi_IMS_k;
    //delete [] gPurity;
    delete fInput;
    delete fOutput;
}


//do a prefit
//MEflag 0 is se, 1 is me
//SIGflag 1-2-3 is single/double/triple Gaussian
void Purity_vs_kstar_2(const int& MEflag, const int& SIGflag, const unsigned& FirstBoot, const unsigned& LastBoot){
//PurityTree->Draw("Purity:Momentum>>h2(20,0,480,64,0.94, 0.98)", "Scenario==2", "COLZ")
    //const double ME_yield_bias = 17;
    const double kMax = 6000;
    //const double kMax = 960;
    //const double kMax = 480;
    //const unsigned BootstrapIter = 128;
    const double pdgLambda = 1115.683;
    const double imWindow = 4.0;

    const int RebinKstar = 6*2;
    //const int RebinIMS = MEflag?4/4:8;
    const int RebinIMS = 4;

    TString OutputFileName = TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/"
                            "pLambda_1/Purity_vs_kstar/Test/fWeightedPurity%s_G%i_%u_%u.root",
                            MEflag?"ME":"SE",SIGflag,FirstBoot,LastBoot);


/*
    TString InputFileName = "/mnt/Ubuntu_Data/CernBox/Sync/CatsFiles_Dimi/pLambda/AnalysisResults_Purity01062020.root";
    TString DirectoryName = "HMDimiResultsQA0";
    TString ListName1 = "HMDimiResultsQA0";
    TString ListName2 = "PairQA";
    TString ListNamePP = "QA_Particle0_Particle2";
    TString ListNameAPAP = "QA_Particle1_Particle3";
    TString HistoNamePP = "MassQA_Particle2_2";
    TString HistoNameAPAP = "MassQA_Particle3_2";
*/

    TString InputFileName;
    TString DirectoryName;
    TString ListName1;
    TString ListName2;
    TString ListNamePP;
    TString ListNameAPAP;
    TString HistoNamePP;
    TString HistoNameAPAP;
    //used to get the number of single lambdas
    TString ListNameV0_1;
    TString ListNameV0_2;
    TString HistoName_V0IMS;
    TString ListNameAV0_1;
    TString ListNameAV0_2;
    TString HistoName_AV0IMS;
    InputFileName = "/mnt/Ubuntu_Data/CernBox/Sync/CatsFiles_Dimi/pLambda/AnalysisResults_MEpurity.root";
    DirectoryName = "HMDimiResultsQAS40";
    ListName1 = "HMDimiResultsQAS40";
    ListName2 = "PairQA";
    ListNamePP = "QA_Particle0_Particle2";
    ListNameAPAP = "QA_Particle1_Particle3";
    if(MEflag==0){
        HistoNamePP = "MassQA_Particle2_2";
        HistoNameAPAP = "MassQA_Particle3_2";
    }
    else{
        HistoNamePP = "MEMassQA_Particle2_2";
        HistoNameAPAP = "MEMassQA_Particle3_2";
    }
    ListNameV0_1 = "HMDimiv0CutsS40";
    ListNameV0_2 = "v0Cuts";
    HistoName_V0IMS = "InvMassPt";
    ListNameAV0_1 = "HMDimiAntiv0CutsS40";
    ListNameAV0_2 = "v0Cuts";
    HistoName_AV0IMS = "InvMassPt";

    TFile* fInput = new TFile(InputFileName,"read");
    printf("fInput = %p\n",fInput);
    //fInput->Print();
    TDirectoryFile* dInput = (TDirectoryFile*)(fInput->FindObjectAny(DirectoryName));
    printf("dInput = %p\n",dInput);
    //dInput->ls();
    TList* lInput1 = NULL;
    dInput->GetObject(ListName1,lInput1);
    printf("lInput1 = %p\n",lInput1);
    //lInput1->ls();
    TList* lInput2 = (TList*)lInput1->FindObject(ListName2);
    printf("lInput2 = %p\n",lInput2);
    TList* lInputPP = (TList*)lInput2->FindObject(ListNamePP);
    printf("lInputPP = %p\n",lInputPP);
    lInputPP->ls();
    TList* lInputAPAP = (TList*)lInput2->FindObject(ListNameAPAP);
    printf("lInputAPAP = %p\n",lInputAPAP);
    lInputAPAP->ls();
    TH2F* hPP_IMS_k = (TH2F*)lInputPP->FindObject(HistoNamePP);
    printf("hPP_IMS_k = %p\n",hPP_IMS_k);
    TH2F* hAPAP_IMS_k = (TH2F*)lInputAPAP->FindObject(HistoNameAPAP);
    printf("hAPAP_IMS_k = %p\n",hAPAP_IMS_k);

    dInput = NULL;
    lInput1 = NULL;
    lInput2 = NULL;
    dInput = (TDirectoryFile*)(fInput->FindObjectAny(ListNameV0_1));
    dInput->GetObject(ListNameV0_1,lInput1);
    lInput2 = (TList*)lInput1->FindObject(ListNameV0_2);
    TH1F* hV0IMS = (TH1F*)lInput2->FindObject(HistoName_V0IMS);
    const double NumV0 = hV0IMS->Integral();

    dInput = NULL;
    lInput1 = NULL;
    lInput2 = NULL;
    dInput = (TDirectoryFile*)(fInput->FindObjectAny(ListNameAV0_1));
    dInput->GetObject(ListNameAV0_1,lInput1);
    lInput2 = (TList*)lInput1->FindObject(ListNameAV0_2);
    TH1F* hAV0IMS = (TH1F*)lInput2->FindObject(HistoName_AV0IMS);
    const double NumAV0 = hAV0IMS->Integral();

    //const double ME_yield_bias = (hPP_IMS_k->Integral()+hAPAP_IMS_k->Integral())/(NumV0+NumAV0);
const double ME_yield_bias = 1;
    printf("ME_yield_bias = %f\n",ME_yield_bias);

    hPP_IMS_k->Sumw2();
    hAPAP_IMS_k->Sumw2();

    hPP_IMS_k->Rebin2D(RebinIMS,RebinKstar);
    hAPAP_IMS_k->Rebin2D(RebinIMS,RebinKstar);

    DLM_Histo<float>* dlm_PP_IMS_k = Convert_TH2F_DlmHisto(hPP_IMS_k);
//return;
    DLM_Histo<float>* dlm_APAP_IMS_k = Convert_TH2F_DlmHisto(hAPAP_IMS_k);
    DLM_Histo<float>* dlm_IMS_k = new DLM_Histo<float>();

    dlm_PP_IMS_k->RescaleAxis(0,1000,false);
    dlm_PP_IMS_k->RescaleAxis(1,1000,false);
    dlm_APAP_IMS_k->RescaleAxis(0,1000,false);
    dlm_APAP_IMS_k->RescaleAxis(1,1000,false);

    dlm_IMS_k->Copy(*dlm_PP_IMS_k);
    *dlm_IMS_k += *dlm_APAP_IMS_k;

    TH2F* hdimi_PP_IMS_k;
    hdimi_PP_IMS_k = Convert_DlmHisto_TH2F(dlm_PP_IMS_k,"hdimi_PP_IMS_k");

    TH2F* hdimi_APAP_IMS_k;
    hdimi_APAP_IMS_k = Convert_DlmHisto_TH2F(dlm_APAP_IMS_k,"hdimi_APAP_IMS_k");

    TH2F* hdimi_IMS_k;
    hdimi_IMS_k = Convert_DlmHisto_TH2F(dlm_IMS_k,"hdimi_IMS_k");



    double DefaultValue;
    double Error;
    double RandomValue;

    //TGraphErrors* gPurity = new TGraphErrors [3];
    TH2F* hInput;
    TH1D* hIMS;

    TFile* fOutput = new TFile(OutputFileName,"recreate");

    float Momentum;
    float Purity;
    float Chi2;
    float NDF;
    //the fraction of SBleft or right
    //in the mass window of the lambda
    float SBleft;
    float SBleftOut;

    float Purity_p2;
    float Chi2_p2;
    float NDF_p2;

    unsigned Scenario;
    TGraph gFit;
    gFit.SetName("gFit");
    TGraph gBgr;
    gBgr.SetName("gBgr");
    TGraphErrors gData;
    gData.SetName("gData");

    TGraph gFit_p2;
    gFit_p2.SetName("gFit_p2");
    TGraph gBgr_p2;
    gBgr_p2.SetName("gBgr_p2");

    TTree* PurityTree = new TTree("PurityTree","PurityTree");
    PurityTree->Branch("Scenario", &Scenario, "Scenario/i");
    PurityTree->Branch("gData","TGraphErrors",&gData,32000,0);//
    PurityTree->Branch("gFit","TGraph",&gFit,32000,0);//
    PurityTree->Branch("gBgr","TGraph",&gBgr,32000,0);//
    PurityTree->Branch("Momentum", &Momentum, "Momentum/F");//
    PurityTree->Branch("Purity", &Purity, "Purity/F");//
    PurityTree->Branch("SBleft", &SBleft, "SBleft/F");//
    PurityTree->Branch("SBleftOut", &SBleftOut, "SBleftOut/F");//
    PurityTree->Branch("Chi2", &Chi2, "Chi2/F");//
    PurityTree->Branch("NDF", &NDF, "NDF/F");//

    PurityTree->Branch("gFit_p2","TGraph",&gFit_p2,32000,0);//
    PurityTree->Branch("gBgr_p2","TGraph",&gBgr_p2,32000,0);//
    PurityTree->Branch("Purity_p2", &Purity_p2, "Purity/F");//
    PurityTree->Branch("Chi2_p2", &Chi2_p2, "Chi2/F");//
    PurityTree->Branch("NDF_p2", &NDF_p2, "NDF/F");//


    const unsigned NumMomBins = TMath::Nint(kMax/hdimi_IMS_k->GetYaxis()->GetBinWidth(1));
    const double MinPur = 0.9;
    const double MaxPur = 1.0;
    const unsigned NumPurBins = 1024;
    const double MinSBl = 0.0;
    const double MaxSBl = 1.0;
    const unsigned NumSBlBins = 1024;
    TH2F* hPurMom = new TH2F("hPurMom","hPurMom",NumMomBins,0,kMax,NumPurBins,MinPur,MaxPur);
    //the fraction of left SB taken by integrating half of the peak region
    TH2F* hSBlMom = new TH2F("hSBlMom","hSBlMom",NumMomBins,0,kMax,NumSBlBins,MinSBl,MaxSBl);
    //the fraction of left SB by integrating just outside the peak region (4-8 MeV)
    TH2F* hSBlOutMom = new TH2F("hSBlOutMom","hSBlOutMom",NumMomBins,0,kMax,NumSBlBins,MinSBl,MaxSBl);

    TH1F* hWeightedPurity = new TH1F("hWeightedPurity","hWeightedPurity",NumPurBins,MinPur,MaxPur);

    fInput->cd();

    //3 scenarios: PP, APAP, Full
    for(unsigned uSce=2; uSce<3; uSce++){
printf("uSce=%u\n",uSce);
        Scenario = uSce;
        fInput->cd();
        if(uSce==0) hInput = hdimi_PP_IMS_k;
        else if(uSce==1) hInput = hdimi_APAP_IMS_k;
        else hInput = hdimi_IMS_k;

        double WeightedPurity = 0;
        double TotalWeight = 0;

        for(unsigned uBin=0; uBin<hInput->GetYaxis()->GetNbins(); uBin++){
            fInput->cd();
            Momentum = hInput->GetYaxis()->GetBinCenter(uBin+1);
            if(Momentum>kMax) continue;
printf(" uBin=%u k*=%.0f MeV\n",uBin,Momentum);
            for(unsigned uBoot=FirstBoot; uBoot<LastBoot; uBoot++){
                TRandom3 rangen(uBoot+1);
                fInput->cd();
                hIMS = hInput->ProjectionX("hIMS",uBin+1,uBin+1);
                //adjust the uncertainty based on the discrepancy between the number of single particles and particles in the ME
                for(unsigned uBinIMS=0; uBinIMS<hIMS->GetNbinsX(); uBinIMS++){
                    hIMS->SetBinError(uBinIMS+1,hIMS->GetBinError(uBinIMS+1)*sqrt(ME_yield_bias));
                }
                const double dataIntegral = hIMS->Integral(hIMS->FindBin(pdgLambda-2.*imWindow),hIMS->FindBin(pdgLambda+2.*imWindow));
                hIMS->Scale(1./dataIntegral,"width");
    ///delete
                TH1D* hIMS_Prefit = (TH1D*)hIMS->Clone("hIMS_Prefit");
                for(unsigned uBinIMS=0; uBinIMS<hIMS_Prefit->GetNbinsX(); uBinIMS++){
                    if(hIMS_Prefit->GetBinCenter(uBinIMS+1)>1104&&hIMS_Prefit->GetBinCenter(uBinIMS+1)<1128){
                        hIMS_Prefit->SetBinError(uBinIMS+1,hIMS_Prefit->GetBinContent(uBinIMS+1)*2);
                    }
                }

                ///used for the PAG presentation
                //double IM_min = 1080;
                //double IM_max = hIMS->GetXaxis()->GetBinUpEdge(hIMS->GetNbinsX());
                double IM_min = 1088;
                double IM_max = 1144;

                for(unsigned uBinIMS=0; uBinIMS<hIMS->GetNbinsX(); uBinIMS++){
                    fInput->cd();
                    if(hIMS->GetBinCenter(uBinIMS+1)<IM_min||hIMS->GetBinCenter(uBinIMS+1)>IM_max){
                        gData.SetPoint(uBinIMS,hIMS->GetBinCenter(uBinIMS+1),0);
                        gData.SetPointError(uBinIMS,0,0);
                        continue;
                    }
                    if(uBoot){
                        DefaultValue = hIMS->GetBinContent(uBinIMS+1);
                        Error = hIMS->GetBinError(uBinIMS+1);
                        RandomValue = 0;
                        while(RandomValue<=0&&DefaultValue>0){
                            RandomValue = rangen.Gaus(DefaultValue,Error);
                        }
                        hIMS->SetBinContent(uBinIMS+1,RandomValue);
                    }
                    gData.SetPoint(uBinIMS,hIMS->GetBinCenter(uBinIMS+1),hIMS->GetBinContent(uBinIMS+1));
                    gData.SetPointError(uBinIMS,0,hIMS->GetBinError(uBinIMS+1));
                }

//1084-1092-1098-1104-1110-1116-1122-1128-1134-1140-1146
//1084-1098-1110-1122-1134-1146
//1104-1128
//1088-1094-1100-1106-1126-1132-1138-1144
//1084-1092-1100-1108-1124-1134-1142-1150
                const int NumKnots = 8;
                double* Nodes_x = new double [NumKnots];
                //Nodes_x[0] = IM_min;
                //Nodes_x[1] = 1098;
                //Nodes_x[2] = 1104;
                //Nodes_x[3] = 1128;
                //Nodes_x[4] = 1134;
                //Nodes_x[5] = IM_max;
                Nodes_x[0] = IM_min;
                Nodes_x[1] = 1094;
                Nodes_x[2] = 1100;
                Nodes_x[3] = 1106;
                Nodes_x[4] = 1126;
                Nodes_x[5] = 1132;
                Nodes_x[6] = 1138;
                Nodes_x[7] = IM_max;

                TF1* fLinearLow = new TF1("fLinearLow","[0]+[1]*x",IM_min-3,IM_min+3);
                TF1* fLinearUp = new TF1("fLinearUp","[0]+[1]*x",IM_max-3,IM_max+3);
                hIMS->Fit(fLinearLow,"Q, S, N, R, M");
                hIMS->Fit(fLinearUp,"Q, S, N, R, M");

                TF1* fData;

                TF1* fData_def = new TF1("fData_def",Fit_IMS_Splines3,IM_min,IM_max,9+3+NumKnots*2);
                fData_def->FixParameter(0,0);
                fData_def->FixParameter(1,pdgLambda);
                fData_def->FixParameter(2,1.);
                fData_def->FixParameter(3,0);
                fData_def->FixParameter(4,0);
                fData_def->FixParameter(5,5.);
                fData_def->FixParameter(6,0);
                fData_def->FixParameter(7,0);
                fData_def->FixParameter(8,5.);

                fData_def->FixParameter(9,NumKnots);
                //derivative at the firs knot
                //fData_def->SetParameter(10,-0.01);
                //fData_def->SetParLimits(10,-10,10);
                //fData_def->SetParameter(11,-0.01);
                //fData_def->SetParLimits(11,-10,10);

                fData_def->FixParameter(10,fLinearLow->GetParameter(1));
                fData_def->FixParameter(11,fLinearUp->GetParameter(1));


                for(unsigned uKnot=0; uKnot<NumKnots; uKnot++){
                    fInput->cd();
                    double HistVal = hIMS->GetBinContent(hIMS->FindBin(Nodes_x[uKnot]));
                    double HistErr = hIMS->GetBinError(hIMS->FindBin(Nodes_x[uKnot]));
                    //if(MEflag==1) HistErr *= sqrt(ME_yield_bias);
                    fData_def->FixParameter(12+uKnot,Nodes_x[uKnot]);

                    fData_def->SetParameter(12+NumKnots+uKnot,HistVal);
                    fData_def->SetParLimits(12+NumKnots+uKnot,0,HistVal*2);
                }
                fData_def->SetNpx(1024);

                TF1* fit_Pol2_DG = new TF1("fit_Pol2_DG",Fit_IMS_Pol2,pdgLambda-25,pdgLambda+25,9);
                fit_Pol2_DG->SetParameter(0,1.*0.5);
                fit_Pol2_DG->SetParLimits(0,0,1.);
                fit_Pol2_DG->SetParameter(1,pdgLambda);
                fit_Pol2_DG->SetParLimits(1,pdgLambda-1.,pdgLambda+1.);
                fit_Pol2_DG->SetParameter(2,1.);
                fit_Pol2_DG->SetParLimits(2,0,6);

                fit_Pol2_DG->SetParameter(3,1.*0.5);
                fit_Pol2_DG->SetParLimits(3,0,1.);
                fit_Pol2_DG->FixParameter(4,0.);
                fit_Pol2_DG->SetParameter(5,5.);
                fit_Pol2_DG->SetParLimits(5,0,6);

                fit_Pol2_DG->SetParameter(6,0.);
                fit_Pol2_DG->SetParameter(7,0.);
                fit_Pol2_DG->SetParameter(8,0.);

                hIMS_Prefit->Fit(fData_def,"Q, S, N, R, M");


                double PAR_VAL;
                double PAR_ERR;
                for(unsigned uPar=0; uPar<9+3+NumKnots*2; uPar++){
                    PAR_VAL = fData_def->GetParameter(uPar);
                    fData_def->FixParameter(uPar,PAR_VAL);
                }


                double PeakFrac = 1./double(SIGflag);
                fData_def->SetParameter(0,1.0*PeakFrac);
                fData_def->SetParLimits(0,0,1.0);
                fData_def->SetParameter(1,pdgLambda);
                fData_def->SetParLimits(1,pdgLambda-1.,pdgLambda+1.);
                fData_def->SetParameter(2,1.5);
                fData_def->SetParLimits(2,0.75,8);

                fData_def->SetParameter(3,1.0*PeakFrac);
                fData_def->SetParLimits(3,0,1.0);
                fData_def->SetParameter(4,pdgLambda);
                //fData_def->SetParLimits(4,pdgLambda-1.,pdgLambda+1.);
                fData_def->FixParameter(4,0);
                fData_def->SetParameter(5,3.);
                fData_def->SetParLimits(5,0.75,8);

                fData_def->SetParameter(6,1.0*PeakFrac);
                fData_def->SetParLimits(6,0,1.0);
                fData_def->SetParameter(7,pdgLambda);
                //fData_def->SetParLimits(7,pdgLambda-1.,pdgLambda+1.);
                fData_def->FixParameter(7,0);
                fData_def->SetParameter(8,4.5);
                fData_def->SetParLimits(8,0.75,8);

                if(SIGflag<3){
                    fData_def->FixParameter(6,0);
                    fData_def->FixParameter(7,0);
                    fData_def->FixParameter(8,1);
                }
                if(SIGflag<2){
                    fData_def->FixParameter(3,0);
                    fData_def->FixParameter(4,0);
                    fData_def->FixParameter(5,1);
                }

            //fData_def->FixParameter(6,0);
            //fData_def->FixParameter(7,0);
            //fData_def->FixParameter(8,1);

                //we keep the first and last knot fixed from the prefit
                for(unsigned uKnot=2; uKnot<NumKnots-2; uKnot++){
                    fInput->cd();

                    PAR_VAL = fData_def->GetParameter(12+NumKnots+uKnot);
                    PAR_ERR = fabs(PAR_VAL)*0.2;
                    fData_def->SetParameter(12+NumKnots+uKnot,PAR_VAL);
                    fData_def->SetParLimits(12+NumKnots+uKnot,PAR_VAL-PAR_ERR,PAR_VAL+PAR_ERR);
                }

                hIMS->Fit(fit_Pol2_DG,"Q, S, N, R, M");
                hIMS->Fit(fData_def,"Q, S, N, R, M");
                Chi2 = fData_def->GetChisquare();
                NDF = fData_def->GetNDF();
                fData = fData_def;

                fData->SetName("fData");

                TF1* fSig = new TF1("fSig",Fit_IMS_Signal3,IM_min,IM_max,9);
                for(unsigned uPar=0; uPar<9; uPar++) fSig->FixParameter(uPar,fData->GetParameter(uPar));

                TF1* fSigPol2 = new TF1("fSigPol2",Fit_IMS_Signal,pdgLambda-25,pdgLambda+25,6);
                for(unsigned uPar=0; uPar<6; uPar++) fSigPol2->FixParameter(uPar,fit_Pol2_DG->GetParameter(uPar));

                TF1* fPol2 = new TF1("fPol2","[0]+x*[1]+x*x*[2]",pdgLambda-25,pdgLambda+25);
                for(unsigned uPar=0; uPar<3; uPar++) fPol2->FixParameter(uPar,fit_Pol2_DG->GetParameter(6+uPar));

                for(unsigned uBinIMS=0; uBinIMS<hIMS->GetNbinsX(); uBinIMS++){
                    fInput->cd();
                    if(hIMS->GetBinCenter(uBinIMS+1)<IM_min||hIMS->GetBinCenter(uBinIMS+1)>IM_max){
                        gFit.SetPoint(uBinIMS,hIMS->GetBinCenter(uBinIMS+1),0);
                        gBgr.SetPoint(uBinIMS,hIMS->GetBinCenter(uBinIMS+1),0);
                        gFit_p2.SetPoint(uBinIMS,hIMS->GetBinCenter(uBinIMS+1),0);
                        gBgr_p2.SetPoint(uBinIMS,hIMS->GetBinCenter(uBinIMS+1),0);
                        continue;
                    }
                    gFit.SetPoint(uBinIMS,hIMS->GetBinCenter(uBinIMS+1),fData->Eval(hIMS->GetBinCenter(uBinIMS+1)));
                    gBgr.SetPoint(uBinIMS,hIMS->GetBinCenter(uBinIMS+1),fData->Eval(hIMS->GetBinCenter(uBinIMS+1))-fSig->Eval(hIMS->GetBinCenter(uBinIMS+1)));
                    gFit_p2.SetPoint(uBinIMS,hIMS->GetBinCenter(uBinIMS+1),fit_Pol2_DG->Eval(hIMS->GetBinCenter(uBinIMS+1)));
                    gBgr_p2.SetPoint(uBinIMS,hIMS->GetBinCenter(uBinIMS+1),fit_Pol2_DG->Eval(hIMS->GetBinCenter(uBinIMS+1))-fSigPol2->Eval(hIMS->GetBinCenter(uBinIMS+1)));
                }

                double Signal = fSig->Integral(pdgLambda-imWindow,pdgLambda+imWindow);
                double SignalLeft = fSig->Integral(pdgLambda-imWindow,pdgLambda);
                double SignalLeftOut = fSig->Integral(pdgLambda-2.*imWindow,pdgLambda-imWindow);
                double SignalRightOut = fSig->Integral(pdgLambda+imWindow,pdgLambda+2.*imWindow);
                double Total = fData->Integral(pdgLambda-imWindow,pdgLambda+imWindow);
                double TotalLeft = fData->Integral(pdgLambda-imWindow,pdgLambda);
                double TotalLeftOut = fData->Integral(pdgLambda-2.*imWindow,pdgLambda-imWindow);
                double TotalRightOut = fData->Integral(pdgLambda+imWindow,pdgLambda+2.*imWindow);

                double SignalPol2 = fSigPol2->Integral(pdgLambda-imWindow,pdgLambda+imWindow);
                double TotalPol2 = fit_Pol2_DG->Integral(pdgLambda-imWindow,pdgLambda+imWindow);

                Purity = Signal/Total;
                Purity_p2 = SignalPol2/TotalPol2;
                SBleft = (TotalLeft-SignalLeft)/(Total-Signal);
                SBleftOut = (TotalLeftOut-SignalLeftOut)/(TotalLeftOut-SignalLeftOut+TotalRightOut-SignalRightOut);
                fOutput->cd();
                PurityTree->Fill();
                hPurMom->Fill(Momentum,Purity);
                hSBlMom->Fill(Momentum,SBleft);
                hSBlOutMom->Fill(Momentum,SBleftOut);

                if(uBoot==FirstBoot){
                    WeightedPurity += dataIntegral*Purity;
                    TotalWeight += dataIntegral;
                }


//printf("  Chi2/NDF = %.0f/%.0f = %.2f\n",Chi2,NDF,Chi2/NDF);
//printf("  Purity = %.2f\n",Purity*100.);
//printf("   For pol2:\n");
//printf("    Chi2/NDF = %.0f/%i = %.2f\n",fit_Pol2_DG->GetChisquare(),fit_Pol2_DG->GetNDF(),fit_Pol2_DG->GetChisquare()/double(fit_Pol2_DG->GetNDF()));
//printf("    Purity = %.2f\n",SignalPol2/TotalPol2*100.);

                hIMS->Write();
                hIMS_Prefit->Write();
                fData->Write();
                fSig->Write();
                fit_Pol2_DG->Write();
                fSigPol2->Write();
                fPol2->Write();
                gFit.Write();
                gBgr.Write();

                fInput->cd();

                delete fData_def;
                delete fSig;
                delete hIMS;
                delete hIMS_Prefit;
                delete fit_Pol2_DG;
                delete fPol2;
                delete fSigPol2;
                delete fLinearLow;
                delete fLinearUp;
            }
        }

        WeightedPurity /= TotalWeight;
        hWeightedPurity->Fill(WeightedPurity);

    }

    TH1F* hPurity = new TH1F("hPurity","hPurity",NumMomBins,0,kMax);
    TH1F* hSBl = new TH1F("hSBl","hSBl",NumMomBins,0,kMax);
    TH1F* hSBlOut = new TH1F("hSBlOut","hSBlOut",NumMomBins,0,kMax);
    for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
        TH1D* hProj = hPurMom->ProjectionY("hProj",uMomBin+1,uMomBin+1);
        double yMean = hProj->GetMean();
        double yStddev = hProj->GetStdDev();
        hPurity->SetBinContent(uMomBin+1,yMean);
        hPurity->SetBinError(uMomBin+1,yStddev);
        delete hProj;

        hProj = hSBlMom->ProjectionY("hProj",uMomBin+1,uMomBin+1);
        yMean = hProj->GetMean();
        yStddev = hProj->GetStdDev();
        hSBl->SetBinContent(uMomBin+1,yMean);
        hSBl->SetBinError(uMomBin+1,yStddev);
        delete hProj;

        hProj = hSBlOutMom->ProjectionY("hProj",uMomBin+1,uMomBin+1);
        yMean = hProj->GetMean();
        yStddev = hProj->GetStdDev();
        hSBlOut->SetBinContent(uMomBin+1,yMean);
        hSBlOut->SetBinError(uMomBin+1,yStddev);
        delete hProj;
    }

    fOutput->cd();
    PurityTree->Write();
    hPurity->Write();
    hPurMom->Write();
    hWeightedPurity->Write();
    hSBl->Write();
    hSBlMom->Write();
    hSBlOut->Write();
    hSBlOutMom->Write();

    delete dlm_PP_IMS_k;
    delete dlm_APAP_IMS_k;
    delete dlm_IMS_k;
    delete hdimi_PP_IMS_k;
    delete hdimi_APAP_IMS_k;
    delete hdimi_IMS_k;
    //delete [] gPurity;
    delete hPurity;
    delete hSBl;
    delete hSBlOut;
    delete hPurMom;
    delete hWeightedPurity;
    delete hSBlMom;
    delete hSBlOutMom;
    delete fInput;
    delete fOutput;
}
/*
void DimiPurityCorrected_pL(){

    TString NameSeFile;
    TString NameMeFile;
    TString NamePurHisto;

    TString NameDataFile;

    TH1F* hPurSe;
    TH1F* hPurMe;
    TH1F* hCorrCs;
    TH1F* hCorrCb;
    TH1F* hData;
    TH1F* hCs;
    TH1F* hCb;

}
*/

//0 is flat
//1 is chiral
//2 is esc
void FitMC_CompareToData_pL(int pSigma0_flag){

    TString InputFileNameData = "/mnt/Ubuntu_Data/CernBox/Sync/CatsFiles/ExpData/ALICE_pp_13TeV_HM/DimiJun20/Norm240_340/CompareDataMC/Ck_pL_Data.root";
    TString InputDirNameData = "Binning_12";
    TString InputHistoNameData = "hCkmult_SUM";
    TFile* fInputData = new TFile(InputFileNameData,"read");
    TDirectoryFile* fDirData = (TDirectoryFile*)(fInputData->FindObjectAny(InputDirNameData));
    printf("fDirData = %p\n",fDirData);
    //fDirData->ls();
    TH1D* hInputData = NULL;
    fDirData->GetObject(InputHistoNameData,hInputData);
    hInputData->SetName("hInputData");
    printf("hInputData = %p\n",hInputData);

    TString InputFileNameMc = "/mnt/Ubuntu_Data/CernBox/Sync/CatsFiles/ExpData/ALICE_pp_13TeV_HM/DimiJun20/Norm240_340/CompareDataMC/Ck_pL_Mc.root";
    TString InputDirNameMc = "Binning_60";
    TString InputHistoNameMc = "hCkmult_SUM";
    TFile* fInputMc = new TFile(InputFileNameMc,"read");
    TDirectoryFile* fDirMc = (TDirectoryFile*)(fInputMc->FindObjectAny(InputDirNameMc));
    printf("fDirMc = %p\n",fDirMc);
    //fDirMc->ls();
    TH1D* hInputMc = NULL;
    fDirMc->GetObject(InputHistoNameMc,hInputMc);
    hInputMc->SetName("hInputMc");
    printf("hInputMc = %p\n",hInputMc);

    TF1* fBaselineData = new TF1("fBaselineData",dimi_pL_May2020_FitterBl,0,500,5);

    //Chiral
    if(pSigma0_flag==1){
        fBaselineData->FixParameter(0,1.001);
        fBaselineData->FixParameter(1,0);
        fBaselineData->FixParameter(2,164.5);
        fBaselineData->FixParameter(3,3.206e-10);
        fBaselineData->FixParameter(4,-1e6);
    }
    //ESC16
    else if(pSigma0_flag==2){
        fBaselineData->FixParameter(0,0.9904);
//fBaselineData->FixParameter(0,1.000);
        fBaselineData->FixParameter(1,0);
        fBaselineData->FixParameter(2,2.864e-12);
        fBaselineData->FixParameter(3,2.504e-10);
        fBaselineData->FixParameter(4,-1e6);
    }
    //Flat
    else{
        fBaselineData->FixParameter(0,0.9913);
//fBaselineData->FixParameter(0,1.000);
        fBaselineData->FixParameter(1,0);
        fBaselineData->FixParameter(2,1.885e-6);
        fBaselineData->FixParameter(3,2.354e-10);
        fBaselineData->FixParameter(4,-1e6);
    }


    TF1* fBaselineMc = new TF1("fBaselineMc",dimi_pL_May2020_FitterBl,0,1000,5);
    fBaselineMc->SetParameter(0,1.001);
    fBaselineMc->SetParLimits(0,0.95,1.05);
    fBaselineMc->FixParameter(1,0);
    fBaselineMc->SetParameter(2,164.5);
    fBaselineMc->SetParLimits(2,100,700);
//fBaselineMc->FixParameter(2,1.885e-6);
    fBaselineMc->SetParameter(3,3.206e-10);
    fBaselineMc->SetParLimits(3,1e-11,1e-9);
    fBaselineMc->FixParameter(4,-1e6);

    TF1* fNormMc = new TF1("fNormMc","[0]",0,300);
    fNormMc->SetParameter(0,1.001);
    fNormMc->SetParLimits(0,0.95,1.05);

    hInputData->Fit(fBaselineData, "S, N, R, M");
    hInputMc->Fit(fBaselineMc, "S, N, R, M");
    hInputMc->Fit(fNormMc, "S, N, R, M");

    TFile* fOutput = new TFile(
    TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/FitMC_CompareToData_pL/fOutput_%i.root",pSigma0_flag)
                               ,"recreate");
    hInputData->Write();
    fBaselineData->Write();
    hInputMc->Write();
    fBaselineMc->Write();
    fNormMc->Write();

    delete fBaselineData;

    delete fBaselineMc;

    delete fInputData;
    delete fInputMc;

    delete fOutput;
}

void LambdaPars(TString DataSet){
    DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder("/mnt/Ubuntu_Data/CernBox/Sync/CatsFiles");

    double Purity_p[5];
    AnalysisObject.GetPurities_p(DataSet, 0, Purity_p);

    double Purity_L[5];
    AnalysisObject.GetPurities_L(DataSet, 0, Purity_L);

    double Frac_p[5];
    AnalysisObject.GetFractions_p(DataSet, 0, Frac_p);

    double Frac_L[5];
    AnalysisObject.GetFractions_L(DataSet, 0, Frac_L);

    double lambda_pars[5];
    AnalysisObject.SetUpLambdaPars_pL(DataSet, 0, 0, lambda_pars);

    printf("Protons:\n");
    printf(" pur p = %.2f\n", Purity_p[0]*100.);

    printf("Λs:\n");
    printf(" pur Λ = %.2f\n", Purity_L[0]*100.);

    printf("p-Λ:\n");
    printf(" p-Λ  -> p-Λ = %.2f\n", lambda_pars[0]*100.);
    printf(" p-Σ0 -> p-Λ = %.2f\n", lambda_pars[1]*100.);
    printf(" p-Ξ- -> p-Λ = %.2f\n", lambda_pars[2]*100.);
    printf(" p-Χ  -> p-Λ = %.2f\n", lambda_pars[3]*100.);
    printf(" fakes (m)   = %.2f\n", lambda_pars[4]*100.);
    printf(" m(Λ) any(p) = %.2f\n", Purity_L[0]*100.);
    printf(" m(p) t(Λ)   = %.2f\n", Purity_p[3]*Purity_L[0]*100.);
    printf(" m(p) m(Λ)   = %.2f\n", Purity_p[3]*Purity_L[4]*100.);


}



void EffectOfLambdaKaon(){
    const double MaxSize = 0.55/100.;

    const double Ref = -0.60;
    const double Imf = 0.51;
    const double EffR = 0.83;

    const unsigned NumMomBins = 30;
    const double kMin = 0;
    const double kMax = 360;

    DLM_Ck Ck_LambdaKaon(1,3,NumMomBins,kMin,kMax,ComplexLednicky_Singlet);
    Ck_LambdaKaon.SetSourcePar(0,1.3);
    Ck_LambdaKaon.SetPotPar(0,Ref);
    Ck_LambdaKaon.SetPotPar(1,Imf);
    Ck_LambdaKaon.SetPotPar(2,EffR);
    Ck_LambdaKaon.Update();
    DLM_CkDecomposition CkDec_LambdaKaon("LambdaKaon",0,Ck_LambdaKaon,NULL);
    CkDec_LambdaKaon.Update();


    DLM_Ck Dummy_pL(0,0,NumMomBins,kMin,kMax,Flat_Residual);
    Dummy_pL.Update();
    DLM_CkDecomposition CkDec_pL("pLambda",1,Dummy_pL,NULL);
    CkDec_pL.AddContribution(0,MaxSize,DLM_CkDecomposition::cFake,&CkDec_LambdaKaon);
    CkDec_pL.Update();

    TGraph g_KL;
    g_KL.SetName("g_KL");
    TGraph g_KL_pL;
    g_KL_pL.SetName("g_KL_pL");

    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        double MOM = Ck_LambdaKaon.GetBinCenter(0,uBin);
        g_KL.SetPoint(uBin,MOM,Ck_LambdaKaon.GetBinContent(uBin));
        g_KL_pL.SetPoint(uBin,MOM,CkDec_pL.EvalCk(MOM));
    }

    TFile fOutput("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/EffectOfLambdaKaon/fOutput.root","recreate");
    g_KL.Write();
    g_KL_pL.Write();
}

void FeedDown_pXi0_pL(const int& SEED, const int& NumIter){
    DLM_DecayMatrix DECM;
    DECM.SetFileName(TString::Format("/mnt/Ubuntu_Data/CernBox/Sync/CatsFiles/DecaySmear/pXi0_pL_RS%i.root",SEED));
    DECM.SetHistoName("pXi0_pL");
    DECM.SetBins(1024,0,1024);
    DECM.SetNumDaughters1(1);
    DECM.SetMotherMass1(Mass_p);
    DECM.SetNumDaughters2(2);
    DECM.SetDaughterMass2(0,Mass_L);
    DECM.SetDaughterMass2(1,Mass_pi0);
    DECM.SetMotherMass2(Mass_Xi0);
    DECM.SetMeanMomentum(0);
    DECM.SetMomentumSpread(350);
    DECM.Run(SEED,NumIter);
}

void EffectOfpXi0(TString DataSet){
    const unsigned NumMomBins = 30;
    const double kMin = 0;
    const double kMax = 360;

    CATS Kitty_pXi0;

    CATSparameters cPars(CATSparameters::tSource,1,true);
    cPars.SetParameter(0,1.0);
    Kitty_pXi0.SetAnaSource(GaussSource, cPars);
    Kitty_pXi0.SetUseAnalyticSource(true);
    Kitty_pXi0.SetAutoNormSource(false);

    Kitty_pXi0.SetMomBins(NumMomBins,kMin,kMax);
    Kitty_pXi0.SetExcludeFailedBins(false);
    Kitty_pXi0.SetQ1Q2(0);
    Kitty_pXi0.SetPdgId(2212, 3312);
    Kitty_pXi0.SetRedMass( (Mass_p*Mass_Xi0)/(Mass_p+Mass_Xi0) );

    //0: S=0 I=1
    //1: S=1 I=1
    Kitty_pXi0.SetNumChannels(2);
    Kitty_pXi0.SetNumPW(0,1);
    Kitty_pXi0.SetNumPW(1,1);
    Kitty_pXi0.SetSpin(0,0);
    Kitty_pXi0.SetSpin(1,1);
    Kitty_pXi0.SetChannelWeight(0, 1./4.);
    Kitty_pXi0.SetChannelWeight(1, 3./4.);

    CATSparameters cPars_pXim_HalI1S0(CATSparameters::tPotential,8,true);
    cPars_pXim_HalI1S0.SetParameter(0,pXim_HALQCDPaper2020);
    cPars_pXim_HalI1S0.SetParameter(1,12);
    cPars_pXim_HalI1S0.SetParameter(2,1);
    cPars_pXim_HalI1S0.SetParameter(3,1);
    cPars_pXim_HalI1S0.SetParameter(4,1);
    cPars_pXim_HalI1S0.SetParameter(5,0);
    cPars_pXim_HalI1S0.SetParameter(6,0);
    cPars_pXim_HalI1S0.SetParameter(7,0);

    CATSparameters cPars_pXim_HalI1S1(CATSparameters::tPotential,8,true);
    cPars_pXim_HalI1S1.SetParameter(0,pXim_HALQCDPaper2020);
    cPars_pXim_HalI1S1.SetParameter(1,12);
    cPars_pXim_HalI1S1.SetParameter(2,1);
    cPars_pXim_HalI1S1.SetParameter(3,1);
    cPars_pXim_HalI1S1.SetParameter(4,1);
    cPars_pXim_HalI1S1.SetParameter(5,1);
    cPars_pXim_HalI1S1.SetParameter(6,0);
    cPars_pXim_HalI1S1.SetParameter(7,1);

    Kitty_pXi0.SetShortRangePotential(0,0,fDlmPot,cPars_pXim_HalI1S0);
    Kitty_pXi0.SetShortRangePotential(1,0,fDlmPot,cPars_pXim_HalI1S1);

    DLM_Ck Ck_pXi0(Kitty_pXi0.GetNumSourcePars(),0,Kitty_pXi0);
    Ck_pXi0.SetSourcePar(0,1.0);
    Ck_pXi0.Update();

    DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder("/mnt/Ubuntu_Data/CernBox/Sync/CatsFiles");
    TH2F* hRes_pXi0 = AnalysisObject.GetResidualMatrix("pLambda","pXi0");
    double lam_pXim[5];
    AnalysisObject.SetUpLambdaPars_pXim(DataSet,0,0,lam_pXim);
    //printf("lam_pXim[0] = %f\n",lam_pXim[0]);
    //printf("lam_pXim[1] = %f\n",lam_pXim[1]);
    //printf("lam_pXim[2] = %f\n",lam_pXim[2]);
    //printf("lam_pXim[3] = %f\n",lam_pXim[3]);
    //printf("lam_pXim[4] = %f\n",lam_pXim[4]);

    DLM_CkDecomposition CkDec_pXi0("pXi0",2,Ck_pXi0,NULL);
    CkDec_pXi0.AddContribution(0,lam_pXim[1]+lam_pXim[2]+lam_pXim[3],DLM_CkDecomposition::cFeedDown);
    CkDec_pXi0.AddContribution(1,lam_pXim[4],DLM_CkDecomposition::cFake);
    CkDec_pXi0.Update();

    DLM_Ck Dummy_pL(0,0,NumMomBins,kMin,kMax,Flat_Residual);
    Dummy_pL.Update();
    DLM_CkDecomposition CkDec_pL("pLambda",1,Dummy_pL,NULL);
    CkDec_pL.AddContribution(0,lam_pXim[1]*2.,DLM_CkDecomposition::cFeedDown,&CkDec_pXi0,hRes_pXi0);
    CkDec_pL.Update();

    TGraph g_pXi0;
    g_pXi0.SetName("g_pXi0");
    TGraph g_pXi0_pL;
    g_pXi0_pL.SetName("g_pXi0_pL");

    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        double MOM = Ck_pXi0.GetBinCenter(0,uBin);
        g_pXi0.SetPoint(uBin,MOM,Ck_pXi0.GetBinContent(uBin));
        g_pXi0_pL.SetPoint(uBin,MOM,CkDec_pL.EvalCk(MOM));
    }

    TFile fOutput("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/EffectOfpXi0/fOutput.root","recreate");
    g_pXi0.Write();
    g_pXi0_pL.Write();

}

void CreateFineBinnedResoMatrix_pp(){
    const TString InputFileName = "/mnt/Ubuntu_Data/CernBox/Sync/CatsFiles/ExpData/ALICE_pp_13TeV_HM_MCnano/AnalysisResults.root";
    const TString OutputFileName = "/mnt/Ubuntu_Data/CernBox/Sync/CatsFiles/MomentumSmear/ALICE_pp_13TeV_MEpp.root";
    const TString DirectoryName = "HMBBarResultsQA8";
    const TString ListName1 = "HMBBarResultsQA8";
    const TString ListName2 = "PairQA";
    const TString ListNamePP = "QA_Particle0_Particle0";
    const TString ListNameAPAP = "QA_Particle1_Particle1";
    const TString HistoNamePP = "MomentumResolutionME_Particle0_Particle0";
    const TString HistoNameAPAP = "MomentumResolutionME_Particle1_Particle1";

    TFile* fInput = new TFile(InputFileName,"read");
    printf("fInput = %p\n",fInput);
    TDirectoryFile* dInput = (TDirectoryFile*)(fInput->FindObjectAny(DirectoryName));
    printf("dInput = %p\n",dInput);
    //dInput->ls();
    TList* lInput1 = NULL;
    dInput->GetObject(ListName1,lInput1);
    printf("lInput1 = %p\n",lInput1);
    //lInput1->ls();
    TList* lInput2 = (TList*)lInput1->FindObject(ListName2);
    printf("lInput2 = %p\n",lInput2);


    TList* lInputPP = (TList*)lInput2->FindObject(ListNamePP);
    printf("lInputPP = %p\n",lInputPP);
    lInputPP->ls();

    TList* lInputAPAP = (TList*)lInput2->FindObject(ListNameAPAP);
    printf("lInputAPAP = %p\n",lInputAPAP);
    lInputAPAP->ls();

    TH2F* hPP_RESO = (TH2F*)lInputPP->FindObject(HistoNamePP);
    printf("hPP_RESO = %p\n",hPP_RESO);
    TH2F* hAPAP_RESO = (TH2F*)lInputAPAP->FindObject(HistoNameAPAP);
    printf("hAPAP_RESO = %p\n",hAPAP_RESO);

    const unsigned NumMomBins = hPP_RESO->GetXaxis()->GetNbins();
    double kMin_MeV = hPP_RESO->GetXaxis()->GetBinLowEdge(1)*1000.;
    double kMax_MeV = hPP_RESO->GetXaxis()->GetBinUpEdge(NumMomBins)*1000.;

    TFile* f_RESO_pp_MeV = new TFile(OutputFileName,"recreate");

    TH2F* h_RESO_pp_MeV = new TH2F("h_RESO_pp_MeV","h_RESO_pp_MeV",
                                   NumMomBins,kMin_MeV,kMax_MeV,NumMomBins,kMin_MeV,kMax_MeV);

    double Val_PP, Val_APAP, Val;
    double Err_PP, Err_APAP, Err;

    for(unsigned uMomBinX=0; uMomBinX<NumMomBins; uMomBinX++){
        for(unsigned uMomBinY=0; uMomBinY<NumMomBins; uMomBinY++){
            Val_PP = hPP_RESO->GetBinContent(uMomBinX+1,uMomBinY+1);
            Val_APAP = hAPAP_RESO->GetBinContent(uMomBinX+1,uMomBinY+1);
            Err_PP = sqrt(Val_PP);
            Err_APAP = sqrt(Val_APAP);
            h_RESO_pp_MeV->SetBinContent(uMomBinX+1,uMomBinY+1,Val_PP+Val_APAP);
            h_RESO_pp_MeV->SetBinError(uMomBinX+1,uMomBinY+1,sqrt(Err_PP*Err_PP+Err_APAP*Err_APAP));
        }
    }

    h_RESO_pp_MeV->Write();

    delete h_RESO_pp_MeV;
    delete f_RESO_pp_MeV;
}

void CreateFineBinnedResoMatrix_pL(){
    const TString InputFileName = "/mnt/Ubuntu_Data/CernBox/Sync/CatsFiles/ExpData/ALICE_pp_13TeV_HM_MCnano/AnalysisResults.root";
    const TString OutputFileName = "/mnt/Ubuntu_Data/CernBox/Sync/CatsFiles/MomentumSmear/ALICE_pp_13TeV_MEpL.root";
    const TString DirectoryName = "HMBBarResultsQA8";
    const TString ListName1 = "HMBBarResultsQA8";
    const TString ListName2 = "PairQA";
    const TString ListNamePP = "QA_Particle0_Particle2";
    const TString ListNameAPAP = "QA_Particle1_Particle3";
    const TString HistoNamePP = "MomentumResolutionME_Particle0_Particle2";
    const TString HistoNameAPAP = "MomentumResolutionME_Particle1_Particle3";

    TFile* fInput = new TFile(InputFileName,"read");
    printf("fInput = %p\n",fInput);
    TDirectoryFile* dInput = (TDirectoryFile*)(fInput->FindObjectAny(DirectoryName));
    printf("dInput = %p\n",dInput);
    //dInput->ls();
    TList* lInput1 = NULL;
    dInput->GetObject(ListName1,lInput1);
    printf("lInput1 = %p\n",lInput1);
    //lInput1->ls();
    TList* lInput2 = (TList*)lInput1->FindObject(ListName2);
    printf("lInput2 = %p\n",lInput2);


    TList* lInputPP = (TList*)lInput2->FindObject(ListNamePP);
    printf("lInputPP = %p\n",lInputPP);
    lInputPP->ls();

    TList* lInputAPAP = (TList*)lInput2->FindObject(ListNameAPAP);
    printf("lInputAPAP = %p\n",lInputAPAP);
    lInputAPAP->ls();

    TH2F* hPP_RESO = (TH2F*)lInputPP->FindObject(HistoNamePP);
    printf("hPP_RESO = %p\n",hPP_RESO);
    TH2F* hAPAP_RESO = (TH2F*)lInputAPAP->FindObject(HistoNameAPAP);
    printf("hAPAP_RESO = %p\n",hAPAP_RESO);

    const unsigned NumMomBins = hPP_RESO->GetXaxis()->GetNbins();
    double kMin_MeV = hPP_RESO->GetXaxis()->GetBinLowEdge(1)*1000.;
    double kMax_MeV = hPP_RESO->GetXaxis()->GetBinUpEdge(NumMomBins)*1000.;

    TFile* f_RESO_pL_MeV = new TFile(OutputFileName,"recreate");

    TH2F* h_RESO_pL_MeV = new TH2F("h_RESO_pL_MeV","h_RESO_pL_MeV",
                                   NumMomBins,kMin_MeV,kMax_MeV,NumMomBins,kMin_MeV,kMax_MeV);

    double Val_PP, Val_APAP, Val;
    double Err_PP, Err_APAP, Err;

    for(unsigned uMomBinX=0; uMomBinX<NumMomBins; uMomBinX++){
        for(unsigned uMomBinY=0; uMomBinY<NumMomBins; uMomBinY++){
            Val_PP = hPP_RESO->GetBinContent(uMomBinX+1,uMomBinY+1);
            Val_APAP = hAPAP_RESO->GetBinContent(uMomBinX+1,uMomBinY+1);
            Err_PP = sqrt(Val_PP);
            Err_APAP = sqrt(Val_APAP);
            h_RESO_pL_MeV->SetBinContent(uMomBinX+1,uMomBinY+1,Val_PP+Val_APAP);
            h_RESO_pL_MeV->SetBinError(uMomBinX+1,uMomBinY+1,sqrt(Err_PP*Err_PP+Err_APAP*Err_APAP));
        }
    }

    h_RESO_pL_MeV->Write();

    delete h_RESO_pL_MeV;
    delete f_RESO_pL_MeV;
}


void pLambda_DummyCk_DifferentRadii(){

  //0 HM and OO, 1 is only OO
  const int PlotType = 0;
  const double lambdapar = 0.47;
  unsigned NumRadii = 2;
  double* Radii = new double [NumRadii];
  Radii[0] = 1.25;
  Radii[1] = 2.5;
  //Radii[2] = 3.0;
  //Radii[3] = 4.0;

  //how large should the new uncertainties be?
  double* UncertaintyScale = new double [NumRadii];
  UncertaintyScale[0] = 1.0;
  UncertaintyScale[1] = 0.7;
  //UncertaintyScale[2] = 1.0;
  //UncertaintyScale[3] = 1.0;
  //get the exp. correlation, use the stat uncertainties from there

  const unsigned NumModels = 2;
  int* ChiralType = new int [NumModels];
  ChiralType[0] = 1600;
  ChiralType[1] = 11600;

  const unsigned NumRatios = (NumModels*(NumModels-1))/2;

  unsigned NumMomBins = 30;
  double kMax = 360;

  TGraphErrors** gCk = new TGraphErrors* [NumModels];
  TGraphErrors** gCkFine = new TGraphErrors* [NumModels];
  for(unsigned uMod=0; uMod<NumModels; uMod++){
    gCk[uMod] = new TGraphErrors[NumRadii];
    gCkFine[uMod] = new TGraphErrors[NumRadii];
    for(unsigned uRad=0; uRad<NumRadii; uRad++){
      gCk[uMod][uRad].SetName(TString::Format("gCk_%.1f_%i",Radii[uRad],ChiralType[uMod]));
      gCkFine[uMod][uRad].SetName(TString::Format("gCkFine_%.1f_%i",Radii[uRad],ChiralType[uMod]));
    }
  }

  TGraphErrors** gCkRat = new TGraphErrors* [NumRadii];
  for(unsigned uRad=0; uRad<NumRadii; uRad++){
    gCkRat[uRad] = new TGraphErrors[NumRatios];
    unsigned uRat=0;
    for(unsigned uMod1=0; uMod1<NumModels; uMod1++){
      for(unsigned uMod2=uMod1+1; uMod2<NumModels; uMod2++){
        gCkRat[uRad][uRat].SetName(TString::Format("gCkRat_%.1f_%i_%i",Radii[uRad],ChiralType[uMod1],ChiralType[uMod2]));
        uRat++;
      }
    }
  }
  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());
  for(unsigned uMod=0; uMod<NumModels; uMod++){
    CATS AB_BIN;
    AB_BIN.SetMomBins(NumMomBins,0,kMax);
    CATS AB_pL;
    AB_pL.SetNotifications(CATS::nWarning);
    AB_pL.SetMomBins(NumMomBins*3,0,kMax);
    //AnalysisObject.SetUpCats_pL(AB_pL,"Chiral_Coupled_SPD","McGauss_ResoTM",ChiralType[uMod],202);//NLO_Coupled_S
    AnalysisObject.SetUpCats_pL(AB_pL,"Chiral_Coupled_SPD","Gauss",ChiralType[uMod],202);//NLO_Coupled_S
    const double CuspWeight = 0.4;//0.54
    AB_pL.SetChannelWeight(7,1./4.*CuspWeight);//1S0 SN(s) -> LN(s)
    AB_pL.SetChannelWeight(8,3./4.*CuspWeight);//3S1 SN(s) -> LN(s)
    AB_pL.SetChannelWeight(10,3./4.*CuspWeight);//3S1 SN(d) -> LN(s)
    AB_pL.SetChannelWeight(13,3./20.*CuspWeight);//3D1 SN(d) -> LN(d)
    AB_pL.SetChannelWeight(15,3./20.*CuspWeight);//3D1 SN(s) -> LN(d)
    for(unsigned uRad=0; uRad<NumRadii; uRad++){
      double Prog = double(uMod*NumRadii+uRad)/double(NumRadii*NumModels);
      unsigned NumStars = TMath::Nint(20.*Prog);
      unsigned NumDash = 20-NumStars;
      printf("\r\033[K Progress %.0f%% ",100.*Prog);
      for(unsigned ustar=0; ustar<NumStars; ustar++) printf("*");
      for(unsigned udash=0; udash<NumDash; udash++) printf("-");
      cout << flush;
      AB_pL.SetAnaSource(0,Radii[uRad]);
      AB_pL.KillTheCat();
      for(unsigned uMom=0; uMom<AB_pL.GetNumMomBins(); uMom++){
        gCkFine[uMod][uRad].SetPoint(uMom,AB_pL.GetMomentum(uMom),lambdapar*AB_pL.GetCorrFun(uMom)+1-lambdapar);
      }
      for(unsigned uMom=0; uMom<NumMomBins; uMom++){
        //if(uMom==23){
          double c1,c2,c3,avg;
          c1 = AB_pL.EvalCorrFun(AB_BIN.GetMomentum(uMom)-6.);
          c2 = AB_pL.EvalCorrFun(AB_BIN.GetMomentum(uMom));
          c3 = AB_pL.EvalCorrFun(AB_BIN.GetMomentum(uMom)+6.);
          avg = (c1+c2+c3)/3.;
          gCk[uMod][uRad].SetPoint(uMom,AB_BIN.GetMomentum(uMom),lambdapar*avg+1-lambdapar);
        //}
        //else
          //gCk[uMod][uRad].SetPoint(uMom,AB_BIN.GetMomentum(uMom),lambdapar*AB_pL.EvalCorrFun(AB_BIN.GetMomentum(uMom))+1-lambdapar);
      }
    }
  }

  for(unsigned uRad=0; uRad<NumRadii; uRad++){
    unsigned uRat=0;
    for(unsigned uMod1=0; uMod1<NumModels; uMod1++){
      for(unsigned uMod2=uMod1+1; uMod2<NumModels; uMod2++){
        double Val1,Val2,mom1,mom2;
        for(unsigned uMom=0; uMom<NumMomBins; uMom++){
          gCk[uMod1][uRad].GetPoint(uMom,mom1,Val1);
          gCk[uMod2][uRad].GetPoint(uMom,mom2,Val2);
          gCkRat[uRad][uRat].SetPoint(uMom,mom1,Val1/Val2);
        }
        uRat++;
      }
    }
  }

  printf("\r\033[K Progress %.0f%% ",100.);
  for(unsigned ustar=0; ustar<20; ustar++) printf("*");
  printf("\n");

  TString OutputFile = TString::Format("%s/pLambda/pLambda_DummyCk_DifferentRadii.root",GetFemtoOutputFolder());
  printf("Saving the output into %s\n",OutputFile.Data());
  TFile fOutput(OutputFile,"recreate");
  for(unsigned uRad=0; uRad<NumRadii; uRad++){
    for(unsigned uMod=0; uMod<NumModels; uMod++){
      gCk[uMod][uRad].Write();
    }
    unsigned uRat=0;
    for(unsigned uMod1=0; uMod1<NumModels; uMod1++){
      for(unsigned uMod2=uMod1+1; uMod2<NumModels; uMod2++){
        gCkRat[uRad][uRat].Write();
        uRat++;
      }
    }
  }



//PLOT

  TRandom3 rangen(7);

  double MOM;
  double CKVAL;
  double CKERR;
  int CKBIN;
  TH1F* hData_pL_Stat = AnalysisObject.GetAliceExpCorrFun("pp13TeV_HM_DimiJul20","pLambda","L53_SL4_SR6_P96_0",2,false,-1);
  fOutput.cd();
  TH1F* hProjection = (TH1F*)hData_pL_Stat->Clone("hProjection");
  TH1F* hGenHM = (TH1F*)hData_pL_Stat->Clone("hGenHM");
  //hData_pL_Stat->Scale(1.01);
  gStyle->SetCanvasPreferGL(1);
  SetStyle();

  hData_pL_Stat->Write();
  hProjection->Write();

  TGraphErrors** gDummyData = new TGraphErrors* [NumModels];
  for(unsigned uMod=0; uMod<NumModels; uMod++){
    gDummyData[uMod] = new TGraphErrors[NumRadii];
    for(unsigned uRad=0; uRad<NumRadii; uRad++){
      gDummyData[uMod][uRad].SetName(TString::Format("gDummyData_%.1f_%i",Radii[uRad],ChiralType[uMod]));
      for(unsigned uMom=0; uMom<NumMomBins; uMom++){
        gCk[uMod][uRad].GetPoint(uMom,MOM,CKVAL);
        //MOM = AB_pL.GetMomentum(uMom);
        CKBIN = hData_pL_Stat->FindBin(MOM);
        CKERR = hData_pL_Stat->GetBinError(CKBIN)*UncertaintyScale[uRad];
        CKVAL = rangen.Gaus(CKVAL,CKERR);
        gDummyData[uMod][uRad].SetPoint(uMom,MOM,CKVAL);
        gDummyData[uMod][uRad].SetPointError(uMom,0,CKERR);
        if(uMod==1&&uRad==1){
          hProjection->SetBinContent(uMom+1,CKVAL);
          hProjection->SetBinError(uMom+1,CKERR);
        }
        if(uMod==1&&uRad==0){
          hGenHM->SetBinContent(uMom+1,CKVAL);
          hGenHM->SetBinError(uMom+1,CKERR);
        }
      }
    }
  }
  const float right = 0.025;
  const float top = 0.025;

  gCkFine[1][0].SetLineColor(kBlack);
  gCkFine[1][0].SetLineWidth(2.5);

  gCkFine[1][1].SetLineColor(kAzure);
  gCkFine[1][1].SetLineWidth(2.5);

  gDummyData[1][0].SetLineColor(kRed);
  gDummyData[1][0].SetLineWidth(5);

  gDummyData[1][1].SetLineColor(kBlue);
  gDummyData[1][1].SetLineWidth(5);

  DLM_SubPads DlmPad(1024,720);
  DlmPad.AddSubPad(0,1,0.33,1);
  DlmPad.AddSubPad(0,1,0,0.33);
  DlmPad.SetMargin(0,0.12,0.02,0.0,0.02);
  DlmPad.SetMargin(1,0.12,0.02,0.09,0.0);
  DlmPad.cd(0);

  hData_pL_Stat->SetTitle("; #it{k*} (MeV/#it{c}); #it{C}(#it{k*})");
  hData_pL_Stat->GetXaxis()->SetRangeUser(0, kMax);
  hData_pL_Stat->GetXaxis()->SetNdivisions(505);
  hData_pL_Stat->GetYaxis()->SetRangeUser(0.85, 2.1);
  hData_pL_Stat->SetFillColor(kGray+1);
  SetStyleHisto2(hData_pL_Stat,2,0);
  SetStyleHisto2(hGenHM,2,0);
  SetStyleHisto2(hProjection,2,0);
  hData_pL_Stat->SetFillColor(kWhite);
  hGenHM->SetFillColor(kWhite);
  hProjection->SetLineColor(kAzure);
  hProjection->SetMarkerColor(kAzure);
  hProjection->SetMarkerStyle(26);
  //hData_pL_Stat->GetYaxis()->SetTitleOffset(1.0);
  hData_pL_Stat->Draw("axis");
  //gDummyData[1][0].Draw("same,e");
  if(PlotType==0) hGenHM->Draw("same");
  hProjection->Draw("same");
  if(PlotType==0) gCkFine[1][0].Draw("same,c");
  gCkFine[1][1].Draw("same,c");
  //gDummyData[1][1].Draw("same,e");

  TString LegendSource_line1 = "Gaussian source";
  //TString LegendSource_line2 = TString::Format("n#sigma#in(%.1f,%.1f)",nSigma_Best[NumBlTypes],nSigma_Worst[NumBlTypes]);

  unsigned NumRows=2;
  TLegend *legend = new TLegend(0.44,0.73-0.055-0.060*NumRows,0.78,0.73-0.055);//lbrt
  legend->SetBorderSize(0);
  legend->SetTextFont(42);
  legend->SetTextSize(gStyle->GetTextSize()*0.90);
  legend->AddEntry(&gDummyData[1][0],TString::Format("#chiEFT NLO19, r_{eff} = %.2f fm (pp HM)",Radii[0]),"l");
  legend->AddEntry(&gDummyData[1][1],TString::Format("#chiEFT NLO19, r_{eff} = %.2f fm (O+O)",Radii[1]),"l");
  //TH1F* hCk_Fake;
  //hCk_Fake = (TH1F*)hData_pL_Stat->Clone("hCk_Fake");
  //hCk_Fake->SetName("hCk_Fake");
  //hCk_Fake->SetLineColor(hCk_Fake->GetFillColor());

  //legend->AddEntry(hCk_Fake, "p#minus #Lambda #oplus #bar{p}#minus #bar{#Lambda} pairs", "fpe");
  //if(!DataOnly){
  //legend->AddEntry(&graph_FitMAX[NumBlTypes],"Femtoscopic fit (NLO)","l");
  //legend->AddEntry(&graph_BlMAX[NumBlTypes],"Baseline","l");
  //}
  //legend->Draw("same");
  TLatex BeamText;
  BeamText.SetTextSize(gStyle->GetTextSize()*1.00);
  BeamText.SetNDC(kTRUE);
  BeamText.DrawLatex(0.25, 0.895, "ALICE Update Projection");
  //BeamText.DrawLatex(0.45, 0.830, "high-mult. (0#minus 0.17% INEL>0) pp #sqrt{#it{s}} = 13 TeV");
  TLegend *legHM = new TLegend(0.22,0.830-0.055*2,0.95,0.830);//lbrt
  legHM->SetBorderSize(0);
  legHM->SetTextFont(42);
  legHM->SetTextSize(gStyle->GetTextSize()*0.90);
  legHM->AddEntry(hData_pL_Stat,TString::Format("Projected data for high-mult. (0#minus 0.17%% INEL>0) pp #sqrt{#it{s}} = 13 TeV"));
  legHM->AddEntry(&gCkFine[1][0],TString::Format("#chiEFT NLO19, #it{r_{0}}=%.2f fm",Radii[0]));
  //BeamText.DrawLatex(0.45, 0.720, "O#minusO #sqrt{#it{s}_{NN}}=6.35 TeV, #it{L}_{int}=1 nb^{-1}");
  TLegend *legOO = new TLegend(0.22,0.700-0.055*2,0.95,0.700);//lbrt
  legOO->SetBorderSize(0);
  legOO->SetTextFont(42);
  legOO->SetTextSize(gStyle->GetTextSize()*0.90);
  legOO->AddEntry(hProjection,TString::Format("Projected data for O#minusO #sqrt{#it{s}_{NN}}=6.35 TeV, #it{L}_{int}=1 nb^{-1}"));
  legOO->AddEntry(&gCkFine[1][1],TString::Format("#chiEFT NLO19, #it{r_{0}}=%.2f fm",Radii[1]));

  if(PlotType==0) legHM->Draw("same");
  legOO->Draw("same");

  //BeamText.DrawLatex(0.45, 0.850, "Data generated from theory with");
  //BeamText.DrawLatex(0.45, 0.795, "reallistic statistical uncertainties");


  TLatex BeamTextSource;
  BeamTextSource.SetTextSize(gStyle->GetTextSize()*1.00);
  BeamTextSource.SetNDC(kTRUE);
  //BeamTextSource.DrawLatex(0.45, 0.785-0.055, LegendSource_line1);
  //BeamTextSource.DrawLatex(0.40, 0.750, LegendSource_line2);

//INLET -------------------------------------------------------------------------------------------------------------------

  DlmPad.cd(1);
  TH1F* DataHisto_Inlet = (TH1F*)hData_pL_Stat->Clone("DataHisto_Inlet");
  DataHisto_Inlet->SetMarkerSize(hData_pL_Stat->GetMarkerSize()*0.67);
  DataHisto_Inlet->SetLineWidth(hData_pL_Stat->GetLineWidth()*0.67);
  DataHisto_Inlet->GetXaxis()->SetTitleSize(hData_pL_Stat->GetXaxis()->GetTitleSize()*2.0);
  DataHisto_Inlet->GetXaxis()->SetLabelSize(hData_pL_Stat->GetXaxis()->GetLabelSize()*2.0);
  DataHisto_Inlet->GetXaxis()->SetRangeUser(0, kMax);
  DataHisto_Inlet->GetXaxis()->SetNdivisions(505);

  DataHisto_Inlet->GetYaxis()->SetTitleSize(hData_pL_Stat->GetYaxis()->GetTitleSize()*2.0);
  DataHisto_Inlet->GetYaxis()->SetLabelSize(hData_pL_Stat->GetYaxis()->GetLabelSize()*2.0);
  DataHisto_Inlet->GetYaxis()->SetTitleOffset(hData_pL_Stat->GetYaxis()->GetTitleOffset()*0.5);
  DataHisto_Inlet->GetYaxis()->SetRangeUser(0.975, 1.0225);
  DataHisto_Inlet->GetYaxis()->SetNdivisions(504);

  DataHisto_Inlet->Draw("axis");
  //gDummyData[1][0].Draw("same,e");
  if(PlotType==0) hGenHM->Draw("same");
  hProjection->Draw("same");
  if(PlotType==0) gCkFine[1][0].Draw("same,c");
  gCkFine[1][1].Draw("same,c");
  //gDummyData[1][1].Draw("same,e");
  DlmPad.GetCanvas()->SaveAs(TString::Format("%s/pLambda/pLambda_DummyCk_DifferentRadii.pdf",GetFemtoOutputFolder()));


  for(unsigned uMod=0; uMod<NumModels; uMod++){
    delete [] gDummyData[uMod];
    delete [] gCk[uMod];
  }
  delete [] gCk;
  delete [] gDummyData;
  delete [] Radii;
  delete [] UncertaintyScale;
  delete [] ChiralType;
}





int PLAMBDA_1_MAIN(int argc, char *argv[]){
printf("PLAMBDA_1_MAIN\n");

pLambda_DummyCk_DifferentRadii();
return 0;

//FitMC_CompareToData_pL(0);
//FitMC_CompareToData_pL(1);
//FitMC_CompareToData_pL(2);

//EffectOfLambdaKaon();
//EffectOfpXi0("pp13TeV_HM_Dec19");
//FeedDown_pXi0_pL(atoi(argv[1]),atoi(argv[2]));
//LambdaPars("pp13TeV_HM_DimiJun20");
//LambdaPars("pp13TeV_HM_Dec19");
//Purity_vs_kstar_2(atoi(argv[1]),atoi(argv[2]),atoi(argv[3]),atoi(argv[4]));
//CompareTo_pp();
//pLambda_Study_Hypotheses();
//ParametrizeSmearMatrix1_pL();

//CreateFineBinnedResoMatrix_pp();
//CreateFineBinnedResoMatrix_pL();
//pLambda_Spline_Fit_Test();
//pLambda_Spline_Fit_Test2();
//pLambda_Spline_Fit_Unfold2(12,"L53_SL4_SR6_P96_0","/mnt/Ubuntu_Data/CernBox/Sync/CatsFiles/ExpData/ALICE_pp_13TeV_HM/DimiJun20/Norm240_340/DataSignal/Unfolded/ppMatrix/");
//const unsigned& SEEDmin, const unsigned& NumIter,
//                                const double& BinWidth, const TString& DataVariation,
//                                const char* CatsFileFolder, const TString& OutputFolder
//pLambda_Spline_Fit_Unfold2(atoi(argv[1]),atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),"L53_SL4_SR6_P96_0","/mnt/Ubuntu_Data/CernBox/Sync/CatsFiles",
//                           "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Unfolding/Test4/");
//UpdateUnfoldFile("/mnt/Ubuntu_Data/CernBox/Sync/CatsFiles",
//                 "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Unfolding/Test4/QA.root",
//                 12,"L53_SL4_SR6_P96_0");
//cout << argv[1] << endl;
//cout << argv[2] << endl;
//cout << argv[3] << endl;
//cout << argv[4] << endl;
//UpdateUnfoldFile(argv[1],argv[2],argv[3],argv[4],atoi(argv[5]));

//const int& SIGMA_FEED,const int& WhichBaseline,
//const int& WhichPotential,const float& ValSourceAlpha,
//TString InputFolder, TString InputFileName, TString OutputFolder
//POT BL SIG ALPHA(20)


//Quick_pLambda_plotter();
//Quick_pLambda_plotter_NLO13_vs_LO13();
//return 0;
//printf("hello\n");

Plot_pL_SystematicsMay2020_2(atoi(argv[3]),atoi(argv[2]),atoi(argv[1]),double(atoi(argv[4]))/10.,
                            ///home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/pL_SystematicsMay2020/BatchFarm/100720_Unfolded/
                            TString::Format("%s/pLambda/100720_Unfolded/",GetCernBoxDimi()),
                            TString::Format("Merged_pp13TeV_HM_DimiJul20_POT%i_BL%i_SIG%i.root",
                            //TString::Format("Merged_pp13TeV_HM_Dec19_POT%i_BL%i_SIG%i.root",
                            atoi(argv[1]),atoi(argv[2]),atoi(argv[3])),
                            //"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/pL_SystematicsMay2020/Test/",
                            //"UnfoldRefine_pp13TeV_HM_DimiJul20_POT11600_BL10_SIG1.root",
                            TString::Format("%s/pLambda/100720_Unfolded/PaperPlotsUpdate1/",GetCernBoxDimi()),
                            //"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/pL_SystematicsMay2020/Test/"
                            atoi(argv[5])///REMOVE FOR THE OLD PLOTS
                          );

//MakeLATEXtable(TString::Format("%s/pLambda/231020_NoBootstrap/Plots091220/",GetCernBoxDimi()),true);
return 0;

//Plot_pL_SystematicsMay2020_2(2,10,1500,2.0,
//        "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/pL_SystematicsMay2020/BatchFarm/040620_Gauss/",
//        "Merged_pp13TeV_HM_Dec19_POT1500_BL10_SIG2.root",
//        "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/pL_SystematicsMay2020/BatchFarm/040620_Gauss/Plots/");
//Plot_pL_SystematicsMay2020_2(1,10,11600,2.0,
//        "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/pL_SystematicsMay2020/BatchFarm/040620_Gauss/",
//        "Merged_pp13TeV_HM_Dec19_POT11600_BL10_SIG1.root",
//        "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/pL_SystematicsMay2020/BatchFarm/040620_Gauss/Plots/");


//Plot_pL_SystematicsMay2020(atoi(argv[1]),atoi(argv[2]),argv[3],argv[4]);

//Plot_pL_SystematicsMay2020(-1,0,
//"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/pL_SystematicsMay2020/BatchFarm/220520_Levy/GaussOnly/",
//"Levy_8192.root");

//Plot_pL_SystematicsMay2020(-1,1,
//"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/pL_SystematicsMay2020/BatchFarm/220520_Levy/Levy_18/",
//"Levy_8192.root");

//Plot_pL_SystematicsMay2020(-1,2,
//"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/pL_SystematicsMay2020/BatchFarm/220520_Levy/Levy_16/",
//"Levy_8192.root");

//Plot_pL_SystematicsMay2020(1000);
//Plot_pL_SystematicsMay2020(-1);
//Plot_pL_SystematicsMay2020(0);
//Plot_pL_SystematicsMay2020(1);
//Plot_pL_SystematicsMay2020(2);
//STUPED_TEST();
//CompareChiralNLO_pLambda();


//unsigned SEED, unsigned BASELINE_VAR, int POT_VAR, int Sigma0_Feed, int Data_Type,
                           //bool DataSyst, bool FitSyst, bool Bootstrap, unsigned NumIter,
                          // const char* CatsFileFolder, const char* OutputFolder
//pL_SystematicsMay2020(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7]), atoi(argv[8]), atoi(argv[9]),
//TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data(),
//TString::Format("%s/pLambda/Dump/",GetCernBoxDimi()).Data());

//pL_SystematicsMay2020(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7]),
//                      argv[8],argv[9]);
//pL_SystematicsMay2020(1, 9, 11600, 0, 64,
//"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/pL_SystematicsMay2020/Test/");
//pL_SystematicsMay2020(1, 6, 11600, false, 3,
//"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/pL_SystematicsMay2020/Test/");
//pL_SystematicsMay2020(1, 8, 11600, false, 3,
//"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/pL_SystematicsMay2020/Test/");
//pL_SystematicsMay2020(1, 9, 11600, false, 3,
//"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/pL_SystematicsMay2020/Test/");
return 0;

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
    //pp_SystematicsHM("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pp/Test/",
    //                 atoi(argv[1]),atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atoi(argv[5]),atoi(argv[6]));
//return 0;
//run for 30 and 32
    //pL_SystematicsHM("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pL/Systematics_020919/",
    //                 30,0,19845,200,1,0);
    //pL_SystematicsHM("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pL/Test/",
    //                 atoi(argv[1]),atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atoi(argv[5]),atoi(argv[6]));
    //pL_SystematicsHM("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pL/SystematicsAdd_300419/",
    //                 atoi(argv[1])+2,atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atoi(argv[5]),atoi(argv[6]));
    //pL_SystematicsHM("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pL/SystematicsAdd_300419/",
    //                 atoi(argv[1])+1,atoi(argv[2])*3,(atoi(argv[3])+1)*3-1,atoi(argv[4])*3,atoi(argv[5]),atoi(argv[6]));

    //pL_SystematicsHM("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pL/SystematicsAdd_070419/",
    //                 0,0,4,1,1,0);

    //the fast plots for Schleching
    //goal: fit pL with systematic variations
    //pL_SchlechingSystematics(11,4);

//const unsigned SEED, const unsigned NumVars_PerBL, const int& POT_VAR
    //pL_SchlechingSystematics(atoi(argv[1]),atoi(argv[2]),atoi(argv[3]));
    //pL_SchlechingPlots("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/pL_AllHaide_May2020/Output_pp13TeV_HM_Dec19_POT1600_1.root");

    //pL_SchlechingSystematics(atoi(argv[1]),atoi(argv[2]));
    //pL_SchlechingPlots("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/pL_SchlechingSystematics/Output_pp13TeV_HM_Dec19_1002.root");
    //pL_SchlechingPlots("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/pL_SchlechingSystematics/Output_pp13TeV_HM_Dec19_POT0_3.root");
    //pL_SchlechingPlots("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/pL_SchlechingSystematics/Output_pp13TeV_HM_Dec19_POT10_4.root");


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

    DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder("/mnt/Ubuntu_Data/CernBox/Sync/CatsFiles");
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
                 1, 0, "/mnt/Ubuntu_Data/CernBox/Sync/pPb/Systematics/C2totalsysPL.root"
                 );

    Fit_pp_pL("pPb5TeV_Run2paper","McGauss_Reso","Norm","Norm");
    Plot_pL_FAST("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pp_pL/250319/NLO_Coupled/",
                 "pPb5TeV_Run2paper","McGauss_Reso","Norm","Norm",
                 "LO_Coupled_S", "NLO_Coupled_S",
                 "Gauss core fixed to p-p", "",
                 1, 0, "/mnt/Ubuntu_Data/CernBox/Sync/pPb/Systematics/C2totalsysPL.root"
                 );

    Fit_pp_pL("pPb5TeV_Run2paper","McLevyNolan_Reso","Norm","Norm");
    Plot_pL_FAST("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pp_pL/250319/NLO_Coupled/",
                 "pPb5TeV_Run2paper","McLevyNolan_Reso","Norm","Norm",
                 "LO_Coupled_S", "NLO_Coupled_S",
                 "Levy core fixed to p-p", "",
                 1, 0, "/mnt/Ubuntu_Data/CernBox/Sync/pPb/Systematics/C2totalsysPL.root"
                 );
*/
/*
    Fit_pp_pL("pPb5TeV_CPR_Mar19","Gauss","Norm","Norm");
    Plot_pL_FAST("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pp_pL/020419/",
                 "pPb5TeV_CPR_Mar19","Gauss","Norm","Norm",
                 "LO_Coupled_S", "NLO_Coupled_S",
                 "Gauss source fixed to p-p (1.41 fm)", "",
                 1, 0, "/mnt/Ubuntu_Data/CernBox/Sync/SystematicsAndCalib/pPbRun2_MB_StreamFIX/C2totalsysPL.root"
                 );
*/
/*
    Fit_pp_pL("pPb5TeV_CPR_Mar19","McGauss_Reso","Norm","Norm");
    Plot_pL_FAST("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pp_pL/020419/",
                 "pPb5TeV_CPR_Mar19","McGauss_Reso","Norm","Norm",
                 "LO_Coupled_S", "NLO_Coupled_S",
                 "Gauss core fixed to p-p", "",
                 1, 0, "/mnt/Ubuntu_Data/CernBox/Sync/SystematicsAndCalib/pPbRun2_MB_StreamFIX/C2totalsysPL.root"
                 );
*/
/*
    Fit_pp_pL("pPb5TeV_CPR_Mar19","McLevyNolan_Reso","Norm","Norm");
    Plot_pL_FAST("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Fit_pp_pL/250319/NLO_Coupled/",
                 "pPb5TeV_CPR_Mar19","McLevyNolan_Reso","Norm","Norm",
                 "LO_Coupled_S", "NLO_Coupled_S",
                 "Levy core fixed to p-p", "",
                 1, 0, "/mnt/Ubuntu_Data/CernBox/Sync/SystematicsAndCalib/pPbRun2_MB_StreamFIX/C2totalsysPL.root"
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
