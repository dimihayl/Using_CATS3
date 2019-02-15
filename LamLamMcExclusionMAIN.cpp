
#include "LamLamMcExclusionMAIN.h"
#include "LamLamMcExclusion.h"

#include "DLM_Sort.h"

#include <iostream>
#include <fstream>

#include "TNtuple.h"
#include "TGraph.h"
#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TLegend.h"

using namespace std;
//the first version to run on the batch farm
//SepBL - 1 fit the BL separately
//SepBL - 0 do not fit the BL separately (i.e. global femto+bl fit)
//SepBL - -1 no bl at all
//BinWidthSyst - do we take as systematics the bin width
//RadiusSyst -> 0 = no variations, 1 only up variation, 2 up and down variations
//FitSyst - do we take as systematics the standard fit variations (i.e. femto and bl region)
//DataSets: "pp13TeV" "pPb5TeV" "pp7TeV" or some combinations of those
//the computer is either XMG or NX1
//Chi2Map: 0 do not save the chi2 map, 1 save the chi2 in the main output file, -1 save the chi2 in a common file
//          if -2 then the common file is separate for each systematic variation
//FitFun = -1 is no BL and no renorm
//          0 only renorm
//          1 full BL
//          2 full BL (2nd order)
//          3 femto+BL are next to one another (i.e. kf=kl), we have only norm
//          4 femto+BL are next to one another (i.e. kf=kl), we have the full BL
//          5 femto+BL are next to one another (i.e. kf=kl), we have the full BL, second order
//          6 BL and femto region are the same (norm)
//          7 BL and femto region are the same (1st order)
//          8 BL and femto region are the same (2st order)
//RandSyst if false => run a specific systematic variation, else at random
void LamLamConfidence_22Oct2018(const TString Computer,
                                const TString OutputFolder, const int Chi2Map, const bool Bootstrap,
                                const int NumBootIter, const double TimeLimit, const int JOBID, const int SystID, const bool RandSyst, const int SystJobID,
                                const int SepBL, const int FitFun,
                                const int BinWidthSyst, const int RadiusSyst, const int FitSyst, const TString DataSets){
/*
    ofstream myfile;
    myfile.open ((OutputFolder+TString::Format("SystId%u_SystJobId%u_JobId%u_OutputLamLamConfidenceINFO.txt",SystID,SystJobID,JOBID)).Data());
    myfile << "Computer = " << Computer << endl;
    myfile << "OutputFolder = " << OutputFolder << endl;
    myfile << "Chi2Map = " << Chi2Map << endl;
    myfile << "Bootstrap = " << Bootstrap << endl;
    //myfile << "NumSystVar = " << NumSystVar << endl;
    myfile << "NumBootIter = " << NumBootIter << endl;
    myfile << "JOBID = " << JOBID << endl;
    myfile << "SystID = " << SystID << endl;
    myfile << "SystJobID = " << SystJobID << endl;
    myfile << "SepBL = " << SepBL << endl;
    myfile << "FitFun = " << FitFun << endl;
    myfile << "BinWidthSyst = " << BinWidthSyst << endl;
    myfile << "RadiusSyst = " << RadiusSyst << endl;
    myfile << "FitSyst = " << FitSyst << endl;
    myfile << "DataSets = " << DataSets.Data() << endl;
    myfile.close();
*/
    LamLamMcExclusion ConfidenceMaker;
    //printf("Constructor done\n");
    ConfidenceMaker.SetOutputFile( (OutputFolder+TString::Format("SystId%u_SystJobId%u_JobId%u_OutputLamLamConfidence.root",SystID,SystJobID,JOBID)).Data() );
    if(Chi2Map%10==0){ConfidenceMaker.SetSaveChi2Map(false); ConfidenceMaker.SetCommonChi2MapFile(NULL);}
    else if(Chi2Map%10==1){ConfidenceMaker.SetSaveChi2Map(true); ConfidenceMaker.SetCommonChi2MapFile(NULL);}
    else if(Chi2Map%10==-1){ConfidenceMaker.SetSaveChi2Map(true); ConfidenceMaker.SetCommonChi2MapFile((OutputFolder+TString::Format("CommonFile.root")).Data());}
    else if(Chi2Map%10==-2){
        if(SystJobID==0){
            ConfidenceMaker.SetSaveChi2Map(true);
            ConfidenceMaker.SetCommonChi2MapFile((OutputFolder+TString::Format("Chi2Map_SystId%u.root",SystID)).Data());
        }
        else{
            ConfidenceMaker.SetSaveChi2Map(false);
            ConfidenceMaker.SetCommonChi2MapFile(NULL);
        }
    }
    else {printf("\033[1;31mERROR: \033[0mUnknown value %i for Chi2Map\n", Chi2Map); return;}
    if(Chi2Map/10==0) ConfidenceMaker.SetUpChi2Map(1024,-2,5,768,0,18);
    else if(abs(Chi2Map/10)==1) ConfidenceMaker.SetUpChi2Map(512,-2,5,384,0,18);
    else if(abs(Chi2Map/10)==2) ConfidenceMaker.SetUpChi2Map(256,-2,5,192,0,18);
    else if(abs(Chi2Map/10)==3) ConfidenceMaker.SetUpChi2Map(128,-2,5,96,0,18);
    else ConfidenceMaker.SetUpChi2Map(64,-2,5,48,0,18);
    ConfidenceMaker.SetRandomSeed(abs(SystID+1),abs(JOBID+1));
    ConfidenceMaker.SetBootstrap(Bootstrap);
    ConfidenceMaker.SetNumSystVariations(1);
    if(!RandSyst){
        ConfidenceMaker.RunSpecificSyst(true,SystID,SystID);
//ConfidenceMaker.RunSpecificSyst(true,0,53);
    }

    //run the default systematics with the first job
    ConfidenceMaker.SetRunDefaultSyst(SystID==0);
    //ConfidenceMaker.SetRunDefaultSyst(false);
    ConfidenceMaker.SetNumBootIter(NumBootIter);
    ConfidenceMaker.SetTimeLimit(TimeLimit);

    TString DataFileName_pp13TeV_20MeV;
    TString DataHistoName_pp13TeV_20MeV;
    TString ResFileName_pp13TeV_20MeV;
    TString ResHistoName_pp13TeV_20MeV;

    TString DataFileName_pp13TeV_16MeV;
    TString DataHistoName_pp13TeV_16MeV;
    TString ResFileName_pp13TeV_16MeV;
    TString ResHistoName_pp13TeV_16MeV;

    TString DataFileName_pp13TeV_12MeV;
    TString DataHistoName_pp13TeV_12MeV;
    TString ResFileName_pp13TeV_12MeV;
    TString ResHistoName_pp13TeV_12MeV;

    TString DataFileName_pPb5TeV_20MeV;
    TString DataHistoName_pPb5TeV_20MeV;
    TString ResFileName_pPb5TeV_20MeV;
    TString ResHistoName_pPb5TeV_20MeV;

    TString DataFileName_pPb5TeV_16MeV;
    TString DataHistoName_pPb5TeV_16MeV;
    TString ResFileName_pPb5TeV_16MeV;
    TString ResHistoName_pPb5TeV_16MeV;

    TString DataFileName_pPb5TeV_12MeV;
    TString DataHistoName_pPb5TeV_12MeV;
    TString ResFileName_pPb5TeV_12MeV;
    TString ResHistoName_pPb5TeV_12MeV;

    TString DataFileName_pp7TeV_20MeV;
    TString DataHistoName_pp7TeV_20MeV;
    TString ResFileName_pp7TeV_20MeV;
    TString ResHistoName_pp7TeV_20MeV;

    if(Computer=="XMG"){
        DataFileName_pp13TeV_20MeV="/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample7/Data/CFOutput_LLMEReweighted_Rebin_5.root";
        DataHistoName_pp13TeV_20MeV="hCkTotNormWeight";
        ResFileName_pp13TeV_20MeV="/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample7/MC/hResMatr_LL_Rebin_5.root";
        ResHistoName_pp13TeV_20MeV="hResMatrix_5";

        DataFileName_pp13TeV_16MeV="/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample7/Data/CFOutput_LLMEReweighted_Rebin_4.root";
        DataHistoName_pp13TeV_16MeV="hCkTotNormWeight";
        ResFileName_pp13TeV_16MeV="/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample7/MC/hResMatr_LL_Rebin_4.root";
        ResHistoName_pp13TeV_16MeV="hResMatrix_4";

        DataFileName_pp13TeV_12MeV="/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample7/Data/CFOutput_LLMEReweighted_Rebin_3.root";
        DataHistoName_pp13TeV_12MeV="hCkTotNormWeight";
        ResFileName_pp13TeV_12MeV="/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample7/MC/hResMatr_LL_Rebin_3.root";
        ResHistoName_pp13TeV_12MeV="hResMatrix_3";

        DataFileName_pPb5TeV_20MeV="/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pPb_5TeV/Sample7/Data/CFOutput_LLMEReweighted_Rebin_5.root";
        DataHistoName_pPb5TeV_20MeV="hCkTotNormWeight";
        ResFileName_pPb5TeV_20MeV="/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pPb_5TeV/Sample7/MC/hResMatr_LL_Rebin_5.root";
        ResHistoName_pPb5TeV_20MeV="hResMatrix_5";

        DataFileName_pPb5TeV_16MeV="/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pPb_5TeV/Sample7/Data/CFOutput_LLMEReweighted_Rebin_4.root";
        DataHistoName_pPb5TeV_16MeV="hCkTotNormWeight";
        ResFileName_pPb5TeV_16MeV="/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pPb_5TeV/Sample7/MC/hResMatr_LL_Rebin_4.root";
        ResHistoName_pPb5TeV_16MeV="hResMatrix_4";

        DataFileName_pPb5TeV_12MeV="/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pPb_5TeV/Sample7/Data/CFOutput_LLMEReweighted_Rebin_3.root";
        DataHistoName_pPb5TeV_12MeV="hCkTotNormWeight";
        ResFileName_pPb5TeV_12MeV="/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pPb_5TeV/Sample7/MC/hResMatr_LL_Rebin_3.root";
        ResHistoName_pPb5TeV_12MeV="hResMatrix_3";

        //the old data
        //DataFileName_pp7TeV_20MeV="/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_7TeV/Range200_400/LHC10bcde_merge_CkLambdaLambda_0.root";
        //DataHistoName_pp7TeV_20MeV="hCkTotNormWeight";
        //ResFileName_pp7TeV_20MeV="/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_7TeV/ResolutionMatrices/LHC14j4bcde_MeV_compact.root";
        //ResHistoName_pp7TeV_20MeV="hSigma_pL_InMeV";

        //the new data, i.e. with multiplicity re-weight
        DataFileName_pp7TeV_20MeV="/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_7TeV/301018/CFOutput_LL.root";
        DataHistoName_pp7TeV_20MeV="hCkTotReweightMeV";
        ResFileName_pp7TeV_20MeV="/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_7TeV/ResolutionMatrices/LHC14j4bcde_MeV_compact.root";
        ResHistoName_pp7TeV_20MeV="hSigma_pL_InMeV";

        //DataFileName_pp7TeV_20MeV="/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample7/Data/CFOutput_LLMEReweighted_Rebin_5.root";
        //DataHistoName_pp7TeV_20MeV="hCkTotNormWeight";
        //ResFileName_pp7TeV_20MeV="/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample7/MC/hResMatr_LL_Rebin_5.root";
        //ResHistoName_pp7TeV_20MeV="hResMatrix_5";
    }
    else if(Computer=="NX1"){
        DataFileName_pp13TeV_20MeV="/scratch6/dmihaylov/CorrelationFiles_2018/ALICE_pp_13TeV/Sample7/Data/CFOutput_LLMEReweighted_Rebin_5.root";
        DataHistoName_pp13TeV_20MeV="hCkTotNormWeight";
        ResFileName_pp13TeV_20MeV="/scratch6/dmihaylov/CorrelationFiles_2018/ALICE_pp_13TeV/Sample7/MC/hResMatr_LL_Rebin_5.root";
        ResHistoName_pp13TeV_20MeV="hResMatrix_5";

        DataFileName_pp13TeV_16MeV="/scratch6/dmihaylov/CorrelationFiles_2018/ALICE_pp_13TeV/Sample7/Data/CFOutput_LLMEReweighted_Rebin_4.root";
        DataHistoName_pp13TeV_16MeV="hCkTotNormWeight";
        ResFileName_pp13TeV_16MeV="/scratch6/dmihaylov/CorrelationFiles_2018/ALICE_pp_13TeV/Sample7/MC/hResMatr_LL_Rebin_4.root";
        ResHistoName_pp13TeV_16MeV="hResMatrix_4";

        DataFileName_pp13TeV_12MeV="/scratch6/dmihaylov/CorrelationFiles_2018/ALICE_pp_13TeV/Sample7/Data/CFOutput_LLMEReweighted_Rebin_3.root";
        DataHistoName_pp13TeV_12MeV="hCkTotNormWeight";
        ResFileName_pp13TeV_12MeV="/scratch6/dmihaylov/CorrelationFiles_2018/ALICE_pp_13TeV/Sample7/MC/hResMatr_LL_Rebin_3.root";
        ResHistoName_pp13TeV_12MeV="hResMatrix_3";

        DataFileName_pPb5TeV_20MeV="/scratch6/dmihaylov/CorrelationFiles_2018/ALICE_pPb_5TeV/Sample7/Data/CFOutput_LLMEReweighted_Rebin_5.root";
        DataHistoName_pPb5TeV_20MeV="hCkTotNormWeight";
        ResFileName_pPb5TeV_20MeV="/scratch6/dmihaylov/CorrelationFiles_2018/ALICE_pPb_5TeV/Sample7/MC/hResMatr_LL_Rebin_5.root";
        ResHistoName_pPb5TeV_20MeV="hResMatrix_5";

        DataFileName_pPb5TeV_16MeV="/scratch6/dmihaylov/CorrelationFiles_2018/ALICE_pPb_5TeV/Sample7/Data/CFOutput_LLMEReweighted_Rebin_4.root";
        DataHistoName_pPb5TeV_16MeV="hCkTotNormWeight";
        ResFileName_pPb5TeV_16MeV="/scratch6/dmihaylov/CorrelationFiles_2018/ALICE_pPb_5TeV/Sample7/MC/hResMatr_LL_Rebin_4.root";
        ResHistoName_pPb5TeV_16MeV="hResMatrix_4";

        DataFileName_pPb5TeV_12MeV="/scratch6/dmihaylov/CorrelationFiles_2018/ALICE_pPb_5TeV/Sample7/Data/CFOutput_LLMEReweighted_Rebin_3.root";
        DataHistoName_pPb5TeV_12MeV="hCkTotNormWeight";
        ResFileName_pPb5TeV_12MeV="/scratch6/dmihaylov/CorrelationFiles_2018/ALICE_pPb_5TeV/Sample7/MC/hResMatr_LL_Rebin_3.root";
        ResHistoName_pPb5TeV_12MeV="hResMatrix_3";

        //the old data
        //DataFileName_pp7TeV_20MeV="/scratch6/dmihaylov/CorrelationFiles_2018/ALICE_pp_7TeV/Range200_400/LHC10bcde_merge_CkLambdaLambda_0.root";
        //DataHistoName_pp7TeV_20MeV="hCkTotNormWeight";
        //ResFileName_pp7TeV_20MeV="/scratch6/dmihaylov/CorrelationFiles_2018/ALICE_pp_7TeV/ResolutionMatrices/LHC14j4bcde_MeV_compact.root";
        //ResHistoName_pp7TeV_20MeV="hSigma_pL_InMeV";

        //the new data, i.e. with multiplicity re-weight
        DataFileName_pp7TeV_20MeV="/scratch6/dmihaylov/CorrelationFiles_2018/ALICE_pp_7TeV/301018/CFOutput_LL.root";
        DataHistoName_pp7TeV_20MeV="hCkTotReweightMeV";
        ResFileName_pp7TeV_20MeV="/scratch6/dmihaylov/CorrelationFiles_2018/ALICE_pp_7TeV/ResolutionMatrices/LHC14j4bcde_MeV_compact.root";
        ResHistoName_pp7TeV_20MeV="hSigma_pL_InMeV";
    }
    else{
        printf("\033[1;31mERROR:\033[0m %s is an unknown machine!\n",Computer.Data());
        return;
    }


    TFile* fd_pp13TeV_20MeV = new TFile(DataFileName_pp13TeV_20MeV, "read");
    TH1F* hd_pp13TeV_20MeV = (TH1F*)fd_pp13TeV_20MeV->Get(DataHistoName_pp13TeV_20MeV);
    TFile* fr_pp13TeV_20MeV = new TFile(ResFileName_pp13TeV_20MeV, "read");
    TH2F* hr_pp13TeV_20MeV = (TH2F*)fr_pp13TeV_20MeV->Get(ResHistoName_pp13TeV_20MeV);

    TFile* fd_pp13TeV_16MeV = new TFile(DataFileName_pp13TeV_16MeV, "read");
    TH1F* hd_pp13TeV_16MeV = (TH1F*)fd_pp13TeV_16MeV->Get(DataHistoName_pp13TeV_16MeV);
    TFile* fr_pp13TeV_16MeV = new TFile(ResFileName_pp13TeV_16MeV, "read");
    TH2F* hr_pp13TeV_16MeV = (TH2F*)fr_pp13TeV_16MeV->Get(ResHistoName_pp13TeV_16MeV);

    TFile* fd_pp13TeV_12MeV = new TFile(DataFileName_pp13TeV_12MeV, "read");
    TH1F* hd_pp13TeV_12MeV = (TH1F*)fd_pp13TeV_12MeV->Get(DataHistoName_pp13TeV_12MeV);
    TFile* fr_pp13TeV_12MeV = new TFile(ResFileName_pp13TeV_12MeV, "read");
    TH2F* hr_pp13TeV_12MeV = (TH2F*)fr_pp13TeV_12MeV->Get(ResHistoName_pp13TeV_12MeV);

    TFile* fd_pPb5TeV_20MeV = new TFile(DataFileName_pPb5TeV_20MeV, "read");
    TH1F* hd_pPb5TeV_20MeV = (TH1F*)fd_pPb5TeV_20MeV->Get(DataHistoName_pPb5TeV_20MeV);
    TFile* fr_pPb5TeV_20MeV = new TFile(ResFileName_pPb5TeV_20MeV, "read");
    TH2F* hr_pPb5TeV_20MeV = (TH2F*)fr_pPb5TeV_20MeV->Get(ResHistoName_pPb5TeV_20MeV);

    TFile* fd_pPb5TeV_16MeV = new TFile(DataFileName_pPb5TeV_16MeV, "read");
    TH1F* hd_pPb5TeV_16MeV = (TH1F*)fd_pPb5TeV_16MeV->Get(DataHistoName_pPb5TeV_16MeV);
    TFile* fr_pPb5TeV_16MeV = new TFile(ResFileName_pPb5TeV_16MeV, "read");
    TH2F* hr_pPb5TeV_16MeV = (TH2F*)fr_pPb5TeV_16MeV->Get(ResHistoName_pPb5TeV_16MeV);

    TFile* fd_pPb5TeV_12MeV = new TFile(DataFileName_pPb5TeV_12MeV, "read");
    TH1F* hd_pPb5TeV_12MeV = (TH1F*)fd_pPb5TeV_12MeV->Get(DataHistoName_pPb5TeV_12MeV);
    TFile* fr_pPb5TeV_12MeV = new TFile(ResFileName_pPb5TeV_12MeV, "read");
    TH2F* hr_pPb5TeV_12MeV = (TH2F*)fr_pPb5TeV_12MeV->Get(ResHistoName_pPb5TeV_12MeV);

    TFile* fd_pp7TeV_20MeV = new TFile(DataFileName_pp7TeV_20MeV, "read");
    TH1F* hd_pp7TeV_20MeV = (TH1F*)fd_pp7TeV_20MeV->Get(DataHistoName_pp7TeV_20MeV);
    TFile* fr_pp7TeV_20MeV = new TFile(ResFileName_pp7TeV_20MeV, "read");
    TH2F* hr_pp7TeV_20MeV = (TH2F*)fr_pp7TeV_20MeV->Get(ResHistoName_pp7TeV_20MeV);

    // 0 = RUN2 13 TeV
    // 1 = RUN2 5 TeV
    // 2 = RUN1 7 TeV
    unsigned NumDataSets = 0;
    //bool UseDataSet[3];
    //later WhichDataSet[uData] will give me the actual ID when I iterate over NumDataSets
    unsigned WhichDataSet[3];
    if(DataSets.Contains("pp13TeV")){
        WhichDataSet[NumDataSets]=0;
        NumDataSets++;
        //UseDataSet[0]=true;
    }
    else{
        //UseDataSet[0]=false;
    }

    if(DataSets.Contains("pPb5TeV")){
        WhichDataSet[NumDataSets]=1;
        NumDataSets++;
        //UseDataSet[1]=true;
    }
    else{
        //UseDataSet[1]=false;
    }

    if(DataSets.Contains("pp7TeV")){
        WhichDataSet[NumDataSets]=2;
        NumDataSets++;
        //UseDataSet[2]=true;
    }
    else{
        //UseDataSet[2]=false;
    }

    ConfidenceMaker.SetNumDataSets(NumDataSets);
    // 0 = 20 MeV bins
    // 1 = 16 MeV bins
    // 2 = 12 MeV bins

    //unsigned DataSetID=0;
    unsigned NumDataVar=0;
    unsigned NumFemtoVar=0;
    unsigned NumBlVar=0;
    unsigned NumRadVar=0;
//pp13TeV

    for(unsigned uData=0; uData<NumDataSets; uData++){
        printf("uData=%u\n",uData);
        //no variation for the 7 TeV sample
        if(BinWidthSyst&&WhichDataSet[uData]!=2) NumDataVar=3;
        else NumDataVar=1;
        ConfidenceMaker.SetNumDataVar(uData,NumDataVar);
        ConfidenceMaker.SetNumFitFunVar(uData,1);
        ConfidenceMaker.SetFitFunction(uData,0,FitFun);
        //ConfidenceMaker.SetFitFunction(uData,0,uData==2?1:-1);
        ConfidenceMaker.SetSeparateBaseline(uData,SepBL);

        if(FitSyst){
            if(FitFun<3){
                NumFemtoVar=3; NumBlVar=3;
            }
            //if femto and BL are glued together
            else if(FitFun<6){
                NumFemtoVar=1; NumBlVar=3;
            }
            else{
                NumFemtoVar=3; NumBlVar=1;
            }
        }
        else {NumFemtoVar=1; NumBlVar=1;}

        if(RadiusSyst==0) NumRadVar=1;
        else if(RadiusSyst==1) NumRadVar=2;
        else if(RadiusSyst==2) NumRadVar=3;

        for(unsigned uDataVar=0; uDataVar<NumDataVar; uDataVar++){
            ConfidenceMaker.SetNumFemtoReg(uData,uDataVar,NumFemtoVar);
            for(unsigned uFemtoVar=0; uFemtoVar<NumFemtoVar; uFemtoVar++){
                     if(uDataVar==0 && uFemtoVar==0) ConfidenceMaker.SetFemtoRegion(uData,uDataVar,uFemtoVar,0,460);//220
                else if(uDataVar==0 && uFemtoVar==1) ConfidenceMaker.SetFemtoRegion(uData,uDataVar,uFemtoVar,0,420);//180
                else if(uDataVar==0 && uFemtoVar==2) ConfidenceMaker.SetFemtoRegion(uData,uDataVar,uFemtoVar,0,500);//260
                else if(uDataVar==1 && uFemtoVar==0) ConfidenceMaker.SetFemtoRegion(uData,uDataVar,uFemtoVar,0,464);//224
                else if(uDataVar==1 && uFemtoVar==1) ConfidenceMaker.SetFemtoRegion(uData,uDataVar,uFemtoVar,0,416);//176
                else if(uDataVar==1 && uFemtoVar==2) ConfidenceMaker.SetFemtoRegion(uData,uDataVar,uFemtoVar,0,496);//256
                else if(uDataVar==2 && uFemtoVar==0) ConfidenceMaker.SetFemtoRegion(uData,uDataVar,uFemtoVar,0,456);//216
                else if(uDataVar==2 && uFemtoVar==1) ConfidenceMaker.SetFemtoRegion(uData,uDataVar,uFemtoVar,0,420);//180
                else if(uDataVar==2 && uFemtoVar==2) ConfidenceMaker.SetFemtoRegion(uData,uDataVar,uFemtoVar,0,504);//264
                else printf("VERY VERY BAD!\n");
            }//uFemtoVar
            ConfidenceMaker.SetNumBlReg(uData,uDataVar,NumBlVar);
            for(unsigned uBlVar=0; uBlVar<NumBlVar; uBlVar++){
                if(FitFun<3){
                         if(uDataVar==0 && uBlVar==0) ConfidenceMaker.SetBlRegion(uData,uDataVar,uBlVar,300,500);
                    else if(uDataVar==0 && uBlVar==1) ConfidenceMaker.SetBlRegion(uData,uDataVar,uBlVar,320,480);
                    else if(uDataVar==0 && uBlVar==2) ConfidenceMaker.SetBlRegion(uData,uDataVar,uBlVar,300,540);
                    else if(uDataVar==1 && uBlVar==0) ConfidenceMaker.SetBlRegion(uData,uDataVar,uBlVar,304,496);
                    else if(uDataVar==1 && uBlVar==1) ConfidenceMaker.SetBlRegion(uData,uDataVar,uBlVar,320,480);
                    else if(uDataVar==1 && uBlVar==2) ConfidenceMaker.SetBlRegion(uData,uDataVar,uBlVar,304,544);
                    else if(uDataVar==2 && uBlVar==0) ConfidenceMaker.SetBlRegion(uData,uDataVar,uBlVar,300,504);
                    else if(uDataVar==2 && uBlVar==1) ConfidenceMaker.SetBlRegion(uData,uDataVar,uBlVar,324,480);
                    else if(uDataVar==2 && uBlVar==2) ConfidenceMaker.SetBlRegion(uData,uDataVar,uBlVar,300,540);
                    else printf("VERY VERY BAD!\n");
                }
                else if(FitFun<6){
                         if(uDataVar==0 && uBlVar==0) ConfidenceMaker.SetBlRegion(uData,uDataVar,uBlVar,220,460);
                    else if(uDataVar==0 && uBlVar==1) ConfidenceMaker.SetBlRegion(uData,uDataVar,uBlVar,220,420);
                    else if(uDataVar==0 && uBlVar==2) ConfidenceMaker.SetBlRegion(uData,uDataVar,uBlVar,220,500);
                    else if(uDataVar==1 && uBlVar==0) ConfidenceMaker.SetBlRegion(uData,uDataVar,uBlVar,224,464);
                    else if(uDataVar==1 && uBlVar==1) ConfidenceMaker.SetBlRegion(uData,uDataVar,uBlVar,224,416);
                    else if(uDataVar==1 && uBlVar==2) ConfidenceMaker.SetBlRegion(uData,uDataVar,uBlVar,224,496);
                    else if(uDataVar==2 && uBlVar==0) ConfidenceMaker.SetBlRegion(uData,uDataVar,uBlVar,216,456);
                    else if(uDataVar==2 && uBlVar==1) ConfidenceMaker.SetBlRegion(uData,uDataVar,uBlVar,216,420);
                    else if(uDataVar==2 && uBlVar==2) ConfidenceMaker.SetBlRegion(uData,uDataVar,uBlVar,216,504);
                    else printf("VERY VERY BAD!\n");
                }
                else{
                    if(uDataVar==0 && uBlVar==0) ConfidenceMaker.SetBlRegion(uData,uDataVar,uBlVar,500,500);
                    else if(uDataVar==1 && uBlVar==0) ConfidenceMaker.SetBlRegion(uData,uDataVar,uBlVar,496,496);
                    else if(uDataVar==2 && uBlVar==0) ConfidenceMaker.SetBlRegion(uData,uDataVar,uBlVar,504,504);
                    else printf("VERY VERY BAD!\n");
                }

            }//uBlVar

                 if(WhichDataSet[uData]==0 && uDataVar==0) ConfidenceMaker.SetData(uData,uDataVar,hd_pp13TeV_20MeV);
            else if(WhichDataSet[uData]==0 && uDataVar==1) ConfidenceMaker.SetData(uData,uDataVar,hd_pp13TeV_16MeV);
            else if(WhichDataSet[uData]==0 && uDataVar==2) ConfidenceMaker.SetData(uData,uDataVar,hd_pp13TeV_12MeV);
            else if(WhichDataSet[uData]==1 && uDataVar==0) ConfidenceMaker.SetData(uData,uDataVar,hd_pPb5TeV_20MeV);
            else if(WhichDataSet[uData]==1 && uDataVar==1) ConfidenceMaker.SetData(uData,uDataVar,hd_pPb5TeV_16MeV);
            else if(WhichDataSet[uData]==1 && uDataVar==2) ConfidenceMaker.SetData(uData,uDataVar,hd_pPb5TeV_12MeV);
            else if(WhichDataSet[uData]==2 && uDataVar==0) ConfidenceMaker.SetData(uData,uDataVar,hd_pp7TeV_20MeV);
            //just in case for the future, anyway now for 7 TeV we have uDataVar<=0
            //else if(WhichDataSet[uData]==2 && uDataVar==1) ConfidenceMaker.SetData(uData,uDataVar,hd_pp7TeV_16MeV);
            //else if(WhichDataSet[uData]==2 && uDataVar==2) ConfidenceMaker.SetData(uData,uDataVar,hd_pp7TeV_12MeV);
            else printf("VERY VERY BAD!\n");

                 if(WhichDataSet[uData]==0 && uDataVar==0) ConfidenceMaker.SetResolution(uData,uDataVar,hr_pp13TeV_20MeV);
            else if(WhichDataSet[uData]==0 && uDataVar==1) ConfidenceMaker.SetResolution(uData,uDataVar,hr_pp13TeV_16MeV);
            else if(WhichDataSet[uData]==0 && uDataVar==2) ConfidenceMaker.SetResolution(uData,uDataVar,hr_pp13TeV_12MeV);
            else if(WhichDataSet[uData]==1 && uDataVar==0) ConfidenceMaker.SetResolution(uData,uDataVar,hr_pPb5TeV_20MeV);
            else if(WhichDataSet[uData]==1 && uDataVar==1) ConfidenceMaker.SetResolution(uData,uDataVar,hr_pPb5TeV_16MeV);
            else if(WhichDataSet[uData]==1 && uDataVar==2) ConfidenceMaker.SetResolution(uData,uDataVar,hr_pPb5TeV_12MeV);
            else if(WhichDataSet[uData]==2 && uDataVar==0) ConfidenceMaker.SetResolution(uData,uDataVar,hr_pp7TeV_20MeV);
            //just in case for the future, anyway now for 7 TeV we have uDataVar<=0
            //else if(WhichDataSet[uData]==2 && uDataVar==1) ConfidenceMaker.SetResolution(uData,uDataVar,hr_pp7TeV_16MeV);
            //else if(WhichDataSet[uData]==2 && uDataVar==2) ConfidenceMaker.SetResolution(uData,uDataVar,hr_pp7TeV_12MeV);
            else printf("VERY VERY BAD!\n");

            ConfidenceMaker.SetNumLamPar(uData,uDataVar,1);
                 if(WhichDataSet[uData]==0 && uDataVar==0) ConfidenceMaker.SetLambdaPar(uData,uDataVar,0,0.338239,0.598166);
            else if(WhichDataSet[uData]==0 && uDataVar==1) ConfidenceMaker.SetLambdaPar(uData,uDataVar,0,0.338239,0.598166);
            else if(WhichDataSet[uData]==0 && uDataVar==2) ConfidenceMaker.SetLambdaPar(uData,uDataVar,0,0.338239,0.598166);
            else if(WhichDataSet[uData]==1 && uDataVar==0) ConfidenceMaker.SetLambdaPar(uData,uDataVar,0,0.239101,0.640295);
            else if(WhichDataSet[uData]==1 && uDataVar==1) ConfidenceMaker.SetLambdaPar(uData,uDataVar,0,0.239101,0.640295);
            else if(WhichDataSet[uData]==1 && uDataVar==2) ConfidenceMaker.SetLambdaPar(uData,uDataVar,0,0.239101,0.640295);
            else if(WhichDataSet[uData]==2 && uDataVar==0) ConfidenceMaker.SetLambdaPar(uData,uDataVar,0,0.2994,0.6222);
        //else if(WhichDataSet[uData]==2 && uDataVar==0) ConfidenceMaker.SetLambdaPar(uData,uDataVar,0,0.338239,0.598166);
            else printf("VERY VERY BAD!\n");
        }//uDataVar
        ConfidenceMaker.SetNumRadVar(uData,NumRadVar);
        for(unsigned uRad=0; uRad<NumRadVar; uRad++){
                 if(WhichDataSet[uData]==0 && uRad==0) ConfidenceMaker.SetRadius(uData,uRad,1.152);//stat+/-0.009; syst+0.025 syst-0.001 and the +5%(+0.058) from resonances
            else if(WhichDataSet[uData]==0 && uRad==1) ConfidenceMaker.SetRadius(uData,uRad,1.244);
            else if(WhichDataSet[uData]==0 && uRad==2) ConfidenceMaker.SetRadius(uData,uRad,1.142);//difference from min to max is 8.9%
            else if(WhichDataSet[uData]==1 && uRad==0) ConfidenceMaker.SetRadius(uData,uRad,1.401);//stat+/-0.008; syst+0.013 syst-0.006 and the +5%(+0.070) from resonances
            else if(WhichDataSet[uData]==1 && uRad==1) ConfidenceMaker.SetRadius(uData,uRad,1.492);
            else if(WhichDataSet[uData]==1 && uRad==2) ConfidenceMaker.SetRadius(uData,uRad,1.387);//difference from min to max is 8.1%
            else if(WhichDataSet[uData]==2 && uRad==0) ConfidenceMaker.SetRadius(uData,uRad,1.125);//stat+/-0.018; syst+0.058 syst-0.035 and the +5%(+0.056) from resonances
        //else if(WhichDataSet[uData]==2 && uRad==0) ConfidenceMaker.SetRadius(uData,uRad,1.152);
            else if(WhichDataSet[uData]==2 && uRad==1) ConfidenceMaker.SetRadius(uData,uRad,1.257);
            else if(WhichDataSet[uData]==2 && uRad==2) ConfidenceMaker.SetRadius(uData,uRad,1.072);//difference from min to max is 17.2%
            else printf("VERY VERY BAD!\n");
        }
    }//uData

    printf("Num Iter = %u\n",ConfidenceMaker.GetTotNumIter());

ConfidenceMaker.SetSaveFitResults(true);//!

    ConfidenceMaker.GoBabyGo();

    delete fd_pp13TeV_20MeV;
    delete fr_pp13TeV_20MeV;
    delete fd_pp13TeV_16MeV;
    delete fr_pp13TeV_16MeV;
    delete fd_pp13TeV_12MeV;
    delete fr_pp13TeV_12MeV;

    delete fd_pPb5TeV_20MeV;
    delete fr_pPb5TeV_20MeV;
    delete fd_pPb5TeV_16MeV;
    delete fr_pPb5TeV_16MeV;
    delete fd_pPb5TeV_12MeV;
    delete fr_pPb5TeV_12MeV;

    delete fd_pp7TeV_20MeV;
    delete fr_pp7TeV_20MeV;
//printf(" SamplingHisto&&SamplingHisto[1] = %p && %p\n",ConfidenceMaker.SamplingHisto,ConfidenceMaker.SamplingHisto?ConfidenceMaker.SamplingHisto[1]:NULL);
    //printf("Ending it\n");
}



void FinalExclusionPlot(const TString Computer,const int Chi2Map,const TString OutputFolder,const unsigned NumSyst){

    unsigned Num_f0_bins;
    unsigned Num_d0_bins;
    if(Chi2Map/10==0) {Num_f0_bins=1024; Num_d0_bins=768;}
    else if(abs(Chi2Map/10)==1) {Num_f0_bins=512; Num_d0_bins=384;}
    else if(abs(Chi2Map/10)==2) {Num_f0_bins=256; Num_d0_bins=192;}
    else if(abs(Chi2Map/10)==3) {Num_f0_bins=128; Num_d0_bins=96;}
    else {Num_f0_bins=64; Num_d0_bins=48;}
    const double f0_min = -2;
    const double f0_max = 5;
    const double d0_min = 0;
    const double d0_max = 18;

    TFile* OutFile = new TFile(OutputFolder+"FinalEP.root","RECREATE");

    //hFinalEP->Write();
    //OutFile->Close();
    //delete hFinalEP;
    //delete OutFile;

    double* dChi2Val;
    TH2F** hConfLev = new TH2F* [NumSyst];
    TH2F** hConfidence = new TH2F* [NumSyst];

    double PVAL;
    double NSIGMA;

    for(unsigned uSyst=0; uSyst<NumSyst; uSyst++){
        printf("uSyst = %u (%u)\n",uSyst,NumSyst);
        TFile* File_dChi2 = new TFile(OutputFolder+TString::Format("MergedSystId%u.root",uSyst),"read");
        TNtuple* nt_dChi2 = (TNtuple*)File_dChi2->Get("ntParDensity");
        unsigned NumEntries = nt_dChi2->GetEntries();
        TGraph gCDF;
        gCDF.SetName("gCDF");
        gCDF.Set(NumEntries+2);
        dChi2Val = new double [NumEntries];
        Float_t dChi2Ndf;
        nt_dChi2->SetBranchAddress("dChi2Ndf",&dChi2Ndf);
        for(unsigned uEntry=0; uEntry<NumEntries; uEntry++){
            nt_dChi2->GetEntry(uEntry);
            dChi2Val[uEntry] = dChi2Ndf;
        }
        //delete nt_dChi2;
        delete File_dChi2;

        DLM_Sort < double, unsigned > SortTool;
        SortTool.SetData(dChi2Val,NumEntries);
        SortTool.MergeSort(false);
        SortTool.GetSortedData(dChi2Val,dChi2Val);

        //unsigned TotEntries=0;
        double CDF_VAL;
        double dChi2Ndf_Min = dChi2Val[0];
        double dChi2Ndf_Max = dChi2Val[NumEntries-1];
        gCDF.SetPoint(0,dChi2Ndf_Min-0.1,0);
        for(unsigned uEntry=0; uEntry<NumEntries; uEntry++){
            CDF_VAL = double(uEntry+1)/double(NumEntries);
            gCDF.SetPoint(uEntry+1,dChi2Val[uEntry],CDF_VAL);
            //if(CDF_VAL<0) printf("CDF_VAL=%f\n",CDF_VAL);
        }
        gCDF.SetPoint(NumEntries+1,dChi2Ndf_Max+0.1,1);

        TFile* File_Chi2Map = new TFile(OutputFolder+TString::Format("Chi2Map_SystId%u.root",uSyst),"read");
        TH2F* hChi2NdfMap = (TH2F*)File_Chi2Map->Get("hChi2NdfMap");

        OutFile->cd();
        hConfLev[uSyst] = new TH2F(TString::Format("hConfLev_%u",uSyst),TString::Format("hConfLev_%u",uSyst),Num_f0_bins,f0_min,f0_max,Num_d0_bins,d0_min,d0_max);
        hConfidence[uSyst] = new TH2F(TString::Format("hConfidence_%u",uSyst),TString::Format("hConfidence_%u",uSyst),Num_f0_bins,f0_min,f0_max,Num_d0_bins,d0_min,d0_max);

        double ConfidenceLevel;
        for(unsigned uBinX=0; uBinX<Num_f0_bins; uBinX++){
            //printf("uBinX=%u\n",uBinX);
            for(unsigned uBinY=0; uBinY<Num_d0_bins; uBinY++){
                dChi2Ndf = hChi2NdfMap->GetBinContent(uBinX+1,uBinY+1);
                if(dChi2Ndf<dChi2Ndf_Min){
                    //printf("dChi2Ndf=%f (%f)\n",dChi2Ndf,dChi2Ndf_Min);
                    dChi2Ndf=dChi2Ndf_Min;
                }
                if(dChi2Ndf>dChi2Ndf_Max){
                    //printf("dChi2Ndf=%f (%f)\n",dChi2Ndf,dChi2Ndf_Max);
                    dChi2Ndf=dChi2Ndf_Max;
                }
                ConfidenceLevel = gCDF.Eval(dChi2Ndf);
                if(ConfidenceLevel<0) printf("ConfidenceLevel=%f\n",ConfidenceLevel);
                hConfLev[uSyst]->SetBinContent(uBinX+1,uBinY+1,ConfidenceLevel);

                PVAL = 1.-ConfidenceLevel;
                if(PVAL>0 && PVAL<1)NSIGMA = sqrt(2)*TMath::ErfcInverse(PVAL);
                else if(PVAL<=0) NSIGMA = 10;
                else NSIGMA = 0;
                if(NSIGMA>10 || NSIGMA!=NSIGMA) NSIGMA = 10;

                hConfidence[uSyst]->SetBinContent(uBinX+1,uBinY+1,NSIGMA);
            }
        }

        //TFile* Chi2MapFile = new TFile(Chi2MapFileName,"UPDATE");
        delete File_Chi2Map;
        delete [] dChi2Val;
        OutFile->cd();
        hConfLev[uSyst]->Write();
        hConfidence[uSyst]->Write();
    }

    OutFile->cd();
    TH2F* hFinalEP_Cdf = new TH2F("hFinalEP_Cdf","hFinalEP_Cdf",Num_f0_bins,f0_min,f0_max,Num_d0_bins,d0_min,d0_max);
    TH2F* hFinalEP_Confidence = new TH2F("hFinalEP_Confidence","hFinalEP_Confidence",Num_f0_bins,f0_min,f0_max,Num_d0_bins,d0_min,d0_max);
    for(unsigned uBinX=0; uBinX<Num_f0_bins; uBinX++){
        for(unsigned uBinY=0; uBinY<Num_d0_bins; uBinY++){
            hFinalEP_Cdf->SetBinContent(uBinX+1,uBinY+1,1e6);
            hFinalEP_Confidence->SetBinContent(uBinX+1,uBinY+1,10);
        }
    }

    for(unsigned uSyst=0; uSyst<NumSyst; uSyst++){
        for(unsigned uBinX=0; uBinX<Num_f0_bins; uBinX++){
            for(unsigned uBinY=0; uBinY<Num_d0_bins; uBinY++){
                if(hFinalEP_Cdf->GetBinContent(uBinX+1,uBinY+1)>hConfLev[uSyst]->GetBinContent(uBinX+1,uBinY+1)){
                    PVAL = 1.-hConfLev[uSyst]->GetBinContent(uBinX+1,uBinY+1);
                    if(PVAL>0 && PVAL<1)NSIGMA = sqrt(2)*TMath::ErfcInverse(PVAL);
                    else if(PVAL<=0) NSIGMA = 10;
                    else NSIGMA = 0;
                    if(NSIGMA>10 || NSIGMA!=NSIGMA) NSIGMA = 10;
                    //if(NSIGMA==0){printf("PVAL=%f\n",PVAL);}
                    hFinalEP_Cdf->SetBinContent(uBinX+1,uBinY+1,hConfLev[uSyst]->GetBinContent(uBinX+1,uBinY+1));
                    hFinalEP_Confidence->SetBinContent(uBinX+1,uBinY+1,NSIGMA);
                    //sqrt(2)*TMath::ErfcInverse();
                }
            }
        }
    }
    hFinalEP_Cdf->Write();
    hFinalEP_Confidence->Write();

    for(unsigned uSyst=0; uSyst<NumSyst; uSyst++){
        delete hConfLev[uSyst];
        delete hConfidence[uSyst];
    }
    delete [] hConfLev;
    delete [] hConfidence;
    delete hFinalEP_Cdf;
    delete hFinalEP_Confidence;
    delete OutFile;
}

void ComputeBindingEnergy(const TString FolderName){

    const unsigned NumBins_Eb=128;
    const unsigned NumBins_Eb_Y=128;
    const double PlotMax_d0=4.5;
    const double PlotMax_f0Inv=0.6;
    const double PlotMax_Eb=9;
    TH2F* hEbin_d0 = new TH2F("hEbin_d0","hEbin_d0",NumBins_Eb,0,PlotMax_Eb,NumBins_Eb_Y,0,PlotMax_d0);
    TH2F* hEbin_f0Inv = new TH2F("hEbin_f0Inv","hEbin_f0Inv",NumBins_Eb,0,PlotMax_Eb,NumBins_Eb_Y,-PlotMax_f0Inv,0);
    for(unsigned uBinX=0; uBinX<NumBins_Eb; uBinX++){
        for(unsigned uBinY=0; uBinY<NumBins_Eb_Y; uBinY++){
            hEbin_d0->SetBinContent(uBinX+1,uBinY+1,0);
            hEbin_f0Inv->SetBinContent(uBinX+1,uBinY+1,0);
        }
    }

    double BinEn_Min=1e6;
    double BinEn_Max=0;
    double BinEn_Best=0;
    double d0_Best=0;
    double f0Inv_Best=0;
    double Pval_Best=0;

    TFile* InFile = new TFile(FolderName+"FinalEP.root","read");
    TH2F* hExclPlot = (TH2F*)InFile->Get("hFinalEP_Cdf");
    //TH2F* hExclPlot = (TH2F*)InFile->Get("hConfLev_0");

    unsigned NumBins_f0Inv = hExclPlot->GetXaxis()->GetNbins();
    unsigned NumBins_d0 = hExclPlot->GetYaxis()->GetNbins();

    double Pval;
    double f0Inv_Val;
    double d0_Val;
    double kappa, Ebinding;
    const double hbarc=197.327;
    const double MassLambda = 1115.683;
    for(unsigned uBinX=0; uBinX<NumBins_f0Inv; uBinX++){
        for(unsigned uBinY=0; uBinY<NumBins_d0; uBinY++){
            Pval = 1.-hExclPlot->GetBinContent(uBinX+1, uBinY+1);
            if(Pval<0.3173) continue;
            f0Inv_Val = hExclPlot->GetXaxis()->GetBinCenter(uBinX+1);
            d0_Val = hExclPlot->GetYaxis()->GetBinCenter(uBinY+1);
            if(f0Inv_Val<-0.7 || f0Inv_Val>=0 || d0_Val>4) continue;
            //printf("f0Inv,d0 = %.2f, %.2f; sqrt = %.2f\n",f0Inv_Val,d0_Val,2.*d0_Val*f0Inv_Val+1);
            if( 2.*d0_Val*f0Inv_Val+1<0 ) continue;
            //printf("f0Inv,d0 = %.2f, %.2f; BE = %.2f\n",f0Inv_Val,d0_Val,Ebinding);

            kappa=(1.-sqrt(1.+2.*d0_Val*f0Inv_Val))/d0_Val*hbarc;
            Ebinding=kappa*kappa/MassLambda;


            if(Pval>Pval_Best) {
                Pval_Best=Pval; BinEn_Best=Ebinding; d0_Best=d0_Val; f0Inv_Best=f0Inv_Val;
            }
            if(Ebinding<BinEn_Min) {BinEn_Min=Ebinding;}
            if(Ebinding>BinEn_Max) {
                printf("f0Inv,d0 = %.2f, %.2f; BE = %.2f\n",f0Inv_Val,d0_Val,Ebinding);
                BinEn_Max=Ebinding;
            }

            hEbin_d0->SetBinContent(hEbin_d0->FindBin(Ebinding,d0_Val),1);
            hEbin_d0->SetBinError(hEbin_d0->FindBin(Ebinding,d0_Val),0);
            hEbin_f0Inv->SetBinContent(hEbin_f0Inv->FindBin(Ebinding,f0Inv_Val),1);
            hEbin_f0Inv->SetBinError(hEbin_f0Inv->FindBin(Ebinding,f0Inv_Val),0);
        }
    }

    //fill up empty gaps in the histos, by seeing the first and last filled bins along x and y axis
    unsigned FirstBin;
    unsigned LastBin;
    for(unsigned uBinX=1; uBinX<=hEbin_d0->GetNbinsX(); uBinX++){
        FirstBin=0;
        LastBin=0;
        for(unsigned uBinY=1; uBinY<=hEbin_d0->GetNbinsY(); uBinY++){
            if(!FirstBin && hEbin_d0->GetBinContent(uBinX,uBinY)){
                FirstBin = uBinY;
            }
            if(hEbin_d0->GetBinContent(uBinX,uBinY) && LastBin<uBinY){
                LastBin=uBinY;
            }
        }
        for(unsigned uBinY=FirstBin; uBinY<=LastBin; uBinY++){
            hEbin_d0->SetBinContent(uBinX,uBinY,1);
        }
    }
    for(unsigned uBinY=1; uBinY<=hEbin_d0->GetNbinsY(); uBinY++){
        FirstBin=0;
        LastBin=0;
        for(unsigned uBinX=1; uBinX<=hEbin_d0->GetNbinsX(); uBinX++){
            if(!FirstBin && hEbin_d0->GetBinContent(uBinX,uBinY)){
                FirstBin = uBinX;
            }
            if(hEbin_d0->GetBinContent(uBinX,uBinY) && LastBin<uBinX){
                LastBin=uBinX;
            }
        }
        for(unsigned uBinX=FirstBin; uBinX<=LastBin; uBinX++){
            hEbin_d0->SetBinContent(uBinX,uBinY,1);
        }
    }

    for(unsigned uBinX=1; uBinX<=hEbin_f0Inv->GetNbinsX(); uBinX++){
        FirstBin=0;
        LastBin=0;
        for(unsigned uBinY=1; uBinY<=hEbin_f0Inv->GetNbinsY(); uBinY++){
            if(!FirstBin && hEbin_f0Inv->GetBinContent(uBinX,uBinY)){
                FirstBin = uBinY;
            }
            if(hEbin_f0Inv->GetBinContent(uBinX,uBinY) && LastBin<uBinY){
                LastBin=uBinY;
            }
        }
        for(unsigned uBinY=FirstBin; uBinY<=LastBin; uBinY++){
            hEbin_f0Inv->SetBinContent(uBinX,uBinY,1);
        }
    }
    for(unsigned uBinY=1; uBinY<=hEbin_f0Inv->GetNbinsY(); uBinY++){
        FirstBin=0;
        LastBin=0;
        for(unsigned uBinX=1; uBinX<=hEbin_f0Inv->GetNbinsX(); uBinX++){
            if(!FirstBin && hEbin_f0Inv->GetBinContent(uBinX,uBinY)){
                FirstBin = uBinX;
            }
            if(hEbin_f0Inv->GetBinContent(uBinX,uBinY) && LastBin<uBinX){
                LastBin=uBinX;
            }
        }
        for(unsigned uBinX=FirstBin; uBinX<=LastBin; uBinX++){
            hEbin_f0Inv->SetBinContent(uBinX,uBinY,1);
        }
    }

    printf("Pval_Best = %.2f\n",Pval_Best);
    printf("BinEn_Best = %.2f\n",BinEn_Best);
    printf("d0_Best = %.2f\n",d0_Best);
    printf("f0Inv_Best = %.2f\n",f0Inv_Best);

    printf("BinEn_Min = %.2f\n",BinEn_Min);
    printf("BinEn_Max = %.2f\n",BinEn_Max);

    TGraph grBestFit_d0(1);
    grBestFit_d0.SetName("grBestFit_d0");
    grBestFit_d0.SetPoint(0,BinEn_Best,d0_Best);

    TGraph grBestFit_f0Inv(1);
    grBestFit_f0Inv.SetName("grBestFit_f0Inv");
    grBestFit_f0Inv.SetPoint(0,BinEn_Best,f0Inv_Best);

    TFile* OutFile = new TFile(FolderName+"BindingEnergy.root","recreate");

    //const unsigned nContours=4;
    //double contours[nContours] = {0,0.75,1.5,2.25};
    //hEbin_d0->SetContour(nContours,contours);

    hEbin_d0->Write();
    hEbin_f0Inv->Write();
    grBestFit_d0.Write();
    grBestFit_f0Inv.Write();

    delete hEbin_d0;
    delete hEbin_f0Inv;
    delete InFile;
    delete OutFile;
}


Double_t ChiSquareDistr(Double_t *x,Double_t *par)
{
// Chisquare density distribution for nrFree degrees of freedom

Double_t nrFree = par[0];
Double_t chi2 = x[0]*nrFree;

if (chi2 > 0) {
Double_t lambda = nrFree/2.;
Double_t norm = TMath::Gamma(lambda)*TMath::Power(2.,lambda)/nrFree;
return TMath::Power(chi2,lambda-1)*TMath::Exp(-0.5*chi2)/norm;
} else
return 0.0;
}

void ShowChi2Plots(const TString InputFolderName, const TString InputFileName){

    double MIN_CHI2=0;
    const double Chi2From=0;
    const double Chi2To=22;
    const unsigned NumBins = 128;

    TH1F* hDeltaChi2 = new TH1F("hDeltaChi2","hDeltaChi2",NumBins,Chi2From,Chi2To);
    TH1F* hDeltaChi2_CDF = new TH1F("hDeltaChi2_CDF","hDeltaChi2_CDF",NumBins,Chi2From,Chi2To);
    TGraph gDeltaChi2_CDF;
    gDeltaChi2_CDF.Set(NumBins);
    gDeltaChi2_CDF.SetLineWidth(4);
    gDeltaChi2_CDF.SetFillColor(kWhite);
    gDeltaChi2_CDF.SetMarkerSize(0.1);

    TFile* InputFile = new TFile(InputFolderName+InputFileName,"read");
    TNtuple* ntParDensity = (TNtuple*)InputFile->Get("ntParDensity");
    Float_t dChi2Ndf;
    Float_t f0Inv;
    Float_t d0_val;
    ntParDensity->SetBranchAddress("dChi2",&dChi2Ndf);
    ntParDensity->SetBranchAddress("f0Inv",&f0Inv);
    ntParDensity->SetBranchAddress("d0",&d0_val);
    unsigned NumEntries = ntParDensity->GetEntries();
    for(unsigned uEntry=0; uEntry<NumEntries; uEntry++){
        ntParDensity->GetEntry(uEntry);
        if(dChi2Ndf<MIN_CHI2) MIN_CHI2=dChi2Ndf;
    }
    unsigned ActualEntries=0;
    for(unsigned uEntry=0; uEntry<NumEntries; uEntry++){
        ntParDensity->GetEntry(uEntry);
        //bring the starting point to zero
        dChi2Ndf += MIN_CHI2;
        //SOME CONDITION??

//if( f0Inv>0 ) continue;

        //if( f0Inv<-0.799 ) continue;
        //if( f0Inv>-0.101 && f0Inv<0 ) continue;
        //if( f0Inv<0.201 && f0Inv>0 ) continue;
        //if( f0Inv>31.999 ) continue;

        //if( d0_val<0.01001 ) continue;
        //if( d0_val>63.999 ) continue;

        ActualEntries++;
        hDeltaChi2->Fill(dChi2Ndf);
    }
    hDeltaChi2->Scale(1./double(ActualEntries),"width");
    double CDF=0;
    for(unsigned uBin=1; uBin<=NumBins; uBin++){
        CDF += hDeltaChi2->GetBinContent(uBin)*hDeltaChi2->GetBinWidth(uBin);
        hDeltaChi2_CDF->SetBinContent(uBin,1.-CDF);
        gDeltaChi2_CDF.SetPoint(uBin,hDeltaChi2_CDF->GetBinCenter(uBin),1.-CDF);
    }

    TF1* fChi2Ndf = new TF1("fChi2Ndf",ChiSquareDistr,0.001,10,1);
    fChi2Ndf->FixParameter(0,2);
    printf("Integral = %.3f\n",fChi2Ndf->Integral(0.001,10));
//FitKitty = new TF1("FitKitty", Dlm_ppFitter1, hOliBinShift->GetBinCenter(DataSet==2?3:1), hOliBinShift->GetBinCenter(NumMomBins),4);
    TCanvas* cDeltaChi2 = new TCanvas(TString::Format("cDeltaChi2"));
    cDeltaChi2->SetCanvasSize(1920, 1080); cDeltaChi2->SetMargin(0.14,0.05,0.2,0.05);//lrbt
    hDeltaChi2->Draw();
    fChi2Ndf->Draw("same");
    cDeltaChi2->SaveAs(InputFolderName+"cDeltaChi2.png");

    TF1* fChi2Ndf_CDF = new TF1("fChi2Ndf_CDF","TMath::Prob(x,[0])",0,Chi2To);
    fChi2Ndf_CDF->FixParameter(0,fChi2Ndf->GetParameter(0));
    fChi2Ndf_CDF->SetLineWidth(4);
    fChi2Ndf_CDF->SetLineColor(kRed);
    fChi2Ndf_CDF->SetLineStyle(5);

    TF1* f1Sigma = new TF1("f1Sigma","[0]",0,Chi2To);
    TF1* f2Sigma = new TF1("f2Sigma","[0]",0,Chi2To);
    TF1* f3Sigma = new TF1("f3Sigma","[0]",0,Chi2To);
    TF1* f4Sigma = new TF1("f4Sigma","[0]",0,Chi2To);

    f1Sigma->FixParameter(0, 1.-0.682689492137086);
    f1Sigma->SetLineColor(kGray);
    f1Sigma->SetLineWidth(2);
    f1Sigma->SetLineStyle(2);
    f2Sigma->FixParameter(0, 1.-0.954499736103642);
    f2Sigma->SetLineColor(kGray);
    f2Sigma->SetLineWidth(2);
    f2Sigma->SetLineStyle(2);
    f3Sigma->FixParameter(0, 1.-0.997300203936740);
    f3Sigma->SetLineColor(kGray);
    f3Sigma->SetLineWidth(2);
    f3Sigma->SetLineStyle(2);
    f4Sigma->FixParameter(0, 1.-0.999936657516334);
    f4Sigma->SetLineColor(kGray);
    f4Sigma->SetLineWidth(2);
    f4Sigma->SetLineStyle(2);


    hDeltaChi2_CDF->SetStats(false);
    hDeltaChi2_CDF->SetTitle("");
    hDeltaChi2_CDF->GetXaxis()->SetLabelSize(0.065);
    hDeltaChi2_CDF->GetXaxis()->SetTitle("#Delta#chi^{2}");
    hDeltaChi2_CDF->GetXaxis()->CenterTitle();
    hDeltaChi2_CDF->GetXaxis()->SetTitleOffset(1.15);
    hDeltaChi2_CDF->GetXaxis()->SetLabelOffset(0.02);
    hDeltaChi2_CDF->GetXaxis()->SetTitleSize(0.075);
    hDeltaChi2_CDF->GetYaxis()->SetLabelSize(0.065);
    hDeltaChi2_CDF->GetYaxis()->SetTitle("p-val");
    hDeltaChi2_CDF->GetYaxis()->CenterTitle();
    hDeltaChi2_CDF->GetYaxis()->SetTitleOffset(0.8);
    hDeltaChi2_CDF->GetYaxis()->SetTitleSize(0.075);

    hDeltaChi2_CDF->GetYaxis()->SetRangeUser(3e-5,2);

//FitKitty = new TF1("FitKitty", Dlm_ppFitter1, hOliBinShift->GetBinCenter(DataSet==2?3:1), hOliBinShift->GetBinCenter(NumMomBins),4);
    TCanvas* cDeltaChi2_CDF = new TCanvas(TString::Format("cDeltaChi2_CDF"));
    cDeltaChi2_CDF->SetCanvasSize(1920/3, 1080/3); cDeltaChi2_CDF->SetMargin(0.14,0.05,0.2,0.05);//lrbt

    hDeltaChi2_CDF->Draw("axis");
    gDeltaChi2_CDF.Draw("same,C");
    fChi2Ndf_CDF->Draw("same");
    f1Sigma->Draw("same");
    f2Sigma->Draw("same");
    f3Sigma->Draw("same");
    f4Sigma->Draw("same");
    cDeltaChi2_CDF->SetGrid();
    cDeltaChi2_CDF->SetLogy();

    TLatex TextSigma;
    TextSigma.SetTextSize(gStyle->GetTextSize()*1);
    TextSigma.SetNDC(kTRUE);
    TextSigma.DrawLatex(0.45, 0.835, "1#sigma");
    TextSigma.DrawLatex(0.45, 0.705, "2#sigma");
    TextSigma.DrawLatex(0.45, 0.515, "3#sigma");
    TextSigma.DrawLatex(0.45, 0.260, "4#sigma");
    TLegend* myLegend = new TLegend(0.6,0.75,0.95,0.95);//lbrt
    myLegend->SetName("myLegend");
    myLegend->SetTextSize(0.05);
    myLegend->AddEntry(&gDeltaChi2_CDF, "Bootstrap");
    myLegend->AddEntry(fChi2Ndf_CDF, "Normal uncertainties");
    myLegend->Draw("same");
    cDeltaChi2_CDF->SaveAs(InputFolderName+"cDeltaChi2_CDF.png");
    cDeltaChi2_CDF->SaveAs(InputFolderName+"cDeltaChi2_CDF.pdf");


    delete hDeltaChi2;
    delete hDeltaChi2_CDF;
    delete fChi2Ndf;
    delete fChi2Ndf_CDF;
    delete f1Sigma;
    delete f2Sigma;
    delete f3Sigma;
    delete f4Sigma;
    delete myLegend;
    delete cDeltaChi2;
    delete cDeltaChi2_CDF;
    delete InputFile;

}
int LamLamMcExclusionMAIN(int narg, char** ARGS){
//1 = Computer
//2 = OutputFolder
//3 = Chi2Map (0 no, +1 individual, -1 common for all, -2 common per systemtics) (1024x768)
//              +10 => 512x384
//              +20 => 256x192
//              +30 => 128x96
//              +40 => 64x48
//4 = Bootstrap
//5 = NumBootIter
//6 = TimeLimit
//7 = JOBID
//8 = SystID
//9 = RandSyst
//10 = SystJobID
//11 = SepBL
//12 = FitFun
//13 = BinWidthSyst
//14 = RadiusSyst
//15 = FitSyst
//16 = DataSets

    const TString Computer = "XMG";
    const int Chi2Map = -12;
    const bool Bootstrap = false;
    const int NumBootIter = 2;
    const double TimeLimit = 9000;
    const int JobId = 0;
    const int SystId = 0;
    const bool RandSyst = false;
    const int BootId = 0;
    const int SepBL = 0;//-1 is no BL
    const int FitFun = 7;

    LamLamConfidence_22Oct2018(Computer.Data(),"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/LamLam/EB_Answers/pp13TeV/",Chi2Map,Bootstrap,NumBootIter,TimeLimit,0,0,false,0,SepBL,FitFun,0,2,0,"pp13TeV");
    LamLamConfidence_22Oct2018(Computer.Data(),"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/LamLam/EB_Answers/pPb5TeV/",Chi2Map,Bootstrap,NumBootIter,TimeLimit,0,0,false,0,SepBL,FitFun,0,2,0,"pPb5TeV");
    LamLamConfidence_22Oct2018(Computer.Data(),"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/LamLam/EB_Answers/Full/",Chi2Map,Bootstrap,NumBootIter,TimeLimit,0,0,false,0,SepBL,FitFun,0,2,0,"pp13TeV,pPb5TeV,pp7TeV");
    //for(unsigned uSystId=0; uSystId<27; uSystId++){
        //printf("uSystId=%u\n",uSystId);
        //LamLamConfidence_22Oct2018(Computer.Data(),"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Temp/",Chi2Map,Bootstrap,NumBootIter,TimeLimit,0,uSystId,false,uSystId,SepBL,FitFun,0,2,0,"pp13TeV,pPb5TeV,pp7TeV");
    //}

    //FinalExclusionPlot(Computer,-2,"/home/dmihaylov/Temp/Output/170119/MC_Default_BL_AllDataSets/",1);
    //FinalExclusionPlot(Computer,-2,"/home/dmihaylov/Temp/Output/170119/MC_FullSyst_BL_AllDataSets/",121);
    //ComputeBindingEnergy("/home/dmihaylov/Temp/Output/170119/MC_Default_BL_AllDataSets/");
    //ComputeBindingEnergy("/home/dmihaylov/Temp/Output/170119/MC_FullSyst_BL_AllDataSets/");

    return 0;
}
