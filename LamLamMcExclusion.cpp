
#include "LamLamMcExclusion.h"

#include <iostream>
#include <unistd.h>
#include <fcntl.h>
#include <cerrno>
#include <cstdio>
#include <string.h>
//#include <omp.h>

#include "DLM_Fitters.h"
#include "DLM_CkDecomposition.h"
#include "DLM_CkModels.h"
#include "DLM_CppTools.h"
#include "TFile.h"
#include "TH1F.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TMath.h"

using namespace std;

LamLamMcExclusion::LamLamMcExclusion(){

    OutputFile = NULL;
    OutputFileName = NULL;
    SaveChi2Map = true;
    SaveFitResults = false;
    Chi2MapFileName = NULL;
    NumBins_f0Inv=256;
    Min_f0Inv=-2;
    Max_f0Inv=-5;
    NumBins_d0=196;
    Min_d0=0;
    Max_d0=18;
    RanGenSyst = new TRandom3();
    RanGenMC = new TRandom3();
    NumJobs=1;
    //JobID=0;
    //FixSeed=false;
    NumDataSets=0;
    NumDataVar=NULL;//[NumDataSets]
    NumFitFunctions=NULL;//[NumDataSets]
    FitFunction=NULL;//[NumDataSets][NumFitFunctions]
    SeparateBaseline = NULL;//[NumDataSets]
    NumFemtoReg = NULL;//[NumDataSets][NumDataVar]
    FemtoRegion = NULL;//[NumDataSets][NumDataVar][NumFemtoRegMin][min/max]
    NumBlReg = NULL;//[NumDataSets][NumDataVar]
    BlRegion = NULL;//[NumDataSets][NumDataVar][NumFemtoRegMin][min/max]

    hData = NULL;//[NumDataSets][NumDataVar]
    hResMatrix = NULL;//[NumDataSets][NumDataVar]
    NumLamPar = NULL;//[NumDataSets][NumDataVar]
    LambdaPar0 = NULL;//[NumDataSets][NumDataVar][NumLamPar]
    LambdaParF = NULL;//[NumDataSets][NumDataVar][NumLamPar]

    NumRadVar = NULL;//[NumDataSets]
    Radius = NULL;//[NumDataSets][NumRadVar]

    Bootstrap = true;
    NumSystVariations = 0;
    SpecificSyst = false;
    FirstSyst = 0;
    LastSyst = unsigned(-1);
    RunDefaultSyst = false;
    NumBootIter = 0;
    TimeLimit = 1e64;
    TIMER = new DLM_Timer();

    ntParDensity = NULL;
    SamplingHisto = NULL;//[NumDataSets]
    SamplingGraph = NULL;//[NumDataSets]
    Chi2min = 0;//[NumDataSets]
    NDFmin = 0;
    f0Inv_min = 0;
    d0_min = 0;

}

LamLamMcExclusion::~LamLamMcExclusion(){

    delete RanGenSyst; RanGenSyst=NULL;
    delete RanGenMC; RanGenMC=NULL;
    ClearMem();
    if(ntParDensity){delete ntParDensity; ntParDensity=NULL;}
    if(OutputFile){delete OutputFile; OutputFile=NULL;}
    if(OutputFileName){delete[]OutputFileName; OutputFileName=NULL;}
    if(Chi2MapFileName){delete[]Chi2MapFileName; Chi2MapFileName=NULL;}
    delete TIMER; TIMER=NULL;

}

//clears memory related to each data set
//if OnlyDataVar we only delete entries, which depend on the data variations.
//the last case is needed in case we do not change the data sets, only the data variations
void LamLamMcExclusion::ClearMemDataSet(const unsigned& uData, const bool& OnlyDataVar){

    if(uData>=NumDataSets){printf("\033[1;33mWARNING:\033[0m Attempting to access uData=%u, but NumDataSets=%u!\n",uData,NumDataSets);return;}
    if(!OnlyDataVar){
        if(NumFitFunctions){NumFitFunctions[uData]=0;}
        if(FitFunction&&FitFunction[uData]){delete[]FitFunction[uData];FitFunction[uData]=NULL;}
    }

    if(FemtoRegion&&FemtoRegion[uData]){
        for(unsigned uDataVar=0; uDataVar<NumDataVar[uData]; uDataVar++){
            for(unsigned uVar=0; uVar<NumFemtoReg[uData][uDataVar]; uVar++){
                delete[]FemtoRegion[uData][uDataVar][uVar];FemtoRegion[uData][uDataVar][uVar]=NULL;
            }
            delete[]FemtoRegion[uData][uDataVar];FemtoRegion[uData][uDataVar]=NULL;
        }
        delete[]FemtoRegion[uData];FemtoRegion[uData]=NULL;
    }
    if(NumFemtoReg){delete[]NumFemtoReg[uData];NumFemtoReg[uData]=NULL;}
    if(BlRegion&&BlRegion[uData]){
        for(unsigned uDataVar=0; uDataVar<NumDataVar[uData]; uDataVar++){
            for(unsigned uVar=0; uVar<NumBlReg[uData][uDataVar]; uVar++){
                delete[]BlRegion[uData][uDataVar][uVar];BlRegion[uData][uDataVar][uVar]=NULL;
            }
            delete[]BlRegion[uData][uDataVar];BlRegion[uData][uDataVar]=NULL;
        }
        delete[]BlRegion[uData];BlRegion[uData]=NULL;
    }
    if(NumBlReg){delete[]NumBlReg[uData];NumBlReg[uData]=NULL;}
    if(hData&&hData[uData]){delete[]hData[uData];hData[uData]=NULL;}
    if(hResMatrix&&hResMatrix[uData]){delete[]hResMatrix[uData];hResMatrix[uData]=NULL;}
    if(NumLamPar&&NumLamPar[uData]){
        delete[]NumLamPar[uData];NumLamPar[uData]=NULL;
    }
    if(LambdaPar0&&LambdaPar0[uData]){
        for(unsigned uVar=0; uVar<NumDataVar[uData]; uVar++){delete[]LambdaPar0[uData][uVar];LambdaPar0[uData][uVar]=NULL;}
        delete[]LambdaPar0[uData];LambdaPar0[uData]=NULL;
    }
    if(LambdaParF&&LambdaParF[uData]){
        for(unsigned uVar=0; uVar<NumDataVar[uData]; uVar++){delete[]LambdaParF[uData][uVar];LambdaParF[uData][uVar]=NULL;}
        delete[]LambdaParF[uData];LambdaParF[uData]=NULL;
    }

    if(!OnlyDataVar){
        if(NumRadVar){NumRadVar[uData]=0;}
        if(Radius&&Radius[uData]){delete[]Radius[uData];Radius[uData]=NULL;}
        if(SamplingHisto&&SamplingHisto[uData]){
            delete SamplingHisto[uData];SamplingHisto[uData]=NULL;
        }
        if(NumDataVar){NumDataVar[uData]=0;}
    }

}

void LamLamMcExclusion::ClearMem(){

    for(unsigned uData=0; uData<NumDataSets; uData++){
        ClearMemDataSet(uData,false);
    }
    NumDataSets=0;
    if(NumDataVar){delete[]NumDataVar;NumDataVar=NULL;}
    if(NumFitFunctions){delete[]NumFitFunctions;NumFitFunctions=NULL;}
    if(FitFunction){delete[]FitFunction;FitFunction=NULL;}
    if(SeparateBaseline){delete[]SeparateBaseline;SeparateBaseline=NULL;}
    if(NumFemtoReg){delete[]NumFemtoReg;NumFemtoReg=NULL;}
    if(FemtoRegion){delete[]FemtoRegion;FemtoRegion=NULL;}
    if(NumBlReg){delete[]NumBlReg;NumBlReg=NULL;}
    if(BlRegion){delete[]BlRegion;BlRegion=NULL;}
    if(hData){delete[]hData;hData=NULL;}
    if(hResMatrix){delete[]hResMatrix;hResMatrix=NULL;}
    if(NumLamPar){delete[]NumLamPar;NumLamPar=NULL;}
    if(LambdaPar0){delete[]LambdaPar0;LambdaPar0=NULL;}
    if(LambdaParF){delete[]LambdaParF;LambdaParF=NULL;}
    if(NumRadVar){delete[]NumRadVar;NumRadVar=NULL;}
    if(Radius){delete[]Radius;Radius=NULL;}
    if(SamplingHisto){delete[]SamplingHisto;SamplingHisto=NULL;}
    if(SamplingGraph){delete[]SamplingGraph;SamplingGraph=NULL;}

}

void LamLamMcExclusion::SetOutputFile(const char* outfilename){

    if(!outfilename){
        if(OutputFileName){delete[]OutputFileName;OutputFileName=NULL;}
        if(OutputFile){delete OutputFile;OutputFile=NULL;}
        return;
    }
    if(OutputFile && strcmp(OutputFileName,outfilename)==0) return;
    if(OutputFileName) delete [] OutputFileName;
    OutputFileName = new char [int(strlen(outfilename))+1];
    strcpy(OutputFileName,outfilename);
    if(OutputFile){delete OutputFile;}
    OutputFile = new TFile(OutputFileName,"recreate");

}
void LamLamMcExclusion::SetSaveChi2Map(const bool& savechi2map){

    SaveChi2Map = savechi2map;

}
void LamLamMcExclusion::SetSaveFitResults(const bool& savefitresults){

    SaveFitResults = savefitresults;

}
void LamLamMcExclusion::SetUpChi2Map(const unsigned& nf0inv, const double& f0invmin, const double& f0invmax,
                                     const unsigned& nd0, const double& d0min, const double& d0max){

    NumBins_f0Inv=nf0inv;
    Min_f0Inv=f0invmin;
    Max_f0Inv=f0invmax;
    NumBins_d0=nd0;
    Min_d0=d0min;
    Max_d0=d0max;

}
void LamLamMcExclusion::SetCommonChi2MapFile(const char* chi2filename){

    if(!chi2filename){
        if(Chi2MapFileName){delete[]Chi2MapFileName;Chi2MapFileName=NULL;}
        return;
    }
    if(Chi2MapFileName) delete [] Chi2MapFileName;
    Chi2MapFileName = new char [int(strlen(chi2filename))+1];
    strcpy(Chi2MapFileName,chi2filename);

}
void LamLamMcExclusion::SetRandomSeed(const unsigned& rs_syst, const unsigned& rs_mc){

    delete RanGenSyst;
    RanGenSyst = new TRandom3(rs_syst);

    delete RanGenMC;
    RanGenMC = new TRandom3(rs_mc);

}
void LamLamMcExclusion::SetNumDataSets(const unsigned& NumDS){

    if(NumDS==NumDataSets) return;
    ClearMem();
    NumDataSets = NumDS;
    NumDataVar = new unsigned [NumDataSets];
    NumFitFunctions = new unsigned [NumDataSets];
    FitFunction = new int* [NumDataSets];
    SeparateBaseline = new int [NumDataSets];
    NumFemtoReg = new unsigned* [NumDataSets];
    FemtoRegion = new double*** [NumDataSets];
    NumBlReg = new unsigned* [NumDataSets];
    BlRegion = new double*** [NumDataSets];
    hData = new TH1F** [NumDataSets];
    hResMatrix = new TH2F** [NumDataSets];
    NumLamPar = new unsigned* [NumDataSets];
    LambdaPar0 = new double** [NumDataSets];
    LambdaParF = new double** [NumDataSets];
    NumRadVar = new unsigned [NumDataSets];
    Radius = new double* [NumDataSets];
    SamplingHisto = new TH1F* [NumDataSets];
    SamplingGraph = new TGraph [NumDataSets];
    for(unsigned uData=0; uData<NumDataSets; uData++){
        NumDataVar[uData]=0;
        NumFitFunctions[uData]=0;
        FitFunction[uData]=NULL;
        SeparateBaseline[uData]=false;
        NumFemtoReg[uData]=NULL;
        FemtoRegion[uData]=NULL;
        NumBlReg[uData]=NULL;
        BlRegion[uData]=NULL;
        hData[uData]=NULL;
        hResMatrix[uData]=NULL;
        NumLamPar[uData]=NULL;
        LambdaPar0[uData]=NULL;
        LambdaParF[uData]=NULL;
        NumRadVar[uData]=0;
        Radius[uData]=NULL;
        SamplingHisto[uData]=NULL;
    }

}
void LamLamMcExclusion::SetNumDataVar(const unsigned& uData, const unsigned& val){

    if(uData>=NumDataSets){printf("\033[1;33mWARNING:\033[0m Attempting to access uData=%u, but NumDataSets=%u!\n",uData,NumDataSets);return;}
    if(NumDataVar[uData]==val) return;
    ClearMemDataSet(uData,true);
    NumDataVar[uData] = val;

    NumFemtoReg[uData] = new unsigned [NumDataVar[uData]];
    FemtoRegion[uData] = new double** [NumDataVar[uData]];
    for(unsigned uDataVar=0; uDataVar<NumDataVar[uData]; uDataVar++){
        FemtoRegion[uData][uDataVar] = NULL;
    }

    NumBlReg[uData] = new unsigned [NumDataVar[uData]];
    BlRegion[uData] = new double** [NumDataVar[uData]];
    for(unsigned uDataVar=0; uDataVar<NumDataVar[uData]; uDataVar++){
        BlRegion[uData][uDataVar] = NULL;
    }

    hData[uData] = new TH1F* [NumDataVar[uData]];

    hResMatrix[uData] = new TH2F* [NumDataVar[uData]];
    NumLamPar[uData] = new unsigned [NumDataVar[uData]];
    LambdaPar0[uData] = new double* [NumDataVar[uData]];
    LambdaParF[uData] = new double* [NumDataVar[uData]];
    for(unsigned uDataVar=0; uDataVar<NumDataVar[uData]; uDataVar++){
        hData[uData][uDataVar]=NULL;
        hResMatrix[uData][uDataVar]=NULL;
        NumLamPar[uData][uDataVar] = 0;
		LambdaPar0[uData][uDataVar]=NULL;
	    LambdaParF[uData][uDataVar]=NULL;
	}

}

void LamLamMcExclusion::SetNumFitFunVar(const unsigned& uData, const unsigned& val){

    if(uData>=NumDataSets){printf("\033[1;33mWARNING:\033[0m Attempting to access uData=%u, but NumDataSets=%u!\n",uData,NumDataSets);return;}
    if(NumFitFunctions[uData]==val) return;
    NumFitFunctions[uData] = val;

    if(FitFunction[uData]){delete[]FitFunction;}
    FitFunction[uData] = new int [NumFitFunctions[uData]];
    for(unsigned uFit=0; uFit<NumFitFunctions[uData]; uFit++){
        FitFunction[uData][uFit] = 1;
    }

}
void LamLamMcExclusion::SetFitFunction(const unsigned& uData, const unsigned& uVar, const unsigned& val){

    if(uData>=NumDataSets){printf("\033[1;33mWARNING:\033[0m Attempting to access uData=%u, but NumDataSets=%u!\n",uData,NumDataSets);return;}
    if(uVar>=NumFitFunctions[uData]){printf("\033[1;33mWARNING:\033[0m Only %u fit functions are allowed for data set %u!\n",NumFitFunctions[uData],uData);return;}
    FitFunction[uData][uVar] = val;

}
void LamLamMcExclusion::SetSeparateBaseline(const unsigned& uData, const int& val){

    if(uData>=NumDataSets){printf("\033[1;33mWARNING:\033[0m Attempting to access uData=%u, but NumDataSets=%u!\n",uData,NumDataSets);return;}
    if(val<-1 || val>1){printf("\033[1;33mWARNING:\033[0m SeparateBaseline should be in [-1,1]\n");return;}
    SeparateBaseline[uData] = val;

}
void LamLamMcExclusion::SetNumFemtoReg(const unsigned& uData, const unsigned& uDataVar, const unsigned& val){

    if(uData>=NumDataSets){printf("\033[1;33mWARNING:\033[0m Attempting to access uData=%u, but NumDataSets=%u!\n",uData,NumDataSets);return;}
    if(uDataVar>=NumDataVar[uData]){printf("\033[1;33mWARNING:\033[0m Only %u variations are allowed for data set %u!\n",NumDataVar[uData],uData);return;}
    if(NumFemtoReg[uData][uDataVar]==val) return;
    if(FemtoRegion[uData][uDataVar]){
        for(unsigned uFemReg=0; uFemReg<NumFemtoReg[uData][uDataVar]; uFemReg++){
            delete[]FemtoRegion[uData][uDataVar][uFemReg];
            FemtoRegion[uData][uDataVar][uFemReg]=NULL;
        }
        delete[]FemtoRegion[uData][uDataVar];
        FemtoRegion[uData][uDataVar] = NULL;
    }
    NumFemtoReg[uData][uDataVar] = val;
    FemtoRegion[uData][uDataVar] = new double* [NumFemtoReg[uData][uDataVar]];
    for(unsigned uVar=0; uVar<NumFemtoReg[uData][uDataVar]; uVar++){
        FemtoRegion[uData][uDataVar][uVar]=new double[2];
        FemtoRegion[uData][uDataVar][uVar][0] = 0;
        FemtoRegion[uData][uDataVar][uVar][1] = 0;
    }

}
void LamLamMcExclusion::SetFemtoRegion(const unsigned& uData, const unsigned& uDataVar, const unsigned& uVar, const double& MIN, const double& MAX){

    if(uData>=NumDataSets){printf("\033[1;33mWARNING:\033[0m Attempting to access uData=%u, but NumDataSets=%u!\n",uData,NumDataSets);return;}
    if(uDataVar>=NumDataVar[uData]){printf("\033[1;33mWARNING:\033[0m Only %u variations are allowed for data set %u!\n",NumDataVar[uData],uData);return;}
    if(uVar>=NumFemtoReg[uData][uDataVar]){printf("\033[1;33mWARNING:\033[0m Only %u femto regions are allowed for data set %u, variation %u!\n",
                                                NumFemtoReg[uData][uDataVar],uData,uDataVar);return;}
    FemtoRegion[uData][uDataVar][uVar][0] = MIN;
    FemtoRegion[uData][uDataVar][uVar][1] = MAX;

}
void LamLamMcExclusion::SetNumBlReg(const unsigned& uData, const unsigned& uDataVar, const unsigned& val){

    if(uData>=NumDataSets){printf("\033[1;33mWARNING:\033[0m Attempting to access uData=%u, but NumDataSets=%u!\n",uData,NumDataSets);return;}
    if(uDataVar>=NumDataVar[uData]){printf("\033[1;33mWARNING:\033[0m Only %u variations are allowed for data set %u!\n",NumDataVar[uData],uData);return;}
    if(NumBlReg[uData][uDataVar]==val) return;
    if(BlRegion[uData][uDataVar]){
        for(unsigned uBlReg=0; uBlReg<NumBlReg[uData][uDataVar]; uBlReg++){
            delete[]BlRegion[uData][uDataVar][uBlReg];
            BlRegion[uData][uDataVar][uBlReg]=NULL;
        }
        delete[]BlRegion[uData][uDataVar];
        BlRegion[uData][uDataVar] = NULL;
    }
    NumBlReg[uData][uDataVar] = val;
    BlRegion[uData][uDataVar] = new double* [NumBlReg[uData][uDataVar]];
    for(unsigned uVar=0; uVar<NumBlReg[uData][uDataVar]; uVar++){
        BlRegion[uData][uDataVar][uVar]=new double[2];
        BlRegion[uData][uDataVar][uVar][0] = 0;
        BlRegion[uData][uDataVar][uVar][1] = 0;
    }

}
void LamLamMcExclusion::SetBlRegion(const unsigned& uData, const unsigned& uDataVar, const unsigned& uVar, const double& MIN, const double& MAX){

    if(uData>=NumDataSets){printf("\033[1;33mWARNING:\033[0m Attempting to access uData=%u, but NumDataSets=%u!\n",uData,NumDataSets);return;}
    if(uDataVar>=NumDataVar[uData]){printf("\033[1;33mWARNING:\033[0m Only %u variations are allowed for data set %u!\n",NumDataVar[uData],uData);return;}
    if(uVar>=NumBlReg[uData][uDataVar]){printf("\033[1;33mWARNING:\033[0m Only %u baseline regions are allowed for data set %u, variation %u!\n",
                                                NumBlReg[uData][uDataVar],uData,uDataVar);return;}
    BlRegion[uData][uDataVar][uVar][0] = MIN;
    BlRegion[uData][uDataVar][uVar][1] = MAX;

}
void LamLamMcExclusion::SetData(const unsigned& uData, const unsigned& uVar, TH1F* hdata){

    if(uData>=NumDataSets){printf("\033[1;33mWARNING:\033[0m Attempting to access uData=%u, but NumDataSets=%u!\n",uData,NumDataSets);return;}
    if(uVar>=NumDataVar[uData]){printf("\033[1;33mWARNING:\033[0m Only %u variations are allowed for data set %u!\n",NumDataVar[uData],uData);return;}
    hData[uData][uVar] = hdata;

}
void LamLamMcExclusion::SetResolution(const unsigned& uData, const unsigned& uVar, TH2F* hresolution){

    if(uData>=NumDataSets){printf("\033[1;33mWARNING:\033[0m Attempting to access uData=%u, but NumDataSets=%u!\n",uData,NumDataSets);return;}
    if(uVar>=NumDataVar[uData]){printf("\033[1;33mWARNING:\033[0m Only %u variations are allowed for data set %u!\n",NumDataVar[uData],uData);return;}
    hResMatrix[uData][uVar] = hresolution;

}
void LamLamMcExclusion::SetNumLamPar(const unsigned& uData, const unsigned& uDataVar, const unsigned& val){

    if(uData>=NumDataSets){printf("\033[1;33mWARNING:\033[0m Attempting to access uData=%u, but NumDataSets=%u!\n",uData,NumDataSets);return;}
    if(uDataVar>=NumDataVar[uData]){printf("\033[1;33mWARNING:\033[0m Only %u variations are allowed for data set %u!\n",NumDataVar[uData],uData);return;}
    if(NumLamPar[uData][uDataVar]==val) return;
    NumLamPar[uData][uDataVar] = val;
    if(LambdaPar0[uData][uDataVar]){delete[]LambdaPar0[uData][uDataVar];}
    if(LambdaParF[uData][uDataVar]){delete[]LambdaParF[uData][uDataVar];}
    LambdaPar0[uData][uDataVar] = new double [NumLamPar[uData][uDataVar]];
    LambdaParF[uData][uDataVar] = new double [NumLamPar[uData][uDataVar]];
    for(unsigned uLamVar=0; uLamVar<NumLamPar[uData][uDataVar]; uLamVar++){
		LambdaPar0[uData][uDataVar][uLamVar]=1;
	    LambdaParF[uData][uDataVar][uLamVar]=0;
	}

}
void LamLamMcExclusion::SetLambdaPar(const unsigned& uData, const unsigned& uDataVar, const unsigned& uLamVar, const double& lam0, const double& lamF){

    if(uData>=NumDataSets){printf("\033[1;33mWARNING:\033[0m Attempting to access uData=%u, but NumDataSets=%u!\n",uData,NumDataSets);return;}
    if(uDataVar>=NumDataVar[uData]){printf("\033[1;33mWARNING:\033[0m Only %u variations are allowed for data set %u!\n",NumDataVar[uData],uData);return;}
    if(uLamVar>=NumLamPar[uData][uDataVar]){printf("\033[1;33mWARNING:\033[0m Only %u lambda parameters are allowed for data set %u, variation %u!\n",
                                            NumLamPar[uData][uDataVar],uData,uDataVar);return;}
    LambdaPar0[uData][uDataVar][uLamVar] = lam0;
    LambdaParF[uData][uDataVar][uLamVar] = lamF;

}
void LamLamMcExclusion::SetNumRadVar(const unsigned& uData, const unsigned& val){

    if(uData>=NumDataSets){printf("\033[1;33mWARNING:\033[0m Attempting to access uData=%u, but NumDataSets=%u!\n",uData,NumDataSets);return;}
    if(NumRadVar[uData]==val) return;
    if(Radius[uData]){delete[]Radius[uData];}
    NumRadVar[uData] = val;
    Radius[uData] = new double [NumRadVar[uData]];
    for(unsigned uRadVar=0; uRadVar<NumRadVar[uData]; uRadVar++){
        Radius[uData][uRadVar] = 2;
    }

}
void LamLamMcExclusion::SetRadius(const unsigned& uData, const unsigned& uRadVar, const double& val){

    if(uData>=NumDataSets){printf("\033[1;33mWARNING:\033[0m Attempting to access uData=%u, but NumDataSets=%u!\n",uData,NumDataSets);return;}
    if(uRadVar>=NumRadVar[uData]){printf("\033[1;33mWARNING:\033[0m Only %u radii are allowed for data set %u!\n",NumRadVar[uData],uData);return;}
    Radius[uData][uRadVar] = val;

}
void LamLamMcExclusion::SetBootstrap(const bool& val){

    Bootstrap = val;

}
void LamLamMcExclusion::SetNumSystVariations(const unsigned& nrandom){

    NumSystVariations = nrandom;

}
void LamLamMcExclusion::RunSpecificSyst(const bool& yesno, const unsigned& sfirst, const unsigned& slast){

    SpecificSyst = yesno;
    FirstSyst = sfirst;
    LastSyst = slast;

}
void LamLamMcExclusion::SetRunDefaultSyst(const bool& rds){

    RunDefaultSyst = rds;

}
void LamLamMcExclusion::SetNumBootIter(const unsigned& nrandom){

    NumBootIter = nrandom;

}
void LamLamMcExclusion::SetTimeLimit(const double& timelimit){

    TimeLimit = timelimit;

}
unsigned LamLamMcExclusion::GetTotNumIter(){

    unsigned TOTITER=1;
    for(unsigned uData=0; uData<NumDataSets; uData++){
        unsigned dataIter=0;
        for(unsigned uDataVar=0; uDataVar<NumDataVar[uData]; uDataVar++){
            for(unsigned uFitVar=0; uFitVar<NumFitFunctions[uData]; uFitVar++){
                for(unsigned uFemtoVar=0; uFemtoVar<NumFemtoReg[uData][uDataVar]; uFemtoVar++){
                    for(unsigned uBlVar=0; uBlVar<NumBlReg[uData][uDataVar]; uBlVar++){
                        for(unsigned uLamPar=0; uLamPar<NumLamPar[uData][uDataVar]; uLamPar++){
                            for(unsigned uRad=0; uRad<NumRadVar[uData]; uRad++){
                                dataIter++;
                            }
                        }
                    }
                }
            }
        }
        TOTITER*=dataIter;
    }
    return TOTITER;

}

void LamLamMcExclusion::GoBabyGo(){

    TIMER->Start();

    if(ntParDensity){
        delete ntParDensity;
        ntParDensity = NULL;
    }

    //UITERfirst = UITER%100000, i.e. the first 5 digits
    //UITERlast = UITER/100000, i.e. the digits after.
    //done so as to be able to save the full precision on UITER
    //the Chi2fit is the the Chi2 of the fit to the random sample
    //dChi2 is the Chi2-Chi2min, where Chi2 is with respect to the original data, Chi2min is the best fit to the original data


    if(!OutputFile){
        printf("\033[1;33mWARNING:\033[0m There is not output file defined, defining LamLamMcExclusion_OUTPUT.root!\n");
        OutputFile = new TFile("LamLamMcExclusion_OUTPUT.root","recreate");
    }

    OutputFile->cd();
    ntParDensity = new TNtuple("ntParDensity","ntParDensity",
                                "UITERfirst:UITERlast:f0Inv:d0:Chi2fit:NDFfit:Chi2refit:dChi2:dChi2Ndf");


    unsigned* DataVar = new unsigned [NumDataSets];
    unsigned* FitFunVar = new unsigned [NumDataSets];
    unsigned* FemRegVar = new unsigned [NumDataSets];
    unsigned* BlRegVar = new unsigned [NumDataSets];
    unsigned* LamParVar = new unsigned [NumDataSets];
    unsigned* RadVar = new unsigned [NumDataSets];

    if(!NumSystVariations || SpecificSyst){

        long int Counter=0;
        unsigned* CurrentDataVar = new unsigned [NumDataSets];
        unsigned* CurrentFitFunVar = new unsigned [NumDataSets];
        unsigned* CurrentFemRegVar = new unsigned [NumDataSets];
        unsigned* CurrentBlRegVar = new unsigned [NumDataSets];
        unsigned* CurrentLamParVar = new unsigned [NumDataSets];
        unsigned* CurrentRadVar = new unsigned [NumDataSets];
        for(unsigned uData=0; uData<NumDataSets; uData++){
            CurrentDataVar[uData]=0;
            CurrentFitFunVar[uData]=0;
            CurrentFemRegVar[uData]=0;
            CurrentBlRegVar[uData]=0;
            CurrentLamParVar[uData]=0;
            CurrentRadVar[uData]=0;
        }

        bool CONTINUE=true;
        while(CONTINUE){
            bool Updated = false;
            for(unsigned uData=0; uData<NumDataSets; uData++){
                DataVar[uData]=CurrentDataVar[uData];
                FitFunVar[uData]=CurrentFitFunVar[uData];
                FemRegVar[uData]=CurrentFemRegVar[uData];
                BlRegVar[uData]=CurrentBlRegVar[uData];
                LamParVar[uData]=CurrentLamParVar[uData];
                RadVar[uData]=CurrentRadVar[uData];
            }
            //do the computation only in case we are inside the range of desired systematics
            if((Counter)>=(long int)(FirstSyst) && (long int)((Counter)<=LastSyst)){
                PerformSingleFit(-1-Counter,DataVar,FitFunVar,FemRegVar,BlRegVar,LamParVar,RadVar);
                //MapOutExclusionPlot(DataVar,FitFunVar,FemRegVar,BlRegVar,LamParVar,RadVar);
                //#pragma omp parallel
                for(int iBoot=0; iBoot<NumBootIter; iBoot++){
                    PerformSingleFit(Counter,DataVar,FitFunVar,FemRegVar,BlRegVar,LamParVar,RadVar);
                    if(double(TIMER->Stop())>TimeLimit*1e6) break;
                }
            }
            for(unsigned uData=0; uData<NumDataSets; uData++){
                DataVar[uData]=CurrentDataVar[uData];
                FitFunVar[uData]=CurrentFitFunVar[uData];
                FemRegVar[uData]=CurrentFemRegVar[uData];
                BlRegVar[uData]=CurrentBlRegVar[uData];
                LamParVar[uData]=CurrentLamParVar[uData];
                RadVar[uData]=CurrentRadVar[uData];

                if(!Updated){
                    unsigned WhichData = uData;
                    while(!Updated){
                        if(WhichData>=NumDataSets){
                            CONTINUE=false;
                            break;
                        }
                        CurrentRadVar[WhichData]++;
                        if(CurrentRadVar[WhichData]>=NumRadVar[WhichData]){
                            CurrentRadVar[WhichData]=0;
                            CurrentLamParVar[WhichData]++;
                        }
                        if(CurrentLamParVar[WhichData]>=NumLamPar[WhichData][CurrentDataVar[WhichData]]){
                            CurrentLamParVar[WhichData]=0;
                            CurrentBlRegVar[WhichData]++;
                        }
                        if(CurrentBlRegVar[WhichData]>=NumBlReg[WhichData][CurrentDataVar[WhichData]]){
                            CurrentBlRegVar[WhichData]=0;
                            CurrentFemRegVar[WhichData]++;
                        }
                        if(CurrentFemRegVar[WhichData]>=NumFemtoReg[WhichData][CurrentDataVar[WhichData]]){
                            CurrentFemRegVar[WhichData]=0;
                            CurrentFitFunVar[WhichData]++;
                        }
                        if(CurrentFitFunVar[WhichData]>=NumFitFunctions[WhichData]){
                            CurrentFitFunVar[WhichData]=0;
                            CurrentDataVar[WhichData]++;
                        }
                        if(CurrentDataVar[WhichData]>=NumDataVar[WhichData]){
                            CurrentDataVar[WhichData]=0;
                            WhichData++;
                        }
                        else{
                            Updated = true;
                        }
                    }
                }
                if(!CONTINUE) break;
            }
            Counter++;
        }

        delete [] CurrentDataVar;
        delete [] CurrentFitFunVar;
        delete [] CurrentFemRegVar;
        delete [] CurrentBlRegVar;
        delete [] CurrentLamParVar;
        delete [] CurrentRadVar;
    }
    else{
        for(int iSystVar=0; iSystVar<NumSystVariations; iSystVar++){
            for(unsigned uData=0; uData<NumDataSets; uData++){
                if(!iSystVar&&RunDefaultSyst){
                    DataVar[uData]=0;
                    FitFunVar[uData]=0;
                    FemRegVar[uData]=0;
                    BlRegVar[uData]=0;
                    LamParVar[uData]=0;
                    RadVar[uData]=0;
                }
                else{
                    DataVar[uData]=RanGenSyst->Integer(NumDataVar[uData]);
                    FitFunVar[uData]=RanGenSyst->Integer(NumFitFunctions[uData]);
                    FemRegVar[uData]=RanGenSyst->Integer(NumFemtoReg[uData][DataVar[uData]]);
                    BlRegVar[uData]=RanGenSyst->Integer(NumBlReg[uData][DataVar[uData]]);
                    LamParVar[uData]=RanGenSyst->Integer(NumLamPar[uData][DataVar[uData]]);
                    RadVar[uData]=RanGenSyst->Integer(NumRadVar[uData]);
                }
            }
            PerformSingleFit(-1-iSystVar,DataVar,FitFunVar,FemRegVar,BlRegVar,LamParVar,RadVar);
            //#pragma omp parallel
            for(int iBoot=0; iBoot<NumBootIter; iBoot++){
                PerformSingleFit(iSystVar,DataVar,FitFunVar,FemRegVar,BlRegVar,LamParVar,RadVar);
                if(double(TIMER->Stop())>TimeLimit*1e6) break;
            }
        }
    }

    OutputFile->cd();
    ntParDensity->Write();

    delete [] DataVar;
    delete [] FitFunVar;
    delete [] FemRegVar;
    delete [] BlRegVar;
    delete [] LamParVar;
    delete [] RadVar;

}

//if uIter==-1 => we simply compute the reference, given this systematic variations
void LamLamMcExclusion::PerformSingleFit(const long int& iIter,
                                        const unsigned* DataVar, const unsigned* FitFunVar, const unsigned* FemRegVar,
                                        const unsigned* BlRegVar, const unsigned* LamParVar, const unsigned* RadVar){

    DLM_Ck*** Ck_Dlm;
    Ck_Dlm = new DLM_Ck** [NumDataSets];

    //in the end when we fix the scat. pars we want to fit each data set individually due to CPU time
    unsigned NumFitters=2+(iIter>=0||SaveChi2Map)*NumDataSets;

    //fit the original with the scattering parameters obtained form fitter0 and 1
    //the fits for the bootstrap. 0,1 only differ by the starting point of the scattering parameters
    //0 biased towards an attractive potential
    //1 biased towards a bound state
    //at the very and the best fit between 0,1 is selected
    //the fitters 2 is used to refit the original data with the best result from 0,1
    DLM_Fitter1** Fitter = new DLM_Fitter1* [NumFitters];
    for(unsigned uFit=0; uFit<NumFitters; uFit++){
        if(uFit<2) Fitter[uFit] = new DLM_Fitter1(NumDataSets);
        else Fitter[uFit] = new DLM_Fitter1(1);
        //Fitter[uFit]->RemoveNegativeCk(true);
    }

    DLM_CkDecomposition*** CkDec = new DLM_CkDecomposition** [NumDataSets];

    TH1F** hBootstrap = new TH1F* [NumDataSets];
    TH1F** hOriginal = new TH1F* [NumDataSets];
    TH1F*** hToFit = new TH1F** [NumDataSets];

    //define all needed inputs for the Fitter, which depend on the data set
    for(unsigned uData=0; uData<NumDataSets; uData++){//uData
        hToFit[uData] = new TH1F* [NumFitters];
        hOriginal[uData] = hData[uData][DataVar[uData]];
        double FemtoMax = FemtoRegion[uData][DataVar[uData]][FemRegVar[uData]][1];
        double FemtoMin = FemtoRegion[uData][DataVar[uData]][FemRegVar[uData]][0];
        double BlMax = BlRegion[uData][DataVar[uData]][BlRegVar[uData]][1];
        double BlMin = BlRegion[uData][DataVar[uData]][BlRegVar[uData]][0];
        unsigned BinFemMax = hOriginal[uData]->GetXaxis()->FindBin(FemtoMax);
        unsigned BinFemMin = hOriginal[uData]->GetXaxis()->FindBin(FemtoMin);
        unsigned BinBlMax = hOriginal[uData]->GetXaxis()->FindBin(BlMax);
        if(hOriginal[uData]->GetXaxis()->GetBinLowEdge(BinFemMax)==FemtoMax) BinFemMax--;
        if(hOriginal[uData]->GetXaxis()->GetBinLowEdge(BinBlMax)==BlMax) BinBlMax--;
        unsigned NumFemtoBins = BinFemMax-BinFemMin+1;
        unsigned NumTotBins = BinBlMax-BinFemMin+1;

        //sample random data using the Bootstrap algorithm
        if(iIter>=0){
            double CkVal;
            double CkErr;
            double CkBoot;
            //double CkErrBoot;
            double Momentum;
            unsigned OriginalBin;
            double EstNumSE;
            double EstNumME;
            double SEboot;
            hBootstrap[uData] = new TH1F(TString::Format("hBootstrap_%u",uData),TString::Format("hBootstrap_%u",uData),
                                         NumTotBins,FemtoMin,BlMax);
            for(unsigned uBin=1; uBin<=NumTotBins; uBin++){
                Momentum = hBootstrap[uData]->GetBinCenter(uBin);
                OriginalBin = hOriginal[uData]->FindBin(Momentum);
                if(Bootstrap) CkVal = hOriginal[uData]->GetBinContent(OriginalBin);
                else CkVal = SamplingGraph[uData].Eval(Momentum);
///N.B. SamplingHisto[uData] is the same as hOriginal[uData], maybe I can optimize how I handle this, to be less confusing
                CkErr = hOriginal[uData]->GetBinError(OriginalBin);
                EstNumSE = round(pow(CkVal/CkErr,2.));
                EstNumME = round(EstNumSE/CkVal);
                SEboot = RanGenMC->Poisson(EstNumSE);
                CkBoot = SEboot/EstNumME;
                hBootstrap[uData]->SetBinContent(uBin,CkBoot);
                hBootstrap[uData]->SetBinError(uBin,CkErr); }
        }
        else{
            hBootstrap[uData] = NULL;
        }

        for(unsigned uFit=0; uFit<NumFitters; uFit++){
            if(hBootstrap[uData] && uFit<2) hToFit[uData][uFit] = hBootstrap[uData];
            else if(uFit<2) hToFit[uData][uFit] = hOriginal[uData];
            else if(!SamplingHisto[uData]) hToFit[uData][uFit] = hOriginal[uData];
            else hToFit[uData][uFit] = SamplingHisto[uData];

        }

        Ck_Dlm[uData] = new DLM_Ck* [NumFitters];
        CkDec[uData] = new DLM_CkDecomposition* [NumFitters];

        for(unsigned uFit=0; uFit<NumFitters; uFit++){
            //so we want to fit each data set individually, when we make the chi2 map => we only need 1 data set per fitter
            if( uFit>=2 && (uFit-2)!=uData ){
                Ck_Dlm[uData][uFit] = NULL;
                CkDec[uData][uFit] = NULL;
                continue;
            }
            Ck_Dlm[uData][uFit] = new DLM_Ck(1,2,NumFemtoBins,FemtoMin,FemtoMax,Lednicky_Identical_Singlet_InvScatLen);//_InvScatLen
            //! careful with implementing the variations on the radius
            Ck_Dlm[uData][uFit]->SetSourcePar(0,Radius[uData][RadVar[uData]]);
            //! also consider fitting with two different starting parameters
            if(uFit%2==0){
                Ck_Dlm[uData][uFit]->SetPotPar(0,1);
                Ck_Dlm[uData][uFit]->SetPotPar(1,4);
            }
            else{
                Ck_Dlm[uData][uFit]->SetPotPar(0,-1);
                Ck_Dlm[uData][uFit]->SetPotPar(1,1);
            }

            Ck_Dlm[uData][uFit]->Update();

            CkDec[uData][uFit] = new DLM_CkDecomposition(TString::Format("LL_%u",uData),2,*Ck_Dlm[uData][uFit],hResMatrix[uData][DataVar[uData]]);
            CkDec[uData][uFit]->AddContribution(0,LambdaParF[uData][DataVar[uData]][LamParVar[uData]],DLM_CkDecomposition::cFeedDown);
            if(LambdaParF[uData][DataVar[uData]][LamParVar[uData]]-LambdaPar0[uData][DataVar[uData]][LamParVar[uData]]>1){
                printf("\033[1;31mERROR: \033[0mWrong lambda parameters!\n");
                return;
            }
            CkDec[uData][uFit]->AddContribution(1,1.-LambdaParF[uData][DataVar[uData]][LamParVar[uData]]-LambdaPar0[uData][DataVar[uData]][LamParVar[uData]],DLM_CkDecomposition::cFake);
            CkDec[uData][uFit]->Update();

            unsigned uDataFitter = uFit<2?uData:0;

            if(FitFunction[uData][FitFunVar[uData]]>=6){
                Fitter[uFit]->SetSystem(uDataFitter,*hToFit[uData][uFit],1,*CkDec[uData][uFit],
                                        FemtoMin,FemtoMax,FemtoMax,FemtoMax);
                Fitter[uFit]->SetSeparateBL(uDataFitter,false);
                Fitter[uFit]->SetParameter(TString::Format("LL_%u",uData),DLM_Fitter1::p_a,1.0,0.7,1.3);
                if(FitFunction[uData][FitFunVar[uData]]>=7) Fitter[uFit]->SetParameter(TString::Format("LL_%u",uData),DLM_Fitter1::p_b,1e-4,-2e-3,2e-3);
                else Fitter[uFit]->FixParameter(TString::Format("LL_%u",uData),DLM_Fitter1::p_b,0);
                if(FitFunction[uData][FitFunVar[uData]]>=8) Fitter[uFit]->SetParameter(TString::Format("LL_%u",uData),DLM_Fitter1::p_c,0,-2e-4,2e-4);
                else Fitter[uFit]->FixParameter(TString::Format("LL_%u",uData),DLM_Fitter1::p_c,0);
                Fitter[uFit]->FixParameter(TString::Format("LL_%u",uData),DLM_Fitter1::p_Cl,-1);
            }
            //femto and BL region are glued together, we always make a global fit, as else it is a fucking mess
            //the BL minimum value is ignored
            else if(FitFunction[uData][FitFunVar[uData]]>=3){
                Fitter[uFit]->SetSystem(uDataFitter,*hToFit[uData][uFit],1,*CkDec[uData][uFit],
                                        FemtoMin,FemtoMax,FemtoMax,BlMax);
                Fitter[uFit]->SetSeparateBL(uDataFitter,false);
                Fitter[uFit]->SetParameter(TString::Format("LL_%u",uData),DLM_Fitter1::p_a,1.0,0.7,1.3);
                if(FitFunction[uData][FitFunVar[uData]]>=4) Fitter[uFit]->SetParameter(TString::Format("LL_%u",uData),DLM_Fitter1::p_b,1e-4,-2e-3,2e-3);
                else Fitter[uFit]->FixParameter(TString::Format("LL_%u",uData),DLM_Fitter1::p_b,0);
                if(FitFunction[uData][FitFunVar[uData]]>=5) Fitter[uFit]->SetParameter(TString::Format("LL_%u",uData),DLM_Fitter1::p_c,0,-2e-4,2e-4);
                else Fitter[uFit]->FixParameter(TString::Format("LL_%u",uData),DLM_Fitter1::p_c,0);

                Fitter[uFit]->FixParameter(TString::Format("LL_%u",uData),DLM_Fitter1::p_kc,FemtoMax+8);
            }
            else{
                if(SeparateBaseline[uData]==-1){
                    Fitter[uFit]->SetSystem(uDataFitter,*hToFit[uData][uFit],1,*CkDec[uData][uFit],
                                        FemtoMin,FemtoMax,FemtoMax,FemtoMax);
                    Fitter[uFit]->SetSeparateBL(uDataFitter,false);
                }
                else{
                    Fitter[uFit]->SetSystem(uDataFitter,*hToFit[uData][uFit],1,*CkDec[uData][uFit],
                                        FemtoMin,FemtoMax,BlMin,BlMax);
                    Fitter[uFit]->SetSeparateBL(uDataFitter,SeparateBaseline[uData]==1);
                }


                if(FitFunction[uData][FitFunVar[uData]]==-1){
                    Fitter[uFit]->SetSystem(uDataFitter,*hToFit[uData][uFit],1,*CkDec[uData][uFit],
                                        FemtoMin,FemtoMax,FemtoMax,FemtoMax);
                    Fitter[uFit]->SetSeparateBL(uDataFitter,false);
                    Fitter[uFit]->FixParameter(TString::Format("LL_%u",uData),DLM_Fitter1::p_a,1.0);
                    Fitter[uFit]->FixParameter(TString::Format("LL_%u",uData),DLM_Fitter1::p_b,0);
                    Fitter[uFit]->FixParameter(TString::Format("LL_%u",uData),DLM_Fitter1::p_c,0);
                }
                else{
                    Fitter[uFit]->SetParameter(TString::Format("LL_%u",uData),DLM_Fitter1::p_a,1.0,0.7,1.3);
                    if(FitFunction[uData][FitFunVar[uData]]>=1) Fitter[uFit]->SetParameter(TString::Format("LL_%u",uData),DLM_Fitter1::p_b,1e-4,-2e-3,2e-3);
                    else Fitter[uFit]->FixParameter(TString::Format("LL_%u",uData),DLM_Fitter1::p_b,0);
                    if(FitFunction[uData][FitFunVar[uData]]>=2) Fitter[uFit]->SetParameter(TString::Format("LL_%u",uData),DLM_Fitter1::p_c,0,-2e-4,2e-4);
                    else Fitter[uFit]->FixParameter(TString::Format("LL_%u",uData),DLM_Fitter1::p_c,0);
                }
                Fitter[uFit]->FixParameter(TString::Format("LL_%u",uData),DLM_Fitter1::p_Cl,-1);
            }

            if(uFit%2==0){
                Fitter[uFit]->SetParameter(TString::Format("LL_%u",uData),DLM_Fitter1::p_pot0,0.5,0.2,32);//2 0.005,32
                Fitter[uFit]->SetParameter(TString::Format("LL_%u",uData),DLM_Fitter1::p_pot1,4.5,0.01,64);//8 0.01 64
            }
            else{
                Fitter[uFit]->SetParameter(TString::Format("LL_%u",uData),DLM_Fitter1::p_pot0,-0.25,-0.8,-0.1);//-0.6,-0.005
                Fitter[uFit]->SetParameter(TString::Format("LL_%u",uData),DLM_Fitter1::p_pot1,0.5,0.01,64);
            }

            Fitter[uFit]->FixParameter(TString::Format("LL_%u",uData),DLM_Fitter1::p_sor0,Radius[uData][RadVar[uData]]);

            if(uData&&uFit<2){
                Fitter[uFit]->AddSamePotential(TString::Format("LL_%u",uData),"LL_0",2);
            }
        }//uFit
    }//uData

    Fitter[0]->GoBabyGo();
    Fitter[1]->GoBabyGo();
    unsigned BestFitter = (Fitter[0]->GetChi2()<Fitter[1]->GetChi2())?0:1;
    //fix the scattering parameters in the fitters that will compare the result to the baseline data
    if(iIter>=0){
        for(unsigned uFit=2; uFit<NumFitters; uFit++){
            Fitter[uFit]->FixParameter(TString::Format("LL_%u",uFit-2),DLM_Fitter1::p_pot0,Fitter[BestFitter]->GetParameter(0,DLM_Fitter1::p_pot0));
            Fitter[uFit]->FixParameter(TString::Format("LL_%u",uFit-2),DLM_Fitter1::p_pot1,Fitter[BestFitter]->GetParameter(0,DLM_Fitter1::p_pot1));
            Fitter[uFit]->GoBabyGo();
        }
    }

    //the default case => save the sampling histo
    if(iIter<=-1){
        for(unsigned uData=0; uData<NumDataSets; uData++){
            SamplingHisto[uData] = (TH1F*)hOriginal[uData]->Clone(TString::Format("SamplingHisto_%u",uData));
            SamplingHisto[uData] = new TH1F(*hOriginal[uData]);
            SamplingHisto[uData]->SetName(TString::Format("SamplingHisto_%u",uData));
            SamplingHisto[uData]->SetTitle(TString::Format("SamplingHisto_%u",uData));

            Fitter[BestFitter]->GetFitGraph(uData,SamplingGraph[uData]);
            SamplingGraph[uData].SetName(TString::Format("gLL_%u",uData));
        }
    }

    double Chi2_ToOriginal=0;
    if(iIter>=0){
        for(unsigned uFit=2; uFit<NumFitters; uFit++){
            Chi2_ToOriginal += Fitter[uFit]->GetChi2();
        }
    }

    if(iIter<=-1){
        Chi2min=Fitter[BestFitter]->GetChi2();
        NDFmin=Fitter[BestFitter]->GetNdf();
        f0Inv_min=Fitter[BestFitter]->GetParameter(TString::Format("LL_0"),DLM_Fitter1::p_pot0);
        d0_min=Fitter[BestFitter]->GetParameter(TString::Format("LL_0"),DLM_Fitter1::p_pot1);
        printf("Chi2min=%.2f; NDFmin=%u\n",Chi2min,NDFmin);

        if(SaveFitResults){
            TGraph* gFitResult = new TGraph [NumDataSets];
            for(unsigned uData=0; uData<NumDataSets; uData++){
                gFitResult[uData].SetName(TString::Format("grFitRes_%u_%ld",uData,-iIter-1));
                Fitter[BestFitter]->GetFitGraph(uData,gFitResult[uData]);
                //printf("Data set %u: a=%.3e; b=%.3e\n",uData,Fitter[BestFitter]->GetParameter(uData,DLM_Fitter1::p_a),Fitter[BestFitter]->GetParameter(uData,DLM_Fitter1::p_b));
                //printf("         %.2f; %.2f (r=%.2f)\n",Fitter[BestFitter]->GetParameter(uData,DLM_Fitter1::p_pot0),Fitter[BestFitter]->GetParameter(uData,DLM_Fitter1::p_pot1),
                 //      Fitter[BestFitter]->GetParameter(uData,DLM_Fitter1::p_sor0));

                //printf("  at %.2f ck is %.2f\n",100.,Ck_Dlm[uData][BestFitter]->Eval(100));
                //printf("  pp0 = %.2f; pp1 = %.2f\n",Ck_Dlm[uData][BestFitter]->GetPotPar(0),Ck_Dlm[uData][BestFitter]->GetPotPar(1));

                OutputFile->cd();
                gFitResult[uData].Write();
            }
            delete [] gFitResult;
        }

        double Val_f0Inv;
        double Val_d0;

        if(SaveChi2Map){
            TH2F* hExclPlot = new TH2F(TString::Format("hExclPlot_%li",-iIter-1),TString::Format("hExclPlot_%li",-iIter-1),NumBins_f0Inv,Min_f0Inv,Max_f0Inv,NumBins_d0,Min_d0,Max_d0);
            TH2F* hNegCk = new TH2F(TString::Format("hNegCk_%li",-iIter-1),TString::Format("hNegCk_%li",-iIter-1),NumBins_f0Inv,Min_f0Inv,Max_f0Inv,NumBins_d0,Min_d0,Max_d0);
            for(unsigned uBin_f0=0; uBin_f0<NumBins_f0Inv; uBin_f0++){
                Val_f0Inv = hExclPlot->GetXaxis()->GetBinCenter(uBin_f0+1);
                for(unsigned uBin_d0=0; uBin_d0<NumBins_d0; uBin_d0++){
                    Val_d0 = hExclPlot->GetYaxis()->GetBinCenter(uBin_d0+1);
                    //Fitter[BestFitter]->FixParameter(TString::Format("LL_0"),DLM_Fitter1::p_pot0,Val_f0Inv);
                    //Fitter[BestFitter]->FixParameter(TString::Format("LL_0"),DLM_Fitter1::p_pot1,Val_d0);
                    //Fitter[BestFitter]->GoBabyGo();
                    //hExclPlot->SetBinContent(uBin_f0+1,uBin_d0+1,(Fitter[BestFitter]->GetChi2()-Chi2min)/double(NDFmin));
                    double GlobChi2=0;
                    hNegCk->SetBinContent(uBin_f0+1,uBin_d0+1,0);
                    for(unsigned uFit=2; uFit<NumFitters; uFit++){
                        Fitter[uFit]->FixParameter(TString::Format("LL_%u",uFit-2),DLM_Fitter1::p_pot0,Val_f0Inv);
                        Fitter[uFit]->FixParameter(TString::Format("LL_%u",uFit-2),DLM_Fitter1::p_pot1,Val_d0);
                        Fitter[uFit]->GoBabyGo();
                        GlobChi2+=Fitter[uFit]->GetChi2();
                        if(hNegCk->GetBinContent(uBin_f0+1,uBin_d0+1)==0 &&
                           Fitter[uFit]->CheckNegativeCk()) {hNegCk->SetBinContent(uBin_f0+1,uBin_d0+1,1);}
                    }
                    hExclPlot->SetBinContent(uBin_f0+1,uBin_d0+1,(GlobChi2-Chi2min)/double(NDFmin));
                    //hMinCk->SetBinContent(uBin_f0+1,uBin_d0+1,MinCk);
                }
            }

            //save the output in the default output file
            if(!Chi2MapFileName){
                OutputFile->cd();
                hExclPlot->Write();
                hNegCk->Write();
            }
            //save the output in a file common to all jobs => make sure the file is not opened when you modify it
            else{
                bool FileIsOpen=false;
                double WaitedInSec=0;
                do{
                    if(FileIsOpen) {usleep(1000); WaitedInSec+=1e-3;}//sleep for a millisecond
                    if(WaitedInSec>60){
                        printf("\033[1;31mERROR: \033[0mCannot access the file %s\n", Chi2MapFileName);
                        printf("\033[1;31m ABORTING THE EXECUTABLE!\033[0m\n");
                        abort();
                    }
                    int fd = open(Chi2MapFileName, O_RDONLY);
                    if (fd < 0) {
                        //perror("open");
                        //file does not exist
                        FileIsOpen = false;
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
                TFile* Chi2MapFile = new TFile(Chi2MapFileName,"UPDATE");

                TH2F* hChi2NdfMap;
                hChi2NdfMap = (TH2F*)Chi2MapFile->Get("hChi2NdfMap");
                TH2F* hNegativeLedni;
                hNegativeLedni = (TH2F*)Chi2MapFile->Get("hNegativeLedni");
                bool Create_Chi2Map = !hChi2NdfMap;
                //if the file does not exist, than
                if(Create_Chi2Map){
                    hChi2NdfMap = new TH2F(TString::Format("hChi2NdfMap"),TString::Format("hChi2NdfMap"),NumBins_f0Inv,Min_f0Inv,Max_f0Inv,NumBins_d0,Min_d0,Max_d0);
                    hNegativeLedni = new TH2F(TString::Format("hNegativeLedni"),TString::Format("hNegativeLedni"),NumBins_f0Inv,Min_f0Inv,Max_f0Inv,NumBins_d0,Min_d0,Max_d0);
                    for(unsigned uBin_f0=0; uBin_f0<NumBins_f0Inv; uBin_f0++){
                        for(unsigned uBin_d0=0; uBin_d0<NumBins_d0; uBin_d0++){
                            hChi2NdfMap->SetBinContent(uBin_f0+1,uBin_d0+1,hExclPlot->GetBinContent(uBin_f0+1,uBin_d0+1));
                            hNegativeLedni->SetBinContent(uBin_f0+1,uBin_d0+1,hNegCk->GetBinContent(uBin_f0+1,uBin_d0+1));
                        }
                    }
                }
                else{
                    double Val_Chi2ndf_Old;
                    double Val_Chi2ndf_New;
                    for(unsigned uBin_f0=0; uBin_f0<NumBins_f0Inv; uBin_f0++){
                        for(unsigned uBin_d0=0; uBin_d0<NumBins_d0; uBin_d0++){
                            Val_Chi2ndf_Old = hExclPlot->GetBinContent(uBin_f0+1,uBin_d0+1);
                            Val_Chi2ndf_New = hChi2NdfMap->GetBinContent(uBin_f0+1,uBin_d0+1);
                            //set the lower value
                            hChi2NdfMap->SetBinContent(uBin_f0+1,uBin_d0+1,Val_Chi2ndf_New<Val_Chi2ndf_Old?Val_Chi2ndf_New:Val_Chi2ndf_Old);
                            if(hNegativeLedni->GetBinContent(uBin_f0+1,uBin_d0+1)==0 && hNegCk->GetBinContent(uBin_f0+1,uBin_d0+1)==0)
                                hNegativeLedni->SetBinContent(uBin_f0+1,uBin_d0+1,1);
                        }
                    }
                }

                Chi2MapFile->cd();
                hChi2NdfMap->Write(0,2,0);
                hNegativeLedni->Write(0,2,0);
                if(Create_Chi2Map) {delete hChi2NdfMap; delete hNegativeLedni;}
                Chi2MapFile->Close();
                delete Chi2MapFile;
            }

            delete hExclPlot;
            delete hNegCk;
        }
    }
    else{
            //UITERfirst:UITERlast:f0Inv:d0:Chi2fit:NDFfit:dChi2
            Float_t CONTAINTER[9];
            CONTAINTER[0] = iIter/100000;
            CONTAINTER[1] = iIter%100000;
            CONTAINTER[2] = Fitter[BestFitter]->GetParameter("LL_0",DLM_Fitter1::p_pot0);
            CONTAINTER[3] = Fitter[BestFitter]->GetParameter("LL_0",DLM_Fitter1::p_pot1);
            CONTAINTER[4] = Fitter[BestFitter]->GetChi2();
            CONTAINTER[5] = Fitter[BestFitter]->GetNdf();
            CONTAINTER[6] = Chi2_ToOriginal;
            CONTAINTER[7] = Chi2_ToOriginal-Chi2min;
            CONTAINTER[8] = (Chi2_ToOriginal-Chi2min)/double(NDFmin);
            if(CONTAINTER[7]<0){
                printf("\033[1;33mWARNING:\033[0m Chi2_ToOriginal-Chi2min = %.2e\n",Chi2_ToOriginal-Chi2min);
                printf(" This happened at 1/f0=%.2f and d0=%.2f (best are 1/f0=%.2f and d0=%.2f)\n",
                       Fitter[BestFitter]->GetParameter("LL_0",DLM_Fitter1::p_pot0),Fitter[BestFitter]->GetParameter("LL_0",DLM_Fitter1::p_pot1),
                       f0Inv_min,d0_min);
            }
            ntParDensity->Fill(CONTAINTER);
        }

        for(unsigned uData=0; uData<NumDataSets; uData++){
            for(unsigned uFit=0; uFit<NumFitters; uFit++){
                if(Ck_Dlm[uData][uFit]){delete Ck_Dlm[uData][uFit]; Ck_Dlm[uData][uFit]=NULL;}
                if(CkDec[uData][uFit]){delete CkDec[uData][uFit]; CkDec[uData][uFit]=NULL;}
            }
            delete Ck_Dlm[uData]; Ck_Dlm[uData]=NULL;
            delete CkDec[uData]; CkDec[uData]=NULL;
            //delete hMCdata[uData]; hMCdata[uData]=NULL;
        }

        delete [] Ck_Dlm;
        delete [] CkDec;
        for(unsigned uData=0; uData<NumDataSets; uData++){
            if(hBootstrap[uData]){
                delete hBootstrap[uData];
                hBootstrap[uData] = NULL;
            }
            delete[]hToFit[uData]; hToFit[uData]=NULL;
        }
        delete [] hBootstrap;
        delete [] hOriginal;
        delete [] hToFit;

        for(unsigned uFit=0; uFit<NumFitters; uFit++){
            delete Fitter[uFit];
        }
        delete [] Fitter;

}

TNtuple* LamLamMcExclusion::GetOutputNtuple(){

    return ntParDensity;

}

