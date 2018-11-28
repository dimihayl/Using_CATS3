
class TH1F;
class TH2F;
class TGraph;
class TNtuple;
class TRandom3;
class DLM_Fitter1;
class TFile;
class TString;
class DLM_Timer;

class LamLamMcExclusion{

public:
    LamLamMcExclusion();
    ~LamLamMcExclusion();

    void SetOutputFile(const char* outfilename);
    void SetSaveChi2Map(const bool& savechi2map);
    void SetSaveFitResults(const bool& savefitresults);
    void SetUpChi2Map(const unsigned& nf0inv, const double& f0invmin, const double& f0invmax,
                                     const unsigned& nd0, const double& d0min, const double& d0max);
    void SetCommonChi2MapFile(const char* chi2filename);
    void SetRandomSeed(const unsigned& rs_syst, const unsigned& rs_mc);
    //N.B. if fixseed, the random seed is automatically assigned based on the jobid, i.e. SetRandomSeed is NOT used any more!!!
    //void SetParallelJobs(const unsigned& numjobs, const unsigned& jobid, const bool& fixseed);
    void SetNumDataSets(const unsigned& NumDS);
    void SetNumDataVar(const unsigned& uData, const unsigned& val);
    void SetNumFitFunVar(const unsigned& uData, const unsigned& val);
    void SetFitFunction(const unsigned& uData, const unsigned& uVar, const unsigned& val);
    void SetSeparateBaseline(const unsigned& uData, const int& val);
    void SetNumFemtoReg(const unsigned& uData, const unsigned& uDataVar, const unsigned& val);
    void SetFemtoRegion(const unsigned& uData, const unsigned& uDataVar, const unsigned& uVar, const double& MIN, const double& MAX);
    void SetNumBlReg(const unsigned& uData, const unsigned& uDataVar, const unsigned& val);
    void SetBlRegion(const unsigned& uData, const unsigned& uDataVar, const unsigned& uVar, const double& MIN, const double& MAX);
    void SetData(const unsigned& uData, const unsigned& uVar, TH1F* hdata);
    void SetResolution(const unsigned& uData, const unsigned& uVar, TH2F* hresolution);
    void SetNumLamPar(const unsigned& uData, const unsigned& uDataVar, const unsigned& val);
    void SetLambdaPar(const unsigned& uData, const unsigned& uDataVar, const unsigned& uLamVar, const double& lam0, const double& lamF);
    void SetNumRadVar(const unsigned& uData, const unsigned& val);
    void SetRadius(const unsigned& uData, const unsigned& uRadVar, const double& val);
    void SetBootstrap(const bool& val);
    void SetNumSystVariations(const unsigned& nrandom);
    void RunSpecificSyst(const bool& yesno, const unsigned& sfirst=0, const unsigned& slast=unsigned(-1));
    void SetRunDefaultSyst(const bool& rds);
    void SetNumBootIter(const unsigned& nrandom);
    //the time limit is given in seconds
    void SetTimeLimit(const double& timelimit);

    unsigned GetTotNumIter();

    //void SetFemtoRegionMax(const unsigned& uData, const unsigned& uVar, const double& val);
    //void SetNumBlRegions(const unsigned& uData, const unsigned& val);
    //void SetBlRegionMin(const unsigned& uData, const double& val);
    //void SetBlRegionMax(const unsigned& uData, const double& val);


    //void SetBins_f0(const unsigned& val, const double& f0inv_min, const double& f0inv_max);
    //void SetBins_d0(const unsigned& val, const double& d0_min, const double& d0_max);


    void GoBabyGo();


    TNtuple* GetOutputNtuple();

//private:

    TRandom3* RanGenSyst;
    TRandom3* RanGenMC;
    unsigned NumJobs;
    //unsigned JobID;
    //bool FixSeed;

    //we can have different data sets (e.g. pp 13 TeV, pPb 5 TeV)
    unsigned NumDataSets;
    //per data set, we can have as systematics different input
    unsigned* NumDataVar;//[DataSet]

    //[DataSet] each data set can have a different fit function.
    //we also allow to perform systematic variations of those functions
    unsigned* NumFitFunctions;//[DataSet]
    //-1 => a=1; b=0; Cl=-1 (no norm and bl, only femto region is fitted)
    //0 => a=free; b=0; Cl=-1 (only norm, bl region can be included if set up)
    //1 => a=free; b=free; Cl=-1 (with bl fit)
    int** FitFunction;//[DataSet][FitFunction]
    //[DataSet] if the baseline is separately fitted or not (at the moment w/o systematic variation)
    //if -1 => no bl at all
    int* SeparateBaseline;
    //[DataSet][DataVar] per data set, including systematic variations
    unsigned** NumFemtoReg;
    //[DataSet][DataVar][Variation][min/max]
    double**** FemtoRegion;
    unsigned** NumBlReg;
    double**** BlRegion;

    //[DataSet][DataVar]
    TH1F*** hData;
    TH2F*** hResMatrix;
    unsigned** NumLamPar;
    //[DataSet][DataVar][LamVar]
    double*** LambdaPar0;
    double*** LambdaParF;
    //[DataSet]: each data set is characterized by a radius, here we allow to have systematic variations
    unsigned* NumRadVar;
    //[DataSet][RadiusVariation]
    double** Radius;

    bool Bootstrap;
    //the number of syst variations. If -1 => we do all syst vars
    int NumSystVariations;
    bool SpecificSyst;
    unsigned FirstSyst;
    unsigned LastSyst;
    bool RunDefaultSyst;
    //number of random bootstrap samples per systematic variation
    int NumBootIter;
    double TimeLimit;
    DLM_Timer* TIMER;

    //TH2F* hParDensity;
    //unsigned NumBins_f0;
    //double f0Inv_Min;
    //double f0Inv_Max;
    //unsigned NumBins_d0;
    //double d0_Min;
    //double d0_Max;
    TNtuple* ntParDensity;
    TFile* OutputFile;
    char* OutputFileName;
    bool SaveChi2Map;
    bool SaveFitResults;
    char* Chi2MapFileName;
    unsigned NumBins_f0Inv;
    double Min_f0Inv;
    double Max_f0Inv;
    unsigned NumBins_d0;
    double Min_d0;
    double Max_d0;

    //this is the histo from which we should sample (and our default baseline corresponding to chi2min, GIVEN a certain systematic variation)
    TH1F** SamplingHisto;
    TGraph* SamplingGraph;
    double Chi2min;
    unsigned NDFmin;
    double f0Inv_min;
    double d0_min;

    //void FitTheReference(const unsigned* DataVar, const unsigned* FitFunVar, const unsigned* FemRegVar,
    //                                    const unsigned* BlRegVar, const unsigned* LamParVar, const unsigned* RadVar);
    void PerformSingleFit(const long int& iIter, const unsigned* DataVar, const unsigned* FitFunVar, const unsigned* FemRegVar,
                                        const unsigned* BlRegVar, const unsigned* LamParVar, const unsigned* RadVar);
    //void SingleInteraction(const unsigned& uIter, const unsigned* DataVar, const unsigned* FitFunVar, const unsigned* FemRegVar,
    //                                    const unsigned* BlRegVar, const unsigned* LamParVar, const unsigned* RadVar);

    void ClearMemDataSet(const unsigned& uData, const bool& OnlyDataVar);
    void ClearMem();

    //void ClearMem();
    //void ClearMem(const unsigned& uData);
    //void ClearMemDataVar(const unsigned& uData);
    //void ClearMemRadVar(const unsigned& uRad);
};

