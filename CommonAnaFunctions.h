#ifndef COMMONANAFUNCTIONS_H
#define COMMONANAFUNCTIONS_H

class TString;
class TH1F;
class TH2F;
class CATS;
class DLM_CleverLevy;
class DLM_CleverMcLevyReso;
//class MS_GaussExp_mT_Simple;

const double Mass_pi0 = 134.9766;
const double Mass_pic = 139.57018;
const double Mass_p = 938.272;
const double Mass_L = 1115.683;
const double Mass_Xim = 1321.7;


class DLM_CommonAnaFunctions{

public:

    DLM_CommonAnaFunctions();
    ~DLM_CommonAnaFunctions();

    //! ALWAY CALL THESE FUNCTIONS AFTER YOU HAVE DEFINED THE MOMENTUM BINS!
    //SOURCE:
    //"Gauss"
    //"Cauchy"
    //"Levy_Nolan"
    //"Levy_Single"
    //"Levy_Diff"
    //"CleverLevy_Nolan"
    //"CleverLevy_Single"
    //"CleverLevy_Diff"
    //"GaussExpTotSimple_2body" (the first version)
    //"McLevyNolan_Reso" (the Monte-Carlo version without mT scaling)
    //"EPOS"
    //"EPOSrescaled" -> starts with a basis rescaling of 1.5
    //"Levy_mT_Reso" (the MC version created for pLambda analysis)
    //POT:
    //  "AV18"
    void SetUpCats_pp(CATS& Kitty, const TString& POT, const TString& SOURCE);
    //POT:
    //  "LO"
    //  "NLO"
    //  "NLO_Coupled_S"
    //  "Usmani"
    void SetUpCats_pL(CATS& Kitty, const TString& POT, const TString& SOURCE);
    //POT:
    //  "pXim_Lattice" (the first version)
    //  "pXim_HALQCD1" (the second version, THE ONE TO USE)
    void SetUpCats_pXim(CATS& Kitty, const TString& POT, const TString& SOURCE);

    void SetUpBinning_pp(const TString& DataSample, unsigned& NumMomBins, double*& MomBins, double*& FitRegion);
    void SetUpBinning_pL(const TString& DataSample, unsigned& NumMomBins, double*& MomBins, double*& FitRegion);

    //DataSamples: SystemEnergy_Trigger_Version
    //the version is there to mark the different versions based on our own analysis, it can be some short description
    //Versions:
    //  Run2paper: as used for all of the first Run2 papers (LL, pXim etc)
    //DataSamples:
    //  pp13TeV_MB_Run2paper
    //  pp13TeV_HM_Run2paper
    //  pPb5TeV_Run2paper
    //The Variation flag is there for the systematics, refer to the functions themselves for more information
    void GetPurities_p(const TString& DataSample, const int& Variation, double* Purities);
    void GetPurities_L(const TString& DataSample, const int& Variation, double* Purities);
    void GetPurities_Xim(const TString& DataSample, const int& Variation, double* Purities);
    void GetFractions_p(const TString& DataSample, const int& Variation, double* Fractions);
    void GetFractions_L(const TString& DataSample, const int& Variation, double* Fractions);
    void GetFractions_Xim(const TString& DataSample, const int& Variation, double* Fractions);
    //primary, pL->pp, XX->pp, pp fake
    void SetUpLambdaPars_pp(const TString& DataSample, const int& Variation_p, double* lambda_pars);
    //primary, pS0->pL, pXim->pL, XX->pL, pp fake
    void SetUpLambdaPars_pL(const TString& DataSample, const int& Variation_p, const int& Variation_L, double* lambda_pars);
    void SetUpLambdaPars_pXim(const TString& DataSample, const int& Variation_p, const int& Variation_Xim, double* lambda_pars);

    TH2F* GetResolutionMatrix(const TString& DataSample,const TString&& System);
    TH2F* GetResidualMatrix(const TString&& FinalSystem, const TString& InitialSystem);
    TH1F* GetAliceExpCorrFun(const TString& DataSample,const TString& System,const int& iReb);

private:
    void Clean_CommonAnaFunctions();
    //MS_GaussExp_mT_Simple* Simple_Reso;
    DLM_CleverLevy* CleverLevy;
    DLM_CleverMcLevyReso* CleverMcLevyReso;
    const unsigned NumCleverLevyObjects;
};

/*
class DLM_Analyzer{

public:

    DLM_Analyzer();
    ~DLM_Analyzer();


private:
    TH1F* hData;
    DLM_Fitter1* fitter;
    TString System;

};
*/


#endif


