#ifndef CECAPAPER_H
#define CECAPAPER_H

#include "TString.h"

class CATS;
class DLM_Ck;
class DLM_CkDecomposition;
class TH1F;
class TH2F;
class TF1;
class DLM_CecaSource_v0;
class DLM_MultiFit;
class DLM_CommonAnaFunctions;
class TString;

int CECA_PAPER(int narg, char** ARGS);

class CecaAnalysis1{
public:
    CecaAnalysis1(TString AnalysisVersion, TString SourceVersion, TString DLM_AnalysisFolder);
    ~CecaAnalysis1();

    //these functions can be executed only one time (for set up)

    //these parameters are passed to the common functions

    void SetUp_pp(const TString &POT, const int &PotVar);
    void SetUp_pL(const TString &POT, const int &PotVar);
    void SetUp_pS0(const TString &POT, const int &PotVar);
    void SetUp_pXim(const TString &POT, const int &PotVar, const double src_scale=1.06);
    void SetUp_pXi0(const TString &POT, const int &PotVar, const double src_scale=1.06);
    //the Xim1530 is automatically set up, where the source is Gaussian and set to the value of pXim source
    //evaluated at an mT of 1.8
    //void SetUp_pXim1530(const TString &POT, const int &PotVar);
    void SetUp_Decomposition(const int &Variation_p, const int &Variation_L);

    //void GaussSource();
    //void CecaSource(bool Global=true);

    TF1* GetFit(TString system, unsigned uMt);
    //double GetChi2();
    //unsigned GetNDF();
    //unsigned GetNDP();

    void DumpCurrentCk(TString OutFileName);
    //Ceca or Gauss (Global vs non-global)
    void SetUp_Fits(TString fittype = "Ceca", bool AutoSrcPars=true, bool AutoBlPars=true);
    void GoBabyGo(bool print_info=false);

    DLM_CecaSource_v0* GetSrc_pp(){return Src_pp?Src_pp:NULL;}
    DLM_CecaSource_v0* GetSrc_pL(){return Src_pL?Src_pL:NULL;}

//private:

//IGNORE mT dependence (almost flat anyways)
  enum EnumLamPar { pp_gen, pp_flt, pp_pL, pp_fake,
                    pL_gen, pL_flt, pL_pS0, pL_pXi0, pL_pXim, pL_fake,
                    pS0_gen, pS0_flt, pXi_gen, pXi_pXi1530, pXi_flt };

  //Model_pS0; 0 = chiral, 1 = flat
  //int Model_pL; the standard weird notation if we deal with Usmani => negative values
  //double FemtoRegMax_pp;//evaluate C(k) up to this value
  //double FitMax_pp;//fit up to this value
  //double CkConv_pp;//the extrapolation thingy
  enum AnaSettings {  Model_pS0, Model_pL, FemtoRegMax_pp, FitMax_pp, CkConv_pp,
                      FemtoRegMax_pL, FitMax_pL, CkConv_pL};
  double* lam;
  double* settings;

  const unsigned NumMtBins_pp;
  const unsigned NumMtBins_pL;

  const TString catsfiles_fld;
  const TString ana_ver;

  //the experimental pp/pL correlations for each [mT] bin
  TH1F** hCkExp_pp;
  TH1F** hCkExp_pL;

  //the fit function for each mT bins of pp and pL
  TF1** fCk_pp;
  TF1** fCk_pL;

  const double pp_cats_max;
  const double pp_ck_max;

  const double pL_cats_max;
  const double pL_ck_max;

  //0 - Ceca (global fit)
  //1 - Ceca (individual fits)
  //
  //11 - Gauss (individual fits)
  int FittingScenario;

  //we can have only one object for all mT bins,
  //as the mT itself is a parameter of the source that we can change

  CATS* Kitty_pp;
  CATS* Kitty_pL;
  CATS* Kitty_pS0;
  CATS* Kitty_pXi0;
  CATS* Kitty_pXim;
  CATS* Kitty_pXim1530;

  DLM_Ck* Ck_pp;
  DLM_Ck* Ck_pL;
  DLM_Ck* Ck_pS0;
  DLM_Ck* Ck_pXi0;
  DLM_Ck* Ck_pXim;
  DLM_Ck* Ck_pXim1530;

  DLM_Ck** CkMt_pp;
  DLM_Ck** CkMt_pL;
  DLM_Ck** CkMt_pS0;
  DLM_Ck** CkMt_pXi0;
  DLM_Ck** CkMt_pXim;
  //DLM_Ck** CkMt_pXim1530;



  DLM_CkDecomposition* CkDec_pp;
  DLM_CkDecomposition* CkDec_pL;
  DLM_CkDecomposition* CkDec_pS0;
  DLM_CkDecomposition* CkDec_pXi0;
  DLM_CkDecomposition* CkDec_pXim;
  DLM_CkDecomposition* CkDec_pXim1530;

  TH2F* hReso_pp;
  TH1F* hPs_pp;
  TH2F* hReso_pL;
  TH1F* hPs_pL;

  TH2F* hFeed_pp_pL;
  TH2F* hFeed_pL_pS0;
  TH2F* hFeed_pL_pXim;
  TH2F* hFeed_pL_pXi0;
  TH2F* hFeed_pXi_pXim1530;

  DLM_CecaSource_v0* Src_pp;
  DLM_CecaSource_v0* Src_pL;

  DLM_MultiFit* GrandeFitter;
  DLM_CommonAnaFunctions* AnalysisObject;

  void LoadData();
  void UpdateDecomp();

  //the parameters will be:
  //[0] = mT (will be fixed)
  //[1] = d
  //[2] = ht
  //[3] = hz or tau
  //[4] = scale
  //BL = DLM_Baseline, which goes in theory to pol4, so we have 5 pars
  //in practice:  [5] is the normalization
  //              [6] is fixed to zero (to ensure zero der at k*=0)
  //              [7] is the position (in k*) of the extrama
  //              [8] is p3 of the polynomial
  //              [9] is fixed to zero (to ensure pol3)
  //after that we allow to have some extra pars, in case we need
  //to include something else, e.g. interaction etc.
  const unsigned BaseFitPars;
  const unsigned ExtraFitPars;
  const unsigned TotFitPars;
  double FitFun_pp(double* xval, double* par);
  double FitFun_pL(double* xval, double* par);
  //in these arrays we keep track of the last parameters that the fitter wanted,
  //and compare them to the new ones. we only update our theory curves in case
  //d, ht, hz, tau or ExtraPars for the interraction changed.
  //i.e. we do not update for mT or BL changes
  double* OldPar;
  double* NewPar;

  double* mT_BinCenter_pp;
  double* mT_BinCenter_pL;

  unsigned FunctionCallCounter;
  //Fit results
  //... some TGraphs

};


int Ceca_pp_or_pL(const TString FileBase, const TString InputFolder, const TString OutputFolder, const TString LogFolder,
                  const int ParID, const int JobID, const int NumCPU=1);
int dadd_f(int argc, char *argv[]);
int dlmhst_root(int argc, char *argv[]);
int dlmhst_ceca_fit_wc0(int argc, char *argv[]);
int dlmhst_ceca_kdpfit_wc0(int argc, char *argv[]);

#endif
