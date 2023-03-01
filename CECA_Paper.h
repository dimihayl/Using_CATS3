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
    CecaAnalysis1(TString SourceVersion, TString DLM_AnalysisFolder);
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


//private:

  const unsigned NumMtBins_pp;
  const unsigned NumMtBins_pL;

  const TString cernbox_fld;

  //the experimental pp/pL correlations for each [mT] bin
  TH1F** hCkExp_pp;
  TH1F** hCkExp_pL;

  //the fit function for each mT bins of pp and pL
  TF1** fCk_pp;
  TF1** fCk_pL;

  //the pp genuine/flat/pL->pp/misid.
  //IGNORE mT dependence (almost flat anyways)
  double lam_pp_gen;
  double lam_pp_flt;
  double lam_pp_pL;
  double lam_pp_fake;

  double lam_pL_gen;
  double lam_pL_flt;
  double lam_pL_pS0;
  double lam_pL_pXi0;
  double lam_pL_pXim;
  double lam_pL_fake;

  //we model charged Xi1530 with Coulomb only
  double lam_pXi_gen;
  double lam_pXi_pXi1530;
  double lam_pXi_flt;

  double lam_pS0_gen;
  double lam_pS0_flt;

  //0 = chiral, 1 = flat
  int Model_pS0;

  //the standard weird notation
  //if we deal with Usmani => negative values
  int Model_pL;

  //evaluate C(k) up to this value
  double FemtoRegMax_pp;
  //fit up to this value
  double FitMax_pp;
  //the extrapolation thingy
  //double CutOff_pp;
  double CkConv_pp;
  const double pp_cats_max;
  const double pp_ck_max;

  //evaluate C(k) up to this value
  double FemtoRegMax_pL;
  //fit up to this value
  double FitMax_pL;
  //the extrapolation thingy
  double CkConv_pL;
  const double pL_cats_max;
  const double pL_ck_max;

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
