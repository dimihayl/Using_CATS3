#ifndef CECAPAPER_H
#define CECAPAPER_H

#include "TString.h"

class CATS;
class DLM_Ck;
class DLM_CkDecomposition;
class TH1F;
class DLM_CecaSource_v0;
class DLM_MultiFit;

int CECA_PAPER(int narg, char** ARGS);

class CecaAnalysis1{
public:
    CecaAnalysis1(TString SourceVersion, TString SourceFolder);
    ~CecaAnalysis1();

private:
  //the experimental pp/pL correlations for each [mT] bin
  TH1F** hCkExp_pp;
  //the pp genuine/flat/pL->pp/misid in [mT]
  double* lam_pp_gen;
  double* lam_pp_flt;
  double* lam_pp_pL;
  double* lam_pp_fake;

  double* lam_pL_gen;
  double* lam_pL_flt;
  double* lam_pL_pS0;
  double* lam_pL_pXi0;
  double* lam_pL_pXim;
  double* lam_pL_fake;

  //we model charged Xi1530 with Coulomb only
  double* lam_Xi_gen;
  double* lam_Xi_Xi1530;
  double* lam_Xi_flt;

  double* lam_pS0_gen;
  double* lam_pS0_flt;

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
  double FitCutOff_pp;

  //evaluate C(k) up to this value
  double FemtoRegMax_pL;
  //fit up to this value
  double FitMax_pL;
  //the extrapolation thingy
  double FitCutOff_pL;

  //[mT]
  CATS* Kitty_pp;
  CATS* Kitty_pL;
  CATS* Kitty_pS0;
  CATS* Kitty_pXi0;
  CATS* Kitty_pXim;

  DLM_Ck** Ck_pp;
  DLM_Ck** Ck_pL;
  DLM_Ck** Ck_pS0;
  DLM_Ck** Ck_pXi0;
  DLM_Ck** Ck_pXim;

  DLM_CkDecomposition** CkDec_pp;
  DLM_CkDecomposition** CkDec_pL;
  DLM_CkDecomposition** CkDec_pS0;
  DLM_CkDecomposition** CkDec_pXi0;
  DLM_CkDecomposition** CkDec_pXim;

  TH1F* hReso_pp;
  TH1F* hPs_pp;
  TH1F* hReso_pL;
  TH1F* hPs_pL;

  DLM_CecaSource_v0* Src_pp;
  DLM_CecaSource_v0* Src_pL;

  DLM_MultiFit* GrandeFitter;

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
