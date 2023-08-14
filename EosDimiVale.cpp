#include <iostream>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>

#include "CATS.h"
#include "CommonAnaFunctions.h"
#include "DLM_Potentials.h"
#include "DLM_Random.h"
#include "DLM_CkDecomposition.h"
#include "DLM_Ck.h"
#include "DLM_Source.h"
#include "DLM_Potentials.h"
#include "DLM_WfModel.h"
#include "DLM_CkModels.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TNtuple.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TNtuple.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TFractionFitter.h"
#include "TGenPhaseSpace.h"
#include "TString.h"
#include "DLM_Fitters.h"
#include "DLM_Histo.h"
#include "TROOT.h"
#include "TMath.h"
#include "TSystem.h"
#include "TLatex.h"
#include "DLM_Integration.h"
#include "CATSconstants.h"
#include "DLM_SubPads.h"
#include "CommonAnaFunctions.h"
#include "TLine.h"
#include "TColor.h"

// make script that evaluates the cross-section using the CATS object already setup and returns the chi-square
const unsigned NumMomBins = 40;
const double kMin = 0;
const double kMax = 160;
const double kChi2 = 135;
const double Eps = 1e-9;

void CrossSectionFit_pL(CATS& Kitty, double& chi2, int& ndp)
{
  // 1) Reading the scattering data
  // 0 = take all
  // 1 = take Haidenbauer (we take the second set from alex, where bins are chosen to have same yield)
  // 2 = take all published
  int DataFlag = 1;

  TGraphErrors *g_CS_Data = new TGraphErrors();
  g_CS_Data->SetName("pLambda_ScatteringData");
  g_CS_Data->SetLineColor(kCyan + 1);
  g_CS_Data->SetLineWidth(2);
  g_CS_Data->SetMarkerColor(kCyan + 1);
  g_CS_Data->SetMarkerSize(1.25);
  g_CS_Data->SetMarkerStyle(20);

  g_CS_Data->GetXaxis()->SetLabelSize(0.055);
  g_CS_Data->GetXaxis()->SetTitleSize(0.06);
  g_CS_Data->GetXaxis()->SetLabelOffset(0.01);
  g_CS_Data->GetXaxis()->SetTitleOffset(1.0);
  g_CS_Data->GetXaxis()->SetLabelFont(42);
  g_CS_Data->GetYaxis()->SetLabelSize(0.055);
  g_CS_Data->GetYaxis()->SetTitleSize(0.06);

  g_CS_Data->SetTitle("; #it{k*} (MeV/#it{c}); #it{#sigma} (mb)");
  // g_CS_Data->GetYaxis()->SetRangeUser(0, 450);
  g_CS_Data->GetXaxis()->SetNdivisions(505);
  g_CS_Data->GetYaxis()->SetLabelOffset(0.01);
  g_CS_Data->GetYaxis()->SetTitleOffset(0.8);

  int NumPts = 0;

  // PhysRev.173.1452 = Alexander
  // PhysRev.175.1735 = B. Sechi-Zorn
  // Nucl. Phys. B 27 = F. Eisele
  // Nucl.Phys.B 125 = J.M. Hauptman
  // Phys. Rev. 159, 853 (R.C. Herndon, Y.C. Tang) --> not bublished

  if (DataFlag == 0 || DataFlag == 1 || DataFlag == 2)
  {
    // PhysRev.175.1735
    g_CS_Data->SetPoint(NumPts, pLab_pCm(135.1, Mass_L, Mass_p), 209);
    g_CS_Data->SetPointError(NumPts++, 0, 58);

    // PhysRev.173.1452
    g_CS_Data->SetPoint(NumPts, pLab_pCm(145, Mass_L, Mass_p), 180);
    g_CS_Data->SetPointError(NumPts++, 0, 22);

    // PhysRev.175.1735
    g_CS_Data->SetPoint(NumPts, pLab_pCm(165.1, Mass_L, Mass_p), 177);
    g_CS_Data->SetPointError(NumPts++, 0, 38);

    // PhysRev.173.1452
    g_CS_Data->SetPoint(NumPts, pLab_pCm(185, Mass_L, Mass_p), 130);
    g_CS_Data->SetPointError(NumPts++, 0, 17);

    // PhysRev.175.1735
    g_CS_Data->SetPoint(NumPts, pLab_pCm(194, Mass_L, Mass_p), 153);
    g_CS_Data->SetPointError(NumPts++, 0, 27);

    // PhysRev.173.1452
    g_CS_Data->SetPoint(NumPts, pLab_pCm(210, Mass_L, Mass_p), 118);
    g_CS_Data->SetPointError(NumPts++, 0, 16);

    // PhysRev.175.1735
    g_CS_Data->SetPoint(NumPts, pLab_pCm(226, Mass_L, Mass_p), 111);
    g_CS_Data->SetPointError(NumPts++, 0, 18);

    // PhysRev.173.1452
    g_CS_Data->SetPoint(NumPts, pLab_pCm(230, Mass_L, Mass_p), 101);
    g_CS_Data->SetPointError(NumPts++, 0, 12);

    // PhysRev.173.1452
    g_CS_Data->SetPoint(NumPts, pLab_pCm(250, Mass_L, Mass_p), 83);
    g_CS_Data->SetPointError(NumPts++, 0, 9);

    // PhysRev.175.1735
    g_CS_Data->SetPoint(NumPts, pLab_pCm(252, Mass_L, Mass_p), 87);
    g_CS_Data->SetPointError(NumPts++, 0, 13);

    // PhysRev.173.1452
    g_CS_Data->SetPoint(NumPts, pLab_pCm(290, Mass_L, Mass_p), 57);
    g_CS_Data->SetPointError(NumPts++, 0, 9);

    // PhysRev.175.1735
    g_CS_Data->SetPoint(NumPts, pLab_pCm(293, Mass_L, Mass_p), 46);
    g_CS_Data->SetPointError(NumPts++, 0, 11);
  }

  //printf("pLab<=300 MeV (k*<=%.0f): #%u\n", pLab_pCm(300, Mass_L, Mass_p), NumPts);
  //printf("pLab<=345 MeV (k*<=%.0f): #%u\n", pLab_pCm(345, Mass_L, Mass_p), NumPts);

  //printf("pLab<=360 MeV (k*<=%.0f): #%u\n", pLab_pCm(360, Mass_L, Mass_p), NumPts);
  //printf("pLab<=380 MeV (k*<=%.0f): #%u\n", pLab_pCm(380, Mass_L, Mass_p), NumPts);
  //printf("pLab<=400 MeV (k*<=%.0f): #%u\n", pLab_pCm(400, Mass_L, Mass_p), NumPts);
  //printf("pLab<=420 MeV (k*<=%.0f): #%u\n", pLab_pCm(420, Mass_L, Mass_p), NumPts);

  //printf("pLab<=450 MeV (k*<=%.0f): #%u\n", pLab_pCm(450, Mass_L, Mass_p), NumPts);

  //TFile fOutput(TString::Format("/Users/sartozza/cernbox/EoSFemto/EoSPaperPheno/pL_CS_Data.root"), "recreate");
  //g_CS_Data->Write();

  TH1F *hPsSin_UsmFemto = new TH1F("hPsSin_UsmFemto", "hPsSin_UsmFemto", NumMomBins, kMin, kMax);
  TH1F *hPsTri_UsmFemto = new TH1F("hPsTri_UsmFemto", "hPsTri_UsmFemto", NumMomBins, kMin, kMax);
  TGraph *gCs_UsmFemto = new TGraph();
  gCs_UsmFemto->SetName("gCs_UsmFemto");

  double CrossSection;
  double PhaseShift;
  double kstar;

  for (unsigned uMom = 0; uMom < NumMomBins; uMom++)
  {

    kstar = Kitty.GetMomentum(uMom);
    hPsSin_UsmFemto->SetBinContent(uMom + 1, Kitty.GetPhaseShift(uMom, 0, 0));
    hPsTri_UsmFemto->SetBinContent(uMom + 1, Kitty.GetPhaseShift(uMom, 1, 0));

    CrossSection = 0;
    PhaseShift = Kitty.GetPhaseShift(uMom, 0, 0);
    CrossSection += 0.25 * 4. * Pi * NuToFm * NuToFm * 10. / (kstar * kstar) * (pow(sin(PhaseShift), 2));
    PhaseShift = Kitty.GetPhaseShift(uMom, 1, 0);
    CrossSection += 0.75 * 4. * Pi * NuToFm * NuToFm * 10. / (kstar * kstar) * (pow(sin(PhaseShift), 2));
    gCs_UsmFemto->SetPoint(uMom, kstar, CrossSection);
  }

  Kitty.SetNotifications(CATS::nWarning);

  //3) Fit to scatt. data
  //TH1F *hFit;
  //TF1 *fitSP;
  //double ScatLen, EffRan;

  //printf("UsmFemto:\n");
  //GetScattParameters(Kitty, ScatLen, EffRan, hFit, fitSP, 3, false, false, 0);
  //printf("   S=0: f0 = %.2f; d0 = %.2f\n", ScatLen, EffRan);
  //delete hFit;
  //delete fitSP;
  //GetScattParameters(Kitty, ScatLen, EffRan, hFit, fitSP, 3, false, false, 1);
  //printf("   S=1: f0 = %.2f; d0 = %.2f\n", ScatLen, EffRan);
  //delete hFit;
  //delete fitSP;

  //hPsSin_UsmFemto->Write();
  //hPsTri_UsmFemto->Write();
  //gCs_UsmFemto->Write();

  //4) Evaluate chi-square
  double CsData;
  double CsErr;

  double Chi2_UsmFemto = 0;

  unsigned NDF = 0;
  for (unsigned uData = 0; uData < g_CS_Data->GetN(); uData++)
  {
    g_CS_Data->GetPoint(uData, kstar, CsData);
    if (kstar > kChi2)
      break;
    CsErr = g_CS_Data->GetErrorY(uData);

    Chi2_UsmFemto += pow((CsData - gCs_UsmFemto->Eval(kstar)) / CsErr, 2.);
    NDF++;
  }

  //printf("chi2/ndf:\n");
  //printf("   UsmFemto: %.1f/%u\n", Chi2_UsmFemto, NDF);

  chi2 = Chi2_UsmFemto;
  ndp = NDF;

  delete g_CS_Data;
  delete hPsSin_UsmFemto;
  delete hPsTri_UsmFemto;
  delete gCs_UsmFemto;
}
int EOSDIMIVALE(int argc, char *argv[])
{

  // 2) TEMPORARY Setting up of the CATS objce (will be setted outside in the final script)
  DLM_CommonAnaFunctions AnalysisObject;

  // CATS set up with usmani, but fixed to femto data
  CATS Kitty_Charlotte;
  Kitty_Charlotte.SetMomBins(NumMomBins, kMin, kMax);
  AnalysisObject.SetUpCats_pL(Kitty_Charlotte, "UsmaniFit", "Gauss", 0, 0);
  Kitty_Charlotte.SetEpsilonConv(Eps);
  Kitty_Charlotte.SetEpsilonProp(Eps);
  Kitty_Charlotte.SetShortRangePotential(1,0,1,2279-13.5);
  Kitty_Charlotte.SetShortRangePotential(1,0,2,0.3394);
  Kitty_Charlotte.SetShortRangePotential(1,0,3,0.2614);
  Kitty_Charlotte.KillTheCat();

  double chi2;
  int ndp;
  CrossSectionFit_pL(Kitty_Charlotte,chi2,ndp);

  return 0;
}
