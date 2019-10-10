#include "GentleDimi.h"

#include <string.h>
#include <iostream>
//#include <stdio.h>
//#include <stdlib.h>

#include "DreamKayTee.h"
#include "ReadDreamFile.h"
#include "DreamPlot.h"

#include "TROOT.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TPaveText.h"

#include "CommonAnaFunctions.h"
#include "CATS.h"
#include "DLM_CkDecomposition.h"
#include "DLM_Fitters.h"
#include "DLM_Source.h"

using namespace std;

void GetCorrelations(const char* foldername, const char* filename, const char* prefix,
                     const char* addon, double& Dev_pL) {
  //gStyle->SetOptStat(0);
  ReadDreamFile* DreamFile = new ReadDreamFile(6, 6);
  DreamFile->SetAnalysisFile(filename, prefix, addon);

  DreamCF* CF_pp = new DreamCF();
  DreamPair* pp = new DreamPair("Part",0.2,0.4);
  DreamPair* ApAp = new DreamPair("AntiPart",0.2,0.4);

  DreamCF* CF_pL = new DreamCF();
  DreamPair* pL = new DreamPair("Part",0.2,0.4);
  DreamPair* ApAL = new DreamPair("AntiPart",0.2,0.4);

  DreamCF* CF_LL = new DreamCF();
  DreamPair* LL = new DreamPair("Part",0.2,0.4);
  DreamPair* ALAL = new DreamPair("AntiPart",0.2,0.4);

  DreamCF* CF_pXi = new DreamCF();
  DreamPair* pXi = new DreamPair("Part",0.2,0.4);
  DreamPair* ApAXi = new DreamPair("AntiPart",0.2,0.4);

  std::cout << "=========================" << std::endl;
  std::cout << "========Pair Set=========" << std::endl;
  std::cout << "=========================" << std::endl;
  pp->SetPair(DreamFile->GetPairDistributions(0, 0, ""));
  ApAp->SetPair(DreamFile->GetPairDistributions(1, 1, ""));

  pL->SetPair(DreamFile->GetPairDistributions(0, 2, ""));
  ApAL->SetPair(DreamFile->GetPairDistributions(1, 3, ""));
  TH1F* hSE = DreamFile->GetPairDistributions(0, 2, "")->GetSEDist();
  int NumPairs = hSE->Integral(1,hSE->FindBin(0.312-1e-6));
  const static int DefNumPairs = hSE->Integral(1,hSE->FindBin(0.312-1e-6));
  static int VariationNr = -1;
  VariationNr++;
  Dev_pL = double(NumPairs-DefNumPairs)/double(DefNumPairs);
  printf("Addon: %s\n",addon);
  printf(" NumPairs = %u k\n",NumPairs/1000);
  printf(" Deviation: ");
  if(Dev_pL<0.2) printf("%.1f%%\n",100.*Dev_pL);
  else printf("\033[1;31m%.1f%%\033[0m\n",100.*Dev_pL);
  //usleep(2000e3);

  LL->SetPair(DreamFile->GetPairDistributions(2, 2, ""));
  ALAL->SetPair(DreamFile->GetPairDistributions(3, 3, ""));

  pXi->SetPair(DreamFile->GetPairDistributions(0, 4, ""));
  ApAXi->SetPair(DreamFile->GetPairDistributions(1, 5, ""));

  std::cout << "=========================" << std::endl;
  std::cout << "======Pair Shifted=======" << std::endl;
  std::cout << "=========================" << std::endl;
  pp->ShiftForEmpty(pp->GetPair());
  ApAp->ShiftForEmpty(ApAp->GetPair());

  pL->ShiftForEmpty(pL->GetPair());
  ApAL->ShiftForEmpty(ApAL->GetPair());

  LL->ShiftForEmpty(LL->GetPair());
  ALAL->ShiftForEmpty(ALAL->GetPair());

  pXi->ShiftForEmpty(pXi->GetPair());
  ApAXi->ShiftForEmpty(ApAXi->GetPair());

  std::cout << "=========================" << std::endl;
  std::cout << "====Pair Fix Shifted=====" << std::endl;
  std::cout << "=========================" << std::endl;
  pp->FixShift(pp->GetPairShiftedEmpty(0), ApAp->GetPairShiftedEmpty(0),
               ApAp->GetFirstBin());
  ApAp->FixShift(ApAp->GetPairShiftedEmpty(0), pp->GetPairShiftedEmpty(0),
                 pp->GetFirstBin());

  pL->FixShift(pL->GetPairShiftedEmpty(0), ApAL->GetPairShiftedEmpty(0),
               ApAL->GetFirstBin());
  ApAL->FixShift(ApAL->GetPairShiftedEmpty(0), pL->GetPairShiftedEmpty(0),
                 pL->GetFirstBin());
  LL->FixShift(LL->GetPairShiftedEmpty(0), ALAL->GetPairShiftedEmpty(0),
               ALAL->GetFirstBin());
  ALAL->FixShift(ALAL->GetPairShiftedEmpty(0), LL->GetPairShiftedEmpty(0),
                 LL->GetFirstBin());

  pXi->FixShift(pXi->GetPairShiftedEmpty(0), ApAXi->GetPairShiftedEmpty(0),
                ApAXi->GetFirstBin());
  ApAXi->FixShift(ApAXi->GetPairShiftedEmpty(0), pXi->GetPairShiftedEmpty(0),
                  pXi->GetFirstBin());

  std::cout << "=========================" << std::endl;
  std::cout << "==Rebinning & Weighting==" << std::endl;
  std::cout << "=========================" << std::endl;

  std::vector<int> rebinVec = {{1,2,3,4,5}};
  for (size_t iReb = 0; iReb < rebinVec.size(); ++iReb) {
    std::cout << "==Rebinning==" << std::endl;
    pL->Rebin(pL->GetPairFixShifted(0), rebinVec[iReb]);
    LL->Rebin(LL->GetPairFixShifted(0), rebinVec[iReb]);
    pXi->Rebin(pXi->GetPairFixShifted(0), rebinVec[iReb]);
    std::cout << "==Weighting==" << std::endl;
    pL->ReweightMixedEvent(pL->GetPairRebinned(iReb), 0.2, 0.9);
    LL->ReweightMixedEvent(LL->GetPairRebinned(iReb), 0.2, 0.9);
    pXi->ReweightMixedEvent(pXi->GetPairRebinned(iReb), 0.2, 0.9);
    std::cout << "==Rebinning==" << std::endl;
    ApAL->Rebin(ApAL->GetPairFixShifted(0), rebinVec[iReb]);
    ALAL->Rebin(ALAL->GetPairFixShifted(0), rebinVec[iReb]);
    ApAXi->Rebin(ApAXi->GetPairFixShifted(0), rebinVec[iReb]);
    std::cout << "==Weighting==" << std::endl;
    ApAL->ReweightMixedEvent(ApAL->GetPairRebinned(iReb), 0.2, 0.9);
    ALAL->ReweightMixedEvent(ALAL->GetPairRebinned(iReb), 0.2, 0.9);
    ApAXi->ReweightMixedEvent(ApAXi->GetPairRebinned(iReb), 0.2, 0.9);
  }
  pp->ReweightMixedEvent(pp->GetPairFixShifted(0), 0.2, 0.9);
  ApAp->ReweightMixedEvent(ApAp->GetPairFixShifted(0), 0.2, 0.9);

  pL->Rebin(pL->GetPair(), 4);
  pL->Rebin(pL->GetPair(), 5);
  ApAL->Rebin(ApAL->GetPair(), 4);
  ApAL->Rebin(ApAL->GetPair(), 5);
  LL->Rebin(LL->GetPair(), 4);
  LL->Rebin(LL->GetPair(), 5);
  ALAL->Rebin(ALAL->GetPair(), 4);
  ALAL->Rebin(ALAL->GetPair(), 5);
  pXi->Rebin(pXi->GetPair(), 4);
  pXi->Rebin(pXi->GetPair(), 5);
  ApAXi->Rebin(ApAXi->GetPair(), 4);
  ApAXi->Rebin(ApAXi->GetPair(), 5);

  std::cout << "=========================" << std::endl;
  std::cout << "=========CFs=============" << std::endl;
  std::cout << "=========================" << std::endl;

  //TString foldername = filename;
  //foldername.ReplaceAll("AnalysisResults.root", "");

  CF_pp->SetPairs(pp, ApAp);
  CF_pp->GetCorrelations();
  if(strcmp(addon,"")==0) CF_pp->WriteOutput(Form("%sCFOutput_pp.root", foldername));
  else CF_pp->WriteOutput(Form("%sCFOutput_pp_%s.root", foldername, addon));

  CF_pL->SetPairs(pL, ApAL);
  CF_pL->GetCorrelations();
  if(strcmp(addon,"")==0) CF_pL->WriteOutput(Form("%sCFOutput_pL.root", foldername));
  else CF_pL->WriteOutput(Form("%sCFOutput_pL_%s.root", foldername, addon));

  CF_LL->SetPairs(LL, ALAL);
  CF_LL->GetCorrelations();
  if(strcmp(addon,"")==0) CF_LL->WriteOutput(Form("%sCFOutput_LL.root", foldername));
  else CF_LL->WriteOutput(Form("%sCFOutput_LL_%s.root", foldername, addon));

  CF_pXi->SetPairs(pXi, ApAXi);
  CF_pXi->GetCorrelations();
  if(strcmp(addon,"")==0) CF_pXi->WriteOutput(Form("%sCFOutput_pXi.root", foldername));
  else CF_pXi->WriteOutput(Form("%sCFOutput_pXi_%s.root", foldername, addon));
}

void GetCorrelationsBbarB(const char* foldername, const char* filename, const char* prefix,
                     const char* addon = "") {
  //gStyle->SetOptStat(0);
  ReadDreamFile* DreamFile = new ReadDreamFile(6, 6);
  DreamFile->SetAnalysisFile(filename, prefix, addon);

  DreamCF* CF_pAp_App = new DreamCF();
  DreamPair* pAp = new DreamPair("PartAntiPart",0.2,0.4);
//  DreamPair* App = new DreamPair("AntiPartPart");

  DreamCF* CF_pAL_ApL = new DreamCF();
  DreamPair* pAL = new DreamPair("PartAntiPart",0.2,0.4);
  DreamPair* ApL = new DreamPair("AntiPartPart",0.2,0.4);

  DreamCF* CF_ALL_LAL = new DreamCF();
  DreamPair* LAL = new DreamPair("PartAntiPart",0.2,0.4);
//  DreamPair* ALL = new DreamPair("AntiPartPart");

  DreamCF* CF_pAXi_ApXi = new DreamCF();
  DreamPair* pAXi = new DreamPair("PartAntiPart",0.2,0.4);
  DreamPair* ApXi = new DreamPair("AntiPartPart",0.2,0.4);

  std::cout << "=========================" << std::endl;
  std::cout << "========Pair Set=========" << std::endl;
  std::cout << "=========================" << std::endl;
  pAp->SetPair(DreamFile->GetPairDistributions(0, 1, ""));
//  App->SetPair(DreamFile->GetPairDistributions(1, 0, ""));

  pAL->SetPair(DreamFile->GetPairDistributions(0, 3, ""));
  ApL->SetPair(DreamFile->GetPairDistributions(1, 2, ""));

  LAL->SetPair(DreamFile->GetPairDistributions(2, 3, ""));
//  ALL->SetPair(DreamFile->GetPairDistributions(3, 2, ""));

  pAXi->SetPair(DreamFile->GetPairDistributions(0, 5, ""));
  ApXi->SetPair(DreamFile->GetPairDistributions(1, 4, ""));

  std::cout << "=========================" << std::endl;
  std::cout << "======Pair Shifted=======" << std::endl;
  std::cout << "=========================" << std::endl;
  pAp->ShiftForEmpty(pAp->GetPair());
//  App->ShiftForEmpty(App->GetPair());

  pAL->ShiftForEmpty(pAL->GetPair());
  ApL->ShiftForEmpty(ApL->GetPair());

  LAL->ShiftForEmpty(LAL->GetPair());
//  ALL->ShiftForEmpty(ALL->GetPair());

  pAXi->ShiftForEmpty(pAXi->GetPair());
  ApXi->ShiftForEmpty(ApXi->GetPair());

  std::cout << "=========================" << std::endl;
  std::cout << "====Pair Fix Shifted=====" << std::endl;
  std::cout << "=========================" << std::endl;
  //Fix shift singe pair cfs anyway to ensure competability!
  pAp->FixShift(pAp->GetPairShiftedEmpty(0), pAp->GetPairShiftedEmpty(0),
                pAp->GetFirstBin());
  pAL->FixShift(pAL->GetPairShiftedEmpty(0), ApL->GetPairShiftedEmpty(0),
               ApL->GetFirstBin());
  ApL->FixShift(ApL->GetPairShiftedEmpty(0), pAL->GetPairShiftedEmpty(0),
                 pAL->GetFirstBin());
  LAL->FixShift(LAL->GetPairShiftedEmpty(0), LAL->GetPairShiftedEmpty(0),
                LAL->GetFirstBin());

  pAXi->FixShift(pAXi->GetPairShiftedEmpty(0), ApXi->GetPairShiftedEmpty(0),
                ApXi->GetFirstBin());
  ApXi->FixShift(ApXi->GetPairShiftedEmpty(0), pAXi->GetPairShiftedEmpty(0),
                  pAXi->GetFirstBin());

  std::cout << "=========================" << std::endl;
  std::cout << "==Rebinning & Weighting==" << std::endl;
  std::cout << "=========================" << std::endl;

  for (int iReb = 4; iReb < 6; ++iReb) {
    std::cout << "==Rebinning==" << std::endl;
    pAL->Rebin(pAL->GetPairFixShifted(0), iReb);
    pAXi->Rebin(pAXi->GetPairFixShifted(0), iReb);
    LAL->Rebin(LAL->GetPairFixShifted(0), iReb);
    std::cout << "==Weighting==" << std::endl;
    pAL->ReweightMixedEvent(pAL->GetPairRebinned(iReb - 4), 0.2, 0.9);
    LAL->ReweightMixedEvent(LAL->GetPairRebinned(iReb - 4), 0.2, 0.9);
    pAXi->ReweightMixedEvent(pAXi->GetPairRebinned(iReb - 4), 0.2, 0.9);
    std::cout << "==Rebinning==" << std::endl;
    ApL->Rebin(ApL->GetPairFixShifted(0), iReb);
    ApXi->Rebin(ApXi->GetPairFixShifted(0), iReb);
    std::cout << "==Weighting==" << std::endl;
    ApL->ReweightMixedEvent(ApL->GetPairRebinned(iReb - 4), 0.2, 0.9);
    ApXi->ReweightMixedEvent(ApXi->GetPairRebinned(iReb - 4), 0.2, 0.9);
  }
  pAp->ReweightMixedEvent(pAp->GetPairShiftedEmpty(0), 0.2, 0.9);

  pAL->Rebin(pAL->GetPair(), 4);
  pAL->Rebin(pAL->GetPair(), 5);
  ApL->Rebin(ApL->GetPair(), 4);
  ApL->Rebin(ApL->GetPair(), 5);
  LAL->Rebin(LAL->GetPair(), 4);
  LAL->Rebin(LAL->GetPair(), 5);
  pAXi->Rebin(pAXi->GetPair(), 4);
  pAXi->Rebin(pAXi->GetPair(), 5);
  ApXi->Rebin(ApXi->GetPair(), 4);
  ApXi->Rebin(ApXi->GetPair(), 5);

  std::cout << "=========================" << std::endl;
  std::cout << "=========CFs=============" << std::endl;
  std::cout << "=========================" << std::endl;

  //TString foldername = filename;
  //foldername.ReplaceAll("AnalysisResults.root", "");
  std::cout << "$PWD " << foldername << std::endl;
  std::cout << "pp CF \n";
  std::cout << "Set Pair \n";
  CF_pAp_App->SetPairs(pAp, nullptr);
  std::cout << "Get CF \n";
  CF_pAp_App->GetCorrelations();
  std::cout << "Write Output \n";
  if(strcmp(addon,"")==0) CF_pAp_App->WriteOutput(Form("%sCFOutput_pAp_App.root", foldername));
  else CF_pAp_App->WriteOutput(Form("%sCFOutput_pAp_App_%s.root", foldername, addon));

  std::cout << "pL CF \n";
  CF_pAL_ApL->SetPairs(pAL, ApL);
  CF_pAL_ApL->GetCorrelations();
  if(strcmp(addon,"")==0) CF_pAL_ApL->WriteOutput(Form("%sCFOutput_pAL_ApL.root", foldername));
  else CF_pAL_ApL->WriteOutput(Form("%sCFOutput_pAL_ApL_%s.root", foldername, addon));

  std::cout << "LL CF \n";
  CF_ALL_LAL->SetPairs(LAL, nullptr);
  CF_ALL_LAL->GetCorrelations();
  if(strcmp(addon,"")==0) CF_ALL_LAL->WriteOutput(Form("%sCFOutput_LAL_ALL.root", foldername));
  else CF_ALL_LAL->WriteOutput(Form("%sCFOutput_LAL_ALL_%s.root", foldername, addon));

  std::cout << "pXi CF \n";
  CF_pAXi_ApXi->SetPairs(pAXi, ApXi);
  CF_pAXi_ApXi->GetCorrelations();
  if(strcmp(addon,"")==0) CF_pAXi_ApXi->WriteOutput(Form("%sCFOutput_pAXi_ApXi.root", foldername));
  else CF_pAXi_ApXi->WriteOutput(Form("%sCFOutput_pAXi_ApXi_%s.root", foldername, addon));
}

void GetQADistributions(const char* PairName, DreamDist* PairOrg,
                        DreamDist* PairRew) {
  gStyle->SetOptStat(false);
  DreamPlot::SetStyle();

  TString CanvasName = Form("Can%s", PairName);
  TCanvas* c1 = new TCanvas(CanvasName.Data(), CanvasName.Data(), 2000, 1000);
  c1->Divide(2, 2);

  TH1F* MERaw = PairOrg->GetMEDist();
  TH2F* SEMultRaw = PairOrg->GetSEMultDist();
  TH2F* MEMultRaw = PairOrg->GetMEMultDist();
  TH1F* CFRaw = PairOrg->GetCF();

  TH1F* MERew = PairRew->GetMEDist();
  TH2F* MEMultRew = PairRew->GetMEMultDist();
  TH1F* CFRew = PairRew->GetCF();
  c1->cd(1);
  TLegend leg = TLegend();
  TH1F* SEMultProjRaw = (TH1F*) SEMultRaw->ProjectionY("SEMultProjRaw",
                                                       SEMultRaw->FindBin(0.2),
                                                       SEMultRaw->FindBin(0.4));
  SEMultProjRaw->SetTitle("; Multiplicity Bin; Normalized counts");
  leg.AddEntry(SEMultProjRaw, "SE", "lep");
  double Scale_SEMultProjRaw = 1. / SEMultProjRaw->Integral();
  SEMultProjRaw->Scale(Scale_SEMultProjRaw);
  TH1F* MEMultProjRaw = (TH1F*) MEMultRaw->ProjectionY("MEMultProjRaw",
                                                       MEMultRaw->FindBin(0.2),
                                                       MEMultRaw->FindBin(0.4));
  double Scale_MEMultProjRaw = 1. / MEMultProjRaw->Integral();
  MEMultProjRaw->Scale(Scale_MEMultProjRaw);
  leg.AddEntry(MEMultProjRaw, "ME Unweighted", "lep");

  TH1F* MEMultProjRew = (TH1F*) MEMultRew->ProjectionY("MEMultProjRew",
                                                       MEMultRew->FindBin(0.2),
                                                       MEMultRew->FindBin(0.4));
  double Scale_MEMultProjRew = 1. / MEMultProjRew->Integral();
  MEMultProjRew->Scale(Scale_MEMultProjRew);
  leg.AddEntry(MEMultProjRew, "ME Weighted", "lep");

  DreamPlot::SetStyleHisto(SEMultProjRaw, 21, kGreen + 2);
  DreamPlot::SetStyleHisto(MEMultProjRaw, 22, kRed + 2);
  DreamPlot::SetStyleHisto(MEMultProjRew, 23, kBlue + 2);
  SEMultProjRaw->DrawCopy();
  SEMultProjRaw->SetStats(false);
  MEMultProjRaw->DrawCopy("same");
  MEMultProjRew->DrawCopy("Same");
  leg.Draw("same");

  c1->cd(2);
  TH1F* MERawNorm = (TH1F*) MERaw->Clone(Form("%s_Norm", MERaw->GetName()));
  double Scale_MERawNorm = 1. / MERawNorm->Integral();
  MERawNorm->Scale(Scale_MERawNorm);
  MERawNorm->SetTitle("; #it{k}* (GeV/#it{c}); Normalized counts");
  MERawNorm->Draw();
  MERawNorm->SetStats(false);
  TH1F* MERewNorm = (TH1F*) MERew->Clone(Form("%s_Norm", MERew->GetName()));
  double Scale_MERewNorm = 1. / MERewNorm->Integral();
  MERewNorm->Scale(Scale_MERewNorm);
  MERewNorm->Draw("same");
  DreamPlot::SetStyleHisto(MERawNorm, 22, kRed + 2);
  DreamPlot::SetStyleHisto(MERewNorm, 23, kBlue + 2);
  TLegend leg2 = TLegend();
  leg2.AddEntry(MERawNorm, "ME Unweighted", "lep");
  leg2.AddEntry(MERewNorm, "ME Weighted", "lep");
  leg2.Draw("same");

  c1->cd(3);
  CFRaw->SetTitle(";#it{k}* (GeV/#it{c}); C(#it{k}*)");
  CFRaw->Draw();
  CFRaw->SetStats(false);
  CFRew->Draw("same");
  DreamPlot::SetStyleHisto(CFRaw, 22, kRed + 2);
  DreamPlot::SetStyleHisto(CFRew, 23, kBlue + 2);
  TLegend leg3 = TLegend();
  leg3.AddEntry(CFRaw, "CF Unweighted", "lep");
  leg3.AddEntry(CFRew, "CF Weighted", "lep");
  leg3.Draw("same");

  c1->cd(4);

  TH1F* Ratio = (TH1F*) CFRaw->Clone(Form("SEMEReweightingRatio_%s",PairName));
  if (Ratio->Divide(CFRew)) {
    Ratio->Draw();
    Ratio->SetStats(false);
    DreamPlot::SetStyleHisto(Ratio, 20, kBlue + 2);
    Ratio->SetTitle(
        "; #it{k}* (GeV/#it{c}); C(#it{k}*)_{raw} / C(#it{k}*)_{reweighted}");
    Ratio->Write();
  }
  c1->Write();
  delete c1;
  delete SEMultProjRaw;
  delete MEMultProjRaw;
  delete MEMultProjRew;
  delete MERawNorm;
  delete MERewNorm;
  return;
}
void ReweightingQA(TList* PairList) {
  TList* UntouchedPairList = (TList*) PairList->FindObject("PairFixShifted");
  TList* PairRebinnedList = (TList*) PairList->FindObject("PairRebinned");
  TList* PairReweightedList = (TList*) PairList->FindObject("PairReweighted");
  const int nEntries = PairReweightedList->GetEntries() / 5;
  for (int iQA = 0; iQA < nEntries; ++iQA) {
    DreamDist* PairRew = new DreamDist();
    PairRew->SetSEDist((TH1F*) PairReweightedList->At(5 * iQA), "_");
    PairRew->SetSEMultDist((TH2F*) PairReweightedList->At(5 * iQA + 1), "_");
    PairRew->SetMEDist((TH1F*) PairReweightedList->At(5 * iQA + 2), "_");
    PairRew->SetMEMultDist((TH2F*) PairReweightedList->At(5 * iQA + 3), "_");
    PairRew->Calculate_CF(0.2,0.4);
    //Get the corresponding unmodified pair, which was just reweighted!
    DreamDist* Pair=nullptr;
    TString SEName = PairRew->GetSEDist()->GetName();
    TString PartName = SEName(SEName.Index("Particle"), 19);
    TString CanName = PartName;

    TString test = "_Rebinned_";
    int iStart = SEName.Index(test.Data());

    std::cout << PartName.Data() << std::endl;
    if (iStart > 0) {
      iStart += test.Length();
      std::cout << SEName.Data() << std::endl;
      TString RebinName = SEName(iStart, 1);
      CanName += "_Rebinned_";
      CanName += RebinName;
      if (PairRebinnedList->FindObject(
          Form("SEDist_%s_clone_Shifted_FixShifted_Rebinned_%s", PartName.Data(),
               RebinName.Data()))) {
        Pair = new DreamDist();
        Pair->SetSEDist(
            (TH1F*) PairRebinnedList->FindObject(
                Form("SEDist_%s_clone_Shifted_FixShifted_Rebinned_%s",
                     PartName.Data(), RebinName.Data())),
            "_");
        Pair->SetSEMultDist(
            (TH2F*) PairRebinnedList->FindObject(
                Form("SEMultDist_%s_clone_Shifted_FixShifted_Rebinned_%s",
                     PartName.Data(), RebinName.Data())),
            "_");
        Pair->SetMEDist(
            (TH1F*) PairRebinnedList->FindObject(
                Form("MEDist_%s_clone_Shifted_FixShifted_Rebinned_%s",
                     PartName.Data(), RebinName.Data())),
            "_");
        Pair->SetMEMultDist(
            (TH2F*) PairRebinnedList->FindObject(
                Form("MEMultDist_%s_clone_Shifted_FixShifted_Rebinned_%s",
                     PartName.Data(), RebinName.Data())),
            "_");
        Pair->Calculate_CF(0.2,0.4);
      } else {
        std::cout
            << "===========" << '\n' << "==Missing==" << '\n' << "==========="
            << '\n'
            << Form("SEDist_%s_clone_Rebinned_%s", PartName.Data(), RebinName.Data())
            << std::endl;
      }
    } else {
      //else we don't rebinn, e.g. pp case, then take the original CF
      std::cout << Form("SEDist_%s", PartName.Data()) << std::endl;
      if (UntouchedPairList->FindObject(
          Form("SEDist_%s_clone_Shifted_FixShifted", PartName.Data()))) {
        Pair = new DreamDist();
        Pair->SetSEDist(
            (TH1F*) UntouchedPairList->FindObject(
                Form("SEDist_%s_clone_Shifted_FixShifted", PartName.Data())),
            "_");
        Pair->SetSEMultDist(
            (TH2F*) UntouchedPairList->FindObject(
                Form("SEMultDist_%s_clone_Shifted_FixShifted", PartName.Data())),
            "_");
        Pair->SetMEDist(
            (TH1F*) UntouchedPairList->FindObject(
                Form("MEDist_%s_clone_Shifted_FixShifted", PartName.Data())),
            "_");
        Pair->SetMEMultDist(
            (TH2F*) UntouchedPairList->FindObject(
                Form("MEMultDist_%s_clone_Shifted_FixShifted", PartName.Data())),
            "_");
        Pair->Calculate_CF(0.2,0.4);
      } else {
        std::cout << "===========" << '\n' << "==Missing==" << '\n'
                  << "===========" << '\n'
                  << Form("SEDist_%s_clone_Shifted_FixShifted", PartName.Data())
                  << std::endl;
      }
    }
    std::cout << CanName.Data() << std::endl;
    if (Pair && PairRew) {
      GetQADistributions(CanName.Data(), Pair, PairRew);
    }
    if (Pair)
      delete Pair;
    if (PairRew)
      delete PairRew;
  }
  return;
}

void METoSEReweighting(const char* foldername, const char* addon) {
  const char* filenames[8] = { "pp", "pXi", "pL", "LL", "pAp_App", "pAL_ApL",
      "LAL_ALL", "pAXi_ApXi" };
  for (int iFile = 0; iFile < 8; ++iFile) {
    TString FileName;
    if(strcmp(addon,"")!=0) FileName = Form("%sCFOutput_%s_%s.root", foldername, filenames[iFile], addon);
    else FileName = Form("%sCFOutput_%s.root", foldername, filenames[iFile]);
    std::cout << FileName.Data() << std::endl;
    TFile* file = TFile::Open(FileName, "update");
    TList* PairDist = (TList*) file->Get("PairDist");
    if (PairDist) {
      ReweightingQA(PairDist);
    } else {
      file->ls();
    }
    TList* AntiPairDist = (TList*) file->Get("AntiPairDist");
    if (AntiPairDist) {
      ReweightingQA(AntiPairDist);
    }
//    TH1F* CFRaw;
//    TH1F* CFRew;
//    if (iFile == 0) {
//      CFRaw=(TH1F*)file->Get("hCkTotNormWeight");
//      CFRew = (TH1F*)file->Get("hCkReweighted_0");
//    }
//    else{
//      CFRaw=(TH1F*)file->Get("hCkReweighted_3");
//      CFRew = (TH1F*)file->Get("hCkReweighted_1");
//    }
//
//    TH1F* Ratio = (TH1F*)CFRaw->Clone(Form("Ratio_%i",iFile));
//    Ratio->Divide(CFRew);
//    TCanvas *c1=new TCanvas(Form("CanRatio%i",iFile),Form("CanRatio%i",iFile),2000,1000);
//    c1->Divide(2,1);
//    c1->cd(1);
//    CFRaw->SetLineColor(3);
//    CFRaw->DrawCopy();
//    CFRew->SetLineColor(4);
//    CFRew->DrawCopy("same");
//    c1->cd(2);
//    Ratio->DrawCopy();
//    c1->Write();
//    delete Ratio;
  }
}

int ExecuteCFmT_pp(int argc, char* argv[]) {
    const char* foldername = argv[0];
  const char* filename = argv[1];
  const char* CalibName = argv[2];
  const char* prefix = argv[3];
  const char* addon = (argv[4]) ? argv[4] : "";

  ReadDreamFile* DreamFile = new ReadDreamFile(6, 6);

  //TString foldername = filename;
  //foldername.ReplaceAll("AnalysisResults.root", "");
  DreamKayTee* mTppDists;
  DreamFile->ReadmTHistos(filename, prefix, addon);

  mTppDists = DreamFile->GetmTPairDistributions(0, 0, 1, 1);

  //std::vector<float> mTppBins = TString(filename)=="/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample9HM/AnalysisResults.root"?{ 0.9 ,1.07,1.13,1.19,1.25,1.37,1.55,1.97,4.5 }:
  //                              TString(filename)=="/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample9/AnalysisResults.root"?{ 0.9 ,1.1,1.25,1.5,4.5 }:
  //                                  {0,1,2,3,4,5};
    std::vector<float> mTppBins;
    if(TString(filename)=="/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample10HM/AnalysisResults.root"||
       TString(filename)=="/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample10HM_ver2/AnalysisResults.root"){
        //mTppBins = { 0.9 ,1.07,1.13,1.19,1.25,1.37,1.55,1.97,4.5 };
        mTppBins = { 1.02,1.14,1.199,1.26,1.38,1.56,1.86,4.5 };
        //mTppBins = { 1.02,1.14,1.2,1.26,1.38,1.56,1.86,4.5 };
    }
    else if(TString(filename)=="/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample11HM/AnalysisResults.root"){
        mTppBins = { 0.9 ,1.07,1.13,1.19,1.25,1.37,1.55,1.97,4.5 };
        //mTppBins = { 1.02,1.14,1.2,1.26,1.38,1.56,1.86,4.5 };
    }
    else if(TString(filename)=="/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample9/AnalysisResults.root"){
        mTppBins = { 0.9 ,1.1,1.25,1.5,4.5 };
    }
    else if(TString(filename)=="/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pPb_5TeV/Sample11/AnalysisResults.root"){
        //mTppBins = { 0.9 ,1.1,1.25,1.5,4.5 };
        mTppBins = { 0.9 ,1.1,1.2,1.4,1.9,4.5 };
    }
    else{
        mTppBins = { 0,1,2,3,4,5 };
    }

//  std::vector<float> mTBins = { 1.01,1.13,1.25,1.55,4.5 };
  TString CalibPP = Form("%s/CFOutput_pp.root",CalibName);
  mTppDists->SetSEMEReweightingRatio(CalibPP,"pp");
  mTppDists->SetKayTeeBins(mTppBins);
  mTppDists->SetNormalization(0.2,0.4);
  mTppDists->ObtainTheCorrelationFunction(gSystem->pwd(), prefix, "pp");

  DreamKayTee* mTpLDists;
  DreamFile->ReadmTHistos(filename, prefix, addon);
  mTpLDists = DreamFile->GetmTPairDistributions(0, 2, 1, 3);

  std::vector<float> mTpLBins =

  { 1.08, 1.26, 1.32, 1.44, 1.62, 1.68, 4.5 };

  //TString CalibPL = Form("%s/CFOutput_pL.root",CalibName);
  //mTpLDists->SetSEMEReweightingRatio(CalibPL, "pL");
  //mTpLDists->SetKayTeeBins(mTpLBins);
  //mTpLDists->SetNormalization(0.2,0.4);
  //mTpLDists->ObtainTheCorrelationFunction(gSystem->pwd(), prefix, "pL");

  return 1;
}

int ExecuteCFmT_pL(int argc, char* argv[]) {
    const char* foldername = argv[0];
  const char* filename = argv[1];
  const char* CalibName = argv[2];
  const char* prefix = argv[3];
  const char* addon = (argv[4]) ? argv[4] : "";

  ReadDreamFile* DreamFile = new ReadDreamFile(6, 6);

  DreamKayTee* mTpLDists;
  DreamFile->ReadmTHistos(filename, prefix, addon);
  mTpLDists = DreamFile->GetmTPairDistributions(0, 2, 1, 3);

    std::vector<float> mTpLBins;
    if(TString(filename)=="/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample12HM/AnalysisResults.root"||
       TString(filename)=="/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample13HM/AnalysisResults.root"){
        //mTpLBins = { 1.08, 1.26, 1.32, 1.44, 1.68, 4.5 };
        mTpLBins = { 1.08, 1.26, 1.32, 1.44, 1.65, 1.9, 4.5 };
    }
    else if(TString(filename)=="/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample9/AnalysisResults.root"){
        mTpLBins = { 1.08, 1.29, 1.44, 1.65, 4.5 };
    }
    else if(TString(filename)=="/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pPb_5TeV/Sample11/AnalysisResults.root"){
        mTpLBins = { 1.08, 1.29, 1.44, 1.7, 4.5 };
    }
    else{
        mTpLBins = { 0,1,2,3,4,5 };
    }

//  std::vector<float> mTBins = { 1.01,1.13,1.25,1.55,4.5 };
  TString CalibPP;
  if(strcmp(addon,"")!=0) CalibPP = Form("%sCFOutput_pL_%s.root",CalibName,addon);
  else CalibPP = Form("%sCFOutput_pL.root",CalibName);
  printf("CalibPP=%s\n",CalibPP.Data());
  mTpLDists->SetSEMEReweightingRatio(CalibPP,"pL");
  mTpLDists->SetKayTeeBins(mTpLBins);
  mTpLDists->SetNormalization(0.2,0.4);
  mTpLDists->SetRebin(3);
  mTpLDists->ObtainTheCorrelationFunction(gSystem->pwd(), prefix, "pL", addon);

  return 1;
}

int ExecuteCFDream(int argc, char* argv[], double& Dev_pL) {
  const char* foldername = argv[0];
  const char* filename = argv[1];
  const char* prefix = argv[3];
  const char* addon = (argv[4]) ? argv[4] : "";
//printf("gSystem->pwd()=%s\n",gSystem->pwd());
  GetCorrelations(foldername, filename, prefix, addon, Dev_pL);
  GetCorrelationsBbarB(foldername, filename, prefix, addon);
  //TString foldername = filename;
  //foldername.ReplaceAll("AnalysisResults.root", "");
  METoSEReweighting(foldername,addon);

  return 1;
}


void QA_for_pp_mT(){

    //const TString TotalCF_FileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/DimiExecuteCFmT/CFOutput_pp.root";
    const TString InputFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/DimiExecuteCFmT/CFOutputALL_mT_pp_HM.root";
    const TString HistoNameWeights = "hWeights_0_ME";
    const unsigned Num_mT_Bins = 8;
    const TString BaseName = "hCkTotNormWeightMeV";

    const TString OutputFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/DimiExecuteCFmT/QA_for_pp_mT.root";

    TFile* fIn = new TFile(InputFileName,"read");

    TH2F* hWeightsME = (TH2F*)fIn->Get(HistoNameWeights);
    TH1F** hCk = new TH1F* [Num_mT_Bins];

    for(unsigned umT=0; umT<Num_mT_Bins; umT++){
        hCk[umT] = (TH1F*)fIn->Get(BaseName+"_mTBin_"+umT);
    }

    TH1F* hTotCk = (TH1F*) hCk[0]->Clone("hTotCk");

    double BinContent=0;
    double BinError=0;
    double Weight;
    for(unsigned uMom=0; uMom<hTotCk->GetNbinsX(); uMom++){
        BinContent=0;
        BinError=0;
        double TotalWeight = 0;
        for(unsigned umT=0; umT<Num_mT_Bins; umT++){

            Weight = hWeightsME->GetBinContent(uMom+1,umT+1);
            TotalWeight+=Weight;
            BinContent += Weight*hCk[umT]->GetBinContent(uMom+1);
            BinError += pow(Weight*hCk[umT]->GetBinError(uMom+1),2.);
        }
        printf("uMom=%u; TotalWeight=%.3f\n",uMom,TotalWeight);
        hTotCk->SetBinContent(uMom+1,BinContent);
        hTotCk->SetBinError(uMom+1,sqrt(BinError));
    }

    TFile* fOut = new TFile(OutputFileName,"recreate");
    hTotCk->Write();
    //for(unsigned umT=0; umT<Num_mT_Bins; umT++){
        //delete hCk[umT];
    //}
    delete [] hCk;
    delete hTotCk;
    delete fIn;
    delete fOut;
}



void pp_FitDiff_CompareToTotal(char* argv[], const TString& DataSample, const TString& SourceType, const TString& FittingMode_pp){

    TString InputFileNameDiff = TString(argv[2])+"CFOutputALL_mT_pp_"+argv[3]+".root";
    TString InputFileNameTot = TString(argv[2])+"CFOutput_pp.root";
    const TString HistoNameWeights = "hWeights_0_ME";
    const TString HistoNameAvg_mT = "AveragemT";
    const TString HistoNameAvg_mT_k = "hAvgKayTee_0_SE";

    const TString BaseName = "hCkTotNormWeightMeV";
    const TString MainHistoName = "hCk_ReweightedMeV_0";

    const TString OutputFolder = argv[0];

    TFile* fInDiff = new TFile(InputFileNameDiff,"read");
//printf("InputFileNameDiff=%s fInDiff=%p\n",InputFileNameDiff.Data(),fInDiff);

    TH2F* hWeightsME = (TH2F*)fInDiff->Get(HistoNameWeights);
    TGraph* Avg_mT = (TGraph*)fInDiff->Get(HistoNameAvg_mT);
    TH1F* hAvgKayTee = (TH1F*)fInDiff->Get(HistoNameAvg_mT_k);
    const unsigned Num_mT_Bins = Avg_mT->GetN();

    TH1F** hCk = new TH1F* [Num_mT_Bins];
    for(unsigned umT=0; umT<Num_mT_Bins; umT++){
        hCk[umT] = (TH1F*)fInDiff->Get(BaseName+"_mTBin_"+umT);
        hCk[umT]->SetName(TString::Format("hCk_%u",umT));
    }

    TH1F* hTotCk;
    TFile* fInTot = new TFile(InputFileNameTot,"read");
    hTotCk = (TH1F*)fInTot->Get(MainHistoName);
    hTotCk->SetName("hTotCk");

    DLM_CommonAnaFunctions AnalysisObject;

    double lambda_pp[4];
    AnalysisObject.SetUpLambdaPars_pp(DataSample,0,lambda_pp);
    double lambda_pL[5];
    AnalysisObject.SetUpLambdaPars_pL(DataSample,0,0,lambda_pL);

    double* MomBins_pp = NULL;
    double* FitRegion_pp = NULL;
    unsigned NumMomBins_pp;

    AnalysisObject.SetUpBinning_pp(DataSample,NumMomBins_pp,MomBins_pp,FitRegion_pp);
    TH2F* hResolution_pp = AnalysisObject.GetResolutionMatrix(DataSample,"pp");
    TH2F* hResidual_pp_pL = AnalysisObject.GetResidualMatrix("pp","pLambda");

    double lam_pp[4];
    double lam_pL[5];

    AnalysisObject.SetUpLambdaPars_pp(DataSample,0,lam_pp);
    AnalysisObject.SetUpLambdaPars_pL(DataSample,0,0,lam_pL);

    TFile* fOut = new TFile(OutputFolder+TString::Format("Output_%s_%s_%s.root",DataSample.Data(),SourceType.Data(),FittingMode_pp.Data()),"recreate");

    for(unsigned umT=0; umT<Num_mT_Bins; umT++) hCk[umT]->Write();

    CATS AB_pp;
    AB_pp.SetMaxNumThreads(4);
    AB_pp.SetMomBins(NumMomBins_pp,MomBins_pp);
    AnalysisObject.SetUpCats_pp(AB_pp,"AV18",SourceType);
    AB_pp.SetNotifications(CATS::nWarning);
    AB_pp.SetEpsilonConv(1e-9);
    AB_pp.SetEpsilonProp(1e-9);
    AB_pp.KillTheCat();
    DLM_Ck* Ck_pp = new DLM_Ck(AB_pp.GetNumSourcePars(),0,AB_pp);

    CATS AB_pL;
    AB_pL.SetMaxNumThreads(4);
    AB_pL.SetMomBins(NumMomBins_pp,MomBins_pp);
AnalysisObject.SetUpCats_pL(AB_pL,"Usmani","Gauss");
AB_pL.SetAnaSource(0,1.2);
    AB_pL.SetNotifications(CATS::nWarning);
    AB_pL.KillTheCat();
    DLM_Ck* Ck_pL = new DLM_Ck(AB_pL.GetNumSourcePars(),0,AB_pL);

    Ck_pp->Update();
    Ck_pL->Update();

    DLM_CkDecomposition CkDec_pp("pp",3,*Ck_pp,hResolution_pp);
    DLM_CkDecomposition CkDec_pL("pLambda",2,*Ck_pL,NULL);

    CkDec_pp.AddContribution(0,lam_pp[1],DLM_CkDecomposition::cFeedDown,&CkDec_pL,hResidual_pp_pL);
    CkDec_pp.AddContribution(1,lam_pp[2],DLM_CkDecomposition::cFeedDown);
    CkDec_pp.AddContribution(2,lam_pp[3],DLM_CkDecomposition::cFake);

    CkDec_pL.AddContribution(1,lam_pL[1]+lam_pL[2]+lam_pL[3],DLM_CkDecomposition::cFeedDown);
    CkDec_pL.AddContribution(2,lam_pL[4],DLM_CkDecomposition::cFake);

    CkDec_pp.Update();
    CkDec_pL.Update();

    TGraph* FitResult_pp = new TGraph[Num_mT_Bins];
    TGraphErrors Radius_pp;Radius_pp.SetName("Radius_pp");Radius_pp.Set(Num_mT_Bins);
    TGraphErrors Alpha_pp;Alpha_pp.SetName("Alpha_pp");Alpha_pp.Set(Num_mT_Bins);
    TGraph Chi2ndf_pp;Chi2ndf_pp.SetName("Chi2ndf_pp");Chi2ndf_pp.Set(Num_mT_Bins);

    unsigned NumSourcePars=0;

    double Value_mT;
    double Dummy;

    TGraph* gCkTh_pp = new TGraph[Num_mT_Bins];

    for(unsigned umT=0; umT<Num_mT_Bins; umT++){
        gCkTh_pp[umT].SetName(TString::Format("gCkTh_pp_%u",umT));
        gCkTh_pp[umT].Set(AB_pp.GetNumMomBins());
        DLM_Fitter1* fitter = new DLM_Fitter1(1);

        if(FittingMode_pp=="Norm"||FittingMode_pp.Contains("Baseline")){
            fitter->SetSystem(0,*hCk[umT],1,CkDec_pp,
                            FitRegion_pp[0],FitRegion_pp[1],FitRegion_pp[1],FitRegion_pp[1]);
        }
        else if(FittingMode_pp.Contains("Longbaseline")){
            fitter->SetSystem(0,*hCk[umT],1,CkDec_pp,
                            FitRegion_pp[0],FitRegion_pp[1],FitRegion_pp[2],FitRegion_pp[3]);
        }
        else{
            printf("\033[1;31mERROR:\033[0m Unknown fitting mode '%s'",FittingMode_pp.Data());
        }
        fitter->SetSeparateBL(0,false);

        fitter->SetParameter("pp",DLM_Fitter1::p_a,1.0,0.5,2.0);
        if(FittingMode_pp.Contains("Baseline")||FittingMode_pp.Contains("Longbaseline"))
            fitter->SetParameter("pp",DLM_Fitter1::p_b,0,-2e-3,2e-3);
        else
            fitter->FixParameter("pp",DLM_Fitter1::p_b,0);
        if(FittingMode_pp.Contains("Baseline2")||FittingMode_pp.Contains("Longbaseline2"))
            fitter->SetParameter("pp",DLM_Fitter1::p_c,0,-2e-4,2e-4);
        else
            fitter->FixParameter("pp",DLM_Fitter1::p_c,0);

        fitter->FixParameter("pp",DLM_Fitter1::p_Cl,-1);

        if(SourceType=="Gauss"){
            NumSourcePars = 1;
            fitter->SetParameter("pp",DLM_Fitter1::p_sor0,1.1,0.6,1.6);
        }
        else if(SourceType=="CleverLevy_Nolan"||SourceType=="CleverLevy_Single"||SourceType=="CleverLevy_Diff"){
            NumSourcePars = 2;
            fitter->SetParameter("pp",DLM_Fitter1::p_sor0,1.0,0.3,1.7);
            //fitter->SetParameter("pp",DLM_Fitter1::p_sor1,1.6,1.0,2.0);
            fitter->FixParameter("pp",DLM_Fitter1::p_sor1,2.0);
        }
        else if(SourceType=="McLevyNolan_Reso"){
            NumSourcePars = 2;
            fitter->SetParameter("pp",DLM_Fitter1::p_sor0,0.8,0.2,1.6);
            fitter->SetParameter("pp",DLM_Fitter1::p_sor1,1.6,1.0,2.0);
        }
        else if(SourceType=="McGauss_Reso"){
            NumSourcePars = 1;
            fitter->SetParameter("pp",DLM_Fitter1::p_sor0,0.8,0.2,1.6);
            //fitter->FixParameter("pp",DLM_Fitter1::p_sor1,2.0);
        }
        else if(SourceType=="GaussExpTotSimple_2body"){
            NumSourcePars = 11;
            //fitter->SetParameter("pp",DLM_Fitter1::p_sor0,0.8,0.2,1.6);
            fitter->FixParameter("pp",DLM_Fitter1::p_sor0,0.8);
        }
        else{
            printf("\033[1;31mERROR:\033[0m '%s' does not exist",SourceType.Data());
        }
        //fitter->AddSameSource("pLambda","pp",NumSourcePars);
        CkDec_pp.Update();
        CkDec_pL.Update();

        printf("Fitting mT bin %u\n",umT);
        fitter->GoBabyGo(false);

        FitResult_pp[umT].SetName(TString::Format("ppFit_mT_%u",umT));
        fitter->GetFitGraph(0, FitResult_pp[umT]);

        //printf("Fitted mT bin %u\n",umT);
        printf(" R = %.3f +/- %.3f\n",fitter->GetParameter("pp",DLM_Fitter1::p_sor0),fitter->GetParError("pp",DLM_Fitter1::p_sor0));

        Avg_mT->GetPoint(umT,Dummy,Value_mT);
        Radius_pp.SetPoint(umT,Value_mT,fitter->GetParameter("pp",DLM_Fitter1::p_sor0));
        Radius_pp.SetPointError(umT,0,fitter->GetParError("pp",DLM_Fitter1::p_sor0));
        if(NumSourcePars>1){
            printf(" α = %.3f +/- %.3f\n",fitter->GetParameter("pp",DLM_Fitter1::p_sor1),fitter->GetParError("pp",DLM_Fitter1::p_sor1));
            Alpha_pp.SetPoint(umT,Value_mT,fitter->GetParameter("pp",DLM_Fitter1::p_sor1));
            Alpha_pp.SetPointError(umT,0,fitter->GetParError("pp",DLM_Fitter1::p_sor1));
        }
        Chi2ndf_pp.SetPoint(umT,Value_mT,fitter->GetChi2Ndf());

        for(unsigned uMom=0; uMom<AB_pp.GetNumMomBins(); uMom++){
            gCkTh_pp[umT].SetPoint(uMom,AB_pp.GetMomentum(uMom),CkDec_pp.EvalCk(AB_pp.GetMomentum(uMom)));
        }

        fOut->cd();
        FitResult_pp[umT].Write();
        gCkTh_pp[umT].Write();

        delete fitter;
    }

    TGraph gCkTot_Th;
    gCkTot_Th.SetName("gCkTot_Th");
    gCkTot_Th.Set(AB_pp.GetNumMomBins());

    for(unsigned uMom=0; uMom<AB_pp.GetNumMomBins(); uMom++){
        double Momentum=AB_pp.GetMomentum(uMom);
        double CkVal=0;
        double Weight;

        //if(uMom==0) printf("Momentum[%u] = %f\n",uMom,Momentum);
        for(unsigned umT=0; umT<Num_mT_Bins; umT++){
            Avg_mT->GetPoint(umT,Dummy,Value_mT);
            Weight=hWeightsME->GetBinContent(hWeightsME->GetXaxis()->FindBin(Momentum*0.001),umT+1);//hWeightsME->GetYaxis()->FindBin(Value_mT)
            if(uMom==0) printf("Value_mT = %.3f, Bin_mT = %i\n",Value_mT,umT);
            CkVal += Weight*gCkTh_pp[umT].Eval(Momentum);
            if(uMom==0) printf(" CkVal[%u] += %.3f(%.3f)*%.3f\n",umT,Weight,Value_mT,gCkTh_pp[umT].Eval(Momentum));
        }
        //if(uMom==0) printf(" CkVal = %f\n",CkVal);
        gCkTot_Th.SetPoint(uMom,Momentum,CkVal);
    }

    TGraph gCkTot_Exp;
    gCkTot_Exp.SetName("gCkTot_Exp");
    gCkTot_Exp.Set(AB_pp.GetNumMomBins());

    TH1F* hCkTot_Exp = new TH1F("hCkTot_Exp","hCkTot_Exp",hTotCk->GetNbinsX(),hTotCk->GetBinLowEdge(1),hTotCk->GetXaxis()->GetBinUpEdge(hTotCk->GetNbinsX()));

    for(unsigned uMom=0; uMom<AB_pp.GetNumMomBins(); uMom++){
        double Momentum=AB_pp.GetMomentum(uMom);
        double CkVal=0;
        double Weight;
        for(unsigned umT=0; umT<Num_mT_Bins; umT++){
            Avg_mT->GetPoint(umT,Dummy,Value_mT);
            Weight=hWeightsME->GetBinContent(hWeightsME->GetXaxis()->FindBin(Momentum*0.001),umT+1);
            CkVal += Weight*FitResult_pp[umT].Eval(Momentum);
        }
        gCkTot_Exp.SetPoint(uMom,Momentum,CkVal);
    }

    for(unsigned uMom=0; uMom<hTotCk->GetNbinsX(); uMom++){
        double Momentum=hTotCk->GetBinCenter(uMom+1);
        double Weight;
        double CkVal_Exp=0;
        double CkErr_Exp=0;
        for(unsigned umT=0; umT<Num_mT_Bins; umT++){
            Avg_mT->GetPoint(umT,Dummy,Value_mT);
            Weight=hWeightsME->GetBinContent(hWeightsME->GetXaxis()->FindBin(Momentum*0.001),umT+1);
            CkVal_Exp += Weight*hCk[umT]->GetBinContent(uMom+1);
            CkErr_Exp += pow(Weight*hCk[umT]->GetBinError(uMom+1),2.);
        }
        hCkTot_Exp->SetBinContent(uMom+1,CkVal_Exp);
        hCkTot_Exp->SetBinError(uMom+1,sqrt(CkErr_Exp));
    }

    TF1* fit_CkTot = new TF1("fit_CkTot",[&](double*x, double *p){ return (p[0]+p[1]*x[0]+p[2]*x[0]*x[0])*gCkTot_Th.Eval(x[0]); }, AB_pp.GetMomBinLowEdge(0), AB_pp.GetMomBinUpEdge(AB_pp.GetNumMomBins()-1), 3);
    if(FittingMode_pp=="Norm"){
        fit_CkTot->FixParameter(1,0);
        fit_CkTot->FixParameter(2,0);
    }
    else if(FittingMode_pp=="Baseline"||FittingMode_pp=="Longbaseline"){
        fit_CkTot->FixParameter(2,0);
    }
    hTotCk->Fit(fit_CkTot,"S, N, R, M");

    double AVG_MT = 0;
    unsigned NUM_AVG_PTS = 0;
    for(double dMom=50; dMom<100; dMom+=4){
        AVG_MT+=hAvgKayTee->GetBinContent(hAvgKayTee->FindBin(dMom*0.001));
        NUM_AVG_PTS++;
    }
    AVG_MT /= double(NUM_AVG_PTS);
    printf("AVG_MT = %.3f\n",AVG_MT);
    double AVG_RAD = Radius_pp.Eval(AVG_MT);
    double AVG_ALFA = Alpha_pp.Eval(AVG_MT);
    printf(" AVG_RAD = %.3f\n",AVG_RAD);
    printf(" AVG_ALFA = %.3f\n",AVG_ALFA);

    DLM_Fitter1* fitter = new DLM_Fitter1(1);

    if(FittingMode_pp=="Norm"||FittingMode_pp.Contains("Baseline")){
        fitter->SetSystem(0,*hTotCk,1,CkDec_pp,
                        FitRegion_pp[0],FitRegion_pp[1],FitRegion_pp[1],FitRegion_pp[1]);
    }
    else if(FittingMode_pp.Contains("Longbaseline")){
        fitter->SetSystem(0,*hTotCk,1,CkDec_pp,
                        FitRegion_pp[0],FitRegion_pp[1],FitRegion_pp[2],FitRegion_pp[3]);
    }
    else{
        printf("\033[1;31mERROR:\033[0m Unknown fitting mode '%s'",FittingMode_pp.Data());
    }
    fitter->SetSeparateBL(0,false);

    fitter->SetParameter("pp",DLM_Fitter1::p_a,1.0,0.5,2.0);
    if(FittingMode_pp.Contains("Baseline")||FittingMode_pp.Contains("Longbaseline"))
        fitter->SetParameter("pp",DLM_Fitter1::p_b,0,-2e-3,2e-3);
    else
        fitter->FixParameter("pp",DLM_Fitter1::p_b,0);
    if(FittingMode_pp.Contains("Baseline2")||FittingMode_pp.Contains("Longbaseline2"))
        fitter->SetParameter("pp",DLM_Fitter1::p_c,0,-2e-4,2e-4);
    else
        fitter->FixParameter("pp",DLM_Fitter1::p_c,0);

    fitter->FixParameter("pp",DLM_Fitter1::p_Cl,-1);

    if(SourceType=="Gauss"){
        NumSourcePars = 1;
        fitter->FixParameter("pp",DLM_Fitter1::p_sor0,AVG_RAD);
    }
    else if(SourceType=="CleverLevy_Nolan"||SourceType=="CleverLevy_Single"||SourceType=="CleverLevy_Diff"){
        NumSourcePars = 2;
        fitter->FixParameter("pp",DLM_Fitter1::p_sor0,AVG_RAD);
        fitter->FixParameter("pp",DLM_Fitter1::p_sor1,AVG_ALFA);
    }
    else if(SourceType=="McLevyNolan_Reso"){
        NumSourcePars = 2;
        fitter->FixParameter("pp",DLM_Fitter1::p_sor0,AVG_RAD);
        fitter->FixParameter("pp",DLM_Fitter1::p_sor1,AVG_ALFA);
    }
    else if(SourceType=="McGauss_Reso"){
        NumSourcePars = 1;
        fitter->FixParameter("pp",DLM_Fitter1::p_sor0,AVG_RAD);
        //fitter->FixParameter("pp",DLM_Fitter1::p_sor1,2.0);
    }
    else if(SourceType=="GaussExpTotSimple_2body"){
        NumSourcePars = 11;
        fitter->FixParameter("pp",DLM_Fitter1::p_sor0,AVG_RAD);
    }
    else{
        printf("\033[1;31mERROR:\033[0m '%s' does not exist",SourceType.Data());
    }
    //fitter->AddSameSource("pLambda","pp",NumSourcePars);
    CkDec_pp.Update();
    CkDec_pL.Update();

    printf("Fitting TOTAL C(k)\n");
    fitter->GoBabyGo(false);

    TGraph ppFit_TOTAL;
    ppFit_TOTAL.SetName(TString::Format("ppFit_TOTAL"));
    fitter->GetFitGraph(0, ppFit_TOTAL);

    //printf("Fitted mT bin %u\n",umT);
    printf(" R = %.3f +/- %.3f\n",fitter->GetParameter("pp",DLM_Fitter1::p_sor0),fitter->GetParError("pp",DLM_Fitter1::p_sor0));

    //printf("AB_pp.EvalCorrFun(25)=%.4f\n",AB_pp.EvalCorrFun(25));
    //printf("AB_pp.EvaluateTheSource(25,2,0)=%.4f\n",AB_pp.EvaluateTheSource(25,2,0));
    //printf("pars: %.3f, %.3f\n",AB_pp.AnaSourcePar->GetParameter(0),AB_pp.AnaSourcePar->GetParameter(1));

    fOut->cd();
    ppFit_TOTAL.Write();

    delete fitter;

    //Radius_pp.GetXaxis()->SetRangeUser(0,3);
    TH1F* FUCK_ROOT = new TH1F("ROOT_IS_GREAT","ROOT_IS_GREAT",300,0,3);

    fOut->cd();
    Radius_pp.Write();
    Alpha_pp.Write();
    Chi2ndf_pp.Write();

    hTotCk->Write();
    hCkTot_Exp->Write();
    gCkTot_Th.Write();
    gCkTot_Exp.Write();
    fit_CkTot->Write();
    FUCK_ROOT->Write();

    delete [] hCk;
    delete hTotCk;
    delete hResolution_pp;
    delete hResidual_pp_pL;
    delete hCkTot_Exp;
    delete fInDiff;
    delete fInTot;
    delete [] MomBins_pp;
    delete [] FitRegion_pp;
    delete Ck_pp;
    delete Ck_pL;
    delete [] FitResult_pp;
    delete fit_CkTot;
    delete fOut;
}


//the Longbaseline2 does not work properly
void pL_FitDiff_CompareToTotal(char* argv[], const TString& DataSample, const TString& SourceType, const TString& FittingMode_pL){

    TString InputFileNameDiff = TString(argv[2])+"CFOutputALL_mT_pL_"+argv[3]+".root";
    TString InputFileNameTot = TString(argv[2])+"CFOutput_pL.root";
    const TString HistoNameWeights = "hWeights_0_ME";
    const TString HistoNameAvg_mT = "AveragemT";
    const TString HistoNameAvg_mT_k = "hAvgKayTee_0_SE";

    const TString BaseName = "hCk_RebinnedMeV_0";
    const TString MainHistoName = "hCk_ReweightedMeV_0";

    const TString OutputFolder = argv[0];

    TFile* fInDiff = new TFile(InputFileNameDiff,"read");
printf("fInDiff=%p\n",fInDiff);
//printf("InputFileNameDiff=%s fInDiff=%p\n",InputFileNameDiff.Data(),fInDiff);

    TH2F* hWeightsME = (TH2F*)fInDiff->Get(HistoNameWeights);
    TGraph* Avg_mT = (TGraph*)fInDiff->Get(HistoNameAvg_mT);
    TH1F* hAvgKayTee = (TH1F*)fInDiff->Get(HistoNameAvg_mT_k);
    const unsigned Num_mT_Bins = Avg_mT->GetN();
printf("hWeightsME=%p\n",hWeightsME);
printf("hAvgKayTee=%p\n",hAvgKayTee);
    TH1F** hCk = new TH1F* [Num_mT_Bins];
    printf("Num_mT_Bins=%u\n",Num_mT_Bins);
    for(unsigned umT=0; umT<Num_mT_Bins; umT++){
        hCk[umT] = (TH1F*)fInDiff->Get(BaseName+"_mTBin_"+umT);
printf("hCk[%u]=%p\n",umT,hCk[umT]);
        hCk[umT]->SetName(TString::Format("hCk_%u",umT));
    }

    TH1F* hTotCk;
    TFile* fInTot = new TFile(InputFileNameTot,"read");
    hTotCk = (TH1F*)fInTot->Get(MainHistoName);
    hTotCk->SetName("hTotCk");

    DLM_CommonAnaFunctions AnalysisObject;

    //double lambda_pp[4];
    //AnalysisObject.SetUpLambdaPars_pp(DataSample,0,lambda_pp);
    double lambda_pL[5];
    AnalysisObject.SetUpLambdaPars_pL(DataSample,0,0,lambda_pL);

    double* MomBins_pL = NULL;
    double* FitRegion_pL = NULL;
    unsigned NumMomBins_pL;

    AnalysisObject.SetUpBinning_pL(DataSample,NumMomBins_pL,MomBins_pL,FitRegion_pL,0,0);
    TH2F* hResolution_pL = AnalysisObject.GetResolutionMatrix(DataSample,"pLambda");

    double lam_pL[5];
    AnalysisObject.SetUpLambdaPars_pL(DataSample,0,0,lam_pL);

    TFile* fOut = new TFile(OutputFolder+TString::Format("Output_%s_%s_%s.root",DataSample.Data(),SourceType.Data(),FittingMode_pL.Data()),"recreate");

    for(unsigned umT=0; umT<Num_mT_Bins; umT++) hCk[umT]->Write();

    CATS AB_pL;
    AB_pL.SetMaxNumThreads(4);
    AB_pL.SetMomBins(NumMomBins_pL,MomBins_pL);
    AnalysisObject.SetUpCats_pL(AB_pL,"NLO_Coupled_S",SourceType);
    AB_pL.SetNotifications(CATS::nWarning);
    AB_pL.KillTheCat();
    DLM_Ck* Ck_pL = new DLM_Ck(AB_pL.GetNumSourcePars(),0,AB_pL);
    Ck_pL->Update();

    DLM_CkDecomposition CkDec_pL("pLambda",2,*Ck_pL,hResolution_pL);

    CkDec_pL.AddContribution(1,lam_pL[1]+lam_pL[2]+lam_pL[3],DLM_CkDecomposition::cFeedDown);
    CkDec_pL.AddContribution(2,lam_pL[4],DLM_CkDecomposition::cFake);

    CkDec_pL.Update();

    TGraph* FitResult_pL = new TGraph[Num_mT_Bins];
    TGraphErrors Radius_pL;Radius_pL.SetName("Radius_pL");Radius_pL.Set(Num_mT_Bins);
    TGraphErrors Alpha_pL;Alpha_pL.SetName("Alpha_pL");Alpha_pL.Set(Num_mT_Bins);
    TGraph Chi2ndf_pL;Chi2ndf_pL.SetName("Chi2ndf_pL");Chi2ndf_pL.Set(Num_mT_Bins);

    unsigned NumSourcePars=0;

    double Value_mT;
    double Dummy;

    TGraph* gCkTh_pL = new TGraph[Num_mT_Bins];

    for(unsigned umT=0; umT<Num_mT_Bins; umT++){
        gCkTh_pL[umT].SetName(TString::Format("gCkTh_pL_%u",umT));
        gCkTh_pL[umT].Set(AB_pL.GetNumMomBins());
        DLM_Fitter1* fitter = new DLM_Fitter1(1);

        if(FittingMode_pL=="Norm"||FittingMode_pL=="Baseline"){
            fitter->SetSystem(0,*hCk[umT],1,CkDec_pL,
                            FitRegion_pL[0],FitRegion_pL[1],FitRegion_pL[1],FitRegion_pL[1]);
        }
        else if(FittingMode_pL=="Longbaseline"){
            fitter->SetSystem(0,*hCk[umT],1,CkDec_pL,
                            FitRegion_pL[0],FitRegion_pL[1],FitRegion_pL[2],FitRegion_pL[3]);
        }
        else if(FittingMode_pL=="Longbaseline2"){
            fitter->SetSystem(0,*hCk[umT],1,CkDec_pL,
                            FitRegion_pL[0],FitRegion_pL[1],FitRegion_pL[2],FitRegion_pL[3]);
        }
        else{
            printf("\033[1;31mERROR:\033[0m Unknown fitting mode '%s'",FittingMode_pL.Data());
        }
        fitter->SetSeparateBL(0,false);

        fitter->SetParameter("pLambda",DLM_Fitter1::p_a,1.0,0.5,2.0);
        if(FittingMode_pL.Contains("Baseline")||FittingMode_pL.Contains("Longbaseline"))
            fitter->SetParameter("pLambda",DLM_Fitter1::p_b,0,-2e-3,2e-3);
        else
            fitter->FixParameter("pLambda",DLM_Fitter1::p_b,0);
        if(FittingMode_pL.Contains("Baseline2")||FittingMode_pL.Contains("Longbaseline2"))
            fitter->SetParameter("pLambda",DLM_Fitter1::p_c,0,-2e-3,2e-3);
        else
            fitter->FixParameter("pLambda",DLM_Fitter1::p_c,0);
        fitter->FixParameter("pLambda",DLM_Fitter1::p_Cl,-1);

        if(SourceType=="Gauss"){
            NumSourcePars = 1;
            fitter->SetParameter("pLambda",DLM_Fitter1::p_sor0,1.1,0.6,1.6);
        }
        else if(SourceType=="CleverLevy_Nolan"||SourceType=="CleverLevy_Single"||SourceType=="CleverLevy_Diff"){
            NumSourcePars = 2;
            fitter->SetParameter("pLambda",DLM_Fitter1::p_sor0,1.0,0.3,1.7);
            //fitter->SetParameter("pLambda",DLM_Fitter1::p_sor1,1.6,1.0,2.0);
            fitter->FixParameter("pLambda",DLM_Fitter1::p_sor1,2.0);
        }
        else if(SourceType=="McLevyNolan_Reso"){
            NumSourcePars = 2;
            fitter->SetParameter("pLambda",DLM_Fitter1::p_sor0,0.8,0.2,1.6);
            fitter->SetParameter("pLambda",DLM_Fitter1::p_sor1,1.6,1.0,2.0);
        }
        else if(SourceType=="McGauss_Reso"){
            NumSourcePars = 1;
            fitter->SetParameter("pLambda",DLM_Fitter1::p_sor0,0.8,0.2,1.6);
            //fitter->FixParameter("pLambda",DLM_Fitter1::p_sor1,2.0);
        }
        else if(SourceType=="GaussExpTotSimple_2body"){
            NumSourcePars = 11;
            //fitter->SetParameter("pLambda",DLM_Fitter1::p_sor0,0.8,0.2,1.6);
            fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,0.8);
        }
        else{
            printf("\033[1;31mERROR:\033[0m '%s' does not exist",SourceType.Data());
        }
        //fitter->AddSameSource("pLambda","pLambda",NumSourcePars);
        CkDec_pL.Update();

        printf("Fitting mT bin %u\n",umT);
        fitter->GoBabyGo(false);

        FitResult_pL[umT].SetName(TString::Format("pLambdaFit_mT_%u",umT));
        fitter->GetFitGraph(0, FitResult_pL[umT]);

        //printf("Fitted mT bin %u\n",umT);
        printf(" R = %.3f +/- %.3f\n",fitter->GetParameter("pLambda",DLM_Fitter1::p_sor0),fitter->GetParError("pLambda",DLM_Fitter1::p_sor0));

        Avg_mT->GetPoint(umT,Dummy,Value_mT);
        Radius_pL.SetPoint(umT,Value_mT,fitter->GetParameter("pLambda",DLM_Fitter1::p_sor0));
        Radius_pL.SetPointError(umT,0,fitter->GetParError("pLambda",DLM_Fitter1::p_sor0));
        if(NumSourcePars>1){
            printf(" α = %.3f +/- %.3f\n",fitter->GetParameter("pLambda",DLM_Fitter1::p_sor1),fitter->GetParError("pLambda",DLM_Fitter1::p_sor1));
            Alpha_pL.SetPoint(umT,Value_mT,fitter->GetParameter("pLambda",DLM_Fitter1::p_sor1));
            Alpha_pL.SetPointError(umT,0,fitter->GetParError("pLambda",DLM_Fitter1::p_sor1));
        }
        Chi2ndf_pL.SetPoint(umT,Value_mT,fitter->GetChi2Ndf());

        for(unsigned uMom=0; uMom<AB_pL.GetNumMomBins(); uMom++){
            gCkTh_pL[umT].SetPoint(uMom,AB_pL.GetMomentum(uMom),CkDec_pL.EvalCk(AB_pL.GetMomentum(uMom)));
        }

        fOut->cd();
        FitResult_pL[umT].Write();
        gCkTh_pL[umT].Write();

        delete fitter;
    }

    TGraph gCkTot_Th;
    gCkTot_Th.SetName("gCkTot_Th");
    gCkTot_Th.Set(AB_pL.GetNumMomBins());

    for(unsigned uMom=0; uMom<AB_pL.GetNumMomBins(); uMom++){
        double Momentum=AB_pL.GetMomentum(uMom);
        double CkVal=0;
        double Weight;

        //if(uMom==0) printf("Momentum[%u] = %f\n",uMom,Momentum);
        for(unsigned umT=0; umT<Num_mT_Bins; umT++){
            Avg_mT->GetPoint(umT,Dummy,Value_mT);
            Weight=hWeightsME->GetBinContent(hWeightsME->GetXaxis()->FindBin(Momentum*0.001),umT+1);
            CkVal += Weight*gCkTh_pL[umT].Eval(Momentum);
            //if(uMom==0) printf(" CkVal[%u] += %f*%f\n",umT,Weight,gCkTh_pL[umT].Eval(Momentum));
        }
        //if(uMom==0) printf(" CkVal = %f\n",CkVal);
        gCkTot_Th.SetPoint(uMom,Momentum,CkVal);
    }

    TGraph gCkTot_Exp;
    gCkTot_Exp.SetName("gCkTot_Exp");
    gCkTot_Exp.Set(AB_pL.GetNumMomBins());

    for(unsigned uMom=0; uMom<AB_pL.GetNumMomBins(); uMom++){
        double Momentum=AB_pL.GetMomentum(uMom);
        double CkVal=0;
        double Weight;
        for(unsigned umT=0; umT<Num_mT_Bins; umT++){
            Avg_mT->GetPoint(umT,Dummy,Value_mT);
            Weight=hWeightsME->GetBinContent(hWeightsME->GetXaxis()->FindBin(Momentum*0.001),umT+1);
            CkVal += Weight*FitResult_pL[umT].Eval(Momentum);
        }
        gCkTot_Exp.SetPoint(uMom,Momentum,CkVal);
    }

    TF1* fit_CkTot = new TF1("fit_CkTot",[&](double*x, double *p){ return (p[0]+p[1]*x[0]+p[2]*x[0]*x[0])*gCkTot_Th.Eval(x[0]); }, AB_pL.GetMomBinLowEdge(0), AB_pL.GetMomBinUpEdge(AB_pL.GetNumMomBins()-1), 3);
    if(FittingMode_pL=="Norm"){
        fit_CkTot->FixParameter(1,0);
        fit_CkTot->FixParameter(2,0);
    }
    else if(FittingMode_pL=="Baseline"||FittingMode_pL=="Longbaseline"){
        fit_CkTot->FixParameter(2,0);
    }
    else{

    }
    hTotCk->Fit(fit_CkTot,"S, N, R, M");

    double AVG_MT = 0;
    unsigned NUM_AVG_PTS = 0;
    for(double dMom=50; dMom<100; dMom+=4){
        AVG_MT+=hAvgKayTee->GetBinContent(hAvgKayTee->FindBin(dMom*0.001));
        NUM_AVG_PTS++;
    }
    AVG_MT /= double(NUM_AVG_PTS);
    printf("AVG_MT = %.3f\n",AVG_MT);
    double AVG_RAD = Radius_pL.Eval(AVG_MT);
    double AVG_ALFA = Alpha_pL.Eval(AVG_MT);
    printf(" AVG_RAD = %.3f\n",AVG_RAD);
    printf(" AVG_ALFA = %.3f\n",AVG_ALFA);

    DLM_Fitter1* fitter = new DLM_Fitter1(1);

    if(FittingMode_pL=="Norm"||FittingMode_pL=="Baseline"){
        fitter->SetSystem(0,*hTotCk,1,CkDec_pL,
                        FitRegion_pL[0],FitRegion_pL[1],FitRegion_pL[1],FitRegion_pL[1]);
    }
    else if(FittingMode_pL=="Longbaseline"||FittingMode_pL=="Longbaseline2"){
        fitter->SetSystem(0,*hTotCk,1,CkDec_pL,
                        FitRegion_pL[0],FitRegion_pL[1],FitRegion_pL[2],FitRegion_pL[3]);
    }
    else{
        printf("\033[1;31mERROR:\033[0m Unknown fitting mode '%s'",FittingMode_pL.Data());
    }
    fitter->SetSeparateBL(0,false);

    fitter->SetParameter("pLambda",DLM_Fitter1::p_a,1.0,0.5,2.0);
    if(FittingMode_pL.Contains("Baseline")||FittingMode_pL.Contains("Longbaseline"))
        fitter->SetParameter("pLambda",DLM_Fitter1::p_b,0,-2e-3,2e-3);
    else
        fitter->FixParameter("pLambda",DLM_Fitter1::p_b,0);
    if(FittingMode_pL.Contains("Baseline2")||FittingMode_pL.Contains("Longbaseline2"))
        fitter->SetParameter("pLambda",DLM_Fitter1::p_c,0,-2e-4,2e-4);
    else
        fitter->FixParameter("pLambda",DLM_Fitter1::p_c,0);
    fitter->FixParameter("pLambda",DLM_Fitter1::p_Cl,-1);

    if(SourceType=="Gauss"){
        NumSourcePars = 1;
        fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,AVG_RAD);
    }
    else if(SourceType=="CleverLevy_Nolan"||SourceType=="CleverLevy_Single"||SourceType=="CleverLevy_Diff"){
        NumSourcePars = 2;
        fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,AVG_RAD);
        fitter->FixParameter("pLambda",DLM_Fitter1::p_sor1,AVG_ALFA);
    }
    else if(SourceType=="McLevyNolan_Reso"){
        NumSourcePars = 2;
        fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,AVG_RAD);
        fitter->FixParameter("pLambda",DLM_Fitter1::p_sor1,AVG_ALFA);
    }
    else if(SourceType=="McGauss_Reso"){
        NumSourcePars = 1;
        fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,AVG_RAD);
        //fitter->FixParameter("pLambda",DLM_Fitter1::p_sor1,2.0);
    }
    else if(SourceType=="GaussExpTotSimple_2body"){
        NumSourcePars = 11;
        fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,AVG_RAD);
    }
    else{
        printf("\033[1;31mERROR:\033[0m '%s' does not exist",SourceType.Data());
    }
    //fitter->AddSameSource("pLambda","pLambda",NumSourcePars);
    CkDec_pL.Update();

    printf("Fitting TOTAL C(k)\n");
    fitter->GoBabyGo(false);

    TGraph pLambdaFit_TOTAL;
    pLambdaFit_TOTAL.SetName(TString::Format("pLambdaFit_TOTAL"));
    fitter->GetFitGraph(0, pLambdaFit_TOTAL);

    //printf("Fitted mT bin %u\n",umT);
    printf(" R = %.3f +/- %.3f\n",fitter->GetParameter("pLambda",DLM_Fitter1::p_sor0),fitter->GetParError("pLambda",DLM_Fitter1::p_sor0));

    //printf("AB_pL.EvalCorrFun(25)=%.4f\n",AB_pL.EvalCorrFun(25));
    //printf("AB_pL.EvaluateTheSource(25,2,0)=%.4f\n",AB_pL.EvaluateTheSource(25,2,0));
    //printf("pars: %.3f, %.3f\n",AB_pL.AnaSourcePar->GetParameter(0),AB_pL.AnaSourcePar->GetParameter(1));

    fOut->cd();
    pLambdaFit_TOTAL.Write();

    delete fitter;

    TH1F* FUCK_ROOT = new TH1F("ROOT_IS_GREAT","ROOT_IS_GREAT",300,0.5,2.5);
    //Radius_pL.GetXaxis()->SetRangeUser(0.5,2.5);

    fOut->cd();
    Radius_pL.Write();
    Alpha_pL.Write();
    Chi2ndf_pL.Write();

    hTotCk->Write();
    gCkTot_Th.Write();
    gCkTot_Exp.Write();
    fit_CkTot->Write();
    FUCK_ROOT->Write();

    delete [] hCk;
    delete hTotCk;
    delete hResolution_pL;
    delete fInDiff;
    delete fInTot;
    delete [] MomBins_pL;
    delete [] FitRegion_pL;
    delete Ck_pL;
    delete [] FitResult_pL;
    delete fit_CkTot;
    delete fOut;
}




void mT_Scaled_Source_Prelim2019(){

    //file with the mT scaling
    const TString InputFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pp_FitDiff_CompareToTotal/250319/HM/Output_pp13TeV_HM_March19_McGauss_Reso_Norm.root";
    const TString OutputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/mT_Scaled_Source_Prelim2019/";
    const TString mT_GraphName = "Radius_pp";

    const unsigned NumPartPairs = 5;
    TString* PairName = new TString [NumPartPairs];
    enum PrevNext { pp, pL, pS0, pXi, pOmega };
    PairName[0] = "pp";
    PairName[1] = "pL";
    PairName[2] = "pS0";
    PairName[3] = "pXim";
    PairName[4] = "pOmega";

    TFile* fInput = new TFile(InputFileName,"read");
    TGraph* Radius_pp = (TGraph*)fInput->Get(mT_GraphName);

    double* Pair_Avg_mT = new double[NumPartPairs];
    double* Part0_Mass = new double[NumPartPairs];
    double* Part0_MassReso = new double[NumPartPairs];
    double* Part0_TauReso = new double[NumPartPairs];
    double* Part0_FracPrim = new double[NumPartPairs];
    double* Part1_Mass = new double[NumPartPairs];
    double* Part1_MassReso = new double[NumPartPairs];
    double* Part1_TauReso = new double[NumPartPairs];
    double* Part1_FracPrim = new double[NumPartPairs];

//VALUE FOR pp HM 13 TEV
    Pair_Avg_mT[pp] = 1.35;
    Part0_Mass[pp] = 938.272;
    Part0_MassReso[pp] = 1361.52;
    Part0_TauReso[pp] = 1.65;
    Part0_FracPrim[pp] = 0.3578;
    Part1_Mass[pp] = 938.272;
    Part1_MassReso[pp] = 1361.52;
    Part1_TauReso[pp] = 1.65;
    Part1_FracPrim[pp] = 0.3578;

    //<mT> p-Xi = 1.84 Units of mT xD

    Pair_Avg_mT[pL] = 1.55;
    Part0_Mass[pL] = 938.272;
    Part0_MassReso[pL] = 1361.52;
    Part0_TauReso[pL] = 1.65;
    Part0_FracPrim[pL] = 0.3578;
    Part1_Mass[pL] = 1115.683;
    Part1_MassReso[pL] = 1462.93;
    Part1_TauReso[pL] = 4.69;
    Part1_FracPrim[pL] = 0.3562;

    Pair_Avg_mT[pS0] = 2.07;
    Part0_Mass[pS0] = 938.272;
    Part0_MassReso[pS0] = 1361.52;
    Part0_TauReso[pS0] = 1.65;
    Part0_FracPrim[pS0] = 0.3578;
    Part1_Mass[pS0] = 1192.642;
    Part1_MassReso[pS0] = 1581.73;
    Part1_TauReso[pS0] = 4.28;
    Part1_FracPrim[pS0] = 0.3735;

    Pair_Avg_mT[pXi] = 1.85;
    Part0_Mass[pXi] = 938.272;
    Part0_MassReso[pXi] = 1361.52;
    Part0_TauReso[pXi] = 1.65;
    Part0_FracPrim[pXi] = 0.3578;
    Part1_Mass[pXi] = 1321.7;
    Part1_MassReso[pXi] = 1361.52;
    Part1_TauReso[pXi] = 0;
    Part1_FracPrim[pXi] = 1;

    Pair_Avg_mT[pOmega] = 2.17;
    Part0_Mass[pOmega] = 938.272;
    Part0_MassReso[pOmega] = 1361.52;
    Part0_TauReso[pOmega] = 1.65;
    Part0_FracPrim[pOmega] = 0.3578;
    Part1_Mass[pOmega] = 1672.45;
    Part1_MassReso[pOmega] = 1672.45;
    Part1_TauReso[pOmega] = 0;
    Part1_FracPrim[pOmega] = 1;

    TString* PairRootName = new TString [NumPartPairs];
    PairRootName[pp] = TString::Format("p#minusp <m_{T}>=%.2f MeV",Pair_Avg_mT[pp]);
    PairRootName[pL] = TString::Format("p#minus#Lambda <m_{T}>=%.2f MeV",Pair_Avg_mT[pL]);
    PairRootName[pS0] = TString::Format("p#minus#Sigma^{0} <m_{T}>=%.2f MeV",Pair_Avg_mT[pS0]);
    PairRootName[pXi] = TString::Format("p#minus#Xi^{#minus} <m_{T}>=%.2f MeV",Pair_Avg_mT[pXi]);
    PairRootName[pOmega] = TString::Format("p#minus#Omega^{#minus} <m_{T}>=%.2f MeV",Pair_Avg_mT[pOmega]);

    double* SourceScale = new double[NumPartPairs];
    for(unsigned uPair=0; uPair<NumPartPairs; uPair++){
        SourceScale[uPair] = Radius_pp->Eval(Pair_Avg_mT[uPair]);
    }

    TGraph* gRad = new TGraph[NumPartPairs];
    const double RadMin=0;
    const double RadMax=32;
    const unsigned NumRadSteps = 256;

    DLM_CleverMcLevyReso* CleverMcLevyReso = new DLM_CleverMcLevyReso [NumPartPairs];
    //#pragma omp parallel for
    for(unsigned uPair=0; uPair<NumPartPairs; uPair++){

        const double RadBinWidth = (RadMax-RadMin)/(NumRadSteps);
        CleverMcLevyReso[uPair].InitStability(1,2-1e-6,2+1e-6);
        CleverMcLevyReso[uPair].InitScale(35,0.25,2.0);
        CleverMcLevyReso[uPair].InitRad(NumRadSteps,RadMin,RadMax);
        CleverMcLevyReso[uPair].InitType(2);
        CleverMcLevyReso[uPair].InitReso(0,1);
        CleverMcLevyReso[uPair].InitReso(1,1);
        CleverMcLevyReso[uPair].SetUpReso(0,0,1.-Part0_FracPrim[uPair],Part0_MassReso[uPair],Part0_TauReso[uPair],Part0_Mass[uPair],139.57);
        CleverMcLevyReso[uPair].SetUpReso(1,0,1.-Part1_FracPrim[uPair],Part1_MassReso[uPair],Part1_TauReso[uPair],Part1_Mass[uPair],139.57);
        CleverMcLevyReso[uPair].InitNumMcIter(2000000);

        #pragma critical
        {
        gRad[uPair].SetName("hRad_"+PairName[uPair]);
        gRad[uPair].Set(NumRadSteps);
        }

        double EvalValue;
        unsigned Counter=0;
        double PARS[5];
        PARS[3]=SourceScale[uPair];
        PARS[4]=2;
        for(double dRad=RadBinWidth*0.5; dRad<RadMax; dRad+=RadBinWidth){
            PARS[1]=dRad;
            EvalValue = CleverMcLevyReso[uPair].Eval(PARS);
            //EvalValue = CleverMcLevyReso[uPair].Eval(0,dRad,0);
            #pragma critical
            {
            gRad[uPair].SetPoint(Counter,dRad,EvalValue);
            }
            Counter++;
        }
    }

    TFile* fOutput = new TFile(OutputFolder+"fOutput.root","recreate");
    for(unsigned uPair=0; uPair<NumPartPairs; uPair++){
        gRad[uPair].Write();
    }

    TH1F* hAxisRad = new TH1F("hAxisRad", "hAxisRad", NumRadSteps,RadMin,RadMax);
    hAxisRad->SetStats(false);
    hAxisRad->SetTitle("");
    hAxisRad->GetXaxis()->SetLabelSize(0.065);
    hAxisRad->GetXaxis()->SetTitle("r (fm)");
    hAxisRad->GetXaxis()->CenterTitle();
    hAxisRad->GetXaxis()->SetTitleOffset(1.2);
    hAxisRad->GetXaxis()->SetLabelOffset(0.02);
    hAxisRad->GetXaxis()->SetTitleSize(0.075);
    hAxisRad->GetYaxis()->SetLabelSize(0.065);
    hAxisRad->GetYaxis()->SetTitle("4#pir^{2}S(r) (fm^{-1})");
    hAxisRad->GetYaxis()->CenterTitle();
    hAxisRad->GetYaxis()->SetTitleOffset(1.00);
    hAxisRad->GetYaxis()->SetTitleSize(0.075);

    hAxisRad->GetYaxis()->SetRangeUser(0, 0.5);
    //hAxisRad->GetXaxis()->SetNdivisions(506);
    hAxisRad->GetXaxis()->SetRangeUser(0, 8);

    int COLORS[NumPartPairs] = {kBlue,kRed+1,kRed+3,kCyan+2,kPink+9};
    TCanvas* cSource = new TCanvas("cSource", "cSource", 1);
    cSource->cd(0); cSource->SetCanvasSize(1920, 1080); cSource->SetMargin(0.15,0.05,0.2,0.05);//lrbt
    hAxisRad->Draw();

    TPaveText* PT1 = new TPaveText(0.4,0.9,0.975,0.975, "blNDC");//lbrt
    PT1->SetName("PT1");
    PT1->SetBorderSize(1);
    PT1->SetTextSize(0.04);
    PT1->SetFillColor(kWhite);
    PT1->SetTextFont(22);
    PT1->AddText("Source using Gaussian core, based on <m_{T}>, plus resonances");

    TLegend* lInfo1 = new TLegend(0.52,0.6,0.975,0.9);//lbrt
    lInfo1->SetName(TString::Format("lInfo1"));
    lInfo1->SetTextSize(0.04);
    for(unsigned uPair=0; uPair<NumPartPairs; uPair++){
        gRad[uPair].SetLineColor(COLORS[uPair]);
        gRad[uPair].SetMarkerColor(COLORS[uPair]);
        gRad[uPair].SetMarkerSize(0);
        gRad[uPair].SetLineWidth(6);
        lInfo1->AddEntry(&gRad[uPair],TString::Format("%s (R_{G} = %.2f fm)",PairRootName[uPair].Data(),SourceScale[uPair]));
        gRad[uPair].Draw("C,same");
    }
    lInfo1->Draw("same");
    PT1->Draw("same");

    cSource->Write();
    cSource->SaveAs(OutputFolder+"cSource.png");


    delete [] PairName;
    delete [] PairRootName;
    delete [] Pair_Avg_mT;
    delete [] Part0_Mass;
    delete [] Part0_MassReso;
    delete [] Part0_TauReso;
    delete [] Part0_FracPrim;
    delete [] Part1_Mass;
    delete [] Part1_MassReso;
    delete [] Part1_TauReso;
    delete [] Part1_FracPrim;
    delete [] gRad;
    delete hAxisRad;
    delete PT1;
    delete lInfo1;
    delete cSource;
    delete fInput;
    delete fOutput;
    delete [] CleverMcLevyReso;
    delete [] SourceScale;
}

void Plot_pL_Vars(char* argv[], const unsigned NumVar, double* Dev_pL){
    TGraph gVars_pL;
    gVars_pL.SetName("gVars_pL");
    for(unsigned uVar=1; uVar<NumVar; uVar++){
        gVars_pL.SetPoint(uVar-1,uVar,Dev_pL[uVar]*100.);
    }
    gVars_pL.SetMarkerColor(kBlack);
    gVars_pL.SetMarkerSize(1);
    gVars_pL.SetMarkerStyle(20);


    TH1F* hAxis = new TH1F("hAxis", "hAxis", 128, 0, NumVar);
    hAxis->SetStats(false);
    hAxis->SetTitle("");
    hAxis->GetXaxis()->SetTitle("Variation");
    hAxis->GetXaxis()->SetTitleSize(0.06);
    hAxis->GetXaxis()->SetLabelSize(0.06);
    hAxis->GetXaxis()->CenterTitle();
    hAxis->GetXaxis()->SetTitleOffset(1.3);
    hAxis->GetXaxis()->SetLabelOffset(0.02);

    hAxis->GetYaxis()->SetTitle("Rel. variation (%)");
    hAxis->GetYaxis()->SetTitleSize(0.06);
    hAxis->GetYaxis()->SetLabelSize(0.06);
    hAxis->GetYaxis()->CenterTitle();
    hAxis->GetYaxis()->SetTitleOffset(1.1);
    //hAxis->GetXaxis()->SetNdivisions(506);
    hAxis->GetYaxis()->SetRangeUser(-25,25);
    //hAxis->GetYaxis()->SetRangeUser(0.90, 1.80);

    TCanvas* cVar = new TCanvas("cVar", "cVar", 1);
    cVar->cd(0); cVar->SetCanvasSize(1920/2, 1080/2); cVar->SetMargin(0.15,0.05,0.2,0.05);//lrbt
    hAxis->Draw("axis");
    gVars_pL.Draw("same,P");

    TString FileName = argv[0];
    FileName+="cVar_pL.png";
    cVar->SaveAs(FileName.Data());

    delete hAxis;
    delete cVar;
}

void DimiExecuteCFmT(){
    char** ARGV = new char* [8];
    for(int i=0; i<8; i++) ARGV[i]=new char [512];

//mT_Scaled_Source_Prelim2019();

//home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample10HM/Systematics
    //strcpy(ARGV[0],"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample13HM/");
    //strcpy(ARGV[1],"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample13HM/AnalysisResults.root");
    //strcpy(ARGV[2],"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample13HM/");
    //strcpy(ARGV[3],"PL");
    //strcpy(ARGV[4],"");

    strcpy(ARGV[0],"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample10HM_ver2/");
    strcpy(ARGV[1],"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample10HM_ver2/AnalysisResults.root");
    strcpy(ARGV[2],"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample10HM_ver2/");
    strcpy(ARGV[3],"HM");
    strcpy(ARGV[4],"");

    //strcpy(ARGV[0],"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample9/");
    //strcpy(ARGV[1],"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample9/AnalysisResults.root");
    //strcpy(ARGV[2],"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample9/");
    //strcpy(ARGV[3],"MB");
    //strcpy(ARGV[4],"");

    //strcpy(ARGV[0],"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pPb_5TeV/Sample10/");
    //strcpy(ARGV[1],"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pPb_5TeV/Sample10/AnalysisResults.root");
    //strcpy(ARGV[2],"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pPb_5TeV/Sample10/");
    //strcpy(ARGV[3],"MB");
    //strcpy(ARGV[4],"");

    //strcpy(ARGV[0],"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pPb_5TeV/Sample11/");
    //strcpy(ARGV[1],"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pPb_5TeV/Sample11/AnalysisResults.root");
    //strcpy(ARGV[2],"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pPb_5TeV/Sample11/");
    //strcpy(ARGV[3],"MB");
    //strcpy(ARGV[4],"");


    gSystem->cd(ARGV[0]);


    const unsigned NumVar=45;//45
    double* Dev_pL = new double[NumVar];
    //ExecuteCFDream(5,ARGV,Dev_pL[0]);
    for(int iVar=0; iVar<NumVar; iVar++){
        char buffer [8];
        //sprintf(buffer, "%i", iVar);
        //strcpy(ARGV[4],buffer);
        //ExecuteCFDream(5,ARGV,Dev_pL[iVar]);
        //ExecuteCFmT_pL(5,ARGV);
    }
    //Plot_pL_Vars(ARGV,NumVar,Dev_pL);
    delete [] Dev_pL;

    ExecuteCFmT_pp(5,ARGV);
    //QA_for_pp_mT();
    //ExecuteCFmT_pL(5,ARGV);
return;
    //strcpy(ARGV[0],"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pp_FitDiff_CompareToTotal/250319/HM/pL/NLO/");
    //pp_FitDiff_CompareToTotal(ARGV,"pp13TeV_HM_March19","Gauss","Norm");
    //pp_FitDiff_CompareToTotal(ARGV,"pp13TeV_HM_March19","McGauss_Reso","Norm");
    //pp_FitDiff_CompareToTotal(ARGV,"pp13TeV_HM_March19","McLevyNolan_Reso","Norm");

    //pL_FitDiff_CompareToTotal(ARGV,"pp13TeV_HM_March19","Gauss","Norm");
    //pL_FitDiff_CompareToTotal(ARGV,"pp13TeV_HM_March19","McGauss_Reso","Norm");

    //pL_FitDiff_CompareToTotal(ARGV,"pp13TeV_HM_March19","Gauss","Longbaseline2");
    //pL_FitDiff_CompareToTotal(ARGV,"pp13TeV_HM_March19","McGauss_Reso","Longbaseline2");

    strcpy(ARGV[0],"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pp_FitDiff_CompareToTotal/Test/");
    //pp_FitDiff_CompareToTotal(ARGV,"pp13TeV_HM_March19","Gauss","Baseline");
    //pp_FitDiff_CompareToTotal(ARGV,"pp13TeV_HM_March19","McGauss_Reso","Norm");
    //pp_FitDiff_CompareToTotal(ARGV,"pp13TeV_HM_March19","McGauss_Reso","Baseline");
    //pp_FitDiff_CompareToTotal(ARGV,"pp13TeV_HM_March19","McGauss_Reso","Baseline2");
    //pp_FitDiff_CompareToTotal(ARGV,"pp13TeV_HM_March19","McLevyNolan_Reso","Norm");
    //pp_FitDiff_CompareToTotal(ARGV,"pp13TeV_HM_March19","McLevyNolan_Reso","Baseline");
    //pp_FitDiff_CompareToTotal(ARGV,"pp13TeV_HM_March19","McLevyNolan_Reso","Baseline2");

    //strcpy(ARGV[0],"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pp_FitDiff_CompareToTotal/");
    //pp_FitDiff_CompareToTotal(ARGV,"pPb5TeV_Run2paper","Gauss","Norm");
    //pp_FitDiff_CompareToTotal(ARGV,"pPb5TeV_Run2paper","McGauss_Reso","Norm");
    //pp_FitDiff_CompareToTotal("ARGV,pPb5TeV_Run2paper","McLevyNolan_Reso","Norm");

//return;
    //strcpy(ARGV[0],"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pp_FitDiff_CompareToTotal/270319/pPb/");
    //pp_FitDiff_CompareToTotal(ARGV,"pp13TeV_MB_Run2paper","Gauss","Norm");
    //pp_FitDiff_CompareToTotal(ARGV,"pp13TeV_MB_Run2paper","GaussExpTotSimple_2body","Norm");
    //pp_FitDiff_CompareToTotal(ARGV,"pp13TeV_MB_Run2paper","McGauss_Reso","Norm");
    //pp_FitDiff_CompareToTotal(ARGV,"pp13TeV_MB_Run2paper","McLevyNolan_Reso","Norm");

    //pp_FitDiff_CompareToTotal(ARGV,"pp13TeV_MB_Run2paper","Gauss","Baseline");
    //pp_FitDiff_CompareToTotal(ARGV,"pp13TeV_MB_Run2paper","GaussExpTotSimple_2body","Norm");
    //pp_FitDiff_CompareToTotal(ARGV,"pp13TeV_MB_Run2paper","McGauss_Reso","Baseline");
    //pp_FitDiff_CompareToTotal(ARGV,"pp13TeV_MB_Run2paper","McLevyNolan_Reso","Baseline");

    //pp_FitDiff_CompareToTotal(ARGV,"pPb5TeV_Run2paper","Gauss","Norm");
    //pp_FitDiff_CompareToTotal(ARGV,"pPb5TeV_Run2paper","GaussExpTotSimple_2body","Norm");
    //pp_FitDiff_CompareToTotal(ARGV,"pPb5TeV_Run2paper","McGauss_Reso","Norm");
    //pp_FitDiff_CompareToTotal(ARGV,"pPb5TeV_Run2paper","McLevyNolan_Reso","Norm");

    //strcpy(ARGV[0],"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pp_FitDiff_CompareToTotal/250319/pPb/");
    //pp_FitDiff_CompareToTotal(ARGV,"pPb5TeV_CPR_Mar19","Gauss","Norm");
    //pp_FitDiff_CompareToTotal(ARGV,"pPb5TeV_CPR_Mar19","McGauss_Reso","Norm");
    //pp_FitDiff_CompareToTotal(ARGV,"pPb5TeV_CPR_Mar19","McLevyNolan_Reso","Norm");
    //pp_FitDiff_CompareToTotal(ARGV,"pPb5TeV_CPR_Mar19","Gauss","Baseline");
    //pp_FitDiff_CompareToTotal(ARGV,"pPb5TeV_CPR_Mar19","McGauss_Reso","Baseline");
    //pp_FitDiff_CompareToTotal(ARGV,"pPb5TeV_CPR_Mar19","McLevyNolan_Reso","Baseline");

    //strcpy(ARGV[0],"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pp_FitDiff_CompareToTotal/250319/pPb/pL/LO/");
    //pL_FitDiff_CompareToTotal(ARGV,"pPb5TeV_CPR_Mar19","McGauss_Reso","Norm");

    for(int i=0; i<8; i++) delete[]ARGV[i];
    delete [] ARGV;
}
