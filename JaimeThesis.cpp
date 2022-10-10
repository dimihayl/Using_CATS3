
#include "JaimeThesis.h"
#include "EnvVars.h"
#include "FemtoBoyzScripts.h"
#include "DLM_SubPads.h"

#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TLegend.h"


//plots the correlation functions (based on data) and the final fit that
//Jaime gave me
//mode == 0; default
void Plot_CorrFunFits(const int& mode, TString FitSettings){

  const unsigned NumDataVar = 43;
  const unsigned NumMtBins = 6;
  double kMin = 0;
  double kMax = 432;
  unsigned NumMomBins = 36;

  TString InDirName = TString::Format("%s/Jaime/Thesis/Input/",GetCernBoxDimi());
  TString InFileName = "ck-fit_"+FitSettings+".root";
  TString OutDirName = TString::Format("%s/Jaime/Thesis/Output/",GetCernBoxDimi());
  TString OutFileName = "ck-plot_"+FitSettings+".root";


  TGraphErrors* gFit = new TGraphErrors[NumMtBins];
  TGraphErrors* gBl = new TGraphErrors[NumMtBins];
  TGraphErrors* gNsig = new TGraphErrors[NumMtBins];
  TH1F** hData = new TH1F*[NumMtBins];
  TH1F** hDataInlet = new TH1F*[NumMtBins];

  TRandom3 rangen(11);

  TFile InFile(InDirName+InFileName,"read");
  for(unsigned uMt=0; uMt<NumMtBins; uMt++){
    gFit[uMt].SetName(TString::Format("gFit_%u",uMt));
    gBl[uMt].SetName(TString::Format("gBl_%u",uMt));
    gNsig[uMt].SetName(TString::Format("gNsig_%u",uMt));
    //gData[uMt].SetName(TString::Format("gData_%u",uMt));
    TH2F* gFitMt = (TH2F*)InFile.Get(TString::Format("gFit mT%u",uMt));
    TH2F* gBlMt = (TH2F*)InFile.Get(TString::Format("gBl mT%u",uMt));
    hData[uMt] = new TH1F(TString::Format("hData_%u",uMt),TString::Format("hData_%u",uMt),NumMomBins,kMin,kMax);
    for(unsigned uMom=0; uMom<NumMomBins; uMom++){
      TH1F* hProj = (TH1F*)gFitMt->ProjectionY(TString::Format("hProj"),uMom+1,uMom+1);
      gFit[uMt].SetPoint(uMom,gFitMt->GetXaxis()->GetBinCenter(uMom+1),hProj->GetMean());
      gFit[uMt].SetPointError(uMom,0,hProj->GetStdDev());
      delete hProj;

      hProj = (TH1F*)gBlMt->ProjectionY(TString::Format("hProj"),uMom+1,uMom+1);
      gBl[uMt].SetPoint(uMom,gBlMt->GetXaxis()->GetBinCenter(uMom+1),hProj->GetMean());
      gBl[uMt].SetPointError(uMom,0,hProj->GetStdDev());
      delete hProj;
    }

    TH2F* hData2D = new TH2F("hData2D","hData2D",NumMomBins,kMin,kMax,2048,0,4);
    TH1F* hDataDef = NULL;
    for(unsigned uVar=0; uVar<NumDataVar; uVar++){
      TFile DataFile(TString::Format("%s/CatsFiles/ExpData/Bernie_Source/pLData/mTBin_%u/CFOutput_mT_pLVar%u_HM_%u.root",GetCernBoxDimi(),uMt+1,uVar,uMt),"read");
      TString HistoName = TString::Format("hCk_RebinnedpLVar%uMeV_0",uVar);
      TH1F* hCkData = (TH1F*)DataFile.Get(HistoName);
      if(uVar==0){
        hDataDef = (TH1F*)DataFile.Get(HistoName);
      }
      for(unsigned uMom=0; uMom<NumMomBins; uMom++){

        const unsigned BootstrapN = 256;
        double CkVal,CkRelErr,CkErr;
        //CkVal = hCkData->GetBinContent(uMom+1);
        CkRelErr = hDataDef->GetBinError(uMom+1)/hDataDef->GetBinContent(uMom+1);
        CkErr = CkVal*CkRelErr;
        //printf("%f\n",CkRelErr);
        for(unsigned uBoot=0; uBoot<BootstrapN; uBoot++){
          CkVal = rangen.Gaus(hCkData->GetBinContent(uMom+1),CkErr);
          hData2D->Fill(hData2D->GetXaxis()->GetBinCenter(uMom+1),CkVal);
        }
        //hData2D->Fill(hData2D->GetXaxis()->GetBinCenter(uMom+1),hCkData->GetBinContent(uMom+1));
      }
    }
    for(unsigned uMom=0; uMom<NumMomBins; uMom++){
      TH1F* hProj = (TH1F*)hData2D->ProjectionY(TString::Format("hProj"),uMom+1,uMom+1);
      //gData[uMt].SetPoint(uMom,hData2D->GetXaxis()->GetBinCenter(uMom+1),hProj->GetMean());
      //gData[uMt].SetPointError(uMom,0,hProj->GetStdDev());

      double DataVal = hProj->GetMean();
      double DataErr = hProj->GetStdDev();
      double FitVal,FitMom;
      gFit[uMt].GetPoint(uMom,FitMom,FitVal);

      hData[uMt]->SetBinContent(uMom+1,DataVal);
      hData[uMt]->SetBinError(uMom+1,DataErr);

      double nsig = (DataVal-FitVal)/DataErr;
      gNsig[uMt].SetPoint(uMom,FitMom,nsig);

      delete hProj;
    }

    hDataInlet[uMt] = (TH1F*)hData[uMt]->Clone(TString::Format("hDataInlet_%u",uMt));
    delete hData2D;
  }//uMt

  int ColorNLO13 = kRed;//+1
  int ColorNLO19 = kCyan;
  //int ColorFT = kViolet+1;
  int ColorLO13 = kGreen;//+1
  int ColorNLO19_Red3S1 = kViolet+1;
  int ColorSigmaNLO = kAzure;//kBlue-7 or //kTeal
  int ColorXiLattice = kPink+1;//kOrange+2 or +1
  int ColorBaseline = kGray+3;//kGray+1

  int ColorInteraction;
  if(FitSettings.Contains("p11600-")){
    ColorInteraction = ColorNLO19;
  }
  else if(FitSettings.Contains("p1600-")){
    ColorInteraction = kRed;
  }
  else if(FitSettings.Contains("p131600-")){
    ColorInteraction = ColorNLO19_Red3S1;
  }


  for(unsigned uMt=0; uMt<NumMtBins; uMt++){
    gFit[uMt].SetFillColorAlpha(ColorInteraction,0.70);
    gFit[uMt].SetLineColor(ColorInteraction);
    gFit[uMt].SetLineWidth(0);

    gBl[uMt].SetFillColorAlpha(ColorBaseline,0.70);
    gBl[uMt].SetLineColor(ColorBaseline);
    gBl[uMt].SetLineWidth(0);

    gNsig[uMt].SetFillColorAlpha(ColorInteraction,0.70);
    gNsig[uMt].SetLineColor(ColorInteraction);
    gNsig[uMt].SetLineWidth(5);

    hData[uMt]->SetStats(false);
    hData[uMt]->SetTitle("; #it{k*} (MeV/#it{c}); #it{C}(#it{k*})");
    hData[uMt]->GetXaxis()->SetRangeUser(kMin,kMax);
    hData[uMt]->GetXaxis()->SetNdivisions(505);
    hData[uMt]->GetYaxis()->SetRangeUser(0.95, 2.1);
    //hData[uMt]->SetFillColor(kGray+1);
    SetStyleHisto_pLambda_PLB(hData[uMt],2,0,1,false);

    hDataInlet[uMt]->SetStats(false);
    hDataInlet[uMt]->SetTitle("; #it{k*} (MeV/#it{c}); #it{C}(#it{k*})");
    hDataInlet[uMt]->SetMarkerSize(hData[uMt]->GetMarkerSize()*0.67);
    hDataInlet[uMt]->SetLineWidth(hData[uMt]->GetLineWidth()*0.67);
    hDataInlet[uMt]->GetXaxis()->SetRangeUser(kMin,kMax);
    SetStyleHisto_pLambda_PLB(hDataInlet[uMt],2,0,2,false);
    hDataInlet[uMt]->GetXaxis()->SetNdivisions(505);
    hDataInlet[uMt]->GetYaxis()->SetNdivisions(505);
    hDataInlet[uMt]->GetYaxis()->SetRangeUser(0.975, 1.085);
  }


  TH1F* hAxisNsig = new TH1F("hAxisNsig", "hAxisNsig", NumMomBins,kMin,kMax);
  hAxisNsig->SetStats(false);
  hAxisNsig->GetXaxis()->SetRangeUser(kMin,kMax);
  hAxisNsig->SetTitle("; #it{k*} (MeV/#it{c}); #it{n_{#sigma}}");
  hAxisNsig->GetYaxis()->SetRangeUser(-3.7,3.7);
  hAxisNsig->GetXaxis()->SetNdivisions(505);
  SetStyleHisto_pLambda_PLB(hAxisNsig,2,0,2,false);
  hAxisNsig->GetYaxis()->SetNdivisions(505);
  //hAxisNsig->SetGridy(true);

  TLegend* legend = new TLegend(0.55,0.6,0.95,0.9);//lbrt
  legend->AddEntry(hData[0],"p#Lambda data (stat+syst)");
  legend->AddEntry(&gFit[0],"Fit");
  legend->AddEntry(&gBl[0],"Baseline");

  for(unsigned uMt=0; uMt<NumMtBins; uMt++){

    DLM_SubPads* DlmPad;
    //PLB
    DlmPad = new DLM_SubPads(1080,1080);
    DlmPad->AddSubPad(0,1,0.5,1);
    DlmPad->AddSubPad(0,1,0.25,0.5);
    DlmPad->AddSubPad(0,1,0.,0.25);
    DlmPad->SetMargin(0,0.15,0.02,0.0,0.02);
    DlmPad->SetMargin(1,0.15,0.02,0.0,0.0);
    DlmPad->SetMargin(2,0.15,0.02,0.1,0.0);

    DlmPad->cd(0);
    hData[uMt]->Draw();
    gFit[uMt].Draw("3L same");
    gBl[uMt].Draw("3L same");
    legend->Draw("same");

    DlmPad->cd(1);
    hDataInlet[uMt]->Draw();
    gFit[uMt].Draw("3L same");
    gBl[uMt].Draw("3L same");

    DlmPad->cd(2);
    DlmPad->SetGridy(2);
    hAxisNsig->Draw();
    gNsig[uMt].Draw("L same");

    DlmPad->GetCanvas()->SaveAs(OutDirName+TString::Format("DlmPad_%u.pdf",uMt));
  }

  TFile fOutput(OutDirName+OutFileName,"recreate");
  for(unsigned uMt=0; uMt<NumMtBins; uMt++){
    hData[uMt]->Write();
    gFit[uMt].Write();
    gBl[uMt].Write();
    gNsig[uMt].Write();
  }
  //hData->Write();

  delete hAxisNsig;
  delete legend;
  delete [] gFit;
  delete [] gBl;
  delete [] gNsig;
  for(unsigned uMt=0; uMt<NumMtBins; uMt++){
    delete hData[uMt];
  }
  delete [] hData;
  //delete [] gData;
  //delete hData;
}


//plots the source distribution out of CECA, using Jaime's parameters
//can be done for pp and pL, we compare the cases of Gaussian, RSM, CECA (<100 MeV), CECA (kstar dep)
//we further plot the differences onto the correlation function
void Plot_SourceFunMt(){
  
}


int JAIME_THESIS(int argc, char *argv[]){

  //Plot_CorrFunFits(0,"cl2741_p11600-a1201-e11-f0");
  //Plot_CorrFunFits(0,"cl2741_p131600-a1201-e11-f2");
  Plot_SourceFunMt();

  return 0;
}
