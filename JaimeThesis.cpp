
#include "JaimeThesis.h"
#include "EnvVars.h"

#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"


//plots the correlation functions (based on data) and the final fit that
//Jaime gave me
//mode == 0; default
void Plot_CorrFunFits(const int& mode, TString FitSettings){

  const unsigned NumDataVar = 43;
  const unsigned NumMtBins = 6;
  double kMin = 0;
  double kMax = 432;
  unsigned NumMomBins = 36;

  TString InDirName = TString::Format("%s/Thesis/Input/",GetCernBoxDimi());
  TString InFileName = "ck-fit_"+FitSettings+".root";
  TString OutFileName = "ck-plot_"+FitSettings+".root";

  TGraphErrors* gFit = new TGraphErrors[NumMtBins];
  TGraphErrors* gBl = new TGraphErrors[NumMtBins];
  TGraphErrors* gData = new TGraphErrors[NumMtBins];

  TFile InFile(InDirName+InFileName,"read");
  for(unsigned uMt=0; uMt<NumMtBins; uMt++){
    gFit[uMt].SetName(TString::Format("gFit_%u",uMt));
    gBl[uMt].SetName(TString::Format("gBl_%u",uMt));
    gData[uMt].SetName(TString::Format("gData_%u",uMt));
    TH2F* gFitMt = (TH2F*)InFile.Get(TString::Format("gFit mT%u",uMt));
    for(unsigned uMom=0; uMom<NumMomBins; uMom++){
      TH1F* hProj = (TH1F*)gFitMt->ProjectionY(TString::Format("hProj"),uMom+1,uMom+1);
      gFit[uMt].SetPoint(uMom,hProj->GetMean());
      gFit[uMt].SetPointError(uMom,0,hProj->GetStdDev());
      delete hProj;

      hProj = (TH1F*)gBlMt->ProjectionY(TString::Format("hProj"),uMom+1,uMom+1);
      gBl[uMt].SetPoint(uMom,hProj->GetMean());
      gBl[uMt].SetPointError(uMom,0,hProj->GetStdDev());
      delete hProj;
    }

    for(unsigned uVar=0; uVar<NumDataVar; uVar++){
      TFile DataFile(TString::Format("%s/CatsFiles/ExpData/Bernie_Source/pLData/mTBin_%u/CFOutput_mT_pLVar%u_HM_%u.root",GetCernBoxDimi(),uMt+1,uVar,uMt),"read");
      for(unsigned uMom=0; uMom<NumMomBins; uMom++){
        
      }
    }
  }

  delete [] gFit;
  delete [] gBl;
  delete [] gData;

}


//plots the source distribution out of CECA, using Jaime's parameters
//can be done for pp and pL, we compare the cases of Gaussian, RSM, CECA (<100 MeV), CECA (kstar dep)
//we further plot the differences onto the correlation function
void Plot_SourceFunMt(){

}


int JAIME_THESIS(int argc, char *argv[]){

  Plot_CorrFunFits(0,"cl2741_p11600-a1201-e11-f0");
  Plot_SourceFunMt();

  return 0;
}
