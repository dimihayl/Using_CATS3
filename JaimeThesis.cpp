
#include "JaimeThesis.h"
#include "EnvVars.h"
#include "FemtoBoyzScripts.h"
#include "DLM_SubPads.h"
#include "CommonAnaFunctions.h"
#include "DLM_Source.h"
#include "CECA.h"
#include "TREPNI.h"
#include "CATSconstants.h"
#include "DLM_RootWrapper.h"
#include "DLM_Histo.h"

#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TLegend.h"
#include "TF1.h"


//plots the correlation functions (based on data) and the final fit that
//Jaime gave me
//mode == 0; default
void Plot_CorrFunFits(const int& mode, TString FitSettings){

  const unsigned NumDataVar = 43;
  const unsigned NumMtBins = 6;
  double kMin = 0;
  double kMax = 432;
  unsigned NumMomBins = 36;
  if(FitSettings.Contains("p0-")){
    kMax = 204;
    NumMomBins = 17;
  }

  TString InDirName = TString::Format("%s/Jaime/Thesis/Input/v1/",GetCernBoxDimi());
  TString InFileName = "ck-fit_"+FitSettings+".root";
  TString OutDirName = TString::Format("%s/Jaime/Thesis/Output/v1/",GetCernBoxDimi());
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
  else if(FitSettings.Contains("pot0-")){
    ColorInteraction = ColorNLO19_Red3S1;
  }
  else if(FitSettings.Contains("p0-")){
    ColorInteraction = kRed;
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


void WhatAreOurLambdaPars(){



  const bool SourcePaper = true;
  std::vector<int> lampar_p_var;
  std::vector<int> lampar_l_var;
  if(SourcePaper){
    //actually for source paper it is mt diff, but we ignore to spare details, so it is 0/1/2 within syst (0 is def)
    lampar_p_var.push_back(0);
    lampar_p_var.push_back(1);
    lampar_p_var.push_back(2);

    //0 is for the source paper, with systematics 1 and 2
    //lampar_l_var.push_back(0);
    //lampar_l_var.push_back(1);
    //lampar_l_var.push_back(2);
    std::vector<int> var_A = {0,1,2};
    std::vector<int> var_B = {0,1,2};
    for (int & iA : var_A)
        for (int & iB : var_B)
          lampar_l_var.push_back(iA + 10 * iB);
  }
  else{
    lampar_p_var.push_back(0);
    lampar_p_var.push_back(1);
    lampar_p_var.push_back(2);

    // CBA: A:  is the S0:L, 0 is 1:3, 2 is 40% (kind of our default
  	//	    looking at best fits) 1 is 27%
  	//      B:  the Xi feed-down, 1 and 2 are some old settings,
  	//	    now working with either 0 (fractions averaged over pT)
  	//          or 3, which is our new default, where the fraction of
  	//	    lambdas is lowered based on the avg pT of the femto pairs
  	//      C: the purity, 4 is 95.3%, 5 is 96.3%, these values are based on the PLB result
  	std::vector<int> var_A = {0,1,2};
  	std::vector<int> var_B = {0,3};
  	std::vector<int> var_C = {5,4};

  	for (int & iA : var_A)
  	 for (int & iB : var_B)
  	  for (int & iC : var_C)
  		  lampar_l_var.push_back(iA + 10*iB + 100*iC);
  }


  for(unsigned ulp=0; ulp<lampar_p_var.size(); ulp++){
    for(unsigned ull=0; ull<lampar_l_var.size(); ull++){
      if(!ulp && !ull){
        printf("DEFAULT:\n");
      }
      DLM_CommonAnaFunctions AnalysisObject;
      double lam_pars_pl[5];
      AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());
      AnalysisObject.SetUpLambdaPars_pL("pp13TeV_HM_Dec19",lampar_p_var.at(ulp),lampar_l_var.at(ull),lam_pars_pl);

      double SUM=0;
      for(unsigned uLam=0; uLam<5; uLam++){
          //printf("λ(pΛ)_%u = %.1f\n",uLam,lam_pars_pl[uLam]*100.);
          SUM+=lam_pars_pl[uLam]*100.;
      }
      //printf("SUM: %.1f\n------------\n",SUM);

      //the thing in brackets it the value to compare with pL paper, we divide by the misid lambda, but multiply by the proton pur, as we dont correct it
      if(SourcePaper){
        printf(" gen = %.2f%%; S0 = %.2f%%; Xim=%.2f%%; flt=%.2f%%; fake=%.2f%%; totflt=%.2f%%\n",
        lam_pars_pl[0]*100.,lam_pars_pl[1]*100.,lam_pars_pl[2]*100.,
        lam_pars_pl[3]*100.,lam_pars_pl[4]*100.,lam_pars_pl[4]*100.+lam_pars_pl[3]*100.);
      }
      else{
        printf(" gen = %.2f%% (%.2f); S0 = %.2f%%; Xim/Xi0=%.2f%%; flt=%.2f%%; fake=%.2f%%; totflt=%.2f%%\n",
        lam_pars_pl[0]*100.,lam_pars_pl[0]/(1.-lam_pars_pl[4])*99.4,lam_pars_pl[1]*100.,lam_pars_pl[2]*100.,
        lam_pars_pl[3]*100.-lam_pars_pl[2]*100.,lam_pars_pl[4]*100.,lam_pars_pl[4]*100.+lam_pars_pl[3]*100.-lam_pars_pl[2]*100.);
      }

    }
  }





}

void Plot_EffSources(){
  TString OutDirName = TString::Format("%s/Jaime/Thesis/Output/v1/",GetCernBoxDimi());

  const double r_pp = 1.28;
  const double r_pL = 1.41;

  TF1* fSource_pp = new TF1("fSource_pp",GaussSourceTF1,0,8,1);
  fSource_pp->SetParameter(0,r_pp);
  fSource_pp->SetLineWidth(6);
  fSource_pp->SetLineColor(kBlue+1);
  fSource_pp->SetNpx(1024);


  TF1* fSource_pL = new TF1("fSource_pL",GaussSourceTF1,0,8,1);
  fSource_pL->SetParameter(0,r_pL);
  fSource_pL->SetLineWidth(5);
  fSource_pL->SetLineColor(kRed+1);
  fSource_pL->SetNpx(1024);

  TH1F* hAxisS = new TH1F("hAxisS", "hAxisS", 128, 0, 8);
  hAxisS->SetStats(false);
  hAxisS->SetTitle("");
  hAxisS->GetXaxis()->SetTitle("r (fm)");
  hAxisS->GetXaxis()->SetTitleSize(0.06);
  hAxisS->GetXaxis()->SetLabelSize(0.06);
  hAxisS->GetXaxis()->SetTitleOffset(1.3);
  hAxisS->GetXaxis()->SetLabelOffset(0.02);

  hAxisS->GetYaxis()->SetTitle("S_{4#pi} (1/fm)");
  hAxisS->GetYaxis()->SetTitleSize(0.06);
  hAxisS->GetYaxis()->SetLabelSize(0.06);
  hAxisS->GetYaxis()->SetTitleOffset(1.10);

  //hAxisS->GetXaxis()->SetNdivisions(506);
  hAxisS->GetXaxis()->SetRangeUser(0,8);
  hAxisS->GetYaxis()->SetRangeUser(0,0.4);

  TLegend* leg = new TLegend(0.65,0.75,0.95,0.95);//lbrt
  leg->SetName("leg");
  leg->SetTextSize(0.045);
  leg->AddEntry(fSource_pp, "p#minusp source (1.28 fm)");
  leg->AddEntry(fSource_pL, "p#minus#Lambda source (1.41 fm)");

  TCanvas* c_Sr = new TCanvas("c_Sr", "c_Sr", 1);
  c_Sr->cd(0);
  c_Sr->SetCanvasSize(1280, 720);
  c_Sr->SetMargin(0.15,0.05,0.2,0.05);//lrbt

  hAxisS->Draw("axis");
  fSource_pp->Draw("same");
  fSource_pL->Draw("same");
  leg->Draw("same");
  c_Sr->SaveAs(OutDirName+TString::Format("c_Sr.pdf"));
}

//for pp
void CompareSources_pp(){
  TString OutDirName = TString::Format("%s/Jaime/Thesis/Output/v1/",GetCernBoxDimi());
  const double rcore = 1.14;
  const double rSP = 1.0;
  const double MomSpread = 850;
  const double TIMEOUT = 15*8;

  printf("rcore(rSP) = %.2f (%.2f)\n",rcore,rSP);

  TF1* fSP = new TF1("fSP",GaussSourceTF1,0,8,1);
  fSP->SetParameter(0,rSP);
  fSP->SetLineWidth(6);
  fSP->SetLineColor(kBlack);
  fSP->SetLineStyle(2);
  fSP->SetNpx(1024);

  TF1* fcore = new TF1("fcore",GaussSourceTF1,0,8,1);
  fcore->SetParameter(0,rcore);
  fcore->SetLineWidth(6);
  fcore->SetLineColor(kBlue+1);
  fcore->SetNpx(1024);

  DLM_CleverMcLevyResoTM RsmFlat;
  SetUp_RSMflat_pp(RsmFlat);
  double reff_RsmFlat = GetReff(RsmFlat,rcore);
  TF1* freff_RsmFlat = new TF1("freff_RsmFlat",GaussSourceTF1,0,8,1);
  freff_RsmFlat->SetParameter(0,reff_RsmFlat);
  freff_RsmFlat->SetLineWidth(6);
  freff_RsmFlat->SetLineColor(kBlue+1);
  freff_RsmFlat->SetNpx(1024);

  TH1F* h_reff_RsmFlat = new TH1F("h_reff_RsmFlat","h_reff_RsmFlat",128,0,32);
  for(unsigned uBin=0; uBin<h_reff_RsmFlat->GetNbinsX(); uBin++){
    double RAD = h_reff_RsmFlat->GetBinCenter(uBin+1);
    double parameters[2];
    parameters[0] = rcore;
    parameters[1] = 2.0;
    double SourceValue = RsmFlat.RootEval(&RAD, parameters);
    h_reff_RsmFlat->SetBinContent(uBin+1,SourceValue);
  }
  h_reff_RsmFlat->SetLineWidth(6);
  h_reff_RsmFlat->SetLineColor(kBlue+1);
  double reff_RsmFlat_2 = Get_reff(h_reff_RsmFlat);
  printf("reff_RsmFlat = %.2f (%.2f)\n",reff_RsmFlat,reff_RsmFlat_2);

  TREPNI Database(0);
  Database.SetSeed(11);
  std::vector<TreParticle*> ParticleList;
  ParticleList.push_back(Database.NewParticle("Proton"));
  ParticleList.push_back(Database.NewParticle("Pion"));
  ParticleList.push_back(Database.NewParticle("ProtonReso"));
  for(TreParticle* prt : ParticleList){
    if(prt->GetName()=="Proton"){
      prt->SetMass(Mass_p);
      prt->SetAbundance(35.78);
    }
    else if(prt->GetName()=="Pion"){
      prt->SetMass(Mass_pic);
      prt->SetAbundance(0);
    }
    else if(prt->GetName()=="ProtonReso"){
      prt->SetMass(1362);
      prt->SetAbundance(64.22);
      prt->SetWidth(hbarc/1.65);

      prt->NewDecay();
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Proton"));
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
      prt->GetDecay(0)->SetBranching(100);
    }
    prt->SetPtPz(prt->GetMass()*MomSpread*0.001,prt->GetMass()*MomSpread*0.001);
  }
  std::vector<std::string> ListOfParticles;
  ListOfParticles.push_back("Proton");
  ListOfParticles.push_back("Proton");
  CECA Ivana(Database,ListOfParticles);
  Ivana.Ghetto_NumRadBins = 512;
  Ivana.Ghetto_NumRadBins = 512;
  Ivana.SetDisplacementT(rSP);
  Ivana.SetDisplacementZ(rSP);
  Ivana.SetTargetStatistics(10);
  Ivana.SetEventMult(2);
  Ivana.SetSourceDim(2);
  Ivana.SetThreadTimeout(TIMEOUT);
  Ivana.SetFemtoRegion(100);
  Ivana.GHETTO_EVENT = true;
  Ivana.GoBabyGo(0);

  Ivana.GhettoFemto_rstar->ComputeError();
  Ivana.GhettoFemto_rstar->ScaleToIntegral();
  Ivana.GhettoFemto_rstar->ScaleToBinSize();
  TH1F* h_reff_Ceca = Convert_DlmHisto_TH1F(Ivana.GhettoFemto_rstar,"h_reff_Ceca");
  h_reff_Ceca->SetLineColor(kMagenta);
  h_reff_Ceca->SetLineWidth(6);

  Ivana.GhettoFemto_rcore->ComputeError();
  Ivana.GhettoFemto_rcore->ScaleToIntegral();
  Ivana.GhettoFemto_rcore->ScaleToBinSize();
  TH1F* h_rcore_Ceca = Convert_DlmHisto_TH1F(Ivana.GhettoFemto_rcore,"h_rcore_Ceca");
  h_rcore_Ceca->SetLineColor(kMagenta);
  h_rcore_Ceca->SetLineWidth(6);

  double rcore_Ceca = Get_reff(h_rcore_Ceca);
  double reff_Ceca = Get_reff(h_reff_Ceca);
  printf("rcore_Ceca = %.2f\n",rcore_Ceca);
  printf("reff_Ceca = %.2f\n",reff_Ceca);


  TH1F* hAxisC = new TH1F("hAxisC", "hAxisC", 128, 0, 8);
  hAxisC->SetStats(false);
  hAxisC->SetTitle("");
  hAxisC->GetXaxis()->SetTitle("r (fm)");
  hAxisC->GetXaxis()->SetTitleSize(0.055);
  hAxisC->GetXaxis()->SetLabelSize(0.055);
  hAxisC->GetXaxis()->SetTitleOffset(1.25);
  hAxisC->GetXaxis()->SetLabelOffset(0.02);

  hAxisC->GetYaxis()->SetTitle("S_{core} (1/fm)");
  hAxisC->GetYaxis()->SetTitleSize(0.055);
  hAxisC->GetYaxis()->SetLabelSize(0.055);
  hAxisC->GetYaxis()->SetTitleOffset(1.25);

  //hAxisC->GetXaxis()->SetNdivisions(506);
  hAxisC->GetXaxis()->SetRangeUser(0,8);
  hAxisC->GetYaxis()->SetRangeUser(0,0.45);

  TLegend* legC = new TLegend(0.5,0.775,0.975,0.975);//lbrt
  legC->SetName("leg");
  legC->SetTextSize(0.04);
  legC->AddEntry(fSP, "Gaussian of 1.00 fm");
  legC->AddEntry(fcore, "Gaussian of 1.14 fm");
  legC->AddEntry(h_rcore_Ceca, "CECA core with r_{d}=1.0 fm");

  TCanvas* c_SrC = new TCanvas("c_SrC", "c_SrC", 1);
  c_SrC->cd(0);
  c_SrC->SetCanvasSize(960, 720);
  c_SrC->SetMargin(0.16,0.05,0.15,0.05);//lrbt

  hAxisC->Draw("axis");
  legC->Draw("same");
  fSP->Draw("same");
  fcore->Draw("same");
  h_rcore_Ceca->Draw("same");
  c_SrC->SaveAs(OutDirName+"CompareSources_pp_core.pdf");



  TH1F* hAxisEff = new TH1F("hAxisEff", "hAxisEff", 128, 0, 8);
  hAxisEff->SetStats(false);
  hAxisEff->SetTitle("");
  hAxisEff->GetXaxis()->SetTitle("r (fm)");
  hAxisEff->GetXaxis()->SetTitleSize(0.055);
  hAxisEff->GetXaxis()->SetLabelSize(0.055);
  hAxisEff->GetXaxis()->SetTitleOffset(1.25);
  hAxisEff->GetXaxis()->SetLabelOffset(0.02);

  hAxisEff->GetYaxis()->SetTitle("S (1/fm)");
  hAxisEff->GetYaxis()->SetTitleSize(0.055);
  hAxisEff->GetYaxis()->SetLabelSize(0.055);
  hAxisEff->GetYaxis()->SetTitleOffset(1.25);

  //hAxisEff->GetXaxis()->SetNdivisions(506);
  hAxisEff->GetXaxis()->SetRangeUser(0,8);
  hAxisEff->GetYaxis()->SetRangeUser(0,0.45);

  TLegend* legEff = new TLegend(0.5,0.775,0.975,0.975);//lbrt
  legEff->SetName("leg");
  legEff->SetTextSize(0.04);
  legEff->AddEntry(h_reff_RsmFlat, "Total source from RSM");
  legEff->AddEntry(h_reff_Ceca, "Total source from CECA");

  TCanvas* c_SrEff = new TCanvas("c_SrEff", "c_SrEff", 1);
  c_SrEff->cd(0);
  c_SrEff->SetCanvasSize(960, 720);
  c_SrEff->SetMargin(0.16,0.05,0.15,0.05);//lrbt
  hAxisEff->Draw("axis");
  legEff->Draw("same");
  h_reff_RsmFlat->Draw("same");
  h_reff_Ceca->Draw("same");

  c_SrEff->SaveAs(OutDirName+"CompareSources_pp_eff.pdf");



  TFile fOutput(OutDirName+"CompareSources_pp.root","recreate");
  fcore->Write();
  fSP->Write();
  h_reff_RsmFlat->Write();
  freff_RsmFlat->Write();
  h_rcore_Ceca->Write();
  h_reff_Ceca->Write();
  c_SrC->Write();
  c_SrEff->Write();
}

int JAIME_THESIS(int argc, char *argv[]){

  //Plot_CorrFunFits(0,"cl2741_p11600-a1201-e11-f0");
  //Plot_CorrFunFits(0,"cl2741_p131600-a1201-e11-f2");
  //Plot_CorrFunFits(0,"cl2741_p0-a0-e11-f0");
  //Plot_SourceFunMt();
  //WhatAreOurLambdaPars();
  //Plot_EffSources();
  CompareSources_pp();

  return 0;
}
