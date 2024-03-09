
#include <stdio.h>

#include "CATS.h"
#include "CATStools.h"
#include "DLM_Ck.h"
#include "DLM_CkDecomposition.h"
#include "CommonAnaFunctions.h"
#include "DLM_WfModel.h"
#include "DLM_Source.h"
#include "CECA.h"
#include "CECA_Paper.h"
#include "CECA.h"
#include "TREPNI.h"
#include "CATSconstants.h"
#include "DLM_RootWrapper.h"

#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TGraph.h"
#include "TString.h"
#include "TList.h"
#include "TROOT.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TLegend.h"

#include "EnvVars.h"
#include "HighMtProblem.h"




//study the effect on Ck of the interaction / source(k*),Levy
//with the RC in p-wave, we see effect at large k*, the bump ~200 MeV can dissapear if we have enough RC
//however, we see this bump in pretty much all mT bins, theoretically we expect to see it for all sources. So all is fine in p-waves.
//for the s-wave the effect of RC is at low k*, it changes peak amplitude and even slightly the position.
//at first look it seems also the s-wave effects all radii, so we should have seen effects also in the other bins
void PlayWithProtons(){

  std::vector<float> rad;
  rad.push_back(0.95);
  rad.push_back(1.3);
  const int NumMomBins = 160;
  const double kMin = 0;
  const double kMax = 320;
  CATS KittyAV18;
  KittyAV18.SetMomBins(NumMomBins,kMin,kMax);

  CATS Kitty;
  Kitty.SetMomBins(NumMomBins,kMin,kMax);

  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());
  AnalysisObject.SetUpCats_pp(Kitty,"pp_AV18_Toy","Gauss",0,0);
  AnalysisObject.SetUpCats_pp(KittyAV18,"AV18","Gauss",0,0);

  //Kitty.SetShortRangePotential(0,0,1,-2000);
  //Kitty.SetShortRangePotential(0,0,2,0.24);

  //Kitty.SetShortRangePotential(0,0,1,-2700);
  //Kitty.SetShortRangePotential(0,0,2,0.25);

  Kitty.SetShortRangePotential(0,0,1,-2700);
  Kitty.SetShortRangePotential(0,0,2,0.75);

  //Kitty.SetShortRangePotential(3,1,1,-4000);
  //Kitty.SetShortRangePotential(3,1,2,0.4);

  //Kitty.SetShortRangePotential(3,1,1,-1500);
  //Kitty.SetShortRangePotential(3,1,2,0.4);

  //Kitty.SetQ1Q2(0);
  //Kitty.SetQuantumStatistics(false);
  //KittyAV18.SetQ1Q2(0);
  //KittyAV18.SetQuantumStatistics(false);

  Kitty.KillTheCat();
  KittyAV18.KillTheCat();

  Kitty.SetNotifications(CATS::nWarning);
  KittyAV18.SetNotifications(CATS::nWarning);

  const int NumCk = 1;
  TGraph* gAV18 = new TGraph [rad.size()];
  TGraph** gCk = new TGraph*[NumCk];
  TString* CkInfo = new TString [NumCk];

  TGraph g1S0_AV18;
  TGraph* g1S0_MODY = new TGraph [NumCk];
  TGraph g3P2_AV18;
  TGraph* g3P2_MODY = new TGraph [NumCk];

  TGraph gPs1S0_AV18;
  TGraph* gPs1S0_MODY = new TGraph [NumCk];
  TGraph gPs3P2_AV18;
  TGraph* gPs3P2_MODY = new TGraph [NumCk];

  CkInfo[0] = "Modified";

  g1S0_AV18.SetName("g1S0_AV18");
  g3P2_AV18.SetName("g3P2_AV18");
  gPs1S0_AV18.SetName("gPs1S0_AV18");
  gPs3P2_AV18.SetName("gPs3P2_AV18");
  for(unsigned uRad=0; uRad<rad.size(); uRad++)
    gAV18[uRad].SetName(TString::Format("gAV18_%.2f",rad.at(uRad)));

  for(unsigned uRstar=0; uRstar<1024; uRstar++){
    double rstar = double(uRstar+1)*4./1000.;
    double V1S0 = KittyAV18.EvaluateThePotential(0,0,5,rstar);
    double V3P2 = KittyAV18.EvaluateThePotential(3,1,5,rstar);
    g1S0_AV18.SetPoint(uRstar,rstar,V1S0);
    g3P2_AV18.SetPoint(uRstar,rstar,V3P2);
  }
  for(unsigned uMom=0; uMom<NumMomBins; uMom++){
    double MOM = Kitty.GetMomentum(uMom);
    double Ps1S0 = KittyAV18.GetPhaseShift(uMom,0,0);
    double Ps3P2 = KittyAV18.GetPhaseShift(uMom,3,1);
    gPs1S0_AV18.SetPoint(uMom,MOM,Ps1S0);
    gPs3P2_AV18.SetPoint(uMom,MOM,Ps3P2);
  }

  for(unsigned uCk=0; uCk<NumCk; uCk++){
    gCk[uCk] = new TGraph [rad.size()];
    g1S0_MODY[uCk].SetName(TString::Format("g1S0_MODY_%s",CkInfo[uCk].Data()));
    g3P2_MODY[uCk].SetName(TString::Format("g3P2_MODY_%s",CkInfo[uCk].Data()));
    gPs1S0_MODY[uCk].SetName(TString::Format("gPs1S0_MODY_%s",CkInfo[uCk].Data()));
    gPs3P2_MODY[uCk].SetName(TString::Format("gPs3P2_MODY_%s",CkInfo[uCk].Data()));
    for(unsigned uRad=0; uRad<rad.size(); uRad++){
      gCk[uCk][uRad].SetName(TString::Format("gCk_%s_%.2f",CkInfo[uCk].Data(),rad.at(uRad)));
      Kitty.SetAnaSource(0,rad.at(uRad));
      KittyAV18.SetAnaSource(0,rad.at(uRad));
      Kitty.KillTheCat();
      KittyAV18.KillTheCat();
      for(unsigned uMom=0; uMom<NumMomBins; uMom++){
        double MOM = Kitty.GetMomentum(uMom);
        double Ck_val = Kitty.GetCorrFun(uMom);
        gCk[uCk][uRad].SetPoint(uMom,MOM,Ck_val);
        if(uCk==0){
          Ck_val = KittyAV18.GetCorrFun(uMom);
          gAV18[uRad].SetPoint(uMom,MOM,Ck_val);
        }
      }
    }//uRad
    for(unsigned uRstar=0; uRstar<1024; uRstar++){
      double rstar = double(uRstar+1)*4./1000.;
      double V1S0 = Kitty.EvaluateThePotential(0,0,5,rstar);
      double V3P2 = Kitty.EvaluateThePotential(3,1,5,rstar);
      g1S0_MODY[uCk].SetPoint(uRstar,rstar,V1S0);
      g3P2_MODY[uCk].SetPoint(uRstar,rstar,V3P2);
    }
    for(unsigned uMom=0; uMom<NumMomBins; uMom++){
      double MOM = Kitty.GetMomentum(uMom);
      double Ps1S0 = Kitty.GetPhaseShift(uMom,0,0);
      double Ps3P2 = Kitty.GetPhaseShift(uMom,3,1);
      gPs1S0_MODY[uCk].SetPoint(uMom,MOM,Ps1S0);
      gPs3P2_MODY[uCk].SetPoint(uMom,MOM,Ps3P2);
    }
  }

  TString OutputFileName = TString::Format("%s/HighMtProblem/PlayWithProtons/fOutput.root",GetFemtoOutputFolder());
  TFile fOutput(OutputFileName,"recreate");

  for(unsigned uRad=0; uRad<rad.size(); uRad++){
    gAV18[uRad].Write();
    for(unsigned uCk=0; uCk<NumCk; uCk++){
      gCk[uCk][uRad].Write();
    }
  }

  g1S0_AV18.Write();
  gPs1S0_AV18.Write();
  g3P2_AV18.Write();
  gPs3P2_AV18.Write();
  for(unsigned uCk=0; uCk<NumCk; uCk++){
    g1S0_MODY[uCk].Write();
    gPs1S0_MODY[uCk].Write();
    g3P2_MODY[uCk].Write();
    gPs3P2_MODY[uCk].Write();
  }


  for(unsigned uRad=0; uRad<rad.size(); uRad++){
    delete [] gCk[uRad];
  }

  delete [] g1S0_MODY;
  delete [] g3P2_MODY;
  delete [] gPs1S0_MODY;
  delete [] gPs3P2_MODY;
  delete [] gAV18;
  delete [] gCk;
  delete [] CkInfo;
}





double pp_simple_fitter(double* x, double* pars){
    return 0;
}

//gets the data from the corresponding mt bin, creates a simple cats fitter (no feed-down, lam_genuine is a fit parameter.
//The CATS object has to be already set up !!!
//all passed pointers will be used to create objects that need to be cleaned up !!!
//WhichMtBin counts from zero
void SetUpProtonProtonAna(int WhichMtBin, CATS& Cat, DLM_Ck* DlmCk, DLM_CkDecomposition* CkDecomp, TH1F* hData, TF1* fData){

    DLM_CommonAnaFunctions AnalysisObject;
    AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());
    AnalysisObject.GetAliceExpCorrFun("pp13TeV_HM_BernieSource","pp","0",0,1,WhichMtBin);

    double MinMom = Cat.GetMomBinLowEdge(0);
    double MaxMom = Cat.GetMomBinUpEdge(Cat.GetNumMomBins()-1);
    int NumMomBins = Cat.GetNumMomBins();
    double dMom = (MaxMom-MinMom)/double(NumMomBins);
    NumMomBins *= 1.5;
    MaxMom = MinMom + dMom*double(NumMomBins);
    DlmCk = new DLM_Ck(Cat.GetNumSourcePars(),0,Cat,NumMomBins,MinMom,MaxMom);
    DlmCk->SetCutOff(Cat.GetMomBinUpEdge(Cat.GetNumMomBins()-1),700);

    TH2F* hReso_pp = AnalysisObject.GetResolutionMatrix("pp13TeV_HM_DimiJun20","pp");
    CkDecomp = new DLM_CkDecomposition("pp",0,*DlmCk,hReso_pp);

    //get the phase space and add it. The add phase space copies the input into a DLM_Histo, so no issues deleting
    for(unsigned uMt=0; uMt<7; uMt++){
      TString FileName = TString::Format("%s/CatsFiles/ExpData/Bernie_Source/ppData/mTBin_%i/CFOutput_mT_ppVar0_HM_%i.root",
                                          GetCernBoxDimi(),uMt+1,uMt);
      TFile* inFile = new TFile(FileName,"read");
      //PARTICLES
      TList* list1_tmp = (TList*)inFile->Get("PairDist");
      TList* list2_tmp = (TList*)list1_tmp->FindObject("PairFixShifted");//there is also PairShifted ?? which one ??
      TH1F* hME_PP = (TH1F*)list2_tmp->FindObject(TString::Format("MEPart_mT_%i_FixShifted",uMt));
      list1_tmp = (TList*)inFile->Get("AntiPairDist");
      list2_tmp = (TList*)list1_tmp->FindObject("PairFixShifted");
      TH1F* hME_APAP = (TH1F*)list2_tmp->FindObject(TString::Format("MEAntiPart_mT_%i_FixShifted",uMt));
      gROOT->cd();
      TH1F* hPs_pp = (TH1F*)hME_PP->Clone("hPs_pp");
      hPs_pp->Add(hME_APAP);
      hPs_pp->GetXaxis()->SetLimits(hPs_pp->GetXaxis()->GetXmin()*1000.,hPs_pp->GetXaxis()->GetXmax()*1000.);
      CkDecomp->AddPhaseSpace(hPs_pp);
      delete hPs_pp;
      delete inFile;
    }

    //parameters:
    //[0] -> normalization
    //[1] -> lambda par
    //[2] -> source par 1 (percentage, how much of Gauss 1 we take)
    //[3] -> source par 2 (Gauss size 1)
    //[4] -> source par 3 (Gauss size 2)
    fData = new TF1("fData",pp_simple_fitter,Cat.GetMomBinLowEdge(0),Cat.GetMomBinUpEdge(Cat.GetNumMomBins()-1),5);

    fData->SetParameter(0,1);
    fData->SetParLimits(0,0.9,1.1);

    fData->SetParameter(1,0.67);
    fData->SetParLimits(1,0.5,0.9);

    fData->SetParameter(2,1);
    fData->SetParLimits(2,0.0,1);

    fData->SetParameter(3,1.2);
    fData->SetParLimits(3,0.9,1.5);

    fData->SetParameter(4,0.7);
    fData->SetParLimits(4,0.4,0.9);

}

    //1 = s only
    //11 = sp
    //1XXX -> AV18 with the magic trick
void Compare_pp_AV18_Haide_SingleGauss(int PotVar){
    std::vector<float> SrcSize = {0.6,0.8,1.0,1.2,1.4,1.6};
    const unsigned NumMomBins = 60;//does not work for 24 ???
    const double kMin = 0;
    const double kMax = 240;

    DLM_CommonAnaFunctions AnalysisObject;
    AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());

    CATS pp_AV18;
    pp_AV18.SetMomBins(NumMomBins,kMin,kMax);
    AnalysisObject.SetUpCats_pp(pp_AV18,"AV18","Gauss",PotVar/1000,0);
    if(PotVar%1000==1){
        pp_AV18.RemoveShortRangePotential(0,2);
        pp_AV18.RemoveShortRangePotential(1,1);
        pp_AV18.RemoveShortRangePotential(2,1);
        pp_AV18.RemoveShortRangePotential(3,1);
    }
    else if(PotVar%1000==11){
        pp_AV18.RemoveShortRangePotential(0,2);
    }


    printf("Johann time\n");
    CATS pp_John;
    pp_John.SetMomBins(NumMomBins,kMin,kMax);
    //Init_pp_Haidenbauer(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data(),pp_John,0);
    //
    AnalysisObject.SetUpCats_pp(pp_John,"Johann","Gauss",PotVar%1000,0);
    printf("...\n");

    TFile fOutput(TString::Format("%s/HighMtProblem/Compare_pp_AV18_Haide_SingleGauss/Compare_pp_AV18_Haide_SingleGauss_pv%i.root",GetFemtoOutputFolder(),PotVar),"recreate");
    for(float R0 : SrcSize){

        pp_AV18.SetAnaSource(0,R0);
        pp_John.SetAnaSource(0,R0);

        pp_AV18.KillTheCat();
        pp_John.KillTheCat();

        
        TGraph gAV18;
        gAV18.SetName(TString::Format("gAV18_%.2f",R0));
        TGraph gJohn;
        gJohn.SetName(TString::Format("gJohn_%.2f",R0));

        TGraph gRatio;
        gRatio.SetName(TString::Format("gRatio_%.2f",R0));

        for(unsigned uMom=0; uMom<NumMomBins; uMom++){
            double kstar = pp_AV18.GetMomentum(uMom);
            gAV18.SetPoint(uMom,kstar,pp_AV18.GetCorrFun(uMom));
            gJohn.SetPoint(uMom,kstar,pp_John.GetCorrFun(uMom));
            gRatio.SetPoint(uMom,kstar,pp_John.GetCorrFun(uMom)/pp_AV18.GetCorrFun(uMom));
        }

        gAV18.Write();
        gJohn.Write();
        gRatio.Write();

    }
}


CATS* DGE_Cat = NULL;
double DGE_Fitter(double* x, double* pars){
    DGE_Cat->SetAnaSource(0,pars[0],true);
    DGE_Cat->KillTheCat();
    return pars[1]*DGE_Cat->EvalCorrFun(*x);
}

//const + par*k^2
double Kstar_Source(double* pars){
    double& kstar = pars[0];//
    double& rad = pars[1];
    double ssize = pars[3];//the source size at 0 MeV
    double& slope = pars[4];
    ssize += slope*0.000001*kstar*kstar;
    return 4.*TMath::Pi()*rad*rad*pow(4.*TMath::Pi()*ssize*ssize,-1.5)*exp(-(rad*rad)/(4.*ssize*ssize));
    //return 4.*Pi*Radius*Radius*pow(4.*Pi*Size*Size,-1.5)*exp(-(Radius*Radius)/(4.*Size*Size));
}
DLM_Histo<float>* KstarDepSource;


//steps:
//generate a correlation function from two gaussians, one small one large.
//loop over several different ratios of the weight of each source
//each result for Ck fit with a normal Gauss source, see the result. To make the fit reallistic,
//generate dummy data with large errors at small k*. Do only theoretical stuff to save time
void Effects(){
    const double CoreSize = 0.5;
    const double EffSize = 1.5;
    std::vector<float> CoreWeight = {0.1, 0.3, 0.5};

    std::vector<float> KstarSlope = {10, 40, 70};


    const unsigned NumMomBins = 50;
    const double kMin = 0;
    const double kMax = 200;

    DLM_CommonAnaFunctions AnalysisObject;
    AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());


    CATS SingleCat;
    SingleCat.SetMomBins(NumMomBins,kMin,kMax);
    AnalysisObject.SetUpCats_pp(SingleCat,"AV18","Gauss",0,0);
    SingleCat.SetAnaSource(0,EffSize);
    SingleCat.SetNotifications(CATS::nWarning);
    DGE_Cat = &SingleCat;

    CATS DoubleCat;
    DoubleCat.SetMomBins(NumMomBins,kMin,kMax);
    AnalysisObject.SetUpCats_pp(DoubleCat,"AV18","DoubleGauss",0,0);


    KstarDepSource = new DLM_Histo<float>();
    KstarDepSource->SetUp(2);
    KstarDepSource->SetUp(0,60,0,240);//kstar
    KstarDepSource->SetUp(1,100,0.4,5.4);//source size

    CATS KstarCat;
    KstarCat.SetMomBins(NumMomBins,kMin,kMax);
    AnalysisObject.SetUpCats_pp(KstarCat,"AV18","NULL",0,0);
    CATSparameters cPars(CATSparameters::tSource, 2, true);
    KstarCat.SetAnaSource(Kstar_Source, cPars);
    //KstarCat.SetAnaSource(GaussSource, cPars);
    KstarCat.SetNormalizedSource(false);
    KstarCat.SetAutoNormSource(true);
    

    TRandom3 rangen(23);
    TFile fOutput(TString::Format("%s/HighMtProblem/Effects/Effects.root",GetFemtoOutputFolder()),"recreate");
    
    for(unsigned uW=0; uW<CoreWeight.size(); uW++){
        printf("w = %.0f\n",CoreWeight.at(uW)*100);
        DoubleCat.SetAnaSource(0,CoreSize);
        DoubleCat.SetAnaSource(1,EffSize);
        DoubleCat.SetAnaSource(2,CoreWeight.at(uW));
        DoubleCat.KillTheCat();
        
        TGraph* gCkDg = new TGraph();
        gCkDg->SetName(TString::Format("gCkDg_%.0f",CoreWeight.at(uW)*100));
        TH1F* hCkDg = new TH1F(TString::Format("hCkDg_%.0f",CoreWeight.at(uW)*100),TString::Format("hCkDg_%.0f",CoreWeight.at(uW)*100),NumMomBins,kMin,kMax);
        TF1* fCkSg = new TF1(TString::Format("fCkSg_%.0f",CoreWeight.at(uW)*100),DGE_Fitter,kMin,kMax,2);
        fCkSg->SetParameter(0,CoreWeight.at(uW)*CoreSize + (1.-CoreWeight.at(uW))*EffSize);
        fCkSg->SetParLimits(0, CoreSize, EffSize);
        fCkSg->SetParameter(1,1);
        fCkSg->SetParLimits(1, 0.9, 1.05);
    fCkSg->FixParameter(1,1);
        SingleCat.SetAnaSource(0,CoreWeight.at(uW)*CoreSize + (1.-CoreWeight.at(uW))*EffSize);
        SingleCat.KillTheCat();
        

        for(unsigned uBin=0; uBin<NumMomBins; uBin++){
            double Kstar = DoubleCat.GetMomentum(uBin);
            double CkVal = DoubleCat.GetCorrFun(uBin);
            double CkErr = 0.005*(1+60*exp(-pow(Kstar/75.,2.)));
            hCkDg->SetBinContent(uBin+1,fabs(rangen.Gaus(CkVal,CkErr)));
            hCkDg->SetBinError(uBin+1,CkErr);
            gCkDg->SetPoint(uBin, Kstar, CkVal);
        }

        hCkDg->Fit(fCkSg, "S, N, R, M");

        hCkDg->SetStats(false);
        hCkDg->SetTitle(TString::Format("%.0f %% Core",CoreWeight.at(uW)*100));
        hCkDg->GetXaxis()->SetTitle("#it{k*} (MeV/#it{c})");
        hCkDg->GetXaxis()->SetTitleSize(0.06);
        hCkDg->GetXaxis()->SetLabelSize(0.06);
        hCkDg->GetXaxis()->CenterTitle();
        hCkDg->GetXaxis()->SetTitleOffset(1.25);
        hCkDg->GetXaxis()->SetLabelOffset(0.02);
        hCkDg->GetXaxis()->SetRangeUser(0.0, 160);

        hCkDg->GetYaxis()->SetTitle("C(k*)");
        hCkDg->GetYaxis()->SetTitleSize(0.06);
        hCkDg->GetYaxis()->SetLabelSize(0.06);
        hCkDg->GetYaxis()->CenterTitle();
        hCkDg->GetYaxis()->SetTitleOffset(1.05);
        //hCkDg->GetXaxis()->SetNdivisions(506);
        //hCkDg->GetYaxis()->SetRangeUser(0.90, 1.80);
        hCkDg->SetLineWidth(5);

        fCkSg->SetLineColor(kRed);
        fCkSg->SetLineWidth(5);
        fCkSg->SetLineStyle(7);

        gCkDg->SetLineColor(kBlack);
        gCkDg->SetLineWidth(5);

        TLegend* lLegend = new TLegend(0.35,0.75,0.975,0.94);//lbrt
        lLegend->SetName(TString::Format("lLegend"));
        lLegend->SetTextSize(0.0425);
        lLegend->AddEntry(gCkDg,"C(k*) using Double Gauss");
        lLegend->AddEntry(hCkDg,"Dummy C(k*) sampled from black");
        lLegend->AddEntry(fCkSg,"Fit with effective single Gauss");
        
        TCanvas* plt = new TCanvas("plt","plt",1);
        plt->cd(0); plt->SetCanvasSize(1920, 1440); plt->SetMargin(0.13,0.05,0.16,0.06);//lrbt
        hCkDg->Draw("");
        gCkDg->Draw("same");
        fCkSg->Draw("same");
        lLegend->Draw("same");
        plt->SaveAs(TString::Format("%s/HighMtProblem/Effects/EffectsDG_%.0f.png",GetFemtoOutputFolder(),CoreWeight.at(uW)*100));


        gCkDg->Write();
        hCkDg->Write();
        fCkSg->Write();

        delete gCkDg;
        delete hCkDg;
        delete fCkSg;
        delete lLegend;
        delete plt;
    }

    for(unsigned uW=0; uW<KstarSlope.size(); uW++){
        printf("wks = %.0f\n",KstarSlope.at(uW));
        KstarCat.SetAnaSource(0,CoreSize);
        KstarCat.SetAnaSource(1,KstarSlope.at(uW));
        KstarCat.KillTheCat();
        
        TGraph* gCkKs = new TGraph();
        gCkKs->SetName(TString::Format("gCkKs_%.0f",KstarSlope.at(uW)));
        TH1F* hCkKs = new TH1F(TString::Format("hCkKs_%.0f",KstarSlope.at(uW)),TString::Format("hCkKs_%.0f",KstarSlope.at(uW)),NumMomBins,kMin,kMax);
        TF1* fCkKs = new TF1(TString::Format("fCkKs_%.0f",KstarSlope.at(uW)),DGE_Fitter,kMin,kMax,2);
        fCkKs->SetParameter(0,CoreSize);
        fCkKs->SetParLimits(0, CoreSize*0.8, CoreSize*5);
        fCkKs->SetParameter(1,1);
        fCkKs->SetParLimits(1, 0.9, 1.05);
        SingleCat.SetAnaSource(0,CoreSize*1.1);
        SingleCat.KillTheCat();
        
        for(unsigned uBin=0; uBin<NumMomBins; uBin++){
            double Kstar = KstarCat.GetMomentum(uBin);
            double CkVal = KstarCat.GetCorrFun(uBin);
            //printf("C(%.0f) = %.3e\n",Kstar,CkVal);
            double CkErr = 0.005*(1+60*exp(-pow(Kstar/75.,2.)));
            hCkKs->SetBinContent(uBin+1,fabs(rangen.Gaus(CkVal,CkErr)));
            hCkKs->SetBinError(uBin+1,CkErr);
            gCkKs->SetPoint(uBin, Kstar, CkVal);
        }

        hCkKs->Fit(fCkKs, "S, N, R, M");
        gCkKs->Write();
        hCkKs->Write();
        fCkKs->Write();

        delete gCkKs;
        delete hCkKs;
        delete fCkKs;
    }
}

//the goal here is to see if CATS can handle this
//to achieve it: we first evaluate C(k*) by changing the source size of a CATS object, for different momenta
//we evaluate the CATS object separately for each momentum. This we do to get the "true" correlation
//after that, we try to reproduce the result with a normal CATS object
void KstarSource_Test1(){
  double BaseSource = 0.9;
  double KstarLin = 0;//3
  double KstarPar = 60;//60
  double kMin = 0;
  double kMax = 200;
  unsigned kSteps = 50;

  TH1F* hCkTrue = new TH1F("hCkTrue","hCkTrue",kSteps,kMin,kMax);
  TH1F* hCkMin = new TH1F("hCkMin","hCkMin",kSteps,kMin,kMax);
  TH1F* hCkMax = new TH1F("hCkMax","hCkMax",kSteps,kMin,kMax);

  //MIN SOURCE
  CATS MouseJerry;
  MouseJerry.SetMomBins(kSteps,kMin,kMax);
  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());
  AnalysisObject.SetUpCats_pp(MouseJerry,"AV18","Gauss",0,0);
  MouseJerry.SetAnaSource(0,BaseSource);
  MouseJerry.KillTheCat();

  //MAX SOURCE
  CATS CatTom;
  CatTom.SetMomBins(kSteps,kMin,kMax);
  AnalysisObject.SetUpCats_pp(CatTom,"AV18","Gauss",0,0);
  CatTom.SetAnaSource(0,BaseSource + KstarLin*1e-3*kMax + KstarPar*1e-6*kMax*kMax);
  CatTom.KillTheCat();

  //A MIDDLE THINGY
  CATS DogSpike;
  DogSpike.SetMomBins(kSteps,kMin,kMax);
  AnalysisObject.SetUpCats_pp(DogSpike,"AV18","Gauss",0,0);
  DogSpike.SetAnaSource(0,BaseSource + KstarLin*1e-3*56 + KstarPar*1e-6*56*56);
  DogSpike.KillTheCat();

  //for the fit
  CATS DogTyke;
  DogTyke.SetMomBins(kSteps,kMin,kMax);
  AnalysisObject.SetUpCats_pp(DogTyke,"AV18","Gauss",0,0);
  DogTyke.SetAnaSource(0,BaseSource + KstarLin*1e-3*100 + KstarPar*1e-6*100*100);
  DogTyke.KillTheCat();

  //ATTEMPT TO DO IT RIGHT
  CATS MouseNibbles;
  MouseNibbles.SetMomBins(kSteps,kMin,kMax);
  AnalysisObject.SetUpCats_pp(MouseNibbles,"AV18","",0,0);
  CATSparameters cPars(CATSparameters::tSource, 3, true);
  MouseNibbles.SetUseAnalyticSource(true);
  MouseNibbles.SetMomentumDependentSource(true);
  MouseNibbles.SetAnaSource(GaussSourceKstarPol2, cPars);
  MouseNibbles.SetAnaSource(0,BaseSource);
  MouseNibbles.SetAnaSource(1,KstarLin);
  MouseNibbles.SetAnaSource(2,KstarPar);
  MouseNibbles.SetAutoNormSource(false);
  MouseNibbles.SetNormalizedSource(true);
  MouseNibbles.KillTheCat();

  TGraph gCkTrue;
  gCkTrue.SetName("gCkTrue");
  gCkTrue.SetLineColor(kGreen+1);
  gCkTrue.SetLineWidth(5);

  TGraph gCkJerry;
  gCkJerry.SetName("gCkJerry");
  gCkJerry.SetLineColor(kOrange+4);
  gCkJerry.SetLineWidth(5);

  TGraph gCkTom;
  gCkTom.SetName("gCkTom");
  gCkTom.SetLineColor(kBlue-4);
  gCkTom.SetLineWidth(5);

  TGraph gCkSpike;
  gCkSpike.SetName("gCkSpike");
  gCkSpike.SetLineColor(kGray);
  gCkSpike.SetLineWidth(5);

  TGraph gCkNibbles;
  gCkNibbles.SetName("gCkNibbles");
  gCkNibbles.SetLineColor(kPink);
  gCkNibbles.SetLineWidth(5);

  TH1F* hCkNibbles = new TH1F("hCkNibbles","hCkNibbles",kSteps,kMin,kMax);
  hCkNibbles->SetLineColor(kPink);
  hCkNibbles->SetLineWidth(5);

  TGraph* gSrNibbles = new TGraph [kSteps];
  TGraph gSkNibbles;
  gSkNibbles.SetName("gSkNibbles");

  for(unsigned uMom=0; uMom<kSteps; uMom++){
    //TRUE
    CATS CatMom;
    CatMom.SetMomBins(1,hCkTrue->GetXaxis()->GetBinLowEdge(uMom+1),hCkTrue->GetXaxis()->GetBinUpEdge(uMom+1));
    AnalysisObject.SetUpCats_pp(CatMom,"AV18","Gauss",0,0);
    double MOM = hCkTrue->GetBinCenter(uMom+1);
    CatMom.SetAnaSource(0,BaseSource + KstarLin*1e-3*MOM + KstarPar*1e-6*MOM*MOM);
    CatMom.KillTheCat();
    gCkTrue.SetPoint(uMom,MOM,CatMom.GetCorrFun(0));

    gCkJerry.SetPoint(uMom,MOM,MouseJerry.GetCorrFun(uMom));
    gCkTom.SetPoint(uMom,MOM,CatTom.GetCorrFun(uMom));
    gCkSpike.SetPoint(uMom,MOM,DogSpike.GetCorrFun(uMom));
    gCkNibbles.SetPoint(uMom,MOM,MouseNibbles.GetCorrFun(uMom));
    hCkNibbles->SetBinContent(uMom+1,MouseNibbles.GetCorrFun(uMom));
    double CkErr = 0.5*exp(-pow(MOM,sqrt(2))/180.)+0.01;
    hCkNibbles->SetBinError(uMom+1,CkErr);

    gSrNibbles[uMom].SetName(TString::Format("gSrNibbles_%.0f",MOM));
    double Smean = 0;
    double Sint = 0;
    unsigned NumRadPts = 0;
    double SPARS[6];
    SPARS[0] = MOM;
    SPARS[2] = 0;
    SPARS[3] = BaseSource;
    SPARS[4] = KstarLin;   
    SPARS[5] = KstarPar;   
    for(double dRad=0.05; dRad<16; dRad+=0.05){
      SPARS[1] = dRad;
      double Sr_value = GaussSourceKstarPol2(SPARS);
      Smean += Sr_value*dRad;
      Sint += Sr_value;
      gSrNibbles[uMom].SetPoint(NumRadPts, dRad, Sr_value);
      NumRadPts++;
    }
    Smean /= Sint;
    gSkNibbles.SetPoint(uMom, MOM, Smean);
  }
  DGE_Cat = &DogTyke;
  DogTyke.SetNotifications(CATS::nWarning);
  TF1* fCkNibbles = new TF1(TString::Format("fCkNibbles"),DGE_Fitter,kMin,kMax,2);
  fCkNibbles->SetParameter(0, BaseSource + KstarLin*1e-3*(kMin+kMax)*0.5 + KstarPar*1e-6*pow((kMin+kMax)*0.5,2.));
  fCkNibbles->SetParLimits(0, BaseSource, BaseSource + KstarLin*1e-3*(kMax) + KstarPar*1e-6*pow(kMax,2.));
  fCkNibbles->SetParameter(1, 1);
  fCkNibbles->SetParLimits(1, 0.8, 1.2);
  hCkNibbles->Fit(fCkNibbles, "S, N, R, M");



  TString OutputFileName = TString::Format("%s/HighMtProblem/KstarSource_Test1/JerryMouse_l%.2f_q%.2f.root",GetFemtoOutputFolder(),KstarLin,KstarPar);
  TFile fOutput(OutputFileName,"recreate");
  gCkTrue.Write();
  gCkJerry.Write();
  gCkTom.Write();
  gCkSpike.Write();
  gCkNibbles.Write();
  hCkNibbles->Write();
  fCkNibbles->Write();
  gSkNibbles.Write();

  //SLIDES FOR LAURA
  TCanvas* cCorrelations = new TCanvas("cCorrelations", "cCorrelations", 1);
  cCorrelations->cd(0);
  cCorrelations->SetCanvasSize(1280, 720);
  cCorrelations->SetMargin(0.15,0.05,0.2,0.05);//lrbt
  TH1F* hAxisPotCk = new TH1F("hAxisPotCk", "hAxisPotCk", 90, 0, 180);
  hAxisPotCk->SetStats(false);
  hAxisPotCk->SetTitle("");
  hAxisPotCk->GetXaxis()->SetTitle("k* (MeV)");
  hAxisPotCk->GetXaxis()->SetTitleSize(0.06);
  hAxisPotCk->GetXaxis()->SetLabelSize(0.06);
  hAxisPotCk->GetXaxis()->CenterTitle();
  hAxisPotCk->GetXaxis()->SetTitleOffset(1.3);
  hAxisPotCk->GetXaxis()->SetLabelOffset(0.02);
  hAxisPotCk->GetYaxis()->SetTitle("C(k*)");
  hAxisPotCk->GetYaxis()->SetTitleSize(0.06);
  hAxisPotCk->GetYaxis()->SetLabelSize(0.06);
  hAxisPotCk->GetYaxis()->CenterTitle();
  hAxisPotCk->GetYaxis()->SetTitleOffset(1.10);
  hAxisPotCk->GetYaxis()->SetRangeUser(0, 5.5);
  hAxisPotCk->Draw("axis");
  gCkNibbles.SetLineWidth(5);
  gCkNibbles.SetLineColor(kGreen+1);
  gCkNibbles.Draw("same,L");
  fCkNibbles->SetLineWidth(5);
  fCkNibbles->SetLineColor(kRed+1);
  fCkNibbles->SetLineStyle(2);
  fCkNibbles->Draw("same");
  TLegend* lCk = new TLegend(0.65,0.75,0.95,0.95);//lbrt
  lCk->SetName("lCk");
  lCk->SetTextSize(0.045);
  lCk->AddEntry(&gCkNibbles, "True correlation");
  lCk->AddEntry(fCkNibbles, "Fit with single Gauss");
  lCk->Draw("same");
  cCorrelations->Write();
  cCorrelations->SaveAs(TString::Format("%s/HighMtProblem/KstarSource_Test1/JerryMouseCk_l%.2f_q%.2f.png",GetFemtoOutputFolder(),KstarLin,KstarPar));

  TCanvas* cCorrelationsZoom = new TCanvas("cCorrelationsZoom", "cCorrelationsZoom", 1);
  cCorrelationsZoom->cd(0);
  cCorrelationsZoom->SetCanvasSize(1280, 720);
  cCorrelationsZoom->SetMargin(0.15,0.05,0.2,0.05);//lrbt
  TH1F* hAxisCkZoom = new TH1F("hAxisCkZoom", "hAxisCkZoom", 60, 60, 180);
  hAxisCkZoom->SetStats(false);
  hAxisCkZoom->SetTitle("");
  hAxisCkZoom->GetXaxis()->SetTitle("k* (MeV)");
  hAxisCkZoom->GetXaxis()->SetTitleSize(0.06);
  hAxisCkZoom->GetXaxis()->SetLabelSize(0.06);
  hAxisCkZoom->GetXaxis()->CenterTitle();
  hAxisCkZoom->GetXaxis()->SetTitleOffset(1.3);
  hAxisCkZoom->GetXaxis()->SetLabelOffset(0.02);
  hAxisCkZoom->GetYaxis()->SetTitle("C(k*)");
  hAxisCkZoom->GetYaxis()->SetTitleSize(0.06);
  hAxisCkZoom->GetYaxis()->SetLabelSize(0.06);
  hAxisCkZoom->GetYaxis()->CenterTitle();
  hAxisCkZoom->GetYaxis()->SetTitleOffset(1.10);
  hAxisCkZoom->GetYaxis()->SetRangeUser(0.85, 1.1);
  hAxisCkZoom->Draw("axis");
  gCkNibbles.Draw("same,L");
  fCkNibbles->Draw("same");
  cCorrelationsZoom->Write();
  cCorrelationsZoom->SaveAs(TString::Format("%s/HighMtProblem/KstarSource_Test1/JerryMouseCkZoom_l%.2f_q%.2f.png",GetFemtoOutputFolder(),KstarLin,KstarPar));


  TCanvas* cSk = new TCanvas("cSk", "cSk", 1);
  cSk->cd(0);
  cSk->SetCanvasSize(1280, 720);
  cSk->SetMargin(0.15,0.05,0.2,0.05);//lrbt
  TH1F* hAxisSk = new TH1F("hAxisSk", "hAxisSk", 90, 0, 180);
  hAxisSk->SetStats(false);
  hAxisSk->SetTitle("");
  hAxisSk->GetXaxis()->SetTitle("k* (MeV)");
  hAxisSk->GetXaxis()->SetTitleSize(0.06);
  hAxisSk->GetXaxis()->SetLabelSize(0.06);
  hAxisSk->GetXaxis()->CenterTitle();
  hAxisSk->GetXaxis()->SetTitleOffset(1.3);
  hAxisSk->GetXaxis()->SetLabelOffset(0.02);
  hAxisSk->GetYaxis()->SetTitle("<S(r)> fm");
  hAxisSk->GetYaxis()->SetTitleSize(0.06);
  hAxisSk->GetYaxis()->SetLabelSize(0.06);
  hAxisSk->GetYaxis()->CenterTitle();
  hAxisSk->GetYaxis()->SetTitleOffset(1.10);
  hAxisSk->GetYaxis()->SetRangeUser(1.7, 3.5);
  hAxisSk->Draw("axis");
  gSkNibbles.SetLineWidth(5);
  gSkNibbles.SetLineColor(kBlack);
  gSkNibbles.Draw("same,L");
  cSk->Write();
  cSk->SaveAs(TString::Format("%s/HighMtProblem/KstarSource_Test1/JerryMouseSk_l%.2f_q%.2f.png",GetFemtoOutputFolder(),KstarLin,KstarPar));


  for(unsigned uMom=0; uMom<kSteps; uMom++){
    gSrNibbles[uMom].Write();
  }
  delete [] gSrNibbles;

  delete hCkNibbles;
  delete fCkNibbles;
}

void CECA_kstar_vs_integrated(const float frac_proton_reso, const float hours, const int SEED=11){

  const double disp = 0.288;
  const double hT = 3.23;
  const double tau = 3.26;

  const unsigned target_yield = 1000*1000*1000;
  const unsigned multiplicity = 2;

  const double TIMEOUT = 30;
  const double GLOB_TIMEOUT = hours*3600;

  const double femto_region = 100;



  TREPNI Database(0);
  Database.SetSeed(SEED);
  std::vector<TreParticle*> ParticleList;
  ParticleList.push_back(Database.NewParticle("Proton"));
  ParticleList.push_back(Database.NewParticle("Pion"));
  ParticleList.push_back(Database.NewParticle("ProtonReso"));

  DLM_Histo<float>* dlm_pT_eta_p;
  dlm_pT_eta_p = GetPtEta_13TeV(
    TString::Format("%s/CatsFiles/Source/CECA/proton_pT/p_dist_13TeV_ClassI.root",GetCernBoxDimi()),
    "Graph1D_y1", 500, 4050, 0.8);

  //const double frac_proton_reso = 64.22;
  //const double frac_proton_reso = 100;
  const double tau_proton_reso = 1.65;
  const double m_proton_reso = 1362;

  for(TreParticle* prt : ParticleList){
    if(prt->GetName()=="Proton"){
      prt->SetMass(Mass_p);
      prt->SetAbundance(100.-frac_proton_reso);
      prt->SetRadius(0);
      prt->SetRadiusSlope(0);
      prt->SetDelayTau(0);
      if(dlm_pT_eta_p) prt->SetPtEtaPhi(*dlm_pT_eta_p);
      else prt->SetPtPz(0.85*prt->GetMass(),0.85*prt->GetMass());
    }
    else if(prt->GetName()=="ProtonReso"){
      prt->SetMass(m_proton_reso);
      prt->SetAbundance(frac_proton_reso);
      prt->SetWidth(hbarc/tau_proton_reso);
      prt->SetDelayTau(0);

      prt->NewDecay();
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Proton"));
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
      prt->GetDecay(0)->SetBranching(100);

      if(dlm_pT_eta_p) prt->SetPtEtaPhi(*dlm_pT_eta_p);
      else prt->SetPtPz(0.85*prt->GetMass(),0.85*prt->GetMass());
    }
    else if(prt->GetName()=="Pion"){
      prt->SetMass(Mass_pic);
      prt->SetAbundance(0);
      prt->SetRadius(0);
      prt->SetRadiusSlope(0);
    }
  }//ParticleList


  std::vector<std::string> ListOfParticles;
  ListOfParticles.push_back("Proton");
  ListOfParticles.push_back("Proton");

  CECA Ivana(Database,ListOfParticles);

  Ivana.SetDisplacementX(disp);
  Ivana.SetDisplacementY(disp);
  Ivana.SetDisplacementZ(disp);

  Ivana.SetHadronizationX(hT);
  Ivana.SetHadronizationY(hT);
  Ivana.SetHadronizationZ(0);
  Ivana.SetHadrFluctuation(0);

  Ivana.SetTau(tau);
  //Ivana.SetTauFluct(tau_fct);

  //Ivana.SetThermalKick(th_kick);
  //Ivana.SetFixedHadr(fixed_hdr);
  //Ivana.SetFragmentBeta(frag_beta);

  Ivana.SetTargetStatistics(target_yield);
  Ivana.SetEventMult(multiplicity);
  Ivana.SetSourceDim(2);
  Ivana.SetDebugMode(false);

  Ivana.SetThreadTimeout(TIMEOUT);
  Ivana.SetGlobalTimeout(GLOB_TIMEOUT);
  Ivana.EqualizeFsiTime(true);
  Ivana.SetFemtoRegion(femto_region);

  Ivana.GHETTO_EVENT = true;

  Ivana.Ghetto_NumMtBins = 5;
  Ivana.Ghetto_MtBins = new double [Ivana.Ghetto_NumMtBins+1];
  Ivana.Ghetto_MtBins[0] = 930;
  Ivana.Ghetto_MtBins[1] = 1080;
  Ivana.Ghetto_MtBins[2] = 1200;
  Ivana.Ghetto_MtBins[3] = 1380;
  Ivana.Ghetto_MtBins[4] = 1840;
  Ivana.Ghetto_MtBins[5] = 4500;

  Ivana.Ghetto_NumMomBins = 100;
  Ivana.Ghetto_MomMin = 0;
  Ivana.Ghetto_MomMax = 400;

  Ivana.Ghetto_NumRadBins = 1024;
  Ivana.Ghetto_RadMin = 0;
  Ivana.Ghetto_RadMax = 64;

  Ivana.GoBabyGo(24);

  Ivana.GhettoFemto_mT_rstar->ComputeError();
  TH2F* h_GhettoFemto_mT_rstar = Convert_DlmHisto_TH2F(Ivana.GhettoFemto_mT_rstar,"GhettoFemto_mT_rstar");


  TString BaseFileName = TString::Format("%s/HighMtProblem/CECA_kstar_vs_integrated/",GetFemtoOutputFolder());
  Ivana.Ghetto_kstar_rstar_mT->QuickWrite(BaseFileName+"Ghetto_kstar_rstar_mT.dlm",true);

  TFile fOutput(BaseFileName+TString::Format("fOutput_FR%.0f_S%i.root",frac_proton_reso,SEED),"recreate");
  TGraph gfm_mT;
  gfm_mT.SetName("gfm_mT");
  gfm_mT.SetLineWidth(3);
  TGraph* gfm_kstar = new TGraph [Ivana.GhettoFemto_mT_rstar->GetNbins(0)];
  for(unsigned umT=0; umT<Ivana.GhettoFemto_mT_rstar->GetNbins(0); umT++){
    TH1D* hProj = h_GhettoFemto_mT_rstar->ProjectionY(TString::Format("hProj"),umT+1,umT+1);
    double Mean = hProj->GetMean();
    delete hProj;
    gfm_mT.SetPoint(umT,Ivana.GhettoFemto_mT_rstar->GetBinCenter(0,umT),Mean);

    //  //how the radius changes as a function of kstar
    gfm_kstar[umT].SetName(TString::Format("gfm_kstar_mT%.0f",Ivana.GhettoFemto_mT_rstar->GetBinCenter(0,umT)));
    gfm_kstar[umT].SetLineColor(kRed+umT);
    gfm_kstar[umT].SetLineWidth(3);
    gfm_kstar[umT].SetMarkerStyle(2);
    gfm_kstar[umT].SetMarkerColor(kBlack);
    gfm_kstar[umT].SetMarkerSize(1);
    TH2F* h_Ghetto_kstar_rstar_mT = new TH2F(TString::Format("h_Ghetto_kstar_rstar_mT%.0f",Ivana.GhettoFemto_mT_rstar->GetBinCenter(0,umT)),
            TString::Format("h_Ghetto_kstar_rstar_mT%.0f",Ivana.GhettoFemto_mT_rstar->GetBinCenter(0,umT)),
            Ivana.Ghetto_kstar_rstar_mT->GetNbins(0),Ivana.Ghetto_kstar_rstar_mT->GetLowEdge(0),Ivana.Ghetto_kstar_rstar_mT->GetUpEdge(0),
            Ivana.Ghetto_kstar_rstar_mT->GetNbins(1),Ivana.Ghetto_kstar_rstar_mT->GetLowEdge(1),Ivana.Ghetto_kstar_rstar_mT->GetUpEdge(1));
    for(unsigned ukstar=0; ukstar<Ivana.Ghetto_kstar_rstar_mT->GetNbins(0); ukstar++){
      for(unsigned urstar=0; urstar<Ivana.Ghetto_kstar_rstar_mT->GetNbins(1); urstar++){
        h_Ghetto_kstar_rstar_mT->SetBinContent(ukstar+1,urstar+1,Ivana.Ghetto_kstar_rstar_mT->GetBinContent(ukstar,urstar,umT));
        //printf("mt%u ks%u rs%u --> %.2f\n",umT,ukstar,urstar,Ivana.Ghetto_kstar_rstar_mT->GetBinContent(ukstar,urstar,umT));
        //usleep(1e3);
      }
    }
    for(unsigned ukstar=0; ukstar<Ivana.Ghetto_kstar_rstar_mT->GetNbins(0); ukstar++){
      hProj = h_Ghetto_kstar_rstar_mT->ProjectionY(TString::Format("hProj"),ukstar+1,ukstar+1);
      Mean = hProj->GetMean();
      delete hProj;
      gfm_kstar[umT].SetPoint(ukstar,Ivana.Ghetto_kstar_rstar_mT->GetBinCenter(0,ukstar),Mean);
    }

    fOutput.cd();
    h_Ghetto_kstar_rstar_mT->Write();
    gfm_kstar[umT].Write();
    //delete h_Ghetto_kstar_rstar_mT;
  }
}


DLM_Histo<float>* dlmCecaSource1 = NULL;

void Make_Ck_from_Ceca(){
  TString BaseFileName = TString::Format("%s/HighMtProblem/CECA_kstar_vs_integrated/",GetFemtoOutputFolder());
  TString IntputFileName = BaseFileName + "fOutput_FR64.root";

  const double kUpperLimit = 240;

  const unsigned NumMtBins = 5;
  const double MtMin = 0.5;
  const double MtMax = 5.5;

  TString* IntputHistName = new TString [NumMtBins];
  IntputHistName[0] = "h_Ghetto_kstar_rstar_mT1005";
  IntputHistName[1] = "h_Ghetto_kstar_rstar_mT1140";
  IntputHistName[2] = "h_Ghetto_kstar_rstar_mT1290";
  IntputHistName[3] = "h_Ghetto_kstar_rstar_mT1610";
  IntputHistName[4] = "h_Ghetto_kstar_rstar_mT3170";

  TFile fInput(IntputFileName,"read");

  TH2F* hTemp = (TH2F*)fInput.Get(IntputHistName[0]);
  DLM_Histo<float>* dlmTemp = Convert_TH2F_DlmHisto(hTemp);
  //hTemp->Rebin2D(2,2);
  unsigned NumMomBins = dlmTemp->GetNbins(0);
  double* kRange = dlmTemp->GetBinRange(0);

  unsigned NumRadBins = dlmTemp->GetNbins(1);
  double* rRange = dlmTemp->GetBinRange(1);

  delete dlmTemp;

  for(unsigned uMom=0; uMom<NumMomBins; uMom++){
    double kstar = (kRange[uMom]+kRange[uMom+1])*0.5;
    if(kstar>kUpperLimit){
      NumMomBins = uMom;
      break;
    }
  }
  //unsigned NumMomBinsCats = NumMomBins*2;
  printf("NumMomBins,kmin,kmax = %u, %.0f, %.0f\n",NumMomBins,kRange[0],kRange[NumMomBins]);

  if(dlmCecaSource1) {delete [] dlmCecaSource1; dlmCecaSource1 = NULL;}
  dlmCecaSource1 = new DLM_Histo<float> [NumMtBins];

//we could rebin later to do our life easier
  dlmCecaSource1->SetUp(3);
  dlmCecaSource1->SetUp(0,NumMtBins,MtMin,MtMax);
  dlmCecaSource1->SetUp(1,NumMomBins,kRange);
  dlmCecaSource1->SetUp(2,NumRadBins,rRange);
  dlmCecaSource1->Initialize();

  DLM_Histo<float>* KstarIntegrSource;
  KstarIntegrSource = new DLM_Histo<float> [NumMtBins];


  for(unsigned uMt=0; uMt<NumMtBins; uMt++){
    KstarIntegrSource[uMt].SetUp(1);
    KstarIntegrSource[uMt].SetUp(0,NumRadBins,rRange);
    TH2F* histo = (TH2F*)fInput.Get(IntputHistName[uMt]);
    for(unsigned uMom=0; uMom<NumMomBins; uMom++){
      for(unsigned uRad=0; uRad<NumRadBins; uRad++){
        dlmCecaSource1->SetBinContent(uMt,uMom,uRad,histo->GetBinContent(uMom+1,uRad+1));
        dlmCecaSource1->SetBinError(uMt,uMom,uRad,histo->GetBinContent(uMom+1,uRad+1));
      }
    }
  }


  //takes a 3D histo of Mt Kstar Rstar and returns a 2D kdp histo of Mt Kstar
  DLM_Histo<KdpPars>* dlmKdpPars = Convert_3Dsource_Kdp(*dlmCecaSource1);
  DLM_MtKstar_KdpSource MyKdpSource(dlmKdpPars);

  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());
  
  CATS* Cat_pp = new CATS [NumMtBins];

  CATS** SillyCat_pp = new CATS* [NumMtBins];
  CATS* FittyCat_pp = new CATS [NumMtBins];

  DLM_HistoSource** HistoSourceKstar = new DLM_HistoSource* [NumMtBins];

  TH1F** hCat_pp = new TH1F* [NumMtBins];
  TH1F** hSillyCat_pp = new TH1F* [NumMtBins];
  TH1F*** hSrc_pp = new TH1F** [NumMtBins];

  TF1** fCkSg = new TF1* [NumMtBins];


  for(unsigned uMt=0; uMt<NumMtBins; uMt++){
    SillyCat_pp[uMt] = new CATS [NumMomBins];
    hCat_pp[uMt] = new TH1F(TString::Format("hCat_pp_uMt%u",uMt),TString::Format("hCat_pp_uMt%u",uMt),
                                NumMomBins, kRange);
    hSillyCat_pp[uMt] = new TH1F(TString::Format("hSillyCat_pp_uMt%u",uMt),TString::Format("hSillyCat_pp_pp_uMt%u",uMt),
                                NumMomBins, kRange);

    
    printf("uMt = %u\n",uMt);
    hSrc_pp[uMt] = new TH1F* [NumMomBins];
    Cat_pp[uMt].SetMomBins(NumMomBins,kRange);
    AnalysisObject.SetUpCats_pp(Cat_pp[uMt],"AV18","",0,0);
    Cat_pp[uMt].SetAnaSource(CatsSourceForwarder, &MyKdpSource, 1);
    Cat_pp[uMt].SetAnaSource(0, dlmCecaSource1->GetBinCenter(0, uMt));
    Cat_pp[uMt].SetNormalizedSource(false);
    Cat_pp[uMt].SetAutoNormSource(false);
    Cat_pp[uMt].SetMomentumDependentSource(true);
    Cat_pp[uMt].SetNotifications(CATS::nWarning);
    Cat_pp[uMt].KillTheCat();
    //Cat_pp[uMt].SetName(TString::Format("gCat_pp_uMt%u",uMt));
    double src_pars[4];
    src_pars[3] = dlmCecaSource1->GetBinCenter(0,uMt);
    for(unsigned uMom=0; uMom<NumMomBins; uMom++){
      double Kstar = Cat_pp[uMt].GetMomentum(uMom);
      
      hCat_pp[uMt]->SetBinContent(uMom+1,Kstar,Cat_pp[uMt].GetCorrFun(uMom));
      //0.3∙exp(−(x^1.25)/180)+0.01
      double CkErr = 0.5*exp(-pow(Kstar,sqrt(2))/180.)+0.01;
      //double CkErr = 0.001*(1+60*exp(-pow(Cat_pp[uMt].GetMomentum(uMom)/75.,2.)));
      hCat_pp[uMt]->SetBinError(uMom+1,CkErr);

      TString hsname = TString::Format("hSrc_pp_uMt%u_k%.0f",uMt,Kstar);
      hSrc_pp[uMt][uMom] = new TH1F(hsname,hsname,1024,0,16);   
      src_pars[0] = Kstar;
      for(unsigned uRad=0; uRad<hSrc_pp[uMt][uMom]->GetNbinsX(); uRad++){
        src_pars[1] = hSrc_pp[uMt][uMom]->GetBinCenter(uRad+1);
        double val = MyKdpSource.Eval(src_pars);
        //printf("up %u %u %u = %.3e\n",uMt,uMom,uRad,val);
        //usleep(100e3);
        hSrc_pp[uMt][uMom]->SetBinContent(uRad+1,val);
        //usleep(10e3);
      }   

      SillyCat_pp[uMt][uMom].SetMomBins(1,dlmCecaSource1->GetBinLowEdge(1,uMom),dlmCecaSource1->GetBinUpEdge(1,uMom));
      AnalysisObject.SetUpCats_pp(SillyCat_pp[uMt][uMom],"AV18","",0,0);
      DLM_Histo<float> CurrentSource;
      CurrentSource.SetUp(1);
      CurrentSource.SetUp(0,NumRadBins,rRange);
      CurrentSource.Initialize();
      for(unsigned uRad=0; uRad<NumRadBins;uRad++){
        src_pars[1] = CurrentSource.GetBinCenter(0,uRad);
        double val = MyKdpSource.Eval(src_pars);
        CurrentSource.SetBinContent(uRad,val);
      }
      DLM_HistoSource HistoSource(CurrentSource);
      SillyCat_pp[uMt][uMom].SetAnaSource(CatsSourceForwarder, &HistoSource, 0);
      SillyCat_pp[uMt][uMom].SetNotifications(CATS::nWarning);
      SillyCat_pp[uMt][uMom].KillTheCat();

      hSillyCat_pp[uMt]->SetBinContent(uMom+1, SillyCat_pp[uMt][uMom].GetCorrFun(0));

    }
    double AvgSrc_at100 = hSrc_pp[uMt][hCat_pp[uMt]->FindBin(100)-1]->GetMean();
    fCkSg[uMt] = new TF1(TString::Format("fCkSg_uMt%u",uMt),DGE_Fitter,0,kUpperLimit,2);
    fCkSg[uMt]->SetParameter(0, AvgSrc_at100/2.3);
    fCkSg[uMt]->SetParLimits(0, AvgSrc_at100/4.6, AvgSrc_at100/1.15);
    fCkSg[uMt]->SetParameter(1, 1);
    fCkSg[uMt]->SetParLimits(1, 0.8, 1.2);    


    FittyCat_pp[uMt].SetMomBins(NumMomBins, kRange);
    AnalysisObject.SetUpCats_pp(FittyCat_pp[uMt],"AV18","Gauss",0,0);
    FittyCat_pp[uMt].SetAnaSource(0, AvgSrc_at100/2.3);
    FittyCat_pp[uMt].SetNotifications(CATS::nWarning);
    FittyCat_pp[uMt].KillTheCat();
    DGE_Cat = &FittyCat_pp[uMt];

    hCat_pp[uMt]->Fit(fCkSg[uMt], "S, N, R, M");

  }

  TFile fOutput(BaseFileName+"Make_Ck_from_Ceca.root","recreate");
  for(unsigned uMt=0; uMt<NumMtBins; uMt++){
    DGE_Cat = &FittyCat_pp[uMt];
    hCat_pp[uMt]->Write();
    hSillyCat_pp[uMt]->Write();
    fCkSg[uMt]->Write();
    delete hCat_pp[uMt];
    delete hSillyCat_pp[uMt];
    delete fCkSg[uMt];
  }

  for(unsigned uMt=0; uMt<NumMtBins; uMt++){
    for(unsigned uMom=0; uMom<NumMomBins; uMom++){
      hSrc_pp[uMt][uMom]->Write();
      delete hSrc_pp[uMt][uMom];
    }
    delete [] hSrc_pp[uMt];
  }
  delete [] hSrc_pp;
  
  delete [] IntputHistName;
  delete [] dlmCecaSource1; dlmCecaSource1 = NULL;
  delete [] Cat_pp; Cat_pp = NULL;
  delete [] HistoSourceKstar; HistoSourceKstar = NULL;
  delete [] hCat_pp; hCat_pp = NULL;
  delete [] hSillyCat_pp; hSillyCat_pp = NULL;
  delete [] fCkSg; fCkSg = NULL;
  delete [] kRange;
  delete [] rRange;
  delete [] FittyCat_pp;
  
}


double KDDG_v1(double* pars){
  double& kstar = pars[0];
  double& rad = pars[1];
  //base size of the first gauss
  double& SZ1 = pars[3];//[0]
  //linear slope of the first gauss
  double& LSL1 = pars[4];//[1]
  //quadratic slope of the first gauss
  double& QSL1 = pars[5];//[2]
  //base size of the second gauss
  double& SZ2 = pars[6];//[3]
  //linear slope of the second gauss
  double& LSL2 = pars[7];//[4]
  //quadratic slope of the second gauss
  double& QSL2 = pars[8];//[5]

  //min weight of the first Gauss
  double& minW1 = pars[9];//[6]
  //max weight of the first Gauss
  double& maxW1 = pars[10];//[7]
  //the position where the two source are most balanced
  double& midway = pars[11];//[8]
  //the slope of the transition
  double& slope = pars[12];//[9]


  double Size1 = SZ1 + LSL1*kstar*1e-3 + QSL1*kstar*kstar*1e-6;
  double Size2 = SZ2 + LSL2*kstar*1e-3 + QSL2*kstar*kstar*1e-6;

  double Gauss1 = 4.*Pi*rad*rad*pow(4.*Pi*Size1*Size1,-1.5)*exp(-(rad*rad)/(4.*Size1*Size1));
  double Gauss2 = 4.*Pi*rad*rad*pow(4.*Pi*Size2*Size2,-1.5)*exp(-(rad*rad)/(4.*Size2*Size2));

  double GentleStep = minW1 + (maxW1-minW1)/(1.+exp((kstar-midway)/slope));

  return GentleStep*Gauss1 + (1.-GentleStep)*Gauss2;
}

//here, we will use a custom source KDDG_v1
void KstarDepDoublGauss_pp_pL_pipi(){

  unsigned NumMomBins = 50;
  double kMin = 0;
  double kMax = 200;

  double BaseCoreSize = 0.8;//0.8
  double GaussSize = 1.5;//1.3
  double LinearKstar = 0;//0
  double QuadraticKstar = 0;//6

  double MinW1 = 0;//0
  double MaxW1 = 1;//1
  double MidWay = 36;//36
  double Slope = 12;//12

  //double BaseHaloSize = 1.4;

  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());

  CATS Kitty_pp;
  Kitty_pp.SetMomBins(NumMomBins,kMin,kMax);
  AnalysisObject.SetUpCats_pp(Kitty_pp,"AV18","Gauss",0,0);


  CATS Kitty_pL;
  Kitty_pL.SetMomBins(NumMomBins,kMin,kMax);
  AnalysisObject.SetUpCats_pL(Kitty_pL,"Usmani","Gauss",0,0);

  CATS Kitty_pipi;
  Kitty_pipi.SetMomBins(NumMomBins,kMin,kMax);
  AnalysisObject.SetUpCats_pipi(Kitty_pipi,"Gauss",0);


  Kitty_pp.SetAnaSource(0,GaussSize);
  Kitty_pL.SetAnaSource(0,GaussSize);
  Kitty_pipi.SetAnaSource(0,GaussSize);

  Kitty_pp.SetMomentumDependentSource(true);
  Kitty_pL.SetMomentumDependentSource(true);
  Kitty_pipi.SetMomentumDependentSource(true);
  //Kitty_pp.SetAutoNormSource(false);
  //Kitty_pp.SetNormalizedSource(true);

  Kitty_pp.KillTheCat();
  Kitty_pL.KillTheCat();
  Kitty_pipi.KillTheCat();

  TGraph gCk_pp_G;
  gCk_pp_G.SetName("gCk_pp_G");
  TGraph gCk_pL_G;
  gCk_pL_G.SetName("gCk_pL_G");
  TGraph gCk_pipi_G;
  gCk_pipi_G.SetName("gCk_pipi_G");

  TGraph gCk_pp_KDDG;
  gCk_pp_KDDG.SetName("gCk_pp_KDDG");
  TGraph gCk_pL_KDDG;
  gCk_pL_KDDG.SetName("gCk_pL_KDDG");
  TGraph gCk_pipi_KDDG;
  gCk_pipi_KDDG.SetName("gCk_pipi_KDDG");

  for(unsigned uBin=0; uBin<NumMomBins; uBin++){
    gCk_pp_G.SetPoint(uBin, Kitty_pp.GetMomentum(uBin), Kitty_pp.GetCorrFun(uBin));
    gCk_pL_G.SetPoint(uBin, Kitty_pL.GetMomentum(uBin), Kitty_pL.GetCorrFun(uBin));
    gCk_pipi_G.SetPoint(uBin, Kitty_pipi.GetMomentum(uBin), Kitty_pipi.GetCorrFun(uBin));
  }

  CATSparameters src_pp(CATSparameters::tSource,10,true);
  CATSparameters src_pL(CATSparameters::tSource,10,true);
  CATSparameters src_pipi(CATSparameters::tSource,10,true);
  
  Kitty_pp.SetAnaSource(KDDG_v1, src_pp);
  Kitty_pL.SetAnaSource(KDDG_v1, src_pL);
  Kitty_pipi.SetAnaSource(KDDG_v1, src_pipi);

  Kitty_pp.SetAnaSource(0, BaseCoreSize);
  Kitty_pp.SetAnaSource(1, 0);
  Kitty_pp.SetAnaSource(2, 0);
  Kitty_pp.SetAnaSource(3, GaussSize);
  Kitty_pp.SetAnaSource(4, LinearKstar);
  Kitty_pp.SetAnaSource(5, QuadraticKstar);
  Kitty_pp.SetAnaSource(6, MinW1);
  Kitty_pp.SetAnaSource(7, MaxW1);
  Kitty_pp.SetAnaSource(8, MidWay);
  Kitty_pp.SetAnaSource(9, Slope);

  Kitty_pL.SetAnaSource(0, BaseCoreSize);
  Kitty_pL.SetAnaSource(1, 0);
  Kitty_pL.SetAnaSource(2, 0);
  Kitty_pL.SetAnaSource(3, GaussSize);
  Kitty_pL.SetAnaSource(4, LinearKstar);
  Kitty_pL.SetAnaSource(5, QuadraticKstar);
  Kitty_pL.SetAnaSource(6, MinW1);
  Kitty_pL.SetAnaSource(7, MaxW1);
  Kitty_pL.SetAnaSource(8, MidWay);
  Kitty_pL.SetAnaSource(9, Slope);

  Kitty_pipi.SetAnaSource(0, BaseCoreSize);
  Kitty_pipi.SetAnaSource(1, 0);
  Kitty_pipi.SetAnaSource(2, 0);
  Kitty_pipi.SetAnaSource(3, GaussSize);
  Kitty_pipi.SetAnaSource(4, LinearKstar);
  Kitty_pipi.SetAnaSource(5, QuadraticKstar);
  Kitty_pipi.SetAnaSource(6, MinW1);
  Kitty_pipi.SetAnaSource(7, MaxW1);
  Kitty_pipi.SetAnaSource(8, MidWay);
  Kitty_pipi.SetAnaSource(9, Slope);

  Kitty_pp.KillTheCat();
  Kitty_pL.KillTheCat();
  Kitty_pipi.KillTheCat();

  for(unsigned uBin=0; uBin<NumMomBins; uBin++){
    gCk_pp_KDDG.SetPoint(uBin, Kitty_pp.GetMomentum(uBin), Kitty_pp.GetCorrFun(uBin));
    gCk_pL_KDDG.SetPoint(uBin, Kitty_pL.GetMomentum(uBin), Kitty_pL.GetCorrFun(uBin));
    gCk_pipi_KDDG.SetPoint(uBin, Kitty_pipi.GetMomentum(uBin), Kitty_pipi.GetCorrFun(uBin));
  }

  TGraph gMeanSrk;
  TGraph gMeanLow;
  TGraph gMeanUp;
  gMeanSrk.SetName("gMeanSrk");
  gMeanLow.SetName("gMeanLow");
  gMeanUp.SetName("gMeanUp");
  for(unsigned uBin=0; uBin<NumMomBins; uBin++){
    double kstar = Kitty_pp.GetMomentum(uBin);

    unsigned NumRadBins = 2048;
    double RadMax = 64;
    double src_prs[13];
    src_prs[0] = kstar;
    src_prs[3] = BaseCoreSize;
    src_prs[4] = 0;
    src_prs[5] = 0;
    src_prs[6] = GaussSize;
    src_prs[7] = LinearKstar;
    src_prs[8] = QuadraticKstar;
    src_prs[11] = MidWay;
    src_prs[12] = Slope;

    TH1F* hRad = new TH1F("hRad","hRad",NumRadBins,0,RadMax);
  
    src_prs[9] = 1;
    src_prs[10] = 1;
    for(unsigned uRad=0; uRad<NumRadBins; uRad++){
      src_prs[1] = hRad->GetBinCenter(uRad+1);
      hRad->SetBinContent(uRad+1, KDDG_v1(src_prs));
    }
    gMeanLow.SetPoint(uBin, kstar, hRad->GetMean());

    src_prs[9] = 0;
    src_prs[10] = 0;
    for(unsigned uRad=0; uRad<NumRadBins; uRad++){
      src_prs[1] = hRad->GetBinCenter(uRad+1);
      hRad->SetBinContent(uRad+1, KDDG_v1(src_prs));
    }
    gMeanUp.SetPoint(uBin, kstar, hRad->GetMean());

    src_prs[9] = MinW1;
    src_prs[10] = MaxW1;
    for(unsigned uRad=0; uRad<NumRadBins; uRad++){
      src_prs[1] = hRad->GetBinCenter(uRad+1);
      hRad->SetBinContent(uRad+1, KDDG_v1(src_prs));
    }
    gMeanSrk.SetPoint(uBin, kstar, hRad->GetMean()); 

    delete hRad;
  }

  TString OutputFileName = TString::Format("%s/HighMtProblem/KstarSource_Test1/JuliNevestuli.root",GetFemtoOutputFolder());
  TFile fOutput(OutputFileName, "recreate");
  gCk_pp_G.Write();
  gCk_pp_KDDG.Write();
  gCk_pL_G.Write();
  gCk_pL_KDDG.Write();
  gCk_pipi_G.Write();
  gCk_pipi_KDDG.Write();
  gMeanSrk.Write();
  gMeanLow.Write();
  gMeanUp.Write();

}



int HIGHMTPROBLEM(int argc, char *argv[]){
    printf("HIGHMTPROBLEM, HELLO!!!\n\n");
    //Compare_pp_AV18_Haide_SingleGauss(atoi(argv[1]));
    //Effects();//double sources
    //PlayWithProtons();

    //KstarSource_Test1();
    //fraction_reso (64.22), hours, seed
    //CECA_kstar_vs_integrated(atof(argv[1]),atof(argv[2]),atoi(argv[3]));
    //Make_Ck_from_Ceca();
    //KstarDepDoublGauss_pp_pL_pipi();

    return 0;
}