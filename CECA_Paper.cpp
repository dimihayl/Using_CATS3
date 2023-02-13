
#include "CECA_Paper.h"
#include "CommonAnaFunctions.h"
#include "EnvVars.h"
#include "DLM_CppTools.h"
#include "DLM_OmpTools.h"
#include "DLM_Random.h"
#include "DLM_RootWrapper.h"
#include "DLM_HistoAnalysis.h"
#include "DLM_Source.h"
#include "DLM_Potentials.h"
#include "DLM_Histo.h"
#include "DLM_MathFunctions.h"
#include "DLM_RootFit.h"

#include <iostream>
#include <unistd.h>
#include <omp.h>

#include "TREPNI.h"
#include "CATS.h"
#include "CATSconstants.h"
#include "CATStools.h"
#include "CECA.h"


#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TF1.h"
#include "TGenPhaseSpace.h"
#include "TRandom3.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TNtuple.h"
#include "TSystem.h"
#include "TROOT.h"

#include <boost/algorithm/string.hpp>

#include<fstream>

using namespace std;

void TestSaveStuctToFile(){
  DoubleLevy SourcePars;
  SourcePars.alpha1 = 1.5;
  SourcePars.sigma1 = 1.4;
  SourcePars.alpha2 = 2.0;
  SourcePars.sigma2 = 1.3;
  SourcePars.wght1 = 0.4;

  SourcePars.Print();

  printf("Size of a float: %lu\n",sizeof(float));
  printf("Size of the stucture: %lu\n",sizeof(DoubleLevy));

  TString OutputFileName = TString::Format("%s/CECA_Paper/TestSaveStuctToFile/DoubleLevy.bin",GetFemtoOutputFolder());
  ofstream wf(OutputFileName.Data(), ios::out | ios::binary);
  if(!wf) {
     cout << "Cannot open file!" << endl;
     return;
  }
  else{
    cout << "File opened!" << endl;
  }

  wf.write((char *) &SourcePars, sizeof(DoubleLevy));

  wf.close();
  if(!wf.good()) {
    cout << "Error occurred at writing time!" << endl;
    return;
  }
  else{
    cout << "Writing to file done!" << endl;
  }


  DoubleLevy SourceParsFromFile;
  ifstream is;
  is.open (OutputFileName.Data(), ios::binary );
  is.read ((char*)&SourceParsFromFile,20);
  is.close();

  SourceParsFromFile.Print();
}

void TestDoubleSourceOperation(){
  DoubleLevy SourcePars1;
  SourcePars1.alpha1 = 1.5;
  SourcePars1.sigma1 = 1.4;
  SourcePars1.alpha2 = 2.0;
  SourcePars1.sigma2 = 1.3;
  SourcePars1.wght1 = 0.4;

  DoubleLevy SourcePars2;
  SourcePars2.alpha1 = 1.7;
  SourcePars2.sigma1 = 1.2;
  SourcePars2.alpha2 = 1.8;
  SourcePars2.sigma2 = 1.5;
  SourcePars2.wght1 = 0.6;

  printf("SourcePars1 ----------------------------\n");
  SourcePars1.Print();
  printf("SourcePars2 ----------------------------\n");
  SourcePars2.Print();

  DLM_Histo<DoubleLevy> SHst;
  SHst.SetUp(1);
  SHst.SetUp(0,2,0,2);
  SHst.Initialize();

  SHst.SetBinContent(unsigned(0),SourcePars1);
  SHst.SetBinContent(unsigned(1),SourcePars2);

  DoubleLevy s1;
  DoubleLevy s2;
  DoubleLevy a12;

  s1 = SHst.Eval(0.5);
  s2 = SHst.Eval(1.5);
  a12 = SHst.Eval(1.0);

  printf("s1 ----------------------------\n");
  s1.Print();
  printf("s2 ----------------------------\n");
  s2.Print();
  printf("a12 ----------------------------\n");
  a12.Print();

  TString HistoFileName = TString::Format("%s/CECA_Paper/TestSaveStuctToFile/Histo1.dlm.hst",GetFemtoOutputFolder());
  SHst.QuickWrite(HistoFileName.Data(),true);

  DLM_Histo<DoubleLevy> SHst_FromFile;
  SHst_FromFile.QuickLoad(HistoFileName.Data());

  DoubleLevy fs1;
  DoubleLevy fs2;
  DoubleLevy fa12;

  fs1 = SHst_FromFile.Eval(0.5);
  fs2 = SHst_FromFile.Eval(1.5);
  fa12 = SHst_FromFile.Eval(1.0);

  printf("fs1 ----------------------------\n");
  fs1.Print();
  printf("fs2 ----------------------------\n");
  fs2.Print();
  printf("fa12 ----------------------------\n");
  fa12.Print();
}

//test the DLM_Histo and statistics for pp
//fit the distos with Levy
void Test_pp_Statistics_1(){
  //12h of running on 6 core old AMD Desktop
  TString InputHistoFile = TString::Format("%s/FunWithCeca/Ceca_pp_EffFix/12h/KstarDist_pp_ET1_PR1_DD0.0_EF-1.Ghetto_kstar_rstar_mT",GetFemtoOutputFolder());
  DLM_Histo<float> kstar_rstar_mT;
  kstar_rstar_mT.QuickLoad(InputHistoFile.Data());

  const unsigned NumMomBins = kstar_rstar_mT.GetNbins(0);
  const unsigned NumRadBins = kstar_rstar_mT.GetNbins(1);
  const unsigned NumMtBins = kstar_rstar_mT.GetNbins(2);

  const double MaxKstar = 200;

  printf("NumMomBins = %u\n",NumMomBins);
  printf("NumRadBins = %u\n",NumRadBins);
  printf("NumMtBins = %u\n",NumMtBins);

  TH1F*** hMtKstar_Rad = new TH1F** [NumMtBins];
  TF1*** fDoubleGauss = new TF1** [NumMtBins];
  //TF1*** fStupidGauss = new TF1** [NumMtBins];

  for(unsigned uMt=0; uMt<NumMtBins; uMt++){
    printf("uMt = %u\n",uMt);
    double Mt = kstar_rstar_mT.GetBinCenter(2,uMt);
    hMtKstar_Rad[uMt] = new TH1F* [NumMomBins];
    fDoubleGauss[uMt] = new TF1* [NumMomBins];
    //fStupidGauss[uMt] = new TF1* [NumMomBins];
    for(unsigned uMom=0; uMom<NumMomBins; uMom++){
      double Momentum = kstar_rstar_mT.GetBinCenter(0,uMom);
      if(Momentum>MaxKstar){
        hMtKstar_Rad[uMt][uMom] = NULL;
        fDoubleGauss[uMt][uMom] = NULL;
        //fStupidGauss[uMt][uMom] = NULL;
        continue;
      }
      TString HistoName = TString::Format("hMtKstar_Rad_%.0f_%.0f",Mt,Momentum);
      hMtKstar_Rad[uMt][uMom] = new TH1F(HistoName,HistoName,
        kstar_rstar_mT.GetNbins(1), kstar_rstar_mT.GetLowEdge(1), kstar_rstar_mT.GetUpEdge(1));
      double Rad,upperlimit,lowerlimit;
      for(unsigned uRad=0; uRad<NumRadBins; uRad++){
        Rad = kstar_rstar_mT.GetBinContent(uMom,uRad,uMt);
        hMtKstar_Rad[uMt][uMom]->SetBinContent(uRad+1,Rad);
      }
      hMtKstar_Rad[uMt][uMom]->Sumw2();
      hMtKstar_Rad[uMt][uMom]->Scale(1./hMtKstar_Rad[uMt][uMom]->Integral(), "width");

      GetCentralInterval(*hMtKstar_Rad[uMt][uMom], 0.98, lowerlimit, upperlimit, true);
      if(lowerlimit>5) lowerlimit = 5;
      if(upperlimit>20) upperlimit = 20;

      TString FitName = TString::Format("fDoubleGauss_%.0f_%.0f",Mt,Momentum);
      fDoubleGauss[uMt][uMom] = new TF1(FitName,NormTripleShiftedGaussTF1,0.,20.,9);

      //like chi2, but not normalized to error
      double Dist2 = 0;
      double NDPts = 0;
      double Chi2 = 0;
      double NDPts_chi2 = 0;
      double Dist_Max = 0;
      double Nsig_AtDistMax = 0;
      double Rad_AtDistMax = 0;

      const double Dist_Limit = 0.015;
      const double Nsig_Limit = 3.0;
      const unsigned Patience = 32;//increase limit after X fits
      const double BadFitWarning = 0.06;
      double Dist_CurLim = Dist_Limit;
      unsigned StuckCount = 0;
      unsigned ResetCount = 0;

      TRandom3 rangen(11);

      do{
        Dist2 = 0;
        NDPts = 0;
        Chi2 = 0;
        NDPts_chi2 = 0;
        Dist_Max = 0;
        Nsig_AtDistMax = 0;

        //NORM
        fDoubleGauss[uMt][uMom]->SetParameter(0,rangen.Uniform(0.9,1.0));
        fDoubleGauss[uMt][uMom]->SetParLimits(0,0.,1.);

        //G1
        fDoubleGauss[uMt][uMom]->SetParameter(1,rangen.Uniform(0.2,0.4));//sigma
        fDoubleGauss[uMt][uMom]->SetParLimits(1,0.,0.8);
        fDoubleGauss[uMt][uMom]->SetParameter(2,rangen.Uniform(0.0,0.5));//shift
        fDoubleGauss[uMt][uMom]->SetParLimits(2,0.,2.0);
        fDoubleGauss[uMt][uMom]->SetParameter(3,rangen.Uniform(0.2,0.4));//weight
        fDoubleGauss[uMt][uMom]->SetParLimits(3,0.,0.95);

        if(hMtKstar_Rad[uMt][uMom]->GetEntries()>200 || ResetCount>=1){
          //G2
          //sigma
          fDoubleGauss[uMt][uMom]->SetParameter(4,rangen.Uniform(2.0,4.0));
          fDoubleGauss[uMt][uMom]->SetParLimits(4,0,20.);
          //shift
          fDoubleGauss[uMt][uMom]->SetParameter(5,rangen.Uniform(0,0.5));
          fDoubleGauss[uMt][uMom]->SetParLimits(5,0.,10.);
          //weight
          fDoubleGauss[uMt][uMom]->SetParameter(6,rangen.Uniform(0.4,0.6));
          fDoubleGauss[uMt][uMom]->SetParLimits(6,0.,1.0);
        }
        else{
          //G2
          //sigma
          fDoubleGauss[uMt][uMom]->FixParameter(4,1);
          //shift
          fDoubleGauss[uMt][uMom]->FixParameter(5,0);
          //weight
          fDoubleGauss[uMt][uMom]->FixParameter(6,0);
        }


        //if we have enough data
        if(hMtKstar_Rad[uMt][uMom]->GetEntries()>1000 || ResetCount>=2){
          //G3
          //sigma
          fDoubleGauss[uMt][uMom]->SetParameter(7,rangen.Uniform(0.4,0.7));
          fDoubleGauss[uMt][uMom]->SetParLimits(7,0,20.);
          //shift
          fDoubleGauss[uMt][uMom]->SetParameter(8,rangen.Uniform(0,0.5));
          fDoubleGauss[uMt][uMom]->SetParLimits(8,0.,10.);
        }
        else{
          //G3
          //sigma
          fDoubleGauss[uMt][uMom]->FixParameter(7,1);
          //shift
          fDoubleGauss[uMt][uMom]->FixParameter(8,0);
        }



        hMtKstar_Rad[uMt][uMom]->Fit(fDoubleGauss[uMt][uMom],"Q, S, N, R, M","",lowerlimit,upperlimit);


        //up to 10 fm
        const double RadDistMax = 8;
        for(unsigned uRad=0; uRad<NumRadBins; uRad++){
          double Rad = hMtKstar_Rad[uMt][uMom]->GetBinCenter(uRad+1);
          if(Rad>RadDistMax) break;
          double dst = hMtKstar_Rad[uMt][uMom]->GetBinContent(uRad+1)-fDoubleGauss[uMt][uMom]->Eval(Rad);
          double err;
          if(hMtKstar_Rad[uMt][uMom]->GetBinContent(uRad+1)){
            err = hMtKstar_Rad[uMt][uMom]->GetBinError(uRad+1);
          }
          else{
            err = fabs(dst)*1000;
          }

          Dist2 += dst*dst;
          NDPts++;
          if(hMtKstar_Rad[uMt][uMom]->GetBinError(uRad+1)){
            Chi2 += (dst*dst)/(err*err);
            NDPts_chi2++;
          }

          if(Dist_Max<fabs(dst)){
            Dist_Max = fabs(dst);
            Nsig_AtDistMax = Dist_Max/err;
            Rad_AtDistMax = Rad;
          }
        }
        Dist2 /= NDPts;
        Chi2 /= NDPts_chi2;
        //printf("uMt_%u uMom_%u Dist2=%e (%e); Dist_Max = %f (%.2f)\n",uMt,uMom,Dist2,Chi2,Dist_Max,Nsig_AtDistMax);
        //usleep(100e3);
        StuckCount++;
        if(StuckCount>Patience){
          Dist_CurLim += Dist_Limit;
          StuckCount = 0;
          ResetCount++;
        }
      }
      while(Dist_Max>Dist_CurLim && Nsig_AtDistMax>Nsig_Limit);
      if(Dist_Max>BadFitWarning && Nsig_AtDistMax>Nsig_Limit){
        printf("WARNING: uMt_%u uMom_%u (r = %.3f)\n",uMt,uMom,Rad_AtDistMax);
      }


/*
      TString FitName = TString::Format("fStupidGauss_%.0f_%.0f",Mt,Momentum);
      fStupidGauss[uMt][uMom] = new TF1(FitName,StupidGaussSumTF1,0.,20.,13);

      fStupidGauss[uMt][uMom]->FixParameter(0,4);

      for(unsigned uG=0; uG<4; uG++){
        fStupidGauss[uMt][uMom]->SetParameter(1+uG*3,0.5);//W
        fStupidGauss[uMt][uMom]->SetParLimits(1+uG*3,0.,1.);
        fStupidGauss[uMt][uMom]->SetParameter(2+uG*3,1.0*pow(1.+double(uG),1.5));//M
        fStupidGauss[uMt][uMom]->SetParLimits(2+uG*3,0.,4.*pow(1.+double(uG),1.5));
        fStupidGauss[uMt][uMom]->SetParameter(3+uG*3,1.0*pow(1.+double(uG),1.5));//S
        fStupidGauss[uMt][uMom]->SetParLimits(3+uG*3,0.,2.*pow(1.+double(uG),1.5));

        //fStupidGauss[uMt][uMom]->FixParameter(1+uG*3,0.5);
        //fStupidGauss[uMt][uMom]->FixParameter(2+uG*3,1.0*pow(1.+double(uG),1.5));//M
        //fStupidGauss[uMt][uMom]->FixParameter(3+uG*3,1.0*pow(1.+double(uG),1.5));//S
        //fStupidGauss[uMt][uMom]->FixParameter(1+uG*3,0.5);
        //fStupidGauss[uMt][uMom]->FixParameter(2+uG*3,1);//M
        //fStupidGauss[uMt][uMom]->FixParameter(3+uG*3,3);//S

      }
*/


    }
  }

  TFile fOutput(TString::Format("%s/CECA_Paper/Test_pp_Statistics_1/HistoAndFits.root",GetFemtoOutputFolder()),"recreate");

  for(unsigned uMt=0; uMt<NumMtBins; uMt++){
    for(unsigned uMom=0; uMom<NumMomBins; uMom++){
      if(hMtKstar_Rad[uMt][uMom]){
        hMtKstar_Rad[uMt][uMom]->Write();
        delete hMtKstar_Rad[uMt][uMom];
        hMtKstar_Rad[uMt][uMom]=NULL;
      }

      if(fDoubleGauss[uMt][uMom]){
        fDoubleGauss[uMt][uMom]->Write();
        delete fDoubleGauss[uMt][uMom];
        fDoubleGauss[uMt][uMom]=NULL;
      }
      //if(fStupidGauss[uMt][uMom]){
      //  fStupidGauss[uMt][uMom]->Write();
      //  delete fStupidGauss[uMt][uMom];
      //  fStupidGauss[uMt][uMom]=NULL;
      //}
    }
    if(hMtKstar_Rad[uMt]){delete [] hMtKstar_Rad[uMt]; hMtKstar_Rad[uMt]=NULL;}
    if(fDoubleGauss[uMt]){delete [] fDoubleGauss[uMt]; fDoubleGauss[uMt]=NULL;}
    //if(fStupidGauss[uMt]){delete [] fStupidGauss[uMt]; fStupidGauss[uMt]=NULL;}
  }
  delete [] hMtKstar_Rad;
  delete [] fDoubleGauss;
  //delete [] fStupidGauss;
}



//parameters to control:
//MUST:
//  GLOB_TIMEOUT, multiplicity, target_yield, femto_region,
//  d_x, d_y, d_z, h_x, h_y, h_z,tau;
//  hdr_size, hdr_slope, th_kick, fixed_hdr, frag_beta;
//  a flag for the momentum distribution;
//  a flag for resonance variations;
//  18 numbers in total
//
// OUTPUT:
//  DLM_Histo: Ghetto_kstar_rstar_mT
//  A settings file, which will again be a histo (or a custum class?),
//  with the values of all parameters above.
//  the names of the file should be unique, e.g. SEED id or whatever

//QUESTIONS: do we care about h_xyz differentially
// -> it leads to 1.5% error -> take as systematics :D
//      actually, if a take the hz = 2x hT, than we have +4.3% on the radius, but the scaling remains the same
//      i.e. we can comment that the choice of hz can simply scale the thing up and down by a small amount

void Ceca_pp_or_pL(const TString type, const double GLOB_TIMEOUT, const double target_yield, const double femto_region,
                    const double d_x, const double d_y, const double d_z,
                    const double h_x, const double h_y, const double h_z, const double h_fct,
                    const double tau, const double tau_fct, const double tau_prp){

  if(type!="pp"&&type!="pL"){
    printf("WHAT IS THIS: %s\n", type.Data());
    return;
  }

  const double TIMEOUT = 30;
  const double EtaCut = 0.8;
  const bool PROTON_RESO = true;
  const bool EQUALIZE_TAU = true;



  TREPNI Database(0);
  Database.SetSeed(11);
  std::vector<TreParticle*> ParticleList;
  ParticleList.push_back(Database.NewParticle("Proton"));
  ParticleList.push_back(Database.NewParticle("ProtonReso"));
  ParticleList.push_back(Database.NewParticle("Lambda"));
  ParticleList.push_back(Database.NewParticle("LambdaReso"));
  ParticleList.push_back(Database.NewParticle("Pion"));
/*
  TString OutputFolderName;
  if(EffFix) OutputFolderName = "FunWithCeca/Ceca_"+type+"_EffFix";
  else OutputFolderName = "FunWithCeca/Ceca_"+type+"_CoreFix";
  //OutputFolderName = "FunWithCeca";
  TString BaseFileName = TString::Format("%s/%s/KstarDist_%s_ET%i_PR%i_DD%.1f_EF%i",
  GetFemtoOutputFolder(),OutputFolderName.Data(),type.Data(),EQUALIZE_TAU,PROTON_RESO,d_delay,EffFix);
  TFile fOutput(BaseFileName+".root","recreate");

  DLM_Histo<float>* dlm_pT_p = NULL;
  DLM_Histo<float>* dlm_pT_d = NULL;
  DLM_Histo<float> dlm_pT_eta_p;
  DLM_Histo<float> dlm_pT_eta_d;
  TH1F* h_pT_p_all = NULL;
  TH1F* h_pT_d_all = NULL;
  TH1F* h_pT_p = NULL;
  TH1F* h_pT_ap = NULL;
  TH1F* h_pT_d = NULL;
  TH1F* h_pT_ad = NULL;

  TString FilePath;
  if(type=="pd"){
    FilePath = TString::Format("%s/CatsFiles/Source/CECA/pd/pT_spectra/",GetCernBoxDimi());
    TString FileNameP1,FileNameAP1,FileNameP2,FileNameAP2;
    FileNameP1 = "protonSpectra.root";
    FileNameAP1 = "AntiProtonSpectra.root";
    FileNameP2 = "deuteronSpectra.root";
    FileNameAP2 = "AntiDeuteronSpectra.root";

    TFile file_p(FilePath+FileNameP1,"read");
    h_pT_p = (TH1F*)file_p.Get("pTDist_after");
    if(!h_pT_p) printf("ISSUE with h_pT_p\n");
    fOutput.cd();
    h_pT_p_all = (TH1F*)h_pT_p->Clone("h_pT_p_all");

    TFile file_ap(FilePath+FileNameAP1,"read");
    h_pT_ap = (TH1F*)file_ap.Get("pTDist_after");
    if(!h_pT_ap) printf("ISSUE with h_pT_ap\n");
    h_pT_p_all->Add(h_pT_ap);

    TFile file_d(FilePath+FileNameP2,"read");
    h_pT_d = (TH1F*)file_d.Get("pTDist_after");
    if(!h_pT_d) printf("ISSUE with h_pT_d\n");
    fOutput.cd();
    h_pT_d_all = (TH1F*)h_pT_d->Clone("h_pT_d_all");

    TFile file_ad(FilePath+FileNameAP2,"read");
    h_pT_ad = (TH1F*)file_ad.Get("pTDist_after");
    if(!h_pT_ad) printf("ISSUE with h_pT_ad\n");
    h_pT_d_all->Add(h_pT_ad);
  }
  if(type=="pp"||type=="pP"){
    FilePath = TString::Format("%s/CatsFiles/Source/CECA/pd/pT_spectra/",GetCernBoxDimi());
    TString FileNameP1,FileNameAP1,FileNameP2,FileNameAP2;
    FileNameP1 = "protonSpectra.root";
    FileNameAP1 = "AntiProtonSpectra.root";
    FileNameP2 = "protonSpectra.root";
    FileNameAP2 = "AntiProtonSpectra.root";

    TFile file_p(FilePath+FileNameP1,"read");
    h_pT_p = (TH1F*)file_p.Get("pTDist_after");
    if(!h_pT_p) printf("ISSUE with h_pT_p\n");
    fOutput.cd();
    h_pT_p_all = (TH1F*)h_pT_p->Clone("h_pT_p_all");

    TFile file_ap(FilePath+FileNameAP1,"read");
    h_pT_ap = (TH1F*)file_ap.Get("pTDist_after");
    if(!h_pT_ap) printf("ISSUE with h_pT_ap\n");
    h_pT_p_all->Add(h_pT_ap);

    TFile file_d(FilePath+FileNameP2,"read");
    h_pT_d = (TH1F*)file_d.Get("pTDist_after");
    if(!h_pT_d) printf("ISSUE with h_pT_d\n");
    fOutput.cd();
    h_pT_d_all = (TH1F*)h_pT_d->Clone("h_pT_d_all");

    TFile file_ad(FilePath+FileNameAP2,"read");
    h_pT_ad = (TH1F*)file_ad.Get("pTDist_after");
    if(!h_pT_ad) printf("ISSUE with h_pT_ad\n");
    h_pT_d_all->Add(h_pT_ad);
  }
  //N.B. here we set the Kaons to the proton histo
  if(type=="Kd"){
    FilePath = TString::Format("%s/CatsFiles/Source/CECA/Kd/",GetCernBoxDimi());
    TFile file(FilePath+"outKD_mydeuteronsOpenSyst.root","read");
    TDirectoryFile *dir=(TDirectoryFile*)(file.FindObjectAny(Form("PtEtaPhi")));

    dir->GetObject(Form("h_kp_pT"),h_pT_p);
    dir->GetObject(Form("h_km_pT"),h_pT_ap);
    fOutput.cd();
    h_pT_p_all = (TH1F*)h_pT_p->Clone("h_pT_p_all");
    h_pT_p_all->Add(h_pT_ap);

    dir->GetObject(Form("h_d_pT"),h_pT_d);
    dir->GetObject(Form("h_ad_pT"),h_pT_ad);
    fOutput.cd();
    h_pT_d_all = (TH1F*)h_pT_d->Clone("h_pT_d_all");
    h_pT_d_all->Add(h_pT_ad);
  }
  //N.B. here we set the Kaons to the deuteron histo
  if(type=="pK"){
    FilePath = TString::Format("%s/CatsFiles/Source/CECA/Kd/",GetCernBoxDimi());
    TFile fileK(FilePath+"outKD_mydeuteronsOpenSyst.root","read");
    TDirectoryFile *dir=(TDirectoryFile*)(fileK.FindObjectAny(Form("PtEtaPhi")));
    dir->GetObject(Form("h_kp_pT"),h_pT_d);
    dir->GetObject(Form("h_km_pT"),h_pT_ad);
    fOutput.cd();
    h_pT_d_all = (TH1F*)h_pT_d->Clone("h_pT_d_all");
    h_pT_d_all->Add(h_pT_ad);


    FilePath = TString::Format("%s/CatsFiles/Source/CECA/pd/pT_spectra/",GetCernBoxDimi());
    TString FileNameP1,FileNameAP1,FileNameP2,FileNameAP2;
    FileNameP1 = "protonSpectra.root";
    FileNameAP1 = "AntiProtonSpectra.root";
    FileNameP2 = "protonSpectra.root";
    FileNameAP2 = "AntiProtonSpectra.root";

    TFile file_p(FilePath+FileNameP1,"read");
    h_pT_p = (TH1F*)file_p.Get("pTDist_after");
    if(!h_pT_p) printf("ISSUE with h_pT_p\n");
    fOutput.cd();
    h_pT_p_all = (TH1F*)h_pT_p->Clone("h_pT_p_all");

    TFile file_ap(FilePath+FileNameAP1,"read");
    h_pT_ap = (TH1F*)file_ap.Get("pTDist_after");
    if(!h_pT_ap) printf("ISSUE with h_pT_ap\n");
    h_pT_p_all->Add(h_pT_ap);

  }

  if(h_pT_p_all){
    dlm_pT_p = Convert_TH1F_DlmHisto(h_pT_p_all);
    dlm_pT_p->RescaleAxis(0,1000,false);
  }
  if(h_pT_d_all){
    dlm_pT_d = Convert_TH1F_DlmHisto(h_pT_d_all);
    dlm_pT_d->RescaleAxis(0,1000,false);
  }

  double* BinRange = NULL;
  double axis[2];

  if(dlm_pT_p){
    dlm_pT_eta_p.SetUp(2);
    BinRange = dlm_pT_p->GetBinRange(0);
    dlm_pT_eta_p.SetUp(0,dlm_pT_p->GetNbins(),BinRange);
    delete [] BinRange;
    dlm_pT_eta_p.SetUp(1,1,-EtaCut,EtaCut);
    dlm_pT_eta_p.Initialize();
    for(unsigned uBin=0; uBin<dlm_pT_p->GetNbins(); uBin++){
      dlm_pT_eta_p.SetBinContent(uBin,0,dlm_pT_p->GetBinContent(uBin));
      //printf("b%u %.3e\n",uBin,dlm_pT_p->GetBinContent(uBin));
    }
  }

  if(dlm_pT_d){
    dlm_pT_eta_d.SetUp(2);
    BinRange = dlm_pT_d->GetBinRange(0);
    dlm_pT_eta_d.SetUp(0,dlm_pT_d->GetNbins(),BinRange);
    delete [] BinRange;
    dlm_pT_eta_d.SetUp(1,1,-EtaCut,EtaCut);
    dlm_pT_eta_d.Initialize();
    for(unsigned uBin=0; uBin<dlm_pT_d->GetNbins(); uBin++){
      dlm_pT_eta_d.SetBinContent(uBin,0,dlm_pT_d->GetBinContent(uBin));
    }
  }

  TH2F* hSampleQA_p = new TH2F("hSampleQA_p","hSampleQA_p",64,-1,8,64,-1,1);
  if(dlm_pT_p){
    for(unsigned uIter=0; uIter<100*1000; uIter++){
      double axisValues[2];
      dlm_pT_eta_p.Sample(axisValues);
      hSampleQA_p->Fill(axisValues[0],axisValues[1]);
    }
  }


  fOutput.cd();
  hSampleQA_p->Write();
  if(h_pT_p_all)h_pT_p_all->Write();
  if(h_pT_d_all)h_pT_d_all->Write();

  for(TreParticle* prt : ParticleList){
    if(prt->GetName()=="Proton"||prt->GetName()=="PrimProton"){
      //if(prt->GetName()=="PrimProton")prt->SetMass(2.*Mass_p);
      //else prt->SetMass(Mass_p);
      prt->SetMass(Mass_p);
      if(type=="pd"||type=="pp"||type=="pK") prt->SetAbundance(35.78+64.22*(!PROTON_RESO));
      else if(type=="pP") prt->SetAbundance(100.);
      else prt->SetAbundance(0);
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);
      if(dlm_pT_p) prt->SetPtEtaPhi(dlm_pT_eta_p);
      else
        prt->SetPtPz(0.85*prt->GetMass(),0.85*prt->GetMass());
    }
    else if(prt->GetName()=="Kaon"){
      prt->SetMass(Mass_Kch);
      if(type=="Kd"||type=="pK") prt->SetAbundance(53.4+47.6*(!PROTON_RESO));
    //if(type=="Kd"||type=="pK") prt->SetAbundance(100);
      else prt->SetAbundance(0);
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);

      if(type=="Kd"){
        if(dlm_pT_p) prt->SetPtEtaPhi(dlm_pT_eta_p);
        else prt->SetPtPz(0.85*prt->GetMass(),0.85*prt->GetMass());
      }
      else if(type=="pK"){
        if(dlm_pT_d) prt->SetPtEtaPhi(dlm_pT_eta_d);
        else prt->SetPtPz(0.85*prt->GetMass(),0.85*prt->GetMass());
      }
      else{
        prt->SetPtPz(0.85*prt->GetMass(),0.85*prt->GetMass());
      }
    }
    else if(prt->GetName()=="Deuteron"){
      prt->SetMass(Mass_d);
      if(type=="pd"||type=="Kd") prt->SetAbundance(100.0);
      else prt->SetAbundance(0);
      prt->SetRadius(HadronSize*1.5);
      prt->SetRadiusSlope(HadronSlope);
      prt->SetDelayTau(d_delay);
      if(dlm_pT_d)prt->SetPtEtaPhi(dlm_pT_eta_d);
      else prt->SetPtPz(0.85*prt->GetMass(),0.85*prt->GetMass());
    }
    else if(prt->GetName()=="Pion"){
      prt->SetMass(Mass_pic);
      prt->SetAbundance(0);
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);
    }
    else if(prt->GetName()=="PionFSI"){
      prt->SetMass(Mass_pic);
      //59.8 is the amount of reso below 5 fm vs 28 prim, here we normed this to 100
      if(type=="pipi") prt->SetAbundance(31.8+68.2*(!PROTON_RESO));
      else if(type=="pipi_core") prt->SetAbundance(100);
      else if(type=="Dpi"){
        prt->SetAbundance(31.8);
        prt->SetDelayTau(d_delay);
      }
      else prt->SetAbundance(0);
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);
      prt->SetPtPz(0.85*prt->GetMass(),0.85*prt->GetMass());
      prt->SetAcceptance_pT(140,100000);
      prt->SetAcceptance_Eta(-EtaCut,EtaCut);
    }
    else if(prt->GetName()=="Dch"){
      prt->SetMass(Mass_Dch);
      if(type=="Dpi") prt->SetAbundance(70.5);
      else prt->SetAbundance(0);
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);
      prt->SetPtPz(0.85*prt->GetMass(),0.85*prt->GetMass());
    }
    else if(prt->GetName()=="Dch_star"){
      prt->SetMass(Mass_Dch_star);
      //if(type=="Dpi") prt->SetAbundance(29.5);
      prt->SetAbundance(0);
      prt->SetWidth(0.0834);
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);
      prt->SetPtPz(0.85*prt->GetMass(),0.85*prt->GetMass());

      prt->NewDecay();
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Dch"));
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
      prt->GetDecay(0)->SetBranching(100);
    }
    else if(prt->GetName()=="ProtonReso"){
      prt->SetMass(1362);
      if(type=="pd"||type=="pp"||type=="pP"||type=="pK") prt->SetAbundance(64.22*PROTON_RESO);
      else prt->SetAbundance(0);
      prt->SetWidth(hbarc/1.65);
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);
      if(dlm_pT_p)prt->SetPtEtaPhi(dlm_pT_eta_p);
      else
        prt->SetPtPz(0.85*prt->GetMass(),0.85*prt->GetMass());


      prt->NewDecay();
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Proton"));
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
      prt->GetDecay(0)->SetBranching(100);
    }
    else if(prt->GetName()=="KaonReso"){
      prt->SetMass(1054);
      if(type=="Kd"||type=="pK") prt->SetAbundance(47.6*PROTON_RESO);
    //if(type=="Kd"||type=="pK") prt->SetAbundance(0);
      else prt->SetAbundance(0);
      prt->SetWidth(hbarc/3.66);
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);

      if(type=="Kd"){
        if(dlm_pT_p) prt->SetPtEtaPhi(dlm_pT_eta_p);
        else prt->SetPtPz(0.85*prt->GetMass(),0.85*prt->GetMass());
      }
      else if(type=="pK"){
        if(dlm_pT_d) prt->SetPtEtaPhi(dlm_pT_eta_d);
        else prt->SetPtPz(0.85*prt->GetMass(),0.85*prt->GetMass());
      }
      else{
        prt->SetPtPz(0.85*prt->GetMass(),0.85*prt->GetMass());
      }

      prt->NewDecay();
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Kaon"));
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
      prt->GetDecay(0)->SetBranching(100);
    }
    else if(prt->GetName()=="PionReso"){
      prt->SetMass(1180);
      if(type=="pipi") prt->SetAbundance(68.2*(PROTON_RESO));
      else if(type=="pipi_core") prt->SetAbundance(0);
      else if(type=="Dpi") prt->SetAbundance(68.2);
      else prt->SetAbundance(0);
      prt->SetWidth(hbarc/1.50);
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);
      prt->SetPtPz(0.85*prt->GetMass(),0.85*prt->GetMass());
//  const double Tau_pi = 1.50;
  //const double Mass_ProtonReso = 1362;
  //const double Mass_piReso = 1180;
      prt->NewDecay();
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("PionFSI"));
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
      prt->GetDecay(0)->SetBranching(100);

      prt->SetAcceptance_pT(140,100000);
      prt->SetAcceptance_Eta(-EtaCut,EtaCut);
    }


  }

  std::vector<std::string> ListOfParticles;
  if(type=="pd"||type=="pK") ListOfParticles.push_back("Proton");
  if(type=="Kd"||type=="pK") ListOfParticles.push_back("Kaon");
  if(type=="pd"||type=="Kd") ListOfParticles.push_back("Deuteron");

  if(type=="Dpi"){
    ListOfParticles.push_back("Dch");
    ListOfParticles.push_back("PionFSI");
  }
  if(type=="pp"){
    ListOfParticles.push_back("Proton");
    ListOfParticles.push_back("Proton");
  }
  if(type=="pP"){
    ListOfParticles.push_back("Proton");
    ListOfParticles.push_back("PrimProton");
  }
  if(type=="pipi"||type=="pipi_core"){
    ListOfParticles.push_back("PionFSI");
    ListOfParticles.push_back("PionFSI");
  }

//printf("%s %s\n", ListOfParticles.at(0).c_str(),ListOfParticles.at(1).c_str());
//return;
  CECA Ivana(Database,ListOfParticles);

  //Ivana.SetDisplacement(rSP_core);
  //Ivana.SetHadronization(0.0);
  //Ivana.SetHadrFluctuation(0.0);

  Ivana.SetDisplacementZ(rSP_dispZ);
  Ivana.SetDisplacementT(rSP_core);
  Ivana.SetHadronizationZ(0);
  Ivana.SetHadronizationT(rSP_hadr);
  Ivana.SetHadrFluctuation(rSP_hflc);
  Ivana.SetTau(rSP_tau,tau_prp);
  Ivana.SetTauFluct(rSP_tflc);
  Ivana.SetThermalKick(rSP_ThK);
  Ivana.SetFixedHadr(rSP_FixedHadr);
  Ivana.SetFragmentBeta(rSP_FragBeta);

  //Ivana.SetDisplacement(3.55*0.25);
  //Ivana.SetHadronizationT(3.55);
  //Ivana.SetHadrFluctuation(0.0);
  //Ivana.SetTau(0);

  //Ivana.SetDisplacement(1);
  //Ivana.SetHadronizationT(0);
  //Ivana.SetTau(0);

  Ivana.SetTargetStatistics(target_yield);
  Ivana.SetEventMult(Multiplicity);
  Ivana.SetSourceDim(2);
  Ivana.SetDebugMode(true);
  Ivana.SetThreadTimeout(TIMEOUT);
  Ivana.SetGlobalTimeout(GLOB_TIMEOUT);
  Ivana.EqualizeFsiTime(EQUALIZE_TAU);
  Ivana.SetFemtoRegion(femto_region);
  Ivana.GHETTO_EVENT = true;
  if(type=="pp"){
    //BinCenter_pp[0] = 1.1077;
    //BinCenter_pp[1] = 1.1683;
    //BinCenter_pp[2] = 1.2284;
    //BinCenter_pp[3] = 1.3156;
    //BinCenter_pp[4] = 1.4628;
    //BinCenter_pp[5] = 1.6872;
    //BinCenter_pp[6] = 2.2116;

    Ivana.Ghetto_NumMtBins = 9;
    Ivana.Ghetto_MtBins = new double [Ivana.Ghetto_NumMtBins+1];
    Ivana.Ghetto_MtBins[0] = 938;
    Ivana.Ghetto_MtBins[1] = 1055;
    Ivana.Ghetto_MtBins[2] = 1135;
    Ivana.Ghetto_MtBins[3] = 1190;
    Ivana.Ghetto_MtBins[4] = 1270;
    Ivana.Ghetto_MtBins[5] = 1390;
    Ivana.Ghetto_MtBins[6] = 1570;
    Ivana.Ghetto_MtBins[7] = 1940;
    Ivana.Ghetto_MtBins[8] = 2500;
    Ivana.Ghetto_MtBins[9] = 4000;

    Ivana.Ghetto_NumMomBins = 150;
    Ivana.Ghetto_MomMin = 0;
    Ivana.Ghetto_MomMax = 600;
  }

  Ivana.GoBabyGo();

  Ivana.Ghetto_kstar_rstar_mT->QuickWrite(BaseFileName+".Ghetto_kstar_rstar_mT",true);

//return;
  double TotPairs = Ivana.GhettoPrimReso[0]+Ivana.GhettoPrimReso[1]+Ivana.GhettoPrimReso[2]+Ivana.GhettoPrimReso[3];
  double TotPP = double(Ivana.GhettoPrimReso[0])/TotPairs;
  double TotPR = double(Ivana.GhettoPrimReso[1])/TotPairs;
  double TotRP = double(Ivana.GhettoPrimReso[2])/TotPairs;
  double TotRR = double(Ivana.GhettoPrimReso[3])/TotPairs;

  double FemtoPairs = Ivana.GhettoFemtoPrimReso[0]+Ivana.GhettoFemtoPrimReso[1]+Ivana.GhettoFemtoPrimReso[2]+Ivana.GhettoFemtoPrimReso[3];
  double FemtoPP = double(Ivana.GhettoFemtoPrimReso[0])/FemtoPairs;
  double FemtoPR = double(Ivana.GhettoFemtoPrimReso[1])/FemtoPairs;
  double FemtoRP = double(Ivana.GhettoFemtoPrimReso[2])/FemtoPairs;
  double FemtoRR = double(Ivana.GhettoFemtoPrimReso[3])/FemtoPairs;

  printf("     Total  Femto\n");
  printf("PP%6.2f%% %6.2f\n",TotPP*100.,FemtoPP*100.);
  printf("PR%6.2f%% %6.2f\n",TotPR*100.,FemtoPR*100.);
  printf("RP%6.2f%% %6.2f\n",TotRP*100.,FemtoRP*100.);
  printf("RR%6.2f%% %6.2f\n",TotRR*100.,FemtoRR*100.);

  Ivana.GhettoFemto_rstar->ComputeError();
  Ivana.GhettoFemto_rstar->ScaleToIntegral();
  Ivana.GhettoFemto_rstar->ScaleToBinSize();

  CATS Kitty_pp;
  Kitty_pp.SetMomBins(80,0,320);
  Kitty_pp.SetNotifications(CATS::nWarning);
  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());
  AnalysisObject.SetUpCats_pp(Kitty_pp,"AV18","NULL",0,0);
  DLM_HistoSource ppHistoSource(*Ivana.GhettoFemto_rstar);
  Kitty_pp.SetAnaSource(CatsSourceForwarder, &ppHistoSource, 0);
  Kitty_pp.SetUseAnalyticSource(true);
  Kitty_pp.SetAutoNormSource(false);
  Kitty_pp.SetNormalizedSource(true);
  Kitty_pp.KillTheCat();
  TGraph Ck_pp;
  Ck_pp.SetName("Ck_pp");
  Ck_pp.SetLineColor(kBlue);
  Ck_pp.SetLineWidth(5);
  for(unsigned uBin=0; uBin<Kitty_pp.GetNumMomBins(); uBin++){
    Ck_pp.SetPoint(uBin,Kitty_pp.GetMomentum(uBin),Kitty_pp.GetCorrFun(uBin));
  }
  fOutput.cd();
  Ck_pp.Write();

  TH1F* h_GhettoFemto_rstar = Convert_DlmHisto_TH1F(Ivana.GhettoFemto_rstar,"GhettoFemto_rstar");
  fOutput.cd();
  h_GhettoFemto_rstar->SetLineWidth(3);
  h_GhettoFemto_rstar->SetLineColor(kAzure);
  h_GhettoFemto_rstar->Write();
  TF1* fit_rstar = new TF1("fit_rstar","[0]*4.*TMath::Pi()*x*x*pow(4.*TMath::Pi()*[1]*[1],-1.5)*exp(-(x*x)/(4.*[1]*[1]))+1.-[0]",0,16);
  fit_rstar->SetLineColor(kAzure);
  fit_rstar->SetLineWidth(3);
  //plot the fit to the total dist as an example
  //below we find the same for each mT, there we will just get out the reff_Ceca
  double reff_Ceca = Get_reff(h_GhettoFemto_rstar);
  printf("reff_Ceca = %f\n",reff_Ceca);
  fOutput.cd();
  fit_rstar->FixParameter(0,1);
  fit_rstar->FixParameter(1,reff_Ceca);
  fit_rstar->Write();




  Ivana.GhettoFemto_rcore->ComputeError();
  Ivana.GhettoFemto_rcore->ScaleToIntegral();
  Ivana.GhettoFemto_rcore->ScaleToBinSize();
  TH1F* h_GhettoFemto_rcore = Convert_DlmHisto_TH1F(Ivana.GhettoFemto_rcore,"GhettoFemto_rcore");
  fOutput.cd();
  h_GhettoFemto_rcore->SetLineWidth(3);
  h_GhettoFemto_rcore->SetLineColor(kBlack);
  h_GhettoFemto_rcore->Write();
  TF1* fit_rcore = new TF1("fit_rcore","[0]*4.*TMath::Pi()*x*x*pow(4.*TMath::Pi()*[1]*[1],-1.5)*exp(-(x*x)/(4.*[1]*[1]))+1.-[0]",0,16);
  fit_rcore->SetLineColor(kBlack);
  fit_rcore->SetLineWidth(3);
  double rcore_Ceca = Get_reff(h_GhettoFemto_rcore);
  printf("rcore_Ceca = %f\n",rcore_Ceca);
  fOutput.cd();
  fit_rcore->FixParameter(0,1);
  fit_rcore->FixParameter(1,rcore_Ceca);
  fit_rcore->Write();


  Ivana.Ghetto_kstar->ComputeError();
  Ivana.Ghetto_kstar->ScaleToIntegral();
  Ivana.Ghetto_kstar->ScaleToBinSize();
  TH1F* h_Ghetto_kstar = Convert_DlmHisto_TH1F(Ivana.Ghetto_kstar,"Ghetto_kstar");
  Ivana.Ghetto_kstar_rstar->ComputeError();
  TH2F* h_Ghetto_kstar_rstar = Convert_DlmHisto_TH2F(Ivana.Ghetto_kstar_rstar,"Ghetto_kstar_rstar");

  Ivana.Ghetto_kstar_rstar_PP->ComputeError();
  TH2F* h_Ghetto_kstar_rstar_PP = Convert_DlmHisto_TH2F(Ivana.Ghetto_kstar_rstar_PP,"Ghetto_kstar_rstar_PP");
  Ivana.Ghetto_kstar_rstar_PR->ComputeError();
  TH2F* h_Ghetto_kstar_rstar_PR = Convert_DlmHisto_TH2F(Ivana.Ghetto_kstar_rstar_PR,"Ghetto_kstar_rstar_PR");
  Ivana.Ghetto_kstar_rstar_RP->ComputeError();
  TH2F* h_Ghetto_kstar_rstar_RP = Convert_DlmHisto_TH2F(Ivana.Ghetto_kstar_rstar_RP,"Ghetto_kstar_rstar_RP");
  Ivana.Ghetto_kstar_rstar_RR->ComputeError();
  TH2F* h_Ghetto_kstar_rstar_RR = Convert_DlmHisto_TH2F(Ivana.Ghetto_kstar_rstar_RR,"Ghetto_kstar_rstar_RR");

  TH2F* h_Ghetto_mT_rstar = Convert_DlmHisto_TH2F(Ivana.Ghetto_mT_rstar,"Ghetto_mT_rstar");
  Ivana.GhettoFemto_mT_rstar->ComputeError();
  TH2F* h_GhettoFemto_mT_rstar = Convert_DlmHisto_TH2F(Ivana.GhettoFemto_mT_rstar,"GhettoFemto_mT_rstar");

//fOutput.cd();
//h_GhettoFemto_mT_rstar->Write();
  Ivana.GhettoFemto_mT_rcore->ComputeError();
  TH2F* h_GhettoFemto_mT_rcore = Convert_DlmHisto_TH2F(Ivana.GhettoFemto_mT_rcore,"GhettoFemto_mT_rcore");
  TH2F* h_GhettoFemto_mT_kstar = Convert_DlmHisto_TH2F(Ivana.GhettoFemto_mT_kstar,"GhettoFemto_mT_kstar");
  TH2F* h_Ghetto_mT_costh = Convert_DlmHisto_TH2F(Ivana.Ghetto_mT_costh,"Ghetto_mT_costh");
  TH2F* h_GhettoSP_pT_th = Convert_DlmHisto_TH2F(Ivana.GhettoSP_pT_th,"GhettoSP_pT_th");
  TH1F* h_GhettoSP_pT_1 = Convert_DlmHisto_TH1F(Ivana.GhettoSP_pT_1,"GhettoSP_pT_1");
  TH1F* h_GhettoSP_pT_2 = Convert_DlmHisto_TH1F(Ivana.GhettoSP_pT_2,"GhettoSP_pT_2");
  //Ivana.GhettoFemto_rstar->ComputeError();
  //TH1F* h_GhettoFemto_rstar = Convert_DlmHisto_TH1F(Ivana.GhettoFemto_rstar,"GhettoFemto_rstar");
  Ivana.Ghetto_PP_AngleRcP1->ComputeError();
  Ivana.Ghetto_PP_AngleRcP2->ComputeError();
  Ivana.Ghetto_PP_AngleP1P2->ComputeError();
  Ivana.Ghetto_RP_AngleRcP1->ComputeError();
  Ivana.Ghetto_PR_AngleRcP2->ComputeError();
  Ivana.Ghetto_RR_AngleRcP1->ComputeError();
  Ivana.Ghetto_RR_AngleRcP2->ComputeError();
  Ivana.Ghetto_RR_AngleP1P2->ComputeError();
  TH1F* h_Ghetto_PP_AngleRcP1 = Convert_DlmHisto_TH1F(Ivana.Ghetto_PP_AngleRcP1,"Ghetto_PP_AngleRcP1");
  TH1F* h_Ghetto_PP_AngleRcP2 = Convert_DlmHisto_TH1F(Ivana.Ghetto_PP_AngleRcP2,"Ghetto_PP_AngleRcP2");
  TH1F* h_Ghetto_PP_AngleP1P2 = Convert_DlmHisto_TH1F(Ivana.Ghetto_PP_AngleP1P2,"Ghetto_PP_AngleP1P2");
  TH1F* h_Ghetto_RP_AngleRcP1 = Convert_DlmHisto_TH1F(Ivana.Ghetto_RP_AngleRcP1,"Ghetto_RP_AngleRcP1");
  TH1F* h_Ghetto_PR_AngleRcP2 = Convert_DlmHisto_TH1F(Ivana.Ghetto_PR_AngleRcP2,"Ghetto_PR_AngleRcP2");
  TH1F* h_Ghetto_RR_AngleRcP1 = Convert_DlmHisto_TH1F(Ivana.Ghetto_RR_AngleRcP1,"Ghetto_RR_AngleRcP1");
  TH1F* h_Ghetto_RR_AngleRcP2 = Convert_DlmHisto_TH1F(Ivana.Ghetto_RR_AngleRcP2,"Ghetto_RR_AngleRcP2");
  TH1F* h_Ghetto_RR_AngleP1P2 = Convert_DlmHisto_TH1F(Ivana.Ghetto_RR_AngleP1P2,"Ghetto_RR_AngleP1P2");

  Ivana.GhettoSPr_X->ComputeError();
  Ivana.GhettoSPr_X->ScaleToIntegral();
  Ivana.GhettoSPr_X->ScaleToBinSize();
  TH1F* h_GhettoSPr_X = Convert_DlmHisto_TH1F(Ivana.GhettoSPr_X,"GhettoSPr_X");

  Ivana.GhettoSPr_Y->ComputeError();
  Ivana.GhettoSPr_Y->ScaleToIntegral();
  Ivana.GhettoSPr_Y->ScaleToBinSize();
  TH1F* h_GhettoSPr_Y = Convert_DlmHisto_TH1F(Ivana.GhettoSPr_Y,"GhettoSPr_Y");

  Ivana.GhettoSPr_Z->ComputeError();
  Ivana.GhettoSPr_Z->ScaleToIntegral();
  Ivana.GhettoSPr_Z->ScaleToBinSize();
  TH1F* h_GhettoSPr_Z = Convert_DlmHisto_TH1F(Ivana.GhettoSPr_Z,"GhettoSPr_Z");

  Ivana.GhettoSPr_Rho->ComputeError();
  Ivana.GhettoSPr_Rho->ScaleToIntegral();
  Ivana.GhettoSPr_Rho->ScaleToBinSize();
  TH1F* h_GhettoSPr_Rho = Convert_DlmHisto_TH1F(Ivana.GhettoSPr_Rho,"GhettoSPr_Rho");

  Ivana.GhettoSPr_R->ComputeError();
  Ivana.GhettoSPr_R->ScaleToIntegral();
  Ivana.GhettoSPr_R->ScaleToBinSize();
  TH1F* h_GhettoSPr_R = Convert_DlmHisto_TH1F(Ivana.GhettoSPr_R,"GhettoSPr_R");

  Ivana.GhettoSP_X->ComputeError();
  Ivana.GhettoSP_X->ScaleToIntegral();
  Ivana.GhettoSP_X->ScaleToBinSize();
  TH1F* h_GhettoSP_X = Convert_DlmHisto_TH1F(Ivana.GhettoSP_X,"GhettoSP_X");

  Ivana.GhettoSP_Y->ComputeError();
  Ivana.GhettoSP_Y->ScaleToIntegral();
  Ivana.GhettoSP_Y->ScaleToBinSize();
  TH1F* h_GhettoSP_Y = Convert_DlmHisto_TH1F(Ivana.GhettoSP_Y,"GhettoSP_Y");

  Ivana.GhettoSP_Z->ComputeError();
  Ivana.GhettoSP_Z->ScaleToIntegral();
  Ivana.GhettoSP_Z->ScaleToBinSize();
  TH1F* h_GhettoSP_Z = Convert_DlmHisto_TH1F(Ivana.GhettoSP_Z,"GhettoSP_Z");

  Ivana.GhettoSP_Rho->ComputeError();
  Ivana.GhettoSP_Rho->ScaleToIntegral();
  Ivana.GhettoSP_Rho->ScaleToBinSize();
  TH1F* h_GhettoSP_Rho = Convert_DlmHisto_TH1F(Ivana.GhettoSP_Rho,"GhettoSP_Rho");

  Ivana.GhettoSP_R->ComputeError();
  Ivana.GhettoSP_R->ScaleToIntegral();
  Ivana.GhettoSP_R->ScaleToBinSize();
  TH1F* h_GhettoSP_R = Convert_DlmHisto_TH1F(Ivana.GhettoSP_R,"GhettoSP_R");

  TGraphErrors g_GhettoFemto_mT_rstar;
  g_GhettoFemto_mT_rstar.SetName("g_GhettoFemto_mT_rstar");
  g_GhettoFemto_mT_rstar.SetMarkerStyle(20);
  g_GhettoFemto_mT_rstar.SetMarkerSize(1);
  g_GhettoFemto_mT_rstar.SetLineWidth(3);

  TGraphErrors g_GhettoFemto_mT_rstar_G;
  g_GhettoFemto_mT_rstar_G.SetName("g_GhettoFemto_mT_rstar_G");
  g_GhettoFemto_mT_rstar_G.SetMarkerStyle(20);
  g_GhettoFemto_mT_rstar_G.SetMarkerSize(1);
  g_GhettoFemto_mT_rstar_G.SetLineWidth(6);
  g_GhettoFemto_mT_rstar_G.SetLineColor(kAzure);

  TGraphErrors g_GhettoFemto_mT_rcore;
  g_GhettoFemto_mT_rcore.SetName("g_GhettoFemto_mT_rcore");
  g_GhettoFemto_mT_rcore.SetMarkerStyle(20);
  g_GhettoFemto_mT_rcore.SetMarkerSize(1);
  g_GhettoFemto_mT_rcore.SetLineWidth(3);

  TGraphErrors g_GhettoFemto_mT_rcore_G;
  g_GhettoFemto_mT_rcore_G.SetName("g_GhettoFemto_mT_rcore_G");
  g_GhettoFemto_mT_rcore_G.SetMarkerStyle(20);
  g_GhettoFemto_mT_rcore_G.SetMarkerSize(1);
  g_GhettoFemto_mT_rcore_G.SetLineWidth(6);
  g_GhettoFemto_mT_rcore_G.SetLineColor(kBlack);

  unsigned uPointRS = 0;
  unsigned uPointRC = 0;
  for(unsigned uBin=0; uBin<h_GhettoFemto_mT_rstar->GetXaxis()->GetNbins(); uBin++){
    TH1F* hProj = (TH1F*)h_GhettoFemto_mT_rstar->ProjectionY(TString::Format("hProj"),uBin+1,uBin+1);
    double Mean = hProj->GetMean();
    double Err = hProj->GetStdDev();
    double mT = h_GhettoFemto_mT_rstar->GetXaxis()->GetBinCenter(uBin+1);
    if(Mean&&Err&&hProj->GetEntries()>256){
      hProj->Scale(1./hProj->Integral(),"width");
      g_GhettoFemto_mT_rstar.SetPoint(uPointRS,mT*0.001,Mean);
      //g_GhettoFemto_mT_rstar.SetPointError(uPointRS,0,Err);
      g_GhettoFemto_mT_rstar.SetPointError(uPointRS,0,0);

      reff_Ceca = Get_reff(hProj);
      //printf("mT%u reff_Ceca=%.3f\n",uBin,reff_Ceca);
      g_GhettoFemto_mT_rstar_G.SetPoint(uPointRS,mT*0.001,reff_Ceca);
      g_GhettoFemto_mT_rstar_G.SetPointError(uPointRS,0,0);

      uPointRS++;
    }
    delete hProj;

    //do the same for
    //GhettoFemto_mT_rcore

    hProj = (TH1F*)h_GhettoFemto_mT_rcore->ProjectionY(TString::Format("hProj"),uBin+1,uBin+1);
    Mean = hProj->GetMean();
    Err = hProj->GetStdDev();
    mT = h_GhettoFemto_mT_rcore->GetXaxis()->GetBinCenter(uBin+1);
    if(Mean&&Err&&hProj->GetEntries()>256){
      hProj->Scale(1./hProj->Integral(),"width");
      g_GhettoFemto_mT_rcore.SetPoint(uPointRC,mT*0.001,Mean);
      g_GhettoFemto_mT_rcore.SetPointError(uPointRC,0,0);

      rcore_Ceca = Get_reff(hProj);
      //printf("    rcore_Ceca=%.3f\n",rcore_Ceca);
      g_GhettoFemto_mT_rcore_G.SetPoint(uPointRC,mT*0.001,rcore_Ceca);
      g_GhettoFemto_mT_rcore_G.SetPointError(uPointRC,0,0);

      uPointRC++;
    }
    delete hProj;
  }

  const double MeanToGauss = 1./2.3;
  double mT_first = g_GhettoFemto_mT_rstar.GetPointY(0);
  printf("r(%.2f GeV) = %.2f\n",g_GhettoFemto_mT_rstar.GetPointX(0)*0.001,mT_first*MeanToGauss);
  double mT_2GeV = g_GhettoFemto_mT_rstar.Eval(2000);
  printf("r(%.2f GeV) = %.2f\n",2.,mT_2GeV*MeanToGauss);
  //printf("----FemTUM----\n");
  //printf("CORE: <rcore> = %.3f; rcore = %.3f\n",);

  h_GhettoFemto_rstar->Scale(1./h_GhettoFemto_rstar->Integral(),"width");

  double lowerlimit;
  double upperlimit;
  GetCentralInterval(*h_GhettoFemto_rstar, 0.9, lowerlimit, upperlimit, true);
  unsigned lowerbin = h_GhettoFemto_rstar->FindBin(lowerlimit);
  unsigned upperbin = h_GhettoFemto_rstar->FindBin(upperlimit);
  printf("The fit will be performed in the range [%.2f, %.2f] fm\n",lowerlimit,upperlimit);
  TF1* fSource = new TF1("fSource",ScaledGauss,0,32,2);
  fSource->SetParameter(0,h_GhettoFemto_rstar->GetMean()/2.3);
  fSource->SetParLimits(0,0.5,5.0);
  //fSource->SetParameter(1,0.5);
  //fSource->SetParLimits(1,0.1,1.0);
  fSource->FixParameter(1,1.0);

  h_GhettoFemto_rstar->Fit(fSource,"Q, S, N, R, M","",lowerlimit,upperlimit);
  printf("The effective Gaussian size is %.3f +/- %.3f fm, λ = %.3f +/- %.3f\n",
  fSource->GetParameter(0),fSource->GetParError(0),fSource->GetParameter(1),fSource->GetParError(1));

  GetCentralInterval(*h_GhettoFemto_rstar, 0.9, lowerlimit, upperlimit, true);
  //TF1* fitDG_rstar = new TF1("fitDG_rstar",NormDoubleGaussSourceTF1,lowerlimit,upperlimit,4);
  TF1* fitDG_rstar = new TF1("fitDG_rstar",NormDoubleGaussSourceTF1,0.4,5.0,4);
  //fitDG_rstar->SetParameter(3,1);
  //fitDG_rstar->SetParLimits(3,0.8,1.0);//norm
  fitDG_rstar->FixParameter(3,1);

  //fitDG_rstar->SetParameter(0,0.36);//3.30932e-01
  //fitDG_rstar->SetParLimits(0,0.0,0.6);//first source
  fitDG_rstar->FixParameter(0,fit_rcore->GetParameter(1));//3.30932e-01

  fitDG_rstar->SetParameter(1,1.30);//1.23015e+00
  fitDG_rstar->SetParLimits(1,0.9,2.0);//second source

  //fitDG_rstar->SetParameter(2,0.24);//weight of first source
  //fitDG_rstar->SetParLimits(2,0.0,0.4);
  fitDG_rstar->FixParameter(2,FemtoPP);

  h_GhettoFemto_rstar->Fit(fitDG_rstar,"S, N, R, M");
  printf("DG = %.3f x %.3ffm + %.3f x %.3ffm\n",
    fitDG_rstar->GetParameter(2),
    fitDG_rstar->GetParameter(0),
    1.-fitDG_rstar->GetParameter(2),
    fitDG_rstar->GetParameter(1)
  );

  TGraph gRadKstar;
  gRadKstar.SetName("gRadKstar");
  gRadKstar.SetLineColor(kRed+1);
  gRadKstar.SetLineWidth(3);

  TGraph gMeanRadKstar;
  gMeanRadKstar.SetName("gMeanRadKstar");
  gMeanRadKstar.SetLineColor(kRed+1);
  gMeanRadKstar.SetLineWidth(3);


  //evaluate the radius based on the mean of the source distribution,
  //setting the actual r value to that of a Gauss source with the same mean
  TGraph gGhettoRadKstar;
  gGhettoRadKstar.SetName("gGhettoRadKstar");
  gGhettoRadKstar.SetLineColor(kRed+1);
  gGhettoRadKstar.SetLineWidth(3);

  TH1D** hkstar_rstar;
  hkstar_rstar = new TH1D* [h_Ghetto_kstar_rstar->GetXaxis()->GetNbins()];

  for(unsigned uMom=0; uMom<h_Ghetto_kstar_rstar->GetXaxis()->GetNbins(); uMom++){
    double kstar = h_Ghetto_kstar_rstar->GetXaxis()->GetBinCenter(uMom+1);
    hkstar_rstar[uMom] = NULL;
    if(kstar>780) continue;
    hkstar_rstar[uMom] = h_Ghetto_kstar_rstar->ProjectionY(TString::Format("hkstar_rstar_%.0f",kstar),uMom+1,uMom+1);
    hkstar_rstar[uMom]->Scale(1./hkstar_rstar[uMom]->Integral(),"width");
    GetCentralInterval(*hkstar_rstar[uMom], 0.9, lowerlimit, upperlimit, true);
    //printf("--- k* = %.0f ---\n",kstar);
    //printf(" The fit will be performed in the range [%.2f, %.2f] fm\n",lowerlimit,upperlimit);
    fSource->SetParameter(0,hkstar_rstar[uMom]->GetMean()/2.3);
    fSource->SetParLimits(0,hkstar_rstar[uMom]->GetMean()/4.,hkstar_rstar[uMom]->GetMean());
    hkstar_rstar[uMom]->Fit(fSource,"Q, S, N, R, M","",lowerlimit,upperlimit);
    //printf(" The effective Gaussian size is %.3f +/- %.3f fm, λ = %.3f +/- %.3f\n",
    //fSource->GetParameter(0),fSource->GetParError(0),fSource->GetParameter(1),fSource->GetParError(1));
    gRadKstar.SetPoint(uMom,kstar,fSource->GetParameter(0));
    gMeanRadKstar.SetPoint(uMom,kstar,hkstar_rstar[uMom]->GetMean());
    double gfm = GaussFromMean(hkstar_rstar[uMom]->GetMean());
    //printf(" The mean is %f; The correponding r0 is %f\n",hkstar_rstar[uMom]->GetMean(),gfm);
    gGhettoRadKstar.SetPoint(uMom,kstar,gfm);
  }

  h_Ghetto_PP_AngleRcP1->Scale(1./h_Ghetto_PP_AngleRcP1->Integral(),"width");
  h_Ghetto_PP_AngleRcP2->Scale(1./h_Ghetto_PP_AngleRcP2->Integral(),"width");
  h_Ghetto_PP_AngleP1P2->Scale(1./h_Ghetto_PP_AngleP1P2->Integral(),"width");
  h_Ghetto_RP_AngleRcP1->Scale(1./h_Ghetto_RP_AngleRcP1->Integral(),"width");
  h_Ghetto_PR_AngleRcP2->Scale(1./h_Ghetto_PR_AngleRcP2->Integral(),"width");
  h_Ghetto_RR_AngleRcP1->Scale(1./h_Ghetto_RR_AngleRcP1->Integral(),"width");
  h_Ghetto_RR_AngleRcP2->Scale(1./h_Ghetto_RR_AngleRcP2->Integral(),"width");
  h_Ghetto_RR_AngleP1P2->Scale(1./h_Ghetto_RR_AngleP1P2->Integral(),"width");






  fOutput.cd();
  //h_Ghetto_rstar->Write();
  //fit_rstar->Write();
  //h_Ghetto_rcore->Write();
  //fit_rcore->Write();

  TH1F* hAxisSource = new TH1F("hAxisSource","hAxisSource",128,0,8);
  hAxisSource->SetStats(false);
  hAxisSource->SetTitle("");
  hAxisSource->GetXaxis()->SetTitle("r (fm)");
  hAxisSource->GetXaxis()->SetTitleSize(0.06);
  hAxisSource->GetXaxis()->SetLabelSize(0.06);
  hAxisSource->GetXaxis()->SetTitleOffset(1.3);
  hAxisSource->GetXaxis()->SetLabelOffset(0.02);
  hAxisSource->GetYaxis()->SetTitle("4#pir*^{2}S(r) (1/fm)");
  hAxisSource->GetYaxis()->SetTitleSize(0.06);
  hAxisSource->GetYaxis()->SetLabelSize(0.06);
  hAxisSource->GetYaxis()->SetTitleOffset(1.00);
  hAxisSource->GetYaxis()->SetRangeUser(0.001,1.50);
  hAxisSource->Write();

  TH1F* hAxisMt = new TH1F("hAxisMt","hAxisMt",128,1.05,2.25);
  hAxisMt->SetStats(false);
  hAxisMt->SetTitle("");
  hAxisMt->GetXaxis()->SetTitle("m_{T} (GeV)");
  hAxisMt->GetXaxis()->SetTitleSize(0.06);
  hAxisMt->GetXaxis()->SetLabelSize(0.06);
  hAxisMt->GetXaxis()->SetTitleOffset(1.3);
  hAxisMt->GetXaxis()->SetLabelOffset(0.02);
  hAxisMt->GetYaxis()->SetTitle("r_{Gauss} (fm)");
  hAxisMt->GetYaxis()->SetTitleSize(0.06);
  hAxisMt->GetYaxis()->SetLabelSize(0.06);
  hAxisMt->GetYaxis()->SetTitleOffset(1.00);
  hAxisMt->GetYaxis()->SetRangeUser(0.4,1.85);
  hAxisMt->Write();

  h_Ghetto_kstar_rstar->GetXaxis()->SetRangeUser(0,1200);
  h_Ghetto_kstar_rstar_PP->GetXaxis()->SetRangeUser(0,1200);
  h_Ghetto_kstar_rstar_PR->GetXaxis()->SetRangeUser(0,1200);
  h_Ghetto_kstar_rstar_RP->GetXaxis()->SetRangeUser(0,1200);
  h_Ghetto_kstar_rstar_RR->GetXaxis()->SetRangeUser(0,1200);

  double UpRU = 100;
  if(type=="pK") UpRU = 40;
  h_Ghetto_kstar_rstar->GetYaxis()->SetRangeUser(0,UpRU);
  h_Ghetto_kstar_rstar_PP->GetYaxis()->SetRangeUser(0,UpRU);
  h_Ghetto_kstar_rstar_PR->GetYaxis()->SetRangeUser(0,UpRU);
  h_Ghetto_kstar_rstar_RP->GetYaxis()->SetRangeUser(0,UpRU);
  h_Ghetto_kstar_rstar_RR->GetYaxis()->SetRangeUser(0,UpRU);

  h_Ghetto_kstar_rstar->GetXaxis()->SetTitle("k* (MeV)");
  h_Ghetto_kstar_rstar_PP->GetXaxis()->SetTitle("k* (MeV)");
  h_Ghetto_kstar_rstar_PR->GetXaxis()->SetTitle("k* (MeV)");
  h_Ghetto_kstar_rstar_RP->GetXaxis()->SetTitle("k* (MeV)");
  h_Ghetto_kstar_rstar_RR->GetXaxis()->SetTitle("k* (MeV)");

  h_Ghetto_kstar_rstar->GetYaxis()->SetTitle("r* (fm)");
  h_Ghetto_kstar_rstar_PP->GetYaxis()->SetTitle("r* (fm)");
  h_Ghetto_kstar_rstar_PR->GetYaxis()->SetTitle("r* (fm)");
  h_Ghetto_kstar_rstar_RP->GetYaxis()->SetTitle("r* (fm)");
  h_Ghetto_kstar_rstar_RR->GetYaxis()->SetTitle("r* (fm)");

  TCanvas* cSource = new TCanvas("cSource", "cSource", 1);
  cSource->cd(0); cSource->SetCanvasSize(1280, 720); cSource->SetMargin(0.15,0.05,0.2,0.05);//lrbt
  cSource->SetLogy(true);
  hAxisSource->Draw("axis");
  h_GhettoFemto_rcore->Draw("same");
  h_GhettoFemto_rstar->Draw("same");
  fit_rcore->Draw("same");
  fit_rstar->Draw("same");


  TCanvas* cMt = new TCanvas("cMt", "cMt", 1);
  cMt->cd(0); cMt->SetCanvasSize(1280, 720); cMt->SetMargin(0.15,0.05,0.2,0.05);//lrbt
  hAxisMt->Draw("axis");
  g_GhettoFemto_mT_rcore_G.Draw("same");
  g_GhettoFemto_mT_rstar_G.Draw("same");

  fOutput.cd();
  cSource->Write();
  cMt->Write();
  //h_GhettoFemto_rstar->Write();
  fSource->Write();
  fitDG_rstar->Write();
  h_Ghetto_kstar->Write();
  h_Ghetto_kstar_rstar->Write();
  h_Ghetto_kstar_rstar_PP->Write();
  h_Ghetto_kstar_rstar_PR->Write();
  h_Ghetto_kstar_rstar_RP->Write();
  h_Ghetto_kstar_rstar_RR->Write();
  h_Ghetto_mT_rstar->Write();
  h_GhettoFemto_mT_rstar->Write();
  g_GhettoFemto_mT_rstar.Write();
  g_GhettoFemto_mT_rstar_G.Write();
  g_GhettoFemto_mT_rcore.Write();
  g_GhettoFemto_mT_rcore_G.Write();
  h_GhettoFemto_mT_kstar->Write();
  h_Ghetto_mT_costh->Write();
//return;
  h_GhettoSP_pT_th->Write();
  h_GhettoSP_pT_1->Write();
  h_GhettoSP_pT_2->Write();
  h_Ghetto_PP_AngleRcP1->Write();
  h_Ghetto_PP_AngleRcP2->Write();
  h_Ghetto_PP_AngleP1P2->Write();
  h_Ghetto_RP_AngleRcP1->Write();
  h_Ghetto_PR_AngleRcP2->Write();
  h_Ghetto_RR_AngleRcP1->Write();
  h_Ghetto_RR_AngleRcP2->Write();
  h_Ghetto_RR_AngleP1P2->Write();
  gRadKstar.Write();
  gMeanRadKstar.Write();
  gGhettoRadKstar.Write();
  h_GhettoSPr_X->Write();
  h_GhettoSPr_Y->Write();
  h_GhettoSPr_Z->Write();
  h_GhettoSPr_Rho->Write();
  h_GhettoSPr_R->Write();
  h_GhettoSP_X->Write();
  h_GhettoSP_Y->Write();
  h_GhettoSP_Z->Write();
  h_GhettoSP_Rho->Write();
  h_GhettoSP_R->Write();
  for(unsigned uMom=0; uMom<h_Ghetto_kstar_rstar->GetXaxis()->GetNbins(); uMom++){
    if(hkstar_rstar[uMom]){
      hkstar_rstar[uMom]->Write();
      delete hkstar_rstar[uMom];
    }
  }


  delete h_pT_p_all;
  delete h_pT_d_all;
  delete hSampleQA_p;

  //delete h_Ghetto_rstar;
  delete h_GhettoFemto_rstar;
  delete fSource;
  delete h_Ghetto_kstar;
  delete h_Ghetto_mT_rstar;
  delete h_GhettoFemto_mT_rstar;
  //delete g_GhettoFemto_mT_rstar;
  delete h_GhettoFemto_mT_kstar;
  delete h_Ghetto_mT_costh;
  delete h_GhettoSP_pT_th;
  delete h_GhettoSP_pT_1;
  delete h_GhettoSP_pT_2;
  delete h_Ghetto_PP_AngleRcP1;
  delete h_Ghetto_PP_AngleRcP2;
  delete h_Ghetto_PP_AngleP1P2;
  delete h_Ghetto_RP_AngleRcP1;
  delete h_Ghetto_PR_AngleRcP2;
  delete h_Ghetto_RR_AngleRcP1;
  delete h_Ghetto_RR_AngleRcP2;
  delete h_Ghetto_RR_AngleP1P2;
  delete h_Ghetto_kstar_rstar;
  delete h_GhettoSPr_X;
  delete h_GhettoSPr_Y;
  delete h_GhettoSPr_Z;
  delete h_GhettoSPr_Rho;
  delete h_GhettoSPr_R;
  delete h_GhettoSP_X;
  delete h_GhettoSP_Y;
  delete h_GhettoSP_Z;
  delete h_GhettoSP_Rho;
  delete h_GhettoSP_R;
  delete h_GhettoFemto_rcore;
  delete h_GhettoFemto_mT_rcore;
  delete fit_rstar;
  delete fit_rcore;
  delete fitDG_rstar;
  delete cSource;
  delete cMt;
*/
}






int CECA_PAPER(int argc, char *argv[]){

  //how to read/write the Levy pars into a file
  //TestSaveStuctToFile();

  //a test that extrapolation of Levy pars works
  //TestDoubleSourceOperation();

  Test_pp_Statistics_1();


  return 0;
}
