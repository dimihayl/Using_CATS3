
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
#include <stdio.h>
#include <stdlib.h>

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
#include "TLegend.h"

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




DLM_Histo<float>* GetPtEta(TString FileNameP, TString FileNameAP,
  TString HistoNameP="pTDist_after", TString HistoNameAP="pTDist_after", const double EtaCut = 0.8){

  TH1F* h_pT_p;
  TH1F* h_pT_ap;
  TH1F* h_pT_all;

  TFile file_p(FileNameP,"read");
  h_pT_p = (TH1F*)file_p.Get(HistoNameP);
  if(!h_pT_p) printf("ISSUE with h_pT_p\n");
  gROOT->cd();
  h_pT_all = (TH1F*)h_pT_p->Clone("h_pT_all");

  TFile file_ap(FileNameAP,"read");
  h_pT_ap = (TH1F*)file_ap.Get(HistoNameAP);
  if(!h_pT_ap) printf("ISSUE with h_pT_ap\n");
  h_pT_all->Add(h_pT_ap);

  DLM_Histo<float>* dlm_pT_p = Convert_TH1F_DlmHisto(h_pT_all);
  dlm_pT_p->RescaleAxis(0,1000,false);

  double* BinRange = NULL;
  double axis[2];

  DLM_Histo<float>* dlm_pT_eta = new DLM_Histo<float>();
  dlm_pT_eta->SetUp(2);
  BinRange = dlm_pT_p->GetBinRange(0);
  dlm_pT_eta->SetUp(0,dlm_pT_p->GetNbins(),BinRange);
  delete [] BinRange;
  dlm_pT_eta->SetUp(1,1,-EtaCut,EtaCut);
  dlm_pT_eta->Initialize();
  for(unsigned uBin=0; uBin<dlm_pT_p->GetNbins(); uBin++){
    dlm_pT_eta->SetBinContent(uBin,0,dlm_pT_p->GetBinContent(uBin));
  }

  file_p.Close();
  file_ap.Close();
  delete h_pT_all;
  return dlm_pT_eta;
}







//parameters to control:
//MUST:
//  SEED, GLOB_TIMEOUT, multiplicity, target_yield, femto_region,
//  d_x, d_y, d_z, h_x, h_y, h_z, h_fct, tau, tau_prp, tau_fct;
//  hdr_size, hdr_slope, th_kick, frag_beta, fixed_hdr;
//  a flag for the momentum distribution;
//  a flag for resonance variations;
//  a flag for the type (1 is pp, 2 is pL)
//  WILD_CARD flag for the future
//  24 numbers in total -> input from file
//
// OUTPUT:
//  DLM_Histo: Ghetto_kstar_rstar_mT
//  A settings file (perhaps not, integrate as input), which will again be a histo (or a custum class?),
//  with the values of all parameters above.
//  the names of the file should be unique, e.g. SEED id or whatever

//QUESTIONS: do we care about h_xyz differentially
// -> it leads to 1.5% error -> take as systematics :D
//      actually, if a take the hz = 2x hT, than we have +4.3% on the radius, but the scaling remains the same
//      i.e. we can comment that the choice of hz can simply scale the thing up and down by a small amount



//the names should be given without extension. They should also have the FULL path!!
//the assumed extension is *.txt for the Input and .dlm.hst for the Output
//WILD_CARD: if -1: the mT binning is super fine, done for testing (figure out what binning to use)
int Ceca_pp_or_pL(const TString FileBase, const TString InputFolder, const TString OutputFolder, const TString LogFolder,
                  const int ParID, const int JobID, const int NumCPU){

  printf("FileBase = %s\n",FileBase.Data());
  printf("InputFolder = %s\n",InputFolder.Data());
  printf("OutputFolder = %s\n",OutputFolder.Data());
  printf("LogFolder = %s\n",LogFolder.Data());
  printf("ParID = %i\n",ParID);
  printf("JobID = %i\n",JobID);
  printf("NumCPU = %i\n",NumCPU);

  const double TIMEOUT = 30;
  const double EtaCut = 0.8;
  const bool PROTON_RESO = true;
  const bool EQUALIZE_TAU = true;

  int SEED = (ParID+1)*10. + JobID;

  TString InputFileName = InputFolder+FileBase+TString::Format(".%i.%i.dlm.job",ParID,JobID);
  TString OutputFileNameFull = OutputFolder+FileBase+TString::Format(".%i.%i.full.dlm.hst",ParID,JobID);
  TString OutputFileNameCore = OutputFolder+FileBase+TString::Format(".%i.%i.core.dlm.hst",ParID,JobID);
  //a binary file that contains information on the statistics we have collected so far
  //this program will search for that file, and if it exists it will read it and add the current yield
  //to the total yield
  TString LogFileName = LogFolder+FileBase+TString::Format(".%i.%i.dlm.log",ParID,JobID);

  printf("\n");
  printf("InputFileName = %s\n",InputFileName.Data());
  printf("OutputFileNameFull = %s\n",OutputFileNameFull.Data());
  printf("OutputFileNameCore = %s\n",OutputFileNameCore.Data());
  printf("LogFileName = %s\n",LogFileName.Data());

  TREPNI Database(0);
  Database.SetSeed(11);
  std::vector<TreParticle*> ParticleList;
  ParticleList.push_back(Database.NewParticle("Proton"));
  ParticleList.push_back(Database.NewParticle("ProtonReso"));
  ParticleList.push_back(Database.NewParticle("Lambda"));
  ParticleList.push_back(Database.NewParticle("LambdaReso"));
  ParticleList.push_back(Database.NewParticle("Pion"));


  //DLM_Histo<double> SettingsHisto;
  //SettingsHisto.QuickLoad(InputFileName.Data());

  double GLOB_TIMEOUT = 0;
  unsigned multiplicity = 0;
  unsigned target_yield = 0;
  float femto_region = 100;
  float d_x = 0;
  float d_y = 0;
  float d_z = 0;
  float h_x = 0;
  float h_y = 0;
  float h_z = 0;
  float h_fct = 0;
  float tau = 0;
  float tau_fct = 0;
  bool tau_prp = true;
  float hdr_size = 0;
  float hdr_slope = 0;
  float th_kick = 0;
  float frag_beta = 0;
  float fixed_hdr = 1;
  int momdst_flag = 1;
  int reso_flag = 1;//default, read out masses etc
  //int type_flag = 0;
  int wildcard_flag = 0;
  TString type = "";
  TString AnaVersion = "";

  double m_proton_reso;
  double frac_proton_reso;
  double tau_proton_reso;

  double m_lambda_reso;
  double frac_lambda_reso;
  double tau_lambda_reso;

  char* cline = new char [512];
  char* cdscr = new char [128];
  char* cval = new char [128];
  double read_value;
  FILE *InFile;
  InFile = fopen(InputFileName.Data(), "r");
  if(!InFile){
      printf("\033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", InputFileName.Data());
      return 0;
  }
  fseek ( InFile , 0 , SEEK_END );
  long EndPos;
  EndPos = ftell (InFile);
  fseek ( InFile , 0 , SEEK_SET );
  long CurPos;
  while(!feof(InFile)){
    if(!fgets(cline, 511, InFile)){
      //printf("\033[1;31mERROR:\033[0m The file\033[0m %s cannot be properly read (%s)!\n", InputFileName.Data(),cline);
    }
    sscanf(cline, "%s %s",cdscr,cval);
    if(strcmp(cdscr,"type")==0){
      type = TString(cval);
    }
    else if(strcmp(cdscr,"AnaVersion")==0){
      AnaVersion = TString(cval);
    }
    else{
      read_value = stod(cval);
            if(strcmp(cdscr,"GLOB_TIMEOUT")==0) {GLOB_TIMEOUT = read_value;}
      else if(strcmp(cdscr,"multiplicity")==0) {multiplicity = unsigned(read_value);}
      else if(strcmp(cdscr,"target_yield")==0) {target_yield = unsigned(read_value);}
      else if(strcmp(cdscr,"femto_region")==0) {femto_region = read_value;}
      else if(strcmp(cdscr,"d_x")==0) {d_x = read_value;}
      else if(strcmp(cdscr,"d_y")==0) {d_y = read_value;}
      else if(strcmp(cdscr,"d_z")==0) {d_z = read_value;}
      else if(strcmp(cdscr,"h_x")==0) {h_x = read_value;}
      else if(strcmp(cdscr,"h_y")==0) {h_y = read_value;}
      else if(strcmp(cdscr,"h_z")==0) {h_z = read_value;}
      else if(strcmp(cdscr,"h_fct")==0) {h_fct = read_value;}
      else if(strcmp(cdscr,"tau")==0) {tau = read_value;}
      else if(strcmp(cdscr,"tau_fct")==0) {tau_fct = read_value;}
      else if(strcmp(cdscr,"tau_prp")==0) {tau_prp = bool(read_value);}
      else if(strcmp(cdscr,"hdr_size")==0) {hdr_size = read_value;}
      else if(strcmp(cdscr,"hdr_slope")==0) {hdr_slope = read_value;}
      else if(strcmp(cdscr,"th_kick")==0) {th_kick = read_value;}
      else if(strcmp(cdscr,"frag_beta")==0) {frag_beta = read_value;}
      else if(strcmp(cdscr,"fixed_hdr")==0) {fixed_hdr = read_value;}
      else if(strcmp(cdscr,"momdst_flag")==0) {momdst_flag = int(read_value);}
      else if(strcmp(cdscr,"reso_flag")==0) {reso_flag = int(read_value);}
      else if(strcmp(cdscr,"wildcard_flag")==0) {wildcard_flag = int(read_value);}
      else if(strcmp(cdscr,"m_proton_reso")==0) {m_proton_reso = read_value;}
      else if(strcmp(cdscr,"tau_proton_reso")==0) {tau_proton_reso = read_value;}
      else if(strcmp(cdscr,"frac_proton_reso")==0) {frac_proton_reso = read_value;}
      else if(strcmp(cdscr,"m_lambda_reso")==0) {m_lambda_reso = read_value;}
      else if(strcmp(cdscr,"tau_lambda_reso")==0) {tau_lambda_reso = read_value;}
      else if(strcmp(cdscr,"frac_lambda_reso")==0) {frac_lambda_reso = read_value;}
    }
  }
  delete [] cval;
  delete [] cline;
  delete [] cdscr;

  printf("AnaVersion = %s\n",AnaVersion.Data());
  printf("GLOB_TIMEOUT = %f\n",GLOB_TIMEOUT);
  printf("multiplicity = %u\n",multiplicity);
  printf("target_yield = %u\n",target_yield);
  //printf("current_yield = %u\n",current_yield);
  printf("femto_region = %f\n",femto_region);
  printf("d_x = %f\n",d_x);
  printf("d_y = %f\n",d_y);
  printf("d_z = %f\n",d_z);
  printf("h_x = %f\n",h_x);
  printf("h_y = %f\n",h_y);
  printf("h_z = %f\n",h_z);
  printf("h_fct = %f\n",h_fct);
  printf("tau = %f\n",tau);
  printf("tau_fct = %f\n",tau_fct);
  printf("tau_prp = %i\n",tau_prp);
  printf("hdr_size = %f\n",hdr_size);
  printf("hdr_slope = %f\n",hdr_slope);
  printf("th_kick = %f\n",th_kick);
  printf("frag_beta = %f\n",frag_beta);
  printf("fixed_hdr = %f\n",fixed_hdr);
  printf("momdst_flag = %i\n",momdst_flag);
  printf("reso_flag = %i\n",reso_flag);
  printf("wildcard_flag = %i\n",wildcard_flag);
  printf("frac_proton_reso = %f\n",frac_proton_reso);
  printf("m_proton_reso = %f\n",m_proton_reso);
  printf("tau_proton_reso = %f\n",tau_proton_reso);
  printf("frac_lambda_reso = %f\n",frac_lambda_reso);
  printf("m_lambda_reso = %f\n",m_lambda_reso);
  printf("tau_lambda_reso = %f\n",tau_lambda_reso);
  printf("type = %s\n",type.Data());

  //if(current_yield>=target_yield){
  //  printf("\033[1;31mERROR:\033[0m current_yield>=target_yield, this should NOT happen!\n");
  //  return 0;
  //}

//return;

  //here use the flags if needed
  //frac_proton_reso = 64.22;
  //m_proton_reso = 1362;
  //tau_proton_reso = 1.65;

  //frac_lambda_reso = 64.38;
  //m_lambda_reso = 1463;
  //tau_lambda_reso = 4.69;

  //basic QA: too lazy to doo it
  if(multiplicity<1){
    printf("Bad multiplicity!\n");
    return 0;
  }
  if(multiplicity==1){
    printf("Potenitally bad multiplicity (1)!\n");
  }

  if(type!="pp"&&type!="pL"){
    printf("WHAT IS THIS: type = %s\n", type.Data());
    return 0;
  }

  DLM_Histo<float>* dlm_pT_eta_p = GetPtEta(
    TString::Format("%s/Jaime/p_pT.root",GetCernBoxDimi()),
    TString::Format("%s/Jaime/ap_pT.root",GetCernBoxDimi()),
    "pTDist_after", "pTDist_after", EtaCut);
  DLM_Histo<float>* dlm_pT_eta_L = GetPtEta(
    TString::Format("%s/Jaime/L_pT.root",GetCernBoxDimi()),
    TString::Format("%s/Jaime/aL_pT.root",GetCernBoxDimi()),
    "pTDist_after", "pTDist_after", EtaCut);

  for(TreParticle* prt : ParticleList){
    if(prt->GetName()=="Proton"){
      prt->SetMass(Mass_p);
      prt->SetAbundance(100.-frac_proton_reso);
      prt->SetRadius(hdr_size);
      prt->SetRadiusSlope(hdr_slope);
      if(dlm_pT_eta_p) prt->SetPtEtaPhi(*dlm_pT_eta_p);
      else prt->SetPtPz(0.85*prt->GetMass(),0.85*prt->GetMass());
    }
    else if(prt->GetName()=="ProtonReso"){
      prt->SetMass(m_proton_reso);
      prt->SetAbundance(frac_proton_reso);
      prt->SetWidth(hbarc/tau_proton_reso);

      prt->NewDecay();
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Proton"));
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
      prt->GetDecay(0)->SetBranching(100);

      if(dlm_pT_eta_p) prt->SetPtEtaPhi(*dlm_pT_eta_p);
      else prt->SetPtPz(0.85*prt->GetMass(),0.85*prt->GetMass());
    }
    else if(prt->GetName()=="Lambda"){
      prt->SetMass(Mass_L);
      prt->SetAbundance(100.-frac_lambda_reso);
      prt->SetRadius(hdr_size);
      prt->SetRadiusSlope(hdr_slope);
      if(dlm_pT_eta_L) prt->SetPtEtaPhi(*dlm_pT_eta_L);
      else prt->SetPtPz(0.85*prt->GetMass(),0.85*prt->GetMass());
    }
    else if(prt->GetName()=="LambdaReso"){
      prt->SetMass(m_lambda_reso);
      prt->SetAbundance(frac_lambda_reso);
      prt->SetWidth(hbarc/tau_lambda_reso);

      prt->NewDecay();
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Lambda"));
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
      prt->GetDecay(0)->SetBranching(100);

      if(dlm_pT_eta_L) prt->SetPtEtaPhi(*dlm_pT_eta_L);
      else prt->SetPtPz(0.85*prt->GetMass(),0.85*prt->GetMass());
    }
    else if(prt->GetName()=="Pion"){
      prt->SetMass(Mass_pic);
      prt->SetAbundance(0);
      prt->SetRadius(hdr_size);
      prt->SetRadiusSlope(hdr_slope);
    }
  }//ParticleList

  std::vector<std::string> ListOfParticles;
  if(type=="pp"){
    ListOfParticles.push_back("Proton");
    ListOfParticles.push_back("Proton");
  }
  else if(type=="pL"){
    ListOfParticles.push_back("Proton");
    ListOfParticles.push_back("Lambda");
  }
  else{
    printf("Issue with the type!\n");
    return 0;
  }

  CECA Ivana(Database,ListOfParticles);

  Ivana.SetDisplacementX(d_x);
  Ivana.SetDisplacementY(d_y);
  Ivana.SetDisplacementZ(d_z);

  Ivana.SetHadronizationX(h_x);
  Ivana.SetHadronizationY(h_y);
  Ivana.SetHadronizationZ(h_z);
  Ivana.SetHadrFluctuation(h_fct);

  Ivana.SetTau(tau,tau_prp);
  Ivana.SetTauFluct(tau_fct);

  Ivana.SetThermalKick(th_kick);
  Ivana.SetFixedHadr(fixed_hdr);
  Ivana.SetFragmentBeta(frag_beta);

  Ivana.SetTargetStatistics(target_yield);
  Ivana.SetEventMult(multiplicity);
  Ivana.SetSourceDim(2);
  Ivana.SetDebugMode(false);

  Ivana.SetThreadTimeout(TIMEOUT);
  Ivana.SetGlobalTimeout(GLOB_TIMEOUT);
  Ivana.EqualizeFsiTime(true);
  Ivana.SetFemtoRegion(femto_region);

  Ivana.GHETTO_EVENT = true;

  if(type=="pp"){
    if(wildcard_flag==-1){
      //Ivana.Ghetto_NumMtBins = 206;
      //Ivana.Ghetto_MtBins = new double [Ivana.Ghetto_NumMtBins+1];
      //for(unsigned uMt=0; uMt<=156; uMt++){
      //  Ivana.Ghetto_MtBins[uMt] = 940. + double(uMt)*10;
      //}
      //for(unsigned uMt=157; uMt<=Ivana.Ghetto_NumMtBins; uMt++){
      //  Ivana.Ghetto_MtBins[uMt] = 2500. + double(uMt-156)*50;
      //}
      Ivana.Ghetto_NumMtBins = 360;
      Ivana.Ghetto_MtBins = new double [Ivana.Ghetto_NumMtBins+1];
      for(unsigned uMt=0; uMt<=360; uMt++){
        Ivana.Ghetto_MtBins[uMt] = 940. + double(uMt)*10;
      }

      Ivana.Ghetto_NumMomBins = 25;
      Ivana.Ghetto_MomMin = 0;
      Ivana.Ghetto_MomMax = 100;

      Ivana.Ghetto_NumRadBins = 192;//twice the bin width compared to default
      Ivana.Ghetto_RadMin = 0;
      Ivana.Ghetto_RadMax = 48;
    }
    else{
      Ivana.Ghetto_NumMtBins = 10;
      Ivana.Ghetto_MtBins = new double [Ivana.Ghetto_NumMtBins+1];
      Ivana.Ghetto_MtBins[0] = 930; //avg  983 ( 985)
      Ivana.Ghetto_MtBins[1] = 1020;//avg 1054 (1055)
      Ivana.Ghetto_MtBins[2] = 1080;//avg 1110 (1110)
      Ivana.Ghetto_MtBins[3] = 1140;//avg 1168 (1170)
      Ivana.Ghetto_MtBins[4] = 1200;//avg 1228 (1230)
      Ivana.Ghetto_MtBins[5] = 1260;//avg 1315 (1315)
      Ivana.Ghetto_MtBins[6] = 1380;//avg 1463 (1460)
      Ivana.Ghetto_MtBins[7] = 1570;//avg 1681 (1680)
      Ivana.Ghetto_MtBins[8] = 1840;//avg 1923 (1920)
      Ivana.Ghetto_MtBins[9] = 2030;//avg 2303 (2300)
      Ivana.Ghetto_MtBins[10] = 4500;

      Ivana.Ghetto_NumMomBins = 150;
      Ivana.Ghetto_MomMin = 0;
      Ivana.Ghetto_MomMax = 600;
    }

  }
  else if(type=="pL"){
//NOT_DONE_YET
    if(wildcard_flag==-1){
      Ivana.Ghetto_NumMtBins = 360;
      Ivana.Ghetto_MtBins = new double [Ivana.Ghetto_NumMtBins+1];
      for(unsigned uMt=0; uMt<=360; uMt++){
        Ivana.Ghetto_MtBins[uMt] = 1000. + double(uMt)*10;
      }

      Ivana.Ghetto_NumMomBins = 25;
      Ivana.Ghetto_MomMin = 0;
      Ivana.Ghetto_MomMax = 100;

      Ivana.Ghetto_NumRadBins = 192;//twice the bin width compared to default
      Ivana.Ghetto_RadMin = 0;
      Ivana.Ghetto_RadMax = 48;
    }
    else{
      Ivana.Ghetto_NumMtBins = 8;
      Ivana.Ghetto_MtBins = new double [Ivana.Ghetto_NumMtBins+1];
      Ivana.Ghetto_MtBins[0] = 1000;//avg 1121 (1120)
      Ivana.Ghetto_MtBins[1] = 1170;//avg 1210 (1210)
      Ivana.Ghetto_MtBins[2] = 1250;//avg 1288 (1290)
      Ivana.Ghetto_MtBins[3] = 1330;//avg 1377 (1380)
      Ivana.Ghetto_MtBins[4] = 1430;//avg 1536 (1540)
      Ivana.Ghetto_MtBins[5] = 1680;//avg 1753 (1750)
      Ivana.Ghetto_MtBins[6] = 1840;//avg 1935 (1935)
      Ivana.Ghetto_MtBins[7] = 2060;//avg 2334 (2330)
      Ivana.Ghetto_MtBins[8] = 4800;

      Ivana.Ghetto_NumMomBins = 150;
      Ivana.Ghetto_MomMin = 0;
      Ivana.Ghetto_MomMax = 600;
    }
  }

  if(NumCPU>1){
    Ivana.SetDebugMode(true);
    for(unsigned uTh=0; uTh<NumCPU; uTh++){
      Ivana.SetSeed(uTh,SEED*(NumCPU)+uTh);
    }
  }
  else{
    Ivana.SetDebugMode(false);
    Ivana.SetSeed(0,SEED);
  }

  Ivana.GoBabyGo(NumCPU);

  Ivana.Ghetto_kstar_rstar_mT->QuickWrite(OutputFileNameFull,true);
  Ivana.Ghetto_kstar_rcore_mT->QuickWrite(OutputFileNameCore,true);

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

  ofstream log_file (LogFileName.Data(), ios::out);
  if(log_file.fail()){
    printf("\033[1;31mERROR:\033[0m The file\033[0m The LOG file %s cannot be created!\n", LogFileName.Data());
  }
  else{
    log_file << Ivana.GetStatistics();
    log_file.close();
  }


  delete dlm_pT_eta_p;
  delete dlm_pT_eta_L;
  //printf("Returning %i\n",int(Ivana.GetStatistics()));
  return Ivana.GetStatistics();
}


//project the output into 2D histos of kstar_vs_rstar for each mT bin
void ReadDlmHst(){
  //system("echo -n '1. Current Directory is '; pwd");
  //TString HistoFileName = TString::Format("%s/CECA_Paper/dadd_f/testout1.dlm.hst",GetFemtoOutputFolder());
  //TString HistoFileName = "./Output/CECA_Paper/dadd_f/Output/TEST1_3_full.dlm.hst";
  TString HistoFileName = "./Output/CECA_Paper/Ceca_pp_or_pL/TEST3_pp_full.dlm.hst";

  DLM_Histo<float> dlmHisto;
  dlmHisto.QuickLoad(HistoFileName.Data());

  const unsigned NumMt = dlmHisto.GetNbins(2);
  const unsigned NumRad = dlmHisto.GetNbins(1);
  const unsigned NumMom = dlmHisto.GetNbins(0);

  double* BinRangeMt = dlmHisto.GetBinRange(2);
  double* BinRangeRad = dlmHisto.GetBinRange(1);
  double* BinRangeMom = dlmHisto.GetBinRange(0);


  //get how many entries we have per mT_kstar bin to build our source
  TH2F* hYield = new TH2F("hYield","hYield",NumMom,BinRangeMom,NumMt,BinRangeMt);
  TH2F** h_kstar_rstar = new TH2F* [NumMt];
  for(unsigned uMt=0; uMt<NumMt; uMt++){
    TString hName = TString::Format("h_kstar_rstar_mT%u",uMt);
    h_kstar_rstar[uMt] = new TH2F(hName,hName,NumMom,BinRangeMom,NumRad,BinRangeRad);
    for(unsigned uMom=0; uMom<NumMom; uMom++){
      hYield->SetBinContent(uMom+1,uMt+1,0);
      for(unsigned uRad=0; uRad<NumRad; uRad++){
        h_kstar_rstar[uMt]->SetBinContent(uMom+1,uRad+1,dlmHisto.GetBinContent(uMom,uRad,uMt));
        hYield->SetBinContent(uMom+1,uMt+1,hYield->GetBinContent(uMom+1,uMt+1)+dlmHisto.GetBinContent(uMom,uRad,uMt));
      }
    }
    h_kstar_rstar[uMt]->Sumw2();
    //printf("Integral %u: %u\n",uMt,unsigned(h_kstar_rstar[uMt]->Integral()));
  }
  hYield->Sumw2();

  //TFile fOutput(TString::Format("%s/CECA_Paper/dadd_f/testout1.root",GetFemtoOutputFolder()),"recreate");
  //TFile fOutput(TString::Format("%s/CECA_Paper/dadd_f/Output/TEST1_3_full.root",GetFemtoOutputFolder()),"recreate");
  TFile fOutput(TString::Format("%s/CECA_Paper/Ceca_pp_or_pL/TEST3_pp_full.root",GetFemtoOutputFolder()),"recreate");
  for(unsigned uMt=0; uMt<NumMt; uMt++){
    h_kstar_rstar[uMt]->Write();
  }
  hYield->Write();

  for(unsigned uMt=0; uMt<NumMt; uMt++){
    delete h_kstar_rstar[uMt];
  }
  delete [] h_kstar_rstar;
  delete hYield;
  delete [] BinRangeRad;
  delete [] BinRangeMom;
  delete [] BinRangeMt;

}

//what is the optimal binning in Ceca
void OptimalBinning(TString syst, TString InputFile){

  double* BinCenter;
  //the lower acceptable limit
  double* BinLowLimit;
  //the upper acceptabe limit
  double* BinUpLimit;
  //the current averaged bin center based on these limits
  double* BinAvg;
  //

  //BinCenter_pp[0] = 1.1077;
  //BinCenter_pp[1] = 1.1683;
  //BinCenter_pp[2] = 1.2284;
  //BinCenter_pp[3] = 1.3156;
  //BinCenter_pp[4] = 1.4628;
  //BinCenter_pp[5] = 1.6872;
  //BinCenter_pp[6] = 2.2116;

  //BinCenter_pL[0] = 1.2124;
  //BinCenter_pL[1] = 1.2896;
  //BinCenter_pL[2] = 1.376;
  //BinCenter_pL[3] = 1.5407;
  //BinCenter_pL[4] = 1.756;
  //BinCenter_pL[5] = 2.2594;


  //read the histo, integral to get total yield
  //get optimal yield and start expanding the bin range to accomodate it
  //if you cannot achieve it, getting up to the bin center of next bin (say within some pre-defined limit), reduce the optimal yield 2x try again
  //

  //unsigned TotalYield = ;

}



//hadd but for my histos (n.b. following bash 0 is okay, 1 is error)
int dadd_f(int argc, char *argv[]){
  //printf("Hello\n");
  //system("echo -n '1. Current Directory is '; pwd");
  if(argc<3) {printf("\033[1;31mERROR:\033[0m Not enough arguments [dadd Output Input]\n"); return 1;}
  //1 default, 1 for the output
  const int NumInputFiles = argc-2;
  printf("Attempting to merge %i number of files\n",NumInputFiles);
  DLM_Histo<float> MergedHisto;
  //DLM_Histo<float>* HISTO;
  //DLM_Histo<float>& uHisto = HISTO[0];
  printf("Final file: %s\n",argv[1]);
  int Status = 0;
  for(unsigned uIn=0; uIn<NumInputFiles; uIn++){
    DLM_Histo<float> uHisto;
    printf(" -> %s\n",argv[2+uIn]);
    if(!uHisto.QuickLoad(argv[2+uIn])) return 1;
    printf("   done\n");
    if(uIn==0){
      Status = !(MergedHisto = uHisto);
    }
    else{
      Status = !(MergedHisto += uHisto);
    }

    //delete HISTO;
    if(Status) return Status;
  }
  MergedHisto.QuickWrite(argv[1],true);
  return 0;
}

//project the output into 2D histos of kstar_vs_rstar for each mT bin
//the input args should be the dlm.hst input file and the .root output file
int dlmhst_root(TString DADD_OUT, TString ROOT_OUT){

  DLM_Histo<float> dlmHisto;
  if(!dlmHisto.QuickLoad(DADD_OUT.Data())) return 1;

  const unsigned NumMt = dlmHisto.GetNbins(2);
  const unsigned NumRad = dlmHisto.GetNbins(1);
  const unsigned NumMom = dlmHisto.GetNbins(0);

  double* BinRangeMt = dlmHisto.GetBinRange(2);
  double* BinRangeRad = dlmHisto.GetBinRange(1);
  double* BinRangeMom = dlmHisto.GetBinRange(0);


  //get how many entries we have per mT_kstar bin to build our source
  TH2F* hYield = new TH2F("hYield","hYield",NumMom,BinRangeMom,NumMt,BinRangeMt);
  TH2F** h_kstar_rstar = new TH2F* [NumMt];
  for(unsigned uMt=0; uMt<NumMt; uMt++){
    TString hName = TString::Format("h_kstar_rstar_mT%u",uMt);
    h_kstar_rstar[uMt] = new TH2F(hName,hName,NumMom,BinRangeMom,NumRad,BinRangeRad);
    for(unsigned uMom=0; uMom<NumMom; uMom++){
      hYield->SetBinContent(uMom+1,uMt+1,0);
      for(unsigned uRad=0; uRad<NumRad; uRad++){
        h_kstar_rstar[uMt]->SetBinContent(uMom+1,uRad+1,dlmHisto.GetBinContent(uMom,uRad,uMt));
        hYield->SetBinContent(uMom+1,uMt+1,hYield->GetBinContent(uMom+1,uMt+1)+dlmHisto.GetBinContent(uMom,uRad,uMt));
      }
    }
    h_kstar_rstar[uMt]->Sumw2();
  }
  hYield->Sumw2();

  TFile fOutput(ROOT_OUT,"recreate");
  if(!fOutput.IsOpen()) return 1;
  for(unsigned uMt=0; uMt<NumMt; uMt++){
    h_kstar_rstar[uMt]->Write();
  }
  hYield->Write();

  for(unsigned uMt=0; uMt<NumMt; uMt++){
    delete h_kstar_rstar[uMt];
  }
  delete [] h_kstar_rstar;
  delete hYield;
  delete [] BinRangeRad;
  delete [] BinRangeMom;
  delete [] BinRangeMt;
  //printf("true\n");
  return 0;
}

int dlmhst_root(int argc, char *argv[]){
  if(argc!=3) {printf("\033[1;31mERROR:\033[0m Wrong args [DlmHst_Root Input Output]\n"); return 1;}
  TString DADD_OUT = TString(argv[1]);
  TString ROOT_OUT = TString(argv[2]);
  return dlmhst_root(DADD_OUT,ROOT_OUT);
}






//FitMode:
// 1 - core (single non-shifted, fully normalized, Gauss). Saves the pars within the 1st Gauss in SrcPar
// 2 - full: takes the pars of the first Gauss in SrcPar to fix the first Gaussian,
//            and fits by adding a second Gaussian that can be shifted
//            G3 is OFF!! if we are stuck, we change a bit the core gauss
TF1* Fit_wc0(int FitMode, TH1F* hSrc, TriGauss& SrcPar, double& Dist_Max, double& Chi2){
  double lowerlimit, upperlimit;
  GetCentralInterval(*hSrc, 0.98, lowerlimit, upperlimit, true);
//printf("SrcPar.shift1 = %f\n",SrcPar.shift1);
//usleep(100e3);
  if(lowerlimit>5) lowerlimit = 5;
  if(upperlimit>20) upperlimit = 20;

  TF1* fSrc = new TF1(TString::Format("fSrc"),NormTripleShiftedGaussTF1,0.,20.,9);

  //like chi2, but not normalized to error
  double Dist2 = 0;
  double NDPts = 0;
  Chi2 = 0;
  //nsig = 0;
  double NDPts_chi2 = 0;
  //double Dist_Max = 0;
  Dist_Max = 0;
  double Nsig_AtDistMax = 0;
  double Rad_AtDistMax = 0;

  const double Dist_Limit = 0.015;
  const double NsigMax_Limit = 3.0;
  const double Chi2_Limit = 2.0;
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


    //if we are very stuck, we allow to change the norm
    if(ResetCount<4) fSrc->FixParameter(0,1);
    else{
      fSrc->SetParameter(0,rangen.Uniform(0.9,1.0));
      fSrc->SetParLimits(0,0.9,1.0);
    }

    if(FitMode==1){
      //G1
      fSrc->SetParameter(1,rangen.Uniform(0.0,1.0));
      fSrc->SetParLimits(1,0,3.);

      fSrc->FixParameter(2,0);

      fSrc->FixParameter(3,1);

      fSrc->FixParameter(4,1);
      fSrc->FixParameter(5,0);
      fSrc->FixParameter(6,1);

      fSrc->FixParameter(7,1);
      fSrc->FixParameter(8,0);
    }
    else if(FitMode==2){
      //if we are stuck, slight variation of the G1 are allowed
      //in steps: first weight, than sigma, finally shift
      //of all fails: also norm (above)
      //as final resort: reduce the acceptance limit (see below)

      //default
      fSrc->FixParameter(1,SrcPar.sigma1);
      fSrc->FixParameter(2,SrcPar.shift1);
      fSrc->FixParameter(3,SrcPar.wght1);

      //sigma g2
      fSrc->SetParameter(4,rangen.Uniform(SrcPar.sigma1,SrcPar.sigma1*2));
      fSrc->SetParLimits(4,0,10.);
      //shift g2
      fSrc->SetParameter(5,0);
      fSrc->SetParLimits(5,0,3.);
      //weight g2
      fSrc->FixParameter(6,1);

      fSrc->FixParameter(7,1);
      fSrc->FixParameter(8,1);

      if(SrcPar.wght1){
        if(ResetCount>=1){
          double low_lim = SrcPar.wght1*0.9;
          double up_lim = SrcPar.wght1*1.1;
          if(up_lim>1) up_lim = 1;
          if(up_lim<0) up_lim = 0;
          if(low_lim>1) low_lim = 1;
          if(low_lim<0) low_lim = 0;
          fSrc->SetParameter(3,SrcPar.wght1);
          //printf("%f %f\n",low_lim,up_lim);
          fSrc->SetParLimits(3,low_lim,up_lim);
        }


        if(ResetCount>=2){
          fSrc->SetParameter(1,SrcPar.sigma1);
          fSrc->SetParLimits(1,SrcPar.sigma1*0.9,SrcPar.sigma1*1.1);
        }

        if(ResetCount>=3){
          fSrc->SetParameter(2,0);
          fSrc->SetParLimits(2,0,0.5);
        }

        if(ResetCount>=4){
          fSrc->SetParameter(6,0.9);
          fSrc->SetParLimits(6,0.8,1.0);

          fSrc->SetParameter(7,5);
          fSrc->SetParLimits(7,0,20.);
          //shift g2
          fSrc->SetParameter(8,0);
          fSrc->SetParLimits(8,0,10.);
        }

      }
      else{
        if(ResetCount>=1){
          fSrc->SetParameter(1,0.5);
          fSrc->SetParLimits(1,0,1.0);
          fSrc->SetParameter(3,0.1);
          fSrc->SetParLimits(3,0,0.2);
        }

        if(ResetCount>=3){
          fSrc->SetParameter(2,0);
          fSrc->SetParLimits(2,0,0.5);
        }

        if(ResetCount>=4){
          fSrc->SetParameter(6,0.9);
          fSrc->SetParLimits(6,0.8,1.0);

          fSrc->SetParameter(7,5);
          fSrc->SetParLimits(7,0,20.);
          //shift g2
          fSrc->SetParameter(8,0);
          fSrc->SetParLimits(8,0,10.);
        }

      }


    }
    else{
      printf("FitMode ERROR\n");
      return NULL;
    }


    hSrc->Fit(fSrc,"Q, S, N, R, M","",lowerlimit,upperlimit);

    //up to 10 fm
    const double RadDistMax = 8;
    for(unsigned uRad=0; uRad<hSrc->GetNbinsX(); uRad++){
      double Rad = hSrc->GetBinCenter(uRad+1);
      if(Rad>RadDistMax) break;
      double dst = hSrc->GetBinContent(uRad+1)-fSrc->Eval(Rad);
      double err;
      if(hSrc->GetBinContent(uRad+1)){
        err = hSrc->GetBinError(uRad+1);
      }
      else{
        err = fabs(dst)*1000;
      }

      Dist2 += dst*dst;
      NDPts++;
      if(hSrc->GetBinContent(uRad+1)){
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
    //if(Chi2/NDPts_chi2 < 4) nsig = sqrt(2)*TMath::ErfcInverse(TMath::Prob(Chi2,NDPts_chi2));
    //else nsig = 10;
    //printf("nsig = %f\n",nsig);
    Chi2 /= NDPts_chi2;
    //printf("Chi2 = %f\n",Chi2);
    StuckCount++;
    if(StuckCount>Patience){
      //if(ResetCount%2) Dist_CurLim += Dist_Limit;
      if(ResetCount>5) Dist_CurLim += Dist_Limit;
      StuckCount = 0;
      ResetCount++;
    }
  }
  while(Dist_Max>Dist_CurLim && Nsig_AtDistMax>NsigMax_Limit && Chi2>Chi2_Limit);
  if(Dist_Max>BadFitWarning && Nsig_AtDistMax>NsigMax_Limit && Chi2>Chi2_Limit){
    printf("WARNING: BadFitWarning (Chi2/ndf = %.2f)\n",Chi2);
  }

  if(FitMode==1){
    SrcPar.norm = fSrc->GetParameter(0);
    SrcPar.sigma1 = fSrc->GetParameter(1);
    SrcPar.shift1 = fSrc->GetParameter(2);
    if(SrcPar.shift1<0){
      printf("wtf1\n");
    }
    SrcPar.wght1 = fSrc->GetParameter(3);
  }
  else if(FitMode==2){
    SrcPar.norm = fSrc->GetParameter(0);
    SrcPar.sigma1 = fSrc->GetParameter(1);
    SrcPar.shift1 = fSrc->GetParameter(2);
    if(SrcPar.shift1<0){
      printf("wtf2\n");
    }
    SrcPar.wght1 = fSrc->GetParameter(3);
    SrcPar.sigma2 = fSrc->GetParameter(4);
    SrcPar.shift2 = fSrc->GetParameter(5);
    SrcPar.wght2 = fSrc->GetParameter(6);
    SrcPar.sigma3 = fSrc->GetParameter(7);
    SrcPar.shift3 = fSrc->GetParameter(8);
  }

  return fSrc;
}










//reads my batch output and fits it, saving the output + info of all pars from the master file
//works only for wildcard_flag==0
int dlmhst_ceca_fit_wc0(int argc, char *argv[]){
//printf("hi...\n");
  if(argc!=5){
    printf("\033[1;31mERROR:\033[0m dlmhst_ceca_fit_wc0 needs 5 args!\n");
    return 1;
  }

  //TString DADD_IN=TString(argv[1]);

  //TString DADD_OUT=TString(argv[1]);
  //TString ROOT_OUT=TString(argv[2]);
  //TString MASTER_FILE=TString(argv[3]);
  //TString FINAL_OUT=TString(argv[4]);

  TString BASE_NAME_PAR = TString(argv[1]);
  TString OUT_FOLDER = TString(argv[2]);//the dlm.hst files
  TString MASTER_FOLDER = TString(argv[3]);
  TString OUTPUT_FOLDER = TString(argv[4]);

  TString DADD_CORE_OUT = OUT_FOLDER+BASE_NAME_PAR+TString(".core.dlm.hst");
  TString DADD_FULL_OUT = OUT_FOLDER+BASE_NAME_PAR+TString(".full.dlm.hst");

  TString ROOT_CORE_OUT = OUT_FOLDER+BASE_NAME_PAR+TString(".core.root");
  TString ROOT_FULL_OUT = OUT_FOLDER+BASE_NAME_PAR+TString(".full.root");

  TString MASTER_FILE = MASTER_FOLDER+BASE_NAME_PAR+TString(".dlm.master");
  TString FINAL_OUT = OUTPUT_FOLDER+BASE_NAME_PAR+TString(".ceca.source");

//printf("dlmhst_root...\n");
  if(dlmhst_root(DADD_CORE_OUT,ROOT_CORE_OUT)){
    return 1;
  }
  if(dlmhst_root(DADD_FULL_OUT,ROOT_FULL_OUT)){
    return 1;
  }
//printf(" --> done\n");

  char* cline = new char [512];
  char* cdscr = new char [128];
  char* cval = new char [128];
  double read_value;
  double mom_bin_width = -1;
  int wildcard_flag;
  TString type = "";

  const unsigned MaxNumLines = 256;
  unsigned NumLines = 0;
  //we save all the info from the master file, that we want to save in our final output
  //this is everything apart from the JOB_YIELD stuff (batch output)
  char** LinesToSave = new char* [MaxNumLines];
  for(unsigned uLine=0; uLine<MaxNumLines; uLine++){
    LinesToSave[uLine] = new char [256];
  }
//printf("InFile...\n");
  FILE *InFile;
  InFile = fopen(MASTER_FILE.Data(), "r");
  if(!InFile){
      printf("\033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", MASTER_FILE.Data());
      return 1;
  }
  fseek ( InFile , 0 , SEEK_END );
  long EndPos;
  EndPos = ftell (InFile);
  fseek ( InFile , 0 , SEEK_SET );
  long CurPos;
  while(!feof(InFile)){
    if(!fgets(cline, 511, InFile)){
      //printf("\033[1;31mERROR:\033[0m The file\033[0m %s cannot be properly read (%s)!\n", InputFileName.Data(),cline);
    }
    sscanf(cline, "%s %s",cdscr,cval);
    if(strcmp(cdscr,"mom_bin_width")==0){
      mom_bin_width = atof(cval);
    }
    if(strcmp(cdscr,"wildcard_flag")==0){
      wildcard_flag = atoi(cval);
      if(wildcard_flag){
        printf("\033[1;31mERROR:\033[0m dlmhst_ceca_fit_wc0 works only for wildcard_flag==0!\n");
        return 1;
      }
    }
    if(strcmp(cdscr,"type")==0){
      if(strcmp(cval,"pp")&&strcmp(cval,"pL")){
        printf("\033[1;31mERROR:\033[0m dlmhst_ceca_fit_wc0 works only for pp and pL!\n");
        return 1;
      }
      type = TString(cval);
    }

    if(strcmp(cdscr,"JOB_ACTIVE")&&strcmp(cdscr,"JOB_YIELD")){
      strcpy(LinesToSave[NumLines],cline);
      NumLines++;
    }
  }//InFile
  fclose(InFile);

  //printf(" --> done\n");

  unsigned UpKstarBin;

  TFile root_file_core(ROOT_CORE_OUT,"update");
  TFile root_file_full(ROOT_FULL_OUT,"update");

  if(!root_file_core.IsOpen()){
    printf("\033[1;31mERROR:\033[0m Cannot open %s\n",ROOT_CORE_OUT.Data());
    return 1;
  }
  if(!root_file_full.IsOpen()){
    printf("\033[1;31mERROR:\033[0m Cannot open %s\n",ROOT_FULL_OUT.Data());
    return 1;
  }

  root_file_core.cd();
  TH2F* hYield_core = (TH2F*)root_file_core.Get("hYield");
  if(!hYield_core){
    printf("\033[1;31mERROR:\033[0m Cannot open hYield in %s\n",ROOT_CORE_OUT.Data());
    root_file_core.Close();
    return 1;
  }
  root_file_full.cd();
  TH2F* hYield_full = (TH2F*)root_file_full.Get("hYield");
  if(!hYield_full){
    printf("\033[1;31mERROR:\033[0m Cannot open hYield in %s\n",ROOT_FULL_OUT.Data());
    root_file_full.Close();
    return 1;
  }


  const unsigned NumMtBins = hYield_full->GetYaxis()->GetNbins();
  if(type=="pp"&&NumMtBins!=10){
    printf("\033[1;31mERROR:\033[0m We have to have 10 mT bins in pp\n");
    root_file_core.Close();
    root_file_full.Close();
    return 1;
  }
  if(type=="pL"&&NumMtBins!=8){
    printf("\033[1;31mERROR:\033[0m We have to have 8 mT bins in pL\n");
    root_file_core.Close();
    root_file_full.Close();
    return 1;
  }

  float* MtBinCenter = new float [NumMtBins];
  TriGauss* SrcPar = new TriGauss[NumMtBins];
  //the loopwhole ghetto: we dont check what is the mean mT
  //here we take the values based on jaime pars ran with high statistics
  //hardcoded for wildcard_flag==0
  if(type=="pp"){
    MtBinCenter[0] = 983;
    MtBinCenter[1] = 1054;
    MtBinCenter[2] = 1110;
    MtBinCenter[3] = 1168;
    MtBinCenter[4] = 1228;
    MtBinCenter[5] = 1315;
    MtBinCenter[6] = 1463;
    MtBinCenter[7] = 1681;
    MtBinCenter[8] = 1923;
    MtBinCenter[9] = 2303;
  }
  //pL
  else{
    MtBinCenter[0] = 1121;
    MtBinCenter[1] = 1210;
    MtBinCenter[2] = 1288;
    MtBinCenter[3] = 1377;
    MtBinCenter[4] = 1536;
    MtBinCenter[5] = 1753;
    MtBinCenter[6] = 1935;
    MtBinCenter[7] = 2334;
  }

  double* Dist_Max = new double [NumMtBins];
  double* Chi2Ndf = new double [NumMtBins];

  for(unsigned uMt=0; uMt<NumMtBins; uMt++){
    root_file_core.cd();
    TH2F* h_kstar_rstar_core = (TH2F*)root_file_core.Get(TString::Format("h_kstar_rstar_mT%u",uMt));
    root_file_full.cd();
    TH2F* h_kstar_rstar_full = (TH2F*)root_file_full.Get(TString::Format("h_kstar_rstar_mT%u",uMt));
    if(!h_kstar_rstar_core || !h_kstar_rstar_full){
      printf("\033[1;31mERROR:\033[0m Cannot open h_kstar_rstar_mT%u in %s or %s\n",uMt,ROOT_CORE_OUT.Data(),ROOT_FULL_OUT.Data());
      root_file_core.Close();
      root_file_full.Close();
      delete [] MtBinCenter;
      return 1;
    }

    gROOT->cd();

    if(mom_bin_width==-1)
      mom_bin_width=100;

    UpKstarBin = h_kstar_rstar_full->GetXaxis()->FindBin(mom_bin_width);

    TH1F* hSrc_core = (TH1F*)h_kstar_rstar_core->ProjectionY(TString::Format("hSrc_core_%u",uMt),1,UpKstarBin);
    double Int_Src_core = hSrc_core->Integral();
    if(Int_Src_core) hSrc_core->Scale(1./Int_Src_core,"width");

    TH1F* hSrc_full = (TH1F*)h_kstar_rstar_full->ProjectionY(TString::Format("hSrc_full_%u",uMt),1,UpKstarBin);
    double Int_Src_full = hSrc_full->Integral();
    if(Int_Src_full) hSrc_full->Scale(1./Int_Src_full,"width");

    double WeightCore = 0;
    if(Int_Src_core+Int_Src_full) WeightCore = Int_Src_core/(Int_Src_core+Int_Src_full);


    //printf("so far so core %u\n",uMt);
    //CORE FIT

    TF1* fSrc_core = NULL;
    //, Dist_Max[uMt], Nsig[uMt] will be overwritten later, but its fine, we only want to QA final result
    if(Int_Src_core) fSrc_core = Fit_wc0(1, hSrc_core,SrcPar[uMt], Dist_Max[uMt], Chi2Ndf[uMt]);
    else{
      SrcPar[uMt].sigma1 = 0.5;
      SrcPar[uMt].shift1 = 0;
      SrcPar[uMt].wght1 = 0;
    }
    if(fSrc_core) fSrc_core->SetName(TString::Format("fSrc_core_%u",uMt));
    //printf("SrcPar[uMt].shift1 = %f\n",SrcPar[uMt].shift1);

    //printf("so far so full %u\n",uMt);
    //FULL FIT
    TF1* fSrc_full = NULL;
    SrcPar[uMt].wght1 = WeightCore;
    if(Int_Src_core+Int_Src_full) fSrc_full = Fit_wc0(2,hSrc_full,SrcPar[uMt], Dist_Max[uMt], Chi2Ndf[uMt]);
    if(fSrc_full) fSrc_full->SetName(TString::Format("fSrc_full_%u",uMt));

    //printf("..write c\n");
    root_file_core.cd();
    if(hSrc_core) hSrc_core->Write();
    if(fSrc_core) fSrc_core->Write();
    //printf("..del c\n");
    if(fSrc_core) delete fSrc_core;
    if(hSrc_core) delete hSrc_core;

    //printf("..write f\n");
    root_file_full.cd();
    if(hSrc_full) hSrc_full->Write();
    if(fSrc_full) fSrc_full->Write();
    //printf("..del f\n");
    if(hSrc_full) delete hSrc_full;
    if(fSrc_full) delete fSrc_full;
  }//uMt

//TF1* fSrc2 = new TF1(TString::Format("fSrc2"),NormTripleShiftedGaussTF1,0.,20.,9);
//fSrc2->FixParameter(0,(1.0000+0.9838)*0.5);
//fSrc2->FixParameter(1,(0.897+0.688)*0.5);
//fSrc2->FixParameter(2,(0.0000+0.029)*0.5);
//fSrc2->FixParameter(3,(0.1500+0.0538)*0.5);
//fSrc2->FixParameter(4,(1.841+1.396)*0.5);
//fSrc2->FixParameter(5,(0.0000+0.0)*0.5);
//fSrc2->FixParameter(6,(1.0000+1.0)*0.5);
//fSrc2->FixParameter(7,(6.500+6.765)*0.5);
//fSrc2->FixParameter(8,(0.500+0.765)*0.5);
//root_file_full.cd();
//fSrc2->Write();

  FILE *fptr;
  fptr=fopen(FINAL_OUT.Data(),"w");
  //NEW LINE NEEDED???
  for(unsigned uLine=0; uLine<NumLines; uLine++){
    //printf("%s\n",LinesToSave[uLine]);
    fprintf(fptr,"%s",LinesToSave[uLine]);
  }
  fprintf(fptr,"SOURCE:\n");
  for(unsigned uMt=0; uMt<NumMtBins; uMt++){
    //norm,sig1,shift1,w1,sig2,shift2,w2,sig3,shift3
    fprintf(fptr,"%.0f %.4f %.3f %.3f %.4f %.3f %.3f %.4f %.3f %.3f %.3f %.1f\n",
        MtBinCenter[uMt],
        SrcPar[uMt].norm,
        SrcPar[uMt].sigma1,
        SrcPar[uMt].shift1,
        SrcPar[uMt].wght1,
        SrcPar[uMt].sigma2,
        SrcPar[uMt].shift2,
        SrcPar[uMt].wght2,
        SrcPar[uMt].sigma3,
        SrcPar[uMt].shift3,
        Dist_Max[uMt],
        Chi2Ndf[uMt]
      );
  }


  for(unsigned uLine=0; uLine<MaxNumLines; uLine++){
    delete [] LinesToSave[uLine];
  }
  delete [] LinesToSave;

  delete [] cline;
  delete [] cdscr;
  delete [] cval;

  root_file_core.Close();
  root_file_full.Close();
  delete [] MtBinCenter;
  delete [] SrcPar;
  delete [] Dist_Max;
  delete [] Chi2Ndf;
  return 0;
}








//FitMode:
// 1 - core (single non-shifted, fully normalized, Gauss). Saves the pars within the 1st Gauss in SrcPar
// 2 - full: takes the pars of the first Gauss in SrcPar to fix the first Gaussian,
//            and fits by adding a second Gaussian that can be shifted
//            G3 is OFF!! if we are stuck, we change a bit the core gauss
TF1* Fit_spl_wc0(int FitMode, TH1F* hSrc, SplPars& SrcPar, double& Chi2, double& Correction){
  double lowerlimit, upperlimit;
  GetCentralInterval(*hSrc, 0.98, lowerlimit, upperlimit, true);
//printf("SrcPar.shift1 = %f\n",SrcPar.shift1);
//usleep(100e3);
  if(lowerlimit>5) lowerlimit = 5;
  if(upperlimit>20) upperlimit = 20;


  //double DlmTSplineFit(double* xVal, double* pars){
      //[0] = NumKnots
      //[1] = der at 0
      //[2] = der at last
      //[3]... posX
      //[...]... poxY
  const unsigned NumKnots = 10;
  const unsigned NumFitPars = 3+NumKnots*2;
  const double SplFitMin = 0;
  const double SplFitMax = 8;//HARDCODED FOR THIS CLASS !!!
  const double Chi2_Limit = 2;

  TF1* fSrc;
  SetUpSplPars(fSrc);
  for(unsigned uKnot=0; uKnot<NumKnots; uKnot++){
    double x_val = fSrc->GetParameter(3+uKnot);
    double y_val = hSrc->GetBinContent(hSrc->FindBin(x_val));
    double y_err = hSrc->GetBinError(hSrc->FindBin(x_val));
    if(y_val==0 || x_val==0){
      fSrc->FixParameter(3+NumKnots+uKnot,0);
    }
    //get as close to zero as possible
    else if(x_val==SplFitMax){
      fSrc->SetParameter(3+NumKnots+uKnot,y_val*0.25);
      fSrc->SetParLimits(3+NumKnots+uKnot,0,y_val*0.5);
      if(y_val==0){
        fSrc->SetParLimits(3+NumKnots+uKnot,0,0.01);
      }
    }
    else{
      fSrc->SetParameter(3+NumKnots+uKnot,y_val);
      //to reduce prob of getting large integral (0.0*y_err)
      fSrc->SetParLimits(3+NumKnots+uKnot,0,y_val+0.0*y_err);
    }
  }


  double Integral = fSrc->Integral(SplFitMin,SplFitMax);
  Correction = 1./Integral;
  //we only worry if the integral is > 1, otherwise the "missing" yield will be obsorbed by cats
  if(Correction>1) Correction=1;

/*
  unsigned Counter=0;
  double SmallBias;
  while(fabs(Correction-1)>0.01){
    SmallBias = 0.1*Correction*double(Counter)+0.9;
    for(unsigned uKnot=0; uKnot<NumKnots; uKnot++){
      fSrc->SetParameter(3+NumKnots+uKnot,fSrc->GetParameter(3+NumKnots+uKnot)*Correction);
      fSrc->SetParLimits(3+NumKnots+uKnot,0,fSrc->GetParameter(3+NumKnots+uKnot)*SmallBias);
    }
    hSrc->Fit(fSrc,"Q, S, N, R, M","",lowerlimit,upperlimit);
    Integral = fSrc->Integral(SplFitMin,SplFitMax);
    Correction = 1./Integral;
    Counter++;
  }
*/

  for(unsigned uKnot=0; uKnot<NumKnots; uKnot++){
    fSrc->FixParameter(3+NumKnots+uKnot,fSrc->GetParameter(3+NumKnots+uKnot)*Correction);
    SrcPar.KnotY[uKnot] = fSrc->GetParameter(3+NumKnots+uKnot);
  }
  double IntegralCorr = fSrc->Integral(SplFitMin,SplFitMax);
  if(IntegralCorr>1+1e-6){
    printf("SUPER BIG BUG WITH THE SPLINES!!!\n");
  }

  Chi2 = 0;
  double NDPts_chi2 = 0;

  //up to 8 fm
  const double RadDistMax = SplFitMax;
  for(unsigned uRad=0; uRad<hSrc->GetNbinsX(); uRad++){
    double Rad = hSrc->GetBinCenter(uRad+1);
    if(Rad>RadDistMax) break;
    double dst = hSrc->GetBinContent(uRad+1)-fSrc->Eval(Rad);
    double err;
    if(hSrc->GetBinContent(uRad+1)){
      err = hSrc->GetBinError(uRad+1);
    }
    else{
      err = fabs(dst)*1000;
    }

    if(hSrc->GetBinContent(uRad+1)){
      Chi2 += (dst*dst)/(err*err);
      NDPts_chi2++;
    }

  }
  Chi2 /= NDPts_chi2;

  if(Chi2>Chi2_Limit){
    printf("WARNING: BadFitWarning (Chi2/ndf = %.2f)\n",Chi2);
  }

  return fSrc;
}









//reads my batch output and fits it, saving the output + info of all pars from the master file
//works only for wildcard_flag==0
int dlmhst_ceca_splfit_wc0(int argc, char *argv[]){
//printf("hi...\n");
  if(argc!=5){
    printf("\033[1;31mERROR:\033[0m dlmhst_ceca_fit_wc0 needs 5 args!\n");
    return 1;
  }

  TString BASE_NAME_PAR = TString(argv[1]);
  TString OUT_FOLDER = TString(argv[2]);//the dlm.hst files
  TString MASTER_FOLDER = TString(argv[3]);
  TString OUTPUT_FOLDER = TString(argv[4]);

  TString DADD_CORE_OUT = OUT_FOLDER+BASE_NAME_PAR+TString(".core.dlm.hst");
  TString DADD_FULL_OUT = OUT_FOLDER+BASE_NAME_PAR+TString(".full.dlm.hst");

  TString ROOT_CORE_OUT = OUT_FOLDER+BASE_NAME_PAR+TString(".core.root");
  TString ROOT_FULL_OUT = OUT_FOLDER+BASE_NAME_PAR+TString(".full.root");

  TString MASTER_FILE = MASTER_FOLDER+BASE_NAME_PAR+TString(".dlm.master");
  TString FINAL_OUT = OUTPUT_FOLDER+BASE_NAME_PAR+TString(".ceca.source");

  if(dlmhst_root(DADD_CORE_OUT,ROOT_CORE_OUT)){
    return 1;
  }
  if(dlmhst_root(DADD_FULL_OUT,ROOT_FULL_OUT)){
    return 1;
  }

  char* cline = new char [512];
  char* cdscr = new char [128];
  char* cval = new char [128];
  double read_value;
  double mom_bin_width = -1;
  int wildcard_flag;
  TString type = "";

  const unsigned MaxNumLines = 256;
  unsigned NumLines = 0;
  //we save all the info from the master file, that we want to save in our final output
  //this is everything apart from the JOB_YIELD stuff (batch output)
  char** LinesToSave = new char* [MaxNumLines];
  for(unsigned uLine=0; uLine<MaxNumLines; uLine++){
    LinesToSave[uLine] = new char [256];
  }
  FILE *InFile;
  InFile = fopen(MASTER_FILE.Data(), "r");
  if(!InFile){
      printf("\033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", MASTER_FILE.Data());
      return 1;
  }
  fseek ( InFile , 0 , SEEK_END );
  long EndPos;
  EndPos = ftell (InFile);
  fseek ( InFile , 0 , SEEK_SET );
  long CurPos;
  while(!feof(InFile)){
    if(!fgets(cline, 511, InFile)){
      //printf("\033[1;31mERROR:\033[0m The file\033[0m %s cannot be properly read (%s)!\n", InputFileName.Data(),cline);
    }
    sscanf(cline, "%s %s",cdscr,cval);
    if(strcmp(cdscr,"mom_bin_width")==0){
      mom_bin_width = atof(cval);
    }
    if(strcmp(cdscr,"wildcard_flag")==0){
      wildcard_flag = atoi(cval);
      if(wildcard_flag){
        printf("\033[1;31mERROR:\033[0m dlmhst_ceca_fit_wc0 works only for wildcard_flag==0!\n");
        return 1;
      }
    }
    if(strcmp(cdscr,"type")==0){
      if(strcmp(cval,"pp")&&strcmp(cval,"pL")){
        printf("\033[1;31mERROR:\033[0m dlmhst_ceca_fit_wc0 works only for pp and pL!\n");
        return 1;
      }
      type = TString(cval);
    }

    if(strcmp(cdscr,"JOB_ACTIVE")&&strcmp(cdscr,"JOB_YIELD")){
      strcpy(LinesToSave[NumLines],cline);
      NumLines++;
    }
  }//InFile
  fclose(InFile);

  unsigned UpKstarBin;

  TFile root_file_core(ROOT_CORE_OUT,"update");
  TFile root_file_full(ROOT_FULL_OUT,"update");

  if(!root_file_core.IsOpen()){
    printf("\033[1;31mERROR:\033[0m Cannot open %s\n",ROOT_CORE_OUT.Data());
    return 1;
  }
  if(!root_file_full.IsOpen()){
    printf("\033[1;31mERROR:\033[0m Cannot open %s\n",ROOT_FULL_OUT.Data());
    return 1;
  }

  root_file_core.cd();
  TH2F* hYield_core = (TH2F*)root_file_core.Get("hYield");
  if(!hYield_core){
    printf("\033[1;31mERROR:\033[0m Cannot open hYield in %s\n",ROOT_CORE_OUT.Data());
    root_file_core.Close();
    return 1;
  }
  root_file_full.cd();
  TH2F* hYield_full = (TH2F*)root_file_full.Get("hYield");
  if(!hYield_full){
    printf("\033[1;31mERROR:\033[0m Cannot open hYield in %s\n",ROOT_FULL_OUT.Data());
    root_file_full.Close();
    return 1;
  }


  const unsigned NumMtBins = hYield_full->GetYaxis()->GetNbins();
  if(type=="pp"&&NumMtBins!=10){
    printf("\033[1;31mERROR:\033[0m We have to have 10 mT bins in pp\n");
    root_file_core.Close();
    root_file_full.Close();
    return 1;
  }
  if(type=="pL"&&NumMtBins!=8){
    printf("\033[1;31mERROR:\033[0m We have to have 8 mT bins in pL\n");
    root_file_core.Close();
    root_file_full.Close();
    return 1;
  }

  float* MtBinCenter = new float [NumMtBins];



  SplPars* SrcPar = new SplPars[NumMtBins];


  //the loopwhole ghetto: we dont check what is the mean mT
  //here we take the values based on jaime pars ran with high statistics
  //hardcoded for wildcard_flag==0
  if(type=="pp"){
    MtBinCenter[0] = 983;
    MtBinCenter[1] = 1054;
    MtBinCenter[2] = 1110;
    MtBinCenter[3] = 1168;
    MtBinCenter[4] = 1228;
    MtBinCenter[5] = 1315;
    MtBinCenter[6] = 1463;
    MtBinCenter[7] = 1681;
    MtBinCenter[8] = 1923;
    MtBinCenter[9] = 2303;
  }
  //pL
  else{
    MtBinCenter[0] = 1121;
    MtBinCenter[1] = 1210;
    MtBinCenter[2] = 1288;
    MtBinCenter[3] = 1377;
    MtBinCenter[4] = 1536;
    MtBinCenter[5] = 1753;
    MtBinCenter[6] = 1935;
    MtBinCenter[7] = 2334;
  }

  double* Chi2Ndf = new double [NumMtBins];
  double* Correction = new double [NumMtBins];

  for(unsigned uMt=0; uMt<NumMtBins; uMt++){
    root_file_core.cd();
    TH2F* h_kstar_rstar_core = (TH2F*)root_file_core.Get(TString::Format("h_kstar_rstar_mT%u",uMt));
    root_file_full.cd();
    TH2F* h_kstar_rstar_full = (TH2F*)root_file_full.Get(TString::Format("h_kstar_rstar_mT%u",uMt));
    if(!h_kstar_rstar_core || !h_kstar_rstar_full){
      printf("\033[1;31mERROR:\033[0m Cannot open h_kstar_rstar_mT%u in %s or %s\n",uMt,ROOT_CORE_OUT.Data(),ROOT_FULL_OUT.Data());
      root_file_core.Close();
      root_file_full.Close();
      delete [] MtBinCenter;
      return 1;
    }

    gROOT->cd();

    if(mom_bin_width==-1)
      mom_bin_width=100;

    UpKstarBin = h_kstar_rstar_full->GetXaxis()->FindBin(mom_bin_width);

    TH1F* hSrc_core = (TH1F*)h_kstar_rstar_core->ProjectionY(TString::Format("hSrc_core_%u",uMt),1,UpKstarBin);
    double Int_Src_core = hSrc_core->Integral();
    if(Int_Src_core) hSrc_core->Scale(1./Int_Src_core,"width");

    TH1F* hSrc_full = (TH1F*)h_kstar_rstar_full->ProjectionY(TString::Format("hSrc_full_%u",uMt),1,UpKstarBin);
    double Int_Src_full = hSrc_full->Integral();
    if(Int_Src_full) hSrc_full->Scale(1./Int_Src_full,"width");

    double WeightCore = 0;
    if(Int_Src_core+Int_Src_full) WeightCore = Int_Src_core/(Int_Src_core+Int_Src_full);

    //CORE FIT
    TF1* fSrc_core = NULL;

    if(Int_Src_core) fSrc_core = Fit_spl_wc0(1, hSrc_core,SrcPar[uMt],Chi2Ndf[uMt],Correction[uMt]);
    else{
      SrcPar[uMt] = 0;
    }
    if(fSrc_core) fSrc_core->SetName(TString::Format("fSrc_core_%u",uMt));

    //FULL FIT
    TF1* fSrc_full = NULL;
    if(Int_Src_core+Int_Src_full) fSrc_full = Fit_spl_wc0(2,hSrc_full,SrcPar[uMt], Chi2Ndf[uMt],Correction[uMt]);
    if(fSrc_full) fSrc_full->SetName(TString::Format("fSrc_full_%u",uMt));

    root_file_core.cd();
    if(hSrc_core) hSrc_core->Write();
    if(fSrc_core) fSrc_core->Write();
    if(fSrc_core) delete fSrc_core;
    if(hSrc_core) delete hSrc_core;

    root_file_full.cd();
    if(hSrc_full) hSrc_full->Write();
    if(fSrc_full) fSrc_full->Write();
    if(hSrc_full) delete hSrc_full;
    if(fSrc_full) delete fSrc_full;
  }//uMt



  FILE *fptr;
  fptr=fopen(FINAL_OUT.Data(),"w");
  //NEW LINE NEEDED???
  for(unsigned uLine=0; uLine<NumLines; uLine++){
    //printf("%s\n",LinesToSave[uLine]);
    fprintf(fptr,"%s",LinesToSave[uLine]);
  }
  fprintf(fptr,"SOURCE:\n");
  for(unsigned uMt=0; uMt<NumMtBins; uMt++){
    //norm,sig1,shift1,w1,sig2,shift2,w2,sig3,shift3
    fprintf(fptr,"%.0f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.2f %.4f\n",
        MtBinCenter[uMt],
        SrcPar[uMt].KnotY[0],
        SrcPar[uMt].KnotY[1],
        SrcPar[uMt].KnotY[2],
        SrcPar[uMt].KnotY[3],
        SrcPar[uMt].KnotY[4],
        SrcPar[uMt].KnotY[5],
        SrcPar[uMt].KnotY[6],
        SrcPar[uMt].KnotY[7],
        SrcPar[uMt].KnotY[8],
        SrcPar[uMt].KnotY[9],
        Chi2Ndf[uMt],
        Correction[uMt]
      );
  }


  for(unsigned uLine=0; uLine<MaxNumLines; uLine++){
    delete [] LinesToSave[uLine];
  }
  delete [] LinesToSave;

  delete [] cline;
  delete [] cdscr;
  delete [] cval;

  root_file_core.Close();
  root_file_full.Close();
  delete [] MtBinCenter;
  delete [] SrcPar;
  delete [] Chi2Ndf;
  delete [] Correction;
  return 0;
}












//FitMode:
// 1 - core
// 2 - full
TF1* Fit_kdp_wc0(int FitMode, TH1F* hSrc, KdpPars& SrcPar, double& Chi2){
  double lowerlimit, upperlimit;
  GetCentralInterval(*hSrc, 0.98, lowerlimit, upperlimit, true);
  //if(lowerlimit>5) lowerlimit = 5;
  lowerlimit = 0;
  if(upperlimit>10) upperlimit = 10;


  const double Chi2_Limit = 3;
  const double KdpFitMax = 8;

  TF1* fSrc;
  SetUpKdpPars(fSrc,2);

  //for(unsigned uP=0; uP<KdpPars::NumDistos; uP++){
  //  fSrc->SetParLimits(1+uP*3,0.1,fSrc->GetParameter(0+uP*3));
  //}
  //fSrc->FixParameter(0,1);
  //fSrc->FixParameter(1,0.5);
  //fSrc->FixParameter(2,0.5);

  //fSrc->FixParameter(3,6);
  //fSrc->FixParameter(4,2);
  //fSrc->FixParameter(5,1);

  hSrc->Fit(fSrc,"Q, S, N, R, M","",lowerlimit,upperlimit);

  double Integral = fSrc->Integral(0,KdpFitMax*4);
  if(fabs(Integral-1)>1e-2){
    printf("SUPER BIG BUG WITH THE KDP (%f)!!!\n",Integral);
  }

  Chi2 = 0;
  double NDPts_chi2 = 0;

  //up to 8 fm
  for(unsigned uRad=0; uRad<hSrc->GetNbinsX(); uRad++){
    double Rad = hSrc->GetBinCenter(uRad+1);
    if(Rad>KdpFitMax) break;
    double dst = hSrc->GetBinContent(uRad+1)-fSrc->Eval(Rad);
    double err = hSrc->GetBinError(uRad+1);
    if(hSrc->GetBinContent(uRad+1)){
      Chi2 += (dst*dst)/(err*err);
      NDPts_chi2++;
    }
  }
  Chi2 /= NDPts_chi2;

  if(Chi2>Chi2_Limit){
    printf("WARNING: BadFitWarning (Chi2/ndf = %.2f)\n",Chi2);
  }

  for(unsigned uP=0; uP<KdpPars::NumDistos; uP++){
    SrcPar.mean[uP] = fSrc->GetParameter(0+uP*3);
    SrcPar.stdv[uP] = fSrc->GetParameter(1+uP*3);
    if(uP!=KdpPars::NumDistos-1)
      SrcPar.wght[uP] = fSrc->GetParameter(2+uP*3);
  }

  return fSrc;
}





//reads my batch output and fits it, saving the output + info of all pars from the master file
//works only for wildcard_flag==0
//kdp = Kernel Density Poisson
int dlmhst_ceca_kdpfit_wc0(int argc, char *argv[]){
//printf("hi...\n");
  if(argc!=5){
    printf("\033[1;31mERROR:\033[0m dlmhst_ceca_fit_wc0 needs 5 args!\n");
    return 1;
  }

  TString BASE_NAME_PAR = TString(argv[1]);
  TString OUT_FOLDER = TString(argv[2]);//the dlm.hst files
  TString MASTER_FOLDER = TString(argv[3]);
  TString OUTPUT_FOLDER = TString(argv[4]);

  TString DADD_CORE_OUT = OUT_FOLDER+BASE_NAME_PAR+TString(".core.dlm.hst");
  TString DADD_FULL_OUT = OUT_FOLDER+BASE_NAME_PAR+TString(".full.dlm.hst");

  TString ROOT_CORE_OUT = OUT_FOLDER+BASE_NAME_PAR+TString(".core.root");
  TString ROOT_FULL_OUT = OUT_FOLDER+BASE_NAME_PAR+TString(".full.root");

  TString MASTER_FILE = MASTER_FOLDER+BASE_NAME_PAR+TString(".dlm.master");
  TString FINAL_OUT = OUTPUT_FOLDER+BASE_NAME_PAR+TString(".ceca.source");

  if(dlmhst_root(DADD_CORE_OUT,ROOT_CORE_OUT)){
    return 1;
  }
  if(dlmhst_root(DADD_FULL_OUT,ROOT_FULL_OUT)){
    return 1;
  }

  char* cline = new char [512];
  char* cdscr = new char [128];
  char* cval = new char [128];
  double read_value;
  double mom_bin_width = -1;
  int wildcard_flag;
  TString type = "";

  const unsigned MaxNumLines = 256;
  unsigned NumLines = 0;
  //we save all the info from the master file, that we want to save in our final output
  //this is everything apart from the JOB_YIELD stuff (batch output)
  char** LinesToSave = new char* [MaxNumLines];
  for(unsigned uLine=0; uLine<MaxNumLines; uLine++){
    LinesToSave[uLine] = new char [256];
  }
  FILE *InFile;
  InFile = fopen(MASTER_FILE.Data(), "r");
  if(!InFile){
      printf("\033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", MASTER_FILE.Data());
      return 1;
  }
  fseek ( InFile , 0 , SEEK_END );
  long EndPos;
  EndPos = ftell (InFile);
  fseek ( InFile , 0 , SEEK_SET );
  long CurPos;
  while(!feof(InFile)){
    if(!fgets(cline, 511, InFile)){
      //printf("\033[1;31mERROR:\033[0m The file\033[0m %s cannot be properly read (%s)!\n", InputFileName.Data(),cline);
    }
    sscanf(cline, "%s %s",cdscr,cval);
    if(strcmp(cdscr,"mom_bin_width")==0){
      mom_bin_width = atof(cval);
    }
    if(strcmp(cdscr,"wildcard_flag")==0){
      wildcard_flag = atoi(cval);
      if(wildcard_flag){
        printf("\033[1;31mERROR:\033[0m dlmhst_ceca_fit_wc0 works only for wildcard_flag==0!\n");
        return 1;
      }
    }
    if(strcmp(cdscr,"type")==0){
      if(strcmp(cval,"pp")&&strcmp(cval,"pL")){
        printf("\033[1;31mERROR:\033[0m dlmhst_ceca_fit_wc0 works only for pp and pL!\n");
        return 1;
      }
      type = TString(cval);
    }

    if(strcmp(cdscr,"JOB_ACTIVE")&&strcmp(cdscr,"JOB_YIELD")){
      strcpy(LinesToSave[NumLines],cline);
      NumLines++;
    }
  }//InFile
  fclose(InFile);

  unsigned UpKstarBin;

  TFile root_file_core(ROOT_CORE_OUT,"update");
  TFile root_file_full(ROOT_FULL_OUT,"update");

  if(!root_file_core.IsOpen()){
    printf("\033[1;31mERROR:\033[0m Cannot open %s\n",ROOT_CORE_OUT.Data());
    return 1;
  }
  if(!root_file_full.IsOpen()){
    printf("\033[1;31mERROR:\033[0m Cannot open %s\n",ROOT_FULL_OUT.Data());
    return 1;
  }

  root_file_core.cd();
  TH2F* hYield_core = (TH2F*)root_file_core.Get("hYield");
  if(!hYield_core){
    printf("\033[1;31mERROR:\033[0m Cannot open hYield in %s\n",ROOT_CORE_OUT.Data());
    root_file_core.Close();
    return 1;
  }
  root_file_full.cd();
  TH2F* hYield_full = (TH2F*)root_file_full.Get("hYield");
  if(!hYield_full){
    printf("\033[1;31mERROR:\033[0m Cannot open hYield in %s\n",ROOT_FULL_OUT.Data());
    root_file_full.Close();
    return 1;
  }


  const unsigned NumMtBins = hYield_full->GetYaxis()->GetNbins();
  if(type=="pp"&&NumMtBins!=10){
    printf("\033[1;31mERROR:\033[0m We have to have 10 mT bins in pp\n");
    root_file_core.Close();
    root_file_full.Close();
    return 1;
  }
  if(type=="pL"&&NumMtBins!=8){
    printf("\033[1;31mERROR:\033[0m We have to have 8 mT bins in pL\n");
    root_file_core.Close();
    root_file_full.Close();
    return 1;
  }

  float* MtBinCenter = new float [NumMtBins];

  KdpPars* SrcPar = new KdpPars[NumMtBins];

  //the loopwhole ghetto: we dont check what is the mean mT
  //here we take the values based on jaime pars ran with high statistics
  //hardcoded for wildcard_flag==0
  if(type=="pp"){
    MtBinCenter[0] = 983;
    MtBinCenter[1] = 1054;
    MtBinCenter[2] = 1110;
    MtBinCenter[3] = 1168;
    MtBinCenter[4] = 1228;
    MtBinCenter[5] = 1315;
    MtBinCenter[6] = 1463;
    MtBinCenter[7] = 1681;
    MtBinCenter[8] = 1923;
    MtBinCenter[9] = 2303;
  }
  //pL
  else{
    MtBinCenter[0] = 1121;
    MtBinCenter[1] = 1210;
    MtBinCenter[2] = 1288;
    MtBinCenter[3] = 1377;
    MtBinCenter[4] = 1536;
    MtBinCenter[5] = 1753;
    MtBinCenter[6] = 1935;
    MtBinCenter[7] = 2334;
  }

  double* Chi2Ndf = new double [NumMtBins];

  for(unsigned uMt=0; uMt<NumMtBins; uMt++){
    root_file_core.cd();
    TH2F* h_kstar_rstar_core = (TH2F*)root_file_core.Get(TString::Format("h_kstar_rstar_mT%u",uMt));
    root_file_full.cd();
    TH2F* h_kstar_rstar_full = (TH2F*)root_file_full.Get(TString::Format("h_kstar_rstar_mT%u",uMt));
    if(!h_kstar_rstar_core || !h_kstar_rstar_full){
      printf("\033[1;31mERROR:\033[0m Cannot open h_kstar_rstar_mT%u in %s or %s\n",uMt,ROOT_CORE_OUT.Data(),ROOT_FULL_OUT.Data());
      root_file_core.Close();
      root_file_full.Close();
      delete [] MtBinCenter;
      return 1;
    }

    gROOT->cd();

    if(mom_bin_width==-1)
      mom_bin_width=100;

    UpKstarBin = h_kstar_rstar_full->GetXaxis()->FindBin(mom_bin_width);

    double ErrAtZero;

    TH1F* hSrc_core = (TH1F*)h_kstar_rstar_core->ProjectionY(TString::Format("hSrc_core_%u",uMt),1,UpKstarBin);
    double Int_Src_core = hSrc_core->Integral();
    ErrAtZero = 1;
    ErrAtZero /= Int_Src_core;
    ErrAtZero /= hSrc_core->GetBinWidth(1);
    if(Int_Src_core) hSrc_core->Scale(1./Int_Src_core,"width");
    for(unsigned uBin=0; uBin<hSrc_core->GetNbinsX(); uBin++){
      if(hSrc_core->GetBinContent(uBin+1)==0){
        hSrc_core->SetBinError(uBin+1,ErrAtZero);
      }
    }


    TH1F* hSrc_full = (TH1F*)h_kstar_rstar_full->ProjectionY(TString::Format("hSrc_full_%u",uMt),1,UpKstarBin);
    double Int_Src_full = hSrc_full->Integral();
    ErrAtZero = 1;
    ErrAtZero /= Int_Src_full;
    ErrAtZero /= hSrc_full->GetBinWidth(1);
    if(Int_Src_full) hSrc_full->Scale(1./Int_Src_full,"width");
    for(unsigned uBin=0; uBin<hSrc_full->GetNbinsX(); uBin++){
      if(hSrc_full->GetBinContent(uBin+1)==0){
        hSrc_full->SetBinError(uBin+1,ErrAtZero);
      }
    }


    //double WeightCore = 0;
    //if(Int_Src_core+Int_Src_full) WeightCore = Int_Src_core/(Int_Src_core+Int_Src_full);

    //CORE FIT
    TF1* fSrc_core = NULL;

    if(Int_Src_core) fSrc_core = Fit_kdp_wc0(1, hSrc_core,SrcPar[uMt],Chi2Ndf[uMt]);
    else{
      SrcPar[uMt] = 0;
    }
    if(fSrc_core) fSrc_core->SetName(TString::Format("fSrc_core_%u",uMt));

    //FULL FIT
    TF1* fSrc_full = NULL;
    if(Int_Src_core+Int_Src_full) fSrc_full = Fit_kdp_wc0(2,hSrc_full,SrcPar[uMt], Chi2Ndf[uMt]);
    if(fSrc_full) fSrc_full->SetName(TString::Format("fSrc_full_%u",uMt));

    root_file_core.cd();
    if(hSrc_core) hSrc_core->Write();
    if(fSrc_core) fSrc_core->Write();
    if(fSrc_core) delete fSrc_core;
    if(hSrc_core) delete hSrc_core;

    root_file_full.cd();
    if(hSrc_full) hSrc_full->Write();
    if(fSrc_full) fSrc_full->Write();
    if(hSrc_full) delete hSrc_full;
    if(fSrc_full) delete fSrc_full;
  }//uMt



  FILE *fptr;
  fptr=fopen(FINAL_OUT.Data(),"w");
  //NEW LINE NEEDED???
  for(unsigned uLine=0; uLine<NumLines; uLine++){
    //printf("%s\n",LinesToSave[uLine]);
    fprintf(fptr,"%s",LinesToSave[uLine]);
  }
  fprintf(fptr,"SOURCE:\n");
  for(unsigned uMt=0; uMt<NumMtBins; uMt++){
    //the first number is the FitMode of the KDP, just that we do not forget!
    fprintf(fptr,"%i %.0f %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e "
                      "%.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e "
                      "%.3e %.3e %.3e %.3e %.3e %.3e %.3e %.2f\n",
        2, MtBinCenter[uMt],
        SrcPar[uMt].mean[0],SrcPar[uMt].stdv[0],SrcPar[uMt].wght[0],
        SrcPar[uMt].mean[1],SrcPar[uMt].stdv[1],SrcPar[uMt].wght[1],
        SrcPar[uMt].mean[2],SrcPar[uMt].stdv[2],SrcPar[uMt].wght[2],
        SrcPar[uMt].mean[3],SrcPar[uMt].stdv[3],SrcPar[uMt].wght[3],
        SrcPar[uMt].mean[4],SrcPar[uMt].stdv[4],SrcPar[uMt].wght[4],
        SrcPar[uMt].mean[5],SrcPar[uMt].stdv[5],SrcPar[uMt].wght[5],
        SrcPar[uMt].mean[6],SrcPar[uMt].stdv[6],SrcPar[uMt].wght[6],
        SrcPar[uMt].mean[7],SrcPar[uMt].stdv[7],SrcPar[uMt].wght[7],
        SrcPar[uMt].mean[8],SrcPar[uMt].stdv[8],SrcPar[uMt].wght[8],
        SrcPar[uMt].mean[9],SrcPar[uMt].stdv[9],
        Chi2Ndf[uMt]
      );
  }


  for(unsigned uLine=0; uLine<MaxNumLines; uLine++){
    delete [] LinesToSave[uLine];
  }
  delete [] LinesToSave;

  delete [] cline;
  delete [] cdscr;
  delete [] cval;

  root_file_core.Close();
  root_file_full.Close();
  delete [] MtBinCenter;
  delete [] SrcPar;
  delete [] Chi2Ndf;
  return 0;
}







void TestReadWriteBinary(){
  TString FileName = GetFemtoOutputFolder();
  FileName += "/BinTest.bin";

  unsigned ThatMany = 7;
  unsigned* Element = new unsigned[ThatMany];

  ofstream obfile (FileName, ios::app | ios::binary);
  for(unsigned uEl=0; uEl<ThatMany; uEl++){
    Element[uEl] = uEl*3+1;
    obfile.write((char *) &Element[uEl], sizeof(unsigned));
  }
  obfile.close();

  unsigned HowMany = 0;
  FileName += "rubbish";
  ifstream ibfile (FileName, ios::in | ios::binary);
  if(ibfile.fail()){
    printf("RUBBISH\n");
    return;
  }
  unsigned element;

  ibfile.seekg (0, ibfile.end);
  int length = ibfile.tellg();
  ibfile.seekg (0, ibfile.beg);
  if(length%sizeof(unsigned)){
    printf("wtf\n");
  }
  HowMany = length/sizeof(unsigned);

  printf("length = %i\n",length);
  printf("HowMany = %u\n",HowMany);
  //while(!ibfile.eof()){
  //while(ibfile){
  for(unsigned uEl=0; uEl<HowMany; uEl++){
    //HowMany++;
    ibfile.read((char *) &element, sizeof(unsigned));
    printf("%u: %u (%i)\n",uEl,element,ibfile.eof());
  }

  delete [] Element;
}

void Test_src_read(){
  //'Cigar2_ds'+str(d_num)+'_hts'+str(ht_num)+'_hzs'+str(hz_num)
  DLM_CecaSource_v0 dlmsrc1("pp","Cigar2_ds24_hts36_hzs36","fld");

}

void LnGammaTest(){
  double InputVal = 3.18;
  double root_val = TMath::LnGamma(InputVal);
  double dlm_val = gammln(InputVal);
  printf("root vs dlm = %.6e vs %.6e\n", root_val, dlm_val);

  double root_poiss = TMath::Poisson(3.1,3.8);
  double dlm_poiss = DLM_Poisson(3.1,3.8);
  printf("root_poiss vs dlm_poiss = %.6e vs %.6e\n", root_poiss, dlm_poiss);
}

void TestReadCigar2(){
  DLM_CecaSource_v0 pp_src("pp","Cigar2_ds24_hts36_hzs36","/home/dimihayl/Software/LocalFemto/Output/CECA_Paper/Output_260223/CecaSource/");

  //double rad = 1.0;
  double pars[4];
  //pars[0] = (1110.+1168.)*0.5;
  pars[0] = 1110;
  pars[1] = 0.1867*2.0;
  pars[2] = 4.067*1.05;
  pars[3] = 10.167*1.05;
  //KdpPars SrcPars = pp_src.RootEval(&rad,pars);
  //SrcPars.Print();
//return;
  TF1* mySource = new TF1("pp_src",&pp_src,&DLM_CecaSource_v0::RootEval,0,16,4,"DLM_CecaSource_v0","RootEval");

  mySource->SetParameter(0,pars[0]);
  mySource->SetParameter(1,pars[1]);
  mySource->SetParameter(2,pars[2]);
  mySource->SetParameter(3,pars[3]);

  TFile fOutput(TString::Format("%s/CECA_Paper/TestReadCigar2/fOutput.root",GetFemtoOutputFolder()),"recreate");
  mySource->SetNpx(4096);
  mySource->Write();

  delete mySource;
}


void TestCigar2_Ck(){
  DLM_CecaSource_v0 pp_src("pp","Cigar2_ds24_hts36_hzs36","/home/dimihayl/Software/LocalFemto/Output/CECA_Paper/Output_260223/CecaSource/");
  DLM_CecaSource_v0 pL_src("pL","Cigar2_ds24_hts36_hzs36","/home/dimihayl/Software/LocalFemto/Output/CECA_Paper/Output_260223/CecaSource/");


  const double mT_val = 1350;
  const double Dist = 0.5;
  const double HadrT = 2.0;
  const double HadrZ = 10.0;

  const double rG_pp = 1.2;
  const double rG_pL = 1.4;

  const unsigned NumMomBins = 80;
  const double kMin = 0;
  const double kMax = 320;

  CATS KittyG_pp;
  CATS KittyG_pL;
  CATS KittyC_pp;
  CATS KittyC_pL;

  KittyG_pp.SetMomBins(NumMomBins,kMin,kMax);
  KittyG_pL.SetMomBins(NumMomBins,kMin,kMax);
  KittyC_pp.SetMomBins(NumMomBins,kMin,kMax);
  KittyC_pL.SetMomBins(NumMomBins,kMin,kMax);

  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());


  AnalysisObject.SetUpCats_pp(KittyG_pp,"AV18","Gauss",0,202);
  AnalysisObject.SetUpCats_pL(KittyG_pL,"Chiral_Coupled_SPD","Gauss",11600,202);

  AnalysisObject.SetUpCats_pp(KittyC_pp,"AV18","",0,202);
  AnalysisObject.SetUpCats_pL(KittyC_pL,"Chiral_Coupled_SPD","",11600,202);


  KittyG_pp.SetAnaSource(0, rG_pp);
  KittyG_pL.SetAnaSource(0, rG_pL);

  KittyC_pp.SetUseAnalyticSource(true);
  KittyC_pp.SetAnaSource(CatsSourceForwarder, &pp_src, 4);
  KittyC_pp.SetAnaSource(0, mT_val);
  KittyC_pp.SetAnaSource(1, Dist);
  KittyC_pp.SetAnaSource(2, HadrT);
  KittyC_pp.SetAnaSource(3, HadrZ);

  KittyC_pL.SetUseAnalyticSource(true);
  KittyC_pL.SetAnaSource(CatsSourceForwarder, &pL_src, 4);
  KittyC_pL.SetAnaSource(0, mT_val);
  KittyC_pL.SetAnaSource(1, Dist);
  KittyC_pL.SetAnaSource(2, HadrT);
  KittyC_pL.SetAnaSource(3, HadrZ);

  KittyG_pp.KillTheCat();
  KittyG_pL.KillTheCat();
  KittyC_pp.KillTheCat();
  KittyC_pL.KillTheCat();

  TFile fOutput(TString::Format("%s/CECA_Paper/TestCigar2_Ck/fOutput.root",GetFemtoOutputFolder()),"recreate");

  TGraph gCk_pp_Gauss;
  TGraph gCk_pL_Gauss;
  TGraph gCk_pp_Ceca;
  TGraph gCk_pL_Ceca;


  gCk_pp_Gauss.SetName("Ck_pp_Gauss");
  gCk_pL_Gauss.SetName("Ck_pL_Gauss");
  gCk_pp_Ceca.SetName("Ck_pp_Ceca");
  gCk_pL_Ceca.SetName("Ck_pL_Ceca");

  for(unsigned uBin=0; uBin<NumMomBins; uBin++){
    double Momentum = KittyG_pp.GetMomentum(uBin);
    gCk_pp_Gauss.SetPoint(uBin,Momentum,KittyG_pp.GetCorrFun(uBin));
    gCk_pL_Gauss.SetPoint(uBin,Momentum,KittyG_pL.GetCorrFun(uBin));
    gCk_pp_Ceca.SetPoint(uBin,Momentum,KittyC_pp.GetCorrFun(uBin));
    gCk_pL_Ceca.SetPoint(uBin,Momentum,KittyC_pL.GetCorrFun(uBin));
  }

  TF1* fSrcC_pp = new TF1("fSrcC_pp",&pp_src,&DLM_CecaSource_v0::RootEval,0,16,4,"DLM_CecaSource_v0","RootEval");
  fSrcC_pp->SetParameter(0,mT_val);
  fSrcC_pp->SetParameter(1,Dist);
  fSrcC_pp->SetParameter(2,HadrT);
  fSrcC_pp->SetParameter(3,HadrZ);

  TF1* fSrcC_pL = new TF1("fSrcC_pL",&pL_src,&DLM_CecaSource_v0::RootEval,0,16,4,"DLM_CecaSource_v0","RootEval");
  fSrcC_pL->SetParameter(0,mT_val);
  fSrcC_pL->SetParameter(1,Dist);
  fSrcC_pL->SetParameter(2,HadrT);
  fSrcC_pL->SetParameter(3,HadrZ);

  TF1* fSrcG_pp = new TF1("fSrcG_pp",GaussSourceTF1,0,16,1);
  fSrcG_pp->SetParameter(0,rG_pp);
  TF1* fSrcG_pL = new TF1("fSrcG_pL",GaussSourceTF1,0,16,1);
  fSrcG_pL->SetParameter(0,rG_pL);

  fSrcC_pp->SetNpx(4096);
  fSrcC_pL->SetNpx(4096);
  fSrcG_pp->SetNpx(4096);
  fSrcG_pL->SetNpx(4096);


  gCk_pp_Gauss.Write();
  fSrcG_pp->Write();
  gCk_pp_Ceca.Write();
  fSrcC_pp->Write();
  gCk_pL_Gauss.Write();
  fSrcG_pL->Write();
  gCk_pL_Ceca.Write();
  fSrcC_pL->Write();

}


void PNG_ANIM_SOURCE(){
  DLM_CecaSource_v0 pp_src("pp","Cigar2_ds24_hts36_hzs36","/home/dimihayl/Software/LocalFemto/Output/CECA_Paper/Output_260223/CecaSource/");
  const double min_mt = 1000;
  const double max_mt = 2300;
  const double step_mt = (max_mt-min_mt)/52.;
  double mt_val = min_mt;

  const double Dist = 0.5;
  const double HadrT = 2.0;
  const double HadrZ = 10.0;

  while(mt_val<=max_mt+0.1*step_mt){

    TF1* fSrcC_pp = new TF1("fSrcC_pp",&pp_src,&DLM_CecaSource_v0::RootEval,0,8,4,"DLM_CecaSource_v0","RootEval");
    fSrcC_pp->SetParameter(0,mt_val);
    fSrcC_pp->SetParameter(1,Dist);
    fSrcC_pp->SetParameter(2,HadrT);
    fSrcC_pp->SetParameter(3,HadrZ);
    fSrcC_pp->SetLineColor(kBlue);
    fSrcC_pp->SetNpx(2048);

    TCanvas* can1 = new TCanvas("can1", "can1", 1);
    can1->cd(0); can1->SetCanvasSize(1280, 720); can1->SetMargin(0.15,0.05,0.2,0.05);//lrbt

    TH1F* hAxis = new TH1F("hAxis", "hAxis", 128, 0, 8);
    hAxis->SetStats(false);
    hAxis->SetTitle("");
    hAxis->GetXaxis()->SetTitle("r* (fm)");
    hAxis->GetXaxis()->SetTitleSize(0.06);
    hAxis->GetXaxis()->SetLabelSize(0.06);
    hAxis->GetXaxis()->CenterTitle();
    hAxis->GetXaxis()->SetTitleOffset(1.3);
    hAxis->GetXaxis()->SetLabelOffset(0.02);

    hAxis->GetYaxis()->SetTitle("4#pir^{2}S(r*)");
    hAxis->GetYaxis()->SetTitleSize(0.06);
    hAxis->GetYaxis()->SetLabelSize(0.06);
    hAxis->GetYaxis()->CenterTitle();
    hAxis->GetYaxis()->SetTitleOffset(1.10);

    //hAxisPot->GetXaxis()->SetNdivisions(506);
    hAxis->GetYaxis()->SetRangeUser(0, 0.4);

    TLegend* leg1 = new TLegend(0.65,0.75,0.95,0.95);//lbrt
    leg1->SetName("leg1");
    leg1->SetTextSize(0.055);
    leg1->AddEntry(fSrcC_pp, TString::Format("m_{T} = %.2f GeV",mt_val*0.001));
    leg1->Draw("same");

    hAxis->Draw("axis");
    fSrcC_pp->Draw("same");
    leg1->Draw("same");

    if(mt_val==min_mt){
      can1->SaveAs(TString::Format("%s/CECA_Paper/TestCigar2_Ck/anim.gif",GetFemtoOutputFolder()));
    }
    else{
      can1->SaveAs(TString::Format("%s/CECA_Paper/TestCigar2_Ck/anim.gif+8",GetFemtoOutputFolder()));
    }

    mt_val += step_mt;

    delete hAxis;
    delete fSrcC_pp;
    delete leg1;
    delete can1;

  }
}













//SmearStrategy == 1 is the new one
void CecaPaper_pp(const TString SourceVar, const int& SmearStrategy, const unsigned DataVar, const TString OutputFolder){

  gROOT->ProcessLine("gErrorIgnoreLevel = 2001;");
  const bool Silent = true;
  //for the pp source, the error on the sign
  const bool CorrectSign = true;
  const unsigned NumMtBins = 7;

  std::vector<int> pp_lam_var = {0,1,2};//done
  std::vector<int> pL_lam_var = {0};

  //0 is the old one, 1 is the new one with the folded ME
  //std::vector<int> SmearStrategy = {0,1};//done
  std::vector<float> FemtoRegion = {376,352,400};//done
  const unsigned NumMomBins_pp = TMath::Nint(FemtoRegion.back()/4.);
  const unsigned NumMomBins_feed = TMath::Nint(FemtoRegion.back()/10.);
  std::vector<float> BaselineRegion = {500};
  std::vector<float> CkCutOff = {700};//done
  std::vector<int> pL_pot_var = {11600,-11600};//done

  //pp13TeV_HM_DimiJun20 -> the ME reso matrix
  //pp13TeV_HM_BernieSource -> Bernie's old smearing matrix
  //pp13TeV_HM_BernieSource -> contains the data
  std::vector<TString> MomSmearVar = {"pp13TeV_HM_DimiJun20"};//done

  //and an additional variation will be done on if pS0 is included as a feed or not
  enum BLTYPE { pol0s,pol1s,pol2s,pol3s,dpol2s,dpol3s,dpol4s,pol2e,pol3e,dpol2e,dpol3e,dpol4e,spl1 };
  //std::vector<int> BaselineVar = {pol0s,pol1s,dpol3e};
  std::vector<int> BaselineVar = {dpol3e};
  //perhaps change to the new one
  const bool pS0_Var = true;

  //all histos to fit
  TH1F** hCk = new TH1F* [NumMtBins];
  //all fit functions
  TF1** fCk = new TF1* [NumMtBins];

  for(unsigned uMt=0; uMt<NumMtBins; uMt++){

  }


/*
  //the region for which the DLM_Ck objects will be defined
  //we put some extra to get away from the edge effects of the smearing
  const double kMin = 0;
  const double kMax = 600;
  const unsigned TotMomBins = TMath::Nint(kMax);
  const unsigned NumMtBins = 7;
  if(imTbin>=NumMtBins){
    printf("\033[1;31mERROR:\033[0m Only %u mT bins are available!\n",NumMtBins);
    return;
  }
  //approximate, based on the output we kind of know we will get
  //std::vector<float> ExpectedRadii = {1.35,1.3,1.25,1.2,1.1,1.05,0.95};
  //const double Scale_pL = 1.1;
  //const double Scale_pS0 = 1.15;
  //const double Scale_pXi = 0.97;
  //const double Scale_core = 0.94;
  std::vector<float> ExpectedRadii = { 1.55, 1.473, 1.421, 1.368, 1.295, 1.220,1.124 };
  //std::vector<float> ExpectedRadii = {1.35,1.3,1.25,1.2,1.1,1.05,0.95};
  const double Scale_pL = 1.0;
  const double Scale_pS0 = 1.0;
  const double Scale_pXi = 1.0;
  const double Scale_core = 1.0;
  //the difference in the effectiv Gaussian compered to pp
//printf("1\n");
  const unsigned NumSourcePars = 1;
  //std::vector<float> pSigma0Radii = {1.55,1.473,1.421,1.368,1.295,1.220,1.124};

  DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());

  //gROOT->cd();
  TFile* fOutputFile = new TFile(OutputFolder+TString::Format("/fOut_%s_SS%i_D%u_mT%i_ppSign%i.root",SourceVar.Data(),SmearStrategy,DataVar,imTbin,CorrectSign),"recreate");

  TH1F* HDATA = NULL;
  //GDATA->Set(TotMomBins);
  TGraph* GFIT = new TGraph();
  //GFIT->Set(TotMomBins);
  TGraph* GBL = new TGraph();
  //GBL->Set(TotMomBins);
  TGraph* GFEMTO = new TGraph();
  //GFEMTO->Set(TotMomBins);
  Int_t BASELINEVAR;
  Float_t FEMTOREGION;
  Int_t PS0;
  Int_t SMEARSTRATEGY;
  Int_t PP_LAM_VAR;
  TString MOMSMEARVAR;
  Float_t CKCUTOFF;
  Int_t PL_POT_VAR;
  TString SOURCEVAR;
  Float_t LAM_PP;
  Float_t LAM_PPL;
  Float_t NSIG;
  Float_t RADIUS;
  Float_t RADERR;

  TTree* ppTree = new TTree("ppTree","ppTree");
  ppTree->Branch("hData","TH1F",&HDATA,32000,0);//
  ppTree->Branch("gFit","TGraph",&GFIT,32000,0);//
  ppTree->Branch("gBl","TGraph",&GBL,32000,0);//
  ppTree->Branch("gFemto","TGraph",&GFEMTO,32000,0);//
  ppTree->Branch("BaselineVar",&BASELINEVAR,"BaselineVar/I");//
  ppTree->Branch("FemtoRegion",&FEMTOREGION,"FemtoRegion/F");//
  ppTree->Branch("pS0",&PS0,"pS0/I");//
  ppTree->Branch("SmearStrategy",&SMEARSTRATEGY,"SmearStrategy/I");//
  ppTree->Branch("pp_lam_var",&PP_LAM_VAR,"pp_lam_var/I");//
  ppTree->Branch("MomSmearVar","TString",&MOMSMEARVAR,8000,0);//
  ppTree->Branch("CkCutOff",&CKCUTOFF,"CkCutOff/F");//
  ppTree->Branch("pL_pot_var",&PL_POT_VAR,"pL_pot_var/I");//
  ppTree->Branch("SourceVar","TString",&SOURCEVAR,8000,0);//
  ppTree->Branch("lam_pp",&LAM_PP,"lam_pp/F");//
  ppTree->Branch("lam_ppL",&LAM_PPL,"lam_ppL/F");//
  ppTree->Branch("nsig",&NSIG,"nsig/F");//
  ppTree->Branch("rad",&RADIUS,"rad/F");//
  ppTree->Branch("raderr",&RADERR,"raderr/F");//

//printf("2\n");
  //for(TString varSource : SourceVar){
    CATS AB_pp;
    AB_pp.SetMomBins(NumMomBins_pp,0,FemtoRegion.back());
    AnalysisObject.SetUpCats_pp(AB_pp,"AV18",SourceVar,0,CorrectSign?-202:202);
    AB_pp.SetAnaSource(0,ExpectedRadii.at(imTbin));
    if(SourceVar.Contains("Levy")) AB_pp.SetAnaSource(1,1.7);
    else AB_pp.SetAnaSource(1,2.0);
    AB_pp.SetNotifications(CATS::nWarning);
    AB_pp.KillTheCat();
//printf("3\n");
    CATS AB_pXim;
    //same binning as pL, as we only use pXim as feed-down
    AB_pXim.SetMomBins(NumMomBins_feed,0,FemtoRegion.back());
    AnalysisObject.SetUpCats_pXim(AB_pXim,"pXim_HALQCDPaper2020","Gauss");
    AB_pXim.SetAnaSource(0,ExpectedRadii.at(imTbin)*Scale_pXi);
    AB_pXim.SetNotifications(CATS::nWarning);
    AB_pXim.KillTheCat();
//printf("4\n");
    CATS AB_pXi1530;
    AB_pXi1530.SetMomBins(NumMomBins_feed,0,FemtoRegion.back());
    AB_pXi1530.SetNotifications(CATS::nWarning);
    AnalysisObject.SetUpCats_pXim(AB_pXi1530,"pXim1530","Gauss");//McLevyNolan_Reso
    AB_pXi1530.SetAnaSource(0,ExpectedRadii.at(imTbin)*Scale_pXi);
    AB_pXi1530.KillTheCat();
//printf("5\n");
    CATS AB_pS0_Chiral;
    //the minus one is to to avoid going above 350 MeV, since we do not have the WF there
    AB_pS0_Chiral.SetMomBins(NumMomBins_feed,0,FemtoRegion.back());
    AnalysisObject.SetUpCats_pS0(AB_pS0_Chiral,"Chiral","Gauss");
    AB_pS0_Chiral.SetAnaSource(0,ExpectedRadii.at(imTbin)*Scale_pS0);
    AB_pS0_Chiral.SetNotifications(CATS::nWarning);
    AB_pS0_Chiral.KillTheCat();
//printf("6\n");
    for(int varPL : pL_pot_var){
      CATS AB_pL;
      AB_pL.SetMomBins(NumMomBins_feed,0,FemtoRegion.back());
      AnalysisObject.SetUpCats_pL(AB_pL,"Chiral_Coupled_SPD","Gauss",varPL,202);//NLO_Coupled_S
      const double CuspWeight = 0.33;//0.54
      if(abs(varPL)>1000){
          AB_pL.SetChannelWeight(7,1./4.*CuspWeight);//1S0 SN(s) -> LN(s)
          AB_pL.SetChannelWeight(8,3./4.*CuspWeight);//3S1 SN(s) -> LN(s)
          AB_pL.SetChannelWeight(10,3./4.*CuspWeight);//3S1 SN(d) -> LN(s)
          AB_pL.SetChannelWeight(13,3./20.*CuspWeight);//3D1 SN(d) -> LN(d)
          AB_pL.SetChannelWeight(15,3./20.*CuspWeight);//3D1 SN(s) -> LN(d)
      }
      AB_pL.SetAnaSource(0,ExpectedRadii.at(imTbin)*Scale_pL);
      AB_pL.SetNotifications(CATS::nError);
      AB_pL.KillTheCat();

      for(float varCutOff : CkCutOff){
        DLM_Ck* Ck_pp = new DLM_Ck(NumSourcePars, 0, AB_pp, TotMomBins, kMin, kMax);
        Ck_pp->SetCutOff(FemtoRegion.back(),varCutOff);
        DLM_Ck* Ck_pL = new DLM_Ck(NumSourcePars, 0, AB_pL, TotMomBins, kMin, kMax);
        Ck_pL->SetCutOff(FemtoRegion.at(0),varCutOff);
        DLM_Ck* Ck_pS0 = new DLM_Ck(NumSourcePars, 0, AB_pS0_Chiral, TotMomBins, kMin, kMax);
        Ck_pS0->SetCutOff(FemtoRegion.at(0),varCutOff);
        DLM_Ck* Ck_pXim = new DLM_Ck(NumSourcePars, 0, AB_pXim, TotMomBins, kMin, kMax);
        Ck_pXim->SetCutOff(FemtoRegion.at(0),varCutOff);
        DLM_Ck* Ck_pXim1530 = new DLM_Ck(NumSourcePars, 0, AB_pXi1530, TotMomBins, kMin, kMax);
        Ck_pXim1530->SetCutOff(FemtoRegion.at(0),varCutOff);

        for(TString varSmear : MomSmearVar){
          TH2F* hResolution_pp = AnalysisObject.GetResolutionMatrix(varSmear,"pp");
          TH2F* hResidual_pp_pL = AnalysisObject.GetResidualMatrix("pp","pLambda");
          TH2F* hResidual_pL_pSigma0 = AnalysisObject.GetResidualMatrix("pLambda","pSigma0");
          TH2F* hResidual_pL_pXim = AnalysisObject.GetResidualMatrix("pLambda","pXim");
          TH2F* hResidual_pXi_pXi1530 = AnalysisObject.GetResidualMatrix("pXim","pXim1530");

          for(int varLam : pp_lam_var){
            double lambda_pp[4];
            AnalysisObject.SetUpLambdaPars_pp("pp13TeV_HM_BernieSource",varLam+imTbin*10,lambda_pp);
            double lambda_pL[5];
            AnalysisObject.SetUpLambdaPars_pL("pp13TeV_HM_Dec19",varLam+imTbin*10,0,lambda_pL);
            double lambda_pXim[5];
            AnalysisObject.SetUpLambdaPars_pXim("pp13TeV_HM_Dec19",varLam+imTbin*10,0,lambda_pXim);
//printf("SetUp\n");
            //for(int varSS : SmearStrategy){
              TH1F* hPhaseSpace_pp=NULL;
              if(SmearStrategy==1){
                TList* list1_tmp;
                TList* list2_tmp;
                TString FileName = TString::Format("%s/CatsFiles/ExpData/Bernie_Source/ppData/mTBin_%i/CFOutput_mT_ppVar%u_HM_%i.root",
                                                    GetCernBoxDimi(),imTbin+1,DataVar,imTbin);
                TFile* inFile = new TFile(FileName,"read");
                //PARTICLES
                list1_tmp = (TList*)inFile->Get("PairDist");
                list2_tmp = (TList*)list1_tmp->FindObject("PairFixShifted");//there is also PairShifted ?? which one ??
                TH1F* hME_PP = (TH1F*)list2_tmp->FindObject(TString::Format("MEPart_mT_%i_FixShifted",imTbin));
                list1_tmp = (TList*)inFile->Get("AntiPairDist");
                list2_tmp = (TList*)list1_tmp->FindObject("PairFixShifted");
                TH1F* hME_APAP = (TH1F*)list2_tmp->FindObject(TString::Format("MEAntiPart_mT_%i_FixShifted",imTbin));
                gROOT->cd();
                hPhaseSpace_pp = (TH1F*)hME_PP->Clone("hPhaseSpace_pp");
                hPhaseSpace_pp->Add(hME_APAP);
                delete inFile;
              }
//printf("GotData\n");
              for(unsigned ipS0=!(pS0_Var); ipS0<2; ipS0++){
              //for(unsigned ipS0=0; ipS0<2; ipS0++){
                DLM_CkDecomposition CkDec_pp("pp",3,*Ck_pp,hResolution_pp);
                DLM_CkDecomposition CkDec_pL("pLambda", 4,*Ck_pL,NULL);
                DLM_CkDecomposition CkDec_pSigma0("pSigma0",0,*Ck_pS0,NULL);
                DLM_CkDecomposition CkDec_pXim("pXim",3,*Ck_pXim,NULL);
                DLM_CkDecomposition CkDec_pXim1530("pXim1530",0,*Ck_pXim1530,NULL);

                CkDec_pp.AddContribution(0,lambda_pp[1],DLM_CkDecomposition::cFeedDown,&CkDec_pL,hResidual_pp_pL);
                CkDec_pp.AddContribution(1,lambda_pp[2],DLM_CkDecomposition::cFeedDown);
                CkDec_pp.AddContribution(2,lambda_pp[3],DLM_CkDecomposition::cFake);
                if(hPhaseSpace_pp){
                  CkDec_pp.AddPhaseSpace(hPhaseSpace_pp);
                  CkDec_pp.AddPhaseSpace(0, hPhaseSpace_pp);
                }

                if(ipS0==0) CkDec_pL.AddContribution(0,lambda_pL[1],DLM_CkDecomposition::cFeedDown);
                else CkDec_pL.AddContribution(0,lambda_pL[1],DLM_CkDecomposition::cFeedDown,&CkDec_pSigma0,hResidual_pL_pSigma0);
                CkDec_pL.AddContribution(1, lambda_pL[2],DLM_CkDecomposition::cFeedDown,&CkDec_pXim,hResidual_pL_pXim);
                CkDec_pL.AddContribution(2, lambda_pL[3],DLM_CkDecomposition::cFeedDown);
                CkDec_pL.AddContribution(3, lambda_pL[4],DLM_CkDecomposition::cFake);
                if(hPhaseSpace_pp){
                  if(ipS0==1) CkDec_pL.AddPhaseSpace(0,hPhaseSpace_pp);
                  CkDec_pL.AddPhaseSpace(1,hPhaseSpace_pp);
                }

                CkDec_pXim.AddContribution(0, lambda_pXim[1],DLM_CkDecomposition::cFeedDown,&CkDec_pXim1530,hResidual_pXi_pXi1530);  //from Xi-(1530)
                CkDec_pXim.AddContribution(1, lambda_pXim[2]+lambda_pXim[3],DLM_CkDecomposition::cFeedDown);  //other feed-down (flat)
                CkDec_pXim.AddContribution(2, lambda_pXim[4], DLM_CkDecomposition::cFake);
                if(hPhaseSpace_pp){
                  CkDec_pXim.AddPhaseSpace(0,hPhaseSpace_pp);
                }

                CkDec_pp.Update();
                CkDec_pL.Update();
                CkDec_pSigma0.Update();
                CkDec_pXim.Update();
                CkDec_pXim1530.Update();
//printf("Updated\n");
                gROOT->cd();
                TH1F* hData = AnalysisObject.GetAliceExpCorrFun("pp13TeV_HM_BernieSource","pp",TString::Format("%u",DataVar),0,0,imTbin);
//printf("GetAliceExpCorrFun %p\n",hData);
                for(float varFit : FemtoRegion){
                  for(int varBL : BaselineVar){
                    double FitRange = varFit;
                    if(varBL==dpol3e) FitRange=BaselineRegion.at(0);
                    //gROOT->cd();
                    fOutputFile->cd();
                    //2 femto and 5 BL fit pars (the BL are norm, pol1,2,3,4)
                    TF1* fData = new TF1("fData",Fit_BernieSource,0,FitRange,7);
                    TF1* fBl = new TF1("fBl",Baseline_BernieSource,0,FitRange,5);
                    fData->SetParameter(0,SourceVar=="Gauss"?ExpectedRadii.at(imTbin):ExpectedRadii.at(imTbin)*Scale_core);
                    fData->SetParLimits(0,fData->GetParameter(0)*0.5,fData->GetParameter(0)*2.0);
                    //fData->SetParLimits(0,fData->GetParameter(0)*0.75,fData->GetParameter(0)*1.25);
//fData->FixParameter(0,1.3);
//fData->SetParameter(0,1.366);
//fData->SetParLimits(0,1.1,1.5);
                    fData->FixParameter(1,0);

                    //BL
                    fData->SetParameter(2,1);

                    //pol0s,pol1s,dpol3e
                    if(varBL==pol1s){
                      fData->SetParameter(3,0);
                      fData->SetParLimits(3,-1e-2,1e-2);
                      fData->FixParameter(4,0);
                      fData->FixParameter(5,0);
                      fData->FixParameter(6,0);
                    }
                    else if(varBL==dpol3e){
                      fData->SetParameter(3,0);
                      fData->SetParLimits(3,-100000,100);
                      fData->SetParameter(4,100);
                      fData->SetParLimits(4,0,400);
                      fData->SetParameter(5,0);
                      fData->SetParLimits(5,-1e-6,1e-6);
                      fData->FixParameter(6,0);
                    }
                    else{
                      fData->FixParameter(3,0);
                      fData->FixParameter(4,0);
                      fData->FixParameter(5,0);
                      fData->FixParameter(6,0);
                    }
                    SOURCE_FIT = &CkDec_pp;
                    printf("BL=%i FIT=%.0f PS0=%i SS=%i LAM=%i SMR=%s PL=%i SRC=%s lam_pp=%.1f lam_ppl=%.1f\n",
                            varBL,varFit,ipS0,SmearStrategy,varLam,varSmear.Data(),varPL,SourceVar.Data(),lambda_pp[0]*100.,lambda_pp[1]*100.);
                    hData->Fit(fData,"Q, S, N, R, M");
//fData->FixParameter(0,1.4);
//fData->FixParameter(2,1);
//TFile TempF(OutputFolder+"/TempF.root","recreate");
//hData->Write();
//fData->Write();
//CkDec_pp.Update(true,true);
//printf("C(20) = %f; cats %f; Ck_pp %f; CkDec_pp %f\n",fData->Eval(20),AB_pp.EvalCorrFun(20),Ck_pp->Eval(20),CkDec_pp.EvalCk(20));
//return;
//printf("fitted\n");
                    fData->SetNpx(TotMomBins);

                    fBl->FixParameter(0,fData->GetParameter(2));
                    fBl->FixParameter(1,fData->GetParameter(3));
                    fBl->FixParameter(2,fData->GetParameter(4));
                    fBl->FixParameter(3,fData->GetParameter(5));
                    fBl->FixParameter(4,fData->GetParameter(6));

printf(" r = %.3f +/- %.3f\n",fData->GetParameter(0),fData->GetParError(0));
                    double Chi2 = fData->GetChisquare();
                    double NDF = fData->GetNDF();
                    double pval = TMath::Prob(Chi2,NDF);
                    double nsig = sqrt(2)*TMath::ErfcInverse(pval);
                    printf(" chi2/ndf = %.2f\n",Chi2/NDF);
                    printf(" nsig = %.2f\n",nsig);

                    //gROOT->cd();
                    fOutputFile->cd();

                    BASELINEVAR = varBL;
                    FEMTOREGION = varFit;
                    PS0 = ipS0;
                    SMEARSTRATEGY = SmearStrategy;
                    PP_LAM_VAR = varLam;
                    MOMSMEARVAR = varSmear;
                    CKCUTOFF = varCutOff;
                    PL_POT_VAR = varPL;
                    SOURCEVAR = SourceVar;
                    LAM_PP = lambda_pp[0];
                    LAM_PPL = lambda_pp[1];
                    NSIG = nsig;
                    RADIUS = fData->GetParameter(0);
                    RADERR = fData->GetParError(0);

                    if(HDATA)delete HDATA;
                    HDATA = (TH1F*)hData->Clone("HDATA");
                    delete GFIT;
                    GFIT = new TGraph();
                    GFIT->SetName("GFIT");
                    delete GBL;
                    GBL = new TGraph();
                    GBL->SetName("GBL");
                    delete GFEMTO;
                    GFEMTO = new TGraph();
                    GFEMTO->SetName("GFEMTO");
                    for(unsigned uBin=0; uBin<NumMomBins_pp; uBin++){
                      double MOM = AB_pp.GetMomentum(uBin);
                      if(MOM>FitRange) break;
                      //GDATA->SetPoint(uBin,MOM,hData->GetBinContent(hData->FindBin(MOM)));
                      GFIT->SetPoint(uBin,MOM,fData->Eval(MOM));
                      GBL->SetPoint(uBin,MOM,fBl->Eval(MOM));
                      GFEMTO->SetPoint(uBin,MOM,fData->Eval(MOM)/fBl->Eval(MOM));
                    }

                    ppTree->Fill();
                    //ppTree->Branch("gData","TH1F",&GDATA,32000,0);//
                    //ppTree->Branch("gFit","TGraph",&GFIT,32000,0);//
                    //ppTree->Branch("gBl","TGraph",&GBL,32000,0);//
                    //ppTree->Branch("gFemto","TGraph",&GFEMTO,32000,0);//

//save output here: OutputFolder

                    TFile outTest(OutputFolder+"/outTest_OLD.root","recreate");

                    TGraph* GFEMTO_TH = new TGraph();
                    TGraph* GFEMTO_SF = new TGraph();
                    GFEMTO_TH->SetName("GFEMTO_TH");
                    GFEMTO_SF->SetName("GFEMTO_SF");
                    CkDec_pp.GetCk()->SetSourcePar(0,RADIUS);
                    CkDec_pp.Update(true,true);
                    for(unsigned uBin=0; uBin<NumMomBins_pp; uBin++){
                      double MOM = AB_pp.GetMomentum(uBin);
                      if(MOM>FitRange) break;
                      GFEMTO_TH->SetPoint(uBin,MOM,CkDec_pp.EvalSignal(MOM)+1);
                      GFEMTO_SF->SetPoint(uBin,MOM,SOURCE_FIT->EvalCk(MOM)+1);
                    }

                    //printf("File opened\n");
                    hData->Write();
                    //printf("Histo written\n");
                    GFIT->Write();
                    GBL->Write();
                    GFEMTO->Write();
                    GFEMTO_TH->Write();
                    GFEMTO_SF->Write();
                    delete GFEMTO_TH;


                    //return;
                    delete fData;
                    delete fBl;
//break;
                  }//varBL (3x)
                }//varFit (3x)
                delete hData;
//break;
              }//ipS0 (2x)
            //}//varSS (2x)
//break;
          }//varLam (3x)
//break;
        }//varSmear (2x)
        delete Ck_pp;
        delete Ck_pL;
        delete Ck_pS0;
        delete Ck_pXim;
        delete Ck_pXim1530;
//break;
      }//varCutOff (1x)
//break;
    }//varPL (2x)
  //}//varSource (2x)

//const unsigned DataVar, const int imTbin
  fOutputFile->cd();
  ppTree->Write();




  for(unsigned uMt=0; uMt<NumMtBins; uMt++){
    if(hCk[uMt]){delete hCk[uMt]; hCk[uMt]=NULL;}
    if(fCk[uMt]){delete fCk[uMt]; fCk[uMt]=NULL;}
  }
  delete [] hCk;
  delete [] fCk;




  delete ppTree;
  delete fOutputFile;

*/

}






CecaAnalysis1::CecaAnalysis1(TString SourceVersion, TString SourceFolder){
  hCkExp_pp = NULL;

  lam_pp_gen = NULL;
  lam_pp_flt = NULL;
  lam_pp_pL = NULL;
  lam_pp_fake = NULL;

  lam_pL_gen = NULL;
  lam_pL_flt = NULL;
  lam_pL_pS0 = NULL;
  lam_pL_pXi0 = NULL;
  lam_pL_pXim = NULL;
  lam_pL_fake = NULL;

  lam_Xi_gen = NULL;
  lam_Xi_Xi1530 = NULL;
  lam_Xi_flt = NULL;

  lam_pS0_gen = NULL;
  lam_pS0_flt = NULL;

  Model_pS0 = 0;
  Model_pL = 11600;

  FemtoRegMax_pp = 376;
  FitMax_pp = 376;
  FitCutOff_pp = 700;

  FemtoRegMax_pL = 336;
  FitMax_pL = 456;
  FitCutOff_pL = 700;

  Kitty_pp = NULL;
  Kitty_pL = NULL;
  Kitty_pS0 = NULL;
  Kitty_pXi0 = NULL;
  Kitty_pXim = NULL;

  Ck_pp = NULL;
  Ck_pL = NULL;
  Ck_pS0 = NULL;
  Ck_pXi0 = NULL;
  Ck_pXim = NULL;

  CkDec_pp = NULL;
  CkDec_pL = NULL;
  CkDec_pS0 = NULL;
  CkDec_pXi0 = NULL;
  CkDec_pXim = NULL;

  hReso_pp = NULL;
  hPs_pp = NULL;
  hReso_pL = NULL;
  hPs_pL = NULL;

  Src_pp = NULL;
  Src_pL = NULL;

  GrandeFitter = NULL;

  Src_pp = new DLM_CecaSource_v0("pp",SourceVersion.Data(),SourceFolder.Data());
  Src_pL = new DLM_CecaSource_v0("pp",SourceVersion.Data(),SourceFolder.Data());

  if(Src_pp->InErrorState()){
    printf("\033[1;31mERROR:\033[0m CecaAnalysis1 failed to load the pp source %s from %s.\n",SourceVersion.Data(),SourceFolder.Data());
    delete Src_pp;
    Src_pp = NULL;
  }
  if(Src_pL->InErrorState()){
    printf("\033[1;31mERROR:\033[0m CecaAnalysis1 failed to load the pL source %s from %s.\n",SourceVersion.Data(),SourceFolder.Data());
    delete Src_pL;
    Src_pL = NULL;
  }

}


CecaAnalysis1::~CecaAnalysis1(){

}






int CECA_PAPER(int argc, char *argv[]){
  printf("CECA_PAPER\n");
  //how to read/write the Levy pars into a file
  //TestSaveStuctToFile();

  //a test that extrapolation of Levy pars works
  //TestDoubleSourceOperation();

  //Test_pp_Statistics_1();

  //return Ceca_pp_or_pL("TEST4_pp",TString::Format("%s/CECA_Paper/Ceca_pp_or_pL/TEST4/Job/",GetFemtoOutputFolder()),
  //TString::Format("%s/CECA_Paper/Ceca_pp_or_pL/TEST4/Out/",GetFemtoOutputFolder()),
  //TString::Format("%s/CECA_Paper/Ceca_pp_or_pL/TEST4/Log/",GetFemtoOutputFolder()),
  //0,0,1);
  //ReadDlmHst();

  //dadd_f(argc,argv);

  //TestReadWriteBinary();

  //TString DADD_IN=TString(argv[1]);
  //TString DADD_OUT=TString(argv[2]);
  //TString ROOT_OUT=TString(argv[3]);
  //TString MASTER_FILE=TString(argv[4]);
  //TString FINAL_OUT=TString(argv[5]);
  char** myinput = new char* [6];
  for(unsigned uch=0; uch<6; uch++){
    myinput[uch] = new char [511];
  }
  //strcpy(myinput[1],"/home/dimihayl/Software/LocalFemto/Output/CECA_Paper/TestBatchOutput/test1/IN/Cigar1.26888.*.full.dlm.hst");
  //strcpy(myinput[1],"Cigar1.26888");
  strcpy(myinput[1],"Cigar1xs.18");
  strcpy(myinput[2],"/home/dimihayl/Software/LocalFemto/Output/CECA_Paper/TestBatchOutput/test1/IN/");
  strcpy(myinput[3],"/home/dimihayl/Software/LocalFemto/Output/CECA_Paper/TestBatchOutput/test1/IN/");
  strcpy(myinput[4],"/home/dimihayl/Software/LocalFemto/Output/CECA_Paper/TestBatchOutput/test1/");
  printf("Lets go\n");
  //return dlmhst_ceca_kdpfit_wc0(5,myinput);
  //Test_src_read();
  //LnGammaTest();
  //CreateDlmHistoCigar2("pp",TString::Format("%s/CECA_Paper/Output_260223/",GetFemtoOutputFolder()),"Cigar2");
  //CreateDlmHistoCigar2("pL",TString::Format("%s/CECA_Paper/Output_260223/",GetFemtoOutputFolder()),"Cigar2");
  //TestReadCigar2();
  //TestCigar2_Ck();
  PNG_ANIM_SOURCE();

  return 0;
}
