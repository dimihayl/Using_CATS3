
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
#include "DLM_CkDecomposition.h"
#include "DLM_Ck.h"
#include "CATS.h"
#include "DLM_MultiFit.h"

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
#include "TFitResultPtr.h"
#include "TFitResult.h"

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
  DLM_CecaSource_v0 pp_src("pp","Cigar2_ds24_hts36_hzs36",TString::Format("%s/CatsFiles/Source/CECA/LookUpSource/Cigar2_ds24_hts36_hzs36/pp/",GetCernBoxDimi()).Data());

  //double rad = 1.0;
  double pars[5];
  //pars[0] = (1110.+1168.)*0.5;
  pars[0] = 1110;
  pars[1] = 0.1867*2.0;
  pars[2] = 4.067*1.05;
  pars[3] = 10.167*1.05;
  pars[4] = 1;
  //KdpPars SrcPars = pp_src.RootEval(&rad,pars);
  //SrcPars.Print();
//return;
  TF1* mySource = new TF1("pp_src",&pp_src,&DLM_CecaSource_v0::RootEval,0,16,5,"DLM_CecaSource_v0","RootEval");

  mySource->SetParameter(0,pars[0]);
  mySource->SetParameter(1,pars[1]);
  mySource->SetParameter(2,pars[2]);
  mySource->SetParameter(3,pars[3]);
  mySource->SetParameter(4,pars[4]);

  TFile fOutput(TString::Format("%s/CECA_Paper/TestReadCigar2/fOutput.root",GetFemtoOutputFolder()),"recreate");
  mySource->SetNpx(4096);
  mySource->Write();

  delete mySource;
}

void TestReadJaime1(){
  DLM_CecaSource_v0 pp_src("pp","Jaime1_ds24_hts36_hzs36",TString::Format("%s/CatsFiles/Source/CECA/LookUpSource/Jaime1_ds24_hts36_hzs36/pp/",GetCernBoxDimi()).Data());

  //double rad = 1.0;
  double pars[5];
  //pars[0] = (1110.+1168.)*0.5;
  pars[0] = 1200;
  pars[1] = 0.2;
  pars[2] = 2.35;
  pars[3] = 3.55;
  pars[4] = 1;
  //KdpPars SrcPars = pp_src.RootEval(&rad,pars);
  //SrcPars.Print();
//return;
  TF1* mySource = new TF1("pp_src",&pp_src,&DLM_CecaSource_v0::RootEval,0,16,5,"DLM_CecaSource_v0","RootEval");

  mySource->SetParameter(0,pars[0]);
  mySource->SetParameter(1,pars[1]);
  mySource->SetParameter(2,pars[2]);
  mySource->SetParameter(3,pars[3]);
  mySource->SetParameter(4,pars[4]);

  TFile fOutput(TString::Format("%s/CECA_Paper/TestReadJaime1/fOutput.root",GetFemtoOutputFolder()),"recreate");
  mySource->SetNpx(4096);
  mySource->Write();

  delete mySource;
}



void TestCigar2_Ck(){
  DLM_CecaSource_v0 pp_src("pp","Cigar2_ds24_hts36_hzs36",TString::Format("%s/CatsFiles/Source/CECA/LookUpSource/Cigar2_ds24_hts36_hzs36/pp/",GetCernBoxDimi()).Data());
  DLM_CecaSource_v0 pL_src("pL","Cigar2_ds24_hts36_hzs36",TString::Format("%s/CatsFiles/Source/CECA/LookUpSource/Cigar2_ds24_hts36_hzs36/pL/",GetCernBoxDimi()).Data());


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
  KittyC_pp.SetAnaSource(CatsSourceForwarder, &pp_src, 5);
  KittyC_pp.SetAnaSource(0, mT_val);
  KittyC_pp.SetAnaSource(1, Dist);
  KittyC_pp.SetAnaSource(2, HadrT);
  KittyC_pp.SetAnaSource(3, HadrZ);
  KittyC_pp.SetAnaSource(4, 1);

  KittyC_pL.SetUseAnalyticSource(true);
  KittyC_pL.SetAnaSource(CatsSourceForwarder, &pL_src, 5);
  KittyC_pL.SetAnaSource(0, mT_val);
  KittyC_pL.SetAnaSource(1, Dist);
  KittyC_pL.SetAnaSource(2, HadrT);
  KittyC_pL.SetAnaSource(3, HadrZ);
  KittyC_pL.SetAnaSource(4, 1);

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

  TF1* fSrcC_pp = new TF1("fSrcC_pp",&pp_src,&DLM_CecaSource_v0::RootEval,0,16,5,"DLM_CecaSource_v0","RootEval");
  fSrcC_pp->SetParameter(0,mT_val);
  fSrcC_pp->SetParameter(1,Dist);
  fSrcC_pp->SetParameter(2,HadrT);
  fSrcC_pp->SetParameter(3,HadrZ);
  fSrcC_pp->FixParameter(4,1);

  TF1* fSrcC_pL = new TF1("fSrcC_pL",&pL_src,&DLM_CecaSource_v0::RootEval,0,16,5,"DLM_CecaSource_v0","RootEval");
  fSrcC_pL->SetParameter(0,mT_val);
  fSrcC_pL->SetParameter(1,Dist);
  fSrcC_pL->SetParameter(2,HadrT);
  fSrcC_pL->SetParameter(3,HadrZ);
  fSrcC_pL->FixParameter(4,1);

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




void TestJaime1_Ck(){
  DLM_CecaSource_v0 pp_src("pp","Jaime1_ds24_hts36_hzs36",TString::Format("%s/CatsFiles/Source/CECA/LookUpSource/Jaime1_ds24_hts36_hzs36/pp/",GetCernBoxDimi()).Data());
  DLM_CecaSource_v0 pL_src("pL","Jaime1_ds24_hts36_hzs36",TString::Format("%s/CatsFiles/Source/CECA/LookUpSource/Jaime1_ds24_hts36_hzs36/pL/",GetCernBoxDimi()).Data());


  const double mT_val = 1350;
  const double Dist = 0.2;
  const double HadrT = 2.3;
  const double Tau = 3.55;

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
  KittyC_pp.SetAnaSource(CatsSourceForwarder, &pp_src, 5);
  KittyC_pp.SetAnaSource(0, mT_val);
  KittyC_pp.SetAnaSource(1, Dist);
  KittyC_pp.SetAnaSource(2, HadrT);
  KittyC_pp.SetAnaSource(3, Tau);
  KittyC_pp.SetAnaSource(4, 1);

  KittyC_pL.SetUseAnalyticSource(true);
  KittyC_pL.SetAnaSource(CatsSourceForwarder, &pL_src, 5);
  KittyC_pL.SetAnaSource(0, mT_val);
  KittyC_pL.SetAnaSource(1, Dist);
  KittyC_pL.SetAnaSource(2, HadrT);
  KittyC_pL.SetAnaSource(3, Tau);
  KittyC_pL.SetAnaSource(4, 1);

  KittyG_pp.KillTheCat();
  KittyG_pL.KillTheCat();
  KittyC_pp.KillTheCat();
  KittyC_pL.KillTheCat();

  TFile fOutput(TString::Format("%s/CECA_Paper/TestJaime1_Ck/fOutput.root",GetFemtoOutputFolder()),"recreate");

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

  TF1* fSrcC_pp = new TF1("fSrcC_pp",&pp_src,&DLM_CecaSource_v0::RootEval,0,16,5,"DLM_CecaSource_v0","RootEval");
  fSrcC_pp->SetParameter(0,mT_val);
  fSrcC_pp->SetParameter(1,Dist);
  fSrcC_pp->SetParameter(2,HadrT);
  fSrcC_pp->SetParameter(3,Tau);
  fSrcC_pp->FixParameter(4,1);

  TF1* fSrcC_pL = new TF1("fSrcC_pL",&pL_src,&DLM_CecaSource_v0::RootEval,0,16,5,"DLM_CecaSource_v0","RootEval");
  fSrcC_pL->SetParameter(0,mT_val);
  fSrcC_pL->SetParameter(1,Dist);
  fSrcC_pL->SetParameter(2,HadrT);
  fSrcC_pL->SetParameter(3,Tau);
  fSrcC_pL->FixParameter(4,1);

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
  //DLM_CecaSource_v0 pp_src("pp","Cigar2_ds24_hts36_hzs36","/home/dimihayl/Software/LocalFemto/Output/CECA_Paper/Output_260223/CecaSource/");
  DLM_CecaSource_v0 pp_src("pp","Cigar2_ds24_hts36_hzs36",TString::Format("%s/CatsFiles/Source/CECA/LookUpSource/Cigar2_ds24_hts36_hzs36/pp/",GetCernBoxDimi()).Data());
  const double min_mt = 1000;
  const double max_mt = 2300;
  const double step_mt = (max_mt-min_mt)/52.;
  double mt_val = min_mt;

  const double Dist = 0.5;
  const double HadrT = 2.0;
  const double HadrZ = 10.0;

  while(mt_val<=max_mt+0.1*step_mt){

    TF1* fSrcC_pp = new TF1("fSrcC_pp",&pp_src,&DLM_CecaSource_v0::RootEval,0,8,5,"DLM_CecaSource_v0","RootEval");
    fSrcC_pp->SetParameter(0,mt_val);
    fSrcC_pp->SetParameter(1,Dist);
    fSrcC_pp->SetParameter(2,HadrT);
    fSrcC_pp->SetParameter(3,HadrZ);
    fSrcC_pp->FixParameter(4,1);
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






CecaAnalysis1::CecaAnalysis1(TString AnalysisVersion, TString SourceVersion, TString DLM_AnalysisFolder):
                            NumMtBins_pp(7),NumMtBins_pL(6),catsfiles_fld(DLM_AnalysisFolder),ana_ver(AnalysisVersion),
                            pp_cats_max(400),pp_ck_max(600),pL_cats_max(360),pL_ck_max(600),
                            SourceFitPars(5),BaselineFitPars(5),
                            ExtraFitPars(AnalysisVersion=="Simple"?0:AnalysisVersion=="Reduced"?0:0),
                            TotFitPars(SourceFitPars+BaselineFitPars+ExtraFitPars)
                            {

                              printf("\033[41mDON'T FORGET TO CHECK IF THE DIFFERENT BINNING FUCKS UP THE FEED-DOWNS!!!!!!!\033[0m\n");
  printf("Activating CecaAnalysis1:\n");
  printf("  DLM_AnalysisFolder = %s\n",catsfiles_fld.Data());
  printf("  SourceVersion = %s\n",SourceVersion.Data());
  printf("  AnalysisVersion = %s\n",AnalysisVersion.Data());
  if(AnalysisVersion!="Simple"&&AnalysisVersion!="Reduced"){
    printf("\033[1;31mERROR:\033[0m Unkown analysis version %s.\n",AnalysisVersion.Data());
    return;
  }

  TString SourceFolder_pp = catsfiles_fld+"/Source/CECA/LookUpSource/"+SourceVersion+"/pp/";
  //printf("  SourceFolder (pp)  = %s\n",SourceFolder_pp.Data());
  TString SourceFolder_pL = catsfiles_fld+"/Source/CECA/LookUpSource/"+SourceVersion+"/pL/";
  //printf("  SourceFolder (pL)  = %s\n",SourceFolder_pL.Data());

  lam = new double [15];
  settings = new double [8];

  FunctionCallCounter = 0;
  FittingScenario = 0;

  mT_BinCenter_pp = new double [NumMtBins_pp];
  mT_BinCenter_pp[0] = 1107.7;
  mT_BinCenter_pp[1] = 1168.3;
  mT_BinCenter_pp[2] = 1228.4;
  mT_BinCenter_pp[3] = 1315.6;
  mT_BinCenter_pp[4] = 1462.8;
  mT_BinCenter_pp[5] = 1687.2;
  mT_BinCenter_pp[6] = 2211.6;

  mT_BinCenter_pL = new double [NumMtBins_pL];
  mT_BinCenter_pL[0] = 1212.4;
  mT_BinCenter_pL[1] = 1289.6;
  mT_BinCenter_pL[2] = 1376.0;
  mT_BinCenter_pL[3] = 1540.7;
  mT_BinCenter_pL[4] = 1756.0;
  mT_BinCenter_pL[5] = 2259.4;

  hCkExp_pp = NULL;
  hCkExp_pL = NULL;

  fCk_pp = NULL;
  fCk_pL = NULL;

  fSrc_pp = NULL;
  fSrc_pL = NULL;


  if(ana_ver=="Simple"){
    settings[Model_pS0] = 0;
    settings[Model_pL] = 11600;

    settings[FemtoRegMax_pp] = 376;
    settings[FitMax_pp] = 376;
    //CutOff_pp = pp_cats_max;
    settings[CkConv_pp] = 700;

    settings[FemtoRegMax_pL] = 336;
    settings[FitMax_pL] = 456;
    //CutOff_pL = pL_cats_max;
    settings[CkConv_pL] = 700;
  }
  else if(ana_ver=="Reduced"){
    settings[Model_pS0] = 0;
    settings[Model_pL] = 11600;

    settings[FemtoRegMax_pp] = 184;
    settings[FitMax_pp] = 184;
    settings[CkConv_pp] = 700;

    settings[FemtoRegMax_pL] = 180;
    settings[FitMax_pL] = 180;
    settings[CkConv_pL] = 700;
  }


  Kitty_pp = NULL;
  Kitty_pL = NULL;
  Kitty_pS0 = NULL;
  Kitty_pXi0 = NULL;
  Kitty_pXim = NULL;
  Kitty_pXim1530 = NULL;

  Ck_pp = NULL;
  Ck_pL = NULL;
  Ck_pS0 = NULL;
  Ck_pXi0 = NULL;
  Ck_pXim = NULL;
  Ck_pXim1530 = NULL;

  //CkMt_pp = NULL;
  //CkMt_pL = NULL;
  //CkMt_pS0 = NULL;
  //CkMt_pXi0 = NULL;
  //CkMt_pXim = NULL;
  ////CkMt_pXim1530 = NULL;
  MaxHistoryEntries = 8196;
  HistoryEntries = 0;
  HistoryCk_pp = NULL;
  HistoryCk_pL = NULL;
  HistoryPars = NULL;
  HistoryPrecision = NULL;
  PresentDay = 0;


  CkDec_pp = NULL;
  CkDec_pL = NULL;
  CkDec_pS0 = NULL;
  CkDec_pXi0 = NULL;
  CkDec_pXim = NULL;
  CkDec_pXim1530 = NULL;

  hReso_pp = NULL;
  hPs_pp = NULL;
  hReso_pL = NULL;
  hPs_pL = NULL;

  hFeed_pp_pL = NULL;
  hFeed_pL_pS0 = NULL;
  hFeed_pL_pXim = NULL;
  hFeed_pL_pXi0 = NULL;
  hFeed_pXi_pXim1530 = NULL;

  Src_pp = NULL;//done
  Src_pL = NULL;//done

  GrandeFitter = NULL;
  AnalysisObject = NULL;//done

  printf("\n");
  printf("...Loading pp source...\n");
  Src_pp = new DLM_CecaSource_v0("pp",SourceVersion.Data(),SourceFolder_pp.Data());
  printf("--> Completed\n");
  printf("...Loading pL source...\n");
  Src_pL = new DLM_CecaSource_v0("pL",SourceVersion.Data(),SourceFolder_pL.Data());
  printf("--> Completed\n");

  if(Src_pp->InErrorState()){
    printf("\033[1;31mERROR:\033[0m CecaAnalysis1 failed to load the pp source %s from %s.\n",SourceVersion.Data(),SourceFolder_pp.Data());
    delete Src_pp;
    Src_pp = NULL;
  }
  if(Src_pL->InErrorState()){
    printf("\033[1;31mERROR:\033[0m CecaAnalysis1 failed to load the pL source %s from %s.\n",SourceVersion.Data(),SourceFolder_pL.Data());
    delete Src_pL;
    Src_pL = NULL;
  }

  AnalysisObject = new DLM_CommonAnaFunctions();
  AnalysisObject->SetCatsFilesFolder(catsfiles_fld.Data());

  Kitty_pp = new CATS();
  Kitty_pL = new CATS();
  Kitty_pS0 = new CATS();
  Kitty_pXi0 = new CATS();
  Kitty_pXim = new CATS();
  Kitty_pXim1530 = new CATS();

  //OldPar = new double[TotFitPars];
  //NewPar = new double[TotFitPars];
  PresentPars =  new double[TotFitPars];

  for(unsigned uPar=0; uPar<TotFitPars; uPar++){
    //OldPar[uPar] = 0;
    //NewPar[uPar] = 0;
  }

  LoadData();

}


CecaAnalysis1::~CecaAnalysis1(){
//printf("~CecaAnalysis1\n");
  DelHistory();
  if(mT_BinCenter_pp) {delete [] mT_BinCenter_pp; mT_BinCenter_pp = NULL;}
  if(mT_BinCenter_pL) {delete [] mT_BinCenter_pL; mT_BinCenter_pL = NULL;}
  if(lam) {delete [] lam; lam = NULL;}
  if(settings) {delete [] settings; settings = NULL;}
  //if(OldPar) {delete [] OldPar; OldPar = NULL;}
  //if(NewPar) {delete [] NewPar; NewPar = NULL;}
  if(PresentPars) {delete [] PresentPars; PresentPars = NULL;}
  if(hCkExp_pp){
    for(unsigned uMt=0; uMt<NumMtBins_pp; uMt++){
      if(hCkExp_pp[uMt]){
        delete hCkExp_pp[uMt];
        hCkExp_pp[uMt] = NULL;
      }
    }
    delete [] hCkExp_pp;
    hCkExp_pp = NULL;
  }
  if(hCkExp_pL){
    for(unsigned uMt=0; uMt<NumMtBins_pL; uMt++){
      if(hCkExp_pL[uMt]){
        delete hCkExp_pL[uMt];
        hCkExp_pL[uMt] = NULL;
      }
    }
    delete [] hCkExp_pL;
    hCkExp_pL = NULL;
  }

  if(fCk_pp){
    for(unsigned uMt=0; uMt<NumMtBins_pp; uMt++){
      if(fCk_pp[uMt]){
        delete fCk_pp[uMt];
        fCk_pp[uMt] = NULL;
      }
    }
    delete [] fCk_pp;
    fCk_pp = NULL;
  }
  if(fCk_pL){
    for(unsigned uMt=0; uMt<NumMtBins_pL; uMt++){
      if(fCk_pL[uMt]){
        delete fCk_pL[uMt];
        fCk_pL[uMt] = NULL;
      }
    }
    delete [] fCk_pL;
    fCk_pL = NULL;
  }

  if(fSrc_pp){
    for(unsigned uMt=0; uMt<NumMtBins_pp; uMt++){
      if(fSrc_pp[uMt]){
        delete fSrc_pp[uMt];
        fSrc_pp[uMt] = NULL;
      }
    }
    delete [] fSrc_pp;
    fSrc_pp = NULL;
  }
  if(fSrc_pL){
    for(unsigned uMt=0; uMt<NumMtBins_pL; uMt++){
      if(fSrc_pL[uMt]){
        delete fSrc_pL[uMt];
        fSrc_pL[uMt] = NULL;
      }
    }
    delete [] fSrc_pL;
    fSrc_pL = NULL;
  }

  if(Kitty_pp){
    delete Kitty_pp;
    Kitty_pp = NULL;
  }
  if(Kitty_pL){
    delete Kitty_pL;
    Kitty_pL = NULL;
  }
  if(Kitty_pS0){
    delete Kitty_pS0;
    Kitty_pS0 = NULL;
  }
  if(Kitty_pXi0){
    delete Kitty_pXi0;
    Kitty_pXi0 = NULL;
  }
  if(Kitty_pXim){
    delete Kitty_pXim;
    Kitty_pXim = NULL;
  }
  if(Kitty_pXim1530){
    delete Kitty_pXim1530;
    Kitty_pXim1530 = NULL;
  }
  if(Ck_pp){
    delete Ck_pp;
    Ck_pp = NULL;
  }
  if(Ck_pL){
    delete Ck_pL;
    Ck_pL = NULL;
  }
  if(Ck_pS0){
    delete Ck_pS0;
    Ck_pS0 = NULL;
  }
  if(Ck_pXi0){
    delete Ck_pXi0;
    Ck_pXi0 = NULL;
  }
  if(Ck_pXim){
    delete Ck_pXim;
    Ck_pXim = NULL;
  }
  if(Ck_pXim1530){
    delete Ck_pXim1530;
    Ck_pXim1530 = NULL;
  }
//  if(CkMt_pp){
//    for(unsigned uMt=0; uMt<NumMtBins_pp; uMt++){
//      if(CkMt_pp[uMt]){delete CkMt_pp[uMt];CkMt_pp[uMt]=NULL;}
//    }
//    delete [] CkMt_pp;
//    CkMt_pp = NULL;
//  }
//  if(CkMt_pL){
//    for(unsigned uMt=0; uMt<NumMtBins_pL; uMt++){
//      if(CkMt_pL[uMt]){delete CkMt_pL[uMt];CkMt_pL[uMt]=NULL;}
//    }
//    delete [] CkMt_pL;
//    CkMt_pL = NULL;
//  }
//  if(CkMt_pS0){
//    for(unsigned uMt=0; uMt<NumMtBins_pL; uMt++){
//      if(CkMt_pS0[uMt]){delete CkMt_pS0[uMt];CkMt_pS0[uMt]=NULL;}
//    }
//    delete [] CkMt_pS0;
//    CkMt_pS0 = NULL;
//  }
//  if(CkMt_pXi0){
//    for(unsigned uMt=0; uMt<NumMtBins_pL; uMt++){
//      if(CkMt_pXi0[uMt]){delete CkMt_pXi0[uMt];CkMt_pXi0[uMt]=NULL;}
//    }
//    delete [] CkMt_pXi0;
//    CkMt_pXi0 = NULL;
//  }
//  if(CkMt_pXim){
//    for(unsigned uMt=0; uMt<NumMtBins_pL; uMt++){
//      if(CkMt_pXim[uMt]){delete CkMt_pXim[uMt];CkMt_pXim[uMt]=NULL;}
//    }
//    delete [] CkMt_pXim;
//    CkMt_pXim = NULL;
//  }
  //if(CkMt_pXim1530){
  //  for(unsigned uMt=0; uMt<NumMtBins_pL; uMt++){
  //    if(CkMt_pXim1530[uMt]){delete CkMt_pXim1530[uMt];CkMt_pXim1530[uMt]=NULL;}
  //  }
  //  delete [] CkMt_pXim1530;
  //  CkMt_pXim1530 = NULL;
  //}
  if(CkDec_pp){
    delete CkDec_pp;
    CkDec_pp = NULL;
  }
  if(CkDec_pL){
    delete CkDec_pL;
    CkDec_pL = NULL;
  }
  if(CkDec_pS0){
    delete CkDec_pS0;
    CkDec_pS0 = NULL;
  }
  if(CkDec_pXi0){
    delete CkDec_pXi0;
    CkDec_pXi0 = NULL;
  }
  if(CkDec_pXim){
    delete CkDec_pXim;
    CkDec_pXim = NULL;
  }
  if(CkDec_pXim1530){
    delete CkDec_pXim1530;
    CkDec_pXim1530 = NULL;
  }
  if(hReso_pp){
    delete hReso_pp;
    hReso_pp = NULL;
  }
  if(hPs_pp){
    delete hPs_pp;
    hPs_pp = NULL;
  }
  if(hReso_pL){
    delete hReso_pL;
    hReso_pL = NULL;
  }
  if(hPs_pL){
    delete hPs_pL;
    hPs_pL = NULL;
  }
  if(hFeed_pp_pL){
    delete hFeed_pp_pL;
    hFeed_pp_pL = NULL;
  }
  if(hFeed_pL_pS0){
    delete hFeed_pL_pS0;
    hFeed_pL_pS0 = NULL;
  }
  if(hFeed_pL_pXim){
    delete hFeed_pL_pXim;
    hFeed_pL_pXim = NULL;
  }
  if(hFeed_pL_pXi0){
    delete hFeed_pL_pXi0;
    hFeed_pL_pXi0 = NULL;
  }
  if(hFeed_pXi_pXim1530){
    delete hFeed_pXi_pXim1530;
    hFeed_pXi_pXim1530 = NULL;
  }
  if(Src_pp){
    delete Src_pp;
    Src_pp = NULL;
  }
  if(Src_pL){
    delete Src_pL;
    Src_pL = NULL;
  }
  if(GrandeFitter){
    delete GrandeFitter;
    GrandeFitter = NULL;
  }
  if(AnalysisObject){
    delete AnalysisObject;
    AnalysisObject = NULL;
  }
}


void CecaAnalysis1::DelHistory(){
  if(HistoryCk_pp){
    for(unsigned uHst=0; uHst<MaxHistoryEntries; uHst++){
      if(HistoryCk_pp[uHst]){
        for(unsigned uMt=0; uMt<NumMtBins_pp; uMt++){
            if(HistoryCk_pp[uHst][uMt]){
              delete HistoryCk_pp[uHst][uMt];
              HistoryCk_pp[uHst][uMt] = NULL;
            }
        }
        delete [] HistoryCk_pp[uHst];
        HistoryCk_pp[uHst] = NULL;
      }
    }
    delete [] HistoryCk_pp;
    HistoryCk_pp = NULL;
  }
  if(HistoryCk_pL){
    for(unsigned uHst=0; uHst<MaxHistoryEntries; uHst++){
      if(HistoryCk_pL[uHst]){
        for(unsigned uMt=0; uMt<NumMtBins_pL; uMt++){
            if(HistoryCk_pL[uHst][uMt]){
              delete HistoryCk_pL[uHst][uMt];
              HistoryCk_pL[uHst][uMt] = NULL;
            }
        }
        delete [] HistoryCk_pL[uHst];
        HistoryCk_pL[uHst] = NULL;
      }
    }
    delete [] HistoryCk_pL;
    HistoryCk_pL = NULL;
  }
  if(HistoryPars){
    for(unsigned uHst=0; uHst<MaxHistoryEntries; uHst++){
      if(HistoryPars[uHst]){
        delete HistoryPars[uHst];
        HistoryPars[uHst] = NULL;
      }
    }
    delete [] HistoryPars;
    HistoryPars = NULL;
  }
  if(HistoryPrecision){
    delete [] HistoryPrecision;
    HistoryPrecision = NULL;
  }
}

void CecaAnalysis1::InitHistory(){
//printf("hey\n");
  DelHistory();
  PresentDay = 0;
//printf("deled\n");
  HistoryCk_pp = new DLM_Ck** [MaxHistoryEntries];
  for(unsigned uHst=0; uHst<MaxHistoryEntries; uHst++){
    HistoryCk_pp[uHst] = new DLM_Ck* [NumMtBins_pp];
    for(unsigned uMt=0; uMt<NumMtBins_pp; uMt++){
      HistoryCk_pp[uHst][uMt] = new DLM_Ck(*Ck_pp,false);
    }
  }
//printf("1\n");
  HistoryCk_pL = new DLM_Ck** [MaxHistoryEntries];
  for(unsigned uHst=0; uHst<MaxHistoryEntries; uHst++){
    HistoryCk_pL[uHst] = new DLM_Ck* [NumMtBins_pL];
    for(unsigned uMt=0; uMt<NumMtBins_pL; uMt++){
      HistoryCk_pL[uHst][uMt] = new DLM_Ck(*Ck_pL,false);
    }
  }
//printf("2\n");
  HistoryPars = new float* [MaxHistoryEntries];
  for(unsigned uHst=0; uHst<MaxHistoryEntries; uHst++){
    HistoryPars[uHst] = new float [TotFitPars];
  }
//printf("3\n");
  HistoryPrecision = new float [TotFitPars];
  for(unsigned uPar=0; uPar<TotFitPars; uPar++){
    HistoryPrecision[uPar] = 0;
  }
  //set up the defaults for the ceca source
  HistoryPrecision[1] = 0.0005;
  HistoryPrecision[2] = 0.005;
  HistoryPrecision[3] = 0.005;
//printf("done\n");
}

void CecaAnalysis1::SetNumHistEntries(unsigned numetries){
  MaxHistoryEntries = numetries;
}

void CecaAnalysis1::LoadData(){
  printf("... Loading the data...\n");
  if(hCkExp_pp || hCkExp_pL){
    printf("\033[1;31mERROR:\033[0m CecaAnalysis1::LoadDatas can load the data only one time!\n");
    return;
  }

  hCkExp_pp = new TH1F* [NumMtBins_pp];
  hCkExp_pL = new TH1F* [NumMtBins_pL];

  for(unsigned uMt=0; uMt<NumMtBins_pp; uMt++){
    hCkExp_pp[uMt] = AnalysisObject->GetAliceExpCorrFun("pp13TeV_HM_BernieSource","pp","0",0,1,uMt);
    hCkExp_pp[uMt]->SetName(TString::Format("hCkExp_pp[%u]",uMt));
  }
  for(unsigned uMt=0; uMt<NumMtBins_pL; uMt++){
    hCkExp_pL[uMt] = AnalysisObject->GetAliceExpCorrFun("pp13TeV_HM_BernieSource","pLambda","0",0,1,uMt);
    hCkExp_pL[uMt]->SetName(TString::Format("hCkExp_pL[%u]",uMt));
  }
  printf("---> Completed\n");
}


DLM_Ck* CecaAnalysis1::GetCkHistory(TString system, unsigned WhichMt){
  if(!CkDec_pp || !CkDec_pL || !CkDec_pS0 || !CkDec_pXim || !CkDec_pXi0 || !CkDec_pXim1530){
    return NULL;
  }

//  FunctionCallCounter++;
//  if(FunctionCallCounter%1000==0){
//    //printf("\r\033[K GetCkHistory has been called %u times",FunctionCallCounter);
//    //cout << flush;
//    printf("GetCkHistory has been called %uk times\n",FunctionCallCounter/1000);
//  }
//printf("GetCkHistory has been called %u times\n",FunctionCallCounter);
//printf("  system = %s\n",system.Data());
//printf("  present pars: %.3e %.3e %.3e %.3e %.3e\n",PresentPars[0],PresentPars[1],PresentPars[2],PresentPars[3],PresentPars[4]);
//printf("  WhichMt=%u\n",WhichMt);
//printf("  PresentDay=%u\n",PresentDay);
//printf("  HistoryEntries=%u\n",HistoryEntries);

  for(int iHst=0; iHst<HistoryEntries; iHst++){
    //unsigned uPresent = (PresentDay+uHst)%(MaxHistoryEntries);
    //we iterate backwards, as we hope to find the best solution from the last few iterations
    int iHistory = (int(PresentDay)-iHst-1);
    if(iHistory<0) iHistory = int(HistoryEntries)+iHistory;

    bool WeHaveIt = true;
    for(unsigned uPar=0; uPar<TotFitPars; uPar++){
      //we only check the non-baseline parameters
      if(uPar>=SourceFitPars && uPar<SourceFitPars+BaselineFitPars) continue;
      //we dont care about mT, as we fill all mT bins any times for a given set of pars
      if(uPar==0) continue;
      //we also do not care about the scaling parameter, as we will not change it for pp/pL
      if(uPar==4) continue;
      //a Gaussian source has only uPar==1
      if(PresentPars[0]==-1 && uPar==2|| uPar==3) continue;
      WeHaveIt *= (fabs(HistoryPars[iHistory][uPar]-PresentPars[uPar])<=HistoryPrecision[uPar]);
      //printf("    cmp: %.3e vs %.3e: %.3e < %.3e\n",HistoryPars[iHistory][uPar],PresentPars[uPar],
      //            fabs(HistoryPars[iHistory][uPar]-PresentPars[uPar]), HistoryPrecision[uPar]);
    }
//    printf("    Compare with %i: %.3e %.3e %.3e %.3e %.3e\n",iHistory,
//      HistoryPars[iHistory][0],HistoryPars[iHistory][1],HistoryPars[iHistory][2],
//      HistoryPars[iHistory][3],HistoryPars[iHistory][4]);
    if(WeHaveIt){
//printf("  -> we have it!\n");
//usleep(1000e3);
      if(system=="pp"){
        return HistoryCk_pp[iHistory][WhichMt];
      }
      else if(system=="pL"){
        //printf("iHistory = %i (%u)\n",iHistory,WhichMt);
        return HistoryCk_pL[iHistory][WhichMt];
      }
      break;//the iHst
    }
  }//iHst
//printf("  -> we DONT have it\n");
//usleep(1000e3);

  FunctionCallCounter++;
  if(FunctionCallCounter%100==0){
    //printf("\r\033[K GetCkHistory has been executed %u times",FunctionCallCounter);
    //cout << flush;
    //printf("GetCkHistory has been executed %u times\n",FunctionCallCounter);
    //printf("  system = %s\n",system.Data());
    //printf("  present pars: %.0f %.5e %.4e %.4e %.1f\n",PresentPars[0],PresentPars[1],PresentPars[2],PresentPars[3],PresentPars[4]);
    //printf("  WhichMt=%u\n",WhichMt);
    //usleep(500e3);
    //printf("  PresentDay=%u\n",PresentDay);
    //printf("  HistoryEntries=%u\n",HistoryEntries);
  }


  //if do not have a suitable entry:
  //evaluate it and plug it into the history
  if(CkDec_pp->GetCk()->GetNumSourcePar()>SourceFitPars){
    printf("\033[1;31mERROR:\033[0m CecaAnalysis1::UpdateDecomp got %u parameters !?!?\n",CkDec_pp->GetCk()->GetNumSourcePar());
    return NULL;
  }
  if(CkDec_pL->GetCk()->GetNumSourcePar()!=CkDec_pp->GetCk()->GetNumSourcePar()){
    printf("\033[1;31mERROR:\033[0m CecaAnalysis1::UpdateDecomp got non-identical sources for pp and pL\n");
    return NULL;
  }



//printf("PresentDay=%u\n",PresentDay);
  //for pp ONLY (later pL)
  for(unsigned uMt=0; uMt<NumMtBins_pp; uMt++){
    double MtVal;
    if(PresentPars[0]==-1) MtVal = -1;
    else MtVal = GetAnalysisMt_pp(uMt);
//printf("MtVal = %.0f\n",MtVal);
    for(unsigned uPar=0; uPar<SourceFitPars; uPar++){
      //this is the scaling parameter, this as assumed to be fixed for each contribution and is
      //NOT a common parameter, hence we do not set it here (neither have an option to fit it!)
      if(uPar==4) continue;
      //mT
      if(uPar==0){
//DEB Check your notes what and how this should be set up
        CkDec_pp->GetCk()->SetSourcePar(uPar,MtVal);
        CkDec_pL->GetCk()->SetSourcePar(uPar,MtVal);
        CkDec_pS0->GetCk()->SetSourcePar(uPar,MtVal);
        CkDec_pXim->GetCk()->SetSourcePar(uPar,MtVal);
        CkDec_pXi0->GetCk()->SetSourcePar(uPar,MtVal);
      }
      else{
        //for a Gaussian source, set only
        //if(uPar==0 && PresentPars[uPar]==-1){
        //  continue;
        //}
        CkDec_pp->GetCk()->SetSourcePar(uPar,PresentPars[uPar]);
        CkDec_pL->GetCk()->SetSourcePar(uPar,PresentPars[uPar]);
        CkDec_pS0->GetCk()->SetSourcePar(uPar,PresentPars[uPar]);
        CkDec_pXim->GetCk()->SetSourcePar(uPar,PresentPars[uPar]);
        CkDec_pXi0->GetCk()->SetSourcePar(uPar,PresentPars[uPar]);
      }
    }//uPar
    //we only update ones if we have a Gaussian source (same for all mt bins)
    if(PresentPars[0]!=-1 || uMt==0){
      CkDec_pXim->Update();
      CkDec_pXi0->Update();
      CkDec_pS0->Update();
      CkDec_pL->Update();
      CkDec_pp->Update();
    }


    for(unsigned uMom=0; uMom<HistoryCk_pp[PresentDay][uMt]->GetNbins(); uMom++){
      double MOM = HistoryCk_pp[PresentDay][uMt]->GetBinCenter(0,uMom);
      HistoryCk_pp[PresentDay][uMt]->SetBinContent(uMom,CkDec_pp->EvalCk(MOM));
    }
  }//uMt pp


  //for pL ONLY
  for(unsigned uMt=0; uMt<NumMtBins_pL; uMt++){
    double MtVal;
    if(PresentPars[0]==-1) MtVal = -1;
    else MtVal = GetAnalysisMt_pL(uMt);
    for(unsigned uPar=0; uPar<SourceFitPars; uPar++){
      //this is the scaling parameter, this as assumed to be fixed for each contribution and is
      //NOT a common parameter, hence we do not set it here (neither have an option to fit it!)
      if(uPar==4) continue;
      //mT
      if(uPar==0){
//DEB Check your notes what and how this should be set up
        //CkDec_pp->GetCk()->SetSourcePar(uPar,MtVal);
        CkDec_pL->GetCk()->SetSourcePar(uPar,MtVal);
        CkDec_pS0->GetCk()->SetSourcePar(uPar,MtVal);
        CkDec_pXim->GetCk()->SetSourcePar(uPar,MtVal);
        CkDec_pXi0->GetCk()->SetSourcePar(uPar,MtVal);
      }
      else{
        //CkDec_pp->GetCk()->SetSourcePar(uPar,PresentPars[uPar]);
        CkDec_pL->GetCk()->SetSourcePar(uPar,PresentPars[uPar]);
        CkDec_pS0->GetCk()->SetSourcePar(uPar,PresentPars[uPar]);
        CkDec_pXim->GetCk()->SetSourcePar(uPar,PresentPars[uPar]);
        CkDec_pXi0->GetCk()->SetSourcePar(uPar,PresentPars[uPar]);
      }
    }//uPar
    //we only update ones if we have a Gaussian source (same for all mt bins)
    if(PresentPars[0]!=-1 || uMt==0){
      CkDec_pXim->Update();
      CkDec_pXi0->Update();
      CkDec_pS0->Update();
      CkDec_pL->Update();
      //CkDec_pp->Update();
    }

    for(unsigned uMom=0; uMom<HistoryCk_pL[PresentDay][uMt]->GetNbins(); uMom++){
      double MOM = HistoryCk_pL[PresentDay][uMt]->GetBinCenter(0,uMom);
      HistoryCk_pL[PresentDay][uMt]->SetBinContent(uMom,CkDec_pL->EvalCk(MOM));
    }
  }//uMt pp


  for(unsigned uPar=0; uPar<TotFitPars; uPar++){
    HistoryPars[PresentDay][uPar] = PresentPars[uPar];
  }

  unsigned PRESENT = PresentDay;
  //we increase the number of entries in history
  if(HistoryEntries<MaxHistoryEntries) HistoryEntries++;
  PresentDay = (PresentDay+1)%(MaxHistoryEntries);

  if(system=="pp"){
    return HistoryCk_pp[PRESENT][WhichMt];
  }
  else{
    return HistoryCk_pL[PRESENT][WhichMt];
  }
}

unsigned CecaAnalysis1::GetAnalysisMtBin_pp(double Mt){
       if(Mt<1150){return 0;}
  else if(Mt<1200){return 1;}
  else if(Mt<1300){return 2;}
  else if(Mt<1400){return 3;}
  else if(Mt<1500){return 4;}
  else if(Mt<1800){return 5;}
  else {return 6;}
  return 0;
}
double CecaAnalysis1::GetAnalysisMt_pp(unsigned WhichMt){
  if(WhichMt>=NumMtBins_pp) return 0;
  return mT_BinCenter_pp[WhichMt];
}

unsigned CecaAnalysis1::GetAnalysisMtBin_pL(double Mt){
       if(Mt<1250){return 0;}
  else if(Mt<1300){return 1;}
  else if(Mt<1400){return 2;}
  else if(Mt<1600){return 3;}
  else if(Mt<1800){return 4;}
  else {return 5;}
  return 0;
}
double CecaAnalysis1::GetAnalysisMt_pL(unsigned WhichMt){
  if(WhichMt>=NumMtBins_pL) return 0;
  return mT_BinCenter_pL[WhichMt];
}

double CecaAnalysis1::FitFun_pp(double* xval, double* par){
  double Baseline = DLM_Baseline(xval,&par[5]);
  double Femto=1;
  double kstar = *xval;
  unsigned WhichMt = GetAnalysisMtBin_pp(par[0]);

  for(unsigned uPar=0; uPar<TotFitPars; uPar++){
    PresentPars[uPar] = par[uPar];
  }

  DLM_Ck* CkSaved = GetCkHistory("pp",WhichMt);
  Femto = CkSaved->Eval(kstar);

  return Baseline*Femto;
}

//N.B. we will reject the data points around the cusp, not to overcomplicate things
double CecaAnalysis1::FitFun_pL(double* xval, double* par){
  double Baseline = DLM_Baseline(xval,&par[5]);
  double Femto=1;
  double kstar = *xval;
  unsigned WhichMt = GetAnalysisMtBin_pL(par[0]);

  if(kstar>260&&kstar<300){
    TF1::RejectPoint();
  }

  for(unsigned uPar=0; uPar<TotFitPars; uPar++){
    PresentPars[uPar] = par[uPar];
  }

  DLM_Ck* CkSaved = GetCkHistory("pL",WhichMt);
  Femto = CkSaved->Eval(kstar);

  return Baseline*Femto;
}

void CecaAnalysis1::SetUp_pp(const TString &POT, const int &PotVar){
  printf("...Setting up pp...\n");
  printf("  POT = %s (%i)\n",POT.Data(),PotVar);

  if(!AnalysisObject){
    printf("\033[1;33mWARNING!\033[0m CecaAnalysis1::SetUp_pp needs an AnalysisObject. BUG??\n");
    return;
  }
  if(Ck_pp||CkDec_pp){
    printf("\033[1;31mERROR:\033[0m At the moment SetUp_pp can be called only ones!");
    return;
  }

  Kitty_pp->SetMomBins(pp_cats_max/4,0,pp_cats_max);
  //we set up without a source, set it up below
  AnalysisObject->SetUpCats_pp(*Kitty_pp,POT,"",PotVar,0);
  Kitty_pp->SetUseAnalyticSource(true);
  Kitty_pp->SetAnaSource(CatsSourceForwarder, Src_pp, 5);
  //this is a Gaussian by default
  Kitty_pp->SetAnaSource(0, -1);
  Kitty_pp->SetAnaSource(1, 1);
  Kitty_pp->SetAnaSource(2, 0);
  Kitty_pp->SetAnaSource(3, 0);
  Kitty_pp->SetAnaSource(4, 1);

  Ck_pp = new DLM_Ck(Kitty_pp->GetNumSourcePars(),0,*Kitty_pp,pp_ck_max/4,0,pp_ck_max);
  for(unsigned uSP=0; uSP<Kitty_pp->GetNumSourcePars(); uSP++){
    Ck_pp->SetSourcePar(uSP,Kitty_pp->GetAnaSourcePar(uSP));
  }
  Ck_pp->SetCutOff(pp_cats_max,CkConv_pp);

  printf("--> Completed\n");
}
void CecaAnalysis1::SetUp_pL(const TString &POT, const int &PotVar){
  printf("...Setting up pL...\n");
  printf("  POT = %s (%i)\n",POT.Data(),PotVar);

  if(!AnalysisObject){
    printf("\033[1;33mWARNING!\033[0m CecaAnalysis1::SetUp_pL needs an AnalysisObject. BUG??\n");
    return;
  }
  if(Ck_pL||CkDec_pL){
    printf("\033[1;31mERROR:\033[0m At the moment SetUp_pL can be called only ones!");
    return;
  }

  Kitty_pL->SetMomBins(pL_cats_max/12,0,pL_cats_max);
  AnalysisObject->SetUpCats_pL(*Kitty_pL,POT,"",PotVar,0);
  Kitty_pL->SetUseAnalyticSource(true);
  Kitty_pL->SetAnaSource(CatsSourceForwarder, Src_pL, 5);
  Kitty_pL->SetAnaSource(0, -1);
  Kitty_pL->SetAnaSource(1, 1);
  Kitty_pL->SetAnaSource(2, 0);
  Kitty_pL->SetAnaSource(3, 0);
  Kitty_pL->SetAnaSource(4, 1);

  if(Kitty_pL->GetNumChannels()>7){
    const double CUSP_WEIGHT = 0.33;
    Kitty_pL->SetChannelWeight(7,1./4.*CUSP_WEIGHT);//1S0 SN(s) -> LN(s)
    Kitty_pL->SetChannelWeight(8,3./4.*CUSP_WEIGHT);//3S1 SN(s) -> LN(s)
    Kitty_pL->SetChannelWeight(10,3./4.*CUSP_WEIGHT);//3S1 SN(d) -> LN(s)
    Kitty_pL->SetChannelWeight(13,3./20.*CUSP_WEIGHT);//3D1 SN(d) -> LN(d)
    Kitty_pL->SetChannelWeight(15,3./20.*CUSP_WEIGHT);//3D1 SN(s) -> LN(d)
  }


  Ck_pL = new DLM_Ck(Kitty_pL->GetNumSourcePars(),0,*Kitty_pL,pL_ck_max/12,0,pL_ck_max);
  for(unsigned uSP=0; uSP<Kitty_pL->GetNumSourcePars(); uSP++){
    Ck_pL->SetSourcePar(uSP,Kitty_pL->GetAnaSourcePar(uSP));
  }
  Ck_pL->SetCutOff(pL_cats_max,CkConv_pL);

  printf("--> Completed\n");
}
void CecaAnalysis1::SetUp_pS0(const TString &POT, const int &PotVar){
  printf("...Setting up pS0...\n");
  printf("  POT = %s (%i)\n",POT.Data(),PotVar);

  if(!AnalysisObject){
    printf("\033[1;33mWARNING!\033[0m CecaAnalysis1::SetUp_pS0 needs an AnalysisObject. BUG??\n");
    return;
  }
  if(Ck_pS0||CkDec_pS0){
    printf("\033[1;31mERROR:\033[0m At the moment SetUp_pS0 can be called only ones!");
    return;
  }

  Kitty_pS0->SetMomBins(pL_cats_max/12,0,pL_cats_max);
  AnalysisObject->SetUpCats_pS0(*Kitty_pS0,POT,"",PotVar,0);
  Kitty_pS0->SetUseAnalyticSource(true);
  //N.B. we use the proton-proton source as it is
  Kitty_pS0->SetAnaSource(CatsSourceForwarder, Src_pp, 5);
  Kitty_pS0->SetAnaSource(0, -1);
  Kitty_pS0->SetAnaSource(1, 1);
  Kitty_pS0->SetAnaSource(2, 0);
  Kitty_pS0->SetAnaSource(3, 0);
  Kitty_pS0->SetAnaSource(4, 1);

  Ck_pS0 = new DLM_Ck(Kitty_pS0->GetNumSourcePars(),0,*Kitty_pS0,pL_ck_max/12,0,pL_ck_max);
  for(unsigned uSP=0; uSP<Kitty_pS0->GetNumSourcePars(); uSP++){
    Ck_pS0->SetSourcePar(uSP,Kitty_pS0->GetAnaSourcePar(uSP));
  }
  Ck_pS0->SetCutOff(pL_cats_max,CkConv_pL);



  printf("--> Completed\n");
}
void CecaAnalysis1::SetUp_pXim(const TString &POT, const int &PotVar, const double src_scale){
  printf("...Setting up pXim...\n");
  printf("  POT = %s (%i)\n",POT.Data(),PotVar);

  if(!AnalysisObject){
    printf("\033[1;33mWARNING!\033[0m CecaAnalysis1::SetUp_pXim needs an AnalysisObject. BUG??\n");
    return;
  }
  if(Ck_pXim||CkDec_pXim||Ck_pXim1530||CkDec_pXim1530){
    printf("\033[1;31mERROR:\033[0m At the moment SetUp_pXim can be called only ones!");
    return;
  }

  Kitty_pXim->SetMomBins(pL_cats_max/12,0,pL_cats_max);
  AnalysisObject->SetUpCats_pXim(*Kitty_pXim,POT,"",PotVar,0);
  Kitty_pXim->SetUseAnalyticSource(true);
  //N.B. we use the proton-proton source scaled by 1.06 (by def, one can change it)
  Kitty_pXim->SetAnaSource(CatsSourceForwarder, Src_pp, 5);
  Kitty_pXim->SetAnaSource(0, -1);
  Kitty_pXim->SetAnaSource(1, 1);
  Kitty_pXim->SetAnaSource(2, 0);
  Kitty_pXim->SetAnaSource(3, 0);
  Kitty_pXim->SetAnaSource(4, src_scale);

  Ck_pXim = new DLM_Ck(Kitty_pXim->GetNumSourcePars(),0,*Kitty_pXim,pL_ck_max/12,0,pL_ck_max);
  for(unsigned uSP=0; uSP<Kitty_pXim->GetNumSourcePars(); uSP++){
    Ck_pXim->SetSourcePar(uSP,Kitty_pXim->GetAnaSourcePar(uSP));
  }
  Ck_pXim->SetCutOff(pL_cats_max,CkConv_pL);

  Kitty_pXim1530->SetMomBins(pL_cats_max/12,0,pL_cats_max);
  AnalysisObject->SetUpCats_pXim(*Kitty_pXim1530,"pXim1530","Gauss",0,0);
  Kitty_pXim1530->SetAnaSource(0, 1);

  Ck_pXim1530 = new DLM_Ck(Kitty_pXim1530->GetNumSourcePars(),0,*Kitty_pXim1530,pL_ck_max/12,0,pL_ck_max);
  for(unsigned uSP=0; uSP<Kitty_pXim1530->GetNumSourcePars(); uSP++){
    Ck_pXim1530->SetSourcePar(uSP,Kitty_pXim1530->GetAnaSourcePar(uSP));
  }
  Ck_pXim1530->SetCutOff(pL_cats_max,CkConv_pL);

  printf("--> Completed\n");
}
void CecaAnalysis1::SetUp_pXi0(const TString &POT, const int &PotVar, const double src_scale){
  printf("...Setting up pXi0...\n");
  printf("  POT = %s (%i)\n",POT.Data(),PotVar);

  if(!AnalysisObject){
    printf("\033[1;33mWARNING!\033[0m CecaAnalysis1::SetUp_pXi0 needs an AnalysisObject. BUG??\n");
    return;
  }
  if(Ck_pXi0||CkDec_pXi0){
    printf("\033[1;31mERROR:\033[0m At the moment SetUp_pXi0 can be called only ones!");
    return;
  }

  Kitty_pXi0->SetMomBins(pL_cats_max/12,0,pL_cats_max);
  AnalysisObject->SetUpCats_pXi0(*Kitty_pXi0,POT,"",PotVar,0);
  Kitty_pXi0->SetUseAnalyticSource(true);
  //N.B. we use the proton-proton source scaled by 1.06 (by def, one can change it)
  Kitty_pXi0->SetAnaSource(CatsSourceForwarder, Src_pp, 5);
  Kitty_pXi0->SetAnaSource(0, -1);
  Kitty_pXi0->SetAnaSource(1, 0);
  Kitty_pXi0->SetAnaSource(2, 0);
  Kitty_pXi0->SetAnaSource(3, 0);
  Kitty_pXi0->SetAnaSource(4, src_scale);

  Ck_pXi0 = new DLM_Ck(Kitty_pXi0->GetNumSourcePars(),0,*Kitty_pXi0,pL_ck_max/12,0,pL_ck_max);
  for(unsigned uSP=0; uSP<Kitty_pXi0->GetNumSourcePars(); uSP++){
    Ck_pXi0->SetSourcePar(uSP,Kitty_pXi0->GetAnaSourcePar(uSP));
  }
  Ck_pXi0->SetCutOff(pL_cats_max,CkConv_pL);

  printf("--> Completed\n");
}


void CecaAnalysis1::SetUp_Decomposition(const int &Variation_p, const int &Variation_L){
  printf("...Setting up the decomposition...\n");
  //set up
  if(!CkDec_pp && !CkDec_pL && !CkDec_pS0 && !CkDec_pXim && !CkDec_pXi0 && !CkDec_pXim1530){

    //CREATE ALL OBJECTS  /////////////////////////////////////////////////////////
    printf("   Create all objects...\n");
    //N.B. for whatever stupid reason, pp13TeV_HM_DimiJun20 has to be used for the smearing matrix (high statistics etc.)
    //although for lambda pars it has to be the older pp13TeV_HM_Dec19, as the new one assumes 100% LambdaPurity (SB corrected)
    hReso_pp = AnalysisObject->GetResolutionMatrix("pp13TeV_HM_DimiJun20","pp");
    CkDec_pp = new DLM_CkDecomposition("pp",3,*Ck_pp,hReso_pp);

    hReso_pL = AnalysisObject->GetResolutionMatrix("pp13TeV_HM_DimiJun20","pLambda");
    CkDec_pL = new DLM_CkDecomposition("pL",5,*Ck_pL,hReso_pL);

    CkDec_pS0 = new DLM_CkDecomposition("pS0",0,*Ck_pS0,NULL);

    CkDec_pXim = new DLM_CkDecomposition("pXim",2,*Ck_pXim,NULL);
    CkDec_pXi0 = new DLM_CkDecomposition("pXi0",2,*Ck_pXi0,NULL);

    CkDec_pXim1530 = new DLM_CkDecomposition("Xim1530",0,*Ck_pXim1530,NULL);

    //GET ALL FEED-DOWN MATRICES
    printf("   Get the feed-down matrices...\n");
    hFeed_pp_pL = AnalysisObject->GetResidualMatrix("pp","pLambda");
    hFeed_pL_pS0 = AnalysisObject->GetResidualMatrix("pLambda","pSigma0");
    hFeed_pL_pXim = AnalysisObject->GetResidualMatrix("pLambda","pXim");
    hFeed_pL_pXi0 = AnalysisObject->GetResidualMatrix("pLambda","pXi0");
    hFeed_pXi_pXim1530 = AnalysisObject->GetResidualMatrix("pXim","pXim1530");
     //////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    //Get the mixed events for the phase space factor of the smearing /////////////////////////////////////////////////////////
    TList* list1_tmp;
    TList* list2_tmp;
    TFile* inFile;
    TH1F* hME_PP;
    TH1F* hME_APAP;
    TString FileName;
    //pp
    printf("   Get the pp phase space...\n");
    for(unsigned uMt=0; uMt<NumMtBins_pp; uMt++){
      FileName = TString::Format("%s/ExpData/Bernie_Source/ppData/mTBin_%i/CFOutput_mT_ppVar0_HM_%i.root",
                                          catsfiles_fld.Data(),uMt+1,uMt);
      inFile = new TFile(FileName,"read");
      //PARTICLES
      list1_tmp = (TList*)inFile->Get("PairDist");
      list2_tmp = (TList*)list1_tmp->FindObject("PairFixShifted");//there is also PairShifted ?? which one ??
      hME_PP = (TH1F*)list2_tmp->FindObject(TString::Format("MEPart_mT_%i_FixShifted",uMt));
      list1_tmp = (TList*)inFile->Get("AntiPairDist");
      list2_tmp = (TList*)list1_tmp->FindObject("PairFixShifted");
      hME_APAP = (TH1F*)list2_tmp->FindObject(TString::Format("MEAntiPart_mT_%i_FixShifted",uMt));
      gROOT->cd();
      hPs_pp = (TH1F*)hME_PP->Clone("hPs_pp");
      hPs_pp->Add(hME_APAP);
      hPs_pp->GetXaxis()->SetLimits(hPs_pp->GetXaxis()->GetXmin()*1000.,hPs_pp->GetXaxis()->GetXmax()*1000.);
      delete inFile;
    }

    //pLambda
    printf("   Get the pL phase space...\n");
    for(unsigned uMt=0; uMt<NumMtBins_pL; uMt++){
      FileName = TString::Format("%s/ExpData/Bernie_Source/pLData/mTBin_%i/CFOutput_mT_pLVar0_HM_%i.root",
                                          catsfiles_fld.Data(),uMt+1,uMt);
      inFile = new TFile(FileName,"read");
      //PARTICLES
      list1_tmp = (TList*)inFile->Get("PairDist");
      list2_tmp = (TList*)list1_tmp->FindObject("PairFixShifted");//there is also PairShifted ?? which one ??
      hME_PP = (TH1F*)list2_tmp->FindObject(TString::Format("MEPart_mT_%i_FixShifted",uMt));
      list1_tmp = (TList*)inFile->Get("AntiPairDist");
      list2_tmp = (TList*)list1_tmp->FindObject("PairFixShifted");
      hME_APAP = (TH1F*)list2_tmp->FindObject(TString::Format("MEAntiPart_mT_%i_FixShifted",uMt));
      gROOT->cd();
      hPs_pL = (TH1F*)hME_PP->Clone("hPs_pL");
      hPs_pL->Add(hME_APAP);
      hPs_pL->GetXaxis()->SetLimits(hPs_pL->GetXaxis()->GetXmin()*1000.,hPs_pL->GetXaxis()->GetXmax()*1000.);
      delete inFile;
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //LAMBDA PARAMETERS ////////////////////////////////////////////////////////////////////////
    printf("   Set up the pp lambda parameters...\n");
    double lambda_pars[5];
    double lambda_sum;

    AnalysisObject->SetUpLambdaPars_pp("pp13TeV_HM_Dec19",Variation_p,lambda_pars);
    lam[pp_gen] = lambda_pars[0];
    lam[pp_pL] = lambda_pars[1];
    lam[pp_flt] = lambda_pars[2];
    lam[pp_fake] = lambda_pars[3];
    lambda_sum = lam[pp_gen]+lam[pp_pL]+lam[pp_flt]+lam[pp_fake];
    if(fabs(1.-lambda_sum)>1e-6){
      printf("\033[1;33mWARNING!\033[0m sum(lambda_pp) = %f\n",lambda_sum);
    }

    printf("   Set up the pL lambda parameters...\n");
    AnalysisObject->SetUpLambdaPars_pL("pp13TeV_HM_Dec19",Variation_p,Variation_L,lambda_pars);
    lam[pL_gen] = lambda_pars[0];
    lam[pL_pS0] = lambda_pars[1];
    lam[pL_pXim] = lambda_pars[2];
    lam[pL_pXi0] = lambda_pars[2];
    lam[pL_flt] = lambda_pars[3]-lambda_pars[2];
    lam[pL_fake] = lambda_pars[4];
    //printf("lam_pL_gen = %.2f\n",lam_pL_gen*100.);
    //printf("lam_pL_pS0 = %.2f\n",lam_pL_pS0*100.);
    //printf("lam_pL_pXim = %.2f\n",lam_pL_pXim*100.);
    //printf("lam_pL_pXi0 = %.2f\n",lam_pL_pXi0*100.);
    //printf("lam_pL_flt = %.2f\n",lam_pL_flt*100.);
    //printf("lam_pL_fake = %.2f\n",lam_pL_fake*100.);

    lambda_sum = lam[pL_gen]+lam[pL_pS0]+lam[pL_pXim]+lam[pL_pXi0]+lam[pL_flt]+lam[pL_fake];
    if(fabs(1.-lambda_sum)>1e-6){
      printf("\033[1;33mWARNING!\033[0m sum(lambda_pL) = %f\n",lambda_sum);
    }

    printf("   Set up the pXi lambda parameters...\n");
    AnalysisObject->SetUpLambdaPars_pXim("pp13TeV_HM_Dec19",Variation_p,0,lambda_pars);
    //here we normalize to the purity, as we do not analyse experimental data and we could not care less
    //the lambda pars for Xim and Xi0 are identical
    lam[pXi_gen] = lambda_pars[0]/(1.-lambda_pars[4]);
    //the charged case
    lam[pXi_pXi1530] = lambda_pars[1]/(1.-lambda_pars[4]);
    //flat plus neutral 1530
    lam[pXi_flt] = (lambda_pars[2]+lambda_pars[3])/(1.-lambda_pars[4]);
    lambda_sum = lam[pXi_gen]+lam[pXi_pXi1530]+lam[pXi_flt];
    if(fabs(1.-lambda_sum)>1e-6){
      printf("\033[1;33mWARNING!\033[0m sum(lambda_pXi) = %f\n",lambda_sum);
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //LINK ALL OBJECTS + PHASE SPACE /////////////////////////////////////////////////////////

    printf("   Link the pp contributions...\n");

    CkDec_pp->AddContribution(0,lam[pp_pL],DLM_CkDecomposition::cFeedDown,CkDec_pL,hFeed_pp_pL);
//CkDec_pp->AddContribution(0,lam[pp_flt],DLM_CkDecomposition::cFeedDown);
    CkDec_pp->AddContribution(1,lam[pp_flt],DLM_CkDecomposition::cFeedDown);
    CkDec_pp->AddContribution(2,lam[pp_fake],DLM_CkDecomposition::cFake);
    CkDec_pp->AddPhaseSpace(hPs_pp);
    CkDec_pp->AddPhaseSpace(0,hPs_pp);

    printf("   Link the pL contributions...\n");
    CkDec_pL->AddContribution(0,lam[pL_pS0],DLM_CkDecomposition::cFeedDown,CkDec_pS0,hFeed_pL_pS0);
    CkDec_pL->AddContribution(1,lam[pL_pXim],DLM_CkDecomposition::cFeedDown,CkDec_pXim,hFeed_pL_pXim);
    CkDec_pL->AddContribution(2,lam[pL_pXi0],DLM_CkDecomposition::cFeedDown,CkDec_pXi0,hFeed_pL_pXi0);
    CkDec_pL->AddContribution(3,lam[pL_flt],DLM_CkDecomposition::cFeedDown);
    CkDec_pL->AddContribution(4,lam[pL_fake],DLM_CkDecomposition::cFake);
    CkDec_pL->AddPhaseSpace(hPs_pL);
    CkDec_pL->AddPhaseSpace(0,hPs_pL);
    CkDec_pL->AddPhaseSpace(1,hPs_pL);
    CkDec_pL->AddPhaseSpace(2,hPs_pL);

    printf("   Link the pXi contributions...\n");
    CkDec_pXim->AddContribution(0,lam[pXi_pXi1530],DLM_CkDecomposition::cFeedDown,CkDec_pXim1530,hFeed_pXi_pXim1530);
    CkDec_pXim->AddContribution(1,lam[pXi_flt],DLM_CkDecomposition::cFeedDown);
    CkDec_pXim->AddPhaseSpace(0,hPs_pL);

    CkDec_pXi0->AddContribution(0,lam[pXi_pXi1530],DLM_CkDecomposition::cFeedDown,CkDec_pXim1530,hFeed_pXi_pXim1530);
    CkDec_pXi0->AddContribution(1,lam[pXi_flt],DLM_CkDecomposition::cFeedDown);
    CkDec_pXi0->AddPhaseSpace(0,hPs_pL);

  }
  else if(!CkDec_pp || !CkDec_pL || !CkDec_pS0 || !CkDec_pXim || !CkDec_pXi0 || !CkDec_pXim1530){
    printf("\033[1;31mERROR:\033[0m SetUp_Decomposition has an impossible set of pointers!!! BUG!?\n");
    return;
  }
  printf("---> Completed\n");
}

TF1* CecaAnalysis1::GetFit(TString system, unsigned uMt){
  if(system=="pp"){
    if(uMt>=7) return NULL;
    return fCk_pp[uMt];
  }
  else if(system=="pL"||system=="pLambda"){
    if(uMt>=6) return NULL;
    return fCk_pL[uMt];
  }
  else{
    printf("\033[1;31mERROR:\033[0m CecaAnalysis1::GetFit got an unkown system (%s)\n",system.Data());
    return NULL;
  }
}

TH1F* CecaAnalysis1::GetData(TString system, unsigned uMt){
  if(system=="pp"){
    if(uMt>=7) return NULL;
    return hCkExp_pp[uMt];
  }
  else if(system=="pL"||system=="pLambda"){
    if(uMt>=6) return NULL;
    return hCkExp_pL[uMt];
  }
  else{
    printf("\033[1;31mERROR:\033[0m CecaAnalysis1::GeData got an unkown system (%s)\n",system.Data());
    return NULL;
  }
}

void CecaAnalysis1::DumpCurrentCk(TString OutFileBase, int plots){

  TFile fOutput(OutFileBase+".root","recreate");

  TH1F* hAxisCk_pp = new TH1F("hAxisCk_pp","hAxisCk_pp",128,0,settings[FemtoRegMax_pp]);
  hAxisCk_pp->SetStats(false);
  hAxisCk_pp->SetTitle("");
  hAxisCk_pp->GetXaxis()->SetTitle("k* (MeV)");
  hAxisCk_pp->GetXaxis()->SetTitleSize(0.06);
  hAxisCk_pp->GetXaxis()->SetLabelSize(0.06);
  hAxisCk_pp->GetXaxis()->CenterTitle();
  hAxisCk_pp->GetXaxis()->SetTitleOffset(1.3);
  hAxisCk_pp->GetXaxis()->SetLabelOffset(0.02);

  hAxisCk_pp->GetYaxis()->SetTitle("C(k*)");
  hAxisCk_pp->GetYaxis()->SetTitleSize(0.06);
  hAxisCk_pp->GetYaxis()->SetLabelSize(0.06);
  hAxisCk_pp->GetYaxis()->CenterTitle();
  hAxisCk_pp->GetYaxis()->SetTitleOffset(1.10);
  //hAxisCk->GetXaxis()->SetNdivisions(506);
  if(plots==1) hAxisCk_pp->GetYaxis()->SetRangeUser(0.8, 4.8);
  else hAxisCk_pp->GetYaxis()->SetRangeUser(0.9, 1.15);

  TCanvas** cCk_pp = new TCanvas* [NumMtBins_pp];
  for(unsigned uMt=0; uMt<NumMtBins_pp; uMt++){
    cCk_pp[uMt] = new TCanvas(TString::Format("cCk_pp[%u]",uMt), TString::Format("cCk_pp[%u]",uMt), 1);
    cCk_pp[uMt]->cd(0);
    cCk_pp[uMt]->SetCanvasSize(1280, 720);
    cCk_pp[uMt]->SetMargin(0.15,0.05,0.2,0.05);//lrbt
  }

  TH1F* hAxisCk_pL = new TH1F("hAxisCk_pL","hAxisCk_pL",128,0,settings[FemtoRegMax_pL]);
  hAxisCk_pL->SetStats(false);
  hAxisCk_pL->SetTitle("");
  hAxisCk_pL->GetXaxis()->SetTitle("k* (MeV)");
  hAxisCk_pL->GetXaxis()->SetTitleSize(0.06);
  hAxisCk_pL->GetXaxis()->SetLabelSize(0.06);
  hAxisCk_pL->GetXaxis()->CenterTitle();
  hAxisCk_pL->GetXaxis()->SetTitleOffset(1.3);
  hAxisCk_pL->GetXaxis()->SetLabelOffset(0.02);

  hAxisCk_pL->GetYaxis()->SetTitle("C(k*)");
  hAxisCk_pL->GetYaxis()->SetTitleSize(0.06);
  hAxisCk_pL->GetYaxis()->SetLabelSize(0.06);
  hAxisCk_pL->GetYaxis()->CenterTitle();
  hAxisCk_pL->GetYaxis()->SetTitleOffset(1.10);
  //hAxisCk->GetXaxis()->SetNdivisions(506);
  if(plots==1) hAxisCk_pL->GetYaxis()->SetRangeUser(0.9, 2.1);
  else hAxisCk_pL->GetYaxis()->SetRangeUser(0.95, 1.1);

  TCanvas** cCk_pL = new TCanvas* [NumMtBins_pL];
  for(unsigned uMt=0; uMt<NumMtBins_pL; uMt++){
    cCk_pL[uMt] = new TCanvas(TString::Format("cCk_pL[%u]",uMt), TString::Format("cCk_pL[%u]",uMt), 1);
    cCk_pL[uMt]->cd(0);
    cCk_pL[uMt]->SetCanvasSize(1280, 720);
    cCk_pL[uMt]->SetMargin(0.15,0.05,0.2,0.05);//lrbt
  }




  TH1F* hAxisMt = new TH1F("hAxisMt","hAxisMt",128,1.0,2.4);
  hAxisMt->SetStats(false);
  hAxisMt->SetTitle("");
  hAxisMt->GetXaxis()->SetTitle("m_{T} (GeV)");
  hAxisMt->GetXaxis()->SetTitleSize(0.06);
  hAxisMt->GetXaxis()->SetLabelSize(0.06);
  hAxisMt->GetXaxis()->CenterTitle();
  hAxisMt->GetXaxis()->SetTitleOffset(1.3);
  hAxisMt->GetXaxis()->SetLabelOffset(0.02);

  hAxisMt->GetYaxis()->SetTitle("<r*> #approx 2.3 x r_{eff} (fm)");
  hAxisMt->GetYaxis()->SetTitleSize(0.06);
  hAxisMt->GetYaxis()->SetLabelSize(0.06);
  hAxisMt->GetYaxis()->CenterTitle();
  hAxisMt->GetYaxis()->SetTitleOffset(1.10);
  //hAxisMt->GetXaxis()->SetNdivisions(506);
  hAxisMt->GetYaxis()->SetRangeUser(2.0, 3.6);


  TGraphErrors mT_Scaling_pp;
  mT_Scaling_pp.SetName("mT_Scaling_pp");

  TGraphErrors mT_Scaling_pL;
  mT_Scaling_pL.SetName("mT_Scaling_pL");

  for(unsigned uMt=0; uMt<NumMtBins_pp; uMt++){
    if(hCkExp_pp){
      if(hCkExp_pp[uMt]){

        hCkExp_pp[uMt]->SetLineWidth(4);
        hCkExp_pp[uMt]->SetLineColor(kBlack);

        hCkExp_pp[uMt]->Write();

        cCk_pp[uMt]->cd();
        hAxisCk_pp->Draw("axis");
        hCkExp_pp[uMt]->Draw("same");

      }
    }
    if(fCk_pp){
      if(fCk_pp[uMt]){
        fCk_pp[uMt]->SetLineWidth(3);
        fCk_pp[uMt]->SetLineColor(kBlue+1);

        fSrc_pp[uMt]->SetLineWidth(3);
        fSrc_pp[uMt]->SetLineColor(kBlue+1);

        fCk_pp[uMt]->Write();
        cCk_pp[uMt]->cd();
        fCk_pp[uMt]->Draw("same");
        //TF1* fSrcC_pp = new TF1("fSrcC_pp",Src_pp,&DLM_CecaSource_v0::RootEval,0,16,5,"DLM_CecaSource_v0","RootEval");

        for(unsigned uP=0; uP<5; uP++){
          fSrc_pp[uMt]->FixParameter(uP,fCk_pp[uMt]->GetParameter(uP));
        }
        fSrc_pp[uMt]->Write();
        double Mean = fSrc_pp[uMt]->Mean(0, 16);
        mT_Scaling_pp.SetPoint(uMt,mT_BinCenter_pp[uMt]*0.001,Mean);
        //delete fSrcC_pp;
      }

      //SOURCE
      //fSrc_pp[uMt]->SetParameter(0,MtVal);
      //fSrc_pp[uMt]->SetParameter(1,fCk_pp[uMt]->GetParameter(1));
      //fSrc_pp[uMt]->SetParameter(2,fCk_pp[uMt]->GetParameter(2));
      //fSrc_pp[uMt]->SetParameter(3,fCk_pp[uMt]->GetParameter(3));
      //fSrc_pp[uMt]->FixParameter(4,fCk_pp[uMt]->GetParameter(4));
      //fSrc_pp[uMt]->SetLineColor(kBlue);
      //fSrc_pp[uMt]->SetNpx(2048);

    }
    if(plots) cCk_pp[uMt]->SaveAs(TString::Format(OutFileBase+"_Ck_pp_Mt%u_p%i.pdf",uMt,plots));
  }//NumMtBins_pp
  for(unsigned uMt=0; uMt<NumMtBins_pL; uMt++){
    if(hCkExp_pL){
      if(hCkExp_pL[uMt]){
        hCkExp_pL[uMt]->SetLineWidth(4);
        hCkExp_pL[uMt]->SetLineColor(kBlack);

        hCkExp_pL[uMt]->Write();

        cCk_pL[uMt]->cd();
        hAxisCk_pL->Draw("axis");
        hCkExp_pL[uMt]->Draw("same");
      }
    }
    if(fCk_pL){
      if(fCk_pL[uMt]){
        fCk_pL[uMt]->SetLineWidth(3);
        fCk_pL[uMt]->SetLineColor(kRed+2);

        fCk_pL[uMt]->Write();
        cCk_pL[uMt]->cd();
        fCk_pL[uMt]->Draw("same");
        //TF1* fSrcC_pL = new TF1("fSrcC_pL",Src_pL,&DLM_CecaSource_v0::RootEval,0,16,5,"DLM_CecaSource_v0","RootEval");
        fSrc_pL[uMt]->SetLineWidth(3);
        fSrc_pL[uMt]->SetLineColor(kRed+2);
        for(unsigned uP=0; uP<5; uP++){
          fSrc_pL[uMt]->FixParameter(uP,fCk_pL[uMt]->GetParameter(uP));
        }
        fSrc_pL[uMt]->Write();
        double Mean = fSrc_pL[uMt]->Mean(0, 16);
        mT_Scaling_pL.SetPoint(uMt,mT_BinCenter_pL[uMt]*0.001,Mean);
        //delete fSrcC_pL;
      }
    }
    if(plots) cCk_pL[uMt]->SaveAs(TString::Format(OutFileBase+"_Ck_pL_Mt%u_p%i.pdf",uMt,plots));
  }//NumMtBins_pL


  TGraphErrors mT_Published_pp;
  mT_Published_pp.SetName("mT_Published_pp");
  TGraphErrors mT_Published_pL;
  mT_Published_pL.SetName("mT_Published_pL");

  TFile fInput(catsfiles_fld+"/Source/SourcePaper_Published.root","read");
  TGraphErrors* g_pp = (TGraphErrors*)fInput.Get("g_reff_pp");
  TGraphErrors* g_pL = (TGraphErrors*)fInput.Get("g_reff_pL");

  for(unsigned uMt=0; uMt<g_pp->GetN(); uMt++){
    double MtVal,Reff;
    g_pp->GetPoint(uMt,MtVal,Reff);
    double dReff = g_pp->GetErrorY(uMt);

    double UpReff = Reff+dReff;
    double LowReff = Reff-dReff;

    double UpMean = ConvertGaussMean(UpReff);
    double LowMean = ConvertGaussMean(LowReff);
    double MeanMean = (UpMean+LowMean)*0.5;

    mT_Published_pp.SetPoint(uMt,MtVal,MeanMean);
    mT_Published_pp.SetPointError(uMt,0,(UpMean-LowMean)*0.5);
  }

  for(unsigned uMt=0; uMt<g_pL->GetN(); uMt++){
    double MtVal,Reff;
    g_pL->GetPoint(uMt,MtVal,Reff);
    double dReff = g_pL->GetErrorY(uMt);

    double UpReff = Reff+dReff;
    double LowReff = Reff-dReff;

    double UpMean = ConvertGaussMean(UpReff);
    double LowMean = ConvertGaussMean(LowReff);
    double MeanMean = (UpMean+LowMean)*0.5;

    mT_Published_pL.SetPoint(uMt,MtVal,MeanMean);
    mT_Published_pL.SetPointError(uMt,0,(UpMean-LowMean)*0.5);
  }

  mT_Published_pp.SetLineWidth(4);
  mT_Published_pp.SetLineColor(kBlue+1);

  mT_Scaling_pp.SetLineWidth(4);
  mT_Scaling_pp.SetLineColor(kAzure+1);

  mT_Published_pL.SetLineWidth(4);
  mT_Published_pL.SetLineColor(kRed+2);

  mT_Scaling_pL.SetLineWidth(4);
  mT_Scaling_pL.SetLineColor(kRed);


  fOutput.cd();
  mT_Scaling_pp.Write();
  mT_Published_pp.Write();

  mT_Scaling_pL.Write();
  mT_Published_pL.Write();


  TLegend lMt(0.65,0.75,0.95,0.95);//lbrt
  lMt.SetName("lMt");
  lMt.SetTextSize(0.045);
  lMt.AddEntry(&mT_Published_pp, "Published pp");
  if(OutFileBase.Contains("Ceca")) lMt.AddEntry(&mT_Scaling_pp, "CECA pp");
  else if(OutFileBase.Contains("Gauss")) lMt.AddEntry(&mT_Scaling_pp, "Gauss pp");
  lMt.AddEntry(&mT_Published_pL, "Published p#Lambda");
  if(OutFileBase.Contains("Ceca")) lMt.AddEntry(&mT_Scaling_pL, "CECA p#Lambda");
  else if(OutFileBase.Contains("Gauss")) lMt.AddEntry(&mT_Scaling_pL, "Gauss p#Lambda");














  TCanvas cMt("cMt","cMt");
  cMt.cd(0);
  cMt.SetCanvasSize(1280, 720);
  cMt.SetMargin(0.15,0.05,0.2,0.05);//lrbt
  hAxisMt->Draw("axis");
  mT_Published_pp.Draw("same");
  mT_Published_pL.Draw("same");
  mT_Scaling_pp.Draw("same");
  mT_Scaling_pL.Draw("same");
  lMt.Draw("same");
  cMt.SaveAs(OutFileBase+"_Mt.pdf");


  delete hAxisCk_pp;
  delete hAxisCk_pL;
  delete hAxisMt;
  for(unsigned uMt=0; uMt<NumMtBins_pp; uMt++){
    delete cCk_pp[uMt];
  }
  delete [] cCk_pp;
  for(unsigned uMt=0; uMt<NumMtBins_pL; uMt++){
    delete cCk_pL[uMt];
  }
  delete [] cCk_pL;
}


void CecaAnalysis1::SetUp_Fits(TString fittype, bool AutoSrcPars, bool AutoBlPars){

  if(fittype=="Ceca"||fittype=="CECA"||fittype=="ceca"){
    FittingScenario = 0;
  }
  else if(fittype=="SingleCeca"||fittype=="SINGLECECA"||fittype=="singlececa"){
    FittingScenario = 1;
  }
  else if(fittype=="Gauss"||fittype=="GAUSS"||fittype=="gauss"||fittype=="Gaus"||fittype=="GAUS"||fittype=="gaus"){
    FittingScenario = 11;
  }
  else{
    printf("\033[1;31mERROR:\033[0m SetUp_Fits got an unkown fit type (%s)\n",fittype.Data());
    return;
  }

  printf("   Setting up the pp fits...\n");
  if(fCk_pp){
    for(unsigned uMt=0; uMt<NumMtBins_pp; uMt++){
      if(fCk_pp[uMt]){
        delete fCk_pp[uMt];
        fCk_pp[uMt] = NULL;
      }
    }
    delete [] fCk_pp;
    fCk_pp = NULL;
  }
  if(fSrc_pp){
    for(unsigned uMt=0; uMt<NumMtBins_pp; uMt++){
      if(fSrc_pp[uMt]){
        delete fSrc_pp[uMt];
        fSrc_pp[uMt] = NULL;
      }
    }
    delete [] fSrc_pp;
    fSrc_pp = NULL;
  }
  fCk_pp = new TF1* [NumMtBins_pp];
  fSrc_pp = new TF1* [NumMtBins_pp];
  for(unsigned uMt=0; uMt<NumMtBins_pp; uMt++){
    fCk_pp[uMt] = new TF1(TString::Format("fCk_pp[%u]",uMt),this,&CecaAnalysis1::FitFun_pp,hCkExp_pp[uMt]->GetBinLowEdge(1),settings[FitMax_pp],TotFitPars,"CecaAnalysis1","FitFun_pp");
    fSrc_pp[uMt] = new TF1(TString::Format("fSrc_pp[%u]",uMt),Src_pp,&DLM_CecaSource_v0::RootEval,0,16,5,"DLM_CecaSource_v0","RootEval");

//printf("FIT%u %f %f\n",uMt,hCkExp_pp[uMt]->GetBinLowEdge(1),settings[FitMax_pp]);
    //the mT is fixed
    if(FittingScenario==11){
      fCk_pp[uMt]->FixParameter(0,-1);
    }
    else{
      fCk_pp[uMt]->FixParameter(0,mT_BinCenter_pp[uMt]);
    }

    //the source parameters
    if(AutoSrcPars){
      if(FittingScenario==11){
        fCk_pp[uMt]->SetParameter(1,1.0);
        fCk_pp[uMt]->SetParLimits(1,0.5,2.0);
      }
      else{
        for(unsigned uP=1; uP<4; uP++){
          double Low = Src_pp->Low_par(uP,true);
          double Up = Src_pp->Up_par(uP,true);
          fCk_pp[uMt]->SetParameter(uP,(Low+Up)*0.5);
          fCk_pp[uMt]->SetParLimits(uP,Low,Up);
        }
      }
    }
    //the scale parameter is fixed for the fit
    fCk_pp[uMt]->FixParameter(4,1);

    //the baseline
    if(AutoBlPars){
      //norm
      fCk_pp[uMt]->SetParameter(5,1);
      fCk_pp[uMt]->SetParLimits(5,0.9,1.1);

      //flat at zero
      fCk_pp[uMt]->FixParameter(6,0);

      //the position of the extremum
      fCk_pp[uMt]->SetParameter(7,100);
      fCk_pp[uMt]->SetParLimits(7,20,800);
      if(ana_ver=="Reduced") fCk_pp[uMt]->FixParameter(7,0);

      //the pol3
      fCk_pp[uMt]->SetParameter(8,0);
      fCk_pp[uMt]->SetParLimits(8,-1e-10,1e-10);
      if(ana_ver=="Reduced") fCk_pp[uMt]->FixParameter(8,0);

      //remove the 4th order
      fCk_pp[uMt]->FixParameter(9,0);
    }

  }//uMt_pp

  printf("   Setting up the pL fits...\n");
  if(fCk_pL){
    for(unsigned uMt=0; uMt<NumMtBins_pL; uMt++){
      if(fCk_pL[uMt]){
        delete fCk_pL[uMt];
        fCk_pL[uMt] = NULL;
      }
    }
    delete [] fCk_pL;
    fCk_pL = NULL;
  }
  if(fSrc_pL){
    for(unsigned uMt=0; uMt<NumMtBins_pL; uMt++){
      if(fSrc_pL[uMt]){
        delete fSrc_pL[uMt];
        fSrc_pL[uMt] = NULL;
      }
    }
    delete [] fSrc_pL;
    fSrc_pL = NULL;
  }
  fCk_pL = new TF1* [NumMtBins_pL];
  fSrc_pL = new TF1* [NumMtBins_pL];
  for(unsigned uMt=0; uMt<NumMtBins_pL; uMt++){
    if(fCk_pL[uMt]){
      delete fCk_pL[uMt];
      fCk_pL[uMt] = NULL;
    }

    fCk_pL[uMt] = new TF1(TString::Format("fCk_pL[%u]",uMt),this,&CecaAnalysis1::FitFun_pL,hCkExp_pL[uMt]->GetBinLowEdge(1),settings[FitMax_pL],TotFitPars,"CecaAnalysis1","FitFun_pL");
    fSrc_pL[uMt] = new TF1(TString::Format("fSrc_pL[%u]",uMt),Src_pL,&DLM_CecaSource_v0::RootEval,0,16,5,"DLM_CecaSource_v0","RootEval");

    //the mT is fixed
    if(FittingScenario==11){
      fCk_pL[uMt]->FixParameter(0,-1);
    }
    else{
      fCk_pL[uMt]->FixParameter(0,mT_BinCenter_pL[uMt]);
    }

    //the source parameters
    if(AutoSrcPars){
      if(FittingScenario==11){
        fCk_pL[uMt]->SetParameter(1,1.0);
        fCk_pL[uMt]->SetParLimits(1,0.5,2.0);
      }
      else{
        for(unsigned uP=1; uP<4; uP++){
          double Low = Src_pL->Low_par(uP,true);
          double Up = Src_pL->Up_par(uP,true);
          fCk_pL[uMt]->SetParameter(uP,(Low+Up)*0.5);
          fCk_pL[uMt]->SetParLimits(uP,Low,Up);
        }
      }
    }
    //the scale parameter is fixed for the fit
    fCk_pL[uMt]->FixParameter(4,1);

    //the baseline
    if(AutoBlPars){
      //norm
      fCk_pL[uMt]->SetParameter(5,1);
      fCk_pL[uMt]->SetParLimits(5,0.9,1.1);

      //flat at zero
      fCk_pL[uMt]->FixParameter(6,0);

      //the position of the extremum
      fCk_pL[uMt]->SetParameter(7,100);
      fCk_pL[uMt]->SetParLimits(7,20,800);
      if(ana_ver=="Reduced") fCk_pL[uMt]->FixParameter(7,0);

      //the pol3
      fCk_pL[uMt]->SetParameter(8,0);
      fCk_pL[uMt]->SetParLimits(8,-1e-10,1e-10);
      if(ana_ver=="Reduced") fCk_pL[uMt]->FixParameter(8,0);

      //remove the 4th order
      fCk_pL[uMt]->FixParameter(9,0);
    }

  }//uMt_pL

}

void CecaAnalysis1::GoBabyGo(bool print_info){
  if(print_info) printf("<<< Go baby go >>>\n");
  InitHistory();
  //GLOBAL FIT
  if(FittingScenario==0){
    if(GrandeFitter){
      delete GrandeFitter;
      GrandeFitter = NULL;
    }

    if(print_info) printf("   Setting up the global fit...\n");
    GrandeFitter = new DLM_MultiFit();
    for(unsigned uMt=0; uMt<NumMtBins_pp; uMt++){
      GrandeFitter->AddSpectrum(hCkExp_pp[uMt],fCk_pp[uMt]);
    }
    for(unsigned uMt=0; uMt<NumMtBins_pL; uMt++){
      GrandeFitter->AddSpectrum(hCkExp_pL[uMt],fCk_pL[uMt]);
    }
  //DEB
    for(unsigned uMt=1; uMt<NumMtBins_pp+NumMtBins_pL; uMt++){
      GrandeFitter->SetEqualPar(uMt,1,0,1);//d
      GrandeFitter->SetEqualPar(uMt,2,0,2);//ht
      GrandeFitter->SetEqualPar(uMt,3,0,3);//hz or tau
    }

    if(print_info) printf("Well... let us fit!!!\n");
    if(print_info) printf(" <a long long wait ahead>\n");
    FunctionCallCounter = 0;
    ROOT::Fit::FitResult fitresptr = GrandeFitter->PerformGlobalFit(print_info);
/*
    TRandom3 rangen(11);
    while(true){
      double drnd = rangen.Uniform(1.1,1.3);

      Ck_pp->SetSourcePar(0,-1);
      Ck_pp->SetSourcePar(1,drnd);

      Ck_pL->SetSourcePar(0,-1);
      Ck_pL->SetSourcePar(1,drnd);

      Ck_pS0->SetSourcePar(0,-1);
      Ck_pS0->SetSourcePar(1,drnd);

      Ck_pXim->SetSourcePar(0,-1);
      Ck_pXim->SetSourcePar(1,drnd);

      Ck_pXi0->SetSourcePar(0,-1);
      Ck_pXi0->SetSourcePar(1,drnd);

      Ck_pXim1530->SetSourcePar(0,-1);
      Ck_pXim1530->SetSourcePar(1,drnd);

      //Ck_pp->Update(true);NO LEAK
      //CkDec_pp->Update(true,true);//LEAK
      //CkDec_pp->Update(true,false);//
      //CkDec_pp->Update(false,true);//
      //CkDec_pp->Update(false,false);//LEAK
      //CkDec_pXim1530->Update(true,true);//NO LEAK

      NewPar[0] = -1;
      NewPar[1] = drnd;
      NewPar[2] = 2;
      NewPar[3] = 6;
      UpdateDecomp();
    }
*/
  }
  else{
    for(unsigned uMt=0; uMt<NumMtBins_pp; uMt++){
      if(print_info) printf(" fitting pp[%u]\n",uMt);
      if(print_info) hCkExp_pp[uMt]->Fit(fCk_pp[uMt],"S, N, R, M");
      else hCkExp_pp[uMt]->Fit(fCk_pp[uMt],"Q, S, N, R, M");
    }
    for(unsigned uMt=0; uMt<NumMtBins_pL; uMt++){
      if(print_info) printf(" fitting pL[%u]\n",uMt);
      if(print_info) hCkExp_pL[uMt]->Fit(fCk_pL[uMt],"S, N, R, M");
      else hCkExp_pL[uMt]->Fit(fCk_pL[uMt],"Q, S, N, R, M");
    }
  }
}

//With a totally free Gauss, chi2 of 748.

void TestSetUpAna(){


  TString Description = "J1_Reduced_USM_Ceca";
  //TString Description = "J1_Reduced_USM_Gauss";
  //TString Description = "NLO_Gauss";//NLO19-600
  //TString Description = "Usmani_Gauss";//usmani 1:1 as NLO19-600
  //TString Description = "C2_Reduced_USM_Ceca";

  TString CecaAnaSettings;
  TString FileBase;
  if(Description=="C2_Reduced_USM_Ceca"){
    FileBase = Description+"/Ceca";
    CecaAnaSettings = "Cigar2_ds24_hts36_hzs36";
  }
  else if(Description=="J1_Reduced_USM_Ceca"){
    FileBase = Description+"/Ceca";
    CecaAnaSettings = "Jaime1_ds24_hts36_hzs36";
  }
  else if(Description=="J1_Reduced_USM_Gauss"){
    FileBase = Description+"/Gauss";
    CecaAnaSettings = "Jaime1_ds24_hts36_hzs36";
  }
  else if(Description=="NLO_Gauss"){
    FileBase = Description+"/NLO_Gauss";
    CecaAnaSettings = "Jaime1_ds24_hts36_hzs36";//should not matter which
  }
  else if(Description=="Usmani_Gauss"){
    FileBase = Description+"/Usmani_Gauss";
    CecaAnaSettings = "Jaime1_ds24_hts36_hzs36";//should not matter which
  }
  else{
    return;
  }

  CecaAnalysis1 CECA_ANA("Reduced",CecaAnaSettings,GetCernBoxDimi()+TString("/CatsFiles/"));



  CECA_ANA.SetUp_pp("AV18",0);

  //CECA_ANA.SetUp_pL("Chiral_Coupled_SPD",9411602);
  //CECA_ANA.SetUp_pL("Chiral_Coupled_SPD",11600);

//Wc=2279.0; Rc=0.3394; Sc=0.2614; f1=1.41; d1=2.53; tDev=0.003
  if(Description=="NLO_Gauss"){
    CECA_ANA.SetUp_pL("Chiral_Coupled_SPD",11600);
  }
  else{
    CECA_ANA.SetUp_pL("UsmaniFit",0);

    //1:1 as NLO600
    if(Description=="Usmani_Gauss"){
      CECA_ANA.Kitty_pL->SetShortRangePotential(1,0,1,2279.0);
      CECA_ANA.Kitty_pL->SetShortRangePotential(1,0,2,0.3394);
      CECA_ANA.Kitty_pL->SetShortRangePotential(1,0,3,0.2614);
    }
    //cigar2 pars
    else if(Description=="C2_Reduced_USM_Ceca"){
      //at the moment the same as best jaime pars
      CECA_ANA.Kitty_pL->SetShortRangePotential(1,0,1,2279.1);
      CECA_ANA.Kitty_pL->SetShortRangePotential(1,0,2,0.34770);
      CECA_ANA.Kitty_pL->SetShortRangePotential(1,0,3,0.26002);
    }
    //best jaime pars
    else{
      CECA_ANA.Kitty_pL->SetShortRangePotential(1,0,1,2279.1);
      CECA_ANA.Kitty_pL->SetShortRangePotential(1,0,2,0.34770);
      CECA_ANA.Kitty_pL->SetShortRangePotential(1,0,3,0.26002);
    }


  }


  CECA_ANA.SetUp_pS0("Chiral",0);
  CECA_ANA.SetUp_pXim("pXim_HALQCDPaper2020",0);
  CECA_ANA.SetUp_pXi0("pXim_HALQCDPaper2020",0);
  CECA_ANA.SetUp_Decomposition(0,0);
  if(Description=="J1_Reduced_USM_Ceca"||Description=="C2_Reduced_USM_Ceca")
    CECA_ANA.SetUp_Fits("SingleCeca");//each mT bin by itself
  else if(Description=="J1_Reduced_USM_Gauss")
    CECA_ANA.SetUp_Fits("Gauss");
  else if(Description=="NLO_Gauss"||Description=="Usmani_Gauss")
    CECA_ANA.SetUp_Fits("Gauss");
  else return;
  //CECA_ANA.SetUp_Fits("Ceca");//Global fit
  //CECA_ANA.SetUp_Fits("Gauss");//each mT bin by itself


  CECA_ANA.settings[CecaAnalysis1::FemtoRegMax_pp] = 184;
  CECA_ANA.settings[CecaAnalysis1::FitMax_pp] = 184;
  CECA_ANA.settings[CecaAnalysis1::CkConv_pp] = 700;

  CECA_ANA.settings[CecaAnalysis1::FemtoRegMax_pL] = 180;
  CECA_ANA.settings[CecaAnalysis1::FitMax_pL] = 180;
  CECA_ANA.settings[CecaAnalysis1::CkConv_pL] = 700;

  //CECA_ANA.GetFit("pp",0)->FixParameter(0,1.10770e+03);
  //CECA_ANA.GetFit("pp",1)->FixParameter(0,1.16830e+03);
  //CECA_ANA.GetFit("pp",2)->FixParameter(0,1.22840e+03);
  //CECA_ANA.GetFit("pp",3)->FixParameter(0,1.31560e+03);
  //CECA_ANA.GetFit("pp",4)->FixParameter(0,1.46280e+03);
  //CECA_ANA.GetFit("pp",5)->FixParameter(0,1.68720e+03);
  //CECA_ANA.GetFit("pp",6)->FixParameter(0,2.21160e+03);
  //CECA_ANA.GetFit("pL",0)->FixParameter(0,1.21240e+03);
  //CECA_ANA.GetFit("pL",1)->FixParameter(0,1.28960e+03);
  //CECA_ANA.GetFit("pL",2)->FixParameter(0,1.37600e+03);
  //CECA_ANA.GetFit("pL",3)->FixParameter(0,1.54070e+03);
  //CECA_ANA.GetFit("pL",4)->FixParameter(0,1.75600e+03);
  //CECA_ANA.GetFit("pL",5)->FixParameter(0,2.25940e+03);

  if(false){//Jaime pars
    for(unsigned uMt=0; uMt<7; uMt++){
      CECA_ANA.GetFit("pp",uMt)->FixParameter(1,0.2);
      CECA_ANA.GetFit("pp",uMt)->FixParameter(2,2.3);
      CECA_ANA.GetFit("pp",uMt)->FixParameter(3,3.5);
      CECA_ANA.GetFit("pp",uMt)->FixParameter(4,1);

      //CECA_ANA.GetFit("pp",uMt)->FixParameter(5,1);
      //CECA_ANA.GetFit("pp",uMt)->FixParameter(6,0);
      //CECA_ANA.GetFit("pp",uMt)->FixParameter(7,0);
      //CECA_ANA.GetFit("pp",uMt)->FixParameter(8,0);
      CECA_ANA.GetFit("pp",uMt)->FixParameter(9,0);
    }
    for(unsigned uMt=0; uMt<6; uMt++){
      CECA_ANA.GetFit("pL",uMt)->FixParameter(1,0.2*1.0);
      CECA_ANA.GetFit("pL",uMt)->FixParameter(2,2.3*1.0);
      CECA_ANA.GetFit("pL",uMt)->FixParameter(3,3.5*1.0);
      CECA_ANA.GetFit("pL",uMt)->FixParameter(4,1);

      //CECA_ANA.GetFit("pL",uMt)->FixParameter(5,1);
      //CECA_ANA.GetFit("pL",uMt)->FixParameter(6,0);
      //CECA_ANA.GetFit("pL",uMt)->FixParameter(7,0);
      //CECA_ANA.GetFit("pL",uMt)->FixParameter(8,0);
      CECA_ANA.GetFit("pL",uMt)->FixParameter(9,0);
    }
  }
  else if(Description=="J1_Reduced_USM_Ceca"){//the first best sol I found (Jaime1)
    for(unsigned uMt=0; uMt<7; uMt++){
      CECA_ANA.GetFit("pp",uMt)->FixParameter(1,0.31218);
      CECA_ANA.GetFit("pp",uMt)->FixParameter(2,3.5202);
      CECA_ANA.GetFit("pp",uMt)->FixParameter(3,2.64406);
      CECA_ANA.GetFit("pp",uMt)->FixParameter(4,1);

      CECA_ANA.GetFit("pp",uMt)->FixParameter(9,0);
    }
    for(unsigned uMt=0; uMt<6; uMt++){
      CECA_ANA.GetFit("pL",uMt)->FixParameter(1,0.31218);
      CECA_ANA.GetFit("pL",uMt)->FixParameter(2,3.5202);
      CECA_ANA.GetFit("pL",uMt)->FixParameter(3,2.64406);
      CECA_ANA.GetFit("pL",uMt)->FixParameter(4,1);

      CECA_ANA.GetFit("pL",uMt)->FixParameter(9,0);
    }
  }
  else if(false){//the first best sol I found (Jaime1), but for pp alone
    for(unsigned uMt=0; uMt<7; uMt++){
      CECA_ANA.GetFit("pp",uMt)->FixParameter(1,0.2319);
      CECA_ANA.GetFit("pp",uMt)->FixParameter(2,2.754);
      CECA_ANA.GetFit("pp",uMt)->FixParameter(3,3.403);
      CECA_ANA.GetFit("pp",uMt)->FixParameter(4,1);

      CECA_ANA.GetFit("pp",uMt)->FixParameter(9,0);
    }
    for(unsigned uMt=0; uMt<6; uMt++){
      CECA_ANA.GetFit("pL",uMt)->FixParameter(1,0.2319);
      CECA_ANA.GetFit("pL",uMt)->FixParameter(2,2.754);
      CECA_ANA.GetFit("pL",uMt)->FixParameter(3,3.403);
      CECA_ANA.GetFit("pL",uMt)->FixParameter(4,1);

      CECA_ANA.GetFit("pL",uMt)->FixParameter(9,0);
    }
  }
  else if(false){//the first best sol I found (Cigar2)
    for(unsigned uMt=0; uMt<7; uMt++){
      CECA_ANA.GetFit("pp",uMt)->FixParameter(1,0.4742);
      CECA_ANA.GetFit("pp",uMt)->FixParameter(2,3.285);
      CECA_ANA.GetFit("pp",uMt)->FixParameter(3,9.758);
      CECA_ANA.GetFit("pp",uMt)->FixParameter(4,1);

      CECA_ANA.GetFit("pp",uMt)->FixParameter(9,0);
    }
    for(unsigned uMt=0; uMt<6; uMt++){
      CECA_ANA.GetFit("pL",uMt)->FixParameter(1,0.4742);
      CECA_ANA.GetFit("pL",uMt)->FixParameter(2,3.285);
      CECA_ANA.GetFit("pL",uMt)->FixParameter(3,9.758);
      CECA_ANA.GetFit("pL",uMt)->FixParameter(4,1);

      CECA_ANA.GetFit("pL",uMt)->FixParameter(9,0);
    }
  }
  else if(false){//best cigar2 pars (954 chi2 from the scan)
    for(unsigned uMt=0; uMt<7; uMt++){
      CECA_ANA.GetFit("pp",uMt)->FixParameter(1,5.081e-01);
      CECA_ANA.GetFit("pp",uMt)->FixParameter(2,3.913e+00);
      CECA_ANA.GetFit("pp",uMt)->FixParameter(3,9.112e+00);
      CECA_ANA.GetFit("pp",uMt)->FixParameter(4,1);

      CECA_ANA.GetFit("pp",uMt)->FixParameter(9,0);
    }
    for(unsigned uMt=0; uMt<6; uMt++){
      CECA_ANA.GetFit("pL",uMt)->FixParameter(1,5.081e-01);
      CECA_ANA.GetFit("pL",uMt)->FixParameter(2,3.913e+00);
      CECA_ANA.GetFit("pL",uMt)->FixParameter(3,9.112e+00);
      CECA_ANA.GetFit("pL",uMt)->FixParameter(4,1);

      CECA_ANA.GetFit("pL",uMt)->FixParameter(9,0);
    }
  }
  else if(false){
    for(unsigned uMt=0; uMt<7; uMt++){
      CECA_ANA.GetFit("pp",uMt)->SetParameter(1,4.533e-01);
      CECA_ANA.GetFit("pp",uMt)->SetParLimits(1,0.9*4.533e-01,1.1*4.533e-01);
      CECA_ANA.GetFit("pp",uMt)->SetParameter(2,4.333e+00);
      CECA_ANA.GetFit("pp",uMt)->SetParLimits(2,0.9*4.333e+00,1.1*4.333e+00);
      CECA_ANA.GetFit("pp",uMt)->SetParameter(3,9.500e+00);
      CECA_ANA.GetFit("pp",uMt)->SetParLimits(3,0.9*9.500e+00,1.1*9.500e+00);
      //CECA_ANA.GetFit("pp",uMt)->SetParameter(3,8.0);
      //CECA_ANA.GetFit("pp",uMt)->SetParLimits(3,2.0,11.9);
    }
    for(unsigned uMt=0; uMt<6; uMt++){
      CECA_ANA.GetFit("pL",uMt)->SetParameter(1,4.533e-01);
      CECA_ANA.GetFit("pL",uMt)->SetParLimits(1,0.9*4.533e-01,1.1*4.533e-01);
      CECA_ANA.GetFit("pL",uMt)->SetParameter(2,4.333e+00);
      CECA_ANA.GetFit("pL",uMt)->SetParLimits(2,0.9*4.333e+00,1.1*4.333e+00);
      CECA_ANA.GetFit("pL",uMt)->SetParameter(3,9.500e+00);
      CECA_ANA.GetFit("pL",uMt)->SetParLimits(3,0.9*9.500e+00,1.1*9.500e+00);
      //CECA_ANA.GetFit("pL",uMt)->SetParameter(3,8.0);
      //CECA_ANA.GetFit("pL",uMt)->SetParLimits(3,2.0,11.9);
    }
  }
  else if(false){
    for(unsigned uMt=0; uMt<7; uMt++){
      CECA_ANA.GetFit("pp",uMt)->FixParameter(1,8.800e-01);
      CECA_ANA.GetFit("pp",uMt)->FixParameter(2,6.667e-02);
      CECA_ANA.GetFit("pp",uMt)->FixParameter(3,4.167e+00);
    }
    for(unsigned uMt=0; uMt<6; uMt++){
      CECA_ANA.GetFit("pL",uMt)->FixParameter(1,8.800e-01);
      CECA_ANA.GetFit("pL",uMt)->FixParameter(2,6.667e-02);
      CECA_ANA.GetFit("pL",uMt)->FixParameter(3,4.167e+00);
    }
  }
  else if(false){
    CECA_ANA.GetFit("pp",0)->FixParameter(0,1.10770e+03);
    CECA_ANA.GetFit("pp",0)->FixParameter(1,6.16591e-01);
    CECA_ANA.GetFit("pp",0)->FixParameter(2,2.86667e+00);
    CECA_ANA.GetFit("pp",0)->FixParameter(3,8.16667e+00);
    CECA_ANA.GetFit("pp",0)->FixParameter(4,1.0);

    CECA_ANA.GetFit("pp",1)->FixParameter(0,1.16830e+03);
    CECA_ANA.GetFit("pp",1)->FixParameter(1,7.73333e-01);
    CECA_ANA.GetFit("pp",1)->FixParameter(2,2.86764e+00);
    CECA_ANA.GetFit("pp",1)->FixParameter(3,6.50000e+00);
    CECA_ANA.GetFit("pp",1)->FixParameter(4,1.0);

    CECA_ANA.GetFit("pp",2)->FixParameter(0,1.22840e+03);
    CECA_ANA.GetFit("pp",2)->FixParameter(1,6.56113e-01);
    CECA_ANA.GetFit("pp",2)->FixParameter(2,2.46668e+00);
    CECA_ANA.GetFit("pp",2)->FixParameter(3,8.24924e+00);
    CECA_ANA.GetFit("pp",2)->FixParameter(4,1.0);

    CECA_ANA.GetFit("pp",3)->FixParameter(0,1.31560e+03);
    CECA_ANA.GetFit("pp",3)->FixParameter(1,6.39995e-01);
    CECA_ANA.GetFit("pp",3)->FixParameter(2,4.12169e+00);
    CECA_ANA.GetFit("pp",3)->FixParameter(3,6.00061e+00);
    CECA_ANA.GetFit("pp",3)->FixParameter(4,1.0);

    CECA_ANA.GetFit("pp",4)->FixParameter(0,1.46280e+03);
    CECA_ANA.GetFit("pp",4)->FixParameter(1,7.46824e-01);
    CECA_ANA.GetFit("pp",4)->FixParameter(2,2.46667e+00);
    CECA_ANA.GetFit("pp",4)->FixParameter(3,5.49974e+00);
    CECA_ANA.GetFit("pp",4)->FixParameter(4,1.0);

    CECA_ANA.GetFit("pp",5)->FixParameter(0,1.68720e+03);
    CECA_ANA.GetFit("pp",5)->FixParameter(1,7.06734e-01);
    CECA_ANA.GetFit("pp",5)->FixParameter(2,2.40001e+00);
    CECA_ANA.GetFit("pp",5)->FixParameter(3,5.57319e+00);
    CECA_ANA.GetFit("pp",5)->FixParameter(4,1.0);

    CECA_ANA.GetFit("pp",6)->FixParameter(0,2.21160e+03);
    CECA_ANA.GetFit("pp",6)->FixParameter(1,5.06662e-01);
    CECA_ANA.GetFit("pp",6)->FixParameter(2,2.33333e+00);
    CECA_ANA.GetFit("pp",6)->FixParameter(3,6.29996e+00);
    CECA_ANA.GetFit("pp",6)->FixParameter(4,1.0);




    CECA_ANA.GetFit("pL",0)->FixParameter(0,1.21240e+03);
    CECA_ANA.GetFit("pL",0)->FixParameter(1,8.20771e-01);
    CECA_ANA.GetFit("pL",0)->FixParameter(2,2.73333e+00);
    CECA_ANA.GetFit("pL",0)->FixParameter(3,7.97590e+00);
    CECA_ANA.GetFit("pL",0)->FixParameter(4,1.0);

    CECA_ANA.GetFit("pL",1)->FixParameter(0,1.28960e+03);
    CECA_ANA.GetFit("pL",1)->FixParameter(1,8.61206e-01);
    CECA_ANA.GetFit("pL",1)->FixParameter(2,6.66673e-02);
    CECA_ANA.GetFit("pL",1)->FixParameter(3,6.83333e+00);
    CECA_ANA.GetFit("pL",1)->FixParameter(4,1.0);

    CECA_ANA.GetFit("pL",2)->FixParameter(0,1.37600e+03);
    CECA_ANA.GetFit("pL",2)->FixParameter(1,8.89911e-01);
    CECA_ANA.GetFit("pL",2)->FixParameter(2,2.39784e+00);
    CECA_ANA.GetFit("pL",2)->FixParameter(3,6.83406e+00);
    CECA_ANA.GetFit("pL",2)->FixParameter(4,1.0);

    CECA_ANA.GetFit("pL",3)->FixParameter(0,1.54070e+03);
    CECA_ANA.GetFit("pL",3)->FixParameter(1,6.85418e-01);
    CECA_ANA.GetFit("pL",3)->FixParameter(2,3.40000e+00);
    CECA_ANA.GetFit("pL",3)->FixParameter(3,9.26255e+00);
    CECA_ANA.GetFit("pL",3)->FixParameter(4,1.0);

    CECA_ANA.GetFit("pL",4)->FixParameter(0,1.75600e+03);
    CECA_ANA.GetFit("pL",4)->FixParameter(1,6.95112e-01);
    CECA_ANA.GetFit("pL",4)->FixParameter(2,2.36499e+00);
    CECA_ANA.GetFit("pL",4)->FixParameter(3,6.16666e+00);
    CECA_ANA.GetFit("pL",4)->FixParameter(4,1.0);

    CECA_ANA.GetFit("pL",5)->FixParameter(0,2.25940e+03);
    CECA_ANA.GetFit("pL",5)->FixParameter(1,4.91680e-01);
    CECA_ANA.GetFit("pL",5)->FixParameter(2,2.20000e+00);
    CECA_ANA.GetFit("pL",5)->FixParameter(3,6.50003e+00);
    CECA_ANA.GetFit("pL",5)->FixParameter(4,1.0);
  }


  CECA_ANA.GoBabyGo(true);

  double TotChi2 = 0;
  double Chi2_pp = 0;
  double Chi2_pL = 0;
  printf("pp:");
  for(unsigned uMt=0; uMt<7; uMt++){
    double Chi2 = CECA_ANA.GetFit("pp",uMt)->GetChisquare();
    TotChi2 += Chi2;
    Chi2_pp += Chi2;
    printf(" %.0f",Chi2);
  }
  printf("\n");
  printf("pL:");
  for(unsigned uMt=0; uMt<6; uMt++){
    double Chi2 = CECA_ANA.GetFit("pL",uMt)->GetChisquare();
    TotChi2 += Chi2;
    Chi2_pL += Chi2;
    printf(" %.0f",Chi2);
  }
  printf("\n");

  printf("TotChi2 = %.0f\n",TotChi2);
  printf("Chi2_pp = %.0f\n",Chi2_pp);
  printf("Chi2_pL = %.0f\n",Chi2_pL);

  CECA_ANA.DumpCurrentCk(TString::Format("%s/CECA_Paper/TestSetUpAna/",GetFemtoOutputFolder())+FileBase,2);
  CECA_ANA.DumpCurrentCk(TString::Format("%s/CECA_Paper/TestSetUpAna/",GetFemtoOutputFolder())+FileBase,1);


}

//
//Cigar2
//Best_Chi2_tot = 1177
// Best_Dsp_tot = 4.533e-01
//  Best_Ht_tot = 4.333e+00
//  Best_Hz_tot = 9.500e+00
//The script terminated after: 10 min 0 s 439 ms
//confirmed with the new code

//finer binning:
//Best_Chi2_tot = 1070
// Best_Dsp_tot = 4.533e-01
//  Best_Ht_tot = 4.067e+00
//  Best_Hz_tot = 9.500e+00

//d=0.51 ht=3.9 hz=9.2
 //pp 139 82 52 58 68 93 116
 //pL 38 29 89 79 58 54
 //Tot 956 (956.48)
 //Best_Chi2_tot = 956
 //Best_Dsp_tot = 5.067e-01
  //Best_Ht_tot = 3.933e+00
//  Best_Hz_tot = 9.167e+00

//for chiral
//Best_Chi2_tot = 954
// Best_Dsp_tot = 5.081e-01
//  Best_Ht_tot = 3.913e+00
//  Best_Hz_tot = 9.112e+00

//for 9411602
//Best_Chi2_tot = 876
// Best_Dsp_tot = 4.600e-01
//  Best_Ht_tot = 3.478e+00
//  Best_Hz_tot = 9.491e+00

//Usmani
//d=0.51 ht=3.9 hz=9.1
// pp 130 78 51 55 61 89 125
// pL 49 37 92 93 64 59
// Tot 983 (982.59)



void ScanPs(){

  //for the Chi2
  const double FemtoLimit = 180;

  CecaAnalysis1 CECA_ANA("Reduced","Cigar2_ds24_hts36_hzs36",GetCernBoxDimi()+TString("/CatsFiles/"));

  CECA_ANA.SetUp_pp("AV18",0);

  //CECA_ANA.SetUp_pL("Chiral_Coupled_SPD",11600);
  //CECA_ANA.SetUp_pL("Chiral_Coupled_SPD",9411602);

////Wc=2279.0; Rc=0.3394; Sc=0.2614; f1=1.41; d1=2.53; tDev=0.003
///Wc=1940.7; Rc=0.4005; Sc=0.2565; f1=1.44; d1=2.54; tDev=0.036
  CECA_ANA.SetUp_pL("UsmaniFit",0);
  CECA_ANA.Kitty_pL->SetShortRangePotential(1,0,1,2279.0);
  CECA_ANA.Kitty_pL->SetShortRangePotential(1,0,2,0.3394);
  CECA_ANA.Kitty_pL->SetShortRangePotential(1,0,3,0.2614*1.015);
  double f0,d0;
  double f1,d1;
  TH1F* hPsFit;
  TF1* fPsFit;
  GetScattParameters(*CECA_ANA.Kitty_pL,f0,d0,hPsFit,fPsFit,2,false,false,0);
  GetScattParameters(*CECA_ANA.Kitty_pL,f1,d1,hPsFit,fPsFit,2,false,false,1);
  printf("S0: f0=%.3f d0=%.3f fm\n",f0,d0);
  printf("S1: f1=%.3f d1=%.3f fm\n",f1,d1);

  //TString SourceType = "Gauss";

  CECA_ANA.SetUp_pS0("Chiral",0);
  CECA_ANA.SetUp_pXim("pXim_HALQCDPaper2020",0);
  CECA_ANA.SetUp_pXi0("pXim_HALQCDPaper2020",0);
  CECA_ANA.SetUp_Decomposition(0,0);
  CECA_ANA.SetUp_Fits("SingleCeca",false);
  //CECA_ANA.SetUp_Fits(SourceType);

  unsigned Nbins_Dsp = CECA_ANA.GetSrc_pp()->GetNbins(1);
  double* BinCenter_Dsp = CECA_ANA.GetSrc_pp()->GetBinCenters(1);
  unsigned Nbins_Ht = CECA_ANA.GetSrc_pp()->GetNbins(2);
  double* BinCenter_Ht = CECA_ANA.GetSrc_pp()->GetBinCenters(2);
  unsigned Nbins_Hz = CECA_ANA.GetSrc_pp()->GetNbins(3);
  double* BinCenter_Hz = CECA_ANA.GetSrc_pp()->GetBinCenters(3);

  //original best sol
  //TH1F* hDsp = new TH1F("hDsp","hDsp",16,0.48,0.54);
  //TH1F* hHt = new TH1F("hHt","hHt",16,3.7,4.1);
  //TH1F* hHz = new TH1F("hHz","hHz",16,9.0,9.4);

  TH1F* hDsp;
  //if(SourceType=="Gauss"){
  //  hDsp = new TH1F("hDsp","hDsp",64,1,2);
  //}
  //else{
    hDsp = new TH1F("hDsp","hDsp",1,0.54,0.57);
  //}

  TH1F* hHt = new TH1F("hHt","hHt",1,3.6,4.8);
  TH1F* hHz = new TH1F("hHz","hHz",1,8.0,12.0);

  //double* Chi2_pp = new double [7];
  double* Best_Chi2_pp = new double [7];
  double* Best_Dsp_pp = new double [7];
  double* Best_Ht_pp = new double [7];
  double* Best_Hz_pp = new double [7];
  for(unsigned uMt=0; uMt<7; uMt++){
    Best_Chi2_pp[uMt] = 1e16;
    Best_Dsp_pp[uMt] = 0;
    Best_Ht_pp[uMt] = 0;
    Best_Hz_pp[uMt] = 0;
  }

  double* Best_Chi2_pL = new double [6];
  double* Best_Dsp_pL = new double [6];
  double* Best_Ht_pL = new double [6];
  double* Best_Hz_pL = new double [6];
  for(unsigned uMt=0; uMt<6; uMt++){
    Best_Chi2_pL[uMt] = 1e16;
    Best_Dsp_pL[uMt] = 0;
    Best_Ht_pL[uMt] = 0;
    Best_Hz_pL[uMt] = 0;
  }

  double Best_Chi2_tot = 1e16;
  double Best_Dsp_tot = 0;
  double Best_Ht_tot = 0;
  double Best_Hz_tot = 0;

  //[tot/pp/pL][tot/up to 180]
  double** Chi2 = new double* [3];
  double** Best_Chi2 = new double* [3];
  double** Best_Dsp = new double* [3];
  double** Best_Ht = new double* [3];
  double** Best_Hz = new double* [3];
  unsigned** Ndp = new unsigned* [3];

  for(unsigned uSce=0; uSce<3; uSce++){
    Chi2[uSce] = new double [2];
    Best_Chi2[uSce] = new double [2];
    Best_Dsp[uSce] = new double [2];
    Best_Ht[uSce] = new double [2];
    Best_Hz[uSce] = new double [2];
    Ndp[uSce] = new unsigned [2];

    for(unsigned uRng=0; uRng<2; uRng++){
      Best_Chi2[uSce][uRng] = 1e16;
      Best_Dsp[uSce][uRng] = 0;
      Best_Ht[uSce][uRng] = 0;
      Best_Hz[uSce][uRng] = 0;
      Ndp[uSce][uRng] = 0;
    }
  }

  //double* NDP_pp = new double [7];
  //double* Chi2_pL = new double [6];
  //double* NDP_pL = new double [6];

  //for(unsigned uDsp=0; uDsp<Nbins_Dsp; uDsp++){
  for(unsigned uDsp=0; uDsp<hDsp->GetNbinsX(); uDsp++){
    //if(uDsp%4) continue;
    //if(uDsp>=4) break;
    printf("uDsp = %u\n", uDsp);
    for(unsigned uHt=0; uHt<hHt->GetNbinsX(); uHt++){
      //if(uHt%4) continue;
      printf(" uHt = %u\n", uHt);
      for(unsigned uHz=0; uHz<hHz->GetNbinsX(); uHz++){
        /*
        if(SourceType=="Gauss"){
          for(unsigned uMt=0; uMt<7; uMt++){
            CECA_ANA.GetFit("pp",uMt)->FixParameter(0,-1);
            CECA_ANA.GetFit("pp",uMt)->FixParameter(1,hDsp->GetBinCenter(uDsp+1));
            CECA_ANA.GetFit("pp",uMt)->FixParameter(2,0);
            CECA_ANA.GetFit("pp",uMt)->FixParameter(3,0);
            CECA_ANA.GetFit("pp",uMt)->FixParameter(4,1.0);
          }

          for(unsigned uMt=0; uMt<6; uMt++){
            CECA_ANA.GetFit("pL",uMt)->FixParameter(0,-1);
            CECA_ANA.GetFit("pL",uMt)->FixParameter(1,hDsp->GetBinCenter(uDsp+1));
            CECA_ANA.GetFit("pL",uMt)->FixParameter(2,0);
            CECA_ANA.GetFit("pL",uMt)->FixParameter(3,0);
            CECA_ANA.GetFit("pL",uMt)->FixParameter(4,1.0);
          }
        }
        */
        //else{
          CECA_ANA.GetFit("pp",0)->FixParameter(0,1.10770e+03);
          CECA_ANA.GetFit("pp",1)->FixParameter(0,1.16830e+03);
          CECA_ANA.GetFit("pp",2)->FixParameter(0,1.22840e+03);
          CECA_ANA.GetFit("pp",3)->FixParameter(0,1.31560e+03);
          CECA_ANA.GetFit("pp",4)->FixParameter(0,1.46280e+03);
          CECA_ANA.GetFit("pp",5)->FixParameter(0,1.68720e+03);
          CECA_ANA.GetFit("pp",6)->FixParameter(0,2.21160e+03);
          for(unsigned uMt=0; uMt<7; uMt++){
            CECA_ANA.GetFit("pp",uMt)->FixParameter(1,hDsp->GetBinCenter(uDsp+1));
            CECA_ANA.GetFit("pp",uMt)->FixParameter(2,hHt->GetBinCenter(uHt+1));
            CECA_ANA.GetFit("pp",uMt)->FixParameter(3,hHz->GetBinCenter(uHz+1));
            CECA_ANA.GetFit("pp",uMt)->FixParameter(4,1.0);
          }

          CECA_ANA.GetFit("pL",0)->FixParameter(0,1.21240e+03);
          CECA_ANA.GetFit("pL",1)->FixParameter(0,1.28960e+03);
          CECA_ANA.GetFit("pL",2)->FixParameter(0,1.37600e+03);
          CECA_ANA.GetFit("pL",3)->FixParameter(0,1.54070e+03);
          CECA_ANA.GetFit("pL",4)->FixParameter(0,1.75600e+03);
          CECA_ANA.GetFit("pL",5)->FixParameter(0,2.25940e+03);
          for(unsigned uMt=0; uMt<6; uMt++){
            CECA_ANA.GetFit("pL",uMt)->FixParameter(1,hDsp->GetBinCenter(uDsp+1));
            CECA_ANA.GetFit("pL",uMt)->FixParameter(2,hHt->GetBinCenter(uHt+1));
            CECA_ANA.GetFit("pL",uMt)->FixParameter(3,hHz->GetBinCenter(uHz+1));
            CECA_ANA.GetFit("pL",uMt)->FixParameter(4,1.0);
          }
        //}

        CECA_ANA.GoBabyGo(false);

        //double TotChi2 = 0;
        //double Chi2_pp = 0;
        double Chi2fl_pp = 0;
        //double Chi2_pL = 0;
        //double Chi2fl_pL = 0;
        //unsigned NDP_pp = 0;
        //unsigned NDP_pL = 0;
        Ndp[0][0] = 0;
        Ndp[1][0] = 0;
        Ndp[2][0] = 0;
        Ndp[0][1] = 0;
        Ndp[1][1] = 0;
        Ndp[2][1] = 0;

        Chi2[0][0] = 0;
        Chi2[1][0] = 0;
        Chi2[2][0] = 0;
        Chi2[0][1] = 0;
        Chi2[1][1] = 0;
        Chi2[2][1] = 0;

        printf("d=%.2f ht=%.1f hz=%.1f\n",hDsp->GetBinCenter(uDsp+1),hHt->GetBinCenter(uHt+1),hHz->GetBinCenter(uHz+1));
        printf(" pp");
        for(unsigned uMt=0; uMt<7; uMt++){
          double Chi2_val = CECA_ANA.GetFit("pp",uMt)->GetChisquare();
          Ndp[0][0] += CECA_ANA.GetFit("pp",uMt)->GetNumberFitPoints();
          Ndp[1][0] += CECA_ANA.GetFit("pp",uMt)->GetNumberFitPoints();
          printf(" %.0f",Chi2_val);
          Chi2[0][0] += Chi2_val;
          Chi2[1][0] += Chi2_val;

          for(unsigned uMom=0; uMom<CECA_ANA.GetData("pp",uMt)->GetNbinsX(); uMom++){
            double MOM = CECA_ANA.GetData("pp",uMt)->GetBinCenter(uMom+1);
            if(MOM>FemtoLimit) break;
            double Chi2_fmt;
            if(CECA_ANA.GetData("pp",uMt)->GetBinContent(uMom+1)==0 || CECA_ANA.GetData("pp",uMt)->GetBinError(uMom+1)==0){
              Chi2_fmt = 0;
            }
            else{
              Chi2_fmt = CECA_ANA.GetData("pp",uMt)->GetBinContent(uMom+1)-CECA_ANA.GetFit("pp",uMt)->Eval(MOM);
              Chi2_fmt /= CECA_ANA.GetData("pp",uMt)->GetBinError(uMom+1);
              Chi2_fmt = Chi2_fmt*Chi2_fmt;
            }
            Chi2[0][1] += Chi2_fmt;
            Chi2[1][1] += Chi2_fmt;
            Ndp[0][1]++;
            Ndp[1][1]++;
          }

        }
        printf("\n");

        printf(" pL");
        for(unsigned uMt=0; uMt<6; uMt++){
          double Chi2_val = CECA_ANA.GetFit("pL",uMt)->GetChisquare();
          Ndp[0][0] += CECA_ANA.GetFit("pL",uMt)->GetNumberFitPoints();
          Ndp[2][0] += CECA_ANA.GetFit("pL",uMt)->GetNumberFitPoints();
          printf(" %.0f",Chi2_val);
          Chi2[0][0] += Chi2_val;
          Chi2[2][0] += Chi2_val;

          for(unsigned uMom=0; uMom<CECA_ANA.GetData("pL",uMt)->GetNbinsX(); uMom++){
            double MOM = CECA_ANA.GetData("pL",uMt)->GetBinCenter(uMom+1);
            if(MOM>FemtoLimit) break;
            if(MOM>260&&MOM<300){
              continue;
            }
            double Chi2_fmt;
            if(CECA_ANA.GetData("pL",uMt)->GetBinContent(uMom+1)==0 || CECA_ANA.GetData("pL",uMt)->GetBinError(uMom+1)==0){
              Chi2_fmt = 0;
            }
            else{
              Chi2_fmt = CECA_ANA.GetData("pL",uMt)->GetBinContent(uMom+1)-CECA_ANA.GetFit("pL",uMt)->Eval(MOM);
              Chi2_fmt /= CECA_ANA.GetData("pL",uMt)->GetBinError(uMom+1);
              Chi2_fmt = Chi2_fmt*Chi2_fmt;
            }
            Chi2[0][1] += Chi2_fmt;
            Chi2[2][1] += Chi2_fmt;
            Ndp[0][1]++;
            Ndp[2][1]++;
          }
        }

        //total
        if(Chi2[0][0]<Best_Chi2[0][0]){
          Best_Chi2[0][0] = Chi2[0][0];
          //all fits are the same, so it does not matter which mT bin or pp/pL we take
          Best_Dsp[0][0] =  CECA_ANA.GetFit("pp",0)->GetParameter(1);
          Best_Ht[0][0] =  CECA_ANA.GetFit("pp",0)->GetParameter(2);
          Best_Hz[0][0] =  CECA_ANA.GetFit("pp",0)->GetParameter(3);
        }

        //pp
        if(Chi2[1][0]<Best_Chi2[1][0]){
          Best_Chi2[1][0] = Chi2[1][0];
          Best_Dsp[1][0] =  CECA_ANA.GetFit("pp",0)->GetParameter(1);
          Best_Ht[1][0] =  CECA_ANA.GetFit("pp",0)->GetParameter(2);
          Best_Hz[1][0] =  CECA_ANA.GetFit("pp",0)->GetParameter(3);
        }

        //pL
        if(Chi2[2][0]<Best_Chi2[2][0]){
          Best_Chi2[2][0] = Chi2[2][0];
          Best_Dsp[2][0] =  CECA_ANA.GetFit("pp",0)->GetParameter(1);
          Best_Ht[2][0] =  CECA_ANA.GetFit("pp",0)->GetParameter(2);
          Best_Hz[2][0] =  CECA_ANA.GetFit("pp",0)->GetParameter(3);
        }

        //total in femto region
        if(Chi2[0][1]<Best_Chi2[0][1]){
          Best_Chi2[0][1] = Chi2[0][1];
          Best_Dsp[0][1] =  CECA_ANA.GetFit("pp",0)->GetParameter(1);
          Best_Ht[0][1] =  CECA_ANA.GetFit("pp",0)->GetParameter(2);
          Best_Hz[0][1] =  CECA_ANA.GetFit("pp",0)->GetParameter(3);
        }

        //pp in femto region
        if(Chi2[1][1]<Best_Chi2[1][1]){
          Best_Chi2[1][1] = Chi2[1][1];
          Best_Dsp[1][1] =  CECA_ANA.GetFit("pp",0)->GetParameter(1);
          Best_Ht[1][1] =  CECA_ANA.GetFit("pp",0)->GetParameter(2);
          Best_Hz[1][1] =  CECA_ANA.GetFit("pp",0)->GetParameter(3);
        }

        //pL in femto region
        if(Chi2[2][1]<Best_Chi2[2][1]){
          Best_Chi2[2][1] = Chi2[2][1];
          Best_Dsp[2][1] =  CECA_ANA.GetFit("pp",0)->GetParameter(1);
          Best_Ht[2][1] =  CECA_ANA.GetFit("pp",0)->GetParameter(2);
          Best_Hz[2][1] =  CECA_ANA.GetFit("pp",0)->GetParameter(3);
        }

        printf("\n");
        printf(" Tot %.0f (%.2f)\n",Chi2[0][0],Best_Chi2[0][0]);
        usleep(100e3);
      }//Hz
    }//Ht
  }//Dsp


  printf("\nTHE BIG SUMMARY:\n");
  printf("\n--- for the total fit ---\n");
  printf("Best_Chi2_tot = %.0f / %u\n",Best_Chi2[0][0],Ndp[0][0]);
  printf(" Best_Dsp_tot = %.3e\n",Best_Dsp[0][0]);
  printf("  Best_Ht_tot = %.3e\n",Best_Ht[0][0]);
  printf("  Best_Hz_tot = %.3e\n",Best_Hz[0][0]);

  printf("Best_Chi2_pp = %.0f / %u\n",Best_Chi2[1][0],Ndp[1][0]);
  printf(" Best_Dsp_pp = %.3e\n",Best_Dsp[1][0]);
  printf("  Best_Ht_pp = %.3e\n",Best_Ht[1][0]);
  printf("  Best_Hz_pp = %.3e\n",Best_Hz[1][0]);

  printf("Best_Chi2_pL = %.0f / %u\n",Best_Chi2[2][0],Ndp[2][0]);
  printf(" Best_Dsp_pL = %.3e\n",Best_Dsp[2][0]);
  printf("  Best_Ht_pL = %.3e\n",Best_Ht[2][0]);
  printf("  Best_Hz_pL = %.3e\n",Best_Hz[2][0]);

  printf("\n--- for k*<%.0f ---\n",FemtoLimit);
  printf("Best_Chi2_tot = %.0f / %u\n",Best_Chi2[0][1],Ndp[0][1]);
  printf(" Best_Dsp_tot = %.3e\n",Best_Dsp[0][1]);
  printf("  Best_Ht_tot = %.3e\n",Best_Ht[0][1]);
  printf("  Best_Hz_tot = %.3e\n",Best_Hz[0][1]);

  printf("Best_Chi2_pp = %.0f / %u\n",Best_Chi2[1][1],Ndp[1][1]);
  printf(" Best_Dsp_pp = %.3e\n",Best_Dsp[1][1]);
  printf("  Best_Ht_pp = %.3e\n",Best_Ht[1][1]);
  printf("  Best_Hz_pp = %.3e\n",Best_Hz[1][1]);

  printf("Best_Chi2_pL = %.0f / %u\n",Best_Chi2[2][1],Ndp[2][1]);
  printf(" Best_Dsp_pL = %.3e\n",Best_Dsp[2][1]);
  printf("  Best_Ht_pL = %.3e\n",Best_Ht[2][1]);
  printf("  Best_Hz_pL = %.3e\n",Best_Hz[2][1]);

  delete [] BinCenter_Dsp;
  delete [] BinCenter_Ht;
  delete [] BinCenter_Hz;
  //delete [] Chi2_pp;
  //delete [] Chi2_pL;
}

















void ScanPsUsmani(TString AnaType, TString SourceDescription, TString cern_box, TString out_folder,
                  double min_d, double max_d,
                  double min_ht, double max_ht,
                  double min_hz, double max_hz,
                  double min_wc, double max_wc,
                  double min_rc, double max_rc,
                  double min_sc, double max_sc,
                  int VAR_L,
                  double Minutes, int SEED
){

  TRandom3 rangen(SEED);

  const long long JobTimeLimit = (long long)(Minutes*60.);
  DLM_Timer JobTime;

  //for the Chi2
  const double FemtoLimit = 180;
  CecaAnalysis1 CECA_ANA(AnaType,SourceDescription,cern_box+TString("/CatsFiles/"));

  CECA_ANA.SetUp_pp("AV18",0);
////Wc=2279.0; Rc=0.3394; Sc=0.2614; f1=1.41; d1=2.53; tDev=0.003
  CECA_ANA.SetUp_pL("UsmaniFit",0);
  CECA_ANA.Kitty_pL->SetShortRangePotential(1,0,1,2279.0);
  CECA_ANA.Kitty_pL->SetShortRangePotential(1,0,2,0.3394);
  CECA_ANA.Kitty_pL->SetShortRangePotential(1,0,3,0.2614);
  CECA_ANA.Kitty_pL->SetNotifications(CATS::nWarning);

  TNtuple* ntResult;
  //16
  ntResult = new TNtuple("ntResult", "ntResult","VAR_L:d:ht:hz:wc:rc:sc:f1:d1:chi2_pp:ndp_pp:chi2_pL:ndp_pL:chi2_pp_fmt:ndp_pp_fmt:chi2_pL_fmt:ndp_pL_fmt");
  Float_t* Entries = new Float_t [17];

  double f0,d0;
  double f1,d1;
  TH1F* hPsFit;
  TF1* fPsFit;

  CECA_ANA.SetUp_pS0("Chiral",0);
  CECA_ANA.SetUp_pXim("pXim_HALQCDPaper2020",0);
  CECA_ANA.SetUp_pXi0("pXim_HALQCDPaper2020",0);
  CECA_ANA.SetUp_Decomposition(0,VAR_L);
  CECA_ANA.SetUp_Fits("SingleCeca",false);

  double* Best_Chi2_pp = new double [7];
  double* Best_Dsp_pp = new double [7];
  double* Best_Ht_pp = new double [7];
  double* Best_Hz_pp = new double [7];
  for(unsigned uMt=0; uMt<7; uMt++){
    Best_Chi2_pp[uMt] = 1e16;
    Best_Dsp_pp[uMt] = 0;
    Best_Ht_pp[uMt] = 0;
    Best_Hz_pp[uMt] = 0;
  }

  double* Best_Chi2_pL = new double [6];
  double* Best_Dsp_pL = new double [6];
  double* Best_Ht_pL = new double [6];
  double* Best_Hz_pL = new double [6];
  for(unsigned uMt=0; uMt<6; uMt++){
    Best_Chi2_pL[uMt] = 1e16;
    Best_Dsp_pL[uMt] = 0;
    Best_Ht_pL[uMt] = 0;
    Best_Hz_pL[uMt] = 0;
  }

  double Best_Chi2_tot = 1e16;
  double Best_Dsp_tot = 0;
  double Best_Ht_tot = 0;
  double Best_Hz_tot = 0;

  //[tot/pp/pL][tot/up to 180]
  double** Chi2 = new double* [3];
  double** Best_Chi2 = new double* [3];
  double** Best_Dsp = new double* [3];
  double** Best_Ht = new double* [3];
  double** Best_Hz = new double* [3];
  unsigned** Ndp = new unsigned* [3];

  for(unsigned uSce=0; uSce<3; uSce++){
    Chi2[uSce] = new double [2];
    Best_Chi2[uSce] = new double [2];
    Best_Dsp[uSce] = new double [2];
    Best_Ht[uSce] = new double [2];
    Best_Hz[uSce] = new double [2];
    Ndp[uSce] = new unsigned [2];

    for(unsigned uRng=0; uRng<2; uRng++){
      Best_Chi2[uSce][uRng] = 1e16;
      Best_Dsp[uSce][uRng] = 0;
      Best_Ht[uSce][uRng] = 0;
      Best_Hz[uSce][uRng] = 0;
      Ndp[uSce][uRng] = 0;
    }
  }

  long long TotalTime = 0;
  int Iter = 0;
  while(TotalTime<JobTimeLimit){
    double ran_d = rangen.Uniform(min_d,max_d);
    double ran_ht = rangen.Uniform(min_ht,max_ht);
    double ran_hz = rangen.Uniform(min_hz,max_hz);
    double ran_wc = rangen.Uniform(min_wc,max_wc);
    double ran_rc = rangen.Uniform(min_rc,max_rc);
    double ran_sc = rangen.Uniform(min_sc,max_sc);

    CECA_ANA.Kitty_pL->SetShortRangePotential(1,0,1,ran_wc);
    CECA_ANA.Kitty_pL->SetShortRangePotential(1,0,2,ran_rc);
    CECA_ANA.Kitty_pL->SetShortRangePotential(1,0,3,ran_sc);
    CECA_ANA.Kitty_pL->KillTheCat();

    GetScattParameters(*CECA_ANA.Kitty_pL,f1,d1,hPsFit,fPsFit,2,false,false,1);
    delete hPsFit;
    delete fPsFit;

    CECA_ANA.GetFit("pp",0)->FixParameter(0,1.10770e+03);
    CECA_ANA.GetFit("pp",1)->FixParameter(0,1.16830e+03);
    CECA_ANA.GetFit("pp",2)->FixParameter(0,1.22840e+03);
    CECA_ANA.GetFit("pp",3)->FixParameter(0,1.31560e+03);
    CECA_ANA.GetFit("pp",4)->FixParameter(0,1.46280e+03);
    CECA_ANA.GetFit("pp",5)->FixParameter(0,1.68720e+03);
    CECA_ANA.GetFit("pp",6)->FixParameter(0,2.21160e+03);
    for(unsigned uMt=0; uMt<7; uMt++){
      CECA_ANA.GetFit("pp",uMt)->FixParameter(1,ran_d);
      CECA_ANA.GetFit("pp",uMt)->FixParameter(2,ran_ht);
      CECA_ANA.GetFit("pp",uMt)->FixParameter(3,ran_hz);
      CECA_ANA.GetFit("pp",uMt)->FixParameter(4,1.0);
    }

    CECA_ANA.GetFit("pL",0)->FixParameter(0,1.21240e+03);
    CECA_ANA.GetFit("pL",1)->FixParameter(0,1.28960e+03);
    CECA_ANA.GetFit("pL",2)->FixParameter(0,1.37600e+03);
    CECA_ANA.GetFit("pL",3)->FixParameter(0,1.54070e+03);
    CECA_ANA.GetFit("pL",4)->FixParameter(0,1.75600e+03);
    CECA_ANA.GetFit("pL",5)->FixParameter(0,2.25940e+03);
    for(unsigned uMt=0; uMt<6; uMt++){
      CECA_ANA.GetFit("pL",uMt)->FixParameter(1,ran_d);
      CECA_ANA.GetFit("pL",uMt)->FixParameter(2,ran_ht);
      CECA_ANA.GetFit("pL",uMt)->FixParameter(3,ran_hz);
      CECA_ANA.GetFit("pL",uMt)->FixParameter(4,1.0);
    }

    CECA_ANA.GoBabyGo(false);

    for(unsigned uSce=0; uSce<3; uSce++){for(unsigned uRng=0; uRng<2; uRng++){
        Ndp[uSce][uRng] = 0;
        Chi2[uSce][uRng] = 0;
    }}

    for(unsigned uMt=0; uMt<7; uMt++){
      double Chi2_val = CECA_ANA.GetFit("pp",uMt)->GetChisquare();
      Ndp[0][0] += CECA_ANA.GetFit("pp",uMt)->GetNumberFitPoints();
      Ndp[1][0] += CECA_ANA.GetFit("pp",uMt)->GetNumberFitPoints();
      Chi2[0][0] += Chi2_val;
      Chi2[1][0] += Chi2_val;

      for(unsigned uMom=0; uMom<CECA_ANA.GetData("pp",uMt)->GetNbinsX(); uMom++){
        double MOM = CECA_ANA.GetData("pp",uMt)->GetBinCenter(uMom+1);
        if(MOM>FemtoLimit) break;
        double Chi2_fmt;
        if(CECA_ANA.GetData("pp",uMt)->GetBinContent(uMom+1)==0 || CECA_ANA.GetData("pp",uMt)->GetBinError(uMom+1)==0){
          Chi2_fmt = 0;
        }
        else{
          Chi2_fmt = CECA_ANA.GetData("pp",uMt)->GetBinContent(uMom+1)-CECA_ANA.GetFit("pp",uMt)->Eval(MOM);
          Chi2_fmt /= CECA_ANA.GetData("pp",uMt)->GetBinError(uMom+1);
          Chi2_fmt = Chi2_fmt*Chi2_fmt;
        }
        Chi2[0][1] += Chi2_fmt;
        Chi2[1][1] += Chi2_fmt;
        Ndp[0][1]++;
        Ndp[1][1]++;
      }
    }//uMt pp

    for(unsigned uMt=0; uMt<6; uMt++){
      double Chi2_val = CECA_ANA.GetFit("pL",uMt)->GetChisquare();
      Ndp[0][0] += CECA_ANA.GetFit("pL",uMt)->GetNumberFitPoints();
      Ndp[2][0] += CECA_ANA.GetFit("pL",uMt)->GetNumberFitPoints();
      Chi2[0][0] += Chi2_val;
      Chi2[2][0] += Chi2_val;

      for(unsigned uMom=0; uMom<CECA_ANA.GetData("pL",uMt)->GetNbinsX(); uMom++){
        double MOM = CECA_ANA.GetData("pL",uMt)->GetBinCenter(uMom+1);
        if(MOM>FemtoLimit) break;
        if(MOM>260&&MOM<300){
          continue;
        }
        double Chi2_fmt;
        if(CECA_ANA.GetData("pL",uMt)->GetBinContent(uMom+1)==0 || CECA_ANA.GetData("pL",uMt)->GetBinError(uMom+1)==0){
          Chi2_fmt = 0;
        }
        else{
          Chi2_fmt = CECA_ANA.GetData("pL",uMt)->GetBinContent(uMom+1)-CECA_ANA.GetFit("pL",uMt)->Eval(MOM);
          Chi2_fmt /= CECA_ANA.GetData("pL",uMt)->GetBinError(uMom+1);
          Chi2_fmt = Chi2_fmt*Chi2_fmt;
        }
        Chi2[0][1] += Chi2_fmt;
        Chi2[2][1] += Chi2_fmt;
        Ndp[0][1]++;
        Ndp[2][1]++;
      }
    }//uMt pL

    for(unsigned uSce=0; uSce<3; uSce++){for(unsigned uRng=0; uRng<2; uRng++){
      if(Chi2[uSce][uRng]<Best_Chi2[uSce][uRng]){
        //all fits are the same, so it does not matter which mT bin or pp/pL we take
        Best_Chi2[uSce][uRng] = Chi2[uSce][uRng];
        Best_Dsp[uSce][uRng] =  CECA_ANA.GetFit("pp",0)->GetParameter(1);
        Best_Ht[uSce][uRng] =  CECA_ANA.GetFit("pp",0)->GetParameter(2);
        Best_Hz[uSce][uRng] =  CECA_ANA.GetFit("pp",0)->GetParameter(3);
      }

    }}

    //d:ht:hz:wc:rc:sc:f1:d1:chi2_pp:ndp_pp:chi2_pL:ndp_pL:chi2_pp_fmt:ndp_pp_fmt:chi2_pL_fmt:ndp_pL_fmt
    Entries[0] = VAR_L;
    Entries[1] = ran_d;
    Entries[2] = ran_ht;
    Entries[3] = ran_hz;
    Entries[4] = ran_wc;
    Entries[5] = ran_rc;
    Entries[6] = ran_sc;
    Entries[7] = f1;
    Entries[8] = d1;
    Entries[9] = Chi2[1][0];
    Entries[10] = Ndp[1][0];
    Entries[11] = Chi2[2][0];
    Entries[12] = Ndp[2][0];
    Entries[13] = Chi2[1][1];
    Entries[14] = Ndp[1][1];
    Entries[15] = Chi2[2][1];
    Entries[16] = Ndp[2][1];

    ntResult->Fill(Entries);

    TotalTime = JobTime.Stop()/1000000.;
    Iter++;
    if(Iter%10==0) printf("iIter = %i\n",Iter);
  }

  TFile fOutput(out_folder+TString::Format("/US_%s_%s_s%i.root",AnaType.Data(),SourceDescription.Data(),SEED),"recreate");
  ntResult->Write();

  delete [] Best_Chi2_pp;
  delete [] Best_Dsp_pp;
  delete [] Best_Ht_pp;
  delete [] Best_Hz_pp;

  delete [] Best_Chi2_pL;
  delete [] Best_Dsp_pL;
  delete [] Best_Ht_pL;
  delete [] Best_Hz_pL;

  for(unsigned uSce=0; uSce<3; uSce++){
    delete [] Ndp[uSce];
    delete [] Chi2[uSce];
    delete [] Best_Chi2[uSce];
    delete [] Best_Dsp[uSce];
    delete [] Best_Ht[uSce];
    delete [] Best_Hz[uSce];
  }
  delete [] Ndp;
  delete [] Chi2;
  delete [] Best_Chi2;
  delete [] Best_Dsp;
  delete [] Best_Ht;
  delete [] Best_Hz;

  delete [] Entries;
  delete ntResult;
}











//Wc=2288; Rc=0.337; Sc=0.261; f1=1.49; d1=2.46; tDev=0.105
//Wc=2288; Rc=0.346; Sc=0.258; f1=1.40; d1=2.65; tDev=0.124
//Wc=1813; Rc=0.421; Sc=0.258; f1=1.31; d1=2.55; tDev=0.100
void UsmaniFirstLook(int SEED){
  TRandom3 rangen(SEED);

  double Wc_Min = rangen.Uniform(2137-400,2137+600);
  double Wc_Max = Wc_Min;
  unsigned Wc_Steps = 1;

  double Rc_Min = rangen.Uniform(0.3,0.8);
  double Rc_Max = Rc_Min*1.1;
  unsigned Rc_Steps = 8;

  double Sc_Min = rangen.Uniform(0.1,0.3);
  double Sc_Max = Sc_Min*1.1;
  unsigned Sc_Steps = 8;

  TH2F* hSL = new TH2F("hSL","hSL",Wc_Steps,Wc_Min,Wc_Max,Rc_Steps,Rc_Min,Rc_Max);
  TH2F* hER = new TH2F("hER","hER",Wc_Steps,Wc_Min,Wc_Max,Rc_Steps,Rc_Min,Rc_Max);

  TH1F* hSc = new TH1F("hSc","hSc",Sc_Steps,Sc_Min,Sc_Max);

  const double fGoal = 1.41;
  const double dGoal = 2.53;
  double fDev, dDev, tDev;

  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());
  CATS Kitty;
  Kitty.SetMomBins(30,0,120);
  Kitty.SetEpsilonConv(4e-8);
  Kitty.SetEpsilonProp(4e-8);
  Kitty.SetNotifications(CATS::nWarning);
  AnalysisObject.SetUpCats_pL(Kitty,"UsmaniFit","Gauss",0,0);//NLO_Coupled_S

  for(unsigned uWc=0; uWc<Wc_Steps; uWc++){
    for(unsigned uRc=0; uRc<Rc_Steps; uRc++){
      for(unsigned uSc=0; uSc<Sc_Steps; uSc++){
        double Wc_Val = hSL->GetXaxis()->GetBinCenter(uWc+1);
        double Rc_Val = hSL->GetYaxis()->GetBinCenter(uRc+1);
        double Sc_Val = hSc->GetXaxis()->GetBinCenter(uSc+1);
        Kitty.SetShortRangePotential(1,0,1,Wc_Val);
        Kitty.SetShortRangePotential(1,0,2,Rc_Val);
        Kitty.SetShortRangePotential(1,0,3,Sc_Val);
        Kitty.KillTheCat();
        double f1,d1;
        TH1F* hPsFit;
        TF1* fPsFit;
        GetScattParameters(Kitty,f1,d1,hPsFit,fPsFit,3,false,false,1);
        fDev = f1-fGoal;
        dDev = d1-dGoal;
        tDev = sqrt(fDev*fDev+dDev*dDev);
        delete hPsFit;
        delete fPsFit;

        //tDev /= 16;
        if(tDev<1./1.){
          if(tDev>1./2.) printf("\033[41m ");
          else if(tDev>1./4.) printf("\033[43m ");
          else if(tDev>1./8.) printf("\033[42m ");
          else if(tDev>1./16.) printf("\033[46m ");
          else printf("\033[44m ");
          printf("Wc=%.0f; Rc=%.3f; Sc=%.3f; f1=%.2f; d1=%.2f; tDev=%.3f\033[0m\n",Wc_Val,Rc_Val,Sc_Val,f1,d1,tDev);
        }

        //hSL->SetBinContent(uWc+1,uRc+1,f1);
        //hER->SetBinContent(uWc+1,uRc+1,d1);
      }
    }
  }

  TFile fOutput(TString::Format("%s/CECA_Paper/UsmaniFirstLook/fOutput.root",GetFemtoOutputFolder()),"recreate");
  hSL->Write();
  hER->Write();

  delete hSL;
  delete hER;
  delete hSc;
}

//first look:
//Wc=2288; Rc=0.337; Sc=0.261; f1=1.49; d1=2.46; tDev=0.105 !!
//Wc=2288; Rc=0.346; Sc=0.258; f1=1.40; d1=2.65; tDev=0.124
//Wc=1813; Rc=0.421; Sc=0.258; f1=1.31; d1=2.55; tDev=0.100 !!

//Wc=1813; Rc=0.414; Sc=0.261; f1=1.33; d1=2.40; tDev=0.158
//Wc=1813; Rc=0.430; Sc=0.252; f1=1.55; d1=2.55; tDev=0.145
//Wc=2350; Rc=0.325; Sc=0.263; f1=1.58; d1=2.37; tDev=0.233 !
//Wc=2350; Rc=0.334; Sc=0.260; f1=1.36; d1=2.65; tDev=0.135 !!
//Wc=1960; Rc=0.403; Sc=0.254; f1=1.46; d1=2.63; tDev=0.115
void UsmaniSecondLook(int Wc_Steps, int SEED){
  TRandom3 rangen(SEED);

  const double fGoal = 1.41;
  const double dGoal = 2.53;
  double fDev, dDev, tDev;

  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());
  CATS Kitty;
  Kitty.SetMomBins(30,0,120);
  Kitty.SetEpsilonConv(4e-8);
  Kitty.SetEpsilonProp(4e-8);
  Kitty.SetNotifications(CATS::nWarning);
  AnalysisObject.SetUpCats_pL(Kitty,"UsmaniFit","Gauss",0,0);//NLO_Coupled_S


  for(unsigned uWc=0; uWc<Wc_Steps; uWc++){
    double Wc_Val = rangen.Uniform(2137-400,2137+600);
    double Rc_Min = rangen.Uniform(0.3,0.8);
    double Rc_Step = Rc_Min*0.01;
    for(unsigned uRc=0; uRc<8; uRc++){
      double Sc_Min = rangen.Uniform(0.1,0.3);
      double Sc_Step = Sc_Min*0.01;
      double Rc_Val = Rc_Min + double(uRc)*Rc_Step;
      for(unsigned uSc=0; uSc<8; uSc++){
        double Sc_Val = Sc_Min + double(uSc)*Sc_Step;

        Kitty.SetShortRangePotential(1,0,1,Wc_Val);
        Kitty.SetShortRangePotential(1,0,2,Rc_Val);
        Kitty.SetShortRangePotential(1,0,3,Sc_Val);
        Kitty.KillTheCat();
        double f1,d1;
        TH1F* hPsFit;
        TF1* fPsFit;
        GetScattParameters(Kitty,f1,d1,hPsFit,fPsFit,3,false,false,1);
        fDev = f1-fGoal;
        dDev = d1-dGoal;
        tDev = sqrt(fDev*fDev+dDev*dDev);
        delete hPsFit;
        delete fPsFit;

        //tDev /= 16;
        if(tDev<1./1.){
          if(tDev>1./2.) printf("\033[41m ");
          else if(tDev>1./4.) printf("\033[43m ");
          else if(tDev>1./8.) printf("\033[42m ");
          else if(tDev>1./16.) printf("\033[46m ");
          else printf("\033[44m ");
          printf("Wc=%.0f; Rc=%.3f; Sc=%.3f; f1=%.2f; d1=%.2f; tDev=%.3f\033[0m\n",Wc_Val,Rc_Val,Sc_Val,f1,d1,tDev);
        }

      }
    }
  }
}


//Color code by 1/4
//Wc=2277; Rc=0.339; Sc=0.262; f1=1.41; d1=2.52; tDev=0.014
//Wc=2300; Rc=0.336; Sc=0.261; f1=1.43; d1=2.51; tDev=0.028
//Wc=2282; Rc=0.340; Sc=0.261; f1=1.41; d1=2.54; tDev=0.010
//Wc=2291; Rc=0.337; Sc=0.262; f1=1.40; d1=2.53; tDev=0.009
//Wc=2278.4; Rc=0.3397; Sc=0.2613; f1=1.41; d1=2.53; tDev=0.004
//Wc=2285.1; Rc=0.3384; Sc=0.2614; f1=1.41; d1=2.53; tDev=0.003
//Wc=2288.3; Rc=0.3380; Sc=0.2614; f1=1.41; d1=2.53; tDev=0.004
//Wc=2279.0; Rc=0.3394; Sc=0.2614; f1=1.41; d1=2.53; tDev=0.003 //use this
//Wc=1940.7; Rc=0.4005; Sc=0.2565; f1=1.44; d1=2.54; tDev=0.036 //or this, as you decide, both QAd with UsmaniMatch to NLO

//try: Wc=1415.99; Rc=0.339256; Sc=0.254946;

void UsmaniFineCheck(){
  const double Wc_Start = 1960;
  const double Rc_Start = 0.403;
  const double Sc_Start = 0.254;
  const double Fineness = 1.01;

  const unsigned IterPerPar = 6;

  TRandom3 rangen(13);

  const double fGoal = 1.41;
  const double dGoal = 2.53;
  double fDev, dDev, tDev;

  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());
  CATS Kitty;
  Kitty.SetMomBins(30,0,120);
  Kitty.SetEpsilonConv(4e-8);
  Kitty.SetEpsilonProp(4e-8);
  Kitty.SetNotifications(CATS::nWarning);
  AnalysisObject.SetUpCats_pL(Kitty,"UsmaniFit","Gauss",0,0);//NLO_Coupled_S


  for(unsigned uWc=0; uWc<IterPerPar; uWc++){
    double Wc_Val = rangen.Uniform(Wc_Start/Fineness,Wc_Start*Fineness);
    for(unsigned uRc=0; uRc<IterPerPar; uRc++){
      double Rc_Val = rangen.Uniform(Rc_Start/Fineness,Rc_Start*Fineness);
      for(unsigned uSc=0; uSc<IterPerPar; uSc++){
        double Sc_Val = rangen.Uniform(Sc_Start/Fineness,Sc_Start*Fineness);

        Kitty.SetShortRangePotential(1,0,1,Wc_Val);
        Kitty.SetShortRangePotential(1,0,2,Rc_Val);
        Kitty.SetShortRangePotential(1,0,3,Sc_Val);
        Kitty.KillTheCat();
        double f1,d1;
        TH1F* hPsFit;
        TF1* fPsFit;
        GetScattParameters(Kitty,f1,d1,hPsFit,fPsFit,3,false,false,1);
        fDev = f1-fGoal;
        dDev = d1-dGoal;
        tDev = sqrt(fDev*fDev+dDev*dDev);
        delete hPsFit;
        delete fPsFit;

        //tDev /= 16;
        if(tDev<1./4.){
          if(tDev>1./8.) printf("\033[41m ");
          else if(tDev>1./16.) printf("\033[43m ");
          else if(tDev>1./32.) printf("\033[42m ");
          else if(tDev>1./64.) printf("\033[46m ");
          else printf("\033[44m ");
          printf("Wc=%.1f; Rc=%.4f; Sc=%.4f; f1=%.2f; d1=%.2f; tDev=%.3f\033[0m\n",Wc_Val,Rc_Val,Sc_Val,f1,d1,tDev);
        }

      }
    }
  }
}


void Usmani_vs_NLO19(){

}

int CECA_PAPER(int argc, char *argv[]){
  printf("CECA_PAPER\n\n");
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
  //printf("Lets go\n");
  //return dlmhst_ceca_kdpfit_wc0(5,myinput);
  //Test_src_read();
  //LnGammaTest();
  //CreateDlmHistoCigar2("pp",TString::Format("%s/CECA_Paper/Output_260223/",GetFemtoOutputFolder()),"Cigar2");
  //CreateDlmHistoCigar2("pL",TString::Format("%s/CECA_Paper/Output_260223/",GetFemtoOutputFolder()),"Cigar2");
  //TestReadCigar2();
  //TestReadJaime1();
  //TestCigar2_Ck();
  //TestJaime1_Ck();
  //PNG_ANIM_SOURCE();
  //TestSetUpAna();
  //ScanPs();
  //UsmaniSecondLook(atoi(argv[1]),atoi(argv[2]));
  //UsmaniFineCheck();
  //Usmani_vs_NLO19();

////Wc=2279.0; Rc=0.3394; Sc=0.2614; f1=1.41; d1=2.53; tDev=0.003


/*
ScanPsUsmani(
                  "Reduced","Cigar2_ds24_hts36_hzs36",
                  TString(GetCernBoxDimi()), TString::Format("%s/CECA_Paper/ScanPsUsmani/",GetFemtoOutputFolder()),
                  0.3, 0.7,//d
                  0.5,4.8,//ht
                  4.5,12.0,//hz
                  2279.0*1.030, 2279.0/1.020,
                  0.3394*1.020, 0.3394/1.015,
                  0.2614*1.015, 0.2614/1.010,
                  atof(argv[1]), atoi(argv[2]));
*/
/*
  ScanPsUsmani(
                    "Reduced","Jaime1_ds24_hts36_hzs36",
                    TString(GetCernBoxDimi()), TString::Format("%s/CECA_Paper/ScanPsUsmani/",GetFemtoOutputFolder()),
                    0.2/1.2, 0.2*1.2,//d
                    2.3/1.2,2.3*1.2,//ht
                    3.55/1.2,3.55*1.2,//hz
                    2279.0*1.030, 2279.0/1.020,
                    0.3394*1.020, 0.3394/1.015,
                    0.2614*1.015, 0.2614/1.010,
                    atof(argv[1]), atoi(argv[2]));
*/
/*
ScanPsUsmani(
                  "Reduced","Jaime1_ds24_hts36_hzs36",
                  TString(GetCernBoxDimi()), TString::Format("%s/CECA_Paper/ScanPsUsmani/",GetFemtoOutputFolder()),
                  0,0.8,//d
                  0,4.8,//ht
                  0,6.0,//hz
                  2279.0*1.030, 2279.0/1.020,
                  0.3394*1.020, 0.3394/1.015,
                  0.2614*1.015, 0.2614/1.010,
                  atof(argv[1]), atoi(argv[2]));
*/
/*
ScanPsUsmani(
                  "Reduced","Jaime1_ds24_hts36_hzs36",
                  TString(GetCernBoxDimi()), TString::Format("%s/CECA_Paper/ScanPsUsmani/",GetFemtoOutputFolder()),
                  0.274-0.050, 0.274+0.050,//d
                  3.63-0.25, 3.63+0.25,//ht
                  2.68-0.2, 2.68+0.2,//hz
                  2311-40, 2311+40,
                  0.3453-0.006, 0.3453+0.006,
                  0.2595-0.001, 0.2595+0.001,
                  atof(argv[1]), atoi(argv[2]));
*/

//USE THIS ONE FOR THE PAPER
ScanPsUsmani(
                  "Reduced","Jaime1_ds24_hts36_hzs36",
                  TString(GetCernBoxDimi()), TString::Format("%s/CECA_Paper/ScanPsUsmani/",GetFemtoOutputFolder()),
                  0.274-0.070*2, 0.274+0.070*2,//d
                  3.63-0.3*2, 3.63+0.3*2,//ht
                  2.68-0.3*2, 2.68+0.3*2,//hz
                  2311-50, 2311+50,
                  0.3453-0.008, 0.3453+0.008,
                  0.2595-0.002, 0.2595+0.002,
                  atoi(argv[1]),//lambda vars (400,401,402,500,501,502)
                  atof(argv[2]), atoi(argv[3]));//mins and seed

//thse parameters are what makes usmani fixed to NLO19
//Wc=2279.0; Rc=0.3394; Sc=0.2614; f1=1.41; d1=2.53; tDev=0.003
/*
ScanPsUsmani(
                  "Reduced","Cigar2_ds24_hts36_hzs36",
                  TString(GetCernBoxDimi()), TString::Format("%s/CECA_Paper/ScanPsUsmani/",GetFemtoOutputFolder()),
                  0.2,0.7,//d
                  0.4,4.8,//ht
                  4.0,12.0,//hz
                  2250, 2400,
                  0.3394-0.001, 0.3394+0.001,
                  0.2614-0.0005, 0.2614+0.0005,
                  atof(argv[1]), atoi(argv[2]));
  */
  return 0;
}
