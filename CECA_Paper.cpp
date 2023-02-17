
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
      Ivana.Ghetto_MtBins[0] = 930;
      Ivana.Ghetto_MtBins[1] = 1020;
      Ivana.Ghetto_MtBins[2] = 1080;
      Ivana.Ghetto_MtBins[3] = 1140;
      Ivana.Ghetto_MtBins[4] = 1200;
      Ivana.Ghetto_MtBins[5] = 1260;
      Ivana.Ghetto_MtBins[6] = 1380;
      Ivana.Ghetto_MtBins[7] = 1570;
      Ivana.Ghetto_MtBins[8] = 1840;
      Ivana.Ghetto_MtBins[9] = 2030;
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
      Ivana.Ghetto_MtBins[0] = 1000;
      Ivana.Ghetto_MtBins[1] = 1170;
      Ivana.Ghetto_MtBins[2] = 1250;
      Ivana.Ghetto_MtBins[3] = 1330;
      Ivana.Ghetto_MtBins[4] = 1430;
      Ivana.Ghetto_MtBins[5] = 1680;
      Ivana.Ghetto_MtBins[6] = 1840;
      Ivana.Ghetto_MtBins[7] = 2060;
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



//hadd but for my histos
bool dadd_f(int argc, char *argv[]){
  //printf("Hello\n");
  system("echo -n '1. Current Directory is '; pwd");
  if(argc<3) {printf("\033[1;31mERROR:\033[0m Not enough arguments [dadd Output Input]\n"); return false;}
  //1 default, 1 for the output
  const int NumInputFiles = argc-2;
  printf("Attempting to merge %i number of files\n",NumInputFiles);
  DLM_Histo<float> MergedHisto;
  //DLM_Histo<float>* HISTO;
  //DLM_Histo<float>& uHisto = HISTO[0];
  printf("Final file: %s\n",argv[1]);
  bool Status = true;
  for(unsigned uIn=0; uIn<NumInputFiles; uIn++){
    DLM_Histo<float> uHisto;
    printf(" -> %s\n",argv[2+uIn]);
    uHisto.QuickLoad(argv[2+uIn]);
    printf("   done\n");
    if(uIn==0){
      Status = (MergedHisto = uHisto);
    }
    else{
      Status = (MergedHisto += uHisto);
    }

    //delete HISTO;
    if(!Status) return Status;
  }
  MergedHisto.QuickWrite(argv[1],true);
  return true;
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

int CECA_PAPER(int argc, char *argv[]){

  //how to read/write the Levy pars into a file
  //TestSaveStuctToFile();

  //a test that extrapolation of Levy pars works
  //TestDoubleSourceOperation();

  //Test_pp_Statistics_1();

  return Ceca_pp_or_pL("TEST4_pp",TString::Format("%s/CECA_Paper/Ceca_pp_or_pL/TEST4/Job/",GetFemtoOutputFolder()),
  TString::Format("%s/CECA_Paper/Ceca_pp_or_pL/TEST4/Out/",GetFemtoOutputFolder()),
  TString::Format("%s/CECA_Paper/Ceca_pp_or_pL/TEST4/Log/",GetFemtoOutputFolder()),
  0,0,1);
  //ReadDlmHst();

  //dadd_f(argc,argv);

  //TestReadWriteBinary();

  return 0;
}
