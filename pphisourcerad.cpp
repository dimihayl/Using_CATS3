// R__ADD_INCLUDE_PATH("Users/emma/Software/CATS/install/include")

#include <iostream>
#include <stdio.h>
#include <string.h>
//#include <omp.h>
#include <complex>

// R__ADD_INCLUDE_PATH(Users/emma/Software/CATS/install/include)
// R__ADD_INCLUDE_PATH("/Users/emma/Software/CATS/source/CATS_Extentions")

// R__ADD_INCLUDE_PATH("/Users/emma/Software/CATS/source/CATS_MathTools")
// R__ADD_INCLUDE_PATH("/Users/emma/Software/CATS/source/CATS_CppTools")

// gSystem->AddIncludePath("/Users/emma/Software/CATS/source/CATS")

#include "CATS.h"
#include "CATSconstants.h"
#include "CATStools.h"
#include "DLM_Potentials.h"
#include "DLM_Source.h"
#include "DLM_Fitters.h"
#include "DLM_CppTools.h"
#include "DLM_CkDecomposition.h"
#include "DLM_Random.h"
#include "DLM_Bessel.h"
#include "DLM_Integration.h"
#include "DLM_WfModel.h"
#include "DLM_Histo.h"
#include "DLM_CkModels.h"
#include "DLM_HistoAnalysis.h"
#include "EnvVars.h"

#include "TROOT.h"
#include "TGraph.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TNtuple.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TNtuple.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TLorentzVector.h"
#include "TVector3.h"

using namespace std;



void GetCentralInterval(TH1* histo, float level, double &lowlim, double &uplim){
float mean=histo->GetMean();
cout<<"mean: "<<mean<<endl;
float sigma=histo->GetStdDev();
cout<<"sigm: "<<sigma<<endl;

float N=histo->GetEntries();
cout<<"N: "<<N<<endl;

float dev=level*sigma/sqrt(N);
cout<<"dev: "<<dev<<endl;

lowlim=mean-dev;
uplim=mean+dev;
}


void EffectiveGaussianpphi(){

//  const double CoreSize = 1.2;
    const double CoreSize = 0.94;

  //DLM_CleverMcLevyResoTM* MagicSource = new DLM_CleverMcLevyResoTM ();
  DLM_CleverMcLevyResoTM MagicSource;

  //DO NOT CHANGE !!! Sets up numerical bullshit, tuned for a Gaussian source
  MagicSource.InitStability(1,2-1e-6,2+1e-6);
  MagicSource.InitScale(38,0.15,2.0);
  MagicSource.InitRad(257*2,0,64);
  MagicSource.InitType(2);
  ///////////////////

  //for p-Xi, set up the amount of secondaries
  //first for the protons (64.22%)
  MagicSource.SetUpReso(0,0.6422);
  //than for the Xis, here its 0% (we have ONLY primordials)
  MagicSource.SetUpReso(1,0.0);

  //the cut off scale in k*, for which the angular distributions from EPOS
  //are evaluated. 200 MeV works okay, you can go up to 300 MeV for systematic checks
  const double k_CutOff = 200;

  //to be used for the NTuple later on
  Float_t k_D;
  Float_t fP1;
  Float_t fP2;
  Float_t fM1;
  Float_t fM2;
  Float_t Tau1;
  Float_t Tau2;
  Float_t AngleRcP1;
  Float_t AngleRcP2;
  Float_t AngleP1P2;
  //random generator dimi style. The input is incompatible with the ROOT random generator,
  //do not mix and match, do not ask me how I know this. Ask Bernie.
  //11 is the seed, you can change that to you favorite number
  DLM_Random RanGen(11);
  //dummies to save random shit
  double RanVal1;
  double RanVal2;
  double RanVal3;

  //open the magic file from dimi with the angular distributions.
//  TFile* F_EposDisto_pReso_Xim = new TFile("/home/emma/pphiRADIUS/files/EposDisto_p_pReso.root");
  TFile* F_EposDisto_pReso_Xim = new TFile(TString::Format("%s/CatsFiles/Source/EposAngularDist/Epos_pReso_Phi.root",GetCernBoxDimi()));

  //set up the ntuple, do not change anything unless told so by dimi
  TNtuple* T_EposDisto_pReso_Xim = (TNtuple*)F_EposDisto_pReso_Xim->Get("InfoTuple_ClosePairs");
  unsigned N_EposDisto_pReso_Xim = T_EposDisto_pReso_Xim->GetEntries();
  T_EposDisto_pReso_Xim->SetBranchAddress("k_D",&k_D);
  T_EposDisto_pReso_Xim->SetBranchAddress("P1",&fP1);
  T_EposDisto_pReso_Xim->SetBranchAddress("P2",&fP2);
  T_EposDisto_pReso_Xim->SetBranchAddress("M1",&fM1);
  T_EposDisto_pReso_Xim->SetBranchAddress("M2",&fM2);
  T_EposDisto_pReso_Xim->SetBranchAddress("Tau1",&Tau1);
  T_EposDisto_pReso_Xim->SetBranchAddress("Tau2",&Tau2);
  T_EposDisto_pReso_Xim->SetBranchAddress("AngleRcP1",&AngleRcP1);
  T_EposDisto_pReso_Xim->SetBranchAddress("AngleRcP2",&AngleRcP2);
  T_EposDisto_pReso_Xim->SetBranchAddress("AngleP1P2",&AngleP1P2);
  //iterate over the ntuple
  for(unsigned uEntry=0; uEntry<N_EposDisto_pReso_Xim; uEntry++){
      //get each entry
      T_EposDisto_pReso_Xim->GetEntry(uEntry);
      //disregard the entry of you are outside the desired k*
      if(k_D>k_CutOff) continue;
      //overwrite the value for the lifetime. This is computed from the
      //stat. hadronization model (Vale) or thermal fist (Max)
      //this is the value for the secondary protons
      Tau1 = 1.65;
      //for primoridials (the Xis) we put 0
      Tau2 = 0;
      //put in the average mass of the resonances (again from SHM or TF)
      //this is the value for protons
      fM1 = 1362;
      //generate a random path length for the propagation of the resonances
      //nothing to change!
      RanVal1 = RanGen.Exponential(fM1/(fP1*Tau1));
      //adds a single entry into the PDF for the angular distribution to be used
      MagicSource.AddBGT_RP(RanVal1,cos(AngleRcP1));
//      MagicSource.AddBGT_RP(RanVal1,-cos(AngleRcP1));


  }
  delete F_EposDisto_pReso_Xim;
  //if you have resonances contributing to both particles, we need to repeat the above procedure
  //for the prim-reso (AddBGT_PR) and reso-reso (AddBGT_RR) cases

  const unsigned NumSourceBins = 128;
  const double rMin = 0;
  const double rMax = 16;
  TFile* fOutput = new TFile("fOutput.root","recreate");
  TH1F* hSource = new TH1F("hSource","hSource",NumSourceBins,rMin,rMax);

  //fill the histo fro the source
  for(unsigned uBin=0; uBin<NumSourceBins; uBin++){
    //get the x-axis (r value) of the current bin
    double xaxis = hSource->GetBinCenter(uBin+1);
    //an array for the parameters, [0] is source size, [1] is == 2 (for a Gaussian)
    double parameters[2];
    parameters[0] = CoreSize;
    parameters[1] = 2.0;
    double SourceValue = MagicSource.RootEval(&xaxis, parameters);
    hSource->SetBinContent(uBin+1,SourceValue);
    //infinite errors for now
    hSource->SetBinError(uBin+1,1000.);
  }

  //idea: fit the source distribution only in a range around its peak
  //to do this: silly idea: put very large uncertainties in the bins outside of this range
  //we can get this range automatically, by evaluating the central (median) integral of the source distribution
  //with this set up, we fit the 68% most central yield of the source distribution
//  double lowerlimit=1.2;
//  double upperlimit=2.5;
  double lowerlimit=1.0;
  double upperlimit=3.5;
  //GetCentralInterval(hSource, 10, lowerlimit, upperlimit);

  unsigned lowerbin = hSource->FindBin(lowerlimit);
  unsigned upperbin = hSource->FindBin(upperlimit);
  for(unsigned uBin=lowerbin; uBin<=upperbin; uBin++){
    hSource->SetBinError(uBin+1,0.01);
  }

  printf("Core size of %.3f fm\n",CoreSize);
  printf("The fit will be performed in the range [%.2f, %.2f] fm\n",lowerlimit,upperlimit);
  //fyi, GaussSourceTF1 is in DLM_Source.h if you want to check it out.
  TF1* fSource = new TF1("fSource",GaussSourceTF1,rMin,rMax,1);
  fSource->SetParameter(0,CoreSize);
  fSource->SetParLimits(0,CoreSize*0.5,CoreSize*2.0);

  int binmax = hSource->GetMaximumBin();
  double maxval = hSource->GetBinContent(binmax);
  fSource->SetParameter(1,maxval);
  fSource->SetParLimits(1,maxval*0.95,maxval*1.05);

  hSource->Fit(fSource,"S, N, R, M");
  printf("The effective Gaussian size is %.3f +/- %.3f fm\n",fSource->GetParameter(0),fSource->GetParError(0));

  //get rid of weird plotting
  for(unsigned uBin=0; uBin<NumSourceBins; uBin++){
    hSource->SetBinError(uBin+1,0.01);
  }
  hSource->Write();
  fSource->Write();

  delete hSource;
  delete fSource;
  delete fOutput;
}



//first attempt: fit 0 - 16
//for the EB jun 2024: 0 - 4
void EmmaDaniel_piD(const double CoreSize, double *weight, double *r1, double *r2, double *chi2ndf, const double fMin, const double fMax)
{
    const int &SmoothSampling = 0;
    bool IncludeBug = false;
    const unsigned SmoothEntires = 1024 * 128;

    // DLM_CleverMcLevyResoTM* MagicSource = new DLM_CleverMcLevyResoTM ();
    DLM_CleverMcLevyResoTM MagicSource;

    // DO NOT CHANGE !!! Sets up numerical bullshit, tuned for a Gaussian source
    MagicSource.InitStability(1, 2 - 1e-6, 2 + 1e-6);
    MagicSource.InitScale(38, 0.15, 2.0);
    MagicSource.InitRad(257 * 2, 0, 64);
    MagicSource.InitType(2);
    ///////////////////

    // the amount of secondaries
    // pions (64.22%)
    MagicSource.SetUpReso(0, 0.0);
    // than for the Xis, here its 0% (we have ONLY primordials)
    MagicSource.SetUpReso(1, 0.682);
    MagicSource.InitNumMcIter(1000000);

    // the cut off scale in k*, for which the angular distributions from EPOS
    // are evaluated. 200 MeV works okay, you can go up to 300 MeV for systematic checks
    const double k_CutOff = 300;

    // to be used for the NTuple later on
    Float_t k_D;
    Float_t fP1;
    Float_t fP2;
    Float_t fM1;
    Float_t fM2;
    Float_t Tau1;
    Float_t Tau2;
    Float_t AngleRcP1;
    Float_t AngleRcP2;
    Float_t AngleP1P2;
    // random generator dimi style. The input is incompatible with the ROOT random generator,
    // do not mix and match, do not ask me how I know this. Ask Bernie.
    // 11 is the seed, you can change that to you favorite number
    DLM_Random RanGen(11);
    // dummies to save random shit
    double RanVal1;
    double RanVal2;
    double RanVal3;
    double RanCos;
    double MeanP2 = 0;
    // printf("1\n");
    // open the magic file from dimi with the angular distributions.
    //TFile *F_EposDisto = new TFile("/Users/emma/Femtoscopy/DMesonLight/oton_selections/radius/files/ALL_D_piReso.root");
    TFile *F_EposDisto = new TFile(TString::Format("%s/CatsFiles/Source/EposAngularDist/TEMP/SecondTry/ALL_D_piReso.root",GetCernBoxDimi()));
    // set up the ntuple, do not change anything unless told so by dimi
    TNtuple *T_EposDisto = (TNtuple *)F_EposDisto->Get("InfoTuple_ClosePairs");
    unsigned N_EposDisto = T_EposDisto->GetEntries();
    T_EposDisto->SetBranchAddress("k_D", &k_D);
    T_EposDisto->SetBranchAddress("P1", &fP1);
    T_EposDisto->SetBranchAddress("P2", &fP2);
    T_EposDisto->SetBranchAddress("M1", &fM1);
    T_EposDisto->SetBranchAddress("M2", &fM2);
    T_EposDisto->SetBranchAddress("Tau1", &Tau1);
    T_EposDisto->SetBranchAddress("Tau2", &Tau2);
    T_EposDisto->SetBranchAddress("AngleRcP1", &AngleRcP1);
    T_EposDisto->SetBranchAddress("AngleRcP2", &AngleRcP2);
    T_EposDisto->SetBranchAddress("AngleP1P2", &AngleP1P2);

    // TH1F* hAngleRcP1 = new TH1F("hAngleRcP1","hAngleRcP1",1024,0,TMath::Pi());

    gROOT->cd();
    TH1F *hAngle = new TH1F("hAngle", "hAngle", 32, 0, TMath::Pi());
    TH1F *hCos = new TH1F("hCos", "hCos", 32, -1., 1.);
    TH1F *hFinalAngle = new TH1F("hFinalAngle", "hFinalAngle", 32, 0, TMath::Pi());
    TH1F *hFinalCos = new TH1F("hFinalCos", "hFinalCos", 32, -1., 1.);
    F_EposDisto->cd();
    int NumUsefulEntries = 0;
    for (unsigned uEntry = 0; uEntry < N_EposDisto; uEntry++)
    {
        T_EposDisto->GetEntry(uEntry);
        if (k_D > k_CutOff)
            continue;
        hAngle->Fill(AngleRcP2);
        hCos->Fill(cos(AngleRcP2));
        MeanP2 += fP2;
        NumUsefulEntries++;
    }
    MeanP2 /= double(NumUsefulEntries);

    // hAngle->SetBinContent(hAngle->FindBin(0.275),0.31);
    hCos->Scale(1. / hAngle->Integral(), "width");
    // printf("2\n");

    if (!SmoothSampling)
    {
        // iterate over the ntuple
        for (unsigned uEntry = 0; uEntry < N_EposDisto; uEntry++)
        {
            // get each entry
            T_EposDisto->GetEntry(uEntry);
            // disregard the entry of you are outside the desired k*
            if (k_D > k_CutOff)
                continue;
            // overwrite the value for the lifetime. This is computed from the
            // stat. hadronization model (Vale) or thermal fist (Max)
            // this is the value for the secondary protons
            Tau1 = 0;
            // for primoridials (the Xis) we put 0
            Tau2 = 1.50;
            // put in the average mass of the resonances (again from SHM or TF)
            // this is the value for protons
            fM2 = 1180;
            // generate a random path length for the propagation of the resonances
            // nothing to change!
            RanVal2 = RanGen.Exponential(fM2 / (fP2 * Tau2));
            // adds a single entry into the PDF for the angular distribution to be used
            if (!IncludeBug)
                MagicSource.AddBGT_PR(RanVal2, cos(AngleRcP2));
            else
                MagicSource.AddBGT_PR(RanVal2, -cos(AngleRcP2));
            hFinalAngle->Fill(AngleRcP2);
            hFinalCos->Fill(cos(AngleRcP2));
        }
    }
    else
        cout << "LOL" << endl;
    // else if(SmoothSampling==1){
    //   gROOT->cd();
    //   TF1* fRandomAngle = MakeSmoothAngularSourceDisto(hAngle);
    //   //TRandom3* ran_gen = new TRandom3(11);
    //   for(unsigned uEntry=0; uEntry<SmoothEntires; uEntry++){
    //       Tau1 = 0;
    //       Tau2 = 1.50;
    //       fM2 = 1180;
    //       fP2 = MeanP2;
    //       RanVal2 = RanGen.Exponential(fM2/(fP2*Tau2));
    //       AngleRcP2 = fRandomAngle->GetRandom(0.,TMath::Pi());
    //       if(!IncludeBug) MagicSource.AddBGT_PR(RanVal2,cos(AngleRcP2));
    //       else MagicSource.AddBGT_PR(RanVal2,-cos(AngleRcP2));
    //       hFinalAngle->Fill(AngleRcP2);
    //       hFinalCos->Fill(cos(AngleRcP2));
    //   }
    //   gROOT->cd();
    //   //delete ran_gen;
    //   delete fRandomAngle;
    // }
    // else{
    //   gROOT->cd();
    //   TF1* fRandomCos = MakeSmoothAngularCosine(hCos);
    //   //TRandom3* ran_gen = new TRandom3(11);
    //   for(unsigned uEntry=0; uEntry<SmoothEntires; uEntry++){
    //       Tau1 = 0;
    //       Tau2 = 1.50;
    //       fM2 = 1180;
    //       fP2 = MeanP2;
    //       RanVal2 = RanGen.Exponential(fM2/(fP2*Tau2));
    //       RanCos = fRandomCos->GetRandom(-1.,1.);
    //       if(!IncludeBug) MagicSource.AddBGT_PR(RanVal2,RanCos);
    //       else MagicSource.AddBGT_PR(RanVal2,-RanCos);
    //       hFinalAngle->Fill(acos(RanCos));
    //       hFinalCos->Fill(RanCos);
    //   }
    //   gROOT->cd();
    //   delete fRandomCos;
    // }

    delete F_EposDisto;

    // if you have resonances contributing to both particles, we need to repeat the above procedure
    // for the prim-reso (AddBGT_PR) and reso-reso (AddBGT_RR) cases

    const unsigned NumSourceBins = 128;
    const double rMin = 0;
    const double rMax = 16;
    TFile *fOutput = new TFile(TString::Format("fSecond_piD_%.2f.root", CoreSize), "recreate");
    TH1F *hSource = new TH1F("hSource", "hSource", NumSourceBins, rMin, rMax);

    // fill the histo fro the source
    for (unsigned uBin = 0; uBin < NumSourceBins; uBin++)
    {
        // get the x-axis (r value) of the current bin
        double xaxis = hSource->GetBinCenter(uBin + 1);
        // an array for the parameters, [0] is source size, [1] is == 2 (for a Gaussian)
        double parameters[2];
        parameters[0] = CoreSize;
        parameters[1] = 2.0;
        double SourceValue = MagicSource.RootEval(&xaxis, parameters);
        hSource->SetBinContent(uBin + 1, SourceValue);
        // infinite errors for now
        hSource->SetBinError(uBin + 1, 1000.);
    }
    // printf("4\n");
    // idea: fit the source distribution only in a range around its peak
    // to do this: silly idea: put very large uncertainties in the bins outside of this range
    // we can get this range automatically, by evaluating the central (median) integral of the source distribution
    // with this set up, we fit the 68% most central yield of the source distribution
    double lowerlimit;
    double upperlimit;
    GetCentralInterval(*hSource, 0.84, lowerlimit, upperlimit, true);
    unsigned lowerbin = hSource->FindBin(lowerlimit);
    unsigned upperbin = hSource->FindBin(upperlimit);
    for (unsigned uBin = lowerbin; uBin <= upperbin; uBin++)
    {
        hSource->SetBinError(uBin + 1, 0.01);
    }


//D(*)K: 0.86+0.09-0.07, D(*)pi: 0.97+0.09-0.08
    printf("Core size of %.3f fm\n", CoreSize);
    printf("The fit will be performed in the range [%.2f, %.2f] fm\n", lowerlimit, upperlimit);
    // fyi, GaussSourceTF1 is in DLM_Source.h if you want to check it out.

    // par0 = weight
    // par1 = radius 1
    // par2 = radius 2

    TF1 *fSource = new TF1("fSource",
                           "[0] * 4 * TMath::Pi() * x ** 2 * 1. / ((4 * TMath::Pi() *[1] ** 2) ** 1.5) * TMath::E() ** "
                           "(- x ** 2 / 4 /[1] ** 2) + (1 - [0]) * 4 * TMath::Pi() * x ** 2 * 1. / ((4 * TMath::Pi() * "
                           "[2] ** 2) ** 1.5) * TMath::E() ** (- x ** 2 / 4 / [2] ** 2)",
                           fMin, fMax);
    fSource->SetParameter(0, 0.5);
    fSource->SetParLimits(0, 0, 1);
    fSource->SetParameter(1, CoreSize);
    // fSource->SetParLimits(1,CoreSize*0.5,CoreSize*2.0);
    fSource->SetParLimits(1, 0, 10);
    fSource->SetParLimits(2, fSource->GetParameter(1), 10);

//GHETTO DIMI STUFF
    //reff1 0.970 --> 0.931
    //reff2 2.520 --> 2.419

    fSource->FixParameter(0, 0.66);
    fSource->FixParameter(1, 0.931);
    fSource->FixParameter(2, 2.419);

    hSource->Fit(fSource, "Q, S, N, R, M");
    printf("The weight size is %.3f +/- %.3f fm\n", fSource->GetParameter(0), fSource->GetParError(0));
    printf("The effective Gaussian size 1 is %.3f +/- %.3f fm\n", fSource->GetParameter(1), fSource->GetParError(1));
    printf("The effective Gaussian size 2 is %.3f +/- %.3f fm\n", fSource->GetParameter(2), fSource->GetParError(2));

    *weight = fSource->GetParameter(0);
    *r1 = fSource->GetParameter(1);
    *r2 = fSource->GetParameter(2);
    *chi2ndf = fSource->GetChisquare() / fSource->GetNDF();

    TF1 *fSourceLambda = new TF1("fSourceLambda",
                                 "[3]* ( [0] * 4 * TMath::Pi() * x ** 2 * 1. / ((4 * TMath::Pi() *[1] ** 2) ** 1.5) * "
                                 "TMath::E() ** (- x ** 2 / 4 /[1] ** 2) + (1 - [0]) * 4 * TMath::Pi() * x ** 2 * 1. / "
                                 "((4 * TMath::Pi() * [2] ** 2) ** 1.5) * TMath::E() ** (- x ** 2 / 4 / [2] ** 2))",
                                 fMin, fMax);
    fSourceLambda->SetParameter(0, 0.5);
    fSourceLambda->SetParLimits(0, 0, 1);
    fSourceLambda->SetParameter(1, CoreSize);
    // fSourceLambda->SetParLimits(1,CoreSize*0.5,CoreSize*2.0);
    fSourceLambda->SetParLimits(1, 0, 10);
    fSourceLambda->SetParLimits(2, fSourceLambda->GetParameter(1), 10);
    fSourceLambda->SetParLimits(3, 0, 1);

    hSource->Fit(fSourceLambda, "Q, S, N, R, M");

    // get rid of weird plotting
    for (unsigned uBin = 0; uBin < NumSourceBins; uBin++)
    {
        hSource->SetBinError(uBin + 1, 0.01);
    }
    fOutput->cd();
    hSource->Write();
    fSource->Write();
    fSourceLambda->Write();
    hAngle->Scale(1. / hAngle->Integral(), "width");
    hAngle->Write();
    hCos->Write();
    hFinalAngle->Scale(1. / hFinalAngle->Integral(), "width");
    hFinalAngle->Write();
    hFinalCos->Scale(1. / hFinalCos->Integral(), "width");
    hFinalCos->Write();

    // gROOT->cd();
    // delete hAngle;

    delete hSource;
    delete fSource;
    delete hAngle;
    delete hCos;
    delete hFinalAngle;
    delete hFinalCos;
    delete fOutput;
    printf("-------------------\n");
    // return 0;
}

void EmmaDaniel_piD_error(const double CoreSize, double *weight, double *weight_error, double *r1, double *r1_error,
                          double *r2, double *r2_error, double *chi2ndf)
{
    const int &SmoothSampling = 0;
    bool IncludeBug = false;
    const unsigned SmoothEntires = 1024 * 128;

    // DLM_CleverMcLevyResoTM* MagicSource = new DLM_CleverMcLevyResoTM ();
    DLM_CleverMcLevyResoTM MagicSource;

    // DO NOT CHANGE !!! Sets up numerical bullshit, tuned for a Gaussian source
    MagicSource.InitStability(1, 2 - 1e-6, 2 + 1e-6);
    MagicSource.InitScale(38, 0.15, 2.0);
    MagicSource.InitRad(257 * 2, 0, 64);
    MagicSource.InitType(2);
    ///////////////////

    // the amount of secondaries
    // pions (64.22%)
    MagicSource.SetUpReso(0, 0.0);
    // than for the Xis, here its 0% (we have ONLY primordials)
    MagicSource.SetUpReso(1, 0.682);
    MagicSource.InitNumMcIter(1000000);

    // the cut off scale in k*, for which the angular distributions from EPOS
    // are evaluated. 200 MeV works okay, you can go up to 300 MeV for systematic checks
    const double k_CutOff = 300;

    // to be used for the NTuple later on
    Float_t k_D;
    Float_t fP1;
    Float_t fP2;
    Float_t fM1;
    Float_t fM2;
    Float_t Tau1;
    Float_t Tau2;
    Float_t AngleRcP1;
    Float_t AngleRcP2;
    Float_t AngleP1P2;
    // random generator dimi style. The input is incompatible with the ROOT random generator,
    // do not mix and match, do not ask me how I know this. Ask Bernie.
    // 11 is the seed, you can change that to you favorite number
    DLM_Random RanGen(11);
    // dummies to save random shit
    double RanVal1;
    double RanVal2;
    double RanVal3;
    double RanCos;
    double MeanP2 = 0;
    // printf("1\n");
    // open the magic file from dimi with the angular distributions.
    TFile *F_EposDisto = new TFile("/Users/emma/Femtoscopy/DMesonLight/oton_selections/radius/files/ALL_D_piReso.root");
    // set up the ntuple, do not change anything unless told so by dimi
    TNtuple *T_EposDisto = (TNtuple *)F_EposDisto->Get("InfoTuple_ClosePairs");
    unsigned N_EposDisto = T_EposDisto->GetEntries();
    T_EposDisto->SetBranchAddress("k_D", &k_D);
    T_EposDisto->SetBranchAddress("P1", &fP1);
    T_EposDisto->SetBranchAddress("P2", &fP2);
    T_EposDisto->SetBranchAddress("M1", &fM1);
    T_EposDisto->SetBranchAddress("M2", &fM2);
    T_EposDisto->SetBranchAddress("Tau1", &Tau1);
    T_EposDisto->SetBranchAddress("Tau2", &Tau2);
    T_EposDisto->SetBranchAddress("AngleRcP1", &AngleRcP1);
    T_EposDisto->SetBranchAddress("AngleRcP2", &AngleRcP2);
    T_EposDisto->SetBranchAddress("AngleP1P2", &AngleP1P2);

    // TH1F* hAngleRcP1 = new TH1F("hAngleRcP1","hAngleRcP1",1024,0,TMath::Pi());

    gROOT->cd();
    TH1F *hAngle = new TH1F("hAngle", "hAngle", 32, 0, TMath::Pi());
    TH1F *hCos = new TH1F("hCos", "hCos", 32, -1., 1.);
    TH1F *hFinalAngle = new TH1F("hFinalAngle", "hFinalAngle", 32, 0, TMath::Pi());
    TH1F *hFinalCos = new TH1F("hFinalCos", "hFinalCos", 32, -1., 1.);
    F_EposDisto->cd();
    int NumUsefulEntries = 0;
    for (unsigned uEntry = 0; uEntry < N_EposDisto; uEntry++)
    {
        T_EposDisto->GetEntry(uEntry);
        if (k_D > k_CutOff)
            continue;
        hAngle->Fill(AngleRcP2);
        hCos->Fill(cos(AngleRcP2));
        MeanP2 += fP2;
        NumUsefulEntries++;
    }
    MeanP2 /= double(NumUsefulEntries);

    // hAngle->SetBinContent(hAngle->FindBin(0.275),0.31);
    hCos->Scale(1. / hAngle->Integral(), "width");
    // printf("2\n");

    if (!SmoothSampling)
    {
        // iterate over the ntuple
        for (unsigned uEntry = 0; uEntry < N_EposDisto; uEntry++)
        {
            // get each entry
            T_EposDisto->GetEntry(uEntry);
            // disregard the entry of you are outside the desired k*
            if (k_D > k_CutOff)
                continue;
            // overwrite the value for the lifetime. This is computed from the
            // stat. hadronization model (Vale) or thermal fist (Max)
            // this is the value for the secondary protons
            Tau1 = 0;
            // for primoridials (the Xis) we put 0
            Tau2 = 1.50;
            // put in the average mass of the resonances (again from SHM or TF)
            // this is the value for protons
            fM2 = 1180;
            // generate a random path length for the propagation of the resonances
            // nothing to change!
            RanVal2 = RanGen.Exponential(fM2 / (fP2 * Tau2));
            // adds a single entry into the PDF for the angular distribution to be used
            if (!IncludeBug)
                MagicSource.AddBGT_PR(RanVal2, cos(AngleRcP2));
            else
                MagicSource.AddBGT_PR(RanVal2, -cos(AngleRcP2));
            hFinalAngle->Fill(AngleRcP2);
            hFinalCos->Fill(cos(AngleRcP2));
        }
    }
    else
        cout << "LOL" << endl;
    // else if(SmoothSampling==1){
    //   gROOT->cd();
    //   TF1* fRandomAngle = MakeSmoothAngularSourceDisto(hAngle);
    //   //TRandom3* ran_gen = new TRandom3(11);
    //   for(unsigned uEntry=0; uEntry<SmoothEntires; uEntry++){
    //       Tau1 = 0;
    //       Tau2 = 1.50;
    //       fM2 = 1180;
    //       fP2 = MeanP2;
    //       RanVal2 = RanGen.Exponential(fM2/(fP2*Tau2));
    //       AngleRcP2 = fRandomAngle->GetRandom(0.,TMath::Pi());
    //       if(!IncludeBug) MagicSource.AddBGT_PR(RanVal2,cos(AngleRcP2));
    //       else MagicSource.AddBGT_PR(RanVal2,-cos(AngleRcP2));
    //       hFinalAngle->Fill(AngleRcP2);
    //       hFinalCos->Fill(cos(AngleRcP2));
    //   }
    //   gROOT->cd();
    //   //delete ran_gen;
    //   delete fRandomAngle;
    // }
    // else{
    //   gROOT->cd();
    //   TF1* fRandomCos = MakeSmoothAngularCosine(hCos);
    //   //TRandom3* ran_gen = new TRandom3(11);
    //   for(unsigned uEntry=0; uEntry<SmoothEntires; uEntry++){
    //       Tau1 = 0;
    //       Tau2 = 1.50;
    //       fM2 = 1180;
    //       fP2 = MeanP2;
    //       RanVal2 = RanGen.Exponential(fM2/(fP2*Tau2));
    //       RanCos = fRandomCos->GetRandom(-1.,1.);
    //       if(!IncludeBug) MagicSource.AddBGT_PR(RanVal2,RanCos);
    //       else MagicSource.AddBGT_PR(RanVal2,-RanCos);
    //       hFinalAngle->Fill(acos(RanCos));
    //       hFinalCos->Fill(RanCos);
    //   }
    //   gROOT->cd();
    //   delete fRandomCos;
    // }

    delete F_EposDisto;

    // if you have resonances contributing to both particles, we need to repeat the above procedure
    // for the prim-reso (AddBGT_PR) and reso-reso (AddBGT_RR) cases

    const unsigned NumSourceBins = 128;
    const double rMin = 0;
    const double rMax = 16;
    TFile *fOutput = new TFile(TString::Format("fSecond_piD_%.2f.root", CoreSize), "recreate");
    TH1F *hSource = new TH1F("hSource", "hSource", NumSourceBins, rMin, rMax);

    // fill the histo fro the source
    for (unsigned uBin = 0; uBin < NumSourceBins; uBin++)
    {
        // get the x-axis (r value) of the current bin
        double xaxis = hSource->GetBinCenter(uBin + 1);
        // an array for the parameters, [0] is source size, [1] is == 2 (for a Gaussian)
        double parameters[2];
        parameters[0] = CoreSize;
        parameters[1] = 2.0;
        double SourceValue = MagicSource.RootEval(&xaxis, parameters);
        hSource->SetBinContent(uBin + 1, SourceValue);
        // infinite errors for now
        hSource->SetBinError(uBin + 1, 1000.);
    }
    // printf("4\n");
    // idea: fit the source distribution only in a range around its peak
    // to do this: silly idea: put very large uncertainties in the bins outside of this range
    // we can get this range automatically, by evaluating the central (median) integral of the source distribution
    // with this set up, we fit the 68% most central yield of the source distribution
    double lowerlimit;
    double upperlimit;
    GetCentralInterval(*hSource, 0.84, lowerlimit, upperlimit, true);
    unsigned lowerbin = hSource->FindBin(lowerlimit);
    unsigned upperbin = hSource->FindBin(upperlimit);
    for (unsigned uBin = lowerbin; uBin <= upperbin; uBin++)
    {
        hSource->SetBinError(uBin + 1, 0.01);
    }

    printf("Core size of %.3f fm\n", CoreSize);
    printf("The fit will be performed in the range [%.2f, %.2f] fm\n", lowerlimit, upperlimit);
    // fyi, GaussSourceTF1 is in DLM_Source.h if you want to check it out.

    // par0 = weight
    // par1 = radius 1
    // par2 = radius 2

    TF1 *fSource = new TF1("fSource",
                           "[0] * 4 * TMath::Pi() * x ** 2 * 1. / ((4 * TMath::Pi() *[1] ** 2) ** 1.5) * TMath::E() ** "
                           "(- x ** 2 / 4 /[1] ** 2) + (1 - [0]) * 4 * TMath::Pi() * x ** 2 * 1. / ((4 * TMath::Pi() * "
                           "[2] ** 2) ** 1.5) * TMath::E() ** (- x ** 2 / 4 / [2] ** 2)",
                           rMin, rMax);
    fSource->SetParameter(0, 0.5);
    fSource->SetParLimits(0, 0.4, 1);
    fSource->SetParameter(1, CoreSize);
    // fSource->SetParLimits(1,CoreSize*0.5,CoreSize*2.0);
    fSource->SetParLimits(1, 0, 10);
    fSource->SetParLimits(2, fSource->GetParameter(1), 10);

    hSource->Fit(fSource, "Q, S, N, R, M");
    printf("The weight size is %.3f +/- %.3f fm\n", fSource->GetParameter(0), fSource->GetParError(0));
    printf("The effective Gaussian size 1 is %.3f +/- %.3f fm\n", fSource->GetParameter(1), fSource->GetParError(1));
    printf("The effective Gaussian size 2 is %.3f +/- %.3f fm\n", fSource->GetParameter(2), fSource->GetParError(2));

    *weight = fSource->GetParameter(0);
    *weight_error = fSource->GetParError(0);

    *r1 = fSource->GetParameter(1);
    *r1_error = fSource->GetParError(1);

    *r2 = fSource->GetParameter(2);
    *r2_error = fSource->GetParError(2);

    *chi2ndf = fSource->GetChisquare() / fSource->GetNDF();

    printf("The chi2/ndf is %.2f fm\n", fSource->GetChisquare() / fSource->GetNDF());

    TF1 *fSourceLambda = new TF1("fSourceLambda",
                                 "[3]* ( [0] * 4 * TMath::Pi() * x ** 2 * 1. / ((4 * TMath::Pi() *[1] ** 2) ** 1.5) * "
                                 "TMath::E() ** (- x ** 2 / 4 /[1] ** 2) + (1 - [0]) * 4 * TMath::Pi() * x ** 2 * 1. / "
                                 "((4 * TMath::Pi() * [2] ** 2) ** 1.5) * TMath::E() ** (- x ** 2 / 4 / [2] ** 2))",
                                 rMin, rMax);
    fSourceLambda->SetParameter(0, 0.5);
    fSourceLambda->SetParLimits(0, 0, 1);
    fSourceLambda->SetParameter(1, CoreSize);
    // fSourceLambda->SetParLimits(1,CoreSize*0.5,CoreSize*2.0);
    fSourceLambda->SetParLimits(1, 0, 10);
    fSourceLambda->SetParLimits(2, fSourceLambda->GetParameter(1), 10);
    fSourceLambda->SetParLimits(3, 0, 1);

    hSource->Fit(fSourceLambda, "Q, S, N, R, M");

    // get rid of weird plotting
    for (unsigned uBin = 0; uBin < NumSourceBins; uBin++)
    {
        hSource->SetBinError(uBin + 1, 0.01);
    }
    fOutput->cd();
    hSource->Write();
    fSource->Write();
    fSourceLambda->Write();
    hAngle->Scale(1. / hAngle->Integral(), "width");
    hAngle->Write();
    hCos->Write();
    hFinalAngle->Scale(1. / hFinalAngle->Integral(), "width");
    hFinalAngle->Write();
    hFinalCos->Scale(1. / hFinalCos->Integral(), "width");
    hFinalCos->Write();

    // gROOT->cd();
    // delete hAngle;

    delete hSource;
    delete fSource;
    delete hAngle;
    delete hCos;
    delete hFinalAngle;
    delete hFinalCos;
    delete fOutput;
    printf("-------------------\n");
    // return 0;
}

void EmmaDaniel_KD(const double CoreSize, double *weight, double *r1, double *r2, double *chi2ndf, const double fMin, const double fMax)
{
    const int &SmoothSampling = 0;
    bool IncludeBug = false;
    const unsigned SmoothEntires = 1024 * 128;

    // DLM_CleverMcLevyResoTM* MagicSource = new DLM_CleverMcLevyResoTM ();
    DLM_CleverMcLevyResoTM MagicSource;

    // DO NOT CHANGE !!! Sets up numerical bullshit, tuned for a Gaussian source
    MagicSource.InitStability(1, 2 - 1e-6, 2 + 1e-6);
    MagicSource.InitScale(38, 0.15, 2.0);
    MagicSource.InitRad(257 * 2, 0, 64);
    MagicSource.InitType(2);
    ///////////////////

    // the amount of secondaries
    // pions (64.22%)
    MagicSource.SetUpReso(0, 0.0);
    // than for the Xis, here its 0% (we have ONLY primordials)
    MagicSource.SetUpReso(1, 0.476);
    MagicSource.InitNumMcIter(1000000);

    // the cut off scale in k*, for which the angular distributions from EPOS
    // are evaluated. 200 MeV works okay, you can go up to 300 MeV for systematic checks
    const double k_CutOff = 300;

    // to be used for the NTuple later on
    Float_t k_D;
    Float_t fP1;
    Float_t fP2;
    Float_t fM1;
    Float_t fM2;
    Float_t Tau1;
    Float_t Tau2;
    Float_t AngleRcP1;
    Float_t AngleRcP2;
    Float_t AngleP1P2;
    // random generator dimi style. The input is incompatible with the ROOT random generator,
    // do not mix and match, do not ask me how I know this. Ask Bernie.
    // 11 is the seed, you can change that to you favorite number
    DLM_Random RanGen(11);
    // dummies to save random shit
    double RanVal1;
    double RanVal2;
    double RanVal3;
    double RanCos;
    double MeanP2 = 0;
    // printf("1\n");
    // open the magic file from dimi with the angular distributions.
    //TFile *F_EposDisto =
    //    new TFile("/Users/emma/Femtoscopy/DMesonLight/oton_selections/radius/files/ALL_D_KaonReso.root");
    TFile *F_EposDisto = new TFile(TString::Format("%s/CatsFiles/Source/EposAngularDist/TEMP/SecondTry/ALL_D_KaonReso.root",GetCernBoxDimi()));
    // set up the ntuple, do not change anything unless told so by dimi
    TNtuple *T_EposDisto = (TNtuple *)F_EposDisto->Get("InfoTuple_ClosePairs");
    unsigned N_EposDisto = T_EposDisto->GetEntries();
    T_EposDisto->SetBranchAddress("k_D", &k_D);
    T_EposDisto->SetBranchAddress("P1", &fP1);
    T_EposDisto->SetBranchAddress("P2", &fP2);
    T_EposDisto->SetBranchAddress("M1", &fM1);
    T_EposDisto->SetBranchAddress("M2", &fM2);
    T_EposDisto->SetBranchAddress("Tau1", &Tau1);
    T_EposDisto->SetBranchAddress("Tau2", &Tau2);
    T_EposDisto->SetBranchAddress("AngleRcP1", &AngleRcP1);
    T_EposDisto->SetBranchAddress("AngleRcP2", &AngleRcP2);
    T_EposDisto->SetBranchAddress("AngleP1P2", &AngleP1P2);

    // TH1F* hAngleRcP1 = new TH1F("hAngleRcP1","hAngleRcP1",1024,0,TMath::Pi());

    gROOT->cd();
    TH1F *hAngle = new TH1F("hAngle", "hAngle", 32, 0, TMath::Pi());
    TH1F *hCos = new TH1F("hCos", "hCos", 32, -1., 1.);
    TH1F *hFinalAngle = new TH1F("hFinalAngle", "hFinalAngle", 32, 0, TMath::Pi());
    TH1F *hFinalCos = new TH1F("hFinalCos", "hFinalCos", 32, -1., 1.);
    F_EposDisto->cd();
    int NumUsefulEntries = 0;
    for (unsigned uEntry = 0; uEntry < N_EposDisto; uEntry++)
    {
        T_EposDisto->GetEntry(uEntry);
        if (k_D > k_CutOff)
            continue;
        hAngle->Fill(AngleRcP2);
        hCos->Fill(cos(AngleRcP2));
        MeanP2 += fP2;
        NumUsefulEntries++;
    }
    MeanP2 /= double(NumUsefulEntries);

    // hAngle->SetBinContent(hAngle->FindBin(0.275),0.31);
    hCos->Scale(1. / hAngle->Integral(), "width");
    // printf("2\n");

    if (!SmoothSampling)
    {
        // iterate over the ntuple
        for (unsigned uEntry = 0; uEntry < N_EposDisto; uEntry++)
        {
            // get each entry
            T_EposDisto->GetEntry(uEntry);
            // disregard the entry of you are outside the desired k*
            if (k_D > k_CutOff)
                continue;
            // overwrite the value for the lifetime. This is computed from the
            // stat. hadronization model (Vale) or thermal fist (Max)
            // this is the value for the secondary protons
            Tau1 = 0;
            // for primoridials (the Xis) we put 0
            Tau2 = 3.66;
            // put in the average mass of the resonances (again from SHM or TF)
            // this is the value for protons
            fM2 = 1054;
            // generate a random path length for the propagation of the resonances
            // nothing to change!
            RanVal2 = RanGen.Exponential(fM2 / (fP2 * Tau2));
            // adds a single entry into the PDF for the angular distribution to be used
            if (!IncludeBug)
                MagicSource.AddBGT_PR(RanVal2, cos(AngleRcP2));
            else
                MagicSource.AddBGT_PR(RanVal2, -cos(AngleRcP2));
            hFinalAngle->Fill(AngleRcP2);
            hFinalCos->Fill(cos(AngleRcP2));
        }
    }
    else
        cout << "LOL KAON" << endl;
    // else if(SmoothSampling==1){
    //   gROOT->cd();
    //   TF1* fRandomAngle = MakeSmoothAngularSourceDisto(hAngle);
    //   //TRandom3* ran_gen = new TRandom3(11);
    //   for(unsigned uEntry=0; uEntry<SmoothEntires; uEntry++){
    //       Tau1 = 0;
    //       Tau2 = 3.66;
    //       fM2 = 1054;
    //       fP2 = MeanP2;
    //       RanVal2 = RanGen.Exponential(fM2/(fP2*Tau2));
    //       AngleRcP2 = fRandomAngle->GetRandom(0.,TMath::Pi());
    //       if(!IncludeBug) MagicSource.AddBGT_PR(RanVal2,cos(AngleRcP2));
    //       else MagicSource.AddBGT_PR(RanVal2,-cos(AngleRcP2));
    //       hFinalAngle->Fill(AngleRcP2);
    //       hFinalCos->Fill(cos(AngleRcP2));
    //   }
    //   gROOT->cd();
    //   //delete ran_gen;
    //   delete fRandomAngle;
    // }
    // else{
    //   gROOT->cd();
    //   TF1* fRandomCos = MakeSmoothAngularCosine(hCos);
    //   //TRandom3* ran_gen = new TRandom3(11);
    //   for(unsigned uEntry=0; uEntry<SmoothEntires; uEntry++){
    //       Tau1 = 0;
    //       Tau2 = 3.66;
    //       fM2 = 1054;
    //       fP2 = MeanP2;
    //       RanVal2 = RanGen.Exponential(fM2/(fP2*Tau2));
    //       RanCos = fRandomCos->GetRandom(-1.,1.);
    //       if(!IncludeBug) MagicSource.AddBGT_PR(RanVal2,RanCos);
    //       else MagicSource.AddBGT_PR(RanVal2,RanCos);
    //       hFinalAngle->Fill(acos(RanCos));
    //       hFinalCos->Fill(RanCos);
    //   }
    //   gROOT->cd();
    //   delete fRandomCos;
    // }

    delete F_EposDisto;

    // if you have resonances contributing to both particles, we need to repeat the above procedure
    // for the prim-reso (AddBGT_PR) and reso-reso (AddBGT_RR) cases

    const unsigned NumSourceBins = 128;
    const double rMin = 0;
    const double rMax = 16;
    TFile *fOutput = new TFile(TString::Format("fSecond_KaonD_%.2f.root", CoreSize), "recreate");
    TH1F *hSource = new TH1F("hSource", "hSource", NumSourceBins, rMin, rMax);

    // fill the histo fro the source
    for (unsigned uBin = 0; uBin < NumSourceBins; uBin++)
    {
        // get the x-axis (r value) of the current bin
        double xaxis = hSource->GetBinCenter(uBin + 1);
        // an array for the parameters, [0] is source size, [1] is == 2 (for a Gaussian)
        double parameters[2];
        parameters[0] = CoreSize;
        parameters[1] = 2.0;
        double SourceValue = MagicSource.RootEval(&xaxis, parameters);
        hSource->SetBinContent(uBin + 1, SourceValue);
        // infinite errors for now
        hSource->SetBinError(uBin + 1, 1000.);
    }
    // printf("4\n");
    // idea: fit the source distribution only in a range around its peak
    // to do this: silly idea: put very large uncertainties in the bins outside of this range
    // we can get this range automatically, by evaluating the central (median) integral of the source distribution
    // with this set up, we fit the 68% most central yield of the source distribution
    double lowerlimit;
    double upperlimit;
    GetCentralInterval(*hSource, 0.84, lowerlimit, upperlimit, true);
    unsigned lowerbin = hSource->FindBin(lowerlimit);
    unsigned upperbin = hSource->FindBin(upperlimit);
    for (unsigned uBin = lowerbin; uBin <= upperbin; uBin++)
    {
        hSource->SetBinError(uBin + 1, 0.01);
    }

    printf("Core size of %.3f fm\n", CoreSize);
    printf("The fit will be performed in the range [%.2f, %.2f] fm\n", lowerlimit, upperlimit);
    // fyi, GaussSourceTF1 is in DLM_Source.h if you want to check it out.
    TF1 *fSource = new TF1("fSource",
                           "[0] * 4 * TMath::Pi() * x ** 2 * 1. / ((4 * TMath::Pi() *[1] ** 2) ** 1.5) * TMath::E() ** "
                           "(- x ** 2 / 4 /[1] ** 2) + (1 - [0]) * 4 * TMath::Pi() * x ** 2 * 1. / ((4 * TMath::Pi() * "
                           "[2] ** 2) ** 1.5) * TMath::E() ** (- x ** 2 / 4 / [2] ** 2)",
                           fMin, fMax);
    fSource->SetParameter(0, 0.5);
    fSource->SetParLimits(0, 0, 1);
    fSource->SetParameter(1, CoreSize);
    // fSource->SetParLimits(1,CoreSize*0.5,CoreSize*2.0);
    fSource->SetParLimits(1, 0, 10);
    fSource->SetParLimits(2, fSource->GetParameter(1), 10);

if(CoreSize==0.81){
    fSource->FixParameter(0, 0.78);
fSource->FixParameter(1, 0.86);
fSource->FixParameter(2, 2.03);
}
else{
    fSource->FixParameter(0, 0.78);
fSource->FixParameter(1, 0.814);
fSource->FixParameter(2, 1.921);
}

    hSource->Fit(fSource, "Q, S, N, R, M");
    printf("The weight size is %.3f +/- %.3f fm\n", fSource->GetParameter(0), fSource->GetParError(0));
    printf("The effective Gaussian size 1 is %.3f +/- %.3f fm\n", fSource->GetParameter(1), fSource->GetParError(1));
    printf("The effective Gaussian size 2 is %.3f +/- %.3f fm\n", fSource->GetParameter(2), fSource->GetParError(2));

    *weight = fSource->GetParameter(0);
    *r1 = fSource->GetParameter(1);
    *r2 = fSource->GetParameter(2);
    *chi2ndf = fSource->GetChisquare() / fSource->GetNDF();

    TF1 *fSourceLambda = new TF1("fSourceLambda",
                                 "[3]* ( [0] * 4 * TMath::Pi() * x ** 2 * 1. / ((4 * TMath::Pi() *[1] ** 2) ** 1.5) * "
                                 "TMath::E() ** (- x ** 2 / 4 /[1] ** 2) + (1 - [0]) * 4 * TMath::Pi() * x ** 2 * 1. / "
                                 "((4 * TMath::Pi() * [2] ** 2) ** 1.5) * TMath::E() ** (- x ** 2 / 4 / [2] ** 2))",
                                 rMin, rMax);
    fSourceLambda->SetParameter(0, 0.5);
    fSourceLambda->SetParLimits(0, 0, 1);
    fSourceLambda->SetParameter(1, CoreSize);
    // fSourceLambda->SetParLimits(1,CoreSize*0.5,CoreSize*2.0);
    fSourceLambda->SetParLimits(1, 0, 10);
    fSourceLambda->SetParLimits(2, fSourceLambda->GetParameter(1), 10);
    fSourceLambda->SetParLimits(3, 0, 1);

    hSource->Fit(fSourceLambda, "Q, S, N, R, M");

    // get rid of weird plotting
    for (unsigned uBin = 0; uBin < NumSourceBins; uBin++)
    {
        hSource->SetBinError(uBin + 1, 0.01);
    }
    hSource->Write();
    fSource->Write();
    fSourceLambda->Write();
    hAngle->Scale(1. / hAngle->Integral(), "width");
    hAngle->Write();
    hCos->Write();
    hFinalAngle->Scale(1. / hFinalAngle->Integral(), "width");
    hFinalAngle->Write();
    hFinalCos->Scale(1. / hFinalCos->Integral(), "width");
    hFinalCos->Write();

    // gROOT->cd();
    // delete hAngle;

    delete hSource;
    delete fSource;
    delete hAngle;
    delete hCos;
    delete hFinalAngle;
    delete hFinalCos;
    delete fOutput;

    printf("-------------------\n");
}

void EmmaDaniel_KD_error(const double CoreSize, double *weight, double *weight_error, double *r1, double *r1_error,
                         double *r2, double *r2_error, double *chi2ndf)
{
    const int &SmoothSampling = 0;
    bool IncludeBug = false;
    const unsigned SmoothEntires = 1024 * 128;

    // DLM_CleverMcLevyResoTM* MagicSource = new DLM_CleverMcLevyResoTM ();
    DLM_CleverMcLevyResoTM MagicSource;

    // DO NOT CHANGE !!! Sets up numerical bullshit, tuned for a Gaussian source
    MagicSource.InitStability(1, 2 - 1e-6, 2 + 1e-6);
    MagicSource.InitScale(38, 0.15, 2.0);
    MagicSource.InitRad(257 * 2, 0, 64);
    MagicSource.InitType(2);
    ///////////////////

    // the amount of secondaries
    // pions (64.22%)
    MagicSource.SetUpReso(0, 0.0);
    // than for the Xis, here its 0% (we have ONLY primordials)
    MagicSource.SetUpReso(1, 0.476);
    MagicSource.InitNumMcIter(1000000);

    // the cut off scale in k*, for which the angular distributions from EPOS
    // are evaluated. 200 MeV works okay, you can go up to 300 MeV for systematic checks
    const double k_CutOff = 300;

    // to be used for the NTuple later on
    Float_t k_D;
    Float_t fP1;
    Float_t fP2;
    Float_t fM1;
    Float_t fM2;
    Float_t Tau1;
    Float_t Tau2;
    Float_t AngleRcP1;
    Float_t AngleRcP2;
    Float_t AngleP1P2;
    // random generator dimi style. The input is incompatible with the ROOT random generator,
    // do not mix and match, do not ask me how I know this. Ask Bernie.
    // 11 is the seed, you can change that to you favorite number
    DLM_Random RanGen(11);
    // dummies to save random shit
    double RanVal1;
    double RanVal2;
    double RanVal3;
    double RanCos;
    double MeanP2 = 0;
    // printf("1\n");
    // open the magic file from dimi with the angular distributions.
    TFile *F_EposDisto =
        new TFile("/Users/emma/Femtoscopy/DMesonLight/oton_selections/radius/files/ALL_D_KaonReso.root");
    // set up the ntuple, do not change anything unless told so by dimi
    TNtuple *T_EposDisto = (TNtuple *)F_EposDisto->Get("InfoTuple_ClosePairs");
    unsigned N_EposDisto = T_EposDisto->GetEntries();
    T_EposDisto->SetBranchAddress("k_D", &k_D);
    T_EposDisto->SetBranchAddress("P1", &fP1);
    T_EposDisto->SetBranchAddress("P2", &fP2);
    T_EposDisto->SetBranchAddress("M1", &fM1);
    T_EposDisto->SetBranchAddress("M2", &fM2);
    T_EposDisto->SetBranchAddress("Tau1", &Tau1);
    T_EposDisto->SetBranchAddress("Tau2", &Tau2);
    T_EposDisto->SetBranchAddress("AngleRcP1", &AngleRcP1);
    T_EposDisto->SetBranchAddress("AngleRcP2", &AngleRcP2);
    T_EposDisto->SetBranchAddress("AngleP1P2", &AngleP1P2);

    // TH1F* hAngleRcP1 = new TH1F("hAngleRcP1","hAngleRcP1",1024,0,TMath::Pi());

    gROOT->cd();
    TH1F *hAngle = new TH1F("hAngle", "hAngle", 32, 0, TMath::Pi());
    TH1F *hCos = new TH1F("hCos", "hCos", 32, -1., 1.);
    TH1F *hFinalAngle = new TH1F("hFinalAngle", "hFinalAngle", 32, 0, TMath::Pi());
    TH1F *hFinalCos = new TH1F("hFinalCos", "hFinalCos", 32, -1., 1.);
    F_EposDisto->cd();
    int NumUsefulEntries = 0;
    for (unsigned uEntry = 0; uEntry < N_EposDisto; uEntry++)
    {
        T_EposDisto->GetEntry(uEntry);
        if (k_D > k_CutOff)
            continue;
        hAngle->Fill(AngleRcP2);
        hCos->Fill(cos(AngleRcP2));
        MeanP2 += fP2;
        NumUsefulEntries++;
    }
    MeanP2 /= double(NumUsefulEntries);

    // hAngle->SetBinContent(hAngle->FindBin(0.275),0.31);
    hCos->Scale(1. / hAngle->Integral(), "width");
    // printf("2\n");

    if (!SmoothSampling)
    {
        // iterate over the ntuple
        for (unsigned uEntry = 0; uEntry < N_EposDisto; uEntry++)
        {
            // get each entry
            T_EposDisto->GetEntry(uEntry);
            // disregard the entry of you are outside the desired k*
            if (k_D > k_CutOff)
                continue;
            // overwrite the value for the lifetime. This is computed from the
            // stat. hadronization model (Vale) or thermal fist (Max)
            // this is the value for the secondary protons
            Tau1 = 0;
            // for primoridials (the Xis) we put 0
            Tau2 = 3.66;
            // put in the average mass of the resonances (again from SHM or TF)
            // this is the value for protons
            fM2 = 1054;
            // generate a random path length for the propagation of the resonances
            // nothing to change!
            RanVal2 = RanGen.Exponential(fM2 / (fP2 * Tau2));
            // adds a single entry into the PDF for the angular distribution to be used
            if (!IncludeBug)
                MagicSource.AddBGT_PR(RanVal2, cos(AngleRcP2));
            else
                MagicSource.AddBGT_PR(RanVal2, -cos(AngleRcP2));
            hFinalAngle->Fill(AngleRcP2);
            hFinalCos->Fill(cos(AngleRcP2));
        }
    }
    else
        cout << "LOL KAON" << endl;
    // else if(SmoothSampling==1){
    //   gROOT->cd();
    //   TF1* fRandomAngle = MakeSmoothAngularSourceDisto(hAngle);
    //   //TRandom3* ran_gen = new TRandom3(11);
    //   for(unsigned uEntry=0; uEntry<SmoothEntires; uEntry++){
    //       Tau1 = 0;
    //       Tau2 = 3.66;
    //       fM2 = 1054;
    //       fP2 = MeanP2;
    //       RanVal2 = RanGen.Exponential(fM2/(fP2*Tau2));
    //       AngleRcP2 = fRandomAngle->GetRandom(0.,TMath::Pi());
    //       if(!IncludeBug) MagicSource.AddBGT_PR(RanVal2,cos(AngleRcP2));
    //       else MagicSource.AddBGT_PR(RanVal2,-cos(AngleRcP2));
    //       hFinalAngle->Fill(AngleRcP2);
    //       hFinalCos->Fill(cos(AngleRcP2));
    //   }
    //   gROOT->cd();
    //   //delete ran_gen;
    //   delete fRandomAngle;
    // }
    // else{
    //   gROOT->cd();
    //   TF1* fRandomCos = MakeSmoothAngularCosine(hCos);
    //   //TRandom3* ran_gen = new TRandom3(11);
    //   for(unsigned uEntry=0; uEntry<SmoothEntires; uEntry++){
    //       Tau1 = 0;
    //       Tau2 = 3.66;
    //       fM2 = 1054;
    //       fP2 = MeanP2;
    //       RanVal2 = RanGen.Exponential(fM2/(fP2*Tau2));
    //       RanCos = fRandomCos->GetRandom(-1.,1.);
    //       if(!IncludeBug) MagicSource.AddBGT_PR(RanVal2,RanCos);
    //       else MagicSource.AddBGT_PR(RanVal2,RanCos);
    //       hFinalAngle->Fill(acos(RanCos));
    //       hFinalCos->Fill(RanCos);
    //   }
    //   gROOT->cd();
    //   delete fRandomCos;
    // }

    delete F_EposDisto;

    // if you have resonances contributing to both particles, we need to repeat the above procedure
    // for the prim-reso (AddBGT_PR) and reso-reso (AddBGT_RR) cases

    const unsigned NumSourceBins = 128;
    const double rMin = 0;
    const double rMax = 16;
    TFile *fOutput = new TFile(TString::Format("fSecond_KaonD_%.2f.root", CoreSize), "recreate");
    TH1F *hSource = new TH1F("hSource", "hSource", NumSourceBins, rMin, rMax);

    // fill the histo fro the source
    for (unsigned uBin = 0; uBin < NumSourceBins; uBin++)
    {
        // get the x-axis (r value) of the current bin
        double xaxis = hSource->GetBinCenter(uBin + 1);
        // an array for the parameters, [0] is source size, [1] is == 2 (for a Gaussian)
        double parameters[2];
        parameters[0] = CoreSize;
        parameters[1] = 2.0;
        double SourceValue = MagicSource.RootEval(&xaxis, parameters);
        hSource->SetBinContent(uBin + 1, SourceValue);
        // infinite errors for now
        hSource->SetBinError(uBin + 1, 1000.);
    }
    // printf("4\n");
    // idea: fit the source distribution only in a range around its peak
    // to do this: silly idea: put very large uncertainties in the bins outside of this range
    // we can get this range automatically, by evaluating the central (median) integral of the source distribution
    // with this set up, we fit the 68% most central yield of the source distribution
    double lowerlimit;
    double upperlimit;
    GetCentralInterval(*hSource, 0.84, lowerlimit, upperlimit, true);
    unsigned lowerbin = hSource->FindBin(lowerlimit);
    unsigned upperbin = hSource->FindBin(upperlimit);
    for (unsigned uBin = lowerbin; uBin <= upperbin; uBin++)
    {
        hSource->SetBinError(uBin + 1, 0.01);
    }

    printf("Core size of %.3f fm\n", CoreSize);
    printf("The fit will be performed in the range [%.2f, %.2f] fm\n", lowerlimit, upperlimit);
    // fyi, GaussSourceTF1 is in DLM_Source.h if you want to check it out.
    TF1 *fSource = new TF1("fSource",
                           "[0] * 4 * TMath::Pi() * x ** 2 * 1. / ((4 * TMath::Pi() *[1] ** 2) ** 1.5) * TMath::E() ** "
                           "(- x ** 2 / 4 /[1] ** 2) + (1 - [0]) * 4 * TMath::Pi() * x ** 2 * 1. / ((4 * TMath::Pi() * "
                           "[2] ** 2) ** 1.5) * TMath::E() ** (- x ** 2 / 4 / [2] ** 2)",
                           rMin, rMax);
    fSource->SetParameter(0, 0.5);
    fSource->SetParLimits(0, 0.4, 1);
    fSource->SetParameter(1, CoreSize);
    // fSource->SetParLimits(1,CoreSize*0.5,CoreSize*2.0);
    fSource->SetParLimits(1, 0, 10);
    fSource->SetParLimits(2, fSource->GetParameter(1), 10);

    hSource->Fit(fSource, "Q, S, N, R, M");
    printf("The weight size is %.3f +/- %.3f fm\n", fSource->GetParameter(0), fSource->GetParError(0));
    printf("The effective Gaussian size 1 is %.3f +/- %.3f fm\n", fSource->GetParameter(1), fSource->GetParError(1));
    printf("The effective Gaussian size 2 is %.3f +/- %.3f fm\n", fSource->GetParameter(2), fSource->GetParError(2));

    *weight = fSource->GetParameter(0);
    *weight_error = fSource->GetParError(0);

    *r1 = fSource->GetParameter(1);
    *r1_error = fSource->GetParError(1);

    *r2 = fSource->GetParameter(2);
    *r2_error = fSource->GetParError(2);

    *chi2ndf = fSource->GetChisquare() / fSource->GetNDF();

    printf("The chi2/ndf is %.2f fm\n", fSource->GetChisquare() / fSource->GetNDF());

    TF1 *fSourceLambda = new TF1("fSourceLambda",
                                 "[3]* ( [0] * 4 * TMath::Pi() * x ** 2 * 1. / ((4 * TMath::Pi() *[1] ** 2) ** 1.5) * "
                                 "TMath::E() ** (- x ** 2 / 4 /[1] ** 2) + (1 - [0]) * 4 * TMath::Pi() * x ** 2 * 1. / "
                                 "((4 * TMath::Pi() * [2] ** 2) ** 1.5) * TMath::E() ** (- x ** 2 / 4 / [2] ** 2))",
                                 rMin, rMax);
    fSourceLambda->SetParameter(0, 0.5);
    fSourceLambda->SetParLimits(0, 0, 1);
    fSourceLambda->SetParameter(1, CoreSize);
    // fSourceLambda->SetParLimits(1,CoreSize*0.5,CoreSize*2.0);
    fSourceLambda->SetParLimits(1, 0, 10);
    fSourceLambda->SetParLimits(2, fSourceLambda->GetParameter(1), 10);
    fSourceLambda->SetParLimits(3, 0, 1);

    hSource->Fit(fSourceLambda, "Q, S, N, R, M");

    // get rid of weird plotting
    for (unsigned uBin = 0; uBin < NumSourceBins; uBin++)
    {
        hSource->SetBinError(uBin + 1, 0.01);
    }
    hSource->Write();
    fSource->Write();
    fSourceLambda->Write();
    hAngle->Scale(1. / hAngle->Integral(), "width");
    hAngle->Write();
    hCos->Write();
    hFinalAngle->Scale(1. / hFinalAngle->Integral(), "width");
    hFinalAngle->Write();
    hFinalCos->Scale(1. / hFinalCos->Integral(), "width");
    hFinalCos->Write();

    // gROOT->cd();
    // delete hAngle;

    delete hSource;
    delete fSource;
    delete hAngle;
    delete hCos;
    delete hFinalAngle;
    delete hFinalCos;
    delete fOutput;

    printf("-------------------\n");
}










int PPHISOURCERAD(int argc, char *argv[]){
    printf("Suck my cock\n");
    DLM_Timer TIMER;

    //EffectiveGaussianpphi();

    //EmmaDaniel_piD(const double CoreSize, double *weight, double *r1, double *r2, double *chi2ndf)
    double weight[10], r1[10], r2[10], chi2ndf[10];
    //EmmaDaniel_piD(0.82, weight, r1, r2, chi2ndf,0,16);
    //printf("AS D. DID FOR PAPER\n");
    //EmmaDaniel_piD(0.82, weight, r1, r2, chi2ndf,0,4);
    //EmmaDaniel_piD(0.76, weight, r1, r2, chi2ndf,0,16);
    //EmmaDaniel_piD(0.76, weight, r1, r2, chi2ndf,0,4);
    EmmaDaniel_KD(0.81, weight, r1, r2, chi2ndf, 0, 16);
    EmmaDaniel_KD(0.75, weight, r1, r2, chi2ndf, 0, 16);
    printf("w=%f; r1=%f; r2=%f; chi2ndf=%f\n", weight[0],r1[0],r2[0],chi2ndf[0]);

    long long ExeTime = TIMER.Stop()/1000.;
    char* strtime = new char [128];
    ShowTime(ExeTime,strtime,0,true,6);
    printf("The script terminated after: %s\n",strtime);
    delete [] strtime;
    return 0;
}
