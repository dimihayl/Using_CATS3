#include "Smearing.h"
#include "EnvVars.h"
#include "CommonAnaFunctions.h"


#include "omp.h"

#include "DLM_Random.h"
#include "DLM_Histo.h"
#include "DLM_RootWrapper.h"
#include "CATS.h"
#include "DLM_Ck.h"
#include "DLM_CkDecomposition.h"
#include "DLM_CkDecomp.h"
#include "DLM_Unfold.h"
#include "DLM_MultiFit.h"

#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TFitResult.h"

using namespace std;

//we work for pLambda for now (ME)
void FitSmearMatrix(){
  //for the smear matrix
  const TString InputFileName = TString::Format("%s/CatsFiles/MomentumSmear/ALICE_pp_13TeV_MEpL.root",GetCernBoxDimi());
  const TString OutputFolder = TString::Format("%s/Smearing/FitSmearMatrix/",GetFemtoOutputFolder());
  const int REBIN = 4;
  TFile fInput(InputFileName,"read");
  TH2F* hMomSmear_pL_ME = (TH2F*)fInput.Get("h_RESO_pL_MeV");
  hMomSmear_pL_ME->Rebin2D(REBIN,REBIN);
  const unsigned NbinsX = hMomSmear_pL_ME->GetNbinsX();
  const unsigned NbinsY = hMomSmear_pL_ME->GetNbinsY();
  const double kMin = hMomSmear_pL_ME->GetXaxis()->GetBinLowEdge(1);
  const double kMax = hMomSmear_pL_ME->GetXaxis()->GetBinUpEdge(NbinsX);
  const unsigned NbinsX_MC = 1500;
  const double kMin_MC = 0;
  const double kMax_MC = 750;
  const unsigned NbinsY_MC = 2000;
  const double deltaY = 100;

  //compare to some basic MC
  const long unsigned NumIter = (long unsigned)(10*1000)*(long unsigned)(1000*1000);
  const double Reso = 4.7;

  //TH2F* hMomSmear_Rotated = (TH2F*)hMomSmear_pL_ME->Clone("histoCopy");
  TH2F* hMomSmear_Rotated = new TH2F("hMomSmear_Rotated","hMomSmear_Rotated",
    NbinsX,kMin,kMax,NbinsY,-kMax*0.5,kMax*0.5);
  for(unsigned uBinX=0; uBinX<NbinsX; uBinX++){
    for(unsigned uBinY=0; uBinY<NbinsY; uBinY++){
      hMomSmear_Rotated->SetBinContent(uBinX+1,uBinY+1,0);
    }
  }
  DLM_Histo<float> dlmMomSmear_Rotated_MC; //= Convert_TH2F_DlmHisto(hMomSmear_Rotated);
  dlmMomSmear_Rotated_MC.SetUp(2);
  dlmMomSmear_Rotated_MC.SetUp(0,NbinsX_MC,kMin_MC,kMax_MC);
  dlmMomSmear_Rotated_MC.SetUp(1,NbinsY_MC,-deltaY,deltaY);
  dlmMomSmear_Rotated_MC.Initialize();

  DLM_Histo<float> dlmMomSmear_MC; //= Convert_TH2F_DlmHisto(hMomSmear_Rotated);
  dlmMomSmear_MC.SetUp(2);
  dlmMomSmear_MC.SetUp(0,NbinsX_MC,kMin_MC,kMax_MC);
  dlmMomSmear_MC.SetUp(1,NbinsX_MC,kMin_MC,kMax_MC);
  dlmMomSmear_MC.Initialize();

  double kx_val,ky_val,delta_k,old_val;
  int WhichBin;
  for(unsigned uBinX=0; uBinX<NbinsX; uBinX++){
    kx_val = hMomSmear_pL_ME->GetXaxis()->GetBinCenter(uBinX+1);
    for(unsigned uBinY=0; uBinY<NbinsY; uBinY++){
      ky_val = hMomSmear_pL_ME->GetYaxis()->GetBinCenter(uBinY+1);
      delta_k = kx_val-ky_val;
      WhichBin = hMomSmear_Rotated->FindBin(kx_val,delta_k);
      old_val = hMomSmear_Rotated->GetBinContent(WhichBin);
      hMomSmear_Rotated->SetBinContent(WhichBin,old_val+hMomSmear_pL_ME->GetBinContent(uBinX+1,uBinY+1));
    }
  }



  unsigned NumThreads = omp_get_num_procs();
  printf("NumThreads = %u\n",NumThreads);
  DLM_Random** dlmrangen;
  dlmrangen = new DLM_Random* [NumThreads];
  for(unsigned uThr=0; uThr<NumThreads; uThr++){
    dlmrangen[uThr] = new DLM_Random(uThr+1);
  }
  long unsigned omp_u;
  #pragma omp parallel for
  for(omp_u=0; omp_u<NumIter; omp_u++){
    unsigned thr = omp_get_thread_num();

    //double xMom = dlmrangen[thr]->Uniform(-0.5*kMax,0.5*kMax);
    //double yMom = dlmrangen[thr]->Uniform(-0.5*kMax,0.5*kMax);
    //double zMom = dlmrangen[thr]->Uniform(-0.5*kMax,0.5*kMax);
    //double xMom = dlmrangen[thr]->Gauss(0,0.25*kMax);
    //double yMom = dlmrangen[thr]->Gauss(0,0.25*kMax);
    //double zMom = dlmrangen[thr]->Gauss(0,0.25*kMax);
    double xMom = dlmrangen[thr]->Exponential(1./(0.12*kMax_MC));
    double yMom = dlmrangen[thr]->Exponential(1./(0.12*kMax_MC));
    double zMom = dlmrangen[thr]->Exponential(1./(0.12*kMax_MC));


    double Mom = sqrt(xMom*xMom+yMom*yMom+zMom*zMom);
    double xSmr = dlmrangen[thr]->Gauss(xMom,Reso);
    double ySmr = dlmrangen[thr]->Gauss(yMom,Reso);
    double zSmr = dlmrangen[thr]->Gauss(zMom,Reso);
    double MomSmr = sqrt(xSmr*xSmr+ySmr*ySmr+zSmr*zSmr);
    double axis_val[2];
    axis_val[0] = Mom;
    axis_val[1] = Mom-MomSmr;
    double axis_mom[2];
    axis_mom[0] = Mom;
    axis_mom[1] = MomSmr;
    //#pragma omp critical
    {
    dlmMomSmear_Rotated_MC.AddAt(axis_val);
    dlmMomSmear_MC.AddAt(axis_mom);
    }

  }
  TH2F* hMomSmear_Rotated_MC = Convert_DlmHisto_TH2F(&dlmMomSmear_Rotated_MC,"hMomSmear_Rotated_MC");
  TH2F* hMomSmear_MC = Convert_DlmHisto_TH2F(&dlmMomSmear_MC,"hMomSmear_MC");

  TFile fOutput(OutputFolder+"fOutput.root","recreate");
  hMomSmear_pL_ME->Write();
  hMomSmear_Rotated->Write();
  hMomSmear_Rotated_MC->Write();
  hMomSmear_MC->Write();

  //delete hMomSmear_pL_ME;
  delete hMomSmear_Rotated;
  delete hMomSmear_Rotated_MC;
  delete hMomSmear_MC;
  for(unsigned uThr=0; uThr<NumThreads; uThr++){
    delete dlmrangen[uThr];
  }
  delete [] dlmrangen;

}


//par 0 NORM
//par 1 WIDTH
//par 2 Mass1
//par 3 Mass2
double TrueMixedEvent(double* x, double* par){
  double& xVal = *x;
  return par[0]*xVal*xVal*exp(-pow(xVal/par[1],2.))/sqrt(par[2]*par[2]+xVal*xVal)/sqrt(par[3]*par[3]+xVal*xVal);
}

//SystemName = pp, pLambda, pXim
void FiniteBinningStudy(const TString SystemName, const double SourceSize){
  const double kMin = 0;
  const double kMax = 384;
  const TString OutputFolder = TString::Format("%s/Smearing/FiniteBinningStudy/",GetFemtoOutputFolder());

/*
  const unsigned NumBinnings = 7;
  unsigned* NumBins = new unsigned[NumBinnings];
  NumBins[0] = 768;//0.5 MeV
  NumBins[1] = 384;//1 MeV
  NumBins[2] = 192;//2 MeV
  NumBins[3] = 96;//4 MeV
  NumBins[4] = 48;//8 MeV
  NumBins[5] = 24;//16 MeV
  NumBins[6] = 12;//32 MeV
*/
  const unsigned NumBinnings = 6;
  unsigned* NumBins = new unsigned[NumBinnings];
  NumBins[0] = 384;//1 MeV
  NumBins[1] = 192;//2 MeV
  NumBins[2] = 96;//4 MeV
  NumBins[3] = 48;//8 MeV
  NumBins[4] = 24;//16 MeV
  NumBins[5] = 12;//32 MeV

  const unsigned NumMomSmear = 2;
  TString* MomSmearMatrixFile = new TString[NumMomSmear];
  TString* MomSmearMatrixHisto = new TString[NumMomSmear];

  MomSmearMatrixFile[0] = TString::Format("%s/Smearing/FitSmearMatrix/fVeryFine10B_4.7MeV.root",GetFemtoOutputFolder());
  MomSmearMatrixHisto[0] = "hMomSmear_MC";

  MomSmearMatrixFile[1] = TString::Format("%s/CatsFiles/MomentumSmear/ALICE_pp_13TeV_MEpL.root",GetCernBoxDimi());
  MomSmearMatrixHisto[1] = "h_RESO_pL_MeV";

  double PSpar[4];
  unsigned DesiredPairsBelow32;
  if(SystemName=="pp"){
    PSpar[0] = 1;
    PSpar[1] = 1000;
    PSpar[2] = 938;
    PSpar[3] = 938;
    DesiredPairsBelow32 = 2000;
  }
  else if(SystemName=="pLambda"){
    PSpar[0] = 1;
    PSpar[1] = 1000;
    PSpar[2] = 938;
    PSpar[3] = 1116;
    DesiredPairsBelow32 = 500;
  }
  else if(SystemName=="pXim"){
    PSpar[0] = 1;
    PSpar[1] = 1000;
    PSpar[2] = 938;
    PSpar[3] = 1320;
    DesiredPairsBelow32 = 100;
  }

  DLM_Histo<float>* dlmPhaseSpace = new DLM_Histo<float>[NumBinnings];
  for(unsigned uBW=0; uBW<NumBinnings; uBW++){
    dlmPhaseSpace[uBW].SetUp(1);
    dlmPhaseSpace[uBW].SetUp(0,NumBins[uBW],kMin,kMax);
    dlmPhaseSpace[uBW].Initialize();
    for(unsigned uBin=0; uBin<NumBins[uBW]; uBin++){
      double MOM = dlmPhaseSpace[uBW].GetBinCenter(0,uBin);
      dlmPhaseSpace[uBW].SetBinContent(uBin,TrueMixedEvent(&MOM,PSpar));
    }
    double k32=32;
    double Rescale = double(DesiredPairsBelow32)/dlmPhaseSpace[uBW].Integral(&kMin,&k32);
    dlmPhaseSpace[uBW] *= Rescale;
  }


  for(unsigned uSmear=0; uSmear<NumMomSmear; uSmear++){
    TFile fMatrixInput(MomSmearMatrixFile[uSmear],"read");
    TH2F* hMomSmear = (TH2F*)fMatrixInput.Get(MomSmearMatrixHisto[uSmear]);
    DLM_Histo<float>* dlmMomSmear;
    dlmMomSmear = new DLM_Histo<float>[NumBinnings];
    DLM_CommonAnaFunctions* AnalysisObject = new DLM_CommonAnaFunctions[NumBinnings];
    CATS* Kitty = new CATS[NumBinnings];
    double lambda[5];
    for(unsigned uBW=0; uBW<NumBinnings; uBW++){
      AnalysisObject[uBW].SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());
      dlmMomSmear[uBW].SetUp(2);
      dlmMomSmear[uBW].SetUp(0,NumBins[uBW],kMin,kMax);
      dlmMomSmear[uBW].SetUp(1,NumBins[uBW],kMin,kMax);
      dlmMomSmear[uBW].Initialize();
      dlmMomSmear[uBW].SetBinContentAll(0);
      unsigned uBin[2];
      for(uBin[0]=0; uBin[0]<NumBins[uBW]; uBin[0]++){
        double minX = dlmMomSmear[uBW].GetBinLowEdge(0,uBin[0]);
        double maxX = dlmMomSmear[uBW].GetBinUpEdge(0,uBin[0]);
        for(uBin[1]=0; uBin[1]<NumBins[uBW]; uBin[1]++){
          double minY = dlmMomSmear[uBW].GetBinLowEdge(1,uBin[1])*(1+1e-6);//avoid silly stuff at the bin edge, like getting into the next bin
          double maxY = dlmMomSmear[uBW].GetBinUpEdge(1,uBin[1])*(1-1e-6);
          unsigned FirstBinX = hMomSmear->GetXaxis()->FindBin(minX);
          unsigned LastBinX = hMomSmear->GetXaxis()->FindBin(maxX);
          unsigned FirstBinY = hMomSmear->GetYaxis()->FindBin(minY);
          unsigned LastBinY = hMomSmear->GetYaxis()->FindBin(maxY);
          dlmMomSmear[uBW].SetBinContent(uBin,hMomSmear->Integral(FirstBinX,LastBinX,FirstBinY,LastBinY));
        }
      }

      Kitty[uBW].SetMomBins(NumBins[uBW],kMin,kMax);
      if(SystemName=="pp")            {AnalysisObject[uBW].SetUpCats_pp(Kitty[uBW],"AV18","Gauss",0,202);
                                      if(!uBW){AnalysisObject[uBW].SetUpLambdaPars_pp("pp13TeV_HM_Dec19",0,lambda);}}
      else if(SystemName=="pLambda")  {AnalysisObject[uBW].SetUpCats_pL(Kitty[uBW],"Chiral_Coupled_SPD","Gauss",11600,202);
                                      if(!uBW){AnalysisObject[uBW].SetUpLambdaPars_pL("pp13TeV_HM_Dec19",0,0,lambda);}}
      else if(SystemName=="pXim")     {AnalysisObject[uBW].SetUpCats_pXim(Kitty[uBW],"pXim_HALQCD1","Gauss");
                                      if(!uBW){AnalysisObject[uBW].SetUpLambdaPars_pXim("pp13TeV_HM_Dec19",0,0,lambda);}}
      else {printf("The system %s is not available\n",SystemName.Data()); continue;}
      Kitty[uBW].SetAnaSource(0,SourceSize);
      Kitty[uBW].SetNotifications(CATS::nWarning);
      //Kitty[uBW].SetEpsilonConv(4e-8);
      //Kitty[uBW].SetEpsilonProp(4e-8);
    }//uBW

    DLM_Histo<float>* dlmPhaseSpaceExp = new DLM_Histo<float>[NumBinnings];

    unsigned omp_ubw;
    #pragma omp parallel for
    for(omp_ubw=0; omp_ubw<NumBinnings; omp_ubw++){
      DLM_Unfold dlmFold;
      dlmFold.SetResponse(&dlmMomSmear[omp_ubw]);
      dlmFold.SetSilentMode(false);
      dlmFold.SetData(&dlmPhaseSpace[omp_ubw]);
      DLM_Histo<float>* temp = dlmFold.Fold();
      dlmPhaseSpaceExp[omp_ubw] = *temp;
      delete temp;
    }

    #pragma omp parallel for
    for(omp_ubw=0; omp_ubw<NumBinnings; omp_ubw++){
      Kitty[omp_ubw].KillTheCat();
    }

    DLM_Ck** CkHisto = new DLM_Ck* [NumBinnings];
    DLM_CkDecomp** CkDecDirect = new DLM_CkDecomp* [NumBinnings];
    DLM_CkDecomp** CkDecPS = new DLM_CkDecomp* [NumBinnings];
    #pragma parallel for
    for(omp_ubw=0; omp_ubw<NumBinnings; omp_ubw++){
      CkHisto[omp_ubw] = new DLM_Ck(Kitty[omp_ubw].GetNumSourcePars(),0,Kitty[omp_ubw]);
      CkHisto[omp_ubw]->Update();
      CkDecDirect[omp_ubw] = new DLM_CkDecomp("Main",1,*CkHisto[omp_ubw],dlmMomSmear);
      CkDecDirect[omp_ubw]->AddContribution(0,1.-lambda[0],DLM_CkDecomp::cFeedDown);
      CkDecDirect[omp_ubw]->Update();
      CkDecPS[omp_ubw] = new DLM_CkDecomp("Main",1,*CkHisto[omp_ubw],dlmMomSmear);
      CkDecPS[omp_ubw]->AddContribution(0,1.-lambda[0],DLM_CkDecomp::cFeedDown);
      CkDecPS[omp_ubw]->AddPhaseSpace(&dlmPhaseSpace[omp_ubw]);
      CkDecPS[omp_ubw]->Update();
    }//omp_ubw

    TFile fOutput(OutputFolder+TString::Format("%s_%.0f_SM%u.root",SystemName.Data(),SourceSize*100.,uSmear),"recreate");
    //the true correlation
    TH1F** hCk_true = new TH1F* [NumBinnings];
    //current procedure
    TH1F** hCk_direct = new TH1F* [NumBinnings];
    //including phase space correction
    TH1F** hCk_ps = new TH1F* [NumBinnings];
    //including ps correction and correction for the bin width
    TH1F** hCk_ps_bin = new TH1F* [NumBinnings];
    TH1F** hPhaseSpace = new TH1F* [NumBinnings];
    TH1F** hPhaseSpaceExp = new TH1F* [NumBinnings];

    for(unsigned uBW=0; uBW<NumBinnings; uBW++){
      hCk_true[uBW] = new TH1F(TString::Format("hCk_true_%.1f",kMax/double(NumBins[uBW])),TString::Format("hCk_true_%.1f",kMax/double(NumBins[uBW])),NumBins[uBW],kMin,kMax);
      hCk_direct[uBW] = new TH1F(TString::Format("hCk_direct_%.1f",kMax/double(NumBins[uBW])),TString::Format("hCk_direct_%.1f",kMax/double(NumBins[uBW])),NumBins[uBW],kMin,kMax);
      hCk_ps[uBW] = new TH1F(TString::Format("hCk_ps_%.1f",kMax/double(NumBins[uBW])),TString::Format("hCk_ps_%.1f",kMax/double(NumBins[uBW])),NumBins[uBW],kMin,kMax);
      hCk_ps_bin[uBW] = new TH1F(TString::Format("hCk_ps_bin_%.1f",kMax/double(NumBins[uBW])),TString::Format("hCk_ps_bin_%.1f",kMax/double(NumBins[uBW])),NumBins[uBW],kMin,kMax);

      //hPhaseSpace[uBW] = new TH1F(TString::Format("hPhaseSpace_%.1f",kMax/double(NumBins[uBW])),TString::Format("hPhaseSpace_%.1f",kMax/double(NumBins[uBW])),NumBins[uBW],kMin,kMax);
      //hPhaseSpaceExp[uBW] = new TH1F(TString::Format("hPhaseSpaceExp_%.1f",kMax/double(NumBins[uBW])),TString::Format("hPhaseSpaceExp_%.1f",kMax/double(NumBins[uBW])),NumBins[uBW],kMin,kMax);
      hPhaseSpace[uBW] = Convert_DlmHisto_TH1F(&dlmPhaseSpace[uBW],TString::Format("hPhaseSpace_%.1f",kMax/double(NumBins[uBW])));
      hPhaseSpaceExp[uBW] = Convert_DlmHisto_TH1F(&dlmPhaseSpaceExp[uBW],TString::Format("hPhaseSpaceExp_%.1f",kMax/double(NumBins[uBW])));

      //uBinDaughters (Y)
      for(unsigned uBin=0; uBin<NumBins[uBW]; uBin++){
        double MOM = Kitty[uBW].GetMomentum(uBin);
        hCk_true[uBW]->SetBinContent(uBin+1,lambda[0]*Kitty[uBW].GetCorrFun(uBin)+1.-lambda[0]);
        hCk_direct[uBW]->SetBinContent(uBin+1,CkDecDirect[uBW]->EvalCk(MOM));
        hCk_ps[uBW]->SetBinContent(uBin+1,CkDecPS[uBW]->EvalCk(MOM));
      }//uBin


      //(Y)
      for(unsigned uBin=0; uBin<NumBins[uBW]; uBin++){
        //perhaps in parallel
        double TotIntegral=0;
        double NormIntegral=0;
        //(X)
        for(unsigned uBinTrue=0; uBinTrue<NumBins[uBW]; uBinTrue++){
          //iter over the bin itself
          unsigned FineBinSteps = NumBins[0]/NumBins[uBW];
          double t,mp,M,Cp;
          unsigned whichfinebin[2];
          for(unsigned udX=0; udX<FineBinSteps; udX++){
            for(unsigned udY=0; udY<FineBinSteps; udY++){
              whichfinebin[0] = uBinTrue*FineBinSteps+udX;
              whichfinebin[1] = uBin*FineBinSteps+udY;
              t = dlmMomSmear[0].GetBinContent(whichfinebin);
              mp = dlmPhaseSpace[0].GetBinContent(whichfinebin[0]);
              M = dlmPhaseSpace[0].GetBinContent(whichfinebin[1]);
              Cp = hCk_true[0]->GetBinContent(whichfinebin[0]+1);
              TotIntegral += fabs(t*mp/M*Cp);
              NormIntegral += fabs(t*mp/M);
            }
            //TotIntegral += dlmMomSmear[uBW].GetBinContent();
          }
        }
        TotIntegral /= NormIntegral;
        hCk_ps_bin[uBW]->SetBinContent(uBin+1,TotIntegral);
      }





      fOutput.cd();
      hCk_true[uBW]->Write();
      hCk_direct[uBW]->Write();
      hCk_ps[uBW]->Write();
      hCk_ps_bin[uBW]->Write();
      hPhaseSpace[uBW]->Write();
      hPhaseSpaceExp[uBW]->Write();
    }//uBW

    delete [] Kitty;
    delete [] AnalysisObject;
    delete [] dlmMomSmear;
    delete [] dlmPhaseSpaceExp;

    for(unsigned uBW=0; uBW<NumBinnings; uBW++){
      delete hCk_true[uBW];
      delete hCk_direct[uBW];
      delete hCk_ps[uBW];
      delete hCk_ps_bin[uBW];
      delete hPhaseSpace[uBW];
      delete hPhaseSpaceExp[uBW];
      delete CkHisto[uBW];
      delete CkDecDirect[uBW];
      delete CkDecPS[uBW];
    }
    delete [] hCk_true;
    delete [] hCk_direct;
    delete [] hCk_ps;
    delete [] hCk_ps_bin;
    delete [] hPhaseSpace;
    delete [] hPhaseSpaceExp;
    delete [] CkHisto;
    delete [] CkDecDirect;
    delete [] CkDecPS;

    fOutput.Close();
//break;
  }//uSmear

  delete [] MomSmearMatrixFile;
  delete [] MomSmearMatrixHisto;
  delete [] NumBins;
  delete [] dlmPhaseSpace;
}

double MSMP1(double*x,double*par){
  double& MOM = *x;
  double& NG = par[0];
  double& meanG = par[1];
  double& sigG = par[2];
  double& limit = par[3];
  double& NB = par[4];
  double& meanB = par[5];
  double& sigB = par[6];
  //below 3sigma -> gauss
  //after that: cauchy
  if(fabs(MOM-meanG)/meanG<limit){
    return 0;
  }
  else{
    return 0;
  }

}

void ParametrizeSmearMatrix_MC(){
  const TString OutputFolder = TString::Format("%s/Smearing/ParametrizeSmearMatrix_MC/",GetFemtoOutputFolder());
  const double kMin=0;
  const double kMax=120;
  const double expected_sigma = 5;
  const unsigned NumMomSmear = 2;
  TString* MomSmearMatrixFile = new TString[NumMomSmear];
  TString* MomSmearMatrixHisto = new TString[NumMomSmear];

  MomSmearMatrixFile[0] = TString::Format("%s/Smearing/FitSmearMatrix/fVeryFine10B_5MeV.root",GetFemtoOutputFolder());
  MomSmearMatrixHisto[0] = "hMomSmear_MC";

  MomSmearMatrixFile[1] = TString::Format("%s/CatsFiles/MomentumSmear/ALICE_pp_13TeV_MEpL.root",GetCernBoxDimi());
  MomSmearMatrixHisto[1] = "h_RESO_pL_MeV";

  for(unsigned uSmear=0; uSmear<NumMomSmear; uSmear++){
    TFile fMatrixInput(MomSmearMatrixFile[uSmear],"read");
    TH2F* hMomSmear = (TH2F*)fMatrixInput.Get(MomSmearMatrixHisto[uSmear]);
    const unsigned NumBinsY = hMomSmear->GetYaxis()->GetNbins();
    double kstar,fit_low,fit_up;
    DLM_MultiFit MultiFitter;
    TFile fOutput(OutputFolder+TString::Format("fOutput_%u.root",uSmear),"recreate");
    TF1** fitFunctions = new TF1* [NumBinsY];
    TH1D** histoToFit = new TH1D* [NumBinsY];
    unsigned histoCounter=0;
    for(unsigned uBinY=0; uBinY<NumBinsY; uBinY++){
      fitFunctions[uBinY] = NULL;
      histoToFit[uBinY] = NULL;
      kstar = hMomSmear->GetYaxis()->GetBinCenter(uBinY+1);
      if(kstar<kMin || kstar>kMax) continue;
      fit_low = kstar-8*expected_sigma;
      fit_up = kstar+8*expected_sigma;
      if(fit_low<kMin)fit_low=kMin;
      if(fit_up>kMax)fit_up=kMax;

      TH1D* hProj = hMomSmear->ProjectionX(TString::Format("hProj%u",uBinY),uBinY+1,uBinY+1);
      fOutput.cd();
      histoToFit[uBinY] = (TH1D*)hProj->Clone(TString::Format("histoToFit%u",uBinY));
      //histoToFit[uBinY]->Sumw2();
      for(unsigned uBin=0; uBin<histoToFit[uBinY]->GetNbinsX(); uBin++){
        double BC = histoToFit[uBinY]->GetBinContent(uBin+1);
        double BE = sqrt(BC);
        if(BE==0) BE = 1;
        histoToFit[uBinY]->SetBinError(uBin+1,BE);
      }
      histoToFit[uBinY]->Scale(1./histoToFit[uBinY]->Integral(),"width");
      delete hProj;
      fOutput.cd();
      //fitFunctions[uBinY] = new TF1(TString::Format("fitFunction%u",uBinY),"[0]/([1]*[1])*TMath::Poisson(x,[1]*[1])",kMin,kMax);
      //fitFunctions[uBinY]->FixParameter(0,kstar);
      //fitFunctions[uBinY]->SetParameter(1,expected_sigma);
      //fitFunctions[uBinY]->SetParLimits(1,0.25*expected_sigma,4.*expected_sigma);


      fitFunctions[uBinY] = new TF1(TString::Format("fitFunction%u",uBinY),
      "[0]*TMath::Gaus(x,[1],[2],1)+[3]*TMath::Landau(x,[4],[5],1)+[6]*TMath::Landau(2.*[7]-x,[7],[8],1)",fit_low,fit_up);
      fitFunctions[uBinY]->SetParameter(0,1);
      fitFunctions[uBinY]->SetParameter(1,kstar);
      fitFunctions[uBinY]->SetParameter(2,expected_sigma);
      fitFunctions[uBinY]->SetParameter(3,1);
      fitFunctions[uBinY]->SetParameter(4,kstar*1.5);
      fitFunctions[uBinY]->SetParameter(5,expected_sigma*0.5);
      fitFunctions[uBinY]->SetParameter(6,1);
      fitFunctions[uBinY]->SetParameter(7,kstar/1.5);
      fitFunctions[uBinY]->SetParameter(8,expected_sigma*0.5);




      //MultiFitter.AddSpectrum(histoToFit[uBinY],fitFunctions[uBinY]);
      //if(histoCounter){
      //  MultiFitter.SetEqualPar(histoCounter,1,histoCounter-1,1);
      //}
      histoCounter++;
      histoToFit[uBinY]->Fit(fitFunctions[uBinY],"Q, S, N, R, M");
    }
    //ROOT::Fit::FitResult fitResult = MultiFitter.PerformGlobalFit();


    for(unsigned uBinY=0; uBinY<NumBinsY; uBinY++){
      if(histoToFit[uBinY]){histoToFit[uBinY]->Write();}
      if(fitFunctions[uBinY]){fitFunctions[uBinY]->Write();}

      if(fitFunctions[uBinY]){delete fitFunctions[uBinY];fitFunctions[uBinY]=NULL;}
      if(histoToFit[uBinY]){delete histoToFit[uBinY];histoToFit[uBinY]=NULL;}
    }
    delete [] fitFunctions;
    delete [] histoToFit;
    fOutput.Close();
//break;
  }



}


int SMEARING(int argc, char *argv[]){

  //FitSmearMatrix();
  //FiniteBinningStudy("pp",1.2);
  //FiniteBinningStudy("pLambda",1.2);
  //FiniteBinningStudy("pXim",1.2);
  //FiniteBinningStudy("pXim",1.43);
  ParametrizeSmearMatrix_MC();

  return 0;
}
