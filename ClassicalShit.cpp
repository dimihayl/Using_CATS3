

#include "omp.h"
#include <unistd.h>
#include <thread>


#include "DLM_Random.h"
#include "DLM_Potentials.h"
#include "CommonAnaFunctions.h"
#include "EnvVars.h"
#include "CATS.h"
#include "DLM_Source.h"
#include "DLM_Potentials.h"
#include "CATSconstants.h"
#include "DLM_Histo.h"
#include "DLM_RootWrapper.h"

#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TF1.h"


//test the idea, so get two particles, say proton and Lambda in spin 0, lets see how C(k) would look like
void Toy_MC_1(){
  const int NumIter = 200*1000;
  const double dt = 0.002*FmToNu*2;
  const double tmax = 20*FmToNu;

  //used for the eval of the derivative
  const double dr = 0.001*FmToNu;

  const double SourceSize = 1.0*FmToNu;
  const double red_mass = (Mass_p*Mass_L)/(Mass_p+Mass_L);

  const double kstar_width = 450/2;

  const double Vr_cutoff = 10;

  //do not allow to produce particles at a lower relative distance
  const double source_cutoff = 0.7*FmToNu;

  int NumThreads = std::thread::hardware_concurrency();

  DLM_Random* dlmran = new DLM_Random [NumThreads];
  for(unsigned uTh=0; uTh<NumThreads; uTh++){
    dlmran[uTh].SetSeed(uTh+1);
  }

  double kCats_min = 0;
  double kCats_max = 300;
  unsigned kCats_step = 25;
  CATS Kitty;
  Kitty.SetMomBins(kCats_step,kCats_min,kCats_max);
  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());
  AnalysisObject.SetUpCats_pL(Kitty,"Usmani","Gauss",0,0);
  Kitty.SetAnaSource(0,SourceSize*NuToFm);
  Kitty.SetChannelWeight(0,1);
  Kitty.SetChannelWeight(1,0);
  Kitty.KillTheCat();
  TH1F* hCkCats = new TH1F("hCkCats","hCkCats",kCats_step,kCats_min,kCats_max);
  for(unsigned uMom=0; uMom<kCats_step; uMom++){
    hCkCats->SetBinContent(uMom+1, Kitty.GetCorrFun(uMom));
  }

  double kstar_min = 0;
  double kstar_max = 4800;
  unsigned kstar_step = 400;

  TH1F* hSame = new TH1F("hSame","hSame",kstar_step,kstar_min,kstar_max);
  TH1F* hRef = new TH1F("hRef","hRef",kstar_step,kstar_min,kstar_max);

  #pragma omp parallel for
  for(unsigned uIter=0; uIter<NumIter; uIter++){
    unsigned ThId = omp_get_thread_num();

    
    double kstar0 = dlmran[ThId].GaussDiffR(3,0,kstar_width);
    double src = 0;
    do{
      src = dlmran[ThId].GaussDiffR(3,0,SourceSize);
    }
    while(src<source_cutoff);

    double pot_pars[3];
    pot_pars[0] = src;
    pot_pars[1] = 0;
    pot_pars[2] = 0;//spin
    double dist = pot_pars[0];
    double kstar = kstar0;
    double rho = dist*kstar;
    double drho = dr*kstar;
    for(double ftime=0; ftime<=tmax; ftime+=dt){
      pot_pars[0] = dist*NuToFm;
      //double pot = UsmaniPotentialCats(pot_pars);
      double pot = 90.*exp(-pow((rho)/0.5,2.)) - 20.*exp(-pow((rho)/1.,2.));
      pot_pars[0] = (dist+dr)*NuToFm;
      //double poth = UsmaniPotentialCats(pot_pars)
      double poth = 90.*exp(-pow((rho+drho)/0.5,2.)) - 20.*exp(-pow((rho+drho)/1.,2.));
      //in this convension, (+) force means the particles fly away from each other, i.e. repulsion
      double force = -(poth-pot)/dr;
      //in relativity, F = dp/dt, or for us dk = F*dt
      double dkstar = force*dt;
      kstar += dkstar;
      double beta = kstar/sqrt(kstar*kstar + red_mass*red_mass);
      dist += beta*dt;

      if(dist>=Vr_cutoff){
        break;
      }

      //printf("kstar = %.2f; rstar = %.2f\n", kstar, dist*NuToFm);
      //printf("   V = %.2f; dkstar = %.4f; dr = %.4f\n", pot, dkstar, beta*dt);
      //usleep(100e3);
    }
    #pragma omp critical
    {
    hSame->Fill(fabs(kstar));
    hRef->Fill(fabs(kstar0));
    }

  }

  TH1F* hCk = (TH1F*)hSame->Clone("hCk");
  hCk->Divide(hRef);

  TH1F* hCkNorm = (TH1F*)hCk->Clone("hCkNorm");
  double kstar_norm = hCkCats->GetBinCenter(kCats_step);
  double norm = hCkCats->GetBinContent(kCats_step)/hCkNorm->GetBinContent(hCkNorm->FindBin(kstar_norm));
  hCkNorm->Scale(norm);
  for(unsigned uMom=0; uMom<kstar_step; uMom++){
    hCkNorm->SetBinError(uMom+1, 0);
  }


  TFile fOutput(TString::Format("%s/ClassicalShit/Toy_MC_1/fOutput.root",GetCernBoxDimi()),"recreate");
  hSame->Write();
  hRef->Write();
  hCk->Write();
  hCkNorm->Write();
  hCkCats->Write();

  delete hSame;
  delete hRef;
  delete hCk;
  delete hCkNorm;
  delete hCkCats;
  delete [] dlmran;
}




int CLASSICALSHIT(int argc, char *argv[]){
  Toy_MC_1();
  return 0;
}
