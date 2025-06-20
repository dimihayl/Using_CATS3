
#include "FunWithCeca.h"
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
#include "DLM_MultiFit.h"
#include "CECA_Paper.h"

#include <iostream>
#include <unistd.h>
#include <omp.h>
#include <fstream>
#include <sstream>
#include <string>
#include <memory>

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
#include "TFitResultPtr.h"
#include "TFitResult.h"


#include <boost/algorithm/string.hpp>


void ceca_test(){
    const double Sigma = 0.6;
    const double LambdaFraction = 0.6;
    //you may put to zero for now
    const double HadronSize = 0.75;
    const double HadronSlope = 0.2;
    //const double HadronSize = 0;
    //const double HadronSlope = 0;
    TREPNI Database(0);
    Database.SetSeed(11);
    std::vector<TreParticle*> ParticleList;
    ParticleList.push_back(Database.NewParticle("Proton"));
    ParticleList.push_back(Database.NewParticle("Lambda"));
    ParticleList.push_back(Database.NewParticle("Pion"));
    ParticleList.push_back(Database.NewParticle("ProtonReso"));
    ParticleList.push_back(Database.NewParticle("LambdaReso"));
    double ExpPP,ExpPR,ExpRP,ExpRR;
    ExpPP = 0.3578*0.3562;
    ExpPR = 0.3578*0.6438;
    ExpRP = 0.6422*0.3562;
    ExpRR = 0.6422*0.6438;
    /*          Particle List          */
    for(TreParticle* prt : ParticleList){
        if(prt->GetName()=="Proton"){
            prt->SetMass(Mass_p);
            prt->SetAbundance(35.78);
            prt->SetRadius(HadronSize);
            prt->SetRadiusSlope(HadronSlope);
        }
        else if(prt->GetName()=="Lambda"){
            prt->SetMass(Mass_L);
            prt->SetAbundance(35.62*LambdaFraction);
            prt->SetRadius(HadronSize);
            prt->SetRadiusSlope(HadronSlope);
        }
        else if(prt->GetName()=="Pion"){
            prt->SetMass(Mass_pic);
            prt->SetAbundance(0);
            prt->SetRadius(HadronSize);
            prt->SetRadiusSlope(HadronSlope);
        }
        else if(prt->GetName()=="ProtonReso"){
            prt->SetMass(1362);
            prt->SetAbundance(64.22*1);
            prt->SetWidth(hbarc/1.65);
            prt->SetRadius(HadronSize);
            prt->SetRadiusSlope(HadronSlope);
            prt->NewDecay();
            prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Proton"));
            prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
            prt->GetDecay(0)->SetBranching(100);
        }
        else if(prt->GetName()=="LambdaReso"){
            prt->SetMass(1463);
            prt->SetAbundance(64.38*1);
            prt->SetWidth(hbarc/4.69);
            prt->SetRadius(HadronSize);
            prt->SetRadiusSlope(HadronSlope);
            prt->NewDecay();
            prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Lambda"));
            prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
            prt->GetDecay(0)->SetBranching(100);
        }
        //prt->SetMomPDF(Sigma*prt->GetMass());
    }
    printf("1\n");
    std::vector<std::string> ListOfParticles;
    ListOfParticles.push_back("Proton");
    //ListOfParticles.push_back("Proton");
    ListOfParticles.push_back("Lambda");
    CECA Ivana(Database,ListOfParticles);
    Ivana.SetDisplacement(0.0);
    Ivana.SetHadronization(1.0);
    Ivana.SetTau(0.0);
    Ivana.SetThermalKick(98);
    Ivana.SetTargetStatistics(10);
    Ivana.SetEventMult(2);
    Ivana.SetSourceDim(2);
    Ivana.SetDebugMode(true);
    Ivana.SetThreadTimeout(30);
    Ivana.GHETTO_EVENT = true;
    printf("2\n");
    Ivana.GoBabyGo(1);
    printf("3\n");
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
    printf("   Expected  Total  Femto\n");
    printf("PP %6.2f%% %6.2f%% %6.2f\n",ExpPP*100.,TotPP*100.,FemtoPP*100.);
    printf("PR %6.2f%% %6.2f%% %6.2f\n",ExpPR*100.,TotPR*100.,FemtoPR*100.);
    printf("RP %6.2f%% %6.2f%% %6.2f\n",ExpRP*100.,TotRP*100.,FemtoRP*100.);
    printf("RR %6.2f%% %6.2f%% %6.2f\n",ExpRR*100.,TotRR*100.,FemtoRR*100.);
    /*          Histograms          */
    TH1F* h_Ghetto_rstar            = Convert_DlmHisto_TH1F(Ivana.Ghetto_rstar,"Ghetto_rstar");
    TH1F* h_Ghetto_kstar            = Convert_DlmHisto_TH1F(Ivana.Ghetto_kstar,"Ghetto_kstar");
    TH2F* h_Ghetto_kstar_rstar      = Convert_DlmHisto_TH2F(Ivana.Ghetto_kstar_rstar,"Ghetto_kstar_rstar");
    TH2F* h_Ghetto_mT_rstar         = Convert_DlmHisto_TH2F(Ivana.Ghetto_mT_rstar,"Ghetto_mT_rstar");
    TH2F* h_GhettoFemto_mT_rstar    = Convert_DlmHisto_TH2F(Ivana.GhettoFemto_mT_rstar,"GhettoFemto_mT_rstar");
    TH2F* h_GhettoFemto_mT_kstar    = Convert_DlmHisto_TH2F(Ivana.GhettoFemto_mT_kstar,"GhettoFemto_mT_kstar");
    TH2F* h_Ghetto_mT_costh         = Convert_DlmHisto_TH2F(Ivana.Ghetto_mT_costh,"Ghetto_mT_costh");
    TH2F* h_GhettoSP_pT_th          = Convert_DlmHisto_TH2F(Ivana.GhettoSP_pT_th,"GhettoSP_pT_th");
    Ivana.GhettoFemto_rstar->ComputeError();
    TH1F* h_GhettoFemto_rstar       = Convert_DlmHisto_TH1F(Ivana.GhettoFemto_rstar,"GhettoFemto_rstar");
    TH1F* h_Ghetto_RP_AngleRcP1     = Convert_DlmHisto_TH1F(Ivana.Ghetto_RP_AngleRcP1,"Ghetto_RP_AngleRcP1");
    TH1F* h_Ghetto_PR_AngleRcP2     = Convert_DlmHisto_TH1F(Ivana.Ghetto_PR_AngleRcP2,"Ghetto_PR_AngleRcP2");
    TH1F* h_Ghetto_RR_AngleRcP1     = Convert_DlmHisto_TH1F(Ivana.Ghetto_RR_AngleRcP1,"Ghetto_RR_AngleRcP1");
    TH1F* h_Ghetto_RR_AngleRcP2     = Convert_DlmHisto_TH1F(Ivana.Ghetto_RR_AngleRcP2,"Ghetto_RR_AngleRcP2");
    TH1F* h_Ghetto_RR_AngleP1P2     = Convert_DlmHisto_TH1F(Ivana.Ghetto_RR_AngleP1P2,"Ghetto_RR_AngleP1P2");
    TGraphErrors g_GhettoFemto_mT_rstar;
    g_GhettoFemto_mT_rstar.SetName("g_GhettoFemto_mT_rstar");
    g_GhettoFemto_mT_rstar.SetMarkerStyle(20);
    g_GhettoFemto_mT_rstar.SetMarkerSize(1);
    g_GhettoFemto_mT_rstar.SetLineWidth(3);
    unsigned uPoint = 0;
    for(unsigned uBin=0; uBin<h_GhettoFemto_mT_rstar->GetXaxis()->GetNbins(); uBin++){
        TH1D* hProj = h_GhettoFemto_mT_rstar->ProjectionY(TString::Format("hProj"),uBin+1,uBin+1);
        double Mean = hProj->GetMean();
        double Err = hProj->GetStdDev();
        double mT = h_GhettoFemto_mT_rstar->GetXaxis()->GetBinCenter(uBin+1);
        if(Mean&&Err){
            g_GhettoFemto_mT_rstar.SetPoint(uPoint,mT,Mean);
            g_GhettoFemto_mT_rstar.SetPointError(uPoint,0,Err);
            uPoint++;
        }
    }
    const double MeanToGauss = 1./2.3;
    double mT_first = g_GhettoFemto_mT_rstar.GetPointY(0);
    printf("r(%.2f GeV) = %.2f\n",g_GhettoFemto_mT_rstar.GetPointX(0)*0.001,mT_first*MeanToGauss);
    double mT_2GeV = g_GhettoFemto_mT_rstar.Eval(2000);
    printf("r(%.2f GeV) = %.2f\n",2.,mT_2GeV*MeanToGauss);
    h_GhettoFemto_rstar->Scale(1./h_GhettoFemto_rstar->Integral(),"width");
    TF1* fSource = new TF1("fSource",GaussSourceTF1,0,5,1);
    fSource->SetParameter(0,1.2);
    fSource->SetParLimits(0,0.4,2.0);
    h_GhettoFemto_rstar->Fit(fSource,"Q, S, N, R, M");
    printf("The effective Gaussian size is %.3f +/- %.3f fm\n",fSource->GetParameter(0),fSource->GetParError(0));
    TFile fOutput("ceca_test.root", "recreate");
    h_Ghetto_rstar->Write();
    h_GhettoFemto_rstar->Write();
    fSource->Write();
    h_Ghetto_kstar->Write();
    h_Ghetto_kstar_rstar->Write();
    h_Ghetto_mT_rstar->Write();
    h_GhettoFemto_mT_rstar->Write();
    g_GhettoFemto_mT_rstar.Write();
    h_GhettoFemto_mT_kstar->Write();
    h_Ghetto_mT_costh->Write();
    h_GhettoSP_pT_th->Write();
    h_Ghetto_RP_AngleRcP1->Write();
    h_Ghetto_PR_AngleRcP2->Write();
    h_Ghetto_RR_AngleRcP1->Write();
    h_Ghetto_RR_AngleRcP2->Write();
    h_Ghetto_RR_AngleP1P2->Write();
}


void Test1(){
  unsigned char ch_1 = 255;
  std::cout << unsigned(ch_1) << std::endl;
  std::cout << (ch_1==255) << std::endl;
  std::cout << (ch_1==-1) << std::endl;

  /*
  TREPNI Trepni(0);
  //Trepni.SetPrintLevel(2,false);

  std::string str ;
  //str = Trepni.GetParticleName(0);
  //str = Trepni.GetParticleName(1);

  Trepni.SetParticle("Test", 100, 100, 300, 300);

  Trepni.GetParticleName(0);
  str = Trepni.GetParticleName(0);
  str = Trepni.GetParticleName(0);

  str = Trepni.GetParticleName(1);
  printf("PN1 = %s\n",str.c_str());

  str = Trepni.GetParticleName(-1);
  printf("PNm1 = %s\n",str.c_str());

  printf("ID of Test: %i\n",Trepni.GetParticleId("Test"));
  printf("ID of anti_Test: %i\n",Trepni.GetParticleId("anti_Test"));
  printf("ID of John: %i\n",Trepni.GetParticleId("John"));

  Trepni.SetParticle("anti_Silly", 200, 200, 320, 320);
  Trepni.SetParticle("", 200, 200, 320, 320);
  Trepni.SetParticle("Silly", 200, 200, 320, 320);
  Trepni.SetParticle("SillyNess_-", 200, 200, 320, 320);
  Trepni.SetParticle("test1", -200, 200, 320, 320);
  Trepni.SetParticle("test1", 200, 200, 350, 320);
  */
}

void Test2(){

  const double Lambda_width = (0.0789)*pow(10,-15)*FmToNu;

  TREPNI Database(0);
  TreParticle* p_pim = Database.NewParticle("pi-");
  TreParticle* p_pi0 = Database.NewParticle();
  TreParticle* p_proton = Database.NewParticle();
  TreParticle* p_neutron = Database.NewParticle();
  TreParticle* p_Lambda = Database.NewParticle();

  //TreParticle p_pim(Database);
  //TreParticle p_pi0(Database);
  //TreParticle p_proton(Database);
  //TreParticle p_neutron(Database);
  //TreParticle p_Lambda(Database);


  //p_pim->SetName("pi-");
  Database.GetParticle("pi-")->SetMass(1.*Mass_pic);
  p_pim->SetMassLimit(1.*Mass_pic,1.*Mass_pic);
  p_pim->SetWidth(0);
  p_pim->SetWidthLimit(0,0);
  p_pim->SetAbundance(0);
  p_pim->SetAbundanceLimit(0,0);

  p_pi0->SetName("pi0");
  p_pi0->SetMass(1.*Mass_pi0);
  p_pi0->SetMassLimit(1.*Mass_pi0,1.*Mass_pi0);
  p_pi0->SetWidth(0);
  p_pi0->SetWidthLimit(0,0);
  p_pi0->SetAbundance(0);
  p_pi0->SetAbundanceLimit(0,0);

  p_proton->SetName("proton");
  p_proton->SetMass(Mass_p);
  p_proton->SetMassLimit(Mass_p,Mass_p);
  p_proton->SetWidth(0);
  p_proton->SetWidthLimit(0,0);
  p_proton->SetAbundance(0);
  p_proton->SetAbundanceLimit(0,0);

  p_neutron->SetName("neutron");
  p_neutron->SetMass(Mass_n);
  p_neutron->SetMassLimit(Mass_n,Mass_n);
  p_neutron->SetWidth(0);
  p_neutron->SetWidthLimit(0,0);
  p_neutron->SetAbundance(0);
  p_neutron->SetAbundanceLimit(0,0);

  p_Lambda->SetName("Lambda");
  p_Lambda->SetMass(Mass_L);
  p_Lambda->SetMassLimit(Mass_L,Mass_L);
  p_Lambda->SetWidth(Lambda_width);
  p_Lambda->SetWidthLimit(Lambda_width,Lambda_width);
  p_Lambda->SetAbundance(100);
  p_Lambda->SetAbundanceLimit(50,150);

  p_Lambda->NewDecay();
  p_Lambda->GetDecay(0)->AddDaughter(*p_pim);
  p_Lambda->GetDecay(0)->AddDaughter(*p_proton);
  p_Lambda->GetDecay(0)->SetBranching(67);
  //p_Lambda->GetDecay(0)->SetBranchingLimit(1, 1);

  p_Lambda->NewDecay();
  p_Lambda->GetDecay(1)->AddDaughter(*p_pi0);
  p_Lambda->GetDecay(1)->AddDaughter(*p_neutron);
  p_Lambda->GetDecay(1)->SetBranching(33);
  p_Lambda->GetDecay(1)->SetBranchingLimit(33, 33);
  //
  p_pim->Print();
  p_proton->Print();
  p_Lambda->Print();

  Database.SetTotalYield(100);
  bool database_qa = Database.QA();
  printf("Database QA: %i\n",database_qa);

}

//test: we assume we have
void TestingFractions(){
  const int NumEvents = 1000000;
  const int EventMult = 3;
  std::vector<float> ResoFrac;
  ResoFrac.push_back(0.24);//0   decay to A
  ResoFrac.push_back(0.12);//1   decay to B
  ResoFrac.push_back(0.11);//2  decay to C
  ResoFrac.push_back(0.1);//3    decay to AB
  ResoFrac.push_back(0.14);//4   decay to AC
  ResoFrac.push_back(0.2);//5    decay to BC
  ResoFrac.push_back(0.05);//6    decay to ABC
  ResoFrac.push_back(0.04);//7   decay into CC

  std::vector<float> ResoFracCum;
  for(unsigned uR=0; uR<ResoFrac.size(); uR++){
    if(uR) ResoFracCum.push_back(ResoFracCum.at(uR-1)+ResoFrac.at(uR));
    else ResoFracCum.push_back(ResoFrac.at(uR));
  }


  TRandom3 RanGen(11);
  TH1F* hFracReso = new TH1F("hFracReso","hFracReso",ResoFrac.size(),-0.5,-0.5+ResoFrac.size());
  TH1F* hFracPart = new TH1F("hFracPart","hFracPart",EventMult,-0.5,-0.5+EventMult);
  TH1F* hTrueMult = new TH1F("hTrueMult","hTrueMult",11,-0.5,10.5);
  TH1F* hFracPart3 = new TH1F("hFracPart3","hFracPart3",EventMult,-0.5,-0.5+EventMult);
  TH1F* hFracReso3 = new TH1F("hFracReso3","hFracReso3",ResoFrac.size(),-0.5,-0.5+ResoFrac.size());
  TH1F* hFracPart3f = new TH1F("hFracPart3f","hFracPart3f",EventMult,-0.5,-0.5+EventMult);
  TH1F* hFracReso3f = new TH1F("hFracReso3f","hFracReso3f",ResoFrac.size(),-0.5,-0.5+ResoFrac.size());
  TH1F** hFracResoT_P = new TH1F* [3];
  TH1F** hFracReso_P = new TH1F* [3];
  TH1F** hFracReso3_P = new TH1F* [3];
  TH1F** hFracReso3f_P = new TH1F* [3];
  for(unsigned uP=0; uP<3; uP++){
    hFracReso_P[uP] = new TH1F(TString::Format("hFracReso_P%u",uP),TString::Format("hFracReso_P%u",uP)
    ,ResoFrac.size(),-0.5,-0.5+ResoFrac.size());
    hFracReso3_P[uP] = new TH1F(TString::Format("hFracReso3_P%u",uP),TString::Format("hFracReso3_P%u",uP)
    ,ResoFrac.size(),-0.5,-0.5+ResoFrac.size());
    hFracReso3f_P[uP] = new TH1F(TString::Format("hFracReso3f_P%u",uP),TString::Format("hFracReso3f_P%u",uP)
    ,ResoFrac.size(),-0.5,-0.5+ResoFrac.size());
    hFracResoT_P[uP] = new TH1F(TString::Format("hFracResoT_P%u",uP),TString::Format("hFracResoT_P%u",uP)
    ,ResoFrac.size(),-0.5,-0.5+ResoFrac.size());
  }

  //expected fractions of resonances feeding into particular species (0,1,2)
  double** TheoMammaFrac = new double* [3];
  double* TheoMammaNorm = new double [3];
  for(unsigned uP=0; uP<3; uP++){
    TheoMammaNorm[uP] = 0;
    TheoMammaFrac[uP] = new double [ResoFrac.size()];
    for(unsigned uR=0; uR<ResoFrac.size(); uR++) TheoMammaFrac[uP][uR]=0;
    if(uP==0){
      TheoMammaFrac[uP][0] = ResoFrac.at(0)*1.;
      TheoMammaFrac[uP][3] = ResoFrac.at(3)*1.;
      TheoMammaFrac[uP][4] = ResoFrac.at(4)*1.;
      TheoMammaFrac[uP][6] = ResoFrac.at(6)*1.;
    }
    if(uP==1){
      TheoMammaFrac[uP][1] = ResoFrac.at(1)*1.;
      TheoMammaFrac[uP][3] = ResoFrac.at(3)*1.;
      TheoMammaFrac[uP][5] = ResoFrac.at(5)*1.;
      TheoMammaFrac[uP][6] = ResoFrac.at(6)*1.;
    }
    if(uP==2){
      TheoMammaFrac[uP][2] = ResoFrac.at(2)*1.;
      TheoMammaFrac[uP][4] = ResoFrac.at(4)*1.;
      TheoMammaFrac[uP][5] = ResoFrac.at(5)*1.;
      TheoMammaFrac[uP][6] = ResoFrac.at(6)*1.;
      TheoMammaFrac[uP][7] = ResoFrac.at(7)*2.;
    }
    for(unsigned uR=0; uR<ResoFrac.size(); uR++)
      TheoMammaNorm[uP] += TheoMammaFrac[uP][uR];
    for(unsigned uR=0; uR<ResoFrac.size(); uR++)
      TheoMammaFrac[uP][uR] /= TheoMammaNorm[uP];
  }

  for(unsigned uP=0; uP<3; uP++){
    for(unsigned uR=0; uR<ResoFrac.size(); uR++){
      hFracResoT_P[uP]->SetBinContent(uR+1,TheoMammaFrac[uP][uR]);
      hFracResoT_P[uP]->SetBinError(uR+1,0);
    }
    hFracResoT_P[uP]->SetMarkerStyle(1);
    hFracResoT_P[uP]->SetMarkerSize(1.5);
  }
  for(unsigned uP=0; uP<3; uP++) delete [] TheoMammaFrac[uP];
  delete [] TheoMammaFrac;
  delete [] TheoMammaNorm;

  for(unsigned uEvent=0; uEvent<NumEvents; uEvent++){
    std::vector<int> ResonanceID;
    std::vector<int> ParticleID;
    std::vector<int> MotherID;
    for(unsigned uMult=0; uMult<EventMult; uMult++){
      double sample = RanGen.Uniform();
      unsigned WhichReso = ResoFrac.size()-1;
      for(unsigned uR=0; uR<ResoFrac.size()-1; uR++){
        if(sample<ResoFracCum.at(uR)){
          WhichReso = uR;
          break;
        }
      }
      //the decays can be into particles 0 1 2, and the decay is
      //either into 1,2 or all 3 of the studied species
      ResonanceID.push_back(WhichReso);
      switch (WhichReso) {
        case 0: ParticleID.push_back(0);
                MotherID.push_back(0);
                break;
        case 1: ParticleID.push_back(1);
                MotherID.push_back(1);
                break;
        case 2: ParticleID.push_back(2);
                MotherID.push_back(2);
                break;
        case 3: ParticleID.push_back(0);
                ParticleID.push_back(1);
                MotherID.push_back(3);
                MotherID.push_back(3);
                break;
        case 4: ParticleID.push_back(0);
                ParticleID.push_back(2);
                MotherID.push_back(4);
                MotherID.push_back(4);
                break;
        case 5: ParticleID.push_back(1);
                ParticleID.push_back(2);
                MotherID.push_back(5);
                MotherID.push_back(5);
                break;
        case 6: ParticleID.push_back(0);
                ParticleID.push_back(1);
                ParticleID.push_back(2);
                MotherID.push_back(6);
                MotherID.push_back(6);
                MotherID.push_back(6);
                break;
        case 7: ParticleID.push_back(2);
                ParticleID.push_back(2);
                MotherID.push_back(7);
                MotherID.push_back(7);
                break;
        default : printf("!error!\n"); break;
      }
      hFracReso->Fill(WhichReso);
    }//uMult


    hTrueMult->Fill(ParticleID.size());
    for(unsigned uPart=0; uPart<ParticleID.size(); uPart++){
      hFracPart->Fill(ParticleID.at(uPart));
      for(unsigned uP=0; uP<3; uP++)
        if(ParticleID.at(uPart)==uP) hFracReso_P[uP]->Fill(MotherID.at(uPart));
    }

    for(unsigned uP0=0; uP0<ParticleID.size(); uP0++){
      for(unsigned uP1=uP0; uP1<ParticleID.size(); uP1++){
        for(unsigned uP2=uP1; uP2<ParticleID.size(); uP2++){
          bool P0 = (ParticleID.at(uP0)==0||ParticleID.at(uP1)==0||ParticleID.at(uP2)==0);
          bool P1 = (ParticleID.at(uP0)==1||ParticleID.at(uP1)==1||ParticleID.at(uP2)==1);
          bool P2 = (ParticleID.at(uP0)==2||ParticleID.at(uP1)==2||ParticleID.at(uP2)==2);
          //all species of interest
          bool GoodTriplet = (P0&&P1&&P2);
          //all specieas of interest + their mothers are different
          //the latter is probably somewhat true for low k* / Q3 pairs (is it true for Q3 ??)
          bool FemtoTriplet = GoodTriplet;
          FemtoTriplet *= (MotherID.at(uP0)!=MotherID.at(uP1));
          FemtoTriplet *= (MotherID.at(uP0)!=MotherID.at(uP2));
          FemtoTriplet *= (MotherID.at(uP1)!=MotherID.at(uP2));
          if(GoodTriplet){
            hFracPart3->Fill(ParticleID.at(uP0));
            hFracPart3->Fill(ParticleID.at(uP1));
            hFracPart3->Fill(ParticleID.at(uP2));
            hFracReso3->Fill(MotherID.at(uP0));
            hFracReso3->Fill(MotherID.at(uP1));
            hFracReso3->Fill(MotherID.at(uP2));
            for(unsigned uP=0; uP<3; uP++){
              if(ParticleID.at(uP)==0) hFracReso3_P[0]->Fill(MotherID.at(uP));
              if(ParticleID.at(uP)==1) hFracReso3_P[1]->Fill(MotherID.at(uP));
              if(ParticleID.at(uP)==2) hFracReso3_P[2]->Fill(MotherID.at(uP));
            }
          }
          if(FemtoTriplet){
            hFracPart3f->Fill(ParticleID.at(uP0));
            hFracPart3f->Fill(ParticleID.at(uP1));
            hFracPart3f->Fill(ParticleID.at(uP2));
            hFracReso3f->Fill(MotherID.at(uP0));
            hFracReso3f->Fill(MotherID.at(uP1));
            hFracReso3f->Fill(MotherID.at(uP2));
            for(unsigned uP=0; uP<3; uP++){
              if(ParticleID.at(uP)==0) hFracReso3f_P[0]->Fill(MotherID.at(uP));
              if(ParticleID.at(uP)==1) hFracReso3f_P[1]->Fill(MotherID.at(uP));
              if(ParticleID.at(uP)==2) hFracReso3f_P[2]->Fill(MotherID.at(uP));
            }
          }
        }
      }
    }
  }//uEvent
  hFracReso->Scale(1./hFracReso->Integral(),"width");
  hFracPart->Scale(1./hFracPart->Integral(),"width");
  hTrueMult->Scale(1./hTrueMult->Integral(),"width");
  hFracPart3->Scale(1./hFracPart3->Integral(),"width");
  hFracReso3->Scale(1./hFracReso3->Integral(),"width");
  hFracPart3f->Scale(1./hFracPart3f->Integral(),"width");
  hFracReso3f->Scale(1./hFracReso3f->Integral(),"width");
  for(unsigned uP=0; uP<3; uP++){
    hFracResoT_P[uP]->Scale(1./hFracResoT_P[uP]->Integral(),"width");
    hFracReso_P[uP]->Scale(1./hFracReso_P[uP]->Integral(),"width");
    hFracReso3_P[uP]->Scale(1./hFracReso3_P[uP]->Integral(),"width");
    hFracReso3f_P[uP]->Scale(1./hFracReso3f_P[uP]->Integral(),"width");
  }

  TFile fOutput(TString::Format("%s/FunWithCeca/TestingFractions.root",GetFemtoOutputFolder()),"recreate");
  for(unsigned uP=0; uP<3; uP++) hFracResoT_P[uP]->Write();
  for(unsigned uP=0; uP<3; uP++) hFracReso_P[uP]->Write();
  hFracReso->Write();
  for(unsigned uP=0; uP<3; uP++) hFracReso3_P[uP]->Write();
  hFracReso3->Write();
  for(unsigned uP=0; uP<3; uP++) hFracReso3f_P[uP]->Write();
  hFracReso3f->Write();
  hFracPart->Write();
  hFracPart3->Write();
  hFracPart3f->Write();
  hTrueMult->Write();


  delete hFracReso;
  delete hFracPart;
  delete hTrueMult;
  delete hFracPart3;
  delete hFracReso3;
  delete hFracPart3f;
  delete hFracReso3f;
  for(unsigned uP=0; uP<3; uP++){
    delete hFracReso_P[uP];
    delete hFracReso3_P[uP];
    delete hFracReso3f_P[uP];
  }
  delete [] hFracReso_P;
  delete [] hFracReso3_P;
  delete [] hFracReso3f_P;
}


void TestTime(){
  DLM_Timer Timer;
  Timer.Start();
  usleep(1000e3);
  long long value = Timer.Stop();
  double dval = value;
  printf("dval = %f\n",dval);

  printf("Nth = %u\n",GetNumThreads());
  printf("Oth = %u\n",GetOptimalThreads());
}

void TestDecay(){
  const unsigned NumIter = 100000;
  const double MotherMass = 1116;
  const double MotherMom = 1000;
  const double DaughterMass1 = 938;
  const double DaughterMass2 = 140;

  TH1F* hGenBod_Mom1_CM = new TH1F("hGenBod_Mom1_CM","hGenBod_Mom1_CM",128,0,1024);
  TH1F* hGenBod_CosTheta1_CM = new TH1F("hGenBod_CosTheta1_CM","hGenBod_CosTheta1_CM",128,-1,1);
  TH1F* hGenBod_Phi1_CM = new TH1F("hGenBod_Phi1_CM","hGenBod_Phi1_CM",128,-2.*TMath::Pi(),2.*TMath::Pi());

  TH1F* hGenBod_Mom2_CM = new TH1F("hGenBod_Mom2_CM","hGenBod_Mom2_CM",128,0,1024);
  TH1F* hGenBod_CosTheta2_CM = new TH1F("hGenBod_CosTheta2_CM","hGenBod_CosTheta2_CM",128,-1,1);
  TH1F* hGenBod_Phi2_CM = new TH1F("hGenBod_Phi2_CM","hGenBod_Phi2_CM",128,-2.*TMath::Pi(),2.*TMath::Pi());

  TH1F* hGenBod_Mom1_LAB = new TH1F("hGenBod_Mom1_LAB","hGenBod_Mom1_LAB",128,0,1024);
  TH1F* hGenBod_CosTheta1_LAB = new TH1F("hGenBod_CosTheta1_LAB","hGenBod_CosTheta1_LAB",128,-1,1);
  TH1F* hGenBod_Phi1_LAB = new TH1F("hGenBod_Phi1_LAB","hGenBod_Phi1_LAB",128,-2.*TMath::Pi(),2.*TMath::Pi());

  TH1F* hGenBod_Mom2_LAB = new TH1F("hGenBod_Mom2_LAB","hGenBod_Mom2_LAB",128,0,1024);
  TH1F* hGenBod_CosTheta2_LAB = new TH1F("hGenBod_CosTheta2_LAB","hGenBod_CosTheta2_LAB",128,-1,1);
  TH1F* hGenBod_Phi2_LAB = new TH1F("hGenBod_Phi2_LAB","hGenBod_Phi2_LAB",128,-2.*TMath::Pi(),2.*TMath::Pi());


  TH1F* hDimi_Mom1_CM = new TH1F("hDimi_Mom1_CM","hDimi_Mom1_CM",128,0,1024);
  TH1F* hDimi_CosTheta1_CM = new TH1F("hDimi_CosTheta1_CM","hDimi_CosTheta1_CM",128,-1,1);
  TH1F* hDimi_Phi1_CM = new TH1F("hDimi_Phi1_CM","hDimi_Phi1_CM",128,-2.*TMath::Pi(),2.*TMath::Pi());

  TH1F* hDimi_Mom2_CM = new TH1F("hDimi_Mom2_CM","hDimi_Mom2_CM",128,0,1024);
  TH1F* hDimi_CosTheta2_CM = new TH1F("hDimi_CosTheta2_CM","hDimi_CosTheta2_CM",128,-1,1);
  TH1F* hDimi_Phi2_CM = new TH1F("hDimi_Phi2_CM","hDimi_Phi2_CM",128,-2.*TMath::Pi(),2.*TMath::Pi());

  TH1F* hDimi_Mom1_LAB = new TH1F("hDimi_Mom1_LAB","hDimi_Mom1_LAB",128,0,1024);
  TH1F* hDimi_CosTheta1_LAB = new TH1F("hDimi_CosTheta1_LAB","hDimi_CosTheta1_LAB",128,-1,1);
  TH1F* hDimi_Phi1_LAB = new TH1F("hDimi_Phi1_LAB","hDimi_Phi1_LAB",128,-2.*TMath::Pi(),2.*TMath::Pi());

  TH1F* hDimi_Mom2_LAB = new TH1F("hDimi_Mom2_LAB","hDimi_Mom2_LAB",128,0,1024);
  TH1F* hDimi_CosTheta2_LAB = new TH1F("hDimi_CosTheta2_LAB","hDimi_CosTheta2_LAB",128,-1,1);
  TH1F* hDimi_Phi2_LAB = new TH1F("hDimi_Phi2_LAB","hDimi_Phi2_LAB",128,-2.*TMath::Pi(),2.*TMath::Pi());

  double MassArray[2];
  MassArray[0] = DaughterMass1;
  MassArray[1] = DaughterMass2;

  DLM_Random DimiGen(11);

  CatsParticle MotherAtRest;
  MotherAtRest.SetMomXYZ(0,0,0);
  MotherAtRest.SetMass(MotherMass);
  MotherAtRest.SetDecayRanGen(DimiGen);

  CatsParticle MotherMoving;
  MotherMoving.SetMomXYZ(0,0,MotherMom);
  MotherMoving.SetMass(MotherMass);
  MotherMoving.SetDecayRanGen(DimiGen);

  TGenPhaseSpace psgenRest;
  TLorentzVector MammaRest;
  MammaRest.SetXYZM(0,0,0,MotherMass);
  psgenRest.SetDecay(MammaRest,2,MassArray);

  TGenPhaseSpace psgenMoving;
  TLorentzVector MammaMoving;
  MammaMoving.SetXYZM(0,0,MotherMom,MotherMass);
  psgenMoving.SetDecay(MammaMoving,2,MassArray);

  for(unsigned uIter=0; uIter<NumIter; uIter++){
//printf("u = %u\n", uIter);
    CatsParticle* DaughtersCM = MotherAtRest.Decay(DaughterMass1,DaughterMass2);
    CatsParticle* DaughtersLAB = MotherMoving.Decay(DaughterMass1,DaughterMass2);

    psgenRest.Generate();
    psgenMoving.Generate();

    TLorentzVector* GenDaughter1_CM = psgenRest.GetDecay(0);
    TLorentzVector* GenDaughter2_CM = psgenRest.GetDecay(1);

    TLorentzVector* GenDaughter1_LAB = psgenMoving.GetDecay(0);
    TLorentzVector* GenDaughter2_LAB = psgenMoving.GetDecay(1);

    hGenBod_Mom1_CM->Fill(GenDaughter1_CM->P());
    hGenBod_CosTheta1_CM->Fill(GenDaughter1_CM->CosTheta());
    hGenBod_Phi1_CM->Fill(GenDaughter1_CM->Phi());

    hGenBod_Mom2_CM->Fill(GenDaughter2_CM->P());
    hGenBod_CosTheta2_CM->Fill(GenDaughter2_CM->CosTheta());
    hGenBod_Phi2_CM->Fill(GenDaughter2_CM->Phi());

    hGenBod_Mom1_LAB->Fill(GenDaughter1_LAB->P());
    hGenBod_CosTheta1_LAB->Fill(GenDaughter1_LAB->CosTheta());
    hGenBod_Phi1_LAB->Fill(GenDaughter1_LAB->Phi());

    hGenBod_Mom2_LAB->Fill(GenDaughter2_LAB->P());
    hGenBod_CosTheta2_LAB->Fill(GenDaughter2_LAB->CosTheta());
    hGenBod_Phi2_LAB->Fill(GenDaughter2_LAB->Phi());

//printf("  p1 = %f\n", DaughtersCM[0].GetP());
    hDimi_Mom1_CM->Fill(DaughtersCM[0].GetP());
    hDimi_CosTheta1_CM->Fill(cos(DaughtersCM[0].GetPtheta()));
    hDimi_Phi1_CM->Fill(DaughtersCM[0].GetPphi());

    hDimi_Mom2_CM->Fill(DaughtersCM[1].GetP());
    hDimi_CosTheta2_CM->Fill(cos(DaughtersCM[1].GetPtheta()));
    hDimi_Phi2_CM->Fill(DaughtersCM[1].GetPphi());

    hDimi_Mom1_LAB->Fill(DaughtersLAB[0].GetP());
    hDimi_CosTheta1_LAB->Fill(cos(DaughtersLAB[0].GetPtheta()));
    hDimi_Phi1_LAB->Fill(DaughtersLAB[0].GetPphi());

    hDimi_Mom2_LAB->Fill(DaughtersLAB[1].GetP());
    hDimi_CosTheta2_LAB->Fill(cos(DaughtersLAB[1].GetPtheta()));
    hDimi_Phi2_LAB->Fill(DaughtersLAB[1].GetPphi());

    delete [] DaughtersCM;
    delete [] DaughtersLAB;

    //delete GenDaughter1_CM;
    //delete GenDaughter2_CM;
    //delete GenDaughter1_LAB;
    //delete GenDaughter2_LAB;

  }

  TFile fOutput(TString::Format("%s/FunWithCeca/TestDecay.root",GetFemtoOutputFolder()),"recreate");

  hGenBod_Mom1_CM->Write();
  hGenBod_CosTheta1_CM->Write();
  hGenBod_Phi1_CM->Write();
  hGenBod_Mom2_CM->Write();
  hGenBod_CosTheta2_CM->Write();
  hGenBod_Phi2_CM->Write();

  hGenBod_Mom1_LAB->Write();
  hGenBod_CosTheta1_LAB->Write();
  hGenBod_Phi1_LAB->Write();
  hGenBod_Mom2_LAB->Write();
  hGenBod_CosTheta2_LAB->Write();
  hGenBod_Phi2_LAB->Write();

  hDimi_Mom1_CM->Write();
  hDimi_CosTheta1_CM->Write();
  hDimi_Phi1_CM->Write();
  hDimi_Mom2_CM->Write();
  hDimi_CosTheta2_CM->Write();
  hDimi_Phi2_CM->Write();

  hDimi_Mom1_LAB->Write();
  hDimi_CosTheta1_LAB->Write();
  hDimi_Phi1_LAB->Write();
  hDimi_Mom2_LAB->Write();
  hDimi_CosTheta2_LAB->Write();
  hDimi_Phi2_LAB->Write();

}

//conclusion
//the time coordinate changes the spacial ones after the boost,
//however ALL relative coordinates (space and time) remain unchanged
//i.e. it does not matter what time we choose in our LAB, as long as it is
//consistent for all particles.
//N.B. same time in LAB is NOT same time in CM, should we propagate ??
void TestBoost(){

  CatsLorentzVector BoostV;
  BoostV.Set(0,0,0,0,500,200,0,0);

  CatsLorentzVector SV1;
  SV1.Set(0,0.01,-0.02,0.01,0,0,0,0);

  CatsLorentzVector SV1A;
  SV1A.Set(0,0.02,-0.02,0.01,0,0,0,0);

  CatsLorentzVector SV2;
  SV2.Set(sqrt(0.01*0.01+0.02*0.02+0.01*0.01),0.01,-0.02,0.01,0,0,0,0);

  CatsLorentzVector SV2A;
  SV2A.Set(sqrt(0.01*0.01+0.02*0.02+0.01*0.01),0.02,-0.02,0.01,0,0,0,0);

  CatsLorentzVector SV3;
  SV3.Set(sqrt(0.02*0.02+0.01*0.01+0.02*0.02+0.01*0.01),0.01,-0.02,0.01,0,0,0,0);

  printf("\n");
  SV1.Print();
  printf("\n");
  SV2.Print();
  printf("\n");
  SV3.Print();

  SV1.Boost(BoostV);
  SV1A.Boost(BoostV);
  SV2.Boost(BoostV);
  SV2A.Boost(BoostV);
  SV3.Boost(BoostV);


  printf("\n");
  printf("Gamma = %.3f\n",BoostV.Gamma());
  printf("Beta  = %.3f\n",BoostV.Beta());
  printf("BetaX = %.3f\n",BoostV.BetaX());

  printf("\n");
  SV1.Print();
  printf("\n");
  SV2.Print();
  printf("\n");
  SV3.Print();

  printf("Relative distance (x) in the final system:\n");
  printf(" SV1 - SV1A: %.2e %.2e\n",SV1.GetX()-SV1A.GetX(),SV1.GetT()-SV1A.GetT());
  printf(" SV2 - SV2A: %.2e %.2e\n",SV2.GetX()-SV2A.GetX(),SV2.GetT()-SV2A.GetT());

}

void Ghetto_Test1(){
  const unsigned NumPairs = 32*1000*1000;
  //const unsigned NumPairs = 1000;
  const double r_SP = 1.0;//1.35
  const double r_TR = 1.0;//1.1
  const double p_SP = 600;
  const double FemtoLimit = 255;

  TREPNI Database(0);
  std::vector<std::string> ListOfParticles;
  ListOfParticles.push_back("proton");
  ListOfParticles.push_back("proton");
  ListOfParticles.push_back("Lambda");
  CECA GT1(Database, ListOfParticles);
  GT1.SetFemtoRegion(FemtoLimit);
  GT1.GhettoTest1(NumPairs,r_SP,p_SP);

  TFile fOutput(TString::Format("%s/FunWithCeca/Ghetto_Test1_1.root",GetFemtoOutputFolder()),"recreate");

  TH2F* h_Ghetto_rstar = Convert_DlmHisto_TH2F(GT1.Ghetto_rstar,"Ghetto_rstar");
  TH2F* h_Ghetto_rcore = Convert_DlmHisto_TH2F(GT1.Ghetto_rcore,"Ghetto_rcore");

  TH1F* hFemto_Ghetto_rstar = (TH1F*)h_Ghetto_rstar->ProjectionY(TString::Format("hFemto_Ghetto_rstar"),1,h_Ghetto_rstar->GetXaxis()->FindBin(FemtoLimit));
  //hFemto_Ghetto_rstar->Sumw2();
  hFemto_Ghetto_rstar->Scale(1./hFemto_Ghetto_rstar->Integral(),"width");

  TH1F* hFemto_Ghetto_rcore = (TH1F*)h_Ghetto_rcore->ProjectionY(TString::Format("hFemto_Ghetto_rcore"),1,h_Ghetto_rcore->GetXaxis()->FindBin(FemtoLimit));
  //hFemto_Ghetto_rcore->Sumw2();
  hFemto_Ghetto_rcore->Scale(1./hFemto_Ghetto_rcore->Integral(),"width");

  TH1D* h_kstar = h_Ghetto_rcore->ProjectionX(TString::Format("h_kstar"),1,h_Ghetto_rcore->GetYaxis()->FindBin(h_Ghetto_rcore->GetYaxis()->GetNbins()));
  h_kstar->Scale(1./h_kstar->Integral(),"width");

  double rMin,rMax;

  GetCentralInterval(*hFemto_Ghetto_rstar, 0.95, rMin, rMax, true);
  TF1* fSource_star = new TF1("fSource_star",GaussSourceTF1,rMin,rMax,1);
  fSource_star->SetParameter(0,r_SP);
  fSource_star->SetParLimits(0,r_SP*0.25,r_SP*4.0);
  hFemto_Ghetto_rstar->Fit(fSource_star,"Q, S, N, R, M");

  GetCentralInterval(*hFemto_Ghetto_rcore, 0.95, rMin, rMax, true);
  TF1* fSource_core = new TF1("fSource_core",GaussSourceTF1,rMin,rMax,1);
  fSource_core->SetParameter(0,r_SP);
  fSource_core->SetParLimits(0,r_SP*0.25,r_SP*4.0);
  hFemto_Ghetto_rcore->Fit(fSource_core,"Q, S, N, R, M");

  printf("r_star = %.3f +/- %.3f\n",fSource_star->GetParameter(0),fSource_star->GetParError(0));
  printf("r_core = %.3f +/- %.3f\n",fSource_core->GetParameter(0),fSource_core->GetParError(0));

  TH1F* hFemto_Old_rstar = Convert_DlmHisto_TH1F(GT1.Old_rstar,"Old_rstar");
  TH1F* hFemto_Old_rcore = Convert_DlmHisto_TH1F(GT1.Old_rcore,"Old_rcore");

  TH1F* hFemto_Old_CosRcP1 = Convert_DlmHisto_TH1F(GT1.Old_CosRcP1,"Old_CosRcP1");
  TH1F* hFemto_Old_CosRcP2 = Convert_DlmHisto_TH1F(GT1.Old_CosRcP2,"Old_CosRcP2");
  TH1F* hFemto_Old_CosP1P2 = Convert_DlmHisto_TH1F(GT1.Old_CosP1P2,"Old_CosP1P2");

  TH1F* hFemto_Old_RcP1 = Convert_DlmHisto_TH1F(GT1.Old_RcP1,"Old_RcP1");
  TH1F* hFemto_Old_RcP2 = Convert_DlmHisto_TH1F(GT1.Old_RcP2,"Old_RcP2");
  TH1F* hFemto_Old_P1P2 = Convert_DlmHisto_TH1F(GT1.Old_P1P2,"Old_P1P2");

  TH1F* hGhettOld_rstar = Convert_DlmHisto_TH1F(GT1.GhettOld_rstar,"GhettOld_rstar");
  hGhettOld_rstar->Scale(1./hGhettOld_rstar->Integral(),"width");

  TH1F* hFemto_OldNew_rstar = (TH1F*)hFemto_Old_rstar->Clone("OldNew_rstar");
  TH1F* hFemto_OldNew_rcore = (TH1F*)hFemto_Old_rcore->Clone("OldNew_rcore");
  for(unsigned uBin=0; uBin<hFemto_Old_rstar->GetNbinsX(); uBin++){
    double rad = hFemto_Old_rstar->GetBinCenter(uBin+1);
    double pars[2];
    pars[0] = fSource_core->GetParameter(0);
    pars[1] = 2.;
    hFemto_OldNew_rstar->SetBinContent(uBin+1,GT1.Old_source->RootEval(&rad,pars));
    hFemto_OldNew_rcore->SetBinContent(uBin+1,GaussSourceTF1(&rad,pars));
  }

  DLM_Histo<float>* dlmFemto_Ghetto_rstar = Convert_TH1F_DlmHisto(hFemto_Ghetto_rstar);
  dlmFemto_Ghetto_rstar->RebinAll(4);
  *dlmFemto_Ghetto_rstar/=4.;
  //printf("%u\n",dlmFemto_Ghetto_rstar->GetNbins());
  DLM_HistoSource HistoSource(dlmFemto_Ghetto_rstar);


  //#,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
  double PotPars1S0[8]={NN_AV18,v18_Coupled3P2,1,1,1,0,0,0};
  double PotPars3P0[8]={NN_AV18,v18_Coupled3P2,1,1,1,1,1,0};
  double PotPars3P1[8]={NN_AV18,v18_Coupled3P2,1,1,1,1,1,1};
  double PotPars3P2[8]={NN_AV18,v18_Coupled3P2,1,1,1,1,1,2};
  double PotPars1D2[8]={NN_AV18,v18_Coupled3P2,1,1,1,0,2,2};
  CATSparameters cPotPars1S0(CATSparameters::tPotential,8,true); cPotPars1S0.SetParameters(PotPars1S0);
  CATSparameters cPotPars3P0(CATSparameters::tPotential,8,true); cPotPars3P0.SetParameters(PotPars3P0);
  CATSparameters cPotPars3P1(CATSparameters::tPotential,8,true); cPotPars3P1.SetParameters(PotPars3P1);
  CATSparameters cPotPars3P2(CATSparameters::tPotential,8,true); cPotPars3P2.SetParameters(PotPars3P2);
  CATSparameters cPotPars1D2(CATSparameters::tPotential,8,true); cPotPars1D2.SetParameters(PotPars1D2);
  


  const unsigned NumMomBins_pp = 100;
  const double kMin_pp = 0;
  const double kMax_pp = 400;
  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles/",GetCernBoxDimi()));

  CATS pp_traditional;
  pp_traditional.SetMomBins(NumMomBins_pp,kMin_pp,kMax_pp);
  AnalysisObject.SetUpCats_pp(pp_traditional,"AV18","McGauss_ResoTM",0,202);
  pp_traditional.SetAnaSource(0,r_TR);
  pp_traditional.KillTheCat();

  CATS pp_new;
  pp_new.SetMomBins(NumMomBins_pp,kMin_pp,kMax_pp);
  pp_new.SetThetaDependentSource(false);
  pp_new.SetMomentumDependentSource(false);
  pp_new.SetAnaSource(CatsSourceForwarder, &HistoSource, 0);
  pp_new.SetUseAnalyticSource(true);
  pp_new.SetExcludeFailedBins(false);
  pp_new.SetQ1Q2(1);
  pp_new.SetPdgId(2212, 2212);
  pp_new.SetRedMass( 0.5*Mass_p );
  pp_new.SetNumChannels(4);
  pp_new.SetNumPW(0,3);
  pp_new.SetNumPW(1,2);
  pp_new.SetNumPW(2,2);
  pp_new.SetNumPW(3,2);
  pp_new.SetSpin(0,0);
  pp_new.SetSpin(1,1);
  pp_new.SetSpin(2,1);
  pp_new.SetSpin(3,1);
  pp_new.SetChannelWeight(0, 3./12.);
  pp_new.SetChannelWeight(1, 1./12.);
  pp_new.SetChannelWeight(2, 3./12.);
  pp_new.SetChannelWeight(3, 5./12.);
  pp_new.SetShortRangePotential(0,0,fDlmPot,cPotPars1S0);
  pp_new.SetShortRangePotential(0,2,fDlmPot,cPotPars1D2);
  pp_new.SetShortRangePotential(1,1,fDlmPot,cPotPars3P0);
  pp_new.SetShortRangePotential(2,1,fDlmPot,cPotPars3P1);
  pp_new.SetShortRangePotential(3,1,fDlmPot,cPotPars3P2);
  pp_new.KillTheCat();


  CATS pp_old;
  pp_old.SetMomBins(NumMomBins_pp,kMin_pp,kMax_pp);
  pp_old.SetThetaDependentSource(false);
  pp_old.SetMomentumDependentSource(false);
  pp_old.SetAnaSource(CatsSourceForwarder, GT1.Old_source, 2);
  pp_old.SetAnaSource(0, fSource_core->GetParameter(0));
  pp_old.SetAnaSource(1, 2.0);
  pp_old.SetUseAnalyticSource(true);
  pp_old.SetExcludeFailedBins(false);
  pp_old.SetQ1Q2(1);
  pp_old.SetPdgId(2212, 2212);
  pp_old.SetRedMass( 0.5*Mass_p );
  pp_old.SetNumChannels(4);
  pp_old.SetNumPW(0,3);
  pp_old.SetNumPW(1,2);
  pp_old.SetNumPW(2,2);
  pp_old.SetNumPW(3,2);
  pp_old.SetSpin(0,0);
  pp_old.SetSpin(1,1);
  pp_old.SetSpin(2,1);
  pp_old.SetSpin(3,1);
  pp_old.SetChannelWeight(0, 3./12.);
  pp_old.SetChannelWeight(1, 1./12.);
  pp_old.SetChannelWeight(2, 3./12.);
  pp_old.SetChannelWeight(3, 5./12.);
  pp_old.SetShortRangePotential(0,0,fDlmPot,cPotPars1S0);
  pp_old.SetShortRangePotential(0,2,fDlmPot,cPotPars1D2);
  pp_old.SetShortRangePotential(1,1,fDlmPot,cPotPars3P0);
  pp_old.SetShortRangePotential(2,1,fDlmPot,cPotPars3P1);
  pp_old.SetShortRangePotential(3,1,fDlmPot,cPotPars3P2);
  pp_old.KillTheCat();


  TH1F* hTraditional_rstar = new TH1F("hTraditional_rstar","hTraditional_rstar",512,0,32);
  for(unsigned uBin=0; uBin<hTraditional_rstar->GetNbinsX(); uBin++){
    double rad = hTraditional_rstar->GetBinCenter(uBin+1);
    hTraditional_rstar->SetBinContent(uBin+1,pp_traditional.EvaluateTheSource(10,rad,0));
  }

  TH1F* Ck_pp_trad = new TH1F("Ck_pp_trad","Ck_pp_trad",NumMomBins_pp,kMin_pp,kMax_pp);
  TH1F* Ck_pp_new = new TH1F("Ck_pp_new","Ck_pp_new",NumMomBins_pp,kMin_pp,kMax_pp);
  TH1F* Ck_pp_old = new TH1F("Ck_pp_old","Ck_pp_old",NumMomBins_pp,kMin_pp,kMax_pp);
  for(unsigned uBin=0; uBin<NumMomBins_pp; uBin++){
    Ck_pp_trad->SetBinContent(uBin+1,pp_traditional.GetCorrFun(uBin));
    Ck_pp_new->SetBinContent(uBin+1,pp_new.GetCorrFun(uBin));
    Ck_pp_old->SetBinContent(uBin+1,pp_old.GetCorrFun(uBin));
  }
  TH1F* Ck_pp_ratio = (TH1F*)Ck_pp_new->Clone("Ck_pp_ratio");
  Ck_pp_ratio->Divide(Ck_pp_old);

  fOutput.cd();
  h_kstar->Write();
  h_Ghetto_rstar->Write();
  h_Ghetto_rcore->Write();
  hFemto_Ghetto_rstar->Write();
  fSource_star->Write();
  hTraditional_rstar->Write();
  hGhettOld_rstar->Write();
  hFemto_Old_rstar->Write();
  hFemto_OldNew_rstar->Write();
  hFemto_Ghetto_rcore->Write();
  fSource_core->Write();
  hFemto_Old_rcore->Write();
  hFemto_OldNew_rcore->Write();
  hFemto_Old_CosRcP1->Write();
  hFemto_Old_CosRcP2->Write();
  hFemto_Old_CosP1P2->Write();
  hFemto_Old_RcP1->Write();
  hFemto_Old_RcP2->Write();
  hFemto_Old_P1P2->Write();
  Ck_pp_trad->Write();
  Ck_pp_new->Write();
  Ck_pp_old->Write();
  Ck_pp_ratio->Write();


  delete h_Ghetto_rstar;
  delete h_Ghetto_rcore;

  delete hFemto_Ghetto_rstar;
  delete hFemto_Ghetto_rcore;

  delete fSource_star;
  delete fSource_core;

  delete hFemto_Old_CosRcP1;
  delete hFemto_Old_CosRcP2;
  delete hFemto_Old_CosP1P2;
  //delete dlmFemto_Ghetto_rstar;

}

void RandomParallel(const unsigned& NumThr){

  const unsigned MaxThreads = 4;

  TH1F** hDist = new TH1F* [MaxThreads+1];
  for(unsigned uThr=0; uThr<=MaxThreads; uThr++){
    hDist[uThr] = new TH1F(TString::Format("hDist_%u_%u",uThr,NumThr),TString::Format("hDist_%u_%u",uThr,NumThr),2048,0,1000);
  }

  omp_set_dynamic(0);
  omp_set_num_threads(NumThr);
  DLM_Random** RanGen = new DLM_Random* [MaxThreads];
  for(unsigned uThr=0; uThr<MaxThreads; uThr++){
    RanGen[uThr] = new DLM_Random(uThr+1);
  }
  std::vector<float>* rnd_num = new std::vector<float>[MaxThreads];

  const unsigned NumIter = 1000*1000*200;
  const unsigned PlotIter = 16;
  #pragma omp parallel for
  for(unsigned uIter=0; uIter<NumIter; uIter++){
    unsigned ThId = omp_get_thread_num();
    unsigned ThN = omp_get_num_threads();
    if(ThN==1) {ThId = uIter%MaxThreads;}
    rnd_num[ThId].push_back(RanGen[ThId]->Gauss(500,50));
  }

  printf("#  ");
  for(unsigned uThr=0; uThr<MaxThreads; uThr++){
    if(uThr<10) printf("   T%u",uThr);
    else printf("  T%u",uThr);
  }
  printf("\n");
  for(unsigned pIter=0; pIter<PlotIter; pIter++){
    printf("%2.0u:",pIter);
    for(unsigned uThr=0; uThr<MaxThreads; uThr++){
      printf("%5.0f",rnd_num[uThr].at(pIter));
    }
    printf("\n");
  }

  for(unsigned uThr=0; uThr<MaxThreads; uThr++){
    for(float rnd : rnd_num[uThr]){
      hDist[uThr]->Fill(rnd);
      hDist[MaxThreads]->Fill(rnd);
    }
  }

  TFile fOutput(TString::Format("%s/FunWithCeca/RandomParallel_%u.root",GetFemtoOutputFolder(),NumThr),"recreate");
  for(unsigned uThr=0; uThr<=MaxThreads; uThr++){
    hDist[uThr]->Write();
    delete hDist[uThr];
  }
  delete [] hDist;

  for(unsigned uThr=0; uThr<MaxThreads; uThr++){
    delete RanGen[uThr];
  }
  delete [] RanGen;
  delete [] rnd_num;
}

void Test3(){
  std::string ListOfParticles;
  ListOfParticles = "proton,proton;,;Lambda";
  std::vector<std::string> Particles = ParseString(ListOfParticles," ");
  std::vector<std::string> Particles2;
  for(std::string str : Particles){
    std::cout << "|" << str << "|" << std::endl;
    std::vector<std::string> tmp = ParseString(str,",");
    for(std::string str2 : tmp){
      Particles2.push_back(str2);
      std::cout << "   " << str2 << std::endl;
    }
  }
  for(std::string str : Particles2){
    std::cout << "*" << str << "*" << std::endl;
  }

  printf("BOOST\n");
  //std::string inputString("One!Two,Three:Four");
  std::string delimiters(",;");
  std::vector<std::string> parts;
  boost::split(parts, ListOfParticles, boost::is_any_of(delimiters));
  for(std::string str : parts){
    std::cout << "*" << str << "*" << std::endl;
  }

}

//for the fun, each pT will be sigma*mass (Fig_spectra_piKp_pp13_ANA.png in Pictures)
void CecaTest_0(){

  const double Sigma = 1.0;

  TREPNI Database(0);
  Database.SetSeed(11);
  std::vector<TreParticle*> ParticleList;
  ParticleList.push_back(Database.NewParticle("pim"));
  ParticleList.push_back(Database.NewParticle("pip"));
  ParticleList.push_back(Database.NewParticle("pi0"));
  ParticleList.push_back(Database.NewParticle("p"));
  //ParticleList.push_back(Database.NewParticle("a_p"));
  ParticleList.push_back(Database.NewParticle("n"));
  //ParticleList.push_back(Database.NewParticle("a_n"));

  ParticleList.push_back(Database.NewParticle("rho"));//pi pi
  ParticleList.push_back(Database.NewParticle("D0_1232"));//N pi
  ParticleList.push_back(Database.NewParticle("Dp_1232"));//N pi
  ParticleList.push_back(Database.NewParticle("Dpp_1232"));//N pi
  //ParticleList.push_back(Database.NewParticle("a_rho"));
  //ParticleList.push_back(Database.NewParticle("a_D1232"));


  for(TreParticle* prt : ParticleList){
    if(prt->GetName()=="pim"||prt->GetName()=="pip"||prt->GetName()=="pi0"){
      prt->SetMass(140);
      prt->SetAbundance(100);
    }
    else if(prt->GetName()=="p"||prt->GetName()=="n"){
      prt->SetMass(940);
      prt->SetAbundance(20);
    }
    else if(prt->GetName()=="rho"){
      prt->SetMass(770);
      prt->SetAbundance(25);
      prt->SetWidth(149);

      prt->NewDecay();
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("pim"));
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("pip"));
      prt->GetDecay(0)->SetBranching(50);

      prt->NewDecay();
      prt->GetDecay(1)->AddDaughter(*Database.GetParticle("pi0"));
      prt->GetDecay(1)->AddDaughter(*Database.GetParticle("pi0"));
      prt->GetDecay(1)->SetBranching(50);
    }
    else if(prt->GetName()=="D0_1232"||prt->GetName()=="Dp_1232"||prt->GetName()=="Dpp_1232"){
      prt->SetMass(1232);
      prt->SetAbundance(15);
      prt->SetWidth(117);
    }


    if(prt->GetName()=="D0_1232"){
      prt->NewDecay();
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("p"));
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("pim"));
      prt->GetDecay(0)->SetBranching(50);

      prt->NewDecay();
      prt->GetDecay(1)->AddDaughter(*Database.GetParticle("n"));
      prt->GetDecay(1)->AddDaughter(*Database.GetParticle("pi0"));
      prt->GetDecay(1)->SetBranching(50);
    }
    else if(prt->GetName()=="Dp_1232"){
      prt->NewDecay();
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("p"));
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("pi0"));
      prt->GetDecay(0)->SetBranching(50);

      prt->NewDecay();
      prt->GetDecay(1)->AddDaughter(*Database.GetParticle("n"));
      prt->GetDecay(1)->AddDaughter(*Database.GetParticle("pip"));
      prt->GetDecay(1)->SetBranching(50);
    }
    else if(prt->GetName()=="Dpp_1232"){
      prt->NewDecay();
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("p"));
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("pip"));
      prt->GetDecay(0)->SetBranching(100);
    }

    prt->SetPtPz(Sigma*prt->GetMass(),Sigma*prt->GetMass());
  }

  std::vector<std::string> ListOfParticles;
  ListOfParticles.push_back("p");
  ListOfParticles.push_back("p");
  ListOfParticles.push_back("pim");
  CECA Ivana(Database,ListOfParticles);
  Ivana.SetDisplacement(0.1);
  Ivana.SetHadronization(0.3);
  Ivana.SetTau(0.5);
  Ivana.SetTargetStatistics(10);
  Ivana.SetEventMult(3);
  Ivana.SetSourceDim(3);
  Ivana.SetDebugMode(true);
  Ivana.SetThreadTimeout(5);
  Ivana.GHETTO_EVENT = true;
  Ivana.GoBabyGo(1);

}


void CecaTest_1(){

  const double Sigma = 0.6;
  const double LambdaFraction = 1.0;
  const double HadronSize = 0.75*0;
  const double HadronSlope = 0.2*0;

  TREPNI Database(0);
  Database.SetSeed(11);
  std::vector<TreParticle*> ParticleList;
  ParticleList.push_back(Database.NewParticle("Proton"));
  ParticleList.push_back(Database.NewParticle("Lambda"));
  ParticleList.push_back(Database.NewParticle("Pion"));

  ParticleList.push_back(Database.NewParticle("ProtonReso"));
  ParticleList.push_back(Database.NewParticle("LambdaReso"));

  double ExpPP,ExpPR,ExpRP,ExpRR;
  ExpPP = 0.3578*0.3562;
  ExpPR = 0.3578*0.6438;
  ExpRP = 0.6422*0.3562;
  ExpRR = 0.6422*0.6438;

  for(TreParticle* prt : ParticleList){
    if(prt->GetName()=="Proton"){
      prt->SetMass(Mass_p);
      prt->SetAbundance(35.78);
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);
    }
    else if(prt->GetName()=="Lambda"){
      prt->SetMass(Mass_L);
      prt->SetAbundance(35.62*LambdaFraction);
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);
    }
    else if(prt->GetName()=="Pion"){
      prt->SetMass(Mass_pic);
      prt->SetAbundance(0);
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);
    }
    else if(prt->GetName()=="ProtonReso"){
      prt->SetMass(1362);
      prt->SetAbundance(64.22*1);
      prt->SetWidth(hbarc/1.65);
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);

      prt->NewDecay();
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Proton"));
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
      prt->GetDecay(0)->SetBranching(100);
    }
    else if(prt->GetName()=="LambdaReso"){
      prt->SetMass(1463);
      prt->SetAbundance(64.38*1);
      prt->SetWidth(hbarc/4.69);
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);

      prt->NewDecay();
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Lambda"));
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
      prt->GetDecay(0)->SetBranching(100);
    }
    prt->SetPtPz(Sigma*prt->GetMass(),Sigma*prt->GetMass());
  }


/*
  HI Volume = 5000 fm^3; r = 10.6 fm
  HI Temp = 154 MeV
  lets assume (TF) pp Temp = 171 MeV
  if we look at engy density vs temp plot, we see that this
  difference in Temp results in c.a. 1.5x difference in energy density (ED)
  since ED = Energy/Volume, and for HI E = 208 * 5.02 GeV = 1000 TeV,
  while for pp it is 2*6.5 = 13 TeV, we can find the ratio between the volume/radius
  in HI and pp collisions. It turns out to be c.a. 14x, so we have 70 fm^3 volume in pp collisions
  this is a radius of 2.55. Now, for an ellipsoid we need 4/3pi 4.8*4.8*0.7 = 70, where the 0.7 is found
  by dividing 4.8 by the gamma factor of the proton.
  Finally, we give the particle a 98.7 MeV kick (in xyz), which corresponds to 171 MeV directional kick
*/

const double VolumeHI = 5000;
const double TempHI = 154;
const double TempPP = 171;//171
const double ThKick = TempPP/sqrt(3.);
const double EngyHI = 208.*5.02e6;
const double ProtonEngy = 6.5e6;
const double EngyPP = 2.*ProtonEngy;
const double EngyDenst_PP_HI = 0.029412*TempPP-3.5294;//ratio
const double VolumePP = EngyPP/EngyHI*VolumeHI/EngyDenst_PP_HI;
const double ProtonGamma = ProtonEngy/Database.GetParticle("Proton")->GetMass();
const double ExpansionZ = 0.3;//will determine the z-axis fireball size
const double ExpansionDiameter = 2.*ExpansionZ;
//const double ExpansionDiameter = 1.2;

//const double a=ExpansionZ*ProtonGamma;
//const double b=VolumePP*ProtonGamma;
//the ellipsoid has a size of a x a x (a/γ+ExpansionZ) fm
//const double EllipsoidSize = 1./3.*(-a+(pow(2.,1./3.)*a*a)/pow(-2.*a*a*a+27.*b+3.*sqrt(3.)*sqrt(-4.*a*a*a*b+27.*b*b),1./3.)+pow(-2.*a*a*a+27.*b+3.*sqrt(3.)*sqrt(-4.*a*a*a*b+27.*b*b),1./3.)/pow(2.,1./3.));;
const double EllipsoidSize = sqrt(VolumePP/4.*3./Pi/ExpansionDiameter);//approximation for very large gamma

const double PancakeT = pow(3.*VolumePP*ProtonGamma/4./Pi,1./3.);
const double PancakeZ = PancakeT/ProtonGamma;

//
printf("The pp volume is %.2f fm^3\n",VolumePP);
printf("γ = %.2f\n",ProtonGamma);
printf("Ellipsoid: %.2f x %.2f x %.2f\n", EllipsoidSize,EllipsoidSize,ExpansionDiameter);
printf("Pancake  : %.2f x %.2f x %.2f\n", PancakeT,PancakeT,PancakeZ);

  std::vector<std::string> ListOfParticles;
  ListOfParticles.push_back("Proton");
  ListOfParticles.push_back("Proton");
  //ListOfParticles.push_back("Lambda");
  CECA Ivana(Database,ListOfParticles);
  //Ivana.SetDisplacementT(0.85*0.0);
  //Ivana.SetDisplacementZ(0.059*0.0);
  //const double HadronEnergy = 6500;
  //const double EmissionScale = 1.0;
  //printf("rh_T = %.2f fm\n",EmissionScale*HadronSize);
  //printf("rh_Z = %.2f fm\n",EmissionScale*HadronSize/ProtonGamma);
  //usleep(1000e3);
  //Ivana.SetHadronizationT(EmissionScale*HadronSize);
  //Ivana.SetHadronizationZ(EmissionScale*HadronSize/ProtonGamma);
  //Ivana.SetHadronizationT(EllipsoidSize);//4.8
  //Ivana.SetHadronizationZ(ExpansionDiameter);//0.7
  //Ivana.SetHadronizationT(8.5);//4.8
  //Ivana.SetHadronizationZ(0.59);//0.7
  //Ivana.SetHadronizationT(0.8);
  //Ivana.SetHadronizationT(0.8/ProtonGamma*2.);
  //Ivana.SetTau(0.8);
  //Ivana.SetThermalKick(ThKick);


//this works perfectly! In terms of getting our angular distributions,
//and source sizes for core and reso. mT scaling sucks
  Ivana.SetDisplacement(1.0);
  Ivana.SetHadronization(0.0);
  Ivana.SetTau(0.0);
  //Ivana.SetThermalKick(98);
  Ivana.SetThermalKick(0);
  Ivana.EqualizeFsiTime(true);

/*
Ivana.SetDisplacementT(0.0);
Ivana.SetDisplacementZ(0.0/ProtonGamma);
Ivana.SetHadronizationT(3.6);
Ivana.SetHadronizationZ(0.77);
Ivana.SetHadrFluctuation(0.00);
Ivana.SetTau(0.0);
Ivana.SetThermalKick(98);
*/
/*
//here the conclusion is: the asymmetry in the source and the HadrFluct are the best parameters to fine tune
//perhaps also the temperature of the system
Ivana.SetDisplacementT(0.0);
Ivana.SetDisplacementZ(0.0/ProtonGamma);
//a larger total volume increases the size of the source.
//this corresponds to lower freeze-out temperature
//the larger the asymmetry here, the larger the source size and the slope of the mT
Ivana.SetHadronizationT(EllipsoidSize);
Ivana.SetHadronizationZ(ExpansionDiameter);
//this can reduce the slope of the mT
//5-15 % are reasonable values, if we want to compare them to the
//ratio between tempreture and particle energy
Ivana.SetHadrFluctuation(0.12);
Ivana.SetTau(0.0);
//this has little influence
Ivana.SetThermalKick(ThKick);
*/
  Ivana.SetTargetStatistics(10);
  Ivana.SetEventMult(2);
  Ivana.SetSourceDim(2);
  Ivana.SetDebugMode(true);
  Ivana.SetThreadTimeout(15);
  Ivana.GHETTO_EVENT = true;


  //!!! some small mem leak might still be there,
  //noticable when running 8 cores
  Ivana.GoBabyGo();

  //Ivana.SetSeed(0,1);
  //Ivana.SetSeed(1,2);
  //Ivana.GoBabyGo(2);

  //Ivana.SetSeed(0,2);
  //Ivana.GoBabyGo(1);



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

  printf("   Expected  Total  Femto\n");
  printf("PP %6.2f%% %6.2f%% %6.2f\n",ExpPP*100.,TotPP*100.,FemtoPP*100.);
  printf("PR %6.2f%% %6.2f%% %6.2f\n",ExpPR*100.,TotPR*100.,FemtoPR*100.);
  printf("RP %6.2f%% %6.2f%% %6.2f\n",ExpRP*100.,TotRP*100.,FemtoRP*100.);
  printf("RR %6.2f%% %6.2f%% %6.2f\n",ExpRR*100.,TotRR*100.,FemtoRR*100.);

  TH1F* h_Ghetto_rstar = Convert_DlmHisto_TH1F(Ivana.Ghetto_rstar,"Ghetto_rstar");
  TH1F* h_Ghetto_kstar = Convert_DlmHisto_TH1F(Ivana.Ghetto_kstar,"Ghetto_kstar");
  TH2F* h_Ghetto_kstar_rstar = Convert_DlmHisto_TH2F(Ivana.Ghetto_kstar_rstar,"Ghetto_kstar_rstar");
  TH2F* h_Ghetto_mT_rstar = Convert_DlmHisto_TH2F(Ivana.Ghetto_mT_rstar,"Ghetto_mT_rstar");
  TH2F* h_GhettoFemto_mT_rstar = Convert_DlmHisto_TH2F(Ivana.GhettoFemto_mT_rstar,"GhettoFemto_mT_rstar");
  TH2F* h_GhettoFemto_mT_kstar = Convert_DlmHisto_TH2F(Ivana.GhettoFemto_mT_kstar,"GhettoFemto_mT_kstar");
  TH2F* h_Ghetto_mT_costh = Convert_DlmHisto_TH2F(Ivana.Ghetto_mT_costh,"Ghetto_mT_costh");
  TH2F* h_GhettoSP_pT_th = Convert_DlmHisto_TH2F(Ivana.GhettoSP_pT_th,"GhettoSP_pT_th");
  Ivana.GhettoFemto_rstar->ComputeError();
  TH1F* h_GhettoFemto_rstar = Convert_DlmHisto_TH1F(Ivana.GhettoFemto_rstar,"GhettoFemto_rstar");
  TH1F* h_Ghetto_RP_AngleRcP1 = Convert_DlmHisto_TH1F(Ivana.Ghetto_RP_AngleRcP1,"Ghetto_RP_AngleRcP1");
  TH1F* h_Ghetto_PR_AngleRcP2 = Convert_DlmHisto_TH1F(Ivana.Ghetto_PR_AngleRcP2,"Ghetto_PR_AngleRcP2");
  TH1F* h_Ghetto_RR_AngleRcP1 = Convert_DlmHisto_TH1F(Ivana.Ghetto_RR_AngleRcP1,"Ghetto_RR_AngleRcP1");
  TH1F* h_Ghetto_RR_AngleRcP2 = Convert_DlmHisto_TH1F(Ivana.Ghetto_RR_AngleRcP2,"Ghetto_RR_AngleRcP2");
  TH1F* h_Ghetto_RR_AngleP1P2 = Convert_DlmHisto_TH1F(Ivana.Ghetto_RR_AngleP1P2,"Ghetto_RR_AngleP1P2");

  TGraphErrors g_GhettoFemto_mT_rstar;
  g_GhettoFemto_mT_rstar.SetName("g_GhettoFemto_mT_rstar");
  g_GhettoFemto_mT_rstar.SetMarkerStyle(20);
  g_GhettoFemto_mT_rstar.SetMarkerSize(1);
  g_GhettoFemto_mT_rstar.SetLineWidth(3);

  unsigned uPoint = 0;
  for(unsigned uBin=0; uBin<h_GhettoFemto_mT_rstar->GetXaxis()->GetNbins(); uBin++){
    TH1D* hProj = h_GhettoFemto_mT_rstar->ProjectionY(TString::Format("hProj"),uBin+1,uBin+1);
    double Mean = hProj->GetMean();
    double Err = hProj->GetStdDev();
    double mT = h_GhettoFemto_mT_rstar->GetXaxis()->GetBinCenter(uBin+1);
    if(Mean&&Err){
      g_GhettoFemto_mT_rstar.SetPoint(uPoint,mT,Mean);
      g_GhettoFemto_mT_rstar.SetPointError(uPoint,0,Err);
      uPoint++;
    }
  }

  const double MeanToGauss = 1./2.3;
  double mT_first = g_GhettoFemto_mT_rstar.GetPointY(0);
  printf("r(%.2f GeV) = %.2f\n",g_GhettoFemto_mT_rstar.GetPointX(0)*0.001,mT_first*MeanToGauss);
  double mT_2GeV = g_GhettoFemto_mT_rstar.Eval(2000);
  printf("r(%.2f GeV) = %.2f\n",2.,mT_2GeV*MeanToGauss);

  h_GhettoFemto_rstar->Scale(1./h_GhettoFemto_rstar->Integral(),"width");

  TF1* fSource = new TF1("fSource",GaussSourceTF1,0,5,1);
  fSource->SetParameter(0,1.2);
  fSource->SetParLimits(0,0.4,2.0);
  h_GhettoFemto_rstar->Fit(fSource,"Q, S, N, R, M");
  printf("The effective Gaussian size is %.3f +/- %.3f fm\n",fSource->GetParameter(0),fSource->GetParError(0));

  TFile fOutput(TString::Format("%s/FunWithCeca/CecaTest_1.root",GetFemtoOutputFolder()),"recreate");
  h_Ghetto_rstar->Write();
  h_GhettoFemto_rstar->Write();
  fSource->Write();
  h_Ghetto_kstar->Write();
  h_Ghetto_kstar_rstar->Write();
  h_Ghetto_mT_rstar->Write();
  h_GhettoFemto_mT_rstar->Write();
  g_GhettoFemto_mT_rstar.Write();
  h_GhettoFemto_mT_kstar->Write();
  h_Ghetto_mT_costh->Write();
  h_GhettoSP_pT_th->Write();
  h_Ghetto_RP_AngleRcP1->Write();
  h_Ghetto_PR_AngleRcP2->Write();
  h_Ghetto_RR_AngleRcP1->Write();
  h_Ghetto_RR_AngleRcP2->Write();
  h_Ghetto_RR_AngleP1P2->Write();

}



void CecaTest_pipi_1(){

  const double Sigma = 0.6;
  const double HadronSize = 0.75;
  const double HadronSlope = 0.2;

  TREPNI Database(0);
  Database.SetSeed(11);
  std::vector<TreParticle*> ParticleList;
  ParticleList.push_back(Database.NewParticle("Pion"));

  ParticleList.push_back(Database.NewParticle("PionReso"));

  double ExpPP,ExpPR,ExpRP,ExpRR;
  ExpPP = 0.318*0.318;
  ExpPR = 0.318*0.682;
  ExpRP = 0.682*0.318;
  ExpRR = 0.682*0.682;

  for(TreParticle* prt : ParticleList){
    if(prt->GetName()=="Pion"){
      prt->SetMass(Mass_pic);
      prt->SetAbundance(31.8);
      prt->SetRadius(HadronSize*2./3.);
      prt->SetRadiusSlope(HadronSlope);
    }
    else if(prt->GetName()=="PionReso"){
      prt->SetMass(1180);
      prt->SetAbundance(68.2*1);
      prt->SetWidth(hbarc/1.50);
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);

      prt->NewDecay();
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
      prt->GetDecay(0)->SetBranching(100);
    }
    prt->SetPtPz(Sigma*prt->GetMass(),Sigma*prt->GetMass());
  }


/*
  HI Volume = 5000 fm^3; r = 10.6 fm
  HI Temp = 154 MeV
  lets assume (TF) pp Temp = 171 MeV
  if we look at engy density vs temp plot, we see that this
  difference in Temp results in c.a. 1.5x difference in energy density (ED)
  since ED = Energy/Volume, and for HI E = 208 * 5.02 GeV = 1000 TeV,
  while for pp it is 2*6.5 = 13 TeV, we can find the ratio between the volume/radius
  in HI and pp collisions. It turns out to be c.a. 14x, so we have 70 fm^3 volume in pp collisions
  this is a radius of 2.55. Now, for an ellipsoid we need 4/3pi 4.8*4.8*0.7 = 70, where the 0.7 is found
  by dividing 4.8 by the gamma factor of the proton.
  Finally, we give the particle a 98.7 MeV kick (in xyz), which corresponds to 171 MeV directional kick
*/

const double VolumeHI = 5000;
const double TempHI = 154;
const double TempPP = 171;//171
const double ThKick = TempPP/sqrt(3.);
const double EngyHI = 208.*5.02e6;
const double ProtonEngy = 6.5e6;
const double EngyPP = 2.*ProtonEngy;
const double EngyDenst_PP_HI = 0.029412*TempPP-3.5294;//ratio
const double VolumePP = EngyPP/EngyHI*VolumeHI/EngyDenst_PP_HI;
const double ProtonGamma = ProtonEngy/938.;
const double ExpansionZ = 0.3;//will determine the z-axis fireball size
const double ExpansionDiameter = 2.*ExpansionZ;
const double EllipsoidSize = sqrt(VolumePP/4.*3./Pi/ExpansionDiameter);//approximation for very large gamma

const double PancakeT = pow(3.*VolumePP*ProtonGamma/4./Pi,1./3.);
const double PancakeZ = PancakeT/ProtonGamma;

printf("The pp volume is %.2f fm^3\n",VolumePP);
printf("γ = %.2f\n",ProtonGamma);
printf("Ellipsoid: %.2f x %.2f x %.2f\n", EllipsoidSize,EllipsoidSize,ExpansionDiameter);
printf("Pancake  : %.2f x %.2f x %.2f\n", PancakeT,PancakeT,PancakeZ);

  std::vector<std::string> ListOfParticles;
  ListOfParticles.push_back("Pion");
  ListOfParticles.push_back("Pion");
  CECA Ivana(Database,ListOfParticles);

/*
  Ivana.SetDisplacementT(0.0);
  Ivana.SetDisplacementZ(0.0/ProtonGamma);
  //Ivana.SetHadronizationT(EllipsoidSize);
  //Ivana.SetHadronizationZ(ExpansionDiameter);
  Ivana.SetHadronizationT(2.15);
  Ivana.SetHadronizationZ(2.15);
  Ivana.SetHadrFluctuation(0.12);
  Ivana.SetTau(0.0);
  Ivana.SetThermalKick(ThKick);
*/

const double Vol = 41.5*1.0;
const double Asymm = 1.0;
const double HadrT = pow(3.*Vol*Asymm/4./Pi,1./3.);
const double HadrZ = HadrT/Asymm;
const double V_QA = 4./3.*Pi*HadrT*HadrT*HadrZ;
printf("----------------------\n");
printf("aT = %.3f fm\n",HadrT);
printf("aZ = %.3f fm\n",HadrZ);
printf(" V = %.3f fm^3\n",V_QA);

Ivana.SetDisplacement(0.0);
Ivana.SetHadronizationT(HadrT);
Ivana.SetHadronizationZ(HadrZ);
Ivana.SetHadrFluctuation(0.0);
Ivana.SetTau(0.0);
Ivana.SetThermalKick(ThKick);

/*
  Ivana.SetDisplacement(1.0);
  Ivana.SetHadronization(0.0);
  Ivana.SetTau(0.0);
  Ivana.SetThermalKick(0);
*/

  Ivana.SetTargetStatistics(10);
  Ivana.SetEventMult(2);
  Ivana.SetSourceDim(2);
  Ivana.SetDebugMode(true);
  Ivana.SetThreadTimeout(30);
  Ivana.GHETTO_EVENT = true;

  Ivana.GoBabyGo(1);

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

  printf("   Expected  Total  Femto\n");
  printf("PP %6.2f%% %6.2f%% %6.2f\n",ExpPP*100.,TotPP*100.,FemtoPP*100.);
  printf("PR %6.2f%% %6.2f%% %6.2f\n",ExpPR*100.,TotPR*100.,FemtoPR*100.);
  printf("RP %6.2f%% %6.2f%% %6.2f\n",ExpRP*100.,TotRP*100.,FemtoRP*100.);
  printf("RR %6.2f%% %6.2f%% %6.2f\n",ExpRR*100.,TotRR*100.,FemtoRR*100.);

  TH1F* h_Ghetto_rstar = Convert_DlmHisto_TH1F(Ivana.Ghetto_rstar,"Ghetto_rstar");
  TH1F* h_Ghetto_kstar = Convert_DlmHisto_TH1F(Ivana.Ghetto_kstar,"Ghetto_kstar");
  TH2F* h_Ghetto_kstar_rstar = Convert_DlmHisto_TH2F(Ivana.Ghetto_kstar_rstar,"Ghetto_kstar_rstar");
  TH2F* h_Ghetto_mT_rstar = Convert_DlmHisto_TH2F(Ivana.Ghetto_mT_rstar,"Ghetto_mT_rstar");
  TH2F* h_GhettoFemto_mT_rstar = Convert_DlmHisto_TH2F(Ivana.GhettoFemto_mT_rstar,"GhettoFemto_mT_rstar");
  TH2F* h_GhettoFemto_mT_kstar = Convert_DlmHisto_TH2F(Ivana.GhettoFemto_mT_kstar,"GhettoFemto_mT_kstar");
  TH2F* h_Ghetto_mT_costh = Convert_DlmHisto_TH2F(Ivana.Ghetto_mT_costh,"Ghetto_mT_costh");
  TH2F* h_GhettoSP_pT_th = Convert_DlmHisto_TH2F(Ivana.GhettoSP_pT_th,"GhettoSP_pT_th");
  Ivana.GhettoFemto_rstar->ComputeError();
  TH1F* h_GhettoFemto_rstar = Convert_DlmHisto_TH1F(Ivana.GhettoFemto_rstar,"GhettoFemto_rstar");
  TH1F* h_Ghetto_RP_AngleRcP1 = Convert_DlmHisto_TH1F(Ivana.Ghetto_RP_AngleRcP1,"Ghetto_RP_AngleRcP1");
  TH1F* h_Ghetto_PR_AngleRcP2 = Convert_DlmHisto_TH1F(Ivana.Ghetto_PR_AngleRcP2,"Ghetto_PR_AngleRcP2");
  TH1F* h_Ghetto_RR_AngleRcP1 = Convert_DlmHisto_TH1F(Ivana.Ghetto_RR_AngleRcP1,"Ghetto_RR_AngleRcP1");
  TH1F* h_Ghetto_RR_AngleRcP2 = Convert_DlmHisto_TH1F(Ivana.Ghetto_RR_AngleRcP2,"Ghetto_RR_AngleRcP2");
  TH1F* h_Ghetto_RR_AngleP1P2 = Convert_DlmHisto_TH1F(Ivana.Ghetto_RR_AngleP1P2,"Ghetto_RR_AngleP1P2");

  TGraphErrors g_GhettoFemto_mT_rstar;
  g_GhettoFemto_mT_rstar.SetName("g_GhettoFemto_mT_rstar");
  g_GhettoFemto_mT_rstar.SetMarkerStyle(20);
  g_GhettoFemto_mT_rstar.SetMarkerSize(1);
  g_GhettoFemto_mT_rstar.SetLineWidth(3);

  unsigned uPoint = 0;
  for(unsigned uBin=0; uBin<h_GhettoFemto_mT_rstar->GetXaxis()->GetNbins(); uBin++){
    TH1D* hProj = h_GhettoFemto_mT_rstar->ProjectionY(TString::Format("hProj"),uBin+1,uBin+1);
    double Mean = hProj->GetMean();
    double Err = hProj->GetStdDev();
    double mT = h_GhettoFemto_mT_rstar->GetXaxis()->GetBinCenter(uBin+1);
    if(Mean&&Err){
      g_GhettoFemto_mT_rstar.SetPoint(uPoint,mT,Mean);
      g_GhettoFemto_mT_rstar.SetPointError(uPoint,0,Err);
      uPoint++;
    }
  }

  const double MeanToGauss = 1./2.3;
  double mT_first = g_GhettoFemto_mT_rstar.GetPointY(0);
  printf("r(%.2f GeV) = %.2f\n",g_GhettoFemto_mT_rstar.GetPointX(0)*0.001,mT_first*MeanToGauss);
  double mT_2GeV = g_GhettoFemto_mT_rstar.Eval(1000);
  printf("r(%.2f GeV) = %.2f\n",1.,mT_2GeV*MeanToGauss);

  h_GhettoFemto_rstar->Scale(1./h_GhettoFemto_rstar->Integral(),"width");

  TF1* fSource = new TF1("fSource",GaussSourceTF1,0,5,1);
  fSource->SetParameter(0,1.2);
  fSource->SetParLimits(0,0.4,3.0);
  h_GhettoFemto_rstar->Fit(fSource,"Q, S, N, R, M");
  printf("The effective Gaussian size is %.3f +/- %.3f fm\n",fSource->GetParameter(0),fSource->GetParError(0));

  TFile fOutput(TString::Format("%s/FunWithCeca/CecaTest_pipi_1.root",GetFemtoOutputFolder()),"recreate");
  h_Ghetto_rstar->Write();
  h_GhettoFemto_rstar->Write();
  fSource->Write();
  h_Ghetto_kstar->Write();
  h_Ghetto_kstar_rstar->Write();
  h_Ghetto_mT_rstar->Write();
  h_GhettoFemto_mT_rstar->Write();
  g_GhettoFemto_mT_rstar.Write();
  h_GhettoFemto_mT_kstar->Write();
  h_Ghetto_mT_costh->Write();
  h_GhettoSP_pT_th->Write();
  h_Ghetto_RP_AngleRcP1->Write();
  h_Ghetto_PR_AngleRcP2->Write();
  h_Ghetto_RR_AngleRcP1->Write();
  h_Ghetto_RR_AngleRcP2->Write();
  h_Ghetto_RR_AngleP1P2->Write();

}


double ScaledGauss(double* x, double* par){
  return par[1]*GaussSourceTF1(x,par);
}

void CecaTest_pd_1(){

  const double Sigma = 0.6;
  const double HadronSize = 0.75;
  const double HadronSlope = 0.2;

  TREPNI Database(0);
  Database.SetSeed(11);
  std::vector<TreParticle*> ParticleList;
  ParticleList.push_back(Database.NewParticle("Proton"));
  ParticleList.push_back(Database.NewParticle("Deuteron"));
  ParticleList.push_back(Database.NewParticle("Pion"));
  ParticleList.push_back(Database.NewParticle("Photon"));
  ParticleList.push_back(Database.NewParticle("ProtonReso"));
  ParticleList.push_back(Database.NewParticle("DeuteronReso"));

  double ExpPP,ExpPR,ExpRP,ExpRR;
  //ExpPP = 0.3578*0.3562;
  //ExpPR = 0.3578*0.6438;
  //ExpRP = 0.6422*0.3562;
  //ExpRR = 0.6422*0.6438;

  for(TreParticle* prt : ParticleList){
    if(prt->GetName()=="Proton"){
      prt->SetMass(Mass_p);
      prt->SetAbundance(35.78);
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);
    }
    else if(prt->GetName()=="Deuteron"){
      prt->SetMass(Mass_d);
      //prt->SetAbundance(35.78);
      prt->SetAbundance(100.0);
      prt->SetRadius(HadronSize*1.5);
      prt->SetRadiusSlope(HadronSlope);
      prt->SetDelayTau(0);
    }
    else if(prt->GetName()=="Pion"){
      prt->SetMass(Mass_pic);
      prt->SetAbundance(0);
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);
    }
    else if(prt->GetName()=="ProtonReso"){
      prt->SetMass(1362);
      prt->SetAbundance(64.22*1);
      prt->SetWidth(hbarc/1.65);
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);

      prt->NewDecay();
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Proton"));
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
      prt->GetDecay(0)->SetBranching(100);
    }
    else if(prt->GetName()=="DeuteronReso"){
      prt->SetMass(Mass_d+0.01);
      prt->SetAbundance(0.0);
      prt->SetWidth(hbarc/100.);
      prt->SetRadius(HadronSize*1.5);
      prt->SetRadiusSlope(HadronSlope);

      prt->NewDecay();
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Deuteron"));
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Photon"));
      prt->GetDecay(0)->SetBranching(100);
    }
    prt->SetPtPz(Sigma*prt->GetMass(),Sigma*prt->GetMass());
  }

  std::vector<std::string> ListOfParticles;
  ListOfParticles.push_back("Proton");
  ListOfParticles.push_back("Deuteron");
  CECA Ivana(Database,ListOfParticles);

//1.0 = r_core of 1.15 fm
//0.8 = r_core of 0.92 fm
  Ivana.SetDisplacement(0.905);
  Ivana.SetHadronization(0.0);
  Ivana.SetHadrFluctuation(0.0);
  Ivana.SetTau(0.0);
  Ivana.SetThermalKick(0);

  Ivana.SetTargetStatistics(10);
  Ivana.SetEventMult(2);
  Ivana.SetSourceDim(2);
  Ivana.SetDebugMode(true);
  Ivana.SetThreadTimeout(15);
  Ivana.EqualizeFsiTime(false);
  Ivana.GHETTO_EVENT = true;

  Ivana.GoBabyGo();

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

  printf("   Expected  Total  Femto\n");
  printf("PP %6.2f%% %6.2f%% %6.2f\n",ExpPP*100.,TotPP*100.,FemtoPP*100.);
  printf("PR %6.2f%% %6.2f%% %6.2f\n",ExpPR*100.,TotPR*100.,FemtoPR*100.);
  printf("RP %6.2f%% %6.2f%% %6.2f\n",ExpRP*100.,TotRP*100.,FemtoRP*100.);
  printf("RR %6.2f%% %6.2f%% %6.2f\n",ExpRR*100.,TotRR*100.,FemtoRR*100.);

  TH1F* h_Ghetto_rstar = Convert_DlmHisto_TH1F(Ivana.Ghetto_rstar,"Ghetto_rstar");
  TH1F* h_Ghetto_kstar = Convert_DlmHisto_TH1F(Ivana.Ghetto_kstar,"Ghetto_kstar");
  Ivana.Ghetto_kstar_rstar->ComputeError();
  TH2F* h_Ghetto_kstar_rstar = Convert_DlmHisto_TH2F(Ivana.Ghetto_kstar_rstar,"Ghetto_kstar_rstar");
  TH2F* h_Ghetto_mT_rstar = Convert_DlmHisto_TH2F(Ivana.Ghetto_mT_rstar,"Ghetto_mT_rstar");
  TH2F* h_GhettoFemto_mT_rstar = Convert_DlmHisto_TH2F(Ivana.GhettoFemto_mT_rstar,"GhettoFemto_mT_rstar");
  TH2F* h_GhettoFemto_mT_kstar = Convert_DlmHisto_TH2F(Ivana.GhettoFemto_mT_kstar,"GhettoFemto_mT_kstar");
  TH2F* h_Ghetto_mT_costh = Convert_DlmHisto_TH2F(Ivana.Ghetto_mT_costh,"Ghetto_mT_costh");
  TH2F* h_GhettoSP_pT_th = Convert_DlmHisto_TH2F(Ivana.GhettoSP_pT_th,"GhettoSP_pT_th");
  TH1F* h_GhettoSP_pT_1 = Convert_DlmHisto_TH1F(Ivana.GhettoSP_pT_1,"GhettoSP_pT_1");
  TH1F* h_GhettoSP_pT_2 = Convert_DlmHisto_TH1F(Ivana.GhettoSP_pT_2,"GhettoSP_pT_2");
  Ivana.GhettoFemto_rstar->ComputeError();
  TH1F* h_GhettoFemto_rstar = Convert_DlmHisto_TH1F(Ivana.GhettoFemto_rstar,"GhettoFemto_rstar");
  TH1F* h_Ghetto_RP_AngleRcP1 = Convert_DlmHisto_TH1F(Ivana.Ghetto_RP_AngleRcP1,"Ghetto_RP_AngleRcP1");
  TH1F* h_Ghetto_PR_AngleRcP2 = Convert_DlmHisto_TH1F(Ivana.Ghetto_PR_AngleRcP2,"Ghetto_PR_AngleRcP2");
  TH1F* h_Ghetto_RR_AngleRcP1 = Convert_DlmHisto_TH1F(Ivana.Ghetto_RR_AngleRcP1,"Ghetto_RR_AngleRcP1");
  TH1F* h_Ghetto_RR_AngleRcP2 = Convert_DlmHisto_TH1F(Ivana.Ghetto_RR_AngleRcP2,"Ghetto_RR_AngleRcP2");
  TH1F* h_Ghetto_RR_AngleP1P2 = Convert_DlmHisto_TH1F(Ivana.Ghetto_RR_AngleP1P2,"Ghetto_RR_AngleP1P2");

  TGraphErrors g_GhettoFemto_mT_rstar;
  g_GhettoFemto_mT_rstar.SetName("g_GhettoFemto_mT_rstar");
  g_GhettoFemto_mT_rstar.SetMarkerStyle(20);
  g_GhettoFemto_mT_rstar.SetMarkerSize(1);
  g_GhettoFemto_mT_rstar.SetLineWidth(3);

  unsigned uPoint = 0;
  for(unsigned uBin=0; uBin<h_GhettoFemto_mT_rstar->GetXaxis()->GetNbins(); uBin++){
    TH1D* hProj = h_GhettoFemto_mT_rstar->ProjectionY(TString::Format("hProj"),uBin+1,uBin+1);
    double Mean = hProj->GetMean();
    double Err = hProj->GetStdDev();
    double mT = h_GhettoFemto_mT_rstar->GetXaxis()->GetBinCenter(uBin+1);
    if(Mean&&Err){
      g_GhettoFemto_mT_rstar.SetPoint(uPoint,mT,Mean);
      g_GhettoFemto_mT_rstar.SetPointError(uPoint,0,Err);
      uPoint++;
    }
  }

  const double MeanToGauss = 1./2.3;
  double mT_first = g_GhettoFemto_mT_rstar.GetPointY(0);
  printf("r(%.2f GeV) = %.2f\n",g_GhettoFemto_mT_rstar.GetPointX(0)*0.001,mT_first*MeanToGauss);
  double mT_2GeV = g_GhettoFemto_mT_rstar.Eval(2000);
  printf("r(%.2f GeV) = %.2f\n",2.,mT_2GeV*MeanToGauss);

  h_GhettoFemto_rstar->Scale(1./h_GhettoFemto_rstar->Integral(),"width");

  double lowerlimit;
  double upperlimit;
  GetCentralInterval(*h_GhettoFemto_rstar, 0.9, lowerlimit, upperlimit, true);
  unsigned lowerbin = h_GhettoFemto_rstar->FindBin(lowerlimit);
  unsigned upperbin = h_GhettoFemto_rstar->FindBin(upperlimit);
  printf("The fit will be performed in the range [%.2f, %.2f] fm\n",lowerlimit,upperlimit);
  TF1* fSource = new TF1("fSource",ScaledGauss,lowerlimit,upperlimit,2);
  fSource->SetParameter(0,h_GhettoFemto_rstar->GetMean()/2.3);
  fSource->SetParLimits(0,0.5,5.0);
  //fSource->SetParameter(1,0.5);
  //fSource->SetParLimits(1,0.1,1.0);
  fSource->FixParameter(1,1.0);

  h_GhettoFemto_rstar->Fit(fSource,"Q, S, N, R, M");
  printf("The effective Gaussian size is %.3f +/- %.3f fm, λ = %.3f +/- %.3f\n",
  fSource->GetParameter(0),fSource->GetParError(0),fSource->GetParameter(1),fSource->GetParError(1));

  TGraph gRadKstar;
  gRadKstar.SetName("gRadKstar");
  gRadKstar.SetLineColor(kRed+1);
  gRadKstar.SetLineWidth(3);

  TH1D** hkstar_rstar;
  hkstar_rstar = new TH1D* [h_Ghetto_kstar_rstar->GetXaxis()->GetNbins()];
  for(unsigned uMom=0; uMom<h_Ghetto_kstar_rstar->GetXaxis()->GetNbins(); uMom++){
    double kstar = h_Ghetto_kstar_rstar->GetXaxis()->GetBinCenter(uMom+1);
    hkstar_rstar[uMom] = NULL;
    if(kstar>150) continue;
    hkstar_rstar[uMom] = h_Ghetto_kstar_rstar->ProjectionY(TString::Format("hkstar_rstar_%.0f",kstar),uMom+1,uMom+1);
    hkstar_rstar[uMom]->Scale(1./hkstar_rstar[uMom]->Integral(),"width");
    GetCentralInterval(*hkstar_rstar[uMom], 0.9, lowerlimit, upperlimit, true);
    //printf("--- k* = %.0f ---\n",kstar);
    //printf(" The fit will be performed in the range [%.2f, %.2f] fm\n",lowerlimit,upperlimit);
    fSource->SetParameter(0,hkstar_rstar[uMom]->GetMean()/2.3);
    fSource->SetParLimits(0,hkstar_rstar[uMom]->GetMean()/4.,hkstar_rstar[uMom]->GetMean());
    hkstar_rstar[uMom]->Fit(fSource,"Q, S, N, R, M");
    //printf(" The effective Gaussian size is %.3f +/- %.3f fm, λ = %.3f +/- %.3f\n",
    //fSource->GetParameter(0),fSource->GetParError(0),fSource->GetParameter(1),fSource->GetParError(1));
    //gRadKstar.SetPoint(uMom,kstar,fSource->GetParameter(0));
  }

  TFile fOutput(TString::Format("%s/FunWithCeca/CecaTest_pd_1.root",GetFemtoOutputFolder()),"recreate");
  h_Ghetto_rstar->Write();
  h_GhettoFemto_rstar->Write();
  fSource->Write();
  h_Ghetto_kstar->Write();
  h_Ghetto_kstar_rstar->Write();
  h_Ghetto_mT_rstar->Write();
  h_GhettoFemto_mT_rstar->Write();
  g_GhettoFemto_mT_rstar.Write();
  h_GhettoFemto_mT_kstar->Write();
  h_Ghetto_mT_costh->Write();
  h_GhettoSP_pT_th->Write();
  h_GhettoSP_pT_1->Write();
  h_GhettoSP_pT_2->Write();
  h_Ghetto_RP_AngleRcP1->Write();
  h_Ghetto_PR_AngleRcP2->Write();
  h_Ghetto_RR_AngleRcP1->Write();
  h_Ghetto_RR_AngleRcP2->Write();
  h_Ghetto_RR_AngleP1P2->Write();
  gRadKstar.Write();
  for(unsigned uMom=0; uMom<h_Ghetto_kstar_rstar->GetXaxis()->GetNbins(); uMom++){
    if(hkstar_rstar[uMom]){
      hkstar_rstar[uMom]->Write();
    }
  }

}

double GaussSourceMean(double* x, double* Pars){
    double& size = Pars[0];
    TF1 f_gauss("f_gauss","[0]*4.*TMath::Pi()*x*x*pow(4.*TMath::Pi()*[1]*[1],-1.5)*exp(-(x*x)/(4.*[1]*[1]))+1.-[0]",0,256);
    f_gauss.FixParameter(0,1);
    f_gauss.FixParameter(1,size);
    return f_gauss.Mean(0,256);
}

double GaussFromMean(const double mean){
  TH1F hHist("hHist","hHist",1,0,1);
  hHist.SetBinContent(1,mean);
  hHist.SetBinError(1,mean*0.01);
  TF1 fHist("fHist",GaussSourceMean,0,1,1);
  fHist.SetParameter(0,mean);
  fHist.SetParLimits(0,mean*0.1,mean*2.);
  hHist.Fit(&fHist,"Q, S, N, R, M");
  return fHist.GetParameter(0);
}

//0.97 -> 1.059 (pd) + 0.04 error => 1.10
//1.04 -> 1.216 (Kd)  upper => 1.314
//we can use the same for pd pr Kd, type == pd or Kd
//N.B. KdReso is as Kd, but the deuteron has an effective resonances, to kind of mimic Delta -> pi p (+n +K) = Kd
//EffFix - for Dpi we can set it to negative value, in which case we have Jaime's core
//pP is essentially that we compare two protons, but one of them is primordial only (useful for pd stuff)
void Ceca_pd_1(const double& d_delay, const int& EffFix, const TString type="pd"){

  if(type!="pp"&&type!="pL"&&type!="pd"&&type!="Kd"&&type!="KdReso"&&type!="Dpi"&&type!="pipi"&&type!="pipi_core"&&type!="pP"&&type!="pK"){
    printf("WHAT IS THIS: %s\n", type.Data());
    return;
  }
  double HadronSize = 0;//0.75
  double HadronSlope = 0;//0.2
  const double EtaCut = 0.8;
  const bool PROTON_RESO = true;
  const double frac_protons = 35.78;
  const double frac_kaons = 52.4*1.1;
  const bool EQUALIZE_TAU = true;
  const double TIMEOUT = 20;
  const double GLOB_TIMEOUT = 60*60;
  const unsigned Multiplicity=2;
  const double femto_region = 100;
  const unsigned target_yield = 512*1000;//originally 4M
  const unsigned NUM_CPU = 6;
  //TString BaseName = "KstarDist";
  TString BaseName = TString::Format("Eta%.1f",EtaCut);

  //we run to either reproduce the core of 0.97,
  //or the upper limit of reff = 1.06+0.04
  //this leads to a 10% difference in the SP core source
  double rSP_core=0;
  double rSP_dispZ=0;
  double rSP_hadr=0;
  double rSP_hadrZ=0;
  double rSP_hflc=0;
  double rSP_tau=0;
  bool tau_prp = true;
  double rSP_tflc=0;
  double rSP_ThK=0;
  //default is true
  bool rSP_FixedHadr = true;
  //default is 0
  float rSP_FragBeta = 0;
  if(type=="pd"){
    rSP_core = EffFix?0.915*1.100:0.915;
    rSP_dispZ = rSP_core;
    if(EffFix==-1){
      rSP_core = 0.2;
      rSP_dispZ = 0.2;
      //rSP_hadr = 2.7;
      //rSP_tau = 1.7;
      //rSP_hadr = 3.55;
      //rSP_tau = 2.3;
      rSP_hadr = 2.3;
      rSP_hadrZ = 0;
      rSP_tau = 3.55;
      rSP_hflc = 0.00;
    }
    else if(EffFix==-1402){
      //this so happens to be exactly as the usmani fit in the ceca paper
      rSP_core = 0.176;
      rSP_dispZ = 0.176;
      rSP_hadr = 2.68;
      rSP_tau = 3.76;
    }
    else if(EffFix<0){
      rSP_core = 0.25;
      rSP_dispZ = 0.25;
      rSP_hadr = 3.55;
      rSP_tau = 2.3;
    }

  }
  if(type=="pp"||type=="pipi"||type=="pipi_core"||type=="pP"||type=="pK"||type=="Kd"||type=="KdReso"||type=="pL"){
    //jaime issue, to get his stuff, rSP is 0, rSP_hflc = 0.15
    rSP_core = EffFix?0.915*1.100:0.915;
    rSP_core = 1.12;
    rSP_dispZ = rSP_core;
    if((type=="pK"||type=="pipi"||type=="pipi_core")&&EffFix==1){
      rSP_core = 0.8;
      rSP_dispZ = 0.8;
    }
    if(EffFix==-1){
      rSP_core = 0.2;
      rSP_dispZ = 0.2;
      //rSP_hadr = 2.7;
      //rSP_tau = 1.7;
      //rSP_hadr = 3.55;
      //rSP_tau = 2.3;
      rSP_hadr = 2.3;
      rSP_tau = 3.55;
      rSP_hflc = 0.00;
    }
    else if(EffFix==-2){
      rSP_core = 0.25;
      rSP_dispZ = 0.25;
      rSP_hadr = 3.55;
      rSP_tau = 2.3;
      rSP_hflc = 0.1;
    }
    else if(EffFix==-3){
      //this works well
      //rSP_core = 0.35;
      //rSP_hadr = 4.50;
      //rSP_tau = 1.0;
      //rSP_hflc = 0.0;
      //HadronSize = 0.8;
      //HadronSlope = 0.2;
//what works is all zero, but this is 0.25 and the tau is 3.10 exponential
//what kind of works (but not perfect) is to have tau of c.a. 4.5, 15% fluctuations
//the former leads to exp source (also the core), while for the latter it is a Gauss
      rSP_core = 0.0;
      rSP_dispZ = 0.0;
      rSP_hadr = 0.0;
      rSP_tau = 4.5;
      rSP_tflc = 0.15;
      tau_prp = false;
      rSP_hflc = 0.0;
      HadronSize = 0.0;
      HadronSlope = 0.0;
      rSP_ThK = 0.0;
      //small increase of scaling at large mT values
      rSP_FixedHadr = true;
      //can also increase the radius.
      //Big effect on Tau, where it creates strong scaling (at low mT values)
      rSP_FragBeta = 1.0;

      //rSP_core = 0.25;
      //rSP_hadr = 3.55;
      //rSP_tau = 1.0;
    }
    //some random tests
    else if(EffFix==-900){
      rSP_core = 0.32;
      rSP_dispZ = 0.32;
      rSP_hadr = 2.3/2.;
      rSP_hadrZ = 2.3*4.;
      rSP_tau = 3.55*0;
      rSP_hflc = 0.00;
    }
    //for FEMTUM
    else if(EffFix==-1001){
      //rSP_core = 0.80;
      //rSP_dispZ = 0.80;
      //these values are for the ceca paper, work well with 1 fm effective
      rSP_core = 0.85;
      rSP_dispZ = 0.85;
      rSP_hadr = 0.00;
      rSP_tau = 0.00;
    }
    else if(EffFix==-1002){
      rSP_core = 1.20;
      rSP_dispZ = 1.20;
      rSP_hadr = 0.00;
      rSP_tau = 0.00;
    }
    else if(EffFix==-1010){
      rSP_core = 0.00;
      rSP_dispZ = 0.00;
      //rSP_hadr = 4.00;
      //these values are for the ceca paper, work well with 1 fm effective
      rSP_hadr = 4.2;
      rSP_tau = 0.00;
    }
    else if(EffFix==-1020){
      rSP_core = 0.00;
      rSP_dispZ = 0.00;
      rSP_hadr = 6.00;
      rSP_tau = 0.00;
    }
    else if(EffFix==-1100){
      rSP_core = 0.00;
      rSP_dispZ = 0.00;
      rSP_hadr = 0.00;
      //rSP_tau = 4.00;
      //these values are for the ceca paper, work well with 1 fm effective
      rSP_tau = 4.2;
    }
    else if(EffFix==-1200){
      rSP_core = 0.00;
      rSP_dispZ = 0.00;
      rSP_hadr = 0.00;
      rSP_tau = 6.00;
    }
    else if(EffFix==-1300){
      rSP_core = 0.25;
      rSP_dispZ = 0.00;
      rSP_hadr = 7.00;
      rSP_tau = 0.00;
      //rSP_FragBeta = 1.0;
      //tau_prp = false;
      rSP_ThK = 100;
    }
    //examples for the animation for hadron 2023
    else if(EffFix==-1400){
      rSP_core = 0.176;
      rSP_dispZ = 0.176;
      rSP_hadr = 0.00;
      rSP_tau = 0.00;
    }
    else if(EffFix==-1401){
      rSP_core = 0.176;
      rSP_dispZ = 0.176;
      rSP_hadr = 2.68;
      rSP_tau = 0.00;
    }
    else if(EffFix==-1402){
      //this so happens to be exactly as the usmani fit in the ceca paper
      rSP_core = 0.176;
      rSP_dispZ = 0.176;
      rSP_hadr = 2.68;
      rSP_tau = 3.76;
    }
    else if(EffFix==-1403){
      rSP_core = 0.176;
      rSP_dispZ = 0.176;
      rSP_hadr = 2.70;
      rSP_tau = 3.78;
    }
    else if(EffFix==-1404){
      rSP_core = 0.176;
      rSP_dispZ = 0.176;
      rSP_hadr = 2.24;
      rSP_tau = 3.13;
    }
    else if(EffFix==-1405){
      rSP_core = 0.176;
      rSP_dispZ = 0.176;
      rSP_hadr = 1.96;
      rSP_tau = 2.74;
    }
    else if(EffFix==-1406){
      rSP_core = 0.176;
      rSP_dispZ = 0.176;
      rSP_hadr = 1.44;
      rSP_tau = 2.02;
    }
    else if(EffFix==-1407){
      rSP_core = 0.176;
      rSP_dispZ = 0.176;
      rSP_hadr = 2.00;
      rSP_tau = 3.78;
    }
    //these values are for the ceca paper -> with NLO
    else if(EffFix==-1408){
      rSP_core = 0.288;
      rSP_dispZ = 0.288;
      rSP_hadr = 3.23;
      rSP_tau = 3.26;
    }
  }
  //if(type=="Kd") {rSP_core = EffFix?0.950*1.075:0.950;rSP_dispZ = rSP_core;}
  if(type=="Dpi"){
    rSP_core = EffFix?0.9*1.075:0.9;
    rSP_dispZ = rSP_core;
    if(EffFix<0){
      rSP_core = 0.25;
      rSP_dispZ = 0.00;
      rSP_hadr = 3.55;
      rSP_tau = 2.3;
      if(EffFix==-3){
        rSP_core = 1.0;
        rSP_dispZ = 0.0;
        rSP_hadr = 0.0;
        rSP_tau = 4.3*0;
        rSP_tflc = 0.25*0;
        tau_prp = false;
        rSP_hflc = 0.0;
        HadronSize = 0.0;
        HadronSlope = 0.0;
        rSP_ThK = 0.0;
        //small increase of scaling at large mT values
        rSP_FixedHadr = true;
        //can also increase the radius.
        //Big effect on Tau, where it creates strong scaling (at low mT values)
        rSP_FragBeta = 1.0;
      }
    }
  }

  TREPNI Database(0);
  Database.SetSeed(11);
  std::vector<TreParticle*> ParticleList;
  ParticleList.push_back(Database.NewParticle("Proton"));
  ParticleList.push_back(Database.NewParticle("PrimProton"));
  ParticleList.push_back(Database.NewParticle("Lambda"));
  ParticleList.push_back(Database.NewParticle("Kaon"));
  ParticleList.push_back(Database.NewParticle("Deuteron"));
  ParticleList.push_back(Database.NewParticle("Pion"));
  ParticleList.push_back(Database.NewParticle("PionReso"));
  ParticleList.push_back(Database.NewParticle("PionFSI"));
  ParticleList.push_back(Database.NewParticle("Photon"));
  ParticleList.push_back(Database.NewParticle("ProtonReso"));
  ParticleList.push_back(Database.NewParticle("LambdaReso"));
  ParticleList.push_back(Database.NewParticle("KaonReso"));
  ParticleList.push_back(Database.NewParticle("DeuteronReso"));
  ParticleList.push_back(Database.NewParticle("Dch"));
  ParticleList.push_back(Database.NewParticle("Dch_star"));

  TString OutputFolderName;
  if(EffFix) OutputFolderName = "FunWithCeca/Ceca_"+type+"_Frascati_v1";
  else OutputFolderName = "FunWithCeca/Ceca_"+type+"_CoreFix";
OutputFolderName = "dKaon/Frascati/Sources/Raw";
  //OutputFolderName = "FunWithCeca";
  TString BaseFileName = TString::Format("%s/%s/%s_%s_ET%i_PR%i_DD%.1f_EF%i",
  GetCernBoxDimi(),OutputFolderName.Data(),BaseName.Data(),type.Data(),EQUALIZE_TAU,PROTON_RESO,d_delay,EffFix);
  if(type=="Kd"||type=="KdReso"){
   BaseFileName = TString::Format("%s/%s/%s_%s_PR%.2f_KR%.2f_DD%.2f",
    GetCernBoxDimi(),OutputFolderName.Data(),BaseName.Data(),type.Data(),PROTON_RESO*frac_protons,frac_kaons,d_delay);   
  }
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
  if(type=="Kd"||type=="KdReso"){
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

  if(type=="pL"){
    DLM_Histo<float>* dlmTemp = GetPtEta_13TeV(
      TString::Format("%s/CatsFiles/Source/CECA/proton_pT/p_dist_13TeV_ClassI.root",GetCernBoxDimi()),
      "Graph1D_y1", 500, 4050, EtaCut);
    dlm_pT_eta_p = dlmTemp[0];
    delete dlmTemp;
    dlmTemp = GetPtEta_13TeV(
      TString::Format("%s/CatsFiles/Source/CECA/Lambda_pT/L_dist_13TeV_ClassI.root",GetCernBoxDimi()),
      "Graph1D_y1", 400, 8000, EtaCut);
    dlm_pT_eta_d = dlmTemp[0];
    delete dlmTemp;
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
      dlm_pT_eta_p.SampleYield(axisValues);
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
      if(type=="pd"||type=="pp"||type=="pK"||type=="pL") prt->SetAbundance(frac_protons+(100.-frac_protons)*(!PROTON_RESO));
      else if(type=="pP") prt->SetAbundance(100.);
      else prt->SetAbundance(0);
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);
      if(dlm_pT_p || type=="pL") prt->SetPtEtaPhi(dlm_pT_eta_p);
      else
        prt->SetPtPz(0.85*prt->GetMass(),0.85*prt->GetMass());
    }
    else if(prt->GetName()=="Lambda"){
      prt->SetMass(Mass_L);
      if(type=="pL") prt->SetAbundance((100.-64.38));
      else prt->SetAbundance(0);
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);
      if(dlm_pT_d || type=="pL") prt->SetPtEtaPhi(dlm_pT_eta_d);
      else
        prt->SetPtPz(0.85*prt->GetMass(),0.85*prt->GetMass());
    }    
    else if(prt->GetName()=="Kaon"){
      prt->SetMass(Mass_Kch);
      if(type=="Kd"||type=="KdReso"||type=="pK") prt->SetAbundance(frac_kaons+(100.-frac_kaons)*(!PROTON_RESO));
    //if(type=="Kd"||type=="pK") prt->SetAbundance(100);
      else prt->SetAbundance(0);
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);

      if(type=="Kd"||type=="KdReso"){
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
      else if(type=="KdReso") prt->SetAbundance(100.*pow(frac_protons*0.01,2.));
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
      if(type=="pd"||type=="pp"||type=="pP"||type=="pK"||type=="pL") prt->SetAbundance((100.-frac_protons)*PROTON_RESO);
      else prt->SetAbundance(0);
      prt->SetWidth(hbarc/1.65);
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);
      if(dlm_pT_p || type=="pL")prt->SetPtEtaPhi(dlm_pT_eta_p);
      else
        prt->SetPtPz(0.85*prt->GetMass(),0.85*prt->GetMass());


      prt->NewDecay();
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Proton"));
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
      prt->GetDecay(0)->SetBranching(100);
    }
    else if(prt->GetName()=="DeuteronReso"){
      //so we assume this effective resonances has the same mass difference
      //to the detueron as the mass difference between protons and their resonances
      prt->SetMass(Mass_d + 1362. - Mass_p);
      if(type=="KdReso") prt->SetAbundance(100.-100.*pow(frac_protons*0.01,2.));
      else prt->SetAbundance(0);
      prt->SetWidth(hbarc/1.65);
      prt->SetRadius(HadronSize*1.5);
      prt->SetRadiusSlope(HadronSlope);
      //the resonance is also delayed. While this makes no sence, it is kind of equivelent as delaying the deuteron formation.
      prt->SetDelayTau(d_delay);
      if(dlm_pT_d)prt->SetPtEtaPhi(dlm_pT_eta_d);
      else prt->SetPtPz(0.85*prt->GetMass(),0.85*prt->GetMass());

      prt->NewDecay();
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Deuteron"));
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
      prt->GetDecay(0)->SetBranching(100);
    }
    else if(prt->GetName()=="LambdaReso"){
      prt->SetMass(1462);
      if(type=="pL") prt->SetAbundance(64.38);
      else prt->SetAbundance(0);
      prt->SetWidth(hbarc/4.69);
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);
      if(dlm_pT_d || type=="pL")prt->SetPtEtaPhi(dlm_pT_eta_d);
      else
        prt->SetPtPz(0.85*prt->GetMass(),0.85*prt->GetMass());


      prt->NewDecay();
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Lambda"));
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
      prt->GetDecay(0)->SetBranching(100);
    }
    else if(prt->GetName()=="KaonReso"){
      prt->SetMass(1054);
      if(type=="Kd"||type=="KdReso"||type=="pK") prt->SetAbundance((100.-frac_kaons)*PROTON_RESO);
    //if(type=="Kd"||type=="pK") prt->SetAbundance(0);
      else prt->SetAbundance(0);
      prt->SetWidth(hbarc/3.66);
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);

      if(type=="Kd"||type=="KdReso"){
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

/*
prt->SetPtPz(0.85*1000,0.85*1000);
double TMPMASS;
if(prt->GetName()=="Proton"||prt->GetName()=="PionFSI"||prt->GetName()=="Dch"){
  TMPMASS = 1000;
  prt->SetMassLimit(TMPMASS,TMPMASS);
  prt->SetMass(TMPMASS);
  prt->SetAbundance(100);
}
if(prt->GetName()=="ProtonReso"||prt->GetName()=="PionReso"){
  TMPMASS = 1400;
  prt->SetMassLimit(TMPMASS,TMPMASS);
  prt->SetMass(TMPMASS);
  prt->SetAbundance(300);
  prt->SetWidthLimit(hbarc/1.50,hbarc/1.50);
  prt->SetWidth(hbarc/1.50);
}


if(prt->GetName()=="Dch"){
  //TMPMASS = 1500;
  //prt->SetMassLimit(TMPMASS,TMPMASS);
  //prt->SetMass(TMPMASS);
  //prt->SetPtPz(1.7*1000,1.7*1000);
  //prt->SetAcceptance_pT(0,1000);
  //prt->SetAcceptance_pT(1000,1e6);
}
if(prt->GetName()=="PionFSI"){
  TMPMASS = 100;
  prt->SetMassLimit(TMPMASS,TMPMASS);
  prt->SetMass(TMPMASS);
  //prt->SetPtPz(0.25*1000,0.25*1000);
  //prt->SetAcceptance_pT(0,1000);
  //prt->SetAcceptance_pT(1000,1e6);
}
*/

  }

  std::vector<std::string> ListOfParticles;
  if(type=="pd"||type=="pK") ListOfParticles.push_back("Proton");
  if(type=="Kd"||type=="pK"||type=="KdReso") ListOfParticles.push_back("Kaon");
  if(type=="pd"||type=="Kd"||type=="KdReso") ListOfParticles.push_back("Deuteron");

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
  if(type=="pL"){
    ListOfParticles.push_back("Proton");
    ListOfParticles.push_back("Lambda");
  }


//printf("%s %s\n", ListOfParticles.at(0).c_str(),ListOfParticles.at(1).c_str());
//return;
  CECA Ivana(Database,ListOfParticles);

  //Ivana.SetDisplacement(rSP_core);
  //Ivana.SetHadronization(0.0);
  //Ivana.SetHadrFluctuation(0.0);

printf("d %.3f %.3f\n",rSP_core,rSP_dispZ);
printf("h %.3f %.3f\n",rSP_hadr,rSP_hadrZ);
printf("t %.3f %.3f\n",rSP_tau,float(tau_prp));
printf("tf %.3f\n",rSP_tflc);
printf("tk %.3f\n",rSP_ThK);
printf("fh %.3f\n",float(rSP_FixedHadr));
printf("fb %.3f\n",rSP_FragBeta);

  Ivana.SetDisplacementZ(rSP_dispZ);
  Ivana.SetDisplacementT(rSP_core);
  Ivana.SetHadronizationZ(rSP_hadrZ);//0
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
/*
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
*/


    //ceca paper
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

    Ivana.Ghetto_NumRadBins = 2048;
    Ivana.Ghetto_RadMin = 0;
    Ivana.Ghetto_RadMax = 64;

  }
  if(type=="Kd"||type=="KdReso"){

    //ceca paper
    Ivana.Ghetto_NumMtBins = 5;
    Ivana.Ghetto_MtBins = new double [Ivana.Ghetto_NumMtBins+1];
    Ivana.Ghetto_MtBins[0] = 1180; //avg  983 ( 985)
    Ivana.Ghetto_MtBins[1] = 1300;//avg 1054 (1055)
    Ivana.Ghetto_MtBins[2] = 1420;//avg 1110 (1110)
    Ivana.Ghetto_MtBins[3] = 1600;//avg 1168 (1170)
    Ivana.Ghetto_MtBins[4] = 1800;//avg 1228 (1230)
    Ivana.Ghetto_MtBins[5] = 2400;//avg 1315 (1315)

    Ivana.Ghetto_NumMomBins = 20;
    Ivana.Ghetto_MomMin = 0;
    Ivana.Ghetto_MomMax = 400;

    Ivana.Ghetto_NumRadBins = 2048;
    Ivana.Ghetto_RadMin = 0;
    Ivana.Ghetto_RadMax = 96;

  }

  if(type=="pL"){

    Ivana.Ghetto_NumMtBins = 16;
    Ivana.Ghetto_MtBins = new double [Ivana.Ghetto_NumMtBins+1];
    Ivana.Ghetto_MtBins[0] = 1000;
    Ivana.Ghetto_MtBins[1] = 1100;
    Ivana.Ghetto_MtBins[2] = 1200;
    Ivana.Ghetto_MtBins[3] = 1300;
    Ivana.Ghetto_MtBins[4] = 1350;
    Ivana.Ghetto_MtBins[5] = 1400;
    Ivana.Ghetto_MtBins[6] = 1500;
    Ivana.Ghetto_MtBins[7] = 1600;
    Ivana.Ghetto_MtBins[8] = 1700;
    Ivana.Ghetto_MtBins[9] = 1800;
    Ivana.Ghetto_MtBins[10] = 1900;
    Ivana.Ghetto_MtBins[11] = 2000;
    Ivana.Ghetto_MtBins[12] = 2100;
    Ivana.Ghetto_MtBins[13] = 2200;
    Ivana.Ghetto_MtBins[14] = 2300;
    Ivana.Ghetto_MtBins[15] = 2400;
    Ivana.Ghetto_MtBins[16] = 2500;

    Ivana.Ghetto_NumMomBins = 150;
    Ivana.Ghetto_MomMin = 0;
    Ivana.Ghetto_MomMax = 600;

    Ivana.Ghetto_NumRadBins = 2048;
    Ivana.Ghetto_RadMin = 0;
    Ivana.Ghetto_RadMax = 64;

  }
  else if(type=="pd"){
    Ivana.Ghetto_NumMtBins = 64;
    Ivana.Ghetto_MtMin  = 1000;
    Ivana.Ghetto_MtMax  = 4200;

    //Ivana.Ghetto_NumMtBins = 10;
    //Ivana.Ghetto_MtBins = new double [Ivana.Ghetto_NumMtBins+1];
    //Ivana.Ghetto_MtBins[0] = 930; //avg  983 ( 985)
    //Ivana.Ghetto_MtBins[1] = 1020;//avg 1054 (1055)
    //Ivana.Ghetto_MtBins[2] = 1080;//avg 1110 (1110)
    //Ivana.Ghetto_MtBins[3] = 1140;//avg 1168 (1170)
    //Ivana.Ghetto_MtBins[4] = 1200;//avg 1228 (1230)
    //Ivana.Ghetto_MtBins[5] = 1260;//avg 1315 (1315)
    //Ivana.Ghetto_MtBins[6] = 1380;//avg 1463 (1460)
    //Ivana.Ghetto_MtBins[7] = 1570;//avg 1681 (1680)
    //Ivana.Ghetto_MtBins[8] = 1840;//avg 1923 (1920)
    //Ivana.Ghetto_MtBins[9] = 2030;//avg 2303 (2300)
    //Ivana.Ghetto_MtBins[10] = 4500;


  }

  Ivana.GoBabyGo(NUM_CPU);

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

  Ivana.GhettoFemto_mT_mTwrong->ComputeError();
  TH2F* h_GhettoFemto_mT_mTwrong = Convert_DlmHisto_TH2F(Ivana.GhettoFemto_mT_mTwrong,"GhettoFemto_mT_mTwrong");

  Ivana.Ghetto_mT_mTwrong->ComputeError();
  TH2F* h_Ghetto_mT_mTwrong = Convert_DlmHisto_TH2F(Ivana.Ghetto_mT_mTwrong,"Ghetto_mT_mTwrong");

  Ivana.GhettoFemto_pT1_pT2->ComputeError();
  TH2F* h_GhettoFemto_pT1_pT2 = Convert_DlmHisto_TH2F(Ivana.GhettoFemto_pT1_pT2,"GhettoFemto_pT1_pT2");

  Ivana.GhettoFemto_pT1_div_pT->ComputeError();
  TH1F* h_GhettoFemto_pT1_div_pT = Convert_DlmHisto_TH1F(Ivana.GhettoFemto_pT1_div_pT,"GhettoFemto_pT1_div_pT");

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
    if(Mean&&Err&&hProj->GetEntries()>128){
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
  TF1* fitDG_rstar = new TF1("fitDG_rstar",NormDoubleGaussSourceTF1,0.4,16.0,4);
  //fitDG_rstar->SetParameter(3,1);
  //fitDG_rstar->SetParLimits(3,0.8,1.0);//norm
  fitDG_rstar->FixParameter(3,1);

  //fitDG_rstar->SetParameter(0,0.36);//3.30932e-01
  //fitDG_rstar->SetParLimits(0,0.0,0.6);//first source
  fitDG_rstar->FixParameter(0,fit_rcore->GetParameter(1));//3.30932e-01

  fitDG_rstar->SetParameter(1,1.30);//1.23015e+00
  fitDG_rstar->SetParLimits(1,0.9,5.0);//second source

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

  h_Ghetto_mT_mTwrong->Write();
  h_GhettoFemto_mT_mTwrong->Write();
  h_GhettoFemto_pT1_pT2->Write();
  h_GhettoFemto_pT1_div_pT->Write();

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
  delete h_Ghetto_mT_mTwrong;
  delete h_GhettoFemto_mT_mTwrong;
  delete h_GhettoFemto_pT1_pT2;
  delete h_GhettoFemto_pT1_div_pT;
  delete fit_rstar;
  delete fit_rcore;
  delete fitDG_rstar;
  delete cSource;
  delete cMt;
}


//the ModeAB: XCBA
//A = Mass; 1 = 1 GeV, 2 = 0.1 GeV
//B = MomDist; 1 = 0.85*Mass GeV (full), 2 : <1*Mass GeV (low), 3 : >1*Mass GeV (high)
//C = 1 no reso, 2: with reso
//the Mode CBA:
//CBA: 1 = random source
//CBA: 2 = hadr only
//CBA: 3 = tau only, light fragments
void Ceca_AB(const int& ModeA, const int& ModeB, const int Mode){
  printf("\n-------------------------------------\n");
  printf("%i %i %i\n",ModeA,ModeB,Mode);
  const double TIMEOUT = 40;

  double DispT=0;
  double DispZ=0;

  double HadrT=0;
  double HadrZ=0;
  double HadrF=0;

  double Tau=0;
  double TauF=0;
  double FragBeta=0;

  double MassA;
  double MassB;
  double MassDA;
  double MassDB;
  double MassRA;
  double MassRB;
  double FracRA;
  double FracRB;

  double EtaLowA;
  double EtaUpA;
  double EtaLowB;
  double EtaUpB;

  double TauRA = 2.0;
  double TauRB = 2.0;

  float PtLowA;
  float PtUpA;
  float PtLowB;
  float PtUpB;

  float PtDistA = 0.85;
  float PtDistB = 0.85;

  TREPNI Database(0);
  Database.SetSeed(11);
  std::vector<TreParticle*> ParticleList;
  ParticleList.push_back(Database.NewParticle("A"));
  ParticleList.push_back(Database.NewParticle("B"));
  ParticleList.push_back(Database.NewParticle("RA"));
  ParticleList.push_back(Database.NewParticle("RB"));
  ParticleList.push_back(Database.NewParticle("DA"));
  ParticleList.push_back(Database.NewParticle("DB"));

  TString OutputFolderName = "FunWithCeca/Ceca_AB";
  TFile fOutput(TString::Format("%s/%s/Ceca_AB_MA%i_MB%i_M%i.root",
  GetFemtoOutputFolder(),OutputFolderName.Data(),ModeA,ModeB,Mode),"recreate");

//N.B. this is 45 deg, since we only apply acceptance to the primordials, to have
//50:50 prim vs second we need to make the ration more or less 1:2 (180 of 360 deg covarage)
  EtaLowA = -0.88;
  EtaUpA = 0.88;
  EtaLowB = -0.88;
  EtaUpB = 0.88;
  switch (ModeA%10) {
    case 1: MassA = 1000; MassDA = 100; MassRA = 1500; break;
    case 2: MassA = 100; MassDA = 100; MassRA = 1500; break;
    default : printf("MassA==???\n"); return;
  }
  switch ((ModeA/10)%10) {
    case 1: PtLowA = 0; PtUpA = 1e6; break;
    case 2: PtLowA = 0; PtUpA = 1*MassA; break;
    case 3: PtLowA = 1*MassA; PtUpA = 1e6; break;
    default : printf("PtA==???\n"); return;
  }
  switch ((ModeA/100)%10) {
    case 1: FracRA = 0; break;
    case 2: FracRA = 66; break;
    default : printf("FracRA==???\n"); return;
  }

  switch (ModeB%10) {
    case 1: MassB = 1000; MassDB = 100; MassRB = 1500; break;
    case 2: MassB = 100; MassDB = 100; MassRB = 1500; break;
    default : printf("MassB==???\n"); return;
  }
  switch ((ModeB/10)%10) {
    case 1: PtLowB = 0; PtUpB = 1e6; break;
    case 2: PtLowB = 0; PtUpB = 1*MassB; break;
    case 3: PtLowB = 1*MassB; PtUpB = 1e6; break;
    default : printf("PtB==???\n"); return;
  }
  switch ((ModeB/100)%10) {
    case 1: FracRB = 0; break;
    case 2: FracRB = 66; break;
    default : printf("FracRB==???\n"); return;
  }


  switch (Mode) {
    case 1:
      DispT = 1.0;
      DispZ = 1.0;
      break;
    case 2:
      HadrT = 4.0;
      HadrZ = 4.0;
      break;
    case 3:
      Tau = 4.0;
      TauF = 0.25;
      FragBeta = 1.0;
      break;
    default : printf("MassB==???\n"); return;
  }


  for(TreParticle* prt : ParticleList){
    if(prt->GetName()=="A"){
      prt->SetMass(MassA);
      prt->SetAbundance(100.-FracRA);
      prt->SetPtPz(PtDistA*prt->GetMass(),PtDistA*prt->GetMass());
      prt->SetAcceptance_pT(PtLowA,PtUpA);
      prt->SetAcceptance_Eta(EtaLowA,EtaUpA);
    }
    else if(prt->GetName()=="B"){
      prt->SetMass(MassB);
      prt->SetAbundance(100.-FracRB);
      prt->SetPtPz(PtDistB*prt->GetMass(),PtDistB*prt->GetMass());
      prt->SetAcceptance_pT(PtLowB,PtUpB);
      prt->SetAcceptance_Eta(EtaLowB,EtaUpB);
    }
    //!!! N.B. The daughters and resonances WILL NOT have an ACCEPTANCE cut.
    //i.e. we assume they can be produced anywhere, it is just A and B that are
    //subject to acceptence, i.e. if they can be measured or not
    else if(prt->GetName()=="DA"){
      prt->SetMass(MassDA);
      prt->SetAbundance(0);
      prt->SetPtPz(PtDistA*prt->GetMass(),PtDistA*prt->GetMass());
    }
    else if(prt->GetName()=="DB"){
      prt->SetMass(MassDB);
      prt->SetAbundance(0);
      prt->SetPtPz(PtDistB*prt->GetMass(),PtDistB*prt->GetMass());
    }
    else if(prt->GetName()=="RA"){
      prt->SetMass(MassRA);
      prt->SetAbundance(FracRA);
      prt->SetPtPz(PtDistA*prt->GetMass(),PtDistA*prt->GetMass());

      prt->NewDecay();
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("A"));
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("DA"));
      prt->GetDecay(0)->SetBranching(100);
      prt->SetWidth(hbarc/TauRA);
    }
    else if(prt->GetName()=="RB"){
      prt->SetMass(MassRB);
      prt->SetAbundance(FracRB);
      prt->SetPtPz(PtDistB*prt->GetMass(),PtDistB*prt->GetMass());

      prt->NewDecay();
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("B"));
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("DB"));
      prt->GetDecay(0)->SetBranching(100);
      prt->SetWidth(hbarc/TauRB);
    }
    else{
      printf("WTF!?!?!?!");
      return;
    }
  }//ParticleList

  std::vector<std::string> ListOfParticles;
  ListOfParticles.push_back("A");
  ListOfParticles.push_back("B");

  CECA Ivana(Database,ListOfParticles);

  Ivana.SetDisplacementZ(DispZ);
  Ivana.SetDisplacementT(DispT);
  Ivana.SetHadronizationZ(HadrZ);
  Ivana.SetHadronizationT(HadrT);
  Ivana.SetHadrFluctuation(HadrF);
  Ivana.SetTau(Tau,false);
  Ivana.SetTauFluct(TauF);
  Ivana.SetThermalKick(0);
  Ivana.SetFixedHadr(true);
  Ivana.SetFragmentBeta(FragBeta);

  Ivana.SetTargetStatistics(10);
  Ivana.SetEventMult(2);
  Ivana.SetSourceDim(2);
  Ivana.SetDebugMode(true);
  Ivana.SetThreadTimeout(TIMEOUT);
  Ivana.EqualizeFsiTime(true);
  Ivana.SetFemtoRegion(100);
  Ivana.GHETTO_EVENT = true;

  Ivana.GoBabyGo();

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
  printf("reff_Ceca = %f\n",reff_Ceca);
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
  h_Ghetto_kstar->Write();
  h_Ghetto_kstar_rstar->Write();
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
  delete cSource;
  delete cMt;
}




//code for Dmeson pion and Dmeson proton
//tdel: delay of Dmesons,LightMesons,LightMesonReso
void Ceca_Dmeson_LightMeson_1(const TString type, std::vector<float>& disp, std::vector<float>& hadr, const float& hfluct,
  std::vector<float>& tdel, const float& TimeDel, const float TIMEOUT=60){

  if(type!="Dpi"&&type!="DK"&&type!="pp"){
    printf("WHAT IS THIS: %s\n", type.Data());
    return;
  }
  TREPNI Database(0);
  Database.SetSeed(11);
  std::vector<TreParticle*> ParticleList;
  ParticleList.push_back(Database.NewParticle("Dmeson"));
  ParticleList.push_back(Database.NewParticle("Kaon"));
  ParticleList.push_back(Database.NewParticle("Pion"));
  ParticleList.push_back(Database.NewParticle("Daughter"));//decay daughters (pions)
  ParticleList.push_back(Database.NewParticle("KaonReso"));
  ParticleList.push_back(Database.NewParticle("PionReso"));

  TString OutputFolderName;
  OutputFolderName = "FunWithCeca/Ceca_Dmeson_LightMeson_1";
  TFile fOutput(TString::Format("%s/%s/Dmesons_%s.root",
  GetFemtoOutputFolder(),OutputFolderName.Data(),type.Data()),"recreate");

  DLM_Histo<float>* dlm_pT = NULL;
  DLM_Histo<float>* dlm_eta = NULL;
  DLM_Histo<float> dlm_pT_eta_D;
  DLM_Histo<float> dlm_pT_eta_LM;

  TH1F* h_pT = NULL;
  TH1F* h_eta = NULL;
  double* BinRange;

  TString FilePath;
  if(type=="Dpi"||type=="DK"){
    FilePath = TString::Format("%s/CatsFiles/Source/CECA/Dmeson_LightMeson/",GetCernBoxDimi());
    TFile* RootFile = NULL;
    double Norm;
    //for the Dmesons//////////////////////////////////////////////////
    RootFile = new TFile(FilePath+"pT_Dmeson.root","read");
    h_pT = (TH1F*)RootFile->Get("slice_px_of_fHistDplusInvMassPtSel");
    Norm = h_pT->Integral()*h_pT->GetNbinsX();
    dlm_pT = Convert_TH1F_DlmHisto(h_pT);
    *dlm_pT /= Norm;
    dlm_pT->RescaleAxis(0,1000,false);
    delete RootFile; RootFile=NULL;

    RootFile = new TFile(FilePath+"eta_Dmeson.root","read");
    h_eta = (TH1F*)RootFile->Get("fHistDplusEta");
    Norm = h_eta->Integral()*h_eta->GetNbinsX();
    dlm_eta = Convert_TH1F_DlmHisto(h_pT);
    *dlm_eta /= Norm;
    delete RootFile; RootFile=NULL;

    dlm_pT_eta_D.SetUp(2);
    BinRange = dlm_pT->GetBinRange(0);
    dlm_pT_eta_D.SetUp(0,dlm_pT->GetNbins(),BinRange);
    delete[]BinRange;
    BinRange = dlm_eta->GetBinRange(0);
    dlm_pT_eta_D.SetUp(1,dlm_eta->GetNbins(),BinRange);
    delete[]BinRange;
    dlm_pT_eta_D.Initialize();

    for(unsigned upT=0; upT<dlm_pT->GetNbins(); upT++){
      for(unsigned ueta=0; ueta<dlm_eta->GetNbins(); ueta++){
        dlm_pT_eta_D.SetBinContent(upT,ueta,dlm_pT->GetBinContent(upT)*dlm_eta->GetBinContent(ueta));
      }
    }

    dlm_pT->ScaleToIntegral();
    dlm_pT->ScaleToBinSize();
    fOutput.cd();
    TH1F* H_pT = Convert_DlmHisto_TH1F(dlm_pT,"dlm_pT_Dmeson");
    fOutput.cd();
    H_pT->Write();
    delete H_pT;

    delete dlm_pT;dlm_pT=NULL;
    delete dlm_eta;dlm_eta=NULL;
    //////////////////////////////////////////////////////////////////

    //for the pions/Kaons//////////////////////////////////////////////////
    RootFile = new TFile(FilePath+TString(type=="Dpi"?"pT_pions.root":"pT_Kaons.root"),"read");
    h_pT = (TH1F*)RootFile->Get("pTDist_after");
    Norm = h_pT->Integral()*h_pT->GetNbinsX();
    dlm_pT = Convert_TH1F_DlmHisto(h_pT);
    *dlm_pT /= Norm;
    dlm_pT->RescaleAxis(0,1000,false);
    delete RootFile; RootFile=NULL;

    RootFile = new TFile(FilePath+TString(type=="Dpi"?"eta_pions.root":"eta_Kaons.root"),"read");
    h_eta = (TH1F*)RootFile->Get("EtaDist_after");
    Norm = h_eta->Integral()*h_eta->GetNbinsX();
    dlm_eta = Convert_TH1F_DlmHisto(h_pT);
    *dlm_eta /= Norm;
    delete RootFile; RootFile=NULL;

    dlm_pT_eta_LM.SetUp(2);
    BinRange = dlm_pT->GetBinRange(0);
    dlm_pT_eta_LM.SetUp(0,dlm_pT->GetNbins(),BinRange);
    delete[]BinRange;
    BinRange = dlm_eta->GetBinRange(0);
    dlm_pT_eta_LM.SetUp(1,dlm_eta->GetNbins(),BinRange);
    delete[]BinRange;
    dlm_pT_eta_LM.Initialize();

    for(unsigned upT=0; upT<dlm_pT->GetNbins(); upT++){
      for(unsigned ueta=0; ueta<dlm_eta->GetNbins(); ueta++){
        dlm_pT_eta_LM.SetBinContent(upT,ueta,dlm_pT->GetBinContent(upT)*dlm_eta->GetBinContent(ueta));
      }
    }

    dlm_pT->ScaleToIntegral();
    dlm_pT->ScaleToBinSize();
    fOutput.cd();
    H_pT = Convert_DlmHisto_TH1F(dlm_pT,"dlm_pT_LightMeson");
    fOutput.cd();
    H_pT->Write();
    delete H_pT;

    delete dlm_pT;dlm_pT=NULL;
    delete dlm_eta;dlm_eta=NULL;
    //////////////////////////////////////////////////////////////////

  }


  fOutput.cd();
  TH2F* h_pT_eta_D = Convert_DlmHisto_TH2F(&dlm_pT_eta_D,"h_pT_eta_D");
  TH2F* h_pT_eta_LM = Convert_DlmHisto_TH2F(&dlm_pT_eta_LM,"h_pT_eta_LM");
  fOutput.cd();
  h_pT_eta_D->Write();
  h_pT_eta_LM->Write();

//  ParticleList.push_back(Database.NewParticle("Dmeson"));
//  ParticleList.push_back(Database.NewParticle("Kaon"));
//  ParticleList.push_back(Database.NewParticle("Pion"));
//  ParticleList.push_back(Database.NewParticle("KaonReso"));
//  ParticleList.push_back(Database.NewParticle("Daughter"));//decay daughters (pions)
//  ParticleList.push_back(Database.NewParticle("PionReso"));
//  ParticleList.push_back(Database.NewParticle("PionDaught"));
  for(TreParticle* prt : ParticleList){
    if(prt->GetName()=="Dmeson"){
      prt->SetMass(Mass_Dch);
      prt->SetAbundance(100);
      prt->SetPtEtaPhi(dlm_pT_eta_D);
      prt->SetDelayTau(tdel.at(0));
    }
    else if(prt->GetName()=="Daughter"){
      prt->SetMass(Mass_pic);
      prt->SetAbundance(0);
    }
    else if(prt->GetName()=="Kaon"){
      prt->SetMass(Mass_Kch);
      if(type=="DK") prt->SetAbundance(53.4);
      else prt->SetAbundance(0);
      prt->SetPtEtaPhi(dlm_pT_eta_LM);
      prt->SetDelayTau(tdel.at(1));
    }
    else if(prt->GetName()=="KaonReso"){
      prt->SetMass(1054);
      if(type=="DK") prt->SetAbundance(47.6);
      else prt->SetAbundance(0);
      prt->SetWidth(hbarc/3.66);
      prt->SetPtEtaPhi(dlm_pT_eta_LM);
      prt->SetDelayTau(tdel.at(2));

      prt->NewDecay();
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Kaon"));
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Daughter"));
      prt->GetDecay(0)->SetBranching(100);
    }
    else if(prt->GetName()=="Pion"){
      prt->SetMass(Mass_pic);
      if(type=="Dpi") prt->SetAbundance(31.8);
      else prt->SetAbundance(0);
      prt->SetPtEtaPhi(dlm_pT_eta_LM);
      prt->SetDelayTau(tdel.at(1));
    }
    else if(prt->GetName()=="PionReso"){
      prt->SetMass(1180);
      if(type=="Dpi") prt->SetAbundance(68.2);
      else prt->SetAbundance(0);
      prt->SetWidth(hbarc/1.50);
      prt->SetPtEtaPhi(dlm_pT_eta_LM);
      prt->SetDelayTau(tdel.at(2));

      prt->NewDecay();
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Daughter"));
      prt->GetDecay(0)->SetBranching(100);
    }
    //prt->SetRadius(1.8);
    //prt->SetRadiusSlope(0.1);

  }

  std::vector<std::string> ListOfParticles;
  if(type=="Dpi"){
    ListOfParticles.push_back("Dmeson");
    ListOfParticles.push_back("Pion");
  }
  if(type=="DK"){
    ListOfParticles.push_back("Dmeson");
    ListOfParticles.push_back("Kaon");
  }

  printf("%s %s\n", ListOfParticles.at(0).c_str(),ListOfParticles.at(1).c_str());

  CECA Ivana(Database,ListOfParticles);
  Ivana.SetDisplacementX(disp.at(0));
  Ivana.SetDisplacementY(disp.at(1));
  Ivana.SetDisplacementZ(disp.at(2));

  Ivana.SetHadronizationX(hadr.at(0));
  Ivana.SetHadronizationY(hadr.at(1));
  Ivana.SetHadronizationZ(hadr.at(2));

  Ivana.SetHadrFluctuation(hfluct);

  Ivana.SetTau(TimeDel);

  Ivana.SetTargetStatistics(10);
  Ivana.SetEventMult(2);
  Ivana.SetSourceDim(2);
  Ivana.SetDebugMode(true);
  Ivana.SetThreadTimeout(TIMEOUT);
  Ivana.GHETTO_EVENT = true;
  Ivana.GoBabyGo();


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


  Ivana.GhettoSP_pT_1->ComputeError();
  Ivana.GhettoSP_pT_1->ScaleToIntegral();
  Ivana.GhettoSP_pT_1->ScaleToBinSize();
  TH1F* h_GhettoSP_pT_1 = Convert_DlmHisto_TH1F(Ivana.GhettoSP_pT_1,"GhettoSP_pT_1");

  Ivana.GhettoSP_pT_2->ComputeError();
  Ivana.GhettoSP_pT_2->ScaleToIntegral();
  Ivana.GhettoSP_pT_2->ScaleToBinSize();
  TH1F* h_GhettoSP_pT_2 = Convert_DlmHisto_TH1F(Ivana.GhettoSP_pT_2,"GhettoSP_pT_2");

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

  Ivana.GhettoFemto_mT_rstar->ComputeError();
  TH2F* h_GhettoFemto_mT_rstar = Convert_DlmHisto_TH2F(Ivana.GhettoFemto_mT_rstar,"GhettoFemto_mT_rstar");

  TGraph gfm_mT;
  gfm_mT.SetName("gfm_mT");
  gfm_mT.SetLineWidth(3);
  TGraph* gfm_kstar = new TGraph [Ivana.GhettoFemto_mT_rstar->GetNbins(0)];
  for(unsigned umT=0; umT<Ivana.GhettoFemto_mT_rstar->GetNbins(0); umT++){
    TH1D* hProj = h_GhettoFemto_mT_rstar->ProjectionY(TString::Format("hProj"),umT+1,umT+1);
    double Mean = hProj->GetMean();
    delete hProj;
    double GFM = GaussFromMean(Mean);
    gfm_mT.SetPoint(umT,Ivana.GhettoFemto_mT_rstar->GetBinCenter(0,umT),GFM);

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
      GFM = GaussFromMean(Mean);
      gfm_kstar[umT].SetPoint(ukstar,Ivana.Ghetto_kstar_rstar_mT->GetBinCenter(0,ukstar),GFM);
    }

    fOutput.cd();
    h_Ghetto_kstar_rstar_mT->Write();
    //delete h_Ghetto_kstar_rstar_mT;
  }


  fOutput.cd();
  h_GhettoFemto_mT_rstar->Write();
  gfm_mT.Write();
  for(unsigned umT=0; umT<Ivana.GhettoFemto_mT_rstar->GetNbins(0); umT++){
    gfm_kstar[umT].Write();
  }
  h_GhettoSP_pT_1->Write();
  h_GhettoSP_pT_2->Write();
  h_GhettoSPr_X->Write();
  h_GhettoSPr_Y->Write();
  h_GhettoSPr_Z->Write();
  h_GhettoSPr_Rho->Write();
  h_GhettoSPr_R->Write();

  delete h_GhettoFemto_mT_rstar;
  delete h_GhettoSP_pT_1;
  delete h_GhettoSP_pT_2;
  delete h_GhettoSPr_X;
  delete h_GhettoSPr_Y;
  delete h_GhettoSPr_Z;
  delete h_GhettoSPr_Rho;
  delete h_GhettoSPr_R;
/*
  TH1F* h_Ghetto_rstar = Convert_DlmHisto_TH1F(Ivana.Ghetto_rstar,"Ghetto_rstar");
  TH1F* h_Ghetto_kstar = Convert_DlmHisto_TH1F(Ivana.Ghetto_kstar,"Ghetto_kstar");
  Ivana.Ghetto_kstar_rstar->ComputeError();
  TH2F* h_Ghetto_kstar_rstar = Convert_DlmHisto_TH2F(Ivana.Ghetto_kstar_rstar,"Ghetto_kstar_rstar");
  TH2F* h_Ghetto_mT_rstar = Convert_DlmHisto_TH2F(Ivana.Ghetto_mT_rstar,"Ghetto_mT_rstar");
  TH2F* h_GhettoFemto_mT_rstar = Convert_DlmHisto_TH2F(Ivana.GhettoFemto_mT_rstar,"GhettoFemto_mT_rstar");
  TH2F* h_GhettoFemto_mT_kstar = Convert_DlmHisto_TH2F(Ivana.GhettoFemto_mT_kstar,"GhettoFemto_mT_kstar");
  TH2F* h_Ghetto_mT_costh = Convert_DlmHisto_TH2F(Ivana.Ghetto_mT_costh,"Ghetto_mT_costh");
  TH2F* h_GhettoSP_pT_th = Convert_DlmHisto_TH2F(Ivana.GhettoSP_pT_th,"GhettoSP_pT_th");
  TH1F* h_GhettoSP_pT_1 = Convert_DlmHisto_TH1F(Ivana.GhettoSP_pT_1,"GhettoSP_pT_1");
  TH1F* h_GhettoSP_pT_2 = Convert_DlmHisto_TH1F(Ivana.GhettoSP_pT_2,"GhettoSP_pT_2");
  Ivana.GhettoFemto_rstar->ComputeError();
  TH1F* h_GhettoFemto_rstar = Convert_DlmHisto_TH1F(Ivana.GhettoFemto_rstar,"GhettoFemto_rstar");
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

  TGraphErrors g_GhettoFemto_mT_rstar;
  g_GhettoFemto_mT_rstar.SetName("g_GhettoFemto_mT_rstar");
  g_GhettoFemto_mT_rstar.SetMarkerStyle(20);
  g_GhettoFemto_mT_rstar.SetMarkerSize(1);
  g_GhettoFemto_mT_rstar.SetLineWidth(3);

  unsigned uPoint = 0;
  for(unsigned uBin=0; uBin<h_GhettoFemto_mT_rstar->GetXaxis()->GetNbins(); uBin++){
    TH1D* hProj = h_GhettoFemto_mT_rstar->ProjectionY(TString::Format("hProj"),uBin+1,uBin+1);
    double Mean = hProj->GetMean();
    double Err = hProj->GetStdDev();
    double mT = h_GhettoFemto_mT_rstar->GetXaxis()->GetBinCenter(uBin+1);
    if(Mean&&Err){
      g_GhettoFemto_mT_rstar.SetPoint(uPoint,mT,Mean);
      g_GhettoFemto_mT_rstar.SetPointError(uPoint,0,Err);
      uPoint++;
    }
  }

  const double MeanToGauss = 1./2.3;
  double mT_first = g_GhettoFemto_mT_rstar.GetPointY(0);
  printf("r(%.2f GeV) = %.2f\n",g_GhettoFemto_mT_rstar.GetPointX(0)*0.001,mT_first*MeanToGauss);
  double mT_2GeV = g_GhettoFemto_mT_rstar.Eval(2000);
  printf("r(%.2f GeV) = %.2f\n",2.,mT_2GeV*MeanToGauss);

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
    //if(kstar>780) continue;
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
  h_Ghetto_rstar->Write();
  h_GhettoFemto_rstar->Write();
  fSource->Write();
  h_Ghetto_kstar->Write();
  h_Ghetto_kstar_rstar->Write();
  h_Ghetto_mT_rstar->Write();
  h_GhettoFemto_mT_rstar->Write();
  g_GhettoFemto_mT_rstar.Write();
  h_GhettoFemto_mT_kstar->Write();
  h_Ghetto_mT_costh->Write();

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
  for(unsigned uMom=0; uMom<h_Ghetto_kstar_rstar->GetXaxis()->GetNbins(); uMom++){
    if(hkstar_rstar[uMom]){
      hkstar_rstar[uMom]->Write();
      delete hkstar_rstar[uMom];
    }
  }


  delete h_pT_p_all;
  delete h_pT_d_all;
  delete hSampleQA_p;

  delete h_Ghetto_rstar;
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
  */
}






void pd_rstar_vs_KstarTau_CreateTH2F(const TString type = "pd"){
  TString BaseName1 = TString::Format("%s/FunWithCeca/Ceca_pd_CoreFix/Ceca_pd_ET1_PR1_",GetFemtoOutputFolder());
  TString BaseName2 = TString::Format("%s/FunWithCeca/Ceca_pd_EffFix/Ceca_pd_ET1_PR1_",GetFemtoOutputFolder());
  TString OutputFile = TString::Format("%s/FunWithCeca/Ceca_pd_rstar_vs_KstarTau_040322.root",GetFemtoOutputFolder());
  if(type=="Kd"){
    BaseName1 = TString::Format("%s/FunWithCeca/Ceca_Kd_CoreFix/Ceca_Kd_ET1_PR1_",GetFemtoOutputFolder());
    BaseName2 = TString::Format("%s/FunWithCeca/Ceca_Kd_EffFix/Ceca_Kd_ET1_PR1_",GetFemtoOutputFolder());
    OutputFile = TString::Format("%s/FunWithCeca/Ceca_Kd_rstar_vs_KstarTau_040322.root",GetFemtoOutputFolder());
  }


  const unsigned NumMomBins = 60;
  const double kMin = 0;
  const double kMax = 800;
  const unsigned NumTauBins = 25;
  const double tMin = -2.5;
  const double tMax = 122.5;
  //const unsigned NumTauBins = 7;
  //const double tMin = -10;
  //const double tMax = 130.0;
  TH2F* rstar_vs_KstarTau = new TH2F("rstar_vs_KstarTau","rstar_vs_KstarTau",NumMomBins,kMin,kMax,NumTauBins,tMin,tMax);

  TFile* fInput1 = NULL;
  TGraph* gGhettoRadKstar1 = NULL;
  TGraph* gRadKstar1 = NULL;
  TFile* fInput2 = NULL;
  TGraph* gGhettoRadKstar2 = NULL;
  TGraph* gRadKstar2 = NULL;

  double tau,kstar,rad,rad_min,rad_max;
  for(unsigned uTau=0; uTau<NumTauBins; uTau++){
    tau = rstar_vs_KstarTau->GetYaxis()->GetBinCenter(uTau+1);

    if(fInput1) {delete fInput1; fInput1=NULL;}
    fInput1 = new TFile(BaseName1+TString::Format("DD%i.root",TMath::Nint(tau)),"read");
    gGhettoRadKstar1 = (TGraph*)fInput1->Get("gGhettoRadKstar");
    gRadKstar1 = (TGraph*)fInput1->Get("gRadKstar");

    if(fInput2) {delete fInput2; fInput2=NULL;}
    fInput2 = new TFile(BaseName2+TString::Format("DD%i.root",TMath::Nint(tau)),"read");
    gGhettoRadKstar2 = (TGraph*)fInput2->Get("gGhettoRadKstar");
    gRadKstar2 = (TGraph*)fInput2->Get("gRadKstar");

    for(unsigned uKstar=0; uKstar<NumMomBins; uKstar++){
      kstar = rstar_vs_KstarTau->GetXaxis()->GetBinCenter(uKstar+1);
      //printf("k*=%.0f; t=%.0f\n",kstar,tau);
      rad_min = 1e64;
      rad_max = 0;

      rad = gRadKstar1->Eval(kstar);
      if(rad>rad_max)rad_max=rad;
      if(rad<rad_min)rad_min=rad;
      //printf("  %.2f %.2f %.2f\n",rad,rad_min,rad_max);

      rad = gGhettoRadKstar1->Eval(kstar);
      if(rad>rad_max)rad_max=rad;
      if(rad<rad_min)rad_min=rad;
      //printf("  %.2f %.2f %.2f\n",rad,rad_min,rad_max);

      rad = gRadKstar2->Eval(kstar);
      if(rad>rad_max)rad_max=rad;
      if(rad<rad_min)rad_min=rad;
      //printf("  %.2f %.2f %.2f\n",rad,rad_min,rad_max);

      rad = gGhettoRadKstar2->Eval(kstar);
      if(rad>rad_max)rad_max=rad;
      if(rad<rad_min)rad_min=rad;
      //printf("  %.2f %.2f %.2f\n",rad,rad_min,rad_max);
      //printf("\n");

      rstar_vs_KstarTau->SetBinContent(uKstar+1,uTau+1,(rad_max+rad_min)*0.5);
      rstar_vs_KstarTau->SetBinError(uKstar+1,uTau+1,(rad_max-rad_min)*0.5);
    }
  }

  TFile* fOutput = new TFile(OutputFile,"recreate");
  rstar_vs_KstarTau->Write();

  delete rstar_vs_KstarTau;
  if(fInput1) {delete fInput1; fInput1=NULL;}
  if(fInput2) {delete fInput2; fInput2=NULL;}
  delete fOutput;

}

DLM_Histo<float>* dlmHist_pd_fit=NULL;
double Fit_pd_rstar_kstar(double* x, double* par){
  if(!dlmHist_pd_fit) return 0;
  if(dlmHist_pd_fit->GetDim()!=2) return 0;
  double kstar_tau[2];
  kstar_tau[0] = *x;
  kstar_tau[1] = par[0];
  return dlmHist_pd_fit->Eval(kstar_tau);
}

TGraphErrors* pd_rstar_kstar(const double& tau){
  TGraphErrors* gResult = new TGraphErrors();
  gResult->SetName(TString::Format("pd_rstar_kstar_t%.0f",tau));
  double rad;
  double err;
  double kstar_tau[2];
  kstar_tau[1] = tau;
  for(unsigned uKstar=0; uKstar<dlmHist_pd_fit->GetNbins(0); uKstar++){
    kstar_tau[0] = dlmHist_pd_fit->GetBinCenter(0,uKstar);
    rad = dlmHist_pd_fit->Eval(kstar_tau);
    err = dlmHist_pd_fit->Eval(kstar_tau,true);
    gResult->SetPoint(uKstar,kstar_tau[0],rad);
    gResult->SetPointError(uKstar,0,err);
  }
  return gResult;
}


void Ceca_mT_1(){

  const double HadronSize = 0.0;
  const double HadronSlope = 0.0;
  const double EtaCut = 0.8;
  const bool PROTON_RESO = true;
  const bool EQUALIZE_TAU = true;
  const double TIMEOUT = 120;
  const double TempPP = 154;//171
  const double Asym = 6;
  const double PancakeFluct = 0.00;
  const double DisplT = 0.8;
  const double DisplZ = 0.0;
  unsigned THREADS = 8;


  const double VolumeHI = 5000;
  const double TempHI = 154;
  const double ThKick = TempPP/sqrt(3.);
  const double EngyHI = 208.*5.02e6;
  const double ProtonEngy = 6.5e6;
  const double EngyPP = 2.*ProtonEngy;
  const double EngyDenst_PP_HI = 0.029412*TempPP-3.5294;//ratio
  const double VolumePP = EngyPP/EngyHI*VolumeHI/EngyDenst_PP_HI;
  const double ProtonGamma = ProtonEngy/938.;

  const double PancakeT = pow(3.*VolumePP*Asym/4./Pi,1./3.);
  const double PancakeZ = PancakeT/Asym;

  printf("The pp volume is %.2f fm^3\n",VolumePP);
  printf("Pancake  : %.2f x %.2f x %.2f\n", PancakeT,PancakeT,PancakeZ);

  TREPNI Database(0);
  Database.SetSeed(11);
  std::vector<TreParticle*> ParticleList;
  ParticleList.push_back(Database.NewParticle("Proton"));
  ParticleList.push_back(Database.NewParticle("Pion"));
  ParticleList.push_back(Database.NewParticle("ProtonReso"));

  TFile fOutput(TString::Format("%s/FunWithCeca/Ceca_mT_ET%i_PR%i.root",GetFemtoOutputFolder(),EQUALIZE_TAU,PROTON_RESO),"recreate");

  DLM_Histo<float>* dlm_pT_p = NULL;
  DLM_Histo<float> dlm_pT_eta_p;
  TH1F* h_pT_p_all = NULL;
  TH1F* h_pT_p = NULL;
  TH1F* h_pT_ap = NULL;

  const TString FilePath = TString::Format("%s/CatsFiles/Source/CECA/pd/pT_spectra/",GetCernBoxDimi());

  TFile file_p(FilePath+"protonSpectra.root","read");
  h_pT_p = (TH1F*)file_p.Get("pTDist_after");
  if(!h_pT_p) printf("ISSUE with h_pT_p\n");
  fOutput.cd();
  h_pT_p_all = (TH1F*)h_pT_p->Clone("h_pT_p_all");

  TFile file_ap(FilePath+"AntiProtonSpectra.root","read");
  h_pT_ap = (TH1F*)file_ap.Get("pTDist_after");
  if(!h_pT_ap) printf("ISSUE with h_pT_ap\n");
  h_pT_p_all->Add(h_pT_ap);

  dlm_pT_p = Convert_TH1F_DlmHisto(h_pT_p_all);
  dlm_pT_p->RescaleAxis(0,1000,false);

  double* BinRange = NULL;
  double axis[2];

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

  fOutput.cd();
  h_pT_p_all->Write();

  for(TreParticle* prt : ParticleList){
    if(prt->GetName()=="Proton"){
      prt->SetMass(Mass_p);
      prt->SetAbundance(35.78+64.22*(!PROTON_RESO));
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);
      prt->SetPtEtaPhi(dlm_pT_eta_p);
    }
    else if(prt->GetName()=="Pion"){
      prt->SetMass(Mass_pic);
      prt->SetAbundance(0);
    }
    else if(prt->GetName()=="ProtonReso"){
      prt->SetMass(1362);
      prt->SetAbundance(64.22*PROTON_RESO);
      prt->SetWidth(hbarc/1.65);
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);
      prt->SetPtEtaPhi(dlm_pT_eta_p);

      prt->NewDecay();
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Proton"));
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
      prt->GetDecay(0)->SetBranching(100);
    }

  }

  std::vector<std::string> ListOfParticles;
  ListOfParticles.push_back("Proton");
  ListOfParticles.push_back("Proton");
  CECA Ivana(Database,ListOfParticles);

/*
  Ivana.SetDisplacement(rSP_core);
  Ivana.SetHadronization(0.0);
  Ivana.SetHadrFluctuation(0.0);
  Ivana.SetTau(0.0);
  Ivana.SetThermalKick(0);
*/

  Ivana.SetDisplacementT(DisplT);
  Ivana.SetDisplacementZ(DisplZ);
  Ivana.SetHadronizationT(PancakeT);
  Ivana.SetHadronizationZ(PancakeZ);
  Ivana.SetHadrFluctuation(PancakeFluct);
  Ivana.SetTau(0.0);
  Ivana.SetThermalKick(0.0);

  Ivana.SetTargetStatistics(10);
  Ivana.SetEventMult(2);
  Ivana.SetSourceDim(2);
  Ivana.SetDebugMode(true);
  Ivana.SetThreadTimeout(TIMEOUT);
  Ivana.EqualizeFsiTime(EQUALIZE_TAU);
  Ivana.GHETTO_EVENT = true;

  Ivana.GoBabyGo(THREADS);

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

  TH1F* h_Ghetto_rstar = Convert_DlmHisto_TH1F(Ivana.Ghetto_rstar,"Ghetto_rstar");
  TH1F* h_Ghetto_kstar = Convert_DlmHisto_TH1F(Ivana.Ghetto_kstar,"Ghetto_kstar");
  Ivana.Ghetto_kstar_rstar->ComputeError();
  TH2F* h_Ghetto_kstar_rstar = Convert_DlmHisto_TH2F(Ivana.Ghetto_kstar_rstar,"Ghetto_kstar_rstar");
  TH2F* h_Ghetto_mT_rstar = Convert_DlmHisto_TH2F(Ivana.Ghetto_mT_rstar,"Ghetto_mT_rstar");
  Ivana.GhettoFemto_mT_rstar->ComputeError();
  TH2F* h_GhettoFemto_mT_rstar = Convert_DlmHisto_TH2F(Ivana.GhettoFemto_mT_rstar,"GhettoFemto_mT_rstar");
  TH2F* h_GhettoFemto_mT_kstar = Convert_DlmHisto_TH2F(Ivana.GhettoFemto_mT_kstar,"GhettoFemto_mT_kstar");
  TH2F* h_Ghetto_mT_costh = Convert_DlmHisto_TH2F(Ivana.Ghetto_mT_costh,"Ghetto_mT_costh");
  TH2F* h_GhettoSP_pT_th = Convert_DlmHisto_TH2F(Ivana.GhettoSP_pT_th,"GhettoSP_pT_th");
  TH1F* h_GhettoSP_pT_1 = Convert_DlmHisto_TH1F(Ivana.GhettoSP_pT_1,"GhettoSP_pT_1");
  TH1F* h_GhettoSP_pT_2 = Convert_DlmHisto_TH1F(Ivana.GhettoSP_pT_2,"GhettoSP_pT_2");
  Ivana.GhettoFemto_rstar->ComputeError();
  TH1F* h_GhettoFemto_rstar = Convert_DlmHisto_TH1F(Ivana.GhettoFemto_rstar,"GhettoFemto_rstar");
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

  TGraphErrors g_GhettoFemto_mT_rstar;
  g_GhettoFemto_mT_rstar.SetName("g_GhettoFemto_mT_rstar");
  g_GhettoFemto_mT_rstar.SetMarkerStyle(20);
  g_GhettoFemto_mT_rstar.SetMarkerSize(2);
  g_GhettoFemto_mT_rstar.SetLineWidth(3);

  TFile fMt(TString::Format("%s/SourceStudies/SourcePaper_pp/Systematics/gMt_pp.root",GetFemtoOutputFolder()),"read");
  TGraphErrors* gMt_pp = (TGraphErrors*)fMt.Get("gMtStatGauss_pp");

  unsigned uPoint = 0;
  for(unsigned uBin=0; uBin<h_GhettoFemto_mT_rstar->GetXaxis()->GetNbins(); uBin++){
    TH1D* hProj = h_GhettoFemto_mT_rstar->ProjectionY(TString::Format("hProj"),uBin+1,uBin+1);
    if(hProj->GetEntries()<200) {delete hProj; continue;}
    hProj->Scale(1./hProj->Integral(),"width");
    double lowerlimit;
    double upperlimit;
    GetCentralInterval(*h_GhettoFemto_rstar, 0.9, lowerlimit, upperlimit, true);
    TF1* fSource_mT = new TF1("fSource_mT",ScaledGauss,lowerlimit,upperlimit,2);
    fSource_mT->SetParameter(0,hProj->GetMean()/2.3);
    fSource_mT->SetParLimits(0,0.1,10.0);
    fSource_mT->FixParameter(1,1);

    hProj->Fit(fSource_mT,"Q, S, N, R, M");
    double mT = h_GhettoFemto_mT_rstar->GetXaxis()->GetBinCenter(uBin+1);
    g_GhettoFemto_mT_rstar.SetPoint(uPoint,mT*0.001,fSource_mT->GetParameter(0));
    uPoint++;
    delete hProj;
    delete fSource_mT;
  }

  double chi2=0;
  for(unsigned uBin=0; uBin<gMt_pp->GetN(); uBin++){
    double mt,r0,dr0,rCeca;
    gMt_pp->GetPoint(uBin,mt,r0);
    dr0 = gMt_pp->GetErrorY(uBin);
    rCeca = g_GhettoFemto_mT_rstar.Eval(mt);
    chi2 += pow((r0-rCeca)/dr0,2.);
  }
  double nsigma = sqrt(2)*TMath::ErfcInverse(TMath::Prob(chi2,gMt_pp->GetN()));
  printf("The deviation from the data r(mT) is %.1f σ (χ2/ndf=%.1f)\n",nsigma,chi2/double(gMt_pp->GetN()));

  const double MeanToGauss = 1./2.3;
  double mT_first = g_GhettoFemto_mT_rstar.GetPointY(0);
  printf("r(%.2f GeV) = %.2f\n",g_GhettoFemto_mT_rstar.GetPointX(0)*0.001,mT_first*MeanToGauss);
  double mT_2GeV = g_GhettoFemto_mT_rstar.Eval(2000);
  printf("r(%.2f GeV) = %.2f\n",2.,mT_2GeV*MeanToGauss);

  h_GhettoFemto_rstar->Scale(1./h_GhettoFemto_rstar->Integral(),"width");

  double lowerlimit;
  double upperlimit;
  GetCentralInterval(*h_GhettoFemto_rstar, 0.9, lowerlimit, upperlimit, true);
  unsigned lowerbin = h_GhettoFemto_rstar->FindBin(lowerlimit);
  unsigned upperbin = h_GhettoFemto_rstar->FindBin(upperlimit);
  printf("The fit will be performed in the range [%.2f, %.2f] fm\n",lowerlimit,upperlimit);
  TF1* fSource = new TF1("fSource",ScaledGauss,lowerlimit,upperlimit,2);
  fSource->SetParameter(0,h_GhettoFemto_rstar->GetMean()/2.3);
  fSource->SetParLimits(0,0.5,5.0);
  //fSource->SetParameter(1,0.5);
  //fSource->SetParLimits(1,0.1,1.0);
  fSource->FixParameter(1,1.0);

  h_GhettoFemto_rstar->Fit(fSource,"Q, S, N, R, M");
  printf("The effective Gaussian size is %.3f +/- %.3f fm, λ = %.3f +/- %.3f\n",
  fSource->GetParameter(0),fSource->GetParError(0),fSource->GetParameter(1),fSource->GetParError(1));

  TGraph gRadKstar;
  gRadKstar.SetName("gRadKstar");
  gRadKstar.SetLineColor(kRed+1);
  gRadKstar.SetLineWidth(3);

  TH1D** hkstar_rstar;
  hkstar_rstar = new TH1D* [h_Ghetto_kstar_rstar->GetXaxis()->GetNbins()];

  for(unsigned uMom=0; uMom<h_Ghetto_kstar_rstar->GetXaxis()->GetNbins(); uMom++){
    double kstar = h_Ghetto_kstar_rstar->GetXaxis()->GetBinCenter(uMom+1);
    hkstar_rstar[uMom] = NULL;
    if(kstar>200) continue;
    hkstar_rstar[uMom] = h_Ghetto_kstar_rstar->ProjectionY(TString::Format("hkstar_rstar_%.0f",kstar),uMom+1,uMom+1);
    hkstar_rstar[uMom]->Scale(1./hkstar_rstar[uMom]->Integral(),"width");
    GetCentralInterval(*hkstar_rstar[uMom], 0.9, lowerlimit, upperlimit, true);
    //printf("--- k* = %.0f ---\n",kstar);
    //printf(" The fit will be performed in the range [%.2f, %.2f] fm\n",lowerlimit,upperlimit);
    fSource->SetParameter(0,hkstar_rstar[uMom]->GetMean()/2.3);
    fSource->SetParLimits(0,hkstar_rstar[uMom]->GetMean()/4.,hkstar_rstar[uMom]->GetMean());
    hkstar_rstar[uMom]->Fit(fSource,"Q, S, N, R, M");
    //printf(" The effective Gaussian size is %.3f +/- %.3f fm, λ = %.3f +/- %.3f\n",
    //fSource->GetParameter(0),fSource->GetParError(0),fSource->GetParameter(1),fSource->GetParError(1));
    //gRadKstar.SetPoint(uMom,kstar,fSource->GetParameter(0));
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
  h_Ghetto_rstar->Write();
  h_GhettoFemto_rstar->Write();
  fSource->Write();
  h_Ghetto_kstar->Write();
  h_Ghetto_kstar_rstar->Write();
  h_Ghetto_mT_rstar->Write();
  h_GhettoFemto_mT_rstar->Write();
  g_GhettoFemto_mT_rstar.Write();
  gMt_pp->Write();
  h_GhettoFemto_mT_kstar->Write();
  h_Ghetto_mT_costh->Write();

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
  for(unsigned uMom=0; uMom<h_Ghetto_kstar_rstar->GetXaxis()->GetNbins(); uMom++){
    if(hkstar_rstar[uMom]){
      hkstar_rstar[uMom]->Write();
    }
  }


  delete h_pT_p_all;
}



DLM_Histo<float>* GetPtEtaDist(TString FileNameP, TString FileNameAP,
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

  return dlm_pT_eta;
}




void Ceca_vs_RSM_1(const std::string part1, const std::string part2){
  const double reff = sqrt(2);
  const double HadronSize = 0.75*0;
  const double HadronSlope = 0.2*0;
  const double EtaCut = 0.8;
  const bool EQUALIZE_TAU = true;
  const double TIMEOUT = 30;
  const double TempPP = 171;//171 or 154
  const double Asym = 40;
  //const double PancakeT = 7.06*1.0;
  //const double PancakeZ = 0.5*3;
  const double PancakeT = 2.3*1.0;
  const double PancakeZ = 3.85*0.0;
  const double Tau = 3.55*1.0;
  const double PancakeFluct = 0;//in % !!!
  //const double DisplT = 0.18*0;
  //const double DisplZ = 0.10*0;
  const double DisplT = 0.2*1.0;
  const double DisplZ = 0.2*1.0;
  //const double DisplZ = DisplT/Asym;
  unsigned THREADS = 1;
  const double MomSpread = 1061;
  const bool THERMAL_KICK = false;

  const double Sensitivity = 0.0033*10000000;

  const double VolumeHI = 5000;
  const double TempHI = 154;
  const double ThKick = TempPP/sqrt(3.);
  const double EngyHI = 208.*5.02e6;
  const double ProtonEngy = 6.5e6;
  const double EngyPP = 2.*ProtonEngy;
  const double EngyDenst_PP_HI = 0.029412*TempPP-3.5294;//ratio
  const double VolumePP = TempPP?EngyPP/EngyHI*VolumeHI/EngyDenst_PP_HI:0;
  const double ProtonGamma = ProtonEngy/938.;
  //const double PancakeT = pow(3.*VolumePP*Asym/4./Pi,1./3.);
  //const double PancakeZ = PancakeT/Asym;

  //const double PancakeT = pow(3.*VolumePP*ProtonGamma/4./Pi,1./3.);
  //const double PancakeZ = PancakeT/ProtonGamma;

  const double q_CutOff = 200;

  double PancakeV = PancakeT*PancakeT*PancakeZ*4./3.*Pi;
  printf("The pp volume is %.2f fm^3\n",VolumePP);
  //printf("Pancake  : %.2f x %.2f x %.2f\n", PancakeT,PancakeT,PancakeZ);
  printf("Displacement: %.2f %.2f %.2f fm\n",DisplT,DisplT,DisplZ);
  printf("Pancake     : %.2f %.2f %.2f fm\n",PancakeT,PancakeT,PancakeZ);
  printf("PancakeV    : %.2f fm^3\n",PancakeV);


  const double Tau_Proton = 1.65;
  const double Tau_Lambda = 4.69;//4.69

  const double Mass_ProtonReso = 1362;
  const double Mass_LambdaReso = 1462;

  const double FracProtonReso = 0.6422;
  const double FracLambdaReso = 0.6438;

  const unsigned NumRadBins = 256;
  const double rMin = 0;
  const double rMax = 32;



  TString FilePath = TString::Format("%s/Jaime/",GetCernBoxDimi());
  TString FileNameP1,FileNameAP1,FileNameP2,FileNameAP2;
  FileNameP1 = "p_pT.root";
  FileNameAP1 = "ap_pT.root";
  FileNameP2 = "L_pT.root";
  FileNameAP2 = "aL_pT.root";


  TH1F* h_pT_p;
  TH1F* h_pT_ap;
  TH1F* h_pT_p_all;
  TH1F* h_pT_d;
  TH1F* h_pT_ad;
  TH1F* h_pT_d_all;


  TFile file_p(FilePath+FileNameP1,"read");
  h_pT_p = (TH1F*)file_p.Get("pTDist_after");
  if(!h_pT_p) printf("ISSUE with h_pT_p\n");
  //fOutput.cd();
  gROOT->cd();
  h_pT_p_all = (TH1F*)h_pT_p->Clone("h_pT_p_all");

  TFile file_ap(FilePath+FileNameAP1,"read");
  h_pT_ap = (TH1F*)file_ap.Get("pTDist_after");
  if(!h_pT_ap) printf("ISSUE with h_pT_ap\n");
  h_pT_p_all->Add(h_pT_ap);

  TFile file_d(FilePath+FileNameP2,"read");
  h_pT_d = (TH1F*)file_d.Get("pTDist_after");
  if(!h_pT_d) printf("ISSUE with h_pT_d\n");
  //fOutput.cd();
  gROOT->cd();
  h_pT_d_all = (TH1F*)h_pT_d->Clone("h_pT_d_all");

  TFile file_ad(FilePath+FileNameAP2,"read");
  h_pT_ad = (TH1F*)file_ad.Get("pTDist_after");
  if(!h_pT_ad) printf("ISSUE with h_pT_ad\n");
  h_pT_d_all->Add(h_pT_ad);

  DLM_Histo<float>* dlm_pT_p = Convert_TH1F_DlmHisto(h_pT_p_all);
  DLM_Histo<float>* dlm_pT_d = Convert_TH1F_DlmHisto(h_pT_d_all);
  dlm_pT_p->RescaleAxis(0,1000,false);
  dlm_pT_d->RescaleAxis(0,1000,false);

  double* BinRange = NULL;
  double axis[2];

  DLM_Histo<float> dlm_pT_eta_p;
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

  DLM_Histo<float> dlm_pT_eta_d;
  dlm_pT_eta_d.SetUp(2);
  BinRange = dlm_pT_d->GetBinRange(0);
  dlm_pT_eta_d.SetUp(0,dlm_pT_d->GetNbins(),BinRange);
  delete [] BinRange;
  dlm_pT_eta_d.SetUp(1,1,-EtaCut,EtaCut);
  dlm_pT_eta_d.Initialize();
  for(unsigned uBin=0; uBin<dlm_pT_d->GetNbins(); uBin++){
    dlm_pT_eta_d.SetBinContent(uBin,0,dlm_pT_d->GetBinContent(uBin));
  }

  gROOT->cd();

  TH1F* h_reff = new TH1F("h_reff","h_reff",NumRadBins,rMin,rMax);
  //TH1F* h_rstar_Ceca = new TH1F("h_rstar_Ceca","h_rstar_Ceca",NumRadBins,rMin,rMax);
  TH1F* h_rstar_Ceca = NULL;
  TH1F* h_rstar_Rsm = new TH1F("h_rstar_Rsm","h_rstar_Rsm",NumRadBins,rMin,rMax);
  TH1F* h_rcore_Rsm = new TH1F("h_rcore_Rsm","h_rcore_Rsm",NumRadBins,rMin,rMax);
  TH1F* h_rcore_Ceca;// = new TH1F("h_rcore_Ceca","h_rcore_Ceca",NumRadBins,rMin,rMax);

  TH1F* h_PR_AngleRcP2_Rsm = new TH1F("h_PR_AngleRcP2_Rsm","h_PR_AngleRcP2_Rsm",64,0,Pi);
  TH1F* h_RP_AngleRcP1_Rsm = new TH1F("h_RP_AngleRcP1_Rsm","h_RP_AngleRcP1_Rsm",64,0,Pi);
  TH1F* h_RR_AngleRcP2_Rsm = new TH1F("h_RR_AngleRcP2_Rsm","h_RR_AngleRcP2_Rsm",64,0,Pi);
  TH1F* h_RR_AngleRcP1_Rsm = new TH1F("h_RR_AngleRcP1_Rsm","h_RR_AngleRcP1_Rsm",64,0,Pi);
  TH1F* h_RR_AngleP1P2_Rsm = new TH1F("h_RR_AngleP1P2_Rsm","h_RR_AngleP1P2_Rsm",64,0,Pi);

  TF1* f_reff = new TF1("f_reff","[0]*4.*TMath::Pi()*x*x*pow(4.*TMath::Pi()*[1]*[1],-1.5)*exp(-(x*x)/(4.*[1]*[1]))+1.-[0]",rMin,rMax);
  TF1* f_rstar_Ceca = new TF1("f_rstar_Ceca","[0]*4.*TMath::Pi()*x*x*pow(4.*TMath::Pi()*[1]*[1],-1.5)*exp(-(x*x)/(4.*[1]*[1]))+1.-[0]",rMin,rMax);
  TF1* f_rstar_Rsm = new TF1("f_rstar_Rsm","[0]*4.*TMath::Pi()*x*x*pow(4.*TMath::Pi()*[1]*[1],-1.5)*exp(-(x*x)/(4.*[1]*[1]))+1.-[0]",rMin,rMax);
  TF1* f_rcore_Ceca = new TF1("f_rcore_Ceca","[0]*4.*TMath::Pi()*x*x*pow(4.*TMath::Pi()*[1]*[1],-1.5)*exp(-(x*x)/(4.*[1]*[1]))+1.-[0]",rMin,rMax);
  TF1* f_rcore_Rsm = new TF1("f_rcore_Rsm","[0]*4.*TMath::Pi()*x*x*pow(4.*TMath::Pi()*[1]*[1],-1.5)*exp(-(x*x)/(4.*[1]*[1]))+1.-[0]",rMin,rMax);


  f_reff->FixParameter(0,1);
  f_reff->FixParameter(1,reff);
  for(unsigned uBin=0; uBin<NumRadBins; uBin++){
    double rad = h_reff->GetBinCenter(uBin+1);
    double val = f_reff->Eval(rad);
    double err = 1e-3;
    h_reff->SetBinContent(uBin+1,val);
    h_reff->SetBinError(uBin+1,err);
  }


  //set up the RSM
  DLM_CleverMcLevyResoTM MagicSource;
  MagicSource.InitStability(1,2-1e-6,2+1e-6);
  MagicSource.InitScale(38,0.15,2.0);
  MagicSource.InitRad(257*2,0,64);
  MagicSource.InitType(2);
  MagicSource.SetUpReso(0,FracProtonReso);
  MagicSource.SetUpReso(1,FracLambdaReso);
  MagicSource.InitNumMcIter(1000000);


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
  DLM_Random RanGen(11);
  double RanVal1;
  double RanVal2;
  double RanVal3;

  TFile* F_EposDisto_p_LamReso = new TFile(TString::Format("%s/CatsFiles/Source/EposAngularDist/EposDisto_p_LamReso.root",GetCernBoxDimi()));
  TNtuple* T_EposDisto_p_LamReso = (TNtuple*)F_EposDisto_p_LamReso->Get("InfoTuple_ClosePairs");
  unsigned N_EposDisto_p_LamReso = T_EposDisto_p_LamReso->GetEntries();
  T_EposDisto_p_LamReso->SetBranchAddress("k_D",&k_D);
  T_EposDisto_p_LamReso->SetBranchAddress("P1",&fP1);
  T_EposDisto_p_LamReso->SetBranchAddress("P2",&fP2);
  T_EposDisto_p_LamReso->SetBranchAddress("M1",&fM1);
  T_EposDisto_p_LamReso->SetBranchAddress("M2",&fM2);
  T_EposDisto_p_LamReso->SetBranchAddress("Tau1",&Tau1);
  T_EposDisto_p_LamReso->SetBranchAddress("Tau2",&Tau2);
  T_EposDisto_p_LamReso->SetBranchAddress("AngleRcP1",&AngleRcP1);
  T_EposDisto_p_LamReso->SetBranchAddress("AngleRcP2",&AngleRcP2);
  T_EposDisto_p_LamReso->SetBranchAddress("AngleP1P2",&AngleP1P2);
  for(unsigned uEntry=0; uEntry<N_EposDisto_p_LamReso; uEntry++){
      T_EposDisto_p_LamReso->GetEntry(uEntry);
      Tau1 = 0;
      Tau2 = Tau_Lambda;
      fM2 = Mass_LambdaReso;
      if(k_D>q_CutOff) continue;
      RanVal2 = RanGen.Exponential(fM2/(fP2*Tau2));
      MagicSource.AddBGT_PR(RanVal2,cos(AngleRcP2));
      h_PR_AngleRcP2_Rsm->Fill(AngleRcP2);
  }
  delete F_EposDisto_p_LamReso;

  TFile* F_EposDisto_pReso_Lam = new TFile(TString::Format("%s/CatsFiles/Source/EposAngularDist/EposDisto_pReso_Lam.root",GetCernBoxDimi()));
  TNtuple* T_EposDisto_pReso_Lam = (TNtuple*)F_EposDisto_pReso_Lam->Get("InfoTuple_ClosePairs");
  unsigned N_EposDisto_pReso_Lam = T_EposDisto_pReso_Lam->GetEntries();
  T_EposDisto_pReso_Lam->SetBranchAddress("k_D",&k_D);
  T_EposDisto_pReso_Lam->SetBranchAddress("P1",&fP1);
  T_EposDisto_pReso_Lam->SetBranchAddress("P2",&fP2);
  T_EposDisto_pReso_Lam->SetBranchAddress("M1",&fM1);
  T_EposDisto_pReso_Lam->SetBranchAddress("M2",&fM2);
  T_EposDisto_pReso_Lam->SetBranchAddress("Tau1",&Tau1);
  T_EposDisto_pReso_Lam->SetBranchAddress("Tau2",&Tau2);
  T_EposDisto_pReso_Lam->SetBranchAddress("AngleRcP1",&AngleRcP1);
  T_EposDisto_pReso_Lam->SetBranchAddress("AngleRcP2",&AngleRcP2);
  T_EposDisto_pReso_Lam->SetBranchAddress("AngleP1P2",&AngleP1P2);
  for(unsigned uEntry=0; uEntry<N_EposDisto_pReso_Lam; uEntry++){
      T_EposDisto_pReso_Lam->GetEntry(uEntry);
      Tau1 = Tau_Proton;
      Tau2 = 0;
      fM1 = Mass_ProtonReso;
      if(k_D>q_CutOff) continue;
      RanVal1 = RanGen.Exponential(fM1/(fP1*Tau1));
      MagicSource.AddBGT_RP(RanVal1,cos(AngleRcP1));
      h_RP_AngleRcP1_Rsm->Fill(AngleRcP1);
  }
  delete F_EposDisto_pReso_Lam;

  TFile* F_EposDisto_pReso_LamReso = new TFile(TString::Format("%s/CatsFiles/Source/EposAngularDist/EposDisto_pReso_LamReso.root",GetCernBoxDimi()));
  TNtuple* T_EposDisto_pReso_LamReso = (TNtuple*)F_EposDisto_pReso_LamReso->Get("InfoTuple_ClosePairs");
  unsigned N_EposDisto_pReso_LamReso = T_EposDisto_pReso_LamReso->GetEntries();
  T_EposDisto_pReso_LamReso->SetBranchAddress("k_D",&k_D);
  T_EposDisto_pReso_LamReso->SetBranchAddress("P1",&fP1);
  T_EposDisto_pReso_LamReso->SetBranchAddress("P2",&fP2);
  T_EposDisto_pReso_LamReso->SetBranchAddress("M1",&fM1);
  T_EposDisto_pReso_LamReso->SetBranchAddress("M2",&fM2);
  T_EposDisto_pReso_LamReso->SetBranchAddress("Tau1",&Tau1);
  T_EposDisto_pReso_LamReso->SetBranchAddress("Tau2",&Tau2);
  T_EposDisto_pReso_LamReso->SetBranchAddress("AngleRcP1",&AngleRcP1);
  T_EposDisto_pReso_LamReso->SetBranchAddress("AngleRcP2",&AngleRcP2);
  T_EposDisto_pReso_LamReso->SetBranchAddress("AngleP1P2",&AngleP1P2);
  for(unsigned uEntry=0; uEntry<N_EposDisto_pReso_LamReso; uEntry++){
      T_EposDisto_pReso_LamReso->GetEntry(uEntry);
      Tau1 = Tau_Proton;
      Tau2 = Tau_Lambda;
      fM1 = Mass_ProtonReso;
      fM2 = Mass_LambdaReso;
      if(k_D>q_CutOff) continue;
      RanVal1 = RanGen.Exponential(fM1/(fP1*Tau1));
      RanVal2 = RanGen.Exponential(fM2/(fP2*Tau2));
      MagicSource.AddBGT_RR(RanVal1,cos(AngleRcP1),RanVal2,cos(AngleRcP2),cos(AngleP1P2));
      h_RR_AngleRcP1_Rsm->Fill(AngleRcP1);
      h_RR_AngleRcP2_Rsm->Fill(AngleRcP2);
      h_RR_AngleP1P2_Rsm->Fill(AngleP1P2);
  }
  delete F_EposDisto_pReso_LamReso;

  h_PR_AngleRcP2_Rsm->Sumw2();
  h_RP_AngleRcP1_Rsm->Sumw2();
  h_RR_AngleRcP1_Rsm->Sumw2();
  h_RR_AngleRcP2_Rsm->Sumw2();
  h_RR_AngleP1P2_Rsm->Sumw2();

  h_PR_AngleRcP2_Rsm->Scale(1./h_PR_AngleRcP2_Rsm->Integral(),"width");
  h_RP_AngleRcP1_Rsm->Scale(1./h_RP_AngleRcP1_Rsm->Integral(),"width");
  h_RR_AngleRcP1_Rsm->Scale(1./h_RR_AngleRcP1_Rsm->Integral(),"width");
  h_RR_AngleRcP2_Rsm->Scale(1./h_RR_AngleRcP2_Rsm->Integral(),"width");
  h_RR_AngleP1P2_Rsm->Scale(1./h_RR_AngleP1P2_Rsm->Integral(),"width");

  ///////////////////////////

  TREPNI Database(0);
  Database.SetSeed(11);
  std::vector<TreParticle*> ParticleList;
  ParticleList.push_back(Database.NewParticle("Proton"));
  ParticleList.push_back(Database.NewParticle("Lambda"));
  ParticleList.push_back(Database.NewParticle("Pion"));

  ParticleList.push_back(Database.NewParticle("ProtonReso"));
  ParticleList.push_back(Database.NewParticle("LambdaReso"));

  double ExpPP,ExpPR,ExpRP,ExpRR;
  ExpPP = (1.-FracProtonReso)*(1.-FracLambdaReso);
  ExpPR = (1.-FracProtonReso)*FracLambdaReso;
  ExpRP = FracProtonReso*(1.-FracLambdaReso);
  ExpRR = FracProtonReso*FracLambdaReso;

  for(TreParticle* prt : ParticleList){
    if(prt->GetName()=="Proton"){
      prt->SetMass(Mass_p);
      prt->SetAbundance((1.-FracProtonReso));
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);

      //prt->SetPtPz(prt->GetMass()*MomSpread*0.001,prt->GetMass()*MomSpread*0.001);
      prt->SetPtEtaPhi(dlm_pT_eta_p);
    }
    else if(prt->GetName()=="Lambda"){
      prt->SetMass(Mass_L);
      prt->SetAbundance((1.-FracLambdaReso));
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);

      //prt->SetPtPz(prt->GetMass()*MomSpread*0.001,prt->GetMass()*MomSpread*0.001);
      prt->SetPtEtaPhi(dlm_pT_eta_d);
    }
    else if(prt->GetName()=="Pion"){
      prt->SetMass(Mass_pic);
      prt->SetAbundance(0);
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);

      prt->SetPtPz(prt->GetMass()*MomSpread*0.001,prt->GetMass()*MomSpread*0.001);
    }
    else if(prt->GetName()=="ProtonReso"){
      prt->SetMass(Mass_ProtonReso);
      prt->SetAbundance(FracProtonReso);
      prt->SetWidth(hbarc/Tau_Proton);
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);

      prt->NewDecay();
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Proton"));
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
      prt->GetDecay(0)->SetBranching(100);

      //prt->SetPtPz(prt->GetMass()*MomSpread*0.001,prt->GetMass()*MomSpread*0.001);
      prt->SetPtEtaPhi(dlm_pT_eta_p);
    }
    else if(prt->GetName()=="LambdaReso"){
      prt->SetMass(Mass_LambdaReso);
      prt->SetAbundance(FracLambdaReso);
      prt->SetWidth(hbarc/Tau_Lambda);
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);

      prt->NewDecay();
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Lambda"));
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
      prt->GetDecay(0)->SetBranching(100);

      //prt->SetPtPz(prt->GetMass()*MomSpread*0.001,prt->GetMass()*MomSpread*0.001);
      prt->SetPtEtaPhi(dlm_pT_eta_d);
    }

    prt->SetAcceptance_Eta(-EtaCut,EtaCut);
    //prt->SetPtPz(prt->GetMass()*MomSpread*0.001,prt->GetMass()*MomSpread*0.001);
    //prt->SetPtEtaPhi(MomSpread*);
    //
    //prt->SetAcceptance_pT(500,1e6);
  }
  std::vector<std::string> ListOfParticles;
  ListOfParticles.push_back(part1);
  ListOfParticles.push_back(part2);
  CECA Ivana(Database,ListOfParticles);
  Ivana.SetDisplacementT(DisplT);
  Ivana.SetDisplacementZ(DisplZ);
  Ivana.SetHadronizationT(PancakeT);
  Ivana.SetHadronizationZ(PancakeZ);
  Ivana.SetHadrFluctuation(PancakeFluct*0.01);
  Ivana.SetTau(Tau);
  Ivana.SetThermalKick(THERMAL_KICK*ThKick);
  Ivana.SetTargetStatistics(10);
  Ivana.SetEventMult(2);
  Ivana.SetSourceDim(2);
  Ivana.SetThreadTimeout(TIMEOUT);
  Ivana.SetFemtoRegion(q_CutOff*0.5);
  Ivana.EqualizeFsiTime(EQUALIZE_TAU);
  Ivana.GHETTO_EVENT = true;
  printf("Goal for reff = %.4f fm; <r*> = %.4f fm\n",reff,h_reff->GetMean());
  printf("--- CECA source ---\n");
  Ivana.GoBabyGo(THREADS);

  double r_fit;
  double adjust;

  Ivana.GhettoFemto_rstar->ComputeError();
  h_rstar_Ceca = Convert_DlmHisto_TH1F(Ivana.GhettoFemto_rstar,"h_rstar_Ceca");
  r_fit = Get_reff(h_rstar_Ceca);
  printf(" Starting parameters: Disp = %.2f %.2f fm; Hadr = %.2f %.2f fm\n",
  Ivana.GetDisplacementT(),Ivana.GetDisplacementZ(),Ivana.GetHadronizationT(),Ivana.GetHadronizationZ());
  double DEV = fabs(r_fit-reff)/fabs(reff);
  while(DEV>Sensitivity){
    adjust = reff/r_fit;
    Ivana.SetDisplacementT(Ivana.GetDisplacementT()*adjust);
    Ivana.SetDisplacementZ(Ivana.GetDisplacementZ()*adjust);
    Ivana.SetHadronizationT(Ivana.GetHadronizationT()*adjust);
    Ivana.SetHadronizationZ(Ivana.GetHadronizationZ()*adjust);
    Ivana.SetTau(Ivana.GetTau()*adjust);
    printf(" Redo with (DEV=%.2f%%): Disp = %.2f %.2f; Hadr = %.2f %.2f; Tau = %.2f\n",
    DEV*100.,Ivana.GetDisplacementT(),Ivana.GetDisplacementZ(),Ivana.GetHadronizationT(),Ivana.GetHadronizationZ(), Ivana.GetTau());
    Ivana.GoBabyGo(THREADS);
    delete h_rstar_Ceca;
    Ivana.GhettoFemto_rstar->ComputeError();
    h_rstar_Ceca = Convert_DlmHisto_TH1F(Ivana.GhettoFemto_rstar,"h_rstar_Ceca");
    r_fit = Get_reff(h_rstar_Ceca);
    DEV = fabs(r_fit-reff)/fabs(reff);
  }
  printf(" Final parameters: Disp = %.2f %.2f fm; Hadr = %.2f %.2f fm => reff = %.4f <r*> = %.4f\n",
  Ivana.GetDisplacementT(),Ivana.GetDisplacementZ(),Ivana.GetHadronizationT(),Ivana.GetHadronizationZ(),r_fit,h_rstar_Ceca->GetMean());
  h_rstar_Ceca->Scale(1./h_rstar_Ceca->Integral(),"width");

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
  printf("   Expected  Total  Femto\n");
  printf("PP %6.2f%% %6.2f%% %6.2f\n",ExpPP*100.,TotPP*100.,FemtoPP*100.);
  printf("PR %6.2f%% %6.2f%% %6.2f\n",ExpPR*100.,TotPR*100.,FemtoPR*100.);
  printf("RP %6.2f%% %6.2f%% %6.2f\n",ExpRP*100.,TotRP*100.,FemtoRP*100.);
  printf("RR %6.2f%% %6.2f%% %6.2f\n",ExpRR*100.,TotRR*100.,FemtoRR*100.);


  Ivana.Ghetto_RP_AngleRcP1->ComputeError();
  Ivana.Ghetto_PR_AngleRcP2->ComputeError();
  Ivana.Ghetto_RR_AngleRcP1->ComputeError();
  Ivana.Ghetto_RR_AngleRcP2->ComputeError();
  Ivana.Ghetto_RR_AngleP1P2->ComputeError();
  TH1F* h_RP_AngleRcP1_Ceca = Convert_DlmHisto_TH1F(Ivana.Ghetto_RP_AngleRcP1,"h_RP_AngleRcP1_Ceca");
  TH1F* h_PR_AngleRcP2_Ceca = Convert_DlmHisto_TH1F(Ivana.Ghetto_PR_AngleRcP2,"h_PR_AngleRcP2_Ceca");
  TH1F* h_RR_AngleRcP1_Ceca = Convert_DlmHisto_TH1F(Ivana.Ghetto_RR_AngleRcP1,"h_RR_AngleRcP1_Ceca");
  TH1F* h_RR_AngleRcP2_Ceca = Convert_DlmHisto_TH1F(Ivana.Ghetto_RR_AngleRcP2,"h_RR_AngleRcP2_Ceca");
  TH1F* h_RR_AngleP1P2_Ceca = Convert_DlmHisto_TH1F(Ivana.Ghetto_RR_AngleP1P2,"h_RR_AngleP1P2_Ceca");

  h_RP_AngleRcP1_Ceca->Scale(1./h_RP_AngleRcP1_Ceca->Integral(),"width");
  h_PR_AngleRcP2_Ceca->Scale(1./h_PR_AngleRcP2_Ceca->Integral(),"width");
  h_RR_AngleRcP1_Ceca->Scale(1./h_RR_AngleRcP1_Ceca->Integral(),"width");
  h_RR_AngleRcP2_Ceca->Scale(1./h_RR_AngleRcP2_Ceca->Integral(),"width");
  h_RR_AngleP1P2_Ceca->Scale(1./h_RR_AngleP1P2_Ceca->Integral(),"width");

  Ivana.GhettoFemto_mT_kstar->ComputeError();
  TH2F* h_GhettoFemto_mT_kstar = Convert_DlmHisto_TH2F(Ivana.GhettoFemto_mT_kstar,"GhettoFemto_mT_kstar");

  Ivana.GhettoFemto_mT_rstar->ComputeError();
  TH2F* h_GhettoFemto_mT_rstar = Convert_DlmHisto_TH2F(Ivana.GhettoFemto_mT_rstar,"GhettoFemto_mT_rstar");

  TGraphErrors g_GhettoFemto_mT_rstar;
  g_GhettoFemto_mT_rstar.SetName("g_GhettoFemto_mT_rstar");
  g_GhettoFemto_mT_rstar.SetMarkerStyle(20);
  g_GhettoFemto_mT_rstar.SetMarkerSize(2);
  g_GhettoFemto_mT_rstar.SetLineWidth(3);

  unsigned uPoint = 0;
  for(unsigned uBin=0; uBin<h_GhettoFemto_mT_rstar->GetXaxis()->GetNbins(); uBin++){
    TH1F* hProj = (TH1F*)h_GhettoFemto_mT_rstar->ProjectionY(TString::Format("hProj"),uBin+1,uBin+1);
    if(hProj->GetEntries()<200) {delete hProj; continue;}
    hProj->Scale(1./hProj->Integral(),"width");
    double mT = h_GhettoFemto_mT_rstar->GetXaxis()->GetBinCenter(uBin+1);
    g_GhettoFemto_mT_rstar.SetPoint(uPoint,mT*0.001,Get_reff(hProj));
    uPoint++;
    delete hProj;
  }

  TFile fMtEff_pp(TString::Format("%s/SourceStudies/SourcePaper_pp/Systematics/gMt_pp.root",GetFemtoOutputFolder()),"read");
  TGraphErrors* gMtEff_pp = (TGraphErrors*)fMtEff_pp.Get("gMtStatGauss_pp");
  //TGraphErrors* gMtCore_pp = (TGraphErrors*)fMtEff.Get("gMtStatCore_pp");
  TFile fMtEff_pL(TString::Format("%s/SourceStudies/SourcePaper_pL/Systematics/gMt_pL.root",GetFemtoOutputFolder()),"read");
  TGraphErrors* gMtEff_pL = (TGraphErrors*)fMtEff_pL.Get("gMtStatGauss_pL_NLO");
  //TGraphErrors* gMtCore_pp = (TGraphErrors*)fMtEff_pL.Get("gMtStatCore_pL_NLO");

  printf("--- CECA core ---\n");
  Database.GetParticle("ProtonReso")->SetAbundance(0);
  Database.GetParticle("LambdaReso")->SetAbundance(0);
  Ivana.SetThreadTimeout(TIMEOUT/4);
  Ivana.GoBabyGo(THREADS);
  Ivana.GhettoFemto_rstar->ComputeError();
  h_rcore_Ceca = Convert_DlmHisto_TH1F(Ivana.GhettoFemto_rstar,"h_rcore_Ceca");
  h_rcore_Ceca->Scale(1./h_rcore_Ceca->Integral(),"width");
  printf(" Final r_core = %.4f fm; reff = %.4f fm; <r*> = %.4f\n",Get_reff(h_rcore_Ceca),r_fit,h_rcore_Ceca->GetMean());

  Ivana.GhettoFemto_mT_rstar->ComputeError();
  TH2F* h_GhettoFemto_mT_core = Convert_DlmHisto_TH2F(Ivana.GhettoFemto_mT_rstar,"h_GhettoFemto_mT_core");

  TGraphErrors g_GhettoFemto_mT_core;
  g_GhettoFemto_mT_core.SetName("g_GhettoFemto_mT_core");
  g_GhettoFemto_mT_core.SetMarkerStyle(20);
  g_GhettoFemto_mT_core.SetMarkerSize(2);
  g_GhettoFemto_mT_core.SetLineWidth(3);

  uPoint = 0;
  for(unsigned uBin=0; uBin<h_GhettoFemto_mT_core->GetXaxis()->GetNbins(); uBin++){
    TH1F* hProj = (TH1F*)h_GhettoFemto_mT_core->ProjectionY(TString::Format("hProj"),uBin+1,uBin+1);
    if(hProj->GetEntries()<200) {delete hProj; continue;}
    hProj->Scale(1./hProj->Integral(),"width");
    double mT = h_GhettoFemto_mT_core->GetXaxis()->GetBinCenter(uBin+1);
    g_GhettoFemto_mT_core.SetPoint(uPoint,mT*0.001,Get_reff(hProj));
    uPoint++;
    delete hProj;
  }


  printf("--- RSM source ---\n");
  //double rcore = reff*0.819;
  //double rcore = reff*1.46;
  double rcore = reff;
  printf(" Starting parameter: rcore = %.4f fm\n",rcore);
  for(unsigned uBin=0; uBin<NumRadBins; uBin++){
    double rstar = h_rstar_Rsm->GetBinCenter(uBin+1);
    double parameters[2];
    parameters[0] = rcore;
    parameters[1] = 2.0;
    double val = MagicSource.RootEval(&rstar,parameters);
    h_rstar_Rsm->SetBinContent(uBin+1,val);
    h_rstar_Rsm->SetBinError(uBin+1,1e-3);
  }

  r_fit = Get_reff(h_rstar_Rsm);
  DEV = fabs(r_fit-reff)/fabs(reff);
  while(DEV>Sensitivity){
    rcore *= reff/r_fit;
    printf(" Redo with (DEV=%.2f%%): rcore = %.3f fm\n",DEV*100.,rcore);
    for(unsigned uBin=0; uBin<NumRadBins; uBin++){
      double rstar = h_rstar_Rsm->GetBinCenter(uBin+1);
      double parameters[2];
      parameters[0] = rcore;
      parameters[1] = 2.0;
      double val = MagicSource.RootEval(&rstar,parameters);
      h_rstar_Rsm->SetBinContent(uBin+1,val);
      h_rstar_Rsm->SetBinError(uBin+1,1e-3);
    }
    r_fit = Get_reff(h_rstar_Rsm);
    DEV = fabs(r_fit-reff)/fabs(reff);
  }

  printf(" Final rcore = %.4f fm; reff = %.4f fm; <r*> = %.4f\n",rcore,r_fit,h_rstar_Rsm->GetMean());
  h_rstar_Rsm->Scale(1./h_rstar_Rsm->Integral(),"width");

  printf("--- RSM core ---\n");
  f_rcore_Rsm->FixParameter(0,1);
  f_rcore_Rsm->FixParameter(1,rcore);
  for(unsigned uBin=0; uBin<NumRadBins; uBin++){
    double rstar = h_rcore_Rsm->GetBinCenter(uBin+1);
    h_rcore_Rsm->SetBinContent(uBin+1,f_rcore_Rsm->Eval(rstar));
    h_rcore_Rsm->SetBinError(uBin+1,1e-3);
  }

  TFile fOutput(TString::Format("%s/FunWithCeca/Ceca_%s_%s_vs_RSM_EqTau%i_Tau%.0f_Dsp%.0f-%.0f_Pnk%.0f-%.0f_Fct%.0f_Kick%.0f_Hds%.0f-%.0f.root",
  GetFemtoOutputFolder(),ListOfParticles.at(0).c_str(),ListOfParticles.at(1).c_str(),
  EQUALIZE_TAU,Tau*100.,Ivana.GetDisplacementT()*100.,Ivana.GetDisplacementZ()*100.,
  Ivana.GetHadronizationT()*100.,Ivana.GetHadronizationZ()*100.,Ivana.GetHadrFluctuation()*100.,
  ThKick*THERMAL_KICK,HadronSize*.100,HadronSlope*100.),"recreate");

  fOutput.cd();
  h_reff->Write();
  f_reff->Write();
  h_rstar_Ceca->Write();
  h_rstar_Rsm->Write();
  h_rcore_Ceca->Write();
  f_rcore_Ceca->Write();
  h_rcore_Rsm->Write();
  f_rcore_Rsm->Write();

  h_PR_AngleRcP2_Ceca->Write();
  h_PR_AngleRcP2_Rsm->Write();

  h_RP_AngleRcP1_Ceca->Write();
  h_RP_AngleRcP1_Rsm->Write();

  h_RR_AngleRcP1_Ceca->Write();
  h_RR_AngleRcP1_Rsm->Write();

  h_RR_AngleRcP2_Ceca->Write();
  h_RR_AngleRcP2_Rsm->Write();

  h_RR_AngleP1P2_Ceca->Write();
  h_RR_AngleP1P2_Rsm->Write();

  h_GhettoFemto_mT_rstar->Write();
  g_GhettoFemto_mT_rstar.Write();
  gMtEff_pp->Write();
  gMtEff_pL->Write();
  h_GhettoFemto_mT_core->Write();
  g_GhettoFemto_mT_core.Write();
  //gMtCore_pp->Write();

  h_GhettoFemto_mT_kstar->Write();

  delete h_reff;
  if(h_rstar_Ceca){delete h_rstar_Ceca; h_rstar_Ceca=NULL;}
  delete h_rstar_Rsm;
  delete h_rcore_Rsm;
  delete h_rcore_Ceca;
  delete h_PR_AngleRcP2_Rsm;
  delete h_RP_AngleRcP1_Rsm;
  delete h_RR_AngleRcP1_Rsm;
  delete h_RR_AngleRcP2_Rsm;
  delete h_RR_AngleP1P2_Rsm;
  delete h_PR_AngleRcP2_Ceca;
  delete h_RP_AngleRcP1_Ceca;
  delete h_RR_AngleRcP1_Ceca;
  delete h_RR_AngleRcP2_Ceca;
  delete h_RR_AngleP1P2_Ceca;
  delete h_GhettoFemto_mT_rstar;
  delete h_GhettoFemto_mT_kstar;
  delete f_reff;
  delete f_rstar_Ceca;
  delete f_rstar_Rsm;
  delete f_rcore_Ceca;
  delete f_rcore_Rsm;
}

void Ceca_vs_RSM_2(const std::string part1, const std::string part2, const bool PropMother, 
      const double mean_mt=0, const double reso_var=1){
  const double reff = sqrt(2);
  const double rSP = 0.176;//1.0,0.25*1.04
  //const double rSP = 0.92;
  const double FORCED_RCORE = 0;//0.98;

  //only concider this mt range in the S(r) disto.
  //To exclude the fit, put the range to be very large
  //Raffa: For p-p --->  mT = 1.30 GeV
  //Raffa: For p-Lambda ---> mT = 1.43 GeV
  double MeanMt = 1.43;
  if(mean_mt) MeanMt = mean_mt;
  //const double MeanMtRange = 1e16;
  const double MeanMtRange = 0.075;

  const double EtaCut = 0.8;
  const bool EQUALIZE_TAU = true;
  const double TIMEOUT = 30.*(1.);
  const double FracProtonReso = 0.6422*1 * reso_var;
  const double FracLambdaReso = 0.6438*1 * reso_var;
  const double PancakeT = 2.68;//2.3*1.04;
  //const double PancakeT = 7.5;
  const double PancakeZ = 0.0;
  const double Tau = 3.76;//3.55*1.04;
  const double PancakeFluct = 0;//in % !!!
  unsigned THREADS = 24;
  const double MomSpread = 1061;
  const bool THERMAL_KICK = false;
  const double q_CutOff = 100*2;
  const double DisplT = rSP;
  const double DisplZ = rSP;

  int Type = 0;
  if(part1=="Proton"&&part2=="Proton") Type = 1;
  else if(part1=="Proton"&&part2=="Lambda") Type = 2;
  else if(part1=="Proton"&&part2=="Kaon") Type = 3;
  else {printf("Cannot do it\n"); return;}

  double PancakeV = PancakeT*PancakeT*PancakeZ*4./3.*Pi;
  printf("Displacement: %.2f %.2f %.2f fm\n",DisplT,DisplT,DisplZ);
  printf("Pancake     : %.2f %.2f %.2f fm\n",PancakeT,PancakeT,PancakeZ);
  printf("PancakeV    : %.2f fm^3\n",PancakeV);

  const double Tau_Proton = 1.65;
  const double Tau_Lambda = 4.69;//4.69
  const double Tau_Kaon = 3.66;

  const double Mass_ProtonReso = 1362;
  const double Mass_LambdaReso = 1462;
  const double Mass_KaonReso = 1054;

  const unsigned NumRadBins = 256;
  const double rMin = 0;
  const double rMax = 32;

  DLM_Histo<float>* dlm_pT_eta_p = GetPtEtaDist(
    TString::Format("%s/Jaime/p_pT.root",GetCernBoxDimi()),
    TString::Format("%s/Jaime/ap_pT.root",GetCernBoxDimi()));

  DLM_Histo<float>* dlm_pT_eta_L = GetPtEtaDist(
    TString::Format("%s/Jaime/L_pT.root",GetCernBoxDimi()),
    TString::Format("%s/Jaime/aL_pT.root",GetCernBoxDimi()));

  TH1F* h_reff = new TH1F("h_reff","h_reff",NumRadBins,rMin,rMax);
  //TH1F* h_rstar_Ceca = new TH1F("h_rstar_Ceca","h_rstar_Ceca",NumRadBins,rMin,rMax);
  TH2F* h_kstar_rstar_Ceca;
  TH1F* h_rstar_CecaAll = NULL;
  TH1F* h_rstar_Ceca = NULL;
  TH1F* h_rstar_Rsm = new TH1F("h_rstar_Rsm","h_rstar_Rsm",NumRadBins,rMin,rMax);
  TH1F* h_rstar_RsmWS = new TH1F("h_rstar_RsmWS","h_rstar_RsmWS",NumRadBins,rMin,rMax);
  TH1F* h_rstar_RsmCeca = new TH1F("h_rstar_RsmCeca","h_rstar_RsmCeca",NumRadBins,rMin,rMax);
  TH1F* h_rcore_Rsm = new TH1F("h_rcore_Rsm","h_rcore_Rsm",NumRadBins,rMin,rMax);
  TH1F* h_rcore_Ceca;// = new TH1F("h_rcore_Ceca","h_rcore_Ceca",NumRadBins,rMin,rMax);
  TH1F* h_rcore_CecaAll;

  TH1F* h_PR_AngleRcP2_Rsm = new TH1F("h_PR_AngleRcP2_Rsm","h_PR_AngleRcP2_Rsm",64,0,Pi);
  TH1F* h_RP_AngleRcP1_Rsm = new TH1F("h_RP_AngleRcP1_Rsm","h_RP_AngleRcP1_Rsm",64,0,Pi);
  TH1F* h_RR_AngleRcP2_Rsm = new TH1F("h_RR_AngleRcP2_Rsm","h_RR_AngleRcP2_Rsm",64,0,Pi);
  TH1F* h_RR_AngleRcP1_Rsm = new TH1F("h_RR_AngleRcP1_Rsm","h_RR_AngleRcP1_Rsm",64,0,Pi);
  TH1F* h_RR_AngleP1P2_Rsm = new TH1F("h_RR_AngleP1P2_Rsm","h_RR_AngleP1P2_Rsm",64,0,Pi);

  TH1F* h_PR_AngleRcP2_RsmWS = new TH1F("h_PR_AngleRcP2_RsmWS","h_PR_AngleRcP2_RsmWS",64,0,Pi);
  TH1F* h_RP_AngleRcP1_RsmWS = new TH1F("h_RP_AngleRcP1_RsmWS","h_RP_AngleRcP1_RsmWS",64,0,Pi);

  TF1* f_reff = new TF1("f_reff","[0]*4.*TMath::Pi()*x*x*pow(4.*TMath::Pi()*[1]*[1],-1.5)*exp(-(x*x)/(4.*[1]*[1]))+1.-[0]",rMin,rMax);
  TF1* f_rstar_Ceca = new TF1("f_rstar_Ceca","[0]*4.*TMath::Pi()*x*x*pow(4.*TMath::Pi()*[1]*[1],-1.5)*exp(-(x*x)/(4.*[1]*[1]))+1.-[0]",rMin,rMax);
  TF1* f_rstar_Rsm = new TF1("f_rstar_Rsm","[0]*4.*TMath::Pi()*x*x*pow(4.*TMath::Pi()*[1]*[1],-1.5)*exp(-(x*x)/(4.*[1]*[1]))+1.-[0]",rMin,rMax);
  TF1* f_rcore_Ceca = new TF1("f_rcore_Ceca","[0]*4.*TMath::Pi()*x*x*pow(4.*TMath::Pi()*[1]*[1],-1.5)*exp(-(x*x)/(4.*[1]*[1]))+1.-[0]",rMin,rMax);
  TF1* f_rcore_Rsm = new TF1("f_rcore_Rsm","[0]*4.*TMath::Pi()*x*x*pow(4.*TMath::Pi()*[1]*[1],-1.5)*exp(-(x*x)/(4.*[1]*[1]))+1.-[0]",rMin,rMax);

  f_reff->FixParameter(0,1);
  f_reff->FixParameter(1,reff);
  for(unsigned uBin=0; uBin<NumRadBins; uBin++){
    double rad = h_reff->GetBinCenter(uBin+1);
    double val = f_reff->Eval(rad);
    double err = 1e-3;
    h_reff->SetBinContent(uBin+1,val);
    h_reff->SetBinError(uBin+1,err);
  }

  //set up the RSM
  DLM_CleverMcLevyResoTM MagicSource;
  MagicSource.InitStability(1,2-1e-6,2+1e-6);
  MagicSource.InitScale(38,0.15,2.0);
  MagicSource.InitRad(257*2,0,64);
  MagicSource.InitType(2);
  MagicSource.SetUpReso(0,FracProtonReso);
  MagicSource.SetUpReso(1,Type==1?FracProtonReso:FracLambdaReso);
  MagicSource.InitNumMcIter(1000000);

  DLM_CleverMcLevyResoTM MagicSourceWrongSign;
  MagicSourceWrongSign.InitStability(1,2-1e-6,2+1e-6);
  MagicSourceWrongSign.InitScale(38,0.15,2.0);
  MagicSourceWrongSign.InitRad(257*2,0,64);
  MagicSourceWrongSign.InitType(2);
  MagicSourceWrongSign.SetUpReso(0,FracProtonReso);
  MagicSourceWrongSign.SetUpReso(1,Type==1?FracProtonReso:FracLambdaReso);
  MagicSourceWrongSign.InitNumMcIter(1000000);

  DLM_CleverMcLevyResoTM MagicSourceCeca;
  MagicSourceCeca.InitStability(1,2-1e-6,2+1e-6);
  MagicSourceCeca.InitScale(38,0.15,2.0);
  MagicSourceCeca.InitRad(257*2,0,64);
  MagicSourceCeca.InitType(2);
  MagicSourceCeca.SetUpReso(0,FracProtonReso);
  MagicSourceCeca.SetUpReso(1,Type==1?FracProtonReso:FracLambdaReso);
  MagicSourceCeca.InitNumMcIter(1000000);

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
  DLM_Random RanGen(11);
  double RanVal1;
  double RanVal2;
  double RanVal3;

  if(Type==1){
    TFile* F_EposDisto_p_pReso;
    F_EposDisto_p_pReso = new TFile(TString::Format("%s/CatsFiles/Source/EposAngularDist/EposDisto_p_pReso.root",GetCernBoxDimi()));
    TNtuple* T_EposDisto_p_pReso = (TNtuple*)F_EposDisto_p_pReso->Get("InfoTuple_ClosePairs");
    unsigned N_EposDisto_p_pReso = T_EposDisto_p_pReso->GetEntries();
    T_EposDisto_p_pReso->SetBranchAddress("k_D",&k_D);
    T_EposDisto_p_pReso->SetBranchAddress("P1",&fP1);
    T_EposDisto_p_pReso->SetBranchAddress("P2",&fP2);
    T_EposDisto_p_pReso->SetBranchAddress("M1",&fM1);
    T_EposDisto_p_pReso->SetBranchAddress("M2",&fM2);
    T_EposDisto_p_pReso->SetBranchAddress("Tau1",&Tau1);
    T_EposDisto_p_pReso->SetBranchAddress("Tau2",&Tau2);
    T_EposDisto_p_pReso->SetBranchAddress("AngleRcP1",&AngleRcP1);
    T_EposDisto_p_pReso->SetBranchAddress("AngleRcP2",&AngleRcP2);
    T_EposDisto_p_pReso->SetBranchAddress("AngleP1P2",&AngleP1P2);
    for(unsigned uEntry=0; uEntry<N_EposDisto_p_pReso; uEntry++){
        T_EposDisto_p_pReso->GetEntry(uEntry);
        Tau1 = 0;
        Tau2 = Tau_Proton;
        fM2 = Mass_ProtonReso;
        if(k_D>q_CutOff) continue;
        RanVal1 = RanGen.Exponential(fM2/(fP2*Tau2));
        //this is with the correct sign
        MagicSource.AddBGT_PR(RanVal1,cos(AngleRcP2));
        MagicSource.AddBGT_RP(RanVal1,-cos(AngleRcP2));
        MagicSourceWrongSign.AddBGT_PR(RanVal1,-cos(AngleRcP2));
        MagicSourceWrongSign.AddBGT_RP(RanVal1,cos(AngleRcP2));

        h_PR_AngleRcP2_Rsm->Fill(AngleRcP2);
        h_RP_AngleRcP1_Rsm->Fill(Pi-AngleRcP2);

        h_PR_AngleRcP2_RsmWS->Fill(Pi-AngleRcP2);
        h_RP_AngleRcP1_RsmWS->Fill(AngleRcP2);
    }
    delete F_EposDisto_p_pReso;

    TFile* F_EposDisto_pReso_pReso;
    F_EposDisto_pReso_pReso = new TFile(TString::Format("%s/CatsFiles/Source/EposAngularDist/EposDisto_pReso_pReso.root",GetCernBoxDimi()));
    TNtuple* T_EposDisto_pReso_pReso = (TNtuple*)F_EposDisto_pReso_pReso->Get("InfoTuple_ClosePairs");
    unsigned N_EposDisto_pReso_pReso = T_EposDisto_pReso_pReso->GetEntries();
    T_EposDisto_pReso_pReso->SetBranchAddress("k_D",&k_D);
    T_EposDisto_pReso_pReso->SetBranchAddress("P1",&fP1);
    T_EposDisto_pReso_pReso->SetBranchAddress("P2",&fP2);
    T_EposDisto_pReso_pReso->SetBranchAddress("M1",&fM1);
    T_EposDisto_pReso_pReso->SetBranchAddress("M2",&fM2);
    T_EposDisto_pReso_pReso->SetBranchAddress("Tau1",&Tau1);
    T_EposDisto_pReso_pReso->SetBranchAddress("Tau2",&Tau2);
    T_EposDisto_pReso_pReso->SetBranchAddress("AngleRcP1",&AngleRcP1);
    T_EposDisto_pReso_pReso->SetBranchAddress("AngleRcP2",&AngleRcP2);
    T_EposDisto_pReso_pReso->SetBranchAddress("AngleP1P2",&AngleP1P2);
    for(unsigned uEntry=0; uEntry<N_EposDisto_pReso_pReso; uEntry++){
        T_EposDisto_pReso_pReso->GetEntry(uEntry);
        Tau1 = Tau_Proton;
        Tau2 = Tau_Proton;
        fM1 = Mass_ProtonReso;
        fM2 = Mass_ProtonReso;
        if(k_D>q_CutOff) continue;
        RanVal1 = RanGen.Exponential(fM1/(fP1*Tau1));
        RanVal2 = RanGen.Exponential(fM2/(fP2*Tau2));
        MagicSource.AddBGT_RR(RanVal1,cos(AngleRcP1),RanVal2,cos(AngleRcP2),cos(AngleP1P2));
        MagicSourceWrongSign.AddBGT_RR(RanVal1,cos(AngleRcP1),RanVal2,cos(AngleRcP2),cos(AngleP1P2));
        h_RR_AngleRcP1_Rsm->Fill(AngleRcP1);
        h_RR_AngleRcP2_Rsm->Fill(AngleRcP2);
        h_RR_AngleP1P2_Rsm->Fill(AngleP1P2);
    }
    delete F_EposDisto_pReso_pReso;

  }
  else{
    TFile* F_EposDisto_p_LamReso = new TFile(TString::Format("%s/CatsFiles/Source/EposAngularDist/EposDisto_p_LamReso.root",GetCernBoxDimi()));
    TNtuple* T_EposDisto_p_LamReso = (TNtuple*)F_EposDisto_p_LamReso->Get("InfoTuple_ClosePairs");
    unsigned N_EposDisto_p_LamReso = T_EposDisto_p_LamReso->GetEntries();
    T_EposDisto_p_LamReso->SetBranchAddress("k_D",&k_D);
    T_EposDisto_p_LamReso->SetBranchAddress("P1",&fP1);
    T_EposDisto_p_LamReso->SetBranchAddress("P2",&fP2);
    T_EposDisto_p_LamReso->SetBranchAddress("M1",&fM1);
    T_EposDisto_p_LamReso->SetBranchAddress("M2",&fM2);
    T_EposDisto_p_LamReso->SetBranchAddress("Tau1",&Tau1);
    T_EposDisto_p_LamReso->SetBranchAddress("Tau2",&Tau2);
    T_EposDisto_p_LamReso->SetBranchAddress("AngleRcP1",&AngleRcP1);
    T_EposDisto_p_LamReso->SetBranchAddress("AngleRcP2",&AngleRcP2);
    T_EposDisto_p_LamReso->SetBranchAddress("AngleP1P2",&AngleP1P2);
    for(unsigned uEntry=0; uEntry<N_EposDisto_p_LamReso; uEntry++){
        T_EposDisto_p_LamReso->GetEntry(uEntry);
        Tau1 = 0;
        Tau2 = Tau_Lambda;
        fM2 = Mass_LambdaReso;
        if(k_D>q_CutOff) continue;
        RanVal2 = RanGen.Exponential(fM2/(fP2*Tau2));
        MagicSource.AddBGT_PR(RanVal2,cos(AngleRcP2));
        MagicSourceWrongSign.AddBGT_PR(RanVal2,cos(AngleRcP2));
        h_PR_AngleRcP2_Rsm->Fill(AngleRcP2);
        h_PR_AngleRcP2_RsmWS->Fill(AngleRcP2);
    }
    delete F_EposDisto_p_LamReso;

    TFile* F_EposDisto_pReso_Lam = new TFile(TString::Format("%s/CatsFiles/Source/EposAngularDist/EposDisto_pReso_Lam.root",GetCernBoxDimi()));
    TNtuple* T_EposDisto_pReso_Lam = (TNtuple*)F_EposDisto_pReso_Lam->Get("InfoTuple_ClosePairs");
    unsigned N_EposDisto_pReso_Lam = T_EposDisto_pReso_Lam->GetEntries();
    T_EposDisto_pReso_Lam->SetBranchAddress("k_D",&k_D);
    T_EposDisto_pReso_Lam->SetBranchAddress("P1",&fP1);
    T_EposDisto_pReso_Lam->SetBranchAddress("P2",&fP2);
    T_EposDisto_pReso_Lam->SetBranchAddress("M1",&fM1);
    T_EposDisto_pReso_Lam->SetBranchAddress("M2",&fM2);
    T_EposDisto_pReso_Lam->SetBranchAddress("Tau1",&Tau1);
    T_EposDisto_pReso_Lam->SetBranchAddress("Tau2",&Tau2);
    T_EposDisto_pReso_Lam->SetBranchAddress("AngleRcP1",&AngleRcP1);
    T_EposDisto_pReso_Lam->SetBranchAddress("AngleRcP2",&AngleRcP2);
    T_EposDisto_pReso_Lam->SetBranchAddress("AngleP1P2",&AngleP1P2);
    for(unsigned uEntry=0; uEntry<N_EposDisto_pReso_Lam; uEntry++){
        T_EposDisto_pReso_Lam->GetEntry(uEntry);
        Tau1 = Tau_Proton;
        Tau2 = 0;
        fM1 = Mass_ProtonReso;
        if(k_D>q_CutOff) continue;
        RanVal1 = RanGen.Exponential(fM1/(fP1*Tau1));
        MagicSource.AddBGT_RP(RanVal1,cos(AngleRcP1));
        MagicSourceWrongSign.AddBGT_RP(RanVal1,cos(AngleRcP1));
        h_RP_AngleRcP1_Rsm->Fill(AngleRcP1);
        h_RP_AngleRcP1_RsmWS->Fill(AngleRcP1);
    }
    delete F_EposDisto_pReso_Lam;

    TFile* F_EposDisto_pReso_LamReso = new TFile(TString::Format("%s/CatsFiles/Source/EposAngularDist/EposDisto_pReso_LamReso.root",GetCernBoxDimi()));
    TNtuple* T_EposDisto_pReso_LamReso = (TNtuple*)F_EposDisto_pReso_LamReso->Get("InfoTuple_ClosePairs");
    unsigned N_EposDisto_pReso_LamReso = T_EposDisto_pReso_LamReso->GetEntries();
    T_EposDisto_pReso_LamReso->SetBranchAddress("k_D",&k_D);
    T_EposDisto_pReso_LamReso->SetBranchAddress("P1",&fP1);
    T_EposDisto_pReso_LamReso->SetBranchAddress("P2",&fP2);
    T_EposDisto_pReso_LamReso->SetBranchAddress("M1",&fM1);
    T_EposDisto_pReso_LamReso->SetBranchAddress("M2",&fM2);
    T_EposDisto_pReso_LamReso->SetBranchAddress("Tau1",&Tau1);
    T_EposDisto_pReso_LamReso->SetBranchAddress("Tau2",&Tau2);
    T_EposDisto_pReso_LamReso->SetBranchAddress("AngleRcP1",&AngleRcP1);
    T_EposDisto_pReso_LamReso->SetBranchAddress("AngleRcP2",&AngleRcP2);
    T_EposDisto_pReso_LamReso->SetBranchAddress("AngleP1P2",&AngleP1P2);
    for(unsigned uEntry=0; uEntry<N_EposDisto_pReso_LamReso; uEntry++){
        T_EposDisto_pReso_LamReso->GetEntry(uEntry);
        Tau1 = Tau_Proton;
        Tau2 = Tau_Lambda;
        fM1 = Mass_ProtonReso;
        fM2 = Mass_LambdaReso;
        if(k_D>q_CutOff) continue;
        RanVal1 = RanGen.Exponential(fM1/(fP1*Tau1));
        RanVal2 = RanGen.Exponential(fM2/(fP2*Tau2));
        MagicSource.AddBGT_RR(RanVal1,cos(AngleRcP1),RanVal2,cos(AngleRcP2),cos(AngleP1P2));
        MagicSourceWrongSign.AddBGT_RR(RanVal1,cos(AngleRcP1),RanVal2,cos(AngleRcP2),cos(AngleP1P2));
        h_RR_AngleRcP1_Rsm->Fill(AngleRcP1);
        h_RR_AngleRcP2_Rsm->Fill(AngleRcP2);
        h_RR_AngleP1P2_Rsm->Fill(AngleP1P2);
    }
    delete F_EposDisto_pReso_LamReso;
  }

  h_PR_AngleRcP2_Rsm->Sumw2();
  h_RP_AngleRcP1_Rsm->Sumw2();
  h_RR_AngleRcP1_Rsm->Sumw2();
  h_RR_AngleRcP2_Rsm->Sumw2();
  h_RR_AngleP1P2_Rsm->Sumw2();

  h_PR_AngleRcP2_RsmWS->Sumw2();
  h_RP_AngleRcP1_RsmWS->Sumw2();

  h_PR_AngleRcP2_Rsm->Scale(1./h_PR_AngleRcP2_Rsm->Integral(),"width");
  h_RP_AngleRcP1_Rsm->Scale(1./h_RP_AngleRcP1_Rsm->Integral(),"width");
  h_RR_AngleRcP1_Rsm->Scale(1./h_RR_AngleRcP1_Rsm->Integral(),"width");
  h_RR_AngleRcP2_Rsm->Scale(1./h_RR_AngleRcP2_Rsm->Integral(),"width");
  h_RR_AngleP1P2_Rsm->Scale(1./h_RR_AngleP1P2_Rsm->Integral(),"width");

  h_PR_AngleRcP2_RsmWS->Scale(1./h_PR_AngleRcP2_Rsm->Integral(),"width");
  h_RP_AngleRcP1_RsmWS->Scale(1./h_RP_AngleRcP1_Rsm->Integral(),"width");
  ///////////////////////////

  TREPNI Database(0);
  Database.SetSeed(11);
  std::vector<TreParticle*> ParticleList;
  ParticleList.push_back(Database.NewParticle("Proton"));
  ParticleList.push_back(Database.NewParticle("Lambda"));
  ParticleList.push_back(Database.NewParticle("Pion"));

  ParticleList.push_back(Database.NewParticle("ProtonReso"));
  ParticleList.push_back(Database.NewParticle("LambdaReso"));

  double ExpPP,ExpPR,ExpRP,ExpRR;
  ExpPP = (1.-FracProtonReso)*(1.-FracLambdaReso);
  ExpPR = (1.-FracProtonReso)*FracLambdaReso;
  ExpRP = FracProtonReso*(1.-FracLambdaReso);
  ExpRR = FracProtonReso*FracLambdaReso;

  for(TreParticle* prt : ParticleList){
    if(prt->GetName()=="Proton"){
      prt->SetMass(Mass_p);
      prt->SetAbundance((1.-FracProtonReso));
      prt->SetPtEtaPhi(*dlm_pT_eta_p);
    }
    else if(prt->GetName()=="Lambda"){
      prt->SetMass(Mass_L);
      if(Type==1) prt->SetAbundance(0);
      else prt->SetAbundance((1.-FracLambdaReso));
      prt->SetPtEtaPhi(*dlm_pT_eta_L);
    }
    else if(prt->GetName()=="Pion"){
      prt->SetMass(Mass_pic);
      prt->SetAbundance(0);
    }
    else if(prt->GetName()=="ProtonReso"){
      prt->SetMass(Mass_ProtonReso);
      prt->SetAbundance(FracProtonReso);
      prt->SetWidth(hbarc/Tau_Proton);

      prt->NewDecay();
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Proton"));
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
      prt->GetDecay(0)->SetBranching(100);

      prt->SetPtEtaPhi(*dlm_pT_eta_p);
    }
    else if(prt->GetName()=="LambdaReso"){
      prt->SetMass(Mass_LambdaReso);
      if(Type==1) prt->SetAbundance(0);
      else prt->SetAbundance(FracLambdaReso);
      prt->SetWidth(hbarc/Tau_Lambda);

      prt->NewDecay();
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Lambda"));
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
      prt->GetDecay(0)->SetBranching(100);

      prt->SetPtEtaPhi(*dlm_pT_eta_L);
    }

    //prt->SetPtPz(prt->GetMass()*MomSpread*0.001,prt->GetMass()*MomSpread*0.001);
    //prt->SetAcceptance_Eta(-EtaCut,EtaCut);
    //prt->SetAcceptance_pT(500,1e6);

  }
  std::vector<std::string> ListOfParticles;
  ListOfParticles.push_back(part1);
  ListOfParticles.push_back(part2);
  CECA Ivana(Database,ListOfParticles);
  Ivana.SetDisplacementT(DisplT);
  Ivana.SetDisplacementZ(DisplZ);
  Ivana.SetHadronizationT(PancakeT);
  Ivana.SetHadronizationZ(PancakeZ);
  Ivana.SetHadrFluctuation(PancakeFluct*0.01);
  Ivana.SetTau(Tau);
  //Ivana.SetThermalKick(THERMAL_KICK*ThKick);
  Ivana.SetTargetStatistics(10);
  Ivana.SetEventMult(2);
  Ivana.SetSourceDim(2);
  Ivana.SetThreadTimeout(TIMEOUT);
  Ivana.SetFemtoRegion(q_CutOff*0.5);
  Ivana.EqualizeFsiTime(EQUALIZE_TAU);
  Ivana.SetPropagateMother(PropMother);
  Ivana.GHETTO_EVENT = true;
  Ivana.SetUp_RSM = &MagicSourceCeca;
  Ivana.Ghetto_MinMt = 1000.*(MeanMt-MeanMtRange);
  Ivana.Ghetto_MaxMt = 1000.*(MeanMt+MeanMtRange);
  printf("--- CECA source ---\n");
  Ivana.GoBabyGo(THREADS);
  Ivana.SetUp_RSM=NULL;

  double r_fit;
  double r_eff_Rsm;
  double r_eff_RsmCeca;
  double r_eff_RsmWS;
  double r_eff_Ceca;
  double adjust;

  Ivana.Ghetto_kstar->ComputeError();
  Ivana.Ghetto_kstar->ScaleToIntegral();
  Ivana.Ghetto_kstar->ScaleToBinSize();
  TH1F* h_Ghetto_kstar = Convert_DlmHisto_TH1F(Ivana.Ghetto_kstar,"Ghetto_kstar");

  Ivana.GhettoFemto_rstar->ComputeError();
  h_rstar_Ceca = Convert_DlmHisto_TH1F(Ivana.GhettoFemto_rstar,"h_rstar_Ceca");
  r_eff_Ceca = Get_reff(h_rstar_Ceca);

  Ivana.Ghetto_kstar_rstar->ComputeError();
  h_kstar_rstar_Ceca = Convert_DlmHisto_TH2F(Ivana.Ghetto_kstar_rstar,"h_kstar_rstar_Ceca");
  Ivana.Ghetto_rstar->ComputeError();
  h_rstar_CecaAll = Convert_DlmHisto_TH1F(Ivana.Ghetto_rstar,"h_rstar_CecaAll");
  h_rstar_CecaAll->Scale(1./h_rstar_CecaAll->Integral(),"width");

  Ivana.Ghetto_ScatteringAngle->ComputeError();
  TH1F* h_fsi_ScatAngle = Convert_DlmHisto_TH1F(Ivana.Ghetto_ScatteringAngle,"h_fsi_ScatAngle");
  h_fsi_ScatAngle->Scale(1./h_fsi_ScatAngle->Integral(),"width");


  printf(" Final parameters: Disp = %.2f %.2f fm; Hadr = %.2f %.2f fm => reff = %.4f <r*> = %.4f\n",
  Ivana.GetDisplacementT(),Ivana.GetDisplacementZ(),Ivana.GetHadronizationT(),Ivana.GetHadronizationZ(),r_eff_Ceca,h_rstar_Ceca->GetMean());
  h_rstar_Ceca->Scale(1./h_rstar_Ceca->Integral(),"width");


  Ivana.GhettoFemto_rcore->ComputeError();
  TH1F* h_rcore_Ceca_alt = Convert_DlmHisto_TH1F(Ivana.GhettoFemto_rcore,"h_rcore_Ceca_alt");
  h_rcore_Ceca_alt->Scale(1./h_rcore_Ceca_alt->Integral(),"width");
  double rcore_alt = Get_reff(h_rcore_Ceca_alt);

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
  printf("   Expected  Total  Femto\n");
  printf("PP %6.2f%% %6.2f%% %6.2f\n",ExpPP*100.,TotPP*100.,FemtoPP*100.);
  printf("PR %6.2f%% %6.2f%% %6.2f\n",ExpPR*100.,TotPR*100.,FemtoPR*100.);
  printf("RP %6.2f%% %6.2f%% %6.2f\n",ExpRP*100.,TotRP*100.,FemtoRP*100.);
  printf("RR %6.2f%% %6.2f%% %6.2f\n",ExpRR*100.,TotRR*100.,FemtoRR*100.);

  Ivana.Ghetto_RP_AngleRcP1->ComputeError();
  Ivana.Ghetto_PR_AngleRcP2->ComputeError();
  Ivana.Ghetto_RR_AngleRcP1->ComputeError();
  Ivana.Ghetto_RR_AngleRcP2->ComputeError();
  Ivana.Ghetto_RR_AngleP1P2->ComputeError();
  TH1F* h_RP_AngleRcP1_Ceca = Convert_DlmHisto_TH1F(Ivana.Ghetto_RP_AngleRcP1,"h_RP_AngleRcP1_Ceca");
  TH1F* h_PR_AngleRcP2_Ceca = Convert_DlmHisto_TH1F(Ivana.Ghetto_PR_AngleRcP2,"h_PR_AngleRcP2_Ceca");
  TH1F* h_RR_AngleRcP1_Ceca = Convert_DlmHisto_TH1F(Ivana.Ghetto_RR_AngleRcP1,"h_RR_AngleRcP1_Ceca");
  TH1F* h_RR_AngleRcP2_Ceca = Convert_DlmHisto_TH1F(Ivana.Ghetto_RR_AngleRcP2,"h_RR_AngleRcP2_Ceca");
  TH1F* h_RR_AngleP1P2_Ceca = Convert_DlmHisto_TH1F(Ivana.Ghetto_RR_AngleP1P2,"h_RR_AngleP1P2_Ceca");

  h_RP_AngleRcP1_Ceca->Scale(1./h_RP_AngleRcP1_Ceca->Integral(),"width");
  h_PR_AngleRcP2_Ceca->Scale(1./h_PR_AngleRcP2_Ceca->Integral(),"width");
  h_RR_AngleRcP1_Ceca->Scale(1./h_RR_AngleRcP1_Ceca->Integral(),"width");
  h_RR_AngleRcP2_Ceca->Scale(1./h_RR_AngleRcP2_Ceca->Integral(),"width");
  h_RR_AngleP1P2_Ceca->Scale(1./h_RR_AngleP1P2_Ceca->Integral(),"width");

  Ivana.GhettoFemto_mT_kstar->ComputeError();
  TH2F* h_GhettoFemto_mT_kstar = Convert_DlmHisto_TH2F(Ivana.GhettoFemto_mT_kstar,"GhettoFemto_mT_kstar");

  Ivana.GhettoFemto_mT_rstar->ComputeError();
  TH2F* h_GhettoFemto_mT_rstar = Convert_DlmHisto_TH2F(Ivana.GhettoFemto_mT_rstar,"GhettoFemto_mT_rstar");

  TGraphErrors g_GhettoFemto_mT_rstar;
  g_GhettoFemto_mT_rstar.SetName("g_GhettoFemto_mT_rstar");
  g_GhettoFemto_mT_rstar.SetMarkerStyle(20);
  g_GhettoFemto_mT_rstar.SetMarkerSize(2);
  g_GhettoFemto_mT_rstar.SetLineWidth(3);

  unsigned uPoint = 0;
  for(unsigned uBin=0; uBin<h_GhettoFemto_mT_rstar->GetXaxis()->GetNbins(); uBin++){
    TH1F* hProj = (TH1F*)h_GhettoFemto_mT_rstar->ProjectionY(TString::Format("hProj"),uBin+1,uBin+1);
    if(hProj->GetEntries()<200) {delete hProj; continue;}
    hProj->Scale(1./hProj->Integral(),"width");
    double mT = h_GhettoFemto_mT_rstar->GetXaxis()->GetBinCenter(uBin+1);
    g_GhettoFemto_mT_rstar.SetPoint(uPoint,mT*0.001,Get_reff(hProj));
    uPoint++;
    delete hProj;
  }

  TFile fMtEff_pp(TString::Format("%s/SourceStudies/SourcePaper_pp/Systematics/gMt_pp.root",GetFemtoOutputFolder()),"read");
  TGraphErrors* gMtEff_pp = (TGraphErrors*)fMtEff_pp.Get("gMtStatGauss_pp");
  //TGraphErrors* gMtCore_pp = (TGraphErrors*)fMtEff.Get("gMtStatCore_pp");
  TFile fMtEff_pL(TString::Format("%s/SourceStudies/SourcePaper_pL/Systematics/gMt_pL.root",GetFemtoOutputFolder()),"read");
  TGraphErrors* gMtEff_pL = (TGraphErrors*)fMtEff_pL.Get("gMtStatGauss_pL_NLO");
  //TGraphErrors* gMtCore_pp = (TGraphErrors*)fMtEff_pL.Get("gMtStatCore_pL_NLO");

  printf("--- CECA core ---\n");
  Database.GetParticle("ProtonReso")->SetAbundance(0);
  Database.GetParticle("LambdaReso")->SetAbundance(0);
  Ivana.SetThreadTimeout(TIMEOUT);
  Ivana.GoBabyGo(THREADS);
  Ivana.GhettoFemto_rstar->ComputeError();
  h_rcore_Ceca = Convert_DlmHisto_TH1F(Ivana.GhettoFemto_rstar,"h_rcore_Ceca");
  h_rcore_Ceca->Scale(1./h_rcore_Ceca->Integral(),"width");
  double rcore = Get_reff(h_rcore_Ceca);
  printf(" Final r_core = %.4f fm; reff = %.4f fm; <r*> = %.4f\n",rcore,r_eff_Ceca,h_rcore_Ceca->GetMean());
  if(FORCED_RCORE){rcore = FORCED_RCORE; printf("  -- reset rcore to: %.4f\n",rcore);}
  Ivana.Ghetto_rstar->ComputeError();
  h_rcore_CecaAll = Convert_DlmHisto_TH1F(Ivana.Ghetto_rstar,"h_rcore_CecaAll");
  h_rcore_CecaAll->Scale(1./h_rcore_CecaAll->Integral(),"width");

  Ivana.GhettoFemto_mT_rstar->ComputeError();
  TH2F* h_GhettoFemto_mT_core = Convert_DlmHisto_TH2F(Ivana.GhettoFemto_mT_rstar,"h_GhettoFemto_mT_core");

  Ivana.Ghetto_ScatteringAngle->ComputeError();
  TH1F* h_core_ScatAngle = Convert_DlmHisto_TH1F(Ivana.Ghetto_ScatteringAngle,"h_core_ScatAngle");
  h_core_ScatAngle->Scale(1./h_core_ScatAngle->Integral(),"width");



  TGraphErrors g_GhettoFemto_mT_core;
  g_GhettoFemto_mT_core.SetName("g_GhettoFemto_mT_core");
  g_GhettoFemto_mT_core.SetMarkerStyle(20);
  g_GhettoFemto_mT_core.SetMarkerSize(2);
  g_GhettoFemto_mT_core.SetLineWidth(3);

  uPoint = 0;
  for(unsigned uBin=0; uBin<h_GhettoFemto_mT_core->GetXaxis()->GetNbins(); uBin++){
    TH1F* hProj = (TH1F*)h_GhettoFemto_mT_core->ProjectionY(TString::Format("hProj"),uBin+1,uBin+1);
    if(hProj->GetEntries()<200) {delete hProj; continue;}
    hProj->Scale(1./hProj->Integral(),"width");
    double mT = h_GhettoFemto_mT_core->GetXaxis()->GetBinCenter(uBin+1);
    g_GhettoFemto_mT_core.SetPoint(uPoint,mT*0.001,Get_reff(hProj));
    uPoint++;
    delete hProj;
  }

  printf("--- RSM source ---\n");
  printf(" Starting parameter: rcore = %.4f fm\n",rcore);
  for(unsigned uBin=0; uBin<NumRadBins; uBin++){
    double rstar = h_rstar_Rsm->GetBinCenter(uBin+1);
    double parameters[2];
    parameters[0] = rcore;
    parameters[1] = 2.0;

    double val = MagicSource.RootEval(&rstar,parameters);
    h_rstar_Rsm->SetBinContent(uBin+1,val);
    h_rstar_Rsm->SetBinError(uBin+1,1e-3);

    val = MagicSourceWrongSign.RootEval(&rstar,parameters);
    h_rstar_RsmWS->SetBinContent(uBin+1,val);
    h_rstar_RsmWS->SetBinError(uBin+1,1e-3);

    val = MagicSourceCeca.RootEval(&rstar,parameters);
    h_rstar_RsmCeca->SetBinContent(uBin+1,val);
    h_rstar_RsmCeca->SetBinError(uBin+1,1e-3);
  }

  r_eff_Rsm = Get_reff(h_rstar_Rsm);
  r_eff_RsmCeca = Get_reff(h_rstar_RsmCeca);
  r_eff_RsmWS = Get_reff(h_rstar_RsmWS);
  printf(" Final rcore = %.4f fm; reff = %.4f fm; <r*> = %.4f\n",rcore,r_eff_Rsm,h_rstar_Rsm->GetMean());
  printf(" Ceca  rcore = %.4f fm\n",rcore_alt);
  printf("         With wrong sign: reff = %.4f fm; <r*> = %.4f\n",r_eff_RsmWS,h_rstar_RsmWS->GetMean());
  printf("          CECA based RSM: reff = %.4f fm; <r*> = %.4f\n",r_eff_RsmCeca,h_rstar_RsmCeca->GetMean());
  printf("               From CECA: reff = %.4f fm; <r*> = %.4f\n",r_eff_Ceca,h_rstar_Ceca->GetMean());
  h_rstar_Rsm->Scale(1./h_rstar_Rsm->Integral(),"width");
  h_rstar_RsmWS->Scale(1./h_rstar_RsmWS->Integral(),"width");
  h_rstar_RsmCeca->Scale(1./h_rstar_Rsm->Integral(),"width");

  printf("--- RSM core ---\n");
  f_rcore_Rsm->FixParameter(0,1);
  f_rcore_Rsm->FixParameter(1,rcore);
  for(unsigned uBin=0; uBin<NumRadBins; uBin++){
    double rstar = h_rcore_Rsm->GetBinCenter(uBin+1);
    h_rcore_Rsm->SetBinContent(uBin+1,f_rcore_Rsm->Eval(rstar));
    h_rcore_Rsm->SetBinError(uBin+1,1e-3);
  }

  TFile fOutput(TString::Format("%s/FunWithCeca/Ceca_%s_%s_vs2_RSM_EqTau%i_Tau%.0f_Dsp%.0f-%.0f_Pnk%.0f-%.0f_Fct%.0f_PM%i.root",
  GetFemtoOutputFolder(),ListOfParticles.at(0).c_str(),ListOfParticles.at(1).c_str(),
  EQUALIZE_TAU,Tau*100.,Ivana.GetDisplacementT()*100.,Ivana.GetDisplacementZ()*100.,
  Ivana.GetHadronizationT()*100.,Ivana.GetHadronizationZ()*100.,Ivana.GetHadrFluctuation()*100.,PropMother),"recreate");

  fOutput.cd();
  h_reff->Write();
  f_reff->Write();
  h_Ghetto_kstar->Write();
  h_kstar_rstar_Ceca->Write();
  h_rstar_Ceca->Write();
  h_rstar_CecaAll->Write();
  h_rstar_Rsm->Write();
  h_rstar_RsmCeca->Write();
  h_rcore_Ceca->Write();
  h_rcore_CecaAll->Write();
  f_rcore_Ceca->Write();
  h_rcore_Rsm->Write();
  f_rcore_Rsm->Write();

  h_PR_AngleRcP2_Ceca->Write();
  h_PR_AngleRcP2_Rsm->Write();

  h_RP_AngleRcP1_Ceca->Write();
  h_RP_AngleRcP1_Rsm->Write();

  h_RR_AngleRcP1_Ceca->Write();
  h_RR_AngleRcP1_Rsm->Write();

  h_RR_AngleRcP2_Ceca->Write();
  h_RR_AngleRcP2_Rsm->Write();

  h_RR_AngleP1P2_Ceca->Write();
  h_RR_AngleP1P2_Rsm->Write();

  h_fsi_ScatAngle->Write();
  h_core_ScatAngle->Write();

  h_GhettoFemto_mT_rstar->Write();
  g_GhettoFemto_mT_rstar.Write();
  gMtEff_pp->Write();
  gMtEff_pL->Write();
  h_GhettoFemto_mT_core->Write();
  g_GhettoFemto_mT_core.Write();
  //gMtCore_pp->Write();

  h_GhettoFemto_mT_kstar->Write();

  delete h_Ghetto_kstar;
  delete h_reff;
  if(h_rstar_Ceca){delete h_rstar_Ceca; h_rstar_Ceca=NULL;}
  delete h_rstar_Rsm;
  delete h_rstar_RsmWS;
  delete h_rstar_RsmCeca;
  delete h_rcore_Rsm;
  delete h_rcore_Ceca;
  delete h_PR_AngleRcP2_Rsm;
  delete h_RP_AngleRcP1_Rsm;
  delete h_PR_AngleRcP2_RsmWS;
  delete h_RP_AngleRcP1_RsmWS;
  delete h_RR_AngleRcP1_Rsm;
  delete h_RR_AngleRcP2_Rsm;
  delete h_RR_AngleP1P2_Rsm;
  delete h_PR_AngleRcP2_Ceca;
  delete h_RP_AngleRcP1_Ceca;
  delete h_RR_AngleRcP1_Ceca;
  delete h_RR_AngleRcP2_Ceca;
  delete h_RR_AngleP1P2_Ceca;
  delete h_GhettoFemto_mT_rstar;
  delete h_GhettoFemto_mT_kstar;
  delete h_kstar_rstar_Ceca;
  delete h_rstar_CecaAll;
  delete h_rcore_CecaAll;
  delete h_fsi_ScatAngle;
  delete h_core_ScatAngle;
  delete h_GhettoFemto_mT_core;
  delete h_rcore_Ceca_alt;
  delete f_reff;
  delete f_rstar_Ceca;
  delete f_rstar_Rsm;
  delete f_rcore_Ceca;
  delete f_rcore_Rsm;
  delete dlm_pT_eta_p;
  delete dlm_pT_eta_L;
}


void CecaTest_p_pi_1(){

  const double reff = 1.3*1.00;
  const double HadronSize = 0.75;
  const double HadronSlope = 0.2;
  const double EtaCut = 0.8;
  const bool EQUALIZE_TAU = true;
  const double TIMEOUT = 40;
  const double TempPP = 171;//171 or 154
  const double Asym = 6;
  const double PancakeT = 0.0;
  const double PancakeZ = 0.0;
  const double Tau = 0.0;
  const double PancakeFluct = 0;//in % !!!
  const double DisplT = 1.06;
  const double DisplZ = 1.06;
  unsigned THREADS = 0;
  const double MomSpread = 850;
  const bool THERMAL_KICK = false;

  const double Sensitivity = 0.005;

  const double VolumeHI = 5000;
  const double TempHI = 154;
  const double ThKick = TempPP/sqrt(3.);
  const double EngyHI = 208.*5.02e6;
  const double ProtonEngy = 6.5e6;
  const double EngyPP = 2.*ProtonEngy;
  const double EngyDenst_PP_HI = 0.029412*TempPP-3.5294;//ratio
  const double VolumePP = TempPP?EngyPP/EngyHI*VolumeHI/EngyDenst_PP_HI:0;
  const double ProtonGamma = ProtonEngy/938.;
  //const double PancakeT = pow(3.*VolumePP*Asym/4./Pi,1./3.);
  //const double PancakeZ = PancakeT/Asym;

  const double q_CutOff = 200;

  printf("The pp volume is %.2f fm^3\n",VolumePP);
  //printf("Pancake  : %.2f x %.2f x %.2f\n", PancakeT,PancakeT,PancakeZ);
  printf("Displacement: %.2f %.2f %.2f fm\n",DisplT,DisplT,DisplZ);
  printf("Pancake     : %.2f %.2f %.2f fm\n",PancakeT,PancakeT,PancakeZ);

  const double Tau_Proton = 1.65;
  const double Tau_pi = 1.50;

  const double Mass_ProtonReso = 1362;
  const double Mass_piReso = 1180;

  const double FracProtonReso = 0.6422;
  const double FracpiReso = 0.682;

  const unsigned NumRadBins = 256;
  const double rMin = 0;
  const double rMax = 32;

  TH1F* h_reff = new TH1F("h_reff","h_reff",NumRadBins,rMin,rMax);
  //TH1F* h_rstar_Ceca = new TH1F("h_rstar_Ceca","h_rstar_Ceca",NumRadBins,rMin,rMax);
  TH1F* h_rstar_Ceca = NULL;
  TH1F* h_rstar_Rsm = new TH1F("h_rstar_Rsm","h_rstar_Rsm",NumRadBins,rMin,rMax);
  TH1F* h_rstar_RsmCeca = new TH1F("h_rstar_RsmCeca","h_rstar_RsmCeca",NumRadBins,rMin,rMax);
  TH1F* h_rcore_Rsm = new TH1F("h_rcore_Rsm","h_rcore_Rsm",NumRadBins,rMin,rMax);
  TH1F* h_rcore_Ceca;// = new TH1F("h_rcore_Ceca","h_rcore_Ceca",NumRadBins,rMin,rMax);

  TH1F* h_PR_AngleRcP2_Rsm = new TH1F("h_PR_AngleRcP2_Rsm","h_PR_AngleRcP2_Rsm",64,0,Pi);
  TH1F* h_RP_AngleRcP1_Rsm = new TH1F("h_RP_AngleRcP1_Rsm","h_RP_AngleRcP1_Rsm",64,0,Pi);
  TH1F* h_RR_AngleRcP2_Rsm = new TH1F("h_RR_AngleRcP2_Rsm","h_RR_AngleRcP2_Rsm",64,0,Pi);
  TH1F* h_RR_AngleRcP1_Rsm = new TH1F("h_RR_AngleRcP1_Rsm","h_RR_AngleRcP1_Rsm",64,0,Pi);
  TH1F* h_RR_AngleP1P2_Rsm = new TH1F("h_RR_AngleP1P2_Rsm","h_RR_AngleP1P2_Rsm",64,0,Pi);

  TF1* f_reff = new TF1("f_reff","[0]*4.*TMath::Pi()*x*x*pow(4.*TMath::Pi()*[1]*[1],-1.5)*exp(-(x*x)/(4.*[1]*[1]))+1.-[0]",rMin,rMax);
  TF1* f_rstar_Ceca = new TF1("f_rstar_Ceca","[0]*4.*TMath::Pi()*x*x*pow(4.*TMath::Pi()*[1]*[1],-1.5)*exp(-(x*x)/(4.*[1]*[1]))+1.-[0]",rMin,rMax);
  TF1* f_rstar_Rsm = new TF1("f_rstar_Rsm","[0]*4.*TMath::Pi()*x*x*pow(4.*TMath::Pi()*[1]*[1],-1.5)*exp(-(x*x)/(4.*[1]*[1]))+1.-[0]",rMin,rMax);
  TF1* f_rcore_Ceca = new TF1("f_rcore_Ceca","[0]*4.*TMath::Pi()*x*x*pow(4.*TMath::Pi()*[1]*[1],-1.5)*exp(-(x*x)/(4.*[1]*[1]))+1.-[0]",rMin,rMax);
  TF1* f_rcore_Rsm = new TF1("f_rcore_Rsm","[0]*4.*TMath::Pi()*x*x*pow(4.*TMath::Pi()*[1]*[1],-1.5)*exp(-(x*x)/(4.*[1]*[1]))+1.-[0]",rMin,rMax);


  f_reff->FixParameter(0,1);
  f_reff->FixParameter(1,reff);
  for(unsigned uBin=0; uBin<NumRadBins; uBin++){
    double rad = h_reff->GetBinCenter(uBin+1);
    double val = f_reff->Eval(rad);
    double err = 1e-3;
    h_reff->SetBinContent(uBin+1,val);
    h_reff->SetBinError(uBin+1,err);
  }


  //set up the RSM
  DLM_CleverMcLevyResoTM MagicSource;
  MagicSource.InitStability(1,2-1e-6,2+1e-6);
  MagicSource.InitScale(38,0.15,2.0);
  MagicSource.InitRad(257*2,0,64);
  MagicSource.InitType(2);
  MagicSource.SetUpReso(0,FracProtonReso);
  MagicSource.SetUpReso(1,FracpiReso);
  MagicSource.InitNumMcIter(1000000);

  DLM_CleverMcLevyResoTM MagicSourceCeca;
  MagicSourceCeca.InitStability(1,2-1e-6,2+1e-6);
  MagicSourceCeca.InitScale(38,0.15,2.0);
  MagicSourceCeca.InitRad(257*2,0,64);
  MagicSourceCeca.InitType(2);
  MagicSourceCeca.SetUpReso(0,FracProtonReso);
  MagicSourceCeca.SetUpReso(1,FracpiReso);
  MagicSourceCeca.InitNumMcIter(1000000);

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
  DLM_Random RanGen(11);
  double RanVal1;
  double RanVal2;
  double RanVal3;

  gROOT->cd();
  //type treated as 3 digit integer (XYZ), where:
  //X00,X01,X10,X11 are the options, each stating if we have a resonance. E.g. 00 would be two primordials.
  //the first digit can be 0 (EPOS) or 1 (CECA)
  TNtuple* ntPiPi = new TNtuple("nt_p_pi","nt_p_pi","Type:k_D:P1:P2:M1:M2:Tau1:Tau2:AngleRcP1:AngleRcP2:AngleP1P2");
  float nt_array[11];

  TFile* F_EposDisto_pi_LamReso = new TFile(TString::Format("%s/CatsFiles/Source/EposAngularDist/ForMaxRamona_pi_SigReso.root",GetCernBoxDimi()));
  TNtuple* T_EposDisto_pi_LamReso = (TNtuple*)F_EposDisto_pi_LamReso->Get("InfoTuple_ClosePairs");
  unsigned N_EposDisto_pi_LamReso = T_EposDisto_pi_LamReso->GetEntries();
  T_EposDisto_pi_LamReso->SetBranchAddress("k_D",&k_D);
  T_EposDisto_pi_LamReso->SetBranchAddress("P1",&fP1);
  T_EposDisto_pi_LamReso->SetBranchAddress("P2",&fP2);
  T_EposDisto_pi_LamReso->SetBranchAddress("M1",&fM1);
  T_EposDisto_pi_LamReso->SetBranchAddress("M2",&fM2);
  T_EposDisto_pi_LamReso->SetBranchAddress("Tau1",&Tau1);
  T_EposDisto_pi_LamReso->SetBranchAddress("Tau2",&Tau2);
  T_EposDisto_pi_LamReso->SetBranchAddress("AngleRcP1",&AngleRcP1);
  T_EposDisto_pi_LamReso->SetBranchAddress("AngleRcP2",&AngleRcP2);
  T_EposDisto_pi_LamReso->SetBranchAddress("AngleP1P2",&AngleP1P2);
  for(unsigned uEntry=0; uEntry<N_EposDisto_pi_LamReso; uEntry++){
      T_EposDisto_pi_LamReso->GetEntry(uEntry);
      Tau1 = Tau_Proton;
      Tau2 = 0;
      fM1 = Mass_ProtonReso;
      if(k_D>q_CutOff) continue;
      RanVal1 = RanGen.Exponential(fM1/(fP2*Tau1));
      MagicSource.AddBGT_RP(RanVal1,cos(Pi-AngleRcP2));
      //MagicSource.AddBGT_RP(1,0.0);
      //printf("RP %f : %f\n",RanVal1,cos(Pi-AngleRcP2));
      h_RP_AngleRcP1_Rsm->Fill(Pi-AngleRcP2);
      nt_array[0]=10; nt_array[1]=k_D; nt_array[2]=fP2; nt_array[3]=0;
      nt_array[4]=fM1; nt_array[5]=0; nt_array[6]=Tau1; nt_array[7]=0;
      nt_array[8]=Pi-AngleRcP2; nt_array[9]=0; nt_array[10]=0;
      ntPiPi->Fill(nt_array);
  }
  delete F_EposDisto_pi_LamReso;

  TFile* F_EposDisto_piReso_Lam = new TFile(TString::Format("%s/CatsFiles/Source/EposAngularDist/ForMaxRamona_piReso_Sig.root",GetCernBoxDimi()));
  TNtuple* T_EposDisto_piReso_Lam = (TNtuple*)F_EposDisto_piReso_Lam->Get("InfoTuple_ClosePairs");
  unsigned N_EposDisto_piReso_Lam = T_EposDisto_piReso_Lam->GetEntries();
  T_EposDisto_piReso_Lam->SetBranchAddress("k_D",&k_D);
  T_EposDisto_piReso_Lam->SetBranchAddress("P1",&fP1);
  T_EposDisto_piReso_Lam->SetBranchAddress("P2",&fP2);
  T_EposDisto_piReso_Lam->SetBranchAddress("M1",&fM1);
  T_EposDisto_piReso_Lam->SetBranchAddress("M2",&fM2);
  T_EposDisto_piReso_Lam->SetBranchAddress("Tau1",&Tau1);
  T_EposDisto_piReso_Lam->SetBranchAddress("Tau2",&Tau2);
  T_EposDisto_piReso_Lam->SetBranchAddress("AngleRcP1",&AngleRcP1);
  T_EposDisto_piReso_Lam->SetBranchAddress("AngleRcP2",&AngleRcP2);
  T_EposDisto_piReso_Lam->SetBranchAddress("AngleP1P2",&AngleP1P2);
  for(unsigned uEntry=0; uEntry<N_EposDisto_piReso_Lam; uEntry++){
      T_EposDisto_piReso_Lam->GetEntry(uEntry);
      Tau1 = 0;
      Tau2 = Tau_pi;
      fM2 = Mass_pic;
      if(k_D>q_CutOff) continue;
      RanVal2 = RanGen.Exponential(fM2/(fP1*Tau2));
      MagicSource.AddBGT_PR(RanVal2,cos(Pi-AngleRcP1));
      //MagicSource.AddBGT_PR(1,0.0);
      //printf("PR %f : %f\n",RanVal2,cos(Pi-AngleRcP1));
      h_PR_AngleRcP2_Rsm->Fill(Pi-AngleRcP1);
      nt_array[0]=1; nt_array[1]=k_D; nt_array[2]=0; nt_array[3]=fP1;
      nt_array[4]=0; nt_array[5]=fM2; nt_array[6]=0; nt_array[7]=Tau2;
      nt_array[8]=0; nt_array[9]=Pi-AngleRcP1; nt_array[10]=0;
      ntPiPi->Fill(nt_array);
  }
  delete F_EposDisto_piReso_Lam;

  TFile* F_EposDisto_piReso_LamReso = new TFile(TString::Format("%s/CatsFiles/Source/EposAngularDist/ForMaxRamona_piReso_SigReso.root",GetCernBoxDimi()));
  TNtuple* T_EposDisto_piReso_LamReso = (TNtuple*)F_EposDisto_piReso_LamReso->Get("InfoTuple_ClosePairs");
  unsigned N_EposDisto_piReso_LamReso = T_EposDisto_piReso_LamReso->GetEntries();
  T_EposDisto_piReso_LamReso->SetBranchAddress("k_D",&k_D);
  T_EposDisto_piReso_LamReso->SetBranchAddress("P1",&fP1);
  T_EposDisto_piReso_LamReso->SetBranchAddress("P2",&fP2);
  T_EposDisto_piReso_LamReso->SetBranchAddress("M1",&fM1);
  T_EposDisto_piReso_LamReso->SetBranchAddress("M2",&fM2);
  T_EposDisto_piReso_LamReso->SetBranchAddress("Tau1",&Tau1);
  T_EposDisto_piReso_LamReso->SetBranchAddress("Tau2",&Tau2);
  T_EposDisto_piReso_LamReso->SetBranchAddress("AngleRcP1",&AngleRcP1);
  T_EposDisto_piReso_LamReso->SetBranchAddress("AngleRcP2",&AngleRcP2);
  T_EposDisto_piReso_LamReso->SetBranchAddress("AngleP1P2",&AngleP1P2);
  for(unsigned uEntry=0; uEntry<N_EposDisto_piReso_LamReso; uEntry++){
      T_EposDisto_piReso_LamReso->GetEntry(uEntry);
      Tau1 = Tau_Proton;
      Tau2 = Tau_pi;
      fM1 = Mass_ProtonReso;
      fM2 = Mass_pic;
      if(k_D>q_CutOff) continue;
      //RanVal1 = RanGen.Exponential(fM1/(fP2*Tau1));
      //RanVal2 = RanGen.Exponential(fM2/(fP1*Tau2));
      RanVal1 = fP2*Tau1/fM1;
      RanVal2 = fP1*Tau2/fM2;
      //MagicSource.AddBGT_RR(RanVal1,cos(Pi-AngleRcP2),RanVal2,cos(Pi-AngleRcP1),cos(AngleP1P2));
      MagicSource.AddBGT_RR(RanVal1,cos(Pi-AngleRcP2),RanVal2,cos(Pi-AngleRcP1),cos(AngleP1P2));

      //MagicSource.AddBGT_RR(1,0.0,1,0,0);
      //MagicSource.AddBGT_RR(RanVal1,cos(Pi-AngleRcP2),RanVal2,cos(Pi-AngleRcP1),RanGen.Uniform(-1,1));
      //printf("RR %f : %f : %f : %f : %f\n",RanVal1,cos(Pi-AngleRcP2),RanVal2,cos(Pi-AngleRcP1),cos(Pi-AngleP1P2));
      //usleep(200e3);
      h_RR_AngleRcP1_Rsm->Fill(Pi-AngleRcP2);
      h_RR_AngleRcP2_Rsm->Fill(Pi-AngleRcP1);
      h_RR_AngleP1P2_Rsm->Fill(AngleP1P2);
      nt_array[0]=11; nt_array[1]=k_D; nt_array[2]=fP2; nt_array[3]=fP1;
      nt_array[4]=fM1; nt_array[5]=fM2; nt_array[6]=Tau1; nt_array[7]=Tau2;
      nt_array[8]=Pi-AngleRcP2; nt_array[9]=Pi-AngleRcP1; nt_array[10]=AngleP1P2;
      ntPiPi->Fill(nt_array);
  }
  delete F_EposDisto_piReso_LamReso;

  h_PR_AngleRcP2_Rsm->Sumw2();
  h_RP_AngleRcP1_Rsm->Sumw2();
  h_RR_AngleRcP1_Rsm->Sumw2();
  h_RR_AngleRcP2_Rsm->Sumw2();
  h_RR_AngleP1P2_Rsm->Sumw2();

  h_PR_AngleRcP2_Rsm->Scale(1./h_PR_AngleRcP2_Rsm->Integral(),"width");
  h_RP_AngleRcP1_Rsm->Scale(1./h_RP_AngleRcP1_Rsm->Integral(),"width");
  h_RR_AngleRcP1_Rsm->Scale(1./h_RR_AngleRcP1_Rsm->Integral(),"width");
  h_RR_AngleRcP2_Rsm->Scale(1./h_RR_AngleRcP2_Rsm->Integral(),"width");
  h_RR_AngleP1P2_Rsm->Scale(1./h_RR_AngleP1P2_Rsm->Integral(),"width");

  ///////////////////////////

  TREPNI Database(0);
  Database.SetSeed(11);
  std::vector<TreParticle*> ParticleList;
  ParticleList.push_back(Database.NewParticle("Proton"));
  ParticleList.push_back(Database.NewParticle("Pion"));
  ParticleList.push_back(Database.NewParticle("DecayPion"));

  ParticleList.push_back(Database.NewParticle("ProtonReso"));
  ParticleList.push_back(Database.NewParticle("PionReso"));

  double ExpPP,ExpPR,ExpRP,ExpRR;
  ExpPP = (1.-FracProtonReso)*(1.-FracpiReso);
  ExpPR = (1.-FracProtonReso)*FracpiReso;
  ExpRP = FracProtonReso*(1.-FracpiReso);
  ExpRR = FracProtonReso*FracpiReso;

  for(TreParticle* prt : ParticleList){
    if(prt->GetName()=="Proton"){
      prt->SetMass(Mass_p);
      prt->SetAbundance((1.-FracProtonReso));
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);
    }
    else if(prt->GetName()=="Pion"){
      prt->SetMass(Mass_pic);
      prt->SetAbundance((1.-FracpiReso));
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);
    }
    else if(prt->GetName()=="DecayPion"){
      prt->SetMass(Mass_pic);
      prt->SetAbundance(0);
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);
    }
    else if(prt->GetName()=="ProtonReso"){
      prt->SetMass(Mass_ProtonReso);
      prt->SetAbundance(FracProtonReso);
      prt->SetWidth(hbarc/Tau_Proton);
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);

      prt->NewDecay();
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Proton"));
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("DecayPion"));
      prt->GetDecay(0)->SetBranching(100);
    }
    else if(prt->GetName()=="PionReso"){
      prt->SetMass(Mass_piReso);
      prt->SetAbundance(FracpiReso);
      prt->SetWidth(hbarc/Tau_pi);
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);

      prt->NewDecay();
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
      prt->GetDecay(0)->SetBranching(100);
    }

    //prt->SetPtEtaPhi(MomSpread*);
    //prt->SetPtEtaPhi(prt->GetMass()*MomSpread*0.001);
    prt->SetPtPz(prt->GetMass()*MomSpread*0.001,prt->GetMass()*MomSpread*0.001);
  }
  std::vector<std::string> ListOfParticles;
  ListOfParticles.push_back("Proton");
  ListOfParticles.push_back("Pion");
  CECA Ivana(Database,ListOfParticles);
  Ivana.SetDisplacementT(DisplT);
  Ivana.SetDisplacementZ(DisplZ);
  Ivana.SetHadronizationT(PancakeT);
  Ivana.SetHadronizationZ(PancakeZ);
  Ivana.SetHadrFluctuation(PancakeFluct);
  Ivana.SetTau(Tau);
  Ivana.SetThermalKick(THERMAL_KICK*ThKick);
  Ivana.SetTargetStatistics(10);
  Ivana.SetEventMult(2);
  Ivana.SetSourceDim(2);
  Ivana.SetThreadTimeout(TIMEOUT);
  Ivana.SetFemtoRegion(q_CutOff*0.5);
  Ivana.GHETTO_EVENT = true;
  printf("Goal for reff = %.4f fm; <r*> = %.4f fm\n",reff,h_reff->GetMean());
  printf("--- CECA source ---\n");
  Ivana.GoBabyGo(THREADS);

  double r_fit;
  double adjust;

  Ivana.GhettoFemto_rstar->ComputeError();
  h_rstar_Ceca = Convert_DlmHisto_TH1F(Ivana.GhettoFemto_rstar,"h_rstar_Ceca");
  r_fit = Get_reff(h_rstar_Ceca);
  printf(" Starting parameters: Disp = %.2f %.2f fm; Hadr = %.2f %.2f fm\n",
  Ivana.GetDisplacementT(),Ivana.GetDisplacementZ(),Ivana.GetHadronizationT(),Ivana.GetHadronizationZ());
  double DEV = fabs(r_fit-reff)/fabs(reff);
  while(DEV>Sensitivity){
    adjust = reff/r_fit;
    Ivana.SetDisplacementT(Ivana.GetDisplacementT()*adjust);
    Ivana.SetDisplacementZ(Ivana.GetDisplacementZ()*adjust);
    Ivana.SetHadronizationT(Ivana.GetHadronizationT()*adjust);
    Ivana.SetHadronizationZ(Ivana.GetHadronizationZ()*adjust);
    Ivana.SetTau(Ivana.GetTau()*adjust);
    printf(" Redo with (DEV=%.2f%%): Disp = %.2f %.2f; Hadr = %.2f %.2f; Tau = %.2f\n",
    DEV*100.,Ivana.GetDisplacementT(),Ivana.GetDisplacementZ(),Ivana.GetHadronizationT(),Ivana.GetHadronizationZ(), Ivana.GetTau());
    Ivana.GoBabyGo(THREADS);
    delete h_rstar_Ceca;
    Ivana.GhettoFemto_rstar->ComputeError();
    h_rstar_Ceca = Convert_DlmHisto_TH1F(Ivana.GhettoFemto_rstar,"h_rstar_Ceca");
    r_fit = Get_reff(h_rstar_Ceca);
    DEV = fabs(r_fit-reff)/fabs(reff);
  }

  Ivana.SetUp_RSM = &MagicSourceCeca;
  std::vector<float*> buffer_ntPiPi;
  Ivana.Buffer_RSM = &buffer_ntPiPi;
  Ivana.GoBabyGo(THREADS);
  Ivana.SetUp_RSM = NULL;
  Ivana.Buffer_RSM = NULL;
  for(float* entry : buffer_ntPiPi){
    ntPiPi->Fill(entry);
    delete [] entry;
    entry = NULL;
  }

  printf(" Final parameters: Disp = %.2f %.2f fm; Hadr = %.2f %.2f fm => reff = %.4f <r*> = %.4f\n",
  Ivana.GetDisplacementT(),Ivana.GetDisplacementZ(),Ivana.GetHadronizationT(),Ivana.GetHadronizationZ(),r_fit,h_rstar_Ceca->GetMean());
  h_rstar_Ceca->Scale(1./h_rstar_Ceca->Integral(),"width");

  Ivana.Ghetto_RP_AngleRcP1->ComputeError();
  Ivana.Ghetto_PR_AngleRcP2->ComputeError();
  Ivana.Ghetto_RR_AngleRcP1->ComputeError();
  Ivana.Ghetto_RR_AngleRcP2->ComputeError();
  Ivana.Ghetto_RR_AngleP1P2->ComputeError();
  TH1F* h_RP_AngleRcP1_Ceca = Convert_DlmHisto_TH1F(Ivana.Ghetto_RP_AngleRcP1,"h_RP_AngleRcP1_Ceca");
  TH1F* h_PR_AngleRcP2_Ceca = Convert_DlmHisto_TH1F(Ivana.Ghetto_PR_AngleRcP2,"h_PR_AngleRcP2_Ceca");
//for(unsigned uBin=0; uBin<h_PR_AngleRcP2_Ceca->GetNbinsX(); uBin++){
//  printf("%5.2e vs %5.2e\n",h_PR_AngleRcP2_Ceca->GetBinContent(uBin+1),Ivana.Ghetto_PR_AngleRcP2->GetBinContent(uBin));
//}
  TH1F* h_RR_AngleRcP1_Ceca = Convert_DlmHisto_TH1F(Ivana.Ghetto_RR_AngleRcP1,"h_RR_AngleRcP1_Ceca");
  TH1F* h_RR_AngleRcP2_Ceca = Convert_DlmHisto_TH1F(Ivana.Ghetto_RR_AngleRcP2,"h_RR_AngleRcP2_Ceca");
  TH1F* h_RR_AngleP1P2_Ceca = Convert_DlmHisto_TH1F(Ivana.Ghetto_RR_AngleP1P2,"h_RR_AngleP1P2_Ceca");

  h_RP_AngleRcP1_Ceca->Scale(1./h_RP_AngleRcP1_Ceca->Integral(),"width");
  h_PR_AngleRcP2_Ceca->Scale(1./h_PR_AngleRcP2_Ceca->Integral(),"width");
  h_RR_AngleRcP1_Ceca->Scale(1./h_RR_AngleRcP1_Ceca->Integral(),"width");
  h_RR_AngleRcP2_Ceca->Scale(1./h_RR_AngleRcP2_Ceca->Integral(),"width");
  h_RR_AngleP1P2_Ceca->Scale(1./h_RR_AngleP1P2_Ceca->Integral(),"width");

  Ivana.GhettoFemto_mT_kstar->ComputeError();
  TH2F* h_GhettoFemto_mT_kstar = Convert_DlmHisto_TH2F(Ivana.GhettoFemto_mT_kstar,"GhettoFemto_mT_kstar");

  Ivana.GhettoFemto_mT_rstar->ComputeError();
  TH2F* h_GhettoFemto_mT_rstar = Convert_DlmHisto_TH2F(Ivana.GhettoFemto_mT_rstar,"GhettoFemto_mT_rstar");

  TGraphErrors g_GhettoFemto_mT_rstar;
  g_GhettoFemto_mT_rstar.SetName("g_GhettoFemto_mT_rstar");
  g_GhettoFemto_mT_rstar.SetMarkerStyle(20);
  g_GhettoFemto_mT_rstar.SetMarkerSize(2);
  g_GhettoFemto_mT_rstar.SetLineWidth(3);

  unsigned uPoint = 0;
  for(unsigned uBin=0; uBin<h_GhettoFemto_mT_rstar->GetXaxis()->GetNbins(); uBin++){
    TH1F* hProj = (TH1F*)h_GhettoFemto_mT_rstar->ProjectionY(TString::Format("hProj"),uBin+1,uBin+1);
    if(hProj->GetEntries()<200) {delete hProj; continue;}
    hProj->Scale(1./hProj->Integral(),"width");
    double mT = h_GhettoFemto_mT_rstar->GetXaxis()->GetBinCenter(uBin+1);
    g_GhettoFemto_mT_rstar.SetPoint(uPoint,mT*0.001,Get_reff(hProj));
    uPoint++;
    delete hProj;
  }

  //TFile fMtEff(TString::Format("%s/SourceStudies/SourcePaper_pp/Systematics/gMt_pp.root",GetFemtoOutputFolder()),"read");
  //TGraphErrors* gMtEff_pp = (TGraphErrors*)fMtEff.Get("gMtStatGauss_pp");
  //TGraphErrors* gMtCore_pp = (TGraphErrors*)fMtEff.Get("gMtStatCore_pp");
  TFile fMtEff(TString::Format("%s/SourceStudies/SourcePaper_pL/Systematics/gMt_pL.root",GetFemtoOutputFolder()),"read");
  TGraphErrors* gMtEff_pp = (TGraphErrors*)fMtEff.Get("gMtStatGauss_pL_NLO");
  TGraphErrors* gMtCore_pp = (TGraphErrors*)fMtEff.Get("gMtStatCore_pL_NLO");

  printf("--- CECA core ---\n");
  Database.GetParticle("ProtonReso")->SetAbundance(0);
  Database.GetParticle("PionReso")->SetAbundance(0);
  Ivana.GoBabyGo(THREADS);
  Ivana.GhettoFemto_rstar->ComputeError();
  h_rcore_Ceca = Convert_DlmHisto_TH1F(Ivana.GhettoFemto_rstar,"h_rcore_Ceca");
  h_rcore_Ceca->Scale(1./h_rcore_Ceca->Integral(),"width");
  printf(" Final r_core = %.4f fm; reff = %.4f fm; <r*> = %.4f\n",Get_reff(h_rcore_Ceca),r_fit,h_rcore_Ceca->GetMean());

  Ivana.GhettoFemto_mT_rstar->ComputeError();
  TH2F* h_GhettoFemto_mT_core = Convert_DlmHisto_TH2F(Ivana.GhettoFemto_mT_rstar,"h_GhettoFemto_mT_core");

  TGraphErrors g_GhettoFemto_mT_core;
  g_GhettoFemto_mT_core.SetName("g_GhettoFemto_mT_core");
  g_GhettoFemto_mT_core.SetMarkerStyle(20);
  g_GhettoFemto_mT_core.SetMarkerSize(2);
  g_GhettoFemto_mT_core.SetLineWidth(3);

  uPoint = 0;
  for(unsigned uBin=0; uBin<h_GhettoFemto_mT_core->GetXaxis()->GetNbins(); uBin++){
    TH1F* hProj = (TH1F*)h_GhettoFemto_mT_core->ProjectionY(TString::Format("hProj"),uBin+1,uBin+1);
    if(hProj->GetEntries()<200) {delete hProj; continue;}
    hProj->Scale(1./hProj->Integral(),"width");
    double mT = h_GhettoFemto_mT_core->GetXaxis()->GetBinCenter(uBin+1);
    g_GhettoFemto_mT_core.SetPoint(uPoint,mT*0.001,Get_reff(hProj));
    uPoint++;
    delete hProj;
  }


  printf("--- RSM source ---\n");

  double rcore = reff*0.819;
  //double rcore = reff*1.46;
  printf(" Starting parameter: rcore = %.4f fm\n",rcore);
  for(unsigned uBin=0; uBin<NumRadBins; uBin++){
    double rstar = h_rstar_Rsm->GetBinCenter(uBin+1);
    double parameters[2];
    parameters[0] = rcore;
    parameters[1] = 2.0;
    //printf("rstar = %f; rcore = %f\n",rstar,rcore);
    double val = MagicSource.RootEval(&rstar,parameters);

    h_rstar_Rsm->SetBinContent(uBin+1,val);
    h_rstar_Rsm->SetBinError(uBin+1,1e-3);

    val = MagicSourceCeca.RootEval(&rstar,parameters);
    h_rstar_RsmCeca->SetBinContent(uBin+1,val);
    h_rstar_RsmCeca->SetBinError(uBin+1,1e-3);
  }

  r_fit = Get_reff(h_rstar_Rsm);
  DEV = fabs(r_fit-reff)/fabs(reff);
  /*
  while(DEV>Sensitivity){
    rcore *= reff/r_fit;
    printf(" Redo with (DEV=%.2f%%): rcore = %.3f fm\n",DEV*100.,rcore);
    for(unsigned uBin=0; uBin<NumRadBins; uBin++){
      double rstar = h_rstar_Rsm->GetBinCenter(uBin+1);
      double parameters[2];
      parameters[0] = rcore;
      parameters[1] = 2.0;
      double val = MagicSource.RootEval(&rstar,parameters);
      h_rstar_Rsm->SetBinContent(uBin+1,val);
      h_rstar_Rsm->SetBinError(uBin+1,1e-3);
      //printf(" vale(%.4f) = %.4f\n",rstar,val);

      val = MagicSourceCeca.RootEval(&rstar,parameters);
      //printf(" valc(%.4f) = %.4f\n",rstar,val);
      h_rstar_RsmCeca->SetBinContent(uBin+1,val);
      h_rstar_RsmCeca->SetBinError(uBin+1,1e-3);
    }
    r_fit = Get_reff(h_rstar_Rsm);
    DEV = fabs(r_fit-reff)/fabs(reff);
  }
  */

  printf(" Final rcore = %.4f fm; reff = %.4f fm; <r*> = %.4f\n",rcore,r_fit,h_rstar_Rsm->GetMean());
  printf("   N.B. reff Rsm vs RsmCeca: %.3f fm %.3f fm\n",r_fit,Get_reff(h_rstar_RsmCeca));
  h_rstar_Rsm->Scale(1./h_rstar_Rsm->Integral(),"width");
  h_rstar_RsmCeca->Scale(1./h_rstar_RsmCeca->Integral(),"width");

  printf("--- RSM core ---\n");
  f_rcore_Rsm->FixParameter(0,1);
  f_rcore_Rsm->FixParameter(1,rcore);
  for(unsigned uBin=0; uBin<NumRadBins; uBin++){
    double rstar = h_rcore_Rsm->GetBinCenter(uBin+1);
    h_rcore_Rsm->SetBinContent(uBin+1,f_rcore_Rsm->Eval(rstar));
    h_rcore_Rsm->SetBinError(uBin+1,1e-3);
  }

  TFile fOutput(TString::Format("%s/FunWithCeca/CecaTest_p_pi_EqTau%i_Tau%.0f_Dsp%.0f-%.0f_Pnk%.0f-%.0f_Fct%.0f_Kick%.0f_Hds%.0f-%.0f.root",
  GetFemtoOutputFolder(),EQUALIZE_TAU,Tau*100.,Ivana.GetDisplacementT()*100.,Ivana.GetDisplacementZ()*100.,
  Ivana.GetHadronizationT()*100.,Ivana.GetHadronizationZ()*100.,Ivana.GetHadrFluctuation()*100.,
  ThKick*THERMAL_KICK,HadronSize*.100,HadronSlope*100.),"recreate");

  fOutput.cd();
  h_reff->Write();
  f_reff->Write();
  h_rstar_Ceca->Write();
  h_rstar_Rsm->Write();
  h_rstar_RsmCeca->Write();
  h_rcore_Ceca->Write();
  f_rcore_Ceca->Write();
  h_rcore_Rsm->Write();
  f_rcore_Rsm->Write();

  h_PR_AngleRcP2_Ceca->Write();
  h_PR_AngleRcP2_Rsm->Write();

  h_RP_AngleRcP1_Ceca->Write();
  h_RP_AngleRcP1_Rsm->Write();

  h_RR_AngleRcP1_Ceca->Write();
  h_RR_AngleRcP1_Rsm->Write();

  h_RR_AngleRcP2_Ceca->Write();
  h_RR_AngleRcP2_Rsm->Write();

  h_RR_AngleP1P2_Ceca->Write();
  h_RR_AngleP1P2_Rsm->Write();

  h_GhettoFemto_mT_rstar->Write();
  g_GhettoFemto_mT_rstar.Write();
  gMtEff_pp->Write();
  h_GhettoFemto_mT_core->Write();
  g_GhettoFemto_mT_core.Write();
  gMtCore_pp->Write();

  h_GhettoFemto_mT_kstar->Write();
  ntPiPi->Write();

  delete h_reff;
  if(h_rstar_Ceca){delete h_rstar_Ceca; h_rstar_Ceca=NULL;}
  delete h_rstar_Rsm;
  delete h_rstar_RsmCeca;
  delete h_rcore_Rsm;
  delete h_rcore_Ceca;
  delete h_PR_AngleRcP2_Rsm;
  delete h_RP_AngleRcP1_Rsm;
  delete h_RR_AngleRcP1_Rsm;
  delete h_RR_AngleRcP2_Rsm;
  delete h_RR_AngleP1P2_Rsm;
  delete h_PR_AngleRcP2_Ceca;
  delete h_RP_AngleRcP1_Ceca;
  delete h_RR_AngleRcP1_Ceca;
  delete h_RR_AngleRcP2_Ceca;
  delete h_RR_AngleP1P2_Ceca;
  delete h_GhettoFemto_mT_rstar;
  delete h_GhettoFemto_mT_kstar;
  delete f_reff;
  delete f_rstar_Ceca;
  delete f_rstar_Rsm;
  delete f_rcore_Ceca;
  delete f_rcore_Rsm;
  delete ntPiPi;
}





void Ceca_pp_pL_source(){

  const double EtaCut = 0.8;
  const double PtCut = 500;
  const bool EQUALIZE_TAU = true;
  const double TIMEOUT = 20;
  const double FracProtonReso = 0.6422*1;
  const double FracLambdaReso = 0.6438*1;
  const double DisplT = 1.2;
  const double DisplZ = 1.2;
  const double PancakeT = 0.0;
  const double PancakeZ = 0;
  const double Tau = 0.0;
  const double TauFluct = 0;
  const double PancakeFluct = 0;//in % !!!
  unsigned THREADS = 0;
  const double MomSpread = 1061;
  const double THERMAL_KICK = 0;
  const double q_CutOff = 200;

  double PancakeV = PancakeT*PancakeT*PancakeZ*4./3.*Pi;
  printf("Displacement: %.2f %.2f %.2f fm\n",DisplT,DisplT,DisplZ);
  printf("Pancake     : %.2f %.2f %.2f fm\n",PancakeT,PancakeT,PancakeZ);
  printf("PancakeV    : %.2f fm^3\n",PancakeV);

  const double Tau_Proton = 1.65;
  const double Tau_Lambda = 4.69;//4.69

  const double Mass_ProtonReso = 1362;
  const double Mass_LambdaReso = 1462;

  const unsigned NumRadBins = 256;
  const double rMin = 0;
  const double rMax = 32;

  DLM_Histo<float>* dlm_pT_eta_p = GetPtEtaDist(
    TString::Format("%s/Jaime/p_pT.root",GetCernBoxDimi()),
    TString::Format("%s/Jaime/ap_pT.root",GetCernBoxDimi()));

  DLM_Histo<float>* dlm_pT_eta_L = GetPtEtaDist(
    TString::Format("%s/Jaime/L_pT.root",GetCernBoxDimi()),
    TString::Format("%s/Jaime/aL_pT.root",GetCernBoxDimi()));

  TREPNI Database(0);
  Database.SetSeed(11);
  std::vector<TreParticle*> ParticleList;
  ParticleList.push_back(Database.NewParticle("Proton"));
  ParticleList.push_back(Database.NewParticle("Lambda"));
  ParticleList.push_back(Database.NewParticle("Pion"));

  ParticleList.push_back(Database.NewParticle("ProtonReso"));
  ParticleList.push_back(Database.NewParticle("LambdaReso"));

  double ExpPP,ExpPR,ExpRP,ExpRR;
  ExpPP = (1.-FracProtonReso)*(1.-FracLambdaReso);
  ExpPR = (1.-FracProtonReso)*FracLambdaReso;
  ExpRP = FracProtonReso*(1.-FracLambdaReso);
  ExpRR = FracProtonReso*FracLambdaReso;

  for(TreParticle* prt : ParticleList){
    if(prt->GetName()=="Proton"){
      prt->SetMass(Mass_p);
      prt->SetAbundance((1.-FracProtonReso));
      prt->SetPtEtaPhi(*dlm_pT_eta_p);
    }
    else if(prt->GetName()=="Lambda"){
      prt->SetMass(Mass_L);
      prt->SetAbundance((1.-FracLambdaReso));
      prt->SetPtEtaPhi(*dlm_pT_eta_L);
    }
    else if(prt->GetName()=="Pion"){
      prt->SetMass(Mass_pic);
      prt->SetAbundance(0);
    }
    else if(prt->GetName()=="ProtonReso"){
      prt->SetMass(Mass_ProtonReso);
      prt->SetAbundance(FracProtonReso);
      prt->SetWidth(hbarc/Tau_Proton);

      prt->NewDecay();
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Proton"));
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
      prt->GetDecay(0)->SetBranching(100);

      prt->SetPtEtaPhi(*dlm_pT_eta_p);
    }
    else if(prt->GetName()=="LambdaReso"){
      prt->SetMass(Mass_LambdaReso);
      prt->SetAbundance(FracLambdaReso);
      prt->SetWidth(hbarc/Tau_Lambda);

      prt->NewDecay();
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Lambda"));
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
      prt->GetDecay(0)->SetBranching(100);

      prt->SetPtEtaPhi(*dlm_pT_eta_L);
    }

    //prt->SetPtPz(prt->GetMass()*MomSpread*0.001,prt->GetMass()*MomSpread*0.001);
    prt->SetAcceptance_Eta(-EtaCut,EtaCut);
    prt->SetAcceptance_pT(PtCut,1e6);

  }


  DLM_CleverMcLevyResoTM MagicSource_Uni_pp;
  DLM_CleverMcLevyResoTM MagicSource_EposWS_pp;
  DLM_CleverMcLevyResoTM MagicSource_Epos_pp;
  DLM_CleverMcLevyResoTM MagicSource_BB_pp;
  DLM_CleverMcLevyResoTM MagicSource_Ceca_pp;

  DLM_CleverMcLevyResoTM MagicSource_Uni_pL;
  DLM_CleverMcLevyResoTM MagicSource_Epos_pL;
  DLM_CleverMcLevyResoTM& MagicSource_EposWS_pL = MagicSource_Epos_pL;
  DLM_CleverMcLevyResoTM MagicSource_BB_pL;
  DLM_CleverMcLevyResoTM MagicSource_Ceca_pL;

  BasicSetUp_MS(MagicSource_Uni_pp,FracProtonReso,FracProtonReso);
  BasicSetUp_MS(MagicSource_EposWS_pp,FracProtonReso,FracProtonReso);
  BasicSetUp_MS(MagicSource_Epos_pp,FracProtonReso,FracProtonReso);
  BasicSetUp_MS(MagicSource_BB_pp,FracProtonReso,FracProtonReso);
  BasicSetUp_MS(MagicSource_Ceca_pp,FracProtonReso,FracProtonReso);

  BasicSetUp_MS(MagicSource_Uni_pL,FracProtonReso,FracLambdaReso);
  BasicSetUp_MS(MagicSource_Epos_pL,FracProtonReso,FracLambdaReso);
  BasicSetUp_MS(MagicSource_BB_pL,FracProtonReso,FracLambdaReso);
  BasicSetUp_MS(MagicSource_Ceca_pL,FracProtonReso,FracLambdaReso);


  std::vector<std::string> ListOfParticles;


  printf("--- CECA pp ---\n");
  ListOfParticles.clear();
  ListOfParticles.push_back("Proton");
  ListOfParticles.push_back("Proton");
  CECA Ivana_pp(Database,ListOfParticles);
  Ivana_pp.SetDisplacementT(DisplT);
  Ivana_pp.SetDisplacementZ(DisplZ);
  Ivana_pp.SetHadronizationT(PancakeT);
  Ivana_pp.SetHadronizationZ(PancakeZ);
  Ivana_pp.SetHadrFluctuation(PancakeFluct*0.01);
  Ivana_pp.SetTau(Tau,TauFluct*0.01);
  Ivana_pp.SetTargetStatistics(10);
  Ivana_pp.SetEventMult(2);
  Ivana_pp.SetSourceDim(2);
  Ivana_pp.SetThreadTimeout(TIMEOUT);
  Ivana_pp.SetFemtoRegion(q_CutOff*0.5);
  Ivana_pp.EqualizeFsiTime(EQUALIZE_TAU);
  Ivana_pp.SetThermalKick(THERMAL_KICK);
  //Ivana_pp.SetPropagateMother(PropMother);
  Ivana_pp.GHETTO_EVENT = true;
  Database.GetParticle("Proton")->SetAbundance(1.-FracProtonReso);
  Database.GetParticle("Lambda")->SetAbundance(0);
  Database.GetParticle("ProtonReso")->SetAbundance(FracProtonReso);
  Database.GetParticle("LambdaReso")->SetAbundance(0);
  Ivana_pp.SetUp_RSM = &MagicSource_Ceca_pp;
  Ivana_pp.SetUp_RSM_UNI = &MagicSource_Uni_pp;
  Ivana_pp.SetUp_RSM_BB = &MagicSource_BB_pp;
  Ivana_pp.GoBabyGo(THREADS);
  Ivana_pp.GhettoFemto_mT_rstar->ComputeError();
  TH2F* h_GhettoFemto_mT_rstar_pp = Convert_DlmHisto_TH2F(Ivana_pp.GhettoFemto_mT_rstar,"h_GhettoFemto_mT_rstar_pp");
  Ivana_pp.GhettoFemto_mT_rcore->ComputeError();
  TH2F* h_GhettoFemto_mT_rcore_pp = Convert_DlmHisto_TH2F(Ivana_pp.GhettoFemto_mT_rcore,"h_GhettoFemto_mT_rcore_pp");
  TGraphErrors pp_GhettoFemto_mT_rstar;
  pp_GhettoFemto_mT_rstar.SetName("pp_GhettoFemto_mT_rstar");
  pp_GhettoFemto_mT_rstar.SetMarkerStyle(20);
  pp_GhettoFemto_mT_rstar.SetMarkerSize(2);
  pp_GhettoFemto_mT_rstar.SetLineWidth(3);
  unsigned uPoint = 0;
  for(unsigned uBin=0; uBin<h_GhettoFemto_mT_rstar_pp->GetXaxis()->GetNbins(); uBin++){
    TH1F* hProj = (TH1F*)h_GhettoFemto_mT_rstar_pp->ProjectionY(TString::Format("hProj"),uBin+1,uBin+1);
    if(hProj->GetEntries()<200) {delete hProj; continue;}
    hProj->Scale(1./hProj->Integral(),"width");
    double mT = h_GhettoFemto_mT_rstar_pp->GetXaxis()->GetBinCenter(uBin+1);
    pp_GhettoFemto_mT_rstar.SetPoint(uPoint,mT*0.001,Get_reff(hProj));
    uPoint++;
    delete hProj;
  }
  TGraphErrors pp_GhettoFemto_mT_rcore;
  pp_GhettoFemto_mT_rcore.SetName("pp_GhettoFemto_mT_rcore");
  pp_GhettoFemto_mT_rcore.SetMarkerStyle(20);
  pp_GhettoFemto_mT_rcore.SetMarkerSize(2);
  pp_GhettoFemto_mT_rcore.SetLineWidth(3);
  uPoint = 0;
  for(unsigned uBin=0; uBin<h_GhettoFemto_mT_rcore_pp->GetXaxis()->GetNbins(); uBin++){
    TH1F* hProj = (TH1F*)h_GhettoFemto_mT_rcore_pp->ProjectionY(TString::Format("hProj"),uBin+1,uBin+1);
    if(hProj->GetEntries()<200) {delete hProj; continue;}
    hProj->Scale(1./hProj->Integral(),"width");
    double mT = h_GhettoFemto_mT_rcore_pp->GetXaxis()->GetBinCenter(uBin+1);
    pp_GhettoFemto_mT_rcore.SetPoint(uPoint,mT*0.001,Get_reff(hProj));
    uPoint++;
    delete hProj;
  }


  printf("--- CECA pL ---\n");
  ListOfParticles.clear();
  ListOfParticles.push_back("Proton");
  ListOfParticles.push_back("Lambda");
  CECA Ivana_pL(Database,ListOfParticles);
  Ivana_pL.SetDisplacementT(DisplT);
  Ivana_pL.SetDisplacementZ(DisplZ);
  Ivana_pL.SetHadronizationT(PancakeT);
  Ivana_pL.SetHadronizationZ(PancakeZ);
  Ivana_pL.SetHadrFluctuation(PancakeFluct*0.01);
  Ivana_pL.SetTau(Tau);
  Ivana_pL.SetTargetStatistics(10);
  Ivana_pL.SetEventMult(2);
  Ivana_pL.SetSourceDim(2);
  Ivana_pL.SetThreadTimeout(TIMEOUT);
  Ivana_pL.SetFemtoRegion(q_CutOff*0.5);
  Ivana_pL.EqualizeFsiTime(EQUALIZE_TAU);
  Ivana_pL.SetThermalKick(THERMAL_KICK);
  //Ivana_pL.SetPropagateMother(PropMother);
  Ivana_pL.GHETTO_EVENT = true;
  Database.GetParticle("Proton")->SetAbundance(1.-FracProtonReso);
  Database.GetParticle("Lambda")->SetAbundance(1.-FracLambdaReso);
  Database.GetParticle("ProtonReso")->SetAbundance(FracProtonReso);
  Database.GetParticle("LambdaReso")->SetAbundance(FracLambdaReso);
  Ivana_pL.SetUp_RSM = &MagicSource_Ceca_pL;
  Ivana_pL.SetUp_RSM_UNI = &MagicSource_Uni_pL;
  Ivana_pL.SetUp_RSM_BB = &MagicSource_BB_pL;
  Ivana_pL.GoBabyGo(THREADS);
  Ivana_pL.GhettoFemto_mT_rstar->ComputeError();
  TH2F* h_GhettoFemto_mT_rstar_pL = Convert_DlmHisto_TH2F(Ivana_pL.GhettoFemto_mT_rstar,"h_GhettoFemto_mT_rstar_pL");
  Ivana_pL.GhettoFemto_mT_rcore->ComputeError();
  TH2F* h_GhettoFemto_mT_rcore_pL = Convert_DlmHisto_TH2F(Ivana_pL.GhettoFemto_mT_rcore,"h_GhettoFemto_mT_rcore_pL");
  TGraphErrors pL_GhettoFemto_mT_rstar;
  pL_GhettoFemto_mT_rstar.SetName("pL_GhettoFemto_mT_rstar");
  pL_GhettoFemto_mT_rstar.SetMarkerStyle(20);
  pL_GhettoFemto_mT_rstar.SetMarkerSize(2);
  pL_GhettoFemto_mT_rstar.SetLineWidth(3);
  uPoint = 0;
  for(unsigned uBin=0; uBin<h_GhettoFemto_mT_rstar_pL->GetXaxis()->GetNbins(); uBin++){
    TH1F* hProj = (TH1F*)h_GhettoFemto_mT_rstar_pL->ProjectionY(TString::Format("hProj"),uBin+1,uBin+1);
    if(hProj->GetEntries()<200) {delete hProj; continue;}
    hProj->Scale(1./hProj->Integral(),"width");
    double mT = h_GhettoFemto_mT_rstar_pL->GetXaxis()->GetBinCenter(uBin+1);
    pL_GhettoFemto_mT_rstar.SetPoint(uPoint,mT*0.001,Get_reff(hProj));
    uPoint++;
    delete hProj;
  }
  TGraphErrors pL_GhettoFemto_mT_rcore;
  pL_GhettoFemto_mT_rcore.SetName("pL_GhettoFemto_mT_rcore");
  pL_GhettoFemto_mT_rcore.SetMarkerStyle(20);
  pL_GhettoFemto_mT_rcore.SetMarkerSize(2);
  pL_GhettoFemto_mT_rcore.SetLineWidth(3);
  uPoint = 0;
  for(unsigned uBin=0; uBin<h_GhettoFemto_mT_rcore_pL->GetXaxis()->GetNbins(); uBin++){
    TH1F* hProj = (TH1F*)h_GhettoFemto_mT_rcore_pL->ProjectionY(TString::Format("hProj"),uBin+1,uBin+1);
    if(hProj->GetEntries()<200) {delete hProj; continue;}
    hProj->Scale(1./hProj->Integral(),"width");
    double mT = h_GhettoFemto_mT_rcore_pL->GetXaxis()->GetBinCenter(uBin+1);
    pL_GhettoFemto_mT_rcore.SetPoint(uPoint,mT*0.001,Get_reff(hProj));
    uPoint++;
    delete hProj;
  }


  TFile fMtEff_pp(TString::Format("%s/SourceStudies/SourcePaper_pp/Systematics/gMt_pp.root",GetFemtoOutputFolder()),"read");
  TGraphErrors* gMtEff_pp = (TGraphErrors*)fMtEff_pp.Get("gMtStatGauss_pp");
  TFile fMtEff_pL(TString::Format("%s/SourceStudies/SourcePaper_pL/Systematics/gMt_pL.root",GetFemtoOutputFolder()),"read");
  TGraphErrors* gMtEff_pL = (TGraphErrors*)fMtEff_pL.Get("gMtStatGauss_pL_NLO");


  const double avg_mT = 0.001*0.5*(h_GhettoFemto_mT_rstar_pp->GetMean(1)+h_GhettoFemto_mT_rstar_pL->GetMean(1));
  const int mt_bin_pp = h_GhettoFemto_mT_rstar_pp->GetXaxis()->FindBin(avg_mT*1000.);
  const int mt_bin_pL = h_GhettoFemto_mT_rstar_pL->GetXaxis()->FindBin(avg_mT*1000.);
  TH1F* h_rstar_avgmt_pp = (TH1F*)h_GhettoFemto_mT_rstar_pp->ProjectionY(TString::Format("h_rstar_avgmt_pp"),mt_bin_pp,mt_bin_pp);
  TH1F* h_rstar_avgmt_pL = (TH1F*)h_GhettoFemto_mT_rstar_pL->ProjectionY(TString::Format("h_rstar_avgmt_pL"),mt_bin_pL,mt_bin_pL);
  TH1F* h_rcore_avgmt_pp = (TH1F*)h_GhettoFemto_mT_rcore_pp->ProjectionY(TString::Format("h_rcore_avgmt_pp"),mt_bin_pp,mt_bin_pp);
  TH1F* h_rcore_avgmt_pL = (TH1F*)h_GhettoFemto_mT_rcore_pL->ProjectionY(TString::Format("h_rcore_avgmt_pL"),mt_bin_pL,mt_bin_pL);

  double rcore_pp = Get_reff(h_rcore_avgmt_pp);
  double rcore_pL = Get_reff(h_rcore_avgmt_pL);
  double rstar_pp = Get_reff(h_rstar_avgmt_pp);
  double rstar_pL = Get_reff(h_rstar_avgmt_pL);

  TH1F* h_rstar_RsmCeca_pp = new TH1F("h_rstar_RsmCeca_pp","h_rstar_RsmCeca_pp",NumRadBins,rMin,rMax);
  TH1F* h_rstar_RsmUni_pp = new TH1F("h_rstar_RsmUni_pp","h_rstar_RsmUni_pp",NumRadBins,rMin,rMax);
  TH1F* h_rstar_RsmBB_pp = new TH1F("h_rstar_RsmBB_pp","h_rstar_RsmBB_pp",NumRadBins,rMin,rMax);
  TH1F* h_rstar_RsmCeca_pL = new TH1F("h_rstar_RsmCeca_pL","h_rstar_RsmCeca_pL",NumRadBins,rMin,rMax);
  TH1F* h_rstar_RsmUni_pL = new TH1F("h_rstar_RsmUni_pL","h_rstar_RsmUni_pL",NumRadBins,rMin,rMax);
  TH1F* h_rstar_RsmBB_pL = new TH1F("h_rstar_RsmBB_pL","h_rstar_RsmBB_pL",NumRadBins,rMin,rMax);
  for(unsigned uBin=0; uBin<NumRadBins; uBin++){
    double rstar = h_rstar_RsmCeca_pp->GetBinCenter(uBin+1);
    double parameters[2];
    parameters[0] = rcore_pp;
    parameters[1] = 2.0;
    double val;

    val = MagicSource_Ceca_pp.RootEval(&rstar,parameters);
    h_rstar_RsmCeca_pp->SetBinContent(uBin+1,val);
    h_rstar_RsmCeca_pp->SetBinError(uBin+1,1e-3);

    val = MagicSource_Uni_pp.RootEval(&rstar,parameters);
    h_rstar_RsmUni_pp->SetBinContent(uBin+1,val);
    h_rstar_RsmUni_pp->SetBinError(uBin+1,1e-3);

    val = MagicSource_BB_pp.RootEval(&rstar,parameters);
    h_rstar_RsmBB_pp->SetBinContent(uBin+1,val);
    h_rstar_RsmBB_pp->SetBinError(uBin+1,1e-3);

    val = MagicSource_Ceca_pL.RootEval(&rstar,parameters);
    h_rstar_RsmCeca_pL->SetBinContent(uBin+1,val);
    h_rstar_RsmCeca_pL->SetBinError(uBin+1,1e-3);

    val = MagicSource_Uni_pL.RootEval(&rstar,parameters);
    h_rstar_RsmUni_pL->SetBinContent(uBin+1,val);
    h_rstar_RsmUni_pL->SetBinError(uBin+1,1e-3);

    val = MagicSource_BB_pL.RootEval(&rstar,parameters);
    h_rstar_RsmBB_pL->SetBinContent(uBin+1,val);
    h_rstar_RsmBB_pL->SetBinError(uBin+1,1e-3);
  }

  double reff_RsmCeca_pp = Get_reff(h_rstar_RsmCeca_pp);
  double reff_RsmUni_pp = Get_reff(h_rstar_RsmUni_pp);
  double reff_RsmBB_pp = Get_reff(h_rstar_RsmBB_pp);
  double reff_RsmCeca_pL = Get_reff(h_rstar_RsmCeca_pL);
  double reff_RsmUni_pL = Get_reff(h_rstar_RsmUni_pL);
  double reff_RsmBB_pL = Get_reff(h_rstar_RsmBB_pL);

  printf("<mT> = %.2f GeV\n",avg_mT);
  printf(" CECA rcore pp  = %.2f fm\n", rcore_pp);
  printf(" CECA reff (pp) = %.2f (%.2f) fm\n", rstar_pp, pp_GhettoFemto_mT_rstar.Eval(avg_mT));
  printf("  RsmCeca       = %.2f fm\n", reff_RsmCeca_pp);
  printf("  RsmUni        = %.2f fm\n", reff_RsmUni_pp);
  printf("  RsmBB         = %.2f fm\n", reff_RsmBB_pp);
  printf(" Data reff (pp) = %.2f fm\n", gMtEff_pp->Eval(avg_mT));
  printf("-------------------\n");
  printf(" CECA rcore pL  = %.2f fm\n", rcore_pL);
  printf(" CECA reff (pL) = %.2f (%.2f) fm\n", rstar_pL, pL_GhettoFemto_mT_rstar.Eval(avg_mT));
  printf("  RsmCeca       = %.2f fm\n", reff_RsmCeca_pL);
  printf("  RsmUni        = %.2f fm\n", reff_RsmUni_pL);
  printf("  RsmBB         = %.2f fm\n", reff_RsmBB_pL);
  printf(" Data reff (pL) = %.2f fm\n", gMtEff_pL->Eval(avg_mT));
  printf("  CECA Δr = %.2f fm\n",pL_GhettoFemto_mT_rstar.Eval(avg_mT)-pp_GhettoFemto_mT_rstar.Eval(avg_mT));
  printf("  Data Δr = %.2f fm\n",gMtEff_pL->Eval(avg_mT)-gMtEff_pp->Eval(avg_mT));



  TFile fOutput(TString::Format("%s/FunWithCeca/Ceca_pp_pL_EqTau%i_Tau%.0f_Dsp%.0f-%.0f_Pnk%.0f-%.0f_Fct%.0f.root",
  GetFemtoOutputFolder(),
  EQUALIZE_TAU,Tau*100.,Ivana_pp.GetDisplacementT()*100.,Ivana_pp.GetDisplacementZ()*100.,
  Ivana_pp.GetHadronizationT()*100.,Ivana_pp.GetHadronizationZ()*100.,Ivana_pp.GetHadrFluctuation()*100.),"recreate");

  fOutput.cd();
  h_GhettoFemto_mT_rcore_pp->Write();
  h_rcore_avgmt_pp->Write();
  h_GhettoFemto_mT_rstar_pp->Write();
  h_rstar_avgmt_pp->Write();
  pp_GhettoFemto_mT_rstar.Write();
  gMtEff_pp->Write();
  h_GhettoFemto_mT_rcore_pL->Write();
  h_rcore_avgmt_pL->Write();
  h_GhettoFemto_mT_rstar_pL->Write();
  h_rstar_avgmt_pL->Write();
  pL_GhettoFemto_mT_rstar.Write();
  gMtEff_pL->Write();



  delete h_GhettoFemto_mT_rstar_pp;
  delete h_GhettoFemto_mT_rstar_pL;
  delete h_rstar_avgmt_pp;
  delete h_rstar_avgmt_pL;
  delete h_rcore_avgmt_pp;
  delete h_rcore_avgmt_pL;
  delete h_rstar_RsmCeca_pp;
  delete h_rstar_RsmUni_pp;
  delete h_rstar_RsmBB_pp;
  delete h_rstar_RsmCeca_pL;
  delete h_rstar_RsmUni_pL;
  delete h_rstar_RsmBB_pL;
  delete dlm_pT_eta_p;
  delete dlm_pT_eta_L;
}

//ramona
//https://indico.cern.ch/event/850345/contributions/3671333/attachments/1960666/3258436/pk_Femtoscopy_UpdatepPbandHMpp.pdf
//mT   reff
//0.85  1.8
//0.95  1.7
//1.05  1.62
//1.15  1.56
//for LK -> rcore = 1.11 +/- 0.04
void Ceca_pK_1(const bool Swap_p_with_Lambda=false){
  const double reff = 1.2*1.00;
  const double HadronSize = 0.0;
  const double HadronSlope = 0.0;
  const double EtaCut = 0.8;
  const double PtCut_p = 400;
  const double PtCut_K = 150;
  const double PtCutMax_p = 1800;
  const double PtCutMax_K = 1400;
  const bool EQUALIZE_TAU = true;
  const double TIMEOUT = 60;
  const double TempPP = 171;//171 or 154
  const double Asym = 6;
  const double PancakeT = 3.55;
  const double PancakeZ = 0.0;
  const double Tau = 2.3;
  const double PancakeFluct = 0;//in % !!!
  const double DisplT = 0.88/4.;
  const double DisplZ = 0.88/4.;
  unsigned THREADS = 0;
  const double MomSpread = 850;
  const bool THERMAL_KICK = false;

  const double Sensitivity = 0.005*1000;

  const double VolumeHI = 5000;
  const double TempHI = 154;
  const double ThKick = TempPP/sqrt(3.);
  const double EngyHI = 208.*5.02e6;
  const double ProtonEngy = 6.5e6;
  const double EngyPP = 2.*ProtonEngy;
  const double EngyDenst_PP_HI = 0.029412*TempPP-3.5294;//ratio
  const double VolumePP = TempPP?EngyPP/EngyHI*VolumeHI/EngyDenst_PP_HI:0;
  const double ProtonGamma = ProtonEngy/938.;
  //const double PancakeT = pow(3.*VolumePP*Asym/4./Pi,1./3.);
  //const double PancakeZ = PancakeT/Asym;

  const double q_CutOff = 200;

  printf("The pp volume is %.2f fm^3\n",VolumePP);
  //printf("Pancake  : %.2f x %.2f x %.2f\n", PancakeT,PancakeT,PancakeZ);
  printf("Displacement: %.2f %.2f %.2f fm\n",DisplT,DisplT,DisplZ);
  printf("Pancake     : %.2f %.2f %.2f fm\n",PancakeT,PancakeT,PancakeZ);

  const double Tau_Proton = Swap_p_with_Lambda?4.69:1.65;
  const double Tau_Kaon = 3.66;

  const double Mass_ProtonReso =  Swap_p_with_Lambda?1463:1362;
  const double Mass_KaonReso = 1054;

  const double FracProtonReso = Swap_p_with_Lambda?0.6438:0.6422;
  const double FracKaonReso = 0.476;

  const unsigned NumRadBins = 256;
  const double rMin = 0;
  const double rMax = 32;



  TH1F* h_reff = new TH1F("h_reff","h_reff",NumRadBins,rMin,rMax);
  //TH1F* h_rstar_Ceca = new TH1F("h_rstar_Ceca","h_rstar_Ceca",NumRadBins,rMin,rMax);
  TH1F* h_rstar_Ceca = NULL;
  TH1F* h_rcore_Ceca;// = new TH1F("h_rcore_Ceca","h_rcore_Ceca",NumRadBins,rMin,rMax);
  TH1F* h_rstar_Rsm = new TH1F("h_rstar_Rsm","h_rstar_Rsm",NumRadBins,rMin,rMax);
  TH1F* h_rcore_Rsm = new TH1F("h_rcore_Rsm","h_rcore_Rsm",NumRadBins,rMin,rMax);
  TH1F* h_rstar_RsmCeca = new TH1F("h_rstar_RsmCeca","h_rstar_RsmCeca",NumRadBins,rMin,rMax);

  TH1F* h_PR_AngleRcP2_Rsm = new TH1F("h_PR_AngleRcP2_Rsm","h_PR_AngleRcP2_Rsm",64,0,Pi);
  TH1F* h_RP_AngleRcP1_Rsm = new TH1F("h_RP_AngleRcP1_Rsm","h_RP_AngleRcP1_Rsm",64,0,Pi);
  TH1F* h_RR_AngleRcP2_Rsm = new TH1F("h_RR_AngleRcP2_Rsm","h_RR_AngleRcP2_Rsm",64,0,Pi);
  TH1F* h_RR_AngleRcP1_Rsm = new TH1F("h_RR_AngleRcP1_Rsm","h_RR_AngleRcP1_Rsm",64,0,Pi);
  TH1F* h_RR_AngleP1P2_Rsm = new TH1F("h_RR_AngleP1P2_Rsm","h_RR_AngleP1P2_Rsm",64,0,Pi);

  TF1* f_reff = new TF1("f_reff","[0]*4.*TMath::Pi()*x*x*pow(4.*TMath::Pi()*[1]*[1],-1.5)*exp(-(x*x)/(4.*[1]*[1]))+1.-[0]",rMin,rMax);
  TF1* f_rstar_Ceca = new TF1("f_rstar_Ceca","[0]*4.*TMath::Pi()*x*x*pow(4.*TMath::Pi()*[1]*[1],-1.5)*exp(-(x*x)/(4.*[1]*[1]))+1.-[0]",rMin,rMax);
  TF1* f_rstar_Rsm = new TF1("f_rstar_Rsm","[0]*4.*TMath::Pi()*x*x*pow(4.*TMath::Pi()*[1]*[1],-1.5)*exp(-(x*x)/(4.*[1]*[1]))+1.-[0]",rMin,rMax);
  TF1* f_rcore_Ceca = new TF1("f_rcore_Ceca","[0]*4.*TMath::Pi()*x*x*pow(4.*TMath::Pi()*[1]*[1],-1.5)*exp(-(x*x)/(4.*[1]*[1]))+1.-[0]",rMin,rMax);
  TF1* f_rcore_Rsm = new TF1("f_rcore_Rsm","[0]*4.*TMath::Pi()*x*x*pow(4.*TMath::Pi()*[1]*[1],-1.5)*exp(-(x*x)/(4.*[1]*[1]))+1.-[0]",rMin,rMax);


  f_reff->FixParameter(0,1);
  f_reff->FixParameter(1,reff);
  for(unsigned uBin=0; uBin<NumRadBins; uBin++){
    double rad = h_reff->GetBinCenter(uBin+1);
    double val = f_reff->Eval(rad);
    double err = 1e-3;
    h_reff->SetBinContent(uBin+1,val);
    h_reff->SetBinError(uBin+1,err);
  }

  //set up the RSM
  DLM_CleverMcLevyResoTM MagicSource;
  MagicSource.InitStability(1,2-1e-6,2+1e-6);
  MagicSource.InitScale(38,0.15,2.0);
  MagicSource.InitRad(257*2,0,64);
  MagicSource.InitType(2);
  MagicSource.SetUpReso(0,FracProtonReso);
  MagicSource.SetUpReso(1,FracKaonReso);
  MagicSource.InitNumMcIter(1000000);

  DLM_CleverMcLevyResoTM MagicSourceCeca;
  MagicSourceCeca.InitStability(1,2-1e-6,2+1e-6);
  MagicSourceCeca.InitScale(38,0.15,2.0);
  MagicSourceCeca.InitRad(257*2,0,64);
  MagicSourceCeca.InitType(2);
  MagicSourceCeca.SetUpReso(0,FracProtonReso);
  MagicSourceCeca.SetUpReso(1,FracKaonReso);
  MagicSourceCeca.InitNumMcIter(1000000);

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
  DLM_Random RanGen(11);
  double RanVal1;
  double RanVal2;
  double RanVal3;

  gROOT->cd();
  TString NKname = Swap_p_with_Lambda?"LK":"pK";
  TNtuple* ntNK = new TNtuple("nt_"+NKname,"nt_"+NKname,"Type:k_D:P1:P2:M1:M2:Tau1:Tau2:AngleRcP1:AngleRcP2:AngleP1P2");
  float nt_array[11];

  TFile* F_EposDisto_p_KaonReso = new TFile(TString::Format("%s/CatsFiles/Source/EposAngularDist/ForRamona_p_KaonReso.root",GetCernBoxDimi()));
  TNtuple* T_EposDisto_p_KaonReso = (TNtuple*)F_EposDisto_p_KaonReso->Get("InfoTuple_ClosePairs");
  unsigned N_EposDisto_p_KaonReso = T_EposDisto_p_KaonReso->GetEntries();
  T_EposDisto_p_KaonReso->SetBranchAddress("k_D",&k_D);
  T_EposDisto_p_KaonReso->SetBranchAddress("P1",&fP1);
  T_EposDisto_p_KaonReso->SetBranchAddress("P2",&fP2);
  T_EposDisto_p_KaonReso->SetBranchAddress("M1",&fM1);
  T_EposDisto_p_KaonReso->SetBranchAddress("M2",&fM2);
  T_EposDisto_p_KaonReso->SetBranchAddress("Tau1",&Tau1);
  T_EposDisto_p_KaonReso->SetBranchAddress("Tau2",&Tau2);
  T_EposDisto_p_KaonReso->SetBranchAddress("AngleRcP1",&AngleRcP1);
  T_EposDisto_p_KaonReso->SetBranchAddress("AngleRcP2",&AngleRcP2);
  T_EposDisto_p_KaonReso->SetBranchAddress("AngleP1P2",&AngleP1P2);
  for(unsigned uEntry=0; uEntry<N_EposDisto_p_KaonReso; uEntry++){
      T_EposDisto_p_KaonReso->GetEntry(uEntry);
      Tau1 = 0;
      Tau2 = Tau_Kaon;
      fM2 = Mass_KaonReso;
      if(k_D>q_CutOff) continue;
      RanVal2 = RanGen.Exponential(fM2/(fP2*Tau2));
      MagicSource.AddBGT_PR(RanVal2,cos(AngleRcP2));
      h_PR_AngleRcP2_Rsm->Fill(AngleRcP2);
      nt_array[0]=1; nt_array[1]=k_D; nt_array[2]=0; nt_array[3]=fP2;
      nt_array[4]=0; nt_array[5]=fM2; nt_array[6]=0; nt_array[7]=Tau2;
      nt_array[8]=0; nt_array[9]=AngleRcP2; nt_array[10]=0;
      gROOT->cd();
      ntNK->Fill(nt_array);
      F_EposDisto_p_KaonReso->cd();
  }
  delete F_EposDisto_p_KaonReso;

  TFile* F_EposDisto_pReso_Kaon = new TFile(TString::Format("%s/CatsFiles/Source/EposAngularDist/ForRamona_pReso_Kaon.root",GetCernBoxDimi()));
  TNtuple* T_EposDisto_pReso_Kaon = (TNtuple*)F_EposDisto_pReso_Kaon->Get("InfoTuple_ClosePairs");
  unsigned N_EposDisto_pReso_Kaon = T_EposDisto_pReso_Kaon->GetEntries();
  T_EposDisto_pReso_Kaon->SetBranchAddress("k_D",&k_D);
  T_EposDisto_pReso_Kaon->SetBranchAddress("P1",&fP1);
  T_EposDisto_pReso_Kaon->SetBranchAddress("P2",&fP2);
  T_EposDisto_pReso_Kaon->SetBranchAddress("M1",&fM1);
  T_EposDisto_pReso_Kaon->SetBranchAddress("M2",&fM2);
  T_EposDisto_pReso_Kaon->SetBranchAddress("Tau1",&Tau1);
  T_EposDisto_pReso_Kaon->SetBranchAddress("Tau2",&Tau2);
  T_EposDisto_pReso_Kaon->SetBranchAddress("AngleRcP1",&AngleRcP1);
  T_EposDisto_pReso_Kaon->SetBranchAddress("AngleRcP2",&AngleRcP2);
  T_EposDisto_pReso_Kaon->SetBranchAddress("AngleP1P2",&AngleP1P2);
  for(unsigned uEntry=0; uEntry<N_EposDisto_pReso_Kaon; uEntry++){
      T_EposDisto_pReso_Kaon->GetEntry(uEntry);
      Tau1 = Tau_Proton;
      Tau2 = 0;
      fM1 = Mass_ProtonReso;
      if(k_D>q_CutOff) continue;
      RanVal1 = RanGen.Exponential(fM1/(fP1*Tau1));
      MagicSource.AddBGT_RP(RanVal1,cos(AngleRcP1));
      h_RP_AngleRcP1_Rsm->Fill(AngleRcP1);
      nt_array[0]=10; nt_array[1]=k_D; nt_array[2]=fP1; nt_array[3]=0;
      nt_array[4]=fM1; nt_array[5]=0; nt_array[6]=Tau1; nt_array[7]=0;
      nt_array[8]=AngleRcP1; nt_array[9]=0; nt_array[10]=0;
      gROOT->cd();
      ntNK->Fill(nt_array);
      F_EposDisto_pReso_Kaon->cd();
  }
  delete F_EposDisto_pReso_Kaon;

  TFile* F_EposDisto_pReso_KaonReso = new TFile(TString::Format("%s/CatsFiles/Source/EposAngularDist/ForRamona_pReso_KaonReso.root",GetCernBoxDimi()));
  TNtuple* T_EposDisto_pReso_KaonReso = (TNtuple*)F_EposDisto_pReso_KaonReso->Get("InfoTuple_ClosePairs");
  unsigned N_EposDisto_pReso_KaonReso = T_EposDisto_pReso_KaonReso->GetEntries();
  T_EposDisto_pReso_KaonReso->SetBranchAddress("k_D",&k_D);
  T_EposDisto_pReso_KaonReso->SetBranchAddress("P1",&fP1);
  T_EposDisto_pReso_KaonReso->SetBranchAddress("P2",&fP2);
  T_EposDisto_pReso_KaonReso->SetBranchAddress("M1",&fM1);
  T_EposDisto_pReso_KaonReso->SetBranchAddress("M2",&fM2);
  T_EposDisto_pReso_KaonReso->SetBranchAddress("Tau1",&Tau1);
  T_EposDisto_pReso_KaonReso->SetBranchAddress("Tau2",&Tau2);
  T_EposDisto_pReso_KaonReso->SetBranchAddress("AngleRcP1",&AngleRcP1);
  T_EposDisto_pReso_KaonReso->SetBranchAddress("AngleRcP2",&AngleRcP2);
  T_EposDisto_pReso_KaonReso->SetBranchAddress("AngleP1P2",&AngleP1P2);
  for(unsigned uEntry=0; uEntry<N_EposDisto_pReso_KaonReso; uEntry++){
      T_EposDisto_pReso_KaonReso->GetEntry(uEntry);
      Tau1 = Tau_Proton;
      Tau2 = Tau_Kaon;
      fM1 = Mass_ProtonReso;
      fM2 = Mass_KaonReso;
      if(k_D>q_CutOff) continue;
      RanVal1 = RanGen.Exponential(fM1/(fP1*Tau1));
      RanVal2 = RanGen.Exponential(fM2/(fP2*Tau2));
      MagicSource.AddBGT_RR(RanVal1,cos(AngleRcP1),RanVal2,cos(AngleRcP2),cos(AngleP1P2));
      h_RR_AngleRcP1_Rsm->Fill(AngleRcP1);
      h_RR_AngleRcP2_Rsm->Fill(AngleRcP2);
      h_RR_AngleP1P2_Rsm->Fill(AngleP1P2);
      nt_array[0]=11; nt_array[1]=k_D; nt_array[2]=fP1; nt_array[3]=fP2;
      nt_array[4]=fM1; nt_array[5]=fM2; nt_array[6]=Tau1; nt_array[7]=Tau2;
      nt_array[8]=AngleRcP1; nt_array[9]=AngleRcP2; nt_array[10]=AngleP1P2;
      gROOT->cd();
      ntNK->Fill(nt_array);
      F_EposDisto_pReso_KaonReso->cd();
  }
  delete F_EposDisto_pReso_KaonReso;

  h_PR_AngleRcP2_Rsm->Sumw2();
  h_RP_AngleRcP1_Rsm->Sumw2();
  h_RR_AngleRcP1_Rsm->Sumw2();
  h_RR_AngleRcP2_Rsm->Sumw2();
  h_RR_AngleP1P2_Rsm->Sumw2();

  h_PR_AngleRcP2_Rsm->Scale(1./h_PR_AngleRcP2_Rsm->Integral(),"width");
  h_RP_AngleRcP1_Rsm->Scale(1./h_RP_AngleRcP1_Rsm->Integral(),"width");
  h_RR_AngleRcP1_Rsm->Scale(1./h_RR_AngleRcP1_Rsm->Integral(),"width");
  h_RR_AngleRcP2_Rsm->Scale(1./h_RR_AngleRcP2_Rsm->Integral(),"width");
  h_RR_AngleP1P2_Rsm->Scale(1./h_RR_AngleP1P2_Rsm->Integral(),"width");

  ///////////////////////////

  TREPNI Database(0);
  Database.SetSeed(11);
  std::vector<TreParticle*> ParticleList;
  ParticleList.push_back(Database.NewParticle("Proton"));
  ParticleList.push_back(Database.NewParticle("Kaon"));
  ParticleList.push_back(Database.NewParticle("Pion"));

  ParticleList.push_back(Database.NewParticle("ProtonReso"));
  ParticleList.push_back(Database.NewParticle("KaonReso"));

  double ExpPP,ExpPR,ExpRP,ExpRR;
  ExpPP = (1.-FracProtonReso)*(1.-FracKaonReso);
  ExpPR = (1.-FracProtonReso)*FracKaonReso;
  ExpRP = FracProtonReso*(1.-FracKaonReso);
  ExpRR = FracProtonReso*FracKaonReso;

  for(TreParticle* prt : ParticleList){
    if(prt->GetName()=="Proton"){
      prt->SetMass(Swap_p_with_Lambda?Mass_L:Mass_p);
      prt->SetAbundance((1.-FracProtonReso));
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);
      //prt->SetDelayTau(0.5);
      prt->SetAcceptance_pT(PtCut_p,PtCutMax_p);
    }
    else if(prt->GetName()=="Kaon"){
      prt->SetMass(Mass_Kch);
      prt->SetAbundance((1.-FracKaonReso));
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);
      //prt->SetDelayTau(1.0);
      prt->SetAcceptance_pT(PtCut_K,PtCutMax_K);
    }
    else if(prt->GetName()=="Pion"){
      prt->SetMass(Mass_pic);
      prt->SetAbundance(0);
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);
    }
    else if(prt->GetName()=="ProtonReso"){
      prt->SetMass(Mass_ProtonReso);
      prt->SetAbundance(FracProtonReso);
      prt->SetWidth(hbarc/Tau_Proton);
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);
      //prt->SetDelayTau(0.6);
      prt->SetAcceptance_pT(PtCut_p,PtCutMax_p);

      prt->NewDecay();
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Proton"));
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
      prt->GetDecay(0)->SetBranching(100);
    }
    else if(prt->GetName()=="KaonReso"){
      prt->SetMass(Mass_KaonReso);
      prt->SetAbundance(FracKaonReso);
      prt->SetWidth(hbarc/Tau_Kaon);
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);
      //prt->SetDelayTau(0.0);
      prt->SetAcceptance_pT(PtCut_K,PtCutMax_K);

      prt->NewDecay();
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Kaon"));
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
      prt->GetDecay(0)->SetBranching(100);
    }

    //prt->SetPtEtaPhi(MomSpread*);
    //prt->SetPtEtaPhi(prt->GetMass()*MomSpread*0.001);
    prt->SetPtPz(prt->GetMass()*MomSpread*0.001,prt->GetMass()*MomSpread*0.001);

    prt->SetAcceptance_Eta(-EtaCut,EtaCut);
  }
  std::vector<std::string> ListOfParticles;
  ListOfParticles.push_back("Proton");
//ListOfParticles.push_back("Proton");
  ListOfParticles.push_back("Kaon");
  CECA Ivana(Database,ListOfParticles);
  Ivana.SetDisplacementT(DisplT);
  Ivana.SetDisplacementZ(DisplZ);
  Ivana.SetHadronizationT(PancakeT);
  Ivana.SetHadronizationZ(PancakeZ);
  Ivana.SetHadrFluctuation(PancakeFluct);
  Ivana.SetTau(Tau);
  Ivana.SetThermalKick(THERMAL_KICK*ThKick);
  Ivana.SetTargetStatistics(10);
  Ivana.SetEventMult(2);
  Ivana.SetSourceDim(2);
  Ivana.SetThreadTimeout(TIMEOUT);
  Ivana.SetFemtoRegion(q_CutOff*0.5);
  Ivana.GHETTO_EVENT = true;
  printf("Goal for reff = %.4f fm; <r*> = %.4f fm\n",reff,h_reff->GetMean());
  printf("--- CECA source ---\n");
  Ivana.SetDebugMode(false);
  Ivana.GoBabyGo(THREADS);

  double r_fit;
  double adjust;

  Ivana.GhettoFemto_rstar->ComputeError();
  h_rstar_Ceca = Convert_DlmHisto_TH1F(Ivana.GhettoFemto_rstar,"h_rstar_Ceca");
  r_fit = Get_reff(h_rstar_Ceca);
  printf(" Starting parameters: Disp = %.2f %.2f fm; Hadr = %.2f %.2f fm\n",
  Ivana.GetDisplacementT(),Ivana.GetDisplacementZ(),Ivana.GetHadronizationT(),Ivana.GetHadronizationZ());
  double DEV = fabs(r_fit-reff)/fabs(reff);
  while(DEV>Sensitivity){
    adjust = reff/r_fit;
    Ivana.SetDisplacementT(Ivana.GetDisplacementT()*adjust);
    Ivana.SetDisplacementZ(Ivana.GetDisplacementZ()*adjust);
    Ivana.SetHadronizationT(Ivana.GetHadronizationT()*adjust);
    Ivana.SetHadronizationZ(Ivana.GetHadronizationZ()*adjust);
    Ivana.SetTau(Ivana.GetTau()*adjust);
    printf(" Redo with (DEV=%.2f%%): Disp = %.2f %.2f; Hadr = %.2f %.2f; Tau = %.2f\n",
    DEV*100.,Ivana.GetDisplacementT(),Ivana.GetDisplacementZ(),Ivana.GetHadronizationT(),Ivana.GetHadronizationZ(), Ivana.GetTau());
    Ivana.GoBabyGo(THREADS);
    delete h_rstar_Ceca;
    Ivana.GhettoFemto_rstar->ComputeError();
    h_rstar_Ceca = Convert_DlmHisto_TH1F(Ivana.GhettoFemto_rstar,"h_rstar_Ceca");
    r_fit = Get_reff(h_rstar_Ceca);
    DEV = fabs(r_fit-reff)/fabs(reff);
  }

  Ivana.SetUp_RSM = &MagicSourceCeca;
  std::vector<float*> buffer_ntNK;
  Ivana.Buffer_RSM = &buffer_ntNK;
  Ivana.GoBabyGo(THREADS);
  Ivana.SetUp_RSM = NULL;
  Ivana.Buffer_RSM = NULL;
  for(float* entry : buffer_ntNK){
  ntNK->Fill(entry);
  delete [] entry;
  entry = NULL;
  }

  printf(" Final parameters: Disp = %.2f %.2f fm; Hadr = %.2f %.2f fm => reff = %.4f <r*> = %.4f\n",
  Ivana.GetDisplacementT(),Ivana.GetDisplacementZ(),Ivana.GetHadronizationT(),Ivana.GetHadronizationZ(),r_fit,h_rstar_Ceca->GetMean());
  h_rstar_Ceca->Scale(1./h_rstar_Ceca->Integral(),"width");
  printf("Creating the MagicSourceCeca\n");
  Ivana.SetUp_RSM = &MagicSourceCeca;
  Ivana.GoBabyGo(THREADS);
  Ivana.SetUp_RSM = NULL;
  printf("  --> done\n");


  Ivana.Ghetto_RP_AngleRcP1->ComputeError();
  Ivana.Ghetto_PR_AngleRcP2->ComputeError();
  Ivana.Ghetto_RR_AngleRcP1->ComputeError();
  Ivana.Ghetto_RR_AngleRcP2->ComputeError();
  Ivana.Ghetto_RR_AngleP1P2->ComputeError();
  TH1F* h_RP_AngleRcP1_Ceca = Convert_DlmHisto_TH1F(Ivana.Ghetto_RP_AngleRcP1,"h_RP_AngleRcP1_Ceca");
  TH1F* h_PR_AngleRcP2_Ceca = Convert_DlmHisto_TH1F(Ivana.Ghetto_PR_AngleRcP2,"h_PR_AngleRcP2_Ceca");
  TH1F* h_RR_AngleRcP1_Ceca = Convert_DlmHisto_TH1F(Ivana.Ghetto_RR_AngleRcP1,"h_RR_AngleRcP1_Ceca");
  TH1F* h_RR_AngleRcP2_Ceca = Convert_DlmHisto_TH1F(Ivana.Ghetto_RR_AngleRcP2,"h_RR_AngleRcP2_Ceca");
  TH1F* h_RR_AngleP1P2_Ceca = Convert_DlmHisto_TH1F(Ivana.Ghetto_RR_AngleP1P2,"h_RR_AngleP1P2_Ceca");

  h_RP_AngleRcP1_Ceca->Scale(1./h_RP_AngleRcP1_Ceca->Integral(),"width");
  h_PR_AngleRcP2_Ceca->Scale(1./h_PR_AngleRcP2_Ceca->Integral(),"width");
  h_RR_AngleRcP1_Ceca->Scale(1./h_RR_AngleRcP1_Ceca->Integral(),"width");
  h_RR_AngleRcP2_Ceca->Scale(1./h_RR_AngleRcP2_Ceca->Integral(),"width");
  h_RR_AngleP1P2_Ceca->Scale(1./h_RR_AngleP1P2_Ceca->Integral(),"width");

  Ivana.GhettoFemto_mT_kstar->ComputeError();
  TH2F* h_GhettoFemto_mT_kstar = Convert_DlmHisto_TH2F(Ivana.GhettoFemto_mT_kstar,"GhettoFemto_mT_kstar");

  Ivana.GhettoFemto_mT_rstar->ComputeError();
  TH2F* h_GhettoFemto_mT_rstar = Convert_DlmHisto_TH2F(Ivana.GhettoFemto_mT_rstar,"GhettoFemto_mT_rstar");

  TGraphErrors g_GhettoFemto_mT_rstar;
  g_GhettoFemto_mT_rstar.SetName("g_GhettoFemto_mT_rstar");
  g_GhettoFemto_mT_rstar.SetMarkerStyle(20);
  g_GhettoFemto_mT_rstar.SetMarkerSize(2);
  g_GhettoFemto_mT_rstar.SetLineWidth(3);

  unsigned uPoint = 0;
  for(unsigned uBin=0; uBin<h_GhettoFemto_mT_rstar->GetXaxis()->GetNbins(); uBin++){
    TH1F* hProj = (TH1F*)h_GhettoFemto_mT_rstar->ProjectionY(TString::Format("hProj"),uBin+1,uBin+1);
    if(hProj->GetEntries()<200) {delete hProj; continue;}
    hProj->Scale(1./hProj->Integral(),"width");
    double mT = h_GhettoFemto_mT_rstar->GetXaxis()->GetBinCenter(uBin+1);
    g_GhettoFemto_mT_rstar.SetPoint(uPoint,mT*0.001,Get_reff(hProj));
    uPoint++;
    delete hProj;
  }

  //TFile fMtEff(TString::Format("%s/SourceStudies/SourcePaper_pp/Systematics/gMt_pp.root",GetFemtoOutputFolder()),"read");
  //TGraphErrors* gMtEff_pp = (TGraphErrors*)fMtEff.Get("gMtStatGauss_pp");
  //TGraphErrors* gMtCore_pp = (TGraphErrors*)fMtEff.Get("gMtStatCore_pp");
  TFile fMtEff(TString::Format("%s/SourceStudies/SourcePaper_pL/Systematics/gMt_pL.root",GetFemtoOutputFolder()),"read");
  TGraphErrors* gMtEff_pp = (TGraphErrors*)fMtEff.Get("gMtStatGauss_pL_NLO");
  TGraphErrors* gMtCore_pp = (TGraphErrors*)fMtEff.Get("gMtStatCore_pL_NLO");

  printf("--- CECA core ---\n");
  Database.GetParticle("ProtonReso")->SetAbundance(0);
  Database.GetParticle("KaonReso")->SetAbundance(0);
  Ivana.GoBabyGo(THREADS);
  Ivana.GhettoFemto_rstar->ComputeError();
  h_rcore_Ceca = Convert_DlmHisto_TH1F(Ivana.GhettoFemto_rstar,"h_rcore_Ceca");
  h_rcore_Ceca->Scale(1./h_rcore_Ceca->Integral(),"width");
  printf(" Final r_core = %.4f fm; reff = %.4f fm; <r*> = %.4f\n",Get_reff(h_rcore_Ceca),r_fit,h_rcore_Ceca->GetMean());

  Ivana.GhettoFemto_mT_rstar->ComputeError();
  TH2F* h_GhettoFemto_mT_core = Convert_DlmHisto_TH2F(Ivana.GhettoFemto_mT_rstar,"h_GhettoFemto_mT_core");

  TGraphErrors g_GhettoFemto_mT_core;
  g_GhettoFemto_mT_core.SetName("g_GhettoFemto_mT_core");
  g_GhettoFemto_mT_core.SetMarkerStyle(20);
  g_GhettoFemto_mT_core.SetMarkerSize(2);
  g_GhettoFemto_mT_core.SetLineWidth(3);

  uPoint = 0;
  for(unsigned uBin=0; uBin<h_GhettoFemto_mT_core->GetXaxis()->GetNbins(); uBin++){
    TH1F* hProj = (TH1F*)h_GhettoFemto_mT_core->ProjectionY(TString::Format("hProj"),uBin+1,uBin+1);
    if(hProj->GetEntries()<200) {delete hProj; continue;}
    hProj->Scale(1./hProj->Integral(),"width");
    double mT = h_GhettoFemto_mT_core->GetXaxis()->GetBinCenter(uBin+1);
    g_GhettoFemto_mT_core.SetPoint(uPoint,mT*0.001,Get_reff(hProj));
    uPoint++;
    delete hProj;
  }


  printf("--- RSM source ---\n");
  double rcore = reff*0.819;
  //double rcore = reff*1.46;
  printf(" Starting parameter: rcore = %.4f fm\n",rcore);
  for(unsigned uBin=0; uBin<NumRadBins; uBin++){
    double rstar = h_rstar_Rsm->GetBinCenter(uBin+1);
    double parameters[2];
    parameters[0] = rcore;
    parameters[1] = 2.0;
    double val = MagicSource.RootEval(&rstar,parameters);
    h_rstar_Rsm->SetBinContent(uBin+1,val);
    h_rstar_Rsm->SetBinError(uBin+1,1e-3);

    val = MagicSourceCeca.RootEval(&rstar,parameters);
    h_rstar_RsmCeca->SetBinContent(uBin+1,val);
    h_rstar_RsmCeca->SetBinError(uBin+1,1e-3);
  }

  r_fit = Get_reff(h_rstar_Rsm);
  DEV = fabs(r_fit-reff)/fabs(reff);
  while(DEV>Sensitivity){
    rcore *= reff/r_fit;
    printf(" Redo with (DEV=%.2f%%): rcore = %.3f fm\n",DEV*100.,rcore);
    for(unsigned uBin=0; uBin<NumRadBins; uBin++){
      double rstar = h_rstar_Rsm->GetBinCenter(uBin+1);
      double parameters[2];
      parameters[0] = rcore;
      parameters[1] = 2.0;
      double val = MagicSource.RootEval(&rstar,parameters);
      h_rstar_Rsm->SetBinContent(uBin+1,val);
      h_rstar_Rsm->SetBinError(uBin+1,1e-3);

      val = MagicSourceCeca.RootEval(&rstar,parameters);
      h_rstar_RsmCeca->SetBinContent(uBin+1,val);
      h_rstar_RsmCeca->SetBinError(uBin+1,1e-3);
    }
    r_fit = Get_reff(h_rstar_Rsm);
    DEV = fabs(r_fit-reff)/fabs(reff);
  }

  printf(" Final rcore = %.4f fm; reff = %.4f fm; <r*> = %.4f\n",rcore,r_fit,h_rstar_Rsm->GetMean());
  printf("   N.B. reff Rsm vs RsmCeca: %.3f fm %.3f fm\n",r_fit,Get_reff(h_rstar_RsmCeca));
  h_rstar_Rsm->Scale(1./h_rstar_Rsm->Integral(),"width");
  h_rstar_RsmCeca->Scale(1./h_rstar_RsmCeca->Integral(),"width");

  printf("--- RSM core ---\n");
  f_rcore_Rsm->FixParameter(0,1);
  f_rcore_Rsm->FixParameter(1,rcore);
  for(unsigned uBin=0; uBin<NumRadBins; uBin++){
    double rstar = h_rcore_Rsm->GetBinCenter(uBin+1);
    h_rcore_Rsm->SetBinContent(uBin+1,f_rcore_Rsm->Eval(rstar));
    h_rcore_Rsm->SetBinError(uBin+1,1e-3);
  }

  TFile fOutput(TString::Format("%s/FunWithCeca/CecaRSM_%s_EqTau%i_Tau%.0f_Dsp%.0f-%.0f_Pnk%.0f-%.0f_Fct%.0f_Kick%.0f_Hds%.0f-%.0f.XXX.root",
  GetFemtoOutputFolder(),Swap_p_with_Lambda?"LK":"pK",EQUALIZE_TAU,Tau*100.,Ivana.GetDisplacementT()*100.,Ivana.GetDisplacementZ()*100.,
  Ivana.GetHadronizationT()*100.,Ivana.GetHadronizationZ()*100.,Ivana.GetHadrFluctuation()*100.,
  ThKick*THERMAL_KICK,HadronSize*.100,HadronSlope*100.),"recreate");

  fOutput.cd();
  h_reff->Write();
  f_reff->Write();
  h_rstar_Ceca->Write();
  h_rstar_Rsm->Write();
  h_rstar_RsmCeca->Write();
  h_rcore_Ceca->Write();
  f_rcore_Ceca->Write();
  h_rcore_Rsm->Write();
  f_rcore_Rsm->Write();

  h_PR_AngleRcP2_Ceca->Write();
  h_PR_AngleRcP2_Rsm->Write();

  h_RP_AngleRcP1_Ceca->Write();
  h_RP_AngleRcP1_Rsm->Write();

  h_RR_AngleRcP1_Ceca->Write();
  h_RR_AngleRcP1_Rsm->Write();

  h_RR_AngleRcP2_Ceca->Write();
  h_RR_AngleRcP2_Rsm->Write();

  h_RR_AngleP1P2_Ceca->Write();
  h_RR_AngleP1P2_Rsm->Write();

  h_GhettoFemto_mT_rstar->Write();
  g_GhettoFemto_mT_rstar.Write();
  gMtEff_pp->Write();
  h_GhettoFemto_mT_core->Write();
  g_GhettoFemto_mT_core.Write();
  gMtCore_pp->Write();

  h_GhettoFemto_mT_kstar->Write();
  ntNK->Write();

  delete h_reff;
  if(h_rstar_Ceca){delete h_rstar_Ceca; h_rstar_Ceca=NULL;}
  delete h_rstar_Rsm;
  delete h_rstar_RsmCeca;
  delete h_rcore_Rsm;
  delete h_rcore_Ceca;
  delete h_PR_AngleRcP2_Rsm;
  delete h_RP_AngleRcP1_Rsm;
  delete h_RR_AngleRcP1_Rsm;
  delete h_RR_AngleRcP2_Rsm;
  delete h_RR_AngleP1P2_Rsm;
  delete h_PR_AngleRcP2_Ceca;
  delete h_RP_AngleRcP1_Ceca;
  delete h_RR_AngleRcP1_Ceca;
  delete h_RR_AngleRcP2_Ceca;
  delete h_RR_AngleP1P2_Ceca;
  delete h_GhettoFemto_mT_rstar;
  delete h_GhettoFemto_mT_kstar;
  delete f_reff;
  delete f_rstar_Ceca;
  delete f_rstar_Rsm;
  delete f_rcore_Ceca;
  delete f_rcore_Rsm;
  delete h_GhettoFemto_mT_core;
  delete ntNK;
}

void LK_Source(const double& CoreSize){
  //this one is c.a. 2% smaller than Rsm (based on pK comparison)
  TH1F* hSr_Ceca = new TH1F("hSr_Ceca","hSr_Ceca",512,0,64);
  TH1F* hSr_Epos = new TH1F("hSr_Epos","hSr_Epos",512,0,64);

  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());

  DLM_CleverMcLevyResoTM* Source_Epos = AnalysisObject.GaussCoreRsm_LK(200);
  DLM_CleverMcLevyResoTM* Source_Ceca = AnalysisObject.GaussCoreRsm_LK(201);


  for(unsigned uRad=0; uRad<hSr_Ceca->GetNbinsX(); uRad++){
    double rad = hSr_Ceca->GetBinCenter(uRad+1);
    double pars[2];
    pars[0] = CoreSize;
    pars[1] = 2.0;
    hSr_Epos->SetBinContent(uRad+1,Source_Epos->RootEval(&rad,pars));
    hSr_Ceca->SetBinContent(uRad+1,Source_Ceca->RootEval(&rad,pars));

  }

  double lowerlimit;
  double upperlimit;
  double redmean;

  GetCentralInterval(*hSr_Epos, 0.95, lowerlimit, upperlimit, true);
  TF1* f_Epos = new TF1("f_Epos",NormDoubleGaussSourceTF1,lowerlimit,upperlimit,4);
  redmean = hSr_Epos->GetMean()/2.3;
  f_Epos->SetParameter(0,redmean*0.5);
  f_Epos->SetParLimits(0,redmean*0.125,redmean*2);
  f_Epos->SetParameter(1,redmean*1.5);
  f_Epos->SetParLimits(1,redmean*0.5,redmean*8);
  f_Epos->SetParameter(2,0.5);//weight
  f_Epos->SetParLimits(2,0,1);
  f_Epos->SetParameter(3,0.8);//norm
  f_Epos->SetParLimits(3,0,1);
  hSr_Epos->Fit(f_Epos,"S, N, R, M");

  GetCentralInterval(*hSr_Ceca, 0.95, lowerlimit, upperlimit, true);
  TF1* f_Ceca = new TF1("f_Ceca",NormDoubleGaussSourceTF1,lowerlimit,upperlimit,4);
  redmean = hSr_Ceca->GetMean()/2.3;
  f_Ceca->SetParameter(0,redmean*0.5);
  f_Ceca->SetParLimits(0,redmean*0.125,redmean*2);
  f_Ceca->SetParameter(1,redmean*1.5);
  f_Ceca->SetParLimits(1,redmean*0.5,redmean*8);
  f_Ceca->SetParameter(2,0.5);//weight
  f_Ceca->SetParLimits(2,0,1);
  f_Ceca->SetParameter(3,0.8);//norm
  f_Ceca->SetParLimits(3,0,1);
  //f_Ceca->FixParameter(3,1.0);
  hSr_Ceca->Fit(f_Ceca,"S, N, R, M");

  TFile fOutput(TString::Format("%s/FunWithCeca/LK_Source_%.2ffm.root",GetFemtoOutputFolder(),CoreSize),"recreate");
  hSr_Epos->Write();
  f_Epos->Write();
  hSr_Ceca->Write();
  f_Ceca->Write();

  delete Source_Ceca;
  delete Source_Epos;
  delete hSr_Epos;
  delete f_Epos;
  delete hSr_Ceca;
  delete f_Ceca;
}

void pd_paper_sources(const bool& anchored){

  const unsigned NumMomBins = 50;
  const double kMin = 0;
  const double kMax = 250;

  const unsigned NumRadBins = 1024;
  const double rMin = 0;
  const double rMax = 64;

  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles/",GetCernBoxDimi()));

  CATS Cat_pp;
  Cat_pp.SetMomBins(NumMomBins,kMin,kMax);
  AnalysisObject.SetUpCats_pp(Cat_pp,"AV18","McGauss_ResoTM",0,202);

  CATS Cat_pL;
  Cat_pL.SetMomBins(NumMomBins,kMin,kMax);
  AnalysisObject.SetUpCats_pL(Cat_pL,"Chiral_Coupled_SPD","McGauss_ResoTM",11600,202);

  DLM_CleverMcLevyResoTM* s_pp = AnalysisObject.GetCleverMcLevyResoTM_pp();
  DLM_CleverMcLevyResoTM* s_pL = AnalysisObject.GetCleverMcLevyResoTM_pL();
  DLM_CleverMcLevyResoTM* s_pK = AnalysisObject.GaussCoreRsm_pK(200);

  TH1F* hs_pp = new TH1F("hs_pp","hs_pp",NumRadBins,rMin,rMax);
  TH1F* hs_pL = new TH1F("hs_pL","hs_pL",NumRadBins,rMin,rMax);
  TH1F* hs_pK = new TH1F("hs_pK","hs_pK",NumRadBins,rMin,rMax);

  for(unsigned uRad=0; uRad<NumRadBins; uRad++){
    double rad = hs_pp->GetBinCenter(uRad+1);
    double pars[2];
    pars[0] = anchored?1.10:1.10;
    pars[1] = 2.0;
    hs_pp->SetBinContent(uRad+1,s_pp->RootEval(&rad,pars));
    pars[0] = anchored?1.02:1.10;
    hs_pL->SetBinContent(uRad+1,s_pL->RootEval(&rad,pars));
    pars[0] = anchored?1.30:1.10;
    hs_pK->SetBinContent(uRad+1,s_pK->RootEval(&rad,pars));
  }

  TFile fOutput(TString::Format("%s/FunWithCeca/pd_paper_sources_%i.root",GetFemtoOutputFolder(),anchored),"recreate");
  hs_pp->Write();
  hs_pL->Write();
  hs_pK->Write();


  delete s_pK;
  delete hs_pp;
  delete hs_pL;
  delete hs_pK;
}


//last digit: propmother*2 + displ flag
//the dist flag 0: we fit disp. if it is 1 we fit hadr only
void Ceca_rcore_vs_reff(const double TIMEOUT, const int flag, const TString system, DLM_Histo<float>& rcore_vs_reff, DLM_Histo<float>& disp_vs_reff,
                          DLM_CleverMcLevyResoTM* MagicSourceCeca=NULL){

  const double EtaCut = 0.8;
  const double PtCut = 500;
  const bool EQUALIZE_TAU = true;
  const double FracProtonReso = 0.6422*1;
  const double FracLambdaReso = 0.6438*1;
  unsigned THREADS = 0;
  const double q_CutOff = 200;
  const bool PropMother = (flag%10)%2;
  printf("Ceca_rcore_vs_reff called:\n");
  printf(" TIMEOUT = %.0f\n",TIMEOUT);
  printf(" flag = %i\n",flag);
  printf(" system = %s\n",system.Data());
  printf(" PropMother = %i\n",PropMother);

  double disp_min = rcore_vs_reff.GetLowEdge(0)/1.5;
  double disp_max = rcore_vs_reff.GetUpEdge(0);
  //2 or 3
  if((flag%10)/2==1){
    disp_min = 3.0;
    disp_max = 10.0;
  }

  //const unsigned disp_bins = (rcore_vs_reff.GetNbins()*3)/2;
  const unsigned disp_bins = rcore_vs_reff.GetNbins();
  //const unsigned disp_bins = 8;
  const double disp_step = (disp_max-disp_min)/double(disp_bins);
  TGraph g_disp_reff;
  g_disp_reff.SetName("g_disp_reff");
  TGraph g_rcore_disp;
  g_rcore_disp.SetName("g_rcore_disp");

  int Type = 0;
  if(system=="pp") Type = 1;
  else if(system=="pL") Type = 2;
  else {printf("Cannot do it\n"); return;}

  const double Tau_Proton = 1.65;
  const double Tau_Lambda = 4.69;//4.69

  const double Mass_ProtonReso = 1362;
  const double Mass_LambdaReso = 1462;

  const unsigned NumRadBins = 256;
  const double rMin = 0;
  const double rMax = 32;

  DLM_Histo<float>* dlm_pT_eta_p = GetPtEtaDist(
    TString::Format("%s/Jaime/p_pT.root",GetCernBoxDimi()),
    TString::Format("%s/Jaime/ap_pT.root",GetCernBoxDimi()));

  DLM_Histo<float>* dlm_pT_eta_L = GetPtEtaDist(
    TString::Format("%s/Jaime/L_pT.root",GetCernBoxDimi()),
    TString::Format("%s/Jaime/aL_pT.root",GetCernBoxDimi()));


  TREPNI Database(0);
  Database.SetSeed(11);

  std::vector<TreParticle*> ParticleList;
  ParticleList.push_back(Database.NewParticle("Proton"));
  ParticleList.push_back(Database.NewParticle("Lambda"));
  ParticleList.push_back(Database.NewParticle("Pion"));
  ParticleList.push_back(Database.NewParticle("ProtonReso"));
  ParticleList.push_back(Database.NewParticle("LambdaReso"));

  for(TreParticle* prt : ParticleList){
    if(prt->GetName()=="Proton"){
      prt->SetMass(Mass_p);
      prt->SetAbundance((1.-FracProtonReso));
      prt->SetPtEtaPhi(*dlm_pT_eta_p);
    }
    else if(prt->GetName()=="Lambda"){
      prt->SetMass(Mass_L);
      if(Type==1) prt->SetAbundance(0);
      else prt->SetAbundance((1.-FracLambdaReso));
      prt->SetPtEtaPhi(*dlm_pT_eta_L);
    }
    else if(prt->GetName()=="Pion"){
      prt->SetMass(Mass_pic);
      prt->SetAbundance(0);
    }
    else if(prt->GetName()=="ProtonReso"){
      prt->SetMass(Mass_ProtonReso);
      prt->SetAbundance(FracProtonReso);
      prt->SetWidth(hbarc/Tau_Proton);

      prt->NewDecay();
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Proton"));
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
      prt->GetDecay(0)->SetBranching(100);

      prt->SetPtEtaPhi(*dlm_pT_eta_p);
    }
    else if(prt->GetName()=="LambdaReso"){
      prt->SetMass(Mass_LambdaReso);
      if(Type==1) prt->SetAbundance(0);
      else prt->SetAbundance(FracLambdaReso);
      prt->SetWidth(hbarc/Tau_Lambda);

      prt->NewDecay();
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Lambda"));
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
      prt->GetDecay(0)->SetBranching(100);

      prt->SetPtEtaPhi(*dlm_pT_eta_L);
    }

    //prt->SetPtPz(prt->GetMass()*MomSpread*0.001,prt->GetMass()*MomSpread*0.001);
    prt->SetAcceptance_Eta(-EtaCut,EtaCut);
    prt->SetAcceptance_pT(PtCut,1e6);
  }

  std::vector<std::string> ListOfParticles;
  ListOfParticles.push_back("Proton");
  ListOfParticles.push_back(system=="pp"?"Proton":"Lambda");
  CECA Ivana(Database,ListOfParticles);
  Ivana.SetTargetStatistics(10);
  Ivana.SetEventMult(2);
  Ivana.SetSourceDim(2);
  Ivana.SetThreadTimeout(TIMEOUT);
  Ivana.SetFemtoRegion(q_CutOff*0.5);
  Ivana.EqualizeFsiTime(EQUALIZE_TAU);
  Ivana.SetPropagateMother(PropMother);
  Ivana.GHETTO_EVENT = true;

  for(unsigned uDisp=0; uDisp<disp_bins; uDisp++){
    //printf("%u / %u\n",uDisp,disp_bins);
    double disp = disp_min+disp_step*double(uDisp);
    if((flag%10)/2==0){
      Ivana.SetDisplacement(disp);
      Ivana.SetHadronization(0);
    }
    else{
      Ivana.SetDisplacement(0);
      Ivana.SetHadronization(disp);
    }
    if(uDisp==disp_bins/2) Ivana.SetUp_RSM = MagicSourceCeca;
    //if(uDisp==0) Ivana.SetUp_RSM = MagicSourceCeca;
    else Ivana.SetUp_RSM = NULL;
    Ivana.GoBabyGo(THREADS);

    Ivana.GhettoFemto_rstar->ComputeError();
    TH1F* h_rstar_Ceca = Convert_DlmHisto_TH1F(Ivana.GhettoFemto_rstar,"h_rstar_Ceca");
    h_rstar_Ceca->Scale(1./h_rstar_Ceca->Integral(),"width");
    g_disp_reff.SetPoint(uDisp,Get_reff(h_rstar_Ceca),disp);

    Ivana.GhettoFemto_rcore->ComputeError();
    TH1F* h_rcore_Ceca = Convert_DlmHisto_TH1F(Ivana.GhettoFemto_rcore,"h_rcore_Ceca");
    h_rcore_Ceca->Scale(1./h_rcore_Ceca->Integral(),"width");
    g_rcore_disp.SetPoint(uDisp,disp,Get_reff(h_rcore_Ceca));

    delete h_rstar_Ceca;
    delete h_rcore_Ceca;
  }

  for(unsigned uBin=0; uBin<rcore_vs_reff.GetNbins(); uBin++){
    double reff = rcore_vs_reff.GetBinCenter(0,uBin);
    double disp = g_disp_reff.Eval(reff);
    double rcore = g_rcore_disp.Eval(disp);
    rcore_vs_reff.SetBinContent(uBin,rcore);
    disp_vs_reff.SetBinContent(uBin,disp);
  }

}




//FLAG==0 Ceca is set up normally
//FLAG==1 Ceca propagates the mothers
void HowWrongIsOurSource(const int FLAG){

  //TFile fMtEff_pp(TString::Format("%s/SourceStudies/SourcePaper_pp/Systematics/gMt_pp.root",GetFemtoOutputFolder()),"read");
  //TGraphErrors* gMtEff_pp = (TGraphErrors*)fMtEff_pp.Get("gMtStatGauss_pp");
  //TFile fMtEff_pL(TString::Format("%s/SourceStudies/SourcePaper_pL/Systematics/gMt_pL.root",GetFemtoOutputFolder()),"read");
  //TGraphErrors* gMtEff_pL = (TGraphErrors*)fMtEff_pL.Get("gMtStatGauss_pL_NLO");

  const double TIMEOUT = 16;
  const double mT_min = 1.2124;
  const double mT_max = 2.2116;
  const unsigned mT_bins = 16;
  const double mT_step = (mT_max-mT_min)/double(mT_bins);

  //const unsigned NumRadBins = 256;
  //const double rMin = 0;
  //const double rMax = 32;

  const unsigned NumMtBins_pp = 7;
  double* BinRange_pp = new double[NumMtBins_pp+1];
  double* BinCenter_pp = new double[NumMtBins_pp];
  BinCenter_pp[0] = 1.1077;
  BinCenter_pp[1] = 1.1683;
  BinCenter_pp[2] = 1.2284;
  BinCenter_pp[3] = 1.3156;
  BinCenter_pp[4] = 1.4628;
  BinCenter_pp[5] = 1.6872;
  BinCenter_pp[6] = 2.2116;
  BinRange_pp[0] = BinCenter_pp[0]-BinCenter_pp[1]+BinCenter_pp[0];
  for(int i=1; i<=6; i++) BinRange_pp[i] = (BinCenter_pp[i-1]+BinCenter_pp[i])*0.5;
  BinRange_pp[7] = BinCenter_pp[6]+BinCenter_pp[6]-BinCenter_pp[5];

  const unsigned NumMtBins_pL = 6;
  double* BinRange_pL = new double[NumMtBins_pL];
  double* BinCenter_pL = new double[NumMtBins_pL];
  BinCenter_pL[0] = 1.2124;
  BinCenter_pL[1] = 1.2896;
  BinCenter_pL[2] = 1.376;
  BinCenter_pL[3] = 1.5407;
  BinCenter_pL[4] = 1.756;
  BinCenter_pL[5] = 2.2594;
  BinRange_pL[0] = BinCenter_pL[0]-BinCenter_pL[1]+BinCenter_pL[0];
  for(int i=1; i<=5; i++) BinRange_pL[i] = (BinCenter_pL[i-1]+BinCenter_pL[i])*0.5;
  BinRange_pL[6] = BinCenter_pL[5]+BinCenter_pL[5]-BinCenter_pL[4];

  DLM_Histo<float> dlm_reff_pp;
  dlm_reff_pp.SetUp(1);
  dlm_reff_pp.SetUp(0,NumMtBins_pp,BinRange_pp,BinCenter_pp);
  dlm_reff_pp.Initialize();

  DLM_Histo<float> dlm_reff_pL;
  dlm_reff_pL.SetUp(1);
  dlm_reff_pL.SetUp(0,NumMtBins_pL,BinRange_pL,BinCenter_pL);
  dlm_reff_pL.Initialize();

  dlm_reff_pp.SetBinContent(unsigned(0),1.3717);
  dlm_reff_pp.SetBinError(unsigned(0),0.02648);
  dlm_reff_pp.SetBinContent(1,1.3005);
  dlm_reff_pp.SetBinError(1,0.02477);
  dlm_reff_pp.SetBinContent(2,1.2709);
  dlm_reff_pp.SetBinError(2,0.02654);
  dlm_reff_pp.SetBinContent(3,1.2138);
  dlm_reff_pp.SetBinError(3,0.02491);
  dlm_reff_pp.SetBinContent(4,1.1423);
  dlm_reff_pp.SetBinError(4,0.02936);
  dlm_reff_pp.SetBinContent(5,1.0482);
  dlm_reff_pp.SetBinError(5,0.02433);
  dlm_reff_pp.SetBinContent(6,0.95678);
  dlm_reff_pp.SetBinError(6,0.02014);

  dlm_reff_pL.SetBinContent(unsigned(0),1.45835);
  dlm_reff_pL.SetBinError(unsigned(0),0.12951);
  dlm_reff_pL.SetBinContent(1,1.39835);
  dlm_reff_pL.SetBinError(1,0.12552);
  dlm_reff_pL.SetBinContent(2,1.38725);
  dlm_reff_pL.SetBinError(2,0.09747);
  dlm_reff_pL.SetBinContent(3,1.3022);
  dlm_reff_pL.SetBinError(3,0.11191);
  dlm_reff_pL.SetBinContent(4,1.19395);
  dlm_reff_pL.SetBinError(4,0.10734);
  dlm_reff_pL.SetBinContent(5,1.017725);
  dlm_reff_pL.SetBinError(5,0.11363);





  DLM_Histo<float> dlm_rcore_pp;
  dlm_rcore_pp.SetUp(1);
  dlm_rcore_pp.SetUp(0,NumMtBins_pp,BinRange_pp,BinCenter_pp);
  dlm_rcore_pp.Initialize();

  DLM_Histo<float> dlm_rcore_pL;
  dlm_rcore_pL.SetUp(1);
  dlm_rcore_pL.SetUp(0,NumMtBins_pL,BinRange_pL,BinCenter_pL);
  dlm_rcore_pL.Initialize();

  dlm_rcore_pp.SetBinContent(unsigned(0),1.3064);
  dlm_rcore_pp.SetBinError(unsigned(0),sqrt(pow(0.027805,2.)+pow(0.0085539,2.)));
  dlm_rcore_pp.SetBinContent(1,1.2316);
  dlm_rcore_pp.SetBinError(1,sqrt(pow(0.022773,2.)+pow(0.010501,2.)));
  dlm_rcore_pp.SetBinContent(2,1.2006);
  dlm_rcore_pp.SetBinError(2,sqrt(pow(0.022552,2.)+pow(0.014732,2.)));
  dlm_rcore_pp.SetBinContent(3,1.1402);
  dlm_rcore_pp.SetBinError(3,sqrt(pow(0.025014,2.)+pow(0.011086,2.)));
  dlm_rcore_pp.SetBinContent(4,1.0628);
  dlm_rcore_pp.SetBinError(4,sqrt(pow(0.025221,2.)+pow(0.010027,2.)));
  dlm_rcore_pp.SetBinContent(5,0.96238);
  dlm_rcore_pp.SetBinError(5,sqrt(pow(0.025234,2.)+pow(0.0098228,2.)));
  dlm_rcore_pp.SetBinContent(6,0.86503);
  dlm_rcore_pp.SetBinError(6,sqrt(pow(0.020369,2.)+pow(0.010841,2.)));


  for(unsigned uBin=0; uBin<6; uBin++){
    double nlo_val,nlo_stat,nlo_syst;
    double lo_val,lo_stat,lo_syst;
    switch (uBin) {
      case 0:   nlo_val = 1.2634;
                nlo_syst = 0.0876;
                nlo_stat = 0.041731;
                lo_val = 1.1816;
                lo_syst = 0.060069;
                lo_stat = 0.042102;
                break;
      case 1:   nlo_val = 1.1926;
                nlo_syst = 0.10035;
                nlo_stat = 0.050217;
                lo_val = 1.1381;
                lo_syst = 0.057701;
                lo_stat = 0.048758;
                break;
      case 2:   nlo_val = 1.1804;
                nlo_syst = 0.065382;
                nlo_stat = 0.042162;
                lo_val = 1.1279;
                lo_syst = 0.070412;
                lo_stat = 0.041525;
                break;
      case 3:   nlo_val = 1.0902;
                nlo_syst = 0.092647;
                nlo_stat = 0.035456;
                lo_val = 1.0516;
                lo_syst = 0.055087;
                lo_stat = 0.033123;
                break;
      case 4:   nlo_val = 0.9717;
                nlo_syst = 0.0825;
                nlo_stat = 0.041914;
                lo_val = 0.94642;
                lo_syst = 0.055;
                lo_stat = 0.038272;
                break;
      case 5:   nlo_val = 0.73832;
                nlo_syst = 0.18077;
                nlo_stat = 0.03668;
                lo_val = 0.82714;
                lo_syst = 0.050213;
                lo_stat = 0.033311;
                break;
      default:  break;
    }
    double avg_val = 0.5*(nlo_val+lo_val);
    double avg_stat = 0.5*(nlo_stat+lo_stat);
    double avg_syst = 0.5*(nlo_syst+lo_syst);
    double diff = nlo_val-lo_val;
    double tot_err = sqrt(avg_stat*avg_stat+avg_syst*avg_syst+diff*diff);
    dlm_rcore_pL.SetBinContent(uBin,avg_val);
    dlm_rcore_pL.SetBinError(uBin,tot_err);
  }


  DLM_Histo<float> dlm_reff_diff(dlm_reff_pL);
  dlm_reff_diff = dlm_reff_pL;
  dlm_reff_diff -= dlm_reff_pp;

  TGraphErrors g_reff_pp;
  g_reff_pp.SetName("g_reff_pp");
  g_reff_pp.SetFillColorAlpha(kBlue+1,0.3);
  g_reff_pp.SetLineColor(kBlue+1);
  g_reff_pp.SetLineWidth(3);
  for(unsigned uBin=0; uBin<dlm_reff_pp.GetNbins(); uBin++){
    g_reff_pp.SetPoint(uBin,dlm_reff_pp.GetBinCenter(0,uBin),dlm_reff_pp.GetBinContent(uBin));
    g_reff_pp.SetPointError(uBin,0,dlm_reff_pp.GetBinError(uBin));
  }

  TGraphErrors g_reff_pL;
  g_reff_pL.SetName("g_reff_pL");
  g_reff_pL.SetFillColorAlpha(kRed+1,0.3);
  g_reff_pL.SetLineColor(kRed+1);
  g_reff_pL.SetLineWidth(3);
  for(unsigned uBin=0; uBin<dlm_reff_pL.GetNbins(); uBin++){
    g_reff_pL.SetPoint(uBin,dlm_reff_pL.GetBinCenter(0,uBin),dlm_reff_pL.GetBinContent(uBin));
    g_reff_pL.SetPointError(uBin,0,dlm_reff_pL.GetBinError(uBin));
  }


  TGraphErrors g_rcore_pp;
  g_rcore_pp.SetName("g_core_pp");
  g_rcore_pp.SetFillColorAlpha(kBlue+1,0.3);
  g_rcore_pp.SetLineColor(kBlue+1);
  g_rcore_pp.SetLineWidth(3);
  for(unsigned uBin=0; uBin<dlm_rcore_pp.GetNbins(); uBin++){
    g_rcore_pp.SetPoint(uBin,dlm_rcore_pp.GetBinCenter(0,uBin),dlm_rcore_pp.GetBinContent(uBin));
    g_rcore_pp.SetPointError(uBin,0,dlm_rcore_pp.GetBinError(uBin));
  }

  TGraphErrors g_rcore_pL;
  g_rcore_pL.SetName("g_rcore_pL");
  g_rcore_pL.SetFillColorAlpha(kRed+1,0.3);
  g_rcore_pL.SetLineColor(kRed+1);
  g_rcore_pL.SetLineWidth(3);
  for(unsigned uBin=0; uBin<dlm_rcore_pL.GetNbins(); uBin++){
    g_rcore_pL.SetPoint(uBin,dlm_rcore_pL.GetBinCenter(0,uBin),dlm_rcore_pL.GetBinContent(uBin));
    g_rcore_pL.SetPointError(uBin,0,dlm_rcore_pL.GetBinError(uBin));
  }


  TGraphErrors gDiffHalf;
  gDiffHalf.SetName("gDiffHalf");
  gDiffHalf.SetFillColorAlpha(kGreen+1,0.3);
  gDiffHalf.SetLineColor(kGreen+1);
  gDiffHalf.SetLineWidth(3);
  TGraphErrors gDiff;
  gDiff.SetName("gDiff");
  gDiff.SetFillColorAlpha(kOrange+1,0.3);
  gDiff.SetLineColor(kOrange+1);
  gDiff.SetLineWidth(3);
  for(unsigned uBin=0; uBin<mT_bins; uBin++){
    double mT = double(uBin)*mT_step+mT_min;
    gDiffHalf.SetPoint(uBin,mT,dlm_reff_diff.Eval(&mT));
    gDiffHalf.SetPointError(uBin,0,dlm_reff_diff.Eval(&mT,true)*0.5);
    gDiff.SetPoint(uBin,mT,dlm_reff_diff.Eval(&mT));
    gDiff.SetPointError(uBin,0,dlm_reff_diff.Eval(&mT,true));
  }

  //hos
  DLM_Histo<float> dlm_Ceca_rcore_vs_reff_pp;
  dlm_Ceca_rcore_vs_reff_pp.SetUp(1);
  dlm_Ceca_rcore_vs_reff_pp.SetUp(0,16,0.92,1.42);
  dlm_Ceca_rcore_vs_reff_pp.Initialize();

  DLM_Histo<float> dlm_Ceca_disp_vs_reff_pp;
  dlm_Ceca_disp_vs_reff_pp.SetUp(1);
  dlm_Ceca_disp_vs_reff_pp.SetUp(0,16,0.92,1.42);
  dlm_Ceca_disp_vs_reff_pp.Initialize();

  DLM_Histo<float> dlm_Ceca_rcore_vs_reff_pL;
  dlm_Ceca_rcore_vs_reff_pL.SetUp(1);
  dlm_Ceca_rcore_vs_reff_pL.SetUp(0,24,0.88,1.64);
  dlm_Ceca_rcore_vs_reff_pL.Initialize();

  DLM_Histo<float> dlm_Ceca_disp_vs_reff_pL;
  dlm_Ceca_disp_vs_reff_pL.SetUp(1);
  dlm_Ceca_disp_vs_reff_pL.SetUp(0,24,0.88,1.64);
  dlm_Ceca_disp_vs_reff_pL.Initialize();

  DLM_CleverMcLevyResoTM MagicSource_Epos_pp;
  DLM_CleverMcLevyResoTM MagicSource_EposWs_pp;
  DLM_CleverMcLevyResoTM MagicSource_Flat_pp;
  DLM_CleverMcLevyResoTM MagicSource_Epos_pL;
  DLM_CleverMcLevyResoTM MagicSource_Flat_pL;
  SetUp_RSM_pp(MagicSource_Epos_pp,GetCernBoxDimi(),false);
  SetUp_RSM_pp(MagicSource_Flat_pp,GetCernBoxDimi(),10);
  SetUp_RSM_pp(MagicSource_EposWs_pp,GetCernBoxDimi(),true);
  SetUp_RSM_pL(MagicSource_Epos_pL,GetCernBoxDimi());
  SetUp_RSM_pL(MagicSource_Flat_pL,GetCernBoxDimi(),10);

  DLM_CleverMcLevyResoTM MagicSource_RsmCeca_pp;
  DLM_CleverMcLevyResoTM MagicSource_RsmCeca_pL;
  BasicSetUp_MS(MagicSource_RsmCeca_pp,0.6422,0.6422);
  BasicSetUp_MS(MagicSource_RsmCeca_pL,0.6422,0.6438);

  //settings flag, pp or pL, histo rcore, disp,
  Ceca_rcore_vs_reff(TIMEOUT,FLAG,"pp",dlm_Ceca_rcore_vs_reff_pp,dlm_Ceca_disp_vs_reff_pp,&MagicSource_RsmCeca_pp);
  Ceca_rcore_vs_reff(TIMEOUT,FLAG,"pL",dlm_Ceca_rcore_vs_reff_pL,dlm_Ceca_disp_vs_reff_pL,&MagicSource_RsmCeca_pL);

  TH1F* h_Ceca_rcore_vs_reff_pp = Convert_DlmHisto_TH1F(&dlm_Ceca_rcore_vs_reff_pp,"h_Ceca_rcore_vs_reff_pp");
  TH1F* h_Ceca_rcore_vs_reff_pL = Convert_DlmHisto_TH1F(&dlm_Ceca_rcore_vs_reff_pL,"h_Ceca_rcore_vs_reff_pL");

  TH1F* h_Ceca_disp_vs_reff_pp = Convert_DlmHisto_TH1F(&dlm_Ceca_disp_vs_reff_pp,"h_Ceca_disp_vs_reff_pp");
  TH1F* h_Ceca_disp_vs_reff_pL = Convert_DlmHisto_TH1F(&dlm_Ceca_disp_vs_reff_pL,"h_Ceca_disp_vs_reff_pL");

  TGraphErrors mT_rcore_pp_Ceca;
  mT_rcore_pp_Ceca.SetName("mT_rcore_pp_Ceca");
  TGraphErrors mT_rcore_pL_Ceca;
  mT_rcore_pL_Ceca.SetName("mT_rcore_pL_Ceca");

  TGraphErrors mT_rcore_pp_EposWs;
  mT_rcore_pp_EposWs.SetName("mT_rcore_pp_EposWs");
  TGraphErrors mT_rcore_pL_EposWs;
  mT_rcore_pL_EposWs.SetName("mT_rcore_pL_EposWs");

  TGraphErrors mT_rcore_pp_Epos;
  mT_rcore_pp_Epos.SetName("mT_rcore_pp_Epos");
  TGraphErrors mT_rcore_pL_Epos;
  mT_rcore_pL_Epos.SetName("mT_rcore_pL_Epos");

  TGraphErrors mT_reff_pL_Epos;
  mT_reff_pL_Epos.SetName("mT_reff_pL_Epos");

  TGraphErrors mT_rcore_pp_Flat;
  mT_rcore_pp_Flat.SetName("mT_rcore_pp_Flat");
  TGraphErrors mT_rcore_pL_Flat;
  mT_rcore_pL_Flat.SetName("mT_rcore_pL_Flat");

  TGraphErrors mT_rcore_pp_RsmCeca;
  mT_rcore_pp_RsmCeca.SetName("mT_rcore_pp_RsmCeca");
  TGraphErrors mT_rcore_pL_RsmCeca;
  mT_rcore_pL_RsmCeca.SetName("mT_rcore_pL_RsmCeca");

  double Chi2_rcore_Epos = 0;
  double Chi2_rcore_Flat = 0;
  double Chi2_rcore_EposWs = 0;
  double Chi2_rcore_Ceca = 0;
  double Chi2_rcore_RsmCeca = 0;
  //do this for all but the last pL bin (we skip the last, too large error and no pp data)
  //for(unsigned umT=0; umT<NumMtBins_pL-1; umT++){
  for(unsigned umT=0; umT<NumMtBins_pL; umT++){
    //double mT = double(umT)*mT_step+mT_min;
    double mT = BinCenter_pL[umT];

    double reff_pp = dlm_reff_pp.Eval(&mT);
    double refferr_pp = dlm_reff_pp.Eval(&mT,true);
    double reff_pL = dlm_reff_pL.Eval(&mT);
    double refferr_pL = dlm_reff_pL.Eval(&mT,true);

    double min_reff_pp = reff_pp-refferr_pp;
    double max_reff_pp = reff_pp+refferr_pp;
    double min_reff_pL = reff_pL-refferr_pL;
    double max_reff_pL = reff_pL+refferr_pL;

    //goal: for each mT bin, based on the reff in pp and pL,
    //evaluate the corresponding rcore (pp and pL) and the correponding errors
    //for the errors: simply evaluate the +/-1sigma limits and assume the rcore
    //is the average of those two
    //than: evaluate the descrepancy for each mT bin for the hypothesis of a common core
    //combine the descrepancy for all bins to get a total significance. Hope to stay within 1 sigma, 2 max

    double min_rcore_pp, max_rcore_pp, min_rcore_pL, max_rcore_pL, rcore_pp, rcore_pL, rcore_diff, rcore_err_pp, rcore_err_pL, rcore_err_diff;
    double rcore_nsig;

    ///////////////

    min_rcore_pp = dlm_Ceca_rcore_vs_reff_pp.Eval(&min_reff_pp);
    max_rcore_pp = dlm_Ceca_rcore_vs_reff_pp.Eval(&max_reff_pp);
    min_rcore_pL = dlm_Ceca_rcore_vs_reff_pL.Eval(&min_reff_pL);
    max_rcore_pL = dlm_Ceca_rcore_vs_reff_pL.Eval(&max_reff_pL);

    rcore_pp = (max_rcore_pp+min_rcore_pp)*0.5;
    rcore_err_pp = (max_rcore_pp-min_rcore_pp)*0.5;
    rcore_pL = (max_rcore_pL+min_rcore_pL)*0.5;
    rcore_err_pL = (max_rcore_pL-min_rcore_pL)*0.5;
    rcore_diff = rcore_pL-rcore_pp;
    rcore_err_diff = sqrt(rcore_err_pp*rcore_err_pp+rcore_err_pL*rcore_err_pL);
    rcore_nsig = fabs(rcore_diff/rcore_err_diff);
    Chi2_rcore_Ceca += rcore_nsig*rcore_nsig;

    mT_rcore_pp_Ceca.SetPoint(umT,mT,rcore_pp);
    mT_rcore_pp_Ceca.SetPointError(umT,0,rcore_err_pp);

    mT_rcore_pL_Ceca.SetPoint(umT,mT,rcore_pL);
    mT_rcore_pL_Ceca.SetPointError(umT,0,rcore_err_pL);

    ///////////////

    max_rcore_pp = GetRcore(MagicSource_EposWs_pp,max_reff_pp);
    min_rcore_pp = GetRcore(MagicSource_EposWs_pp,min_reff_pp);
    max_rcore_pL = GetRcore(MagicSource_Epos_pL,max_reff_pL);
    min_rcore_pL = GetRcore(MagicSource_Epos_pL,min_reff_pL);

    rcore_pp = (max_rcore_pp+min_rcore_pp)*0.5;
    rcore_err_pp = (max_rcore_pp-min_rcore_pp)*0.5;
    rcore_pL = (max_rcore_pL+min_rcore_pL)*0.5;
    rcore_err_pL = (max_rcore_pL-min_rcore_pL)*0.5;
    rcore_diff = rcore_pL-rcore_pp;
    rcore_err_diff = sqrt(rcore_err_pp*rcore_err_pp+rcore_err_pL*rcore_err_pL);
    rcore_nsig = fabs(rcore_diff/rcore_err_diff);
    Chi2_rcore_EposWs += rcore_nsig*rcore_nsig;

    mT_rcore_pp_EposWs.SetPoint(umT,mT,rcore_pp);
    mT_rcore_pp_EposWs.SetPointError(umT,0,rcore_err_pp);

    mT_rcore_pL_EposWs.SetPoint(umT,mT,rcore_pL);
    mT_rcore_pL_EposWs.SetPointError(umT,0,rcore_err_pL);
    ///////////////

    max_rcore_pp = GetRcore(MagicSource_Epos_pp,max_reff_pp);
    min_rcore_pp = GetRcore(MagicSource_Epos_pp,min_reff_pp);
    max_rcore_pL = GetRcore(MagicSource_Epos_pL,max_reff_pL);
    min_rcore_pL = GetRcore(MagicSource_Epos_pL,min_reff_pL);

    rcore_pp = (max_rcore_pp+min_rcore_pp)*0.5;
    rcore_err_pp = (max_rcore_pp-min_rcore_pp)*0.5;
    rcore_pL = (max_rcore_pL+min_rcore_pL)*0.5;
    rcore_err_pL = (max_rcore_pL-min_rcore_pL)*0.5;
    rcore_diff = rcore_pL-rcore_pp;
    rcore_err_diff = sqrt(rcore_err_pp*rcore_err_pp+rcore_err_pL*rcore_err_pL);
    rcore_nsig = fabs(rcore_diff/rcore_err_diff);
    Chi2_rcore_Epos += rcore_nsig*rcore_nsig;

    mT_rcore_pp_Epos.SetPoint(umT,mT,rcore_pp);
    mT_rcore_pp_Epos.SetPointError(umT,0,rcore_err_pp);

    mT_rcore_pL_Epos.SetPoint(umT,mT,rcore_pL);
    mT_rcore_pL_Epos.SetPointError(umT,0,rcore_err_pL);

//printf("%f %f\n",max_rcore_pp,min_rcore_pp);
//printf(" %f %f\n",GetReff(MagicSource_Epos_pp,max_rcore_pp),GetReff(MagicSource_Epos_pp,min_rcore_pp));
//printf(" %f %f\n",GetReff(MagicSource_Epos_pL,max_rcore_pp),GetReff(MagicSource_Epos_pL,min_rcore_pp));
    mT_reff_pL_Epos.SetPoint(umT,mT,0.5*(GetReff(MagicSource_Epos_pL,max_rcore_pp)+GetReff(MagicSource_Epos_pL,min_rcore_pp)));
    mT_reff_pL_Epos.SetPointError(umT,0,0.5*(GetReff(MagicSource_Epos_pL,max_rcore_pp)-GetReff(MagicSource_Epos_pL,min_rcore_pp)));

    ///////////////

    max_rcore_pp = GetRcore(MagicSource_Flat_pp,max_reff_pp);
    min_rcore_pp = GetRcore(MagicSource_Flat_pp,min_reff_pp);
    max_rcore_pL = GetRcore(MagicSource_Flat_pL,max_reff_pL);
    min_rcore_pL = GetRcore(MagicSource_Flat_pL,min_reff_pL);

    rcore_pp = (max_rcore_pp+min_rcore_pp)*0.5;
    rcore_err_pp = (max_rcore_pp-min_rcore_pp)*0.5;
    rcore_pL = (max_rcore_pL+min_rcore_pL)*0.5;
    rcore_err_pL = (max_rcore_pL-min_rcore_pL)*0.5;
    rcore_diff = rcore_pL-rcore_pp;
    rcore_err_diff = sqrt(rcore_err_pp*rcore_err_pp+rcore_err_pL*rcore_err_pL);
    rcore_nsig = fabs(rcore_diff/rcore_err_diff);
    Chi2_rcore_Flat += rcore_nsig*rcore_nsig;

    mT_rcore_pp_Flat.SetPoint(umT,mT,rcore_pp);
    mT_rcore_pp_Flat.SetPointError(umT,0,rcore_err_pp);

    mT_rcore_pL_Flat.SetPoint(umT,mT,rcore_pL);
    mT_rcore_pL_Flat.SetPointError(umT,0,rcore_err_pL);

    ///////////////

    max_rcore_pp = GetRcore(MagicSource_RsmCeca_pp,max_reff_pp);
    min_rcore_pp = GetRcore(MagicSource_RsmCeca_pp,min_reff_pp);
    max_rcore_pL = GetRcore(MagicSource_RsmCeca_pL,max_reff_pL);
    min_rcore_pL = GetRcore(MagicSource_RsmCeca_pL,min_reff_pL);

    rcore_pp = (max_rcore_pp+min_rcore_pp)*0.5;
    rcore_err_pp = (max_rcore_pp-min_rcore_pp)*0.5;
    rcore_pL = (max_rcore_pL+min_rcore_pL)*0.5;
    rcore_err_pL = (max_rcore_pL-min_rcore_pL)*0.5;
    rcore_diff = rcore_pL-rcore_pp;
    rcore_err_diff = sqrt(rcore_err_pp*rcore_err_pp+rcore_err_pL*rcore_err_pL);
    rcore_nsig = fabs(rcore_diff/rcore_err_diff);
    Chi2_rcore_RsmCeca += rcore_nsig*rcore_nsig;

    mT_rcore_pp_RsmCeca.SetPoint(umT,mT,rcore_pp);
    mT_rcore_pp_RsmCeca.SetPointError(umT,0,rcore_err_pp);

    mT_rcore_pL_RsmCeca.SetPoint(umT,mT,rcore_pL);
    mT_rcore_pL_RsmCeca.SetPointError(umT,0,rcore_err_pL);
  }
  //double Nsigma_rcore = GetNsigma(Chi2_rcore,NumMtBins_pL-1);


  TFile fOutput(TString::Format("%s/FunWithCeca/HowWrongIsOurSource/HowWrongIsOurSource_FLAG%i.root",GetFemtoOutputFolder(),FLAG),"recreate");
  //gDiffHalf.Write();
  //gDiff.Write();
  g_reff_pp.Write();
  g_reff_pL.Write();
  g_rcore_pp.Write();
  g_rcore_pL.Write();
  h_Ceca_rcore_vs_reff_pp->Write();
  h_Ceca_rcore_vs_reff_pL->Write();
  h_Ceca_disp_vs_reff_pp->Write();
  h_Ceca_disp_vs_reff_pL->Write();
  mT_rcore_pp_Ceca.Write();
  mT_rcore_pL_Ceca.Write();
  mT_rcore_pp_EposWs.Write();
  mT_rcore_pL_EposWs.Write();
  mT_rcore_pp_Epos.Write();
  mT_rcore_pL_Epos.Write();
  mT_rcore_pp_Flat.Write();
  mT_rcore_pL_Epos.Write();
  mT_rcore_pp_RsmCeca.Write();
  mT_rcore_pL_RsmCeca.Write();
  mT_reff_pL_Epos.Write();

  delete [] BinRange_pp;
  delete [] BinCenter_pp;
  delete [] BinRange_pL;
  delete [] BinCenter_pL;

  delete h_Ceca_rcore_vs_reff_pp;
  delete h_Ceca_rcore_vs_reff_pL;
  delete h_Ceca_disp_vs_reff_pp;
  delete h_Ceca_disp_vs_reff_pL;
}


//a parallel cross check of Jaime's analysis
//i.e. the extraction of the mT scaling for pp/pL, where we do save the
//S(r,mT) and the corresponding fit. We also compare the shape of the
//theoretical C(k,mT) based on the "true" source and the Gaussian approximation.
//regarding the true source, there we only take pairs below 100 MeV, at some point we should
//look at the k* dependence.
void Sources_In_SourcePaper(const TString WhichSystem, const int mode){

  const float TIMEOUT = 60*16;
  const float k_CutOff = 108;
  const unsigned NumThreads = 0;
  const float EtaCut = 0.8;
  const float PrimProt = 0.3578;
  const float PrimLamb = 0.3562;
  const double TauProt = 1.65;
  const double TauLamb = 4.69;//4.69
  const double Mass_ProtR = 1362;
  const double Mass_LambR = 1462;
  //multiply by that factor (for CATS with RSM)
  const double BiasCore = 0.96;

  //const double MeanMt = 1.55;
  //const double MeanMtRange = 0.075;


  const unsigned NumMomBins = 27*1;
  const double kMin = 0;
  const double kMax = 324;

  std::vector<float> mT_val;
  std::vector<float> mT_reff;
  TString Name_mT_reff;
  TString Name_mT_rcore;
  if(WhichSystem=="pp"){
    Name_mT_reff = "g_reff_pp";
    Name_mT_rcore = "g_rcore_pp";
  }
  else if(WhichSystem=="pL"){
    Name_mT_reff = "g_reff_pL";
    Name_mT_rcore = "g_rcore_pL";
  }
  else{
    printf("WhichSystem == WTF ???\n");
    return;
  }
  TFile fMt(TString::Format("%s/Jaime/SourcePaper_Published.root",GetCernBoxDimi()),"read");
  TGraph* g_mT_reff = (TGraph*)fMt.Get(Name_mT_reff);
  TGraph* g_mT_rcore = (TGraph*)fMt.Get(Name_mT_rcore);

  //const unsigned NumMt = g_mT_reff->GetN();

  float dispT,dispZ;
  float hadrT,hadrZ,hadr_flc;
  float tau,tau_flc,tau_prp;
  float hdr_size,hdr_slope;
  float thm_kick;
  bool fixed_hadron;
  float frag_beta;

  switch (mode) {
    //jaime parameters
    case 1:
      dispT = 0.25*1.03;
      dispZ = 0.25*1.03;
      hadrT = 2.3*1.03;
      hadrZ = 0;
      hadr_flc = 0;
      tau = 3.55*1.03;
      tau_flc = 0;
      tau_prp = true;
      hdr_size = 0.0;
      hdr_slope = 0.0;
      thm_kick = 0;
      fixed_hadron = true;
      frag_beta = 0;
      break;
    default: printf("Unknown mode\n"); return;
  }


  TString FilePathMom = TString::Format("%s/Jaime/",GetCernBoxDimi());
  TString FileNameMomP1,FileNameMomAP1,FileNameMomP2,FileNameMomAP2;
  FileNameMomP1 = "p_pT.root";
  FileNameMomAP1 = "ap_pT.root";
  FileNameMomP2 = "L_pT.root";
  FileNameMomAP2 = "aL_pT.root";
/*
  TGraph gCkRelError;
  gCkRelError.SetName("gCkRelError");
  //this is Ck in the 0th mT bin (def). We use it to get reallisitc errors on C(k)
  //note that in all mT bins we have approx the same yield, hence similar uncertainties,
  //i.e. we asssume that those errors are constant with mT
  TFile fData_def(FilePathMom+"pL_Data_Def_TotUnc.root","read");
  TH1F* hData_def = (TH1F*)fData_def.Get("hData_0");
  for(unsigned uMomBin=0; uMomBin<hData_def->GetNbinsX(); uMomBin++){
    gCkRelError.SetPoint(uMomBin,hData_def->GetBinCenter(uMomBin+1),
                      hData_def->GetBinError(uMomBin+1)/hData_def->GetBinContent(uMomBin+1));
  }
*/
  TH1F* h_pT_p;
  TH1F* h_pT_ap;
  TH1F* h_pT_p_all;
  TH1F* h_pT_L;
  TH1F* h_pT_aL;
  TH1F* h_pT_L_all;

  TFile file_p(FilePathMom+FileNameMomP1,"read");
  h_pT_p = (TH1F*)file_p.Get("pTDist_after");
  if(!h_pT_p) printf("ISSUE with h_pT_p\n");
  //fOutput.cd();
  gROOT->cd();
  h_pT_p_all = (TH1F*)h_pT_p->Clone("h_pT_p_all");

  TFile file_ap(FilePathMom+FileNameMomAP1,"read");
  h_pT_ap = (TH1F*)file_ap.Get("pTDist_after");
  if(!h_pT_ap) printf("ISSUE with h_pT_ap\n");
  h_pT_p_all->Add(h_pT_ap);

  TFile file_L(FilePathMom+FileNameMomP2,"read");
  h_pT_L = (TH1F*)file_L.Get("pTDist_after");
  if(!h_pT_L) printf("ISSUE with h_pT_L\n");
  //fOutput.cd();
  gROOT->cd();
  h_pT_L_all = (TH1F*)h_pT_L->Clone("h_pT_L_all");

  TFile file_aL(FilePathMom+FileNameMomAP2,"read");
  h_pT_aL = (TH1F*)file_aL.Get("pTDist_after");
  if(!h_pT_aL) printf("ISSUE with h_pT_aL\n");
  h_pT_L_all->Add(h_pT_aL);

  DLM_Histo<float>* dlm_pT_p = Convert_TH1F_DlmHisto(h_pT_p_all);
  DLM_Histo<float>* dlm_pT_L = Convert_TH1F_DlmHisto(h_pT_L_all);
  dlm_pT_p->RescaleAxis(0,1000,false);
  dlm_pT_L->RescaleAxis(0,1000,false);

  double* BinRange = NULL;
  double axis[2];

  DLM_Histo<float> dlm_pT_eta_p;
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

  DLM_Histo<float> dlm_pT_eta_L;
  dlm_pT_eta_L.SetUp(2);
  BinRange = dlm_pT_L->GetBinRange(0);
  dlm_pT_eta_L.SetUp(0,dlm_pT_L->GetNbins(),BinRange);
  delete [] BinRange;
  dlm_pT_eta_L.SetUp(1,1,-EtaCut,EtaCut);
  dlm_pT_eta_L.Initialize();
  for(unsigned uBin=0; uBin<dlm_pT_L->GetNbins(); uBin++){
    dlm_pT_eta_L.SetBinContent(uBin,0,dlm_pT_L->GetBinContent(uBin));
  }

  gROOT->cd();

  TFile fOutput(TString::Format("%s/FunWithCeca/Sources_In_SourcePaper/SSP_191222_%s_%i.root",
  GetFemtoOutputFolder(),WhichSystem.Data(),mode),"recreate");
  h_pT_p_all->Write();
  h_pT_L_all->Write();
  g_mT_reff->Write();
  g_mT_rcore->Write();

  gROOT->cd();

  double ExpPP,ExpPR,ExpRP,ExpRR;
  if(WhichSystem=="pp"){
    ExpRR = (1.-PrimProt)*(1.-PrimLamb);
    ExpRP = (1.-PrimProt)*PrimLamb;
    ExpPR = PrimProt*(1.-PrimLamb);
    ExpPP = PrimProt*PrimLamb;
  }
  else{
    ExpRR = (1.-PrimProt)*(1.-PrimProt);
    ExpRP = (1.-PrimProt)*PrimProt;
    ExpPR = PrimProt*(1.-PrimProt);
    ExpPP = PrimProt*PrimProt;
  }




  TREPNI Database(0);
  Database.SetSeed(11);
  std::vector<TreParticle*> ParticleList;
  ParticleList.push_back(Database.NewParticle("Proton"));
  ParticleList.push_back(Database.NewParticle("Lambda"));
  ParticleList.push_back(Database.NewParticle("Pion"));

  ParticleList.push_back(Database.NewParticle("ProtonReso"));
  ParticleList.push_back(Database.NewParticle("LambdaReso"));

  for(TreParticle* prt : ParticleList){
    if(prt->GetName()=="Proton"){
      prt->SetMass(Mass_p);
      prt->SetAbundance(PrimProt);
      prt->SetRadius(hdr_size);
      prt->SetRadiusSlope(hdr_slope);
      prt->SetPtEtaPhi(dlm_pT_eta_p);
      }
      else if(prt->GetName()=="Lambda"){
        prt->SetMass(Mass_L);
        if(WhichSystem=="pL") prt->SetAbundance(PrimLamb);
        else prt->SetAbundance(0);
        prt->SetRadius(hdr_size);
        prt->SetRadiusSlope(hdr_slope);
        //prt->SetPtPz(prt->GetMass()*MomSpread*0.001,prt->GetMass()*MomSpread*0.001);
        prt->SetPtEtaPhi(dlm_pT_eta_L);
      }
      else if(prt->GetName()=="Pion"){
        prt->SetMass(Mass_pic);
        prt->SetAbundance(0);
        prt->SetRadius(hdr_size);
        prt->SetRadiusSlope(hdr_slope);
        //prt->SetPtPz(prt->GetMass()*MomSpread*0.001,prt->GetMass()*MomSpread*0.001);
      }
      else if(prt->GetName()=="ProtonReso"){
        prt->SetMass(Mass_ProtR);
        prt->SetAbundance(1.-PrimProt);
        prt->SetWidth(hbarc/TauProt);
        prt->SetRadius(hdr_size);
        prt->SetRadiusSlope(hdr_slope);

        prt->NewDecay();
        prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Proton"));
        prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
        prt->GetDecay(0)->SetBranching(100);

        //prt->SetPtPz(prt->GetMass()*MomSpread*0.001,prt->GetMass()*MomSpread*0.001);
        prt->SetPtEtaPhi(dlm_pT_eta_p);
      }
      else if(prt->GetName()=="LambdaReso"){
        prt->SetMass(Mass_LambR);
        if(WhichSystem=="pL") prt->SetAbundance(1.-PrimLamb);
        else prt->SetAbundance(0);
        prt->SetWidth(hbarc/TauLamb);
        prt->SetRadius(hdr_size);
        prt->SetRadiusSlope(hdr_slope);

        prt->NewDecay();
        prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Lambda"));
        prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
        prt->GetDecay(0)->SetBranching(100);

        //prt->SetPtPz(prt->GetMass()*MomSpread*0.001,prt->GetMass()*MomSpread*0.001);
        prt->SetPtEtaPhi(dlm_pT_eta_L);
      }
      prt->SetAcceptance_Eta(-EtaCut,EtaCut);
    }

    std::string part1;
    std::string part2;

    if(WhichSystem=="pp"){
      part1 = "Proton";
      part2 = "Proton";
    }
    else{
      part1 = "Proton";
      part2 = "Lambda";
    }

    std::vector<std::string> ListOfParticles;
    ListOfParticles.push_back(part1);
    ListOfParticles.push_back(part2);
    CECA Ivana(Database,ListOfParticles);
    Ivana.SetDisplacementT(dispT);
    Ivana.SetDisplacementZ(dispZ);
    Ivana.SetHadronizationT(hadrT);
    Ivana.SetHadronizationZ(hadrZ);
    Ivana.SetHadrFluctuation(hadr_flc);
    Ivana.SetTau(tau,tau_prp);
    Ivana.SetTauFluct(tau_flc);
    Ivana.SetThermalKick(thm_kick);
    Ivana.SetTargetStatistics(10);
    Ivana.SetEventMult(2);
    Ivana.SetSourceDim(2);
    Ivana.SetThreadTimeout(TIMEOUT);
    Ivana.SetFemtoRegion(k_CutOff);
    Ivana.EqualizeFsiTime(true);
    Ivana.SetFixedHadr(fixed_hadron);
    Ivana.SetFragmentBeta(frag_beta);
    Ivana.GHETTO_EVENT = true;
    //Ivana.Ghetto_NumMtBins = 9;
    //Ivana.Ghetto_MtMin = 1000;
    //Ivana.Ghetto_MtMax = 2350;
    Ivana.Ghetto_NumMtBins = 9;
    Ivana.Ghetto_MtMin = 1025;
    Ivana.Ghetto_MtMax = 2375;
    Ivana.Ghetto_NumRadBins = 128*2;
    Ivana.Ghetto_RadMin = 0;
    Ivana.Ghetto_RadMax = 32;
    Ivana.Ghetto_NumMomBins = NumMomBins*8;
    Ivana.Ghetto_MomMin = kMin;
    Ivana.Ghetto_MomMax = kMin+(kMax-kMin)*8;

    //Ivana.Ghetto_MinMt = 1000.*(MeanMt-MeanMtRange);
    //Ivana.Ghetto_MaxMt = 1000.*(MeanMt+MeanMtRange);
    //Ivana.SetDebugMode(true);
    //printf("Goal for reff = %.4f fm; <r*> = %.4f fm\n",reff,h_reff->GetMean());
    //printf("--- CECA source ---\n");
    Ivana.GoBabyGo(NumThreads);

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
    printf("   Expected  Total  Femto\n");
    printf("PP %6.2f%% %6.2f%% %6.2f\n",ExpPP*100.,TotPP*100.,FemtoPP*100.);
    printf("PR %6.2f%% %6.2f%% %6.2f\n",ExpPR*100.,TotPR*100.,FemtoPR*100.);
    printf("RP %6.2f%% %6.2f%% %6.2f\n",ExpRP*100.,TotRP*100.,FemtoRP*100.);
    printf("RR %6.2f%% %6.2f%% %6.2f\n",ExpRR*100.,TotRR*100.,FemtoRR*100.);


    Ivana.GhettoFemto_mT_rstar->ComputeError();
    //DEL?
    TH2F* h_GhettoFemto_mT_rstar = Convert_DlmHisto_TH2F(Ivana.GhettoFemto_mT_rstar,"GhettoFemto_mT_rstar");
    //TH2F* h_GhettoFemto_mT_rcore = Convert_DlmHisto_TH2F(Ivana.GhettoFemto_mT_rcore,"GhettoFemto_mT_rcore");
    TGraph gfm_mT;
    Ivana.Ghetto_kstar_rstar_mT->ComputeError();
    Ivana.Ghetto_kstar_rcore_mT->ComputeError();

    //DEL?
    TGraph* gfm_kstar = new TGraph [Ivana.GhettoFemto_mT_rstar->GetNbins(0)];

    TH2F** h_Ghetto_kstar_rstar_mT = new TH2F* [Ivana.GhettoFemto_mT_rstar->GetNbins(0)];
    DLM_Histo<float>** dlm_Ghetto_kstar_rstar_mT = new DLM_Histo<float>* [Ivana.GhettoFemto_mT_rstar->GetNbins(0)];

    TH2F** h_Ghetto_kstar_rcore_mT = new TH2F* [Ivana.GhettoFemto_mT_rcore->GetNbins(0)];
    DLM_Histo<float>** dlm_Ghetto_kstar_rcore_mT = new DLM_Histo<float>* [Ivana.GhettoFemto_mT_rcore->GetNbins(0)];

    for(unsigned umT=0; umT<Ivana.GhettoFemto_mT_rstar->GetNbins(0); umT++){
      TH1D* hProj = h_GhettoFemto_mT_rstar->ProjectionY(TString::Format("hProj"),umT+1,umT+1);
      ////this spits out tons of TF1 errors
      double Mean = hProj->GetMean();
      delete hProj;
      if(Mean==0) continue;
      double GFM = GaussFromMean(Mean);
      gfm_mT.SetPoint(umT,Ivana.GhettoFemto_mT_rstar->GetBinCenter(0,umT),GFM);

      //  //how the radius changes as a function of kstar
      gfm_kstar[umT].SetName(TString::Format("gfm_kstar_mT%.0f",Ivana.GhettoFemto_mT_rstar->GetBinCenter(0,umT)));
      gfm_kstar[umT].SetLineColor(kRed+umT);
      gfm_kstar[umT].SetLineWidth(3);
      gfm_kstar[umT].SetMarkerStyle(2);
      gfm_kstar[umT].SetMarkerColor(kBlack);
      gfm_kstar[umT].SetMarkerSize(1);
      //h_Ghetto_kstar_rstar_mT[umT] = new TH2F(TString::Format("h_Ghetto_kstar_rstar_mT%.0f",Ivana.GhettoFemto_mT_rstar->GetBinCenter(0,umT)),
      //        TString::Format("h_Ghetto_kstar_rstar_mT%.0f",Ivana.GhettoFemto_mT_rstar->GetBinCenter(0,umT)),
      //        Ivana.Ghetto_kstar_rstar_mT->GetNbins(0),Ivana.Ghetto_kstar_rstar_mT->GetLowEdge(0),Ivana.Ghetto_kstar_rstar_mT->GetUpEdge(0),
      //        Ivana.Ghetto_kstar_rstar_mT->GetNbins(1),Ivana.Ghetto_kstar_rstar_mT->GetLowEdge(1),Ivana.Ghetto_kstar_rstar_mT->GetUpEdge(1));
      dlm_Ghetto_kstar_rstar_mT[umT] = new DLM_Histo<float>();
      dlm_Ghetto_kstar_rstar_mT[umT]->SetUp(2);
      dlm_Ghetto_kstar_rstar_mT[umT]->SetUp(0,Ivana.Ghetto_kstar_rstar_mT->GetNbins(0),Ivana.Ghetto_kstar_rstar_mT->GetLowEdge(0),Ivana.Ghetto_kstar_rstar_mT->GetUpEdge(0));
      dlm_Ghetto_kstar_rstar_mT[umT]->SetUp(1,Ivana.Ghetto_kstar_rstar_mT->GetNbins(1),Ivana.Ghetto_kstar_rstar_mT->GetLowEdge(1),Ivana.Ghetto_kstar_rstar_mT->GetUpEdge(1));
      dlm_Ghetto_kstar_rstar_mT[umT]->Initialize();
      for(unsigned ukstar=0; ukstar<Ivana.Ghetto_kstar_rstar_mT->GetNbins(0); ukstar++){
        for(unsigned urstar=0; urstar<Ivana.Ghetto_kstar_rstar_mT->GetNbins(1); urstar++){
          //h_Ghetto_kstar_rstar_mT[umT]->SetBinContent(ukstar+1,urstar+1,Ivana.Ghetto_kstar_rstar_mT->GetBinContent(ukstar,urstar,umT));
          //printf("mt%u ks%u rs%u --> %.2f\n",umT,ukstar,urstar,Ivana.Ghetto_kstar_rstar_mT->GetBinContent(ukstar,urstar,umT));
          //usleep(1e3);
          dlm_Ghetto_kstar_rstar_mT[umT]->SetBinContent(ukstar,urstar,Ivana.Ghetto_kstar_rstar_mT->GetBinContent(ukstar,urstar,umT));
        }
      }
      dlm_Ghetto_kstar_rstar_mT[umT]->ComputeError();
      NormalizeSource_rk(dlm_Ghetto_kstar_rstar_mT[umT]);
      h_Ghetto_kstar_rstar_mT[umT] = Convert_DlmHisto_TH2F(dlm_Ghetto_kstar_rstar_mT[umT],TString::Format("h_Ghetto_kstar_rstar_mT%.0f",Ivana.GhettoFemto_mT_rstar->GetBinCenter(0,umT)));

      for(unsigned ukstar=0; ukstar<Ivana.Ghetto_kstar_rstar_mT->GetNbins(0); ukstar++){
        hProj = h_Ghetto_kstar_rstar_mT[umT]->ProjectionY(TString::Format("hProj"),ukstar+1,ukstar+1);
        Mean = hProj->GetMean();
        delete hProj;
        if(Mean==0) continue;
        ////this spits out tons of TF1 errors
        GFM = GaussFromMean(Mean);
        gfm_kstar[umT].SetPoint(ukstar,Ivana.Ghetto_kstar_rstar_mT->GetBinCenter(0,ukstar),GFM);
      }

      fOutput.cd();
      h_Ghetto_kstar_rstar_mT[umT]->Write();
      gfm_kstar[umT].Write();
      //delete h_Ghetto_kstar_rstar_mT;
    }




    for(unsigned umT=0; umT<Ivana.GhettoFemto_mT_rcore->GetNbins(0); umT++){
      double GFM,Mean;
      TH1D* hProj;
      dlm_Ghetto_kstar_rcore_mT[umT] = new DLM_Histo<float>();
      dlm_Ghetto_kstar_rcore_mT[umT]->SetUp(2);
      dlm_Ghetto_kstar_rcore_mT[umT]->SetUp(0,Ivana.Ghetto_kstar_rcore_mT->GetNbins(0),Ivana.Ghetto_kstar_rcore_mT->GetLowEdge(0),Ivana.Ghetto_kstar_rcore_mT->GetUpEdge(0));
      dlm_Ghetto_kstar_rcore_mT[umT]->SetUp(1,Ivana.Ghetto_kstar_rcore_mT->GetNbins(1),Ivana.Ghetto_kstar_rcore_mT->GetLowEdge(1),Ivana.Ghetto_kstar_rcore_mT->GetUpEdge(1));
      dlm_Ghetto_kstar_rcore_mT[umT]->Initialize();
      for(unsigned ukstar=0; ukstar<Ivana.Ghetto_kstar_rcore_mT->GetNbins(0); ukstar++){
        for(unsigned urcore=0; urcore<Ivana.Ghetto_kstar_rcore_mT->GetNbins(1); urcore++){
          //h_Ghetto_kstar_rcore_mT[umT]->SetBinContent(ukstar+1,urcore+1,Ivana.Ghetto_kstar_rcore_mT->GetBinContent(ukstar,urcore,umT));
          //printf("mt%u ks%u rs%u --> %.2f\n",umT,ukstar,urcore,Ivana.Ghetto_kstar_rcore_mT->GetBinContent(ukstar,urcore,umT));
          //usleep(1e3);
          dlm_Ghetto_kstar_rcore_mT[umT]->SetBinContent(ukstar,urcore,Ivana.Ghetto_kstar_rcore_mT->GetBinContent(ukstar,urcore,umT));
        }
      }
      dlm_Ghetto_kstar_rcore_mT[umT]->ComputeError();
      NormalizeSource_rk(dlm_Ghetto_kstar_rcore_mT[umT]);
      h_Ghetto_kstar_rcore_mT[umT] = Convert_DlmHisto_TH2F(dlm_Ghetto_kstar_rcore_mT[umT],TString::Format("h_Ghetto_kstar_rcore_mT%.0f",Ivana.GhettoFemto_mT_rcore->GetBinCenter(0,umT)));

      fOutput.cd();
      h_Ghetto_kstar_rcore_mT[umT]->Write();
      //delete h_Ghetto_kstar_rcore_mT;
    }





    Ivana.GhettoFemto_rstar->ComputeError();
    Ivana.GhettoFemto_rstar->ScaleToIntegral();
    Ivana.GhettoFemto_rstar->ScaleToBinSize();
    TH1F* h_GhettoFemto_rstar = Convert_DlmHisto_TH1F(Ivana.GhettoFemto_rstar,"GhettoFemto_rstar");
    fOutput.cd();
    h_GhettoFemto_rstar->SetLineWidth(3);
    h_GhettoFemto_rstar->SetLineColor(kAzure);
    h_GhettoFemto_rstar->Write();

    //Ivana.GhettoFemto_mT_rstar->ComputeError();
    //TH2F* h_GhettoFemto_mT_rstar = Convert_DlmHisto_TH2F(Ivana.GhettoFemto_mT_rstar,"GhettoFemto_mT_rstar");
    TGraphErrors g_GhettoFemto_mT_rstar_G;
    g_GhettoFemto_mT_rstar_G.SetName("g_GhettoFemto_mT_rstar_G");
    g_GhettoFemto_mT_rstar_G.SetMarkerStyle(20);
    g_GhettoFemto_mT_rstar_G.SetMarkerSize(1);
    g_GhettoFemto_mT_rstar_G.SetLineWidth(6);
    g_GhettoFemto_mT_rstar_G.SetLineColor(kAzure);

    Ivana.GhettoFemto_mT_rcore->ComputeError();
    TH2F* h_GhettoFemto_mT_rcore = Convert_DlmHisto_TH2F(Ivana.GhettoFemto_mT_rcore,"GhettoFemto_mT_rcore");
    TGraphErrors g_GhettoFemto_mT_rcore_G;
    g_GhettoFemto_mT_rcore_G.SetName("g_GhettoFemto_mT_rcore_G");
    g_GhettoFemto_mT_rcore_G.SetMarkerStyle(20);
    g_GhettoFemto_mT_rcore_G.SetMarkerSize(1);
    g_GhettoFemto_mT_rcore_G.SetLineWidth(6);
    g_GhettoFemto_mT_rcore_G.SetLineColor(kBlack);

    unsigned uPointRS = 0;
    unsigned uPointRC = 0;
    double reff_Ceca, rcore_Ceca;

    CATS Kitty_G;
    Kitty_G.SetMomBins(NumMomBins,kMin,kMax);
    DLM_CommonAnaFunctions AnalysisObject;
    AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());
    if(WhichSystem=="pp"){
      AnalysisObject.SetUpCats_pp(Kitty_G,"AV18","Gauss",0,0);
    }
    else if(WhichSystem=="pL"){
      AnalysisObject.SetUpCats_pL(Kitty_G,"Chiral_Coupled_SPD","Gauss",11600,0);
    }
    Kitty_G.SetNotifications(CATS::nWarning);

    CATS Kitty_G2;
    Kitty_G2.SetMomBins(NumMomBins,kMin,kMax);
    if(WhichSystem=="pp"){
      AnalysisObject.SetUpCats_pp(Kitty_G2,"AV18","Gauss",0,0);
    }
    else if(WhichSystem=="pL"){
      AnalysisObject.SetUpCats_pL(Kitty_G2,"Chiral_Coupled_SPD","Gauss",11600,0);
    }
    Kitty_G2.SetNotifications(CATS::nWarning);

    CATS Kitty_RSM;
    Kitty_RSM.SetMomBins(NumMomBins,kMin,kMax);
    if(WhichSystem=="pp"){
      AnalysisObject.SetUpCats_pp(Kitty_RSM,"AV18","McGauss_ResoTM",0,202);
    }
    else if(WhichSystem=="pL"){
      AnalysisObject.SetUpCats_pL(Kitty_RSM,"Chiral_Coupled_SPD","McGauss_ResoTM",11600,202);
    }

    Kitty_RSM.SetNotifications(CATS::nWarning);
    //DLM_CleverMcLevyResoTM MagicSource;
    //not set up for pp!!!
    //SetUp_RSM_pL(MagicSource,GetCernBoxDimi());

    CATS Kitty_C;
    Kitty_C.SetMomBins(NumMomBins,kMin,kMax);
    if(WhichSystem=="pp"){
      AnalysisObject.SetUpCats_pp(Kitty_C,"AV18","Gauss",0,0);
    }
    else if(WhichSystem=="pL"){
      AnalysisObject.SetUpCats_pL(Kitty_C,"Chiral_Coupled_SPD","Gauss",11600,0);
    }
    Kitty_C.SetNotifications(CATS::nWarning);

    CATS Kitty_Ckstar;
    Kitty_Ckstar.SetMomBins(NumMomBins,kMin,kMax);
    if(WhichSystem=="pp"){
      AnalysisObject.SetUpCats_pp(Kitty_Ckstar,"AV18","Gauss",0,0);
    }
    else if(WhichSystem=="pL"){
      AnalysisObject.SetUpCats_pL(Kitty_Ckstar,"Chiral_Coupled_SPD","Gauss",11600,0);
    }
    Kitty_Ckstar.SetNotifications(CATS::nWarning);

    for(unsigned uBin=0; uBin<h_GhettoFemto_mT_rstar->GetXaxis()->GetNbins(); uBin++){
      //printf("uBin=%u\n",uBin);
      TH1F* hProj = (TH1F*)h_GhettoFemto_mT_rstar->ProjectionY(TString::Format("hProj"),uBin+1,uBin+1);
      //printf("--hProj\n");
      TF1* fmT;
      TF1* fmT_full = NULL;

      TGraph* gCk_G;
      TGraph* gCk_G2;
      TGraph* gCk_C;
      TGraph* gCk_Ckstar;
      TGraph* gCk_RSM;
      TGraph* gCk_RatGC;
      TGraph* gCk_RatG2C;
      TGraph* gCk_RatRC;
      double Mean = hProj->GetMean();
      double Err = hProj->GetStdDev();
      double mT = h_GhettoFemto_mT_rstar->GetXaxis()->GetBinCenter(uBin+1);
      //printf("--BEF IF\n");
      if(Mean&&Err&&hProj->GetEntries()>320){
        //printf("--AFT IF\n");
        hProj->Scale(1./hProj->Integral(),"width");
        //printf("0\n");
        reff_Ceca = Get_reff_TF1(hProj,fmT,1,0.95);

        fmT_full = new TF1("fmT_full","[0]*4.*TMath::Pi()*x*x*pow(4.*TMath::Pi()*[1]*[1],-1.5)*exp(-(x*x)/(4.*[1]*[1]))+1.-[0]",0,16);
        fmT_full->SetNpx(1024);
        fmT_full->SetParameter(0,fmT->GetParameter(0));
        fmT_full->SetParameter(1,fmT->GetParameter(1));
      //reff_Ceca = g_mT_reff->Eval(mT*0.001);
        printf("reff_Ceca @ %.2f = %f\n",mT*0.001,reff_Ceca);
        g_GhettoFemto_mT_rstar_G.SetPoint(uPointRS,mT*0.001,reff_Ceca);
        g_GhettoFemto_mT_rstar_G.SetPointError(uPointRS,0,0);
        //printf("2\n");
        gCk_C = new TGraph();
        gCk_C->SetName(TString::Format("gCk_C_mT%u",uBin));
        gCk_C->SetLineWidth(6);
        gCk_C->SetLineColor(kRed+1);

        gCk_Ckstar = new TGraph();
        gCk_Ckstar->SetName(TString::Format("gCk_Ckstar_mT%u",uBin));
        gCk_Ckstar->SetLineWidth(6);
        gCk_Ckstar->SetLineColor(kOrange+1);

        gCk_RSM = new TGraph();
        gCk_RSM->SetName(TString::Format("gCk_RSM_mT%u",uBin));
        gCk_RSM->SetLineWidth(4);
        gCk_RSM->SetLineColor(kGreen+1);
        //printf("3\n");
        gCk_G = new TGraph();
        gCk_G->SetName(TString::Format("gCk_G_mT%u",uBin));
        gCk_G->SetLineWidth(4);
        gCk_G->SetLineColor(kBlue+1);

        gCk_G2 = new TGraph();
        gCk_G2->SetName(TString::Format("gCk_G2_mT%u",uBin));
        gCk_G2->SetLineWidth(4);
        gCk_G2->SetLineColor(kBlue-3);

        //printf("4\n");
        gCk_RatGC = new TGraph();
        gCk_RatGC->SetName(TString::Format("gCk_RatGC_mT%u",uBin));
        gCk_RatGC->SetLineWidth(6);
        gCk_RatGC->SetLineColor(kBlack);

        gCk_RatG2C = new TGraph();
        gCk_RatG2C->SetName(TString::Format("gCk_RatG2C_mT%u",uBin));
        gCk_RatG2C->SetLineWidth(6);
        gCk_RatG2C->SetLineColor(kBlack);

        gCk_RatRC = new TGraph();
        gCk_RatRC->SetName(TString::Format("gCk_RatRC_mT%u",uBin));
        gCk_RatRC->SetLineWidth(6);
        gCk_RatRC->SetLineColor(kBlack);

        Kitty_G.SetAnaSource(0,reff_Ceca);
        Kitty_G2.SetAnaSource(0,reff_Ceca*0.98);

        double DeltaREFF = g_mT_reff->Eval(mT*0.001)-reff_Ceca;
        double rcore = g_mT_rcore->Eval(mT*0.001)-DeltaREFF;
        rcore *= BiasCore;
        //here we want the core, given the effective CECA source
        //since the rcore is evaluated based on the measured reff, which is larger,
        //we need to subtract the difference to get to the CECA core
      //double rcore = GetRcore(MagicSource,reff_Ceca);
      //double rcore = g_mT_rcore->Eval(mT*0.001);
      printf("mT=%.2f; reff vs rcore = %.2f %.2f\n",mT,reff_Ceca,rcore);
        //printf("%f vs %f vs %f\n",rcore,g_mT_rcore->Eval(mT*0.001)-DeltaREFF,g_mT_rcore->Eval(mT*0.001));
        //Kitty_RSM.SetAnaSource(0,g_mT_rcore->Eval(mT*0.001)-DeltaREFF);
        Kitty_RSM.SetAnaSource(0,rcore*1.0);
        Kitty_RSM.SetAnaSource(1,2.);

        DLM_Histo<float>* dlmProj;
        dlmProj = Convert_TH1F_DlmHisto(hProj);
        //printf("Converted\n");
        dlmProj->RebinAll(2);
        //printf("Rebinned\n");
        *dlmProj/=2.;
        //printf("Rescaled\n");
        DLM_HistoSource HistoSource(*dlmProj);
        //printf("Created DLM_Histo\n");
        Kitty_C.SetAnaSource(CatsSourceForwarder, &HistoSource, 0);
        //printf("SetAnaSource done\n");

        //DLM_Histo<float>* dlmSkstar = Convert_TH2F_DlmHisto(h_Ghetto_kstar_rstar_mT[uBin]);
        DLM_HistoSource HistoSourceKstar(*dlm_Ghetto_kstar_rstar_mT[uBin]);

        //NormalizeSource_rk
  Kitty_Ckstar.SetAnaSource(CatsSourceForwarder, &HistoSourceKstar, 0);
  //Kitty_Ckstar.SetAnaSource(CatsSourceForwarder, &HistoSource, 0);
        //printf("dlmSkstar done\n");

        #pragma omp parallel sections
        {
          #pragma omp section
          {
          Kitty_G.KillTheCat();
          }
          #pragma omp section
          {
          Kitty_G2.KillTheCat();
          }
          #pragma omp section
          {
          Kitty_RSM.KillTheCat();
          }
          #pragma omp section
          {
          Kitty_C.KillTheCat(uPointRS?CATS::kSourceChanged:CATS::kAllChanged);
          }
          #pragma omp section
          {
          Kitty_Ckstar.KillTheCat(uPointRS?CATS::kSourceChanged:CATS::kAllChanged);
          }
        }

        //printf("KillTheCat done\n");
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
          gCk_G->SetPoint(uMomBin,Kitty_G.GetMomentum(uMomBin),Kitty_G.GetCorrFun(uMomBin));
          gCk_G2->SetPoint(uMomBin,Kitty_G2.GetMomentum(uMomBin),Kitty_G2.GetCorrFun(uMomBin));
          gCk_C->SetPoint(uMomBin,Kitty_C.GetMomentum(uMomBin),Kitty_C.GetCorrFun(uMomBin));
          gCk_Ckstar->SetPoint(uMomBin,Kitty_Ckstar.GetMomentum(uMomBin),Kitty_Ckstar.GetCorrFun(uMomBin));
          gCk_RSM->SetPoint(uMomBin,Kitty_RSM.GetMomentum(uMomBin),Kitty_RSM.GetCorrFun(uMomBin));
          gCk_RatGC->SetPoint(uMomBin,Kitty_C.GetMomentum(uMomBin),Kitty_G.GetCorrFun(uMomBin)/Kitty_C.GetCorrFun(uMomBin));
          gCk_RatG2C->SetPoint(uMomBin,Kitty_C.GetMomentum(uMomBin),Kitty_G2.GetCorrFun(uMomBin)/Kitty_C.GetCorrFun(uMomBin));
          gCk_RatRC->SetPoint(uMomBin,Kitty_C.GetMomentum(uMomBin),Kitty_RSM.GetCorrFun(uMomBin)/Kitty_C.GetCorrFun(uMomBin));
        }

        hProj->SetName(TString::Format("h_reff_mT%u",uBin));
        fmT->SetName(TString::Format("f_reff_mT%u",uBin));
        if(fmT_full) fmT_full->SetName(TString::Format("full_reff_mT%u",uBin));
        fOutput.cd();
        hProj->Write();
        fmT->Write();
        if(fmT_full) fmT_full->Write();
        gCk_G->Write();
        gCk_G2->Write();
        gCk_RSM->Write();
        gCk_C->Write();
        gCk_Ckstar->Write();
        gCk_RatGC->Write();
        gCk_RatG2C->Write();
        gCk_RatRC->Write();

        uPointRS++;
        //printf("Attempting to delete\n");
        delete dlmProj;
        //delete dlmSkstar;
        delete fmT;
        if(fmT_full) {delete fmT_full;fmT_full=NULL;}
        delete gCk_G;
        delete gCk_G2;
        delete gCk_RSM;
        delete gCk_C;
        delete gCk_Ckstar;
        delete gCk_RatGC;
        delete gCk_RatG2C;
        delete gCk_RatRC;
        //printf("Deleted\n");
      }

      gROOT->cd();
      delete hProj;

      //do the same for
      //GhettoFemto_mT_rcore

      hProj = (TH1F*)h_GhettoFemto_mT_rcore->ProjectionY(TString::Format("hProj"),uBin+1,uBin+1);
      Mean = hProj->GetMean();
      Err = hProj->GetStdDev();
      mT = h_GhettoFemto_mT_rcore->GetXaxis()->GetBinCenter(uBin+1);
      fmT_full = NULL;
      //printf("A\n");
      if(Mean&&Err&&hProj->GetEntries()>160){
        hProj->Scale(1./hProj->Integral(),"width");
        //printf("AA\n");
        rcore_Ceca = Get_reff_TF1(hProj,fmT,1,0.95);
        fmT_full = new TF1("fmT_full","[0]*4.*TMath::Pi()*x*x*pow(4.*TMath::Pi()*[1]*[1],-1.5)*exp(-(x*x)/(4.*[1]*[1]))+1.-[0]",0,16);
        fmT_full->SetNpx(1024);
        fmT_full->SetParameter(0,fmT->GetParameter(0));
        fmT_full->SetParameter(1,fmT->GetParameter(1));
        printf("rcore_Ceca @ %.2f = %f\n",mT*0.001,rcore_Ceca);
        g_GhettoFemto_mT_rcore_G.SetPoint(uPointRC,mT*0.001,rcore_Ceca);
        g_GhettoFemto_mT_rcore_G.SetPointError(uPointRC,0,0);

        hProj->SetName(TString::Format("h_rcore_mT%u",uBin));
        fmT->SetName(TString::Format("f_rcore_mT%u",uBin));
        if(fmT_full) fmT_full->SetName(TString::Format("full_rcore_mT%u",uBin));
        fOutput.cd();
        hProj->Write();
        fmT->Write();
        fmT_full->Write();

        uPointRC++;
        //printf("DEL\n");
        delete fmT;
        //printf("DELED\n");
      }

      gROOT->cd();
      //if(fmT) {delete fmT; fmT=NULL;}
      if(fmT_full) {delete fmT_full;fmT_full=NULL;}
      //printf("end\n");
      delete hProj;
      //printf("enddel\n");
    }

    fOutput.cd();
    g_GhettoFemto_mT_rstar_G.Write();
    g_GhettoFemto_mT_rcore_G.Write();
}


//flag =  1 usmani fit results
//        2 NLO fit results
void kstar_source_simple_1(const int& flag, const TString type="pp", int SEED=1){

  double hdr_size = 0;//0.75
  double hdr_slope = 0;//0.2

  double frac_proton_reso = 64.22;
  double m_proton_reso = 1362;
  double tau_proton_reso = 1.65;
  double frac_lambda_reso = 64.38;
  double m_lambda_reso = 1463;
  double tau_lambda_reso = 4.69;

  const double EtaCut = 0.8;
  const double TIMEOUT = 60*1;
  const double GLOB_TIMEOUT = 30*60;
  const unsigned Multiplicity=2;
  const double femto_region = 100;
  const unsigned target_yield = 1331*100;
  const unsigned NumCPU = 7;

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
  float th_kick = 0;
  float fixed_hdr = 1;
  float frag_beta = 0;
  float del_proton = 0;
  float del_proton_reso = 0;
  float del_lambda = 0;
  float del_lambda_reso = 0;
  int multiplicity = 2;

  TString BaseName = TString::Format("kstarSrc");

  if(type!="pp"&&type!="pL"){
    printf("WHAT IS THIS: %s\n", type.Data());
    return;
  }

  if(type=="pp"||type=="pL"){
    //for FEMTUM
    if(flag==1){
      //these values are for the ceca paper -> with usmani fit
      d_x = 0.176;
      d_x = d_y; d_x = d_z;
      h_x = 2.68;
      h_y = h_x;
      tau = 3.76;
    }
    else if(flag==2){
      //these values are for the ceca paper -> with NLO
      d_x = 0.288;
      d_x = d_y; d_x = d_z;
      h_x = 3.23;
      h_y = h_x;
      tau = 3.26;
    }
  }

  TREPNI Database(0);
  Database.SetSeed(SEED);
  std::vector<TreParticle*> ParticleList;
  ParticleList.push_back(Database.NewParticle("Proton"));
  ParticleList.push_back(Database.NewParticle("ProtonReso"));
  ParticleList.push_back(Database.NewParticle("Lambda"));
  ParticleList.push_back(Database.NewParticle("LambdaReso"));
  ParticleList.push_back(Database.NewParticle("Pion"));

  TString OutputFolderName;
  OutputFolderName = "FunWithCeca/kstar_source_simple_1";
  TString BaseFileName = TString::Format("%s/%s/%s_%s_flag%i_sd%i",
  GetFemtoOutputFolder(),OutputFolderName.Data(),BaseName.Data(),type.Data(),flag,SEED);
  TFile fOutput(BaseFileName+".root","recreate");

  DLM_Histo<float>* dlm_pT_eta_p = NULL;
  DLM_Histo<float>* dlm_pT_eta_L = NULL;
  //1 = Jaime (measured)
  //2 = full ALICE
  const int momdst_flag = 2;
  if(momdst_flag%100==1){
    dlm_pT_eta_p = GetPtEta(
     TString::Format("%s/Jaime/p_pT.root",GetCernBoxDimi()),
     TString::Format("%s/Jaime/ap_pT.root",GetCernBoxDimi()),
     "pTDist_after", "pTDist_after", EtaCut);
    dlm_pT_eta_L = GetPtEta(
       TString::Format("%s/Jaime/L_pT.root",GetCernBoxDimi()),
       TString::Format("%s/Jaime/aL_pT.root",GetCernBoxDimi()),
       "pTDist_after", "pTDist_after", EtaCut);
  }
  else if(momdst_flag%100==2){
    dlm_pT_eta_p = GetPtEta_13TeV(
      TString::Format("%s/CatsFiles/Source/CECA/proton_pT/p_dist_13TeV_ClassI.root",GetCernBoxDimi()),
      "Graph1D_y1", 500, 4050, EtaCut);
    dlm_pT_eta_L = GetPtEta_13TeV(
      TString::Format("%s/CatsFiles/Source/CECA/Lambda_pT/L_dist_13TeV_ClassI.root",GetCernBoxDimi()),
      "Graph1D_y1", 400, 8000, EtaCut);
  }
  else{
    printf("ERROR momdst_flag\n");
    return;
  }

  for(TreParticle* prt : ParticleList){
    if(prt->GetName()=="Proton"){
      prt->SetMass(Mass_p);
      prt->SetAbundance(100.-frac_proton_reso);
      prt->SetRadius(hdr_size);
      prt->SetRadiusSlope(hdr_slope);
      prt->SetDelayTau(del_proton);
      if(dlm_pT_eta_p) prt->SetPtEtaPhi(*dlm_pT_eta_p);
      else prt->SetPtPz(0.85*prt->GetMass(),0.85*prt->GetMass());
    }
    else if(prt->GetName()=="ProtonReso"){
      prt->SetMass(m_proton_reso);
      prt->SetAbundance(frac_proton_reso);
      prt->SetWidth(hbarc/tau_proton_reso);
      prt->SetDelayTau(del_proton_reso);

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
      prt->SetDelayTau(del_lambda);
      if(dlm_pT_eta_L) prt->SetPtEtaPhi(*dlm_pT_eta_L);
      else prt->SetPtPz(0.85*prt->GetMass(),0.85*prt->GetMass());
    }
    else if(prt->GetName()=="LambdaReso"){
      prt->SetMass(m_lambda_reso);
      prt->SetAbundance(frac_lambda_reso);
      prt->SetWidth(hbarc/tau_lambda_reso);
      prt->SetDelayTau(del_lambda_reso);

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
    return;
  }

  printf("Set up IVANA\n");

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

      Ivana.Ghetto_NumMomBins = 128;
      Ivana.Ghetto_MomMin = 0;
      Ivana.Ghetto_MomMax = 512;
  }
  else if(type=="pL"){
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

      Ivana.Ghetto_NumMomBins = 128;
      Ivana.Ghetto_MomMin = 0;
      Ivana.Ghetto_MomMax = 512;
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

  printf("    Total  Femto\n");
  printf("PP %6.2f%% %6.2f\n",TotPP*100.,FemtoPP*100.);
  printf("PR %6.2f%% %6.2f\n",TotPR*100.,FemtoPR*100.);
  printf("RP %6.2f%% %6.2f\n",TotRP*100.,FemtoRP*100.);
  printf("RR %6.2f%% %6.2f\n",TotRR*100.,FemtoRR*100.);

  TH1F* h_Ghetto_rstar = Convert_DlmHisto_TH1F(Ivana.Ghetto_rstar,"Ghetto_rstar");
  TH1F* h_Ghetto_kstar = Convert_DlmHisto_TH1F(Ivana.Ghetto_kstar,"Ghetto_kstar");
  Ivana.Ghetto_kstar_rstar->ComputeError();
  TH2F* h_Ghetto_kstar_rstar = Convert_DlmHisto_TH2F(Ivana.Ghetto_kstar_rstar,"Ghetto_kstar_rstar");
  TH2F* h_Ghetto_mT_rstar = Convert_DlmHisto_TH2F(Ivana.Ghetto_mT_rstar,"Ghetto_mT_rstar");
  TH2F* h_GhettoFemto_mT_rstar = Convert_DlmHisto_TH2F(Ivana.GhettoFemto_mT_rstar,"GhettoFemto_mT_rstar");
  TH2F* h_GhettoFemto_mT_kstar = Convert_DlmHisto_TH2F(Ivana.GhettoFemto_mT_kstar,"GhettoFemto_mT_kstar");
  TH2F* h_Ghetto_mT_costh = Convert_DlmHisto_TH2F(Ivana.Ghetto_mT_costh,"Ghetto_mT_costh");
  TH2F* h_GhettoSP_pT_th = Convert_DlmHisto_TH2F(Ivana.GhettoSP_pT_th,"GhettoSP_pT_th");
  TH1F* h_GhettoSP_pT_1 = Convert_DlmHisto_TH1F(Ivana.GhettoSP_pT_1,"GhettoSP_pT_1");
  TH1F* h_GhettoSP_pT_2 = Convert_DlmHisto_TH1F(Ivana.GhettoSP_pT_2,"GhettoSP_pT_2");
  Ivana.GhettoFemto_rstar->ComputeError();
  TH1F* h_GhettoFemto_rstar = Convert_DlmHisto_TH1F(Ivana.GhettoFemto_rstar,"GhettoFemto_rstar");

  fOutput.cd();

  h_Ghetto_kstar_rstar->Write();
  h_GhettoFemto_rstar->Write();

  TGraph gMean_mT;
  gMean_mT.SetName("gMean_mT");
  TGraph* gMean_kstar = new TGraph [Ivana.GhettoFemto_mT_rstar->GetNbins(0)];
  TH1F** hPP_kstar = new TH1F* [Ivana.GhettoFemto_mT_rstar->GetNbins(0)];
  TH1F** hPR_kstar = new TH1F* [Ivana.GhettoFemto_mT_rstar->GetNbins(0)];
  TH1F** hRP_kstar = new TH1F* [Ivana.GhettoFemto_mT_rstar->GetNbins(0)];
  TH1F** hRR_kstar = new TH1F* [Ivana.GhettoFemto_mT_rstar->GetNbins(0)];

  for(unsigned umT=0; umT<Ivana.GhettoFemto_mT_rstar->GetNbins(0); umT++){
    TH1D* hProj = h_GhettoFemto_mT_rstar->ProjectionY(TString::Format("hProj"),umT+1,umT+1);
    double Mean = hProj->GetMean();
    delete hProj;
    //double GFM = GaussFromMean(Mean);
    gMean_mT.SetPoint(umT,Ivana.GhettoFemto_mT_rstar->GetBinCenter(0,umT),Mean);

    //  //how the radius changes as a function of kstar
    gMean_kstar[umT].SetName(TString::Format("gMean_kstar_mT%.0f",Ivana.GhettoFemto_mT_rstar->GetBinCenter(0,umT)));
    if(umT<4) gMean_kstar[umT].SetLineColor(kRed+umT);
    else if(umT<8)  gMean_kstar[umT].SetLineColor(kGreen+umT-4);
    else if(umT<12) gMean_kstar[umT].SetLineColor(kBlue+umT-8);
    else if(umT<16) gMean_kstar[umT].SetLineColor(kMagenta+umT-12);
    else            gMean_kstar[umT].SetLineColor(kBlack);
    gMean_kstar[umT].SetLineWidth(3);
    gMean_kstar[umT].SetMarkerStyle(2);
    gMean_kstar[umT].SetMarkerColor(kViolet);
    gMean_kstar[umT].SetMarkerSize(1);
    TH2F* h_Ghetto_kstar_rstar_mT = new TH2F(TString::Format("h_Ghetto_kstar_rstar_mT%.0f",Ivana.GhettoFemto_mT_rstar->GetBinCenter(0,umT)),
            TString::Format("h_Ghetto_kstar_rstar_mT%.0f",Ivana.GhettoFemto_mT_rstar->GetBinCenter(0,umT)),
            Ivana.Ghetto_kstar_rstar_mT->GetNbins(0),Ivana.Ghetto_kstar_rstar_mT->GetLowEdge(0),Ivana.Ghetto_kstar_rstar_mT->GetUpEdge(0),
            Ivana.Ghetto_kstar_rstar_mT->GetNbins(1),Ivana.Ghetto_kstar_rstar_mT->GetLowEdge(1),Ivana.Ghetto_kstar_rstar_mT->GetUpEdge(1));
    hPP_kstar[umT] = new TH1F(TString::Format("hPP_kstar_mT%.0f",Ivana.GhettoFemto_mT_rstar->GetBinCenter(0,umT)),
            TString::Format("hPP_kstar_mT%.0f",Ivana.GhettoFemto_mT_rstar->GetBinCenter(0,umT)),
            Ivana.Ghetto_kstar_rstar_mT->GetNbins(0),Ivana.Ghetto_kstar_rstar_mT->GetLowEdge(0),Ivana.Ghetto_kstar_rstar_mT->GetUpEdge(0)
            );
    hPR_kstar[umT] = new TH1F(TString::Format("hPR_kstar_mT%.0f",Ivana.GhettoFemto_mT_rstar->GetBinCenter(0,umT)),
            TString::Format("hPR_kstar_mT%.0f",Ivana.GhettoFemto_mT_rstar->GetBinCenter(0,umT)),
            Ivana.Ghetto_kstar_rstar_mT->GetNbins(0),Ivana.Ghetto_kstar_rstar_mT->GetLowEdge(0),Ivana.Ghetto_kstar_rstar_mT->GetUpEdge(0)
            );
    hRP_kstar[umT] = new TH1F(TString::Format("hRP_kstar_mT%.0f",Ivana.GhettoFemto_mT_rstar->GetBinCenter(0,umT)),
            TString::Format("hRP_kstar_mT%.0f",Ivana.GhettoFemto_mT_rstar->GetBinCenter(0,umT)),
            Ivana.Ghetto_kstar_rstar_mT->GetNbins(0),Ivana.Ghetto_kstar_rstar_mT->GetLowEdge(0),Ivana.Ghetto_kstar_rstar_mT->GetUpEdge(0)
            );
    hRR_kstar[umT] = new TH1F(TString::Format("hRR_kstar_mT%.0f",Ivana.GhettoFemto_mT_rstar->GetBinCenter(0,umT)),
            TString::Format("hRR_kstar_mT%.0f",Ivana.GhettoFemto_mT_rstar->GetBinCenter(0,umT)),
            Ivana.Ghetto_kstar_rstar_mT->GetNbins(0),Ivana.Ghetto_kstar_rstar_mT->GetLowEdge(0),Ivana.Ghetto_kstar_rstar_mT->GetUpEdge(0)
            );
    for(unsigned ukstar=0; ukstar<Ivana.Ghetto_kstar_rstar_mT->GetNbins(0); ukstar++){
      for(unsigned urstar=0; urstar<Ivana.Ghetto_kstar_rstar_mT->GetNbins(1); urstar++){
        h_Ghetto_kstar_rstar_mT->SetBinContent(ukstar+1,urstar+1,Ivana.Ghetto_kstar_rstar_mT->GetBinContent(ukstar,urstar,umT));
        //printf("mt%u ks%u rs%u --> %.2f\n",umT,ukstar,urstar,Ivana.Ghetto_kstar_rstar_mT->GetBinContent(ukstar,urstar,umT));
        //usleep(1e3);
      }
      double TOT_PAIRS = Ivana.Ghetto_kstar_reso_mT->GetBinContent(ukstar,0,umT);
      TOT_PAIRS += Ivana.Ghetto_kstar_reso_mT->GetBinContent(ukstar,1,umT);
      TOT_PAIRS += Ivana.Ghetto_kstar_reso_mT->GetBinContent(ukstar,2,umT);
      TOT_PAIRS += Ivana.Ghetto_kstar_reso_mT->GetBinContent(ukstar,3,umT);
      if(!TOT_PAIRS) TOT_PAIRS = 1;
      hPP_kstar[umT]->SetBinContent(ukstar+1,double(Ivana.Ghetto_kstar_reso_mT->GetBinContent(ukstar,0,umT))/TOT_PAIRS);
      hPR_kstar[umT]->SetBinContent(ukstar+1,double(Ivana.Ghetto_kstar_reso_mT->GetBinContent(ukstar,1,umT))/TOT_PAIRS);
      hRP_kstar[umT]->SetBinContent(ukstar+1,double(Ivana.Ghetto_kstar_reso_mT->GetBinContent(ukstar,2,umT))/TOT_PAIRS);
      hRR_kstar[umT]->SetBinContent(ukstar+1,double(Ivana.Ghetto_kstar_reso_mT->GetBinContent(ukstar,3,umT))/TOT_PAIRS);
    }
    for(unsigned ukstar=0; ukstar<Ivana.Ghetto_kstar_rstar_mT->GetNbins(0); ukstar++){
      hProj = h_Ghetto_kstar_rstar_mT->ProjectionY(TString::Format("hProj"),ukstar+1,ukstar+1);
      Mean = hProj->GetMean();
      delete hProj;
      //GFM = GaussFromMean(Mean);
      gMean_kstar[umT].SetPoint(ukstar,Ivana.Ghetto_kstar_rstar_mT->GetBinCenter(0,ukstar),Mean);
    }

    fOutput.cd();
    h_Ghetto_kstar_rstar_mT->Write();
    gMean_kstar[umT].Write();
    hPP_kstar[umT]->Write();
    hPR_kstar[umT]->Write();
    hRP_kstar[umT]->Write();
    hRR_kstar[umT]->Write();
    //delete h_Ghetto_kstar_rstar_mT;
  }

  gMean_mT.Write();

  for(unsigned umT=0; umT<Ivana.GhettoFemto_mT_rstar->GetNbins(0); umT++){
    delete hPP_kstar[umT];
    delete hPR_kstar[umT];
    delete hRP_kstar[umT];
    delete hRR_kstar[umT];
  }
  delete [] hRR_kstar;

  delete [] gMean_kstar;
}


void kstar_source_ck_TEST1(){

  const double SourceSize = 1.0;
  const double kMin = 0;
  const double kMax = 320;
  const unsigned NumMomBins = 80;
  const double KstarStrength = 0.9;

  //1D histo
  DLM_Histo<float> HistoSourceR;
  //2D histo (k,r in that order)
  DLM_Histo<float> HistoSourceRK;
  //the same as above, but the k relation is constant
  //qa to reproduce HistoSourceR
  DLM_Histo<float> HistoSourceRKconst;

  const unsigned NumRadBins = 8192;
  const double rMin = 0;
  const double rMax = 96;

  HistoSourceR.SetUp(1);
  HistoSourceR.SetUp(0,NumRadBins,rMin,rMax);
  HistoSourceR.Initialize();

  HistoSourceRK.SetUp(2);
  HistoSourceRK.SetUp(0,NumMomBins,kMin,kMax);
  HistoSourceRK.SetUp(1,NumRadBins,rMin,rMax);
  HistoSourceRK.Initialize();

  HistoSourceRKconst.SetUp(2);
  HistoSourceRKconst.SetUp(0,NumMomBins,kMin,kMax);
  HistoSourceRKconst.SetUp(1,NumRadBins,rMin,rMax);
  HistoSourceRKconst.Initialize();

  double PARS[3];

  for(unsigned uRad=0; uRad<NumRadBins; uRad++){
    double rad = HistoSourceR.GetBinCenter(0,uRad);
    PARS[1] = rad;
    PARS[3] = SourceSize;
    double Sr = GaussSource(PARS);
    HistoSourceR.SetBinContent(uRad,Sr);
    for(unsigned uMom=0; uMom<NumMomBins; uMom++){
      double Kstar = HistoSourceRK.GetBinCenter(0,uMom);
      HistoSourceRKconst.SetBinContent(uMom,uRad,Sr);
      //with this relation, the original source size is restored at 100 MeV
      //while at kstar = 0 we have KstarStrength percent of the sources
      PARS[3] = KstarStrength*SourceSize + (SourceSize - KstarStrength*SourceSize)*Kstar/100.;
      SourceSize*Kstar/1000.;
      double SrK = GaussSource(PARS);
      HistoSourceRK.SetBinContent(uMom,uRad,SrK);
    }
  }

  TString BaseFileName = TString::Format("%s/FunWithCeca/kstar_source_ck_TEST1/test1.root",GetFemtoOutputFolder());
  TFile fOutput(BaseFileName,"recreate");

  TH1F* hHistoSourceR = Convert_DlmHisto_TH1F(&HistoSourceR,"HistoSourceR");
  TH2F* hHistoSourceRK = Convert_DlmHisto_TH2F(&HistoSourceRK,"HistoSourceRK");
  TH2F* hHistoSourceRKconst = Convert_DlmHisto_TH2F(&HistoSourceRKconst,"HistoSourceRKconst");

  fOutput.cd();
  hHistoSourceR->Write();
  hHistoSourceRK->Write();
  hHistoSourceRKconst->Write();

  CATS Kitty_Gr;
  CATS Kitty_HGr;
  CATS Kitty_HGrK;
  CATS Kitty_HGrKconst;

  Kitty_Gr.SetMomBins(NumMomBins,kMin,kMax);

  Kitty_HGr.SetMomBins(NumMomBins,kMin,kMax);
  Kitty_HGrK.SetMomBins(NumMomBins,kMin,kMax);
  Kitty_HGrKconst.SetMomBins(NumMomBins,kMin,kMax);

  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles/",GetCernBoxDimi()));

  AnalysisObject.SetUpCats_pp(Kitty_Gr,"AV18","Gauss",0,0);
  AnalysisObject.SetUpCats_pp(Kitty_HGr,"AV18","",0,0);
  AnalysisObject.SetUpCats_pp(Kitty_HGrK,"AV18","",0,0);
  AnalysisObject.SetUpCats_pp(Kitty_HGrKconst,"AV18","",0,0);

  Kitty_Gr.SetAnaSource(0,SourceSize);

  DLM_HistoSource src_HGr(HistoSourceR);
  Kitty_HGr.SetAnaSource(CatsSourceForwarder, &src_HGr, 0);
  Kitty_HGr.SetUseAnalyticSource(true);
  Kitty_HGr.SetAutoNormSource(false);
  Kitty_HGr.SetNormalizedSource(true);

  DLM_HistoSource src_HGrKconst(HistoSourceRKconst);
  Kitty_HGrKconst.SetAnaSource(CatsSourceForwarder, &src_HGrKconst, 0);
  Kitty_HGrKconst.SetUseAnalyticSource(true);
  Kitty_HGrKconst.SetAutoNormSource(false);
  Kitty_HGrKconst.SetNormalizedSource(true);
  Kitty_HGrKconst.SetMomentumDependentSource(true);

  DLM_HistoSource src_HGrK(HistoSourceRK);
  Kitty_HGrK.SetAnaSource(CatsSourceForwarder, &src_HGrK, 0);
  Kitty_HGrK.SetUseAnalyticSource(true);
  Kitty_HGrK.SetAutoNormSource(false);
  Kitty_HGrK.SetNormalizedSource(true);
  Kitty_HGrK.SetMomentumDependentSource(true);

  Kitty_Gr.KillTheCat();
  Kitty_HGr.KillTheCat();
  Kitty_HGrK.KillTheCat();
  Kitty_HGrKconst.KillTheCat();

  TGraph Ck_Gr;
  TGraph Ck_HGr;
  TGraph Ck_HGrK;
  TGraph Ck_HGrKconst;

  Ck_Gr.SetName("Ck_Gr");
  Ck_HGr.SetName("Ck_HGr");
  Ck_HGrK.SetName("Ck_HGrK");
  Ck_HGrKconst.SetName("Ck_HGrKconst");

  for(unsigned uMom=0; uMom<NumMomBins; uMom++){
    double kstar = Kitty_Gr.GetMomentum(uMom);
    Ck_Gr.SetPoint(uMom,kstar,Kitty_Gr.GetCorrFun(uMom));
    Ck_HGr.SetPoint(uMom,kstar,Kitty_HGr.GetCorrFun(uMom));
    Ck_HGrK.SetPoint(uMom,kstar,Kitty_HGrK.GetCorrFun(uMom));
    Ck_HGrKconst.SetPoint(uMom,kstar,Kitty_HGrKconst.GetCorrFun(uMom));
  }

  Ck_Gr.Write();
  Ck_HGr.Write();
  Ck_HGrK.Write();
  Ck_HGrKconst.Write();



}

//study the distribution of primordials as mT and kstar
void CECA_primoridal_disto(int SEED, unsigned NUM_CPU){

  double EtaCut = 0.8;
  DLM_Histo<float>* dlm_pT_eta_p;
  dlm_pT_eta_p = GetPtEta_13TeV(
      TString::Format("%s/CatsFiles/Source/CECA/proton_pT/p_dist_13TeV_ClassI.root",GetCernBoxDimi()),
      "Graph1D_y1", 500, 4050, EtaCut);

  double HadronSize = 0;//0.75
  double HadronSlope = 0;//0.2
  const bool PROTON_RESO = true;
  const bool EQUALIZE_TAU = true;
  const double TIMEOUT = 30;
  const double GLOB_TIMEOUT = 2*60;
  const unsigned Multiplicity=2;
  const double femto_region = 100;
  const unsigned target_yield = 1000*1000;
  TString BaseName = TString::Format("Eta%.1f",EtaCut);

  //we run to either reproduce the core of 0.97,
  //or the upper limit of reff = 1.06+0.04
  //this leads to a 10% difference in the SP core source
  double rSP_core=0;
  double rSP_dispZ=0;
  double rSP_hadr=0;
  double rSP_hadrZ=0;
  double rSP_hflc=0;
  double rSP_tau=0;
  bool tau_prp = true;
  double rSP_tflc=0;
  double rSP_ThK=0;
  //default is true
  bool rSP_FixedHadr = true;
  //default is 0
  float rSP_FragBeta = 0;

  //this so happens to be exactly the mean of the results from the ceca paper
  //the errors are the half of the difference of the two ceca results + in ^2 the two errors from the paper
  rSP_core = (0.176+0.288)*0.5;
  rSP_dispZ = rSP_core;
  rSP_hadr = (2.68+3.23)*0.5;
  rSP_tau = (3.76+3.26)*0.5;

  TREPNI Database(0);
  Database.SetSeed(11);
  std::vector<TreParticle*> ParticleList;
  ParticleList.push_back(Database.NewParticle("Proton"));
  ParticleList.push_back(Database.NewParticle("Pion"));
  ParticleList.push_back(Database.NewParticle("ProtonReso"));

  TString OutputFolderName;
  OutputFolderName = TString::Format("%s/FunWithCeca/CECA_primoridal_disto",GetFemtoOutputFolder());

  for(TreParticle* prt : ParticleList){
    if(prt->GetName()=="Proton"){
      prt->SetMass(Mass_p);
      prt->SetAbundance(35.78+64.22*(!PROTON_RESO));
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);
      prt->SetPtEtaPhi(*dlm_pT_eta_p);
    }
    else if(prt->GetName()=="Pion"){
      prt->SetMass(Mass_pic);
      prt->SetAbundance(0);
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);
    }
    else if(prt->GetName()=="ProtonReso"){
      prt->SetMass(1362);
      prt->SetAbundance(64.22*PROTON_RESO);
      prt->SetWidth(hbarc/1.65);
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);
      prt->SetPtEtaPhi(*dlm_pT_eta_p);

      prt->NewDecay();
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Proton"));
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
      prt->GetDecay(0)->SetBranching(100);
    }
  }

  std::vector<std::string> ListOfParticles;

  ListOfParticles.push_back("Proton");
  ListOfParticles.push_back("Proton");

  CECA Ivana(Database,ListOfParticles);

  Ivana.SetDisplacementZ(rSP_dispZ);
  Ivana.SetDisplacementT(rSP_core);
  Ivana.SetHadronizationZ(rSP_hadrZ);//0
  Ivana.SetHadronizationT(rSP_hadr);
  Ivana.SetHadrFluctuation(rSP_hflc);
  Ivana.SetTau(rSP_tau,tau_prp);
  Ivana.SetTauFluct(rSP_tflc);
  Ivana.SetThermalKick(rSP_ThK);
  Ivana.SetFixedHadr(rSP_FixedHadr);
  Ivana.SetFragmentBeta(rSP_FragBeta);

  Ivana.SetTargetStatistics(target_yield);
  Ivana.SetEventMult(Multiplicity);
  Ivana.SetSourceDim(2);
  Ivana.SetDebugMode(true);
  Ivana.SetThreadTimeout(TIMEOUT);
  Ivana.SetGlobalTimeout(GLOB_TIMEOUT);
  Ivana.EqualizeFsiTime(EQUALIZE_TAU);
  Ivana.SetFemtoRegion(femto_region);
  Ivana.GHETTO_EVENT = true;

  if(NUM_CPU>1){
    Ivana.SetDebugMode(true);
    for(unsigned uTh=0; uTh<NUM_CPU; uTh++){
      Ivana.SetSeed(uTh,SEED*(NUM_CPU)+uTh);
    }
  }
  else{
    Ivana.SetDebugMode(false);
    Ivana.SetSeed(0,SEED);
  }

  Ivana.GoBabyGo(NUM_CPU);

  Ivana.GhettoFemtoPrimordial_mT_kstar->ComputeError();
  TH2F* h_GhettoFemtoPrimordial_mT_kstar = Convert_DlmHisto_TH2F(Ivana.GhettoFemtoPrimordial_mT_kstar,"GhettoFemtoPrimordial_mT_kstar");
  Ivana.GhettoFemto_mT_kstar->ComputeError();
  TH2F* h_GhettoFemto_mT_kstar = Convert_DlmHisto_TH2F(Ivana.GhettoFemto_mT_kstar,"GhettoFemto_mT_kstar");

  TH2F* h_Frac_Prim = (TH2F*)h_GhettoFemtoPrimordial_mT_kstar->Clone("h_Frac_Prim");
  h_Frac_Prim->Divide(h_GhettoFemto_mT_kstar);

  // Save TNtuple to ROOT file
  TString pair_root_name(TString::Format(OutputFolderName+"/prim_info_S%i_N%i.root",SEED,NUM_CPU).Data());
  TFile* outFile = new TFile(pair_root_name, "RECREATE");
  h_GhettoFemtoPrimordial_mT_kstar->Write();
  h_GhettoFemto_mT_kstar->Write();
  h_Frac_Prim->Write();
  outFile->Close();


  // Clean up memory
  delete outFile;
  delete dlm_pT_eta_p;
}



//a first attempt to simulate CECA events for pp, to be used by the AI people
void CECA_sim_for_AI_pp_v0(int SEED, unsigned NUM_CPU){

  double EtaCut = 0.8;
  DLM_Histo<float>* dlm_pT_eta_p;
  dlm_pT_eta_p = GetPtEta_13TeV(
      TString::Format("%s/CatsFiles/Source/CECA/proton_pT/p_dist_13TeV_ClassI.root",GetCernBoxDimi()),
      "Graph1D_y1", 500, 4050, EtaCut);

  double HadronSize = 0;//0.75
  double HadronSlope = 0;//0.2
  const bool PROTON_RESO = true;
  const bool EQUALIZE_TAU = true;
  const double TIMEOUT = 30;
  const double GLOB_TIMEOUT = 30*60;
  const unsigned Multiplicity=2;
  const double femto_region = 100;
  const unsigned target_yield = 1000*1000;
  TString BaseName = TString::Format("Eta%.1f",EtaCut);

  //we run to either reproduce the core of 0.97,
  //or the upper limit of reff = 1.06+0.04
  //this leads to a 10% difference in the SP core source
  double rSP_core=0;
  double rSP_dispZ=0;
  double rSP_hadr=0;
  double rSP_hadrZ=0;
  double rSP_hflc=0;
  double rSP_tau=0;
  bool tau_prp = true;
  double rSP_tflc=0;
  double rSP_ThK=0;
  //default is true
  bool rSP_FixedHadr = true;
  //default is 0
  float rSP_FragBeta = 0;

  //this so happens to be exactly the mean of the results from the ceca paper
  //the errors are the half of the difference of the two ceca results + in ^2 the two errors from the paper
  rSP_core = (0.176+0.288)*0.5;
  double rSP_core_SIG = fabs(0.176-0.288)*0.5;
  rSP_core_SIG = sqrt(pow(rSP_core_SIG,2.)+pow(0.013,2.)+pow(0.005,2.));

  rSP_dispZ = rSP_core;
  double rSP_dispZ_SIG = rSP_core_SIG;

  rSP_hadr = (2.68+3.23)*0.5;
  double rSP_hadr_SIG = fabs(2.68-3.23)*0.5;
  rSP_hadr_SIG = sqrt(pow(rSP_hadr_SIG,2.)+pow(0.175,2.)+pow(0.05,2.));

  rSP_tau = (3.76+3.26)*0.5;
  double rSP_tau_SIG = fabs(3.76-3.26)*0.5;
  rSP_tau_SIG = sqrt(pow(rSP_tau_SIG,2.)+pow(0.1,2.)+pow(0.04,2.));

  printf("rSP_core = %.4e +/- %.4e\n",rSP_core, rSP_core_SIG);
  printf("rSP_hadr = %.4e +/- %.4e\n",rSP_hadr, rSP_hadr_SIG);
  printf("rSP_tau = %.4e +/- %.4e\n",rSP_tau, rSP_tau_SIG);


  TREPNI Database(0);
  Database.SetSeed(11);
  std::vector<TreParticle*> ParticleList;
  ParticleList.push_back(Database.NewParticle("Proton"));
  ParticleList.push_back(Database.NewParticle("Pion"));
  ParticleList.push_back(Database.NewParticle("ProtonReso"));


  TString OutputFolderName;
  OutputFolderName = TString::Format("%s/FunWithCeca/CECA_sim_for_AI_pp",GetFemtoOutputFolder());
  //TFile fOutput(BaseFileName+".root","recreate");


  for(TreParticle* prt : ParticleList){
    if(prt->GetName()=="Proton"){
      prt->SetMass(Mass_p);
      prt->SetAbundance(35.78+64.22*(!PROTON_RESO));
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);
      prt->SetPtEtaPhi(*dlm_pT_eta_p);
    }
    else if(prt->GetName()=="Pion"){
      prt->SetMass(Mass_pic);
      prt->SetAbundance(0);
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);
    }
    else if(prt->GetName()=="ProtonReso"){
      prt->SetMass(1362);
      prt->SetAbundance(64.22*PROTON_RESO);
      prt->SetWidth(hbarc/1.65);
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);
      prt->SetPtEtaPhi(*dlm_pT_eta_p);

      prt->NewDecay();
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Proton"));
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
      prt->GetDecay(0)->SetBranching(100);
    }
  }

  std::vector<std::string> ListOfParticles;

  ListOfParticles.push_back("Proton");
  ListOfParticles.push_back("Proton");

  CECA Ivana(Database,ListOfParticles);

  Ivana.SetDisplacementZ(rSP_dispZ);
  Ivana.SetDisplacementT(rSP_core);
  Ivana.SetDisplacementEbe(rSP_core_SIG/rSP_core);
  Ivana.SetHadronizationZ(rSP_hadrZ);//0
  Ivana.SetHadronizationT(rSP_hadr);
  Ivana.SetHadronizationEbe(rSP_hadr_SIG/rSP_hadr);
  Ivana.SetHadrFluctuation(rSP_hflc);
  Ivana.SetTau(rSP_tau,tau_prp);
  Ivana.SetTauEbe(rSP_tau_SIG/rSP_tau);
  Ivana.SetTauFluct(rSP_tflc);
  Ivana.SetThermalKick(rSP_ThK);
  Ivana.SetFixedHadr(rSP_FixedHadr);
  Ivana.SetFragmentBeta(rSP_FragBeta);

  Ivana.SetTargetStatistics(target_yield);
  Ivana.SetEventMult(Multiplicity);
  Ivana.SetSourceDim(2);
  Ivana.SetDebugMode(true);
  Ivana.SetThreadTimeout(TIMEOUT);
  Ivana.SetGlobalTimeout(GLOB_TIMEOUT);
  Ivana.EqualizeFsiTime(EQUALIZE_TAU);
  Ivana.SetFemtoRegion(femto_region);
  Ivana.GHETTO_EVENT = true;
  std::string pair_file_name(TString::Format(OutputFolderName+"/pair_info_S%i_N%i.data",SEED,NUM_CPU).Data());
  Ivana.SetExportPairs(2,pair_file_name);

  if(NUM_CPU>1){
    Ivana.SetDebugMode(true);
    for(unsigned uTh=0; uTh<NUM_CPU; uTh++){
      Ivana.SetSeed(uTh,SEED*(NUM_CPU)+uTh);
    }
  }
  else{
    Ivana.SetDebugMode(false);
    Ivana.SetSeed(0,SEED);
  }

  Ivana.GoBabyGo(NUM_CPU);

  // Open the text file
  std::ifstream inFile(pair_file_name);

  if (!inFile.is_open()) {
      std::cerr << "Error opening input file." << std::endl;
      return;
  }

  // Create TNtuple
  TNtuple* ntuple = new TNtuple("CECA_Data", "CECA_Data", "kstar:rstar:mT:rd:h:tau");

  // Skip the first line
  std::string line;
  std::getline(inFile, line);

  // Read and process the remaining lines
  while (std::getline(inFile, line)) {
      std::istringstream iss(line);
      float kstar, rstar, mT, rd, h, tau;

      // Read values from the line
      iss >> kstar >> rstar >> mT >> rd >> h >> tau;

      // Fill the TNtuple with the read values
      ntuple->Fill(kstar, rstar, mT, rd, h, tau);
  }

  // Close the input file
  inFile.close();

  // Save TNtuple to ROOT file
  TString pair_root_name(TString::Format(OutputFolderName+"/pair_info_S%i_N%i.root",SEED,NUM_CPU).Data());
  TFile* outFile = new TFile(pair_root_name, "RECREATE");
  ntuple->Write();
  outFile->Close();

  // Clean up memory
  delete outFile;
  delete ntuple;

//rd,h,tau --> get a list of pairs with r* and mT
//sample rd,h,tau randomly from a gauss
//for the output, we need not a histo, but an ASCII or binary file with:
//kstar,rstar,mT,rd,h,tau

  delete dlm_pT_eta_p;
}


//a second attempt to simulate CECA events for pp, to be used by the AI people
//now, we have a fixed amount of yield per configurations simulated
//NUM_CONFIGS -> random configurations that we simulate
//the flag, so far is:
//0 -> clean output into a root file
//2 -> also ASCII files, with rd, hT, tau
//do NOT use other flags, as the output will be even more
void CECA_sim_for_AI_pp_v1(int SEED, unsigned NUM_CPU, unsigned NUM_CONFIGS, unsigned K_YIELD_PER_CONFIG, int flag, TString OutputFolderName){

  unsigned YIELD_PER_CONFIG = K_YIELD_PER_CONFIG*1000;
  double EtaCut = 0.8;
  DLM_Histo<float>* dlm_pT_eta_p;
  dlm_pT_eta_p = GetPtEta_13TeV(
      TString::Format("%s/CatsFiles/Source/CECA/proton_pT/p_dist_13TeV_ClassI.root",GetCernBoxDimi()),
      "Graph1D_y1", 500, 4050, EtaCut);

  double HadronSize = 0;//0.75
  double HadronSlope = 0;//0.2
  const bool PROTON_RESO = true;
  const bool EQUALIZE_TAU = true;
  const double TIMEOUT = 15;
  const double GLOB_TIMEOUT = 30*60;
  const unsigned Multiplicity=2;
  const double femto_region = 100;
  const unsigned target_yield = YIELD_PER_CONFIG;
  TString BaseName = TString::Format("Eta%.1f",EtaCut);

  //we run to either reproduce the core of 0.97,
  //or the upper limit of reff = 1.06+0.04
  //this leads to a 10% difference in the SP core source
  double rSP_core=0;
  double rSP_dispZ=0;
  double rSP_hadr=0;
  double rSP_hadrZ=0;
  double rSP_hflc=0;
  double rSP_tau=0;
  bool tau_prp = true;
  double rSP_tflc=0;
  double rSP_ThK=0;
  //default is true
  bool rSP_FixedHadr = true;
  //default is 0
  float rSP_FragBeta = 0;

  TRandom3 rangen(SEED);

  //TString OutputFolderName;
  //OutputFolderName = TString::Format("%s/FunWithCeca/CECA_sim_for_AI_pp_v3",GetFemtoOutputFolder());
  //TFile fOutput(BaseFileName+".root","recreate");

  // Save TNtuple to ROOT file
  TString pair_root_name(TString::Format(OutputFolderName+"/pair_info_S%i_N%i_YPC%iK.root",SEED,NUM_CPU,K_YIELD_PER_CONFIG).Data());
  TFile* outFile = new TFile(pair_root_name, "RECREATE");

  // Create TNtuple
  TNtuple* ntuple = new TNtuple("CECA_Data", "CECA_Data", "config_id:kstar:rstar:mT:rd:h:tau");


  for(unsigned uConfig=0; uConfig<NUM_CONFIGS; uConfig++){
    //this so happens to be exactly the mean of the results from the ceca paper
    //the errors are the half of the difference of the two ceca results + in ^2 the two errors from the paper
    rSP_core = (0.176+0.288)*0.5;
    double rSP_core_AVG = rSP_core;
    double rSP_core_SIG = fabs(0.176-0.288)*0.5;
    rSP_core_SIG = sqrt(pow(rSP_core_SIG,2.)+pow(0.013,2.)+pow(0.005,2.));
    rSP_core = rangen.Uniform(rSP_core_AVG-rSP_core_SIG,rSP_core_AVG+rSP_core_SIG);

    rSP_dispZ = rSP_core;
    double rSP_dispZ_SIG = rSP_core_SIG;

    rSP_hadr = (2.68+3.23)*0.5;
    double rSP_hadr_AVG = rSP_hadr;
    double rSP_hadr_SIG = fabs(2.68-3.23)*0.5;
    rSP_hadr_SIG = sqrt(pow(rSP_hadr_SIG,2.)+pow(0.175,2.)+pow(0.05,2.));
    rSP_hadr = rangen.Uniform(rSP_hadr_AVG-rSP_hadr_SIG,rSP_hadr_AVG+rSP_hadr_SIG);


    rSP_tau = (3.76+3.26)*0.5;
    double rSP_tau_AVG = rSP_tau;
    double rSP_tau_SIG = fabs(3.76-3.26)*0.5;
    rSP_tau_SIG = sqrt(pow(rSP_tau_SIG,2.)+pow(0.1,2.)+pow(0.04,2.));
    rSP_tau = rangen.Uniform(rSP_hadr_AVG-rSP_tau_SIG,rSP_tau_AVG+rSP_tau_SIG);

    printf("rSP_core = %.4e +/- %.4e\n",rSP_core, rSP_core_SIG);
    printf("rSP_hadr = %.4e +/- %.4e\n",rSP_hadr, rSP_hadr_SIG);
    printf("rSP_tau = %.4e +/- %.4e\n",rSP_tau, rSP_tau_SIG);   

    TREPNI Database(0);
    Database.SetSeed(11);
    std::vector<TreParticle*> ParticleList;
    ParticleList.push_back(Database.NewParticle("Proton"));
    ParticleList.push_back(Database.NewParticle("Pion"));
    ParticleList.push_back(Database.NewParticle("ProtonReso"));    

    for(TreParticle* prt : ParticleList){
      if(prt->GetName()=="Proton"){
        prt->SetMass(Mass_p);
        prt->SetAbundance(35.78+64.22*(!PROTON_RESO));
        prt->SetRadius(HadronSize);
        prt->SetRadiusSlope(HadronSlope);
        prt->SetPtEtaPhi(*dlm_pT_eta_p);
      }
      else if(prt->GetName()=="Pion"){
        prt->SetMass(Mass_pic);
        prt->SetAbundance(0);
        prt->SetRadius(HadronSize);
        prt->SetRadiusSlope(HadronSlope);
      }
      else if(prt->GetName()=="ProtonReso"){
        prt->SetMass(1362);
        prt->SetAbundance(64.22*PROTON_RESO);
        prt->SetWidth(hbarc/1.65);
        prt->SetRadius(HadronSize);
        prt->SetRadiusSlope(HadronSlope);
        prt->SetPtEtaPhi(*dlm_pT_eta_p);

        prt->NewDecay();
        prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Proton"));
        prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
        prt->GetDecay(0)->SetBranching(100);
      }
    }

    std::vector<std::string> ListOfParticles;

    ListOfParticles.push_back("Proton");
    ListOfParticles.push_back("Proton");

    CECA Ivana(Database,ListOfParticles);

    Ivana.SetDisplacementZ(rSP_dispZ);
    Ivana.SetDisplacementT(rSP_core);
    //Ivana.SetDisplacementEbe(rSP_core_SIG/rSP_core);
    Ivana.SetHadronizationZ(rSP_hadrZ);//0
    Ivana.SetHadronizationT(rSP_hadr);
    //Ivana.SetHadronizationEbe(rSP_hadr_SIG/rSP_hadr);
    Ivana.SetHadrFluctuation(rSP_hflc);
    Ivana.SetTau(rSP_tau,tau_prp);
    //Ivana.SetTauEbe(rSP_tau_SIG/rSP_tau);
    Ivana.SetTauFluct(rSP_tflc);
    Ivana.SetThermalKick(rSP_ThK);
    Ivana.SetFixedHadr(rSP_FixedHadr);
    Ivana.SetFragmentBeta(rSP_FragBeta);

    Ivana.SetTargetStatistics(target_yield);
    Ivana.SetEventMult(Multiplicity);
    Ivana.SetSourceDim(2);
    Ivana.SetDebugMode(true);
    Ivana.SetThreadTimeout(TIMEOUT);
    Ivana.SetGlobalTimeout(GLOB_TIMEOUT);
    Ivana.EqualizeFsiTime(EQUALIZE_TAU);
    Ivana.SetFemtoRegion(femto_region);
    Ivana.GHETTO_EVENT = true;
    std::string pair_file_name(TString::Format(OutputFolderName+"/pair_info_S%i_N%i_YPC%iK_C%u.data",SEED,NUM_CPU,K_YIELD_PER_CONFIG,uConfig).Data());
    Ivana.SetExportPairs(flag%10,pair_file_name);

    if(NUM_CPU>1){
      Ivana.SetDebugMode(true);
      for(unsigned uTh=0; uTh<NUM_CPU; uTh++){
        Ivana.SetSeed(uTh,SEED*(NUM_CPU)+uTh);
      }
    }
    else{
      Ivana.SetDebugMode(false);
      Ivana.SetSeed(0,SEED);
    }

    Ivana.GoBabyGo(NUM_CPU);

    // Open the text file
    std::ifstream inFile(pair_file_name);

    if (!inFile.is_open()) {
        std::cerr << "Error opening input file." << std::endl;
        return;
    }

    // Skip the first line
    std::string line;
    std::getline(inFile, line);

    outFile->cd();
    // Read and process the remaining lines
    while (std::getline(inFile, line)) {
        std::istringstream iss(line);
        float kstar, rstar, mT, rd, h, tau;

        // Read values from the line
        iss >> kstar >> rstar >> mT >> rd >> h >> tau;

        // Fill the TNtuple with the read values
        ntuple->Fill(uConfig, kstar, rstar, mT, rd, h, tau);
    }

    // Close the input file
    inFile.close();


  }//uConfig

  outFile->cd();
  ntuple->Write();

  // Clean up memory
  delete ntuple;
  delete outFile;
  
//rd,h,tau --> get a list of pairs with r* and mT
//sample rd,h,tau randomly from a gauss
//for the output, we need not a histo, but an ASCII or binary file with:
//kstar,rstar,mT,rd,h,tau

  delete dlm_pT_eta_p;
}

//a third attempt to simulate CECA events for pp, to be used by the AI/Carla people
//now, we have a fixed amount of yield per configurations simulated, plus flexibility of the set up (phase space to explore)
//NUM_CONFIGS -> random configurations that we simulate
//the flag, so far is:
//0 -> clean output into a root file
//2 -> also ASCII files, with rd, hT, tau
//do NOT use other flags, as the output will be even more
void CECA_sim_for_AI_pp_v2(int SEED, unsigned NUM_CPU, unsigned NUM_CONFIGS, unsigned K_YIELD_PER_CONFIG, 
  float min_rd, float max_rd, float min_ht, float max_ht, float min_tau, float max_tau,
  int flag, TString OutputFolderName){

  unsigned YIELD_PER_CONFIG = K_YIELD_PER_CONFIG*1000;
  double EtaCut = 0.8;
  DLM_Histo<float>* dlm_pT_eta_p;
  dlm_pT_eta_p = GetPtEta_13TeV(
      TString::Format("%s/CatsFiles/Source/CECA/proton_pT/p_dist_13TeV_ClassI.root",GetCernBoxDimi()),
      "Graph1D_y1", 500, 4050, EtaCut);

  double HadronSize = 0;//0.75
  double HadronSlope = 0;//0.2
  const bool PROTON_RESO = true;
  const bool EQUALIZE_TAU = true;
  const double TIMEOUT = 30;
  const double GLOB_TIMEOUT = 10*60;
  const unsigned Multiplicity=2;
  const double femto_region = 100;
  const unsigned target_yield = YIELD_PER_CONFIG;
  TString BaseName = TString::Format("Eta%.1f",EtaCut);

  //we run to either reproduce the core of 0.97,
  //or the upper limit of reff = 1.06+0.04
  //this leads to a 10% difference in the SP core source
  double rSP_core=0;
  double rSP_dispZ=0;
  double rSP_hadr=0;
  double rSP_hadrZ=0;
  double rSP_hflc=0;
  double rSP_tau=0;
  bool tau_prp = true;
  double rSP_tflc=0;
  double rSP_ThK=0;
  //default is true
  bool rSP_FixedHadr = true;
  //default is 0
  float rSP_FragBeta = 0;

  TRandom3 rangen(SEED);

  //TString OutputFolderName;
  //OutputFolderName = TString::Format("%s/FunWithCeca/CECA_sim_for_AI_pp_v3",GetFemtoOutputFolder());
  //TFile fOutput(BaseFileName+".root","recreate");

  // Save TNtuple to ROOT file
  TString pair_root_name(TString::Format(OutputFolderName+"/pair_info_S%i_N%i_YPC%iK.root",SEED,NUM_CPU,K_YIELD_PER_CONFIG).Data());
  TFile* outFile = new TFile(pair_root_name, "RECREATE");

  // Create TNtuple
  TNtuple* ntuple = new TNtuple("CECA_Data", "CECA_Data", "seed:config_id:kstar:rstar:mT:rd:h:tau");


  for(unsigned uConfig=0; uConfig<NUM_CONFIGS; uConfig++){

    rSP_core = rangen.Uniform(min_rd, max_rd);
    rSP_dispZ = rSP_core;

    rSP_hadr = rangen.Uniform(min_ht, max_ht);

    rSP_tau = rangen.Uniform(min_tau, max_tau);

    TREPNI Database(0);
    Database.SetSeed(11);
    std::vector<TreParticle*> ParticleList;
    ParticleList.push_back(Database.NewParticle("Proton"));
    ParticleList.push_back(Database.NewParticle("Pion"));
    ParticleList.push_back(Database.NewParticle("ProtonReso"));    

    for(TreParticle* prt : ParticleList){
      if(prt->GetName()=="Proton"){
        prt->SetMass(Mass_p);
        prt->SetAbundance(35.78+64.22*(!PROTON_RESO));
        prt->SetRadius(HadronSize);
        prt->SetRadiusSlope(HadronSlope);
        prt->SetPtEtaPhi(*dlm_pT_eta_p);
      }
      else if(prt->GetName()=="Pion"){
        prt->SetMass(Mass_pic);
        prt->SetAbundance(0);
        prt->SetRadius(HadronSize);
        prt->SetRadiusSlope(HadronSlope);
      }
      else if(prt->GetName()=="ProtonReso"){
        prt->SetMass(1362);
        prt->SetAbundance(64.22*PROTON_RESO);
        prt->SetWidth(hbarc/1.65);
        prt->SetRadius(HadronSize);
        prt->SetRadiusSlope(HadronSlope);
        prt->SetPtEtaPhi(*dlm_pT_eta_p);

        prt->NewDecay();
        prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Proton"));
        prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
        prt->GetDecay(0)->SetBranching(100);
      }
    }

    std::vector<std::string> ListOfParticles;

    ListOfParticles.push_back("Proton");
    ListOfParticles.push_back("Proton");

    CECA Ivana(Database,ListOfParticles);

    Ivana.SetDisplacementZ(rSP_dispZ);
    Ivana.SetDisplacementT(rSP_core);
    //Ivana.SetDisplacementEbe(rSP_core_SIG/rSP_core);
    Ivana.SetHadronizationZ(rSP_hadrZ);//0
    Ivana.SetHadronizationT(rSP_hadr);
    //Ivana.SetHadronizationEbe(rSP_hadr_SIG/rSP_hadr);
    Ivana.SetHadrFluctuation(rSP_hflc);
    Ivana.SetTau(rSP_tau,tau_prp);
    //Ivana.SetTauEbe(rSP_tau_SIG/rSP_tau);
    Ivana.SetTauFluct(rSP_tflc);
    Ivana.SetThermalKick(rSP_ThK);
    Ivana.SetFixedHadr(rSP_FixedHadr);
    Ivana.SetFragmentBeta(rSP_FragBeta);

    Ivana.SetTargetStatistics(target_yield);
    Ivana.SetEventMult(Multiplicity);
    Ivana.SetSourceDim(2);
    Ivana.SetDebugMode(true);
    Ivana.SetThreadTimeout(TIMEOUT);
    Ivana.SetGlobalTimeout(GLOB_TIMEOUT);
    Ivana.EqualizeFsiTime(EQUALIZE_TAU);
    Ivana.SetFemtoRegion(femto_region);
    Ivana.GHETTO_EVENT = true;
    std::string pair_file_name(TString::Format(OutputFolderName+"/pair_info_S%i_N%i_YPC%iK_C%u.data",SEED,NUM_CPU,K_YIELD_PER_CONFIG,uConfig).Data());
    Ivana.SetExportPairs(flag%10,pair_file_name);

    if(NUM_CPU>1){
      Ivana.SetDebugMode(true);
      for(unsigned uTh=0; uTh<NUM_CPU; uTh++){
        Ivana.SetSeed(uTh,SEED*(NUM_CPU)+uTh);
      }
    }
    else{
      Ivana.SetDebugMode(true);
      Ivana.SetSeed(0,SEED);
    }

    Ivana.GoBabyGo(NUM_CPU);

    // Open the text file
    std::ifstream inFile(pair_file_name);

    if (!inFile.is_open()) {
        std::cerr << "Error opening input file." << std::endl;
        return;
    }

    // Skip the first line
    std::string line;
    std::getline(inFile, line);

    outFile->cd();
    // Read and process the remaining lines
    while (std::getline(inFile, line)) {
        std::istringstream iss(line);
        float kstar, rstar, mT, rd, h, tau;

        // Read values from the line
        iss >> kstar >> rstar >> mT >> rd >> h >> tau;

        // Fill the TNtuple with the read values
        ntuple->Fill(SEED, uConfig, kstar, rstar, mT, rd, h, tau);
    }

    // Close the input file
    inFile.close();

    if (std::remove(pair_file_name.c_str()) == 0) {
      std::cout << "File deleted successfully.\n";
    } 
    else{
      std::perror("Error deleting file");
    }
  
  

  }//uConfig

  outFile->cd();
  ntuple->Write();

  // Clean up memory
  delete ntuple;
  delete outFile;
  
//rd,h,tau --> get a list of pairs with r* and mT
//sample rd,h,tau randomly from a gauss
//for the output, we need not a histo, but an ASCII or binary file with:
//kstar,rstar,mT,rd,h,tau

  delete dlm_pT_eta_p;
}




// rho0 = 2.5, RM = 1.25, based on the two-body source
//it seems raffa used for the plot M = 3mp, just as it comes out from the code below, as it stands on 27th March 2024
// the best fit of the 3B points to rho0 = 1.8 to 2.0 (RM = 0.9 - 1)
void ThreeBody_test1(){

  double HadronSize = 0;//0.75
  double HadronSlope = 0;//0.2
  const double EtaCut = 0.8;
  const bool PROTON_RESO = true;
  const bool EQUALIZE_TAU = true;
  //const double TIMEOUT = 60;
  //const double GLOB_TIMEOUT = 30*60;
  const double TIMEOUT = 20;
  const double GLOB_TIMEOUT = 1*60;
  const unsigned Multiplicity=3;
  const bool ThreeBody = true;
  const double femto_region = ThreeBody?400:100;
  const unsigned target_yield = 4096*1000;
  const unsigned NUM_CPU = 8;
  //2 or 3- body
  

  DLM_Histo<float>* dlm_pT_eta_p;
  DLM_Histo<float>* dlm_pT_eta_L;

  dlm_pT_eta_p = GetPtEta_13TeV(
    TString::Format("%s/CatsFiles/Source/CECA/proton_pT/p_dist_13TeV_ClassI.root",GetCernBoxDimi()),
    "Graph1D_y1", 500, 4050, EtaCut);
  dlm_pT_eta_L = GetPtEta_13TeV(
    TString::Format("%s/CatsFiles/Source/CECA/Lambda_pT/L_dist_13TeV_ClassI.root",GetCernBoxDimi()),
    "Graph1D_y1", 400, 8000, EtaCut);

  double rSP_core=0;
  double rSP_dispZ=0;
  double rSP_hadr=0;
  double rSP_hadrZ=0;
  double rSP_hflc=0;
  double rSP_tau=0;
  bool tau_prp = true;
  double rSP_tflc=0;
  double rSP_ThK=0;
  //default is true
  bool rSP_FixedHadr = true;
  //default is 0
  float rSP_FragBeta = 0;

  TString type = "ppp";

  if(type=="ppp"){
    //this so happens to be exactly as the usmani fit in the ceca paper
    rSP_core = 0.176;
    rSP_dispZ = 0.176;
    rSP_hadr = 2.68;
    rSP_tau = 3.76;

    //a silly test with a gauss only
    //rSP_core = 1.10;
    //rSP_dispZ = 1.10;
    //rSP_hadr = 0.0;
    //rSP_tau = 0.0;   
  }
  else{
    printf("BAD STUFF!\n");
    return;
  }

  TREPNI Database(0);
  Database.SetSeed(11);
  std::vector<TreParticle*> ParticleList;
  ParticleList.push_back(Database.NewParticle("Proton"));
  ParticleList.push_back(Database.NewParticle("Pion"));
  ParticleList.push_back(Database.NewParticle("ProtonReso"));

  TString OutputFolderName;
  OutputFolderName = "FunWithCeca/ThreeBody_test1";
  TString BaseFileName = TString::Format("%s/%s/fOutput_%uM_%uB",GetFemtoOutputFolder(),OutputFolderName.Data(), Multiplicity, unsigned(ThreeBody)+2);

  TFile fOutput(BaseFileName+".root","recreate");

  for(TreParticle* prt : ParticleList){
    if(prt->GetName()=="Proton"){
      prt->SetMass(Mass_p);
      if(type=="pp"||type=="ppp") prt->SetAbundance(35.78+64.22*(!PROTON_RESO));
      //if(type=="pp"||type=="ppp") prt->SetAbundance(0);
      else if(type=="pP") prt->SetAbundance(100.);
      else prt->SetAbundance(0);
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);
      if(dlm_pT_eta_p) prt->SetPtEtaPhi(*dlm_pT_eta_p);
      else
        prt->SetPtPz(0.85*prt->GetMass(),0.85*prt->GetMass());
    }
    else if(prt->GetName()=="Pion"){
        prt->SetMass(Mass_pic);
        prt->SetAbundance(0);
        prt->SetRadius(HadronSize);
        prt->SetRadiusSlope(HadronSlope);
    }
    else if(prt->GetName()=="ProtonReso"){
      prt->SetMass(1362);
      prt->SetAbundance(64.22*PROTON_RESO);
      //prt->SetAbundance(100);
      prt->SetWidth(hbarc/1.65);
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);
      if(dlm_pT_eta_p) prt->SetPtEtaPhi(*dlm_pT_eta_p);
      else
        prt->SetPtPz(0.85*prt->GetMass(),0.85*prt->GetMass());

      prt->NewDecay();
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Proton"));
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
      prt->GetDecay(0)->SetBranching(100);
    }
  }

  std::vector<std::string> ListOfParticles;

  ListOfParticles.push_back("Proton");
  ListOfParticles.push_back("Proton");
  if(ThreeBody)
    ListOfParticles.push_back("Proton");

  CECA Ivana(Database,ListOfParticles);

  Ivana.SetDisplacementZ(rSP_dispZ);
  Ivana.SetDisplacementT(rSP_core);
  Ivana.SetHadronizationZ(rSP_hadrZ);//0
  Ivana.SetHadronizationT(rSP_hadr);
  Ivana.SetHadrFluctuation(rSP_hflc);
  Ivana.SetTau(rSP_tau,tau_prp);
  Ivana.SetTauFluct(rSP_tflc);
  Ivana.SetThermalKick(rSP_ThK);
  Ivana.SetFixedHadr(rSP_FixedHadr);
  Ivana.SetFragmentBeta(rSP_FragBeta);

  Ivana.SetTargetStatistics(target_yield);
  Ivana.SetEventMult(Multiplicity);
  Ivana.SetSourceDim(2+ThreeBody);
  Ivana.SetDebugMode(true);
  Ivana.SetThreadTimeout(TIMEOUT);
  Ivana.SetGlobalTimeout(GLOB_TIMEOUT);
  Ivana.EqualizeFsiTime(EQUALIZE_TAU);
  Ivana.SetFemtoRegion(femto_region);
  Ivana.GHETTO_EVENT = true;

  Ivana.Ghetto_NumRadBins *= 2;

  printf("GoBabyGo\n");
  Ivana.GoBabyGo(NUM_CPU);

  fOutput.cd();
  Ivana.Ghetto_kstar_rstar->ComputeError();
  TH2F* h_Ghetto_kstar_rstar = Convert_DlmHisto_TH2F(Ivana.Ghetto_kstar_rstar,"Ghetto_kstar_rstar");

  Ivana.GhettoFemto_mT_rstar->ComputeError();
  TH2F* h_GhettoFemto_mT_rstar = Convert_DlmHisto_TH2F(Ivana.GhettoFemto_mT_rstar,"GhettoFemto_mT_rstar");

  fOutput.cd();
  h_Ghetto_kstar_rstar->Write();
  h_GhettoFemto_mT_rstar->Write();
}

void ThreeBody_test1_plots(){

  TString OutputFolderName;
  OutputFolderName = "FunWithCeca/ThreeBody_test1";
  TString InFile_2M_2B = TString::Format("%s/%s/2025-03-14_CECA/fOutput_2M_2B.root",GetFemtoOutputFolder(),OutputFolderName.Data());
  TString InFile_3M_2B = TString::Format("%s/%s/2025-03-14_CECA/fOutput_3M_2B.root",GetFemtoOutputFolder(),OutputFolderName.Data());
  TString InFile_3M_3B = TString::Format("%s/%s/2025-03-14_CECA/fOutput_3M_3B.root",GetFemtoOutputFolder(),OutputFolderName.Data());
  TString OutFileName = TString::Format("%s/%s/2025-03-14_CECA/fPlots.root",GetFemtoOutputFolder(),OutputFolderName.Data());


  const double Q3_femto = 400;
  const double kstar_femto = 100;

  TFile fIn3M_3B(InFile_3M_3B,"read");
  TH2F* h3M_3B_kstar_rstar = (TH2F*)fIn3M_3B.Get("Ghetto_kstar_rstar");
  TH2F* h3M_3B_mT_rstar = (TH2F*)fIn3M_3B.Get("GhettoFemto_mT_rstar");

  TFile fIn3M_2B(InFile_3M_2B,"read");
  TH2F* h3M_2B_kstar_rstar = (TH2F*)fIn3M_2B.Get("Ghetto_kstar_rstar");
  TH2F* h3M_2B_mT_rstar = (TH2F*)fIn3M_2B.Get("GhettoFemto_mT_rstar");
  
  TFile fIn2M_2B(InFile_2M_2B,"read");
  TH2F* h2M_2B_kstar_rstar = (TH2F*)fIn2M_2B.Get("Ghetto_kstar_rstar");
  TH2F* h2M_2B_mT_rstar = (TH2F*)fIn2M_2B.Get("GhettoFemto_mT_rstar");

  unsigned kstar_bin_min;
  unsigned kstar_bin_max;
  unsigned NumMt;
  TH1F* hSrc_2M_2B;
  TH1F* hSrc_3M_2B;
  TH1F* hSrc_3M_3B;

  TFile fOutput(OutFileName, "recreate");

  //3B source, mT integrated
  kstar_bin_min = h3M_3B_kstar_rstar->GetXaxis()->FindBin(0.);
  kstar_bin_max = h3M_3B_kstar_rstar->GetXaxis()->FindBin(Q3_femto);
  fOutput.cd();
  hSrc_3M_3B = (TH1F*)h3M_3B_kstar_rstar->ProjectionY(TString::Format("hSrc_3M_3B"),kstar_bin_min,kstar_bin_max);

  hSrc_3M_3B->SetStats(false);
  hSrc_3M_3B->SetTitle("");
  hSrc_3M_3B->GetXaxis()->SetTitle("#rho (fm)");
  hSrc_3M_3B->GetXaxis()->SetTitleSize(0.06);
  hSrc_3M_3B->GetXaxis()->SetLabelSize(0.06);
  hSrc_3M_3B->GetXaxis()->SetTitleOffset(1.3);
  hSrc_3M_3B->GetXaxis()->SetLabelOffset(0.02);
  hSrc_3M_3B->GetXaxis()->SetRangeUser(0.0,5.0);
  hSrc_3M_3B->GetYaxis()->SetTitle("#rho^{5}S(#rho) (1/fm)");
  hSrc_3M_3B->GetYaxis()->SetTitleSize(0.06);
  hSrc_3M_3B->GetYaxis()->SetLabelSize(0.06);
  hSrc_3M_3B->GetYaxis()->SetTitleOffset(1.00);
  hSrc_3M_3B->GetYaxis()->SetRangeUser(0.001,1.50);
  hSrc_3M_3B->SetLineColor(kOrange+2);
  hSrc_3M_3B->SetLineWidth(6);
  hSrc_3M_3B->Scale(1./hSrc_3M_3B->Integral(),"width");
  hSrc_3M_3B->Write();
  printf("3M_3B mean %.3f\n",hSrc_3M_3B->GetMean());

  TH1F* hSrc_3M_3B_red = new TH1F("hSrc_3M_3B_red","hSrc_3M_3B_red",hSrc_3M_3B->GetNbinsX(),hSrc_3M_3B->GetXaxis()->GetXmin()*0.5,hSrc_3M_3B->GetXaxis()->GetXmax()*0.5);
  for(unsigned uBin=0; uBin<hSrc_3M_3B_red->GetNbinsX(); uBin++){
    hSrc_3M_3B_red->SetBinContent(uBin+1,hSrc_3M_3B->GetBinContent(uBin+1)*2);
  }
  hSrc_3M_3B_red->Write();

  TH1F* hSrc_Ana_3B = (TH1F*)hSrc_3M_3B->Clone("hSrc_Ana_3B");
  const double RM = 1.25;
  const double rho0 = sqrt(2./3.)*RM;
  for(unsigned uRho=0; uRho<hSrc_Ana_3B->GetNbinsX(); uRho++){
    double rho = hSrc_Ana_3B->GetBinCenter(uRho+1);
    
    double src3b = pow(rho,5)*exp(-rho*rho/rho0/rho0)/(pow(TMath::Pi(),3)*pow(rho0,6.));
    hSrc_Ana_3B->SetBinContent(uRho+1, src3b);
    hSrc_Ana_3B->SetBinError(uRho+1, 0);
  }
  hSrc_Ana_3B->SetLineColor(kGreen+2);
  hSrc_Ana_3B->SetLineWidth(4);
  hSrc_Ana_3B->Scale(1./hSrc_Ana_3B->Integral(),"width");
  hSrc_Ana_3B->Write();
  printf("Ana_3B mean %.3f\n",hSrc_Ana_3B->GetMean());


  //3B source, mT scaling
  NumMt = h3M_3B_mT_rstar->GetXaxis()->GetNbins();
  fOutput.cd();
  TH1F* hmT_3M_3B = (TH1F*)h3M_3B_mT_rstar->ProjectionX(TString::Format("hmT_3M_3B"),1,h3M_3B_mT_rstar->GetXaxis()->GetNbins());

  hmT_3M_3B->SetStats(false);
  hmT_3M_3B->SetTitle("");
  hmT_3M_3B->GetXaxis()->SetTitle("m_{T} (MeV)");
  hmT_3M_3B->GetXaxis()->SetTitleSize(0.06);
  hmT_3M_3B->GetXaxis()->SetLabelSize(0.06);
  hmT_3M_3B->GetXaxis()->SetTitleOffset(1.3);
  hmT_3M_3B->GetXaxis()->SetLabelOffset(0.02);
  hmT_3M_3B->GetXaxis()->SetRangeUser(1000,2400);
  hmT_3M_3B->GetYaxis()->SetTitle("<#rho> or <r*> (fm)");
  hmT_3M_3B->GetYaxis()->SetTitleSize(0.06);
  hmT_3M_3B->GetYaxis()->SetLabelSize(0.06);
  hmT_3M_3B->GetYaxis()->SetTitleOffset(1.00);
  //hmT_3M_3B->GetYaxis()->SetRangeUser(1000,2400);
  hmT_3M_3B->SetLineColor(kOrange+2);
  hmT_3M_3B->SetLineWidth(8);

  for(unsigned uMt=0; uMt<NumMt; uMt++){
    //printf("uMt %u\n",uMt);
    gROOT->cd();
    TH1F* hProj = (TH1F*)h3M_3B_mT_rstar->ProjectionY(TString::Format("hProj"),uMt+1,uMt+1);
    double rad_mean = hProj->GetMean();
    //printf(" %f\n",rad_mean);
    hmT_3M_3B->SetBinContent(uMt+1, rad_mean);
    hmT_3M_3B->SetBinError(uMt+1, 0);
    delete hProj;
  }
  fOutput.cd();
  hmT_3M_3B->Write();


  //3B sim, 2B source, mT integrated
  kstar_bin_min = h3M_2B_kstar_rstar->GetXaxis()->FindBin(0.);
  kstar_bin_max = h3M_2B_kstar_rstar->GetXaxis()->FindBin(kstar_femto);
  fOutput.cd();
  hSrc_3M_2B = (TH1F*)h3M_2B_kstar_rstar->ProjectionY(TString::Format("hSrc_3M_2B"),kstar_bin_min,kstar_bin_max);
  hSrc_3M_2B->SetStats(false);
  hSrc_3M_2B->SetTitle("");
  hSrc_3M_2B->GetXaxis()->SetTitle("r (fm)");
  hSrc_3M_2B->GetXaxis()->SetTitleSize(0.06);
  hSrc_3M_2B->GetXaxis()->SetLabelSize(0.06);
  hSrc_3M_2B->GetXaxis()->SetTitleOffset(1.3);
  hSrc_3M_2B->GetXaxis()->SetLabelOffset(0.02);
  hSrc_3M_2B->GetXaxis()->SetRangeUser(0.0,7.0);
  hSrc_3M_2B->GetYaxis()->SetTitle("4#pir*^{2}S(r) (1/fm)");
  hSrc_3M_2B->GetYaxis()->SetTitleSize(0.06);
  hSrc_3M_2B->GetYaxis()->SetLabelSize(0.06);
  hSrc_3M_2B->GetYaxis()->SetTitleOffset(1.00);
  hSrc_3M_2B->GetYaxis()->SetRangeUser(0.001,1.50);
  hSrc_3M_2B->SetLineColor(kBlue);
  hSrc_3M_2B->SetLineWidth(4);
  hSrc_3M_2B->Scale(1./hSrc_3M_2B->Integral(),"width");
  hSrc_3M_2B->Write();
  printf("3M_2B mean %.3f\n",hSrc_3M_2B->GetMean());


  //3B sim, 2B source, mT scaling
  NumMt = h3M_2B_mT_rstar->GetXaxis()->GetNbins();
  fOutput.cd();
  TH1F* hmT_3M_2B = (TH1F*)h3M_2B_mT_rstar->ProjectionX(TString::Format("hmT_3M_2B"),1,h3M_2B_mT_rstar->GetXaxis()->GetNbins());

  hmT_3M_2B->SetStats(false);
  hmT_3M_2B->SetTitle("");
  hmT_3M_2B->GetXaxis()->SetTitle("m_{T} (MeV)");
  hmT_3M_2B->GetXaxis()->SetTitleSize(0.06);
  hmT_3M_2B->GetXaxis()->SetLabelSize(0.06);
  hmT_3M_2B->GetXaxis()->SetTitleOffset(1.3);
  hmT_3M_2B->GetXaxis()->SetLabelOffset(0.02);
  hmT_3M_2B->GetXaxis()->SetRangeUser(1000,2400);
  hmT_3M_2B->GetYaxis()->SetTitle("<#rho> or <r*> (fm)");
  hmT_3M_2B->GetYaxis()->SetTitleSize(0.06);
  hmT_3M_2B->GetYaxis()->SetLabelSize(0.06);
  hmT_3M_2B->GetYaxis()->SetTitleOffset(1.00);
  //hmT_3M_2B->GetYaxis()->SetRangeUser(1000,2400);
  hmT_3M_2B->SetLineColor(kBlue);
  hmT_3M_2B->SetLineWidth(8);

  for(unsigned uMt=0; uMt<NumMt; uMt++){
    //printf("uMt %u\n",uMt);
    gROOT->cd();
    TH1F* hProj = (TH1F*)h3M_2B_mT_rstar->ProjectionY(TString::Format("hProj"),uMt+1,uMt+1);
    double rad_mean = hProj->GetMean();
    //printf(" %f\n",rad_mean);
    hmT_3M_2B->SetBinContent(uMt+1, rad_mean);
    hmT_3M_2B->SetBinError(uMt+1, 0);
    delete hProj;
  }
  fOutput.cd();
  hmT_3M_2B->Write();


  //2B source, mT integrated
  kstar_bin_min = h2M_2B_kstar_rstar->GetXaxis()->FindBin(0.);
  kstar_bin_max = h2M_2B_kstar_rstar->GetXaxis()->FindBin(kstar_femto);
  fOutput.cd();
  hSrc_2M_2B = (TH1F*)h2M_2B_kstar_rstar->ProjectionY(TString::Format("hSrc_2M_2B"),kstar_bin_min,kstar_bin_max);
  
  
  TH1F* hAxisSource = new TH1F("hAxisSource","hAxisSource",128,0,8);
  hSrc_2M_2B->SetStats(false);
  hSrc_2M_2B->SetTitle("");
  hSrc_2M_2B->GetXaxis()->SetTitle("r (fm)");
  hSrc_2M_2B->GetXaxis()->SetTitleSize(0.06);
  hSrc_2M_2B->GetXaxis()->SetLabelSize(0.06);
  hSrc_2M_2B->GetXaxis()->SetTitleOffset(1.3);
  hSrc_2M_2B->GetXaxis()->SetLabelOffset(0.02);
  hSrc_2M_2B->GetXaxis()->SetRangeUser(0.0,7.0);
  hSrc_2M_2B->GetYaxis()->SetTitle("4#pir*^{2}S(r) (1/fm)");
  hSrc_2M_2B->GetYaxis()->SetTitleSize(0.06);
  hSrc_2M_2B->GetYaxis()->SetLabelSize(0.06);
  hSrc_2M_2B->GetYaxis()->SetTitleOffset(1.00);
  hSrc_2M_2B->GetYaxis()->SetRangeUser(0.001,1.50);
  hSrc_2M_2B->SetLineColor(kViolet);
  hSrc_2M_2B->SetLineWidth(6);
  hSrc_2M_2B->Scale(1./hSrc_2M_2B->Integral(),"width");
  hSrc_2M_2B->Write();
  printf("2M_2B mean %.3f\n",hSrc_2M_2B->GetMean());

  TH1F* hSrc_Ana_2B = (TH1F*)hSrc_2M_2B->Clone("hSrc_Ana_2B");
  for(unsigned uRad=0; uRad<hSrc_Ana_2B->GetNbinsX(); uRad++){
    double rad = hSrc_Ana_2B->GetBinCenter(uRad+1);
    double rad0 = RM;
    double src2b = 4.*TMath::Pi()*rad*rad/(pow(4.*TMath::Pi()*rad0*rad0,1.5))*exp(-rad*rad/(4.*rad0*rad0));
    hSrc_Ana_2B->SetBinContent(uRad+1, src2b);
    hSrc_Ana_2B->SetBinError(uRad+1, 0);
  }
  hSrc_Ana_2B->SetLineColor(kGreen+2);
  hSrc_Ana_2B->SetLineWidth(4);
  hSrc_Ana_2B->Scale(1./hSrc_Ana_2B->Integral(),"width");
  hSrc_Ana_2B->Write();
  printf("Ana_2B mean %.3f\n",hSrc_Ana_2B->GetMean());



  //2B source, mT scaling
  NumMt = h2M_2B_mT_rstar->GetXaxis()->GetNbins();
  fOutput.cd();
  TH1F* hmT_2M_2B = (TH1F*)h2M_2B_mT_rstar->ProjectionX(TString::Format("hmT_2M_2B"),1,h2M_2B_mT_rstar->GetXaxis()->GetNbins());

  hmT_2M_2B->SetStats(false);
  hmT_2M_2B->SetTitle("");
  hmT_2M_2B->GetXaxis()->SetTitle("m_{T} (MeV)");
  hmT_2M_2B->GetXaxis()->SetTitleSize(0.06);
  hmT_2M_2B->GetXaxis()->SetLabelSize(0.06);
  hmT_2M_2B->GetXaxis()->SetTitleOffset(1.3);
  hmT_2M_2B->GetXaxis()->SetLabelOffset(0.02);
  hmT_2M_2B->GetXaxis()->SetRangeUser(1000,2400);
  hmT_2M_2B->GetYaxis()->SetTitle("<#rho> or <r*> (fm)");
  hmT_2M_2B->GetYaxis()->SetTitleSize(0.06);
  hmT_2M_2B->GetYaxis()->SetLabelSize(0.06);
  hmT_2M_2B->GetYaxis()->SetTitleOffset(1.00);
  //hmT_2M_2B->GetYaxis()->SetRangeUser(1000,2400);
  hmT_2M_2B->SetLineColor(kViolet);
  hmT_2M_2B->SetLineWidth(8);

  for(unsigned uMt=0; uMt<NumMt; uMt++){
    //printf("uMt %u\n",uMt);
    gROOT->cd();
    TH1F* hProj = (TH1F*)h2M_2B_mT_rstar->ProjectionY(TString::Format("hProj"),uMt+1,uMt+1);
    double rad_mean = hProj->GetMean();
    //printf(" %f\n",rad_mean);
    hmT_2M_2B->SetBinContent(uMt+1, rad_mean);
    hmT_2M_2B->SetBinError(uMt+1, 0);
    delete hProj;
  }
  fOutput.cd();
  hmT_2M_2B->Write();

  //delete hSrc_3M_3B;
  //delete hmT_3M_3B;
  delete hSrc_Ana_3B;
}


//generate a source, and a create a dummy correlation with errors. Fit it using CECA source.
//give it to Max, see what he gets with ML
void MC_closure_for_Max(int SEED){
  //std::vector<double> input_rd = {0.21, 0.23, 0.25, 0.27, 0.29};
  //std::vector<double> input_ht = {3.9, 4.2, 4.5, 4.8, 5.1};
  //std::vector<double> input_tau = {2.1, 2.3, 2.5, 2.7, 2.9};

  //std::vector<double> input_rd = {0.22, 0.25, 0.28};
  //std::vector<double> input_ht = {4.0, 4.5, 5.0};
  //std::vector<double> input_tau = {2.2, 2.5, 2.8};

  std::vector<double> input_rd = {0.22, 0.28};
  std::vector<double> input_ht = {4.0, 5.0};
  std::vector<double> input_tau = {2.2, 2.8};


  std::vector<float> mT_values = {1100, 1200, 1300, 1400, 1600, 1900};

  double kMin = 0;
  double kMax = 160;
  unsigned kSteps = 20;
  //in bins of 8 MeV, errors given in %
  TH1F* hErrors = new TH1F("hErrors", "hErrors", kSteps, kMin, kMax);
  hErrors->SetBinContent(1, 1);
  hErrors->SetBinError(1, 20./100.);
  for(unsigned uMom=1; uMom<kSteps; uMom++){
    hErrors->SetBinContent(uMom+1, 1);
    double old_error = hErrors->GetBinError(uMom);
    double new_error = old_error / pow(2,0.5);
    if(new_error<0.08) new_error = old_error / pow(2,0.25);
    if(new_error<0.02) new_error = old_error / pow(2,0.125);
    hErrors->SetBinError(uMom+1, new_error);
    printf("%u %.3f\n", uMom, new_error*100);
  }

  double EtaCut = 0.8;
  DLM_Histo<float>* dlm_pT_eta_p;
  dlm_pT_eta_p = GetPtEta_13TeV(
      TString::Format("%s/CatsFiles/Source/CECA/proton_pT/p_dist_13TeV_ClassI.root",GetCernBoxDimi()),
      "Graph1D_y1", 500, 4050, EtaCut);

  double HadronSize = 0;//0.75
  double HadronSlope = 0;//0.2
  const bool PROTON_RESO = true;
  const bool EQUALIZE_TAU = true;
  const double TIMEOUT = 15;
  const double GLOB_TIMEOUT = 24*60*60;
  const unsigned Multiplicity=2;
  const double femto_region = 100;
  const unsigned target_yield = 11*23*1000;
  //const unsigned target_yield = 3*23*1000;
  const unsigned NUM_CPU = 8;
  TString BaseName = TString::Format("Eta%.1f",EtaCut);

  //we run to either reproduce the core of 0.97,
  //or the upper limit of reff = 1.06+0.04
  //this leads to a 10% difference in the SP core source
  double rSP_core=0;
  double rSP_dispZ=0;
  double rSP_hadr=0;
  double rSP_hadrZ=0;
  double rSP_hflc=0;
  double rSP_tau=0;
  bool tau_prp = true;
  double rSP_tflc=0;
  double rSP_ThK=0;
  //default is true
  bool rSP_FixedHadr = true;
  //default is 0
  float rSP_FragBeta = 0;

  TRandom3 rangen(SEED);

  TString OutputFolderName;
  OutputFolderName = TString::Format("%s/FunWithCeca/MC_closure_for_Max",GetFemtoOutputFolder());
  //TFile fOutput(BaseFileName+".root","recreate");

  // Save TNtuple to ROOT file
  TString pair_root_name(TString::Format(OutputFolderName+"/pp_S%i_TY%iK.root",SEED,target_yield/1000).Data());
  TFile* outFile = new TFile(pair_root_name, "RECREATE");

  // Create TNtuple
  //TNtuple* ntuple = new TNtuple("CECA_Data", "CECA_Data", "config_id:kstar:rstar:mT:rd:h:tau");

  for(unsigned urd=0; urd<input_rd.size(); urd++){
  for(unsigned uht=0; uht<input_ht.size(); uht++){
  for(unsigned utau=0; utau<input_tau.size(); utau++){

  //for(unsigned uConfig=0; uConfig<NUM_CONFIGS; uConfig++){
    //this so happens to be exactly the mean of the results from the ceca paper
    //the errors are the half of the difference of the two ceca results + in ^2 the two errors from the paper
    //rSP_core = (0.176+0.288)*0.5;
    //double rSP_core_AVG = rSP_core;
    //double rSP_core_SIG = fabs(0.176-0.288)*0.5;
    //rSP_core_SIG = sqrt(pow(rSP_core_SIG,2.)+pow(0.013,2.)+pow(0.005,2.));
    rSP_core = input_rd.at(urd);

    rSP_dispZ = rSP_core;
    //double rSP_dispZ_SIG = rSP_core_SIG;

    //rSP_hadr = (2.68+3.23)*0.5;
    //double rSP_hadr_AVG = rSP_hadr;
    //double rSP_hadr_SIG = fabs(2.68-3.23)*0.5;
    //rSP_hadr_SIG = sqrt(pow(rSP_hadr_SIG,2.)+pow(0.175,2.)+pow(0.05,2.));
    rSP_hadr = input_ht.at(uht);

    //rSP_tau = (3.76+3.26)*0.5;
    //double rSP_tau_AVG = rSP_tau;
    //double rSP_tau_SIG = fabs(3.76-3.26)*0.5;
    //rSP_tau_SIG = sqrt(pow(rSP_tau_SIG,2.)+pow(0.1,2.)+pow(0.04,2.));
    rSP_tau = input_tau.at(utau);

    printf("rSP_core = %.4e\n",rSP_core);
    printf("rSP_hadr = %.4e\n",rSP_hadr);
    printf("rSP_tau = %.4e\n",rSP_tau);   

    TREPNI Database(0);
    Database.SetSeed(11);
    std::vector<TreParticle*> ParticleList;
    ParticleList.push_back(Database.NewParticle("Proton"));
    ParticleList.push_back(Database.NewParticle("Pion"));
    ParticleList.push_back(Database.NewParticle("ProtonReso"));    

    for(TreParticle* prt : ParticleList){
      if(prt->GetName()=="Proton"){
        prt->SetMass(Mass_p);
        prt->SetAbundance(35.78+64.22*(!PROTON_RESO));
        prt->SetRadius(HadronSize);
        prt->SetRadiusSlope(HadronSlope);
        prt->SetPtEtaPhi(*dlm_pT_eta_p);
      }
      else if(prt->GetName()=="Pion"){
        prt->SetMass(Mass_pic);
        prt->SetAbundance(0);
        prt->SetRadius(HadronSize);
        prt->SetRadiusSlope(HadronSlope);
      }
      else if(prt->GetName()=="ProtonReso"){
        prt->SetMass(1362);
        prt->SetAbundance(64.22*PROTON_RESO);
        prt->SetWidth(hbarc/1.65);
        prt->SetRadius(HadronSize);
        prt->SetRadiusSlope(HadronSlope);
        prt->SetPtEtaPhi(*dlm_pT_eta_p);

        prt->NewDecay();
        prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Proton"));
        prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
        prt->GetDecay(0)->SetBranching(100);
      }
    }

    std::vector<std::string> ListOfParticles;

    ListOfParticles.push_back("Proton");
    ListOfParticles.push_back("Proton");

    CECA Ivana(Database,ListOfParticles);

    Ivana.SetDisplacementZ(rSP_dispZ);
    Ivana.SetDisplacementT(rSP_core);
    //Ivana.SetDisplacementEbe(rSP_core_SIG/rSP_core);
    Ivana.SetHadronizationZ(rSP_hadrZ);//0
    Ivana.SetHadronizationT(rSP_hadr);
    //Ivana.SetHadronizationEbe(rSP_hadr_SIG/rSP_hadr);
    Ivana.SetHadrFluctuation(rSP_hflc);
    Ivana.SetTau(rSP_tau,tau_prp);
    //Ivana.SetTauEbe(rSP_tau_SIG/rSP_tau);
    Ivana.SetTauFluct(rSP_tflc);
    Ivana.SetThermalKick(rSP_ThK);
    Ivana.SetFixedHadr(rSP_FixedHadr);
    Ivana.SetFragmentBeta(rSP_FragBeta);

    Ivana.SetTargetStatistics(target_yield);
    Ivana.SetEventMult(Multiplicity);
    Ivana.SetSourceDim(2);
    Ivana.SetDebugMode(true);
    Ivana.SetThreadTimeout(TIMEOUT);
    Ivana.SetGlobalTimeout(GLOB_TIMEOUT);
    Ivana.EqualizeFsiTime(EQUALIZE_TAU);
    Ivana.SetFemtoRegion(femto_region);
    Ivana.GHETTO_EVENT = true;
    std::string pair_file_name(TString::Format(OutputFolderName+"/pair_info_S%i_TY%iK.data",SEED,target_yield/1000).Data());
    Ivana.SetExportPairs(2,pair_file_name);

    if(NUM_CPU>1){
      Ivana.SetDebugMode(true);
      for(unsigned uTh=0; uTh<NUM_CPU; uTh++){
        Ivana.SetSeed(uTh,SEED*(NUM_CPU)+uTh);
      }
    }
    else{
      Ivana.SetDebugMode(false);
      Ivana.SetSeed(0,SEED);
    }

    Ivana.Ghetto_NumRadBins = 4096;
    Ivana.Ghetto_RadMax = 64;
    Ivana.Ghetto_NumMtBins = 192;
    Ivana.GoBabyGo(NUM_CPU);

    Ivana.GhettoFemto_mT_rstar->ComputeError();


    Ivana.GhettoFemto_mT_kstar->ComputeError();
    outFile->cd();
    TH2F* h_GhettoFemto_mT_kstar = Convert_DlmHisto_TH2F(Ivana.GhettoFemto_mT_kstar,TString::Format("GhettoFemto_mT_kstar_%.2f_%.2f_%.2f",rSP_core, rSP_hadr, rSP_tau));
    TH2F* h_GhettoFemto_mT_rstar = Convert_DlmHisto_TH2F(Ivana.GhettoFemto_mT_rstar,TString::Format("GhettoFemto_mT_rstar_%.2f_%.2f_%.2f",rSP_core, rSP_hadr, rSP_tau));
    

    for(unsigned uMt=0; uMt<mT_values.size(); uMt++){

      double mT_val = mT_values.at(uMt);
      
      

      DLM_Histo<float> HistoSourceR;
      HistoSourceR.SetUp(1);
      HistoSourceR.SetUp(0,Ivana.Ghetto_NumRadBins,Ivana.Ghetto_RadMin,Ivana.Ghetto_RadMax);
      HistoSourceR.Initialize();

      TH1F* hSource = new TH1F(TString::Format("hSource_%.0f_%.2f_%.2f_%.2f", mT_val, rSP_core, rSP_hadr, rSP_tau), 
        TString::Format("hSource_%.0f_%.2f_%.2f_%.2f", mT_val, rSP_core, rSP_hadr, rSP_tau), 
        Ivana.Ghetto_NumRadBins,Ivana.Ghetto_RadMin,Ivana.Ghetto_RadMax);
      
      //TH1F* hCk = new TH1F(TString::Format("hCk_%.0f", mT_val), TString::Format("hCk_%.0f", mT_val), 
      //  kSteps, kMin, kMax); 
      for(unsigned uRad=0; uRad<Ivana.Ghetto_NumRadBins; uRad++){
        //mT,rstar
        double rad_val = HistoSourceR.GetBinCenter(0, uRad);
        //HistoSourceR.SetBinContent(uRad, Ivana.GhettoFemto_mT_rstar->Eval2D(mT_val, rad_val));
        //HistoSourceR.SetBinError(uRad, Ivana.GhettoFemto_mT_rstar->EvalError2D(mT_val, rad_val));
        hSource->SetBinContent(uRad+1, Ivana.GhettoFemto_mT_rstar->Eval2D(mT_val, rad_val));
        hSource->SetBinError(uRad+1, Ivana.GhettoFemto_mT_rstar->EvalError2D(mT_val, rad_val));
      }
      hSource->Scale(1./hSource->Integral(), "width");
      for(unsigned uRad=0; uRad<Ivana.Ghetto_NumRadBins; uRad++){
        HistoSourceR.SetBinContent(uRad, hSource->GetBinContent(uRad+1));
        HistoSourceR.SetBinError(uRad, hSource->GetBinError(uRad+1));      
      }

      //CATS Kitty_Dummy;
      //Kitty_Dummy.SetMomBins(kSteps,kMin,kMax);
      //DLM_CommonAnaFunctions AnalysisObject;
      //AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles/",GetCernBoxDimi()));
      //AnalysisObject.SetUpCats_pp(Kitty_Dummy,"AV18","",0,0);
      //DLM_HistoSource src_CECA(HistoSourceR);
      //Kitty_Dummy.SetAnaSource(CatsSourceForwarder, &src_CECA, 0);
      //Kitty_Dummy.SetUseAnalyticSource(true);
      //Kitty_Dummy.SetAutoNormSource(false);
      //Kitty_Dummy.SetNormalizedSource(true);

      //Kitty_Dummy.KillTheCat();

      //for(unsigned uMom=0; uMom<kSteps; uMom++){
      //  double Ck_Mean = Kitty_Dummy.GetCorrFun(uMom);
      //  double Ck_Err = Ck_Mean*hErrors->GetBinError(uMom+1);
      //  if(uMom==0) Ck_Err = hErrors->GetBinError(uMom+1);
      //  double Ck_Val = -1;
      //  while(Ck_Val<0){
      //    Ck_Val = rangen.Gaus(Ck_Mean, Ck_Err);
      //  }
      //  printf("%u %.3f %.3f\n", uMom, Ck_Mean, Ck_Err);
      //  hCk->SetBinContent(uMom+1, Ck_Val);
      //  hCk->SetBinError(uMom+1, Ck_Err);
      //}

      outFile->cd();
      hSource->Write();
      
      //hCk->Write();
      delete hSource;
      
      //delete hCk;
    }

    outFile->cd();
    h_GhettoFemto_mT_kstar->Write();
    h_GhettoFemto_mT_rstar->Write();
    
    //ntuple->Write();

    // Clean up memory
    //delete ntuple;
    //delete outFile;
    delete h_GhettoFemto_mT_kstar;
    delete h_GhettoFemto_mT_rstar;
  
  }}} 

  delete outFile;
//rd,h,tau --> get a list of pairs with r* and mT
//sample rd,h,tau randomly from a gauss
//for the output, we need not a histo, but an ASCII or binary file with:
//kstar,rstar,mT,rd,h,tau

  delete dlm_pT_eta_p;
}

//we read a file of TH1s corresponding to different mT, rd, ht, tau
//we save as an output the sources and corresponding correlation functions
//for the central values of the parameters, we also generate a random Cks for each mt  
//which we will end up using for the fitting procedure
void MC_closure_FIT_for_Max(){

  //int default_values = {2,2,2};
  std::vector<double> input_rd = {0.21,0.23, 0.25, 0.27, 0.29};
  std::vector<double> input_ht = {3.9, 4.2, 4.5, 4.8, 5.1};
  std::vector<double> input_tau = {2.1, 2.3, 2.5, 2.7, 2.9};
  std::vector<float> input_mT = {1100, 1200, 1300, 1400, 1600, 1900};

  //std::vector<double> input_rd = {0.21, 0.25};
  //std::vector<double> input_ht = {3.9, 4.5};
  //std::vector<double> input_tau = {2.1, 2.5};
  //std::vector<float> input_mT = {1300, 1400};

  TString BaseFileName = TString::Format("%s/FunWithCeca/MC_closure_for_Max/Official_v2/pp_S69_TY253K",GetFemtoOutputFolder());
  TString InputFileName = TString::Format(BaseFileName+".root",GetFemtoOutputFolder());
  TString OutputFileName = TString::Format(BaseFileName+"_kdp.root",GetFemtoOutputFolder());
  TString OutputFileNameDlm = TString::Format(BaseFileName+"_kdp.dlm",GetFemtoOutputFolder());
  TFile InputFile(InputFileName,"read");
  TFile OutputFile(OutputFileName,"recreate");

  TRandom3 rangen(23);

  unsigned NumBins_mT = input_mT.size();
  unsigned NumBins_rd = input_rd.size();
  unsigned NumBins_ht = input_ht.size();
  unsigned NumBins_tau = input_tau.size();

  // Create a unique pointer to the array
  double* BinRange_mT = new double[NumBins_mT+1];
  double* BinRange_rd = new double[NumBins_rd+1];
  double* BinRange_ht = new double[NumBins_ht+1];
  double* BinRange_tau = new double[NumBins_tau+1];

  double* BinCenter_mT = new double[NumBins_mT];
  double* BinCenter_rd = new double[NumBins_rd];
  double* BinCenter_ht = new double[NumBins_ht];
  double* BinCenter_tau = new double[NumBins_tau];

  BinRange_mT[0] = input_mT.at(0) - 0.5*(input_mT.at(1)-input_mT.at(0));
  BinCenter_mT[0] = input_mT.at(0);
  for(unsigned umT=0; umT<NumBins_mT-1; umT++){
    BinRange_mT[umT+1] = 0.5*(input_mT.at(umT+1)+input_mT.at(umT));
    BinCenter_mT[umT+1] = input_mT.at(umT+1);
  }
  BinRange_mT[NumBins_mT] = input_mT.at(NumBins_mT-1) + 0.5*(input_mT.at(NumBins_mT-1)-input_mT.at(NumBins_mT-2));
  

  BinRange_rd[0] = input_rd.at(0) - 0.5*(input_rd.at(1)-input_rd.at(0));
  BinCenter_rd[0] = input_rd.at(0);
  for(unsigned urd=0; urd<NumBins_rd-1; urd++){
    BinRange_rd[urd+1] = 0.5*(input_rd.at(urd+1)+input_rd.at(urd));
    BinCenter_rd[urd+1] = input_rd.at(urd+1);
  }
  BinRange_rd[NumBins_rd] = input_rd.at(NumBins_rd-1) + 0.5*(input_rd.at(NumBins_rd-1)-input_rd.at(NumBins_rd-2));

  BinRange_ht[0] = input_ht.at(0) - 0.5*(input_ht.at(1)-input_ht.at(0));
  BinCenter_ht[0] = input_ht.at(0);
  for(unsigned uht=0; uht<NumBins_ht-1; uht++){
    BinRange_ht[uht+1] = 0.5*(input_ht.at(uht+1)+input_ht.at(uht));
    BinCenter_ht[uht+1] = input_ht.at(uht+1);
  }
  BinRange_ht[NumBins_ht] = input_ht.at(NumBins_ht-1) + 0.5*(input_ht.at(NumBins_ht-1)-input_ht.at(NumBins_ht-2));

  BinRange_tau[0] = input_tau.at(0) - 0.5*(input_tau.at(1)-input_tau.at(0));
  BinCenter_tau[0] = input_tau.at(0);
  for(unsigned utau=0; utau<NumBins_tau-1; utau++){
    BinRange_tau[utau+1] = 0.5*(input_tau.at(utau+1)+input_tau.at(utau));
    BinCenter_tau[utau+1] = input_tau.at(utau+1);
  }
  BinRange_tau[NumBins_tau] = input_tau.at(NumBins_tau-1) + 0.5*(input_tau.at(NumBins_tau-1)-input_tau.at(NumBins_tau-2));
  for(unsigned utau=0; utau<NumBins_tau; utau++){
    //printf("%.2f %.2f %.2f\n", BinRange_tau[utau], input_tau.at(utau), BinRange_tau[utau+1]);
  }

  //mT, rd, ht, tau
  DLM_Histo<KdpPars> dlm_kdp_pars;
  dlm_kdp_pars.SetUp(4);
  dlm_kdp_pars.SetUp(0, NumBins_mT, BinRange_mT, BinCenter_mT);
  dlm_kdp_pars.SetUp(1, NumBins_rd, BinRange_rd, BinCenter_rd);
  dlm_kdp_pars.SetUp(2, NumBins_ht, BinRange_ht, BinCenter_ht);
  dlm_kdp_pars.SetUp(3, NumBins_tau, BinRange_tau, BinCenter_tau);
  dlm_kdp_pars.Initialize();

  double kMin = 0;
  double kMax = 160;
  unsigned kSteps = 20;

  TH1F* hErrors = new TH1F("hErrors", "hErrors", kSteps, kMin, kMax);
  hErrors->SetBinContent(1, 1);
  hErrors->SetBinError(1, 20./100.);
  for(unsigned uMom=1; uMom<kSteps; uMom++){
    hErrors->SetBinContent(uMom+1, 1);
    double old_error = hErrors->GetBinError(uMom);
    double new_error = old_error / pow(2,0.5);
    if(new_error<0.08) new_error = old_error / pow(2,0.25);
    if(new_error<0.02) new_error = old_error / pow(2,0.125);
    hErrors->SetBinError(uMom+1, new_error);
    //printf("%u %.3f\n", uMom, new_error*100);
  }

  TH1F** hCk_Synthetic = new TH1F*[NumBins_mT];

  for(unsigned umT=0; umT<NumBins_mT; umT++){
    double mT_val = input_mT.at(umT);
    OutputFile.cd();
    hCk_Synthetic[umT] = new TH1F(TString::Format("hCk_Synthetic_mT%.0f",mT_val), 
    TString::Format("hCk_Synthetic_mT%.0f",mT_val), kSteps, kMin, kMax);
    for(unsigned urd=0; urd<NumBins_rd; urd++){
      double rd_val = input_rd.at(urd);
      for(unsigned uht=0; uht<NumBins_ht; uht++){
        double ht_val = input_ht[uht];
        for(unsigned utau=0; utau<NumBins_tau; utau++){
          double tau_val = input_tau[utau];
          printf("%.0f %.2f %.2f %.2f\n", mT_val, rd_val, ht_val, tau_val);
          TString histo_name = TString::Format("hSource_%.0f_%.2f_%.2f_%.2f", mT_val, rd_val, ht_val, tau_val);
          TString histo_name_mt = TString::Format("GhettoFemto_mT_kstar_%.2f_%.2f_%.2f", rd_val, ht_val, tau_val);
          InputFile.cd();
          TH1F* hSource = (TH1F*)InputFile.Get(histo_name);
          TH2F* hMtKstar = NULL;
          if(umT==0) {hMtKstar = (TH2F*)InputFile.Get(histo_name_mt);}
          hSource->Scale(1./hSource->Integral(),"width");
          OutputFile.cd();
          //TH1F* hSourceCat = (TH1F*)hSource->Clone(TString::Format("hSourceCat_%.0f_%.2f_%.2f_%.2f", mT_val, rd_val, ht_val, tau_val));
          DLM_Histo<float>* dlmSource = Convert_TH1F_DlmHisto(hSource);
          double Chi2=0;
          KdpPars my_kdp;
          OutputFile.cd();
          //mode 2 or 3?
          TF1* fit_ptr = fit_source_kdp(hSource, my_kdp, Chi2, 3,96,2);
          unsigned WhichBin[4] = {umT,urd,uht,utau};
          dlm_kdp_pars.SetBinContent(WhichBin,my_kdp);
          fit_ptr->SetName(TString::Format("fSource_%.0f_%.2f_%.2f_%.2f", mT_val, rd_val, ht_val, tau_val));
          OutputFile.cd();
          if(hMtKstar) hMtKstar->Write();
          hSource->Write();
          fit_ptr->Write();
          
          TH1F* hCk = new TH1F(TString::Format("hCk_%.0f_%.2f_%.2f_%.2f", mT_val, rd_val, ht_val, tau_val),
                  TString::Format("hCk_%.0f_%.2f_%.2f_%.2f", mT_val, rd_val, ht_val, tau_val),
                  kSteps, kMin, kMax);

          DLM_KdpSource my_kdp_source(my_kdp);
          //
          double dummy_r = 2;
          double dummy_kxc[3] = {0,dummy_r,0};
          DLM_HistoSource dlmSRC(dlmSource);


  CATS Kitty_Dummy;
  Kitty_Dummy.SetMomBins(kSteps,kMin,kMax);
  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles/",GetCernBoxDimi()));
  AnalysisObject.SetUpCats_pp(Kitty_Dummy,"AV18","",0,0);          
          //printf("S(2) = %.3f vs %.3f vs %.3f\n",fit_ptr->Eval(dummy_r),my_kdp_source.Eval(&dummy_r), dlmSRC.Eval(dummy_kxc));
          //Kitty_Dummy.SetAnaSource(CatsSourceForwarder, &dlmSRC, 0);
          Kitty_Dummy.SetAnaSource(CatsSourceForwarder, &my_kdp_source, 0);
          Kitty_Dummy.SetUseAnalyticSource(true);
          Kitty_Dummy.SetAutoNormSource(false);
          Kitty_Dummy.SetNormalizedSource(true);
          Kitty_Dummy.SetNotifications(CATS::nError);
          //I am not sure if the source is always properly listed as not update, just in case we always do it
          //if(Kitty_Dummy.SourceStatus()==true){
          //  printf("Update source\n");
          //  Kitty_Dummy.KillTheCat(CATS::kSourceChanged);
          //}
          //else{
            //printf("Normal update");
            Kitty_Dummy.KillTheCat();

            //for(unsigned uRad=0; uRad<hSourceCat->GetNbinsX(); uRad++){
            //  double RAD = hSourceCat->GetBinCenter(uRad+1);
            //  hSourceCat->SetBinContent(uRad+1, Kitty_Dummy.EvaluateTheSource(1,RAD,0));
            //}

          //}
          for(unsigned uMom=0; uMom<kSteps; uMom++){
            hCk->SetBinContent(uMom+1, Kitty_Dummy.GetCorrFun(uMom));
          }

          OutputFile.cd();
          //hSourceCat->Write();
          hCk->Write();
          
          
          if(urd==NumBins_rd/2 && uht==NumBins_ht/2 && utau==NumBins_tau/2){
            OutputFile.cd();
            for(unsigned uMom=0; uMom<kSteps; uMom++){
              double Ck_Mean = Kitty_Dummy.GetCorrFun(uMom);
              double Ck_Err = Ck_Mean*hErrors->GetBinError(uMom+1);
              if(uMom==0) Ck_Err = hErrors->GetBinError(uMom+1);
              double Ck_Val = -1;
              while(Ck_Val<0){
                Ck_Val = rangen.Gaus(Ck_Mean, Ck_Err);
              }
            //  printf("%u %.3f %.3f\n", uMom, Ck_Mean, Ck_Err);
              hCk_Synthetic[umT]->SetBinContent(uMom+1, Ck_Val);
              hCk_Synthetic[umT]->SetBinError(uMom+1, Ck_Err);
            }            
          }

          delete dlmSource;
          delete fit_ptr;
          delete hCk;
        }
      }
    }
    //hCk_Synthetic[umT]->Write();
  }

  OutputFile.cd();
  for(unsigned umT=0; umT<NumBins_mT; umT++){
    hCk_Synthetic[umT]->Write();
    delete hCk_Synthetic[umT];
  }
  


  dlm_kdp_pars.QuickWrite(OutputFileNameDlm.Data(), true);

  delete [] BinRange_mT;
  delete [] BinRange_rd;
  delete [] BinRange_ht;
  delete [] BinRange_tau;

  delete [] BinCenter_mT;
  delete [] BinCenter_rd;
  delete [] BinCenter_ht;
  delete [] BinCenter_tau; 

  delete [] hCk_Synthetic;

    //CATS Kitty_Dummy;
    //Kitty_Dummy.SetMomBins(kSteps,kMin,kMax);
    //DLM_CommonAnaFunctions AnalysisObject;
    //AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles/",GetCernBoxDimi()));
    //AnalysisObject.SetUpCats_pp(Kitty_Dummy,"AV18","",0,0);
    //DLM_HistoSource src_CECA(HistoSourceR);
    //Kitty_Dummy.SetAnaSource(CatsSourceForwarder, &src_CECA, 0);
    //Kitty_Dummy.SetUseAnalyticSource(true);
    //Kitty_Dummy.SetAutoNormSource(false);
    //Kitty_Dummy.SetNormalizedSource(true);

}


CATS* MC_closure_Cat=NULL;
DLM_Histo<KdpPars>* MC_closure_kdp = NULL;


//create custom source function HERE


//pars[0,1,2] = rd,ht,tau
double MC_closure_CK_for_Max_FF(double* x, double* pars){
  double where_to_eval[3] = {pars[0],pars[1],pars[2]};
  KdpPars current_kdp = MC_closure_kdp->Eval(where_to_eval);
  MC_closure_Cat->SetAnaSource(CatsSourceForwarder, &current_kdp, 0);
  MC_closure_Cat->KillTheCat();
  return MC_closure_Cat->EvalCorrFun(*x);
}

//here we read the KDP file and finally make and fit the correlation
//as a start, make QA to see if you can read everything alright from the kdp file
void MC_closure_CK_for_Max(){
/*
  TString BaseFileName = TString::Format("%s/FunWithCeca/MC_closure_for_Max/Official_v1/pp_S23_TY253K",GetFemtoOutputFolder());
  TString InputFileName = TString::Format(BaseFileName+"_kdp.root",GetFemtoOutputFolder());
  TString OutputFileName = TString::Format(BaseFileName+"_kdp_fits.root",GetFemtoOutputFolder());
  TString InputFileNameDlm = TString::Format(BaseFileName+"_kdp.dlm",GetFemtoOutputFolder());

  TFile fOutput(OutputFileName, "recreate");
  TFile fInput(InputFileName, "read");
  

  //DLM_Histo<KdpPars> dlm_kdp_pars;
  if(!MC_closure_kdp) MC_closure_kdp = new DLM_Histo<KdpPars>();
  MC_closure_kdp->QuickLoad(InputFileNameDlm.Data());

  unsigned NumBins_mT = dlm_kdp_pars.GetNbins(0);
  unsigned NumBins_rd = dlm_kdp_pars.GetNbins(1);
  unsigned NumBins_ht = dlm_kdp_pars.GetNbins(2);
  unsigned NumBins_tau = dlm_kdp_pars.GetNbins(3);

  DLM_MultiFit GrandeFitter;

  TH1F** hToFit = new TH1F*[NumBins_mT];
  TF1** fitFun = new TF1*[NumBins_mT];
  CATS* Kitty = new CATS[NumBins_mT];
  for(unsigned umT=0; umT<NumBins_mT; umT++){
    double mT_val = dlm_kdp_pars.GetBinCenter(0,umT);
    fInput.cd();
    fInput.ls();
    TH1F* hTemp = (TH1F*)fInput.Get(TString::Format("hCk_Synthetic_mT%.0f",mT_val));
    double kMin = hTemp->GetXaxis()->GetBinLowEdge(1);
    double kMax = hTemp->GetXaxis()->GetBinUpEdge(hTemp->GetNbinsX());
    fOutput.cd();
    hToFit[umT] = (TH1F*)hTemp->Clone(TString::Format("hCk_Fit_mT%.0f",mT_val));
    fitFun[umT] = new TF1(TString::Format("fitFun_mT%.0f",mT_val),MC_closure_CK_for_Max_FF,kMin,kMax,3);
    //printf("add %u\n",umT);
    GrandeFitter.AddSpectrum(hToFit[umT],fitFun[umT]);
    //printf("done\n");
  }

  for(unsigned umT=1; umT<NumBins_mT; umT++){
    GrandeFitter.SetEqualPar(umT,1,0,1);//d
    GrandeFitter.SetEqualPar(umT,2,0,2);//ht
    GrandeFitter.SetEqualPar(umT,3,0,3);//hz or tau
  } 

  ROOT::Fit::FitResult fitresptr = GrandeFitter.PerformGlobalFit(true);



  for(unsigned umT=0; umT<NumBins_mT; umT++){
    delete hToFit[umT];
    delete fitFun[umT];
  }
  delete [] hToFit;
  delete [] fitFun;
  delete [] Kitty;
*/

/*
  for(unsigned umT=0; umT<dlm_kdp_pars.GetNbins(0); umT++){
    double mT_val = dlm_kdp_pars.GetBinCenter(0,umT);
    for(unsigned urd=0; urd<dlm_kdp_pars.GetNbins(1); urd++){
      double rd_val = dlm_kdp_pars.GetBinCenter(1,urd);
      for(unsigned uht=0; uht<dlm_kdp_pars.GetNbins(2); uht++){
        double ht_val = dlm_kdp_pars.GetBinCenter(2,uht);
        for(unsigned utau=0; utau<dlm_kdp_pars.GetNbins(3); utau++){
          double tau_val = dlm_kdp_pars.GetBinCenter(3,utau);
          //printf("%.0f %.2f %.2f %.2f\n", mT_val, rd_val, ht_val, tau_val);
          fOutput.cd();
          
          TH1F* hSource = new TH1F(TString::Format("hKDP_%.0f %.2f %.2f %.2f", mT_val, rd_val, ht_val, tau_val), 
                        TString::Format("hKDP_%.0f %.2f %.2f %.2f", mT_val, rd_val, ht_val, tau_val),
                        2048,0,32);
          for(unsigned uRad=0; uRad<hSource->GetNbinsX(); uRad++){
            double RAD = hSource->GetBinCenter(uRad+1);
            unsigned WhichBin[4] = {umT,urd,uht,utau};
            //printf("uRad = %.3f -> %.3f\n", RAD, PoissonSum(RAD, dlm_kdp_pars.GetBinContent(WhichBin))); usleep(10e3);
            hSource->SetBinContent(uRad+1, PoissonSum(RAD, dlm_kdp_pars.GetBinContent(WhichBin)));
          }
          fOutput.cd();
          hSource->SetLineColor(kGreen);
          hSource->SetLineWidth(3);
          




          //hSource->Write();
          delete hSource;
        }
      }
    }
  }  
*/
  //double PoissonSum(const double& xVal, const KdpPars& kdppars)
}

int FUN_WITH_CECA(int argc, char *argv[]){
  //Sources_In_SourcePaper("pp",1);
  //Sources_In_SourcePaper("pL",1);
  //return 0;
//DLM_Random rangen(1);
//for(int i=0; i<20; i++) printf("%i\n", rangen.Int(1));


  //MC_closure_for_Max(17); return 0;
  //MC_closure_FIT_for_Max(); return 0;
  //MC_closure_CK_for_Max(); return 0;


  //Test1();
  //Test2();
  //Test3();
  //TestingFractions();
  //TestDecay();
  //TestBoost();

  //CecaTest_0();
  //CecaTest_1();
  //CecaTest_pd_1();
  //CecaTest_p_pi_1();
  //CecaTest_p_pi_1(true);

  //CECA_sim_for_AI_pp_v0(atoi(argv[1]), atoi(argv[2]));
  //void CECA_sim_for_AI_pp_v1(int SEED, unsigned NUM_CPU, unsigned NUM_CONFIGS, unsigned K_YIELD_PER_CONFIG, int flag, TString OutputFolderName)
  //CECA_sim_for_AI_pp_v1(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), TString::Format("%s/FunWithCeca/CECA_sim_for_AI_pp_v3",GetFemtoOutputFolder())); return 0;
  //int SEED, unsigned NUM_CPU, unsigned NUM_CONFIGS, unsigned K_YIELD_PER_CONFIG, 
  //  float min_rd, float max_rd, float min_ht, float max_ht, float min_tau, float max_tau,
  //  int flag, TString OutputFolderName
  CECA_sim_for_AI_pp_v2(1,1,4,16,
    0,1,0,6,0,6,
    2,TString::Format("%s/FunWithCeca/CECA_sim_for_AI_pp_v3",GetFemtoOutputFolder())); return 0;
  //CECA_primoridal_disto(atoi(argv[1]), atoi(argv[2]));

//printf("GaussFromMean 2.97 = %f\n",GaussFromMean(2.97));

  //for(int TAU=0; TAU<=120; TAU+=5){
  //  if(TMath::Nint(TAU)%20==0) continue;
  //for(int TAU=0; TAU<=120; TAU+=20){
  //std::vector<double> tau_list = {0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5};
  //std::vector<double> tau_list = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0};
  //std::vector<double> tau_list = {2.25, 4.75};
  //Ceca_pd_1(0,-1402,"KdReso");
  //for(double TAU : tau_list){
    //Ceca_pd_1(TAU,true,"Kd");
  //  Ceca_pd_1(TAU,false,"Kd");
    //Ceca_pd_1(TAU,-1408,"Kd");
    //Ceca_pd_1(TAU,-1408,"KdReso");
  //}return 0;
  //Ceca_pd_1(0,true,"pd");
  //Ceca_pd_1(30,true,"pd");
  //Ceca_pd_1(60,true,"pd");
  //Ceca_pd_1(90,true,"pd");
  //Ceca_pd_1(3,true,"pd");
  //Ceca_pd_1(0.0,1,"Dpi");
  //Ceca_pd_1(0.0,-1,"Dpi");
//Ceca_pd_1(0.0,-3,"Dpi");
  //Ceca_pd_1(0.0,1,"pd");
  //Ceca_pd_1(0.0,-1,"pd");
  //Ceca_pd_1(0.0,1,"pp");
  //Ceca_pd_1(0.0,-2,"pp");
  //Ceca_pd_1(0.0,-1402,"pL");//usm fit
  //Ceca_pd_1(0.0,-1408,"pL");//nlo19
  //return 0;

//Ceca_pd_1(0.0,-3,"pp");
  //Ceca_pd_1(0.0,1,"pp");
  //Ceca_pd_1(0.0,-1,"pp");
  //Ceca_pd_1(0.0,-1001,"pp");
  //Ceca_pd_1(0.0,-1002,"pp");
  //Ceca_pd_1(0.0,-1010,"pp");
  //Ceca_pd_1(0.0,-1020,"pp");
  //Ceca_pd_1(0.0,-1100,"pp");
  //Ceca_pd_1(0.0,-1200,"pp");
  //Ceca_pd_1(0.0,-900,"pp");
  //Ceca_pd_1(0.0,atoi(argv[1]),"pp");

  //Ceca_pd_1(0.0,-1402,"pd");
  //Ceca_pd_1(0.0,-1402,"pp");


//TFile fOutput(,"recreate");
  //kstar_source_simple_1(1,"pp",12);
  //kstar_source_ck_1(TString::Format("%s/FunWithCeca/kstar_source_simple_1/kstarSrc_pp_flag1_sd11and12.root",GetFemtoOutputFolder()));
  //kstar_source_ck_TEST1();

  //Ceca_pd_1(0.0,-1300,"pp");
  //Ceca_pd_1(0.0,-1,"pipi_core");
  //Ceca_pd_1(0.0,1,"pipi_core");

  //Ceca_pd_1(0.0,-1,"pd");
  //Ceca_pd_1(0.0,-1,"pP");
  //Ceca_pd_1(0.0,-1,"pipi");

  //Ceca_pd_1(0.0,-1,"pK");
  //Ceca_pd_1(0.0,1,"pK");
  //Ceca_pd_1(0.0,-3,"pK");
  //Ceca_pd_1(0.0,-1100,"pK");

//return 0;

/*
  Ceca_AB(211,211,1);
  Ceca_AB(212,212,1);
  Ceca_AB(211,212,1);

  Ceca_AB(221,221,1);
  Ceca_AB(222,222,1);
  Ceca_AB(221,222,1);

  Ceca_AB(231,231,1);
  Ceca_AB(232,232,1);
  Ceca_AB(231,232,1);

/////////////////////////////////
  Ceca_AB(211,211,3);

  Ceca_AB(212,212,3);
  Ceca_AB(211,212,3);

  Ceca_AB(221,221,3);
  Ceca_AB(222,222,3);
  Ceca_AB(221,222,3);

  Ceca_AB(231,231,3);
  Ceca_AB(232,232,3);
  Ceca_AB(231,232,3);
*/
//return 0;
  //Ceca_pd_1(65,true);
  //pd_rstar_vs_KstarTau_CreateTH2F("Kd");

  //Ceca_pK_1(false);
  //Ceca_pK_1(true);
  //LK_Source(1.07);
  //LK_Source(1.11);
  //LK_Source(1.15);

//const TString type, std::vector<float>& disp, std::vector<float>& hadr,
//  std::vector<float>& tdel, const float& TimeDel, const float TIMEOUT=60
  //std::vector<float> disp = {0.0,0.0,0.0};
  float hfluct = 0.25*0;
  std::vector<float> disp = {0.25,0.25,0.25};
  std::vector<float> hadr = {3.55,3.55,0.0};
  //delay of Dmesons,LightMesons,LightMesonReso
  std::vector<float> tdel = {0.0,0.0,0.0};
  float TimeDel = 2.3;
  //std::vector<float> disp = {0.82,0.82,0.82};
  //std::vector<float> hadr = {0.0,0.0,0.0};
  //Ceca_Dmeson_LightMeson_1("Dpi",disp,hadr,hfluct,tdel,TimeDel,20);
  //Ceca_Dmeson_LightMeson_1("DK",disp,hadr,tdel,TimeDel,20);

  //pd_paper_sources(true);
  //pd_paper_sources(false);

  //Ceca_mT_1();
  //Ceca_vs_RSM_1("Proton","Proton");
  //Ceca_vs_RSM_1("Proton","Lambda");
  //ceca_test();

//concentrate on this one, in particular the prop mother ceca based rsm, seems to almost work
  //Ceca_vs_RSM_2("Proton","Proton",true,1.3,1.1);
  //Ceca_vs_RSM_2("Proton","Lambda",true,1.43,1.1);
  //Ceca_vs_RSM_2("Proton","Proton",false);
  //Ceca_vs_RSM_2("Proton","Lambda",false);

  //Ceca_vs_RSM_2("Proton","Kaon",false);

  //HowWrongIsOurSource(0);
  //HowWrongIsOurSource(1);

  //Ceca_pp_pL_source();

  //Ghetto_Test1();

  //RandomParallel(1);
  //RandomParallel(4);

  //TestTime();

  ThreeBody_test1();
  //ThreeBody_test1_plots();

  return 0;
}
