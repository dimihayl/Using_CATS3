
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

#include <iostream>
#include <unistd.h>
#include <omp.h>

#include "TREPNI.h"
#include "CATS.h"
#include "CATSconstants.h"
#include "CATStools.h"
#include "CECA.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TF1.h"
#include "TGenPhaseSpace.h"
#include "TRandom3.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TNtuple.h"

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
  dlmFemto_Ghetto_rstar->Rebin(4);
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

    prt->SetPtEtaPhi(Sigma*prt->GetMass());
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
    prt->SetPtEtaPhi(Sigma*prt->GetMass());
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
    prt->SetPtEtaPhi(Sigma*prt->GetMass());
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
    prt->SetPtEtaPhi(Sigma*prt->GetMass());
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

//0.97 -> 1.059
//we can use the same for pd pr Kd, type == pd or Kd
void Ceca_pd_1(const double& d_delay, const bool& EffFix, const TString type="pd"){

  if(type!="pd"&&type!="Kd"){
    printf("WHAT IS THIS: %s\n", type.Data());
    return;
  }
  const double HadronSize = 0.75*0;
  const double HadronSlope = 0.2*0;
  const double EtaCut = 0.8;
  const bool PROTON_RESO = true;
  const bool EQUALIZE_TAU = true;
  const double TIMEOUT = 180/6;
  //we run to either reproduce the core of 0.97,
  //or the upper limit of reff = 1.06+0.04
  //this leads to a 10% difference in the SP core source
  double rSP_core;
  if(type=="pd") rSP_core = EffFix?0.815*1.10:0.815;
  if(type=="Kd") rSP_core = EffFix?0.815*1.10:0.815;

  TREPNI Database(0);
  Database.SetSeed(11);
  std::vector<TreParticle*> ParticleList;
  ParticleList.push_back(Database.NewParticle("Proton"));
  ParticleList.push_back(Database.NewParticle("Kaon"));
  ParticleList.push_back(Database.NewParticle("Deuteron"));
  ParticleList.push_back(Database.NewParticle("Pion"));
  ParticleList.push_back(Database.NewParticle("Photon"));
  ParticleList.push_back(Database.NewParticle("ProtonReso"));
  ParticleList.push_back(Database.NewParticle("KaonReso"));
  ParticleList.push_back(Database.NewParticle("DeuteronReso"));

  TString OutputFolderName;
  if(EffFix) OutputFolderName = "FunWithCeca/Ceca_"+type+"_EffFix";
  else OutputFolderName = "FunWithCeca/Ceca_"+type+"_CoreFix";
  //OutputFolderName = "FunWithCeca";
  TFile fOutput(TString::Format("%s/%s/Ceca_%s_ET%i_PR%i_DD%.0f.root",
  GetFemtoOutputFolder(),OutputFolderName.Data(),type.Data(),EQUALIZE_TAU,PROTON_RESO,d_delay),"recreate");

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
  TString FileNameP1,FileNameAP1,FileNameP2,FileNameAP2;
  if(type=="pd"){
    FilePath = TString::Format("%s/CatsFiles/Source/CECA/pd/pT_spectra/",GetCernBoxDimi());
    FileNameP1 = "protonSpectra.root";
    FileNameAP1 = "AntiProtonSpectra.root";
    FileNameP2 = "deuteronSpectra.root";
    FileNameAP2 = "AntiDeuteronSpectra.root";
  }
  if(type=="Kd"){
    printf("NB The Kd path not set yet!!!\n");
    FilePath = TString::Format("%s/CatsFiles/Source/CECA/pd/pT_spectra/",GetCernBoxDimi());
    FileNameP1 = "protonSpectra.root";
    FileNameAP1 = "AntiProtonSpectra.root";
    FileNameP2 = "deuteronSpectra.root";
    FileNameAP2 = "AntiDeuteronSpectra.root";
  }

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

  dlm_pT_p = Convert_TH1F_DlmHisto(h_pT_p_all);
  dlm_pT_d = Convert_TH1F_DlmHisto(h_pT_d_all);
  dlm_pT_p->RescaleAxis(0,1000,false);
  dlm_pT_d->RescaleAxis(0,1000,false);

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

  dlm_pT_eta_d.SetUp(2);
  BinRange = dlm_pT_d->GetBinRange(0);
  dlm_pT_eta_d.SetUp(0,dlm_pT_d->GetNbins(),BinRange);
  delete [] BinRange;
  dlm_pT_eta_d.SetUp(1,1,-EtaCut,EtaCut);
  dlm_pT_eta_d.Initialize();
  for(unsigned uBin=0; uBin<dlm_pT_d->GetNbins(); uBin++){
    dlm_pT_eta_d.SetBinContent(uBin,0,dlm_pT_d->GetBinContent(uBin));
  }

  TH2F* hSampleQA_p = new TH2F("hSampleQA_p","hSampleQA_p",64,-1,8,64,-1,1);
  for(unsigned uIter=0; uIter<100*1000; uIter++){
    double axisValues[2];
    dlm_pT_eta_p.Sample(axisValues);
    hSampleQA_p->Fill(axisValues[0],axisValues[1]);
  }

  fOutput.cd();
  hSampleQA_p->Write();
  h_pT_p_all->Write();
  h_pT_d_all->Write();

  for(TreParticle* prt : ParticleList){
    if(prt->GetName()=="Proton"){
      prt->SetMass(Mass_p);
      if(type=="pd") prt->SetAbundance(35.78+64.22*(!PROTON_RESO));
      else prt->SetAbundance(0);
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);
      prt->SetPtEtaPhi(dlm_pT_eta_p);
    }
    else if(prt->GetName()=="Kaon"){
      prt->SetMass(Mass_Kch);
      if(type=="Kd") prt->SetAbundance(53.4+47.6*(!PROTON_RESO));
      else prt->SetAbundance(0);
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);
      prt->SetPtEtaPhi(dlm_pT_eta_p);
    }
    else if(prt->GetName()=="Deuteron"){
      prt->SetMass(Mass_d);
      prt->SetAbundance(100.0);
      prt->SetRadius(HadronSize*1.5);
      prt->SetRadiusSlope(HadronSlope);
      prt->SetDelayTau(d_delay);
      prt->SetPtEtaPhi(dlm_pT_eta_d);
    }
    else if(prt->GetName()=="Pion"){
      prt->SetMass(Mass_pic);
      prt->SetAbundance(0);
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);
    }
    else if(prt->GetName()=="ProtonReso"){
      prt->SetMass(1362);
      if(type=="pd") prt->SetAbundance(64.22*PROTON_RESO);
      else prt->SetAbundance(0);
      prt->SetWidth(hbarc/1.65);
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);
      prt->SetPtEtaPhi(dlm_pT_eta_p);

      prt->NewDecay();
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Proton"));
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
      prt->GetDecay(0)->SetBranching(100);
    }
    else if(prt->GetName()=="KaonReso"){
      prt->SetMass(1054);
      if(type=="Kd") prt->SetAbundance(47.6*PROTON_RESO);
      else prt->SetAbundance(0);
      prt->SetWidth(hbarc/3.66);
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);
      prt->SetPtEtaPhi(dlm_pT_eta_p);

      prt->NewDecay();
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Kaon"));
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
      prt->GetDecay(0)->SetBranching(100);
    }
  }

  std::vector<std::string> ListOfParticles;
  if(type=="pd") ListOfParticles.push_back("Proton");
  if(type=="Kd") ListOfParticles.push_back("Kaon");
  ListOfParticles.push_back("Deuteron");
//printf("%s %s\n", ListOfParticles.at(0).c_str(),ListOfParticles.at(1).c_str());
//return;
  CECA Ivana(Database,ListOfParticles);

  Ivana.SetDisplacement(rSP_core);
  Ivana.SetHadronization(0.0);
  Ivana.SetHadrFluctuation(0.0);
  Ivana.SetTau(0.0);
  Ivana.SetThermalKick(0);

  Ivana.SetTargetStatistics(10);
  Ivana.SetEventMult(2);
  Ivana.SetSourceDim(2);
  Ivana.SetDebugMode(true);
  Ivana.SetThreadTimeout(TIMEOUT);
  Ivana.EqualizeFsiTime(EQUALIZE_TAU);
  Ivana.GHETTO_EVENT = true;

  Ivana.GoBabyGo();
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
}


void pd_rstar_vs_KstarTau_CreateTH2F(){
  const TString BaseName1 = TString::Format("%s/FunWithCeca/Ceca_pd_CoreFix/Ceca_pd_ET1_PR1_",GetFemtoOutputFolder());
  const TString BaseName2 = TString::Format("%s/FunWithCeca/Ceca_pd_EffFix/Ceca_pd_ET1_PR1_",GetFemtoOutputFolder());
  const TString OutputFile = TString::Format("%s/FunWithCeca/Ceca_pd_rstar_vs_KstarTau_040322.root",GetFemtoOutputFolder());

  const unsigned NumMomBins = 60;
  const double kMin = 0;
  const double kMax = 800;
  const unsigned NumTauBins = 25;
  const double tMin = -2.5;
  const double tMax = 122.5;
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



double Get_reff(TH1F* hsource, const float lambda=1, const float CEI=0.9){
  TH1F* hfit4325 = (TH1F*)hsource->Clone("hfit4325");
  hfit4325->Scale(1./hfit4325->Integral(),"width");

  double lowerlimit;
  double upperlimit;
  GetCentralInterval(*hfit4325, CEI, lowerlimit, upperlimit, true);

  TF1* fit4325 = new TF1("fit4325","[0]*4.*TMath::Pi()*x*x*pow(4.*TMath::Pi()*[1]*[1],-1.5)*exp(-(x*x)/(4.*[1]*[1]))+1.-[0]",lowerlimit,upperlimit);
  fit4325->FixParameter(1,lambda);
  fit4325->SetParameter(1,hfit4325->GetMean()/2.3);
  fit4325->SetParLimits(1,hfit4325->GetMean()/10.,hfit4325->GetMean()*2.);

  hfit4325->Fit(fit4325,"Q, S, N, R, M");

  return fit4325->GetParameter(1);
}


void Ceca_vs_RSM_1(){
  const double reff = 1.32*1.00;
  const double HadronSize = 0.75;
  const double HadronSlope = 0.2;
  const double EtaCut = 0.8;
  const bool EQUALIZE_TAU = true;
  const double TIMEOUT = 60;
  const double TempPP = 171;//171 or 154
  const double Asym = 6;
  const double PancakeT = 2.09;
  const double PancakeZ = 0.42;
  const double Tau = 0.0;
  const double PancakeFluct = 0;//in % !!!
  const double DisplT = 1.05;
  const double DisplZ = 0.21;
  unsigned THREADS = 0;
  const double MomSpread = 850;
  const bool THERMAL_KICK = false;

  const double Sensitivity = 0.0033*10;

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
  const double Tau_Lambda = 4.69;

  const double Mass_ProtonReso = 1362;
  const double Mass_LambdaReso = 1462;

  const double FracProtonReso = 0.6422;
  const double FracLambdaReso = 0.6438;

  const unsigned NumRadBins = 256;
  const double rMin = 0;
  const double rMax = 32;



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
    }
    else if(prt->GetName()=="Lambda"){
      prt->SetMass(Mass_L);
      prt->SetAbundance((1.-FracLambdaReso));
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
      prt->SetMass(Mass_ProtonReso);
      prt->SetAbundance(FracProtonReso);
      prt->SetWidth(hbarc/Tau_Proton);
      prt->SetRadius(HadronSize);
      prt->SetRadiusSlope(HadronSlope);

      prt->NewDecay();
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Proton"));
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
      prt->GetDecay(0)->SetBranching(100);
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
    }

    //prt->SetPtEtaPhi(MomSpread*);
    prt->SetPtEtaPhi(prt->GetMass()*MomSpread*0.001);
  }
  std::vector<std::string> ListOfParticles;
  ListOfParticles.push_back("Proton");
  ListOfParticles.push_back("Lambda");
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
  Database.GetParticle("LambdaReso")->SetAbundance(0);
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

  TFile fOutput(TString::Format("%s/FunWithCeca/Ceca_vs_RSM_EqTau%i_Tau%.0f_Dsp%.0f-%.0f_Pnk%.0f-%.0f_Fct%.0f_Kick%.0f_Hds%.0f-%.0f.root",
  GetFemtoOutputFolder(),EQUALIZE_TAU,Tau*100.,Ivana.GetDisplacementT()*100.,Ivana.GetDisplacementZ()*100.,
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
  h_GhettoFemto_mT_core->Write();
  g_GhettoFemto_mT_core.Write();
  gMtCore_pp->Write();

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

void CecaTest_p_pi_1(const bool& CecaAngles=false){

  const double reff = 1.32*1.00;
  const double HadronSize = 0.75;
  const double HadronSlope = 0.2;
  const double EtaCut = 0.8;
  const bool EQUALIZE_TAU = true;
  const double TIMEOUT = 20;
  const double TempPP = 171;//171 or 154
  const double Asym = 6;
  const double PancakeT = 0.0;
  const double PancakeZ = 0.0;
  const double Tau = 0.0;
  const double PancakeFluct = 0;//in % !!!
  const double DisplT = 1.35;
  const double DisplZ = 0.25;
  unsigned THREADS = 0;
  const double MomSpread = 850;
  const bool THERMAL_KICK = false;

  const double Sensitivity = 0.0033*30;

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
      RanVal1 = RanGen.Exponential(fM1/(fP1*Tau1));
      if(!CecaAngles) MagicSource.AddBGT_RP(RanVal1,cos(Pi-AngleRcP2));
      //MagicSource.AddBGT_RP(1,0.5);
      //printf("RP %f : %f\n",RanVal1,cos(Pi-AngleRcP2));
      h_RP_AngleRcP1_Rsm->Fill(Pi-AngleRcP2);
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
      RanVal2 = RanGen.Exponential(fM2/(fP2*Tau2));
      if(!CecaAngles) MagicSource.AddBGT_PR(RanVal2,cos(Pi-AngleRcP1));
      //printf("PR %f : %f\n",RanVal2,cos(Pi-AngleRcP1));
      h_PR_AngleRcP2_Rsm->Fill(Pi-AngleRcP1);
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
      RanVal1 = RanGen.Exponential(fM1/(fP1*Tau1));
      RanVal2 = RanGen.Exponential(fM2/(fP2*Tau2));
      if(!CecaAngles) MagicSource.AddBGT_RR(RanVal1,cos(Pi-AngleRcP2),RanVal2,cos(Pi-AngleRcP1),cos(AngleP1P2));
      //MagicSource.AddBGT_RR(RanVal1,cos(Pi-AngleRcP2),RanVal2,cos(Pi-AngleRcP1),RanGen.Uniform(-1,1));
      //printf("RR %f : %f : %f : %f : %f\n",RanVal1,cos(Pi-AngleRcP2),RanVal2,cos(Pi-AngleRcP1),cos(Pi-AngleP1P2));
      //usleep(200e3);
      h_RR_AngleRcP1_Rsm->Fill(Pi-AngleRcP2);
      h_RR_AngleRcP2_Rsm->Fill(Pi-AngleRcP1);
      h_RR_AngleP1P2_Rsm->Fill(Pi-AngleP1P2);
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
    prt->SetPtEtaPhi(prt->GetMass()*MomSpread*0.001);
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
  if(CecaAngles){
    Ivana.SetUp_RSM = &MagicSource;
    Ivana.GoBabyGo(THREADS);
    Ivana.SetUp_RSM = NULL;
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
    //printf(" val = %f\n",val);
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

  TFile fOutput(TString::Format("%s/FunWithCeca/CecaTest_p_pi_EqTau%i_Tau%.0f_Dsp%.0f-%.0f_Pnk%.0f-%.0f_Fct%.0f_Kick%.0f_Hds%.0f-%.0f_%s.root",
  GetFemtoOutputFolder(),EQUALIZE_TAU,Tau*100.,Ivana.GetDisplacementT()*100.,Ivana.GetDisplacementZ()*100.,
  Ivana.GetHadronizationT()*100.,Ivana.GetHadronizationZ()*100.,Ivana.GetHadrFluctuation()*100.,
  ThKick*THERMAL_KICK,HadronSize*.100,HadronSlope*100.,CecaAngles?"CECA":"EPOS"),"recreate");

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
  h_GhettoFemto_mT_core->Write();
  g_GhettoFemto_mT_core.Write();
  gMtCore_pp->Write();

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


int FUN_WITH_CECA(int argc, char *argv[]){

//DLM_Random rangen(1);
//for(int i=0; i<20; i++) printf("%i\n", rangen.Int(1));

  //Test1();
  //Test2();
  //Test3();
  //TestingFractions();
  //TestDecay();
  //TestBoost();

  //CecaTest_0();
  //CecaTest_1();
  //CecaTest_pd_1();
  //CecaTest_p_pi_1(false);
  //CecaTest_p_pi_1(true);

//printf("GaussFromMean 2.97 = %f\n",GaussFromMean(2.97));
  for(int TAU=0; TAU<=120; TAU+=40){
    Ceca_pd_1(TAU,true,"pd");
    Ceca_pd_1(TAU,false,"pd");
  }
  //Ceca_pd_1(65,true);
  //pd_rstar_vs_KstarTau_CreateTH2F();

  //Ceca_mT_1();
  //Ceca_vs_RSM_1();
  //ceca_test();


  //Ghetto_Test1();

  //RandomParallel(1);
  //RandomParallel(4);


  //TestTime();

  return 0;
}
