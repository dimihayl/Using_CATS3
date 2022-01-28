
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

#include <boost/algorithm/string.hpp>

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

int FUN_WITH_CECA(int argc, char *argv[]){

//DLM_Random rangen(1);
//for(int i=0; i<20; i++) printf("%i\n", rangen.Int(1));

  //Test1();
  //Test2();
  Test3();
  //TestingFractions();
  //TestDecay();
  //TestBoost();

  //Ghetto_Test1();

  //RandomParallel(1);
  //RandomParallel(4);


  //TestTime();

  return 0;
}
