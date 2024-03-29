
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
#include <thread>

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

#include "TMath.h"
#include "TF1.h"
#include "TCanvas.h"

#include "TLegend.h"

//#include "Math/SpecFuncMathMore.h"

unsigned NumThreads;

void MC1_SingleEvent(){

}

//action at two places
//the momentum of each particle is modified by adding the relative momentum modified by some factor,
//which depends on the relative momentum itself
//N.B. positive strength corresponds to attraction, 0 no interaction
float MC1_FSI_XY( const float& rel_mom,
                  const float& pos0, const float& strength0, const float& range0,
                  const float& pos1=0, const float& strength1=0, const float& range1=1){
  return strength0*exp(-pow((rel_mom-pos0)/range0,2.))+strength1*exp(-pow((rel_mom-pos1)/range1,2.));
}

//modifies the momenta (of single particles) by that factor, based on the proximity of other particles
//N.B. positive strength corresponds to attraction
float MC1_FSI_XYZ(  const float& Q3, const float& strength, const float& range){
  return strength*exp(-pow(Q3/range,2.));
}

float MC1_P(const float* particle){
  return sqrt(particle[0]*particle[0]+particle[1]*particle[1]+particle[2]*particle[2]);
}

float MC1_k(const float* p1, const float* p2){
  return sqrt(pow(p1[0]-p2[0],2.)+pow(p1[1]-p2[1],2.)+pow(p1[2]-p2[2],2.));
}

float MC1_Q3(const float* p1, const float* p2, const float* p3){
  return sqrt(pow(MC1_k(p1,p2),2.)+pow(MC1_k(p1,p3),2.)+pow(MC1_k(p2,p3),2.));
}

struct MC1_Particle {
  float xyz[3];
} ;

class MC1_MixBuffer{
public:
  MC1_MixBuffer(const unsigned& depth=1000):Depth(depth),ran_gen(0){
    Part = new float* [Depth];
    for(unsigned ud=0; ud<Depth; ud++){
      Part[ud] = new float [3];
    }
    Next = 0;
    NumEl = 0;
  }
  ~MC1_MixBuffer(){
    for(unsigned ud; ud<Depth; ud++){
      delete [] Part[ud];
    }
    delete [] Part;
  }

  void push_back(float* part){
    for(int i=0; i<3; i++){
      Part[Next][i]=part[i];
    }
    Next++;
    if(Next>=Depth)Next=0;
    if(NumEl<Depth)NumEl++;
  }

  void get_random(float* part){
    int index = ran_gen.Uniform(0,NumEl);
    for(int i=0; i<3; i++) part[i]=Part[index][i];
  }

  float fill_status(){
    return double(NumEl)/double(Depth);
  }

private:
  const unsigned Depth;
  DLM_Random ran_gen;
  float** Part;
  unsigned Next;
  unsigned NumEl;
};


//mode 0: the original, with 0 3-body (genuine)
//mode 1: the original, where I put some correlation based on kstar/Q3
//mode 10: same as 0, but only interaction between 12
//mode 11: same as 1, but only interaction between 12
//mode 2: only QS and no interaction
//mode 1213: we dont have any 3 body, and we only have the 12 and 13 two body (no 23)
void MC1(const unsigned NumEvents, const int Mode){
  const float p_sigma = 150*4;
  const unsigned NumMomBins = 128*1;
  const float MinMom = 0;
  const float MaxMom = p_sigma*4;
  const float MinQ3 = 0;
  const float MaxQ3 = p_sigma*6;
  const unsigned MixingDepth = 1000;
  const float QS_width = 100;
  unsigned CurrentDepth=0;

  //cases are:
  //no final state interaction (FSI)
  //only FSI between particles 1 and 2
  //only FSI between 1 and 3
  //only FSI between 2 and 3
  //FSI between all 2-body pairs
  //no 2-body, only genuine 3 body
  //all 2-body and genuine 3 body
  //the FSIM are like FSI_ALL, only that the unlisted particle
  //has been substituted with a random particle.
  //Still, the listed 2 particles stem from event where they had 2 and 3 body interaction
  //LSM is Laura Serksnyte Mixing, i.e. a set of 3 random particles taken from the mixed events,
  //where each of the 3 particles is taken from an FSI_ALL triplet
printf("Change this to 7 to run faster, or simply do not random generate >FSI_ALL\n");
  const unsigned NumScenario = 11;
  bool QuickMix = true;//if true, we only mix the FSI_ALL case
  const unsigned MaxThreads = std::thread::hardware_concurrency();
  enum Scenario { noFSI,  FSI_12, FSI_13, FSI_23, FSI_2B, FSI_3G, FSI_ALL
                          ,FSIM_12, FSIM_13, FSIM_23, LSM
                        };
  TString* StrSce = new TString [NumScenario];
  StrSce[0] = "noFSI";
  StrSce[1] = "FSI_12";
  StrSce[2] = "FSI_13";
  StrSce[3] = "FSI_23";
  StrSce[4] = "FSI_2B";
  StrSce[5] = "FSI_3G";
  StrSce[6] = "FSI_ALL";
  if(NumScenario>7) StrSce[7] = "FSIM_12";
  if(NumScenario>7) StrSce[8] = "FSIM_13";
  if(NumScenario>7) StrSce[9] = "FSIM_23";
  if(NumScenario>7) StrSce[10] = "LSM";

  //mixed event buffer
//I created MC1_MixBuffer, have a look
  //std::vector<MC1_Particle>** MB1 = new std::vector<MC1_Particle>* [MaxThreads];
  //std::vector<MC1_Particle>** MB2 = new std::vector<MC1_Particle>* [MaxThreads];
  //std::vector<MC1_Particle>** MB3 = new std::vector<MC1_Particle>* [MaxThreads];
  MC1_MixBuffer** MB1 = new MC1_MixBuffer* [MaxThreads];
  MC1_MixBuffer** MB2 = new MC1_MixBuffer* [MaxThreads];
  MC1_MixBuffer** MB3 = new MC1_MixBuffer* [MaxThreads];
  for(unsigned uTh=0; uTh<MaxThreads; uTh++){
    //MB1[uTh] = new std::vector<MC1_Particle> [NumScenario];
    //MB2[uTh] = new std::vector<MC1_Particle> [NumScenario];
    //MB3[uTh] = new std::vector<MC1_Particle> [NumScenario];
    MB1[uTh] = new MC1_MixBuffer [NumScenario];
    MB2[uTh] = new MC1_MixBuffer [NumScenario];
    MB3[uTh] = new MC1_MixBuffer [NumScenario];
  }

  DLM_Histo<float>** S1_p = new DLM_Histo<float>* [MaxThreads+1];
  DLM_Histo<float>** S2_p = new DLM_Histo<float>* [MaxThreads+1];
  DLM_Histo<float>** S3_p = new DLM_Histo<float>* [MaxThreads+1];
  DLM_Histo<float>** S12_k = new DLM_Histo<float>* [MaxThreads+1];
  DLM_Histo<float>** S13_k = new DLM_Histo<float>* [MaxThreads+1];
  DLM_Histo<float>** S23_k = new DLM_Histo<float>* [MaxThreads+1];
  DLM_Histo<float>** S_Q3 = new DLM_Histo<float>* [MaxThreads+1];
  for(unsigned uTh=0; uTh<MaxThreads+1; uTh++){
    S1_p[uTh] = new DLM_Histo<float> [NumScenario];
    S2_p[uTh] = new DLM_Histo<float> [NumScenario];
    S3_p[uTh] = new DLM_Histo<float> [NumScenario];
    S12_k[uTh] = new DLM_Histo<float> [NumScenario];
    S13_k[uTh] = new DLM_Histo<float> [NumScenario];
    S23_k[uTh] = new DLM_Histo<float> [NumScenario];
    S_Q3[uTh] = new DLM_Histo<float> [NumScenario];
    for(unsigned uSce=0; uSce<NumScenario; uSce++){
      S1_p[uTh][uSce].SetUp(1);
      S1_p[uTh][uSce].SetUp(0,NumMomBins,MinMom,MaxMom);
      S1_p[uTh][uSce].Initialize();

      S2_p[uTh][uSce].SetUp(1);
      S2_p[uTh][uSce].SetUp(0,NumMomBins,MinMom,MaxMom);
      S2_p[uTh][uSce].Initialize();

      S3_p[uTh][uSce].SetUp(1);
      S3_p[uTh][uSce].SetUp(0,NumMomBins,MinMom,MaxMom);
      S3_p[uTh][uSce].Initialize();

      S12_k[uTh][uSce].SetUp(1);
      S12_k[uTh][uSce].SetUp(0,NumMomBins,MinMom,MaxMom);
      S12_k[uTh][uSce].Initialize();

      S13_k[uTh][uSce].SetUp(1);
      S13_k[uTh][uSce].SetUp(0,NumMomBins,MinMom,MaxMom);
      S13_k[uTh][uSce].Initialize();

      S23_k[uTh][uSce].SetUp(1);
      S23_k[uTh][uSce].SetUp(0,NumMomBins,MinMom,MaxMom);
      S23_k[uTh][uSce].Initialize();

      S_Q3[uTh][uSce].SetUp(1);
      S_Q3[uTh][uSce].SetUp(0,NumMomBins,MinQ3,MaxQ3);
      S_Q3[uTh][uSce].Initialize();
    }
  }

  DLM_Random** RanGen = new DLM_Random* [MaxThreads];
  for(unsigned uTh=0; uTh<MaxThreads; uTh++){
    RanGen[uTh] = new DLM_Random(uTh+1);
  }

  float** Part1 = new float* [MaxThreads];
  float** Part2 = new float* [MaxThreads];
  float** Part3 = new float* [MaxThreads];
  //a mixed event particle
  //float** PartM1 = new float* [MaxThreads];
  //float** PartM2 = new float* [MaxThreads];
  //float** PartM3 = new float* [MaxThreads];

  float*** Part1_fsi = new float** [MaxThreads];
  float*** Part2_fsi = new float** [MaxThreads];
  float*** Part3_fsi = new float** [MaxThreads];

  for(unsigned uTh=0; uTh<MaxThreads; uTh++){
    Part1[uTh] = new float [3];
    Part2[uTh] = new float [3];
    Part3[uTh] = new float [3];
    //PartM1[uTh] = new float [3];
    //PartM2[uTh] = new float [3];
    //PartM3[uTh] = new float [3];

    Part1_fsi[uTh] = new float* [NumScenario];
    Part2_fsi[uTh] = new float* [NumScenario];
    Part3_fsi[uTh] = new float* [NumScenario];

    for(unsigned uSce=0; uSce<NumScenario; uSce++){
      Part1_fsi[uTh][uSce] = new float [3];
      Part2_fsi[uTh][uSce] = new float [3];
      Part3_fsi[uTh][uSce] = new float [3];
    }
  }

  long long NUM_EVENTS = NumEvents;
  NUM_EVENTS *= 1000000;
  #pragma omp parallel for
  for(long long uEvent=0; uEvent<NUM_EVENTS; uEvent++){
//printf("EV %u\n",uEvent);
    int ThId = omp_get_thread_num();
    static bool FirstIter=true;
    float p1;
    float p2;
    float p3;
    float k12;
    float k13;
    float k23;
    float Q3;

    bool GoodToGo = true;
    float Condition = Mode==2?RanGen[ThId]->Uniform():1;
    //do{
    for(int i=0; i<3; i++){
      Part1[ThId][i] = RanGen[ThId]->Gauss(0,p_sigma);
      Part2[ThId][i] = RanGen[ThId]->Gauss(0,p_sigma);
      Part3[ThId][i] = RanGen[ThId]->Gauss(0,p_sigma);
      if(FirstIter){
        //PartM1[ThId][i] = RanGen[ThId]->Gauss(0,p_sigma);
      }
    }
    p1 = MC1_P(Part1[ThId]);
    p2 = MC1_P(Part2[ThId]);
    p3 = MC1_P(Part3[ThId]);
    k12 = MC1_k(Part1[ThId],Part2[ThId]);
    k13 = MC1_k(Part1[ThId],Part3[ThId]);
    k23 = MC1_k(Part2[ThId],Part3[ThId]);
    Q3 = MC1_Q3(Part1[ThId],Part2[ThId],Part3[ThId]);

    //while(!GoodToGo);


    for(unsigned uSce=0; uSce<NumScenario; uSce++){
      for(int i=0; i<3; i++){
        Part1_fsi[ThId][uSce][i] = Part1[ThId][i];
        Part2_fsi[ThId][uSce][i] = Part2[ThId][i];
        Part3_fsi[ThId][uSce][i] = Part3[ThId][i];
      }
    }
    for(unsigned uSce=0; uSce<NumScenario; uSce++){
      if(Mode==0||Mode==1||Mode==10||Mode==11||Mode==1213){
        //say 1,2 are protons, 3 is pion
        if(uSce==FSI_12 || uSce==FSI_2B || uSce>=FSI_ALL){
          for(int i=0; i<3; i++){
            Part1_fsi[ThId][uSce][i] += (Part2[ThId][i]-Part1[ThId][i])*MC1_FSI_XY(k12,0,-0.6,100);//-0.4 and 80
            Part2_fsi[ThId][uSce][i] += (Part1[ThId][i]-Part2[ThId][i])*MC1_FSI_XY(k12,0,-0.6,100);
            //Part1_fsi[ThId][uSce][i] += (Part2[ThId][i]-Part1[ThId][i])*MC1_FSI_XY(k12,0,0.15,50,0,0.05,250);
            //Part2_fsi[ThId][uSce][i] += (Part1[ThId][i]-Part2[ThId][i])*MC1_FSI_XY(k12,0,0.15,50,0,0.05,250);
            //Part1_fsi[ThId][uSce][i] += (Part2_fsi[ThId][i]-Part1_fsi[ThId][i])*MC1_FSI_XY(k12,0,0.5,80);
            //Part2_fsi[ThId][uSce][i] += (Part1_fsi[ThId][i]-Part2_fsi[ThId][i])*MC1_FSI_XY(k12,0,0.5,80);
            //Part1_fsi[ThId][uSce][i] += (Part2_fsi[ThId][i]-Part1_fsi[ThId][i])*MC1_FSI_XY(k12,0,-0.2,120,20,0.1,10);
            //Part2_fsi[ThId][uSce][i] += (Part1_fsi[ThId][i]-Part2_fsi[ThId][i])*MC1_FSI_XY(k12,0,-0.2,120,20,0.1,10);
          }
        }
        if(uSce==FSI_13 || uSce==FSI_2B || uSce>=FSI_ALL){
          for(int i=0; i<3; i++){
            if(Mode==0||Mode==1||Mode==1213){
              //Part1_fsi[ThId][uSce][i] += (Part3[ThId][i]-Part1[ThId][i])*MC1_FSI_XY(k13,0,-0.5,40);
              //Part3_fsi[ThId][uSce][i] += (Part1[ThId][i]-Part3[ThId][i])*MC1_FSI_XY(k13,0,-0.5,40);
              Part1_fsi[ThId][uSce][i] += (Part3[ThId][i]-Part1[ThId][i])*MC1_FSI_XY(k13,0,-0.6,75);
              Part3_fsi[ThId][uSce][i] += (Part1[ThId][i]-Part3[ThId][i])*MC1_FSI_XY(k13,0,-0.6,75);
            }
            //Part1_fsi[ThId][uSce][i] += (Part3[ThId][i]-Part1[ThId][i])*MC1_FSI_XY(k13,0,0.3,400,320,0.1,40);//p pi-
            //Part3_fsi[ThId][uSce][i] += (Part1[ThId][i]-Part3[ThId][i])*MC1_FSI_XY(k13,0,0.3,400,320,0.1,40);
            //Part1_fsi[ThId][uSce][i] += (Part3[ThId][i]-Part1[ThId][i])*MC1_FSI_XY(k13,0,-0.05,300,500,0.05,250);//p pi+
            //Part3_fsi[ThId][uSce][i] += (Part1[ThId][i]-Part3[ThId][i])*MC1_FSI_XY(k13,0,-0.05,300,320,0.05,250);
            //Part1_fsi[ThId][uSce][i] += (Part3_fsi[ThId][i]-Part1_fsi[ThId][i])*MC1_FSI_XY(k13,0,-0.2,60);
            //Part3_fsi[ThId][uSce][i] += (Part1_fsi[ThId][i]-Part3_fsi[ThId][i])*MC1_FSI_XY(k13,0,-0.2,60);
            //Part1_fsi[ThId][uSce][i] += (Part3_fsi[ThId][i]-Part1_fsi[ThId][i])*MC1_FSI_XY(k13,0,-0.2,120,20,0.1,10);
            //Part3_fsi[ThId][uSce][i] += (Part1_fsi[ThId][i]-Part3_fsi[ThId][i])*MC1_FSI_XY(k13,0,-0.2,120,20,0.1,10);
          }
        }
        if(uSce==FSI_23 || uSce==FSI_2B || uSce>=FSI_ALL){
          for(int i=0; i<3; i++){
            if(Mode==0||Mode==1){
              //Part2_fsi[ThId][uSce][i] += (Part3[ThId][i]-Part2[ThId][i])*MC1_FSI_XY(k23,0,-0.3,80,200,-0.05,80);
              //Part3_fsi[ThId][uSce][i] += (Part2[ThId][i]-Part3[ThId][i])*MC1_FSI_XY(k23,0,-0.3,80,200,-0.05,80);
            Part2_fsi[ThId][uSce][i] += (Part3[ThId][i]-Part2[ThId][i])*MC1_FSI_XY(k23,0,-1.00,20);
            Part3_fsi[ThId][uSce][i] += (Part2[ThId][i]-Part3[ThId][i])*MC1_FSI_XY(k23,0,-1.00,20);
            //Part2_fsi[ThId][uSce][i] += (Part3[ThId][i]-Part2[ThId][i])*MC1_FSI_XY(k23,0,0.3,400,320,0.1,40);//p pi-
            //Part3_fsi[ThId][uSce][i] += (Part2[ThId][i]-Part3[ThId][i])*MC1_FSI_XY(k23,0,0.3,400,320,0.1,40);
            //Part2_fsi[ThId][uSce][i] += (Part3[ThId][i]-Part2[ThId][i])*MC1_FSI_XY(k23,0,-0.05,300,500,0.05,250);//p pi+
            //Part3_fsi[ThId][uSce][i] += (Part2[ThId][i]-Part3[ThId][i])*MC1_FSI_XY(k23,0,-0.05,300,320,0.05,250);
            //Part2_fsi[ThId][uSce][i] += (Part3[ThId][i]-Part2[ThId][i])*MC1_FSI_XY(k23,0,-0.4,80,200,-0.07,80);
            //Part3_fsi[ThId][uSce][i] += (Part2[ThId][i]-Part3[ThId][i])*MC1_FSI_XY(k23,0,-0.4,80,200,-0.07,80);
            //Part2_fsi[ThId][uSce][i] += (Part3_fsi[ThId][i]-Part2_fsi[ThId][i])*MC1_FSI_XY(k23,0,-0.4,80,200,-0.07,80);
            //Part3_fsi[ThId][uSce][i] += (Part2_fsi[ThId][i]-Part3_fsi[ThId][i])*MC1_FSI_XY(k23,0,-0.4,80,200,-0.07,80);
            //Part2_fsi[ThId][uSce][i] += (Part3_fsi[ThId][i]-Part2_fsi[ThId][i])*MC1_FSI_XY(k23,0,-0.2,120,20,0.1,10);
            //Part3_fsi[ThId][uSce][i] += (Part2_fsi[ThId][i]-Part3_fsi[ThId][i])*MC1_FSI_XY(k23,0,-0.2,120,20,0.1,10);
            }
          }
        }
        if( (Mode==1 || Mode==11 || Mode==2) && (uSce==FSI_3G || uSce>=FSI_ALL) ){
          for(int i=0; i<3; i++){
            Part1_fsi[ThId][uSce][i] -= Part1[ThId][i]*MC1_FSI_XYZ(Q3,0.15,300);//0.4
            Part2_fsi[ThId][uSce][i] -= Part2[ThId][i]*MC1_FSI_XYZ(Q3,0.15,300);
            Part3_fsi[ThId][uSce][i] -= Part3[ThId][i]*MC1_FSI_XYZ(Q3,0.15,300);
          }
        }
      }
      else if(Mode==2){
        if(uSce==FSI_12 || uSce==FSIM_12 || uSce==FSI_2B || uSce==FSI_3G || uSce>=FSI_ALL){
          GoodToGo *= Condition<=(1.-exp(-pow(k12/QS_width,2.)));
        }
        else if(uSce==FSI_13 || uSce==FSIM_13 || uSce==FSI_2B || uSce==FSI_3G || uSce>=FSI_ALL){
          GoodToGo *= Condition<=(1.-exp(-pow(k13/QS_width,2.)));
        }
        else if(uSce==FSI_23 || uSce==FSIM_23 || uSce==FSI_2B || uSce==FSI_3G || uSce>=FSI_ALL){
          GoodToGo *= Condition<=(1.-exp(-pow(k23/QS_width,2.)));
        }
      }
      else{
        printf("Wrong mode\n");
      }

      if(!GoodToGo){
        break;
      }

      //at this point, the particles are fully modifed by the interaction
      //now we do event mixing, where we switch one particle by an uncorrelated particle
      //from the mixed event. For the event mixing, the particles are also set after fully interacting
      const int MIX_SCE = FSI_ALL;
      if(uSce==FSIM_12){
        //printf("fill_status = %f\n",MB3[ThId][uSce].fill_status());
        if(MB3[ThId][MIX_SCE].fill_status()>0.5){
          for(int i=0; i<3; i++){
            Part1_fsi[ThId][uSce][i] = Part1_fsi[ThId][FSI_ALL][i];
            Part2_fsi[ThId][uSce][i] = Part2_fsi[ThId][FSI_ALL][i];
          }
          MB3[ThId][MIX_SCE].get_random(Part3_fsi[ThId][uSce]);
        }
        else{
          for(int i=0; i<3; i++){
            Part1_fsi[ThId][uSce][i] = 0;
            Part2_fsi[ThId][uSce][i] = 0;
            Part3_fsi[ThId][uSce][i] = 0;
          }
        }
      }

      if(uSce==FSIM_13){
        if(MB2[ThId][MIX_SCE].fill_status()>0.5){
          for(int i=0; i<3; i++){
            Part1_fsi[ThId][uSce][i] = Part1_fsi[ThId][FSI_ALL][i];
            Part3_fsi[ThId][uSce][i] = Part3_fsi[ThId][FSI_ALL][i];
          }
          MB2[ThId][MIX_SCE].get_random(Part2_fsi[ThId][uSce]);
        }
        else{
          for(int i=0; i<3; i++){
            Part1_fsi[ThId][uSce][i] = 0;
            Part2_fsi[ThId][uSce][i] = 0;
            Part3_fsi[ThId][uSce][i] = 0;
          }
        }
      }

      if(uSce==FSIM_23){
        if(MB1[ThId][MIX_SCE].fill_status()>0.5){
          for(int i=0; i<3; i++){
            Part2_fsi[ThId][uSce][i] = Part2_fsi[ThId][FSI_ALL][i];
            Part3_fsi[ThId][uSce][i] = Part3_fsi[ThId][FSI_ALL][i];
          }
          MB1[ThId][MIX_SCE].get_random(Part1_fsi[ThId][uSce]);
        }
        else{
          for(int i=0; i<3; i++){
            Part1_fsi[ThId][uSce][i] = 0;
            Part2_fsi[ThId][uSce][i] = 0;
            Part3_fsi[ThId][uSce][i] = 0;
          }
        }
      }

      if(uSce==LSM){
        if( MB1[ThId][MIX_SCE].fill_status()>0.5 &&
            MB2[ThId][MIX_SCE].fill_status()>0.5 &&
            MB3[ThId][MIX_SCE].fill_status()>0.5){
          MB1[ThId][MIX_SCE].get_random(Part1_fsi[ThId][uSce]);
          MB2[ThId][MIX_SCE].get_random(Part2_fsi[ThId][uSce]);
          MB3[ThId][MIX_SCE].get_random(Part3_fsi[ThId][uSce]);
        }
        else{
          for(int i=0; i<3; i++){
            Part1_fsi[ThId][uSce][i] = 0;
            Part2_fsi[ThId][uSce][i] = 0;
            Part3_fsi[ThId][uSce][i] = 0;
          }
        }
      }

      //if(uSce==FSIM_12){
        //for(int i=0; i<3; i++) Part3_fsi[ThId][uSce][i] = PartM1[ThId][i];
      //}
      //if(uSce==FSIM_13){
        //for(int i=0; i<3; i++) Part2_fsi[ThId][uSce][i] = PartM1[ThId][i];
      //}
      //if(uSce==FSIM_23){
        //for(int i=0; i<3; i++) Part1_fsi[ThId][uSce][i] = PartM1[ThId][i];
      //}
    }//uSce

    for(unsigned uSce=0; uSce<NumScenario; uSce++){
      double value;
      //we dont fill if we have tracks with zero momentum (dummy above for useless event)
      if(MC1_P(Part1_fsi[ThId][uSce])&&MC1_P(Part2_fsi[ThId][uSce])&&MC1_P(Part3_fsi[ThId][uSce])){
        value = MC1_P(Part1_fsi[ThId][uSce]);
        S1_p[ThId][uSce].AddAt(&value);

        value = MC1_P(Part2_fsi[ThId][uSce]);
        S2_p[ThId][uSce].AddAt(&value);

        value = MC1_P(Part3_fsi[ThId][uSce]);
        S3_p[ThId][uSce].AddAt(&value);

        value = MC1_k(Part1_fsi[ThId][uSce],Part2_fsi[ThId][uSce]);
        S12_k[ThId][uSce].AddAt(&value);

        value = MC1_k(Part1_fsi[ThId][uSce],Part3_fsi[ThId][uSce]);
        S13_k[ThId][uSce].AddAt(&value);

        value = MC1_k(Part2_fsi[ThId][uSce],Part3_fsi[ThId][uSce]);
        S23_k[ThId][uSce].AddAt(&value);

        value = MC1_Q3(Part1_fsi[ThId][uSce],Part2_fsi[ThId][uSce],Part3_fsi[ThId][uSce]);
        S_Q3[ThId][uSce].AddAt(&value);
      }

    }

    //event mixing
    for(unsigned uSce=0; uSce<NumScenario; uSce++){
      if(QuickMix&&uSce!=FSI_ALL) continue;
      MC1_Particle tmp;

      int WhichPart = RanGen[ThId]->Integer(0,3);

//we need to replace the first element
//perhaps a small class to do that?
//we can have fixed amount of particles, and an iterator that saves the info
//about which is the newest and oldest member
      if(WhichPart==0){
        for(int i=0; i<3; i++) tmp.xyz[i]=Part1_fsi[ThId][uSce][i];
        //if(MB1[ThId][uSce].size()<MixingDepth)
        MB1[ThId][uSce].push_back(Part1_fsi[ThId][uSce]);
      }
      else if(WhichPart==1){
        //for(int i=0; i<3; i++) tmp.xyz[i]=Part2_fsi[ThId][uSce][i];
        MB2[ThId][uSce].push_back(Part2_fsi[ThId][uSce]);
      }
      else if(WhichPart==2){
        //for(int i=0; i<3; i++) tmp.xyz[i]=Part3_fsi[ThId][uSce][i];
        MB3[ThId][uSce].push_back(Part3_fsi[ThId][uSce]);
      }
      else printf("????!!!!!?????!!!!!\n");
    }

    //PartM1[ThId][0] = Part1[ThId][0];
    //PartM1[ThId][1] = Part1[ThId][1];
    //PartM1[ThId][2] = Part1[ThId][2];
    FirstIter = false;
  }//the parallel loop

printf("Out of loop\n");
  for(unsigned uTh=0; uTh<MaxThreads; uTh++){
    for(unsigned uSce=0; uSce<NumScenario; uSce++){
      S1_p[MaxThreads][uSce] += S1_p[uTh][uSce];
      S2_p[MaxThreads][uSce] += S2_p[uTh][uSce];
      S3_p[MaxThreads][uSce] += S3_p[uTh][uSce];
      S12_k[MaxThreads][uSce] += S12_k[uTh][uSce];
      S13_k[MaxThreads][uSce] += S13_k[uTh][uSce];
      S23_k[MaxThreads][uSce] += S23_k[uTh][uSce];
      S_Q3[MaxThreads][uSce] += S_Q3[uTh][uSce];
    }
  }
  for(unsigned uSce=0; uSce<NumScenario; uSce++){
    S1_p[MaxThreads][uSce].ComputeError();
    S2_p[MaxThreads][uSce].ComputeError();
    S3_p[MaxThreads][uSce].ComputeError();
    S12_k[MaxThreads][uSce].ComputeError();
    S13_k[MaxThreads][uSce].ComputeError();
    S23_k[MaxThreads][uSce].ComputeError();
    S_Q3[MaxThreads][uSce].ComputeError();
  }


  DLM_Histo<float> C12_exclusive = (S12_k[MaxThreads][FSI_12]/S12_k[MaxThreads][noFSI]);
  DLM_Histo<float> C12_open = (S12_k[MaxThreads][FSI_2B]/S12_k[MaxThreads][noFSI]);
  DLM_Histo<float> C12_orgy = (S12_k[MaxThreads][FSI_ALL]/S12_k[MaxThreads][noFSI]);
  DLM_Histo<float> C12_mix = C12_orgy;
  if(NumScenario>FSIM_12) C12_mix = (S12_k[MaxThreads][FSIM_12]/S12_k[MaxThreads][noFSI]);
  DLM_Histo<float> C12_LSM = C12_orgy;
  if(NumScenario>FSIM_12) C12_LSM = (S12_k[MaxThreads][FSIM_12]/S12_k[MaxThreads][LSM]);

  DLM_Histo<float> C13_exclusive = (S13_k[MaxThreads][FSI_13]/S13_k[MaxThreads][noFSI]);
  DLM_Histo<float> C13_open = (S13_k[MaxThreads][FSI_2B]/S13_k[MaxThreads][noFSI]);
  DLM_Histo<float> C13_orgy = (S13_k[MaxThreads][FSI_ALL]/S13_k[MaxThreads][noFSI]);
  DLM_Histo<float> C13_mix = C13_orgy;
  if(NumScenario>FSIM_13) C13_mix = (S13_k[MaxThreads][FSIM_13]/S13_k[MaxThreads][noFSI]);
  DLM_Histo<float> C13_LSM = C13_orgy;
  if(NumScenario>FSIM_13) C13_LSM = (S13_k[MaxThreads][FSIM_13]/S13_k[MaxThreads][LSM]);

  DLM_Histo<float> C23_exclusive = (S23_k[MaxThreads][FSI_23]/S23_k[MaxThreads][noFSI]);
  DLM_Histo<float> C23_open = (S23_k[MaxThreads][FSI_2B]/S23_k[MaxThreads][noFSI]);
  DLM_Histo<float> C23_orgy = (S23_k[MaxThreads][FSI_ALL]/S23_k[MaxThreads][noFSI]);
  DLM_Histo<float> C23_mix = C23_orgy;
  if(NumScenario>FSIM_23) C23_mix = (S23_k[MaxThreads][FSIM_23]/S23_k[MaxThreads][noFSI]);
  DLM_Histo<float> C23_LSM = C23_orgy;
  if(NumScenario>FSIM_23) C23_LSM = (S23_k[MaxThreads][FSIM_23]/S23_k[MaxThreads][LSM]);

  DLM_Histo<float> C12_voyeur = (S_Q3[MaxThreads][FSI_12]/S_Q3[MaxThreads][noFSI]);
  DLM_Histo<float> C13_voyeur = (S_Q3[MaxThreads][FSI_13]/S_Q3[MaxThreads][noFSI]);
  DLM_Histo<float> C23_voyeur = (S_Q3[MaxThreads][FSI_23]/S_Q3[MaxThreads][noFSI]);


printf("basic\n");


  //this is the analytical computation
  //original
  DLM_Histo<float> C123_open_Kubo = (C12_voyeur+C13_voyeur+C23_voyeur)-2.;
  DLM_Histo<float> C123_open_Dimi = C12_voyeur*C13_voyeur*C23_voyeur;
  //DLM_Histo<float> C123_open_DimiM = C12_mixed*C13_mixed*C23_mixed;
  DLM_Histo<float> R1R2R3_R = (S_Q3[MaxThreads][FSI_12]+S_Q3[MaxThreads][FSI_13]+S_Q3[MaxThreads][FSI_23])/S_Q3[MaxThreads][FSI_ALL];
  R1R2R3_R /= 3.;


  DLM_Histo<float> C123_open = (S_Q3[MaxThreads][FSI_2B]/S_Q3[MaxThreads][noFSI]);
  DLM_Histo<float> C123_threesome = (S_Q3[MaxThreads][FSI_3G]/S_Q3[MaxThreads][noFSI]);
  DLM_Histo<float> C123_Cum = (S_Q3[MaxThreads][FSI_3G]/S_Q3[MaxThreads][noFSI])-1.;
  DLM_Histo<float> C123_DimiCum = (S_Q3[MaxThreads][FSI_3G]/S_Q3[MaxThreads][noFSI])-1.;
  DLM_Histo<float> C123_orgy = (S_Q3[MaxThreads][FSI_ALL]/S_Q3[MaxThreads][noFSI]);

  //here we assume the correct 3-body, just for the fun of it
  DLM_Histo<float> C123_KuboOrgy = C123_Cum+C123_open_Kubo;
  DLM_Histo<float> C123_KuboDimiOrgy = C123_Cum+C123_open_Kubo;
  DLM_Histo<float> C123_DimiOrgy = (C123_Cum+1.)*C12_voyeur*C13_voyeur*C23_voyeur;

  DLM_Histo<float> C123_mix = C123_orgy;
  if(NumScenario>FSIM_12) C123_mix = (S_Q3[MaxThreads][FSIM_12]/S_Q3[MaxThreads][noFSI]);
  DLM_Histo<float> C132_mix = C123_orgy;
  if(NumScenario>FSIM_13) C132_mix = (S_Q3[MaxThreads][FSIM_13]/S_Q3[MaxThreads][noFSI]);
  DLM_Histo<float> C231_mix = C123_orgy;
  if(NumScenario>FSIM_23) C231_mix = (S_Q3[MaxThreads][FSIM_23]/S_Q3[MaxThreads][noFSI]);

  DLM_Histo<float> C123_LSM = C123_orgy;
  if(NumScenario>FSIM_12) C123_LSM = (S_Q3[MaxThreads][FSIM_12]/S_Q3[MaxThreads][LSM]);
  DLM_Histo<float> C132_LSM = C123_orgy;
  if(NumScenario>FSIM_13) C132_LSM = (S_Q3[MaxThreads][FSIM_13]/S_Q3[MaxThreads][LSM]);
  DLM_Histo<float> C231_LSM = C123_orgy;
  if(NumScenario>FSIM_23) C231_LSM = (S_Q3[MaxThreads][FSIM_23]/S_Q3[MaxThreads][LSM]);

  DLM_Histo<float> C123_3Gx2B = C123_threesome*C123_open;


  DLM_Histo<float> C123_Cum_Kubo = C123_orgy-C12_voyeur-C13_voyeur-C23_voyeur+2.;
  DLM_Histo<float> C123_threesome_Kubo = C123_orgy-C12_voyeur-C13_voyeur-C23_voyeur+3.;
  DLM_Histo<float> C123_Cum_Dimi = C123_orgy/(C12_voyeur*C13_voyeur*C23_voyeur)-1.;
  DLM_Histo<float> C123_threesome_Dimi = C123_orgy/(C12_voyeur*C13_voyeur*C23_voyeur);

  //assume that each pair gets a modification of its distribution, that ONLY depends on Q3 and not on 2-body
  //correct for that by summing (with equal weight) the all XX(MIXED) and correcting their modification with respect
  //to their sum using 2-body only

  DLM_Histo<float> S_Q3_AvgMix = S_Q3[MaxThreads][FSI_ALL];
  if(NumScenario>FSIM_23) S_Q3_AvgMix = (S_Q3[MaxThreads][FSIM_12]+S_Q3[MaxThreads][FSIM_13]+S_Q3[MaxThreads][FSIM_23])/3.;
  DLM_Histo<float> S_Q3_Avg = (S_Q3[MaxThreads][FSI_12]+S_Q3[MaxThreads][FSI_13]+S_Q3[MaxThreads][FSI_23])/3.;
  DLM_Histo<float> C123_threesome_Dimitros = C123_orgy/(C12_voyeur*C13_voyeur*C23_voyeur)*S_Q3_AvgMix/S_Q3_Avg;


  //DLM_Histo<float> C123_Cum_Kubo3 = C123_threesome-C12_orgy-C13_orgy-C23_orgy-1.;
  DLM_Histo<float> C123_Cum_Kubo3 = C123_orgy-C12_orgy-C13_orgy-C23_orgy+2.;
  DLM_Histo<float> C123_Cum_Dimi3 = C123_orgy/(C12_orgy*C13_orgy*C23_orgy)-1.;

  DLM_Histo<float> C123_Cum_Dimi_3Gx2B = C123_orgy/(C123_open)-1.;

  DLM_Histo<float> C2B_KuboTrue = C123_orgy-C123_threesome+3.;
  DLM_Histo<float> C2B_DimiTrue = C123_orgy/C123_threesome;

  DLM_Histo<float> C2B_Kubo = C12_voyeur+C13_voyeur+C23_voyeur;
  DLM_Histo<float> C2B_Dimi = C12_voyeur*C13_voyeur*C23_voyeur;

  DLM_Histo<float> C2B_KuboOrgy = C12_orgy+C13_orgy+C23_orgy;
  DLM_Histo<float> C2B_DimiOrgy = C12_orgy*C13_orgy*C23_orgy;


  DLM_Histo<float> C2B_KuboMix = C123_mix+C132_mix+C231_mix;
  DLM_Histo<float> C2B_DimiMix = C123_mix*C132_mix*C231_mix;
printf("all\n");




//for the MEETING

DLM_Histo<float> C_2B = (S_Q3[MaxThreads][FSI_2B]/S_Q3[MaxThreads][noFSI]);
DLM_Histo<float> C_3G = (S_Q3[MaxThreads][FSI_3G]/S_Q3[MaxThreads][noFSI]);
DLM_Histo<float> C_3G_LSM = (S_Q3[MaxThreads][FSI_3G]/S_Q3[MaxThreads][LSM]);
DLM_Histo<float> C_3Gc = (S_Q3[MaxThreads][FSI_3G]/S_Q3[MaxThreads][noFSI])-1.;
DLM_Histo<float> C_ALL = (S_Q3[MaxThreads][FSI_ALL]/S_Q3[MaxThreads][noFSI]);
DLM_Histo<float> C_ALL_LSM = (S_Q3[MaxThreads][FSI_ALL]/S_Q3[MaxThreads][LSM]);
//DLM_Histo<float> C_3Gc_new = C_ALL-C_2B+2.;

DLM_Histo<float> C_2B_Dimi = C12_voyeur*C13_voyeur*C23_voyeur;
DLM_Histo<float> C_2B_DimiLSM = C123_LSM*C132_LSM*C231_LSM;
DLM_Histo<float> C_2B_Kubo = C12_voyeur+C13_voyeur+C23_voyeur-2.;

DLM_Histo<float> C_3G_Dimi = C_ALL/(C12_voyeur*C13_voyeur*C23_voyeur);
DLM_Histo<float> C_3G_DimiLSM = C_ALL_LSM/(C123_LSM*C132_LSM*C231_LSM);
//if this equals C_3G, the 3-body genuine factorizes with respect to the 2-body
DLM_Histo<float> C_3G_DimiTrue = C_ALL/C_2B;
DLM_Histo<float> C_3G_DimiTrueRatio = C_3G_DimiTrue/C_3G;

DLM_Histo<float> C_3Gx2B = C_3G*C_2B;

DLM_Histo<float> C_3Gc_Kubo = C_ALL-C12_voyeur-C13_voyeur-C23_voyeur+2.;
DLM_Histo<float> C_3Gc_KuboNew = (C_3Gc_Kubo+C12_voyeur+C13_voyeur+C23_voyeur-2.)/(C12_voyeur*C13_voyeur*C23_voyeur)-1.;
DLM_Histo<float> C_3G_Kubo = C_3Gc_Kubo+1;
DLM_Histo<float> C_3G_KuboNew = C_3Gc_KuboNew+1;

///////////////////////////////////
printf("meeting\n");




//printf("Creating file\n");
  TFile fOutput(TString::Format("%s/MCUM/2B_is_12_13/MC1_m%u_%uM.root",GetFemtoOutputFolder(),Mode,NumEvents),"recreate");
  TDirectory* Singles = fOutput.mkdir(TString::Format("Singles"));
  TDirectory* Pairs = fOutput.mkdir(TString::Format("Pairs"));
  TDirectory* Threesome = fOutput.mkdir(TString::Format("Threesome"));
  TDirectory* Ck = fOutput.mkdir(TString::Format("Ck"));
  TDirectory* CQ3 = fOutput.mkdir(TString::Format("CQ3"));
//printf("tdir\n");
  TH1F** h_S1_p = new TH1F* [NumScenario];
  TH1F** h_S2_p = new TH1F* [NumScenario];
  TH1F** h_S3_p = new TH1F* [NumScenario];
  TH1F** h_S12_k = new TH1F* [NumScenario];
  TH1F** h_S13_k = new TH1F* [NumScenario];
  TH1F** h_S23_k = new TH1F* [NumScenario];
  TH1F** h_S_Q3 = new TH1F* [NumScenario];
//printf("bh\n");
  TH1F* h_C12_exclusive = Convert_DlmHisto_TH1F(&C12_exclusive,"h_C12_exclusive");
//printf("h_C12_exclusive\n");
  TH1F* h_C12_open = Convert_DlmHisto_TH1F(&C12_open,"h_C12_open");
  TH1F* h_C12_orgy = Convert_DlmHisto_TH1F(&C12_orgy,"h_C12_orgy");
//printf("a\n");
  TH1F* h_C12_mix = Convert_DlmHisto_TH1F(&C12_mix,"h_C12_mix");
  TH1F* h_C12_LSM = Convert_DlmHisto_TH1F(&C12_LSM,"h_C12_LSM");
  TH1F* h_C13_exclusive = Convert_DlmHisto_TH1F(&C13_exclusive,"h_C13_exclusive");
  TH1F* h_C13_open = Convert_DlmHisto_TH1F(&C13_open,"h_C13_open");
  TH1F* h_C13_orgy = Convert_DlmHisto_TH1F(&C13_orgy,"h_C13_orgy");
//printf("b\n");
  TH1F* h_C13_mix = Convert_DlmHisto_TH1F(&C13_mix,"h_C13_mix");
  TH1F* h_C13_LSM = Convert_DlmHisto_TH1F(&C13_LSM,"h_C13_LSM");
  TH1F* h_C23_exclusive = Convert_DlmHisto_TH1F(&C23_exclusive,"h_C23_exclusive");
  TH1F* h_C23_open = Convert_DlmHisto_TH1F(&C23_open,"h_C23_open");
  TH1F* h_C23_orgy = Convert_DlmHisto_TH1F(&C23_orgy,"h_C23_orgy");
  TH1F* h_C23_mix = Convert_DlmHisto_TH1F(&C23_mix,"h_C23_mix");
  TH1F* h_C23_LSM = Convert_DlmHisto_TH1F(&C23_LSM,"h_C23_LSM");
  TH1F* h_C12_voyeur = Convert_DlmHisto_TH1F(&C12_voyeur,"h_C12_voyeur");
  TH1F* h_C13_voyeur = Convert_DlmHisto_TH1F(&C13_voyeur,"h_C13_voyeur");
  TH1F* h_C23_voyeur = Convert_DlmHisto_TH1F(&C23_voyeur,"h_C23_voyeur");
  //TH1F* h_C12_mixed = Convert_DlmHisto_TH1F(&C12_mixed,"h_C12_mixed");
  //TH1F* h_C13_mixed = Convert_DlmHisto_TH1F(&C13_mixed,"h_C13_mixed");
  //TH1F* h_C23_mixed = Convert_DlmHisto_TH1F(&C23_mixed,"h_C23_mixed");
  TH1F* h_C123_open = Convert_DlmHisto_TH1F(&C123_open,"h_C123_open");
  TH1F* h_C123_threesome = Convert_DlmHisto_TH1F(&C123_threesome,"h_C123_threesome");
  TH1F* h_C123_Cum = Convert_DlmHisto_TH1F(&C123_Cum,"h_C123_Cum");
  TH1F* h_C123_orgy = Convert_DlmHisto_TH1F(&C123_orgy,"h_C123_orgy");
  TH1F* h_C123_KuboOrgy = Convert_DlmHisto_TH1F(&C123_KuboOrgy,"h_C123_KuboOrgy");
  TH1F* h_C123_KuboDimiOrgy = Convert_DlmHisto_TH1F(&C123_KuboDimiOrgy,"h_C123_KuboDimiOrgy");
  TH1F* h_C123_DimiOrgy = Convert_DlmHisto_TH1F(&C123_DimiOrgy,"h_C123_DimiOrgy");
//printf("1\n");
  TH1F* h_C123_mix = Convert_DlmHisto_TH1F(&C123_mix,"h_C123_mix");
  TH1F* h_C132_mix = Convert_DlmHisto_TH1F(&C132_mix,"h_C132_mix");
  TH1F* h_C231_mix = Convert_DlmHisto_TH1F(&C231_mix,"h_C231_mix");
  TH1F* h_C123_LSM = Convert_DlmHisto_TH1F(&C123_LSM,"h_C123_LSM");
  TH1F* h_C132_LSM = Convert_DlmHisto_TH1F(&C132_LSM,"h_C132_LSM");
  TH1F* h_C231_LSM = Convert_DlmHisto_TH1F(&C231_LSM,"h_C231_LSM");
//printf("2\n");

  TH1F* h_C123_3Gx2B = Convert_DlmHisto_TH1F(&C123_3Gx2B,"h_C123_3Gx2B");

  TH1F* h_C123_open_Kubo = Convert_DlmHisto_TH1F(&C123_open_Kubo,"h_C123_open_Kubo");
  TH1F* h_C123_open_Dimi = Convert_DlmHisto_TH1F(&C123_open_Dimi,"h_C123_open_Dimi");
  //TH1F* h_C123_open_KuboM = Convert_DlmHisto_TH1F(&C123_open_KuboM,"h_C123_open_KuboM");
  //TH1F* h_C123_open_DimiM = Convert_DlmHisto_TH1F(&C123_open_DimiM,"h_C123_open_DimiM");
  TH1F* h_C123_Cum_Kubo = Convert_DlmHisto_TH1F(&C123_Cum_Kubo,"h_C123_Cum_Kubo");
  TH1F* h_C123_Cum_Dimi = Convert_DlmHisto_TH1F(&C123_Cum_Dimi,"h_C123_Cum_Dimi");
  TH1F* h_C123_threesome_Kubo = Convert_DlmHisto_TH1F(&C123_threesome_Kubo,"h_C123_threesome_Kubo");
  TH1F* h_C123_threesome_Dimi = Convert_DlmHisto_TH1F(&C123_threesome_Dimi,"h_C123_threesome_Dimi");
  TH1F* h_C123_threesome_Dimitros = Convert_DlmHisto_TH1F(&C123_threesome_Dimitros,"h_C123_threesome_Dimitros");
  TH1F* h_C123_Cum_Kubo3 = Convert_DlmHisto_TH1F(&C123_Cum_Kubo3,"h_C123_Cum_Kubo3");
  TH1F* h_C123_Cum_Dimi3 = Convert_DlmHisto_TH1F(&C123_Cum_Dimi3,"h_C123_Cum_Dimi3");
  TH1F* h_C123_Cum_Dimi_3Gx2B = Convert_DlmHisto_TH1F(&C123_Cum_Dimi_3Gx2B,"h_C123_Cum_Dimi_3Gx2B");
  TH1F* h_C2B_KuboTrue = Convert_DlmHisto_TH1F(&C2B_KuboTrue,"h_C2B_KuboTrue");
  TH1F* h_C2B_DimiTrue = Convert_DlmHisto_TH1F(&C2B_DimiTrue,"h_C2B_DimiTrue");
  TH1F* h_C2B_Kubo = Convert_DlmHisto_TH1F(&C2B_Kubo,"h_C2B_Kubo");
  TH1F* h_C2B_Dimi = Convert_DlmHisto_TH1F(&C2B_Dimi,"h_C2B_Dimi");
  TH1F* h_C2B_KuboOrgy = Convert_DlmHisto_TH1F(&C2B_KuboOrgy,"h_C2B_KuboOrgy");
  TH1F* h_C2B_DimiOrgy = Convert_DlmHisto_TH1F(&C2B_DimiOrgy,"h_C2B_DimiOrgy");
  TH1F* h_C2B_KuboMix = Convert_DlmHisto_TH1F(&C2B_KuboMix,"h_C2B_KuboMix");
  TH1F* h_C2B_DimiMix = Convert_DlmHisto_TH1F(&C2B_DimiMix,"h_C2B_DimiMix");
//printf("3\n");

  TH1F* h_C_2B = Convert_DlmHisto_TH1F(&C_2B,"h_C_2B");
  TH1F* h_C_3G = Convert_DlmHisto_TH1F(&C_3G,"h_C_3G");
  TH1F* h_C_3G_LSM = Convert_DlmHisto_TH1F(&C_3G_LSM,"h_C_3G_LSM");
  TH1F* h_C_3Gc = Convert_DlmHisto_TH1F(&C_3Gc,"h_C_3Gc");
  TH1F* h_C_ALL = Convert_DlmHisto_TH1F(&C_ALL,"h_C_ALL");
  TH1F* h_C_ALL_LSM = Convert_DlmHisto_TH1F(&C_ALL_LSM,"h_C_ALL_LSM");


  TH1F* h_C_2B_Dimi = Convert_DlmHisto_TH1F(&C_2B_Dimi,"h_C_2B_Dimi");
  TH1F* h_C_2B_Kubo = Convert_DlmHisto_TH1F(&C_2B_Kubo,"h_C_2B_Kubo");

  TH1F* h_C_3G_Dimi = Convert_DlmHisto_TH1F(&C_3G_Dimi,"h_C_3G_Dimi");
  TH1F* h_C_3G_DimiLSM = Convert_DlmHisto_TH1F(&C_3G_DimiLSM,"h_C_3G_DimiLSM");
  TH1F* h_C_3G_DimiTrue = Convert_DlmHisto_TH1F(&C_3G_DimiTrue,"h_C_3G_DimiTrue");
  TH1F* h_C_3G_DimiTrueRatio = Convert_DlmHisto_TH1F(&C_3G_DimiTrueRatio,"h_C_3G_DimiTrueRatio");


  TH1F* h_C_3Gx2B = Convert_DlmHisto_TH1F(&C_3Gx2B,"h_C_3Gx2B");
  TH1F* h_C_3Gc_Kubo = Convert_DlmHisto_TH1F(&C_3Gc_Kubo,"h_C_3Gc_Kubo");
  TH1F* h_C_3Gc_KuboNew = Convert_DlmHisto_TH1F(&C_3Gc_KuboNew,"h_C_3Gc_KuboNew");
  TH1F* h_C_3G_Kubo = Convert_DlmHisto_TH1F(&C_3G_Kubo,"h_C_3G_Kubo");
  TH1F* h_C_3G_KuboNew = Convert_DlmHisto_TH1F(&C_3G_KuboNew,"h_C_3G_KuboNew");

//printf("4\n");

  TH1F* h_R1R2R3_R = Convert_DlmHisto_TH1F(&R1R2R3_R,"h_R1R2R3_R");
  for(unsigned uSce=0; uSce<NumScenario; uSce++){
    h_S1_p[uSce] = Convert_DlmHisto_TH1F(&S1_p[MaxThreads][uSce],TString::Format("h_S1_p_%s",StrSce[uSce].Data()));
    h_S2_p[uSce] = Convert_DlmHisto_TH1F(&S2_p[MaxThreads][uSce],TString::Format("h_S2_p_%s",StrSce[uSce].Data()));
    h_S3_p[uSce] = Convert_DlmHisto_TH1F(&S3_p[MaxThreads][uSce],TString::Format("h_S3_p_%s",StrSce[uSce].Data()));
    h_S12_k[uSce] = Convert_DlmHisto_TH1F(&S12_k[MaxThreads][uSce],TString::Format("h_S12_k_%s",StrSce[uSce].Data()));
    h_S13_k[uSce] = Convert_DlmHisto_TH1F(&S13_k[MaxThreads][uSce],TString::Format("h_S13_k_%s",StrSce[uSce].Data()));
    h_S23_k[uSce] = Convert_DlmHisto_TH1F(&S23_k[MaxThreads][uSce],TString::Format("h_S23_k_%s",StrSce[uSce].Data()));
    h_S_Q3[uSce] = Convert_DlmHisto_TH1F(&S_Q3[MaxThreads][uSce],TString::Format("h_S_Q3_%s",StrSce[uSce].Data()));
  }

printf("Creating histos\n");
  fOutput.cd();
  Singles->cd();
  for(unsigned uSce=0; uSce<NumScenario; uSce++) Singles->Add(h_S1_p[uSce]);
  for(unsigned uSce=0; uSce<NumScenario; uSce++) Singles->Add(h_S2_p[uSce]);
  for(unsigned uSce=0; uSce<NumScenario; uSce++) Singles->Add(h_S3_p[uSce]);
  Singles->Write();
  Pairs->cd();
  for(unsigned uSce=0; uSce<NumScenario; uSce++) Pairs->Add(h_S12_k[uSce]);
  for(unsigned uSce=0; uSce<NumScenario; uSce++) Pairs->Add(h_S13_k[uSce]);
  for(unsigned uSce=0; uSce<NumScenario; uSce++) Pairs->Add(h_S23_k[uSce]);
  Pairs->Write();
  Threesome->cd();
  for(unsigned uSce=0; uSce<NumScenario; uSce++) Threesome->Add(h_S_Q3[uSce]);
  Threesome->Write();


  Ck->cd();
  Ck->Add(h_C12_exclusive);
  Ck->Add(h_C12_open);
  Ck->Add(h_C12_orgy);
  Ck->Add(h_C12_mix);
  Ck->Add(h_C12_LSM);
  Ck->Add(h_C13_exclusive);
  Ck->Add(h_C13_open);
  Ck->Add(h_C13_orgy);
  Ck->Add(h_C13_mix);
  Ck->Add(h_C13_LSM);
  Ck->Add(h_C23_exclusive);
  Ck->Add(h_C23_open);
  Ck->Add(h_C23_orgy);
  Ck->Add(h_C23_mix);
  Ck->Add(h_C23_LSM);
  Ck->Write();

  CQ3->cd();
  CQ3->Add(h_C12_voyeur);
  CQ3->Add(h_C13_voyeur);
  CQ3->Add(h_C23_voyeur);
  CQ3->Add(h_C123_LSM);
  CQ3->Add(h_C132_LSM);
  CQ3->Add(h_C231_LSM);
  CQ3->Add(h_C_2B);
  CQ3->Add(h_C_3G);
  CQ3->Add(h_C_3G_LSM);
  CQ3->Add(h_C_3Gc);
  CQ3->Add(h_C_ALL);
  CQ3->Add(h_C_ALL_LSM);
  CQ3->Add(h_C_2B_Dimi);
  CQ3->Add(h_C_2B_Kubo);
  CQ3->Add(h_C_3G_Dimi);
  CQ3->Add(h_C_3G_DimiLSM);
  CQ3->Add(h_C_3G_DimiTrue);
  CQ3->Add(h_C_3G_DimiTrueRatio);
  CQ3->Add(h_C_3Gx2B);
  CQ3->Add(h_C_3Gc_Kubo);
  CQ3->Add(h_C_3Gc_KuboNew);
  CQ3->Add(h_C_3G_Kubo);
  CQ3->Add(h_C_3G_KuboNew);

  /*
  CQ3->Add(h_C123_open_Kubo);
  CQ3->Add(h_C123_open_Dimi);
  CQ3->Add(h_R1R2R3_R);
  CQ3->Add(h_C123_open);
  CQ3->Add(h_C123_threesome);
  CQ3->Add(h_C123_threesome_Kubo);
  CQ3->Add(h_C123_threesome_Dimi);
  CQ3->Add(h_C123_threesome_Dimitros);
  CQ3->Add(h_C123_Cum);
  CQ3->Add(h_C123_Cum_Kubo);
  CQ3->Add(h_C123_Cum_Dimi);
  CQ3->Add(h_C123_Cum_Dimi_3Gx2B);
  CQ3->Add(h_C123_orgy);
  CQ3->Add(h_C123_KuboOrgy);
  CQ3->Add(h_C123_KuboDimiOrgy);
  CQ3->Add(h_C123_DimiOrgy);

  CQ3->Add(h_C123_3Gx2B);
  CQ3->Add(h_C123_mix);
  CQ3->Add(h_C132_mix);
  CQ3->Add(h_C231_mix);
  CQ3->Add(h_C2B_KuboTrue);
  CQ3->Add(h_C2B_DimiTrue);
  CQ3->Add(h_C2B_Kubo);
  CQ3->Add(h_C2B_Dimi);
  CQ3->Add(h_C2B_KuboOrgy);
  CQ3->Add(h_C2B_DimiOrgy);
  CQ3->Add(h_C2B_KuboMix);
  CQ3->Add(h_C2B_DimiMix);
  */

  CQ3->Write();

//printf("Writing histos\n");
  for(unsigned uTh=0; uTh<MaxThreads+1; uTh++){
    delete [] S1_p[uTh];
    delete [] S2_p[uTh];
    delete [] S3_p[uTh];
    delete [] S12_k[uTh];
    delete [] S13_k[uTh];
    delete [] S23_k[uTh];
    delete [] S_Q3[uTh];
  }
  delete [] S1_p;
  delete [] S2_p;
  delete [] S3_p;
  delete [] S12_k;
  delete [] S13_k;
  delete [] S23_k;
  delete [] S_Q3;
//printf("Deleted DLM histos\n");

  for(unsigned uSce=0; uSce<NumScenario; uSce++){
    delete h_S1_p[uSce];
    delete h_S2_p[uSce];
    delete h_S3_p[uSce];
    delete h_S12_k[uSce];
    delete h_S13_k[uSce];
    delete h_S23_k[uSce];
    delete h_S_Q3[uSce];
  }
  delete [] h_S1_p;
  delete [] h_S2_p;
  delete [] h_S3_p;
  delete [] h_S12_k;
  delete [] h_S13_k;
  delete [] h_S23_k;
  delete [] h_S_Q3;

  delete h_C12_exclusive;
  delete h_C12_open;
  delete h_C12_orgy;
  delete h_C12_mix;
  delete h_C13_exclusive;
  delete h_C13_open;
  delete h_C13_orgy;
  delete h_C13_mix;
  delete h_C23_exclusive;
  delete h_C23_open;
  delete h_C23_orgy;
  delete h_C23_mix;
  delete h_C12_voyeur;
  delete h_C13_voyeur;
  delete h_C23_voyeur;
  //delete h_C12_mixed;
  //delete h_C13_mixed;
  //delete h_C23_mixed;
  delete h_C123_open;
  delete h_C123_threesome;
  delete h_C123_Cum;
  delete h_C123_orgy;
  delete h_C123_KuboOrgy;
  delete h_C123_KuboDimiOrgy;
  delete h_C123_DimiOrgy;
  delete h_C123_mix;
  delete h_C132_mix;
  delete h_C231_mix;
  delete h_C123_3Gx2B;
  delete h_C123_open_Kubo;
  delete h_C123_open_Dimi;
  //delete h_C123_open_KuboM;
  //delete h_C123_open_DimiM;
  delete h_R1R2R3_R;
  delete h_C123_threesome_Kubo;
  delete h_C123_threesome_Dimi;
  delete h_C123_threesome_Dimitros;
  delete h_C123_Cum_Kubo;
  delete h_C123_Cum_Dimi;
  delete h_C123_Cum_Kubo3;
  delete h_C123_Cum_Dimi3;
  delete h_C123_Cum_Dimi_3Gx2B;
  delete h_C2B_KuboTrue;
  delete h_C2B_DimiTrue;
  delete h_C2B_Kubo;
  delete h_C2B_Dimi;
  delete h_C2B_KuboOrgy;
  delete h_C2B_DimiOrgy;
  delete h_C2B_KuboMix;
  delete h_C2B_DimiMix;
  delete h_C123_LSM;
  delete h_C132_LSM;
  delete h_C231_LSM;
  delete h_C12_LSM;
  delete h_C13_LSM;
  delete h_C23_LSM;

  delete h_C_2B;
  delete h_C_3G;
  delete h_C_3G_LSM;
  delete h_C_3Gc;
  delete h_C_ALL;
  delete h_C_ALL_LSM;
  delete h_C_2B_Dimi;
  delete h_C_2B_Kubo;
  delete h_C_3G_Dimi;
  delete h_C_3G_DimiLSM;
  delete h_C_3G_DimiTrue;
  delete h_C_3G_DimiTrueRatio;
  delete h_C_3Gx2B;
  delete h_C_3Gc_Kubo;
  delete h_C_3Gc_KuboNew;
  delete h_C_3G_Kubo;
  delete h_C_3G_KuboNew;

//printf("Deleted histos\n");
  for(unsigned uTh=0; uTh<MaxThreads; uTh++){
    delete RanGen[uTh];

    for(unsigned uSce=0; uSce<NumScenario; uSce++){
      delete [] Part1_fsi[uTh][uSce];
      delete [] Part2_fsi[uTh][uSce];
      delete [] Part3_fsi[uTh][uSce];
    }
    delete [] Part1_fsi[uTh];
    delete [] Part2_fsi[uTh];
    delete [] Part3_fsi[uTh];

    delete [] Part1[uTh];
    delete [] Part2[uTh];
    delete [] Part3[uTh];
    //delete [] PartM1[uTh];
  }
  delete [] RanGen;
  delete [] Part1_fsi;
  delete [] Part2_fsi;
  delete [] Part3_fsi;
  delete [] Part1;
  delete [] Part2;
  delete [] Part3;
  //delete [] PartM1;
  delete [] StrSce;

  delete Singles;
  delete Pairs;
  delete Threesome;
  delete Ck;
  delete CQ3;
  for(unsigned uTh=0; uTh<MaxThreads; uTh++){
    delete [] MB1[uTh];
    delete [] MB2[uTh];
    delete [] MB3[uTh];
  }
  delete [] MB1;
  delete [] MB2;
  delete [] MB3;

}

int MCUM(int argc, char *argv[]){
  NumThreads = std::thread::hardware_concurrency()==1?0:std::thread::hardware_concurrency()-1;
  omp_set_dynamic(0);
  omp_set_num_threads(NumThreads);
  printf("MCUM #%u\n",NumThreads);


  //MC1(atoi(argv[1]), atoi(argv[2]));
  //MC1(atoi(argv[1]), 1);
  //MC1(atoi(argv[1]), 0);
  //MC1(atoi(argv[1]), 10);
  //MC1(atoi(argv[1]), 11);
  MC1(atoi(argv[1]), 1213);

  //MC2(atoi(argv[1])*1000000);

  return 0;
}
