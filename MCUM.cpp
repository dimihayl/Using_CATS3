
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
    for(unsigned ud; ud<Depth; ud++){
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
    for(int i=0; i<3; i++)Part[Next][i]=part[i];
    Next++;
    if(Next>=Depth)Next=0;
    if(NumEl<Depth)NumEl++;
  }

  void get_random(float* part){
    int index = ran_gen.Uniform(0,NumEl);
    for(int i=0; i<3; i++) part[i]=Part[index][i];
  }

private:
  const unsigned Depth;
  DLM_Random ran_gen;
  float** Part;
  unsigned Next;
  unsigned NumEl;
};


void MC1(const int NumEvents){
  const float p_sigma = 200;
  const unsigned NumMomBins = 256;
  const float MinMom = 0;
  const float MaxMom = p_sigma*5;
  const float MinQ3 = 0;
  const float MaxQ3 = p_sigma*8;
  const unsigned MixingDepth = 1000;
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
printf("Change this to 7 to run faster, or simply do not random generate >FSI_ALL\n");
  const unsigned NumScenario = 7;
  const unsigned MaxThreads = std::thread::hardware_concurrency();
  enum Scenario { noFSI,  FSI_12, FSI_13, FSI_23, FSI_2B, FSI_3G, FSI_ALL
                          //,FSIM_12, FSIM_13, FSIM_23
                        };
  TString* StrSce = new TString [NumScenario];
  StrSce[0] = "noFSI";
  StrSce[1] = "FSI_12";
  StrSce[2] = "FSI_13";
  StrSce[3] = "FSI_23";
  StrSce[4] = "FSI_2B";
  StrSce[5] = "FSI_3G";
  StrSce[6] = "FSI_ALL";
  //if(NumScenario>7) StrSce[7] = "FSIM_12";
  //if(NumScenario>7) StrSce[8] = "FSIM_13";
  //if(NumScenario>7) StrSce[9] = "FSIM_23";

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

  #pragma omp parallel for
  for(unsigned uEvent=0; uEvent<NumEvents; uEvent++){
//printf("EV %u\n",uEvent);
    int ThId = omp_get_thread_num();
    static bool FirstIter=true;

    for(int i=0; i<3; i++){
      Part1[ThId][i] = RanGen[ThId]->Gauss(0,p_sigma);
      Part2[ThId][i] = RanGen[ThId]->Gauss(0,p_sigma);
      Part3[ThId][i] = RanGen[ThId]->Gauss(0,p_sigma);
      if(FirstIter){
        //PartM1[ThId][i] = RanGen[ThId]->Gauss(0,p_sigma);
      }
    }

    float p1 = MC1_P(Part1[ThId]);
    float p2 = MC1_P(Part2[ThId]);
    float p3 = MC1_P(Part3[ThId]);
    float k12 = MC1_k(Part1[ThId],Part2[ThId]);
    float k13 = MC1_k(Part1[ThId],Part3[ThId]);
    float k23 = MC1_k(Part2[ThId],Part3[ThId]);
    float Q3 = MC1_Q3(Part1[ThId],Part2[ThId],Part3[ThId]);
    for(unsigned uSce=0; uSce<NumScenario; uSce++){
      for(int i=0; i<3; i++){
        Part1_fsi[ThId][uSce][i] = Part1[ThId][i];
        Part2_fsi[ThId][uSce][i] = Part2[ThId][i];
        Part3_fsi[ThId][uSce][i] = Part3[ThId][i];
      }
    }
    for(unsigned uSce=0; uSce<NumScenario; uSce++){
      /*
      std::vector<float*> PARTICLES;
      if(uEvent%6==0){
        PARTICLES.push_back(Part1[ThId][i]);
        PARTICLES.push_back(Part2[ThId][i]);
        PARTICLES.push_back(Part3[ThId][i]);
      }
      else if(uEvent%6==1){
        PARTICLES.push_back(Part1[ThId][i]);
        PARTICLES.push_back(Part3[ThId][i]);
        PARTICLES.push_back(Part2[ThId][i]);
      }
      else if(uEvent%6==2){
        PARTICLES.push_back(Part2[ThId][i]);
        PARTICLES.push_back(Part1[ThId][i]);
        PARTICLES.push_back(Part3[ThId][i]);
      }
      else if(uEvent%6==3){
        PARTICLES.push_back(Part2[ThId][i]);
        PARTICLES.push_back(Part3[ThId][i]);
        PARTICLES.push_back(Part1[ThId][i]);
      }
      else if(uEvent%6==4){
        PARTICLES.push_back(Part3[ThId][i]);
        PARTICLES.push_back(Part1[ThId][i]);
        PARTICLES.push_back(Part2[ThId][i]);
      }
      else if(uEvent%6==5){
        PARTICLES.push_back(Part3[ThId][i]);
        PARTICLES.push_back(Part2[ThId][i]);
        PARTICLES.push_back(Part1[ThId][i]);
      }
      */

      if(uSce==FSI_12 || uSce==FSI_2B || uSce>=FSI_ALL){
        for(int i=0; i<3; i++){
          //Part1_fsi[ThId][uSce][i] += (Part2[ThId][i]-Part1[ThId][i])*MC1_FSI_XY(k12,0,0.5,80);
          //Part2_fsi[ThId][uSce][i] += (Part1[ThId][i]-Part2[ThId][i])*MC1_FSI_XY(k12,0,0.5,80);
          //Part1_fsi[ThId][uSce][i] += (Part2_fsi[ThId][i]-Part1_fsi[ThId][i])*MC1_FSI_XY(k12,0,0.5,80);
          //Part2_fsi[ThId][uSce][i] += (Part1_fsi[ThId][i]-Part2_fsi[ThId][i])*MC1_FSI_XY(k12,0,0.5,80);
          Part1_fsi[ThId][uSce][i] += (Part2[ThId][i]-Part1[ThId][i])*MC1_FSI_XY(k12,0,-0.2,120,20,0.1,10);
          Part2_fsi[ThId][uSce][i] += (Part1[ThId][i]-Part2[ThId][i])*MC1_FSI_XY(k12,0,-0.2,120,20,0.1,10);
        }
      }
      if(uSce==FSI_13 || uSce==FSI_2B || uSce>=FSI_ALL){
        for(int i=0; i<3; i++){
          //Part1_fsi[ThId][uSce][i] += (Part3[ThId][i]-Part1[ThId][i])*MC1_FSI_XY(k13,0,-0.2,60);
          //Part3_fsi[ThId][uSce][i] += (Part1[ThId][i]-Part3[ThId][i])*MC1_FSI_XY(k13,0,-0.2,60);
          //Part1_fsi[ThId][uSce][i] += (Part3_fsi[ThId][i]-Part1_fsi[ThId][i])*MC1_FSI_XY(k13,0,-0.2,60);
          //Part3_fsi[ThId][uSce][i] += (Part1_fsi[ThId][i]-Part3_fsi[ThId][i])*MC1_FSI_XY(k13,0,-0.2,60);
          Part1_fsi[ThId][uSce][i] += (Part3[ThId][i]-Part1[ThId][i])*MC1_FSI_XY(k13,0,-0.2,120,20,0.1,10);
          Part3_fsi[ThId][uSce][i] += (Part1[ThId][i]-Part3[ThId][i])*MC1_FSI_XY(k13,0,-0.2,120,20,0.1,10);
        }
      }
      if(uSce==FSI_23 || uSce==FSI_2B || uSce>=FSI_ALL){
        for(int i=0; i<3; i++){
          //Part2_fsi[ThId][uSce][i] += (Part3[ThId][i]-Part2[ThId][i])*MC1_FSI_XY(k23,0,0.4,80,200,-0.07,80);
          //Part3_fsi[ThId][uSce][i] += (Part2[ThId][i]-Part3[ThId][i])*MC1_FSI_XY(k23,0,0.4,80,200,-0.07,80);
          //Part2_fsi[ThId][uSce][i] += (Part3[ThId][i]-Part2[ThId][i])*MC1_FSI_XY(k23,0,-0.4,80,200,-0.07,80);
          //Part3_fsi[ThId][uSce][i] += (Part2[ThId][i]-Part3[ThId][i])*MC1_FSI_XY(k23,0,-0.4,80,200,-0.07,80);
          //Part2_fsi[ThId][uSce][i] += (Part3_fsi[ThId][i]-Part2_fsi[ThId][i])*MC1_FSI_XY(k23,0,-0.4,80,200,-0.07,80);
          //Part3_fsi[ThId][uSce][i] += (Part2_fsi[ThId][i]-Part3_fsi[ThId][i])*MC1_FSI_XY(k23,0,-0.4,80,200,-0.07,80);
          Part2_fsi[ThId][uSce][i] += (Part3[ThId][i]-Part2[ThId][i])*MC1_FSI_XY(k23,0,-0.2,120,20,0.1,10);
          Part3_fsi[ThId][uSce][i] += (Part2[ThId][i]-Part3[ThId][i])*MC1_FSI_XY(k23,0,-0.2,120,20,0.1,10);
        }
      }

      if(uSce==FSI_3G || uSce>=FSI_ALL){
        for(int i=0; i<3; i++){
          Part1_fsi[ThId][uSce][i] -= Part1[ThId][i]*MC1_FSI_XYZ(Q3,0.0,300);//0.4
          Part2_fsi[ThId][uSce][i] -= Part2[ThId][i]*MC1_FSI_XYZ(Q3,0.0,300);
          Part3_fsi[ThId][uSce][i] -= Part3[ThId][i]*MC1_FSI_XYZ(Q3,0.0,300);
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
    }



    for(unsigned uSce=0; uSce<NumScenario; uSce++){
      double value;

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

    //event mixing
    for(unsigned uSce=0; uSce<NumScenario; uSce++){
      MC1_Particle tmp;

      int WhichPart = RanGen[ThId]->Integer(0,3);

//we need to replace the first element
//perhaps a small class to do that?
//we can have fixed amount of particles, and an iterator that saves the info
//about which is the newest and oldest member
      if(WhichPart==0){
        //for(int i=0; i<3; i++) tmp.xyz[i]=Part1_fsi[ThId][uSce][i];
        //if(MB1[ThId][uSce].size()<MixingDepth)
        //MB1[ThId][uSce].push_back(Part1_fsi[ThId][uSce]);
      }
      else if(WhichPart==1){
        //for(int i=0; i<3; i++) tmp.xyz[i]=Part2_fsi[ThId][uSce][i];
        //MB2[ThId][uSce].push_back(Part2_fsi[ThId][uSce]);
      }
      else if(WhichPart==2){
        //for(int i=0; i<3; i++) tmp.xyz[i]=Part3_fsi[ThId][uSce][i];
        //MB3[ThId][uSce].push_back(Part3_fsi[ThId][uSce]);
      }
      else printf("????!!!!!?????!!!!!\n");
    }




    //PartM1[ThId][0] = Part1[ThId][0];
    //PartM1[ThId][1] = Part1[ThId][1];
    //PartM1[ThId][2] = Part1[ThId][2];
    FirstIter = false;
  }//the parallel loop

//printf("Out of loop\n");
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

  DLM_Histo<float> C13_exclusive = (S13_k[MaxThreads][FSI_13]/S13_k[MaxThreads][noFSI]);
  DLM_Histo<float> C13_open = (S13_k[MaxThreads][FSI_2B]/S13_k[MaxThreads][noFSI]);
  DLM_Histo<float> C13_orgy = (S13_k[MaxThreads][FSI_ALL]/S13_k[MaxThreads][noFSI]);

  DLM_Histo<float> C23_exclusive = (S23_k[MaxThreads][FSI_23]/S23_k[MaxThreads][noFSI]);
  DLM_Histo<float> C23_open = (S23_k[MaxThreads][FSI_2B]/S23_k[MaxThreads][noFSI]);
  DLM_Histo<float> C23_orgy = (S23_k[MaxThreads][FSI_ALL]/S23_k[MaxThreads][noFSI]);


  DLM_Histo<float> C12_voyeur = (S_Q3[MaxThreads][FSI_12]/S_Q3[MaxThreads][noFSI]);
  DLM_Histo<float> C13_voyeur = (S_Q3[MaxThreads][FSI_13]/S_Q3[MaxThreads][noFSI]);
  DLM_Histo<float> C23_voyeur = (S_Q3[MaxThreads][FSI_23]/S_Q3[MaxThreads][noFSI]);


  //DLM_Histo<float> C12_mixed = (S_Q3[MaxThreads][FSIM_12]/S_Q3[MaxThreads][noFSI]);
  //DLM_Histo<float> C13_mixed = (S_Q3[MaxThreads][FSIM_13]/S_Q3[MaxThreads][noFSI]);
  //DLM_Histo<float> C23_mixed = (S_Q3[MaxThreads][FSIM_23]/S_Q3[MaxThreads][noFSI]);

  //this is the analytical computation
  //original
  DLM_Histo<float> C123_open_Kubo = (C12_voyeur+C13_voyeur+C23_voyeur)-2.;
  //DLM_Histo<float> C123_open_KuboM = (C12_mixed+C13_mixed+C23_mixed)-2.;


  //dimi style
  /*
  DLM_Histo<float> C123_open_Dimi =
    (S_Q3[MaxThreads][FSI_12]/S_Q3[MaxThreads][FSI_ALL])*C12_voyeur/3. +
    (S_Q3[MaxThreads][FSI_13]/S_Q3[MaxThreads][FSI_ALL])*C13_voyeur/3. +
    (S_Q3[MaxThreads][FSI_23]/S_Q3[MaxThreads][FSI_ALL])*C23_voyeur/3.;
  */
  //DLM_Histo<float> C123_open_Dimi = C12_voyeur/3. + C13_voyeur/3. + C23_voyeur/3.;
  DLM_Histo<float> C123_open_Dimi = C12_voyeur*C13_voyeur*C23_voyeur;
  //DLM_Histo<float> C123_open_DimiM = C12_mixed*C13_mixed*C23_mixed;
  DLM_Histo<float> R1R2R3_R = (S_Q3[MaxThreads][FSI_12]+S_Q3[MaxThreads][FSI_13]+S_Q3[MaxThreads][FSI_23])/S_Q3[MaxThreads][FSI_ALL];
  R1R2R3_R /= 3.;


  DLM_Histo<float> C123_open = (S_Q3[MaxThreads][FSI_2B]/S_Q3[MaxThreads][noFSI]);
  DLM_Histo<float> C123_threesome = (S_Q3[MaxThreads][FSI_3G]/S_Q3[MaxThreads][noFSI]);
  DLM_Histo<float> C123_Cum = (S_Q3[MaxThreads][FSI_3G]/S_Q3[MaxThreads][noFSI])-1.;
  DLM_Histo<float> C123_orgy = (S_Q3[MaxThreads][FSI_ALL]/S_Q3[MaxThreads][noFSI]);

  DLM_Histo<float> C123_Cum_Kubo = C123_orgy-C12_voyeur-C13_voyeur-C23_voyeur+2.;
  DLM_Histo<float> C123_Cum_Dimi = C123_orgy/(C12_voyeur*C13_voyeur*C23_voyeur)-1.;

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


//printf("Creating file\n");
  TFile fOutput(TString::Format("%s/MCUM/MC1.root",GetFemtoOutputFolder()),"recreate");
  TDirectory* Singles = fOutput.mkdir(TString::Format("Singles"));
  TDirectory* Pairs = fOutput.mkdir(TString::Format("Pairs"));
  TDirectory* Threesome = fOutput.mkdir(TString::Format("Threesome"));
  TDirectory* Ck = fOutput.mkdir(TString::Format("Ck"));
  TDirectory* CQ3 = fOutput.mkdir(TString::Format("CQ3"));
  TH1F** h_S1_p = new TH1F* [NumScenario];
  TH1F** h_S2_p = new TH1F* [NumScenario];
  TH1F** h_S3_p = new TH1F* [NumScenario];
  TH1F** h_S12_k = new TH1F* [NumScenario];
  TH1F** h_S13_k = new TH1F* [NumScenario];
  TH1F** h_S23_k = new TH1F* [NumScenario];
  TH1F** h_S_Q3 = new TH1F* [NumScenario];
  TH1F* h_C12_exclusive = Convert_DlmHisto_TH1F(&C12_exclusive,"h_C12_exclusive");
  TH1F* h_C12_open = Convert_DlmHisto_TH1F(&C12_open,"h_C12_open");
  TH1F* h_C12_orgy = Convert_DlmHisto_TH1F(&C12_orgy,"h_C12_orgy");
  TH1F* h_C13_exclusive = Convert_DlmHisto_TH1F(&C13_exclusive,"h_C13_exclusive");
  TH1F* h_C13_open = Convert_DlmHisto_TH1F(&C13_open,"h_C13_open");
  TH1F* h_C13_orgy = Convert_DlmHisto_TH1F(&C13_orgy,"h_C13_orgy");
  TH1F* h_C23_exclusive = Convert_DlmHisto_TH1F(&C23_exclusive,"h_C23_exclusive");
  TH1F* h_C23_open = Convert_DlmHisto_TH1F(&C23_open,"h_C23_open");
  TH1F* h_C23_orgy = Convert_DlmHisto_TH1F(&C23_orgy,"h_C23_orgy");
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
  TH1F* h_C123_open_Kubo = Convert_DlmHisto_TH1F(&C123_open_Kubo,"h_C123_open_Kubo");
  TH1F* h_C123_open_Dimi = Convert_DlmHisto_TH1F(&C123_open_Dimi,"h_C123_open_Dimi");
  //TH1F* h_C123_open_KuboM = Convert_DlmHisto_TH1F(&C123_open_KuboM,"h_C123_open_KuboM");
  //TH1F* h_C123_open_DimiM = Convert_DlmHisto_TH1F(&C123_open_DimiM,"h_C123_open_DimiM");
  TH1F* h_C123_Cum_Kubo = Convert_DlmHisto_TH1F(&C123_Cum_Kubo,"h_C123_Cum_Kubo");
  TH1F* h_C123_Cum_Dimi = Convert_DlmHisto_TH1F(&C123_Cum_Dimi,"h_C123_Cum_Dimi");
  TH1F* h_C123_Cum_Kubo3 = Convert_DlmHisto_TH1F(&C123_Cum_Kubo3,"h_C123_Cum_Kubo3");
  TH1F* h_C123_Cum_Dimi3 = Convert_DlmHisto_TH1F(&C123_Cum_Dimi3,"h_C123_Cum_Dimi3");
  TH1F* h_C123_Cum_Dimi_3Gx2B = Convert_DlmHisto_TH1F(&C123_Cum_Dimi_3Gx2B,"h_C123_Cum_Dimi_3Gx2B");
  TH1F* h_C2B_KuboTrue = Convert_DlmHisto_TH1F(&C2B_KuboTrue,"h_C2B_KuboTrue");
  TH1F* h_C2B_DimiTrue = Convert_DlmHisto_TH1F(&C2B_DimiTrue,"h_C2B_DimiTrue");
  TH1F* h_C2B_Kubo = Convert_DlmHisto_TH1F(&C2B_Kubo,"h_C2B_Kubo");
  TH1F* h_C2B_Dimi = Convert_DlmHisto_TH1F(&C2B_Dimi,"h_C2B_Dimi");
  TH1F* h_C2B_KuboOrgy = Convert_DlmHisto_TH1F(&C2B_KuboOrgy,"h_C2B_KuboOrgy");
  TH1F* h_C2B_DimiOrgy = Convert_DlmHisto_TH1F(&C2B_DimiOrgy,"h_C2B_DimiOrgy");


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


//printf("Creating histos\n");
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
  Ck->Add(h_C13_exclusive);
  Ck->Add(h_C13_open);
  Ck->Add(h_C13_orgy);
  Ck->Add(h_C23_exclusive);
  Ck->Add(h_C23_open);
  Ck->Add(h_C23_orgy);
  Ck->Write();

  CQ3->cd();
  CQ3->Add(h_C12_voyeur);
  CQ3->Add(h_C13_voyeur);
  CQ3->Add(h_C23_voyeur);
  //CQ3->Add(h_C12_mixed);
  //CQ3->Add(h_C13_mixed);
  //CQ3->Add(h_C23_mixed);
  CQ3->Add(h_C123_open_Kubo);
  CQ3->Add(h_C123_open_Dimi);
  //CQ3->Add(h_C123_open_KuboM);
  //CQ3->Add(h_C123_open_DimiM);
  CQ3->Add(h_R1R2R3_R);
  CQ3->Add(h_C123_open);
  CQ3->Add(h_C123_threesome);
  CQ3->Add(h_C123_Cum);
  CQ3->Add(h_C123_Cum_Kubo);
  CQ3->Add(h_C123_Cum_Dimi);
  CQ3->Add(h_C123_Cum_Kubo3);
  CQ3->Add(h_C123_Cum_Dimi3);
  CQ3->Add(h_C123_Cum_Dimi_3Gx2B);
  CQ3->Add(h_C2B_KuboTrue);
  CQ3->Add(h_C2B_DimiTrue);
  CQ3->Add(h_C2B_Kubo);
  CQ3->Add(h_C2B_Dimi);
  CQ3->Add(h_C2B_KuboOrgy);
  CQ3->Add(h_C2B_DimiOrgy);
  CQ3->Add(h_C123_orgy);
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
  delete h_C13_exclusive;
  delete h_C13_open;
  delete h_C13_orgy;
  delete h_C23_exclusive;
  delete h_C23_open;
  delete h_C23_orgy;
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
  delete h_C123_open_Kubo;
  delete h_C123_open_Dimi;
  //delete h_C123_open_KuboM;
  //delete h_C123_open_DimiM;
  delete h_R1R2R3_R;
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
  MC1(atoi(argv[1]));
  return 0;
}
