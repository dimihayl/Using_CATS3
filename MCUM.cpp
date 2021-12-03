
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

void MC1(const int NumEvents){
  //cases are:
  //no final state interaction (FSI)
  //only FSI between particles 1 and 2
  //only FSI between 1 and 3
  //only FSI between 2 and 3
  //FSI between all 2-body pairs
  //no 2-body, only genuine 3 body
  //all 2-body and genuine 3 body
  const unsigned NumScenario = 7;
  enum Scenario { noFSI, FSI_12, FSI_13, FSI_23, FSI_2B, FSI_G3, FSI_ALL };
  //single particle momenta
  float** S1_p = new float* [NumScenario];
  float** S2_p = new float* [NumScenario];
  float** S3_p = new float* [NumScenario];

  float** S12_k = new float* [NumScenario];
  float** S13_k = new float* [NumScenario];
  float** S23_k = new float* [NumScenario];

  float** S12_Q3 = new float* [NumScenario];
  float** S13_Q3 = new float* [NumScenario];
  float** S23_Q3 = new float* [NumScenario];

  float** S123_Q3 = new float* [NumScenario];

  for(unsigned uSce=0; uSce<NumScenario; uSce++){
    S1_p[uSce] = new float [NumEvents];
    S2_p[uSce] = new float [NumEvents];
    S3_p[uSce] = new float [NumEvents];
    S12_k[uSce] = new float [NumEvents];
    S13_k[uSce] = new float [NumEvents];
    S23_k[uSce] = new float [NumEvents];
    S12_Q3[uSce] = new float [NumEvents];
    S13_Q3[uSce] = new float [NumEvents];
    S23_Q3[uSce] = new float [NumEvents];
    S123_Q3[uSce] = new float [NumEvents];
  }

  #pragma omp parallel for
  for(unsigned uEvent=0; uEvent<NumEvents; uEvent++){

  }


  for(unsigned uSce=0; uSce<NumScenario; uSce++){
    delete [] S1_p[uSce];
    delete [] S2_p[uSce];
    delete [] S3_p[uSce];
    delete [] S12_k[uSce];
    delete [] S13_k[uSce];
    delete [] S23_k[uSce];
    delete [] S12_Q3[uSce];
    delete [] S13_Q3[uSce];
    delete [] S23_Q3[uSce];
    delete [] S123_Q3[uSce];
  }
  delete [] S1_p;
  delete [] S2_p;
  delete [] S3_p;
  delete [] S12_k;
  delete [] S13_k;
  delete [] S23_k;
  delete [] S12_Q3;
  delete [] S13_Q3;
  delete [] S23_Q3;
  delete [] S123_Q3;
}

int MCUM(int argc, char *argv[]){
  NumThreads = std::thread::hardware_concurrency()==1?0:std::thread::hardware_concurrency()-1;
  printf("MCUM #%u\n",NumThreads);

  MC1(atoi(argv[1]));

  return 0;
}
