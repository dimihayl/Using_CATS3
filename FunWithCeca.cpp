
#include "FunWithCeca.h"

#include <iostream>

#include "TREPNI.h"

void Test1(){
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
}

int FUN_WITH_CECA(int argc, char *argv[]){

  Test1();

  return 0;
}
