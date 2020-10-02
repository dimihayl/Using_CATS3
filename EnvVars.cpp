#include "EnvVars.h"

#include <iostream>
#include <string>
#include <cstdlib>

using namespace std;

string  strCernBoxDimi = "CERNBOX_DIMI";
string  strFemtoOutputFolder = "FEMTO_OUTPUT";

char* GetCernBoxDimi(){
  char* loc;
  loc = getenv(strCernBoxDimi.c_str());
  return loc;
}

char* GetFemtoOutputFolder(){
  char* loc;
  loc = getenv(strFemtoOutputFolder.c_str());
  return loc;
}
