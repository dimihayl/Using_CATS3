#include "EnvVars.h"

#include <iostream>
#include <string>
#include <cstdlib>

using namespace std;

string  strCernBoxDimi = "CERNBOX_DIMI";
string  strFemtoOutputFolder = "FEMTO_OUTPUT";
string  strNx1Folder = "PATH_NX1";
string  strNx2Folder = "PATH_NX2";
string  strNx3Folder = "PATH_NX3";

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

char* GetNx1Folder(){
  char* loc;
  loc = getenv(strNx1Folder.c_str());
  return loc;
}
char* GetNx2Folder(){
  char* loc;
  loc = getenv(strNx2Folder.c_str());
  return loc;
}
char* GetNx3Folder(){
  char* loc;
  loc = getenv(strNx3Folder.c_str());
  return loc;
}
