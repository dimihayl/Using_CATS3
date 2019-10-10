#ifndef KAONPROTON_H
#define KAONPROTON_H

class TString;

int KAONPROTON_MAIN(int argc, char *argv[]);
double* FindPotentialPars(const TString& InputFileName, const double& ScatLen, const double& EffRan, const unsigned WhichSol=0);
#endif
