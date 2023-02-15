#ifndef CECAPAPER_H
#define CECAPAPER_H

#include "TString.h"

int CECA_PAPER(int narg, char** ARGS);

void Ceca_pp_or_pL(const TString InputFileBase, const TString OutputFileBase, const int SEED, const int NumCPU=1);
bool dadd_f(int argc, char *argv[]);

#endif
