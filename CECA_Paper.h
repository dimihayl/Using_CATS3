#ifndef CECAPAPER_H
#define CECAPAPER_H

#include "TString.h"

int CECA_PAPER(int narg, char** ARGS);

int Ceca_pp_or_pL(const TString FileBase, const TString InputFolder, const TString OutputFolder, const TString LogFolder,
                  const int ParID, const int JobID, const int NumCPU=1);
int dadd_f(int argc, char *argv[]);
int dlmhst_root(int argc, char *argv[]);
int dlmhst_ceca_fit_wc0(int argc, char *argv[]);

#endif
