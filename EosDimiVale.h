#ifndef EOSDIMIVALE_H
#define EOSDIMIVALE_H

class CATS;
class TFile;

void CrossSectionFit_pL(CATS& Kitty, double& chi2, int& ndp, TFile* fOutput=NULL);
int EOSDIMIVALE(int argc, char *argv[]);

#endif
