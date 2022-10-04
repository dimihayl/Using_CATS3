#ifndef OTHERTASKS_H
#define OTHERTASKS_H

class CATS;
class TH1F;
class TF1;


void pp_CompareToNorfolk();
bool Eval_ScattParameters(CATS& Kitty, double& ScatLen, double& EffRan, TH1F*& hFit, TF1*& fitSP,
  const int& Nterms=2, const bool& Fixf0=false, const bool& Fixd0=false, const unsigned short& usCh=0);
int OTHERTASKS(int narg, char** ARGS);

#endif
