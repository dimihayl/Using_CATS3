#ifndef PION_ANALYSIS_H
#define PION_ANALYSIS_H

#include <iostream>
#include <vector>

//
//class DLM_FemtoBootstrap{
//public:
//    DLM_Bootstrap();
//    ~DLM_Bootstrap();

//protected:

//    TH2F*

//};

class CATS;
class DLM_CleverMcLevyResoTM;

std::vector<double> CalculateLambdaParam(double pur1, double pur2, bool SamePart, bool extinput, bool verbose);
std::vector<double> CalculateLambdaParam2(double pur1, double pur2, bool SamePart, bool extinput, bool verbose);

void SetUpReso_pipi_CUSTOM(CATS& Kitty, DLM_CleverMcLevyResoTM& CleverMcLevyResoTM, const double& cTau, const double& avgmass, const double& fraction, const double& fractionomega, const double& fractionlong, const bool& randomangle=false);


int PION_ANA(int narg, char** ARGS);

#endif
