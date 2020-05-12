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

std::vector<double> CalculateLambdaParam(double pur1, double pur2, bool SamePart, bool extinput, bool verbose);
std::vector<double> CalculateLambdaParam2(double pur1, double pur2, bool SamePart, bool extinput, bool verbose);

int PION_ANA(int narg, char** ARGS);

#endif
