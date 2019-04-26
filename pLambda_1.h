
#ifndef PLAMBDA_1_H
#define PLAMBDA_1_H

class TString;
class DLM_CommonAnaFunctions;

void Fit_pL(DLM_CommonAnaFunctions& AnalysisObject, const TString& OutputFolder, const TString& OutputFileName,
            const TString& DataSample, const TString& SourceType, const double& SourceScale, const double& SourceStability,
            const TString& pL_Pot,const TString& FittingMode_pL, const int* VARIATIONS, const int& UniqueID, const TString ntFileName, const TString ntName);
void pL_SystematicsHM(const TString& OutputFolder, const int& WhichConfiguration, const int& FirstIter, const int& LastIter, const unsigned& OnlyFraction,
                      const int& RANDOMSEED, const int& Mode);


int PLAMBDA_1_MAIN(int argc, char *argv[]);



#endif

