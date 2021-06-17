
#ifndef PLAMBDA_1_H
#define PLAMBDA_1_H

class TString;
class DLM_CommonAnaFunctions;

void Fit_pL(DLM_CommonAnaFunctions& AnalysisObject, const TString& OutputFolder, const TString& OutputFileName,
            const TString& DataSample, const TString& SourceType, const double& SourceScale, const double& SourceStability,
            const TString& pL_Pot,const TString& FittingMode_pL, const int* VARIATIONS, const int& UniqueID, const TString ntFileName, const TString ntName);
void pL_SystematicsHM(const TString& OutputFolder, const int& WhichConfiguration, const int& FirstIter, const int& LastIter, const unsigned& OnlyFraction,
                      const int& RANDOMSEED, const int& Mode);


void pL_SystematicsMay2020(unsigned SEED, unsigned BASELINE_VAR, int POT_VAR, int Sigma0_Feed, int Data_Type,
                           bool DataSyst, bool FitSyst, bool Bootstrap, unsigned NumIter,
                           const char* CatsFileFolder, const char* OutputFolder);
void pLambda_Spline_Fit_Unfold2(const unsigned& SEEDmin, const unsigned& NumIter, const unsigned& TimeLimit,
                                const double& Perfect_chi2ndf, const double& VeryGood_chi2ndf, const double& Unacceptable_chi2ndf,
                                const double& BinWidth, const TString& DataVariation, const int& WhichPS,
                                const char* CatsFileFolder, const TString& OutputFolder);
void UpdateUnfoldFile(const char* CatsFileFolder, const TString& InputFolderName,
                      const TString& InputFileName, const TString& DataVariation,
                      const int& BinWidth);

int PLAMBDA_1_MAIN(int argc, char *argv[]);



#endif
