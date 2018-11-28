
class TString;

void LamLamConfidence_22Oct2018(const TString Computer,
                                const TString OutputFolder, const int Chi2Map, const bool Bootstrap,
                                const int NumBootIter, const double TimeLimit, const int JOBID, const int SystID, const bool RandSyst, const int BootstrapID,
                                const int SepBL, const int FitFun,
                                const int BinWidthSyst, const int RadiusSyst, const int FitSyst, const TString DataSets);

void FinalExclusionPlot(const TString Computer,const int Chi2Map,const TString OutputFolder,const unsigned NumSyst);

int LamLamMcExclusionMAIN(int narg, char** ARGS);
