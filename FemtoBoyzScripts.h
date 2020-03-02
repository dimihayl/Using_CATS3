
#ifndef FEMTOBOYZSCRIPTS_H
#define FEMTOBOYZSCRIPTS_H

#include "TString.h"
#include "TH1.h"

class TGraphErrors;


void SetStyle(bool graypalette=false, bool title=false);
void SetStyle2(bool graypalette=false, bool title=false);
void plotMorePreliminaries(const TString system);
void SetStyleHisto(TH1 *histo, int marker=-1, int color=-1);
void SetStyleHisto2(TH1 *histo, int marker, int color, double factor=1);
void PlotDimiExclusion_ver1(const TString FolderName, const TString HistoNameEP, const TString FileNameLED, const TString HistoNameLED,
                            const unsigned RebinX, const unsigned RebinY,
                            const TString DataSets);
void PlotDimiExclusionBE_ver1(const TString OutputFolder,
                              const TString FileNameStat, const TString FileNameSyst, const TString f0Inv_or_d0);
TGraphErrors *DrawSystematicError_FAST(TH1F* histexp,TH1F *histerr,TF1* ferr,double errorwidth);

void Plot_pL_FAST(
                        //const TString& WorkFolder,
                        //const TString& DataSample, const TString& SourceType, const TString& FittingMode_pp, const TString& FittingMode_pL,
                        //const TString& pL_Model1, const TString& pL_Model2,
                        //const TString& LegendSource_line1, const TString& LegendSource_line2,
                        //const double parA, const double parB, const TString DataSystFileName
                        );

int FEMTOBOYZ_MAIN(int narg, char** ARGS);

#endif // FEMTOBOYZSCRIPTS_H
