
#ifndef CATSPAPER_H
#define CATSPAPER_H

#include "CATS.h"
#include "TString.h"

class CATSparameters;

class CATSpaper{
public:
   CATSpaper(const TString& outdir);
   ~CATSpaper();

    //CATS* Kitty_pp;
    //CATS* Kitty_pL;
    CATS* FitKitty;

    void PhaseShiftPlot_pp();
    void PlotCS_pLambda();
    void PlotCk_pp();
    void PlotCk_pp_ver2();
    void BormioPlotCk_pp();
    void PlotCk_pL();
    void PlotCk_pL_ver2();
    void BormioPlotCk_pL();
    void PlotCk_EPOS();
    void PlotCk_LL();
    void PlotCk_LL_GIF();
    void PlotCk_LL_ver2();
    //void PlotCk_LL_ver2_Corrected1();
    void PlotCk_KminProton();
    void PlotCk_KplusProton();
    void PlotCk_XiMinusProton(const int WhichChannel, const int WhichAxis);
    void PlotCk_XiMinusProton_Corrected1();
    void CompareCatsCrab(const int& DataSet);
    void ExampleFit();
    void FitHadesData_pp();
    void FitHadesData_pL();

    void LamLamStudy1();

private:

    const TString OutputDir;

    const double Mass_p;
    const double Mass_L;

    CATSparameters* ppPars1S0;
    CATSparameters* ppPars3P0;
    CATSparameters* ppPars3P1;
    CATSparameters* ppPars3P2;
    double* ppSourcePars;
    unsigned NumMomBins_pp;
    double kMin_pp;
    double kMax_pp;

    CATSparameters* pLamPars1S0;
    CATSparameters* pLamPars3S1;
    double* pLamSourcePars;
    unsigned NumMomBins_pL;
    double kMin_pL;
    double kMax_pL;

    void StandardKitty_pp();
    void StandardKitty_pL();

    double Fitter1(double* mom, double* par);

    //pLab to pCm, Mass2 is the mass of the particle at rest
    double pLab_pCm(const double& pLab, const double& Mass1, const double& Mass2);
    //tLab to kCm, Mass2 is the mass of the particle at rest
    double tLab_pCm(const double& tLab, const double& Mass1, const double& Mass2);
    double pCm_pLab(const double& pCm, const double& Mass1, const double& Mass2);
    double pCm_tLab(const double& pCm, const double& Mass1, const double& Mass2);

};


#endif // CATSTOOLS_H
