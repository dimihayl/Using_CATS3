
#include "ForBernie.h"

#include "CATStools.h"

#include "DLM_Source.h"
#include "DLM_Potentials.h"
#include "DLM_CkModels.h"
#include "CATS.h"
#include "DLM_CkDecomposition.h"
#include "DLM_Fitters.h"

#include "TRandom3.h"
#include "TH2F.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TGraph.h"
#include "TF1.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TMath.h"
#include "TCanvas.h"

#include "DLM_WfModel.h"

#include "math.h"

//#include "Ugly_mT_Resonance_Source1.h"

//so this one includes kT dependence and resonances, notice that everything is FIXED here, instead of giving huge arrays as input pars.
//[3] = linear parameter of the mT scaling
//[4] = slope parameter of the mT scaling
//[5] = which system we look at, 0,1,2,3 = pp pL pXim LL
//[6] = parameter flag 1 (resonances):
//      0 = the very first values from Vale (19th Sep 2018)
//      1 = the modified values from thermal model (21th Sep 2018)
//[7] = parameter flag 2 (mT binning):
//      0 = the very first values from Andi (19th Sep 2018), pp and pL are both fitted
void UglySource_MB13TeV_1(double* Pars){
    enum eSyst { pp, pL, pXim, LL };

    const unsigned NumSystems = 4;

    double ResMass_p;
    double ResMass_L;
    double ResMass_Xim;

    double ResTau_p;
    double ResTau_L;
    double ResTau_Xim;

    double Prim_p;
    double Prim_L;
    double Prim_Xim;

    const double Mass_p = 938.272;
    const double Mass_pi = 138; //kind of avg between pi+ and pi0
    const double Mass_L = 1115.683;
    const double Mass_Xim = 1321.7;

    if(round(Pars[6]==0)){

    }
    else if(round(Pars[6]==1)){

    }
    else{

    }

    unsigned* Num_mT = new unsigned [NumSystems];

    double** Mean_mT = new double* [NumSystems];
    double** Weight_mT = new double* [NumSystems];
    double** Radius_mT = new double* [NumSystems];
    double LinearTerm;
    double SlopeTerm;

    if(round(Pars[7]==0)){
        Num_mT[pp] = 4;
        Num_mT[pL] = 4;
        Num_mT[pXim] = 4;
        Num_mT[LL] = 4;

        for(unsigned uSyst=0; uSyst<NumSystems; uSyst++){
            Mean_mT[uSyst] = new double [Num_mT[uSyst]];
            Weight_mT[uSyst] = new double [Num_mT[uSyst]];
            Radius_mT[uSyst] = new double [Num_mT[uSyst]];
        }

        Mean_mT[0][0] = 1047;   Weight_mT[0][0] = 0.2902592257;
        Mean_mT[0][1] = 1165;   Weight_mT[0][1] = 0.3337219832;
        Mean_mT[0][2] = 1356;   Weight_mT[0][2] = 0.2568598596;
        Mean_mT[0][3] = 1850;   Weight_mT[0][3] = 0.1191589314;

        Mean_mT[1][0] = 1080;   Weight_mT[1][0] = 0.0031967451;
        Mean_mT[1][1] = 1168;   Weight_mT[1][1] = 0.2513975076;
        Mean_mT[1][2] = 1363;   Weight_mT[1][2] = 0.4929569038;
        Mean_mT[1][3] = 1911;   Weight_mT[1][3] = 0.2524488435;

        Mean_mT[2][0] = 1110;   Weight_mT[2][0] = 0;
        Mean_mT[2][1] = 1181;   Weight_mT[2][1] = 0.0252263126;
        Mean_mT[2][2] = 1369;   Weight_mT[2][2] = 0.5030778515;
        Mean_mT[2][3] = 1966;   Weight_mT[2][3] = 0.4716958358;

        Mean_mT[3][0] = 1110;   Weight_mT[3][0] = 0;
        Mean_mT[3][1] = 1193;   Weight_mT[3][1] = 0.0015342129;
        Mean_mT[3][2] = 1376;   Weight_mT[3][2] = 0.4470696533;
        Mean_mT[3][3] = 1991;   Weight_mT[3][3] = 0.5513961338;

        LinearTerm = 1.16320e+00;
        SlopeTerm = -2.52473e-04;
    }
    else{

    }

    //[system][paricle1/2]
    double** ResMass = new double* [NumSystems];
    double** ResTau = new double* [NumSystems];
    double** PrimFrac = new double* [NumSystems];
    double** Mass0 = new double* [NumSystems];
    double** Mass1 = new double* [NumSystems];
    for(unsigned uSyst=0; uSyst<NumSystems; uSyst++){
        ResMass[uSyst] = new double [2];
        ResTau[uSyst] = new double [2];
        PrimFrac[uSyst] = new double [2];
        Mass0[uSyst] = new double [2];
        Mass1[uSyst] = new double [2];
    }
    ResMass[pp][0] = ResMass_p;
    ResTau[pp][0] = ResTau_p;
    PrimFrac[pp][0] = Prim_p;
    Mass0[pp][0] = Mass_p;
    Mass1[pp][0] = Mass_pi;
    ResMass[pp][1] = ResMass_p;
    ResTau[pp][1] = ResTau_p;
    PrimFrac[pp][1] = Prim_p;
    Mass0[pp][1] = Mass_p;
    Mass1[pp][1] = Mass_pi;

    ResMass[pL][0] = ResMass_p;
    ResTau[pL][0] = ResTau_p;
    PrimFrac[pL][0] = Prim_p;
    Mass0[pL][0] = Mass_p;
    Mass1[pL][0] = Mass_pi;
    ResMass[pL][1] = ResMass_L;
    ResTau[pL][1] = ResTau_L;
    PrimFrac[pL][1] = Prim_L;
    Mass0[pL][1] = Mass_L;
    Mass1[pL][1] = Mass_pi;

    ResMass[pXim][0] = ResMass_p;
    ResTau[pXim][0] = ResTau_p;
    PrimFrac[pXim][0] = Prim_p;
    Mass0[pXim][0] = Mass_p;
    Mass1[pXim][0] = Mass_pi;
    ResMass[pXim][1] = ResMass_Xim;
    ResTau[pXim][1] = ResTau_Xim;
    PrimFrac[pXim][1] = Prim_Xim;
    Mass0[pXim][1] = Mass_Xim;
    Mass1[pXim][1] = Mass_pi;

    ResMass[LL][0] = ResMass_L;
    ResTau[LL][0] = ResTau_L;
    PrimFrac[LL][0] = Prim_L;
    Mass0[LL][0] = Mass_L;
    Mass1[LL][0] = Mass_pi;
    ResMass[LL][1] = ResMass_L;
    ResTau[LL][1] = ResTau_L;
    PrimFrac[LL][1] = Prim_L;
    Mass0[LL][1] = Mass_L;
    Mass1[LL][1] = Mass_pi;

    const unsigned NumRadBins = 256;
    const double rMin = 0;
    const double rMax = 8;
    //used for the fake random source
    const double RelError = 0.000;
    const double AbsError = 0.001;

    TF1* f_mT_R = new TF1("f_mT_R","[0]+x*[1]",rMin,rMax);
    f_mT_R->FixParameter(0,LinearTerm);
    f_mT_R->FixParameter(1,SlopeTerm);

    for(unsigned uSyst=0; uSyst<NumSystems; uSyst++){
        for(unsigned umT=0; umT<Num_mT[uSyst]; umT++){
            Radius_mT[uSyst][umT] = f_mT_R->Eval(Mean_mT[uSyst][umT]);
        }
    }



}

void RUN2_main(const unsigned& NumIter, const unsigned& NumJobs, const unsigned& JobID){

    TRandom3 rangen(1+JobID);

    //!SETTING THAT YOU MIGHT NEED
    //What source to use: 0 = Gauss; 1=Cauchy; 2=DoubleGauss; 3=EPOS; 4=Gauss+Cauchy; 5=MC_GAUSS_OutSideLong; 6=mT scaling + simple resonances
    //7=stable, 8=EPOS fitted with a double Gauss; 9=analytic levy; clever levy
    //int vSource = rangen.Integer(2); if(vSource==1) vSource=3; //chose Gauss/EPOS at random
    int vSource=10;

//(it worked with 1.1,1.75,1.55)
//(it worked with 1.0,1.60,1.40)
    const double Stable_Scale = 1.00;
    const double pp_Stability = 1.6;//
    const double pL_Stability = 1.4;//(for pL)
    const double pXim_Stability = 1.9;//

    const bool FAST_PLOT = true;
    const bool FULL_PLOT = false;
    //some bullshit that you might need (and ask me) when using EPOS
    const int EPOS_MAX_PAIRS = 64e6/8;
    const int EPOS_MIX_DEP = 8;
    const bool EPOS_THETA_DEP = false;
    const bool TAU_CORR = false;
    const double TRANSP_RENORM = 1.5;

    //pp13TeV,pPb5TeV
    const TString DataPeriod = "pp13TeV";
    //const TString DataPeriod = "pPb5TeV";
    //0 => the first analysis
    //1 => as we did it for the LL paper
    //2 => the new event mixing introduced (NOT FINISHED)
    const int DataFlag = 2;
    //0 initial
    //1 papers
    const int FitRegion = 1;

    //just some temp folder for temp files. Really not important, I will try to get rid of this soon
    const TString OutputDir = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/ForBernie/";

    //perform the systematics by using the different cut variations as further permutations, do NOT add the systematic errors in quadrature
    //TString SystematicsType = "CutVarIterate";
    //perform the systematics by adding the systematics errors of the bin quadratically to the data points
    TString SystematicsType = "CutVarAdd";
    const bool ExcludeXi = false; //(fit only pp/pL/LL just as in run1)



    double PurityProton;
    double PurityLambda;
    double PurityXi;
    double pp_f0;
    double pp_f1;
    double pL_f0;
    double pL_f1;
    double pL_f2;

    //(single particle quantities)
    if(DataPeriod=="pp13TeV"){
        //as used initially
        if(DataFlag==0){
            PurityProton = 0.991213;//991213
            PurityLambda = 0.965964;
            PurityXi = 0.956;

            pp_f0 = 0.874808;
            pp_f1 = 0.0876342;//fraction of Lambda

            pL_f0 = 0.619493;//fraction of primary Lambdas
            pL_f1 = 0.206498;//fraction of Sigma0
            pL_f2 = 0.0870044;//fractions of Xi0/m
        }
        else if(DataFlag==1||DataFlag==2){
            PurityProton = 0.989859;
            PurityLambda = 0.96768;
            PurityXi = 0.956;//is this correct?

            pp_f0 = 0.87397;
            pp_f1 = 0.0882211;//fraction of Lambda

            pL_f0 = 0.601008;//fraction of primary Lambdas
            pL_f1 = 0.200336;//fraction of Sigma0
            pL_f2 = 0.0993283;//fractions of Xi0/m
        }
    }
    else if(DataPeriod=="pPb5TeV"){
        if(DataFlag==0){
            PurityProton = 0.984177 ;//pPb 5 TeV
            PurityLambda = 0.93171;
            PurityXi = 0.9;//new cuts

            pp_f0 = 0.862431;
            pp_f1 = 0.0962984;

            pL_f0 = 0.744516;//fraction of primary Lambdas
            pL_f1 = 0.248172;//fraction of Sigma0
            pL_f2 = 0.003656;//fractions of Xi0/m



            PurityProton = 0.984265 ;//pPb 5 TeV
            PurityLambda = 0.937761;
            PurityXi = 0.88;//new cuts

            pp_f0 = 0.862814;
            pp_f1 = 0.09603;

            pL_f0 = 0.521433;
            pL_f1 = 0.173811;
            pL_f2 = 0.152378;
        }
        else if(DataFlag==1||DataFlag==2){
            //from 20th Aug 2018
            PurityProton = 0.984265 ;//pPb 5 TeV
            PurityLambda = 0.937761;
            PurityXi = 0.88;//new cuts

            pp_f0 = 0.862814;
            pp_f1 = 0.09603;

            pL_f0 = 0.521433;
            pL_f1 = 0.173811;
            pL_f2 = 0.152378;
        }
    }
    else{
        printf("FUCK YOU! Stirb du Wicht! Fat incular! Еби се в гъза!\n");
        return;
    }

//ADVANCED***
    double ProtonPrim = pp_f0;
    double arrayPercLamProton[3]={pp_f1/(1.-pp_f0)*0.8,pp_f1/(1.-pp_f0),pp_f1/(1.-pp_f0)*1.2};//+/- 20%

    const double SigLambdaPrimDir = pL_f0+pL_f1;
    double arrayPercSigLambda[3]={0.8*pL_f1/pL_f0,pL_f1/pL_f0,1.2*pL_f1/pL_f0};//1/3 +/- 20%
    double arrayPercXiLambda[3]={pL_f2/(1.-pL_f0-pL_f1)*0.8,pL_f2/(1.-pL_f0-pL_f1),pL_f2/(1.-pL_f0-pL_f1)*1.2};//+/- 20%

    //ratio Xi-(1530) to Xi-
    const double Xim1530_to_Xim = 0.32*(1./3.);
    //ratio Xi0(1530) to Xi0 (n=neutral)
    const double Xin1530_to_Xim = 0.32*(2./3.);
    const double Omegam_to_Xim = 0.1;
    const double OmegamXim_BR = 0.086;

    //following my lambda pars with the 3 possible modifications
    //for the proton:
    //0 = primary
    //1 = from Lambda
    //2 = other feeddown (flat)
    //3 = missidentified
    const unsigned NumChannels_p = 4;
    double** Purities_p = new double* [3];
    double** Fraction_p = new double* [3];
    for(unsigned uVar=0; uVar<3; uVar++){
        Purities_p[uVar] = new double [NumChannels_p];
        Fraction_p[uVar] = new double [NumChannels_p];

        Purities_p[uVar][0] = PurityProton;
        Purities_p[uVar][1] = PurityProton;
        Purities_p[uVar][2] = PurityProton;
        Purities_p[uVar][3] = 1.-PurityProton;

        Fraction_p[uVar][0] = ProtonPrim;
        Fraction_p[uVar][1] = (1.-ProtonPrim)*(arrayPercLamProton[uVar]);
        Fraction_p[uVar][2] = (1.-ProtonPrim)*(1.-arrayPercLamProton[uVar]);
        Fraction_p[uVar][3] = 1.;
    }

    //for the Lambda:
    //0 = primary
    //1 = from Sigma0
    //2 = from Xim
    //3 = from Xi0
    //4 = missidentified
    const unsigned NumChannels_L = 5;
    double*** Purities_L = new double** [3];
    double*** Fraction_L = new double** [3];
    for(unsigned uVarSL=0; uVarSL<3; uVarSL++){
        Purities_L[uVarSL] = new double* [3];
        Fraction_L[uVarSL] = new double* [3];
        for(unsigned uVarXi=0; uVarXi<3; uVarXi++){
            Purities_L[uVarSL][uVarXi] = new double [NumChannels_L];
            Fraction_L[uVarSL][uVarXi] = new double [NumChannels_L];

            Purities_L[uVarSL][uVarXi][0] = PurityLambda;
            Purities_L[uVarSL][uVarXi][1] = PurityLambda;
            Purities_L[uVarSL][uVarXi][2] = PurityLambda;
            Purities_L[uVarSL][uVarXi][3] = PurityLambda;
            Purities_L[uVarSL][uVarXi][4] = 1.-PurityLambda;

            //the array is r = S/L, and S+L=1 are the fractions of Sigmas and Lambdas
            double FracOfLambda = 1./(1.+arrayPercSigLambda[uVarSL]);
            Fraction_L[uVarSL][uVarXi][0] = SigLambdaPrimDir*FracOfLambda;
            Fraction_L[uVarSL][uVarXi][1] = SigLambdaPrimDir*(1.-FracOfLambda);
            Fraction_L[uVarSL][uVarXi][2] = (1.-SigLambdaPrimDir)*(arrayPercXiLambda[uVarXi]);
            Fraction_L[uVarSL][uVarXi][3] = (1.-SigLambdaPrimDir)*(1.-arrayPercXiLambda[uVarXi]);
            Fraction_L[uVarSL][uVarXi][4] = 1.;
        }
    }

    //for the Xi:
    //0 = primary
    //1 = from Xi-(1530)
    //2 = from Xi0(1530)
    //3 = from Omega
    //4 = missidentified
    const unsigned NumChannels_Xim = 5;
    double** Purities_Xim = new double* [3];
    double** Fraction_Xim = new double* [3];
    for(unsigned uVar=0; uVar<3; uVar++){
        Purities_Xim[uVar] = new double [NumChannels_Xim];
        Fraction_Xim[uVar] = new double [NumChannels_Xim];

        Purities_Xim[uVar][0] = PurityXi;
        Purities_Xim[uVar][1] = PurityXi;
        Purities_Xim[uVar][2] = PurityXi;
        Purities_Xim[uVar][3] = PurityXi;
        Purities_Xim[uVar][4] = 1.-PurityXi;

        //the ratios that we have for Xis are referred to the total number of Xi particles (which already include all contributions)
        //hence Xi1530_to_Xi indeed is simply the number of Xis that stem from a Xi1530
        Fraction_Xim[uVar][0] = 1.-Xim1530_to_Xim-Xin1530_to_Xim-Omegam_to_Xim*OmegamXim_BR;
        Fraction_Xim[uVar][1] = Xim1530_to_Xim;
        Fraction_Xim[uVar][2] = Xin1530_to_Xim;
        Fraction_Xim[uVar][3] = Omegam_to_Xim*OmegamXim_BR;
        Fraction_Xim[uVar][4] = 1.;
    }
//***

    //!Binning
    //This is for the CATS objects, make sure it covers the full femto range
    const unsigned NumMomBins_pp = 95;
    const double kMin_pp = 0;
    const double kMax_pp = kMin_pp+4*NumMomBins_pp;//(4 is the bin width)

    const unsigned NumMomBins_pL = 20;
    const double kMin_pL = 0;
    const double kMax_pL = kMin_pL+20*NumMomBins_pL;

    const unsigned NumMomBins_LL = 24;
    const double kMin_LL = 0;
    const double kMax_LL = kMin_LL+20*NumMomBins_LL;

    const unsigned NumMomBins_pXim = 20;
    const double kMin_pXim = 0;
    const double kMax_pXim = kMin_pXim+20*NumMomBins_pXim;

    //starting value, do not worry about it too much
    const double GaussSourceSize = 1.2;

//ADVANCED***
    //#,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
    double PotPars1S0[10]={0,0,NN_AV18,v18_Coupled3P2,1,1,1,0,0,0};
    double PotPars3P0[10]={0,0,NN_AV18,v18_Coupled3P2,1,1,1,1,1,0};
    double PotPars3P1[10]={0,0,NN_AV18,v18_Coupled3P2,1,1,1,1,1,1};
    double PotPars3P2[10]={0,0,NN_AV18,v18_Coupled3P2,1,1,1,1,1,2};

    CATSparameters cPotPars1S0(CATSparameters::tPotential,8,true);
    cPotPars1S0.SetParameters(&PotPars1S0[2]);
    CATSparameters cPotPars3P0(CATSparameters::tPotential,8,true);
    cPotPars3P0.SetParameters(&PotPars3P0[2]);
    CATSparameters cPotPars3P1(CATSparameters::tPotential,8,true);
    cPotPars3P1.SetParameters(&PotPars3P1[2]);
    CATSparameters cPotPars3P2(CATSparameters::tPotential,8,true);
    cPotPars3P2.SetParameters(&PotPars3P2[2]);

    const double Weight1S0 = 3./12.;
    const double Weight3P0 = 1./12.;
    const double Weight3P1 = 3./12.;
    const double Weight3P2 = 5./12.;

    MS_GaussExp_mT_Simple MSGEMTS_pp;
    MS_GaussExp_mT_Simple MSGEMTS_pL;
    MS_GaussExp_mT_Simple MSGEMTS_LL;
    MS_GaussExp_mT_Simple MSGEMTS_pXim;

    DLM_StableDistribution Stable_pp;
    DLM_StableDistribution Stable_pL;
    DLM_StableDistribution Stable_LL;
    DLM_StableDistribution Stable_pXim;

    DLM_CleverLevy CleverLevy_pp;
    DLM_CleverLevy CleverLevy_pL;
    DLM_CleverLevy CleverLevy_LL;
    DLM_CleverLevy CleverLevy_pXim;

    //0 is single Levy
    //1 is Levy for a difference
    //2 is the Nolan parametrization
    const int ClevLevType = 2;

    CleverLevy_pp.InitStability(20,1,2);
    CleverLevy_pp.InitScale(30,0.3,1.7);
    CleverLevy_pp.InitRad(256,0,64);
    CleverLevy_pp.InitType(ClevLevType);

    CleverLevy_pL.InitStability(20,1,2);
    CleverLevy_pL.InitScale(30,0.3,1.7);
    CleverLevy_pL.InitRad(256,0,64);
    CleverLevy_pL.InitType(ClevLevType);

    CleverLevy_LL.InitStability(20,1,2);
    CleverLevy_LL.InitScale(30,0.3,1.7);
    CleverLevy_LL.InitRad(256,0,64);
    CleverLevy_LL.InitType(ClevLevType);

    CleverLevy_pXim.InitStability(20,1,2);
    CleverLevy_pXim.InitScale(30,0.3,1.7);
    CleverLevy_pXim.InitRad(256,0,64);
    CleverLevy_pXim.InitType(ClevLevType);

    const double Mass_p=938.272; const double Mass_L=1115.683;
//const TString THERMALFILE = "/media/dmihaylov/My\ Book/TempCalc/thermalTest2_PRS=1.10_LRS=1.10_PTS=1.10_LTS=1.10.f19";
//const TString THERMALFILE = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/LambdaProtonPotentials/CRAB/CRAB_DLM2/ThermalSpectrum12/thermalTest12_PRS=1.60_1.60_1.10_PTS=1.00.f19";
//const TString THERMALFILE = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/LambdaProtonPotentials/CRAB/CRAB_DLM2/ThermalSpectrum_pXi_2/ThermalSpectrum_pXi_2_P=0.80_1.40_0.50_Xi=1.10_1.10_0.00_PTS=1.00.f19";

    CATS AB_pp;
    AB_pp.SetMaxNumThreads(1);
    double Pars_pp[9] = {0,0,0,GaussSourceSize*1.2,GaussSourceSize/1.2,0.5,1,1};
    CATSparameters Source_pp(CATSparameters::tSource,6,true);
    Source_pp.SetParameters(&Pars_pp[3]);
    //double Pars_pp[9] = {0,0,0,1.0,1.63,0.4253,1335.5,Mass_p,140};
    if(vSource==0){
        AB_pp.SetAnaSource(GaussSource, Source_pp);
        //AB_pp.SetAnaSource(GaussExpTotIdenticalSimple_2body, Pars_pp);
        //AB_pp.SetAnaSource((umtrs1_1.*(&UMTRS1::Source)),Pars_pp);
//(umtrs1_1.*Source)()
//&UMTRS1::Source
        AB_pp.SetUseAnalyticSource(true);
        AB_pp.SetThetaDependentSource(false);
    }
    else if(vSource==1){
        AB_pp.SetAnaSource(CauchySource, Source_pp);
        AB_pp.SetUseAnalyticSource(true);
        AB_pp.SetThetaDependentSource(false);
    }
    else if(vSource==2){
        AB_pp.SetAnaSource(DoubleGaussSource, Source_pp);
        AB_pp.SetUseAnalyticSource(true);
        AB_pp.SetThetaDependentSource(false);
    }
    else if(vSource==3){
        AB_pp.SetUseAnalyticSource(false);
        //AB_pp.SetInputFileName("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/Scratch9_OSCAR1997_AllLrzProtonFiles.f19");
        AB_pp.SetInputFileName("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/EPOS_LBF_pp91/pp_nist0_030219.f19");
        //AB_pp.SetInputFileName(THERMALFILE.Data());
        if(EPOS_MAX_PAIRS) AB_pp.SetMaxPairsToRead(EPOS_MAX_PAIRS);
        AB_pp.SetMaxPairsPerBin(16000/2);
        AB_pp.SetMixingDepth(EPOS_MIX_DEP);
        AB_pp.SetThetaDependentSource(EPOS_THETA_DEP);
        AB_pp.SetTransportRenorm(TRANSP_RENORM);
        if(EPOS_THETA_DEP){
            AB_pp.SetGridEpsilon(1./8192);
            AB_pp.SetGridMaxDepth(8);
        }
        AB_pp.SetTauCorrection(TAU_CORR);
    }
    else if(vSource==4){
        AB_pp.SetAnaSource(GaussCauchySource, Source_pp);
        AB_pp.SetUseAnalyticSource(true);
        AB_pp.SetThetaDependentSource(false);
    }
    else if(vSource==5){
        AB_pp.SetAnaSource(GaussOSL_MC, Source_pp);
        AB_pp.SetUseAnalyticSource(true);
        AB_pp.SetThetaDependentSource(false);
    }
    else if(vSource==6){
        MSGEMTS_pp.SetNum_mT(4);
        MSGEMTS_pp.SetMean_mT(0,1.047325*1000.);
        MSGEMTS_pp.SetMean_mT(1,1.16534*1000.);
        MSGEMTS_pp.SetMean_mT(2,1.35613*1000.);
        MSGEMTS_pp.SetMean_mT(3,1.849765*1000.);
        MSGEMTS_pp.SetWeight_mT(0,0.2902592257);
        MSGEMTS_pp.SetWeight_mT(1,0.3337219832);
        MSGEMTS_pp.SetWeight_mT(2,0.2568598596);
        MSGEMTS_pp.SetWeight_mT(3,0.1191589314);
        //Initial values
        //MSGEMTS_pp.SetLinear_mT(1.10576);
        //MSGEMTS_pp.SetSlope_mT(-0.253250e-3);
        //November revision
        MSGEMTS_pp.SetLinear_mT(1.11830e+00*1.1);
        MSGEMTS_pp.SetSlope_mT(-0.253753e-03*1.1);

    }
    else if(vSource==7){
        AB_pp.SetAnaSource(CatsSourceForwarder, &Stable_pp, 4);
        AB_pp.SetAnaSource(0,pp_Stability,false);
        AB_pp.SetAnaSource(1,0,false);
        AB_pp.SetAnaSource(2,Stable_Scale,false);
        AB_pp.SetAnaSource(3,0,false);
        AB_pp.SetUseAnalyticSource(true);
        AB_pp.SetThetaDependentSource(false);
    }
    else if(vSource==8){
        //use the mT scaling as a dummy way of including two Gaussians into the fit

        //function trough two points.
        //we set dummy values for the mT (1 and 2), while the values
        //for the radii are set based on our fit results
        double fun0=0.515*1.00;
        double fun1=1.105*1.00;
        double xVal0=1;
        double xVal1=2;
        double parA = (fun1-fun0)/(xVal1-xVal0);
        double parB = (fun0*xVal1-fun1*xVal0)/(xVal1-xVal0);

        MSGEMTS_pp.SetNum_mT(2);
        MSGEMTS_pp.SetMean_mT(0,1.);
        MSGEMTS_pp.SetMean_mT(1,2.);

        //these are the weights of the two radii
        MSGEMTS_pp.SetWeight_mT(0,0.563);
        MSGEMTS_pp.SetWeight_mT(1,0.437);

        MSGEMTS_pp.SetLinear_mT(parA);
        MSGEMTS_pp.SetSlope_mT(parB);
    }
    if(vSource==6||vSource==8){
        MSGEMTS_pp.SetMass(0,Mass_p);
        MSGEMTS_pp.SetMass(1,Mass_p);
        //Initial values
        //MSGEMTS_pp.SetMassR(0,1368.92);
        //MSGEMTS_pp.SetMassR(1,1368.92);
        //MSGEMTS_pp.SetMassD(0,138);
        //MSGEMTS_pp.SetMassD(1,138);
        //MSGEMTS_pp.SetTau(0,1.6034212583);
        //MSGEMTS_pp.SetTau(1,1.6034212583);
        //MSGEMTS_pp.SetResonanceWeight(0,1.-0.3261);
        //MSGEMTS_pp.SetResonanceWeight(1,1.-0.3261);
        //November revision
        MSGEMTS_pp.SetMassR(0,1361.52);
        MSGEMTS_pp.SetMassR(1,1361.52);
        MSGEMTS_pp.SetMassD(0,138);
        MSGEMTS_pp.SetMassD(1,138);
        MSGEMTS_pp.SetTau(0,1.65);
        MSGEMTS_pp.SetTau(1,1.65);
        MSGEMTS_pp.SetResonanceWeight(0,1.-0.3578);
        MSGEMTS_pp.SetResonanceWeight(1,1.-0.3578);
        AB_pp.SetAnaSource(CatsSourceForwarder, &MSGEMTS_pp);
        AB_pp.SetUseAnalyticSource(true);
        AB_pp.SetThetaDependentSource(false);
    }
    else if(vSource==9){
        Pars_pp[3] = GaussSourceSize; Pars_pp[4] = 2.0; Pars_pp[5] = 3;
        Source_pp.SetParameters(&Pars_pp[3]);
        AB_pp.SetAnaSource(LevySource3D_2particle, Source_pp);
        AB_pp.SetUseAnalyticSource(true);
        AB_pp.SetThetaDependentSource(false);
        AB_pp.SetMomentumDependentSource(false);
    }
    else if(vSource==10){
        AB_pp.SetAnaSource(CatsSourceForwarder, &CleverLevy_pp, 2);
        AB_pp.SetAnaSource(0,Stable_Scale,false);
        AB_pp.SetAnaSource(1,1,false);
        AB_pp.SetUseAnalyticSource(true);
        AB_pp.SetThetaDependentSource(false);
        AB_pp.SetMomentumDependentSource(false);
    }

    AB_pp.SetExcludeFailedBins(false);
    AB_pp.SetMomBins(NumMomBins_pp,kMin_pp,kMax_pp);

    AB_pp.SetNumChannels(4);
    AB_pp.SetNumPW(0,2);
    AB_pp.SetNumPW(1,2);
    AB_pp.SetNumPW(2,2);
    AB_pp.SetNumPW(3,2);
    AB_pp.SetSpin(0,0);
    AB_pp.SetSpin(1,1);
    AB_pp.SetSpin(2,1);
    AB_pp.SetSpin(3,1);
    AB_pp.SetChannelWeight(0, Weight1S0);
    AB_pp.SetChannelWeight(1, Weight3P0);
    AB_pp.SetChannelWeight(2, Weight3P1);
    AB_pp.SetChannelWeight(3, Weight3P2);

    AB_pp.SetQ1Q2(1);
    AB_pp.SetPdgId(2212, 2212);
    AB_pp.SetRedMass( 0.5*Mass_p );

    AB_pp.SetShortRangePotential(0,0,fDlmPot,cPotPars1S0);
    AB_pp.SetShortRangePotential(1,1,fDlmPot,cPotPars3P0);
    AB_pp.SetShortRangePotential(2,1,fDlmPot,cPotPars3P1);
    AB_pp.SetShortRangePotential(3,1,fDlmPot,cPotPars3P2);
//AB_pp.SetGamow(true);
//AB_pp.SetMaxNumThreads(1);

    AB_pp.KillTheCat();

    CATS AB_pL;
    AB_pL.SetMaxNumThreads(1);
    double Pars_pL[8] = {0,0,0,GaussSourceSize*1.2,GaussSourceSize/1.2,0.5,1,1};
    CATSparameters Source_pL(CATSparameters::tSource,6,true);
    Source_pL.SetParameters(&Pars_pL[3]);
    //double Pars_pL[14] = {0,0,0,1.0,1.63,0.4253,1335.5,Mass_p,140,5.49,0.3972,1338.6,Mass_L,140};
    if(vSource==0){
        AB_pL.SetAnaSource(GaussSource, Source_pL);
        //AB_pL.SetAnaSource(GaussExpTotSimple_2body, Pars_pL);
        AB_pL.SetUseAnalyticSource(true);
        AB_pL.SetThetaDependentSource(false);
    }
    else if(vSource==1){
        AB_pL.SetAnaSource(CauchySource, Source_pL);
        AB_pL.SetUseAnalyticSource(true);
        AB_pL.SetThetaDependentSource(false);
    }
    else if(vSource==2){
        AB_pL.SetAnaSource(DoubleGaussSource, Source_pL);
        AB_pL.SetUseAnalyticSource(true);
        AB_pL.SetThetaDependentSource(false);
    }
    else if(vSource==3){
        AB_pL.SetUseAnalyticSource(false);
        //AB_pL.SetInputFileName("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/Scratch9_OSCAR1997_AllLrzLambdaFiles.f19");
        AB_pL.SetInputFileName("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/EPOS_LBF_pp91/pL_nist0_030219.f19");
        //AB_pL.SetInputFileName(THERMALFILE.Data());
        if(EPOS_MAX_PAIRS) AB_pL.SetMaxPairsToRead(EPOS_MAX_PAIRS);
        AB_pL.SetMaxPairsPerBin(16000/2);
        AB_pL.SetMixingDepth(EPOS_MIX_DEP);
        AB_pL.SetThetaDependentSource(EPOS_THETA_DEP);
        AB_pL.SetTransportRenorm(TRANSP_RENORM);
        if(EPOS_THETA_DEP){
            AB_pL.SetGridEpsilon(1./8192);
            AB_pL.SetGridMaxDepth(8);
        }
        AB_pL.SetTauCorrection(TAU_CORR);
    }
    else if(vSource==4){
        AB_pL.SetAnaSource(GaussCauchySource, Source_pL);
        AB_pL.SetUseAnalyticSource(true);
        AB_pL.SetThetaDependentSource(false);
    }
    else if(vSource==5){
        AB_pL.SetAnaSource(GaussOSL_MC, Source_pL);
        AB_pL.SetUseAnalyticSource(true);
        AB_pL.SetThetaDependentSource(false);
    }
    else if(vSource==6){
        MSGEMTS_pL.SetNum_mT(4);
        MSGEMTS_pL.SetMean_mT(0,1.047325*1000.);
        MSGEMTS_pL.SetMean_mT(1,1.16534*1000.);
        MSGEMTS_pL.SetMean_mT(2,1.35613*1000.);
        MSGEMTS_pL.SetMean_mT(3,1.849765*1000.);
        MSGEMTS_pL.SetWeight_mT(0,0.0031967451);
        MSGEMTS_pL.SetWeight_mT(1,0.2513975076);
        MSGEMTS_pL.SetWeight_mT(2,0.4929569038);
        MSGEMTS_pL.SetWeight_mT(3,0.2524488435);
        //Initial values
        //MSGEMTS_pL.SetLinear_mT(1.10576);
        //MSGEMTS_pL.SetSlope_mT(-0.253250e-3);
        //November revision
        MSGEMTS_pL.SetLinear_mT(1.11830e+00*1.1);
        MSGEMTS_pL.SetSlope_mT(-0.253753e-03*1.1);
    }
    else if(vSource==7){
        AB_pL.SetAnaSource(CatsSourceForwarder, &Stable_pL, 4);
        AB_pL.SetAnaSource(0,pL_Stability,false);
        AB_pL.SetAnaSource(1,0,false);
        AB_pL.SetAnaSource(2,Stable_Scale,false);
        AB_pL.SetAnaSource(3,0,false);
        AB_pL.SetUseAnalyticSource(true);
        AB_pL.SetThetaDependentSource(false);
    }
    else if(vSource==8){
        //use the mT scaling as a dummy way of including two Gaussians into the fit

        //function trough two points.
        //we set dummy values for the mT (1 and 2), while the values
        //for the radii are set based on our fit results
        double fun0=0.515*1.00;
        double fun1=1.105*1.00;
        double xVal0=1;
        double xVal1=2;
        double parA = (fun1-fun0)/(xVal1-xVal0);
        double parB = (fun0*xVal1-fun1*xVal0)/(xVal1-xVal0);

        MSGEMTS_pL.SetNum_mT(2);
        MSGEMTS_pL.SetMean_mT(0,1.);
        MSGEMTS_pL.SetMean_mT(1,2.);

        //these are the weights of the two radii
        MSGEMTS_pL.SetWeight_mT(0,0.563);
        MSGEMTS_pL.SetWeight_mT(1,0.437);

        MSGEMTS_pL.SetLinear_mT(parA);
        MSGEMTS_pL.SetSlope_mT(parB);
    }
    if(vSource==6||vSource==8){
        //Initial values
        //MSGEMTS_pL.SetMass(0,Mass_p);
        //MSGEMTS_pL.SetMass(1,Mass_L);
        //MSGEMTS_pL.SetMassR(0,1368.92);
        //MSGEMTS_pL.SetMassR(1,1461.460557);
        //MSGEMTS_pL.SetMassD(0,138);
        //MSGEMTS_pL.SetMassD(1,138);
        //MSGEMTS_pL.SetTau(0,1.6034212583);
        //MSGEMTS_pL.SetTau(1,4.777240795);
        //MSGEMTS_pL.SetResonanceWeight(0,1.-0.3261);
        //MSGEMTS_pL.SetResonanceWeight(1,1.-0.3449);
        //November revision
        MSGEMTS_pL.SetMass(0,Mass_p);
        MSGEMTS_pL.SetMass(1,Mass_L);
        MSGEMTS_pL.SetMassR(0,1361.52);
        MSGEMTS_pL.SetMassR(1,1462.93);
        MSGEMTS_pL.SetMassD(0,138);
        MSGEMTS_pL.SetMassD(1,138);
        MSGEMTS_pL.SetTau(0,1.65);
        MSGEMTS_pL.SetTau(1,4.69);
        MSGEMTS_pL.SetResonanceWeight(0,1.-0.3578);
        MSGEMTS_pL.SetResonanceWeight(1,1.-0.3562);

        AB_pL.SetAnaSource(CatsSourceForwarder, &MSGEMTS_pL);
        AB_pL.SetUseAnalyticSource(true);
        AB_pL.SetThetaDependentSource(false);
    }
    else if(vSource==9){
        Pars_pL[3] = GaussSourceSize; Pars_pL[4] = 2.0; Pars_pL[5] = 3;
        Source_pL.SetParameters(&Pars_pL[3]);
        AB_pL.SetAnaSource(LevySource3D_2particle, Source_pL);
        AB_pL.SetUseAnalyticSource(true);
        AB_pL.SetThetaDependentSource(false);
        AB_pL.SetMomentumDependentSource(false);
    }
    else if(vSource==10){
        AB_pL.SetAnaSource(CatsSourceForwarder, &CleverLevy_pL, 2);
        AB_pL.SetAnaSource(0,Stable_Scale,false);
        AB_pL.SetAnaSource(1,1,false);
        AB_pL.SetUseAnalyticSource(true);
        AB_pL.SetThetaDependentSource(false);
        AB_pL.SetMomentumDependentSource(false);
    }
    AB_pL.SetMomBins(NumMomBins_pL,kMin_pL,kMax_pL);

    AB_pL.SetNumChannels(2);
    AB_pL.SetNumPW(0,1);
    AB_pL.SetNumPW(1,1);
    AB_pL.SetSpin(0,0);
    AB_pL.SetSpin(1,1);
    AB_pL.SetChannelWeight(0, 0.25);
    AB_pL.SetChannelWeight(1, 0.75);

    AB_pL.SetQ1Q2(0);
    AB_pL.SetPdgId(2212, 3122);
    AB_pL.SetRedMass( (Mass_p*Mass_L)/(Mass_p+Mass_L) );

    //!TEMPORARY: WE USE USMANI, ASK DIMI ABOUT FURTHER CHANGES (LIKE USE THE NLO)
    //#,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
    //double pLamPotPars1S0[10]={0,0,pL_UsmaniOli,0,0,0,0,0,0,0};
    //double pLamPotPars3S1[10]={0,0,pL_UsmaniOli,0,0,0,0,1,0,1};
    //AB_pL.SetShortRangePotential(0,0,fDlmPot,pLamPotPars1S0);
    //AB_pL.SetShortRangePotential(1,0,fDlmPot,pLamPotPars3S1);

    DLM_Histo<complex<double>>*** ExternalWF=NULL;
    ExternalWF = Init_pL_Haidenbauer(   "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/Haidenbauer/pLambdaNLO/",
                                AB_pL, 1, 600);
    AB_pL.SetExternalWaveFunction(ExternalWF);

//AB_pL.SetMaxNumThreads(1);
    AB_pL.KillTheCat();

    double pXimPotParsI0S0[11]={0,0,pXim_Lattice,12,0,-1,1,0,0,0,0};
    double pXimPotParsI0S1[11]={0,0,pXim_Lattice,12,0,-1,1,1,0,1,0};
    double pXimPotParsI1S0[11]={0,0,pXim_Lattice,6,1,1,1,0,0,0,0};
    double pXimPotParsI1S1[11]={0,0,pXim_Lattice,6,1,1,1,1,0,1,0};

    CATSparameters c_pXimPotParsI0S0(CATSparameters::tPotential,9,true);
    c_pXimPotParsI0S0.SetParameters(&pXimPotParsI0S0[2]);
    CATSparameters c_pXimPotParsI0S1(CATSparameters::tPotential,9,true);
    c_pXimPotParsI0S1.SetParameters(&pXimPotParsI0S1[2]);
    CATSparameters c_pXimPotParsI1S0(CATSparameters::tPotential,9,true);
    c_pXimPotParsI1S0.SetParameters(&pXimPotParsI1S0[2]);
    CATSparameters c_pXimPotParsI1S1(CATSparameters::tPotential,9,true);
    c_pXimPotParsI1S1.SetParameters(&pXimPotParsI1S1[2]);

    CATS AB_pXim;
    AB_pXim.SetMaxNumThreads(1);
    const double Mass_Xim = 1321.7;
    double Pars_pXi[8] = {0,0,0,GaussSourceSize*1.2,GaussSourceSize/1.2,0.5,1,1};
    CATSparameters Source_pXi(CATSparameters::tSource,6,true);
    Source_pXi.SetParameters(&Pars_pXi[3]);
    if(vSource==0){
        AB_pXim.SetAnaSource(GaussSource, Source_pXi);
        AB_pXim.SetUseAnalyticSource(true);
        AB_pXim.SetThetaDependentSource(false);
    }
    else if(vSource==1){
        AB_pXim.SetAnaSource(CauchySource, Source_pXi);
        AB_pXim.SetUseAnalyticSource(true);
        AB_pXim.SetThetaDependentSource(false);
    }
    else if(vSource==2){
        AB_pXim.SetAnaSource(DoubleGaussSource, Source_pXi);
        AB_pXim.SetUseAnalyticSource(true);
        AB_pXim.SetThetaDependentSource(false);
    }
    else if(vSource==3){
        AB_pXim.SetUseAnalyticSource(false);
        //AB_pXim.SetInputFileName("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/Scratch9_OSCAR1997_AllLrzLambdaFiles.f19");
        AB_pXim.SetInputFileName("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/EPOS_LBF_pp91/pL_nist0_030219.f19");
        //AB_pXim.SetInputFileName(THERMALFILE.Data());
        if(EPOS_MAX_PAIRS) AB_pXim.SetMaxPairsToRead(EPOS_MAX_PAIRS);
        AB_pXim.SetMaxPairsPerBin(16000/2);
        AB_pXim.SetMixingDepth(EPOS_MIX_DEP);
        AB_pXim.SetThetaDependentSource(EPOS_THETA_DEP);
        AB_pXim.SetTransportRenorm(TRANSP_RENORM);
        if(EPOS_THETA_DEP){
            AB_pXim.SetGridEpsilon(1./8192);
            AB_pXim.SetGridMaxDepth(8);
        }
        AB_pXim.SetTauCorrection(TAU_CORR);
    }
    else if(vSource==4){
        AB_pXim.SetAnaSource(GaussCauchySource, Source_pXi);
        AB_pXim.SetUseAnalyticSource(true);
        AB_pXim.SetThetaDependentSource(false);
    }
    else if(vSource==5){
        AB_pXim.SetAnaSource(GaussOSL_MC, Source_pXi);
        AB_pXim.SetUseAnalyticSource(true);
        AB_pXim.SetThetaDependentSource(false);
    }
    else if(vSource==6){
        MSGEMTS_pXim.SetNum_mT(4);
        MSGEMTS_pXim.SetMean_mT(0,1.047325*1000.);
        MSGEMTS_pXim.SetMean_mT(1,1.16534*1000.);
        MSGEMTS_pXim.SetMean_mT(2,1.35613*1000.);
        MSGEMTS_pXim.SetMean_mT(3,1.849765*1000.);
        MSGEMTS_pXim.SetWeight_mT(0,0.0);
        MSGEMTS_pXim.SetWeight_mT(1,0.0015342129);
        MSGEMTS_pXim.SetWeight_mT(2,0.4470696533);
        MSGEMTS_pXim.SetWeight_mT(3,0.5513961338);
        //initial values
        //MSGEMTS_pXim.SetLinear_mT(1.10576);
        //MSGEMTS_pXim.SetSlope_mT(-0.253250e-3);
        //November revision
        MSGEMTS_pXim.SetLinear_mT(1.11830e+00*1.1);
        MSGEMTS_pXim.SetSlope_mT(-0.253753e-03*1.1);
    }
    else if(vSource==7){
        AB_pXim.SetAnaSource(CatsSourceForwarder, &Stable_pXim, 4);
        AB_pXim.SetAnaSource(0,pXim_Stability,false);
        AB_pXim.SetAnaSource(1,0,false);
        AB_pXim.SetAnaSource(2,Stable_Scale,false);
        AB_pXim.SetAnaSource(3,0,false);
        AB_pXim.SetUseAnalyticSource(true);
        AB_pXim.SetThetaDependentSource(false);
    }
    else if(vSource==8){
        //use the mT scaling as a dummy way of including two Gaussians into the fit

        //function trough two points.
        //we set dummy values for the mT (1 and 2), while the values
        //for the radii are set based on our fit results
        double fun0=0.515*1.00;
        double fun1=1.105*1.00;
        double xVal0=1;
        double xVal1=2;
        double parA = (fun1-fun0)/(xVal1-xVal0);
        double parB = (fun0*xVal1-fun1*xVal0)/(xVal1-xVal0);

        MSGEMTS_pXim.SetNum_mT(2);
        MSGEMTS_pXim.SetMean_mT(0,1.);
        MSGEMTS_pXim.SetMean_mT(1,2.);

        //these are the weights of the two radii
        MSGEMTS_pXim.SetWeight_mT(0,0.563);
        MSGEMTS_pXim.SetWeight_mT(1,0.437);

        MSGEMTS_pXim.SetLinear_mT(parA);
        MSGEMTS_pXim.SetSlope_mT(parB);
    }
    if(vSource==6||vSource==8){
        //initial values
        //MSGEMTS_pXim.SetMass(0,Mass_p);
        //MSGEMTS_pXim.SetMass(1,Mass_Xim);
        //MSGEMTS_pXim.SetMassR(0,1368.92);
        //MSGEMTS_pXim.SetMassR(1,1823);
        //MSGEMTS_pXim.SetMassD(0,138);
        //MSGEMTS_pXim.SetMassD(1,138);
        //MSGEMTS_pXim.SetTau(0,1.6034212583);
        //MSGEMTS_pXim.SetTau(1,8.220833333);
        //MSGEMTS_pXim.SetResonanceWeight(0,1.-0.3261);
        //MSGEMTS_pXim.SetResonanceWeight(1,1.-0.9545);
        //November revision
        MSGEMTS_pXim.SetMass(0,Mass_p);
        MSGEMTS_pXim.SetMass(1,Mass_Xim);
        MSGEMTS_pXim.SetMassR(0,1361.52);
        MSGEMTS_pXim.SetMassR(1,1823);
        MSGEMTS_pXim.SetMassD(0,138);
        MSGEMTS_pXim.SetMassD(1,138);
        MSGEMTS_pXim.SetTau(0,1.65);
        MSGEMTS_pXim.SetTau(1,8.220833333);
        MSGEMTS_pXim.SetResonanceWeight(0,1.-0.3261);
        MSGEMTS_pXim.SetResonanceWeight(1,1.-0.9545);

        AB_pXim.SetAnaSource(CatsSourceForwarder, &MSGEMTS_pXim);
        AB_pXim.SetUseAnalyticSource(true);
        AB_pXim.SetThetaDependentSource(false);
    }
    else if(vSource==9){
        Pars_pXi[3] = GaussSourceSize; Pars_pXi[4] = 2.0; Pars_pXi[5] = 3;
        Source_pXi.SetParameters(&Pars_pXi[3]);
        AB_pXim.SetAnaSource(LevySource3D_2particle, Source_pXi);
        AB_pXim.SetUseAnalyticSource(true);
        AB_pXim.SetThetaDependentSource(false);
        AB_pXim.SetMomentumDependentSource(false);
    }
    else if(vSource==10){
        AB_pXim.SetAnaSource(CatsSourceForwarder, &CleverLevy_pXim, 2);
        AB_pXim.SetAnaSource(0,Stable_Scale,false);
        AB_pXim.SetAnaSource(1,1,false);
        AB_pXim.SetUseAnalyticSource(true);
        AB_pXim.SetThetaDependentSource(false);
        AB_pXim.SetMomentumDependentSource(false);
    }
    AB_pXim.SetExcludeFailedBins(false);
    AB_pXim.SetMomBins(NumMomBins_pXim,kMin_pXim,kMax_pXim);

    AB_pXim.SetNumChannels(4);
    AB_pXim.SetNumPW(0,1);
    AB_pXim.SetNumPW(1,1);
    AB_pXim.SetNumPW(2,1);
    AB_pXim.SetNumPW(3,1);
    AB_pXim.SetSpin(0,0);//I=0; S=0
    AB_pXim.SetSpin(1,1);//I=0; S=1
    AB_pXim.SetSpin(2,0);//I=1; S=0
    AB_pXim.SetSpin(3,1);//I=1; S=1
    AB_pXim.SetChannelWeight(0, 1./8.);
    AB_pXim.SetChannelWeight(1, 3./8.);
    AB_pXim.SetChannelWeight(2, 1./8.);
    AB_pXim.SetChannelWeight(3, 3./8.);

    AB_pXim.SetQ1Q2(-1);
    //AB_pXim.SetPdgId(2212, 3312);

    AB_pXim.SetPdgId(2212, 3122);//same as Lambda, in case we want to use EPOS pL source

    //const double Mass_p = 938.272;
    AB_pXim.SetRedMass( (Mass_p*Mass_Xim)/(Mass_p+Mass_Xim) );

    AB_pXim.SetShortRangePotential(0,0,fDlmPot,c_pXimPotParsI0S0);
    AB_pXim.SetShortRangePotential(1,0,fDlmPot,c_pXimPotParsI0S1);
    AB_pXim.SetShortRangePotential(2,0,fDlmPot,c_pXimPotParsI1S0);
    AB_pXim.SetShortRangePotential(3,0,fDlmPot,c_pXimPotParsI1S1);

    AB_pXim.SetMaxRad(64);
    AB_pXim.SetMaxRho(32);
//AB_pXim.SetMaxNumThreads(1);
    AB_pXim.KillTheCat();

    CATS AB_pXim1530;
    AB_pXim1530.SetMaxNumThreads(1);
    double Pars_pXim1530[8] = {0,0,0,GaussSourceSize*1.2,GaussSourceSize/1.2,0.5,1,1};
    CATSparameters Source_pXi1530(CATSparameters::tSource,6,true);
    Source_pXi1530.SetParameters(&Pars_pXim1530[3]);
    if(vSource==0){
        AB_pXim1530.SetAnaSource(GaussSource, Source_pXi1530);
        AB_pXim1530.SetUseAnalyticSource(true);
        AB_pXim1530.SetThetaDependentSource(false);
    }
    else if(vSource==1){
        AB_pXim1530.SetAnaSource(CauchySource, Source_pXi1530);
        AB_pXim1530.SetUseAnalyticSource(true);
        AB_pXim1530.SetThetaDependentSource(false);
    }
    else if(vSource==2){
        AB_pXim1530.SetAnaSource(DoubleGaussSource, Source_pXi1530);
        AB_pXim1530.SetUseAnalyticSource(true);
        AB_pXim1530.SetThetaDependentSource(false);
    }
    else if(vSource==3){
        AB_pXim1530.SetUseAnalyticSource(false);
        //AB_pXim1530.SetInputFileName("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/Scratch9_OSCAR1997_AllLrzLambdaFiles.f19");
        AB_pXim1530.SetInputFileName("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/EPOS_LBF_pp91/pL_nist0_030219.f19");
        //AB_pXim1530.SetInputFileName(THERMALFILE.Data());
        if(EPOS_MAX_PAIRS) AB_pXim1530.SetMaxPairsToRead(EPOS_MAX_PAIRS);
        AB_pXim1530.SetMaxPairsPerBin(16000/2);
        AB_pXim1530.SetMixingDepth(EPOS_MIX_DEP);
        AB_pXim1530.SetThetaDependentSource(EPOS_THETA_DEP);
        AB_pXim1530.SetTransportRenorm(TRANSP_RENORM);
        if(EPOS_THETA_DEP){
            AB_pXim1530.SetGridEpsilon(1./8192);
            AB_pXim1530.SetGridMaxDepth(8);
        }
        AB_pXim1530.SetTauCorrection(TAU_CORR);
    }
    else if(vSource==4){
        AB_pXim1530.SetAnaSource(GaussCauchySource, Source_pXi1530);
        AB_pXim1530.SetUseAnalyticSource(true);
        AB_pXim1530.SetThetaDependentSource(false);
    }
    else if(vSource==5){
        AB_pXim1530.SetAnaSource(GaussOSL_MC, Source_pXi1530);
        AB_pXim1530.SetUseAnalyticSource(true);
        AB_pXim1530.SetThetaDependentSource(false);
    }
    else if(vSource==6){
        AB_pXim1530.SetAnaSource(GaussSource, Source_pXi1530);
        AB_pXim1530.SetAnaSource(0,1.0);
        AB_pXim1530.SetUseAnalyticSource(true);
        AB_pXim1530.SetThetaDependentSource(false);
    }
    else if(vSource==7){
        AB_pXim1530.SetAnaSource(CatsSourceForwarder, &Stable_pXim, 4);
        AB_pXim1530.SetAnaSource(0,pp_Stability,false);
        AB_pXim1530.SetAnaSource(1,0,false);
        AB_pXim1530.SetAnaSource(2,Stable_Scale,false);
        AB_pXim1530.SetAnaSource(3,0,false);
        AB_pXim1530.SetUseAnalyticSource(true);
        AB_pXim1530.SetThetaDependentSource(false);
    }
    else if(vSource==8){
        AB_pXim1530.SetAnaSource(GaussSource, Source_pXi1530);
        AB_pXim1530.SetAnaSource(0,1.0);
        AB_pXim1530.SetUseAnalyticSource(true);
        AB_pXim1530.SetThetaDependentSource(false);
    }
    else if(vSource==9){
        Pars_pXim1530[3] = GaussSourceSize; Pars_pXim1530[4] = 2.0; Pars_pXim1530[5] = 3;
        Source_pXi1530.SetParameters(&Pars_pXim1530[3]);
        AB_pXim1530.SetAnaSource(LevySource3D_2particle, Source_pXi1530);
        AB_pXim1530.SetUseAnalyticSource(true);
        AB_pXim1530.SetThetaDependentSource(false);
        AB_pXim1530.SetMomentumDependentSource(false);
    }
    else if(vSource==10){
        AB_pXim1530.SetAnaSource(CatsSourceForwarder, &CleverLevy_pXim, 2);
        AB_pXim1530.SetAnaSource(0,Stable_Scale,false);
        AB_pXim1530.SetAnaSource(1,2,false);
        AB_pXim1530.SetUseAnalyticSource(true);
        AB_pXim1530.SetThetaDependentSource(false);
        AB_pXim1530.SetMomentumDependentSource(false);
    }
    AB_pXim1530.SetExcludeFailedBins(false);
    AB_pXim1530.SetMomBins(NumMomBins_pXim,kMin_pXim,kMax_pXim);

    AB_pXim1530.SetNumChannels(1);
    AB_pXim1530.SetNumPW(0,1);
    AB_pXim1530.SetSpin(0,0);
    AB_pXim1530.SetChannelWeight(0, 1.);

    AB_pXim1530.SetQ1Q2(-1);
    //AB_pXim1530.SetPdgId(2212, 3312);
    AB_pXim1530.SetPdgId(2212, 3122);

    const double Mass_Xim1530 = 1535;
    AB_pXim1530.SetRedMass( (Mass_p*Mass_Xim1530)/(Mass_p+Mass_Xim1530) );
//AB_pXim1530.SetMaxNumThreads(1);
    AB_pXim1530.KillTheCat();

//***

    //location of the file containing the feed-down smearing matrices
    const TString ResMatrixFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/DecayMatrices/run2_decay_matrices_old.root";
    //location of the file containing the momentum smearing matrices
    const TString SigmaMatrixFileName =
        DataPeriod=="pp13TeV"?  "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/ResolutionMatrices/Sample6_MeV_compact.root":
                                "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pPb_5TeV/ResolutionMatrices/Sample3_MeV_compact.root";

//ADVANCED***
    //!DO NOT TOUCH UNLESS YOU CHANGE THE ABOVE FILES, IN WHICH CASE ASK DIMI FOR HELP
    //1/FractionOfBins th number of original bins of the correction matrix are taken into account
    //originally: 1000 MeV => 1/2 should be good most of the time
    const int Fraction_Res = 2;
    const int Fraction_Sig = 1;
    const double UnitConv_Res = 1;
    const double UnitConv_Sig = 1;

    TH2F* hRes_pp_pL;
    TH2F* hRes_pL_pSigma0;
    TH2F* hRes_pL_pXim;
    TH2F* hRes_pXim_pXim1530;
    TH2F* hSigma_pp;
    TH2F* hSigma_pL;
    TH2F* hSigma_LL;
    TH2F* hSigma_pXim;

    TFile* FileRes = new TFile(ResMatrixFileName, "read");
    TFile* FileSigma = new TFile(SigmaMatrixFileName, "read");

    //for the Xi we will make the assumption that all residuals are flat since we do not know better
    FileRes->cd();
    hRes_pp_pL = (TH2F*)FileRes->Get("hRes_pp_pL");
    hRes_pL_pSigma0 = (TH2F*)FileRes->Get("hRes_pL_pSigma0");
    hRes_pL_pXim = (TH2F*)FileRes->Get("hRes_pL_pXim");
    hRes_pXim_pXim1530 = ExcludeXi?NULL:(TH2F*)FileRes->Get("hRes_pXim_pXim1530");

    FileSigma->cd();

    hSigma_pp = (TH2F*)FileSigma->Get("hSigmaMeV_Proton_Proton");
    hSigma_pL = (TH2F*)FileSigma->Get("hSigmaMeV_Proton_Lambda");
    hSigma_LL = (TH2F*)FileSigma->Get("hSigmaMeV_Lambda_Lambda");
    hSigma_pXim = (TH2F*)FileSigma->Get("hSigmaMeV_Proton_Xim");

    TH2F* hRes_pp_pL_MeV = new TH2F("hRes_pp_pL_MeV", "hRes_pp_pL_MeV",
        hRes_pp_pL->GetNbinsX()/Fraction_Res, hRes_pp_pL->GetXaxis()->GetBinLowEdge(1)*UnitConv_Res,hRes_pp_pL->GetXaxis()->GetBinUpEdge(hRes_pp_pL->GetNbinsX()/Fraction_Res)*UnitConv_Res,
        hRes_pp_pL->GetNbinsY()/Fraction_Res, hRes_pp_pL->GetYaxis()->GetBinLowEdge(1)*UnitConv_Res,hRes_pp_pL->GetXaxis()->GetBinUpEdge(hRes_pp_pL->GetNbinsY()/Fraction_Res)*UnitConv_Res);
    TH2F* hRes_pL_pSigma0_MeV = new TH2F("hRes_pL_pSigma0_MeV", "hRes_pL_pSigma0_MeV",
        hRes_pL_pSigma0->GetNbinsX()/Fraction_Res, hRes_pL_pSigma0->GetXaxis()->GetBinLowEdge(1)*UnitConv_Res,hRes_pL_pSigma0->GetXaxis()->GetBinUpEdge(hRes_pL_pSigma0->GetNbinsX()/Fraction_Res)*UnitConv_Res,
        hRes_pL_pSigma0->GetNbinsY()/Fraction_Res, hRes_pL_pSigma0->GetYaxis()->GetBinLowEdge(1)*UnitConv_Res,hRes_pL_pSigma0->GetXaxis()->GetBinUpEdge(hRes_pL_pSigma0->GetNbinsY()/Fraction_Res)*UnitConv_Res);
    TH2F* hRes_pL_pXim_MeV = new TH2F("hRes_pL_pXim_MeV", "hRes_pL_pXim_MeV",
        hRes_pL_pXim->GetNbinsX()/Fraction_Res, hRes_pL_pXim->GetXaxis()->GetBinLowEdge(1)*UnitConv_Res,hRes_pL_pXim->GetXaxis()->GetBinUpEdge(hRes_pL_pXim->GetNbinsX()/Fraction_Res)*UnitConv_Res,
        hRes_pL_pXim->GetNbinsY()/Fraction_Res, hRes_pL_pXim->GetYaxis()->GetBinLowEdge(1)*UnitConv_Res,hRes_pL_pXim->GetXaxis()->GetBinUpEdge(hRes_pL_pXim->GetNbinsY()/Fraction_Res)*UnitConv_Res);
    TH2F* hRes_pXim_pXim1530_MeV = hRes_pXim_pXim1530==NULL?NULL:new TH2F("hRes_pXim_pXim1530_MeV", "hRes_pXim_pXim1530_MeV",
        hRes_pXim_pXim1530->GetNbinsX()/Fraction_Res, hRes_pXim_pXim1530->GetXaxis()->GetBinLowEdge(1)*UnitConv_Res,hRes_pXim_pXim1530->GetXaxis()->GetBinUpEdge(hRes_pXim_pXim1530->GetNbinsX()/Fraction_Res)*UnitConv_Res,
        hRes_pXim_pXim1530->GetNbinsY()/Fraction_Res, hRes_pXim_pXim1530->GetYaxis()->GetBinLowEdge(1)*UnitConv_Res,hRes_pXim_pXim1530->GetXaxis()->GetBinUpEdge(hRes_pXim_pXim1530->GetNbinsY()/Fraction_Res)*UnitConv_Res);

    TH2F* hSigma_pp_MeV = new TH2F("hSigma_pp_MeV", "hSigma_pp_MeV",
        hSigma_pp->GetNbinsX()/Fraction_Sig, hSigma_pp->GetXaxis()->GetBinLowEdge(1)*UnitConv_Sig,hSigma_pp->GetXaxis()->GetBinUpEdge(hSigma_pp->GetNbinsX()/Fraction_Sig)*UnitConv_Sig,
        hSigma_pp->GetNbinsY()/Fraction_Sig, hSigma_pp->GetYaxis()->GetBinLowEdge(1)*UnitConv_Sig,hSigma_pp->GetXaxis()->GetBinUpEdge(hSigma_pp->GetNbinsY()/Fraction_Sig)*UnitConv_Sig);
    TH2F* hSigma_pL_MeV = new TH2F("hSigma_pL_MeV", "hSigma_pL_MeV",
        hSigma_pL->GetNbinsX()/Fraction_Sig, hSigma_pL->GetXaxis()->GetBinLowEdge(1)*UnitConv_Sig,hSigma_pL->GetXaxis()->GetBinUpEdge(hSigma_pL->GetNbinsX()/Fraction_Sig)*UnitConv_Sig,
        hSigma_pL->GetNbinsY()/Fraction_Sig, hSigma_pL->GetYaxis()->GetBinLowEdge(1)*UnitConv_Sig,hSigma_pL->GetXaxis()->GetBinUpEdge(hSigma_pL->GetNbinsY()/Fraction_Sig)*UnitConv_Sig);
    TH2F* hSigma_LL_MeV = hSigma_LL==NULL?NULL:new TH2F("hSigma_LL_MeV", "hSigma_LL_MeV",
        hSigma_LL->GetNbinsX()/Fraction_Sig, hSigma_LL->GetXaxis()->GetBinLowEdge(1)*UnitConv_Sig,hSigma_LL->GetXaxis()->GetBinUpEdge(hSigma_LL->GetNbinsX()/Fraction_Sig)*UnitConv_Sig,
        hSigma_LL->GetNbinsY()/Fraction_Sig, hSigma_LL->GetYaxis()->GetBinLowEdge(1)*UnitConv_Sig,hSigma_LL->GetXaxis()->GetBinUpEdge(hSigma_LL->GetNbinsY()/Fraction_Sig)*UnitConv_Sig);
    TH2F* hSigma_pXim_MeV = hSigma_pXim==NULL?NULL:new TH2F("hSigma_pXi_MeV", "hSigma_pXi_MeV",
        hSigma_pXim->GetNbinsX()/Fraction_Sig, hSigma_pXim->GetXaxis()->GetBinLowEdge(1)*UnitConv_Sig,hSigma_pXim->GetXaxis()->GetBinUpEdge(hSigma_pXim->GetNbinsX()/Fraction_Sig)*UnitConv_Sig,
        hSigma_pXim->GetNbinsY()/Fraction_Sig, hSigma_pXim->GetYaxis()->GetBinLowEdge(1)*UnitConv_Sig,hSigma_pXim->GetXaxis()->GetBinUpEdge(hSigma_pXim->GetNbinsY()/Fraction_Sig)*UnitConv_Sig);

    for(int iBinX=1; iBinX<=hRes_pp_pL->GetNbinsX()/Fraction_Res; iBinX++){
        for(int iBinY=1; iBinY<=hRes_pp_pL->GetNbinsY()/Fraction_Res; iBinY++){
            hRes_pp_pL_MeV->SetBinContent(iBinX, iBinY, hRes_pp_pL->GetBinContent(iBinX, iBinY));
        }
    }
    for(int iBinX=1; iBinX<=hRes_pL_pSigma0->GetNbinsX()/Fraction_Res; iBinX++){
        for(int iBinY=1; iBinY<=hRes_pL_pSigma0->GetNbinsY()/Fraction_Res; iBinY++){
            hRes_pL_pSigma0_MeV->SetBinContent(iBinX, iBinY, hRes_pL_pSigma0->GetBinContent(iBinX, iBinY));
        }
    }
    for(int iBinX=1; iBinX<=hRes_pL_pXim->GetNbinsX()/Fraction_Res; iBinX++){
        for(int iBinY=1; iBinY<=hRes_pL_pXim->GetNbinsY()/Fraction_Res; iBinY++){
            hRes_pL_pXim_MeV->SetBinContent(iBinX, iBinY, hRes_pL_pXim->GetBinContent(iBinX, iBinY));
        }
    }
    if(hRes_pXim_pXim1530 && hRes_pXim_pXim1530_MeV){
        for(int iBinX=1; iBinX<=hRes_pXim_pXim1530->GetNbinsX()/Fraction_Res; iBinX++){
            for(int iBinY=1; iBinY<=hRes_pXim_pXim1530->GetNbinsY()/Fraction_Res; iBinY++){
                hRes_pXim_pXim1530_MeV->SetBinContent(iBinX, iBinY, hRes_pXim_pXim1530->GetBinContent(iBinX, iBinY));
            }
        }
    }
    for(int iBinX=1; iBinX<=hSigma_pp->GetNbinsX()/Fraction_Sig; iBinX++){
        for(int iBinY=1; iBinY<=hSigma_pp->GetNbinsY()/Fraction_Sig; iBinY++){
            hSigma_pp_MeV->SetBinContent(iBinX, iBinY, hSigma_pp->GetBinContent(iBinX, iBinY));
        }
    }
    for(int iBinX=1; iBinX<=hSigma_pL->GetNbinsX()/Fraction_Sig; iBinX++){
        for(int iBinY=1; iBinY<=hSigma_pL->GetNbinsY()/Fraction_Sig; iBinY++){
            hSigma_pL_MeV->SetBinContent(iBinX, iBinY, hSigma_pL->GetBinContent(iBinX, iBinY));
        }
    }
    if(hSigma_LL && hSigma_LL_MeV){
        for(int iBinX=1; iBinX<=hSigma_LL->GetNbinsX()/Fraction_Sig; iBinX++){
            for(int iBinY=1; iBinY<=hSigma_LL->GetNbinsY()/Fraction_Sig; iBinY++){
                hSigma_LL_MeV->SetBinContent(iBinX, iBinY, hSigma_LL->GetBinContent(iBinX, iBinY));
            }
        }
    }
    if(hSigma_pXim && hSigma_pXim_MeV){
        for(int iBinX=1; iBinX<=hSigma_pXim->GetNbinsX()/Fraction_Sig; iBinX++){
            for(int iBinY=1; iBinY<=hSigma_pXim->GetNbinsY()/Fraction_Sig; iBinY++){
                hSigma_pXim_MeV->SetBinContent(iBinX, iBinY, hSigma_pXim->GetBinContent(iBinX, iBinY));
            }
        }
    }

    int vCutID;//which data file (cut combination) should you take. 0 = default
    //int vSource;//which source we use, see above
    int vFemReg_pp;//which femto region we use for pp (1 = default)
    int vFemReg_pL;//which femto region we use for pL (1 = default)
    int vFemReg_LL;//which femto region we use for LL (1 = default)
    int vFemReg_pXim;//which femto region we use for pXim (1 = default)
    int vBlReg;//which baseline region to use (1 = default)
    int vMod_pL;//which pL function to use: //0=exact NLO; 1=Ledni NLO; 2=Ledni LO; 3=ESC08
    int vFrac_pp_pL;//fraction of protons coming from Lambda variation (1 = default)
    int vFrac_pL_pSigma0;//fraction of Lambdas coming from Sigma0 variation (1 = default)
    int vFrac_pL_pXim;//fraction of Lambdas coming from Xim (1 = default)
    int vFit;//the shape of the baseline: 0=no BL (only renorm), 1=1st pol BL, 2=2nd pot BL
    int vStartPar_LL;//not important (starting fit values of the LL scattering parameters)
    int vSameRad;//0 = use different radii in each systems; 1 = fit with a common radius

    //each JOB produces a separate output file
    TFile* OutFile = new TFile(TString::Format("%sOutFile_%s_Iter%u_JOBS%u_ID%u.root",OutputDir.Data(),SystematicsType.Data(),NumIter,NumJobs,JobID), "recreate");
    //you save a lot of stuff in an NTuple
    TNtuple* ntResult = new TNtuple("ntResult", "ntResult", "IterID:vCutID:vSource:vFemReg_pp:vFemReg_pL:vFemReg_LL:vFemReg_pXim:"
                                    "vBlReg:vMod_pL:vFrac_pp_pL:vFrac_pL_pSigma0:vFrac_pL_pXim:vFit:vStartPar_LL:vSameRad:"
                                    "Radius_pp:RadiusErr_pp:Radius_pL:RadiusErr_pL:Radius_LL:RadiusErr_LL:Radius_pXim:RadiusErr_pXim:"
                                    "a0:a0Err:rEff:rEffErr:Chi2Ndf:pval");
    Float_t ntBuffer[29];

    unsigned WhichPart = JobID+1;
    unsigned NumJobsPart = NumIter;

    unsigned iSplitInto = NumJobs;
    unsigned FirstIter=0; unsigned LastIter=0;
    Printf("JobID = %i", JobID);
    while(iSplitInto>0 && WhichPart>0){
        FirstIter = LastIter + bool(LastIter);
        LastIter += NumJobsPart/iSplitInto-1+bool(LastIter);
        NumJobsPart = NumIter - LastIter - 1;
        iSplitInto--; WhichPart--;
    }
    if(NumIter==NumJobs) {FirstIter = JobID; LastIter=JobID;}

    Printf(" FirstIter = %i", FirstIter);
    Printf(" LastIter = %i", LastIter);
    Printf("  Nruns = %i", LastIter-FirstIter+1);
//***

//the 0 iter is always the default!
    for(unsigned uIter=FirstIter; uIter<=LastIter; uIter++){
        vCutID=(SystematicsType=="CutVarIterate")?rangen.Integer(31):0;
        vFemReg_pp=rangen.Integer(3);
        vFemReg_pL=rangen.Integer(3);
        vFemReg_LL=rangen.Integer(3);
        vFemReg_pXim=rangen.Integer(3);
        vBlReg=rangen.Integer(3);
        vFrac_pp_pL=rangen.Integer(3);
        vFrac_pL_pSigma0=rangen.Integer(3);
        vFrac_pL_pXim=rangen.Integer(3);
        vFit=rangen.Integer(2);
        vStartPar_LL=rangen.Integer(2);
        vSameRad=rangen.Integer(2);

        //The defaults
        if(uIter==0){
            vCutID = 0;
            vFemReg_pp = 1;
            vFemReg_pL = 1;
            vFemReg_LL = 1;
            vFemReg_pXim = 1;
            vBlReg = 1;
            vMod_pL = 2;
            vFrac_pp_pL = 1;
            vFrac_pL_pSigma0 = 1;
            vFrac_pL_pXim = 1;
            vFit = 1;
            vStartPar_LL = 1;
            vSameRad = 1;
        }
vMod_pL=0;
vSameRad=0;
vFit=0;
if(FitRegion==1) vFit=0;

if(vSource==6||vSource==8) vSameRad=0;
        printf("uIter=%u\n",uIter);
        printf("vCutID=%u\n",vCutID);
        printf("vSource=%u\n",vSource);
        printf("vFemReg_pp=%u\n",vFemReg_pp);
        printf("vFemReg_pL=%u\n",vFemReg_pL);
        printf("vFemReg_LL=%u\n",vFemReg_LL);
        printf("vFemReg_pXim=%u\n",vFemReg_pXim);
        printf("vBlReg=%u\n",vBlReg);
        printf("vMod_pL=%u\n",vMod_pL);
        printf("vFrac_pp_pL=%u\n",vFrac_pp_pL);
        printf("vFrac_pL_pSigma0=%u\n",vFrac_pL_pSigma0);
        printf("vFrac_pL_pXim=%u\n",vFrac_pL_pXim);
        printf("vFit=%u\n",vFit);
        printf("vStartPar_LL=%u\n",vStartPar_LL);
        printf("vSameRad=%u\n",vSameRad);


        //true => do the BL separately (as an RUN1), false => fit femto and BL region together (RUN2)
        bool SEPARATE_BL = false;
        //if true renormalization is NOT allowed
        bool FIX_CL = false;
FIX_CL = true;
        //note that for the time being the systematics folder is NOT updated automatically
        TString pp13Folder;TString pp13SystFolder;
        TString pPb5Folder;TString pPb5SystFolder;
        if(DataFlag==0){
            pp13Folder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample6/";
            pp13SystFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample6/";
            pPb5Folder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pPb_5TeV/Sample3/";
            pPb5SystFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pPb_5TeV/Sample3/";

        }
        else if(DataFlag==1){
            pp13Folder = "/home/dmihaylov/CernBox/Femto_pp13/data/";
            pp13SystFolder = "/home/dmihaylov/CernBox/Femto_pp13/";
            pPb5Folder = "/home/dmihaylov/CernBox/pPb/DataFullest/";
            pPb5SystFolder = "/home/dmihaylov/CernBox/SystematicsAndCalib/pPbRun2_MB/";
        }
        else{
            pp13Folder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/January2019_PaperDrafts/";
            pp13SystFolder = "/home/dmihaylov/CernBox/Femto_pp13/";
            pPb5Folder = "/home/dmihaylov/CernBox/pPb/DataFullest/";
            pPb5SystFolder = "/home/dmihaylov/CernBox/SystematicsAndCalib/pPbRun2_MB/";
        }

        //! DATA FILE
        TString OliFileName_pp;
        if(DataPeriod=="pp13TeV"&&DataFlag<2){
            OliFileName_pp =
                TString::Format("%sAnalysisResults_AndiBernieSystME_CkProtonProton_%u.root",pp13Folder.Data(),vCutID);
        }
        else if(DataPeriod=="pp13TeV"&&DataFlag==2){
            OliFileName_pp =
                TString::Format("%sCFOutput_pp.root",pp13Folder.Data());
        }
        else if(DataPeriod=="pPb5TeV"){
            OliFileName_pp =
                TString::Format("%sAnalysisResults_AndiBernieSystME_CkProtonProton_%u.root",pPb5Folder.Data(),vCutID);
        }

        TString OliHistoName_pp = "hCkTotNormWeight";
        if(DataPeriod=="pPb5TeV" && DataFlag==1) OliHistoName_pp = "hCk_ReweightedMeV_0";
        if(DataPeriod=="pp13TeV" && DataFlag==2) OliHistoName_pp = "hCk_ReweightedMeV_0";
        TFile* OliFile_pp = OliFileName_pp!=""?new TFile(OliFileName_pp, "read"):NULL;
        TH1F* OliHisto_pp = OliFile_pp?(TH1F*)OliFile_pp->Get(OliHistoName_pp):NULL;

        //!CHANGE PATH HERE
        TString SystErrFileName_pp =
            DataPeriod=="pp13TeV"?TString::Format("%sC2totalsysPP.root",pp13SystFolder.Data()):
            TString::Format("%sC2totalsysPP.root",pPb5SystFolder.Data());
        TString SystErrHistoName_pp = "C2totalsysPP";
        TFile* SystErrFile_pp = SystErrFileName_pp!=""?new TFile(SystErrFileName_pp, "read"):NULL;
        TH1F* SystErrHisto_pp = SystErrFile_pp?(TH1F*)SystErrFile_pp->Get(SystErrHistoName_pp):NULL;
        int NumSEB_pp = SystErrHisto_pp==NULL?0:SystErrHisto_pp->GetNbinsX();
        if(NumSEB_pp>OliHisto_pp->GetNbinsX()) NumSEB_pp=OliHisto_pp->GetNbinsX();
        for(int iBin=0; iBin<NumSEB_pp; iBin++){
            if(SystematicsType=="CutVarIterate") continue;
            OliHisto_pp->SetBinError(iBin+1, sqrt(pow(OliHisto_pp->GetBinError(iBin+1),2.)+pow(SystErrHisto_pp->GetBinContent(iBin+1),2.)));
        }

        //!CHANGE PATH HERE
        TString OliFileName_pL;
        if(DataPeriod=="pp13TeV"&&DataFlag<2){
            OliFileName_pL =
                TString::Format("%sAnalysisResults_AndiBernieSystME_CkProtonLambda_%u.root",pp13Folder.Data(),vCutID);
        }
        else if(DataPeriod=="pp13TeV"&&DataFlag==2){
            OliFileName_pL =
                TString::Format("%sCFOutput_pL.root",pp13Folder.Data());
        }
        else if(DataPeriod=="pPb5TeV"){
            OliFileName_pL =
                TString::Format("%sAnalysisResults_AndiBernieSystME_CkProtonLambda_%u.root",pPb5Folder.Data(),vCutID);
        }

        TString OliHistoName_pL = "hCkTotNormWeight";
        if(DataPeriod=="pPb5TeV" && DataFlag==1) OliHistoName_pL = "hCk_ReweightedMeV_1";
        if(DataPeriod=="pp13TeV" && DataFlag==2) OliHistoName_pL = "hCk_ReweightedMeV_1";
        TFile* OliFile_pL;
        if(OliFileName_pp==OliFileName_pL) OliFile_pL=OliFile_pp;
        else OliFile_pL = OliFileName_pL!=""?new TFile(OliFileName_pL, "read"):NULL;
        TH1F* OliHisto_pL = OliFile_pL?(TH1F*)OliFile_pL->Get(OliHistoName_pL):NULL;

        //!CHANGE PATH HERE
        TString SystErrFileName_pL =
            DataPeriod=="pp13TeV"?TString::Format("%sC2totalsysPL.root",pp13SystFolder.Data()):
            TString::Format("%sC2totalsysPL.root",pPb5SystFolder.Data());
        TString SystErrHistoName_pL = "C2totalsysPL";
        TFile* SystErrFile_pL = SystErrFileName_pL!=""?new TFile(SystErrFileName_pL, "read"):NULL;
        TH1F* SystErrHisto_pL = SystErrFile_pL?(TH1F*)SystErrFile_pL->Get(SystErrHistoName_pL):NULL;
        int NumSEB_pL = SystErrHisto_pL==NULL?0:SystErrHisto_pL->GetNbinsX();
        if(NumSEB_pL>OliHisto_pL->GetNbinsX()) NumSEB_pL=OliHisto_pL->GetNbinsX();
        for(int iBin=0; iBin<NumSEB_pL; iBin++){
            if(SystematicsType=="CutVarIterate") continue;
            OliHisto_pL->SetBinError(iBin+1, sqrt(pow(OliHisto_pL->GetBinError(iBin+1),2.)+pow(SystErrHisto_pL->GetBinContent(iBin+1),2.)));
            //printf("  Doing it to pLambda as well!\n");
        }

        //!CHANGE PATH HERE
        TString OliFileName_LL;
        if(DataPeriod=="pp13TeV"&&DataFlag<2){
            OliFileName_LL =
                TString::Format("%sAnalysisResults_AndiBernieSystME_CkLambdaLambda_%u.root",pp13Folder.Data(),vCutID);
        }
        else if(DataPeriod=="pp13TeV"&&DataFlag==2){
            OliFileName_LL =
                TString::Format("%sCFOutput_LL.root",pp13Folder.Data());
        }
        else if(DataPeriod=="pPb5TeV"){
            OliFileName_LL =
                TString::Format("%sAnalysisResults_AndiBernieSystME_CkLambdaLambda_%u.root",pPb5Folder.Data(),vCutID);
        }

        TString OliHistoName_LL = "hCkTotNormWeight";
        if(DataPeriod=="pPb5TeV" && DataFlag==1) OliHistoName_LL = "hCk_ReweightedMeV_1";
        if(DataPeriod=="pp13TeV" && DataFlag==2) OliHistoName_LL = "hCk_ReweightedMeV_1";
        TFile* OliFile_LL;
        if(OliFileName_pp==OliFileName_LL) OliFile_LL=OliFile_pp;
        else OliFile_LL = OliFileName_LL!=""?new TFile(OliFileName_LL, "read"):NULL;
        TH1F* OliHisto_LL = OliFile_LL?(TH1F*)OliFile_LL->Get(OliHistoName_LL):NULL;

        //!CHANGE PATH HERE
        TString SystErrFileName_LL =
            DataPeriod=="pp13TeV"?TString::Format("%sC2totalsysLL.root",pp13SystFolder.Data()):
            TString::Format("%sC2totalsysLL.root",pPb5SystFolder.Data());
        TString SystErrHistoName_LL = "C2totalsysLL";
        TFile* SystErrFile_LL = SystErrFileName_LL!=""?new TFile(SystErrFileName_LL, "read"):NULL;
        TH1F* SystErrHisto_LL = SystErrFile_LL?(TH1F*)SystErrFile_LL->Get(SystErrHistoName_LL):NULL;
        int NumSEB_LL = SystErrHisto_LL==NULL?0:SystErrHisto_LL->GetNbinsX();
        if(NumSEB_LL>OliHisto_LL->GetNbinsX()) NumSEB_LL=OliHisto_LL->GetNbinsX();
        for(int iBin=0; iBin<NumSEB_LL; iBin++){
            if(SystematicsType=="CutVarIterate") continue;
            OliHisto_LL->SetBinError(iBin+1, sqrt(pow(OliHisto_LL->GetBinError(iBin+1),2.)+pow(SystErrHisto_LL->GetBinContent(iBin+1),2.)));
        }

        //!CHANGE PATH HERE
        TString OliFileName_pXim;
        if(DataPeriod=="pp13TeV"&&DataFlag<2){
            OliFileName_pXim =
                TString::Format("%sAnalysisResults_AndiBernieSystME_CkProtonXim_%u.root",pp13Folder.Data(),vCutID);
        }
        else if(DataPeriod=="pp13TeV"&&DataFlag==2){
            OliFileName_pXim =
                TString::Format("%sCFOutput_pXi.root",pp13Folder.Data());
        }
        else if(DataPeriod=="pPb5TeV"){
            OliFileName_pXim =
                TString::Format("%sAnalysisResults_AndiBernieSystME_CkProtonXim_%u.root",pPb5Folder.Data(),vCutID);
        }

        if(ExcludeXi) OliFileName_pXim="";
        TString OliHistoName_pXim = "hCkTotNormWeight";
        if(DataPeriod=="pPb5TeV" && DataFlag==1) OliHistoName_pXim = "hCk_ReweightedMeV_1";
        if(DataPeriod=="pp13TeV" && DataFlag==2) OliHistoName_pXim = "hCk_ReweightedMeV_1";
//TString OliHistoName_pXim = "hCkPartNorm";
        TFile* OliFile_pXim;
        if(OliFileName_pp==OliFileName_pXim) OliFile_pXim=OliFile_pp;
        else OliFile_pXim = OliFileName_pXim!=""?new TFile(OliFileName_pXim, "read"):NULL;
        TH1F* OliHisto_pXim = OliFile_pXim?(TH1F*)OliFile_pXim->Get(OliHistoName_pXim):NULL;

        //!CHANGE PATH HERE
        TString SystErrFileName_pXim =
            DataPeriod=="pp13TeV"?TString::Format("%sC2totalsysPXi.root",pp13SystFolder.Data()):
            TString::Format("%sC2totalsysPXi.root",pPb5SystFolder.Data());
        TString SystErrHistoName_pXim = "C2totalsysPXi";
        if(ExcludeXi) SystErrFileName_pXim="";
        TFile* SystErrFile_pXim = SystErrFileName_pXim!=""?new TFile(SystErrFileName_pXim, "read"):NULL;
        TH1F* SystErrHisto_pXim = SystErrFile_pXim?(TH1F*)SystErrFile_pXim->Get(SystErrHistoName_pXim):NULL;
        int NumSEB_pXim = SystErrHisto_pXim==NULL?0:SystErrHisto_pXim->GetNbinsX();
        if(!ExcludeXi){
            if(NumSEB_pXim>OliHisto_pXim->GetNbinsX()) NumSEB_pXim=OliHisto_pXim->GetNbinsX();
            for(int iBin=0; iBin<NumSEB_pXim; iBin++){
                if(SystematicsType=="CutVarIterate") continue;
                OliHisto_pXim->SetBinError(iBin+1, sqrt(pow(OliHisto_pXim->GetBinError(iBin+1),2.)+pow(SystErrHisto_pXim->GetBinContent(iBin+1),2.)));
                //printf("  Doing it to pLambda as well!\n");
            }
        }


        //!The Femtoregions for pp/pL/LL/pXim
        //if you modify you may need to change the CATS ranges somewhere below
        double FemtoRegion_pp[3][2];
        double FemtoRegion_pL[3][2];
        double FemtoRegion_LL[3][2];
        double FemtoRegion_pXim[3][2];
        double BlRegion[3][2];
        //initial fit
        if(FitRegion==0){
            FemtoRegion_pp[0][0]=OliHisto_pp->GetBinLowEdge(1); FemtoRegion_pp[0][1]=OliHisto_pp->GetBinLowEdge(32);
            FemtoRegion_pp[1][0]=OliHisto_pp->GetBinLowEdge(1); FemtoRegion_pp[1][1]=OliHisto_pp->GetBinLowEdge(42);
            FemtoRegion_pp[2][0]=OliHisto_pp->GetBinLowEdge(1); FemtoRegion_pp[2][1]=OliHisto_pp->GetBinLowEdge(52);
        //printf("FemtoRegion_pp[0][1]=%f\n",FemtoRegion_pp[0][1]);
        //return;

            FemtoRegion_pL[0][0]=OliHisto_pL->GetBinLowEdge(1); FemtoRegion_pL[0][1]=OliHisto_pL->GetBinLowEdge(10);
            FemtoRegion_pL[1][0]=OliHisto_pL->GetBinLowEdge(1); FemtoRegion_pL[1][1]=OliHisto_pL->GetBinLowEdge(12);
            FemtoRegion_pL[2][0]=OliHisto_pL->GetBinLowEdge(1); FemtoRegion_pL[2][1]=OliHisto_pL->GetBinLowEdge(14);


            FemtoRegion_LL[0][0]=OliHisto_LL->GetBinLowEdge(1); FemtoRegion_LL[0][1]=OliHisto_LL->GetBinLowEdge(10);
            FemtoRegion_LL[1][0]=OliHisto_LL->GetBinLowEdge(1); FemtoRegion_LL[1][1]=OliHisto_LL->GetBinLowEdge(12);
            FemtoRegion_LL[2][0]=OliHisto_LL->GetBinLowEdge(1); FemtoRegion_LL[2][1]=OliHisto_LL->GetBinLowEdge(14);


            FemtoRegion_pXim[0][0]=OliHisto_pXim->GetBinLowEdge(1); FemtoRegion_pXim[0][1]=OliHisto_pXim->GetBinLowEdge(10);
            FemtoRegion_pXim[1][0]=OliHisto_pXim->GetBinLowEdge(1); FemtoRegion_pXim[1][1]=OliHisto_pXim->GetBinLowEdge(12);
            FemtoRegion_pXim[2][0]=OliHisto_pXim->GetBinLowEdge(1); FemtoRegion_pXim[2][1]=OliHisto_pXim->GetBinLowEdge(14);

            //!The baseline region (the same for all systems)
            BlRegion[0][0]=320; BlRegion[0][1]=480;
            BlRegion[1][0]=300; BlRegion[1][1]=500;
            BlRegion[2][0]=300; BlRegion[2][1]=540;
        }
        //as for the papers
        else{
            FemtoRegion_pp[0][0]=OliHisto_pp->GetBinLowEdge(1); FemtoRegion_pp[0][1]=OliHisto_pp->GetBinLowEdge(94);
            FemtoRegion_pp[1][0]=OliHisto_pp->GetBinLowEdge(1); FemtoRegion_pp[1][1]=OliHisto_pp->GetBinLowEdge(94);
            FemtoRegion_pp[2][0]=OliHisto_pp->GetBinLowEdge(1); FemtoRegion_pp[2][1]=OliHisto_pp->GetBinLowEdge(94);
        //printf("FemtoRegion_pp[0][1]=%f\n",FemtoRegion_pp[0][1]);
        //return;

            FemtoRegion_pL[0][0]=OliHisto_pL->GetBinLowEdge(1); FemtoRegion_pL[0][1]=OliHisto_pL->GetBinLowEdge(12);
            FemtoRegion_pL[1][0]=OliHisto_pL->GetBinLowEdge(1); FemtoRegion_pL[1][1]=OliHisto_pL->GetBinLowEdge(12);
            FemtoRegion_pL[2][0]=OliHisto_pL->GetBinLowEdge(1); FemtoRegion_pL[2][1]=OliHisto_pL->GetBinLowEdge(12);


            FemtoRegion_LL[0][0]=OliHisto_LL->GetBinLowEdge(1); FemtoRegion_LL[0][1]=OliHisto_LL->GetBinLowEdge(23);
            FemtoRegion_LL[1][0]=OliHisto_LL->GetBinLowEdge(1); FemtoRegion_LL[1][1]=OliHisto_LL->GetBinLowEdge(23);
            FemtoRegion_LL[2][0]=OliHisto_LL->GetBinLowEdge(1); FemtoRegion_LL[2][1]=OliHisto_LL->GetBinLowEdge(23);


            FemtoRegion_pXim[0][0]=OliHisto_pXim->GetBinLowEdge(1); FemtoRegion_pXim[0][1]=OliHisto_pXim->GetBinLowEdge(15);
            FemtoRegion_pXim[1][0]=OliHisto_pXim->GetBinLowEdge(1); FemtoRegion_pXim[1][1]=OliHisto_pXim->GetBinLowEdge(15);
            FemtoRegion_pXim[2][0]=OliHisto_pXim->GetBinLowEdge(1); FemtoRegion_pXim[2][1]=OliHisto_pXim->GetBinLowEdge(15);

            //!The baseline region (the same for all systems)
            BlRegion[0][0]=320; BlRegion[0][1]=480;
            BlRegion[1][0]=300; BlRegion[1][1]=500;
            BlRegion[2][0]=300; BlRegion[2][1]=540;

        }



        const unsigned NumSourcePars =  vSource==0?1:
                                        vSource==1?1:
                                        vSource==2?3:
                                        vSource==3?1:
                                        vSource==4?3:
                                        vSource==5?5:
                                        vSource==6?0:
                                        vSource==7?0:
                                        vSource==8?0:
                                        vSource==9?2:
                                                    2;

        //this way you define a correlation function using a CATS object.
        //needed inputs: num source/pot pars, CATS obj
        DLM_Ck* Ck_pp = new DLM_Ck(NumSourcePars,0,AB_pp);
        DLM_Ck* Ck_pL = vMod_pL==0?new DLM_Ck(NumSourcePars,0,AB_pL):new DLM_Ck(1,4,NumMomBins_pL,kMin_pL,kMax_pL,Lednicky_SingletTriplet);
        //this way you define a correlation function using Lednicky.
        //needed inputs: num source/pot pars, mom. binning, pointer to a function which computes C(k)
        DLM_Ck* Ck_LL = new DLM_Ck(1,2,NumMomBins_LL,kMin_LL,kMax_LL,Lednicky_Identical_Singlet);
        DLM_Ck* Ck_pSigma0 = new DLM_Ck(1,0,NumMomBins_pL,kMin_pL,kMax_pL,Lednicky_gauss_Sigma0);
        Ck_pSigma0->SetSourcePar(0,GaussSourceSize);
        //DLM_Ck* Ck_pXiMinus = new DLM_Ck(1,0,NumMomBins_pL,kMin_pL,kMax_pL,pXi_pheno);
        //Ck_pXiMinus->SetSourcePar(0,3.88);
        DLM_Ck* Ck_pXim = new DLM_Ck(NumSourcePars,0,AB_pXim);
        DLM_Ck* Ck_pXim1530 = new DLM_Ck(NumSourcePars,0,AB_pXim1530);

        Ck_pL->SetSourcePar(0,GaussSourceSize);
        if(vMod_pL==1){
            Ck_pL->SetPotPar(0,2.91);
            Ck_pL->SetPotPar(1,2.78);
            Ck_pL->SetPotPar(2,1.54);
            Ck_pL->SetPotPar(3,2.72);
        }
        else if(vMod_pL==2){
            Ck_pL->SetPotPar(0,1.91);
            Ck_pL->SetPotPar(1,1.4);
            Ck_pL->SetPotPar(2,1.23);
            Ck_pL->SetPotPar(3,2.13);
        }
        else if(vMod_pL==3){
            Ck_pL->SetPotPar(0,2.7);
            Ck_pL->SetPotPar(1,2.97);
            Ck_pL->SetPotPar(2,1.65);
            Ck_pL->SetPotPar(3,3.63);
        }

        Ck_LL->SetSourcePar(0,GaussSourceSize);
        Ck_LL->SetPotPar(0,1.2);
        Ck_LL->SetPotPar(1,4.5);
        Ck_pp->Update();
        Ck_pL->Update();
        Ck_LL->Update();
        Ck_pSigma0->Update();
        Ck_pXim->Update();
        Ck_pXim1530->Update();

        DLM_CkDecomposition CkDec_pp("pp",3,*Ck_pp,hSigma_pp_MeV);
        DLM_CkDecomposition CkDec_pL("pLambda",4,*Ck_pL,hSigma_pL_MeV);
        DLM_CkDecomposition CkDec_LL("LambdaLambda",2,*Ck_LL,hSigma_LL_MeV);
        DLM_CkDecomposition CkDec_pSigma0("pSigma0",0,*Ck_pSigma0,NULL);
        DLM_CkDecomposition CkDec_pXim("pXim",3,*Ck_pXim,ExcludeXi?NULL:hSigma_pXim_MeV);
        DLM_CkDecomposition CkDec_pXim1530("pXim1530",0,*Ck_pXim1530,NULL);

        double lam_pp =   Purities_p[vFrac_pp_pL][0]*Fraction_p[vFrac_pp_pL][0]*
                                Purities_p[vFrac_pp_pL][0]*Fraction_p[vFrac_pp_pL][0];
        double lam_pp_pL =    Purities_p[vFrac_pp_pL][0]*Fraction_p[vFrac_pp_pL][0]*
                                    Purities_p[vFrac_pp_pL][1]*Fraction_p[vFrac_pp_pL][1]*2;
        double lam_pp_fake =  Purities_p[vFrac_pp_pL][3]*Purities_p[vFrac_pp_pL][0]+
                                    Purities_p[vFrac_pp_pL][0]*Purities_p[vFrac_pp_pL][3]+
                                    Purities_p[vFrac_pp_pL][3]*Purities_p[vFrac_pp_pL][3];

        printf("lam_pp = %.3f\n", lam_pp);
        printf("lam_pp_pL = %.3f\n", lam_pp_pL);
        printf("lam_pp_fake = %.3f\n", lam_pp_fake);
        printf("\n");
        CkDec_pp.AddContribution(0,lam_pp_pL,DLM_CkDecomposition::cFeedDown,&CkDec_pL,hRes_pp_pL_MeV);
        CkDec_pp.AddContribution(1,1.-lam_pp-lam_pp_pL-lam_pp_fake,DLM_CkDecomposition::cFeedDown);
        CkDec_pp.AddContribution(2,lam_pp_fake,DLM_CkDecomposition::cFake);//0.02

        double lam_pL =   Purities_p[vFrac_pp_pL][0]*Fraction_p[vFrac_pp_pL][0]*
                                Purities_L[vFrac_pL_pSigma0][vFrac_pL_pXim][0]*Fraction_L[vFrac_pL_pSigma0][vFrac_pL_pXim][0];
        double lam_pL_pS0 =   Purities_p[vFrac_pp_pL][0]*Fraction_p[vFrac_pp_pL][0]*
                                    Purities_L[vFrac_pL_pSigma0][vFrac_pL_pXim][1]*Fraction_L[vFrac_pL_pSigma0][vFrac_pL_pXim][1];
        double lam_pL_pXm =   Purities_p[vFrac_pp_pL][0]*Fraction_p[vFrac_pp_pL][0]*
                                    Purities_L[vFrac_pL_pSigma0][vFrac_pL_pXim][2]*Fraction_L[vFrac_pL_pSigma0][vFrac_pL_pXim][2];
        double lam_pL_fake =  Purities_p[vFrac_pp_pL][3]*Purities_L[vFrac_pL_pSigma0][vFrac_pL_pXim][0]+
                                    Purities_p[vFrac_pp_pL][0]*Purities_L[vFrac_pL_pSigma0][vFrac_pL_pXim][4]+
                                    Purities_p[vFrac_pp_pL][3]*Purities_L[vFrac_pL_pSigma0][vFrac_pL_pXim][4];

        printf("lam_pL=%.3f\n",lam_pL);
        printf("lam_pL_pS0=%.3f\n",lam_pL_pS0);
        printf("lam_pL_pXm=%.3f\n",lam_pL_pXm);
        printf("lam_pL_fake=%.3f\n",lam_pL_fake);
        printf("\n");

        CkDec_pL.AddContribution(0,lam_pL_pS0,DLM_CkDecomposition::cFeedDown,&CkDec_pSigma0,hRes_pL_pSigma0_MeV);
        CkDec_pL.AddContribution(1,lam_pL_pXm,DLM_CkDecomposition::cFeedDown,&CkDec_pXim,hRes_pL_pXim_MeV);
        CkDec_pL.AddContribution(2,1.-lam_pL-lam_pL_pS0-lam_pL_pXm-lam_pL_fake,DLM_CkDecomposition::cFeedDown);
        CkDec_pL.AddContribution(3,lam_pL_fake,DLM_CkDecomposition::cFake);//0.03

        double lam_LL =   Purities_L[vFrac_pL_pSigma0][vFrac_pL_pXim][0]*Fraction_L[vFrac_pL_pSigma0][vFrac_pL_pXim][0]*
                                Purities_L[vFrac_pL_pSigma0][vFrac_pL_pXim][0]*Fraction_L[vFrac_pL_pSigma0][vFrac_pL_pXim][0];
        double lam_LL_fake =  Purities_L[vFrac_pL_pSigma0][vFrac_pL_pXim][4]*Purities_L[vFrac_pL_pSigma0][vFrac_pL_pXim][0]+
                                    Purities_L[vFrac_pL_pSigma0][vFrac_pL_pXim][0]*Purities_L[vFrac_pL_pSigma0][vFrac_pL_pXim][4]+
                                    Purities_L[vFrac_pL_pSigma0][vFrac_pL_pXim][4]*Purities_L[vFrac_pL_pSigma0][vFrac_pL_pXim][4];

        printf("lam_LL=%.6f\n",lam_LL);
        printf("lam_LL_fake=%.6f\n",lam_LL_fake);
        printf("\n");

        CkDec_LL.AddContribution(0,1.-lam_LL-lam_LL_fake,DLM_CkDecomposition::cFeedDown);//0.65
        CkDec_LL.AddContribution(1,lam_LL_fake,DLM_CkDecomposition::cFake);//0.03

        const double lam_pXim =   Purities_p[vFrac_pp_pL][0]*Fraction_p[vFrac_pp_pL][0]*
                                Purities_Xim[0][0]*Fraction_Xim[0][0];
        const double lam_pXim_pXim1530 =    Purities_p[vFrac_pp_pL][0]*Fraction_p[vFrac_pp_pL][0]*
                                    Purities_Xim[0][1]*Fraction_Xim[0][1];
        const double lam_pXim_fake =  Purities_p[vFrac_pp_pL][3]*Purities_Xim[0][0]+
                                    Purities_p[vFrac_pp_pL][0]*Purities_Xim[0][4]+
                                    Purities_p[vFrac_pp_pL][3]*Purities_Xim[0][4];

        printf("lam_pXim = %.3f\n", lam_pXim);
        printf("lam_pXim_pXim1530 = %.3f\n", lam_pXim_pXim1530);
        printf("lam_pXim_fake = %.3f\n", lam_pXim_fake);
        printf("\n");

        if(!ExcludeXi){
            CkDec_pXim.AddContribution(0,lam_pXim_pXim1530,DLM_CkDecomposition::cFeedDown,&CkDec_pXim1530,hRes_pXim_pXim1530_MeV);//from Xi-(1530)
            CkDec_pXim.AddContribution(1,1.-lam_pXim-lam_pXim_pXim1530-lam_pXim_fake,DLM_CkDecomposition::cFeedDown);//other feed-down (flat)
            CkDec_pXim.AddContribution(2,lam_pXim_fake,DLM_CkDecomposition::cFake);
        }

        DLM_Fitter1* fitter = new DLM_Fitter1(ExcludeXi?3:4);
        fitter->SetOutputDir("./OutputRun2/RUN2_SYSTEMATICS_Express1/Temp/");

        if(vFit==0){
            fitter->SetSystem(0,*OliHisto_pp,1,CkDec_pp,
                              FemtoRegion_pp[vFemReg_pp][0],FemtoRegion_pp[vFemReg_pp][1],
                              FemtoRegion_pp[vFemReg_pp][1],FemtoRegion_pp[vFemReg_pp][1]);
//printf("FEM REG: %f %f\n",FemtoRegion_pp[vFemReg_pp][0],FemtoRegion_pp[vFemReg_pp][1]);
            fitter->SetSystem(1,*OliHisto_pL,1,CkDec_pL,
                              FemtoRegion_pL[vFemReg_pL][0],FemtoRegion_pL[vFemReg_pL][1],
                              FemtoRegion_pL[vFemReg_pL][1],FemtoRegion_pL[vFemReg_pL][1]);

            fitter->SetSystem(2,*OliHisto_LL,1,CkDec_LL,
                              FemtoRegion_LL[vFemReg_LL][0],FemtoRegion_LL[vFemReg_LL][1],
                              FemtoRegion_LL[vFemReg_LL][1],FemtoRegion_LL[vFemReg_LL][1]);

            if(!ExcludeXi){
                fitter->SetSystem(3,*OliHisto_pXim,1,CkDec_pXim,
                                  FemtoRegion_pXim[vFemReg_pXim][0],FemtoRegion_pXim[vFemReg_pXim][1],
                                  FemtoRegion_pXim[vFemReg_pXim][1],FemtoRegion_pXim[vFemReg_pXim][1]);
            }
        }
        else{
            fitter->SetSystem(0,*OliHisto_pp,1,CkDec_pp,
                              FemtoRegion_pp[vFemReg_pp][0],FemtoRegion_pp[vFemReg_pp][1],
                              BlRegion[vBlReg][0],BlRegion[vBlReg][1]);

            fitter->SetSystem(1,*OliHisto_pL,1,CkDec_pL,
                              FemtoRegion_pL[vFemReg_pL][0],FemtoRegion_pL[vFemReg_pL][1],
                              BlRegion[vBlReg][0],BlRegion[vBlReg][1]);
            fitter->SetSystem(2,*OliHisto_LL,1,CkDec_LL,
                              FemtoRegion_LL[vFemReg_LL][0],FemtoRegion_LL[vFemReg_LL][1],
                              BlRegion[vBlReg][0],BlRegion[vBlReg][1]);

            if(!ExcludeXi){
                fitter->SetSystem(3,*OliHisto_pXim,1,CkDec_pXim,
                                  FemtoRegion_pXim[vFemReg_pXim][0],FemtoRegion_pXim[vFemReg_pXim][1],
                                  BlRegion[vBlReg][0],BlRegion[vBlReg][1]);
            }
        }

        if(SEPARATE_BL){
            fitter->SetSeparateBL(0,vFit!=0);
            fitter->SetSeparateBL(1,vFit!=0);
            fitter->SetSeparateBL(2,vFit!=0);
            if(!ExcludeXi) fitter->SetSeparateBL(3,vFit!=0);
        }
        else{
            fitter->SetSeparateBL(0,false);
            fitter->SetSeparateBL(1,false);
            fitter->SetSeparateBL(2,false);
            if(!ExcludeXi) fitter->SetSeparateBL(3,false);
        }

        if(vSameRad==0){
            if(vSource==6||vSource==7||vSource==8){

            }
            else{
                if(vSource==0) fitter->AddSameSource("LambdaLambda","pp",1);
                fitter->AddSameSource("pSigma0","pLambda",1);
                if(!ExcludeXi) fitter->AddSameSource("pXim1530","pp",1);
            }
        }
        else{
            fitter->AddSameSource("pLambda","pp",1);
            if(vSource==0) fitter->AddSameSource("LambdaLambda","pp",1);
            fitter->AddSameSource("pSigma0","pp",1);
            fitter->AddSameSource("pXim","pp",1);
            if(!ExcludeXi) fitter->AddSameSource("pXim1530","pp",1);
        }

        fitter->SetParameter("pp",DLM_Fitter1::p_a,1.0,0.7,1.3);
        if(vFit<1) fitter->FixParameter("pp",DLM_Fitter1::p_b,0);
        else fitter->SetParameter("pp",DLM_Fitter1::p_b,1e-4,0,2e-3);
        if(vFit<2) fitter->FixParameter("pp",DLM_Fitter1::p_c,0);
        else fitter->SetParameter("pp",DLM_Fitter1::p_c,0,-2e-4,2e-4);
        if(vSource==0){
            fitter->SetParameter("pp",DLM_Fitter1::p_sor0,1.2,0.8,1.8);
//fitter->SetParameter("pp",DLM_Fitter1::p_sor0,0.9,0.5,1.8);
//fitter->FixParameter("pp",DLM_Fitter1::p_sor0,1.08);
        }
        else if(vSource==1){
            fitter->SetParameter("pp",DLM_Fitter1::p_sor0,GaussSourceSize/1.4,0.5,1.6);
        }
        else if(vSource==2){
            fitter->SetParameter("pp",DLM_Fitter1::p_sor0,0.8,0.4,8);
            fitter->SetParameter("pp",DLM_Fitter1::p_sor1,3.2,0.4,8);
            fitter->SetParameter("pp",DLM_Fitter1::p_sor2,0.5,0,1);
        }
        else if(vSource==3){
            fitter->SetParameter("pp",DLM_Fitter1::p_sor0,1.0,1.0,1.3);
//fitter->FixParameter("pp",DLM_Fitter1::p_sor0,1.0);
        }
        else if(vSource==4){
            fitter->SetParameter("pp",DLM_Fitter1::p_sor0,1.2,0.4,8);
//fitter->FixParameter("pp",DLM_Fitter1::p_sor0,2.0);
            fitter->SetParameter("pp",DLM_Fitter1::p_sor1,1.2,0.4,8);
            fitter->SetParameter("pp",DLM_Fitter1::p_sor2,0.5,0,1);
        }
        else if(vSource==5){
            fitter->SetParameter("pp",DLM_Fitter1::p_sor0,1.5*1.5,0.4*1.5,3.6*1.5);
            fitter->SetParameter("pp",DLM_Fitter1::p_sor1,1.5*1.5,0.4*1.5,3.6*1.5);
            fitter->SetParameter("pp",DLM_Fitter1::p_sor2,0.8*1.5,0.4*1.5,3.6*1.5);
//fitter->FixParameter("pp",DLM_Fitter1::p_sor0,1.6*1.5);
//fitter->FixParameter("pp",DLM_Fitter1::p_sor1,1.6*1.5);
//fitter->FixParameter("pp",DLM_Fitter1::p_sor2,0.8*1.5);
            fitter->FixParameter("pp",DLM_Fitter1::p_sor3,1);
            fitter->FixParameter("pp",DLM_Fitter1::p_sor4,140);
        }
        else if(vSource==9){
            fitter->SetParameter("pp",DLM_Fitter1::p_sor0,1.2,0.4,2);
            fitter->SetParameter("pp",DLM_Fitter1::p_sor1,1.5,1.0,2.0);
            fitter->FixParameter("pp",DLM_Fitter1::p_sor2,3);
        }
        else if(vSource==10){
            fitter->SetParameter("pp",DLM_Fitter1::p_sor0,1.0,0.3,1.7);
            fitter->SetParameter("pp",DLM_Fitter1::p_sor1,1.5,1.0,2.0);
        }

        if(FIX_CL) fitter->FixParameter("pp",DLM_Fitter1::p_Cl,-1);
        else fitter->SetParameter("pp",DLM_Fitter1::p_Cl,-0.9,-1.2,-0.8);

        fitter->SetParameter("pLambda",DLM_Fitter1::p_a,1.0,0.7,1.3);
        if(vFit<1) fitter->FixParameter("pLambda",DLM_Fitter1::p_b,0);
        else fitter->SetParameter("pLambda",DLM_Fitter1::p_b,1e-4,0,2e-3);
        if(vFit<2) fitter->FixParameter("pLambda",DLM_Fitter1::p_c,0);
        else fitter->SetParameter("pLambda",DLM_Fitter1::p_c,0,-2e-4,2e-4);
        if(vSameRad==false){
            if(vSource==0){
                fitter->SetParameter("pLambda",DLM_Fitter1::p_sor0,GaussSourceSize,0.8,1.8);
//fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,1.12);
            }
            else if(vSource==1){
                fitter->SetParameter("pLambda",DLM_Fitter1::p_sor0,GaussSourceSize/1.4,0.5,1.6);
            }
            else if(vSource==2){
                fitter->SetParameter("pLambda",DLM_Fitter1::p_sor0,0.8,0.4,3.6);
                fitter->SetParameter("pLambda",DLM_Fitter1::p_sor1,3.2,0.4,3.6);
                fitter->SetParameter("pLambda",DLM_Fitter1::p_sor2,0.5,0,1);
            }
            else if(vSource==3){
                fitter->SetParameter("pLambda",DLM_Fitter1::p_sor0,1,0.75,2);
//fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,1);
            }
            else if(vSource==4){
                fitter->SetParameter("pLambda",DLM_Fitter1::p_sor0,1.2,0.4,3.6);
//fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,2.0);
                fitter->SetParameter("pLambda",DLM_Fitter1::p_sor1,1.2,0.4,3.6);
                fitter->SetParameter("pLambda",DLM_Fitter1::p_sor2,0.5,0,1);
            }
            else if(vSource==5){
                fitter->SetParameter("pLambda",DLM_Fitter1::p_sor0,1.2,0.4,3.6);
                fitter->SetParameter("pLambda",DLM_Fitter1::p_sor1,1.2,0.4,3.6);
                fitter->SetParameter("pLambda",DLM_Fitter1::p_sor2,1.2,0.4,3.6);
                fitter->FixParameter("pLambda",DLM_Fitter1::p_sor3,1);
                fitter->FixParameter("pLambda",DLM_Fitter1::p_sor4,140);
            }
            else if(vSource==9){
                fitter->SetParameter("pLambda",DLM_Fitter1::p_sor0,1.2,0.4,2);
                fitter->SetParameter("pLambda",DLM_Fitter1::p_sor1,1.5,1.0,2.0);
                fitter->FixParameter("pLambda",DLM_Fitter1::p_sor2,3);
            }
            else if(vSource==10){
                fitter->SetParameter("pLambda",DLM_Fitter1::p_sor0,1.2,0.4,2);
                fitter->SetParameter("pLambda",DLM_Fitter1::p_sor1,1.5,1.0,2.0);
            }
        }
        if(FIX_CL) fitter->FixParameter("pLambda",DLM_Fitter1::p_Cl,-1);
        else fitter->SetParameter("pLambda",DLM_Fitter1::p_Cl,-0.9,-1.2,-0.8);
        if(vMod_pL==1){
            fitter->FixParameter("pLambda",DLM_Fitter1::p_pot0,2.91);
            fitter->FixParameter("pLambda",DLM_Fitter1::p_pot1,2.78);
            fitter->FixParameter("pLambda",DLM_Fitter1::p_pot2,1.54);
            fitter->FixParameter("pLambda",DLM_Fitter1::p_pot3,2.72);
        }
        else if(vMod_pL==2){
            fitter->FixParameter("pLambda",DLM_Fitter1::p_pot0,1.91);
            fitter->FixParameter("pLambda",DLM_Fitter1::p_pot1,1.4);
            fitter->FixParameter("pLambda",DLM_Fitter1::p_pot2,1.23);
            fitter->FixParameter("pLambda",DLM_Fitter1::p_pot3,2.13);
        }
        else if(vMod_pL==3){
            fitter->FixParameter("pLambda",DLM_Fitter1::p_pot0,2.7);
            fitter->FixParameter("pLambda",DLM_Fitter1::p_pot1,2.97);
            fitter->FixParameter("pLambda",DLM_Fitter1::p_pot2,1.65);
            fitter->FixParameter("pLambda",DLM_Fitter1::p_pot3,3.63);
        }

        fitter->SetParameter("LambdaLambda",DLM_Fitter1::p_a,1.0,0.7,1.3);
        if(vFit<1) fitter->FixParameter("LambdaLambda",DLM_Fitter1::p_b,0);
        else fitter->SetParameter("LambdaLambda",DLM_Fitter1::p_b,1e-4,0,2e-3);
        if(vFit<2) fitter->FixParameter("LambdaLambda",DLM_Fitter1::p_c,0);
        else fitter->SetParameter("LambdaLambda",DLM_Fitter1::p_c,0,-2e-4,2e-4);
        if(vSource!=0) fitter->FixParameter("LambdaLambda",DLM_Fitter1::p_sor0,1.2);//!
        //if(FIX_CL) fitter->FixParameter("LambdaLambda",DLM_Fitter1::p_Cl,-1);
        //else fitter->SetParameter("LambdaLambda",DLM_Fitter1::p_Cl,-0.9,-1.2,-0.8);
        //!fix Cl for now, since the LambdaLambda is very shallow end prob. there is very little deviation from unity.
        fitter->FixParameter("LambdaLambda",DLM_Fitter1::p_Cl,-1);
        if(vStartPar_LL) fitter->SetParameter("LambdaLambda",DLM_Fitter1::p_pot0,4,-8,24);
        else fitter->SetParameter("LambdaLambda",DLM_Fitter1::p_pot0,-4,-8,24);
        fitter->SetParameter("LambdaLambda",DLM_Fitter1::p_pot1,8,0,24);

        if(!ExcludeXi){
            fitter->SetParameter("pXim",DLM_Fitter1::p_a,1.0,0.7,1.3);
            if(vFit<1) fitter->FixParameter("pXim",DLM_Fitter1::p_b,0);
            else fitter->SetParameter("pXim",DLM_Fitter1::p_b,1e-4,0,2e-3);
            if(vFit<2) fitter->FixParameter("pXim",DLM_Fitter1::p_c,0);
            else fitter->SetParameter("pXim",DLM_Fitter1::p_c,0,-2e-4,2e-4);
            if(vSameRad==false){
                if(vSource==0){
                    fitter->SetParameter("pXim",DLM_Fitter1::p_sor0,GaussSourceSize,0.8,1.8);
                }
                else if(vSource==1){
                    fitter->SetParameter("pXim",DLM_Fitter1::p_sor0,GaussSourceSize/1.4,0.5,1.6);
                }
                else if(vSource==2){
                    fitter->SetParameter("pXim",DLM_Fitter1::p_sor0,0.8,0.4,3.6);
                    fitter->SetParameter("pXim",DLM_Fitter1::p_sor1,3.2,0.4,3.6);
                    fitter->SetParameter("pXim",DLM_Fitter1::p_sor2,0.5,0,1);
                }
                else if(vSource==3){
                    fitter->SetParameter("pXim",DLM_Fitter1::p_sor0,1,0.75,2);
//fitter->FixParameter("pXim",DLM_Fitter1::p_sor0,1);
                }
                else if(vSource==4){
                    fitter->SetParameter("pXim",DLM_Fitter1::p_sor0,1.2,0.4,3.6);
//fitter->FixParameter("pXim",DLM_Fitter1::p_sor0,2.0);
                    fitter->SetParameter("pXim",DLM_Fitter1::p_sor1,1.2,0.4,3.6);
                    fitter->SetParameter("pXim",DLM_Fitter1::p_sor2,0.5,0,1);
                }
                else if(vSource==5){
                    fitter->SetParameter("pXim",DLM_Fitter1::p_sor0,1.2,0.4,3.6);
                    fitter->SetParameter("pXim",DLM_Fitter1::p_sor1,1.2,0.4,3.6);
                    fitter->SetParameter("pXim",DLM_Fitter1::p_sor2,1.2,0.4,3.6);
                    fitter->FixParameter("pXim",DLM_Fitter1::p_sor3,1);
                    fitter->FixParameter("pXim",DLM_Fitter1::p_sor4,140);
                }
                else if(vSource==9){
                    fitter->SetParameter("pXim",DLM_Fitter1::p_sor0,1.2,0.4,2);
                    fitter->SetParameter("pXim",DLM_Fitter1::p_sor1,1.5,1.0,2.0);
                    fitter->FixParameter("pXim",DLM_Fitter1::p_sor2,3);
                }
                else if(vSource==10){
                    fitter->SetParameter("pXim",DLM_Fitter1::p_sor0,1.2,0.4,2);
                    fitter->SetParameter("pXim",DLM_Fitter1::p_sor1,1.5,1.0,2.0);
                }
            }
            if(FIX_CL) fitter->FixParameter("pXim",DLM_Fitter1::p_Cl,-1);
            else fitter->SetParameter("pXim",DLM_Fitter1::p_Cl,-0.9,-1.2,-0.8);
        }

//fitter->FixParameter("pp",DLM_Fitter1::p_sor0,1.0);
//fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,1.0);
//fitter->FixParameter("LambdaLambda",DLM_Fitter1::p_sor0,1.2);
//fitter->FixParameter("pXim",DLM_Fitter1::p_sor0,1.0);

//fitter->FixParameter("pp",DLM_Fitter1::p_sor1,1.80);
//fitter->FixParameter("pLambda",DLM_Fitter1::p_sor1,1.55);

//! introducing this breaks the cat after 5.22 M function calls
//! I have added some additional output in CATS at line 2711 and 2774, take it out when you fix it
//! possible places to look for a fix: the whole DLM_Fitters.cpp, how you fix the parameters, pass them on to the Ck objects etc.
//okay, I believe it is the stable parameter that goes insane! it goes to 0.4 or so, so probably you should set some limits on that!
fitter->FixParameter("pp",DLM_Fitter1::p_sor0,1.21);
fitter->FixParameter("pp",DLM_Fitter1::p_sor1,1.6);
fitter->FixParameter("pLambda",DLM_Fitter1::p_sor1,1.3);
fitter->FixParameter("pXim",DLM_Fitter1::p_sor1,1.8);

//AddSameParameter(const TString& System, const unsigned& WhichPar, const TString& ParentSystem, const unsigned& ParentPar)
fitter->AddSameParameter("pLambda",DLM_Fitter1::p_sor0,"pp",DLM_Fitter1::p_sor0);
fitter->AddSameParameter("pXim",DLM_Fitter1::p_sor0,"pp",DLM_Fitter1::p_sor0);

/*
const double BigSource = 2.1;
const double SmallSource = 1.1;
const double DilutedProtons = 0.4;
fitter->FixParameter("pp",DLM_Fitter1::p_sor0,BigSource);
fitter->FixParameter("pp",DLM_Fitter1::p_sor1,SmallSource);
fitter->FixParameter("pp",DLM_Fitter1::p_sor2,DilutedProtons);
fitter->FixParameter("pLambda",DLM_Fitter1::p_sor0,BigSource);
fitter->FixParameter("pLambda",DLM_Fitter1::p_sor1,SmallSource);
fitter->FixParameter("pLambda",DLM_Fitter1::p_sor2,DilutedProtons*0.5);
fitter->FixParameter("pXim",DLM_Fitter1::p_sor0,BigSource);
fitter->FixParameter("pXim",DLM_Fitter1::p_sor1,SmallSource);
fitter->FixParameter("pXim",DLM_Fitter1::p_sor2,DilutedProtons*0.5);
*/
        CkDec_pp.Update();
        CkDec_pL.Update();
        CkDec_LL.Update();
        if(!ExcludeXi) CkDec_pXim.Update();
        fitter->GoBabyGo();

        ntBuffer[0]=uIter;
        ntBuffer[1]=vCutID;
        ntBuffer[2]=vSource;
        ntBuffer[3]=vFemReg_pp;
        ntBuffer[4]=vFemReg_pL;
        ntBuffer[5]=vFemReg_LL;
        ntBuffer[6]=vFemReg_pXim;
        ntBuffer[7]=vBlReg;
        ntBuffer[8]=vMod_pL;
        ntBuffer[9]=vFrac_pp_pL;
        ntBuffer[10]=vFrac_pL_pSigma0;
        ntBuffer[11]=vFrac_pL_pXim;
        ntBuffer[12]=vFit;
        ntBuffer[13]=vStartPar_LL;
        ntBuffer[14]=vSameRad;
        ntBuffer[15]=fitter->GetParameter("pp",DLM_Fitter1::p_sor0);
        ntBuffer[16]=fitter->GetParError("pp",DLM_Fitter1::p_sor0);
        ntBuffer[17]=fitter->GetParameter("pLambda",DLM_Fitter1::p_sor0);
        ntBuffer[18]=fitter->GetParError("pLambda",DLM_Fitter1::p_sor0);
        ntBuffer[19]=fitter->GetParameter("LambdaLambda",DLM_Fitter1::p_sor0);
        ntBuffer[20]=fitter->GetParError("LambdaLambda",DLM_Fitter1::p_sor0);
        ntBuffer[21]=ExcludeXi?0:fitter->GetParameter("pXim",DLM_Fitter1::p_sor0);
        ntBuffer[22]=ExcludeXi?0:fitter->GetParError("pXim",DLM_Fitter1::p_sor0);
        ntBuffer[23]=fitter->GetParameter("LambdaLambda",DLM_Fitter1::p_pot0);
        ntBuffer[24]=fitter->GetParError("LambdaLambda",DLM_Fitter1::p_pot0);
        ntBuffer[25]=fitter->GetParameter("LambdaLambda",DLM_Fitter1::p_pot1);
        ntBuffer[26]=fitter->GetParError("LambdaLambda",DLM_Fitter1::p_pot1);
        ntBuffer[27]=fitter->GetChi2Ndf();
        ntBuffer[28]=fitter->GetPval();
        ntResult->Fill(ntBuffer);

        TFile* GraphFile = new TFile(TString::Format("%sGraphFile_%s_Iter%u_uIter%u.root",OutputDir.Data(),SystematicsType.Data(),NumIter,uIter), "recreate");

        GraphFile->cd();

        TGraph FitResult_pp;
        FitResult_pp.SetName(TString::Format("FitResult_pp_%u",uIter));
        fitter->GetFitGraph(0, FitResult_pp);

        TGraph FitResult_pL;
        FitResult_pL.SetName(TString::Format("FitResult_pL_%u",uIter));
        fitter->GetFitGraph(1, FitResult_pL);

        TGraph FitResult_LL;
        FitResult_LL.SetName(TString::Format("FitResult_LL_%u",uIter));
        fitter->GetFitGraph(2, FitResult_LL);

        TGraph FitResult_pXim;
        FitResult_pXim.SetName(TString::Format("FitResult_pXim_%u",uIter));
        if(!ExcludeXi) fitter->GetFitGraph(3, FitResult_pXim);

        double Chi2 = fitter->GetChi2();
        unsigned NDF = fitter->GetNdf();
        double RadiusResult = fitter->GetParameter("pp",DLM_Fitter1::p_sor0);
        if(vSource==3) RadiusResult*=TRANSP_RENORM;
        double RadiusError = fitter->GetParError("pp",DLM_Fitter1::p_sor0);
        double ClResult_pp = fitter->GetParameter("pp",DLM_Fitter1::p_Cl);
        double ClError_pp = fitter->GetParError("pp",DLM_Fitter1::p_Cl);
        double RadiusResult2 = fitter->GetParameter("pp",DLM_Fitter1::p_sor1);
        double RadiusError2 = fitter->GetParError("pp",DLM_Fitter1::p_sor1);
        double RadiusResult3 = fitter->GetParameter("pp",DLM_Fitter1::p_sor2);
        double RadiusError3 = fitter->GetParError("pp",DLM_Fitter1::p_sor2);
        double RadiusResult_pL = fitter->GetParameter("pLambda",DLM_Fitter1::p_sor0);
        if(vSource==3) RadiusResult_pL*=TRANSP_RENORM;
        double RadiusError_pL = fitter->GetParError("pLambda",DLM_Fitter1::p_sor0);
        double ClResult_pL = fitter->GetParameter("pLambda",DLM_Fitter1::p_Cl);
        double ClError_pL = fitter->GetParError("pLambda",DLM_Fitter1::p_Cl);
        double RadiusResult2_pL = fitter->GetParameter("pLambda",DLM_Fitter1::p_sor1);
        double RadiusError2_pL = fitter->GetParError("pLambda",DLM_Fitter1::p_sor1);
        double RadiusResult3_pL = fitter->GetParameter("pLambda",DLM_Fitter1::p_sor2);
        double RadiusError3_pL = fitter->GetParError("pLambda",DLM_Fitter1::p_sor2);
        double RadiusResult_LL = fitter->GetParameter("LambdaLambda",DLM_Fitter1::p_sor0);
        double RadiusError_LL = fitter->GetParError("LambdaLambda",DLM_Fitter1::p_sor0);
        double ClResult_LL = fitter->GetParameter("LambdaLambda",DLM_Fitter1::p_Cl);
        double ClError_LL = fitter->GetParError("LambdaLambda",DLM_Fitter1::p_Cl);
        double RadiusResult_pXim = fitter->GetParameter("pXim",DLM_Fitter1::p_sor0);
        if(vSource==3) RadiusResult_pXim*=TRANSP_RENORM;
        double RadiusError_pXim = fitter->GetParError("pXim",DLM_Fitter1::p_sor0);
        double RadiusResult2_pXim = fitter->GetParameter("pXim",DLM_Fitter1::p_sor1);
        double RadiusError2_pXim = fitter->GetParError("pXim",DLM_Fitter1::p_sor1);
        double RadiusResult3_pXim = fitter->GetParameter("pXim",DLM_Fitter1::p_sor2);
        double RadiusError3_pXim = fitter->GetParError("pXim",DLM_Fitter1::p_sor2);
        double ClResult_pXim = fitter->GetParameter("pXim",DLM_Fitter1::p_Cl);
        double ClError_pXim = fitter->GetParError("pXim",DLM_Fitter1::p_Cl);

        if(Chi2/double(NDF)!=fitter->GetChi2Ndf()){printf("Oh boy...\n");}

        double Chi2_pp=0;
        unsigned EffNumBins_pp=0;
        for(unsigned uBin=0; uBin<NumMomBins_pp; uBin++){

            double mom = AB_pp.GetMomentum(uBin);
            double dataY;
            double dataErr;
            double theoryX;
            double theoryY;

            if(mom>FemtoRegion_pp[vFemReg_pp][1]) continue;

            FitResult_pp.GetPoint(uBin,theoryX,theoryY);
            //if(mom!=theoryX) printf("PROBLEM! (%f <-> %f)\n",mom,theoryX);
            dataY = OliHisto_pp->GetBinContent(uBin+1);
            dataErr = OliHisto_pp->GetBinError(uBin+1);

            Chi2_pp += (dataY-theoryY)*(dataY-theoryY)/(dataErr*dataErr);
            EffNumBins_pp++;
        }

        double Chi2_pL=0;
        unsigned EffNumBins_pL=0;
        for(unsigned uBin=0; uBin<NumMomBins_pL; uBin++){

            double mom = AB_pL.GetMomentum(uBin);
            //double dataX;
            double dataY;
            double dataErr;
            double theoryX;
            double theoryY;

            if(mom>FemtoRegion_pL[vFemReg_pL][1]) continue;

            FitResult_pL.GetPoint(uBin,theoryX,theoryY);
            //if(mom!=theoryX) printf("PROBLEM! (%f <-> %f)\n",mom,theoryX);
            dataY = OliHisto_pL->GetBinContent(uBin+1);
            dataErr = OliHisto_pL->GetBinError(uBin+1);

            Chi2_pL += (dataY-theoryY)*(dataY-theoryY)/(dataErr*dataErr);
            EffNumBins_pL++;
        }

        double Chi2_LL=0;
        unsigned EffNumBins_LL=0;
        for(unsigned uBin=0; uBin<NumMomBins_LL; uBin++){

            double mom = Ck_LL->GetBinCenter(0,uBin);
            //double dataX;
            double dataY;
            double dataErr;
            double theoryX;
            double theoryY;

            if(mom>FemtoRegion_LL[vFemReg_LL][1]) continue;

            FitResult_LL.GetPoint(uBin,theoryX,theoryY);
            //if(mom!=theoryX) printf("PROBLEM! (%f <-> %f)\n",mom,theoryX);
            dataY = OliHisto_LL->GetBinContent(uBin+1);
            dataErr = OliHisto_LL->GetBinError(uBin+1);

            Chi2_LL += (dataY-theoryY)*(dataY-theoryY)/(dataErr*dataErr);
            EffNumBins_LL++;
        }
        double a0_LL = fitter->GetParameter("LambdaLambda",DLM_Fitter1::p_pot0);
        double a0_err_LL = fitter->GetParError("LambdaLambda",DLM_Fitter1::p_pot0);
        double reff_LL = fitter->GetParameter("LambdaLambda",DLM_Fitter1::p_pot1);
        double reff_err_LL = fitter->GetParError("LambdaLambda",DLM_Fitter1::p_pot1);

        double Chi2_pXim=0;
        unsigned EffNumBins_pXim=0;
        for(unsigned uBin=0; uBin<NumMomBins_pXim; uBin++){

            double mom = AB_pXim.GetMomentum(uBin);
            //double dataX;
            double dataY;
            double dataErr;
            double theoryX;
            double theoryY;

            if(mom>FemtoRegion_pXim[vFemReg_pXim][1]) continue;

            FitResult_pXim.GetPoint(uBin,theoryX,theoryY);
            //if(mom!=theoryX) printf("PROBLEM!\n");
            dataY = OliHisto_pXim->GetBinContent(uBin+1);
            dataErr = OliHisto_pXim->GetBinError(uBin+1);

            Chi2_pXim += (dataY-theoryY)*(dataY-theoryY)/(dataErr*dataErr);
            EffNumBins_pXim++;
        }

        printf("Chi2_pp/bins = %.2f/%u = %.2f\n", Chi2_pp, EffNumBins_pp, Chi2_pp/double(EffNumBins_pp));
        printf("Chi2_pL/bins = %.2f/%u = %.2f\n", Chi2_pL, EffNumBins_pL, Chi2_pL/double(EffNumBins_pL));
        printf("Chi2_LL/bins = %.2f/%u = %.2f\n", Chi2_LL, EffNumBins_LL, Chi2_LL/double(EffNumBins_LL));
        printf("Chi2_pXim/bins = %.2f/%u = %.2f\n", Chi2_pXim, EffNumBins_pXim, Chi2_pXim/double(EffNumBins_pXim));

        GraphFile->cd();
        FitResult_pp.Write();
        FitResult_pL.Write();
        FitResult_LL.Write();
        FitResult_pXim.Write();
        fitter->GetFit()->SetName(TString::Format("GlobalFit_%u",uIter));
        fitter->GetFit()->Write();

        //Fix all parameters to their current values
        for(unsigned uPar=0; uPar<fitter->GetNumParPerSyst(); uPar++){
            fitter->FixParameter("pp",uPar);
            fitter->FixParameter("pLambda",uPar);
            fitter->FixParameter("LambdaLambda",uPar);
            fitter->FixParameter("pXim",uPar);
        }

        //change the pLambda to LO (will have an effect only if fitting with Lednicky)
        Ck_pL->SetPotPar(0,1.91);
        Ck_pL->SetPotPar(1,1.4);
        Ck_pL->SetPotPar(2,1.23);
        Ck_pL->SetPotPar(3,2.13);
        Ck_pL->SetSourcePar(0,vSameRad==true?fitter->GetParameter("pp",DLM_Fitter1::p_sor0):fitter->GetParameter("pLambda",DLM_Fitter1::p_sor0));
        Ck_pL->Update();

        fitter->FixParameter("pLambda",DLM_Fitter1::p_pot0,1.91);
        fitter->FixParameter("pLambda",DLM_Fitter1::p_pot1,1.4);
        fitter->FixParameter("pLambda",DLM_Fitter1::p_pot2,1.23);
        fitter->FixParameter("pLambda",DLM_Fitter1::p_pot3,2.13);

        CkDec_pp.Update();
        CkDec_pL.Update();
        CkDec_LL.Update();
        if(!ExcludeXi) CkDec_pXim.Update();
        fitter->GoBabyGo();
        GraphFile->cd();
        TGraph FitResult_pL_LO;
        FitResult_pL_LO.SetName(TString::Format("FitResult_pL_LO_%u",uIter));
        fitter->GetFitGraph(1, FitResult_pL_LO);

        //change the pLambda to NLO (will have an effect only if fitting with Lednicky)
        Ck_pL->SetPotPar(0,2.91);
        Ck_pL->SetPotPar(1,2.78);
        Ck_pL->SetPotPar(2,1.54);
        Ck_pL->SetPotPar(3,2.72);
        Ck_pL->SetSourcePar(0,vSameRad==true?fitter->GetParameter("pp",DLM_Fitter1::p_sor0):fitter->GetParameter("pLambda",DLM_Fitter1::p_sor0));
        Ck_pL->Update();
        fitter->FixParameter("pLambda",DLM_Fitter1::p_pot0,2.91);
        fitter->FixParameter("pLambda",DLM_Fitter1::p_pot1,2.78);
        fitter->FixParameter("pLambda",DLM_Fitter1::p_pot2,1.54);
        fitter->FixParameter("pLambda",DLM_Fitter1::p_pot3,2.72);

        CkDec_pp.Update();
        CkDec_pL.Update();
        CkDec_LL.Update();
        if(!ExcludeXi) CkDec_pXim.Update();
        fitter->GoBabyGo();

        GraphFile->cd();
        TGraph FitResult_pL_NLO;
        FitResult_pL_NLO.SetName(TString::Format("FitResult_pL_NLO_%u",uIter));
        fitter->GetFitGraph(1, FitResult_pL_NLO);

        //we get the STAR result
        fitter->FixParameter("LambdaLambda",DLM_Fitter1::p_pot0,-1.1);
        fitter->FixParameter("LambdaLambda",DLM_Fitter1::p_pot1,8.52);
        CkDec_pp.Update();
        CkDec_pL.Update();
        CkDec_LL.Update();
        if(!ExcludeXi) CkDec_pXim.Update();
        fitter->GoBabyGo();

        GraphFile->cd();
        TGraph FitResult_LL_STAR;
        FitResult_LL_STAR.SetName(TString::Format("FitResult_LL_STAR_%u",uIter));
        fitter->GetFitGraph(2, FitResult_LL_STAR);

        //remove the pXi strong interaction
        if(!ExcludeXi){
            printf("AB_pXim[0] = %.2f\n", AB_pXim.GetCorrFun(0));
            AB_pXim.RemoveShortRangePotential(0,0);
            AB_pXim.RemoveShortRangePotential(1,0);
            AB_pXim.RemoveShortRangePotential(2,0);
            AB_pXim.RemoveShortRangePotential(3,0);
            AB_pXim.KillTheCat(CATS::kPotentialChanged);
            //AB_pXim.KillTheCat(CATS::kAllChanged);
            printf("NEW AB_pXim[0] = %.2f\n", AB_pXim.GetCorrFun(0));
        }

        CkDec_pp.Update(true);
        CkDec_pL.Update(true);
        CkDec_LL.Update(true);
        if(!ExcludeXi){
            CkDec_pXim.Update(true);
        }

        GraphFile->cd();
        TGraph FitResult_pXim_COULOMB;
        FitResult_pXim_COULOMB.SetName(TString::Format("FitResult_pXim_COULOMB_%u",uIter));
        if(!ExcludeXi) fitter->GetFitGraph(3, FitResult_pXim_COULOMB);

        FitResult_pL_LO.Write();
        FitResult_pL_NLO.Write();
        FitResult_LL_STAR.Write();
        FitResult_pXim_COULOMB.Write();

        double Chi2_pL_LO=0;
        unsigned EffNumBins_pL_LO=0;
        for(unsigned uBin=0; uBin<NumMomBins_pL; uBin++){

            double mom = AB_pL.GetMomentum(uBin);
            double dataY;
            double dataErr;
            double theoryX;
            double theoryY;

            if(mom>FemtoRegion_pL[vFemReg_pL][1]) continue;

            FitResult_pL_LO.GetPoint(uBin,theoryX,theoryY);
            //if(mom!=theoryX) printf("PROBLEM!\n");
            dataY = OliHisto_pL->GetBinContent(uBin+1);
            dataErr = OliHisto_pL->GetBinError(uBin+1);

            Chi2_pL_LO += (dataY-theoryY)*(dataY-theoryY)/(dataErr*dataErr);
            EffNumBins_pL_LO++;
        }

        double Chi2_pL_NLO=0;
        unsigned EffNumBins_pL_NLO=0;
        for(unsigned uBin=0; uBin<NumMomBins_pL; uBin++){

            double mom = AB_pL.GetMomentum(uBin);
            double dataY;
            double dataErr;
            double theoryX;
            double theoryY;

            if(mom>FemtoRegion_pL[vFemReg_pL][1]) continue;

            FitResult_pL_NLO.GetPoint(uBin,theoryX,theoryY);
            //if(mom!=theoryX) printf("PROBLEM!\n");
            dataY = OliHisto_pL->GetBinContent(uBin+1);
            dataErr = OliHisto_pL->GetBinError(uBin+1);

            Chi2_pL_NLO += (dataY-theoryY)*(dataY-theoryY)/(dataErr*dataErr);
            EffNumBins_pL_NLO++;
        }

        double Chi2_LL_STAR=0;
        unsigned EffNumBins_LL_STAR=0;
        for(unsigned uBin=0; uBin<NumMomBins_LL; uBin++){

            double mom = Ck_LL->GetBinCenter(0,uBin);
            double dataY;
            double dataErr;
            double theoryX;
            double theoryY;

            if(mom>FemtoRegion_LL[vFemReg_LL][1]) continue;

            FitResult_LL_STAR.GetPoint(uBin,theoryX,theoryY);
            //if(mom!=theoryX) printf("PROBLEM!\n");
            dataY = OliHisto_LL->GetBinContent(uBin+1);
            dataErr = OliHisto_LL->GetBinError(uBin+1);

            Chi2_LL_STAR += (dataY-theoryY)*(dataY-theoryY)/(dataErr*dataErr);
            EffNumBins_LL_STAR++;
        }

        double Chi2_pXim_COULOMB=0;
        unsigned EffNumBins_pXim_COULOMB=0;
        for(unsigned uBin=0; uBin<NumMomBins_pXim; uBin++){

            double mom = AB_pXim.GetMomentum(uBin);
            double dataY;
            double dataErr;
            double theoryX;
            double theoryY;

            if(mom>FemtoRegion_pXim[vFemReg_pXim][1]) continue;

            FitResult_pXim_COULOMB.GetPoint(uBin,theoryX,theoryY);
            //if(mom!=theoryX) printf("PROBLEM!\n");
            dataY = OliHisto_pXim->GetBinContent(uBin+1);
            dataErr = OliHisto_pXim->GetBinError(uBin+1);

            Chi2_pXim_COULOMB += (dataY-theoryY)*(dataY-theoryY)/(dataErr*dataErr);
            EffNumBins_pXim_COULOMB++;
        }

        if(FAST_PLOT){
            TPaveText* info1 = new TPaveText(0.45,0.65,0.9,0.95, "blNDC");//lbrt
            info1->SetName("info1");
            info1->SetBorderSize(1);
            info1->SetTextSize(0.04);
            info1->SetFillColor(kWhite);
            info1->SetTextFont(22);
            TString SOURCE_NAME;
            if(vSource==0) SOURCE_NAME="Gauss";
            if(vSource==1) SOURCE_NAME="Cauchy";
            if(vSource==2) SOURCE_NAME="DoubleGauss";
            if(vSource==3) SOURCE_NAME="EPOS";
            if(vSource==4) SOURCE_NAME="Gauss+Cauchy";
            if(vSource==5) SOURCE_NAME="3D Gauss";
            if(vSource==9) SOURCE_NAME="Levy";
            if(vSource==10) SOURCE_NAME="CleverLevy";
            if(vSource==6){
                info1->AddText("GaussCore+Reso (fixed)");
            }
            else if(vSource==8){
                info1->AddText("DoubleGaussCore+Reso (fixed)");
            }
            else{
                info1->AddText(TString::Format("R(%s)=%.3f#pm%.3f",SOURCE_NAME.Data(),
                        RadiusResult,RadiusError));
            }
            if(vSource==2 || vSource==4){
                info1->AddText(TString::Format("R2=%.3f#pm%.3f",
                        RadiusResult2,RadiusError2));
                info1->AddText(TString::Format("w=%.3f#pm%.3f",
                        RadiusResult3,RadiusError3));
            }
            else if(vSource==5){
                info1->AddText(TString::Format("R_{side}=%.3f#pm%.3f",
                        RadiusResult2,RadiusError2));
                info1->AddText(TString::Format("R_{long}=%.3f#pm%.3f",
                        RadiusResult3,RadiusError3));
            }
            else if(vSource==9||vSource==10){
                info1->AddText(TString::Format("#alpha=%.3f#pm%.3f",
                        RadiusResult2,RadiusError2));
            }
            info1->AddText(TString::Format("C_{l}=%.3f#pm%.3f",
                            ClResult_pp,ClError_pp));
            info1->AddText(TString::Format("Global #chi^{2}_{ndf}=%.0f/%u=%.2f, p_{val}=%.3f",
                        Chi2,NDF,Chi2/double(NDF),TMath::Prob(Chi2,round(NDF))));
            info1->AddText(TString::Format("Local #chi^{2}/ndf=%.0f/%u=%.2f, p_{val}=%.3f",
                        Chi2_pp,EffNumBins_pp,Chi2_pp/double(EffNumBins_pp),TMath::Prob(Chi2_pp,round(EffNumBins_pp))));

            TPaveText* info2 = new TPaveText(0.35,0.63,0.9,0.95, "blNDC");//lbrt
            info2->SetName("info2");
            info2->SetBorderSize(1);
            info2->SetTextSize(0.04);
            info2->SetFillColor(kWhite);
            info2->SetTextFont(22);
            if(vSource==0) SOURCE_NAME="Gauss";
            if(vSource==1) SOURCE_NAME="Cauchy";
            if(vSource==2) SOURCE_NAME="DoubleGauss";
            if(vSource==3) SOURCE_NAME="EPOS";
            if(vSource==4) SOURCE_NAME="Gauss+Cauchy";
            if(vSource==5) SOURCE_NAME="3D Gauss";
            if(vSource==9) SOURCE_NAME="Levy";
            if(vSource==10) SOURCE_NAME="CleverLevy";
            if(vSource==6){
                info2->AddText("GaussCore+Reso (fixed)");
            }
            else if(vSource==8){
                info2->AddText("DoubleGaussCore+Reso (fixed)");
            }
            else{
                info2->AddText(TString::Format("R(%s)=%.3f#pm%.3f",SOURCE_NAME.Data(),
                        RadiusResult_pL,RadiusError_pL));
            }
            if((vSource==2 || vSource==4) && vSameRad==false){
                info2->AddText(TString::Format("R2=%.3f#pm%.3f",
                        RadiusResult2_pL,RadiusError2_pL));
                info2->AddText(TString::Format("w=%.3f#pm%.3f",
                        RadiusResult3_pL,RadiusError3_pL));
            }
            else if(vSource==5 && vSameRad==false ){
                info2->AddText(TString::Format("R_{side}=%.3f#pm%.3f",
                        RadiusResult2_pL,RadiusError2_pL));
                info2->AddText(TString::Format("R_{long}=%.3f#pm%.3f",
                        RadiusResult3_pL,RadiusError3_pL));
            }
            else if((vSource==9||vSource==10) && vSameRad==false){
                info2->AddText(TString::Format("#alpha=%.3f#pm%.3f",
                        RadiusResult2_pL,RadiusError2_pL));
            }
            info2->AddText(TString::Format("C_{l}=%.3f#pm%.3f",
                            ClResult_pL,ClError_pL));
            info2->AddText(TString::Format("Local #chi^{2}/ndf=%.0f/%u=%.2f, p_{val}=%.3f, n#sigma=%.2f",
                                           Chi2_pL,EffNumBins_pL,Chi2_pL/double(EffNumBins_pL),TMath::Prob(Chi2_pL,round(EffNumBins_pL)),sqrt(2)*TMath::ErfcInverse(TMath::Prob(Chi2_pL,round(EffNumBins_pL)))));
            info2->AddText(TString::Format("LO #chi^{2}/ndf=%.0f/%u=%.2f, p_{val}=%.3f, n#sigma=%.2f",
                                           Chi2_pL_LO,EffNumBins_pL,Chi2_pL_LO/double(EffNumBins_pL),TMath::Prob(Chi2_pL_LO,round(EffNumBins_pL)),sqrt(2)*TMath::ErfcInverse(TMath::Prob(Chi2_pL_LO,round(EffNumBins_pL)))));
            info2->AddText(TString::Format("NLO #chi^{2}/ndf=%.0f/%u=%.2f, p_{val}=%.3f, n#sigma=%.2f",
                                           Chi2_pL_NLO,EffNumBins_pL,Chi2_pL_NLO/double(EffNumBins_pL),TMath::Prob(Chi2_pL_NLO,round(EffNumBins_pL)),sqrt(2)*TMath::ErfcInverse(TMath::Prob(Chi2_pL_NLO,round(EffNumBins_pL)))));
            TPaveText* info3 = new TPaveText(0.35,0.63,0.9,0.95, "blNDC");//lbrt
            info3->SetName("info3");
            info3->SetBorderSize(1);
            info3->SetTextSize(0.04);
            info3->SetFillColor(kWhite);
            info3->SetTextFont(22);
            SOURCE_NAME="Gauss";
            if(vSource==0)
            info3->AddText(TString::Format("R(%s)=%.3f#pm%.3f",SOURCE_NAME.Data(),
                   RadiusResult_LL,RadiusError_LL));
            else info3->AddText("R(Gauss) fixed to 1.2");
            info3->AddText(TString::Format("a_{0}=%.2f#pm%.2f",a0_LL,a0_err_LL));
            info3->AddText(TString::Format("r_{eff}=%.2f#pm%.2f",reff_LL,reff_err_LL));
            info3->AddText(TString::Format("C_{l}=%.2f#pm%.2f",
                            ClResult_LL,ClError_LL));
            info3->AddText(TString::Format("Local #chi^{2}/ndf=%.0f/%u=%.2f, p_{val}=%.3f, n#sigma=%.2f",
                                           Chi2_LL,EffNumBins_LL,Chi2_LL/double(EffNumBins_LL),TMath::Prob(Chi2_LL,round(EffNumBins_LL)),sqrt(2)*TMath::ErfcInverse(TMath::Prob(Chi2_LL,round(EffNumBins_LL)))));
            info3->AddText(TString::Format("STAR #chi^{2}/ndf=%.0f/%u=%.2f, p_{val}=%.3f, n#sigma=%.2f",
                                           Chi2_LL_STAR,EffNumBins_LL,Chi2_LL_STAR/double(EffNumBins_LL),TMath::Prob(Chi2_LL_STAR,round(EffNumBins_LL)),sqrt(2)*TMath::ErfcInverse(TMath::Prob(Chi2_LL_STAR,round(EffNumBins_LL)))));

            TPaveText* info4 = new TPaveText(0.35,0.68,0.9,0.95, "blNDC");//lbrt
            info4->SetName("info4");
            info4->SetBorderSize(1);
            info4->SetTextSize(0.04);
            info4->SetFillColor(kWhite);
            info4->SetTextFont(22);
            if(vSource==0) SOURCE_NAME="Gauss";
            if(vSource==1) SOURCE_NAME="Cauchy";
            if(vSource==2) SOURCE_NAME="DoubleGauss";
            if(vSource==3) SOURCE_NAME="EPOS";
            if(vSource==4) SOURCE_NAME="Gauss+Cauchy";
            if(vSource==5) SOURCE_NAME="3D Gauss";
            if(vSource==9) SOURCE_NAME="Levy";
            if(vSource==10) SOURCE_NAME="CleverLevy";
            if(ExcludeXi==false){
                if(vSource==6){
                    info4->AddText("GaussCore+Reso (fixed)");
                }
                else if(vSource==8){
                    info4->AddText("DoubleGaussCore+Reso (fixed)");
                }
                else{
                    info4->AddText(TString::Format("R(%s)=%.3f#pm%.3f",SOURCE_NAME.Data(),
                            RadiusResult_pXim,RadiusError_pXim));
                }
                if((vSource==2 || vSource==4) && vSameRad==false){
                    info4->AddText(TString::Format("R2=%.3f#pm%.3f",
                            fitter->GetParameter("pXim",DLM_Fitter1::p_sor1),fitter->GetParError("pXim",DLM_Fitter1::p_sor1)));
                    info4->AddText(TString::Format("w=%.3f#pm%.3f",
                            fitter->GetParameter("pXim",DLM_Fitter1::p_sor2),fitter->GetParError("pXim",DLM_Fitter1::p_sor2)));
                }
                else if(vSource==5 && vSameRad==false ){
                    info4->AddText(TString::Format("R_{side}=%.3f#pm%.3f",
                            RadiusResult2_pXim,RadiusError2_pXim));
                    info4->AddText(TString::Format("R_{long}=%.3f#pm%.3f",
                            RadiusResult3_pXim,RadiusError3_pXim));
                }
                else if((vSource==9||vSource==10) && vSameRad==false){
                    info4->AddText(TString::Format("#alpha=%.3f#pm%.3f",
                            RadiusResult2_pXim,RadiusError2_pXim));
                }
                info4->AddText(TString::Format("C_{l}=%.3f#pm%.3f",
                                ClResult_pXim,ClError_pXim));
                info4->AddText(TString::Format("Local #chi^{2}/ndf=%.0f/%u=%.2f, p_{val}=%.3f, n#sigma=%.2f",
                                               Chi2_pXim,EffNumBins_pXim,Chi2_pXim/double(EffNumBins_pXim),TMath::Prob(Chi2_pXim,round(EffNumBins_pXim)),sqrt(2)*TMath::ErfcInverse(TMath::Prob(Chi2_pXim,round(EffNumBins_pXim)))));
                info4->AddText(TString::Format("Coulomb #chi^{2}/ndf=%.0f/%u=%.2f, p_{val}=%.3f, n#sigma=%.2f",
                                               Chi2_pXim_COULOMB,EffNumBins_pXim,Chi2_pXim_COULOMB/double(EffNumBins_pXim),TMath::Prob(Chi2_pXim_COULOMB,round(EffNumBins_pXim)),sqrt(2)*TMath::ErfcInverse(TMath::Prob(Chi2_pXim_COULOMB,round(EffNumBins_pXim)))));
            }


            double Yoffset = 1.2;
            TH1F* hAxis_pp = new TH1F("hAxis_pp", "hAxis_pp", 600, 0, 600);
            hAxis_pp->SetStats(false);
            hAxis_pp->SetTitle("");
            hAxis_pp->GetXaxis()->SetLabelSize(0.065);
            hAxis_pp->GetXaxis()->CenterTitle();
            hAxis_pp->GetXaxis()->SetTitleOffset(1.35);
            hAxis_pp->GetXaxis()->SetLabelOffset(0.02);
            hAxis_pp->GetXaxis()->SetTitleSize(0.075);
            hAxis_pp->GetYaxis()->SetLabelSize(0.065);
            hAxis_pp->GetYaxis()->CenterTitle();
            hAxis_pp->GetYaxis()->SetTitleOffset(Yoffset);
            hAxis_pp->GetYaxis()->SetTitleSize(0.075);
            hAxis_pp->GetXaxis()->SetRangeUser(0,kMax_pp);
            if(DataPeriod=="pp13TeV") hAxis_pp->GetYaxis()->SetRangeUser(0.5, 4);//pp
            else hAxis_pp->GetYaxis()->SetRangeUser(0.5, 3);//pPb

            TH1F* hAxis_pL = new TH1F("hAxis_pL", "hAxis_pL", 600, 0, 600);
            hAxis_pL->SetStats(false);
            hAxis_pL->SetTitle("");
            hAxis_pL->GetXaxis()->SetLabelSize(0.065);
            hAxis_pL->GetXaxis()->CenterTitle();
            hAxis_pL->GetXaxis()->SetTitleOffset(1.35);
            hAxis_pL->GetXaxis()->SetLabelOffset(0.02);
            hAxis_pL->GetXaxis()->SetTitleSize(0.075);
            hAxis_pL->GetYaxis()->SetLabelSize(0.065);
            hAxis_pL->GetYaxis()->CenterTitle();
            hAxis_pL->GetYaxis()->SetTitleOffset(Yoffset);
            hAxis_pL->GetYaxis()->SetTitleSize(0.075);
            hAxis_pL->GetXaxis()->SetRangeUser(0,kMax_pL);
            if(DataPeriod=="pp13TeV") hAxis_pL->GetYaxis()->SetRangeUser(0.8, 2.2);//pp
            else hAxis_pL->GetYaxis()->SetRangeUser(0.8, 2.0);//pPb

            TH1F* hAxis_LL = new TH1F("hAxis_LL", "hAxis_LL", 600, 0, 600);
            hAxis_LL->SetStats(false);
            hAxis_LL->SetTitle("");
            hAxis_LL->GetXaxis()->SetLabelSize(0.065);
            //hAxis_LL->GetXaxis()->SetTitle("k (MeV)");
            hAxis_LL->GetXaxis()->CenterTitle();
            hAxis_LL->GetXaxis()->SetTitleOffset(1.35);
            hAxis_LL->GetXaxis()->SetLabelOffset(0.02);
            hAxis_LL->GetXaxis()->SetTitleSize(0.075);
            hAxis_LL->GetYaxis()->SetLabelSize(0.065);
            //hAxis_LL->GetYaxis()->SetTitle("C(k)");
            hAxis_LL->GetYaxis()->CenterTitle();
            hAxis_LL->GetYaxis()->SetTitleOffset(Yoffset);
            hAxis_LL->GetYaxis()->SetTitleSize(0.075);
            //hAxis_LL->GetXaxis()->SetRangeUser(0,kMax_LL);
            hAxis_LL->GetYaxis()->SetRangeUser(0.5, 1.4);

            TH1F* hAxis_pXim = new TH1F("hAxis_pXim", "hAxis_pXim", 600, 0, 600);
            hAxis_pXim->SetStats(false);
            hAxis_pXim->SetTitle("");
            hAxis_pXim->GetXaxis()->SetLabelSize(0.065);
            hAxis_pXim->GetXaxis()->CenterTitle();
            hAxis_pXim->GetXaxis()->SetTitleOffset(1.35);
            hAxis_pXim->GetXaxis()->SetLabelOffset(0.02);
            hAxis_pXim->GetXaxis()->SetTitleSize(0.075);
            hAxis_pXim->GetYaxis()->SetLabelSize(0.065);
            hAxis_pXim->GetYaxis()->CenterTitle();
            hAxis_pXim->GetYaxis()->SetTitleOffset(Yoffset);
            hAxis_pXim->GetYaxis()->SetTitleSize(0.075);
            hAxis_pXim->GetXaxis()->SetRangeUser(0,kMax_pXim);
            if(DataPeriod=="pp13TeV") hAxis_pXim->GetYaxis()->SetRangeUser(0.7, 3.5);//pp
            else hAxis_pXim->GetYaxis()->SetRangeUser(0.7, 4.5);//pPb

            TCanvas* cfast = new TCanvas("cfast","cfast",1);
            cfast->cd(0); cfast->SetCanvasSize(1920, 1280); cfast->SetMargin(0.15,0.05,0.2,0.05);//lrbt
            cfast->Divide(2,2);

            cfast->cd(1);
            OliHisto_pp->SetStats(false);
            OliHisto_pp->SetTitle("pp");
            OliHisto_pp->SetLineWidth(2);
            OliHisto_pp->SetLineColor(kBlack);
            FitResult_pp.SetLineWidth(2);
            FitResult_pp.SetLineColor(kRed);
            FitResult_pp.SetMarkerStyle(24);
            FitResult_pp.SetMarkerColor(kRed);
            FitResult_pp.SetMarkerSize(1);

            hAxis_pp->Draw("axis");
            OliHisto_pp->Draw("same");
            FitResult_pp.Draw("CP,same");
            info1->Draw("same");

            cfast->cd(2);
            OliHisto_pL->SetTitle("p#Lambda");
            OliHisto_pL->SetLineWidth(2);
            OliHisto_pL->SetLineColor(kBlack);
            FitResult_pL.SetLineWidth(2);
            FitResult_pL.SetLineColor(kRed);
            FitResult_pL.SetMarkerStyle(24);
            FitResult_pL.SetMarkerColor(kRed);
            FitResult_pL.SetMarkerSize(1);

            hAxis_pL->Draw("axis");
            OliHisto_pL->Draw("same");
            FitResult_pL.Draw("CP,same");
            info2->Draw("same");

            cfast->cd(3);
            OliHisto_LL->SetTitle("#Lambda#Lambda");
            OliHisto_LL->SetLineWidth(2);
            OliHisto_LL->SetLineColor(kBlack);
            FitResult_LL.SetLineWidth(2);
            FitResult_LL.SetLineColor(kRed);
            FitResult_LL.SetMarkerStyle(24);
            FitResult_LL.SetMarkerColor(kRed);
            FitResult_LL.SetMarkerSize(1);

            hAxis_LL->Draw("axis");
            OliHisto_LL->Draw("same");
            FitResult_LL.Draw("CP,same");
            info3->Draw("same");

            if(ExcludeXi==false){
                cfast->cd(4);
                OliHisto_pXim->SetTitle("p#Xi^{#minus}");
                OliHisto_pXim->SetLineWidth(2);
                OliHisto_pXim->SetLineColor(kBlack);
                FitResult_pXim.SetLineWidth(2);
                FitResult_pXim.SetLineColor(kRed);
                FitResult_pXim.SetMarkerStyle(24);
                FitResult_pXim.SetMarkerColor(kRed);
                FitResult_pXim.SetMarkerSize(1);

                hAxis_pXim->Draw("axis");
                OliHisto_pXim->Draw("same");
                FitResult_pXim.Draw("CP,same");
                info4->Draw("same");
            }

            cfast->Write();
            cfast->SaveAs(TString::Format("%scfast.png", OutputDir.Data()));

            AB_pXim.SetShortRangePotential(0,0,fDlmPot,c_pXimPotParsI0S0);
            AB_pXim.SetShortRangePotential(1,0,fDlmPot,c_pXimPotParsI0S1);
            AB_pXim.SetShortRangePotential(2,0,fDlmPot,c_pXimPotParsI1S0);
            AB_pXim.SetShortRangePotential(3,0,fDlmPot,c_pXimPotParsI1S1);
            AB_pXim.KillTheCat(CATS::kPotentialChanged);

            if(FULL_PLOT){
                TGraph gFeed_To_pp_from_pL;
                gFeed_To_pp_from_pL.SetName("gFeed_To_pp_from_pL");
                TGraph gFeed_To_pL_from_pSigma0;
                gFeed_To_pL_from_pSigma0.SetName("gFeed_To_pL_from_pSigma0");
                TGraph gFeed_To_pL_from_pXim;
                gFeed_To_pL_from_pXim.SetName("gFeed_To_pL_from_pXim");
                TGraph gFeed_To_pXim_from_pXim1530;
                gFeed_To_pXim_from_pXim1530.SetName("gFeed_To_pXim_from_pXim1530");

                TGraph g_pp;
                g_pp.SetName("g_pp");
                TGraph g_pL;
                g_pL.SetName("g_pL");
                TGraph g_LL;
                g_LL.SetName("g_LL");
                TGraph g_pXim;
                g_pXim.SetName("g_pXim");
                TGraph g_pXim1530;
                g_pXim1530.SetName("g_pXim1530");
                TGraph g_pSigma0;
                g_pSigma0.SetName("g_pSigma0");

                TGraph gSmear_pp;
                gSmear_pp.SetName("gSmear_pp");
                TGraph gSmear_pL;
                gSmear_pL.SetName("gSmear_pL");
                TGraph gSmear_LL;
                gSmear_LL.SetName("gSmear_LL");
                TGraph gSmear_pXim;
                gSmear_pXim.SetName("gSmear_pXim");
                TGraph gSmear_pXim1530;
                gSmear_pXim1530.SetName("gSmear_pXim1530");

            }//FULL_PLOT

            delete info1;
            delete info2;
            delete info3;
            delete info4;
            delete cfast;

        }//FAST_PLOT

        delete fitter;
        delete Ck_pp;
        delete Ck_pL;
        delete Ck_LL;
        delete Ck_pSigma0;
        delete Ck_pXim;
        delete Ck_pXim1530;

        if(OliFile_pp) {delete OliFile_pp; OliFile_pp=NULL;}
        if(OliFile_pL) {delete OliFile_pL; OliFile_pL=NULL;}
        if(OliFile_LL) {delete OliFile_LL; OliFile_LL=NULL;}
        if(OliFile_pXim) {delete OliFile_pXim; OliFile_pXim=NULL;}
        if(SystErrFile_pp) {delete SystErrFile_pp; SystErrFile_pp=NULL;}
        if(SystErrFile_pL) {delete SystErrFile_pL; SystErrFile_pL=NULL;}
        if(SystErrFile_LL) {delete SystErrFile_LL; SystErrFile_LL=NULL;}
        if(SystErrFile_pXim) {delete SystErrFile_pXim; SystErrFile_pXim=NULL;}

        delete GraphFile;
    }//END OF THE BIG FOR LOOP OVER ITER

    OutFile->cd();
    ntResult->Write();

    delete hRes_pp_pL_MeV;
    delete hRes_pL_pSigma0_MeV;
    delete hRes_pL_pXim_MeV;
    delete hRes_pXim_pXim1530_MeV;
    if(hSigma_pp_MeV){delete hSigma_pp_MeV; hSigma_pp_MeV=NULL;}
    if(hSigma_pL_MeV){delete hSigma_pL_MeV; hSigma_pL_MeV=NULL;}
    if(hSigma_LL_MeV){delete hSigma_LL_MeV; hSigma_LL_MeV=NULL;}
    if(hSigma_pXim_MeV){delete hSigma_pXim_MeV; hSigma_pXim_MeV=NULL;}
    //delete hSigma_pXim_MeV;
    delete ntResult;
    delete OutFile;
    delete FileRes;
    delete FileSigma;

    for(unsigned uVar=0; uVar<3; uVar++){
        delete [] Purities_p[uVar];
        delete [] Fraction_p[uVar];
        delete [] Purities_L[uVar];
        delete [] Fraction_L[uVar];
    }
    delete [] Purities_p;
    delete [] Fraction_p;
    delete [] Purities_L;
    delete [] Fraction_L;

    if(ExternalWF) CleanUpWfHisto(AB_pL,ExternalWF);
}


void CALL_BERNIE_AND_VALE(){
    //printf("Calling for Берни и Вале\n");
    //printf("Cock\n");
    RUN2_main(1,1,0);
}

