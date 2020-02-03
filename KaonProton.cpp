#include <stdlib.h>

#include "KaonProton.h"

#include "CATS.h"
#include "DLM_Potentials.h"
#include "DLM_Source.h"
#include "DLM_DrawingTools.h"
#include "DLM_CkDecomposition.h"
#include "DLM_CkModels.h"
#include "DLM_Fitters.h"
#include "DLM_WfModel.h"

#include "DLM_SubPads.h"

#include "TNtuple.h"
#include "TString.h"
#include "TGraph.h"
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TRandom3.h"

//without the coulomb interaction!
void AnitKaonProton_PhaseShiftsAndWaveFunctions(){
/*
    const double kMin=5;
    const double kMax=305;
    const unsigned NumMomBins=30;

    const TString OutputFolder = "./OutputKaonProton/AKP_PS_WF/070518/";


    //0 Haide
    //1 Tetsuo pots
    //2 Tetsuo WF
    const unsigned NumScenarios = 3;
    TString ScenarioName[NumScenarios];

    ScenarioName[0] = "Haidenbauer";
    ScenarioName[1] = "TetsuoPot";
    ScenarioName[2] = "TetsuoWf";

    double SourcePars_Kp[4] = {0,0,0,1.19};

    CATS MyCat;
    MyCat.SetUseAnalyticSource(true);
    MyCat.SetAnaSource(GaussSource, SourcePars_Kp);

    MyCat.SetExcludeFailedBins(false);
    MyCat.SetMomBins(NumMomBins,kMin,kMax);

    MyCat.SetNumChannels(3);
    MyCat.SetNumPW(0,1);
    MyCat.SetNumPW(1,1);
    MyCat.SetNumPW(2,1);
    //change to double?
    MyCat.SetSpin(0,1);//I=0
    MyCat.SetSpin(1,1);//I=1
    MyCat.SetSpin(2,1);//the no SI channel

    //this is deliberately, since we would like to study only the SI
    MyCat.SetQ1Q2(1);
    MyCat.SetPdgId(2212, 222);//check the pdg id of the kaon
    const double Mass_p = 938.272;
    const double Mass_Km = 493.677;
    MyCat.SetRedMass( (Mass_p*Mass_Km)/(Mass_p+Mass_Km) );

    complex<double>**** H_WaveFunctionU; double*** H_PhaseShifts; double* H_RadBins; unsigned H_NumRadBins=0;
    InitHaidenbauerKaonMinus("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Haidenbauer/KaonMinus10MeV/",MyCat,&H_WaveFunctionU,&H_PhaseShifts,&H_RadBins,H_NumRadBins);

    complex<double>**** T_WaveFunctionU; double*** T_PhaseShifts; double* T_RadBins; unsigned T_NumRadBins=0;
    InitTetsuoKaonMinus("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Tetsuo/KminusProton/",MyCat,&T_WaveFunctionU,&T_PhaseShifts,&T_RadBins,T_NumRadBins, 1);

    double pKmPotParsI0[10]={0,0,pKm_Tetsuo,0,0,1,1,1,0,1};
    double pKmPotParsI1[10]={0,0,pKm_Tetsuo,0,1,1,1,1,0,1};

    const double MaxRadius = 32;
    const unsigned NumRadiusBins = 384;
    TH1F*** hRadiusI0 = new TH1F** [NumScenarios];
    TH1F*** hRadiusI1 = new TH1F** [NumScenarios];

    TH1F** hPsI0 = new TH1F* [NumScenarios];
    TH1F** hPsI1 = new TH1F* [NumScenarios];

    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        hRadiusI0[uSce] = new TH1F* [NumMomBins];
        hRadiusI1[uSce] = new TH1F* [NumMomBins];

        hPsI0[uSce] = new TH1F( TString::Format("hPsI0_%s",ScenarioName[uSce].Data()),
                                TString::Format("hPsI0_%s",ScenarioName[uSce].Data()),
                                NumMomBins, kMin, kMax);
        hPsI1[uSce] = new TH1F( TString::Format("hPsI1_%s",ScenarioName[uSce].Data()),
                                TString::Format("hPsI1_%s",ScenarioName[uSce].Data()),
                                NumMomBins, kMin, kMax);

        //reset the external function
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            hRadiusI0[uSce][uMomBin] = new TH1F(    TString::Format("hRadiusI0_%s_%.0f",ScenarioName[uSce].Data(),MyCat.GetMomentum(uMomBin)),
                                                    TString::Format("hRadiusI0_%s_%.0f",ScenarioName[uSce].Data(),MyCat.GetMomentum(uMomBin)),
                                                    NumRadiusBins, 0, MaxRadius);
            hRadiusI1[uSce][uMomBin] = new TH1F(    TString::Format("hRadiusI1_%s_%.0f",ScenarioName[uSce].Data(),MyCat.GetMomentum(uMomBin)),
                                                    TString::Format("hRadiusI1_%s_%.0f",ScenarioName[uSce].Data(),MyCat.GetMomentum(uMomBin)),
                                                    NumRadiusBins, 0, MaxRadius);
            MyCat.UseExternalWaveFunction(uMomBin,0,0);//I=0
            MyCat.UseExternalWaveFunction(uMomBin,1,0);//I=1
        }
        switch(uSce){
        case 0 :
            for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
                MyCat.UseExternalWaveFunction(uMomBin,0,0,H_WaveFunctionU[uMomBin][0][0],H_NumRadBins,H_RadBins,H_PhaseShifts[uMomBin][0][0]);//I=0
                MyCat.UseExternalWaveFunction(uMomBin,1,0,H_WaveFunctionU[uMomBin][2][0],H_NumRadBins,H_RadBins,H_PhaseShifts[uMomBin][2][0]);//I=1
            }
            MyCat.SetChannelWeight(0, 1./2.);
            MyCat.SetChannelWeight(1, 1./2.);
            break;
        case 1 :
            MyCat.SetShortRangePotential(0,0,fDlmPot,pKmPotParsI0);
            MyCat.SetShortRangePotential(1,0,fDlmPot,pKmPotParsI1);
            MyCat.SetChannelWeight(0, 1./2.);
            MyCat.SetChannelWeight(1, 1./2.);
            break;
        case 2  :
            for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
                MyCat.UseExternalWaveFunction(uMomBin,0,0,T_WaveFunctionU[uMomBin][0][0],T_NumRadBins,T_RadBins,T_PhaseShifts[uMomBin][0][0]);
            }
            MyCat.SetChannelWeight(0, 1);
            MyCat.SetChannelWeight(1, 0);
            break;
        default :
            break;
        }

        MyCat.KillTheCat();

        double Radius;
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            for(unsigned uRad=0; uRad<NumRadiusBins; uRad++){
                Radius = hRadiusI0[uSce][uMomBin]->GetBinCenter(uRad+1);
                hRadiusI0[uSce][uMomBin]->SetBinContent(uRad+1, real(MyCat.EvalRadialWaveFunction(uMomBin,0,0,Radius,false)));
                hRadiusI1[uSce][uMomBin]->SetBinContent(uRad+1, real(MyCat.EvalRadialWaveFunction(uMomBin,1,0,Radius,false)));
            }

            hPsI0[uSce]->SetBinContent(uMomBin+1, MyCat.GetPhaseShift(uMomBin,0,0));
            hPsI1[uSce]->SetBinContent(uMomBin+1, MyCat.GetPhaseShift(uMomBin,1,0));
        }
    }

    TFile* OutputFile = new TFile(TString::Format("%sOutputFile.root",OutputFolder.Data()),"recreate");
    for(unsigned uSce=0; uSce<NumScenarios; uSce++) for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++) hRadiusI0[uSce][uMomBin]->Write();
    for(unsigned uSce=0; uSce<NumScenarios; uSce++) for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++) hRadiusI1[uSce][uMomBin]->Write();
    for(unsigned uSce=0; uSce<NumScenarios; uSce++) hPsI0[uSce]->Write();
    for(unsigned uSce=0; uSce<NumScenarios; uSce++) hPsI1[uSce]->Write();


    //char* fBuffer = new char [512];
    FILE * pFileTetWf;
    pFileTetWf = fopen (TString::Format("%sKminusProton_Wf_WithCoulomb_Tetsuo.txt",OutputFolder.Data()),"w");
    fprintf (pFileTetWf, "%12s","k (MeV)");
    fprintf (pFileTetWf, "%12s","r (fm)");
    fprintf (pFileTetWf, "%12s","u (1/fm^2)");
    fprintf(pFileTetWf,"\n");
    for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
        for(unsigned uRad=0; uRad<NumRadiusBins; uRad++){
            fprintf (pFileTetWf, "%12.3e",MyCat.GetMomentum(uMomBin));
            fprintf (pFileTetWf, "%12.3e",hRadiusI0[2][uMomBin]->GetBinCenter(uRad+1));
            fprintf (pFileTetWf, "%12.3e",hRadiusI0[2][uMomBin]->GetBinContent(uRad+1));
            fprintf(pFileTetWf,"\n");
        }
    }

    //char* fBuffer = new char [512];
    FILE * pFileTetI0;
    pFileTetI0 = fopen (TString::Format("%sKminusProton_I0_WithCoulomb_Tetsuo.txt",OutputFolder.Data()),"w");
    fprintf (pFileTetI0, "%12s","k (MeV)");
    fprintf (pFileTetI0, "%12s","r (fm)");
    fprintf (pFileTetI0, "%12s","u (1/fm^2)");
    fprintf(pFileTetI0,"\n");
    for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
        for(unsigned uRad=0; uRad<NumRadiusBins; uRad++){
            fprintf (pFileTetI0, "%12.3e",MyCat.GetMomentum(uMomBin));
            fprintf (pFileTetI0, "%12.3e",hRadiusI0[1][uMomBin]->GetBinCenter(uRad+1));
            fprintf (pFileTetI0, "%12.3e",hRadiusI0[1][uMomBin]->GetBinContent(uRad+1));
            fprintf(pFileTetI0,"\n");
        }
    }
    FILE * pFileTetI1;
    pFileTetI1 = fopen (TString::Format("%sKminusProton_I1_WithCoulomb_Tetsuo.txt",OutputFolder.Data()),"w");
    fprintf (pFileTetI1, "%12s","k (MeV)");
    fprintf (pFileTetI1, "%12s","r (fm)");
    fprintf (pFileTetI1, "%12s","u (1/fm^2)");
    fprintf(pFileTetI1,"\n");
    for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
        for(unsigned uRad=0; uRad<NumRadiusBins; uRad++){
            fprintf (pFileTetI1, "%12.3e",MyCat.GetMomentum(uMomBin));
            fprintf (pFileTetI1, "%12.3e",hRadiusI1[1][uMomBin]->GetBinCenter(uRad+1));
            fprintf (pFileTetI1, "%12.3e",hRadiusI1[1][uMomBin]->GetBinContent(uRad+1));
            fprintf(pFileTetI1,"\n");
        }
    }

    FILE * pFileHaideI0;
    pFileHaideI0 = fopen (TString::Format("%sKminusProton_I0_WithCoulomb_Haidenbauer.txt",OutputFolder.Data()),"w");
    fprintf (pFileHaideI0, "%12s","k (MeV)");
    fprintf (pFileHaideI0, "%12s","r (fm)");
    fprintf (pFileHaideI0, "%12s","u (1/fm^2)");
    fprintf(pFileHaideI0,"\n");
    for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
        for(unsigned uRad=0; uRad<NumRadiusBins; uRad++){
            fprintf (pFileHaideI0, "%12.3e",MyCat.GetMomentum(uMomBin));
            fprintf (pFileHaideI0, "%12.3e",hRadiusI0[0][uMomBin]->GetBinCenter(uRad+1));
            fprintf (pFileHaideI0, "%12.3e",hRadiusI0[0][uMomBin]->GetBinContent(uRad+1));
            fprintf(pFileHaideI0,"\n");
        }
    }
    FILE * pFileHaideI1;
    pFileHaideI1 = fopen (TString::Format("%sKminusProton_I1_WithCoulomb_Haidenbauer.txt",OutputFolder.Data()),"w");
    fprintf (pFileHaideI1, "%12s","k (MeV)");
    fprintf (pFileHaideI1, "%12s","r (fm)");
    fprintf (pFileHaideI1, "%12s","u (1/fm^2)");
    fprintf(pFileHaideI1,"\n");
    for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
        for(unsigned uRad=0; uRad<NumRadiusBins; uRad++){
            fprintf (pFileHaideI1, "%12.3e",MyCat.GetMomentum(uMomBin));
            fprintf (pFileHaideI1, "%12.3e",hRadiusI1[0][uMomBin]->GetBinCenter(uRad+1));
            fprintf (pFileHaideI1, "%12.3e",hRadiusI1[0][uMomBin]->GetBinContent(uRad+1));
            fprintf(pFileHaideI1,"\n");
        }
    }

    CleanHaidenbauer(MyCat, &H_WaveFunctionU,&H_PhaseShifts,&H_RadBins);
    CleanHaidenbauer(MyCat, &T_WaveFunctionU,&T_PhaseShifts,&T_RadBins);

    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            delete hRadiusI0[uSce][uMomBin];
            delete hRadiusI1[uSce][uMomBin];
        }
        delete [] hRadiusI0[uSce];
        delete [] hRadiusI1[uSce];
        delete hPsI0[uSce];
        delete hPsI1[uSce];
    }
    delete [] hRadiusI0;
    delete [] hRadiusI1;
    delete [] hPsI0;
    delete [] hPsI1;

    delete OutputFile;
}

void KminProton_ImagWF(){
//! CHECK THE INPUT FILES FOR HAIDENBAUER VER 2
//! we have the Kp wave function
//! ws1 is s-wave in the 2J=1 case (J=1/2)
//! wp1 is the p-wave in the 2J=1 case (J=1/2)
//! wp3 is the p-wave in the 2J=3 case (J=3/2)
    const TString OutputDir = "./OutputKaonProton/KminProton_ImagWF/060718/";

    //0 = Haidenbauer (Juelich)
    //1 = Tetsuo (Kyoto) I=0 + I=1
    //2 = Tetsuo Avg in Isospin
    //3 = Tetsuo (Kyoto) I=0 + I=1, exact wave function (taken with a square)
    //4 = Tetsuo (Kyoto) I=0 + I=1, exact wave function (taken with a square + the dirty trick for sigh flipping)
    //5 = Juelich + mass difference
    const int WHICHPOT = 5;

    const double Plot_kMin = 5;
    const double Plot_kMax = 305;
    const unsigned NumPlotPts = 30;
    //CATShisto<double> CatHist(NumPlotPts,Plot_kMin,Plot_kMax);

    unsigned NumRad = 2;
    double Radius[NumRad] = {0.85,1.2};
//const double Radius[NumRad] = {3.0,2.0,1.0};

    //const unsigned NumRad = 41;
    //double Radius[NumRad];
    //for(unsigned uRad=0; uRad<NumRad; uRad++){
    //    Radius[uRad] = 0.5 + double(uRad)*0.01;
    //}

//for ramona with error bands
    //const unsigned NumRad = 6;
    //const double Radius[NumRad] = {1.113,1.144,1.162,1.188,1.213,1.232};

    //double CkVal;
    //double Momentum;
    DLM_DtColor DlmCol;

    const unsigned NumScenarios = 2;
    TString ScenarioName[NumScenarios];
    ScenarioName[0] = "True PWF";
    ScenarioName[1] = "Re PWF";
    //ScenarioName[2] = "Reid SC as in CRAB";

    TString ScenarioShortName[NumScenarios];
    ScenarioShortName[0] = "NoSI";
    ScenarioShortName[1] = "WithSI";

    //0 = |PW|
    //2 = Imaginary WF
    const int TYPE[NumScenarios] = {2,0};

    double SourcePars_Kp[4] = {0,0,0,Radius[0]};

    CATS MyCat;
    MyCat.SetUseAnalyticSource(true);
    MyCat.SetAnaSource(GaussSource, SourcePars_Kp);

    MyCat.SetExcludeFailedBins(false);
    MyCat.SetMomBins(NumPlotPts,Plot_kMin,Plot_kMax);

    if(WHICHPOT==0){
        MyCat.SetNumChannels(5);
        MyCat.SetNumPW(0,2);
        MyCat.SetNumPW(1,2);
        MyCat.SetNumPW(2,2);
        MyCat.SetNumPW(3,2);
        MyCat.SetNumPW(4,1);
        //change to double?
        MyCat.SetSpin(0,1);//s01+p01; I=0; J(p-wave)=1/2
        MyCat.SetSpin(1,1);//s01+p03; I=0; J(p-wave)=3/2
        MyCat.SetSpin(2,1);//s11+p01; I=1; J(p-wave)=1/2
        MyCat.SetSpin(3,1);//s11+p03; I=1; J(p-wave)=3/2
        MyCat.SetSpin(4,1);//the no SI channel
    }
    else if(WHICHPOT==1){
        MyCat.SetNumChannels(3);
        MyCat.SetNumPW(0,1);
        MyCat.SetNumPW(1,1);
        MyCat.SetNumPW(2,1);
        //change to double?
        MyCat.SetSpin(0,1);//I=0
        MyCat.SetSpin(1,1);//I=1
        MyCat.SetSpin(2,1);//the no SI channel
    }
    else if(WHICHPOT==2){
        MyCat.SetNumChannels(2);
        MyCat.SetNumPW(0,1);
        MyCat.SetNumPW(1,1);
        //change to double?
        MyCat.SetSpin(0,1);
        MyCat.SetSpin(1,1);//the no SI channel
    }
    else if(WHICHPOT==3){
        MyCat.SetNumChannels(2);
        MyCat.SetNumPW(0,1);
        MyCat.SetNumPW(1,1);
        //change to double?
        MyCat.SetSpin(0,1);
        MyCat.SetSpin(1,1);//the no SI channel
    }
    else if(WHICHPOT==4){
        MyCat.SetNumChannels(2);
        MyCat.SetNumPW(0,1);
        MyCat.SetNumPW(1,1);
        //change to double?
        MyCat.SetSpin(0,1);
        MyCat.SetSpin(1,1);//the no SI channel
    }
    else if(WHICHPOT==5){
        MyCat.SetNumChannels(3);
        MyCat.SetNumPW(0,1);
        MyCat.SetNumPW(1,1);
        MyCat.SetNumPW(2,1);
        //change to double?
        MyCat.SetSpin(0,1);//I=0
        MyCat.SetSpin(1,1);//I=1
        MyCat.SetSpin(2,1);//the no SI channel
    }
    else{
        printf("BAD!\n");
    }

    MyCat.SetQ1Q2(-1);
    MyCat.SetPdgId(2212, 222);//check the pdg id of the kaon
    //const double Mass_p = 938.272;
    const double Mass_Km = 493.677;
    const double Mass_p = 938.272;
    MyCat.SetRedMass( (Mass_p*Mass_Km)/(Mass_p+Mass_Km) );

    //in bins of Scenario/Radius
    TGraph** grCk;
    grCk = new TGraph* [NumScenarios];

    TGraph** grRatio;
    grRatio = new TGraph* [NumScenarios];

    TGraph* grWf2 = new TGraph [NumScenarios];

    TGraph* grPwRe = new TGraph [NumScenarios];
    TGraph* grPwIm = new TGraph [NumScenarios];
    TGraph* grPw2 = new TGraph [NumScenarios];

    const double AtMom = 100;
    const unsigned NumRadPlotPts = 128;
    const double MaxPadPlot = 4;

    complex<double>***** WaveFunctionU = new complex<double>**** [NumScenarios];
    double**** PhaseShifts = new double*** [2];
    double** RadBins = new double* [2];
    unsigned* NumRadBins = new unsigned [2];
    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        NumRadBins[uSce]=0;
        if(WHICHPOT==0){
            InitHaidenbauerKaonMinus("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/Haidenbauer/KaonMinus10MeV/",
                                     MyCat,&WaveFunctionU[uSce],&PhaseShifts[uSce],&RadBins[uSce],NumRadBins[uSce],TYPE[uSce]);
        }
        else if(WHICHPOT==3){
            InitTetsuoKaonMinus("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/Tetsuo/KminusProton/"
                                ,MyCat,&WaveFunctionU[uSce],&PhaseShifts[uSce],&RadBins[uSce],NumRadBins[uSce], TYPE[uSce]);
        }
        else{
            InitHaidenbauerKaonMinus_ver2("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/Haidenbauer/KaonMinus10MeV/ver2/",
                                          MyCat,&WaveFunctionU[uSce],&PhaseShifts[uSce],&RadBins[uSce],NumRadBins[uSce], TYPE[uSce]);
        }
    }

//    double Kp_potPars[2];

    //double pKmPotParsI0[10]={0,0,pKm_Tetsuo,0,0,1,1,1,0,1};
    //double pKmPotParsI1[10]={0,0,pKm_Tetsuo,0,1,1,1,1,0,1};
    double pKmPotParsI0[10]={0,0,pKm_Tetsuo,0,0,1,1,1,0,1};
    double pKmPotParsI1[10]={0,0,pKm_Tetsuo,0,1,1,1,1,0,1};
    double pKmPotParsIavg[10]={0,0,pKm_Tetsuo,2,0,1,1,1,0,1};

    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        grCk[uSce] = new TGraph[NumRad];
        grRatio[uSce] = new TGraph[NumRad];

        grWf2[uSce].Set(NumRadPlotPts);
        grPwRe[uSce].Set(NumRadPlotPts);
        grPwIm[uSce].Set(NumRadPlotPts);
        grPw2[uSce].Set(NumRadPlotPts);

        grWf2[uSce].SetName(TString::Format("grWf2_%u",uSce));
        grPwRe[uSce].SetName(TString::Format("grPwRe_%u",uSce));
        grPwIm[uSce].SetName(TString::Format("grPwIm_%u",uSce));
        grPw2[uSce].SetName(TString::Format("grPw2_%u",uSce));

        switch(uSce){
        case 0 :
            if(WHICHPOT==0){
                MyCat.SetChannelWeight(0, 1./6.);
                MyCat.SetChannelWeight(1, 1./3.);
                MyCat.SetChannelWeight(2, 1./6.);
                MyCat.SetChannelWeight(3, 1./3.);
                MyCat.SetChannelWeight(4, 0.);
                for(unsigned uMomBin=0; uMomBin<NumPlotPts; uMomBin++){
                    MyCat.UseExternalWaveFunction(uMomBin,0,0,WaveFunctionU[uSce][uMomBin][0][0],NumRadBins[uSce],RadBins[uSce],PhaseShifts[uSce][uMomBin][0][0]);
                    MyCat.UseExternalWaveFunction(uMomBin,0,1,WaveFunctionU[uSce][uMomBin][0][1],NumRadBins[uSce],RadBins[uSce],PhaseShifts[uSce][uMomBin][0][1]);
                    MyCat.UseExternalWaveFunction(uMomBin,1,0,WaveFunctionU[uSce][uMomBin][1][0],NumRadBins[uSce],RadBins[uSce],PhaseShifts[uSce][uMomBin][1][0]);
                    MyCat.UseExternalWaveFunction(uMomBin,1,1,WaveFunctionU[uSce][uMomBin][1][1],NumRadBins[uSce],RadBins[uSce],PhaseShifts[uSce][uMomBin][1][1]);
                    MyCat.UseExternalWaveFunction(uMomBin,2,0,WaveFunctionU[uSce][uMomBin][2][0],NumRadBins[uSce],RadBins[uSce],PhaseShifts[uSce][uMomBin][2][0]);
                    MyCat.UseExternalWaveFunction(uMomBin,2,1,WaveFunctionU[uSce][uMomBin][2][1],NumRadBins[uSce],RadBins[uSce],PhaseShifts[uSce][uMomBin][2][1]);
                    MyCat.UseExternalWaveFunction(uMomBin,3,0,WaveFunctionU[uSce][uMomBin][3][0],NumRadBins[uSce],RadBins[uSce],PhaseShifts[uSce][uMomBin][3][0]);
                    MyCat.UseExternalWaveFunction(uMomBin,3,1,WaveFunctionU[uSce][uMomBin][3][1],NumRadBins[uSce],RadBins[uSce],PhaseShifts[uSce][uMomBin][3][1]);
                }
            }
            else if(WHICHPOT==1){
                //MyCat.SetChannelWeight(0, 1./2.*0.);
                //MyCat.SetChannelWeight(1, 1./2.*2.);
                MyCat.SetChannelWeight(0, 1./2.);
                MyCat.SetChannelWeight(1, 1./2.);
                MyCat.SetChannelWeight(2, 0.);
                MyCat.SetShortRangePotential(0,0,fDlmPot,pKmPotParsI0);
                MyCat.SetShortRangePotential(1,0,fDlmPot,pKmPotParsI1);
            }
            else if(WHICHPOT==2){
                MyCat.SetChannelWeight(0, 1.);
                MyCat.SetChannelWeight(1, 0.);
                MyCat.SetShortRangePotential(0,0,fDlmPot,pKmPotParsIavg);
            }
            else if(WHICHPOT==3){
                MyCat.SetChannelWeight(0, 1.);
                MyCat.SetChannelWeight(1, 0.);
                for(unsigned uMomBin=0; uMomBin<NumPlotPts; uMomBin++){
                    MyCat.UseExternalWaveFunction(uMomBin,0,0,WaveFunctionU[uSce][uMomBin][0][0],NumRadBins[uSce],RadBins[uSce],PhaseShifts[uSce][uMomBin][0][0]);
                }
            }
            else if(WHICHPOT==4){
                MyCat.SetChannelWeight(0, 1.);
                MyCat.SetChannelWeight(1, 0.);
                for(unsigned uMomBin=0; uMomBin<NumPlotPts; uMomBin++){
                    MyCat.UseExternalWaveFunction(uMomBin,0,0,WaveFunctionU[uSce][uMomBin][0][0],NumRadBins[uSce],RadBins[uSce],PhaseShifts[uSce][uMomBin][0][0]);
                }
            }
            else if(WHICHPOT==5){
                MyCat.SetChannelWeight(0, 1./2.);
                MyCat.SetChannelWeight(1, 1./2.);
                MyCat.SetChannelWeight(2, 0.);
                for(unsigned uMomBin=0; uMomBin<NumPlotPts; uMomBin++){
                    MyCat.UseExternalWaveFunction(uMomBin,0,0,WaveFunctionU[uSce][uMomBin][0][0],NumRadBins[uSce],RadBins[uSce],PhaseShifts[uSce][uMomBin][0][0]);
                    MyCat.UseExternalWaveFunction(uMomBin,1,0,WaveFunctionU[uSce][uMomBin][1][0],NumRadBins[uSce],RadBins[uSce],PhaseShifts[uSce][uMomBin][1][0]);
                }
            }
            else{
                printf("Default K- potential fuck!\n");
            }
            break;
        case 1 :
            if(WHICHPOT==0){
                MyCat.SetChannelWeight(0, 1./6.);
                MyCat.SetChannelWeight(1, 1./3.);
                MyCat.SetChannelWeight(2, 1./6.);
                MyCat.SetChannelWeight(3, 1./3.);
                MyCat.SetChannelWeight(4, 0.);
                for(unsigned uMomBin=0; uMomBin<NumPlotPts; uMomBin++){
                    MyCat.UseExternalWaveFunction(uMomBin,0,0,WaveFunctionU[uSce][uMomBin][0][0],NumRadBins[uSce],RadBins[uSce],PhaseShifts[uSce][uMomBin][0][0]);
                    MyCat.UseExternalWaveFunction(uMomBin,0,1,WaveFunctionU[uSce][uMomBin][0][1],NumRadBins[uSce],RadBins[uSce],PhaseShifts[uSce][uMomBin][0][1]);
                    MyCat.UseExternalWaveFunction(uMomBin,1,0,WaveFunctionU[uSce][uMomBin][1][0],NumRadBins[uSce],RadBins[uSce],PhaseShifts[uSce][uMomBin][1][0]);
                    MyCat.UseExternalWaveFunction(uMomBin,1,1,WaveFunctionU[uSce][uMomBin][1][1],NumRadBins[uSce],RadBins[uSce],PhaseShifts[uSce][uMomBin][1][1]);
                    MyCat.UseExternalWaveFunction(uMomBin,2,0,WaveFunctionU[uSce][uMomBin][2][0],NumRadBins[uSce],RadBins[uSce],PhaseShifts[uSce][uMomBin][2][0]);
                    MyCat.UseExternalWaveFunction(uMomBin,2,1,WaveFunctionU[uSce][uMomBin][2][1],NumRadBins[uSce],RadBins[uSce],PhaseShifts[uSce][uMomBin][2][1]);
                    MyCat.UseExternalWaveFunction(uMomBin,3,0,WaveFunctionU[uSce][uMomBin][3][0],NumRadBins[uSce],RadBins[uSce],PhaseShifts[uSce][uMomBin][3][0]);
                    MyCat.UseExternalWaveFunction(uMomBin,3,1,WaveFunctionU[uSce][uMomBin][3][1],NumRadBins[uSce],RadBins[uSce],PhaseShifts[uSce][uMomBin][3][1]);
                }
            }
            else if(WHICHPOT==1){
                //MyCat.SetChannelWeight(0, 1./2.*0.);
                //MyCat.SetChannelWeight(1, 1./2.*2.);
                MyCat.SetChannelWeight(0, 1./2.);
                MyCat.SetChannelWeight(1, 1./2.);
                MyCat.SetChannelWeight(2, 0.);
                MyCat.SetShortRangePotential(0,0,fDlmPot,pKmPotParsI0);
                MyCat.SetShortRangePotential(1,0,fDlmPot,pKmPotParsI1);
            }
            else if(WHICHPOT==2){
                MyCat.SetChannelWeight(0, 1.);
                MyCat.SetChannelWeight(1, 0.);
                MyCat.SetShortRangePotential(0,0,fDlmPot,pKmPotParsIavg);
            }
            else if(WHICHPOT==3){
                MyCat.SetChannelWeight(0, 1.);
                MyCat.SetChannelWeight(1, 0.);
                for(unsigned uMomBin=0; uMomBin<NumPlotPts; uMomBin++){
                    MyCat.UseExternalWaveFunction(uMomBin,0,0,WaveFunctionU[uSce][uMomBin][0][0],NumRadBins[uSce],RadBins[uSce],PhaseShifts[uSce][uMomBin][0][0]);
                }
            }
            else if(WHICHPOT==4){
                MyCat.SetChannelWeight(0, 1.);
                MyCat.SetChannelWeight(1, 0.);
                for(unsigned uMomBin=0; uMomBin<NumPlotPts; uMomBin++){
                    MyCat.UseExternalWaveFunction(uMomBin,0,0,WaveFunctionU[uSce][uMomBin][0][0],NumRadBins[uSce],RadBins[uSce],PhaseShifts[uSce][uMomBin][0][0]);
                }
            }
            else if(WHICHPOT==5){
                MyCat.SetChannelWeight(0, 1./2.);
                MyCat.SetChannelWeight(1, 1./2.);
                MyCat.SetChannelWeight(2, 0.);
                for(unsigned uMomBin=0; uMomBin<NumPlotPts; uMomBin++){
                    MyCat.UseExternalWaveFunction(uMomBin,0,0,WaveFunctionU[uSce][uMomBin][0][0],NumRadBins[uSce],RadBins[uSce],PhaseShifts[uSce][uMomBin][0][0]);
                    MyCat.UseExternalWaveFunction(uMomBin,1,0,WaveFunctionU[uSce][uMomBin][1][0],NumRadBins[uSce],RadBins[uSce],PhaseShifts[uSce][uMomBin][1][0]);
                }
            }
            else{
                printf("Default K- potential fuck!\n");
            }

            break;
        default :
            printf("Default K- fuck!\n");
            break;
        }

        MyCat.KillTheCat(CATS::kAllChanged);

        //grWf[uSce].Set(NumRadPlotPts);
        //grWf2[uSce].Set(NumRadPlotPts);
        //grPwRe[uSce].Set(NumRadPlotPts);
        //grPw2[uSce].Set(NumRadPlotPts);
        double RAD;
        unsigned WhichMomBin = MyCat.GetMomBin(AtMom);
        complex <double> ComplexNum;
        for(unsigned uRadPlot=0; uRadPlot<NumRadPlotPts; uRadPlot++){
            RAD = (double(uRadPlot)+0.5)*MaxPadPlot/double(NumRadPlotPts);
            grWf2[uSce].SetPoint(uRadPlot,RAD,MyCat.EvalWaveFun2(WhichMomBin,RAD,0));
            ComplexNum = MyCat.EvalRadialWaveFunction(WhichMomBin,0,0,RAD,false);
            grPwRe[uSce].SetPoint(uRadPlot,RAD,real(ComplexNum));
            grPwIm[uSce].SetPoint(uRadPlot,RAD,imag(ComplexNum));
            printf("z = %f +i*%f\n",real(ComplexNum),imag(ComplexNum));
            grPw2[uSce].SetPoint(uRadPlot,RAD,pow(abs(ComplexNum),2));
        }


        for(unsigned uRad=0; uRad<NumRad; uRad++){
            MyCat.SetAnaSource(0, Radius[uRad]);
            MyCat.KillTheCat();

            grCk[uSce][uRad].Set(NumPlotPts);
            grCk[uSce][uRad].SetName(TString::Format("grCk_%u_%.2f",uSce,Radius[uRad]));
            grCk[uSce][uRad].SetMarkerColor(DlmCol.GetColor(1));
            grCk[uSce][uRad].SetMarkerStyle(20);
            grCk[uSce][uRad].SetMarkerSize(0);
            grCk[uSce][uRad].SetLineColor(DlmCol.GetColor(1+uSce));
            //the Coulomb only case
            //grCk[uSce][uRad].SetLineWidth(8.-4.*double(uSce)/double(NumScenarios-1));
            grCk[uSce][uRad].SetLineWidth(6);
            grCk[uSce][uRad].SetFillColor(kWhite);

            for(unsigned uMomBin=0; uMomBin<NumPlotPts; uMomBin++){
                grCk[uSce][uRad].SetPoint(uMomBin, MyCat.GetMomentum(uMomBin), MyCat.GetCorrFun(uMomBin));
            }
            if(uSce==0){
                //grCk[uSce][uRad].SetLineColor(kGray+1);
                grCk[uSce][uRad].SetLineStyle(2);
                grCk[uSce][uRad].SetLineWidth(9);
            }

            grRatio[uSce][uRad].Set(NumPlotPts);
            grRatio[uSce][uRad].SetName(TString::Format("grRatio_%u_%.2f",uSce,Radius[uRad]));
            grRatio[uSce][uRad].SetMarkerColor(DlmCol.GetColor(1));
            grRatio[uSce][uRad].SetMarkerStyle(20);
            grRatio[uSce][uRad].SetMarkerSize(0);
            grRatio[uSce][uRad].SetLineColor(DlmCol.GetColor(1+uSce));
            //the Coulomb only case
            //grCk[uSce][uRad].SetLineWidth(8.-4.*double(uSce)/double(NumScenarios-1));
            grRatio[uSce][uRad].SetLineWidth(9);
            grRatio[uSce][uRad].SetFillColor(kWhite);

            double xVal;
            double yVal;

            for(unsigned uMomBin=0; uMomBin<NumPlotPts; uMomBin++){
                grCk[0][uRad].GetPoint(uMomBin,xVal,yVal);
                grRatio[uSce][uRad].SetPoint(uMomBin, MyCat.GetMomentum(uMomBin), MyCat.GetCorrFun(uMomBin)/yVal);
            }
        }
    }

    char* fBuffer = new char [512];
    FILE * pFile;
    pFile = fopen (TString::Format("%sPlotCk_KmProton.txt",OutputDir.Data()),"w");
    fprintf (pFile, "%18s","k (MeV)");
    for(unsigned uRad=0; uRad<NumRad; uRad++){
        for(unsigned uSce=0; uSce<NumScenarios; uSce++){
            sprintf(fBuffer,"R%.2f_%s",Radius[uRad],ScenarioShortName[uSce].Data());
            fprintf (pFile,"%18s",fBuffer);
        }
    }

    double xVal,yVal;
    for(unsigned uMomBin=0; uMomBin<NumPlotPts; uMomBin++){
        fprintf (pFile, "\n%18.3f",MyCat.GetMomentum(uMomBin));
        for(unsigned uRad=0; uRad<NumRad; uRad++){
            for(unsigned uSce=0; uSce<NumScenarios; uSce++){
                grCk[uSce][uRad].GetPoint(uMomBin,xVal,yVal);
                if(xVal!=MyCat.GetMomentum(uMomBin)) {printf("WTF MAN!\n");}
                fprintf (pFile, "%18.3f",yVal);
            }
        }
    }

    const double Yoffset = 0.65;
    TH1F* hAxis = new TH1F("hAxis", "hAxis", NumPlotPts,Plot_kMin,Plot_kMax);
    hAxis->SetStats(false);
    hAxis->SetTitle("");
    hAxis->GetXaxis()->SetLabelSize(0.065);
    hAxis->GetXaxis()->SetTitle("k (MeV)");
    hAxis->GetXaxis()->CenterTitle();
    hAxis->GetXaxis()->SetTitleOffset(1.35);
    hAxis->GetXaxis()->SetLabelOffset(0.02);
    hAxis->GetXaxis()->SetTitleSize(0.075);
    hAxis->GetYaxis()->SetLabelSize(0.065);
    hAxis->GetYaxis()->SetTitle("C(k)");
    hAxis->GetYaxis()->CenterTitle();
    hAxis->GetYaxis()->SetTitleOffset(Yoffset);
    hAxis->GetYaxis()->SetTitleSize(0.075);

    hAxis->GetYaxis()->SetRangeUser(0.5, 2.5);
    hAxis->GetYaxis()->SetLimits(0.5, 2.5);
//hAxis->GetYaxis()->SetRangeUser(0.85, 1.3);
//hAxis->GetYaxis()->SetLimits(0.85, 1.3);

    hAxis->GetXaxis()->SetNdivisions(504);

    TH1F* hAxisZoom = new TH1F("hAxisZoom", "hAxisZoom", NumPlotPts,Plot_kMin,Plot_kMax);
    hAxisZoom->SetStats(false);
    hAxisZoom->SetTitle("");
    hAxisZoom->GetXaxis()->SetLabelSize(0.065);
    hAxisZoom->GetXaxis()->SetTitle("k (MeV)");
    hAxisZoom->GetXaxis()->CenterTitle();
    hAxisZoom->GetXaxis()->SetTitleOffset(1.05);
    hAxisZoom->GetXaxis()->SetLabelOffset(0.02);
    hAxisZoom->GetXaxis()->SetTitleSize(0.075);
    hAxisZoom->GetYaxis()->SetLabelSize(0.065);
    hAxisZoom->GetYaxis()->SetTitle("Ratio");
    hAxisZoom->GetYaxis()->CenterTitle();
    hAxisZoom->GetYaxis()->SetTitleOffset(Yoffset/1.5);
    hAxisZoom->GetYaxis()->SetTitleSize(0.075);

    hAxisZoom->GetXaxis()->SetNdivisions(504);
    hAxisZoom->GetYaxis()->SetRangeUser(0.999, 1.001);
    hAxisZoom->GetYaxis()->SetLimits(0.999, 1.001);
    hAxisZoom->GetYaxis()->SetNdivisions(603);

    TPaveText** myPT = new TPaveText* [NumRad];
    for(unsigned uRad=0; uRad<NumRad; uRad++){
        myPT[uRad] = new TPaveText(0.5,0.80,0.95,0.95, "blNDC");//lbrt
        myPT[uRad]->SetName("PT12");
        myPT[uRad]->SetBorderSize(1);
        myPT[uRad]->SetTextSize(0.08);
        myPT[uRad]->SetFillColor(kWhite);
        myPT[uRad]->SetTextFont(22);
        //myPT[uRad]->AddText("Gauss source");
        myPT[uRad]->AddText(TString::Format("R = %.2f fm", Radius[uRad]));
    }

    TF1* UnitLine = new TF1("UnitLine", "1",Plot_kMin,Plot_kMax);
    UnitLine->SetLineColor(kGray+1);
    UnitLine->SetLineWidth(6);

    DLM_SubPads DrawBoard_Ck_pp(2560,1080);
    double PadLength = 1./double(NumRad);
    double CurrentPos = 0;
    for(unsigned uRad=0; uRad<NumRad; uRad++){
        DrawBoard_Ck_pp.AddSubPadTL(CurrentPos,CurrentPos+PadLength,0.08,0.63);//lrtb
        DrawBoard_Ck_pp.AddSubPadTL(CurrentPos,CurrentPos+PadLength,0.63,1);//lrtb
        CurrentPos+=PadLength;
//printf("uRad=%u; %.2f <--> %.2f\n", uRad, );
    }
    DrawBoard_Ck_pp.AddSubPadTL(0,1,0,0.075);//lrtb

    for(unsigned uRad=0; uRad<NumRad; uRad++){
        DrawBoard_Ck_pp.SetMargin(uRad*2, 0.07, 0.01, 0, 0.04);//lrbt
        DrawBoard_Ck_pp.SetMargin(uRad*2+1, 0.07, 0.01, 0.14, 0);//lrbt
    }
    DrawBoard_Ck_pp.SetMargin(2*NumRad, 0, 0, 0, 0);//lrbt

    for(unsigned uRad=0; uRad<NumRad; uRad++){
        DrawBoard_Ck_pp.GetCanvas()->SetName(TString::Format("DrawBoard_Ck_KmProton_%u",uRad));
        DrawBoard_Ck_pp.cd(uRad*2);
        hAxis->Draw("AXIS");
        UnitLine->Draw("same");
        for(unsigned uSce=0; uSce<NumScenarios; uSce++) grCk[uSce][uRad].Draw("same,C");
        for(unsigned uSce=0; uSce<NumScenarios; uSce++) myPT[uRad]->Draw("same");

        DrawBoard_Ck_pp.SetLabelSize(uRad*2, hAxis->GetXaxis(), 17);
        DrawBoard_Ck_pp.SetLabelSize(uRad*2, hAxis->GetYaxis(), 17);
        DrawBoard_Ck_pp.SetTitleSize(uRad*2, hAxis->GetXaxis(), 17);
        DrawBoard_Ck_pp.SetTitleSize(uRad*2, hAxis->GetYaxis(), 17);
    }
    for(unsigned uRad=0; uRad<NumRad; uRad++){
        DrawBoard_Ck_pp.GetCanvas()->SetName(TString::Format("DrawBoard_Ck_KmProton_%u",uRad));
        DrawBoard_Ck_pp.cd(uRad*2+1);
        hAxisZoom->Draw("AXIS");
        UnitLine->Draw("same");
        //for(unsigned uSce=0; uSce<NumScenarios; uSce++) grCk[uSce][uRad].Draw("same,C");
        for(unsigned uSce=1; uSce<NumScenarios; uSce++) grRatio[uSce][uRad].Draw("same,C");

        DrawBoard_Ck_pp.SetLabelSize(uRad*2+1, hAxisZoom->GetXaxis(), 17);
        DrawBoard_Ck_pp.SetLabelSize(uRad*2+1, hAxisZoom->GetYaxis(), 17);
        DrawBoard_Ck_pp.SetTitleSize(uRad*2+1, hAxisZoom->GetXaxis(), 17);
        DrawBoard_Ck_pp.SetTitleSize(uRad*2+1, hAxisZoom->GetYaxis(), 17);
    }

    TLegend* myLegend = new TLegend(0.12,0,0.95,0.96);//lbrt
    myLegend->SetName("myLegend");
    myLegend->SetTextSize(0.65);
    myLegend->SetNColumns(NumScenarios);
    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        myLegend->AddEntry(&grCk[uSce][0],  ScenarioName[uSce]);
    }

    DrawBoard_Ck_pp.cd(NumRad*2);
    myLegend->Draw();

    gStyle->SetLineWidth(2.5);
    DrawBoard_Ck_pp.GetCanvas()->SaveAs(TString::Format("%sPlotCk_KmProton.png",OutputDir.Data()));
    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        for(unsigned uRad=0; uRad<NumRad; uRad++){
            grCk[uSce][uRad].SetLineWidth(grCk[uSce][uRad].GetLineWidth()/2.5);
            grRatio[uSce][uRad].SetLineWidth(grRatio[uSce][uRad].GetLineWidth()/2.5);
        }
    }
    hAxis->SetLineWidth(hAxis->GetLineWidth()/2.5);
    UnitLine->SetLineWidth(UnitLine->GetLineWidth()/2.5);
    gStyle->SetLineWidth(1);
    DrawBoard_Ck_pp.GetCanvas()->SaveAs(TString::Format("%sPlotCk_KmProton.pdf",OutputDir.Data()));

    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        for(unsigned uRad=0; uRad<NumRad; uRad++){
            grCk[uSce][uRad].SetLineWidth(grCk[uSce][uRad].GetLineWidth()*2.5);
            grRatio[uSce][uRad].SetLineWidth(grRatio[uSce][uRad].GetLineWidth()*2.5);
        }
    }
    hAxis->SetLineWidth(hAxis->GetLineWidth()*2.5);
    UnitLine->SetLineWidth(UnitLine->GetLineWidth()*2.5);

    TFile* RootFile = new TFile(TString::Format("%sKmProton.root",OutputDir.Data()),"recreate");
    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        for(unsigned uRad=0; uRad<NumRad; uRad++){
            grCk[uSce][uRad].Write();
            grRatio[uSce][uRad].Write();
        }
        grWf2[uSce].Write();
        grPwRe[uSce].Write();
        grPwIm[uSce].Write();
        grPw2[uSce].Write();
    }


    DrawBoard_Ck_pp.GetCanvas()->Write();






    delete hAxis;
    delete hAxisZoom;
    delete UnitLine;
    delete myLegend;

    for(unsigned uRad=0; uRad<NumRad; uRad++) delete myPT[uRad];
    delete [] myPT;

    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        delete [] grCk[uSce];
        delete [] grRatio[uSce];
    }
    delete [] grCk;
    delete [] grRatio;
    delete RootFile;

    for(unsigned uSce=0; uSce<NumScenarios; uSce++) CleanHaidenbauer(MyCat, &WaveFunctionU[uSce],&PhaseShifts[uSce],&RadBins[uSce]);
    delete [] WaveFunctionU;
    delete [] PhaseShifts;
    delete [] RadBins;
    delete [] NumRadBins;

    delete [] grWf2;
    delete [] grPwRe;
    delete [] grPwIm;
    delete [] grPw2;
*/
}

void TestHaideJuly2018(){
/*
    const unsigned NumModels = 7;
    const bool Gamow = true;
    TGraph* grCk = new TGraph [NumModels];

    DLM_DtColor DlmCol;

    for(unsigned uType=0; uType<NumModels; uType++){

        printf("hi\n");

        complex<double>**** WaveFunctionU; double*** PhaseShifts; double* RadBins; unsigned NumRadBins=0;
        CATS cat1;
        if(uType<6) InitCatForHaidenbauerKaonProton1("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/Haidenbauer/OriginalJuelich/",
                cat1,&WaveFunctionU,&PhaseShifts,&RadBins,NumRadBins,uType);
        else{

            cat1.SetMomBins(20, 5, 205);
            cat1.SetNumChannels(1);
            cat1.SetChannelWeight(0,1);
            cat1.SetNumPW(0,1);
            cat1.SetSpin(0,1);

            const double Mass_p = 938.272;
            //const double Mass_n = 939.565;
            const double Mass_Km = 493.677;

            cat1.SetQ1Q2(-1);
            cat1.SetPdgId(2212, -321);
            cat1.SetRedMass( (Mass_p*Mass_Km)/(Mass_p+Mass_Km) );

            InitHaidenbauerKaonMinus_ver2("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/Haidenbauer/KaonMinus10MeV/ver2/",cat1,&WaveFunctionU,&PhaseShifts,&RadBins,NumRadBins,1);
        }

        double SourcePars_Kp[4] = {0,0,0,1.2};
        cat1.SetUseAnalyticSource(true);
        cat1.SetAnaSource(GaussSource, SourcePars_Kp);
        cat1.SetExcludeFailedBins(false);

        unsigned NumMomBins = cat1.GetNumMomBins();
        printf("NumMomBins = %u\n",NumMomBins);
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            cat1.UseExternalWaveFunction(uMomBin,0,0,WaveFunctionU[uMomBin][0][0],NumRadBins,RadBins,PhaseShifts[uMomBin][0][0]);
        }
        cat1.SetGamow(Gamow);
        cat1.KillTheCat(CATS::kAllChanged);

        grCk[uType].Set(NumMomBins);
        grCk[uType].SetName(TString::Format("grCk[%u]",uType));
        grCk[uType].SetMarkerColor(DlmCol.GetColor(uType));
        grCk[uType].SetMarkerStyle(20);
        grCk[uType].SetMarkerSize(0);
        grCk[uType].SetLineColor(DlmCol.GetColor(uType));
        grCk[uType].SetLineWidth(6-0.5*double(uType));
        grCk[uType].SetFillColor(kWhite);
        for(unsigned uBin=0; uBin<NumMomBins; uBin++){
            printf("B%u: C(%.2f)=%.2f\n",uBin,cat1.GetMomentum(uBin),cat1.GetCorrFun(uBin));
            //complex<double> RWF = cat1.EvalRadialWaveFunction(uBin,0,0,2,true);
            //printf("B%u: C(%.2f)=%.2f\n",uBin,cat1.GetMomentum(uBin),real(RWF));
            grCk[uType].SetPoint(uBin,cat1.GetMomentum(uBin),cat1.GetCorrFun(uBin));
        }

        CleanHaidenbauer(cat1,&WaveFunctionU,&PhaseShifts,&RadBins);

    }

    TFile* fout = new TFile("./OutputKaonProton/TestHaideJuly2018/fout.root", "recreate");
    for(unsigned uType=0; uType<NumModels; uType++){
        grCk[uType].Write();
    }

    delete [] grCk;
    delete fout;
*/
}

void TestHaideSiddharta(const double& RADIUS){
    //1.03 - 1.14 - 1.25
    //1.07 - 1.19 - 1.31
/*
    const unsigned NumModels = 7;
    const bool Gamow = true;
    TGraph* grCk = new TGraph [NumModels];
    TGraph* grCkSmeared = new TGraph [NumModels];

    DLM_DtColor DlmCol;

    TString ResolutionFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/ResolutionMatrices/Sample6_MeV_compact.root";
    TFile* ResolutionFile = new TFile(ResolutionFileName, "read");
    ResolutionFile->cd();
    TH2F* hResolution_pKmin = (TH2F*)ResolutionFile->Get("hSigmaMeV_Proton_Proton");

    for(unsigned uType=0; uType<NumModels; uType++){
        //printf("hi\n");
        complex<double>**** WaveFunctionU; double*** PhaseShifts; double* RadBins; unsigned NumRadBins=0;
        CATS cat1;
        if(uType<6) InitCatForHaidenbauerKaonProton2("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/Haidenbauer/Siddharta1/",
                cat1,&WaveFunctionU,&PhaseShifts,&RadBins,NumRadBins,uType);
        else{

            cat1.SetMomBins(20, 5, 205);
            cat1.SetNumChannels(1);
            cat1.SetChannelWeight(0,1);
            cat1.SetNumPW(0,1);
            cat1.SetSpin(0,1);

            const double Mass_p = 938.272;
            //const double Mass_n = 939.565;
            const double Mass_Km = 493.677;

            cat1.SetQ1Q2(-1);
            cat1.SetPdgId(2212, -321);
            cat1.SetRedMass( (Mass_p*Mass_Km)/(Mass_p+Mass_Km) );

            InitHaidenbauerKaonMinus_ver2("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/Haidenbauer/KaonMinus10MeV/ver2/",cat1,&WaveFunctionU,&PhaseShifts,&RadBins,NumRadBins,1);
        }
//cat1.SetQ1Q2(0);
        cat1.SetUseAnalyticSource(true);
        CATSparameters SourcePars_pKmin(CATSparameters::tSource,1,true);
        SourcePars_pKmin.SetParameter(0,RADIUS);
        cat1.SetAnaSource(GaussSource, SourcePars_pKmin);
        cat1.SetExcludeFailedBins(false);

        unsigned NumMomBins = cat1.GetNumMomBins();
        //printf("NumMomBins = %u\n",NumMomBins);
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            cat1.UseExternalWaveFunction(uMomBin,0,0,WaveFunctionU[uMomBin][0][0],NumRadBins,RadBins,PhaseShifts[uMomBin][0][0]);
        }
        cat1.SetGamow(Gamow);
        if(uType==0) cat1.SetMaxPw(256);
        else cat1.SetMaxPw(1);//we add ONLY the s-wave to the total wave function
        cat1.KillTheCat(CATS::kAllChanged);

        DLM_Ck Ck_pKmin(0,0,cat1);
        Ck_pKmin.Update();

        DLM_CkDecomposition CkDec_pKmin("pKmin",0,Ck_pKmin,hResolution_pKmin);
        CkDec_pKmin.Update();

        grCk[uType].Set(NumMomBins);
        grCk[uType].SetName(TString::Format("grCk[%u]",uType));
        grCk[uType].SetMarkerColor(DlmCol.GetColor(uType));
        grCk[uType].SetMarkerStyle(20);
        grCk[uType].SetMarkerSize(0);
        grCk[uType].SetLineColor(DlmCol.GetColor(uType));
        grCk[uType].SetLineWidth(6-0.5*double(uType));
        grCk[uType].SetFillColor(kWhite);

        grCkSmeared[uType].Set(NumMomBins);
        grCkSmeared[uType].SetName(TString::Format("grCkSmeared[%u]",uType));
        grCkSmeared[uType].SetMarkerColor(DlmCol.GetColor(uType));
        grCkSmeared[uType].SetMarkerStyle(20);
        grCkSmeared[uType].SetMarkerSize(0);
        grCkSmeared[uType].SetLineColor(DlmCol.GetColor(uType));
        grCkSmeared[uType].SetLineWidth(6-0.5*double(uType));
        grCkSmeared[uType].SetFillColor(kWhite);

        for(unsigned uBin=0; uBin<NumMomBins; uBin++){
            //printf("B%u: C(%.2f)=%.2f\n",uBin,cat1.GetMomentum(uBin),cat1.GetCorrFun(uBin));
            //complex<double> RWF = cat1.EvalRadialWaveFunction(uBin,0,0,2,true);
            //printf("B%u: C(%.2f)=%.2f\n",uBin,cat1.GetMomentum(uBin),real(RWF));
            grCk[uType].SetPoint(uBin,cat1.GetMomentum(uBin),cat1.GetCorrFun(uBin));
            grCkSmeared[uType].SetPoint(uBin,cat1.GetMomentum(uBin),CkDec_pKmin.EvalCk(cat1.GetMomentum(uBin)));
        }

        CleanHaidenbauer(cat1,&WaveFunctionU,&PhaseShifts,&RadBins);

    }

    unsigned uNumBins=grCk[0].GetN();
    TGraph gTest1;
    gTest1.Set(uNumBins);
    gTest1.SetName(TString::Format("gTest1"));
    gTest1.SetMarkerColor(DlmCol.GetColor(7));
    gTest1.SetMarkerStyle(20);
    gTest1.SetMarkerSize(0);
    gTest1.SetLineColor(DlmCol.GetColor(7));
    gTest1.SetLineWidth(6-0.5*double(7));
    gTest1.SetFillColor(kWhite);
    double* xVal = new double [6];
    double* yVal = new double [6];
    double yValTot;

    //for(unsigned uBin=0; uBin<uNumBins; uBin++){
    //    grCk[0].GetPoint(uBin,xVal1,yVal1);
    //    grCk[1].GetPoint(uBin,xVal2,yVal2);
    //    gTest1.SetPoint(uBin,xVal1,(yVal1+yVal2));
    //}

    for(unsigned uBin=0; uBin<uNumBins; uBin++){
        yValTot = 0;
        for(unsigned uType=0; uType<6; uType++){
            grCk[uType].GetPoint(uBin,xVal[uType],yVal[uType]);
            yValTot += yVal[uType];
            //printf("UT=%u; xVal=%.2f; yVal=%.2f\n",uType,xVal[uType],yVal[uType]);
        }
       // yValTot *= 1.25; //same norm as data
        //yValTot = 1.2*(1.+0.68*(yValTot-1.));//bring it all on the level of the data
        //printf(" yValTot=%.2f\n",yValTot);
        gTest1.SetPoint(uBin,xVal[0],yValTot);
    }


    TGraph gTest1Smeared;
    gTest1Smeared.Set(uNumBins);
    gTest1Smeared.SetName(TString::Format("gTest1Smeared"));
    gTest1Smeared.SetMarkerColor(DlmCol.GetColor(7));
    gTest1Smeared.SetMarkerStyle(20);
    gTest1Smeared.SetMarkerSize(0);
    gTest1Smeared.SetLineColor(DlmCol.GetColor(7));
    gTest1Smeared.SetLineWidth(6-0.5*double(7));
    gTest1Smeared.SetFillColor(kWhite);
    for(unsigned uBin=0; uBin<uNumBins; uBin++){
        yValTot = 0;
        for(unsigned uType=0; uType<6; uType++){
            grCkSmeared[uType].GetPoint(uBin,xVal[uType],yVal[uType]);
            yValTot += yVal[uType];
        }
        gTest1Smeared.SetPoint(uBin,xVal[0],yValTot);
    }

    TFile* fout = new TFile(TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/KaonProton/TestHaideSiddharta/190219/fout_R%.2f.root", RADIUS),
                            "recreate");
    for(unsigned uType=0; uType<NumModels; uType++){
        grCk[uType].Write();
        grCkSmeared[uType].Write();
    }
    gTest1.Write();
    gTest1Smeared.Write();

    delete [] xVal;
    delete [] yVal;
    delete [] grCk;
    delete [] grCkSmeared;
    delete ResolutionFile;
    delete fout;
*/
}

void GamowCheck(){
/*
    const TString OutputFolder = "./OutputKaonProton/GamowCheck/";

    CATS CoulombCat;
    CATS StrongCat;
    CATS** cat = new CATS* [2];
    cat[0] = &CoulombCat;//Coulomb + Strong
    cat[1] = &StrongCat;//Strong

    const double kMin = 5;
    const double kMax = 305;
    const unsigned NumMomBins = 30;

    double SourcePars_Kp[4] = {0,0,0,1.2};

    complex<double>**** cWaveFunctionU; double*** cPhaseShifts; double* cRadBins; unsigned cNumRadBins=0;
    complex<double>**** sWaveFunctionU; double*** sPhaseShifts; double* sRadBins; unsigned sNumRadBins=0;

    const double Mass_p = 938.272;
    const double Mass_Km = 493.677;
    const double RedMass = (Mass_p*Mass_Km)/(Mass_p+Mass_Km);

    for(unsigned uInter=0; uInter<2; uInter++){
        cat[uInter]->SetUseAnalyticSource(true);
        cat[uInter]->SetAnaSource(GaussSource, SourcePars_Kp);
        cat[uInter]->SetExcludeFailedBins(false);
        cat[uInter]->SetMomBins(NumMomBins,kMin,kMax);

        cat[uInter]->SetNumChannels(2);
        cat[uInter]->SetNumPW(0,1);
        cat[uInter]->SetNumPW(1,1);

        cat[uInter]->SetSpin(0,1);
        cat[uInter]->SetSpin(0,1);

        cat[uInter]->SetQ1Q2(-1);
        cat[uInter]->SetPdgId(2212, 222);//check the pdg id of the kaon
        cat[uInter]->SetRedMass( RedMass );

        if(uInter==0){
            InitHaidenbauerKaonMinus("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Haidenbauer/KaonMinus10MeV/",CoulombCat,&cWaveFunctionU,&cPhaseShifts,&cRadBins,cNumRadBins);
            for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
                cat[uInter]->UseExternalWaveFunction(uMomBin,0,0,cWaveFunctionU[uMomBin][0][0],cNumRadBins,cRadBins,cPhaseShifts[uMomBin][0][0]);
                cat[uInter]->UseExternalWaveFunction(uMomBin,1,0,cWaveFunctionU[uMomBin][2][0],cNumRadBins,cRadBins,cPhaseShifts[uMomBin][0][0]);
            }
        }
        else{
            InitHaidenbauerKaonMinus("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Haidenbauer/KaonMinus_SIonly/",StrongCat,&sWaveFunctionU,&sPhaseShifts,&sRadBins,sNumRadBins);
            for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
                cat[uInter]->UseExternalWaveFunction(uMomBin,0,0,sWaveFunctionU[uMomBin][0][0],sNumRadBins,sRadBins,sPhaseShifts[uMomBin][0][0]);
                cat[uInter]->UseExternalWaveFunction(uMomBin,1,0,sWaveFunctionU[uMomBin][2][0],sNumRadBins,sRadBins,sPhaseShifts[uMomBin][0][0]);
            }
        }

        cat[uInter]->SetChannelWeight(0,0.5);
        cat[uInter]->SetChannelWeight(1,0.5);

        cat[uInter]->KillTheCat();
    }

    if(sNumRadBins!=cNumRadBins){
        printf("Something bad has happened\n");
        return;
    }

    TGraph* grGamowTh = new TGraph [NumMomBins];
    TGraph* grGamowHaide = new TGraph [NumMomBins];
    TGraph* grRatio = new TGraph [NumMomBins];

    TGraph grGamowHaideMom;
    grGamowHaideMom.Set(NumMomBins);
    grGamowHaideMom.SetName("grGamowHaideMom");

    TGraph grGamowThMom;
    grGamowThMom.Set(NumMomBins);
    grGamowThMom.SetName("grGamowThMom");

    TGraph grGamowRatioMom;
    grGamowRatioMom.Set(NumMomBins);
    grGamowRatioMom.SetName("grGamowRatioMom");

    TFile* OutFile = new TFile(OutputFolder+"OutFile.root","recreate");

    for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
        grGamowTh[uMomBin].Set(sNumRadBins);
        grGamowHaide[uMomBin].Set(sNumRadBins);
        grRatio[uMomBin].Set(sNumRadBins);

        grGamowTh[uMomBin].SetName(TString::Format("grGamowTh_%u",uMomBin));
        grGamowHaide[uMomBin].SetName(TString::Format("grGamowHaide_%u",uMomBin));
        grRatio[uMomBin].SetName(TString::Format("grRatio_%u",uMomBin));

        grGamowTh[uMomBin].SetLineColor(kRed);
        grGamowHaide[uMomBin].SetLineColor(kBlue);

        double Radius;
        double Momentum;
        complex<double> sWAVEFUN;
        complex<double> cWAVEFUN;
        complex<double> GAMOW;

        Momentum = StrongCat.GetMomentum(uMomBin);

        double CPF = CoulombPenetrationFactor(CoulombEta(Momentum,RedMass,-1));

        //printf("C(%.2f)=%.2f\n",Momentum,CoulombCat.GetCorrFun(uMomBin));

        for(unsigned uRad=0; uRad<sNumRadBins; uRad++){
            Radius = (sRadBins[uRad]+sRadBins[uRad+1])*0.5;
            sWAVEFUN = StrongCat.EvalRadialWaveFunction(uMomBin,0,0,Radius);
            cWAVEFUN = CoulombCat.EvalRadialWaveFunction(uMomBin,0,0,Radius);
            GAMOW = cWAVEFUN/sWAVEFUN;
            grGamowHaide[uMomBin].SetPoint(uRad,Radius,abs(GAMOW)*abs(GAMOW));
            grGamowTh[uMomBin].SetPoint(uRad,Radius,CPF);
            grRatio[uMomBin].SetPoint(uRad,Radius,abs(GAMOW)*abs(GAMOW)/CPF);

            if(uMomBin==6){
                printf("R=%.2f; H %e; T %e; R %e\n",Radius, abs(GAMOW)*abs(GAMOW),CPF,abs(GAMOW)*abs(GAMOW)/CPF);
            }
        }

        grGamowHaideMom.SetPoint(uMomBin,Momentum,abs(GAMOW)*abs(GAMOW));
        grGamowThMom.SetPoint(uMomBin,Momentum,CPF);
        grGamowRatioMom.SetPoint(uMomBin,Momentum,abs(GAMOW)*abs(GAMOW)/CPF);

        grGamowHaide[uMomBin].Write();
        grGamowTh[uMomBin].Write();
        grRatio[uMomBin].Write();
    }

    grGamowHaideMom.Write();
    grGamowThMom.Write();
    grGamowRatioMom.Write();

    CleanHaidenbauer(CoulombCat, &cWaveFunctionU,&cPhaseShifts,&cRadBins);
    CleanHaidenbauer(StrongCat, &sWaveFunctionU,&sPhaseShifts,&sRadBins);

    delete [] grGamowTh;
    delete [] grGamowHaide;
    delete [] grRatio;
    delete [] cat;

    delete OutFile;
*/
}

void GamowCheck2(){
/*
    const TString OutputFolder = "./OutputKaonProton/GamowCheck2/";

    CATS catHaide;
    CATS catGamow;
    CATS catHaide2;

    const unsigned NumSce=3;

    CATS** cat = new CATS* [NumSce];
    cat[0] = &catHaide;
    cat[1] = &catGamow;
    cat[2] = &catHaide2;

    TFile* OutFile = new TFile(OutputFolder+"OutFile.root","recreate");

    TGraph* grCk = new TGraph[NumSce];

    const double kMin = 5;
    const double kMax = 305;
    const unsigned NumMomBins = 30;

    double SourcePars_Kp[4] = {0,0,0,1.2};

    complex<double>**** cWaveFunctionU; double*** cPhaseShifts; double* cRadBins; unsigned cNumRadBins=0;
    complex<double>**** sWaveFunctionU; double*** sPhaseShifts; double* sRadBins; unsigned sNumRadBins=0;
    complex<double>**** h2WaveFunctionU; double*** h2PhaseShifts; double* h2RadBins; unsigned h2NumRadBins=0;

    const double Mass_p = 938.272;
    const double Mass_Km = 493.677;
    const double RedMass = (Mass_p*Mass_Km)/(Mass_p+Mass_Km);

    for(unsigned uInter=0; uInter<NumSce; uInter++){
        cat[uInter]->SetUseAnalyticSource(true);
        cat[uInter]->SetAnaSource(GaussSource, SourcePars_Kp);
        cat[uInter]->SetExcludeFailedBins(false);
        cat[uInter]->SetMomBins(NumMomBins,kMin,kMax);

        cat[uInter]->SetNumChannels(2);
        cat[uInter]->SetNumPW(0,1);
        cat[uInter]->SetNumPW(1,1);

        cat[uInter]->SetSpin(0,1);
        cat[uInter]->SetSpin(0,1);

        cat[uInter]->SetQ1Q2(-1);
        cat[uInter]->SetPdgId(2212, 222);//check the pdg id of the kaon
        cat[uInter]->SetRedMass( RedMass );

        cat[uInter]->SetChannelWeight(0,0.5);
        cat[uInter]->SetChannelWeight(1,0.5);

        if(uInter==0){
            InitHaidenbauerKaonMinus("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Haidenbauer/KaonMinus10MeV/",catHaide,&cWaveFunctionU,&cPhaseShifts,&cRadBins,cNumRadBins);
            for(unsigned uMomBin=0; uMomBin<catHaide.GetNumMomBins(); uMomBin++){
                cat[uInter]->UseExternalWaveFunction(uMomBin,0,0,cWaveFunctionU[uMomBin][0][0],cNumRadBins,cRadBins,cPhaseShifts[uMomBin][0][0]);
                cat[uInter]->UseExternalWaveFunction(uMomBin,1,0,cWaveFunctionU[uMomBin][2][0],cNumRadBins,cRadBins,cPhaseShifts[uMomBin][0][0]);
            }
        }
        else if(uInter==1){
            InitHaidenbauerKaonMinus("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Haidenbauer/KaonMinus_SIonly/",catGamow,&sWaveFunctionU,&sPhaseShifts,&sRadBins,sNumRadBins);
            for(unsigned uMomBin=0; uMomBin<catGamow.GetNumMomBins(); uMomBin++){
                cat[uInter]->UseExternalWaveFunction(uMomBin,0,0,sWaveFunctionU[uMomBin][0][0],sNumRadBins,sRadBins,sPhaseShifts[uMomBin][0][0]);
                cat[uInter]->UseExternalWaveFunction(uMomBin,1,0,sWaveFunctionU[uMomBin][2][0],sNumRadBins,sRadBins,sPhaseShifts[uMomBin][0][0]);
                cat[uInter]->SetGamow(true);
            }
        }
        else{
            InitHaidenbauerKaonMinus_ver2("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/Haidenbauer/KaonMinus10MeV/ver2/",
                                          catHaide2,&h2WaveFunctionU,&h2PhaseShifts,&h2RadBins,h2NumRadBins, 1);
            for(unsigned uMomBin=0; uMomBin<catHaide2.GetNumMomBins(); uMomBin++){
                cat[uInter]->UseExternalWaveFunction(uMomBin,0,0,h2WaveFunctionU[uMomBin][0][0],h2NumRadBins,h2RadBins,h2PhaseShifts[uMomBin][0][0]);
                cat[uInter]->SetGamow(false);
            }
            cat[uInter]->SetChannelWeight(0,1.);
            cat[uInter]->SetChannelWeight(1,0.);
        }


        cat[uInter]->KillTheCat();

        grCk[uInter].SetName(TString::Format("grCk_%u",uInter));
        grCk[uInter].Set(NumMomBins);
        for(unsigned uBin=0; uBin<NumMomBins; uBin++){
            grCk[uInter].SetPoint(uBin,cat[uInter]->GetMomentum(uBin),cat[uInter]->GetCorrFun(uBin));
        }
        grCk[uInter].Write();
    }
    TGraph grRatio;
    grRatio.SetName("grRatio");
    grRatio.Set(NumMomBins);
    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        grRatio.SetPoint(uBin,catHaide.GetMomentum(uBin),catHaide.GetCorrFun(uBin)/catGamow.GetCorrFun(uBin));
    }
    grRatio.Write();

    CleanHaidenbauer(catHaide, &cWaveFunctionU,&cPhaseShifts,&cRadBins);
    CleanHaidenbauer(catGamow, &sWaveFunctionU,&sPhaseShifts,&sRadBins);
    CleanHaidenbauer(catHaide2, &h2WaveFunctionU,&h2PhaseShifts,&h2RadBins);
    delete [] cat;
    delete [] grCk;
    delete OutFile;
*/
}

void Toy_pKplus(){
    const unsigned NumMomBins = 150;
    const double kMin = 0;
    const double kMax = 300;
    const double Mass_p = 938.272;
    const double Mass_KaonCh = 493.677;
    const TString OutputDir = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/KaonProton/Toy_pKplus/";
    CATSparameters cPars(CATSparameters::tSource,1,true);
    cPars.SetParameter(0,1.2);

    CATS KittyCoulomb;
    KittyCoulomb.SetMomBins(NumMomBins,kMin,kMax);
    KittyCoulomb.SetAnaSource(GaussSource, cPars);
    KittyCoulomb.SetUseAnalyticSource(true);
    KittyCoulomb.SetMomentumDependentSource(false);
    KittyCoulomb.SetThetaDependentSource(false);
    //should you include in the result any bins, where the Schroedinger solver failed
    KittyCoulomb.SetExcludeFailedBins(false);
    KittyCoulomb.SetQ1Q2(1);
    KittyCoulomb.SetPdgId(2212, 321);
    KittyCoulomb.SetRedMass( (Mass_p*Mass_KaonCh)/(Mass_p+Mass_KaonCh) );
    KittyCoulomb.SetNumChannels(1);
    KittyCoulomb.SetNumPW(0,1);
    KittyCoulomb.SetSpin(0,0);
    KittyCoulomb.SetChannelWeight(0,1.);
    KittyCoulomb.SetMaxNumThreads(4);
    KittyCoulomb.KillTheCat();


    CATS KittyToyI;
    KittyToyI.SetMomBins(NumMomBins,kMin,kMax);
    KittyToyI.SetAnaSource(GaussSource, cPars);
    KittyToyI.SetUseAnalyticSource(true);
    KittyToyI.SetMomentumDependentSource(false);
    KittyToyI.SetThetaDependentSource(false);
    //should you include in the result any bins, where the Schroedinger solver failed
    KittyToyI.SetExcludeFailedBins(false);
    KittyToyI.SetQ1Q2(1);
    KittyToyI.SetPdgId(2212, 321);
    KittyToyI.SetRedMass( (Mass_p*Mass_KaonCh)/(Mass_p+Mass_KaonCh) );
    KittyToyI.SetNumChannels(1);
    KittyToyI.SetNumPW(0,1);
    //KittyToyI.SetNumPW(1,1);
    KittyToyI.SetSpin(0,0);
    //KittyToyI.SetSpin(1,0);
    KittyToyI.SetChannelWeight(0,1.);
    //KittyToyI.SetChannelWeight(1,1./2.);

    double* SOL_POT_I0 = FindPotentialPars("/home/dmihaylov/Temp/pKplus_1M.root",0.057,0.373,0);
    CATSparameters POT_PARS_I0(CATSparameters::tPotential,4,true);
    POT_PARS_I0.SetParameters(&SOL_POT_I0[2]);

    //double* SOL_POT_I1 = FindPotentialPars("/home/dmihaylov/Temp/pKplus_1M.root",-0.316,0.373,0);
    //double* SOL_POT_I1 = FindPotentialPars("/home/dmihaylov/Temp/pKplus_1M.root",-0.304,0.261,0);
    double* SOL_POT_I1 = FindPotentialPars("/home/dmihaylov/Temp/pKplus_1M.root",-0.310,0.373,0);
    CATSparameters POT_PARS_I1(CATSparameters::tPotential,4,true);
    POT_PARS_I1.SetParameters(&SOL_POT_I1[2]);

    KittyToyI.SetShortRangePotential(0,0,DoubleGaussSum,POT_PARS_I1);
    //KittyToyI.SetShortRangePotential(1,0,DoubleGaussSum,POT_PARS_I0);
    KittyToyI.SetMaxNumThreads(4);
    KittyToyI.SetEpsilonConv(1e-8);
    KittyToyI.SetEpsilonProp(1e-8);
    KittyToyI.KillTheCat();

    CATS KittyToyI_SI;
    KittyToyI_SI.SetMomBins(NumMomBins,kMin,kMax);
    KittyToyI_SI.SetAnaSource(GaussSource, cPars);
    KittyToyI_SI.SetUseAnalyticSource(true);
    KittyToyI_SI.SetMomentumDependentSource(false);
    KittyToyI_SI.SetThetaDependentSource(false);
    //should you include in the result any bins, where the Schroedinger solver failed
    KittyToyI_SI.SetExcludeFailedBins(false);
    KittyToyI_SI.SetQ1Q2(0);
    KittyToyI_SI.SetPdgId(2212, 321);
    KittyToyI_SI.SetRedMass( (Mass_p*Mass_KaonCh)/(Mass_p+Mass_KaonCh) );
    KittyToyI_SI.SetNumChannels(1);
    KittyToyI_SI.SetNumPW(0,1);
    //KittyToyI_SI.SetNumPW(1,1);
    KittyToyI_SI.SetSpin(0,0);
    //KittyToyI_SI.SetSpin(1,0);
    KittyToyI_SI.SetChannelWeight(0,1.);
    //KittyToyI_SI.SetChannelWeight(1,1./2.);

    KittyToyI_SI.SetShortRangePotential(0,0,DoubleGaussSum,POT_PARS_I1);
    //KittyToyI_SI.SetShortRangePotential(1,0,DoubleGaussSum,POT_PARS_I0);
    KittyToyI_SI.SetMaxNumThreads(4);
    KittyToyI_SI.SetEpsilonConv(1e-8);
    KittyToyI_SI.SetEpsilonProp(1e-8);
    KittyToyI_SI.KillTheCat();

    DLM_Histo<complex<double>>*** ExternalWF_0=NULL;
    CATS KittyHaide;
    KittyHaide.SetMomBins(30,5,305);
    KittyHaide.SetUseAnalyticSource(true);
    KittyHaide.SetAnaSource(GaussSource,cPars);
    ExternalWF_0 = InitHaidenbauerKaonPlus("/home/dmihaylov/CernBox/CATS_potentials/Haidenbauer/KaonPlus_SIonly/",KittyHaide,0);
    KittyHaide.SetExternalWaveFunction(0,0,ExternalWF_0[0][0][0],ExternalWF_0[1][0][0]);
    KittyHaide.KillTheCat();

    DLM_Histo<complex<double>>*** ExternalWF_1=NULL;
    CATS KittyHaideGamow;
    KittyHaideGamow.SetMomBins(30,5,305);
    KittyHaideGamow.SetUseAnalyticSource(true);
    KittyHaideGamow.SetAnaSource(GaussSource,cPars);
    ExternalWF_1 = InitHaidenbauerKaonPlus("/home/dmihaylov/CernBox/CATS_potentials/Haidenbauer/KaonPlus10MeV/",KittyHaideGamow,1);
    KittyHaideGamow.SetExternalWaveFunction(0,0,ExternalWF_1[0][0][0],ExternalWF_1[1][0][0]);
    KittyHaideGamow.KillTheCat();

    TGraph grCoulomb;
    grCoulomb.SetName("grCoulomb");
    grCoulomb.Set(NumMomBins);
    for(unsigned uMom=0; uMom<NumMomBins; uMom++){
        grCoulomb.SetPoint(uMom,KittyCoulomb.GetMomentum(uMom),KittyCoulomb.GetCorrFun(uMom));
    }

    TGraph grToyI;
    grToyI.SetName("grToyI");
    grToyI.Set(NumMomBins);
    for(unsigned uMom=0; uMom<NumMomBins; uMom++){
        grToyI.SetPoint(uMom,KittyToyI.GetMomentum(uMom),KittyToyI.GetCorrFun(uMom));
    }

    TGraph grToyI_SI;
    grToyI_SI.SetName("grToyI_SI");
    grToyI_SI.Set(NumMomBins);
    for(unsigned uMom=0; uMom<NumMomBins; uMom++){
        grToyI_SI.SetPoint(uMom,KittyToyI_SI.GetMomentum(uMom),KittyToyI_SI.GetCorrFun(uMom));
    }

    TGraph grHaideGamow;
    grHaideGamow.SetName("grHaideGamow");
    grHaideGamow.Set(30);
    for(unsigned uMom=0; uMom<30; uMom++){
        grHaideGamow.SetPoint(uMom,KittyHaideGamow.GetMomentum(uMom),KittyHaideGamow.GetCorrFun(uMom));
    }

    TGraph grHaide;
    grHaide.SetName("grHaide");
    grHaide.Set(30);
    for(unsigned uMom=0; uMom<30; uMom++){
        grHaide.SetPoint(uMom,KittyHaide.GetMomentum(uMom),KittyHaide.GetCorrFun(uMom));
    }

    TFile fOutput(OutputDir+"fOutput.root","recreate");
    grCoulomb.Write();
    grToyI.Write();
    grToyI_SI.Write();
    grHaideGamow.Write();
    grHaide.Write();

    delete [] SOL_POT_I0;
    delete [] SOL_POT_I1;
}

void Toy_pKplus_2(){
    const unsigned NumMomBins = 31;
    const double kMin = 0;
    const double kMax = 310;
    const double Mass_p = 938.272;
    const double Mass_KaonCh = 493.677;

    const double min_f0 = -0.390;
    const double max_f0 = -0.210;
    const unsigned NumSteps = 19;
    const double f0_step = (max_f0-min_f0)/double(NumSteps-1);

    const TString OutputDir = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/KaonProton/Toy_pKplus/";
    CATSparameters cPars(CATSparameters::tSource,1,true);
    cPars.SetParameter(0,1.2);

    DLM_Histo<complex<double>>*** ExternalWF_0=NULL;
    CATS KittyHaide;
    KittyHaide.SetMomBins(30,5,305);
    KittyHaide.SetUseAnalyticSource(true);
    KittyHaide.SetAnaSource(GaussSource,cPars);
    ExternalWF_0 = InitHaidenbauerKaonPlus("/home/dmihaylov/CernBox/CATS_potentials/Haidenbauer/KaonPlus_SIonly/",KittyHaide,0);
    KittyHaide.SetExternalWaveFunction(0,0,ExternalWF_0[0][0][0],ExternalWF_0[1][0][0]);
    KittyHaide.KillTheCat();

    DLM_Histo<complex<double>>*** ExternalWF_1=NULL;
    CATS KittyHaideGamow;
    KittyHaideGamow.SetMomBins(30,5,305);
    KittyHaideGamow.SetUseAnalyticSource(true);
    KittyHaideGamow.SetAnaSource(GaussSource,cPars);
    ExternalWF_1 = InitHaidenbauerKaonPlus("/home/dmihaylov/CernBox/CATS_potentials/Haidenbauer/KaonPlus10MeV/",KittyHaideGamow,1);
    KittyHaideGamow.SetExternalWaveFunction(0,0,ExternalWF_1[0][0][0],ExternalWF_1[1][0][0]);
    KittyHaideGamow.KillTheCat();

    TGraph grHaideGamow;
    grHaideGamow.SetName("grHaideGamow");
    grHaideGamow.Set(30);
    for(unsigned uMom=0; uMom<30; uMom++){
        grHaideGamow.SetPoint(uMom,KittyHaideGamow.GetMomentum(uMom),KittyHaideGamow.GetCorrFun(uMom));
    }

    TGraph grHaide;
    grHaide.SetName("grHaide");
    grHaide.Set(30);
    for(unsigned uMom=0; uMom<30; uMom++){
        grHaide.SetPoint(uMom,KittyHaide.GetMomentum(uMom),KittyHaide.GetCorrFun(uMom));
    }

    CATS KittyToyI;
    KittyToyI.SetMomBins(NumMomBins,kMin,kMax);
    KittyToyI.SetAnaSource(GaussSource, cPars);
    KittyToyI.SetUseAnalyticSource(true);
    KittyToyI.SetMomentumDependentSource(false);
    KittyToyI.SetThetaDependentSource(false);
    //should you include in the result any bins, where the Schroedinger solver failed
    KittyToyI.SetExcludeFailedBins(false);
    KittyToyI.SetQ1Q2(1);
    KittyToyI.SetPdgId(2212, 321);
    KittyToyI.SetRedMass( (Mass_p*Mass_KaonCh)/(Mass_p+Mass_KaonCh) );
    KittyToyI.SetNumChannels(1);
    KittyToyI.SetNumPW(0,1);
    KittyToyI.SetSpin(0,0);
    KittyToyI.SetChannelWeight(0,1.);

    CATS KittyToyI_SI;
    KittyToyI_SI.SetMomBins(NumMomBins,kMin,kMax);
    KittyToyI_SI.SetAnaSource(GaussSource, cPars);
    KittyToyI_SI.SetUseAnalyticSource(true);
    KittyToyI_SI.SetMomentumDependentSource(false);
    KittyToyI_SI.SetThetaDependentSource(false);
    //should you include in the result any bins, where the Schroedinger solver failed
    KittyToyI_SI.SetExcludeFailedBins(false);
    KittyToyI_SI.SetQ1Q2(0);
    KittyToyI_SI.SetPdgId(2212, 321);
    KittyToyI_SI.SetRedMass( (Mass_p*Mass_KaonCh)/(Mass_p+Mass_KaonCh) );
    KittyToyI_SI.SetNumChannels(1);
    KittyToyI_SI.SetNumPW(0,1);
    KittyToyI_SI.SetSpin(0,0);
    KittyToyI_SI.SetChannelWeight(0,1.);

    TGraph* grToyI = new TGraph[NumSteps];
    TGraph* grToyI_SI = new TGraph[NumSteps];

    for(unsigned uStep=0; uStep<NumSteps; uStep++){
        double* SOL_POT_I1 = FindPotentialPars("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/KaonProton/ScatParsFromRandPotential/Test/N11000.root",min_f0+f0_step*double(uStep),0.373,0);
        printf("%.3f, %.3f: V0=%.1f, mu0=%.3f, V1=%.1f, mu1=%.3f\n",
               SOL_POT_I1[0],SOL_POT_I1[1],SOL_POT_I1[2],SOL_POT_I1[3],SOL_POT_I1[4],SOL_POT_I1[5]);
        CATSparameters POT_PARS_I1(CATSparameters::tPotential,4,true);
        POT_PARS_I1.SetParameters(&SOL_POT_I1[2]);

        KittyToyI.SetShortRangePotential(0,0,DoubleGaussSum,POT_PARS_I1);
        KittyToyI.SetMaxNumThreads(4);
        KittyToyI.SetEpsilonConv(1e-8);
        KittyToyI.SetEpsilonProp(1e-8);
        KittyToyI.SetNotifications(CATS::nError);
        KittyToyI.KillTheCat();

        KittyToyI_SI.SetShortRangePotential(0,0,DoubleGaussSum,POT_PARS_I1);
        KittyToyI_SI.SetMaxNumThreads(4);
        KittyToyI_SI.SetEpsilonConv(1e-8);
        KittyToyI_SI.SetEpsilonProp(1e-8);
        KittyToyI_SI.SetNotifications(CATS::nError);
        KittyToyI_SI.KillTheCat();

        grToyI[uStep].SetName(TString::Format("grToyI_%.2f_%.2f",SOL_POT_I1[0],SOL_POT_I1[1]));
        grToyI[uStep].Set(NumMomBins);
        for(unsigned uMom=0; uMom<NumMomBins; uMom++){
            grToyI[uStep].SetPoint(uMom,KittyToyI.GetMomentum(uMom),KittyToyI.GetCorrFun(uMom));
        }

        grToyI_SI[uStep].SetName(TString::Format("grToyI_SI_%.2f_%.2f",SOL_POT_I1[0],SOL_POT_I1[1]));
        grToyI_SI[uStep].Set(NumMomBins);
        for(unsigned uMom=0; uMom<NumMomBins; uMom++){
            grToyI_SI[uStep].SetPoint(uMom,KittyToyI_SI.GetMomentum(uMom),KittyToyI_SI.GetCorrFun(uMom));
        }

        delete [] SOL_POT_I1;
    }

    TFile fOutput(OutputDir+"fOutput2.root","recreate");
    grHaideGamow.Write();
    grHaide.Write();
    for(unsigned uStep=0; uStep<NumSteps; uStep++){
        grToyI[uStep].Write();
        grToyI_SI[uStep].Write();
    }

}

//femto convention
double* FindPotentialPars(const TString& InputFileName, const double& ScatLen, const double& EffRan, const unsigned WhichSol){
    TFile fInput(InputFileName,"read");
    TNtuple* ntuple=NULL;
    ntuple = (TNtuple*)fInput.Get("Table");
    if(!ntuple){
        printf("Problem with the ntuple!\n");
        abort();
    }
    unsigned NumNtEntries = ntuple->GetEntries();
    Float_t par_V0;
    Float_t par_V1;
    Float_t par_V2;
    Float_t par_mu0;
    Float_t par_mu1;
    Float_t par_mu2;
    Float_t par_f0Inv;
    Float_t par_d0;

    const Float_t DevLimit = 0.02;
    Float_t* Deviation = new Float_t[NumNtEntries];
    Float_t* Solution_V0 = new Float_t[NumNtEntries];
    Float_t* Solution_V1 = new Float_t[NumNtEntries];
    Float_t* Solution_V2 = new Float_t[NumNtEntries];
    Float_t* Solution_mu0 = new Float_t[NumNtEntries];
    Float_t* Solution_mu1 = new Float_t[NumNtEntries];
    Float_t* Solution_mu2 = new Float_t[NumNtEntries];
    Float_t* Solution_f0Inv = new Float_t[NumNtEntries];
    Float_t* Solution_d0 = new Float_t[NumNtEntries];

    ntuple->SetBranchAddress("V0",&par_V0);
    ntuple->SetBranchAddress("V1",&par_V1);
    ntuple->SetBranchAddress("V2",&par_V2);
    ntuple->SetBranchAddress("mu0",&par_mu0);
    ntuple->SetBranchAddress("mu1",&par_mu1);
    ntuple->SetBranchAddress("mu2",&par_mu2);
    ntuple->SetBranchAddress("f0Inv",&par_f0Inv);
    ntuple->SetBranchAddress("d0",&par_d0);

    unsigned PossibleSolutions = 0;

    for(unsigned uEntry=0; uEntry<NumNtEntries; uEntry++){
        ntuple->GetEntry(uEntry);
        //we put twice the weight to the scattering length
        Float_t Delta_f0 = (ScatLen-1./par_f0Inv);
        Float_t Delta_d0 = (EffRan-par_d0);
        Deviation[PossibleSolutions] = sqrt(4*Delta_f0*Delta_f0+Delta_d0*Delta_d0);
        if(Deviation[PossibleSolutions]>DevLimit) continue;
        Solution_V0[PossibleSolutions] = par_V0;
        Solution_V1[PossibleSolutions] = par_V1;
        Solution_V2[PossibleSolutions] = par_V2;
        Solution_mu0[PossibleSolutions] = par_mu0;
        Solution_mu1[PossibleSolutions] = par_mu1;
        Solution_mu2[PossibleSolutions] = par_mu2;
        Solution_f0Inv[PossibleSolutions] = par_f0Inv;
        Solution_d0[PossibleSolutions] = par_d0;
        PossibleSolutions++;
    }

    DLM_Sort < Float_t, unsigned > SortTool;
    SortTool.SetData(Deviation,PossibleSolutions);
    SortTool.MergeSort();
    SortTool.GetSortedData(Deviation,Deviation);

    Float_t* Temp;
    Temp = new Float_t[PossibleSolutions];

    for(unsigned uSol=0; uSol<PossibleSolutions; uSol++){
        Temp[uSol] = Solution_V0[SortTool.GetKey()[uSol]];
    }
    for(unsigned uEl=0; uEl<PossibleSolutions; uEl++){
        Solution_V0[uEl] = Temp[uEl];
    }

    for(unsigned uSol=0; uSol<PossibleSolutions; uSol++){
        Temp[uSol] = Solution_V1[SortTool.GetKey()[uSol]];
    }
    for(unsigned uEl=0; uEl<PossibleSolutions; uEl++){
        Solution_V1[uEl] = Temp[uEl];
    }

    for(unsigned uSol=0; uSol<PossibleSolutions; uSol++){
        Temp[uSol] = Solution_V2[SortTool.GetKey()[uSol]];
    }
    for(unsigned uEl=0; uEl<PossibleSolutions; uEl++){
        Solution_V2[uEl] = Temp[uEl];
    }

    for(unsigned uSol=0; uSol<PossibleSolutions; uSol++){
        Temp[uSol] = Solution_mu0[SortTool.GetKey()[uSol]];
    }
    for(unsigned uEl=0; uEl<PossibleSolutions; uEl++){
        Solution_mu0[uEl] = Temp[uEl];
    }

    for(unsigned uSol=0; uSol<PossibleSolutions; uSol++){
        Temp[uSol] = Solution_mu1[SortTool.GetKey()[uSol]];
    }
    for(unsigned uEl=0; uEl<PossibleSolutions; uEl++){
        Solution_mu1[uEl] = Temp[uEl];
    }

    for(unsigned uSol=0; uSol<PossibleSolutions; uSol++){
        Temp[uSol] = Solution_mu2[SortTool.GetKey()[uSol]];
    }
    for(unsigned uEl=0; uEl<PossibleSolutions; uEl++){
        Solution_mu2[uEl] = Temp[uEl];
    }

    for(unsigned uSol=0; uSol<PossibleSolutions; uSol++){
        Temp[uSol] = Solution_f0Inv[SortTool.GetKey()[uSol]];
    }
    for(unsigned uEl=0; uEl<PossibleSolutions; uEl++){
        Solution_f0Inv[uEl] = Temp[uEl];
    }

    for(unsigned uSol=0; uSol<PossibleSolutions; uSol++){
        Temp[uSol] = Solution_d0[SortTool.GetKey()[uSol]];
    }
    for(unsigned uEl=0; uEl<PossibleSolutions; uEl++){
        Solution_d0[uEl] = Temp[uEl];
    }

    delete [] Temp;


//    for(unsigned uSol=0; uSol<PossibleSolutions; uSol++){
//        printf("%u: %f f0=%.3f; d0=%.3f; V0=%.2f; mu0=%.4f; V1=%.2f; mu1=%.4f\n",uSol,Deviation[uSol],
//               1./Solution_f0Inv[uSol],Solution_d0[uSol],Solution_V0[uSol],Solution_mu0[uSol],Solution_V1[uSol],Solution_mu1[uSol]);
//    }

    double* SOLUTION = NULL;
    if(PossibleSolutions){
        SOLUTION=new double [8];
        SOLUTION[0]=1./Solution_f0Inv[WhichSol];
        SOLUTION[1]=Solution_d0[WhichSol];
        SOLUTION[2]=Solution_V0[WhichSol];
        SOLUTION[3]=Solution_mu0[WhichSol];
        SOLUTION[4]=Solution_V1[WhichSol];
        SOLUTION[5]=Solution_mu1[WhichSol];
        SOLUTION[6]=Solution_V2[WhichSol];
        SOLUTION[7]=Solution_mu2[WhichSol];
    }


    delete [] Deviation;
    delete [] Solution_V0;
    delete [] Solution_V1;
    delete [] Solution_V2;
    delete [] Solution_mu0;
    delete [] Solution_mu1;
    delete [] Solution_mu2;
    delete [] Solution_f0Inv;
    delete [] Solution_d0;

    return SOLUTION;
}

void ScatParsFromRandPotential(const TString OutputFolder,
                               const unsigned RandomSeed, const unsigned& NumIter,
                               const double& V0_min, const double& V0_max,
                               const double& mu0_min, const double& mu0_max,
                               const double& V1_min, const double& V1_max,
                               const double& mu1_min, const double& mu1_max,
                               const double& V2_min, const double& V2_max,
                               const double& mu2_min, const double& mu2_max){



    const double Mass_p = 938.272;
    const double Mass_KaonCh = 493.677;

    DLM_Timer TIMER;
    TIMER.Start();

    TFile* OutputFile = new TFile(OutputFolder+TString::Format("N%u_RS%u_Table.root",NumIter,RandomSeed),"RECREATE");
    TNtuple* Table = new TNtuple("Table","Table","V0:mu0:V1:mu1:V2:mu2:f0Inv:d0");

    TRandom3 RanGen(RandomSeed);
    double Ampl0,Ampl1,Ampl2,mu0,mu1,mu2;

    //Kitty.KillTheCat();

    CATSparameters POT_PARS_I1(CATSparameters::tPotential,6,true);
    CATSparameters cPars(CATSparameters::tSource,1,true);
    cPars.SetParameter(0,1.2);
    //!Femto sign convention

    CATS Kitty;
    const double kFrom = 4;
    const double kTo = 68;
    const unsigned kNumStep = 32;
    const double kStep = (kTo-kFrom)/(kNumStep-1);
    //double SorPars[4] = {0,0,0,2};

    Kitty.SetUseAnalyticSource(true);
    Kitty.SetAnaSource(GaussSource, cPars);
    Kitty.SetMaxRad(100);
    Kitty.SetMaxRho(50);
    Kitty.SetEpsilonProp(2e-9);
    Kitty.SetEpsilonConv(1e-8);
    Kitty.SetExcludeFailedBins(false);
    Kitty.SetMomBins(kNumStep,kFrom,kTo);
    Kitty.SetNumChannels(1);
    Kitty.SetNumPW(0,1);
    Kitty.SetSpin(0,0);
    Kitty.SetChannelWeight(0, 1);
    Kitty.SetQ1Q2(0);
    Kitty.SetPdgId(2212, 321);
    Kitty.SetRedMass( (Mass_p*Mass_KaonCh)/(Mass_p+Mass_KaonCh) );
    //Kitty.SetMaxNumThreads(4);
    Kitty.KillTheCat();
    Kitty.SetNotifications(CATS::nSilent);

    printf("Ready to iterate!\n");

    const double rMin = 0;
    const double rMax = 4;
    const int rSteps = 256;

    for(unsigned uIter=0; uIter<NumIter; uIter++){

        TH1F* hPhaseShift0 = new TH1F(TString::Format("hPhaseShift0_%u",uIter), TString::Format("hPhaseShift0_%u",uIter), kNumStep, kFrom-0.5*kStep, kTo+0.5*kStep);
        TF1* fPhaseShift0 = new TF1(TString::Format("fPhaseShift0%u",uIter), "0.5*[1]/197.327*x*x+[0]*197.327", kFrom, kTo);
        TH1F* hPotential = new TH1F(TString::Format("hPotential_%u",uIter), TString::Format("hPotential_%u",uIter), rSteps, rMin, rMax);

        if(V0_min<=V0_max){
            Ampl0 = RanGen.Uniform(V0_min,V0_max);
            if(mu0_min<=mu0_max){
                mu0 = RanGen.Uniform(mu0_min,mu0_max);
            }
            else{
                mu0 = 1;
            }
        }
        else{
            Ampl0=0;
            mu0=1;
        }

        if(V1_min<=V1_max){
            Ampl1 = RanGen.Uniform(V1_min,V1_max);
            if(mu1_min<=mu1_max){
                mu1 = RanGen.Uniform(mu1_min,mu1_max);
            }
            else{
                mu1 = 1;
            }
        }
        else{
            Ampl1=0;
            mu1=1;
        }

        if(V2_min<=V2_max){
            Ampl2 = RanGen.Uniform(V2_min,V2_max);
            if(mu2_min<=mu2_max){
                mu2 = RanGen.Uniform(mu2_min,mu2_max);
            }
            else{
                mu2 = 1;
            }
        }
        else{
            Ampl2=0;
            mu2=1;
        }

        POT_PARS_I1.SetParameter(0,Ampl0);
        POT_PARS_I1.SetParameter(1,mu0);
        POT_PARS_I1.SetParameter(2,Ampl1);
        POT_PARS_I1.SetParameter(3,mu1);
        POT_PARS_I1.SetParameter(4,Ampl2);
        POT_PARS_I1.SetParameter(5,mu2);
        Kitty.SetShortRangePotential(0,0,GaussExpSum,POT_PARS_I1);
        //
        Kitty.ComputeTheRadialWaveFunction();

        for(unsigned uBin=0; uBin<kNumStep; uBin++){
            hPhaseShift0->SetBinContent(uBin+1,Kitty.GetMomentum(uBin)/tan(Kitty.GetPhaseShift(uBin,0,0)));
            hPhaseShift0->SetBinError(uBin+1, hPhaseShift0->GetBinContent(uBin+1)*0.0001);
        }
        for(unsigned uBin=0; uBin<rSteps; uBin++){
            hPotential->SetBinContent(uBin+1,Kitty.EvaluateThePotential(0,0,(kTo+kFrom)*0.5,hPotential->GetBinCenter(uBin+1)));
        }
        fPhaseShift0->SetParameter(0,0);
        fPhaseShift0->SetParLimits(0,-1000,1000);
        fPhaseShift0->SetParameter(1,0);
        fPhaseShift0->SetParLimits(1,-1000,1000);
        hPhaseShift0->Fit(fPhaseShift0, "Q, S, N, R, M");
        printf("%.1f*exp(-r^2/%.2f)",Ampl0,mu0);
        if(Ampl1>0) printf("+%.1f*exp(-r^2/%.2f)",Ampl1,mu1);
        else printf("-%.1f*exp(-r^2/%.2f)",-Ampl1,mu1);
        if(Ampl2>0) printf("+%.1f*exp(-r^2/%.2f)",Ampl2,mu2);
        else printf("-%.1f*exp(-r^2/%.2f)",-Ampl2,mu2);
        printf("\n 1/f0 = %.2f",fPhaseShift0->GetParameter(0));
        printf(" d0 = %.2f\n",fPhaseShift0->GetParameter(1));
        Float_t CONTAINER[8];
        CONTAINER[0] = Ampl0;
        CONTAINER[1] = mu0;
        CONTAINER[2] = Ampl1;
        CONTAINER[3] = mu1;
        CONTAINER[4] = Ampl2;
        CONTAINER[5] = mu2;
        CONTAINER[6] = fPhaseShift0->GetParameter(0);
        CONTAINER[7] = fPhaseShift0->GetParameter(1);
        Table->Fill(CONTAINER);

        //hPhaseShift0->Write();
        //fPhaseShift0->Write();
        //hPotential->Write();

        delete hPhaseShift0;
        delete fPhaseShift0;
        delete hPotential;

        long long ExecTime = TIMER.Stop();
        //break after 2 h 15 mins
        if(double(ExecTime)*1e-6/60.>30) break;
    }
    Table->Write();



    delete Table;
    delete OutputFile;
    //Table->Draw("d0:f0Inv>>h2(128,-1000,1000,96,-1000,1000)","","colz")
}

void TestKyoto2019(const double& RADIUS){
    const unsigned NumChannels = 6;

    //(tSource or tPotential, Number of parameters, keep it true)
    CATSparameters cPars(CATSparameters::tSource,1,true);
    cPars.SetParameter(0,RADIUS);

    CATS KittyTemp;
    //KittyTemp.SetMomBins(150,1,301);
    DLM_Histo<complex<double>>*** BINNING=Init_pKminus_Kyoto2019("/home/dmihaylov/CernBox/CATS_potentials/Tetsuo/Kyoto2019/",KittyTemp,0);
    unsigned NumMomBins = BINNING[0][0][0].GetNbins(0);
    double* MomBins = BINNING[0][0][0].GetBinRange(0);
    double* MomBinsCenter = BINNING[0][0][0].GetBinCenters(0);
    printf("NumMomBins=%u\n",NumMomBins);
    for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
        printf("%u: %3.f %.3f %.3f\n",uMomBin,MomBins[uMomBin],MomBinsCenter[uMomBin],MomBins[uMomBin+1]);
    }

    CATS KittyStrong;
    KittyStrong.SetMomBins(NumMomBins,MomBins,MomBinsCenter);
    KittyStrong.SetAnaSource(GaussSource, cPars);
    KittyStrong.SetUseAnalyticSource(true);
    DLM_Histo<complex<double>>*** ExternalWF_Strong=Init_pKminus_Kyoto2019("/home/dmihaylov/CernBox/CATS_potentials/Tetsuo/Kyoto2019/",KittyStrong,0);
    //printf("NumMomBins=%u\n",ExternalWF_Strong[0][0][0].GetNbins());
    for(unsigned uCh=0; uCh<NumChannels; uCh++){
        KittyStrong.SetExternalWaveFunction(uCh,0,ExternalWF_Strong[0][uCh][0],ExternalWF_Strong[1][uCh][0]);
    }
    KittyStrong.KillTheCat();

    CATS KittyGamow;
    KittyGamow.SetMomBins(NumMomBins,MomBins,MomBinsCenter);
    KittyGamow.SetAnaSource(GaussSource, cPars);
    KittyGamow.SetUseAnalyticSource(true);
    DLM_Histo<complex<double>>*** ExternalWF_Gamow=Init_pKminus_Kyoto2019("/home/dmihaylov/CernBox/CATS_potentials/Tetsuo/Kyoto2019/",KittyGamow,0);
    for(unsigned uCh=0; uCh<NumChannels; uCh++){
        KittyGamow.SetExternalWaveFunction(uCh,0,ExternalWF_Gamow[0][uCh][0],ExternalWF_Gamow[1][uCh][0]);
    }
    KittyGamow.SetQ1Q2(-1);
    KittyGamow.SetGamow(true);
    KittyGamow.KillTheCat();

    CATS KittyFull;
    KittyFull.SetMomBins(NumMomBins,MomBins,MomBinsCenter);
    KittyFull.SetAnaSource(GaussSource, cPars);
    KittyFull.SetUseAnalyticSource(true);
    DLM_Histo<complex<double>>*** ExternalWF_Full=Init_pKminus_Kyoto2019("/home/dmihaylov/CernBox/CATS_potentials/Tetsuo/Kyoto2019/",KittyFull,1);
    for(unsigned uCh=0; uCh<NumChannels; uCh++){
        KittyFull.SetExternalWaveFunction(uCh,0,ExternalWF_Full[0][uCh][0],ExternalWF_Full[1][uCh][0]);
    }
    //for channels 1-5 these are the omega weights
    //KittyFull.SetChannelWeight(WHICH_CHANNEL, WEIGHT);
    //if you change the weights, kill the cat again
    KittyFull.KillTheCat();

    TGraph* gStrong = new TGraph [NumChannels+1];
    TGraph* gGamow = new TGraph [NumChannels+1];
    TGraph* gFull = new TGraph [NumChannels+1];
    TH1F** hWfStrong = new TH1F* [NumChannels+1];
    TH1F** hWfFull = new TH1F* [NumChannels+1];
    for(unsigned uCh=0; uCh<=NumChannels; uCh++){
        hWfStrong[uCh] = new TH1F(TString::Format("hWfStrong_%u",uCh),TString::Format("hWfStrong_%u",uCh),16384,0,128);
        hWfFull[uCh] = new TH1F(TString::Format("hWfFull_%u",uCh),TString::Format("hWfFull_%u",uCh),16384,0,128);
    }

    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        gStrong[NumChannels].SetPoint(uBin,KittyStrong.GetMomentum(uBin),KittyStrong.GetCorrFun(uBin));
        gGamow[NumChannels].SetPoint(uBin,KittyGamow.GetMomentum(uBin),KittyGamow.GetCorrFun(uBin));
        gFull[NumChannels].SetPoint(uBin,KittyFull.GetMomentum(uBin),KittyFull.GetCorrFun(uBin));
    }

    KittyStrong.SetNotifications(CATS::nError);
    KittyGamow.SetNotifications(CATS::nError);
    KittyFull.SetNotifications(CATS::nError);

    for(unsigned uCh=0; uCh<NumChannels; uCh++){
        for(unsigned uCh2=0; uCh2<NumChannels; uCh2++){
            KittyStrong.SetChannelWeight(uCh2,0);
            KittyGamow.SetChannelWeight(uCh2,0);
            KittyFull.SetChannelWeight(uCh2,0);
        }
        KittyStrong.SetChannelWeight(uCh,1);
        KittyGamow.SetChannelWeight(uCh,1);
        KittyFull.SetChannelWeight(uCh,1);
        KittyStrong.KillTheCat();
        KittyGamow.KillTheCat();
        KittyFull.KillTheCat();
        for(unsigned uBin=0; uBin<NumMomBins; uBin++){
            gStrong[uCh].SetPoint(uBin,KittyStrong.GetMomentum(uBin),KittyStrong.GetCorrFun(uBin));
            gGamow[uCh].SetPoint(uBin,KittyGamow.GetMomentum(uBin),KittyGamow.GetCorrFun(uBin));
            gFull[uCh].SetPoint(uBin,KittyFull.GetMomentum(uBin),KittyFull.GetCorrFun(uBin));
        }
        for(unsigned uBin=1; uBin<=hWfStrong[uCh]->GetNbinsX(); uBin++){
            hWfStrong[uCh]->SetBinContent(uBin,std::abs(KittyStrong.EvalRadialWaveFunction(40,uCh,0,hWfStrong[uCh]->GetBinCenter(uBin),true)));
            hWfFull[uCh]->SetBinContent(uBin,std::abs(KittyFull.EvalRadialWaveFunction(40,uCh,0,hWfFull[uCh]->GetBinCenter(uBin),true)));
        }
    }

    TFile fOutput(TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/TestKyoto2019/fOutput_%.2f.root",
                                  RADIUS),"recreate");

    for(unsigned uCh=0; uCh<=NumChannels; uCh++){
        TString AddOn;
        if(uCh<NumChannels) AddOn = TString::Format("%u",uCh);
        else AddOn = "ALLw1";
        gStrong[uCh].SetName(TString::Format("gStrong_%s",AddOn.Data()));
        gGamow[uCh].SetName(TString::Format("gGamow_%s",AddOn.Data()));
        gFull[uCh].SetName(TString::Format("gFull_%s",AddOn.Data()));
        gStrong[uCh].Write();
        //gGamow[uCh].Write();
        gFull[uCh].Write();
        //hWfStrong[uCh]->Write();
        //hWfFull[uCh]->Write();
    }

    FILE * CkFile;
    CkFile = fopen (TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/TestKyoto2019/CkProtonKaon.txt"),"w");
    fprintf (CkFile, "%10s%10s%10s%10s\n","k* (MeV)","C_S(k*)","C_G(k*)","C_C(k*)");
    double kVal_S,CkVal_S;
	double kVal_G,CkVal_G;
	double kVal_C,CkVal_C;
	for(unsigned uBin=0; uBin<NumMomBins; uBin++){
		gStrong->GetPoint(uBin,kVal_S,CkVal_S);
		gGamow->GetPoint(uBin,kVal_G,CkVal_G);
		gFull->GetPoint(uBin,kVal_C,CkVal_C);
		fprintf (CkFile, "%10.1f%10.3f%10.3f%10.3f\n",kVal_S,CkVal_S,CkVal_G,CkVal_C);
	}
	fclose(CkFile);

    delete [] gStrong;
    delete [] gGamow;
    delete [] gFull;
    for(unsigned uCh=0; uCh<=NumChannels; uCh++){
        delete hWfStrong[uCh];
        delete hWfFull[uCh];
    }
    delete [] hWfStrong;
    delete [] hWfFull;
    delete [] MomBins;
}

void SmearTest_pKminus(){
    const unsigned NumMomBins = 100;
    const double MinMomentum = 0;
    const double MaxMomentum = 250;
    const double RADIUS = 1.0;
    TString OutputFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/KaonProton/SmearTest_pKminus/fOutput.root";

    //(tSource or tPotential, Number of parameters, keep it true)
    CATSparameters cPars(CATSparameters::tSource,1,true);
    cPars.SetParameter(0,RADIUS);

    CATS KittyFull;
    KittyFull.SetMomBins(NumMomBins,MinMomentum,MaxMomentum);
    KittyFull.SetAnaSource(GaussSource, cPars);
    KittyFull.SetUseAnalyticSource(true);
    DLM_Histo<complex<double>>*** ExternalWF_Full=Init_pKminus_Kyoto2019("/home/dmihaylov/CernBox/CATS_potentials/Tetsuo/Kyoto2019/",KittyFull,1);
    for(unsigned uCh=0; uCh<6; uCh++){
        KittyFull.SetExternalWaveFunction(uCh,0,ExternalWF_Full[0][uCh][0],ExternalWF_Full[1][uCh][0]);
    }
    //for channels 1-5 these are the omega weights
    //KittyFull.SetChannelWeight(WHICH_CHANNEL, WEIGHT);
    //if you change the weights, kill the cat again
    KittyFull.KillTheCat();

    //in essence, this is a custom made histogram class designed to work directly with CATS
    DLM_Ck CkHisto(KittyFull.GetNumSourcePars(),0,KittyFull);
    //fills all bins with the contents from KittyFull
    CkHisto.Update();

    //the next lines are getting the momentum-smear matrix from a file, saving a local copy to work with
    TString FileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/ResolutionMatrices/Sample6_MeV_compact.root";
    TString HistoName = "hSigmaMeV_Proton_Proton";
    ///ROOT SUCKS!
    //so we need to copy our histogram, as else we lose it when we delete the file
    //and we need to change to the "central" root directory, as else hResolution will also be lost
    //and we need to play with the name a little bit, else we are fucked!
    TFile* FileROOT = new TFile(FileName, "read");
    TH2F* histoTemp = (TH2F*)FileROOT->Get(HistoName);
    if(!histoTemp){printf("\033[1;31mERROR:\033[0m The histoTemp '%s' if file '%s' does not exist\n",HistoName.Data(),FileName.Data());}
    TString Name = histoTemp->GetName();
    gROOT->cd();
    TH2F *hResolution = (TH2F*)histoTemp->Clone("hResolution");
    delete FileROOT;
    hResolution->SetName(Name);

    //this class is used to apply feed-down and resolution corrections onto the theory curve
    //in this example we only use it to smear the Ck for the momentum resolution
    //arguments: "Arbitrary name of this object", number of feed-down contributions, Ck histogram, resolution matrix
    DLM_CkDecomposition CkCorrected("pKminus",0,CkHisto,hResolution);
    //applies the smearing and computes the corresponding histogram
    CkCorrected.Update();

    TFile fOutput(OutputFileName,"recreate");
    //finally we save the output, both the theoretical and smeared correlations
    TGraph gCkTheory;
    gCkTheory.SetName("gCkTheory");
    TGraph gCkCorrected;
    gCkCorrected.SetName("gCkCorrected");

    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        double Momentum;
        double Value_CkTh;
        double Value_CkCorr;
        Momentum = CkHisto.GetBinCenter(0,uBin);
        Value_CkTh = CkHisto.Eval(Momentum);
        Value_CkCorr = CkCorrected.EvalCk(Momentum);
        gCkTheory.SetPoint(uBin,Momentum,Value_CkTh);
        gCkCorrected.SetPoint(uBin,Momentum,Value_CkCorr);
    }

    gCkTheory.Write();
    gCkCorrected.Write();

}

int KAONPROTON_MAIN(int argc, char *argv[]){

    SmearTest_pKminus();

    //TestKyoto2019(1.2);
    //for(double rad=1; rad<=7.5; rad+=0.5){
    //    TestKyoto2019(rad);
    //}

    //Toy_pKplus();
    //Toy_pKplus_2();
    //ScatParsFromRandPotential(argv[1],atoi(argv[2]),atoi(argv[3]),atof(argv[4]),atof(argv[5]),
    //                          atof(argv[6]),atof(argv[7]),atof(argv[8]),atof(argv[9]),atof(argv[10]),atof(argv[11]),
    //                          atof(argv[12]),atof(argv[13]),atof(argv[14]),atof(argv[15]));
    //ScatParsFromRandPotential("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/KaonProton/ScatParsFromRandPotential/Test/",
    //                          1,10,0,2500,0,1,0,2500,0,1);

    //FindPotentialPars("/home/dmihaylov/Temp/pKplus_1M.root",-0.316,0.373);

    //AnitKaonProton_PhaseShiftsAndWaveFunctions();
    //KminProton_ImagWF();
    //TestHaideJuly2018();
    //TestHaideSiddharta(1.072);
    //TestHaideSiddharta(1.125);
    //TestHaideSiddharta(1.201);
    //TestHaideSiddharta(1.173);
    //TestHaideSiddharta(1.182);
    //TestHaideSiddharta(1.195);
    //GamowCheck();
    //GamowCheck2();
    return 0;
}
