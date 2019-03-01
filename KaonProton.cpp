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

#include "TString.h"
#include "TGraph.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TStyle.h"

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
    /*
    for(unsigned uBin=0; uBin<uNumBins; uBin++){
        grCk[0].GetPoint(uBin,xVal1,yVal1);
        grCk[1].GetPoint(uBin,xVal2,yVal2);
        gTest1.SetPoint(uBin,xVal1,(yVal1+yVal2));
    }
    */
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

int KAONPROTON_MAIN(int argc, char *argv[]){

    //AnitKaonProton_PhaseShiftsAndWaveFunctions();
    //KminProton_ImagWF();
    //TestHaideJuly2018();
    TestHaideSiddharta(1.072);
    TestHaideSiddharta(1.125);
    TestHaideSiddharta(1.201);
    TestHaideSiddharta(1.173);
    TestHaideSiddharta(1.182);
    TestHaideSiddharta(1.195);
    //GamowCheck();
    //GamowCheck2();
    return 0;
}
