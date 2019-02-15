#include "CATSpaper.h"
#include "CATStools.h"
//#include "DLM_CATSresults2.h"
#include "DLM_CkModels.h"
#include "DLM_WfModel.h"

#include "DLM_Potentials.h"
#include "DLM_Source.h"
#include "DLM_DrawingTools.h"
#include "DLM_SubPads.h"
#include "DLM_CkDecomposition.h"

#include "TH1F.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TStyle.h"

#include <vector>
#include <boost/algorithm/string.hpp>
#include <fstream>

using namespace boost;
using namespace std;

CATSpaper::CATSpaper(const TString& outdir):OutputDir(outdir),Mass_p(938.272),Mass_L(1115.683){
    //Kitty_pp = NULL;
    //Kitty_pL = NULL;
    ppPars1S0 = NULL;
    ppPars3P0 = NULL;
    ppPars3P1 = NULL;
    ppPars3P2 = NULL;
    ppSourcePars = NULL;
    pLamPars1S0 = NULL;
    pLamPars3S1 = NULL;
    pLamSourcePars = NULL;

    StandardKitty_pp();
    StandardKitty_pL();
}
CATSpaper::~CATSpaper(){
    //if(Kitty_pp) {delete Kitty_pp; Kitty_pp=NULL;}
    //if(Kitty_pL) {delete Kitty_pL; Kitty_pL=NULL;}

    if(ppPars1S0) {delete ppPars1S0; ppPars1S0=NULL;}
    if(ppPars3P0) {delete ppPars3P0; ppPars3P0=NULL;}
    if(ppPars3P1) {delete ppPars3P1; ppPars3P1=NULL;}
    if(ppPars3P2) {delete ppPars3P2; ppPars3P2=NULL;}
    if(ppSourcePars) {delete ppSourcePars; ppSourcePars=NULL;}
    if(pLamPars1S0) {delete pLamPars1S0; pLamPars1S0=NULL;}
    if(pLamPars3S1) {delete pLamPars3S1; pLamPars3S1=NULL;}
    if(pLamSourcePars) {delete pLamSourcePars; pLamSourcePars=NULL;}
}

void CATSpaper::StandardKitty_pp(){
    //if(Kitty_pp) {delete Kitty_pp; Kitty_pp=NULL;}
    //Kitty_pp = new CATS();

    const double RadLarge = 1.2;
    //const double RadSmall = 0.85;

    if(!ppPars1S0) ppPars1S0 = new CATSparameters(CATSparameters::tPotential,8,true);
    ppPars1S0->SetParameter(0,NN_AV18);
    ppPars1S0->SetParameter(1,v18_Coupled3P2);
    ppPars1S0->SetParameter(2,1);
    ppPars1S0->SetParameter(3,1);
    ppPars1S0->SetParameter(4,1);
    ppPars1S0->SetParameter(5,0);
    ppPars1S0->SetParameter(6,0);
    ppPars1S0->SetParameter(7,0);

    if(!ppPars3P0) ppPars3P0 = new CATSparameters(CATSparameters::tPotential,8,true);
    ppPars3P0->SetParameter(0,NN_AV18);
    ppPars3P0->SetParameter(1,v18_Coupled3P2);
    ppPars3P0->SetParameter(2,1);
    ppPars3P0->SetParameter(3,1);
    ppPars3P0->SetParameter(4,1);
    ppPars3P0->SetParameter(5,1);
    ppPars3P0->SetParameter(6,1);
    ppPars3P0->SetParameter(7,0);

    if(!ppPars3P1) ppPars3P1 = new CATSparameters(CATSparameters::tPotential,8,true);
    ppPars3P1->SetParameter(0,NN_AV18);
    ppPars3P1->SetParameter(1,v18_Coupled3P2);
    ppPars3P1->SetParameter(2,1);
    ppPars3P1->SetParameter(3,1);
    ppPars3P1->SetParameter(4,1);
    ppPars3P1->SetParameter(5,1);
    ppPars3P1->SetParameter(6,1);
    ppPars3P1->SetParameter(7,1);

    if(!ppPars3P2) ppPars3P2 = new CATSparameters(CATSparameters::tPotential,8,true);
    ppPars3P2->SetParameter(0,NN_AV18);
    ppPars3P2->SetParameter(1,v18_Coupled3P2);
    ppPars3P2->SetParameter(2,1);
    ppPars3P2->SetParameter(3,1);
    ppPars3P2->SetParameter(4,1);
    ppPars3P2->SetParameter(5,1);
    ppPars3P2->SetParameter(6,1);
    ppPars3P2->SetParameter(7,2);

    if(!ppSourcePars) ppSourcePars = new double [4];
    ppSourcePars[3] = RadLarge;

    NumMomBins_pp = 125;
    kMin_pp = 0;
    kMax_pp = 500;
/*
    Kitty_pp->SetUseAnalyticSource(true);
    Kitty_pp->SetAnaSource(GaussSource, ppSourcePars);

    Kitty_pp->SetExcludeFailedBins(false);
    Kitty_pp->SetMomBins(NumMomBins_pp,kMin_pp,kMax_pp);

    Kitty_pp->SetNumChannels(4);
    Kitty_pp->SetNumPW(0,2);
    Kitty_pp->SetNumPW(1,2);
    Kitty_pp->SetNumPW(2,2);
    Kitty_pp->SetNumPW(3,2);
    Kitty_pp->SetSpin(0,0);
    Kitty_pp->SetSpin(1,1);
    Kitty_pp->SetSpin(2,1);
    Kitty_pp->SetSpin(3,1);
    Kitty_pp->SetChannelWeight(0, 3./12.);
    Kitty_pp->SetChannelWeight(1, 1./12.);
    Kitty_pp->SetChannelWeight(2, 3./12.);
    Kitty_pp->SetChannelWeight(3, 5./12.);
    Kitty_pp->SetShortRangePotential(0,0,fDlmPot,ppPars1S0);
    Kitty_pp->SetShortRangePotential(1,1,fDlmPot,ppPars3P0);
    Kitty_pp->SetShortRangePotential(2,1,fDlmPot,ppPars3P1);
    Kitty_pp->SetShortRangePotential(3,1,fDlmPot,ppPars3P2);

    Kitty_pp->SetQ1Q2(1);
    Kitty_pp->SetPdgId(2212, 2212);

    Kitty_pp->SetRedMass( 0.5*Mass_p );
*/
}
void CATSpaper::StandardKitty_pL(){
    //if(Kitty_pL) {delete Kitty_pL; Kitty_pL=NULL;}

    const double RadLarge = 1.2;
    //const double RadSmall = 0.85;

    if(!pLamPars1S0) pLamPars1S0 = new CATSparameters(CATSparameters::tPotential,8,true);
    pLamPars1S0->SetParameter(0,pL_UsmaniOli);
    pLamPars1S0->SetParameter(1,0);
    pLamPars1S0->SetParameter(2,1);
    pLamPars1S0->SetParameter(3,1);
    pLamPars1S0->SetParameter(4,1);
    pLamPars1S0->SetParameter(5,0);
    pLamPars1S0->SetParameter(6,0);
    pLamPars1S0->SetParameter(7,0);

    if(!pLamPars3S1) pLamPars3S1 = new CATSparameters(CATSparameters::tPotential,8,true);
    pLamPars1S0->SetParameter(0,pL_UsmaniOli);
    pLamPars1S0->SetParameter(1,0);
    pLamPars1S0->SetParameter(2,1);
    pLamPars1S0->SetParameter(3,1);
    pLamPars1S0->SetParameter(4,1);
    pLamPars1S0->SetParameter(5,1);
    pLamPars1S0->SetParameter(6,0);
    pLamPars1S0->SetParameter(7,1);

    if(!pLamSourcePars) pLamSourcePars = new double [4];
    pLamSourcePars[3] = RadLarge;

    NumMomBins_pL = 25;
    kMin_pL = 0;
    kMax_pL = 500;
}

void CATSpaper::PhaseShiftPlot_pp(){
    const double Plot_kMin = 0;
    const double Plot_kMax = 420;
    const double NumPlotPts = 70;
    CATShisto<double> CatHist(NumPlotPts,Plot_kMin,Plot_kMax);

    const unsigned NumDataPts = 11;
    const unsigned NumChannels_pp = 4;
    const unsigned NumPotentials = 2;
    const int Potentials[NumPotentials] = {NN_AV18, pp_ReidSC};
    //const int Potentials[NumPotentials] = {NN_AV18, pp_ReidVale};
    const int PotFlags[NumPotentials] = {v18_Coupled3P2, 0};
    const double DataPoints1S0[NumDataPts] = {32.68,54.74,55.09,48.51,38.78,25.01,15,6.99,0.23,-5.64,-10.86};//1S0
    const double DataPoints3P0[NumDataPts] = {0.134,1.582,3.729,8.575,11.47,9.45,4.74,-0.37,-5.43,-10.39,-15.30};//3P0
    const double DataPoints3P1[NumDataPts] = {-0.081,-0.902,-2.060,-4.932,-8.317,-13.258,-17.434,-21.25,-24.77,-27.99,-30.89};//3P1
    const double DataPoints3P2[NumDataPts] = {0.014,0.214,0.651,2.491,5.855,11.013,13.982,15.63,16.59,17.17,17.54};//3P2
    const double DataPtsT[NumDataPts] = {1,5,10,25,50,100,150,200,250,300,350};

    const double** DataPoints = new const double* [NumChannels_pp];
    DataPoints[0] = DataPoints1S0;
    DataPoints[1] = DataPoints3P0;
    DataPoints[2] = DataPoints3P1;
    DataPoints[3] = DataPoints3P2;

    TString ChannelShortName[NumChannels_pp];
    ChannelShortName[0] = "1S0";
    ChannelShortName[1] = "3P0";
    ChannelShortName[2] = "3P1";
    ChannelShortName[3] = "3P2";

    TString PotShortName[NumPotentials];
    PotShortName[0] = "AV18";
    PotShortName[1] = "CastratedRSC";

    CATSparameters cSourcePars(CATSparameters::tSource,1,true);
    cSourcePars.SetParameters(&ppSourcePars[3]);

    CATS MyCat;
    MyCat.SetUseAnalyticSource(true);
    MyCat.SetAnaSource(GaussSource, cSourcePars);

    MyCat.SetExcludeFailedBins(false);
    MyCat.SetMomBins(NumPlotPts,Plot_kMin,Plot_kMax);

    MyCat.SetNumChannels(4);
    MyCat.SetNumPW(0,2);
    MyCat.SetNumPW(1,2);
    MyCat.SetNumPW(2,2);
    MyCat.SetNumPW(3,2);
    MyCat.SetSpin(0,0);
    MyCat.SetSpin(1,1);
    MyCat.SetSpin(2,1);
    MyCat.SetSpin(3,1);
    MyCat.SetChannelWeight(0, 3./12.);
    MyCat.SetChannelWeight(1, 1./12.);
    MyCat.SetChannelWeight(2, 3./12.);
    MyCat.SetChannelWeight(3, 5./12.);
    MyCat.SetShortRangePotential(0,0,fDlmPot,*ppPars1S0);
    MyCat.SetShortRangePotential(1,1,fDlmPot,*ppPars3P0);
    MyCat.SetShortRangePotential(2,1,fDlmPot,*ppPars3P1);
    MyCat.SetShortRangePotential(3,1,fDlmPot,*ppPars3P2);

    MyCat.SetQ1Q2(1);
    MyCat.SetPdgId(2212, 2212);

    MyCat.SetRedMass( 0.5*Mass_p );

    MyCat.SetEpsilonProp(1e-7);
    MyCat.SetEpsilonConv(1e-7);

    double Tlab;
    double pCM;
    DLM_DtColor DlmCol;
    TGraph** gCatsPs = new TGraph* [NumPotentials];
    TGraph* gData_k = new TGraph [NumChannels_pp];

    TPaveText** myPT = new TPaveText* [NumChannels_pp];

    DLM_SubPads DrawBoard_k(1920,1080);
    DrawBoard_k.AddSubPadTL(0,0.5,0.08,0.5);//lrtb
    DrawBoard_k.AddSubPadTL(0.5,1,0.08,0.5);//lrtb
    DrawBoard_k.AddSubPadTL(0,0.5,0.5,1);//lrtb
    DrawBoard_k.AddSubPadTL(0.5,1,0.5,1);//lrtb
    DrawBoard_k.AddSubPadTL(0,1,0,0.07);//lrtb

    DrawBoard_k.SetMargin(0, 0.08, 0.01, 0, 0.02);//lrbt
    DrawBoard_k.SetMargin(1, 0.08, 0.01, 0, 0.02);//lrbt
    DrawBoard_k.SetMargin(2, 0.08, 0.01, 0.11, 0);//lrbt
    DrawBoard_k.SetMargin(3, 0.08, 0.01, 0.11, 0);//lrbt
    DrawBoard_k.SetMargin(4, 0,0,0,0);//lrbt

    TH1F** hAxis = new TH1F* [NumChannels_pp];
    for(unsigned uCh=0; uCh<NumChannels_pp; uCh++){
        hAxis[uCh] = new TH1F (TString::Format("hAxis%u",uCh), TString::Format("hAxis%u",uCh), NumDataPts, Plot_kMin, Plot_kMax);
        hAxis[uCh]->SetStats(false);
        hAxis[uCh]->SetTitle("");
        hAxis[uCh]->GetXaxis()->SetLabelSize(0.06*1.1);
        hAxis[uCh]->GetXaxis()->SetTitle("k (MeV)");
        hAxis[uCh]->GetXaxis()->CenterTitle();
        hAxis[uCh]->GetXaxis()->SetTitleOffset(1.2);
        hAxis[uCh]->GetXaxis()->SetLabelOffset(0.02);
        hAxis[uCh]->GetXaxis()->SetTitleSize(0.07*1.1);
        hAxis[uCh]->GetXaxis()->SetNdivisions(506);

        hAxis[uCh]->GetYaxis()->SetLabelSize(0.06*1.1);
        hAxis[uCh]->GetYaxis()->SetTitle("#delta(k) (deg)");
        hAxis[uCh]->GetYaxis()->CenterTitle();
        hAxis[uCh]->GetYaxis()->SetTitleOffset(uCh<2?0.7:0.83);
        hAxis[uCh]->GetYaxis()->SetTitleSize(0.07*1.1);
        hAxis[uCh]->GetYaxis()->SetNdivisions(406);
        switch(uCh){
        case 0: hAxis[uCh]->GetYaxis()->SetRangeUser(-19,70);
                hAxis[uCh]->GetYaxis()->SetLimits(-19,70);
                break;
        case 1: hAxis[uCh]->GetYaxis()->SetRangeUser(-25,20);
                hAxis[uCh]->GetYaxis()->SetLimits(-25,20);
                break;
        case 2: hAxis[uCh]->GetYaxis()->SetRangeUser(-40,11);
                hAxis[uCh]->GetYaxis()->SetLimits(-40,11);
                break;
        case 3: hAxis[uCh]->GetYaxis()->SetRangeUser(-15,25);//-10, 40 original plots
                hAxis[uCh]->GetYaxis()->SetLimits(-15,25);//-10, 40 original plots
                break;
        default: break;
        }
        for(unsigned uData=0; uData<NumDataPts; uData++){
            Tlab = DataPtsT[uData];
            pCM = tLab_pCm(Tlab, Mass_p, Mass_p);
            gData_k[uCh].SetPoint(uData, pCM, DataPoints[uCh][uData]);
        }
        gData_k[uCh].SetLineColor(kRed);
        gData_k[uCh].SetLineWidth(6);
        gData_k[uCh].SetMarkerStyle(34);
        gData_k[uCh].SetMarkerSize(3.6);
        gData_k[uCh].SetMarkerColor(kRed);
        gData_k[uCh].SetFillColor(kWhite);

        switch(uCh){
            case 0:
                myPT[uCh] = new TPaveText(0.45,0.85,0.6,0.98, "blNDC");//lbrt
                myPT[uCh]->SetTextSize(0.08*1.1);
                myPT[uCh]->AddText("a) ^{1}S_{0}");
                break;
            case 1:
                myPT[uCh] = new TPaveText(0.25,0.85,0.4,0.98, "blNDC");//lbrt
                myPT[uCh]->SetTextSize(0.08*1.1);
                myPT[uCh]->AddText("b) ^{3}P_{0}");
                break;
            case 2:
                myPT[uCh] = new TPaveText(0.45,0.85,0.6,0.98, "blNDC");//lbrt
                myPT[uCh]->SetTextSize(0.08*0.84*1.1);
                myPT[uCh]->AddText("c) ^{3}P_{1}");
                break;
            case 3:
                myPT[uCh] = new TPaveText(0.25,0.85,0.4,0.98, "blNDC");//lbrt
                myPT[uCh]->SetTextSize(0.08*0.84*1.1);
                myPT[uCh]->AddText("d) ^{3}P_{2}");
                break;
            default: break;
        }
        myPT[uCh]->SetName(TString::Format("myPT_%u",uCh));
        myPT[uCh]->SetBorderSize(1);
        myPT[uCh]->SetFillColor(kWhite);
        myPT[uCh]->SetTextFont(22);
        //myPT[uCh]->Draw("same");

    }

    TF1* ZeroLine = new TF1("ZeroLine", "0", Plot_kMin, Plot_kMax);
    ZeroLine->SetLineColor(kGray);
    ZeroLine->SetLineWidth(4);

    for(unsigned uPot=0; uPot<NumPotentials; uPot++){
        ppPars1S0->SetParameter(0,Potentials[uPot]);
        ppPars1S0->SetParameter(1,PotFlags[uPot]);
        ppPars3P0->SetParameter(0,Potentials[uPot]);
        ppPars3P0->SetParameter(1,PotFlags[uPot]);
        ppPars3P1->SetParameter(0,Potentials[uPot]);
        ppPars3P1->SetParameter(1,PotFlags[uPot]);
        ppPars3P2->SetParameter(0,Potentials[uPot]);
        ppPars3P2->SetParameter(1,PotFlags[uPot]);

        MyCat.KillTheCat(CATS::kPotentialChanged);

        gCatsPs[uPot] = new TGraph [NumChannels_pp];

        for(unsigned uCh=0; uCh<NumChannels_pp; uCh++){
            gCatsPs[uPot][uCh].Set(NumPlotPts);

            for(unsigned uPlot=0; uPlot<NumPlotPts; uPlot++){
                //Tlab = DataPtsT[uData];
                //pCM = tLab_pCm(Tlab, Mass_p, Mass_p);
                pCM = CatHist.GetBinCenter(uPlot);
                //gCatsPs[uPot][uCh].SetPoint(uPlot, pCM, MyCat.EvalPhaseShift(pCM,uCh,bool(uCh))*180./3.14159);
                gCatsPs[uPot][uCh].SetPoint(uPlot, pCM, MyCat.GetPhaseShift(uPlot,uCh,bool(uCh))*180./3.14159);
            }

            gCatsPs[uPot][uCh].SetLineColor(DlmCol.GetColor(uPot+1));
            gCatsPs[uPot][uCh].SetLineWidth(9.-4.*double(uPot)/double(NumPotentials-1));
            //gCatsPs[uPot][uCh].SetMarkerStyle(34);
            gCatsPs[uPot][uCh].SetMarkerSize(4.2-1.*double(uPot)/double(NumPotentials-1));
            gCatsPs[uPot][uCh].SetMarkerColor(DlmCol.GetColor(uPot+1));
            gCatsPs[uPot][uCh].SetFillColor(kWhite);

            DrawBoard_k.cd(uCh);
            DrawBoard_k.SetLabelSize(uCh, hAxis[uCh]->GetXaxis(), 16);
            DrawBoard_k.SetLabelSize(uCh, hAxis[uCh]->GetYaxis(), 16);
            DrawBoard_k.SetTitleSize(uCh, hAxis[uCh]->GetXaxis(), 16);
            DrawBoard_k.SetTitleSize(uCh, hAxis[uCh]->GetYaxis(), 16);

            if(uPot==0){
                hAxis[uCh]->Draw("AXIS");
                ZeroLine->Draw("same");
                gData_k[uCh].Draw("same,P");
                myPT[uCh]->Draw("same");
            }
            gCatsPs[uPot][uCh].Draw("same,L");
        }

    }

    char* fBuffer = new char [512];

    FILE * CatsFile;
    CatsFile = fopen (TString::Format("%sPhaseShiftPlot_pp_Cats.txt",OutputDir.Data()),"w");
    fprintf (CatsFile, "%18s","k (MeV)");
    for(unsigned uCh=0; uCh<NumChannels_pp; uCh++){
        for(unsigned uPot=0; uPot<NumPotentials; uPot++){
            sprintf(fBuffer,"%s_%s",ChannelShortName[uCh].Data(),PotShortName[uPot].Data());
            fprintf (CatsFile,"%18s",fBuffer);
        }
    }

    double xVal,yVal;
    for(unsigned uMomBin=0; uMomBin<NumPlotPts; uMomBin++){
        fprintf (CatsFile, "\n%18.3f",MyCat.GetMomentum(uMomBin));
        for(unsigned uCh=0; uCh<NumChannels_pp; uCh++){
            for(unsigned uPot=0; uPot<NumPotentials; uPot++){
                gCatsPs[uPot][uCh].GetPoint(uMomBin,xVal,yVal);
                if(xVal!=MyCat.GetMomentum(uMomBin)) {printf("WTF MAN!\n");}
                fprintf (CatsFile, "%18.3f",yVal);
            }
        }
    }

    FILE * DataFile;
    DataFile = fopen (TString::Format("%sPhaseShiftPlot_pp_Data.txt",OutputDir.Data()),"w");
    fprintf (DataFile, "%18s","k (MeV)");
    for(unsigned uCh=0; uCh<NumChannels_pp; uCh++){
        sprintf(fBuffer,"%s_%s",ChannelShortName[uCh].Data(),"Data");
        fprintf (DataFile,"%18s",fBuffer);
    }


    for(unsigned uData=0; uData<NumDataPts; uData++){
        gData_k[0].GetPoint(uData,xVal,yVal);
        fprintf (DataFile, "\n%18.3f",xVal);
        for(unsigned uCh=0; uCh<NumChannels_pp; uCh++){
            gData_k[uCh].GetPoint(uData,xVal,yVal);
            fprintf (DataFile, "%18.3f",yVal);
        }
    }


    TLegend* myLegend = new TLegend(0.08,0,0.99,0.96);//lbrt
    myLegend->SetName("myLegend");
    myLegend->SetTextSize(0.6);
    myLegend->SetNColumns(3);
    myLegend->AddEntry(gData_k,    "Data          ");
    myLegend->AddEntry(gCatsPs[0], "Argonne V18   ");
    myLegend->AddEntry(gCatsPs[1], "Reid SC as in CRAB");

    DrawBoard_k.cd(4);

    myLegend->Draw();

    Width_t width_ZeroLine=ZeroLine->GetLineWidth();
    Width_t* width_gData_k = new Width_t [NumChannels_pp];
    for(unsigned uCh=0; uCh<NumChannels_pp; uCh++){
        width_gData_k[uCh] = gData_k[uCh].GetLineWidth();
    }
    Width_t** width_gCatsPs = new Width_t* [NumPotentials];
    for(unsigned uPot=0; uPot<NumPotentials; uPot++){
        width_gCatsPs[uPot] = new Width_t[NumChannels_pp];
        for(unsigned uCh=0; uCh<NumChannels_pp; uCh++){
            width_gCatsPs[uPot][uCh] = gCatsPs[uPot][uCh].GetLineWidth();
        }
    }

    DrawBoard_k.GetCanvas()->SaveAs(TString::Format("%sPhaseShiftPlot_pp.png",OutputDir.Data()));
    ZeroLine->SetLineWidth(width_ZeroLine/2.5);
    for(unsigned uCh=0; uCh<NumChannels_pp; uCh++){
        gData_k[uCh].SetLineWidth(width_gData_k[uCh]/2.5);
    }
    for(unsigned uPot=0; uPot<NumPotentials; uPot++){
        for(unsigned uCh=0; uCh<NumChannels_pp; uCh++){
            gCatsPs[uPot][uCh].SetLineWidth(width_gCatsPs[uPot][uCh]/2.5);
        }
    }
    DrawBoard_k.GetCanvas()->SaveAs(TString::Format("%sPhaseShiftPlot_pp.pdf",OutputDir.Data()));
    ZeroLine->SetLineWidth(width_ZeroLine);
    for(unsigned uCh=0; uCh<NumChannels_pp; uCh++){
        gData_k[uCh].SetLineWidth(width_gData_k[uCh]);
    }
    for(unsigned uPot=0; uPot<NumPotentials; uPot++){
        width_gCatsPs[uPot] = new Width_t[NumChannels_pp];
        for(unsigned uCh=0; uCh<NumChannels_pp; uCh++){
            gCatsPs[uPot][uCh].SetLineWidth(width_gCatsPs[uPot][uCh]);
        }
    }

    MyCat.SetEpsilonProp(5e-6);
    MyCat.SetEpsilonConv(5e-6);

    ppPars1S0->SetParameter(0,NN_AV18);
    ppPars1S0->SetParameter(1,v18_Coupled3P2);
    ppPars3P0->SetParameter(0,NN_AV18);
    ppPars3P0->SetParameter(1,v18_Coupled3P2);
    ppPars3P1->SetParameter(0,NN_AV18);
    ppPars3P1->SetParameter(1,v18_Coupled3P2);
    ppPars3P2->SetParameter(0,NN_AV18);
    ppPars3P2->SetParameter(1,v18_Coupled3P2);

    for(unsigned uPot=0; uPot<NumPotentials; uPot++){
        delete [] gCatsPs[uPot];
        delete [] width_gCatsPs[uPot];
    }
    delete [] gCatsPs;
    delete [] gData_k;
    delete [] width_gCatsPs;
    delete [] width_gData_k;
    for(unsigned uCh=0; uCh<NumChannels_pp; uCh++){
        delete hAxis[uCh];
        delete myPT[uCh];
    }
    delete [] myPT;
    delete [] hAxis;
    delete [] DataPoints;
    delete ZeroLine;
    delete myLegend;
    delete [] fBuffer;
}

void CATSpaper::PlotCS_pLambda(){
    const double kFrom = 5;
    const double kTo = 145;
    const unsigned kNumStep = 14;

    const unsigned NumDataPts = 12;
    const double DataPoints[NumDataPts] = {209,177,153,111,87,46,180,130,118,101,83,57};
    const double DataPointsErr[NumDataPts] = {58,38,27,18,13,11,22,17,16,12,9,9};
    const double DataPtsPlab[NumDataPts] = {135,165,194,226,252,293,145,185,210,230,250,290};

    CATS Afterburner;
    double Pars[4] = {0,0,0,1.3};
    CATSparameters cSourcePars(CATSparameters::tSource,1,true);
    cSourcePars.SetParameters(&Pars[3]);
    Afterburner.SetUseAnalyticSource(true);
    Afterburner.SetAnaSource(GaussSource, cSourcePars);
    //Afterburner.SetMaxRad(100);
    //Afterburner.SetMaxRho(50);
    Afterburner.SetEpsilonProp(1e-7);
    Afterburner.SetEpsilonConv(1e-7);
    Afterburner.SetExcludeFailedBins(false);
    Afterburner.SetMomBins(kNumStep,kFrom,kTo);

    Afterburner.SetNumChannels(2);
    Afterburner.SetNumPW(0,2);
    Afterburner.SetNumPW(1,2);
    Afterburner.SetSpin(0,0);
    Afterburner.SetSpin(1,1);
    Afterburner.SetChannelWeight(0, 0.25);
    Afterburner.SetChannelWeight(1, 0.75);

    Afterburner.SetQ1Q2(0);
    Afterburner.SetPdgId(2212, 3122);
    double Mass1=1115.683; double Mass2=938.272;
    Afterburner.SetRedMass( (Mass1*Mass2)/(Mass1+Mass2) );;

    double pCM;
    double pLAB;
    double CStot;

    //#,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
    double pLamPotPars1S0[10]={0,0,pL_UsmaniOli,0,0,0,0,0,0,0};
    double pLamPotPars3S1[10]={0,0,pL_UsmaniOli,0,0,0,0,1,0,1};
    CATSparameters c_pLamPotPars1S0(CATSparameters::tPotential,8,true);
    c_pLamPotPars1S0.SetParameters(&pLamPotPars1S0[2]);
    CATSparameters c_pLamPotPars3S1(CATSparameters::tPotential,8,true);
    c_pLamPotPars3S1.SetParameters(&pLamPotPars3S1[2]);
    //Afterburner.SetShortRangePotential(UsmaniPotential, PotPars);
    Afterburner.SetShortRangePotential(0,0,fDlmPot,c_pLamPotPars1S0);
    Afterburner.SetShortRangePotential(1,0,fDlmPot,c_pLamPotPars3S1);

    Afterburner.KillTheCat();
    TH1F* hCs_Pot1_pCm = new TH1F("hCs_Pot1_pCm", "hCs_Pot1_pCm", kNumStep, kFrom, kTo);
    TGraphErrors gCs_Pot1_pCm;
    gCs_Pot1_pCm.Set(kNumStep);
    //TH1F* hDiffPot1Pot2_k = new TH1F("hDiffPot1Pot2_k", "hDiffPot1Pot2_k", kNumStep, kFrom, kTo);
    for(unsigned iBin=0; iBin<kNumStep; iBin++){
        pCM = Afterburner.GetMomentum(iBin);
        CStot = 0;
        CStot += Afterburner.GetChannelWeight(0)*(pow(sin(Afterburner.GetPhaseShift(iBin,0,0)),2));//singlet
        CStot += Afterburner.GetChannelWeight(1)*(pow(sin(Afterburner.GetPhaseShift(iBin,1,0)),2));//triplet
        CStot *= (4.*3.14159/pCM/pCM);//units are in 1/MeV^2
        CStot *= pow(Afterburner.NuFm(),2);//convert to fm^2
        CStot *= 10.; //fm to mbarn

        hCs_Pot1_pCm->SetBinContent(iBin+1, CStot);
        hCs_Pot1_pCm->SetBinError(iBin+1, CStot*0.0001);

        gCs_Pot1_pCm.SetPoint(iBin, pCM, CStot);
        gCs_Pot1_pCm.SetPointError(iBin, 0, 0);
    }

    complex<double>**** WaveFunctionU; double*** PhaseShifts; double* RadBins; unsigned NumRadBins=0;
    InitHaidenbauerNLO("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/Haidenbauer/pLambdaNLO/",Afterburner,&WaveFunctionU,&PhaseShifts,&RadBins,NumRadBins);

//[0][uMomBin][usCh][usPw]
    TGraphErrors gCs_NLO_pCm;
    gCs_NLO_pCm.Set(kNumStep);
    for(unsigned iBin=0; iBin<kNumStep; iBin++){
        pCM = Afterburner.GetMomentum(iBin);
        CStot = 0;
        CStot += Afterburner.GetChannelWeight(0)*(pow(sin(PhaseShifts[iBin][0][0]),2));//singlet
        CStot += Afterburner.GetChannelWeight(1)*(pow(sin(PhaseShifts[iBin][1][0]),2));//triplet
        CStot *= (4.*3.14159/pCM/pCM);//units are in 1/MeV^2
        CStot *= pow(Afterburner.NuFm(),2);//convert to fm^2
        CStot *= 10.; //fm to mbarn

        hCs_Pot1_pCm->SetBinContent(iBin+1, CStot);
        hCs_Pot1_pCm->SetBinError(iBin+1, CStot*0.0001);

        gCs_NLO_pCm.SetPoint(iBin, pCM, CStot);
        gCs_NLO_pCm.SetPointError(iBin, 0, 0);
    }

    TGraphErrors gData_pCm;
    gData_pCm.Set(NumDataPts);
    for(unsigned uData=0; uData<NumDataPts; uData++){
        pLAB = DataPtsPlab[uData];
        pCM = pLab_pCm(pLAB, Mass1, Mass2);
        gData_pCm.SetPoint(uData, pCM, DataPoints[uData]);
        gData_pCm.SetPointError(uData, 0, DataPointsErr[uData]);
    }

    TH1F* hRatioPot1_pCm = new TH1F("hRatioPot1_pCm", "hRatioPot1_pCm", kNumStep, kFrom, kTo);
    TGraphErrors gRatioPot1_pCm;
    gRatioPot1_pCm.Set(NumDataPts);

    TGraphErrors gRatioNLO_pCm;
    gRatioNLO_pCm.Set(NumDataPts);

    //double Ratio;
    double Numerator;
    double Denumerator;
    double DenumError;
    double DummyX;

    for(unsigned uData=0; uData<NumDataPts; uData++){
        pLAB = DataPtsPlab[uData];
        pCM = pLab_pCm(pLAB, Mass1, Mass2);
        Numerator = gCs_Pot1_pCm.Eval(pCM);
        gData_pCm.GetPoint(uData, DummyX, Denumerator);
        DenumError = gData_pCm.GetErrorY(uData);
        gRatioPot1_pCm.SetPoint(uData, pCM, Numerator/Denumerator);
        gRatioPot1_pCm.SetPoint(uData, pCM, gCs_Pot1_pCm.Eval(pCM)/gData_pCm.Eval(pCM));
        gRatioPot1_pCm.SetPointError(uData, 0, pow(Denumerator,-2.)*sqrt(Numerator*Numerator*DenumError*DenumError));
    }
    for(unsigned uData=0; uData<NumDataPts; uData++){
        pLAB = DataPtsPlab[uData];
        pCM = pLab_pCm(pLAB, Mass1, Mass2);
        Numerator = gCs_NLO_pCm.Eval(pCM);
        gData_pCm.GetPoint(uData, DummyX, Denumerator);
        DenumError = gData_pCm.GetErrorY(uData);
        gRatioNLO_pCm.SetPoint(uData, pCM, Numerator/Denumerator);
        gRatioNLO_pCm.SetPoint(uData, pCM, gCs_NLO_pCm.Eval(pCM)/gData_pCm.Eval(pCM));
        gRatioNLO_pCm.SetPointError(uData, 0, pow(Denumerator,-2.)*sqrt(Numerator*Numerator*DenumError*DenumError));
    }

    DLM_DtColor DlmCol;

    char* fBuffer = new char [512];
    FILE * CatsFile;
    CatsFile = fopen (TString::Format("%sPlotCS_pLambda_Cats.txt",OutputDir.Data()),"w");
    fprintf (CatsFile, "%18s","k (MeV)");
    fprintf (CatsFile,"%18s","Usmani");
    fprintf (CatsFile,"%18s","NLO");

    double xVal,yVal;
    for(unsigned iBin=0; iBin<kNumStep; iBin++){
        fprintf (CatsFile, "\n%18.3f",Afterburner.GetMomentum(iBin));
        gCs_Pot1_pCm.GetPoint(iBin,xVal,yVal);
        if(xVal!=Afterburner.GetMomentum(iBin)) {printf("WTF MAN!\n");}
        fprintf (CatsFile, "%18.3f",yVal);
        gCs_NLO_pCm.GetPoint(iBin,xVal,yVal);
        if(xVal!=Afterburner.GetMomentum(iBin)) {printf("WTF MAN!\n");}
        fprintf (CatsFile, "%18.3f",yVal);
    }
//STILL LEFT TO DO IS THE DATA AND RATIO

    FILE * DataFile;
    DataFile = fopen (TString::Format("%sPlotCS_pLambda_Data.txt",OutputDir.Data()),"w");
    fprintf (DataFile, "%18s","k (MeV)");
    fprintf (DataFile,"%18s","Data");
    fprintf (DataFile,"%18s","DataErr");
    fprintf (DataFile,"%18s","RatioUsm");
    fprintf (DataFile,"%18s","RatioUsmErr");
    fprintf (DataFile,"%18s","RatioNlo");
    fprintf (DataFile,"%18s","RatioNloErr");

    for(unsigned uData=0; uData<NumDataPts; uData++){
        gData_pCm.GetPoint(uData,xVal,yVal);
        fprintf (DataFile, "\n%18.3f",xVal);
        fprintf (DataFile, "%18.3f",yVal);
        yVal = gData_pCm.GetErrorY(uData);
        fprintf (DataFile, "%18.3f",yVal);

        gRatioPot1_pCm.GetPoint(uData,xVal,yVal);
        fprintf (DataFile, "%18.3f",yVal);
        yVal = gRatioPot1_pCm.GetErrorY(uData);
        fprintf (DataFile, "%18.3f",yVal);

        gRatioNLO_pCm.GetPoint(uData,xVal,yVal);
        fprintf (DataFile, "%18.3f",yVal);
        yVal = gRatioNLO_pCm.GetErrorY(uData);
        fprintf (DataFile, "%18.3f",yVal);
    }

    TF1* UnitLine_pCm = new TF1("UnitLine_pCm", "1", kFrom,kTo);
    UnitLine_pCm->SetLineColor(kGray);
    UnitLine_pCm->SetLineWidth(5);

    hCs_Pot1_pCm->SetStats(false);
    hCs_Pot1_pCm->SetTitle("");
    hCs_Pot1_pCm->GetXaxis()->SetLabelSize(0.06*1.25);
    hCs_Pot1_pCm->GetXaxis()->SetTitle("k (MeV)");
    hCs_Pot1_pCm->GetXaxis()->CenterTitle();
    hCs_Pot1_pCm->GetXaxis()->SetTitleOffset(1.12);
    hCs_Pot1_pCm->GetXaxis()->SetLabelOffset(0.02);
    hCs_Pot1_pCm->GetXaxis()->SetTitleSize(0.07*1.25);
    hCs_Pot1_pCm->GetYaxis()->SetLabelSize(0.06*1.25);
    hCs_Pot1_pCm->GetYaxis()->SetTitle("#sigma_{p#Lambda}(k) (mb)");
    hCs_Pot1_pCm->GetYaxis()->CenterTitle();
    hCs_Pot1_pCm->GetYaxis()->SetTitleOffset(0.6);
    hCs_Pot1_pCm->GetYaxis()->SetTitleSize(0.07*1.25);

    hCs_Pot1_pCm->GetYaxis()->SetRangeUser(0, 600);
    hCs_Pot1_pCm->GetYaxis()->SetLimits(0, 600);

    hCs_Pot1_pCm->SetLineColor(DlmCol.GetColor(1));
    hCs_Pot1_pCm->SetLineWidth(6);
    hCs_Pot1_pCm->SetMarkerStyle(20);
    hCs_Pot1_pCm->SetMarkerSize(2.6);
    hCs_Pot1_pCm->SetMarkerColor(DlmCol.GetColor(1));

    gCs_Pot1_pCm.SetLineColor(DlmCol.GetColor(1));
    gCs_Pot1_pCm.SetLineWidth(9);
    gCs_Pot1_pCm.SetMarkerStyle(20);
    gCs_Pot1_pCm.SetMarkerSize(2.6);
    gCs_Pot1_pCm.SetMarkerColor(DlmCol.GetColor(1));
    gCs_Pot1_pCm.SetFillColor(kWhite);

    gCs_NLO_pCm.SetLineColor(DlmCol.GetColor(2));
    gCs_NLO_pCm.SetLineWidth(6);
    gCs_NLO_pCm.SetMarkerStyle(20);
    gCs_NLO_pCm.SetMarkerSize(2.6);
    gCs_NLO_pCm.SetMarkerColor(DlmCol.GetColor(2));
    gCs_NLO_pCm.SetFillColor(kWhite);

    hRatioPot1_pCm->SetStats(false);
    hRatioPot1_pCm->SetTitle("");
    hRatioPot1_pCm->GetXaxis()->SetLabelSize(0.06*1.25);
    hRatioPot1_pCm->GetXaxis()->SetTitle("k (MeV)");
    hRatioPot1_pCm->GetXaxis()->CenterTitle();
    hRatioPot1_pCm->GetXaxis()->SetTitleOffset(1.12);
    hRatioPot1_pCm->GetXaxis()->SetLabelOffset(0.02);
    hRatioPot1_pCm->GetXaxis()->SetTitleSize(0.07*1.25);
    hRatioPot1_pCm->GetYaxis()->SetLabelSize(0.06*1.25);
    hRatioPot1_pCm->GetYaxis()->SetTitle("Ratio");
    hRatioPot1_pCm->GetYaxis()->CenterTitle();
    hRatioPot1_pCm->GetYaxis()->SetTitleOffset(0.4);
    hRatioPot1_pCm->GetYaxis()->SetTitleSize(0.07*1.25);

    hRatioPot1_pCm->GetYaxis()->SetRangeUser(0.63, 1.6);
    hRatioPot1_pCm->GetYaxis()->SetLimits(0.63, 1.6);
    hRatioPot1_pCm->GetYaxis()->SetNdivisions(404);

    hRatioPot1_pCm->SetLineColor(kBlack);
    hRatioPot1_pCm->SetLineWidth(9);
    hRatioPot1_pCm->SetMarkerStyle(21);
    hRatioPot1_pCm->SetMarkerSize(2.8);
    hRatioPot1_pCm->SetMarkerColor(kBlack);

    gRatioPot1_pCm.SetLineColor(DlmCol.GetColor(1));
    gRatioPot1_pCm.SetLineWidth(9);
    gRatioPot1_pCm.SetMarkerStyle(21);
    gRatioPot1_pCm.SetMarkerSize(2.8);
    gRatioPot1_pCm.SetMarkerColor(DlmCol.GetColor(1));
    gRatioPot1_pCm.SetFillColor(kWhite);

    gRatioNLO_pCm.SetLineColor(DlmCol.GetColor(2));
    gRatioNLO_pCm.SetLineWidth(6);
    gRatioNLO_pCm.SetMarkerStyle(21);
    gRatioNLO_pCm.SetMarkerSize(2.1);
    gRatioNLO_pCm.SetMarkerColor(DlmCol.GetColor(2));
    gRatioNLO_pCm.SetFillColor(kWhite);

    gData_pCm.SetLineColor(kRed);
    gData_pCm.SetLineWidth(9);
    gData_pCm.SetMarkerStyle(20);
    gData_pCm.SetMarkerSize(2.8);
    gData_pCm.SetMarkerColor(kRed);
    gData_pCm.SetFillColor(kWhite);


    DLM_SubPads DrawBoard_pCm(1920,1080);
    DrawBoard_pCm.AddSubPadTL(0,1,0,0.6);//lrtb
    DrawBoard_pCm.AddSubPadTL(0,1,0.6,1);//lrtb

    DrawBoard_pCm.SetMargin(0, 0.14, 0.04, 0, 0.04);//lrbt
    DrawBoard_pCm.SetMargin(1, 0.14, 0.04, 0.14, 0);//lrbt

    DrawBoard_pCm.GetCanvas()->SetName(TString::Format("DrawBoard_pCm"));
    //DrawBoard_pCm.SetLogy(0, true);
//! originally I was plotting gCs_Pot1_pCm/gData_pCm (symbols), but eventually changed to a line
    DrawBoard_pCm.cd(0);
    DrawBoard_pCm.SetGrid(0, true);
    hCs_Pot1_pCm->Draw("AXIS");
    gCs_Pot1_pCm.Draw("same,C");
    gCs_NLO_pCm.Draw("same,C");
    gData_pCm.Draw("same,P");
    //UnitLine->Draw("same");
    DrawBoard_pCm.SetLabelSize(0, hCs_Pot1_pCm->GetXaxis(), 18*1.2);
    DrawBoard_pCm.SetLabelSize(0, hCs_Pot1_pCm->GetYaxis(), 18*1.2);
    DrawBoard_pCm.SetTitleSize(0, hCs_Pot1_pCm->GetXaxis(), 18*1.2);
    DrawBoard_pCm.SetTitleSize(0, hCs_Pot1_pCm->GetYaxis(), 18*1.2);

    DrawBoard_pCm.cd(1);
    DrawBoard_pCm.SetGrid(1, true);
    hRatioPot1_pCm->Draw("AXIS");
    //gRatioPot1_pCm.Draw("AXIS");
    UnitLine_pCm->Draw("same");
    gRatioPot1_pCm.Draw("same,P");
    gRatioNLO_pCm.Draw("same,P");
    DrawBoard_pCm.SetLabelSize(1, hRatioPot1_pCm->GetXaxis(), 18*1.2);
    DrawBoard_pCm.SetLabelSize(1, hRatioPot1_pCm->GetYaxis(), 18*1.2);
    DrawBoard_pCm.SetTitleSize(1, hRatioPot1_pCm->GetXaxis(), 18*1.2);
    DrawBoard_pCm.SetTitleSize(1, hRatioPot1_pCm->GetYaxis(), 18*1.2);

    DrawBoard_pCm.GetCanvas()->cd(0);
    TLegend* myLegend_pCm = new TLegend(0.7,0.7,1,1);//lbrt
    myLegend_pCm->SetName("myLegend_pCm");
    myLegend_pCm->SetTextSize(0.05);
    myLegend_pCm->AddEntry(&gCs_Pot1_pCm, "CATS+Usmani");
    myLegend_pCm->AddEntry(&gCs_NLO_pCm, "NLO");
    myLegend_pCm->AddEntry(&gData_pCm, "Exp. data");
    myLegend_pCm->AddEntry(hRatioPot1_pCm, "Ratio to data");
    myLegend_pCm->Draw("same");

    DrawBoard_pCm.GetCanvas()->SaveAs(TString::Format("%sPlotCS_pLambda.png",OutputDir.Data()));
    hCs_Pot1_pCm->SetLineWidth(hCs_Pot1_pCm->GetLineWidth()/2.5);
    hRatioPot1_pCm->SetLineWidth(hRatioPot1_pCm->GetLineWidth()/2.5);
    gCs_Pot1_pCm.SetLineWidth(gCs_Pot1_pCm.GetLineWidth()/2.5);
    gCs_NLO_pCm.SetLineWidth(gCs_NLO_pCm.GetLineWidth()/2.5);
    gData_pCm.SetLineWidth(gData_pCm.GetLineWidth()/2.5);
    gRatioPot1_pCm.SetLineWidth(gRatioPot1_pCm.GetLineWidth()/2.5);
    gRatioNLO_pCm.SetLineWidth(gRatioNLO_pCm.GetLineWidth()/2.5);
    UnitLine_pCm->SetLineWidth(UnitLine_pCm->GetLineWidth()/2.5);
    DrawBoard_pCm.GetCanvas()->SaveAs(TString::Format("%sPlotCS_pLambda.pdf",OutputDir.Data()));
    hCs_Pot1_pCm->SetLineWidth(hCs_Pot1_pCm->GetLineWidth()*2.5);
    hRatioPot1_pCm->SetLineWidth(hRatioPot1_pCm->GetLineWidth()*2.5);
    gCs_Pot1_pCm.SetLineWidth(gCs_Pot1_pCm.GetLineWidth()*2.5);
    gCs_NLO_pCm.SetLineWidth(gCs_NLO_pCm.GetLineWidth()*2.5);
    gData_pCm.SetLineWidth(gData_pCm.GetLineWidth()*2.5);
    gRatioPot1_pCm.SetLineWidth(gRatioPot1_pCm.GetLineWidth()*2.5);
    gRatioNLO_pCm.SetLineWidth(gRatioNLO_pCm.GetLineWidth()*2.5);
    UnitLine_pCm->SetLineWidth(UnitLine_pCm->GetLineWidth()*2.5);

    CleanHaidenbauer(Afterburner,&WaveFunctionU,&PhaseShifts,&RadBins);

    delete [] fBuffer;
}

void CATSpaper::PlotCk_pp(){
    const double Plot_kMin = 0;
    const double Plot_kMax = 200;
    const double NumPlotPts = 50;
    //CATShisto<double> CatHist(NumPlotPts,Plot_kMin,Plot_kMax);

    const double RadLarge = 1.2;
    const double RadSmall = 0.85;

    double CkVal;
    double Momentum;
    DLM_DtColor DlmCol;

    CATSparameters cSourcePars(CATSparameters::tSource,1,true);
    cSourcePars.SetParameters(&ppSourcePars[3]);

    CATS MyCat;
    MyCat.SetUseAnalyticSource(true);
    MyCat.SetAnaSource(GaussSource, cSourcePars);

    MyCat.SetExcludeFailedBins(false);
    MyCat.SetMomBins(NumPlotPts,Plot_kMin,Plot_kMax);

    MyCat.SetNumChannels(4);
    MyCat.SetNumPW(0,2);
    MyCat.SetNumPW(1,2);
    MyCat.SetNumPW(2,2);
    MyCat.SetNumPW(3,2);
    MyCat.SetSpin(0,0);
    MyCat.SetSpin(1,1);
    MyCat.SetSpin(2,1);
    MyCat.SetSpin(3,1);
    MyCat.SetChannelWeight(0, 3./12.);
    MyCat.SetChannelWeight(1, 1./12.);
    MyCat.SetChannelWeight(2, 3./12.);
    MyCat.SetChannelWeight(3, 5./12.);
    MyCat.SetShortRangePotential(0,0,fDlmPot,*ppPars1S0);
    MyCat.SetShortRangePotential(1,1,fDlmPot,*ppPars3P0);
    MyCat.SetShortRangePotential(2,1,fDlmPot,*ppPars3P1);
    MyCat.SetShortRangePotential(3,1,fDlmPot,*ppPars3P2);
    MyCat.SetQ1Q2(1);
    MyCat.SetPdgId(2212, 2212);
    MyCat.SetRedMass( 0.5*Mass_p );

    ppSourcePars[3] = RadLarge;
    MyCat.KillTheCat();

    TGraph gAV18_G12;
    gAV18_G12.Set(NumPlotPts);
    gAV18_G12.SetName("gAV18_G12");
    gAV18_G12.SetMarkerColor(DlmCol.GetColor(1));
    gAV18_G12.SetMarkerStyle(20);
    gAV18_G12.SetMarkerSize(0);
    gAV18_G12.SetLineColor(DlmCol.GetColor(1));
    gAV18_G12.SetLineWidth(10);
    gAV18_G12.SetFillColor(kWhite);
    for(unsigned uBin=0; uBin<NumPlotPts; uBin++){
        gAV18_G12.SetPoint(uBin, MyCat.GetMomentum(uBin), MyCat.GetCorrFun(uBin));
    }

    MyCat.SetAnaSource(0, RadSmall);
    MyCat.KillTheCat();
    TGraph gAV18_G08;
    gAV18_G08.Set(NumPlotPts);
    gAV18_G08.SetName("gAV18_G08");
    gAV18_G08.SetMarkerColor(DlmCol.GetColor(1));
    gAV18_G08.SetMarkerStyle(20);
    gAV18_G08.SetMarkerSize(0);
    gAV18_G08.SetLineColor(DlmCol.GetColor(1));
    gAV18_G08.SetLineWidth(10);
    gAV18_G08.SetFillColor(kWhite);
    for(unsigned uBin=0; uBin<NumPlotPts; uBin++){
        gAV18_G08.SetPoint(uBin, MyCat.GetMomentum(uBin), MyCat.GetCorrFun(uBin));
    }
/*
    MyCat.SetUseAnalyticSource(false);
    MyCat.SetMaxPairsPerBin(8000);
    //Afterburner.SetMaxPairsToRead(8e6);
    //MyCat.SetMixingDepth(8);
    MyCat.SetThetaDependentSource(true);
    MyCat.SetInputFileName("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/Scratch9_OSCAR1997_AllLrzProtonFiles.f19");
    MyCat.KillTheCat(CATS::kSourceChanged);

    TGraphErrors gAV18_EPOS;
    gAV18_EPOS.Set(NumPlotPts);
    gAV18_EPOS.SetName("gAV18_EPOS");
    gAV18_EPOS.SetMarkerColor(kOrange+1);
    gAV18_EPOS.SetMarkerStyle(20);
    gAV18_EPOS.SetMarkerSize(0);
    gAV18_EPOS.SetLineColor(kOrange+1);
    gAV18_EPOS.SetLineWidth(10);
    gAV18_EPOS.SetFillColor(kWhite);
    for(unsigned uBin=0; uBin<NumPlotPts; uBin++){
        gAV18_EPOS.SetPoint(uBin, MyCat.GetMomentum(uBin), MyCat.GetCorrFun(uBin));
        gAV18_EPOS.SetPointError(uBin, 0, MyCat.GetCorrFunErr(uBin));
    }
*/
    ppPars1S0[2] = pp_ReidSC;
    ppPars1S0[3] = 0;
    ppPars3P0[2] = pp_ReidSC;
    ppPars3P0[3] = 0;
    ppPars3P1[2] = pp_ReidSC;
    ppPars3P1[3] = 0;
    ppPars3P2[2] = pp_ReidSC;
    ppPars3P2[3] = 0;
/*
    ppPars1S0[2] = NN_ReidV8;
    ppPars1S0[3] = v18_Coupled3P2;
    ppPars3P0[2] = NN_ReidV8;
    ppPars3P0[3] = v18_Coupled3P2;
    ppPars3P1[2] = NN_ReidV8;
    ppPars3P1[3] = v18_Coupled3P2;
    ppPars3P2[2] = NN_ReidV8;
    ppPars3P2[3] = v18_Coupled3P2;
*/
    //MyCat.SetChannelWeight(0, 1./4.);
    //MyCat.SetChannelWeight(1, 0);
    //MyCat.SetChannelWeight(2, 0);
    //MyCat.SetChannelWeight(3, 3./4.);

    MyCat.SetUseAnalyticSource(true);
    ppSourcePars[3] = RadLarge;
    MyCat.KillTheCat(CATS::kAllChanged);

    TGraph gRSC_G12;
    gRSC_G12.Set(NumPlotPts);
    gRSC_G12.SetName("gRSC_G12");
    gRSC_G12.SetMarkerColor(DlmCol.GetColor(2));
    gRSC_G12.SetMarkerStyle(20);
    gRSC_G12.SetMarkerSize(0);
    gRSC_G12.SetLineColor(DlmCol.GetColor(2));
    gRSC_G12.SetLineWidth(6);
    gRSC_G12.SetFillColor(kWhite);

    TGraph gRatio_G12;
    gRatio_G12.Set(NumPlotPts);
    gRatio_G12.SetName("gRatio_G12");
    gRatio_G12.SetMarkerColor(kBlack);
    gRatio_G12.SetMarkerStyle(20);
    gRatio_G12.SetMarkerSize(0);
    gRatio_G12.SetLineColor(kBlack);
    gRatio_G12.SetLineWidth(10);
    gRatio_G12.SetFillColor(kWhite);

    for(unsigned uBin=0; uBin<NumPlotPts; uBin++){
        gRSC_G12.SetPoint(uBin, MyCat.GetMomentum(uBin), MyCat.GetCorrFun(uBin));

        gAV18_G12.GetPoint(uBin, Momentum, CkVal);
        gRatio_G12.SetPoint(uBin, Momentum, MyCat.GetCorrFun(uBin)/CkVal);
    }

    MyCat.SetAnaSource(0, RadSmall);
    MyCat.KillTheCat();
    TGraph gRSC_G08;
    gRSC_G08.Set(NumPlotPts);
    gRSC_G08.SetName("gRSC_G08");
    gRSC_G08.SetMarkerColor(DlmCol.GetColor(2));
    gRSC_G08.SetMarkerStyle(20);
    gRSC_G08.SetMarkerSize(0);
    gRSC_G08.SetLineColor(DlmCol.GetColor(2));
    gRSC_G08.SetLineWidth(6);
    gRSC_G08.SetFillColor(kWhite);

    TGraph gRatio_G08;
    gRatio_G08.Set(NumPlotPts);
    gRatio_G08.SetName("gRatio_G08");
    gRatio_G08.SetMarkerColor(kBlack);
    gRatio_G08.SetMarkerStyle(20);
    gRatio_G08.SetMarkerSize(0);
    gRatio_G08.SetLineColor(kBlack);
    gRatio_G08.SetLineWidth(10);
    gRatio_G08.SetFillColor(kWhite);

    for(unsigned uBin=0; uBin<NumPlotPts; uBin++){
        gRSC_G08.SetPoint(uBin, MyCat.GetMomentum(uBin), MyCat.GetCorrFun(uBin));

        gAV18_G08.GetPoint(uBin, Momentum, CkVal);
        gRatio_G08.SetPoint(uBin, Momentum, MyCat.GetCorrFun(uBin)/CkVal);
    }

    const double Yoffset = 1;
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

    hAxis->GetYaxis()->SetRangeUser(0, 6.2);
    hAxis->GetYaxis()->SetLimits(0, 6.2);

    hAxis->GetXaxis()->SetNdivisions(504);

    TH1F* hAxisRatio = new TH1F("hAxisRatio", "hAxisRatio", NumPlotPts,Plot_kMin,Plot_kMax);
    hAxisRatio->SetStats(false);
    hAxisRatio->SetTitle("");
    hAxisRatio->GetXaxis()->SetLabelSize(0.065);
    hAxisRatio->GetXaxis()->SetTitle("k (MeV)");
    hAxisRatio->GetXaxis()->CenterTitle();
    hAxisRatio->GetXaxis()->SetTitleOffset(1.25);
    hAxisRatio->GetXaxis()->SetLabelOffset(0.02);
    hAxisRatio->GetXaxis()->SetTitleSize(0.075);
    hAxisRatio->GetYaxis()->SetLabelSize(0.065);
    hAxisRatio->GetYaxis()->SetTitle("Ratio");
    hAxisRatio->GetYaxis()->CenterTitle();
    hAxisRatio->GetYaxis()->SetTitleOffset(Yoffset/1.5);
    hAxisRatio->GetYaxis()->SetTitleSize(0.075);

    hAxisRatio->GetXaxis()->SetNdivisions(504);
    hAxisRatio->GetYaxis()->SetRangeUser(0.85, 1.09);
    hAxisRatio->GetYaxis()->SetLimits(0.85, 1.09);
    hAxisRatio->GetYaxis()->SetNdivisions(604);

    TPaveText* PT12 = new TPaveText(0.5,0.68,0.95,0.95, "blNDC");//lbrt
    PT12->SetName("PT12");
    PT12->SetBorderSize(1);
    PT12->SetTextSize(0.075);
    PT12->SetFillColor(kWhite);
    PT12->SetTextFont(22);
    PT12->AddText("Gauss source");
    PT12->AddText(TString::Format("R = %.2f fm", RadLarge));

    TPaveText* PT08 = new TPaveText(0.5,0.68,0.95,0.95, "blNDC");//lbrt
    PT08->SetName("PT08");
    PT08->SetBorderSize(1);
    PT08->SetTextSize(0.075);
    PT08->SetFillColor(kWhite);
    PT08->SetTextFont(22);
    PT08->AddText("Gauss source");
    PT08->AddText(TString::Format("R = %.2f fm", RadSmall));

    TF1* UnitLine = new TF1("UnitLine", "1",Plot_kMin,Plot_kMax);
    UnitLine->SetLineColor(kGray+1);
    UnitLine->SetLineWidth(6);

    DLM_SubPads DrawBoard_Ck_pp(1920,1080);
    DrawBoard_Ck_pp.AddSubPadTL(0,0.5,0.08,0.63);//lrtb
    DrawBoard_Ck_pp.AddSubPadTL(0,0.5,0.63,1);//lrtb
    DrawBoard_Ck_pp.AddSubPadTL(0.5,1,0.08,0.63);//lrtb
    DrawBoard_Ck_pp.AddSubPadTL(0.5,1,0.63,1);//lrtb
    DrawBoard_Ck_pp.AddSubPadTL(0,1,0,0.07);//lrtb

    DrawBoard_Ck_pp.SetMargin(0, 0.09, 0.04, 0, 0.04);//lrbt
    DrawBoard_Ck_pp.SetMargin(1, 0.09, 0.04, 0.14, 0);//lrbt
    DrawBoard_Ck_pp.SetMargin(2, 0.09, 0.04, 0, 0.04);//lrbt
    DrawBoard_Ck_pp.SetMargin(3, 0.09, 0.04, 0.14, 0);//lrbt
    DrawBoard_Ck_pp.SetMargin(4, 0, 0, 0, 0);//lrbt

    DrawBoard_Ck_pp.GetCanvas()->SetName(TString::Format("DrawBoard_Ck_pp"));
    DrawBoard_Ck_pp.cd(0);
    hAxis->Draw("AXIS");
    gAV18_G12.Draw("same,C");
    gRSC_G12.Draw("same,C");
    UnitLine->Draw("same");
    PT12->Draw("same");
    DrawBoard_Ck_pp.SetLabelSize(0, hAxis->GetXaxis(), 18);
    DrawBoard_Ck_pp.SetLabelSize(0, hAxis->GetYaxis(), 18);
    DrawBoard_Ck_pp.SetTitleSize(0, hAxis->GetXaxis(), 18);
    DrawBoard_Ck_pp.SetTitleSize(0, hAxis->GetYaxis(), 18);

    DrawBoard_Ck_pp.cd(1);
    hAxisRatio->Draw("AXIS");
    gRatio_G12.Draw("same,C");
    UnitLine->Draw("same");
    DrawBoard_Ck_pp.SetLabelSize(1, hAxisRatio->GetXaxis(), 18);
    DrawBoard_Ck_pp.SetLabelSize(1, hAxisRatio->GetYaxis(), 18);
    DrawBoard_Ck_pp.SetTitleSize(1, hAxisRatio->GetXaxis(), 18);
    DrawBoard_Ck_pp.SetTitleSize(1, hAxisRatio->GetYaxis(), 18);

    DrawBoard_Ck_pp.cd(2);
    hAxis->Draw("AXIS");
    gAV18_G08.Draw("same,C");
    gRSC_G08.Draw("same,C");
    UnitLine->Draw("same");
    PT08->Draw("same");
    DrawBoard_Ck_pp.SetLabelSize(0, hAxis->GetXaxis(), 18);
    DrawBoard_Ck_pp.SetLabelSize(0, hAxis->GetYaxis(), 18);
    DrawBoard_Ck_pp.SetTitleSize(0, hAxis->GetXaxis(), 18);
    DrawBoard_Ck_pp.SetTitleSize(0, hAxis->GetYaxis(), 18);

    DrawBoard_Ck_pp.cd(3);
    hAxisRatio->Draw("AXIS");
    gRatio_G08.Draw("same,C");
    UnitLine->Draw("same");
    DrawBoard_Ck_pp.SetLabelSize(1, hAxisRatio->GetXaxis(), 18);
    DrawBoard_Ck_pp.SetLabelSize(1, hAxisRatio->GetYaxis(), 18);
    DrawBoard_Ck_pp.SetTitleSize(1, hAxisRatio->GetXaxis(), 18);
    DrawBoard_Ck_pp.SetTitleSize(1, hAxisRatio->GetYaxis(), 18);


    TLegend* myLegend = new TLegend(0.12,0,0.95,0.96);//lbrt
    myLegend->SetName("myLegend");
    myLegend->SetTextSize(0.6);
    myLegend->SetNColumns(3);
    myLegend->AddEntry(&gAV18_G12,  "Argonne V18   ");
    myLegend->AddEntry(&gRSC_G12,   "Reid-Soft-Core");
    myLegend->AddEntry(&gRatio_G12, "Ratio to AV18   ");

    DrawBoard_Ck_pp.cd(4);
    myLegend->Draw();

    DrawBoard_Ck_pp.GetCanvas()->SaveAs(TString::Format("%sPlotCk_pp.png",OutputDir.Data()));
    hAxis->SetLineWidth(hAxis->GetLineWidth()/2.5);
    gAV18_G08.SetLineWidth(gAV18_G08.GetLineWidth()/2.5);
    gAV18_G12.SetLineWidth(gAV18_G12.GetLineWidth()/2.5);
    gRSC_G08.SetLineWidth(gRSC_G08.GetLineWidth()/2.5);
    gRSC_G12.SetLineWidth(gRSC_G12.GetLineWidth()/2.5);
    gRatio_G08.SetLineWidth(gRatio_G08.GetLineWidth()/2.5);
    gRatio_G12.SetLineWidth(gRatio_G12.GetLineWidth()/2.5);
    UnitLine->SetLineWidth(UnitLine->GetLineWidth()/2.5);
    DrawBoard_Ck_pp.GetCanvas()->SaveAs(TString::Format("%sPlotCk_pp.pdf",OutputDir.Data()));
    gAV18_G08.SetLineWidth(gAV18_G08.GetLineWidth()*2.5);
    gAV18_G12.SetLineWidth(gAV18_G12.GetLineWidth()*2.5);
    gRSC_G08.SetLineWidth(gRSC_G08.GetLineWidth()*2.5);
    gRSC_G12.SetLineWidth(gRSC_G12.GetLineWidth()*2.5);
    gRatio_G08.SetLineWidth(gRatio_G08.GetLineWidth()*2.5);
    gRatio_G12.SetLineWidth(gRatio_G12.GetLineWidth()*2.5);
    UnitLine->SetLineWidth(UnitLine->GetLineWidth()*2.5);

    delete hAxis;
    delete hAxisRatio;
    delete UnitLine;
    delete myLegend;
    delete PT12;
    delete PT08;
}

void CATSpaper::PlotCk_pp_ver2(){
    const double Plot_kMin = 0;
    const double Plot_kMax = 250;
    const unsigned NumPlotPts = 50;
    //CATShisto<double> CatHist(NumPlotPts,Plot_kMin,Plot_kMax);

    const unsigned NumRad = 3;
    const double Radius[NumRad] = {2.0,1.2,0.85};

    //double CkVal;
    //double Momentum;
    DLM_DtColor DlmCol;

    const unsigned NumScenarios = 3;
    TString ScenarioName[NumScenarios];
    ScenarioName[0] = "AV18 (s+p)";
    ScenarioName[1] = "Reid SC (s)";
    ScenarioName[2] = "Reid SC as in CRAB";

    TString ScenarioShortName[NumScenarios];
    ScenarioShortName[0] = "AV18sp";
    ScenarioShortName[1] = "RSCs";
    ScenarioShortName[2] = "RSCcrab";

    CATSparameters cSourcePars(CATSparameters::tSource,1,true);
    cSourcePars.SetParameters(&ppSourcePars[3]);

    CATS MyCat;
    MyCat.SetUseAnalyticSource(true);
    MyCat.SetAnaSource(GaussSource, cSourcePars);

    MyCat.SetExcludeFailedBins(false);
    MyCat.SetMomBins(NumPlotPts,Plot_kMin,Plot_kMax);

    MyCat.SetNumChannels(5);
    MyCat.SetNumPW(0,2);
    MyCat.SetNumPW(1,2);
    MyCat.SetNumPW(2,2);
    MyCat.SetNumPW(3,2);
    MyCat.SetNumPW(4,2);
    MyCat.SetSpin(0,0);
    MyCat.SetSpin(1,1);
    MyCat.SetSpin(2,1);
    MyCat.SetSpin(3,1);
    MyCat.SetSpin(4,1);//zero channel

    MyCat.SetQ1Q2(1);
    MyCat.SetPdgId(2212, 2212);
    MyCat.SetRedMass( 0.5*Mass_p );

    MyCat.SetShortRangePotential(0,0,fDlmPot,ppPars1S0);
    MyCat.SetShortRangePotential(1,1,fDlmPot,ppPars3P0);
    MyCat.SetShortRangePotential(2,1,fDlmPot,ppPars3P1);
    MyCat.SetShortRangePotential(3,1,fDlmPot,ppPars3P2);

    MyCat.KillTheCat();

    //in bins of Scenario/Radius
    TGraph** grCk;
    grCk = new TGraph* [NumScenarios];
    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        grCk[uSce] = new TGraph[NumRad];
        switch(uSce){
        //AV18 a la Vale
        case 0 :
            MyCat.SetChannelWeight(0, 3./12.);
            MyCat.SetChannelWeight(1, 1./12.);
            MyCat.SetChannelWeight(2, 3./12.);
            MyCat.SetChannelWeight(3, 5./12.);
            MyCat.SetChannelWeight(4, 0.);
            ppPars1S0[2] = NN_AV18;
            ppPars1S0[3] = v18_Coupled3P2;
            ppPars3P0[2] = NN_AV18;
            ppPars3P0[3] = v18_Coupled3P2;
            ppPars3P1[2] = NN_AV18;
            ppPars3P1[3] = v18_Coupled3P2;
            ppPars3P2[2] = NN_AV18;
            ppPars3P2[3] = v18_Coupled3P2;

            //ppPars1S0[2] = pp_ReidVale;
            //ppPars1S0[3] = 0;
            //ppPars3P0[2] = pp_ReidVale;
            //ppPars3P0[3] = 0;
            //ppPars3P1[2] = pp_ReidVale;
            //ppPars3P1[3] = 0;
            //ppPars3P2[2] = pp_ReidVale;
            //ppPars3P2[3] = 0;

            //MyCat.KillTheCat(CATS::kAllChanged);
            break;
        //AV18 s wave only
        case 1 :
            MyCat.SetChannelWeight(0, 1./4.);
            MyCat.SetChannelWeight(1, 0.);
            MyCat.SetChannelWeight(2, 0.);
            MyCat.SetChannelWeight(3, 0.);
            MyCat.SetChannelWeight(4, 3./4.);
            //MyCat.KillTheCat(CATS::kAllChanged);
            break;
        //Reid as in CRAB
        case 2 :
            MyCat.SetChannelWeight(0, 1./4.);
            MyCat.SetChannelWeight(1, 0.);
            MyCat.SetChannelWeight(2, 0.);
            MyCat.SetChannelWeight(3, 3./4.);
            MyCat.SetChannelWeight(4, 0.);
            ppPars1S0[2] = pp_ReidSC;
            ppPars1S0[3] = 0;
            ppPars3P0[2] = pp_ReidSC;
            ppPars3P0[3] = 0;
            ppPars3P1[2] = pp_ReidSC;
            ppPars3P1[3] = 0;
            ppPars3P2[2] = pp_ReidSC;
            ppPars3P2[3] = 0;
            //MyCat.KillTheCat(CATS::kAllChanged);
            break;
        default :
            printf("Default fuck!\n");
            break;
        }

        MyCat.KillTheCat(CATS::kAllChanged);

        for(unsigned uRad=0; uRad<NumRad; uRad++){
            MyCat.SetAnaSource(0, Radius[uRad]);
            MyCat.KillTheCat();

            grCk[uSce][uRad].Set(NumPlotPts);
            grCk[uSce][uRad].SetName(TString::Format("grCk_%u_%u",uSce,uRad));
            grCk[uSce][uRad].SetMarkerColor(DlmCol.GetColor(1));
            grCk[uSce][uRad].SetMarkerStyle(20);
            grCk[uSce][uRad].SetMarkerSize(0);
            grCk[uSce][uRad].SetLineColor(DlmCol.GetColor(1+uSce));
            grCk[uSce][uRad].SetLineWidth(8.-4.*double(uSce)/double(NumScenarios-1));
            grCk[uSce][uRad].SetFillColor(kWhite);

            for(unsigned uMomBin=0; uMomBin<NumPlotPts; uMomBin++){
                grCk[uSce][uRad].SetPoint(uMomBin, MyCat.GetMomentum(uMomBin), MyCat.GetCorrFun(uMomBin));
            }

        }
    }

    char* fBuffer = new char [512];
    FILE * pFile;
    pFile = fopen (TString::Format("%sPlotCk_pp_ver2.txt",OutputDir.Data()),"w");
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

    const double Yoffset = 0.9;
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

    hAxis->GetYaxis()->SetRangeUser(0, 6.2);
    hAxis->GetYaxis()->SetLimits(0, 6.2);

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
    hAxisZoom->GetYaxis()->SetTitle("C(k)");
    hAxisZoom->GetYaxis()->CenterTitle();
    hAxisZoom->GetYaxis()->SetTitleOffset(Yoffset/1.5);
    hAxisZoom->GetYaxis()->SetTitleSize(0.075);

    hAxisZoom->GetXaxis()->SetNdivisions(504);
    hAxisZoom->GetYaxis()->SetRangeUser(0.80, 1.05);
    hAxisZoom->GetYaxis()->SetLimits(0.80, 1.15);
    hAxisZoom->GetYaxis()->SetNdivisions(604);

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
        DrawBoard_Ck_pp.GetCanvas()->SetName(TString::Format("DrawBoard_Ck_pp_%u",uRad));
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
        DrawBoard_Ck_pp.GetCanvas()->SetName(TString::Format("DrawBoard_Ck_pp_%u",uRad));
        DrawBoard_Ck_pp.cd(uRad*2+1);
        hAxisZoom->Draw("AXIS");
        UnitLine->Draw("same");
        for(unsigned uSce=0; uSce<NumScenarios; uSce++) grCk[uSce][uRad].Draw("same,C");

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
    DrawBoard_Ck_pp.GetCanvas()->SaveAs(TString::Format("%sPlotCk_pp_ver2.png",OutputDir.Data()));
    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        for(unsigned uRad=0; uRad<NumRad; uRad++){
            grCk[uSce][uRad].SetLineWidth(grCk[uSce][uRad].GetLineWidth()/2.5);
        }
    }
    hAxis->SetLineWidth(hAxis->GetLineWidth()/2.5);
    UnitLine->SetLineWidth(UnitLine->GetLineWidth()/2.5);
    gStyle->SetLineWidth(1);
    DrawBoard_Ck_pp.GetCanvas()->SaveAs(TString::Format("%sPlotCk_pp_ver2.pdf",OutputDir.Data()));

    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        for(unsigned uRad=0; uRad<NumRad; uRad++){
            grCk[uSce][uRad].SetLineWidth(grCk[uSce][uRad].GetLineWidth()*2.5);
        }
    }
    hAxis->SetLineWidth(hAxis->GetLineWidth()*2.5);
    UnitLine->SetLineWidth(UnitLine->GetLineWidth()*2.5);

    delete hAxis;
    delete hAxisZoom;
    delete UnitLine;
    delete myLegend;

    for(unsigned uRad=0; uRad<NumRad; uRad++) delete myPT[uRad];
    delete [] myPT;

    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        delete [] grCk[uSce];
    }
    delete [] grCk;

    delete [] fBuffer;
}



void CATSpaper::BormioPlotCk_pp(){
    const double Plot_kMin = 0;
    const double Plot_kMax = 160;
    const double NumPlotPts = 40;
    //CATShisto<double> CatHist(NumPlotPts,Plot_kMin,Plot_kMax);

    const double RadLarge = 2;
    const double RadSmall = 1;

    const double ScaleResolution=2.2;
    //gStyle->SetLineWidth(gStyle->GetLineWidth()*ScaleResolution);

    double CkVal;
    double Momentum;
    DLM_DtColor DlmCol;

    CATSparameters cSourcePars(CATSparameters::tSource,1,true);
    cSourcePars.SetParameters(&ppSourcePars[3]);

    CATS MyCat;
    MyCat.SetUseAnalyticSource(true);
    MyCat.SetAnaSource(GaussSource, cSourcePars);

    MyCat.SetExcludeFailedBins(false);
    MyCat.SetMomBins(NumPlotPts,Plot_kMin,Plot_kMax);

    MyCat.SetNumChannels(4);
    MyCat.SetNumPW(0,2);
    MyCat.SetNumPW(1,2);
    MyCat.SetNumPW(2,2);
    MyCat.SetNumPW(3,2);
    MyCat.SetSpin(0,0);
    MyCat.SetSpin(1,1);
    MyCat.SetSpin(2,1);
    MyCat.SetSpin(3,1);
    MyCat.SetChannelWeight(0, 3./12.);
    MyCat.SetChannelWeight(1, 1./12.);
    MyCat.SetChannelWeight(2, 3./12.);
    MyCat.SetChannelWeight(3, 5./12.);
    MyCat.SetShortRangePotential(0,0,fDlmPot,ppPars1S0);
    MyCat.SetShortRangePotential(1,1,fDlmPot,ppPars3P0);
    MyCat.SetShortRangePotential(2,1,fDlmPot,ppPars3P1);
    MyCat.SetShortRangePotential(3,1,fDlmPot,ppPars3P2);
    MyCat.SetQ1Q2(1);
    MyCat.SetPdgId(2212, 2212);
    MyCat.SetRedMass( 0.5*Mass_p );

    ppSourcePars[3] = RadLarge;
    MyCat.KillTheCat();

    TGraph gAV18_G12;
    gAV18_G12.Set(NumPlotPts);
    gAV18_G12.SetName("gAV18_G12");
    gAV18_G12.SetMarkerColor(DlmCol.GetColor(1));
    gAV18_G12.SetMarkerStyle(20);
    gAV18_G12.SetMarkerSize(0);
    gAV18_G12.SetLineColor(DlmCol.GetColor(1));
    gAV18_G12.SetLineWidth(10*ScaleResolution);
    gAV18_G12.SetFillColor(kWhite);
    for(unsigned uBin=0; uBin<NumPlotPts; uBin++){
        gAV18_G12.SetPoint(uBin, MyCat.GetMomentum(uBin), MyCat.GetCorrFun(uBin));
    }

    MyCat.SetAnaSource(0, RadSmall);
    MyCat.KillTheCat();
    TGraph gAV18_G08;
    gAV18_G08.Set(NumPlotPts);
    gAV18_G08.SetName("gAV18_G08");
    gAV18_G08.SetMarkerColor(DlmCol.GetColor(1));
    gAV18_G08.SetMarkerStyle(20);
    gAV18_G08.SetMarkerSize(0);
    gAV18_G08.SetLineColor(DlmCol.GetColor(1));
    gAV18_G08.SetLineWidth(10*ScaleResolution);
    gAV18_G08.SetFillColor(kWhite);
    for(unsigned uBin=0; uBin<NumPlotPts; uBin++){
        gAV18_G08.SetPoint(uBin, MyCat.GetMomentum(uBin), MyCat.GetCorrFun(uBin));
    }

    ppPars1S0[2] = pp_ReidSC;
    ppPars1S0[3] = 0;
    ppPars3P0[2] = pp_ReidSC;
    ppPars3P0[3] = 0;
    ppPars3P1[2] = pp_ReidSC;
    ppPars3P1[3] = 0;
    ppPars3P2[2] = pp_ReidSC;
    ppPars3P2[3] = 0;
    MyCat.SetChannelWeight(0, 3./12.);
    MyCat.SetChannelWeight(1, 1./12.);
    MyCat.SetChannelWeight(2, 3./12.);
    MyCat.SetChannelWeight(3, 5./12.);

    MyCat.SetUseAnalyticSource(true);
    ppSourcePars[3] = RadLarge;
    MyCat.KillTheCat(CATS::kAllChanged);

    TGraph gRSC_G12;
    gRSC_G12.Set(NumPlotPts);
    gRSC_G12.SetName("gRSC_G12");
    gRSC_G12.SetMarkerColor(DlmCol.GetColor(2));
    gRSC_G12.SetMarkerStyle(20);
    gRSC_G12.SetMarkerSize(0);
    gRSC_G12.SetLineColor(DlmCol.GetColor(2));
    gRSC_G12.SetLineWidth(6*ScaleResolution);
    gRSC_G12.SetFillColor(kWhite);

    TGraph gRatio_G12;
    gRatio_G12.Set(NumPlotPts);
    gRatio_G12.SetName("gRatio_G12");
    gRatio_G12.SetMarkerColor(kBlack);
    gRatio_G12.SetMarkerStyle(20);
    gRatio_G12.SetMarkerSize(0);
    gRatio_G12.SetLineColor(kBlack);
    gRatio_G12.SetLineWidth(10*ScaleResolution);
    gRatio_G12.SetFillColor(kWhite);

    for(unsigned uBin=0; uBin<NumPlotPts; uBin++){
        gRSC_G12.SetPoint(uBin, MyCat.GetMomentum(uBin), MyCat.GetCorrFun(uBin));

        gAV18_G12.GetPoint(uBin, Momentum, CkVal);
        gRatio_G12.SetPoint(uBin, Momentum, MyCat.GetCorrFun(uBin)/CkVal);
    }

    MyCat.SetAnaSource(0, RadSmall);
    MyCat.KillTheCat();
    TGraph gRSC_G08;
    gRSC_G08.Set(NumPlotPts);
    gRSC_G08.SetName("gRSC_G08");
    gRSC_G08.SetMarkerColor(DlmCol.GetColor(2));
    gRSC_G08.SetMarkerStyle(20);
    gRSC_G08.SetMarkerSize(0);
    gRSC_G08.SetLineColor(DlmCol.GetColor(2));
    gRSC_G08.SetLineWidth(6*ScaleResolution);
    gRSC_G08.SetFillColor(kWhite);

    TGraph gRatio_G08;
    gRatio_G08.Set(NumPlotPts);
    gRatio_G08.SetName("gRatio_G08");
    gRatio_G08.SetMarkerColor(kBlack);
    gRatio_G08.SetMarkerStyle(20);
    gRatio_G08.SetMarkerSize(0);
    gRatio_G08.SetLineColor(kBlack);
    gRatio_G08.SetLineWidth(10*ScaleResolution);
    gRatio_G08.SetFillColor(kWhite);

    for(unsigned uBin=0; uBin<NumPlotPts; uBin++){
        gRSC_G08.SetPoint(uBin, MyCat.GetMomentum(uBin), MyCat.GetCorrFun(uBin));

        gAV18_G08.GetPoint(uBin, Momentum, CkVal);
        gRatio_G08.SetPoint(uBin, Momentum, MyCat.GetCorrFun(uBin)/CkVal);
    }

    const double Yoffset = 1;
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

    hAxis->GetYaxis()->SetRangeUser(0, 5.5);
    hAxis->GetYaxis()->SetLimits(0, 5.5);

    hAxis->GetXaxis()->SetNdivisions(504);

    TH1F* hAxisRatio = new TH1F("hAxisRatio", "hAxisRatio", NumPlotPts,Plot_kMin,Plot_kMax);
    hAxisRatio->SetStats(false);
    hAxisRatio->SetTitle("");
    hAxisRatio->GetXaxis()->SetLabelSize(0.065);
    hAxisRatio->GetXaxis()->SetTitle("k (MeV)");
    hAxisRatio->GetXaxis()->CenterTitle();
    hAxisRatio->GetXaxis()->SetTitleOffset(1.25);
    hAxisRatio->GetXaxis()->SetLabelOffset(0.02);
    hAxisRatio->GetXaxis()->SetTitleSize(0.075);
    hAxisRatio->GetYaxis()->SetLabelSize(0.065);
    hAxisRatio->GetYaxis()->SetTitle("Ratio");
    hAxisRatio->GetYaxis()->CenterTitle();
    hAxisRatio->GetYaxis()->SetTitleOffset(Yoffset/1.5);
    hAxisRatio->GetYaxis()->SetTitleSize(0.075);

    hAxisRatio->GetXaxis()->SetNdivisions(504);
    hAxisRatio->GetYaxis()->SetRangeUser(0.85, 1.15);
    hAxisRatio->GetYaxis()->SetLimits(0.85, 1.15);
    hAxisRatio->GetYaxis()->SetNdivisions(603);

    TPaveText* PT12 = new TPaveText(0.5,0.68,0.95,0.95, "blNDC");//lbrt
    PT12->SetName("PT12");
    PT12->SetBorderSize(1);
    PT12->SetTextSize(0.075);
    PT12->SetFillColor(kWhite);
    PT12->SetTextFont(22);
    PT12->AddText("Gauss source");
    PT12->AddText(TString::Format("R = %.2f fm", RadLarge));

    TPaveText* PT08 = new TPaveText(0.5,0.68,0.95,0.95, "blNDC");//lbrt
    PT08->SetName("PT08");
    PT08->SetBorderSize(1);
    PT08->SetTextSize(0.075);
    PT08->SetFillColor(kWhite);
    PT08->SetTextFont(22);
    PT08->AddText("Gauss source");
    PT08->AddText(TString::Format("R = %.2f fm", RadSmall));

    TF1* UnitLine = new TF1("UnitLine", "1",Plot_kMin,Plot_kMax);
    UnitLine->SetLineColor(kGray+1);
    UnitLine->SetLineWidth(6*ScaleResolution);

    DLM_SubPads DrawBoard_Ck_pp(4620,2599);
    DrawBoard_Ck_pp.AddSubPadTL(0,0.5,0.08,0.63);//lrtb
    DrawBoard_Ck_pp.AddSubPadTL(0,0.5,0.63,1);//lrtb
    DrawBoard_Ck_pp.AddSubPadTL(0.5,1,0.08,0.63);//lrtb
    DrawBoard_Ck_pp.AddSubPadTL(0.5,1,0.63,1);//lrtb
    DrawBoard_Ck_pp.AddSubPadTL(0,1,0,0.07);//lrtb

    DrawBoard_Ck_pp.SetMargin(0, 0.09, 0.04, 0, 0.04);//lrbt
    DrawBoard_Ck_pp.SetMargin(1, 0.09, 0.04, 0.14, 0);//lrbt
    DrawBoard_Ck_pp.SetMargin(2, 0.09, 0.04, 0, 0.04);//lrbt
    DrawBoard_Ck_pp.SetMargin(3, 0.09, 0.04, 0.14, 0);//lrbt
    DrawBoard_Ck_pp.SetMargin(4, 0, 0, 0, 0);//lrbt




    DrawBoard_Ck_pp.GetCanvas()->SetName(TString::Format("DrawBoard_Ck_pp"));
    DrawBoard_Ck_pp.cd(0);
    hAxis->Draw("AXIS");
    gAV18_G12.Draw("same,C");
    gRSC_G12.Draw("same,C");
    UnitLine->Draw("same");
    PT12->Draw("same");
    DrawBoard_Ck_pp.SetLabelSize(0, hAxis->GetXaxis(), 18);
    DrawBoard_Ck_pp.SetLabelSize(0, hAxis->GetYaxis(), 18);
    DrawBoard_Ck_pp.SetTitleSize(0, hAxis->GetXaxis(), 18);
    DrawBoard_Ck_pp.SetTitleSize(0, hAxis->GetYaxis(), 18);

    DrawBoard_Ck_pp.cd(1);
    hAxisRatio->Draw("AXIS");
    gRatio_G12.Draw("same,C");
    UnitLine->Draw("same");
    DrawBoard_Ck_pp.SetLabelSize(1, hAxisRatio->GetXaxis(), 18);
    DrawBoard_Ck_pp.SetLabelSize(1, hAxisRatio->GetYaxis(), 18);
    DrawBoard_Ck_pp.SetTitleSize(1, hAxisRatio->GetXaxis(), 18);
    DrawBoard_Ck_pp.SetTitleSize(1, hAxisRatio->GetYaxis(), 18);

    DrawBoard_Ck_pp.cd(2);
    hAxis->Draw("AXIS");
    gAV18_G08.Draw("same,C");
    gRSC_G08.Draw("same,C");
    UnitLine->Draw("same");
    PT08->Draw("same");
    DrawBoard_Ck_pp.SetLabelSize(0, hAxis->GetXaxis(), 18);
    DrawBoard_Ck_pp.SetLabelSize(0, hAxis->GetYaxis(), 18);
    DrawBoard_Ck_pp.SetTitleSize(0, hAxis->GetXaxis(), 18);
    DrawBoard_Ck_pp.SetTitleSize(0, hAxis->GetYaxis(), 18);

    DrawBoard_Ck_pp.cd(3);
    hAxisRatio->Draw("AXIS");
    gRatio_G08.Draw("same,C");
    UnitLine->Draw("same");
    DrawBoard_Ck_pp.SetLabelSize(1, hAxisRatio->GetXaxis(), 18);
    DrawBoard_Ck_pp.SetLabelSize(1, hAxisRatio->GetYaxis(), 18);
    DrawBoard_Ck_pp.SetTitleSize(1, hAxisRatio->GetXaxis(), 18);
    DrawBoard_Ck_pp.SetTitleSize(1, hAxisRatio->GetYaxis(), 18);


    TLegend* myLegend = new TLegend(0.12,0,0.95,0.96);//lbrt
    myLegend->SetName("myLegend");
    myLegend->SetTextSize(0.6);
    myLegend->SetNColumns(3);
    myLegend->AddEntry(&gAV18_G12,  "Argonne V18   ");
    myLegend->AddEntry(&gRSC_G12,   "Reid-Soft-Core");
    myLegend->AddEntry(&gRatio_G12, "Ratio to AV18   ");

    DrawBoard_Ck_pp.cd(4);
    myLegend->Draw();

    DrawBoard_Ck_pp.GetCanvas()->SaveAs(TString::Format("%sBormioPlotCk_pp.png",OutputDir.Data()));
    hAxis->SetLineWidth(hAxis->GetLineWidth()/2.5);
    gAV18_G08.SetLineWidth(gAV18_G08.GetLineWidth()/2.5);
    gAV18_G12.SetLineWidth(gAV18_G12.GetLineWidth()/2.5);
    gRSC_G08.SetLineWidth(gRSC_G08.GetLineWidth()/2.5);
    gRSC_G12.SetLineWidth(gRSC_G12.GetLineWidth()/2.5);
    gRatio_G08.SetLineWidth(gRatio_G08.GetLineWidth()/2.5);
    gRatio_G12.SetLineWidth(gRatio_G12.GetLineWidth()/2.5);
    UnitLine->SetLineWidth(UnitLine->GetLineWidth()/2.5);
    DrawBoard_Ck_pp.GetCanvas()->SaveAs(TString::Format("%sBormioPlotCk_pp.pdf",OutputDir.Data()));
    gAV18_G08.SetLineWidth(gAV18_G08.GetLineWidth()*2.5);
    gAV18_G12.SetLineWidth(gAV18_G12.GetLineWidth()*2.5);
    gRSC_G08.SetLineWidth(gRSC_G08.GetLineWidth()*2.5);
    gRSC_G12.SetLineWidth(gRSC_G12.GetLineWidth()*2.5);
    gRatio_G08.SetLineWidth(gRatio_G08.GetLineWidth()*2.5);
    gRatio_G12.SetLineWidth(gRatio_G12.GetLineWidth()*2.5);
    UnitLine->SetLineWidth(UnitLine->GetLineWidth()*2.5);

    delete hAxis;
    delete hAxisRatio;
    delete UnitLine;
    delete myLegend;
    delete PT12;
    delete PT08;
}

void CATSpaper::PlotCk_pL(){
    const double Plot_kMin = 5;
    const double Plot_kMax = 245;
    const double NumPlotPts = 24;

    const double RadLarge = 1.2;
    const double RadSmall = 0.85;

    double CkVal;
    double Momentum;
    DLM_DtColor DlmCol;

    CATSparameters cSourcePars(CATSparameters::tSource,1,true);
    cSourcePars.SetParameters(&pLamSourcePars[3]);

    CATS MyCat;
    MyCat.SetUseAnalyticSource(true);
    MyCat.SetAnaSource(GaussSource, cSourcePars);

    MyCat.SetExcludeFailedBins(false);
    MyCat.SetMomBins(NumPlotPts,Plot_kMin,Plot_kMax);

    MyCat.SetNumChannels(2);
    MyCat.SetNumPW(0,1);
    MyCat.SetNumPW(1,1);
    MyCat.SetSpin(0,0);
    MyCat.SetSpin(1,1);
    MyCat.SetChannelWeight(0, 1./4.);
    MyCat.SetChannelWeight(1, 3./4.);
    MyCat.SetShortRangePotential(0,0,fDlmPot,pLamPars1S0);
    MyCat.SetShortRangePotential(1,0,fDlmPot,pLamPars3S1);
    MyCat.SetQ1Q2(0);
    MyCat.SetPdgId(2212, 3122);
    MyCat.SetRedMass( (Mass_p*Mass_L)/(Mass_p+Mass_L) );

    pLamSourcePars[3] = RadLarge;
    MyCat.KillTheCat();

    TGraph gUsm_G12;
    gUsm_G12.Set(NumPlotPts);
    gUsm_G12.SetName("gUsm_G12");
    gUsm_G12.SetMarkerColor(DlmCol.GetColor(2));
    gUsm_G12.SetMarkerStyle(20);
    gUsm_G12.SetMarkerSize(0);
    gUsm_G12.SetLineColor(DlmCol.GetColor(2));
    gUsm_G12.SetLineWidth(6);
    gUsm_G12.SetFillColor(kWhite);
    for(unsigned uBin=0; uBin<NumPlotPts; uBin++){
        gUsm_G12.SetPoint(uBin, MyCat.GetMomentum(uBin), MyCat.GetCorrFun(uBin));
    }

    MyCat.SetAnaSource(0, RadSmall);
    MyCat.KillTheCat();
    TGraph gUsm_G08;
    gUsm_G08.Set(NumPlotPts);
    gUsm_G08.SetName("gUsm_G08");
    gUsm_G08.SetMarkerColor(DlmCol.GetColor(2));
    gUsm_G08.SetMarkerStyle(20);
    gUsm_G08.SetMarkerSize(0);
    gUsm_G08.SetLineColor(DlmCol.GetColor(2));
    gUsm_G08.SetLineWidth(6);
    gUsm_G08.SetFillColor(kWhite);
    for(unsigned uBin=0; uBin<NumPlotPts; uBin++){
        gUsm_G08.SetPoint(uBin, MyCat.GetMomentum(uBin), MyCat.GetCorrFun(uBin));
    }


    complex<double>**** WaveFunctionU; double*** PhaseShifts; double* RadBins; unsigned NumRadBins=0;
    InitHaidenbauerNLO("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/Haidenbauer/pLambdaNLO/",MyCat,&WaveFunctionU,&PhaseShifts,&RadBins,NumRadBins);

    for(unsigned uBin=0; uBin<NumPlotPts; uBin++){
        MyCat.UseExternalWaveFunction(uBin,0,0,WaveFunctionU[uBin][0][0], NumRadBins, RadBins, PhaseShifts[uBin][0][0]);
        MyCat.UseExternalWaveFunction(uBin,1,0,WaveFunctionU[uBin][1][0], NumRadBins, RadBins, PhaseShifts[uBin][1][0]);
    }

    pLamSourcePars[3] = RadLarge;
    MyCat.KillTheCat(CATS::kAllChanged);

    TGraph gHaide_G12;
    gHaide_G12.Set(NumPlotPts);
    gHaide_G12.SetName("gHaide_G12");
    gHaide_G12.SetMarkerColor(DlmCol.GetColor(1));
    gHaide_G12.SetMarkerStyle(20);
    gHaide_G12.SetMarkerSize(0);
    gHaide_G12.SetLineColor(DlmCol.GetColor(1));
    gHaide_G12.SetLineWidth(10);
    gHaide_G12.SetFillColor(kWhite);

    TGraph gRatio_G12;
    gRatio_G12.Set(NumPlotPts);
    gRatio_G12.SetName("gRatio_G12");
    gRatio_G12.SetMarkerColor(kBlack);
    gRatio_G12.SetMarkerStyle(20);
    gRatio_G12.SetMarkerSize(0);
    gRatio_G12.SetLineColor(kBlack);
    gRatio_G12.SetLineWidth(10);
    gRatio_G12.SetFillColor(kWhite);

    for(unsigned uBin=0; uBin<NumPlotPts; uBin++){
        gHaide_G12.SetPoint(uBin, MyCat.GetMomentum(uBin), MyCat.GetCorrFun(uBin));

        gUsm_G12.GetPoint(uBin, Momentum, CkVal);
        gRatio_G12.SetPoint(uBin, Momentum, CkVal/MyCat.GetCorrFun(uBin));
    }

    MyCat.SetAnaSource(0, RadSmall);
    MyCat.KillTheCat();
    TGraph gHaide_G08;
    gHaide_G08.Set(NumPlotPts);
    gHaide_G08.SetName("gHaide_G08");
    gHaide_G08.SetMarkerColor(DlmCol.GetColor(1));
    gHaide_G08.SetMarkerStyle(20);
    gHaide_G08.SetMarkerSize(0);
    gHaide_G08.SetLineColor(DlmCol.GetColor(1));
    gHaide_G08.SetLineWidth(10);
    gHaide_G08.SetFillColor(kWhite);

    TGraph gRatio_G08;
    gRatio_G08.Set(NumPlotPts);
    gRatio_G08.SetName("gRatio_G08");
    gRatio_G08.SetMarkerColor(kBlack);
    gRatio_G08.SetMarkerStyle(20);
    gRatio_G08.SetMarkerSize(0);
    gRatio_G08.SetLineColor(kBlack);
    gRatio_G08.SetLineWidth(10);
    gRatio_G08.SetFillColor(kWhite);

    for(unsigned uBin=0; uBin<NumPlotPts; uBin++){
        gHaide_G08.SetPoint(uBin, MyCat.GetMomentum(uBin), MyCat.GetCorrFun(uBin));

        gUsm_G08.GetPoint(uBin, Momentum, CkVal);
        gRatio_G08.SetPoint(uBin, Momentum, CkVal/MyCat.GetCorrFun(uBin));
    }

    const double Yoffset = 1;
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

    hAxis->GetYaxis()->SetRangeUser(0.5, 4);
    hAxis->GetYaxis()->SetLimits(0.5, 4);

    hAxis->GetXaxis()->SetNdivisions(505);

    TH1F* hAxisRatio = new TH1F("hAxisRatio", "hAxisRatio", NumPlotPts,Plot_kMin,Plot_kMax);
    hAxisRatio->SetStats(false);
    hAxisRatio->SetTitle("");
    hAxisRatio->GetXaxis()->SetLabelSize(0.065);
    hAxisRatio->GetXaxis()->SetTitle("k (MeV)");
    hAxisRatio->GetXaxis()->CenterTitle();
    hAxisRatio->GetXaxis()->SetTitleOffset(1.25);
    hAxisRatio->GetXaxis()->SetLabelOffset(0.02);
    hAxisRatio->GetXaxis()->SetTitleSize(0.075);
    hAxisRatio->GetYaxis()->SetLabelSize(0.065);
    hAxisRatio->GetYaxis()->SetTitle("Ratio");
    hAxisRatio->GetYaxis()->CenterTitle();
    hAxisRatio->GetYaxis()->SetTitleOffset(Yoffset/1.5);
    hAxisRatio->GetYaxis()->SetTitleSize(0.075);

    hAxisRatio->GetXaxis()->SetNdivisions(505);
    hAxisRatio->GetYaxis()->SetRangeUser(0.85, 1.09);
    hAxisRatio->GetYaxis()->SetLimits(0.85, 1.09);
    hAxisRatio->GetYaxis()->SetNdivisions(604);

    TPaveText* PT12 = new TPaveText(0.5,0.68,0.95,0.95, "blNDC");//lbrt
    PT12->SetName("PT12");
    PT12->SetBorderSize(1);
    PT12->SetTextSize(0.075);
    PT12->SetFillColor(kWhite);
    PT12->SetTextFont(22);
    PT12->AddText("Gauss source");
    PT12->AddText(TString::Format("R = %.2f fm", RadLarge));

    TPaveText* PT08 = new TPaveText(0.5,0.68,0.95,0.95, "blNDC");//lbrt
    PT08->SetName("PT08");
    PT08->SetBorderSize(1);
    PT08->SetTextSize(0.075);
    PT08->SetFillColor(kWhite);
    PT08->SetTextFont(22);
    PT08->AddText("Gauss source");
    PT08->AddText(TString::Format("R = %.2f fm", RadSmall));

    TF1* UnitLine = new TF1("UnitLine", "1",Plot_kMin,Plot_kMax);
    UnitLine->SetLineColor(kGray+1);
    UnitLine->SetLineWidth(6);

    DLM_SubPads DrawBoard_Ck_pp(1920,1080);
    DrawBoard_Ck_pp.AddSubPadTL(0,0.5,0.08,0.63);//lrtb
    DrawBoard_Ck_pp.AddSubPadTL(0,0.5,0.63,1);//lrtb
    DrawBoard_Ck_pp.AddSubPadTL(0.5,1,0.08,0.63);//lrtb
    DrawBoard_Ck_pp.AddSubPadTL(0.5,1,0.63,1);//lrtb
    DrawBoard_Ck_pp.AddSubPadTL(0,1,0,0.07);//lrtb

    DrawBoard_Ck_pp.SetMargin(0, 0.09, 0.04, 0, 0.04);//lrbt
    DrawBoard_Ck_pp.SetMargin(1, 0.09, 0.04, 0.14, 0);//lrbt
    DrawBoard_Ck_pp.SetMargin(2, 0.09, 0.04, 0, 0.04);//lrbt
    DrawBoard_Ck_pp.SetMargin(3, 0.09, 0.04, 0.14, 0);//lrbt
    DrawBoard_Ck_pp.SetMargin(4, 0, 0, 0, 0);//lrbt

    DrawBoard_Ck_pp.GetCanvas()->SetName(TString::Format("DrawBoard_Ck_pp"));
    DrawBoard_Ck_pp.cd(0);
    hAxis->Draw("AXIS");
    gHaide_G12.Draw("same,C");
    gUsm_G12.Draw("same,C");
    UnitLine->Draw("same");
    PT12->Draw("same");
    DrawBoard_Ck_pp.SetLabelSize(0, hAxis->GetXaxis(), 18);
    DrawBoard_Ck_pp.SetLabelSize(0, hAxis->GetYaxis(), 18);
    DrawBoard_Ck_pp.SetTitleSize(0, hAxis->GetXaxis(), 18);
    DrawBoard_Ck_pp.SetTitleSize(0, hAxis->GetYaxis(), 18);

    DrawBoard_Ck_pp.cd(1);
    hAxisRatio->Draw("AXIS");
    gRatio_G12.Draw("same,C");
    UnitLine->Draw("same");
    DrawBoard_Ck_pp.SetLabelSize(1, hAxisRatio->GetXaxis(), 18);
    DrawBoard_Ck_pp.SetLabelSize(1, hAxisRatio->GetYaxis(), 18);
    DrawBoard_Ck_pp.SetTitleSize(1, hAxisRatio->GetXaxis(), 18);
    DrawBoard_Ck_pp.SetTitleSize(1, hAxisRatio->GetYaxis(), 18);

    DrawBoard_Ck_pp.cd(2);
    hAxis->Draw("AXIS");
    gHaide_G08.Draw("same,C");
    gUsm_G08.Draw("same,C");
    UnitLine->Draw("same");
    PT08->Draw("same");
    DrawBoard_Ck_pp.SetLabelSize(0, hAxis->GetXaxis(), 18);
    DrawBoard_Ck_pp.SetLabelSize(0, hAxis->GetYaxis(), 18);
    DrawBoard_Ck_pp.SetTitleSize(0, hAxis->GetXaxis(), 18);
    DrawBoard_Ck_pp.SetTitleSize(0, hAxis->GetYaxis(), 18);

    DrawBoard_Ck_pp.cd(3);
    hAxisRatio->Draw("AXIS");
    gRatio_G08.Draw("same,C");
    UnitLine->Draw("same");
    DrawBoard_Ck_pp.SetLabelSize(1, hAxisRatio->GetXaxis(), 18);
    DrawBoard_Ck_pp.SetLabelSize(1, hAxisRatio->GetYaxis(), 18);
    DrawBoard_Ck_pp.SetTitleSize(1, hAxisRatio->GetXaxis(), 18);
    DrawBoard_Ck_pp.SetTitleSize(1, hAxisRatio->GetYaxis(), 18);


    TLegend* myLegend = new TLegend(0.12,0,0.95,0.96);//lbrt
    myLegend->SetName("myLegend");
    myLegend->SetTextSize(0.6);
    myLegend->SetNColumns(3);
    myLegend->AddEntry(&gHaide_G12,   "NLO   ");
    myLegend->AddEntry(&gUsm_G12,  "Usmani   ");
    myLegend->AddEntry(&gRatio_G12, "Ratio to NLO   ");

    DrawBoard_Ck_pp.cd(4);
    myLegend->Draw();

    DrawBoard_Ck_pp.GetCanvas()->SaveAs(TString::Format("%sPlotCk_pL.png",OutputDir.Data()));
    hAxis->SetLineWidth(hAxis->GetLineWidth()/2.5);
    gUsm_G08.SetLineWidth(gUsm_G08.GetLineWidth()/2.5);
    gUsm_G12.SetLineWidth(gUsm_G12.GetLineWidth()/2.5);
    gHaide_G08.SetLineWidth(gHaide_G08.GetLineWidth()/2.5);
    gHaide_G12.SetLineWidth(gHaide_G12.GetLineWidth()/2.5);
    gRatio_G08.SetLineWidth(gRatio_G08.GetLineWidth()/2.5);
    gRatio_G12.SetLineWidth(gRatio_G12.GetLineWidth()/2.5);
    UnitLine->SetLineWidth(UnitLine->GetLineWidth()/2.5);
    DrawBoard_Ck_pp.GetCanvas()->SaveAs(TString::Format("%sPlotCk_pL.pdf",OutputDir.Data()));
    gUsm_G08.SetLineWidth(gUsm_G08.GetLineWidth()*2.5);
    gUsm_G12.SetLineWidth(gUsm_G12.GetLineWidth()*2.5);
    gHaide_G08.SetLineWidth(gHaide_G08.GetLineWidth()*2.5);
    gHaide_G12.SetLineWidth(gHaide_G12.GetLineWidth()*2.5);
    gRatio_G08.SetLineWidth(gRatio_G08.GetLineWidth()*2.5);
    gRatio_G12.SetLineWidth(gRatio_G12.GetLineWidth()*2.5);
    UnitLine->SetLineWidth(UnitLine->GetLineWidth()*2.5);

    delete hAxis;
    delete hAxisRatio;
    delete UnitLine;
    delete myLegend;
    delete PT12;
    delete PT08;

    CleanHaidenbauer(MyCat,&WaveFunctionU,&PhaseShifts,&RadBins);
}


void CATSpaper::PlotCk_pL_ver2(){
    const double Plot_kMin = 5;
    const double Plot_kMax = 245;
    const double NumPlotPts = 24;

    const unsigned NumRad = 3;
    const double Radius[NumRad] = {2.0,1.2,0.85};

    double CkVal0;
    double Momentum0;
    double CkVal;
    double Momentum;
    DLM_DtColor DlmCol;

    const unsigned NumScenarios = 3;
    TString ScenarioName[NumScenarios];
    ScenarioName[0] = "NLO";
    ScenarioName[1] = "Usmani";
    ScenarioName[2] = "NLO (Lednicky)";

    TString ScenarioShortName[NumScenarios];
    ScenarioShortName[0] = "NLO";
    ScenarioShortName[1] = "Usmani";
    ScenarioShortName[2] = "LedniNLO";

    CATSparameters cSourcePars(CATSparameters::tSource,1,true);
    cSourcePars.SetParameters(&pLamSourcePars[3]);

    CATS MyCat;
    MyCat.SetUseAnalyticSource(true);
    MyCat.SetAnaSource(GaussSource, cSourcePars);

    MyCat.SetExcludeFailedBins(false);
    MyCat.SetMomBins(NumPlotPts,Plot_kMin,Plot_kMax);

    MyCat.SetNumChannels(2);
    MyCat.SetNumPW(0,1);
    MyCat.SetNumPW(1,1);
    MyCat.SetSpin(0,0);
    MyCat.SetSpin(1,1);

    MyCat.SetQ1Q2(0);
    MyCat.SetPdgId(2212, 3122);
    MyCat.SetRedMass( (Mass_p*Mass_L)/(Mass_p+Mass_L) );

    MyCat.SetChannelWeight(0, 1./4.);
    MyCat.SetChannelWeight(1, 3./4.);


    //const double UsmaniScatPar[4] = {2.88,2.92,1.66,3.78};
    const double NloScatPar[4] = {2.91,2.78,1.53,2.72};

    complex<double>**** WaveFunctionU; double*** PhaseShifts; double* RadBins; unsigned NumRadBins=0;
    InitHaidenbauerNLO("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/Haidenbauer/pLambdaNLO/",MyCat,&WaveFunctionU,&PhaseShifts,&RadBins,NumRadBins);

    //in bins of Scenario/Radius
    TGraph** grCk;
    grCk = new TGraph* [NumScenarios];
    TGraph** grCkRatio;
    grCkRatio = new TGraph* [NumScenarios];
    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        grCk[uSce] = new TGraph[NumRad];
        grCkRatio[uSce] = new TGraph[NumRad];
        switch(uSce){
        //NLO
        case 0 :
            for(unsigned uBin=0; uBin<NumPlotPts; uBin++){
                MyCat.UseExternalWaveFunction(uBin,0,0,WaveFunctionU[uBin][0][0], NumRadBins, RadBins, PhaseShifts[uBin][0][0]);
                MyCat.UseExternalWaveFunction(uBin,1,0,WaveFunctionU[uBin][1][0], NumRadBins, RadBins, PhaseShifts[uBin][1][0]);
            }
            break;
        //Usmani
        case 1 :
            for(unsigned uBin=0; uBin<NumPlotPts; uBin++){
                MyCat.UseExternalWaveFunction(uBin,0,0);
                MyCat.UseExternalWaveFunction(uBin,1,0);
            }
            MyCat.SetShortRangePotential(0,0,fDlmPot,pLamPars1S0);
            MyCat.SetShortRangePotential(1,0,fDlmPot,pLamPars3S1);
            break;
        //Ledni NLO
        case 2 :
            break;
        default :
            printf("Default pΛ fuck!\n");
            break;
        }

        if(uSce<2) MyCat.KillTheCat(CATS::kAllChanged);
        for(unsigned uRad=0; uRad<NumRad; uRad++){
            if(uSce<2){
                MyCat.SetAnaSource(0, Radius[uRad]);
                MyCat.KillTheCat();
            }

            grCk[uSce][uRad].Set(NumPlotPts);
            grCk[uSce][uRad].SetName(TString::Format("grCk_%u_%u",uSce,uRad));
            grCk[uSce][uRad].SetMarkerColor(DlmCol.GetColor(1));
            grCk[uSce][uRad].SetMarkerStyle(20);
            grCk[uSce][uRad].SetMarkerSize(0);
            grCk[uSce][uRad].SetLineColor(DlmCol.GetColor(1+uSce));
            grCk[uSce][uRad].SetLineWidth(8.-2.*double(uSce)/double(NumScenarios-1));
            grCk[uSce][uRad].SetFillColor(kWhite);

            for(unsigned uMomBin=0; uMomBin<NumPlotPts; uMomBin++){
                if(uSce<2) grCk[uSce][uRad].SetPoint(uMomBin, MyCat.GetMomentum(uMomBin), MyCat.GetCorrFun(uMomBin));
                else{
                    Momentum = MyCat.GetMomentum(uMomBin);
                    grCk[uSce][uRad].SetPoint(uMomBin, Momentum, Lednicky_SingletTriplet(Momentum, &Radius[uRad], NloScatPar));
                }
            }

        }
    }

    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        for(unsigned uRad=0; uRad<NumRad; uRad++){
            grCkRatio[uSce][uRad].Set(NumPlotPts);
            grCkRatio[uSce][uRad].SetName(TString::Format("grCkRatio_%u_%u",uSce,uRad));
            grCkRatio[uSce][uRad].SetMarkerColor(DlmCol.GetColor(1+uSce));
            grCkRatio[uSce][uRad].SetMarkerStyle(20);
            grCkRatio[uSce][uRad].SetMarkerSize(0);
            grCkRatio[uSce][uRad].SetLineColor(DlmCol.GetColor(1+uSce));
            grCkRatio[uSce][uRad].SetLineWidth(8.-2.*double(uSce)/double(NumScenarios-1));
            grCkRatio[uSce][uRad].SetLineStyle(2.);
            grCkRatio[uSce][uRad].SetFillColor(kWhite);
            for(unsigned uMomBin=0; uMomBin<NumPlotPts; uMomBin++){
                grCk[0][uRad].GetPoint(uMomBin,Momentum0,CkVal0);
                grCk[uSce][uRad].GetPoint(uMomBin,Momentum,CkVal);
                grCkRatio[uSce][uRad].SetPoint(uMomBin, Momentum, CkVal/CkVal0);
            }
        }
    }

    char* fBuffer = new char [512];
    FILE * pFile;
    pFile = fopen (TString::Format("%sPlotCk_pL_ver2.txt",OutputDir.Data()),"w");
    fprintf (pFile, "%18s","k (MeV)");
    for(unsigned uRad=0; uRad<NumRad; uRad++){
        for(unsigned uSce=0; uSce<NumScenarios; uSce++){
            sprintf(fBuffer,"R%.2f_%s",Radius[uRad],ScenarioShortName[uSce].Data());
            fprintf (pFile,"%18s",fBuffer);
            fprintf (pFile,"%18s","<-RatioNLO");
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
                grCkRatio[uSce][uRad].GetPoint(uMomBin,xVal,yVal);
                fprintf (pFile, "%18.3f",yVal);
            }
        }
    }


    TGraph grDummyRatio;
    grDummyRatio.SetName(TString::Format("grDummyRatio"));
    grDummyRatio.SetMarkerColor(kBlack);
    grDummyRatio.SetMarkerStyle(20);
    grDummyRatio.SetMarkerSize(0);
    grDummyRatio.SetLineColor(kBlack);
    grDummyRatio.SetLineWidth(8.);
    grDummyRatio.SetLineStyle(2.);
    grDummyRatio.SetFillColor(kWhite);


    const double Yoffset = 0.9;
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

    hAxis->GetYaxis()->SetRangeUser(0.5, 4);
    hAxis->GetYaxis()->SetLimits(0.5, 4);

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
    hAxisZoom->GetYaxis()->SetRangeUser(0.9, 1.15);
    hAxisZoom->GetYaxis()->SetLimits(0.9, 1.15);
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

    DLM_SubPads DrawBoard_Ck_pL(2560,1080);
    double PadLength = 1./double(NumRad);
    double CurrentPos = 0;
    for(unsigned uRad=0; uRad<NumRad; uRad++){
        DrawBoard_Ck_pL.AddSubPadTL(CurrentPos,CurrentPos+PadLength,0.08,0.63);//lrtb
        DrawBoard_Ck_pL.AddSubPadTL(CurrentPos,CurrentPos+PadLength,0.63,1);//lrtb
        CurrentPos+=PadLength;
//printf("uRad=%u; %.2f <--> %.2f\n", uRad, );
    }
    DrawBoard_Ck_pL.AddSubPadTL(0,1,0,0.075);//lrtb

    for(unsigned uRad=0; uRad<NumRad; uRad++){
        DrawBoard_Ck_pL.SetMargin(uRad*2, 0.07, 0.01, 0, 0.04);//lrbt
        DrawBoard_Ck_pL.SetMargin(uRad*2+1, 0.07, 0.01, 0.14, 0);//lrbt
    }
    DrawBoard_Ck_pL.SetMargin(2*NumRad, 0, 0, 0, 0);//lrbt

    for(unsigned uRad=0; uRad<NumRad; uRad++){
        DrawBoard_Ck_pL.GetCanvas()->SetName(TString::Format("DrawBoard_Ck_pL_%u",uRad));
        DrawBoard_Ck_pL.cd(uRad*2);
        hAxis->Draw("AXIS");
        UnitLine->Draw("same");
        for(unsigned uSce=0; uSce<NumScenarios; uSce++) grCk[uSce][uRad].Draw("same,C");
        for(unsigned uSce=0; uSce<NumScenarios; uSce++) myPT[uRad]->Draw("same");

        DrawBoard_Ck_pL.SetLabelSize(uRad*2, hAxis->GetXaxis(), 17);
        DrawBoard_Ck_pL.SetLabelSize(uRad*2, hAxis->GetYaxis(), 17);
        DrawBoard_Ck_pL.SetTitleSize(uRad*2, hAxis->GetXaxis(), 17);
        DrawBoard_Ck_pL.SetTitleSize(uRad*2, hAxis->GetYaxis(), 17);
    }
    for(unsigned uRad=0; uRad<NumRad; uRad++){
        DrawBoard_Ck_pL.GetCanvas()->SetName(TString::Format("DrawBoard_Ck_pL_%u",uRad));
        DrawBoard_Ck_pL.cd(uRad*2+1);
        hAxisZoom->Draw("AXIS");
        UnitLine->Draw("same");
        for(unsigned uSce=0; uSce<NumScenarios; uSce++) if(uSce!=0) grCkRatio[uSce][uRad].Draw("same,C");
        //grCkRatio[1][uRad].Draw("same,C");

        DrawBoard_Ck_pL.SetLabelSize(uRad*2+1, hAxisZoom->GetXaxis(), 17);
        DrawBoard_Ck_pL.SetLabelSize(uRad*2+1, hAxisZoom->GetYaxis(), 17);
        DrawBoard_Ck_pL.SetTitleSize(uRad*2+1, hAxisZoom->GetXaxis(), 17);
        DrawBoard_Ck_pL.SetTitleSize(uRad*2+1, hAxisZoom->GetYaxis(), 17);
    }

    TLegend* myLegend = new TLegend(0.12,0,0.95,0.96);//lbrt
    myLegend->SetName("myLegend");
    myLegend->SetTextSize(0.65);
    myLegend->SetNColumns(NumScenarios+1);
    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        myLegend->AddEntry(&grCk[uSce][0],  ScenarioName[uSce]);
    }
    myLegend->AddEntry(&grDummyRatio,  TString::Format("Ratios to %s",ScenarioName[0].Data()));

    DrawBoard_Ck_pL.cd(NumRad*2);
    myLegend->Draw();

    gStyle->SetLineWidth(2.5);
    DrawBoard_Ck_pL.GetCanvas()->SaveAs(TString::Format("%sPlotCk_pL_ver2.png",OutputDir.Data()));
    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        for(unsigned uRad=0; uRad<NumRad; uRad++){
            grCk[uSce][uRad].SetLineWidth(grCk[uSce][uRad].GetLineWidth()/2.5);
            grCkRatio[uSce][uRad].SetLineWidth(grCkRatio[uSce][uRad].GetLineWidth()/2.5);
        }
    }
    hAxis->SetLineWidth(hAxis->GetLineWidth()/2.5);
    UnitLine->SetLineWidth(UnitLine->GetLineWidth()/2.5);
    gStyle->SetLineWidth(1);
    DrawBoard_Ck_pL.GetCanvas()->SaveAs(TString::Format("%sPlotCk_pL_ver2.pdf",OutputDir.Data()));

    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        for(unsigned uRad=0; uRad<NumRad; uRad++){
            grCk[uSce][uRad].SetLineWidth(grCk[uSce][uRad].GetLineWidth()*2.5);
            grCkRatio[uSce][uRad].SetLineWidth(grCkRatio[uSce][uRad].GetLineWidth()*2.5);
        }
    }
    hAxis->SetLineWidth(hAxis->GetLineWidth()*2.5);
    UnitLine->SetLineWidth(UnitLine->GetLineWidth()*2.5);

    delete hAxis;
    delete hAxisZoom;
    delete UnitLine;
    delete myLegend;

    for(unsigned uRad=0; uRad<NumRad; uRad++) delete myPT[uRad];
    delete [] myPT;

    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        delete [] grCk[uSce];
        delete [] grCkRatio[uSce];
    }
    delete [] grCk;
    delete [] grCkRatio;
    delete [] fBuffer;

    CleanHaidenbauer(MyCat,&WaveFunctionU,&PhaseShifts,&RadBins);
}


void CATSpaper::BormioPlotCk_pL(){
    const double Plot_kMin = 0;
    const double Plot_kMax = 200;
    const double NumPlotPts = 50;

    const double RadLarge = 2;
    const double RadSmall = 1;

    const double ScaleResolution = 2.2;
    //gStyle->SetLineWidth(gStyle->GetLineWidth()*ScaleResolution);

    double CkVal;
    double Momentum;
    DLM_DtColor DlmCol;

    CATSparameters cSourcePars(CATSparameters::tSource,1,true);
    cSourcePars.SetParameters(&pLamSourcePars[3]);

    CATS MyCat;
    MyCat.SetUseAnalyticSource(true);
    MyCat.SetAnaSource(GaussSource, cSourcePars);

    MyCat.SetExcludeFailedBins(false);
    MyCat.SetMomBins(NumPlotPts,Plot_kMin,Plot_kMax);

    MyCat.SetNumChannels(2);
    MyCat.SetNumPW(0,1);
    MyCat.SetNumPW(1,1);
    MyCat.SetSpin(0,0);
    MyCat.SetSpin(1,1);
    MyCat.SetChannelWeight(0, 1./4.);
    MyCat.SetChannelWeight(1, 3./4.);
    MyCat.SetShortRangePotential(0,0,fDlmPot,pLamPars1S0);
    MyCat.SetShortRangePotential(1,0,fDlmPot,pLamPars3S1);
    MyCat.SetQ1Q2(0);
    MyCat.SetPdgId(2212, 3122);
    MyCat.SetRedMass( (Mass_p*Mass_L)/(Mass_p+Mass_L) );

    pLamSourcePars[3] = RadLarge;
    MyCat.KillTheCat();

    TGraph gUsm_G12;
    gUsm_G12.Set(NumPlotPts);
    gUsm_G12.SetName("gUsm_G12");
    gUsm_G12.SetMarkerColor(DlmCol.GetColor(1));
    gUsm_G12.SetMarkerStyle(20);
    gUsm_G12.SetMarkerSize(0);
    gUsm_G12.SetLineColor(DlmCol.GetColor(1));
    gUsm_G12.SetLineWidth(10*ScaleResolution);
    gUsm_G12.SetFillColor(kWhite);
    for(unsigned uBin=0; uBin<NumPlotPts; uBin++){
        gUsm_G12.SetPoint(uBin, MyCat.GetMomentum(uBin), MyCat.GetCorrFun(uBin));
    }

    MyCat.SetAnaSource(0, RadSmall);
    MyCat.KillTheCat();
    TGraph gUsm_G08;
    gUsm_G08.Set(NumPlotPts);
    gUsm_G08.SetName("gUsm_G08");
    gUsm_G08.SetMarkerColor(DlmCol.GetColor(1));
    gUsm_G08.SetMarkerStyle(20);
    gUsm_G08.SetMarkerSize(0);
    gUsm_G08.SetLineColor(DlmCol.GetColor(1));
    gUsm_G08.SetLineWidth(10*ScaleResolution);
    gUsm_G08.SetFillColor(kWhite);
    for(unsigned uBin=0; uBin<NumPlotPts; uBin++){
        gUsm_G08.SetPoint(uBin, MyCat.GetMomentum(uBin), MyCat.GetCorrFun(uBin));
    }

    TGraph gLedni_G12;
    gLedni_G12.Set(NumPlotPts);
    gLedni_G12.SetName("gLedni_G12");
    gLedni_G12.SetMarkerColor(DlmCol.GetColor(2));
    gLedni_G12.SetMarkerStyle(20);
    gLedni_G12.SetMarkerSize(0);
    gLedni_G12.SetLineColor(DlmCol.GetColor(2));
    gLedni_G12.SetLineWidth(6*ScaleResolution);
    gLedni_G12.SetFillColor(kWhite);

    TGraph gRatio_G12;
    gRatio_G12.Set(NumPlotPts);
    gRatio_G12.SetName("gRatio_G12");
    gRatio_G12.SetMarkerColor(kBlack);
    gRatio_G12.SetMarkerStyle(20);
    gRatio_G12.SetMarkerSize(0);
    gRatio_G12.SetLineColor(kBlack);
    gRatio_G12.SetLineWidth(10*ScaleResolution);
    gRatio_G12.SetFillColor(kWhite);

    const double UsmaniScatPar[4] = {2.88,2.92,1.66,3.78};

    for(unsigned uBin=0; uBin<NumPlotPts; uBin++){
        Momentum = MyCat.GetMomentum(uBin);
        gLedni_G12.SetPoint(uBin, Momentum, Lednicky_SingletTriplet(Momentum, &RadLarge, UsmaniScatPar));

        gUsm_G12.GetPoint(uBin, Momentum, CkVal);
        gRatio_G12.SetPoint(uBin, Momentum, Lednicky_SingletTriplet(Momentum, &RadLarge, UsmaniScatPar)/CkVal);
    }

    MyCat.SetAnaSource(0, RadSmall);
    MyCat.KillTheCat();
    TGraph gLedni_G08;
    gLedni_G08.Set(NumPlotPts);
    gLedni_G08.SetName("gLedni_G08");
    gLedni_G08.SetMarkerColor(DlmCol.GetColor(2));
    gLedni_G08.SetMarkerStyle(20);
    gLedni_G08.SetMarkerSize(0);
    gLedni_G08.SetLineColor(DlmCol.GetColor(2));
    gLedni_G08.SetLineWidth(6*ScaleResolution);
    gLedni_G08.SetFillColor(kWhite);

    TGraph gRatio_G08;
    gRatio_G08.Set(NumPlotPts);
    gRatio_G08.SetName("gRatio_G08");
    gRatio_G08.SetMarkerColor(kBlack);
    gRatio_G08.SetMarkerStyle(20);
    gRatio_G08.SetMarkerSize(0);
    gRatio_G08.SetLineColor(kBlack);
    gRatio_G08.SetLineWidth(10*ScaleResolution);
    gRatio_G08.SetFillColor(kWhite);

    for(unsigned uBin=0; uBin<NumPlotPts; uBin++){
        Momentum = MyCat.GetMomentum(uBin);
        gLedni_G08.SetPoint(uBin, Momentum, Lednicky_SingletTriplet(Momentum, &RadSmall, UsmaniScatPar));

        gUsm_G08.GetPoint(uBin, Momentum, CkVal);
        gRatio_G08.SetPoint(uBin, Momentum, Lednicky_SingletTriplet(Momentum, &RadSmall, UsmaniScatPar)/CkVal);
    }

    const double Yoffset = 1;
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

    hAxis->GetYaxis()->SetRangeUser(0.5, 4);
    hAxis->GetYaxis()->SetLimits(0.5, 4);

    hAxis->GetXaxis()->SetNdivisions(505);

    TH1F* hAxisRatio = new TH1F("hAxisRatio", "hAxisRatio", NumPlotPts,Plot_kMin,Plot_kMax);
    hAxisRatio->SetStats(false);
    hAxisRatio->SetTitle("");
    hAxisRatio->GetXaxis()->SetLabelSize(0.065);
    hAxisRatio->GetXaxis()->SetTitle("k (MeV)");
    hAxisRatio->GetXaxis()->CenterTitle();
    hAxisRatio->GetXaxis()->SetTitleOffset(1.25);
    hAxisRatio->GetXaxis()->SetLabelOffset(0.02);
    hAxisRatio->GetXaxis()->SetTitleSize(0.075);
    hAxisRatio->GetYaxis()->SetLabelSize(0.065);
    hAxisRatio->GetYaxis()->SetTitle("Ratio");
    hAxisRatio->GetYaxis()->CenterTitle();
    hAxisRatio->GetYaxis()->SetTitleOffset(Yoffset/1.5);
    hAxisRatio->GetYaxis()->SetTitleSize(0.075);

    hAxisRatio->GetXaxis()->SetNdivisions(505);
    hAxisRatio->GetYaxis()->SetRangeUser(0.88, 1.06);
    hAxisRatio->GetYaxis()->SetLimits(0.88, 1.06);
    hAxisRatio->GetYaxis()->SetNdivisions(603);

    TPaveText* PT12 = new TPaveText(0.5,0.68,0.95,0.95, "blNDC");//lbrt
    PT12->SetName("PT12");
    PT12->SetBorderSize(1);
    PT12->SetTextSize(0.075);
    PT12->SetFillColor(kWhite);
    PT12->SetTextFont(22);
    PT12->AddText("Gauss source");
    PT12->AddText(TString::Format("R = %.2f fm", RadLarge));

    TPaveText* PT08 = new TPaveText(0.5,0.68,0.95,0.95, "blNDC");//lbrt
    PT08->SetName("PT08");
    PT08->SetBorderSize(1);
    PT08->SetTextSize(0.075);
    PT08->SetFillColor(kWhite);
    PT08->SetTextFont(22);
    PT08->AddText("Gauss source");
    PT08->AddText(TString::Format("R = %.2f fm", RadSmall));

    TF1* UnitLine = new TF1("UnitLine", "1",Plot_kMin,Plot_kMax);
    UnitLine->SetLineColor(kGray+1);
    UnitLine->SetLineWidth(6*ScaleResolution);

    DLM_SubPads DrawBoard_Ck_pp(4620,2599);
    DrawBoard_Ck_pp.AddSubPadTL(0,0.5,0.08,0.63);//lrtb
    DrawBoard_Ck_pp.AddSubPadTL(0,0.5,0.63,1);//lrtb
    DrawBoard_Ck_pp.AddSubPadTL(0.5,1,0.08,0.63);//lrtb
    DrawBoard_Ck_pp.AddSubPadTL(0.5,1,0.63,1);//lrtb
    DrawBoard_Ck_pp.AddSubPadTL(0,1,0,0.07);//lrtb

    DrawBoard_Ck_pp.SetMargin(0, 0.09, 0.04, 0, 0.04);//lrbt
    DrawBoard_Ck_pp.SetMargin(1, 0.09, 0.04, 0.14, 0);//lrbt
    DrawBoard_Ck_pp.SetMargin(2, 0.09, 0.04, 0, 0.04);//lrbt
    DrawBoard_Ck_pp.SetMargin(3, 0.09, 0.04, 0.14, 0);//lrbt
    DrawBoard_Ck_pp.SetMargin(4, 0, 0, 0, 0);//lrbt

    DrawBoard_Ck_pp.GetCanvas()->SetName(TString::Format("DrawBoard_Ck_pp"));
    DrawBoard_Ck_pp.cd(0);
    hAxis->Draw("AXIS");
    gUsm_G12.Draw("same,C");
    gLedni_G12.Draw("same,C");
    UnitLine->Draw("same");
    PT12->Draw("same");
    DrawBoard_Ck_pp.SetLabelSize(0, hAxis->GetXaxis(), 18);
    DrawBoard_Ck_pp.SetLabelSize(0, hAxis->GetYaxis(), 18);
    DrawBoard_Ck_pp.SetTitleSize(0, hAxis->GetXaxis(), 18);
    DrawBoard_Ck_pp.SetTitleSize(0, hAxis->GetYaxis(), 18);

    DrawBoard_Ck_pp.cd(1);
    hAxisRatio->Draw("AXIS");
    gRatio_G12.Draw("same,C");
    UnitLine->Draw("same");
    DrawBoard_Ck_pp.SetLabelSize(1, hAxisRatio->GetXaxis(), 18);
    DrawBoard_Ck_pp.SetLabelSize(1, hAxisRatio->GetYaxis(), 18);
    DrawBoard_Ck_pp.SetTitleSize(1, hAxisRatio->GetXaxis(), 18);
    DrawBoard_Ck_pp.SetTitleSize(1, hAxisRatio->GetYaxis(), 18);

    DrawBoard_Ck_pp.cd(2);
    hAxis->Draw("AXIS");
    gUsm_G08.Draw("same,C");
    gLedni_G08.Draw("same,C");
    UnitLine->Draw("same");
    PT08->Draw("same");
    DrawBoard_Ck_pp.SetLabelSize(0, hAxis->GetXaxis(), 18);
    DrawBoard_Ck_pp.SetLabelSize(0, hAxis->GetYaxis(), 18);
    DrawBoard_Ck_pp.SetTitleSize(0, hAxis->GetXaxis(), 18);
    DrawBoard_Ck_pp.SetTitleSize(0, hAxis->GetYaxis(), 18);

    DrawBoard_Ck_pp.cd(3);
    hAxisRatio->Draw("AXIS");
    gRatio_G08.Draw("same,C");
    UnitLine->Draw("same");
    DrawBoard_Ck_pp.SetLabelSize(1, hAxisRatio->GetXaxis(), 18);
    DrawBoard_Ck_pp.SetLabelSize(1, hAxisRatio->GetYaxis(), 18);
    DrawBoard_Ck_pp.SetTitleSize(1, hAxisRatio->GetXaxis(), 18);
    DrawBoard_Ck_pp.SetTitleSize(1, hAxisRatio->GetYaxis(), 18);


    TLegend* myLegend = new TLegend(0.12,0,0.95,0.96);//lbrt
    myLegend->SetName("myLegend");
    myLegend->SetTextSize(0.6);
    myLegend->SetNColumns(3);
    myLegend->AddEntry(&gUsm_G12,   "Usmani            ");
    myLegend->AddEntry(&gLedni_G12,  "Lednicky         ");
    myLegend->AddEntry(&gRatio_G12, "Ratio to Usmani   ");

    DrawBoard_Ck_pp.cd(4);
    myLegend->Draw();

    DrawBoard_Ck_pp.GetCanvas()->SaveAs(TString::Format("%sBormioPlotCk_pL.png",OutputDir.Data()));
    hAxis->SetLineWidth(hAxis->GetLineWidth()/2.5);
    gUsm_G08.SetLineWidth(gUsm_G08.GetLineWidth()/2.5);
    gUsm_G12.SetLineWidth(gUsm_G12.GetLineWidth()/2.5);
    gLedni_G08.SetLineWidth(gLedni_G08.GetLineWidth()/2.5);
    gLedni_G12.SetLineWidth(gLedni_G12.GetLineWidth()/2.5);
    gRatio_G08.SetLineWidth(gRatio_G08.GetLineWidth()/2.5);
    gRatio_G12.SetLineWidth(gRatio_G12.GetLineWidth()/2.5);
    UnitLine->SetLineWidth(UnitLine->GetLineWidth()/2.5);
    DrawBoard_Ck_pp.GetCanvas()->SaveAs(TString::Format("%sBormioPlotCk_pL.pdf",OutputDir.Data()));
    gUsm_G08.SetLineWidth(gUsm_G08.GetLineWidth()*2.5);
    gUsm_G12.SetLineWidth(gUsm_G12.GetLineWidth()*2.5);
    gLedni_G08.SetLineWidth(gLedni_G08.GetLineWidth()*2.5);
    gLedni_G12.SetLineWidth(gLedni_G12.GetLineWidth()*2.5);
    gRatio_G08.SetLineWidth(gRatio_G08.GetLineWidth()*2.5);
    gRatio_G12.SetLineWidth(gRatio_G12.GetLineWidth()*2.5);
    UnitLine->SetLineWidth(UnitLine->GetLineWidth()*2.5);

    delete hAxis;
    delete hAxisRatio;
    delete UnitLine;
    delete myLegend;
    delete PT12;
    delete PT08;

}

void CATSpaper::PlotCk_EPOS(){

    const double Plot_kMin_pp = 0;
    const double Plot_kMax_pp = 240;
    const double NumPlotPts_pp = 60;

    const double Plot_kMin_pL = 5;
    const double Plot_kMax_pL = 245;
    const double NumPlotPts_pL = 24;

    const double RadLarge = 1.2;
    const double RadSmall = 0.85;
    //const double RadLarge = 0.8;
    //const double RadSmall = 0.4;

    const TString SourceType = "Gauss";
    const bool ThetaDepSource = true;
    const bool TauCorrection = false;
    const double EposRenorm = 1.0;

    double CkVal;
    double Momentum;
    DLM_DtColor DlmCol;

    CATSparameters cSourcePars_pp(CATSparameters::tSource,1,true);
    cSourcePars_pp.SetParameters(&ppSourcePars[3]);

    CATS MyCat_pp;
    MyCat_pp.SetThetaDependentSource(ThetaDepSource);
    MyCat_pp.SetTauCorrection(TauCorrection);
    MyCat_pp.SetUseAnalyticSource(true);
    if(MyCat_pp.GetThetaDependentSource()){
        if(SourceType=="Gauss") MyCat_pp.SetAnaSource(GaussSourceTheta, cSourcePars_pp);
        else if(SourceType=="Cauchy") MyCat_pp.SetAnaSource(CauchySourceTheta, cSourcePars_pp);
        else {printf("FUCK!!!\n"); return;}
    }
    else{
        if(SourceType=="Gauss") MyCat_pp.SetAnaSource(GaussSource, cSourcePars_pp);
        else if(SourceType=="Cauchy") MyCat_pp.SetAnaSource(CauchySource, cSourcePars_pp);
        else {printf("FUCK!!!\n"); return;}
    }

    MyCat_pp.SetExcludeFailedBins(false);
    MyCat_pp.SetMomBins(NumPlotPts_pp,Plot_kMin_pp,Plot_kMax_pp);

    MyCat_pp.SetNumChannels(4);
    MyCat_pp.SetNumPW(0,2);
    MyCat_pp.SetNumPW(1,2);
    MyCat_pp.SetNumPW(2,2);
    MyCat_pp.SetNumPW(3,2);
    MyCat_pp.SetSpin(0,0);
    MyCat_pp.SetSpin(1,1);
    MyCat_pp.SetSpin(2,1);
    MyCat_pp.SetSpin(3,1);
    MyCat_pp.SetChannelWeight(0, 3./12.);
    MyCat_pp.SetChannelWeight(1, 1./12.);
    MyCat_pp.SetChannelWeight(2, 3./12.);
    MyCat_pp.SetChannelWeight(3, 5./12.);
    MyCat_pp.SetShortRangePotential(0,0,fDlmPot,ppPars1S0);
    MyCat_pp.SetShortRangePotential(1,1,fDlmPot,ppPars3P0);
    MyCat_pp.SetShortRangePotential(2,1,fDlmPot,ppPars3P1);
    MyCat_pp.SetShortRangePotential(3,1,fDlmPot,ppPars3P2);
    MyCat_pp.SetQ1Q2(1);
    MyCat_pp.SetPdgId(2212, 2212);
    MyCat_pp.SetRedMass( 0.5*Mass_p );

    ppSourcePars[3] = RadLarge;
    MyCat_pp.KillTheCat();

    TGraph g_pp_G12;
    g_pp_G12.Set(NumPlotPts_pp);
    g_pp_G12.SetName("g_pp_G12");
    g_pp_G12.SetMarkerColor(DlmCol.GetColor(0));
    g_pp_G12.SetMarkerStyle(20);
    g_pp_G12.SetMarkerSize(0);
    g_pp_G12.SetLineColor(DlmCol.GetColor(0));
    g_pp_G12.SetLineWidth(10);
    g_pp_G12.SetFillColor(kWhite);
    for(unsigned uBin=0; uBin<NumPlotPts_pp; uBin++){
        g_pp_G12.SetPoint(uBin, MyCat_pp.GetMomentum(uBin), MyCat_pp.GetCorrFun(uBin));
    }

    MyCat_pp.SetAnaSource(0, RadSmall);
    MyCat_pp.KillTheCat();

    TGraph g_pp_G08;
    g_pp_G08.Set(NumPlotPts_pp);
    g_pp_G08.SetName("g_pp_G08");
    g_pp_G08.SetMarkerColor(DlmCol.GetColor(4));
    g_pp_G08.SetMarkerStyle(20);
    g_pp_G08.SetMarkerSize(0);
    g_pp_G08.SetLineColor(DlmCol.GetColor(4));
    g_pp_G08.SetLineWidth(10);
    g_pp_G08.SetFillColor(kWhite);
    for(unsigned uBin=0; uBin<NumPlotPts_pp; uBin++){
        g_pp_G08.SetPoint(uBin, MyCat_pp.GetMomentum(uBin), MyCat_pp.GetCorrFun(uBin));
    }

    MyCat_pp.SetUseAnalyticSource(false);
    MyCat_pp.SetMaxPairsPerBin(16000);
    MyCat_pp.SetMaxPairsToRead(256e6);
    MyCat_pp.SetMixingDepth(8);
    MyCat_pp.SetInputFileName("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/Scratch9_OSCAR1997_AllLrzProtonFiles.f19");
    MyCat_pp.SetTransportRenorm(EposRenorm);
    MyCat_pp.KillTheCat(CATS::kSourceChanged);

    TGraphErrors g_pp_EPOS;
    g_pp_EPOS.Set(NumPlotPts_pp);
    g_pp_EPOS.SetName("g_pp_EPOS");
    g_pp_EPOS.SetMarkerColor(DlmCol.GetColor(6));
    g_pp_EPOS.SetMarkerStyle(20);
    g_pp_EPOS.SetMarkerSize(0);
    g_pp_EPOS.SetLineColor(DlmCol.GetColor(6));
    g_pp_EPOS.SetLineWidth(6);
    g_pp_EPOS.SetFillColor(kWhite);

    TGraph gRatio_pp_G12;
    gRatio_pp_G12.Set(NumPlotPts_pp);
    gRatio_pp_G12.SetName("gRatio_pp_G12");
    gRatio_pp_G12.SetMarkerColor(DlmCol.GetColor(0));
    gRatio_pp_G12.SetMarkerStyle(20);
    gRatio_pp_G12.SetMarkerSize(0);
    gRatio_pp_G12.SetLineColor(DlmCol.GetColor(0));
    gRatio_pp_G12.SetLineWidth(10);
    gRatio_pp_G12.SetFillColor(kWhite);

    TGraph gRatio_pp_G08;
    gRatio_pp_G08.Set(NumPlotPts_pp);
    gRatio_pp_G08.SetName("gRatio_pp_G08");
    gRatio_pp_G08.SetMarkerColor(DlmCol.GetColor(4));
    gRatio_pp_G08.SetMarkerStyle(20);
    gRatio_pp_G08.SetMarkerSize(0);
    gRatio_pp_G08.SetLineColor(DlmCol.GetColor(4));
    gRatio_pp_G08.SetLineWidth(10);
    gRatio_pp_G08.SetFillColor(kWhite);

    for(unsigned uBin=0; uBin<NumPlotPts_pp; uBin++){
        g_pp_EPOS.SetPoint(uBin, MyCat_pp.GetMomentum(uBin), MyCat_pp.GetCorrFun(uBin));
        //g_pp_EPOS.SetPointError(uBin, 0, MyCat_pp.GetCorrFunErr(uBin));

        g_pp_G12.GetPoint(uBin, Momentum, CkVal);
        gRatio_pp_G12.SetPoint(uBin, Momentum, CkVal/MyCat_pp.GetCorrFun(uBin));

        g_pp_G08.GetPoint(uBin, Momentum, CkVal);
        gRatio_pp_G08.SetPoint(uBin, Momentum, CkVal/MyCat_pp.GetCorrFun(uBin));
    }

    CATSparameters cSourcePars_pL(CATSparameters::tSource,1,true);
    cSourcePars_pL.SetParameters(&pLamSourcePars[3]);

    CATS MyCat_pL;
    MyCat_pL.SetThetaDependentSource(ThetaDepSource);
    MyCat_pL.SetTauCorrection(TauCorrection);
    MyCat_pL.SetUseAnalyticSource(true);
    if(MyCat_pL.GetThetaDependentSource()){
        if(SourceType=="Gauss") MyCat_pL.SetAnaSource(GaussSourceTheta, cSourcePars_pL);
        else if(SourceType=="Cauchy") MyCat_pL.SetAnaSource(CauchySourceTheta, cSourcePars_pL);
        else {printf("FUCK!!!\n"); return;}
    }
    else{
        if(SourceType=="Gauss") MyCat_pL.SetAnaSource(GaussSource, cSourcePars_pL);
        else if(SourceType=="Cauchy") MyCat_pL.SetAnaSource(CauchySource, cSourcePars_pL);
        else {printf("FUCK!!!\n"); return;}
    }

    MyCat_pL.SetExcludeFailedBins(false);
    MyCat_pL.SetMomBins(NumPlotPts_pL,Plot_kMin_pL,Plot_kMax_pL);

    MyCat_pL.SetNumChannels(4);
    MyCat_pL.SetNumPW(0,2);
    MyCat_pL.SetNumPW(1,2);
    MyCat_pL.SetNumPW(2,2);
    MyCat_pL.SetNumPW(3,2);
    MyCat_pL.SetSpin(0,0);
    MyCat_pL.SetSpin(1,1);
    MyCat_pL.SetSpin(2,1);
    MyCat_pL.SetSpin(3,1);
    MyCat_pL.SetChannelWeight(0, 1./4.);
    MyCat_pL.SetChannelWeight(1, 1./12.);
    MyCat_pL.SetChannelWeight(2, 3./12.);
    MyCat_pL.SetChannelWeight(3, 5./12.);
    //MyCat_pL.SetShortRangePotential(0,0,fDlmPot,pLamPars1S0);
    //MyCat_pL.SetShortRangePotential(1,0,fDlmPot,pLamPars3S1);
    MyCat_pL.SetQ1Q2(0);
    MyCat_pL.SetPdgId(2212, 3122);
    MyCat_pL.SetRedMass( (Mass_p*Mass_L)/(Mass_p+Mass_L) );

    pLamSourcePars[3] = RadLarge;
    complex<double>**** WaveFunctionU; double*** PhaseShifts; double* RadBins; unsigned NumRadBins=0;
    InitHaidenbauerNLO("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/Haidenbauer/pLambdaNLO/",MyCat_pL,&WaveFunctionU,&PhaseShifts,&RadBins,NumRadBins);
    for(unsigned uBin=0; uBin<NumPlotPts_pL; uBin++){
        MyCat_pL.UseExternalWaveFunction(uBin,0,0,WaveFunctionU[uBin][0][0], NumRadBins, RadBins, PhaseShifts[uBin][0][0]);
        MyCat_pL.UseExternalWaveFunction(uBin,0,1,WaveFunctionU[uBin][0][1], NumRadBins, RadBins, PhaseShifts[uBin][0][1]);
        MyCat_pL.UseExternalWaveFunction(uBin,1,0,WaveFunctionU[uBin][1][0], NumRadBins, RadBins, PhaseShifts[uBin][1][0]);
        MyCat_pL.UseExternalWaveFunction(uBin,1,1,WaveFunctionU[uBin][1][1], NumRadBins, RadBins, PhaseShifts[uBin][1][1]);
        MyCat_pL.UseExternalWaveFunction(uBin,2,0,WaveFunctionU[uBin][1][0], NumRadBins, RadBins, PhaseShifts[uBin][1][0]);
        MyCat_pL.UseExternalWaveFunction(uBin,2,1,WaveFunctionU[uBin][2][1], NumRadBins, RadBins, PhaseShifts[uBin][2][1]);
        MyCat_pL.UseExternalWaveFunction(uBin,3,0,WaveFunctionU[uBin][1][0], NumRadBins, RadBins, PhaseShifts[uBin][1][0]);
        MyCat_pL.UseExternalWaveFunction(uBin,3,1,WaveFunctionU[uBin][3][1], NumRadBins, RadBins, PhaseShifts[uBin][3][1]);
    }
    MyCat_pL.KillTheCat();


    TGraph g_pL_G12;
    g_pL_G12.Set(NumPlotPts_pL);
    g_pL_G12.SetName("g_pL_G12");
    g_pL_G12.SetMarkerColor(DlmCol.GetColor(0));
    g_pL_G12.SetMarkerStyle(20);
    g_pL_G12.SetMarkerSize(0);
    g_pL_G12.SetLineColor(DlmCol.GetColor(0));
    g_pL_G12.SetLineWidth(10);
    g_pL_G12.SetFillColor(kWhite);
    for(unsigned uBin=0; uBin<NumPlotPts_pL; uBin++){
        g_pL_G12.SetPoint(uBin, MyCat_pL.GetMomentum(uBin), MyCat_pL.GetCorrFun(uBin));
    }

    MyCat_pL.SetAnaSource(0, RadSmall);
    MyCat_pL.KillTheCat();

    TGraph g_pL_G08;
    g_pL_G08.Set(NumPlotPts_pL);
    g_pL_G08.SetName("g_pL_G08");
    g_pL_G08.SetMarkerColor(DlmCol.GetColor(4));
    g_pL_G08.SetMarkerStyle(20);
    g_pL_G08.SetMarkerSize(0);
    g_pL_G08.SetLineColor(DlmCol.GetColor(4));
    g_pL_G08.SetLineWidth(10);
    g_pL_G08.SetFillColor(kWhite);
    for(unsigned uBin=0; uBin<NumPlotPts_pL; uBin++){
        g_pL_G08.SetPoint(uBin, MyCat_pL.GetMomentum(uBin), MyCat_pL.GetCorrFun(uBin));
    }

    MyCat_pL.SetUseAnalyticSource(false);
    MyCat_pL.SetMaxPairsPerBin(16000);
    MyCat_pL.SetMaxPairsToRead(256e6);
    MyCat_pL.SetMixingDepth(8);
    MyCat_pL.SetInputFileName("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/Scratch9_OSCAR1997_AllLrzLambdaFiles.f19");
    MyCat_pL.SetTransportRenorm(EposRenorm);
    MyCat_pL.KillTheCat(CATS::kSourceChanged);

    TGraphErrors g_pL_EPOS;
    g_pL_EPOS.Set(NumPlotPts_pL);
    g_pL_EPOS.SetName("g_pL_EPOS");
    g_pL_EPOS.SetMarkerColor(DlmCol.GetColor(6));
    g_pL_EPOS.SetMarkerStyle(20);
    g_pL_EPOS.SetMarkerSize(0);
    g_pL_EPOS.SetLineColor(DlmCol.GetColor(6));
    g_pL_EPOS.SetLineWidth(6);
    g_pL_EPOS.SetFillColor(kWhite);

    TGraph gRatio_pL_G12;
    gRatio_pL_G12.Set(NumPlotPts_pL);
    gRatio_pL_G12.SetName("gRatio_pL_G12");
    gRatio_pL_G12.SetMarkerColor(DlmCol.GetColor(0));
    gRatio_pL_G12.SetMarkerStyle(20);
    gRatio_pL_G12.SetMarkerSize(0);
    gRatio_pL_G12.SetLineColor(DlmCol.GetColor(0));
    gRatio_pL_G12.SetLineWidth(10);
    gRatio_pL_G12.SetFillColor(kWhite);

    TGraph gRatio_pL_G08;
    gRatio_pL_G08.Set(NumPlotPts_pL);
    gRatio_pL_G08.SetName("gRatio_pL_G08");
    gRatio_pL_G08.SetMarkerColor(DlmCol.GetColor(4));
    gRatio_pL_G08.SetMarkerStyle(20);
    gRatio_pL_G08.SetMarkerSize(0);
    gRatio_pL_G08.SetLineColor(DlmCol.GetColor(4));
    gRatio_pL_G08.SetLineWidth(10);
    gRatio_pL_G08.SetFillColor(kWhite);

    for(unsigned uBin=0; uBin<NumPlotPts_pL; uBin++){
        g_pL_EPOS.SetPoint(uBin, MyCat_pL.GetMomentum(uBin), MyCat_pL.GetCorrFun(uBin));
        //g_pL_EPOS.SetPointError(uBin, 0, MyCat_pL.GetCorrFunErr(uBin));

        g_pL_G12.GetPoint(uBin, Momentum, CkVal);
        gRatio_pL_G12.SetPoint(uBin, Momentum, CkVal/MyCat_pL.GetCorrFun(uBin));

        g_pL_G08.GetPoint(uBin, Momentum, CkVal);
        gRatio_pL_G08.SetPoint(uBin, Momentum, CkVal/MyCat_pL.GetCorrFun(uBin));
    }

    char* fBuffer = new char [512];

    FILE * ppFile;
    ppFile = fopen (TString::Format("%sPlotCk_EPOS_pp.txt",OutputDir.Data()),"w");
    fprintf (ppFile, "%18s","k (MeV)");

    sprintf(fBuffer,"R%.2f_%s",RadLarge,"AV18sp");
    fprintf (ppFile,"%18s",fBuffer);

    sprintf(fBuffer,"R%.2f_%s",RadSmall,"AV18sp");
    fprintf (ppFile,"%18s",fBuffer);

    sprintf(fBuffer,"EPOS_%s","AV18sp");
    fprintf (ppFile,"%18s",fBuffer);

    sprintf(fBuffer,"RatioR%.2f",RadLarge);
    fprintf (ppFile,"%18s",fBuffer);

    sprintf(fBuffer,"RatioR%.2f",RadSmall);
    fprintf (ppFile,"%18s",fBuffer);

    double xVal,yVal;
    for(unsigned uMomBin=0; uMomBin<NumPlotPts_pp; uMomBin++){
        fprintf (ppFile, "\n%18.3f",MyCat_pp.GetMomentum(uMomBin));

        g_pp_G12.GetPoint(uMomBin,xVal,yVal);
        if(xVal!=MyCat_pp.GetMomentum(uMomBin)) {printf("WTF MAN pp!\n");}
        fprintf (ppFile, "%18.3f",yVal);

        g_pp_G08.GetPoint(uMomBin,xVal,yVal);
        if(xVal!=MyCat_pp.GetMomentum(uMomBin)) {printf("WTF MAN pp!\n");}
        fprintf (ppFile, "%18.3f",yVal);

        g_pp_EPOS.GetPoint(uMomBin,xVal,yVal);
        if(xVal!=MyCat_pp.GetMomentum(uMomBin)) {printf("WTF MAN pp!\n");}
        fprintf (ppFile, "%18.3f",yVal);

        gRatio_pp_G12.GetPoint(uMomBin,xVal,yVal);
        if(xVal!=MyCat_pp.GetMomentum(uMomBin)) {printf("WTF MAN pp!\n");}
        fprintf (ppFile, "%18.3f",yVal);

        gRatio_pp_G08.GetPoint(uMomBin,xVal,yVal);
        if(xVal!=MyCat_pp.GetMomentum(uMomBin)) {printf("WTF MAN pp!\n");}
        fprintf (ppFile, "%18.3f",yVal);
    }


    FILE * pLamFile;
    pLamFile = fopen (TString::Format("%sPlotCk_EPOS_pL.txt",OutputDir.Data()),"w");
    fprintf (pLamFile, "%18s","k (MeV)");

    sprintf(fBuffer,"R%.2f_%s",RadLarge,"AV18sp");
    fprintf (pLamFile,"%18s",fBuffer);

    sprintf(fBuffer,"R%.2f_%s",RadSmall,"AV18sp");
    fprintf (pLamFile,"%18s",fBuffer);

    sprintf(fBuffer,"EPOS_%s","AV18sp");
    fprintf (pLamFile,"%18s",fBuffer);

    sprintf(fBuffer,"RatioR%.2f",RadLarge);
    fprintf (pLamFile,"%18s",fBuffer);

    sprintf(fBuffer,"RatioR%.2f",RadSmall);
    fprintf (pLamFile,"%18s",fBuffer);

    for(unsigned uMomBin=0; uMomBin<NumPlotPts_pL; uMomBin++){
        fprintf (pLamFile, "\n%18.3f",MyCat_pL.GetMomentum(uMomBin));

        g_pL_G12.GetPoint(uMomBin,xVal,yVal);
        if(xVal!=MyCat_pL.GetMomentum(uMomBin)) {printf("WTF MAN pp!\n");}
        fprintf (pLamFile, "%18.3f",yVal);

        g_pL_G08.GetPoint(uMomBin,xVal,yVal);
        if(xVal!=MyCat_pL.GetMomentum(uMomBin)) {printf("WTF MAN pp!\n");}
        fprintf (pLamFile, "%18.3f",yVal);

        g_pL_EPOS.GetPoint(uMomBin,xVal,yVal);
        if(xVal!=MyCat_pL.GetMomentum(uMomBin)) {printf("WTF MAN pp!\n");}
        fprintf (pLamFile, "%18.3f",yVal);

        gRatio_pL_G12.GetPoint(uMomBin,xVal,yVal);
        if(xVal!=MyCat_pL.GetMomentum(uMomBin)) {printf("WTF MAN pp!\n");}
        fprintf (pLamFile, "%18.3f",yVal);

        gRatio_pL_G08.GetPoint(uMomBin,xVal,yVal);
        if(xVal!=MyCat_pL.GetMomentum(uMomBin)) {printf("WTF MAN pp!\n");}
        fprintf (pLamFile, "%18.3f",yVal);
    }




    TF1* UnitLine = new TF1("UnitLine", "1",0,800);
    UnitLine->SetLineColor(kGray+1);
    UnitLine->SetLineWidth(6);

    const double Yoffset_pp = 1;
    TH1F* hAxis_pp = new TH1F("hAxis_pp", "hAxis_pp", NumPlotPts_pp,Plot_kMin_pp,Plot_kMax_pp);
    hAxis_pp->SetStats(false);
    hAxis_pp->SetTitle("");
    hAxis_pp->GetXaxis()->SetLabelSize(0.065);
    hAxis_pp->GetXaxis()->SetTitle("k (MeV)");
    hAxis_pp->GetXaxis()->CenterTitle();
    hAxis_pp->GetXaxis()->SetTitleOffset(1.35);
    hAxis_pp->GetXaxis()->SetLabelOffset(0.02);
    hAxis_pp->GetXaxis()->SetTitleSize(0.075);
    hAxis_pp->GetYaxis()->SetLabelSize(0.065);
    hAxis_pp->GetYaxis()->SetTitle("C(k)");
    hAxis_pp->GetYaxis()->CenterTitle();
    hAxis_pp->GetYaxis()->SetTitleOffset(Yoffset_pp);
    hAxis_pp->GetYaxis()->SetTitleSize(0.075);

    hAxis_pp->GetYaxis()->SetRangeUser(0, 6.2);
    hAxis_pp->GetYaxis()->SetLimits(0, 6.2);

    hAxis_pp->GetXaxis()->SetNdivisions(504);

    TH1F* hAxisRatio_pp = new TH1F("hAxisRatio_pp", "hAxisRatio_pp", NumPlotPts_pp,Plot_kMin_pp,Plot_kMax_pp);
    hAxisRatio_pp->SetStats(false);
    hAxisRatio_pp->SetTitle("");
    hAxisRatio_pp->GetXaxis()->SetLabelSize(0.065);
    hAxisRatio_pp->GetXaxis()->SetTitle("k (MeV)");
    hAxisRatio_pp->GetXaxis()->CenterTitle();
    hAxisRatio_pp->GetXaxis()->SetTitleOffset(1.25);
    hAxisRatio_pp->GetXaxis()->SetLabelOffset(0.02);
    hAxisRatio_pp->GetXaxis()->SetTitleSize(0.075);
    hAxisRatio_pp->GetYaxis()->SetLabelSize(0.065);
    hAxisRatio_pp->GetYaxis()->SetTitle("Ratio");
    hAxisRatio_pp->GetYaxis()->CenterTitle();
    hAxisRatio_pp->GetYaxis()->SetTitleOffset(Yoffset_pp/1.5);
    hAxisRatio_pp->GetYaxis()->SetTitleSize(0.075);

    hAxisRatio_pp->GetXaxis()->SetNdivisions(504);
    hAxisRatio_pp->GetYaxis()->SetRangeUser(0.65, 1.32);
    hAxisRatio_pp->GetYaxis()->SetLimits(0.65, 1.32);
    hAxisRatio_pp->GetYaxis()->SetNdivisions(604);

    const double Yoffset_pL = 1;
    TH1F* hAxis_pL = new TH1F("hAxis_pL", "hAxis_pL", NumPlotPts_pL,Plot_kMin_pL,Plot_kMax_pL);
    hAxis_pL->SetStats(false);
    hAxis_pL->SetTitle("");
    hAxis_pL->GetXaxis()->SetLabelSize(0.065);
    hAxis_pL->GetXaxis()->SetTitle("k (MeV)");
    hAxis_pL->GetXaxis()->CenterTitle();
    hAxis_pL->GetXaxis()->SetTitleOffset(1.35);
    hAxis_pL->GetXaxis()->SetLabelOffset(0.02);
    hAxis_pL->GetXaxis()->SetTitleSize(0.075);
    hAxis_pL->GetYaxis()->SetLabelSize(0.065);
    hAxis_pL->GetYaxis()->SetTitle("C(k)");
    hAxis_pL->GetYaxis()->CenterTitle();
    hAxis_pL->GetYaxis()->SetTitleOffset(Yoffset_pL);
    hAxis_pL->GetYaxis()->SetTitleSize(0.075);

    hAxis_pL->GetYaxis()->SetRangeUser(0.5, 4);
    hAxis_pL->GetYaxis()->SetLimits(0.5, 4);

    hAxis_pL->GetXaxis()->SetNdivisions(505);

    TH1F* hAxisRatio_pL = new TH1F("hAxisRatio_pL", "hAxisRatio_pL", NumPlotPts_pL,Plot_kMin_pL,Plot_kMax_pL);
    hAxisRatio_pL->SetStats(false);
    hAxisRatio_pL->SetTitle("");
    hAxisRatio_pL->GetXaxis()->SetLabelSize(0.065);
    hAxisRatio_pL->GetXaxis()->SetTitle("k (MeV)");
    hAxisRatio_pL->GetXaxis()->CenterTitle();
    hAxisRatio_pL->GetXaxis()->SetTitleOffset(1.25);
    hAxisRatio_pL->GetXaxis()->SetLabelOffset(0.02);
    hAxisRatio_pL->GetXaxis()->SetTitleSize(0.075);
    hAxisRatio_pL->GetYaxis()->SetLabelSize(0.065);
    hAxisRatio_pL->GetYaxis()->SetTitle("Ratio");
    hAxisRatio_pL->GetYaxis()->CenterTitle();
    hAxisRatio_pL->GetYaxis()->SetTitleOffset(Yoffset_pL/1.5);
    hAxisRatio_pL->GetYaxis()->SetTitleSize(0.075);

    hAxisRatio_pL->GetXaxis()->SetNdivisions(505);
    hAxisRatio_pL->GetYaxis()->SetRangeUser(0.8, 1.32);
    hAxisRatio_pL->GetYaxis()->SetLimits(0.8, 1.32);
    hAxisRatio_pL->GetYaxis()->SetNdivisions(604);

    TPaveText* PTpp = new TPaveText(0.5,0.68,0.95,0.95, "blNDC");//lbrt
    PTpp->SetName("PTpp");
    PTpp->SetBorderSize(1);
    PTpp->SetTextSize(0.075);
    PTpp->SetFillColor(kWhite);
    PTpp->SetTextFont(22);
    PTpp->AddText("pp correlation");
    PTpp->AddText("Argonne V18");

    TPaveText* PTpL = new TPaveText(0.5,0.68,0.95,0.95, "blNDC");//lbrt
    PTpL->SetName("PTpL");
    PTpL->SetBorderSize(1);
    PTpL->SetTextSize(0.075);
    PTpL->SetFillColor(kWhite);
    PTpL->SetTextFont(22);
    PTpL->AddText("p#Lambda correlation");
    PTpL->AddText("NLO");

    DLM_SubPads DrawBoard_Ck_EPOS(1920,1080);
    DrawBoard_Ck_EPOS.AddSubPadTL(0,0.5,0.08,0.63);//lrtb
    DrawBoard_Ck_EPOS.AddSubPadTL(0,0.5,0.63,1);//lrtb
    DrawBoard_Ck_EPOS.AddSubPadTL(0.5,1,0.08,0.63);//lrtb
    DrawBoard_Ck_EPOS.AddSubPadTL(0.5,1,0.63,1);//lrtb
    DrawBoard_Ck_EPOS.AddSubPadTL(0,1,0,0.07);//lrtb

    DrawBoard_Ck_EPOS.SetMargin(0, 0.09, 0.04, 0, 0.04);//lrbt
    DrawBoard_Ck_EPOS.SetMargin(1, 0.09, 0.04, 0.14, 0);//lrbt
    DrawBoard_Ck_EPOS.SetMargin(2, 0.09, 0.04, 0, 0.04);//lrbt
    DrawBoard_Ck_EPOS.SetMargin(3, 0.09, 0.04, 0.14, 0);//lrbt
    DrawBoard_Ck_EPOS.SetMargin(4, 0, 0, 0, 0);//lrbt

    DrawBoard_Ck_EPOS.GetCanvas()->SetName(TString::Format("DrawBoard_Ck_EPOS"));
    DrawBoard_Ck_EPOS.cd(0);
    hAxis_pp->Draw("AXIS");
    g_pp_G12.Draw("same,C");
    g_pp_G08.Draw("same,C");
    g_pp_EPOS.Draw("same,C");
    UnitLine->Draw("same");
    PTpp->Draw("same");
    DrawBoard_Ck_EPOS.SetLabelSize(0, hAxis_pp->GetXaxis(), 18);
    DrawBoard_Ck_EPOS.SetLabelSize(0, hAxis_pp->GetYaxis(), 18);
    DrawBoard_Ck_EPOS.SetTitleSize(0, hAxis_pp->GetXaxis(), 18);
    DrawBoard_Ck_EPOS.SetTitleSize(0, hAxis_pp->GetYaxis(), 18);

    DrawBoard_Ck_EPOS.cd(1);
    hAxisRatio_pp->Draw("AXIS");
    gRatio_pp_G12.Draw("same,C");
    gRatio_pp_G08.Draw("same,C");
    UnitLine->Draw("same");
    DrawBoard_Ck_EPOS.SetLabelSize(1, hAxisRatio_pp->GetXaxis(), 18);
    DrawBoard_Ck_EPOS.SetLabelSize(1, hAxisRatio_pp->GetYaxis(), 18);
    DrawBoard_Ck_EPOS.SetTitleSize(1, hAxisRatio_pp->GetXaxis(), 18);
    DrawBoard_Ck_EPOS.SetTitleSize(1, hAxisRatio_pp->GetYaxis(), 18);

    DrawBoard_Ck_EPOS.cd(2);
    hAxis_pL->Draw("AXIS");
    g_pL_G12.Draw("same,C");
    g_pL_G08.Draw("same,C");
    g_pL_EPOS.Draw("same,C");
    UnitLine->Draw("same");
    PTpL->Draw("same");
    DrawBoard_Ck_EPOS.SetLabelSize(0, hAxis_pL->GetXaxis(), 18);
    DrawBoard_Ck_EPOS.SetLabelSize(0, hAxis_pL->GetYaxis(), 18);
    DrawBoard_Ck_EPOS.SetTitleSize(0, hAxis_pL->GetXaxis(), 18);
    DrawBoard_Ck_EPOS.SetTitleSize(0, hAxis_pL->GetYaxis(), 18);

    DrawBoard_Ck_EPOS.cd(3);
    hAxisRatio_pL->Draw("AXIS");
    gRatio_pL_G12.Draw("same,C");
    gRatio_pL_G08.Draw("same,C");
    UnitLine->Draw("same");
    DrawBoard_Ck_EPOS.SetLabelSize(1, hAxisRatio_pL->GetXaxis(), 18);
    DrawBoard_Ck_EPOS.SetLabelSize(1, hAxisRatio_pL->GetYaxis(), 18);
    DrawBoard_Ck_EPOS.SetTitleSize(1, hAxisRatio_pL->GetXaxis(), 18);
    DrawBoard_Ck_EPOS.SetTitleSize(1, hAxisRatio_pL->GetYaxis(), 18);

    TLegend* myLegend = new TLegend(0.12,0,0.95,0.96);//lbrt
    myLegend->SetName("myLegend");
    myLegend->SetTextSize(0.6);
    myLegend->SetNColumns(3);
    myLegend->AddEntry(&g_pp_G12,  TString::Format("%s %.2f fm", SourceType.Data(), RadLarge));
    myLegend->AddEntry(&g_pp_G08,  TString::Format("%s %.2f fm", SourceType.Data(), RadSmall));
    if(EposRenorm==1) myLegend->AddEntry(&g_pp_EPOS,  "EPOS   ");
    else myLegend->AddEntry(&g_pp_EPOS,  TString::Format("EPOS x %.2f", EposRenorm));

    DrawBoard_Ck_EPOS.cd(4);
    myLegend->Draw();

    DrawBoard_Ck_EPOS.GetCanvas()->SaveAs(TString::Format("%sPlotCk_EPOS.png",OutputDir.Data()));
    hAxis_pp->SetLineWidth(hAxis_pp->GetLineWidth()/2.5);
    g_pp_G12.SetLineWidth(g_pp_G12.GetLineWidth()/2.5);
    g_pp_G08.SetLineWidth(g_pp_G08.GetLineWidth()/2.5);
    g_pp_EPOS.SetLineWidth(g_pp_EPOS.GetLineWidth()/2.5);
    gRatio_pp_G12.SetLineWidth(gRatio_pp_G12.GetLineWidth()/2.5);
    gRatio_pp_G08.SetLineWidth(gRatio_pp_G08.GetLineWidth()/2.5);
    hAxis_pL->SetLineWidth(hAxis_pL->GetLineWidth()/2.5);
    g_pL_G12.SetLineWidth(g_pL_G12.GetLineWidth()/2.5);
    g_pL_G08.SetLineWidth(g_pL_G08.GetLineWidth()/2.5);
    g_pL_EPOS.SetLineWidth(g_pL_EPOS.GetLineWidth()/2.5);
    gRatio_pL_G12.SetLineWidth(gRatio_pL_G12.GetLineWidth()/2.5);
    gRatio_pL_G08.SetLineWidth(gRatio_pL_G08.GetLineWidth()/2.5);
    UnitLine->SetLineWidth(UnitLine->GetLineWidth()/2.5);
    DrawBoard_Ck_EPOS.GetCanvas()->SaveAs(TString::Format("%sPlotCk_EPOS.pdf",OutputDir.Data()));
    g_pp_G12.SetLineWidth(g_pp_G12.GetLineWidth()*2.5);
    g_pp_G08.SetLineWidth(g_pp_G08.GetLineWidth()*2.5);
    g_pp_EPOS.SetLineWidth(g_pp_EPOS.GetLineWidth()*2.5);
    gRatio_pp_G12.SetLineWidth(gRatio_pp_G12.GetLineWidth()*2.5);
    gRatio_pp_G08.SetLineWidth(gRatio_pp_G08.GetLineWidth()*2.5);
    g_pL_G12.SetLineWidth(g_pL_G12.GetLineWidth()*2.5);
    g_pL_G08.SetLineWidth(g_pL_G08.GetLineWidth()*2.5);
    g_pL_EPOS.SetLineWidth(g_pL_EPOS.GetLineWidth()*2.5);
    gRatio_pL_G12.SetLineWidth(gRatio_pL_G12.GetLineWidth()*2.5);
    gRatio_pL_G08.SetLineWidth(gRatio_pL_G08.GetLineWidth()*2.5);
    UnitLine->SetLineWidth(UnitLine->GetLineWidth()*2.5);

    delete UnitLine;
    delete hAxis_pp;
    delete hAxisRatio_pp;
    delete hAxis_pL;
    delete hAxisRatio_pL;
    delete myLegend;
    delete PTpp;
    delete PTpL;

    delete [] fBuffer;

    CleanHaidenbauer(MyCat_pL,&WaveFunctionU,&PhaseShifts,&RadBins);

}

void CATSpaper::PlotCk_LL(){

    const double Plot_kMin_LL = 0;
    const double Plot_kMax_LL = 240;
    const double NumPlotPts_LL = 48;

    const unsigned NumPotPts = 128;
    const double rMin = 0.2;
    const double rMax = 1.5;
    const double rStep = (rMax-rMin)/double(NumPotPts-1);

    const double RadLarge = 1.2;
    //const double RadSmall = 0.85;

    const TString SourceType = "Gauss";
    const bool ThetaDepSource = false;

    const unsigned NumPots = 4;
    //const TString PotNames[NumPots] = {"HKMYY","NF42","NF50","ND56"};
    const TString PotNames[NumPots] = {"NF42", "NF50","ND56","fss2"};

    //const double PotScatLen[NumPots] = {-0.575, 3.659, -0.772, -1.179};
    //const double PotEffRan[NumPots] = {6.45,0.975,4.271,4.656};

    double** PotPars = new double* [NumPots];
    for(unsigned uPot=0; uPot<NumPots; uPot++){
        PotPars[uPot] = new double [8];
        for(unsigned uEl=0; uEl<8; uEl++) PotPars[uPot][uEl]=0;

        if(PotNames[uPot]=="HKMYY")
        {PotPars[uPot][2]=-10.96; PotPars[uPot][3]=1.342; PotPars[uPot][4]=-141.75; PotPars[uPot][5]=0.777; PotPars[uPot][6]=2136.6; PotPars[uPot][7]=0.35;}

        if(PotNames[uPot]=="NF42")
        {PotPars[uPot][2]=-878.97; PotPars[uPot][3]=0.6; PotPars[uPot][4]=1048.58; PotPars[uPot][5]=0.45;}

        if(PotNames[uPot]=="NF44")
        {PotPars[uPot][2]=-1066.98; PotPars[uPot][3]=0.6; PotPars[uPot][4]=1646.65; PotPars[uPot][5]=0.45;}

        if(PotNames[uPot]=="NF50")
        {PotPars[uPot][2]=-2007.35; PotPars[uPot][3]=0.6; PotPars[uPot][4]=5678.97; PotPars[uPot][5]=0.45;}

        if(PotNames[uPot]=="ND56")
        {PotPars[uPot][2]=-144.26; PotPars[uPot][3]=1.0; PotPars[uPot][4]=1413.75; PotPars[uPot][5]=0.45;}

        if(PotNames[uPot]=="fss2")
        {PotPars[uPot][2]=-103.9; PotPars[uPot][3]=0.92; PotPars[uPot][4]=658.2; PotPars[uPot][5]=0.41;}

    }
    //PotPars[0][2]=-10.96; PotPars[0][3]=1.342; PotPars[0][4]=-141.75; PotPars[0][5]=0.777; PotPars[0][6]=2136.6; PotPars[0][7]=0.35;
    //PotPars[1][2]=-878.97; PotPars[1][3]=0.6; PotPars[1][4]=1048.58; PotPars[1][5]=0.45;
    //PotPars[2][2]=-2007.35; PotPars[2][3]=0.6; PotPars[2][4]=5678.97; PotPars[2][5]=0.45;
    //PotPars[3][2]=-144.26; PotPars[3][3]=1.0; PotPars[3][4]=1413.75; PotPars[3][5]=0.45;

    double LamLamSourcePars[4] = {0,0,0,RadLarge};

    CATSparameters cSourcePars_LL(CATSparameters::tSource,1,true);
    cSourcePars_LL.SetParameters(&LamLamSourcePars[3]);

    CATS Nala;
    Nala.SetEpsilonConv(1e-7);
    Nala.SetEpsilonProp(1e-7);
    Nala.SetThetaDependentSource(ThetaDepSource);
    Nala.SetUseAnalyticSource(true);
    if(Nala.GetThetaDependentSource()){
        if(SourceType=="Gauss") Nala.SetAnaSource(GaussSourceTheta, cSourcePars_LL);
        else if(SourceType=="Cauchy") Nala.SetAnaSource(CauchySourceTheta, cSourcePars_LL);
        else {printf("FUCK!!!\n"); return;}
    }
    else{
        if(SourceType=="Gauss") Nala.SetAnaSource(GaussSource, cSourcePars_LL);
        else if(SourceType=="Cauchy") Nala.SetAnaSource(CauchySource, cSourcePars_LL);
        else {printf("FUCK!!!\n"); return;}
    }

    Nala.SetExcludeFailedBins(false);
    Nala.SetMomBins(NumPlotPts_LL,Plot_kMin_LL,Plot_kMax_LL);

    Nala.SetNumChannels(2);
    Nala.SetNumPW(0,1);
    Nala.SetNumPW(1,1);
    Nala.SetSpin(0,0);
    Nala.SetSpin(1,1);
    Nala.SetChannelWeight(0, 0.25);
    Nala.SetChannelWeight(1, 0.75);
    Nala.SetQ1Q2(0);
    Nala.SetPdgId(3122, 3122);
    Nala.SetRedMass( 0.5*Mass_L );

    LamLamSourcePars[3] = RadLarge;

    Nala.KillTheCat();
    TGraph gCkBaseline_LL;
    gCkBaseline_LL.Set(NumPlotPts_LL);
    gCkBaseline_LL.SetName(TString::Format("gCkBaseline_LL"));
    gCkBaseline_LL.SetLineColor(kGray+1);
    gCkBaseline_LL.SetLineWidth(6);
    gCkBaseline_LL.SetLineStyle(7);
    gCkBaseline_LL.SetMarkerStyle(1);
    gCkBaseline_LL.SetMarkerSize(0);
    gCkBaseline_LL.SetMarkerColor(kWhite);
    gCkBaseline_LL.SetFillColor(kWhite);
    for(unsigned uBin=0; uBin<NumPlotPts_LL; uBin++){
        gCkBaseline_LL.SetPoint(uBin,Nala.GetMomentum(uBin),Nala.GetCorrFun(uBin));
    }

    TGraph* gCk_LL = new TGraph [NumPots];
    DLM_DtColor DlmCol;
    for(unsigned uPot=0; uPot<NumPots; uPot++){

        FILE * pFile;
        pFile = fopen (TString::Format("%spLam_WF_%s.txt",OutputDir.Data(),PotNames[uPot].Data()).Data(),"w");
        fprintf (pFile, "         k         r      u(r)\n");


        gCk_LL[uPot].Set(NumPlotPts_LL);
        gCk_LL[uPot].SetName(TString::Format("gCk_LL_%s",PotNames[uPot].Data()));
        gCk_LL[uPot].SetLineColor(DlmCol.GetColor(uPot));
        gCk_LL[uPot].SetLineWidth(12.-4.*double(uPot)/double(NumPots-1));
        gCk_LL[uPot].SetMarkerStyle(1);
        gCk_LL[uPot].SetMarkerSize(0);
        gCk_LL[uPot].SetMarkerColor(DlmCol.GetColor(uPot));
        gCk_LL[uPot].SetFillColor(kWhite);
        Nala.SetShortRangePotential(0,0,TripleGaussSum,PotPars[uPot]);
        Nala.KillTheCat();
        //printf("POT %u\n", uPot);
        for(unsigned uBin=0; uBin<NumPlotPts_LL; uBin++){
            //printf(" C(%.2f) = %.3f\n", Nala.GetMomentum(uBin), Nala.GetCorrFun(uBin));
            gCk_LL[uPot].SetPoint(uBin,Nala.GetMomentum(uBin),Nala.GetCorrFun(uBin));

            double WF_rMin=0.1;
            double WF_rMax=16;
            double WF_rStep=0.1;
            for(double WF_r=WF_rMin; WF_r<=WF_rMax; WF_r+=WF_rStep){
                fprintf (pFile, "%.4e %.4e %.4e\n",Nala.GetMomentum(uBin),WF_r,abs(Nala.EvalRadialWaveFunction(uBin,0,0,WF_r,false)));
            }

        }
        fclose(pFile);
    }

    TH1F* hCkAxis = new TH1F("hCkAxis", "hCkAxis", NumPlotPts_LL,Plot_kMin_LL,Plot_kMax_LL);
    hCkAxis->SetStats(false);
    hCkAxis->SetTitle("");
    hCkAxis->GetXaxis()->SetLabelSize(0.065);
    hCkAxis->GetXaxis()->SetTitle("k (MeV)");
    hCkAxis->GetXaxis()->CenterTitle();
    hCkAxis->GetXaxis()->SetTitleOffset(1.15);
    hCkAxis->GetXaxis()->SetLabelOffset(0.02);
    hCkAxis->GetXaxis()->SetTitleSize(0.075);
    hCkAxis->GetYaxis()->SetLabelSize(0.065);
    hCkAxis->GetYaxis()->SetTitle("C(k)");
    hCkAxis->GetYaxis()->CenterTitle();
    hCkAxis->GetYaxis()->SetTitleOffset(0.75);
    hCkAxis->GetYaxis()->SetTitleSize(0.075);
    hCkAxis->GetYaxis()->SetRangeUser(0.4, 1.3);
    hCkAxis->GetYaxis()->SetLimits(0.4, 1.3);
    hCkAxis->GetXaxis()->SetNdivisions(506);

    TLegend* CkLegend = new TLegend(0.7,0.3,0.975,0.58);//lbrt
    CkLegend->SetName("CkLegend");
    CkLegend->SetTextSize(0.05);
    for(unsigned uPot=0; uPot<NumPots; uPot++){
        CkLegend->AddEntry(&gCk_LL[uPot], TString::Format("%s", PotNames[uPot].Data()));
    }
    CkLegend->AddEntry(&gCkBaseline_LL, "W/o strong inter.");

    TF1* UnitLine = new TF1("UnitLine", "1", Plot_kMin_LL, Plot_kMax_LL);
    UnitLine->SetLineColor(kGray);
    UnitLine->SetLineWidth(6);

    TCanvas* cCk_LL = new TCanvas("cCk_LL", "cCk_LL", 1);
    cCk_LL->cd(0); cCk_LL->SetCanvasSize(1920, 1080); cCk_LL->SetMargin(0.12,0.025,0.18,0.05);//lrbt
    hCkAxis->Draw("AXIS");
    gCkBaseline_LL.Draw("same");
    for(unsigned uPot=0; uPot<NumPots; uPot++){
        gCk_LL[uPot].Draw("C,same");
    }
    //UnitLine->Draw("same");
    CkLegend->Draw("same");
    cCk_LL->SaveAs(TString::Format("%sPlotCk_LL.png",OutputDir.Data()));
    for(unsigned uPot=0; uPot<NumPots; uPot++){
        gCk_LL[uPot].SetLineWidth(gCk_LL[uPot].GetLineWidth()/2.5);
    }
    UnitLine->SetLineWidth(UnitLine->GetLineWidth()/2.5);
    cCk_LL->SaveAs(TString::Format("%sPlotCk_LL.pdf",OutputDir.Data()));
    for(unsigned uPot=0; uPot<NumPots; uPot++){
        gCk_LL[uPot].SetLineWidth(gCk_LL[uPot].GetLineWidth()*2.5);
    }
    UnitLine->SetLineWidth(UnitLine->GetLineWidth()*2.5);


    //PLOT THE POTENTIALS

    TGraph* gPot_LL = new TGraph [NumPots];
    for(unsigned uPot=0; uPot<NumPots; uPot++){
        gPot_LL[uPot].Set(NumPotPts);
        gPot_LL[uPot].SetName(TString::Format("gPot_LL_%s",PotNames[uPot].Data()));
        gPot_LL[uPot].SetLineColor(DlmCol.GetColor(uPot));
        gPot_LL[uPot].SetLineWidth(12.-4.*double(uPot)/double(NumPots-1));
        gPot_LL[uPot].SetMarkerStyle(1);
        gPot_LL[uPot].SetMarkerSize(0);
        gPot_LL[uPot].SetMarkerColor(DlmCol.GetColor(uPot));
        gPot_LL[uPot].SetFillColor(kWhite);
        double Radius;
        for(unsigned uRad=0; uRad<NumPotPts; uRad++){
            Radius = rMin + double(uRad)*rStep;
            PotPars[uPot][0] = Radius;
            gPot_LL[uPot].SetPoint(uRad,Radius,TripleGaussSum(PotPars[uPot]));
        }
    }

    TH1F* hPotAxis = new TH1F("hPotAxis", "hPotAxis", NumPotPts,rMin,rMax);
    hPotAxis->SetStats(false);
    hPotAxis->SetTitle("");
    hPotAxis->GetXaxis()->SetLabelSize(0.065);
    hPotAxis->GetXaxis()->SetTitle("r (fm)");
    hPotAxis->GetXaxis()->CenterTitle();
    hPotAxis->GetXaxis()->SetTitleOffset(1.15);
    hPotAxis->GetXaxis()->SetLabelOffset(0.02);
    hPotAxis->GetXaxis()->SetTitleSize(0.075);
    hPotAxis->GetYaxis()->SetLabelSize(0.065);
    hPotAxis->GetYaxis()->SetTitle("V (MeV)");
    hPotAxis->GetYaxis()->CenterTitle();
    hPotAxis->GetYaxis()->SetTitleOffset(0.75);
    hPotAxis->GetYaxis()->SetTitleSize(0.075);
    hPotAxis->GetYaxis()->SetRangeUser(-200, 300);
    hPotAxis->GetYaxis()->SetLimits(-200, 300);
    hPotAxis->GetYaxis()->SetNdivisions(405,false);

    TLegend* PotLegend = new TLegend(0.8,0.6,0.975,0.88);//lbrt
    PotLegend->SetName("PotLegend");
    PotLegend->SetTextSize(0.05);
    for(unsigned uPot=0; uPot<NumPots; uPot++){
        PotLegend->AddEntry(&gCk_LL[uPot], TString::Format("%s", PotNames[uPot].Data()));
    }

    TF1* ZeroLine = new TF1("ZeroLine", "0", rMin, rMax);
    ZeroLine->SetLineColor(kGray);
    ZeroLine->SetLineWidth(6);

    TCanvas* cPot_LL = new TCanvas("cPot_LL", "cPot_LL", 1);
    cPot_LL->cd(0); cPot_LL->SetCanvasSize(1920, 1080); cPot_LL->SetMargin(0.12,0.025,0.18,0.05);//lrbt
    hPotAxis->Draw("AXIS");
    for(unsigned uPot=0; uPot<NumPots; uPot++){
        gPot_LL[uPot].Draw("C,same");
    }
    ZeroLine->Draw("same");
    PotLegend->Draw("same");
    cPot_LL->SaveAs(TString::Format("%sPlotPot_LL.png",OutputDir.Data()));
    for(unsigned uPot=0; uPot<NumPots; uPot++){
        gPot_LL[uPot].SetLineWidth(gPot_LL[uPot].GetLineWidth()/2.5);
    }
    gCkBaseline_LL.SetLineWidth(gCkBaseline_LL.GetLineWidth()/2.5);
    ZeroLine->SetLineWidth(ZeroLine->GetLineWidth()/2.5);
    cPot_LL->SaveAs(TString::Format("%sPlotPot_LL.pdf",OutputDir.Data()));
    for(unsigned uPot=0; uPot<NumPots; uPot++){
        gPot_LL[uPot].SetLineWidth(gPot_LL[uPot].GetLineWidth()*2.5);
    }
    gCkBaseline_LL.SetLineWidth(gCkBaseline_LL.GetLineWidth()*2.5);
    ZeroLine->SetLineWidth(ZeroLine->GetLineWidth()*2.5);


    for(unsigned uPot=0; uPot<NumPots; uPot++){
        delete [] PotPars[uPot];
    }
    delete [] PotPars;
    delete [] gCk_LL;
    delete [] gPot_LL;
    delete hCkAxis;
    delete hPotAxis;
    delete UnitLine;
    delete ZeroLine;
    delete CkLegend;
    delete PotLegend;
    delete cCk_LL;
    delete cPot_LL;

}

void CATSpaper::PlotCk_LL_ver2(){
    const double Plot_kMin = 5;
    const double Plot_kMax = 245;
    const double NumPlotPts = 24;

    const unsigned NumRad = 3;
    const double Radius[NumRad] = {2.0,1.2,0.85};

    const double NumPotentialPts = 128;
    const double rMin = 0.2;
    const double rMax = 1.5;
    const double rStep = (rMax-rMin)/double(NumPotentialPts-1);

    const double lambdaPar = 0.32;
    const double lambdaPurity = 0.03;
    const double lambdaFeed = 1.-lambdaPar-lambdaPurity;

    const TString SourceType = "Gauss";

//    double CkVal0;
//    double Momentum0;
//    double CkVal;
//    double Momentum;
    DLM_DtColor DlmCol;

    const unsigned NumScenarios = 5;
    TString ScenarioName[NumScenarios];
    ScenarioName[0] = "NF42";
    ScenarioName[1] = "NF50";
    ScenarioName[2] = "ND56";
    ScenarioName[3] = "fss2";
    ScenarioName[4] = "Without strong interaction";

    TString ScenarioShortName[NumScenarios];
    ScenarioShortName[0] = "NF42";
    ScenarioShortName[1] = "NF50";
    ScenarioShortName[2] = "ND56";
    ScenarioShortName[3] = "fss2";
    ScenarioShortName[4] = "NoSI";

    double** PotPars = new double* [NumScenarios];
    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        PotPars[uSce] = new double [8];
        for(unsigned uEl=0; uEl<8; uEl++) PotPars[uSce][uEl]=0;

        if(ScenarioName[uSce]=="HKMYY")
        {PotPars[uSce][2]=-10.96; PotPars[uSce][3]=1.342; PotPars[uSce][4]=-141.75; PotPars[uSce][5]=0.777; PotPars[uSce][6]=2136.6; PotPars[uSce][7]=0.35;}

        if(ScenarioName[uSce]=="NF42")
        {PotPars[uSce][2]=-878.97; PotPars[uSce][3]=0.6; PotPars[uSce][4]=1048.58; PotPars[uSce][5]=0.45;}

        if(ScenarioName[uSce]=="NF44")
        {PotPars[uSce][2]=-1066.98; PotPars[uSce][3]=0.6; PotPars[uSce][4]=1646.65; PotPars[uSce][5]=0.45;}

        if(ScenarioName[uSce]=="NF50")
        {PotPars[uSce][2]=-2007.35; PotPars[uSce][3]=0.6; PotPars[uSce][4]=5678.97; PotPars[uSce][5]=0.45;}

        if(ScenarioName[uSce]=="ND56")
        {PotPars[uSce][2]=-144.26; PotPars[uSce][3]=1.0; PotPars[uSce][4]=1413.75; PotPars[uSce][5]=0.45;}

        if(ScenarioName[uSce]=="fss2")
        {PotPars[uSce][2]=-103.9; PotPars[uSce][3]=0.92; PotPars[uSce][4]=658.2; PotPars[uSce][5]=0.41;}

        if(ScenarioName[uSce]=="Without strong interaction")
        {PotPars[uSce][2]=0; PotPars[uSce][3]=0; PotPars[uSce][4]=0; PotPars[uSce][5]=0; PotPars[uSce][6]=0; PotPars[uSce][7]=0;}
    }

    double LamLamSourcePars[4] = {0,0,0,Radius[0]};

    CATSparameters cSourcePars_LL(CATSparameters::tSource,1,true);
    cSourcePars_LL.SetParameters(&LamLamSourcePars[3]);

    CATS Nala;
    Nala.SetEpsilonConv(1e-7);
    Nala.SetEpsilonProp(1e-7);
    //Nala.SetThetaDependentSource(ThetaDepSource);
    Nala.SetUseAnalyticSource(true);
    if(Nala.GetThetaDependentSource()){
        if(SourceType=="Gauss") Nala.SetAnaSource(GaussSourceTheta, cSourcePars_LL);
        else if(SourceType=="Cauchy") Nala.SetAnaSource(CauchySourceTheta, cSourcePars_LL);
        else {printf("FUCK!!!\n"); return;}
    }
    else{
        if(SourceType=="Gauss") Nala.SetAnaSource(GaussSource, cSourcePars_LL);
        else if(SourceType=="Cauchy") Nala.SetAnaSource(CauchySource, cSourcePars_LL);
        else {printf("FUCK!!!\n"); return;}
    }
    Nala.SetExcludeFailedBins(false);
    Nala.SetMomBins(NumPlotPts,Plot_kMin,Plot_kMax);

    Nala.SetNumChannels(2);
    Nala.SetNumPW(0,1);
    Nala.SetNumPW(1,1);
    Nala.SetSpin(0,0);
    Nala.SetSpin(1,1);
    Nala.SetChannelWeight(0,1./4.);
    Nala.SetChannelWeight(1,3./4.);
    Nala.SetQ1Q2(0);
    Nala.SetPdgId(3122, 3122);
    Nala.SetRedMass( 0.5*Mass_L );

    //in bins of Scenario/Radius
    TGraph** grCk;
    grCk = new TGraph* [NumScenarios];

    const TString SigmaMatrixFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_7TeV/ResolutionMatrices/LHC14j4bcde_MeV.root";
    const double TrMatUnitConv = 1;//1000

    TH2F* hSigma_pL;
    TFile* FileSigma = new TFile(SigmaMatrixFileName, "read");
    FileSigma->cd();
    hSigma_pL = (TH2F*)FileSigma->Get("hSigma_pL_InMeV");

    //1/FractionOfBins th number of original bins of the correction matrix are taken into account
    //originally: 3000 MeV => 1/10 should be good most of the time
    const int FractionOfBins = 6;

    TH2F* hSigma_pL_InFuckingMeV = new TH2F("hSigma_pL_InFuckingMeV", "hSigma_pL_InFuckingMeV",
        hSigma_pL->GetNbinsX()/FractionOfBins, hSigma_pL->GetXaxis()->GetBinLowEdge(1)*TrMatUnitConv,hSigma_pL->GetXaxis()->GetBinUpEdge(hSigma_pL->GetNbinsX()/FractionOfBins)*TrMatUnitConv,
        hSigma_pL->GetNbinsY()/FractionOfBins, hSigma_pL->GetYaxis()->GetBinLowEdge(1)*TrMatUnitConv,hSigma_pL->GetXaxis()->GetBinUpEdge(hSigma_pL->GetNbinsY()/FractionOfBins)*TrMatUnitConv);

    for(int iBinX=1; iBinX<=hSigma_pL->GetNbinsX()/FractionOfBins; iBinX++){
        for(int iBinY=1; iBinY<=hSigma_pL->GetNbinsY()/FractionOfBins; iBinY++){
            hSigma_pL_InFuckingMeV->SetBinContent(iBinX, iBinY, hSigma_pL->GetBinContent(iBinX, iBinY));
        }
    }

    //in bins of Scenario/Radius
    TGraph** grCkCorr;
    grCkCorr = new TGraph* [NumScenarios];

    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        grCk[uSce] = new TGraph[NumRad];
        grCkCorr[uSce] = new TGraph[NumRad];
        Nala.SetShortRangePotential(0,0,TripleGaussSum,PotPars[uSce]);
        Nala.KillTheCat(CATS::kAllChanged);

        for(unsigned uRad=0; uRad<NumRad; uRad++){
            Nala.SetAnaSource(0, Radius[uRad]);
            Nala.KillTheCat();
//Save the uncorrected correlation
            grCk[uSce][uRad].Set(NumPlotPts);
            grCk[uSce][uRad].SetName(TString::Format("grCk_%u_%u",uSce,uRad));
            grCk[uSce][uRad].SetMarkerColor(DlmCol.GetColor(1+uSce));
            grCk[uSce][uRad].SetMarkerStyle(20);
            grCk[uSce][uRad].SetMarkerSize(0);
            grCk[uSce][uRad].SetLineColor(DlmCol.GetColor(1+uSce));
            grCk[uSce][uRad].SetLineWidth(8.-2.*double(uSce)/double(NumScenarios-1));
            if(uSce==NumScenarios-1){
                grCk[uSce][uRad].SetLineWidth(8.);
                grCk[uSce][uRad].SetLineColor(kGray+1);
                grCk[uSce][uRad].SetLineStyle(2);
            }
            grCk[uSce][uRad].SetFillColor(kWhite);

            for(unsigned uMomBin=0; uMomBin<NumPlotPts; uMomBin++){
                grCk[uSce][uRad].SetPoint(uMomBin, Nala.GetMomentum(uMomBin), Nala.GetCorrFun(uMomBin));
            }
//Apply the correlations
            DLM_Ck* Ck_LL = new DLM_Ck(1,0,Nala);
            Ck_LL->Update();
            DLM_CkDecomposition CkDec_LL("LambdaLambda",2,*Ck_LL,hSigma_pL_InFuckingMeV);
            CkDec_LL.AddContribution(0,lambdaFeed,DLM_CkDecomposition::cFeedDown,NULL);
            CkDec_LL.AddContribution(1,lambdaPurity,DLM_CkDecomposition::cFake);
            CkDec_LL.Update();

            grCkCorr[uSce][uRad].Set(NumPlotPts);
            grCkCorr[uSce][uRad].SetName(TString::Format("grCkCorr_%u_%u",uSce,uRad));
            grCkCorr[uSce][uRad].SetMarkerColor(DlmCol.GetColor(1+uSce));
            grCkCorr[uSce][uRad].SetMarkerStyle(20);
            grCkCorr[uSce][uRad].SetMarkerSize(0);
            grCkCorr[uSce][uRad].SetLineColor(DlmCol.GetColor(1+uSce));
            grCkCorr[uSce][uRad].SetLineWidth(8.-2.*double(uSce)/double(NumScenarios-1));
            if(uSce==NumScenarios-1){
                grCkCorr[uSce][uRad].SetLineWidth(8.);
                grCkCorr[uSce][uRad].SetLineColor(kGray+1);
                grCkCorr[uSce][uRad].SetLineStyle(2);
            }
            grCkCorr[uSce][uRad].SetFillColor(kWhite);

            for(unsigned uMomBin=0; uMomBin<NumPlotPts; uMomBin++){
                grCkCorr[uSce][uRad].SetPoint(uMomBin, Nala.GetMomentum(uMomBin), CkDec_LL.EvalCk(Nala.GetMomentum(uMomBin)));
            }

            delete Ck_LL;
        }
    }


    char* fBuffer = new char [512];
    FILE * CkFile;
    CkFile = fopen (TString::Format("%sPlotCk_LL_ver2.txt",OutputDir.Data()),"w");
    fprintf (CkFile, "%18s","k (MeV)");
    for(unsigned uRad=0; uRad<NumRad; uRad++){
        for(unsigned uSce=0; uSce<NumScenarios; uSce++){
            sprintf(fBuffer,"R%.2f_%s",Radius[uRad],ScenarioShortName[uSce].Data());
            fprintf (CkFile,"%18s",fBuffer);
            //fprintf (CkFile,"%18s","<-RatioNLO");
        }
    }

    double xVal,yVal;
    for(unsigned uMomBin=0; uMomBin<NumPlotPts; uMomBin++){
        fprintf (CkFile, "\n%18.3f",Nala.GetMomentum(uMomBin));
        for(unsigned uRad=0; uRad<NumRad; uRad++){
            for(unsigned uSce=0; uSce<NumScenarios; uSce++){
                grCk[uSce][uRad].GetPoint(uMomBin,xVal,yVal);
                if(xVal!=Nala.GetMomentum(uMomBin)) {printf("WTF MAN LL!\n");}
                fprintf (CkFile, "%18.3f",yVal);
                //grCkRatio[uSce][uRad].GetPoint(uMomBin,xVal,yVal);
                //fprintf (CkFile, "%18.3f",yVal);
            }
        }
    }

    FILE * CkFileCorr;
    CkFileCorr = fopen (TString::Format("%sPlotCk_LL_ver2_Corrected1.txt",OutputDir.Data()),"w");
    fprintf (CkFileCorr, "%18s","k (MeV)");
    for(unsigned uRad=0; uRad<NumRad; uRad++){
        for(unsigned uSce=0; uSce<NumScenarios; uSce++){
            sprintf(fBuffer,"R%.2f_%s",Radius[uRad],ScenarioShortName[uSce].Data());
            fprintf (CkFileCorr,"%18s",fBuffer);
        }
    }

    for(unsigned uMomBin=0; uMomBin<NumPlotPts; uMomBin++){
        fprintf (CkFileCorr, "\n%18.3f",Nala.GetMomentum(uMomBin));
        for(unsigned uRad=0; uRad<NumRad; uRad++){
            for(unsigned uSce=0; uSce<NumScenarios; uSce++){
                grCkCorr[uSce][uRad].GetPoint(uMomBin,xVal,yVal);
                if(xVal!=Nala.GetMomentum(uMomBin)) {printf("WTF MAN LL!\n");}
                fprintf (CkFileCorr, "%18.3f",yVal);
            }
        }
    }

    const double Yoffset = 1.25;
    TH1F* hAxis = new TH1F("hAxis", "hAxis", NumPlotPts,Plot_kMin,Plot_kMax);
    hAxis->SetStats(false);
    hAxis->SetTitle("");
    hAxis->GetXaxis()->SetLabelSize(0.065);
    hAxis->GetXaxis()->SetTitle("k (MeV)");
    hAxis->GetXaxis()->CenterTitle();
    hAxis->GetXaxis()->SetTitleOffset(1.05);
    hAxis->GetXaxis()->SetLabelOffset(0.02);
    hAxis->GetXaxis()->SetTitleSize(0.075);
    hAxis->GetYaxis()->SetLabelSize(0.065);
    hAxis->GetYaxis()->SetTitle("C(k)");
    hAxis->GetYaxis()->CenterTitle();
    hAxis->GetYaxis()->SetTitleOffset(Yoffset);
    hAxis->GetYaxis()->SetTitleSize(0.075);

    //hAxis->GetYaxis()->SetRangeUser(0.2, 1.3);
    //hAxis->GetYaxis()->SetLimits(0.2, 1.3);
hAxis->GetYaxis()->SetRangeUser(0.8, 1.45);
hAxis->GetYaxis()->SetLimits(0.8, 1.45);

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
    hAxisZoom->GetYaxis()->SetRangeUser(0.85, 1.15);
    hAxisZoom->GetYaxis()->SetLimits(0.85, 1.15);
    hAxisZoom->GetYaxis()->SetNdivisions(603);

    TPaveText** myPT = new TPaveText* [NumRad];
    for(unsigned uRad=0; uRad<NumRad; uRad++){
        myPT[uRad] = new TPaveText(0.5,0.885,0.95,0.975, "blNDC");//lbrt
        myPT[uRad]->SetName("PT12");
        myPT[uRad]->SetBorderSize(1);
        myPT[uRad]->SetTextSize(0.08/1.4);
        myPT[uRad]->SetFillColor(kWhite);
        myPT[uRad]->SetTextFont(22);
        //myPT[uRad]->AddText("Gauss source");
        myPT[uRad]->AddText(TString::Format("R = %.2f fm", Radius[uRad]));
    }

    TF1* UnitLine = new TF1("UnitLine", "1",Plot_kMin,Plot_kMax);
    UnitLine->SetLineColor(kGray+1);
    UnitLine->SetLineWidth(6);

    DLM_SubPads DrawBoard_Ck_LL(2560,1080);
    double PadLength = 1./double(NumRad);
    double CurrentPos = 0;
    for(unsigned uRad=0; uRad<NumRad; uRad++){
        DrawBoard_Ck_LL.AddSubPadTL(CurrentPos,CurrentPos+PadLength,0.08,1.);//lrtb
        //DrawBoard_Ck_LL.AddSubPadTL(CurrentPos,CurrentPos+PadLength,0.63,1);//lrtb
        CurrentPos+=PadLength;
//printf("uRad=%u; %.2f <--> %.2f\n", uRad, );
    }
    DrawBoard_Ck_LL.AddSubPadTL(0,1,0,0.075);//lrtb

    for(unsigned uRad=0; uRad<NumRad; uRad++){
        DrawBoard_Ck_LL.SetMargin(uRad, 0.07, 0.01, 0.15, 0.04);//lrbt
    }
    DrawBoard_Ck_LL.SetMargin(NumRad, 0, 0, 0, 0);//lrbt

    for(unsigned uRad=0; uRad<NumRad; uRad++){
        DrawBoard_Ck_LL.GetCanvas()->SetName(TString::Format("DrawBoard_Ck_LL_%u",uRad));
        DrawBoard_Ck_LL.cd(uRad);
        hAxis->Draw("AXIS");
        //UnitLine->Draw("same");
        for(unsigned uSce=0; uSce<NumScenarios; uSce++) grCk[uSce][uRad].Draw("same,C");
        for(unsigned uSce=0; uSce<NumScenarios; uSce++) myPT[uRad]->Draw("same");

        DrawBoard_Ck_LL.SetLabelSize(uRad, hAxis->GetXaxis(), 17);
        DrawBoard_Ck_LL.SetLabelSize(uRad, hAxis->GetYaxis(), 17);
        DrawBoard_Ck_LL.SetTitleSize(uRad, hAxis->GetXaxis(), 17);
        DrawBoard_Ck_LL.SetTitleSize(uRad, hAxis->GetYaxis(), 17);
    }

    TLegend* myLegend = new TLegend(0.12,0,0.95,0.96);//lbrt
    myLegend->SetName("myLegend");
    myLegend->SetTextSize(0.65);
    myLegend->SetNColumns(NumScenarios);
    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        myLegend->AddEntry(&grCk[uSce][0],  ScenarioName[uSce]);
    }

    DrawBoard_Ck_LL.cd(NumRad);
    myLegend->Draw();

    gStyle->SetLineWidth(2.5);
    DrawBoard_Ck_LL.GetCanvas()->SaveAs(TString::Format("%sPlotCk_LL_ver2.png",OutputDir.Data()));
    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        for(unsigned uRad=0; uRad<NumRad; uRad++){
            grCk[uSce][uRad].SetLineWidth(grCk[uSce][uRad].GetLineWidth()/2.5);
        }
    }
    hAxis->SetLineWidth(hAxis->GetLineWidth()/2.5);
    UnitLine->SetLineWidth(UnitLine->GetLineWidth()/2.5);
    gStyle->SetLineWidth(1);
    DrawBoard_Ck_LL.GetCanvas()->SaveAs(TString::Format("%sPlotCk_LL_ver2.pdf",OutputDir.Data()));

    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        for(unsigned uRad=0; uRad<NumRad; uRad++){
            grCk[uSce][uRad].SetLineWidth(grCk[uSce][uRad].GetLineWidth()*2.5);
        }
    }
    hAxis->SetLineWidth(hAxis->GetLineWidth()*2.5);
    UnitLine->SetLineWidth(UnitLine->GetLineWidth()*2.5);


    DLM_SubPads DrawBoard_Ck_LL_Corrected1(2560,1080);
    CurrentPos=0;
    for(unsigned uRad=0; uRad<NumRad; uRad++){
        DrawBoard_Ck_LL_Corrected1.AddSubPadTL(CurrentPos,CurrentPos+PadLength,0.08,1.);//lrtb
        CurrentPos+=PadLength;
    }
    DrawBoard_Ck_LL_Corrected1.AddSubPadTL(0,1,0,0.075);//lrtb

    for(unsigned uRad=0; uRad<NumRad; uRad++){
        DrawBoard_Ck_LL_Corrected1.SetMargin(uRad, 0.07, 0.01, 0.15, 0.04);//lrbt
    }
    DrawBoard_Ck_LL_Corrected1.SetMargin(NumRad, 0, 0, 0, 0);//lrbt

    for(unsigned uRad=0; uRad<NumRad; uRad++){
        DrawBoard_Ck_LL_Corrected1.GetCanvas()->SetName(TString::Format("DrawBoard_Ck_LL_%u",uRad));
        DrawBoard_Ck_LL_Corrected1.cd(uRad);
        hAxis->Draw("AXIS");
        //UnitLine->Draw("same");
        for(unsigned uSce=0; uSce<NumScenarios; uSce++) grCkCorr[uSce][uRad].Draw("same,C");
        for(unsigned uSce=0; uSce<NumScenarios; uSce++) myPT[uRad]->Draw("same");

        DrawBoard_Ck_LL_Corrected1.SetLabelSize(uRad, hAxis->GetXaxis(), 17);
        DrawBoard_Ck_LL_Corrected1.SetLabelSize(uRad, hAxis->GetYaxis(), 17);
        DrawBoard_Ck_LL_Corrected1.SetTitleSize(uRad, hAxis->GetXaxis(), 17);
        DrawBoard_Ck_LL_Corrected1.SetTitleSize(uRad, hAxis->GetYaxis(), 17);
    }

    DrawBoard_Ck_LL_Corrected1.cd(NumRad);
    myLegend->Draw();



    TPaveText* ptCorrected0 = new TPaveText(0.4,0.47,0.95,0.60, "blNDC");//lbrt
    ptCorrected0->SetName("ptCorrected0");
    ptCorrected0->SetBorderSize(1);
    ptCorrected0->SetTextSize(0.075/1.4);
    ptCorrected0->SetFillColor(kWhite);
    ptCorrected0->SetTextFont(22);
    ptCorrected0->AddText("Corrections were");
    ptCorrected0->AddText("applied for:");

    TPaveText* ptCorrected = new TPaveText(0.4,0.20,0.95,0.47, "blNDC");//lbrt
    ptCorrected->SetName("ptCorrected");
    ptCorrected->SetBorderSize(1);
    ptCorrected->SetTextSize(0.065/1.4);
    ptCorrected->SetFillColor(kWhite);
    ptCorrected->SetTextFont(22);
    ptCorrected->SetTextAlign(11);
    ptCorrected->AddText("* Momentum resolution");
    ptCorrected->AddText("* Residuals assuming:");
    ptCorrected->AddText(TString::Format("   %.0f%% primary #Lambda#Lambda  pairs",lambdaPar*100));
    ptCorrected->AddText(TString::Format("   %.0f%% flat residuals",(1.-lambdaPar)*100));

    DrawBoard_Ck_LL_Corrected1.cd(NumRad/2);
    ptCorrected0->Draw("same");
    ptCorrected->Draw("same");



    gStyle->SetLineWidth(2.5);
    DrawBoard_Ck_LL_Corrected1.GetCanvas()->SaveAs(TString::Format("%sPlotCk_LL_ver2_Corrected1.png",OutputDir.Data()));
    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        for(unsigned uRad=0; uRad<NumRad; uRad++){
            grCk[uSce][uRad].SetLineWidth(grCk[uSce][uRad].GetLineWidth()/2.5);
            grCkCorr[uSce][uRad].SetLineWidth(grCkCorr[uSce][uRad].GetLineWidth()/2.5);
        }
    }
    hAxis->SetLineWidth(hAxis->GetLineWidth()/2.5);
    UnitLine->SetLineWidth(UnitLine->GetLineWidth()/2.5);
    gStyle->SetLineWidth(1);
    DrawBoard_Ck_LL_Corrected1.GetCanvas()->SaveAs(TString::Format("%sPlotCk_LL_ver2_Corrected1.pdf",OutputDir.Data()));

    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        for(unsigned uRad=0; uRad<NumRad; uRad++){
            grCk[uSce][uRad].SetLineWidth(grCk[uSce][uRad].GetLineWidth()*2.5);
            grCkCorr[uSce][uRad].SetLineWidth(grCkCorr[uSce][uRad].GetLineWidth()*2.5);
        }
    }
    hAxis->SetLineWidth(hAxis->GetLineWidth()*2.5);
    UnitLine->SetLineWidth(UnitLine->GetLineWidth()*2.5);


    //!PLOT THE POTENTIALS
    double RAD;
    TGraph* gPot_LL = new TGraph [NumScenarios];
    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        gPot_LL[uSce].Set(NumPotentialPts);
        gPot_LL[uSce].SetName(TString::Format("gPot_LL_%s",ScenarioName[uSce].Data()));
        gPot_LL[uSce].SetLineColor(DlmCol.GetColor(uSce+1));
        gPot_LL[uSce].SetLineWidth(12.-4.*double(uSce)/double(NumScenarios-1));
        gPot_LL[uSce].SetMarkerStyle(1);
        gPot_LL[uSce].SetMarkerSize(0);
        gPot_LL[uSce].SetMarkerColor(DlmCol.GetColor(uSce+1));
        gPot_LL[uSce].SetFillColor(kWhite);
        if(uSce==NumScenarios-1){
            gPot_LL[uSce].SetLineWidth(8.);
            gPot_LL[uSce].SetLineColor(kGray+1);
            gPot_LL[uSce].SetLineStyle(2);
        }

        for(unsigned uRad=0; uRad<NumPotentialPts; uRad++){
            RAD = rMin + double(uRad)*rStep;
            PotPars[uSce][0] = RAD;
            gPot_LL[uSce].SetPoint(uRad,RAD,TripleGaussSum(PotPars[uSce]));
        }
    }

    FILE * PotFile;
    PotFile = fopen (TString::Format("%sPlotPot_LL_ver2.txt",OutputDir.Data()),"w");
    fprintf (PotFile, "%18s","r (fm)");
    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        sprintf(fBuffer,"%s",ScenarioShortName[uSce].Data());
        fprintf (PotFile,"%18s",fBuffer);
    }

    for(unsigned uRad=0; uRad<NumPotentialPts; uRad++){
        RAD = rMin + double(uRad)*rStep;
        fprintf (PotFile, "\n%18.3f",RAD);
        for(unsigned uSce=0; uSce<NumScenarios; uSce++){
            gPot_LL[uSce].GetPoint(uRad,xVal,yVal);
            if(xVal!=RAD) {printf("WTF MAN LL pot!\n");}
            fprintf (PotFile, "%18.3f",yVal);
        }
    }


    TH1F* hPotAxis = new TH1F("hPotAxis", "hPotAxis", NumPotentialPts,rMin,rMax);
    hPotAxis->SetStats(false);
    hPotAxis->SetTitle("");
    hPotAxis->GetXaxis()->SetLabelSize(0.065);
    hPotAxis->GetXaxis()->SetTitle("r (fm)");
    hPotAxis->GetXaxis()->CenterTitle();
    hPotAxis->GetXaxis()->SetTitleOffset(1.15);
    hPotAxis->GetXaxis()->SetLabelOffset(0.02);
    hPotAxis->GetXaxis()->SetTitleSize(0.075);
    hPotAxis->GetYaxis()->SetLabelSize(0.065);
    hPotAxis->GetYaxis()->SetTitle("V (MeV)");
    hPotAxis->GetYaxis()->CenterTitle();
    hPotAxis->GetYaxis()->SetTitleOffset(0.75);
    hPotAxis->GetYaxis()->SetTitleSize(0.075);
    hPotAxis->GetYaxis()->SetRangeUser(-200, 300);
    hPotAxis->GetYaxis()->SetLimits(-200, 300);
    hPotAxis->GetYaxis()->SetNdivisions(405,false);

    TLegend* PotLegend = new TLegend(0.8,0.6,0.975,0.88);//lbrt
    PotLegend->SetName("PotLegend");
    PotLegend->SetTextSize(0.05);
    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        if(uSce!=NumScenarios-1) PotLegend->AddEntry(&gPot_LL[uSce], TString::Format("%s", ScenarioName[uSce].Data()));
    }


    //TF1* ZeroLine = new TF1("ZeroLine", "0", rMin, rMax);
    //ZeroLine->SetLineColor(kGray);
    //ZeroLine->SetLineWidth(6);

    TCanvas* cPot_LL = new TCanvas("cPot_LL", "cPot_LL", 1);
    cPot_LL->cd(0); cPot_LL->SetCanvasSize(1920, 1080); cPot_LL->SetMargin(0.12,0.025,0.18,0.05);//lrbt
    hPotAxis->Draw("AXIS");
    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        gPot_LL[uSce].Draw("C,same");
    }
    //ZeroLine->Draw("same");
    PotLegend->Draw("same");
    cPot_LL->SaveAs(TString::Format("%sPlotPot_LL_ver2.png",OutputDir.Data()));
    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        gPot_LL[uSce].SetLineWidth(gPot_LL[uSce].GetLineWidth()/2.5);
    }
    //gCkBaseline_LL.SetLineWidth(gCkBaseline_LL.GetLineWidth()/2.5);
    //ZeroLine->SetLineWidth(ZeroLine->GetLineWidth()/2.5);
    cPot_LL->SaveAs(TString::Format("%sPlotPot_LL_ver2.pdf",OutputDir.Data()));
    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        gPot_LL[uSce].SetLineWidth(gPot_LL[uSce].GetLineWidth()*2.5);
    }
    //gCkBaseline_LL.SetLineWidth(gCkBaseline_LL.GetLineWidth()*2.5);
    //ZeroLine->SetLineWidth(ZeroLine->GetLineWidth()*2.5);


    delete hAxis;
    delete hAxisZoom;
    delete UnitLine;
    delete myLegend;
    delete ptCorrected0;
    delete ptCorrected;

    delete hSigma_pL_InFuckingMeV;

    for(unsigned uRad=0; uRad<NumRad; uRad++) delete myPT[uRad];
    delete [] myPT;

    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        delete [] grCk[uSce];
        delete [] grCkCorr[uSce];
    }
    delete [] grCk;
    delete [] grCkCorr;
    delete [] gPot_LL;
    delete hPotAxis;
    delete PotLegend;
    //delete ZeroLine;
    delete cPot_LL;

    delete [] fBuffer;

    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        delete [] PotPars[uSce];
    }
    delete [] PotPars;

}



void CATSpaper::LamLamStudy1(){
    const double Plot_kMin = 5;
    const double Plot_kMax = 205;
    const double NumPlotPts = 20;

    //const unsigned NumRad = 3;
    //const double Radius[NumRad] = {0.8,1.15,3.0};
    //const unsigned NumRad = 2;
    //const double Radius[NumRad] = {1.15,3.0};

    const unsigned NumRad = 1;
    const double Radius[NumRad] = {1.19};

    const double NumPotentialPts = 128;
    const double rMin = 0.2;
    const double rMax = 2.4;
    const double rStep = (rMax-rMin)/double(NumPotentialPts-1);

    const double lambdaPar = 0.37;
    const double lambdaPurity = 0.03;
    const double lambdaFeed = 1.-lambdaPar-lambdaPurity;

    const TString SourceType = "Gauss";

//    double CkVal0;
//    double Momentum0;
//    double CkVal;
//    double Momentum;
    DLM_DtColor DlmCol;

/*
    const unsigned NumScenarios = 8;
    TString ScenarioName[NumScenarios];
    ScenarioName[0] = "ND46";
    ScenarioName[1] = "ND48";
    ScenarioName[2] = "ND50";
    ScenarioName[3] = "ND52";
    ScenarioName[4] = "ND54";
    ScenarioName[5] = "ND56";
    ScenarioName[6] = "ND58";
    ScenarioName[7] = "QS only";

    TString ScenarioShortName[NumScenarios];
    ScenarioShortName[0] = "ND46";
    ScenarioShortName[1] = "ND48";
    ScenarioShortName[2] = "ND50";
    ScenarioShortName[3] = "ND52";
    ScenarioShortName[4] = "ND54";
    ScenarioShortName[5] = "ND56";
    ScenarioShortName[6] = "ND58";
    ScenarioShortName[7] = "NoSI";
*/
/*
    const unsigned NumScenarios = 7;
    TString ScenarioName[NumScenarios];
    ScenarioName[0] = "NSC97a";
    ScenarioName[1] = "NSC97b";
    ScenarioName[2] = "NSC97c";
    ScenarioName[3] = "NSC97d";
    ScenarioName[4] = "NSC97e";
    ScenarioName[5] = "NSC97f";
    ScenarioName[6] = "QS only";

    TString ScenarioShortName[NumScenarios];
    ScenarioShortName[0] = "NSC97a";
    ScenarioShortName[1] = "NSC97b";
    ScenarioShortName[2] = "NSC97c";
    ScenarioShortName[3] = "NSC97d";
    ScenarioShortName[4] = "NSC97e";
    ScenarioShortName[5] = "NSC97f";
    ScenarioShortName[6] = "NoSI";
*/
/*
    const unsigned NumScenarios = 3;
    TString ScenarioName[NumScenarios];
    ScenarioName[0] = "Toy1";
    ScenarioName[1] = "Toy2";
    ScenarioName[2] = "Without strong interaction";

    TString ScenarioShortName[NumScenarios];
    ScenarioShortName[0] = "Toy1";
    ScenarioShortName[1] = "Toy2";
    ScenarioShortName[2] = "NoSI";
*/

    const unsigned NumScenarios = 5;
    TString ScenarioName[NumScenarios];
    ScenarioName[0] = "ND46";
    ScenarioName[1] = "NF44";
    ScenarioName[2] = "Ehime";
    ScenarioName[3] = "ESC08";
    ScenarioName[4] = "HKMYY";

    TString ScenarioShortName[NumScenarios];
    ScenarioShortName[0] = "ND46";
    ScenarioShortName[1] = "NF44";
    ScenarioShortName[2] = "Ehime";
    ScenarioShortName[3] = "ESC08";
    ScenarioShortName[4] = "HKMYY";

    double** PotPars = new double* [NumScenarios];
    double** ScatPars = new double* [NumScenarios];
    double** ScatParErr = new double* [NumScenarios];
    double** ScatParsMorita = new double* [NumScenarios];

    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        PotPars[uSce] = new double [8];
        ScatPars[uSce] = new double [2];
        ScatParErr[uSce] = new double [2];
        ScatParsMorita[uSce] = new double [2];
        for(unsigned uEl=0; uEl<8; uEl++) PotPars[uSce][uEl]=0;

        if(ScenarioName[uSce]=="Toy1")
        {//PotPars[uSce][2]=-140.0; PotPars[uSce][3]=2.0; PotPars[uSce][4]=340.0; PotPars[uSce][5]=0.45;
         PotPars[uSce][2]=-144.5; PotPars[uSce][3]=2.11; PotPars[uSce][4]=520.0; PotPars[uSce][5]=0.54;//the ones to use
         ScatParsMorita[uSce][0]=0.73; ScatParsMorita[uSce][1]=7.72;}

        if(ScenarioName[uSce]=="Toy2")
        {//PotPars[uSce][2]=-140.0; PotPars[uSce][3]=2.0; PotPars[uSce][4]=340.0; PotPars[uSce][5]=0.45;
         //PotPars[uSce][2]=-105; PotPars[uSce][3]=1.0; PotPars[uSce][4]=2550; PotPars[uSce][5]=0.46;
         //PotPars[uSce][2]=100; PotPars[uSce][3]=1.3; PotPars[uSce][4]=500; PotPars[uSce][5]=0.45;//THE ONES TO USE
         PotPars[uSce][2]=-30.0; PotPars[uSce][3]=2.0; PotPars[uSce][4]=1080.0; PotPars[uSce][5]=0.55;
         ScatParsMorita[uSce][0]=1.09; ScatParsMorita[uSce][1]=0.42;}

        if(ScenarioName[uSce]=="HKMYY")
        {PotPars[uSce][2]=-10.96; PotPars[uSce][3]=1.342; PotPars[uSce][4]=-141.75; PotPars[uSce][5]=0.777; PotPars[uSce][6]=2136.6; PotPars[uSce][7]=0.35;}

        if(ScenarioName[uSce]=="NF42")
        {PotPars[uSce][2]=-878.97; PotPars[uSce][3]=0.6; PotPars[uSce][4]=1048.58; PotPars[uSce][5]=0.45;}

        if(ScenarioName[uSce]=="NF44")
        {PotPars[uSce][2]=-1066.98; PotPars[uSce][3]=0.6; PotPars[uSce][4]=1646.65; PotPars[uSce][5]=0.45;
         ScatParsMorita[uSce][0]=23.956; ScatParsMorita[uSce][1]=1.258;}

        if(ScenarioName[uSce]=="NF50")
        {PotPars[uSce][2]=-2007.35; PotPars[uSce][3]=0.6; PotPars[uSce][4]=5678.97; PotPars[uSce][5]=0.45;
         ScatParsMorita[uSce][0]=-0.772; ScatParsMorita[uSce][1]=4.271;}

        if(ScenarioName[uSce]=="NF52")
        {PotPars[uSce][2]=-2276.73; PotPars[uSce][3]=0.6; PotPars[uSce][4]=7415.56; PotPars[uSce][5]=0.45;
         ScatParsMorita[uSce][0]=-0.406; ScatParsMorita[uSce][1]=8.868;}

        if(ScenarioName[uSce]=="ND46")
        {PotPars[uSce][2]=-144.89; PotPars[uSce][3]=1.0; PotPars[uSce][4]=127.87; PotPars[uSce][5]=0.45;
         ScatParsMorita[uSce][0]=4.621; ScatParsMorita[uSce][1]=1.3;}

        if(ScenarioName[uSce]=="ND48")
        {PotPars[uSce][2]=-150.83; PotPars[uSce][3]=1.0; PotPars[uSce][4]=355.09; PotPars[uSce][5]=0.45;
         ScatParsMorita[uSce][0]=14.394; ScatParsMorita[uSce][1]=1.633;}

        if(ScenarioName[uSce]=="ND50")
        {PotPars[uSce][2]=-151.54; PotPars[uSce][3]=1.0; PotPars[uSce][4]=587.21; PotPars[uSce][5]=0.45;
         ScatParsMorita[uSce][0]=-10.629; ScatParsMorita[uSce][1]=2.042;}

        if(ScenarioName[uSce]=="ND52")
        {PotPars[uSce][2]=-150.29; PotPars[uSce][3]=1.0; PotPars[uSce][4]=840.55; PotPars[uSce][5]=0.45;
         ScatParsMorita[uSce][0]=-3.483; ScatParsMorita[uSce][1]=2.592;}

        if(ScenarioName[uSce]=="ND54")
        {PotPars[uSce][2]=-147.65; PotPars[uSce][3]=1.0; PotPars[uSce][4]=1114.72; PotPars[uSce][5]=0.45;
         ScatParsMorita[uSce][0]=-1.893; ScatParsMorita[uSce][1]=3.389;}

        if(ScenarioName[uSce]=="ND56")
        {PotPars[uSce][2]=-144.26; PotPars[uSce][3]=1.0; PotPars[uSce][4]=1413.75; PotPars[uSce][5]=0.45;
         ScatParsMorita[uSce][0]=-1.179; ScatParsMorita[uSce][1]=4.656;}

        if(ScenarioName[uSce]=="ND58")
        {PotPars[uSce][2]=-137.74; PotPars[uSce][3]=1.0; PotPars[uSce][4]=1666.78; PotPars[uSce][5]=0.45;
         ScatParsMorita[uSce][0]=-0.764; ScatParsMorita[uSce][1]=6.863;}
         //ScatParsMorita[uSce][0]=1.1; ScatParsMorita[uSce][1]=8.52;}

        if(ScenarioName[uSce]=="NSC89-1020")
        {PotPars[uSce][2]=-22.89; PotPars[uSce][3]=1.0; PotPars[uSce][4]=67.45; PotPars[uSce][5]=0.45;
         ScatParsMorita[uSce][0]=-0.250; ScatParsMorita[uSce][1]=7.200;}

        if(ScenarioName[uSce]=="NSC89-920")
        {PotPars[uSce][2]=-1080.35; PotPars[uSce][3]=0.6; PotPars[uSce][4]=2039.54; PotPars[uSce][5]=0.45;
         ScatParsMorita[uSce][0]=-2.100; ScatParsMorita[uSce][1]=1.900;}

        if(ScenarioName[uSce]=="NSC89-820")
        {PotPars[uSce][2]=-1904.41; PotPars[uSce][3]=0.6; PotPars[uSce][4]=4996.93; PotPars[uSce][5]=0.45;
         ScatParsMorita[uSce][0]=-1.110; ScatParsMorita[uSce][1]=3.200;}

        if(ScenarioName[uSce]=="NSC97a")
        {PotPars[uSce][2]=-69.45; PotPars[uSce][3]=1.0; PotPars[uSce][4]=653.86; PotPars[uSce][5]=0.45;
         ScatParsMorita[uSce][0]=-0.329; ScatParsMorita[uSce][1]=12.370;}

        if(ScenarioName[uSce]=="NSC97b")
        {PotPars[uSce][2]=-78.42; PotPars[uSce][3]=1.0; PotPars[uSce][4]=741.76; PotPars[uSce][5]=0.45;
         ScatParsMorita[uSce][0]=-0.397; ScatParsMorita[uSce][1]=10.360;}

        if(ScenarioName[uSce]=="NSC97c")
        {PotPars[uSce][2]=-91.80; PotPars[uSce][3]=1.0; PotPars[uSce][4]=914.67; PotPars[uSce][5]=0.45;
         ScatParsMorita[uSce][0]=-0.476; ScatParsMorita[uSce][1]=9.130;}

        if(ScenarioName[uSce]=="NSC97d")
        {PotPars[uSce][2]=-445.77; PotPars[uSce][3]=0.4; PotPars[uSce][4]=373.64; PotPars[uSce][5]=0.30;
         ScatParsMorita[uSce][0]=-0.401; ScatParsMorita[uSce][1]=1.150;}

        if(ScenarioName[uSce]=="NSC97e")
        {PotPars[uSce][2]=-110.45; PotPars[uSce][3]=1.0; PotPars[uSce][4]=1309.55; PotPars[uSce][5]=0.45;
         ScatParsMorita[uSce][0]=-0.501; ScatParsMorita[uSce][1]=9.840;}

        if(ScenarioName[uSce]=="NSC97f")
        {PotPars[uSce][2]=-106.53; PotPars[uSce][3]=1.0; PotPars[uSce][4]=1469.33; PotPars[uSce][5]=0.45;
         ScatParsMorita[uSce][0]=-0.350; ScatParsMorita[uSce][1]=16.330;}

        if(ScenarioName[uSce]=="fss2")
        {PotPars[uSce][2]=-103.9; PotPars[uSce][3]=0.92; PotPars[uSce][4]=658.2; PotPars[uSce][5]=0.41;
         ScatParsMorita[uSce][0]=-0.810; ScatParsMorita[uSce][1]=3.99;}

        if(ScenarioName[uSce]=="Ehime")
        {PotPars[uSce][2]=-146.6; PotPars[uSce][3]=1.0; PotPars[uSce][4]=720.9; PotPars[uSce][5]=0.45;
         ScatParsMorita[uSce][0]=-4.21; ScatParsMorita[uSce][1]=2.41;}

        if(ScenarioName[uSce]=="ESC08")
        {PotPars[uSce][2]=-293.66; PotPars[uSce][3]=0.8; PotPars[uSce][4]=1429.27; PotPars[uSce][5]=0.45;
         ScatParsMorita[uSce][0]=-0.97; ScatParsMorita[uSce][1]=3.86;}

        if(ScenarioName[uSce]=="HKMYY")
        {PotPars[uSce][2]=-10.96; PotPars[uSce][3]=1.342; PotPars[uSce][4]=-141.75; PotPars[uSce][5]=0.777; PotPars[uSce][6]=2136.6; PotPars[uSce][7]=0.35;
         ScatParsMorita[uSce][0]=-0.575; ScatParsMorita[uSce][1]=6.45;}

        if(ScenarioName[uSce]=="Without strong interaction")
        {PotPars[uSce][2]=0; PotPars[uSce][3]=0; PotPars[uSce][4]=0; PotPars[uSce][5]=0; PotPars[uSce][6]=0; PotPars[uSce][7]=0;}
    }

    double LamLamSourcePars[4] = {0,0,0,Radius[0]};

    CATSparameters cSourcePars_LL(CATSparameters::tSource,1,true);
    cSourcePars_LL.SetParameters(&LamLamSourcePars[3]);

    CATS Nala;
    Nala.SetEpsilonConv(1e-7);
    Nala.SetEpsilonProp(1e-7);
    //Nala.SetThetaDependentSource(ThetaDepSource);
    Nala.SetUseAnalyticSource(true);
    if(Nala.GetThetaDependentSource()){
        if(SourceType=="Gauss") Nala.SetAnaSource(GaussSourceTheta, cSourcePars_LL);
        else if(SourceType=="Cauchy") Nala.SetAnaSource(CauchySourceTheta, cSourcePars_LL);
        else {printf("FUCK!!!\n"); return;}
    }
    else{
        if(SourceType=="Gauss") Nala.SetAnaSource(GaussSource, cSourcePars_LL);
        else if(SourceType=="Cauchy") Nala.SetAnaSource(CauchySource, cSourcePars_LL);
        else {printf("FUCK!!!\n"); return;}
    }
    Nala.SetExcludeFailedBins(false);
    Nala.SetMomBins(NumPlotPts,Plot_kMin,Plot_kMax);

    Nala.SetNumChannels(2);
    Nala.SetNumPW(0,1);
    Nala.SetNumPW(1,1);
    Nala.SetSpin(0,0);
    Nala.SetSpin(1,1);
    Nala.SetChannelWeight(0,1./4.);
    Nala.SetChannelWeight(1,3./4.);
    Nala.SetQ1Q2(0);
    Nala.SetPdgId(3122, 3122);
    Nala.SetRedMass( 0.5*Mass_L );

    unsigned NumFineSteps = 100;//100
    double kMinFine = 2;//2
    double kMaxFine = 40;//40
    CATS NalaFine;
    NalaFine.SetEpsilonConv(1e-9);
    NalaFine.SetEpsilonProp(1e-9);
    //Nala.SetThetaDependentSource(ThetaDepSource);
    NalaFine.SetUseAnalyticSource(true);
    if(NalaFine.GetThetaDependentSource()){
        if(SourceType=="Gauss") NalaFine.SetAnaSource(GaussSourceTheta, cSourcePars_LL);
        else if(SourceType=="Cauchy") NalaFine.SetAnaSource(CauchySourceTheta, cSourcePars_LL);
        else {printf("FUCK!!!\n"); return;}
    }
    else{
        if(SourceType=="Gauss") NalaFine.SetAnaSource(GaussSource, cSourcePars_LL);
        else if(SourceType=="Cauchy") NalaFine.SetAnaSource(CauchySource, cSourcePars_LL);
        else {printf("FUCK!!!\n"); return;}
    }
    NalaFine.SetExcludeFailedBins(false);
    NalaFine.SetMomBins(NumFineSteps,kMinFine,kMaxFine);

    NalaFine.SetNumChannels(2);
    NalaFine.SetNumPW(0,1);
    NalaFine.SetNumPW(1,1);
    NalaFine.SetSpin(0,0);
    NalaFine.SetSpin(1,1);
    NalaFine.SetChannelWeight(0,1./4.);
    NalaFine.SetChannelWeight(1,3./4.);
    NalaFine.SetQ1Q2(0);
    NalaFine.SetPdgId(3122, 3122);
    NalaFine.SetRedMass( 0.5*Mass_L );


    //in bins of Scenario/Radius
    TGraph** grCk;
    grCk = new TGraph* [NumScenarios];

    TGraph** grCkLedni;
    grCkLedni = new TGraph* [NumScenarios];

    const TString SigmaMatrixFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_7TeV/ResolutionMatrices/LHC14j4bcde_MeV.root";
    const double TrMatUnitConv = 1;//1000

    TH2F* hSigma_pL;
    TFile* FileSigma = new TFile(SigmaMatrixFileName, "read");
    FileSigma->cd();
    hSigma_pL = (TH2F*)FileSigma->Get("hSigma_pL_InMeV");

    //1/FractionOfBins th number of original bins of the correction matrix are taken into account
    //originally: 3000 MeV => 1/10 should be good most of the time
    const int FractionOfBins = 6;

    TH2F* hSigma_pL_InFuckingMeV = new TH2F("hSigma_pL_InFuckingMeV", "hSigma_pL_InFuckingMeV",
        hSigma_pL->GetNbinsX()/FractionOfBins, hSigma_pL->GetXaxis()->GetBinLowEdge(1)*TrMatUnitConv,hSigma_pL->GetXaxis()->GetBinUpEdge(hSigma_pL->GetNbinsX()/FractionOfBins)*TrMatUnitConv,
        hSigma_pL->GetNbinsY()/FractionOfBins, hSigma_pL->GetYaxis()->GetBinLowEdge(1)*TrMatUnitConv,hSigma_pL->GetXaxis()->GetBinUpEdge(hSigma_pL->GetNbinsY()/FractionOfBins)*TrMatUnitConv);

    for(int iBinX=1; iBinX<=hSigma_pL->GetNbinsX()/FractionOfBins; iBinX++){
        for(int iBinY=1; iBinY<=hSigma_pL->GetNbinsY()/FractionOfBins; iBinY++){
            hSigma_pL_InFuckingMeV->SetBinContent(iBinX, iBinY, hSigma_pL->GetBinContent(iBinX, iBinY));
        }
    }

    //in bins of Scenario/Radius
    TGraph** grCkCorr;
    grCkCorr = new TGraph* [NumScenarios];

    TGraph** grCkLedniCorr;
    grCkLedniCorr = new TGraph* [NumScenarios];

    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        grCk[uSce] = new TGraph[NumRad];
        grCkLedni[uSce] = new TGraph[NumRad];
        grCkCorr[uSce] = new TGraph[NumRad];
        grCkLedniCorr[uSce] = new TGraph[NumRad];
        Nala.SetShortRangePotential(0,0,TripleGaussSum,PotPars[uSce]);
        Nala.KillTheCat(CATS::kAllChanged);

        NalaFine.SetShortRangePotential(0,0,TripleGaussSum,PotPars[uSce]);
        NalaFine.KillTheCat(CATS::kAllChanged);

        for(unsigned uRad=0; uRad<NumRad; uRad++){
            Nala.SetAnaSource(0, Radius[uRad]);
            Nala.KillTheCat();

            DLM_Ck* Ck_LL = new DLM_Ck(1,0,Nala);
            Ck_LL->Update();

            DLM_Ck* CkLedni_LL = new DLM_Ck(1,2,NumPlotPts,Plot_kMin,Plot_kMax,Lednicky_Identical_Singlet);
            CkLedni_LL->SetPotPar(0,-ScatParsMorita[uSce][0]);
            CkLedni_LL->SetPotPar(1,ScatParsMorita[uSce][1]);
            CkLedni_LL->SetSourcePar(0,Radius[uRad]);
            CkLedni_LL->Update();

//Save the uncorrected correlation
            grCk[uSce][uRad].Set(NumPlotPts);
            grCk[uSce][uRad].SetName(TString::Format("grCk_%u_%u",uSce,uRad));
            grCk[uSce][uRad].SetMarkerColor(DlmCol.GetColor(1+uSce));
            grCk[uSce][uRad].SetMarkerStyle(20);
            grCk[uSce][uRad].SetMarkerSize(0);
            grCk[uSce][uRad].SetLineColor(DlmCol.GetColor(1+uSce));
            grCk[uSce][uRad].SetLineWidth(8.-2.*double(uSce)/double(NumScenarios-1));
            if(uSce==NumScenarios-1){
                grCk[uSce][uRad].SetLineWidth(8.);
                grCk[uSce][uRad].SetLineColor(kGray+1);
                grCk[uSce][uRad].SetLineStyle(2);
            }
            grCk[uSce][uRad].SetFillColor(kWhite);

            for(unsigned uMomBin=0; uMomBin<NumPlotPts; uMomBin++){
                grCk[uSce][uRad].SetPoint(uMomBin, Nala.GetMomentum(uMomBin), Nala.GetCorrFun(uMomBin));
            }

            grCkLedni[uSce][uRad].Set(NumPlotPts);
            grCkLedni[uSce][uRad].SetName(TString::Format("grCkLedni_%u_%u",uSce,uRad));
            grCkLedni[uSce][uRad].SetMarkerColor(DlmCol.GetColor(1+uSce));
            grCkLedni[uSce][uRad].SetMarkerStyle(20);
            grCkLedni[uSce][uRad].SetMarkerSize(0);
            grCkLedni[uSce][uRad].SetLineColor(DlmCol.GetColor(1+uSce));
            grCkLedni[uSce][uRad].SetLineWidth(8.-2.*double(uSce)/double(NumScenarios-1));
            grCkLedni[uSce][uRad].SetLineStyle(2);
            if(uSce==NumScenarios-1){
                grCkLedni[uSce][uRad].SetLineWidth(8.);
                grCkLedni[uSce][uRad].SetLineColor(kGray+1);
                grCkLedni[uSce][uRad].SetLineStyle(2);
            }
            grCkLedni[uSce][uRad].SetFillColor(kWhite);

            for(unsigned uMomBin=0; uMomBin<NumPlotPts; uMomBin++){
                grCkLedni[uSce][uRad].SetPoint(uMomBin, Nala.GetMomentum(uMomBin), CkLedni_LL->Eval(Nala.GetMomentum(uMomBin)));
            }

//Apply the correlations

            DLM_CkDecomposition CkDec_LL("LambdaLambda",2,*Ck_LL,hSigma_pL_InFuckingMeV);
            CkDec_LL.AddContribution(0,lambdaFeed,DLM_CkDecomposition::cFeedDown,NULL);
            CkDec_LL.AddContribution(1,lambdaPurity,DLM_CkDecomposition::cFake);
            CkDec_LL.Update();

            DLM_CkDecomposition CkDecLedni_LL("LambdaLambda",2,*CkLedni_LL,hSigma_pL_InFuckingMeV);
            CkDecLedni_LL.AddContribution(0,lambdaFeed,DLM_CkDecomposition::cFeedDown,NULL);
            CkDecLedni_LL.AddContribution(1,lambdaPurity,DLM_CkDecomposition::cFake);
            CkDecLedni_LL.Update();

            grCkCorr[uSce][uRad].Set(NumPlotPts);
            grCkCorr[uSce][uRad].SetName(TString::Format("grCkCorr_%u_%u",uSce,uRad));
            grCkCorr[uSce][uRad].SetMarkerColor(DlmCol.GetColor(1+uSce));
            grCkCorr[uSce][uRad].SetMarkerStyle(20);
            grCkCorr[uSce][uRad].SetMarkerSize(0);
            grCkCorr[uSce][uRad].SetLineColor(DlmCol.GetColor(1+uSce));
            grCkCorr[uSce][uRad].SetLineWidth(8.-2.*double(uSce)/double(NumScenarios-1));
            if(uSce==NumScenarios-1){
                grCkCorr[uSce][uRad].SetLineWidth(8.);
                grCkCorr[uSce][uRad].SetLineColor(kGray+1);
                grCkCorr[uSce][uRad].SetLineStyle(2);
            }
            grCkCorr[uSce][uRad].SetFillColor(kWhite);

            for(unsigned uMomBin=0; uMomBin<NumPlotPts; uMomBin++){
                grCkCorr[uSce][uRad].SetPoint(uMomBin, Nala.GetMomentum(uMomBin), CkDec_LL.EvalCk(Nala.GetMomentum(uMomBin)));
            }

            grCkLedniCorr[uSce][uRad].Set(NumPlotPts);
            grCkLedniCorr[uSce][uRad].SetName(TString::Format("grCkLedniCorr_%u_%u",uSce,uRad));
            grCkLedniCorr[uSce][uRad].SetMarkerColor(DlmCol.GetColor(1+uSce));
            grCkLedniCorr[uSce][uRad].SetMarkerStyle(20);
            grCkLedniCorr[uSce][uRad].SetMarkerSize(0);
            grCkLedniCorr[uSce][uRad].SetLineColor(DlmCol.GetColor(1+uSce));
            grCkLedniCorr[uSce][uRad].SetLineWidth(8.-2.*double(uSce)/double(NumScenarios-1));
            grCkLedniCorr[uSce][uRad].SetLineStyle(2);
            if(uSce==NumScenarios-1){
                grCkLedniCorr[uSce][uRad].SetLineWidth(8.);
                grCkLedniCorr[uSce][uRad].SetLineColor(kGray+1);
                grCkLedniCorr[uSce][uRad].SetLineStyle(2);
            }
            grCkLedniCorr[uSce][uRad].SetFillColor(kWhite);

            for(unsigned uMomBin=0; uMomBin<NumPlotPts; uMomBin++){
                grCkLedniCorr[uSce][uRad].SetPoint(uMomBin, Nala.GetMomentum(uMomBin), CkDecLedni_LL.EvalCk(Nala.GetMomentum(uMomBin)));
            }

            delete Ck_LL;
            delete CkLedni_LL;


        }

        double kFineStep = (kMaxFine-kMinFine)/double(NumFineSteps-1);
        TH1F* hPhaseShift0 = new TH1F("hPhaseShift0", "hPhaseShift0", NumFineSteps, kMinFine-0.5*kFineStep, kMaxFine+0.5*kFineStep);
        for(unsigned iBin=0; iBin<NumFineSteps; iBin++){
            hPhaseShift0->SetBinContent(iBin+1, NalaFine.GetMomentum(iBin)/tan(NalaFine.GetPhaseShift(iBin,0,0)));
            //printf("k = %.2f; delta = %.3f; whatever/tan=%.3f\n",NalaFine.GetMomentum(iBin),NalaFine.GetPhaseShift(iBin,0,0)*180./3.14,
            //       NalaFine.GetMomentum(iBin)/tan(NalaFine.GetPhaseShift(iBin,0,0)));
            hPhaseShift0->SetBinError(iBin+1, hPhaseShift0->GetBinContent(iBin+1)*0.0001);
        }


        TF1* fPhaseShift0 = new TF1("fPhaseShift0", "[2]*x*x*x*x+0.5*[1]*x*x-1./[0]", kMinFine, kMaxFine);
        fPhaseShift0->SetParameter(0, ScatParsMorita[uSce][0]/197.327);
        fPhaseShift0->SetParLimits(0, -30./197.327, 30./197.327);
        fPhaseShift0->SetParameter(1, ScatParsMorita[uSce][1]/197.327);
        fPhaseShift0->SetParLimits(1, 0./197.327, 64./197.327);
        fPhaseShift0->FixParameter(2, 0);

        //fPhaseShift0->SetLineColor(kRed);
        //fPhaseShift0->SetLineWidth(8);

        hPhaseShift0->Fit(fPhaseShift0, "S, N, R, M");
        ScatPars[uSce][0] = fPhaseShift0->GetParameter(0)*197.327;
        ScatPars[uSce][1] = fPhaseShift0->GetParameter(1)*197.327;
        ScatParErr[uSce][0] = fPhaseShift0->GetParError(0)*197.327;
        ScatParErr[uSce][1] = fPhaseShift0->GetParError(1)*197.327;

        delete hPhaseShift0;
        delete fPhaseShift0;
    }

    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        printf("Morita vs CATS:\n");
        printf(" a0(%s) = %.2f vs %.2f+/-%.2f\n",ScenarioShortName[uSce].Data(), ScatParsMorita[uSce][0], ScatPars[uSce][0], ScatParErr[uSce][0]);
        printf(" reff(%s) = %.2f vs %.2f+/-%.2f\n",ScenarioShortName[uSce].Data(), ScatParsMorita[uSce][1], ScatPars[uSce][1], ScatParErr[uSce][1]);

    }

    char* fBuffer = new char [512];
    FILE * CkFile;
    CkFile = fopen (TString::Format("%sPlotCk_LL_ver2.txt",OutputDir.Data()),"w");
    fprintf (CkFile, "%18s","k (MeV)");
    for(unsigned uRad=0; uRad<NumRad; uRad++){
        for(unsigned uSce=0; uSce<NumScenarios; uSce++){
            sprintf(fBuffer,"R%.2f_%s",Radius[uRad],ScenarioShortName[uSce].Data());
            fprintf (CkFile,"%18s",fBuffer);
            //fprintf (CkFile,"%18s","<-RatioNLO");
        }
    }

    double xVal,yVal;
    for(unsigned uMomBin=0; uMomBin<NumPlotPts; uMomBin++){
        fprintf (CkFile, "\n%18.3f",Nala.GetMomentum(uMomBin));
        for(unsigned uRad=0; uRad<NumRad; uRad++){
            for(unsigned uSce=0; uSce<NumScenarios; uSce++){
                grCk[uSce][uRad].GetPoint(uMomBin,xVal,yVal);
                if(xVal!=Nala.GetMomentum(uMomBin)) {printf("WTF MAN LL!\n");}
                fprintf (CkFile, "%18.3f",yVal);
                //grCkRatio[uSce][uRad].GetPoint(uMomBin,xVal,yVal);
                //fprintf (CkFile, "%18.3f",yVal);
            }
        }
    }

    FILE * CkFileCorr;
    CkFileCorr = fopen (TString::Format("%sPlotCk_LL_ver2_Corrected1.txt",OutputDir.Data()),"w");
    fprintf (CkFileCorr, "%18s","k (MeV)");
    for(unsigned uRad=0; uRad<NumRad; uRad++){
        for(unsigned uSce=0; uSce<NumScenarios; uSce++){
            sprintf(fBuffer,"R%.2f_%s",Radius[uRad],ScenarioShortName[uSce].Data());
            fprintf (CkFileCorr,"%18s",fBuffer);
        }
    }

    for(unsigned uMomBin=0; uMomBin<NumPlotPts; uMomBin++){
        fprintf (CkFileCorr, "\n%18.3f",Nala.GetMomentum(uMomBin));
        for(unsigned uRad=0; uRad<NumRad; uRad++){
            for(unsigned uSce=0; uSce<NumScenarios; uSce++){
                grCkCorr[uSce][uRad].GetPoint(uMomBin,xVal,yVal);
                if(xVal!=Nala.GetMomentum(uMomBin)) {printf("WTF MAN LL!\n");}
                fprintf (CkFileCorr, "%18.3f",yVal);
            }
        }
    }

    const double Yoffset = 1.15;
    TH1F* hAxis = new TH1F("hAxis", "hAxis", NumPlotPts,Plot_kMin,Plot_kMax);
    hAxis->SetStats(false);
    hAxis->SetTitle("");
    hAxis->GetXaxis()->SetLabelSize(0.065);
    hAxis->GetXaxis()->SetTitle("k (MeV)");
    hAxis->GetXaxis()->CenterTitle();
    hAxis->GetXaxis()->SetTitleOffset(1.05);
    hAxis->GetXaxis()->SetLabelOffset(0.02);
    hAxis->GetXaxis()->SetTitleSize(0.075);
    hAxis->GetYaxis()->SetLabelSize(0.065);
    hAxis->GetYaxis()->SetTitle("C(k)");
    hAxis->GetYaxis()->CenterTitle();
    hAxis->GetYaxis()->SetTitleOffset(Yoffset);
    hAxis->GetYaxis()->SetTitleSize(0.075);

    //hAxis->GetYaxis()->SetRangeUser(0.2, 1.3);
    //hAxis->GetYaxis()->SetLimits(0.2, 1.3);
hAxis->GetYaxis()->SetRangeUser(0.4, 1.8);
hAxis->GetYaxis()->SetLimits(0.4, 1.8);
//hAxis->GetYaxis()->SetRangeUser(0.5, 1.1);
//hAxis->GetYaxis()->SetLimits(0.5, 1.1);


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
    hAxisZoom->GetYaxis()->SetRangeUser(0.85, 1.15);
    hAxisZoom->GetYaxis()->SetLimits(0.85, 1.15);
    hAxisZoom->GetYaxis()->SetNdivisions(603);

    TPaveText** myPT = new TPaveText* [NumRad];
    for(unsigned uRad=0; uRad<NumRad; uRad++){
        myPT[uRad] = new TPaveText(0.5,0.885,0.95,0.975, "blNDC");//lbrt
        myPT[uRad]->SetName("PT12");
        myPT[uRad]->SetBorderSize(1);
        myPT[uRad]->SetTextSize(0.08/1.4);
        myPT[uRad]->SetFillColor(kWhite);
        myPT[uRad]->SetTextFont(22);
        //myPT[uRad]->AddText("Gauss source");
        myPT[uRad]->AddText(TString::Format("R = %.2f fm", Radius[uRad]));
    }

    TF1* UnitLine = new TF1("UnitLine", "1",Plot_kMin,Plot_kMax);
    UnitLine->SetLineColor(kGray+1);
    UnitLine->SetLineWidth(6);

    DLM_SubPads DrawBoard_Ck_LL(2560,1080);
    double PadLength = 1./double(NumRad);
    double CurrentPos = 0;
    for(unsigned uRad=0; uRad<NumRad; uRad++){
        DrawBoard_Ck_LL.AddSubPadTL(CurrentPos,CurrentPos+PadLength,0.08,1.);//lrtb
        //DrawBoard_Ck_LL.AddSubPadTL(CurrentPos,CurrentPos+PadLength,0.63,1);//lrtb
        CurrentPos+=PadLength;
//printf("uRad=%u; %.2f <--> %.2f\n", uRad, );
    }
    DrawBoard_Ck_LL.AddSubPadTL(0,1,0,0.075);//lrtb

    for(unsigned uRad=0; uRad<NumRad; uRad++){
        DrawBoard_Ck_LL.SetMargin(uRad, 0.07, 0.01, 0.15, 0.04);//lrbt
    }
    DrawBoard_Ck_LL.SetMargin(NumRad, 0, 0, 0, 0);//lrbt

    for(unsigned uRad=0; uRad<NumRad; uRad++){
        DrawBoard_Ck_LL.GetCanvas()->SetName(TString::Format("DrawBoard_Ck_LL_%u",uRad));
        DrawBoard_Ck_LL.cd(uRad);
        hAxis->Draw("AXIS");
        //UnitLine->Draw("same");
        for(unsigned uSce=0; uSce<NumScenarios; uSce++) grCk[uSce][uRad].Draw("same,C");
        for(unsigned uSce=0; uSce<NumScenarios; uSce++) grCkLedni[uSce][uRad].Draw("same,C");
        for(unsigned uSce=0; uSce<NumScenarios; uSce++) myPT[uRad]->Draw("same");

        DrawBoard_Ck_LL.SetLabelSize(uRad, hAxis->GetXaxis(), 17);
        DrawBoard_Ck_LL.SetLabelSize(uRad, hAxis->GetYaxis(), 17);
        DrawBoard_Ck_LL.SetTitleSize(uRad, hAxis->GetXaxis(), 17);
        DrawBoard_Ck_LL.SetTitleSize(uRad, hAxis->GetYaxis(), 17);
    }

    TLegend* myLegend = new TLegend(0.12,0,0.95,0.96);//lbrt
    myLegend->SetName("myLegend");
    myLegend->SetTextSize(0.65);
    myLegend->SetNColumns(NumScenarios);
    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        myLegend->AddEntry(&grCk[uSce][0],  ScenarioName[uSce]);
    }

    DrawBoard_Ck_LL.cd(NumRad);
    myLegend->Draw();

    gStyle->SetLineWidth(2.5);
    DrawBoard_Ck_LL.GetCanvas()->SaveAs(TString::Format("%sPlotCk_LL_ver2.png",OutputDir.Data()));
    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        for(unsigned uRad=0; uRad<NumRad; uRad++){
            grCk[uSce][uRad].SetLineWidth(grCk[uSce][uRad].GetLineWidth()/2.5);
            grCkLedni[uSce][uRad].SetLineWidth(grCkLedni[uSce][uRad].GetLineWidth()/2.5);
        }
    }
    hAxis->SetLineWidth(hAxis->GetLineWidth()/2.5);
    UnitLine->SetLineWidth(UnitLine->GetLineWidth()/2.5);
    gStyle->SetLineWidth(1);
    DrawBoard_Ck_LL.GetCanvas()->SaveAs(TString::Format("%sPlotCk_LL_ver2.pdf",OutputDir.Data()));

    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        for(unsigned uRad=0; uRad<NumRad; uRad++){
            grCk[uSce][uRad].SetLineWidth(grCk[uSce][uRad].GetLineWidth()*2.5);
            grCkLedni[uSce][uRad].SetLineWidth(grCkLedni[uSce][uRad].GetLineWidth()*2.5);
        }
    }
    hAxis->SetLineWidth(hAxis->GetLineWidth()*2.5);
    UnitLine->SetLineWidth(UnitLine->GetLineWidth()*2.5);


    DLM_SubPads DrawBoard_Ck_LL_Corrected1(2560,1080);
    CurrentPos=0;
    for(unsigned uRad=0; uRad<NumRad; uRad++){
        DrawBoard_Ck_LL_Corrected1.AddSubPadTL(CurrentPos,CurrentPos+PadLength,0.08,1.);//lrtb
        CurrentPos+=PadLength;
    }
    DrawBoard_Ck_LL_Corrected1.AddSubPadTL(0,1,0,0.075);//lrtb

    for(unsigned uRad=0; uRad<NumRad; uRad++){
        DrawBoard_Ck_LL_Corrected1.SetMargin(uRad, 0.07, 0.01, 0.15, 0.04);//lrbt
    }
    DrawBoard_Ck_LL_Corrected1.SetMargin(NumRad, 0, 0, 0, 0);//lrbt

    for(unsigned uRad=0; uRad<NumRad; uRad++){
        DrawBoard_Ck_LL_Corrected1.GetCanvas()->SetName(TString::Format("DrawBoard_Ck_LL_%u",uRad));
        DrawBoard_Ck_LL_Corrected1.cd(uRad);
        hAxis->Draw("AXIS");
        //UnitLine->Draw("same");
        for(unsigned uSce=0; uSce<NumScenarios; uSce++) grCkCorr[uSce][uRad].Draw("same,C");
        for(unsigned uSce=0; uSce<NumScenarios; uSce++) grCkLedniCorr[uSce][uRad].Draw("same,C");
        for(unsigned uSce=0; uSce<NumScenarios; uSce++) myPT[uRad]->Draw("same");

        DrawBoard_Ck_LL_Corrected1.SetLabelSize(uRad, hAxis->GetXaxis(), 17);
        DrawBoard_Ck_LL_Corrected1.SetLabelSize(uRad, hAxis->GetYaxis(), 17);
        DrawBoard_Ck_LL_Corrected1.SetTitleSize(uRad, hAxis->GetXaxis(), 17);
        DrawBoard_Ck_LL_Corrected1.SetTitleSize(uRad, hAxis->GetYaxis(), 17);
    }

    DrawBoard_Ck_LL_Corrected1.cd(NumRad);
    myLegend->Draw();

    TPaveText* ptCorrected0 = new TPaveText(0.4,0.47,0.95,0.60, "blNDC");//lbrt
    ptCorrected0->SetName("ptCorrected0");
    ptCorrected0->SetBorderSize(1);
    ptCorrected0->SetTextSize(0.075/1.4);
    ptCorrected0->SetFillColor(kWhite);
    ptCorrected0->SetTextFont(22);
    ptCorrected0->AddText("Corrections were");
    ptCorrected0->AddText("applied for:");

    TPaveText* ptCorrected = new TPaveText(0.4,0.19,0.95,0.4, "blNDC");//lbrt
    ptCorrected->SetName("ptCorrected");
    ptCorrected->SetBorderSize(1);
    ptCorrected->SetTextSize(0.065/1.4);
    ptCorrected->SetFillColor(kWhite);
    ptCorrected->SetTextFont(22);
    ptCorrected->SetTextAlign(11);
    //ptCorrected->AddText("* Momentum resolution");
    //ptCorrected->AddText("* Residuals assuming:");
    ptCorrected->AddText(TString::Format("   %.0f%% primary #Lambda#Lambda  pairs",lambdaPar*100));
    ptCorrected->AddText(TString::Format("   %.0f%% flat residuals",(1.-lambdaPar)*100));

    DrawBoard_Ck_LL_Corrected1.cd(NumRad/2);
    //ptCorrected0->Draw("same");
    ptCorrected->Draw("same");



    gStyle->SetLineWidth(2.5);
    DrawBoard_Ck_LL_Corrected1.GetCanvas()->SaveAs(TString::Format("%sPlotCk_LL_ver2_Corrected1.png",OutputDir.Data()));
    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        for(unsigned uRad=0; uRad<NumRad; uRad++){
            grCk[uSce][uRad].SetLineWidth(grCk[uSce][uRad].GetLineWidth()/2.5);
            grCkLedni[uSce][uRad].SetLineWidth(grCkLedni[uSce][uRad].GetLineWidth()/2.5);
            grCkCorr[uSce][uRad].SetLineWidth(grCkCorr[uSce][uRad].GetLineWidth()/2.5);
            grCkLedniCorr[uSce][uRad].SetLineWidth(grCkLedniCorr[uSce][uRad].GetLineWidth()/2.5);
        }
    }
    hAxis->SetLineWidth(hAxis->GetLineWidth()/2.5);
    UnitLine->SetLineWidth(UnitLine->GetLineWidth()/2.5);
    gStyle->SetLineWidth(1);
    DrawBoard_Ck_LL_Corrected1.GetCanvas()->SaveAs(TString::Format("%sPlotCk_LL_ver2_Corrected1.pdf",OutputDir.Data()));

    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        for(unsigned uRad=0; uRad<NumRad; uRad++){
            grCk[uSce][uRad].SetLineWidth(grCk[uSce][uRad].GetLineWidth()*2.5);
            grCkLedni[uSce][uRad].SetLineWidth(grCkLedni[uSce][uRad].GetLineWidth()*2.5);
            grCkCorr[uSce][uRad].SetLineWidth(grCkCorr[uSce][uRad].GetLineWidth()*2.5);
            grCkLedniCorr[uSce][uRad].SetLineWidth(grCkLedniCorr[uSce][uRad].GetLineWidth()*2.5);
        }
    }
    hAxis->SetLineWidth(hAxis->GetLineWidth()*2.5);
    UnitLine->SetLineWidth(UnitLine->GetLineWidth()*2.5);



    for(unsigned uRad=0; uRad<NumRad; uRad++){
        DLM_SubPads DrawBoard_Ck_LL_Compare(2560,1080);

        DrawBoard_Ck_LL_Compare.AddSubPadTL(0,0.5,0.08,1.);//lrtb
        DrawBoard_Ck_LL_Compare.AddSubPadTL(0.5,1,0.08,1.);//lrtb
        DrawBoard_Ck_LL_Compare.AddSubPadTL(0,1,0,0.075);//lrtb

        for(unsigned uPlot=0; uPlot<2; uPlot++){
            DrawBoard_Ck_LL_Compare.SetMargin(uPlot, 0.07, 0.01, 0.15, 0.04);//lrbt
        }
        DrawBoard_Ck_LL_Compare.SetMargin(2, 0, 0, 0, 0);//lrbt


        DrawBoard_Ck_LL_Compare.GetCanvas()->SetName(TString::Format("DrawBoard_Ck_LL_Compare_%u",uRad));
        DrawBoard_Ck_LL_Compare.cd(0);
        hAxis->Draw("AXIS");
        for(unsigned uSce=0; uSce<NumScenarios; uSce++) grCk[uSce][uRad].Draw("same,C");
        for(unsigned uSce=0; uSce<NumScenarios; uSce++) grCkLedni[uSce][uRad].Draw("same,C");
        for(unsigned uSce=0; uSce<NumScenarios; uSce++) myPT[uRad]->Draw("same");
        DrawBoard_Ck_LL_Compare.SetLabelSize(0, hAxis->GetXaxis(), 17);
        DrawBoard_Ck_LL_Compare.SetLabelSize(0, hAxis->GetYaxis(), 17);
        DrawBoard_Ck_LL_Compare.SetTitleSize(0, hAxis->GetXaxis(), 17);
        DrawBoard_Ck_LL_Compare.SetTitleSize(0, hAxis->GetYaxis(), 17);

        DrawBoard_Ck_LL_Compare.cd(1);
        hAxis->Draw("AXIS");
        for(unsigned uSce=0; uSce<NumScenarios; uSce++) grCkCorr[uSce][uRad].Draw("same,C");
        for(unsigned uSce=0; uSce<NumScenarios; uSce++) grCkLedniCorr[uSce][uRad].Draw("same,C");
        for(unsigned uSce=0; uSce<NumScenarios; uSce++) myPT[uRad]->Draw("same");
        ptCorrected->Draw("same");
        DrawBoard_Ck_LL_Compare.SetLabelSize(1, hAxis->GetXaxis(), 17);
        DrawBoard_Ck_LL_Compare.SetLabelSize(1, hAxis->GetYaxis(), 17);
        DrawBoard_Ck_LL_Compare.SetTitleSize(1, hAxis->GetXaxis(), 17);
        DrawBoard_Ck_LL_Compare.SetTitleSize(1, hAxis->GetYaxis(), 17);

        DrawBoard_Ck_LL_Compare.cd(2);
        myLegend->Draw();

        gStyle->SetLineWidth(2.5);
        DrawBoard_Ck_LL_Compare.GetCanvas()->SaveAs(TString::Format("%sDrawBoard_Ck_LL_Compare_%u.png",OutputDir.Data(),uRad));
        for(unsigned uSce=0; uSce<NumScenarios; uSce++){
            for(unsigned uRad=0; uRad<NumRad; uRad++){
                grCk[uSce][uRad].SetLineWidth(grCk[uSce][uRad].GetLineWidth()/2.5);
                grCkLedni[uSce][uRad].SetLineWidth(grCkLedni[uSce][uRad].GetLineWidth()/2.5);
                grCkCorr[uSce][uRad].SetLineWidth(grCkCorr[uSce][uRad].GetLineWidth()/2.5);
                grCkLedniCorr[uSce][uRad].SetLineWidth(grCkLedniCorr[uSce][uRad].GetLineWidth()/2.5);
            }
        }
        hAxis->SetLineWidth(hAxis->GetLineWidth()/2.5);
        UnitLine->SetLineWidth(UnitLine->GetLineWidth()/2.5);
        gStyle->SetLineWidth(1);
        DrawBoard_Ck_LL_Compare.GetCanvas()->SaveAs(TString::Format("%sDrawBoard_Ck_LL_Compare_%u.pdf",OutputDir.Data(),uRad));

        for(unsigned uSce=0; uSce<NumScenarios; uSce++){
            for(unsigned uRad=0; uRad<NumRad; uRad++){
                grCk[uSce][uRad].SetLineWidth(grCk[uSce][uRad].GetLineWidth()*2.5);
                grCkLedni[uSce][uRad].SetLineWidth(grCkLedni[uSce][uRad].GetLineWidth()*2.5);
                grCkCorr[uSce][uRad].SetLineWidth(grCkCorr[uSce][uRad].GetLineWidth()*2.5);
                grCkLedniCorr[uSce][uRad].SetLineWidth(grCkLedniCorr[uSce][uRad].GetLineWidth()*2.5);
            }
        }
        hAxis->SetLineWidth(hAxis->GetLineWidth()*2.5);
        UnitLine->SetLineWidth(UnitLine->GetLineWidth()*2.5);

    }



    //!PLOT THE POTENTIALS
    double RAD;
    TGraph* gPot_LL = new TGraph [NumScenarios];
    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        gPot_LL[uSce].Set(NumPotentialPts);
        gPot_LL[uSce].SetName(TString::Format("gPot_LL_%s",ScenarioName[uSce].Data()));
        gPot_LL[uSce].SetLineColor(DlmCol.GetColor(uSce+1));
        gPot_LL[uSce].SetLineWidth(12.-4.*double(uSce)/double(NumScenarios-1));
        gPot_LL[uSce].SetMarkerStyle(1);
        gPot_LL[uSce].SetMarkerSize(0);
        gPot_LL[uSce].SetMarkerColor(DlmCol.GetColor(uSce+1));
        gPot_LL[uSce].SetFillColor(kWhite);
        if(uSce==NumScenarios-1){
            gPot_LL[uSce].SetLineWidth(8.);
            gPot_LL[uSce].SetLineColor(kGray+1);
            gPot_LL[uSce].SetLineStyle(2);
        }

        for(unsigned uRad=0; uRad<NumPotentialPts; uRad++){
            RAD = rMin + double(uRad)*rStep;
            PotPars[uSce][0] = RAD;
            gPot_LL[uSce].SetPoint(uRad,RAD,TripleGaussSum(PotPars[uSce]));
        }
    }

    FILE * PotFile;
    PotFile = fopen (TString::Format("%sPlotPot_LL_ver2.txt",OutputDir.Data()),"w");
    fprintf (PotFile, "%18s","r (fm)");
    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        sprintf(fBuffer,"%s",ScenarioShortName[uSce].Data());
        fprintf (PotFile,"%18s",fBuffer);
    }

    for(unsigned uRad=0; uRad<NumPotentialPts; uRad++){
        RAD = rMin + double(uRad)*rStep;
        fprintf (PotFile, "\n%18.3f",RAD);
        for(unsigned uSce=0; uSce<NumScenarios; uSce++){
            gPot_LL[uSce].GetPoint(uRad,xVal,yVal);
            if(xVal!=RAD) {printf("WTF MAN LL pot!\n");}
            fprintf (PotFile, "%18.3f",yVal);
        }
    }


    TH1F* hPotAxis = new TH1F("hPotAxis", "hPotAxis", NumPotentialPts,rMin,rMax);
    hPotAxis->SetStats(false);
    hPotAxis->SetTitle("");
    hPotAxis->GetXaxis()->SetLabelSize(0.065);
    hPotAxis->GetXaxis()->SetTitle("r (fm)");
    hPotAxis->GetXaxis()->CenterTitle();
    hPotAxis->GetXaxis()->SetTitleOffset(1.15);
    hPotAxis->GetXaxis()->SetLabelOffset(0.02);
    hPotAxis->GetXaxis()->SetTitleSize(0.075);
    hPotAxis->GetYaxis()->SetLabelSize(0.065);
    hPotAxis->GetYaxis()->SetTitle("V (MeV)");
    hPotAxis->GetYaxis()->CenterTitle();
    hPotAxis->GetYaxis()->SetTitleOffset(0.75);
    hPotAxis->GetYaxis()->SetTitleSize(0.075);
    hPotAxis->GetYaxis()->SetRangeUser(-200, 400);
    hPotAxis->GetYaxis()->SetLimits(-200, 400);
    hPotAxis->GetYaxis()->SetNdivisions(405,false);

    TLegend* PotLegend = new TLegend(0.8,0.6,0.975,0.88);//lbrt
    PotLegend->SetName("PotLegend");
    PotLegend->SetTextSize(0.05);
    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        if(uSce!=NumScenarios-1) PotLegend->AddEntry(&gPot_LL[uSce], TString::Format("%s", ScenarioName[uSce].Data()));
    }


    //TF1* ZeroLine = new TF1("ZeroLine", "0", rMin, rMax);
    //ZeroLine->SetLineColor(kGray);
    //ZeroLine->SetLineWidth(6);

    TCanvas* cPot_LL = new TCanvas("cPot_LL", "cPot_LL", 1);
    cPot_LL->cd(0); cPot_LL->SetCanvasSize(1920, 1080); cPot_LL->SetMargin(0.12,0.025,0.18,0.05);//lrbt
    hPotAxis->Draw("AXIS");
    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        gPot_LL[uSce].Draw("C,same");
    }
    //ZeroLine->Draw("same");
    PotLegend->Draw("same");
    cPot_LL->SaveAs(TString::Format("%sPlotPot_LL_ver2.png",OutputDir.Data()));
    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        gPot_LL[uSce].SetLineWidth(gPot_LL[uSce].GetLineWidth()/2.5);
    }
    //gCkBaseline_LL.SetLineWidth(gCkBaseline_LL.GetLineWidth()/2.5);
    //ZeroLine->SetLineWidth(ZeroLine->GetLineWidth()/2.5);
    cPot_LL->SaveAs(TString::Format("%sPlotPot_LL_ver2.pdf",OutputDir.Data()));
    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        gPot_LL[uSce].SetLineWidth(gPot_LL[uSce].GetLineWidth()*2.5);
    }
    //gCkBaseline_LL.SetLineWidth(gCkBaseline_LL.GetLineWidth()*2.5);
    //ZeroLine->SetLineWidth(ZeroLine->GetLineWidth()*2.5);


    delete hAxis;
    delete hAxisZoom;
    delete UnitLine;
    delete myLegend;
    delete ptCorrected0;
    delete ptCorrected;

    delete hSigma_pL_InFuckingMeV;

    for(unsigned uRad=0; uRad<NumRad; uRad++) delete myPT[uRad];
    delete [] myPT;

    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        delete [] grCk[uSce];
        delete [] grCkLedni[uSce];
        delete [] grCkCorr[uSce];
        delete [] grCkLedniCorr[uSce];
    }
    delete [] grCk;
    delete [] grCkLedni;
    delete [] grCkCorr;
    delete [] grCkLedniCorr;
    delete [] gPot_LL;
    delete hPotAxis;
    delete PotLegend;
    //delete ZeroLine;
    delete cPot_LL;

    delete [] fBuffer;

    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        delete [] PotPars[uSce];
        delete [] ScatPars[uSce];
        delete [] ScatParErr[uSce];
        delete [] ScatParsMorita[uSce];
    }
    delete [] PotPars;
    delete [] ScatPars;
    delete [] ScatParErr;
    delete [] ScatParsMorita;
}


void CATSpaper::PlotCk_LL_GIF(){

    const double Plot_kMin_LL = 0;
    const double Plot_kMax_LL = 240;
    const double NumPlotPts_LL = 48;

    const unsigned NumRadPts = 64;
    const double radMin = 0.4;
    const double radMax = 2.0;
    const double radStep = (radMax-radMin)/double(NumRadPts-1);

    const TString SourceType = "Gauss";
    const bool ThetaDepSource = false;

    const unsigned NumPots = 4;
    //const TString PotNames[NumPots] = {"HKMYY","NF42","NF50","ND56"};
    const TString PotNames[NumPots] = {"NF42", "NF50","ND56","fss2"};

    const unsigned VideoLength = 8;
    const unsigned Delay = (100*VideoLength)/NumRadPts;

    //const double PotScatLen[NumPots] = {-0.575, 3.659, -0.772, -1.179};
    //const double PotEffRan[NumPots] = {6.45,0.975,4.271,4.656};

    double** PotPars = new double* [NumPots];
    for(unsigned uPot=0; uPot<NumPots; uPot++){
        PotPars[uPot] = new double [8];
        for(unsigned uEl=0; uEl<8; uEl++) PotPars[uPot][uEl]=0;

        if(PotNames[uPot]=="HKMYY")
        {PotPars[uPot][2]=-10.96; PotPars[uPot][3]=1.342; PotPars[uPot][4]=-141.75; PotPars[uPot][5]=0.777; PotPars[uPot][6]=2136.6; PotPars[uPot][7]=0.35;}

        if(PotNames[uPot]=="NF42")
        {PotPars[uPot][2]=-878.97; PotPars[uPot][3]=0.6; PotPars[uPot][4]=1048.58; PotPars[uPot][5]=0.45;}

        if(PotNames[uPot]=="NF44")
        {PotPars[uPot][2]=-1066.98; PotPars[uPot][3]=0.6; PotPars[uPot][4]=1646.65; PotPars[uPot][5]=0.45;}

        if(PotNames[uPot]=="NF50")
        {PotPars[uPot][2]=-2007.35; PotPars[uPot][3]=0.6; PotPars[uPot][4]=5678.97; PotPars[uPot][5]=0.45;}

        if(PotNames[uPot]=="ND56")
        {PotPars[uPot][2]=-144.26; PotPars[uPot][3]=1.0; PotPars[uPot][4]=1413.75; PotPars[uPot][5]=0.45;}

        if(PotNames[uPot]=="fss2")
        {PotPars[uPot][2]=-103.9; PotPars[uPot][3]=0.92; PotPars[uPot][4]=658.2; PotPars[uPot][5]=0.41;}

    }
    //PotPars[0][2]=-10.96; PotPars[0][3]=1.342; PotPars[0][4]=-141.75; PotPars[0][5]=0.777; PotPars[0][6]=2136.6; PotPars[0][7]=0.35;
    //PotPars[1][2]=-878.97; PotPars[1][3]=0.6; PotPars[1][4]=1048.58; PotPars[1][5]=0.45;
    //PotPars[2][2]=-2007.35; PotPars[2][3]=0.6; PotPars[2][4]=5678.97; PotPars[2][5]=0.45;
    //PotPars[3][2]=-144.26; PotPars[3][3]=1.0; PotPars[3][4]=1413.75; PotPars[3][5]=0.45;

    double LamLamSourcePars[4] = {0,0,0,radMin};

    CATSparameters cSourcePars_LL(CATSparameters::tSource,1,true);
    cSourcePars_LL.SetParameters(&LamLamSourcePars[3]);

    CATS Nala;
    Nala.SetEpsilonConv(1e-7);
    Nala.SetEpsilonProp(1e-7);
    Nala.SetThetaDependentSource(ThetaDepSource);
    Nala.SetUseAnalyticSource(true);
    if(Nala.GetThetaDependentSource()){
        if(SourceType=="Gauss") Nala.SetAnaSource(GaussSourceTheta, cSourcePars_LL);
        else if(SourceType=="Cauchy") Nala.SetAnaSource(CauchySourceTheta, cSourcePars_LL);
        else {printf("FUCK!!!\n"); return;}
    }
    else{
        if(SourceType=="Gauss") Nala.SetAnaSource(GaussSource, cSourcePars_LL);
        else if(SourceType=="Cauchy") Nala.SetAnaSource(CauchySource, cSourcePars_LL);
        else {printf("FUCK!!!\n"); return;}
    }

    Nala.SetExcludeFailedBins(false);
    Nala.SetMomBins(NumPlotPts_LL,Plot_kMin_LL,Plot_kMax_LL);

    Nala.SetNumChannels(2);
    Nala.SetNumPW(0,1);
    Nala.SetNumPW(1,1);
    Nala.SetSpin(0,0);
    Nala.SetSpin(1,1);
    Nala.SetChannelWeight(0, 0.25);
    Nala.SetChannelWeight(1, 0.75);
    Nala.SetQ1Q2(0);
    Nala.SetPdgId(3122, 3122);
    Nala.SetRedMass( 0.5*Mass_L );

    LamLamSourcePars[3] = radMin;

    Nala.KillTheCat();

    double Radius;
    for(unsigned uRad=0; uRad<NumRadPts; uRad++){
        Radius = radMin+double(uRad)*radStep;
printf("Radius=%f\n",Radius);
Nala.SetNotifications(CATS::nError);
        Nala.SetAnaSource(0,Radius);
        Nala.RemoveShortRangePotential();
        Nala.KillTheCat();

        TGraph gCkBaseline_LL;
        gCkBaseline_LL.Set(NumPlotPts_LL);
        gCkBaseline_LL.SetName(TString::Format("gCkBaseline_LL_uRad%u",uRad));
        gCkBaseline_LL.SetLineColor(kGray+1);
        gCkBaseline_LL.SetLineWidth(6);
        gCkBaseline_LL.SetLineStyle(7);
        gCkBaseline_LL.SetMarkerStyle(1);
        gCkBaseline_LL.SetMarkerSize(0);
        gCkBaseline_LL.SetMarkerColor(kWhite);
        gCkBaseline_LL.SetFillColor(kWhite);
        for(unsigned uBin=0; uBin<NumPlotPts_LL; uBin++){
            gCkBaseline_LL.SetPoint(uBin,Nala.GetMomentum(uBin),Nala.GetCorrFun(uBin));
        }

        TGraph* gCk_LL = new TGraph [NumPots];
        DLM_DtColor DlmCol;



        for(unsigned uPot=0; uPot<NumPots; uPot++){
            gCk_LL[uPot].Set(NumPlotPts_LL);
            gCk_LL[uPot].SetName(TString::Format("gCk_LL_%s_uRad%u",PotNames[uPot].Data(),uRad));
            gCk_LL[uPot].SetLineColor(DlmCol.GetColor(uPot));
            gCk_LL[uPot].SetLineWidth(12.-4.*double(uPot)/double(NumPots-1));
            gCk_LL[uPot].SetMarkerStyle(1);
            gCk_LL[uPot].SetMarkerSize(0);
            gCk_LL[uPot].SetMarkerColor(DlmCol.GetColor(uPot));
            gCk_LL[uPot].SetFillColor(kWhite);
            Nala.SetShortRangePotential(0,0,TripleGaussSum,PotPars[uPot]);
            Nala.KillTheCat();
            //printf("POT %u\n", uPot);
            for(unsigned uBin=0; uBin<NumPlotPts_LL; uBin++){
                //printf(" C(%.2f) = %.3f\n", Nala.GetMomentum(uBin), Nala.GetCorrFun(uBin));
                gCk_LL[uPot].SetPoint(uBin,Nala.GetMomentum(uBin),Nala.GetCorrFun(uBin));

            }
        }

        TH1F* hCkAxis = new TH1F("hCkAxis", "hCkAxis", NumPlotPts_LL,Plot_kMin_LL,Plot_kMax_LL);
        hCkAxis->SetStats(false);
        hCkAxis->SetTitle("");
        hCkAxis->GetXaxis()->SetLabelSize(0.065);
        hCkAxis->GetXaxis()->SetTitle("k (MeV)");
        hCkAxis->GetXaxis()->CenterTitle();
        hCkAxis->GetXaxis()->SetTitleOffset(1.15);
        hCkAxis->GetXaxis()->SetLabelOffset(0.02);
        hCkAxis->GetXaxis()->SetTitleSize(0.075);
        hCkAxis->GetYaxis()->SetLabelSize(0.065);
        hCkAxis->GetYaxis()->SetTitle("C(k)");
        hCkAxis->GetYaxis()->CenterTitle();
        hCkAxis->GetYaxis()->SetTitleOffset(0.75);
        hCkAxis->GetYaxis()->SetTitleSize(0.075);
        hCkAxis->GetYaxis()->SetRangeUser(0.0, 3);
        hCkAxis->GetYaxis()->SetLimits(0.0, 3);
        hCkAxis->GetXaxis()->SetNdivisions(506);

        TLegend* CkLegend = new TLegend(0.7,0.7,0.975,0.98);//lbrt
        CkLegend->SetName(TString::Format("CkLegend_uRad%u",uRad));
        CkLegend->SetTextSize(0.05);
        for(unsigned uPot=0; uPot<NumPots; uPot++){
            CkLegend->AddEntry(&gCk_LL[uPot], TString::Format("%s", PotNames[uPot].Data()));
        }
        CkLegend->AddEntry(&gCkBaseline_LL, "W/o strong inter.");

        TPaveText* myPT = new TPaveText(0.5,0.8,0.7,0.98, "blNDC");//lbrt
        myPT->SetName(TString::Format("ptSource_uRad%u",uRad));
        myPT->SetBorderSize(1);
        myPT->SetTextSize(0.05);
        myPT->SetFillColor(kWhite);
        myPT->SetTextFont(22);
        myPT->AddText("Gauss source");
        myPT->AddText(TString::Format("R = %.2f fm", Radius));

        TF1* UnitLine = new TF1(TString::Format("UnitLine_uRad%u",uRad), "1", Plot_kMin_LL, Plot_kMax_LL);
        UnitLine->SetLineColor(kGray);
        UnitLine->SetLineWidth(6);

        TCanvas* cCk_LL = new TCanvas(TString::Format("cCk_LL_uRad%u",uRad), TString::Format("cCk_LL_uRad%u",uRad), 1);
        cCk_LL->cd(0); cCk_LL->SetCanvasSize(1920, 1080); cCk_LL->SetMargin(0.12,0.025,0.18,0.05);//lrbt
        hCkAxis->Draw("AXIS");
        gCkBaseline_LL.Draw("same");
        for(unsigned uPot=0; uPot<NumPots; uPot++){
            gCk_LL[uPot].Draw("C,same");
        }
        //UnitLine->Draw("same");
        CkLegend->Draw("same");
        myPT->Draw("same");

        if(!uRad) cCk_LL->SaveAs(TString::Format("%sAnimationCk_LL.gif",OutputDir.Data()));
        cCk_LL->SaveAs(TString::Format("%sAnimationCk_LL.gif+%u",OutputDir.Data(),Delay));

        delete [] gCk_LL;
        delete hCkAxis;
        delete UnitLine;
        delete CkLegend;
        delete cCk_LL;
        delete myPT;
    }

    for(unsigned uPot=0; uPot<NumPots; uPot++){
        delete [] PotPars[uPot];
    }
    delete [] PotPars;

}


void CATSpaper::PlotCk_KminProton(){

    //0 = Haidenbauer (Juelich)
    //1 = Tetsuo (Kyoto) I=0 + I=1
    //2 = Tetsuo Avg in Isospin
    //3 = Tetsuo (Kyoto) I=0 + I=1, exact wave function (taken with a square)
    //4 = Tetsuo (Kyoto) I=0 + I=1, exact wave function (taken with a square + the dirty trick for sigh flipping)
    //5 = Tetsuo as above, but real+imag
    const int WHICHPOT = 5;

    const double Plot_kMin = 5;
    const double Plot_kMax = 305;
    const unsigned NumPlotPts = 30;
    //CATShisto<double> CatHist(NumPlotPts,Plot_kMin,Plot_kMax);

    unsigned NumRad = 3;
    double Radius[NumRad] = {2.0,1.2,0.85};
    //double Radius[NumRad] = {5,3,2,1};
//const double Radius[NumRad] = {3.0,4.0,5.0};

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
    ScenarioName[0] = "Coulomb only";
    ScenarioName[1] = "Coulomb + Strong";
    //ScenarioName[2] = "Reid SC as in CRAB";

    TString ScenarioShortName[NumScenarios];
    ScenarioShortName[0] = "NoSI";
    ScenarioShortName[1] = "WithSI";

    double SourcePars_Kp[4] = {0,0,0,Radius[0]};

    CATSparameters cSourcePars_Kp(CATSparameters::tSource,1,true);
    cSourcePars_Kp.SetParameters(&SourcePars_Kp[3]);

    CATS MyCat;
    MyCat.SetUseAnalyticSource(true);
    MyCat.SetAnaSource(GaussSource, cSourcePars_Kp);

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
        MyCat.SetNumChannels(2);
        MyCat.SetNumPW(0,1);
        MyCat.SetNumPW(1,1);
        //change to double?
        MyCat.SetSpin(0,1);
        MyCat.SetSpin(1,1);//the no SI channel
    }

    MyCat.SetQ1Q2(-1);
    MyCat.SetPdgId(2212, 222);//check the pdg id of the kaon
    //const double Mass_p = 938.272;
    const double Mass_Km = 493.677;
    MyCat.SetRedMass( (Mass_p*Mass_Km)/(Mass_p+Mass_Km) );

    //in bins of Scenario/Radius
    TGraph** grCk;
    grCk = new TGraph* [NumScenarios];

    complex<double>**** WaveFunctionU; double*** PhaseShifts; double* RadBins; unsigned NumRadBins=0;
    int TetsPotType = int(WHICHPOT)-3;
    if(TetsPotType<0)TetsPotType=0;
    if(WHICHPOT==0)InitHaidenbauerKaonMinus("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Haidenbauer/KaonMinus10MeV/",MyCat,&WaveFunctionU,&PhaseShifts,&RadBins,NumRadBins);
    else InitTetsuoKaonMinus("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Tetsuo/KminusProton/",MyCat,&WaveFunctionU,&PhaseShifts,&RadBins,NumRadBins, TetsPotType);
    double Kp_potPars[2];

//IsoSpin, const int& t2p1, const int& t2p2,
//                     const int& Spin, const int& AngMom, const int& TotMom

    //double pKmPotParsI0[10]={0,0,pKm_Tetsuo,0,0,1,1,1,0,1};
    //double pKmPotParsI1[10]={0,0,pKm_Tetsuo,0,1,1,1,1,0,1};
    double pKmPotParsI0[10]={0,0,pKm_Tetsuo,0,0,1,1,1,0,1};
    double pKmPotParsI1[10]={0,0,pKm_Tetsuo,0,1,1,1,1,0,1};
    double pKmPotParsIavg[10]={0,0,pKm_Tetsuo,2,0,1,1,1,0,1};

    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        grCk[uSce] = new TGraph[NumRad];
        switch(uSce){
        case 0 :
            if(WHICHPOT==0){
                MyCat.SetChannelWeight(0, 0.);
                MyCat.SetChannelWeight(1, 0.);
                MyCat.SetChannelWeight(2, 0.);
                MyCat.SetChannelWeight(3, 0.);
                MyCat.SetChannelWeight(4, 1.);
                MyCat.SetShortRangePotential(4,0,ZeroPotential,Kp_potPars);
            }
            else if(WHICHPOT==1){
                MyCat.SetChannelWeight(0, 0.);
                MyCat.SetChannelWeight(1, 0.);
                MyCat.SetChannelWeight(2, 1.);
                MyCat.SetShortRangePotential(2,0,ZeroPotential,Kp_potPars);
            }
            else if(WHICHPOT==2){
                MyCat.SetChannelWeight(0, 0.);
                MyCat.SetChannelWeight(1, 1.);
                MyCat.SetShortRangePotential(1,0,ZeroPotential,Kp_potPars);
            }
            else if(WHICHPOT==3){
                MyCat.SetChannelWeight(0, 0.);
                MyCat.SetChannelWeight(1, 1.);
                MyCat.SetShortRangePotential(1,0,ZeroPotential,Kp_potPars);
            }
            else if(WHICHPOT==4){
                MyCat.SetChannelWeight(0, 0.);
                MyCat.SetChannelWeight(1, 1.);
                MyCat.SetShortRangePotential(1,0,ZeroPotential,Kp_potPars);
            }
            else if(WHICHPOT==5){
                MyCat.SetChannelWeight(0, 0.);
                MyCat.SetChannelWeight(1, 1.);
                MyCat.SetShortRangePotential(1,0,ZeroPotential,Kp_potPars);
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
                    MyCat.UseExternalWaveFunction(uMomBin,0,0,WaveFunctionU[uMomBin][0][0],NumRadBins,RadBins,PhaseShifts[uMomBin][0][0]);
                    MyCat.UseExternalWaveFunction(uMomBin,0,1,WaveFunctionU[uMomBin][0][1],NumRadBins,RadBins,PhaseShifts[uMomBin][0][1]);
                    MyCat.UseExternalWaveFunction(uMomBin,1,0,WaveFunctionU[uMomBin][1][0],NumRadBins,RadBins,PhaseShifts[uMomBin][1][0]);
                    MyCat.UseExternalWaveFunction(uMomBin,1,1,WaveFunctionU[uMomBin][1][1],NumRadBins,RadBins,PhaseShifts[uMomBin][1][1]);
                    MyCat.UseExternalWaveFunction(uMomBin,2,0,WaveFunctionU[uMomBin][2][0],NumRadBins,RadBins,PhaseShifts[uMomBin][2][0]);
                    MyCat.UseExternalWaveFunction(uMomBin,2,1,WaveFunctionU[uMomBin][2][1],NumRadBins,RadBins,PhaseShifts[uMomBin][2][1]);
                    MyCat.UseExternalWaveFunction(uMomBin,3,0,WaveFunctionU[uMomBin][3][0],NumRadBins,RadBins,PhaseShifts[uMomBin][3][0]);
                    MyCat.UseExternalWaveFunction(uMomBin,3,1,WaveFunctionU[uMomBin][3][1],NumRadBins,RadBins,PhaseShifts[uMomBin][3][1]);
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
                    MyCat.UseExternalWaveFunction(uMomBin,0,0,WaveFunctionU[uMomBin][0][0],NumRadBins,RadBins,PhaseShifts[uMomBin][0][0]);
                }
            }
            else if(WHICHPOT==4){
                MyCat.SetChannelWeight(0, 1.);
                MyCat.SetChannelWeight(1, 0.);
                for(unsigned uMomBin=0; uMomBin<NumPlotPts; uMomBin++){
                    MyCat.UseExternalWaveFunction(uMomBin,0,0,WaveFunctionU[uMomBin][0][0],NumRadBins,RadBins,PhaseShifts[uMomBin][0][0]);
                }
            }
            else if(WHICHPOT==5){
                MyCat.SetChannelWeight(0, 1.);
                MyCat.SetChannelWeight(1, 0.);
                for(unsigned uMomBin=0; uMomBin<NumPlotPts; uMomBin++){
                    MyCat.UseExternalWaveFunction(uMomBin,0,0,WaveFunctionU[uMomBin][0][0],NumRadBins,RadBins,PhaseShifts[uMomBin][0][0]);
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
            grCk[uSce][uRad].SetLineWidth(8.);
            grCk[uSce][uRad].SetFillColor(kWhite);

            for(unsigned uMomBin=0; uMomBin<NumPlotPts; uMomBin++){
                grCk[uSce][uRad].SetPoint(uMomBin, MyCat.GetMomentum(uMomBin), MyCat.GetCorrFun(uMomBin));
            }
            if(uSce==0){
                //grCk[uSce][uRad].SetLineColor(kGray+1);
                grCk[uSce][uRad].SetLineStyle(2);
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

    const double Yoffset = 0.9;
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

    hAxis->GetYaxis()->SetRangeUser(0.75, 3.6);
    hAxis->GetYaxis()->SetLimits(0.75, 3.6);
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
    hAxisZoom->GetYaxis()->SetTitle("C(k)");
    hAxisZoom->GetYaxis()->CenterTitle();
    hAxisZoom->GetYaxis()->SetTitleOffset(Yoffset/1.5);
    hAxisZoom->GetYaxis()->SetTitleSize(0.075);

    hAxisZoom->GetXaxis()->SetNdivisions(504);
    hAxisZoom->GetYaxis()->SetRangeUser(0.95, 1.05);
    hAxisZoom->GetYaxis()->SetLimits(0.95, 1.05);
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
        for(unsigned uSce=0; uSce<NumScenarios; uSce++) grCk[uSce][uRad].Draw("same,C");

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
        }
    }
    hAxis->SetLineWidth(hAxis->GetLineWidth()/2.5);
    UnitLine->SetLineWidth(UnitLine->GetLineWidth()/2.5);
    gStyle->SetLineWidth(1);
    DrawBoard_Ck_pp.GetCanvas()->SaveAs(TString::Format("%sPlotCk_KmProton.pdf",OutputDir.Data()));

    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        for(unsigned uRad=0; uRad<NumRad; uRad++){
            grCk[uSce][uRad].SetLineWidth(grCk[uSce][uRad].GetLineWidth()*2.5);
        }
    }
    hAxis->SetLineWidth(hAxis->GetLineWidth()*2.5);
    UnitLine->SetLineWidth(UnitLine->GetLineWidth()*2.5);

    TFile* RootFile = new TFile(TString::Format("%sKmProton.root",OutputDir.Data()),"recreate");
    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        for(unsigned uRad=0; uRad<NumRad; uRad++){
            grCk[uSce][uRad].Write();
        }
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
    }
    delete [] grCk;
    delete RootFile;

    CleanHaidenbauer(MyCat, &WaveFunctionU,&PhaseShifts,&RadBins);
}

void CATSpaper::PlotCk_KplusProton(){
    const double Plot_kMin = 5;
    const double Plot_kMax = 305;
    const unsigned NumPlotPts = 30;
    //CATShisto<double> CatHist(NumPlotPts,Plot_kMin,Plot_kMax);

    //const unsigned NumRad = 3;
    //const double Radius[NumRad] = {2.0,1.2,0.85};

    const unsigned NumRad = 6;
    const double Radius[NumRad] = {1.113,1.144,1.162,1.188,1.213,1.232};

    //double CkVal;
    //double Momentum;
    DLM_DtColor DlmCol;

    const unsigned NumScenarios = 2;
    TString ScenarioName[NumScenarios];
    ScenarioName[0] = "Coulomb only";
    ScenarioName[1] = "Coulomb + Strong";
    //ScenarioName[2] = "Reid SC as in CRAB";

    TString ScenarioShortName[NumScenarios];
    ScenarioShortName[0] = "NoSI";
    ScenarioShortName[1] = "WithSI";

    double SourcePars_Kp[4] = {0,0,0,Radius[0]};

    CATSparameters cSourcePars_Kp(CATSparameters::tSource,1,true);
    cSourcePars_Kp.SetParameters(&SourcePars_Kp[3]);

    CATS MyCat;
    MyCat.SetUseAnalyticSource(true);
    MyCat.SetAnaSource(GaussSource, cSourcePars_Kp);

    MyCat.SetExcludeFailedBins(false);
    MyCat.SetMomBins(NumPlotPts,Plot_kMin,Plot_kMax);

    MyCat.SetNumChannels(3);
    MyCat.SetNumPW(0,2);
    MyCat.SetNumPW(1,2);
    MyCat.SetNumPW(2,2);
    //change to double?
    MyCat.SetSpin(0,1);//s11+p11
    MyCat.SetSpin(1,1);//s11+p13
    MyCat.SetSpin(2,1);//n/a


    MyCat.SetQ1Q2(1);
    MyCat.SetPdgId(2212, 222);//check the pdg id of the kaon
    //const double Mass_p = 938.272;
    const double Mass_Km = 493.677;
    MyCat.SetRedMass( (Mass_p*Mass_Km)/(Mass_p+Mass_Km) );

    //in bins of Scenario/Radius
    TGraph** grCk;
    grCk = new TGraph* [NumScenarios];

    complex<double>**** WaveFunctionU; double*** PhaseShifts; double* RadBins; unsigned NumRadBins=0;
    InitHaidenbauerKaonPlus("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Haidenbauer/KaonPlus10MeV/",MyCat,&WaveFunctionU,&PhaseShifts,&RadBins,NumRadBins);
    double Kp_potPars[2];

    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        grCk[uSce] = new TGraph[NumRad];
        switch(uSce){
        case 0 :
            MyCat.SetChannelWeight(0, 0.);
            MyCat.SetChannelWeight(1, 0.);
            MyCat.SetChannelWeight(2, 1.);
            MyCat.SetShortRangePotential(4,0,ZeroPotential,Kp_potPars);
            break;
        case 1 :
            MyCat.SetChannelWeight(0, 1./3.);
            MyCat.SetChannelWeight(1, 2./3.);
            MyCat.SetChannelWeight(2, 0.);

            for(unsigned uMomBin=0; uMomBin<NumPlotPts; uMomBin++){
                MyCat.UseExternalWaveFunction(uMomBin,0,0,WaveFunctionU[uMomBin][0][0],NumRadBins,RadBins,PhaseShifts[uMomBin][0][0]);
                MyCat.UseExternalWaveFunction(uMomBin,0,1,WaveFunctionU[uMomBin][0][1],NumRadBins,RadBins,PhaseShifts[uMomBin][0][1]);
                MyCat.UseExternalWaveFunction(uMomBin,1,0,WaveFunctionU[uMomBin][1][0],NumRadBins,RadBins,PhaseShifts[uMomBin][1][0]);
                MyCat.UseExternalWaveFunction(uMomBin,1,1,WaveFunctionU[uMomBin][1][1],NumRadBins,RadBins,PhaseShifts[uMomBin][1][1]);
            }

            break;
        default :
            printf("Default K+ fuck!\n");
            break;
        }

        MyCat.KillTheCat(CATS::kAllChanged);

        for(unsigned uRad=0; uRad<NumRad; uRad++){
            MyCat.SetAnaSource(0, Radius[uRad]);
            MyCat.KillTheCat();

            grCk[uSce][uRad].Set(NumPlotPts);
            //grCk[uSce][uRad].SetName(TString::Format("grCk_%u_%u",uSce,uRad));
            grCk[uSce][uRad].SetName(TString::Format("grCk_%u_%.2f",uSce,Radius[uRad]));
            grCk[uSce][uRad].SetMarkerColor(DlmCol.GetColor(1));
            grCk[uSce][uRad].SetMarkerStyle(20);
            grCk[uSce][uRad].SetMarkerSize(0);
            grCk[uSce][uRad].SetLineColor(DlmCol.GetColor(1+uSce));
            //the Coulomb only case
            //grCk[uSce][uRad].SetLineWidth(8.-4.*double(uSce)/double(NumScenarios-1));
            grCk[uSce][uRad].SetLineWidth(8.);
            grCk[uSce][uRad].SetFillColor(kWhite);

            for(unsigned uMomBin=0; uMomBin<NumPlotPts; uMomBin++){
                grCk[uSce][uRad].SetPoint(uMomBin, MyCat.GetMomentum(uMomBin), MyCat.GetCorrFun(uMomBin));
            }
            if(uSce==0){
                //grCk[uSce][uRad].SetLineColor(kGray+1);
                grCk[uSce][uRad].SetLineStyle(2);
            }
        }
    }


    char* fBuffer = new char [512];
    FILE * pFile;
    pFile = fopen (TString::Format("%sPlotCk_KpProton.txt",OutputDir.Data()),"w");
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

    const double Yoffset = 1.25;
    TH1F* hAxis = new TH1F("hAxis", "hAxis", NumPlotPts,Plot_kMin,Plot_kMax);
    hAxis->SetStats(false);
    hAxis->SetTitle("");
    hAxis->GetXaxis()->SetLabelSize(0.065);
    hAxis->GetXaxis()->SetTitle("k (MeV)");
    hAxis->GetXaxis()->CenterTitle();
    hAxis->GetXaxis()->SetTitleOffset(1.05);
    hAxis->GetXaxis()->SetLabelOffset(0.02);
    hAxis->GetXaxis()->SetTitleSize(0.075);
    hAxis->GetYaxis()->SetLabelSize(0.065);
    hAxis->GetYaxis()->SetTitle("C(k)");
    hAxis->GetYaxis()->CenterTitle();
    hAxis->GetYaxis()->SetTitleOffset(Yoffset);
    hAxis->GetYaxis()->SetTitleSize(0.075);

    hAxis->GetYaxis()->SetRangeUser(0.0, 1.3);
    hAxis->GetYaxis()->SetLimits(0.0, 1.3);

    hAxis->GetXaxis()->SetNdivisions(504);

    TH1F* hAxisZoom = new TH1F("hAxisZoom", "hAxisZoom", NumPlotPts,Plot_kMin,Plot_kMax);
    hAxisZoom->SetStats(false);
    hAxisZoom->SetTitle("");
    hAxisZoom->GetXaxis()->SetLabelSize(0.065);
    hAxisZoom->GetXaxis()->SetTitle("k (MeV)");
    hAxisZoom->GetXaxis()->CenterTitle();
    hAxisZoom->GetXaxis()->SetTitleOffset(1.3);
    hAxisZoom->GetXaxis()->SetLabelOffset(0.02);
    hAxisZoom->GetXaxis()->SetTitleSize(0.075);
    hAxisZoom->GetYaxis()->SetLabelSize(0.065);
    hAxisZoom->GetYaxis()->SetTitle("C(k)");
    hAxisZoom->GetYaxis()->CenterTitle();
    hAxisZoom->GetYaxis()->SetTitleOffset(Yoffset/1.5);
    hAxisZoom->GetYaxis()->SetTitleSize(0.075);

    hAxisZoom->GetXaxis()->SetNdivisions(504);
    hAxisZoom->GetYaxis()->SetRangeUser(0.85, 1.15);
    hAxisZoom->GetYaxis()->SetLimits(0.85, 1.15);
    hAxisZoom->GetYaxis()->SetNdivisions(603);

    TPaveText** myPT = new TPaveText* [NumRad];
    for(unsigned uRad=0; uRad<NumRad; uRad++){
        myPT[uRad] = new TPaveText(0.5,0.885,0.95,0.975, "blNDC");//lbrt
        myPT[uRad]->SetName("PT12");
        myPT[uRad]->SetBorderSize(1);
        myPT[uRad]->SetTextSize(0.08/1.4);
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
        DrawBoard_Ck_pp.AddSubPadTL(CurrentPos,CurrentPos+PadLength,0.08,1);//lrtb
        //DrawBoard_Ck_pp.AddSubPadTL(CurrentPos,CurrentPos+PadLength,0.08,0.63);//lrtb
        //DrawBoard_Ck_pp.AddSubPadTL(CurrentPos,CurrentPos+PadLength,0.63,1);//lrtb
        CurrentPos+=PadLength;
//printf("uRad=%u; %.2f <--> %.2f\n", uRad, );
    }
    DrawBoard_Ck_pp.AddSubPadTL(0,1,0,0.075);//lrtb

    for(unsigned uRad=0; uRad<NumRad; uRad++){
        DrawBoard_Ck_pp.SetMargin(uRad, 0.07, 0.01, 0.15, 0.04);//lrbt
        //DrawBoard_Ck_pp.SetMargin(uRad*2+1, 0.07, 0.01, 0.14, 0);//lrbt
    }
    DrawBoard_Ck_pp.SetMargin(NumRad, 0, 0, 0, 0);//lrbt

    for(unsigned uRad=0; uRad<NumRad; uRad++){
        DrawBoard_Ck_pp.GetCanvas()->SetName(TString::Format("DrawBoard_Ck_KmProton_%u",uRad));
        DrawBoard_Ck_pp.cd(uRad);
        hAxis->Draw("AXIS");
        UnitLine->Draw("same");
        for(unsigned uSce=0; uSce<NumScenarios; uSce++) grCk[uSce][uRad].Draw("same,C");
        for(unsigned uSce=0; uSce<NumScenarios; uSce++) myPT[uRad]->Draw("same");

        DrawBoard_Ck_pp.SetLabelSize(uRad, hAxis->GetXaxis(), 17);
        DrawBoard_Ck_pp.SetLabelSize(uRad, hAxis->GetYaxis(), 17);
        DrawBoard_Ck_pp.SetTitleSize(uRad, hAxis->GetXaxis(), 17);
        DrawBoard_Ck_pp.SetTitleSize(uRad, hAxis->GetYaxis(), 17);
    }
    //for(unsigned uRad=0; uRad<NumRad; uRad++){
    //    DrawBoard_Ck_pp.GetCanvas()->SetName(TString::Format("DrawBoard_Ck_KmProton_%u",uRad));
    //    DrawBoard_Ck_pp.cd(uRad*2+1);
    //    hAxisZoom->Draw("AXIS");
    //    UnitLine->Draw("same");
    //    for(unsigned uSce=0; uSce<NumScenarios; uSce++) grCk[uSce][uRad].Draw("same,C");

    //    DrawBoard_Ck_pp.SetLabelSize(uRad*2+1, hAxisZoom->GetXaxis(), 17);
    //    DrawBoard_Ck_pp.SetLabelSize(uRad*2+1, hAxisZoom->GetYaxis(), 17);
    //    DrawBoard_Ck_pp.SetTitleSize(uRad*2+1, hAxisZoom->GetXaxis(), 17);
    //    DrawBoard_Ck_pp.SetTitleSize(uRad*2+1, hAxisZoom->GetYaxis(), 17);
    //}

    TLegend* myLegend = new TLegend(0.12,0,0.95,0.96);//lbrt
    myLegend->SetName("myLegend");
    myLegend->SetTextSize(0.65);
    myLegend->SetNColumns(NumScenarios);
    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        myLegend->AddEntry(&grCk[uSce][0],  ScenarioName[uSce]);
    }

    DrawBoard_Ck_pp.cd(NumRad);
    myLegend->Draw();

    gStyle->SetLineWidth(2.5);
    DrawBoard_Ck_pp.GetCanvas()->SaveAs(TString::Format("%sPlotCk_KpProton.png",OutputDir.Data()));
    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        for(unsigned uRad=0; uRad<NumRad; uRad++){
            grCk[uSce][uRad].SetLineWidth(grCk[uSce][uRad].GetLineWidth()/2.5);
        }
    }
    hAxis->SetLineWidth(hAxis->GetLineWidth()/2.5);
    UnitLine->SetLineWidth(UnitLine->GetLineWidth()/2.5);
    gStyle->SetLineWidth(1);
    DrawBoard_Ck_pp.GetCanvas()->SaveAs(TString::Format("%sPlotCk_KpProton.pdf",OutputDir.Data()));

    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        for(unsigned uRad=0; uRad<NumRad; uRad++){
            grCk[uSce][uRad].SetLineWidth(grCk[uSce][uRad].GetLineWidth()*2.5);
        }
    }
    hAxis->SetLineWidth(hAxis->GetLineWidth()*2.5);
    UnitLine->SetLineWidth(UnitLine->GetLineWidth()*2.5);

    TFile* RootFile = new TFile(TString::Format("%sKpProton.root",OutputDir.Data()),"recreate");
    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        for(unsigned uRad=0; uRad<NumRad; uRad++){
            grCk[uSce][uRad].Write();
        }
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
    }
    delete [] grCk;
    delete RootFile;

    CleanHaidenbauer(MyCat, &WaveFunctionU,&PhaseShifts,&RadBins);
}

//WhichChannel
//0= I0 S0
//1= I0 S1
//2= I1 S0
//3= I1 S1
//10=AVERAGED
//anything else is the full potential

//WhichAxis
//0 = standard
//1 = zoomed
//2 = very zoomed
void CATSpaper::PlotCk_XiMinusProton(const int WhichChannel, const int WhichAxis){
    const double Plot_kMin = 0;
    const double Plot_kMax = 360;
    const unsigned NumPlotPts = 90;
    //CATShisto<double> CatHist(NumPlotPts,Plot_kMin,Plot_kMax);

    const unsigned NumRad = 3;
    const double Radius[NumRad] = {2.0,1.2,0.85};

    //double CkVal;
    //double Momentum;
    DLM_DtColor DlmCol;

    const unsigned NumScenarios = 2;
    TString ScenarioName[NumScenarios];
    ScenarioName[0] = "Coulomb only";
    if(WhichChannel==0) ScenarioName[1] = "Coulomb + Strong(I=0,S=0)";
    else if(WhichChannel==1) ScenarioName[1] = "Coulomb + Strong(I=0,S=1)";
    else if(WhichChannel==2) ScenarioName[1] = "Coulomb + Strong(I=1,S=0)";
    else if(WhichChannel==3) ScenarioName[1] = "Coulomb + Strong(I=1,S=1)";
    else if(WhichChannel==10) ScenarioName[1] = "Coulomb + Strong(I,S avg)";
    else ScenarioName[1] = "Coulomb + Strong";
    //ScenarioName[2] = "Reid SC as in CRAB";

    TString ScenarioShortName[NumScenarios];
    ScenarioShortName[0] = "NoSI";
    if(WhichChannel==0) ScenarioShortName[1] = "I=0,S=0";
    else if(WhichChannel==1) ScenarioShortName[1] = "I=0,S=1";
    else if(WhichChannel==2) ScenarioShortName[1] = "I=1,S=0";
    else if(WhichChannel==3) ScenarioShortName[1] = "I=1,S=1";
    else if(WhichChannel==10) ScenarioShortName[1] = "I,S avg";
    else ScenarioShortName[1] = "Full";

    double SourcePars_Xip[4] = {0,0,0,Radius[0]};

    CATSparameters cSourcePars_Xip(CATSparameters::tSource,1,true);
    cSourcePars_Xip.SetParameters(&SourcePars_Xip[3]);

    CATS MyCat;
    MyCat.SetUseAnalyticSource(true);
    MyCat.SetAnaSource(GaussSource, cSourcePars_Xip);

    MyCat.SetExcludeFailedBins(false);
    MyCat.SetMomBins(NumPlotPts,Plot_kMin,Plot_kMax);

    MyCat.SetNumChannels(4);
    MyCat.SetNumPW(0,1);
    MyCat.SetNumPW(1,1);
    MyCat.SetNumPW(2,1);
    MyCat.SetNumPW(3,1);
    MyCat.SetSpin(0,0);//I=0; S=0
    MyCat.SetSpin(1,1);//I=0; S=1
    MyCat.SetSpin(2,0);//I=1; S=0
    MyCat.SetSpin(3,1);//I=1; S=1
    if(WhichChannel==0){
        MyCat.SetChannelWeight(0, 1.);
        MyCat.SetChannelWeight(1, 0.);
        MyCat.SetChannelWeight(2, 0.);
        MyCat.SetChannelWeight(3, 0.);
    }
    else if(WhichChannel==1){
        MyCat.SetChannelWeight(0, 0.);
        MyCat.SetChannelWeight(1, 1.);
        MyCat.SetChannelWeight(2, 0.);
        MyCat.SetChannelWeight(3, 0.);
    }
    else if(WhichChannel==2){
        MyCat.SetChannelWeight(0, 0.);
        MyCat.SetChannelWeight(1, 0.);
        MyCat.SetChannelWeight(2, 1.);
        MyCat.SetChannelWeight(3, 0.);
    }
    else if(WhichChannel==3){
        MyCat.SetChannelWeight(0, 0.);
        MyCat.SetChannelWeight(1, 0.);
        MyCat.SetChannelWeight(2, 0.);
        MyCat.SetChannelWeight(3, 1.);
    }
    else if(WhichChannel==10){
        MyCat.SetChannelWeight(0, 1.);
        MyCat.SetChannelWeight(1, 0.);
        MyCat.SetChannelWeight(2, 0.);
        MyCat.SetChannelWeight(3, 0.);
    }
    else{
        MyCat.SetChannelWeight(0, 1./8.);
        MyCat.SetChannelWeight(1, 3./8.);
        MyCat.SetChannelWeight(2, 1./8.);
        MyCat.SetChannelWeight(3, 3./8.);
    }

    MyCat.SetQ1Q2(-1);
    MyCat.SetPdgId(2212, 3312);
    //const double Mass_p = 938.272;
    const double Mass_Xim = 1321.7;
    MyCat.SetRedMass( (Mass_p*Mass_Xim)/(Mass_p+Mass_Xim) );

    MyCat.SetMaxRad(64);
    MyCat.SetMaxRho(32);

    MyCat.SetEpsilonConv(1e-8);
    MyCat.SetEpsilonProp(1e-8);

    //in bins of Scenario/Radius
    TGraph** grCk;
    grCk = new TGraph* [NumScenarios];

    double pXimPotParsISavg[10]={0,0,pXim_LatticeAvg,12,0,-1,1,0,0,0};
    double pXimPotParsI01S0[10]={0,0,pXim_Lattice,12,0,-1,1,0,0,0};
    double pXimPotParsI03S1[10]={0,0,pXim_Lattice,12,0,-1,1,1,0,1};
    double pXimPotParsI11S0[10]={0,0,pXim_Lattice,6,1,1,1,0,0,0};
    double pXimPotParsI13S1[10]={0,0,pXim_Lattice,6,1,1,1,1,0,1};


    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        grCk[uSce] = new TGraph[NumRad];
        switch(uSce){
        case 0 :
            MyCat.SetShortRangePotential(0,0,ZeroPotential,pXimPotParsI01S0);
            MyCat.SetShortRangePotential(1,0,ZeroPotential,pXimPotParsI03S1);
            MyCat.SetShortRangePotential(2,0,ZeroPotential,pXimPotParsI11S0);
            MyCat.SetShortRangePotential(3,0,ZeroPotential,pXimPotParsI13S1);
            break;
        case 1 :
            if(WhichChannel==10){
                MyCat.SetShortRangePotential(0,0,fDlmPot,pXimPotParsISavg);
                MyCat.SetShortRangePotential(1,0,ZeroPotential,pXimPotParsI03S1);
                MyCat.SetShortRangePotential(2,0,ZeroPotential,pXimPotParsI11S0);
                MyCat.SetShortRangePotential(3,0,ZeroPotential,pXimPotParsI13S1);
            }
            else{
                MyCat.SetShortRangePotential(0,0,fDlmPot,pXimPotParsI01S0);
                MyCat.SetShortRangePotential(1,0,fDlmPot,pXimPotParsI03S1);
                MyCat.SetShortRangePotential(2,0,fDlmPot,pXimPotParsI11S0);
                MyCat.SetShortRangePotential(3,0,fDlmPot,pXimPotParsI13S1);
            }
            //MyCat.SetShortRangePotential(2,0,ZeroPotential,pXimPotParsI11S0);
            //MyCat.SetShortRangePotential(3,0,ZeroPotential,pXimPotParsI13S1);
            break;
        default :
            printf("Default Xi- fuck!\n");
            break;
        }
        MyCat.KillTheCat(CATS::kAllChanged);

//double TEMP_X;
//double TEMP_Y;

        for(unsigned uRad=0; uRad<NumRad; uRad++){
            MyCat.SetAnaSource(0, Radius[uRad]);
            MyCat.KillTheCat();

            grCk[uSce][uRad].Set(NumPlotPts);
            grCk[uSce][uRad].SetName(TString::Format("grCk_%u_%u",uSce,uRad));
            grCk[uSce][uRad].SetMarkerColor(DlmCol.GetColor(1));
            grCk[uSce][uRad].SetMarkerStyle(20);
            grCk[uSce][uRad].SetMarkerSize(0);
            grCk[uSce][uRad].SetLineColor(DlmCol.GetColor(1+uSce));
            //the Coulomb only case
            //grCk[uSce][uRad].SetLineWidth(8.-4.*double(uSce)/double(NumScenarios-1));
            grCk[uSce][uRad].SetLineWidth(8.);
            grCk[uSce][uRad].SetFillColor(kWhite);

            for(unsigned uMomBin=0; uMomBin<NumPlotPts; uMomBin++){
//grCk[0][uRad].GetPoint(uMomBin, TEMP_X,TEMP_Y);
//if(uSce==0) TEMP_Y=MyCat.GetCorrFun(uMomBin);
//grCk[uSce][uRad].SetPoint(uMomBin, MyCat.GetMomentum(uMomBin), MyCat.GetCorrFun(uMomBin)/TEMP_Y);
                grCk[uSce][uRad].SetPoint(uMomBin, MyCat.GetMomentum(uMomBin), MyCat.GetCorrFun(uMomBin));
            }
            if(uSce==0){
                //grCk[uSce][uRad].SetLineColor(kGray+1);
                grCk[uSce][uRad].SetLineStyle(2);
            }
        }
    }


    char* fBuffer = new char [512];
    FILE * pFile;
    pFile = fopen (TString::Format("%sPlotCk_XimProton.txt",OutputDir.Data()),"w");
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

    const double Yoffset = 1.25;
    TH1F* hAxis = new TH1F("hAxis", "hAxis", NumPlotPts,Plot_kMin,Plot_kMax);
    hAxis->SetStats(false);
    hAxis->SetTitle("");
    hAxis->GetXaxis()->SetLabelSize(0.065);
    hAxis->GetXaxis()->SetTitle("k (MeV)");
    hAxis->GetXaxis()->CenterTitle();
    hAxis->GetXaxis()->SetTitleOffset(1.05);
    hAxis->GetXaxis()->SetLabelOffset(0.02);
    hAxis->GetXaxis()->SetTitleSize(0.075);
    hAxis->GetYaxis()->SetLabelSize(0.065);
    hAxis->GetYaxis()->SetTitle("C(k)");
    hAxis->GetYaxis()->CenterTitle();
    hAxis->GetYaxis()->SetTitleOffset(Yoffset);
    hAxis->GetYaxis()->SetTitleSize(0.075);

    if(WhichAxis==1){
        hAxis->GetYaxis()->SetRangeUser(0.50, 2.0);
        hAxis->GetYaxis()->SetLimits(0.50, 2.0);
    }
    else if(WhichAxis==2){
        hAxis->GetYaxis()->SetRangeUser(0.85, 1.15);
        hAxis->GetYaxis()->SetLimits(0.85, 1.15);
    }
    else if(WhichAxis==3){
        hAxis->GetYaxis()->SetRangeUser(0.70, 1.1);
        hAxis->GetYaxis()->SetLimits(0.70, 1.1);
    }
    else{
        hAxis->GetYaxis()->SetRangeUser(0.97, 4);
        hAxis->GetYaxis()->SetLimits(0.97, 4);
    }


    hAxis->GetXaxis()->SetNdivisions(504);

    TH1F* hAxisZoom = new TH1F("hAxisZoom", "hAxisZoom", NumPlotPts,Plot_kMin,Plot_kMax);
    hAxisZoom->SetStats(false);
    hAxisZoom->SetTitle("");
    hAxisZoom->GetXaxis()->SetLabelSize(0.065);
    hAxisZoom->GetXaxis()->SetTitle("k (MeV)");
    hAxisZoom->GetXaxis()->CenterTitle();
    hAxisZoom->GetXaxis()->SetTitleOffset(1.3);
    hAxisZoom->GetXaxis()->SetLabelOffset(0.02);
    hAxisZoom->GetXaxis()->SetTitleSize(0.075);
    hAxisZoom->GetYaxis()->SetLabelSize(0.065);
    hAxisZoom->GetYaxis()->SetTitle("C(k)");
    hAxisZoom->GetYaxis()->CenterTitle();
    hAxisZoom->GetYaxis()->SetTitleOffset(Yoffset/1.5);
    hAxisZoom->GetYaxis()->SetTitleSize(0.075);

    hAxisZoom->GetXaxis()->SetNdivisions(504);
    hAxisZoom->GetYaxis()->SetRangeUser(0.85, 1.15);
    hAxisZoom->GetYaxis()->SetLimits(0.85, 1.15);
    hAxisZoom->GetYaxis()->SetNdivisions(603);

    TPaveText** myPT = new TPaveText* [NumRad];
    for(unsigned uRad=0; uRad<NumRad; uRad++){
        myPT[uRad] = new TPaveText(0.5,0.885,0.95,0.975, "blNDC");//lbrt
        myPT[uRad]->SetName("PT12");
        myPT[uRad]->SetBorderSize(1);
        myPT[uRad]->SetTextSize(0.08/1.4);
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
        DrawBoard_Ck_pp.AddSubPadTL(CurrentPos,CurrentPos+PadLength,0.08,1);//lrtb
        //DrawBoard_Ck_pp.AddSubPadTL(CurrentPos,CurrentPos+PadLength,0.08,0.63);//lrtb
        //DrawBoard_Ck_pp.AddSubPadTL(CurrentPos,CurrentPos+PadLength,0.63,1);//lrtb
        CurrentPos+=PadLength;
//printf("uRad=%u; %.2f <--> %.2f\n", uRad, );
    }
    DrawBoard_Ck_pp.AddSubPadTL(0,1,0,0.075);//lrtb

    for(unsigned uRad=0; uRad<NumRad; uRad++){
        DrawBoard_Ck_pp.SetMargin(uRad, 0.07, 0.02, 0.15, 0.04);//lrbt
        //DrawBoard_Ck_pp.SetMargin(uRad*2+1, 0.07, 0.01, 0.14, 0);//lrbt
    }
    DrawBoard_Ck_pp.SetMargin(NumRad, 0, 0, 0, 0);//lrbt

    for(unsigned uRad=0; uRad<NumRad; uRad++){
        DrawBoard_Ck_pp.GetCanvas()->SetName(TString::Format("DrawBoard_Ck_XimProton_%u",uRad));
        DrawBoard_Ck_pp.cd(uRad);
        hAxis->Draw("AXIS");
        UnitLine->Draw("same");
        for(unsigned uSce=0; uSce<NumScenarios; uSce++) grCk[uSce][uRad].Draw("same,C");
        for(unsigned uSce=0; uSce<NumScenarios; uSce++) myPT[uRad]->Draw("same");

        DrawBoard_Ck_pp.SetLabelSize(uRad, hAxis->GetXaxis(), 17);
        DrawBoard_Ck_pp.SetLabelSize(uRad, hAxis->GetYaxis(), 17);
        DrawBoard_Ck_pp.SetTitleSize(uRad, hAxis->GetXaxis(), 17);
        DrawBoard_Ck_pp.SetTitleSize(uRad, hAxis->GetYaxis(), 17);
    }
    //for(unsigned uRad=0; uRad<NumRad; uRad++){
    //    DrawBoard_Ck_pp.GetCanvas()->SetName(TString::Format("DrawBoard_Ck_KmProton_%u",uRad));
    //    DrawBoard_Ck_pp.cd(uRad*2+1);
    //    hAxisZoom->Draw("AXIS");
    //    UnitLine->Draw("same");
    //    for(unsigned uSce=0; uSce<NumScenarios; uSce++) grCk[uSce][uRad].Draw("same,C");

    //    DrawBoard_Ck_pp.SetLabelSize(uRad*2+1, hAxisZoom->GetXaxis(), 17);
    //    DrawBoard_Ck_pp.SetLabelSize(uRad*2+1, hAxisZoom->GetYaxis(), 17);
    //    DrawBoard_Ck_pp.SetTitleSize(uRad*2+1, hAxisZoom->GetXaxis(), 17);
    //    DrawBoard_Ck_pp.SetTitleSize(uRad*2+1, hAxisZoom->GetYaxis(), 17);
    //}

    TLegend* myLegend = new TLegend(0.12,0,0.95,0.96);//lbrt
    myLegend->SetName("myLegend");
    myLegend->SetTextSize(0.65);
    myLegend->SetNColumns(NumScenarios);
    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        myLegend->AddEntry(&grCk[uSce][0],  ScenarioName[uSce]);
    }

    DrawBoard_Ck_pp.cd(NumRad);
    myLegend->Draw();

    gStyle->SetLineWidth(2.5);

    TString ChannelName = "";
    if(WhichChannel==0){
        ChannelName = "_I0S0";
    }
    else if(WhichChannel==1){
        ChannelName = "_I0S1";
    }
    else if(WhichChannel==2){
        ChannelName = "_I1S0";
    }
    else if(WhichChannel==3){
        ChannelName = "_I1S1";
    }
    else if(WhichChannel==10){
        ChannelName = "_ISavg";
    }

    DrawBoard_Ck_pp.GetCanvas()->SaveAs(TString::Format("%sPlotCk_XimProton%s.png",OutputDir.Data(),ChannelName.Data()));
    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        for(unsigned uRad=0; uRad<NumRad; uRad++){
            grCk[uSce][uRad].SetLineWidth(grCk[uSce][uRad].GetLineWidth()/2.5);
        }
    }
    hAxis->SetLineWidth(hAxis->GetLineWidth()/2.5);
    UnitLine->SetLineWidth(UnitLine->GetLineWidth()/2.5);
    gStyle->SetLineWidth(1);
    DrawBoard_Ck_pp.GetCanvas()->SaveAs(TString::Format("%sPlotCk_XimProton%s.pdf",OutputDir.Data(),ChannelName.Data()));

    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        for(unsigned uRad=0; uRad<NumRad; uRad++){
            grCk[uSce][uRad].SetLineWidth(grCk[uSce][uRad].GetLineWidth()*2.5);
        }
    }
    hAxis->SetLineWidth(hAxis->GetLineWidth()*2.5);
    UnitLine->SetLineWidth(UnitLine->GetLineWidth()*2.5);

    delete hAxis;
    delete hAxisZoom;
    delete UnitLine;
    delete myLegend;

    for(unsigned uRad=0; uRad<NumRad; uRad++) delete myPT[uRad];
    delete [] myPT;

    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        delete [] grCk[uSce];
    }
    delete [] grCk;

}

//! The lambdaPar is computed as follows:
//https://arxiv.org/pdf/1406.3206.pdf claims that around 1/3 of Xis come from the Xi1530
//https://arxiv.org/pdf/1204.0282.pdf and https://arxiv.org/pdf/0912.2855.pdf report that the Ω- is c.a. 10x less abundant compared to Xi
//if we assume (for no good reason) 93% purity of the Xi and 99% of the protons, than
//we should come roughly to the value of 46% prompt pXi pairs and 7% of contamination due to missidentifications.
void CATSpaper::PlotCk_XiMinusProton_Corrected1(){

    const double lambdaPar = 0.52;
    const double lambdaPurity = 0.11;
    const double lambdaFeed = 1.-lambdaPar-lambdaPurity;

    const double Plot_kMin = 0;
    const double Plot_kMax = 280;
    const unsigned NumPlotPts = 70;
    //CATShisto<double> CatHist(NumPlotPts,Plot_kMin,Plot_kMax);

    const unsigned NumRad = 3;
    const double Radius[NumRad] = {2.0,1.2,0.85};

//0= I0 S0
//1= I0 S1
//2= I1 S0
//3= I1 S1
//10=AVERAGED
//anything else is the full potential
const int WhichChannel = 100;

    //double CkVal;
    //double Momentum;
    DLM_DtColor DlmCol;

    const unsigned NumScenarios = 2;
    TString ScenarioName[NumScenarios];
    ScenarioName[0] = "Coulomb only";
    if(WhichChannel==0) ScenarioName[1] = "Coulomb + Strong(I=0,S=0)";
    else if(WhichChannel==1) ScenarioName[1] = "Coulomb + Strong(I=0,S=1)";
    else if(WhichChannel==2) ScenarioName[1] = "Coulomb + Strong(I=1,S=0)";
    else if(WhichChannel==3) ScenarioName[1] = "Coulomb + Strong(I=1,S=1)";
    else if(WhichChannel==10) ScenarioName[1] = "Coulomb + Strong(I,S avg)";
    else ScenarioName[1] = "Coulomb + Strong";
    //ScenarioName[2] = "Reid SC as in CRAB";

    TString ScenarioShortName[NumScenarios];
    ScenarioShortName[0] = "NoSI";
    if(WhichChannel==0) ScenarioShortName[1] = "I=0,S=0";
    else if(WhichChannel==1) ScenarioShortName[1] = "I=0,S=1";
    else if(WhichChannel==2) ScenarioShortName[1] = "I=1,S=0";
    else if(WhichChannel==3) ScenarioShortName[1] = "I=1,S=1";
    else if(WhichChannel==10) ScenarioShortName[1] = "I,S avg";
    else ScenarioShortName[1] = "Full";

    double SourcePars_Xip[4] = {0,0,0,Radius[0]};

    CATSparameters cSourcePars_Xip(CATSparameters::tSource,1,true);
    cSourcePars_Xip.SetParameters(&SourcePars_Xip[3]);

    CATS MyCat;
    MyCat.SetUseAnalyticSource(true);
    MyCat.SetAnaSource(GaussSource, cSourcePars_Xip);

    MyCat.SetExcludeFailedBins(false);
    MyCat.SetMomBins(NumPlotPts,Plot_kMin,Plot_kMax);

    MyCat.SetNumChannels(4);
    MyCat.SetNumPW(0,1);
    MyCat.SetNumPW(1,1);
    MyCat.SetNumPW(2,1);
    MyCat.SetNumPW(3,1);
    MyCat.SetSpin(0,0);//I=0; S=0
    MyCat.SetSpin(1,1);//I=0; S=1
    MyCat.SetSpin(2,0);//I=1; S=0
    MyCat.SetSpin(3,1);//I=1; S=1
    if(WhichChannel==0){
        MyCat.SetChannelWeight(0, 1.);
        MyCat.SetChannelWeight(1, 0.);
        MyCat.SetChannelWeight(2, 0.);
        MyCat.SetChannelWeight(3, 0.);
    }
    else if(WhichChannel==1){
        MyCat.SetChannelWeight(0, 0.);
        MyCat.SetChannelWeight(1, 1.);
        MyCat.SetChannelWeight(2, 0.);
        MyCat.SetChannelWeight(3, 0.);
    }
    else if(WhichChannel==2){
        MyCat.SetChannelWeight(0, 0.);
        MyCat.SetChannelWeight(1, 0.);
        MyCat.SetChannelWeight(2, 1.);
        MyCat.SetChannelWeight(3, 0.);
    }
    else if(WhichChannel==3){
        MyCat.SetChannelWeight(0, 0.);
        MyCat.SetChannelWeight(1, 0.);
        MyCat.SetChannelWeight(2, 0.);
        MyCat.SetChannelWeight(3, 1.);
    }
    else if(WhichChannel==10){
        MyCat.SetChannelWeight(0, 1.);
        MyCat.SetChannelWeight(1, 0.);
        MyCat.SetChannelWeight(2, 0.);
        MyCat.SetChannelWeight(3, 0.);
    }
    else{
        MyCat.SetChannelWeight(0, 1./8.);
        MyCat.SetChannelWeight(1, 3./8.);
        MyCat.SetChannelWeight(2, 1./8.);
        MyCat.SetChannelWeight(3, 3./8.);
    }

    MyCat.SetQ1Q2(-1);
    MyCat.SetPdgId(2212, 3312);//check the pdg id of the kaon
    //const double Mass_p = 938.272;
    const double Mass_Xim = 1321.7;
    MyCat.SetRedMass( (Mass_p*Mass_Xim)/(Mass_p+Mass_Xim) );

    MyCat.SetMaxRad(64);
    MyCat.SetMaxRho(32);

    const TString SigmaMatrixFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_7TeV/ResolutionMatrices/LHC14j4bcde_MeV.root";
    const double TrMatUnitConv = 1;//1000

    //TH2F* hSigma_pp;
    TH2F* hSigma_pL;

    TFile* FileSigma = new TFile(SigmaMatrixFileName, "read");
    FileSigma->cd();
    //hSigma_pp = (TH2F*)FileSigma->Get("hSigma_pp_InMeV");
    hSigma_pL = (TH2F*)FileSigma->Get("hSigma_pL_InMeV");

    //1/FractionOfBins th number of original bins of the correction matrix are taken into account
    //originally: 3000 MeV => 1/10 should be good most of the time
    const int FractionOfBins = 6;

    TH2F* hSigma_pL_InFuckingMeV = new TH2F("hSigma_pL_InFuckingMeV", "hSigma_pL_InFuckingMeV",
        hSigma_pL->GetNbinsX()/FractionOfBins, hSigma_pL->GetXaxis()->GetBinLowEdge(1)*TrMatUnitConv,hSigma_pL->GetXaxis()->GetBinUpEdge(hSigma_pL->GetNbinsX()/FractionOfBins)*TrMatUnitConv,
        hSigma_pL->GetNbinsY()/FractionOfBins, hSigma_pL->GetYaxis()->GetBinLowEdge(1)*TrMatUnitConv,hSigma_pL->GetXaxis()->GetBinUpEdge(hSigma_pL->GetNbinsY()/FractionOfBins)*TrMatUnitConv);

    for(int iBinX=1; iBinX<=hSigma_pL->GetNbinsX()/FractionOfBins; iBinX++){
        for(int iBinY=1; iBinY<=hSigma_pL->GetNbinsY()/FractionOfBins; iBinY++){
            hSigma_pL_InFuckingMeV->SetBinContent(iBinX, iBinY, hSigma_pL->GetBinContent(iBinX, iBinY));
        }
    }

    //in bins of Scenario/Radius
    TGraph** grCk;
    grCk = new TGraph* [NumScenarios];
//! DO THE INTEGRATED BIN HERE, JUST FOR CROSS CHECK
    //TGraph** grCk;
    //grCk = new TGraph* [NumScenarios];

    double pXimPotParsISavg[10]={0,0,pXim_LatticeAvg,12,0,-1,1,0,0,0};
    double pXimPotParsI01S0[10]={0,0,pXim_Lattice,12,0,-1,1,0,0,0};
    double pXimPotParsI03S1[10]={0,0,pXim_Lattice,12,0,-1,1,1,0,1};
    double pXimPotParsI11S0[10]={0,0,pXim_Lattice,6,1,1,1,0,0,0};
    double pXimPotParsI13S1[10]={0,0,pXim_Lattice,6,1,1,1,1,0,1};

    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        grCk[uSce] = new TGraph[NumRad];
        switch(uSce){
        case 0 :
            MyCat.SetShortRangePotential(0,0,ZeroPotential,pXimPotParsI01S0);
            MyCat.SetShortRangePotential(1,0,ZeroPotential,pXimPotParsI03S1);
            MyCat.SetShortRangePotential(2,0,ZeroPotential,pXimPotParsI11S0);
            MyCat.SetShortRangePotential(3,0,ZeroPotential,pXimPotParsI13S1);
            break;
        case 1 :
            if(WhichChannel==10){
                MyCat.SetShortRangePotential(0,0,fDlmPot,pXimPotParsISavg);
                MyCat.SetShortRangePotential(1,0,ZeroPotential,pXimPotParsI03S1);
                MyCat.SetShortRangePotential(2,0,ZeroPotential,pXimPotParsI11S0);
                MyCat.SetShortRangePotential(3,0,ZeroPotential,pXimPotParsI13S1);
            }
            else{
                MyCat.SetShortRangePotential(0,0,fDlmPot,pXimPotParsI01S0);
                MyCat.SetShortRangePotential(1,0,fDlmPot,pXimPotParsI03S1);
                MyCat.SetShortRangePotential(2,0,fDlmPot,pXimPotParsI11S0);
                MyCat.SetShortRangePotential(3,0,fDlmPot,pXimPotParsI13S1);
            }
            break;
        default :
            printf("Default Xi- fuck!\n");
            break;
        }
        MyCat.KillTheCat(CATS::kAllChanged);

        for(unsigned uRad=0; uRad<NumRad; uRad++){
            MyCat.SetAnaSource(0, Radius[uRad]);
            MyCat.KillTheCat();

            DLM_Ck* Ck_pXi = new DLM_Ck(1,0,MyCat);
            Ck_pXi->Update();
            DLM_CkDecomposition CkDec_pXi("pXim",2,*Ck_pXi,hSigma_pL_InFuckingMeV);
            CkDec_pXi.AddContribution(0,lambdaFeed,DLM_CkDecomposition::cFeedDown,NULL);
            CkDec_pXi.AddContribution(1,lambdaPurity,DLM_CkDecomposition::cFake);
            CkDec_pXi.Update();

            grCk[uSce][uRad].Set(NumPlotPts);
            grCk[uSce][uRad].SetName(TString::Format("grCk_%u_%u",uSce,uRad));
            grCk[uSce][uRad].SetMarkerColor(DlmCol.GetColor(1));
            grCk[uSce][uRad].SetMarkerStyle(20);
            grCk[uSce][uRad].SetMarkerSize(0);
            grCk[uSce][uRad].SetLineColor(DlmCol.GetColor(1+uSce));
            //the Coulomb only case
            //grCk[uSce][uRad].SetLineWidth(8.-4.*double(uSce)/double(NumScenarios-1));
            grCk[uSce][uRad].SetLineWidth(8.);
            grCk[uSce][uRad].SetFillColor(kWhite);

            for(unsigned uMomBin=0; uMomBin<NumPlotPts; uMomBin++){
                //grCk[uSce][uRad].SetPoint(uMomBin, MyCat.GetMomentum(uMomBin), MyCat.GetCorrFun(uMomBin));
                grCk[uSce][uRad].SetPoint(uMomBin, MyCat.GetMomentum(uMomBin), CkDec_pXi.EvalCk(MyCat.GetMomentum(uMomBin)));
            }
            if(uSce==0){
                //grCk[uSce][uRad].SetLineColor(kGray+1);
                grCk[uSce][uRad].SetLineStyle(2);
            }

            delete Ck_pXi;
        }
    }


    char* fBuffer = new char [512];
    FILE * pFile;
    pFile = fopen (TString::Format("%sPlotCk_XimProton_Corrected1.txt",OutputDir.Data()),"w");
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

    const double Yoffset = 1.25;
    TH1F* hAxis = new TH1F("hAxis", "hAxis", NumPlotPts,Plot_kMin,Plot_kMax);
    hAxis->SetStats(false);
    hAxis->SetTitle("");
    hAxis->GetXaxis()->SetLabelSize(0.065);
    hAxis->GetXaxis()->SetTitle("k (MeV)");
    hAxis->GetXaxis()->CenterTitle();
    hAxis->GetXaxis()->SetTitleOffset(1.05);
    hAxis->GetXaxis()->SetLabelOffset(0.02);
    hAxis->GetXaxis()->SetTitleSize(0.075);
    hAxis->GetYaxis()->SetLabelSize(0.065);
    hAxis->GetYaxis()->SetTitle("C(k)");
    hAxis->GetYaxis()->CenterTitle();
    hAxis->GetYaxis()->SetTitleOffset(Yoffset);
    hAxis->GetYaxis()->SetTitleSize(0.075);

    hAxis->GetYaxis()->SetRangeUser(0.95, 4);
    hAxis->GetYaxis()->SetLimits(0.95, 4);

    hAxis->GetXaxis()->SetNdivisions(504);

    TH1F* hAxisZoom = new TH1F("hAxisZoom", "hAxisZoom", NumPlotPts,Plot_kMin,Plot_kMax);
    hAxisZoom->SetStats(false);
    hAxisZoom->SetTitle("");
    hAxisZoom->GetXaxis()->SetLabelSize(0.065);
    hAxisZoom->GetXaxis()->SetTitle("k (MeV)");
    hAxisZoom->GetXaxis()->CenterTitle();
    hAxisZoom->GetXaxis()->SetTitleOffset(1.3);
    hAxisZoom->GetXaxis()->SetLabelOffset(0.02);
    hAxisZoom->GetXaxis()->SetTitleSize(0.075);
    hAxisZoom->GetYaxis()->SetLabelSize(0.065);
    hAxisZoom->GetYaxis()->SetTitle("C(k)");
    hAxisZoom->GetYaxis()->CenterTitle();
    hAxisZoom->GetYaxis()->SetTitleOffset(Yoffset/1.5);
    hAxisZoom->GetYaxis()->SetTitleSize(0.075);

    hAxisZoom->GetXaxis()->SetNdivisions(504);
    hAxisZoom->GetYaxis()->SetRangeUser(0.85, 1.15);
    hAxisZoom->GetYaxis()->SetLimits(0.85, 1.15);
    hAxisZoom->GetYaxis()->SetNdivisions(603);

    TPaveText** myPT = new TPaveText* [NumRad];
    for(unsigned uRad=0; uRad<NumRad; uRad++){
        myPT[uRad] = new TPaveText(0.5,0.885,0.95,0.975, "blNDC");//lbrt
        myPT[uRad]->SetName("PT12");
        myPT[uRad]->SetBorderSize(1);
        myPT[uRad]->SetTextSize(0.08/1.4);
        myPT[uRad]->SetFillColor(kWhite);
        myPT[uRad]->SetTextFont(22);
        //myPT[uRad]->AddText("Gauss source");
        myPT[uRad]->AddText(TString::Format("R = %.2f fm", Radius[uRad]));
    }

    TPaveText* ptCorrected0 = new TPaveText(0.4,0.72,0.95,0.885, "blNDC");//lbrt
    ptCorrected0->SetName("ptCorrected0");
    ptCorrected0->SetBorderSize(1);
    ptCorrected0->SetTextSize(0.075/1.4);
    ptCorrected0->SetFillColor(kWhite);
    ptCorrected0->SetTextFont(22);
    ptCorrected0->AddText("Corrections were");
    ptCorrected0->AddText("applied for:");


    TPaveText* ptCorrected = new TPaveText(0.4,0.45,0.95,0.72, "blNDC");//lbrt
    ptCorrected->SetName("ptCorrected");
    ptCorrected->SetBorderSize(1);
    ptCorrected->SetTextSize(0.065/1.4);
    ptCorrected->SetFillColor(kWhite);
    ptCorrected->SetTextFont(22);
    //ptCorrected->SetTextColor(kRed+1);
    ptCorrected->SetTextAlign(11);
    //ptCorrected->AddText("'Rule of thumb'");
    //ptCorrected->AddText("corrections for:");
    ptCorrected->AddText("* Momentum resolution");
    ptCorrected->AddText("* Residuals assuming:");
    ptCorrected->AddText(TString::Format("   %.0f%% primary p#Xi^{#minus} pairs",lambdaPar*100));
    ptCorrected->AddText(TString::Format("   %.0f%% flat residuals",(1.-lambdaPar)*100));

    TF1* UnitLine = new TF1("UnitLine", "1",Plot_kMin,Plot_kMax);
    UnitLine->SetLineColor(kGray+1);
    UnitLine->SetLineWidth(6);

    DLM_SubPads DrawBoard_Ck_pp(2560,1080);
    double PadLength = 1./double(NumRad);
    double CurrentPos = 0;
    for(unsigned uRad=0; uRad<NumRad; uRad++){
        DrawBoard_Ck_pp.AddSubPadTL(CurrentPos,CurrentPos+PadLength,0.08,1);//lrtb
        //DrawBoard_Ck_pp.AddSubPadTL(CurrentPos,CurrentPos+PadLength,0.08,0.63);//lrtb
        //DrawBoard_Ck_pp.AddSubPadTL(CurrentPos,CurrentPos+PadLength,0.63,1);//lrtb
        CurrentPos+=PadLength;
//printf("uRad=%u; %.2f <--> %.2f\n", uRad, );
    }
    DrawBoard_Ck_pp.AddSubPadTL(0,1,0,0.075);//lrtb

    for(unsigned uRad=0; uRad<NumRad; uRad++){
        DrawBoard_Ck_pp.SetMargin(uRad, 0.07, 0.02, 0.15, 0.04);//lrbt
        //DrawBoard_Ck_pp.SetMargin(uRad*2+1, 0.07, 0.01, 0.14, 0);//lrbt
    }
    DrawBoard_Ck_pp.SetMargin(NumRad, 0, 0, 0, 0);//lrbt

    for(unsigned uRad=0; uRad<NumRad; uRad++){
        DrawBoard_Ck_pp.GetCanvas()->SetName(TString::Format("DrawBoard_Ck_XimProton_Corrected1_%u",uRad));
        DrawBoard_Ck_pp.cd(uRad);
        hAxis->Draw("AXIS");
        UnitLine->Draw("same");
        for(unsigned uSce=0; uSce<NumScenarios; uSce++) grCk[uSce][uRad].Draw("same,C");
        for(unsigned uSce=0; uSce<NumScenarios; uSce++) myPT[uRad]->Draw("same");

        DrawBoard_Ck_pp.SetLabelSize(uRad, hAxis->GetXaxis(), 17);
        DrawBoard_Ck_pp.SetLabelSize(uRad, hAxis->GetYaxis(), 17);
        DrawBoard_Ck_pp.SetTitleSize(uRad, hAxis->GetXaxis(), 17);
        DrawBoard_Ck_pp.SetTitleSize(uRad, hAxis->GetYaxis(), 17);
    }
    //for(unsigned uRad=0; uRad<NumRad; uRad++){
    //    DrawBoard_Ck_pp.GetCanvas()->SetName(TString::Format("DrawBoard_Ck_KmProton_%u",uRad));
    //    DrawBoard_Ck_pp.cd(uRad*2+1);
    //    hAxisZoom->Draw("AXIS");
    //    UnitLine->Draw("same");
    //    for(unsigned uSce=0; uSce<NumScenarios; uSce++) grCk[uSce][uRad].Draw("same,C");

    //    DrawBoard_Ck_pp.SetLabelSize(uRad*2+1, hAxisZoom->GetXaxis(), 17);
    //    DrawBoard_Ck_pp.SetLabelSize(uRad*2+1, hAxisZoom->GetYaxis(), 17);
    //    DrawBoard_Ck_pp.SetTitleSize(uRad*2+1, hAxisZoom->GetXaxis(), 17);
    //    DrawBoard_Ck_pp.SetTitleSize(uRad*2+1, hAxisZoom->GetYaxis(), 17);
    //}

    TLegend* myLegend = new TLegend(0.12,0,0.95,0.96);//lbrt
    myLegend->SetName("myLegend");
    myLegend->SetTextSize(0.65);
    myLegend->SetNColumns(NumScenarios);
    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        myLegend->AddEntry(&grCk[uSce][0],  ScenarioName[uSce]);
    }

    DrawBoard_Ck_pp.cd(NumRad);
    myLegend->Draw();

    DrawBoard_Ck_pp.cd(NumRad/2);
    ptCorrected0->Draw("same");
    ptCorrected->Draw("same");

    gStyle->SetLineWidth(2.5);

    TString ChannelName = "";
    if(WhichChannel==0){
        ChannelName = "_I0S0";
    }
    else if(WhichChannel==1){
        ChannelName = "_I0S1";
    }
    else if(WhichChannel==2){
        ChannelName = "_I0S1";
    }
    else if(WhichChannel==3){
        ChannelName = "_I1S1";
    }
    else if(WhichChannel==10){
        ChannelName = "_ISavg";
    }

    DrawBoard_Ck_pp.GetCanvas()->SaveAs(TString::Format("%sPlotCk_XimProton%s_Corrected1.png",OutputDir.Data(),ChannelName.Data()));
    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        for(unsigned uRad=0; uRad<NumRad; uRad++){
            grCk[uSce][uRad].SetLineWidth(grCk[uSce][uRad].GetLineWidth()/2.5);
        }
    }
    hAxis->SetLineWidth(hAxis->GetLineWidth()/2.5);
    UnitLine->SetLineWidth(UnitLine->GetLineWidth()/2.5);
    gStyle->SetLineWidth(1);
    DrawBoard_Ck_pp.GetCanvas()->SaveAs(TString::Format("%sPlotCk_XimProton%s_Corrected1.pdf",OutputDir.Data(),ChannelName.Data()));

    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        for(unsigned uRad=0; uRad<NumRad; uRad++){
            grCk[uSce][uRad].SetLineWidth(grCk[uSce][uRad].GetLineWidth()*2.5);
        }
    }
    hAxis->SetLineWidth(hAxis->GetLineWidth()*2.5);
    UnitLine->SetLineWidth(UnitLine->GetLineWidth()*2.5);


    delete hAxis;
    delete hAxisZoom;
    delete ptCorrected;
    delete ptCorrected0;
    delete UnitLine;
    delete myLegend;

    for(unsigned uRad=0; uRad<NumRad; uRad++) delete myPT[uRad];
    delete [] myPT;

    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        delete [] grCk[uSce];
    }
    delete [] grCk;

    delete hSigma_pL_InFuckingMeV;
    delete FileSigma;


}



void CATSpaper::CompareCatsCrab(const int& DataSet){
    //0 = pp ALICE
    //1 = pL ALICE
    //2 = pp HADES
    //3 = pL HADES
    //4 = pp Gauss 3 fm 150 MeV

    double PlotKmin = 0;
    double PlotKmax = 160;


if(DataSet==1 || DataSet==3) {printf("Not tuned to plot that, sorry!\n"); return;}
    const unsigned NumMomBins =
        DataSet==0?60:
        DataSet==1?12:
        DataSet==2?60:
        DataSet==3?24:
                   30;
    const double kMin = 0;
    const double kMax =
        DataSet==0?240:
        DataSet==1?240:
        DataSet==2?240:
        DataSet==3?240:
                   120;

    const TString CrabFileName =
        DataSet==0?"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/LambdaProtonPotentials/CRAB/CRAB_DLM2/CompareWithCATS_2/qinv_pp3_EPOS.dat":
        DataSet==1?"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/LambdaProtonPotentials/CRAB/CRAB_DLM2/CompareWithCATS_2/qinv_pL3_EPOS.dat":
        DataSet==2?"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/LambdaProtonPotentials/CRAB/CRAB_DLM2/CompareWithCATS_2/qinv_pp3_UrQMD.dat":
        DataSet==3?"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/LambdaProtonPotentials/CRAB/CRAB_DLM2/CompareWithCATS_2/qinv_pL3_UrQMD.dat":
                   "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/LambdaProtonPotentials/CRAB/CRAB_DLM2/CompareWithCATS_2/qinv_pp3_TEGv2.dat";

    const unsigned CrabRebin=
        DataSet==0?24:
        DataSet==1?48:
        DataSet==2?24:
        DataSet==3?24:
                   24;

    const double Ymin =
    DataSet==0?0.0:
    DataSet==1?0.5:
    DataSet==2?0.0:
    DataSet==3?0.5:
               0.0;

    const double Ymax =
    DataSet==0?7.0:
    DataSet==1?4.0:
    DataSet==2?3.5:
    DataSet==3?4.0:
               1.6;

    const TString ExpSetInfo1 =
    DataSet==0?"Correlation: pp":
    DataSet==1?"Correlation: p#Lambda":
    DataSet==2?"Correlation: pp":
    DataSet==3?"Correlation: p#Lambda":
               "Correlation: pp"
    ;
    const TString ExpSetInfo2 =
    DataSet==0?"EPOS: pp @ #sqrt{s}=7 TeV":
    DataSet==1?"EPOS: pp @ #sqrt{s}=7 TeV":
    DataSet==2?"UrQMD: pNb @ #sqrt{s}=3.18 GeV":
    DataSet==3?"UrQMD: pNb @ #sqrt{s}=3.18 GeV":
               "PS GEN: 150 MeV, Gauss 3.0 fm"
    ;

    //char* NameBuffer = new char [128];

    //if(DataSet%2==0) strcpy(NameBuffer, "Reid");
    //else strcpy(NameBuffer, "Usmani");

    const TString LegendCats =
    DataSet==0?TString::Format("CATS + Argonne V18"):
    DataSet==1?TString::Format(" "):
    DataSet==2?TString::Format("CATS + Argonne V18"):
    DataSet==3?TString::Format(" "):
               TString::Format(" ")
    ;

    const TString LegendCrab =
    DataSet==0?TString::Format("CRAB + modified Reid SC"):
    DataSet==1?TString::Format(" "):
    DataSet==2?TString::Format("CRAB + modified Reid SC"):
    DataSet==3?TString::Format(" "):
               TString::Format(" ")
    ;

    const double GaussSourceSize =
    DataSet==0?1.2:
    DataSet==1?1.2:
    DataSet==2?2.02:
    DataSet==3?1.62:
               3.00
    ;

    const TString StrCrab = "Crab";
    const TString StrNumDenum = "CATS/CRAB";

    CATS MyCat;
    double Pars[4] = {0,0,0,GaussSourceSize};
    CATSparameters cSourcePars(CATSparameters::tSource,1,true);
    cSourcePars.SetParameters(&Pars[3]);
    MyCat.SetAnaSource(GaussSource, cSourcePars);
    MyCat.SetInputFileName(   DataSet==0?"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/Scratch9_OSCAR1997_AllLrzProtonFiles.f19":
                                    DataSet==1?"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/Scratch9_OSCAR1997_AllLrzLambdaFiles.f19":
                                    DataSet==2? "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/LambdaProtonPotentials/CppLocal/UrQMD_Crab_1/InputFiles/"
                                                "BigProtonSample/BigProtonSample_FM3.f19":
                                    DataSet==3? "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/LambdaProtonPotentials/CppLocal/UrQMD_Crab_1/InputFiles/"
                                                "BigLambdaSample/BigLambdaSample_pL_full.f19":
                                               "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/LambdaProtonPotentials/CRAB/CRAB_DLM2/ThermalSpectrum11_pp1/thermalTest2_PRS=3.00_PTS=3.00.f19");

    MyCat.SetExcludeFailedBins(false);
    MyCat.SetMomBins(NumMomBins,kMin,kMax);


    MyCat.SetNumChannels(4);
    MyCat.SetNumPW(0,2);
    MyCat.SetNumPW(1,2);
    MyCat.SetNumPW(2,2);
    MyCat.SetNumPW(3,2);
    MyCat.SetSpin(0,0);
    MyCat.SetSpin(1,1);
    MyCat.SetSpin(2,1);
    MyCat.SetSpin(3,1);
    MyCat.SetChannelWeight(0, 3./12.);
    MyCat.SetChannelWeight(1, 1./12.);
    MyCat.SetChannelWeight(2, 3./12.);
    MyCat.SetChannelWeight(3, 5./12.);

    MyCat.SetQ1Q2((DataSet%2)?0:1);
    MyCat.SetPdgId(2212, (DataSet%2)?3122:2212);
    double Mass1=938.272; double Mass2=(DataSet%2)?1115.683:938.272;
    MyCat.SetRedMass( (Mass1*Mass2)/(Mass1+Mass2) );
    //#,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
    double pLamPotPars1S0[10]={0,0,pL_UsmaniOli,0,0,0,0,0,0,0};
    double pLamPotPars3S1[10]={0,0,pL_UsmaniOli,0,0,0,0,1,0,1};
    if(DataSet%2){
        MyCat.SetShortRangePotential(0,0,fDlmPot,pLamPotPars1S0);
        MyCat.SetShortRangePotential(1,0,fDlmPot,pLamPotPars3S1);
    }
    else{
        MyCat.SetShortRangePotential(0,0,fDlmPot,ppPars1S0);
        MyCat.SetShortRangePotential(1,1,fDlmPot,ppPars3P0);
        MyCat.SetShortRangePotential(2,1,fDlmPot,ppPars3P1);
        MyCat.SetShortRangePotential(3,1,fDlmPot,ppPars3P2);
    }

    MyCat.SetMixingDepth(2);
    MyCat.SetUseAnalyticSource(false);
    MyCat.SetThetaDependentSource(true);
    MyCat.SetMaxPairsPerBin(16000);
    MyCat.SetMaxPairsToRead(256e6);

    MyCat.KillTheCat();

    DLM_DtColor DlmCol;

    TGraphErrors GrCats;
    GrCats.Set(NumMomBins);
    GrCats.SetName("GrCats");
    GrCats.SetMarkerColor(DlmCol.GetColor(1));
    GrCats.SetMarkerStyle(20);
    GrCats.SetMarkerSize(0);
    GrCats.SetLineColor(DlmCol.GetColor(1));
    GrCats.SetLineWidth(10);
    GrCats.SetFillColor(DlmCol.GetColor(1));
    GrCats.SetFillStyle(2001);//1001 solid fill

    double CkValue;
    double CkError;
    double Momentum;

    double Numerator;
    double NumError;
    double Denumerator;
    double DenumError;
    double NumDenum;
    //used for the asymmetric errors
    double NumMin;
    double NumMax;
    double DenumMin;
    double DenumMax;
    double RatioMin;
    double RatioMax;

    for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
        GrCats.SetPoint(uMomBin, MyCat.GetMomentum(uMomBin), MyCat.GetCorrFun(uMomBin));
        CkError = MyCat.GetCorrFunErr(uMomBin);
        //if(CkError<0.1)CkError=0.1;
        GrCats.SetPointError(uMomBin, 0, CkError);
    }

    TGraphErrors GrCrab;
    GrCrab.Set(NumMomBins);
    GrCrab.SetName("GrCrab");
    GrCrab.SetMarkerColor(DlmCol.GetColor(2));
    GrCrab.SetMarkerStyle(23);
    GrCrab.SetMarkerSize(0);
    GrCrab.SetLineColor(DlmCol.GetColor(2));
    GrCrab.SetLineWidth(8);

    GrCrab.SetFillColor(DlmCol.GetColor(2));
    GrCrab.SetFillStyle(2001);

    TGraphAsymmErrors GrRatioCrab;
    GrRatioCrab.Set(NumMomBins);
    GrRatioCrab.SetName("GrRatioCrab");
    GrRatioCrab.SetMarkerColor(kGray+1);
    GrRatioCrab.SetMarkerStyle(23);
    GrRatioCrab.SetMarkerSize(0);
    GrRatioCrab.SetLineColor(kGray+1);
    GrRatioCrab.SetLineWidth(8);
    GrRatioCrab.SetFillColor(kBlack);
    GrRatioCrab.SetFillStyle(2001);

    //if(ShowCrab){

        vector <string> fields;
        string line;
        //double RelMom;
        double CorrFun;
        double CorrFunErr;
        double Integral;
        double IntegralError;
        //double Chi2;
        bool FieldsCorr;
        unsigned uGet;
        unsigned CurrentBin;
        unsigned UsedBins;

        //hQinv = new TH1F("hQinv", "hQinv", NumBins, kMin, kMax);


        //ifstream CurrentQinvFile((TString("./TempCRAB/qinv_")+BaseFileName+TString(".dat")).Data());
        ifstream CurrentQinvFile(CrabFileName.Data());
        getline (CurrentQinvFile,line);
        getline (CurrentQinvFile,line);
        uGet = 0;
        CurrentBin=0;
        //Chi2=0;
        Integral=0;
        IntegralError=0;
        UsedBins=0;
        //ndof=0;
        while(getline (CurrentQinvFile,line)){
            uGet++;
            if(uGet%CrabRebin == 0){
                CkValue = Integral/double(UsedBins);
                CkError = sqrt(IntegralError)/double(UsedBins);
                //if(CkError<0.1)CkError=0.1;
                Momentum = MyCat.GetMomentum(CurrentBin);

                GrCrab.SetPoint(CurrentBin, Momentum, CkValue);
                //GrCrab.SetPointError(CurrentBin, 0.5*(MyCat.GetMomBinUpEdge(CurrentBin)-MyCat.GetMomBinLowEdge(CurrentBin)), CkError);
                GrCrab.SetPointError(CurrentBin, 0, CkError);

                Integral = 0;
                IntegralError = 0;
                UsedBins = 0;

                CurrentBin++;
            }
            if(CurrentBin>=NumMomBins) break;
            split( fields, line, is_any_of( " \n" ), token_compress_on );
            FieldsCorr = fields[0]=="";
            //RelMom = atof(fields[0+FieldsCorr].c_str());
            CorrFun = atof(fields[1+FieldsCorr].c_str());
            CorrFunErr = atof(fields[2+FieldsCorr].c_str());
            if(CorrFunErr!=-1){
                Integral += CorrFun;
                IntegralError += CorrFunErr*CorrFunErr;
                UsedBins++;
            }
        }
        CurrentQinvFile.close();

        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            GrCats.GetPoint(uMomBin, Momentum, Numerator);
            NumError = GrCats.GetErrorY(uMomBin);
            GrCrab.GetPoint(uMomBin, Momentum, Denumerator);
            DenumError = GrCrab.GetErrorY(uMomBin);
            NumDenum = Denumerator?(Numerator/Denumerator):0;
            NumMin = Numerator-NumError; if(NumMin<=0) NumMin=0;
            NumMax = Numerator+NumError;
            DenumMin = Denumerator-DenumError; if(DenumMin<=0) DenumMin=Denumerator*1e-9;
            DenumMax = Denumerator+DenumError;
            RatioMin = NumMin/DenumMax;
            RatioMax = NumMax/DenumMin;

            GrRatioCrab.SetPoint(uMomBin, Momentum, NumDenum);
            if( (NumDenum<0.2 || NumDenum>5) && fabs(DenumError/Denumerator)<1 ) GrRatioCrab.SetPointError(uMomBin, 0, 0, NumDenum-0, 2-NumDenum);
            else GrRatioCrab.SetPointError(uMomBin, 0, 0, NumDenum-RatioMin, RatioMax-NumDenum);
//printf("Momentum=%f; eup=%f; elow=%f\n",Momentum,NumDenum-RatioMin, RatioMax-NumDenum);
//printf(" Numerator=%.3f+/-%.3f; Denumerator=%.3f+/-%.3f; NumDenum=%f;\n",Numerator,NumError,Denumerator,DenumError,NumDenum);
        }

    //}

    char* fBuffer = new char [512];

    FILE * pFile;
    pFile = fopen (TString::Format("%sCompareCatsCrab_%s.txt",OutputDir.Data(),DataSet==0?"pp_ALICE":"pp_HADES"),"w");
    fprintf (pFile, "%18s","k (MeV)");

    sprintf(fBuffer,"Cats");
    fprintf (pFile,"%18s",fBuffer);

    sprintf(fBuffer,"CatsErr");
    fprintf (pFile,"%18s",fBuffer);

    sprintf(fBuffer,"Crab");
    fprintf (pFile,"%18s",fBuffer);

    sprintf(fBuffer,"CrabErr");
    fprintf (pFile,"%18s",fBuffer);

    sprintf(fBuffer,"Ratio");
    fprintf (pFile,"%18s",fBuffer);

    sprintf(fBuffer,"RatioErr");
    fprintf (pFile,"%18s",fBuffer);

    double xVal,yVal;
    for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
        fprintf (pFile, "\n%18.3f",MyCat.GetMomentum(uMomBin));

        GrCats.GetPoint(uMomBin,xVal,yVal);
        if(xVal!=MyCat.GetMomentum(uMomBin)) {printf("WTF MAN!\n");}
        fprintf (pFile, "%18.3f",yVal);

        yVal = GrCats.GetErrorY(uMomBin);
        fprintf (pFile, "%18.3f",yVal);

        GrCrab.GetPoint(uMomBin,xVal,yVal);
        if(xVal!=MyCat.GetMomentum(uMomBin)) {printf("WTF MAN!\n");}
        fprintf (pFile, "%18.3f",yVal);

        yVal = GrCrab.GetErrorY(uMomBin);
        fprintf (pFile, "%18.3f",yVal);

        GrRatioCrab.GetPoint(uMomBin,xVal,yVal);
        if(xVal!=MyCat.GetMomentum(uMomBin)) {printf("WTF MAN!\n");}
        fprintf (pFile, "%18.3f",yVal);

        yVal = GrRatioCrab.GetErrorY(uMomBin);
        fprintf (pFile, "%18.3f",yVal);
    }


    TH1F* hAxis = new TH1F("hAxis", "hAxis", 1, PlotKmin, PlotKmax);
    hAxis->SetStats(false);
    hAxis->SetTitle("");
    hAxis->GetXaxis()->SetLabelSize(0.065*1.25);
    hAxis->GetXaxis()->SetTitle("k (MeV)");
    hAxis->GetXaxis()->CenterTitle();
    hAxis->GetXaxis()->SetTitleOffset(1.12);
    hAxis->GetXaxis()->SetLabelOffset(0.02);
    hAxis->GetXaxis()->SetTitleSize(0.075*1.25);
    hAxis->GetYaxis()->SetLabelSize(0.06*1.25);
    hAxis->GetYaxis()->SetTitle("C(k)");
    hAxis->GetYaxis()->CenterTitle();
    hAxis->GetYaxis()->SetTitleOffset(0.5);
    hAxis->GetYaxis()->SetTitleSize(0.07*1.25);

    hAxis->GetYaxis()->SetRangeUser(Ymin, Ymax);
    hAxis->GetYaxis()->SetLimits(Ymin, Ymax);


    TH1F* hRatioAxis = new TH1F("hRatioAxis", "hRatioAxis", 1, PlotKmin, PlotKmax);
    hRatioAxis->SetStats(false);
    hRatioAxis->SetTitle("");
    hRatioAxis->GetXaxis()->SetLabelSize(0.06*1.25);
    hRatioAxis->GetXaxis()->SetTitle("k (MeV)");
    hRatioAxis->GetXaxis()->CenterTitle();
    hRatioAxis->GetXaxis()->SetTitleOffset(1.12);
    hRatioAxis->GetXaxis()->SetLabelOffset(0.02);
    hRatioAxis->GetXaxis()->SetTitleSize(0.075*1.25);
    hRatioAxis->GetYaxis()->SetLabelSize(0.06*1.25);
    hRatioAxis->GetYaxis()->SetTitle("Ratio");
    hRatioAxis->GetYaxis()->CenterTitle();
    hRatioAxis->GetYaxis()->SetTitleOffset(0.5/1.5);
    hRatioAxis->GetYaxis()->SetTitleSize(0.07*1.25);
    hRatioAxis->GetYaxis()->SetNdivisions(504);

    hRatioAxis->GetYaxis()->SetRangeUser(0.7, 1.26);
    hRatioAxis->GetYaxis()->SetLimits(0.7, 1.26);

    TF1* UnitLine = new TF1("UnitLine", "1", kMin, kMax);
    UnitLine->SetLineColor(kGray+1);
    UnitLine->SetLineWidth(5);

    const double DownShift = DataSet==(4)?0.2:0;
    const double LeftShift = DataSet==(4)?0.1:0;

    TPaveText* myPT = new TPaveText(0.58,0.85,1,1.00, "blNDC");//lbrt
    myPT->SetName("myPT");
    myPT->SetBorderSize(1);
    myPT->SetTextSize(0.045);
    myPT->SetFillColor(kWhite);
    myPT->SetTextFont(22);
    myPT->AddText(ExpSetInfo1);
    myPT->AddText(ExpSetInfo2);

    TPaveText* myPT3 = new TPaveText(0.78,0.32,0.955,0.39, "blNDC");//lbrt
    myPT3->SetName("myPT3");
    myPT3->SetBorderSize(1);
    myPT3->SetTextSize(0.045);
    myPT3->SetFillColor(kWhite);
    myPT3->SetTextFont(42);
    myPT3->AddText(StrNumDenum);

    TLegend* myLegend = new TLegend(0.58,0.62-DownShift,1.0-LeftShift,0.85-DownShift);//lbrt
    myLegend->SetName("myLegend");
    myLegend->SetTextSize(0.045);
    myLegend->AddEntry(&GrCats, LegendCats);
    myLegend->AddEntry(&GrCrab, LegendCrab);
    myLegend->AddEntry(&GrRatioCrab, "Ratio with errors");
    //myLegend->AddEntry(UnitLine, "Baseline");

    DLM_SubPads DrawBoard1(1920,1080);
    DrawBoard1.AddSubPadTL(0,1,0,0.6);//lrtb
    DrawBoard1.AddSubPadTL(0,1,0.6,1);//lrtb

    DrawBoard1.SetMargin(0, 0.12, 0.04, 0, 0.04);//lrbt
    DrawBoard1.SetMargin(1, 0.12, 0.04, 0.14, 0);//lrbt

    DrawBoard1.SetLogy(0,false);
    DrawBoard1.cd(0);

    hAxis->Draw("AXIS");
    UnitLine->Draw("same");
    GrCats.Draw("0,4,C,same");
    GrCrab.Draw("0,4,C,same");

    DrawBoard1.SetLabelSize(0, hAxis->GetXaxis(), 18*1.2);
    DrawBoard1.SetLabelSize(0, hAxis->GetYaxis(), 18*1.2);
    DrawBoard1.SetTitleSize(0, hAxis->GetXaxis(), 18*1.2);
    DrawBoard1.SetTitleSize(0, hAxis->GetYaxis(), 18*1.2);

    DrawBoard1.cd(1);
    hRatioAxis->Draw("AXIS");
    UnitLine->Draw("same");
    GrRatioCrab.Draw("0,3,L,same");

    DrawBoard1.SetLabelSize(1, hRatioAxis->GetXaxis(), 18*1.2);
    DrawBoard1.SetLabelSize(1, hRatioAxis->GetYaxis(), 18*1.2);
    DrawBoard1.SetTitleSize(1, hRatioAxis->GetXaxis(), 18*1.2);
    DrawBoard1.SetTitleSize(1, hRatioAxis->GetYaxis(), 18*1.2);

    DrawBoard1.GetCanvas()->SetName(TString::Format("CompareCats%s_DS%i", StrCrab.Data(), DataSet));
    DrawBoard1.GetCanvas()->cd(0);

    myPT->Draw("same");
    myPT3->Draw("same");
    myLegend->Draw("same");

    DrawBoard1.GetCanvas()->SaveAs(TString::Format("%sCompareCats%s_%s.png",
                                    OutputDir.Data(), StrCrab.Data(), DataSet==0?"pp_ALICE":"pp_HADES"));
    UnitLine->SetLineWidth(UnitLine->GetLineWidth()/2.);
    GrCats.SetLineWidth(GrCats.GetLineWidth()/2.);
    GrCrab.SetLineWidth(GrCrab.GetLineWidth()/2.);
    GrRatioCrab.SetLineWidth(GrRatioCrab.GetLineWidth()/2.);
    DrawBoard1.GetCanvas()->SaveAs(TString::Format("%sCompareCats%s_%s.pdf",
                                    OutputDir.Data(), StrCrab.Data(), DataSet==0?"pp_ALICE":"pp_HADES"));
    UnitLine->SetLineWidth(UnitLine->GetLineWidth()*2.);
    GrCats.SetLineWidth(GrCats.GetLineWidth()*2.);
    GrCrab.SetLineWidth(GrCrab.GetLineWidth()*2.);
    GrRatioCrab.SetLineWidth(GrRatioCrab.GetLineWidth()*2.);

    delete UnitLine;
    delete myPT;
    delete myPT3;
    delete myLegend;
    delete hAxis;
    delete hRatioAxis;

    delete [] fBuffer;

}


void CATSpaper::ExampleFit(){

    const double kFrom = 0;
    const double kTo = 160;
    const unsigned kNumStep = 40;

    const double DefaultRadius = 1.2;

    const int RanSeed = 22;
    TRandom3 rangen(RanSeed);
    const double MinError = 0.05;
    const double MaxError = 0.4;
    const double MinErrorAt = 120;
    const double SlopePar = (MinError-MaxError)/(MinErrorAt);

    double kVal;
    double CkVal;
    double ranval;
    double Sigma;

    CATSparameters cSourcePars(CATSparameters::tSource,1,true);
    cSourcePars.SetParameters(&ppSourcePars[3]);

    CATS MyCat;
    MyCat.SetUseAnalyticSource(true);
    MyCat.SetAnaSource(GaussSource, cSourcePars);

    MyCat.SetExcludeFailedBins(false);
    MyCat.SetMomBins(kNumStep,kFrom,kTo);

    MyCat.SetNumChannels(4);
    MyCat.SetNumPW(0,2);
    MyCat.SetNumPW(1,2);
    MyCat.SetNumPW(2,2);
    MyCat.SetNumPW(3,2);
    MyCat.SetSpin(0,0);
    MyCat.SetSpin(1,1);
    MyCat.SetSpin(2,1);
    MyCat.SetSpin(3,1);
    MyCat.SetChannelWeight(0, 3./12.);
    MyCat.SetChannelWeight(1, 1./12.);
    MyCat.SetChannelWeight(2, 3./12.);
    MyCat.SetChannelWeight(3, 5./12.);
    MyCat.SetShortRangePotential(0,0,fDlmPot,ppPars1S0);
    MyCat.SetShortRangePotential(1,1,fDlmPot,ppPars3P0);
    MyCat.SetShortRangePotential(2,1,fDlmPot,ppPars3P1);
    MyCat.SetShortRangePotential(3,1,fDlmPot,ppPars3P2);
    MyCat.SetQ1Q2(1);
    MyCat.SetPdgId(2212, 2212);
    MyCat.SetRedMass( 0.5*Mass_p );

    ppSourcePars[3] = DefaultRadius;
    MyCat.KillTheCat();

    double* MomBins = new double [kNumStep+1];
    for(unsigned uBin=0; uBin<=kNumStep; uBin++){
        MomBins[uBin] = MyCat.GetMomBinLowEdge(uBin);
    }

    //TH1F* hAxis = new TH1F("hAxis", "hAxis", kNumStep, MomBins);
    TH1F* hAV18GaussDefault = new TH1F("hAV18GaussDefault", "hAV18GaussDefault", kNumStep, MomBins);
    hAV18GaussDefault->SetStats(false);
    hAV18GaussDefault->SetTitle("");
    hAV18GaussDefault->GetXaxis()->SetLabelSize(0.065);
    hAV18GaussDefault->GetXaxis()->SetTitle("k (MeV)");
    hAV18GaussDefault->GetXaxis()->CenterTitle();
    hAV18GaussDefault->GetXaxis()->SetTitleOffset(1.1);
    hAV18GaussDefault->GetXaxis()->SetLabelOffset(0.02);
    hAV18GaussDefault->GetXaxis()->SetTitleSize(0.075);
    hAV18GaussDefault->GetYaxis()->SetLabelSize(0.065);
    hAV18GaussDefault->GetYaxis()->SetTitle("C(k)");
    hAV18GaussDefault->GetYaxis()->CenterTitle();
    hAV18GaussDefault->GetYaxis()->SetTitleOffset(0.7);
    hAV18GaussDefault->GetYaxis()->SetTitleSize(0.075);

    hAV18GaussDefault->GetYaxis()->SetRangeUser(0, 5);
    hAV18GaussDefault->GetYaxis()->SetLimits(0, 5);

    hAV18GaussDefault->GetXaxis()->SetNdivisions(504);

    hAV18GaussDefault->SetLineWidth(10);
    hAV18GaussDefault->SetLineColor(kBlack);

    for(unsigned uBin=0; uBin<kNumStep; uBin++){
        kVal = hAV18GaussDefault->GetBinCenter(uBin+1);
        CkVal = MyCat.EvalCorrFun(kVal);
        if(kVal>MinErrorAt) Sigma = MinError;
        else Sigma = SlopePar*kVal + MaxError;
        do ranval = rangen.Gaus(CkVal, Sigma);
        while(ranval<0);
        hAV18GaussDefault->SetBinContent(uBin+1, ranval);
        hAV18GaussDefault->SetBinError(uBin+1, Sigma);
    }
    //   TF1 * f = new TF1("f",fptr,&MyFunction::Evaluate,0,1,npar,"MyFunction","Evaluate");   // create TF1 class.
    FitKitty = &MyCat;
    MyCat.SetNotifications(CATS::nWarning);
    TF1* MyFit = new TF1("MyFit",this,&CATSpaper::Fitter1,kFrom,kTo,4,"CATSpaper","Fitter1");
    MyFit->SetParameter(0,1);//a
    MyFit->SetParLimits(0,0.5,2);
    MyFit->FixParameter(1,0);//b
    MyFit->FixParameter(2,1);//λ
    MyFit->SetParameter(3,1);//r
    MyFit->SetParLimits(3,0.5,2);
    TGraph grFit;
    grFit.Set(kNumStep);
    grFit.SetName("grFit");
    grFit.SetLineWidth(8);
    grFit.SetLineColor(kRed+1);
    grFit.SetMarkerSize(0);
    grFit.SetMarkerColor(kRed+1);
    grFit.SetFillColor(kWhite);

    hAV18GaussDefault->Fit(MyFit, "S, N, R, M");
    for(unsigned uBin=0; uBin<kNumStep; uBin++){
        grFit.SetPoint(uBin,MyCat.GetMomentum(uBin),MyFit->Eval(MyCat.GetMomentum(uBin)));
    }

    TF1* UnitLine = new TF1("UnitLine", "1", kFrom,kTo);
    UnitLine->SetLineWidth(8);
    UnitLine->SetLineColor(kGray+1);

    TLegend* myLegend = new TLegend(0.5,0.75,0.975,0.95);//lbrt
    myLegend->SetName("myLegend");
    myLegend->SetTextSize(0.05);
    myLegend->AddEntry(hAV18GaussDefault, TString::Format("Dummy data: R_{true}=%.2f fm", DefaultRadius));
    myLegend->AddEntry(&grFit, TString::Format(           "CATS fit: R_{fit} =%.2f#pm%.2f fm", MyFit->GetParameter(0),MyFit->GetParError(0)));

    TCanvas* cCk_pp = new TCanvas("cCk_pp", "cCk_pp", 1);
    cCk_pp->cd(0); cCk_pp->SetCanvasSize(1920, 1080); cCk_pp->SetMargin(0.12,0.025,0.18,0.05);//lrbt

    hAV18GaussDefault->Draw("axis");
    UnitLine->Draw("same");
    hAV18GaussDefault->Draw("same");
    grFit.Draw("same,C");
    myLegend->Draw("same");

    cCk_pp->SaveAs(TString::Format("%sExampleFit.png",OutputDir.Data()));
    hAV18GaussDefault->SetLineWidth(hAV18GaussDefault->GetLineWidth()/2.5);
    UnitLine->SetLineWidth(UnitLine->GetLineWidth()/2.5);
    grFit.SetLineWidth(grFit.GetLineWidth()/2.5);
    cCk_pp->SaveAs(TString::Format("%sExampleFit.pdf",OutputDir.Data()));
    hAV18GaussDefault->SetLineWidth(hAV18GaussDefault->GetLineWidth()*2.5);
    UnitLine->SetLineWidth(UnitLine->GetLineWidth()*2.5);
    grFit.SetLineWidth(grFit.GetLineWidth()*2.5);

    delete hAV18GaussDefault;
    delete MyFit;
    delete UnitLine;
    delete myLegend;
    delete cCk_pp;

}



void CATSpaper::FitHadesData_pp(){

    const bool ShiftHalfBin = true;

    const unsigned NumMomBins = 30;
    const double BinWidth = 4;
    const double kMin = ShiftHalfBin*BinWidth*0.5;
    const double kMax = kMin+double(NumMomBins)*BinWidth;

    CATS MyCat;
    FitKitty = &MyCat;

    double Pars[4] = {0,0,0,2};
    CATSparameters cSourcePars(CATSparameters::tSource,1,true);
    cSourcePars.SetParameters(&Pars[3]);
    MyCat.SetUseAnalyticSource(true);
    MyCat.SetAnaSource(GaussSource, cSourcePars);

    MyCat.SetExcludeFailedBins(false);
    MyCat.SetMomBins(NumMomBins,kMin,kMax);

    MyCat.SetNumChannels(4);
    MyCat.SetNumPW(0,2);
    MyCat.SetNumPW(1,2);
    MyCat.SetNumPW(2,2);
    MyCat.SetNumPW(3,2);
    MyCat.SetSpin(0,0);
    MyCat.SetSpin(1,1);
    MyCat.SetSpin(2,1);
    MyCat.SetSpin(3,1);
    MyCat.SetChannelWeight(0, 3./12.);
    MyCat.SetChannelWeight(1, 1./12.);
    MyCat.SetChannelWeight(2, 3./12.);
    MyCat.SetChannelWeight(3, 5./12.);

    MyCat.SetQ1Q2(1);
    MyCat.SetPdgId(2212,2212);
    MyCat.SetRedMass( 0.5*Mass_p );

    //ppPars1S0[2]=pp_ReidSC;
    //ppPars1S0[3]=0;
    MyCat.SetShortRangePotential(0,0,fDlmPot,ppPars1S0);
    MyCat.SetShortRangePotential(1,1,fDlmPot,ppPars3P0);
    MyCat.SetShortRangePotential(2,1,fDlmPot,ppPars3P1);
    MyCat.SetShortRangePotential(3,1,fDlmPot,ppPars3P2);

    MyCat.KillTheCat();
    MyCat.SetNotifications(CATS::nSilent);

    const TString OliFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/OliverPlotter/pp_qinv_fullycorrected.root";
    const TString OliHistoName = "C2_pp_fullyCorrected";
    const double UnitsConv = 1;
    TFile* OliFile = OliFileName!=""?new TFile(OliFileName, "read"):NULL;
    TH1F* OliHisto = OliFile?(TH1F*)OliFile->Get(OliHistoName):NULL;

    TH1F* hOliBinShift = new TH1F("hOliBinShift","hOliBinShift",NumMomBins+1,
                                  (OliHisto->GetBinLowEdge(1)-double(ShiftHalfBin)*OliHisto->GetBinWidth(1)*0.5)*UnitsConv,
                                  //(OliHisto->GetXaxis()->GetBinUpEdge(NumMomBins)-double(ShiftHalfBin)*OliHisto->GetBinWidth(NumMomBins)*0.5)*UnitsConv);
                                  (OliHisto->GetXaxis()->GetBinUpEdge(NumMomBins)+double(ShiftHalfBin)*OliHisto->GetBinWidth(NumMomBins)*0.5)*UnitsConv);

    if(hOliBinShift->GetBinWidth(1)!=OliHisto->GetBinWidth(1)*UnitsConv) printf("TROUBLE!\n");
    for(unsigned uBin=1; uBin<=NumMomBins+1; uBin++){
        hOliBinShift->SetBinContent(uBin, OliHisto->GetBinContent(uBin));
        hOliBinShift->SetBinError(uBin, OliHisto->GetBinError(uBin));
    }

    TH1F* hAxis = new TH1F("hAxis","hAxis",NumMomBins,kMin,kMax);
    hAxis->SetStats(false);
    hAxis->SetTitle("");
    hAxis->GetXaxis()->SetLabelSize(0.065);
    hAxis->GetXaxis()->SetTitle("k (MeV)");
    hAxis->GetXaxis()->CenterTitle();
    hAxis->GetXaxis()->SetTitleOffset(1.1);
    hAxis->GetXaxis()->SetLabelOffset(0.02);
    hAxis->GetXaxis()->SetTitleSize(0.075);
    hAxis->GetYaxis()->SetLabelSize(0.065);
    hAxis->GetYaxis()->SetTitle("C(k)");
    hAxis->GetYaxis()->CenterTitle();
    hAxis->GetYaxis()->SetTitleOffset(0.7);
    hAxis->GetYaxis()->SetTitleSize(0.075);

    hAxis->GetYaxis()->SetRangeUser(0, 2.8);
    hAxis->GetYaxis()->SetLimits(0, 2.8);

    //hOliBinShift->GetXaxis()->SetNdivisions(504);

    hOliBinShift->SetLineWidth(10);
    hOliBinShift->SetLineColor(kBlack);

    FitKitty = &MyCat;
    MyCat.SetNotifications(CATS::nWarning);
    TF1* MyFit = new TF1("MyFit",this,&CATSpaper::Fitter1,hOliBinShift->GetBinCenter(3), hOliBinShift->GetBinCenter(NumMomBins),4,"CATSpaper","Fitter1");
    MyFit->SetParameter(0,1);//a
    MyFit->SetParLimits(0,0.5,1.5);
//MyFit->FixParameter(0,1);//a
    MyFit->FixParameter(1,0);//b
    MyFit->FixParameter(2,1);//λ
    MyFit->SetParameter(3,2);//r
    MyFit->SetParLimits(3,1.5,2.5);
    TGraph grFit;
    grFit.Set(NumMomBins);
    grFit.SetName("grFit");
    grFit.SetLineWidth(8);
    grFit.SetLineColor(kRed+1);
    grFit.SetMarkerSize(0);
    grFit.SetMarkerColor(kRed+1);
    grFit.SetFillColor(kWhite);

    hOliBinShift->Fit(MyFit, "S, N, R, M");
    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        grFit.SetPoint(uBin,MyCat.GetMomentum(uBin),MyFit->Eval(MyCat.GetMomentum(uBin)));
    }

    char* fBuffer = new char [512];

    FILE * pFile;
    pFile = fopen (TString::Format("%sFitHadesData_pp.txt",OutputDir.Data()),"w");
    fprintf (pFile, "%18s","k (MeV)");

    sprintf(fBuffer,"Data");
    fprintf (pFile,"%18s",fBuffer);

    sprintf(fBuffer,"DataErr");
    fprintf (pFile,"%18s",fBuffer);

    sprintf(fBuffer,"CatsFit");
    fprintf (pFile,"%18s",fBuffer);

    double xVal,yVal;
    for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
        fprintf (pFile, "\n%18.3f",MyCat.GetMomentum(uMomBin));
//printf("k=%.3f; h=%.3f\n",MyCat.GetMomentum(uMomBin),hOliBinShift->GetBinCenter(uMomBin+2));

        //the +2 comes from the fact, that the first histo bin sits at k=0, while the first CATS bin is at k=4
        yVal = hOliBinShift->GetBinContent(uMomBin+2);
        fprintf (pFile, "%18.3f",yVal);

        yVal = hOliBinShift->GetBinError(uMomBin+2);
        fprintf (pFile, "%18.3f",yVal);

        grFit.GetPoint(uMomBin,xVal,yVal);
        if(xVal!=MyCat.GetMomentum(uMomBin)) {printf("WTF MAN!\n");}
        fprintf (pFile, "%18.3f",yVal);

    }

    TF1* UnitLine = new TF1("UnitLine", "1", kMin,kMax);
    UnitLine->SetLineWidth(8);
    UnitLine->SetLineColor(kGray+1);

    TLegend* myLegend = new TLegend(0.65,0.75,0.975,0.95);//lbrt
    myLegend->SetName("myLegend");
    myLegend->SetTextSize(0.05);
    myLegend->AddEntry(hOliBinShift, TString::Format("HADES data"));
    myLegend->AddEntry(&grFit, TString::Format("CATS fit with AV18"));

    TPaveText* myPT = new TPaveText(0.6,0.55,0.975,0.75, "blNDC");//lbrt
    myPT->SetName(TString::Format("myPT"));
    myPT->SetBorderSize(1);
    myPT->SetTextSize(0.05);
    myPT->SetFillColor(kWhite);
    myPT->SetTextFont(22);
    myPT->AddText(TString::Format("R_{HADES} = %.2f #pm %.2f fm", 2.02, 0.01));
    myPT->AddText(TString::Format("#color[%u]{R_{CATS} = %.2f #pm %.2f fm}", unsigned(kRed+1), MyFit->GetParameter(3),MyFit->GetParError(3)));
    //for(unsigned u=0; u<10; u++){
    //    const unsigned FROM=kRed;
    //    myPT->AddText(TString::Format("%u: #color[%u]{R_{fit} = %.2f #pm %.2f fm}", FROM+u, FROM+u, MyFit->GetParameter(3),MyFit->GetParError(3)));
    //}

    TCanvas* cCk_pp = new TCanvas("cCk_pp", "cCk_pp", 1);
    cCk_pp->cd(0); cCk_pp->SetCanvasSize(1920, 1080); cCk_pp->SetMargin(0.12,0.025,0.18,0.05);//lrbt

    hAxis->Draw("axis");
    UnitLine->Draw("same");
    hOliBinShift->Draw("same");
    grFit.Draw("same,C");
    myLegend->Draw("same");
    myPT->Draw("same");

    cCk_pp->SaveAs(TString::Format("%sFitHadesData_pp.png",OutputDir.Data()));
    hOliBinShift->SetLineWidth(hOliBinShift->GetLineWidth()/2.5);
    UnitLine->SetLineWidth(UnitLine->GetLineWidth()/2.5);
    grFit.SetLineWidth(grFit.GetLineWidth()/2.5);
    cCk_pp->SaveAs(TString::Format("%sFitHadesData_pp.pdf",OutputDir.Data()));
    hOliBinShift->SetLineWidth(hOliBinShift->GetLineWidth()*2.5);
    UnitLine->SetLineWidth(UnitLine->GetLineWidth()*2.5);
    grFit.SetLineWidth(grFit.GetLineWidth()*2.5);

    delete MyFit;
    delete hAxis;
    delete hOliBinShift;
    delete UnitLine;
    delete OliFile;
    delete myLegend;
    delete myPT;
    delete cCk_pp;

    delete [] fBuffer;

}


void CATSpaper::FitHadesData_pL(){

    const bool ShiftHalfBin = false;

    const unsigned NumMomBins = 20;
    const double BinWidth = 10;
    const double kMin = 5;//ShiftHalfBin*BinWidth*0.5;
    const double kMax = kMin+double(NumMomBins)*BinWidth;

    CATS MyCat;
    FitKitty = &MyCat;

    double Pars[4] = {0,0,0,1.62};
    CATSparameters cSourcePars(CATSparameters::tSource,1,true);
    cSourcePars.SetParameters(&Pars[3]);
    MyCat.SetUseAnalyticSource(true);
    MyCat.SetAnaSource(GaussSource, cSourcePars);

    MyCat.SetExcludeFailedBins(false);
    MyCat.SetMomBins(NumMomBins,kMin,kMax);

    MyCat.SetNumChannels(2);
    MyCat.SetNumPW(0,1);
    MyCat.SetNumPW(1,1);
    MyCat.SetSpin(0,0);
    MyCat.SetSpin(1,1);
    MyCat.SetChannelWeight(0, 1./4.);
    MyCat.SetChannelWeight(1, 3./4.);
    MyCat.SetShortRangePotential(0,0,fDlmPot,pLamPars1S0);
    MyCat.SetShortRangePotential(1,0,fDlmPot,pLamPars3S1);
    MyCat.SetQ1Q2(0);
    MyCat.SetPdgId(2212, 3122);
    MyCat.SetRedMass( (Mass_p*Mass_L)/(Mass_p+Mass_L) );

    complex<double>**** WaveFunctionU; double*** PhaseShifts; double* RadBins; unsigned NumRadBins=0;
    InitHaidenbauerNLO("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/Haidenbauer/pLambdaNLO/",MyCat,&WaveFunctionU,&PhaseShifts,&RadBins,NumRadBins);

    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        MyCat.UseExternalWaveFunction(uBin,0,0,WaveFunctionU[uBin][0][0], NumRadBins, RadBins, PhaseShifts[uBin][0][0]);
        MyCat.UseExternalWaveFunction(uBin,1,0,WaveFunctionU[uBin][1][0], NumRadBins, RadBins, PhaseShifts[uBin][1][0]);
    }

    MyCat.KillTheCat();
    MyCat.SetNotifications(CATS::nSilent);

    const TString OliFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/OliverPlotter/lambda_histos.root";
    const TString OliHistoName = "hC2_lambdap_fullycorr";
    const double UnitsConv = 1;
    TFile* OliFile = OliFileName!=""?new TFile(OliFileName, "read"):NULL;
    TH1F* OliHisto = OliFile?(TH1F*)OliFile->Get(OliHistoName):NULL;

    TH1F* hOliBinShift = new TH1F("hOliBinShift","hOliBinShift",NumMomBins+1,
                                  (OliHisto->GetBinLowEdge(1)-double(ShiftHalfBin)*OliHisto->GetBinWidth(1)*0.5)*UnitsConv,
                                  //(OliHisto->GetXaxis()->GetBinUpEdge(NumMomBins)-double(ShiftHalfBin)*OliHisto->GetBinWidth(NumMomBins)*0.5)*UnitsConv);
                                  (OliHisto->GetXaxis()->GetBinUpEdge(NumMomBins)+OliHisto->GetBinWidth(NumMomBins)-
                                   double(ShiftHalfBin)*OliHisto->GetBinWidth(NumMomBins)*0.5)*UnitsConv);

    if(hOliBinShift->GetBinWidth(1)!=OliHisto->GetBinWidth(1)*UnitsConv) printf("TROUBLE!\n");
    for(unsigned uBin=1; uBin<=NumMomBins+1; uBin++){
        hOliBinShift->SetBinContent(uBin, OliHisto->GetBinContent(uBin));
        hOliBinShift->SetBinError(uBin, OliHisto->GetBinError(uBin));
    }

    TH1F* hAxis = new TH1F("hAxis","hAxis",NumMomBins,0,kMax);
    hAxis->SetStats(false);
    hAxis->SetTitle("");
    hAxis->GetXaxis()->SetLabelSize(0.065);
    hAxis->GetXaxis()->SetTitle("k (MeV)");
    hAxis->GetXaxis()->CenterTitle();
    hAxis->GetXaxis()->SetTitleOffset(1.1);
    hAxis->GetXaxis()->SetLabelOffset(0.02);
    hAxis->GetXaxis()->SetTitleSize(0.075);
    hAxis->GetYaxis()->SetLabelSize(0.065);
    hAxis->GetYaxis()->SetTitle("C(k)");
    hAxis->GetYaxis()->CenterTitle();
    hAxis->GetYaxis()->SetTitleOffset(0.7);
    hAxis->GetYaxis()->SetTitleSize(0.075);

    hAxis->GetYaxis()->SetRangeUser(0.5, 3.2);
    hAxis->GetYaxis()->SetLimits(0.5, 3.2);

    //hOliBinShift->GetXaxis()->SetNdivisions(504);

    hOliBinShift->SetLineWidth(10);
    hOliBinShift->SetLineColor(kBlack);

    FitKitty = &MyCat;
    MyCat.SetNotifications(CATS::nWarning);

    TF1* MyFit = new TF1("MyFit",this,&CATSpaper::Fitter1,kMin,kMax,4,"CATSpaper","Fitter1");
    //MyFit->SetParameter(0,1);//a
    //MyFit->SetParLimits(0,0.5,1.5);
MyFit->FixParameter(0,1);//a
    MyFit->FixParameter(1,0);//b
    MyFit->FixParameter(2,1);//λ
    //MyFit->SetParameter(3,1.62);//r
    //MyFit->SetParLimits(3,1,3);
    MyFit->FixParameter(3,1.62);//r
//printf("C(%.2f) = %.2f (%.2f)\n", 10., MyCat.EvalCorrFun(10), MyFit->Eval(10));
    TGraph grFit;
    grFit.Set(NumMomBins);
    grFit.SetName("grFit");
    grFit.SetLineWidth(8);
    grFit.SetLineColor(kRed+1);
    grFit.SetMarkerSize(0);
    grFit.SetMarkerColor(kRed+1);
    grFit.SetFillColor(kWhite);

    hOliBinShift->Fit(MyFit, "S, N, R, M");
    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        //grFit.SetPoint(uBin,MyCat.GetMomentum(uBin),MyFit->Eval(MyCat.GetMomentum(uBin)));
        grFit.SetPoint(uBin,MyCat.GetMomentum(uBin),MyCat.GetCorrFun(uBin));
    }

    char* fBuffer = new char [512];

    FILE * pFile;
    pFile = fopen (TString::Format("%sFitHadesData_pL.txt",OutputDir.Data()),"w");
    fprintf (pFile, "%18s","k (MeV)");

    sprintf(fBuffer,"Data");
    fprintf (pFile,"%18s",fBuffer);

    sprintf(fBuffer,"DataErr");
    fprintf (pFile,"%18s",fBuffer);

    sprintf(fBuffer,"CatsFit");
    fprintf (pFile,"%18s",fBuffer);

    double xVal,yVal;
    for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
        fprintf (pFile, "\n%18.3f",MyCat.GetMomentum(uMomBin));
//printf("k=%.3f; h=%.3f\n",MyCat.GetMomentum(uMomBin),hOliBinShift->GetBinCenter(uMomBin+2));

        //the +2 comes from the fact, that the first histo bin sits at k=0, while the first CATS bin is at k=4
        yVal = hOliBinShift->GetBinContent(uMomBin+2);
        fprintf (pFile, "%18.3f",yVal);

        yVal = hOliBinShift->GetBinError(uMomBin+2);
        fprintf (pFile, "%18.3f",yVal);

        grFit.GetPoint(uMomBin,xVal,yVal);
        if(xVal!=MyCat.GetMomentum(uMomBin)) {printf("WTF MAN!\n");}
        fprintf (pFile, "%18.3f",yVal);

    }

    TF1* UnitLine = new TF1("UnitLine", "1", kMin,kMax);
    UnitLine->SetLineWidth(8);
    UnitLine->SetLineColor(kGray+1);

    TLegend* myLegend = new TLegend(0.65,0.75,0.975,0.95);//lbrt
    myLegend->SetName("myLegend");
    myLegend->SetTextSize(0.05);
    myLegend->AddEntry(hOliBinShift, TString::Format("HADES data"));
    myLegend->AddEntry(&grFit, TString::Format("CATS with NLO"));

    TPaveText* myPT = new TPaveText(0.6,0.55,0.975,0.75, "blNDC");//lbrt
    myPT->SetName(TString::Format("myPT"));
    myPT->SetBorderSize(1);
    myPT->SetTextSize(0.05);
    myPT->SetFillColor(kWhite);
    myPT->SetTextFont(22);
    myPT->AddText(TString::Format("R_{HADES} = %.2f #pm %.2f fm", 1.62, 0.01));
    myPT->AddText(TString::Format("#color[%u]{R_{CATS} = %.2f #pm %.2f fm}", unsigned(kRed+1), MyFit->GetParameter(3),MyFit->GetParError(3)));
    //for(unsigned u=0; u<10; u++){
    //    const unsigned FROM=kRed;
    //    myPT->AddText(TString::Format("%u: #color[%u]{R_{fit} = %.2f #pm %.2f fm}", FROM+u, FROM+u, MyFit->GetParameter(3),MyFit->GetParError(3)));
    //}

    TCanvas* cCk_pL = new TCanvas("cCk_pL", "cCk_pL", 1);
    cCk_pL->cd(0); cCk_pL->SetCanvasSize(1920, 1080); cCk_pL->SetMargin(0.12,0.025,0.18,0.05);//lrbt

    hAxis->Draw("axis");
    UnitLine->Draw("same");
    hOliBinShift->Draw("same");
    grFit.Draw("same,C");
    myLegend->Draw("same");
    //myPT->Draw("same");

    cCk_pL->SaveAs(TString::Format("%sFitHadesData_pL.png",OutputDir.Data()));
    hOliBinShift->SetLineWidth(hOliBinShift->GetLineWidth()/2.5);
    UnitLine->SetLineWidth(UnitLine->GetLineWidth()/2.5);
    grFit.SetLineWidth(grFit.GetLineWidth()/2.5);
    cCk_pL->SaveAs(TString::Format("%sFitHadesData_pL.pdf",OutputDir.Data()));
    hOliBinShift->SetLineWidth(hOliBinShift->GetLineWidth()*2.5);
    UnitLine->SetLineWidth(UnitLine->GetLineWidth()*2.5);
    grFit.SetLineWidth(grFit.GetLineWidth()*2.5);

    delete MyFit;
    delete hAxis;
    delete hOliBinShift;
    delete UnitLine;
    delete OliFile;
    delete myLegend;
    delete myPT;
    delete cCk_pL;

    delete [] fBuffer;

    CleanHaidenbauer(MyCat,&WaveFunctionU,&PhaseShifts,&RadBins);

}


//void CATSpaper::PlotCk_pp(){
/*
    const double Plot_kMin = 0;
    const double Plot_kMax = 200;
    const double NumPlotPts = 50;
    //CATShisto<double> CatHist(NumPlotPts,Plot_kMin,Plot_kMax);

    CATS MyCatAV18;
    MyCatAV18.SetUseAnalyticSource(true);
    MyCatAV18.SetAnaSource(GaussSource, ppSourcePars);

    MyCatAV18.SetExcludeFailedBins(false);
    MyCatAV18.SetMomBins(NumPlotPts,Plot_kMin,Plot_kMax);

    MyCatAV18.SetNumChannels(4);
    MyCatAV18.SetNumPW(0,2);
    MyCatAV18.SetNumPW(1,2);
    MyCatAV18.SetNumPW(2,2);
    MyCatAV18.SetNumPW(3,2);
    MyCatAV18.SetSpin(0,0);
    MyCatAV18.SetSpin(1,1);
    MyCatAV18.SetSpin(2,1);
    MyCatAV18.SetSpin(3,1);
    MyCatAV18.SetChannelWeight(0, 3./12.);
    MyCatAV18.SetChannelWeight(1, 1./12.);
    MyCatAV18.SetChannelWeight(2, 3./12.);
    MyCatAV18.SetChannelWeight(3, 5./12.);
    MyCatAV18.SetShortRangePotential(0,0,fDlmPot,ppPars1S0);
    MyCatAV18.SetShortRangePotential(1,1,fDlmPot,ppPars3P0);
    MyCatAV18.SetShortRangePotential(2,1,fDlmPot,ppPars3P1);
    MyCatAV18.SetShortRangePotential(3,1,fDlmPot,ppPars3P2);
    MyCatAV18.SetQ1Q2(1);
    MyCatAV18.SetPdgId(2212, 2212);
    MyCatAV18.SetRedMass( 0.5*Mass_p );

    ppSourcePars[3] = 1.2;
    MyCat.KillTheCat();
*/
//}

double CATSpaper::Fitter1(double* mom, double* par){
    FitKitty->SetAnaSource(0, par[3], true);
    FitKitty->KillTheCat();
    return (par[0]+par[1]**mom)*(par[2]*FitKitty->EvalCorrFun(*mom)+1-par[2]);
}

//pLab to pCm, Mass2 is the mass of the particle at rest
double CATSpaper::pLab_pCm(const double& pLab, const double& Mass1, const double& Mass2){
    return sqrt(pow(Mass2*pLab,2)/(pow(Mass2,2)+pow(Mass1,2)+2*Mass2*sqrt(pow(Mass1,2)+pow(pLab,2))));
}

//tLab to kCm, Mass2 is the mass of the particle at rest
double CATSpaper::tLab_pCm(const double& tLab, const double& Mass1, const double& Mass2){
    return pLab_pCm(sqrt(tLab*tLab+2*tLab*Mass1), Mass1, Mass2);
}

double CATSpaper::pCm_pLab(const double& pCm, const double& Mass1, const double& Mass2){
    return pCm/Mass2*(sqrt(Mass2*Mass2+pCm*pCm)+sqrt(Mass1*Mass1+pCm*pCm));
}

double CATSpaper::pCm_tLab(const double& pCm, const double& Mass1, const double& Mass2){
    return sqrt(Mass1*Mass1+pow(pCm_pLab(pCm, Mass1, Mass2),2))-Mass1;
}
