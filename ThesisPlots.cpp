#include "ThesisPlots.h"

#include <stdlib.h>

#include "CATS.h"
#include "CATSconstants.h"
#include "DLM_Potentials.h"
#include "DLM_Source.h"
#include "DLM_SubPads.h"
#include "DLM_CkDecomposition.h"
#include "DLM_Ck.h"
#include "DLM_CkDecomp.h"
#include "DLM_RootWrapper.h"
#include "DLM_CkModels.h"
#include "CommonAnaFunctions.h"

#include "TString.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TPaveText.h"


void SetStyleAxis_1(TH1F* hAxis){
    hAxis->SetStats(false);
    hAxis->SetTitle("");

    hAxis->GetXaxis()->SetTitleSize(0.07);
    hAxis->GetXaxis()->SetLabelSize(0.07);
    hAxis->GetXaxis()->SetTitleOffset(1.3);
    hAxis->GetXaxis()->SetLabelOffset(0.02);

    hAxis->GetYaxis()->SetTitleSize(0.07);
    hAxis->GetYaxis()->SetLabelSize(0.07);
    hAxis->GetYaxis()->SetTitleOffset(1.00);
    hAxis->GetYaxis()->SetNdivisions(505);
}

void SetStyleAxis_2(TH1F* hAxis){
    hAxis->SetStats(false);
    hAxis->SetTitle("");

    hAxis->GetXaxis()->SetTitleSize(0.06);
    hAxis->GetXaxis()->SetLabelSize(0.06);
    hAxis->GetXaxis()->SetTitleOffset(1.3);
    hAxis->GetXaxis()->SetLabelOffset(0.02);

    hAxis->GetYaxis()->SetTitleSize(0.06);
    hAxis->GetYaxis()->SetLabelSize(0.06);
    hAxis->GetYaxis()->SetTitleOffset(0.95);
    hAxis->GetYaxis()->SetNdivisions(505);
}

void SetStyleAxis2D_2(TH2F* hAxis){
    hAxis->SetStats(false);
    hAxis->SetTitle("");

    hAxis->GetXaxis()->SetTitleSize(0.06);
    hAxis->GetXaxis()->SetLabelSize(0.06);
    hAxis->GetXaxis()->SetTitleOffset(1.3);
    hAxis->GetXaxis()->SetLabelOffset(0.02);

    hAxis->GetYaxis()->SetTitleSize(0.06);
    hAxis->GetYaxis()->SetLabelSize(0.06);
    hAxis->GetYaxis()->SetTitleOffset(0.95);
    hAxis->GetYaxis()->SetNdivisions(505);

    hAxis->GetZaxis()->SetTitleSize(0.05);
    hAxis->GetZaxis()->SetLabelSize(0.05);
    hAxis->GetZaxis()->SetTitleOffset(1.25);

}

void SetStyleAxis_3(TH1F* hAxis){
    hAxis->SetStats(false);
    hAxis->SetTitle("");

    hAxis->GetXaxis()->SetTitleSize(0.08);
    hAxis->GetXaxis()->SetLabelSize(0.08);
    hAxis->GetXaxis()->SetTitleOffset(1.3);
    hAxis->GetXaxis()->SetLabelOffset(0.02);

    hAxis->GetYaxis()->SetTitleSize(0.08);
    hAxis->GetYaxis()->SetLabelSize(0.08);
    hAxis->GetYaxis()->SetTitleOffset(0.5);
    hAxis->GetYaxis()->SetNdivisions(505);
}

void ToyPotentials_PS_WF_CF(){

    printf("Start: ToyPotentials_PS_WF_CF\n");

    const TString OutputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/PhD_Thesis/Pictures/InteractionExamples/";

    const double MassX = 1116;
    const double MassY = 1116;
    const double MassXY = (MassX*MassY)/(MassX+MassY);

    const double kMin = 0;
    const double kMax = 200;
    const unsigned NumMomBins = 100;

    const double dfit_Min = 2;
    const double dfit_Max = 100;

    const double rMin = 0;
    const double rMax = 48;
    const unsigned NumRadBins = 512;
    const double rWidth = (rMax-rMin)/double(NumRadBins);

    const unsigned NumRad = 3;
    double* RADIUS = new double[NumRad];
    //RADIUS[0] = 0.8;
    //RADIUS[1] = 1.2;
    //RADIUS[2] = 1.6;
    //RADIUS[3] = 3.0;
    RADIUS[0] = 1.0;
    RADIUS[1] = 1.5;
    RADIUS[2] = 4.0;

    const unsigned NumFreeWF = 3;
    double* FreeWaveK = new double[NumFreeWF];
    FreeWaveK[0] = 50;
    FreeWaveK[1] = 100;
    FreeWaveK[2] = 200;

    const unsigned NumPot = 3;
    TString* PotNames = new TString[NumPot];
    PotNames[0] = "Toy1_X";
    PotNames[1] = "ND54_X";
    PotNames[2] = "ND46_X";

    TString PotLegend[NumPot];
    PotLegend[0] = "V_{#minus}";
    PotLegend[1] = "V_{+}";
    PotLegend[2] = "V_{B}";

    int* PotColors = new int[NumPot];
    PotColors[0] = kRed+1;
    PotColors[1] = kGreen+2;
    PotColors[2] = kAzure+10;


    int* FwfColors = new int[NumFreeWF];
    FwfColors[0] = kMagenta+1;
    FwfColors[1] = kCyan+1;
    FwfColors[2] = kYellow+2;


    TGraph** gCk = new TGraph* [NumPot];
    for(unsigned uPot=0; uPot<NumPot; uPot++){
        gCk[uPot] = new TGraph [NumRad];
        for(unsigned uSor=0; uSor<NumRad; uSor++){
            gCk[uPot][uSor].SetName("gCk_"+PotNames[uPot]+TString::Format("_%.1f",RADIUS[uSor]));
        }
    }

    TF1** fPot = new TF1* [NumPot];
    for(unsigned uPot=0; uPot<NumPot; uPot++){
        //V0*exp(-r^2/β0^2)+V1*exp(-r^2/β1^2)+V2*exp(-r^2/β2^2)
        fPot[uPot] = new TF1("fPot_"+PotNames[uPot],"[0]*exp(-pow(x/[1],2.))+[2]*exp(-pow(x/[3],2.))+[4]*exp(-pow(x/[5],2.))",0,6);
    }

    TGraph* gWFref = new TGraph [NumPot];
    TGraph* gWFas = new TGraph [NumPot];
    TGraph* gWF = new TGraph [NumPot];
    TGraph* grWFref = new TGraph [NumPot];
    TGraph* grWFas = new TGraph [NumPot];
    TGraph* grWF = new TGraph [NumPot];
    TGraph* gWFd = new TGraph [NumPot];//diff to ref
    TGraph* gWFda = new TGraph [NumPot];//diff to as
    TGraph* gTotWF = new TGraph [NumPot];

    for(unsigned uPot=0; uPot<NumPot; uPot++){
        gWFref[uPot].SetName("gWFref_"+PotNames[uPot]);
        gWFas[uPot].SetName("gWFas_"+PotNames[uPot]);
        gWF[uPot].SetName("gWF_"+PotNames[uPot]);
        grWFref[uPot].SetName("grWFref_"+PotNames[uPot]);
        grWFas[uPot].SetName("grWFas_"+PotNames[uPot]);
        grWF[uPot].SetName("grWF_"+PotNames[uPot]);
        gWFd[uPot].SetName("gWFd_"+PotNames[uPot]);
        gWFda[uPot].SetName("gWFda_"+PotNames[uPot]);
        gTotWF[uPot].SetName("gTotWF_"+PotNames[uPot]);
    }

    TGraph* gFreeWF = new TGraph [NumFreeWF];
    for(unsigned uFWF=0; uFWF<NumFreeWF; uFWF++){
        gFreeWF[uFWF].SetName(TString::Format("gFreeWF_%.0f",FreeWaveK[uFWF]));
    }

    TGraph* gSource = new TGraph [NumRad];
    for(unsigned uSor=0; uSor<NumRad; uSor++){
        gSource[uSor].SetName("gSource_"+PotNames[uSor]);
    }


    TH1F** hPhaseShiftsPar = new TH1F* [NumPot];
    TH1F** hPhaseShifts = new TH1F* [NumPot];
    TF1** fPhaseShiftsPar = new TF1* [NumPot];
    TF1** fPhaseShifts = new TF1* [NumPot];

    for(unsigned uPot=0; uPot<NumPot; uPot++){
        hPhaseShiftsPar[uPot] = new TH1F("hPhaseShiftsPar_"+PotNames[uPot],"hPhaseShiftsPar_"+PotNames[uPot],NumMomBins,kMin,kMax);
        hPhaseShifts[uPot] = new TH1F("hPhaseShifts_"+PotNames[uPot],"hPhaseShifts_"+PotNames[uPot],NumMomBins,kMin,kMax);
        fPhaseShiftsPar[uPot] = new TF1("fPhaseShiftsPar_"+PotNames[uPot],"[2]*x*x*x*x+0.5*[1]/197.327*x*x+197.327/[0]", dfit_Min, dfit_Max);
        if(uPot==0) fPhaseShifts[uPot] = new TF1("fPhaseShifts_"+PotNames[uPot],"57.2957795*TMath::ATan2(1./([2]*x*x*x+0.5*[1]/197.327*x+197.327/[0]/x),1.)", dfit_Min, dfit_Max);
        else if(uPot==1) fPhaseShifts[uPot] = new TF1("fPhaseShifts_"+PotNames[uPot],"57.2957795*TMath::ATan2(1./([2]*x*x*x+0.5*[1]/197.327*x+197.327/[0]/x),1.)", dfit_Min, dfit_Max);
        else fPhaseShifts[uPot] = new TF1("fPhaseShifts_"+PotNames[uPot],"-57.2957795*TMath::ATan2(1./([2]*x*x*x+0.5*[1]/197.327*x+197.327/[0]/x),-1.)", dfit_Min, dfit_Max);
    }

    CATSparameters ParsSource(CATSparameters::tSource,1,true);
    ParsSource.SetParameter(0,1.2);
    CATSparameters ParsPot(CATSparameters::tPotential,6,true);

    CATS Kitty;
    Kitty.SetMomBins(NumMomBins,kMin,kMax);
    Kitty.SetThetaDependentSource(false);
    Kitty.SetAnaSource(GaussSource, ParsSource);
    Kitty.SetUseAnalyticSource(true);
    Kitty.SetMomentumDependentSource(false);
    Kitty.SetExcludeFailedBins(false);
    Kitty.SetQ1Q2(0);
    Kitty.SetQuantumStatistics(false);
    Kitty.SetPdgId(1,2);
    Kitty.SetRedMass( MassXY );
    Kitty.SetNumChannels(1);
    Kitty.SetNumPW(0,1);
    Kitty.SetSpin(0,0);
    Kitty.SetChannelWeight(0, 1.);
    Kitty.SetEpsilonConv(2e-9);
    Kitty.SetEpsilonProp(2e-9);
    Kitty.SetMaxNumThreads(4);
    for(unsigned uPot=0; uPot<NumPot; uPot++){
        switch(uPot){
        //Toy1_X (modified)
        case 0 :    //ParsPot.SetParameter(0,-144.5);
                    //ParsPot.SetParameter(1,2.11);
                    //ParsPot.SetParameter(2,520.0);
                    //ParsPot.SetParameter(3,0.54);
                    ParsPot.SetParameter(0,-150);
                    ParsPot.SetParameter(1,1.0);
                    ParsPot.SetParameter(2,2000);
                    ParsPot.SetParameter(3,0.7);
                    ParsPot.SetParameter(4,0.0);
                    ParsPot.SetParameter(5,100.0);
                    break;
        //ND54_X (modified)
        case 1 :    ParsPot.SetParameter(0,-150);
                    ParsPot.SetParameter(1,1.0);
                    ParsPot.SetParameter(2,1000);
                    ParsPot.SetParameter(3,0.5);
                    ParsPot.SetParameter(4,0.0);
                    ParsPot.SetParameter(5,100.0);
                    break;
        //ND46_X (modified)
        case 2 :    ParsPot.SetParameter(0,-150);
                    ParsPot.SetParameter(1,1.0);
                    ParsPot.SetParameter(2,125);
                    ParsPot.SetParameter(3,0.5);
                    ParsPot.SetParameter(4,0.0);
                    ParsPot.SetParameter(5,100.0);
                    break;
        default : break;
        }
        Kitty.RemoveShortRangePotential(0,0);
        Kitty.SetShortRangePotential(0,0,TripleGaussSum,ParsPot);
        for(unsigned uSor=0; uSor<NumRad; uSor++){
            Kitty.SetAnaSource(0,RADIUS[uSor]);
            Kitty.KillTheCat();
            for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
                gCk[uPot][uSor].SetPoint(uMomBin,Kitty.GetMomentum(uMomBin),Kitty.GetCorrFun(uMomBin));
                if(uSor==0){
                    hPhaseShiftsPar[uPot]->SetBinContent(uMomBin+1, Kitty.GetMomentum(uMomBin)/tan(Kitty.GetPhaseShift(uMomBin,0,0)));
                    hPhaseShifts[uPot]->SetBinContent(uMomBin+1, Kitty.GetPhaseShift(uMomBin,0,0)*RadToDeg);
                }
            }

            if(uSor==0){
                switch(uPot){
                //Toy1
                case 0 :    fPhaseShiftsPar[uPot]->SetParameter(0, -0.73);
                            fPhaseShiftsPar[uPot]->SetParameter(1, 7.72);
                            fPhaseShifts[uPot]->SetParameter(0, -0.73);
                            fPhaseShifts[uPot]->SetParameter(1, 7.72);
                            break;
                //ND54
                case 1 :    fPhaseShiftsPar[uPot]->SetParameter(0, 1.89);
                            fPhaseShiftsPar[uPot]->SetParameter(1, 3.39);
                            fPhaseShifts[uPot]->SetParameter(0, 1.89);
                            fPhaseShifts[uPot]->SetParameter(1, 3.39);
                            break;
                //ND46
                case 2 :    fPhaseShiftsPar[uPot]->SetParameter(0, -4.62);
                            fPhaseShiftsPar[uPot]->SetParameter(1, 1.3);
                            fPhaseShifts[uPot]->SetParameter(0, -4.72);
                            fPhaseShifts[uPot]->SetParameter(1, 1.31);
                            break;
                default : break;
                }
                fPhaseShiftsPar[uPot]->SetParLimits(0, -10., 10.);
                fPhaseShiftsPar[uPot]->SetParLimits(1, 0., 32.);
                fPhaseShiftsPar[uPot]->FixParameter(2, 0);
                hPhaseShiftsPar[uPot]->Fit(fPhaseShiftsPar[uPot], "S, N, R, M");

                fPhaseShifts[uPot]->FixParameter(0, fPhaseShiftsPar[uPot]->GetParameter(0));
                fPhaseShifts[uPot]->FixParameter(1, fPhaseShiftsPar[uPot]->GetParameter(1));
                //fPhaseShifts[uPot]->SetParLimits(0, fPhaseShiftsPar[uPot]->GetParameter(0)*0.9, fPhaseShiftsPar[uPot]->GetParameter(0)*1.1);
                //fPhaseShifts[uPot]->SetParLimits(1, fPhaseShiftsPar[uPot]->GetParameter(1)*0.9, fPhaseShiftsPar[uPot]->GetParameter(1)*1.1);
                fPhaseShifts[uPot]->FixParameter(2, 0);
                hPhaseShifts[uPot]->Fit(fPhaseShifts[uPot], "S, N, R, M");


                fPot[uPot]->SetParameter(0,ParsPot.GetParameter(0));
                fPot[uPot]->SetParameter(1,ParsPot.GetParameter(1));
                fPot[uPot]->SetParameter(2,ParsPot.GetParameter(2));
                fPot[uPot]->SetParameter(3,ParsPot.GetParameter(3));
                fPot[uPot]->SetParameter(4,ParsPot.GetParameter(4));
                fPot[uPot]->SetParameter(5,ParsPot.GetParameter(5));

                for(unsigned uRadBin=0; uRadBin<NumRadBins; uRadBin++){
                    double RAD = rMin+0.5*rWidth+double(uRadBin)*rWidth;
                    gWFref[uPot].SetPoint(uRadBin,RAD,real(Kitty.EvalReferenceRadialWF(Kitty.GetMomBin(50),0,RAD,true)));
                    grWFref[uPot].SetPoint(uRadBin,RAD,real(Kitty.EvalReferenceRadialWF(Kitty.GetMomBin(50),0,RAD,false)));
                    if(uPot==0){
                        for(unsigned uFWF=0; uFWF<NumFreeWF; uFWF++){
                            gFreeWF[uFWF].SetPoint(uRadBin,RAD,real(Kitty.EvalReferenceRadialWF(Kitty.GetMomBin(FreeWaveK[uFWF]),0,RAD,true)));
                        }
                    }
                    gWFas[uPot].SetPoint(uRadBin,RAD,real(Kitty.EvalAsymptoticRadialWF(Kitty.GetMomBin(50),0,0,RAD,true)));
                    grWFas[uPot].SetPoint(uRadBin,RAD,real(Kitty.EvalAsymptoticRadialWF(Kitty.GetMomBin(50),0,0,RAD,false)));
                    gWF[uPot].SetPoint(uRadBin,RAD,real(Kitty.EvalRadialWaveFunction(Kitty.GetMomBin(50),0,0,RAD,true)));
                    grWF[uPot].SetPoint(uRadBin,RAD,real(Kitty.EvalRadialWaveFunction(Kitty.GetMomBin(50),0,0,RAD,false)));
                    gWFd[uPot].SetPoint(uRadBin,RAD,real(Kitty.EvalRadialWaveFunction(Kitty.GetMomBin(50),0,0,RAD,true))-real(Kitty.EvalReferenceRadialWF(Kitty.GetMomBin(50),0,RAD,true)));
                    gWFda[uPot].SetPoint(uRadBin,RAD,fabs(real(Kitty.EvalRadialWaveFunction(Kitty.GetMomBin(50),0,0,RAD,true))-real(Kitty.EvalAsymptoticRadialWF(Kitty.GetMomBin(50),0,0,RAD,true))));
                    gTotWF[uPot].SetPoint(uRadBin,RAD,Kitty.EvalWaveFun2(Kitty.GetMomBin(50),RAD,0));
                }
            }
            if(uPot==0){
                for(unsigned uRadBin=0; uRadBin<NumRadBins; uRadBin++){
                    double RAD = rMin+0.5*rWidth+double(uRadBin)*rWidth;
                    gSource[uSor].SetPoint(uRadBin,RAD,Kitty.EvaluateTheSource(50,RAD,0));
                }
            }
        }
    }

    //POTENTIAL PLOTS
    TH1F* hPot_Dummy = new TH1F("hPot_Dummy","hPot_Dummy",128,0,2.5);
    SetStyleAxis_1(hPot_Dummy);
    hPot_Dummy->GetXaxis()->SetTitle("r (fm)");
    hPot_Dummy->GetYaxis()->SetTitle("V(r) (MeV)");
    hPot_Dummy->GetYaxis()->SetRangeUser(-100,400);

    TLegend* lLegend_Pot = new TLegend(0.7,0.60,0.85,0.85);//lbrt
    lLegend_Pot->SetName(TString::Format("lLegend_Pot"));
    lLegend_Pot->SetTextSize(0.07);
    for(unsigned uPot=0; uPot<NumPot; uPot++){
        lLegend_Pot->AddEntry(fPot[uPot],PotLegend[uPot]);
    }

    TCanvas* cPotentials = new TCanvas("cPotentials", "cPotentials", 1);
    cPotentials->cd(0); cPotentials->SetCanvasSize(1920, 1080); cPotentials->SetMargin(0.15,0.05,0.2,0.05);//lrbt
    cPotentials->SetGrid(true);
    hPot_Dummy->Draw("axis");
    for(unsigned uPot=0; uPot<NumPot; uPot++){
        fPot[uPot]->SetLineWidth(4.-0.33*double(uPot));
        fPot[uPot]->SetLineColor(PotColors[uPot]);
        fPot[uPot]->Draw("same");
    }
    lLegend_Pot->Draw("same");
    cPotentials->SaveAs(OutputFolder+"cPotentials.png");
    cPotentials->SaveAs(OutputFolder+"cPotentials.pdf");


    //FREE S-WAVE
    TH1F* hFWF_Dummy = new TH1F("hFWF_Dummy","hFWF_Dummy",128,0,32);
    SetStyleAxis_1(hFWF_Dummy);
    hFWF_Dummy->GetXaxis()->SetTitle("r (fm)");
    hFWF_Dummy->GetYaxis()->SetTitle("j_{0}(kr) = a_{k,0}(r)/r");
    hFWF_Dummy->GetYaxis()->SetRangeUser(-0.3,1.0);
    TLegend* lLegend_FWF = new TLegend(0.475,0.575,0.925,0.925);//lbrt
    lLegend_FWF->SetName(TString::Format("lLegend_FWF"));
    lLegend_FWF->SetTextSize(0.07);
    for(unsigned uFWF=0; uFWF<NumFreeWF; uFWF++){
        lLegend_FWF->AddEntry(&gFreeWF[uFWF],TString::Format("j_{0}(kr) (k=%.0f MeV)",FreeWaveK[uFWF]));
    }
    TCanvas* cFWF = new TCanvas("cFWF", "cFWF", 1);
    cFWF->cd(0); cFWF->SetCanvasSize(1920, 1080); cFWF->SetMargin(0.15,0.05,0.2,0.05);//lrbt
    cFWF->SetGrid(true);
    hFWF_Dummy->Draw("axis");
    for(unsigned uFWF=0; uFWF<NumFreeWF; uFWF++){
        gFreeWF[uFWF].SetLineWidth(4.-0.33*double(uFWF));
        gFreeWF[uFWF].SetLineColor(FwfColors[uFWF]);
        gFreeWF[uFWF].Draw("C,same");
    }
    lLegend_FWF->Draw("same");
    cFWF->SaveAs(OutputFolder+"cFWF.png");
    cFWF->SaveAs(OutputFolder+"cFWF.pdf");

    //plot the wave functions

    double WF_Mom = Kitty.GetMomentum(Kitty.GetMomBin(50));
    double WFu_yMax = 1.8;
    double WFtot_yMax = 2.4;
    double WF_rMax = 16;
    double WF_rhoMax = WF_rMax*FmToNu*WF_Mom;
    TH1F* hWFu_Dummy = new TH1F("hWFu_Dummy","hWFu_Dummy",128,0,WF_rMax);
    SetStyleAxis_2(hWFu_Dummy);
    hWFu_Dummy->GetXaxis()->SetTitle("r (fm)");
    hWFu_Dummy->GetYaxis()->SetTitle("u_{k,0}(r)/r");
    hWFu_Dummy->GetYaxis()->SetRangeUser(-0.575,WFu_yMax);
    TF1* fWFu_base = new TF1("fWFu_base","0",0,WF_rMax);
    fWFu_base->SetLineColor(kBlack);
    fWFu_base->SetLineWidth(1.5);
    fWFu_base->SetLineStyle(2);

    TH1F* hWFda_Dummy = new TH1F("hWFda_Dummy","hWFda_Dummy",128,0,WF_rMax);
    SetStyleAxis_2(hWFda_Dummy);
    hWFda_Dummy->GetXaxis()->SetTitle("r (fm)");
    hWFda_Dummy->GetYaxis()->SetTitle("[u_{k,0}(r)-a_{k,0}(r)]/r");
    hWFda_Dummy->GetYaxis()->SetRangeUser(-1.0,1.0);

    TH1F* hWFtot_Dummy = new TH1F("hWFtot_Dummy","hWFtot_Dummy",128,0,WF_rMax);
    SetStyleAxis_2(hWFtot_Dummy);
    hWFtot_Dummy->GetXaxis()->SetTitle("r (fm)");
    hWFtot_Dummy->GetYaxis()->SetTitle("|#Psi_{k}(r)|^{2}");
    hWFtot_Dummy->GetYaxis()->SetRangeUser(-0.1,WFtot_yMax);
    TF1* fWFtot_base = new TF1("fWFtot_base","1",0,WF_rMax);
    fWFtot_base->SetLineColor(kBlack);
    fWFtot_base->SetLineWidth(1.5);
    fWFtot_base->SetLineStyle(2);

    TH1F* hSource_Dummy = new TH1F("hSource_Dummy","hSource_Dummy",128,0,WF_rMax);
    SetStyleAxis_2(hSource_Dummy);
    hSource_Dummy->GetXaxis()->SetTitle("r (fm)");
    hSource_Dummy->GetYaxis()->SetTitle("4#pir^{2}S(r) (1/fm)");
    hSource_Dummy->GetYaxis()->SetRangeUser(0,0.45);

    //DLM_SubPads Pad_WF(1080,1920);
    DLM_SubPads Pad_WF(1080,1440);
    Pad_WF.AddSubPad(0,1,0.65,1);
    Pad_WF.AddSubPad(0,1,0.35,0.65);
    Pad_WF.AddSubPad(0,1,0,0.35);

    Pad_WF.SetMargin(0,0.12,0.02,0.0,0.06);
    Pad_WF.SetMargin(1,0.12,0.02,0.,0.);
    Pad_WF.SetMargin(2,0.12,0.02,0.07,0.0);

    Pad_WF.cd(0);

    hWFu_Dummy->Draw("axis");
    fWFu_base->Draw("same");
    gWFref[0].SetLineWidth(3);
    gWFref[0].SetLineColor(kGray+1);
    gWFref[0].SetLineStyle(4);
    grWFref[0].SetLineWidth(3);
    grWFref[0].SetLineColor(kGray+1);
    grWFref[0].SetLineStyle(4);
    gWFref[0].Draw("same");
    for(unsigned uPot=0; uPot<NumPot; uPot++){
        //gWF[uPot].Draw(uPot?"same":"");
        gWF[uPot].SetLineWidth(4.-0.33*double(uPot));
        gWF[uPot].SetLineColor(PotColors[uPot]);
        gWF[uPot].Draw("same");

        grWF[uPot].SetLineWidth(4.-0.33*double(uPot));
        grWF[uPot].SetLineColor(PotColors[uPot]);
        grWF[uPot].Draw("same");

        gWFas[uPot].SetLineWidth(4.-0.33*double(uPot));
        gWFas[uPot].SetLineColor(PotColors[uPot]);
        gWFas[uPot].SetLineStyle(3);
        gWFas[uPot].Draw("same");

        grWFas[uPot].SetLineWidth(4.-0.33*double(uPot));
        grWFas[uPot].SetLineColor(PotColors[uPot]);
        grWFas[uPot].SetLineStyle(3);
        grWFas[uPot].Draw("same");
    }
    //TGaxis *Arho = new TGaxis(0,2,10,2,"f1",510,"-");
    TGaxis* Arho = new TGaxis(0,WFu_yMax,WF_rMax,WFu_yMax,0,WF_rhoMax,505,"-");
    Arho->SetTitle("#rho=kr");
    Arho->SetLabelSize(hWFu_Dummy->GetXaxis()->GetLabelSize());
    Arho->SetTitleSize(hWFu_Dummy->GetXaxis()->GetTitleSize());
    Arho->SetLabelFont(hWFu_Dummy->GetXaxis()->GetLabelFont());
    Arho->SetTitleFont(hWFu_Dummy->GetXaxis()->GetTitleFont());
    Arho->Draw("same");

    TLegend* lLegend_WF_up = new TLegend(0.65,0.55,0.95,0.80);//lbrt
    lLegend_WF_up->SetName(TString::Format("lLegend_WF_up"));
    lLegend_WF_up->SetTextSize(0.06);
    lLegend_WF_up->AddEntry(&gWFref[0],"Free wave");
    TGraph gSFW_Dummy;
    gSFW_Dummy.SetLineWidth(3.5);
    gSFW_Dummy.SetLineStyle(3);
    gSFW_Dummy.SetLineColor(kBlack);
    lLegend_WF_up->AddEntry(&gSFW_Dummy,"Shifted free wave");
    lLegend_WF_up->Draw("same");

    Pad_WF.cd(1);
    hWFtot_Dummy->Draw("axis");
    for(unsigned uPot=0; uPot<NumPot; uPot++){
        //gWFas[uPot].Draw(uPot?"same":"");
        gTotWF[uPot].SetLineWidth(4.-0.33*double(uPot));
        gTotWF[uPot].SetLineColor(PotColors[uPot]);
        fWFtot_base->Draw("same");
        gTotWF[uPot].Draw("same");
    }
    TLegend* lLegend_WF_low = new TLegend(0.8,0.70,0.95,0.95);//lbrt
    lLegend_WF_low->SetName(TString::Format("lLegend_WF_low"));
    lLegend_WF_low->SetTextSize(0.06);
    for(unsigned uPot=0; uPot<NumPot; uPot++){
        lLegend_WF_low->AddEntry(&gTotWF[uPot],PotLegend[uPot]);
    }
    lLegend_WF_low->Draw("same");

    Pad_WF.cd(2);
    hSource_Dummy->Draw("axis");
    for(unsigned uSor=0; uSor<NumRad; uSor++){
        gSource[uSor].SetLineColor(kBlue+3);
        gSource[uSor].SetLineStyle(uSor+1);
        gSource[uSor].SetLineWidth(3.5);
        gSource[uSor].Draw("same");
    }
    TLegend* lLegend_Source = new TLegend(0.6,0.65,0.95,0.95);//lbrt
    lLegend_Source->SetName(TString::Format("lLegend_Source"));
    lLegend_Source->SetTextSize(0.06);
    for(unsigned uSor=0; uSor<NumRad; uSor++){
        lLegend_Source->AddEntry(&gSource[uSor],TString::Format("Source size: %.1f fm",RADIUS[uSor]));
    }

    lLegend_Source->Draw("same");

    Pad_WF.GetCanvas()->SaveAs(OutputFolder+"Pad_WF.png");
    Pad_WF.GetCanvas()->SaveAs(OutputFolder+"Pad_WF.pdf");


    //PHASE SHIFTS
    TH1F* hPS_Dummy = new TH1F("hPS_Dummy","hPS_Dummy",NumMomBins,kMin,kMax);
    SetStyleAxis_1(hPS_Dummy);
    hPS_Dummy->GetXaxis()->SetTitle("k (MeV)");
    hPS_Dummy->GetYaxis()->SetTitle("#delta (deg)");
    hPS_Dummy->GetYaxis()->SetRangeUser(-90,180);
    hPS_Dummy->GetYaxis()->SetNdivisions(510);

    TLegend* lLegend_PS = new TLegend(0.75,0.65,0.90,0.90);//lbrt
    lLegend_PS->SetName(TString::Format("lLegend_PS"));
    lLegend_PS->SetTextSize(0.07);
    for(unsigned uPot=0; uPot<NumPot; uPot++){
        lLegend_PS->AddEntry(fPot[uPot],PotLegend[uPot]);
    }

    TF1* fPS_base = new TF1("fWFu_base","0",kMin,kMax);
    fPS_base->SetLineColor(kBlack);
    fPS_base->SetLineWidth(1.5);
    fPS_base->SetLineStyle(2);

    TCanvas* cPS = new TCanvas("cPS", "cPS", 1);
    cPS->cd(0); cPS->SetCanvasSize(1920, 1080); cPS->SetMargin(0.15,0.05,0.2,0.05);//lrbt
    cPS->SetGrid(true);
    hPS_Dummy->Draw("axis");
    for(unsigned uPot=0; uPot<NumPot; uPot++){
        hPhaseShifts[uPot]->SetLineWidth(4.-0.33*double(uPot));
        hPhaseShifts[uPot]->SetLineColor(PotColors[uPot]);
        hPhaseShifts[uPot]->Draw("C,same");

        fPhaseShifts[uPot]->SetLineWidth(4.-0.33*double(uPot));
        fPhaseShifts[uPot]->SetLineColor(kBlack);
        fPhaseShifts[uPot]->SetLineStyle(4);
        fPhaseShifts[uPot]->Draw("C,same");
    }
    lLegend_PS->Draw("same");
    fPS_base->Draw("same");
    cPS->SaveAs(OutputFolder+"cPS.png");
    cPS->SaveAs(OutputFolder+"cPS.pdf");


    //CORR FUN

    DLM_SubPads Pad_Ck(1280,1280);
    Pad_Ck.AddSubPad(0,1,0.64,1);
    Pad_Ck.AddSubPad(0,1,0.36,0.64);
    Pad_Ck.AddSubPad(0,1,0,0.36);
    Pad_Ck.SetMargin(0,0.12,0.03,0.0,0.06);
    Pad_Ck.SetMargin(1,0.12,0.03,0.,0.);
    Pad_Ck.SetMargin(2,0.12,0.03,0.08,0.0);

    TH1F** hCk_Dummy = new TH1F* [NumRad];
    double CkCan_Scale = 1.286;
    for(unsigned uSor=0; uSor<NumRad; uSor++){
        hCk_Dummy[uSor] = new TH1F(TString::Format("hCk_Dummy_%u",uSor),TString::Format("hCk_Dummy_%u",uSor),NumMomBins,kMin,kMax);
        SetStyleAxis_3(hCk_Dummy[uSor]);
        hCk_Dummy[uSor]->GetXaxis()->SetTitle("k (MeV)");
        hCk_Dummy[uSor]->GetYaxis()->SetTitle("C(k)");
        hCk_Dummy[uSor]->GetYaxis()->SetRangeUser(-0.15,3.15);
        if(uSor==1) hCk_Dummy[uSor]->GetYaxis()->SetTitleOffset(hCk_Dummy[uSor]->GetYaxis()->GetTitleOffset()/CkCan_Scale);
        Pad_Ck.SetLabelSize(uSor, hCk_Dummy[uSor]->GetXaxis(), 14);
        Pad_Ck.SetTitleSize(uSor, hCk_Dummy[uSor]->GetXaxis(), 14);
        Pad_Ck.SetLabelSize(uSor, hCk_Dummy[uSor]->GetYaxis(), 14);
        Pad_Ck.SetTitleSize(uSor, hCk_Dummy[uSor]->GetYaxis(), 14);
        //hCk_Dummy->GetYaxis()->SetNdivisions(510);
    }


    TLegend* lLegend_Ck = new TLegend(0.75,0.65,0.90,0.90);//lbrt
    lLegend_Ck->SetName(TString::Format("lLegend_Ck"));
    lLegend_Ck->SetTextSize(0.07);
    for(unsigned uPot=0; uPot<NumPot; uPot++){
        lLegend_Ck->AddEntry(fPot[uPot],PotLegend[uPot]);
    }

    double PadCk_l = 0.3;
    double PadCk_r = 0.5;
    //double PadCk_b =
    TPaveText** PT_RAD = new TPaveText* [NumRad];
    for(unsigned uSor=0; uSor<NumRad; uSor++){
        PT_RAD[uSor] = new TPaveText(0.35,0.65,0.65,0.85, "blNDC");//lbrt
        PT_RAD[uSor]->SetName(TString::Format("PT_RAD_%u",uSor));
        PT_RAD[uSor]->SetBorderSize(1);
        PT_RAD[uSor]->SetTextSize(0.055);
        PT_RAD[uSor]->SetFillColor(kWhite);
        PT_RAD[uSor]->SetTextFont(22);
        PT_RAD[uSor]->AddText(TString::Format("Source size: %.1f fm",RADIUS[uSor]));
    }

    TF1* fCk_base = new TF1("fCk_base","1",kMin,kMax);
    fCk_base->SetLineColor(kBlack);
    fCk_base->SetLineWidth(1.5);
    fCk_base->SetLineStyle(2);

    for(unsigned uSor=0; uSor<NumRad; uSor++){
        Pad_Ck.cd(uSor);
        hCk_Dummy[uSor]->Draw("axis");
        for(unsigned uPot=0; uPot<NumPot; uPot++){
            gCk[uPot][uSor].SetLineWidth(4.-0.33*double(uPot));
            gCk[uPot][uSor].SetLineColor(PotColors[uPot]);
            //gCk[uPot][uSor].SetLineStyle(uSor+1);
            gCk[uPot][uSor].Draw("C,same");
        }
        if(uSor==0){
            lLegend_Ck->Draw("same");
        }
        fCk_base->Draw("same");
        PT_RAD[uSor]->Draw("same");
    }

    Pad_Ck.GetCanvas()->SaveAs(OutputFolder+"Pad_Ck.png");
    Pad_Ck.GetCanvas()->SaveAs(OutputFolder+"Pad_Ck.pdf");

    TFile fOutput(OutputFolder+"fOutput.root","recreate");
    for(unsigned uPot=0; uPot<NumPot; uPot++){
        for(unsigned uSor=0; uSor<NumRad; uSor++){
            gCk[uPot][uSor].Write();
        }
        hPhaseShifts[uPot]->Write();
        fPhaseShifts[uPot]->Write();
        hPhaseShiftsPar[uPot]->Write();
        fPhaseShiftsPar[uPot]->Write();
        fPot[uPot]->Write();
        gWFref[uPot].Write();
        gWFas[uPot].Write();
        gWF[uPot].Write();
        grWFref[uPot].Write();
        grWFas[uPot].Write();
        grWF[uPot].Write();
        gWFd[uPot].Write();
        gWFda[uPot].Write();
        gTotWF[uPot].Write();
    }

    for(unsigned uFWF=0; uFWF<NumFreeWF; uFWF++){
        gFreeWF[uFWF].Write();
    }

    for(unsigned uSor=0; uSor<NumRad; uSor++){
        gSource[uSor].Write();
    }



    for(unsigned uPot=0; uPot<NumPot; uPot++){
        delete hPhaseShifts[uPot];
        delete hPhaseShiftsPar[uPot];
        delete fPhaseShiftsPar[uPot];
        delete fPot[uPot];
        delete [] gCk[uPot];
    }
    delete [] hPhaseShifts;
    delete [] hPhaseShiftsPar;
    delete [] fPhaseShiftsPar;
    delete [] fPhaseShifts;
    delete [] gCk;
    delete [] fPot;
    delete [] gWFref;
    delete [] gWFas;
    delete [] gWF;
    delete [] grWFref;
    delete [] grWFas;
    delete [] grWF;
    delete [] gWFd;
    delete [] gWFda;
    delete [] gTotWF;
    delete [] gFreeWF;
    delete [] gSource;

    delete hPot_Dummy;
    delete lLegend_Pot;
    delete cPotentials;

    delete hFWF_Dummy;
    delete lLegend_FWF;
    delete cFWF;

    delete hPS_Dummy;
    delete lLegend_PS;
    delete fPS_base;
    delete cPS;

    delete hSource_Dummy;
    delete lLegend_Source;

    for(unsigned uSor=0; uSor<NumRad; uSor++){
        delete PT_RAD[uSor];
        delete hCk_Dummy[uSor];
    }
    delete [] PT_RAD;
    delete [] hCk_Dummy;
    delete lLegend_Ck;
    delete fCk_base;
    //delete cCk;

    delete hWFu_Dummy;
    delete hWFtot_Dummy;
    delete hWFda_Dummy;
    delete fWFu_base;
    delete fWFtot_base;
}



//pi pi plots for the appendix

//pipi interaction (CATS)
//ExpSource==true => use exponential source
//ExpSource==false => use Gaussian source
TGraph* Basics_PiPiCATS(const bool& ExpSource){

    //constants related to the number of bins and range (in MeV)
    const unsigned NumMomBins = 300;
    const double kMin = 0;
    const double kMax = 300;

    //initialize the CATS object
    CATS PionKitty;
    //(#bins,min,max) or (#bins,BinRangeArray[] as in ROOT)
    PionKitty.SetMomBins(NumMomBins,kMin,kMax);

    //object for the parameters to be used by the source function
    CATSparameters SOURCE_PARS(CATSparameters::tSource,1,true);
    //set the first and only par (source size)
    SOURCE_PARS.SetParameter(0,1.5);
    //say to CATS which Source function to use, and with which parameter set
    PionKitty.SetAnaSource(ExpSource?ExponentialSource:GaussSource,SOURCE_PARS);
    //this step is needed, to set up CATS for a source from a function
    //there is one more possibility, namely to use the source as a
    //direct input from a transport model. This method will not be discussed here
    PionKitty.SetUseAnalyticSource(true);
    //reduces CPU time, for a source that does not has a momentum or an angular dependence
    PionKitty.SetMomentumDependentSource(false);
    PionKitty.SetThetaDependentSource(false);
    //if true, the source is automatically renormalized in the range 0-64 fm. Nice to dummy proof the source, but problematic for sources with large tails
    //for the Gaussian example above, both options should be completely identical
    PionKitty.SetAutoNormSource(false);

    //standard settings for a CATS object which has no strong interaction potential included
    PionKitty.SetNumChannels(1);
    //#which channel, how many PWs
    PionKitty.SetNumPW(0,0);
    //which channel, spin value
    PionKitty.SetSpin(0,0);
    //which channel, weight
    PionKitty.SetChannelWeight(0, 1);

    //include the coulomb interaction. Q1Q2 is the multiplied charge numbers of the two particles
    PionKitty.SetQ1Q2(1);
    //the reduced mass of the two pions
    PionKitty.SetRedMass( 0.5*Mass_pic );
    //for identical particles, we need to include the quantum statistics
    PionKitty.SetQuantumStatistics(true);
    //evaluate the correlation function
    PionKitty.KillTheCat();

    //save the result in a TGraph (ROOT)
    TGraph* grCk = new TGraph();
    grCk->SetName("grCk");
    grCk->Set(NumMomBins);
    //iterate over all momentum bins
    for(unsigned uMom=0; uMom<NumMomBins; uMom++){
        //the GetMomentum() function sets the x-axis
        //the GetCorrFun() function evaluates C(k) (y-axis)
        grCk->SetPoint(uMom,PionKitty.GetMomentum(uMom),PionKitty.GetCorrFun(uMom));
    }

    return grCk;
}

//computes the expectation based on quantum statistics only
//the formula used is from chapter 4.1 in Phys. Rev. C 96 (2017) 064908 (ATLAS paper on pipi correlations in p-Pb)
TGraph* Basics_PiPiTheory(const bool& ExpSource){
    //load from a Mathematica output file
    FILE *InFileCBE;
    const TString CBEname = ExpSource?"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/PhD_Thesis/Mathematica/tab_txCBE_Exp.dat":
        "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/PhD_Thesis/Mathematica/tab_txCBE.dat";
    InFileCBE = fopen(CBEname.Data(), "r");
    if(!InFileCBE){
        printf("\033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", CBEname.Data());
        return NULL;
    }
    char*  cdummy = new char [512];
    float kstar;
    float ckval;
    TGraph gCBE;
    gCBE.SetName("gCBE");
    unsigned NumPointsCBE=0;
    while(!feof(InFileCBE)){
        if(!fgets(cdummy, 511, InFileCBE)) continue;
        sscanf(cdummy, "%f %f",&kstar, &ckval);
        gCBE.SetPoint(NumPointsCBE,kstar,ckval);
        NumPointsCBE++;
    }
    fclose(InFileCBE);

    FILE *InFileK;
    const TString Kname = ExpSource?"/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/PhD_Thesis/Mathematica/tab_txCoulombSame_Exp.dat":
        "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/PhD_Thesis/Mathematica/tab_txCoulombSame.dat";
    InFileK = fopen(Kname.Data(), "r");
    TGraph gK;
    gK.SetName("gK");
    unsigned NumPointsK=0;
    if(!InFileK){
        printf("\033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", Kname.Data());
        return NULL;
    }
    while(!feof(InFileK)){
        if(!fgets(cdummy, 511, InFileK)) continue;
        sscanf(cdummy, "%f %f",&kstar, &ckval);
        gK.SetPoint(NumPointsK,kstar,ckval);
        NumPointsK++;
    }
    fclose(InFileK);

    TGraph* gCk = new TGraph();
    gCk->SetName("gPiPiTheory");
    double dkstar,dckval;
    for(unsigned uBin=0; uBin<NumPointsCBE; uBin++){
        gCBE.GetPoint(uBin,dkstar,dckval);
        gCk->SetPoint(uBin,dkstar,dckval*gK.Eval(dkstar));
    }

    delete [] cdummy;
    return gCk;
}

void ComparePionPion(const bool& ExpSource){
    TGraph* grCATS = Basics_PiPiCATS(ExpSource);
    TGraph* grTheory = Basics_PiPiTheory(ExpSource);

    grTheory->SetLineColor(kRed+1);
    grTheory->SetLineStyle(2);
    grTheory->SetLineWidth(4);

    grCATS->SetLineColor(kBlue+1);
    grCATS->SetLineWidth(2.5);

    TH1F* hAxis = new TH1F("hAxis","hAxis",128,0,300);
    SetStyleAxis_1(hAxis);
    hAxis->GetYaxis()->SetRangeUser(0.9, 1.85);
    hAxis->GetXaxis()->SetTitle("k* (MeV)");
    hAxis->GetYaxis()->SetTitle("C(k*)");

    TF1* UnitLine = new TF1("UnitLine","1",-20,320);
    UnitLine->SetLineColor(kBlack);
    UnitLine->SetLineWidth(2);
    UnitLine->SetLineStyle(1);

    TPaveText* PTinfo = new TPaveText(0.25,0.85,0.5,0.95, "blNDC");//lbrt
    PTinfo->SetName("PTinfo");
    PTinfo->SetBorderSize(1);
    PTinfo->SetTextSize(0.05);
    PTinfo->SetFillColor(kWhite);
    PTinfo->SetTextFont(22);
    PTinfo->AddText(ExpSource?"Cauchy source":"Gaussian source");

    TLegend *legend = new TLegend(0.5,0.8,0.95,0.95);
    legend->SetTextFont(42);
    legend->SetTextSize(0.05);
    legend->AddEntry(grTheory, "Bowler-Sinyukov (approx)");
    legend->AddEntry(grCATS, "CATS (exact)");

    TCanvas* cPiPi = new TCanvas("cPiPi", "cPiPi", 1);
    cPiPi->cd(0);
    cPiPi->SetCanvasSize(1920, 1080);
    cPiPi->SetMargin(0.15,0.05,0.2,0.05);//lrbt
    cPiPi->SetGrid();

    hAxis->Draw();
    UnitLine->Draw("same");
    grTheory->Draw("C,same");
    grCATS->Draw("C,same");
    legend->Draw("same");
    PTinfo->Draw("same");

    TString OutFileName = "ComparePionPion";
    OutFileName += ExpSource?"_Exp":"_Gauss";

    TFile* OutputFile = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/PhD_Thesis/Pictures/PionPion/"+OutFileName+".root","recreate");
    grCATS->Write();
    grTheory->Write();
    cPiPi->Write();

    cPiPi->SaveAs("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/PhD_Thesis/Pictures/PionPion/"+OutFileName+".pdf");

    delete grCATS;
    delete grTheory;
    delete UnitLine;
    delete hAxis;
    delete legend;
    delete PTinfo;
    delete cPiPi;
    delete OutputFile;
}

void Plot_ProtonLambda(TGraph* grCk){

    grCk->SetLineColor(kRed+1);
    grCk->SetLineWidth(4);

    TH1F* hAxis = new TH1F("hAxis","hAxis",128,0,300);
    SetStyleAxis_1(hAxis);
    hAxis->GetYaxis()->SetRangeUser(0.9, 3.35);
    hAxis->GetXaxis()->SetTitle("k* (MeV)");
    hAxis->GetYaxis()->SetTitle("C(k*)");

    TF1* UnitLine = new TF1("UnitLine","1",-20,320);
    UnitLine->SetLineColor(kBlack);
    UnitLine->SetLineWidth(2.0);
    UnitLine->SetLineStyle(1);

    TPaveText* PTinfo = new TPaveText(0.25,0.85,0.5,0.95, "blNDC");//lbrt
    PTinfo->SetName("PTinfo");
    PTinfo->SetBorderSize(1);
    PTinfo->SetTextSize(0.05);
    PTinfo->SetFillColor(kWhite);
    PTinfo->SetTextFont(22);
    PTinfo->AddText("Gaussian source");

    TLegend *legend = new TLegend(0.5,0.8,0.95,0.95);
    legend->SetTextFont(42);
    legend->SetTextSize(0.05);
    legend->AddEntry(grCk, "p#minus#Lambda (Usmani)");

    TCanvas* cpL = new TCanvas("cpL", "cpL", 10, 10, 1920, 1080);
    //cpL->SetCanvasSize(1920, 1080);
    cpL->SetMargin(0.15,0.05,0.2,0.05);//lrbt
    cpL->SetGrid(1,1);
    //cpL->cd(0);

    hAxis->Draw();
    UnitLine->Draw("same");
    grCk->Draw("C,same");
    legend->Draw("same");
    PTinfo->Draw("same");


    //INLET
    TH1F* hAxis_Inlet = (TH1F*)hAxis->Clone("hAxis_Inlet");
    hAxis_Inlet->GetXaxis()->SetTitleSize(0);
    hAxis_Inlet->GetXaxis()->SetLabelSize(0);
    hAxis_Inlet->GetXaxis()->SetTitleOffset(0);
    hAxis_Inlet->GetXaxis()->SetLabelOffset(0);
    hAxis_Inlet->GetXaxis()->SetRangeUser(100, 300);
    hAxis_Inlet->GetXaxis()->SetNdivisions(505);
    hAxis_Inlet->GetXaxis()->SetTitle("");

    hAxis_Inlet->GetYaxis()->SetTitleSize(hAxis->GetYaxis()->GetTitleSize()*1.75);
    hAxis_Inlet->GetYaxis()->SetLabelSize(hAxis->GetYaxis()->GetLabelSize()*1.75);
    hAxis_Inlet->GetYaxis()->SetTitleOffset(hAxis->GetYaxis()->GetTitleOffset()*0.5);
    hAxis_Inlet->GetYaxis()->SetRangeUser(0.965, 1.015);

    const double fXMinInlet=0.190;
    const double fYMinInlet=0.26;
    const double fXMaxInlet=0.989;
    const double fYMaxInlet=0.70;
    TPad *inset_pad = new TPad("inset_pad", "inset_pad", fXMinInlet, fYMinInlet, fXMaxInlet, fYMaxInlet);
    inset_pad->SetTopMargin(0.01);
    inset_pad->SetRightMargin(0.05);
    inset_pad->SetBottomMargin(0.28);
    inset_pad->SetLeftMargin(0.28);
    inset_pad->SetFillStyle(4000);
    inset_pad->Draw();
    inset_pad->cd();
    inset_pad->SetGrid();
    hAxis_Inlet->Draw();
    UnitLine->Draw("same");
    grCk->Draw("C,same");



    TString OutFileName = "Usmani_Gauss";

    TFile* OutputFile = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/PhD_Thesis/Pictures/CATS_examples/"+OutFileName+".root","recreate");
    grCk->Write();
    cpL->Write();

    cpL->SaveAs("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/PhD_Thesis/Pictures/CATS_examples/"+OutFileName+".pdf");

    delete grCk;
    delete UnitLine;
    delete hAxis;
    delete legend;
    delete PTinfo;
    delete cpL;
    delete OutputFile;
}

//initialize the Gauss source for (any) CATS object
void CATS_GaussSource(CATS& Kitty, const double& SourceSize){
    CATSparameters cPars(CATSparameters::tSource,1,true);
    cPars.SetParameter(0,SourceSize);
    Kitty.SetAnaSource(GaussSource, cPars);
    Kitty.SetUseAnalyticSource(true);
    Kitty.SetMomentumDependentSource(false);
    Kitty.SetThetaDependentSource(false);
    //if true, the source is automatically renormalized in the range 0-64 fm.
    //Nice to dummy proof the source, but problematic for sources with large tails
    //for the Gaussian example above, both options should be completely identical
    Kitty.SetAutoNormSource(false);
}

//Basic initialization of a CATS object for pp
void CATS_pp_Basic(CATS& Kitty,
        const unsigned& NumMomBins, const double& kMin, const double& kMax){
    Kitty.SetMomBins(NumMomBins,kMin,kMax);
    Kitty.SetExcludeFailedBins(false);
    Kitty.SetQ1Q2(1);
    Kitty.SetQuantumStatistics(true);
    Kitty.SetRedMass( 0.5*Mass_p );
}
//initialize the interaction for pp, using the AV18 potential.
//the s,l and d waves are included
void CATS_pp_AV18(CATS& Kitty, const bool& pwaves, const bool& dwaves){
    //the 4 channels for pp are:
    //s=0: 1S0 + 3D1
    //s=1: 3P0
    //s=1: 3P1
    //s=1: 3P2
    //note that for s=0 the p-waves are Pauli blocked,
    //for s=1 these are the s and d waves
    if(pwaves){
        Kitty.SetNumChannels(4);
        if(dwaves) Kitty.SetNumPW(0,3);
        else Kitty.SetNumPW(0,1);
        Kitty.SetNumPW(1,2);
        Kitty.SetNumPW(2,2);
        Kitty.SetNumPW(3,2);
        Kitty.SetSpin(0,0);
        Kitty.SetSpin(1,1);
        Kitty.SetSpin(2,1);
        Kitty.SetSpin(3,1);
        Kitty.SetChannelWeight(0, 3./12.);
        Kitty.SetChannelWeight(1, 1./12.);
        Kitty.SetChannelWeight(2, 3./12.);
        Kitty.SetChannelWeight(3, 5./12.);
    }
    else{
        //important: even with the p-waves switched off,
        //physics wise the spin 1 state still exists and
        //the p-waves are there, just in there asymptotic state (free wave).
        //To include this in the computation,
        //CATS still needs a second channel, even if it is `empty`!
        Kitty.SetNumChannels(2);
        if(dwaves) Kitty.SetNumPW(0,3);
        else Kitty.SetNumPW(0,1);
        Kitty.SetNumPW(1,0);
        Kitty.SetSpin(0,0);
        Kitty.SetSpin(1,1);
        Kitty.SetChannelWeight(0, 1./4.);
        Kitty.SetChannelWeight(1, 3./4.);
    }

    //to set up the strong interaction, one can use the
    //predefined functions available in DLM_Potentials.h
    //The input parameters are by default 9 and are defined as follows:
    //0: potential flag
    //1: a second flag, that can be used if needed
    //2: total isospin
    //3: 2 x isospin of particle 1
    //4: 2 x isospin of particle 2
    //5: total spin
    //6: l quantum number
    //7: j quantum number
    CATSparameters cPars_pp_1S0(CATSparameters::tPotential,8,true);
    cPars_pp_1S0.SetParameter(0,NN_ReidV8);//choose the AV18
    cPars_pp_1S0.SetParameter(1,v18_Coupled3P2);
    cPars_pp_1S0.SetParameter(2,1);
    cPars_pp_1S0.SetParameter(3,1);
    cPars_pp_1S0.SetParameter(4,1);
    cPars_pp_1S0.SetParameter(5,0);
    cPars_pp_1S0.SetParameter(6,0);
    cPars_pp_1S0.SetParameter(7,0);

    //copy all settings from cPars_pp_1S0, and just change quantum numbers s,l,j
    CATSparameters cPars_pp_3P0(cPars_pp_1S0);
    cPars_pp_3P0.SetParameter(5,1);
    cPars_pp_3P0.SetParameter(6,1);
    cPars_pp_3P0.SetParameter(7,0);

    CATSparameters cPars_pp_3P1(cPars_pp_1S0);
    cPars_pp_3P1.SetParameter(5,1);
    cPars_pp_3P1.SetParameter(6,1);
    cPars_pp_3P1.SetParameter(7,1);

    CATSparameters cPars_pp_3P2(cPars_pp_1S0);
    cPars_pp_3P2.SetParameter(5,1);
    cPars_pp_3P2.SetParameter(6,1);
    cPars_pp_3P2.SetParameter(7,2);

    CATSparameters cPars_pp_1D2(cPars_pp_1S0);
    cPars_pp_1D2.SetParameter(5,0);
    cPars_pp_1D2.SetParameter(6,2);
    cPars_pp_1D2.SetParameter(7,2);

    //plug in the strong potential for each channel and partial wave
    //the arguments are:
    //#WhichChannel,#WhichPartialWave,#PotentialFunction,#PotentialParameters
    Kitty.SetShortRangePotential(0,0,fDlmPot,cPars_pp_1S0);
    if(pwaves){
        Kitty.SetShortRangePotential(1,1,fDlmPot,cPars_pp_3P0);
        Kitty.SetShortRangePotential(2,1,fDlmPot,cPars_pp_3P1);
        Kitty.SetShortRangePotential(3,1,fDlmPot,cPars_pp_3P2);
    }
    if(dwaves){
        Kitty.SetShortRangePotential(0,2,fDlmPot,cPars_pp_1D2);
    }
    //if later on you would like to switch some contribution off,
    //this can be done with:
    //Kitty.RemoveShortRangePotential(#WhichChannel,#WhichPartialWave);
}

//setup the basic properties for pLambda
void CATS_pL_Basic(CATS& Kitty, const unsigned& NumMomBins, const double& kMin, const double& kMax){
    Kitty.SetMomBins(NumMomBins,kMin,kMax);
    //should you include in the result any bins,
    //where the Schroedinger solver failed
    Kitty.SetExcludeFailedBins(false);
    Kitty.SetQ1Q2(0);
    Kitty.SetQuantumStatistics(false);
    Kitty.SetRedMass( (Mass_p*Mass_L)/(Mass_p+Mass_L) );
}
//setup the CATS interaction for pLambda (Usmani)
void CATS_pL_Usmani(CATS& Kitty){
    Kitty.SetNumChannels(2);
    Kitty.SetNumPW(0,1);
    Kitty.SetNumPW(1,1);
    Kitty.SetSpin(0,0);
    Kitty.SetSpin(1,1);
    //set up the weights of the spin 0 and 1 channel.
    //the singlet/triplet configuration leads to these weights
    Kitty.SetChannelWeight(0, 1./4.);
    Kitty.SetChannelWeight(1, 3./4.);

    //define the parameters for the 1S0 partial wave
    CATSparameters POT_PARS_1S0(CATSparameters::tPotential,1,true);
    //the only parameter is the spin (here 0)
    POT_PARS_1S0.SetParameter(0,0);
    //define the parameters for the 3S1 partial wave
    CATSparameters POT_PARS_3S1(CATSparameters::tPotential,1,true);
    //the only parameter is the spin (here 1)
    POT_PARS_3S1.SetParameter(0,1);
    //For each channel and each partial wave, the strong potential
    //has to be added using this function, taking as input arguments
    //#WhichChannel,#WhichPW,PotentialFunction,Parameters
    Kitty.SetShortRangePotential(0,0,UsmaniPotentialCats,POT_PARS_1S0);
    Kitty.SetShortRangePotential(1,0,UsmaniPotentialCats,POT_PARS_3S1);
}

void CATS_pXim_Basic(CATS& Kitty,
        const unsigned& NumMomBins, const double& kMin, const double& kMax){
    Kitty.SetMomBins(NumMomBins,kMin,kMax);
    Kitty.SetExcludeFailedBins(false);
    Kitty.SetQ1Q2(-1);
    Kitty.SetQuantumStatistics(false);
    Kitty.SetRedMass( (Mass_p*Mass_Xim)/(Mass_p+Mass_Xim) );
}
void CATS_pXim_Hal(CATS& Kitty){
    //0: S=0 I=0
    //1: S=1 I=0
    //2: S=0 I=1
    //3: S=1 I=1
    Kitty.SetNumChannels(4);
    Kitty.SetNumPW(0,1);
    Kitty.SetNumPW(1,1);
    Kitty.SetNumPW(2,1);
    Kitty.SetNumPW(3,1);
    Kitty.SetSpin(0,0);
    Kitty.SetSpin(1,1);
    Kitty.SetSpin(2,0);
    Kitty.SetSpin(3,1);
    Kitty.SetChannelWeight(0, 1./8.);
    Kitty.SetChannelWeight(1, 3./8.);
    Kitty.SetChannelWeight(2, 1./8.);
    Kitty.SetChannelWeight(3, 3./8.);

    CATSparameters cPars_pXim_HalI0S0(CATSparameters::tPotential,8,true);
    cPars_pXim_HalI0S0.SetParameter(0,pXim_HALQCD1);
    //the euclidean time
    cPars_pXim_HalI0S0.SetParameter(1,12);
    cPars_pXim_HalI0S0.SetParameter(2,0);
    cPars_pXim_HalI0S0.SetParameter(3,-1);
    cPars_pXim_HalI0S0.SetParameter(4,1);
    cPars_pXim_HalI0S0.SetParameter(5,0);
    cPars_pXim_HalI0S0.SetParameter(6,0);
    cPars_pXim_HalI0S0.SetParameter(7,0);

    CATSparameters cPars_pXim_HalI0S1(CATSparameters::tPotential,8,true);
    cPars_pXim_HalI0S1.SetParameter(0,pXim_HALQCD1);
    //the euclidean time
    cPars_pXim_HalI0S1.SetParameter(1,12);
    cPars_pXim_HalI0S1.SetParameter(2,0);
    cPars_pXim_HalI0S1.SetParameter(3,-1);
    cPars_pXim_HalI0S1.SetParameter(4,1);
    cPars_pXim_HalI0S1.SetParameter(5,1);
    cPars_pXim_HalI0S1.SetParameter(6,0);
    cPars_pXim_HalI0S1.SetParameter(7,1);

    CATSparameters cPars_pXim_HalI1S0(CATSparameters::tPotential,8,true);
    cPars_pXim_HalI1S0.SetParameter(0,pXim_HALQCD1);
    //the euclidean time
    cPars_pXim_HalI1S0.SetParameter(1,12);
    cPars_pXim_HalI1S0.SetParameter(2,1);
    cPars_pXim_HalI1S0.SetParameter(3,1);
    cPars_pXim_HalI1S0.SetParameter(4,1);
    cPars_pXim_HalI1S0.SetParameter(5,0);
    cPars_pXim_HalI1S0.SetParameter(6,0);
    cPars_pXim_HalI1S0.SetParameter(7,0);

    CATSparameters cPars_pXim_HalI1S1(CATSparameters::tPotential,8,true);
    cPars_pXim_HalI1S1.SetParameter(0,pXim_HALQCD1);
    //the euclidean time
    cPars_pXim_HalI1S1.SetParameter(1,12);
    cPars_pXim_HalI1S1.SetParameter(2,1);
    cPars_pXim_HalI1S1.SetParameter(3,1);
    cPars_pXim_HalI1S1.SetParameter(4,1);
    cPars_pXim_HalI1S1.SetParameter(5,1);
    cPars_pXim_HalI1S1.SetParameter(6,0);
    cPars_pXim_HalI1S1.SetParameter(7,1);

    Kitty.SetShortRangePotential(0,0,fDlmPot,cPars_pXim_HalI0S0);
    Kitty.SetShortRangePotential(1,0,fDlmPot,cPars_pXim_HalI0S1);
    Kitty.SetShortRangePotential(2,0,fDlmPot,cPars_pXim_HalI1S0);
    Kitty.SetShortRangePotential(3,0,fDlmPot,cPars_pXim_HalI1S1);
}

double pXiPhenomenologicalCk(const double& Momentum, const double* SourcePar, const double* PotPar){
    return 1 + exp(-2.*Momentum*SourcePar[0])/(2.*Momentum*SourcePar[0]);
}

TGraph* Basics_ProtonLambda(){

    const unsigned NumMomBins = 150;
    const double kMin = 0;
    const double kMax = 300;

    CATS Kitty_pL;
    CATS_pL_Basic(Kitty_pL,NumMomBins,kMin,kMax);
    CATS_GaussSource(Kitty_pL,1.5);
    CATS_pL_Usmani(Kitty_pL);

    //these are two optional functions, which can make the computing
    //grid finer, at the expense of increased CPU time to solve
    //the Schroedinger equation. The default value is 5e-6. Here we choose
    //a very fine grid, in order to produce a smooth 2 MeV binned correlation
    Kitty_pL.SetEpsilonConv(5e-8);
    Kitty_pL.SetEpsilonProp(5e-8);
    Kitty_pL.KillTheCat();

    //save the result in a TGraph (ROOT)
    TGraph* grCk = new TGraph();
    grCk->SetName("grCk");
    grCk->Set(NumMomBins);
    //iterate over all momentum bins
    for(unsigned uMom=0; uMom<NumMomBins; uMom++){
        //the GetMomentum() function sets the x-axis
        //the GetCorrFun() function evaluates C(k) (y-axis)
        grCk->SetPoint(uMom,Kitty_pL.GetMomentum(uMom),Kitty_pL.GetCorrFun(uMom));
    }

    return grCk;
}

TGraph* ProtonLambda_SimpleXiFeeddown(){
    const unsigned NumMomBins = 150;
    const double kMin = 0;
    const double kMax = 300;
    const double SourceSize = 1.5;

    //the standard set up of the pLambda CATS object
    CATS Kitty_pL;
    CATS_pL_Basic(Kitty_pL,NumMomBins,kMin,kMax);
    CATS_GaussSource(Kitty_pL,SourceSize);
    CATS_pL_Usmani(Kitty_pL);
    Kitty_pL.SetEpsilonConv(5e-8);
    Kitty_pL.SetEpsilonProp(5e-8);
    Kitty_pL.KillTheCat();

    //transfer the result to a histogram of the type Ck_pL
    DLM_Ck Ck_pL(Kitty_pL.GetNumSourcePars(),0,Kitty_pL);
    //Update (i.e. fill) the histogram
    Ck_pL.Update();

    //to take pXim->pLambda feeddown, we need the pXim
    //correlation function. Here we model it with a simple
    //phenomenological function (exponential), which takes a
    //single input parameter (to mimic the source size)
    DLM_Ck Ck_pXim(1,0,NumMomBins,kMin,kMax,pXiPhenomenologicalCk);
    Ck_pXim.SetSourcePar(0,SourceSize);
    Ck_pXim.Update();

    //get the residual matrix from a file
    TString FileName_Feed_pXim_pL = "../Files/DecayMatrices_Oli.root";
    TString HistoName_Feed_pXim_pL = "hRes_pp_pL";
    //TH2F* hFeedDown_pL_pp = GetSmearMatrix(FileName_Feed_pXim_pL,HistoName_Feed_pXim_pL);


    //DLM_CkDecomposition CkDec_pL("pLambda",3,Ck_pL,HistoName_Feed_pXim_pL);
    //DLM_CkDecomposition CkDec_pXim("pXim",2,Ck_pXim,NULL);





}

//an example for which we consider the pLambda
//correlation function including the feeddown
//from pXi. The rest of the feeddown is flat
TGraph* ProtonLambda_XiFeeddown(){
    const unsigned NumMomBins = 150;
    const double kMin = 0;
    const double kMax = 300;

    CATS Kitty_pL;
    CATS_pL_Basic(Kitty_pL,NumMomBins,kMin,kMax);
    CATS_GaussSource(Kitty_pL,1.5);
    CATS_pL_Usmani(Kitty_pL);

    CATS Kitty_pXim;
    CATS_pXim_Basic(Kitty_pL,NumMomBins,kMin,kMax);
    CATS_GaussSource(Kitty_pXim,1.5);


}

void RepulsiveCore_AppendixExample(){

    const unsigned NumPts = 256;
    const double rMin = 0;
    const double rMax = 3;
    const double Vmin = 0;
    const double Vmax = 11;

    double* Parameters1 = new double [5];
    Parameters1[2] = 10;
    Parameters1[3] = 1;
    Parameters1[4] = 0.2;

    double* Parameters2 = new double [5];
    Parameters2[2] = 10;
    Parameters2[3] = 1;
    Parameters2[4] = 0.02;

    TGraph grRepCore1;
    grRepCore1.SetName("grRepCore1");
    grRepCore1.Set(NumPts);
    grRepCore1.SetLineColor(kBlue+2);
    grRepCore1.SetLineWidth(4);

    TGraph grRepCore2;
    grRepCore2.SetName("grRepCore2");
    grRepCore2.Set(NumPts);
    grRepCore2.SetLineColor(kBlue+2);
    grRepCore2.SetLineWidth(4);
    grRepCore2.SetLineStyle(2);

    TH1F* hAxis = new TH1F("hAxis","hAxis",NumPts,rMin,rMax);
    hAxis->GetYaxis()->SetRangeUser(Vmin,Vmax);
    hAxis->GetXaxis()->SetTitle("r (fm)");
    hAxis->GetYaxis()->SetTitle("V (GeV)");
    SetStyleAxis_2(hAxis);

    for(unsigned uPts=0; uPts<NumPts; uPts++){
        Parameters1[0] = hAxis->GetBinCenter(uPts+1);
        Parameters2[0] = hAxis->GetBinCenter(uPts+1);
        grRepCore1.SetPoint(uPts,hAxis->GetBinCenter(uPts+1),RepulsiveCore(Parameters1));
        grRepCore2.SetPoint(uPts,hAxis->GetBinCenter(uPts+1),RepulsiveCore(Parameters2));
    }

    TCanvas* cRepCore = new TCanvas("cRepCore", "cRepCore", 1);
    cRepCore->cd(0); cRepCore->SetCanvasSize(1920, 1080); cRepCore->SetMargin(0.15,0.05,0.2,0.05);//lrbt
    cRepCore->SetGrid(true);

    hAxis->Draw();
    grRepCore1.Draw("C,same");
    grRepCore2.Draw("C,same");

    //TFile fOutput("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/PhD_Thesis/Pictures/CATS_examples/RepulsiveCoreExample.root","recreate");
    cRepCore->SaveAs("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/PhD_Thesis/Pictures/CATS_examples/RepulsiveCoreExample.pdf");

    delete hAxis;
    delete cRepCore;
    delete [] Parameters1;
    delete [] Parameters2;
}

TGraph* Lednicky_ProtonLambda(){
    const unsigned NumMomBins = 150;
    const double kMin = 0;
    const double kMax = 300;
    const double SourceSize = 1.5;
    //anchored to the Usmani potential
    const double ScattLen_1S0 = 2.88;
    const double EffRan_1S0 = 2.92;
    const double ScattLen_3S1 = 1.66;
    const double EffRan_3S1 = 3.78;

    //initialize the histogram object using the Lednicky model
    //the first two arguments represent the number of source parameters
    //and the number of potential parameters required by the function
    DLM_Ck CkProtonLambda_Lednicky(1,4,NumMomBins,kMin,kMax,Lednicky_Triplet);
    //set the source and interaction parameters
    CkProtonLambda_Lednicky.SetSourcePar(0,SourceSize);
    CkProtonLambda_Lednicky.SetPotPar(0,ScattLen_1S0);
    CkProtonLambda_Lednicky.SetPotPar(1,EffRan_1S0);
    CkProtonLambda_Lednicky.SetPotPar(2,ScattLen_3S1);
    CkProtonLambda_Lednicky.SetPotPar(3,EffRan_3S1);
    //fill the histogram
    CkProtonLambda_Lednicky.Update();

    //save the result in a TGraph (ROOT)
    TGraph* grCk = new TGraph();
    grCk->SetName("grCk_Lednicky_ProtonLambda");
    grCk->Set(NumMomBins);
    //iterate over all momentum bins
    for(unsigned uMom=0; uMom<NumMomBins; uMom++){
        //DLM_Ck has similar function naming as a ROOT TH1F.
        //however the bin counting starts from zero, and the GetBinCenter
        //has an extra parameter corresponding to the axis number.
        //DLM_Ck is 1D, thus this parameter can only be zero
        grCk->SetPoint(uMom,CkProtonLambda_Lednicky.GetBinCenter(0,uMom),
                       CkProtonLambda_Lednicky.GetBinContent(uMom));
    }

    return grCk;
}

void PlotCompareUsmaniLednicky(){
    TGraph* CkUsmani = Basics_ProtonLambda();
    TGraph* CkLednicky = Lednicky_ProtonLambda();

    CkUsmani->SetLineColor(kRed+1);
    CkUsmani->SetLineWidth(4);

    CkLednicky->SetLineColor(kCyan+1);
    CkLednicky->SetLineWidth(4);
    CkLednicky->SetLineStyle(2);

    TH1F* hAxis = new TH1F("hAxis","hAxis",128,0,300);
    SetStyleAxis_1(hAxis);
    hAxis->GetYaxis()->SetRangeUser(0.9, 3.35);
    hAxis->GetXaxis()->SetTitle("k* (MeV)");
    hAxis->GetYaxis()->SetTitle("C(k*)");

    TF1* UnitLine = new TF1("UnitLine","1",-20,320);
    UnitLine->SetLineColor(kBlack);
    UnitLine->SetLineWidth(2.0);
    UnitLine->SetLineStyle(1);

    TPaveText* PTinfo = new TPaveText(0.25,0.85,0.5,0.95, "blNDC");//lbrt
    PTinfo->SetName("PTinfo");
    PTinfo->SetBorderSize(1);
    PTinfo->SetTextSize(0.05);
    PTinfo->SetFillColor(kWhite);
    PTinfo->SetTextFont(22);
    PTinfo->AddText("Gaussian source");

    TLegend *legend = new TLegend(0.5,0.8,0.95,0.95);
    legend->SetTextFont(42);
    legend->SetTextSize(0.05);
    legend->AddEntry(CkUsmani, "p#minus#Lambda (Usmani)");
    legend->AddEntry(CkLednicky, "p#minus#Lambda (Lednicky)");

    TCanvas* cpL = new TCanvas("cpL", "cpL", 10, 10, 1280, 720);
    //cpL->SetCanvasSize(1920, 1080);
    cpL->SetMargin(0.15,0.05,0.2,0.05);//lrbt
    cpL->SetGrid(1,1);
    //cpL->cd(0);

    hAxis->Draw();
    UnitLine->Draw("same");
    CkUsmani->Draw("C,same");
    CkLednicky->Draw("C,same");
    legend->Draw("same");
    PTinfo->Draw("same");

    //INLET
    TH1F* hAxis_Inlet = (TH1F*)hAxis->Clone("hAxis_Inlet");
    hAxis_Inlet->GetXaxis()->SetTitleSize(0);
    hAxis_Inlet->GetXaxis()->SetLabelSize(0);
    hAxis_Inlet->GetXaxis()->SetTitleOffset(0);
    hAxis_Inlet->GetXaxis()->SetLabelOffset(0);
    hAxis_Inlet->GetXaxis()->SetRangeUser(100, 300);
    hAxis_Inlet->GetXaxis()->SetNdivisions(505);
    hAxis_Inlet->GetXaxis()->SetTitle("");

    hAxis_Inlet->GetYaxis()->SetTitleSize(hAxis->GetYaxis()->GetTitleSize()*1.75);
    hAxis_Inlet->GetYaxis()->SetLabelSize(hAxis->GetYaxis()->GetLabelSize()*1.75);
    hAxis_Inlet->GetYaxis()->SetTitleOffset(hAxis->GetYaxis()->GetTitleOffset()*0.5);
    hAxis_Inlet->GetYaxis()->SetRangeUser(0.965, 1.015);

    const double fXMinInlet=0.190;
    const double fYMinInlet=0.26;
    const double fXMaxInlet=0.989;
    const double fYMaxInlet=0.70;
    TPad *inset_pad = new TPad("inset_pad", "inset_pad", fXMinInlet, fYMinInlet, fXMaxInlet, fYMaxInlet);
    inset_pad->SetTopMargin(0.01);
    inset_pad->SetRightMargin(0.05);
    inset_pad->SetBottomMargin(0.28);
    inset_pad->SetLeftMargin(0.28);
    inset_pad->SetFillStyle(4000);
    inset_pad->Draw();
    inset_pad->cd();
    inset_pad->SetGrid();
    hAxis_Inlet->Draw();
    UnitLine->Draw("same");
    CkUsmani->Draw("C,same");
    CkLednicky->Draw("C,same");

    TString OutFileName = "UsmaniLednicky_Gauss";

    TFile* OutputFile = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/PhD_Thesis/Pictures/CATS_examples/"+OutFileName+".root","recreate");
    CkUsmani->Write();
    CkLednicky->Write();
    cpL->Write();

    cpL->SaveAs("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/PhD_Thesis/Pictures/CATS_examples/"+OutFileName+".pdf");

    delete CkUsmani;
    delete CkLednicky;
    delete UnitLine;
    delete hAxis;
    delete legend;
    delete PTinfo;
    delete cpL;
    delete OutputFile;
}


DLM_Histo<float>* GetSmearMatrix_pp(){
    DLM_Histo<float>* MATRIX = NULL;
    TString filename = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/ResolutionMatrices/Sample6_MeV_compact.root";
    TString histoname = "hSigmaMeV_Proton_Proton";
    TFile* FileROOT = new TFile(filename, "read");
    TH2F* histo = (TH2F*)FileROOT->Get(histoname);
    if(!histo){printf("\033[1;31mERROR:\033[0m The histo '%s' in file '%s' does not exist\n",histoname.Data(),filename.Data());return NULL;}
    MATRIX = Convert_TH2F_DlmHisto(histo);
    delete FileROOT;
    return MATRIX;
}
DLM_Histo<float>* GetDecayMatrix_pL_pp(){
    DLM_Histo<float>* MATRIX = NULL;
    TString filename = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/DecayMatrices/run2_decay_matrices_old.root";
    TString histoname = "hRes_pp_pL";
    TFile* FileROOT = new TFile(filename, "read");
    TH2F* histo = (TH2F*)FileROOT->Get(histoname);
    if(!histo){printf("\033[1;31mERROR:\033[0m The histo '%s' in file '%s' does not exist\n",histoname.Data(),filename.Data());return NULL;}
    MATRIX = Convert_TH2F_DlmHisto(histo);
    delete FileROOT;
    return MATRIX;
}
DLM_Histo<float>* GetDecayMatrix_pXim_pL(){
    DLM_Histo<float>* MATRIX = NULL;
    TString filename = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/DecayMatrices/run2_decay_matrices_old.root";
    TString histoname = "hRes_pL_pXim";
    TFile* FileROOT = new TFile(filename, "read");
    TH2F* histo = (TH2F*)FileROOT->Get(histoname);
    if(!histo){printf("\033[1;31mERROR:\033[0m The histo '%s' in file '%s' does not exist\n",histoname.Data(),filename.Data());return NULL;}
    MATRIX = Convert_TH2F_DlmHisto(histo);
    delete FileROOT;
    return MATRIX;
}
DLM_Histo<float>* GetDecayMatrix_pSigma0_pL(){
    DLM_Histo<float>* MATRIX = NULL;
    TString filename = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/DecayMatrices/run2_decay_matrices_old.root";
    TString histoname = "hRes_pL_pSigma0";
    TFile* FileROOT = new TFile(filename, "read");
    TH2F* histo = (TH2F*)FileROOT->Get(histoname);
    if(!histo){printf("\033[1;31mERROR:\033[0m The histo '%s' in file '%s' does not exist\n",histoname.Data(),filename.Data());return NULL;}
    MATRIX = Convert_TH2F_DlmHisto(histo);
    delete FileROOT;
    return MATRIX;
}

///an example how compute the pp theoretical correlation
void Ck_pp_Decomposition(){

    //the same source profile and size will be assumed for all
    //contributions. This does NOT need to be the case
    const double SourceSize = 1.5;
    const double NumMomBins = 40;
    const double kMin = 0;
    const double kMax = 160;

    //set up pp
    CATS Kitty_pp;
    CATS_GaussSource(Kitty_pp,SourceSize);
    CATS_pp_Basic(Kitty_pp,NumMomBins,kMin,kMax);
    //setup the pp interaction with s,p,d waves
    CATS_pp_AV18(Kitty_pp,true,true);
    Kitty_pp.KillTheCat();
    //create the C(k) histogram
    DLM_Ck dlmCk_pp(1,0,Kitty_pp);

    //set up pLambda
    CATS Kitty_pL;
    CATS_GaussSource(Kitty_pL,SourceSize);
    CATS_pL_Basic(Kitty_pL,NumMomBins,kMin,kMax);
    CATS_pL_Usmani(Kitty_pL);
    Kitty_pL.KillTheCat();
    DLM_Ck dlmCk_pL(1,0,Kitty_pL);

    //set up pXi
    CATS Kitty_pXim;
    CATS_GaussSource(Kitty_pXim,SourceSize);
    CATS_pXim_Basic(Kitty_pXim,NumMomBins,kMin,kMax);
    CATS_pXim_Hal(Kitty_pXim);
    Kitty_pXim.KillTheCat();
    DLM_Ck dlmCk_pXim(1,0,Kitty_pXim);

    //set up pSigma0
    DLM_Ck dlmCk_pSigma0(1,0,SourceSize,kMin,kMax,Lednicky_gauss_Sigma0);
    dlmCk_pSigma0.SetSourcePar(0,SourceSize);

    //one could also get the histograms in a TH2F format,
    //however than one should use the DLM_CkDecomposition
    //object instead of DLM_CkDecomp

    //momentum smearing matrix.
    //it is applied only to the final correlation (pp)
    DLM_Histo<float>* hResolution_pp = GetSmearMatrix_pp();
    //the decay matrix for the feeddown of pLambda->pp
    DLM_Histo<float>* hFeedDown_pL_pp = GetDecayMatrix_pL_pp();
    //the decay matrix for the feeddown of pXim->pLambda
    DLM_Histo<float>* hFeedDown_pXim_pL = GetDecayMatrix_pXim_pL();
    //the decay matrix for the feeddown of pSigma0->pL
    DLM_Histo<float>* hFeedDown_pSigma0_pL = GetDecayMatrix_pSigma0_pL();

    //the object is initialized by: a name of your choice,
    //number of contributions to C(k) APART from the primary,
    //DLM_Ck object, resolution matrix (NULL to avoid mom. smear)
    DLM_CkDecomp CkDec_pp("pp",3,dlmCk_pp,hResolution_pp);
    DLM_CkDecomp CkDec_pL("pLambda",3,dlmCk_pL,NULL);
    DLM_CkDecomp CkDec_pXim("pXim",1,dlmCk_pXim,NULL);
    DLM_CkDecomp CkDec_pSigma0("pSigma0",1,dlmCk_pSigma0,NULL);

    //example lambda pars (more or less realistic values)
    const double lambda_pL_pp = 0.15;
    const double lambda_pp_flat = 0.08;
    const double lambda_pp_misid = 0.02;
    //the primary lambda is assumed to be:
    //1. - lambda_pL_pp - lambda_pp_flat - lambda_pp_misid
    //to add the different contributions to the correlation, we need to specify:
    //#WhichContribution,#Type(cFeedDown or cFake),
    //OPTIONALLY: DLM_CkDecomp of the feed-down contribution, transformation matrix
    //if the last two entries are not given,
    //it is assumed that the feed-down correlation is flat
    CkDec_pp.AddContribution(0,lambda_pL_pp,DLM_CkDecomp::cFeedDown,
                             &CkDec_pL,hFeedDown_pL_pp);
    CkDec_pp.AddContribution(1,lambda_pp_flat,DLM_CkDecomp::cFeedDown);
    CkDec_pp.AddContribution(2,lambda_pp_misid,DLM_CkDecomp::cFake);

    double lambda_pXim_pL = 0.10;
    double lambda_pSigma0_pL = 0.17;
    double lambda_pL_flat = 0.23;
    CkDec_pL.AddContribution(0,lambda_pXim_pL,DLM_CkDecomp::cFeedDown,
                             &CkDec_pXim,hFeedDown_pXim_pL);
    CkDec_pL.AddContribution(1,lambda_pSigma0_pL,DLM_CkDecomp::cFeedDown,
                             &CkDec_pSigma0,hFeedDown_pSigma0_pL);
    CkDec_pL.AddContribution(2,lambda_pL_flat,DLM_CkDecomp::cFeedDown);

    double lambda_pXim_flat = 0.37;
    CkDec_pXim.AddContribution(0,lambda_pXim_flat,DLM_CkDecomp::cFeedDown);

//based on the paper, I think it should be 0.176
    double lambda_pSigma0_flat = 0.12;
    CkDec_pSigma0.AddContribution(0,lambda_pSigma0_flat,DLM_CkDecomp::cFeedDown);

    //computes everything, including the decomposition of the correlation signal C(k)-1,
    //i.e. how much is each contribution to the total correlation
    CkDec_pp.Update(true,true);
    CkDec_pL.Update(true,true);
    CkDec_pXim.Update(true,true);
    CkDec_pSigma0.Update(true,true);

    TString OutFileName = "Decomposition_pp";
    TFile* OutputFile = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/PhD_Thesis/Pictures/CATS_examples/"+OutFileName+".root","recreate");

    TGraph* gCk_pp = new TGraph();
    gCk_pp->SetName("gCk_pp");

    TGraph* gCk_pp_pp = new TGraph();
    gCk_pp_pp->SetName("gCk_pp_pp");

    TGraph* gCk_pL_pp = new TGraph();
    gCk_pL_pp->SetName("gCk_pL_pp");

    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        //the momentum
        double Momentum = dlmCk_pp.GetBinCenter(0,uBin);
        //evaluates the total pp correlation,
        //including momentum smearing and feed-downs
        double Ck_pp = CkDec_pp.EvalCk(Momentum);
        //the contribution of the primary correlation (smeared)
        double Ck_pp_pp = CkDec_pp.EvalSignalSmearedMain(Momentum);
        //the contribution of the pL->pp correlation (smeared)
        double Ck_pL_pp = CkDec_pp.EvalSignalSmearedChild(0,Momentum);
        gCk_pp->SetPoint(uBin,Momentum,Ck_pp);
        gCk_pp_pp->SetPoint(uBin,Momentum,Ck_pp_pp+1);
        gCk_pL_pp->SetPoint(uBin,Momentum,Ck_pL_pp+1);
    }

    //PLOTING MACROS
    //...

    gCk_pp->SetLineColor(kBlue+1);
    gCk_pp->SetLineWidth(4);

    gCk_pp_pp->SetLineColor(kBlue+1);
    gCk_pp_pp->SetLineWidth(4);
    gCk_pp_pp->SetLineStyle(2);

    gCk_pL_pp->SetLineColor(kRed+1);
    gCk_pL_pp->SetLineWidth(4);
    gCk_pL_pp->SetLineStyle(2);

    OutputFile->cd();

    TH1F* hAxis = new TH1F("hAxis","hAxis",NumMomBins,kMin,kMax);
    SetStyleAxis_1(hAxis);
    hAxis->GetYaxis()->SetRangeUser(0.9, 3.05);
    hAxis->GetXaxis()->SetTitle("k* (MeV)");
    hAxis->GetYaxis()->SetTitle("C(k*)");

    TF1* UnitLine = new TF1("UnitLine","1",-20,320);
    UnitLine->SetLineColor(kBlack);
    UnitLine->SetLineWidth(2.0);
    UnitLine->SetLineStyle(1);

    TPaveText* PTinfo = new TPaveText(0.25,0.85,0.5,0.95, "blNDC");//lbrt
    PTinfo->SetName("PTinfo");
    PTinfo->SetBorderSize(1);
    PTinfo->SetTextSize(0.05);
    PTinfo->SetFillColor(kWhite);
    PTinfo->SetTextFont(22);
    PTinfo->AddText("Gaussian source");

    TLegend *legend = new TLegend(0.5,0.75,0.95,0.95);
    legend->SetTextFont(42);
    legend->SetTextSize(0.05);
    legend->AddEntry(gCk_pp, "p#minusp (total)");
    legend->AddEntry(gCk_pp_pp, "p#minusp#rightarrowp#minusp");
    legend->AddEntry(gCk_pL_pp, "p#minus#Lambda#rightarrowp#minusp");

    TCanvas* canDecomp = new TCanvas("canDecomp", "canDecomp", 1);
    canDecomp->SetCanvasSize(1280, 720);
    canDecomp->SetMargin(0.15,0.05,0.2,0.05);//lrbt
    canDecomp->SetGrid(1,1);
    //canDecomp->cd(0);

    hAxis->Draw();
    UnitLine->Draw("same");
    gCk_pp->Draw("C,same");
    gCk_pp_pp->Draw("C,same");
    gCk_pL_pp->Draw("C,same");
    legend->Draw("same");
    PTinfo->Draw("same");

    //INLET
    TH1F* hAxis_Inlet = (TH1F*)hAxis->Clone("hAxis_Inlet");
    //hAxis_Inlet->GetXaxis()->SetTitleSize(0);
    //hAxis_Inlet->GetXaxis()->SetLabelSize(0);
    //hAxis_Inlet->GetXaxis()->SetTitleOffset(0);
    //hAxis_Inlet->GetXaxis()->SetLabelOffset(0);
    //hAxis_Inlet->GetXaxis()->SetRangeUser(100, 300);
    hAxis_Inlet->GetXaxis()->SetTitleSize(hAxis->GetXaxis()->GetTitleSize()*1.75);
    hAxis_Inlet->GetXaxis()->SetLabelSize(hAxis->GetXaxis()->GetLabelSize()*1.75);
    hAxis_Inlet->GetXaxis()->SetTitleOffset(hAxis->GetXaxis()->GetTitleOffset()*0.75);
    hAxis_Inlet->GetXaxis()->SetNdivisions(505);
    //hAxis_Inlet->GetXaxis()->SetTitle("");

    hAxis_Inlet->GetYaxis()->SetTitleSize(hAxis->GetYaxis()->GetTitleSize()*1.75);
    hAxis_Inlet->GetYaxis()->SetLabelSize(hAxis->GetYaxis()->GetLabelSize()*1.75);
    hAxis_Inlet->GetYaxis()->SetTitleOffset(hAxis->GetYaxis()->GetTitleOffset()*0.5);
    hAxis_Inlet->GetYaxis()->SetRangeUser(0.95, 1.05);

    const double fXMinInlet=0.215;
    const double fYMinInlet=0.31;
    const double fXMaxInlet=0.988;
    const double fYMaxInlet=0.75;
    TPad *inset_pad = new TPad("inset_pad", "inset_pad", fXMinInlet, fYMinInlet, fXMaxInlet, fYMaxInlet);
    inset_pad->SetTopMargin(0.01);
    inset_pad->SetRightMargin(0.05);
    inset_pad->SetBottomMargin(0.28);
    inset_pad->SetLeftMargin(0.28);
    inset_pad->SetFillStyle(4000);
    inset_pad->Draw();
    inset_pad->cd();
    inset_pad->SetGrid();
    hAxis_Inlet->Draw();
    UnitLine->Draw("same");
    gCk_pp->Draw("C,same");
    gCk_pp_pp->Draw("C,same");
    gCk_pL_pp->Draw("C,same");

    gCk_pp->Write();
    gCk_pp_pp->Write();
    gCk_pL_pp->Write();
    canDecomp->Write();

    canDecomp->SaveAs("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/PhD_Thesis/Pictures/CATS_examples/"+OutFileName+".pdf");


    //do not forget to free the memory
    delete hResolution_pp;
    delete hFeedDown_pL_pp;
    delete hFeedDown_pSigma0_pL;
    delete hFeedDown_pXim_pL;

    delete gCk_pp;
    delete gCk_pp_pp;
    delete gCk_pL_pp;

    printf("The end\n");
}
/*
void ModifiedGauss(double* x, double* par){
    par[0]*exp(-pow((x-[par[1])/par[2],2.)*(1.-exp(-pow(,2.)))
}

void ModifiedGauss(double* x, double* par){
    const int Lambda = 4;
    const double Mean = 2.0;
    const double BinProb = Mean/double(Lambda);

    const unsigned NumIter = 1000000;

    TH1F* hPDF = new TH1F("hPDF","hPDF",256,0,16);
    DLM_Random rangen(11);

    for(unsigned uIter=0; uIter<NumIter; uIter++){
        hPDF->Fill(rangen.Binomial(Lambda,BinProb));
    }

    TFile fOutput("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/PhD_Thesis/Tests/ModifiedPoisson/fOutput.root","recreate");

    hPDF->Write();

    delete hPDF;
}
*/
void pL_Feed(){

    TString OutputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/PhD_Thesis/Pictures/Matrices/";
    TString FileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/DecayMatrices/run2_decay_matrices_old.root";
    TString HistoName = "hRes_pp_pL";

    //original histo
    TFile InputFile(FileName,"read");
    TH2F* hOriginal = (TH2F*)InputFile.Get(HistoName);
    if(!hOriginal) {printf("!hOriginal\n"); return;}
    const unsigned NumTotBins = hOriginal->GetXaxis()->GetNbins();
    const unsigned FirstBin = hOriginal->GetXaxis()->FindBin(double(0));
    const unsigned LastBin = hOriginal->GetXaxis()->FindBin(double(499));
    const double kMin = hOriginal->GetXaxis()->GetBinLowEdge(FirstBin);
    const double kMax = hOriginal->GetXaxis()->GetBinUpEdge(LastBin);
    const unsigned NumMomBins = LastBin-FirstBin+1;
    const unsigned Rescale = 5;


    const double SourceSize = 1.5;
    const double NumCatsBins = 80;
    const double kMinCats = 0;
    const double kMaxCats = 200;
    const double kMaxCatsPlot = 160;
    //set up pp
    CATS Kitty_pp;
    CATS_GaussSource(Kitty_pp,SourceSize);
    CATS_pp_Basic(Kitty_pp,NumCatsBins,kMinCats,kMaxCats);
    CATS_pp_AV18(Kitty_pp,false,false);
    Kitty_pp.KillTheCat();
    DLM_Ck dlmCk_pp(1,0,Kitty_pp);

    //set up pLambda
    CATS Kitty_pL;
    CATS_GaussSource(Kitty_pL,SourceSize);
    CATS_pL_Basic(Kitty_pL,NumCatsBins,kMinCats,kMaxCats);
    CATS_pL_Usmani(Kitty_pL);
    Kitty_pL.SetEpsilonConv(1e-8);
    Kitty_pL.SetEpsilonProp(1e-8);
    Kitty_pL.KillTheCat();
    DLM_Ck dlmCk_pL(1,0,Kitty_pL);

    DLM_CkDecomposition CkDec_pp("pp",1,dlmCk_pp,NULL);
    DLM_CkDecomposition CkDec_pL("pLambda",0,dlmCk_pL,NULL);
    CkDec_pp.AddContribution(0,1.-1e-6,DLM_CkDecomp::cFeedDown,&CkDec_pL,hOriginal);
    CkDec_pp.Update(true,true);
    CkDec_pL.Update(true,true);
    TGraph gCk_pL;
    gCk_pL.SetName("gCk_pL");
    gCk_pL.SetLineColor(kRed+2);
    gCk_pL.SetLineWidth(4);
    gCk_pL.SetLineStyle(1);
    TGraph gCk_pL_pp;
    gCk_pL_pp.SetName("gCk_pL_pp");
    gCk_pL_pp.SetLineColor(kRed+1);
    gCk_pL_pp.SetLineWidth(4);
    gCk_pL_pp.SetLineStyle(2);
    TGraph gCkLamda_pL_pp;
    gCkLamda_pL_pp.SetName("gCkLamda_pL_pp");
    gCkLamda_pL_pp.SetLineColor(kAzure+1);
    gCkLamda_pL_pp.SetLineWidth(4);
    gCkLamda_pL_pp.SetLineStyle(5);
    for(unsigned uBin=0; uBin<NumCatsBins; uBin++){
        double MOMENTUM = Kitty_pL.GetMomentum(uBin);
        if(MOMENTUM>kMaxCatsPlot) break;
        gCk_pL.SetPoint(uBin,MOMENTUM,CkDec_pL.EvalSignalMain(MOMENTUM)+1);
        gCk_pL_pp.SetPoint(uBin,MOMENTUM,CkDec_pp.EvalSignalChild(0,MOMENTUM)+1);
        gCkLamda_pL_pp.SetPoint(uBin,MOMENTUM,0.08*CkDec_pp.EvalSignalChild(0,MOMENTUM)+1);
    }

    //renormalized
    TH2F* hRenorm = new TH2F("hRenorm","hRenorm",NumMomBins,kMin,kMax,NumMomBins,kMin,kMax);

    //profiles (at 25, 50, 100 MeV)
    const unsigned NumProf = 5;
    double ProfVal[NumProf];
    ProfVal[0] = 15;
    ProfVal[1] = 30;
    ProfVal[2] = 60;
    ProfVal[3] = 90;
    ProfVal[4] = 120;
    TGraph** gProf = new TGraph* [NumProf];
    TH1F** hProf = new TH1F* [NumProf];

    TH1D* hProj;
    for(unsigned uBinX=0; uBinX<NumMomBins; uBinX++){
        hProj = hOriginal->ProjectionY("hProj",uBinX+1,uBinX+1);
        hProj->Sumw2();
        hProj->Scale(1./hProj->Integral(),"width");
        for(unsigned uBinY=0; uBinY<NumMomBins; uBinY++){
            hRenorm->SetBinContent(uBinX+1,uBinY+1,hProj->GetBinContent(uBinY+1));
            hRenorm->SetBinError(uBinX+1,uBinY+1,hProj->GetBinError(uBinY+1));
        }
        delete hProj;
    }

    for(unsigned uProf=0; uProf<NumProf; uProf++){
        gProf[uProf] = new TGraph();
        gProf[uProf]->SetName(TString::Format("gProf_%.0fMeV",ProfVal[uProf]));
        hProf[uProf] = new TH1F(TString::Format("hProf_%.0fMeV",ProfVal[uProf]),TString::Format("hProf_%.0fMeV",ProfVal[uProf]),
                                NumMomBins,kMin,kMax);
        unsigned WhichXbin = hRenorm->GetXaxis()->FindBin(ProfVal[uProf]);
        hProj = hRenorm->ProjectionY("hProj",WhichXbin,WhichXbin);
        for(unsigned uBinY=0; uBinY<NumMomBins; uBinY++){
            hProf[uProf]->SetBinContent(uBinY+1,hRenorm->GetBinContent(WhichXbin,uBinY+1));
            hProf[uProf]->SetBinError(uBinY+1,hRenorm->GetBinError(WhichXbin,uBinY+1));

        }
        //hProf[uProf]->Sumw2();
        hProj->Rebin(Rescale);
        hProj->Scale(1./double(Rescale));
        for(unsigned uBinY=0; uBinY<hProj->GetNbinsX(); uBinY++){
            double MOM = hProj->GetBinCenter(uBinY+1);
            double VAL = hProj->GetBinContent(uBinY+1);
            gProf[uProf]->SetPoint(uBinY,MOM,VAL);
        }
    }

    TGraph* gProfOrig = new TGraph [NumProf+1];
    TGraph* gProfOrig2 = new TGraph [NumProf+1];

    gProfOrig[0].SetLineColor(kViolet+2);
    gProfOrig[1].SetLineColor(kBlue+2);
    gProfOrig[2].SetLineColor(kGreen+2);
    gProfOrig[3].SetLineColor(kYellow+2);
    gProfOrig[4].SetLineColor(kRed+2);
    gProfOrig[5].SetLineColor(kBlack);

    gProfOrig[0].SetMarkerColor(kViolet+2);
    gProfOrig[1].SetMarkerColor(kBlue+2);
    gProfOrig[2].SetMarkerColor(kGreen+2);
    gProfOrig[3].SetMarkerColor(kYellow+2);
    gProfOrig[4].SetMarkerColor(kRed+2);
    gProfOrig[5].SetMarkerColor(kBlack);

    for(unsigned uProf=0; uProf<=NumProf; uProf++){
        gProfOrig[uProf].SetLineWidth(0);
        gProfOrig[uProf].SetMarkerStyle(uProf<NumProf?29:30);
        gProfOrig[uProf].SetMarkerSize(5);

        gProfOrig2[uProf].SetLineWidth(0);
        gProfOrig2[uProf].SetMarkerStyle(uProf<NumProf?29:30);
        gProfOrig2[uProf].SetMarkerSize(3);
        gProfOrig2[uProf].SetMarkerColor(kWhite);

        double xVal = uProf<NumProf?ProfVal[uProf]:0;
        //double yVal = hProf[uProf]->GetBinContent(hProf[uProf]->FindBin(xVal));
        double yVal = 0.0105;
        //printf("xMean : yVal = %f : %f\n",xMean,yVal);
        gProfOrig[uProf].SetPoint(0,xVal,yVal);
        gProfOrig2[uProf].SetPoint(0,xVal,yVal);
    }


    SetStyleAxis2D_2(hOriginal);
    hOriginal->GetXaxis()->SetTitle("k' (p#minus#Lambda) (MeV)");
    hOriginal->GetXaxis()->SetRangeUser(kMin,kMax);
    hOriginal->GetYaxis()->SetTitle("k* (p#minusp) (MeV)");
    hOriginal->GetYaxis()->SetRangeUser(kMin,kMax);
    hOriginal->GetZaxis()->SetTitle("N");
    TCanvas* cOriginal = new TCanvas("cOriginal", "cOriginal", 1);
    cOriginal->cd(0); cOriginal->SetCanvasSize(1280, 720); cOriginal->SetMargin(0.12,0.2,0.17,0.05);//lrbt
    cOriginal->SetGrid(true);
    hOriginal->Draw("COLZ");

    SetStyleAxis2D_2(hRenorm);
    hRenorm->GetXaxis()->SetTitle("k' (p#minus#Lambda) (MeV)");
    hRenorm->GetXaxis()->SetRangeUser(kMin,kMax);
    hRenorm->GetYaxis()->SetTitle("k* (p#minusp) (MeV)");
    hRenorm->GetYaxis()->SetRangeUser(kMin,kMax);
    hRenorm->GetZaxis()->SetTitle("dN/dk* (1/MeV)");
    hRenorm->GetZaxis()->SetRangeUser(0.0,0.02);
    TCanvas* cRenorm = new TCanvas("cRenorm", "cRenorm", 1);
    cRenorm->cd(0); cRenorm->SetCanvasSize(1280, 720); cRenorm->SetMargin(0.12,0.2,0.17,0.05);//lrbt
    cRenorm->SetGrid(true);
    hRenorm->Draw("COLZ");

    TLegend* lLegend_Prof = new TLegend(0.60,0.60,0.95,0.95);//lbrt
    lLegend_Prof->SetName(TString::Format("lLegend_Prof"));
    lLegend_Prof->SetTextSize(0.05);
    for(unsigned uProf=0; uProf<NumProf; uProf++){
        lLegend_Prof->AddEntry(hProf[uProf],TString::Format("{k'=%.0f MeV} #rightarrow {k*}",ProfVal[uProf]));
    }
    lLegend_Prof->AddEntry(&gProfOrig[NumProf],"k' (p#minus#Lambda)");

    TH1F* hAxisProf = new TH1F("hAxisProf","hAxisProf",128,0,192);
    SetStyleAxis_2(hAxisProf);
    hAxisProf->GetYaxis()->SetRangeUser(0.0, 0.05);
    hAxisProf->GetXaxis()->SetTitle("k* (p#minusp) (MeV)");
    hAxisProf->GetYaxis()->SetTitle("dN/dk* (1/MeV)");
    TCanvas* cProf = new TCanvas(TString::Format("cProf"), TString::Format("cProf"), 1);
    cProf->cd(0); cProf->SetCanvasSize(1280, 720); cProf->SetMargin(0.15,0.05,0.2,0.05);//lrbt
    cProf->SetGrid(true);
    hAxisProf->Draw("axis");
    hProf[0]->SetLineColor(kViolet+2);
    hProf[1]->SetLineColor(kBlue+2);
    hProf[2]->SetLineColor(kGreen+2);
    hProf[3]->SetLineColor(kYellow+2);
    hProf[4]->SetLineColor(kRed+2);
    for(unsigned uProf=0; uProf<NumProf; uProf++){
        hProf[uProf]->SetLineWidth(2.5);
        hProf[uProf]->Draw("same");
        gProfOrig[uProf].Draw("same,P");
        gProfOrig2[uProf].Draw("same,P");
    }
    lLegend_Prof->Draw("same");


    TH1F* hAxisCk = new TH1F("hAxisCk","hAxisCk",NumCatsBins,kMinCats,kMaxCatsPlot);
    SetStyleAxis_1(hAxisCk);
    hAxisCk->GetYaxis()->SetRangeUser(0.9, 3.05);
    hAxisCk->GetXaxis()->SetTitle("k* (MeV)");
    hAxisCk->GetYaxis()->SetTitle("C(k*)");

    TF1* UnitLine = new TF1("UnitLine","1",-20,320);
    UnitLine->SetLineColor(kBlack);
    UnitLine->SetLineWidth(2.0);
    UnitLine->SetLineStyle(1);

    TPaveText* PTinfoCk = new TPaveText(0.25,0.85,0.6,0.95, "blNDC");//lbrt
    PTinfoCk->SetName("PTinfoCk");
    PTinfoCk->SetBorderSize(1);
    PTinfoCk->SetTextSize(0.05);
    PTinfoCk->SetFillColor(kWhite);
    PTinfoCk->SetTextFont(22);
    PTinfoCk->AddText("Gauss source (1.5 fm)");

    TLegend *legendCk = new TLegend(0.6,0.75,0.95,0.95);
    legendCk->SetTextFont(42);
    legendCk->SetTextSize(0.05);
    legendCk->AddEntry(&gCk_pL, "p#minus#Lambda");
    legendCk->AddEntry(&gCk_pL_pp, "p#minus#Lambda#rightarrowp#minusp");
    legendCk->AddEntry(&gCkLamda_pL_pp, "p#minus#Lambda#rightarrowp#minusp (15%)");

    TCanvas* canDecomp = new TCanvas("canDecomp", "canDecomp", 1);
    canDecomp->SetCanvasSize(1280, 720);
    canDecomp->SetMargin(0.15,0.05,0.2,0.05);//lrbt
    canDecomp->SetGrid(1,1);
    //canDecomp->cd(0);

    hAxisCk->Draw();
    UnitLine->Draw("same");
    gCk_pL.Draw("C,same");
    gCk_pL_pp.Draw("C,same");
    gCkLamda_pL_pp.Draw("C,same");
    legendCk->Draw("same");
    PTinfoCk->Draw("same");

    //INLET
    TH1F* hAxisCk_Inlet = (TH1F*)hAxisCk->Clone("hAxisCk_Inlet");
    hAxisCk_Inlet->GetXaxis()->SetTitleSize(0);
    hAxisCk_Inlet->GetXaxis()->SetLabelSize(0);
    hAxisCk_Inlet->GetXaxis()->SetTitleOffset(0);
    hAxisCk_Inlet->GetXaxis()->SetRangeUser(70, kMaxCatsPlot);
    hAxisCk_Inlet->GetXaxis()->SetNdivisions(505);

    hAxisCk_Inlet->GetYaxis()->SetTitleSize(hAxisCk->GetYaxis()->GetTitleSize()*1.75);
    hAxisCk_Inlet->GetYaxis()->SetLabelSize(hAxisCk->GetYaxis()->GetLabelSize()*1.75);
    hAxisCk_Inlet->GetYaxis()->SetTitleOffset(hAxisCk->GetYaxis()->GetTitleOffset()*0.6);
    hAxisCk_Inlet->GetYaxis()->SetRangeUser(0.96, 1.165);

    const double fXMinInlet=0.311;
    const double fYMinInlet=0.31;
    const double fXMaxInlet=0.984;
    const double fYMaxInlet=0.76;
    TPad *inset_pad = new TPad("inset_pad", "inset_pad", fXMinInlet, fYMinInlet, fXMaxInlet, fYMaxInlet);
    inset_pad->SetTopMargin(0.01);
    inset_pad->SetRightMargin(0.05);
    inset_pad->SetBottomMargin(0.28);
    inset_pad->SetLeftMargin(0.28);
    inset_pad->SetFillStyle(4000);
    inset_pad->Draw();
    inset_pad->cd();
    inset_pad->SetGrid();
    hAxisCk_Inlet->Draw();
    UnitLine->Draw("same");
    gCk_pL.Draw("C,same");
    gCk_pL_pp.Draw("C,same");
    gCkLamda_pL_pp.Draw("C,same");

    TFile OutputFile(OutputFolder+"fPLfeed.root","recreate");

    hOriginal->Write();
    hRenorm->Write();
    for(unsigned uProf=0; uProf<NumProf; uProf++){
        hProf[uProf]->Write();
        gProf[uProf]->Write();
    }
    gCk_pL.Write();
    gCk_pL_pp.Write();
    gCkLamda_pL_pp.Write();
    cOriginal->Write();
    cRenorm->Write();
    cProf->Write();
    canDecomp->Write();

    cOriginal->SaveAs(OutputFolder+"c_pL_pp_Matrix.png");
    cRenorm->SaveAs(OutputFolder+"cRenorm_pL_pp_Matrix.png");
    cProf->SaveAs(OutputFolder+"cProf.pdf");
    canDecomp->SaveAs(OutputFolder+"CkFeed_pL_pp.pdf");


    for(unsigned uProf=0; uProf<NumProf; uProf++) {delete hProf[uProf]; delete gProf[uProf];}
    delete [] gProf;
    delete [] hProf;
    delete [] gProfOrig;
    delete [] gProfOrig2;
    delete lLegend_Prof;
    delete hAxisCk;
    delete UnitLine;
    delete PTinfoCk;
    delete legendCk;
    delete inset_pad;
    delete canDecomp;
    delete cProf;
    delete hRenorm;
    delete hAxisProf;
}

double FitSmearMatrix(double* x, double* par){
    double& Amplitude = par[0];
    double& Mean = par[1];
    double& Sigma = par[2];
    double cpi = TMath::Pi();
    //double Norm = (2.*exp(Mean/Sigma)*sqrt(1./Sigma/Sigma))/
    //(sqrt(cpi)*(exp(Mean/Sigma)+
    //            exp(Mean/Sigma)*TMath::Erf(Mean*sqrt(1./Sigma/Sigma))+
    //            exp(0.25)*TMath::Erfc(0.5*sqrt(1./Sigma/Sigma)*(-2.*Mean+Sigma))));
    //return Amplitude*Norm*exp(-pow((*x-Mean)/Sigma,2.))*(1.-exp(-pow(*x/Sigma,2.)));
    return Amplitude*TMath::Gaus(*x,Mean,Sigma,true);
}

void pp_MomReso(){
    TString OutputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/PhD_Thesis/Pictures/Matrices/";
    TString FileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/ResolutionMatrices/Sample6_MeV_compact.root";
    TString HistoName = "hSigmaMeV_Proton_Proton";

    //original histo
    TFile InputFile(FileName,"read");
    TH2F* hOriginal = (TH2F*)InputFile.Get(HistoName);
    if(!hOriginal) {printf("!hOriginal\n"); return;}
    const unsigned NumTotBins = hOriginal->GetXaxis()->GetNbins();
    const unsigned FirstBin = hOriginal->GetXaxis()->FindBin(double(0));
    const unsigned LastBin = hOriginal->GetXaxis()->FindBin(double(499));
    const double kMin = hOriginal->GetXaxis()->GetBinLowEdge(FirstBin);
    const double kMax = hOriginal->GetXaxis()->GetBinUpEdge(LastBin);
    const unsigned NumMomBins = LastBin-FirstBin+1;
    const double Rescale = 5;

    const double SourceSize = 1.5;
    const double NumCatsBins = 80;
    const double kMinCats = 0;
    const double kMaxCats = 200;
    const double kMaxCatsPlot = 160;
    //set up pp
    CATS Kitty_pp;
    CATS_GaussSource(Kitty_pp,SourceSize);
    CATS_pp_Basic(Kitty_pp,NumCatsBins,kMinCats,kMaxCats);
    CATS_pp_AV18(Kitty_pp,true,true);
    Kitty_pp.KillTheCat();
    DLM_Ck dlmCk_pp(1,0,Kitty_pp);

    DLM_CkDecomposition CkDec_pp("pp",0,dlmCk_pp,hOriginal);
    CkDec_pp.Update(true,true);


    TGraph gCk_pp;
    gCk_pp.SetName("gCk_pp");
    gCk_pp.SetLineColor(kBlue+2);
    gCk_pp.SetLineWidth(4);
    gCk_pp.SetLineStyle(1);
    TGraph gCk_pp_Smeared;
    gCk_pp_Smeared.SetName("gCk_pp_Smeared");
    gCk_pp_Smeared.SetLineColor(kBlue+1);
    gCk_pp_Smeared.SetLineWidth(4);
    gCk_pp_Smeared.SetLineStyle(2);
    for(unsigned uBin=0; uBin<NumCatsBins; uBin++){
        double MOMENTUM = Kitty_pp.GetMomentum(uBin);
        if(MOMENTUM>kMaxCatsPlot) break;
        gCk_pp.SetPoint(uBin,MOMENTUM,CkDec_pp.EvalSignal(MOMENTUM)+1);
        gCk_pp_Smeared.SetPoint(uBin,MOMENTUM,CkDec_pp.EvalSignalSmeared(MOMENTUM)+1);
    }

    //renormalized
    TH2F* hRenorm = new TH2F("hRenorm","hRenorm",NumMomBins,kMin,kMax,NumMomBins,kMin,kMax);

    const unsigned NumProf = 5;
    double ProfVal[NumProf];
    ProfVal[0] = 5;
    ProfVal[1] = 20;
    ProfVal[2] = 40;
    ProfVal[3] = 60;
    ProfVal[4] = 80;
    TGraph** gProf = new TGraph* [NumProf];
    TH1F** hProf = new TH1F* [NumProf];
    TF1** fProf = new TF1* [NumProf];

    TH1D* hProj;
    for(unsigned uBinX=0; uBinX<NumMomBins; uBinX++){
        hProj = hOriginal->ProjectionY("hProj",uBinX+1,uBinX+1);
        hProj->Sumw2();
        hProj->Scale(1./hProj->Integral(),"width");
        for(unsigned uBinY=0; uBinY<NumMomBins; uBinY++){
            hRenorm->SetBinContent(uBinX+1,uBinY+1,hProj->GetBinContent(uBinY+1));
            hRenorm->SetBinError(uBinX+1,uBinY+1,hProj->GetBinError(uBinY+1));
        }
        delete hProj;
    }

    for(unsigned uProf=0; uProf<NumProf; uProf++){
        gProf[uProf] = new TGraph();
        gProf[uProf]->SetName(TString::Format("gProf_%.0fMeV",ProfVal[uProf]));
        hProf[uProf] = new TH1F(TString::Format("hProf_%.0fMeV",ProfVal[uProf]),TString::Format("hProf_%.0fMeV",ProfVal[uProf]),
                                NumMomBins,kMin,kMax);
        fProf[uProf] = new TF1(TString::Format("fProf_%.0fMeV",ProfVal[uProf]),FitSmearMatrix,kMinCats,kMaxCatsPlot,3);
        fProf[uProf]->FixParameter(0,1);
        //fProf[uProf]->SetParLimits(0,0.1,10);
        fProf[uProf]->FixParameter(1,ProfVal[uProf]);
        //fProf[uProf]->SetParLimits(1,ProfVal[uProf]*0.9,ProfVal[uProf]*1.1);
        fProf[uProf]->SetParameter(2,10);
        fProf[uProf]->SetParLimits(2,1,20);
        unsigned WhichXbin = hRenorm->GetXaxis()->FindBin(ProfVal[uProf]);
        hProj = hRenorm->ProjectionY("hProj",WhichXbin,WhichXbin);
        for(unsigned uBinY=0; uBinY<NumMomBins; uBinY++){
            hProf[uProf]->SetBinContent(uBinY+1,hRenorm->GetBinContent(WhichXbin,uBinY+1));
            hProf[uProf]->SetBinError(uBinY+1,hRenorm->GetBinError(WhichXbin,uBinY+1));

        }
        //hProf[uProf]->Fit(fProf[uProf],"S, N, R, M");
        //hProf[uProf]->Sumw2();
        hProj->Rebin(Rescale);
        hProj->Scale(1./double(Rescale));
        for(unsigned uBinY=0; uBinY<hProj->GetNbinsX(); uBinY++){
            double MOM = hProj->GetBinCenter(uBinY+1);
            double VAL = hProj->GetBinContent(uBinY+1);
            gProf[uProf]->SetPoint(uBinY,MOM,VAL);
        }
    }

    TGraph* gProfOrig = new TGraph [NumProf+1];
    TGraph* gProfOrig2 = new TGraph [NumProf+1];

    gProfOrig[0].SetLineColor(kViolet+2);
    gProfOrig[1].SetLineColor(kBlue+2);
    gProfOrig[2].SetLineColor(kGreen+2);
    gProfOrig[3].SetLineColor(kYellow+2);
    gProfOrig[4].SetLineColor(kRed+2);
    gProfOrig[5].SetLineColor(kBlack);

    gProfOrig[0].SetMarkerColor(kViolet+2);
    gProfOrig[1].SetMarkerColor(kBlue+2);
    gProfOrig[2].SetMarkerColor(kGreen+2);
    gProfOrig[3].SetMarkerColor(kYellow+2);
    gProfOrig[4].SetMarkerColor(kRed+2);
    gProfOrig[5].SetMarkerColor(kBlack);

    for(unsigned uProf=0; uProf<=NumProf; uProf++){
        gProfOrig[uProf].SetLineWidth(0);
        gProfOrig[uProf].SetMarkerStyle(uProf<NumProf?29:30);
        gProfOrig[uProf].SetMarkerSize(5);

        gProfOrig2[uProf].SetLineWidth(0);
        gProfOrig2[uProf].SetMarkerStyle(uProf<NumProf?29:30);
        gProfOrig2[uProf].SetMarkerSize(3);
        gProfOrig2[uProf].SetMarkerColor(kWhite);

        double xVal = uProf<NumProf?ProfVal[uProf]:0;
        //double yVal = hProf[uProf]->GetBinContent(hProf[uProf]->FindBin(xVal));
        double yVal = 0.0105;
        //printf("xMean : yVal = %f : %f\n",xMean,yVal);
        gProfOrig[uProf].SetPoint(0,xVal,yVal);
        gProfOrig2[uProf].SetPoint(0,xVal,yVal);
    }


    SetStyleAxis2D_2(hOriginal);
    hOriginal->GetXaxis()->SetTitle("k' (true) (MeV)");
    hOriginal->GetXaxis()->SetRangeUser(kMin,kMax);
    hOriginal->GetYaxis()->SetTitle("k* (reconstructed) (MeV)");
    hOriginal->GetYaxis()->SetRangeUser(kMin,kMax);
    hOriginal->GetZaxis()->SetTitle("N");
    hRenorm->GetZaxis()->SetRangeUser(0.9,30000);
    TCanvas* cOriginal = new TCanvas("cOriginal", "cOriginal", 1);
    cOriginal->cd(0); cOriginal->SetCanvasSize(1280, 720); cOriginal->SetMargin(0.12,0.2,0.17,0.05);//lrbt
    cOriginal->SetGrid(true);
    cOriginal->SetLogz(true);
    hOriginal->Draw("COLZ");

    SetStyleAxis2D_2(hRenorm);
    hRenorm->GetXaxis()->SetTitle("k' (true) (MeV)");
    hRenorm->GetXaxis()->SetRangeUser(kMin,kMax);
    hRenorm->GetYaxis()->SetTitle("k* (reconstructed) (MeV)");
    hRenorm->GetYaxis()->SetRangeUser(kMin,kMax);
    hRenorm->GetZaxis()->SetTitle("dN/dk* (1/MeV)");
    hRenorm->GetZaxis()->SetRangeUser(0.000001,0.1);
    TCanvas* cRenorm = new TCanvas("cRenorm", "cRenorm", 1);
    cRenorm->cd(0); cRenorm->SetCanvasSize(1280, 720); cRenorm->SetMargin(0.12,0.2,0.17,0.05);//lrbt
    cRenorm->SetGrid(true);
    cRenorm->SetLogz(true);
    hRenorm->Draw("COLZ");

    TLegend* lLegend_Prof = new TLegend(0.69,0.65,0.95,0.95);//lbrt
    lLegend_Prof->SetName(TString::Format("lLegend_Prof"));
    lLegend_Prof->SetTextSize(0.04);
    for(unsigned uProf=0; uProf<NumProf; uProf++){
        lLegend_Prof->AddEntry(hProf[uProf],TString::Format("{k'=%.0f MeV} #rightarrow {k*}",ProfVal[uProf]));
    }
    //lLegend_Prof->AddEntry(&gProfOrig[NumProf],"k' (p#minus#Lambda)");

    TH1F* hAxisProf = new TH1F("hAxisProf","hAxisProf",128,0,128);
    SetStyleAxis_2(hAxisProf);
    hAxisProf->GetYaxis()->SetRangeUser(0.0, 0.11);
    hAxisProf->GetXaxis()->SetTitle("k* (p#minusp) (MeV)");
    hAxisProf->GetYaxis()->SetTitle("dN/dk* (1/MeV)");
    TCanvas* cProf = new TCanvas(TString::Format("cProf"), TString::Format("cProf"), 1);
    cProf->cd(0); cProf->SetCanvasSize(1280, 720); cProf->SetMargin(0.15,0.05,0.2,0.05);//lrbt
    cProf->SetGrid(true);
    hAxisProf->Draw("axis");
    hProf[0]->SetLineColor(kViolet+2);
    hProf[1]->SetLineColor(kBlue+2);
    hProf[2]->SetLineColor(kGreen+2);
    hProf[3]->SetLineColor(kYellow+2);
    hProf[4]->SetLineColor(kRed+2);
    for(unsigned uProf=0; uProf<NumProf; uProf++){
        fProf[uProf]->SetLineColor(hProf[uProf]->GetLineColor());
        fProf[uProf]->SetLineWidth(3.5);
        hProf[uProf]->SetLineWidth(2.0);
        hProf[uProf]->Draw("same");
        //fProf[uProf]->Draw("same");
        //gProfOrig[uProf].Draw("same,P");
        //gProfOrig2[uProf].Draw("same,P");
    }
    lLegend_Prof->Draw("same");


    TH1F* hAxisCk = new TH1F("hAxisCk","hAxisCk",NumCatsBins,kMinCats,kMaxCatsPlot);
    SetStyleAxis_1(hAxisCk);
    hAxisCk->GetYaxis()->SetRangeUser(0.0, 3.8);
    hAxisCk->GetXaxis()->SetTitle("k* (MeV)");
    hAxisCk->GetYaxis()->SetTitle("C(k*)");

    TF1* UnitLine = new TF1("UnitLine","1",-20,320);
    UnitLine->SetLineColor(kBlack);
    UnitLine->SetLineWidth(2.0);
    UnitLine->SetLineStyle(1);

    TPaveText* PTinfoCk = new TPaveText(0.30,0.85,0.6,0.95, "blNDC");//lbrt
    PTinfoCk->SetName("PTinfoCk");
    PTinfoCk->SetBorderSize(1);
    PTinfoCk->SetTextSize(0.05);
    PTinfoCk->SetFillColor(kWhite);
    PTinfoCk->SetTextFont(22);
    PTinfoCk->AddText("Gauss source (1.5 fm)");

    TLegend *legendCk = new TLegend(0.6,0.75,0.95,0.95);
    legendCk->SetTextFont(42);
    legendCk->SetTextSize(0.05);
    legendCk->AddEntry(&gCk_pp, "p#minusp (theory)");
    legendCk->AddEntry(&gCk_pp_Smeared, "p#minusp (smeared)");

    TCanvas* canDecomp = new TCanvas("canDecomp", "canDecomp", 1);
    canDecomp->SetCanvasSize(1280, 720);
    canDecomp->SetMargin(0.15,0.05,0.2,0.05);//lrbt
    canDecomp->SetGrid(1,1);
    //canDecomp->cd(0);

    hAxisCk->Draw();
    UnitLine->Draw("same");
    gCk_pp.Draw("C,same");
    gCk_pp_Smeared.Draw("C,same");
    legendCk->Draw("same");
    PTinfoCk->Draw("same");

    //INLET
    TH1F* hAxisCk_Inlet = (TH1F*)hAxisCk->Clone("hAxisCk_Inlet");
    hAxisCk_Inlet->GetXaxis()->SetTitleSize(0);
    hAxisCk_Inlet->GetXaxis()->SetLabelSize(0);
    hAxisCk_Inlet->GetXaxis()->SetTitleOffset(0);
    hAxisCk_Inlet->GetXaxis()->SetRangeUser(70, kMaxCatsPlot);
    hAxisCk_Inlet->GetXaxis()->SetNdivisions(505);

    hAxisCk_Inlet->GetYaxis()->SetTitleSize(hAxisCk->GetYaxis()->GetTitleSize()*1.75);
    hAxisCk_Inlet->GetYaxis()->SetLabelSize(hAxisCk->GetYaxis()->GetLabelSize()*1.75);
    hAxisCk_Inlet->GetYaxis()->SetTitleOffset(hAxisCk->GetYaxis()->GetTitleOffset()*0.6);
    hAxisCk_Inlet->GetYaxis()->SetRangeUser(0.93, 1.02);

    const double fXMinInlet=0.311;
    const double fYMinInlet=0.31;
    const double fXMaxInlet=0.984;
    const double fYMaxInlet=0.76;
    TPad *inset_pad = new TPad("inset_pad", "inset_pad", fXMinInlet, fYMinInlet, fXMaxInlet, fYMaxInlet);
    inset_pad->SetTopMargin(0.01);
    inset_pad->SetRightMargin(0.05);
    inset_pad->SetBottomMargin(0.28);
    inset_pad->SetLeftMargin(0.28);
    inset_pad->SetFillStyle(4000);
    inset_pad->Draw();
    inset_pad->cd();
    inset_pad->SetGrid();
    hAxisCk_Inlet->Draw();
    UnitLine->Draw("same");
    gCk_pp.Draw("C,same");
    gCk_pp_Smeared.Draw("C,same");

    TFile OutputFile(OutputFolder+"fppSigma.root","recreate");

    hOriginal->Write();
    hRenorm->Write();
    for(unsigned uProf=0; uProf<NumProf; uProf++){
        hProf[uProf]->Write();
        fProf[uProf]->Write();
        gProf[uProf]->Write();
    }
    gCk_pp.Write();
    gCk_pp_Smeared.Write();
    cOriginal->Write();
    cRenorm->Write();
    cProf->Write();
    canDecomp->Write();

    cOriginal->SaveAs(OutputFolder+"c_pp_SigmaMatrix.png");
    cRenorm->SaveAs(OutputFolder+"cRenorm_pp_SigmaMatrix.png");
    cProf->SaveAs(OutputFolder+"cProfSigma.pdf");
    canDecomp->SaveAs(OutputFolder+"CkSigma_pp.pdf");

    for(unsigned uProf=0; uProf<NumProf; uProf++) {delete fProf[uProf]; delete hProf[uProf]; delete gProf[uProf];}
    delete [] gProf;
    delete [] hProf;
    delete [] fProf;
    delete [] gProfOrig;
    delete [] gProfOrig2;
    delete lLegend_Prof;
    delete hAxisCk;
    delete UnitLine;
    delete PTinfoCk;
    delete legendCk;
    delete inset_pad;
    delete canDecomp;
    delete cProf;
    delete hRenorm;
    delete hAxisProf;
}

int THESIS_PLOTS(int narg, char** ARGS){
    //ToyPotentials_PS_WF_CF();
    //ComparePionPion(0);
    //ComparePionPion(1);
    //Plot_ProtonLambda(Basics_ProtonLambda());
    //PlotCompareUsmaniLednicky();
    //RepulsiveCore_AppendixExample();
    //Ck_pp_Decomposition();
    //ModifiedPoisson();
    //pL_Feed();
    pp_MomReso();
}
