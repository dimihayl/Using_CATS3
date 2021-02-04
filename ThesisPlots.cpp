#include "ThesisPlots.h"

#include <stdlib.h>

#include "FemtoBoyzScripts.h"
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
#include "DLM_Fitters.h"
#include "CommonAnaFunctions.h"
#include "EnvVars.h"

#include "TStyle.h"
#include "TString.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TPaveText.h"
#include "TROOT.h"

DLM_CkDecomposition* CKDEC_FIT_PP;

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

    const TString OutputFolder = TString::Format("%s/Plots/PhDThesis/CATS/",GetCernBoxDimi()).Data();

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


    TLegend* lLegend_Ck = new TLegend(0.75,0.55,0.90,0.90);//lbrt
    lLegend_Ck->SetName(TString::Format("lLegend_Ck"));
    lLegend_Ck->SetTextSize(0.08);
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
        //PT_RAD[uSor]->SetTextSize(0.055);
        if(uSor==0) PT_RAD[uSor]->SetTextSize(0.065);
        else if(uSor==1) PT_RAD[uSor]->SetTextSize(0.08);
        else PT_RAD[uSor]->SetTextSize(0.065);
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
    //TString OutputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/PhD_Thesis/Pictures/Matrices/";
    //TString FileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/ResolutionMatrices/Sample6_MeV_compact.root";
    //TString HistoName = "hSigmaMeV_Proton_Proton";

    //TString OutputFolder = "/home/dmihaylov/CernBox/CatsFiles/MomentumSmear/NicePlots/";
    //TString FileName = "/home/dmihaylov/CernBox/CatsFiles/MomentumSmear/ALICE_pp_13TeV.root";
    //TString HistoName = "hSigmaMeV_Proton_Lambda";

    TString OutputFolder = "/home/dmihaylov/CernBox/CatsFiles/MomentumSmear/NicePlots/";
    TString FileName = "/home/dmihaylov/CernBox/CatsFiles/MomentumSmear/ALICE_pp_13TeV_MEpL.root";
    TString HistoName = "h_RESO_pL_MeV";

    //original histo
    TFile InputFile(FileName,"read");
    TH2F* hOriginal = (TH2F*)InputFile.Get(HistoName);
    if(!hOriginal) {printf("!hOriginal\n"); return;}
    //hOriginal->Rebin2D(10,10);
    const unsigned NumTotBins = hOriginal->GetXaxis()->GetNbins();
    const unsigned FirstBin = hOriginal->GetXaxis()->FindBin(double(0));
    //const unsigned LastBin = hOriginal->GetXaxis()->FindBin(double(499));
    const unsigned LastBin = hOriginal->GetXaxis()->FindBin(double(249));
    const double kMin = hOriginal->GetXaxis()->GetBinLowEdge(FirstBin);
    const double kMax = hOriginal->GetXaxis()->GetBinUpEdge(LastBin);
    const unsigned NumMomBins = LastBin-FirstBin+1;
    //const double Rescale = 5;
    const double Rescale = 1;


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
    //hAxisProf->GetXaxis()->SetTitle("k* (p#minusp) (MeV)");
    hAxisProf->GetXaxis()->SetTitle("k* (p#minus#Lambda) (MeV)");
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


void QuantumBaseline(){

  printf("Executing QuantumBaseline()\n");
  printf("  CernBox folder: %s\n",GetCernBoxDimi());

  TString OutputFolder = TString::Format("%s/Plots/PhDThesis/QuantumBaseline/",GetCernBoxDimi());
  TString OutputFileName = OutputFolder+"QuantumBaseline.root";
  const unsigned NumSpecies = 10;
  TString* SpeciesName = new TString[NumSpecies];
  TString* SpeciesLegend = new TString[NumSpecies];
  SpeciesName[0] = "pp"; SpeciesLegend[0] = "p-p";
  SpeciesName[1] = "pL"; SpeciesLegend[1] = "p-#Lambda";
  SpeciesName[2] = "pXi"; SpeciesLegend[2] = "p-#Xi";
  SpeciesName[3] = "pOmega"; SpeciesLegend[3] = "p-#Omega";
  SpeciesName[4] = "LL"; SpeciesLegend[4] = "#Lambda-#Lambda";
  SpeciesName[5] = "pAp"; SpeciesLegend[5] = "p-#bar{p}";
  SpeciesName[6] = "pAL"; SpeciesLegend[6] = "p-#bar{#Lambda}";
  SpeciesName[7] = "pAXi"; SpeciesLegend[7] = "p-#bar{#Xi}";
  SpeciesName[8] = "pAOmega"; SpeciesLegend[8] = "p-#bar{#Omega}";
  SpeciesName[9] = "LAL"; SpeciesLegend[9] = "#Lambda-#bar{#Lambda}";

  TString* DataFolder = new TString[NumSpecies];
  DataFolder[0] = TString::Format("%s/CatsFiles/ExpData/pp13TeV_HM_Baseline/MyResults_Vale/Fast_BBar/Data/",GetCernBoxDimi());
  //DataFolder[0] = TString::Format("%s/ALICE_pp_13TeV/Sample9/";
  DataFolder[1] = TString::Format("%s/CatsFiles/ExpData/pp13TeV_HM_Baseline/MyResults_Vale/Fast_BBar/Data/",GetCernBoxDimi());
  DataFolder[2] = TString::Format("%s/CatsFiles/ExpData/pp13TeV_HM_Baseline/MyResults_Vale/Fast_BBar/Data/",GetCernBoxDimi());
  DataFolder[3] = TString::Format("%s/CatsFiles/ExpData/pp13TeV_HM_Baseline/pOmega_1/",GetCernBoxDimi());
  DataFolder[4] = TString::Format("%s/CatsFiles/ExpData/pp13TeV_HM_Baseline/MyResults_Vale/Fast_BBar/Data/",GetCernBoxDimi());
  DataFolder[5] = TString::Format("%s/CatsFiles/ExpData/pp13TeV_HM_Baseline/MyResults_Vale/Fast_BBar/Data/",GetCernBoxDimi());
  DataFolder[6] = TString::Format("%s/CatsFiles/ExpData/pp13TeV_HM_Baseline/MyResults_Vale/Fast_BBar/Data/",GetCernBoxDimi());
  DataFolder[7] = TString::Format("%s/CatsFiles/ExpData/pp13TeV_HM_Baseline/MyResults_Vale/Fast_BBar/Data/",GetCernBoxDimi());
  DataFolder[8] = TString::Format("%s/CatsFiles/ExpData/pp13TeV_HM_Baseline/pOmega_1/",GetCernBoxDimi());
  DataFolder[9] = TString::Format("%s/CatsFiles/ExpData/pp13TeV_HM_Baseline/MyResults_Vale/Fast_BBar/Data/",GetCernBoxDimi());

  TString* McFolder = new TString[NumSpecies];
  McFolder[0] = TString::Format("%s/CatsFiles/ExpData/pp13TeV_HM_Baseline/MyResults_Vale/Fast_BBar/Trains_MCAOD/Norm018028/",GetCernBoxDimi());
  McFolder[1] = TString::Format("%s/CatsFiles/ExpData/pp13TeV_HM_Baseline/MyResults_Vale/Fast_BBar/Trains_MCAOD/Norm018028/",GetCernBoxDimi());
  McFolder[2] = TString::Format("%s/CatsFiles/ExpData/pp13TeV_HM_Baseline/MyResults_Vale/Fast_BBar/Trains_MCAOD/Norm018028/",GetCernBoxDimi());
  McFolder[3] = TString::Format("%s/CatsFiles/ExpData/pp13TeV_HM_Baseline/pOmega_1/",GetCernBoxDimi());
  McFolder[4] = TString::Format("%s/CatsFiles/ExpData/pp13TeV_HM_Baseline/MyResults_Vale/Fast_BBar/Trains_MCAOD/Norm018028/",GetCernBoxDimi());
  McFolder[5] = TString::Format("%s/CatsFiles/ExpData/pp13TeV_HM_Baseline/MyResults_Vale/Fast_BBar/Trains_MCAOD/Norm018028/",GetCernBoxDimi());
  McFolder[6] = TString::Format("%s/CatsFiles/ExpData/pp13TeV_HM_Baseline/MyResults_Vale/Fast_BBar/Trains_MCAOD/Norm018028/",GetCernBoxDimi());
  McFolder[7] = TString::Format("%s/CatsFiles/ExpData/pp13TeV_HM_Baseline/MyResults_Vale/Fast_BBar/Trains_MCAOD/Norm018028/",GetCernBoxDimi());
  McFolder[8] = TString::Format("%s/CatsFiles/ExpData/pp13TeV_HM_Baseline/pOmega_1/",GetCernBoxDimi());
  McFolder[9] = TString::Format("%s/CatsFiles/ExpData/pp13TeV_HM_Baseline/MyResults_Vale/Fast_BBar/Trains_MCAOD/Norm018028/",GetCernBoxDimi());

  TString* DataFile = new TString[NumSpecies];
  DataFile[0] = "CFOutput_pp_8.root";
  //DataFile[0] = "CFOutput_pp.root";
  DataFile[1] = "CFOutput_pL_8.root";
  DataFile[2] = "CFOutput_pXi_8.root";
  DataFile[3] = "outmTsyst.root";
  DataFile[4] = "CFOutput_LL_8.root";
  DataFile[5] = "CFOutput_pAp_8.root";
  DataFile[6] = "CFOutput_pAL_8.root";
  DataFile[7] = "CFOutput_pAXi_8.root";
  DataFile[8] = "outmTpantiomega.root";
  DataFile[9] = "CFOutput_LAL_8.root";

  TString* McFile = new TString[NumSpecies];
  McFile[0] = "CFOutput_pp_8.root";
  McFile[1] = "CFOutput_pL_8.root";
  McFile[2] = "CFOutput_pXi_8.root";
  McFile[3] = "outmTsyst.root";
  McFile[4] = "CFOutput_LL_8.root";
  McFile[5] = "CFOutput_pAp_8.root";
  McFile[6] = "CFOutput_pAL_8.root";
  McFile[7] = "CFOutput_pAXi_8.root";
  McFile[8] = "outmTpantiomega.root";
  McFile[9] = "CFOutput_LAL_8.root";

  TString* DataHisto = new TString[NumSpecies];
  DataHisto[0] = "hCk_ReweightedMeV_0";
  DataHisto[1] = "hCk_ReweightedMeV_0";
  DataHisto[2] = "hCk_ReweightedMeV_0";
  DataHisto[3] = "CF";
  DataHisto[4] = "hCk_ReweightedMeV_0";
  DataHisto[5] = "hCk_ReweightedMeV_0";
  DataHisto[6] = "hCk_ReweightedMeV_0";
  DataHisto[7] = "hCk_ReweightedMeV_0";
  DataHisto[8] = "CF";
  DataHisto[9] = "hCk_ReweightedMeV_0";


  TString* McHisto = new TString[NumSpecies];
  McHisto[0] = "hCk_ReweightedMeV_0";
  McHisto[1] = "hCk_ReweightedMeV_0";
  McHisto[2] = "hCk_ReweightedMeV_0";
  McHisto[3] = "CF";
  McHisto[4] = "hCk_ReweightedMeV_0";
  McHisto[5] = "hCk_ReweightedMeV_0";
  McHisto[6] = "hCk_ReweightedMeV_0";
  McHisto[7] = "hCk_ReweightedMeV_0";
  McHisto[8] = "CF";
  McHisto[9] = "hCk_ReweightedMeV_0";

  int* RescaleFactor = new int[NumSpecies];
  RescaleFactor[0] = 1;//pp
  RescaleFactor[1] = 2;//pL
  RescaleFactor[2] = 2;//pXi
  RescaleFactor[3] = 4;//pOmega
  RescaleFactor[4] = 2;//LL
  RescaleFactor[5] = 1;
  RescaleFactor[6] = 1;
  RescaleFactor[7] = 2;
  RescaleFactor[8] = 4;
  RescaleFactor[9] = 2;

  int* RescaleSeMeAxis = new int[NumSpecies];
  RescaleSeMeAxis[0] = 1000;
  RescaleSeMeAxis[1] = 1000;
  RescaleSeMeAxis[2] = 1000;
  RescaleSeMeAxis[3] = 1;
  RescaleSeMeAxis[4] = 1000;
  RescaleSeMeAxis[5] = 1000;
  RescaleSeMeAxis[6] = 1000;
  RescaleSeMeAxis[7] = 1000;
  RescaleSeMeAxis[8] = 1;
  RescaleSeMeAxis[9] = 1000;

  double* NormFrom = new double[NumSpecies];
  double* NormTo = new double[NumSpecies];
  NormFrom[0] = 200; NormTo[0] = 400;
  NormFrom[1] = 200; NormTo[1] = 400;
  NormFrom[2] = 200; NormTo[2] = 400;
  NormFrom[3] = 200; NormTo[3] = 400;
  NormFrom[4] = 200; NormTo[4] = 400;
  NormFrom[5] = 1000; NormTo[5] = 2000;
  NormFrom[6] = 1000; NormTo[6] = 2000;
  NormFrom[7] = 1000; NormTo[7] = 2000;
  NormFrom[8] = 1000; NormTo[8] = 2000;
  NormFrom[9] = 1000; NormTo[9] = 2000;

  const double fitmin = 600;
  const double fitmax = 1600;
  double* FitFrom = new double[NumSpecies];
  double* FitTo = new double[NumSpecies];
  FitFrom[0] = fitmin; FitTo[0] = fitmax;
  FitFrom[1] = fitmin; FitTo[1] = fitmax;
  FitFrom[2] = fitmin; FitTo[2] = fitmax;
  FitFrom[3] = 0; FitTo[3] = 0;
  FitFrom[4] = fitmin; FitTo[4] = fitmax;
  FitFrom[5] = 0; FitTo[5] = 0;
  FitFrom[6] = 0; FitTo[6] = 0;
  FitFrom[7] = 0; FitTo[7] = 0;
  FitFrom[8] = 0; FitTo[8] = 0;
  FitFrom[9] = 0; FitTo[9] = 0;

  double Purity_p[5];
  double Purity_L[5];
  double Purity_Xim[5];
  double Frac_p[5];
  double Frac_L[5];
  double Frac_Xim[5];
  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.GetPurities_p("pp13TeV_HM_Dec19",0,Purity_p);
  AnalysisObject.GetPurities_L("pp13TeV_HM_Dec19",0,Purity_L);
  AnalysisObject.GetPurities_Xim("pp13TeV_HM_Dec19",0,Purity_Xim);

  AnalysisObject.GetFractions_p("pp13TeV_HM_Dec19",0,Frac_p);
  AnalysisObject.GetFractions_L("pp13TeV_HM_Dec19",0,Frac_L);
  AnalysisObject.GetFractions_Xim("pp13TeV_HM_Dec19",0,Frac_Xim);

  TH1F** hCkPlot = new TH1F* [NumSpecies];
  TF1** fQuantumBaseline = new TF1* [NumSpecies];
  //scan the r0 to make it a common fit parameter
  TF1** fClumsyFit = new TF1* [NumSpecies];
  TF1** fBestClumsyFit = new TF1* [NumSpecies];

  for(unsigned uSpec=0; uSpec<NumSpecies; uSpec++){
    printf("Opening file %s%s\n",DataFolder[uSpec].Data(),DataFile[uSpec].Data());
    TFile* datafile = new TFile(DataFolder[uSpec]+DataFile[uSpec],"read");
    TH1F* hCk_tmp = (TH1F*)datafile->Get(DataHisto[uSpec]);
    hCk_tmp->Rebin(RescaleFactor[uSpec]);
    hCk_tmp->Scale(1./double(RescaleFactor[uSpec]));
    TF1* fNorm = new TF1("fNorm", "[0]",NormFrom[uSpec],NormTo[uSpec]);
    hCk_tmp->Fit(fNorm,"Q, S, N, R, M");
    hCk_tmp->Scale(1./fNorm->GetParameter(0));
    //OutputFile->cd();
    gROOT->cd();
    hCkPlot[uSpec] = (TH1F*)hCk_tmp->Clone(TString::Format("hCkPlot_%s",SpeciesName[uSpec].Data()));
    //hCkPlot[uSpec] = NULL;
    //datafile->cd();
    delete fNorm;
    delete datafile;
  }

  gROOT->cd();
  for(unsigned uSpec=0; uSpec<NumSpecies; uSpec++){
    fQuantumBaseline[uSpec] = new TF1(TString::Format("fQuantumBaseline_%s",SpeciesName[uSpec].Data()),
                                "[0]*(1.+[1]*x)*(1.-[2]*exp(-x*x*[3]*[3]/197.327/197.327))",FitFrom[uSpec],FitTo[uSpec]);
    fQuantumBaseline[uSpec]->SetParameter(0,1.);
    fQuantumBaseline[uSpec]->SetParameter(1,-1e-5);
    fQuantumBaseline[uSpec]->SetParameter(2,0.2);
    fQuantumBaseline[uSpec]->SetParameter(3,0.3);
    fQuantumBaseline[uSpec]->SetParLimits(3,0.25,0.35);
    if(FitFrom[uSpec]==FitTo[uSpec]) continue;
    hCkPlot[uSpec]->Fit(fQuantumBaseline[uSpec],"Q, S, N, R, M");
    printf("For the %s system:\n",SpeciesName[uSpec].Data());
    printf("    N = %.3f +/- %.3f\n",fQuantumBaseline[uSpec]->GetParameter(0),fQuantumBaseline[uSpec]->GetParError(0));
    printf("   p1 = %.3e +/- %.3e 1/MeV\n",fQuantumBaseline[uSpec]->GetParameter(1),fQuantumBaseline[uSpec]->GetParError(1));
    double Strength = fQuantumBaseline[uSpec]->GetParameter(2);
    double StrengthErr = fQuantumBaseline[uSpec]->GetParError(2);
    double Purity;
    double PrimFrac;
    if(uSpec==0){Purity=Purity_p[0]*Purity_p[0];PrimFrac=Frac_p[0]*Frac_p[0];}
    if(uSpec==1){Purity=Purity_p[0]*Purity_L[0];PrimFrac=Frac_p[0]*Frac_L[0];}
    if(uSpec==2){Purity=Purity_p[0]*Purity_Xim[0];PrimFrac=Frac_p[0]*Frac_Xim[0];}
    if(uSpec==4){Purity=Purity_L[0]*Purity_L[0];PrimFrac=Frac_L[0]*Frac_L[0];}
    printf("    A = %.3f +/- %.3f\n",Strength,StrengthErr);
    printf("  A/P = %.3f +/- %.3f\n",Strength/Purity,StrengthErr/Purity);
    printf("  A/PF = %.3f +/- %.3f\n",Strength/Purity/PrimFrac,StrengthErr/Purity/PrimFrac);
    printf("   r0 = %.3f +/- %.3f fm\n",fQuantumBaseline[uSpec]->GetParameter(3),fQuantumBaseline[uSpec]->GetParError(3));

    fClumsyFit[uSpec] = new TF1(TString::Format("fClumsyFit_%s",SpeciesName[uSpec].Data()),
                                "[0]*(1.+[1]*x)*(1.-[2]*exp(-x*x*[3]*[3]/197.327/197.327))",FitFrom[uSpec],FitTo[uSpec]);
    fClumsyFit[uSpec]->SetParameter(0,fQuantumBaseline[uSpec]->GetParameter(0));
    fClumsyFit[uSpec]->SetParameter(1,fQuantumBaseline[uSpec]->GetParameter(1));
    fClumsyFit[uSpec]->SetParameter(2,fQuantumBaseline[uSpec]->GetParameter(2));

    fBestClumsyFit[uSpec] = new TF1(TString::Format("fBestClumsyFit_%s",SpeciesName[uSpec].Data()),
                                "[0]*(1.+[1]*x)*(1.-[2]*exp(-x*x*[3]*[3]/197.327/197.327))",FitFrom[uSpec],FitTo[uSpec]);
  }

  const double Min_r0 = 0.20;
  const double Max_r0 = 0.40;
  const double Step_r0 = 0.01;
  double BestChi2 = 1e6;
  double Best_r0 = 0;
  for(double dr0=Min_r0; dr0<=Max_r0+1e-6; dr0+=Step_r0){
    double TotChi2 = 0;
    //printf("CLUMSY FITTER for r0=%.3f:\n",dr0);
    for(unsigned uSpec=0; uSpec<NumSpecies; uSpec++){
      //printf("    Clumsy fitting %s\n",SpeciesName[uSpec].Data());
      if(FitFrom[uSpec]==FitTo[uSpec]) continue;
      fClumsyFit[uSpec]->FixParameter(3,dr0);
      hCkPlot[uSpec]->Fit(fClumsyFit[uSpec],"Q, S, N, R, M");
      TotChi2 += fClumsyFit[uSpec]->GetChisquare();
      //printf("    chi2 = %.0f\n",fClumsyFit[uSpec]->GetChisquare());
    }
    //printf("  TotChi2 = %.0f\n",TotChi2);
    if(TotChi2<BestChi2) {
      BestChi2=TotChi2; Best_r0=dr0;
      for(unsigned uSpec=0; uSpec<NumSpecies; uSpec++){
        if(FitFrom[uSpec]==FitTo[uSpec]) continue;
        fBestClumsyFit[uSpec]->SetParameter(0,fClumsyFit[uSpec]->GetParameter(0));
        fBestClumsyFit[uSpec]->SetParError(0,fClumsyFit[uSpec]->GetParError(0));
        fBestClumsyFit[uSpec]->SetParameter(1,fClumsyFit[uSpec]->GetParameter(1));
        fBestClumsyFit[uSpec]->SetParError(1,fClumsyFit[uSpec]->GetParError(1));
        fBestClumsyFit[uSpec]->SetParameter(2,fClumsyFit[uSpec]->GetParameter(2));
        fBestClumsyFit[uSpec]->SetParError(2,fClumsyFit[uSpec]->GetParError(2));
        fBestClumsyFit[uSpec]->SetParameter(3,fClumsyFit[uSpec]->GetParameter(3));
        fBestClumsyFit[uSpec]->SetParError(3,fClumsyFit[uSpec]->GetParError(3));
      }
    }
  }
  printf("Best r0 = %.2f with chi2 = %.0f\n",Best_r0,BestChi2);
  for(unsigned uSpec=0; uSpec<NumSpecies; uSpec++){
    if(FitFrom[uSpec]==FitTo[uSpec]) continue;

    double Strength = fBestClumsyFit[uSpec]->GetParameter(2);
    double StrengthErr = fBestClumsyFit[uSpec]->GetParError(2);
    double Purity;
    double PrimFrac;
    if(uSpec==0){Purity=Purity_p[0]*Purity_p[0];PrimFrac=Frac_p[0]*Frac_p[0];}
    if(uSpec==1){Purity=Purity_p[0]*Purity_L[0];PrimFrac=Frac_p[0]*Frac_L[0];}
    if(uSpec==2){Purity=Purity_p[0]*Purity_Xim[0];PrimFrac=Frac_p[0]*Frac_Xim[0];}
    if(uSpec==4){Purity=Purity_L[0]*Purity_L[0];PrimFrac=Frac_L[0]*Frac_L[0];}
    printf("    A = %.3f +/- %.3f\n",Strength,StrengthErr);
    printf("  A/P = %.3f +/- %.3f\n",Strength/Purity,StrengthErr/Purity);
    printf("  A/PF = %.3f +/- %.3f\n",Strength/Purity/PrimFrac,StrengthErr/Purity/PrimFrac);

  }


  TFile* OutputFile = new TFile(OutputFileName,"recreate");
  OutputFile->cd();
  for(unsigned uSpec=0; uSpec<NumSpecies; uSpec++){
    if(hCkPlot[uSpec]) hCkPlot[uSpec]->Write();
    if(fQuantumBaseline[uSpec]) fQuantumBaseline[uSpec]->Write();
    if(fBestClumsyFit[uSpec]) fBestClumsyFit[uSpec]->Write();
  }

  for(unsigned uSpec=0; uSpec<NumSpecies; uSpec++){
    if(hCkPlot[uSpec]) delete hCkPlot[uSpec];
    if(fQuantumBaseline[uSpec]) delete fQuantumBaseline[uSpec];
    if(fClumsyFit[uSpec]) delete fClumsyFit[uSpec];
    if(fBestClumsyFit[uSpec]) delete fBestClumsyFit[uSpec];
  }
  delete [] hCkPlot;
  delete [] fQuantumBaseline;
  delete [] fClumsyFit;
  delete [] fBestClumsyFit;
  delete [] NormFrom;
  delete [] NormTo;
  printf("Ending QuantumBaseline\n");
}

void EPOS_QS_PLOT(){

  //strength QS_s of 0,0.2 and 0.4
  TString EposInputFile0 = TString::Format("%s/MixedEvents/dEta_dPhi_Ck_QS/QS0_pp_1.3eta.root",GetFemtoOutputFolder());
  TString EposInputFile2 = TString::Format("%s/MixedEvents/dEta_dPhi_Ck_QS/QS2_pp_1.3eta.root",GetFemtoOutputFolder());
  TString EposInputFile4 = TString::Format("%s/MixedEvents/dEta_dPhi_Ck_QS/QS4_pp_1.3eta.root",GetFemtoOutputFolder());
  TString HistoInputName = "hCphiQsRs";
  TGraphErrors Polish_HM_Data;

  Polish_HM_Data.SetPoint(0,-1.286,0.984);
  Polish_HM_Data.SetPoint(1,-1.077,0.965);
  Polish_HM_Data.SetPoint(2,-0.873,0.938);
  Polish_HM_Data.SetPoint(3,-0.637,0.909);
  Polish_HM_Data.SetPoint(4,-0.44,0.899);
  Polish_HM_Data.SetPoint(5,-0.212,0.888);
  Polish_HM_Data.SetPoint(6,0.0,0.903);
  Polish_HM_Data.SetPoint(7,0.22,0.889);
  Polish_HM_Data.SetPoint(8,0.44,0.896);
  Polish_HM_Data.SetPoint(9,0.64,0.911);
  Polish_HM_Data.SetPoint(10,0.86,0.934);
  Polish_HM_Data.SetPoint(11,1.08,0.958);
  Polish_HM_Data.SetPoint(12,1.3,0.985);
  Polish_HM_Data.SetPoint(13,1.519,1.006);
  Polish_HM_Data.SetPoint(14,1.739,1.026);
  Polish_HM_Data.SetPoint(15,1.96,1.047);
  Polish_HM_Data.SetPoint(16,2.172,1.056);
  Polish_HM_Data.SetPoint(17,2.424,1.069);
  Polish_HM_Data.SetPoint(18,2.613,1.084);
  Polish_HM_Data.SetPoint(19,2.817,1.092);
  Polish_HM_Data.SetPoint(20,3.053,1.092);
  Polish_HM_Data.SetPoint(21,3.234,1.092);
  Polish_HM_Data.SetPoint(22,3.478,1.087);
  Polish_HM_Data.SetPoint(23,3.691,1.084);
  Polish_HM_Data.SetPoint(24,3.895,1.073);
  Polish_HM_Data.SetPoint(25,4.115,1.057);
  Polish_HM_Data.SetPoint(26,4.328,1.045);
  Polish_HM_Data.SetPoint(27,4.54,1.026);
  Polish_HM_Data.SetPoint(28,4.772,1.007);

  Polish_HM_Data.SetPoint(0,-1.286,0.984);
  Polish_HM_Data.SetPoint(1,-1.077,0.965);
  Polish_HM_Data.SetPoint(2,-0.873,0.938);
  Polish_HM_Data.SetPoint(3,-0.637,0.909);
  Polish_HM_Data.SetPoint(4,-0.44,0.899);
  Polish_HM_Data.SetPoint(5,-0.212,0.888);
  Polish_HM_Data.SetPoint(6,0.008,0.903);
  Polish_HM_Data.SetPoint(7,0.213,0.889);
  Polish_HM_Data.SetPoint(8,0.425,0.896);
  Polish_HM_Data.SetPoint(9,0.661,0.911);
  Polish_HM_Data.SetPoint(10,0.874,0.934);
  Polish_HM_Data.SetPoint(11,1.094,0.958);
  Polish_HM_Data.SetPoint(12,1.299,0.985);
  Polish_HM_Data.SetPoint(13,1.519,1.006);
  Polish_HM_Data.SetPoint(14,1.739,1.026);
  Polish_HM_Data.SetPoint(15,1.96,1.047);
  Polish_HM_Data.SetPoint(16,2.172,1.056);
  Polish_HM_Data.SetPoint(17,2.424,1.069);
  Polish_HM_Data.SetPoint(18,2.613,1.084);
  Polish_HM_Data.SetPoint(19,2.817,1.092);
  Polish_HM_Data.SetPoint(20,3.053,1.092);
  Polish_HM_Data.SetPoint(21,3.234,1.092);
  Polish_HM_Data.SetPoint(22,3.478,1.087);
  Polish_HM_Data.SetPoint(23,3.691,1.084);
  Polish_HM_Data.SetPoint(24,3.895,1.073);
  Polish_HM_Data.SetPoint(25,4.115,1.057);
  Polish_HM_Data.SetPoint(26,4.328,1.045);
  Polish_HM_Data.SetPoint(27,4.54,1.026);
  Polish_HM_Data.SetPoint(28,4.772,1.007);

  for(unsigned uPts=0; uPts<Polish_HM_Data.GetN(); uPts++){
    Polish_HM_Data.SetPointError(uPts,0.22,0.015);
  }

  TFile* fEpos0 = new TFile(EposInputFile0,"read");
  TH1F* hEpos0 = (TH1F*)fEpos0->Get(HistoInputName);
  TGraph gEpos0;
  for(unsigned uBin=0; uBin<hEpos0->GetNbinsX(); uBin++){
    gEpos0.SetPoint(uBin,hEpos0->GetBinCenter(uBin+1),hEpos0->GetBinContent(uBin+1));
  }

  TFile* fEpos2 = new TFile(EposInputFile2,"read");
  TH1F* hEpos2 = (TH1F*)fEpos2->Get(HistoInputName);
  TGraph gEpos2;
  for(unsigned uBin=0; uBin<hEpos2->GetNbinsX(); uBin++){
    gEpos2.SetPoint(uBin,hEpos2->GetBinCenter(uBin+1),hEpos2->GetBinContent(uBin+1));
  }

  TFile* fEpos4 = new TFile(EposInputFile4,"read");
  TH1F* hEpos4 = (TH1F*)fEpos4->Get(HistoInputName);
  TGraph gEpos4;
  for(unsigned uBin=0; uBin<hEpos4->GetNbinsX(); uBin++){
    gEpos4.SetPoint(uBin,hEpos4->GetBinCenter(uBin+1),hEpos4->GetBinContent(uBin+1));
  }

  TH1F* hAxis = new TH1F("hAxis", "hAxis", 128, -1.4, 4.8);
  hAxis->SetStats(false);
  hAxis->SetTitle("");
  hAxis->GetXaxis()->SetTitle("#Delta#varphi");
  hAxis->GetXaxis()->SetTitleSize(0.06);
  hAxis->GetXaxis()->SetLabelSize(0.06);
  hAxis->GetXaxis()->SetTitleOffset(1.2);
  hAxis->GetXaxis()->SetLabelOffset(0.02);

  hAxis->GetYaxis()->SetTitle("C(#Delta#varphi)");
  hAxis->GetYaxis()->SetRangeUser(0.85, 1.15);
  hAxis->GetYaxis()->SetTitleSize(0.06);
  hAxis->GetYaxis()->SetLabelSize(0.06);
  hAxis->GetYaxis()->SetTitleOffset(0.90);

  Polish_HM_Data.SetName("Polish_HM_Data");
  Polish_HM_Data.SetFillColorAlpha(kBlue+3,0.3);
  Polish_HM_Data.SetLineColor(kBlue+3);
  Polish_HM_Data.SetLineWidth(0);

  gEpos0.SetName("gEpos0");
  gEpos0.SetLineColor(kBlack);
  gEpos0.SetLineWidth(8);
  hEpos0->SetLineColor(kBlack);
  hEpos0->SetLineWidth(3);
  hEpos0->SetMarkerColor(kBlack);

  gEpos2.SetName("gEpos2");
  gEpos2.SetLineColor(kRed+1);
  gEpos2.SetLineWidth(8);
  hEpos2->SetLineColor(kRed+1);
  hEpos2->SetLineWidth(3);
  hEpos2->SetMarkerColor(kRed+1);

  gEpos4.SetName("gEpos3");
  gEpos4.SetLineColor(kGreen+1);
  gEpos4.SetLineWidth(8);
  hEpos4->SetLineColor(kGreen+1);
  hEpos4->SetLineWidth(3);
  hEpos4->SetMarkerColor(kGreen+1);

  TLegend* lLegend = new TLegend(0.67,0.18,0.90,0.41);//lbrt
  lLegend->SetName(TString::Format("lLegend"));
  lLegend->SetTextSize(0.04);
  lLegend->SetBorderSize(0);
  lLegend->AddEntry(&Polish_HM_Data,"ALICE preliminary");
  lLegend->AddEntry(hEpos0,"EPOS default");
  lLegend->AddEntry(hEpos2,"EPOS A=0.2");
  lLegend->AddEntry(hEpos4,"EPOS A=0.4");

  TCanvas* cEpos = new TCanvas("cEpos", "cEpos", 1);
  cEpos->cd(0);
  cEpos->SetCanvasSize(1280, 720);
  cEpos->SetMargin(0.12,0.05,0.15,0.05);//lrbt
  cEpos->SetGrid();

  hAxis->Draw("axis");
  Polish_HM_Data.Draw("3 same");
  //gEpos4.Draw("C same");
  //gEpos0.Draw("C same");
  //gEpos2.Draw("C same");
  hEpos4->Draw("same");
  hEpos0->Draw("same");
  hEpos2->Draw("same");
  lLegend->Draw("same");

  cEpos->SaveAs(TString::Format("%s/Plots/PhDThesis/QuantumBaseline/cEpos_QS.png",GetCernBoxDimi()));

  delete lLegend;
  delete hEpos0;
  delete hEpos2;
  delete hEpos4;
  delete fEpos0;
  delete fEpos2;
  delete fEpos4;
}

//if par[0]==-1e6 : use splines
//else
//par[0] is an overall normalization
//than we have a pol4 = p0*(1+p1*k+p2*k^2+p3*k^4+p4*k^4), which has 3 free arguments and the following properties
//par4!=0 (pol4 flat at 0)
//	par1,par2 the two extrema, par3 is the p4, par4 is dummy
//par4==0&&par3!=0 (pol3)
//	par1,par2 the two extrema, par3 is p3
//par4==0&&par3==0&&par2!=0 (pol2)
//	par1 is the extrema, par2 is p2
//par4==0&&par3==0&&par2==0&&par1!=0 (pol1)
//	par1 is p1
//to avoid problems with a starting parameter of zero, to switch the order of the par we use -1e6 as a value
//a Mathematica computation of the equations is in your Femto folder
double ThesisFitterBl_pp(double* x, double* par){
    double& k = *x;
    double& p0 = par[0];
    //splines
    if(p0==-1e6){
        //double Spline = DLM_FITTER2_FUNCTION_SPLINE3(x,&par[1]);
        //printf("hi %f\n",Spline);
        return DLM_FITTER2_FUNCTION_SPLINE3(x,&par[1]);
    }
    //constrained polynomials
    else{
        double p1;
        double p2;
        double p3;
        double p4;
        if(par[4]!=-1e6){
            p4 = par[3];
            p3 = -4./3.*(par[1]+par[2])*p4;
            p2 = 2.*par[1]*par[2]*p4;
            p1 = 0;
        }
        else if(par[3]!=-1e6){
            p4 = 0;
            p3 = par[3];
            p2 = -1.5*(par[1]+par[2])*p3;
            p1 = 3.*par[1]*par[2]*p3;
        }
        else if(par[2]!=-1e6){
            p4 = 0;
            p3 = 0;
            p2 = par[2];
            p1 = -2.*par[1]*p2;
        }
        else{
            p4 = 0;
            p3 = 0;
            p2 = 0;
            p1 = par[1];
        }
        return p0*(1.+p1*k+p2*pow(k,2)+p3*pow(k,3)+p4*pow(k,4));
    }
}
double ThesisFitterFemto_pp(double* x, double* par){
    double& MOM = *x;
    CKDEC_FIT_PP->GetCk()->SetSourcePar(0,par[0]);
    CKDEC_FIT_PP->Update(true,true);
    return CKDEC_FIT_PP->EvalCk(MOM);
}
double ThesisFitter_pp(double* x, double* par){
    return ThesisFitterFemto_pp(x,par)*ThesisFitterBl_pp(x,&par[1]);
}
double SimplePol4_ppFit(double* x, double* par){
    double& MOM = *x;
    return ThesisFitterFemto_pp(x,par)*par[1]*(1.+par[2]*MOM+par[3]*MOM*MOM+par[4]*MOM*MOM*MOM+par[5]*MOM*MOM*MOM*MOM);
}


void Fit_pp(){
  //double* MomBins_pp = NULL;
  double* FitRegion_pp = NULL;
  unsigned NumMomBins_pp=100;
  double kMax = 400;
  TString DataSample = "pp13TeV_HM_Dec19";
  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());
  //AnalysisObject.SetUpBinning_pp(DataSample,NumMomBins_pp,MomBins_pp,FitRegion_pp);
  double lam_pp[5];
  double lam_pL[5];
  AnalysisObject.SetUpLambdaPars_pp(DataSample,0,lam_pp);
  AnalysisObject.SetUpLambdaPars_pL(DataSample,0,0,lam_pL);

  TH2F* hResolution_pp = AnalysisObject.GetResolutionMatrix(DataSample,"pp");
  TH2F* hResidual_pp_pL = AnalysisObject.GetResidualMatrix("pp","pLambda");

  TString SourceDescriptionG = "Gauss";
  TString SourceDescriptionR = "McGauss_ResoTM";

  CATS AB_pL;
  DLM_Ck* Ck_pL;
  AB_pL.SetMomBins(NumMomBins_pp,0,kMax);
  AnalysisObject.SetUpCats_pL(AB_pL,"Usmani","Gauss");
  AB_pL.SetAnaSource(0,1.25);
  AB_pL.KillTheCat();
  AB_pL.SetNotifications(CATS::nWarning);
  Ck_pL = new DLM_Ck(AB_pL.GetNumSourcePars(),0,AB_pL);
  Ck_pL->Update();
  DLM_CkDecomposition CkDec_pL("pLambda",2,*Ck_pL,NULL);
  CkDec_pL.AddContribution(2,lam_pL[1]+lam_pL[2]+lam_pL[3],DLM_CkDecomposition::cFeedDown);
  CkDec_pL.AddContribution(3,lam_pL[4],DLM_CkDecomposition::cFake);//0.03

  CATS ABG_pp;
  DLM_Ck* CkG_pp;
  ABG_pp.SetMomBins(NumMomBins_pp,0,kMax);
  AnalysisObject.SetUpCats_pp(ABG_pp,"AV18",SourceDescriptionG,0,0);

  ABG_pp.SetAnaSource(0,1.28);
  if(SourceDescriptionG.Contains("Mc")){
      ABG_pp.SetAnaSource(0,1.20);
      ABG_pp.SetAnaSource(1,2.0);
  }
  //AB_pp.SetAnaSource(1,2.0);
  //AB_pp.SetNotifications(CATS::nWarning);
  ABG_pp.SetEpsilonConv(4e-8);
  ABG_pp.SetEpsilonProp(4e-8);
  ABG_pp.KillTheCat();
  ABG_pp.SetNotifications(CATS::nWarning);
  CkG_pp = new DLM_Ck(ABG_pp.GetNumSourcePars(),0,ABG_pp);
  CkG_pp->Update();

  CATS ABR_pp;
  DLM_Ck* CkR_pp;
  ABR_pp.SetMomBins(NumMomBins_pp,0,kMax);
  //AnalysisObject.SetUpCats_pp(AB_pp,"AV18","Gauss",0,0);//McLevyNolan_Reso
  AnalysisObject.SetUpCats_pp(ABR_pp,"AV18",SourceDescriptionR,0,202);//McLevyNolan_Reso

  ABR_pp.SetAnaSource(0,1.3);
  if(SourceDescriptionR.Contains("Mc")){
      ABR_pp.SetAnaSource(0,1.20);
      ABR_pp.SetAnaSource(1,2.0);
  }
  //AB_pp.SetAnaSource(1,2.0);
  //AB_pp.SetNotifications(CATS::nWarning);
  ABR_pp.SetEpsilonConv(4e-8);
  ABR_pp.SetEpsilonProp(4e-8);
  ABR_pp.KillTheCat();
  ABR_pp.SetNotifications(CATS::nWarning);
  CkR_pp = new DLM_Ck(ABR_pp.GetNumSourcePars(),0,ABR_pp);
  CkR_pp->Update();

  TH1F* hData_pp = AnalysisObject.GetAliceExpCorrFun(DataSample,"pp","_0",0,false,-1);

  TString OutputFolder = TString::Format("%s/Plots/PhDThesis/ppFits/",GetCernBoxDimi());
  TString OutFileName = "ppFit.root";
  TFile* OutputFile = new TFile(OutputFolder+OutFileName,"recreate");
  hData_pp->Write();
/*
  DLM_CkDecomposition CkDecG_pp("pp",3,*CkG_pp,hResolution_pp);

  CkDecG_pp.AddContribution(0,lam_pp[1],DLM_CkDecomposition::cFeedDown,&CkDec_pL,hResidual_pp_pL);
  CkDecG_pp.AddContribution(1,lam_pp[2],DLM_CkDecomposition::cFeedDown);
  CkDecG_pp.AddContribution(2,lam_pp[3],DLM_CkDecomposition::cFake);

  CkDecG_pp.Update();
  CkDec_pL.Update();

  DLM_CkDecomposition CkDecR_pp("pp",3,*CkR_pp,hResolution_pp);

  CkDecR_pp.AddContribution(0,lam_pp[1],DLM_CkDecomposition::cFeedDown,&CkDec_pL,hResidual_pp_pL);
  CkDecR_pp.AddContribution(1,lam_pp[2],DLM_CkDecomposition::cFeedDown);
  CkDecR_pp.AddContribution(2,lam_pp[3],DLM_CkDecomposition::cFake);

  CkDecR_pp.Update();
  CkDec_pL.Update();

  const double FitMin = 0;
  const double FitMax = 376;

  CKDEC_FIT_PP = &CkDecG_pp;
  OutputFile->cd();
  TF1* fitG_pp = new TF1("fitG_pp",SimplePol4_ppFit,FitMin,FitMax,6);//ThesisFitter_pp,SimplePol4_ppFit
  fitG_pp->SetParameter(0,1.25);
  fitG_pp->SetParLimits(0,1.1,1.4);
  fitG_pp->SetParameter(1,1.0);//norm
  fitG_pp->SetParLimits(1,0.8,1.2);
  fitG_pp->SetParameter(2,50);//extremum 1
  fitG_pp->SetParLimits(2,-100000,200);
  fitG_pp->SetParameter(3,300);//extremum 2
  fitG_pp->SetParLimits(3,0,600);
  fitG_pp->SetParameter(4,0);//pol3
  fitG_pp->SetParLimits(4,-1e-6,1e6);
  fitG_pp->FixParameter(5,0);


  //fitG_pp->FixParameter(2,0);
  fitG_pp->SetParameter(2,0);
  fitG_pp->SetParLimits(2,-1e-3,1e-3);

  //fitG_pp->FixParameter(3,0);
  fitG_pp->SetParameter(3,0);
  fitG_pp->SetParLimits(3,-1e-5,1e-5);

  //fitG_pp->FixParameter(4,0);
  fitG_pp->SetParameter(4,0);
  fitG_pp->SetParLimits(4,-1e-7,1e-7);

  //fitG_pp->FixParameter(5,0);
  fitG_pp->SetParameter(5,0);
  fitG_pp->SetParLimits(5,-1e-9,1e-9);


  //hData_pp->Fit(fitG_pp,"S, N, R, M");
  printf("Gaussian source %.3f+/-%.3f fm\n",fitG_pp->GetParameter(0),fitG_pp->GetParError(0));
  printf(" chi2/ndf = %.0f/%i = %.2f\n",fitG_pp->GetChisquare(),fitG_pp->GetNDF(),fitG_pp->GetChisquare()/double(fitG_pp->GetNDF()));

  CKDEC_FIT_PP = &CkDecR_pp;
  OutputFile->cd();
  TF1* fitR_pp = new TF1("fitR_pp",SimplePol4_ppFit,FitMin,FitMax,6);//ThesisFitter_pp,SimplePol4_ppFit
  fitR_pp->SetParameter(0,1.0);
  fitR_pp->SetParLimits(0,0.9,1.1);
  fitR_pp->SetParameter(1,1.0);//norm
  fitR_pp->SetParLimits(1,0.8,1.2);
  fitR_pp->SetParameter(2,50);//extremum 1
  fitR_pp->SetParLimits(2,-100000,200);
  fitR_pp->SetParameter(3,300);//extremum 2
  fitR_pp->SetParLimits(3,0,600);
  fitR_pp->SetParameter(4,0);//pol3
  fitR_pp->SetParLimits(4,-1e-6,1e6);
  fitR_pp->FixParameter(5,0);

  fitR_pp->FixParameter(2,0);
  fitR_pp->FixParameter(3,0);
  fitR_pp->FixParameter(4,0);
  fitR_pp->FixParameter(5,0);

  //hData_pp->Fit(fitR_pp,"S, N, R, M");
  printf("Gaussian core %.3f+/-%.3f fm\n",fitR_pp->GetParameter(0),fitR_pp->GetParError(0));
*/
  TGraph gTheoG;
  gTheoG.SetName("gTheoG");
  for(unsigned uBin=0; uBin<ABG_pp.GetNumMomBins(); uBin++){
    gTheoG.SetPoint(uBin,ABG_pp.GetMomentum(uBin),ABG_pp.GetCorrFun(uBin));
  }

  TGraph gTheoR;
  gTheoR.SetName("gTheoR");
  for(unsigned uBin=0; uBin<ABR_pp.GetNumMomBins(); uBin++){
    gTheoR.SetPoint(uBin,ABR_pp.GetMomentum(uBin),ABR_pp.GetCorrFun(uBin));
  }

  TGraph gTheoGR;
  gTheoGR.SetName("gTheoGR");
  for(unsigned uBin=0; uBin<ABR_pp.GetNumMomBins(); uBin++){
    gTheoGR.SetPoint(uBin,ABR_pp.GetMomentum(uBin),ABG_pp.GetCorrFun(uBin)/ABR_pp.GetCorrFun(uBin));
  }

  gTheoG.SetLineColor(kBlue);
  gTheoG.SetLineWidth(6);
  gTheoR.SetLineColor(kRed);
  gTheoR.SetLineWidth(4);
  gTheoR.SetLineStyle(2);
  gTheoGR.SetLineColor(kBlack);
  gTheoGR.SetLineWidth(6);

  TLegend* legend = new TLegend(0.45,0.5,0.95,0.9);//lbrt
  legend->SetBorderSize(0);
  legend->SetTextFont(42);
  legend->SetTextSize(gStyle->GetTextSize()*1.3);
  legend->AddEntry(&gTheoG, "Gaussian source");
  legend->AddEntry(&gTheoR, "Resonance source model");
  legend->AddEntry(&gTheoGR, "Ratio");

  DLM_SubPads* DlmPad = new DLM_SubPads(720,1080);
  DlmPad->AddSubPad(0,1,0.67,1);
  DlmPad->AddSubPad(0,1,0.34,0.67);
  DlmPad->AddSubPad(0,1,0.0,0.34);
  DlmPad->SetMargin(0,0.14,0.02,0.0,0.02);
  DlmPad->SetMargin(1,0.14,0.02,0.0,0.0);
  DlmPad->SetMargin(2,0.14,0.02,0.05,0.0);
  DlmPad->cd(0);

  TH1F* hCkAxis = new TH1F("hCkAxis","hCkAxis",100,0,400);
  hCkAxis->SetStats(false);
  hCkAxis->SetTitle("; #it{k*} (MeV/#it{c}); #it{C}(#it{k*})");
  hCkAxis->GetXaxis()->SetRangeUser(0, 456);
  hCkAxis->GetXaxis()->SetNdivisions(505);
  hCkAxis->GetYaxis()->SetRangeUser(0.55, 4.5);
  hCkAxis->SetFillColor(kGray+1);
  SetStyleHisto2a(hCkAxis,2,0);
  hCkAxis->GetYaxis()->SetLabelSize(gStyle->GetTextSize()*1.3);
  hCkAxis->GetYaxis()->SetTitleSize(gStyle->GetTextSize()*1.3);


  hCkAxis->Draw("axis");
  gTheoG.Draw("L,same");
  gTheoR.Draw("L,same");
  legend->Draw("same");

  DlmPad->cd(1);
  TH1F* hCkAxisInlet = new TH1F("hCkAxisInlet","hCkAxisInlet",100,0,400);
  hCkAxisInlet->SetStats(false);
  hCkAxisInlet->SetTitle("; #it{k*} (MeV/#it{c}); #it{C}(#it{k*})");
  hCkAxisInlet->GetXaxis()->SetRangeUser(0, 456);
  hCkAxisInlet->GetXaxis()->SetNdivisions(505);
  hCkAxisInlet->GetYaxis()->SetRangeUser(0.935, 1.005);
  hCkAxisInlet->SetFillColor(kGray+1);
  SetStyleHisto2a(hCkAxisInlet,2,0);
  hCkAxisInlet->GetYaxis()->SetLabelSize(gStyle->GetTextSize()*1.3);
  hCkAxisInlet->GetYaxis()->SetTitleSize(gStyle->GetTextSize()*1.3);

  hCkAxisInlet->Draw("axis");
  gTheoG.Draw("L,same");
  gTheoR.Draw("L,same");

  DlmPad->cd(2);

  TH1F* hRatioAxis = new TH1F("hRatioAxis","hRatioAxis",100,0,400);
  hRatioAxis->SetStats(false);
  hRatioAxis->SetTitle("; #it{k*} (MeV/#it{c}); Ratio");
  hRatioAxis->GetXaxis()->SetRangeUser(0, 456);
  hRatioAxis->GetXaxis()->SetNdivisions(505);
  hRatioAxis->GetYaxis()->SetRangeUser(0.990, 1.014);
  hRatioAxis->SetFillColor(kGray+1);
  SetStyleHisto2a(hRatioAxis,2,0);
  hRatioAxis->GetXaxis()->SetLabelSize(gStyle->GetTextSize()*1.3);
  hRatioAxis->GetXaxis()->SetTitleSize(gStyle->GetTextSize()*1.3);
  hRatioAxis->GetYaxis()->SetLabelSize(gStyle->GetTextSize()*1.3);
  hRatioAxis->GetYaxis()->SetTitleSize(gStyle->GetTextSize()*1.3);

  hRatioAxis->Draw("axis");
  gTheoGR.Draw("L,same");

  DlmPad->GetCanvas()->SaveAs(OutputFolder+"JuliHelpsWithResonances.pdf");

  OutputFile->cd();
  //fitG_pp->Write();
  //fitR_pp->Write();
  gTheoG.Write();
  gTheoR.Write();
  gTheoGR.Write();


  delete DlmPad;
}


void Compare_pL_models(){

  TString OutputFolder = GetCernBoxDimi();
  OutputFolder += "/Plots/PhDThesis/pLambda/";

  CATS Kitty_NLO19_600;
  CATS Kitty_NLO19_500;
  CATS Kitty_NLO13_600;
  CATS Kitty_LO13_600;
  CATS Kitty_pSo_Feed;

  const unsigned NumMomBins = 90;
  const double kMax = 360;

  Kitty_NLO19_600.SetMomBins(NumMomBins,0,kMax);
  Kitty_NLO19_500.SetMomBins(NumMomBins,0,kMax);
  Kitty_NLO13_600.SetMomBins(NumMomBins,0,kMax);
  Kitty_LO13_600.SetMomBins(NumMomBins,0,kMax);
  Kitty_pSo_Feed.SetMomBins(NumMomBins,0,kMax);

  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()));
  AnalysisObject.SetUpCats_pL(Kitty_NLO19_600,"Chiral_Coupled_SPD","McGauss_ResoTM",11600,202);
  AnalysisObject.SetUpCats_pL(Kitty_NLO19_500,"Chiral_Coupled_SPD","McGauss_ResoTM",11500,202);
  AnalysisObject.SetUpCats_pL(Kitty_NLO13_600,"Chiral_Coupled_SPD","McGauss_ResoTM",1600,202);
  AnalysisObject.SetUpCats_pL(Kitty_LO13_600,"Chiral_Coupled_SPD","McGauss_ResoTM",-11600,202);
  AnalysisObject.SetUpCats_pS0(Kitty_pSo_Feed,"Chiral","Gauss",0,0);

  const double CuspWeight = 0.33;
  const double CoreRadius = 1.06;
  const double SigmaRadius = 1.25;
  const double lambda_pL = 0.5;
  const double lambda_pSo = 0.2;

  Kitty_NLO19_600.SetChannelWeight(7,1./4.*CuspWeight);//1S0 SN(s) -> LN(s)
  Kitty_NLO19_600.SetChannelWeight(8,3./4.*CuspWeight);//3S1 SN(s) -> LN(s)
  Kitty_NLO19_600.SetChannelWeight(10,3./4.*CuspWeight);//3S1 SN(d) -> LN(s)
  Kitty_NLO19_600.SetChannelWeight(13,3./20.*CuspWeight);//3D1 SN(d) -> LN(d)
  Kitty_NLO19_600.SetChannelWeight(15,3./20.*CuspWeight);//3D1 SN(s) -> LN(d)

  Kitty_NLO19_500.SetChannelWeight(7,1./4.*CuspWeight);//1S0 SN(s) -> LN(s)
  Kitty_NLO19_500.SetChannelWeight(8,3./4.*CuspWeight);//3S1 SN(s) -> LN(s)
  Kitty_NLO19_500.SetChannelWeight(10,3./4.*CuspWeight);//3S1 SN(d) -> LN(s)
  Kitty_NLO19_500.SetChannelWeight(13,3./20.*CuspWeight);//3D1 SN(d) -> LN(d)
  Kitty_NLO19_500.SetChannelWeight(15,3./20.*CuspWeight);//3D1 SN(s) -> LN(d)

  Kitty_NLO13_600.SetChannelWeight(7,1./4.*CuspWeight);//1S0 SN(s) -> LN(s)
  Kitty_NLO13_600.SetChannelWeight(8,3./4.*CuspWeight);//3S1 SN(s) -> LN(s)
  Kitty_NLO13_600.SetChannelWeight(10,3./4.*CuspWeight);//3S1 SN(d) -> LN(s)
  Kitty_NLO13_600.SetChannelWeight(13,3./20.*CuspWeight);//3D1 SN(d) -> LN(d)
  Kitty_NLO13_600.SetChannelWeight(15,3./20.*CuspWeight);//3D1 SN(s) -> LN(d)

  Kitty_LO13_600.SetChannelWeight(7,1./4.*CuspWeight);//1S0 SN(s) -> LN(s)
  Kitty_LO13_600.SetChannelWeight(8,3./4.*CuspWeight);//3S1 SN(s) -> LN(s)
  Kitty_LO13_600.SetChannelWeight(10,3./4.*CuspWeight);//3S1 SN(d) -> LN(s)
  Kitty_LO13_600.SetChannelWeight(13,3./20.*CuspWeight);//3D1 SN(d) -> LN(d)
  Kitty_LO13_600.SetChannelWeight(15,3./20.*CuspWeight);//3D1 SN(s) -> LN(d)

  Kitty_NLO19_600.SetAnaSource(0,CoreRadius);//c.a. 10% smaller compared to p-p due to the mT scaling
  Kitty_NLO19_600.SetAnaSource(1,2.0);
  Kitty_NLO19_600.KillTheCat();

  Kitty_NLO19_500.SetAnaSource(0,CoreRadius);//c.a. 10% smaller compared to p-p due to the mT scaling
  Kitty_NLO19_500.SetAnaSource(1,2.0);
  Kitty_NLO19_500.KillTheCat();

  Kitty_NLO13_600.SetAnaSource(0,CoreRadius);//c.a. 10% smaller compared to p-p due to the mT scaling
  Kitty_NLO13_600.SetAnaSource(1,2.0);
  Kitty_NLO13_600.KillTheCat();

  Kitty_LO13_600.SetAnaSource(0,CoreRadius);//c.a. 10% smaller compared to p-p due to the mT scaling
  Kitty_LO13_600.SetAnaSource(1,2.0);
  Kitty_LO13_600.KillTheCat();

  Kitty_pSo_Feed.SetAnaSource(0,SigmaRadius);
  Kitty_pSo_Feed.KillTheCat();

  DLM_Ck Ck_pL(Kitty_NLO19_600.GetNumSourcePars(),0,Kitty_NLO19_600,NumMomBins,0,kMax);
  Ck_pL.Update();
  DLM_Ck Ck_pS0(Kitty_pSo_Feed.GetNumSourcePars(),0,Kitty_pSo_Feed,NumMomBins,0,kMax);
  Ck_pS0.Update();
  DLM_CkDecomposition CkDec_pS0("pSigma0",0,Ck_pS0,NULL);
  TH2F* hResidual_pL_pSigma0 = AnalysisObject.GetResidualMatrix("pLambda","pSigma0");

  //TH2F* hResolution_pL = AnalysisObject.GetResolutionMatrix("pp13TeV_HM_Dec19","pLambda");
  TH1F* hData = AnalysisObject.GetAliceExpCorrFun("pp13TeV_HM_Dec19","pLambda","_0",2,false,-1);
  DLM_CkDecomposition CkDec_pL("pLambda",2,Ck_pL,NULL);
  CkDec_pL.AddContribution(0,lambda_pSo,DLM_CkDecomposition::cFeedDown,&CkDec_pS0,hResidual_pL_pSigma0);
  CkDec_pL.AddContribution(1,1.-lambda_pSo,DLM_CkDecomposition::cFeedDown);
  CkDec_pS0.Update();
  CkDec_pL.Update();

  int ColorNLO13 = kRed;//+1
  int ColorNLO19 = kCyan;
  int ColorLO13 = kGreen;//+1
  int ColorNLO19_500 = kCyan+2;
  int ColorSigma = kAzure;

  gStyle->SetCanvasPreferGL(1);
  SetStyle();

  TGraphErrors gNLO19_600;
  gNLO19_600.SetName("gNLO19_600");
  gNLO19_600.SetLineColor(ColorNLO19);
  gNLO19_600.SetLineWidth(0);
  gNLO19_600.SetLineStyle(1);
  gNLO19_600.SetFillColorAlpha(ColorNLO19,0.60);

  TGraphErrors gNLO13_600;
  gNLO13_600.SetName("gNLO13_600");
  gNLO13_600.SetLineColor(ColorNLO13);
  gNLO13_600.SetLineWidth(3);
  gNLO13_600.SetLineStyle(1);
  gNLO13_600.SetFillColorAlpha(ColorNLO13,0.60);

  TGraphErrors gLO13_600;
  gLO13_600.SetName("gLO13_600");
  gLO13_600.SetLineColor(ColorLO13);
  gLO13_600.SetLineWidth(3);
  gLO13_600.SetLineStyle(1);
  gLO13_600.SetFillColorAlpha(ColorLO13,0.60);

  TGraphErrors gNLO19_500;
  gNLO19_500.SetName("gNLO19_500");
  gNLO19_500.SetLineColor(ColorNLO19_500);
  gNLO19_500.SetLineWidth(3);
  gNLO19_500.SetLineStyle(2);
  gNLO19_500.SetFillColorAlpha(ColorNLO19_500,0.60);

  TGraphErrors gSigmaFeed;
  gSigmaFeed.SetName("gSigmaFeed");
  gSigmaFeed.SetLineColor(ColorSigma);
  gSigmaFeed.SetLineWidth(3);
  gSigmaFeed.SetLineStyle(3);
  gSigmaFeed.SetFillColorAlpha(ColorSigma,0.60);


  for(unsigned uBin=0; uBin<NumMomBins; uBin++){
    double MOM = Kitty_NLO19_600.GetMomentum(uBin);
    gNLO19_600.SetPoint(uBin,MOM,lambda_pL*Kitty_NLO19_600.GetCorrFun(uBin)+1.-lambda_pL);
    gNLO19_500.SetPoint(uBin,MOM,lambda_pL*Kitty_NLO19_500.GetCorrFun(uBin)+1.-lambda_pL);
    gNLO13_600.SetPoint(uBin,MOM,lambda_pL*Kitty_NLO13_600.GetCorrFun(uBin)+1.-lambda_pL);
    gLO13_600.SetPoint(uBin,MOM,lambda_pL*Kitty_LO13_600.GetCorrFun(uBin)+1.-lambda_pL);
    gSigmaFeed.SetPoint(uBin,MOM,CkDec_pL.EvalSignalChild(0,MOM)+1.);

    gNLO19_600.SetPointError(uBin,0,hData->GetBinError(hData->FindBin(MOM)));
    gNLO19_500.SetPointError(uBin,0,hData->GetBinError(hData->FindBin(MOM)));
    gNLO13_600.SetPointError(uBin,0,hData->GetBinError(hData->FindBin(MOM)));
    gLO13_600.SetPointError(uBin,0,hData->GetBinError(hData->FindBin(MOM)));
    gSigmaFeed.SetPointError(uBin,0,hData->GetBinError(hData->FindBin(MOM)));

  }


  TLegend* myLegend = new TLegend(0.50,0.50,0.95,0.95);//lbrt
  myLegend->SetName(TString::Format("myLegend"));
  myLegend->SetTextSize(0.07);
  myLegend->AddEntry(&gNLO19_600,"p-#Lambda (NLO19-600)");
  myLegend->AddEntry(&gNLO19_500,"p-#Lambda (NLO19-500)");
  myLegend->AddEntry(&gNLO13_600,"p-#Lambda (NLO13-600)");
  myLegend->AddEntry(&gLO13_600,"p-#Lambda (LO-600)");
  myLegend->AddEntry(&gSigmaFeed,"p-#Sigma^{0}#rightarrowp-#Lambda (NLO19)");

  DLM_SubPads* DlmPad = new DLM_SubPads(1080,1080);
  DlmPad->AddSubPad(0,1,0.5,1.0);
  DlmPad->AddSubPad(0,1,0.0,0.5);
  DlmPad->SetMargin(0,0.14,0.02,0.0,0.02);
  DlmPad->SetMargin(1,0.14,0.02,0.09,0.0);

  TH1F* hAxisFull = new TH1F("hAxisFull","hAxisFull",NumMomBins,0,kMax);
  hAxisFull->SetTitle("; #it{k*} (MeV/#it{c}); #it{C}(#it{k*})");
  hAxisFull->GetXaxis()->SetRangeUser(0, kMax);
  hAxisFull->GetXaxis()->SetNdivisions(505);
  hAxisFull->GetYaxis()->SetRangeUser(0.9, 2.1);
  hAxisFull->SetFillColor(kGray+1);
  hAxisFull->SetStats(false);
  SetStyleHisto2(hAxisFull,2,0);
  SetStyleAxis_2(hAxisFull);

  TH1F* hAxisZoom = new TH1F("hAxisZoom","hAxisZoom",NumMomBins,0,kMax);
  hAxisZoom->SetTitle("; #it{k*} (MeV/#it{c}); #it{C}(#it{k*})");
  hAxisZoom->GetXaxis()->SetRangeUser(0, kMax);
  hAxisZoom->GetXaxis()->SetNdivisions(505);
  hAxisZoom->GetYaxis()->SetRangeUser(0.96, 1.04);
  hAxisZoom->SetFillColor(kGray+1);
  hAxisZoom->SetStats(false);
  SetStyleHisto2(hAxisZoom,2,0);
  SetStyleAxis_2(hAxisZoom);

  DlmPad->cd(0);
  hAxisFull->Draw("axis");
  gNLO19_600.Draw("3L,same");
  gNLO13_600.Draw("3L,same");
  gLO13_600.Draw("3L,same");
  gNLO19_500.Draw("3L,same");
  gSigmaFeed.Draw("3L,same");
  myLegend->Draw("same");

  DlmPad->cd(1);
  hAxisZoom->Draw("axis");
  gNLO19_600.Draw("3L,same");
  gNLO13_600.Draw("3L,same");
  gLO13_600.Draw("3L,same");
  gNLO19_500.Draw("3L,same");
  gSigmaFeed.Draw("3L,same");

  DlmPad->GetCanvas()->SaveAs(OutputFolder+"pL_ChiralModels.pdf");

  delete hAxisFull;
  delete hAxisZoom;

}


int THESIS_PLOTS(int narg, char** ARGS){
    ToyPotentials_PS_WF_CF();
    //ComparePionPion(0);
    //ComparePionPion(1);
    //Plot_ProtonLambda(Basics_ProtonLambda());
    //PlotCompareUsmaniLednicky();
    //RepulsiveCore_AppendixExample();
    //Ck_pp_Decomposition();
    //ModifiedPoisson();
    //pL_Feed();
    //pp_MomReso();
    //QuantumBaseline();
    //Fit_pp();
    //EPOS_QS_PLOT();
    //Compare_pL_models();

    return 0;
}
