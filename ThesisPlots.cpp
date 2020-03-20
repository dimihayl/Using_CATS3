#include "ThesisPlots.h"

#include <stdlib.h>

#include "CATS.h"
#include "CATSconstants.h"
#include "DLM_Potentials.h"
#include "DLM_Source.h"
#include "DLM_SubPads.h"

#include "TString.h"
#include "TGraph.h"
#include "TH1F.h"
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
    const double kMax = 400;
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
    RADIUS[2] = 3.0;

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
    TGraph* gWFd = new TGraph [NumPot];//diff to ref
    TGraph* gWFda = new TGraph [NumPot];//diff to as
    TGraph* gTotWF = new TGraph [NumPot];

    for(unsigned uPot=0; uPot<NumPot; uPot++){
        gWFref[uPot].SetName("gWFref_"+PotNames[uPot]);
        gWFas[uPot].SetName("gWFas_"+PotNames[uPot]);
        gWF[uPot].SetName("gWF_"+PotNames[uPot]);
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
                    //ParsPot.SetParameter(4,0.0);
                    //ParsPot.SetParameter(5,100.0);
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
                    if(uPot==0){
                        for(unsigned uFWF=0; uFWF<NumFreeWF; uFWF++){
                            gFreeWF[uFWF].SetPoint(uRadBin,RAD,real(Kitty.EvalReferenceRadialWF(Kitty.GetMomBin(FreeWaveK[uFWF]),0,RAD,true)));
                        }
                    }
                    gWFas[uPot].SetPoint(uRadBin,RAD,real(Kitty.EvalAsymptoticRadialWF(Kitty.GetMomBin(50),0,0,RAD,true)));
                    gWF[uPot].SetPoint(uRadBin,RAD,real(Kitty.EvalRadialWaveFunction(Kitty.GetMomBin(50),0,0,RAD,true)));
                    gWFd[uPot].SetPoint(uRadBin,RAD,real(Kitty.EvalRadialWaveFunction(Kitty.GetMomBin(50),0,0,RAD,true))-real(Kitty.EvalReferenceRadialWF(Kitty.GetMomBin(50),0,RAD,true)));
                    gWFda[uPot].SetPoint(uRadBin,RAD,real(Kitty.EvalRadialWaveFunction(Kitty.GetMomBin(50),0,0,RAD,true))-real(Kitty.EvalAsymptoticRadialWF(Kitty.GetMomBin(50),0,0,RAD,true)));
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

/*
        TH1F* hPhaseShift0 = new TH1F("hPhaseShift0", "hPhaseShift0", NumFineSteps, kMinFine-0.5*kFineStep, kMaxFine+0.5*kFineStep);
        for(unsigned iBin=0; iBin<NumFineSteps; iBin++){
            hPhaseShift0->SetBinContent(iBin+1, NalaFine.GetMomentum(iBin)/tan(NalaFine.GetPhaseShift(iBin,0,0)));
            //printf("k = %.2f; delta = %.3f; whatever/tan=%.3f\n",NalaFine.GetMomentum(iBin),NalaFine.GetPhaseShift(iBin,0,0)*180./3.14,
            //       NalaFine.GetMomentum(iBin)/tan(NalaFine.GetPhaseShift(iBin,0,0)));
            hPhaseShift0->SetBinError(iBin+1, hPhaseShift0->GetBinContent(iBin+1)*0.0001);
        }


        TF1* fPhaseShift0 = new TF1("fPhaseShift0", "[2]*x*x*x*x+0.5*[1]*x*x-1./[0]", kMinFine, kMaxFine);
*/

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
    hWFtot_Dummy->GetYaxis()->SetRangeUser(0,WFtot_yMax);
    TF1* fWFtot_base = new TF1("fWFtot_base","1",0,WF_rMax);
    fWFtot_base->SetLineColor(kBlack);
    fWFtot_base->SetLineWidth(1.5);
    fWFtot_base->SetLineStyle(2);

    //DLM_SubPads Pad_WF(1080,1920);
    DLM_SubPads Pad_WF(1280,1280);
    Pad_WF.AddSubPad(0,1,0.5,1);
    //Pad_WF.AddSubPad(0,1,0.38,0.62);
    Pad_WF.AddSubPad(0,1,0,0.5);

    Pad_WF.SetMargin(0,0.12,0.02,0.0,0.07);
    //Pad_WF.SetMargin(1,0.12,0.02,0.,0.);
    Pad_WF.SetMargin(1,0.12,0.02,0.09,0.0);

    Pad_WF.cd(0);

    hWFu_Dummy->Draw("axis");
    fWFu_base->Draw("same");
    gWFref[0].SetLineWidth(3);
    gWFref[0].SetLineColor(kGray+1);
    gWFref[0].SetLineStyle(4);
    gWFref[0].Draw("same");
    for(unsigned uPot=0; uPot<NumPot; uPot++){
        //gWF[uPot].Draw(uPot?"same":"");
        gWF[uPot].SetLineWidth(4.-0.33*double(uPot));
        gWF[uPot].SetLineColor(PotColors[uPot]);
        gWF[uPot].Draw("same");

        gWFas[uPot].SetLineWidth(4.-0.33*double(uPot));
        gWFas[uPot].SetLineColor(PotColors[uPot]);
        gWFas[uPot].SetLineStyle(3);
        gWFas[uPot].Draw("same");
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

    //Pad_WF.cd(1);
    //hWFda_Dummy->Draw("axis");
    //for(unsigned uPot=0; uPot<NumPot; uPot++){
    //    //gWFas[uPot].Draw(uPot?"same":"");
    //    gWFda[uPot].SetLineWidth(4.-0.33*double(uPot));
    //    gWFda[uPot].SetLineColor(PotColors[uPot]);
    //    gWFda[uPot].Draw("same");
    //}

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



int THESIS_PLOTS(int narg, char** ARGS){
    ToyPotentials_PS_WF_CF();
}
