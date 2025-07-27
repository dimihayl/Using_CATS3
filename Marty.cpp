#include "Marty.h"

#include "EnvVars.h"


#include <iostream>
#include <stdio.h>
#include <string.h>

#include <unistd.h>

#include "TH1F.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLegend.h"

#include "DLM_Source.h"
#include "DLM_RootFit.h"
#include "CATS.h"
#include "CommonAnaFunctions.h"

using namespace std;

void test_mem_a(float a, float* b, float& c){
    printf("a = %f\n",a);
    printf("&a = %p\n",&a);
    printf("b = %p\n",b);
    printf("b[0] = %f\n",b[0]);
    printf("c = %f\n",c);
    printf("&c = %p\n",&c);

    a = 11;
    c = 15;
    printf("\n");  
    printf("a = %f\n",a);
    printf("&a = %p\n",&a);
    printf("b = %p\n",b);
    printf("b[0] = %f\n",b[0]);
    printf("c = %f\n",c);
    printf("&c = %p\n",&c);
    printf("end \n");    
}

void test_mem(){
    float* pointer1 = new float[10];
    float number1;

    pointer1[0] = 1;
    pointer1[1] = 2;
    pointer1[2] = 3;

    printf("pointer1 = %p\n",pointer1);
    printf("pointer1[0] = %f\n",pointer1[0]);
    printf("&pointer1[0] = %p\n",&pointer1[0]);
    printf("pointer1[1] = %f\n",pointer1[1]);
    printf("&pointer1[1] = %p\n",&pointer1[1]);
    printf("pointer1[2] = %f\n",pointer1[2]);
    printf("pointer1[2] = %p\n",&pointer1[2]);
    //printf("&number1 = %p\n",&number1);

    test_mem_a(pointer1[0],pointer1,pointer1[2]);


    printf("pointer1 = %p\n",pointer1);
    printf("pointer1[0] = %f\n",pointer1[0]);
    printf("&pointer1[0] = %p\n",&pointer1[0]);
    printf("pointer1[1] = %f\n",pointer1[1]);
    printf("&pointer1[1] = %p\n",&pointer1[1]);
    printf("pointer1[2] = %f\n",pointer1[2]);
    printf("&pointer1[2] = %p\n",&pointer1[2]);
}


/*
CATS* global_cat;
//*x = kstar
//pars[]
double fit_phase_shifts(double* x, double* pars){
    global_cat->SetShortRangePotential(0,0,0,pars[0]);
    global_cat->SetShortRangePotential(0,0,1,pars[1]);
    global_cat->SetShortRangePotential(0,0,2,pars[2]);
    global_cat->SetShortRangePotential(0,0,3,pars[3]);
    global_cat->KillTheCat();
    return global_cat->EvalPhaseShift(*x,0,0);
}

void my_analysis_function(){

    CATS cat_av18;
    //set up a CATS with AV18
    cat_av18.KillTheCat();
    TH1F* hAV8_PhaseShifts = new TH1F(//blah blah, use same binning as the CATS object//
    );
    //iterate over all bins
    for(){
        hAV8_PhaseShifts->SetBinContent(uBin+1, cat_av18.GetPhaseShift(uBin,0,0));
        //choose some small error, otherwise the fit may fail
        hAV8_PhaseShifts->SetBinError(uBin+1, eps);
    }

    CATS cat;
    //set up your cat here 
    // ....
    //

    //set the global pointer to point to the CATS object you want to use for the fit
    global_cat = &cat;
    //silent the output to avoid flooding the screen while fitting. Warnings/errors will still be shown.
    global_cat->SetNotifications(CATS::nWarning);
    //last argument is the number of parameters, e.g. 4 for a double-Gauss potential
    TF1* custom_fit = new TF1("custom_fit",fit_phase_shifts,0,240,4);
    //set parameters etc., and fit hAV8_PhaseShifts with custom_fit
}
*/



void distroy_the_pc(){
    for(int i=0; i<1000000; i++){
        TH1F* histo;// = new float[8];
        histo = new TH1F("histo","histo",4000,0,100);
        histo->SetBinContent(1,23);

        TH1F histo2("histo2","histo2",100,0,10);

        //printf("1el = %f\n",histo->GetBinContent(1));
        delete histo;
        //printf("1el = %f\n",histo->GetBinContent(1));
    }
    usleep(10e6);
}

void test_root_1(){
    TH1F* histo1 = new TH1F("histo1","histo1",100,0,10);
    TH1F histo2("histo2","histo2",100,0,10);

    

    delete histo1;
}



void plot_WFs(){

    std::vector<float> WF_at = {18, 118};
    const double rmin = 0;
    const double rmax = 8;
    const double source_size = 1.249;

    TString InputFileNameAV18 = TString::Format("%s/Marti/Files_2025-05-12/AV18_histo_1.25.root",GetCernBoxDimi());
    TString InputFileNameSG = TString::Format("%s/Marti/Files_2025-05-12/Gauss_to_AV18_histo_1.25.root",GetCernBoxDimi());
    TString InputFileNameDG = TString::Format("%s/Marti/Files_2025-05-12/Double_Gauss_to_AV18_histo_1.25.root",GetCernBoxDimi());

    TFile fAV18(InputFileNameAV18, "read");
    TH1F* h_Pot_AV18 = (TH1F*)fAV18.Get("spotential");
    TH2F* h_WF2_AV18 = (TH2F*)fAV18.Get("av18_psi2");
    TH2F* h_RWF_AV18 = (TH2F*)fAV18.Get("av18_radial_psi");

    TFile fSG(InputFileNameSG, "read");
    TH1F* h_Pot_SG = (TH1F*)fSG.Get("gauss_potential_h");
    TH2F* h_WF2_SG = (TH2F*)fSG.Get("gauss_to_AV18_psi2");
    TH2F* h_RWF_SG = (TH2F*)fSG.Get("gauss_to_av18_radial_psi");

    TFile fDG(InputFileNameDG, "read");
    TH1F* h_Pot_DG = (TH1F*)fDG.Get("gauss_potential_h");
    TH2F* h_WF2_DG = (TH2F*)fDG.Get("double_gauss_to_AV18_psi2");
    TH2F* h_RWF_DG = (TH2F*)fDG.Get("double_gauss_to_av18_radial_psi");  

    //TH1F* hAxisWF2 = new TH1F("hAxisWF2", "hAxisWF2", 128, 30, 300);
    TF1* fPP_Gauss = new TF1(TString::Format("fPP_Gauss"),GaussSourceScaledTF1,rmin,rmax,2);
    fPP_Gauss->FixParameter(0, source_size);
    fPP_Gauss->FixParameter(1, 1);
    fPP_Gauss->SetLineColor(kGray);
    fPP_Gauss->SetLineWidth(2);
    double max_src = fPP_Gauss->GetMaximum(rmin,rmax);



    for(int iWF=0; iWF<WF_at.size(); iWF++){
        float kstar = WF_at.at(iWF);
        int kstar_bin = h_WF2_AV18->GetYaxis()->FindBin(kstar);

        TH1F* hWF2_AV18 = (TH1F*)h_WF2_AV18->ProjectionX(TString::Format("hWF2_AV18_%.0f",kstar),kstar_bin,kstar_bin);
        TH1F* hRWF_AV18 = (TH1F*)h_RWF_AV18->ProjectionX(TString::Format("hRWF_AV18_%.0f",kstar),kstar_bin,kstar_bin);
        hWF2_AV18->SetLineColor(kBlack);
        hWF2_AV18->SetLineWidth(8);
        hRWF_AV18->SetLineColor(kBlack);
        hRWF_AV18->SetLineWidth(8);

        hWF2_AV18->SetStats(false);
        hWF2_AV18->SetTitle("");
        hWF2_AV18->GetXaxis()->SetTitle("r* (fm)");
        hWF2_AV18->GetXaxis()->SetTitleSize(0.06);
        hWF2_AV18->GetXaxis()->SetLabelSize(0.06);
        hWF2_AV18->GetXaxis()->CenterTitle();
        hWF2_AV18->GetXaxis()->SetTitleOffset(1.3);
        hWF2_AV18->GetXaxis()->SetLabelOffset(0.02);

        hWF2_AV18->GetYaxis()->SetTitle("|#Psi_{k*}(r*)|^{2}");
        hWF2_AV18->GetYaxis()->SetTitleSize(0.06);
        hWF2_AV18->GetYaxis()->SetLabelSize(0.06);
        hWF2_AV18->GetYaxis()->CenterTitle();
        hWF2_AV18->GetYaxis()->SetTitleOffset(1.10);        

        hWF2_AV18->GetXaxis()->SetRangeUser(rmin, rmax);


        hRWF_AV18->SetStats(false);
        hRWF_AV18->SetTitle("");
        hRWF_AV18->GetXaxis()->SetTitle("r* (fm)");
        hRWF_AV18->GetXaxis()->SetTitleSize(0.06);
        hRWF_AV18->GetXaxis()->SetLabelSize(0.06);
        hRWF_AV18->GetXaxis()->CenterTitle();
        hRWF_AV18->GetXaxis()->SetTitleOffset(1.3);
        hRWF_AV18->GetXaxis()->SetLabelOffset(0.02);

        hRWF_AV18->GetYaxis()->SetTitle("R_{k*}(r*)");
        hRWF_AV18->GetYaxis()->SetTitleSize(0.06);
        hRWF_AV18->GetYaxis()->SetLabelSize(0.06);
        hRWF_AV18->GetYaxis()->CenterTitle();
        hRWF_AV18->GetYaxis()->SetTitleOffset(1.10); 

        hRWF_AV18->GetXaxis()->SetRangeUser(rmin, rmax);


        TH1F* hWF2_SG = (TH1F*)h_WF2_SG->ProjectionX(TString::Format("hWF2_SG_%.0f",kstar),kstar_bin,kstar_bin);
        TH1F* hRWF_SG = (TH1F*)h_RWF_SG->ProjectionX(TString::Format("hRWF_SG_%.0f",kstar),kstar_bin,kstar_bin);
        hWF2_SG->SetLineColor(kRed);
        hWF2_SG->SetLineWidth(6);
        hRWF_SG->SetLineColor(kRed);
        hRWF_SG->SetLineWidth(6);

        TH1F* hWF2_DG = (TH1F*)h_WF2_DG->ProjectionX(TString::Format("hWF2_DG_%.0f",kstar),kstar_bin,kstar_bin);
        TH1F* hRWF_DG = (TH1F*)h_RWF_DG->ProjectionX(TString::Format("hRWF_DG_%.0f",kstar),kstar_bin,kstar_bin);
        hWF2_DG->SetLineColor(kGreen+1);
        hWF2_DG->SetLineWidth(4.5);
        hRWF_DG->SetLineColor(kGreen+1);
        hRWF_DG->SetLineWidth(4.5);


        TLegend* lWFs = new TLegend(0.55,0.75,0.95,0.95);//lbrt
        lWFs->SetName("lWFs");
        lWFs->SetTextSize(0.045);
        lWFs->AddEntry(hWF2_AV18, "AV18");
        lWFs->AddEntry(hWF2_SG, "Gauss");
        lWFs->AddEntry(hWF2_DG, "Double Gauss");
        lWFs->AddEntry(fPP_Gauss, "Scaled Gauss (r_{eff}=1.249 fm)");

        TLegend* lRWFs = new TLegend(0.65,0.75,0.95,0.95);//lbrt
        lRWFs->SetName("lRWFs");
        lRWFs->SetTextSize(0.045);
        lRWFs->AddEntry(hWF2_AV18, "AV18");
        lRWFs->AddEntry(hWF2_SG, "Gauss");
        lRWFs->AddEntry(hWF2_DG, "Double Gauss");

        double ymin,ymax;
        ymin = hWF2_AV18->GetMinimum();
        ymax = hWF2_AV18->GetMaximum();
        if(hWF2_SG->GetMinimum()<ymin) ymin=hWF2_SG->GetMinimum();
        if(hWF2_DG->GetMinimum()<ymin) ymin=hWF2_DG->GetMinimum();
        if(hWF2_SG->GetMaximum()>ymax) ymax=hWF2_SG->GetMaximum();
        if(hWF2_DG->GetMaximum()>ymax) ymax=hWF2_DG->GetMaximum();
        hWF2_AV18->GetYaxis()->SetRangeUser(ymin-(ymax-ymin)*0.05, ymax+(ymax-ymin)*0.05);

        //printf("%f %f\n",ymax, max_src);
        fPP_Gauss->FixParameter(1, ymax/max_src*0.9);

        ymin = hRWF_AV18->GetMinimum();
        ymax = hRWF_AV18->GetMaximum();
        if(hRWF_SG->GetMinimum()<ymin) ymin=hRWF_SG->GetMinimum();
        if(hRWF_DG->GetMinimum()<ymin) ymin=hRWF_DG->GetMinimum();
        if(hRWF_SG->GetMaximum()>ymax) ymax=hRWF_SG->GetMaximum();
        if(hRWF_DG->GetMaximum()>ymax) ymax=hRWF_DG->GetMaximum();
        hRWF_AV18->GetYaxis()->SetRangeUser(ymin-(ymax-ymin)*0.05, ymax+(ymax-ymin)*0.05);

 

        TCanvas* c_WF2 = new TCanvas("c_WF2", "c_WF2", 1);
        c_WF2->cd(0);
        c_WF2->SetCanvasSize(1280, 720);
        c_WF2->SetMargin(0.15,0.05,0.2,0.05);//lrbt

        hWF2_AV18->Draw("axis");
        fPP_Gauss->Draw("same");
        hWF2_AV18->Draw("same");
        hWF2_SG->Draw("same");
        hWF2_DG->Draw("same");
        lWFs->Draw("same");

        c_WF2->SaveAs(TString::Format("%s/Marti/Files_2025-05-12/c_WF2_%.0f.png",GetCernBoxDimi(),kstar));

        TCanvas* c_RWF = new TCanvas("c_RWF", "c_RWF", 1);
        c_RWF->cd(0);
        c_RWF->SetCanvasSize(1280, 720);
        c_RWF->SetMargin(0.15,0.05,0.2,0.05);//lrbt

        //hRWF_AV18->Draw("axis");
        hRWF_AV18->Draw();
        hRWF_SG->Draw("same");
        hRWF_DG->Draw("same");
        lRWFs->Draw("same");

        c_RWF->SaveAs(TString::Format("%s/Marti/Files_2025-05-12/c_RWF_%.0f.png",GetCernBoxDimi(),kstar));

        delete hWF2_AV18;
        delete hRWF_AV18;

        delete hWF2_SG;
        delete hRWF_SG;

        delete hWF2_DG;
        delete hRWF_DG;   
        
        //delete fPP_Gauss;

        delete lWFs;
        delete lRWFs;
        delete c_WF2;
        delete c_RWF;
    }

    




    //hAxisWF2->GetYaxis()->SetRangeUser(-0.1,WFtot_yMax);


}



  // Pars[0] Radius in fm
  // Pars[1] Momentum in MeV/c
  // pars[2] = num spline knots
  // pars[3] = der at 0
  // pars[4] = der at large r, good idea to be zero if we want a pot that goes to zero
  // pars[5+iKnot] = x positions
  // pars[5+NumKnots+iKnot] = y positions
double SplinePotentialCats(double* pars){
    //double& rad = pars[0];
    //double& mom = pars[1];
    //int num_knots = round(pars[2]);
    return DlmTSplineFit(pars,&pars[3]);
}
  // pars[0] = num spline knots
  // pars[1] = der at 0
  // pars[2] = der at large r, good idea to be zero if we want a pot that goes to zero
  // pars[3+iKnot] = x positions
  // pars[3+NumKnots+iKnot] = y positions
//double SplinePotential(double*x, double* pars){
//    return DlmTSplineFit(x,pars);
//}

//we try to fit the av18 with a spline potential, where we do two things:
//1: fit the potential itself
//2: fit the phase shifts
void fit_av18_withSpline(){

    unsigned NumKnots = 17;
    //std::vector<float> x_knots = {0,0.3,0.6,1.0,1.3,1.6,3.2,6.4};
    std::vector<float> x_knots;
    for(int iKnot=0; iKnot<NumKnots; iKnot++){
        x_knots.push_back(iKnot*3.2/float(NumKnots-1));
    }

    unsigned kSteps = 50;
    double kMin = 0;
    double kMax = 200;
    double source_size = 1.2;
    CATS smelly_cat_av18;
	smelly_cat_av18.SetMomBins(kSteps, kMin, kMax);

    DLM_CommonAnaFunctions AnalysisObject;
    AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles/",GetCernBoxDimi()));
    AnalysisObject.SetUpCats_pp(smelly_cat_av18,"AV18_s","Gauss",0, 0);
    smelly_cat_av18.SetAnaSource(0, source_size);
    smelly_cat_av18.SetQ1Q2(0);
    smelly_cat_av18.KillTheCat();

    unsigned rSteps = 512;
    double rMin = 0;
    double rMax = x_knots.at(NumKnots-1)+x_knots.at(1)*0.5;
    TH1F* hAV18_V = new TH1F("hAV18_V","hAV18_V",rSteps,rMin,rMax);
    for(unsigned uRad=0; uRad<rSteps; uRad++){
        double rad = hAV18_V->GetBinCenter(uRad+1);
        hAV18_V->SetBinContent(uRad+1, smelly_cat_av18.EvaluateThePotential(0,0,1,rad));
        hAV18_V->SetBinError(uRad+1, 1e-2+smelly_cat_av18.EvaluateThePotential(0,0,1,rad)*1e-3);
    }
    TH1F* hAV18_delta = new TH1F("hAV18_delta","hAV18_delta",kSteps,kMin,kMax);
    for(unsigned uMom=0; uMom<kSteps; uMom++){
        double MOM = hAV18_delta->GetBinCenter(uMom+1);
        hAV18_delta->SetBinContent(uMom+1, smelly_cat_av18.GetPhaseShift(uMom,0,0));
        hAV18_delta->SetBinError(uMom+1, 1e-2+smelly_cat_av18.GetPhaseShift(uMom,0,0)*1e-2);
    }

    TF1* fPot = new TF1("fPot",DlmTSplineFit,rMin,rMax,NumKnots*2+3);
    fPot->FixParameter(0, NumKnots);
    fPot->SetParameter(1, 0);
    fPot->FixParameter(2, 0);
    for(int iKnot=0; iKnot<NumKnots; iKnot++){
        fPot->FixParameter(3+iKnot, x_knots.at(iKnot));
        //fPot->SetParameter(3+NumKnots+iKnot, x_knots.at(iKnot));
        //fPot->SetParLimits(3+NumKnots+iKnot, -1e4, 1e4);

        //double rad = hAV18_V->GetBinCenter(hAV18_V->FindBin(x_knots.at(iKnot)));
        fPot->SetParameter(3+NumKnots+iKnot, hAV18_V->GetBinContent(hAV18_V->FindBin(x_knots.at(iKnot))));
        printf("%i %f %i %f\n",
                    iKnot,
                    x_knots.at(iKnot),
                    hAV18_V->FindBin(x_knots.at(iKnot)),
                    hAV18_V->GetBinContent(hAV18_V->FindBin(x_knots.at(iKnot))));
    }

    //fPot->FixParameter(1, 0);

    //hAV18_V->Fit(fPot, "S, N, R, M");
    fPot->SetNpx(2048);

    CATS smelly_cat_spline_V;
	smelly_cat_spline_V.SetMomBins(kSteps, kMin, kMax);    

    TFile fOutput(TString::Format("%s/Marty/fit_av18_withSpline.root",GetFemtoOutputFolder()), "recreate");
    hAV18_V->Write();
    hAV18_delta->Write();
    fPot->Write();
}


DLM_Histo<float>* dlm_pot = NULL;
//pars[0] = rad
//pars[1] = kstar
double DLM_Potential(double* pars){
    //if(!dlm_pot) return 0;
    //if(dlm_pot->GetTotBin(pars))
    //return dlm_pot->Eval(pars);
}

//the idea is to craft a potential that is defined as a DLM_Histo and interpolates between 
//several points
void DLM_Potential_main1(TString OutputFolder, const unsigned NumRadPts, const unsigned NumIter, const int SEED){
    const double kMin = 0;
    const double kMax = 200;
    const unsigned kSteps = 50;

    const double rMin = 0;
    const double rMax = 8;
    //we will have 2 points at rMin, rMax, NumRadPts are the points inbetween
    //the potential function will return zero if we try to go outside rMin or rMax
    const unsigned rSteps = NumRadPts+2;
    //DLM_Histo<float> my_potential;
    //my_potential.SetUp(2);
    //my_potential.SetUp(0,);
}



int MARTY_MAIN(int argc, char *argv[]){
    printf("Hello Marty\n");
    //plot_WFs();
    //test_mem();
    //test_root_1();
    //distroy_the_pc();
    fit_av18_withSpline();
    return 0;
}