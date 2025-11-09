
#include "OtherTasks.h"
#include "gsl_sf_coulomb.h"
#include <bitset>

#include "CATS.h"
#include "CATSconstants.h"
#include "CommonAnaFunctions.h"
#include "DLM_Potentials.h"
#include "DLM_Random.h"
#include "DLM_CkDecomposition.h"
#include "DLM_Source.h"
#include "DLM_Potentials.h"
#include "DLM_WfModel.h"
#include "DLM_CkModels.h"
#include "DLM_Ck.h"
#include "DLM_RootWrapper.h"
#include "DLM_HistoAnalysis.h"
#include "EnvVars.h"
#include "DLM_Fitters.h"
#include "DLM_RootFit.h"
#include "DLM_Unfold.h"
#include "DLM_Sort.h"
#include "DLM_DrawingTools.h"
#include "DLM_DecayMatrix.h"
#include "DLM_MathFunctions.h"

#include <TApplication.h>
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
//#include "TGraphMultiErrors.h"
#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH3D.h"
#include "TUnfoldDensity.h"
#include "TNtuple.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TFractionFitter.h"
#include "TGenPhaseSpace.h"
#include "TString.h"
#include "TList.h"
#include "TROOT.h"
#include "Math/MinimizerOptions.h"
//./math/mathcore/inc/Math/MinimizerOptions.h
#include "gsl_sf_dawson.h"
#include "Math/Vector4Dfwd.h"

void pp_CompareToNorfolk(){
    const TString OutputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/OtherTasks/pp_CompareToNorfolk/Reid/FAST/";
    const unsigned NumPotentials = 4;
    CATS* Kitty = new CATS [NumPotentials];
    TString* PotNames = new TString [NumPotentials];
    PotNames[0] = "AV18";
    PotNames[1] = "ReidV8";
    PotNames[2] = "ReidSC";
    PotNames[3] = "Norfolk";
    //PotNames[4] = "Norfolk";
    //PotNames[5] = "Norfolk";
    //PotNames[6] = "Norfolk";
    //PotNames[7] = "Norfolk";
    //PotNames[8] = "Norfolk";

    int* PotVars = new int [NumPotentials];
    PotVars[0] = 0;
    PotVars[1] = 0;
    PotVars[2] = 0;
    PotVars[3] = 1090;
    //PotVars[3] = 1050;
    //PotVars[4] = 1060;
    //PotVars[5] = 1090;
    //PotVars[6] = 1091;
    //PotVars[7] = 1092;
    //PotVars[6] = 1100;

    TGraph* gPot_S = new TGraph[NumPotentials];
    TCanvas* cPot_S = new TCanvas[NumPotentials];
    TGraph* gPot_D = new TGraph[NumPotentials];
    TCanvas* cPot_D = new TCanvas[NumPotentials];
    TGraph** gPot_P = new TGraph*[NumPotentials];
    TCanvas** cPot_P = new TCanvas*[NumPotentials];
    TGraph* gCk = new TGraph[NumPotentials];
    TCanvas* cCk = new TCanvas[NumPotentials];
    TGraph** gRatio = new TGraph*[NumPotentials];
    TCanvas** cRatio = new TCanvas*[NumPotentials];
    for(unsigned uPot=0; uPot<NumPotentials; uPot++) gPot_P[uPot] = new TGraph[3];
    for(unsigned uPot=0; uPot<NumPotentials; uPot++){
        cPot_P[uPot] = new TCanvas[NumPotentials];
        cRatio[uPot] = new TCanvas[NumPotentials];
    }
    const unsigned NumMomBins = 80;
    const double kMin = 0;
    const double kMax = 320;

    DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder("/home/dmihaylov/CernBox/CatsFiles");

    for(unsigned uPot=0; uPot<NumPotentials; uPot++){
        Kitty[uPot].SetMomBins(NumMomBins,kMin,kMax);
//printf("SETUP THE CAT\n");
        AnalysisObject.SetUpCats_pp(Kitty[uPot],PotNames[uPot],"Gauss",PotVars[uPot]);
        Kitty[uPot].SetAnaSource(0,0.6);
        Kitty[uPot].RemoveShortRangePotential(0,2);
        Kitty[uPot].RemoveShortRangePotential(1,1);
        Kitty[uPot].RemoveShortRangePotential(2,1);
        Kitty[uPot].RemoveShortRangePotential(3,1);

//printf("1S0: %f\n",Kitty[uPot].EvaluateThePotential(0,0,10,1.0));
//printf("3P0: %f\n",Kitty[uPot].EvaluateThePotential(1,1,10,1.0));
//double PARS[5];
//PARS[0]=1.0; PARS[2]=105; PARS[3]=-1; PARS[4]=310;
//printf("pp_Norfolk: %f\n",pp_Norfolk(PARS));
//printf("KILLING THE CAT\n");
        Kitty[uPot].SetEpsilonConv(4e-8);
        Kitty[uPot].SetEpsilonProp(4e-8);
        Kitty[uPot].KillTheCat();
    }

    const unsigned NumRadPts = 1024;
    const double rMin = 0;
    const double rMax = 16;
    const double rStep = (rMax-rMin)/double(NumRadPts);

    double RadVal;
    double PotVal;

    for(unsigned uPot=0; uPot<NumPotentials; uPot++){
        gPot_S[uPot].SetName("V_1S0_"+PotNames[uPot]+"_"+PotVars[uPot]);
        gPot_D[uPot].SetName("V_1D2_"+PotNames[uPot]+"_"+PotVars[uPot]);
        gPot_P[uPot][0].SetName("V_3P0_"+PotNames[uPot]+"_"+PotVars[uPot]);
        gPot_P[uPot][1].SetName("V_3P1_"+PotNames[uPot]+"_"+PotVars[uPot]);
        gPot_P[uPot][2].SetName("V_3P2_"+PotNames[uPot]+"_"+PotVars[uPot]);
        for(unsigned uRad=0; uRad<NumRadPts; uRad++){
            RadVal = rMin + rStep*0.5 + rStep*double(uRad);

            gPot_S[uPot].SetPoint(uRad,RadVal,Kitty[uPot].EvaluateThePotential(0,0,10,RadVal));
            gPot_D[uPot].SetPoint(uRad,RadVal,Kitty[uPot].EvaluateThePotential(0,2,10,RadVal));
            gPot_P[uPot][0].SetPoint(uRad,RadVal,Kitty[uPot].EvaluateThePotential(1,1,10,RadVal));
            gPot_P[uPot][1].SetPoint(uRad,RadVal,Kitty[uPot].EvaluateThePotential(2,1,10,RadVal));
            gPot_P[uPot][2].SetPoint(uRad,RadVal,Kitty[uPot].EvaluateThePotential(3,1,10,RadVal));
        }

        gCk[uPot].SetName("Ck_"+PotNames[uPot]+"_"+PotVars[uPot]);
        for(unsigned uMom=0; uMom<NumMomBins; uMom++){
            gCk[uPot].SetPoint(uMom,Kitty[uPot].GetMomentum(uMom),Kitty[uPot].GetCorrFun(uMom));
        }
    }

    for(unsigned uPot=0; uPot<NumPotentials; uPot++){
        gRatio[uPot] = new TGraph [NumPotentials];
        for(unsigned uPot2=0; uPot2<NumPotentials; uPot2++){
            gRatio[uPot][uPot2].SetName("CkRatio_"+PotNames[uPot]+"_"+PotVars[uPot]+"_vs_"+PotNames[uPot2]+"_"+PotVars[uPot2]);
            double CkPOT, CkPOT2, Mom;
            for(unsigned uMom=0; uMom<NumMomBins; uMom++){
                gCk[uPot].GetPoint(uMom,Mom,CkPOT);
                gCk[uPot2].GetPoint(uMom,Mom,CkPOT2);
                gRatio[uPot][uPot2].SetPoint(uMom,Mom,CkPOT/CkPOT2);
            }
        }
    }

    TFile fOutput(OutputFolder+"fOutput.root","recreate");
    for(unsigned uPot=0; uPot<NumPotentials; uPot++){gPot_S[uPot].Write();}
    for(unsigned uPot=0; uPot<NumPotentials; uPot++){gPot_P[uPot][0].Write();}
    for(unsigned uPot=0; uPot<NumPotentials; uPot++){gPot_P[uPot][1].Write();}
    for(unsigned uPot=0; uPot<NumPotentials; uPot++){gPot_P[uPot][2].Write();}
    for(unsigned uPot=0; uPot<NumPotentials; uPot++){gPot_D[uPot].Write();}
    for(unsigned uPot=0; uPot<NumPotentials; uPot++){gCk[uPot].Write();}
    for(unsigned uPot=0; uPot<NumPotentials; uPot++)for(unsigned uPot2=0; uPot2<NumPotentials; uPot2++)
    {if(uPot!=uPot2)gRatio[uPot][uPot2].Write();}

    //SLIDES FOR LAURA
    TCanvas* cPotLaura = new TCanvas("cPotLaura", "cPotLaura", 1);
    cPotLaura->cd(0);
    cPotLaura->SetCanvasSize(1280, 720);
    cPotLaura->SetMargin(0.15,0.05,0.2,0.05);//lrbt

    TH1F* hAxisPot = new TH1F("hAxisPot", "hAxisPot", 128, 0, 1.5);
    hAxisPot->SetStats(false);
    hAxisPot->SetTitle("");
    hAxisPot->GetXaxis()->SetTitle("r (fm)");
    hAxisPot->GetXaxis()->SetTitleSize(0.06);
    hAxisPot->GetXaxis()->SetLabelSize(0.06);
    hAxisPot->GetXaxis()->CenterTitle();
    hAxisPot->GetXaxis()->SetTitleOffset(1.3);
    hAxisPot->GetXaxis()->SetLabelOffset(0.02);

    hAxisPot->GetYaxis()->SetTitle("V(r) (MeV)");
    hAxisPot->GetYaxis()->SetTitleSize(0.06);
    hAxisPot->GetYaxis()->SetLabelSize(0.06);
    hAxisPot->GetYaxis()->CenterTitle();
    hAxisPot->GetYaxis()->SetTitleOffset(1.10);

    //hAxisPot->GetXaxis()->SetNdivisions(506);
    hAxisPot->GetYaxis()->SetRangeUser(-200, 3200);

    gPot_S[0].SetLineWidth(4);
    gPot_S[1].SetLineWidth(4);
    gPot_S[3].SetLineWidth(4);

    gPot_S[0].SetLineColor(kBlue);
    gPot_S[1].SetLineColor(kRed);
    gPot_S[3].SetLineColor(kGreen+1);

    gPot_S[0].SetLineStyle(1);
    gPot_S[1].SetLineStyle(2);
    gPot_S[3].SetLineStyle(5);

    TLegend* lPotLaura = new TLegend(0.65,0.75,0.95,0.95);//lbrt
    lPotLaura->SetName("lPotLaura");
    lPotLaura->SetTextSize(0.045);
    lPotLaura->AddEntry(&gPot_S[0], "AV18 (1S0)");
    lPotLaura->AddEntry(&gPot_S[1], "Reid (1S0)");
    lPotLaura->AddEntry(&gPot_S[3], "Norfolk (1S0)");// R=(1.0,0.8) fm E_{lab}<200 MeV
    hAxisPot->Draw("axis");
    gPot_S[0].Draw("C");
    gPot_S[1].Draw("C,same");
    gPot_S[3].Draw("C,same");
    lPotLaura->Draw("same");

    cPotLaura->SaveAs(OutputFolder+"cPotLaura.png");


    TCanvas* cCkLaura = new TCanvas("cCkLaura", "cCkLaura", 1);
    cCkLaura->cd(0);
    cCkLaura->SetCanvasSize(1280, 720);
    cCkLaura->SetMargin(0.15,0.05,0.2,0.05);//lrbt

    TH1F* hAxisCk = new TH1F("hAxisCk", "hAxisCk", 128, 0, 300);
    hAxisCk->SetStats(false);
    hAxisCk->SetTitle("");
    hAxisCk->GetXaxis()->SetTitle("k* (MeV)");
    hAxisCk->GetXaxis()->SetTitleSize(0.06);
    hAxisCk->GetXaxis()->SetLabelSize(0.06);
    hAxisCk->GetXaxis()->CenterTitle();
    hAxisCk->GetXaxis()->SetTitleOffset(1.3);
    hAxisCk->GetXaxis()->SetLabelOffset(0.02);

    hAxisCk->GetYaxis()->SetTitle("C(k*)");
    hAxisCk->GetYaxis()->SetTitleSize(0.06);
    hAxisCk->GetYaxis()->SetLabelSize(0.06);
    hAxisCk->GetYaxis()->CenterTitle();
    hAxisCk->GetYaxis()->SetTitleOffset(1.10);

    //hAxisCk->GetXaxis()->SetNdivisions(506);
    hAxisCk->GetYaxis()->SetRangeUser(0, 5);

    gCk[0].SetLineWidth(4);
    gCk[1].SetLineWidth(4);
    gCk[3].SetLineWidth(4);

    gCk[0].SetLineColor(kBlue);
    gCk[1].SetLineColor(kRed);
    gCk[3].SetLineColor(kGreen+1);

    gCk[0].SetLineStyle(1);
    gCk[1].SetLineStyle(2);
    gCk[3].SetLineStyle(5);

    TLegend* lCkLaura = new TLegend(0.65,0.75,0.95,0.95);//lbrt
    lCkLaura->SetName("lCkLaura");
    lCkLaura->SetTextSize(0.045);
    //lCkLaura->AddEntry(&gCk[0], "AV18 (1S0)");
    //lCkLaura->AddEntry(&gCk[1], "Reid (1S0)");
    //lCkLaura->AddEntry(&gCk[3], "Norfolk (1S0)");// R=(1.0,0.8) fm E_{lab}<200 MeV
    lCkLaura->AddEntry(&gCk[0], "AV18 (1S0)");//s+p waves
    lCkLaura->AddEntry(&gCk[1], "Reid (1S0)");
    lCkLaura->AddEntry(&gCk[3], "Norfolk (1S0)");// R=(1.0,0.8) fm E_{lab}<200 MeV
    hAxisCk->Draw("axis");
    gCk[0].Draw("C");
    gCk[1].Draw("C,same");
    gCk[3].Draw("C,same");
    lCkLaura->Draw("same");

    cCkLaura->SaveAs(OutputFolder+"cCkLaura.png");


    TCanvas* cCkZoomLaura = new TCanvas("cCkZoomLaura", "cCkZoomLaura", 1);
    cCkZoomLaura->cd(0);
    cCkZoomLaura->SetCanvasSize(1280, 720);
    cCkZoomLaura->SetMargin(0.15,0.05,0.2,0.05);//lrbt

    TH1F* hAxisCkZoom = new TH1F("hAxisCkZoom", "hAxisCkZoom", 128, 30, 300);
    hAxisCkZoom->SetStats(false);
    hAxisCkZoom->SetTitle("");
    hAxisCkZoom->GetXaxis()->SetTitle("k* (MeV)");
    hAxisCkZoom->GetXaxis()->SetTitleSize(0.06);
    hAxisCkZoom->GetXaxis()->SetLabelSize(0.06);
    hAxisCkZoom->GetXaxis()->CenterTitle();
    hAxisCkZoom->GetXaxis()->SetTitleOffset(1.3);
    hAxisCkZoom->GetXaxis()->SetLabelOffset(0.02);

    hAxisCkZoom->GetYaxis()->SetTitle("C(k*)");
    hAxisCkZoom->GetYaxis()->SetTitleSize(0.06);
    hAxisCkZoom->GetYaxis()->SetLabelSize(0.06);
    hAxisCkZoom->GetYaxis()->CenterTitle();
    hAxisCkZoom->GetYaxis()->SetTitleOffset(1.10);

    //hAxisCk->GetXaxis()->SetNdivisions(506);
    hAxisCkZoom->GetYaxis()->SetRangeUser(0.87, 1.05);

    hAxisCkZoom->Draw("axis");
    gCk[0].Draw("C");
    gCk[1].Draw("C,same");
    gCk[3].Draw("C,same");
    lCkLaura->Draw("same");

    cCkZoomLaura->SaveAs(OutputFolder+"cCkZoomLaura.png");

    TCanvas* cRatioLaura = new TCanvas("cRatioLaura", "cRatioLaura", 1);
    cRatioLaura->cd(0);
    cRatioLaura->SetCanvasSize(1280, 720);
    cRatioLaura->SetMargin(0.15,0.05,0.2,0.05);//lrbt

    TH1F* hAxisRatio = new TH1F("hAxisRatio", "hAxisRatio", 128, 0, 300);
    hAxisRatio->SetStats(false);
    hAxisRatio->SetTitle("");
    hAxisRatio->GetXaxis()->SetTitle("k* (MeV)");
    hAxisRatio->GetXaxis()->SetTitleSize(0.06);
    hAxisRatio->GetXaxis()->SetLabelSize(0.06);
    hAxisRatio->GetXaxis()->CenterTitle();
    hAxisRatio->GetXaxis()->SetTitleOffset(1.3);
    hAxisRatio->GetXaxis()->SetLabelOffset(0.02);

    hAxisRatio->GetYaxis()->SetTitle("C(k*) / C(k*)");
    hAxisRatio->GetYaxis()->SetTitleSize(0.06);
    hAxisRatio->GetYaxis()->SetLabelSize(0.06);
    hAxisRatio->GetYaxis()->CenterTitle();
    hAxisRatio->GetYaxis()->SetTitleOffset(1.10);

    //hAxisCk->GetXaxis()->SetNdivisions(506);
    hAxisRatio->GetYaxis()->SetRangeUser(0.984, 1.01);

    TLegend* lRatioLaura = new TLegend(0.65,0.75,0.95,0.95);//lbrt
    lRatioLaura->SetName("lRatioLaura");
    lRatioLaura->SetTextSize(0.045);
    lRatioLaura->AddEntry(&gRatio[0][1], "AV18 / Reid");
    lRatioLaura->AddEntry(&gRatio[0][3], "AV18 / Norfolk");

    gRatio[0][1].SetLineWidth(4);
    gRatio[0][3].SetLineWidth(4);

    gRatio[0][1].SetLineColor(kRed);
    gRatio[0][3].SetLineColor(kGreen+1);

    gRatio[0][1].SetLineStyle(2);
    gRatio[0][3].SetLineStyle(5);

    hAxisRatio->Draw("axis");
    gRatio[0][1].Draw("C");
    gRatio[0][3].Draw("C,same");
    lRatioLaura->Draw("same");

    cRatioLaura->SaveAs(OutputFolder+"cRatioLaura.png");

    delete [] Kitty;
    delete [] PotNames;
    delete [] gCk;
    delete [] gPot_S;
    delete [] gPot_D;
    for(unsigned uPot=0; uPot<NumPotentials; uPot++) delete [] gPot_P[uPot];
    delete [] gPot_P;
    for(unsigned uPot=0; uPot<NumPotentials; uPot++) delete [] gRatio[uPot];
    delete [] gRatio;

    delete [] cCk;
    delete [] cPot_S;
    delete [] cPot_D;
    for(unsigned uPot=0; uPot<NumPotentials; uPot++) delete [] cPot_P[uPot];
    delete [] cPot_P;
    for(unsigned uPot=0; uPot<NumPotentials; uPot++) delete [] cRatio[uPot];
    delete [] cRatio;

    delete [] PotVars;
}


void pp_pL_CorrectedMC_EXP(){
    TString DataFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/Data/CF/Norm024034/";
    TString McFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/MC/CF/AOD_Trains/Norm024034/";
    TString OutputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/OtherTasks/pp_pL_CorrectedMC_EXP/";
    TString DataFile_pp = "CFOutput_pp_8.root";
    TString DataFile_pL = "CFOutput_pL_8.root";
    TString McFile_pp = "CFOutput_pp_8.root";
    TString McFile_pL = "CFOutput_pL_8.root";
    TString DataHisto_pp = "hCk_ReweightedMeV_0";
    TString DataHisto_pL = "hCk_ReweightedMeV_0";
    TString McHisto_pp = "hCk_ReweightedMeV_0";
    TString McHisto_pL = "hCk_ReweightedMeV_0";

    TString DataSeList1_pp = "PairDist";
    TString DataSeList2_pp = "PairReweighted";
    TString DataSe_pp = "SEDist_Particle0_Particle0_clone_Shifted_FixShifted_Reweighted";
    TString DataMe_pp = "MEDist_Particle0_Particle0_clone_Shifted_FixShifted_Reweighted";
    TString DataSeList1_pL = "PairDist";
    TString DataSeList2_pL = "PairReweighted";
    TString DataSe_pL = "SEDist_Particle0_Particle2_clone_Shifted_FixShifted_Rebinned_4_Reweighted";
    TString DataMe_pL = "MEDist_Particle0_Particle2_clone_Shifted_FixShifted_Rebinned_4_Reweighted";
    TString McSeList1_pp = "PairDist";
    TString McSeList2_pp = "PairReweighted";
    TString McSe_pp = "SEDist_Particle0_Particle0_clone_Shifted_FixShifted_Reweighted";
    TString McMe_pp = "MEDist_Particle0_Particle0_clone_Shifted_FixShifted_Reweighted";
    TString McSeList1_pL = "PairDist";
    TString McSeList2_pL = "PairReweighted";
    TString McSe_pL = "SEDist_Particle0_Particle2_clone_Shifted_FixShifted_Rebinned_4_Reweighted";
    TString McMe_pL = "MEDist_Particle0_Particle2_clone_Shifted_FixShifted_Rebinned_4_Reweighted";

    DLM_Histo<float> dlm_DataSe_pp;
    DLM_Histo<float> dlm_DataSe_pL;
    DLM_Histo<float> dlm_McSe_pp;
    DLM_Histo<float> dlm_McSe_pL;

    DLM_Histo<float> dlm_DataMe_pp;
    DLM_Histo<float> dlm_DataMe_pL;
    DLM_Histo<float> dlm_McMe_pp;
    DLM_Histo<float> dlm_McMe_pL;

    DLM_Histo<float> dlm_DataCk_pp;
    DLM_Histo<float> dlm_DataCk_pL;
    DLM_Histo<float> dlm_McCk_pp;
    DLM_Histo<float> dlm_McCk_pL;

    DLM_Histo<float> dlm_RatioCk_ppData_ppMc;
    DLM_Histo<float> dlm_RatioCk_pLData_pLMc;
    DLM_Histo<float> dlm_RatioSe_ppData_ppMc;
    DLM_Histo<float> dlm_RatioSe_pLData_pLMc;
    DLM_Histo<float> dlm_RatioMe_ppData_ppMc;
    DLM_Histo<float> dlm_RatioMe_pLData_pLMc;

    unsigned NumBins;
    double kMin;
    double kMax;

    TH1F* hCk_tmp;
    TList* list1_tmp;
    TList* list2_tmp;
    TH1F* hSe_tmp;
    TH1F* hMe_tmp;

    TFile datafile_pp(DataFolder+DataFile_pp);
    hCk_tmp = (TH1F*)datafile_pp.Get(DataHisto_pp);
    //hCk_tmp->Rebin(4); hCk_tmp->Scale(1./4.);
    list1_tmp = (TList*)datafile_pp.Get(DataSeList1_pp);
    list2_tmp = (TList*)list1_tmp->FindObject(DataSeList2_pp);
    hSe_tmp = (TH1F*)list2_tmp->FindObject(DataSe_pp);
    hSe_tmp->Scale(1./hSe_tmp->Integral(),"width");
    hMe_tmp = (TH1F*)list2_tmp->FindObject(DataMe_pp);
    hMe_tmp->Scale(1./hMe_tmp->Integral(),"width");
    NumBins = hCk_tmp->GetNbinsX();
    kMin = hCk_tmp->GetBinLowEdge(1);
    kMax = hCk_tmp->GetXaxis()->GetBinUpEdge(NumBins);
    dlm_DataCk_pp.SetUp(1);
    dlm_DataCk_pp.SetUp(0,NumBins,kMin,kMax);
    dlm_DataCk_pp.Initialize();
    dlm_DataSe_pp.SetUp(1);
    dlm_DataSe_pp.SetUp(0,NumBins,kMin,kMax);
    dlm_DataSe_pp.Initialize();
    dlm_DataMe_pp.SetUp(1);
    dlm_DataMe_pp.SetUp(0,NumBins,kMin,kMax);
    dlm_DataMe_pp.Initialize();
    for(unsigned uBin=0; uBin<NumBins; uBin++){
        dlm_DataCk_pp.SetBinContent(uBin,hCk_tmp->GetBinContent(uBin+1));
        dlm_DataCk_pp.SetBinError(uBin,hCk_tmp->GetBinError(uBin+1));
        dlm_DataSe_pp.SetBinContent(uBin,hSe_tmp->GetBinContent(uBin+1));
        dlm_DataSe_pp.SetBinError(uBin,hSe_tmp->GetBinError(uBin+1));
        dlm_DataMe_pp.SetBinContent(uBin,hMe_tmp->GetBinContent(uBin+1));
        dlm_DataMe_pp.SetBinError(uBin,hMe_tmp->GetBinError(uBin+1));
    }
    datafile_pp.Close();

    TFile datafile_pL(DataFolder+DataFile_pL);
    hCk_tmp = (TH1F*)datafile_pL.Get(DataHisto_pL);
    list1_tmp = (TList*)datafile_pL.Get(DataSeList1_pL);
    list2_tmp = (TList*)list1_tmp->FindObject(DataSeList2_pL);
    hSe_tmp = (TH1F*)list2_tmp->FindObject(DataSe_pL);
    hSe_tmp->Scale(1./hSe_tmp->Integral(),"width");
    hMe_tmp = (TH1F*)list2_tmp->FindObject(DataMe_pL);
    hMe_tmp->Scale(1./hMe_tmp->Integral(),"width");
    NumBins = hCk_tmp->GetNbinsX();
    kMin = hCk_tmp->GetBinLowEdge(1);
    kMax = hCk_tmp->GetXaxis()->GetBinUpEdge(NumBins);
    dlm_DataCk_pL.SetUp(1);
    dlm_DataCk_pL.SetUp(0,NumBins,kMin,kMax);
    dlm_DataCk_pL.Initialize();
    dlm_DataSe_pL.SetUp(1);
    dlm_DataSe_pL.SetUp(0,NumBins,kMin,kMax);
    dlm_DataSe_pL.Initialize();
    dlm_DataMe_pL.SetUp(1);
    dlm_DataMe_pL.SetUp(0,NumBins,kMin,kMax);
    dlm_DataMe_pL.Initialize();
    for(unsigned uBin=0; uBin<NumBins; uBin++){
        dlm_DataCk_pL.SetBinContent(uBin,hCk_tmp->GetBinContent(uBin+1));
        dlm_DataCk_pL.SetBinError(uBin,hCk_tmp->GetBinError(uBin+1));
        dlm_DataSe_pL.SetBinContent(uBin,hSe_tmp->GetBinContent(uBin+1));
        dlm_DataSe_pL.SetBinError(uBin,hSe_tmp->GetBinError(uBin+1));
        dlm_DataMe_pL.SetBinContent(uBin,hMe_tmp->GetBinContent(uBin+1));
        dlm_DataMe_pL.SetBinError(uBin,hMe_tmp->GetBinError(uBin+1));
    }
    datafile_pL.Close();

    TFile mcfile_pp(McFolder+McFile_pp);
    hCk_tmp = (TH1F*)mcfile_pp.Get(McHisto_pp);
    //hCk_tmp->Rebin(4); hCk_tmp->Scale(1./4.);
    list1_tmp = (TList*)mcfile_pp.Get(McSeList1_pp);
    list2_tmp = (TList*)list1_tmp->FindObject(McSeList2_pp);
    hSe_tmp = (TH1F*)list2_tmp->FindObject(McSe_pp);
    hSe_tmp->Scale(1./hSe_tmp->Integral(),"width");
    hMe_tmp = (TH1F*)list2_tmp->FindObject(McMe_pp);
    hMe_tmp->Scale(1./hMe_tmp->Integral(),"width");
    NumBins = hCk_tmp->GetNbinsX();
    kMin = hCk_tmp->GetBinLowEdge(1);
    kMax = hCk_tmp->GetXaxis()->GetBinUpEdge(NumBins);
    dlm_McCk_pp.SetUp(1);
    dlm_McCk_pp.SetUp(0,NumBins,kMin,kMax);
    dlm_McCk_pp.Initialize();
    dlm_McSe_pp.SetUp(1);
    dlm_McSe_pp.SetUp(0,NumBins,kMin,kMax);
    dlm_McSe_pp.Initialize();
    dlm_McMe_pp.SetUp(1);
    dlm_McMe_pp.SetUp(0,NumBins,kMin,kMax);
    dlm_McMe_pp.Initialize();
    for(unsigned uBin=0; uBin<NumBins; uBin++){
        dlm_McCk_pp.SetBinContent(uBin,hCk_tmp->GetBinContent(uBin+1));
        dlm_McCk_pp.SetBinError(uBin,hCk_tmp->GetBinError(uBin+1));
        dlm_McSe_pp.SetBinContent(uBin,hSe_tmp->GetBinContent(uBin+1));
        dlm_McSe_pp.SetBinError(uBin,hSe_tmp->GetBinError(uBin+1));
        dlm_McMe_pp.SetBinContent(uBin,hMe_tmp->GetBinContent(uBin+1));
        dlm_McMe_pp.SetBinError(uBin,hMe_tmp->GetBinError(uBin+1));
    }
    mcfile_pp.Close();

    TFile mcfile_pL(McFolder+McFile_pL);
    hCk_tmp = (TH1F*)mcfile_pL.Get(McHisto_pL);
    list1_tmp = (TList*)mcfile_pL.Get(McSeList1_pL);
    list2_tmp = (TList*)list1_tmp->FindObject(McSeList2_pL);
    hSe_tmp = (TH1F*)list2_tmp->FindObject(McSe_pL);
    hSe_tmp->Scale(1./hSe_tmp->Integral(),"width");
    hMe_tmp = (TH1F*)list2_tmp->FindObject(McMe_pL);
    hMe_tmp->Scale(1./hMe_tmp->Integral(),"width");
    NumBins = hCk_tmp->GetNbinsX();
    kMin = hCk_tmp->GetBinLowEdge(1);
    kMax = hCk_tmp->GetXaxis()->GetBinUpEdge(NumBins);
    dlm_McCk_pL.SetUp(1);
    dlm_McCk_pL.SetUp(0,NumBins,kMin,kMax);
    dlm_McCk_pL.Initialize();
    dlm_McSe_pL.SetUp(1);
    dlm_McSe_pL.SetUp(0,NumBins,kMin,kMax);
    dlm_McSe_pL.Initialize();
    dlm_McMe_pL.SetUp(1);
    dlm_McMe_pL.SetUp(0,NumBins,kMin,kMax);
    dlm_McMe_pL.Initialize();
    for(unsigned uBin=0; uBin<NumBins; uBin++){
        dlm_McCk_pL.SetBinContent(uBin,hCk_tmp->GetBinContent(uBin+1));
        dlm_McCk_pL.SetBinError(uBin,hCk_tmp->GetBinError(uBin+1));
        dlm_McSe_pL.SetBinContent(uBin,hSe_tmp->GetBinContent(uBin+1));
        dlm_McSe_pL.SetBinError(uBin,hSe_tmp->GetBinError(uBin+1));
        dlm_McMe_pL.SetBinContent(uBin,hMe_tmp->GetBinContent(uBin+1));
        dlm_McMe_pL.SetBinError(uBin,hMe_tmp->GetBinError(uBin+1));
    }
    mcfile_pL.Close();

    //dlm_DataSe_pp.RescaleAxis(0,1000.,true);
    //dlm_DataSe_pL.RescaleAxis(0,1000.,true);
    //dlm_McSe_pp.RescaleAxis(0,1000.,true);
    //dlm_McSe_pL.RescaleAxis(0,1000.,true);
    //dlm_DataMe_pp.RescaleAxis(0,1000.,true);
    //dlm_DataMe_pL.RescaleAxis(0,1000.,true);
    //dlm_McMe_pp.RescaleAxis(0,1000.,true);
    //dlm_McMe_pL.RescaleAxis(0,1000.,true);


    dlm_RatioCk_ppData_ppMc = dlm_DataCk_pp;
    dlm_RatioCk_ppData_ppMc /= dlm_McCk_pp;

    dlm_RatioCk_pLData_pLMc = dlm_DataCk_pL;
    dlm_RatioCk_pLData_pLMc /= dlm_McCk_pL;

    dlm_RatioSe_ppData_ppMc = dlm_DataSe_pp;
printf("b1 = %f\n",dlm_RatioSe_ppData_ppMc.GetBinContent(4));
printf(" e1 = %f\n",dlm_RatioSe_ppData_ppMc.GetBinError(4));
printf("b2 = %f\n",dlm_McSe_pp.GetBinContent(4));
printf(" e2 = %f\n",dlm_McSe_pp.GetBinError(4));
    dlm_RatioSe_ppData_ppMc /= dlm_McSe_pp;

    dlm_RatioSe_pLData_pLMc = dlm_DataSe_pL;
    dlm_RatioSe_pLData_pLMc /= dlm_McSe_pL;

    dlm_RatioMe_ppData_ppMc = dlm_DataMe_pp;
    dlm_RatioMe_ppData_ppMc /= dlm_McMe_pp;

    dlm_RatioMe_pLData_pLMc = dlm_DataMe_pL;
    dlm_RatioMe_pLData_pLMc /= dlm_McMe_pL;

    TFile fOutput(OutputFolder+"fOutput.root","recreate");
    TH1F* histo_DataSe_pp = new TH1F("histo_DataSe_pp","histo_DataSe_pp",dlm_DataSe_pp.GetNbins(),dlm_DataSe_pp.GetLowEdge(0),dlm_DataSe_pp.GetUpEdge(0));
    TH1F* histo_DataSe_pL = new TH1F("histo_DataSe_pL","histo_DataSe_pL",dlm_DataSe_pL.GetNbins(),dlm_DataSe_pL.GetLowEdge(0),dlm_DataSe_pL.GetUpEdge(0));
    TH1F* histo_McSe_pp = new TH1F("histo_McSe_pp","histo_McSe_pp",dlm_McSe_pp.GetNbins(),dlm_McSe_pp.GetLowEdge(0),dlm_McSe_pp.GetUpEdge(0));
    TH1F* histo_McSe_pL = new TH1F("histo_McSe_pL","histo_McSe_pL",dlm_McSe_pL.GetNbins(),dlm_McSe_pL.GetLowEdge(0),dlm_McSe_pL.GetUpEdge(0));

    TH1F* histo_DataMe_pp = new TH1F("histo_DataMe_pp","histo_DataMe_pp",dlm_DataMe_pp.GetNbins(),dlm_DataMe_pp.GetLowEdge(0),dlm_DataMe_pp.GetUpEdge(0));
    TH1F* histo_DataMe_pL = new TH1F("histo_DataMe_pL","histo_DataMe_pL",dlm_DataMe_pL.GetNbins(),dlm_DataMe_pL.GetLowEdge(0),dlm_DataMe_pL.GetUpEdge(0));
    TH1F* histo_McMe_pp = new TH1F("histo_McMe_pp","histo_McMe_pp",dlm_McMe_pp.GetNbins(),dlm_McMe_pp.GetLowEdge(0),dlm_McMe_pp.GetUpEdge(0));
    TH1F* histo_McMe_pL = new TH1F("histo_McMe_pL","histo_McMe_pL",dlm_McMe_pL.GetNbins(),dlm_McMe_pL.GetLowEdge(0),dlm_McMe_pL.GetUpEdge(0));

    TH1F* histo_DataCk_pp = new TH1F("histo_DataCk_pp","histo_DataCk_pp",dlm_DataCk_pp.GetNbins(),dlm_DataCk_pp.GetLowEdge(0),dlm_DataCk_pp.GetUpEdge(0));
    TH1F* histo_DataCk_pL = new TH1F("histo_DataCk_pL","histo_DataCk_pL",dlm_DataCk_pL.GetNbins(),dlm_DataCk_pL.GetLowEdge(0),dlm_DataCk_pL.GetUpEdge(0));
    TH1F* histo_McCk_pp = new TH1F("histo_McCk_pp","histo_McCk_pp",dlm_McCk_pp.GetNbins(),dlm_McCk_pp.GetLowEdge(0),dlm_McCk_pp.GetUpEdge(0));
    TH1F* histo_McCk_pL = new TH1F("histo_McCk_pL","histo_McCk_pL",dlm_McCk_pL.GetNbins(),dlm_McCk_pL.GetLowEdge(0),dlm_McCk_pL.GetUpEdge(0));

    TH1F* histo_RatioCk_ppData_ppMc = new TH1F("histo_RatioCk_ppData_ppMc","histo_RatioCk_ppData_ppMc",dlm_RatioCk_ppData_ppMc.GetNbins(),dlm_RatioCk_ppData_ppMc.GetLowEdge(0),dlm_RatioCk_ppData_ppMc.GetUpEdge(0));
    TH1F* histo_RatioCk_pLData_pLMc = new TH1F("histo_RatioCk_pLData_pLMc","histo_RatioCk_pLData_pLMc",dlm_RatioCk_pLData_pLMc.GetNbins(),dlm_RatioCk_pLData_pLMc.GetLowEdge(0),dlm_RatioCk_pLData_pLMc.GetUpEdge(0));
    TH1F* histo_RatioSe_ppData_ppMc = new TH1F("histo_RatioSe_ppData_ppMc","histo_RatioSe_ppData_ppMc",dlm_RatioSe_ppData_ppMc.GetNbins(),dlm_RatioSe_ppData_ppMc.GetLowEdge(0),dlm_RatioSe_ppData_ppMc.GetUpEdge(0));
    TH1F* histo_RatioSe_pLData_pLMc = new TH1F("histo_RatioSe_pLData_pLMc","histo_RatioSe_pLData_pLMc",dlm_RatioSe_pLData_pLMc.GetNbins(),dlm_RatioSe_pLData_pLMc.GetLowEdge(0),dlm_RatioSe_pLData_pLMc.GetUpEdge(0));
    TH1F* histo_RatioMe_ppData_ppMc = new TH1F("histo_RatioMe_ppData_ppMc","histo_RatioMe_ppData_ppMc",dlm_RatioMe_ppData_ppMc.GetNbins(),dlm_RatioMe_ppData_ppMc.GetLowEdge(0),dlm_RatioMe_ppData_ppMc.GetUpEdge(0));
    TH1F* histo_RatioMe_pLData_pLMc = new TH1F("histo_RatioMe_pLData_pLMc","histo_RatioMe_pLData_pLMc",dlm_RatioMe_pLData_pLMc.GetNbins(),dlm_RatioMe_pLData_pLMc.GetLowEdge(0),dlm_RatioMe_pLData_pLMc.GetUpEdge(0));

    for(unsigned uBin=0; uBin<dlm_DataSe_pp.GetNbins(); uBin++){
        histo_DataSe_pp->SetBinContent(uBin+1,dlm_DataSe_pp.GetBinContent(uBin));
        histo_DataSe_pp->SetBinError(uBin+1,dlm_DataSe_pp.GetBinError(uBin));
    }
    histo_DataSe_pp->Write();

    for(unsigned uBin=0; uBin<dlm_DataSe_pL.GetNbins(); uBin++){
        histo_DataSe_pL->SetBinContent(uBin+1,dlm_DataSe_pL.GetBinContent(uBin));
        histo_DataSe_pL->SetBinError(uBin+1,dlm_DataSe_pL.GetBinError(uBin));
    }
    histo_DataSe_pL->Write();

    for(unsigned uBin=0; uBin<dlm_McSe_pp.GetNbins(); uBin++){
        histo_McSe_pp->SetBinContent(uBin+1,dlm_McSe_pp.GetBinContent(uBin));
        histo_McSe_pp->SetBinError(uBin+1,dlm_McSe_pp.GetBinError(uBin));
    }
    histo_McSe_pp->Write();

    for(unsigned uBin=0; uBin<dlm_McSe_pL.GetNbins(); uBin++){
        histo_McSe_pL->SetBinContent(uBin+1,dlm_McSe_pL.GetBinContent(uBin));
        histo_McSe_pL->SetBinError(uBin+1,dlm_McSe_pL.GetBinError(uBin));
    }
    histo_McSe_pL->Write();

    for(unsigned uBin=0; uBin<dlm_DataMe_pp.GetNbins(); uBin++){
        histo_DataMe_pp->SetBinContent(uBin+1,dlm_DataMe_pp.GetBinContent(uBin));
        histo_DataMe_pp->SetBinError(uBin+1,dlm_DataMe_pp.GetBinError(uBin));
    }
    histo_DataMe_pp->Write();

    for(unsigned uBin=0; uBin<dlm_DataMe_pL.GetNbins(); uBin++){
        histo_DataMe_pL->SetBinContent(uBin+1,dlm_DataMe_pL.GetBinContent(uBin));
        histo_DataMe_pL->SetBinError(uBin+1,dlm_DataMe_pL.GetBinError(uBin));
    }
    histo_DataMe_pL->Write();

    for(unsigned uBin=0; uBin<dlm_McMe_pp.GetNbins(); uBin++){
        histo_McMe_pp->SetBinContent(uBin+1,dlm_McMe_pp.GetBinContent(uBin));
        histo_McMe_pp->SetBinError(uBin+1,dlm_McMe_pp.GetBinError(uBin));
    }
    histo_McMe_pp->Write();

    for(unsigned uBin=0; uBin<dlm_McMe_pL.GetNbins(); uBin++){
        histo_McMe_pL->SetBinContent(uBin+1,dlm_McMe_pL.GetBinContent(uBin));
        histo_McMe_pL->SetBinError(uBin+1,dlm_McMe_pL.GetBinError(uBin));
    }
    histo_McMe_pL->Write();

    for(unsigned uBin=0; uBin<dlm_DataCk_pp.GetNbins(); uBin++){
        histo_DataCk_pp->SetBinContent(uBin+1,dlm_DataCk_pp.GetBinContent(uBin));
        histo_DataCk_pp->SetBinError(uBin+1,dlm_DataCk_pp.GetBinError(uBin));
    }
    /*
    TF1* fit_DataCk_pp = new TF1("fit_DataCk_pp","(1-[0]*exp(-pow(x*[1]/197.327,[2])))*[3]*(1+[4]*x+[5]*x*x+[6]*x*x*x+[7]*x*x*x*x+[8]*pow(x,5.))",300,2000);
    fit_DataCk_pp->SetParameter(0,0.2); fit_DataCk_pp->SetParLimits(0,0,1.0);
    fit_DataCk_pp->SetParameter(1,0.3); fit_DataCk_pp->SetParLimits(1,0.1,0.5);
    //fit_DataCk_pp->SetParameter(2,1.5); fit_DataCk_pp->SetParLimits(2,1,2);
    fit_DataCk_pp->FixParameter(2,2);
    fit_DataCk_pp->SetParameter(3,1); fit_DataCk_pp->SetParLimits(3,-10.,10.);
    fit_DataCk_pp->SetParameter(4,-1e-4); fit_DataCk_pp->SetParLimits(4,-1e-3,1e-3);
    fit_DataCk_pp->SetParameter(5,1e-7); fit_DataCk_pp->SetParLimits(5,-1e-6,1e-6);
    fit_DataCk_pp->SetParameter(6,-1e-10); fit_DataCk_pp->SetParLimits(6,-1e-9,1e-9);
    fit_DataCk_pp->SetParameter(7,1e-14); fit_DataCk_pp->SetParLimits(7,-1e-12,1e-12);
    fit_DataCk_pp->FixParameter(8,0);
    //fit_DataCk_pp->SetParameter(8,0); fit_DataCk_pp->SetParLimits(8,-1e-15,1e-15);
    */
    TF1* fit_DataCk_pp = new TF1("fit_DataCk_pp","(1-[0]*exp(-pow(x*[1]/197.327,[2])))*[3]*(1+[4]*x+[5]*x*x+[6]*x*x*x)*",300,2000);

    histo_DataCk_pp->Fit(fit_DataCk_pp,"S, N, R, M");
    histo_DataCk_pp->Write();
    fit_DataCk_pp->Write();

    for(unsigned uBin=0; uBin<dlm_DataCk_pL.GetNbins(); uBin++){
        histo_DataCk_pL->SetBinContent(uBin+1,dlm_DataCk_pL.GetBinContent(uBin));
        histo_DataCk_pL->SetBinError(uBin+1,dlm_DataCk_pL.GetBinError(uBin));
    }
    histo_DataCk_pL->Write();

    for(unsigned uBin=0; uBin<dlm_McCk_pp.GetNbins(); uBin++){
        histo_McCk_pp->SetBinContent(uBin+1,dlm_McCk_pp.GetBinContent(uBin));
        histo_McCk_pp->SetBinError(uBin+1,dlm_McCk_pp.GetBinError(uBin));
    }
    histo_McCk_pp->Write();

    for(unsigned uBin=0; uBin<dlm_McCk_pL.GetNbins(); uBin++){
        histo_McCk_pL->SetBinContent(uBin+1,dlm_McCk_pL.GetBinContent(uBin));
        histo_McCk_pL->SetBinError(uBin+1,dlm_McCk_pL.GetBinError(uBin));
    }
    histo_McCk_pL->Write();

    for(unsigned uBin=0; uBin<dlm_RatioCk_ppData_ppMc.GetNbins(); uBin++){
        histo_RatioCk_ppData_ppMc->SetBinContent(uBin+1,dlm_RatioCk_ppData_ppMc.GetBinContent(uBin));
        histo_RatioCk_ppData_ppMc->SetBinError(uBin+1,dlm_RatioCk_ppData_ppMc.GetBinError(uBin));
        if(dlm_RatioCk_ppData_ppMc.GetBinContent(uBin)>1000||dlm_RatioCk_ppData_ppMc.GetBinError(uBin)>1000){
            histo_RatioCk_ppData_ppMc->SetBinContent(uBin+1,1);
            histo_RatioCk_ppData_ppMc->SetBinError(uBin+1,1);
        }
    }
    histo_RatioCk_ppData_ppMc->Write();

    for(unsigned uBin=0; uBin<dlm_RatioCk_pLData_pLMc.GetNbins(); uBin++){
        histo_RatioCk_pLData_pLMc->SetBinContent(uBin+1,dlm_RatioCk_pLData_pLMc.GetBinContent(uBin));
        histo_RatioCk_pLData_pLMc->SetBinError(uBin+1,dlm_RatioCk_pLData_pLMc.GetBinError(uBin));
        if(dlm_RatioCk_pLData_pLMc.GetBinContent(uBin)>1000||dlm_RatioCk_pLData_pLMc.GetBinError(uBin)>1000){
            histo_RatioCk_pLData_pLMc->SetBinContent(uBin+1,1);
            histo_RatioCk_pLData_pLMc->SetBinError(uBin+1,1);
        }
    }
    histo_RatioCk_pLData_pLMc->Write();

    for(unsigned uBin=0; uBin<dlm_RatioSe_ppData_ppMc.GetNbins(); uBin++){
        histo_RatioSe_ppData_ppMc->SetBinContent(uBin+1,dlm_RatioSe_ppData_ppMc.GetBinContent(uBin));
        histo_RatioSe_ppData_ppMc->SetBinError(uBin+1,dlm_RatioSe_ppData_ppMc.GetBinError(uBin));
        if(dlm_RatioSe_ppData_ppMc.GetBinContent(uBin)>1000||dlm_RatioSe_ppData_ppMc.GetBinError(uBin)>1000){
            histo_RatioSe_ppData_ppMc->SetBinContent(uBin+1,1);
            histo_RatioSe_ppData_ppMc->SetBinError(uBin+1,1);
        }
    }
    histo_RatioSe_ppData_ppMc->Write();

    for(unsigned uBin=0; uBin<dlm_RatioSe_pLData_pLMc.GetNbins(); uBin++){
        histo_RatioSe_pLData_pLMc->SetBinContent(uBin+1,dlm_RatioSe_pLData_pLMc.GetBinContent(uBin));
        histo_RatioSe_pLData_pLMc->SetBinError(uBin+1,dlm_RatioSe_pLData_pLMc.GetBinError(uBin));
        if(dlm_RatioSe_pLData_pLMc.GetBinContent(uBin)>1000||dlm_RatioSe_pLData_pLMc.GetBinError(uBin)>1000){
            histo_RatioSe_pLData_pLMc->SetBinContent(uBin+1,1);
            histo_RatioSe_pLData_pLMc->SetBinError(uBin+1,1);
        }
    }
    histo_RatioSe_pLData_pLMc->Write();

    for(unsigned uBin=0; uBin<dlm_RatioMe_ppData_ppMc.GetNbins(); uBin++){
        histo_RatioMe_ppData_ppMc->SetBinContent(uBin+1,dlm_RatioMe_ppData_ppMc.GetBinContent(uBin));
        histo_RatioMe_ppData_ppMc->SetBinError(uBin+1,dlm_RatioMe_ppData_ppMc.GetBinError(uBin));
        if(dlm_RatioMe_ppData_ppMc.GetBinContent(uBin)>1000||dlm_RatioMe_ppData_ppMc.GetBinError(uBin)>1000){
            histo_RatioMe_ppData_ppMc->SetBinContent(uBin+1,1);
            histo_RatioMe_ppData_ppMc->SetBinError(uBin+1,1);
        }
    }
    histo_RatioMe_ppData_ppMc->Write();

    for(unsigned uBin=0; uBin<dlm_RatioMe_pLData_pLMc.GetNbins(); uBin++){
        histo_RatioMe_pLData_pLMc->SetBinContent(uBin+1,dlm_RatioMe_pLData_pLMc.GetBinContent(uBin));
        histo_RatioMe_pLData_pLMc->SetBinError(uBin+1,dlm_RatioMe_pLData_pLMc.GetBinError(uBin));
        if(dlm_RatioMe_pLData_pLMc.GetBinContent(uBin)>1000||dlm_RatioMe_pLData_pLMc.GetBinError(uBin)>1000){
            histo_RatioMe_pLData_pLMc->SetBinContent(uBin+1,1);
            histo_RatioMe_pLData_pLMc->SetBinError(uBin+1,1);
        }
    }
    histo_RatioMe_pLData_pLMc->Write();

    delete histo_DataSe_pp;
    delete histo_DataSe_pL;
    delete histo_McSe_pp;
    delete histo_McSe_pL;
    delete histo_DataMe_pp;
    delete histo_DataMe_pL;
    delete histo_McMe_pp;
    delete histo_McMe_pL;
    delete histo_DataCk_pp;
    delete histo_DataCk_pL;
    delete histo_McCk_pp;
    delete histo_McCk_pL;
    delete histo_RatioCk_ppData_ppMc;
    delete histo_RatioCk_pLData_pLMc;
    delete histo_RatioSe_ppData_ppMc;
    delete histo_RatioSe_pLData_pLMc;
    delete histo_RatioMe_ppData_ppMc;
    delete histo_RatioMe_pLData_pLMc;

    delete fit_DataCk_pp;

}


DLM_Histo<float>* dlm_ff_t1;
DLM_Histo<float>* dlm_ff_t2;


double dimi_fraction_fitter_T1(double* x, double* par){
    return  par[0]*dlm_ff_t1->Eval(x);
}
double dimi_fraction_fitter_T2(double* x, double* par){
    return  par[0]*dlm_ff_t2->Eval(x);
}
double dimi_fraction_fitter(double* x, double* par){
    return  (par[0]*dlm_ff_t1->Eval(x)+par[1]*dlm_ff_t2->Eval(x))+
            (par[2]+par[3]*0.5*(*x)+par[4]*(*x)*(*x));
}
double dimi_fractionQS_fitter(double* x, double* par){
    double Templates = par[0]*dlm_ff_t1->Eval(x)+par[1]*dlm_ff_t2->Eval(x);
    //double Correction = par[5]+par[6]*0.5*(*x)+par[7]*(*x)*(*x);
    double Correction = par[5]+par[6]*(*x)+par[7]*(*x)*(*x)+par[8]*(*x)*(*x)*(*x)+par[9]*(*x)*(*x)*(*x)*(*x);
    double Pauli = 1.-par[2]*exp(-pow((*x)*par[3]/197.327,par[4]));
    return (Templates+Correction)*Pauli;
}
double dimi_fractionPolQS_fitter(double* x, double* par){
    double k = *x;
    double Baseline = par[6]*(1.+par[7]*k+par[8]*k*k+par[9]*k*k*k);
    double Pauli = 1.-par[0]*exp(-pow((*x)*par[1]/197.327,par[2]));
    double Peak = 1.+par[3]*TMath::Gaus(k,par[4],par[5],0);
    //double Peak = 1.+par[3]*dlm_ff_t2->Eval(x)+par[4];
    return Baseline*Pauli*Peak;
}


void ALL_CorrectedMC_EXP(){
    TString OutputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/OtherTasks/ALL_CorrectedMC_EXP/";
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
    DataFolder[0] = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/Fast_BBar/Data/";
    //DataFolder[0] = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample9/";
    DataFolder[1] = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/Fast_BBar/Data/";
    DataFolder[2] = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/Fast_BBar/Data/";
    DataFolder[3] = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/pOmega_1/";
    DataFolder[4] = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/Fast_BBar/Data/";
    DataFolder[5] = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/Fast_BBar/Data/";
    DataFolder[6] = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/Fast_BBar/Data/";
    DataFolder[7] = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/Fast_BBar/Data/";
    DataFolder[8] = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/pOmega_1/";
    DataFolder[9] = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/Fast_BBar/Data/";

    TString* McFolder = new TString[NumSpecies];
    McFolder[0] = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/Fast_BBar/Trains_MCAOD/Norm018028/";
    McFolder[1] = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/Fast_BBar/Trains_MCAOD/Norm018028/";
    McFolder[2] = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/Fast_BBar/Trains_MCAOD/Norm018028/";
    McFolder[3] = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/pOmega_1/";
    McFolder[4] = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/Fast_BBar/Trains_MCAOD/Norm018028/";
    McFolder[5] = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/Fast_BBar/Trains_MCAOD/Norm018028/";
    McFolder[6] = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/Fast_BBar/Trains_MCAOD/Norm018028/";
    McFolder[7] = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/Fast_BBar/Trains_MCAOD/Norm018028/";
    McFolder[8] = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/pOmega_1/";
    McFolder[9] = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/Fast_BBar/Trains_MCAOD/Norm018028/";

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

    TString* DataSeList1 = new TString[NumSpecies];
    for(unsigned uSpec=0; uSpec<NumSpecies; uSpec++) DataSeList1[uSpec] = "PairDist";
    DataSeList1[3] = "";DataSeList1[8] = "";

    TString* DataSeList2 = new TString[NumSpecies];
    for(unsigned uSpec=0; uSpec<NumSpecies; uSpec++) DataSeList2[uSpec] = "PairReweighted";
    DataSeList2[3] = "";DataSeList2[8] = "";

    TString* McSeList1 = new TString[NumSpecies];
    for(unsigned uSpec=0; uSpec<NumSpecies; uSpec++) McSeList1[uSpec] = "PairDist";
    McSeList1[3] = "";McSeList1[8] = "";

    TString* McSeList2 = new TString[NumSpecies];
    for(unsigned uSpec=0; uSpec<NumSpecies; uSpec++) McSeList2[uSpec] = "PairReweighted";
    McSeList2[3] = "";McSeList2[8] = "";

    TString* DataSe = new TString[NumSpecies];
    DataSe[0] = "SEDist_Particle0_Particle0_clone_Shifted_FixShifted_Reweighted";
    DataSe[1] = "SEDist_Particle0_Particle2_clone_Shifted_FixShifted_Rebinned_4_Reweighted";
    DataSe[2] = "SEDist_Particle0_Particle4_clone_Shifted_FixShifted_Rebinned_4_Reweighted";
    DataSe[3] = "SE";
    DataSe[4] = "SEDist_Particle2_Particle2_clone_Shifted_FixShifted_Rebinned_4_Reweighted";
    DataSe[5] = "SEDist_Particle0_Particle1_clone_Shifted_FixShifted_Reweighted";
    DataSe[6] = "SEDist_Particle0_Particle3_clone_Shifted_FixShifted_Rebinned_4_Reweighted";
    DataSe[7] = "SEDist_Particle0_Particle5_clone_Shifted_FixShifted_Rebinned_4_Reweighted";
    DataSe[8] = "SE";
    DataSe[9] = "SEDist_Particle2_Particle3_clone_Shifted_FixShifted_Rebinned_4_Reweighted";

    TString* McSe = new TString[NumSpecies];
    McSe[0] = "SEDist_Particle0_Particle0_clone_Shifted_FixShifted_Reweighted";
    McSe[1] = "SEDist_Particle0_Particle2_clone_Shifted_FixShifted_Rebinned_4_Reweighted";
    McSe[2] = "SEDist_Particle0_Particle4_clone_Shifted_FixShifted_Rebinned_4_Reweighted";
    McSe[3] = "SE";
    McSe[4] = "SEDist_Particle2_Particle2_clone_Shifted_FixShifted_Rebinned_4_Reweighted";
    McSe[5] = "SEDist_Particle0_Particle1_clone_Shifted_FixShifted_Reweighted";
    McSe[6] = "SEDist_Particle0_Particle3_clone_Shifted_FixShifted_Rebinned_4_Reweighted";
    McSe[7] = "SEDist_Particle0_Particle5_clone_Shifted_FixShifted_Rebinned_4_Reweighted";
    McSe[8] = "SE";
    McSe[9] = "SEDist_Particle2_Particle3_clone_Shifted_FixShifted_Rebinned_4_Reweighted";

    TString* DataMe = new TString[NumSpecies];
    DataMe[0] = "MEDist_Particle0_Particle0_clone_Shifted_FixShifted_Reweighted";
    DataMe[1] = "MEDist_Particle0_Particle2_clone_Shifted_FixShifted_Rebinned_4_Reweighted";
    DataMe[2] = "MEDist_Particle0_Particle4_clone_Shifted_FixShifted_Rebinned_4_Reweighted";
    DataMe[3] = "ME";
    DataMe[4] = "MEDist_Particle2_Particle2_clone_Shifted_FixShifted_Rebinned_4_Reweighted";
    DataMe[5] = "MEDist_Particle0_Particle1_clone_Shifted_FixShifted_Reweighted";
    DataMe[6] = "MEDist_Particle0_Particle3_clone_Shifted_FixShifted_Rebinned_4_Reweighted";
    DataMe[7] = "MEDist_Particle0_Particle5_clone_Shifted_FixShifted_Rebinned_4_Reweighted";
    DataMe[8] = "ME";
    DataMe[9] = "MEDist_Particle2_Particle3_clone_Shifted_FixShifted_Rebinned_4_Reweighted";

    TString* McMe = new TString[NumSpecies];
    McMe[0] = "MEDist_Particle0_Particle0_clone_Shifted_FixShifted_Reweighted";
    McMe[1] = "MEDist_Particle0_Particle2_clone_Shifted_FixShifted_Rebinned_4_Reweighted";
    McMe[2] = "MEDist_Particle0_Particle4_clone_Shifted_FixShifted_Rebinned_4_Reweighted";
    McMe[3] = "ME";
    McMe[4] = "MEDist_Particle2_Particle2_clone_Shifted_FixShifted_Rebinned_4_Reweighted";
    McMe[5] = "MEDist_Particle0_Particle1_clone_Shifted_FixShifted_Reweighted";
    McMe[6] = "MEDist_Particle0_Particle3_clone_Shifted_FixShifted_Rebinned_4_Reweighted";
    McMe[7] = "MEDist_Particle0_Particle5_clone_Shifted_FixShifted_Rebinned_4_Reweighted";
    McMe[8] = "ME";
    McMe[9] = "MEDist_Particle2_Particle3_clone_Shifted_FixShifted_Rebinned_4_Reweighted";

    DLM_Histo<float>* dlm_DataSe = new DLM_Histo<float> [NumSpecies];
    DLM_Histo<float>* dlm_McSe = new DLM_Histo<float> [NumSpecies];
    DLM_Histo<float>* dlm_DataMe = new DLM_Histo<float> [NumSpecies];
    DLM_Histo<float>* dlm_McMe = new DLM_Histo<float> [NumSpecies];
    DLM_Histo<float>* dlm_DataCk = new DLM_Histo<float> [NumSpecies];
    DLM_Histo<float>* dlm_McCk = new DLM_Histo<float> [NumSpecies];

    DLM_Histo<float>** dlm_RatioCk_Data_Mc = new DLM_Histo<float>* [NumSpecies];
    DLM_Histo<float>** dlm_RatioSe_Data_Mc = new DLM_Histo<float>* [NumSpecies];
    DLM_Histo<float>** dlm_RatioMe_Data_Mc = new DLM_Histo<float>* [NumSpecies];

    DLM_Histo<float>** dlm_RatioCk_Data_Data = new DLM_Histo<float>* [NumSpecies];
    DLM_Histo<float>** dlm_RatioSe_Data_Data = new DLM_Histo<float>* [NumSpecies];
    DLM_Histo<float>** dlm_RatioMe_Data_Data = new DLM_Histo<float>* [NumSpecies];

    DLM_Histo<float>** dlm_RatioCk_Mc_Mc = new DLM_Histo<float>* [NumSpecies];
    DLM_Histo<float>** dlm_RatioSe_Mc_Mc = new DLM_Histo<float>* [NumSpecies];
    DLM_Histo<float>** dlm_RatioMe_Mc_Mc = new DLM_Histo<float>* [NumSpecies];

    for(unsigned uSpec=0; uSpec<NumSpecies; uSpec++){
        dlm_RatioCk_Data_Mc[uSpec] = new DLM_Histo<float> [NumSpecies];
        dlm_RatioSe_Data_Mc[uSpec] = new DLM_Histo<float> [NumSpecies];
        dlm_RatioMe_Data_Mc[uSpec] = new DLM_Histo<float> [NumSpecies];
        dlm_RatioCk_Data_Data[uSpec] = new DLM_Histo<float> [NumSpecies];
        dlm_RatioSe_Data_Data[uSpec] = new DLM_Histo<float> [NumSpecies];
        dlm_RatioMe_Data_Data[uSpec] = new DLM_Histo<float> [NumSpecies];
        dlm_RatioCk_Mc_Mc[uSpec] = new DLM_Histo<float> [NumSpecies];
        dlm_RatioSe_Mc_Mc[uSpec] = new DLM_Histo<float> [NumSpecies];
        dlm_RatioMe_Mc_Mc[uSpec] = new DLM_Histo<float> [NumSpecies];
    }

    for(unsigned uSpec=0; uSpec<NumSpecies; uSpec++){
        unsigned NumBins;
        double kMin;
        double kMax;
        TH1F* hCk_tmp;
        TList* list1_tmp=NULL;
        TList* list2_tmp=NULL;
        TH1F* hSe_tmp;
        TH1F* hMe_tmp;
printf("Opening file %s%s\n",DataFolder[uSpec].Data(),DataFile[uSpec].Data());
        TFile datafile(DataFolder[uSpec]+DataFile[uSpec],"read");
        hCk_tmp = (TH1F*)datafile.Get(DataHisto[uSpec]);
//printf("hCk_tmp = %p\n",hCk_tmp);
//printf("hCk_tmp->GetNbinsX()=%u\n",hCk_tmp->GetNbinsX());
        hCk_tmp->Rebin(RescaleFactor[uSpec]); hCk_tmp->Scale(1./double(RescaleFactor[uSpec]));
//printf("hCk_tmp->GetNbinsX()=%u\n",hCk_tmp->GetNbinsX());
        if(DataSeList1[uSpec]!="") list1_tmp = (TList*)datafile.Get(DataSeList1[uSpec]);
        if(list1_tmp) list2_tmp = (TList*)list1_tmp->FindObject(DataSeList2[uSpec]);
        if(list2_tmp) hSe_tmp = (TH1F*)list2_tmp->FindObject(DataSe[uSpec]);
        else hSe_tmp = (TH1F*)datafile.Get(DataSe[uSpec]);
//printf("hSe_tmp = %p\n",hSe_tmp);
//printf("hSe_tmp->Integral() = %f\n",hSe_tmp->Integral());
//printf("hSe_tmp->GetBinWidth() = %f\n",hSe_tmp->GetBinWidth(1));
        hSe_tmp->Rebin(RescaleFactor[uSpec]); hSe_tmp->Scale(1./hSe_tmp->Integral()/RescaleSeMeAxis[uSpec],"width");
//printf(" hSe_tmp->Integral() = %f\n",hSe_tmp->Integral());
        if(list1_tmp) hMe_tmp = (TH1F*)list2_tmp->FindObject(DataMe[uSpec]);
        else hMe_tmp = (TH1F*)datafile.Get(DataMe[uSpec]);
//printf("hMe_tmp = %p\n",hMe_tmp);
        hMe_tmp->Rebin(RescaleFactor[uSpec]); hMe_tmp->Scale(1./hMe_tmp->Integral()/RescaleSeMeAxis[uSpec],"width");
        NumBins = hCk_tmp->GetNbinsX();
        kMin = hCk_tmp->GetBinLowEdge(1);
        kMax = hCk_tmp->GetXaxis()->GetBinUpEdge(NumBins);
        dlm_DataCk[uSpec].SetUp(1);
        dlm_DataCk[uSpec].SetUp(0,NumBins,kMin,kMax);
        dlm_DataCk[uSpec].Initialize();
        dlm_DataSe[uSpec].SetUp(1);
        dlm_DataSe[uSpec].SetUp(0,NumBins,kMin,kMax);
        dlm_DataSe[uSpec].Initialize();
        dlm_DataMe[uSpec].SetUp(1);
        dlm_DataMe[uSpec].SetUp(0,NumBins,kMin,kMax);
        dlm_DataMe[uSpec].Initialize();
        for(unsigned uBin=0; uBin<NumBins; uBin++){
            dlm_DataCk[uSpec].SetBinContent(uBin,hCk_tmp->GetBinContent(uBin+1));
            dlm_DataCk[uSpec].SetBinError(uBin,hCk_tmp->GetBinError(uBin+1));
            dlm_DataSe[uSpec].SetBinContent(uBin,hSe_tmp->GetBinContent(uBin+1));
            dlm_DataSe[uSpec].SetBinError(uBin,hSe_tmp->GetBinError(uBin+1));
            dlm_DataMe[uSpec].SetBinContent(uBin,hMe_tmp->GetBinContent(uBin+1));
            dlm_DataMe[uSpec].SetBinError(uBin,hMe_tmp->GetBinError(uBin+1));
        }
        datafile.Close();

        list1_tmp=NULL;
        list2_tmp=NULL;

        TFile mcfile(McFolder[uSpec]+McFile[uSpec],"read");
        hCk_tmp = (TH1F*)mcfile.Get(McHisto[uSpec]);
        hCk_tmp->Rebin(RescaleFactor[uSpec]); hCk_tmp->Scale(1./double(RescaleFactor[uSpec]));
        if(McSeList1[uSpec]!="") list1_tmp = (TList*)mcfile.Get(McSeList1[uSpec]);
        if(list1_tmp) list2_tmp = (TList*)list1_tmp->FindObject(McSeList2[uSpec]);
        if(list2_tmp) hSe_tmp = (TH1F*)list2_tmp->FindObject(McSe[uSpec]);
        else hSe_tmp = (TH1F*)mcfile.Get(McSe[uSpec]);
//mcfile.ls();
//printf("McSe=%s\n",McSe[uSpec].Data());
//printf("hSe_tmp = %p\n",hSe_tmp);
        hSe_tmp->Rebin(RescaleFactor[uSpec]); hSe_tmp->Scale(1./hSe_tmp->Integral()/RescaleSeMeAxis[uSpec],"width");
        if(list2_tmp) hMe_tmp = (TH1F*)list2_tmp->FindObject(McMe[uSpec]);
        else hMe_tmp = (TH1F*)mcfile.Get(McMe[uSpec]);
        hMe_tmp->Rebin(RescaleFactor[uSpec]); hMe_tmp->Scale(1./hMe_tmp->Integral()/RescaleSeMeAxis[uSpec],"width");
        NumBins = hCk_tmp->GetNbinsX();
        kMin = hCk_tmp->GetBinLowEdge(1);
        kMax = hCk_tmp->GetXaxis()->GetBinUpEdge(NumBins);
        dlm_McCk[uSpec].SetUp(1);
        dlm_McCk[uSpec].SetUp(0,NumBins,kMin,kMax);
        dlm_McCk[uSpec].Initialize();
        dlm_McSe[uSpec].SetUp(1);
        dlm_McSe[uSpec].SetUp(0,NumBins,kMin,kMax);
        dlm_McSe[uSpec].Initialize();
        dlm_McMe[uSpec].SetUp(1);
        dlm_McMe[uSpec].SetUp(0,NumBins,kMin,kMax);
        dlm_McMe[uSpec].Initialize();
        for(unsigned uBin=0; uBin<NumBins; uBin++){
            dlm_McCk[uSpec].SetBinContent(uBin,hCk_tmp->GetBinContent(uBin+1));
            dlm_McCk[uSpec].SetBinError(uBin,hCk_tmp->GetBinError(uBin+1));
            dlm_McSe[uSpec].SetBinContent(uBin,hSe_tmp->GetBinContent(uBin+1));
            dlm_McSe[uSpec].SetBinError(uBin,hSe_tmp->GetBinError(uBin+1));
            dlm_McMe[uSpec].SetBinContent(uBin,hMe_tmp->GetBinContent(uBin+1));
            dlm_McMe[uSpec].SetBinError(uBin,hMe_tmp->GetBinError(uBin+1));
        }
        mcfile.Close();

        list1_tmp=NULL;
        list2_tmp=NULL;
    }

    for(unsigned uSpec=0; uSpec<NumSpecies; uSpec++){
        for(unsigned uSpec2=0; uSpec2<NumSpecies; uSpec2++){
            dlm_RatioCk_Data_Mc[uSpec][uSpec2] = dlm_DataCk[uSpec];
            dlm_RatioCk_Data_Mc[uSpec][uSpec2] /= dlm_McCk[uSpec2];

            dlm_RatioSe_Data_Mc[uSpec][uSpec2] = dlm_DataSe[uSpec];
            dlm_RatioSe_Data_Mc[uSpec][uSpec2] /= dlm_McSe[uSpec2];

            dlm_RatioMe_Data_Mc[uSpec][uSpec2] = dlm_DataMe[uSpec];
            dlm_RatioMe_Data_Mc[uSpec][uSpec2] /= dlm_McMe[uSpec2];

            dlm_RatioCk_Data_Data[uSpec][uSpec2] = dlm_DataCk[uSpec];
            dlm_RatioCk_Data_Data[uSpec][uSpec2] /= dlm_DataCk[uSpec2];

            dlm_RatioSe_Data_Data[uSpec][uSpec2] = dlm_DataSe[uSpec];
            dlm_RatioSe_Data_Data[uSpec][uSpec2] /= dlm_DataSe[uSpec2];

            dlm_RatioMe_Data_Data[uSpec][uSpec2] = dlm_DataMe[uSpec];
            dlm_RatioMe_Data_Data[uSpec][uSpec2] /= dlm_DataMe[uSpec2];

            dlm_RatioCk_Mc_Mc[uSpec][uSpec2] = dlm_McCk[uSpec];
            dlm_RatioCk_Mc_Mc[uSpec][uSpec2] /= dlm_McCk[uSpec2];

            dlm_RatioSe_Mc_Mc[uSpec][uSpec2] = dlm_McSe[uSpec];
            dlm_RatioSe_Mc_Mc[uSpec][uSpec2] /= dlm_McSe[uSpec2];

            dlm_RatioMe_Mc_Mc[uSpec][uSpec2] = dlm_McMe[uSpec];
            dlm_RatioMe_Mc_Mc[uSpec][uSpec2] /= dlm_McMe[uSpec2];
        }
    }

    TFile fOutput(OutputFolder+"fOutput.root","recreate");
    TH1F** histo_DataSe = new TH1F* [NumSpecies];
    TH1F** histo_McSe = new TH1F* [NumSpecies];
    TH1F** histo_DataMe = new TH1F* [NumSpecies];
    TH1F** histo_McMe = new TH1F* [NumSpecies];
    TH1F** histo_DataCk = new TH1F* [NumSpecies];
    TH1F** histo_McCk = new TH1F* [NumSpecies];

    TH1F*** histo_RatioSe_Data_Mc = new TH1F** [NumSpecies];
    TH1F*** histo_RatioMe_Data_Mc = new TH1F** [NumSpecies];
    TH1F*** histo_RatioCk_Data_Mc = new TH1F** [NumSpecies];
    TH1F*** histo_RatioSe_Data_Data = new TH1F** [NumSpecies];
    TH1F*** histo_RatioMe_Data_Data = new TH1F** [NumSpecies];
    TH1F*** histo_RatioCk_Data_Data = new TH1F** [NumSpecies];
    TH1F*** histo_RatioSe_Mc_Mc = new TH1F** [NumSpecies];
    TH1F*** histo_RatioMe_Mc_Mc = new TH1F** [NumSpecies];
    TH1F*** histo_RatioCk_Mc_Mc = new TH1F** [NumSpecies];

    for(unsigned uSpec=0; uSpec<NumSpecies; uSpec++){

        histo_DataSe[uSpec] = new TH1F("histo_DataSe_"+SpeciesName[uSpec],"histo_DataSe_"+SpeciesName[uSpec],
                                       dlm_DataSe[uSpec].GetNbins(),dlm_DataSe[uSpec].GetLowEdge(0),dlm_DataSe[uSpec].GetUpEdge(0));
        histo_McSe[uSpec] = new TH1F("histo_McSe_"+SpeciesName[uSpec],"histo_McSe_"+SpeciesName[uSpec],
                                       dlm_McSe[uSpec].GetNbins(),dlm_McSe[uSpec].GetLowEdge(0),dlm_McSe[uSpec].GetUpEdge(0));
        histo_DataMe[uSpec] = new TH1F("histo_DataMe_"+SpeciesName[uSpec],"histo_DataMe_"+SpeciesName[uSpec],
                                       dlm_DataMe[uSpec].GetNbins(),dlm_DataMe[uSpec].GetLowEdge(0),dlm_DataMe[uSpec].GetUpEdge(0));
        histo_McMe[uSpec] = new TH1F("histo_McMe_"+SpeciesName[uSpec],"histo_McMe_"+SpeciesName[uSpec],
                                       dlm_McMe[uSpec].GetNbins(),dlm_McMe[uSpec].GetLowEdge(0),dlm_McMe[uSpec].GetUpEdge(0));
        histo_DataCk[uSpec] = new TH1F("histo_DataCk_"+SpeciesName[uSpec],"histo_DataCk_"+SpeciesName[uSpec],
                                       dlm_DataCk[uSpec].GetNbins(),dlm_DataCk[uSpec].GetLowEdge(0),dlm_DataCk[uSpec].GetUpEdge(0));
        histo_McCk[uSpec] = new TH1F("histo_McCk_"+SpeciesName[uSpec],"histo_McCk_"+SpeciesName[uSpec],
                                       dlm_McCk[uSpec].GetNbins(),dlm_McCk[uSpec].GetLowEdge(0),dlm_McCk[uSpec].GetUpEdge(0));
        histo_RatioSe_Data_Mc[uSpec] = new TH1F* [NumSpecies];
        histo_RatioMe_Data_Mc[uSpec] = new TH1F* [NumSpecies];
        histo_RatioCk_Data_Mc[uSpec] = new TH1F* [NumSpecies];
        histo_RatioSe_Data_Data[uSpec] = new TH1F* [NumSpecies];
        histo_RatioMe_Data_Data[uSpec] = new TH1F* [NumSpecies];
        histo_RatioCk_Data_Data[uSpec] = new TH1F* [NumSpecies];
        histo_RatioSe_Mc_Mc[uSpec] = new TH1F* [NumSpecies];
        histo_RatioMe_Mc_Mc[uSpec] = new TH1F* [NumSpecies];
        histo_RatioCk_Mc_Mc[uSpec] = new TH1F* [NumSpecies];
        for(unsigned uSpec2=0; uSpec2<NumSpecies; uSpec2++){
        histo_RatioSe_Data_Mc[uSpec][uSpec2] = new TH1F("histo_RatioSe_Data:"+SpeciesName[uSpec]+"_Mc:"+SpeciesName[uSpec2],
                                                        "histo_RatioSe_Data:"+SpeciesName[uSpec]+"_Mc:"+SpeciesName[uSpec2],
                                       dlm_RatioSe_Data_Mc[uSpec][uSpec2].GetNbins(),dlm_RatioSe_Data_Mc[uSpec][uSpec2].GetLowEdge(0),dlm_RatioSe_Data_Mc[uSpec][uSpec2].GetUpEdge(0));
        histo_RatioMe_Data_Mc[uSpec][uSpec2] = new TH1F("histo_RatioMe_Data:"+SpeciesName[uSpec]+"_Mc:"+SpeciesName[uSpec2],
                                                        "histo_RatioMe_Data:"+SpeciesName[uSpec]+"_Mc:"+SpeciesName[uSpec2],
                                       dlm_RatioMe_Data_Mc[uSpec][uSpec2].GetNbins(),dlm_RatioMe_Data_Mc[uSpec][uSpec2].GetLowEdge(0),dlm_RatioMe_Data_Mc[uSpec][uSpec2].GetUpEdge(0));
        histo_RatioCk_Data_Mc[uSpec][uSpec2] = new TH1F("histo_RatioCk_Data:"+SpeciesName[uSpec]+"_Mc:"+SpeciesName[uSpec2],
                                                        "histo_RatioCk_Data:"+SpeciesName[uSpec]+"_Mc:"+SpeciesName[uSpec2],
                                       dlm_RatioCk_Data_Mc[uSpec][uSpec2].GetNbins(),dlm_RatioCk_Data_Mc[uSpec][uSpec2].GetLowEdge(0),dlm_RatioCk_Data_Mc[uSpec][uSpec2].GetUpEdge(0));
        histo_RatioSe_Data_Data[uSpec][uSpec2] = new TH1F("histo_RatioSe_Data:"+SpeciesName[uSpec]+"_Data:"+SpeciesName[uSpec2],
                                                        "histo_RatioSe_Data:"+SpeciesName[uSpec]+"_Data:"+SpeciesName[uSpec2],
                                       dlm_RatioSe_Data_Data[uSpec][uSpec2].GetNbins(),dlm_RatioSe_Data_Data[uSpec][uSpec2].GetLowEdge(0),dlm_RatioSe_Data_Data[uSpec][uSpec2].GetUpEdge(0));
        histo_RatioMe_Data_Data[uSpec][uSpec2] = new TH1F("histo_RatioMe_Data:"+SpeciesName[uSpec]+"_Data:"+SpeciesName[uSpec2],
                                                        "histo_RatioMe_Data:"+SpeciesName[uSpec]+"_Data:"+SpeciesName[uSpec2],
                                       dlm_RatioMe_Data_Data[uSpec][uSpec2].GetNbins(),dlm_RatioMe_Data_Data[uSpec][uSpec2].GetLowEdge(0),dlm_RatioMe_Data_Data[uSpec][uSpec2].GetUpEdge(0));
        histo_RatioCk_Data_Data[uSpec][uSpec2] = new TH1F("histo_RatioCk_Data:"+SpeciesName[uSpec]+"_Data:"+SpeciesName[uSpec2],
                                                        "histo_RatioCk_Data:"+SpeciesName[uSpec]+"_Data:"+SpeciesName[uSpec2],
                                       dlm_RatioCk_Data_Data[uSpec][uSpec2].GetNbins(),dlm_RatioCk_Data_Data[uSpec][uSpec2].GetLowEdge(0),dlm_RatioCk_Data_Data[uSpec][uSpec2].GetUpEdge(0));
        histo_RatioSe_Mc_Mc[uSpec][uSpec2] = new TH1F("histo_RatioSe_Mc:"+SpeciesName[uSpec]+"_Mc:"+SpeciesName[uSpec2],
                                                        "histo_RatioSe_Mc:"+SpeciesName[uSpec]+"_Mc:"+SpeciesName[uSpec2],
                                       dlm_RatioSe_Mc_Mc[uSpec][uSpec2].GetNbins(),dlm_RatioSe_Mc_Mc[uSpec][uSpec2].GetLowEdge(0),dlm_RatioSe_Mc_Mc[uSpec][uSpec2].GetUpEdge(0));
        histo_RatioMe_Mc_Mc[uSpec][uSpec2] = new TH1F("histo_RatioMe_Mc:"+SpeciesName[uSpec]+"_Mc:"+SpeciesName[uSpec2],
                                                        "histo_RatioMe_Mc:"+SpeciesName[uSpec]+"_Mc:"+SpeciesName[uSpec2],
                                       dlm_RatioMe_Mc_Mc[uSpec][uSpec2].GetNbins(),dlm_RatioMe_Mc_Mc[uSpec][uSpec2].GetLowEdge(0),dlm_RatioMe_Mc_Mc[uSpec][uSpec2].GetUpEdge(0));
        histo_RatioCk_Mc_Mc[uSpec][uSpec2] = new TH1F("histo_RatioCk_Mc:"+SpeciesName[uSpec]+"_Mc:"+SpeciesName[uSpec2],
                                                        "histo_RatioCk_Mc:"+SpeciesName[uSpec]+"_Mc:"+SpeciesName[uSpec2],
                                       dlm_RatioCk_Mc_Mc[uSpec][uSpec2].GetNbins(),dlm_RatioCk_Mc_Mc[uSpec][uSpec2].GetLowEdge(0),dlm_RatioCk_Mc_Mc[uSpec][uSpec2].GetUpEdge(0));

        }
        for(unsigned uBin=0; uBin<dlm_DataSe[uSpec].GetNbins(); uBin++){
            histo_DataSe[uSpec]->SetBinContent(uBin+1,dlm_DataSe[uSpec].GetBinContent(uBin));
            histo_DataSe[uSpec]->SetBinError(uBin+1,dlm_DataSe[uSpec].GetBinError(uBin));
        }
        for(unsigned uBin=0; uBin<dlm_McSe[uSpec].GetNbins(); uBin++){
            histo_McSe[uSpec]->SetBinContent(uBin+1,dlm_McSe[uSpec].GetBinContent(uBin));
            histo_McSe[uSpec]->SetBinError(uBin+1,dlm_McSe[uSpec].GetBinError(uBin));
        }
        for(unsigned uBin=0; uBin<dlm_DataMe[uSpec].GetNbins(); uBin++){
            histo_DataMe[uSpec]->SetBinContent(uBin+1,dlm_DataMe[uSpec].GetBinContent(uBin));
            histo_DataMe[uSpec]->SetBinError(uBin+1,dlm_DataMe[uSpec].GetBinError(uBin));
        }
        for(unsigned uBin=0; uBin<dlm_McMe[uSpec].GetNbins(); uBin++){
            histo_McMe[uSpec]->SetBinContent(uBin+1,dlm_McMe[uSpec].GetBinContent(uBin));
            histo_McMe[uSpec]->SetBinError(uBin+1,dlm_McMe[uSpec].GetBinError(uBin));
        }
        for(unsigned uBin=0; uBin<dlm_DataCk[uSpec].GetNbins(); uBin++){
            histo_DataCk[uSpec]->SetBinContent(uBin+1,dlm_DataCk[uSpec].GetBinContent(uBin));
if(dlm_DataCk[uSpec].GetBinCenter(0,uBin)<250) histo_DataCk[uSpec]->SetBinContent(uBin+1,1.003);
if(dlm_DataCk[uSpec].GetBinCenter(0,uBin)<150) histo_DataCk[uSpec]->SetBinContent(uBin+1,1.0015);
if(dlm_DataCk[uSpec].GetBinCenter(0,uBin)<50) histo_DataCk[uSpec]->SetBinContent(uBin+1,1.0);

if(uSpec==1&&dlm_DataCk[uSpec].GetBinCenter(0,uBin)<305&&dlm_DataCk[uSpec].GetBinCenter(0,uBin)>260) histo_DataCk[uSpec]->SetBinError(uBin+1,0.1);

            histo_DataCk[uSpec]->SetBinError(uBin+1,dlm_DataCk[uSpec].GetBinError(uBin));
//if(dlm_DataCk[uSpec].GetBinCenter(0,uBin)<500||dlm_DataCk[uSpec].GetBinCenter(0,uBin)>2500)
//histo_DataCk[uSpec]->SetBinError(uBin+1,dlm_DataCk[uSpec].GetBinError(uBin));
//else
//histo_DataCk[uSpec]->SetBinError(uBin+1,dlm_DataCk[uSpec].GetBinError(uBin)*4.0);
        }
        for(unsigned uBin=0; uBin<dlm_McCk[uSpec].GetNbins(); uBin++){
            histo_McCk[uSpec]->SetBinContent(uBin+1,dlm_McCk[uSpec].GetBinContent(uBin));
            histo_McCk[uSpec]->SetBinError(uBin+1,dlm_McCk[uSpec].GetBinError(uBin));
        }
        for(unsigned uSpec2=0; uSpec2<NumSpecies; uSpec2++){
            for(unsigned uBin=0; uBin<dlm_RatioSe_Data_Mc[uSpec][uSpec2].GetNbins(); uBin++){
                histo_RatioSe_Data_Mc[uSpec][uSpec2]->SetBinContent(uBin+1,dlm_RatioSe_Data_Mc[uSpec][uSpec2].GetBinContent(uBin));
                histo_RatioSe_Data_Mc[uSpec][uSpec2]->SetBinError(uBin+1,dlm_RatioSe_Data_Mc[uSpec][uSpec2].GetBinError(uBin));
                if(dlm_RatioSe_Data_Mc[uSpec][uSpec2].GetBinContent(uBin)>1000||dlm_RatioSe_Data_Mc[uSpec][uSpec2].GetBinError(uBin)>1000){
                    histo_RatioSe_Data_Mc[uSpec][uSpec2]->SetBinContent(uBin+1,1);
                    histo_RatioSe_Data_Mc[uSpec][uSpec2]->SetBinError(uBin+1,1);
                }
            }
        }
        for(unsigned uSpec2=0; uSpec2<NumSpecies; uSpec2++){
            for(unsigned uBin=0; uBin<dlm_RatioMe_Data_Mc[uSpec][uSpec2].GetNbins(); uBin++){
                histo_RatioMe_Data_Mc[uSpec][uSpec2]->SetBinContent(uBin+1,dlm_RatioMe_Data_Mc[uSpec][uSpec2].GetBinContent(uBin));
                histo_RatioMe_Data_Mc[uSpec][uSpec2]->SetBinError(uBin+1,dlm_RatioMe_Data_Mc[uSpec][uSpec2].GetBinError(uBin));
                if(dlm_RatioMe_Data_Mc[uSpec][uSpec2].GetBinContent(uBin)>1000||dlm_RatioMe_Data_Mc[uSpec][uSpec2].GetBinError(uBin)>1000){
                    histo_RatioMe_Data_Mc[uSpec][uSpec2]->SetBinContent(uBin+1,1);
                    histo_RatioMe_Data_Mc[uSpec][uSpec2]->SetBinError(uBin+1,1);
                }
            }
        }
        for(unsigned uSpec2=0; uSpec2<NumSpecies; uSpec2++){
            for(unsigned uBin=0; uBin<dlm_RatioCk_Data_Mc[uSpec][uSpec2].GetNbins(); uBin++){
                histo_RatioCk_Data_Mc[uSpec][uSpec2]->SetBinContent(uBin+1,dlm_RatioCk_Data_Mc[uSpec][uSpec2].GetBinContent(uBin));
                histo_RatioCk_Data_Mc[uSpec][uSpec2]->SetBinError(uBin+1,dlm_RatioCk_Data_Mc[uSpec][uSpec2].GetBinError(uBin));
                if(dlm_RatioCk_Data_Mc[uSpec][uSpec2].GetBinContent(uBin)>1000||dlm_RatioCk_Data_Mc[uSpec][uSpec2].GetBinError(uBin)>1000){
                    histo_RatioCk_Data_Mc[uSpec][uSpec2]->SetBinContent(uBin+1,1);
                    histo_RatioCk_Data_Mc[uSpec][uSpec2]->SetBinError(uBin+1,1);
                }
            }
        }

        for(unsigned uSpec2=0; uSpec2<NumSpecies; uSpec2++){
            for(unsigned uBin=0; uBin<dlm_RatioSe_Data_Data[uSpec][uSpec2].GetNbins(); uBin++){
                histo_RatioSe_Data_Data[uSpec][uSpec2]->SetBinContent(uBin+1,dlm_RatioSe_Data_Data[uSpec][uSpec2].GetBinContent(uBin));
                histo_RatioSe_Data_Data[uSpec][uSpec2]->SetBinError(uBin+1,dlm_RatioSe_Data_Data[uSpec][uSpec2].GetBinError(uBin));
                if(dlm_RatioSe_Data_Data[uSpec][uSpec2].GetBinContent(uBin)>1000||dlm_RatioSe_Data_Data[uSpec][uSpec2].GetBinError(uBin)>1000){
                    histo_RatioSe_Data_Data[uSpec][uSpec2]->SetBinContent(uBin+1,1);
                    histo_RatioSe_Data_Data[uSpec][uSpec2]->SetBinError(uBin+1,1);
                }
            }
        }
        for(unsigned uSpec2=0; uSpec2<NumSpecies; uSpec2++){
            for(unsigned uBin=0; uBin<dlm_RatioMe_Data_Data[uSpec][uSpec2].GetNbins(); uBin++){
                histo_RatioMe_Data_Data[uSpec][uSpec2]->SetBinContent(uBin+1,dlm_RatioMe_Data_Data[uSpec][uSpec2].GetBinContent(uBin));
                histo_RatioMe_Data_Data[uSpec][uSpec2]->SetBinError(uBin+1,dlm_RatioMe_Data_Data[uSpec][uSpec2].GetBinError(uBin));
                if(dlm_RatioMe_Data_Data[uSpec][uSpec2].GetBinContent(uBin)>1000||dlm_RatioMe_Data_Data[uSpec][uSpec2].GetBinError(uBin)>1000){
                    histo_RatioMe_Data_Data[uSpec][uSpec2]->SetBinContent(uBin+1,1);
                    histo_RatioMe_Data_Data[uSpec][uSpec2]->SetBinError(uBin+1,1);
                }
            }
        }
        for(unsigned uSpec2=0; uSpec2<NumSpecies; uSpec2++){
            for(unsigned uBin=0; uBin<dlm_RatioCk_Data_Data[uSpec][uSpec2].GetNbins(); uBin++){
                histo_RatioCk_Data_Data[uSpec][uSpec2]->SetBinContent(uBin+1,dlm_RatioCk_Data_Data[uSpec][uSpec2].GetBinContent(uBin));
                histo_RatioCk_Data_Data[uSpec][uSpec2]->SetBinError(uBin+1,dlm_RatioCk_Data_Data[uSpec][uSpec2].GetBinError(uBin));
                if(dlm_RatioCk_Data_Data[uSpec][uSpec2].GetBinContent(uBin)>1000||dlm_RatioCk_Data_Data[uSpec][uSpec2].GetBinError(uBin)>1000){
                    histo_RatioCk_Data_Data[uSpec][uSpec2]->SetBinContent(uBin+1,1);
                    histo_RatioCk_Data_Data[uSpec][uSpec2]->SetBinError(uBin+1,1);
                }
            }
        }

        for(unsigned uSpec2=0; uSpec2<NumSpecies; uSpec2++){
            for(unsigned uBin=0; uBin<dlm_RatioSe_Mc_Mc[uSpec][uSpec2].GetNbins(); uBin++){
                histo_RatioSe_Mc_Mc[uSpec][uSpec2]->SetBinContent(uBin+1,dlm_RatioSe_Mc_Mc[uSpec][uSpec2].GetBinContent(uBin));
                histo_RatioSe_Mc_Mc[uSpec][uSpec2]->SetBinError(uBin+1,dlm_RatioSe_Mc_Mc[uSpec][uSpec2].GetBinError(uBin));
                if(dlm_RatioSe_Mc_Mc[uSpec][uSpec2].GetBinContent(uBin)>1000||dlm_RatioSe_Mc_Mc[uSpec][uSpec2].GetBinError(uBin)>1000){
                    histo_RatioSe_Mc_Mc[uSpec][uSpec2]->SetBinContent(uBin+1,1);
                    histo_RatioSe_Mc_Mc[uSpec][uSpec2]->SetBinError(uBin+1,1);
                }
            }
        }
        for(unsigned uSpec2=0; uSpec2<NumSpecies; uSpec2++){
            for(unsigned uBin=0; uBin<dlm_RatioMe_Mc_Mc[uSpec][uSpec2].GetNbins(); uBin++){
                histo_RatioMe_Mc_Mc[uSpec][uSpec2]->SetBinContent(uBin+1,dlm_RatioMe_Mc_Mc[uSpec][uSpec2].GetBinContent(uBin));
                histo_RatioMe_Mc_Mc[uSpec][uSpec2]->SetBinError(uBin+1,dlm_RatioMe_Mc_Mc[uSpec][uSpec2].GetBinError(uBin));
                if(dlm_RatioMe_Mc_Mc[uSpec][uSpec2].GetBinContent(uBin)>1000||dlm_RatioMe_Mc_Mc[uSpec][uSpec2].GetBinError(uBin)>1000){
                    histo_RatioMe_Mc_Mc[uSpec][uSpec2]->SetBinContent(uBin+1,1);
                    histo_RatioMe_Mc_Mc[uSpec][uSpec2]->SetBinError(uBin+1,1);
                }
            }
        }
        for(unsigned uSpec2=0; uSpec2<NumSpecies; uSpec2++){
            for(unsigned uBin=0; uBin<dlm_RatioCk_Mc_Mc[uSpec][uSpec2].GetNbins(); uBin++){
                histo_RatioCk_Mc_Mc[uSpec][uSpec2]->SetBinContent(uBin+1,dlm_RatioCk_Mc_Mc[uSpec][uSpec2].GetBinContent(uBin));
                histo_RatioCk_Mc_Mc[uSpec][uSpec2]->SetBinError(uBin+1,dlm_RatioCk_Mc_Mc[uSpec][uSpec2].GetBinError(uBin));
                if(dlm_RatioCk_Mc_Mc[uSpec][uSpec2].GetBinContent(uBin)>1000||dlm_RatioCk_Mc_Mc[uSpec][uSpec2].GetBinError(uBin)>1000){
                    histo_RatioCk_Mc_Mc[uSpec][uSpec2]->SetBinContent(uBin+1,1);
                    histo_RatioCk_Mc_Mc[uSpec][uSpec2]->SetBinError(uBin+1,1);
                }
            }
        }

//histo_DataSe[uSpec]->Rebin(RescaleFactor[uSpec]); histo_DataSe[uSpec]->Scale(1./double(RescaleFactor[uSpec]));
        histo_DataSe[uSpec]->Write();
        histo_McSe[uSpec]->Write();
        histo_DataMe[uSpec]->Write();
        histo_McMe[uSpec]->Write();
        histo_DataCk[uSpec]->Write();
        histo_McCk[uSpec]->Write();
        for(unsigned uSpec2=0; uSpec2<NumSpecies; uSpec2++){
            histo_RatioSe_Data_Mc[uSpec][uSpec2]->Write();
            histo_RatioMe_Data_Mc[uSpec][uSpec2]->Write();
            histo_RatioCk_Data_Mc[uSpec][uSpec2]->Write();
            histo_RatioSe_Data_Data[uSpec][uSpec2]->Write();
            histo_RatioMe_Data_Data[uSpec][uSpec2]->Write();
            histo_RatioCk_Data_Data[uSpec][uSpec2]->Write();
            histo_RatioSe_Mc_Mc[uSpec][uSpec2]->Write();
            histo_RatioMe_Mc_Mc[uSpec][uSpec2]->Write();
            histo_RatioCk_Mc_Mc[uSpec][uSpec2]->Write();
        }

    }


    //OTHER STUFF ///////////////////////////////////////////

    TF1** fitTemplate = new TF1* [NumSpecies];
    TF1** fitTemplatePert = new TF1* [NumSpecies];
    TF1** fitTemplate1 = new TF1* [NumSpecies];
    TF1** fitTemplate2 = new TF1* [NumSpecies];
    TF1** fitTemplateFull = new TF1* [NumSpecies];
    TH1F** hCorrected = new TH1F* [NumSpecies];
    TH1F** hOriginal = new TH1F* [NumSpecies];

    TF1** fitMickeyMouse = new TF1* [NumSpecies];
    TF1** fitMickeyMouseBL = new TF1* [NumSpecies];
    TF1** fitMickeyMouseQS = new TF1* [NumSpecies];
    TF1** fitMickeyMouseJet = new TF1* [NumSpecies];
    TF1** fitMickeyMouseNoJet = new TF1* [NumSpecies];
    TH1F** hMickeyMouseCorrected = new TH1F* [NumSpecies];

    for(unsigned uSpec=0; uSpec<NumSpecies; uSpec++){
        fitTemplate[uSpec]=NULL;
        fitTemplatePert[uSpec]=NULL;
        fitTemplate1[uSpec]=NULL;
        fitTemplate2[uSpec]=NULL;
        fitTemplateFull[uSpec]=NULL;
        hCorrected[uSpec]=NULL;
        hOriginal[uSpec]=NULL;
        fitMickeyMouse[uSpec]=NULL;
        fitMickeyMouseBL[uSpec]=NULL;
        fitMickeyMouseQS[uSpec]=NULL;
        fitMickeyMouseJet[uSpec]=NULL;
        fitMickeyMouseNoJet[uSpec]=NULL;
        hMickeyMouseCorrected[uSpec]=NULL;
    }

    //fit baryon-antibaryon
    //using templates from MC
    for(unsigned uSpec=NumSpecies/2; uSpec<NumSpecies; uSpec++){
        dlm_ff_t1 = &dlm_McCk[uSpec-NumSpecies/2];//baryon baryon template
        dlm_ff_t2 = &dlm_McCk[uSpec];//baryon antibaryon template

        fitTemplate[uSpec] = new TF1("fitTemplate_"+SpeciesName[uSpec],dimi_fraction_fitter,280,3500,5);

        fitTemplate[uSpec]->SetParameter(0,0.7);
        fitTemplate[uSpec]->SetParLimits(0,0.05,1.4);

        fitTemplate[uSpec]->SetParameter(1,0.15);
        fitTemplate[uSpec]->SetParLimits(1,0.01,1.00);

        fitTemplate[uSpec]->SetParameter(2,0.0);
        fitTemplate[uSpec]->SetParLimits(2,-0.15,0.15);

        fitTemplate[uSpec]->SetParameter(3,0.0);
        fitTemplate[uSpec]->SetParLimits(3,-1e-4*100,1e-4*100);

        fitTemplate[uSpec]->SetParameter(4,0.0);
        fitTemplate[uSpec]->SetParLimits(4,-1e-7*100,1e-7*100);

        histo_DataCk[uSpec]->Fit(fitTemplate[uSpec],"S, N, R, M");

        fitTemplate1[uSpec] = new TF1("fitTemplatePP_"+SpeciesName[uSpec],dimi_fraction_fitter_T1,0,4500,1);
        fitTemplate1[uSpec]->SetParameter(0,fitTemplate[uSpec]->GetParameter(0));

        fitTemplate2[uSpec] = new TF1("fitTemplatePAP_"+SpeciesName[uSpec],dimi_fraction_fitter_T2,0,4500,1);
        fitTemplate2[uSpec]->SetParameter(0,fitTemplate[uSpec]->GetParameter(1));

        fitTemplatePert[uSpec] = new TF1("fitTemplatePert_"+SpeciesName[uSpec],"[0]+[1]*0.5*x+[2]*x*x",0,4500);
        fitTemplatePert[uSpec]->SetParameter(0,fitTemplate[uSpec]->GetParameter(2));
        fitTemplatePert[uSpec]->SetParameter(1,fitTemplate[uSpec]->GetParameter(3));
        fitTemplatePert[uSpec]->SetParameter(2,fitTemplate[uSpec]->GetParameter(4));

        fitTemplateFull[uSpec] = new TF1("fitTemplateFull_"+SpeciesName[uSpec],dimi_fraction_fitter,0,4500,5);
        fitTemplateFull[uSpec]->SetParameter(0,fitTemplate[uSpec]->GetParameter(0));
        fitTemplateFull[uSpec]->SetParameter(1,fitTemplate[uSpec]->GetParameter(1));
        fitTemplateFull[uSpec]->SetParameter(2,fitTemplate[uSpec]->GetParameter(2));
        fitTemplateFull[uSpec]->SetParameter(3,fitTemplate[uSpec]->GetParameter(3));
        fitTemplateFull[uSpec]->SetParameter(4,fitTemplate[uSpec]->GetParameter(4));

        fitTemplate[uSpec]->SetNpx(1024);
        fitTemplatePert[uSpec]->SetNpx(1024);
        fitTemplateFull[uSpec]->SetNpx(1024);

        fitTemplate[uSpec]->Write();
        fitTemplate1[uSpec]->Write();
        fitTemplate2[uSpec]->Write();
        fitTemplatePert[uSpec]->Write();
        fitTemplateFull[uSpec]->Write();

        hCorrected[uSpec] = new TH1F("hCorrected_"+SpeciesName[uSpec],"hCorrected_"+SpeciesName[uSpec],histo_DataCk[uSpec]->GetNbinsX(),
                                    histo_DataCk[uSpec]->GetBinLowEdge(1),histo_DataCk[uSpec]->GetXaxis()->GetBinUpEdge(histo_DataCk[uSpec]->GetNbinsX()));

        for(unsigned uBin=1; uBin<=histo_DataCk[uSpec]->GetNbinsX(); uBin++){
            double MOM = hCorrected[uSpec]->GetBinCenter(uBin);
            //hCorrected[uSpec]->SetBinContent(uBin,histo_DataCk[5]->GetBinContent(uBin)/fitTemplateFull[uSpec]->Eval(MOM));
            //hCorrected[uSpec]->SetBinError(uBin,histo_DataCk[uSpec]->GetBinError(uBin)/fitTemplateFull[uSpec]->Eval(MOM));
            hCorrected[uSpec]->SetBinContent(uBin,dlm_DataCk[uSpec].GetBinContent(uBin-1)/fitTemplateFull[uSpec]->Eval(MOM));
            hCorrected[uSpec]->SetBinError(uBin,dlm_DataCk[uSpec].GetBinError(uBin-1)/fitTemplateFull[uSpec]->Eval(MOM));
            if(dlm_DataCk[uSpec].GetBinContent(uBin-1)/fitTemplateFull[uSpec]->Eval(MOM)!=dlm_DataCk[uSpec].GetBinContent(uBin-1)/fitTemplateFull[uSpec]->Eval(MOM)||
               dlm_DataCk[uSpec].GetBinContent(uBin-1)/fitTemplateFull[uSpec]->Eval(MOM)>100){
            hCorrected[uSpec]->SetBinContent(uBin,100);
            hCorrected[uSpec]->SetBinError(uBin,100);
               }
            //
        }

        hCorrected[uSpec]->Write();


    }

    //fit baryon-baryon
    //using Mickey Mouse fit with QS
    for(unsigned uSpec=0; uSpec<NumSpecies; uSpec++){
printf("uSpec = %u\n",uSpec);
        //double Baseline = par[6]*(1.+par[7]*k+par[8]*k*k+par[9]*k*k*k);
        //double Pauli = 1.-par[0]*exp(-pow((*x)*par[1]/197.327,par[2]));
        //double Peak = 1.+par[3]*TMath::Gaus(k,par[4],par[5],1);
        //return Baseline*Pauli*Peak;
        //different fit range for baryon-antibaryon
        if(uSpec>=NumSpecies/2) fitMickeyMouse[uSpec] = new TF1("fitMickeyMouse_"+SpeciesName[uSpec],dimi_fractionPolQS_fitter,280,2500,10);
        //for baryon-baryon
        else fitMickeyMouse[uSpec] = new TF1("fitMickeyMouse_"+SpeciesName[uSpec],dimi_fractionPolQS_fitter,0,2500,10);

        fitMickeyMouse[uSpec]->SetParameter(0,0.1);
        fitMickeyMouse[uSpec]->SetParLimits(0,0.0,1.0);

        fitMickeyMouse[uSpec]->SetParameter(1,0.3);
        fitMickeyMouse[uSpec]->SetParLimits(1,0.1,0.5);

        fitMickeyMouse[uSpec]->FixParameter(2,2.0);

        //no QS for baryon-antibaryon
        if(uSpec>=NumSpecies/2){
            fitMickeyMouse[uSpec]->FixParameter(0,0.0);
            fitMickeyMouse[uSpec]->FixParameter(1,1.0);
            fitMickeyMouse[uSpec]->FixParameter(2,2.0);
        }

        fitMickeyMouse[uSpec]->SetParameter(3,1.0);
        fitMickeyMouse[uSpec]->SetParLimits(3,0.01,10.0);

        fitMickeyMouse[uSpec]->SetParameter(4,200.);
        fitMickeyMouse[uSpec]->SetParLimits(4,60.,400.);
//fitMickeyMouse[uSpec]->SetParameter(4,150.);
//fitMickeyMouse[uSpec]->SetParLimits(4,100.,200.);

        fitMickeyMouse[uSpec]->SetParameter(5,200.);
        fitMickeyMouse[uSpec]->SetParLimits(5,80.,600.);
//fitMickeyMouse[uSpec]->SetParameter(5,300.);
//fitMickeyMouse[uSpec]->SetParLimits(5,200.,400.);

        fitMickeyMouse[uSpec]->SetParameter(6,1.0);
        fitMickeyMouse[uSpec]->SetParLimits(6,0.0,2.0);

        fitMickeyMouse[uSpec]->SetParameter(7,0.0);
        fitMickeyMouse[uSpec]->SetParLimits(7,-1e-2,1e-2);

        fitMickeyMouse[uSpec]->SetParameter(8,0.0);
        fitMickeyMouse[uSpec]->SetParLimits(8,-1e-4,1e-4);

        fitMickeyMouse[uSpec]->FixParameter(9,0.0);

        histo_DataCk[uSpec]->Fit(fitMickeyMouse[uSpec],"S, N, R, M");

        fitMickeyMouseQS[uSpec] = new TF1("fitMickeyMouseQS_"+SpeciesName[uSpec],"1-[0]*exp(-pow(x*[1]/197.327,[2]))",0,4500);
        fitMickeyMouseQS[uSpec]->SetParameter(0,fitMickeyMouse[uSpec]->GetParameter(0));
        fitMickeyMouseQS[uSpec]->SetParameter(1,fitMickeyMouse[uSpec]->GetParameter(1));
        fitMickeyMouseQS[uSpec]->SetParameter(2,fitMickeyMouse[uSpec]->GetParameter(2));

        fitMickeyMouseJet[uSpec] = new TF1("fitMickeyMouseJet_"+SpeciesName[uSpec],"1.+[0]*TMath::Gaus(x,[1],[2],0)",0,4500);
        fitMickeyMouseJet[uSpec]->SetParameter(0,fitMickeyMouse[uSpec]->GetParameter(3));
        fitMickeyMouseJet[uSpec]->SetParameter(1,fitMickeyMouse[uSpec]->GetParameter(4));
        fitMickeyMouseJet[uSpec]->SetParameter(2,fitMickeyMouse[uSpec]->GetParameter(5));

        fitMickeyMouseBL[uSpec] = new TF1("fitMickeyMouseBL_"+SpeciesName[uSpec],"[0]*(1.+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x)",0,4500);
        fitMickeyMouseBL[uSpec]->SetParameter(0,fitMickeyMouse[uSpec]->GetParameter(6));
        fitMickeyMouseBL[uSpec]->SetParameter(1,fitMickeyMouse[uSpec]->GetParameter(7));
        fitMickeyMouseBL[uSpec]->SetParameter(2,fitMickeyMouse[uSpec]->GetParameter(8));
        fitMickeyMouseBL[uSpec]->SetParameter(3,fitMickeyMouse[uSpec]->GetParameter(9));
        fitMickeyMouseBL[uSpec]->SetParameter(4,0);

        fitMickeyMouseNoJet[uSpec] = new TF1("fitMickeyMouseNoJet_"+SpeciesName[uSpec],"(1-[0]*exp(-pow(x*[1]/197.327,[2])))*[3]*(1.+[4]*x+[5]*x*x+[6]*x*x*x+[7]*x*x*x*x)",0,4500);
        fitMickeyMouseNoJet[uSpec]->SetParameter(0,fitMickeyMouse[uSpec]->GetParameter(0));
        fitMickeyMouseNoJet[uSpec]->SetParameter(1,fitMickeyMouse[uSpec]->GetParameter(1));
        fitMickeyMouseNoJet[uSpec]->SetParameter(2,fitMickeyMouse[uSpec]->GetParameter(2));
        fitMickeyMouseNoJet[uSpec]->SetParameter(3,fitMickeyMouse[uSpec]->GetParameter(6));
        fitMickeyMouseNoJet[uSpec]->SetParameter(4,fitMickeyMouse[uSpec]->GetParameter(7));
        fitMickeyMouseNoJet[uSpec]->SetParameter(5,fitMickeyMouse[uSpec]->GetParameter(8));
        fitMickeyMouseNoJet[uSpec]->SetParameter(6,fitMickeyMouse[uSpec]->GetParameter(9));
        fitMickeyMouseNoJet[uSpec]->SetParameter(7,0);


        fitMickeyMouse[uSpec]->SetNpx(1024);
        fitMickeyMouseQS[uSpec]->SetNpx(1024);
        fitMickeyMouseJet[uSpec]->SetNpx(1024);
        fitMickeyMouseBL[uSpec]->SetNpx(1024);
        fitMickeyMouseNoJet[uSpec]->SetNpx(1024);

        fitMickeyMouse[uSpec]->Write();
        fitMickeyMouseQS[uSpec]->Write();
        fitMickeyMouseJet[uSpec]->Write();
        fitMickeyMouseBL[uSpec]->Write();
        fitMickeyMouseNoJet[uSpec]->Write();

        hMickeyMouseCorrected[uSpec] = new TH1F("hMickeyMouseCorrected_"+SpeciesName[uSpec],"hMickeyMouseCorrected_"+SpeciesName[uSpec],histo_DataCk[uSpec]->GetNbinsX(),
                                    histo_DataCk[uSpec]->GetBinLowEdge(1),histo_DataCk[uSpec]->GetXaxis()->GetBinUpEdge(histo_DataCk[uSpec]->GetNbinsX()));
printf("HELLO\n");
        for(unsigned uBin=1; uBin<=histo_DataCk[uSpec]->GetNbinsX(); uBin++){
            double MOM = hMickeyMouseCorrected[uSpec]->GetBinCenter(uBin);
            //hCorrected[uSpec]->SetBinContent(uBin,histo_DataCk[5]->GetBinContent(uBin)/fitTemplateFull[uSpec]->Eval(MOM));
            //hCorrected[uSpec]->SetBinError(uBin,histo_DataCk[uSpec]->GetBinError(uBin)/fitTemplateFull[uSpec]->Eval(MOM));
            hMickeyMouseCorrected[uSpec]->SetBinContent(uBin,dlm_DataCk[uSpec].GetBinContent(uBin-1)/fitMickeyMouse[uSpec]->Eval(MOM));
            hMickeyMouseCorrected[uSpec]->SetBinError(uBin,dlm_DataCk[uSpec].GetBinError(uBin-1)/fitMickeyMouse[uSpec]->Eval(MOM));
            if(dlm_DataCk[uSpec].GetBinContent(uBin-1)/fitMickeyMouse[uSpec]->Eval(MOM)!=dlm_DataCk[uSpec].GetBinContent(uBin-1)/fitMickeyMouse[uSpec]->Eval(MOM)||
               dlm_DataCk[uSpec].GetBinContent(uBin-1)/fitMickeyMouse[uSpec]->Eval(MOM)>100){
            hMickeyMouseCorrected[uSpec]->SetBinContent(uBin,100);
            hMickeyMouseCorrected[uSpec]->SetBinError(uBin,100);
               }
            //
        }
printf("HELLO2\n");
        hMickeyMouseCorrected[uSpec]->Write();

        hOriginal[uSpec] = new TH1F("hOriginal_"+SpeciesName[uSpec],"hOriginal_"+SpeciesName[uSpec],histo_DataCk[uSpec]->GetNbinsX(),
                                    histo_DataCk[uSpec]->GetBinLowEdge(1),histo_DataCk[uSpec]->GetXaxis()->GetBinUpEdge(histo_DataCk[uSpec]->GetNbinsX()));
        for(unsigned uBin=1; uBin<=histo_DataCk[uSpec]->GetNbinsX(); uBin++){
            double MOM = hOriginal[uSpec]->GetBinCenter(uBin);
            if(fitTemplateFull[uSpec]&&fitTemplateFull[uSpec]->Eval(MOM)<1e-8){
                hOriginal[uSpec]->SetBinContent(uBin,100);
                hOriginal[uSpec]->SetBinError(uBin,100);
            }
            else if(fitMickeyMouse[uSpec]&&fitMickeyMouse[uSpec]->Eval(MOM)<1e-8){
                hOriginal[uSpec]->SetBinContent(uBin,100);
                hOriginal[uSpec]->SetBinError(uBin,100);
            }
            else{
                hOriginal[uSpec]->SetBinContent(uBin,dlm_DataCk[uSpec].GetBinContent(uBin-1));
                hOriginal[uSpec]->SetBinError(uBin,dlm_DataCk[uSpec].GetBinError(uBin-1));
            }
        }
        hOriginal[uSpec]->Write();
    }



/*
    //FIT pAp
    dlm_ff_t1 = &dlm_McCk[0];
    dlm_ff_t2 = &dlm_McCk[5];
    TF1* fit_pAp = new TF1("fit_pAp",dimi_fraction_fitter,350,3500,5);
    //TF1* fit_pAp = new TF1("fit_pAp",dimi_fraction_fitter,350,3500,8);

    fit_pAp->SetParameter(0,0.7);
    fit_pAp->SetParLimits(0,0.3,1.1);

    fit_pAp->SetParameter(1,0.15);
    fit_pAp->SetParLimits(1,0.05,0.65);

    fit_pAp->SetParameter(2,0.0);
    fit_pAp->SetParLimits(2,-0.1*5,0.1*5);

    fit_pAp->SetParameter(3,0.0);
    fit_pAp->SetParLimits(3,-1e-4*5,1e-4*5);

    fit_pAp->SetParameter(4,0.0);
    fit_pAp->SetParLimits(4,-1e-7*5,1e-7*5);

    //fit_pAp->FixParameter(0,0.55);
    //fit_pAp->FixParameter(1,0.45);
    //fit_pAp->FixParameter(2,0);
    //fit_pAp->FixParameter(3,0);
    //fit_pAp->FixParameter(4,0);

    histo_DataCk[5]->Fit(fit_pAp,"S, N, R, M");

    TF1* fitPert_pAp = new TF1("fitPert_pAp","[0]+[1]*0.5*x+[2]*x*x",0,4500);
    fitPert_pAp->SetParameter(0,fit_pAp->GetParameter(2));
    fitPert_pAp->SetParameter(1,fit_pAp->GetParameter(3));
    fitPert_pAp->SetParameter(2,fit_pAp->GetParameter(4));

    TF1* fitFull_pAp = new TF1("fitFull_pAp",dimi_fraction_fitter,0,4500,5);
    fitFull_pAp->SetParameter(0,fit_pAp->GetParameter(0));
    fitFull_pAp->SetParameter(1,fit_pAp->GetParameter(1));
    fitFull_pAp->SetParameter(2,fit_pAp->GetParameter(2));
    fitFull_pAp->SetParameter(3,fit_pAp->GetParameter(3));
    fitFull_pAp->SetParameter(4,fit_pAp->GetParameter(4));

    fit_pAp->Write();
    fitPert_pAp->Write();
    fitFull_pAp->Write();
    TH1F* hCorrected_pAp = new TH1F("hCorrected_pAp","hCorrected_pAp",histo_DataCk[5]->GetNbinsX(),
                                    histo_DataCk[5]->GetBinLowEdge(1),histo_DataCk[5]->GetXaxis()->GetBinUpEdge(histo_DataCk[5]->GetNbinsX()));
    for(unsigned uBin=1; uBin<=histo_DataCk[5]->GetNbinsX(); uBin++){
        double MOM = hCorrected_pAp->GetBinCenter(uBin);
        //hCorrected_pAp->SetBinContent(uBin,histo_DataCk[5]->GetBinContent(uBin)/fitFull_pAp->Eval(MOM));
        //hCorrected_pAp->SetBinError(uBin,histo_DataCk[5]->GetBinError(uBin)/fitFull_pAp->Eval(MOM));
        hCorrected_pAp->SetBinContent(uBin,dlm_DataCk[5].GetBinContent(uBin-1)/fitFull_pAp->Eval(MOM));
        hCorrected_pAp->SetBinError(uBin,dlm_DataCk[5].GetBinError(uBin-1)/fitFull_pAp->Eval(MOM));
        if(dlm_DataCk[5].GetBinContent(uBin-1)/fitFull_pAp->Eval(MOM)!=dlm_DataCk[5].GetBinContent(uBin-1)/fitFull_pAp->Eval(MOM)||
           dlm_DataCk[5].GetBinContent(uBin-1)/fitFull_pAp->Eval(MOM)>100){
        hCorrected_pAp->SetBinContent(uBin,100);
        hCorrected_pAp->SetBinError(uBin,100);
           }
        //
    }
    hCorrected_pAp->Write();
*/

    //FIT ppp

TH1F* hWHICHONE = histo_DataCk[0];
DLM_Histo<float>& dWHICHONE = dlm_DataCk[0];
//TH1F* hWHICHONE = histo_McCk[0];
//DLM_Histo<float>& dWHICHONE = dlm_McCk[0];

    dlm_ff_t1 = &dlm_McCk[0];
    dlm_ff_t2 = &dlm_McCk[5];
    //TF1* fit_pp = new TF1("fit_pp",dimi_fractionQS_fitter,200,3500,10);
    TF1* fit_pp = new TF1("fit_pp",dimi_fractionPolQS_fitter,0,2500,10);

/*
    fit_pp->FixParameter(0,0.0);
    fit_pp->FixParameter(1,0.0);
    fit_pp->FixParameter(2,0.25);
    fit_pp->FixParameter(3,0.3);
    fit_pp->FixParameter(4,2);
    fit_pp->FixParameter(5,0);
    fit_pp->FixParameter(6,0);
    fit_pp->FixParameter(7,0);


    fit_pp->SetParameter(0,1.0);
    fit_pp->SetParLimits(0,0.3,3.0);
    fit_pp->SetParameter(1,1.0);
    fit_pp->SetParLimits(1,0.5,3.0);
    fit_pp->SetParameter(2,0.1);
    fit_pp->SetParLimits(2,0,1.0);
    fit_pp->SetParameter(3,0.3);
    fit_pp->SetParLimits(3,0.1,0.5);
    //fit_pp->SetParameter(4,2);
    //fit_pp->SetParLimits(4,1,2);

    fit_pp->SetParameter(5,0.3);
    fit_pp->SetParLimits(5,-1.0,1.0);

    fit_pp->SetParameter(6,0.0);
    fit_pp->SetParLimits(6,-1e-5,1e-5);

    fit_pp->SetParameter(7,0.0);
    fit_pp->SetParLimits(7,-1e-9,1e-9);

    fit_pp->SetParameter(8,0.0);
    fit_pp->SetParLimits(8,-1e-13,1e-13);

    fit_pp->SetParameter(9,0.0);
    fit_pp->SetParLimits(9,-1e-17,1e-17);
*/

    fit_pp->SetParameter(0,0.1);
    fit_pp->SetParLimits(0,0,1.0);
//fit_pp->SetParameter(0,0.67);
    fit_pp->SetParameter(1,0.3);
    fit_pp->SetParLimits(1,0.1,0.5);
    fit_pp->SetParameter(2,2);
    fit_pp->SetParLimits(2,1,2);
fit_pp->FixParameter(2,2);

    fit_pp->SetParameter(3,10);
    fit_pp->SetParLimits(3,0.01,1000);
    fit_pp->SetParameter(4,260);
    fit_pp->SetParLimits(4,60,300);
//fit_pp->FixParameter(4,100);
    fit_pp->SetParameter(5,460);
    fit_pp->SetParLimits(5,80,500);
//fit_pp->FixParameter(5,150);

//fit_pp->SetParameter(3,0.1);
//fit_pp->SetParLimits(3,0,2.0);
//fit_pp->SetParameter(4,0);
//fit_pp->SetParLimits(4,-1,1);

//fit_pp->FixParameter(3,0);
//fit_pp->FixParameter(4,0);
//fit_pp->FixParameter(5,0);

    fit_pp->SetParameter(6,1.0);
    fit_pp->SetParLimits(6,0.0,2.0);

    fit_pp->SetParameter(7,0);
    fit_pp->SetParLimits(7,-1e-2,1e-2);

    fit_pp->SetParameter(8,0);
    fit_pp->SetParLimits(8,-1e-4,1e-4);
//fit_pp->FixParameter(8,0);

    fit_pp->SetParameter(9,0);
    fit_pp->SetParLimits(9,-1e-6,1e-6);
fit_pp->FixParameter(9,0);


fit_pp->FixParameter(0,2.23021e-01);
fit_pp->FixParameter(1,2.93243e-01);
fit_pp->FixParameter(2,2.00000e+00);
fit_pp->FixParameter(3,2.76903e+01);
fit_pp->FixParameter(4,9.76146e+01);
fit_pp->FixParameter(5,1.50618e+02);
fit_pp->FixParameter(6,1.19729e+00);
fit_pp->FixParameter(7,-3.13238e-05);
fit_pp->FixParameter(8,2.03317e-08);
fit_pp->FixParameter(9,0.00000e+00);




    hWHICHONE->Fit(fit_pp,"S, N, R, M");
//histo_McCk[0]->Fit(fit_pp,"S, N, R, M");
//histo_McCk[4]->Fit(fit_pp,"S, N, R, M");

    TF1* fitNonjet_pp = new TF1("fitNonjet_pp","[0]*(1.+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x)",0,4500);
    fitNonjet_pp->SetParameter(0,fit_pp->GetParameter(6));
    fitNonjet_pp->SetParameter(1,fit_pp->GetParameter(7));
    fitNonjet_pp->SetParameter(2,fit_pp->GetParameter(8));
    fitNonjet_pp->SetParameter(3,fit_pp->GetParameter(9));
    fitNonjet_pp->SetParameter(4,0);

    TF1* fitJet_pp = new TF1("fitJet_pp","1.+[0]*TMath::Gaus(x,[1],[2],1);",0,4500);
    fitJet_pp->SetParameter(0,fit_pp->GetParameter(3));
    fitJet_pp->SetParameter(1,fit_pp->GetParameter(4));
    fitJet_pp->SetParameter(2,fit_pp->GetParameter(5));


    TF1* fitQS_pp = new TF1("fitQS_pp","1-[0]*exp(-pow(x*[1]/197.327,[2]))",0,4500);
    fitQS_pp->SetParameter(0,fit_pp->GetParameter(0));
    fitQS_pp->SetParameter(1,fit_pp->GetParameter(1));
    fitQS_pp->SetParameter(2,fit_pp->GetParameter(2));

    TF1* fit_wo_Jet_pp = new TF1("fit_wo_Jet_pp","(1-[0]*exp(-pow(x*[1]/197.327,[2])))*[3]*(1.+[4]*x+[5]*x*x+[6]*x*x*x+[7]*x*x*x*x)",0,4500);
    fit_wo_Jet_pp->SetParameter(0,fit_pp->GetParameter(0));
    fit_wo_Jet_pp->SetParameter(1,fit_pp->GetParameter(1));
    fit_wo_Jet_pp->SetParameter(2,fit_pp->GetParameter(2));
    fit_wo_Jet_pp->SetParameter(3,fit_pp->GetParameter(6));
    fit_wo_Jet_pp->SetParameter(4,fit_pp->GetParameter(7));
    fit_wo_Jet_pp->SetParameter(5,fit_pp->GetParameter(8));
    fit_wo_Jet_pp->SetParameter(6,fit_pp->GetParameter(9));
    fit_wo_Jet_pp->SetParameter(7,0);

    //TF1* fitFull_pp = new TF1("fitFull_pp",dimi_fractionQS_fitter,0,4500,10);
    TF1* fitFull_pp = new TF1("fitFull_pp",dimi_fractionPolQS_fitter,0,4500,10);

    fitFull_pp->SetParameter(0,fit_pp->GetParameter(0));
    fitFull_pp->SetParameter(1,fit_pp->GetParameter(1));
    fitFull_pp->SetParameter(2,fit_pp->GetParameter(2));
    fitFull_pp->SetParameter(3,fit_pp->GetParameter(3));
    fitFull_pp->SetParameter(4,fit_pp->GetParameter(4));
    fitFull_pp->SetParameter(5,fit_pp->GetParameter(5));
    fitFull_pp->SetParameter(6,fit_pp->GetParameter(6));
    fitFull_pp->SetParameter(7,fit_pp->GetParameter(7));
    fitFull_pp->SetParameter(8,fit_pp->GetParameter(8));
    fitFull_pp->SetParameter(9,fit_pp->GetParameter(9));
    fitFull_pp->SetNpx(1024);
/*
    fit_pp->Write();
    fitNonjet_pp->Write();
    fitJet_pp->Write();
    fitQS_pp->Write();
    fit_wo_Jet_pp->Write();
    fitFull_pp->Write();

    TH1F* hOriginal_pp = new TH1F("hOriginal_pp","hOriginal_pp",hWHICHONE->GetNbinsX(),
                                    hWHICHONE->GetBinLowEdge(1),hWHICHONE->GetXaxis()->GetBinUpEdge(hWHICHONE->GetNbinsX()));
    for(unsigned uBin=1; uBin<=hWHICHONE->GetNbinsX(); uBin++){
        double MOM = hOriginal_pp->GetBinCenter(uBin);
        if(fitFull_pp->Eval(MOM)<1e-8){
            hOriginal_pp->SetBinContent(uBin,100);
            hOriginal_pp->SetBinError(uBin,100);
        }
        else{
            //hCorrected_pp->SetBinContent(uBin,hWHICHONE->GetBinContent(uBin)/fitFull_pp->Eval(MOM));
            //hCorrected_pp->SetBinError(uBin,hWHICHONE->GetBinError(uBin)/fitFull_pp->Eval(MOM));
            hOriginal_pp->SetBinContent(uBin,dWHICHONE.GetBinContent(uBin-1));
            hOriginal_pp->SetBinError(uBin,dWHICHONE.GetBinError(uBin-1));
            //printf("%f %e (%e)\n",MOM,hWHICHONE->GetBinContent(uBin),fitFull_pp->Eval(MOM));
        }
    }
    hOriginal_pp->Write();

    TH1F* hCorrected_pp = new TH1F("hCorrected_pp","hCorrected_pp",hWHICHONE->GetNbinsX(),
                                    hWHICHONE->GetBinLowEdge(1),hWHICHONE->GetXaxis()->GetBinUpEdge(hWHICHONE->GetNbinsX()));
    for(unsigned uBin=1; uBin<=hWHICHONE->GetNbinsX(); uBin++){
        double MOM = hCorrected_pp->GetBinCenter(uBin);
        if(fitFull_pp->Eval(MOM)<1e-8){
            hCorrected_pp->SetBinContent(uBin,100);
            hCorrected_pp->SetBinError(uBin,100);
        }
        else{
            //hCorrected_pp->SetBinContent(uBin,hWHICHONE->GetBinContent(uBin)/fitFull_pp->Eval(MOM));
            //hCorrected_pp->SetBinError(uBin,hWHICHONE->GetBinError(uBin)/fitFull_pp->Eval(MOM));
            hCorrected_pp->SetBinContent(uBin,dWHICHONE.GetBinContent(uBin-1)/fitFull_pp->Eval(MOM));
            hCorrected_pp->SetBinError(uBin,dWHICHONE.GetBinError(uBin-1)/fitFull_pp->Eval(MOM));
            //printf("%f %e (%e)\n",MOM,hWHICHONE->GetBinContent(uBin),fitFull_pp->Eval(MOM));
        }
    }
    hCorrected_pp->Write();
*/
    TCanvas* c_pp = new TCanvas("c_pp", "c_pp", 1);
    c_pp->cd(0);
    c_pp->SetCanvasSize(1280, 720);
    c_pp->SetMargin(0.15,0.05,0.2,0.05);//lrbt

    hWHICHONE->GetYaxis()->SetRangeUser(0.6,1.4);
    hWHICHONE->Draw();
    fitFull_pp->Draw("same");
    //fitPert_pp->Draw("same");
    c_pp->SaveAs(OutputFolder+"c_pp.png");

    ////////////////////////////////////////////////////////
}

void ParametrizeTemplates(const TString Particle1, const TString Particle2){

    const TString InputFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/OtherTasks/ALL_CorrectedMC_EXP/fOutput.root";
    const TString OutputFileName =
    "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/OtherTasks/ParametrizeTemplates/fTemplates_"+Particle1+Particle2+".root";

    TString SystemName = Particle1+Particle2;
    TString AntiSystemName = Particle1+"A"+Particle2;
    printf("SystemName = %s (%s)\n",SystemName.Data(),AntiSystemName.Data());

    TString McCk_PP = "histo_McCk_"+SystemName;
    TString McCk_PAP = "histo_McCk_"+AntiSystemName;
    TString DataCk_PP = "histo_DataCk_"+SystemName;
    TString DataCk_PAP = "histo_DataCk_"+AntiSystemName;

    TFile fInput(InputFileName,"read");

    TH1F* Template_PP;
    TH1F* Template_PAP;
    TH1F* Data_PP;
    TH1F* Data_PAP;


    Template_PP = (TH1F*)fInput.Get(McCk_PP);
    Template_PAP = (TH1F*)fInput.Get(McCk_PAP);
    Data_PP = (TH1F*)fInput.Get(DataCk_PP);
    Data_PAP = (TH1F*)fInput.Get(DataCk_PAP);

    TFile fOutput(OutputFileName,"recreate");

    TF1* fTemplate_PP = new TF1("fTemplate_PP","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x",500,3000);
    fTemplate_PP->SetParameter(0,1.03);
    fTemplate_PP->SetParameter(1,-1.21e-4);
    fTemplate_PP->SetParLimits(1,-1e-3,1e-3);
    fTemplate_PP->SetParameter(2,1.08e-7);
    fTemplate_PP->SetParLimits(2,-1e-6,1e-6);
    fTemplate_PP->SetParameter(3,-2.76e-11);
    fTemplate_PP->SetParLimits(3,-1e-9,1e-9);
    fTemplate_PP->SetParameter(4,2.76e-15);
    fTemplate_PP->SetParLimits(4,-1e-12,1e-12);
    Template_PP->Fit(fTemplate_PP,"S, N, R, M");





    Template_PP->Write();
    fTemplate_PP->Write();

    Template_PAP->Write();
    //fTemplate_PP->Write();

}

void ReweightME(){
    TString OutputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/OtherTasks/ReweightME/";
    const unsigned NumSpecies = 8;
    TString* SpeciesName = new TString[NumSpecies];
    TString* SpeciesLegend = new TString[NumSpecies];
    SpeciesName[0] = "pp"; SpeciesLegend[0] = "p-p";
    SpeciesName[1] = "pL"; SpeciesLegend[1] = "p-#Lambda";
    SpeciesName[2] = "pXi"; SpeciesLegend[2] = "p-#Xi";
    SpeciesName[3] = "LL"; SpeciesLegend[3] = "#Lambda-#Lambda";
    SpeciesName[4] = "pAp"; SpeciesLegend[4] = "p-#bar{p}";
    SpeciesName[5] = "pAL"; SpeciesLegend[5] = "p-#bar{#Lambda}";
    SpeciesName[6] = "pAXi"; SpeciesLegend[6] = "p-#bar{#Xi}";
    SpeciesName[7] = "LAL"; SpeciesLegend[7] = "#Lambda-#bar{#Lambda}";


    TString* DataFolder = new TString[NumSpecies];
    for(unsigned uSpec=0; uSpec<NumSpecies; uSpec++)
    DataFolder[uSpec] = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/Data/AnalysisResults/Nano/";

    TString* McFolder = new TString[NumSpecies];
    McFolder[0] = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/Fast_BBar/Trains_MCAOD/Norm018028/";
    McFolder[1] = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/Fast_BBar/Trains_MCAOD/Norm018028/";
    McFolder[2] = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/Fast_BBar/Trains_MCAOD/Norm018028/";
    McFolder[3] = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/Fast_BBar/Trains_MCAOD/Norm018028/";
    McFolder[4] = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/Fast_BBar/Trains_MCAOD/Norm018028/";
    McFolder[5] = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/Fast_BBar/Trains_MCAOD/Norm018028/";
    McFolder[6] = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/Fast_BBar/Trains_MCAOD/Norm018028/";
    McFolder[7] = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/Fast_BBar/Trains_MCAOD/Norm018028/";

    TString* DataFile = new TString[NumSpecies];
    //DataFile[0] = "CFOutput_pp_8.root";
    DataFile[0] = "AnalysisResults.root";
    DataFile[1] = "AnalysisResults.root";
    DataFile[2] = "AnalysisResults.root";
    DataFile[3] = "AnalysisResults.root";
    DataFile[4] = "AnalysisResults.root";
    DataFile[5] = "AnalysisResults.root";
    DataFile[6] = "AnalysisResults.root";
    DataFile[7] = "AnalysisResults.root";

    TString* McFile = new TString[NumSpecies];
    McFile[0] = "CFOutput_pp_8.root";
    McFile[1] = "CFOutput_pL_8.root";
    McFile[2] = "CFOutput_pXi_8.root";
    McFile[3] = "CFOutput_LL_8.root";
    McFile[4] = "CFOutput_pAp_8.root";
    McFile[5] = "CFOutput_pAL_8.root";
    McFile[6] = "CFOutput_pAXi_8.root";
    McFile[7] = "CFOutput_LAL_8.root";

    TString* DataHisto = new TString[NumSpecies];
    DataHisto[0] = "hCk_ReweightedMeV_0";
    DataHisto[1] = "hCk_ReweightedMeV_0";
    DataHisto[2] = "hCk_ReweightedMeV_0";
    DataHisto[3] = "hCk_ReweightedMeV_0";
    DataHisto[4] = "hCk_ReweightedMeV_0";
    DataHisto[5] = "hCk_ReweightedMeV_0";
    DataHisto[6] = "hCk_ReweightedMeV_0";
    DataHisto[7] = "hCk_ReweightedMeV_0";

    TString* McHisto = new TString[NumSpecies];
    McHisto[0] = "hCk_ReweightedMeV_0";
    McHisto[1] = "hCk_ReweightedMeV_0";
    McHisto[2] = "hCk_ReweightedMeV_0";
    McHisto[3] = "hCk_ReweightedMeV_0";
    McHisto[4] = "hCk_ReweightedMeV_0";
    McHisto[5] = "hCk_ReweightedMeV_0";
    McHisto[6] = "hCk_ReweightedMeV_0";
    McHisto[7] = "hCk_ReweightedMeV_0";

    int* RescaleFactor = new int[NumSpecies];
    RescaleFactor[0] = 1;
    RescaleFactor[1] = 1;
    RescaleFactor[2] = 1;
    RescaleFactor[3] = 1;
    RescaleFactor[4] = 1;
    RescaleFactor[5] = 1;
    RescaleFactor[6] = 1;
    RescaleFactor[7] = 1;

    int* RescaleSeMeAxis = new int[NumSpecies];
    RescaleSeMeAxis[0] = 1000;
    RescaleSeMeAxis[1] = 1000;
    RescaleSeMeAxis[2] = 1000;
    RescaleSeMeAxis[3] = 1000;
    RescaleSeMeAxis[4] = 1000;
    RescaleSeMeAxis[5] = 1000;
    RescaleSeMeAxis[6] = 1000;
    RescaleSeMeAxis[7] = 1000;

    TString* DataResultsList = new TString[NumSpecies];
    for(unsigned uSpec=0; uSpec<NumSpecies; uSpec++) DataResultsList[uSpec] = "HMBBarResults8";

    TString* DataPP_List = new TString[NumSpecies];
    DataPP_List[0] = "Particle0_Particle0";
    DataPP_List[1] = "Particle0_Particle2";
    DataPP_List[2] = "Particle0_Particle4";
    DataPP_List[3] = "Particle2_Particle2";
    DataPP_List[4] = "Particle0_Particle1";
    DataPP_List[5] = "Particle0_Particle3";
    DataPP_List[6] = "Particle0_Particle5";
    DataPP_List[7] = "Particle2_Particle3";

    TString* DataAPAP_List = new TString[NumSpecies];
    DataAPAP_List[0] = "Particle1_Particle1";
    DataAPAP_List[1] = "Particle1_Particle3";
    DataAPAP_List[2] = "Particle1_Particle5";
    DataAPAP_List[3] = "Particle3_Particle3";
    DataAPAP_List[4] = "Particle1_Particle1";
    DataAPAP_List[5] = "Particle1_Particle3";
    DataAPAP_List[6] = "Particle3_Particle3";
    DataAPAP_List[7] = "Particle1_Particle3";

    TString* DataPP_SEMultDist = new TString[NumSpecies];
    for(unsigned uSpec=0; uSpec<NumSpecies; uSpec++)DataPP_SEMultDist[uSpec] = "SEMultDist_"+DataPP_List[uSpec];

    TString* DataAPAP_SEMultDist = new TString[NumSpecies];
    for(unsigned uSpec=0; uSpec<NumSpecies; uSpec++)DataAPAP_SEMultDist[uSpec] = "SEMultDist_"+DataAPAP_List[uSpec];

    TString* DataPP_MEMultDist = new TString[NumSpecies];
    for(unsigned uSpec=0; uSpec<NumSpecies; uSpec++)DataPP_MEMultDist[uSpec] = "MEMultDist_"+DataPP_List[uSpec];

    TString* DataAPAP_MEMultDist = new TString[NumSpecies];
    for(unsigned uSpec=0; uSpec<NumSpecies; uSpec++)DataAPAP_MEMultDist[uSpec] = "MEMultDist_"+DataAPAP_List[uSpec];

    DLM_Histo<float>* dlm_SEMult_PP = new DLM_Histo<float>[NumSpecies];
    DLM_Histo<float>* dlm_SEMult_APAP = new DLM_Histo<float>[NumSpecies];
    DLM_Histo<float>* dlm_SEMult = new DLM_Histo<float>[NumSpecies];

    DLM_Histo<float>* dlm_MEMult_PP = new DLM_Histo<float>[NumSpecies];
    DLM_Histo<float>* dlm_MEMult_APAP = new DLM_Histo<float>[NumSpecies];
    DLM_Histo<float>* dlm_MEMult = new DLM_Histo<float>[NumSpecies];

    DLM_Histo<float>* dlm_MEMultW_PP = new DLM_Histo<float>[NumSpecies];
    DLM_Histo<float>* dlm_MEMultW_APAP = new DLM_Histo<float>[NumSpecies];
    DLM_Histo<float>* dlm_MEMultW = new DLM_Histo<float>[NumSpecies];

    DLM_Histo<float>* dlm_CkMult_PP = new DLM_Histo<float>[NumSpecies];
    DLM_Histo<float>* dlm_CkMult_APAP = new DLM_Histo<float>[NumSpecies];
    DLM_Histo<float>* dlm_CkMult = new DLM_Histo<float>[NumSpecies];

    DLM_Histo<float>* dlm_CkMultW_PP = new DLM_Histo<float>[NumSpecies];
    DLM_Histo<float>* dlm_CkMultW_APAP = new DLM_Histo<float>[NumSpecies];
    DLM_Histo<float>* dlm_CkMultW = new DLM_Histo<float>[NumSpecies];

    TH2F** h_CkMult = new TH2F* [NumSpecies];
    TH2F** h_CkMultW = new TH2F* [NumSpecies];

    for(unsigned uSpec=0; uSpec<NumSpecies; uSpec++){
        unsigned NumBins;
        double kMin;
        double kMax;
        TH2F* hSEpp_tmp;
        TH2F* hMEpp_tmp;
        TH2F* hSEapap_tmp;
        TH2F* hMEapap_tmp;
        TDirectoryFile *dirResults=NULL;
        TList* listresults_tmp=NULL;
        TList* listpp_tmp=NULL;
        TList* listapap_tmp=NULL;

printf("Opening file %s%s\n",DataFolder[uSpec].Data(),DataFile[uSpec].Data());

        TFile datafile(DataFolder[uSpec]+DataFile[uSpec],"read");
        dirResults=(TDirectoryFile*)(datafile.FindObjectAny(DataResultsList[uSpec]));
printf("dirResults=%p\n",dirResults);
        dirResults->GetObject(DataResultsList[uSpec],listresults_tmp);
printf("listresults_tmp=%p\n",listresults_tmp);
        listpp_tmp = (TList*)listresults_tmp->FindObject(DataPP_List[uSpec]);
printf("listpp_tmp=%p\n",listpp_tmp);
        if(DataAPAP_List[uSpec]!="")listapap_tmp = (TList*)listresults_tmp->FindObject(DataAPAP_List[uSpec]);
        else listapap_tmp=NULL;
printf("listapap_tmp=%p\n",listapap_tmp);
        hSEpp_tmp = (TH2F*)listpp_tmp->FindObject(DataPP_SEMultDist[uSpec]);
printf("hSEpp_tmp=%p\n",hSEpp_tmp);
        if(listapap_tmp)hSEapap_tmp = (TH2F*)listapap_tmp->FindObject(DataAPAP_SEMultDist[uSpec]);
        else hSEapap_tmp=NULL;
printf("hSEapap_tmp=%p\n",hSEapap_tmp);
        hMEpp_tmp = (TH2F*)listpp_tmp->FindObject(DataPP_MEMultDist[uSpec]);
printf("hMEpp_tmp=%p\n",hMEpp_tmp);
        if(listapap_tmp)hMEapap_tmp = (TH2F*)listapap_tmp->FindObject(DataAPAP_MEMultDist[uSpec]);
        else hMEapap_tmp=NULL;
printf("hMEapap_tmp=%p\n",hMEapap_tmp);

        dlm_SEMult_PP[uSpec].SetUp(2);
        dlm_SEMult_PP[uSpec].SetUp(0,hSEpp_tmp->GetNbinsX(),
                                   hSEpp_tmp->GetXaxis()->GetBinLowEdge(1)*RescaleSeMeAxis[uSpec],
                                   hSEpp_tmp->GetXaxis()->GetBinUpEdge(hSEpp_tmp->GetNbinsX())*RescaleSeMeAxis[uSpec]);
        dlm_SEMult_PP[uSpec].SetUp(1,hSEpp_tmp->GetNbinsY(),
                                   hSEpp_tmp->GetYaxis()->GetBinLowEdge(1),
                                   hSEpp_tmp->GetYaxis()->GetBinUpEdge(hSEpp_tmp->GetNbinsY()));
        dlm_SEMult_PP[uSpec].Initialize();
        for(unsigned uBinX=0; uBinX<hSEpp_tmp->GetNbinsX(); uBinX++){
            for(unsigned uBinY=0; uBinY<hSEpp_tmp->GetNbinsY(); uBinY++){
                unsigned WhichBin[2]; WhichBin[0]=uBinX; WhichBin[1]=uBinY;
                dlm_SEMult_PP[uSpec].SetBinContent(WhichBin,hSEpp_tmp->GetBinContent(uBinX+1,uBinY+1));
                dlm_SEMult_PP[uSpec].SetBinError(WhichBin,hSEpp_tmp->GetBinError(uBinX+1,uBinY+1));
            }
        }

        dlm_SEMult_APAP[uSpec].SetUp(2);
        dlm_SEMult_APAP[uSpec].SetUp(0,hSEapap_tmp->GetNbinsX(),
                                   hSEapap_tmp->GetXaxis()->GetBinLowEdge(1)*RescaleSeMeAxis[uSpec],
                                   hSEapap_tmp->GetXaxis()->GetBinUpEdge(hSEapap_tmp->GetNbinsX())*RescaleSeMeAxis[uSpec]);
        dlm_SEMult_APAP[uSpec].SetUp(1,hSEapap_tmp->GetNbinsY(),
                                   hSEapap_tmp->GetYaxis()->GetBinLowEdge(1),
                                   hSEapap_tmp->GetYaxis()->GetBinUpEdge(hSEapap_tmp->GetNbinsY()));
        dlm_SEMult_APAP[uSpec].Initialize();
        for(unsigned uBinX=0; uBinX<hSEapap_tmp->GetNbinsX(); uBinX++){
            for(unsigned uBinY=0; uBinY<hSEapap_tmp->GetNbinsY(); uBinY++){
                unsigned WhichBin[2]; WhichBin[0]=uBinX; WhichBin[1]=uBinY;
                dlm_SEMult_APAP[uSpec].SetBinContent(WhichBin,hSEapap_tmp->GetBinContent(uBinX+1,uBinY+1));
                dlm_SEMult_APAP[uSpec].SetBinError(WhichBin,hSEapap_tmp->GetBinError(uBinX+1,uBinY+1));
            }
        }

        dlm_MEMult_PP[uSpec].SetUp(2);
        dlm_MEMult_PP[uSpec].SetUp(0,hMEpp_tmp->GetNbinsX(),
                                   hMEpp_tmp->GetXaxis()->GetBinLowEdge(1)*RescaleSeMeAxis[uSpec],
                                   hMEpp_tmp->GetXaxis()->GetBinUpEdge(hMEpp_tmp->GetNbinsX())*RescaleSeMeAxis[uSpec]);
        dlm_MEMult_PP[uSpec].SetUp(1,hMEpp_tmp->GetNbinsY(),
                                   hMEpp_tmp->GetYaxis()->GetBinLowEdge(1),
                                   hMEpp_tmp->GetYaxis()->GetBinUpEdge(hMEpp_tmp->GetNbinsY()));
        dlm_MEMult_PP[uSpec].Initialize();
        for(unsigned uBinX=0; uBinX<hMEpp_tmp->GetNbinsX(); uBinX++){
            for(unsigned uBinY=0; uBinY<hMEpp_tmp->GetNbinsY(); uBinY++){
                unsigned WhichBin[2]; WhichBin[0]=uBinX; WhichBin[1]=uBinY;
                dlm_MEMult_PP[uSpec].SetBinContent(WhichBin,hMEpp_tmp->GetBinContent(uBinX+1,uBinY+1));
                dlm_MEMult_PP[uSpec].SetBinError(WhichBin,hMEpp_tmp->GetBinError(uBinX+1,uBinY+1));
            }
        }

        dlm_MEMult_APAP[uSpec].SetUp(2);
        dlm_MEMult_APAP[uSpec].SetUp(0,hMEapap_tmp->GetNbinsX(),
                                   hMEapap_tmp->GetXaxis()->GetBinLowEdge(1)*RescaleSeMeAxis[uSpec],
                                   hMEapap_tmp->GetXaxis()->GetBinUpEdge(hMEapap_tmp->GetNbinsX())*RescaleSeMeAxis[uSpec]);
        dlm_MEMult_APAP[uSpec].SetUp(1,hMEapap_tmp->GetNbinsY(),
                                   hMEapap_tmp->GetYaxis()->GetBinLowEdge(1),
                                   hMEapap_tmp->GetYaxis()->GetBinUpEdge(hMEapap_tmp->GetNbinsY()));
        dlm_MEMult_APAP[uSpec].Initialize();
        for(unsigned uBinX=0; uBinX<hMEapap_tmp->GetNbinsX(); uBinX++){
            for(unsigned uBinY=0; uBinY<hMEapap_tmp->GetNbinsY(); uBinY++){
                unsigned WhichBin[2]; WhichBin[0]=uBinX; WhichBin[1]=uBinY;
                dlm_MEMult_APAP[uSpec].SetBinContent(WhichBin,hMEapap_tmp->GetBinContent(uBinX+1,uBinY+1));
                dlm_MEMult_APAP[uSpec].SetBinError(WhichBin,hMEapap_tmp->GetBinError(uBinX+1,uBinY+1));
            }
        }

        dlm_CkMult_PP[uSpec] = dlm_SEMult_PP[uSpec]; dlm_CkMult_PP[uSpec] /= dlm_MEMult_PP[uSpec];
        dlm_CkMult_APAP[uSpec] = dlm_SEMult_APAP[uSpec]; dlm_CkMult_APAP[uSpec] /= dlm_MEMult_APAP[uSpec];
        dlm_CkMult[uSpec] = dlm_CkMult_PP[uSpec]; dlm_CkMult[uSpec].AddWeightedHisto(dlm_CkMult_APAP[uSpec]);




        DLM_Histo<float> MultSE_PP_kIntegrated;
        MultSE_PP_kIntegrated.SetUp(1);
        MultSE_PP_kIntegrated.SetUp(0,dlm_SEMult_PP[uSpec].GetNbins(1),dlm_SEMult_PP[uSpec].GetLowEdge(1),dlm_SEMult_PP[uSpec].GetUpEdge(1));
        MultSE_PP_kIntegrated.Initialize();
        double kRange[2];kRange[0]=dlm_SEMult_PP[uSpec].GetLowEdge(0);kRange[1]=dlm_SEMult_PP[uSpec].GetUpEdge(0);
        for(unsigned uBin=0; uBin<MultSE_PP_kIntegrated.GetNbins(); uBin++){
            double MultRange[2];MultRange[0]=MultSE_PP_kIntegrated.GetBinLowEdge(0,uBin);MultRange[1]=MultSE_PP_kIntegrated.GetBinUpEdge(0,uBin);
            MultSE_PP_kIntegrated.SetBinContent(uBin,dlm_SEMult_PP[uSpec].Integral(kRange,MultRange));
        }

        DLM_Histo<float> MultME_PP_kIntegrated;
        MultME_PP_kIntegrated.SetUp(1);
        MultME_PP_kIntegrated.SetUp(0,dlm_MEMult_PP[uSpec].GetNbins(1),dlm_MEMult_PP[uSpec].GetLowEdge(1),dlm_MEMult_PP[uSpec].GetUpEdge(1));
        MultME_PP_kIntegrated.Initialize();
        kRange[0]=dlm_MEMult_PP[uSpec].GetLowEdge(0); kRange[1]=dlm_MEMult_PP[uSpec].GetUpEdge(0);
        for(unsigned uBin=0; uBin<MultME_PP_kIntegrated.GetNbins(); uBin++){
            double MultRange[2];MultRange[0]=MultME_PP_kIntegrated.GetBinLowEdge(0,uBin);MultRange[1]=MultME_PP_kIntegrated.GetBinUpEdge(0,uBin);
            MultME_PP_kIntegrated.SetBinContent(uBin,dlm_MEMult_PP[uSpec].Integral(kRange,MultRange));
        }
        dlm_MEMultW_PP[uSpec] = dlm_MEMult_PP[uSpec];
        for(unsigned uBinX=0; uBinX<dlm_MEMultW_PP[uSpec].GetNbins(0); uBinX++){
            for(unsigned uBinY=0; uBinY<dlm_MEMultW_PP[uSpec].GetNbins(1); uBinY++){
                unsigned WhichBin[2]; WhichBin[0]=uBinX; WhichBin[1]=uBinY;
                double Rescale = MultSE_PP_kIntegrated.GetBinContent(uBinY)/MultME_PP_kIntegrated.GetBinContent(uBinY);
                dlm_MEMultW_PP[uSpec].SetBinContent(WhichBin,dlm_MEMultW_PP[uSpec].GetBinContent(WhichBin)*Rescale);
                dlm_MEMultW_PP[uSpec].SetBinError(WhichBin,dlm_MEMultW_PP[uSpec].GetBinError(WhichBin)*Rescale);
            }
        }
        dlm_CkMultW_PP[uSpec] = dlm_SEMult_PP[uSpec];
        dlm_CkMultW_PP[uSpec] /= dlm_MEMultW_PP[uSpec];



        DLM_Histo<float> MultSE_APAP_kIntegrated;
        MultSE_APAP_kIntegrated.SetUp(1);
        MultSE_APAP_kIntegrated.SetUp(0,dlm_SEMult_APAP[uSpec].GetNbins(1),dlm_SEMult_APAP[uSpec].GetLowEdge(1),dlm_SEMult_APAP[uSpec].GetUpEdge(1));
        MultSE_APAP_kIntegrated.Initialize();
        kRange[2];kRange[0]=dlm_SEMult_APAP[uSpec].GetLowEdge(0);kRange[1]=dlm_SEMult_APAP[uSpec].GetUpEdge(0);
        for(unsigned uBin=0; uBin<MultSE_APAP_kIntegrated.GetNbins(); uBin++){
            double MultRange[2];MultRange[0]=MultSE_APAP_kIntegrated.GetBinLowEdge(0,uBin);MultRange[1]=MultSE_APAP_kIntegrated.GetBinUpEdge(0,uBin);
            MultSE_APAP_kIntegrated.SetBinContent(uBin,dlm_SEMult_APAP[uSpec].Integral(kRange,MultRange));
        }
        DLM_Histo<float> MultME_APAP_kIntegrated;
        MultME_APAP_kIntegrated.SetUp(1);
        MultME_APAP_kIntegrated.SetUp(0,dlm_MEMult_APAP[uSpec].GetNbins(1),dlm_MEMult_APAP[uSpec].GetLowEdge(1),dlm_MEMult_APAP[uSpec].GetUpEdge(1));
        MultME_APAP_kIntegrated.Initialize();
        kRange[0]=dlm_MEMult_APAP[uSpec].GetLowEdge(0); kRange[1]=dlm_MEMult_APAP[uSpec].GetUpEdge(0);
        for(unsigned uBin=0; uBin<MultME_APAP_kIntegrated.GetNbins(); uBin++){
            double MultRange[2];MultRange[0]=MultME_APAP_kIntegrated.GetBinLowEdge(0,uBin);MultRange[1]=MultME_APAP_kIntegrated.GetBinUpEdge(0,uBin);
            MultME_APAP_kIntegrated.SetBinContent(uBin,dlm_MEMult_APAP[uSpec].Integral(kRange,MultRange));
        }
        dlm_MEMultW_APAP[uSpec] = dlm_MEMult_APAP[uSpec];
        for(unsigned uBinX=0; uBinX<dlm_MEMultW_APAP[uSpec].GetNbins(0); uBinX++){
            for(unsigned uBinY=0; uBinY<dlm_MEMultW_APAP[uSpec].GetNbins(1); uBinY++){
                unsigned WhichBin[2]; WhichBin[0]=uBinX; WhichBin[1]=uBinY;
                double Rescale = MultSE_APAP_kIntegrated.GetBinContent(uBinY)/MultME_APAP_kIntegrated.GetBinContent(uBinY);
                dlm_MEMultW_APAP[uSpec].SetBinContent(WhichBin,dlm_MEMultW_APAP[uSpec].GetBinContent(WhichBin)*Rescale);
                dlm_MEMultW_APAP[uSpec].SetBinError(WhichBin,dlm_MEMultW_APAP[uSpec].GetBinError(WhichBin)*Rescale);
            }
        }
        dlm_CkMultW_APAP[uSpec] = dlm_SEMult_APAP[uSpec];
        dlm_CkMultW_APAP[uSpec] /= dlm_MEMultW_APAP[uSpec];

        dlm_CkMultW[uSpec] = dlm_CkMultW_PP[uSpec]; dlm_CkMultW[uSpec].AddWeightedHisto(dlm_CkMultW_APAP[uSpec]);
    }

    TFile fOutput(OutputFolder+"fReweightME.root","recreate");
    for(unsigned uSpec=0; uSpec<NumSpecies; uSpec++){
        h_CkMult[uSpec] = new TH2F("h_CkMult_"+SpeciesName[uSpec],"h_CkMult_"+SpeciesName[uSpec],
                                   dlm_CkMult[uSpec].GetNbins(0),dlm_CkMult[uSpec].GetLowEdge(0),dlm_CkMult[uSpec].GetUpEdge(0),
                                   dlm_CkMult[uSpec].GetNbins(1),dlm_CkMult[uSpec].GetLowEdge(1),dlm_CkMult[uSpec].GetUpEdge(1));
        for(unsigned uBinX=0; uBinX<dlm_CkMult[uSpec].GetNbins(0); uBinX++){
            for(unsigned uBinY=0; uBinY<dlm_CkMult[uSpec].GetNbins(1); uBinY++){
                unsigned WhichBin[2]; WhichBin[0]=uBinX; WhichBin[1]=uBinY;
                h_CkMult[uSpec]->SetBinContent(uBinX+1,uBinY+1,dlm_CkMult[uSpec].GetBinContent(WhichBin));
                h_CkMult[uSpec]->SetBinError(uBinX+1,uBinY+1,dlm_CkMult[uSpec].GetBinError(WhichBin));
                if(dlm_CkMult[uSpec].GetBinContent(WhichBin)>100||dlm_CkMult[uSpec].GetBinContent(WhichBin)!=dlm_CkMult[uSpec].GetBinContent(WhichBin)){
                    h_CkMult[uSpec]->SetBinContent(uBinX+1,uBinY+1,100);
                    h_CkMult[uSpec]->SetBinError(uBinX+1,uBinY+1,100);
                }
            }
        }

        h_CkMultW[uSpec] = new TH2F("h_CkMultW_"+SpeciesName[uSpec],"h_CkMultW_"+SpeciesName[uSpec],
                                   dlm_CkMultW[uSpec].GetNbins(0),dlm_CkMultW[uSpec].GetLowEdge(0),dlm_CkMultW[uSpec].GetUpEdge(0),
                                   dlm_CkMultW[uSpec].GetNbins(1),dlm_CkMultW[uSpec].GetLowEdge(1),dlm_CkMultW[uSpec].GetUpEdge(1));
        for(unsigned uBinX=0; uBinX<dlm_CkMultW[uSpec].GetNbins(0); uBinX++){
            for(unsigned uBinY=0; uBinY<dlm_CkMultW[uSpec].GetNbins(1); uBinY++){
                unsigned WhichBin[2]; WhichBin[0]=uBinX; WhichBin[1]=uBinY;
                h_CkMultW[uSpec]->SetBinContent(uBinX+1,uBinY+1,dlm_CkMultW[uSpec].GetBinContent(WhichBin));
                h_CkMultW[uSpec]->SetBinError(uBinX+1,uBinY+1,dlm_CkMultW[uSpec].GetBinError(WhichBin));
                if(dlm_CkMultW[uSpec].GetBinContent(WhichBin)>100||dlm_CkMultW[uSpec].GetBinContent(WhichBin)!=dlm_CkMultW[uSpec].GetBinContent(WhichBin)){
                    h_CkMultW[uSpec]->SetBinContent(uBinX+1,uBinY+1,100);
                    h_CkMultW[uSpec]->SetBinError(uBinX+1,uBinY+1,100);
                }
            }
        }

        h_CkMult[uSpec]->Write();
        h_CkMultW[uSpec]->Write();
    }


}



void BabyToyJet(){

    DLM_Random dRanGen(11);
    const unsigned NumIter = 100000;

    TH1F* hBoltzmann = new TH1F("hBoltzmann","hBoltzmann",256,0,16);
    TH1F* hMom_PP = new TH1F("hMom_PP","hMom_PP",256,0,1);
    TH1F* hMom_APAP = new TH1F("hMom_APAP","hMom_APAP",256,0,1);

    double Energy;
    const double MASS = 0.938;
    double masses[4] = { MASS,MASS,MASS,MASS } ;
    double Prob4decay;
    TGenPhaseSpace event;

    for(unsigned uIter=0; uIter<NumIter; uIter++){
        Energy = dRanGen.StableR(3,2.0,0,1.5,0);
        if(Energy<2.*MASS) continue;
        if(Energy<4.*MASS) {Prob4decay = 0;}
        else{
            Prob4decay = (Energy-4.*MASS)/(Energy-6.*MASS);
        }
        TLorentzVector String(0,0,Energy,Energy);
        int NumPart;
        if(dRanGen.Uniform(0,1)>Prob4decay) {NumPart=2;}
        else {NumPart=4;}
        event.SetDecay(String, NumPart, masses);
        event.Generate();
        //TLorentzVector** Decay = new TLorentzVector [NumPart];
        for(int iPart=0; iPart<NumPart; iPart++){
            //Decay[iPart] = event.GetDecay(iPart);
        }

    }




}


DLM_CkDecomposition* MM_PL;
CATS* MM_CatPL;
double MickeyFitter_pL(double* x, double* par){
    double& MOM = *x;
//MM_CatPL->SetChannelWeight(2,0.25*par[11]);
//MM_CatPL->SetChannelWeight(3,0.75*par[11]);
    MM_CatPL->SetAnaSource(0,par[12],true);
    if(MM_CatPL->GetNumSourcePars()>1){
        MM_CatPL->SetAnaSource(1,par[13],true);
    }
    MM_CatPL->KillTheCat();
    MM_PL->GetCk()->SetCutOff(320,par[10]);
    //MM_PL->GetCk()->SetSourcePar(0,par[12]);
    MM_PL->Update(true);
    static int COUNTER = 0;
    if(MOM<MM_PL->GetCk()->GetBinUpEdge(0,0)){
    COUNTER++;
    //if(COUNTER%100==0)
    printf("COUNTER=%i\n",COUNTER);
    }

    double CkVal;
    //if(MOM<320) CkVal = MM_PL->EvalCk(MOM);
    //else{
    //    CkVal = (MOM-320-MM_PL->EvalCk(320)*(MOM-par[10]))/(par[10]-320);
    //}
    CkVal = MM_PL->EvalCk(MOM);
    if(MOM<280) CkVal *= par[14];
    double Pauli = 1.-par[0]*exp(-pow(MOM*par[1]/197.327,par[2]));
    double Jet = 1.+par[3]*TMath::Gaus(MOM,par[4],par[5],0);
    double NonJet = par[6]*(1.+par[7]*MOM+par[8]*MOM*MOM+par[9]*MOM*MOM*MOM);
    return CkVal*NonJet*Jet*Pauli;
}

//including QS
void Fit_pL_MickeyMouse(){
    double ResidualSourceSize=1.5;

    double* MomBins_pL = NULL;
    double* FitRegion_pL = NULL;
    const unsigned NumBinsCk = 36;
    const double BinWidthCk = 12;
    const double MaxBinValCk = double(NumBinsCk)*BinWidthCk;
    unsigned NumMomBins_pL;
    //TString DataSample = "pp13TeV_HM_Dec19";
    TString DataSample = "pp13TeV_HM_RotPhiDec19";
    DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder("/home/dmihaylov/CernBox/CatsFiles");
    AnalysisObject.SetUpBinning_pL(DataSample,NumMomBins_pL,MomBins_pL,FitRegion_pL,0,0);

    DLM_Ck* Ck_pSigma0 = new DLM_Ck(1,0,NumMomBins_pL,MomBins_pL,Lednicky_gauss_Sigma0);
    Ck_pSigma0->SetSourcePar(0,ResidualSourceSize);

    TH2F* hResolution_pL = AnalysisObject.GetResolutionMatrix(DataSample,"pLambda");
    TH2F* hResidual_pL_pSigma0 = AnalysisObject.GetResidualMatrix("pLambda","pSigma0");
    TH2F* hResidual_pL_pXim = AnalysisObject.GetResidualMatrix("pLambda","pXim");

    TH1F* hData_pL = AnalysisObject.GetAliceExpCorrFun(DataSample,"pLambda","_0",1,false,-1);

    double lam_pL[5];
    double lam_pXim[5];

    AnalysisObject.SetUpLambdaPars_pL(DataSample,0,0,lam_pL);
    AnalysisObject.SetUpLambdaPars_pXim(DataSample,0,0,lam_pXim);

    //TString SourceDescription = "Gauss";
    TString SourceDescription = "McGauss_ResoTM";
    //TString SourceDescription = "McLevy_ResoTM";

    TString OutputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/OtherTasks/Fit_pL_MickeyMouse/";
    TString OutFileName = "OutSPD.root";
    TFile* OutputFile = new TFile(OutputFolder+OutFileName,"recreate");
    hData_pL->Write();

    CATS AB_pL;
    DLM_Ck* Ck_pL;
    AB_pL.SetMomBins(NumMomBins_pL,MomBins_pL);
    //AnalysisObject.SetUpCats_pL(AB_pL,"NLO_Coupled_S",SourceDescription,0,202);//NLO_Coupled_S
    AnalysisObject.SetUpCats_pL(AB_pL,"NLO_Coupled_SPD",SourceDescription,0,202);//NLO_Coupled_S
//AB_pL.SetChannelWeight(8,AB_pL.GetChannelWeight(8)*0.54);
//AB_pL.SetChannelWeight(13,AB_pL.GetChannelWeight(13)*0.54);
const double CUSP_WEIGHT = 0.33;//0.54
AB_pL.SetChannelWeight(7,1./4.*CUSP_WEIGHT);//1S0 SN(s) -> LN(s)
AB_pL.SetChannelWeight(8,3./4.*CUSP_WEIGHT);//3S1 SN(s) -> LN(s)
AB_pL.SetChannelWeight(10,3./4.*CUSP_WEIGHT);//3S1 SN(d) -> LN(s)
AB_pL.SetChannelWeight(13,3./20.*CUSP_WEIGHT);//3D1 SN(d) -> LN(d)
AB_pL.SetChannelWeight(15,3./20.*CUSP_WEIGHT);//3D1 SN(s) -> LN(d)

    AB_pL.SetAnaSource(0,1.4);
    if(SourceDescription.Contains("Mc")){
        AB_pL.SetAnaSource(0,1.10);//c.a. 10% smaller compared to p-p due to the mT scaling
        AB_pL.SetAnaSource(1,2.0);
    }
    AB_pL.SetNotifications(CATS::nWarning);
    AB_pL.KillTheCat();
    Ck_pL = new DLM_Ck(AB_pL.GetNumSourcePars(),0,AB_pL,NumBinsCk,0,MaxBinValCk);
    Ck_pL->SetSourcePar(0,AB_pL.GetAnaSourcePar(0));
    if(SourceDescription.Contains("Mc")){
        Ck_pL->SetSourcePar(0,AB_pL.GetAnaSourcePar(0));
        Ck_pL->SetSourcePar(1,AB_pL.GetAnaSourcePar(1));
    }
    Ck_pL->SetCutOff(320,500);

    CATS AB_pXim;
    //same binning as pL, as we only use pXim as feed-down
    AB_pXim.SetMomBins(NumMomBins_pL,MomBins_pL);
    AnalysisObject.SetUpCats_pXim(AB_pXim,"pXim_HALQCD1","Gauss");
    AB_pXim.SetAnaSource(0,ResidualSourceSize);
    AB_pXim.SetNotifications(CATS::nWarning);
    AB_pXim.KillTheCat();
    DLM_Ck* Ck_pXim = new DLM_Ck(AB_pXim.GetNumSourcePars(),0,AB_pXim);
    Ck_pL->Update();
    Ck_pSigma0->Update();
    Ck_pXim->Update();

    DLM_CkDecomposition CkDec_pL("pLambda",4,*Ck_pL,hResolution_pL);
    DLM_CkDecomposition CkDec_pSigma0("pSigma0",0,*Ck_pSigma0,NULL);
    DLM_CkDecomposition CkDec_pXim("pXim",2,*Ck_pXim,NULL);

    CkDec_pL.AddContribution(0,lam_pL[1],DLM_CkDecomposition::cFeedDown,&CkDec_pSigma0,hResidual_pL_pSigma0);
    CkDec_pL.AddContribution(1,lam_pL[2],DLM_CkDecomposition::cFeedDown,&CkDec_pXim,hResidual_pL_pXim);
    //CkDec_pL.AddContribution(0,lam_pL[1],DLM_CkDecomposition::cFeedDown);
    //CkDec_pL.AddContribution(1,lam_pL[2],DLM_CkDecomposition::cFeedDown);
    CkDec_pL.AddContribution(2,lam_pL[3],DLM_CkDecomposition::cFeedDown);
    CkDec_pL.AddContribution(3,lam_pL[4],DLM_CkDecomposition::cFake);//0.03

    //for Xim we simplify a bit and take ALL feed-down as flat
    CkDec_pXim.AddContribution(0,lam_pXim[1]+lam_pXim[2]+lam_pXim[3],DLM_CkDecomposition::cFeedDown);
    CkDec_pXim.AddContribution(1,lam_pXim[4],DLM_CkDecomposition::cFake);

    CkDec_pL.Update();

    MM_PL = &CkDec_pL;
    MM_CatPL = &AB_pL;
    OutputFile->cd();
    const double FitMin = 0;
    const double FitMax = 2000;
    TF1* fit_pL = new TF1("fit_pL",MickeyFitter_pL,FitMin,FitMax,15);
    fit_pL->SetParameter(0,0.1);
    fit_pL->SetParLimits(0,0.0,1.0);

    fit_pL->SetParameter(1,0.3);
    fit_pL->SetParLimits(1,0.1,0.5);

    fit_pL->FixParameter(2,2.0);

    fit_pL->SetParameter(3,1.0);
    fit_pL->SetParLimits(3,0.01,10.0);

    fit_pL->SetParameter(4,200.);
    fit_pL->SetParLimits(4,-100.,400.);

    fit_pL->SetParameter(5,200.);
    fit_pL->SetParLimits(5,80.,600.);

    fit_pL->SetParameter(6,1.0);
    fit_pL->SetParLimits(6,0.0,2.0);

    fit_pL->SetParameter(7,0.0);
    fit_pL->SetParLimits(7,-1e-2,1e-2);

    fit_pL->SetParameter(8,0.0);
    fit_pL->SetParLimits(8,-1e-4,1e-4);

    fit_pL->FixParameter(9,0.0);

    fit_pL->SetParameter(10,1000);
    fit_pL->SetParLimits(10,400,1500);

    fit_pL->SetParameter(11,1);
    fit_pL->SetParLimits(11,0.5,1.5);

    fit_pL->SetParameter(12,Ck_pL->GetSourcePar(0));
    fit_pL->SetParLimits(12,0.8*Ck_pL->GetSourcePar(0),1.2*Ck_pL->GetSourcePar(0));

    fit_pL->SetParameter(13,Ck_pL->GetSourcePar(1));
    fit_pL->SetParLimits(13,1.0,2.0);

    fit_pL->SetParameter(14,1);
    fit_pL->SetParLimits(14,0.8,1.2);


/*
//used for the resonances
//fit_pL->FixParameter(0,0.2);
//fit_pL->FixParameter(1,0.32);
//fit_pL->FixParameter(2,2.0);
//fit_pL->FixParameter(3,0.195);
//fit_pL->FixParameter(4,162*0.8);
//fit_pL->FixParameter(5,185*1.15);
//fit_pL->FixParameter(6,1.0);
//fit_pL->FixParameter(7,0.0);
fit_pL->FixParameter(8,0.0);
fit_pL->FixParameter(9,0.0);
fit_pL->FixParameter(10,500);
fit_pL->FixParameter(11,0.67);
fit_pL->FixParameter(12,Ck_pL->GetSourcePar(0));//1.36
fit_pL->FixParameter(13,2.00);//1.36
fit_pL->FixParameter(14,1);
*/
fit_pL->FixParameter(4,0);//jet mean

fit_pL->FixParameter(8,0.0);//k^2
fit_pL->FixParameter(9,0.0);//k^3

fit_pL->FixParameter(11,0.67);
fit_pL->FixParameter(12,1.11);
fit_pL->FixParameter(13,2.00);//1.36
fit_pL->FixParameter(14,1);
    hData_pL->Fit(fit_pL,"S, N, R, M");
    printf("chi2/ndf = %.2f / %i\n",fit_pL->GetChisquare(),fit_pL->GetNDF());
    printf("prob = %.4f\n",fit_pL->GetProb());
    printf("nsigma = %.2f\n",sqrt(2)*TMath::ErfcInverse(fit_pL->GetProb()));
    fit_pL->SetNpx(1024);

    TF1* fit_pL_Jet = new TF1("fit_pL_Jet","1.+[0]*TMath::Gaus(x,[1],[2],0)",0,4500);
    fit_pL_Jet->SetParameter(0,fit_pL->GetParameter(3));
    fit_pL_Jet->SetParameter(1,fit_pL->GetParameter(4));
    fit_pL_Jet->SetParameter(2,fit_pL->GetParameter(5));

    TF1* fit_pL_QS = new TF1("fit_pL_QS","(1.-[0]*exp(-pow(x*[1]/197.327,[2])))",0,4500);
    fit_pL_QS->SetParameter(0,fit_pL->GetParameter(0));
    fit_pL_QS->SetParameter(1,fit_pL->GetParameter(1));
    fit_pL_QS->SetParameter(2,fit_pL->GetParameter(2));

    TF1* fit_pL_NonFemtoNonJet = new TF1("fit_pL_NonFemtoNonJet","(1.-[0]*exp(-pow(x*[1]/197.327,[2])))*(1.+[3]*TMath::Gaus(x,[4],[5],0))*[6]*(1.+[7]*x+[8]*x*x+[9]*x*x*x)",0,4500);
    fit_pL_NonFemtoNonJet->SetParameter(0,fit_pL->GetParameter(0));
    fit_pL_NonFemtoNonJet->SetParameter(1,fit_pL->GetParameter(1));
    fit_pL_NonFemtoNonJet->SetParameter(2,fit_pL->GetParameter(2));
    fit_pL_NonFemtoNonJet->SetParameter(3,0);
    fit_pL_NonFemtoNonJet->SetParameter(4,fit_pL->GetParameter(4));
    fit_pL_NonFemtoNonJet->SetParameter(5,fit_pL->GetParameter(5));
    fit_pL_NonFemtoNonJet->SetParameter(6,fit_pL->GetParameter(6));
    fit_pL_NonFemtoNonJet->SetParameter(7,fit_pL->GetParameter(7));
    fit_pL_NonFemtoNonJet->SetParameter(8,fit_pL->GetParameter(8));
    fit_pL_NonFemtoNonJet->SetParameter(9,fit_pL->GetParameter(9));

    TF1* fit_pL_NonFemto = new TF1("fit_pL_NonFemto","(1.-[0]*exp(-pow(x*[1]/197.327,[2])))*(1.+[3]*TMath::Gaus(x,[4],[5],0))*[6]*(1.+[7]*x+[8]*x*x+[9]*x*x*x)",0,4500);
    fit_pL_NonFemto->SetParameter(0,fit_pL->GetParameter(0));
    fit_pL_NonFemto->SetParameter(1,fit_pL->GetParameter(1));
    fit_pL_NonFemto->SetParameter(2,fit_pL->GetParameter(2));
    fit_pL_NonFemto->SetParameter(3,fit_pL->GetParameter(3));
    fit_pL_NonFemto->SetParameter(4,fit_pL->GetParameter(4));
    fit_pL_NonFemto->SetParameter(5,fit_pL->GetParameter(5));
    fit_pL_NonFemto->SetParameter(6,fit_pL->GetParameter(6));
    fit_pL_NonFemto->SetParameter(7,fit_pL->GetParameter(7));
    fit_pL_NonFemto->SetParameter(8,fit_pL->GetParameter(8));
    fit_pL_NonFemto->SetParameter(9,fit_pL->GetParameter(9));

    TGraph fit_nsigma;
    fit_nsigma.SetName("fit_nsigma");

    TGraphErrors fit_sourcepars;
    fit_sourcepars.SetName("fit_sourcepars");
    fit_sourcepars.SetPoint(0,fit_pL->GetParameter(12),fit_pL->GetParameter(13));
    fit_sourcepars.SetPointError(0,fit_pL->GetParError(12),fit_pL->GetParError(13));


    TH1F* hfit_Ratio = new TH1F("hfit_Ratio","hfit_Ratio",
            hData_pL->GetNbinsX(),hData_pL->GetBinLowEdge(1),hData_pL->GetXaxis()->GetBinUpEdge(hData_pL->GetNbinsX()));

    int NumPts=0;
    double Avg_nsigma = 0;//should be around 0
    for(unsigned uBin=0; uBin<hData_pL->GetNbinsX(); uBin++){
        double MOM = hData_pL->GetBinCenter(uBin+1);
        if(MOM<FitMin || MOM>FitMax) continue;
        double CkData = hData_pL->GetBinContent(uBin+1);
        double CkErr = hData_pL->GetBinError(uBin+1);
        double CkFit = fit_pL->Eval(MOM);
        double nsigma = (CkData-CkFit)/CkErr;
        fit_nsigma.SetPoint(NumPts,MOM,nsigma);
        Avg_nsigma += nsigma;

        hfit_Ratio->SetBinContent(uBin+1,CkData/CkFit);
        hfit_Ratio->SetBinError(uBin+1,CkErr/CkFit);

        NumPts++;
    }
    Avg_nsigma /= double(NumPts);
    printf("Avg_nsigma = %.3f\n",Avg_nsigma);

    OutputFile->cd();
    fit_pL->Write();
    fit_pL_Jet->Write();
    fit_pL_QS->Write();
    fit_pL_NonFemtoNonJet->Write();
    fit_pL_NonFemto->Write();
    fit_nsigma.Write();
    fit_sourcepars.Write();
    hfit_Ratio->Write();
    printf("Closing\n");
    OutputFile->Close();
    printf("Closed\n");

    //delete fit_pL;
}



DLM_CkDecomposition* MM_PP;
CATS* MM_CatPP;
double MickeyFitter_pp(double* x, double* par){
    double& MOM = *x;
    MM_PP->GetCk()->SetCutOff(280,par[10]);
    MM_PP->GetCk()->SetSourcePar(0,par[11]);
    //MM_CatPP->SetAnaSource(0,par[11],true);
//printf(" MM_CatPP->SourceStatus() = %i (%.3f)\n\n",MM_CatPP->SourceStatus(),par[11]);
    if(MM_CatPP->GetNumSourcePars()>1){
        //MM_CatPP->SetAnaSource(1,par[12],true);
        MM_PP->GetCk()->SetSourcePar(1,par[12]);
    }
    //
    //MM_CatPP->KillTheCat();
    MM_PP->Update(false);
    static int COUNTER = 0;
    //if(MOM<MM_PP->GetCk()->GetBinUpEdge(0,0)){
    //if(MOM>MM_PP->GetCk()->GetBinLowEdge(0,10)&&MOM<MM_PP->GetCk()->GetBinUpEdge(0,10)){
    if(MOM>30&&MOM<35){
    COUNTER++;
    printf("COUNTER=%i\n",COUNTER);
    }

    double CkVal;
    CkVal = MM_PP->EvalCk(MOM);
    double Pauli = 1.-par[0]*exp(-pow(MOM*par[1]/197.327,par[2]));
    double Jet = 1.+par[3]*TMath::Gaus(MOM,par[4],par[5],0);
    double NonJet = par[6]*(1.+par[7]*MOM+par[8]*MOM*MOM+par[9]*MOM*MOM*MOM);
    return CkVal*NonJet*Jet*Pauli;
}
//including QS
void Fit_pp_MickeyMouse(){
    double ResidualSourceSize=1.35;

    double* MomBins_pp = NULL;
    double* FitRegion_pp = NULL;
    unsigned NumMomBins_pp;
    TString DataSample = "pp13TeV_HM_Dec19";
    //TString DataSample = "pp13TeV_HM_RotPhiDec19";
    DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder("/home/dmihaylov/CernBox/CatsFiles");
    AnalysisObject.SetUpBinning_pp(DataSample,NumMomBins_pp,MomBins_pp,FitRegion_pp);

    TH2F* hResolution_pp = AnalysisObject.GetResolutionMatrix(DataSample,"pp");
    TH2F* hResidual_pp_pL = AnalysisObject.GetResidualMatrix("pp","pLambda");

    TH1F* hData_pp = AnalysisObject.GetAliceExpCorrFun(DataSample,"pp","_0",0,false,-1);

    double lam_pp[5];
    double lam_pL[5];

    AnalysisObject.SetUpLambdaPars_pp(DataSample,0,lam_pp);
    AnalysisObject.SetUpLambdaPars_pL(DataSample,0,0,lam_pL);

    //TString SourceDescription = "Gauss";
    TString OutputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/OtherTasks/Fit_pp_MickeyMouse/";
    TString OutFileName = "Out.root";
    TFile* OutputFile = new TFile(OutputFolder+OutFileName,"recreate");
    hData_pp->Write();

    //Gauss
    TString SourceDescription = "Gauss";
    //TString SourceDescription = "McGauss_ResoTM";
    //TString SourceDescription = "McLevy_ResoTM";

    CATS AB_pp;
    DLM_Ck* Ck_pp;
    AB_pp.SetMomBins(NumMomBins_pp,MomBins_pp);
    AB_pp.SetNotifications(CATS::nWarning);
    //AnalysisObject.SetUpCats_pp(AB_pp,"AV18","Gauss",0,0);//McLevyNolan_Reso
    AnalysisObject.SetUpCats_pp(AB_pp,"AV18",SourceDescription,0,202);//McLevyNolan_Reso

    AB_pp.SetAnaSource(0,1.3);
    if(SourceDescription.Contains("Mc")){
        AB_pp.SetAnaSource(0,1.20);
        AB_pp.SetAnaSource(1,2.0);
    }
    //AB_pp.SetAnaSource(1,2.0);
    //AB_pp.SetNotifications(CATS::nWarning);
    //AB_pp.SetEpsilonConv(5e-8);
    //AB_pp.SetEpsilonProp(5e-8);
    AB_pp.KillTheCat();
    Ck_pp = new DLM_Ck(AB_pp.GetNumSourcePars(),0,AB_pp);
    Ck_pp->Update();

    CATS AB_pL;
    DLM_Ck* Ck_pL;
    AB_pL.SetMomBins(NumMomBins_pp,MomBins_pp);
    //AnalysisObject.SetUpCats_pL(AB_pL,"NLO_Coupled_S","Gauss");
    AnalysisObject.SetUpCats_pL(AB_pL,"Usmani","Gauss");
    AB_pL.SetAnaSource(0,ResidualSourceSize);
    AB_pL.SetNotifications(CATS::nWarning);
    AB_pL.KillTheCat();
    Ck_pL = new DLM_Ck(AB_pL.GetNumSourcePars(),0,AB_pL);
    Ck_pL->Update();

    DLM_CkDecomposition CkDec_pp("pp",3,*Ck_pp,hResolution_pp);
    DLM_CkDecomposition CkDec_pL("pLambda",2,*Ck_pL,NULL);

    CkDec_pp.AddContribution(0,lam_pp[1],DLM_CkDecomposition::cFeedDown,&CkDec_pL,hResidual_pp_pL);
    CkDec_pp.AddContribution(1,lam_pp[2],DLM_CkDecomposition::cFeedDown);
    CkDec_pp.AddContribution(2,lam_pp[3],DLM_CkDecomposition::cFake);

    CkDec_pL.AddContribution(2,lam_pL[1]+lam_pL[2]+lam_pL[3],DLM_CkDecomposition::cFeedDown);
    CkDec_pL.AddContribution(3,lam_pL[4],DLM_CkDecomposition::cFake);//0.03

    CkDec_pp.Update();
    CkDec_pL.Update();

    MM_PP = &CkDec_pp;
    MM_CatPP = &AB_pp;
    OutputFile->cd();

/*
   1  p0           1.76256e-01   5.85904e-01  -0.00000e+00   8.51297e+03
   2  p1           3.00695e-01   2.84307e-01   0.00000e+00  -3.14685e+03
   3  p2           2.00000e+00     fixed
   4  p3           1.04360e-01   5.71488e+00   0.00000e+00  -2.16095e+03
   5  p4           5.00001e+01   1.53683e+02   0.00000e+00** at limit **
   6  p5           1.65192e+02   1.82562e+02   0.00000e+00  -7.73015e+02
   7  p6           1.13557e+00   1.44802e+00   0.00000e+00  -3.42455e+04
   8  p7           1.05201e-05   1.42098e-02   0.00000e+00  -2.10732e+05
   9  p8           0.00000e+00     fixed
  10  p9           0.00000e+00     fixed
  11  p10          9.99999e+02   3.69205e+02  -0.00000e+00   1.28217e-01
  12  p11          1.27085e+00   2.76548e-01  -0.00000e+00   8.47135e+03
chi2/ndf = 649.11 / 612
nsigma = 1.46

*/

    double FitMin = 0;
    double FitMax = 2500;

    TF1* fit_pp = new TF1("fit_pp",MickeyFitter_pp,FitMin,FitMax,13);
    fit_pp->SetParameter(0,0.18);
    fit_pp->SetParLimits(0,0.0,1.0);

    fit_pp->SetParameter(1,0.3);
    fit_pp->SetParLimits(1,0.1,0.5);

    fit_pp->FixParameter(2,2.0);

    fit_pp->SetParameter(3,0.1);
    fit_pp->SetParLimits(3,0.01,10.0);

    fit_pp->SetParameter(4,50.);
    fit_pp->SetParLimits(4,0.,200.);

    fit_pp->SetParameter(5,160.);
    fit_pp->SetParLimits(5,100.,300.);

    fit_pp->SetParameter(6,1.13);
    fit_pp->SetParLimits(6,0.8,1.3);

    fit_pp->SetParameter(7,1e-5);
    fit_pp->SetParLimits(7,-1e-3,1e-3);

    fit_pp->SetParameter(8,0.0);
    fit_pp->SetParLimits(8,-1e-4,1e-4);

    fit_pp->SetParameter(9,0.0);
    fit_pp->SetParLimits(9,-1e-7,1e-7);

    fit_pp->SetParameter(10,1000);
    fit_pp->SetParLimits(10,400,1600);

    fit_pp->SetParameter(11,Ck_pp->GetSourcePar(0));
    fit_pp->SetParLimits(11,Ck_pp->GetSourcePar(0)*0.9,Ck_pp->GetSourcePar(0)*1.1);

    fit_pp->SetParameter(12,Ck_pp->GetSourcePar(1));
    fit_pp->SetParLimits(12,1.0,2.0);
    if(SourceDescription.Contains("Gauss")) fit_pp->FixParameter(12,2.0);

//fit_pp->FixParameter(0,0.0);
//fit_pp->FixParameter(1,0.3);
//fit_pp->FixParameter(2,2.0);
//fit_pp->FixParameter(3,0.0);
//fit_pp->FixParameter(4,100);
//fit_pp->FixParameter(5,150);
//fit_pp->FixParameter(6,1.0);
//fit_pp->FixParameter(7,0.0);
//fit_pp->FixParameter(8,0.0);
fit_pp->FixParameter(9,0.0);
//fit_pp->FixParameter(10,1000);
fit_pp->FixParameter(11,1.29);
fit_pp->FixParameter(12,2.0);

    hData_pp->Fit(fit_pp,"S, N, R, M");
    printf("chi2/ndf = %.2f / %i\n",fit_pp->GetChisquare(),fit_pp->GetNDF());
    printf("nsigma = %.2f\n",sqrt(2)*TMath::ErfcInverse(fit_pp->GetProb()));
    fit_pp->SetNpx(3000);

    TF1* fit_pp_Jet = new TF1("fit_pp_Jet","1.+[0]*TMath::Gaus(x,[1],[2],0)",0,4500);
    fit_pp_Jet->SetParameter(0,fit_pp->GetParameter(3));
    fit_pp_Jet->SetParameter(1,fit_pp->GetParameter(4));
    fit_pp_Jet->SetParameter(2,fit_pp->GetParameter(5));

    TF1* fit_pp_QS = new TF1("fit_pp_QS","(1.-[0]*exp(-pow(x*[1]/197.327,[2])))",0,4500);
    fit_pp_QS->SetParameter(0,fit_pp->GetParameter(0));
    fit_pp_QS->SetParameter(1,fit_pp->GetParameter(1));
    fit_pp_QS->SetParameter(2,fit_pp->GetParameter(2));

    TF1* fit_pp_Femto = new TF1("fit_pp_Femto",MickeyFitter_pp,0,3000,12);
    fit_pp_Femto->FixParameter(0,0.0);//
    fit_pp_Femto->FixParameter(1,0.3);
    fit_pp_Femto->FixParameter(2,2.0);
    fit_pp_Femto->FixParameter(3,0.0);//
    fit_pp_Femto->FixParameter(4,100);
    fit_pp_Femto->FixParameter(5,150);
    fit_pp_Femto->FixParameter(6,1.0);//
    fit_pp_Femto->FixParameter(7,0.0);
    fit_pp_Femto->FixParameter(8,0.0);
    fit_pp_Femto->FixParameter(9,0.0);
    fit_pp_Femto->FixParameter(10,fit_pp->GetParameter(10));
    fit_pp_Femto->FixParameter(11,fit_pp->GetParameter(11));
    fit_pp_Femto->FixParameter(12,fit_pp->GetParameter(12));
    fit_pp_Femto->SetNpx(3000);

    TF1* fit_pp_NonFemto = new TF1("fit_pp_NonFemto","(1.-[0]*exp(-pow(x*[1]/197.327,[2])))*(1.+[3]*TMath::Gaus(x,[4],[5],0))*[6]*(1.+[7]*x+[8]*x*x+[9]*x*x*x)",0,4500);
    fit_pp_NonFemto->SetParameter(0,fit_pp->GetParameter(0));
    fit_pp_NonFemto->SetParameter(1,fit_pp->GetParameter(1));
    fit_pp_NonFemto->SetParameter(2,fit_pp->GetParameter(2));
    fit_pp_NonFemto->SetParameter(3,fit_pp->GetParameter(3));
    fit_pp_NonFemto->SetParameter(4,fit_pp->GetParameter(4));
    fit_pp_NonFemto->SetParameter(5,fit_pp->GetParameter(5));
    fit_pp_NonFemto->SetParameter(6,1);
    fit_pp_NonFemto->SetParameter(7,fit_pp->GetParameter(7));
    fit_pp_NonFemto->SetParameter(8,fit_pp->GetParameter(8));
    fit_pp_NonFemto->SetParameter(9,fit_pp->GetParameter(9));

    TGraph gCk_Theory;
    TGraph gCk_Theory_pL;
    TGraph gCk_Smeared;
    TGraph gCk_Main;
    TGraph gCk_SmearedMain;
    TGraph gCk_MainFeed;
    TGraph gCk_ChildMainFeed;
    gCk_Theory.SetName("gCk_Theory");
    gCk_Theory_pL.SetName("gCk_Theory_pL");
    gCk_Smeared.SetName("gCk_Smeared");
    gCk_Main.SetName("gCk_Main");
    gCk_SmearedMain.SetName("gCk_SmearedMain");
    gCk_MainFeed.SetName("gCk_MainFeed");
    gCk_ChildMainFeed.SetName("gCk_ChildMainFeed");
    for(unsigned uBin=0; uBin<Ck_pp->GetNbins(); uBin++){
        gCk_Theory.SetPoint(uBin,Ck_pp->GetBinCenter(0,uBin),Ck_pp->GetBinContent(uBin));
        gCk_Theory_pL.SetPoint(uBin,Ck_pL->GetBinCenter(0,uBin),Ck_pL->GetBinContent(uBin));
        gCk_Smeared.SetPoint(uBin,Ck_pp->GetBinCenter(0,uBin),CkDec_pp.EvalCk(Ck_pp->GetBinCenter(0,uBin)));
        gCk_Main.SetPoint(uBin,Ck_pp->GetBinCenter(0,uBin),CkDec_pp.EvalMain(Ck_pp->GetBinCenter(0,uBin)));
        gCk_SmearedMain.SetPoint(uBin,Ck_pp->GetBinCenter(0,uBin),CkDec_pp.EvalSmearedMain(Ck_pp->GetBinCenter(0,uBin)));
        gCk_MainFeed.SetPoint(uBin,Ck_pp->GetBinCenter(0,uBin),CkDec_pp.EvalMainFeed(Ck_pp->GetBinCenter(0,uBin)));
        gCk_ChildMainFeed.SetPoint(uBin,Ck_pp->GetBinCenter(0,uBin),CkDec_pp.GetChild(0)->EvalMainFeed(Ck_pp->GetBinCenter(0,uBin)));
        //if(uBin<AB_pp.GetNumMomBins()){
        //    printf("k=%.0f: %.3f vs %.3f\n", AB_pp.GetMomentum(uBin), AB_pp.GetCorrFun(uBin), Ck_pp->GetBinContent(uBin));
        //}
    }

    TGraph fit_nsigma;
    fit_nsigma.SetName("fit_nsigma");
    TH1F* hfit_Ratio = new TH1F("hfit_Ratio","hfit_Ratio",
            hData_pp->GetNbinsX(),hData_pp->GetBinLowEdge(1),hData_pp->GetXaxis()->GetBinUpEdge(hData_pp->GetNbinsX()));

    int NumPts=0;
    double Avg_nsigma = 0;//should be around 0
    for(unsigned uBin=0; uBin<hData_pp->GetNbinsX(); uBin++){
        double MOM = hData_pp->GetBinCenter(uBin+1);
        if(MOM<FitMin || MOM>FitMax) continue;
        double CkData = hData_pp->GetBinContent(uBin+1);
        double CkErr = hData_pp->GetBinError(uBin+1);
        double CkFit = fit_pp->Eval(MOM);
        double nsigma = (CkData-CkFit)/CkErr;
        fit_nsigma.SetPoint(NumPts,MOM,nsigma);
        Avg_nsigma += nsigma;

        hfit_Ratio->SetBinContent(uBin+1,CkData/CkFit);
        hfit_Ratio->SetBinError(uBin+1,CkErr/CkFit);

        NumPts++;
    }
    Avg_nsigma /= double(NumPts);
    printf("Avg_nsigma = %.3f\n",Avg_nsigma);

    TGraph gSourcePars;
    gSourcePars.SetName("gSourcePars");
    gSourcePars.Set(1);
    gSourcePars.SetMarkerStyle(20);
    gSourcePars.SetMarkerSize(1.5);
    gSourcePars.SetPoint(0,fit_pp->GetParameter(11),fit_pp->GetParameter(12));


    OutputFile->cd();
    fit_pp->Write();
    fit_pp_Jet->Write();
    fit_pp_QS->Write();
    fit_pp_Femto->Write();
    fit_pp_NonFemto->Write();
    fit_nsigma.Write();
    hfit_Ratio->Write();
    gSourcePars.Write();
    gCk_Theory.Write();
    gCk_Theory_pL.Write();
    gCk_Smeared.Write();
    gCk_Main.Write();
    gCk_SmearedMain.Write();
    gCk_MainFeed.Write();
    gCk_ChildMainFeed.Write();

    //delete fit_pp;
    //delete fit_pp_Jet;
    //delete fit_pp_Femto;
}

void Compare_RotPhi_ME(){
    DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder("/home/dmihaylov/CernBox/CatsFiles");

    //pp
    //TH1F* hData_pp_ME = AnalysisObject.GetAliceExpCorrFun("pp13TeV_HM_Dec19","pp","_0",0,false,-1);
    //TH1F* hData_pp_RE = AnalysisObject.GetAliceExpCorrFun("pp13TeV_HM_RotPhiDec19","pp","_0",0,false,-1);
    //pLambda
    //TH1F* hData_pp_ME = AnalysisObject.GetAliceExpCorrFun("pp13TeV_HM_Dec19","pLambda","_0",3,false,-1);
    //TH1F* hData_pp_RE = AnalysisObject.GetAliceExpCorrFun("pp13TeV_HM_RotPhiDec19","pLambda","_0",0,false,-1);

    //MC pp
    TFile* InFile = new TFile(
    "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/Fast_BBar/Trains_MCAOD/Norm018028/CFOutput_pp_8.root"
    ,"read");
    TH1F* hDummy = (TH1F*)InFile->Get("hCk_ReweightedMeV_0");
    gROOT->cd();
    TH1F *hData_pp_ME = (TH1F*)hDummy->Clone("hData_pp_ME");
    delete InFile; InFile=NULL;

    InFile = new TFile(
    "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/MC/CF/NanoMC/CFOutput_pp_3.root"
    ,"read");
    hDummy = (TH1F*)InFile->Get("hCk_ReweightedMeV_0");
    gROOT->cd();
    TH1F *hData_pp_RE = (TH1F*)hDummy->Clone("hData_pp_RE");
    delete InFile; InFile=NULL;

/*
    //MC pLambda
    TFile* InFile = new TFile(
    "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/Fast_BBar/Trains_MCAOD/Norm018028/CFOutput_pL_8.root"
    ,"read");
    TH1F* hDummy = (TH1F*)InFile->Get("hCk_ReweightedMeV_0");
    gROOT->cd();
    TH1F *hData_pp_ME = (TH1F*)hDummy->Clone("hData_pp_ME");
    delete InFile; InFile=NULL;

    InFile = new TFile(
    "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/MC/CF/NanoMC/CFOutput_pL_3.root"
    ,"read");
    hDummy = (TH1F*)InFile->Get("hCk_ReweightedMeV_0");
    gROOT->cd();
    TH1F *hData_pp_RE = (TH1F*)hDummy->Clone("hData_pp_RE");
    delete InFile; InFile=NULL;
*/
    DLM_Histo<float> dlmData_pp_ME;
    DLM_Histo<float> dlmData_pp_RE;
    DLM_Histo<float> dlmData_pp_MEvsRE;

    unsigned NumBins = hData_pp_ME->GetNbinsX()>hData_pp_RE->GetNbinsX()?hData_pp_RE->GetNbinsX():hData_pp_ME->GetNbinsX();
    double kMin = hData_pp_RE->GetBinLowEdge(1);
    double kMax = hData_pp_RE->GetXaxis()->GetBinUpEdge(NumBins);

unsigned RebinFactor = 4;
NumBins /= RebinFactor;
kMax =  hData_pp_RE->GetXaxis()->GetBinUpEdge(NumBins*RebinFactor);

    printf("NumBins=%u; kMin kMax = %f %f\n",NumBins,kMin,kMax);

    dlmData_pp_ME.SetUp(1);
    dlmData_pp_ME.SetUp(0,NumBins,kMin,kMax);
    dlmData_pp_ME.Initialize();

    dlmData_pp_RE.SetUp(1);
    dlmData_pp_RE.SetUp(0,NumBins,kMin,kMax);
    dlmData_pp_RE.Initialize();
    for(unsigned uBin=0; uBin<NumBins; uBin++){

        double BinVal=0; double BinErr=0;
        for(unsigned uReb=0; uReb<RebinFactor; uReb++){
            BinVal += hData_pp_ME->GetBinContent(uBin*RebinFactor+uReb+1);
            BinErr += hData_pp_ME->GetBinError(uBin*RebinFactor+uReb+1)*hData_pp_ME->GetBinError(uBin*RebinFactor+uReb+1);
        }
        BinErr = sqrt(BinErr);
        dlmData_pp_ME.SetBinContent(uBin,BinVal);
        dlmData_pp_ME.SetBinError(uBin,BinErr);

        BinVal=0; BinErr=0;
        for(unsigned uReb=0; uReb<RebinFactor; uReb++){
            BinVal += hData_pp_RE->GetBinContent(uBin*RebinFactor+uReb+1);
            BinErr += hData_pp_RE->GetBinError(uBin*RebinFactor+uReb+1)*hData_pp_RE->GetBinError(uBin*RebinFactor+uReb+1);
        }
        BinVal += 0.001;
        BinErr = sqrt(BinErr);
        dlmData_pp_RE.SetBinContent(uBin,BinVal);
        dlmData_pp_RE.SetBinError(uBin,BinErr);


        //dlmData_pp_ME.SetBinContent(uBin,hData_pp_ME->GetBinContent(uBin+1));
        //dlmData_pp_ME.SetBinError(uBin,hData_pp_ME->GetBinError(uBin+1));
        //dlmData_pp_RE.SetBinContent(uBin,hData_pp_RE->GetBinContent(uBin+1)+0.001);
        //dlmData_pp_RE.SetBinError(uBin,hData_pp_RE->GetBinError(uBin+1));
    }

    dlmData_pp_MEvsRE = dlmData_pp_ME;
    dlmData_pp_MEvsRE /= dlmData_pp_RE;

    TFile fOutput("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/OtherTasks/Compare_RotPhi_ME/fOutput.root","recreate");
    TH1F* hData_pp_MEvsRE = new TH1F("hData_pp_MEvsRE","hData_pp_MEvsRE",NumBins,kMin,kMax);
    for(unsigned uBin=0; uBin<NumBins; uBin++){
        hData_pp_MEvsRE->SetBinContent(uBin+1,dlmData_pp_MEvsRE.GetBinContent(uBin));
        hData_pp_MEvsRE->SetBinError(uBin+1,dlmData_pp_MEvsRE.GetBinError(uBin));
    }
    hData_pp_MEvsRE->Write();

    delete hData_pp_MEvsRE;
}



void ReadDeuteronWF(const char* InputFileName, DLM_Histo<float>& OutputU, DLM_Histo<float>& OutputW){

    FILE *InFile;
    InFile = fopen(InputFileName, "r");
    if(!InFile){
        printf("\033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", InputFileName);
        return;
    }

    char* cdummy = new char [1024];
    unsigned CurrentLine=0;
    float fRadius, fuWF, fwWF;
    double* Radius = new double [1024];
    double* RadiusBins = new double [1025];

    unsigned NumBins = 0;
    RadiusBins[0] = 0;
    //read line-by-line until the end of the file
    while(!feof(InFile)){
        CurrentLine++;
        //read the next line
        if(!fgets(cdummy, 1023, InFile)) continue;
        if(CurrentLine<=3) continue;
        //printf("Line#%u: %s\n",CurrentLine,cdummy);
        sscanf(cdummy, "%f %f %f",&fRadius, &fuWF, &fwWF);
        //printf("Line#%u: %.2e; %.2e; %.2e;\n",CurrentLine,fRadius,fuWF,fwWF);

        //if((CurrentLine)%int(NumRadBins)==0){
         //   sscanf(cdummy, "%f",&fMomentum);
        Radius[NumBins] = fRadius;
        if(NumBins){
            //set the bin range in between the last two bin centers
            RadiusBins[NumBins] = 0.5*(Radius[NumBins]+Radius[NumBins-1]);
        }
        //}

        NumBins++;
    }
    fclose(InFile);

    //set the upper edge of the last bin, where we just add the bin width of the last bin
    //i.e. if we have l(low) c(center) u(up), we have that u=c+(c-l)=2c-l
    RadiusBins[NumBins] = 2.*Radius[NumBins-1]-RadiusBins[NumBins-1];

    //for(unsigned uBin=0; uBin<NumBins; uBin++){
    //    printf("#%u: %f %f %f\n",uBin,RadiusBins[uBin],Radius[uBin],RadiusBins[uBin+1]);
    //}

    //dimension (1D)
    OutputU.SetUp(1);
    OutputU.SetUp(0,NumBins,RadiusBins,Radius);
    OutputU.Initialize();

    OutputW.SetUp(1);
    OutputW.SetUp(0,NumBins,RadiusBins,Radius);
    OutputW.Initialize();

    CurrentLine=0;
    NumBins=0;
    InFile = fopen(InputFileName, "r");
    while(!feof(InFile)){
        CurrentLine++;
        //read the next line
        if(!fgets(cdummy, 1023, InFile)) continue;
        if(CurrentLine<=3) continue;
        sscanf(cdummy, "%f %f %f",&fRadius, &fuWF, &fwWF);
        OutputU.SetBinContent(NumBins,fuWF);
        OutputW.SetBinContent(NumBins,fwWF);
        NumBins++;
    }
    fclose(InFile);

    delete [] cdummy;
    delete [] Radius;
}

double Fitter_MC_SimpleTemplate_pL(double* x, double* par){
    double& MOM = *x;
MM_CatPL->SetChannelWeight(7,1./4.*par[6]);//1S0 SN(s) -> LN(s)
MM_CatPL->SetChannelWeight(8,3./4.*par[6]);//3S1 SN(s) -> LN(s)
MM_CatPL->SetChannelWeight(10,3./4.*par[6]);//3S1 SN(d) -> LN(s)
MM_CatPL->SetChannelWeight(13,3./20.*par[6]);//3D1 SN(d) -> LN(d)
MM_CatPL->SetChannelWeight(15,3./20.*par[6]);//3D1 SN(s) -> LN(d)
    MM_CatPL->SetAnaSource(0,par[7],true);
    if(MM_CatPL->GetNumSourcePars()>1){
        MM_CatPL->SetAnaSource(1,par[8],true);
    }
    MM_CatPL->KillTheCat();
    MM_PL->GetCk()->SetCutOff(340,par[5]);
    //MM_PL->GetCk()->SetSourcePar(0,par[7]);
    MM_PL->Update(true);
    static int COUNTER = 0;
    if(MOM<MM_PL->GetCk()->GetBinUpEdge(0,0)){
    COUNTER++;
    //if(COUNTER%100==0)
    printf("COUNTER=%i\n",COUNTER);
    }

    double CkVal;
    //if(MOM<320) CkVal = MM_PL->EvalCk(MOM);
    //else{
    //    CkVal = (MOM-320-MM_PL->EvalCk(320)*(MOM-par[10]))/(par[10]-320);
    //}
    CkVal = MM_PL->EvalCk(MOM);
    if(MOM<280) CkVal *= par[9];
    //[0]*(1.+[1]*exp(-pow((x-[2])/[3],2.)))
    double Baseline = par[0]*(1.+par[1]*par[4]*exp(-pow((MOM-par[2])/(2.*par[3]),2.)));
    return CkVal*Baseline;
}

//including QS
void Fit_pL_MC_SimpleTemplate(){
    const bool Include_dWaves = true;
    double ResidualSourceSize = 1.5;

    double* MomBins_pL = NULL;
    double* FitRegion_pL = NULL;
    const unsigned NumBinsCk = 36;
    const double BinWidthCk = 12;
    const double MaxBinValCk = double(NumBinsCk)*BinWidthCk;
    unsigned NumMomBins_pL;
    //TString DataSample = "pp13TeV_HM_Dec19";
    TString DataSample = "pp13TeV_HM_RotPhiDec19";

    const double MCfit_Min = 0;
    const double MCfit_Max = 540;
    TString McFile = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/MC/CF/NanoMC/CFOutput_pL_3.root";
    TString McHisto = "hCk_ReweightedMeV_1";

    //so we need to copy our histogram, as else we lose it when we delete the file
    //and we need to change to the "central" root directory, as else histoCopy will also be lost
    //and we need to play with the name a little bit, else we are fucked!
    TFile* mcfile = new TFile(McFile,"read");
    TH1F* hCk_tmp = (TH1F*)mcfile->Get(McHisto);
    if(!hCk_tmp){printf("\033[1;31mERROR:\033[0m The hCk_tmp '%s' if file '%s' does not exist\n",McHisto.Data(),McFile.Data());return;}
    TString Name = hCk_tmp->GetName();
    gROOT->cd();
    TH1F *hCk_MC = (TH1F*)hCk_tmp->Clone("hCk_MC");
    mcfile->Close();
    //delete mcfile;

    TF1* fit_MC = new TF1("fit_MC","[0]*(1.+[1]*exp(-pow((x-[2])/(2.*[3]),2.)))",MCfit_Min,MCfit_Max);
    fit_MC->SetParameter(0,1);
    fit_MC->SetParLimits(0,0.5,2.);
    fit_MC->SetParameter(1,0.03);
    fit_MC->SetParLimits(1,0.005,0.1);
    //fit_MC->SetParameter(2,0);
    //fit_MC->SetParLimits(2,-50,50);
    fit_MC->FixParameter(2,0);
    fit_MC->SetParameter(3,200);
    fit_MC->SetParLimits(3,50,500);

    hCk_MC->Fit(fit_MC,"S, N, R, M");

    DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder("/home/dmihaylov/CernBox/CatsFiles");
    AnalysisObject.SetUpBinning_pL(DataSample,NumMomBins_pL,MomBins_pL,FitRegion_pL,0,0);

    DLM_Ck* Ck_pSigma0 = new DLM_Ck(1,0,NumMomBins_pL,MomBins_pL,Lednicky_gauss_Sigma0);
    Ck_pSigma0->SetSourcePar(0,ResidualSourceSize);

    TH2F* hResolution_pL = AnalysisObject.GetResolutionMatrix(DataSample,"pLambda");
    TH2F* hResidual_pL_pSigma0 = AnalysisObject.GetResidualMatrix("pLambda","pSigma0");
    TH2F* hResidual_pL_pXim = AnalysisObject.GetResidualMatrix("pLambda","pXim");

    TH1F* hData_pL = AnalysisObject.GetAliceExpCorrFun(DataSample,"pLambda","_0",1,false,-1);

    double lam_pL[5];
    double lam_pXim[5];

    AnalysisObject.SetUpLambdaPars_pL(DataSample,0,0,lam_pL);
    AnalysisObject.SetUpLambdaPars_pXim(DataSample,0,0,lam_pXim);

    //TString SourceDescription = "Gauss";
    //TString SourceDescription = "McGauss_ResoTM";
    TString SourceDescription = "McLevy_ResoTM";

    TString OutputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/OtherTasks/Fit_pL_MC_SimpleTemplate/";
    TString OutFileName = "fOutput.root";
    TFile* OutputFile = new TFile(OutputFolder+OutFileName,"recreate");
    hCk_MC->Write();
    fit_MC->Write();
    hData_pL->Write();


    CATS AB_pL;
    DLM_Ck* Ck_pL;
    AB_pL.SetMomBins(NumMomBins_pL,MomBins_pL);
    //AnalysisObject.SetUpCats_pL(AB_pL,"NLO_Coupled_S",SourceDescription,0,202);//NLO_Coupled_S
    AnalysisObject.SetUpCats_pL(AB_pL,"NLO_Coupled_SPD",SourceDescription,0,202);//NLO_Coupled_S
//AB_pL.SetChannelWeight(8,AB_pL.GetChannelWeight(8)*0.54);
//AB_pL.SetChannelWeight(13,AB_pL.GetChannelWeight(13)*0.54);
    AB_pL.SetAnaSource(0,1.4);
    if(SourceDescription.Contains("Mc")){
        AB_pL.SetAnaSource(0,1.15);//c.a. 10% smaller compared to p-p due to the mT scaling
        AB_pL.SetAnaSource(1,2.0);
    }

    for(unsigned short usCh=1; usCh<=6; usCh++){
        //remove the p-waves
        //AB_pL.RemoveExternalWaveFunction(usCh,1);
    }
    for(unsigned short usCh=1; usCh<=3; usCh++){
        //remove the d-waves
        if(!Include_dWaves) AB_pL.RemoveExternalWaveFunction(usCh,2);
    }
    for(unsigned short usCh=8; usCh<=15; usCh++){
        //remove all coupling
        //AB_pL.RemoveExternalWaveFunction(usCh,0);

        //remove all coupling, but the s-wave SN->LN
        //if(usCh!=8&&!Include_dWaves) AB_pL.RemoveExternalWaveFunction(usCh,0);
        //remove all coupling, but the s,d-wave SN->LN
        //else if(usCh!=8&&usCh!=13) AB_pL.RemoveExternalWaveFunction(usCh,0);
    }

    AB_pL.SetNotifications(CATS::nWarning);
    AB_pL.KillTheCat();
    Ck_pL = new DLM_Ck(AB_pL.GetNumSourcePars(),0,AB_pL,NumBinsCk,0,MaxBinValCk);
    Ck_pL->SetSourcePar(0,AB_pL.GetAnaSourcePar(0));
    if(SourceDescription.Contains("Mc")){
        Ck_pL->SetSourcePar(0,AB_pL.GetAnaSourcePar(0));
        Ck_pL->SetSourcePar(1,AB_pL.GetAnaSourcePar(1));
    }
    Ck_pL->SetCutOff(320,500);

    CATS AB_pXim;
    //same binning as pL, as we only use pXim as feed-down
    AB_pXim.SetMomBins(NumMomBins_pL,MomBins_pL);
    AnalysisObject.SetUpCats_pXim(AB_pXim,"pXim_HALQCD1","Gauss");
    AB_pXim.SetAnaSource(0,ResidualSourceSize);
    AB_pXim.SetNotifications(CATS::nWarning);
    AB_pXim.KillTheCat();
    DLM_Ck* Ck_pXim = new DLM_Ck(AB_pXim.GetNumSourcePars(),0,AB_pXim);
    Ck_pL->Update();
    Ck_pSigma0->Update();
    Ck_pXim->Update();

    DLM_CkDecomposition CkDec_pL("pLambda",4,*Ck_pL,hResolution_pL);
    DLM_CkDecomposition CkDec_pSigma0("pSigma0",0,*Ck_pSigma0,NULL);
    DLM_CkDecomposition CkDec_pXim("pXim",2,*Ck_pXim,NULL);

    //CkDec_pL.AddContribution(0,lam_pL[1],DLM_CkDecomposition::cFeedDown,&CkDec_pSigma0,hResidual_pL_pSigma0);
    //CkDec_pL.AddContribution(1,lam_pL[2],DLM_CkDecomposition::cFeedDown,&CkDec_pXim,hResidual_pL_pXim);
    CkDec_pL.AddContribution(0,lam_pL[1],DLM_CkDecomposition::cFeedDown);
    CkDec_pL.AddContribution(1,lam_pL[2],DLM_CkDecomposition::cFeedDown);
    CkDec_pL.AddContribution(2,lam_pL[3],DLM_CkDecomposition::cFeedDown);
    CkDec_pL.AddContribution(3,lam_pL[4],DLM_CkDecomposition::cFake);//0.03

    //for Xim we simplify a bit and take ALL feed-down as flat
    CkDec_pXim.AddContribution(0,lam_pXim[1]+lam_pXim[2]+lam_pXim[3],DLM_CkDecomposition::cFeedDown);
    CkDec_pXim.AddContribution(1,lam_pXim[4],DLM_CkDecomposition::cFake);

    CkDec_pL.Update();

    MM_PL = &CkDec_pL;
    MM_CatPL = &AB_pL;
    OutputFile->cd();
    const double FitMin = 0;
    const double FitMax = 340;
    TF1* fit_pL = new TF1("fit_pL",Fitter_MC_SimpleTemplate_pL,FitMin,FitMax,10);

    //norm
    fit_pL->SetParameter(0,fit_MC->GetParameter(0));
    fit_pL->SetParLimits(0,0.5*fit_MC->GetParameter(0),2.0*fit_MC->GetParameter(0));
//fit_pL->FixParameter(0,fit_MC->GetParameter(0));

    //strength
    fit_pL->FixParameter(1,fit_MC->GetParameter(1));
    //exp mean
    fit_pL->FixParameter(2,fit_MC->GetParameter(2));
    //exp width
    fit_pL->FixParameter(3,fit_MC->GetParameter(3));

    //rescale of the strength
    fit_pL->SetParameter(4,1);
    fit_pL->SetParLimits(4,0.1,2.0);
fit_pL->FixParameter(4,1.0);

    //cutoff
    fit_pL->FixParameter(5,500);

    //weight of NSigma
    fit_pL->SetParameter(6,0.5);
    fit_pL->SetParLimits(6,0.2,0.8);
fit_pL->FixParameter(6,0.33);

    //radius
    fit_pL->SetParameter(7,Ck_pL->GetSourcePar(0));
    fit_pL->SetParLimits(7,0.8*Ck_pL->GetSourcePar(0),1.2*Ck_pL->GetSourcePar(0));
//fit_pL->FixParameter(7,1.35);
    //alpha
//fit_pL->FixParameter(7,1.15);
fit_pL->FixParameter(7,1.2);
    fit_pL->FixParameter(8,1.3);//1.3 gives the best fit
    //fit_pL->SetParameter(8,Ck_pL->GetSourcePar(1));
    //fit_pL->SetParLimits(8,1.0,2.0);

    //some bullshit
    fit_pL->FixParameter(9,1);

    hData_pL->Fit(fit_pL,"S, N, R, M");
    printf("chi2/ndf = %.2f / %i\n",fit_pL->GetChisquare(),fit_pL->GetNDF());
    printf("prob = %.4f\n",fit_pL->GetProb());
    printf("nsigma = %.2f\n",sqrt(2)*TMath::ErfcInverse(fit_pL->GetProb()));
    fit_pL->SetNpx(1024);
/*
    TF1* fit_pL_Jet = new TF1("fit_pL_Jet","1.+[0]*TMath::Gaus(x,[1],[2],0)",0,4500);
    fit_pL_Jet->SetParameter(0,fit_pL->GetParameter(3));
    fit_pL_Jet->SetParameter(1,fit_pL->GetParameter(4));
    fit_pL_Jet->SetParameter(2,fit_pL->GetParameter(5));

    TF1* fit_pL_QS = new TF1("fit_pL_QS","(1.-[0]*exp(-pow(x*[1]/197.327,[2])))",0,4500);
    fit_pL_QS->SetParameter(0,fit_pL->GetParameter(0));
    fit_pL_QS->SetParameter(1,fit_pL->GetParameter(1));
    fit_pL_QS->SetParameter(2,fit_pL->GetParameter(2));

    TF1* fit_pL_NonFemtoNonJet = new TF1("fit_pL_NonFemtoNonJet","(1.-[0]*exp(-pow(x*[1]/197.327,[2])))*(1.+[3]*TMath::Gaus(x,[4],[5],0))*[6]*(1.+[7]*x+[8]*x*x+[9]*x*x*x)",0,4500);
    fit_pL_NonFemtoNonJet->SetParameter(0,fit_pL->GetParameter(0));
    fit_pL_NonFemtoNonJet->SetParameter(1,fit_pL->GetParameter(1));
    fit_pL_NonFemtoNonJet->SetParameter(2,fit_pL->GetParameter(2));
    fit_pL_NonFemtoNonJet->SetParameter(3,0);
    fit_pL_NonFemtoNonJet->SetParameter(4,fit_pL->GetParameter(4));
    fit_pL_NonFemtoNonJet->SetParameter(5,fit_pL->GetParameter(5));
    fit_pL_NonFemtoNonJet->SetParameter(6,fit_pL->GetParameter(6));
    fit_pL_NonFemtoNonJet->SetParameter(7,fit_pL->GetParameter(7));
    fit_pL_NonFemtoNonJet->SetParameter(8,fit_pL->GetParameter(8));
    fit_pL_NonFemtoNonJet->SetParameter(9,fit_pL->GetParameter(9));
*/
    TF1* fit_pL_NonFemto = new TF1("fit_pL_NonFemto","[0]*(1.+[1]*exp(-pow((x-[2])/(2.*[3]),2.)))",MCfit_Min,MCfit_Max);
    fit_pL_NonFemto->SetParameter(0,fit_pL->GetParameter(0));
    fit_pL_NonFemto->SetParameter(1,fit_pL->GetParameter(1)*fit_pL->GetParameter(4));
    fit_pL_NonFemto->SetParameter(2,fit_pL->GetParameter(2));
    fit_pL_NonFemto->SetParameter(3,fit_pL->GetParameter(3));

    TGraph fit_nsigma;
    fit_nsigma.SetName("fit_nsigma");

    TGraphErrors fit_sourcepars;
    fit_sourcepars.SetName("fit_sourcepars");
    fit_sourcepars.SetPoint(0,fit_pL->GetParameter(7),fit_pL->GetParameter(8));
    fit_sourcepars.SetPointError(0,fit_pL->GetParError(7),fit_pL->GetParError(8));


    TH1F* hfit_Ratio = new TH1F("hfit_Ratio","hfit_Ratio",
            hData_pL->GetNbinsX(),hData_pL->GetBinLowEdge(1),hData_pL->GetXaxis()->GetBinUpEdge(hData_pL->GetNbinsX()));

    int NumPts=0;
    double Avg_nsigma = 0;//should be around 0
    for(unsigned uBin=0; uBin<hData_pL->GetNbinsX(); uBin++){
        double MOM = hData_pL->GetBinCenter(uBin+1);
        if(MOM<FitMin || MOM>FitMax) continue;
        double CkData = hData_pL->GetBinContent(uBin+1);
        double CkErr = hData_pL->GetBinError(uBin+1);
        double CkFit = fit_pL->Eval(MOM);
        double nsigma = (CkData-CkFit)/CkErr;
        fit_nsigma.SetPoint(NumPts,MOM,nsigma);
        Avg_nsigma += nsigma;

        hfit_Ratio->SetBinContent(uBin+1,CkData/CkFit);
        hfit_Ratio->SetBinError(uBin+1,CkErr/CkFit);

        NumPts++;
    }
    Avg_nsigma /= double(NumPts);
    printf("Avg_nsigma = %.3f\n",Avg_nsigma);

    OutputFile->cd();
    fit_pL->Write();
    //fit_pL_Jet->Write();
    //fit_pL_QS->Write();
    //fit_pL_NonFemtoNonJet->Write();
    fit_pL_NonFemto->Write();
    fit_nsigma.Write();
    fit_sourcepars.Write();
    hfit_Ratio->Write();
    printf("Closing\n");
    OutputFile->Close();
    printf("Closed\n");

    //delete fit_pL;
}


double Fitter_MC_SimpleTemplate_pp(double* x, double* par){
    double& MOM = *x;
    MM_CatPP->SetAnaSource(0,par[5],true);
    if(MM_CatPP->GetNumSourcePars()>1){
        MM_CatPP->SetAnaSource(1,par[6],true);
    }
    MM_CatPP->KillTheCat();
    MM_PP->Update(true);
    static int COUNTER = 0;
    if(MOM<MM_PP->GetCk()->GetBinUpEdge(0,0)){
    COUNTER++;
    //if(COUNTER%100==0)
    printf("COUNTER=%i\n",COUNTER);
    }

    double CkVal;

    CkVal = MM_PP->EvalCk(MOM);
    //[0]*(1.+[1]*exp(-pow((x-[2])/[3],2.)))
    double Baseline = par[0]*(1.+par[1]*par[4]*exp(-pow((MOM-par[2])/(2.*par[3]),2.)));
    return CkVal*Baseline;
}


void Fit_pp_MC_SimpleTemplate(){
    const bool Include_dWaves = true;
    double ResidualSourceSize = 1.5;

    //double* MomBins_pp = NULL;
    //double* FitRegion_pp = NULL;
    //const unsigned NumBinsCk = 36;
    //const double BinWidthCk = 12;
    //const double MaxBinValCk = double(NumBinsCk)*BinWidthCk;
    unsigned NumMomBins_pp = 90;
    double kMin = 0;
    double kMax = 360;
    TString DataSample = "pp13TeV_HM_Dec19";
    //TString DataSample = "pp13TeV_HM_RotPhiDec19";

    const double MCfit_Min = 0;
    const double MCfit_Max = 600;
    TString McFile = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/MC/CF/NanoMC/CFOutput_pp_3.root";
    //TString McFile = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/MC/CF/NanoMC/CFOutput_pp_3.root";
    TString McHisto = "hCk_ReweightedMeV_0";

    //so we need to copy our histogram, as else we lose it when we delete the file
    //and we need to change to the "central" root directory, as else histoCopy will also be lost
    //and we need to play with the name a little bit, else we are fucked!
    TFile* mcfile = new TFile(McFile,"read");
    TH1F* hCk_tmp = (TH1F*)mcfile->Get(McHisto);
    if(!hCk_tmp){printf("\033[1;31mERROR:\033[0m The hCk_tmp '%s' if file '%s' does not exist\n",McHisto.Data(),McFile.Data());return;}
    TString Name = hCk_tmp->GetName();
    gROOT->cd();
    TH1F *hCk_MC = (TH1F*)hCk_tmp->Clone("hCk_MC");
    mcfile->Close();
    //delete mcfile;

    TF1* fit_MC = new TF1("fit_MC","[0]*(1.+[1]*exp(-pow((x-[2])/(2.*[3]),2.)))",MCfit_Min,MCfit_Max);
    fit_MC->SetParameter(0,1);
    fit_MC->SetParLimits(0,0.5,2.);
    fit_MC->SetParameter(1,0.03);
    fit_MC->SetParLimits(1,0.005,0.1);
    //fit_MC->SetParameter(2,0);
    //fit_MC->SetParLimits(2,-50,50);
    fit_MC->FixParameter(2,0);
    fit_MC->SetParameter(3,200);
    fit_MC->SetParLimits(3,50,500);

    hCk_MC->Fit(fit_MC,"S, N, R, M");


    //Gauss
    //TString SourceDescription = "Gauss";
    TString SourceDescription = "McGauss_ResoTM";
    //TString SourceDescription = "McLevy_ResoTM";

    DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder("/home/dmihaylov/CernBox/CatsFiles");

    CATS AB_pp;
    DLM_Ck* Ck_pp;
    AB_pp.SetMomBins(NumMomBins_pp,kMin,kMax);
    AB_pp.SetNotifications(CATS::nWarning);
    //AnalysisObject.SetUpCats_pp(AB_pp,"AV18","Gauss",0,0);//McLevyNolan_Reso
    AnalysisObject.SetUpCats_pp(AB_pp,"AV18",SourceDescription,0,202);//McLevyNolan_Reso


    TH2F* hResolution_pp = AnalysisObject.GetResolutionMatrix(DataSample,"pp");
    TH2F* hResidual_pp_pL = AnalysisObject.GetResidualMatrix("pp","pLambda");

    TH1F* hData_pp = AnalysisObject.GetAliceExpCorrFun(DataSample,"pp","_0",0,false,-1);

    //TString SourceDescription = "Gauss";
    TString OutputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/OtherTasks/Fit_pp_MC_SimpleTemplate/";
    TString OutFileName = "Out.root";
    TFile* OutputFile = new TFile(OutputFolder+OutFileName,"recreate");
    hData_pp->Write();

    double lam_pp[5];
    double lam_pL[5];

    AnalysisObject.SetUpLambdaPars_pp(DataSample,0,lam_pp);
    AnalysisObject.SetUpLambdaPars_pL(DataSample,0,0,lam_pL);

    AB_pp.SetAnaSource(0,1.3);
    if(SourceDescription.Contains("Mc")){
        AB_pp.SetAnaSource(0,1.20);
        if(SourceDescription.Contains("Levy")) AB_pp.SetAnaSource(1,1.5);
        else AB_pp.SetAnaSource(1,2.0);
    }
    //AB_pp.SetAnaSource(1,2.0);
    //AB_pp.SetNotifications(CATS::nWarning);
    //AB_pp.SetEpsilonConv(5e-8);
    //AB_pp.SetEpsilonProp(5e-8);
    AB_pp.KillTheCat();
    Ck_pp = new DLM_Ck(AB_pp.GetNumSourcePars(),0,AB_pp);
    Ck_pp->Update();

    CATS AB_pL;
    DLM_Ck* Ck_pL;
    AB_pL.SetMomBins(NumMomBins_pp,kMin,kMax);
    //AnalysisObject.SetUpCats_pL(AB_pL,"NLO_Coupled_S","Gauss");
    AnalysisObject.SetUpCats_pL(AB_pL,"Usmani","Gauss");
    AB_pL.SetAnaSource(0,ResidualSourceSize);
    AB_pL.SetNotifications(CATS::nWarning);
    AB_pL.KillTheCat();
    Ck_pL = new DLM_Ck(AB_pL.GetNumSourcePars(),0,AB_pL);
    Ck_pL->Update();

    DLM_CkDecomposition CkDec_pp("pp",3,*Ck_pp,hResolution_pp);
    DLM_CkDecomposition CkDec_pL("pLambda",2,*Ck_pL,NULL);

    CkDec_pp.AddContribution(0,lam_pp[1],DLM_CkDecomposition::cFeedDown,&CkDec_pL,hResidual_pp_pL);
    CkDec_pp.AddContribution(1,lam_pp[2],DLM_CkDecomposition::cFeedDown);
    CkDec_pp.AddContribution(2,lam_pp[3],DLM_CkDecomposition::cFake);

    CkDec_pL.AddContribution(2,lam_pL[1]+lam_pL[2]+lam_pL[3],DLM_CkDecomposition::cFeedDown);
    CkDec_pL.AddContribution(3,lam_pL[4],DLM_CkDecomposition::cFake);//0.03

    CkDec_pp.Update();
    CkDec_pL.Update();

    MM_PP = &CkDec_pp;
    MM_CatPP = &AB_pp;
    OutputFile->cd();
    const double FitMin = 0;
    const double FitMax = 320;
    TF1* fit_pp = new TF1("fit_pp",Fitter_MC_SimpleTemplate_pp,FitMin,FitMax,7);

    //norm
    fit_pp->SetParameter(0,fit_MC->GetParameter(0));
    fit_pp->SetParLimits(0,0.5*fit_MC->GetParameter(0),2.0*fit_MC->GetParameter(0));
//fit_pp->FixParameter(0,fit_MC->GetParameter(0));

    //strength
    fit_pp->FixParameter(1,fit_MC->GetParameter(1));
    //exp mean
    fit_pp->FixParameter(2,fit_MC->GetParameter(2));
    //exp width
    fit_pp->FixParameter(3,fit_MC->GetParameter(3));

    //rescale of the strength
    //fit_pp->SetParameter(4,1);
    //fit_pp->SetParLimits(4,0.1,2.0);
fit_pp->FixParameter(4,1.0);

    //radius
    fit_pp->SetParameter(5,Ck_pp->GetSourcePar(0));
    fit_pp->SetParLimits(5,0.8*Ck_pp->GetSourcePar(0),1.2*Ck_pp->GetSourcePar(0));
//fit_pp->FixParameter(7,1.35);

    //fit_pp->FixParameter(6,2.0);//1.3 gives the best fit
    fit_pp->SetParameter(6,Ck_pp->GetSourcePar(1));
    fit_pp->SetParLimits(6,1.0,2.0);
fit_pp->FixParameter(6,1.3);

    hData_pp->Fit(fit_pp,"S, N, R, M");
    printf("chi2/ndf = %.2f / %i\n",fit_pp->GetChisquare(),fit_pp->GetNDF());
    printf("prob = %.4f\n",fit_pp->GetProb());
    printf("nsigma = %.2f\n",sqrt(2)*TMath::ErfcInverse(fit_pp->GetProb()));
    fit_pp->SetNpx(1024);
/*
    TF1* fit_pL_Jet = new TF1("fit_pL_Jet","1.+[0]*TMath::Gaus(x,[1],[2],0)",0,4500);
    fit_pL_Jet->SetParameter(0,fit_pL->GetParameter(3));
    fit_pL_Jet->SetParameter(1,fit_pL->GetParameter(4));
    fit_pL_Jet->SetParameter(2,fit_pL->GetParameter(5));

    TF1* fit_pL_QS = new TF1("fit_pL_QS","(1.-[0]*exp(-pow(x*[1]/197.327,[2])))",0,4500);
    fit_pL_QS->SetParameter(0,fit_pL->GetParameter(0));
    fit_pL_QS->SetParameter(1,fit_pL->GetParameter(1));
    fit_pL_QS->SetParameter(2,fit_pL->GetParameter(2));

    TF1* fit_pL_NonFemtoNonJet = new TF1("fit_pL_NonFemtoNonJet","(1.-[0]*exp(-pow(x*[1]/197.327,[2])))*(1.+[3]*TMath::Gaus(x,[4],[5],0))*[6]*(1.+[7]*x+[8]*x*x+[9]*x*x*x)",0,4500);
    fit_pL_NonFemtoNonJet->SetParameter(0,fit_pL->GetParameter(0));
    fit_pL_NonFemtoNonJet->SetParameter(1,fit_pL->GetParameter(1));
    fit_pL_NonFemtoNonJet->SetParameter(2,fit_pL->GetParameter(2));
    fit_pL_NonFemtoNonJet->SetParameter(3,0);
    fit_pL_NonFemtoNonJet->SetParameter(4,fit_pL->GetParameter(4));
    fit_pL_NonFemtoNonJet->SetParameter(5,fit_pL->GetParameter(5));
    fit_pL_NonFemtoNonJet->SetParameter(6,fit_pL->GetParameter(6));
    fit_pL_NonFemtoNonJet->SetParameter(7,fit_pL->GetParameter(7));
    fit_pL_NonFemtoNonJet->SetParameter(8,fit_pL->GetParameter(8));
    fit_pL_NonFemtoNonJet->SetParameter(9,fit_pL->GetParameter(9));
*/
    TF1* fit_pp_NonFemto = new TF1("fit_pp_NonFemto","[0]*(1.+[1]*exp(-pow((x-[2])/(2.*[3]),2.)))",MCfit_Min,MCfit_Max);
    fit_pp_NonFemto->SetParameter(0,fit_pp->GetParameter(0));
    fit_pp_NonFemto->SetParameter(1,fit_pp->GetParameter(1)*fit_pp->GetParameter(4));
    fit_pp_NonFemto->SetParameter(2,fit_pp->GetParameter(2));
    fit_pp_NonFemto->SetParameter(3,fit_pp->GetParameter(3));

    TGraph fit_nsigma;
    fit_nsigma.SetName("fit_nsigma");

    TGraphErrors fit_sourcepars;
    fit_sourcepars.SetName("fit_sourcepars");
    fit_sourcepars.SetPoint(0,fit_pp->GetParameter(5),fit_pp->GetParameter(6));
    fit_sourcepars.SetPointError(0,fit_pp->GetParError(5),fit_pp->GetParError(6));


    TH1F* hfit_Ratio = new TH1F("hfit_Ratio","hfit_Ratio",
            hData_pp->GetNbinsX(),hData_pp->GetBinLowEdge(1),hData_pp->GetXaxis()->GetBinUpEdge(hData_pp->GetNbinsX()));

    int NumPts=0;
    double Avg_nsigma = 0;//should be around 0
    for(unsigned uBin=0; uBin<hData_pp->GetNbinsX(); uBin++){
        double MOM = hData_pp->GetBinCenter(uBin+1);
        if(MOM<FitMin || MOM>FitMax) continue;
        double CkData = hData_pp->GetBinContent(uBin+1);
        double CkErr = hData_pp->GetBinError(uBin+1);
        double CkFit = fit_pp->Eval(MOM);
        double nsigma = (CkData-CkFit)/CkErr;
        fit_nsigma.SetPoint(NumPts,MOM,nsigma);
        Avg_nsigma += nsigma;

        hfit_Ratio->SetBinContent(uBin+1,CkData/CkFit);
        hfit_Ratio->SetBinError(uBin+1,CkErr/CkFit);

        NumPts++;
    }
    Avg_nsigma /= double(NumPts);
    printf("Avg_nsigma = %.3f\n",Avg_nsigma);

    OutputFile->cd();
    fit_pp->Write();
    //fit_pL_Jet->Write();
    //fit_pL_QS->Write();
    //fit_pL_NonFemtoNonJet->Write();
    fit_pp_NonFemto->Write();
    fit_nsigma.Write();
    fit_sourcepars.Write();
    hfit_Ratio->Write();
    printf("Closing\n");
    OutputFile->Close();
    printf("Closed\n");

    //delete fit_pp;
}


double dimi_MC_template_fitter_pp(double* x, double* par){
    MM_PP->GetCk()->SetSourcePar(0,par[9]);
    if(MM_CatPP->GetNumSourcePars()>1){
        MM_PP->GetCk()->SetSourcePar(1,par[10]);
    }
    MM_PP->Update(false);
    double Com = par[0]*(1.+par[1]*exp(-pow((*x-par[2])/par[3],2.)));
    double NonCom = par[4]*(1.+par[5]*exp(-pow((*x-par[6])/par[7],2.)));
    double Femto = MM_PP->EvalCk(*x);
    return par[8]*(Com+NonCom)*Femto;
}

void Fit_pp_CommonAncestorTemplate_Ck(){

    const unsigned NumSystIter = 256;

    TString fName_MC_Common = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/pLambda_Rush/Ancestors_YieldNorm/CFOutput_pp_8_Common.root";
    TString fName_MC_NonCommon = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/pLambda_Rush/Ancestors_YieldNorm/CFOutput_pp_8_NonCommon.root";
    TString fName_MC = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/MC/CF/NanoMC/CFOutput_pp_8.root";
    TString fName_Data = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/Data/CF/Norm024034/CFOutput_pp_8.root";
    TString fName_TempInfoData = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/OtherTasks/Fit_pp_CommonAncestorTemplate_dPhi/OutputFile_pp_Data.root";
    TString fName_TempInfoMC = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/OtherTasks/Fit_pp_CommonAncestorTemplate_dPhi/OutputFile_pp_MC.root";

    TString hName_MC_Common = "hCk_FixShiftedMeV_0";
    TString hName_MC_NonCommon = "hCk_FixShiftedMeV_0";
    TString hName_MC = "hCk_FixShiftedMeV_0";
    TString hName_Data = "hCk_FixShiftedMeV_0";

    TH1D* h_MC_Common;
    TH1D* h_MC_NonCommon;
    TH1D* h_MC;
    TH1D* h_MC_TEMPLATE;
    TH1D* h_Data;

    const unsigned NumPlotBins = 200;
    const float PlotMin = 0;
    const float PlotMax = 400;
    const float PlotBinWidth = (PlotMax-PlotMin)/float(NumPlotBins);
    float* BlData_Up = new float [NumPlotBins];
    float* BlData_Down = new float [NumPlotBins];
    float* BlMc_Up = new float [NumPlotBins];
    float* BlMc_Down = new float [NumPlotBins];
    float* Data_Up = new float [NumPlotBins];
    float* Data_Down = new float [NumPlotBins];
    for(unsigned uBin=0; uBin<NumPlotBins; uBin++){
        BlData_Up[uBin] = 0;
        BlMc_Up[uBin] = 0;
        Data_Up[uBin] = 0;

        BlData_Down[uBin] = 1000;
        BlMc_Down[uBin] = 1000;
        Data_Down[uBin] = 1000;
    }

    double WC_MC_mean;// = 0.166;
    double WC_MC_err;
    //double WN_MC_mean;// = 0.834;
    //double WN_MC_err;
    double WC_MC;
    double WN_MC;
    //const double WC_Data = 0.069;
    //const double WN_Data = 0.931;
    double WC_Data_mean;// = 0.166;
    double WC_Data_err;
    //double WN_Data_mean;// = 0.834;
    //double WN_Data_err;
    double WC_Data;
    double WN_Data;
    double ddummy;
    TFile* fInputFile;
    TH1D* InputHisto;
    TGraphErrors* graph_temp;

    fInputFile = new TFile(fName_TempInfoMC,"read");
    graph_temp = (TGraphErrors*)fInputFile->Get("gWC_TotErr");
    graph_temp->GetPoint(0,ddummy,WC_MC_mean);
    WC_MC_err = graph_temp->GetErrorY(0);
    //graph_temp = (TGraphErrors*)fInputFile->Get("gWN_TotErr");
    //graph_temp->GetPoint(0,ddummy,WN_MC_mean);
    //WN_MC_err = graph_temp->GetErrorY(0);
    delete fInputFile;
//return;
    fInputFile = new TFile(fName_TempInfoData,"read");
    graph_temp = (TGraphErrors*)fInputFile->Get("gWC_TotErr");
    graph_temp->GetPoint(0,ddummy,WC_Data_mean);
    WC_Data_err = graph_temp->GetErrorY(0);
    //graph_temp = (TGraphErrors*)fInputFile->Get("gWN_TotErr");
    //graph_temp->GetPoint(0,ddummy,WN_Data_mean);
    //WN_Data_err = graph_temp->GetErrorY(0);
    delete fInputFile;

    TRandom3 rangen(11);

    TFile* OutputFile = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/OtherTasks/Fit_pp_CAT_Ck/OutputFile.root","recreate");
    double FitMin;
    double FitMax;

    TString DataSample = "pp13TeV_HM_Dec19";
    DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder("/home/dmihaylov/CernBox/CatsFiles");

    TH2F* hResolution_pp = AnalysisObject.GetResolutionMatrix(DataSample,"pp");
    TH2F* hResidual_pp_pL = AnalysisObject.GetResidualMatrix("pp","pLambda");
    double lam_pp[5];
    double lam_pL[5];
    AnalysisObject.SetUpLambdaPars_pp(DataSample,0,lam_pp);
    AnalysisObject.SetUpLambdaPars_pL(DataSample,0,0,lam_pL);
    //Gauss
    //TString SourceDescription = "Gauss";
    TString SourceDescription = "McGauss_ResoTM";
    //TString SourceDescription = "McLevy_ResoTM";
    CATS AB_pp;
    DLM_Ck* Ck_pp;
    AB_pp.SetMomBins(105,0,420);
    AB_pp.SetNotifications(CATS::nWarning);
    AnalysisObject.SetUpCats_pp(AB_pp,"AV18",SourceDescription,0,202);//McLevyNolan_Reso
    AB_pp.SetAnaSource(0,1.3);
    double InitialRadius;
    double InitialAlpha = 2.0;
    if(SourceDescription.Contains("Mc")){
        AB_pp.SetAnaSource(0,1.20);
        AB_pp.SetAnaSource(1,2.0);
        InitialAlpha = 2.0;
    }
    InitialRadius = AB_pp.GetAnaSourcePar(0);

    //AB_pp.SetAnaSource(1,2.0);
    //AB_pp.SetNotifications(CATS::nWarning);
    AB_pp.SetEpsilonConv(2e-8);
    AB_pp.SetEpsilonProp(2e-8);
    AB_pp.KillTheCat();
    Ck_pp = new DLM_Ck(AB_pp.GetNumSourcePars(),0,AB_pp);
    Ck_pp->Update();

    CATS AB_pL;
    DLM_Ck* Ck_pL;
    AB_pL.SetMomBins(105,0,420);
    //AnalysisObject.SetUpCats_pL(AB_pL,"NLO_Coupled_S","Gauss");
    AnalysisObject.SetUpCats_pL(AB_pL,"Usmani","Gauss");
    AB_pL.SetAnaSource(0,1.3);
    AB_pL.SetNotifications(CATS::nWarning);
    AB_pL.KillTheCat();
    Ck_pL = new DLM_Ck(AB_pL.GetNumSourcePars(),0,AB_pL);
    Ck_pL->Update();

    DLM_CkDecomposition CkDec_pp("pp",3,*Ck_pp,hResolution_pp);
    DLM_CkDecomposition CkDec_pL("pLambda",2,*Ck_pL,NULL);

    CkDec_pp.AddContribution(0,lam_pp[1],DLM_CkDecomposition::cFeedDown,&CkDec_pL,hResidual_pp_pL);
    CkDec_pp.AddContribution(1,lam_pp[2],DLM_CkDecomposition::cFeedDown);
    CkDec_pp.AddContribution(2,lam_pp[3],DLM_CkDecomposition::cFake);

    CkDec_pL.AddContribution(2,lam_pL[1]+lam_pL[2]+lam_pL[3],DLM_CkDecomposition::cFeedDown);
    CkDec_pL.AddContribution(3,lam_pL[4],DLM_CkDecomposition::cFake);//0.03

    CkDec_pp.Update();
    CkDec_pL.Update();

    MM_PP = &CkDec_pp;
    MM_CatPP = &AB_pp;

    TH1F* hRadius = new TH1F("hRadius","hRadius",1024,0.5,2.0);
    TH1F* hAlpha = new TH1F("hAlpha","hAlpha",1024,0.95,2.05);
    TH1F* hnSigma = new TH1F("hnSigma","hnSigma",1024,0,20.0);

    float kStarVal;
    for(unsigned uSyst=0; uSyst<NumSystIter; uSyst++){
        printf("\r\033[K uSyst=%u",uSyst);
        cout << flush;

        if(uSyst==0) WC_MC = WC_MC_mean;
        else WC_MC = rangen.Gaus(WC_MC_mean,WC_MC_err);
        WN_MC = 1. - WC_MC;

        if(uSyst==0) WC_Data = WC_Data_mean;
        else WC_Data = rangen.Gaus(WC_Data_mean,WC_Data_err);
        WN_Data = 1. - WC_Data;

        fInputFile = new TFile(fName_MC_Common,"read");
        InputHisto = (TH1D*)fInputFile->Get(hName_MC_Common);
        OutputFile->cd();
        h_MC_Common = (TH1D*)InputHisto->Clone("h_MC_Common");
        if(uSyst==0) h_MC_Common->Write();
        else{
            float Value;
            float Error;
            float RanVal;
            for(unsigned uBin=0; uBin<h_MC_Common->GetNbinsX(); uBin++){
                Value = h_MC_Common->GetBinContent(uBin+1);
                Error = h_MC_Common->GetBinError(uBin+1);
                RanVal = rangen.Gaus(Value,Error);
                h_MC_Common->SetBinContent(uBin+1,RanVal);
                h_MC_Common->SetBinError(uBin+1,Error);
            }
        }
        FitMin = 56;
        FitMax = 380;
        TF1* fTemplate_MC_Common = new TF1("fTemplate_MC_Common","[0]*(1.+[1]*exp(-pow((x-[2])/[3],2.)))",FitMin,FitMax);
        fTemplate_MC_Common->SetParameter(0,1);
        fTemplate_MC_Common->SetParameter(1,0.3);
        fTemplate_MC_Common->SetParameter(2,250);
        fTemplate_MC_Common->SetParameter(3,150);
        h_MC_Common->Fit(fTemplate_MC_Common,"Q, S, N, R, M");
        delete fInputFile;

        fInputFile = new TFile(fName_MC_NonCommon,"read");
        InputHisto = (TH1D*)fInputFile->Get(hName_MC_NonCommon);
        OutputFile->cd();
        h_MC_NonCommon = (TH1D*)InputHisto->Clone("h_MC_NonCommon");
        if(uSyst==0) h_MC_NonCommon->Write();
        else{
            float Value;
            float Error;
            float RanVal;
            for(unsigned uBin=0; uBin<h_MC_NonCommon->GetNbinsX(); uBin++){
                Value = h_MC_NonCommon->GetBinContent(uBin+1);
                Error = h_MC_NonCommon->GetBinError(uBin+1);
                RanVal = rangen.Gaus(Value,Error);
                h_MC_NonCommon->SetBinContent(uBin+1,RanVal);
                h_MC_NonCommon->SetBinError(uBin+1,Error);
            }
        }
        FitMin = 20;
        FitMax = 380;
        TF1* fTemplate_MC_NonCommon = new TF1("fTemplate_MC_NonCommon","[0]*(1.+[1]*exp(-pow((x-[2])/[3],2.)))",FitMin,FitMax);
        fTemplate_MC_NonCommon->SetParameter(0,0.9);
        fTemplate_MC_NonCommon->SetParameter(1,0.03);
        fTemplate_MC_NonCommon->FixParameter(2,0);
        fTemplate_MC_NonCommon->SetParameter(3,120);
        h_MC_NonCommon->Fit(fTemplate_MC_NonCommon,"Q, S, N, R, M");
        delete fInputFile;

        fInputFile = new TFile(fName_MC,"read");
        InputHisto = (TH1D*)fInputFile->Get(hName_MC);
        OutputFile->cd();
        h_MC = (TH1D*)InputHisto->Clone("h_MC");
        if(uSyst==0) h_MC->Write();
        else{
            float Value;
            float Error;
            float RanVal;
            for(unsigned uBin=0; uBin<h_MC->GetNbinsX(); uBin++){
                Value = h_MC->GetBinContent(uBin+1);
                Error = h_MC->GetBinError(uBin+1);
                RanVal = rangen.Gaus(Value,Error);
                h_MC->SetBinContent(uBin+1,RanVal);
                h_MC->SetBinError(uBin+1,Error);
            }
        }
        FitMin = 56;
        FitMax = 380;
        TF1* fTemplate_MC = new TF1("fTemplate_MC","[8]*([0]*(1.+[1]*exp(-pow((x-[2])/[3],2.)))+[4]*(1.+[5]*exp(-pow((x-[6])/[7],2.))))",FitMin,FitMax);
        fTemplate_MC->FixParameter(0,fTemplate_MC_Common->GetParameter(0)*WC_MC);
        fTemplate_MC->FixParameter(1,fTemplate_MC_Common->GetParameter(1));
        fTemplate_MC->FixParameter(2,fTemplate_MC_Common->GetParameter(2));
        fTemplate_MC->FixParameter(3,fTemplate_MC_Common->GetParameter(3));
        fTemplate_MC->FixParameter(4,fTemplate_MC_NonCommon->GetParameter(0)*WN_MC);
        fTemplate_MC->FixParameter(5,fTemplate_MC_NonCommon->GetParameter(1));
        fTemplate_MC->FixParameter(6,fTemplate_MC_NonCommon->GetParameter(2));
        fTemplate_MC->FixParameter(7,fTemplate_MC_NonCommon->GetParameter(3));
        fTemplate_MC->SetParameter(8,1.0178);
        h_MC->Fit(fTemplate_MC,"Q, S, N, R, M");
        delete fInputFile;

        OutputFile->cd();

        h_MC_TEMPLATE = (TH1D*)h_MC_Common->Clone("h_MC_TEMPLATE");
        h_MC_TEMPLATE->Scale(WC_MC);
        h_MC_TEMPLATE->Add(h_MC_NonCommon,WN_MC);
        h_MC_TEMPLATE->Scale(fTemplate_MC->GetParameter(8));

        for(unsigned uBin=0; uBin<NumPlotBins; uBin++){
            kStarVal = PlotMin+float(uBin)*PlotBinWidth+PlotBinWidth*0.5;
            if(BlMc_Up[uBin]<fTemplate_MC->Eval(kStarVal)) BlMc_Up[uBin]=fTemplate_MC->Eval(kStarVal);
            if(BlMc_Down[uBin]>fTemplate_MC->Eval(kStarVal)) BlMc_Down[uBin]=fTemplate_MC->Eval(kStarVal);
        }

        fInputFile = new TFile(fName_Data,"read");
        InputHisto = (TH1D*)fInputFile->Get(hName_Data);
        OutputFile->cd();
        h_Data = (TH1D*)InputHisto->Clone("h_Data");
        if(uSyst==0) h_Data->Write();
        FitMin = 0;
        FitMax = 340;
        delete fInputFile;

        TF1* fTemplate_Data = new TF1("fTemplate_Data",dimi_MC_template_fitter_pp,FitMin,FitMax,11);
        fTemplate_Data->FixParameter(0,fTemplate_MC_Common->GetParameter(0)*WC_Data);
        fTemplate_Data->FixParameter(1,fTemplate_MC_Common->GetParameter(1));
        fTemplate_Data->FixParameter(2,fTemplate_MC_Common->GetParameter(2));
        fTemplate_Data->FixParameter(3,fTemplate_MC_Common->GetParameter(3));
        fTemplate_Data->FixParameter(4,fTemplate_MC_NonCommon->GetParameter(0)*WN_Data);
        fTemplate_Data->FixParameter(5,fTemplate_MC_NonCommon->GetParameter(1));
        fTemplate_Data->FixParameter(6,fTemplate_MC_NonCommon->GetParameter(2));
        fTemplate_Data->FixParameter(7,fTemplate_MC_NonCommon->GetParameter(3));
        fTemplate_Data->SetParameter(8,1.0);
        fTemplate_Data->SetParLimits(8,0.5,2.0);
        fTemplate_Data->SetParameter(9,InitialRadius);
        fTemplate_Data->SetParLimits(9,InitialRadius*0.8,InitialRadius*1.2);
        if(AB_pp.GetNumSourcePars()>1){
            fTemplate_Data->SetParameter(10,InitialAlpha);
            fTemplate_Data->SetParLimits(10,1,2);
        }
        else fTemplate_Data->FixParameter(10,2.0);
fTemplate_Data->FixParameter(10,2.0);
        h_Data->Fit(fTemplate_Data,"Q, S, N, R, M");

        for(unsigned uBin=0; uBin<NumPlotBins; uBin++){
            kStarVal = PlotMin+float(uBin)*PlotBinWidth+PlotBinWidth*0.5;
            if(Data_Up[uBin]<fTemplate_Data->Eval(kStarVal)) Data_Up[uBin]=fTemplate_Data->Eval(kStarVal);
            if(Data_Down[uBin]>fTemplate_Data->Eval(kStarVal)) Data_Down[uBin]=fTemplate_Data->Eval(kStarVal);
        }

        TF1* fTemplate_Data_BL = new TF1("fTemplate_Data_BL","[8]*([0]*(1.+[1]*exp(-pow((x-[2])/[3],2.)))+[4]*(1.+[5]*exp(-pow((x-[6])/[7],2.))))",FitMin,FitMax);
        fTemplate_Data_BL->FixParameter(0,fTemplate_Data->GetParameter(0));
        fTemplate_Data_BL->FixParameter(1,fTemplate_Data->GetParameter(1));
        fTemplate_Data_BL->FixParameter(2,fTemplate_Data->GetParameter(2));
        fTemplate_Data_BL->FixParameter(3,fTemplate_Data->GetParameter(3));
        fTemplate_Data_BL->FixParameter(4,fTemplate_Data->GetParameter(4));
        fTemplate_Data_BL->FixParameter(5,fTemplate_Data->GetParameter(5));
        fTemplate_Data_BL->FixParameter(6,fTemplate_Data->GetParameter(6));
        fTemplate_Data_BL->FixParameter(7,fTemplate_Data->GetParameter(7));
        fTemplate_Data_BL->FixParameter(8,fTemplate_Data->GetParameter(8));

        for(unsigned uBin=0; uBin<NumPlotBins; uBin++){
            kStarVal = PlotMin+float(uBin)*PlotBinWidth+PlotBinWidth*0.5;
            if(BlData_Up[uBin]<fTemplate_Data_BL->Eval(kStarVal)) BlData_Up[uBin]=fTemplate_Data_BL->Eval(kStarVal);
            if(BlData_Down[uBin]>fTemplate_Data_BL->Eval(kStarVal)) BlData_Down[uBin]=fTemplate_Data_BL->Eval(kStarVal);
            //if(uBin==10){
            //    printf("at %f; Up = %f; Down = %f; Val = %f\n",kStarVal,BlData_Up[uBin],BlData_Down[uBin],fTemplate_Data_BL->Eval(kStarVal));
            //}
        }

        TH1F* hfit_Ratio = new TH1F("hfit_Ratio","hfit_Ratio",
                h_Data->GetNbinsX(),h_Data->GetBinLowEdge(1),h_Data->GetXaxis()->GetBinUpEdge(h_Data->GetNbinsX()));
        TGraph fit_nsigma;
        fit_nsigma.SetName("fit_nsigma");

        int NumPts=0;
        double Avg_nsigma = 0;//should be around 0
        for(unsigned uBin=0; uBin<h_Data->GetNbinsX(); uBin++){
            double MOM = h_Data->GetBinCenter(uBin+1);
            if(MOM<FitMin || MOM>FitMax) continue;
            double CkData = h_Data->GetBinContent(uBin+1);
            double CkErr = h_Data->GetBinError(uBin+1);
            double CkFit = fTemplate_Data->Eval(MOM);
            double nsigma = (CkData-CkFit)/CkErr;
            fit_nsigma.SetPoint(NumPts,MOM,nsigma);
            Avg_nsigma += nsigma;

            hfit_Ratio->SetBinContent(uBin+1,CkData/CkFit);
            hfit_Ratio->SetBinError(uBin+1,CkErr/CkFit);

            NumPts++;
        }
        Avg_nsigma /= double(NumPts);
        //printf("Avg_nsigma = %.3f\n",Avg_nsigma);

        TGraphErrors fit_sourcepars;
        fit_sourcepars.SetName("fit_sourcepars");
        fit_sourcepars.SetPoint(0,fTemplate_Data->GetParameter(9),fTemplate_Data->GetParameter(10));
        fit_sourcepars.SetPointError(0,fTemplate_Data->GetParError(9),fTemplate_Data->GetParError(10));

        if(uSyst==0){
            OutputFile->cd();
            fTemplate_MC_Common->Write();
            fTemplate_MC_NonCommon->Write();
            fTemplate_MC->Write();
            h_MC_TEMPLATE->Write();
            fTemplate_Data->Write();
            fTemplate_Data_BL->Write();
            hfit_Ratio->Write();
            fit_nsigma.Write();
            fit_sourcepars.Write();
            //printf("nsigma = %.2f\n",sqrt(2)*TMath::ErfcInverse(fTemplate_Data->GetProb()));
        }
//printf("nsigma = %.2f\n",sqrt(2)*TMath::ErfcInverse(fTemplate_Data->GetProb()));

        hRadius->Fill(fTemplate_Data->GetParameter(9));
        hAlpha->Fill(fTemplate_Data->GetParameter(10));
        hnSigma->Fill(sqrt(2)*TMath::ErfcInverse(fTemplate_Data->GetProb()));

        ///DELETES
//printf("0\n");
        delete h_MC_Common;
//printf("1\n");
        delete h_MC_NonCommon;
//printf("2\n");
        delete h_MC;
//printf("3\n");
        delete h_MC_TEMPLATE;
//printf("4\n");
        delete h_Data;
//printf("5\n");
        delete hfit_Ratio;
//printf("6\n");

        delete fTemplate_MC_Common;
//printf("7\n");
        delete fTemplate_MC_NonCommon;
//printf("8\n");
        delete fTemplate_MC;
//printf("9\n");
        delete fTemplate_Data;
//printf("10\n");
        delete fTemplate_Data_BL;
//printf("11\n");
    }

    cout << endl;

    TGraphErrors graph_DataFit;
    graph_DataFit.SetName("graph_DataFit");
    TGraphErrors graph_DataBl;
    graph_DataBl.SetName("graph_DataBl");
    TGraphErrors graph_McBl;
    graph_McBl.SetName("graph_McBl");
    for(unsigned uBin=0; uBin<NumPlotBins; uBin++){
        graph_DataFit.SetPoint(uBin,PlotMin+float(uBin)*PlotBinWidth+PlotBinWidth*0.5,(Data_Up[uBin]+Data_Down[uBin])*0.5);
        graph_DataFit.SetPointError(uBin,0,(Data_Up[uBin]-Data_Down[uBin])*0.5);

        graph_DataBl.SetPoint(uBin,PlotMin+float(uBin)*PlotBinWidth+PlotBinWidth*0.5,(BlData_Up[uBin]+BlData_Down[uBin])*0.5);
        graph_DataBl.SetPointError(uBin,0,(BlData_Up[uBin]-BlData_Down[uBin])*0.5);

        graph_McBl.SetPoint(uBin,PlotMin+float(uBin)*PlotBinWidth+PlotBinWidth*0.5,(BlMc_Up[uBin]+BlMc_Down[uBin])*0.5);
        graph_McBl.SetPointError(uBin,0,(BlMc_Up[uBin]-BlMc_Down[uBin])*0.5);
    }

    hRadius->Write();
    hAlpha->Write();
    hnSigma->Write();
    graph_DataFit.Write();
    graph_DataBl.Write();
    graph_McBl.Write();

printf("END pp ------------------------------------\n");
    delete hRadius;
    delete hAlpha;
    delete hnSigma;
    delete OutputFile;
    delete [] BlMc_Down;
    delete [] BlMc_Up;
    delete [] BlData_Down;
    delete [] BlData_Up;
    delete [] Data_Down;
    delete [] Data_Up;
}


double dimi_MC_template_fitter_pL(double* x, double* par){
    MM_PL->GetCk()->SetSourcePar(0,par[9]);
    if(MM_CatPL->GetNumSourcePars()>1){
        MM_PL->GetCk()->SetSourcePar(1,par[10]);
    }
    MM_PL->Update(false);
    double Com = par[0]*(1.+par[1]*exp(-pow((*x-par[2])/par[3],2.)));
    double NonCom = par[4]*(1.+par[5]*exp(-pow((*x-par[6])/par[7],2.)));
    double Femto = MM_PL->EvalCk(*x);
    return par[8]*(Com+NonCom)*Femto;
}
double dimi_MC_CrazyGauss_fitter_pL(double* x, double* par){
    MM_PL->GetCk()->SetSourcePar(0,par[24]);
    if(MM_CatPL->GetNumSourcePars()>1){
        MM_PL->GetCk()->SetSourcePar(1,par[25]);
    }
    MM_PL->Update(false);
    double Com = par[0]*(1.+par[1]*exp(-pow((*x-par[2])/par[3],2.)))*(1.+par[4]*exp(-pow((*x-par[5])/par[6],2.)))*(1.+par[7]*exp(-pow((*x-par[8])/par[9],2.)))*(1.+par[10]**x);
    double NonCom = par[11]*(1.+par[12]*exp(-pow((*x-par[13])/par[14],2.)))*(1.+par[15]*exp(-pow((*x-par[16])/par[17],2.)))*(1.+par[18]**x);
    double Fermi = 1.-par[22]*exp(-pow(*x*par[23]*FmToNu,2.));
    double Femto = MM_PL->EvalCk(*x);
    return par[19]*(par[20]*Com+(1.-par[20]+par[21])*NonCom)*Fermi*Femto;
}

void Fit_pL_CommonAncestorTemplate_Ck(){

    const unsigned NumSystIter = 256;

    TString fName_MC_Common = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/pLambda_Rush/Ancestors_YieldNorm/CFOutput_pL_8_Common.root";
    TString fName_MC_NonCommon = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/pLambda_Rush/Ancestors_YieldNorm/CFOutput_pL_8_NonCommon.root";
    TString fName_MC = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/MC/CF/NanoMC/CFOutput_pL_8.root";
    TString fName_Data = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/Data/CF/Norm024034/CFOutput_pL_8.root";
    TString fName_TempInfoData = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/OtherTasks/Fit_pp_CommonAncestorTemplate_dPhi/OutputFile_pp_Data.root";
    TString fName_TempInfoMC = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/OtherTasks/Fit_pp_CommonAncestorTemplate_dPhi/OutputFile_pp_MC.root";

    //TString hName_MC_Common = "hCk_RebinnedMeV_1";
    //TString hName_MC_NonCommon = "hCk_RebinnedMeV_1";
    //TString hName_MC = "hCk_RebinnedMeV_1";
    //TString hName_Data = "hCk_RebinnedMeV_1";

    TString hName_MC_Common = "hCkTotNormWeightMeV";
    TString hName_MC_NonCommon = "hCkTotNormWeightMeV";
    TString hName_MC = "hCkTotNormWeightMeV";
    TString hName_Data = "hCkTotNormWeightMeV";
    int REBIN = 4;

    const unsigned NumPlotBins = 100;
    const float PlotMin = 0;
    const float PlotMax = 400;
    const float PlotBinWidth = (PlotMax-PlotMin)/float(NumPlotBins);
    float* BlData_Up = new float [NumPlotBins];
    float* BlData_Down = new float [NumPlotBins];
    float* BlMc_Up = new float [NumPlotBins];
    float* BlMc_Down = new float [NumPlotBins];
    float* Data_Up = new float [NumPlotBins];
    float* Data_Down = new float [NumPlotBins];
    for(unsigned uBin=0; uBin<NumPlotBins; uBin++){
        BlData_Up[uBin] = 0;
        BlMc_Up[uBin] = 0;
        Data_Up[uBin] = 0;

        BlData_Down[uBin] = 1000;
        BlMc_Down[uBin] = 1000;
        Data_Down[uBin] = 1000;
    }

    TH1D* h_MC_Common;
    TH1D* h_MC_NonCommon;
    TH1D* h_MC;
    TH1D* h_MC_TEMPLATE;
    TH1D* h_Data;

    double WC_MC_mean;// = 0.166;
    double WC_MC_err;
    //double WN_MC_mean;// = 0.834;
    //double WN_MC_err;
    double WC_MC;
    double WN_MC;
    //const double WC_Data = 0.069;
    //const double WN_Data = 0.931;
    double WC_Data_mean;// = 0.166;
    double WC_Data_err;
    //double WN_Data_mean;// = 0.834;
    //double WN_Data_err;
    double WC_Data;
    double WN_Data;
    double ddummy;
    TFile* fInputFile;
    TH1D* InputHisto;
    TGraphErrors* graph_temp;

    //const double WC_MC = 0.162;
    //const double WN_MC = 0.838;
    ////const double WC_MC = 0.8;
    ////const double WN_MC = 0.2;

    //const double WC_Data = 0.100;
    //const double WN_Data = 0.900;

    fInputFile = new TFile(fName_TempInfoMC,"read");
    graph_temp = (TGraphErrors*)fInputFile->Get("gWC_TotErr");
    graph_temp->GetPoint(0,ddummy,WC_MC_mean);
    WC_MC_err = graph_temp->GetErrorY(0);
    //graph_temp = (TGraphErrors*)fInputFile->Get("gWN_TotErr");
    //graph_temp->GetPoint(0,ddummy,WN_MC_mean);
    //WN_MC_err = graph_temp->GetErrorY(0);
    delete fInputFile;
//return;
    fInputFile = new TFile(fName_TempInfoData,"read");
    graph_temp = (TGraphErrors*)fInputFile->Get("gWC_TotErr");
    graph_temp->GetPoint(0,ddummy,WC_Data_mean);
    WC_Data_err = graph_temp->GetErrorY(0);
    //graph_temp = (TGraphErrors*)fInputFile->Get("gWN_TotErr");
    //graph_temp->GetPoint(0,ddummy,WN_Data_mean);
    //WN_Data_err = graph_temp->GetErrorY(0);
    delete fInputFile;

    TRandom3 rangen(14);

    TFile* OutputFile = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/OtherTasks/Fit_pL_CAT_Ck/OutputFile.root","recreate");
    double FitMin;
    double FitMax;

    TString DataSample = "pp13TeV_HM_Dec19";
    DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder("/home/dmihaylov/CernBox/CatsFiles");
    const double ResidualSourceSize = 1.4;

    DLM_Ck* Ck_pSigma0 = new DLM_Ck(1,0,24,0,384,Lednicky_gauss_Sigma0);
    Ck_pSigma0->SetSourcePar(0,ResidualSourceSize);

    TH2F* hResolution_pL = AnalysisObject.GetResolutionMatrix(DataSample,"pLambda");
    TH2F* hResidual_pL_pSigma0 = AnalysisObject.GetResidualMatrix("pLambda","pSigma0");
    TH2F* hResidual_pL_pXim = AnalysisObject.GetResidualMatrix("pLambda","pXim");

    TH1F* hData_pL = AnalysisObject.GetAliceExpCorrFun(DataSample,"pLambda","_0",1,false,-1);

    double lam_pL[5];
    double lam_pXim[5];

    AnalysisObject.SetUpLambdaPars_pL(DataSample,0,0,lam_pL);
    AnalysisObject.SetUpLambdaPars_pXim(DataSample,0,0,lam_pXim);

    //TString SourceDescription = "Gauss";
    TString SourceDescription = "McGauss_ResoTM";
    //TString SourceDescription = "McLevy_ResoTM";

    CATS AB_pL;
    DLM_Ck* Ck_pL;
    AB_pL.SetMomBins(22,0,352);
    AnalysisObject.SetUpCats_pL(AB_pL,"NLO_Coupled_SPD",SourceDescription,0,202);//NLO_Coupled_S

    const double CUSP_WEIGHT = 0.33;//0.54
    AB_pL.SetChannelWeight(7,1./4.*CUSP_WEIGHT);//1S0 SN(s) -> LN(s)
    AB_pL.SetChannelWeight(8,3./4.*CUSP_WEIGHT);//3S1 SN(s) -> LN(s)
    AB_pL.SetChannelWeight(10,3./4.*CUSP_WEIGHT);//3S1 SN(d) -> LN(s)
    AB_pL.SetChannelWeight(13,3./20.*CUSP_WEIGHT);//3D1 SN(d) -> LN(d)
    AB_pL.SetChannelWeight(15,3./20.*CUSP_WEIGHT);//3D1 SN(s) -> LN(d)
    double InitialRadius;
    double InitialAlpha = 2.0;
    AB_pL.SetAnaSource(0,1.4);
    if(SourceDescription.Contains("Mc")){
        AB_pL.SetAnaSource(0,1.10);//c.a. 10% smaller compared to p-p due to the mT scaling
        AB_pL.SetAnaSource(1,2.0);
    }
    InitialRadius = AB_pL.GetAnaSourcePar(0);
    //AB_pL.SetNotifications(CATS::nWarning);
AB_pL.SetNotifications(CATS::nError);
    AB_pL.KillTheCat();
    Ck_pL = new DLM_Ck(AB_pL.GetNumSourcePars(),0,AB_pL,24,0,384);
    Ck_pL->SetSourcePar(0,AB_pL.GetAnaSourcePar(0));
    if(SourceDescription.Contains("Mc")){
        Ck_pL->SetSourcePar(0,AB_pL.GetAnaSourcePar(0));
        Ck_pL->SetSourcePar(1,AB_pL.GetAnaSourcePar(1));
    }
    Ck_pL->SetCutOff(350,500);

    CATS AB_pXim;
    //same binning as pL, as we only use pXim as feed-down
    AB_pXim.SetMomBins(24,0,384);
    AnalysisObject.SetUpCats_pXim(AB_pXim,"pXim_HALQCD1","Gauss");
    AB_pXim.SetAnaSource(0,ResidualSourceSize);
    AB_pXim.SetNotifications(CATS::nWarning);
    AB_pXim.KillTheCat();
    DLM_Ck* Ck_pXim = new DLM_Ck(AB_pXim.GetNumSourcePars(),0,AB_pXim);
    Ck_pL->Update();
    Ck_pSigma0->Update();
    Ck_pXim->Update();

    DLM_CkDecomposition CkDec_pL("pLambda",4,*Ck_pL,hResolution_pL);
    DLM_CkDecomposition CkDec_pSigma0("pSigma0",0,*Ck_pSigma0,NULL);
    DLM_CkDecomposition CkDec_pXim("pXim",2,*Ck_pXim,NULL);

    //CkDec_pL.AddContribution(0,lam_pL[1],DLM_CkDecomposition::cFeedDown,&CkDec_pSigma0,hResidual_pL_pSigma0);
    //CkDec_pL.AddContribution(1,lam_pL[2],DLM_CkDecomposition::cFeedDown,&CkDec_pXim,hResidual_pL_pXim);
    CkDec_pL.AddContribution(0,lam_pL[1],DLM_CkDecomposition::cFeedDown);
    CkDec_pL.AddContribution(1,lam_pL[2],DLM_CkDecomposition::cFeedDown);
    CkDec_pL.AddContribution(2,lam_pL[3],DLM_CkDecomposition::cFeedDown);
    CkDec_pL.AddContribution(3,lam_pL[4],DLM_CkDecomposition::cFake);//0.03

    //for Xim we simplify a bit and take ALL feed-down as flat
    CkDec_pXim.AddContribution(0,lam_pXim[1]+lam_pXim[2]+lam_pXim[3],DLM_CkDecomposition::cFeedDown);
    CkDec_pXim.AddContribution(1,lam_pXim[4],DLM_CkDecomposition::cFake);

    CkDec_pL.Update();

    MM_PL = &CkDec_pL;
    MM_CatPL = &AB_pL;

    TH1F* hRadius = new TH1F("hRadius","hRadius",1024,0.5,2.0);
    TH1F* hAlpha = new TH1F("hAlpha","hAlpha",1024,0.95,2.05);
    TH1F* hnSigma = new TH1F("hnSigma","hnSigma",1024,0,20.0);

    float kStarVal;
    for(unsigned uSyst=0; uSyst<NumSystIter; uSyst++){
        printf("\r\033[K uSyst=%u",uSyst);
        cout << flush;

        if(uSyst==0) WC_MC = WC_MC_mean;
        else WC_MC = rangen.Gaus(WC_MC_mean,WC_MC_err);
        WN_MC = 1. - WC_MC;

        if(uSyst==0) WC_Data = WC_Data_mean;
        else WC_Data = rangen.Gaus(WC_Data_mean,WC_Data_err);
        WN_Data = 1. - WC_Data;

        fInputFile = new TFile(fName_MC_Common,"read");
        InputHisto = (TH1D*)fInputFile->Get(hName_MC_Common);
        OutputFile->cd();
        h_MC_Common = (TH1D*)InputHisto->Clone("h_MC_Common");
        h_MC_Common->Rebin(REBIN);
        h_MC_Common->Scale(1./double(REBIN));
        if(uSyst==0) h_MC_Common->Write();
        else{
            float Value;
            float Error;
            float RanVal;
            for(unsigned uBin=0; uBin<h_MC_Common->GetNbinsX(); uBin++){
                Value = h_MC_Common->GetBinContent(uBin+1);
                Error = h_MC_Common->GetBinError(uBin+1);
                RanVal = rangen.Gaus(Value,Error);
                h_MC_Common->SetBinContent(uBin+1,RanVal);
                h_MC_Common->SetBinError(uBin+1,Error);
            }
        }
        FitMin = 0;
        FitMax = 380;
        TF1* fTemplate_MC_Common = new TF1("fTemplate_MC_Common","[0]*(1.+[1]*exp(-pow((x-[2])/[3],2.)))",FitMin,FitMax);
        fTemplate_MC_Common->SetParameter(0,1);
        fTemplate_MC_Common->SetParameter(1,0.3);
        fTemplate_MC_Common->SetParameter(2,250);
        fTemplate_MC_Common->SetParameter(3,150);
        h_MC_Common->Fit(fTemplate_MC_Common,"Q, S, N, R, M");
        delete fInputFile;

        fInputFile = new TFile(fName_MC_NonCommon,"read");
        InputHisto = (TH1D*)fInputFile->Get(hName_MC_NonCommon);
        OutputFile->cd();
        h_MC_NonCommon = (TH1D*)InputHisto->Clone("h_MC_NonCommon");
        h_MC_NonCommon->Rebin(REBIN);
        h_MC_NonCommon->Scale(1./double(REBIN));
        if(uSyst==0) h_MC_NonCommon->Write();
        else{
            float Value;
            float Error;
            float RanVal;
            for(unsigned uBin=0; uBin<h_MC_NonCommon->GetNbinsX(); uBin++){
                Value = h_MC_NonCommon->GetBinContent(uBin+1);
                Error = h_MC_NonCommon->GetBinError(uBin+1);
                RanVal = rangen.Gaus(Value,Error);
                h_MC_NonCommon->SetBinContent(uBin+1,RanVal);
                h_MC_NonCommon->SetBinError(uBin+1,Error);
            }
        }
        FitMin = 0;
        FitMax = 380;
        TF1* fTemplate_MC_NonCommon = new TF1("fTemplate_MC_NonCommon","[0]*(1.+[1]*exp(-pow((x-[2])/[3],2.)))",FitMin,FitMax);
        fTemplate_MC_NonCommon->SetParameter(0,0.9);
        fTemplate_MC_NonCommon->SetParameter(1,0.03);
        fTemplate_MC_NonCommon->FixParameter(2,0);
        fTemplate_MC_NonCommon->SetParameter(3,120);
        h_MC_NonCommon->Fit(fTemplate_MC_NonCommon,"Q, S, N, R, M");
        delete fInputFile;

        fInputFile = new TFile(fName_MC,"read");
        InputHisto = (TH1D*)fInputFile->Get(hName_MC);
        OutputFile->cd();
        h_MC = (TH1D*)InputHisto->Clone("h_MC");
        h_MC->Rebin(REBIN);
        h_MC->Scale(1./double(REBIN));
        if(uSyst==0) h_MC->Write();
        else{
            float Value;
            float Error;
            float RanVal;
            for(unsigned uBin=0; uBin<h_MC->GetNbinsX(); uBin++){
                Value = h_MC->GetBinContent(uBin+1);
                Error = h_MC->GetBinError(uBin+1);
                RanVal = rangen.Gaus(Value,Error);
                h_MC->SetBinContent(uBin+1,RanVal);
                h_MC->SetBinError(uBin+1,Error);
            }
        }
        FitMin = 0;
        FitMax = 380;
        TF1* fTemplate_MC = new TF1("fTemplate_MC","[8]*([0]*(1.+[1]*exp(-pow((x-[2])/[3],2.)))+[4]*(1.+[5]*exp(-pow((x-[6])/[7],2.))))",FitMin,FitMax);
        fTemplate_MC->FixParameter(0,fTemplate_MC_Common->GetParameter(0)*WC_MC);
        fTemplate_MC->FixParameter(1,fTemplate_MC_Common->GetParameter(1));
        fTemplate_MC->FixParameter(2,fTemplate_MC_Common->GetParameter(2));
        fTemplate_MC->FixParameter(3,fTemplate_MC_Common->GetParameter(3));
        fTemplate_MC->FixParameter(4,fTemplate_MC_NonCommon->GetParameter(0)*WN_MC);
        fTemplate_MC->FixParameter(5,fTemplate_MC_NonCommon->GetParameter(1));
        fTemplate_MC->FixParameter(6,fTemplate_MC_NonCommon->GetParameter(2));
        fTemplate_MC->FixParameter(7,fTemplate_MC_NonCommon->GetParameter(3));
        fTemplate_MC->SetParameter(8,1.0178);
        h_MC->Fit(fTemplate_MC,"Q, S, N, R, M");
        delete fInputFile;

        OutputFile->cd();
        h_MC_TEMPLATE = (TH1D*)h_MC_Common->Clone("h_MC_TEMPLATE");
        h_MC_TEMPLATE->Scale(WC_MC);
        h_MC_TEMPLATE->Add(h_MC_NonCommon,WN_MC);
        h_MC_TEMPLATE->Scale(fTemplate_MC->GetParameter(8));

        for(unsigned uBin=0; uBin<NumPlotBins; uBin++){
            kStarVal = PlotMin+float(uBin)*PlotBinWidth+PlotBinWidth*0.5;
            if(BlMc_Up[uBin]<fTemplate_MC->Eval(kStarVal)) BlMc_Up[uBin]=fTemplate_MC->Eval(kStarVal);
            if(BlMc_Down[uBin]>fTemplate_MC->Eval(kStarVal)) BlMc_Down[uBin]=fTemplate_MC->Eval(kStarVal);
        }

        fInputFile = new TFile(fName_Data,"read");
        InputHisto = (TH1D*)fInputFile->Get(hName_Data);
        OutputFile->cd();
        h_Data = (TH1D*)InputHisto->Clone("h_Data");
        h_Data->Rebin(REBIN);
        h_Data->Scale(1./double(REBIN));
        if(uSyst==0) h_Data->Write();
        FitMin = 0;
        FitMax = 380;
        delete fInputFile;

        FitMin = 0;
        FitMax = 340;
        OutputFile->cd();
        TF1* fTemplate_Data = new TF1("fTemplate_Data",dimi_MC_template_fitter_pL,FitMin,FitMax,11);
        fTemplate_Data->FixParameter(0,fTemplate_MC_Common->GetParameter(0)*WC_Data);
        fTemplate_Data->FixParameter(1,fTemplate_MC_Common->GetParameter(1));
        fTemplate_Data->FixParameter(2,fTemplate_MC_Common->GetParameter(2));
        fTemplate_Data->FixParameter(3,fTemplate_MC_Common->GetParameter(3));
        fTemplate_Data->FixParameter(4,fTemplate_MC_NonCommon->GetParameter(0)*WN_Data);
        fTemplate_Data->FixParameter(5,fTemplate_MC_NonCommon->GetParameter(1));
        fTemplate_Data->FixParameter(6,fTemplate_MC_NonCommon->GetParameter(2));
        fTemplate_Data->FixParameter(7,fTemplate_MC_NonCommon->GetParameter(3));
        fTemplate_Data->SetParameter(8,1.0);
        fTemplate_Data->SetParLimits(8,0.5,2.0);
        fTemplate_Data->SetParameter(9,InitialRadius);
        fTemplate_Data->SetParLimits(9,InitialRadius*0.8,InitialRadius*1.2);
        if(AB_pL.GetNumSourcePars()>1){
            fTemplate_Data->SetParameter(10,InitialAlpha);
            fTemplate_Data->SetParLimits(10,1,2);
        }
        else fTemplate_Data->FixParameter(10,2.0);

fTemplate_Data->FixParameter(9,1.02);
    //fTemplate_Data->FixParameter(9,1.18);
    fTemplate_Data->FixParameter(10,2.0);

        h_Data->Fit(fTemplate_Data,"Q, S, N, R, M");

        for(unsigned uBin=0; uBin<NumPlotBins; uBin++){
            kStarVal = PlotMin+float(uBin)*PlotBinWidth+PlotBinWidth*0.5;
            if(Data_Up[uBin]<fTemplate_Data->Eval(kStarVal)) Data_Up[uBin]=fTemplate_Data->Eval(kStarVal);
            if(Data_Down[uBin]>fTemplate_Data->Eval(kStarVal)) Data_Down[uBin]=fTemplate_Data->Eval(kStarVal);
        }

        TF1* fTemplate_Data_BL = new TF1("fTemplate_Data_BL","[8]*([0]*(1.+[1]*exp(-pow((x-[2])/[3],2.)))+[4]*(1.+[5]*exp(-pow((x-[6])/[7],2.))))",FitMin,FitMax);
        fTemplate_Data_BL->FixParameter(0,fTemplate_Data->GetParameter(0));
        fTemplate_Data_BL->FixParameter(1,fTemplate_Data->GetParameter(1));
        fTemplate_Data_BL->FixParameter(2,fTemplate_Data->GetParameter(2));
        fTemplate_Data_BL->FixParameter(3,fTemplate_Data->GetParameter(3));
        fTemplate_Data_BL->FixParameter(4,fTemplate_Data->GetParameter(4));
        fTemplate_Data_BL->FixParameter(5,fTemplate_Data->GetParameter(5));
        fTemplate_Data_BL->FixParameter(6,fTemplate_Data->GetParameter(6));
        fTemplate_Data_BL->FixParameter(7,fTemplate_Data->GetParameter(7));
        fTemplate_Data_BL->FixParameter(8,fTemplate_Data->GetParameter(8));

        for(unsigned uBin=0; uBin<NumPlotBins; uBin++){
            kStarVal = PlotMin+float(uBin)*PlotBinWidth+PlotBinWidth*0.5;
            if(BlData_Up[uBin]<fTemplate_Data_BL->Eval(kStarVal)) BlData_Up[uBin]=fTemplate_Data_BL->Eval(kStarVal);
            if(BlData_Down[uBin]>fTemplate_Data_BL->Eval(kStarVal)) BlData_Down[uBin]=fTemplate_Data_BL->Eval(kStarVal);
        }

        TH1F* hfit_Ratio = new TH1F("hfit_Ratio","hfit_Ratio",
                h_Data->GetNbinsX(),h_Data->GetBinLowEdge(1),h_Data->GetXaxis()->GetBinUpEdge(h_Data->GetNbinsX()));
        TGraph fit_nsigma;
        fit_nsigma.SetName("fit_nsigma");

        int NumPts=0;
        double Avg_nsigma = 0;//should be around 0
        for(unsigned uBin=0; uBin<h_Data->GetNbinsX(); uBin++){
            double MOM = h_Data->GetBinCenter(uBin+1);
            if(MOM<FitMin || MOM>FitMax) continue;
            double CkData = h_Data->GetBinContent(uBin+1);
            double CkErr = h_Data->GetBinError(uBin+1);
            double CkFit = fTemplate_Data->Eval(MOM);
            double nsigma = (CkData-CkFit)/CkErr;
            fit_nsigma.SetPoint(NumPts,MOM,nsigma);
            Avg_nsigma += nsigma;

            hfit_Ratio->SetBinContent(uBin+1,CkData/CkFit);
            hfit_Ratio->SetBinError(uBin+1,CkErr/CkFit);

            NumPts++;
        }
        Avg_nsigma /= double(NumPts);
        //printf("Avg_nsigma = %.3f\n",Avg_nsigma);

        TGraphErrors fit_sourcepars;
        fit_sourcepars.SetName("fit_sourcepars");
        fit_sourcepars.SetPoint(0,fTemplate_Data->GetParameter(9),fTemplate_Data->GetParameter(10));
        fit_sourcepars.SetPointError(0,fTemplate_Data->GetParError(9),fTemplate_Data->GetParError(10));

        if(uSyst==0){
            OutputFile->cd();
            fTemplate_MC_Common->Write();
            fTemplate_MC_NonCommon->Write();
            fTemplate_MC->Write();
            h_MC_TEMPLATE->Write();
            fTemplate_Data->Write();
            fTemplate_Data_BL->Write();
            hfit_Ratio->Write();
            fit_nsigma.Write();
            fit_sourcepars.Write();
        }

        hRadius->Fill(fTemplate_Data->GetParameter(9));
        hAlpha->Fill(fTemplate_Data->GetParameter(10));
        hnSigma->Fill(sqrt(2)*TMath::ErfcInverse(fTemplate_Data->GetProb()));

        //printf("nsigma = %.2f\n",sqrt(2)*TMath::ErfcInverse(fTemplate_Data->GetProb()));

        ///DELETES
        delete h_MC_Common;
        delete h_MC_NonCommon;
        delete h_MC;
        delete h_MC_TEMPLATE;
        delete h_Data;
        delete hfit_Ratio;
        delete fTemplate_MC_Common;
        delete fTemplate_MC_NonCommon;
        delete fTemplate_MC;
        delete fTemplate_Data;
        delete fTemplate_Data_BL;
    }
    cout << endl;

    TGraphErrors graph_DataFit;
    graph_DataFit.SetName("graph_DataFit");
    TGraphErrors graph_DataBl;
    graph_DataBl.SetName("graph_DataBl");
    TGraphErrors graph_McBl;
    graph_McBl.SetName("graph_McBl");
    for(unsigned uBin=0; uBin<NumPlotBins; uBin++){
        graph_DataFit.SetPoint(uBin,PlotMin+float(uBin)*PlotBinWidth+PlotBinWidth*0.5,(Data_Up[uBin]+Data_Down[uBin])*0.5);
        graph_DataFit.SetPointError(uBin,0,(Data_Up[uBin]-Data_Down[uBin])*0.5);

        graph_DataBl.SetPoint(uBin,PlotMin+float(uBin)*PlotBinWidth+PlotBinWidth*0.5,(BlData_Up[uBin]+BlData_Down[uBin])*0.5);
        graph_DataBl.SetPointError(uBin,0,(BlData_Up[uBin]-BlData_Down[uBin])*0.5);

        graph_McBl.SetPoint(uBin,PlotMin+float(uBin)*PlotBinWidth+PlotBinWidth*0.5,(BlMc_Up[uBin]+BlMc_Down[uBin])*0.5);
        graph_McBl.SetPointError(uBin,0,(BlMc_Up[uBin]-BlMc_Down[uBin])*0.5);
    }

    hRadius->Write();
    hAlpha->Write();
    hnSigma->Write();
    graph_DataFit.Write();
    graph_DataBl.Write();
    graph_McBl.Write();

printf("END pL ------------------------------------\n");
    delete hRadius;
    delete hAlpha;
    delete hnSigma;
    delete OutputFile;
    delete [] BlMc_Down;
    delete [] BlMc_Up;
    delete [] BlData_Down;
    delete [] BlData_Up;
    delete [] Data_Down;
    delete [] Data_Up;
}

//Last hope!
void Fit_pL_CommonAncestorTemplate_CrazyGauss(){

    const unsigned NumSystIter = 4;

    TString fName_MC_Common = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/pLambda_Rush/Ancestors_YieldNorm/CFOutput_pL_8_Common.root";
    TString fName_MC_NonCommon = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/pLambda_Rush/Ancestors_YieldNorm/CFOutput_pL_8_NonCommon.root";
    TString fName_MC = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/MC/CF/NanoMC/CFOutput_pL_8.root";
    TString fName_Data = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/Data/CF/Norm024034/CFOutput_pL_8.root";

    //TString hName_MC_Common = "hCk_RebinnedMeV_1";
    //TString hName_MC_NonCommon = "hCk_RebinnedMeV_1";
    //TString hName_MC = "hCk_RebinnedMeV_1";
    //TString hName_Data = "hCk_RebinnedMeV_1";

    TString hName_MC_Common = "hCkTotNormWeightMeV";
    TString hName_MC_NonCommon = "hCkTotNormWeightMeV";
    TString hName_MC = "hCkTotNormWeightMeV";
    TString hName_Data = "hCkTotNormWeightMeV";
    int REBIN = 4;

    const unsigned NumPlotBins = 100;
    const float PlotMin = 0;
    const float PlotMax = 400;
    const float PlotBinWidth = (PlotMax-PlotMin)/float(NumPlotBins);
    float* BlData_Up = new float [NumPlotBins];
    float* BlData_Down = new float [NumPlotBins];
    float* BlMc_Up = new float [NumPlotBins];
    float* BlMc_Down = new float [NumPlotBins];
    float* Data_Up = new float [NumPlotBins];
    float* Data_Down = new float [NumPlotBins];
    for(unsigned uBin=0; uBin<NumPlotBins; uBin++){
        BlData_Up[uBin] = 0;
        BlMc_Up[uBin] = 0;
        Data_Up[uBin] = 0;

        BlData_Down[uBin] = 1000;
        BlMc_Down[uBin] = 1000;
        Data_Down[uBin] = 1000;
    }

    TH1D* h_MC_Common;
    TH1D* h_MC_NonCommon;
    TH1D* h_MC;
    TH1D* h_MC_TEMPLATE;
    TH1D* h_Data;

    //double WC_MC_mean;// = 0.166;
    //double WC_MC_err;
    //double WN_MC_mean;// = 0.834;
    //double WN_MC_err;
    //double WC_MC;
    //double WN_MC;
    //const double WC_Data = 0.069;
    //const double WN_Data = 0.931;
    //double WC_Data_mean;// = 0.166;
    //double WC_Data_err;
    //double WN_Data_mean;// = 0.834;
    //double WN_Data_err;
    //double WC_Data;
    //double WN_Data;
    double ddummy;
    TFile* fInputFile;
    TH1D* InputHisto;

    TRandom3 rangen(14);

    TFile* OutputFile = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/OtherTasks/Fit_pL_CommonAncestorTemplate_CrazyGauss/OutputFile.root","recreate");
    double FitMin;
    double FitMax;

    TString DataSample = "pp13TeV_HM_Dec19";
    DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder("/home/dmihaylov/CernBox/CatsFiles");
    const double ResidualSourceSize = 1.4;

    DLM_Ck* Ck_pSigma0 = new DLM_Ck(1,0,24,0,384,Lednicky_gauss_Sigma0);
    Ck_pSigma0->SetSourcePar(0,ResidualSourceSize);

    TH2F* hResolution_pL = AnalysisObject.GetResolutionMatrix(DataSample,"pLambda");
    TH2F* hResidual_pL_pSigma0 = AnalysisObject.GetResidualMatrix("pLambda","pSigma0");
    TH2F* hResidual_pL_pXim = AnalysisObject.GetResidualMatrix("pLambda","pXim");

    TH1F* hData_pL = AnalysisObject.GetAliceExpCorrFun(DataSample,"pLambda","_0",1,false,-1);

    double lam_pL[5];
    double lam_pXim[5];

    AnalysisObject.SetUpLambdaPars_pL(DataSample,0,0,lam_pL);
    AnalysisObject.SetUpLambdaPars_pXim(DataSample,0,0,lam_pXim);

    TString SourceDescription = "Gauss";
    //TString SourceDescription = "McGauss_ResoTM";
    //TString SourceDescription = "McLevy_ResoTM";

    CATS AB_pL;
    DLM_Ck* Ck_pL;
    AB_pL.SetMomBins(22,0,352);
    AnalysisObject.SetUpCats_pL(AB_pL,"NLO_Coupled_SPD",SourceDescription,0,202);//NLO_Coupled_S

    const double CUSP_WEIGHT = 0.33;//0.54
    AB_pL.SetChannelWeight(7,1./4.*CUSP_WEIGHT);//1S0 SN(s) -> LN(s)
    AB_pL.SetChannelWeight(8,3./4.*CUSP_WEIGHT);//3S1 SN(s) -> LN(s)
    AB_pL.SetChannelWeight(10,3./4.*CUSP_WEIGHT);//3S1 SN(d) -> LN(s)
    AB_pL.SetChannelWeight(13,3./20.*CUSP_WEIGHT);//3D1 SN(d) -> LN(d)
    AB_pL.SetChannelWeight(15,3./20.*CUSP_WEIGHT);//3D1 SN(s) -> LN(d)
    double InitialRadius;
    double InitialAlpha = 2.0;
    AB_pL.SetAnaSource(0,1.4);
    if(SourceDescription.Contains("Mc")){
        AB_pL.SetAnaSource(0,1.10);//c.a. 10% smaller compared to p-p due to the mT scaling
        AB_pL.SetAnaSource(1,2.0);
    }
    InitialRadius = AB_pL.GetAnaSourcePar(0);
    //AB_pL.SetNotifications(CATS::nWarning);
AB_pL.SetNotifications(CATS::nError);
    AB_pL.KillTheCat();
    Ck_pL = new DLM_Ck(AB_pL.GetNumSourcePars(),0,AB_pL,24,0,384);
    Ck_pL->SetSourcePar(0,AB_pL.GetAnaSourcePar(0));
    if(SourceDescription.Contains("Mc")){
        Ck_pL->SetSourcePar(0,AB_pL.GetAnaSourcePar(0));
        Ck_pL->SetSourcePar(1,AB_pL.GetAnaSourcePar(1));
    }
    Ck_pL->SetCutOff(350,500);

    CATS AB_pXim;
    //same binning as pL, as we only use pXim as feed-down
    AB_pXim.SetMomBins(24,0,384);
    AnalysisObject.SetUpCats_pXim(AB_pXim,"pXim_HALQCD1","Gauss");
    AB_pXim.SetAnaSource(0,ResidualSourceSize);
    AB_pXim.SetNotifications(CATS::nWarning);
    AB_pXim.KillTheCat();
    DLM_Ck* Ck_pXim = new DLM_Ck(AB_pXim.GetNumSourcePars(),0,AB_pXim);
    Ck_pL->Update();
    Ck_pSigma0->Update();
    Ck_pXim->Update();

    DLM_CkDecomposition CkDec_pL("pLambda",4,*Ck_pL,hResolution_pL);
    DLM_CkDecomposition CkDec_pSigma0("pSigma0",0,*Ck_pSigma0,NULL);
    DLM_CkDecomposition CkDec_pXim("pXim",2,*Ck_pXim,NULL);

    //CkDec_pL.AddContribution(0,lam_pL[1],DLM_CkDecomposition::cFeedDown,&CkDec_pSigma0,hResidual_pL_pSigma0);
    //CkDec_pL.AddContribution(1,lam_pL[2],DLM_CkDecomposition::cFeedDown,&CkDec_pXim,hResidual_pL_pXim);
    CkDec_pL.AddContribution(0,lam_pL[1],DLM_CkDecomposition::cFeedDown);
    CkDec_pL.AddContribution(1,lam_pL[2],DLM_CkDecomposition::cFeedDown);
    CkDec_pL.AddContribution(2,lam_pL[3],DLM_CkDecomposition::cFeedDown);
    CkDec_pL.AddContribution(3,lam_pL[4],DLM_CkDecomposition::cFake);//0.03

    //for Xim we simplify a bit and take ALL feed-down as flat
    CkDec_pXim.AddContribution(0,lam_pXim[1]+lam_pXim[2]+lam_pXim[3],DLM_CkDecomposition::cFeedDown);
    CkDec_pXim.AddContribution(1,lam_pXim[4],DLM_CkDecomposition::cFake);

    CkDec_pL.Update();

    MM_PL = &CkDec_pL;
    MM_CatPL = &AB_pL;

    TH1F* hRadius = new TH1F("hRadius","hRadius",1024,0.5,2.0);
    TH1F* hAlpha = new TH1F("hAlpha","hAlpha",1024,0.95,2.05);
    TH1F* hnSigma = new TH1F("hnSigma","hnSigma",1024,0,20.0);
    TH1F* hWC = new TH1F("hWC","hWC",512,-1,1);
    TH1F* hWNC = new TH1F("hWNC","hWNC",512,-1,1);

    float kStarVal;
    for(unsigned uSyst=0; uSyst<NumSystIter; uSyst++){
        printf("\r\033[K uSyst=%u",uSyst);
        cout << flush;

        fInputFile = new TFile(fName_MC_Common,"read");
        InputHisto = (TH1D*)fInputFile->Get(hName_MC_Common);
        OutputFile->cd();
        h_MC_Common = (TH1D*)InputHisto->Clone("h_MC_Common");
        h_MC_Common->Rebin(REBIN);
        h_MC_Common->Scale(1./double(REBIN));
        if(uSyst==0) h_MC_Common->Write();
        else{
            float Value;
            float Error;
            float RanVal;
            for(unsigned uBin=0; uBin<h_MC_Common->GetNbinsX(); uBin++){
                Value = h_MC_Common->GetBinContent(uBin+1);
                Error = h_MC_Common->GetBinError(uBin+1);
                RanVal = rangen.Gaus(Value,Error);
                h_MC_Common->SetBinContent(uBin+1,RanVal);
                h_MC_Common->SetBinError(uBin+1,Error);
            }
        }
        FitMin = 0;
        FitMax = 2500;
        TF1* fTemplate_MC_Common = new TF1("fTemplate_MC_Common","[0]*(1.+[1]*exp(-pow((x-[2])/[3],2.)))*(1.+[4]*exp(-pow((x-[5])/[6],2.)))*(1.+[7]*exp(-pow((x-[8])/[9],2.)))*(1.+[10]*x)",FitMin,FitMax);
        fTemplate_MC_Common->SetParameter(0,0.96);
        fTemplate_MC_Common->SetParameter(1,0.47);
        fTemplate_MC_Common->SetParameter(2,200.);
        fTemplate_MC_Common->SetParameter(3,500.);
        fTemplate_MC_Common->FixParameter(4,0);
        fTemplate_MC_Common->FixParameter(5,0.0);
        fTemplate_MC_Common->FixParameter(6,946);
        fTemplate_MC_Common->FixParameter(7,0);
        fTemplate_MC_Common->FixParameter(8,600);
        fTemplate_MC_Common->FixParameter(9,100);
        fTemplate_MC_Common->SetParameter(10,-6e-5);
        h_MC_Common->Fit(fTemplate_MC_Common,"Q, S, N, R, M");
        delete fInputFile;

        fInputFile = new TFile(fName_MC_NonCommon,"read");
        InputHisto = (TH1D*)fInputFile->Get(hName_MC_NonCommon);
        OutputFile->cd();
        h_MC_NonCommon = (TH1D*)InputHisto->Clone("h_MC_NonCommon");
        h_MC_NonCommon->Rebin(REBIN);
        h_MC_NonCommon->Scale(1./double(REBIN));
        if(uSyst==0) h_MC_NonCommon->Write();
        else{
            float Value;
            float Error;
            float RanVal;
            for(unsigned uBin=0; uBin<h_MC_NonCommon->GetNbinsX(); uBin++){
                Value = h_MC_NonCommon->GetBinContent(uBin+1);
                Error = h_MC_NonCommon->GetBinError(uBin+1);
                RanVal = rangen.Gaus(Value,Error);
                h_MC_NonCommon->SetBinContent(uBin+1,RanVal);
                h_MC_NonCommon->SetBinError(uBin+1,Error);
            }
        }
        FitMin = 0;
        FitMax = 2500;
        TF1* fTemplate_MC_NonCommon = new TF1("fTemplate_MC_NonCommon","[0]*(1.+[1]*exp(-pow((x-[2])/[3],2.)))*(1.+[4]*exp(-pow((x-[5])/[6],2.)))*(1.+[7]*x)",FitMin,FitMax);
        fTemplate_MC_NonCommon->SetParameter(0,0.9);
        fTemplate_MC_NonCommon->SetParameter(1,-0.02);
        fTemplate_MC_NonCommon->SetParameter(2,390);
        fTemplate_MC_NonCommon->SetParameter(3,185);
        fTemplate_MC_NonCommon->SetParameter(4,0.1);
        fTemplate_MC_NonCommon->SetParameter(5,-168);
        fTemplate_MC_NonCommon->SetParameter(6,216);
        fTemplate_MC_NonCommon->SetParameter(7,0.0001);
        h_MC_NonCommon->Fit(fTemplate_MC_NonCommon,"Q, S, N, R, M");
        delete fInputFile;

        fInputFile = new TFile(fName_MC,"read");
        InputHisto = (TH1D*)fInputFile->Get(hName_MC);
        OutputFile->cd();
        h_MC = (TH1D*)InputHisto->Clone("h_MC");
        h_MC->Rebin(REBIN);
        h_MC->Scale(1./double(REBIN));
        if(uSyst==0) h_MC->Write();
        else{
            float Value;
            float Error;
            float RanVal;
            for(unsigned uBin=0; uBin<h_MC->GetNbinsX(); uBin++){
                Value = h_MC->GetBinContent(uBin+1);
                Error = h_MC->GetBinError(uBin+1);
                RanVal = rangen.Gaus(Value,Error);
                h_MC->SetBinContent(uBin+1,RanVal);
                h_MC->SetBinError(uBin+1,Error);
            }
        }
        FitMin = 0;
        FitMax = 2500;
        TF1* fTemplate_MC = new TF1("fTemplate_MC","[19]*([20]*[0]*(1.+[1]*exp(-pow((x-[2])/[3],2.)))*(1.+[4]*exp(-pow((x-[5])/[6],2.)))*(1.+[7]*exp(-pow((x-[8])/[9],2.)))*(1.+[10]*x)"
                                    "+(1-[20]+[21])*[11]*(1.+[12]*exp(-pow((x-[13])/[14],2.)))*(1.+[15]*exp(-pow((x-[16])/[17],2.)))*(1.+[18]*x))",FitMin,FitMax);
        for(unsigned uPar=0; uPar<=10; uPar++){
            fTemplate_MC->FixParameter(uPar,fTemplate_MC_Common->GetParameter(uPar));
        }
        for(unsigned uPar=11; uPar<=18; uPar++){
            fTemplate_MC->FixParameter(uPar,fTemplate_MC_NonCommon->GetParameter(uPar-11));
        }
        fTemplate_MC->SetParameter(19,1);
        fTemplate_MC->SetParameter(20,0.1);
        fTemplate_MC->FixParameter(21,0.0);
        h_MC->Fit(fTemplate_MC,"Q, S, N, R, M");
        delete fInputFile;

        OutputFile->cd();
        h_MC_TEMPLATE = (TH1D*)h_MC_Common->Clone("h_MC_TEMPLATE");
        h_MC_TEMPLATE->Scale(fTemplate_MC->GetParameter(20));
        h_MC_TEMPLATE->Add(h_MC_NonCommon,1.-fTemplate_MC->GetParameter(20));
        h_MC_TEMPLATE->Scale(fTemplate_MC->GetParameter(19));

        for(unsigned uBin=0; uBin<NumPlotBins; uBin++){
            kStarVal = PlotMin+float(uBin)*PlotBinWidth+PlotBinWidth*0.5;
            if(BlMc_Up[uBin]<fTemplate_MC->Eval(kStarVal)) BlMc_Up[uBin]=fTemplate_MC->Eval(kStarVal);
            if(BlMc_Down[uBin]>fTemplate_MC->Eval(kStarVal)) BlMc_Down[uBin]=fTemplate_MC->Eval(kStarVal);
        }

        fInputFile = new TFile(fName_Data,"read");
        InputHisto = (TH1D*)fInputFile->Get(hName_Data);
        OutputFile->cd();
        h_Data = (TH1D*)InputHisto->Clone("h_Data");
        h_Data->Rebin(REBIN);
        h_Data->Scale(1./double(REBIN));
        if(uSyst==0) h_Data->Write();
        delete fInputFile;

        FitMin = 0;
        FitMax = 340;//340
        OutputFile->cd();
        TF1* fTemplate_Data = new TF1("fTemplate_Data",dimi_MC_CrazyGauss_fitter_pL,FitMin,FitMax,26);

        for(unsigned uPar=0; uPar<=9; uPar++){
            fTemplate_Data->FixParameter(uPar,fTemplate_MC_Common->GetParameter(uPar));
        }
        for(unsigned uPar=11; uPar<=17; uPar++){
            fTemplate_Data->FixParameter(uPar,fTemplate_MC_NonCommon->GetParameter(uPar-11));
        }
        fTemplate_Data->SetParameter(10,fTemplate_MC_Common->GetParameter(10));
        fTemplate_Data->SetParLimits(10,fTemplate_MC_Common->GetParameter(10)-0.25*fabs(fTemplate_MC_Common->GetParameter(10)),fTemplate_MC_Common->GetParameter(10)+.25*fabs(fTemplate_MC_Common->GetParameter(10)));
        fTemplate_Data->SetParameter(18,fTemplate_MC_NonCommon->GetParameter(7));
        fTemplate_Data->SetParLimits(18,fTemplate_MC_NonCommon->GetParameter(7)-0.25*fabs(fTemplate_MC_NonCommon->GetParameter(7)),fTemplate_MC_NonCommon->GetParameter(7)+0.25*fabs(fTemplate_MC_NonCommon->GetParameter(7)));

/*
        for(unsigned uPar=0; uPar<=10; uPar++){
            fTemplate_Data->SetParameter(uPar,fTemplate_MC_Common->GetParameter(uPar));
            fTemplate_Data->SetParLimits(uPar,
                                         fTemplate_MC_Common->GetParameter(uPar)-0.25*fabs(fTemplate_MC_Common->GetParameter(uPar)),
                                         fTemplate_MC_Common->GetParameter(uPar)+0.25*fabs(fTemplate_MC_Common->GetParameter(uPar)));
        }
        for(unsigned uPar=11; uPar<=18; uPar++){
            fTemplate_Data->SetParameter(uPar,fTemplate_MC_Common->GetParameter(uPar-11));
            fTemplate_Data->SetParLimits(uPar,
                                         fTemplate_MC_NonCommon->GetParameter(uPar-11)-0.25*fabs(fTemplate_MC_NonCommon->GetParameter(uPar-11)),
                                         fTemplate_MC_NonCommon->GetParameter(uPar-11)+0.25*fabs(fTemplate_MC_NonCommon->GetParameter(uPar-11)));
        }
*/
        fTemplate_Data->SetParameter(19,1);
        fTemplate_Data->SetParameter(20,0.1);
        fTemplate_Data->SetParLimits(20,0.0,0.2);
        fTemplate_Data->SetParameter(21,0.0);
        fTemplate_Data->SetParLimits(21,-0.5,0.5);
        fTemplate_Data->SetParameter(22,0.2);
fTemplate_Data->FixParameter(22,0.0);
        fTemplate_Data->SetParameter(23,0.4);
fTemplate_Data->FixParameter(23,0.4);
        fTemplate_Data->SetParameter(24,InitialRadius);
        fTemplate_Data->SetParLimits(24,InitialRadius*0.8,InitialRadius*1.2);
        if(AB_pL.GetNumSourcePars()>1){
            fTemplate_Data->SetParameter(25,InitialAlpha);
            fTemplate_Data->SetParLimits(25,1,2);
        }
        else fTemplate_Data->FixParameter(25,2.0);

    //fTemplate_Data->FixParameter(9,1.08);
    //fTemplate_Data->FixParameter(9,1.18);
    fTemplate_Data->FixParameter(25,2.0);

        if(uSyst==0) h_Data->Fit(fTemplate_Data,"S, N, R, M");
        else h_Data->Fit(fTemplate_Data,"Q, S, N, R, M");

        for(unsigned uBin=0; uBin<NumPlotBins; uBin++){
            kStarVal = PlotMin+float(uBin)*PlotBinWidth+PlotBinWidth*0.5;
            if(Data_Up[uBin]<fTemplate_Data->Eval(kStarVal)) Data_Up[uBin]=fTemplate_Data->Eval(kStarVal);
            if(Data_Down[uBin]>fTemplate_Data->Eval(kStarVal)) Data_Down[uBin]=fTemplate_Data->Eval(kStarVal);
        }

        TF1* fTemplate_Data_BL = new TF1("fTemplate_Data_BL","[19]*([20]*[0]*(1.+[1]*exp(-pow((x-[2])/[3],2.)))*(1.+[4]*exp(-pow((x-[5])/[6],2.)))*(1.+[7]*exp(-pow((x-[8])/[9],2.)))*(1.+[10]*x)"
                                    "+(1-[20]+[21])*[11]*(1.+[12]*exp(-pow((x-[13])/[14],2.)))*(1.+[15]*exp(-pow((x-[16])/[17],2.)))*(1.+[18]*x))",FitMin,FitMax);
        for(unsigned uPar=0; uPar<=21; uPar++){
            fTemplate_Data_BL->SetParameter(uPar,fTemplate_Data->GetParameter(uPar));
        }

        for(unsigned uBin=0; uBin<NumPlotBins; uBin++){
            kStarVal = PlotMin+float(uBin)*PlotBinWidth+PlotBinWidth*0.5;
            if(BlData_Up[uBin]<fTemplate_Data_BL->Eval(kStarVal)) BlData_Up[uBin]=fTemplate_Data_BL->Eval(kStarVal);
            if(BlData_Down[uBin]>fTemplate_Data_BL->Eval(kStarVal)) BlData_Down[uBin]=fTemplate_Data_BL->Eval(kStarVal);
        }

        TH1F* hfit_Ratio = new TH1F("hfit_Ratio","hfit_Ratio",
                h_Data->GetNbinsX(),h_Data->GetBinLowEdge(1),h_Data->GetXaxis()->GetBinUpEdge(h_Data->GetNbinsX()));
        TGraph fit_nsigma;
        fit_nsigma.SetName("fit_nsigma");

        int NumPts=0;
        double Avg_nsigma = 0;//should be around 0
        for(unsigned uBin=0; uBin<h_Data->GetNbinsX(); uBin++){
            double MOM = h_Data->GetBinCenter(uBin+1);
            if(MOM<FitMin || MOM>FitMax) continue;
            double CkData = h_Data->GetBinContent(uBin+1);
            double CkErr = h_Data->GetBinError(uBin+1);
            double CkFit = fTemplate_Data->Eval(MOM);
            double nsigma = (CkData-CkFit)/CkErr;
            fit_nsigma.SetPoint(NumPts,MOM,nsigma);
            Avg_nsigma += nsigma;

            hfit_Ratio->SetBinContent(uBin+1,CkData/CkFit);
            hfit_Ratio->SetBinError(uBin+1,CkErr/CkFit);

            NumPts++;
        }
        Avg_nsigma /= double(NumPts);
        //printf("Avg_nsigma = %.3f\n",Avg_nsigma);

        TGraphErrors fit_sourcepars;
        fit_sourcepars.SetName("fit_sourcepars");
        fit_sourcepars.SetPoint(0,fTemplate_Data->GetParameter(24),fTemplate_Data->GetParameter(25));
        fit_sourcepars.SetPointError(0,fTemplate_Data->GetParError(24),fTemplate_Data->GetParError(25));

        if(uSyst==0){
            OutputFile->cd();
            fTemplate_MC_Common->Write();
            fTemplate_MC_NonCommon->Write();
            fTemplate_MC->Write();
            h_MC_TEMPLATE->Write();
            fTemplate_Data->Write();
            fTemplate_Data_BL->Write();
            hfit_Ratio->Write();
            fit_nsigma.Write();
            fit_sourcepars.Write();
        }

        hRadius->Fill(fTemplate_Data->GetParameter(24));
        hAlpha->Fill(fTemplate_Data->GetParameter(25));
        hnSigma->Fill(sqrt(2)*TMath::ErfcInverse(fTemplate_Data->GetProb()));
        hWC->Fill(fTemplate_Data->GetParameter(20));
        hWNC->Fill(1.-fTemplate_Data->GetParameter(20)+fTemplate_Data->GetParameter(21));

        //printf("nsigma = %.2f\n",sqrt(2)*TMath::ErfcInverse(fTemplate_Data->GetProb()));

        ///DELETES
        delete h_MC_Common;
        delete h_MC_NonCommon;
        delete h_MC;
        delete h_MC_TEMPLATE;
        delete h_Data;
        delete hfit_Ratio;
        delete fTemplate_MC_Common;
        delete fTemplate_MC_NonCommon;
        delete fTemplate_MC;
        delete fTemplate_Data;
        delete fTemplate_Data_BL;
    }
    cout << endl;

    TGraphErrors graph_DataFit;
    graph_DataFit.SetName("graph_DataFit");
    TGraphErrors graph_DataBl;
    graph_DataBl.SetName("graph_DataBl");
    TGraphErrors graph_McBl;
    graph_McBl.SetName("graph_McBl");
    for(unsigned uBin=0; uBin<NumPlotBins; uBin++){
        graph_DataFit.SetPoint(uBin,PlotMin+float(uBin)*PlotBinWidth+PlotBinWidth*0.5,(Data_Up[uBin]+Data_Down[uBin])*0.5);
        graph_DataFit.SetPointError(uBin,0,(Data_Up[uBin]-Data_Down[uBin])*0.5);

        graph_DataBl.SetPoint(uBin,PlotMin+float(uBin)*PlotBinWidth+PlotBinWidth*0.5,(BlData_Up[uBin]+BlData_Down[uBin])*0.5);
        graph_DataBl.SetPointError(uBin,0,(BlData_Up[uBin]-BlData_Down[uBin])*0.5);

        graph_McBl.SetPoint(uBin,PlotMin+float(uBin)*PlotBinWidth+PlotBinWidth*0.5,(BlMc_Up[uBin]+BlMc_Down[uBin])*0.5);
        graph_McBl.SetPointError(uBin,0,(BlMc_Up[uBin]-BlMc_Down[uBin])*0.5);
    }

    hRadius->Write();
    hAlpha->Write();
    hnSigma->Write();
    hWC->Write();
    hWNC->Write();
    graph_DataFit.Write();
    graph_DataBl.Write();
    graph_McBl.Write();

printf("END pL ------------------------------------\n");
    delete hRadius;
    delete hAlpha;
    delete hnSigma;
    delete hWC;
    delete hWNC;
    delete OutputFile;
    delete [] BlMc_Down;
    delete [] BlMc_Up;
    delete [] BlData_Down;
    delete [] BlData_Up;
    delete [] Data_Down;
    delete [] Data_Up;
}


void Fit_pp_CommonAncestorTemplate_dPhi(const bool useMC){

    const unsigned NumSystIter = 256;
    TH1F* hWeight_Common = new TH1F("hWeight_Common","hWeight_Common",4096,0,0.5);
    TGraphErrors gWC_Def;//with stat uncertainty
    gWC_Def.SetName("gWC_Def");
    gWC_Def.SetMarkerStyle(20);
    gWC_Def.SetMarkerSize(1);
    TGraphErrors gWC_Syst;
    gWC_Syst.SetName("gWC_Syst");
    gWC_Syst.SetMarkerStyle(20);
    gWC_Syst.SetMarkerSize(1);
    TGraphErrors gWC_TotErr;
    gWC_TotErr.SetName("gWC_TotErr");
    gWC_TotErr.SetMarkerStyle(20);
    gWC_TotErr.SetMarkerSize(1);

    TString fName_MC_Common = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/pLambda_Rush/Ancestors/dEtadPhiAncestors.root";
    TString fName_MC_NonCommon = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/pLambda_Rush/Ancestors/dEtadPhiAncestors.root";
    TString fName_Data;
    if(useMC) fName_Data = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/pLambda_Rush/dEtadPhi_Data/MC/dEtadPhi.root";
    else fName_Data = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/pLambda_Rush/dEtadPhi_Data/Data/dEtadPhi.root";
    TString l1Name_MC_Common = "ProtonProton";
    TString l1Name_MC_NonCommon = "ProtonProton";
    TString l1Name_Data = "ProtonProton";
    TString l2Name_MC_Common = "ppCommon";
    TString l2Name_MC_NonCommon = "ppNonCommon";
    TString l2Name_Data = "ppIntegrated";
    TString hName_MC_Common = "SEdPhidEtaDistCommon_Particle0_Particle0ShiftedProjection";
    TString hName_MC_NonCommon = "SEdPhidEtaDistNonCommon_Particle0_Particle0AShiftedProjection";
    TString hName_Data = "SEdPhidEtaDist_Particle0_Particle0ShiftedProjection";

    TFile* fInputFile;
    TList* List1;
    TList* List2;
    TH1D* InputHisto;

    TFile* OutputFile = new TFile(
    TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/OtherTasks/Fit_pp_CommonAncestorTemplate_dPhi/OutputFile_pp_%s.root",
                    useMC?"MC":"Data"),"recreate");    double FitMin;
    double FitMax;
    double StatErr;

    TRandom3 rangen(11);

    for(unsigned uSyst=0; uSyst<NumSystIter; uSyst++){
        printf("\r\033[K uSyst=%u",uSyst);
        cout << flush;
        fInputFile = new TFile(fName_MC_Common,"read");
        List1 = (TList*)fInputFile->Get(l1Name_MC_Common);
        List2 = (TList*)List1->FindObject(l2Name_MC_Common);
        InputHisto = (TH1D*)List2->FindObject(hName_MC_Common);
        FitMin = -1.8;
        FitMax = 4.3;//4.3
        OutputFile->cd();
        if(uSyst==0) InputHisto->Write();
        else{
            float Value;
            float Error;
            float RanVal;
            for(unsigned uBin=0; uBin<InputHisto->GetNbinsX(); uBin++){
                Value = InputHisto->GetBinContent(uBin+1);
                Error = InputHisto->GetBinError(uBin+1);
                RanVal = rangen.Gaus(Value,Error);
                InputHisto->SetBinContent(uBin+1,RanVal);
                InputHisto->SetBinError(uBin+1,Error);
            }
        }
        //TF1* fTemplate_MC_Common = new TF1("fTemplate_MC_Common","[0]*(1.+[1]*sin([2]*x+[3])+[4]*exp(-pow((x-[5])/[6],2.)))",FitMin,FitMax);
        //TF1* fTemplate_MC_Common = new TF1("fTemplate_MC_Common",
        //"[0]*(1.+[1]*exp(-pow((x-[2])/[3],2.)))*(1.+[4]*x+[5]*x*x+[6]*pow(x,3.)+[7]*pow(x,4.)+[8]*pow(x,5.)+[9]*pow(x,6.)+[10]*pow(x,7.))",
        //FitMin,FitMax);
        /*
        fTemplate_MC_Common->SetParameter(0,1);
        fTemplate_MC_Common->SetParLimits(0,0.5,1.5);
        fTemplate_MC_Common->SetParameter(1,0.1);
        fTemplate_MC_Common->SetParLimits(1,0.033,0.3);
        fTemplate_MC_Common->SetParameter(2,1);
        fTemplate_MC_Common->SetParLimits(2,0.33,3);
        fTemplate_MC_Common->SetParameter(3,-1.4);
        fTemplate_MC_Common->SetParLimits(3,-3,3);
        fTemplate_MC_Common->SetParameter(4,1);
        fTemplate_MC_Common->SetParLimits(4,0.33,3);
        fTemplate_MC_Common->FixParameter(5,0);
        fTemplate_MC_Common->SetParameter(6,0.4);
        fTemplate_MC_Common->SetParLimits(6,0.1,1.6);
        */
        //fTemplate_MC_Common->SetParameter(0,1);
        //fTemplate_MC_Common->SetParameter(1,1);
        //fTemplate_MC_Common->FixParameter(2,0);
        //fTemplate_MC_Common->SetParameter(3,1);


        TF1* fTemplate_MC_Common = new TF1("fTemplate_MC_Common","[0]+[1]*exp(-pow((fabs(x)-[2])/[3],2.))+[4]*exp(-pow((fabs(x)-[5])/[6],2.))",FitMin,FitMax);
        //TF1* fTemplate_MC_Common = new TF1("fTemplate_MC_Common","[0]+[0]*[1]*exp(-pow((x-[2])/[3],2.))+[0]*[4]*exp(-pow((x-[5])/[6],2.))",FitMin,FitMax);
        //TF1* fTemplate_MC_Common = new TF1("fTemplate_MC_Common","[0]*(1+[1]*exp(-pow((x-[2])/[3],2.)))*(1+[4]*exp(-pow((x-[5])/[6],2.)))",FitMin,FitMax);
        //TF1* fTemplate_MC_Common = new TF1("fTemplate_MC_Common","[0]*([7]+[1]*exp(-pow((x-[2])/[3],2.))+[4]*exp(-pow((x-[5])/[6],2.)))",FitMin,FitMax);
        fTemplate_MC_Common->SetParameter(0,1);
        fTemplate_MC_Common->SetParameter(1,1);
        fTemplate_MC_Common->FixParameter(2,0);
        fTemplate_MC_Common->SetParameter(3,1);
        fTemplate_MC_Common->SetParameter(4,-0.1);
        fTemplate_MC_Common->FixParameter(5,TMath::Pi());
        fTemplate_MC_Common->SetParameter(6,1);


        //fTemplate_MC_Common->SetParameter(7,1);

        InputHisto->Fit(fTemplate_MC_Common,"Q, S, N, R, M");
        delete fInputFile;

        fInputFile = new TFile(fName_MC_NonCommon,"read");
        List1 = (TList*)fInputFile->Get(l1Name_MC_NonCommon);
        List2 = (TList*)List1->FindObject(l2Name_MC_NonCommon);
        InputHisto = (TH1D*)List2->FindObject(hName_MC_NonCommon);
        //FitMin = InputHisto->GetBinLowEdge(1);
        //FitMax = InputHisto->GetBinLowEdge(InputHisto->GetNbinsX());
        FitMin = -1.8;
        FitMax = 4.3;
        OutputFile->cd();
        if(uSyst==0) InputHisto->Write();
        else{
            float Value;
            float Error;
            float RanVal;
            for(unsigned uBin=0; uBin<InputHisto->GetNbinsX(); uBin++){
                Value = InputHisto->GetBinContent(uBin+1);
                Error = InputHisto->GetBinError(uBin+1);
                RanVal = rangen.Gaus(Value,Error);
                InputHisto->SetBinContent(uBin+1,RanVal);
                InputHisto->SetBinError(uBin+1,Error);
            }
        }
        //TF1* fTemplate_MC_NonCommon = new TF1("fTemplate_MC_NonCommon","[0]*(1.+[1]*sin([2]*x+[3]))",FitMin,FitMax);
        //TF1* fTemplate_MC_NonCommon = new TF1("fTemplate_MC_NonCommon","[0]*(1.+[1]*sin([2]*x+[3]))*(1.+[4]*x+[5]*x*x+[6]*pow(x,3.))",FitMin,FitMax);
        /*
        fTemplate_MC_NonCommon->SetParameter(0,1);
        fTemplate_MC_NonCommon->SetParLimits(0,0.5,1.5);
        fTemplate_MC_NonCommon->SetParameter(1,0.1);
        fTemplate_MC_NonCommon->SetParLimits(1,0.033,0.3);
        fTemplate_MC_NonCommon->SetParameter(2,1);
        fTemplate_MC_NonCommon->SetParLimits(2,0.33,3);
        fTemplate_MC_NonCommon->SetParameter(3,-1.4);
        fTemplate_MC_NonCommon->SetParLimits(3,-3,3);
        */
        //fTemplate_MC_NonCommon->SetParameter(0,1);
        //fTemplate_MC_NonCommon->SetParameter(1,0.1);
        //fTemplate_MC_NonCommon->FixParameter(2,1);
        //fTemplate_MC_NonCommon->SetParameter(3,-1.4);


        TF1* fTemplate_MC_NonCommon = new TF1("fTemplate_MC_NonCommon","[0]+[1]*exp(-pow((fabs(x)-[2])/[3],2.))+[4]*exp(-pow((fabs(x)-[5])/[6],2.))",FitMin,FitMax);
        fTemplate_MC_NonCommon->SetParameter(0,1);
        fTemplate_MC_NonCommon->SetParameter(1,-0.1);
        fTemplate_MC_NonCommon->FixParameter(2,0);
        fTemplate_MC_NonCommon->SetParameter(3,1);
        fTemplate_MC_NonCommon->SetParameter(4,0.1);
        fTemplate_MC_NonCommon->FixParameter(5,TMath::Pi());
        fTemplate_MC_NonCommon->SetParameter(6,1);

        InputHisto->Fit(fTemplate_MC_NonCommon,"Q, S, N, R, M");
        delete fInputFile;


        fInputFile = new TFile(fName_Data,"read");
        List1 = (TList*)fInputFile->Get(l1Name_Data);
        List2 = (TList*)List1->FindObject(l2Name_Data);
        InputHisto = (TH1D*)List2->FindObject(hName_Data);
        FitMin = -1.8;
        FitMax = 4.3;
        OutputFile->cd();
        if(uSyst==0) InputHisto->Write();
        else{
            float Value;
            float Error;
            float RanVal;
            for(unsigned uBin=0; uBin<InputHisto->GetNbinsX(); uBin++){
                Value = InputHisto->GetBinContent(uBin+1);
                Error = InputHisto->GetBinError(uBin+1);
                RanVal = rangen.Gaus(Value,Error);
                InputHisto->SetBinContent(uBin+1,RanVal);
                InputHisto->SetBinError(uBin+1,Error);
            }
        }
        //TF1* fTemplate_Data = new TF1("fTemplate_Data","[0]*(1.+[1]*sin([2]*x+[3])+[4]*exp(-pow((x-[5])/[6],2.)))+[7]*(1.+[8]*sin([9]*x+[10]))",
        //                              FitMin,FitMax);
        //TF1* fTemplate_Data = new TF1("fTemplate_Data",
    //"[18]*[0]*(1.+[1]*exp(-pow((x-[2])/[3],2.)))*(1.+[4]*x+[5]*x*x+[6]*pow(x,3.)+[7]*pow(x,4.)+[8]*pow(x,5.)+[9]*pow(x,6.)+[10]*pow(x,7.))+[19]*[11]*(1.+[12]*sin([13]*x+[14]))*(1.+[15]*x+[16]*x*x+[17]*pow(x,3.))",
    //                                  FitMin,FitMax);
        ////fTemplate_Data->SetParameter(0,0.1*fTemplate_MC_Common->GetParameter(0));
        //for(unsigned uPar=0; uPar<=10; uPar++) fTemplate_Data->FixParameter(uPar,fTemplate_MC_Common->GetParameter(uPar));

        ////fTemplate_Data->SetParameter(11,0.9*fTemplate_MC_NonCommon->GetParameter(0));
        //for(unsigned uPar=0; uPar<=6; uPar++) fTemplate_Data->FixParameter(uPar+11,fTemplate_MC_NonCommon->GetParameter(uPar));
        //fTemplate_Data->SetParameter(18,0.2);
        //fTemplate_Data->SetParameter(19,0.8);

        TF1* fTemplate_Data = new TF1("fTemplate_Data","[15]*([14]*([0]+[1]*exp(-pow((fabs(x)-[2])/[3],2.))+[4]*exp(-pow((fabs(x)-[5])/[6],2.)))+"
                                      "(1.-[14])*([7]+[8]*exp(-pow((fabs(x)-[9])/[10],2.))+[11]*exp(-pow((fabs(x)-[12])/[13],2.))))",FitMin,FitMax);
        for(unsigned uPar=0; uPar<=6; uPar++) fTemplate_Data->FixParameter(uPar,fTemplate_MC_Common->GetParameter(uPar));
        for(unsigned uPar=7; uPar<=13; uPar++) fTemplate_Data->FixParameter(uPar,fTemplate_MC_NonCommon->GetParameter(uPar-7));
        fTemplate_Data->SetParameter(14,0.1);
        fTemplate_Data->SetParameter(15,1);

        InputHisto->Fit(fTemplate_Data,"Q, S, N, R, M");

        delete fInputFile;

        OutputFile->cd();
        if(uSyst==0){
            fTemplate_MC_Common->Write();
            fTemplate_MC_NonCommon->Write();
            fTemplate_Data->Write();
            gWC_Def.SetPoint(0,0.9,fTemplate_Data->GetParameter(14));
            gWC_Def.SetPointError(0,0,fTemplate_Data->GetParError(14));
            StatErr = fTemplate_Data->GetParError(14);
        }

        hWeight_Common->Fill(fTemplate_Data->GetParameter(14));

        delete fTemplate_MC_Common;
        delete fTemplate_MC_NonCommon;
        delete fTemplate_Data;
    }
    cout << endl;
    gWC_Syst.SetPoint(0,1.1,hWeight_Common->GetMean());
    gWC_Syst.SetPointError(0,0,hWeight_Common->GetStdDev());

    gWC_TotErr.SetPoint(0,1.2,sqrt(StatErr*StatErr+hWeight_Common->GetMean()*hWeight_Common->GetMean()));
    gWC_TotErr.SetPointError(0,0,sqrt(StatErr*StatErr+hWeight_Common->GetStdDev()*hWeight_Common->GetStdDev()));

    OutputFile->cd();
    hWeight_Common->Write();
    gWC_Def.Write();
    gWC_Syst.Write();
    gWC_TotErr.Write();

    delete hWeight_Common;
    delete OutputFile;

}

void Fit_pL_CommonAncestorTemplate_dPhi(const bool useMC){

    const unsigned NumSystIter = 256;
    TH1F* hWeight_Common = new TH1F("hWeight_Common","hWeight_Common",4096,0,0.5);
    TGraphErrors gWC_Def;//with stat uncertainty
    gWC_Def.SetName("gWC_Def");
    gWC_Def.SetMarkerStyle(20);
    gWC_Def.SetMarkerSize(1);
    TGraphErrors gWC_Syst;
    gWC_Syst.SetName("gWC_Syst");
    gWC_Syst.SetMarkerStyle(20);
    gWC_Syst.SetMarkerSize(1);
    TGraphErrors gWC_TotErr;
    gWC_TotErr.SetName("gWC_TotErr");
    gWC_TotErr.SetMarkerStyle(20);
    gWC_TotErr.SetMarkerSize(1);

    TString fName_MC_Common = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/pLambda_Rush/Ancestors/dEtadPhiAncestors.root";
    TString fName_MC_NonCommon = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/pLambda_Rush/Ancestors/dEtadPhiAncestors.root";
    TString fName_Data;
    if(useMC)fName_Data = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/pLambda_Rush/dEtadPhi_Data/MC/dEtadPhi.root";
    else fName_Data = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/pLambda_Rush/dEtadPhi_Data/Data/dEtadPhi.root";
    TString l1Name_MC_Common = "ProtonLambda";
    TString l1Name_MC_NonCommon = "ProtonLambda";
    TString l1Name_Data = "ProtonLambda";
    TString l2Name_MC_Common = "pLCommon";
    TString l2Name_MC_NonCommon = "pLNonCommon";
    TString l2Name_Data = "pLIntegrated";
    TString hName_MC_Common = "SEdPhidEtaDistCommon_Particle0_Particle2ShiftedProjection";
    TString hName_MC_NonCommon = "SEdPhidEtaDistNonCommon_Particle0_Particle2AShiftedProjection";
    TString hName_Data = "SEdPhidEtaDist_Particle0_Particle2ShiftedProjection";

    TFile* OutputFile = new TFile(
    TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/OtherTasks/Fit_pL_CommonAncestorTemplate_dPhi/OutputFile_pL_%s.root",
                    useMC?"MC":"Data"),"recreate");
    double FitMin;
    double FitMax;
    double StatErr;

    TRandom3 rangen(11);

    for(unsigned uSyst=0; uSyst<NumSystIter; uSyst++){
        printf("\r\033[K uSyst=%u",uSyst);
        cout << flush;
        TFile* fInputFile;
        TList* List1;
        TList* List2;
        TH1D* InputHisto;

        fInputFile = new TFile(fName_MC_Common,"read");
        List1 = (TList*)fInputFile->Get(l1Name_MC_Common);
        List2 = (TList*)List1->FindObject(l2Name_MC_Common);
        InputHisto = (TH1D*)List2->FindObject(hName_MC_Common);
        FitMin = -1.8;
        FitMax = 4.3;
        OutputFile->cd();
        if(uSyst==0) InputHisto->Write();
        else{
            float Value;
            float Error;
            float RanVal;
            for(unsigned uBin=0; uBin<InputHisto->GetNbinsX(); uBin++){
                Value = InputHisto->GetBinContent(uBin+1);
                Error = InputHisto->GetBinError(uBin+1);
                RanVal = rangen.Gaus(Value,Error);
                InputHisto->SetBinContent(uBin+1,RanVal);
                InputHisto->SetBinError(uBin+1,Error);
            }
        }
        TF1* fTemplate_MC_Common = new TF1("fTemplate_MC_Common",
        "[0]+[1]*exp(-pow((fabs(x)-[2])/[3],2.))+[4]*exp(-pow((fabs(x)-[5])/[6],2.))",
        FitMin,FitMax);
        fTemplate_MC_Common->SetParameter(0,1);
        fTemplate_MC_Common->SetParameter(1,1);
        fTemplate_MC_Common->FixParameter(2,0);
        fTemplate_MC_Common->SetParameter(3,1);
        fTemplate_MC_Common->SetParameter(4,-0.1);
        fTemplate_MC_Common->FixParameter(5,TMath::Pi());
        fTemplate_MC_Common->SetParameter(6,1);

        InputHisto->Fit(fTemplate_MC_Common,"Q, S, N, R, M");
        delete fInputFile;

        fInputFile = new TFile(fName_MC_NonCommon,"read");
        List1 = (TList*)fInputFile->Get(l1Name_MC_NonCommon);
        List2 = (TList*)List1->FindObject(l2Name_MC_NonCommon);
        InputHisto = (TH1D*)List2->FindObject(hName_MC_NonCommon);
        //FitMin = InputHisto->GetBinLowEdge(1);
        //FitMax = InputHisto->GetBinLowEdge(InputHisto->GetNbinsX());
        FitMin = -1.8;
        FitMax = 4.3;
        OutputFile->cd();
        if(uSyst==0) {InputHisto->Write();}
        else{
            float Value;
            float Error;
            float RanVal;
            for(unsigned uBin=0; uBin<InputHisto->GetNbinsX(); uBin++){
                Value = InputHisto->GetBinContent(uBin+1);
                Error = InputHisto->GetBinError(uBin+1);
                RanVal = rangen.Gaus(Value,Error);
                InputHisto->SetBinContent(uBin+1,RanVal);
                InputHisto->SetBinError(uBin+1,Error);
            }
        }
        TF1* fTemplate_MC_NonCommon = new TF1("fTemplate_MC_NonCommon","[0]+[1]*exp(-pow((fabs(x)-[2])/[3],2.))+[4]*exp(-pow((fabs(x)-[5])/[6],2.))",FitMin,FitMax);
        fTemplate_MC_NonCommon->SetParameter(0,1);
        fTemplate_MC_NonCommon->SetParameter(1,-0.1);
        fTemplate_MC_NonCommon->FixParameter(2,0);
        fTemplate_MC_NonCommon->SetParameter(3,1);
        fTemplate_MC_NonCommon->SetParameter(4,0.1);
        fTemplate_MC_NonCommon->FixParameter(5,TMath::Pi());
        fTemplate_MC_NonCommon->SetParameter(6,1);

        InputHisto->Fit(fTemplate_MC_NonCommon,"Q, S, N, R, M");
        delete fInputFile;


        fInputFile = new TFile(fName_Data,"read");
        List1 = (TList*)fInputFile->Get(l1Name_Data);
        List2 = (TList*)List1->FindObject(l2Name_Data);
        InputHisto = (TH1D*)List2->FindObject(hName_Data);
        FitMin = -1.8;
        FitMax = 4.3;
        OutputFile->cd();
        if(uSyst==0) {InputHisto->Write();}
        else{
            float Value;
            float Error;
            float RanVal;
            for(unsigned uBin=0; uBin<InputHisto->GetNbinsX(); uBin++){
                Value = InputHisto->GetBinContent(uBin+1);
                Error = InputHisto->GetBinError(uBin+1);
                RanVal = rangen.Gaus(Value,Error);
                InputHisto->SetBinContent(uBin+1,RanVal);
                InputHisto->SetBinError(uBin+1,Error);
            }
        }

        TF1* fTemplate_Data = new TF1("fTemplate_Data","[15]*([14]*([0]+[1]*exp(-pow((fabs(x)-[2])/[3],2.))+[4]*exp(-pow((fabs(x)-[5])/[6],2.)))+"
                                      "(1.-[14])*([7]+[8]*exp(-pow((fabs(x)-[9])/[10],2.))+[11]*exp(-pow((fabs(x)-[12])/[13],2.))))",FitMin,FitMax);
        for(unsigned uPar=0; uPar<=6; uPar++) fTemplate_Data->FixParameter(uPar,fTemplate_MC_Common->GetParameter(uPar));
        for(unsigned uPar=7; uPar<=13; uPar++) fTemplate_Data->FixParameter(uPar,fTemplate_MC_NonCommon->GetParameter(uPar-7));
        fTemplate_Data->SetParameter(14,0.1);
        fTemplate_Data->SetParameter(15,1);

        fTemplate_Data->SetParameter(18,0.2);
        fTemplate_Data->SetParameter(19,0.8);

        InputHisto->Fit(fTemplate_Data,"Q, S, N, R, M");
        delete fInputFile;

        OutputFile->cd();
        if(uSyst==0){
            fTemplate_MC_Common->Write();
            fTemplate_MC_NonCommon->Write();
            fTemplate_Data->Write();
            gWC_Def.SetPoint(0,0.9,fTemplate_Data->GetParameter(14));
            gWC_Def.SetPointError(0,0,fTemplate_Data->GetParError(14));
            StatErr = fTemplate_Data->GetParError(14);
        }

        hWeight_Common->Fill(fTemplate_Data->GetParameter(14));

        delete fTemplate_MC_Common;
        delete fTemplate_MC_NonCommon;
        delete fTemplate_Data;
    }
    cout << endl;
    gWC_Syst.SetPoint(0,1.1,hWeight_Common->GetMean());
    gWC_Syst.SetPointError(0,0,hWeight_Common->GetStdDev());

    gWC_TotErr.SetPoint(0,1.2,sqrt(StatErr*StatErr+hWeight_Common->GetMean()*hWeight_Common->GetMean()));
    gWC_TotErr.SetPointError(0,0,sqrt(StatErr*StatErr+hWeight_Common->GetStdDev()*hWeight_Common->GetStdDev()));

    OutputFile->cd();
    hWeight_Common->Write();
    gWC_Def.Write();
    gWC_Syst.Write();
    gWC_TotErr.Write();

    delete hWeight_Common;
    delete OutputFile;

}


DLM_CommonAnaFunctions* DariuszSource_ppReso=NULL;
double DariuszSourceFunction_pp(double* Pars){
    if(!DariuszSource_ppReso){
        DariuszSource_ppReso = new DLM_CommonAnaFunctions();
        CATS Dummy;
        Dummy.SetMomBins(1,0,10);
        DariuszSource_ppReso->SetUpCats_pp(Dummy,"AV18","McGauss_ResoTM",0,202);
    }
    double CoreResoValue = DariuszSource_ppReso->GetCleverMcLevyResoTM_pp()->Eval(Pars);
    double GaussValue = GaussSource(Pars);
    if(Pars[1]<Pars[3] && CoreResoValue>GaussValue) return GaussValue;
    else return CoreResoValue;
}

void pp_DariuszSetUp(CATS& Kitty){
    Kitty.SetExcludeFailedBins(false);
    Kitty.SetQ1Q2(1);
    Kitty.SetPdgId(2212, 2212);
    Kitty.SetRedMass( 0.5*Mass_p );
    Kitty.SetNumChannels(4);
    Kitty.SetNumPW(0,3);
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

    //#,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
    double PotPars1S0[8]={NN_AV18,v18_Coupled3P2,1,1,1,0,0,0};
    double PotPars3P0[8]={NN_AV18,v18_Coupled3P2,1,1,1,1,1,0};
    double PotPars3P1[8]={NN_AV18,v18_Coupled3P2,1,1,1,1,1,1};
    double PotPars3P2[8]={NN_AV18,v18_Coupled3P2,1,1,1,1,1,2};
    double PotPars1D2[8]={NN_AV18,v18_Coupled3P2,1,1,1,0,2,2};
    CATSparameters cPotPars1S0(CATSparameters::tPotential,8,true); cPotPars1S0.SetParameters(PotPars1S0);
    CATSparameters cPotPars3P0(CATSparameters::tPotential,8,true); cPotPars3P0.SetParameters(PotPars3P0);
    CATSparameters cPotPars3P1(CATSparameters::tPotential,8,true); cPotPars3P1.SetParameters(PotPars3P1);
    CATSparameters cPotPars3P2(CATSparameters::tPotential,8,true); cPotPars3P2.SetParameters(PotPars3P2);
    CATSparameters cPotPars1D2(CATSparameters::tPotential,8,true); cPotPars1D2.SetParameters(PotPars1D2);

    Kitty.SetShortRangePotential(0,0,fDlmPot,cPotPars1S0);
    Kitty.SetShortRangePotential(0,2,fDlmPot,cPotPars1D2);
    Kitty.SetShortRangePotential(1,1,fDlmPot,cPotPars3P0);
    Kitty.SetShortRangePotential(2,1,fDlmPot,cPotPars3P1);
    Kitty.SetShortRangePotential(3,1,fDlmPot,cPotPars3P2);
}

//showing off, that the resonances source that converges towards a Gaussian (i.e. remove the effect of particle coming together) does not matter.
//so we build up a source function, that is GaussCore+Reso in general, but below 2 fm if S_core > S_gauss, we switch back to a Gaussian source
void pp_DariuszSource(){
    TString OutputFileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/OtherTasks/pp_DariuszSource/pp_DariuszSource.root";
    CATS KittyD;
    KittyD.SetMomBins(100,0,400);
    KittyD.SetThetaDependentSource(false);
    CATSparameters cPars(CATSparameters::tSource,1,true);
    cPars.SetParameter(0,1.0);
    KittyD.SetAnaSource(DariuszSourceFunction_pp, cPars);
    KittyD.SetAutoNormSource(true);
    KittyD.SetAnaSource(0,cPars.GetParameter(0));
    KittyD.SetUseAnalyticSource(true);
    KittyD.SetMomentumDependentSource(false);
    pp_DariuszSetUp(KittyD);
    KittyD.KillTheCat();

    CATS KittyR;
    KittyR.SetMomBins(100,0,400);
    DLM_CommonAnaFunctions ANA;
    ANA.SetUpCats_pp(KittyR,"AV18","McGauss_ResoTM",0,202);
    KittyR.KillTheCat();

    RootFile_DlmCk(OutputFileName,"Ck_pp_Dariusz",&KittyD);
    RootFile_DlmCk(OutputFileName,"Ck_pp_Reso",&KittyR);

    RootFile_DlmSource(OutputFileName,"Spdf_pp_Dariusz",&KittyD,128,0,10,1,true);
    RootFile_DlmSource(OutputFileName,"Spdf_pp_Reso",&KittyR,128,0,10,1,true);

    RootFile_DlmSource(OutputFileName,"Sr_pp_Dariusz",&KittyD,128,0,10,1,false);
    RootFile_DlmSource(OutputFileName,"Sr_pp_Reso",&KittyR,128,0,10,1,false);


}


DLM_Histo<float> huWF;
DLM_Histo<float> hwWF;
float Evaluate_d_u(double Radius){
    return huWF.Eval(&Radius);
}
float Evaluate_d_w(double Radius){
    return hwWF.Eval(&Radius);
}
/*
//par[0] is the order of the pol
//par[1/2] are the two limits at which the derivative should be flat
//if on limit is desired, set par[] = -1e6
double PolN_flat_end(double* x, double* par){
    double lim1 = par[1];
    double lim2 = par[2];
    double temp;
    unsigned Order = par[0];
    double* pol = &par[3];
    if(lim2!=-1e6){
        pol[2] = 0;
        for(unsigned uOrd=3; uOrd<=Order; uOrd++){
            pol[2] -= uOrd*pol[uOrd]*(pow(lim1,uOrd-1)+pow(lim2,uOrd-1));
        }
        pol[2] /= (2.*(lim1+lim2));
    }
    if(lim1!=-1e6){
        pol[1] = 0;
        for(unsigned uOrd=2; uOrd<=Order; uOrd++){
            pol[1] += uOrd*pol[uOrd]*(pow(lim1,uOrd-1)+pow(lim2,uOrd-1));
        }
    }
    double Result=0;
    for(unsigned uOrd=0; uOrd<=Order; uOrd++){
        Result += pol[uOrd]*pow(*x,uOrd);
    }

    return Result;
}
*/

//par[0] is the order of the pol
//par[1] is the limit at which the derivative should be flat
//par[2...] are the pol coefficients (0,1...)
double PolN_flat_end(double* x, double* par){
    double& lim = par[1];
    unsigned Order = par[0];
    double* pol = &par[2];
    if(lim!=-1e6){
        pol[1] = 0;
        for(unsigned uOrd=2; uOrd<=Order; uOrd++){
            pol[1] -= uOrd*pol[uOrd]*(pow(lim,uOrd-1));
        }
    }
    double Result=0;
    for(unsigned uOrd=0; uOrd<=Order; uOrd++){
        Result += pol[uOrd]*pow(*x,uOrd);
    }
    return Result;
}


void Test_flat_pol(){

    //pol5
    //fix at 0
    TF1* fun1 = new TF1("fun1",PolN_flat_end,0,500,2+4);
    //fix at 300
    TF1* fun2 = new TF1("fun2",PolN_flat_end,0,500,2+4);

    fun1->FixParameter(0,3);
    fun1->FixParameter(1,0);
    fun1->SetParameter(2,1);//0
    fun1->SetParameter(3,1);//1
    fun1->SetParameter(4,0.01);//2
    fun1->SetParameter(5,0.0001);//3

    fun2->FixParameter(0,3);
    fun2->FixParameter(1,300);
    fun2->SetParameter(2,1);//0
    fun2->SetParameter(3,1);//1
    fun2->SetParameter(4,0.01);//2
    fun2->SetParameter(5,0.0001);//3

    TFile fOutput("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/OtherTasks/Test_flat_pol/fOutput.root","recreate");
    fun1->Write();
    fun2->Write();

}


void TestOmegaOmega(){

    const unsigned NumBins = 200;
    CATS Kitty;
    CATSparameters cPars(CATSparameters::tSource,1,true);
    cPars.SetParameter(0,0.85);

    Kitty.SetMomBins(NumBins,0,400);

    Kitty.SetAnaSource(GaussSource, cPars);
    Kitty.SetUseAnalyticSource(true);

    Kitty.SetMomentumDependentSource(false);
    //Kitty.SetThetaDependentSource(false);
    Kitty.SetExcludeFailedBins(false);

    Kitty.SetQ1Q2(1);
    Kitty.SetPdgId(2212, 2212);
    Kitty.SetRedMass( 0.5*MassOmega );

    Kitty.SetNumChannels(4);
    Kitty.SetNumPW(0,0);
    Kitty.SetNumPW(1,0);
    Kitty.SetNumPW(2,0);
    Kitty.SetNumPW(3,0);
    Kitty.SetSpin(0,0);
    Kitty.SetSpin(1,1);
    Kitty.SetSpin(2,2);
    Kitty.SetSpin(3,3);
    Kitty.SetChannelWeight(0, 1./16.);
    Kitty.SetChannelWeight(1, 3./16.);
    Kitty.SetChannelWeight(2, 5./16.);
    Kitty.SetChannelWeight(3, 7./16.);

    Kitty.KillTheCat();

    TFile* RootFile = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/OtherTasks/TestOmegaOmega/OO_Dimi.root","recreate");
    TGraph graph;
    graph.SetName("Ck_OmegaOmega");
    graph.Set(NumBins);
    for(unsigned uBin=0; uBin<NumBins; uBin++){
        graph.SetPoint(uBin,Kitty.GetMomentum(uBin),Kitty.GetCorrFun(uBin));
    }
    graph.Write("",TObject::kOverwrite);
    delete RootFile;



}

void pXi_Bug_vs_True_Potential(const unsigned WhichChannel, const TString ChannelName){
    const unsigned NumMomBins = 100;
    const double kMin = 0;
    const double kMax = 200;
    const unsigned NumRadBins = 512;
    const double rMin = 0;
    const double rMax = 32;
    DLM_Histo<float> DummyRad;
    DummyRad.SetUp(1);
    DummyRad.SetUp(0,NumRadBins,rMin,rMax);
    DummyRad.Initialize();

    const unsigned NumKstarBins = 4;
    double* WF_Kstar = new double [NumKstarBins];
    WF_Kstar[0] = 10;
    WF_Kstar[1] = 30;
    WF_Kstar[2] = 50;
    WF_Kstar[3] = 70;

    DLM_CommonAnaFunctions AnaObj;

    CATS KittyOld;
    KittyOld.SetMomBins(NumMomBins,kMin,kMax);
    AnaObj.SetUpCats_pXim(KittyOld,"pXim_HALQCD1","Gauss",-12,0);
    KittyOld.SetAnaSource(0,1.0);
    KittyOld.KillTheCat();

    CATS KittyNew;
    KittyNew.SetMomBins(NumMomBins,kMin,kMax);
    AnaObj.SetUpCats_pXim(KittyNew,"pXim_HALQCD1","Gauss",12,0);
    KittyNew.SetAnaSource(0,1.0);
    KittyNew.KillTheCat();

    TString OutputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/OtherTasks/pXi_Bug_vs_True_Potential/";
    OutputFolder += ChannelName;
    OutputFolder += "_";
    TFile fOutput(OutputFolder+"fOutput.root","recreate");
    TGraph* WF_I0S0_Old = new TGraph [NumKstarBins];
    TGraph* WF_I0S0_New = new TGraph [NumKstarBins];
    TCanvas** cWF = new TCanvas* [NumKstarBins];
    TLegend** wfLegend = new TLegend* [NumKstarBins];
    TPaveText** PT_ks = new TPaveText* [NumKstarBins];

    TH1F* hWF_Dummy = new TH1F("hWF_Dummy","hWF_Dummy",128,rMin,rMax);
    hWF_Dummy->SetStats(false);
    hWF_Dummy->SetTitle("");

    hWF_Dummy->GetXaxis()->SetTitleSize(0.07);
    hWF_Dummy->GetXaxis()->SetLabelSize(0.07);
    hWF_Dummy->GetXaxis()->SetTitleOffset(1.3);
    hWF_Dummy->GetXaxis()->SetLabelOffset(0.02);

    hWF_Dummy->GetYaxis()->SetTitleSize(0.07);
    hWF_Dummy->GetYaxis()->SetLabelSize(0.07);
    hWF_Dummy->GetYaxis()->SetTitleOffset(1.00);
    hWF_Dummy->GetYaxis()->SetNdivisions(505);

    hWF_Dummy->GetXaxis()->SetTitle("r (fm)");
    hWF_Dummy->GetYaxis()->SetTitle("u_{0}(kr)/r");
    hWF_Dummy->GetYaxis()->SetRangeUser(-0.5,5.0);

    for(unsigned uKstar=0; uKstar<NumKstarBins; uKstar++){
        WF_I0S0_Old[uKstar].SetName(TString::Format("WF_I0S0_Old_%.0f",WF_Kstar[uKstar]));
        WF_I0S0_New[uKstar].SetName(TString::Format("WF_I0S0_New_%.0f",WF_Kstar[uKstar]));
        WF_I0S0_Old[uKstar].SetLineColor(kOrange+3);
        WF_I0S0_New[uKstar].SetLineColor(kOrange+1);
        WF_I0S0_Old[uKstar].SetLineWidth(4);
        WF_I0S0_New[uKstar].SetLineWidth(4);
        WF_I0S0_Old[uKstar].SetLineStyle(4);
        WF_I0S0_New[uKstar].SetLineStyle(1);

        for(unsigned uRad=0; uRad<NumRadBins; uRad++){
            double RAD = DummyRad.GetBinCenter(0,uRad);
            unsigned WhichMomBin = KittyOld.GetMomBin(WF_Kstar[uKstar]);
            WF_I0S0_Old[uKstar].SetPoint(uRad,RAD,real(KittyOld.EvalRadialWaveFunction(WhichMomBin,WhichChannel,0,RAD,true)));
            WF_I0S0_New[uKstar].SetPoint(uRad,RAD,real(KittyNew.EvalRadialWaveFunction(WhichMomBin,WhichChannel,0,RAD,true)));
        }

        wfLegend[uKstar] = new TLegend(0.6,0.70,0.95,0.95);//lbrt
        wfLegend[uKstar]->SetName(TString::Format("wfLegend_I0S0_%.0f",WF_Kstar[uKstar]));
        wfLegend[uKstar]->SetTextSize(0.06);
        wfLegend[uKstar]->AddEntry(&WF_I0S0_Old[uKstar],"With Bug");
        wfLegend[uKstar]->AddEntry(&WF_I0S0_New[uKstar],"Without Bug");

        PT_ks[uKstar] = new TPaveText(0.25,0.75,0.6,0.95, "blNDC");//lbrt
        PT_ks[uKstar]->SetName(TString::Format("PT_ks_%.0f",WF_Kstar[uKstar]));
        PT_ks[uKstar]->SetBorderSize(1);
        PT_ks[uKstar]->SetTextSize(0.055);
        PT_ks[uKstar]->SetFillColor(kWhite);
        PT_ks[uKstar]->SetTextFont(22);
        PT_ks[uKstar]->AddText(TString::Format("Ch: %s; k*=%.0f MeV",ChannelName.Data(),WF_Kstar[uKstar]));

        cWF[uKstar] = new TCanvas(TString::Format("cWF_I0S0_%.0f",WF_Kstar[uKstar]), TString::Format("cWF_I0S0_%.0f",WF_Kstar[uKstar]), 1);
        cWF[uKstar]->cd(0); cWF[uKstar]->SetCanvasSize(1920, 1080); cWF[uKstar]->SetMargin(0.15,0.05,0.2,0.05);//lrbt
        cWF[uKstar]->SetGrid(true);

        hWF_Dummy->Draw("axis");
        WF_I0S0_Old[uKstar].Draw("C,same");
        WF_I0S0_New[uKstar].Draw("C,same");
        wfLegend[uKstar]->Draw("same");
        PT_ks[uKstar]->Draw("same");

        cWF[uKstar]->SaveAs(OutputFolder+TString::Format("cWF_%.0f.pdf",WF_Kstar[uKstar]));

        WF_I0S0_Old[uKstar].Write();
        WF_I0S0_New[uKstar].Write();
        wfLegend[uKstar]->Write();
        PT_ks[uKstar]->Write();
        cWF[uKstar]->Write();
    }

    delete [] WF_I0S0_Old;
    delete [] WF_I0S0_New;
    for(unsigned uKstar=0; uKstar<NumKstarBins; uKstar++) {
        delete wfLegend[uKstar];
        delete PT_ks[uKstar];
        delete cWF[uKstar];
    }
    delete [] wfLegend;
    delete [] PT_ks;
    delete [] cWF;
    delete [] WF_Kstar;
}

void MemBugHunt(){
    const unsigned NumMomBins = 16;
    const double kMin = 0;
    const double kMax = 320;

    unsigned uDummy;
    double* MomBins_pL = NULL;
    double* FitRegion_pL = NULL;

    DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder("/home/dmihaylov/CernBox/CatsFiles");
    AnalysisObject.SetUpBinning_pL("pp13TeV_HM_Dec19",uDummy,MomBins_pL,FitRegion_pL,0,0);

    CATS AB_pL;
    DLM_Ck* Ck_pL;
    AB_pL.SetMomBins(NumMomBins,kMin,kMax);
    AnalysisObject.SetUpCats_pL(AB_pL,"Chiral_Coupled_SPD","Gauss",11600,202);//NLO_Coupled_S

    double CuspWeight = 0.33;//0.54
    AB_pL.SetChannelWeight(7,1./4.*CuspWeight);//1S0 SN(s) -> LN(s)
    AB_pL.SetChannelWeight(8,3./4.*CuspWeight);//3S1 SN(s) -> LN(s)
    AB_pL.SetChannelWeight(10,3./4.*CuspWeight);//3S1 SN(d) -> LN(s)
    AB_pL.SetChannelWeight(13,3./20.*CuspWeight);//3D1 SN(d) -> LN(d)
    AB_pL.SetChannelWeight(15,3./20.*CuspWeight);//3D1 SN(s) -> LN(d)

    AB_pL.SetAnaSource(0,1.4);
    AB_pL.SetNotifications(CATS::nError);
    AB_pL.KillTheCat();

    CATS AB_pXim;
    //same binning as pL, as we only use pXim as feed-down
    AB_pXim.SetMomBins(NumMomBins,0,kMax);
    AnalysisObject.SetUpCats_pXim(AB_pXim,"pXim_HALQCD1","Gauss");
    AB_pXim.SetAnaSource(0,1.5);
    AB_pXim.SetNotifications(CATS::nError);
    AB_pXim.KillTheCat();

    Ck_pL = new DLM_Ck(AB_pL.GetNumSourcePars(),0,AB_pL,NumMomBins*2,kMin,kMax*2);
    Ck_pL->SetSourcePar(0,AB_pL.GetAnaSourcePar(0));
    Ck_pL->SetCutOff(360,700);

    DLM_Ck* Ck_pSigma0 = new DLM_Ck(1,0,NumMomBins*2,kMin,kMax*2,Lednicky_gauss_Sigma0);
    Ck_pSigma0->SetSourcePar(0,1.5);
    Ck_pSigma0->SetCutOff(360,700);

    DLM_Ck* Ck_pXim = new DLM_Ck(AB_pXim.GetNumSourcePars(),0,AB_pXim,NumMomBins*2,kMin,kMax*2);
    Ck_pXim->SetCutOff(360,700);
    Ck_pL->Update();
    Ck_pSigma0->Update();
    Ck_pXim->Update();

    TH2F* hResolution_pL = AnalysisObject.GetResolutionMatrix("pp13TeV_HM_Dec19","pLambda");
    TH2F* hResidual_pL_pSigma0 = AnalysisObject.GetResidualMatrix("pLambda","pSigma0");
    TH2F* hResidual_pL_pXim = AnalysisObject.GetResidualMatrix("pLambda","pXim");

    DLM_Histo<float>* dlmResolution_pL = Convert_TH2F_DlmHisto(hResolution_pL);
    DLM_Histo<float>* dlmResidual_pL_pSigma0 = Convert_TH2F_DlmHisto(hResidual_pL_pSigma0);
    DLM_Histo<float>* dlmResidual_pL_pXim = Convert_TH2F_DlmHisto(hResidual_pL_pXim);

    const unsigned NumIter = 10000;

    TRandom3 rangen(11);

    unsigned Progress;
    for(unsigned uIter=0; uIter<NumIter; uIter++){
        printf("\r\033[K Progress=%u",uIter);
        cout << flush;
        AnalysisObject.SetUpBinning_pL("pp13TeV_HM_Dec19",uDummy,MomBins_pL,FitRegion_pL,0,0);

        double lam_pL[5];
        double lam_pXim[5];

        AnalysisObject.SetUpLambdaPars_pL("pp13TeV_HM_Dec19",0,0,lam_pL);
        AnalysisObject.SetUpLambdaPars_pXim("pp13TeV_HM_Dec19",0,0,lam_pXim);

/*
        DLM_CkDecomposition* CkDec_pL = new DLM_CkDecomposition("pLambda",4,*Ck_pL,hResolution_pL);
        DLM_CkDecomposition* CkDec_pSigma0 = new DLM_CkDecomposition("pSigma0",0,*Ck_pSigma0,NULL);
        DLM_CkDecomposition* CkDec_pXim = new DLM_CkDecomposition("pXim",2,*Ck_pXim,NULL);
    usleep(10e3);
    delete CkDec_pL;
    delete CkDec_pSigma0;
    delete CkDec_pXim;
    continue;
        CkDec_pL->AddContribution(0,lam_pL[1],DLM_CkDecomposition::cFeedDown,CkDec_pSigma0,hResidual_pL_pSigma0);
        CkDec_pL->AddContribution(1,lam_pL[2],DLM_CkDecomposition::cFeedDown,CkDec_pXim,hResidual_pL_pXim);
        //CkDec_pL->AddContribution(0,lam_pL[1],DLM_CkDecomposition::cFeedDown);
        //CkDec_pL->AddContribution(1,lam_pL[2],DLM_CkDecomposition::cFeedDown);
        CkDec_pL->AddContribution(2,lam_pL[3],DLM_CkDecomposition::cFeedDown);
        CkDec_pL->AddContribution(3,lam_pL[4],DLM_CkDecomposition::cFake);//0.03

        //for Xim we simplify a bit and take ALL feed-down as flat
        CkDec_pXim->AddContribution(0,lam_pXim[1]+lam_pXim[2]+lam_pXim[3],DLM_CkDecomposition::cFeedDown);
        CkDec_pXim->AddContribution(1,lam_pXim[4],DLM_CkDecomposition::cFake);
*/

        DLM_CkDecomp* CkDec_pL = new DLM_CkDecomp("pLambda",4,*Ck_pL,dlmResolution_pL);
        DLM_CkDecomp* CkDec_pSigma0 = new DLM_CkDecomp("pSigma0",0,*Ck_pSigma0,NULL);
        DLM_CkDecomp* CkDec_pXim = new DLM_CkDecomp("pXim",2,*Ck_pXim,NULL);
    usleep(10e3);
    delete CkDec_pL;
    delete CkDec_pSigma0;
    delete CkDec_pXim;
    continue;
        CkDec_pL->AddContribution(0,lam_pL[1],DLM_CkDecomp::cFeedDown,CkDec_pSigma0,dlmResidual_pL_pSigma0);
        CkDec_pL->AddContribution(1,lam_pL[2],DLM_CkDecomp::cFeedDown,CkDec_pXim,dlmResidual_pL_pXim);
        //CkDec_pL->AddContribution(0,lam_pL[1],DLM_CkDecomp::cFeedDown);
        //CkDec_pL->AddContribution(1,lam_pL[2],DLM_CkDecomp::cFeedDown);
        CkDec_pL->AddContribution(2,lam_pL[3],DLM_CkDecomp::cFeedDown);
        CkDec_pL->AddContribution(3,lam_pL[4],DLM_CkDecomp::cFake);//0.03

        //for Xim we simplify a bit and take ALL feed-down as flat
        CkDec_pXim->AddContribution(0,lam_pXim[1]+lam_pXim[2]+lam_pXim[3],DLM_CkDecomp::cFeedDown);
        CkDec_pXim->AddContribution(1,lam_pXim[4],DLM_CkDecomp::cFake);

        CkDec_pL->Update(true,true);

        const unsigned NumFitIter = 0;
        for(unsigned uFitIter=0; uFitIter<NumFitIter; uFitIter++){
            double par[5];
            par[0] = rangen.Uniform(1.3,1.5);
            par[2] = rangen.Uniform(0.2,0.4);
            par[3] = rangen.Uniform(350,370);
            par[4] = rangen.Uniform(600,800);

            AB_pL.SetChannelWeight(7,1./4.*par[2]);//1S0 SN(s) -> LN(s)
            AB_pL.SetChannelWeight(8,3./4.*par[2]);//3S1 SN(s) -> LN(s)
            AB_pL.SetChannelWeight(10,3./4.*par[2]);//3S1 SN(d) -> LN(s)
            AB_pL.SetChannelWeight(13,3./20.*par[2]);//3D1 SN(d) -> LN(d)
            AB_pL.SetChannelWeight(15,3./20.*par[2]);//3D1 SN(s) -> LN(d)
            AB_pL.SetAnaSource(0,par[0],true);

            AB_pL.KillTheCat();
            CkDec_pL->GetCk()->SetCutOff(par[3],par[4]);
            CkDec_pL->Update(true);
        }

        CkDec_pL->Update(true,true);


    }
}

void DlmHistoMemBugHunt(){
    const unsigned NumIter = 100000;
    //DLM_Histo<float>* dlmHisto;
    DLM_Histo<float>* dlmHisto = new DLM_Histo<float>();
    dlmHisto->SetUp(2);
    dlmHisto->SetUp(0,128,0,128);
    dlmHisto->SetUp(1,128,0,128);
    dlmHisto->Initialize();
    TRandom3 rangen(11);
    for(unsigned uIter=0; uIter<NumIter; uIter++){
        //printf("\r\033[K Progress=%u",uIter);
        //cout << flush;
        printf("Progress=%u\n",uIter);

        //dlmHisto = new DLM_Histo<float>();
        dlmHisto->SetUp(2);
        unsigned NumBins = 118+rangen.Integer(21);
        dlmHisto->SetUp(0,NumBins,0,float(NumBins));
        dlmHisto->SetUp(1,NumBins,0,float(NumBins));
        dlmHisto->Initialize();
        usleep(10e3);

        //delete dlmHisto;
    }

    delete dlmHisto;
}

//f(x) = (3+0.001x+(−(0.6∙10^−3.5)x^2+(0.1∙10^−6)∙x^3)∙exp(−(x/200)^2))
double NonFlatCorrelation1(const double &Momentum, const double* SourcePar, const double* PotPar){
    double BASE = (3.+0.001*Momentum+(-pow(Momentum,2.)*0.6*pow(10.,-3.5)+0.1*pow(Momentum,3.)*pow(10.,-6))*exp(-pow(Momentum/200,2.)));
    double PEAK = 2.*exp(-pow((Momentum-200.)/5.,2.));
    return BASE+PEAK;
}
void SmearNonFlatCorrelation(){
    const unsigned NumBinsCk = 250;
    const double kMin = 0;
    const double kMax = 500;
    DLM_Ck dlmCk(0,0,NumBinsCk,kMin,kMax,NonFlatCorrelation1);
    dlmCk.Update();

    DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder("/home/dmihaylov/CernBox/CatsFiles");
    TH2F* hResolution = AnalysisObject.GetResolutionMatrix("pp13TeV_HM_Dec19","pp");

    DLM_CkDecomposition* CkDecomp = new DLM_CkDecomposition("dlmCk",0,dlmCk,hResolution);
    CkDecomp->Update();

    TFile fOutput("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/OtherTasks/SmearNonFlatCorrelation/fOutput1.root","recreate");
    TGraph CkOriginal;
    CkOriginal.SetName("CkOriginal");
    CkOriginal.SetLineColor(kGreen+2);
    CkOriginal.SetLineWidth(3);
    TGraph CkSmeared;
    CkSmeared.SetName("CkSmeared");
    CkSmeared.SetLineColor(kAzure+2);
    CkSmeared.SetLineWidth(3);
    for(unsigned uBin=0; uBin<NumBinsCk; uBin++){
        double MOM = dlmCk.GetBinCenter(0,uBin);
        CkOriginal.SetPoint(uBin,MOM,dlmCk.GetBinContent(uBin));
        CkSmeared.SetPoint(uBin,MOM,CkDecomp->EvalCk(MOM));
    }

    CkOriginal.Write();
    CkSmeared.Write();

    delete CkDecomp;
    delete hResolution;
}

//double Spline_Fitter(double* xVal, double* pars){
//    return DLM_FITTER2_FUNCTION_SPLINE3(xVal,pars);
//}
TF1* MakeSmoothAngularSourceDisto(TH1F* hAngleDisto){

  const double FitMin = 0;
  const double FitMax = TMath::Pi();
  const unsigned NumKnots = 6;
  const double KnotStep = (FitMax-FitMin)/double(NumKnots-1);
  double* Nodes_x = new double [NumKnots];
  //int iKnot = 0;
  Nodes_x[0] = 0;
  for(unsigned uKnot=1; uKnot<NumKnots-1; uKnot++){
    Nodes_x[uKnot] = Nodes_x[uKnot-1]+KnotStep;
  }
  Nodes_x[NumKnots-1] = TMath::Pi();
  for(unsigned uKnot=0; uKnot<NumKnots; uKnot++){
    printf("QA: Nodes_x[%u] = %.4f\n",uKnot,Nodes_x[uKnot]);
  }

  TF1* fit_angle = new TF1("fit_angle",DLM_FITTER2_FUNCTION_SPLINE3,FitMin,FitMax,3+NumKnots*2);
  //number of knots
  fit_angle->FixParameter(0,NumKnots);
  //derivative at the first knot
  fit_angle->SetParameter(1,0.1);
  fit_angle->SetParLimits(1,0,100);
  //derivative at the last knot
  fit_angle->SetParameter(2,-0.1);
  fit_angle->SetParLimits(2,-100,0);
  //knots [3, 3+NumKnots-1] aret the position of the knot on the x-axis (typically fixed)
  //knots [3+NumKnots, 3+2*NumKnots-1] aret the position of the knot on the y-axis (typically fitted)
  for(unsigned uKnot=0; uKnot<NumKnots; uKnot++){
      double HistVal = hAngleDisto->GetBinContent(hAngleDisto->FindBin(Nodes_x[uKnot]));
      fit_angle->FixParameter(3+uKnot,Nodes_x[uKnot]);
      fit_angle->SetParameter(3+NumKnots+uKnot,HistVal);
      fit_angle->SetParLimits(3+NumKnots+uKnot,0,HistVal*2);
  }
  fit_angle->SetNpx(1024);
  hAngleDisto->Fit(fit_angle,"S, N, R, M");

  delete [] Nodes_x;

  return fit_angle;
}
TF1* MakeSmoothAngularCosine(TH1F* hCosDisto){

  const double FitMin = -1;
  const double FitBl = 0.7;
  const double FitMax = 1;

  //hCosDisto->Scale(1./hCosDisto->Integral(),"width");
  const int NumBins = hCosDisto->GetNbinsX();
  TF1* fBaseline = new TF1("fBaseline","[0]",FitBl,FitMax);
  hCosDisto->Fit(fBaseline,"S, N, R, M");
  const double BlPar = fBaseline->GetParameter(0);
  const double BlErr = fBaseline->GetParError(0);
  TF1* fit_cos = new TF1("fit_cos","[0]/(1.-[1]*exp(-[1]*[2]*(x+1.)))",FitMin,FitMax);
  fit_cos->SetParameter(0,BlPar);
  fit_cos->SetParLimits(0,0,BlPar+2.*BlErr);
  fit_cos->SetParameter(1,0.5);
  fit_cos->SetParLimits(1,0.1,1.0);
  fit_cos->SetParameter(2,5);
  fit_cos->SetParLimits(2,0,25.0);
  fit_cos->SetNpx(1024);
  hCosDisto->Fit(fit_cos,"S, N, R, M");

  delete fBaseline;

  return fit_cos;
}

//0 = Original
//1 = fit the angle with splines
//2 = fit the cos with 1/(1-exp)
void Georgios_LXi_ResoTest(const int& SmoothSampling){

    //const double CoreSize = 0.842;
    //const double CoreSize = 0.894;
    const double CoreSize = 0.948;
    const unsigned SmoothEntires = 1024*128;

    //DLM_CleverMcLevyResoTM* MagicSource = new DLM_CleverMcLevyResoTM ();
    DLM_CleverMcLevyResoTM MagicSource;

    //DO NOT CHANGE !!! Sets up numerical bullshit, tuned for a Gaussian source
    MagicSource.InitStability(1,2-1e-6,2+1e-6);
    MagicSource.InitScale(38,0.15,2.0);
    MagicSource.InitRad(257*2,0,64);
    MagicSource.InitType(2);
    ///////////////////

    //for p-Xi, set up the amount of secondaries
    //first for the protons (64.22%)
    //MagicSource.SetUpReso(0,0.6422);
    MagicSource.SetUpReso(0,0.6438);
    //than for the Xis, here its 0% (we have ONLY primordials)
    MagicSource.SetUpReso(1,0.0);
    MagicSource.InitNumMcIter(1000000);

    //the cut off scale in k*, for which the angular distributions from EPOS
    //are evaluated. 200 MeV works okay, you can go up to 300 MeV for systematic checks
    const double k_CutOff = 200;

    //to be used for the NTuple later on
    Float_t k_D;
    Float_t fP1;
    Float_t fP2;
    Float_t fM1;
    Float_t fM2;
    Float_t Tau1;
    Float_t Tau2;
    Float_t AngleRcP1;
    Float_t AngleRcP2;
    Float_t AngleP1P2;
    //random generator dimi style. The input is incompatible with the ROOT random generator,
    //do not mix and match, do not ask me how I know this. Ask Bernie.
    //11 is the seed, you can change that to you favorite number
    DLM_Random RanGen(11);
    //dummies to save random shit
    double RanVal1;
    double RanVal2;
    double RanVal3;
    double RanCos;
    double MeanP1=0;
//printf("1\n");
    //open the magic file from dimi with the angular distributions.
    TFile* F_EposDisto_pReso_Xim = new TFile(TString::Format("%s/CatsFiles/Source/EposAngularDist/Epos_LamReso_Xim.root",GetCernBoxDimi()));
    //set up the ntuple, do not change anything unless told so by dimi
    TNtuple* T_EposDisto_pReso_Xim = (TNtuple*)F_EposDisto_pReso_Xim->Get("InfoTuple_ClosePairs");
    unsigned N_EposDisto_pReso_Xim = T_EposDisto_pReso_Xim->GetEntries();
    T_EposDisto_pReso_Xim->SetBranchAddress("k_D",&k_D);
    T_EposDisto_pReso_Xim->SetBranchAddress("P1",&fP1);
    T_EposDisto_pReso_Xim->SetBranchAddress("P2",&fP2);
    T_EposDisto_pReso_Xim->SetBranchAddress("M1",&fM1);
    T_EposDisto_pReso_Xim->SetBranchAddress("M2",&fM2);
    T_EposDisto_pReso_Xim->SetBranchAddress("Tau1",&Tau1);
    T_EposDisto_pReso_Xim->SetBranchAddress("Tau2",&Tau2);
    T_EposDisto_pReso_Xim->SetBranchAddress("AngleRcP1",&AngleRcP1);
    T_EposDisto_pReso_Xim->SetBranchAddress("AngleRcP2",&AngleRcP2);
    T_EposDisto_pReso_Xim->SetBranchAddress("AngleP1P2",&AngleP1P2);

    //TH1F* hAngleRcP1 = new TH1F("hAngleRcP1","hAngleRcP1",1024,0,TMath::Pi());

    gROOT->cd();
    TH1F* hAngle = new TH1F("hAngle","hAngle",32,0,TMath::Pi());
    TH1F* hCos = new TH1F("hCos","hCos",32,-1.,1.);
    TH1F* hFinalAngle = new TH1F("hFinalAngle","hFinalAngle",32,0,TMath::Pi());
    TH1F* hFinalCos = new TH1F("hFinalCos","hFinalCos",32,-1.,1.);
    F_EposDisto_pReso_Xim->cd();
    int NumUsefulEntries = 0;
    for(unsigned uEntry=0; uEntry<N_EposDisto_pReso_Xim; uEntry++){
      T_EposDisto_pReso_Xim->GetEntry(uEntry);
      if(k_D>k_CutOff) continue;
      hAngle->Fill(AngleRcP1);
      hCos->Fill(cos(AngleRcP1));
      MeanP1 += fP1;
      NumUsefulEntries++;
    }
    MeanP1 /= double(NumUsefulEntries);

    //hAngle->SetBinContent(hAngle->FindBin(0.275),0.31);
    hCos->Scale(1./hAngle->Integral(),"width");
//printf("2\n");

    if(!SmoothSampling){
      //iterate over the ntuple
      for(unsigned uEntry=0; uEntry<N_EposDisto_pReso_Xim; uEntry++){
          //get each entry
          T_EposDisto_pReso_Xim->GetEntry(uEntry);
          //disregard the entry of you are outside the desired k*
          if(k_D>k_CutOff) continue;
          //overwrite the value for the lifetime. This is computed from the
          //stat. hadronization model (Vale) or thermal fist (Max)
          //this is the value for the secondary protons
          //Tau1 = 1.65;
          Tau1 = 4.69;
          //for primoridials (the Xis) we put 0
          Tau2 = 0;
          //put in the average mass of the resonances (again from SHM or TF)
          //this is the value for protons
          //fM1 = 1361.52;
          fM1 = 1462.93;
          //generate a random path length for the propagation of the resonances
          //nothing to change!
          RanVal1 = RanGen.Exponential(fM1/(fP1*Tau1));
          //adds a single entry into the PDF for the angular distribution to be used
          MagicSource.AddBGT_RP(RanVal1,cos(AngleRcP1));
          hFinalAngle->Fill(AngleRcP1);
          hFinalCos->Fill(cos(AngleRcP1));
      }
    }
    else if(SmoothSampling==1){
      gROOT->cd();
      TF1* fRandomAngle = MakeSmoothAngularSourceDisto(hAngle);
      //TRandom3* ran_gen = new TRandom3(11);
      for(unsigned uEntry=0; uEntry<SmoothEntires; uEntry++){
          Tau1 = 4.69;
          Tau2 = 0;
          fM1 = 1462.93;
          fP1 = MeanP1;
          RanVal1 = RanGen.Exponential(fM1/(fP1*Tau1));
          AngleRcP1 = fRandomAngle->GetRandom(0.,TMath::Pi());
          MagicSource.AddBGT_RP(RanVal1,cos(AngleRcP1));
          hFinalAngle->Fill(AngleRcP1);
          hFinalCos->Fill(cos(AngleRcP1));
      }
      gROOT->cd();
      //delete ran_gen;
      delete fRandomAngle;
    }
    else{
      gROOT->cd();
      TF1* fRandomCos = MakeSmoothAngularCosine(hCos);
      //TRandom3* ran_gen = new TRandom3(11);
      for(unsigned uEntry=0; uEntry<SmoothEntires; uEntry++){
          Tau1 = 4.69;
          Tau2 = 0;
          fM1 = 1462.93;
          fP1 = MeanP1;
          RanVal1 = RanGen.Exponential(fM1/(fP1*Tau1));
          RanCos = fRandomCos->GetRandom(-1.,1.);
          MagicSource.AddBGT_RP(RanVal1,RanCos);
          hFinalAngle->Fill(acos(RanCos));
          hFinalCos->Fill(RanCos);
      }
      gROOT->cd();
      delete fRandomCos;
    }

    delete F_EposDisto_pReso_Xim;

    //if you have resonances contributing to both particles, we need to repeat the above procedure
    //for the prim-reso (AddBGT_PR) and reso-reso (AddBGT_RR) cases

    const unsigned NumSourceBins = 128;
    const double rMin = 0;
    const double rMax = 16;
    TFile* fOutput = new TFile(TString::Format("%s/OtherTasks/Georgios_LXi_ResoTest/fOutput_%i_%.2f.root",GetFemtoOutputFolder(),int(SmoothSampling),CoreSize),"recreate");
    TH1F* hSource = new TH1F("hSource","hSource",NumSourceBins,rMin,rMax);

    //fill the histo fro the source
    for(unsigned uBin=0; uBin<NumSourceBins; uBin++){
      //get the x-axis (r value) of the current bin
      double xaxis = hSource->GetBinCenter(uBin+1);
      //an array for the parameters, [0] is source size, [1] is == 2 (for a Gaussian)
      double parameters[2];
      parameters[0] = CoreSize;
      parameters[1] = 2.0;
      double SourceValue = MagicSource.RootEval(&xaxis, parameters);
      hSource->SetBinContent(uBin+1,SourceValue);
      //infinite errors for now
      hSource->SetBinError(uBin+1,1000.);
    }
//printf("4\n");
    //idea: fit the source distribution only in a range around its peak
    //to do this: silly idea: put very large uncertainties in the bins outside of this range
    //we can get this range automatically, by evaluating the central (median) integral of the source distribution
    //with this set up, we fit the 68% most central yield of the source distribution
    double lowerlimit;
    double upperlimit;
    GetCentralInterval(*hSource, 0.84, lowerlimit, upperlimit, true);
    unsigned lowerbin = hSource->FindBin(lowerlimit);
    unsigned upperbin = hSource->FindBin(upperlimit);
    for(unsigned uBin=lowerbin; uBin<=upperbin; uBin++){
      hSource->SetBinError(uBin+1,0.01);
    }

    printf("Core size of %.3f fm\n",CoreSize);
    printf("The fit will be performed in the range [%.2f, %.2f] fm\n",lowerlimit,upperlimit);
    //fyi, GaussSourceTF1 is in DLM_Source.h if you want to check it out.
    TF1* fSource = new TF1("fSource",GaussSourceTF1,rMin,rMax,1);
    fSource->SetParameter(0,CoreSize);
    fSource->SetParLimits(0,CoreSize*0.5,CoreSize*2.0);
    hSource->Fit(fSource,"S, N, R, M");
    printf("The effective Gaussian size is %.3f +/- %.3f fm\n",fSource->GetParameter(0),fSource->GetParError(0));

    //get rid of weird plotting
    for(unsigned uBin=0; uBin<NumSourceBins; uBin++){
      hSource->SetBinError(uBin+1,0.01);
    }
    hSource->Write();
    fSource->Write();
    hAngle->Scale(1./hAngle->Integral(),"width");
    hAngle->Write();
    hCos->Write();
    hFinalAngle->Scale(1./hFinalAngle->Integral(),"width");
    hFinalAngle->Write();
    hFinalCos->Scale(1./hFinalCos->Integral(),"width");
    hFinalCos->Write();

    //gROOT->cd();
    //delete hAngle;

    delete hSource;
    delete fSource;
    delete hAngle;
    delete hCos;
    delete hFinalAngle;
    delete hFinalCos;
    delete fOutput;
}


void EmmaDaniel_piD(const double CoreSize=0.82, const int& SmoothSampling=0, bool IncludeBug = false){
  printf("pi-D\n");
    const unsigned SmoothEntires = 1024*128;

    //DLM_CleverMcLevyResoTM* MagicSource = new DLM_CleverMcLevyResoTM ();
    DLM_CleverMcLevyResoTM MagicSource;

    //DO NOT CHANGE !!! Sets up numerical bullshit, tuned for a Gaussian source
    MagicSource.InitStability(1,2-1e-6,2+1e-6);
    MagicSource.InitScale(38,0.15,2.0);
    MagicSource.InitRad(257*2,0,64);
    MagicSource.InitType(2);
    ///////////////////

    //the amount of secondaries
    //pions (64.22%)
    MagicSource.SetUpReso(0,0.0);
    //than for the Xis, here its 0% (we have ONLY primordials)
    MagicSource.SetUpReso(1,0.682);
    MagicSource.InitNumMcIter(1000000);

    //the cut off scale in k*, for which the angular distributions from EPOS
    //are evaluated. 200 MeV works okay, you can go up to 300 MeV for systematic checks
    const double k_CutOff = 300;

    //to be used for the NTuple later on
    Float_t k_D;
    Float_t fP1;
    Float_t fP2;
    Float_t fM1;
    Float_t fM2;
    Float_t Tau1;
    Float_t Tau2;
    Float_t AngleRcP1;
    Float_t AngleRcP2;
    Float_t AngleP1P2;
    //random generator dimi style. The input is incompatible with the ROOT random generator,
    //do not mix and match, do not ask me how I know this. Ask Bernie.
    //11 is the seed, you can change that to you favorite number
    DLM_Random RanGen(11);
    //dummies to save random shit
    double RanVal1;
    double RanVal2;
    double RanVal3;
    double RanCos;
    double MeanP2=0;
//printf("1\n");
    //open the magic file from dimi with the angular distributions.
    TFile* F_EposDisto = new TFile(TString::Format("%s/CatsFiles/Source/EposAngularDist/TEMP/SecondTry/ALL_D_piReso.root",GetCernBoxDimi()));
    //set up the ntuple, do not change anything unless told so by dimi
    TNtuple* T_EposDisto = (TNtuple*)F_EposDisto->Get("InfoTuple_ClosePairs");
    unsigned N_EposDisto = T_EposDisto->GetEntries();
    T_EposDisto->SetBranchAddress("k_D",&k_D);
    T_EposDisto->SetBranchAddress("P1",&fP1);
    T_EposDisto->SetBranchAddress("P2",&fP2);
    T_EposDisto->SetBranchAddress("M1",&fM1);
    T_EposDisto->SetBranchAddress("M2",&fM2);
    T_EposDisto->SetBranchAddress("Tau1",&Tau1);
    T_EposDisto->SetBranchAddress("Tau2",&Tau2);
    T_EposDisto->SetBranchAddress("AngleRcP1",&AngleRcP1);
    T_EposDisto->SetBranchAddress("AngleRcP2",&AngleRcP2);
    T_EposDisto->SetBranchAddress("AngleP1P2",&AngleP1P2);

    //TH1F* hAngleRcP1 = new TH1F("hAngleRcP1","hAngleRcP1",1024,0,TMath::Pi());

    gROOT->cd();
    TH1F* hAngle = new TH1F("hAngle","hAngle",32,0,TMath::Pi());
    TH1F* hCos = new TH1F("hCos","hCos",32,-1.,1.);
    TH1F* hFinalAngle = new TH1F("hFinalAngle","hFinalAngle",32,0,TMath::Pi());
    TH1F* hFinalCos = new TH1F("hFinalCos","hFinalCos",32,-1.,1.);
    F_EposDisto->cd();
    int NumUsefulEntries = 0;
    for(unsigned uEntry=0; uEntry<N_EposDisto; uEntry++){
      T_EposDisto->GetEntry(uEntry);
      if(k_D>k_CutOff) continue;
      hAngle->Fill(AngleRcP2);
      hCos->Fill(cos(AngleRcP2));
      MeanP2 += fP2;
      NumUsefulEntries++;
    }
    MeanP2 /= double(NumUsefulEntries);

    //hAngle->SetBinContent(hAngle->FindBin(0.275),0.31);
    hCos->Scale(1./hAngle->Integral(),"width");
//printf("2\n");

    if(!SmoothSampling){
      //iterate over the ntuple
      for(unsigned uEntry=0; uEntry<N_EposDisto; uEntry++){
          //get each entry
          T_EposDisto->GetEntry(uEntry);
          //disregard the entry of you are outside the desired k*
          if(k_D>k_CutOff) continue;
          //overwrite the value for the lifetime. This is computed from the
          //stat. hadronization model (Vale) or thermal fist (Max)
          //this is the value for the secondary protons
          Tau1 = 0;
          //for primoridials (the Xis) we put 0
          Tau2 = 1.50;
          //put in the average mass of the resonances (again from SHM or TF)
          //this is the value for protons
          fM2 = 1180;
          //generate a random path length for the propagation of the resonances
          //nothing to change!
          RanVal2 = RanGen.Exponential(fM2/(fP2*Tau2));
          //adds a single entry into the PDF for the angular distribution to be used
          if(!IncludeBug) MagicSource.AddBGT_PR(RanVal2,cos(AngleRcP2));
          else MagicSource.AddBGT_PR(RanVal2,-cos(AngleRcP2));
          hFinalAngle->Fill(AngleRcP2);
          hFinalCos->Fill(cos(AngleRcP2));
      }
    }
    else if(SmoothSampling==1){
      gROOT->cd();
      TF1* fRandomAngle = MakeSmoothAngularSourceDisto(hAngle);
      //TRandom3* ran_gen = new TRandom3(11);
      for(unsigned uEntry=0; uEntry<SmoothEntires; uEntry++){
          Tau1 = 0;
          Tau2 = 1.50;
          fM2 = 1180;
          fP2 = MeanP2;
          RanVal2 = RanGen.Exponential(fM2/(fP2*Tau2));
          AngleRcP2 = fRandomAngle->GetRandom(0.,TMath::Pi());
          if(!IncludeBug) MagicSource.AddBGT_PR(RanVal2,cos(AngleRcP2));
          else MagicSource.AddBGT_PR(RanVal2,-cos(AngleRcP2));
          hFinalAngle->Fill(AngleRcP2);
          hFinalCos->Fill(cos(AngleRcP2));
      }
      gROOT->cd();
      //delete ran_gen;
      delete fRandomAngle;
    }
    else{
      gROOT->cd();
      TF1* fRandomCos = MakeSmoothAngularCosine(hCos);
      //TRandom3* ran_gen = new TRandom3(11);
      for(unsigned uEntry=0; uEntry<SmoothEntires; uEntry++){
          Tau1 = 0;
          Tau2 = 1.50;
          fM2 = 1180;
          fP2 = MeanP2;
          RanVal2 = RanGen.Exponential(fM2/(fP2*Tau2));
          RanCos = fRandomCos->GetRandom(-1.,1.);
          if(!IncludeBug) MagicSource.AddBGT_PR(RanVal2,RanCos);
          else MagicSource.AddBGT_PR(RanVal2,-RanCos);
          hFinalAngle->Fill(acos(RanCos));
          hFinalCos->Fill(RanCos);
      }
      gROOT->cd();
      delete fRandomCos;
    }

    delete F_EposDisto;

    //if you have resonances contributing to both particles, we need to repeat the above procedure
    //for the prim-reso (AddBGT_PR) and reso-reso (AddBGT_RR) cases

    const unsigned NumSourceBins = 128;
    const double rMin = 0;
    const double rMax = 16;
    TFile* fOutput = new TFile(TString::Format("%s/OtherTasks/EmmaDaniel_Dmesons/fSecond%s_piD_%i_%.2f.root",GetFemtoOutputFolder(),IncludeBug?"BUG":"",
    int(SmoothSampling),CoreSize),"recreate");
    TH1F* hSource = new TH1F("hSource","hSource",NumSourceBins,rMin,rMax);

    //fill the histo fro the source
    for(unsigned uBin=0; uBin<NumSourceBins; uBin++){
      //get the x-axis (r value) of the current bin
      double xaxis = hSource->GetBinCenter(uBin+1);
      //an array for the parameters, [0] is source size, [1] is == 2 (for a Gaussian)
      double parameters[2];
      parameters[0] = CoreSize;
      parameters[1] = 2.0;
      double SourceValue = MagicSource.RootEval(&xaxis, parameters);
      hSource->SetBinContent(uBin+1,SourceValue);
      //infinite errors for now
      hSource->SetBinError(uBin+1,1000.);
    }
//printf("4\n");
    //idea: fit the source distribution only in a range around its peak
    //to do this: silly idea: put very large uncertainties in the bins outside of this range
    //we can get this range automatically, by evaluating the central (median) integral of the source distribution
    //with this set up, we fit the 68% most central yield of the source distribution
    double lowerlimit;
    double upperlimit;
    GetCentralInterval(*hSource, 0.84, lowerlimit, upperlimit, true);
    unsigned lowerbin = hSource->FindBin(lowerlimit);
    unsigned upperbin = hSource->FindBin(upperlimit);
    for(unsigned uBin=lowerbin; uBin<=upperbin; uBin++){
      hSource->SetBinError(uBin+1,0.01);
    }

    printf("Core size of %.3f fm\n",CoreSize);
    printf("The fit will be performed in the range [%.2f, %.2f] fm\n",lowerlimit,upperlimit);
    //fyi, GaussSourceTF1 is in DLM_Source.h if you want to check it out.
    TF1* fSource = new TF1("fSource",GaussSourceTF1,rMin,rMax,1);
    fSource->SetParameter(0,CoreSize);
    fSource->SetParLimits(0,CoreSize*0.5,CoreSize*2.0);
    hSource->Fit(fSource,"Q, S, N, R, M");
    printf("The effective Gaussian size is %.3f +/- %.3f fm\n",fSource->GetParameter(0),fSource->GetParError(0));

    //get rid of weird plotting
    for(unsigned uBin=0; uBin<NumSourceBins; uBin++){
      hSource->SetBinError(uBin+1,0.01);
    }
    hSource->Write();
    fSource->Write();
    hAngle->Scale(1./hAngle->Integral(),"width");
    hAngle->Write();
    hCos->Write();
    hFinalAngle->Scale(1./hFinalAngle->Integral(),"width");
    hFinalAngle->Write();
    hFinalCos->Scale(1./hFinalCos->Integral(),"width");
    hFinalCos->Write();

    //gROOT->cd();
    //delete hAngle;

    delete hSource;
    delete fSource;
    delete hAngle;
    delete hCos;
    delete hFinalAngle;
    delete hFinalCos;
    delete fOutput;
    printf("-------------------\n");
}


void EmmaDaniel_KD(const double CoreSize=0.81, const int& SmoothSampling=0, bool IncludeBug = false){
  printf("K-D\n");
    const unsigned SmoothEntires = 1024*128;

    //DLM_CleverMcLevyResoTM* MagicSource = new DLM_CleverMcLevyResoTM ();
    DLM_CleverMcLevyResoTM MagicSource;

    //DO NOT CHANGE !!! Sets up numerical bullshit, tuned for a Gaussian source
    MagicSource.InitStability(1,2-1e-6,2+1e-6);
    MagicSource.InitScale(38,0.15,2.0);
    MagicSource.InitRad(257*2,0,64);
    MagicSource.InitType(2);
    ///////////////////

    //the amount of secondaries
    //pions (64.22%)
    MagicSource.SetUpReso(0,0.0);
    //than for the Xis, here its 0% (we have ONLY primordials)
    MagicSource.SetUpReso(1,0.476);
    MagicSource.InitNumMcIter(1000000);

    //the cut off scale in k*, for which the angular distributions from EPOS
    //are evaluated. 200 MeV works okay, you can go up to 300 MeV for systematic checks
    const double k_CutOff = 300;

    //to be used for the NTuple later on
    Float_t k_D;
    Float_t fP1;
    Float_t fP2;
    Float_t fM1;
    Float_t fM2;
    Float_t Tau1;
    Float_t Tau2;
    Float_t AngleRcP1;
    Float_t AngleRcP2;
    Float_t AngleP1P2;
    //random generator dimi style. The input is incompatible with the ROOT random generator,
    //do not mix and match, do not ask me how I know this. Ask Bernie.
    //11 is the seed, you can change that to you favorite number
    DLM_Random RanGen(11);
    //dummies to save random shit
    double RanVal1;
    double RanVal2;
    double RanVal3;
    double RanCos;
    double MeanP2=0;
//printf("1\n");
    //open the magic file from dimi with the angular distributions.
    TFile* F_EposDisto = new TFile(TString::Format("%s/CatsFiles/Source/EposAngularDist/TEMP/SecondTry/ALL_D_KaonReso.root",GetCernBoxDimi()));
    //set up the ntuple, do not change anything unless told so by dimi
    TNtuple* T_EposDisto = (TNtuple*)F_EposDisto->Get("InfoTuple_ClosePairs");
    unsigned N_EposDisto = T_EposDisto->GetEntries();
    T_EposDisto->SetBranchAddress("k_D",&k_D);
    T_EposDisto->SetBranchAddress("P1",&fP1);
    T_EposDisto->SetBranchAddress("P2",&fP2);
    T_EposDisto->SetBranchAddress("M1",&fM1);
    T_EposDisto->SetBranchAddress("M2",&fM2);
    T_EposDisto->SetBranchAddress("Tau1",&Tau1);
    T_EposDisto->SetBranchAddress("Tau2",&Tau2);
    T_EposDisto->SetBranchAddress("AngleRcP1",&AngleRcP1);
    T_EposDisto->SetBranchAddress("AngleRcP2",&AngleRcP2);
    T_EposDisto->SetBranchAddress("AngleP1P2",&AngleP1P2);

    //TH1F* hAngleRcP1 = new TH1F("hAngleRcP1","hAngleRcP1",1024,0,TMath::Pi());

    gROOT->cd();
    TH1F* hAngle = new TH1F("hAngle","hAngle",32,0,TMath::Pi());
    TH1F* hCos = new TH1F("hCos","hCos",32,-1.,1.);
    TH1F* hFinalAngle = new TH1F("hFinalAngle","hFinalAngle",32,0,TMath::Pi());
    TH1F* hFinalCos = new TH1F("hFinalCos","hFinalCos",32,-1.,1.);
    F_EposDisto->cd();
    int NumUsefulEntries = 0;
    for(unsigned uEntry=0; uEntry<N_EposDisto; uEntry++){
      T_EposDisto->GetEntry(uEntry);
      if(k_D>k_CutOff) continue;
      hAngle->Fill(AngleRcP2);
      hCos->Fill(cos(AngleRcP2));
      MeanP2 += fP2;
      NumUsefulEntries++;
    }
    MeanP2 /= double(NumUsefulEntries);

    //hAngle->SetBinContent(hAngle->FindBin(0.275),0.31);
    hCos->Scale(1./hAngle->Integral(),"width");
//printf("2\n");

    if(!SmoothSampling){
      //iterate over the ntuple
      for(unsigned uEntry=0; uEntry<N_EposDisto; uEntry++){
          //get each entry
          T_EposDisto->GetEntry(uEntry);
          //disregard the entry of you are outside the desired k*
          if(k_D>k_CutOff) continue;
          //overwrite the value for the lifetime. This is computed from the
          //stat. hadronization model (Vale) or thermal fist (Max)
          //this is the value for the secondary protons
          Tau1 = 0;
          //for primoridials (the Xis) we put 0
          Tau2 = 3.66;
          //put in the average mass of the resonances (again from SHM or TF)
          //this is the value for protons
          fM2 = 1054;
          //generate a random path length for the propagation of the resonances
          //nothing to change!
          RanVal2 = RanGen.Exponential(fM2/(fP2*Tau2));
          //adds a single entry into the PDF for the angular distribution to be used
          if(!IncludeBug) MagicSource.AddBGT_PR(RanVal2,cos(AngleRcP2));
          else MagicSource.AddBGT_PR(RanVal2,-cos(AngleRcP2));
          hFinalAngle->Fill(AngleRcP2);
          hFinalCos->Fill(cos(AngleRcP2));
      }
    }
    else if(SmoothSampling==1){
      gROOT->cd();
      TF1* fRandomAngle = MakeSmoothAngularSourceDisto(hAngle);
      //TRandom3* ran_gen = new TRandom3(11);
      for(unsigned uEntry=0; uEntry<SmoothEntires; uEntry++){
          Tau1 = 0;
          Tau2 = 3.66;
          fM2 = 1054;
          fP2 = MeanP2;
          RanVal2 = RanGen.Exponential(fM2/(fP2*Tau2));
          AngleRcP2 = fRandomAngle->GetRandom(0.,TMath::Pi());
          if(!IncludeBug) MagicSource.AddBGT_PR(RanVal2,cos(AngleRcP2));
          else MagicSource.AddBGT_PR(RanVal2,-cos(AngleRcP2));
          hFinalAngle->Fill(AngleRcP2);
          hFinalCos->Fill(cos(AngleRcP2));
      }
      gROOT->cd();
      //delete ran_gen;
      delete fRandomAngle;
    }
    else{
      gROOT->cd();
      TF1* fRandomCos = MakeSmoothAngularCosine(hCos);
      //TRandom3* ran_gen = new TRandom3(11);
      for(unsigned uEntry=0; uEntry<SmoothEntires; uEntry++){
          Tau1 = 0;
          Tau2 = 3.66;
          fM2 = 1054;
          fP2 = MeanP2;
          RanVal2 = RanGen.Exponential(fM2/(fP2*Tau2));
          RanCos = fRandomCos->GetRandom(-1.,1.);
          if(!IncludeBug) MagicSource.AddBGT_PR(RanVal2,RanCos);
          else MagicSource.AddBGT_PR(RanVal2,RanCos);
          hFinalAngle->Fill(acos(RanCos));
          hFinalCos->Fill(RanCos);
      }
      gROOT->cd();
      delete fRandomCos;
    }

    delete F_EposDisto;

    //if you have resonances contributing to both particles, we need to repeat the above procedure
    //for the prim-reso (AddBGT_PR) and reso-reso (AddBGT_RR) cases

    const unsigned NumSourceBins = 128;
    const double rMin = 0;
    const double rMax = 16;
    TFile* fOutput = new TFile(TString::Format("%s/OtherTasks/EmmaDaniel_Dmesons/fSecond%s_KaonD_%i_%.2f.root",GetFemtoOutputFolder(),IncludeBug?"BUG":"",
    int(SmoothSampling),CoreSize),"recreate");
    TH1F* hSource = new TH1F("hSource","hSource",NumSourceBins,rMin,rMax);

    //fill the histo fro the source
    for(unsigned uBin=0; uBin<NumSourceBins; uBin++){
      //get the x-axis (r value) of the current bin
      double xaxis = hSource->GetBinCenter(uBin+1);
      //an array for the parameters, [0] is source size, [1] is == 2 (for a Gaussian)
      double parameters[2];
      parameters[0] = CoreSize;
      parameters[1] = 2.0;
      double SourceValue = MagicSource.RootEval(&xaxis, parameters);
      hSource->SetBinContent(uBin+1,SourceValue);
      //infinite errors for now
      hSource->SetBinError(uBin+1,1000.);
    }
//printf("4\n");
    //idea: fit the source distribution only in a range around its peak
    //to do this: silly idea: put very large uncertainties in the bins outside of this range
    //we can get this range automatically, by evaluating the central (median) integral of the source distribution
    //with this set up, we fit the 68% most central yield of the source distribution
    double lowerlimit;
    double upperlimit;
    GetCentralInterval(*hSource, 0.84, lowerlimit, upperlimit, true);
    unsigned lowerbin = hSource->FindBin(lowerlimit);
    unsigned upperbin = hSource->FindBin(upperlimit);
    for(unsigned uBin=lowerbin; uBin<=upperbin; uBin++){
      hSource->SetBinError(uBin+1,0.01);
    }

    printf("Core size of %.3f fm\n",CoreSize);
    printf("The fit will be performed in the range [%.2f, %.2f] fm\n",lowerlimit,upperlimit);
    //fyi, GaussSourceTF1 is in DLM_Source.h if you want to check it out.
    TF1* fSource = new TF1("fSource",GaussSourceTF1,rMin,rMax,1);
    fSource->SetParameter(0,CoreSize);
    fSource->SetParLimits(0,CoreSize*0.5,CoreSize*2.0);
    hSource->Fit(fSource,"Q, S, N, R, M");
    printf("The effective Gaussian size is %.3f +/- %.3f fm\n",fSource->GetParameter(0),fSource->GetParError(0));

    //get rid of weird plotting
    for(unsigned uBin=0; uBin<NumSourceBins; uBin++){
      hSource->SetBinError(uBin+1,0.01);
    }
    hSource->Write();
    fSource->Write();
    hAngle->Scale(1./hAngle->Integral(),"width");
    hAngle->Write();
    hCos->Write();
    hFinalAngle->Scale(1./hFinalAngle->Integral(),"width");
    hFinalAngle->Write();
    hFinalCos->Scale(1./hFinalCos->Integral(),"width");
    hFinalCos->Write();

    //gROOT->cd();
    //delete hAngle;

    delete hSource;
    delete fSource;
    delete hAngle;
    delete hCos;
    delete hFinalAngle;
    delete hFinalCos;
    delete fOutput;

    printf("-------------------\n");
}




double MakeMaxRamonaFit(double* x, double* pars){
  return pars[1]*GaussSourceTF1(x,pars)+(1.-pars[1]);
}

//0 = Original
//1 = fit the angle with splines
//2 = fit the cos with 1/(1-exp)
void MaxRamona_piSig_ResoTest(){

    const double CoreSize = 0.8;
    const unsigned SmoothEntires = 1024*128;

    //DLM_CleverMcLevyResoTM* MagicSource = new DLM_CleverMcLevyResoTM ();
    DLM_CleverMcLevyResoTM MagicSource;

    //DO NOT CHANGE !!! Sets up numerical bullshit, tuned for a Gaussian source
    MagicSource.InitStability(1,2-1e-6,2+1e-6);
    MagicSource.InitScale(38,0.15,2.0);
    MagicSource.InitRad(257*2,0,64);
    MagicSource.InitType(2);
    ///////////////////

    //for p-Xi, set up the amount of secondaries
    //first for the protons (64.22%)
    //MagicSource.SetUpReso(0,0.6422);
    MagicSource.SetUpReso(0,0.61);
    //MagicSource.SetUpReso(0,0.6422);
    //than for the Xis, here its 0% (we have ONLY primordials)
    MagicSource.SetUpReso(1,0.6438);
    //MagicSource.SetUpReso(1,0.6422);
    MagicSource.InitNumMcIter(1000000);

    //the cut off scale in k*, for which the angular distributions from EPOS
    //are evaluated. 200 MeV works okay, you can go up to 300 MeV for systematic checks
    const double k_CutOff = 250;

    //to be used for the NTuple later on
    Float_t k_D;
    Float_t fP1;
    Float_t fP2;
    Float_t fM1;
    Float_t fM2;
    Float_t Tau1;
    Float_t Tau2;
    Float_t AngleRcP1;
    Float_t AngleRcP2;
    Float_t AngleP1P2;
    //random generator dimi style. The input is incompatible with the ROOT random generator,
    //do not mix and match, do not ask me how I know this. Ask Bernie.
    //11 is the seed, you can change that to you favorite number
    DLM_Random RanGen(11);
    //dummies to save random shit
    double RanVal1;
    double RanVal2;
    double RanVal3;
    double RanCos;
    double MeanP1=0;

    TFile* F_EposDisto_piReso_Sig = new TFile(TString::Format("%s/CatsFiles/Source/EposAngularDist/ForMaxRamona_piReso_Sig.root",GetCernBoxDimi()));
    //TFile* F_EposDisto_piReso_Sig = new TFile(TString::Format("%s/CatsFiles/Source/EposAngularDist/EposDisto_p_pReso.root",GetCernBoxDimi()));
    TNtuple* T_EposDisto_piReso_Sig = (TNtuple*)F_EposDisto_piReso_Sig->Get("InfoTuple_ClosePairs");
    unsigned N_EposDisto_piReso_Sig = T_EposDisto_piReso_Sig->GetEntries();
    T_EposDisto_piReso_Sig->SetBranchAddress("k_D",&k_D);
    T_EposDisto_piReso_Sig->SetBranchAddress("P1",&fP1);
    T_EposDisto_piReso_Sig->SetBranchAddress("P2",&fP2);
    T_EposDisto_piReso_Sig->SetBranchAddress("M1",&fM1);
    T_EposDisto_piReso_Sig->SetBranchAddress("M2",&fM2);
    T_EposDisto_piReso_Sig->SetBranchAddress("Tau1",&Tau1);
    T_EposDisto_piReso_Sig->SetBranchAddress("Tau2",&Tau2);
    T_EposDisto_piReso_Sig->SetBranchAddress("AngleRcP1",&AngleRcP1);
    T_EposDisto_piReso_Sig->SetBranchAddress("AngleRcP2",&AngleRcP2);
    T_EposDisto_piReso_Sig->SetBranchAddress("AngleP1P2",&AngleP1P2);
    gROOT->cd();
    //iterate over the ntuple
    for(unsigned uEntry=0; uEntry<N_EposDisto_piReso_Sig; uEntry++){
        //get each entry
        T_EposDisto_piReso_Sig->GetEntry(uEntry);
        //disregard the entry of you are outside the desired k*
        if(k_D>k_CutOff) continue;
        Tau1 = 1.5;
        Tau2 = 0;
        fM1 = 1124.0;
        RanVal1 = RanGen.Exponential(fM1/(fP1*Tau1));
        MagicSource.AddBGT_RP(RanVal1,cos(AngleRcP1));

        //Tau1 = 0;
        //Tau2 = 1.65;
        //fM2 = 1362;
        //RanVal1 = RanGen.Exponential(fM2/(fP2*Tau2));
        //MagicSource.AddBGT_PR(RanVal1,-cos(AngleRcP2));
        //MagicSource.AddBGT_RP(RanVal1,cos(AngleRcP2));
    }
    delete T_EposDisto_piReso_Sig;
    delete F_EposDisto_piReso_Sig;



    TFile* F_EposDisto_pi_SigReso = new TFile(TString::Format("%s/CatsFiles/Source/EposAngularDist/ForMaxRamona_pi_SigReso.root",GetCernBoxDimi()));
    TNtuple* T_EposDisto_pi_SigReso = (TNtuple*)F_EposDisto_pi_SigReso->Get("InfoTuple_ClosePairs");
    unsigned N_EposDisto_pi_SigReso = T_EposDisto_pi_SigReso->GetEntries();
    T_EposDisto_pi_SigReso->SetBranchAddress("k_D",&k_D);
    T_EposDisto_pi_SigReso->SetBranchAddress("P1",&fP1);
    T_EposDisto_pi_SigReso->SetBranchAddress("P2",&fP2);
    T_EposDisto_pi_SigReso->SetBranchAddress("M1",&fM1);
    T_EposDisto_pi_SigReso->SetBranchAddress("M2",&fM2);
    T_EposDisto_pi_SigReso->SetBranchAddress("Tau1",&Tau1);
    T_EposDisto_pi_SigReso->SetBranchAddress("Tau2",&Tau2);
    T_EposDisto_pi_SigReso->SetBranchAddress("AngleRcP1",&AngleRcP1);
    T_EposDisto_pi_SigReso->SetBranchAddress("AngleRcP2",&AngleRcP2);
    T_EposDisto_pi_SigReso->SetBranchAddress("AngleP1P2",&AngleP1P2);
    gROOT->cd();
    //iterate over the ntuple
    for(unsigned uEntry=0; uEntry<N_EposDisto_pi_SigReso; uEntry++){
        //get each entry
        T_EposDisto_pi_SigReso->GetEntry(uEntry);
        //disregard the entry of you are outside the desired k*
        if(k_D>k_CutOff) continue;
        Tau1 = 0;
        Tau2 = 4.69;
        fM2 = 1462.93;
        RanVal2 = RanGen.Exponential(fM2/(fP2*Tau2));
        MagicSource.AddBGT_PR(RanVal2,-cos(AngleRcP2));
    }
    delete T_EposDisto_pi_SigReso;
    delete F_EposDisto_pi_SigReso;


    TFile* F_EposDisto_piReso_SigReso = new TFile(TString::Format("%s/CatsFiles/Source/EposAngularDist/ForMaxRamona_piReso_SigReso.root",GetCernBoxDimi()));
    //TFile* F_EposDisto_piReso_SigReso = new TFile(TString::Format("%s/CatsFiles/Source/EposAngularDist/EposDisto_pReso_dReso.root",GetCernBoxDimi()));

    TNtuple* T_EposDisto_piReso_SigReso = (TNtuple*)F_EposDisto_piReso_SigReso->Get("InfoTuple_ClosePairs");
    unsigned N_EposDisto_piReso_SigReso = T_EposDisto_piReso_SigReso->GetEntries();
    T_EposDisto_piReso_SigReso->SetBranchAddress("k_D",&k_D);
    T_EposDisto_piReso_SigReso->SetBranchAddress("P1",&fP1);
    T_EposDisto_piReso_SigReso->SetBranchAddress("P2",&fP2);
    T_EposDisto_piReso_SigReso->SetBranchAddress("M1",&fM1);
    T_EposDisto_piReso_SigReso->SetBranchAddress("M2",&fM2);
    T_EposDisto_piReso_SigReso->SetBranchAddress("Tau1",&Tau1);
    T_EposDisto_piReso_SigReso->SetBranchAddress("Tau2",&Tau2);
    T_EposDisto_piReso_SigReso->SetBranchAddress("AngleRcP1",&AngleRcP1);
    T_EposDisto_piReso_SigReso->SetBranchAddress("AngleRcP2",&AngleRcP2);
    T_EposDisto_piReso_SigReso->SetBranchAddress("AngleP1P2",&AngleP1P2);
    gROOT->cd();
    //iterate over the ntuple
    for(unsigned uEntry=0; uEntry<N_EposDisto_piReso_SigReso; uEntry++){
        //get each entry
        T_EposDisto_piReso_SigReso->GetEntry(uEntry);
        //disregard the entry of you are outside the desired k*
        if(k_D>k_CutOff) continue;
        Tau1 = 1.5;
        fM1 = 1124.0;
        RanVal1 = RanGen.Exponential(fM1/(fP1*Tau1));
        Tau2 = 4.69;
        fM2 = 1462.93;
        RanVal2 = RanGen.Exponential(fM2/(fP2*Tau2));
        //check the signs
        MagicSource.AddBGT_RR(RanVal1,cos(AngleRcP1),RanVal2,cos(AngleRcP2),cos(AngleP1P2));

        //Tau1 = 1.65;
        //Tau2 = 1.65;
        //fM1 = 1362;
        //fM2 = 1362;
        //RanVal1 = RanGen.Exponential(fM1/(fP1*Tau1));
        //RanVal2 = RanGen.Exponential(fM2/(fP2*Tau2));
        //MagicSource.AddBGT_RR(RanVal1,cos(AngleRcP1),RanVal2,cos(AngleRcP2),cos(AngleP1P2));

    }
    delete T_EposDisto_piReso_SigReso;
    delete F_EposDisto_piReso_SigReso;


    const unsigned NumSourceBins = 128;
    const double rMin = 0;
    const double rMax = 16;
    TFile* fOutput = new TFile(TString::Format("%s/OtherTasks/MaxRamona_piSig_ResoTest/fOutput.root",GetFemtoOutputFolder()),"recreate");
    TH1F* hSource = new TH1F("hSource","hSource",NumSourceBins,rMin,rMax);

    //fill the histo fro the source
    for(unsigned uBin=0; uBin<NumSourceBins; uBin++){
      //get the x-axis (r value) of the current bin
      double xaxis = hSource->GetBinCenter(uBin+1);
      //an array for the parameters, [0] is source size, [1] is == 2 (for a Gaussian)
      double parameters[2];
      parameters[0] = CoreSize;
      parameters[1] = 2.0;
      double SourceValue = MagicSource.RootEval(&xaxis, parameters);
      hSource->SetBinContent(uBin+1,SourceValue);
      //infinite errors for now
      hSource->SetBinError(uBin+1,1000.);
    }
    //idea: fit the source distribution only in a range around its peak
    //to do this: silly idea: put very large uncertainties in the bins outside of this range
    //we can get this range automatically, by evaluating the central (median) integral of the source distribution
    //with this set up, we fit the 68% most central yield of the source distribution
    double lowerlimit;
    double upperlimit;
    GetCentralInterval(*hSource, 0.84, lowerlimit, upperlimit, true);
    unsigned lowerbin = hSource->FindBin(lowerlimit);
    unsigned upperbin = hSource->FindBin(upperlimit);
    for(unsigned uBin=lowerbin; uBin<=upperbin; uBin++){
      hSource->SetBinError(uBin+1,0.01);
    }

    printf("Core size of %.3f fm\n",CoreSize);
    printf("The fit will be performed in the range [%.2f, %.2f] fm\n",lowerlimit,upperlimit);
    //fyi, GaussSourceTF1 is in DLM_Source.h if you want to check it out.
    TF1* fSource = new TF1("fSource",GaussSourceScaledTF1,rMin,rMax,2);
    //fSource->SetParameter(0,CoreSize);
    //fSource->SetParLimits(0,CoreSize*0.5,CoreSize*2.0);
    fSource->FixParameter(0,1.0);
    fSource->SetParameter(1,0.5);
    fSource->SetParLimits(1,0,1);
    //fSource->FixParameter(1,1.);
    hSource->Fit(fSource,"S, N, R, M");
    printf("The effective Gaussian size is %.3f +/- %.3f fm\n",fSource->GetParameter(0),fSource->GetParError(0));
    printf("The lambda par is %.1f +/- %.1f %%\n",fSource->GetParameter(1)*100.,fSource->GetParError(1)*100.);

    //get rid of weird plotting
    for(unsigned uBin=0; uBin<NumSourceBins; uBin++){
      hSource->SetBinError(uBin+1,0.01);
    }
    hSource->Write();
    fSource->Write();

    delete hSource;
    delete fSource;
    delete fOutput;
}


void Ramona_pK_ResoTest(){

      const double CoreSize = 0.8;

      //DLM_CleverMcLevyResoTM* MagicSource = new DLM_CleverMcLevyResoTM ();
      DLM_CleverMcLevyResoTM MagicSource;

      //DO NOT CHANGE !!! Sets up numerical bullshit, tuned for a Gaussian source
      MagicSource.InitStability(1,2-1e-6,2+1e-6);
      MagicSource.InitScale(38,0.15,2.0);
      MagicSource.InitRad(257*2,0,64);
      MagicSource.InitType(2);
      ///////////////////

      MagicSource.SetUpReso(0,0.6422);
      MagicSource.SetUpReso(1,0.4575);
      MagicSource.InitNumMcIter(1000000);

      //the cut off scale in k*, for which the angular distributions from EPOS
      //are evaluated. 200 MeV works okay, you can go up to 300 MeV for systematic checks
      const double k_CutOff = 200;

      //to be used for the NTuple later on
      Float_t k_D;
      Float_t fP1;
      Float_t fP2;
      Float_t fM1;
      Float_t fM2;
      Float_t Tau1;
      Float_t Tau2;
      Float_t AngleRcP1;
      Float_t AngleRcP2;
      Float_t AngleP1P2;
      //random generator dimi style. The input is incompatible with the ROOT random generator,
      //do not mix and match, do not ask me how I know this. Ask Bernie.
      //11 is the seed, you can change that to you favorite number
      DLM_Random RanGen(11);
      //dummies to save random shit
      double RanVal1;
      double RanVal2;
      double RanVal3;
      double RanCos;
      double MeanP1=0;

      TFile* F_EposDisto_pReso_Kaon = new TFile(TString::Format("%s/CatsFiles/Source/EposAngularDist/ForRamona_pReso_Kaon.root",GetCernBoxDimi()));
      TNtuple* T_EposDisto_pReso_Kaon = (TNtuple*)F_EposDisto_pReso_Kaon->Get("InfoTuple_ClosePairs");
      unsigned N_EposDisto_pReso_Kaon = T_EposDisto_pReso_Kaon->GetEntries();
      T_EposDisto_pReso_Kaon->SetBranchAddress("k_D",&k_D);
      T_EposDisto_pReso_Kaon->SetBranchAddress("P1",&fP1);
      T_EposDisto_pReso_Kaon->SetBranchAddress("P2",&fP2);
      T_EposDisto_pReso_Kaon->SetBranchAddress("M1",&fM1);
      T_EposDisto_pReso_Kaon->SetBranchAddress("M2",&fM2);
      T_EposDisto_pReso_Kaon->SetBranchAddress("Tau1",&Tau1);
      T_EposDisto_pReso_Kaon->SetBranchAddress("Tau2",&Tau2);
      T_EposDisto_pReso_Kaon->SetBranchAddress("AngleRcP1",&AngleRcP1);
      T_EposDisto_pReso_Kaon->SetBranchAddress("AngleRcP2",&AngleRcP2);
      T_EposDisto_pReso_Kaon->SetBranchAddress("AngleP1P2",&AngleP1P2);
      gROOT->cd();
      //iterate over the ntuple
      for(unsigned uEntry=0; uEntry<N_EposDisto_pReso_Kaon; uEntry++){
          //get each entry
          T_EposDisto_pReso_Kaon->GetEntry(uEntry);
          //disregard the entry of you are outside the desired k*
          if(k_D>k_CutOff) continue;
          Tau1 = 1.65;
          Tau2 = 0;
          fM1 = 1362;
          RanVal1 = RanGen.Exponential(fM1/(fP1*Tau1));
          MagicSource.AddBGT_RP(RanVal1,cos(AngleRcP1));
      }
      delete T_EposDisto_pReso_Kaon;
      delete F_EposDisto_pReso_Kaon;



      TFile* F_EposDisto_p_KaonReso = new TFile(TString::Format("%s/CatsFiles/Source/EposAngularDist/ForRamona_p_KaonReso.root",GetCernBoxDimi()));
      TNtuple* T_EposDisto_p_KaonReso = (TNtuple*)F_EposDisto_p_KaonReso->Get("InfoTuple_ClosePairs");
      unsigned N_EposDisto_p_KaonReso = T_EposDisto_p_KaonReso->GetEntries();
      T_EposDisto_p_KaonReso->SetBranchAddress("k_D",&k_D);
      T_EposDisto_p_KaonReso->SetBranchAddress("P1",&fP1);
      T_EposDisto_p_KaonReso->SetBranchAddress("P2",&fP2);
      T_EposDisto_p_KaonReso->SetBranchAddress("M1",&fM1);
      T_EposDisto_p_KaonReso->SetBranchAddress("M2",&fM2);
      T_EposDisto_p_KaonReso->SetBranchAddress("Tau1",&Tau1);
      T_EposDisto_p_KaonReso->SetBranchAddress("Tau2",&Tau2);
      T_EposDisto_p_KaonReso->SetBranchAddress("AngleRcP1",&AngleRcP1);
      T_EposDisto_p_KaonReso->SetBranchAddress("AngleRcP2",&AngleRcP2);
      T_EposDisto_p_KaonReso->SetBranchAddress("AngleP1P2",&AngleP1P2);
      gROOT->cd();
      //iterate over the ntuple
      for(unsigned uEntry=0; uEntry<N_EposDisto_p_KaonReso; uEntry++){
          //get each entry
          T_EposDisto_p_KaonReso->GetEntry(uEntry);
          //disregard the entry of you are outside the desired k*
          if(k_D>k_CutOff) continue;
          Tau1 = 0;
          Tau2 = 3.66;
          fM2 = 1054;
          RanVal2 = RanGen.Exponential(fM2/(fP2*Tau2));
          MagicSource.AddBGT_PR(RanVal2,cos(AngleRcP2));
      }
      delete T_EposDisto_p_KaonReso;
      delete F_EposDisto_p_KaonReso;


      TFile* F_EposDisto_pReso_KaonReso = new TFile(TString::Format("%s/CatsFiles/Source/EposAngularDist/ForRamona_pReso_KaonReso.root",GetCernBoxDimi()));
      //TFile* F_EposDisto_pReso_KaonReso = new TFile(TString::Format("%s/CatsFiles/Source/EposAngularDist/EposDisto_pReso_dReso.root",GetCernBoxDimi()));

      TNtuple* T_EposDisto_pReso_KaonReso = (TNtuple*)F_EposDisto_pReso_KaonReso->Get("InfoTuple_ClosePairs");
      unsigned N_EposDisto_pReso_KaonReso = T_EposDisto_pReso_KaonReso->GetEntries();
      T_EposDisto_pReso_KaonReso->SetBranchAddress("k_D",&k_D);
      T_EposDisto_pReso_KaonReso->SetBranchAddress("P1",&fP1);
      T_EposDisto_pReso_KaonReso->SetBranchAddress("P2",&fP2);
      T_EposDisto_pReso_KaonReso->SetBranchAddress("M1",&fM1);
      T_EposDisto_pReso_KaonReso->SetBranchAddress("M2",&fM2);
      T_EposDisto_pReso_KaonReso->SetBranchAddress("Tau1",&Tau1);
      T_EposDisto_pReso_KaonReso->SetBranchAddress("Tau2",&Tau2);
      T_EposDisto_pReso_KaonReso->SetBranchAddress("AngleRcP1",&AngleRcP1);
      T_EposDisto_pReso_KaonReso->SetBranchAddress("AngleRcP2",&AngleRcP2);
      T_EposDisto_pReso_KaonReso->SetBranchAddress("AngleP1P2",&AngleP1P2);
      gROOT->cd();
      //iterate over the ntuple
      for(unsigned uEntry=0; uEntry<N_EposDisto_pReso_KaonReso; uEntry++){
          //get each entry
          T_EposDisto_pReso_KaonReso->GetEntry(uEntry);
          //disregard the entry of you are outside the desired k*
          if(k_D>k_CutOff) continue;
          Tau1 = 1.65;
          fM1 = 1362;
          RanVal1 = RanGen.Exponential(fM1/(fP1*Tau1));
          Tau2 = 3.66;
          fM2 = 1054;
          RanVal2 = RanGen.Exponential(fM2/(fP2*Tau2));
          //check the signs
          MagicSource.AddBGT_RR(RanVal1,cos(AngleRcP1),RanVal2,cos(AngleRcP2),cos(AngleP1P2));
      }
      delete T_EposDisto_pReso_KaonReso;
      delete F_EposDisto_pReso_KaonReso;

      CATS Kitty;
      Kitty.SetMomBins(50,0,250);
      //Substitue all definitions regaring the source with the following lines:
      Kitty.SetAnaSource(CatsSourceForwarder, &MagicSource, 2);//use as it is, no questions asked
      Kitty.SetAnaSource(0,CoreSize);//size of the core
      Kitty.SetAnaSource(1,2);//important, but you do not care about it (specifies its a Gaussian core)
      Kitty.SetUseAnalyticSource(true);
      Kitty.SetNormalizedSource(false);//to avoid problems in case of a longer tail. Should not really matter for pK
      //... and so on, you continue with the standard definition of your CATS object

//FROM HERE ON AN EXAMPLE HOW TO GET THE EFFECTIVE GAUSSIAN, NOT NEEDED (i.e. remove in case you cannot compile it)
      const unsigned NumSourceBins = 128;
      const double rMin = 0;
      const double rMax = 16;
      TFile* fOutput = new TFile(TString::Format("%s/OtherTasks/Ramona_pK_ResoTest/fOutput.root",GetFemtoOutputFolder()),"recreate");
      TH1F* hSource = new TH1F("hSource","hSource",NumSourceBins,rMin,rMax);

      //fill the histo fro the source
      for(unsigned uBin=0; uBin<NumSourceBins; uBin++){
        //get the x-axis (r value) of the current bin
        double xaxis = hSource->GetBinCenter(uBin+1);
        //an array for the parameters, [0] is source size, [1] is == 2 (for a Gaussian)
        double parameters[2];
        parameters[0] = CoreSize;
        parameters[1] = 2.0;
        double SourceValue = MagicSource.RootEval(&xaxis, parameters);
        hSource->SetBinContent(uBin+1,SourceValue);
        //infinite errors for now
        hSource->SetBinError(uBin+1,1000.);
      }
      //idea: fit the source distribution only in a range around its peak
      //to do this: silly idea: put very large uncertainties in the bins outside of this range
      //we can get this range automatically, by evaluating the central (median) integral of the source distribution
      //with this set up, we fit the 68% (or whatever, 84% in this case) most central yield of the source distribution
      double lowerlimit;
      double upperlimit;
      //you will need to include DLM_HistoAnalysis.h
      GetCentralInterval(*hSource, 0.84, lowerlimit, upperlimit, true);
      unsigned lowerbin = hSource->FindBin(lowerlimit);
      unsigned upperbin = hSource->FindBin(upperlimit);
      for(unsigned uBin=lowerbin; uBin<=upperbin; uBin++){
        hSource->SetBinError(uBin+1,0.01);
      }

      printf("Core size of %.3f fm\n",CoreSize);
      printf("The fit will be performed in the range [%.2f, %.2f] fm\n",lowerlimit,upperlimit);
      //fyi, GaussSourceTF1 is in DLM_Source.h if you want to check it out.
      TF1* fSource = new TF1("fSource",GaussSourceScaledTF1,rMin,rMax,2);
      fSource->SetParameter(0,CoreSize);
      fSource->SetParLimits(0,CoreSize*0.5,CoreSize*2.0);
      fSource->FixParameter(1,1.);
      hSource->Fit(fSource,"S, N, R, M");
      printf("The effective Gaussian size is %.3f +/- %.3f fm\n",fSource->GetParameter(0),fSource->GetParError(0));
      printf("The lambda par is %.1f +/- %.1f %%\n",fSource->GetParameter(1)*100.,fSource->GetParError(1)*100.);

      //get rid of weird plotting
      for(unsigned uBin=0; uBin<NumSourceBins; uBin++){
        hSource->SetBinError(uBin+1,0.01);
      }
      hSource->Write();
      fSource->Write();

      delete hSource;
      delete fSource;
      delete fOutput;
}



void StableDisto_Test(){
  const double CauchySize = 1.0;
  const double GaussSize = CauchySize*2.;
  const unsigned NumBins = 256;
  const double RadMax = 16;
  TH1F* hCauchyApprox = new TH1F("hCauchyApprox","hCauchyApprox",NumBins,0,RadMax);
  TH1F* hCauchyExact = new TH1F("hCauchyExact","hCauchyExact",NumBins,0,RadMax);
  TH1F* hLevy = new TH1F("hLevy","hLevy",NumBins,0,RadMax);
  TH1F* hGauss1 = new TH1F("hGauss1","hGauss1",NumBins,0,RadMax);
  TH1F* hGauss2 = new TH1F("hGauss2","hGauss2",NumBins,0,RadMax);

  double* PARS = new double [4];
  DLM_CleverMcLevyReso Mc_Source;
  Mc_Source.InitStability(1,1-1e-6,1+1e-6);
  Mc_Source.InitScale(1,CauchySize-1e-6,CauchySize+1e-6);
  Mc_Source.InitRad(257*4,0,64);
  Mc_Source.InitType(2);
  Mc_Source.InitNumMcIter(2000000);

  for(unsigned uRad=0; uRad<NumBins; uRad++){
    double& RAD = PARS[1];
    double& SIZE = PARS[3];
    RAD = hCauchyApprox->GetBinCenter(uRad+1);
    SIZE = CauchySize;
    hLevy->SetBinContent(uRad+1,Mc_Source.Eval(PARS));
    hCauchyApprox->SetBinContent(uRad+1,CauchySource(PARS));
    SIZE *= 0.5;hCauchyExact->SetBinContent(uRad+1,CauchySource_v2(PARS));SIZE *= 2.0;
    hGauss1->SetBinContent(uRad+1,GaussSource(PARS));
    SIZE = GaussSize;
    hGauss2->SetBinContent(uRad+1,GaussSource(PARS));
  }

  TFile fOutput(TString::Format("%s/OtherTasks/StableDisto_Test.root",GetFemtoOutputFolder()),"recreate");
  hCauchyApprox->Write();
  hCauchyExact->Write();
  hLevy->Write();
  hGauss1->Write();
  hGauss2->Write();

  delete [] PARS;
  delete hCauchyApprox;
  delete hCauchyExact;
  delete hLevy;
  delete hGauss1;
  delete hGauss2;

}

void Andi_pDminus_1(){
  TString pDminus_Folder = TString::Format("%s/CatsFiles/Interaction/Haidenbauer/pDminus/",GetCernBoxDimi());
  const double kMin=0;
  const double kMax=340;
  const unsigned nMomBins = 68;
  CATSparameters cPars(CATSparameters::tSource,1,true);
  cPars.SetParameter(0,1.4);
  CATS Kitty;
  Kitty.SetMomBins(nMomBins,kMin,kMax);
  Kitty.SetAnaSource(GaussSource, cPars);
  Kitty.SetAnaSource(0, cPars.GetParameter(0));
  DLM_Histo<complex<double>>*** ExternalWF=NULL;
  ExternalWF = Init_pDminus_Haidenbauer(pDminus_Folder,Kitty,0);
  Kitty.SetExternalWaveFunction(0,0,ExternalWF[0][0][0],ExternalWF[1][0][0]);
  Kitty.SetExternalWaveFunction(1,0,ExternalWF[0][1][0],ExternalWF[1][1][0]);
  Kitty.SetChannelWeight(0,1);
  Kitty.SetChannelWeight(1,1);
  Kitty.KillTheCat();

  TGraph gCk;
  gCk.SetName("gCk");
  TGraph* gReWf_u = new TGraph[nMomBins];
  TGraph* gReWf_R = new TGraph[nMomBins];
  TGraph* gImWf_u = new TGraph[nMomBins];
  TGraph* gImWf_R = new TGraph[nMomBins];
  TGraph* gWf2 = new TGraph[nMomBins];

  TGraph* gReCc_u = new TGraph[nMomBins];
  TGraph* gReCc_R = new TGraph[nMomBins];
  TGraph* gImCc_u = new TGraph[nMomBins];
  TGraph* gImCc_R = new TGraph[nMomBins];

  const unsigned nRadBins = 256;
  const double rMin = 0;
  const double rMax = 16;
  const double rWidth = (rMax-rMin)/double(nRadBins);
  for(unsigned uMom=0; uMom<nMomBins; uMom++){
    double MOM = Kitty.GetMomentum(uMom);
    gReWf_u[uMom].SetName(TString::Format("gReWf_u_%.0f",MOM));
    gReWf_R[uMom].SetName(TString::Format("gReWf_R_%.0f",MOM));
    gImWf_u[uMom].SetName(TString::Format("gImWf_u_%.0f",MOM));
    gImWf_R[uMom].SetName(TString::Format("gImWf_R_%.0f",MOM));
    gWf2[uMom].SetName(TString::Format("gWf2_%.0f",MOM));

    gReCc_u[uMom].SetName(TString::Format("gReCc_u_%.0f",MOM));
    gReCc_R[uMom].SetName(TString::Format("gReCc_R_%.0f",MOM));
    gImCc_u[uMom].SetName(TString::Format("gImCc_u_%.0f",MOM));
    gImCc_R[uMom].SetName(TString::Format("gImCc_R_%.0f",MOM));

    gCk.SetPoint(uMom,MOM,Kitty.GetCorrFun(uMom));
    for(unsigned uRad=0; uRad<nRadBins; uRad++){
      double RAD = rWidth*0.5 + rWidth*double(uRad);
      gReWf_u[uMom].SetPoint(uRad,RAD,real(Kitty.EvalRadialWaveFunction(uMom, 0, 0, RAD, false)));
      gImWf_u[uMom].SetPoint(uRad,RAD,imag(Kitty.EvalRadialWaveFunction(uMom, 0, 0, RAD, false)));
      gReWf_R[uMom].SetPoint(uRad,RAD,real(Kitty.EvalRadialWaveFunction(uMom, 0, 0, RAD, true)));
      gImWf_R[uMom].SetPoint(uRad,RAD,imag(Kitty.EvalRadialWaveFunction(uMom, 0, 0, RAD, true)));
      gWf2[uMom].SetPoint(uRad,RAD,Kitty.EvalWaveFun2(uMom,RAD,0));

      gReCc_u[uMom].SetPoint(uRad,RAD,real(Kitty.EvalRadialWaveFunction(uMom, 1, 0, RAD, false)));
      gImCc_u[uMom].SetPoint(uRad,RAD,imag(Kitty.EvalRadialWaveFunction(uMom, 1, 0, RAD, false)));
      gReCc_R[uMom].SetPoint(uRad,RAD,real(Kitty.EvalRadialWaveFunction(uMom, 1, 0, RAD, true)));
      gImCc_R[uMom].SetPoint(uRad,RAD,imag(Kitty.EvalRadialWaveFunction(uMom, 1, 0, RAD, true)));
    }
  }

  TFile OutputFile(TString::Format("%s/OtherTasks/Andi_pDminus_1/OutputFile.root",GetFemtoOutputFolder()),"recreate");
  gCk.Write();
  for(unsigned uMom=0; uMom<nMomBins; uMom++){
    gReWf_u[uMom].Write();
    gImWf_u[uMom].Write();
    gReWf_R[uMom].Write();
    gImWf_R[uMom].Write();

    gReCc_u[uMom].Write();
    gImCc_u[uMom].Write();
    gReCc_R[uMom].Write();
    gImCc_R[uMom].Write();

    gWf2[uMom].Write();
  }

  CleanUpWfHisto(Kitty,ExternalWF);
  delete [] gReWf_u;
  delete [] gReWf_R;
  delete [] gImWf_u;
  delete [] gImWf_R;
  delete [] gReCc_u;
  delete [] gReCc_R;
  delete [] gImCc_u;
  delete [] gImCc_R;
  delete [] gWf2;
}

void Fast_Bootstrap_Example(){
  unsigned NumSystVars = 3;
  double* Constant = new double [NumSystVars];
  double* Shift = new double [NumSystVars];
  //systematic shift of 5% per bin
  Shift[0] = 0;
  Shift[1] = 0.05;
  Shift[2] = -0.05;
  Constant[0] = 1;
  Constant[1] = 1+Shift[1];
  Constant[2] = 1+Shift[2];

  //stat uncertainty of 5% per bin
  const double StatErr = 0.05;

  TRandom3 rangen(11);

  unsigned NumBins = 32;
  TH1F** hCorrFun = new TH1F* [NumSystVars];



  for(unsigned uVar=0; uVar<NumSystVars; uVar++){
    hCorrFun[uVar] = new TH1F(TString::Format("hCorrFun_%u",uVar),TString::Format("hCorrFun_%u",uVar),NumBins,0,256);
  }

  //generate experimental data, default variation
  for(unsigned uBin=0; uBin<NumBins; uBin++){
    hCorrFun[0]->SetBinContent(uBin+1,rangen.Gaus(Constant[0],StatErr));
    hCorrFun[0]->SetBinError(uBin+1,StatErr);
  }

  //generate experimental data, systematically shifted variations
  for(unsigned uVar=1; uVar<NumSystVars; uVar++){
    for(unsigned uBin=0; uBin<NumBins; uBin++){
      hCorrFun[uVar]->SetBinContent(uBin+1,hCorrFun[0]->GetBinContent(uBin+1)+Shift[uVar]);
      hCorrFun[uVar]->SetBinError(uBin+1,StatErr);
    }
  }

  TH1F* hResult_Stat = new TH1F("hResult_Stat","hResult_Stat",128,0.5,1.5);
  TH1F* hResult_Tot = new TH1F("hResult_Tot","hResult_Tot",128,0.5,1.5);
  double ExpectedStatErr;

  const unsigned NumIter = 2048;
  //perform both bootstrap (_Stat) and bootstrap+systematic (_Tot) variations
  for(unsigned uIter=0; uIter<NumIter; uIter++){
    TF1* fit_Tot = new TF1("fit_Tot","[0]",0,256);
    TF1* fit_Stat = new TF1("fit_Stat","[0]",0,256);
    fit_Tot->SetParameter(0,Constant[0]);
    fit_Stat->SetParameter(0,Constant[0]);

    TH1F* h_Tot = new TH1F("h_Tot","h_Tot",NumBins,0,256);
    TH1F* h_Stat = new TH1F("h_Stat","h_Stat",NumBins,0,256);

    for(unsigned uBin=0; uBin<NumBins; uBin++){
      //on iter 0, the default data set (no variations)
      if(uIter==0){
        h_Tot->SetBinContent(uBin+1,hCorrFun[0]->GetBinContent(uBin+1));
        h_Tot->SetBinError(uBin+1,StatErr);
        h_Stat->SetBinContent(uBin+1,hCorrFun[0]->GetBinContent(uBin+1));
        h_Stat->SetBinError(uBin+1,StatErr);
      }
      //generate a bootstrapped sample, which is than fitted to extract the value of the constant function
      else{
        unsigned WhichVar = rangen.Integer(NumSystVars);
        h_Tot->SetBinContent(uBin+1,rangen.Gaus(hCorrFun[WhichVar]->GetBinContent(uBin+1),StatErr));
        h_Tot->SetBinError(uBin+1,StatErr);
        h_Stat->SetBinContent(uBin+1,rangen.Gaus(hCorrFun[0]->GetBinContent(uBin+1),StatErr));
        h_Stat->SetBinError(uBin+1,StatErr);
      }

    }

    //fit, save the stat err on the fit parameter to compare with our result
    if(uIter==0){
      h_Tot->Fit(fit_Tot,"S, N, R, M");
      h_Stat->Fit(fit_Stat,"S, N, R, M");
      ExpectedStatErr = fit_Stat->GetParError(0);
    }
    else{
      h_Tot->Fit(fit_Tot,"Q, S, N, R, M");
      h_Stat->Fit(fit_Stat,"Q, S, N, R, M");
    }

    //save the distribution of fit results (the constant value of our input function)
    hResult_Tot->Fill(fit_Tot->GetParameter(0));
    hResult_Stat->Fill(fit_Stat->GetParameter(0));

    delete fit_Tot;
    delete fit_Stat;
    delete h_Tot;
    delete h_Stat;
  }

  printf(" ExpectedStatErr = %.2e\n",ExpectedStatErr);
  printf("ExtractedStatErr = %.2e\n",hResult_Stat->GetStdDev());
  printf(" ExtractedTotErr = %.2e\n",hResult_Tot->GetStdDev());
  printf("ExtractedSystErr = %.2e\n",sqrt(hResult_Tot->GetStdDev()*hResult_Tot->GetStdDev()-hResult_Stat->GetStdDev()*hResult_Stat->GetStdDev()));


  TFile OutputFile(TString::Format("%s/OtherTasks/Fast_Bootstrap_Example/fOutput.root",GetFemtoOutputFolder()), "recreate");
  for(unsigned uVar=0; uVar<NumSystVars; uVar++){
    hCorrFun[uVar]->Write();
  }
  hResult_Tot->Write();
  hResult_Stat->Write();
  //delete hCorrFun;
  //delete [] Constant;
}

bool Eval_ScattParameters(CATS& Kitty, double& ScatLen, double& EffRan, TH1F*& hFit, TF1*& fitSP,
  const int& Nterms, const bool& Fixf0, const bool& Fixd0, const unsigned short& usCh){
  Kitty.KillTheCat();
  double* MomBins = Kitty.CopyMomBin();
  hFit = new TH1F("hFit","hFit",Kitty.GetNumMomBins(),MomBins);
  double LAST_POINT;
  double CURRENT_POINT;
  for(unsigned uMom=0; uMom<Kitty.GetNumMomBins(); uMom++){
    CURRENT_POINT = Kitty.GetMomentum(uMom)/tan(Kitty.GetPhaseShift(uMom,usCh,0));
    if(uMom){
      if(CURRENT_POINT*LAST_POINT<0&&fabs(CURRENT_POINT-LAST_POINT)>1000&&Kitty.GetMomentum(uMom)<120)
      {fitSP=NULL;delete[]MomBins;return false;}
    }
    hFit->SetBinContent(uMom+1,CURRENT_POINT);
    hFit->SetBinError(uMom+1,1.);
    LAST_POINT = CURRENT_POINT;
  }
  TF1* fitSP2 = new TF1("fitSP2","0.5*[1]/197.327*x*x+197.327*[0]", 10, 90);
  TF1* fitSP4 = new TF1("fitSP4","[2]*pow(x,4.)+0.5*[1]/197.327*x*x+197.327*[0]", 10, 90);
  TF1* fitSP6 = new TF1("fitSP6","[3]*pow(x,6.)+[2]*pow(x,4.)+0.5*[1]/197.327*x*x+197.327*[0]", 10, 90);
  double inv_f0 = ScatLen==0?0:1./ScatLen;
  //EffRan = 0;
//printf("ScatLen = %e\n",ScatLen);
//printf("inv_f0 = %f\n",inv_f0);
//printf("EffRan = %f\n",EffRan);
  if(Fixf0) {fitSP2->FixParameter(0,inv_f0);fitSP4->FixParameter(0,inv_f0);fitSP6->FixParameter(0,inv_f0);}
  else {fitSP2->SetParameter(0,inv_f0);fitSP2->SetParLimits(0,-100,100);
        fitSP4->SetParameter(0,inv_f0);fitSP4->SetParLimits(0,-100,100);
        fitSP6->SetParameter(0,inv_f0);fitSP6->SetParLimits(0,-100,100);}
  if(Fixd0) { fitSP2->FixParameter(1,EffRan);
              fitSP4->FixParameter(1,EffRan);
              fitSP6->FixParameter(1,EffRan);}
  else {fitSP2->SetParameter(1,EffRan);fitSP2->SetParLimits(1,-50,50);
        fitSP4->SetParameter(1,EffRan);fitSP4->SetParLimits(1,-50,50);
        fitSP6->SetParameter(1,EffRan);fitSP6->SetParLimits(1,-50,50);}
  fitSP4->SetParameter(2,0);fitSP6->SetParameter(2,0);
  fitSP6->SetParameter(3,0);


  double Chi2_Old = 1e64;

  //hFit->Fit(fitSP2, "Q, S, N, R, M");
  ROOT::Math::MinimizerOptions MinOpt;
  MinOpt.SetMinimizerType("Minuit2");
  MinOpt.SetPrintLevel(0);
  DLM_FitHisto(hFit, fitSP2, "Q, S, N, R, M", "", &MinOpt);
  //printf("f0 %f\n", 1./fitSP2->GetParameter(0));
  ScatLen = 1./fitSP2->GetParameter(0);
  EffRan = fitSP2->GetParameter(1);
  if(Nterms<=2){delete fitSP4; delete fitSP6; fitSP=fitSP2; delete[]MomBins; return true;}

  hFit->Fit(fitSP4, "S, N, R, M");
  DLM_FitHisto(hFit, fitSP4, "Q, S, N, R, M", "", &MinOpt);
  if(fitSP4->GetChisquare()/fitSP2->GetChisquare()>0.8)
  {delete fitSP4; delete fitSP6; fitSP=fitSP2; delete[]MomBins; return true;}
  ScatLen = 1./fitSP4->GetParameter(0);
  EffRan = fitSP4->GetParameter(1);
  if(Nterms<=3){delete fitSP2; delete fitSP6; fitSP=fitSP4; delete[]MomBins; return true;}

  //hFit->Fit(fitSP6, "Q, S, N, R, M");
  DLM_FitHisto(hFit, fitSP6, "Q, S, N, R, M", "", &MinOpt);
  if(fitSP6->GetChisquare()/fitSP4->GetChisquare()>0.8)
  {delete fitSP2; delete fitSP6; fitSP=fitSP4; delete[]MomBins; return true;}
  ScatLen = 1./fitSP6->GetParameter(0);
  EffRan = fitSP6->GetParameter(1);
  delete fitSP2; delete fitSP4; fitSP=fitSP6; delete[]MomBins; return true;
}

void Ledni_SmallRad(TString PotentialName){
  //const TString PotentialName = "NSC97b";
  //const TString PotentialName = "NF48";
  //const TString PotentialName = "emma";
  //const TString PotentialName = "custom";

  const double kMin = 0;
  const double kMax = 300;
  const double kStep = 3;
  const unsigned nMom = TMath::Nint(kMax/kStep);
  const double Radius = 1.2;

  CATSparameters sPars(CATSparameters::tSource,1,true);
  sPars.SetParameter(0,Radius);
  CATS Kitty_SE;
  Kitty_SE.SetMomBins(nMom,kMin,kMax);
  Kitty_SE.SetAnaSource(GaussSource, sPars);
  Kitty_SE.SetUseAnalyticSource(true);
  Kitty_SE.SetQ1Q2(0);
  Kitty_SE.SetQuantumStatistics(false);
  //Kitty_SE.SetRedMass(Mass_L*0.5);
  Kitty_SE.SetRedMass((Mass_p*Mass_Kch)/(Mass_p+Mass_Kch));
  Kitty_SE.SetNumChannels(1);
  Kitty_SE.SetNumPW(0,1);
  Kitty_SE.SetSpin(0,0);
  Kitty_SE.SetChannelWeight(0, 1.);
  CATSparameters pPars(CATSparameters::tPotential,7,true);
  double c_f0,c_d0;
  if(PotentialName=="NSC97b"){
    pPars.SetParameter(0,-78.42);
    pPars.SetParameter(1,1.0);
    pPars.SetParameter(2,741.76);
    pPars.SetParameter(3,0.45);
    c_f0 = 0.397;
    c_d0 = 10.360;
  }
  else if(PotentialName=="NF48"){
    pPars.SetParameter(0,-1647.40);
    pPars.SetParameter(1,0.6);
    pPars.SetParameter(2,3888.96);
    pPars.SetParameter(3,0.45);
    c_f0 = 1.511;
    c_d0 = 2.549;
  }
  else if(PotentialName=="emma"){
    //for emma -> NSC97f basis
    pPars.SetParameter(0,-106.53*0.85);
    pPars.SetParameter(1,1.0*1.18);
    pPars.SetParameter(2,1469.33);
    pPars.SetParameter(3,0.45*1.1);
    c_f0 = 0.350;
    c_d0 = 16.330;
  }
  else if(PotentialName=="Toy1"){
    pPars.SetParameter(0,-144.5);
    pPars.SetParameter(1,2.11);
    pPars.SetParameter(2,520.0);
    pPars.SetParameter(3,0.54);
    c_f0 = -0.73;
    c_d0 = 7.72;
  }
  else if(PotentialName=="ND46"){
    pPars.SetParameter(0,144.89);
    pPars.SetParameter(1,1.0);
    pPars.SetParameter(2,127.87);
    pPars.SetParameter(3,0.45);
    c_f0 = -4.621;
    c_d0 = 1.3;
  }
  else if(PotentialName=="Yukawa1"){
    pPars.SetParameter(0,1.0);
    pPars.SetParameter(1,1.0);
    pPars.SetParameter(2,100.0);
    pPars.SetParameter(3,1.0);
    pPars.SetParameter(4,0.4);
    c_f0 = 0;
    c_d0 = 1;
  }
  else if(PotentialName=="pKplusI0"){
    pPars.SetParameter(0,0.0);
    c_f0 = 0.03;
    c_d0 = 0.0;
  }
  else if(PotentialName=="pKplusI1"){
    pPars.SetParameter(0,1.0);
    c_f0 = -0.3;
    c_d0 = 0.0;
  }
  else if(PotentialName=="pKplusYuki"){
    pPars.SetParameter(0,0.376);//0.376;0.335
    pPars.SetParameter(1,sqrt(200.*(Mass_p*Mass_Kch)/(Mass_p+Mass_Kch)));
    pPars.SetParameter(2,3);
    pPars.SetParameter(3,2084);
    pPars.SetParameter(4,50.81);
    pPars.SetParameter(5,18.34);
    pPars.SetParameter(6,-1.752);
    c_f0 = -0.3;
    c_d0 = 0.0;
  }
  else{
    pPars.SetParameter(0,-5.50337);
    pPars.SetParameter(1,1./sqrt(2.148805));
    pPars.SetParameter(2,0);
    pPars.SetParameter(3,1);
    //pPars.SetParameter(0,-78.42*0.39*4.5);//0.39,0.4
    //pPars.SetParameter(1,1.0*1.35);
    //pPars.SetParameter(2,741.76*4.5);
    //pPars.SetParameter(3,0.45*1.4);
    //NF46 as a stariting point
    //pPars.SetParameter(0,-1327.26*1.0);
    //pPars.SetParameter(1,0.6);
    //pPars.SetParameter(2,2561.56);
    //pPars.SetParameter(3,0.45);
    c_f0 = 0.02;
    c_d0 = 30.0;
    printf(" Hello\n");

  }
  Kitty_SE.SetEpsilonConv(5e-9);
  Kitty_SE.SetEpsilonProp(5e-9);
  if(PotentialName.Contains("Yukawa"))
    Kitty_SE.SetShortRangePotential(0,0,YukawaDimiCore,pPars);
  else if(PotentialName.Contains("pKplusI")){
    Kitty_SE.SetShortRangePotential(0,0,KpProtonEquivalentPotential,pPars);
    Kitty_SE.SetEpsilonConv(1e-9);
    Kitty_SE.SetEpsilonProp(1e-9);
  }
  else if(PotentialName.Contains("pKplusYuki")){
    Kitty_SE.SetShortRangePotential(0,0,SingleGaussDynamic,pPars);
    Kitty_SE.SetEpsilonConv(1e-9);
    Kitty_SE.SetEpsilonProp(1e-9);
  }
  else Kitty_SE.SetShortRangePotential(0,0,DoubleGaussSum,pPars);
  Kitty_SE.KillTheCat();

  const TString OutputFolder = TString::Format("%s/MMM/Ledni_SmallRad/",GetCernBoxDimi());
  TFile fOutput(OutputFolder+"fOutput_"+PotentialName+".root","recreate");

  TH1F* h_kcotd=NULL;
  TF1* f_kcotd=NULL;
  if(PotentialName.Contains("pKplus")) Eval_ScattParameters(Kitty_SE,c_f0,c_d0,h_kcotd,f_kcotd,2,false,false);
  else Eval_ScattParameters(Kitty_SE,c_f0,c_d0,h_kcotd,f_kcotd);
  //c_f0 = 1./c_f0;
  printf("c_f0 = %.2f fm\n",c_f0);
  printf("c_d0 = %.2f fm\n",c_d0);
  //printf("h_kcotd = %p\n",h_kcotd);
  //printf("f_kcotd = %p\n",f_kcotd);
  fOutput.cd();
  h_kcotd->Write();
  f_kcotd->Write();
  delete h_kcotd;
  delete f_kcotd;
printf("-->\n");
  const double Re_f0 = PotentialName=="NSC97b"?0.397:PotentialName=="NF48"?1.511:c_f0;//NSC97b
  const double Im_f0 = 0;
  double Val_d0 = PotentialName=="NSC97b"?10.36:PotentialName=="NF48"?2.549:c_d0;
  const complex<double> Val_f0(Re_f0*FmToNu,Im_f0*FmToNu);
  const complex<double> Inv_f0 = 1./(Val_f0+1e-64);
  const double Val_r0 = Radius*FmToNu;
  Val_d0 *= FmToNu;

  TGraph gCk_SE;
  gCk_SE.SetName("gCk_"+PotentialName);
  gCk_SE.SetLineWidth(3);
  gCk_SE.SetLineColor(kGreen);
  TGraph gCk_default;
  gCk_default.SetName("gCk_default");
  gCk_default.SetLineWidth(3);
  gCk_default.SetLineColor(kBlack);
  TGraph gCk_csmall;
  gCk_csmall.SetName("gCk_csmall");
  gCk_csmall.SetLineWidth(3);
  gCk_csmall.SetLineColor(kRed+2);
  TGraph g_csmall;
  g_csmall.SetName("g_csmall");
  g_csmall.SetLineWidth(2);
  g_csmall.SetLineColor(kRed);
  TGraph g_f2;
  g_f2.SetName("g_f2");
  g_f2.SetLineWidth(2);
  g_f2.SetLineColor(kGreen+1);
  TGraph g_f2c;
  g_f2c.SetName("g_f2c");
  g_f2c.SetLineWidth(2);
  g_f2c.SetLineColor(kYellow+2);
  TGraph g_F1F2;
  g_F1F2.SetName("g_F1F2");
  g_F1F2.SetLineWidth(2);
  g_F1F2.SetLineColor(kBlue);

  TGraph g_F1F2_f2c;
  g_F1F2_f2c.SetName("g_F1F2_f2c");
  g_F1F2_f2c.SetLineWidth(2);
  g_F1F2_f2c.SetLineColor(kGray+1);

  TGraph g_f0d0;
  g_f0d0.SetName("g_f0d0");
  g_f0d0.SetMarkerStyle(2);
  g_f0d0.SetMarkerSize(2);
  g_f0d0.SetPoint(0,c_f0,c_d0);

  TGraph g_Vmu;
  g_Vmu.SetName("g_Vmu");
  g_Vmu.SetMarkerStyle(3);
  g_Vmu.SetMarkerSize(2);
  g_Vmu.SetPoint(0,pPars.GetParameter(0),pPars.GetParameter(1));
  g_Vmu.SetPoint(1,pPars.GetParameter(2),pPars.GetParameter(3));

  //the true correction
  TGraph g_XXX;
  g_XXX.SetName("g_XXX");
  g_XXX.SetLineWidth(4);
  g_XXX.SetLineStyle(2);
  g_XXX.SetLineColor(kPink+1);

  //attempted correction with MM model
  //TGraph g_MM;
  //g_MM.SetName("g_MM");
  //g_MM.SetLineWidth(4);
  //g_MM.SetLineColor(kPink+3);

  TGraph g_XXXa;
  g_XXXa.SetName("g_XXXa");
  g_XXXa.SetLineWidth(4);
  g_XXXa.SetLineStyle(3);
  g_XXXa.SetLineColor(kPink+1);

  TGraph g_XXXm;
  g_XXXm.SetName("g_XXXm");
  g_XXXm.SetLineWidth(4);
  g_XXXm.SetLineStyle(4);
  g_XXXm.SetLineColor(kPink+1);

  TGraph g_delta;
  g_delta.SetName("g_delta");
  g_delta.SetLineWidth(4);
  g_delta.SetLineStyle(2);
  g_delta.SetLineColor(kBlack);

  TGraph g_kcotdelta;
  g_kcotdelta.SetName("g_kcotdelta");
  g_kcotdelta.SetLineWidth(4);
  g_kcotdelta.SetLineStyle(2);
  g_kcotdelta.SetLineColor(kBlack);

  TGraph g_cotdelta;
  g_cotdelta.SetName("g_cotdelta");
  g_cotdelta.SetLineWidth(4);
  g_cotdelta.SetLineStyle(2);
  g_cotdelta.SetLineColor(kBlack);

  TGraph g_tgdelta;
  g_tgdelta.SetName("g_tgdelta");
  g_tgdelta.SetLineWidth(4);
  g_tgdelta.SetLineStyle(2);
  g_tgdelta.SetLineColor(kBlack);

  TGraph g_sindelta;
  g_sindelta.SetName("g_sindelta");
  g_sindelta.SetLineWidth(4);
  g_sindelta.SetLineStyle(2);
  g_sindelta.SetLineColor(kBlack);

  //real part of scatt. amplitude as a function of sqrt(s)
  TGraph g_s_ReF;
  g_s_ReF.SetName("g_s_ReF");
  g_s_ReF.SetLineWidth(4);
  g_s_ReF.SetLineStyle(2);
  g_s_ReF.SetLineColor(kBlack);

  //imag part of scatt. amplitude as a function of sqrt(s)
  TGraph g_s_ImF;
  g_s_ImF.SetName("g_s_ImF");
  g_s_ImF.SetLineWidth(4);
  g_s_ImF.SetLineStyle(2);
  g_s_ImF.SetLineColor(kBlack);

  TGraph g_MM_B;
  g_MM_B.SetName("g_MM_B");
  g_MM_B.SetLineWidth(4);
  g_MM_B.SetLineColor(kAzure+10);

  TGraph g_MM_A;
  g_MM_A.SetName("g_MM_A");
  g_MM_A.SetLineWidth(4);
  g_MM_A.SetLineColor(kViolet+1);

  TGraph g_MM_XXX;
  g_MM_XXX.SetName("g_MM_XXX");
  g_MM_XXX.SetLineWidth(4);
  g_MM_XXX.SetLineColor(kRed);


  unsigned NumPts = 0;
  double XXX0;
//I am 100% sure it goes as 1/Radius
  double exp_MMB = c_d0/pow(TMath::Pi(),3.)/Radius;
  double exp_MMA;
printf("--->\n");
  for(double MOM=kMin+kStep*0.5; MOM<kMax; MOM+=kStep){
    double F1 = gsl_sf_dawson(2.*MOM*Val_r0)/(2.*MOM*Val_r0);
    double F2 = (1.-exp(-4.*MOM*MOM*Val_r0*Val_r0))/(2.*MOM*Val_r0);
    complex<double> Val_f = pow(Inv_f0+0.5*Val_d0*MOM*MOM-i*MOM,-1.);
    double csmall = (1-(Val_d0)/(2*sqrt(Pi)*Val_r0));
    double v_f2 = 0.5*pow(abs(Val_f)/Val_r0,2);
    double F1F2 = 2*real(Val_f)*F1/(sqrt(Pi)*Val_r0)-imag(Val_f)*F2/Val_r0;

    double Ck_default = 1.+v_f2+F1F2;
    double Ck_csmall = 1.+v_f2*csmall+F1F2;

    gCk_default.SetPoint(NumPts,MOM,Ck_default);
    gCk_csmall.SetPoint(NumPts,MOM,Ck_csmall);
    g_csmall.SetPoint(NumPts,MOM,csmall);
    g_f2.SetPoint(NumPts,MOM,v_f2);
    g_f2c.SetPoint(NumPts,MOM,v_f2*csmall);
    g_F1F2.SetPoint(NumPts,MOM,F1F2);
    g_F1F2_f2c.SetPoint(NumPts,MOM,fabs(F1F2/v_f2/csmall));

    gCk_SE.SetPoint(NumPts,MOM,Kitty_SE.GetCorrFun(NumPts));

    //the correction as in ledni
    double XXX = (Kitty_SE.GetCorrFun(NumPts)-1.-F1F2)/v_f2;
    if(NumPts==0){
      XXX0 = XXX;
      exp_MMA = XXX0-exp_MMB;
    }
    g_XXX.SetPoint(NumPts,MOM,XXX);
    //the correction concidered as a general addition factor to C(k)
    g_XXXa.SetPoint(NumPts,MOM,(Kitty_SE.GetCorrFun(NumPts)-1.-v_f2-F1F2));
    //the correction concidered as a general multiplication factor to C(k)
    g_XXXm.SetPoint(NumPts,MOM,Kitty_SE.GetCorrFun(NumPts)/(1.+v_f2+F1F2));

    //g_MM_XXX.SetPoint(NumPts,MOM,exp_MMA-exp_MMB*sqrt(1.+MOM*MOM*0.5*c_f0*c_d0*FmToNu*FmToNu));

    //wrong still
    double UnderRoot = 1.+0.5*MOM*MOM*c_d0*c_f0*FmToNu*FmToNu;
    double RootSign = -1;
    //bond state
    if(c_f0<0&&(1.+2*c_d0/c_f0)>0) RootSign = 1;
    if(UnderRoot<0) {UnderRoot=-UnderRoot; RootSign=-RootSign;}
    double MM_XXX;
    MM_XXX = 1.+RootSign*c_d0/(2.*sqrt(Pi)*Radius)*sqrt(UnderRoot);
    g_MM_XXX.SetPoint(NumPts,MOM,MM_XXX);

    g_delta.SetPoint(NumPts,MOM,Kitty_SE.GetPhaseShift(NumPts,0,0));
    g_kcotdelta.SetPoint(NumPts,MOM,MOM/tan(Kitty_SE.GetPhaseShift(NumPts,0,0)));
    g_cotdelta.SetPoint(NumPts,MOM,1./tan(Kitty_SE.GetPhaseShift(NumPts,0,0)));
    g_tgdelta.SetPoint(NumPts,MOM,tan(Kitty_SE.GetPhaseShift(NumPts,0,0)));
    g_sindelta.SetPoint(NumPts,MOM,sin(Kitty_SE.GetPhaseShift(NumPts,0,0)));

    complex<double> ScatAmpl = exp(i*Kitty_SE.GetPhaseShift(NumPts,0,0))*sin(Kitty_SE.GetPhaseShift(NumPts,0,0))/(MOM);
    g_s_ReF.SetPoint(NumPts,sqrt(Mass_Kch*Mass_Kch+MOM*MOM)+sqrt(Mass_p*Mass_p+MOM*MOM),real(ScatAmpl)*197.);
    g_s_ImF.SetPoint(NumPts,sqrt(Mass_Kch*Mass_Kch+MOM*MOM)+sqrt(Mass_p*Mass_p+MOM*MOM),imag(ScatAmpl)*197.);

    double MM_B,MM_A;
    //MM_B = (XXX-XXX0)/(sqrt(fabs(c_f0*(1./c_f0+0.5*MOM*MOM*c_d0)))-1.);
    //MM_B = (XXX-XXX0)/(c_f0*(1./c_f0+0.5*MOM*MOM*c_d0)-1.);
    MM_B = (XXX-XXX0)/(sqrt(fabs(1.+MOM*MOM*0.5*c_f0*c_d0*FmToNu*FmToNu))-1.);
    MM_A = XXX0-MM_B;
    //printf("k = %.1f\n",MOM);
    //printf(" X=%.3f\n",XXX);
    //printf(" X-X0=%.3f\n",XXX-XXX0);
    //printf(" C = %.2e\n",0.5*c_f0*c_d0*FmToNu*FmToNu);
    //printf(" prop=%.3f\n",sqrt(fabs(1.+MOM*MOM*0.5*c_f0*c_d0*FmToNu*FmToNu))-1.);
    //printf(" B = %.2e\n",MM_B);
    //printf(" A = %.2e\n",MM_A);
    //printf("----------------------\n");
    //usleep(500e3);
    if(NumPts){
      g_MM_B.SetPoint(NumPts-1,MOM,MM_B);
      g_MM_A.SetPoint(NumPts-1,MOM,MM_A);
    }

    NumPts++;
  }

  TGraph g_pot0;
  g_pot0.SetName("g_pot0");
  g_pot0.SetLineWidth(4);
  g_pot0.SetLineStyle(2);
  g_pot0.SetLineColor(kBlue-1);
  unsigned uRad=0;
  for(double RAD=0.01; RAD<5; RAD+=0.01){
    g_pot0.SetPoint(uRad++,RAD,Kitty_SE.EvaluateThePotential(0,0,0,RAD));
  }
  TGraph g_pot360;
  g_pot360.SetName("g_pot360");
  g_pot360.SetLineWidth(4);
  g_pot360.SetLineStyle(2);
  g_pot360.SetLineColor(kGreen-1);
  uRad=0;
  for(double RAD=0.01; RAD<5; RAD+=0.01){
    g_pot360.SetPoint(uRad++,RAD,Kitty_SE.EvaluateThePotential(0,0,359.7,RAD));
  }


  TGraph g_yukipot0;
  g_yukipot0.SetName("g_yukipot0");
  g_yukipot0.SetLineWidth(2);
  g_yukipot0.SetLineStyle(1);
  g_yukipot0.SetLineColor(kBlue+1);

  TGraph g_yukipot360;
  g_yukipot360.SetName("g_yukipot360");
  g_yukipot360.SetLineWidth(2);
  g_yukipot360.SetLineStyle(1);
  g_yukipot360.SetLineColor(kGreen+1);


  //for pKplus - read yuki file
  float fRad;
  float fPot;
  char* InputFileName = new char[256];
  char* cdummy = new char[256];
  strcpy(InputFileName,GetCernBoxDimi());
  strcat(InputFileName,"/CatsFiles_Dimi/pKplus/YukiPotentials/pote_E_0MeV.dat");
  FILE *InFile;
  InFile = fopen(InputFileName, "r");
  if(!InFile){
      printf("          \033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", InputFileName);
      return;
  }
  fseek ( InFile , 0 , SEEK_END );
  long EndPos;
  EndPos = ftell (InFile);
  fseek ( InFile , 0 , SEEK_SET );
  long CurPos;
  //Read the header lines
  for(unsigned short us=0; us<1; us++){
      if(!fgets(cdummy, 255, InFile)){
          printf("\033[1;33mWARNING!\033[0m Possible bad input-file, error when reading from %s!\n",InputFileName);
          continue;
      }
  }
  uRad=0;
  while(!feof(InFile)){
    if(!fscanf(InFile,"%f %f",&fRad,&fPot)){
        printf("\033[1;33mWARNING!\033[0m Possible bad input-file, error when reading from %s!\n",InputFileName);
        continue;
    }
    g_yukipot0.SetPoint(uRad++,fRad,fPot);
  }
  fclose(InFile);
  //delete InFile;

  strcpy(InputFileName,GetCernBoxDimi());
  strcat(InputFileName,"/CatsFiles_Dimi/pKplus/YukiPotentials/pote_E_200MeV.dat");
  InFile = fopen(InputFileName, "r");
  if(!InFile){
      printf("          \033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", InputFileName);
      return;
  }
  fseek ( InFile , 0 , SEEK_END );
  EndPos = ftell (InFile);
  fseek ( InFile , 0 , SEEK_SET );
  //Read the header lines
  for(unsigned short us=0; us<1; us++){
      if(!fgets(cdummy, 255, InFile)){
          printf("\033[1;33mWARNING!\033[0m Possible bad input-file, error when reading from %s!\n",InputFileName);
          continue;
      }
  }
  uRad=0;
  while(!feof(InFile)){
    if(!fscanf(InFile,"%f %f",&fRad,&fPot)){
        printf("\033[1;33mWARNING!\033[0m Possible bad input-file, error when reading from %s!\n",InputFileName);
        continue;
    }
    g_yukipot360.SetPoint(uRad++,fRad,fPot);
  }
  fclose(InFile);

  //printf("f0 = %.2f\n",c_f0);
  //printf("1/f0 = %.2f\n",1./c_f0);
  //printf("MMA = %.2f\n",exp_MMA);
  //printf("d0 = %.2f\n",c_d0);
  //printf("MMB = %.2f\n",exp_MMB);


  gCk_SE.Write();
  gCk_default.Write();
  gCk_csmall.Write();
  g_csmall.Write();
  g_f2.Write();
  g_f2c.Write();
  g_F1F2.Write();
  g_F1F2_f2c.Write();
  g_XXX.Write();
  g_MM_XXX.Write();
  g_MM_A.Write();
  g_MM_B.Write();
  g_XXXa.Write();
  g_XXXm.Write();
  g_delta.Write();
  g_kcotdelta.Write();
  g_cotdelta.Write();
  g_tgdelta.Write();
  g_sindelta.Write();
  g_s_ReF.Write();
  g_s_ImF.Write();
  g_f0d0.Write();
  g_Vmu.Write();
  g_pot0.Write();
  g_pot360.Write();
  g_yukipot0.Write();
  g_yukipot360.Write();
  //h_kcotd->Write();
  //f_kcotd->Write();

  //delete h_kcotd;
  //delete f_kcotd;
  delete [] InputFileName;
  delete [] cdummy;
}

void Write_MM_Plots(TFile*fOutput, const CATS& Kitty, const double& r0, const double& f0, const double& d0, const TString& suffix){
    fOutput->cd();
    TGraph gCk_SE;
    gCk_SE.SetName("gCk"+suffix);
    gCk_SE.SetLineWidth(3);
    gCk_SE.SetLineColor(kGreen);
    TGraph gCk_LL;
    gCk_LL.SetName("gCk_LL"+suffix);
    gCk_LL.SetLineWidth(3);
    gCk_LL.SetLineColor(kBlack);
    TGraph gCk_LLX;
    gCk_LLX.SetName("gCk_LLX"+suffix);
    gCk_LLX.SetLineWidth(3);
    gCk_LLX.SetLineColor(kRed+2);
    TGraph g_LLX;
    g_LLX.SetName("g_LLX"+suffix);
    g_LLX.SetLineWidth(2);
    g_LLX.SetLineColor(kRed);
    TGraph g_f2;
    g_f2.SetName("g_f2"+suffix);
    g_f2.SetLineWidth(2);
    g_f2.SetLineColor(kGreen+1);
    TGraph g_f2X;
    g_f2X.SetName("g_f2X"+suffix);
    g_f2X.SetLineWidth(2);
    g_f2X.SetLineColor(kYellow+2);
    TGraph g_F1F2;
    g_F1F2.SetName("g_F1F2"+suffix);
    g_F1F2.SetLineWidth(2);
    g_F1F2.SetLineColor(kBlue);

    TGraph g_F1F2_f2X;
    g_F1F2_f2X.SetName("g_F1F2_f2X"+suffix);
    g_F1F2_f2X.SetLineWidth(2);
    g_F1F2_f2X.SetLineColor(kGray+1);

    TGraph g_PS;
    g_PS.SetName("g_PS"+suffix);
    g_PS.SetLineWidth(4);
    g_PS.SetLineStyle(1);
    g_PS.SetLineColor(kPink+2);

    TGraph g_kcot;
    g_kcot.SetName("g_kcot"+suffix);
    g_kcot.SetLineWidth(4);
    g_kcot.SetLineStyle(2);
    g_kcot.SetLineColor(kPink+2);

    TGraph g_sqrtkcot;
    g_sqrtkcot.SetName("g_sqrtkcot"+suffix);
    g_sqrtkcot.SetLineWidth(4);
    g_sqrtkcot.SetLineStyle(3);
    g_sqrtkcot.SetLineColor(kPink+2);


    //the true correction
    TGraph g_MMX;
    g_MMX.SetName("g_MMX"+suffix);
    g_MMX.SetLineWidth(4);
    g_MMX.SetLineStyle(2);
    g_MMX.SetLineColor(kPink+1);

    TGraph g_MMB;
    g_MMB.SetName("g_MMB"+suffix);
    g_MMB.SetLineWidth(4);
    g_MMB.SetLineColor(kAzure+10);

    TGraph g_MMA;
    g_MMA.SetName("g_MMA"+suffix);
    g_MMA.SetLineWidth(4);
    g_MMA.SetLineColor(kViolet+1);

    unsigned NumPts = 0;
    double if0 = 1./f0;;
    for(unsigned uMom=0; uMom<Kitty.GetNumMomBins(); uMom++){
      double MOM = Kitty.GetMomentum(uMom);
      double F1 = gsl_sf_dawson(2.*MOM*r0*FmToNu)/(2.*MOM*r0*FmToNu);
      double F2 = (1.-exp(-4.*MOM*MOM*r0*r0*FmToNu*FmToNu))/(2.*MOM*r0*FmToNu);
      complex<double> Val_f = pow(if0/FmToNu+0.5*d0*FmToNu*MOM*MOM-i*MOM,-1.);
      double LLX = (1-(d0)/(2*sqrt(Pi)*r0));
      double v_f2 = 0.5*pow(abs(Val_f)/r0/FmToNu,2);
      double F1F2 = 2*real(Val_f)*F1/(sqrt(Pi)*r0*FmToNu)-imag(Val_f)*F2/r0/FmToNu;

      double Ck_LL = 1.+v_f2+F1F2;
      double Ck_LLX = 1.+v_f2*LLX+F1F2;

      gCk_LL.SetPoint(NumPts,MOM,Ck_LL);
      gCk_LLX.SetPoint(NumPts,MOM,Ck_LLX);
      g_LLX.SetPoint(NumPts,MOM,LLX);
      g_f2.SetPoint(NumPts,MOM,v_f2);
      g_f2X.SetPoint(NumPts,MOM,v_f2*LLX);
      g_F1F2.SetPoint(NumPts,MOM,F1F2);
      g_F1F2_f2X.SetPoint(NumPts,MOM,fabs(F1F2/v_f2/LLX));

      gCk_SE.SetPoint(NumPts,MOM,Kitty.GetCorrFun(NumPts));
      g_PS.SetPoint(NumPts,MOM,Kitty.GetPhaseShift(NumPts,0,0));
      g_kcot.SetPoint(NumPts,MOM,MOM/tan(Kitty.GetPhaseShift(NumPts,0,0)));
      g_sqrtkcot.SetPoint(NumPts,MOM,sqrt(fabs(MOM/tan(Kitty.GetPhaseShift(NumPts,0,0)))));

      //the correction as in ledni
      double MMX = (Kitty.GetCorrFun(NumPts)-1.-F1F2)/v_f2;
      g_MMX.SetPoint(NumPts,MOM,MMX);

      NumPts++;
    }


    TGraph g_pot;
    g_pot.SetName("g_pot");
    g_pot.SetLineWidth(4);
    g_pot.SetLineStyle(2);
    g_pot.SetLineColor(kBlue-1);
    unsigned uRad=0;
    for(double RAD=0.1; RAD<5; RAD+=0.025){
      g_pot.SetPoint(uRad++,RAD,Kitty.EvaluateThePotential(0,0,50,RAD));
    }

    gCk_SE.Write();
    gCk_LL.Write();
    gCk_LLX.Write();
    g_LLX.Write();
    g_f2.Write();
    g_f2X.Write();
    g_F1F2.Write();
    g_F1F2_f2X.Write();
    g_PS.Write();
    g_kcot.Write();
    g_sqrtkcot.Write();
    g_MMX.Write();
    g_pot.Write();
}

int Get_MMclass(const double& r0, const double& f0, const double& d0){
  //f0>0
  //r0>d0
  //|f0|>2*|d0| (BS cond) -> satisfied for 2 dominant or 3 surpressed
  //r0>|f|/2
  //interesting classes:
  //+/- AB, with A,B = 1,2,3 (r0,f0,d0)
  //A = which one is dominant (e.g. 1 if r0>d0 && r0>|f|/2)
  //B = which one is surpressed (e.g. 3 if d0<r0 && d0<|f0|/2)
  int CLASS = 0;
  if(r0<fabs(d0)&&r0<fabs(f0)*0.5) CLASS = 1;
  else if(fabs(f0)<2.*r0&&fabs(f0)<2.*fabs(d0)) CLASS = 2;
  else if(fabs(d0)<r0&&fabs(d0)<fabs(f0)*0.5) CLASS = 3;

  if(r0>fabs(d0)&&r0>fabs(f0)*0.5) CLASS += 10;
  else if(fabs(f0)>2.*r0&&fabs(f0)>2.*fabs(d0)) CLASS += 20;
  else if(fabs(d0)>r0&&fabs(d0)>fabs(f0)*0.5) CLASS += 30;
  CLASS = f0>0?CLASS:-CLASS;
  return CLASS;
}

//DimiFlag -> any further info. At the moment:
//0 - randomly generated potential
//10 - ManufacturePotential used (random phase)
//11 - ManufacturePotential used (converged phase)
void Fill_ntMM(TNtuple* ntMM, const CATS& Kitty, const double& r0, const double& f0, const double& d0,
                const double& V1, const double& mu1, const double& V2, const double& mu2, const double& s2, const int& DimiFlag){
    //ROOT::Math::MinimizerOptions MinOpt;
    //MinOpt.SetMinimizerType("Minuit2");
    //MinOpt.SetPrintLevel(0);

    const unsigned nMom = Kitty.GetNumMomBins();
    double* MomBins = Kitty.CopyMomBin();
    const double kCorrection = 120;
    TH1F* hMMX = new TH1F("hMMX","hMMX",nMom,MomBins);
    TF1* fMMX = new TF1("fMMX","[0]+[1]*sqrt(1.+[2]*x*x)",Kitty.GetMomBinLowEdge(0),kCorrection);
    fMMX->SetParameter(0,1);
    fMMX->SetParLimits(0,-40,40);
    fMMX->SetParameter(1,0);
    fMMX->SetParLimits(1,-40,40);
    fMMX->FixParameter(2,fabs(f0*d0*FmToNu*FmToNu*0.5));
    //the weight of the f^2 factor with respect to F1F2 (after correction)
    double wf2=0;
    const double if0=f0?1./f0:0;
    Float_t ntBuffer[16];
    for(unsigned uMom=0; uMom<nMom; uMom++){
      double MOM = Kitty.GetMomentum(uMom);
      double F1 = gsl_sf_dawson(2.*MOM*r0*FmToNu)/(2.*MOM*r0*FmToNu);
      double F2 = (1.-exp(-4.*MOM*MOM*r0*r0*FmToNu*FmToNu))/(2.*MOM*r0*FmToNu);
      complex<double> Val_f = pow(if0/FmToNu+0.5*d0*FmToNu*MOM*MOM-i*MOM,-1.);
      double v_f2 = 0.5*pow(abs(Val_f)/(r0*FmToNu),2.);
      double F1F2 = 2*real(Val_f)*F1/(sqrt(Pi)*r0*FmToNu)-imag(Val_f)*F2/(r0*FmToNu);
      //double Ck_Ledni = 1.+v_f2+F1F2;
      double Ck_Cats = Kitty.GetCorrFun(uMom);
      //the correction factor we want to study
      double MMX = (Ck_Cats-1.-F1F2)/v_f2;
      hMMX->SetBinContent(uMom+1,MMX);
      hMMX->SetBinError(uMom+1,0.001);
      wf2 += fabs((Ck_Cats-1.-F1F2)/F1F2);
    }
    wf2 /= double(nMom);
    //hMMX->Fit(fMMX,"Q, S, N, R, M");
    ROOT::Math::MinimizerOptions MinOpt;
    MinOpt.SetMinimizerType("Minuit2");
    MinOpt.SetPrintLevel(0);
    DLM_FitHisto(hMMX, fMMX, "Q, S, N, R, M", "", &MinOpt);
    double csmall = (1-(d0)/(2*sqrt(Pi)*r0));
    ntBuffer[0] = r0;
    ntBuffer[1] = f0;
    ntBuffer[2] = d0;
    ntBuffer[3] = fMMX->GetParameter(0);
    ntBuffer[4] = fMMX->GetParameter(1);
    ntBuffer[5] = hMMX->GetBinContent(1);
    ntBuffer[6] = hMMX->GetBinContent(nMom);
    ntBuffer[7] = csmall;
    ntBuffer[8] = V1;
    ntBuffer[9] = mu1;
    ntBuffer[10] = V2;
    ntBuffer[11] = mu2;
    ntBuffer[12] = s2;
    ntBuffer[13] = wf2;
    ntBuffer[14] = Get_MMclass(r0,f0,d0);
    ntBuffer[15] = DimiFlag;

    ntMM->Fill(ntBuffer);
    delete hMMX; delete fMMX;
    delete [] MomBins;
}

void Ledni_SmallRad_Random(const unsigned SEED, const unsigned NumIter){
//boosts classes 21,23,31,32
// *radius fixed to 1.08 done
// *the boosters only apply to the above 4 classes, never to others done
// *the Gaussian V1 has to be negative done
// *for a Gaussian mu1<1 done
// *for Yukawa V1<80  done
const bool EmmaBias = true;
  const double kMin = 0;
  const double kMax = 150;
  const double kCorrection = 120;
  const double kStep = 3;
  const unsigned nMom = TMath::Nint(kMax/kStep);

  //const double V1_min = -1000;
  //const double V1_max = 0;
  //const double mu1_min = 0.6;
  //const double mu1_max = 1.2;
  //const double V2_min = 0;
  //const double V2_max = 5000;
  //const double mu2_min = 0.3;
  //const double mu2_max = 0.5;

  //const double r_min = 1.00;
  //const double r_max = 1.50;\

  double Radii[8];
  Radii[0] = 0.9;
  Radii[1] = 1.0;
  Radii[2] = 1.2;
  Radii[3] = 1.5;
  Radii[4] = 2.0;
  Radii[5] = 2.5;
  Radii[6] = 3.0;
  Radii[7] = 4.0;

  //separtion based on:
  //f0>0
  //r0>d0
  //|f0|>2*d0 (BS cond) -> satisfied for 2 dominant or 3 surpressed
  //r0>|f|/2

  //interesting classes:
  //+/- AB, with A,B = 1,2,3 (r0,f0,d0)
  //A = which one is dominant (e.g. 1 if r0>d0 && r0>|f|/2)
  //B = which one is surpressed (e.g. 3 if d0<r0 && d0<|f0|/2)
  //N.B.
  //ALL POSSIBLE CLASSES:
  // 12: shallow attraction, large radius. Ledni should be perfect here   (OK)
  // 13: moderate attraction, large radius. Ledni should work well        (OK)
  // 21: strong attraction, small radius. Ledni could fail                (!)
  // 23: strong attraction, small range. Ledni could work                 (?)
  // 31: attraction, small radius, large range. Ledni will fail           (!!)
  // 32: weak attraction, small radius, large range. Ledni could fail     (!)
  //-12: shallow repulsion, large radius. Ledni could work                (?)
  //-13: bound state, large radius. Ledni will work                       (OK)
  //-21: bound state, small radius. Ledni could fail                      (!)
  //-23: deep bound state, moderate raidus. Ledni could work              (?)
  //-31: repulsion, small radius. Ledni will fail                         (!!)
  //-32: shallow repulsion. Ledni will fail                               (!!)
         //--> actually that one works sometimes, in case r0>d0

  TRandom3 rangen(SEED);

  CATSparameters sPars(CATSparameters::tSource,1,true);
  sPars.SetParameter(0,2.5);
  CATS Kitty_SE;
  Kitty_SE.SetMomBins(nMom,kMin,kMax);
  Kitty_SE.SetAnaSource(GaussSource, sPars);
  Kitty_SE.SetUseAnalyticSource(true);
  Kitty_SE.SetQ1Q2(0);
  Kitty_SE.SetQuantumStatistics(false);
  Kitty_SE.SetRedMass(Mass_L*0.5);
  Kitty_SE.SetNumChannels(1);
  Kitty_SE.SetNumPW(0,1);
  Kitty_SE.SetSpin(0,0);
  Kitty_SE.SetChannelWeight(0, 1.);
  Kitty_SE.SetEpsilonConv(5e-9);
  Kitty_SE.SetEpsilonProp(5e-9);
  Kitty_SE.SetMaxRad(64);
  Kitty_SE.SetMaxRho(24);
  Kitty_SE.SetNotifications(CATS::nWarning);
  Kitty_SE.KillTheCat();
  CATSparameters pPars(CATSparameters::tPotential,5,true);
  //Kitty_SE.SetShortRangePotential(0,0,DoubleGaussSum,pPars);
  Kitty_SE.SetShortRangePotential(0,0,YukawaDimiCore,pPars);
  double c_f0,c_if0,c_d0,V_1,V_2,mu_1,mu_2,r_0,s_2,l_f0,l_if0,l_d0;
  const TString OutputFolder = TString::Format("%s/MMM/Emma2/",GetCernBoxDimi());
  TFile fOutput(OutputFolder+TString::Format("fOut_S%u_I%u.root",SEED,NumIter),"recreate");
  TNtuple* ntMM = new TNtuple("ntMM", "ntMM","r0:f0:d0:MMA:MMB:C_X0:C_X150:L_X0:V1:mu1:V2:mu2:s2:wf2:Class");
  Float_t ntBuffer[15];
  TH1F* hClassCount = new TH1F("hClassCount","hClassCount",100,-50.5,49.5);
  TH2F* hif0d0Count = new TH2F("hf0d0Count","hf0d0Count",20,-10,10,20,0,40);
  int LastClass=0;
  int SillyPhaseShifts = 0;
  for(unsigned uIter=0; uIter<NumIter; uIter++){
    printf("\r\033[K Progress %5u from %5u\n",uIter,NumIter);
    //printf("Progress %5u from %5u: ",uIter,NumIter);

    //try to enhance the statistics here
    //if((1.+2.*c_d0/c_f0)<0){
    //  V_1 = rangen.Gaus(V_1,30);
    //  V_2 = rangen.Gaus(V_2,30);
    //  mu_1 = rangen.Gaus(mu_1,0.02);
    //  mu_2 = rangen.Gaus(mu_2,0.02);
    //}
    /*
    if(c_d0>5&&c_f0>0){
      V_1 = rangen.Gaus(V_1,15);
      V_2 = rangen.Gaus(V_2,15);
      mu_1 = rangen.Gaus(mu_1,0.01);
      mu_2 = rangen.Gaus(mu_2,0.01);
    }
    else{
      //try to enhance bound states a bit
      if(rangen.Uniform()<0.2){
        V_1 = rangen.Gaus(-400,200);
        V_2 = rangen.Gaus(600,300);
      }
      else{
        V_1 = rangen.Gaus(-100,75);
        V_2 = rangen.Uniform(V2_min,V2_max);
      }
      mu_1 = rangen.Uniform(mu1_min,mu1_max);
      mu_2 = rangen.Uniform(mu2_min,mu2_max);
    }
    */

//V_1 = -144.89;
//V_2 = 127.87;
//mu_1 = 1.0;
//mu_2 = 0.45;
//V_1 = -144.5;CLASS
//V_2 = 520.0;
//mu_1 = 2.11;
//mu_2 = 0.54;
    int AvgN = uIter/(EmmaBias?6:12);
    int ClassN = hClassCount->GetBinContent(hClassCount->FindBin(double(LastClass)));
    int AvgSP = uIter/400;//based on num bins
    int if0d0N = hif0d0Count->GetBinContent(hif0d0Count->FindBin(l_if0,l_d0));
    bool SUPER_BOOST = (EmmaBias&&l_f0>0.65&&l_f0<1.25&&l_d0>6.2&&l_d0<17.8);
    //printf(" (Class=%i N%i A%i)",LastClass,ClassN,AvgN);
    //in case of EmmaBias, it activates for 31,32 or SUPER_BOOST regardless on the amount of entires for this class
    //in case of EmmaBias, classes 21 and 23 activate in case there is two little amount of them
    if((ClassN<AvgN||if0d0N<AvgSP||SUPER_BOOST||LastClass>30)&&SillyPhaseShifts<3&&(LastClass>20||!EmmaBias)){
      double FactorSP=1;double FactorR=1;
      if(ClassN<AvgN/2){FactorSP*=0.5;FactorR*=0.5;}
      else if(if0d0N<AvgSP/2){FactorSP*=0.5;}
      if(ClassN<AvgN/4){FactorSP*=0.5;FactorR*=0.5;}
      else if(if0d0N<AvgSP/4){FactorSP*=0.5;}
      if(SUPER_BOOST)FactorSP=1./6.;
      //if we are close (30%, 2x for d0 in upper direction) to 0.95 8.9 which is the ledni fit result (real)
      //printf(" -> EB%i lf0=%.2f ld0=%.2f ", EmmaBias,l_f0,l_d0);
      V_1 = rangen.Gaus(V_1,s_2>0?3*FactorSP:12*FactorSP);
      V_2 = rangen.Gaus(V_2,s_2>0?12*FactorSP:16*FactorSP);
      //mu_2 = 0.45;
      //mu_2 = rangen.Gaus(mu_2,0.01*Factor);
      mu_1 = rangen.Gaus(mu_1,0.02*FactorSP);

      if(EmmaBias) r_0 = 1.08;
      else if(ClassN<AvgN) r_0 = r_0;
      else r_0 = Radii[rangen.Integer(8)];
    }
    else{
      if(rangen.Uniform()<0.5) s_2 = 0.1;
      else s_2 = 0.0;
      if(s_2>0){
        //Yukawa
        if(EmmaBias) V_1 = rangen.Uniform(5,80);
        else V_1 = rangen.Uniform(5,150);
        mu_1 = rangen.Uniform(0.6,2.1);
        V_2 = rangen.Uniform(5,600);
        s_2 = 0.1;
      }
      else{
        //V_1 = rangen.Uniform(-2500,500);
        if(EmmaBias) V_1 = rangen.Gaus(-1100,0);
        else V_1 = rangen.Gaus(-1100,700);
        //V_2 = rangen.Uniform(-500,5500);
        V_2 = rangen.Gaus(1400,800);
        if(EmmaBias) mu_1 = rangen.Uniform(0.5,1.0);
        else mu_1 = rangen.Uniform(0.5,1.6);
        s_2 = 0;
      }
      if(rangen.Uniform()<0.75) mu_2 = 0.4;
      else mu_2 = 0.7;
      if(EmmaBias) r_0 = 1.08;
      else r_0 = Radii[rangen.Integer(8)];
    }
         if(V_1>0&&V_2>0)printf("   V1=%5.0f mu1=%.2f V2=%5.0f mu2=%.2f s2=%.1f r0=%.2f", V_1,mu_1,V_2,mu_2,s_2,r_0);
    else if(V_1>0&&V_2<0)printf("   V1=%5.0f mu1=%.2f V2=%4.0f mu2=%.2f s2=%.1f r0=%.2f", V_1,mu_1,V_2,mu_2,s_2,r_0);
    else if(V_1<0&&V_2>0)printf("   V1=%4.0f mu1=%.2f V2=%5.0f mu2=%.2f s2=%.1f r0=%.2f", V_1,mu_1,V_2,mu_2,s_2,r_0);
    else                 printf("   V1=%4.0f mu1=%.2f V2=%4.0f mu2=%.2f s2=%.1f r0=%.2f", V_1,mu_1,V_2,mu_2,s_2,r_0);

    c_f0 = 0.0;
    c_d0 = 0.0;
    //Kitty_SE.SetShortRangePotential(0,0,0,V_1);
    //Kitty_SE.SetShortRangePotential(0,0,1,mu_1);
    //Kitty_SE.SetShortRangePotential(0,0,2,V_2);
    //Kitty_SE.SetShortRangePotential(0,0,3,mu_2);

    if(s_2==0) Kitty_SE.SetShortRangePotential(0,0,DoubleGaussSum,pPars);
    else if(s_2==-1) Kitty_SE.SetShortRangePotential(0,0,GaussExpSum,pPars);
    else Kitty_SE.SetShortRangePotential(0,0,YukawaDimiCore,pPars);
    Kitty_SE.SetShortRangePotential(0,0,0,V_1);
    Kitty_SE.SetShortRangePotential(0,0,1,mu_1);
    Kitty_SE.SetShortRangePotential(0,0,2,V_2);
    //Kitty_SE.SetShortRangePotential(0,0,3,rangen.Uniform(0.3,0.5));
    Kitty_SE.SetShortRangePotential(0,0,3,mu_2);
    //Kitty_SE.SetShortRangePotential(0,0,4,rangen.Uniform(0.05,0.15));
    Kitty_SE.SetShortRangePotential(0,0,4,s_2);

    //Kitty_SE.SetShortRangePotential(0,0,0,50.);
    //Kitty_SE.SetShortRangePotential(0,0,1,1.);
    //Kitty_SE.SetShortRangePotential(0,0,2,200.);
    //Kitty_SE.SetShortRangePotential(0,0,3,0.45);
    //Kitty_SE.SetShortRangePotential(0,0,4,0.1);

    Kitty_SE.SetAnaSource(0,r_0,EmmaBias);
    Kitty_SE.KillTheCat();
    TH1F* hDummy; TF1* fDummy;
    //potential with silly phase shifts
    if(!Eval_ScattParameters(Kitty_SE,c_f0,c_d0,hDummy,fDummy)){
      uIter--;
      if(hDummy) delete hDummy;
      if(fDummy) delete fDummy;
      SillyPhaseShifts++;
      cout << flush;
      cout << "\033[F";
      cout << flush;
//cout << endl;
      continue;
    }
    SillyPhaseShifts=0;
    c_if0 = 1./c_f0;
//hDummy->Write();fDummy->Write();
    delete hDummy; delete fDummy;
    if(fabs(c_if0)>100||fabs(c_f0)>100||fabs(c_d0)>100){
      uIter--;
      cout << flush;
      cout << "\033[F";
      cout << flush;
//cout << endl;
      continue;
    }

    if(c_f0>0)printf(" =uIter> f0=%4.2f d0=%3.2f--------",c_f0,c_d0);
    else printf(" => f0=%3.2f d0=%3.2f--------",c_f0,c_d0);
    if(SUPER_BOOST){
      printf("SUPER_BOOST");
      //usleep(2000e3);
    }
    cout << flush;
    cout << "\033[F";
    cout << flush;
//cout << endl;

    TH1F* hXXX = new TH1F("hXXX","hXXX",nMom,kMin,kMax);
    TF1* fXXX = new TF1("fXXX","[0]+[1]*sqrt(1.+[2]*x*x)",kMin,kCorrection);
    fXXX->SetParameter(0,1);
    fXXX->SetParLimits(0,-40,40);
    fXXX->SetParameter(1,0);
    fXXX->SetParLimits(1,-40,40);
    fXXX->FixParameter(2,fabs(c_f0*c_d0*FmToNu*FmToNu*0.5));
    //the weight of the f^2 factor with respect to F1F2 (after correction)
    double wf2=0;
    for(unsigned uMom=0; uMom<nMom; uMom++){
      double MOM = Kitty_SE.GetMomentum(uMom);
      double F1 = gsl_sf_dawson(2.*MOM*r_0*FmToNu)/(2.*MOM*r_0*FmToNu);
      double F2 = (1.-exp(-4.*MOM*MOM*r_0*r_0*FmToNu*FmToNu))/(2.*MOM*r_0*FmToNu);
      complex<double> Val_f = pow(c_if0/FmToNu+0.5*c_d0*FmToNu*MOM*MOM-i*MOM,-1.);
      double v_f2 = 0.5*pow(abs(Val_f)/(r_0*FmToNu),2.);
      double F1F2 = 2*real(Val_f)*F1/(sqrt(Pi)*r_0*FmToNu)-imag(Val_f)*F2/(r_0*FmToNu);
      //double Ck_Ledni = 1.+v_f2+F1F2;
      double Ck_Cats = Kitty_SE.GetCorrFun(uMom);
      //printf("k=%.0f Cats=%.2f Ledni=%.2f\n", MOM, Ck_Cats, 1.+v_f2+F1F2);
      //usleep(100e3);
      //the correction factor we want to study
      double XXX = (Ck_Cats-1.-F1F2)/v_f2;
      hXXX->SetBinContent(uMom+1,XXX);
      hXXX->SetBinError(uMom+1,0.001);
      wf2 += fabs((Ck_Cats-1.-F1F2)/F1F2);
    }
    wf2 /= double(nMom);
    hXXX->Fit(fXXX,"Q, S, N, R, M");
    double csmall = (1-(c_d0)/(2*sqrt(Pi)*r_0));
    ntBuffer[0] = r_0;
    ntBuffer[1] = c_f0;
    ntBuffer[2] = c_d0;
    ntBuffer[3] = fXXX->GetParameter(0);
    ntBuffer[4] = fXXX->GetParameter(1);
    ntBuffer[5] = hXXX->GetBinContent(1);
    ntBuffer[6] = hXXX->GetBinContent(nMom);
    ntBuffer[7] = csmall;
    ntBuffer[8] = V_1;
    ntBuffer[9] = mu_1;
    ntBuffer[10] = V_2;
    ntBuffer[11] = mu_2;
    ntBuffer[12] = s_2;
    ntBuffer[13] = wf2;


    //f0>0
    //r0>d0
    //|f0|>2*|d0| (BS cond) -> satisfied for 2 dominant or 3 surpressed
    //r0>|f|/2
    //interesting classes:
    //+/- AB, with A,B = 1,2,3 (r0,f0,d0)
    //A = which one is dominant (e.g. 1 if r0>d0 && r0>|f|/2)
    //B = which one is surpressed (e.g. 3 if d0<r0 && d0<|f0|/2)
    int CLASS = 0;
    if(r_0<fabs(c_d0)&&r_0<fabs(c_f0)*0.5) CLASS = 1;
    else if(fabs(c_f0)<2.*r_0&&fabs(c_f0)<2.*fabs(c_d0)) CLASS = 2;
    else if(fabs(c_d0)<r_0&&fabs(c_d0)<fabs(c_f0)*0.5) CLASS = 3;

    if(r_0>fabs(c_d0)&&r_0>fabs(c_f0)*0.5) CLASS += 10;
    else if(fabs(c_f0)>2.*r_0&&fabs(c_f0)>2.*fabs(c_d0)) CLASS += 20;
    else if(fabs(c_d0)>r_0&&fabs(c_d0)>fabs(c_f0)*0.5) CLASS += 30;

    CLASS = c_f0>0?CLASS:-CLASS;
    ntBuffer[14] = CLASS;
    //printf("BEFORE: %f\n",hClassCount->GetBinContent(hClassCount->FindBin(double(CLASS))));
    hClassCount->Fill(ntBuffer[14]);
    //printf("AFTER: %f\n",hClassCount->GetBinContent(hClassCount->FindBin(double(CLASS))));
    LastClass = CLASS;

    hif0d0Count->Fill(c_if0,c_d0);
    l_f0=c_f0;
    l_if0=c_if0;
    l_d0=c_d0;

    ntMM->Fill(ntBuffer);
    //hXXX->Write();
    //fXXX->Write();
    delete hXXX; delete fXXX;
  }
  printf("\n");

  ntMM->Write();
  hClassCount->Write();
  hif0d0Count->Write();
  delete ntMM;
  delete hClassCount;
  delete hif0d0Count;
}

void SelectEmmaPotential(){
  const double f_0 = 0.95;
  const double d_0 = 8.9;
  const double delta = 2./100.;
  const double df0 = fabs(f_0*delta);
  const double dd0 = fabs(d_0*delta);

  TString InputFileName = TString::Format("%s/MMM/Emma2/I4K.root",GetCernBoxDimi());
  TFile InputFile(InputFileName,"read");
  TNtuple* ntMM = (TNtuple*)InputFile.Get("ntMM");
  unsigned NumEntries = ntMM->GetEntries();
  Float_t ft_r0;  ntMM->SetBranchAddress("r0",&ft_r0);
  Float_t ft_f0;  ntMM->SetBranchAddress("f0",&ft_f0);
  Float_t ft_d0;  ntMM->SetBranchAddress("d0",&ft_d0);
  Float_t ft_V1;  ntMM->SetBranchAddress("V1",&ft_V1);
  Float_t ft_mu1; ntMM->SetBranchAddress("mu1",&ft_mu1);
  Float_t ft_V2;  ntMM->SetBranchAddress("V2",&ft_V2);
  Float_t ft_mu2; ntMM->SetBranchAddress("mu2",&ft_mu2);
  Float_t ft_s2;  ntMM->SetBranchAddress("s2",&ft_s2);
  for(unsigned uEntry; uEntry<NumEntries; uEntry++){
    ntMM->GetEntry(uEntry);
    if(ft_f0>f_0-df0&&ft_f0<f_0+df0&&ft_d0>d_0-dd0&&ft_d0<d_0+dd0){
      printf("f0=%.3f; d0=%.2f;\n",ft_f0,ft_d0);
      if(ft_s2>0.01) printf(" Yukawa potential with repulsive core\n");
      else printf(" Double Gaussian potential\n");
      printf("   V1 = %.6e\n",ft_V1);
      printf("   mu1 = %.6e\n",ft_mu1);
      printf("   V2 = %.6e\n",ft_V2);
      printf("   mu2 = %.6e\n",ft_mu2);
      printf("---------------------------------\n");
    }
  }
}

//evaluated w/o coulomb and for red mass of p-phi (so around 490 MeV)
void OkayishStartingPars( const TString Potential, const double& f0, const double d0,
                          double& V1, double& mu1, double& V2, double& mu2, double& s2){
  V1 = 0.1;
  mu1 = 1;
  V2 = 0.1;
  mu2 = 0.5;
  if(f0>0&&fabs(f0)>fabs(2.*d0)){
    if(Potential=="DoubleGaussSum"){

    }
    else if(Potential=="YukawaDimiCore"){
      s2 = mu2*0.2;
    }
    else if(Potential=="Gaussian"){

    }
    else if(Potential=="Yukawa"){

    }
    else{
      V1 = 0;
      mu1 = 0;
      V2 = 0;
      mu2 = 0;
      s2 = 0;
    }
  }
  else if(f0>0&&fabs(f0)<=fabs(2.*d0)){
    if(Potential=="DoubleGaussSum"){
      //f0 = 0.993 fm
      //d0 = 3.97 fm
      V1 = -8.576280e+02;
      mu1 = 6.943480e-01;
      V2 = 7.481622e+03;
      mu2 = 3.903410e-01;
      s2 = 0;
    }
    else if(Potential=="YukawaDimiCore"){
      //f0 = 1.027 fm
      //d0 = 4.19 fm
      V1 = 3.455340e+00;
      mu1 = 3.881920e-01;
      V2 = 2.211957e+02;
      mu2 = 6.202096e-01;
      s2 = mu2*0.2;
    }
    else if(Potential=="Gaussian"){
      //f0 = 0.976 fm
      //d0 = 3.77 fm
      V1 = -2.469257e+01;
      mu1 = 1.301668e+00;
      V2 = 0.000000e+00;
      mu2 = 0.000000e+00;
      s2 = -1;
    }
    else if(Potential=="Yukawa"){
      s2 = -2;
    }
    else{
      V1 = 0;
      mu1 = 0;
      V2 = 0;
      mu2 = 0;
      s2 = 0;
    }
  }
  else if(f0<0&&fabs(f0)>fabs(2.*d0)){
    if(Potential=="DoubleGaussSum"){
      //f0 = -3.976 fm
      //d0 = 1.01 fm
      V1 = -7.414734e+01;
      mu1 = 1.111177e+00;
      V2 = -6.127154e+02;
      mu2 = 3.208735e-01;
      s2 = 0;
    }
    else if(Potential=="YukawaDimiCore"){
      //f0 = -4.094 fm
      //d0 = 1.04 fm
      V1 = 1.246311e+01;
      mu1 = 2.290515e-01;
      V2 = -8.227462e+01;
      mu2 = 9.258949e-01;
      s2 = mu2*0.2;
    }
    else if(Potential=="Gaussian"){
      //f0 = -3.970 fm
      //d0 = 1.05 fm
      V1 = -1.881959e+02;
      mu1 = 8.635644e-01;
      V2 = 0.000000e+00;
      mu2 = 0.000000e+00;
      s2 = -1;
    }
    else if(Potential=="Yukawa"){
      //f0 = -3.837 fm
      //d0 = 1.05 fm
      //V1 = 8.792615e-01;
      //mu1 = 1.045134e+02;
      //V2 = 0.000000e+00;
      //mu2 = 0.000000e+00;
      s2 = -2;
    }
    else{
      V1 = 0;
      mu1 = 0;
      V2 = 0;
      mu2 = 0;
      s2 = 0;
    }

  }
  else if(f0<0&&fabs(f0)<=fabs(2.*d0)){
    if(Potential=="DoubleGaussSum"){
      //f0 = -0.948 fm
      //d0 = 0.95 fm
      V1 = -1.418785e+03;
      mu1 = 1.123364e+00;
      V2 = -1.147570e+03;
      mu2 = 1.031657e+00;
      s2 = 0;
    }
    else if(Potential=="YukawaDimiCore"){
      //f0 = -1.048 fm
      //d0 = 1.04 fm
      V1 = 1.945222e+01;
      mu1 = 4.377242e-01;
      V2 = 1.085863e+03;
      mu2 = 1.111706e+00;
      s2 = mu2*0.2;
    }
    else if(Potential=="Gaussian"){
      //f0 = -1.054 fm
      //d0 = 1.02 fm
      V1 = -1.260122e+02;
      mu1 = 1.991708e+00;
      V2 = 0.000000e+00;
      mu2 = 0.000000e+00;
      s2 = -1;
    }
    else if(Potential=="Yukawa"){
      //f0 = -1.042 fm
      //d0 = 1.04 fm
      V1 = 7.584738e+00;
      mu1 = 5.422729e-01;
      V2 = 0.000000e+00;
      mu2 = 0.000000e+00;
      s2 = -2;
    }
    else{
      V1 = 0;
      mu1 = 0;
      V2 = 0;
      mu2 = 0;
      s2 = 0;
    }
  }
  else{
    printf(" OkayishStartingPars says f*** !?\n");
  }




}

//Potential==Dynamic, means we try to do the double Gaussian, if it fails we repeat it all with Yukawa
//FineTune < 1 to decrease the intrisic minimum step for the current potential
void ManufacturePotential(const double f0, const double df0,
                                const double d0, const double dd0, const double* Radii, const unsigned NumRad,
                                double& V1, double& mu1, double& V2, double& mu2,
                                const TString Potential, const TString OutputFolder, const int& SEED=11, const double* StartPars=NULL,
                                const double RedMass=500, const double FineTune=1){


  //b = best, l = last, d = difference to desired
  TString CurrentPot = Potential;
  if(Potential=="Dynamic") CurrentPot = "DoubleGaussSum";
  double V_1,bV_1;
  double V_2,bV_2;
  double mu_1,bmu_1;
  double mu_2,bmu_2;
  double s_2;
  double f_0,lf_0,bf_0,df_0,bdf_0;
  double d_0,ld_0,bd_0,dd_0,bdd_0;
  bool Starting = true;
  bool Started = false;
  //dist,f0,d0,V1,mu1,V2,mu2: used when VeryStuck>=2 resets the whole thing
  double FallBack[8];FallBack[0]=1e16;
  f_0 = 0; d_0 = 0;
  bV_1 = 0; bV_2 = 0;
  bmu_1 = 1.5; bmu_2 = 0.5;
  bdf_0 = 1e6; bdd_0 = 1e6;
  bf_0 = f0; bd_0 = d0;
  lf_0=f0; ld_0=d0;
  const unsigned MaxIter = 10000;
  unsigned uIter=0;
  TRandom3 rangen(SEED);
  CATSparameters sPars(CATSparameters::tSource,1,true);
  sPars.SetParameter(0,2.5);
  const double kStepFine=3;
  const double kStepCoarse=10;
  const unsigned nMomFine = 40;
  const unsigned nMomCoarse = 18;
  const unsigned nMom = nMomFine+nMomCoarse;
  double* MomBins = new double [nMom+1];
  for(unsigned uMom=0; uMom<nMomFine; uMom++){
    MomBins[uMom] = kStepFine*double(uMom);
  }
  for(unsigned uMom=0; uMom<=nMomCoarse; uMom++){
    MomBins[nMomFine+uMom] = kStepFine*double(nMomFine)+kStepCoarse*double(uMom);
  }
  CATS Kitty_SE;
  Kitty_SE.SetMomBins(nMom,MomBins);
  Kitty_SE.SetAnaSource(GaussSource, sPars);
  Kitty_SE.SetUseAnalyticSource(true);
  Kitty_SE.SetQ1Q2(0);
  Kitty_SE.SetQuantumStatistics(false);
  Kitty_SE.SetRedMass( RedMass );
  //Kitty_SE.SetRedMass(Mass_L*0.5);
//Kitty_SE.SetRedMass(488.6);//pphi
  Kitty_SE.SetNumChannels(1);
  Kitty_SE.SetNumPW(0,1);
  Kitty_SE.SetSpin(0,0);
  Kitty_SE.SetChannelWeight(0, 1.);
  Kitty_SE.SetEpsilonConv(1e-8);
  Kitty_SE.SetEpsilonProp(1e-8);
  //Kitty_SE.SetMaxRad(96);
  //Kitty_SE.SetMaxRho(32);
  Kitty_SE.SetNotifications(CATS::nSilent);
  //Kitty_SE.SetGridEpsilon(1./512.);
  double RadForGrid;
  Kitty_SE.KillTheCat();
  CATSparameters pPars(CATSparameters::tPotential,5,true);
  if(CurrentPot=="YukawaDimiCore") Kitty_SE.SetShortRangePotential(0,0,YukawaDimiCore,pPars);
  else if(CurrentPot=="Gassian") Kitty_SE.SetShortRangePotential(0,0,Gaussian,pPars);
  else if(CurrentPot=="Yukawa") Kitty_SE.SetShortRangePotential(0,0,YukawaDimiSmooth,pPars);
  else Kitty_SE.SetShortRangePotential(0,0,DoubleGaussSum,pPars);
  bool GoodGoing=false;
  double dist=1e16;
  double bdist=1e16;
  //the scale to which 1fm roughly correponds
  const double dist_unit = sqrt(0.5/df0/df0+0.5/dd0/dd0);

  //Fluctuations: should be the same order as bdist
  //if too large -> decrease Convergence
  //if too small -> increse Convergence
  const unsigned fluctN = 8;
  double fluct[fluctN];
  unsigned ufluct = 0;
  double avg_fluct;
  const unsigned dirN = 32;
  //how many of the last dirN iterations moved towards the true solution,
  //which is defined as +/-45 deg angle on the 1/f0 d0 plane (1/4 of the angle)
  //we expect 4 movements on avg, if its less than that we increase Convergence
  bool dir[dirN];
  unsigned udir = 0;
  double avg_dir;
  //double Adjust = 10;
  double Convergence[5];
  for(unsigned uPar=0; uPar<5; uPar++)Convergence[uPar]=1;
  unsigned Stuck = 0;
  unsigned BadPhaseShifts = 0;
  unsigned UnstuckCounter = 0;
  const unsigned StuckLimit = 64;
  //how many iterations to do by changing only single pars
  const unsigned UnstuckPerPar = 32*2;
  //1 is V1, 2 is mu1, 3 is V2 and 4 is mu2. 0 the normal mode
  unsigned SinglePar = 0;
  for(unsigned uf=0; uf<fluctN; uf++) fluct[uf]=-1000*dist_unit;
  for(unsigned ud=0; ud<dirN; ud++) dir[ud]=true;
//TFile fDump(TString::Format("%s/OtherTasks/ManufactureYukawaPotential/fDump.root",GetFemtoOutputFolder()),"recreate");
const bool DEBUG = false;
  unsigned VeryStuck=0;
  TFile fOutput(TString::Format("%s/fOut_%.2f_%.1f_%s.root",OutputFolder.Data(),f0,d0,Potential.Data()),"recreate");
  TNtuple* ntMM = new TNtuple("ntMM", "ntMM","r0:f0:d0:MMA:MMB:C_X0:C_X150:L_X0:V1:mu1:V2:mu2:s2:wf2:Class:Manufactured");
  TH1F* hClassCount = new TH1F("hClassCount","hClassCount",100,-50.5,49.5);
  TH2F* hif0d0Count = new TH2F("hf0d0Count","hf0d0Count",20,-10,10,20,0,40);
  bool FoundIt = false;
  while(uIter<MaxIter&&!FoundIt){
    if(!DEBUG){
      printf("\r\033[K Goal: %s (f0,d0) [Achieved]: (%.3f,%.3f)+/-(%.3f,%.3f) [%3.0f%%,%3.0f%%], Break (%5u/%5u)",Potential.Data(),f0,d0,fabs(df0),fabs(dd0),fabs(df0/bdf_0*100.),fabs(dd0/bdd_0*100.),uIter,MaxIter);
    }
    else printf("\nProgress: Break (%5u/%5u), Precision f0(%3.0f%%) d0(%3.0f%%)",uIter,MaxIter,fabs(df0/bdf_0*100.),fabs(dd0/bdd_0*100.));
    //printf("Progress %5u from %5u, bd=%.1f(%.1f)\n",uIter,MaxIter,bdist,sqrt(0.5*bdf_0*bdf_0/df0/df0+0.5*bdd_0*bdd_0/dd0/dd0));

    avg_fluct = 0;
    for(unsigned uf=0; uf<fluctN; uf++) avg_fluct+=fluct[uf];
    avg_fluct /= double(fluctN);
    double desired_fluct = bdist;

    avg_dir = 0;
    for(unsigned ud=0; ud<dirN; ud++) avg_dir+=double(dir[ud]);
    avg_dir /= double(dirN);

    if(DEBUG){
      printf("\n B: Stuck=%u(%u); UnstuckCounter=%u; SinglePar=%u",Stuck,VeryStuck,UnstuckCounter,SinglePar);
    }

    unsigned Npars = 4;
    if(CurrentPot=="YukawaDimiCore") Npars=4;
    if(CurrentPot=="Gaussian") Npars=2;
    if(CurrentPot=="Yukawa") Npars=2;
    else Npars=4;

    if(Stuck==StuckLimit&&SinglePar==0){
      VeryStuck=1;
      UnstuckCounter=0;
      Stuck=0;
      SinglePar=1;
      for(unsigned uPar=1;uPar<=Npars;uPar++)Convergence[uPar]*=Convergence[0];
      Convergence[0]=1;
      desired_fluct /= 2.;
    }
    if(BadPhaseShifts==10) {VeryStuck=3;Stuck=100000;UnstuckCounter=UnstuckPerPar*Npars;}
    if(UnstuckCounter==UnstuckPerPar*Npars){
      //if VeryStuck goes two times: reset and start from scratch
      if(Stuck<UnstuckCounter){
        VeryStuck=0;
        UnstuckCounter=0;
        Stuck=0;
        SinglePar=0;
        Convergence[0]=0;
        //make it the max value
        for(unsigned uPar=1;uPar<=Npars;uPar++)
          if(Convergence[uPar]>Convergence[0])Convergence[0]=Convergence[uPar];
        for(unsigned uPar=1;uPar<=Npars;uPar++)Convergence[uPar]/=Convergence[0];
      }
      else if(VeryStuck<3){
        VeryStuck = 1+Stuck/UnstuckCounter;
        desired_fluct /= (2.*VeryStuck);
        UnstuckCounter=0;
        //Stuck=0;
        SinglePar=1;
      }
      else{
        if(Potential=="Dynamic"){
          if(CurrentPot=="DoubleGaussSum") CurrentPot="YukawaDimiCore";
          else CurrentPot="DoubleGaussSum";
        }
        if(!DEBUG){
          //cout<<flush;
          printf("\n\n----- RESET (%s) -----\n",CurrentPot.Data());
          printf("\r\033[K Progress: Break (%5u/%5u)",uIter,MaxIter);
        }
        else{
          printf(" --- RESET (%s) ---\n",CurrentPot.Data());
        }

        //make it the max value
        for(unsigned uPar=0;uPar<=Npars;uPar++)
          Convergence[uPar]=1;
        if(bdist<FallBack[0]){
          FallBack[0] = bdist;
          FallBack[1] = bf_0;
          FallBack[2] = bd_0;
          FallBack[3] = bV_1;
          FallBack[4] = bmu_1;
          FallBack[5] = bV_2;
          FallBack[6] = bmu_2;
          FallBack[7] = s_2;
        }
        VeryStuck=0;
        UnstuckCounter=0;
        Stuck=0;
        SinglePar=0;
        bdist = 1e16;
        f_0 = 0; d_0 = 0;
        bV_1 = 0; bV_2 = 0;
        bmu_1 = 1.5; bmu_2 = 0.5;
        bdf_0 = 1e6; bdd_0 = 1e6;
        bf_0 = f0; bd_0 = d0;
        for(unsigned uf=0; uf<fluctN; uf++) fluct[uf]=-1000*dist_unit;
        for(unsigned ud=0; ud<dirN; ud++) dir[ud]=true;
        Starting = true;
        Started = false;
      }
    }
    if(SinglePar){
      SinglePar = 1+UnstuckCounter/UnstuckPerPar;
      if(Stuck<UnstuckPerPar*Npars) printf(" stuck, but working on it...");
      else printf(" very stuck, Jesus is working on it...");
      desired_fluct = bdist/(2.*VeryStuck);
      //desired_fluct = 10.*dist_unit;//10fm... i.e. almost start over
    }
    else{
      desired_fluct = bdist;
    }
    if(DEBUG) printf("\n A: Stuck=%u(%u); UnstuckCounter=%u; SinglePar=%u",Stuck,VeryStuck,UnstuckCounter,SinglePar);
//desired_fluct = 1;

    //if(CurrentPot=="Yukawa"){
    //  Kitty_SE.SetEpsilonConv(5e-8);
    //  Kitty_SE.SetEpsilonProp(5e-8);
    //  Kitty_SE.SetMaxRad(128);
    //  Kitty_SE.SetMaxRho(64);
    //}

    double MaxConv;
    if(CurrentPot=="YukawaDimiCore") MaxConv = 4;
    else MaxConv = 2;
    if(avg_fluct<0||Starting) Convergence[SinglePar]=1;
    else if(avg_fluct/desired_fluct>1) Convergence[SinglePar]/=1.25;
    else if(Convergence[SinglePar]<0.8) Convergence[SinglePar]*=1.25;
    //not good for the potentials to get too scattered values
    else if(Convergence[SinglePar]<MaxConv) Convergence[SinglePar]*=1.05;
    else Convergence[SinglePar] = MaxConv;

    //this will make min:max value Convergence[SinglePar]:1
    //will go towards one if we have no entries getting towards the solution
    //apply only in normal mode, not when already stuck.
    if(Convergence[SinglePar]<1&&SinglePar==0){
      Convergence[SinglePar] *= (1.+(1./Convergence[SinglePar]-1)/(1.+exp((avg_dir-0.25/2.)/(0.25/8.))));
    }



//ADD NUM ITER WITHOUT IMPROVEMENT TO INCREASE THE Convergence
    if(DEBUG){
      printf("\n fluct = %.2f(%.2f)",avg_fluct,desired_fluct);
      printf("\n avg_dir = %.2f(%.2f)",avg_dir,(1.+(1./Convergence[SinglePar]-1)/(1.+exp((avg_dir-0.25/2.)/(0.25/8.)))));
      printf("\n Convergence = %.3f (%.3f)",Convergence[SinglePar],Convergence[0]);
      printf(" %.2e; %.2f; %.2e %.2e\n",bV_1,bmu_1,bV_2,bmu_2);
      //printf("\n bdf_0=%.2f; bdd_0=%.2f",bdf_0,bdd_0);
      //printf("\n bf_0=%.2f; bd_0=%.2f",bf_0,bd_0);
    }

    if(Starting&&!Started){
      if(!StartPars) OkayishStartingPars(CurrentPot,f0,d0,bV_1,bmu_1,bV_2,bmu_2,s_2);
      else {bV_1=StartPars[0];bmu_1=StartPars[1];bV_2=StartPars[2];bmu_2=StartPars[3];}
      V_1 = bV_1;
      mu_1 = bmu_1;
      V_2 = bV_2;
      mu_2 = bmu_2;
      Started = true;
      //printf("STARTING: V1=%.f\n",V_1);
    }
    else{
      if(CurrentPot=="YukawaDimiCore"){
        if(!SinglePar||SinglePar==1){
          do V_1 = rangen.Gaus(bV_1,0.005+5.*Convergence[1]*Convergence[0]*FineTune);
          while(V_1<0);
        }
        if(!SinglePar||SinglePar==2){
          do{mu_1 = rangen.Gaus(bmu_1,0.001+1.*Convergence[2]*Convergence[0]*FineTune);}
          while(mu_1<0.||mu_1>3.0);
        }
        if(!SinglePar||SinglePar==3){
          V_2 = rangen.Gaus(bV_2,0.2+200.*Convergence[3]*Convergence[0]*FineTune);
          //do V_2 = rangen.Gaus(bV_2,0.2+200.*Convergence[3]*Convergence[0]);
          //while(V_2<0);
        }
        if(!SinglePar||SinglePar==4){
          do{mu_2 = rangen.Gaus(bmu_2,0.0005+0.5*Convergence[4]*Convergence[0]*FineTune);}
          while(mu_2<0.1||mu_2>1.5);
        }
      }
      else if(CurrentPot=="Gaussian"){
        if(!SinglePar||SinglePar==1){
          V_1 = rangen.Gaus(bV_1,0.05+40.*Convergence[1]*Convergence[0]*FineTune);
        }
        if(!SinglePar||SinglePar==2){
          do{mu_1 = rangen.Gaus(bmu_1,0.001+1.0*Convergence[2]*Convergence[0]*FineTune);}//----
          while(mu_1<0.||mu_1>3.0);
        }
        V_2 = 0;
        mu_2 = 0;
      }
      else if(CurrentPot=="Yukawa"){
        if(!SinglePar||SinglePar==1){
          do V_1 = rangen.Gaus(bV_1,0.05+20.*Convergence[1]*Convergence[0]*FineTune);
          while(V_1<0);
        }
        if(!SinglePar||SinglePar==2){
          do{mu_1 = rangen.Gaus(bmu_1,0.001+1.0*Convergence[2]*Convergence[0]*FineTune);}
          while(mu_1<0.||mu_1>3.0);
        }
        V_2 = 0;
        mu_2 = 0;
      }
      else{
        //printf("V1=%.f\n",V_1);
        //this is the long range guy, we expect it to be negative for positive f0
        //do {
        if(!SinglePar||SinglePar==1){
          V_1 = rangen.Gaus(bV_1,0.2+500.*Convergence[1]*Convergence[0]*FineTune);
        }
        //}
        //while(f0>0&&V_1>0);
        if(!SinglePar||SinglePar==2){
          do{mu_1 = rangen.Gaus(bmu_1,0.001+1.*Convergence[2]*Convergence[0]*FineTune);}
          while(mu_1<0.||mu_1>2.5);
        }
        //the second gaussian should be short ranged and stronger in amplitude
        //do{
        if(!SinglePar||SinglePar==3){
          V_2 = rangen.Gaus(bV_2,0.2+500.*Convergence[3]*Convergence[0]*FineTune);
        }
        //}
        //while(fabs(V_2)<fabs(V_1));
        if(!SinglePar||SinglePar==4){
          do{mu_2 = rangen.Gaus(bmu_2,0.001+1.*Convergence[4]*Convergence[0]*FineTune);}
          while(mu_2<0.||mu_2>2.5);
        }
        //printf(" ---> V1=%.f\n",V_1);
        //usleep(100e3);
      }
    }


    TH1F* hDummy; TF1* fDummy;
    Kitty_SE.SetShortRangePotential(0,0,0,V_1);
    Kitty_SE.SetShortRangePotential(0,0,1,mu_1);
    Kitty_SE.SetShortRangePotential(0,0,2,V_2);
    Kitty_SE.SetShortRangePotential(0,0,3,mu_2);
    if(CurrentPot=="YukawaDimiCore"){
      if(StartPars&&StartPars[4]){
        s_2 = StartPars[4];
      }
      else s_2 = mu_2*0.2;
    }
    else if(CurrentPot=="Gaussian") s_2 = -1;
    else if(CurrentPot=="Yukawa") s_2 = -2;
    else s_2 = 0;
    Kitty_SE.SetShortRangePotential(0,0,4,s_2);
    //Kitty_SE.SetNotifications(CATS::nAll);
    //printf(" KillTheCat\n");
    Kitty_SE.SetAnaSource(0,Radii[0]);
    Kitty_SE.KillTheCat();

    if(DEBUG) printf("\n    V1=%.1f mu1=%.3f V2=%.1f mu2=%.3f",V_1,mu_1,V_2,mu_2);
    if(DEBUG) printf("\n b: V1=%.1f mu1=%.3f V2=%.1f mu2=%.3f",bV_1,bmu_1,bV_2,bmu_2);
//printf("kfgnslkjglfkgnjdkf\n");
    if(!DEBUG){
      printf("\n Current solution: V1=%.4e  mu1=%.4e  V2=%.4e  mu2=%.4e <--> f0=%.3f  d0=%.3f      ",
      bV_1,bmu_1,bV_2,bmu_2,bf_0,bd_0);
    }

    if(Eval_ScattParameters(Kitty_SE,f_0,d_0,hDummy,fDummy)) BadPhaseShifts = 0;
    else BadPhaseShifts++;
    //hDummy->Write();
    //fDummy->Write();

    if(hDummy) delete hDummy; if(fDummy) delete fDummy;
    if(BadPhaseShifts) {
      //cout<<flush;cout<<"\033[F";cout<<flush;cout<<"\033[F";
      //cout<<flush;cout<<"\033[F";cout<<flush;cout<<"\033[F";
      if(!DEBUG){
        //cout << "JULI" <<endl;
        cout << flush;
        cout << "\e[A";
        cout << flush;
      }
      else printf("\n ISSUE");
      continue;
    }
    df_0 = fabs(f_0-f0);
    dd_0 = fabs(d_0-d0);
    dist = sqrt(0.5*df_0*df_0/df0/df0+0.5*dd_0*dd_0/dd0/dd0);
    if(DEBUG) printf("\n    f_0=%.2f df_0=%.2e; d_0=%.2f dd_0=%.2e",f_0,df_0,d_0,dd_0);
    if(DEBUG) printf("\n b: f_0=%.2f df_0=%.2e; d_0=%.2f dd_0=%.2e",bf_0,bdf_0,bd_0,bdd_0);

    //we get a better result, that has the correct sign
    //GoodGoing = ( (df_0<bdf_0||df_0<df0)&&(dd_0<bdd_0||dd_0<dd0)
    //              &&(f_0*bf_0>0)&&(f_0*f0>0)&&(d_0*bd_0>0)&&(d_0*d0>0));
    //we reduce df_0 or dd_0 (obsolete?)(df_0<bdf_0||dd_0<bdd_0)&&
    //we have the correct sign
    //we reduce the distance
    //GoodGoing = (dist<bdist&&(f_0*bf_0>0)&&(f_0*f0>0)&&(d_0*bd_0>0)&&(d_0*d0>0));
    //the above line files for negative effective range
    GoodGoing = (dist<bdist&&(f_0*bf_0>0)&&(f_0*f0>0));
    //if(bdist<1e15) fluct[ufluct] = fabs(dist-bdist);
    //if(bdist<1e15){
//CHANGED TO BEST, SHOULD I???
      //fluct[ufluct] = sqrt(0.5*fabs(f_0-lf_0)*fabs(f_0-lf_0)/df0/df0+0.5*fabs(d_0-ld_0)*fabs(d_0-ld_0)/dd0/dd0);
      fluct[ufluct] = sqrt(0.5*fabs(f_0-bf_0)*fabs(f_0-bf_0)/df0/df0+0.5*fabs(d_0-bd_0)*fabs(d_0-bd_0)/dd0/dd0);
      //imagine 3 pts, B,N,G (best, new, goal), spanning vectors from BG=g and BN=n
      //than the vectors are given as (nx-bx,ny-by) etc., where x in my case is f0 and y is d0
      //than we compute the angle from the dot product rule. We want it |alpha|<45 deg to move in the right direction
      double nx=(f_0-bf_0); double ny=(d_0-bd_0);
      double gx=(f0-bf_0); double gy=(d0-bd_0);
      double CosAlpha;
      if((nx==0&&ny==0)||(gx==0&&gy==0)) CosAlpha = 1;
      else CosAlpha = (nx*gx+ny*gy)/(sqrt(nx*nx+ny*ny)*sqrt(gx*gx+gy*gy));
      //so true, if we fluctuated towards the soltion, without going more than dist in the other direction
//perhaps stuck mode, where we increase fluct, and change the best solution as long as we go to "the other side"
      dir[udir] = (CosAlpha>(sqrt(2.)/2.))&&(fluct[ufluct]<2.*bdist);

      if(DEBUG){
        printf("\n Current iter: GG=%i; CosAlpha=%.2f",GoodGoing,CosAlpha);
      }
    //}
    //else{
    //  fluct[ufluct] = dist;
    //}
    ufluct++;ufluct=ufluct%fluctN;
    udir++;udir=udir%dirN;

    if(GoodGoing){
      bf_0 = f_0;
      bdf_0 = df_0;
      bd_0 = d_0;
      bdd_0 = dd_0;
      bV_1 = V_1;
      bmu_1 = mu_1;
      bV_2 = V_2;
      bmu_2 = mu_2;
      //Distance = sqrt(bf_0*bf_0+bd_0*bd_0);
      bdist = sqrt(0.5*df_0*df_0/df0/df0+0.5*dd_0*dd_0/dd0/dd0);
      Stuck = 0;
      Starting = false;
    }
    lf_0 = f_0;
    ld_0 = d_0;

    FoundIt = (fabs(df0)>=df_0&&fabs(dd0)>=dd_0);
    for(unsigned uRad=0; uRad<NumRad; uRad++){
      if(uRad){
        Kitty_SE.SetAnaSource(0,Radii[uRad]);
        Kitty_SE.KillTheCat();
      }
      int Class = Get_MMclass(Radii[uRad],f_0,d_0);
      Fill_ntMM(ntMM,Kitty_SE,Radii[uRad],f_0,d_0,V_1,mu_1,V_2,mu_2,s_2,10+FoundIt);
      hClassCount->Fill(Class);
      if(uRad==0) hif0d0Count->Fill(1./f_0,d_0);
    }

    if(SinglePar){
      UnstuckCounter++;
    }
    else{
      uIter++;
    }
    if(bdf_0<0.999e6||bdf_0>1.001e6){
      Stuck++;
    }
    if(DEBUG){
      printf("\n FoundIt=(df0>=df_0&&dd0>=dd_0); %i=(%.3f>=%.3f && %.3f>=%.3f)\n-----------------------------",FoundIt,fabs(df0),df_0,fabs(dd0),dd_0);
      usleep(500e3);
    }
//printf("hi\n");
    //cout<<flush;cout<<"\033[F";cout<<flush;cout<<"\033[F";
    //cout<<flush;cout<<"\033[F";cout<<flush;cout<<"\033[F";
    //cout<<flush;cout<<"\033[F";
    if(!DEBUG){
      cout << flush;
      cout << "\e[A";
      cout << flush;
    }

  }
  if(DEBUG){printf("\ndist: %.2f vs %.2f\n",bdist,FallBack[0]);}
  if(bdist>FallBack[0]){
    bdist = FallBack[0];
    bf_0 = FallBack[1];
    bd_0 = FallBack[2];
    bV_1 = FallBack[3];
    bmu_1 = FallBack[4];
    bV_2 = FallBack[5];
    bmu_2 = FallBack[6];
    s_2 = FallBack[7];
    if(fabs(s_2)<0.001) CurrentPot = "DoubleGaussSum";
    else CurrentPot = "YukawaDimiCore";
  }
  //TH1F* hClassCount = new TH1F("hClassCount","hClassCount",100,-50.5,49.5);
  //TH2F* hif0d0Count = new TH2F("hf0d0Count","hf0d0Count",20,-10,10,20,0,40);
  fOutput.cd();
  ntMM->Write();
  hClassCount->Write();
  hif0d0Count->Write();
  for(unsigned uRad=0; uRad<NumRad; uRad++){
    Kitty_SE.SetAnaSource(0,Radii[uRad]);
    Kitty_SE.KillTheCat();
    //printf("\n Cross check the scattering parameters");
//BINNING_ISSUE->JUST LOOK AT THE KCOTG IN BOTH PLOTS
    //double cc_f0=bf_0; double cc_d0=bd_0; TH1F* hCC; TF1* fCC;
//THE BINNING COMING OUT OF HERE
    //bool statusCC = Eval_ScattParameters(Kitty_SE,cc_f0,cc_d0,hCC,fCC);
    //printf("\n status=%i; f0=%.3f; d0=%.3f",statusCC,cc_f0,cc_d0);
    //hCC->Write();
    //fCC->Write();
    //delete hCC; delete fCC;
//IS DIFFERENT THAN OUT OF THERE
    Write_MM_Plots(&fOutput,Kitty_SE,Radii[uRad],bf_0,bd_0,TString::Format("_r%id%i",int(Radii[uRad]*100.)/100,int(Radii[uRad]*100.)%100));
  }

  printf("\n");
  printf("Suitable %s potential found:\n",CurrentPot.Data());
  printf(" f0 = %.3f fm\n", bf_0);
  printf(" d0 = %.2f fm\n", bd_0);
  printf("  V1 = %.6e\n",bV_1);
  printf("  mu1 = %.6e\n",bmu_1);
  printf("  V2 = %.6e\n",bV_2);
  printf("  mu2 = %.6e\n",bmu_2);

  V1 = bV_1;
  V2 = bV_2;
  mu1 = bmu_1;
  mu2 = bmu_2;

  delete hClassCount;
  delete hif0d0Count;
  delete ntMM;
  delete [] MomBins;
}





//Potential==Dynamic, means we try to do the double Gaussian, if it fails we repeat it all with Yukawa
void ManufacturePotential2(const double f0, const double df0,
                                const double d0, const double dd0, const double* Radii, const unsigned NumRad,
                                double& V1, double& mu1, double& V2, double& mu2,
                                const TString Potential, const TString OutputFolder, const int& SEED=11){

  const double RedMass = 488.6;//pphi
  const bool DEBUG = false;
  //the scale to which 1fm roughly correponds
  //const double dist_unit = df0*df0+dd0*dd0;
  const unsigned MaxIter = 8192;

  //b = best, l = last, d = difference to desired
  TString CurrentPot = Potential;
  if(Potential=="Dynamic") CurrentPot = "DoubleGaussSum";

  bool Starting = true;
  bool Started = false;
  const unsigned NumSavedPars = 8;
  //dist,f0,d0,V1,mu1,V2,mu2,s2:
  enum ManufacturedSolution { ms_dist, ms_f0, ms_d0, ms_V1, ms_mu1, ms_V2, ms_mu2, ms_s2 };
  //the globally best solution
  double* GBS = new double [NumSavedPars];
  for(unsigned usp=0; usp<NumSavedPars; usp++){GBS[usp]=1e16;}
  //the best solution of the current attempt
  double* CAS = new double [NumSavedPars];
  for(unsigned usp=0; usp<NumSavedPars; usp++){CAS[usp]=1e16;}
  //the current set of parameters
  double* CSP = new double [NumSavedPars];
  for(unsigned usp=0; usp<NumSavedPars; usp++){CSP[usp]=1e16;}
  //position of last calibration
  double* PLC = new double [NumSavedPars];
  for(unsigned usp=0; usp<NumSavedPars; usp++){PLC[usp]=1e16;}

  unsigned uIter=0;
  TRandom3 rangen(SEED);
  CATSparameters sPars(CATSparameters::tSource,1,true);
  sPars.SetParameter(0,2.5);
  const double kStepFine=3;
  const double kStepCoarse=10;
  const unsigned nMomFine = 40;
  const unsigned nMomCoarse = 12;
  const unsigned nMom = nMomFine+nMomCoarse;
  double* MomBins = new double [nMom+1];
  for(unsigned uMom=0; uMom<nMomFine; uMom++){
    MomBins[uMom] = kStepFine*double(uMom);
  }
  for(unsigned uMom=0; uMom<=nMomCoarse; uMom++){
    MomBins[nMomFine+uMom] = kStepFine*double(nMomFine)+kStepCoarse*double(uMom);
  }
  CATS Kitty_SE;
  Kitty_SE.SetMomBins(nMom,MomBins);
  Kitty_SE.SetAnaSource(GaussSource, sPars);
  Kitty_SE.SetUseAnalyticSource(true);
  Kitty_SE.SetQ1Q2(0);
  Kitty_SE.SetQuantumStatistics(false);
  //Kitty_SE.SetRedMass(Mass_L*0.5);
  Kitty_SE.SetRedMass(RedMass);//pphi
  Kitty_SE.SetNumChannels(1);
  Kitty_SE.SetNumPW(0,1);
  Kitty_SE.SetSpin(0,0);
  Kitty_SE.SetChannelWeight(0, 1.);
  Kitty_SE.SetEpsilonConv(5e-9);
  Kitty_SE.SetEpsilonProp(5e-9);
  Kitty_SE.SetMaxRad(96);
  Kitty_SE.SetMaxRho(32);
  Kitty_SE.SetNotifications(CATS::nSilent);
  Kitty_SE.SetGridEpsilon(1./512.);
  double RadForGrid;
  Kitty_SE.KillTheCat();

  CATSparameters pPars(CATSparameters::tPotential,5,true);
  if(CurrentPot=="YukawaDimiCore") Kitty_SE.SetShortRangePotential(0,0,YukawaDimiCore,pPars);
  else if(CurrentPot=="Gassian") Kitty_SE.SetShortRangePotential(0,0,SingleGauss,pPars);
  else if(CurrentPot=="Yukawa") Kitty_SE.SetShortRangePotential(0,0,YukawaDimiSmooth,pPars);
  else Kitty_SE.SetShortRangePotential(0,0,DoubleGaussSum,pPars);

  bool GoodGoing=false;
  double dist=1e16;
  double bdist=1e16;
  //Fluctuations: the average distance between our working solution and new solution
  //if too large -> decrease Convergence
  //if too small -> increse Convergence
  //criteria for small/large: ditance between working solution and goal + the sphere
  //formed by the opening angle crossing a line perpendicular to the above distance and going through the goal
  const unsigned fluctN = 32;
  const unsigned fluctNth = 8;
  double fluct[fluctN];
  unsigned ufluct = 0;
  double avg_fluct;

  const unsigned TanN = 32;
  const unsigned TanNth = 8;
  //how many of the last dirN iterations moved towards the true solution, based on the opening angle
  //between the working solution and the goal. The opening angle will be adjusted based on the
  //expectation value based on the angle (assuming a flat anglular disto), and the value
  //measured one. The ratio of the two determines if we will increase or decrease the opening angle
  //the measured one is beeing reset on finding a new working solution, and we will NOT change the default angle
  //between we reach a threshold
  double Tan[TanN];
  unsigned uTan = 0;
  double TanInRange;
  double Exp_TanInRange;

  //the opening angle and its min/max values. Notice, that if we find a new working solution, we start
  //from the default value, which is tuned such that it is not possible to get the next solution further away than the
  //current one. If TanAlpha>1, we can actually end up with a soltion that goes away from our goal (used to avoid beeing stuck)
  //if we cross TanAlphaAdjust (going up), we make a certain number of iterations where we are adjusting the Convergence of the
  //single parameters. During that time, TanAlpha is not allowed to change
  const double TanAlphaDefault = 1.;
  const double TanAlphaStep = TanAlphaDefault/8.;//0.125
  const double TanAlphaMin = TanAlphaStep;//0.125
  const double TanAlphaAdjust = TanAlphaDefault*2.;
  const double TanAlphaMax = TanAlphaDefault*8.;
  double TanAlpha;

  //this is the procedure of fine tuning Convergence for each parameter. It is executed at the start,
  //and later on each time some of the parameters differs either by sign or by more than 2x its value at the last calibration
  //The calibration will be interrupted only if we find a new solution that lies within 1./TanAlphaAdjust or reach our goal.
  //the iterations during calibration do NOT count towards the stuck count.
  bool NeedToCalibrate = true;
  const unsigned SingleParN = 16;
  unsigned uSingleParN = 0;

  const unsigned Npars = 4;
  double Convergence[Npars+1];
  for(unsigned uPar=0; uPar<Npars+1; uPar++)Convergence[uPar]=1;
  //Num Iter since last new working solution, and corresponding limit at which we give up
  unsigned WorkStuck = 0;
  const unsigned MaxWorkStuck = 128;
  //Num Iter since last new best solution
  unsigned GlobalStuck = 0;
  const unsigned MaxGlobalStuck = 1024;

  //Num Iter with consecutive error in the PS calculation
  //if it reacheas a threshold, the TanAlpha is artificially increased by 2x.
  //these iterations DO count towards the stuck
  unsigned BadPS = 0;


  //for(unsigned uf=0; uf<fluctN; uf++) fluct[uf]=0;
  //for(unsigned ud=0; ud<dirN; ud++) dir[ud]=0;
//TFile fDump(TString::Format("%s/OtherTasks/ManufactureYukawaPotential/fDump.root",GetFemtoOutputFolder()),"recreate");

  TFile fOutput(TString::Format("%s/fOut_%.2f_%.1f_%s.root",OutputFolder.Data(),f0,d0,Potential.Data()),"recreate");
  TNtuple* ntMM = new TNtuple("ntMM", "ntMM","r0:f0:d0:MMA:MMB:C_X0:C_X150:L_X0:V1:mu1:V2:mu2:s2:wf2:Class:Manufactured");
  TH1F* hClassCount = new TH1F("hClassCount","hClassCount",100,-50.5,49.5);
  TH2F* hif0d0Count = new TH2F("hf0d0Count","hf0d0Count",20,-10,10,20,0,40);
  bool FoundIt = false;
  bool GiveUp = false;
  bool NewGlobSol = false;
  bool NewWorkSol = false;
/*
  while(!GiveUp&&!FoundIt){
    //PRINT

//after we compute the new guys:
//* if(BadPS)...
//* else:
//* check what is the opening angle
//* check what is the distance from start pt (fluct) and from goal
//* based on those, decide on NewWorkSol,NewGlobSol,FoundIt
//*   -> above apply special rules if we are currently calibrating
//* compute the Convergence for the next step
//* check if we need to calibrate




    if(!NewWorkSol&&!NewGlobSol){
      //if appropriate, evaluate the current avereges
//CHECK ALL OF THIS !!!!!

      //for the distance, we compute the avg distance we moved
      //away from our current working solution. The avg of that
      //is the quantity we need
      avg_fluct = 0;
      if(ufluct>=fluctNth){
        for(unsigned uf=0; uf<ufluct; uf++){
          avg_fluct+=fluct[uf];
        }
        avg_fluct /= double(ufluct);
        fluct[ufluct] = avg_fluct;
        ufluct++;
        if(ufluct>fluctN) ufluct=fluctN;
      }

//CHECK ALL OF THIS !!!!!
      //avg_fluct = 0;
      TanAlpha = TanAlphaDefault;
      if(uTan>=TanNth){
        for(TanAlpha=TanAlphaStep; TanAlpha<=TanAlphaMax; TanAlpha+=TanAlphaStep){
          double Alpha = arctan(TanAlpha);
          Exp_TanInRange = Alpha/(2.*Pi);
          TanInRange=0;
          for(unsigned ut=0; ut<uTan; ut++){
            TanInRange += (fabs(Tan[uTan])<TanAlpha)
          }
          TanInRange /= double(uTan);
          if(TanInRange/Exp_TanInRange>=1) break;
        }
        //uTan++;
        //if(uTan>fluctN) ufluct=fluctN;
      }

      const unsigned TanN = 32;
      const unsigned TanNth = 8;
      //how many of the last dirN iterations moved towards the true solution, based on the opening angle
      //between the working solution and the goal. The opening angle will be adjusted based on the
      //expectation value based on the angle (assuming a flat anglular disto), and the value
      //measured one. The ratio of the two determines if we will increase or decrease the opening angle
      //the measured one is beeing reset on finding a new working solution, and we will NOT change the default angle
      //between we reach a threshold
      double Tan[TanN];
      unsigned uTan = 0;
      double TanInRange;
      double Exp_TanInRange;

    }
    //if we find a new solution, working or global, the avg are reset
    else{
      avg_fluct = 0;
      ufluct = 0;


    }

    if(!NeedToCalibrate){
      uIter++;
      if(NewWorkSol) WorkStuck=0;
      else WorkStuck++;
      if(NewGlobSol) GlobalStuck=0;
      else GlobalStuck++;
    }
    FoundIt = (fabs(CSP[ms_f0]-f0)<fabs(df0) && fabs(CSP[ms_d0]-d0)<fabs(dd0));
    GiveUp = (WorkStuck>=MaxWorkStuck||GlobalStuck>=MaxGlobalStuck||uIter>MaxIter);
  }


  printf("\n");
  printf("Suitable %s potential found:\n",CurrentPot.Data());
  printf(" f0 = %.3f fm\n", bf_0);
  printf(" d0 = %.2f fm\n", bd_0);
  printf("  V1 = %.6e\n",bV_1);
  printf("  mu1 = %.6e\n",bmu_1);
  printf("  V2 = %.6e\n",bV_2);
  printf("  mu2 = %.6e\n",bmu_2);

  V1 = bV_1;
  V2 = bV_2;
  mu1 = bmu_1;
  mu2 = bmu_2;

  delete hClassCount;
  delete hif0d0Count;
  delete ntMM;
  delete [] MomBins;
  delete [] GBS;
  delete [] CAS;
  delete [] CSP;
  delete [] PLC;
*/
}







//!!!!!!!!! I have now introduced the mass as input, but I have set it only for Bhawani
//each time you rerun, CHANGE the masses to the value you desire!
void MakePotentials(int flag){
  printf("Make a potential %i\n",flag);
  double V1,V2,mu1,mu2;
  //ManufactureYukawaPotential(0.5,0.0025,1,0.005,V1,mu1,V2,mu2);
  //ManufactureYukawaPotential(0.5,0.0025,2,0.005*2,V1,mu1,V2,mu2);
  //ManufactureYukawaPotential(0.5,0.0025,4,0.005*4,V1,mu1,V2,mu2);
  //ManufactureYukawaPotential(0.5,0.0025,8,0.005*8,V1,mu1,V2,mu2);

//make the dist and fluct refer to d0 and f0 separately, check if the parametes get closer
//to their goal (not only distance, but direction...)

  unsigned NumR=16;
  double Radii[NumR];

  //NumR=8;
  //Radii[0] = 0.9;Radii[1] = 1.1;Radii[2] = 1.3;Radii[3] = 1.5;
  //Radii[4] = 2.0;Radii[5] = 2.5;Radii[6] = 3.0;Radii[7] = 4.0;

  NumR=1;
  Radii[0]=1.08;Radii[1]=1.08*2;Radii[2]=1.08*4;

  double f0,d0;
  double ef0 = 0.005;
  double ed0 = 0.005;
//NumR=1; f0=25.00; d0=25.00; ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,"Dynamic");
//return;

  //ADD-ONS FOR EACH FLAG
  TString OutputFolder = TString::Format("%s/MMM/ManufacturePotential/",GetCernBoxDimi());
  //flag = last digit is for the potential 0 - dynamic, 1 - DoubleGaussSum, 2 - YukawaDimiCore, 3 - Gaussian, 4 - Yukawa
  TString Potential = "Dynamic";
  if(flag%10==1) Potential = "DoubleGaussSum";
  else if(flag%10==2) Potential = "YukawaDimiCore";
  else if(flag%10==3) Potential = "Gaussian";
  else if(flag%10==4) Potential = "Yukawa";
  if(flag/10==1){
    OutputFolder += "";
    f0=0.25; d0=0.50; ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder);
    f0=0.50; d0=0.50; ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder);
    f0=1.00; d0=0.50; ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder);
    f0=2.00; d0=0.50; ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder);
    f0=4.00; d0=0.50; ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder);
    f0=8.00; d0=0.50; ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder);
  }
  else if(flag/10==2){
    OutputFolder += "";
    f0=0.25; d0=1.00; ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder);
    f0=0.50; d0=1.00; ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder);
    f0=1.00; d0=1.00; ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder);
    f0=2.00; d0=1.00; ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder);
    f0=4.00; d0=1.00; ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder);
    f0=8.00; d0=1.00; ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder);
  }
  else if(flag/10==3){
    OutputFolder += "";
    f0=0.25; d0=2.00; ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder);
    f0=0.50; d0=2.00; ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder);
    f0=1.00; d0=2.00; ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder);
    f0=2.00; d0=2.00; ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder);
    f0=4.00; d0=2.00; ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder);
    f0=8.00; d0=2.00; ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder);
  }
  else if(flag/10==4){
    OutputFolder += "";
    f0=0.25; d0=4.00; ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder);
    f0=0.50; d0=4.00; ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder);
    f0=1.00; d0=4.00; ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder);
    f0=2.00; d0=4.00; ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder);
    f0=4.00; d0=4.00; ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder);
    f0=8.00; d0=4.00; ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder);
  }
  else if(flag/10==5){
    OutputFolder += "";
    f0=0.25; d0=8.00; ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder);
    f0=0.50; d0=8.00; ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder);
    f0=1.00; d0=8.00; ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder);
    f0=2.00; d0=8.00; ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder);
    f0=4.00; d0=8.00; ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder);
    f0=8.00; d0=8.00; ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder);
  }
  else if(flag/10==6){
    OutputFolder += "";
    f0=0.25; d0=16.00; ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder);
    f0=0.50; d0=16.00; ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder);
    f0=1.00; d0=16.00; ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder);
    f0=4.00; d0=16.00; ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder);
    f0=2.00; d0=16.00; ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder);
    f0=8.00; d0=16.00; ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder);
  }
  else if(flag/10==7){
    OutputFolder += "";
    f0=0.25; d0=0.125; ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder);
    f0=0.50; d0=0.250; ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder);
    f0=32.0; d0=16.00; ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder);
    f0=16.0; d0=8.000; ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder);
  }
  //for emma 20X, 21X....
  else if(flag/10>=20&&flag/10<=27){
    //ef0 = 0.01;
    //ed0 = 0.01;
    NumR = 1;
    Radii[0] = 1.08;
    OutputFolder += "Emma2/";
    int VAR_FLAG = (flag/10)%10;
    f0 = 0.45+0.1*double(VAR_FLAG);
    for(unsigned ud0=0; ud0<=5; ud0++){
      d0 = 3.00+3.0*double(ud0);
      ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder);
    }
  }
  //for ramona 40X, 41X....
  else if(flag/10>=40&&flag/10<=44){
    ef0 = 0.01;
    ed0 = 0.01;
    NumR = 2;
    Radii[0] = 0.95;
    Radii[1] = 1.30;
    OutputFolder += "Ramona1/";
    int VAR_FLAG = (flag/10)%10;
    f0 = -0.3-0.05*double(VAR_FLAG);
    for(int ud0=3; ud0>=0; ud0--){
      d0 = 5.00-1.5*double(ud0);
      ManufacturePotential(f0,fabs(f0*ef0)+0.01,d0,fabs(d0*ed0)+0.01,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder);
    }
  }
  //for Prof. S 50X, 51X....
  // 50X,51X - doublet potentials
  // 52X,53,54X - quarted potentials, if other -> perform all of them
  // the 1XXX are for systematics
  // e.g. for the flag 50X, systematics are 150X, 160X, 170X, 180X,
  // for flag 51X we will have 151X, 161X, 171X, 181X etc.
  else if( (flag/10>=50&&flag/10<=59)||(flag/10>=150&&flag/10<=199) ){
    ef0 = 0.1/10.;
    ed0 = 0.05/5.;
    NumR = 1;
    Radii[0] = 1.2;
    OutputFolder += "ProfS/";
    int VAR_FLAG = (flag/10)%10;
    if(flag/10>=150){
      VAR_FLAG += (flag/100)*10;//+150,160,170,180
    }
    //for(unsigned uVar=0; uVar<5; uVar++){
      //if(VAR_FLAG<5&&uVar!=VAR_FLAG) continue;
      switch (VAR_FLAG) {
        case 0: f0=-16.8; d0 = 2.3; break;
        case 150: f0=-21.2; d0 = 2.3; break;
        case 160: f0=-14.4; d0 = 2.3; break;
        case 170: f0=-16.8; d0 = 2.6; break;
        case 180: f0=-16.8; d0 = 2.0; break;
        case 1: f0=-16.3; d0 = 3.2; break;
        case 2: f0=7.6; d0 = 3.6; break;
        case 3: f0=10.8; d0 = 3.8; break;
        case 4: f0=17.3; d0 = 3.6; break;
        default: printf("Weird flags for producing the potentials for Prof. S.\n"); return;
      }
      double StartPars[5];
      double FT = 1;

      if(Potential=="DoubleGaussSum"){
        //this does not work, ones you have result for uVar==1, plug in those here
        if(VAR_FLAG==0){
          //StartPars[0]=-8.999757e+02;
          //StartPars[1]=5.081559e-01;
          //StartPars[2]= -4.021392e+02;
          //StartPars[3]=1.019864e+00;
          //StartPars[4]=0;
          //StartPars[0]=-1.943168e+02;
          //StartPars[1]=1.372106e+00;
          //StartPars[2]= 3.646892e+02;
          //StartPars[3]=9.833930e-01;
          //StartPars[4]=0;
          StartPars[0]=-1.455808e+02;
          StartPars[1]=1.165081e+00;
          StartPars[2]=3.700546e+02;
          StartPars[3]=6.652901e-01;
          StartPars[4]=0;
        }
        else if(VAR_FLAG==150){
          StartPars[0]=-1.454633e+02;
          StartPars[1]=1.156062e+00;
          StartPars[2]=3.649164e+02;
          StartPars[3]=6.643602e-01;
          StartPars[4]=0;
          FT = 1./128.;
  //        Make a potential 1501
  //         Goal: DoubleGaussSum (f0,d0) [Achieved]: (-21.200,2.300)+/-(0.010,0.010) [ 38%,251%], Break (   97/10000) stuck, but working on it...
  //        Suitable DoubleGaussSum potential found:.1559e+00  V2=3.6454e+02  mu2=6.6436e-01 <--> f0=-21.174  d0=2.304
  //         f0 = -21.197 fm
  //         d0 = 2.30 fm
  //          V1 = -1.454633e+02
  //          mu1 = 1.156062e+00
  //          V2 = 3.649164e+02
  //          mu2 = 6.643602e-01


        }
        else if(VAR_FLAG==160){
          StartPars[0]=-1.406252e+02;
          StartPars[1]=1.183071e+00;
          StartPars[2]=3.671337e+02;
          StartPars[3]=6.646221e-01;
          StartPars[4]=0;
          FT = 1./1024.;
          // Current solution: V1=-1.4051e+02  mu1=1.1833e+00  V2=3.6563e+02  mu2=6.6525e-01 <--> f0=-14.397  d0=2.289
        }
        else if(VAR_FLAG==170){
          StartPars[0]=-1.442826e+02;
          StartPars[1]=1.251641e+00;
          StartPars[2]=3.755765e+02;
          StartPars[3]=7.397668e-01;
          StartPars[4]=0;
          FT = 1./256.;
          //ef0 = 1.0;

    //      Make a potential 1701
    //       Goal: DoubleGaussSum (f0,d0) [Achieved]: (-16.800,2.600)+/-(0.010,0.010) [ 15%,139%], Break (   89/10000) stuck, but working on it...
    //      Suitable DoubleGaussSum potential found:.2535e+00  V2=3.7566e+02  mu2=7.4164e-01 <--> f0=-16.867  d0=2.593
    //       f0 = -16.801 fm
    //       d0 = 2.59 fm
    //        V1 = -1.446716e+02
    //        mu1 = 1.253497e+00
    //        V2 = 3.756614e+02
    //        mu2 = 7.415382e-01

        }
        else if(VAR_FLAG==180){
          StartPars[0]=-1.420678e+02;
          StartPars[1]=1.081553e+00;
          StartPars[2]=3.579042e+02;
          StartPars[3]=5.812257e-01;
          StartPars[4]=0;
          FT = 1./2048.;
          //ef0 = 1.0;
        }
        else if(VAR_FLAG==1){
          //StartPars[0]=-8.999757e+02;
          //StartPars[1]=5.081559e-01;
          //StartPars[2]= -4.021392e+02;
          //StartPars[3]=1.019864e+00;
          //StartPars[4]=0;
          StartPars[0]=-1.707025e+02;
          StartPars[1]=1.495874e+00;
          StartPars[2]= 3.626808e+02;
          StartPars[3]=9.711152e-01;
          StartPars[4]=0;
        }
        else if(VAR_FLAG==2){
          StartPars[0]=-1.008636e+01;
          StartPars[1]=1.348757e+00;
          StartPars[2]= -1.410447e+01;
          StartPars[3]=2.147635e+00;
          StartPars[4]=0;
        }
        else if(VAR_FLAG==3){
          StartPars[0]=-4.867921e+02;
          StartPars[1]=1.175196e+00;
          StartPars[2]= -8.110269e+02;
          StartPars[3]=1.169148e-01;
          StartPars[4]=0;
        }
        else{
          StartPars[0]=-4.940004e+02;
          StartPars[1]=1.177133e+00;
          StartPars[2]= -8.133919e+02;
          StartPars[3]=1.007169e-01;
          StartPars[4]=0;
        }
        //if(uVar==1) ManufacturePotential(f0,ef0,d0,ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder,11,StartPars);
        //else ManufacturePotential(f0,ef0,d0,ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder);
        ManufacturePotential(f0,ef0,d0,ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder,11,StartPars,
          (Mass_L*Mass_d)/(Mass_L+Mass_d),FT);

      }
      else{
        ManufacturePotential(f0,ef0,d0,ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder,11,NULL,
          (Mass_L*Mass_d)/(Mass_L+Mass_d),FT);
      }


    //}
  }
  //for fopra S 60X.
  // pLambda fixed to NLO 600 scattering parameters. Effective single channel potential
  else if(flag/10>=60&&flag/10<=69){
    ef0 = 0.01;
    ed0 = 0.01;
    NumR = 1;
    Radii[0] = 1.38;
    OutputFolder += "FemtoFopra/";
    int VAR_FLAG = (flag/10)%10;
    for(unsigned uVar=0; uVar<4; uVar++){
      if(VAR_FLAG<5&&uVar!=VAR_FLAG) continue;
      switch (uVar) {
        case 0: f0=1.88; d0 = 2.74; break;
        case 1: f0=1.9; ef0=0.05; d0 = 2.75;  ed0=0.05; break;
        case 2: f0=1.88; ef0=0.02; d0 = 2.74;  ed0=0.03; break;
        case 3: f0=1.88; ef0=0.02; d0 = 2.75;  ed0=0.05; break;
        default: printf("Weird flags for producing the potentials for the femto fopra\n"); return;
      }
      double StartPars[5];
      if(VAR_FLAG==2){
        StartPars[0] = 2.76705e+00;
        StartPars[1] = 4.46223e-01;
        StartPars[2] = 1.72309e+03;
        StartPars[3] = 3.84851e-01;
        StartPars[4] = 0.1;
      }
      else if(VAR_FLAG==3){
        StartPars[0] = 5.04406e+00;
        StartPars[1] = 3.88958e-01;
        StartPars[2] = 1.83273e+03;
        StartPars[3] = 3.63198e-01;
        StartPars[4] = 1.03943e-01;
      }
      else OkayishStartingPars(Potential,f0,d0,StartPars[0],StartPars[1],StartPars[2],StartPars[3],StartPars[4]);
      ManufacturePotential(f0,ef0,d0,ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder,11,StartPars);
    }
  }
  //for Bhawani, p-d (7XX)
  //0 = Oers, Brockmann et al.(1967) (quartet)
  //1 = Oers, Brockmann et al.(1967) (doublet)
  //2 = Arvieux et al.(1973) (quartet)
  //3 = Arvieux et al.(1973) (doublet)
  //4 = Huttel et al.(1983) (quartet)
  //5 = Huttel et al.(1983) (doublet)
  //6 = Kievsky et al.(1997) (quartet)
  //7 = Kievsky et al.(1997) (doublet)
  //8 = Black et al. (1999) (quartet)
  //9 = Black et al. (1999) (doublet)
  //80X is something from the p-D meson I think, some check Bhawani needed
  else if(flag/10>=70&&flag/10<=89){
    NumR = 1;
    Radii[0] = 1.2;
    OutputFolder += "Bhawani/";
    int VAR_FLAG = (flag/10)%10;
    ef0 = VAR_FLAG%10==0?0.1:0.01;
    ed0 = 1.0;
    d0 = 1.2;
    double FineTune = 1./16.;
    if(flag/10<80){
      for(unsigned uVar=0; uVar<10; uVar++){
        if(VAR_FLAG<10&&uVar!=VAR_FLAG) continue;
        switch (uVar) {
          case 0: f0=11.4; break;
          case 1: f0=1.2; break;
          case 2: f0=11.88; break;
          case 3: f0=2.73; break;
          case 4: f0=11.1; break;
          case 5: f0=4.0; break;
          case 6: f0=13.8; break;
          case 7: f0=0.024; ef0=0.002; break;
          case 8: f0=14.7; break;
          case 9: f0=-0.13; break;
          default: printf("Weird flags for producing the potentials for Bhawani\n"); return;
        }
      }
    }
    else if(flag/10<90){
      for(unsigned uVar=0; uVar<10; uVar++){
        if(VAR_FLAG<10&&uVar!=VAR_FLAG) continue;
        switch (uVar) {
          //    Nucl.Phys.A 221 (1974) 253-268 (case 0 and 1)
          case 0: f0 = -2.73; d0 = 2.27; break;//S = 1/2
          case 1: f0 = -11.88; d0 = 2.63; break;//S = 3/2
          default: printf("Weird flags for producing the potentials for Bhawani\n"); return;
        }
      }
      FineTune = 1./4.;
      ef0 = 0.01*64*FineTune;//0.01 for fine tune = 64, at FT=8 this is 0.08
      ed0 = 0.02*64*FineTune;
    }
    else{
      //NEGATIVE EFFECTIVE RANGE FOR THIS MODEL, WE PUT IT ON HOLD
    }
    double StartPars[5];
    if(VAR_FLAG==0){
      //first good solution, with negative potential ???
      //it should be, as BS is here. It struggled to find good sol with positive pot pars
      //StartPars[0] = -1.548378e+03;
      //StartPars[1] = 5.455562e-01;
      //StartPars[2] = -1.053892e+03;
      //StartPars[3] = 1.560974e+00;
      //ed0 /= 4;
      StartPars[0] = -7.792753e+02;
      StartPars[1] = 1.014193e+00;
      StartPars[2] = 3.425267e+02;
      StartPars[3] = 2.014962e+00;
      ed0/=4;

      //FineTune = 1./11512.;
      //printf("hello\n");
    }
    else if(VAR_FLAG==1){
      //f0 = -11.885 fm
      //d0 = 2.54 fm
      // V1 = 1.613084e+02
      // mu1 = 1.631784e+00
      // V2 = -5.414474e+02
      // mu2 = 1.364351e+00
      StartPars[0] =1.613084e+02;
      StartPars[1] = 1.631784e+00;
      StartPars[2] = -5.414474e+02;
      StartPars[3] = 1.364351e+00;

      //f0 = -11.879 fm
      //d0 = 2.61 fm
      //V1 = 1.607812e+02
      //mu1 = 1.614273e+00
      //V2 = -5.372152e+02
      //mu2 = 1.362657e+00


    }
    //OkayishStartingPars(Potential,f0,d0,StartPars[0],StartPars[1],StartPars[2],StartPars[3],StartPars[4]);
    //ManufacturePotential(f0,ef0,d0,ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder,11,StartPars);
    ManufacturePotential(f0,ef0,d0,ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder,16,StartPars,
      (Mass_p*Mass_d)/(Mass_p+Mass_d),FineTune);
  }
  //901,911,921,... (for double gauss)
  else if(flag/10>=90&&flag/10<=99){
    NumR = 1;
    Radii[0] = 1.2;
    int VAR_FLAG = (flag/10)%10;
    if(VAR_FLAG>=7){
      OutputFolder += "Marcelo/";
    }
    else{
      OutputFolder += "Rafa/";
    }

    double StartPars[5];
    //OkayishStartingPars(Potential,f0,d0,StartPars[0],StartPars[1],StartPars[2],StartPars[3],StartPars[4]);
    //ManufacturePotential(f0,ef0,d0,ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder,11,StartPars);
    double FineTune = 1.;


    for(unsigned uVar=0; uVar<10; uVar++){
      if(VAR_FLAG<10&&uVar!=VAR_FLAG) continue;
      switch (uVar) {
////1) f0 = 0.1248 +/- 0.0011 MeV
//2) f0 = -0.1810 +/- 0.0085 MeV
        case 0: f0=0.1248; ef0=0.0011; d0=-0.22; ed0=0.02; break;
        case 1: f0=-0.1810; ef0=0.0085; d0=-0.22; ed0=0.02; break;
        //p pim (Philipp thesis)
        case 2: f0=0.12306; ef0=f0*0.0025; d0=11.56911; ed0=d0*0.0025; break;
        case 3: f0=0.12306; ef0=f0*0.0025; d0=0; ed0=0.025; break;
        //p pip (Philipp thesis)
        case 4: f0=0.11819; ef0=f0*0.0025; d0=1.90136; ed0=d0*0.0025; break;
        case 5: f0=0.11819; ef0=f0*0.0025; d0=0; ed0=0.025; break;
        case 6: f0=0.11819; ef0=f0*0.0035; d0=21.5; ed0=d0*0.0035; break;
        //Marcelo
        case 7: f0=-0.1246; ef0=f0*0.001; d0=0.1; ed0=0.4; break;
        default: printf("Weird flags for producing the potentials for Rafa\n"); return;
      }
    }

    if(VAR_FLAG==0||VAR_FLAG==3){
      //V0 = 203.1 mu0 = 0.6203 V1 = -1106 mu1 = 0.4156
      StartPars[0] = 200;
      StartPars[1] = 0.6;
      StartPars[2] = -1100;
      StartPars[3] = 0.41;
      FineTune = 1./256.;
    }
    else if(VAR_FLAG==2){
      StartPars[0] = -3.259084e+01;
      StartPars[1] = 1.071914e+00;
      StartPars[2] = -2.225066e+02;
      StartPars[3] = 9.515194e-02;
      FineTune = 1./512.;
    }
    else if(VAR_FLAG==4){
      // Current solution: V1=-4.4841e+02  mu1=1.8420e-01  V2=-1.5662e+03  mu2=2.1340e-01 <--> f0=0.062  d0=1.875
      // Current solution: V1=-4.4897e+02  mu1=9.5603e-02  V2=-1.5418e+03  mu2=2.3632e-01 <--> f0=0.072  d0=1.847
      //f0=0.15; ef0=0.05; d0=2.5; ed0=0.5;
      StartPars[0] = 5.680787e+01;
      StartPars[1] = 7.816589e-01;
      StartPars[2] = -3.371110e+02;
      StartPars[3] = 5.692803e-01;
      FineTune = 1./256.;
    }
    else if(VAR_FLAG==5){
      StartPars[0] = 2.041351e+02;
      StartPars[1] = 5.924833e-01;
      StartPars[2] = -1.105715e+03;
      StartPars[3] = 4.042254e-01;
      FineTune = 1./512.;
    }
    else if(VAR_FLAG==1){
      StartPars[0] = 750;
      StartPars[1] = 0.4;
      StartPars[2] = 650;
      StartPars[3] = 0.45;
      FineTune = 1./256.;
    }
    //V1=7.6505e+02  mu1=3.6658e-01  V2=6.8593e+02  mu2=3.3039e-01
    else if(VAR_FLAG==7){
      //StartPars[0] = 750;
      //StartPars[1] = 0.4;
      //StartPars[2] = 650;
      //StartPars[3] = 0.45;
      StartPars[0] = 7.6505e+02;
      StartPars[1] = 3.6658e-01;
      StartPars[2] = 6.8593e+02;
      StartPars[3] = 3.3039e-01;
      FineTune = 1./256.;
    }
    else if(VAR_FLAG==6){
      StartPars[0] = -1.306545e+01;
      StartPars[1] = 1.455581e+00;
      StartPars[2] = -2.321616e+01;
      StartPars[3] = 1.100759e-01;
      FineTune = 1./256.;
    }
    StartPars[4] = 0;

    ManufacturePotential(f0,ef0,d0,ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder,11,StartPars,
      (Mass_p*Mass_pic)/(Mass_p+Mass_pic),FineTune);
  }



  else if(flag/10>=110&&flag/10<=119){
    NumR = 1;
    Radii[0] = 1.2;
    OutputFolder += "EmmaDani/";
    int VAR_FLAG = (flag/10)%10;
    double FT=1;
    double StartPars[5];
    StartPars[0]=1;
    StartPars[1]=1;
    StartPars[2]=-1;
    StartPars[3]=1;
    //if(flag/10<80){
    double RedMass = 0;
    for(unsigned uVar=0; uVar<10; uVar++){
      if(VAR_FLAG<10&&uVar!=VAR_FLAG) continue;
      switch (uVar) {
        case 0: f0=-4.9e-2; ef0=-0.125*3e-3; d0=7.1644; ed0=0.125*0.164; FT=0.125/64.;
                RedMass = (Mass_Kch*Mass_Dch)/(Mass_Kch+Mass_Dch);
                // Current solution: V1=-5.0966e+01  mu1=8.9633e-01  V2=7.8366e+01  mu2=8.1096e-01 <--> f0=-0.042  d0=7.141
                //Current solution: V1=-5.1056e+01  mu1=8.9607e-01  V2=7.9119e+01  mu2=8.0913e-01 <--> f0=-0.042  d0=7.140
                //Current solution: V1=-5.1048e+01  mu1=8.9693e-01  V2=7.9520e+01  mu2=8.0884e-01 <--> f0=-0.043  d0=7.132
                StartPars[0]=-5.176874e+01; StartPars[1]=9.056358e-01;
                StartPars[2]=8.557196e+01; StartPars[3]=8.044966e-01;
                break;
        case 1: f0=-1.1038; ef0=0.005*0.5; d0=-3.3492; ed0=0.05*0.5; FT=1./64.;
                StartPars[0]=-1.421800e+02; StartPars[1]=1.238405e+00;
                StartPars[2]=7.358287e+01; StartPars[3]=1.845156e+00;
                RedMass = (Mass_Kch*Mass_Dch)/(Mass_Kch+Mass_Dch);
                break;
        case 2: f0=-0.2633; ef0=0.002*1; d0=-2.2207; ed0=0.02*0.5; FT=0.25/64.;
                StartPars[0]=3.017167e+01; StartPars[1]=1.131973e+00;
                StartPars[2]=-1.253010e+00; StartPars[3]=1.777278e+00;
                RedMass = (Mass_Kch*Mass_Dch)/(Mass_Kch+Mass_Dch);
                break;
        case 3:
                //
                //f0=0.3782; ef0=0.02; d0=-0.9374; ed0=0.2; FT=1./8.;
                //StartPars[0]=9.214819e+01; StartPars[1]=1.124058e+00;
                //StartPars[2]=-3.264728e+00; StartPars[3]=1.809686e+00;
                f0=0.3782; ef0=0.001; d0=-0.9374; ed0=0.002; FT=1./256.;
                StartPars[0]=-9.255770e+01; StartPars[1]=9.124695e-01;
                StartPars[2]=1.613572e+01; StartPars[3]=1.406829e+00;
                RedMass = (Mass_Kch*Mass_Dch)/(Mass_Kch+Mass_Dch);
                //Current solution: V1=-8.8926e+01  mu1=1.0927e+00  V2=2.7415e+01  mu2=1.4687e+00 <--> f0=0.453  d0=-0.965
                //f0 = 0.378 fm
                //d0 = -0.94 fm
                // V1 = -9.255770e+01
                // mu1 = 9.124695e-01
                // V2 = 1.613572e+01
                // mu2 = 1.406829e+00

                break;
        case 4: f0=-0.1135; ef0=0.0005; d0=12.9067; ed0=0.2; FT=1./512.;
                //StartPars[0]=-3.704255e+01; StartPars[1]=9.367889e-01;
                //StartPars[2]=1.925511e+02; StartPars[3]=6.310400e-01;
                //V1=-3.5850e+01  mu1=1.0085e+00  V2=2.1450e+02  mu2=7.3332e-01 <--> f0=-0.113  d0=14.003
                StartPars[0]=-6.325317e+01; StartPars[1]=1.274817e+00;
                StartPars[2]=2.986266e+02; StartPars[3]=8.529017e-01;
                // Current solution: V1=-6.4271e+01  mu1=1.2460e+00  V2=2.8766e+02  mu2=8.4920e-01 <--> f0=-0.111  d0=11.222
                //Current solution: V1=-6.8380e+01  mu1=1.2647e+00  V2=2.9255e+02  mu2=8.6874e-01 <--> f0=-0.113  d0=12.336
                RedMass = (Mass_pic*Mass_Dch)/(Mass_pic+Mass_Dch);
                break;
        case 5: f0=0.3294; ef0=0.0004; d0=-3.5534; ed0=0.003; FT=1./512.;
                //StartPars[0]=-3.704255e+01; StartPars[1]=9.367889e-01;
                //StartPars[2]=1.925511e+02; StartPars[3]=6.310400e-01;
                //StartPars[0]=9.214819e+01; StartPars[1]=1.124058e+00;
                //StartPars[2]=-3.264728e+00; StartPars[3]=1.809686e+00;
                StartPars[0]=-1.465727e+02; StartPars[1]=1.175136e+00;
                StartPars[2]=2.718944e+01; StartPars[3]=1.832570e+00;
                RedMass = (Mass_pic*Mass_Dch)/(Mass_pic+Mass_Dch);
                break;
        default: printf("Weird flags for producing the potentials for EmmaDani\n"); return;
      }
    }
    //}
    //else{
      //NEGATIVE EFFECTIVE RANGE FOR THIS MODEL, WE PUT IT ON HOLD
    //}

    //ManufacturePotential(f0,ef0,d0,ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder,11,StartPars);
    //printf("Juli1\n");
    ManufacturePotential(f0,ef0,d0,ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder,11,StartPars,
      RedMass,FT);
      //printf("Juli2\n");
  }
  //toy models to demonstrate effects on cumulant
  //used mass - p Lambda
  else if(flag/10>=120&&flag/10<=139){
    NumR = 1;
    Radii[0] = 1.2;
    OutputFolder += "CumToy/";
    int VAR_FLAG = (flag/10)-120;
    double FT=1;
    double StartPars[5];
    StartPars[0]=100;
    StartPars[1]=0.5;
    StartPars[2]=-100;
    StartPars[3]=0.5;
    //if(flag/10<80){
    double RedMass = (Mass_p*Mass_L)/(Mass_p+Mass_L);
    switch (VAR_FLAG) {
      case 0: f0=0.25; ef0=f0*0.01; d0=0; ed0=1.0; FT=1./32.;//done
              //StartPars[0]=-3.003428e+01; StartPars[1]=6.119355e-01;
              //StartPars[2]=-3.860350e+02; StartPars[3]=2.728766e-01;
              StartPars[0]=-5.471780e+01; StartPars[1]=3.656692e-01;
              StartPars[2]=-5.085896e+02; StartPars[3]=2.699276e-01;
              break;
      case 1: f0=0.50; ef0=f0*0.01; d0=0; ed0=1.0; FT=1./32.;//done
              StartPars[0]=-2.491145e+01; StartPars[1]=6.789944e-01;
              StartPars[2]=-4.202059e+02; StartPars[3]=3.223011e-01;
              break;
      case 2: f0=0.75; ef0=f0*0.01; ed0=1.0; FT=1./32.;//done
      // Current solution: V1=7.5986e+00  mu1=7.6563e-01  V2=-7.6076e+01  mu2=8.0932e-01 <--> f0=0.750  d0=2.090
              StartPars[0]=-4.996459e+02; StartPars[1]=4.707076e-02;
              StartPars[2]=-1.119795e+03; StartPars[3]=2.561758e-01;
              break;
      case 3: f0=1.00; ef0=f0*0.01; d0=0; ed0=1.0; FT=1./32.;//done
              //StartPars[0]=-1.356962e+01; StartPars[1]=7.758964e-01;
              //StartPars[2]=-6.924030e+01; StartPars[3]=8.038734e-01;
              StartPars[0]=-4.819020e+02; StartPars[1]=6.510546e-02;
              StartPars[2]=-1.204676e+03; StartPars[3]=2.554564e-01;
              break;
      case 4: f0=1.5; ef0=f0*0.01; d0=0; ed0=1.0; FT=1./32.;//done
              StartPars[0]=-4.117917e+02; StartPars[1]=3.102178e-02;
              StartPars[2]=-1.100670e+03; StartPars[3]=2.753894e-01;
              break;
      case 5: f0=2.0; ef0=f0*0.01; d0=0; ed0=1.0; FT=1./32.;//done
              StartPars[0]=-4.372201e+02; StartPars[1]=1.757978e-01;
              StartPars[2]=-1.146283e+03; StartPars[3]=2.594236e-01;
              break;
      default: printf("Weird flags for producing the potentials for CumToy\n"); return;
    }

    ManufacturePotential(f0,ef0,d0,ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder,11,StartPars,
      RedMass,FT);
  }






  /*
  Suitable DoubleGaussSum potential found:
   f0 = -0.360 fm
   d0 = 1.98 fm
    V1 = -7.594562e+02
    mu1 = 6.466316e-01
    V2 = -1.200472e+02
    mu2 = 1.076839e+00

  Suitable DoubleGaussSum potential found:
   f0 = -0.353 fm
   d0 = 0.50 fm
    V1 = -1.693235e+03
    mu1 = 4.200525e-01
    V2 = -1.946644e+02
    mu2 = 7.536386e-01

  Suitable DoubleGaussSum potential found:
   f0 = -0.388 fm
   d0 = 3.51 fm
    V1 = -1.388675e+03
    mu1 = 9.657631e-01
    V2 = -1.072519e+03
    mu2 = 3.781579e-01

  Suitable DoubleGaussSum potential found:
   f0 = -0.399 fm
   d0 = 2.01 fm
    V1 = -1.471788e+03
    mu1 = 6.912759e-01
    V2 = -1.086707e+03
    mu2 = 8.669076e-01

  Suitable DoubleGaussSum potential found:
   f0 = -0.386 fm
   d0 = 0.49 fm
    V1 = -8.883082e+02
    mu1 = 4.924323e-01
    V2 = -3.805122e+02
    mu2 = 7.725780e-01
*/

  //Morita
  else if(flag==100){
    OutputFolder += "";
NumR=1;
    ef0 = 0.002;
    ed0 = 0.002;
    //NF48
    f0=1.511; d0=2.549; //ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,"YukawaDimiCore");
                        ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,"DoubleGaussSum",OutputFolder);
  }
  else if(flag<0){
    OutputFolder += "Temp/";
    NumR=1;
    ef0 = 0.05;
    ed0 = 0.05;
    if(flag==-1){
      f0 = 1.0; d0 = 1.0;
      ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,"DoubleGaussSum",OutputFolder);
      ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,"YukawaDimiCore",OutputFolder);
      ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,"Gaussian",OutputFolder);
      ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,"Yukawa",OutputFolder);
    }
    else if(flag==-2){//
      f0 = 1.0; d0 = 4.0;
      //ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,"DoubleGaussSum",OutputFolder);
      ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,"YukawaDimiCore",OutputFolder);
      //ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,"Gaussian",OutputFolder);
      ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,"Yukawa",OutputFolder);

    }
    else if(flag==-3){//
      f0 = -4.0; d0 = 1.0;
      //ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,"DoubleGaussSum",OutputFolder);
      //ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,"YukawaDimiCore",OutputFolder);
      //ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,"Gaussian",OutputFolder);
      ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,"Yukawa",OutputFolder,8);
    }
    else if(flag==-4){//
      f0 = -1.0; d0 = 1.0;
      //ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,"DoubleGaussSum",OutputFolder);
      ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,"YukawaDimiCore",OutputFolder);
      //ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,"Gaussian",OutputFolder);
      ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,"Yukawa",OutputFolder);
    }
    else if(flag==-5){//
      f0 = 4.0; d0 = 1.0;
      ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,"DoubleGaussSum",OutputFolder);
      ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,"YukawaDimiCore",OutputFolder);
      ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,"Gaussian",OutputFolder);
      ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,"Yukawa",OutputFolder);
    }
    else{
      f0 = 0.95; d0 = 8.9; ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,"DoubleGaussSum",OutputFolder);
    }






  }
  //ManufacturePotential(4.0,0.02,4.0,0.02,V1,mu1,V2,mu2,"DoubleGaussSum");

/*
  ManufacturePotential(0.25,0.00125,4.0,0.02,V1,mu1,V2,mu2,"DoubleGaussSum");
  ManufacturePotential(0.5,0.0025,4.0,0.02,V1,mu1,V2,mu2,"DoubleGaussSum");
  ManufacturePotential(1.0,0.005,4.0,0.02,V1,mu1,V2,mu2,"DoubleGaussSum");
  ManufacturePotential(2.0,0.01,4.0,0.02,V1,mu1,V2,mu2,"DoubleGaussSum");
  ManufacturePotential(4.0,0.02,4.0,0.02,V1,mu1,V2,mu2,"DoubleGaussSum");

  ManufacturePotential(0.25,0.00125,1.0,0.02,V1,mu1,V2,mu2,"DoubleGaussSum");
  ManufacturePotential(0.5,0.0025,1.0,0.02,V1,mu1,V2,mu2,"DoubleGaussSum");
  ManufacturePotential(1.0,0.005,1.0,0.02,V1,mu1,V2,mu2,"DoubleGaussSum");
  ManufacturePotential(2.0,0.01,1.0,0.02,V1,mu1,V2,mu2,"DoubleGaussSum");
  ManufacturePotential(4.0,0.02,1.0,0.02,V1,mu1,V2,mu2,"DoubleGaussSum");
*/

}

void pKplus_ScattPars(){
  CATS Kitty_I0;
  Kitty_I0.SetMomBins(41,0,123);


  CATS Kitty_I1;

  //Eval_ScattParameters


}
//double fDlmPot(const int& DlmPot, const int& DlmPotFlag,
//               const int& IsoSpin, const int& t2p1, const int& t2p2, const int& Spin,
//const int& AngMom, const int& TotMom, double* Radius, const double& CutOff, double* OtherPars){
//
/*
void pXi_Pot(){
  TGraph gr[4];
  unsigned nRad=0;
  double PotParsI0S0[9]={pXim_HALQCD1,-12,0,-1,1,0,0,0,0};
  double PotParsI0S1[9]={pXim_HALQCD1,-12,0,-1,1,1,0,1,0};
  double PotParsI1S0[9]={pXim_HALQCD1,-12,1,1,1,0,0,0,0};
  double PotParsI1S1[9]={pXim_HALQCD1,-12,1,1,1,1,0,1,0};
  CATSparameters* cPotParsI0S0 = new CATSparameters(CATSparameters::tPotential,9,true); cPotParsI0S0->SetParameters(PotParsI0S0);
  CATSparameters* cPotParsI0S1 = new CATSparameters(CATSparameters::tPotential,9,true); cPotParsI0S1->SetParameters(PotParsI0S1);
  CATSparameters* cPotParsI1S0 = new CATSparameters(CATSparameters::tPotential,9,true); cPotParsI1S0->SetParameters(PotParsI1S0);
  CATSparameters* cPotParsI1S1 = new CATSparameters(CATSparameters::tPotential,9,true); cPotParsI1S1->SetParameters(PotParsI1S1);
  for(double RAD = 0.01; RAD<3; RAD+=0.02){
    gr[0].SetPoint(nRad,RAD,fDlmPot.Eval(pXim_HALQCD1,-12,0,-1,1,0,0,0));
    gr[1].SetPoint(nRad,RAD,fDlmPot.Eval(pXim_HALQCD1,-12,0,-1,1,1,0,1));
    gr[2].SetPoint(nRad,RAD,fDlmPot.Eval(pXim_HALQCD1,-12,1,1,1,0,0,0));
    gr[3].SetPoint(nRad,RAD,fDlmPot.Eval(pXim_HALQCD1,-12,1,1,1,1,0,1));
    nRad++;
  }

}
*/


void Raffa_Errors(){

  TH1I* hSE1 = new TH1I("hSE1","hSE1",4,0,2);
  TH1I* hME1 = new TH1I("hME1","hME1",4,0,2);
  float Norm1 = 1.;

  TH1I* hSE2 = new TH1I("hSE2","hSE2",4,0,2);
  TH1I* hME2 = new TH1I("hME2","hME2",4,0,2);
  float Norm2 = -1.;

  hSE1->SetBinContent(1,20000);  hSE2->SetBinContent(1,20000+0.1*sqrt(20000.));
  hME1->SetBinContent(1,10000);  hME2->SetBinContent(1,10000);

  hSE1->SetBinContent(2,15000);  hSE2->SetBinContent(2,15000+1.5*sqrt(15000.));
  hME1->SetBinContent(2,10000);  hME2->SetBinContent(2,10000);

  hSE1->SetBinContent(3,10000);  hSE2->SetBinContent(3,10000-3.5*sqrt(10000.));
  hME1->SetBinContent(3,10000);  hME2->SetBinContent(3,10000);

  hSE1->SetBinContent(4,5000);   hSE2->SetBinContent(4,5000+0.5*sqrt(5000.));
  hME1->SetBinContent(4,10000);  hME2->SetBinContent(4,10000);

  HistoAddRatios HAD(2,"two_histos",1024,-0.2,0.2);
  HAD.SetNormalization(0,Norm1);
  HAD.SetNumerator(0,hSE1);
  HAD.SetDenominator(0,hME1);

  HAD.SetNormalization(1,Norm2);
  HAD.SetNumerator(1,hSE2);
  HAD.SetDenominator(1,hME2);
  HAD.SetIgnoreUncertainty(1);

  HAD.SetCompareToNullHypothesis(true,1000);
  TH2F* hResult = HAD.GetResult();

  for(unsigned uBin=0; uBin<hSE1->GetNbinsX(); uBin++){
    printf(" (%u) nsig = %.2f\n",uBin,HAD.GetNsig(uBin));
  }
  printf("nsig fish = %.2f\n",HAD.GetTotNsig(true));
  printf("nsig dimi = %.2f\n",HAD.GetTotNsig(false));

  TFile OutputFile(TString::Format("%s/OtherTasks/Raffa_Errors/OutputFile.root",GetFemtoOutputFolder()),"recreate");
  hResult->Write();

  delete hSE1;
  delete hSE2;
  delete hME1;
  delete hME2;
}

/*
//FeedDown = 0 NO
//FeedDown = 1 Yes, but a constant
//FeedDown = 2 Yes, sample form TH2F
void ppp_errors(int FeedDown, bool Projector){
  //the 3 body thing, I need to scale by 3 and than subrtract 2, related to how we build the cumulant
  const float Xmin = 0;
  const float Xmax = 0.4;
  const float Norm = 1./(0.5*(0.928459+0.929896));
  const float Norm_2 = 1./(0.5*(0.972599+0.972939));
  const float lam_gen = 0.618;
  const float lam_res = 0.196;
  const float lam_flt = 1.-lam_gen-lam_res;

  TString List1_p = "PairDist";
  TString ListREB_p = "PairRebinned";
  TString ListREW_p = "PairReweighted";
  TString histoSe_p = "sameEventDistributionPPPfSE_Shifted_FixShifted_Rebinned_30";
  TString histoMe_p = "mixedEventDistributionPPPfME_Shifted_FixShifted_Rebinned_30";
  TString histoSe_2_p = "sameEventDistributionppSamepMixedfSE_Shifted_FixShifted_Rebinned_30";
  TString histoMe_2_p = "mixedEventDistributionPPPfME_Shifted_FixShifted_Rebinned_30";

  TString List1_ap = "AntiPairDist";
  TString ListREB_ap = "PairRebinned";
  TString ListREW_ap = "PairReweighted";
  TString histoSe_ap = "sameEventDistributionAPAPAPfSE_Shifted_FixShifted_Rebinned_30";
  TString histoMe_ap = "mixedEventDistributionAPAPAPfME_Shifted_FixShifted_Rebinned_30";
  TString histoSe_2_ap = "sameEventDistributionapapSameapMixedfSE_Shifted_FixShifted_Rebinned_30";
  TString histoMe_2_ap = "mixedEventDistributionAPAPAPfME_Shifted_FixShifted_Rebinned_30";

  TString DataFileFolder = TString::Format("%s/3body/LauraPrelim/",GetCernBoxDimi());
  TList* list1_tmp;
  TList* list2_tmp;

  TFile* datafile_3 = new TFile(DataFileFolder+"CFOutput_PPP0Perfect.root");

  //PARTICLES
  list1_tmp = (TList*)datafile_3->Get(List1_p);
  //printf("list1_tmp=%p\n",list1_tmp);
  list2_tmp = (TList*)list1_tmp->FindObject(ListREB_p);
  //printf("list2_tmp=%p\n",list2_tmp);
  TH1F* hSe_3_p = (TH1F*)list2_tmp->FindObject(histoSe_p);
  //printf("hSe_3=%p\n",hSe_3_p);
  TH1F* hMe_3_p = (TH1F*)list2_tmp->FindObject(histoMe_p);
  //printf("hMe_3=%p\n",hMe_3_p);
  list2_tmp = (TList*)list1_tmp->FindObject(ListREW_p);
  //printf("list2_tmp=%p\n",list2_tmp);
  TH1F* hSeRew_3_p = (TH1F*)list2_tmp->FindObject(histoSe_p+"_Reweighted");
  //printf("hSeRew_3=%p\n",hSeRew_3_p);
  TH1F* hMeRew_3_p = (TH1F*)list2_tmp->FindObject(histoMe_p+"_Reweighted");
  //printf("hMeRew_3=%p\n",hMeRew_3_p);

  //ANTI-PARTICLES
  list1_tmp = (TList*)datafile_3->Get(List1_ap);
  //printf("list1_tmp=%p\n",list1_tmp);
  list2_tmp = (TList*)list1_tmp->FindObject(ListREB_ap);
  //printf("list2_tmp=%p\n",list2_tmp);
  TH1F* hSe_3_ap = (TH1F*)list2_tmp->FindObject(histoSe_ap);
  //printf("hSe_3=%p\n",hSe_3_ap);
  TH1F* hMe_3_ap = (TH1F*)list2_tmp->FindObject(histoMe_ap);
  //printf("hMe_3=%p\n",hMe_3_ap);
  list2_tmp = (TList*)list1_tmp->FindObject(ListREW_ap);
  //printf("list2_tmp=%p\n",list2_tmp);
  TH1F* hSeRew_3_ap = (TH1F*)list2_tmp->FindObject(histoSe_ap+"_Reweighted");
  //printf("hSeRew_3=%p\n",hSeRew_3_ap);
  TH1F* hMeRew_3_ap = (TH1F*)list2_tmp->FindObject(histoMe_ap+"_Reweighted");
  //printf("hMeRew_3=%p\n",hMeRew_3_ap);

  TH1F* hSe_3 = (TH1F*)hSe_3_p->Clone("hSe_3");
  hSe_3->Add(hSe_3_ap);

  TH1F* hMe_3 = (TH1F*)hMe_3_p->Clone("hMe_3");
  hMe_3->Add(hMe_3_ap);

  TH1F* hMeRew_3 = (TH1F*)hMeRew_3_p->Clone("hMeRew_3");
  hMeRew_3->Add(hMeRew_3_ap);

  //this is actually 1/weights, as this is the input we need for our code
  TH1F* hMe_3_w = (TH1F*)hMe_3->Clone("hMe_3_w");
  hMe_3_w->Divide(hMeRew_3);
  hMe_3_w->Sumw2();
  hMe_3_w->Scale(Norm);

  TH1F* hCk_3 = (TH1F*)hSe_3->Clone("hCk_3");
  hCk_3->Divide(hMeRew_3);
  hCk_3->Scale(Norm);

  //QA
  TFile* fOutput = new TFile(TString::Format("%s/OtherTasks/ppp_errors/fOutput_%i_%i.root",GetFemtoOutputFolder(),FeedDown,Projector),"recreate");
  hCk_3->Write();
  hMe_3_w->Write();


  //Next: the 2-body bit
  TFile* datafile_2 = new TFile(DataFileFolder+"CFOutput_ppSamepMixed0Perfect.root");

  //PARTICLES
  list1_tmp = (TList*)datafile_2->Get(List1_p);
  //printf("list1_tmp=%p\n",list1_tmp);
  list2_tmp = (TList*)list1_tmp->FindObject(ListREB_p);
  //printf("list2_tmp=%p\n",list2_tmp);
  //list2_tmp->ls();
  TH1F* hSe_2_p = (TH1F*)list2_tmp->FindObject(histoSe_2_p);
  //printf("hSe_2_p=%p\n",hSe_2_p);
  TH1F* hMe_2_p = (TH1F*)list2_tmp->FindObject(histoMe_2_p);
  //printf("hMe_2_p=%p\n",hMe_2_p);
  list2_tmp = (TList*)list1_tmp->FindObject(ListREW_p);
  //printf("list2_tmp=%p\n",list2_tmp);
  TH1F* hSeRew_2_p = (TH1F*)list2_tmp->FindObject(histoSe_2_p+"_Reweighted");
  //printf("hSeRew_2=%p\n",hSeRew_2_p);
  TH1F* hMeRew_2_p = (TH1F*)list2_tmp->FindObject(histoMe_2_p+"_Reweighted");
  //printf("hMeRew_2=%p\n",hMeRew_2_p);

  //ANTI-PARTICLES
  list1_tmp = (TList*)datafile_2->Get(List1_ap);
  //printf("list1_tmp=%p\n",list1_tmp);
  list2_tmp = (TList*)list1_tmp->FindObject(ListREB_ap);
  //printf("list2_tmp=%p\n",list2_tmp);
  //list2_tmp->ls();
  TH1F* hSe_2_ap = (TH1F*)list2_tmp->FindObject(histoSe_2_ap);
  //printf("hSe_2_ap=%p\n",hSe_2_ap);
  TH1F* hMe_2_ap = (TH1F*)list2_tmp->FindObject(histoSe_2_ap);
  //printf("hMe_2_ap=%p\n",hMe_2_ap);
  list2_tmp = (TList*)list1_tmp->FindObject(ListREW_ap);
  //printf("list2_tmp=%p\n",list2_tmp);
  TH1F* hSeRew_2_ap = (TH1F*)list2_tmp->FindObject(histoSe_2_ap+"_Reweighted");
  //printf("hSeRew_2=%p\n",hSeRew_2_ap);
  TH1F* hMeRew_2_ap = (TH1F*)list2_tmp->FindObject(histoMe_2_ap+"_Reweighted");
  //printf("hMeRew_2=%p\n",hMeRew_2_ap);

  //Next: the 2-body bit -> projector
  TFile* datafile_pr2 = new TFile(DataFileFolder+"CpppBINNED.root");
  //PARTICLES
  TGraphAsymmErrors* gr_Ck_2_pr = (TGraphAsymmErrors*)datafile_pr2->Get("htheo_ppp_stat");

  //the cumulant feed down
  TFile* datafile_f = new TFile(DataFileFolder+"DecaysCenteredGaussian2.root");
  TGraphMultiErrors* hFeedCum = (TGraphMultiErrors*)datafile_f->Get(FeedDown==1?"ppLFeedDownScaled":"ppLFeedDownScaled3Sigma");

  //the cumulant feed down (with probability density)
  TFile* datafile_fpdf = new TFile(DataFileFolder+"outputPPLSampled.root");
  list1_tmp = (TList*)list1_tmp->FindObject(ListREW_ap);
  TH2F* hCumPdf = (TH2F*)datafile_fpdf->Get("distributions2DSampledonPPPScaled");


  TH1F* hSe_2 = (TH1F*)hSe_2_p->Clone("hSe_2");
  hSe_2->Add(hSe_2_ap);

  TH1F* hMe_2 = (TH1F*)hMe_2_p->Clone("hMe_2");
  hMe_2->Add(hMe_2_ap);

  TH1F* hMeRew_2 = (TH1F*)hMeRew_2_p->Clone("hMeRew_2");
  hMeRew_2->Add(hMeRew_2_ap);

  //this is actually 1/weights, as this is the input we need for our code
  TH1F* hMe_2_w = (TH1F*)hMe_2->Clone("hMe_2_w");
  hMe_2_w->Divide(hMeRew_2);
  hMe_2_w->Sumw2();
  hMe_2_w->Scale(Norm_2);

  TH1F* hCk_2 = (TH1F*)hSe_2->Clone("hCk_2");
  hCk_2->Divide(hMeRew_2);
  //hCk_2->Sumw2();
  hCk_2->Scale(Norm_2);
  //for(unsigned uBin=0; uBin<hCk_2->GetNbinsX(); uBin++){
  //  hCk_2->SetBinContent(uBin+1,hCk_2->GetBinContent(uBin+1)-2);
  //}

  fOutput->cd();
  hCk_2->Write();
  hMe_2_w->Write();
  //hSe_2_p->Write();
  //hSe_2_ap->Write();
  //hSe_2->Write();
  //hMe_2_p->Write();
  //hMe_2_ap->Write();
  //hMe_2->Write();

  TH1F* hSe_3_exp;// = (TH1F*)hSe_2->Clone("hSe_3_exp");
  TH1F* hCtimes3minus2 = NULL;
  TH1F* hCtimes3minus2_m = NULL;
  //hSe_3_exp->Multiply(hMe_3);
  //hSe_3_exp->Divide(hMe_3_w);
  //hSe_3_exp->Divide(hMe_2);
  //hSe_3_exp->Multiply(hMe_2_w);
  //hSe_3_exp->Scale(3.);
  if(Projector){
    hCtimes3minus2 = (TH1F*)hMe_3->Clone("hCtimes3minus2");
    hSe_3_exp = (TH1F*)hMe_3->Clone("hSe_3_exp");
    hSe_3_exp->Divide(hMe_3_w);
    double Ctimes3minus2;
    double xDummy;
    for(unsigned uBin=0; uBin<gr_Ck_2_pr->GetN(); uBin++){
      gr_Ck_2_pr->GetPoint(uBin,xDummy,Ctimes3minus2);
      hCtimes3minus2->SetBinContent(uBin+1,Ctimes3minus2);
      hCtimes3minus2->SetBinError(uBin+1,gr_Ck_2_pr->GetErrorY(uBin));
      hSe_3_exp->SetBinContent(uBin+1,hSe_3_exp->GetBinContent(uBin+1)*Ctimes3minus2);
    }
    hCtimes3minus2_m = (TH1F*)hCtimes3minus2->Clone("hCtimes3minus2_m");
    hCtimes3minus2_m->Scale(-1.);
  }
  else{
    hSe_3_exp = (TH1F*)hSe_2->Clone("hSe_3_exp");
    hSe_3_exp->Multiply(hMe_2_w);
    hSe_3_exp->Scale(3.);
    hSe_3_exp->Add(hMe_2,-2);
    hSe_3_exp->Divide(hMe_2);
    hSe_3_exp->Multiply(hMe_3);
    hSe_3_exp->Divide(hMe_3_w);
  }
//double XTEST,YTEST;
//gr_Ck_2_pr->GetPoint(2,XTEST,YTEST);
//printf("d-p = %f vs %f\n",hSe_3_exp->GetBinContent(2),YTEST);
  hSe_3->Write();
  hSe_3_exp->Write();
  hMeRew_3->Write();

  TH1F* hMe_3exp_w = (TH1F*)hMe_3_w->Clone("hMe_3exp_w");
  hMe_3exp_w->Scale(-1.);

  TH1F* hMe_2_w_3x = (TH1F*)hMe_2_w->Clone("hMe_2_w_3x");
  hMe_2_w_3x->Scale(3.);
  TH1F* hMe_2_w_3xm = (TH1F*)hMe_2_w->Clone("hMe_2_w_3xm");
  hMe_2_w_3xm->Scale(-3.);
*/
/*
  printf("Last bin:\n");
  printf(" S3 = %.2f\n",hSe_3->GetBinContent(9));
  printf(" M3 = %.2f\n",hMe_3->GetBinContent(9));
  printf(" w3 = %.2f\n",hMe_3_w->GetBinContent(9));
  printf(" C3 = %.2f (%.2f)\n",hSe_3->GetBinContent(9)/hMe_3->GetBinContent(9)*hMe_3_w->GetBinContent(9),hCk_3->GetBinContent(9));
  printf(" S2 = %.2f\n",hSe_3_exp->GetBinContent(9));
  printf(" M2 = %.2f\n",hMe_3->GetBinContent(9));
  printf(" w2 = %.2f\n",hMe_3exp_w->GetBinContent(9));
  printf(" C2 = %.2f (%.2f)\n",hSe_3_exp->GetBinContent(9)/hMe_3->GetBinContent(9)*hMe_3exp_w->GetBinContent(9),hCk_2->GetBinContent(9));
  printf(" C2 (1) = %.2f\n",hMe_2_w->GetBinContent(9)*hSe_2->GetBinContent(9)/hMe_2->GetBinContent(9));
  printf(" Sexp = %.1f (%.1f)\n",hSe_3_exp->GetBinContent(9),3.*hMe_2_w->GetBinContent(9)*hSe_2->GetBinContent(9)/hMe_2->GetBinContent(9)*hMe_3->GetBinContent(9)/hMe_3_w->GetBinContent(9));
  printf("  c = %.2f\n",hSe_3->GetBinContent(9)/hMe_3->GetBinContent(9)*hMe_3_w->GetBinContent(9)+
    hSe_3_exp->GetBinContent(9)/hMe_3->GetBinContent(9)*hMe_3exp_w->GetBinContent(9)+2.);
*/
/*
//1.27-1.07*3+2 = 0.06
  TH1F* hFeed = (TH1F*)hSe_3->Clone("hFeed");
  if(FeedDown==0){
    for(unsigned uBin=0; uBin<hFeed->GetNbinsX(); uBin++){
      hFeed->SetBinContent(uBin+1,0);
    }
  }
  else if(FeedDown==1){
    for(unsigned uBin=0; uBin<hFeedCum->GetN(); uBin++){
      double xVal,yVal;
      //this is already scaled for lamda_res
      hFeedCum->GetPoint(uBin,xVal,yVal);
      hFeed->SetBinContent(uBin+1,-yVal);
    }
  }
  else{
    hFeed = NULL;
  }
  hMe_3_w->Scale(1./lam_gen);
  hMe_3exp_w->Scale(1./lam_gen);
  hCumPdf->Scale(1./lam_gen);
  //hFeed->Scale(1./lam_gen);
  //hMe_2_w->Scale(1./lam_gen);
  //hMe_2_w_minus->Scale(1./lam_gen);
  for(unsigned uBin=0; uBin<hSe_3->GetNbinsX(); uBin++){
    hSe_3->SetBinError(uBin+1,0);
    hMe_3->SetBinError(uBin+1,0);
    //double MuMe3 = hMe_3->GetBinContent(uBin+1);
    double MuSe2 = hSe_2->GetBinContent(uBin+1);
    double MuMe2 = hMe_2->GetBinContent(uBin+1);
    //double ErrMe3 = sqrt();
    //double ErrSe2 = sqrt(MuSe2);
    //double ErrMe2 = sqrt(MuMe2);
    double RelErrCk2;
    //the uncertainty related to the two-body expectiation, that is assumed to be Gaussian and with large number of entries in the ME
    if(Projector){RelErrCk2 = hCtimes3minus2->GetBinError(uBin+1)/hCtimes3minus2->GetBinContent(uBin+1);}
    else{
      RelErrCk2 = sqrt(1./(MuSe2*MuMe2)+1./MuSe2+1./MuMe2);
    }
    hSe_3_exp->SetBinError(uBin+1,hSe_3_exp->GetBinContent(uBin+1)*RelErrCk2);
    //printf(" err %u: %.1f%%\n",uBin,RelErrCk2*100.);
  }
  //TH1F* hNormMinus = (TH1F*)hSe_3->Clone("hSe_3");
  //for(unsigned uBin=0; uBin<hNormMinus->GetNbinsX(); uBin++){
  //  hNormMinus->SetBinContent(uBin+1,-1);
  //  hNormMinus->SetBinError(uBin+1,0);
  //}

  HistoAddRatios HAD(5,"hResult",8192,-24,8);
  HAD.SetNormalization(0,hMe_3_w);
  HAD.SetNumerator(0,hSe_3);
  HAD.SetDenominator(0,hMe_3);
  HAD.SetIgnoreUncertainty(0);

  //HAD.SetNormalization(1,hMe_2_w);
  //HAD.SetNumerator(1,hSe_2);
  //HAD.SetDenominator(1,hMe_2);
  HAD.SetNormalization(1,hMe_3exp_w);
  HAD.SetNumerator(1,hSe_3_exp);
  HAD.SetDenominator(1,hMe_3);

  //HAD.SetRatio(2,hFeed,true);
  if(hFeed){
    HAD.SetRatio(2,hFeed);
    HAD.SetIgnoreUncertainty(2);
  }
  else{
    HAD.SetRatio(2,hCumPdf);
    HAD.SetNormalization(2,-1.);
  }


/*
  if(Projector){
    HAD.SetRatio(3,hCtimes3minus2);

    HAD.SetRatio(4,hCtimes3minus2_m);
    HAD.SetIgnoreUncertainty(4);
  }
  else{
    HAD.SetNormalization(3,hMe_2_w_3x);
    HAD.SetNumerator(3,hSe_2);
    HAD.SetDenominator(3,hMe_2);

    HAD.SetNormalization(4,hMe_2_w_3xm);
    HAD.SetNumerator(4,hSe_2);
    HAD.SetDenominator(4,hMe_2);
    HAD.SetIgnoreUncertainty(4);
  }
*/
/*
  HAD.SetConstant(3,0,0);
  HAD.SetConstant(4,0,0);
  HAD.SetExpectationUncertainty(true);


  HAD.SetCompareToNullHypothesis(true,256);
  HAD.SetRange(Xmin,Xmax);
  HAD.SetNumIter(500000000);
  //HAD.SetNumIter(510000000);
  HAD.SetMinIter(100000);
  //HAD.SetNumIter(1000000);

  fOutput->cd();
  TH2F* hResult = HAD.GetResult();
  double Chi2_3b = 0;
  double Chi2_6b = 0;
  double Chi2_9b = 0;
  TGraphAsymmErrors* gaseResult_1s = new TGraphAsymmErrors();
  gaseResult_1s->SetName("gaseResult_1s");
  TGraphAsymmErrors* gaseResult_3s = new TGraphAsymmErrors();
  gaseResult_3s->SetName("gaseResult_3s");
  //double Chi2_29b = 0;//1 GeV
  for(unsigned uBin=0; uBin<9; uBin++){
    printf(" (%u) nsig (pval) = %.2f (%.5f)\n",uBin,HAD.GetNsig(uBin),HAD.GetPval(uBin));
    double med,low,up;
    double xval = hSe_3->GetBinCenter(uBin+1);
    double xwidth =  hSe_3->GetBinWidth(uBin+1);
    HAD.GetCentralInt(uBin,1,med,low,up);
    gaseResult_1s->SetPoint(uBin,xval,med);
    gaseResult_3s->SetPoint(uBin,xval,med);
    //printf("   mean value = %.2f\n",med);
    //printf("   1 sig = (%.2f, %.2f)\n",low,up);
    gaseResult_1s->SetPointError(uBin,xwidth/2.,xwidth/2.,med-low,up-med);
    HAD.GetCentralInt(uBin,3,med,low,up);
    //printf("   3 sig = (%.2f, %.2f)\n",low,up);
    gaseResult_3s->SetPointError(uBin,xwidth/2.,xwidth/2.,med-low,up-med);
    if(uBin<3) Chi2_3b += HAD.GetNsig(uBin)*HAD.GetNsig(uBin);
    if(uBin<6) Chi2_6b += HAD.GetNsig(uBin)*HAD.GetNsig(uBin);
    if(uBin<9) Chi2_9b += HAD.GetNsig(uBin)*HAD.GetNsig(uBin);
    //Chi2_29b += HAD.GetNsig(uBin)*HAD.GetNsig(uBin);
  }
  printf(" nsig 3 bins = %.2f\n", sqrt(2)*TMath::ErfcInverse(TMath::Prob(Chi2_3b,3)));
  printf(" nsig 6 bins = %.2f\n", sqrt(2)*TMath::ErfcInverse(TMath::Prob(Chi2_6b,6)));
  printf(" nsig 9 bins = %.2f\n", sqrt(2)*TMath::ErfcInverse(TMath::Prob(Chi2_9b,9)));
  //printf(" nsig 29 bins = %.2f\n", sqrt(2)*TMath::ErfcInverse(TMath::Prob(Chi2_29b,29)));

  fOutput->cd();
  hResult->Write();
  gaseResult_1s->Write();
  gaseResult_3s->Write();
}
*/

void ppL_errors(int FeedDown, bool Projector){
  //the 3 body thing, I need to scale by 3 and than subrtract 2, related to how we build the cumulant
  const float Xmin = 0;
  const float Xmax = 1.1;
  const float Norm = 1./(0.5*(0.951354+0.945446));
  const float Norm_2pp = 1./(0.5*(0.97576+0.976197));
  const float Norm_2pL = 1./(0.5*(0.983855+0.979637));

  TString List1_p = "PairDist";
  TString ListREB_p = "PairRebinned";
  TString ListREW_p = "PairReweighted";
  TString histoSe_p = "sameEventDistributionPPLfSE_Shifted_FixShifted_Rebinned_30";
  TString histoMe_p = "mixedEventDistributionPPLfME_Shifted_FixShifted_Rebinned_30";
  TString histoSe_2pp_p = "sameEventDistributionppSameLMixedfSE_Shifted_FixShifted_Rebinned_30";
  TString histoMe_2pp_p = "mixedEventDistributionPPLfME_Shifted_FixShifted_Rebinned_30";
  TString histoSe_2pL_p = "sameEventDistributionpLSamepMixedfSE_Shifted_FixShifted_Rebinned_30";
  TString histoMe_2pL_p = "mixedEventDistributionPPLfME_Shifted_FixShifted_Rebinned_30";

  TString List1_ap = "AntiPairDist";
  TString ListREB_ap = "PairRebinned";
  TString ListREW_ap = "PairReweighted";
  TString histoSe_ap = "sameEventDistributionAPAPALfSE_Shifted_FixShifted_Rebinned_30";
  TString histoMe_ap = "mixedEventDistributionAPAPALfME_Shifted_FixShifted_Rebinned_30";
  TString histoSe_2pp_ap = "sameEventDistributionapapSameaLMixedfSE_Shifted_FixShifted_Rebinned_30";
  TString histoMe_2pp_ap = "mixedEventDistributionAPAPALfME_Shifted_FixShifted_Rebinned_30";
  TString histoSe_2pL_ap = "sameEventDistributionapaLSameapMixedfSE_Shifted_FixShifted_Rebinned_30";
  TString histoMe_2pL_ap = "mixedEventDistributionAPAPALfME_Shifted_FixShifted_Rebinned_30";

  TString DataFileFolder = TString::Format("%s/3body/LauraPrelim/",GetCernBoxDimi());
  TList* list1_tmp;
  TList* list2_tmp;

  TFile* datafile_3 = new TFile(DataFileFolder+"CFOutput_PPL0Perfect.root");

  //PARTICLES
  list1_tmp = (TList*)datafile_3->Get(List1_p);
  printf("list1_tmp=%p\n",list1_tmp);
  list2_tmp = (TList*)list1_tmp->FindObject(ListREB_p);
  printf("list2_tmp=%p\n",list2_tmp);
  TH1F* hSe_3_p = (TH1F*)list2_tmp->FindObject(histoSe_p);
  printf("hSe_3=%p\n",hSe_3_p);
  TH1F* hMe_3_p = (TH1F*)list2_tmp->FindObject(histoMe_p);
  printf("hMe_3=%p\n",hMe_3_p);
  list2_tmp = (TList*)list1_tmp->FindObject(ListREW_p);
  printf("list2_tmp=%p\n",list2_tmp);
  TH1F* hSeRew_3_p = (TH1F*)list2_tmp->FindObject(histoSe_p+"_Reweighted");
  printf("hSeRew_3=%p\n",hSeRew_3_p);
  TH1F* hMeRew_3_p = (TH1F*)list2_tmp->FindObject(histoMe_p+"_Reweighted");
  printf("hMeRew_3=%p\n",hMeRew_3_p);

  //ANTI-PARTICLES
  list1_tmp = (TList*)datafile_3->Get(List1_ap);
  printf("list1_tmp=%p\n",list1_tmp);
  list2_tmp = (TList*)list1_tmp->FindObject(ListREB_ap);
  printf("list2_tmp=%p\n",list2_tmp);
  TH1F* hSe_3_ap = (TH1F*)list2_tmp->FindObject(histoSe_ap);
  printf("hSe_3=%p\n",hSe_3_ap);
  TH1F* hMe_3_ap = (TH1F*)list2_tmp->FindObject(histoMe_ap);
  printf("hMe_3=%p\n",hMe_3_ap);
  list2_tmp = (TList*)list1_tmp->FindObject(ListREW_ap);
  printf("list2_tmp=%p\n",list2_tmp);
  TH1F* hSeRew_3_ap = (TH1F*)list2_tmp->FindObject(histoSe_ap+"_Reweighted");
  printf("hSeRew_3=%p\n",hSeRew_3_ap);
  TH1F* hMeRew_3_ap = (TH1F*)list2_tmp->FindObject(histoMe_ap+"_Reweighted");
  printf("hMeRew_3=%p\n",hMeRew_3_ap);

  TH1F* hSe_3 = (TH1F*)hSe_3_p->Clone("hSe_3");
  hSe_3->Add(hSe_3_ap);

  TH1F* hMe_3 = (TH1F*)hMe_3_p->Clone("hMe_3");
  hMe_3->Add(hMe_3_ap);

  TH1F* hMeRew_3 = (TH1F*)hMeRew_3_p->Clone("hMeRew_3");
  hMeRew_3->Add(hMeRew_3_ap);

  //this is actually 1/weights, as this is the input we need for our code
  TH1F* hMe_3_w = (TH1F*)hMe_3->Clone("hMe_3_w");
  hMe_3_w->Divide(hMeRew_3);
  hMe_3_w->Sumw2();
  hMe_3_w->Scale(Norm);

  TH1F* hCk_3 = (TH1F*)hSe_3->Clone("hCk_3");
  hCk_3->Divide(hMeRew_3);
  hCk_3->Scale(Norm);


  //QA
  TFile* fOutput = new TFile(TString::Format("%s/OtherTasks/ppL_errors/fOutput_%i_%i.root",GetFemtoOutputFolder(),FeedDown,Projector),"recreate");
  hCk_3->Write();
  hMe_3_w->Write();


  //Next: the 2-body bit (pp Lmixed)
  TFile* datafile_2pp = new TFile(DataFileFolder+"CFOutput_ppSameLMixed0Perfect.root");

  //PARTICLES
  list1_tmp = (TList*)datafile_2pp->Get(List1_p);
  printf("list1_tmp=%p\n",list1_tmp);
  list2_tmp = (TList*)list1_tmp->FindObject(ListREB_p);
  printf("list2_tmp=%p\n",list2_tmp);
  //list2_tmp->ls();
  TH1F* hSe_2pp_p = (TH1F*)list2_tmp->FindObject(histoSe_2pp_p);
  printf("hSe_2pp_p=%p\n",hSe_2pp_p);
  TH1F* hMe_2pp_p = (TH1F*)list2_tmp->FindObject(histoMe_2pp_p);
  printf("hMe_2pp_p=%p\n",hMe_2pp_p);
  list2_tmp = (TList*)list1_tmp->FindObject(ListREW_p);
  printf("list2_tmp=%p\n",list2_tmp);
  TH1F* hSeRew_2pp_p = (TH1F*)list2_tmp->FindObject(histoSe_2pp_p+"_Reweighted");
  printf("hSeRew_2pp_p=%p\n",hSeRew_2pp_p);
  TH1F* hMeRew_2pp_p = (TH1F*)list2_tmp->FindObject(histoMe_2pp_p+"_Reweighted");
  printf("hMeRew_2pp_p=%p\n",hMeRew_2pp_p);

  //ANTI-PARTICLES
  list1_tmp = (TList*)datafile_2pp->Get(List1_ap);
  printf("list1_tmp=%p\n",list1_tmp);
  list2_tmp = (TList*)list1_tmp->FindObject(ListREB_ap);
  printf("list2_tmp=%p\n",list2_tmp);
  //list2_tmp->ls();
  TH1F* hSe_2pp_ap = (TH1F*)list2_tmp->FindObject(histoSe_2pp_ap);
  printf("hSe_2pp_ap=%p\n",hSe_2pp_ap);
  TH1F* hMe_2pp_ap = (TH1F*)list2_tmp->FindObject(histoSe_2pp_ap);
  printf("hMe_2pp_ap=%p\n",hMe_2pp_ap);
  list2_tmp = (TList*)list1_tmp->FindObject(ListREW_ap);
  printf("list2_tmp=%p\n",list2_tmp);
  TH1F* hSeRew_2pp_ap = (TH1F*)list2_tmp->FindObject(histoSe_2pp_ap+"_Reweighted");
  printf("hSeRew_2pp_p=%p\n",hSeRew_2pp_ap);
  TH1F* hMeRew_2pp_ap = (TH1F*)list2_tmp->FindObject(histoMe_2pp_ap+"_Reweighted");
  printf("hMeRew_2pp_p=%p\n",hMeRew_2pp_ap);


  TH1F* hSe_2pp = (TH1F*)hSe_2pp_p->Clone("hSe_2pp");
  hSe_2pp->Add(hSe_2pp_ap);

  TH1F* hMe_2pp = (TH1F*)hMe_2pp_p->Clone("hMe_2pp");
  hMe_2pp->Add(hMe_2pp_ap);

  TH1F* hMeRew_2pp = (TH1F*)hMeRew_2pp_p->Clone("hMeRew_2pp");
  hMeRew_2pp->Add(hMeRew_2pp_ap);

  //this is actually 1/weights, as this is the input we need for our code
  TH1F* hMe_2pp_w = (TH1F*)hMe_2pp->Clone("hMe_2pp_w");
  hMe_2pp_w->Divide(hMeRew_2pp);
  hMe_2pp_w->Sumw2();
  hMe_2pp_w->Scale(Norm_2pp);

  TH1F* hCk_2pp = (TH1F*)hSe_2pp->Clone("hCk_2pp");
  hCk_2pp->Divide(hMeRew_2pp);
  hCk_2pp->Scale(Norm_2pp);

  //Next: the 2-body bit (pL pmixed)
  TFile* datafile_2pL = new TFile(DataFileFolder+"CFOutput_pLSamepMixed0Perfect.root");

  //PARTICLES
  list1_tmp = (TList*)datafile_2pL->Get(List1_p);
  printf("list1_tmp=%p\n",list1_tmp);
  list2_tmp = (TList*)list1_tmp->FindObject(ListREB_p);
  printf("list2_tmp=%p\n",list2_tmp);
  //list2_tmp->ls();
  TH1F* hSe_2pL_p = (TH1F*)list2_tmp->FindObject(histoSe_2pL_p);
  printf("hSe_2pL_p=%p\n",hSe_2pL_p);
  TH1F* hMe_2pL_p = (TH1F*)list2_tmp->FindObject(histoMe_2pL_p);
  printf("hMe_2pL_p=%p\n",hMe_2pL_p);
  list2_tmp = (TList*)list1_tmp->FindObject(ListREW_p);
  printf("list2_tmp=%p\n",list2_tmp);
  TH1F* hSeRew_2pL_p = (TH1F*)list2_tmp->FindObject(histoSe_2pL_p+"_Reweighted");
  printf("hSeRew_2pL_p=%p\n",hSeRew_2pL_p);
  TH1F* hMeRew_2pL_p = (TH1F*)list2_tmp->FindObject(histoMe_2pL_p+"_Reweighted");
  printf("hMeRew_2pL_p=%p\n",hMeRew_2pL_p);

  //ANTI-PARTICLES
  list1_tmp = (TList*)datafile_2pL->Get(List1_ap);
  printf("list1_tmp=%p\n",list1_tmp);
  list2_tmp = (TList*)list1_tmp->FindObject(ListREB_ap);
  printf("list2_tmp=%p\n",list2_tmp);
  //list2_tmp->ls();
  TH1F* hSe_2pL_ap = (TH1F*)list2_tmp->FindObject(histoSe_2pL_ap);
  printf("hSe_2pL_ap=%p\n",hSe_2pL_ap);
  TH1F* hMe_2pL_ap = (TH1F*)list2_tmp->FindObject(histoSe_2pL_ap);
  printf("hMe_2pL_ap=%p\n",hMe_2pL_ap);
  list2_tmp = (TList*)list1_tmp->FindObject(ListREW_ap);
  printf("list2_tmp=%p\n",list2_tmp);
  TH1F* hSeRew_2pL_ap = (TH1F*)list2_tmp->FindObject(histoSe_2pL_ap+"_Reweighted");
  printf("hSeRew_2pL_p=%p\n",hSeRew_2pL_ap);
  TH1F* hMeRew_2pL_ap = (TH1F*)list2_tmp->FindObject(histoMe_2pL_ap+"_Reweighted");
  printf("hMeRew_2pL_p=%p\n",hMeRew_2pL_ap);

  //Next: the 2-body bit -> projector
  TFile* datafile_pr2 = new TFile(DataFileFolder+"CppL_BINNED.root");
  //PARTICLES
  TGraphAsymmErrors* gr_Ck_2_pr = (TGraphAsymmErrors*)datafile_pr2->Get("htheo_ppL_stat");

  TH1F* hSe_2pL = (TH1F*)hSe_2pL_p->Clone("hSe_2pL");
  hSe_2pL->Add(hSe_2pL_ap);

  TH1F* hMe_2pL = (TH1F*)hMe_2pL_p->Clone("hMe_2pL");
  hMe_2pL->Add(hMe_2pL_ap);

  TH1F* hMeRew_2pL = (TH1F*)hMeRew_2pL_p->Clone("hMeRew_2pL");
  hMeRew_2pL->Add(hMeRew_2pL_ap);

  //this is actually 1/weights, as this is the input we need for our code
  TH1F* hMe_2pL_w = (TH1F*)hMe_2pL->Clone("hMe_2pL_w");
  hMe_2pL_w->Divide(hMeRew_2pL);
  hMe_2pL_w->Sumw2();
  hMe_2pL_w->Scale(Norm_2pL);

  TH1F* hCk_2pL = (TH1F*)hSe_2pL->Clone("hCk_2pL");
  hCk_2pL->Divide(hMeRew_2pL);
  hCk_2pL->Scale(Norm_2pL);


  //decay matrix, used only for testing
  TFile* datafile_decay = new TFile(DataFileFolder+"DLM_DecayMatrix_ppL_ppp.root");
  TH2F* hDecay_ppL_ppp = (TH2F*)datafile_decay->Get("DLM_DecayMatrix");

  TH1F* hMe_3_decay = (TH1F*)hMe_3->Clone("hMe_3_decay");
  TH1F* hMe_3_decayDimi = (TH1F*)hMe_3->Clone("hMe_3_decayDimi");
  for(unsigned uBinY=0; uBinY<hMe_3->GetNbinsX(); uBinY++){
    hMe_3_decay->SetBinContent(uBinY+1,0);
    double kMother,kDaughter;
    kDaughter = hMe_3->GetBinCenter(uBinY+1);
    if(kDaughter>Xmax) break;
    for(unsigned uBinX=0; uBinX<hMe_3->GetNbinsX(); uBinX++){
      int FirstBin = TMath::Nint(hMe_3->GetXaxis()->GetBinLowEdge(uBinX+1)*1000.);
      int LastBin = TMath::Nint(hMe_3->GetXaxis()->GetBinUpEdge(uBinX+1)*1000.)-1;
      int Rebin = LastBin-FirstBin+1;
      kMother = hMe_3->GetBinCenter(uBinX+1);
      if(kMother>Xmax) break;
      TH1D* hTransform = hDecay_ppL_ppp->ProjectionY(TString::Format("hTransform"),FirstBin,LastBin);
      //TH1F* hTransf = new TH1F("hTransf","hTransf",65,hTransform->GetBinLowEdge(1),hTransform->GetXaxis()->GetBinUpEdge(65));
      //for(unsigned uBin=0; uBin<65; uBin++){
      //  hTransf->SetBinContent(uBin+1,hTransform->GetBinContent(uBin*6+1));
      //  for(unsigned uBinReb=1; uBinReb<6; uBinReb++){
      //    hTransf->SetBinContent(uBin+1,hTransf->GetBinContent(uBin+1)+hTransform->GetBinContent(uBin*6+1+uBinReb));
      //  }
      //}
      //printf(" Rebin=%i\n",Rebin);
      //printf(" b: NBX=%i; l/u=%.3f/%.3f\n",hTransform->GetNbinsX(),hTransform->GetXaxis()->GetBinLowEdge(1),hTransform->GetXaxis()->GetBinUpEdge(1));
      hTransform->Rebin(Rebin);
      hTransform->Scale(1./hTransform->Integral());
      //printf(" a: NBX=%i; l/u=%.3f/%.3f\n",hTransform->GetNbinsX(),hTransform->GetXaxis()->GetBinLowEdge(1),hTransform->GetXaxis()->GetBinUpEdge(1));
      //usleep(200e3);
      unsigned Yield = hMe_3_decay->GetBinContent(uBinX+1);
      unsigned NewYield = hMe_3->GetBinContent(uBinX+1)*hTransform->GetBinContent(uBinX+1);
      hMe_3_decay->SetBinContent(uBinX+1,Yield+NewYield);
      delete hTransform;
      //delete hTransf;
    }
  }

  for(unsigned uBin=0; uBin<hMe_3->GetNbinsX(); uBin++){
    hMe_3_decayDimi->SetBinContent(uBin+1,0);
  }
  for(unsigned uBin=0; uBin<hMe_3->GetNbinsX(); uBin++){
    double kMother,kDaughter;
    unsigned Yield = hMe_3->GetBinContent(uBin+1);
    int FirstBin = TMath::Nint(hMe_3->GetXaxis()->GetBinLowEdge(uBin+1)*1000.);
    int LastBin = TMath::Nint(hMe_3->GetXaxis()->GetBinUpEdge(uBin+1)*1000.)-1;
    int Rebin = LastBin-FirstBin+1;
    kMother = hMe_3->GetBinCenter(uBin+1);
    if(kMother>Xmax) break;
    TH1D* hTransform = hDecay_ppL_ppp->ProjectionY(TString::Format("hTransform"),FirstBin,LastBin);
    //TH1F* hTransf = new TH1F("hTransf","hTransf",65,hTransform->GetBinLowEdge(1),hTransform->GetXaxis()->GetBinUpEdge(65));
    //for(unsigned uBin=0; uBin<65; uBin++){
    //  hTransf->SetBinContent(uBin+1,hTransform->GetBinContent(uBin*6+1));
    //  for(unsigned uBinReb=1; uBinReb<6; uBinReb++){
    //    hTransf->SetBinContent(uBin+1,hTransf->GetBinContent(uBin+1)+hTransform->GetBinContent(uBin*6+1+uBinReb));
    //  }
    //}
    hTransform->Rebin(Rebin);
    hTransform->Scale(1./hTransform->Integral());
    for(unsigned uYield=0; uYield<Yield; uYield++){
      kDaughter = hTransform->GetRandom();
      hMe_3_decayDimi->Fill(kDaughter);
    }
    delete hTransform;
    //delete hTransf;
  }

  TH1F* hMe_3_decayfactor = (TH1F*)hMe_3->Clone("hMe_3_decayfactor");
  hMe_3_decayfactor->Divide(hMe_3_decay);
  TH1F* hMe_3_decayfactorDimi = (TH1F*)hMe_3->Clone("hMe_3_decayfactorDimi");
  hMe_3_decayfactorDimi->Divide(hMe_3_decayDimi);

  fOutput->cd();
  //hCk_2->Write();
  //hMe_2_w->Write();
  //hSe_2_p->Write();
  //hSe_2_ap->Write();
  //hSe_2->Write();
  //hMe_2_p->Write();
  //hMe_2_ap->Write();
  //hMe_2->Write();
  hDecay_ppL_ppp->Write();
  hMe_3->Write();
  hMe_3_decay->Write();
  hMe_3_decayDimi->Write();
  hMe_3_decayfactor->Write();
  hMe_3_decayfactorDimi->Write();

  TH1F* hSe_3_exp; //= (TH1F*)hSe_2pp->Clone("hSe_3_exp");
  TH1F* hCtimes3minus2 = NULL;
  TH1F* hCtimes3minus2_m = NULL;
  //hSe_3_exp->Multiply(hMe_3);
  //hSe_3_exp->Divide(hMe_3_w);
  //hSe_3_exp->Divide(hMe_2);
  //hSe_3_exp->Multiply(hMe_2_w);
  //hSe_3_exp->Scale(3.);
  if(Projector){
    hCtimes3minus2 = (TH1F*)hMe_3->Clone("hCtimes3minus2");
    hSe_3_exp = (TH1F*)hMe_3->Clone("hSe_3_exp");
    hSe_3_exp->Divide(hMe_3_w);
    double Ctimes3minus2;
    double xDummy;
    for(unsigned uBin=0; uBin<gr_Ck_2_pr->GetN(); uBin++){
      gr_Ck_2_pr->GetPoint(uBin,xDummy,Ctimes3minus2);
      /*
      if(uBin==0){
        Ctimes3minus2 += 2.4;
      }
      if(uBin==1){
        Ctimes3minus2 -= 1.9;
      }
      if(uBin==2){
        Ctimes3minus2 -= 1.38;
      }
      if(uBin==3){
        Ctimes3minus2 -= 1.17;
      }
      if(uBin==4){//
        Ctimes3minus2 -= 1.28;
      }
      if(uBin==5){//
        Ctimes3minus2 -= 0.83;
      }
      if(uBin==6){
        Ctimes3minus2 -= 0.57;
      }
      if(uBin==7){
        Ctimes3minus2 -= 0.428;
      }
      if(uBin==8){
        Ctimes3minus2 -= 0.297;
      }
      */
      hCtimes3minus2->SetBinContent(uBin+1,Ctimes3minus2);
      hCtimes3minus2->SetBinError(uBin+1,gr_Ck_2_pr->GetErrorY(uBin));
      hSe_3_exp->SetBinContent(uBin+1,hSe_3_exp->GetBinContent(uBin+1)*Ctimes3minus2);
    }
    hCtimes3minus2_m = (TH1F*)hCtimes3minus2->Clone("hCtimes3minus2_m");
    hCtimes3minus2_m->Scale(-1.);
  }
  else{
    hSe_3_exp = (TH1F*)hSe_2pp->Clone("hSe_3_exp");
    for(unsigned uBin=0; uBin<hSe_3_exp->GetNbinsX(); uBin++){
      if(hSe_3_exp->GetBinCenter(uBin+1)>Xmax) break;
      double BinValue_pp = hSe_2pp->GetBinContent(uBin+1);
      BinValue_pp *= hMe_2pp_w->GetBinContent(uBin+1);
      BinValue_pp /= hMe_2pp->GetBinContent(uBin+1);
      double BinValue_pL = hSe_2pL->GetBinContent(uBin+1);
      BinValue_pL *= hMe_2pL_w->GetBinContent(uBin+1);
      BinValue_pL /= hMe_2pL->GetBinContent(uBin+1);
      double BinValue = 2.*BinValue_pL+BinValue_pp-2.;
      BinValue *= hMe_3->GetBinContent(uBin+1);
      BinValue /= hMe_3_w->GetBinContent(uBin+1);
      hSe_3_exp->SetBinContent(uBin+1,BinValue);
    }
  }
  hSe_3->Write();
  hSe_3_exp->Write();
  hMeRew_3->Write();


  TH1F* hMe_3exp_w = (TH1F*)hMe_3_w->Clone("hMe_3exp_w");
  hMe_3exp_w->Scale(-1);

  TH1F* hMe_2pp_w_1x = (TH1F*)hMe_2pp_w->Clone("hMe_2pp_w_1x");
  //hMe_2_wpp_3x->Scale(1.);
  TH1F* hMe_2pp_w_1xm = (TH1F*)hMe_2pp_w->Clone("hMe_2pp_w_1xm");
  hMe_2pp_w_1xm->Scale(-1.);

  TH1F* hMe_2pL_w_2x = (TH1F*)hMe_2pL_w->Clone("hMe_2pL_w_2x");
  hMe_2pL_w_2x->Scale(2.);
  TH1F* hMe_2pL_w_2xm = (TH1F*)hMe_2pL_w->Clone("hMe_2pL_w_2xm");
  hMe_2pL_w_2xm->Scale(-2.);


  HistoAddRatios HAD(7,"hResult",8192,-8,40);
  HAD.SetNormalization(0,hMe_3_w);
  HAD.SetNumerator(0,hSe_3);
  HAD.SetDenominator(0,hMe_3);
  HAD.SetIgnoreUncertainty(0);

  //HAD.SetNormalization(1,hMe_2_w);
  //HAD.SetNumerator(1,hSe_2);
  //HAD.SetDenominator(1,hMe_2);
  HAD.SetNormalization(1,hMe_3exp_w);
  HAD.SetNumerator(1,hSe_3_exp);
  HAD.SetDenominator(1,hMe_3);

  HAD.SetConstant(2,0,0);

  if(Projector){
    HAD.SetRatio(3,hCtimes3minus2);

    HAD.SetRatio(4,hCtimes3minus2_m);
    HAD.SetIgnoreUncertainty(4);

    HAD.SetConstant(5,0,0);
    HAD.SetConstant(6,0,0);
  }
  else{
    HAD.SetNormalization(3,hMe_2pp_w_1x);
    HAD.SetNumerator(3,hSe_2pp);
    HAD.SetDenominator(3,hMe_2pp);

    HAD.SetNormalization(4,hMe_2pp_w_1xm);
    HAD.SetNumerator(4,hSe_2pp);
    HAD.SetDenominator(4,hMe_2pp);
    HAD.SetIgnoreUncertainty(4);

    HAD.SetNormalization(5,hMe_2pL_w_2x);
    HAD.SetNumerator(5,hSe_2pL);
    HAD.SetDenominator(5,hMe_2pL);

    HAD.SetNormalization(6,hMe_2pL_w_2xm);
    HAD.SetNumerator(6,hSe_2pL);
    HAD.SetDenominator(6,hMe_2pL);
    HAD.SetIgnoreUncertainty(6);
  }

  HAD.SetCompareToNullHypothesis(true,256);
  HAD.SetRange(Xmin,Xmax);
  HAD.SetNumIter(100000000);
  HAD.SetMinIter(10000);


  fOutput->cd();
  TH2F* hResult = HAD.GetResult();
  double Chi2_3b = 0;
  double Chi2_6b = 0;
  for(unsigned uBin=0; uBin<9; uBin++){
    printf(" (%u) nsig (pval) = %.2f (%.5f)\n",uBin,HAD.GetNsig(uBin),HAD.GetPval(uBin));
    if(uBin<3) Chi2_3b += HAD.GetNsig(uBin)*HAD.GetNsig(uBin);
    if(uBin<6) Chi2_6b += HAD.GetNsig(uBin)*HAD.GetNsig(uBin);
  }
  printf(" nsig 3 bins = %.2f\n", sqrt(2)*TMath::ErfcInverse(TMath::Prob(Chi2_3b,3)));
  printf(" nsig 6 bins = %.2f\n", sqrt(2)*TMath::ErfcInverse(TMath::Prob(Chi2_6b,6)));

  fOutput->cd();
  hResult->Write();

}


//Emission: 0, all particles are emitted at dist rSize (point source), otherwise Gaussian source
//PhaseSpace: 0 is uniform, 1 is boltzmann
void Silly2body(const int Emission, const int PhaseSpace){
  const unsigned NumPairs = 100*1000;
  const double rCutOff = 10;
  const double tCutOff = 1e6;
  const double rSize = 1.0;
  //the radius fluctuates (say QM) to avoid zero positions of the potential
  const double rFluct = 0.02;

  //the potential is A*exp(-x^2/R^2)
  const double V_A = -20;
  const double V_R = 2.0;
  const double V_CA = 500;
  const double V_CR = 0.5;
  const double RedMass = 500;
  const double MomMean = 1000;
  const double tEps = 1e-2;
  const double tMin = 1e-2;
  const double tMax = 1e4;

  double Rad;
  double Mom;
  double Force;
  double dt;
  double dr;

  DLM_Random DimiRanGen(11);

  unsigned NBINS = 1024;
  TFile* fOutput = new TFile(TString::Format("%s/OtherTasks/Silly2body/fOutput_%i_%i.root",GetFemtoOutputFolder(),Emission,PhaseSpace),"recreate");
  TH1F* hSource = new TH1F("hSource","hSource",NBINS,0,8);
  TH1F* hRef = new TH1F("hRef","hRef",NBINS,0,4096);
  TH1F* hSig = new TH1F("hSig","hSig",NBINS,0,4096);
  TH1F* hDist = new TH1F("hDist","hDist",NBINS,0,8*rCutOff);
  TH1F* hTime = new TH1F("hTime","hTime",NBINS,0,4*tCutOff);
  TH1F* hPos = new TH1F("hPos","hPos",NBINS,0,4*rCutOff);

  for(unsigned uPair=0; uPair<NumPairs; uPair++){
    double Prog = double(uPair)/double(NumPairs);
    //printf("uPair=%u\n",uPair);
    printf("\r\033[K Progress %.0f%% ",100.*Prog);
    cout << flush;
    if(Emission==0) Rad=rSize;
    else {Rad=DimiRanGen.GaussDiffR(3,0,rSize);}
    double Rad0 = Rad;

    if(PhaseSpace==0) {Mom=DimiRanGen.Uniform(0,2.*MomMean);}
    else {Mom=DimiRanGen.GaussR(3,0,MomMean/sqrt(2));}
    Mom = 100;
    double Mom0 = Mom;

    double DistanceTraveled = 0;
    double TimeTraveled = 0;
    //sign conventions. If Mom>0 we move away from each other, else towards one another
    //if the force is negative => attraction => the absolute value of the Mom should decrease
    //the sign of the radius is arbitrary, in the end we do start off with positive values, and
    //in the end we should only care about |Rad|
    while(fabs(Rad)<rCutOff&&DistanceTraveled<2.*rCutOff&&TimeTraveled<tCutOff){
      Force = V_A*exp(-pow(fabs(Rad)/V_R,2.))+V_CA*exp(-pow(fabs(Rad)/V_CR,2.));
      //if we are stuck around a zero potential value within the range
      //if(fabs(Force)<tEps&&fabs(Rad)<V_R){
      //  Force = tEps*Force/fabs(Force);
      //}
      //dt = fabs(Mom/Force)*(sqrt(1.+fabs(2.*tEps*Rad*RedMass*Force))-1.);
      dt = fabs(Mom/Force*tEps);
      if(dt>tMax) dt=tMax;
      if(dt<tMin) dt=tMin;
      dr = Mom/RedMass*dt+0.5*Force/RedMass*dt*dt;
      /*if(uPair==252){
        printf("Mom = %.2f\n",Mom);
        printf("  F = %.2f\n",Force);
        printf(" dt = %.5f\n",dt);
        printf(" dp = %.2f\n",Force*dt);
        printf("Rad = %.2f\n",Rad);
        printf(" dr = %.2f\n",dt);
        printf("\n");
        usleep(25e3);
      }*/

      Rad += dr;
      Rad += DimiRanGen.Gauss(0,rFluct);
      DistanceTraveled += fabs(dr);
      TimeTraveled += dt;
      if(Mom>0) Mom += Force*dt;
      else Mom -= Force*dt;
    }

    //if(fabs(Mom)<1){
      //printf(" DistanceTraveled = %.2f\n",DistanceTraveled);
      //printf(" TimeTraveled = %.2f\n",TimeTraveled);
      //printf(" Rad = %.2f\n",Rad);
      //printf(" Mom = %.2f\n",Mom);
      //printf("\n");
    //}

    //this is a pair that got stuck with 0 momentum and close position
    //experimentally this is not really detectable, so we remove these pairs
    //if(Mom>1||fabs(Rad)>rCutOff){
      hDist->Fill(DistanceTraveled);
      hSig->Fill(fabs(Mom));
      hTime->Fill(TimeTraveled);
      hPos->Fill(fabs(Rad));
      hSource->Fill(Rad0);
      hRef->Fill(Mom0);
    //}

  }
  printf("\r\033[K Progress %.0f%% ",100.);

  TH1F* hCk = (TH1F*)hSig->Clone("hCk");
  hCk->Divide(hRef);

  hSource->Write();
  hSig->Write();
  hRef->Write();
  hCk->Write();
  hPos->Write();
  hDist->Write();
  hTime->Write();

  delete hSource;
  delete hSig;
  delete hRef;
  delete hCk;
  delete hPos;
  delete hDist;
  delete hTime;
  delete fOutput;
}


void ppL_smear_test(){

}

void KilledByFeedDown(){

  const double kMin = 0;
  const double kMax = 768;
  //8 MeV bins
  const double nBins = 64;
  const double BinWidth = kMax/double(nBins);
  const double lam_S0 = 0.18;
  const double lam_Xi = 0.221;
  const double lam_pL = 0.427;
  const double lam_flt = 1.-lam_S0-lam_Xi-lam_pL;
  double lam_S0_genuine = 1.-0.18;//this is the flat contribution to pSigma0 itself

  TFile* fOutput = new TFile(TString::Format("%s/OtherTasks/KilledByFeedDown/fOutput.root",GetFemtoOutputFolder()),"recreate");
  TH1F* hCk_pS0 = new TH1F("hCk_pS0","hCk_pS0",nBins,kMin,kMax);
  TH1F* hCk_pS0_pL_old = new TH1F("hCk_pS0_pL_old","hCk_pS0_pL_old",nBins,kMin,kMax);
  TH1F* hCk_pS0_pL_new = new TH1F("hCk_pS0_pL_new","hCk_pS0_pL_new",nBins,kMin,kMax);
  TH1F* hCk_pXim = new TH1F("hCk_pXim","hCk_pXim",nBins,kMin,kMax);
  TH1F* hCk_pXi0 = new TH1F("hCk_pXi0","hCk_pXi0",nBins,kMin,kMax);
  TH1F* hCk_pXi = new TH1F("hCk_pXi","hCk_pXi",nBins,kMin,kMax);
  TH1F* hCk_pXi_pL_old = new TH1F("hCk_pXi_pL_old","hCk_pXi_pL_old",nBins,kMin,kMax);
  TH1F* hCk_pXi_pL_new = new TH1F("hCk_pXi_pL_new","hCk_pXi_pL_new",nBins,kMin,kMax);

  TH1F* hCk_pS0_pL_dlm = new TH1F("hCk_pS0_pL_dlm","hCk_pS0_pL_dlm",nBins,kMin,kMax);
  TH1F* hCk_pXi_pL_dlm = new TH1F("hCk_pXi_pL_dlm","hCk_pXi_pL_dlm",nBins,kMin,kMax);

  TFile* fDM = new TFile(TString::Format("%s/CatsFiles/DecaySmear/Decay_matrices_2020.root",GetCernBoxDimi()));
  TH2F* hDM_pS0_pL = (TH2F*)fDM->Get("hRes_pL_pSigma0");
  TH2F* hDM_pXi_pL = (TH2F*)fDM->Get("hRes_pL_pXim");

  TFile* fME = new TFile(TString::Format("%s/CatsFiles/ExpData/ALICE_pp_13TeV_HM/DimiJun20/Norm240_340/DataSignal/CkREW_pL_0.root",GetCernBoxDimi()));
  TDirectoryFile *dirME=(TDirectoryFile*)(fME->FindObjectAny(Form("Binning_%i",TMath::Nint(BinWidth))));
  TH1F *hMe_pL;
  dirME->GetObject(Form("hME_SUM"),hMe_pL);
  TH1F* hMe_pL_S0_smeared = (TH1F*)hMe_pL->Clone("hMe_pL_S0_smeared");
  TH1F* hMe_pL_Xi_smeared = (TH1F*)hMe_pL->Clone("hMe_pL_Xi_smeared");

  DLM_CleverMcLevyResoTM CleverMcLevyReso;
  CleverMcLevyReso.InitStability(21,1,2);
  CleverMcLevyReso.InitScale(91,100,1000);
  CleverMcLevyReso.InitRad(513,0,8192);
  CleverMcLevyReso.InitType(2);
  CleverMcLevyReso.InitNumMcIter(1000000);
  //TF1* fMe = new TF1("fMe",&CleverMcLevyReso,&DLM_CleverMcLevyResoTM::RootEvalNorm,0,1000,3,"DLM_CleverMcLevyResoTM","RootEvalNorm");
  TF1* fitMe = new TF1("fitMe","[2]*(x*x*exp(-pow(x/sqrt(2)/[0],[1])))/pow([0],3.)",0,1000);
  fitMe->SetParameter(0,500);
  fitMe->SetParameter(1,1.5);
  fitMe->SetParLimits(1,0,5);
  fitMe->SetParameter(2,1);
  hMe_pL->Fit(fitMe,"S, N, R, M");
fitMe->SetParameter(0,fitMe->GetParameter(0)*0.9);
//[2]*(x*x*exp(-pow(x/sqrt(2)/[0],[1])))/pow([0],3.)


  TH1F* hGHETTO_PS = new TH1F("hGHETTO_PS","hGHETTO_PS",250,0,3000);
  for(unsigned uBin=0; uBin<250; uBin++){
    double MOM = hGHETTO_PS->GetBinCenter(uBin+1);
    //double VAL = 5.40165e-08-1.06309e-08*pow(MOM,1.)+5.71623e-09*pow(MOM,2.)-
    //-1.15454e-11*pow(MOM,3.)+1.17292e-14*pow(MOM,4.)-7.94266e-18*pow(MOM,5.)+
    //3.70109e-21*pow(MOM,6.)-1.10644e-24*pow(MOM,7.)+1.86899e-28*pow(MOM,8.)-1.34242e-32*pow(MOM,9.);
    double VAL = fitMe->Eval(MOM);
    hGHETTO_PS->SetBinContent(uBin+1,VAL);
  }


  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());

  double lam_pXim[5];
  AnalysisObject.SetUpLambdaPars_pXim("pp13TeV_HM_DimiJun20",0,0,lam_pXim);
  double lam_Xim_genuine = (1.-lam_pXim[1]-lam_pXim[2]-lam_pXim[3]-lam_pXim[4])/(1.-lam_pXim[4]);
  double lam_Xim_Flat = (lam_pXim[1]+lam_pXim[2]+lam_pXim[3])/(1.-lam_pXim[4]);
  printf("lam_Xim_genuine = %.2f\n",lam_Xim_genuine*100.);
  printf("lam_Xim_Flat = %.2f\n",lam_Xim_Flat*100.);

  CATS AB_pL;
  AB_pL.SetMomBins(nBins,kMin,kMax);
  AnalysisObject.SetUpCats_pL(AB_pL,"Chiral_Coupled_SPD","McGauss_ResoTM",11600,202);//NLO_Coupled_S

  const double CuspWeight = 0.4;//0.54
  AB_pL.SetChannelWeight(7,1./4.*CuspWeight);//1S0 SN(s) -> LN(s)
  AB_pL.SetChannelWeight(8,3./4.*CuspWeight);//3S1 SN(s) -> LN(s)
  AB_pL.SetChannelWeight(10,3./4.*CuspWeight);//3S1 SN(d) -> LN(s)
  AB_pL.SetChannelWeight(13,3./20.*CuspWeight);//3D1 SN(d) -> LN(d)
  AB_pL.SetChannelWeight(15,3./20.*CuspWeight);//3D1 SN(s) -> LN(d)

  AB_pL.SetAnaSource(0,1.09);
  AB_pL.SetAnaSource(1,2);
  AB_pL.SetNotifications(CATS::nError);
  AB_pL.KillTheCat();

  DLM_Ck* Ck_pL = new DLM_Ck(AB_pL.GetNumSourcePars(),0,AB_pL,nBins,kMin,kMax);
  Ck_pL->SetSourcePar(0,AB_pL.GetAnaSourcePar(0));
  Ck_pL->SetSourcePar(1,AB_pL.GetAnaSourcePar(1));
  Ck_pL->SetCutOff(330,kMax);
  Ck_pL->Update();

  CATS AB_pXim;
  //same binning as pL, as we only use pXim as feed-down
  AB_pXim.SetMomBins(nBins,kMin,kMax);
  AnalysisObject.SetUpCats_pXim(AB_pXim,"pXim_HALQCDPaper2020","Gauss");
//!
  AB_pXim.SetAnaSource(0,1.0);
  AB_pXim.SetNotifications(CATS::nError);
  AB_pXim.KillTheCat();

  DLM_Ck* Ck_pXim = new DLM_Ck(AB_pXim.GetNumSourcePars(),0,AB_pXim,nBins,kMin,kMax);
  Ck_pXim->SetSourcePar(0,AB_pXim.GetAnaSourcePar(0));
  Ck_pXim->SetSourcePar(1,AB_pXim.GetAnaSourcePar(1));
  Ck_pXim->SetCutOff(330,kMax);
  Ck_pXim->Update();

  CATS AB_pXi0;
  //same binning as pL, as we only use pXim as feed-down
  AB_pXi0.SetMomBins(nBins,kMin,kMax);
  AnalysisObject.SetUpCats_pXi0(AB_pXi0,"pXim_HALQCDPaper2020","Gauss");
//!
  AB_pXi0.SetAnaSource(0,1.0);
  AB_pXi0.SetNotifications(CATS::nError);
  AB_pXi0.KillTheCat();

  DLM_Ck* Ck_pXi0 = new DLM_Ck(AB_pXi0.GetNumSourcePars(),0,AB_pXi0,nBins,kMin,kMax);
  Ck_pXi0->SetSourcePar(0,AB_pXi0.GetAnaSourcePar(0));
  Ck_pXi0->SetSourcePar(1,AB_pXi0.GetAnaSourcePar(1));
  Ck_pXi0->SetCutOff(330,kMax);
  Ck_pXi0->Update();

  CATS AB_pS0_Chiral;
  //the minus one is to to avoid going above 350 MeV, since we do not have the WF there
  AB_pS0_Chiral.SetMomBins(nBins,kMin,kMax);
  AnalysisObject.SetUpCats_pS0(AB_pS0_Chiral,"Chiral","Gauss");
//!
  AB_pS0_Chiral.SetAnaSource(0,1.25);
  AB_pS0_Chiral.SetNotifications(CATS::nError);
  AB_pS0_Chiral.KillTheCat();
  DLM_Ck* Ck_pS0_Chiral = new DLM_Ck(AB_pS0_Chiral.GetNumSourcePars(),0,AB_pS0_Chiral,nBins,kMin,kMax);
  Ck_pS0_Chiral->SetSourcePar(0,1.25);
  Ck_pS0_Chiral->SetCutOff(330,kMax);
  Ck_pS0_Chiral->Update();

  TH2F* hResidual_pL_pSigma0 = AnalysisObject.GetResidualMatrix("pLambda","pSigma0");
  TH2F* hResidual_pL_pXim = AnalysisObject.GetResidualMatrix("pLambda","pXim");
  TH2F* hResidual_pL_pXi0 = AnalysisObject.GetResidualMatrix("pLambda","pXi0");

  DLM_CkDecomposition* CkDec_pL = new DLM_CkDecomposition("pLambda",5,*Ck_pL,NULL);
  DLM_CkDecomposition* CkDec_pS0 = new DLM_CkDecomposition("pSigma0",1,*Ck_pS0_Chiral,NULL);;
  DLM_CkDecomposition* CkDec_pXim = new DLM_CkDecomposition("pXim",2,*Ck_pXim,NULL);
  DLM_CkDecomposition* CkDec_pXi0 = new DLM_CkDecomposition("pXi0",2,*Ck_pXi0,NULL);

  CkDec_pL->AddContribution(0,lam_S0,DLM_CkDecomposition::cFeedDown,CkDec_pS0,hResidual_pL_pSigma0);
  CkDec_pL->AddContribution(1,lam_Xi/2.,DLM_CkDecomposition::cFeedDown,CkDec_pXim,hResidual_pL_pXim);
  CkDec_pL->AddContribution(2,lam_Xi/2.,DLM_CkDecomposition::cFeedDown,CkDec_pXi0,hResidual_pL_pXi0);
  CkDec_pL->AddContribution(3,lam_flt,DLM_CkDecomposition::cFeedDown);
  CkDec_pL->AddContribution(4,0,DLM_CkDecomposition::cFake);

  CkDec_pL->AddPhaseSpace(0,hGHETTO_PS);
  CkDec_pL->AddPhaseSpace(1,hGHETTO_PS);
  CkDec_pL->AddPhaseSpace(2,hGHETTO_PS);
  //CkDec_pL->AddPhaseSpace(0,hMe_pL);
  //CkDec_pL->AddPhaseSpace(1,hMe_pL);
  //CkDec_pL->AddPhaseSpace(2,hMe_pL);
//TF1* fMe = new TF1("fMe",GaussSourceScaledTF1,0,1000,2);


//TF1* fData = new TF1("fData",dimi_pL_May2020_Fitter,FitRegion_pL[0],FitRegion_pL[3],NumFitPars);
//TH1F* testPS0 = Convert_DlmHisto_TH1F(CkDec_pL->PS_Child[0],"testPS0");


  CkDec_pS0->AddContribution(0,1.-lam_S0_genuine,DLM_CkDecomposition::cFeedDown);

  CkDec_pXim->AddContribution(0,lam_pXim[1]+lam_pXim[2]+lam_pXim[3],DLM_CkDecomposition::cFeedDown);
  CkDec_pXim->AddContribution(1,lam_pXim[4],DLM_CkDecomposition::cFake);

  CkDec_pXi0->AddContribution(0,lam_pXim[1]+lam_pXim[2]+lam_pXim[3],DLM_CkDecomposition::cFeedDown);
  CkDec_pXi0->AddContribution(1,lam_pXim[4],DLM_CkDecomposition::cFake);

  CkDec_pL->Update(true,true);

  for(unsigned uBin=0; uBin<nBins; uBin++){
    double MOM = hCk_pXi_pL_dlm->GetBinCenter(uBin+1);
    hCk_pS0_pL_dlm->SetBinContent(uBin+1,CkDec_pL->EvalSignalChild(0,MOM)+1.);
    hCk_pXi_pL_dlm->SetBinContent(uBin+1,CkDec_pL->EvalSignalChild(1,MOM)+CkDec_pL->EvalSignalChild(2,MOM)+1.);
    //hCk_pXi_pL_dlm->SetBinContent(uBin+1,CkDec_pL->EvalCk(MOM));
  }



  for(unsigned uBin=0; uBin<nBins; uBin++){
    hCk_pXim->SetBinContent(uBin+1,AB_pXim.GetCorrFun(uBin));
    hCk_pXi0->SetBinContent(uBin+1,AB_pXi0.GetCorrFun(uBin));
    hCk_pXi->SetBinContent(uBin+1,0.5*(AB_pXim.GetCorrFun(uBin)+AB_pXi0.GetCorrFun(uBin)));
    hCk_pS0->SetBinContent(uBin+1,Ck_pS0_Chiral->GetBinContent(uBin));
  }
//printf("hi\n");
  //smear the shit
  //X = pS0, Y = pL
  for(unsigned uBinY=0; uBinY<nBins; uBinY++){
    double BWDM = hDM_pS0_pL->GetXaxis()->GetBinWidth(1);
    hDM_pS0_pL->Rebin2D(TMath::Nint(BinWidth/BWDM),TMath::Nint(BinWidth/BWDM));
    TH1D* hTrS0 = hDM_pS0_pL->ProjectionX(TString::Format("hTrS0"),uBinY+1,uBinY+1);
    //hDM_pXi_pL
    hTrS0->Scale(1./hTrS0->Integral());

    BWDM = hDM_pXi_pL->GetXaxis()->GetBinWidth(1);
    hDM_pXi_pL->Rebin2D(TMath::Nint(BinWidth/BWDM),TMath::Nint(BinWidth/BWDM));
    TH1D* hTrXi = hDM_pXi_pL->ProjectionX(TString::Format("hTrXi"),uBinY+1,uBinY+1);
    hTrXi->Scale(1./hTrXi->Integral());

    double BinValueS0 = 0;
    double BinValueXi = 0;
    double BinValueMeS0 = 0;
    double BinValueMeXi = 0;
    for(unsigned uBinX=0; uBinX<nBins; uBinX++){
      BinValueS0 += hTrS0->GetBinContent(uBinX+1)*Ck_pS0_Chiral->GetBinContent(uBinX);
      BinValueXi += hTrXi->GetBinContent(uBinX+1)*hCk_pXi->GetBinContent(uBinX+1);
      BinValueMeS0 += hTrS0->GetBinContent(uBinX+1)*hMe_pL->GetBinContent(uBinX+1);
      BinValueMeXi += hTrXi->GetBinContent(uBinX+1)*hMe_pL->GetBinContent(uBinX+1);
    }
    hCk_pS0_pL_old->SetBinContent(uBinY+1,lam_S0*lam_S0_genuine*BinValueS0+1.-lam_S0*lam_S0_genuine);
    hCk_pXi_pL_old->SetBinContent(uBinY+1,lam_Xi*lam_Xim_genuine*BinValueXi+1.-lam_Xi*lam_Xim_genuine);
    hMe_pL_S0_smeared->SetBinContent(uBinY+1,BinValueMeS0);
    hMe_pL_Xi_smeared->SetBinContent(uBinY+1,BinValueMeXi);
    delete hTrS0;
    delete hTrXi;
  }
  //I think the norm is fucked up due to the tail, renorm in tighter region
  for(unsigned uBin=(nBins*3)/4; uBin<nBins; uBin++){
    hMe_pL_S0_smeared->SetBinContent(uBin+1,hMe_pL->GetBinContent(uBin+1));
    hMe_pL_Xi_smeared->SetBinContent(uBin+1,hMe_pL->GetBinContent(uBin+1));
  }
  //hMe_pL_smeared->Scale(hMe_pL->Integral(1,(nBins*2)/3)/hMe_pL_smeared->Integral(1,(nBins*2)/3));
  //printf(" IME=%.4e; ISME=%.4e\n",hMe_pL->Integral(),hMe_pL_smeared->Integral());
  //hMe_pL_smeared->Scale(hMe_pL->Integral()/hMe_pL_smeared->Integral());
//printf("hi2\n");
  //new smearing
  for(unsigned uBinY=0; uBinY<nBins; uBinY++){
    double BWDM = hDM_pS0_pL->GetXaxis()->GetBinWidth(1);
    hDM_pS0_pL->Rebin2D(TMath::Nint(BinWidth/BWDM),TMath::Nint(BinWidth/BWDM));
    TH1D* hTrS0 = hDM_pS0_pL->ProjectionX(TString::Format("hTrS0"),uBinY+1,uBinY+1);
    hTrS0->Scale(1./hTrS0->Integral());

    BWDM = hDM_pXi_pL->GetXaxis()->GetBinWidth(1);
    hDM_pXi_pL->Rebin2D(TMath::Nint(BinWidth/BWDM),TMath::Nint(BinWidth/BWDM));
    TH1D* hTrXi = hDM_pXi_pL->ProjectionX(TString::Format("hTrXi"),uBinY+1,uBinY+1);
    hTrXi->Scale(1./hTrXi->Integral());

    double BinValueS0 = 0;
    double BinValueXi = 0;
    for(unsigned uBinX=0; uBinX<nBins; uBinX++){
      BinValueS0 +=  hTrS0->GetBinContent(uBinX+1)*
                      hMe_pL->GetBinContent(uBinX+1)/
                      hMe_pL_S0_smeared->GetBinContent(uBinY+1)*
                      Ck_pS0_Chiral->GetBinContent(uBinX);
      BinValueXi +=  hTrXi->GetBinContent(uBinX+1)*
                      hMe_pL->GetBinContent(uBinX+1)/
                      hMe_pL_Xi_smeared->GetBinContent(uBinY+1)*
                      hCk_pXi->GetBinContent(uBinX+1);
    }
    hCk_pS0_pL_new->SetBinContent(uBinY+1,lam_S0*BinValueS0+1.-lam_S0);
    hCk_pXi_pL_new->SetBinContent(uBinY+1,lam_Xi*lam_Xim_genuine*BinValueXi+1.-lam_Xi*lam_Xim_genuine);
    delete hTrS0;
    delete hTrXi;
  }
//printf("hi3\n");
  fOutput->cd();
  hGHETTO_PS->Write();
  hMe_pL->Write();
  fitMe->Write();
  hMe_pL_S0_smeared->Write();
  hMe_pL_Xi_smeared->Write();
  hCk_pS0->Write();
  hCk_pS0_pL_old->Write();
  hCk_pS0_pL_new->Write();
  hCk_pS0_pL_dlm->Write();
  hCk_pXim->Write();
  hCk_pXi0->Write();
  hCk_pXi->Write();
  hCk_pXi_pL_old->Write();
  hCk_pXi_pL_new->Write();
  hCk_pXi_pL_dlm->Write();
//testPS0->Write();

/*
  delete hCk_pS0;
  delete hCk_pS0_pL_old;
  delete hCk_pS0_pL_new;
  delete hCk_pXim;
  delete hCk_pXi0;
  delete hCk_pXi;
  delete hCk_pXi_pL_old;
  delete hCk_pXi_pL_new;
  delete Ck_pS0_Chiral;
  delete hMe_pL_S0_smeared;
  delete hMe_pL_Xi_smeared;
  delete fOutput;
  delete fDM;
  delete fME;
  delete hGHETTO_PS;
*/
}


//DLM_CkDecomposition* MM_PP;
//CATS* MM_CatPP;
double VerySimpleFitter_pp(double* x, double* par){
    double& MOM = *x;
    MM_PP->GetCk()->SetSourcePar(0,par[0]);
    MM_PP->Update(false);
    double CkVal;
    CkVal = MM_PP->EvalCk(MOM);
    return CkVal*(par[1]+MOM*par[2]);
}

void SmearEffect_pp(){
  const double ResidualSourceSize = 1.35;
  const double SourceSize = 1.25;

  double* MomBins_pp = NULL;
  double* FitRegion_pp = NULL;
  unsigned NumMomBins_pp;
  TString DataSample = "pp13TeV_HM_Dec19";
  //TString DataSample = "pp13TeV_HM_RotPhiDec19";
  DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());
  AnalysisObject.SetUpBinning_pp(DataSample,NumMomBins_pp,MomBins_pp,FitRegion_pp);
NumMomBins_pp=30;
for(unsigned uBin=0; uBin<=NumMomBins_pp; uBin++){
  MomBins_pp[uBin] = double(uBin)*12.;
}

  TH2F* hResolution_pp = AnalysisObject.GetResolutionMatrix(DataSample,"pp");
  TH2F* hResidual_pp_pL = AnalysisObject.GetResidualMatrix("pp","pLambda");

  TH1F* hData_pp = AnalysisObject.GetAliceExpCorrFun(DataSample,"pp","_0",0,false,-1);

  double lam_pp[5];
  double lam_pL[5];

  AnalysisObject.SetUpLambdaPars_pp(DataSample,0,lam_pp);
  AnalysisObject.SetUpLambdaPars_pL(DataSample,0,0,lam_pL);

  //TString SourceDescription = "Gauss";
  TString OutputFolder = TString::Format("%s/OtherTasks/SmearEffect_pp/",GetFemtoOutputFolder());
  TString OutFileName = "fOutput.root";
  TFile* OutputFile = new TFile(OutputFolder+OutFileName,"recreate");
  hData_pp->Write();

  //Gauss
  TString SourceDescription = "Gauss";
  //TString SourceDescription = "McGauss_ResoTM";
  //TString SourceDescription = "McLevy_ResoTM";

  CATS AB_pp;
  DLM_Ck* Ck_pp;
  AB_pp.SetMomBins(NumMomBins_pp,MomBins_pp);
  AB_pp.SetNotifications(CATS::nWarning);
  //AnalysisObject.SetUpCats_pp(AB_pp,"AV18","Gauss",0,0);//McLevyNolan_Reso
  AnalysisObject.SetUpCats_pp(AB_pp,"AV18",SourceDescription,0,202);//McLevyNolan_Reso

  AB_pp.SetAnaSource(0,SourceSize);
  if(SourceDescription.Contains("Mc")){
      AB_pp.SetAnaSource(0,1.20);
      AB_pp.SetAnaSource(1,2.0);
  }
//AB_pp.RemoveShortRangePotential();
//AB_pp.SetQ1Q2(0);
//AB_pp.SetQuantumStatistics(false);
  //AB_pp.SetAnaSource(1,2.0);
  //AB_pp.SetNotifications(CATS::nWarning);
  //AB_pp.SetEpsilonConv(5e-8);
  //AB_pp.SetEpsilonProp(5e-8);
  AB_pp.KillTheCat();
  Ck_pp = new DLM_Ck(AB_pp.GetNumSourcePars(),0,AB_pp);
  Ck_pp->Update();

  CATS AB_pL;
  DLM_Ck* Ck_pL;
  AB_pL.SetMomBins(NumMomBins_pp,MomBins_pp);
  //AnalysisObject.SetUpCats_pL(AB_pL,"NLO_Coupled_S","Gauss");
  AnalysisObject.SetUpCats_pL(AB_pL,"Usmani","Gauss");
  AB_pL.SetAnaSource(0,ResidualSourceSize);
  AB_pL.SetNotifications(CATS::nWarning);
//AB_pL.RemoveShortRangePotential();
  AB_pL.KillTheCat();
  Ck_pL = new DLM_Ck(AB_pL.GetNumSourcePars(),0,AB_pL);
  Ck_pL->Update();

  DLM_CkDecomposition CkDec_pp("pp",3,*Ck_pp,hResolution_pp);
  //DLM_CkDecomposition CkDec_pp("pp",3,*Ck_pp,NULL);
  DLM_CkDecomposition CkDec_pL("pLambda",2,*Ck_pL,hResolution_pp);

  CkDec_pp.AddContribution(0,lam_pp[1]*0,DLM_CkDecomposition::cFeedDown,&CkDec_pL,hResidual_pp_pL);
  CkDec_pp.AddContribution(1,lam_pp[2]*0,DLM_CkDecomposition::cFeedDown);
  CkDec_pp.AddContribution(2,lam_pp[3]*0,DLM_CkDecomposition::cFake);

  CkDec_pL.AddContribution(2,lam_pL[1]+lam_pL[2]+lam_pL[3],DLM_CkDecomposition::cFeedDown);
  CkDec_pL.AddContribution(3,lam_pL[4],DLM_CkDecomposition::cFake);//0.03

  CkDec_pp.Update();
  CkDec_pL.Update();

  MM_PP = &CkDec_pp;
  TF1* fit_pp_data_old = new TF1("fit_pp_data_old",VerySimpleFitter_pp,0,240,3);
  fit_pp_data_old->SetParameter(0,SourceSize);
  fit_pp_data_old->SetParameter(1,1);
  fit_pp_data_old->FixParameter(2,0);
  hData_pp->Fit(fit_pp_data_old,"S, N, R, M");
  CkDec_pp.GetCk()->SetSourcePar(0,SourceSize);

  TH1F* h_pL_pp_old = new TH1F("h_pL_pp_old","h_pL_pp_old",NumMomBins_pp,MomBins_pp);
  TH1F* h_pL_pp_new = new TH1F("h_pL_pp_new","h_pL_pp_new",NumMomBins_pp,MomBins_pp);
  TH1F* h_pp_old = new TH1F("h_pp_old","h_pp_old",NumMomBins_pp,MomBins_pp);
  TH1F* h_pp_new = new TH1F("h_pp_new","h_pp_new",NumMomBins_pp,MomBins_pp);

  TH1F* h_pL_old = new TH1F("h_pL_old","h_pL_old",NumMomBins_pp,MomBins_pp);
  TH1F* h_pL_new = new TH1F("h_pL_new","h_pL_new",NumMomBins_pp,MomBins_pp);

  printf("lam main = %.2f\n",CkDec_pp.GetLambdaMain()*100.);
  for(unsigned uBin=0; uBin<NumMomBins_pp; uBin++){
    double MOM = h_pL_pp_old->GetBinCenter(uBin+1);
    h_pp_old->SetBinContent(uBin+1,CkDec_pp.EvalCk(MOM));
    h_pp_old->SetBinError(uBin+1,0.005);
    h_pL_old->SetBinContent(uBin+1,CkDec_pL.EvalCk(MOM));
    h_pL_old->SetBinError(uBin+1,0.005);
    if(MOM<50) h_pp_old->SetBinError(uBin+1,0.05);
    h_pL_pp_old->SetBinContent(uBin+1,CkDec_pp.EvalSignalSmearedChild(0,MOM)+1.);

    if(MOM>60) continue;
    printf(" oECK(%.0f) = %.3f = %.3f+1 = %.3f = %.3f+%.3f+1 = %.3f (%.3f)\n",MOM,
            CkDec_pp.EvalCk(MOM),
            CkDec_pp.EvalSignalSmeared(MOM),
            CkDec_pp.EvalSignalSmeared(MOM)+1,
            CkDec_pp.EvalSignalSmearedMain(MOM),CkDec_pp.EvalSignalSmearedChild(0,MOM),
            CkDec_pp.EvalSignalSmearedMain(MOM)+CkDec_pp.EvalSignalSmearedChild(0,MOM)+1,
            CkDec_pp.GetCk()->Eval(MOM)
          );
  }

  TH1F* hGHETTO_PS = new TH1F("hGHETTO_PS","hGHETTO_PS",NumMomBins_pp,MomBins_pp);
  TF1* fitPS = new TF1("fitPS","[2]*(x*x*exp(-pow(x/sqrt(2)/[0],[1])))/pow([0],3.)",0,1000);
  fitPS->SetParameter(0,3.54129e+02);
  fitPS->SetParameter(1,1.03198e+00);
  fitPS->SetParameter(2,2.47853e-01);
  for(unsigned uBin=0; uBin<NumMomBins_pp; uBin++){
    double MOM = hGHETTO_PS->GetBinCenter(uBin+1);
    hGHETTO_PS->SetBinContent(uBin+1,fitPS->Eval(MOM));
    //hGHETTO_PS->SetBinContent(uBin+1,1);
  }

  CkDec_pp.GetCk()->SetSourcePar(0,SourceSize);
  CkDec_pp.AddPhaseSpace(hGHETTO_PS);
  CkDec_pp.AddPhaseSpace(0,hGHETTO_PS);
  CkDec_pp.Update(true,true);

  CkDec_pL.AddPhaseSpace(hGHETTO_PS);
  CkDec_pL.Update(true,true);

  printf("lam main = %.2f\n",CkDec_pp.GetLambdaMain()*100.);
  for(unsigned uBin=0; uBin<NumMomBins_pp; uBin++){
    double MOM = h_pL_pp_old->GetBinCenter(uBin+1);
    h_pp_new->SetBinContent(uBin+1,CkDec_pp.EvalCk(MOM));
    h_pL_pp_new->SetBinContent(uBin+1,CkDec_pp.EvalSignalSmearedChild(0,MOM)+1.);

    h_pL_new->SetBinContent(uBin+1,CkDec_pL.EvalCk(MOM));

    if(MOM>60) continue;
    printf(" nECK(%.0f) = %.3f = %.3f+1 = %.3f = %.3f+%.3f+1 = %.3f (%.3f)\n",MOM,
            CkDec_pp.EvalCk(MOM),
            CkDec_pp.EvalSignalSmeared(MOM),
            CkDec_pp.EvalSignalSmeared(MOM)+1,
            CkDec_pp.EvalSignalSmearedMain(MOM),CkDec_pp.EvalSignalSmearedChild(0,MOM),
            CkDec_pp.EvalSignalSmearedMain(MOM)+CkDec_pp.EvalSignalSmearedChild(0,MOM)+1,
            CkDec_pp.GetCk()->Eval(MOM)
          );
  }

  TF1* fit_pp_data_new = new TF1("fit_pp_data_new",VerySimpleFitter_pp,0,240,3);
  fit_pp_data_new->SetParameter(0,SourceSize);
  fit_pp_data_new->SetParameter(1,1);
  fit_pp_data_new->FixParameter(2,0);
  hData_pp->Fit(fit_pp_data_new,"S, N, R, M");
  CkDec_pp.GetCk()->SetSourcePar(0,SourceSize);

  TF1* fit_pp_Femto = new TF1("fit_pp_Femto",VerySimpleFitter_pp,0,300,3);
  fit_pp_Femto->SetParameter(0,SourceSize);
  fit_pp_Femto->FixParameter(1,1);
  fit_pp_Femto->FixParameter(2,0);
  h_pp_old->Fit(fit_pp_Femto,"S, N, R, M");
  CkDec_pp.GetCk()->SetSourcePar(0,SourceSize);
  printf("True value: r = %.3f fm\n", SourceSize);
  printf(" Fit value: r = %.3f fm\n", fit_pp_Femto->GetParameter(0));
  printf(" Difference: %.1f%%\n",fabs(SourceSize-fit_pp_Femto->GetParameter(0))/(SourceSize)*100.);

  TH1F* h_pp_ratio;
  h_pp_ratio = (TH1F*)h_pp_new->Clone("h_pp_ratio");
  h_pp_ratio->Divide(h_pp_old);
  TH1F* h_pp_diff;
  h_pp_diff = (TH1F*)h_pp_new->Clone("h_pp_diff");
  h_pp_diff->Add(h_pp_old,-1.);
  TH1F* h_pL_pp_diff;
  h_pL_pp_diff = (TH1F*)h_pL_pp_new->Clone("h_pL_pp_diff");
  h_pL_pp_diff->Add(h_pL_pp_old,-1.);

  TGraphErrors* g_pp_ratio_1s = new TGraphErrors();
  g_pp_ratio_1s->SetName("g_pp_ratio_1s");
  g_pp_ratio_1s->SetFillColorAlpha(kBlue,0.4);
  g_pp_ratio_1s->SetLineColor(kBlue);
  g_pp_ratio_1s->SetLineWidth(0);
  TGraphErrors* g_pp_ratio_3s = new TGraphErrors();
  g_pp_ratio_3s->SetName("g_pp_ratio_3s");
  g_pp_ratio_3s->SetFillColorAlpha(kRed,0.4);
  g_pp_ratio_3s->SetLineColor(kRed);
  g_pp_ratio_3s->SetLineWidth(0);
  TGraphErrors* g_pp_ratio_3smT = new TGraphErrors();
  g_pp_ratio_3smT->SetName("g_pp_ratio_3smT");
  g_pp_ratio_3smT->SetFillColorAlpha(kGreen,0.4);
  g_pp_ratio_3smT->SetLineColor(kGreen);
  g_pp_ratio_3smT->SetLineWidth(0);
  for(unsigned uBin=0; uBin<h_pp_ratio->GetNbinsX(); uBin++){
    //h_pp_ratio->SetBinError(uBin+1,0);
    double RelErr = hData_pp->GetBinError(uBin+1)/hData_pp->GetBinContent(uBin+1);
    h_pp_ratio->SetBinError(uBin+1,h_pp_ratio->GetBinContent(uBin+1)*RelErr);
    g_pp_ratio_1s->SetPoint(uBin,h_pp_ratio->GetBinCenter(uBin+1),h_pp_ratio->GetBinContent(uBin+1));
    g_pp_ratio_1s->SetPointError(uBin,0,h_pp_ratio->GetBinError(uBin+1));
    g_pp_ratio_3s->SetPoint(uBin,h_pp_ratio->GetBinCenter(uBin+1),h_pp_ratio->GetBinContent(uBin+1));
    g_pp_ratio_3s->SetPointError(uBin,0,3.*h_pp_ratio->GetBinError(uBin+1));
    g_pp_ratio_3smT->SetPoint(uBin,h_pp_ratio->GetBinCenter(uBin+1),h_pp_ratio->GetBinContent(uBin+1));
    g_pp_ratio_3smT->SetPointError(uBin,0,3.*h_pp_ratio->GetBinError(uBin+1)*sqrt(7.));//7 mT bins of equal yields
  }

  TCanvas* can1 = new TCanvas("can1", "can1", 1);
  TH1F* h_pp_ratio_axis;
  h_pp_ratio_axis = (TH1F*)h_pp_ratio->Clone("h_pp_ratio_axis");
  for(unsigned uBin=0; uBin<h_pp_ratio->GetNbinsX(); uBin++){
    h_pp_ratio_axis->SetBinContent(uBin+1,1);
    h_pp_ratio_axis->SetBinError(uBin+1,0);
  }
  h_pp_ratio_axis->SetLineWidth(3);
  h_pp_ratio_axis->SetLineColor(kBlack);
  h_pp_ratio_axis->GetYaxis()->SetRangeUser(0.92,1.05);
  h_pp_ratio_axis->SetStats(false);
  h_pp_ratio_axis->SetTitle("; #it{k*} (MeV/#it{c}); #it{C}(#it{k*})_{new}/#it{C}(#it{k*})_{old}");
  h_pp_ratio_axis->GetXaxis()->SetTitleSize(0.05);
  h_pp_ratio_axis->GetXaxis()->SetLabelSize(0.05);
  h_pp_ratio_axis->GetXaxis()->SetTitleOffset(1.10);
  h_pp_ratio_axis->GetYaxis()->SetTitleSize(0.05);
  h_pp_ratio_axis->GetYaxis()->SetLabelSize(0.05);
  h_pp_ratio_axis->GetYaxis()->SetTitleOffset(1.10);

  can1->cd(0); can1->SetCanvasSize(1920/3*2, 1080/3*2); can1->SetMargin(0.12,0.05,0.14,0.05);//lrbt
  TLegend *legend = new TLegend(0.6,0.2,0.9,0.45);//lbrt
  legend->SetBorderSize(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->AddEntry(g_pp_ratio_1s, "1#sigma error (data, m_{T} int.)", "fpe");
  legend->AddEntry(g_pp_ratio_3s, "3#sigma error (data, m_{T} int.)", "fpe");
  legend->AddEntry(g_pp_ratio_3smT, "3#sigma error (data, single m_{T} bin)", "fpe");

  h_pp_ratio_axis->Draw();
  g_pp_ratio_3smT->Draw("3,same");
  g_pp_ratio_3s->Draw("3,same");
  g_pp_ratio_1s->Draw("3,same");
  legend->Draw("same");
  can1->SaveAs(OutputFolder+"can1.pdf");

  fit_pp_data_old->Write();
  fit_pp_data_new->Write();
  hGHETTO_PS->Write();
  h_pL_pp_old->Write();
  h_pL_pp_new->Write();
  h_pL_pp_diff->Write();
  h_pp_old->Write();
  fit_pp_Femto->Write();
  h_pp_new->Write();
  h_pp_diff->Write();
  h_pp_ratio->Write();
  g_pp_ratio_1s->Write();
  g_pp_ratio_3s->Write();
  h_pL_old->Write();
  h_pL_new->Write();
  can1->Write();

  delete Ck_pp;
  delete Ck_pL;
  delete h_pL_pp_old;
  delete h_pL_pp_new;
  delete h_pp_old;
  delete h_pp_new;
  delete h_pp_ratio;
  delete fit_pp_Femto;
  delete h_pL_old;
  delete h_pL_new;
  delete hGHETTO_PS;
  delete OutputFile;
}

void Test_new_fold(){
  const double SourceSize = 1.1;
  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles/",GetCernBoxDimi()));
  TH2F* hResolution_pp = AnalysisObject.GetResolutionMatrix("pp13TeV_HM_DimiJun20","pLambda");//AnalysisObject.GetResolutionMatrix("pp13TeV_HM_Dec19","pp");

  printf("Bin width: %.0f\n",hResolution_pp->GetXaxis()->GetBinWidth(1));
  hResolution_pp->Rebin2D(2,2);
  //DLM_Histo<float>* dlmResolution_pp = Convert_TH2F_DlmHisto(hResolution_pp);

  //TH2F* hResidual_pp_pL = AnalysisObject.GetResidualMatrix("pp","pLambda");
  //TH1F* hData_pp = AnalysisObject.GetAliceExpCorrFun(DataSample,"pp","_0",0,false,-1);
  double lam_pp[5];
  //AnalysisObject.SetUpLambdaPars_pp("pp13TeV_HM_Dec19",0,lam_pp);

  TString OutputFolder = TString::Format("%s/OtherTasks/Test_new_fold/",GetFemtoOutputFolder());
  TString OutFileName = "fOutput.root";
  TFile* OutputFile = new TFile(OutputFolder+OutFileName,"recreate");
  //hData_pp->Write();

  TString SourceDescription = "Gauss";
  //TString SourceDescription = "McGauss_ResoTM";
  //TString SourceDescription = "McLevy_ResoTM";

  unsigned NumMomBins_pp = 75;
  double kMin_pp = 0;
  double kMax_pp = 300;
  CATS AB_pp;
  DLM_Ck* Ck_pp;
  AB_pp.SetMomBins(NumMomBins_pp,kMin_pp,kMax_pp);
  AB_pp.SetNotifications(CATS::nWarning);
  AnalysisObject.SetUpCats_pp(AB_pp,"AV18",SourceDescription,0,202);//McLevyNolan_Reso

  AB_pp.SetAnaSource(0,SourceSize);
  AB_pp.KillTheCat();
  printf("Killed\n");
  //DLM_Histo<float> dlm_pp;
  //dlm_pp.SetUp(1);
  //dlm_pp.SetUp(0,NumMomBins_pp,kMin_pp,kMax_pp);
  //dlm_pp.Initialize();
  //printf("Initialized\n");
  //for(unsigned uBin=0; uBin<NumMomBins_pp; uBin++){
  //  dlm_pp.SetBinContent(uBin,AB_pp.GetCorrFun(uBin));
  //}
  //printf("Filled\n");
  TH1F* h_pp = new TH1F("h_pp","h_pp",NumMomBins_pp,kMin_pp,kMax_pp);
  for(unsigned uBin=0; uBin<NumMomBins_pp; uBin++){
    h_pp->SetBinContent(uBin+1,AB_pp.GetCorrFun(uBin));
  }
  DLM_Unfold dlmFold;
  dlmFold.SetData(h_pp);
  dlmFold.SetResponse(hResolution_pp);
  dlmFold.SetSilentMode(false);
  //printf("All set\n");
  //OutputFile->cd();
  DLM_Histo<float>* dlm_pp_fold = dlmFold.Fold();
  printf("Folded\n");
  OutputFile->cd();
  TH1F* h_pp_fold = Convert_DlmHisto_TH1F(dlm_pp_fold,"h_pp_fold");
  printf("Converted\n");


  //test to see what my integral function does if applied to rebin
  DLM_Histo<float> dlm_pp_fold_reb;//(*dlm_pp_fold);

  double* BinRange = dlm_pp_fold->GetBinRange(0);
  double* BinRangeNew = new double[dlm_pp_fold->GetNbins()];
  unsigned BinCounter = 0;
  for(unsigned uBin=0; uBin<dlm_pp_fold->GetNbins(); uBin++){
    //printf("BinRange[%u] = %.2f\n",uBin,BinRange[uBin]);
    if(uBin%2==0){
      BinRangeNew[BinCounter] = BinRange[uBin];
      //printf("  BinRangeNew[%u] = %.2f\n",BinCounter,BinRangeNew[BinCounter]);
      BinCounter++;
    }
    BinRangeNew[BinCounter] = dlm_pp_fold->GetUpEdge(0);
  }
  dlm_pp_fold_reb.SetUp(1);
  dlm_pp_fold_reb.SetUp(0,BinCounter,BinRangeNew);
  dlm_pp_fold_reb.Initialize();

  DLM_Histo<float> dlm_pp_fold_reb_2(dlm_pp_fold_reb);

  delete [] BinRangeNew;
  for(unsigned uBin=0; uBin<dlm_pp_fold_reb.GetNbins(); uBin++){
  //for(unsigned uBin=0; uBin<1; uBin++){
    double xAxisLow = dlm_pp_fold_reb.GetBinLowEdge(0,uBin);
    double xAxisUp = dlm_pp_fold_reb.GetBinUpEdge(0,uBin);
    dlm_pp_fold_reb.SetBinContent(uBin,dlm_pp_fold->Integral(&xAxisLow,&xAxisUp,true));
    //printf("%u : %f\n",uBin,dlm_pp_fold->Integral(&xAxisLow,&xAxisUp));
  }
  TH1F* hlm_pp_fold_reb = Convert_DlmHisto_TH1F(&dlm_pp_fold_reb,"hlm_pp_fold_reb");

  //!!N.B. see below
  dlm_pp_fold->Rebin(dlm_pp_fold_reb_2,true);
  TH1F* hlm_pp_fold_reb_2 = Convert_DlmHisto_TH1F(&dlm_pp_fold_reb_2,"hlm_pp_fold_reb_2");

  DLM_Histo<float> dlm_pp_fold_reb_3(*dlm_pp_fold);
  unsigned REBIN_3 = 2;
  printf("before: %u\n",dlm_pp_fold_reb_3.GetNbins());
  //N.B. there might be an issue with the scaling of this histo, or dlm_pp_fold
  //there were some changes with Rebin in DLM_Histo, if problems, try to
  //make above ,false
  dlm_pp_fold_reb_3.Rebin(&REBIN_3);
  printf("after: %u\n",dlm_pp_fold_reb_3.GetNbins());
  TH1F* hlm_pp_fold_reb_3 = Convert_DlmHisto_TH1F(&dlm_pp_fold_reb_3,"hlm_pp_fold_reb_3");
//return;
  printf("The old stuff\n");
  Ck_pp = new DLM_Ck(AB_pp.GetNumSourcePars(),0,AB_pp);
  Ck_pp->Update();

  DLM_CkDecomposition CkDec_pp("pp",0,*Ck_pp,hResolution_pp);
  CkDec_pp.Update();

  printf("Done\n");

  TH1F* h_pp_cats = new TH1F("h_pp_cats","h_pp_cats",NumMomBins_pp,kMin_pp,kMax_pp);
  for(unsigned uBin=0; uBin<NumMomBins_pp; uBin++){
    h_pp_cats->SetBinContent(uBin+1,CkDec_pp.EvalCk(h_pp_cats->GetBinCenter(uBin+1)));
  }
  OutputFile->cd();
  h_pp->Write();
  h_pp_cats->Write();
  h_pp_fold->Write();
  hlm_pp_fold_reb->Write();
  hlm_pp_fold_reb_2->Write();
  hlm_pp_fold_reb_3->Write();

  printf("Written\n");
  //UNFOLD

  //OutputFile->cd();
  DLM_Unfold dlmUnfold;
  for(unsigned uBin=0; uBin<NumMomBins_pp; uBin++){
    h_pp_fold->SetBinError(uBin+1,0.5);
  }
  dlmUnfold.SetData(h_pp_fold);
  dlmUnfold.SetResponse(hResolution_pp);
  dlmUnfold.SetSilentMode(false);
  dlmUnfold.SetUnfoldPrecision(0.001,0.5);

  DLM_Histo<float>*  dlm_pp_unfolded = dlmUnfold.Unfold();
  printf("Unfolded\n");
  TH1F* h_pp_unfolded = Convert_DlmHisto_TH1F(dlm_pp_unfolded,"h_pp_unfolded");
  for(unsigned uBin=0; uBin<NumMomBins_pp; uBin++){
    if(!h_pp_unfolded) break;
    h_pp_unfolded->SetBinError(uBin+1,0.1);
  }
  printf("Converted\n");

  OutputFile->cd();
  if(h_pp_unfolded)h_pp_unfolded->Write();
  //dlmUnfold.GetFitFoldOriginal()->Write();
  //dlmUnfold.GetFitFoldFinal()->Write();
  //dlmUnfold.GetFitUnfoldFinal()->Write();
  printf("Written\n");

  delete h_pp;
  delete h_pp_cats;
  delete h_pp_fold;
  printf("Del 1\n");
  if(h_pp_unfolded)delete h_pp_unfolded;
  if(dlm_pp_unfolded)delete dlm_pp_unfolded;
  printf("Del 2\n");
  delete Ck_pp;
  delete hResolution_pp;
  printf("Del 3\n");
  //OutputFile->Close();
  delete OutputFile;

}


//we smear pp with some dummy ME and SE
void TestUnfoldSEME(){

  const double SourceSize = 1.2;
  const unsigned NumSE = 100*1000*1000;
  const unsigned NumME = 100*NumSE;
  const double Norm = 1.01;
  const double lambda = 0.67/0.99;//corrected for purity
  const double lambda_pL = 0.20/0.99;
  const bool Flat_pL = true;
  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles/",GetCernBoxDimi()));
  //TH2F* hResolution_pp = AnalysisObject.GetResolutionMatrix("pp13TeV_HM_Dec19","pp");
  TH2F* hResolution_pp = AnalysisObject.GetResolutionMatrix("pp13TeV_HM_DimiJun20","pp");
  TH2F* hResidual_pp_pL = AnalysisObject.GetResidualMatrix("pp","pLambda");
  hResolution_pp->Rebin2D(4,4);
  hResidual_pp_pL->Rebin2D(4,4);
  printf("Bin width: %.0f (%.0f)\n",hResolution_pp->GetXaxis()->GetBinWidth(1),hResidual_pp_pL->GetXaxis()->GetBinWidth(1));

  double lam_pp[5];

  TString OutputFolder = TString::Format("%s/OtherTasks/Test_new_fold/",GetFemtoOutputFolder());
  TString OutFileName = "TestUnfoldSEME.root";
  TFile* OutputFile = new TFile(OutputFolder+OutFileName,"recreate");

  TString SourceDescription = "Gauss";

  unsigned NumMomBins_pp = 80;
  double kMin_pp = 0;
  double kMax_pp = 320;
  CATS AB_pp;
  AB_pp.SetMomBins(NumMomBins_pp,kMin_pp,kMax_pp);
  AB_pp.SetNotifications(CATS::nWarning);
  AnalysisObject.SetUpCats_pp(AB_pp,"AV18",SourceDescription,0,202);//McLevyNolan_Reso
  AB_pp.SetAnaSource(0,SourceSize);
  AB_pp.KillTheCat();
  printf("Killed\n");

  CATS AB_pL;
  AB_pL.SetMomBins(NumMomBins_pp,kMin_pp,kMax_pp);
  AB_pL.SetNotifications(CATS::nWarning);
  AnalysisObject.SetUpCats_pL(AB_pL,"Chiral_Coupled_SPD",SourceDescription,0,202);//McLevyNolan_Reso
  AB_pL.SetAnaSource(0,SourceSize);
  AB_pL.KillTheCat();
  printf("Killed\n");

  DLM_Ck Ck_pp(AB_pp.GetNumSourcePars(),0,AB_pp,200,0,800);
  Ck_pp.SetCutOff(kMax_pp,kMax_pp*2);
  Ck_pp.Update();


  DLM_Ck CkShort_pp(AB_pp.GetNumSourcePars(),0,AB_pp,200,0,800);
  CkShort_pp.SetCutOff(kMax_pp,kMax_pp*2);
  CkShort_pp.Update();

  DLM_Ck CkShort_pL(AB_pL.GetNumSourcePars(),0,AB_pL,200,0,800);
  CkShort_pL.SetCutOff(kMax_pp,kMax_pp*2*2);
  CkShort_pL.Update();

  //DLM_CkDecomposition CkDec_pp_Feed("pp",2,CkShort_pp,hResolution_pp);
  DLM_CkDecomposition CkDec_pL("pLambda",1,CkShort_pL,NULL);
  //CkDec_pp_Feed.AddContribution(0,lambda_pL,DLM_CkDecomposition::cFeedDown,&CkDec_pL,hResidual_pp_pL);
  //CkDec_pp_Feed.AddContribution(1,1.-lambda-lambda_pL,DLM_CkDecomposition::cFeedDown);
  CkDec_pL.AddContribution(0,0.5,DLM_CkDecomposition::cFeedDown);
  //CkDec_pp_Feed.Update();

  TH1F* h_pp_ME_true = new TH1F("h_pp_ME_true","h_pp_ME_true",200,0,800);
  TH1F* h_pp_SE_true = new TH1F("h_pp_SE_true","h_pp_SE_true",200,0,800);
  TH1F* h_pp_Ck_true = new TH1F("h_pp_Ck_true","h_pp_Ck_true",200,0,800);
  TF1* fitPS = new TF1("fitPS","[2]*(x*x*exp(-pow(x/sqrt(2)/[0],[1])))/pow([0],3.)",0,800);
  //i fitted ME of PL with this functions, these were the parameters.
  //we allow -10 to +20% of variations (+20 we expect from ME of pXi etc..., -10 for the fun of it)
  //fitPS->SetParameter(0,3.54129e+02);
  fitPS->SetParameter(0,3.54129e+02);
  fitPS->SetParameter(1,1.03198e+00);
  fitPS->SetParameter(2,2.47853e-01);
  for(unsigned uBin=0; uBin<h_pp_ME_true->GetNbinsX(); uBin++){
    double MOM = h_pp_ME_true->GetBinCenter(uBin+1);
  //  double VAL = 5.40165e-08-1.06309e-08*pow(MOM,1.)+5.71623e-09*pow(MOM,2.)-
  //  -1.15454e-11*pow(MOM,3.)+1.17292e-14*pow(MOM,4.)-7.94266e-18*pow(MOM,5.)+
  //  3.70109e-21*pow(MOM,6.)-1.10644e-24*pow(MOM,7.)+1.86899e-28*pow(MOM,8.)-1.34242e-32*pow(MOM,9.);
    double VAL = fitPS->Eval(MOM);
    h_pp_ME_true->SetBinContent(uBin+1,VAL);
    h_pp_ME_true->SetBinError(uBin+1,0);
    h_pp_SE_true->SetBinContent(uBin+1,VAL*Ck_pp.Eval(MOM));
    h_pp_SE_true->SetBinError(uBin+1,0);
    h_pp_Ck_true->SetBinContent(uBin+1,Ck_pp.Eval(MOM));
    h_pp_Ck_true->SetBinError(uBin+1,0);
  }
  h_pp_SE_true->Scale(1./h_pp_SE_true->Integral());
  h_pp_ME_true->Scale(1./h_pp_ME_true->Integral());

  //fold the correlation,SE and ME
  DLM_Unfold dlmFold;
  dlmFold.SetResponse(hResolution_pp);
  dlmFold.SetSilentMode(false);

  dlmFold.SetData(h_pp_SE_true);
  DLM_Histo<float>* dlm_pp_SE_exp = dlmFold.Fold();
  TH1F* h_pp_SE_exp = Convert_DlmHisto_TH1F(dlm_pp_SE_exp,"h_pp_SE_exp");

  dlmFold.SetData(h_pp_ME_true);
  DLM_Histo<float>* dlm_pp_ME_exp = dlmFold.Fold();
  TH1F* h_pp_ME_exp = Convert_DlmHisto_TH1F(dlm_pp_ME_exp,"h_pp_ME_exp");

  DLM_Histo<float>* dlm_pp_Ck_exp(dlm_pp_SE_exp);
  dlm_pp_Ck_exp[0] /= dlm_pp_ME_exp[0];
  TH1F* h_pp_Ck_exp = Convert_DlmHisto_TH1F(dlm_pp_Ck_exp,"h_pp_Ck_exp");

  dlmFold.SetData(h_pp_Ck_true);
  DLM_Histo<float>* dlm_pp_Ck_truefold = dlmFold.Fold();
  TH1F* h_pp_Ck_truefold = Convert_DlmHisto_TH1F(dlm_pp_Ck_truefold,"h_pp_Ck_truefold");

  ////DLM_CkDecomposition CkDec_pp_Feed("pp",2,CkShort_pp,hResolution_pp);
  //DLM_CkDecomposition CkDec_pL("pLambda",1,CkShort_pL,NULL);
  ////CkDec_pp_Feed.AddContribution(0,lambda_pL,DLM_CkDecomposition::cFeedDown,&CkDec_pL,hResidual_pp_pL);
  ////CkDec_pp_Feed.AddContribution(1,1.-lambda-lambda_pL,DLM_CkDecomposition::cFeedDown);
  //CkDec_pL.AddContribution(0,0.5,DLM_CkDecomposition::cFeedDown);


  //from CATS
  DLM_CkDecomposition CkDec_pp_direct("pp",2,CkShort_pp,hResolution_pp);
  if(Flat_pL) CkDec_pp_direct.AddContribution(0,lambda_pL,DLM_CkDecomposition::cFeedDown);
  else {CkDec_pp_direct.AddContribution(0,lambda_pL,DLM_CkDecomposition::cFeedDown,&CkDec_pL,hResidual_pp_pL);}
  CkDec_pp_direct.AddContribution(1,1.-lambda-lambda_pL,DLM_CkDecomposition::cFeedDown);
  CkDec_pp_direct.Update();
  DLM_CkDecomposition CkDec_pp_PS("pp",2,CkShort_pp,hResolution_pp);
  CkDec_pp_PS.AddPhaseSpace(h_pp_ME_true);
  if(Flat_pL) CkDec_pp_PS.AddContribution(0,lambda_pL,DLM_CkDecomposition::cFeedDown);
  else {CkDec_pp_PS.AddContribution(0,lambda_pL,DLM_CkDecomposition::cFeedDown,&CkDec_pL,hResidual_pp_pL); CkDec_pp_PS.AddPhaseSpace(0,h_pp_ME_true);}
  CkDec_pp_PS.AddContribution(1,1.-lambda-lambda_pL,DLM_CkDecomposition::cFeedDown);
  CkDec_pp_PS.Update();

  DLM_CkDecomposition CkDec_pp_PSexp("pp",2,CkShort_pp,hResolution_pp);
  CkDec_pp_PSexp.AddPhaseSpace(h_pp_ME_exp);
  if(Flat_pL) CkDec_pp_PSexp.AddContribution(0,lambda_pL,DLM_CkDecomposition::cFeedDown);
  else {CkDec_pp_PSexp.AddContribution(0,lambda_pL,DLM_CkDecomposition::cFeedDown,&CkDec_pL,hResidual_pp_pL); CkDec_pp_PS.AddPhaseSpace(0,h_pp_ME_exp);}
  CkDec_pp_PSexp.AddContribution(1,1.-lambda-lambda_pL,DLM_CkDecomposition::cFeedDown);
  CkDec_pp_PSexp.Update();

  TH1F* hCk_CATS_direct = new TH1F("hCk_CATS_direct","hCk_CATS_direct",NumMomBins_pp,kMin_pp,kMax_pp);
  TH1F* hCk_CATS_PS = new TH1F("hCk_CATS_PS","hCk_CATS_PS",NumMomBins_pp,kMin_pp,kMax_pp);
  TH1F* hCk_CATS_PSexp = new TH1F("hCk_CATS_PSexp","hCk_CATS_PSexp",NumMomBins_pp,kMin_pp,kMax_pp);
  //TH1F* hCk_CATS_PSexp_Feed = new TH1F("hCk_CATS_PSexp_Feed","hCk_CATS_PSexp_Feed",NumMomBins_pp,kMin_pp,kMax_pp);

  for(unsigned uBin=0; uBin<NumMomBins_pp; uBin++){
    double MOM = AB_pp.GetMomentum(uBin);
    hCk_CATS_direct->SetBinContent(uBin+1,CkDec_pp_direct.EvalCk(MOM));
    hCk_CATS_direct->SetBinError(uBin+1,0);
    hCk_CATS_PS->SetBinContent(uBin+1,CkDec_pp_PS.EvalCk(MOM));
    hCk_CATS_PS->SetBinError(uBin+1,0);
    hCk_CATS_PSexp->SetBinContent(uBin+1,CkDec_pp_PSexp.EvalCk(MOM));
    hCk_CATS_PSexp->SetBinError(uBin+1,0);
    //hCk_CATS_PSexp_Feed->SetBinContent(uBin+1,CkDec_pp_Feed.EvalCk(MOM));
    //hCk_CATS_PSexp_Feed->SetBinError(uBin+1,0);
    //printf("Ck_pp(%.0f) = %.2f\n",MOM,Ck_pp.GetBinContent(uBin));
    //printf("CkDec_pp_direct(%.0f) = %.2f\n",MOM,CkDec_pp_direct.EvalCk(MOM));
    //printf("CkDec_pp_PS(%.0f) = %.2f\n",MOM,CkDec_pp_PS.EvalCk(MOM));
    //printf("CkDec_pp_PSexp(%.0f) = %.2f\n",MOM,CkDec_pp_PSexp.EvalCk(MOM));
    //printf("\n");
  }

  //unfold the correlation,SE and ME, building the SE/ME to compare to the direct correlation unfold

  for(unsigned uBin=0; uBin<NumMomBins_pp; uBin++){
    double CkVal;
    double CkErr;

    CkVal=h_pp_Ck_true->GetBinContent(uBin+1);
    CkErr=h_pp_Ck_true->GetBinError(uBin+1);
    h_pp_Ck_true->SetBinContent(uBin+1,Norm*(lambda*CkVal+1-lambda));
    h_pp_Ck_true->SetBinError(uBin+1,Norm*lambda*CkErr);

    CkVal=h_pp_Ck_truefold->GetBinContent(uBin+1);
    CkErr=h_pp_Ck_truefold->GetBinError(uBin+1);
    h_pp_Ck_truefold->SetBinContent(uBin+1,Norm*(lambda*CkVal+1-lambda));
    h_pp_Ck_truefold->SetBinError(uBin+1,Norm*lambda*CkErr);

    CkVal=h_pp_Ck_exp->GetBinContent(uBin+1);
    CkErr=h_pp_Ck_exp->GetBinError(uBin+1);
    h_pp_Ck_exp->SetBinContent(uBin+1,Norm*(lambda*CkVal+1-lambda));
    h_pp_Ck_exp->SetBinError(uBin+1,Norm*lambda*CkErr);

    CkVal=hCk_CATS_direct->GetBinContent(uBin+1);
    CkErr=hCk_CATS_direct->GetBinError(uBin+1);
    hCk_CATS_direct->SetBinContent(uBin+1,Norm*CkVal);
    hCk_CATS_direct->SetBinError(uBin+1,Norm*CkErr);

    //CkVal=hCk_CATS_PSexp_Feed->GetBinContent(uBin+1);
    //CkErr=hCk_CATS_PSexp_Feed->GetBinError(uBin+1);
    //hCk_CATS_PSexp_Feed->SetBinContent(uBin+1,Norm*CkVal);
    //hCk_CATS_PSexp_Feed->SetBinError(uBin+1,Norm*CkErr);

    CkVal=hCk_CATS_PS->GetBinContent(uBin+1);
    CkErr=hCk_CATS_PS->GetBinError(uBin+1);
    hCk_CATS_PS->SetBinContent(uBin+1,Norm*CkVal);
    hCk_CATS_PS->SetBinError(uBin+1,Norm*CkErr);

    CkVal=hCk_CATS_PSexp->GetBinContent(uBin+1);
    CkErr=hCk_CATS_PSexp->GetBinError(uBin+1);
    hCk_CATS_PSexp->SetBinContent(uBin+1,Norm*CkVal);
    hCk_CATS_PSexp->SetBinError(uBin+1,Norm*CkErr);
  }

  OutputFile->cd();
  h_pp_SE_true->Write();
  h_pp_ME_true->Write();
  h_pp_Ck_true->Write();
  h_pp_Ck_truefold->Write();
  //h_pp_Ck_perfect->Write();

  h_pp_SE_exp->Write();
  h_pp_ME_exp->Write();
  h_pp_Ck_exp->Write();

  hCk_CATS_direct->Write();
  hCk_CATS_PS->Write();
  hCk_CATS_PSexp->Write();
  //hCk_CATS_PSexp_Feed->Write();

  //h_pp_Ck_fold->Write();


  delete h_pp_SE_true;
  delete h_pp_ME_true;
  delete h_pp_Ck_true;
  //delete h_pp_Ck_perfect;
  delete h_pp_SE_exp;
  delete h_pp_ME_exp;
  delete hCk_CATS_direct;
  delete hCk_CATS_PS;
  delete hCk_CATS_PSexp;
  //delete h_pp_Ck_fold;
  delete fitPS;
  delete OutputFile;
}




void TestUnfoldSEME_pL(){

  const double SourceSize = 1.2;
  const unsigned NumSE = 100*1000*1000;
  const unsigned NumME = 100*NumSE;
  const double Norm = 1;
  const double lambda = 0.47;
  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles/",GetCernBoxDimi()));
  TH2F* hResolution_pL = AnalysisObject.GetResolutionMatrix("pp13TeV_HM_DimiJun20","pLambda");
  hResolution_pL->Rebin2D(12,12);
  printf("BINS %u\n",hResolution_pL->GetXaxis()->GetNbins());
  TH2F* hResolution_pp = AnalysisObject.GetResolutionMatrix("pp13TeV_HM_Dec19","pp");
  printf("Bin width: %.0f\n",hResolution_pL->GetXaxis()->GetBinWidth(1));
  double lam_pL[5];

  TString OutputFolder = TString::Format("%s/OtherTasks/Test_new_fold/",GetFemtoOutputFolder());
  TString OutFileName = "TestUnfoldSEME_pL.root";
  TFile* OutputFile = new TFile(OutputFolder+OutFileName,"recreate");

  TString SourceDescription = "Gauss";

  unsigned NumMomBins_pL = 30;
  double kMin_pL = 0;
  double kMax_pL = 360;
  CATS AB_pL;
  AB_pL.SetMomBins(NumMomBins_pL,kMin_pL,kMax_pL);
  AB_pL.SetNotifications(CATS::nWarning);
  AnalysisObject.SetUpCats_pL(AB_pL,"Chiral_Coupled_SPD",SourceDescription,0,202);//McLevyNolan_Reso
  AB_pL.SetAnaSource(0,SourceSize);
  AB_pL.KillTheCat();
  printf("Killed\n");

  DLM_Ck Ck_pL(AB_pL.GetNumSourcePars(),0,AB_pL,75,0,900);
  Ck_pL.SetCutOff(kMax_pL,kMax_pL*2);
  Ck_pL.Update();

  DLM_Ck CkShort_pL(AB_pL.GetNumSourcePars(),0,AB_pL,60,0,720);
  CkShort_pL.SetCutOff(kMax_pL,kMax_pL*2);
  CkShort_pL.Update();

  printf("DLM_Ck\n");

  TH1F* h_pL_ME_true = new TH1F("h_pL_ME_true","h_pL_ME_true",75,0,900);
  TH1F* h_pL_SE_true = new TH1F("h_pL_SE_true","h_pL_SE_true",75,0,900);
  TH1F* h_pL_Ck_true = new TH1F("h_pL_Ck_true","h_pL_Ck_true",75,0,900);
  TF1* fitPS = new TF1("fitPS","[2]*(x*x*exp(-pow(x/sqrt(2)/[0],[1])))/pow([0],3.)",0,900);
  //i fitted ME of PL with this functions, these were the parameters.
  //we allow -10 to +20% of variations (+20 we expect from ME of pXi etc..., -10 for the fun of it)
  //fitPS->SetParameter(0,3.54129e+02);
  fitPS->SetParameter(0,3.54129e+02);
  fitPS->SetParameter(1,1.03198e+00);
  fitPS->SetParameter(2,2.47853e-01);
  for(unsigned uBin=0; uBin<h_pL_ME_true->GetNbinsX(); uBin++){
    double MOM = h_pL_ME_true->GetBinCenter(uBin+1);
  //  double VAL = 5.40165e-08-1.06309e-08*pow(MOM,1.)+5.71623e-09*pow(MOM,2.)-
  //  -1.15454e-11*pow(MOM,3.)+1.17292e-14*pow(MOM,4.)-7.94266e-18*pow(MOM,5.)+
  //  3.70109e-21*pow(MOM,6.)-1.10644e-24*pow(MOM,7.)+1.86899e-28*pow(MOM,8.)-1.34242e-32*pow(MOM,9.);
    double VAL = fitPS->Eval(MOM);
    h_pL_ME_true->SetBinContent(uBin+1,VAL);
    h_pL_ME_true->SetBinError(uBin+1,0);
    h_pL_SE_true->SetBinContent(uBin+1,VAL*Ck_pL.Eval(MOM));
    h_pL_SE_true->SetBinError(uBin+1,0);
    h_pL_Ck_true->SetBinContent(uBin+1,Ck_pL.Eval(MOM));
    h_pL_Ck_true->SetBinError(uBin+1,0);
  }
  h_pL_SE_true->Scale(1./h_pL_SE_true->Integral());
  h_pL_ME_true->Scale(1./h_pL_ME_true->Integral());

  printf("Folding fun\n");
  //fold the correlation,SE and ME
  DLM_Unfold dlmFold;
  dlmFold.SetResponse(hResolution_pL);
  dlmFold.SetSilentMode(false);

  dlmFold.SetData(h_pL_SE_true);
  DLM_Histo<float>* dlm_pL_SE_exp = dlmFold.Fold();
  printf("Fold dlm_pL_SE_exp\n");
  TH1F* h_pL_SE_exp = Convert_DlmHisto_TH1F(dlm_pL_SE_exp,"h_pL_SE_exp");

  dlmFold.SetData(h_pL_ME_true);
  DLM_Histo<float>* dlm_pL_ME_exp = dlmFold.Fold();
  printf("Fold dlm_pL_ME_exp\n");
  TH1F* h_pL_ME_exp = Convert_DlmHisto_TH1F(dlm_pL_ME_exp,"h_pL_ME_exp");

  DLM_Histo<float>* dlm_pL_Ck_exp(dlm_pL_SE_exp);
  dlm_pL_Ck_exp[0] /= dlm_pL_ME_exp[0];
  TH1F* h_pL_Ck_exp = Convert_DlmHisto_TH1F(dlm_pL_Ck_exp,"h_pL_Ck_exp");

  dlmFold.SetData(h_pL_Ck_true);
  DLM_Histo<float>* dlm_pL_Ck_truefold = dlmFold.Fold();
  printf("Fold dlm_pL_Ck_truefold\n");
  TH1F* h_pL_Ck_truefold = Convert_DlmHisto_TH1F(dlm_pL_Ck_truefold,"h_pL_Ck_truefold");


  //from CATS
  DLM_CkDecomposition CkDec_pL_direct("pp",0,CkShort_pL,hResolution_pL);
  CkDec_pL_direct.Update();
  DLM_CkDecomposition CkDec_pL_PS("pp",0,CkShort_pL,hResolution_pL);
  CkDec_pL_PS.AddPhaseSpace(h_pL_ME_true);
  CkDec_pL_PS.Update();

  DLM_CkDecomposition CkDec_pL_PSexp("pp",0,CkShort_pL,hResolution_pL);
  CkDec_pL_PSexp.AddPhaseSpace(h_pL_ME_exp);
  CkDec_pL_PSexp.Update();

  TH1F* hCk_CATS_direct = new TH1F("hCk_CATS_direct","hCk_CATS_direct",NumMomBins_pL,kMin_pL,kMax_pL);
  TH1F* hCk_CATS_PS = new TH1F("hCk_CATS_PS","hCk_CATS_PS",NumMomBins_pL,kMin_pL,kMax_pL);
  TH1F* hCk_CATS_PSexp = new TH1F("hCk_CATS_PSexp","hCk_CATS_PSexp",NumMomBins_pL,kMin_pL,kMax_pL);

  for(unsigned uBin=0; uBin<NumMomBins_pL; uBin++){
    double MOM = AB_pL.GetMomentum(uBin);
    hCk_CATS_direct->SetBinContent(uBin+1,CkDec_pL_direct.EvalCk(MOM));
    hCk_CATS_PS->SetBinContent(uBin+1,CkDec_pL_PS.EvalCk(MOM));
    hCk_CATS_PSexp->SetBinContent(uBin+1,CkDec_pL_PSexp.EvalCk(MOM));

    hCk_CATS_direct->SetBinError(uBin+1,0);
    hCk_CATS_PS->SetBinError(uBin+1,0);
    hCk_CATS_PSexp->SetBinError(uBin+1,0);
    //printf("Ck_pL(%.0f) = %.2f\n",MOM,Ck_pL.GetBinContent(uBin));
    //printf("CkDec_pL_direct(%.0f) = %.2f\n",MOM,CkDec_pL_direct.EvalCk(MOM));
    //printf("CkDec_pL_PS(%.0f) = %.2f\n",MOM,CkDec_pL_PS.EvalCk(MOM));
    //printf("CkDec_pL_PSexp(%.0f) = %.2f\n",MOM,CkDec_pL_PSexp.EvalCk(MOM));
    //printf("\n");
  }

  TH1F* hData_pL = AnalysisObject.GetAliceExpCorrFun("pp13TeV_HM_Dec19","pLambda","_0",2,false,-1);
  //unfold the correlation,SE and ME, building the SE/ME to compare to the direct correlation unfold

  for(unsigned uBin=0; uBin<NumMomBins_pL; uBin++){
    double CkVal;
    double CkErr;

    CkVal=h_pL_Ck_true->GetBinContent(uBin+1);
    CkErr=h_pL_Ck_true->GetBinError(uBin+1);
    h_pL_Ck_true->SetBinContent(uBin+1,Norm*(lambda*CkVal+1-lambda));
    h_pL_Ck_true->SetBinError(uBin+1,Norm*lambda*CkErr);

    CkVal=h_pL_Ck_truefold->GetBinContent(uBin+1);
    CkErr=h_pL_Ck_truefold->GetBinError(uBin+1);
    h_pL_Ck_truefold->SetBinContent(uBin+1,Norm*(lambda*CkVal+1-lambda));
    h_pL_Ck_truefold->SetBinError(uBin+1,Norm*lambda*CkErr);

    CkVal=h_pL_Ck_exp->GetBinContent(uBin+1);
    CkErr=h_pL_Ck_exp->GetBinError(uBin+1);
    h_pL_Ck_exp->SetBinContent(uBin+1,Norm*(lambda*CkVal+1-lambda));
    h_pL_Ck_exp->SetBinError(uBin+1,Norm*lambda*CkErr);

    CkVal=hCk_CATS_direct->GetBinContent(uBin+1);
    CkErr=hCk_CATS_direct->GetBinError(uBin+1);
    hCk_CATS_direct->SetBinContent(uBin+1,Norm*(lambda*CkVal+1-lambda));
    hCk_CATS_direct->SetBinError(uBin+1,Norm*lambda*CkErr);

    CkVal=hCk_CATS_PS->GetBinContent(uBin+1);
    CkErr=hCk_CATS_PS->GetBinError(uBin+1);
    hCk_CATS_PS->SetBinContent(uBin+1,Norm*(lambda*CkVal+1-lambda));
    hCk_CATS_PS->SetBinError(uBin+1,Norm*lambda*CkErr);

    CkVal=hCk_CATS_PSexp->GetBinContent(uBin+1);
    CkErr=hCk_CATS_PSexp->GetBinError(uBin+1);
    hCk_CATS_PSexp->SetBinContent(uBin+1,Norm*(lambda*CkVal+1-lambda));
    hCk_CATS_PSexp->SetBinError(uBin+1,Norm*lambda*CkErr);
  }

  TH1F* h_pL_Ck_exp_SmallRange = (TH1F*)hCk_CATS_direct->Clone("hCkRatio_CATS_direct");
  for(unsigned uBin=0; uBin<NumMomBins_pL; uBin++){
    h_pL_Ck_exp_SmallRange->SetBinContent(uBin+1,h_pL_Ck_exp->GetBinContent(uBin+1));
    h_pL_Ck_exp_SmallRange->SetBinError(uBin+1,0);
  }

  TH1F* hCkRatio_CATS_direct = (TH1F*)hCk_CATS_direct->Clone("hCkRatio_CATS_direct");
  hCkRatio_CATS_direct->Divide(h_pL_Ck_exp_SmallRange);
  TH1F* hCkRatio_CATS_PS = (TH1F*)hCk_CATS_PS->Clone("hCkRatio_CATS_PS");
  hCkRatio_CATS_PS->Divide(h_pL_Ck_exp_SmallRange);
  TH1F* hCkRatio_CATS_PSexp = (TH1F*)hCk_CATS_PSexp->Clone("hCkRatio_CATS_PSexp");
  hCkRatio_CATS_PSexp->Divide(h_pL_Ck_exp_SmallRange);

  TH1F* hUpperErrors1 = new TH1F("hUpperErrors1","hUpperErrors1",NumMomBins_pL,kMin_pL,kMax_pL);
  TH1F* hLowerErrors1 = new TH1F("hLowerErrors1","hLowerErrors1",NumMomBins_pL,kMin_pL,kMax_pL);
  TH1F* hUpperErrors3 = new TH1F("hUpperErrors3","hUpperErrors3",NumMomBins_pL,kMin_pL,kMax_pL);
  TH1F* hLowerErrors3 = new TH1F("hLowerErrors3","hLowerErrors3",NumMomBins_pL,kMin_pL,kMax_pL);
  TGraph gUpperErrors1;
  gUpperErrors1.SetName("gUpperErrors1");
  TGraph gLowerErrors1;
  gLowerErrors1.SetName("gLowerErrors1");
  TGraph gUpperErrors3;
  gUpperErrors3.SetName("gUpperErrors3");
  TGraph gLowerErrors3;
  gLowerErrors3.SetName("gLowerErrors3");


  hUpperErrors1->SetLineColor(kGreen+1);
  hUpperErrors1->SetLineWidth(3);
  hUpperErrors1->SetLineStyle(1);
  hLowerErrors1->SetLineColor(kGreen+1);
  hLowerErrors1->SetLineWidth(3);
  hLowerErrors1->SetLineStyle(1);

  hUpperErrors3->SetLineColor(kRed);
  hUpperErrors3->SetLineWidth(3);
  hUpperErrors3->SetLineStyle(1);
  hLowerErrors3->SetLineColor(kRed);
  hLowerErrors3->SetLineWidth(3);
  hLowerErrors3->SetLineStyle(1);

  gUpperErrors1.SetLineColor(kGreen+1);
  gUpperErrors1.SetLineWidth(3);
  gUpperErrors1.SetLineStyle(1);
  gLowerErrors1.SetLineColor(kGreen+1);
  gLowerErrors1.SetLineWidth(3);
  gLowerErrors1.SetLineStyle(1);

  gUpperErrors3.SetLineColor(kRed);
  gUpperErrors3.SetLineWidth(3);
  gUpperErrors3.SetLineStyle(1);
  gLowerErrors3.SetLineColor(kRed);
  gLowerErrors3.SetLineWidth(3);
  gLowerErrors3.SetLineStyle(1);

  h_pL_Ck_exp->SetLineColor(kBlack);
  h_pL_Ck_exp->SetLineWidth(7);
  h_pL_Ck_exp->SetLineStyle(1);

  hCk_CATS_direct->SetLineColor(kRed-3);
  hCk_CATS_direct->SetLineWidth(6);
  hCk_CATS_direct->SetLineStyle(2);
  hCkRatio_CATS_direct->SetLineColor(kRed-3);
  hCkRatio_CATS_direct->SetLineWidth(6);
  hCkRatio_CATS_direct->SetLineStyle(2);

  hCk_CATS_PSexp->SetLineColor(kViolet);
  hCk_CATS_PSexp->SetLineWidth(6);
  hCk_CATS_PSexp->SetLineStyle(3);
  hCkRatio_CATS_PSexp->SetLineColor(kViolet);
  hCkRatio_CATS_PSexp->SetLineWidth(6);
  hCkRatio_CATS_PSexp->SetLineStyle(5);

  hCk_CATS_PS->SetLineColor(kAzure);
  hCk_CATS_PS->SetLineWidth(5);
  hCk_CATS_PS->SetLineStyle(1);
  hCkRatio_CATS_PS->SetLineColor(kAzure);
  hCkRatio_CATS_PS->SetLineWidth(5);
  hCkRatio_CATS_PS->SetLineStyle(1);

  for(unsigned uBin=0; uBin<NumMomBins_pL; uBin++){
    double RelErr = hData_pL->GetBinError(uBin+1)/hData_pL->GetBinContent(uBin+1);
    double MOM = hData_pL->GetBinCenter(uBin+1);
    hUpperErrors1->SetBinContent(uBin+1,1.*(1.+RelErr));
    hLowerErrors1->SetBinContent(uBin+1,1./(1.+RelErr));
    hUpperErrors3->SetBinContent(uBin+1,1.*(1.+3.*RelErr));
    hLowerErrors3->SetBinContent(uBin+1,1./(1.+3.*RelErr));

    hUpperErrors1->SetBinError(uBin+1,0);
    hLowerErrors1->SetBinError(uBin+1,0);
    hUpperErrors3->SetBinError(uBin+1,0);
    hLowerErrors3->SetBinError(uBin+1,0);

    gUpperErrors1.SetPoint(uBin,MOM,1.*(1.+RelErr));
    gLowerErrors1.SetPoint(uBin,MOM,1./(1.+RelErr));
    gUpperErrors3.SetPoint(uBin,MOM,1.*(1.+3.*RelErr));
    gLowerErrors3.SetPoint(uBin,MOM,1./(1.+3.*RelErr));
  }

  //for the plots, plug in error h_pL_Ck_exp
  for(unsigned uBin=0; uBin<NumMomBins_pL; uBin++){
    h_pL_Ck_exp->SetBinError(uBin+1,hData_pL->GetBinError(uBin+1));
  }

  OutputFile->cd();
  hData_pL->Write();
  h_pL_SE_true->Write();
  h_pL_ME_true->Write();
  h_pL_Ck_true->Write();
  h_pL_Ck_truefold->Write();
  //h_pL_Ck_perfect->Write();

  h_pL_SE_exp->Write();
  h_pL_ME_exp->Write();
  h_pL_Ck_exp->Write();//the true one

  hCk_CATS_direct->Write();
  hCk_CATS_PS->Write();
  hCk_CATS_PSexp->Write();

  hCkRatio_CATS_direct->Write();
  hCkRatio_CATS_PS->Write();
  hCkRatio_CATS_PSexp->Write();

  hUpperErrors1->Write();
  hLowerErrors1->Write();
  hUpperErrors3->Write();
  hLowerErrors3->Write();

  gUpperErrors1.Write();
  gLowerErrors1.Write();
  gUpperErrors3.Write();
  gLowerErrors3.Write();
  //h_pL_Ck_fold->Write();


  delete h_pL_SE_true;
  delete h_pL_ME_true;
  delete h_pL_Ck_true;
  //delete h_pL_Ck_perfect;
  delete h_pL_SE_exp;
  delete h_pL_ME_exp;
  delete hCk_CATS_direct;
  delete hCk_CATS_PS;
  delete hCk_CATS_PSexp;
  //delete h_pL_Ck_fold;
  delete fitPS;
  delete hData_pL;
  delete OutputFile;
}


//Conclusion on pXi:
//the first bin might change by 10-20% which is on the edge of the uncertainty
//due to the bin shift which was used, actually there is some mimicing of reducing the phase space at low momenta,
//so in fact the effect will be lessened. Given the current binning of the momentum smearing matrix etc it is diffictult to
//really quantify the effect with high precision unless more time is invested.
void TestUnfoldSEME_pXi(){

  const int FINER_BINNING = 5;
  const double SourceSize = 1.43;
  const unsigned NumSE = 100*1000*1000;
  const unsigned NumME = 100*NumSE;
  const double Norm = 1;
  const double lambda_misid = 0.085+0.029;
  const double lambda_gen = 0.513;
  const double lambda_xi1530 = 0.082;
  const double lambda_flat = 0.291;

  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles/",GetCernBoxDimi()));

  //TH2F* hResolution_pXi = AnalysisObject.GetResolutionMatrix("pPb5TeV_Run2paper","pXim");
  //hResolution_pXi->Rebin2D(2,2);

  TH2F* hResolution_pXi = AnalysisObject.GetResolutionMatrix("pp13TeV_HM_DimiJun20","pLambda");
  hResolution_pXi->Rebin2D(20/FINER_BINNING,20/FINER_BINNING);

  TH2F* hResidual_pXi_pXi1530 = AnalysisObject.GetResidualMatrix("pXim","pXim1530");
  hResidual_pXi_pXi1530->Rebin2D(20/FINER_BINNING,20/FINER_BINNING);
  printf("Bin width: %.0f (%.0f)\n",hResolution_pXi->GetXaxis()->GetBinWidth(1),hResidual_pXi_pXi1530->GetXaxis()->GetBinWidth(1));
  //double lam_pL[5];
//return;
  TString OutputFolder = TString::Format("%s/OtherTasks/Test_new_fold/",GetFemtoOutputFolder());
  TString OutFileName = "TestUnfoldSEME_pXi.root";
  TFile* OutputFile = new TFile(OutputFolder+OutFileName,"recreate");

  TString SourceDescription = "Gauss";

  const bool SHIFTED = true;


  unsigned NumMomBins_pXi = 18*FINER_BINNING;
  double kMin_pXi = 0+(SHIFTED?8:0);
  double kMax_pXi = 360+(SHIFTED?8:0);
  CATS AB_pXi;
  AB_pXi.SetMomBins(NumMomBins_pXi,kMin_pXi,kMax_pXi);
  AB_pXi.SetNotifications(CATS::nWarning);
  AnalysisObject.SetUpCats_pXim(AB_pXi,"pXim_Lattice",SourceDescription);//McLevyNolan_Reso
  AB_pXi.SetAnaSource(0,SourceSize);
  AB_pXi.KillTheCat();
  printf("Killed\n");

  CATS AB_pXi1530;
  AB_pXi1530.SetMomBins(NumMomBins_pXi,kMin_pXi,kMax_pXi);
  AB_pXi1530.SetNotifications(CATS::nWarning);
  AnalysisObject.SetUpCats_pXim(AB_pXi1530,"pXim1530",SourceDescription);//McLevyNolan_Reso
  AB_pXi1530.SetAnaSource(0,SourceSize);
  AB_pXi1530.KillTheCat();
  printf("Killed\n");

  DLM_Ck Ck_pXi(AB_pXi.GetNumSourcePars(),0,AB_pXi,40*FINER_BINNING,0+(SHIFTED?8:0),800+(SHIFTED?8:0));
  Ck_pXi.SetCutOff(kMax_pXi,kMax_pXi*2);
  Ck_pXi.Update();

  DLM_Ck Ck_pXi1530(AB_pXi1530.GetNumSourcePars(),0,AB_pXi1530,40*FINER_BINNING,0+(SHIFTED?8:0),800+(SHIFTED?8:0));
  Ck_pXi1530.SetCutOff(kMax_pXi,kMax_pXi*2);
  Ck_pXi1530.Update();

  printf("DLM_Ck\n");

  TH1F* h_pXi_ME_true = new TH1F("h_pXi_ME_true","h_pXi_ME_true",40*FINER_BINNING,0+(SHIFTED?8:0),800+(SHIFTED?8:0));
  TH1F* h_pXi_SE_true = new TH1F("h_pXi_SE_true","h_pXi_SE_true",40*FINER_BINNING,0+(SHIFTED?8:0),800+(SHIFTED?8:0));
  TH1F* h_pXi_Ck_true = new TH1F("h_pXi_Ck_true","h_pXi_Ck_true",40*FINER_BINNING,0+(SHIFTED?8:0),800+(SHIFTED?8:0));
  TF1* fitPS = new TF1("fitPS","[2]*(x*x*exp(-pow(x/sqrt(2)/[0],[1])))/pow([0],3.)",0+(SHIFTED?8:0),800+(SHIFTED?8:0));

  //i fitted ME of PXI
  fitPS->SetParameter(0,4.51e+02);
  fitPS->SetParameter(1,1.597+00);
  fitPS->SetParameter(2,1);
  for(unsigned uBin=0; uBin<h_pXi_ME_true->GetNbinsX(); uBin++){
    double MOM = h_pXi_ME_true->GetBinCenter(uBin+1);
  //  double VAL = 5.40165e-08-1.06309e-08*pow(MOM,1.)+5.71623e-09*pow(MOM,2.)-
  //  -1.15454e-11*pow(MOM,3.)+1.17292e-14*pow(MOM,4.)-7.94266e-18*pow(MOM,5.)+
  //  3.70109e-21*pow(MOM,6.)-1.10644e-24*pow(MOM,7.)+1.86899e-28*pow(MOM,8.)-1.34242e-32*pow(MOM,9.);
    double VAL = fitPS->Eval(MOM);
    h_pXi_ME_true->SetBinContent(uBin+1,VAL);
    h_pXi_ME_true->SetBinError(uBin+1,0);
    h_pXi_SE_true->SetBinContent(uBin+1,VAL*Ck_pXi.Eval(MOM));
    h_pXi_SE_true->SetBinError(uBin+1,0);
    h_pXi_Ck_true->SetBinContent(uBin+1,Ck_pXi.Eval(MOM));
    h_pXi_Ck_true->SetBinError(uBin+1,0);
  }
  h_pXi_SE_true->Scale(1./h_pXi_SE_true->Integral());
  h_pXi_ME_true->Scale(1./h_pXi_ME_true->Integral());

  printf("Folding fun\n");
  //fold the correlation,SE and ME
  DLM_Unfold dlmFold;
  dlmFold.SetResponse(hResolution_pXi);
  dlmFold.SetSilentMode(false);

  dlmFold.SetData(h_pXi_SE_true);
  DLM_Histo<float>* dlm_pXi_SE_exp = dlmFold.Fold();
  printf("Fold dlm_pXi_SE_exp\n");
  TH1F* h_pXi_SE_exp = Convert_DlmHisto_TH1F(dlm_pXi_SE_exp,"h_pXi_SE_exp");

  dlmFold.SetData(h_pXi_ME_true);
  DLM_Histo<float>* dlm_pXi_ME_exp = dlmFold.Fold();
  printf("Fold dlm_pXi_ME_exp\n");
  TH1F* h_pXi_ME_exp = Convert_DlmHisto_TH1F(dlm_pXi_ME_exp,"h_pXi_ME_exp");

  DLM_Histo<float>* dlm_pXi_Ck_exp(dlm_pXi_SE_exp);
  dlm_pXi_Ck_exp[0] /= dlm_pXi_ME_exp[0];
  TH1F* h_pXi_Ck_exp = Convert_DlmHisto_TH1F(dlm_pXi_Ck_exp,"h_pXi_Ck_exp");

  dlmFold.SetData(h_pXi_Ck_true);
  DLM_Histo<float>* dlm_pXi_Ck_truefold = dlmFold.Fold();
  printf("Fold dlm_pXi_Ck_truefold\n");
  TH1F* h_pXi_Ck_truefold = Convert_DlmHisto_TH1F(dlm_pXi_Ck_truefold,"h_pXi_Ck_truefold");


  //from CATS
  DLM_CkDecomposition CkDec_pXi_1530("pXi1530",0,Ck_pXi1530,NULL);

  DLM_CkDecomposition CkDec_pXi_direct("pXi",3,Ck_pXi,hResolution_pXi);
  CkDec_pXi_direct.AddContribution(0,lambda_xi1530,DLM_CkDecomposition::cFeedDown,&CkDec_pXi_1530,hResidual_pXi_pXi1530);
  //CkDec_pXi_direct.AddPhaseSpace(0, h_pXi_ME_true);
  CkDec_pXi_direct.AddContribution(1,lambda_flat,DLM_CkDecomposition::cFeedDown);
  CkDec_pXi_direct.AddContribution(2,lambda_misid,DLM_CkDecomposition::cFake);
  CkDec_pXi_direct.Update();
  DLM_CkDecomposition CkDec_pXi_PS("pXi",3,Ck_pXi,hResolution_pXi);
  CkDec_pXi_PS.AddPhaseSpace(h_pXi_ME_true);
  CkDec_pXi_PS.AddContribution(0,lambda_xi1530,DLM_CkDecomposition::cFeedDown,&CkDec_pXi_1530,hResidual_pXi_pXi1530);
  CkDec_pXi_PS.AddPhaseSpace(0, h_pXi_ME_true);
  CkDec_pXi_PS.AddContribution(1,lambda_flat,DLM_CkDecomposition::cFeedDown);
  CkDec_pXi_PS.AddContribution(2,lambda_misid,DLM_CkDecomposition::cFake);
  CkDec_pXi_PS.Update();

  DLM_CkDecomposition CkDec_pXi_PSexp("pXi",3,Ck_pXi,hResolution_pXi);
  CkDec_pXi_PSexp.AddPhaseSpace(h_pXi_ME_exp);
  CkDec_pXi_PSexp.AddContribution(0,lambda_xi1530,DLM_CkDecomposition::cFeedDown,&CkDec_pXi_1530,hResidual_pXi_pXi1530);
  CkDec_pXi_PSexp.AddPhaseSpace(0, h_pXi_ME_exp);
  CkDec_pXi_PSexp.AddContribution(1,lambda_flat,DLM_CkDecomposition::cFeedDown);
  CkDec_pXi_PSexp.AddContribution(2,lambda_misid,DLM_CkDecomposition::cFake);
  CkDec_pXi_PSexp.Update();

  TH1F* hCk_CATS_direct = new TH1F("hCk_CATS_direct","hCk_CATS_direct",NumMomBins_pXi,kMin_pXi,kMax_pXi);
  TH1F* hCk_CATS_PS = new TH1F("hCk_CATS_PS","hCk_CATS_PS",NumMomBins_pXi,kMin_pXi,kMax_pXi);
  TH1F* hCk_CATS_PSexp = new TH1F("hCk_CATS_PSexp","hCk_CATS_PSexp",NumMomBins_pXi,kMin_pXi,kMax_pXi);

  for(unsigned uBin=0; uBin<NumMomBins_pXi; uBin++){
    double MOM = AB_pXi.GetMomentum(uBin);
    hCk_CATS_direct->SetBinContent(uBin+1,CkDec_pXi_direct.EvalCk(MOM));
    hCk_CATS_PS->SetBinContent(uBin+1,CkDec_pXi_PS.EvalCk(MOM));
    hCk_CATS_PSexp->SetBinContent(uBin+1,CkDec_pXi_PSexp.EvalCk(MOM));

    hCk_CATS_direct->SetBinError(uBin+1,0);
    hCk_CATS_PS->SetBinError(uBin+1,0);
    hCk_CATS_PSexp->SetBinError(uBin+1,0);
    //printf("Ck_pXi(%.0f) = %.2f\n",MOM,Ck_pXi.GetBinContent(uBin));
    //printf("CkDec_pXi_direct(%.0f) = %.2f\n",MOM,CkDec_pXi_direct.EvalCk(MOM));
    //printf("CkDec_pXi_PS(%.0f) = %.2f\n",MOM,CkDec_pXi_PS.EvalCk(MOM));
    //printf("CkDec_pXi_PSexp(%.0f) = %.2f\n",MOM,CkDec_pXi_PSexp.EvalCk(MOM));
    //printf("\n");
  }

  TH1F* hData_pXi = AnalysisObject.GetAliceExpCorrFun("pp13TeV_HM_Dec19","pLambda","_0",2,false,-1);
  //unfold the correlation,SE and ME, building the SE/ME to compare to the direct correlation unfold

  for(unsigned uBin=0; uBin<NumMomBins_pXi; uBin++){
    double CkVal;
    double CkErr;

    CkVal=h_pXi_Ck_true->GetBinContent(uBin+1);
    CkErr=h_pXi_Ck_true->GetBinError(uBin+1);
    h_pXi_Ck_true->SetBinContent(uBin+1,Norm*(lambda_gen*CkVal+1-lambda_gen));
    h_pXi_Ck_true->SetBinError(uBin+1,Norm*lambda_gen*CkErr);

    CkVal=h_pXi_Ck_truefold->GetBinContent(uBin+1);
    CkErr=h_pXi_Ck_truefold->GetBinError(uBin+1);
    h_pXi_Ck_truefold->SetBinContent(uBin+1,Norm*(lambda_gen*CkVal+1-lambda_gen));
    h_pXi_Ck_truefold->SetBinError(uBin+1,Norm*lambda_gen*CkErr);

    CkVal=h_pXi_Ck_exp->GetBinContent(uBin+1);
    CkErr=h_pXi_Ck_exp->GetBinError(uBin+1);
    h_pXi_Ck_exp->SetBinContent(uBin+1,Norm*(lambda_gen*CkVal+1-lambda_gen));
    h_pXi_Ck_exp->SetBinError(uBin+1,Norm*lambda_gen*CkErr);

    CkVal=hCk_CATS_direct->GetBinContent(uBin+1);
    CkErr=hCk_CATS_direct->GetBinError(uBin+1);
    hCk_CATS_direct->SetBinContent(uBin+1,Norm*(lambda_gen*CkVal+1-lambda_gen));
    hCk_CATS_direct->SetBinError(uBin+1,Norm*lambda_gen*CkErr);

    CkVal=hCk_CATS_PS->GetBinContent(uBin+1);
    CkErr=hCk_CATS_PS->GetBinError(uBin+1);
    hCk_CATS_PS->SetBinContent(uBin+1,Norm*(lambda_gen*CkVal+1-lambda_gen));
    hCk_CATS_PS->SetBinError(uBin+1,Norm*lambda_gen*CkErr);

    CkVal=hCk_CATS_PSexp->GetBinContent(uBin+1);
    CkErr=hCk_CATS_PSexp->GetBinError(uBin+1);
    hCk_CATS_PSexp->SetBinContent(uBin+1,Norm*(lambda_gen*CkVal+1-lambda_gen));
    hCk_CATS_PSexp->SetBinError(uBin+1,Norm*lambda_gen*CkErr);
  }

  TH1F* h_pXi_Ck_exp_SmallRange = (TH1F*)hCk_CATS_direct->Clone("hCkRatio_CATS_direct");
  for(unsigned uBin=0; uBin<NumMomBins_pXi; uBin++){
    h_pXi_Ck_exp_SmallRange->SetBinContent(uBin+1,h_pXi_Ck_exp->GetBinContent(uBin+1));
    h_pXi_Ck_exp_SmallRange->SetBinError(uBin+1,0);
  }

  TH1F* hCkRatio_CATS_direct = (TH1F*)hCk_CATS_direct->Clone("hCkRatio_CATS_direct");
  hCkRatio_CATS_direct->Divide(h_pXi_Ck_exp_SmallRange);
  TH1F* hCkRatio_CATS_PS = (TH1F*)hCk_CATS_PS->Clone("hCkRatio_CATS_PS");
  hCkRatio_CATS_PS->Divide(h_pXi_Ck_exp_SmallRange);
  TH1F* hCkRatio_CATS_PSexp = (TH1F*)hCk_CATS_PSexp->Clone("hCkRatio_CATS_PSexp");
  hCkRatio_CATS_PSexp->Divide(h_pXi_Ck_exp_SmallRange);

  TH1F* hUpperErrors1 = new TH1F("hUpperErrors1","hUpperErrors1",NumMomBins_pXi,kMin_pXi,kMax_pXi);
  TH1F* hLowerErrors1 = new TH1F("hLowerErrors1","hLowerErrors1",NumMomBins_pXi,kMin_pXi,kMax_pXi);
  TH1F* hUpperErrors3 = new TH1F("hUpperErrors3","hUpperErrors3",NumMomBins_pXi,kMin_pXi,kMax_pXi);
  TH1F* hLowerErrors3 = new TH1F("hLowerErrors3","hLowerErrors3",NumMomBins_pXi,kMin_pXi,kMax_pXi);
  TGraph gUpperErrors1;
  gUpperErrors1.SetName("gUpperErrors1");
  TGraph gLowerErrors1;
  gLowerErrors1.SetName("gLowerErrors1");
  TGraph gUpperErrors3;
  gUpperErrors3.SetName("gUpperErrors3");
  TGraph gLowerErrors3;
  gLowerErrors3.SetName("gLowerErrors3");


  hUpperErrors1->SetLineColor(kGreen+1);
  hUpperErrors1->SetLineWidth(3);
  hUpperErrors1->SetLineStyle(1);
  hLowerErrors1->SetLineColor(kGreen+1);
  hLowerErrors1->SetLineWidth(3);
  hLowerErrors1->SetLineStyle(1);

  hUpperErrors3->SetLineColor(kRed);
  hUpperErrors3->SetLineWidth(3);
  hUpperErrors3->SetLineStyle(1);
  hLowerErrors3->SetLineColor(kRed);
  hLowerErrors3->SetLineWidth(3);
  hLowerErrors3->SetLineStyle(1);

  gUpperErrors1.SetLineColor(kGreen+1);
  gUpperErrors1.SetLineWidth(3);
  gUpperErrors1.SetLineStyle(1);
  gLowerErrors1.SetLineColor(kGreen+1);
  gLowerErrors1.SetLineWidth(3);
  gLowerErrors1.SetLineStyle(1);

  gUpperErrors3.SetLineColor(kRed);
  gUpperErrors3.SetLineWidth(3);
  gUpperErrors3.SetLineStyle(1);
  gLowerErrors3.SetLineColor(kRed);
  gLowerErrors3.SetLineWidth(3);
  gLowerErrors3.SetLineStyle(1);

  h_pXi_Ck_exp->SetLineColor(kBlack);
  h_pXi_Ck_exp->SetLineWidth(7);
  h_pXi_Ck_exp->SetLineStyle(1);

  hCk_CATS_direct->SetLineColor(kRed-3);
  hCk_CATS_direct->SetLineWidth(6);
  hCk_CATS_direct->SetLineStyle(2);
  hCkRatio_CATS_direct->SetLineColor(kRed-3);
  hCkRatio_CATS_direct->SetLineWidth(6);
  hCkRatio_CATS_direct->SetLineStyle(2);

  hCk_CATS_PSexp->SetLineColor(kViolet);
  hCk_CATS_PSexp->SetLineWidth(6);
  hCk_CATS_PSexp->SetLineStyle(3);
  hCkRatio_CATS_PSexp->SetLineColor(kViolet);
  hCkRatio_CATS_PSexp->SetLineWidth(6);
  hCkRatio_CATS_PSexp->SetLineStyle(5);

  hCk_CATS_PS->SetLineColor(kAzure);
  hCk_CATS_PS->SetLineWidth(5);
  hCk_CATS_PS->SetLineStyle(1);
  hCkRatio_CATS_PS->SetLineColor(kAzure);
  hCkRatio_CATS_PS->SetLineWidth(5);
  hCkRatio_CATS_PS->SetLineStyle(1);

  for(unsigned uBin=0; uBin<NumMomBins_pXi; uBin++){
    double RelErr = hData_pXi->GetBinError(uBin+1)/hData_pXi->GetBinContent(uBin+1);
    double MOM = hData_pXi->GetBinCenter(uBin+1);
    hUpperErrors1->SetBinContent(uBin+1,1.*(1.+RelErr));
    hLowerErrors1->SetBinContent(uBin+1,1./(1.+RelErr));
    hUpperErrors3->SetBinContent(uBin+1,1.*(1.+3.*RelErr));
    hLowerErrors3->SetBinContent(uBin+1,1./(1.+3.*RelErr));

    hUpperErrors1->SetBinError(uBin+1,0);
    hLowerErrors1->SetBinError(uBin+1,0);
    hUpperErrors3->SetBinError(uBin+1,0);
    hLowerErrors3->SetBinError(uBin+1,0);

    gUpperErrors1.SetPoint(uBin,MOM,1.*(1.+RelErr));
    gLowerErrors1.SetPoint(uBin,MOM,1./(1.+RelErr));
    gUpperErrors3.SetPoint(uBin,MOM,1.*(1.+3.*RelErr));
    gLowerErrors3.SetPoint(uBin,MOM,1./(1.+3.*RelErr));
  }

  //for the plots, plug in error h_pXi_Ck_exp
  for(unsigned uBin=0; uBin<NumMomBins_pXi; uBin++){
    h_pXi_Ck_exp->SetBinError(uBin+1,hData_pXi->GetBinError(uBin+1));
  }

  OutputFile->cd();
  hData_pXi->Write();
  h_pXi_SE_true->Write();
  h_pXi_ME_true->Write();
  h_pXi_Ck_true->Write();
  h_pXi_Ck_truefold->Write();
  //h_pXi_Ck_perfect->Write();

  h_pXi_SE_exp->Write();
  h_pXi_ME_exp->Write();
  h_pXi_Ck_exp->Write();//the true one

  hCk_CATS_direct->Write();
  hCk_CATS_PS->Write();
  hCk_CATS_PSexp->Write();

  hCkRatio_CATS_direct->Write();
  hCkRatio_CATS_PS->Write();
  hCkRatio_CATS_PSexp->Write();

  hUpperErrors1->Write();
  hLowerErrors1->Write();
  hUpperErrors3->Write();
  hLowerErrors3->Write();

  gUpperErrors1.Write();
  gLowerErrors1.Write();
  gUpperErrors3.Write();
  gLowerErrors3.Write();
  //h_pXi_Ck_fold->Write();


  delete h_pXi_SE_true;
  delete h_pXi_ME_true;
  delete h_pXi_Ck_true;
  //delete h_pXi_Ck_perfect;
  delete h_pXi_SE_exp;
  delete h_pXi_ME_exp;
  delete hCk_CATS_direct;
  delete hCk_CATS_PS;
  delete hCk_CATS_PSexp;
  //delete h_pXi_Ck_fold;
  delete fitPS;
  delete hData_pXi;
  delete OutputFile;
}


//const double& Momentum, const double* SourcePar, const double* PotPar
double LedniFit(double* x, double* par){
  return par[0]*Lednicky_Identical_Singlet(*x,&par[1],&par[2]);
}

//attempt to find scattering parameters or a potential that
//can describe (approx) a Gaussian shaped correlation function
void Rafa_2body_expCk(){
  TString WorkFolder = TString::Format("%s/OtherTasks/Rafa_2body_expCk/",GetFemtoOutputFolder());
  TString InputFile = WorkFolder+"femto.root";
  TString InputHisto = "CF2PartNorm";
  const double kMin = 0;
  const double kMax = 250;
  TFile fInput(InputFile,"read");
  TH1D* hInput  = (TH1D*)fInput.Get(InputHisto);
  TF1* fLedni = new TF1("fLedni",LedniFit,kMin,kMax,4);
  fLedni->SetParameter(0,1);
  fLedni->SetParameter(1,1.2);
  fLedni->SetParameter(2,2.0);
  fLedni->SetParameter(3,0.1);
  hInput->Fit(fLedni,"S, N, R, M");
  printf("at 50: %f\n",fLedni->Eval(50));
  TH1F* hFit = new TH1F("hFit","hFit",int(kMax-kMin),kMin,kMax);
  for(unsigned uBin=0; uBin<hFit->GetNbinsX(); uBin++){
    hFit->SetBinContent(uBin+1,fLedni->Eval(hFit->GetBinCenter(uBin+1)));
  }

  TFile fOutput(WorkFolder+"FuckThisShit.root","recreate");
  hInput->Write();
  hFit->Write();
  fLedni->Write();

  delete hInput;
}

//0 = norm
//1 = radius
//2,3,4,5 -> double gaussian potential
CATS* KITTY_FITTER=NULL;
double CATS_FITTER(double* x, double* par){
  KITTY_FITTER->SetAnaSource(0,par[1],true);
  KITTY_FITTER->SetShortRangePotential(0,0,0,par[2]);
  KITTY_FITTER->SetShortRangePotential(0,0,1,par[3]);
  KITTY_FITTER->SetShortRangePotential(0,0,2,par[4]);
  KITTY_FITTER->SetShortRangePotential(0,0,3,par[5]);
  KITTY_FITTER->KillTheCat();
  return par[0]*KITTY_FITTER->EvalCorrFun(*x);
}

//same as above, but with a Gaussian potential
void Rafa_2body_expCk_CATS(){
  TString WorkFolder = TString::Format("%s/OtherTasks/Rafa_2body_expCk/",GetFemtoOutputFolder());
  TString InputFile = WorkFolder+"femto.root";
  TString InputHisto = "CF2PartNorm";
  const double kMin = 0;
  const double kMax = 240;
  CATS Kitty;
  Kitty.SetMomBins(TMath::Nint(kMax/2.),kMin,kMax);
  CATSparameters cPars (CATSparameters::tSource,1,true);
  cPars.SetParameter(0,1.2);
  Kitty.SetAnaSource(GaussSource, cPars);
  Kitty.SetUseAnalyticSource(true);
  Kitty.SetMomentumDependentSource(false);
  Kitty.SetExcludeFailedBins(false);
  Kitty.SetQ1Q2(0);
  Kitty.SetQuantumStatistics(false);
  //Kitty.SetRedMass( (938.*1116.)/(938.+1116.) );
  Kitty.SetRedMass( 0.5*938. );

  Kitty.SetNumChannels(1);
  Kitty.SetNumPW(0,1);
  Kitty.SetSpin(0,0);
  Kitty.SetChannelWeight(0, 1.);
  CATSparameters pPars(CATSparameters::tPotential,4,true);
  pPars.SetParameter(0,-50);
  pPars.SetParameter(1,1.0);
  pPars.SetParameter(2,0);
  pPars.SetParameter(3,1);
  Kitty.SetShortRangePotential(0,0,DoubleGaussSum,pPars);
  Kitty.SetNotifications(CATS::nWarning);
  Kitty.SetMaxRad(128);
  Kitty.SetMaxRho(64);
  Kitty.KillTheCat();
  KITTY_FITTER = &Kitty;

  TFile fInput(InputFile,"read");
  TH1D* hInput  = (TH1D*)fInput.Get(InputHisto);
  TF1* fCATS = new TF1("fCATS",CATS_FITTER,kMin,kMax,6);

  fCATS->SetParameter(0,1);
  fCATS->SetParameter(1,1.2);//r
  fCATS->SetParLimits(1,1.0,1.5);
  fCATS->SetParameter(2,-20);//V0
  fCATS->SetParLimits(2,-100,0);
  fCATS->SetParameter(3,1.0);//mu0
  fCATS->SetParLimits(3,0.5,5.0);
  fCATS->SetParameter(4,100.0);//V1
  fCATS->SetParLimits(4,0,1000);
  fCATS->SetParameter(5,0.1);//mu1
  fCATS->SetParLimits(5,0.05,0.5);

  /*
  fCATS->FixParameter(0,1.01191e+00);
  fCATS->FixParameter(1,1.20000e+00);
  fCATS->FixParameter(2,-1.01057e+01);
  fCATS->FixParameter(3,2.54662e+00);
  fCATS->FixParameter(4,8.11437e+02);
  fCATS->FixParameter(5,6.76762e-02);
  */
  hInput->Fit(fCATS,"S, N, R, M");
  //printf("at 50: %f\n",fLedni->Eval(50));
  TH1F* hFitResult = new TH1F("hFitResult","hFitResult",int(kMax-kMin),kMin,kMax);
  for(unsigned uBin=0; uBin<hFitResult->GetNbinsX(); uBin++){
    hFitResult->SetBinContent(uBin+1,fCATS->Eval(hFitResult->GetBinCenter(uBin+1)));
  }

  TH1F* h_kcotd=NULL;
  TF1* f_kcotd=NULL;
  double c_f0,c_d0=1;
  c_f0=0;
  c_d0=1;
  //(CATS& Kitty, double& ScatLen, double& EffRan, TH1F*& hFit, TF1*& fitSP, const int& Nterms=2, const bool& Fixf0=false, const bool& Fixd0=false)
  Kitty.SetEpsilonConv(1e-8);
  Kitty.SetEpsilonProp(1e-8);
  //V0=-10.11; mu=2.55; V1=811.44; mu1=0.07
  //1  p0           1.01191e+00   3.07375e-04   0.00000e+00   3.73857e+03
  //2  p1           1.20000e+00     fixed
  //3  p2          -1.01057e+01   5.58941e-05  -0.00000e+00  -1.21703e+05
  //4  p3           2.54662e+00   1.08960e-06  -0.00000e+00   9.07165e+04
  //5  p4           8.11437e+02   2.70487e-01  -0.00000e+00   2.20305e+03
  //6  p5           6.76762e-02   3.80730e-06   0.00000e+00  -2.79412e+04

  Eval_ScattParameters(Kitty,c_f0,c_d0,h_kcotd,f_kcotd,2,false,false);
  printf("f0 = %.2f; d0 = %.2f\n",c_f0,c_d0);
  printf("V0=%.2f; mu=%.2f; V1=%.2f; mu1=%.2f\n",fCATS->GetParameter(2),fCATS->GetParameter(3),fCATS->GetParameter(4),fCATS->GetParameter(5));

  TFile fOutput(WorkFolder+"FuckThisCat.root","recreate");
  hInput->Write();
  hFitResult->Write();
  fCATS->Write();
  h_kcotd->Write();
  f_kcotd->Write();

  delete fCATS;
  delete h_kcotd;
  delete f_kcotd;
  delete hInput;
}



void rootmathboost_test1(){
/*
  TLorentzVector v1;
  v1.SetXYZM(50,-50,100,1000);
  TLorentzVector v2;
  v2.SetXYZM(100,-20,200,1000);
  ROOT::Math::PtEtaPhiMVector<double> part1(v1.Pt(),v1.Eta(),v1.Phi(),v1.M());
  ROOT::Math::PtEtaPhiMVector<double> part2(v2.Pt(),v2.Eta(),v2.Phi(),v2.M());


  TLorentzVector sum = v1+v2;
  TVector3 bv = -sum.BoostVector();
  v1.Boost(bv);
  v2.Boost(bv);
  printf("x: %.1f, %.1f\n",v1.X(),v2.X());
  printf("y: %.1f, %.1f\n",v1.Y(),v2.Y());
  printf("z: %.1f, %.1f\n",v1.Z(),v2.Z());
*/
}


//we simulate an yield following a Levy distribution (to kind of mimic ME)
//and than random sample from it. We fold it, and than we try to unfold it.
void Test_unfold_yield(){
  const long unsigned Yield = (long unsigned)(10)*(long unsigned)(1000*1000*1000);
  const TString OutputFile = TString::Format("%s/Test_unfold_yield/fOutput.root",GetFemtoOutputFolder());
  const unsigned NumBins = 1024/4;
  const double kMin = 0;
  const double kMax = 1024;
  const double kFold = 360;
  const double kUnfold = 320;
  TRandom3 rangen(11);
  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());
  TH2F* hResolution_pp = AnalysisObject.GetResolutionMatrix("pp13TeV_HM_DimiJun20","pp");
  TFile* fOutput = new TFile(OutputFile,"recreate");
  TF1* fYield = new TF1("fYield","[0]*x*x*exp(-pow(x/[1],[2]))/(sqrt(938.*938.+x*x)*sqrt(1116.*1116.+x*x))",0,8192);
  fYield->SetParameter(0,1);
  fYield->SetParameter(1,1000);
  fYield->SetParameter(2,1.45);
  double fIntegral = fYield->Integral(0,8192);
  fYield->SetParameter(0,double(Yield)/fIntegral);

  TH1F* hExpYield = new TH1F("hExpYield","hExpYield",NumBins,kMin,kMax);
  TH1F* hYield = new TH1F("hYield","hYield",NumBins,kMin,kMax);
  //const unsigned short& dim,const double& stability=2, const double& location=0, const double& scale=1, const double& skewness=0
  for(unsigned uBin=0; uBin<hYield->GetNbinsX(); uBin++){
    double kStar = hYield->GetBinCenter(uBin+1);
    double klow = hYield->GetXaxis()->GetBinLowEdge(uBin+1);
    double kup = hYield->GetXaxis()->GetBinUpEdge(uBin+1);
    double ExpectedYield = fYield->Integral(klow,kup);
    //double RandomYield = rangen.Poisson(EffectiveCounts)*OriginalValue/EffectiveCounts;
    double RandomYield = rangen.Poisson(ExpectedYield);
    hExpYield->SetBinContent(uBin+1,ExpectedYield);
    hYield->SetBinContent(uBin+1,RandomYield);
  }

  DLM_Unfold dlmFold;
  dlmFold.SetData(hYield);
  dlmFold.SetResponse(hResolution_pp);
  dlmFold.SetSilentMode(false);
  dlmFold.SetFoldRange(0,kFold);
  printf("Folding...\n");
  DLM_Histo<float>* dlm_yield_fold = dlmFold.Fold();
  printf("Folded\n");
  fOutput->cd();
  TH1F* h_yield_fold = Convert_DlmHisto_TH1F(dlm_yield_fold,"h_yield_fold");
  for(unsigned uBin=0; uBin<h_yield_fold->GetNbinsX(); uBin++){
    double BinVal = h_yield_fold->GetBinContent(uBin+1);
    double BinErr = h_yield_fold->GetBinError(uBin+1);
    BinVal = TMath::Nint(BinVal);
    BinErr = sqrt(BinVal);
    if(!BinErr) BinErr = 1;
    h_yield_fold->SetBinContent(uBin+1,BinVal);
    h_yield_fold->SetBinError(uBin+1,BinErr);
  }
  //h_yield_fold->Sumw2();


  DLM_Unfold dlmUnfold;
  dlmUnfold.SetData(h_yield_fold);
  dlmUnfold.SetResponse(hResolution_pp);
  dlmUnfold.SetSilentMode(false);
  dlmUnfold.SetUnfoldPrecision(0.1,1.0);
  dlmUnfold.SetUnfoldRange(0,kUnfold);
  dlmUnfold.SetUnfoldMinutes(115);
  printf("Unfolding...\n");
  DLM_Histo<float>*  dlm_yield_unfolded = dlmUnfold.Unfold(1);
  printf("Unfolded\n");
  fOutput->cd();
  TH1F* h_yield_unfolded = Convert_DlmHisto_TH1F(dlm_yield_unfolded,"h_yield_unfolded");
  fOutput->cd();
  hResolution_pp->Write();
  hExpYield->Write();
  hYield->Write();
  fYield->Write();
  h_yield_fold->Write();
  h_yield_unfolded->Write();

  delete hExpYield;
  delete hYield;
  delete fYield;
  delete h_yield_fold;
  delete h_yield_unfolded;
  delete fOutput;
}



void Bhawani_pot(TString PotentialName) {

	double Mass_p = 1116;
	double Mass_Kch = 1116;
	int NumMomBins = 60;
	//const TString PotentialName = "NSC97b";
	//const TString PotentialName = "NF48";
	//const TString PotentialName = "emma";
	//const TString PotentialName = "custom";

	const double kMin = 0;
	const double kMax = 300;
	const double kStep = 3;
	const unsigned nMom = TMath::Nint(kMax / kStep);
	const double Radius = 1.0;

	CATSparameters sPars(CATSparameters::tSource, 1, true);
	sPars.SetParameter(0, Radius);
	CATS Kitty_SE;
	Kitty_SE.SetMomBins(nMom, kMin, kMax);
	Kitty_SE.SetAnaSource(GaussSource, sPars);
	Kitty_SE.SetUseAnalyticSource(true);
	Kitty_SE.SetQ1Q2(0);
	Kitty_SE.SetQuantumStatistics(false);
	//Kitty_SE.SetRedMass(Mass_L*0.5);
	Kitty_SE.SetRedMass((Mass_p * Mass_Kch) / (Mass_p + Mass_Kch));
	Kitty_SE.SetNumChannels(1);
	Kitty_SE.SetNumPW(0, 1);
	Kitty_SE.SetSpin(0, 0);
	Kitty_SE.SetChannelWeight(0, 1.);
	CATSparameters pPars(CATSparameters::tPotential, 7, true);
	double c_f0, c_d0;
	if (PotentialName == "NSC97b") {
		pPars.SetParameter(0, -78.42);
		pPars.SetParameter(1, 1.0);
		pPars.SetParameter(2, 741.76);
		pPars.SetParameter(3, 0.45);
		c_f0 = 0.397;
		c_d0 = 10.360;
	} else if (PotentialName == "NF48") {
		pPars.SetParameter(0, -1647.40);
		pPars.SetParameter(1, 0.6);
		pPars.SetParameter(2, 3888.96);
		pPars.SetParameter(3, 0.45);
		c_f0 = 1.511;
		c_d0 = 2.549;
	} else if (PotentialName == "emma") {
		//for emma -> NSC97f basis
		pPars.SetParameter(0, -106.53 * 0.85);
		pPars.SetParameter(1, 1.0 * 1.18);
		pPars.SetParameter(2, 1469.33);
		pPars.SetParameter(3, 0.45 * 1.1);
		c_f0 = 0.350;
		c_d0 = 16.330;
	} else if (PotentialName == "Toy1") {
		pPars.SetParameter(0, -144.5);
		pPars.SetParameter(1, 2.11);
		pPars.SetParameter(2, 520.0);
		pPars.SetParameter(3, 0.54);
		c_f0 = -0.73;
		c_d0 = 7.72;
	} else if (PotentialName == "ND46") {//bound
		pPars.SetParameter(0, -144.89);
		pPars.SetParameter(1, 1.0);
		pPars.SetParameter(2, 127.87);
		pPars.SetParameter(3, 0.45);
		c_f0 = -4.621;
		c_d0 = 1.3;
	} else if (PotentialName == "Yukawa1") {
		pPars.SetParameter(0, 1.0);
		pPars.SetParameter(1, 1.0);
		pPars.SetParameter(2, 100.0);
		pPars.SetParameter(3, 1.0);
		pPars.SetParameter(4, 0.4);
		c_f0 = 0;
		c_d0 = 1;
	} else if (PotentialName == "pKplusI0") {
		pPars.SetParameter(0, 0.0);
		c_f0 = 0.03;
		c_d0 = 0.0;
	} else if (PotentialName == "pKplusI1") {
		pPars.SetParameter(0, 1.0);
		c_f0 = -0.3;
		c_d0 = 0.0;
	} else if (PotentialName == "pKplusYuki") {
		pPars.SetParameter(0, 0.376); //0.376;0.335
		pPars.SetParameter(1, sqrt(200.*(Mass_p * Mass_Kch) / (Mass_p + Mass_Kch)));
		pPars.SetParameter(2, 3);
		pPars.SetParameter(3, 2084);
		pPars.SetParameter(4, 50.81);
		pPars.SetParameter(5, 18.34);
		pPars.SetParameter(6, -1.752);
		c_f0 = -0.3;
		c_d0 = 0.0;
	} else {
		pPars.SetParameter(0, -5.50337);
		pPars.SetParameter(1, 1. / sqrt(2.148805));
		pPars.SetParameter(2, 0);
		pPars.SetParameter(3, 1);
		//pPars.SetParameter(0,-78.42*0.39*4.5);//0.39,0.4
		//pPars.SetParameter(1,1.0*1.35);
		//pPars.SetParameter(2,741.76*4.5);
		//pPars.SetParameter(3,0.45*1.4);
		//NF46 as a stariting point
		//pPars.SetParameter(0,-1327.26*1.0);
		//pPars.SetParameter(1,0.6);
		//pPars.SetParameter(2,2561.56);
		//pPars.SetParameter(3,0.45);
		c_f0 = 0.02;
		c_d0 = 30.0;
		printf(" Hello\n");

	}
	Kitty_SE.SetEpsilonConv(5e-9);
	Kitty_SE.SetEpsilonProp(5e-9);
	//if (PotentialName.Contains("Yukawa"))
	//	Kitty_SE.SetShortRangePotential(0, 0, YukawaDimiCore, pPars);
	if (PotentialName.Contains("pKplusI")) {
		Kitty_SE.SetShortRangePotential(0, 0, KpProtonEquivalentPotential, pPars);
		Kitty_SE.SetEpsilonConv(1e-9);
		Kitty_SE.SetEpsilonProp(1e-9);
	//} else if (PotentialName.Contains("pKplusYuki")) {
	///	Kitty_SE.SetShortRangePotential(0, 0, SingleGaussDynamic, pPars);
		//Kitty_SE.SetEpsilonConv(1e-9);
		//Kitty_SE.SetEpsilonProp(1e-9);
	} else Kitty_SE.SetShortRangePotential(0, 0, DoubleGaussSum, pPars);

	Kitty_SE.KillTheCat();

	TFile* OutputFile = new TFile(
	    TString::Format("%s/OtherTasks/Bhawani_%s.root",GetFemtoOutputFolder(),PotentialName.Data()), "recreate");
	printf("File Created\n");
	TGraph gKitty;
	gKitty.SetName(TString::Format("gKitty"));
	gKitty.Set(NumMomBins);
	for (unsigned uBin = 0; uBin < NumMomBins; uBin++) {
		printf("C(%.2f) = %.2f\n", Kitty_SE.GetMomentum(uBin), Kitty_SE.GetCorrFun(uBin));
		gKitty.SetPoint(uBin, Kitty_SE.GetMomentum(uBin), Kitty_SE.GetCorrFun(uBin));
	}
	gKitty.Write();
	delete OutputFile;
	//delete cPars;
}


//D0-Dstar_ch
void DDstar_TetsuoTest1(){

  const double Mass_D0 = 1864.84;
  const double Mass_Dch = 1869.5;
  const double Mass_Dstar0 = 2006.85;
  const double Mass_Dstarch = 2010.26;
  const double Mass_pi0 = 134.9768;

  const int NumSource = 6;
  double* SourceSize = new double[NumSource];
  SourceSize[0] = 1.0;
  SourceSize[1] = 1.5;
  SourceSize[2] = 2.0;
  SourceSize[3] = 3.0;
  SourceSize[4] = 4.0;
  SourceSize[5] = 5.0;

  const int NumPots = 2;
  double* par_V0 = new double[NumPots];
  double* par_R = new double[NumPots];

  par_V0[0] = -34.5;
  par_R[0] = NuToFm*(1./Mass_pic);

  par_V0[1] = -33.1;
  par_R[1] = NuToFm*(1./Mass_pic);

  TGraph** gCk = new TGraph* [NumPots];
  for(unsigned uPot=0; uPot<NumPots; uPot++){
    gCk[uPot] = new TGraph[NumSource];
    for(unsigned uSor=0; uSor<NumSource; uSor++){
      gCk[uPot][uSor].SetName(TString::Format("Ck_DDstar_%.1fMeV_%.1ffm",par_V0[uPot],SourceSize[uSor]));
    }
  }
  TGraph* gYuki = new TGraph [NumSource];
  for(unsigned uSor=0; uSor<NumSource; uSor++){
    gYuki[uSor].SetName(TString::Format("Yuki_DDstar_%.1ffm",SourceSize[uSor]));
  }

  for(unsigned uPot=0; uPot<NumPots; uPot++){
    CATS Kitty;
    Kitty.SetMomBins(300,0,300);
    CATSparameters cPars(CATSparameters::tSource,1,true);
    cPars.SetParameter(0,1);
    Kitty.SetAnaSource(GaussSource, cPars);
    Kitty.SetUseAnalyticSource(true);
    Kitty.SetAutoNormSource(false);
    Kitty.SetQ1Q2(0);
    Kitty.SetQuantumStatistics(false);
    Kitty.SetRedMass( (Mass_D0*Mass_Dstarch)/(Mass_D0+Mass_Dstarch) );
    Kitty.SetNumChannels(1);
    Kitty.SetNumPW(0,1);
    Kitty.SetSpin(0,0);
    Kitty.SetChannelWeight(0, 1.);
    Kitty.SetEpsilonConv(5e-9);
    Kitty.SetEpsilonProp(5e-9);
    Kitty.SetMaxRad(96);
    Kitty.SetMaxRho(32);
    Kitty.SetNotifications(CATS::nWarning);
    CATSparameters pPars(CATSparameters::tPotential,2,true);
    pPars.SetParameter(0,par_V0[uPot]);
    pPars.SetParameter(1,par_R[uPot]);
    Kitty.SetShortRangePotential(0,0,SingleGauss,pPars);
    for(unsigned uSor=0; uSor<NumSource; uSor++){
      Kitty.SetAnaSource(0,SourceSize[uSor],false);
      printf("Killing %u %u / %u %u\n",uPot,uSor,NumPots,NumSource);
      Kitty.KillTheCat();
      for(unsigned uMom=0; uMom<Kitty.GetNumMomBins(); uMom++){
        gCk[uPot][uSor].SetPoint(uMom,Kitty.GetMomentum(uMom),1.0*(Kitty.GetCorrFun(uMom)+0.0));
      }
    }
    double ScatLen;
    double EffRan;
    TH1F* hFit;
    TF1* fitSP;

    Eval_ScattParameters(Kitty, ScatLen, EffRan, hFit, fitSP);

    printf("f0, d0 = %.2f %.2f\n", ScatLen, EffRan);
  }

  for(unsigned uSor=0; uSor<NumSource; uSor++){
    if(SourceSize[uSor]==1 || SourceSize[uSor]==2 || SourceSize[uSor]==3 || SourceSize[uSor]==5){
      TString InputFileName = TString::Format("%s/OtherTasks/DDstar/corr_DDstar/%.0ffm.dat",GetFemtoOutputFolder(),SourceSize[uSor]);
      FILE *InFile;
      InFile = fopen(InputFileName.Data(), "r");
      if(!InFile){
          printf("          \033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", InputFileName.Data());
          return;
      }
      fseek ( InFile , 0 , SEEK_END );
      long EndPos;
      EndPos = ftell (InFile);
      fseek ( InFile , 0 , SEEK_SET );
      long CurPos;
      float Mom;
      float CkVal;
      float fDummy;
      unsigned uMom=0;
      while(!feof(InFile)){
        if(!fscanf(InFile,"%f %f %f",&Mom,&CkVal,&fDummy)){
            printf("\033[1;33mWARNING!\033[0m Possible bad input-file, error when reading from %s!\n",InputFileName.Data());
            continue;
        }
        gYuki[uSor].SetPoint(uMom++,Mom,CkVal);
      }
      fclose(InFile);
    }
    else{
      CATS Kitty;
      Kitty.SetMomBins(300,0,300);
      for(unsigned uMom=0; uMom<Kitty.GetNumMomBins(); uMom++){
        gYuki[uSor].SetPoint(uMom,Kitty.GetMomentum(uMom),0);
      }
    }
  }

  TFile fOutput(TString::Format("%s/OtherTasks/DDstar/TetsuoTest1.root",GetFemtoOutputFolder()),"recreate");
  for(unsigned uPot=0; uPot<NumPots; uPot++){
    for(unsigned uSor=0; uSor<NumSource; uSor++){
      gCk[uPot][uSor].Write();
    }
  }
  for(unsigned uSor=0; uSor<NumSource; uSor++){
    gYuki[uSor].Write();
  }

  for(unsigned uPot=0; uPot<NumPots; uPot++){
    delete [] gCk[uPot];
  }
  delete [] gCk;
  delete [] gYuki;
  delete [] SourceSize;
  delete [] par_V0;
  delete [] par_R;
}

//as a function of the reduced mass
void AV18peak(){
  const unsigned NumRedMass = 120;
  const double RedMassMin = 200;
  const double RedMassMax = 800;
  const unsigned kStep = 38;
  double* MomentumBins = new double[kStep+1];
  MomentumBins[0] = 0;
  for(unsigned uMom=1; uMom<=2; uMom++){
    MomentumBins[uMom] = MomentumBins[uMom-1]+4.;
  }
  for(unsigned uMom=3; uMom<=26; uMom++){
    MomentumBins[uMom] = MomentumBins[uMom-1]+1.;
  }
  for(unsigned uMom=27; uMom<=38; uMom++){
    MomentumBins[uMom] = MomentumBins[uMom-1]+4.;
  }
  //for(unsigned uMom=0; uMom<kStep; uMom++){
  //  printf("%u %.1f : %.1f\n",uMom,MomentumBins[uMom],MomentumBins[uMom+1]);
  //}
  //const double kMin = 0;
  //const double kMax = 80;


  TH1F* hRedMass_Peak = new TH1F("hRedMass_Peak","hRedMass_Peak",NumRedMass,RedMassMin,RedMassMax);
  TH1F** hCk = new TH1F* [NumRedMass];

  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());
  for(unsigned uMass=0; uMass<NumRedMass; uMass++){
    printf("\r\033[K  Progress %.0f%%",100.*double(uMass)/double(NumRedMass));
    cout << flush;
    double RedMass = hRedMass_Peak->GetBinCenter(uMass+1);
    hCk[uMass] = new TH1F(TString::Format("hCk_%.0fMeV",RedMass),TString::Format("hCk_%.0fMeV",RedMass),kStep,MomentumBins);
    CATS Kitty;
    Kitty.SetMomBins(kStep,MomentumBins);
    CATSparameters cPars(CATSparameters::tSource,1,true);
    cPars.SetParameter(0,1.2);
    Kitty.SetAnaSource(GaussSource, cPars);
    Kitty.SetUseAnalyticSource(true);
    Kitty.SetAutoNormSource(false);
    Kitty.SetMomentumDependentSource(false);
    Kitty.SetRedMass(RedMass);
    Kitty.SetQ1Q2(1);
    Kitty.SetPdgId(2212, 2212);
    Kitty.SetNumChannels(2);
    Kitty.SetSpin(0,0);
    Kitty.SetSpin(1,1);
    Kitty.SetNumPW(0,1);
    Kitty.SetNumPW(1,0);
    Kitty.SetChannelWeight(0, 0.25);
    Kitty.SetChannelWeight(1, 0.75);
    double PotPars1S0[8]={NN_AV18,v18_Coupled3P2,1,1,1,0,0,0};
    CATSparameters cPotPars1S0(CATSparameters::tPotential,8,true);
    cPotPars1S0.SetParameters(PotPars1S0);
    Kitty.SetShortRangePotential(0,0,fDlmPot,cPotPars1S0);

    Kitty.SetEpsilonConv(1e-8);
    Kitty.SetEpsilonProp(1e-8);
    Kitty.SetMaxRad(96);
    Kitty.SetMaxRho(32);
    Kitty.SetNotifications(CATS::nWarning);
    Kitty.KillTheCat();
    double PeakPosition = -1;
    double PeakValue = 0;
    for(unsigned uMom=0; uMom<kStep; uMom++){
      hCk[uMass]->SetBinContent(uMom+1,Kitty.GetCorrFun(uMom));
      if(Kitty.GetCorrFun(uMom)>PeakValue){
        PeakValue = Kitty.GetCorrFun(uMom);
        PeakPosition = Kitty.GetMomentum(uMom);
      }
    }
    hRedMass_Peak->SetBinContent(uMass+1,PeakPosition);
  }
  printf("\r\033[K  DONE!");
  TFile fOutput(TString::Format("%s/OtherTasks/AV18peak.root",GetFemtoOutputFolder()),"recreate");
  hRedMass_Peak->Write();
  for(unsigned uMass=0; uMass<NumRedMass; uMass++){
    hCk[uMass]->Write();
  }

  delete hRedMass_Peak;
  for(unsigned uMass=0; uMass<NumRedMass; uMass++){
    delete hCk[uMass];
  }
  delete [] hCk;
}

void TestDLM_Sample(){
  const unsigned NumBins = 10;
  const unsigned NumEntries = 1000000;
  DLM_Histo<float> dlmHist;
  dlmHist.SetUp(1);
  dlmHist.SetUp(0,10,0,NumBins);
  dlmHist.Initialize();
  for(unsigned uBin=0; uBin<NumBins; uBin++){
    dlmHist.SetBinContent(uBin,dlmHist.GetBinCenter(0,uBin));
  }
  TH1F* hSample = new TH1F("hSample","hSample",NumBins*10,0,10);
  double axis_val[3];
  for(unsigned uSample=0; uSample<NumEntries; uSample++){
    dlmHist.SampleYield(axis_val);
    hSample->Fill(axis_val[0]);
  }
  TFile fOutput(TString::Format("%s/OtherTasks/TestDLM_Sample.root",GetFemtoOutputFolder()),"recreate");
  hSample->Write();
  delete hSample;

}

void Test_div_any_th1f(){
  TH1F* h1 = new TH1F("h1","h1",64,0,64);
  TH1F* h2 = new TH1F("h2","h2",31,2,64);
  //TH1F* h2 = new TH1F("h2","h2",4,0,4);

  for(unsigned uBin=0; uBin<h1->GetNbinsX(); uBin++){
    h1->SetBinContent(uBin+1,h1->GetBinCenter(uBin+1));
  }
  for(unsigned uBin=0; uBin<h2->GetNbinsX(); uBin++){
    h2->SetBinContent(uBin+1,(h2->GetBinCenter(uBin+1))*1);
  }
  h1->Sumw2();
  h2->Sumw2();

  TH1F* hr = new TH1F("hr","hr",64,0,64);

  DivideAnyTH1Fs(h1,h2,hr,true);
  for(unsigned uBin=0; uBin<hr->GetNbinsX(); uBin++){
    printf("%u %e %e\n",uBin,hr->GetBinContent(uBin+1),hr->GetBinError(uBin+1));
  }

  TFile fOut(TString::Format("%s/OtherTasks/DivideAnyTH1Fs.root",GetFemtoOutputFolder()),"recreate");
  h1->Write();
  h2->Write();
  hr->Write();
}

void Test_align_th1f(){
  TH1F* h1 = new TH1F("h1","h1",128,0,128);
  TH1F* h2 = new TH1F("h2","h2",126,2,128);
  TH1F* h3 = new TH1F("h3","h3",128,1,129);
  TRandom3 rangen(11);
  for(unsigned uBin=0; uBin<h1->GetNbinsX(); uBin++){
    h1->SetBinContent(uBin+1,uBin>=2?(rangen.Integer(10)+10):0);
    if(uBin>=2){
      h2->SetBinContent(uBin+1-2,h1->GetBinContent(uBin+1));
    }
    if(uBin>=1){
      h3->SetBinContent(uBin+1-1,h1->GetBinContent(uBin+1));
    }
  }
  h1->Sumw2();
  h2->Sumw2();
  h3->Sumw2();
  h2->Rebin(2);
  h3->Rebin(2);
  TH1F* h4 = (TH1F*)h3->Clone("h4");

  AlignTH1Fs(h2,h4);

  TFile fOut(TString::Format("%s/OtherTasks/Test_align_th1f.root",GetFemtoOutputFolder()),"recreate");
  h1->Write();
  h2->Write();
  h3->Write();
  h4->Write();

}

void RoughPiPiPotScan(
  TNtuple* result,
  TH1F* hV0, TH1F* hmu0, TH1F* hV1, TH1F* hmu1
){
  double RedMass = (Mass_p*Mass_pic)/(Mass_p+Mass_pic);
  CATSparameters sPars(CATSparameters::tSource,1,true);
  sPars.SetParameter(0,1.0);
  CATS Kitty_SE;
  Kitty_SE.SetMomBins(30,0,120);
  Kitty_SE.SetAnaSource(GaussSource, sPars);
  Kitty_SE.SetUseAnalyticSource(true);
  Kitty_SE.SetQ1Q2(0);
  Kitty_SE.SetQuantumStatistics(false);
  Kitty_SE.SetRedMass( RedMass );
  Kitty_SE.SetNumChannels(1);
  Kitty_SE.SetNumPW(0,1);
  Kitty_SE.SetSpin(0,0);
  Kitty_SE.SetChannelWeight(0, 1.);
  Kitty_SE.SetEpsilonConv(4e-8);
  Kitty_SE.SetEpsilonProp(4e-8);
  //Kitty_SE.SetMaxRad(96);
  //Kitty_SE.SetMaxRho(32);
  Kitty_SE.SetNotifications(CATS::nSilent);
  CATSparameters pPars(CATSparameters::tPotential,4,true);
  Kitty_SE.SetShortRangePotential(0,0,DoubleGaussSum,pPars);
  Kitty_SE.SetShortRangePotential(0,0,0,0);
  Kitty_SE.SetShortRangePotential(0,0,1,1);
  Kitty_SE.SetShortRangePotential(0,0,2,0);
  Kitty_SE.SetShortRangePotential(0,0,3,1);
  Kitty_SE.SetAnaSource(0,1.0);
  Kitty_SE.KillTheCat();

  float BUFFER[6];
  for(unsigned uV0=0; uV0<hV0->GetNbinsX(); uV0++){
    for(unsigned umu0=0; umu0<hmu0->GetNbinsX(); umu0++){
      for(unsigned uV1=0; uV1<hV1->GetNbinsX(); uV1++){
        for(unsigned umu1=0; umu1<hmu1->GetNbinsX(); umu1++){
          BUFFER[0] = hV0->GetBinCenter(uV0+1);
          BUFFER[1] = hmu0->GetBinCenter(umu0+1);
          BUFFER[2] = hV1->GetBinCenter(uV1+1);
          BUFFER[3] = hmu1->GetBinCenter(umu1+1);
          if(uV1==0 && umu1==0) printf("%.3f %.3f\n",BUFFER[0],BUFFER[1]);

          Kitty_SE.SetShortRangePotential(0,0,0,BUFFER[0]);
          Kitty_SE.SetShortRangePotential(0,0,1,BUFFER[1]);
          Kitty_SE.SetShortRangePotential(0,0,2,BUFFER[2]);
          Kitty_SE.SetShortRangePotential(0,0,3,BUFFER[3]);
          Kitty_SE.KillTheCat();

          TH1F* hDummy; TF1* fDummy;
          double df0,dd0;
          //potential with silly phase shifts
          if(!Eval_ScattParameters(Kitty_SE,df0,dd0,hDummy,fDummy)){
            df0=0;
            dd0=0;
          }
          BUFFER[4] = df0;
          BUFFER[5] = dd0;
          result->Fill(BUFFER);
          if(hDummy) delete hDummy;
          if(fDummy) delete fDummy;
        }
      }
    }
  }
}


//The goal is:
//1) f0 = 0.1248 +/- 0.0011 MeV
//2) f0 = -0.1810 +/- 0.0085 MeV

//Solution:
//1): V0 = 203.1 mu0 = 0.6203 V1 = -1106 mu1 = 0.4156


//COUNTING FROM 1 AS ROOT
void RoughPiPiPotScan(int V0_first, int V0_last){

  //printf("The goal is:\n");
  //printf("1) f0 = %.4f +/- %.4f MeV\n",0.0883/Mass_pic*NuToFm,0.0008/Mass_pic*NuToFm);
  //printf("2) f0 = %.4f +/- %.4f MeV\n",-0.128/Mass_pic*NuToFm,0.006/Mass_pic*NuToFm);

  //originally I started with 16
  const int BIN = 16;

  TFile fOutput(TString::Format("%s/OtherTasks/RoughPiPiPotScan/fOutput%i_%i_%i.root",
  GetFemtoOutputFolder(),BIN,V0_first,V0_last),"recreate");
  TNtuple* ntPiPi = new TNtuple("ntPiPi", "ntPiPi",
  "V0:mu0:V1:mu1:f0:d0");

  //initial scan
  /*
  int V0_Bin=BIN; float V0_Min=-500; float V0_Max=1500;
  int mu0_Bin=BIN; float mu0_Min=0.1; float mu0_Max=1.2;
  int V1_Bin=BIN; float V1_Min=-1500; float V1_Max=500;
  int mu1_Bin=BIN; float mu1_Min=0.1; float mu1_Max=1.2;
  */
  /*
  //fine tune (1) for the attractive interaction
  int V0_Bin=BIN; float V0_Min=-500; float V0_Max=500;
  int mu0_Bin=BIN; float mu0_Min=0.1; float mu0_Max=0.7;
  int V1_Bin=BIN; float V1_Min=-1500; float V1_Max=-1000;
  int mu1_Bin=BIN; float mu1_Min=0.1; float mu1_Max=0.5;
  */
/*
  //fine tune (2) for the attractive interaction
  int V0_Bin=BIN; float V0_Min=0; float V0_Max=500;
  int mu0_Bin=BIN; float mu0_Min=0.55; float mu0_Max=0.7;
  int V1_Bin=BIN; float V1_Min=-1500; float V1_Max=-900;
  int mu1_Bin=BIN; float mu1_Min=0.32; float mu1_Max=0.5;
*/
  //fine tune (1) for the repulsive interaction
  int V0_Bin=BIN; float V0_Min=400; float V0_Max=1600;
  int mu0_Bin=BIN; float mu0_Min=0.2; float mu0_Max=0.6;
  int V1_Bin=BIN; float V1_Min=-500; float V1_Max=1000;
  int mu1_Bin=BIN; float mu1_Min=0.0; float mu1_Max=0.8;


  TH1F* hV0_tot = new TH1F("hV0_tot","hV0_tot",V0_Bin,V0_Min,V0_Max);

  TH1F* hV0 = new TH1F("hV0","hV0",V0_last-V0_first+1,
  hV0_tot->GetXaxis()->GetBinLowEdge(V0_first),
  hV0_tot->GetXaxis()->GetBinUpEdge(V0_last));
  TH1F* hmu0 = new TH1F("hmu0","hmu0",mu0_Bin,mu0_Min,mu0_Max);
  TH1F* hV1 = new TH1F("hV1","hV1",V1_Bin,V1_Min,V1_Max);
  TH1F* hmu1 = new TH1F("hmu1","hmu1",mu1_Bin,mu1_Min,mu1_Max);

  RoughPiPiPotScan(ntPiPi,hV0,hmu0,hV1,hmu1);

  ntPiPi->Write();

  delete ntPiPi;
}

void pi_proton(){
  const double Weight0 = 0.5;
  double RedMass = (Mass_p*Mass_pic)/(Mass_p+Mass_pic);

  const unsigned NumBins = 160;
  const double kMin = 0;
  const double kMax = 320;
  const double SourceSize = 1.4;

  CATSparameters sPars(CATSparameters::tSource,1,true);
  sPars.SetParameter(0,SourceSize);
  CATS Kitty_SE;
  Kitty_SE.SetMomBins(NumBins,kMin,kMax);
  Kitty_SE.SetAnaSource(GaussSource, sPars);
  Kitty_SE.SetUseAnalyticSource(true);
  Kitty_SE.SetQ1Q2(-1);
  Kitty_SE.SetQuantumStatistics(false);
  Kitty_SE.SetRedMass( RedMass );
  Kitty_SE.SetNumChannels(2);
  Kitty_SE.SetNumPW(0,1);
  Kitty_SE.SetNumPW(1,1);
  Kitty_SE.SetSpin(0,0);
  Kitty_SE.SetSpin(1,1);
  Kitty_SE.SetChannelWeight(0, Weight0);
  Kitty_SE.SetChannelWeight(1, 1.-Weight0);
  Kitty_SE.SetEpsilonConv(5e-9);
  Kitty_SE.SetEpsilonProp(5e-9);
  Kitty_SE.SetMaxRad(96);
  Kitty_SE.SetMaxRho(32);
  //Kitty_SE.SetNotifications(CATS::nSilent);
  CATSparameters pPars0(CATSparameters::tPotential,4,true);
  pPars0.SetParameter(0,1.973303e+02);
  pPars0.SetParameter(1,6.187704e-01);
  pPars0.SetParameter(2,-1.100288e+03);
  pPars0.SetParameter(3,4.147391e-01);
  Kitty_SE.SetShortRangePotential(0,0,DoubleGaussSum,pPars0);

  CATSparameters pPars1(CATSparameters::tPotential,4,true);
  pPars1.SetParameter(0,7.495875e+02);
  pPars1.SetParameter(1,3.944096e-01);
  pPars1.SetParameter(2,6.549811e+02);
  pPars1.SetParameter(3,4.493379e-01);
  //CHANNEL IS WRONG!!!
  Kitty_SE.SetShortRangePotential(0,0,DoubleGaussSum,pPars1);

  Kitty_SE.KillTheCat();

  TFile fOutput(TString::Format("%s/OtherTasks/RoughPiPiPotScan/pim_proton.root",GetFemtoOutputFolder()),"recreate");

  TH1F* hCk = new TH1F("hCk","hCk",NumBins,kMin,kMax);
  TH1F* hCkGeV = new TH1F("hCkGeV","hCkGeV",NumBins,kMin/1000.,kMax/1000.);
  for(unsigned uBin=0; uBin<NumBins; uBin++){
    hCk->SetBinContent(uBin+1,Kitty_SE.GetCorrFun(uBin));
    hCkGeV->SetBinContent(uBin+1,Kitty_SE.GetCorrFun(uBin));
  }
  hCk->Write();
  hCkGeV->Write();
  delete hCk;
  delete hCkGeV;
}

void Jaime_test1(){
  unsigned numBins = 100;
  double kmin = 0;
  double kmax = 400;
  double mass1 = 938;
  double mass2 = 1116;
  CATS Kitty;
  Kitty.SetMomBins(numBins, kmin, kmax);
  // Create the object as a pointer:
  //CATS* Kitty;
  //Kitty = new CATS();
  //Kitty->SetMomBins(numBins, kmin, kmax);
  // definition of the source:
  CATSparameters source_func(CATSparameters::tSource, 1, true);
  Kitty.SetAnaSource(GaussSource, source_func);
  Kitty.SetAnaSource(0, 1.2);
  Kitty.SetQ1Q2(0);
  Kitty.SetQuantumStatistics(0);
  Kitty.SetRedMass(mass1 * mass2 /(mass1 + mass2));
  Kitty.SetNumChannels(2);
  Kitty.SetNumPW(0,1);
  Kitty.SetNumPW(1,1);
  Kitty.SetSpin(0, 0);
  Kitty.SetSpin(1, 1);
  Kitty.SetChannelWeight(0, 1./4.);
  Kitty.SetChannelWeight(1, 3./4.);
  CATSparameters potential_spin0(CATSparameters::tPotential, 2, true);
  potential_spin0.SetParameter(0, 0);
  Kitty.SetShortRangePotential(0, 0, Gaussian, potential_spin0);
  //Kitty.SetShortRangePotential(0, 0, 0, 0);
  //Kitty.SetShortRangePotential(0, 0, 1, 1);
  //Kitty.SetShortRangePotential(0, 0, 0, 0);
  //Kitty.SetShortRangePotential(0, 0, 0, 0);
  //Kitty.SetShortRangePotential(0, 0, 0, 0);
  //Kitty.SetShortRangePotential(0, 0, 0, 0);
  Kitty.KillTheCat();
  TH1F* his = new TH1F("his","his", numBins, kmin, kmax);
  for (int n=0; n<numBins; ++n){
      double kstar = his->GetBinCenter(n+1);
      his->SetBinContent(n + 1, Kitty.EvalCorrFun(kstar));
  }
  //his->Draw();

  TFile* output = new TFile(TString::Format("%s/OtherTasks/Jaime_test1.root",GetFemtoOutputFolder()), "RECREATE");
  //output->cd();
  his->Write();

}

void pp_for_rock(const double SourceSize){
  const unsigned NumBins = 400;
  const double kMin = 0;
  const double kMax = 400;
  const bool RSM = true;
  CATS AB_pp;
  AB_pp.SetMomBins(NumBins,kMin,kMax);
  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());
  AnalysisObject.SetUpCats_pp(AB_pp,"AV18",RSM?"McGauss_ResoTM":"Gauss",0,RSM?202:0);
  AB_pp.SetAnaSource(0,SourceSize);
  AB_pp.KillTheCat();
  TH1F* hCk_pp = new TH1F("hCk_pp","hCk_pp",NumBins,kMin,kMax);
  TH1F* hSr_pp = new TH1F("hSr_pp","hSr_pp",512,0,64);
  for(unsigned uBin=0; uBin<NumBins; uBin++){
    hCk_pp->SetBinContent(uBin+1,AB_pp.GetCorrFun(uBin));
  }
  for(unsigned uRad=0; uRad<hSr_pp->GetNbinsX(); uRad++){
    double rad = hSr_pp->GetBinCenter(uRad+1);
    hSr_pp->SetBinContent(uRad+1,AB_pp.EvaluateTheSource(0,rad,0));
    hSr_pp->SetBinError(uRad+1,AB_pp.EvaluateTheSource(0,rad,0)*0.01);
  }
  double reff = Get_reff(hSr_pp);
  printf("r_core = %.3f\n",SourceSize);
  printf("r_eff = %.3f\n",reff);

  TFile fOutput(TString::Format("%s/OtherTasks/pp_for_rock_%s%.2f.root",GetFemtoOutputFolder(),RSM?"rcore":"reff",SourceSize), "RECREATE");
  hCk_pp->Write();
  hSr_pp->Write();
  delete hCk_pp;
  delete hSr_pp;
}

void ppSource_bugHunting(const bool IncludeBug){

  const double CoreSize = 0.95;

  DLM_CleverMcLevyResoTM MagicSource;

  MagicSource.InitStability(1,2-1e-6,2+1e-6);
  MagicSource.InitScale(38,0.15,2.0);
  MagicSource.InitRad(257*2,0,64);
  MagicSource.InitType(2);

  MagicSource.SetUpReso(0,0.6422);
  MagicSource.SetUpReso(1,0.6422);
  MagicSource.InitNumMcIter(1000000);

  const double k_CutOff = 200;
  Float_t k_D;
  Float_t fP1;
  Float_t fP2;
  Float_t fM1;
  Float_t fM2;
  Float_t Tau1;
  Float_t Tau2;
  Float_t AngleRcP1;
  Float_t AngleRcP2;
  Float_t AngleP1P2;
  DLM_Random RanGen(11);
  double RanVal1;
  double RanVal2;
  double RanVal3;

  TFile* F_EposDisto_p_pReso;
  F_EposDisto_p_pReso = new TFile(TString::Format("%s/CatsFiles/Source/EposAngularDist/EposDisto_p_pReso.root",GetCernBoxDimi()));
  //printf("F_EposDisto_p_pReso=%p (%s)\n",F_EposDisto_p_pReso,TString(CatsFilesFolder[0]+"/Source/EposAngularDist/Epos3body_p_pReso_3body.root").Data());
  TNtuple* T_EposDisto_p_pReso = (TNtuple*)F_EposDisto_p_pReso->Get("InfoTuple_ClosePairs");
  unsigned N_EposDisto_p_pReso = T_EposDisto_p_pReso->GetEntries();
  T_EposDisto_p_pReso->SetBranchAddress("k_D",&k_D);
  T_EposDisto_p_pReso->SetBranchAddress("P1",&fP1);
  T_EposDisto_p_pReso->SetBranchAddress("P2",&fP2);
  T_EposDisto_p_pReso->SetBranchAddress("M1",&fM1);
  T_EposDisto_p_pReso->SetBranchAddress("M2",&fM2);
  T_EposDisto_p_pReso->SetBranchAddress("Tau1",&Tau1);
  T_EposDisto_p_pReso->SetBranchAddress("Tau2",&Tau2);
  T_EposDisto_p_pReso->SetBranchAddress("AngleRcP1",&AngleRcP1);
  T_EposDisto_p_pReso->SetBranchAddress("AngleRcP2",&AngleRcP2);
  T_EposDisto_p_pReso->SetBranchAddress("AngleP1P2",&AngleP1P2);
  for(unsigned uEntry=0; uEntry<N_EposDisto_p_pReso; uEntry++){
      T_EposDisto_p_pReso->GetEntry(uEntry);
      Tau1 = 0;
      Tau2 = 1.65;
      fM2 = 1362;
      if(k_D>k_CutOff) continue;
      RanVal1 = RanGen.Exponential(fM2/(fP2*Tau2));
      if(IncludeBug){
        MagicSource.AddBGT_PR(RanVal1,-cos(AngleRcP2));
        MagicSource.AddBGT_RP(RanVal1,cos(AngleRcP2));
      }
      else{
        MagicSource.AddBGT_PR(RanVal1,cos(AngleRcP2));
        MagicSource.AddBGT_RP(RanVal1,-cos(AngleRcP2));
      }
  }
  delete F_EposDisto_p_pReso;

  TFile* F_EposDisto_pReso_dReso;
  F_EposDisto_pReso_dReso = new TFile(TString::Format("%s/CatsFiles/Source/EposAngularDist/EposDisto_pReso_dReso.root",GetCernBoxDimi()));
  //printf("F_EposDisto_p_pReso=%p (%s)\n",F_EposDisto_p_pReso,TString(CatsFilesFolder[0]+"/Source/EposAngularDist/Epos3body_pReso_dReso_3body.root").Data());
  TNtuple* T_EposDisto_pReso_dReso = (TNtuple*)F_EposDisto_pReso_dReso->Get("InfoTuple_ClosePairs");
  unsigned N_EposDisto_pReso_dReso = T_EposDisto_pReso_dReso->GetEntries();
  T_EposDisto_pReso_dReso->SetBranchAddress("k_D",&k_D);
  T_EposDisto_pReso_dReso->SetBranchAddress("P1",&fP1);
  T_EposDisto_pReso_dReso->SetBranchAddress("P2",&fP2);
  T_EposDisto_pReso_dReso->SetBranchAddress("M1",&fM1);
  T_EposDisto_pReso_dReso->SetBranchAddress("M2",&fM2);
  T_EposDisto_pReso_dReso->SetBranchAddress("Tau1",&Tau1);
  T_EposDisto_pReso_dReso->SetBranchAddress("Tau2",&Tau2);
  T_EposDisto_pReso_dReso->SetBranchAddress("AngleRcP1",&AngleRcP1);
  T_EposDisto_pReso_dReso->SetBranchAddress("AngleRcP2",&AngleRcP2);
  T_EposDisto_pReso_dReso->SetBranchAddress("AngleP1P2",&AngleP1P2);
  for(unsigned uEntry=0; uEntry<N_EposDisto_pReso_dReso; uEntry++){
      T_EposDisto_pReso_dReso->GetEntry(uEntry);
      Tau1 = 1.65;
      Tau2 = 1.65;
      fM1 = 1362;
      fM2 = 1362;
      if(k_D>k_CutOff) continue;
      RanVal1 = RanGen.Exponential(fM1/(fP1*Tau1));
      RanVal2 = RanGen.Exponential(fM2/(fP2*Tau2));
      MagicSource.AddBGT_RR(RanVal1,cos(AngleRcP1),RanVal2,cos(AngleRcP2),cos(AngleP1P2));
  }
  delete F_EposDisto_pReso_dReso;


  const unsigned NumSourceBins = 128;
  const double rMin = 0;
  const double rMax = 16;
  TFile* fOutput = new TFile(TString::Format("%s/OtherTasks/ppSource_bugHunting/fOutput%s_pp_%.2f.root",GetFemtoOutputFolder(),IncludeBug?"BUG":"",
  CoreSize),"recreate");
  TH1F* hSource = new TH1F("hSource","hSource",NumSourceBins,rMin,rMax);

  //fill the histo fro the source
  for(unsigned uBin=0; uBin<NumSourceBins; uBin++){
    //get the x-axis (r value) of the current bin
    double xaxis = hSource->GetBinCenter(uBin+1);
    //an array for the parameters, [0] is source size, [1] is == 2 (for a Gaussian)
    double parameters[2];
    parameters[0] = CoreSize;
    parameters[1] = 2.0;
    double SourceValue = MagicSource.RootEval(&xaxis, parameters);
    hSource->SetBinContent(uBin+1,SourceValue);
    //infinite errors for now
    hSource->SetBinError(uBin+1,1000.);
  }
//printf("4\n");
  //idea: fit the source distribution only in a range around its peak
  //to do this: silly idea: put very large uncertainties in the bins outside of this range
  //we can get this range automatically, by evaluating the central (median) integral of the source distribution
  //with this set up, we fit the 68% most central yield of the source distribution
  double lowerlimit;
  double upperlimit;
  GetCentralInterval(*hSource, 0.84, lowerlimit, upperlimit, true);
  unsigned lowerbin = hSource->FindBin(lowerlimit);
  unsigned upperbin = hSource->FindBin(upperlimit);
  for(unsigned uBin=lowerbin; uBin<=upperbin; uBin++){
    hSource->SetBinError(uBin+1,0.01);
  }

  printf("Core size of %.3f fm\n",CoreSize);
  printf("The fit will be performed in the range [%.2f, %.2f] fm\n",lowerlimit,upperlimit);
  //fyi, GaussSourceTF1 is in DLM_Source.h if you want to check it out.
  TF1* fSource = new TF1("fSource",GaussSourceTF1,rMin,rMax,1);
  fSource->SetParameter(0,CoreSize);
  fSource->SetParLimits(0,CoreSize*0.5,CoreSize*2.0);
  hSource->Fit(fSource,"S, N, R, M");
  printf("The effective Gaussian size is %.3f +/- %.3f fm\n",fSource->GetParameter(0),fSource->GetParError(0));

  //get rid of weird plotting
  for(unsigned uBin=0; uBin<NumSourceBins; uBin++){
    hSource->SetBinError(uBin+1,0.01);
  }
  hSource->Write();
  fSource->Write();

  delete hSource;
  delete fSource;
  delete fOutput;

}


void PlugInWaveFunction(){
  //to be used for binning. In this example: 5 MeV bins
  const unsigned NumMomBins = 40;
  const double kstar_min = 0;
  const double kstar_max = 200;

  //later to be used for the radial relation of the wave function
  //this is important to be rather fine. The maximum value should be such that
  //the tail of the source function is negligible. As a rule of thumb, for
  //a Gaussian source this corresponds to c.a. 10 times the r0 value.
  const unsigned NumRadBins = 1000;
  const double r_min = 0;
  const double r_max = 50;

  //the source size
  const double r0 = 3.0;

  //some "container" for the source parameters.
  //here trivial, as we will use 1D Gauss source (single parameter)
  CATSparameters cPars(CATSparameters::tSource,1,true);
  CATS Kitty;
  Kitty.SetMomBins(NumMomBins,kstar_min,kstar_max);
  //definition of GaussSource in DLM_Source.h
  Kitty.SetAnaSource(GaussSource, cPars);
  Kitty.SetAnaSource(0,r0);
  Kitty.SetUseAnalyticSource(true);
  Kitty.SetMomentumDependentSource(false);

  //charge, q1*q2
  Kitty.SetQ1Q2(0);
  Kitty.SetQuantumStatistics(false);
  Kitty.SetRedMass( (Mass_p*Mass_Xim)/(Mass_p+Mass_Xim) );

  //in case of different spin/isospin channels, they can be separately defined here
  //in this example, we only define a single effective channel, that contains
  //a single partial wave (PW) that is deviating from the asymptotic form
  Kitty.SetNumChannels(1);
  Kitty.SetNumPW(0,1);
  //dummy value in this case, only relevant in case of identical particle
  Kitty.SetSpin(0,0);
  //as we have 1 channel, we take it with weight of 1 (100%)
  Kitty.SetChannelWeight(0,1.);

  //the wave function has to be passed by using some custom made histograms

  //the histograms should be 2D, with one axis corresponding to momentum k (x)
  //and the other to the distance r (y). The value (z-axis) corresponds to the
  //value of the radial partial wave functions (k,r) for the corresponding channel
  //in this example: single channel and single PW, thus a single wave function is needed
  DLM_Histo<complex<double>> histWF;
  //two dimensions
  histWF.SetUp(2);
  //how many bins/range in each dimension
  histWF.SetUp(0,NumMomBins,kstar_min,kstar_max);
  histWF.SetUp(1,NumRadBins,r_min,r_max);

  histWF.Initialize();

  //we can also provide phase shifts, but this is not needed, so concider
  //this histogram a dummy
  DLM_Histo<complex<double>> histPS(histWF);

  //example how to set the bin values of the wave function histogram
  //in this example, we will simply set it equal to the Coulomb solution (trivial),
  //so we should get a flat correlation
  unsigned axis[2];
  /*
  for(axis[0]=0; axis[0]<histWF.GetNbins(0); axis[0]++){//over kstar
    for(axis[1]=0; axis[1]<histWF.GetNbins(1); axis[1]++){//over radial dependence
      //this is a trick to evaluate the s coulomb partial wave using CATS
      //outside of this example, here is where the wave function needs to be evaluated
      //at a fixed kstar and as a function of its radial relation
      //complex<double> CPW = Kitty.EvalReferenceRadialWF(
      //              axis[0],0,histWF.GetBinCenter(1,axis[1]),false);
//complex<double> CPW = Kitty.ReferencePartialWave(
//                        histWF.GetBinCenter(1,axis[1]),
//                        histWF.GetBinCenter(0,axis[0]), 0, -1);
      //complex<double> RWF = Kitty.GetRadialWaveFunction(axis[0], 0, 0, axis[1]);
      //printf("%e vs %e\n",abs(CPW),abs(RWF));

double Rho = histWF.GetBinCenter(1,axis[1])*histWF.GetBinCenter(0,axis[0]);
complex<double> CPW = histWF.GetBinCenter(1,axis[1])*sin(Rho)/Rho;
//printf("%f %f\n", histWF.GetBinCenter(1,axis[1]), histWF.GetBinCenter(0,axis[0]));

      histWF.SetBinContent(axis,CPW);
    }
  }
*/
  //Kitty.SetExternalWaveFunction(0,0,histWF,histPS);
  Kitty.KillTheCat();

  TH1F* h_CorrFunction = new TH1F("h_CorrFunction","h_CorrFunction",NumMomBins,kstar_min,kstar_max);
  for(unsigned uBin=0; uBin<NumMomBins; uBin++){
    h_CorrFunction->SetBinContent(uBin+1,Kitty.GetCorrFun(uBin));
  }
  //TFile fCatsExample("fCatsExample.root","recreate");
  //h_CorrFunction->Write();
}

void pd_delayed(const double& delay, const bool& reso){

    const double CoreSize = 0.97;

  //  TH1F* hFemtoMom_p = TH1F("hFemtoMom_p","hFemtoMom_p",128,0,4096);
  //  TH1F* hFemtoMom_pReso = TH1F("hFemtoMom_pReso","hFemtoMom_pReso",128,0,4096);
  //  TH1F* hFemtoMom_d = TH1F("hFemtoMom_d","hFemtoMom_d",128,0,4096);

    TH1F* p_Omega_AngleRcP2 = new TH1F("p_Omega_AngleRcP2","p_Omega_AngleRcP2",32,0,Pi);
    TH1F* pReso_Omega_AngleRcP1 = new TH1F("pReso_Omega_AngleRcP1","pReso_Omega_AngleRcP1",32,0,Pi);
    TH1F* pReso_dReso_AngleRcP1 = new TH1F("pReso_dReso_AngleRcP1","pReso_dReso_AngleRcP1",32,0,Pi);
    TH1F* pReso_dReso_AngleRcP2 = new TH1F("pReso_dReso_AngleRcP2","pReso_dReso_AngleRcP2",32,0,Pi);
    TH1F* pReso_dReso_AngleP1P2 = new TH1F("pReso_dReso_AngleP1P2","pReso_dReso_AngleP1P2",32,0,Pi);

    DLM_CleverMcLevyResoTM MagicSource;

    MagicSource.InitStability(1,2-1e-6,2+1e-6);
    MagicSource.InitScale(38,0.15,2.0);
    MagicSource.InitRad(257*2,0,64);
    MagicSource.InitType(2);

    MagicSource.SetUpReso(0,0.6422*reso);
    MagicSource.SetUpReso(1,1.0*reso);
    MagicSource.InitNumMcIter(1000000);

    const double k_CutOff = 200;
    Float_t k_D;
    Float_t fP1;
    Float_t fP2;
    Float_t fM1;
    Float_t fM2;
    Float_t Tau1;
    Float_t Tau2;
    Float_t AngleRcP1;
    Float_t AngleRcP2;
    Float_t AngleP1P2;
    DLM_Random RanGen(11);
    double RanVal1;
    double RanVal2;
    double RanVal3;

    TFile* F_EposDisto_p_dReso;
    F_EposDisto_p_dReso = new TFile(TString::Format("%s/CatsFiles/Source/EposAngularDist/Bhawani_p_Omega.root",GetCernBoxDimi()));
    TNtuple* T_EposDisto_p_dReso = (TNtuple*)F_EposDisto_p_dReso->Get("InfoTuple_ClosePairs");
    unsigned N_EposDisto_p_dReso = T_EposDisto_p_dReso->GetEntries();
    T_EposDisto_p_dReso->SetBranchAddress("k_D",&k_D);
    T_EposDisto_p_dReso->SetBranchAddress("P1",&fP1);
    T_EposDisto_p_dReso->SetBranchAddress("P2",&fP2);
    T_EposDisto_p_dReso->SetBranchAddress("M1",&fM1);
    T_EposDisto_p_dReso->SetBranchAddress("M2",&fM2);
    T_EposDisto_p_dReso->SetBranchAddress("Tau1",&Tau1);
    T_EposDisto_p_dReso->SetBranchAddress("Tau2",&Tau2);
    T_EposDisto_p_dReso->SetBranchAddress("AngleRcP1",&AngleRcP1);
    T_EposDisto_p_dReso->SetBranchAddress("AngleRcP2",&AngleRcP2);
    T_EposDisto_p_dReso->SetBranchAddress("AngleP1P2",&AngleP1P2);
    for(unsigned uEntry=0; uEntry<N_EposDisto_p_dReso; uEntry++){
        T_EposDisto_p_dReso->GetEntry(uEntry);
        Tau1 = 0;
        Tau2 = delay;
        fM2 = Mass_d*1.005;
        if(k_D>k_CutOff) continue;
        RanVal2 = fP2/fM2*Tau2;
        if(RanVal2>1000){
          printf("AddBGT_PR %.2e %.2e\n",RanVal2,cos(AngleRcP2));
          usleep(100e3);
        }
        MagicSource.AddBGT_PR(RanVal2,cos(AngleRcP2));
        p_Omega_AngleRcP2->Fill(AngleRcP2);
    }
    delete F_EposDisto_p_dReso;

    TFile* F_EposDisto_pReso_d;
    F_EposDisto_pReso_d = new TFile(TString::Format("%s/CatsFiles/Source/EposAngularDist/EposDisto_pReso_Omega.root",GetCernBoxDimi()));
    //printf("F_EposDisto_p_pReso=%p (%s)\n",F_EposDisto_p_pReso,TString(CatsFilesFolder[0]+"/Source/EposAngularDist/Epos3body_pReso_d_3body.root").Data());
    TNtuple* T_EposDisto_pReso_d = (TNtuple*)F_EposDisto_pReso_d->Get("InfoTuple_ClosePairs");
    unsigned N_EposDisto_pReso_d = T_EposDisto_pReso_d->GetEntries();
    T_EposDisto_pReso_d->SetBranchAddress("k_D",&k_D);
    T_EposDisto_pReso_d->SetBranchAddress("P1",&fP1);
    T_EposDisto_pReso_d->SetBranchAddress("P2",&fP2);
    T_EposDisto_pReso_d->SetBranchAddress("M1",&fM1);
    T_EposDisto_pReso_d->SetBranchAddress("M2",&fM2);
    T_EposDisto_pReso_d->SetBranchAddress("Tau1",&Tau1);
    T_EposDisto_pReso_d->SetBranchAddress("Tau2",&Tau2);
    T_EposDisto_pReso_d->SetBranchAddress("AngleRcP1",&AngleRcP1);
    T_EposDisto_pReso_d->SetBranchAddress("AngleRcP2",&AngleRcP2);
    T_EposDisto_pReso_d->SetBranchAddress("AngleP1P2",&AngleP1P2);
    for(unsigned uEntry=0; uEntry<N_EposDisto_pReso_d; uEntry++){
        T_EposDisto_pReso_d->GetEntry(uEntry);
        Tau1 = 1.65;
        Tau2 = 0;
        fM1 = 1362;
        if(k_D>k_CutOff) continue;
        RanVal1 = RanGen.Exponential(fM1/(fP1*Tau1));
        MagicSource.AddBGT_RP(RanVal1,cos(AngleRcP1));
        pReso_Omega_AngleRcP1->Fill(AngleRcP1);
        //if(RanVal1>1000){
        //  printf("AddBGT_PR %.2e %.2e\n",RanVal1,cos(RanVal1));
        //  usleep(100e3);
        //}

    }
    delete F_EposDisto_pReso_d;


    TFile* F_EposDisto_pReso_dReso;
    F_EposDisto_pReso_dReso = new TFile(TString::Format("%s/CatsFiles/Source/EposAngularDist/EposDisto_pReso_Omega.root",GetCernBoxDimi()));
    //printf("F_EposDisto_p_pReso=%p (%s)\n",F_EposDisto_p_pReso,TString(CatsFilesFolder[0]+"/Source/EposAngularDist/Epos3body_pReso_dReso_3body.root").Data());
    TNtuple* T_EposDisto_pReso_dReso = (TNtuple*)F_EposDisto_pReso_dReso->Get("InfoTuple_ClosePairs");
    unsigned N_EposDisto_pReso_dReso = T_EposDisto_pReso_dReso->GetEntries();
    T_EposDisto_pReso_dReso->SetBranchAddress("k_D",&k_D);
    T_EposDisto_pReso_dReso->SetBranchAddress("P1",&fP1);
    T_EposDisto_pReso_dReso->SetBranchAddress("P2",&fP2);
    T_EposDisto_pReso_dReso->SetBranchAddress("M1",&fM1);
    T_EposDisto_pReso_dReso->SetBranchAddress("M2",&fM2);
    T_EposDisto_pReso_dReso->SetBranchAddress("Tau1",&Tau1);
    T_EposDisto_pReso_dReso->SetBranchAddress("Tau2",&Tau2);
    T_EposDisto_pReso_dReso->SetBranchAddress("AngleRcP1",&AngleRcP1);
    T_EposDisto_pReso_dReso->SetBranchAddress("AngleRcP2",&AngleRcP2);
    T_EposDisto_pReso_dReso->SetBranchAddress("AngleP1P2",&AngleP1P2);
    for(unsigned uEntry=0; uEntry<N_EposDisto_pReso_dReso; uEntry++){
        T_EposDisto_pReso_dReso->GetEntry(uEntry);
        Tau1 = 1.65;
        Tau2 = delay;
        fM1 = 1362;
        fM2 = Mass_d*1.005;
        if(k_D>k_CutOff) continue;
        RanVal1 = RanGen.Exponential(fM1/(fP1*Tau1));
        RanVal2 = fP2/fM2*Tau2;
        MagicSource.AddBGT_RR(RanVal1,cos(AngleRcP1),RanVal2,cos(AngleRcP2),cos(AngleP1P2));
        pReso_dReso_AngleRcP1->Fill(AngleRcP1);
        pReso_dReso_AngleRcP2->Fill(AngleRcP2);
        pReso_dReso_AngleP1P2->Fill(AngleP1P2);
        //if(RanVal1>1000||RanVal2>1000){
        //  printf("AddBGT_RR %.2e %.2e %.2e %.2e %.2e\n",RanVal1,cos(AngleRcP1),RanVal2,cos(AngleRcP2),cos(AngleP1P2));
        //  usleep(100e3);
        //}
    }
    delete F_EposDisto_pReso_dReso;


    const unsigned NumSourceBins = 128;
    const double rMin = 0;
    const double rMax = 16;
    TFile* fOutput = new TFile(TString::Format("%s/OtherTasks/pd_delayed/fOutput_PR%i_DD%.0f.root",GetFemtoOutputFolder(),reso,delay),"recreate");
    TH1F* hSource = new TH1F("hSource","hSource",NumSourceBins,rMin,rMax);

    //fill the histo fro the source
    for(unsigned uBin=0; uBin<NumSourceBins; uBin++){
      //get the x-axis (r value) of the current bin
      double xaxis = hSource->GetBinCenter(uBin+1);
      //an array for the parameters, [0] is source size, [1] is == 2 (for a Gaussian)
      double parameters[2];
      parameters[0] = CoreSize;
      parameters[1] = 2.0;
      double SourceValue = MagicSource.RootEval(&xaxis, parameters);
      hSource->SetBinContent(uBin+1,SourceValue);
      //infinite errors for now
      hSource->SetBinError(uBin+1,1000.);
    }
  //printf("4\n");
    //idea: fit the source distribution only in a range around its peak
    //to do this: silly idea: put very large uncertainties in the bins outside of this range
    //we can get this range automatically, by evaluating the central (median) integral of the source distribution
    //with this set up, we fit the 68% most central yield of the source distribution
    double lowerlimit;
    double upperlimit;
    GetCentralInterval(*hSource, 0.9, lowerlimit, upperlimit, true);
    unsigned lowerbin = hSource->FindBin(lowerlimit);
    unsigned upperbin = hSource->FindBin(upperlimit);
    for(unsigned uBin=lowerbin; uBin<=upperbin; uBin++){
      hSource->SetBinError(uBin+1,0.01);
    }

    printf("Core size of %.3f fm\n",CoreSize);
    printf("The fit will be performed in the range [%.2f, %.2f] fm\n",lowerlimit,upperlimit);
    //fyi, GaussSourceTF1 is in DLM_Source.h if you want to check it out.
    TF1* fSource = new TF1("fSource",GaussSourceTF1,lowerlimit,upperlimit,1);
    fSource->SetParameter(0,CoreSize);
    fSource->SetParLimits(0,CoreSize*0.5,CoreSize*10.0);
    hSource->Fit(fSource,"S, N, R, M");
    printf("The effective Gaussian size is %.3f +/- %.3f fm\n",fSource->GetParameter(0),fSource->GetParError(0));

    p_Omega_AngleRcP2->Scale(1./p_Omega_AngleRcP2->Integral(),"width");
    pReso_Omega_AngleRcP1->Scale(1./pReso_Omega_AngleRcP1->Integral(),"width");
    pReso_dReso_AngleRcP1->Scale(1./pReso_dReso_AngleRcP1->Integral(),"width");
    pReso_dReso_AngleRcP2->Scale(1./pReso_dReso_AngleRcP2->Integral(),"width");
    pReso_dReso_AngleP1P2->Scale(1./pReso_dReso_AngleP1P2->Integral(),"width");

    //get rid of weird plotting
    for(unsigned uBin=0; uBin<NumSourceBins; uBin++){
      hSource->SetBinError(uBin+1,0.01);
    }
    hSource->Write();
    fSource->Write();
    p_Omega_AngleRcP2->Write();
    pReso_Omega_AngleRcP1->Write();
    pReso_dReso_AngleRcP1->Write();
    pReso_dReso_AngleRcP2->Write();
    pReso_dReso_AngleP1P2->Write();

    delete hSource;
    delete fSource;
    delete p_Omega_AngleRcP2;
    delete pReso_Omega_AngleRcP1;
    delete pReso_dReso_AngleRcP1;
    delete pReso_dReso_AngleRcP2;
    delete pReso_dReso_AngleP1P2;
    delete fOutput;

}



void pp_C_vs_SI(){

  const unsigned NumMomBins = 80;
  const double kMin = 0;
  const double kMax = 320;

  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());
  CATSparameters cPar(CATSparameters::tSource,1,true);

  CATS Cat;
  Cat.SetMomBins(NumMomBins,kMin,kMax);
  AnalysisObject.SetUpCats_pp(Cat,"AV18","Gauss",0,0);
  Cat.SetAnaSource(0,1.2);

  Cat.KillTheCat();
  TH1F* CkFull = new TH1F("CkFull","CkFull",NumMomBins,kMin,kMax);
  CkFull->SetLineColor(kBlack);
  CkFull->SetLineWidth(6);
  for(unsigned uBin=0; uBin<Cat.GetNumMomBins(); uBin++){
    CkFull->SetBinContent(uBin+1,Cat.GetCorrFun(uBin));
  }

  Cat.SetQ1Q2(0);
  Cat.KillTheCat();
  TH1F* CkSI = new TH1F("CkSI","CkSI",NumMomBins,kMin,kMax);
  CkSI->SetLineColor(kRed);
  CkSI->SetLineWidth(5);
  for(unsigned uBin=0; uBin<Cat.GetNumMomBins(); uBin++){
    CkSI->SetBinContent(uBin+1,Cat.GetCorrFun(uBin));
  }

  Cat.SetQ1Q2(1);
  Cat.RemoveShortRangePotential();
  Cat.KillTheCat();
  TH1F* CkC = new TH1F("CkC","CkC",NumMomBins,kMin,kMax);
  CkC->SetLineColor(kBlue);
  CkC->SetLineWidth(4);
  for(unsigned uBin=0; uBin<Cat.GetNumMomBins(); uBin++){
    CkC->SetBinContent(uBin+1,Cat.GetCorrFun(uBin));
  }

  TH1F* Cumulant = (TH1F*)CkFull->Clone("Cumulant");
  TH1F* CumulantRed = (TH1F*)CkFull->Clone("CumulantRed");
  CumulantRed->SetLineColor(kGreen+1);
  CumulantRed->SetLineWidth(6);
  TH1F* CumulantExp = (TH1F*)CkFull->Clone("CumulantExp");
  TH1F* CkExp = (TH1F*)CkFull->Clone("CkExp");
  CumulantRed->SetLineColor(kGreen+1);
  CumulantRed->SetLineWidth(4);
  for(unsigned uBin=0; uBin<Cat.GetNumMomBins(); uBin++){
    Cumulant->SetBinContent(uBin+1,Cumulant->GetBinContent(uBin+1)-1);
    CkExp->SetBinContent(uBin+1,CkSI->GetBinContent(uBin+1)*CkC->GetBinContent(uBin+1));
    CumulantExp->SetBinContent(uBin+1,CkSI->GetBinContent(uBin+1)*CkC->GetBinContent(uBin+1)-1);
    CumulantRed->SetBinContent(uBin+1,Cumulant->GetBinContent(uBin+1)-CumulantExp->GetBinContent(uBin+1));
  }

  TFile fOutput(TString::Format("%s/OtherTasks/pp_C_vs_SI.root",GetFemtoOutputFolder()),"recreate");
  CkFull->Write();
  CkC->Write();
  CkSI->Write();
  CkExp->Write();
  Cumulant->Write();
  CumulantExp->Write();
  CumulantRed->Write();

}


void p_pi_CommonTest1(){
  const unsigned NumBins = 200;
  const double kMin = 0;
  const double kMax = 400;
  bool RSM = true;
  const int source_settings = 201;
  const double SourceSize = 1.2;
  CATS Cat;
  Cat.SetMomBins(NumBins,kMin,kMax);

  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());
  AnalysisObject.SetUpCats_ppic(Cat,"DG_pip_d",RSM?"McGauss_ResoTM":"Gauss",0,RSM?source_settings:0);
  Cat.SetAnaSource(0,SourceSize);
  //if(RSM) Cat.SetAnaSource(1,2.0);
  Cat.KillTheCat();

  TH1F* hCk_ppi = new TH1F("hCk_ppi","hCk_ppi",NumBins,kMin,kMax);
  for(unsigned uBin=0; uBin<NumBins; uBin++){
    hCk_ppi->SetBinContent(uBin+1,Cat.GetCorrFun(uBin));
  }

  TH1F* hSource = new TH1F("hSource","hSource",1024,0,32);
  for(unsigned uBin=0; uBin<hSource->GetNbinsX(); uBin++){
    double rstar = hSource->GetBinCenter(uBin+1);
    hSource->SetBinContent(uBin+1,Cat.EvaluateTheSource(0,rstar,0));
  }


  TFile fOutput(TString::Format("%s/OtherTasks/p_pi_CommonTest1_%s%.2f.root",GetFemtoOutputFolder(),RSM?"rcore":"reff",SourceSize), "RECREATE");
  hCk_ppi->Write();
  hSource->Write();

  delete hCk_ppi;
  delete hSource;
}


void L_pi_CommonSource_based_on_p_pic(){
  const unsigned NumBins = 200;
  const double kMin = 0;
  const double kMax = 400;
  bool RSM = true;
  const int source_settings = -200;
  //r_core = 1.163 +0.047 -0.047 fm
  const double SourceSize = 1.163;
  CATS Cat;
  Cat.SetMomBins(NumBins,kMin,kMax);

  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());
  AnalysisObject.SetUpCats_ppic(Cat,"",RSM?"McGauss_ResoTM":"Gauss",0,RSM?source_settings:0);
  Cat.SetAnaSource(0,SourceSize);
  //if(RSM) Cat.SetAnaSource(1,2.0);
  Cat.KillTheCat();

  TH1F* hCk_Lpi = new TH1F("hCk_Lpi","hCk_Lpi",NumBins,kMin,kMax);
  for(unsigned uBin=0; uBin<NumBins; uBin++){
    hCk_Lpi->SetBinContent(uBin+1,Cat.GetCorrFun(uBin));
  }

  TH1F* hSource = new TH1F("hSource","hSource",1024,0,32);
  for(unsigned uBin=0; uBin<hSource->GetNbinsX(); uBin++){
    double rstar = hSource->GetBinCenter(uBin+1);
    hSource->SetBinContent(uBin+1,Cat.EvaluateTheSource(0,rstar,0)*10000);
  }
  hSource->Sumw2();

  TF1* fSource;
  //double r_eff = Get_reff_TF1(hSource, fsource, 1,  0.9);
  //double r_eff = GetReff(&AnalysisObject.GetCleverMcLevyResoTM_ppic(), SourceSize);

  hSource->Scale(1. / hSource->Integral(), "width");

  double lowerlimit;
  double upperlimit;
  GetCentralInterval(*hSource, 0.9, lowerlimit, upperlimit, true);

  fSource = new TF1("fSource", "[0]*4.*TMath::Pi()*x*x*pow(4.*TMath::Pi()*[1]*[1],-1.5)*exp(-(x*x)/(4.*[1]*[1]))+1.-[0]", 0.5, 6);
  fSource->SetParameter(0, 0.8);
  fSource->SetParLimits(0, 0.6, 1.0);
  fSource->SetParameter(1, hSource->GetMean() / 2.3);
  fSource->SetParLimits(1, hSource->GetMean() / 10., hSource->GetMean() * 2.);

  hSource->Fit(fSource, "S, N, R, M");

  double r_eff = fSource->GetParameter(1);
  double lambda = fSource->GetParameter(0);

  printf("r_core / r_eff = %.3f / %.3f fm with lambda = %.3f\n", SourceSize, r_eff, lambda);


  TFile fOutput(TString::Format("%s/OtherTasks/L_pi_CommonSource_based_on_p_pic_%s%.2f.root",GetFemtoOutputFolder(),RSM?"rcore":"reff",SourceSize), "RECREATE");
  hCk_Lpi->Write();
  hSource->Write();
  fSource->Write();

  delete hCk_Lpi;
  delete hSource;
  delete fSource;
}



//0 1
//2 3 4 5
//6 7 8 9
double TwoDoubleGauss(double* Pars){
  double Pars1[6];
  double Pars2[6];

  for(unsigned u=0; u<2; u++){
    Pars1[u] = Pars[u];
    Pars2[u] = Pars[u];
  }
  for(unsigned u=2; u<6; u++){
    Pars1[u] = Pars[u];
    Pars2[u] = Pars[u+4];
  }
  return DoubleGaussSum(Pars1)+DoubleGaussSum(Pars2);
}

void ShowEffectOfCum(){

  const double kMin = 0;
  const double kMax = 320;
  const unsigned NumMomBins = 160;
  const double SourceSize = 1.2;
  std::vector<double> f0_val;
  f0_val.push_back(0.001);
  f0_val.push_back(0.002);
  f0_val.push_back(0.005);
  f0_val.push_back(0.01);
  f0_val.push_back(0.02);
  f0_val.push_back(0.05);
  TH1F* hCk1 = new TH1F(TString::Format("hCk1_f%.3f",2.0),TString::Format("hCk1_f%.3f",2.0),NumMomBins,kMin,kMax);
  TH1F** hCk2 = new TH1F* [f0_val.size()+1];
  TH1F** hCk12 = new TH1F* [f0_val.size()+1];
  TH1F** hCk12M = new TH1F* [f0_val.size()+1];

  hCk2[0] = new TH1F(TString::Format("hCk2_f%.3f",0.0),TString::Format("hCk2_f%.3f",0.0),NumMomBins,kMin,kMax);
  hCk12[0] = new TH1F(TString::Format("hCk12_f%.3f",0.0),TString::Format("hCk12_f%.3f",0.0),NumMomBins,kMin,kMax);
  hCk12M[0] = new TH1F(TString::Format("hCk12M_f%.3f",0.0),TString::Format("hCk12M_f%.3f",0.0),NumMomBins,kMin,kMax);
  for(unsigned uF=0; uF<f0_val.size(); uF++){
    hCk2[uF+1] = new TH1F(TString::Format("hCk2_f%.3f",f0_val.at(uF)),TString::Format("hCk2_f%.3f",f0_val.at(uF)),NumMomBins,kMin,kMax);
    hCk12[uF+1] = new TH1F(TString::Format("hCk12_f%.3f",f0_val.at(uF)),TString::Format("hCk12_f%.3f",f0_val.at(uF)),NumMomBins,kMin,kMax);
    hCk12M[uF+1] = new TH1F(TString::Format("hCk12M_f%.3f",f0_val.at(uF)),TString::Format("hCk12M_f%.3f",f0_val.at(uF)),NumMomBins,kMin,kMax);
  }

  CATSparameters cSource(CATSparameters::tSource,2,true);
  cSource.SetParameter(0,SourceSize);
  cSource.SetParameter(1,2);

  CATS Cat1;
  Cat1.SetMomBins(NumMomBins,kMin,kMax);
  Cat1.SetThetaDependentSource(false);
  Cat1.SetAnaSource(GaussSource, cSource);
  Cat1.SetUseAnalyticSource(true);
  Cat1.SetMomentumDependentSource(false);
  Cat1.SetExcludeFailedBins(false);
  Cat1.SetQ1Q2(0);
  Cat1.SetQuantumStatistics(false);
  Cat1.SetRedMass( (Mass_p*Mass_L)/(Mass_p+Mass_L) );
  Cat1.SetNumChannels(1);
  Cat1.SetNumPW(0,1);
  Cat1.SetSpin(0,0);
  Cat1.SetChannelWeight(0, 1.);
  Cat1.SetEpsilonConv(2e-8);
  Cat1.SetEpsilonProp(2e-8);
  CATSparameters pPars(CATSparameters::tPotential,8,true);
  pPars.SetParameter(0,-4.372201e+02);
  pPars.SetParameter(1,1.757978e-01);
  pPars.SetParameter(2,-1.146283e+03);
  pPars.SetParameter(3,2.594236e-01);
  Cat1.SetShortRangePotential(0,0,DoubleGaussSum,pPars);
  Cat1.KillTheCat();
  for(unsigned uBin=0; uBin<NumMomBins; uBin++){
    hCk1->SetBinContent(uBin+1,Cat1.GetCorrFun(uBin));
  }

  CATS Cat2;
  Cat2.SetMomBins(NumMomBins,kMin,kMax);
  Cat2.SetThetaDependentSource(false);
  Cat2.SetAnaSource(GaussSource, cSource);
  Cat2.SetUseAnalyticSource(true);
  Cat2.SetMomentumDependentSource(false);
  Cat2.SetExcludeFailedBins(false);
  Cat2.SetQ1Q2(0);
  Cat2.SetQuantumStatistics(false);
  Cat2.SetRedMass( (Mass_p*Mass_L)/(Mass_p+Mass_L) );
  Cat2.SetNumChannels(1);
  Cat2.SetNumPW(0,1);
  Cat2.SetSpin(0,0);
  Cat2.SetChannelWeight(0, 1.);
  Cat2.SetEpsilonConv(2e-8);
  Cat2.SetEpsilonProp(2e-8);
  Cat2.KillTheCat();

  CATS Cat12;
  Cat12.SetMomBins(NumMomBins,kMin,kMax);
  Cat12.SetThetaDependentSource(false);
  Cat12.SetAnaSource(GaussSource, cSource);
  Cat12.SetUseAnalyticSource(true);
  Cat12.SetMomentumDependentSource(false);
  Cat12.SetExcludeFailedBins(false);
  Cat12.SetQ1Q2(0);
  Cat12.SetQuantumStatistics(false);
  Cat12.SetRedMass( (Mass_p*Mass_L)/(Mass_p+Mass_L) );
  Cat12.SetNumChannels(1);
  Cat12.SetNumPW(0,1);
  Cat12.SetSpin(0,0);
  Cat12.SetChannelWeight(0, 1.);
  Cat12.SetEpsilonConv(2e-8);
  Cat12.SetEpsilonProp(2e-8);
  pPars.SetParameter(0,-4.372201e+02);
  pPars.SetParameter(1,1.757978e-01);
  pPars.SetParameter(2,-1.146283e+03);
  pPars.SetParameter(3,2.594236e-01);
  Cat12.SetShortRangePotential(0,0,DoubleGaussSum,pPars);
  Cat12.KillTheCat();

  for(unsigned uBin=0; uBin<NumMomBins; uBin++){
    hCk2[0]->SetBinContent(uBin+1,Cat2.GetCorrFun(uBin));
    hCk12[0]->SetBinContent(uBin+1,Cat12.GetCorrFun(uBin));
    hCk12M[0]->SetBinContent(uBin+1,Cat1.GetCorrFun(uBin)*Cat2.GetCorrFun(uBin));
  }


  for(unsigned uF=0; uF<f0_val.size(); uF++){
    pPars.SetParameter(4,-4.372201e+02);
    pPars.SetParameter(5,1.757978e-01);
    pPars.SetParameter(6,-1.146283e+03);
    pPars.SetParameter(7,2.594236e-01);

    pPars.SetParameter(0,pPars.GetParameter(4)*f0_val.at(uF));
    pPars.SetParameter(1,pPars.GetParameter(5));
    pPars.SetParameter(2,pPars.GetParameter(6)*f0_val.at(uF));
    pPars.SetParameter(3,pPars.GetParameter(7));

    Cat2.SetShortRangePotential(0,0,DoubleGaussSum,pPars);
    Cat2.KillTheCat();

    //pPars.SetParameter(4,0);
    //pPars.SetParameter(5,1);
    //pPars.SetParameter(6,0);
    //pPars.SetParameter(7,1);
    Cat12.SetShortRangePotential(0,0,TwoDoubleGauss,pPars);
    Cat12.KillTheCat();

    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
      hCk2[uF+1]->SetBinContent(uBin+1,Cat2.GetCorrFun(uBin));
      hCk12[uF+1]->SetBinContent(uBin+1,Cat12.GetCorrFun(uBin));
      hCk12M[uF+1]->SetBinContent(uBin+1,Cat1.GetCorrFun(uBin)*Cat2.GetCorrFun(uBin));
    }
  }
printf("OUTPUT\n");
  TFile fOutput(TString::Format("%s/OtherTasks/ShowEffectOfCum.root",GetFemtoOutputFolder()), "RECREATE");
  for(unsigned uF=0; uF<=f0_val.size(); uF++){
    hCk2[uF]->Write();
    hCk12[uF]->Write();
    hCk12M[uF]->Write();
    delete hCk2[uF];
    delete hCk12[uF];
    delete hCk12M[uF];
  }
  delete [] hCk2;
  delete [] hCk12;
  delete [] hCk12M;
}


void pp_QS_Tests(){
  const int NumMomBin = 100;
  const double kMin = 0;
  const double kMax = 200;
  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles/",GetCernBoxDimi()));
  CATS Kitty;
  Kitty.SetMomBins(NumMomBin,kMin,kMax);
  AnalysisObject.SetUpCats_pp(Kitty,"AV18","Gauss",0,0);
  Kitty.SetAnaSource(0,1.2);
  Kitty.SetEpsilonConv(1e-8);
  Kitty.SetEpsilonProp(1e-8);
  Kitty.KillTheCat();

  TGraph gCk_QS;
  gCk_QS.SetName("gCk_QS");
  gCk_QS.SetLineWidth(3);
  gCk_QS.SetLineColor(kBlack);
  TGraph gCk_CO;
  gCk_CO.SetName("gCk_CO");
  gCk_CO.SetLineWidth(3);
  gCk_CO.SetLineColor(kGreen+1);
  TGraph gCk_SI;
  gCk_SI.SetName("gCk_SI");

  TGraph gCk_SI_CO;
  gCk_SI_CO.SetName("gCk_SI_CO");
  gCk_SI_CO.SetLineWidth(3);
  gCk_SI_CO.SetLineColor(kOrange+1);


  TGraph gCk_QS_CO;
  gCk_QS_CO.SetName("gCk_QS_CO");
  gCk_QS_CO.SetLineWidth(3);
  gCk_QS_CO.SetLineColor(kRed);


  TGraph gCk_Full;
  gCk_Full.SetName("gCk_Full");
  gCk_Full.SetLineWidth(3);
  gCk_Full.SetLineColor(kBlue+1);



  for(unsigned uBin=0; uBin<NumMomBin; uBin++){
    gCk_Full.SetPoint(uBin,Kitty.GetMomentum(uBin),Kitty.GetCorrFun(uBin));
  }

  Kitty.SetQuantumStatistics(false);
  Kitty.KillTheCat();
  for(unsigned uBin=0; uBin<NumMomBin; uBin++){
    gCk_SI_CO.SetPoint(uBin,Kitty.GetMomentum(uBin),Kitty.GetCorrFun(uBin));
  }


  Kitty.RemoveShortRangePotential(0,0);
  Kitty.RemoveShortRangePotential(0,2);
  Kitty.RemoveShortRangePotential(1,1);
  Kitty.RemoveShortRangePotential(2,1);
  Kitty.RemoveShortRangePotential(3,1);
  Kitty.SetQuantumStatistics(true);
  Kitty.KillTheCat();
  for(unsigned uBin=0; uBin<NumMomBin; uBin++){
    gCk_QS_CO.SetPoint(uBin,Kitty.GetMomentum(uBin),Kitty.GetCorrFun(uBin));
  }

  Kitty.SetQuantumStatistics(false);
  Kitty.KillTheCat();
  for(unsigned uBin=0; uBin<NumMomBin; uBin++){
    gCk_CO.SetPoint(uBin,Kitty.GetMomentum(uBin),Kitty.GetCorrFun(uBin));
  }

  Kitty.SetQuantumStatistics(true);
  Kitty.SetQ1Q2(0);
  Kitty.KillTheCat();
  for(unsigned uBin=0; uBin<NumMomBin; uBin++){
    gCk_QS.SetPoint(uBin,Kitty.GetMomentum(uBin),Kitty.GetCorrFun(uBin));
  }

  TFile fOutput(TString::Format("%s/OtherTasks/pp_QS_Tests.root",GetFemtoOutputFolder()), "RECREATE");
  gCk_Full.Write();
  gCk_SI_CO.Write();
  gCk_QS_CO.Write();
  gCk_CO.Write();
  gCk_QS.Write();


}

void pn_potential(){
  unsigned NumRadBins = 1024;
  double rMin = 0;
  double rMax = 5;
  TH1F* hPot_3S1 = new TH1F("hPot_3S1","hPot_3S1",NumRadBins,rMin,rMax);

  for(unsigned uRad=0; uRad<NumRadBins; uRad++){
    double rad = hPot_3S1->GetBinCenter(uRad+1);
    double Vrad = fDlmPot(NN_AV18,v18_Coupled3P2,0,1,-1,1,0,1,&rad,0);
    hPot_3S1->SetBinContent(uRad+1,Vrad);
  }
  TFile fOutput(TString::Format("%s/OtherTasks/pn_potential.root",GetFemtoOutputFolder()), "RECREATE");
  hPot_3S1->Write();
  delete hPot_3S1;
}

//compare AV18 vs Hulthen HulthenSmooth (3S1)
void pn_Ck(const double SourceSize=1.3){
  const unsigned NumMomBins = 75;
  const double kMin = 0;
  const double kMax = 300;

  const unsigned NumPots = 3;
  TGraph* gPot = new TGraph[NumPots];
  TGraph* gCk = new TGraph[NumPots];

  CATS* Kitty = new CATS[NumPots];
  TString* PotName = new TString[NumPots];
  CATSparameters sPars(CATSparameters::tSource,1,true);
  CATSparameters pPars(CATSparameters::tPotential,8,true);
  sPars.SetParameter(0,SourceSize);
  for(unsigned uPot=0; uPot<NumPots; uPot++){
    Kitty[uPot].SetMomBins(NumMomBins,kMin,kMax);
    Kitty[uPot].SetAnaSource(GaussSource, sPars);
    Kitty[uPot].SetUseAnalyticSource(true);
    Kitty[uPot].SetAutoNormSource(false);
    Kitty[uPot].SetQ1Q2(0);
    Kitty[uPot].SetQuantumStatistics(false);
    Kitty[uPot].SetRedMass( (Mass_p*Mass_n)/(Mass_p+Mass_n) );
    Kitty[uPot].SetNumChannels(1);
    Kitty[uPot].SetNumPW(0,1);
    Kitty[uPot].SetSpin(0,1);
    Kitty[uPot].SetChannelWeight(0, 1.);
    //Kitty[uPot].SetEpsilonConv(1e-8);
    //Kitty[uPot].SetEpsilonProp(1e-8);
    //Kitty[uPot].SetMaxRad(96);
    //Kitty[uPot].SetMaxRho(32);
    switch (uPot) {
      case 0:
        PotName[uPot] = "AV18";
        pPars.SetParameter(0,NN_AV18);
        pPars.SetParameter(1,v18_Coupled3P2);
        pPars.SetParameter(2,0);
        pPars.SetParameter(3,1);
        pPars.SetParameter(4,-1);
        pPars.SetParameter(5,1);
        pPars.SetParameter(6,0);
        pPars.SetParameter(7,1);
        Kitty[uPot].SetShortRangePotential(0,0,fDlmPot,pPars);
        break;
      case 1:
        PotName[uPot] = "Hulthen";
        pPars.SetParameter(0,32.5);
        pPars.SetParameter(1,226);
        Kitty[uPot].SetShortRangePotential(0,0,Hulthen,pPars);
        break;
      case 2:
        PotName[uPot] = "Hulthen alla Dimi";
        pPars.SetParameter(0,32.5);
        pPars.SetParameter(1,226);
        pPars.SetParameter(2,0.001);
        Kitty[uPot].SetShortRangePotential(0,0,HulthenSmooth,pPars);
        break;
      default:
        PotName[uPot] = "Unknown";
        break;
    }
    if(PotName[uPot]=="Unknown"){
      printf("WARNING: Unknown potential ID. Skipping it!\n");
      continue;
    }
    gPot[uPot].SetName(TString::Format("gPot_%s",PotName[uPot].Data()));
    gPot[uPot].SetLineColor(uPot+2);
    gPot[uPot].SetLineWidth(6-double(uPot)*0.75);
    gCk[uPot].SetName(TString::Format("gCk_%s",PotName[uPot].Data()));
    gCk[uPot].SetLineColor(uPot+2);
    gCk[uPot].SetLineWidth(6-double(uPot)*0.75);
    Kitty[uPot].KillTheCat();

    unsigned Counter=0;
    for(double Rad=0.025; Rad<32; Rad+=0.05){
      double PotVal = Kitty[uPot].EvaluateThePotential(0,0,1,Rad);
      gPot[uPot].SetPoint(Counter,Rad,PotVal);
      Counter++;
    }
    for(double uMom=0; uMom<NumMomBins; uMom++){
      gCk[uPot].SetPoint(uMom,Kitty[uPot].GetMomentum(uMom),Kitty[uPot].GetCorrFun(uMom));
    }
  }
  TFile fOutput(TString::Format("%s/OtherTasks/pn_Ck_r%.2f.root",GetFemtoOutputFolder(),SourceSize), "RECREATE");
  for(unsigned uPot=0; uPot<NumPots; uPot++){
    gPot[uPot].Write();
    gCk[uPot].Write();
  }
  delete [] Kitty;
  delete [] PotName;
  delete [] gPot;
  delete [] gCk;
}

void PlotAv18(){

  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());

  CATS AB_pp;
  AB_pp.SetMomBins(50,0,200);
  AnalysisObject.SetUpCats_pp(AB_pp,"AV18","Gauss",0,0);//McLevyNolan_Reso
  AB_pp.SetAnaSource(0,1.2);
  AB_pp.KillTheCat();

  TGraph gAV18_1S0;
  gAV18_1S0.SetName("gAV18_1S0");
  const unsigned NumRad = 1024;
  const double rMin = 0;
  const double rMax = 4;
  const double kStep = (rMax-rMin)/double(NumRad);

  for(unsigned uRad=0; uRad<NumRad; uRad++){
    double rad = kStep*0.5 + double(uRad)*kStep;
    double pot = AB_pp.EvaluateThePotential(0,0,10,rad);
    gAV18_1S0.SetPoint(uRad,rad,pot);
  }
  TFile fOutput(TString::Format("%s/OtherTasks/PlotAv18.root",GetFemtoOutputFolder()), "RECREATE");
  gAV18_1S0.Write();

}

void pp_at_different_radii(TString SourceDescription){
  std::vector<double> radii;
  radii.push_back(0.2);
  radii.push_back(0.4);
  radii.push_back(0.6);
  radii.push_back(0.8);
  radii.push_back(1.0);
  radii.push_back(1.2);
  radii.push_back(1.4);
  radii.push_back(1.6);

  unsigned NumMomBins = 75;
  double kMin = 0;
  double kMax = 300;

  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());

  std::vector<TGraph*> gCk;

  for(unsigned uRad=0; uRad<radii.size(); uRad++){
    printf("uRad = %u/%u\n",uRad+1,unsigned(radii.size()));
    double SourceSize = radii.at(uRad);
    //TString SourceDescription = "Gauss";
    //TString SourceDescription = "McGauss_ResoTM";
    //TString SourceDescription = "McLevy_ResoTM";

    CATS AB_pp;
    //DLM_Ck* Ck_pp;
    AB_pp.SetMomBins(NumMomBins,kMin,kMax);
    AB_pp.SetNotifications(CATS::nWarning);
    AnalysisObject.SetUpCats_pp(AB_pp,"AV18",SourceDescription,0,0);//McLevyNolan_Reso
    AB_pp.SetAnaSource(0,SourceSize);
    if(SourceDescription.Contains("Mc")){
        AB_pp.SetAnaSource(0,SourceSize);
        AB_pp.SetAnaSource(1,2.0);
    }
    AB_pp.SetEpsilonConv(2e-8);
    AB_pp.SetEpsilonProp(2e-8);
    AB_pp.KillTheCat();

    gCk.push_back(new TGraph());
    gCk.back()->SetName(TString::Format("gCk_%.2f",SourceSize));

    for(unsigned uMom=0; uMom<NumMomBins; uMom++){
      gCk.back()->SetPoint(uMom,AB_pp.GetMomentum(uMom),AB_pp.GetCorrFun(uMom));
    }
  }

  TFile fOutput(TString::Format("%s/OtherTasks/pp_at_different_radii_%s.root",GetFemtoOutputFolder(),SourceDescription.Data()), "RECREATE");
  for(unsigned uRad=0; uRad<radii.size(); uRad++){
    gCk.at(uRad)->Write();
  }
}

void TestMergeSort(){
  unsigned NumElements = 1*1000*1000;
  float* Element = new float[NumElements];
  TRandom3 rangen(11);
  for(unsigned uEl=0; uEl<NumElements; uEl++){
    Element[uEl] = rangen.Uniform(0,1000);
  }

  DLM_Sort < float, unsigned > SortTool;
  SortTool.SetData(Element,NumElements);
  SortTool.MergeSort();

  unsigned* Key1 = new unsigned[NumElements];
  SortTool.CopyKey(Key1);
  float* SortedElement = new float[NumElements];
  SortTool.GetSortedData(Element,SortedElement);

  for(unsigned uEl=0; uEl<NumElements; uEl++){
    //printf("%u @ %u: %.3f (%f)\n",uEl,Key1[uEl],Element[Key1[uEl]],SortedElement[uEl]);
    //usleep(100e3);
  }

  delete [] Element;
  delete [] Key1;
  delete [] SortedElement;
}





void LamLamStudy1(){

    //std::vector<float> SrcSize = {0.8,1.0,1.2,1.4,1.6,1.8,2.0};
    std::vector<float> SrcSize = {1.0,1.2,1.3,1.4,1.6};
    const double lambda_par = 0.34;
    const unsigned NumRad = SrcSize.size();

    const unsigned NumScenarios = 6;
    TString ScenarioName[NumScenarios];
    ScenarioName[0] = "ND46";
    ScenarioName[1] = "NF44";
    ScenarioName[2] = "Ehime";
    ScenarioName[3] = "ESC08";
    ScenarioName[4] = "HKMYY";
    ScenarioName[5] = "Without strong interaction";

    TString ScenarioShortName[NumScenarios];
    ScenarioShortName[0] = "ND46";
    ScenarioShortName[1] = "NF44";
    ScenarioShortName[2] = "Ehime";
    ScenarioShortName[3] = "ESC08";
    ScenarioShortName[4] = "HKMYY";
    ScenarioShortName[5] = "NoFSI";

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




    TString OutFolder = TString::Format("%s/OtherTasks/LamLamStudy1/",GetFemtoOutputFolder());
    printf("OutFolder = %s\n",OutFolder.Data());
    TString OutFile = OutFolder+"femto_LL.root";
    const double kMin = 0;
    const double kMax = 400;
    const unsigned NumMomBins = 100;
    CATS Kitty;
    Kitty.SetMomBins(NumMomBins,kMin,kMax);

    CATSparameters cPars (CATSparameters::tSource,1,true);
    cPars.SetParameter(0,1.2);

    Kitty.SetAnaSource(GaussSource, cPars);
    Kitty.SetUseAnalyticSource(true);
    Kitty.SetMomentumDependentSource(false);
    Kitty.SetExcludeFailedBins(false);
    Kitty.SetQ1Q2(0);
    Kitty.SetQuantumStatistics(true);
    Kitty.SetRedMass( Mass_L*0.5 );

    Kitty.SetNumChannels(2);
    Kitty.SetNumPW(0,1);
    Kitty.SetNumPW(1,0);

    Kitty.SetSpin(0,0);
    Kitty.SetSpin(1,1);

    Kitty.SetChannelWeight(0, 1./4.);
    Kitty.SetChannelWeight(1, 3./4.);


    CATSparameters pPars(CATSparameters::tPotential,6,true);
    //pPars.SetParameter(0,PotPars[uSce][2]);
    //pPars.SetParameter(1,PotPars[uSce][3]);
    //pPars.SetParameter(2,PotPars[uSce][4]);
    //pPars.SetParameter(3,PotPars[uSce][5]);
    //pPars.SetParameter(4,PotPars[uSce][6]);
    //pPars.SetParameter(5,PotPars[uSce][7]);
    Kitty.SetShortRangePotential(0,0,TripleGaussSum,pPars);
    Kitty.SetNotifications(CATS::nWarning);
    Kitty.SetEpsilonConv(1e-8);
    Kitty.SetEpsilonProp(1e-8);

    TGraph** Ck_LL = new TGraph* [NumScenarios];
    TGraph** Ck_LL_exp = new TGraph* [NumScenarios];
    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
      Ck_LL[uSce] = new TGraph [NumRad];
      Ck_LL_exp[uSce] = new TGraph [NumRad];
      for(unsigned uRad=0; uRad<NumRad; uRad++){
        Ck_LL[uSce][uRad].SetName(TString::Format("Ck_LL_%s_r%.2f",ScenarioShortName[uSce].Data(),SrcSize.at(uRad)));
        Ck_LL_exp[uSce][uRad].SetName(TString::Format("CkExp_LL_%s_r%.2f",ScenarioShortName[uSce].Data(),SrcSize.at(uRad)));
      }
    }

    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
      Kitty.SetShortRangePotential(0,0,0,PotPars[uSce][2]);
      Kitty.SetShortRangePotential(0,0,1,PotPars[uSce][3]);
      Kitty.SetShortRangePotential(0,0,2,PotPars[uSce][4]);
      Kitty.SetShortRangePotential(0,0,3,PotPars[uSce][5]);
      Kitty.SetShortRangePotential(0,0,4,PotPars[uSce][6]);
      Kitty.SetShortRangePotential(0,0,5,PotPars[uSce][7]);

      for(unsigned uRad=0; uRad<NumRad; uRad++){
        printf("%s at %.2f fm\n",ScenarioShortName[uSce].Data(),SrcSize.at(uRad));
        Kitty.SetAnaSource(0, SrcSize.at(uRad));
        Kitty.KillTheCat();

        for(unsigned uMom=0; uMom<NumMomBins; uMom++){
          Ck_LL[uSce][uRad].SetPoint(uMom,Kitty.GetMomentum(uMom),Kitty.GetCorrFun(uMom));
          Ck_LL_exp[uSce][uRad].SetPoint(uMom,Kitty.GetMomentum(uMom),Kitty.GetCorrFun(uMom)*lambda_par + 1.-lambda_par);
        }
      }
    }

    TFile fOutput(OutFile,"recreate");
    for(unsigned uRad=0; uRad<NumRad; uRad++){
      for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        Ck_LL[uSce][uRad].Write();
      }
    }

    for(unsigned uRad=0; uRad<NumRad; uRad++){
      for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        Ck_LL_exp[uSce][uRad].Write();
      }
    }

    for(unsigned uSce=0; uSce<NumScenarios; uSce++){
        delete [] PotPars[uSce];
        delete [] ScatPars[uSce];
        delete [] ScatParErr[uSce];
        delete [] ScatParsMorita[uSce];

        delete [] Ck_LL[uSce];
        delete [] Ck_LL_exp[uSce];
    }
    delete [] PotPars;
    delete [] ScatPars;
    delete [] ScatParErr;
    delete [] ScatParsMorita;

    delete [] Ck_LL;
    delete [] Ck_LL_exp;
}



void pXi_BUG_TEST(){
  CATS Cat_pXim_correct;
  CATS Cat_pXim_BUG;

  const unsigned NumMomBins = 60;
  const double kMax = 240;

  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());

  Cat_pXim_correct.SetMomBins(NumMomBins,0,kMax);
  Cat_pXim_BUG.SetMomBins(NumMomBins,0,kMax);

  AnalysisObject.SetUpCats_pXim(Cat_pXim_correct,"pXim_HALQCDPaper2020","Gauss");
  AnalysisObject.SetUpCats_pXim(Cat_pXim_BUG,"pXim_HALQCDPaper2020","Gauss");

  Cat_pXim_correct.SetRedMass((Mass_p * Mass_Xim) / (Mass_p + Mass_Xim));
  Cat_pXim_BUG.SetRedMass((Mass_p * Mass_Xim1530) / (Mass_p + Mass_Xim1530));

  Cat_pXim_correct.KillTheCat();
  Cat_pXim_BUG.KillTheCat();

  TGraph g_pXim_correct;
  g_pXim_correct.SetName("g_pXim_correct");

  TGraph g_pXim_BUG;
  g_pXim_BUG.SetName("g_pXim_BUG");

  for(unsigned uBin=0; uBin<NumMomBins; uBin++){
    double MOM = Cat_pXim_correct.GetMomentum(uBin);
    g_pXim_correct.SetPoint(uBin,MOM,Cat_pXim_correct.GetCorrFun(uBin));
    g_pXim_BUG.SetPoint(uBin,MOM,Cat_pXim_BUG.GetCorrFun(uBin));
  }

  TFile fOutput(TString::Format("%s/OtherTasks/pXi_BUG_TEST.root",GetFemtoOutputFolder()),"recreate");
  g_pXim_correct.Write();
  g_pXim_BUG.Write();

}


void coal_test_1(){

  const double reff_1stMt = 1.36;
  const unsigned NumBins = 5;
  double* BinLims = new double[NumBins+1];
  BinLims[0] = 0.4;
  BinLims[1] = 0.5;
  BinLims[2] = 0.6;
  BinLims[3] = 0.8;
  BinLims[4] = 1.0;
  BinLims[5] = 1.5;
  TH1F* hB2_tot = new TH1F("hB2_tot","hB2_tot",NumBins,BinLims);
  hB2_tot->SetBinContent(1, 0.0095);  hB2_tot->SetBinError(1, 0.001);
  hB2_tot->SetBinContent(2, 0.01);    hB2_tot->SetBinError(2, 0.001);
  hB2_tot->SetBinContent(3, 0.011);   hB2_tot->SetBinError(3, 0.001);
  hB2_tot->SetBinContent(4, 0.012);   hB2_tot->SetBinError(4, 0.001);
  hB2_tot->SetBinContent(5, 0.013);   hB2_tot->SetBinError(5, 0.001);

  TH1F* hB2_jet = new TH1F("hB2_jet","hB2_jet",NumBins,BinLims);
  hB2_jet->SetBinContent(1, 0.512);   hB2_jet->SetBinError(1, 0.19);
  hB2_jet->SetBinContent(2, 0.413);   hB2_jet->SetBinError(2, 0.33);
  hB2_jet->SetBinContent(3, 0.417);   hB2_jet->SetBinError(3, 0.06);
  hB2_jet->SetBinContent(4, 0.322);   hB2_jet->SetBinError(4, 0.03);
  hB2_jet->SetBinContent(5, 0.376);   hB2_jet->SetBinError(5, 0.02);

  //const double PowCoeff = 0.67;
  const double PowCoeff = 3;
  const double Nconv = hB2_tot->GetBinContent(2)*pow(reff_1stMt,PowCoeff);

  TFile fInput(TString::Format("%s/CatsFiles/Source/SourcePaper_Published.root",GetCernBoxDimi()),"read");
  TGraphErrors* g_pp = (TGraphErrors*)fInput.Get("g_reff_pp");
  //TGraphErrors* g_pL = (TGraphErrors*)fInput.Get("g_reff_pL");

  TGraph* gRcoal = new TGraph();
  gRcoal->SetName("gRcoal");

  TGraph* gRjet = new TGraph();
  gRjet->SetName("gRjet");

  for(unsigned uBin=0; uBin<NumBins; uBin++){
    double reff = pow(Nconv/hB2_tot->GetBinContent(uBin+1),1./PowCoeff);
    gRcoal->SetPoint(uBin,hB2_tot->GetBinCenter(uBin+1)*2,reff);

    reff = pow(Nconv/hB2_jet->GetBinContent(uBin+1),1./PowCoeff);
    gRjet->SetPoint(uBin,hB2_jet->GetBinCenter(uBin+1)*2,reff);
  }
  TFile fOutput(TString::Format("%s/OtherTasks/coal_test_1.root",GetFemtoOutputFolder()),"recreate");
  g_pp->Write();
  gRcoal->Write();
  gRjet->Write();

  delete [] BinLims;
}


void test_sqwell(){

  //std::vector<float> SrcSize = {0.7,1.0,1.4,2.0,2.8,4.0};
  std::vector<float> SrcSize = {1.0};

  const double kMin = 0;
  const double kMax = 200;
  const unsigned NumBins = 100;

  const double SourceSize = 1.0;

  const double V_1S0 = -16200./Mass_p*1.;//c.a. 17.27 MeV
  const double d_1S0 = 2.3039*1.;

  TGraph gPsSW;
  gPsSW.SetName(TString::Format("gPsSW"));
  gPsSW.Set(NumBins);

  TGraph gPsAnaSW;
  gPsAnaSW.SetName(TString::Format("gPsAnaSW"));
  gPsAnaSW.Set(NumBins);

  TGraph gPsAV18;
  gPsAV18.SetName(TString::Format("gPsAV18"));
  gPsAV18.Set(NumBins);

  CATS KittySW;
  KittySW.SetMomBins(NumBins,kMin,kMax);
  CATSparameters cPars (CATSparameters::tSource,1,true);
  cPars.SetParameter(0,SourceSize);
  KittySW.SetAnaSource(GaussSource, cPars);
  KittySW.SetUseAnalyticSource(true);
  KittySW.SetMomentumDependentSource(false);
  KittySW.SetExcludeFailedBins(false);
  KittySW.SetQ1Q2(1);
  KittySW.SetQuantumStatistics(true);
  //KittySW.SetRedMass( (938.*1116.)/(938.+1116.) );
  KittySW.SetRedMass( 0.5*Mass_p );
  KittySW.SetNumChannels(2);
  KittySW.SetNumPW(0,1);
  KittySW.SetNumPW(1,0);
  KittySW.SetSpin(0,0);
  KittySW.SetSpin(1,1);
  KittySW.SetChannelWeight(0, 1./4.);
  KittySW.SetChannelWeight(1, 3./4.);
  CATSparameters pPars1S0(CATSparameters::tPotential,2,true);
  pPars1S0.SetParameter(0,V_1S0);
  pPars1S0.SetParameter(1,d_1S0);
  KittySW.SetShortRangePotential(0,0,SquareWell,pPars1S0);
  KittySW.SetNotifications(CATS::nWarning);
  //KittySW.SetMaxRad(128);
  //KittySW.SetMaxRho(64);
  //KittySW.KillTheCat();

  CATS KittyAV18;
  KittyAV18.SetMomBins(NumBins,kMin,kMax);
  KittyAV18.SetAnaSource(GaussSource, cPars);
  KittyAV18.SetUseAnalyticSource(true);
  KittyAV18.SetMomentumDependentSource(false);
  KittyAV18.SetExcludeFailedBins(false);
  KittyAV18.SetQ1Q2(1);
  KittyAV18.SetQuantumStatistics(true);
  KittyAV18.SetRedMass( 0.5*Mass_p );
  KittyAV18.SetNumChannels(2);
  KittyAV18.SetNumPW(0,1);
  KittyAV18.SetNumPW(1,0);
  KittyAV18.SetSpin(0,0);
  KittyAV18.SetSpin(1,1);
  KittyAV18.SetChannelWeight(0, 1./4.);
  KittyAV18.SetChannelWeight(1, 3./4.);
  CATSparameters pPars1S0_AV18(CATSparameters::tPotential,8,true);
  //        double PotPars1S0[8] = {NN_AV18, v18_Coupled3P2, 1, 1, 1, 0, 0, 0};
  pPars1S0_AV18.SetParameter(0,NN_AV18);
  pPars1S0_AV18.SetParameter(1,v18_Coupled3P2);
  pPars1S0_AV18.SetParameter(2,1);
  pPars1S0_AV18.SetParameter(3,1);
  pPars1S0_AV18.SetParameter(4,1);
  pPars1S0_AV18.SetParameter(5,0);
  pPars1S0_AV18.SetParameter(6,0);
  pPars1S0_AV18.SetParameter(7,0);
  KittyAV18.SetShortRangePotential(0,0,fDlmPot,pPars1S0_AV18);
  KittyAV18.SetNotifications(CATS::nWarning);
  //KittyAV18.SetMaxRad(128);
  //KittyAV18.SetMaxRho(64);
  //KittyAV18.KillTheCat();

  TFile* OutputFile = new TFile(
	    TString::Format("%s/OtherTasks/test_sqwell/test1.root",GetFemtoOutputFolder()), "recreate");
	printf("File Created\n");

  for(float src_size : SrcSize){

    KittySW.SetAnaSource(0,src_size);
    KittySW.KillTheCat();
    KittyAV18.SetAnaSource(0,src_size);
    KittyAV18.KillTheCat();

    TGraph gKittySW;
    gKittySW.SetName(TString::Format("gKittySW_r%.1f",src_size));
    gKittySW.Set(NumBins);


    TGraph gKittyAV18;
    gKittyAV18.SetName(TString::Format("gKittyAV18_r%.1f",src_size));
    gKittyAV18.Set(NumBins);





    for (unsigned uBin = 0; uBin < NumBins; uBin++) {
      double kstar = KittySW.GetMomentum(uBin);
      //printf("C(%.2f) = %.2f\n", Kitty_SE.GetMomentum(uBin), Kitty_SE.GetCorrFun(uBin));
      gKittySW.SetPoint(uBin, kstar, KittySW.GetCorrFun(uBin));
      gKittyAV18.SetPoint(uBin, kstar, KittyAV18.GetCorrFun(uBin));
    }//mom
    gKittySW.Write();
    gKittyAV18.Write();


  }//src

  for (unsigned uBin = 0; uBin < NumBins; uBin++) {
    double kstar = KittySW.GetMomentum(uBin);
    gPsSW.SetPoint(uBin, kstar, KittySW.GetPhaseShift(uBin,0,0));
    gPsAV18.SetPoint(uBin, kstar, KittyAV18.GetPhaseShift(uBin,0,0));

    double redmass = KittySW.GetRedMass();
    double q1q2 = KittySW.GetQ1Q2();
    double kbar = sqrt(kstar*kstar-2.*redmass*V_1S0);
    double eta = redmass*q1q2*AlphaFS/kstar;
    double eta_bar = redmass*q1q2*AlphaFS/kbar;
    double F_L;
    double G_L;
    double F1_L;
    double G1_L;
    double F_L_bar;
    double G_L_bar;
    double F1_L_bar;
    double G1_L_bar;
    double f_l;
    double g_l;
    double f_l_bar;
    double g_l_bar;
    double phase_shift;
    double F_O;
    double G_O;

    double diff = 0.05*FmToNu;
    gsl_sf_coulomb_wave_FG_array(0,0,eta,kstar*d_1S0*FmToNu,&F_L,&G_L,&F_O,&G_O);
    gsl_sf_coulomb_wave_FG_array(0,0,eta,kstar*(d_1S0+diff)*FmToNu,&F1_L,&G1_L,&F_O,&G_O);
    gsl_sf_coulomb_wave_FG_array(0,0,eta_bar,kbar*d_1S0*FmToNu,&F_L_bar,&G_L_bar,&F_O,&G_O);
    gsl_sf_coulomb_wave_FG_array(0,0,eta_bar,kbar*(d_1S0+diff)*FmToNu,&F1_L_bar,&G1_L_bar,&F_O,&G_O);

  //G_L = fabs(G_L);
  //G_L_bar = fabs(G_L_bar);
  //G1_L = fabs(G1_L);
  //G1_L_bar = fabs(G1_L_bar);

    f_l = log(F1_L)-log(F_L);
  f_l = log(fabs(F1_L))-log(fabs(F_L));
    f_l /= diff;

    g_l = log(G1_L)-log(G_L);
  g_l = log(fabs(G1_L))-log(fabs(G_L));
    g_l /= diff;
  //g_l = 0;

    f_l_bar = log(F1_L_bar)-log(F_L_bar);
  f_l_bar = log(fabs(F1_L_bar))-log(fabs(F_L_bar));
    f_l_bar /= diff;

    g_l_bar = log(G1_L_bar)-log(G_L_bar);
  g_l_bar = log(fabs(G1_L_bar))-log(fabs(G_L_bar));
    g_l_bar /= diff;

    phase_shift = atan(F_L*(kstar*f_l-kbar*f_l_bar)/G_L/(kbar*f_l_bar-kstar*g_l));//as in slides

    if(kstar>120 && kstar<160){
      printf("kstar = %.2e\n",kstar);
      printf("kbar = %.2e\n",kbar);
      printf("F_L = %.2e\n",F_L);
      printf("F1_L = %.2e\n",F1_L);
      printf("F_L_bar = %.2e\n",F_L_bar);
      printf("G_L = %.2e\n",G_L);
      printf("G1_L = %.2e\n",G1_L);
      printf("G_L_bar = %.2e\n",G_L_bar);
      printf("f_l = %.2e\n",f_l);
      printf("f_l_bar = %.2e\n",f_l_bar);
      printf("g_l = %.2e\n",g_l);
      printf("g_l_bar = %.2e\n",g_l_bar);
      printf("phase_shift = %.2e\n",phase_shift);
    }

    gPsAnaSW.SetPoint(uBin, kstar, phase_shift);
  }


  gPsSW.Write();
  gPsAV18.Write();
  gPsAnaSW.Write();




  double f0_SW;
  double d0_SW;
  double f0_AV18;
  double d0_AV18;
  TH1F* hDummy;
  TF1* fDummy;
  KittySW.SetQ1Q2(0);
  KittySW.KillTheCat();
  KittyAV18.SetQ1Q2(0);
  KittyAV18.KillTheCat();
  Eval_ScattParameters(KittySW,f0_SW,d0_SW,hDummy,fDummy,2,false,false,0);
  delete hDummy; delete fDummy;
  Eval_ScattParameters(KittyAV18,f0_AV18,d0_AV18,hDummy,fDummy,2,false,false,0);
  delete hDummy; delete fDummy;
  printf("SW vs AV18:\n");
  printf(" f0: %.2f %.2f\n",f0_SW,f0_AV18);
  printf(" d0: %.2f %.2f\n",d0_SW,d0_AV18);

  TGraph gPsSW_0;
  gPsSW_0.SetName(TString::Format("gPsSW_0"));
  gPsSW_0.Set(NumBins);
  for (unsigned uBin = 0; uBin < NumBins; uBin++) {
    double kstar = KittySW.GetMomentum(uBin);
    gPsSW_0.SetPoint(uBin,kstar,KittySW.GetPhaseShift(uBin,0,0));
  }
  gPsSW_0.Write();

	delete OutputFile;
}


//par[0] is an overall normalization
//than we have a pol4 = p0*(1+p1*k+p2*k^2+p3*k^3+p4*k^4), which has 3 free arguments and the following properties
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
double DongFang_BL1(double* x, double* par){
    double& k = *x;
    double& p0 = par[0];

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

DLM_CkDecomposition* DongFang_FITTER;
double DongFang_Femto1(double* x, double* par){
  double& kstar = x[0];
  double& source_size = par[0];
  DongFang_FITTER->GetCk()->SetSourcePar(0,source_size);
  DongFang_FITTER->Update();
  return DongFang_FITTER->EvalCk(kstar);
}

//[0-4]: BL pars
//[5]: radius
double DongFang_Fit1(double* x, double* par){
  double Baseline = DongFang_BL1(x,par);
  double Femto = DongFang_Femto1(x,&par[5]);
  return Baseline*Femto;
}

void DongFang_Example1(){


  DLM_CommonAnaFunctions AnalysisObject;
  //change the path to the location where you have downloaded the relevant CernBox folders
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles/",GetCernBoxDimi()));
  TString SourceDescription = "Gauss";
  unsigned NumMomBins_pp = 70;
  double kMin_pp = 0;
  double kMax_pp = 350;
  double CkConv_pp = 700;//the point at which C(k) = 1
  double SourceSize = 5.0;
  const double FitRegion_min = 0;
  const double FitRegion_max = kMax_pp;

  //set up pp
  CATS Kitty_pp;
  Kitty_pp.SetMomBins(NumMomBins_pp,kMin_pp,kMax_pp);
  AnalysisObject.SetUpCats_pp(Kitty_pp,"AV18",SourceDescription,0,0);
  Kitty_pp.SetAnaSource(0,SourceSize);
  Kitty_pp.SetNotifications(CATS::nWarning);//avoid printing stuff
  //evaluate the theoretical correlation function (NO lambda pars, mom resolution etc.)
  Kitty_pp.KillTheCat();

  //something like a histogram
  DLM_Ck* Ck_pp = new DLM_Ck(Kitty_pp.GetNumSourcePars(),0,Kitty_pp,2.*NumMomBins_pp,0,kMax_pp*2);
  for(unsigned uSP=0; uSP<Kitty_pp.GetNumSourcePars(); uSP++){
    Ck_pp->SetSourcePar(uSP,Kitty_pp.GetAnaSourcePar(uSP));
  }
  //extrapolate C(k) with a linear function between C(kMax_pp) --> 1 in the k* range of kMax_pp --> CkConv_pp
  Ck_pp->SetCutOff(kMax_pp,CkConv_pp);

  //set up pLambda (needed for feed-down)
  CATS Kitty_pL;
  Kitty_pL.SetMomBins(NumMomBins_pp,kMin_pp,kMax_pp);
  //we will use the Usmani potential, due to the limited range in r* of the chiral wave functions
  AnalysisObject.SetUpCats_pL(Kitty_pL,"UsmaniFit",SourceDescription,0,0);
  //set the correct potential parameters foro the Usmani potential (tuned to NLO19)
  //https://arxiv.org/pdf/2305.08441.pdf
  Kitty_pL.SetShortRangePotential(1,0,1,2279.0);
  Kitty_pL.SetShortRangePotential(1,0,2,0.3394);
  Kitty_pL.SetShortRangePotential(1,0,3,0.2614);
  Kitty_pL.SetAnaSource(0,SourceSize);
  Kitty_pL.SetNotifications(CATS::nWarning);//avoid printing stuff
  //evaluate the theoretical correlation function (NO lambda pars, mom resolution etc.)
  Kitty_pL.KillTheCat();

  //something like a histogram
  DLM_Ck* Ck_pL = new DLM_Ck(Kitty_pL.GetNumSourcePars(),0,Kitty_pL,2.*NumMomBins_pp,0,kMax_pp*2);
  for(unsigned uSP=0; uSP<Kitty_pL.GetNumSourcePars(); uSP++){
    Ck_pL->SetSourcePar(uSP,Kitty_pL.GetAnaSourcePar(uSP));
  }
  //extrapolate C(k) with a linear function between C(kMax_pp) --> 1 in the k* range of kMax_pp --> CkConv_pp
  Ck_pL->SetCutOff(kMax_pp,CkConv_pp);

  //52.3% feed-down into pLambda, given 0.6% putity (https://inspirehep.net/literature/1857549)
  //since we will not define fake contribution, we need to rescale the feed-down by 100%-0.6%
  //N.B. this is now for pp 13 TeV, parhaps change in the future, no idea how different it is in PbPb
  const double lam_pL_flat = (0.523)/0.994;

  //figure it out from the template fits and purities
  //details in appendix of https://inspirehep.net/literature/1675759
  const double lam_pp_pL = 0.15;
  const double lam_pp_flat = 0.15;
  const double lam_pp_fake = 0.005;


  //use you own, but here you can compare to the pp 13 TeV
  //N.B.!!!!!!! if you use your own histogram, you have to have the axis in MEV !!!!!!!!!!!!!!!
  TH2F* hReso_pp = AnalysisObject.GetResolutionMatrix("pp13TeV_HM_DimiJun20","pp");
  //get the decay smear matrix
  TH2F* hFeed_pp_pL = AnalysisObject.GetResidualMatrix("pp","pLambda");

  DLM_CkDecomposition* CkDec_pL = new DLM_CkDecomposition("pLambda",1,*Ck_pL,NULL);
  CkDec_pL->AddContribution(0,lam_pL_flat,DLM_CkDecomposition::cFeedDown);

  //perhaps interesting to check differences in the hReso taken from pp 13 TeV HM (mixed events)
  //and the one that you will get in PbPb (same events)
  DLM_CkDecomposition* CkDec_pp = new DLM_CkDecomposition("pp",3,*Ck_pp,hReso_pp);;
  CkDec_pp->AddContribution(0,lam_pp_pL,DLM_CkDecomposition::cFeedDown,CkDec_pL,hFeed_pp_pL);
  CkDec_pp->AddContribution(1,lam_pp_flat,DLM_CkDecomposition::cFeedDown);
  CkDec_pp->AddContribution(2,lam_pp_fake,DLM_CkDecomposition::cFake);
  //this is your mixed event distribution (in MEV!!!), needed for the correct application of the smearing matrix
  //https://docs.google.com/presentation/d/1fu88IniMeMGpJxmhl7KAwjby7s3qIwE4hNtsdZwjaRk/edit?usp=sharing
  //CkDec_pp->AddPhaseSpace(hPs_pp);
  //CkDec_pp->AddPhaseSpace(0,hPs_pp);

  //load you data in MEV
  TH1F* hData;


  //[0] = norm constant
  //[1,2] = the positions of the extrema of the polynomial function
  //        typically we fix [1] to 0 and fit only [2] (=> smooth baseline at kstar=0)
  //[3] = the 3rd coeficient of the polynomial function
  //[4] = fix to zero to fit with a pol3 BL (otherwise pol4)
  //[5] = source size
  TF1* fit_pp = new TF1("fit_pp",DongFang_Fit1,FitRegion_min,FitRegion_max,6);
  //set up starting values for the parameters, limits etc...
  //fit your hData
  //.........

}

void TestNewCutOff(){
  const unsigned NumMomBins = 50;
  const double kMin = 0;
  const double kMax = 200;
  const double kNonZero = 100;
  DLM_Ck TestCk(NumMomBins,kMin,kMax);
  for(unsigned uBin=0; uBin<NumMomBins; uBin++){
    double Momentum = TestCk.GetBinCenter(0,uBin);
    if(Momentum<kNonZero){
      TestCk.SetBinContent(uBin,1.+exp(-pow(Momentum/50.,2.))-0.15*exp(-pow((Momentum-100)/30.,2.)));
    }
    else{
      TestCk.SetBinContent(uBin,0);
    }
  }
  TestCk.Update();

  TGraph gTestCk;
  gTestCk.SetName("gTestCk");

  TGraph gEvTestCk;
  gEvTestCk.SetName("gEvTestCk");

  for(unsigned uBin=0; uBin<NumMomBins; uBin++){
    double Momentum = TestCk.GetBinCenter(0,uBin);
    gTestCk.SetPoint(uBin,Momentum,TestCk.GetBinContent(uBin));
    gEvTestCk.SetPoint(uBin,Momentum,TestCk.Eval(Momentum));
  }

  TestCk.SetCutOff(kNonZero,kMax*2);
  TestCk.Update(true);
  TGraph gCutTestCk;
  gCutTestCk.SetName("gCutTestCk");
  TGraph gCutEvTestCk;
  gCutEvTestCk.SetName("gCutEvTestCk");
  for(unsigned uBin=0; uBin<NumMomBins; uBin++){
    double Momentum = TestCk.GetBinCenter(0,uBin);
    gCutTestCk.SetPoint(uBin,Momentum,TestCk.GetBinContent(uBin));
    gCutEvTestCk.SetPoint(uBin,Momentum,TestCk.Eval(Momentum));
  }



  TFile fOutput(TString::Format("%s/OtherTasks/TestNewCutOff1.root",GetFemtoOutputFolder()),"recreate");
  gTestCk.Write();
  gEvTestCk.Write();
  gCutTestCk.Write();
  gCutEvTestCk.Write();


}


//study the effect of the feed-down
void pSigmaPlus_pp(TString Potential, double SourzeSize){

  const int NumMomBins = 160;
  const double kMin = 0;
  const double kMax = 320;
  CATS Kitty_pp;
  Kitty_pp.SetMomBins(NumMomBins,kMin,kMax);

  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());
  AnalysisObject.SetUpCats_pp(Kitty_pp,"AV18","Gauss",1,0);
  Kitty_pp.SetAnaSource(0, SourzeSize);
  Kitty_pp.KillTheCat();

  CATSparameters cPars(CATSparameters::tSource,1,true);
  cPars.SetParameter(0,SourzeSize);

  CATS Kitty_pSpC;
  Kitty_pSpC.SetMomBins(NumMomBins,kMin,kMax);
  Kitty_pSpC.SetThetaDependentSource(false);
  Kitty_pSpC.SetAnaSource(GaussSource, cPars);
  Kitty_pSpC.SetAutoNormSource(false);
  Kitty_pSpC.SetUseAnalyticSource(true);
  Kitty_pSpC.SetMomentumDependentSource(false);

  Kitty_pSpC.SetExcludeFailedBins(false);
  Kitty_pSpC.SetQ1Q2(1);
  Kitty_pSpC.SetQuantumStatistics(false);
  Kitty_pSpC.SetRedMass( (Mass_p*Mass_Sch)/(Mass_p+Mass_Sch) );
  Kitty_pSpC.SetNumChannels(2);
  Kitty_pSpC.SetNumPW(0,1);
  Kitty_pSpC.SetNumPW(1,1);
  Kitty_pSpC.SetChannelWeight(0, 1./4.);
  Kitty_pSpC.SetChannelWeight(1, 3./4.);

  Kitty_pSpC.KillTheCat();


  CATS Kitty_pSpSI;
  Kitty_pSpSI.SetMomBins(NumMomBins,kMin,kMax);
  Kitty_pSpSI.SetThetaDependentSource(false);
  Kitty_pSpSI.SetAnaSource(GaussSource, cPars);
  Kitty_pSpSI.SetAutoNormSource(false);
  Kitty_pSpSI.SetUseAnalyticSource(true);
  Kitty_pSpSI.SetMomentumDependentSource(false);

  Kitty_pSpSI.SetExcludeFailedBins(false);
  Kitty_pSpSI.SetQ1Q2(1);
  Kitty_pSpSI.SetQuantumStatistics(false);
  Kitty_pSpSI.SetRedMass( (Mass_p*Mass_Sch)/(Mass_p+Mass_Sch) );
  Kitty_pSpSI.SetNumChannels(2);
  Kitty_pSpSI.SetNumPW(0,1);
  Kitty_pSpSI.SetNumPW(1,1);
  Kitty_pSpSI.SetChannelWeight(0, 1./4.);
  Kitty_pSpSI.SetChannelWeight(1, 3./4.);

  CATSparameters pPars1S0(CATSparameters::tPotential,4,true);
  CATSparameters pPars3S1(CATSparameters::tPotential,4,true);

  if(Potential=="N2LO"){
    //f_best     3.3900003896617
    //d_best     3.7267469996733706
    Kitty_pSpSI.SetShortRangePotential(0, 0, DoubleGaussSum, pPars1S0);
    Kitty_pSpSI.SetShortRangePotential(0, 0, 0, -79.07298063891983);
    Kitty_pSpSI.SetShortRangePotential(0, 0, 1, 1.3716237629795942);
    Kitty_pSpSI.SetShortRangePotential(0, 0, 2, 2287.907638142889);
    Kitty_pSpSI.SetShortRangePotential(0, 0, 3, 0.43371332814881547);

    //f_best     -0.47641737043176346
    //d_best     -5.495670225918186
    Kitty_pSpSI.SetShortRangePotential(1, 0, DoubleGaussSum, pPars3S1);
    Kitty_pSpSI.SetShortRangePotential(1, 0, 0, 6.6923386649252645);
    Kitty_pSpSI.SetShortRangePotential(1, 0, 1, 1.970250601599477);
    Kitty_pSpSI.SetShortRangePotential(1, 0, 2, -48.30709539991798);
    Kitty_pSpSI.SetShortRangePotential(1, 0, 3, 0.3549784459798298);
  }
  else if(Potential=="NLO19"){
    //    f_best     3.61947029128352
    //    d_best     3.503042374656197
    Kitty_pSpSI.SetShortRangePotential(0, 0, DoubleGaussSum, pPars1S0);
    Kitty_pSpSI.SetShortRangePotential(0, 0, 0, -82.20693060436437);
    Kitty_pSpSI.SetShortRangePotential(0, 0, 1, 1.3392323926534693);
    Kitty_pSpSI.SetShortRangePotential(0, 0, 2, 2274.508479434418);
    Kitty_pSpSI.SetShortRangePotential(0, 0, 3, 0.4182867242441114);

    //f_best     -0.46965864027889725
    //d_best     -5.778925131802048
    Kitty_pSpSI.SetShortRangePotential(1, 0, DoubleGaussSum, pPars3S1);
    Kitty_pSpSI.SetShortRangePotential(1, 0, 0, 6.295180245126864);
    Kitty_pSpSI.SetShortRangePotential(1, 0, 1, 1.995145226219902);
    Kitty_pSpSI.SetShortRangePotential(1, 0, 2, -58.68689895288922);
    Kitty_pSpSI.SetShortRangePotential(1, 0, 3, 0.32017712049353114);
  }

  Kitty_pSpSI.KillTheCat();

  DLM_Ck dlmCk_pp(Kitty_pp.GetNumSourcePars(),0,Kitty_pp);
  DLM_Ck dlmCk_pSpC(Kitty_pSpC.GetNumSourcePars(),0,Kitty_pSpC);
  DLM_Ck dlmCk_pSpSI(Kitty_pSpSI.GetNumSourcePars(),0,Kitty_pSpSI);

  TH2F* hResolution_pp = AnalysisObject.GetResolutionMatrix("pp13TeV_HM_DimiJun20","pp");
  TH2F* hFeed_pp_pSp = AnalysisObject.GetResidualMatrix("pp","pSigmaPlus");
  const double lam_gen = 0.67;
  const double lam_pSp = 0.085;
  const double lam_flat = 1.-lam_gen-lam_pSp;

  DLM_CkDecomposition decCk_pp_flt("pp",1,dlmCk_pp,hResolution_pp);
  DLM_CkDecomposition decCk_pp_pSpC("pp",2,dlmCk_pp,hResolution_pp);
  DLM_CkDecomposition decCk_pp_pSpSI("pp",2,dlmCk_pp,hResolution_pp);

  DLM_CkDecomposition decCk_pSpC("pSpSI",0,dlmCk_pSpC,NULL);
  DLM_CkDecomposition decCk_pSpSI("pSpSI",0,dlmCk_pSpSI,NULL);

  decCk_pp_flt.AddContribution(0,lam_pSp+lam_flat,DLM_CkDecomposition::cFeedDown);

  decCk_pp_pSpC.AddContribution(0,lam_pSp,DLM_CkDecomposition::cFeedDown,&decCk_pSpC,hFeed_pp_pSp);
  decCk_pp_pSpC.AddContribution(1,lam_flat,DLM_CkDecomposition::cFeedDown);

  decCk_pp_pSpSI.AddContribution(0,lam_pSp,DLM_CkDecomposition::cFeedDown,&decCk_pSpC,hFeed_pp_pSp);
  decCk_pp_pSpSI.AddContribution(1,lam_flat,DLM_CkDecomposition::cFeedDown);

  TH1F* hME_PP;
  TH1F* hME_APAP;
  TH1F* hPs_pp;
  TString FileName;
  //for(unsigned uMt=0; uMt<NumMtBins_pp; uMt++){
  unsigned uMt = 3;
  FileName = TString::Format("%s/CatsFiles/ExpData/Bernie_Source/ppData/mTBin_%i/CFOutput_mT_ppVar0_HM_%i.root",
                                      GetCernBoxDimi(),uMt+1,uMt);
  TFile* inFile = new TFile(FileName,"read");
  //PARTICLES
  TList* list1_tmp = (TList*)inFile->Get("PairDist");
  TList* list2_tmp = (TList*)list1_tmp->FindObject("PairFixShifted");//there is also PairShifted ?? which one ??
  hME_PP = (TH1F*)list2_tmp->FindObject(TString::Format("MEPart_mT_%i_FixShifted",uMt));
  list1_tmp = (TList*)inFile->Get("AntiPairDist");
  list2_tmp = (TList*)list1_tmp->FindObject("PairFixShifted");
  hME_APAP = (TH1F*)list2_tmp->FindObject(TString::Format("MEAntiPart_mT_%i_FixShifted",uMt));
  gROOT->cd();
  hPs_pp = (TH1F*)hME_PP->Clone("hPs_pp");
  hPs_pp->Add(hME_APAP);
  hPs_pp->GetXaxis()->SetLimits(hPs_pp->GetXaxis()->GetXmin()*1000.,hPs_pp->GetXaxis()->GetXmax()*1000.);
  delete inFile;
  //}

  decCk_pp_flt.AddPhaseSpace(hPs_pp);
  decCk_pp_pSpC.AddPhaseSpace(hPs_pp);
  decCk_pp_pSpSI.AddPhaseSpace(hPs_pp);

  decCk_pp_pSpC.AddPhaseSpace(0,hPs_pp);
  decCk_pp_pSpSI.AddPhaseSpace(0,hPs_pp);

  decCk_pp_flt.Update(true,true);
  decCk_pp_pSpC.Update(true,true);
  decCk_pp_pSpSI.Update(true,true);

  TGraph g_pSp_Coulomb;
  g_pSp_Coulomb.SetName("g_pSp_Coulomb");
  TGraph g_pSp_CoulombSI;
  g_pSp_CoulombSI.SetName("g_pSp_CoulombSI");

  TGraph g_pp_flt;
  g_pp_flt.SetName("g_pp_flt");
  TGraph g_pp_pSpC;
  g_pp_pSpC.SetName("g_pp_pSpC");
  TGraph g_pp_pSpSI;
  g_pp_pSpSI.SetName("g_pp_pSpSI");



  for(unsigned uMom=0; uMom<NumMomBins; uMom++){
    double MOM = Kitty_pSpC.GetMomentum(uMom);
    g_pSp_Coulomb.SetPoint(uMom,MOM,Kitty_pSpC.GetCorrFun(uMom));
    g_pSp_CoulombSI.SetPoint(uMom,MOM,Kitty_pSpSI.GetCorrFun(uMom));

    g_pp_flt.SetPoint(uMom,MOM,decCk_pp_flt.EvalCk(MOM));
    g_pp_pSpC.SetPoint(uMom,MOM,decCk_pp_pSpC.EvalCk(MOM));
    g_pp_pSpSI.SetPoint(uMom,MOM,decCk_pp_pSpSI.EvalCk(MOM));
  }
  TFile fOutput(TString::Format("%s/OtherTasks/pSigmaPlus_pp/fOutput_%s_r%.2f.root",GetFemtoOutputFolder(),Potential.Data(),SourzeSize),"recreate");
  g_pSp_Coulomb.Write();
  g_pSp_CoulombSI.Write();
  g_pp_flt.Write();
  g_pp_pSpC.Write();
  g_pp_pSpSI.Write();

  delete hPs_pp;
}


//compare this to the actual result from JH wave function
void pSigmaPlus_NLO19_data(){
  TGraph gNLO19_1p2;
  gNLO19_1p2.SetName("gNLO19_1p2");
  gNLO19_1p2.SetPoint(0, 9.773158484268876, 0.5042434679244419);
  gNLO19_1p2.SetPoint(1, 11.273497565473114, 0.6018135129574737);
  gNLO19_1p2.SetPoint(2, 12.457614849470772, 0.7003527968484495);
  gNLO19_1p2.SetPoint(3, 14.586587623917005, 0.8008153188506288);
  gNLO19_1p2.SetPoint(4, 17.031020215335587, 0.9012747929318712);
  gNLO19_1p2.SetPoint(5, 19.15999298978182, 1.0017373149340503);
  gNLO19_1p2.SetPoint(6, 23.49566052256602, 1.1041108833636852);
  gNLO19_1p2.SetPoint(7, 28.510252444051567, 1.1456083269199995);
  gNLO19_1p2.SetPoint(8, 31.03240701936189, 1.1475163254265182);
  gNLO19_1p2.SetPoint(9, 33.56065743654607, 1.1416947964369801);
  gNLO19_1p2.SetPoint(10, 40.21122092092927, 1.1087802982390635);
  gNLO19_1p2.SetPoint(11, 46.87778599023138, 1.0555757903639977);
  gNLO19_1p2.SetPoint(12, 53.229653222795385, 1.0014081394728618);
  gNLO19_1p2.SetPoint(13, 64.94357536365507, 0.9481548648628053);
  gNLO19_1p2.SetPoint(14, 74.12010332451973, 0.9123174104863718);
  gNLO19_1p2.SetPoint(15, 86.12738785555902, 0.8870806251285839);
  gNLO19_1p2.SetPoint(16, 102.85209201673308, 0.8801557487598769);
  gNLO19_1p2.SetPoint(17, 117.04473585955176, 0.8838833560657433);
  gNLO19_1p2.SetPoint(18, 134.7021038274267, 0.8943407728003533);
  gNLO19_1p2.SetPoint(19, 151.0961085669438, 0.9067427630927252);
  gNLO19_1p2.SetPoint(20, 167.8055731234332, 0.91914170546416);
  gNLO19_1p2.SetPoint(21, 188.301317463825, 0.9305378818473445);
  gNLO19_1p2.SetPoint(22, 211.0083284439602, 0.9380479590359423);
  gNLO19_1p2.SetPoint(23, 235.29263850895703, 0.9455427966198553);
  gNLO19_1p2.SetPoint(24, 262.1006271097328, 0.953013250836273);
  gNLO19_1p2.SetPoint(25, 294.2744805199754, 0.9565671266487343);

  TGraph gNLO19_2p5;
  gNLO19_2p5.SetName("gNLO19_2p5");
  gNLO19_2p5.SetPoint(0, 11.907379733799992, 0.500377077993976);
  gNLO19_2p5.SetPoint(1, 13.255414205077791, 0.5989005118040123);
  gNLO19_2p5.SetPoint(2, 16.195383698151794, 0.7032620730772506);
  gNLO19_2p5.SetPoint(3, 19.12915889153009, 0.7998188167339007);
  gNLO19_2p5.SetPoint(4, 24.924700544324104, 0.9022012992543611);
  gNLO19_2p5.SetPoint(5, 29.406276374166687, 0.9489868448560212);
  gNLO19_2p5.SetPoint(6, 37.3543371712182, 0.963543449140928);
  gNLO19_2p5.SetPoint(7, 50.04336009786994, 0.9517123367221314);
  gNLO19_2p5.SetPoint(8, 68.45746451827708, 0.9534839064351031);
  gNLO19_2p5.SetPoint(9, 121.80896779738447, 0.9763780381104289);
  gNLO19_2p5.SetPoint(10, 166.57749456063058, 0.9847217598005436);
  gNLO19_2p5.SetPoint(11, 210.70723416775712, 0.9881936647799863);
  gNLO19_2p5.SetPoint(12, 246.26406299602792, 0.9897979884011738);
  gNLO19_2p5.SetPoint(13, 297.6930105070808, 0.9902718523278953);

  const int NumMomBins = 160;
  const double kMin = 0;
  const double kMax = 320;
  const double SourzeSize = 1.2;


  CATS Kitty_pSp_Haide;
  Kitty_pSp_Haide.SetMomBins(NumMomBins,kMin,kMax);
  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());
  AnalysisObject.SetUpCats_pSp(Kitty_pSp_Haide,"N2LO","Gauss",1,0);
  Kitty_pSp_Haide.SetAnaSource(0, SourzeSize);
  Kitty_pSp_Haide.KillTheCat();

  CATS Kitty_pSp_DG;
  Kitty_pSp_DG.SetMomBins(NumMomBins,kMin,kMax);
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());
  AnalysisObject.SetUpCats_pSp(Kitty_pSp_DG,"DG_N2LO","Gauss",1,0);
  Kitty_pSp_DG.SetAnaSource(0, SourzeSize);
  Kitty_pSp_DG.KillTheCat();

  TGraph gHaide;
  gHaide.SetName("gHaide");
  TGraph gDG;
  gDG.SetName("gDG");
  for(unsigned uMom=0; uMom<NumMomBins; uMom++){
    gHaide.SetPoint(uMom,Kitty_pSp_Haide.GetMomentum(uMom),Kitty_pSp_Haide.GetCorrFun(uMom));
    gDG.SetPoint(uMom,Kitty_pSp_DG.GetMomentum(uMom),Kitty_pSp_DG.GetCorrFun(uMom));
  }


  TFile fOutput(TString::Format("%s/OtherTasks/pSigmaPlus_pp/fData.root",GetFemtoOutputFolder()),"recreate");
  gNLO19_1p2.Write();
  gNLO19_2p5.Write();
  gHaide.Write();
  gDG.Write();

}

void Test_CoulombEffRangeExp(){

  TApplication *myapp = new TApplication("myapp", 0, 0);

  const unsigned NumMomBins = 100;
  const double kMin = 0;
  const double kMax = 200;

  const double f_val = 4.0;
  const double d_val = 3.0;
  double q1q2rm = (1.)*(Mass_p*Mass_Sch)/(Mass_p+Mass_Sch);
  double PARS[5];
  PARS[0] = q1q2rm;
  PARS[1] = 1./f_val;
  PARS[2] = d_val;
  PARS[3] = 0;
  PARS[4] = 0;

  TH1F* hPS_cs = new TH1F("hPS_cs","hPS_cs",NumMomBins,kMin,kMax);
  TH1F* hPS_s = new TH1F("hPS_s","hPS_s",NumMomBins,kMin,kMax);

  for(unsigned uMom=0; uMom<NumMomBins; uMom++){
    double kstar = hPS_cs->GetBinCenter(uMom+1);

    PARS[0] = q1q2rm;
    double ps_val = atan(kstar/EffRangeExp(&kstar,PARS));
    hPS_cs->SetBinContent(uMom+1,ps_val);

    PARS[0] = 0;
    ps_val = atan(kstar/EffRangeExp(&kstar,PARS));
    hPS_s->SetBinContent(uMom+1,ps_val);
  }

  hPS_cs->Draw();
  hPS_s->Draw("same");

  myapp->Run();
}


void pp_pSp_Decay(int SEED, int MilNumIter){
  DLM_DecayMatrix decMat;
  TString OutputFileName = TString::Format("%s/CatsFiles/DecaySmear/Decay_matrix_pp_pSp_RS%i.root",GetCernBoxDimi(),SEED);
  decMat.SetFileName(OutputFileName.Data());
  decMat.SetHistoName("hRes_pp_pSp");
  decMat.SetBins(1024,0,1024);
  decMat.SetNumDaughters1(1);
  decMat.SetMotherMass1(Mass_p);
  decMat.SetNumDaughters2(2);
  decMat.SetDaughterMass2(0,Mass_p);
  decMat.SetDaughterMass2(1,Mass_pi0);
  decMat.SetMotherMass2(Mass_Sch);
  decMat.SetMeanMomentum(0);
  decMat.SetMomentumSpread(350);
  decMat.Run(SEED,MilNumIter*1000*1000);
}

void piXi1530_piXi_Decay(int SEED, int MilNumIter){
  DLM_DecayMatrix decMat;
  TString OutputFileName = TString::Format("%s/CatsFiles/DecaySmear/Decay_matrix_piXim_piXim1530_RS%i.root",GetCernBoxDimi(),SEED);
  decMat.SetFileName(OutputFileName.Data());
  decMat.SetHistoName("hRes_piXim_piXim1530");
  decMat.SetBins(1024,0,1024);

  decMat.SetNumDaughters1(1);
  decMat.SetMotherMass1(Mass_pic);

  decMat.SetNumDaughters2(2);
  decMat.SetDaughterMass2(0,Mass_Xim);
  decMat.SetDaughterMass2(1,Mass_pic);
  decMat.SetMotherMass2(Mass_Xim1530);
  decMat.SetMeanMomentum(0);
  decMat.SetMomentumSpread(350);
  decMat.Run(SEED,MilNumIter*1000*1000);
}

void KplusXi1530_KplusXi_Decay(int SEED, int MilNumIter){
  DLM_DecayMatrix decMat;
  TString OutputFileName = TString::Format("%s/CatsFiles/DecaySmear/Decay_matrix_KplusXim_KplusXim1530_RS%i.root",GetCernBoxDimi(),SEED);
  decMat.SetFileName(OutputFileName.Data());
  decMat.SetHistoName("hRes_KplusXim_KplusXim1530");
  decMat.SetBins(1024,0,1024);

  decMat.SetNumDaughters1(1);
  decMat.SetMotherMass1(Mass_Kch);

  decMat.SetNumDaughters2(2);
  decMat.SetDaughterMass2(0,Mass_Xim);
  decMat.SetDaughterMass2(1,Mass_pic);
  decMat.SetMotherMass2(Mass_Xim1530);
  decMat.SetMeanMomentum(0);
  decMat.SetMomentumSpread(350);
  decMat.Run(SEED,MilNumIter*1000*1000);
}

void pL_LL_Decay(int SEED, int MilNumIter){
  DLM_DecayMatrix decMat;
  //TString OutputFileName = TString::Format("%s/CatsFiles/DecaySmear/Decay_matrix_pL_LL_RS%i.root",GetCernBoxDimi(),SEED);
  TString OutputFileName = TString::Format("%s/OtherTasks/pL_LL_Decay/Decay_matrix_pL_LL_RS%i.root",GetFemtoOutputFolder(),SEED);
  decMat.SetFileName(OutputFileName.Data());
  decMat.SetHistoName("hRes_pL_LL");
  decMat.SetBins(1024,0,1024);
  decMat.SetNumDaughters1(1);
  decMat.SetMotherMass1(Mass_L);
  decMat.SetNumDaughters2(2);
  decMat.SetDaughterMass2(0,Mass_p);
  decMat.SetDaughterMass2(1,Mass_pi0);
  decMat.SetMotherMass2(Mass_L);
  decMat.SetMeanMomentum(0);
  decMat.SetMomentumSpread(350);
  decMat.Run(SEED,MilNumIter*1000*1000);
}

void pp_LL_Decay(int SEED, int MilNumIter){
  DLM_DecayMatrix decMat;
  //TString OutputFileName = TString::Format("%s/CatsFiles/DecaySmear/Decay_matrix_pL_LL_RS%i.root",GetCernBoxDimi(),SEED);
  TString OutputFileName = TString::Format("%s/OtherTasks/pL_LL_Decay/Decay_matrix_pp_LL_RS%i.root",GetFemtoOutputFolder(),SEED);
  decMat.SetFileName(OutputFileName.Data());
  decMat.SetHistoName("hRes_pL_LL");
  decMat.SetBins(1024,0,1024);
  decMat.SetNumDaughters1(2);
  decMat.SetDaughterMass1(0,Mass_p);
  decMat.SetDaughterMass1(1,Mass_pi0);
  decMat.SetMotherMass1(Mass_L);
  decMat.SetNumDaughters2(2);
  decMat.SetDaughterMass2(0,Mass_p);
  decMat.SetDaughterMass2(1,Mass_pi0);
  decMat.SetMotherMass2(Mass_L);
  decMat.SetMeanMomentum(0);
  decMat.SetMomentumSpread(350);
  decMat.Run(SEED,MilNumIter*1000*1000); 
}




void pn_Xchecks_1(){
  const unsigned NumMomBins = 80;
  const double kMin = 0;
  const double kMax = 320;
  const double SourceSize = 1.0;

  DLM_CommonAnaFunctions AnalysisObject;
  //change the path to the location where you have downloaded the relevant CernBox folders
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles/",GetCernBoxDimi()));
  TString SourceDescription = "Gauss";

  //s and p waves
  CATS Kitty_pn;
  Kitty_pn.SetMomBins(NumMomBins,kMin,kMax);
  AnalysisObject.SetUpCats_pp(Kitty_pn,"AV18_pn",SourceDescription,0,0);
  Kitty_pn.SetAnaSource(0, SourceSize);

  //Kitty_pn.RemoveShortRangePotential(0,1);
  Kitty_pn.SetOnlyNumericalPw(0,0);
  Kitty_pn.SetChannelWeight(0,1);
  Kitty_pn.SetChannelWeight(1,0);
  Kitty_pn.SetChannelWeight(2,0);
  Kitty_pn.SetChannelWeight(3,0);

  Kitty_pn.KillTheCat();

  //Kitty_pn.SetThetaDependentSource(false);
  //CATSparameters cPars(CATSparameters::tSource, 1, true);
  //cPars.SetParameter(0, SourceSize);
  //Kitty_pn.SetAnaSource(GaussSource, cPars);
  //Kitty_pn.SetUseAnalyticSource(true);
  //Kitty_pn.SetMomentumDependentSource(false);
  //Kitty_pn.SetExcludeFailedBins(false);
  //Kitty_pn.SetQ1Q2(1);
  //Kitty_pn.SetPdgId(2212, 2212);
  //Kitty_pn.SetRedMass(0.5 * Mass_p);
  //Kitty_pn.SetNumChannels(4);
  //Kitty_pn.SetNumPW(0, 2);
  //Kitty_pn.SetNumPW(1, 2);
  //Kitty_pn.SetNumPW(2, 2);
  //Kitty_pn.SetNumPW(3, 2);
  //Kitty_pn.SetSpin(0, 0);
  //Kitty_pn.SetSpin(1, 1);
  //Kitty_pn.SetSpin(2, 1);
  //Kitty_pn.SetSpin(3, 1);
  //Kitty_pn.SetChannelWeight(0, 3. / 12.);
  //Kitty_pn.SetChannelWeight(1, 1. / 12.);
  //Kitty_pn.SetChannelWeight(2, 3. / 12.);
  //Kitty_pn.SetChannelWeight(3, 5. / 12.);

  TGraph gCk_pn;
  gCk_pn.SetName("gCk_pn");

  TGraph gPs1S0_pn;
  gPs1S0_pn.SetName("gPs1S0_pn");

  for(unsigned uMom=0; uMom<NumMomBins; uMom++){
    double MOM = Kitty_pn.GetMomentum(uMom);
    double Ck_val = Kitty_pn.GetCorrFun(uMom);
    gCk_pn.SetPoint(uMom,MOM,Ck_val);
    gPs1S0_pn.SetPoint(uMom,MOM,Kitty_pn.GetPhaseShift(uMom,0,0));
  }

  TString OutputFileName = TString::Format("%s/OtherTasks/pn_Xchecks_1/fOutput.root",GetFemtoOutputFolder());
  TFile fOutput(OutputFileName,"recreate");

  gCk_pn.Write();
  gPs1S0_pn.Write();

}

void Test_chiral_2023(TString Singlet, TString Triplet){
  CATS Kitty;
  Kitty.SetMomBins(40,0,200);
  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());

  CATSparameters cPars(CATSparameters::tSource, 1, true);
  cPars.SetParameter(0, 1.2);
  Kitty.SetAnaSource(GaussSource, cPars);
  DLM_Histo<complex<double>> ***ExternalWF = NULL;
  ExternalWF = Init_pL_Haidenbauer2023(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data(),Kitty,Singlet.Data(),Triplet.Data());
  for(unsigned uCh=0; uCh<Kitty.GetNumChannels(); uCh++){
    Kitty.SetChannelWeight(uCh,0);
  }
  Kitty.SetChannelWeight(0,1./4.);
  Kitty.SetChannelWeight(1,3./4.);
  Kitty.SetExternalWaveFunction(0, 0, ExternalWF[0][0][0], ExternalWF[1][0][0]);
  Kitty.SetExternalWaveFunction(1, 0, ExternalWF[0][1][0], ExternalWF[1][1][0]);
  Kitty.KillTheCat();

  CATS FatKitty;
  FatKitty.SetMomBins(40,0,320);
  FatKitty.SetAnaSource(GaussSource, cPars);
  AnalysisObject.SetUpCats_pL(FatKitty,"Chiral2023_"+Singlet+"_and_"+Triplet,"Gauss",0,0);
  FatKitty.SetAnaSource(0,1.2);
  FatKitty.KillTheCat();


  TString OutputFileName = TString::Format("%s/OtherTasks/Test_chiral_2023/%s_and_%s.root",GetFemtoOutputFolder(),Singlet.Data(),Triplet.Data());
  TFile fOutput(OutputFileName,"recreate");

  TGraph gCk;
  gCk.SetName(Singlet+"_and_"+Triplet);

  TGraph gCkFat;
  gCkFat.SetName(Singlet+"_FAT_"+Triplet);

  for(unsigned uBin=0; uBin<Kitty.GetNumMomBins(); uBin++){
    gCk.SetPoint(uBin, Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
    gCkFat.SetPoint(uBin, FatKitty.GetMomentum(uBin), FatKitty.GetCorrFun(uBin));
  }

  gCk.Write();
  gCkFat.Write();

  //if(Singlet=="NLO19_600_291_s" && Triplet=="NLO19_600_141_s"){
    CATS Original;
    Original.SetMomBins(40,0,320);
    Original.SetAnaSource(GaussSource, cPars);
    AnalysisObject.SetUpCats_pL(Original,"Chiral_Coupled_SPD","Gauss",11600,0);
    Original.SetAnaSource(0,1.2);
    Original.KillTheCat();
    TGraph gCkO;
    gCkO.SetName(Singlet+"_Paper_"+Triplet);
    for(unsigned uBin=0; uBin<Original.GetNumMomBins(); uBin++){
      gCkO.SetPoint(uBin, Original.GetMomentum(uBin), Original.GetCorrFun(uBin));
    }
    gCkO.Write();
  //}

}

void printFloatBits(float value) {
    // Create a pointer to the float
    unsigned int* floatBits = reinterpret_cast<unsigned int*>(&value);

    // Use std::bitset to print the bits
    std::bitset<sizeof(float) * 8> bits(*floatBits);
    std::cout << "Float: " << value << "\nBits: " << bits << std::endl;
}

void StableDisto_DlmRan_Test(){
  const float srcsize = 1.0;
  std::vector<float> alpha_par;
  alpha_par.push_back(1);
  alpha_par.push_back(1.2);
  alpha_par.push_back(1.4);
  alpha_par.push_back(1.6);
  alpha_par.push_back(1.8);
  alpha_par.push_back(2.0);

  const unsigned NumRndPts = 1000*1000;
  
  DLM_Random RanGen(11);

  TH1F** hLevy = new TH1F* [alpha_par.size()];
  TH1F** hLevy3D = new TH1F* [alpha_par.size()];
  TH1F** hLevy6D = new TH1F* [alpha_par.size()];
  TH1F** hLevySrc1 = new TH1F* [alpha_par.size()];
  const unsigned NumRadBins = 1024;
  const float rMax = 16;

  for(unsigned uLev=0; uLev<alpha_par.size(); uLev++){
    float lev = alpha_par.at(uLev);
    float rx1,ry1,rz1;
    float rx2,ry2,rz2;
    hLevy[uLev] = new TH1F(TString::Format("hLevy_%.2f",lev),TString::Format("hLevy_%.2f",lev),NumRadBins,-rMax*0.5,rMax*0.5);
    hLevy3D[uLev] = new TH1F(TString::Format("hLevy3D_%.2f",lev),TString::Format("hLevy3D_%.2f",lev),NumRadBins,0,rMax);
    hLevy6D[uLev] = new TH1F(TString::Format("hLevy6D_%.2f",lev),TString::Format("hLevy6D_%.2f",lev),NumRadBins,0,rMax);
    hLevySrc1[uLev] = new TH1F(TString::Format("hLevySrc1_%.2f",lev),TString::Format("hLevySrc1_%.2f",lev),NumRadBins,0,rMax);
    for(unsigned uIter=0; uIter<NumRndPts; uIter++){
      rx1 = RanGen.Stable(lev, 0, srcsize);
      ry1 = RanGen.Stable(lev, 0, srcsize);
      rz1 = RanGen.Stable(lev, 0, srcsize);
      rx2 = RanGen.Stable(lev, 0, srcsize);
      ry2 = RanGen.Stable(lev, 0, srcsize);
      rz2 = RanGen.Stable(lev, 0, srcsize);
      hLevy[uLev]->Fill(rx1);
      hLevy[uLev]->Fill(ry1);
      hLevy[uLev]->Fill(rz1);
      hLevy[uLev]->Fill(rx2);
      hLevy[uLev]->Fill(ry2);
      hLevy[uLev]->Fill(rz2);      

      hLevy3D[uLev]->Fill(sqrt(rx1*rx1+ry1*ry1+rz1*rz1));
      hLevy3D[uLev]->Fill(sqrt(rx2*rx2+ry2*ry2+rz2*rz2));

      hLevy6D[uLev]->Fill(sqrt(pow(rx1-rx2,2.)+pow(ry1-ry2,2.)+pow(rz1-rz2,2.)));
    }
    double PARS[5];
    for(unsigned uRad=0; uRad<NumRadBins; uRad++){
      double RAD = hLevySrc1[uLev]->GetBinCenter(uRad+1);
      PARS[1] = RAD;
      PARS[3] = srcsize;
      PARS[4] = lev;
      hLevySrc1[uLev]->SetBinContent(uRad+1, LevySource3D_2particle(PARS));
    }
  }

  

  TString OutputFileName = TString::Format("%s/OtherTasks/StableDisto_DlmRan_Test/fOutput.root",GetFemtoOutputFolder());
  TFile fOutput(OutputFileName,"recreate");
  for(unsigned uLev=0; uLev<alpha_par.size(); uLev++){
    hLevy[uLev]->Sumw2();
    hLevy[uLev]->Scale(1./hLevy[uLev]->Integral(),"width");
    hLevy3D[uLev]->Scale(1./hLevy3D[uLev]->Integral(),"width");
    hLevy6D[uLev]->Scale(hLevySrc1[uLev]->Integral()/hLevy6D[uLev]->Integral());


    hLevy[uLev]->Write();
    hLevy3D[uLev]->Write();
    hLevy6D[uLev]->Write();
    hLevySrc1[uLev]->Write();
  }

}


//from kstar 0 to 200 typically 2x less CA
//at kstar = 0 we have some larger amount, say 70% of CA, say f(x) = 2/3∙exp(−x/300) for CA weight
void pp_Gauss_Cauchy_Mix(){
  const double NonCaSrcSize = 1.2;//gauss
  const double CaSrcSize = 1.2*sqrt(2);//cauchy
  
  const unsigned NumMomBins = 80;
  const double kMin = 0;
  const double kMax = 320;

  const unsigned NumRadBins = 1024;
  const double rMin = 0;
  const double rMax = 16; 

  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles/",GetCernBoxDimi()));


  CATS Gauss_Cat;
  CATS Cauchy_Cat;

  Gauss_Cat.SetMomBins(NumMomBins,kMin,kMax);
  Cauchy_Cat.SetMomBins(NumMomBins,kMin,kMax);

  AnalysisObject.SetUpCats_pp(Gauss_Cat,"AV18","Gauss",0,0);
  AnalysisObject.SetUpCats_pp(Cauchy_Cat,"AV18","Cauchy",0,0);

  Gauss_Cat.SetAnaSource(0, NonCaSrcSize);
  Cauchy_Cat.SetAnaSource(0, CaSrcSize);

  Gauss_Cat.KillTheCat();
  Cauchy_Cat.KillTheCat();

  TGraph gGauss_Cat;
  gGauss_Cat.SetName("gGauss_Cat");
  TGraph gCauchy_Cat;
  gCauchy_Cat.SetName("gCauchy_Cat");

  TGraph gGauss_Src;
  gGauss_Src.SetName("gGauss_Src");
  TGraph gCauchy_Src;
  gCauchy_Src.SetName("gCauchy_Src");

  for(unsigned uBin=0; uBin<NumMomBins; uBin++){
    double kstar = Gauss_Cat.GetMomentum(uBin);
    gGauss_Cat.SetPoint(uBin, kstar, Gauss_Cat.GetCorrFun(uBin));
    gCauchy_Cat.SetPoint(uBin, kstar, Cauchy_Cat.GetCorrFun(uBin));
  }


  for(unsigned uRad=0; uRad<NumRadBins; uRad++){
    double rstar = (rMax-rMin)/double(NumRadBins) * (0.5 + double(uRad));
    gGauss_Src.SetPoint(uRad, rstar, Gauss_Cat.EvaluateTheSource(10, rstar, 0));
    gCauchy_Src.SetPoint(uRad, rstar, Cauchy_Cat.EvaluateTheSource(10, rstar, 0));
  }

  TString OutputFileName = TString::Format("%s/OtherTasks/pp_Gauss_Cauchy_Mix/fOutput.root",GetFemtoOutputFolder());
  TFile fOutput(OutputFileName,"recreate");
  gGauss_Cat.Write();
  gCauchy_Cat.Write();
  gGauss_Src.Write();
  gCauchy_Src.Write();
  
}
/*
void Unfold_test1(){

  unsigned NumMomBins = 80;
  double kMin = 0;
  double kMax = 320;
  double kMaxX = 960;
  unsigned NumMomBinsXX = 240;
  double SourceSize = 1.2;
  double ERROR_FACTOR = 25;

  //ME true
  TH1F* hRef_True;
  //ME CATS smear
  TH1F* hRef_Det;
  //ME brute force smear
  TH1F* hRef_BS;
  //SE true
  TH1F* hSame_True;
  //ME CATS smear
  TH1F* hSame_Det;
  //ME brute force smear
  TH1F* hSame_BS;

  //Ck theo
  TH1F* hCk_Theo;
  //Ck true
  TH1F* hCk_True;
  //Ck CATS smear of SE/ME
  TH1F* hCk_Det;
  //Ck CATS smear with PS
  TH1F* hCk_Decomp;

  const TString ME_FileName = TString::Format("%s/CatsFiles/ExpData/ALICE_pp_13TeV_HM/Sample10HM/AnalysisResults.root",GetCernBoxDimi());
  const TString DirName = "HMResults";
  //this is the reweighted ME sample
  const TString HistoName = "MEDist_Particle0_Particle0";
  TFile* fIn_ME = new TFile(ME_FileName,"read");
  TDirectoryFile* dir_ME=(TDirectoryFile*)(fIn_ME->FindObjectAny(DirName));
  TList* list1_tmp;
  dir_ME->GetObject(DirName,list1_tmp);
  TList* list2_tmp = (TList*)list1_tmp->FindObject("Particle0_Particle0");

  TH1F* hD_ME;
  hD_ME = (TH1F*)list2_tmp->FindObject(HistoName);
  printf("hD_ME %p\n",hD_ME);
  gROOT->cd();
  hRef_True = new TH1F("hRef_True","hRef_True",hD_ME->GetNbinsX(),
  hD_ME->GetXaxis()->GetBinLowEdge(1)*1000,hD_ME->GetXaxis()->GetBinLowEdge(hD_ME->GetNbinsX()+1)*1000);

  hCk_Theo = new TH1F("hCk_Theo","hCk_Theo",hD_ME->GetNbinsX(),
  hD_ME->GetXaxis()->GetBinLowEdge(1)*1000,hD_ME->GetXaxis()->GetBinLowEdge(hD_ME->GetNbinsX()+1)*1000);

  for(unsigned uBin=0; uBin<hD_ME->GetNbinsX(); uBin++){
    hRef_True->SetBinContent(uBin+1,hD_ME->GetBinContent(uBin+1));
    //hRef_True->SetBinError(uBin+1,hD_ME->GetBinError(uBin+1));
  }
  hRef_True->Sumw2();

  hRef_Det = (TH1F*)hRef_True->Clone("hRef_Det");
  hRef_BS = (TH1F*)hRef_True->Clone("hRef_BS");
  TH1F* hRef_X = new TH1F("hRef_X","hRef_X",NumMomBinsXX,kMin,kMaxX);
//fill in this limited distribution, to try and unfold with less binsw
 // for(unsigned uBin=0; uBin<NumMomBinsXX; uBin++){
  //  hRef_X
 // }

  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles/",GetCernBoxDimi()));

  CATS Kitty_pp;
  Kitty_pp.SetMomBins(NumMomBins,kMin,kMax);
  AnalysisObject.SetUpCats_pp(Kitty_pp,"AV18","Gauss",0,0);
  Kitty_pp.SetAnaSource(0,SourceSize);
  Kitty_pp.SetNotifications(CATS::nWarning);//avoid printing stuff
  //evaluate the theoretical correlation function (NO lambda pars, mom resolution etc.)
  Kitty_pp.KillTheCat();

  DLM_Ck DlmCk_pp(Kitty_pp.GetNumSourcePars(),0,Kitty_pp);
  DlmCk_pp.SetCutOff(300,700);

  TH2F* hReso_pp = AnalysisObject.GetResolutionMatrix("pp13TeV_HM_DimiJun20","pp");
  DLM_CkDecomposition CkDecomp_pp("pp",0,DlmCk_pp,hReso_pp);
  CkDecomp_pp.AddPhaseSpace(hRef_True);
  CkDecomp_pp.Update();

  TRandom3 rangen(11);

  hSame_True = (TH1F*)hRef_True->Clone("hSame_True");
  for(unsigned uBin=0; uBin<hSame_True->GetNbinsX(); uBin++){
    double kstar = hSame_True->GetBinCenter(uBin+1);
    double MeVal = hRef_True->GetBinContent(uBin+1);
    double SeVal = -1;
    double SeErr = hRef_True->GetBinError(uBin+1)*ERROR_FACTOR;
    //double CkVal = 0;
    while(SeVal<0){
      SeVal = rangen.Gaus(MeVal*DlmCk_pp.Eval(kstar), SeErr);
    }
    hSame_True->SetBinContent(uBin+1, SeVal);
    hSame_True->SetBinError(uBin+1, SeErr);
    hCk_Theo->SetBinContent(uBin+1, DlmCk_pp.Eval(kstar));
  }

  hSame_Det = (TH1F*)hSame_True->Clone("hSame_Det");
  hSame_BS = (TH1F*)hSame_True->Clone("hSame_BS");

  hCk_True = (TH1F*)hSame_True->Clone("hCk_True");
  hCk_Decomp = new TH1F("hCk_Decomp","hCk_Decomp",NumMomBins,kMin,kMax);
  hCk_True->Divide(hRef_True);

  for(unsigned uBin=0; uBin<NumMomBins; uBin++){
    double kstar = Kitty_pp.GetMomentum(uBin);
    hCk_Decomp->SetBinContent(uBin+1, CkDecomp_pp.EvalCk(kstar));
    //hCk_Decomp->SetBinContent(uBin+1, DlmCk_pp.Eval(kstar));
 
  }

  DLM_Ck DlmSe_pp(NumMomBinsXX, kMin, kMaxX);
  DLM_Ck DlmMe_pp(NumMomBinsXX, kMin, kMaxX);
  for(unsigned uBin=0; uBin<NumMomBinsXX; uBin++){
    DlmSe_pp.SetBinContent(uBin, hSame_True->GetBinContent(uBin+1));
    DlmMe_pp.SetBinContent(uBin, hRef_True->GetBinContent(uBin+1));
  }
  DLM_CkDecomposition SeDecomp_pp("pp",0,DlmSe_pp,hReso_pp);
  DLM_CkDecomposition MeDecomp_pp("pp",0,DlmMe_pp,hReso_pp);

  hCk_Det = new TH1F("hCk_Det","hCk_Det",NumMomBinsXX,kMin,kMaxX);
  for(unsigned uBin=0; uBin<NumMomBinsXX; uBin++){
    double kstar = DlmSe_pp.GetBinCenter(0, uBin);
    hCk_Det->SetBinContent(uBin+1,SeDecomp_pp.EvalCk(kstar)/MeDecomp_pp.EvalCk(kstar));
    hSame_Det->SetBinContent(uBin+1,SeDecomp_pp.EvalCk(kstar));
    hRef_Det->SetBinContent(uBin+1,MeDecomp_pp.EvalCk(kstar));
  }

  double bin_width = hSame_True->GetBinWidth(1);
  double bin_width_reso = hReso_pp->GetXaxis()->GetBinWidth(1);
  int rescale_factor = TMath::Nint(bin_width/bin_width_reso);
  printf("bw, bwr, rf = %.2f, %.2f, %i\n", bin_width, bin_width_reso, rescale_factor);
  TH2F* hReso_pp_rsc = (TH2F*)hReso_pp->Clone("hReso_pp_rsc");
  hReso_pp_rsc->Rebin2D(rescale_factor,rescale_factor);
  TH1F* hCk_BS2 = (TH1F*)hCk_True->Clone("hCk_BS2");
  //brute force smear, should actually be the most accurate
  for(unsigned uBinY=0; uBinY<NumMomBinsXX; uBinY++){
    double Sk_val = 0;
    double Rk_val = 0;
    double Ck_val = 0;
    double Norm_val = 0;
    double NormCk_val = 0;
    for(unsigned uBinX=0; uBinX<NumMomBinsXX; uBinX++){
      Sk_val += hReso_pp_rsc->GetBinContent(uBinX+1, uBinY+1)*hSame_True->GetBinContent(uBinX+1);
      Rk_val += hReso_pp_rsc->GetBinContent(uBinX+1, uBinY+1)*hRef_True->GetBinContent(uBinX+1);
      Ck_val += hReso_pp_rsc->GetBinContent(uBinX+1, uBinY+1)*hRef_True->GetBinContent(uBinX+1)*hCk_True->GetBinContent(uBinX+1);
      Norm_val += hReso_pp_rsc->GetBinContent(uBinX+1, uBinY+1);
      NormCk_val += hReso_pp_rsc->GetBinContent(uBinX+1, uBinY+1)*hRef_True->GetBinContent(uBinX+1);
    }
    Sk_val /= Norm_val;
    Rk_val /= Norm_val;
    Ck_val /= NormCk_val;
    hSame_BS->SetBinContent(uBinY+1, Sk_val);
    hSame_BS->SetBinError(uBinY+1, sqrt(Sk_val)*ERROR_FACTOR);
    hRef_BS->SetBinContent(uBinY+1, Rk_val);
    hRef_BS->SetBinError(uBinY+1, sqrt(Rk_val));
    hCk_BS2->SetBinContent(uBinY+1, Ck_val);
  }
  TH1F* hCk_BS = (TH1F*)hSame_BS->Clone("hCk_BS");
  hCk_BS->Divide(hRef_BS);


  // regularize curvature
  TUnfold::ERegMode regMode =
     TUnfold::kRegModeCurvature;
 // preserve the area
  TUnfold::EConstraint constraintMode=
     TUnfold::kEConstraintArea;
  // bin content is divided by the bin width
  TUnfoldDensity::EDensityMode densityFlags=
     TUnfoldDensity::kDensityModeBinWidth;
  // set up matrix of migrations
  TUnfoldDensity unfold1(hReso_pp,TUnfold::kHistMapOutputHoriz,
                        regMode,constraintMode,densityFlags);

    // define the input vector (the measured data distribution)
  unfold1.SetInput(hRef_BS);
  // run the unfolding
  Int_t nScan=30;
  TSpline *logTauX,*logTauY;
  TGraph *lCurve;
  // this method scans the parameter tau and finds the kink in the L curve
  // finally, the unfolding is done for the best choice of tau
  Int_t iBest=unfold1.ScanLcurve(nScan,0.,0.,&lCurve,&logTauX,&logTauY);
  cout<<"chi**2="<<unfold1.GetChi2A()<<"+"<<unfold1.GetChi2L()
      <<" / "<<unfold1.GetNdf()<<"\n";




  TFile fOutput(TString::Format("%s/OtherTasks/Unfold_test1/Unfold_test1.root",GetFemtoOutputFolder()), "recreate");
  hSame_True->Write();
  hSame_Det->Write();
  hSame_BS->Write();
  hRef_True->Write();
  hRef_Det->Write();
  hRef_BS->Write();
  hCk_Theo->Write();
  hCk_True->Write();
  hReso_pp->Write();
  hCk_Decomp->Write();
  hCk_Det->Write();
  hCk_BS->Write();
  //should actually be the most accurate
  hCk_BS2->Write();
  //hUnfolded1->Write();
  
  delete hRef_True;
  delete hRef_Det;
  delete hRef_BS;
  delete hSame_True;
  delete hSame_Det;
  delete hSame_BS;
  delete hCk_True;
  delete hCk_Det;
  delete hCk_Theo;
  delete hCk_BS;
  delete hCk_BS2;
  //delete hUnfolded1; 
}
*/
void generate_pp_Neelima(){

  double kMin = 1.25;
  double kMax = 301.25;
  unsigned NumMomBins = 120;

  double rMin = 0.7-0.01/2.;
  double rMax = 1.4+0.01/2.;
  unsigned NumRadBins = 71;

  CATS Kitty_pp;
  Kitty_pp.SetMomBins(NumMomBins, kMin, kMax);

  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());

  AnalysisObject.SetUpCats_pp(Kitty_pp,"AV18","Gauss",0,0);
  //
  Kitty_pp.SetNotifications(CATS::nWarning);//avoid printing stuff
  //evaluate the theoretical correlation function (NO lambda pars, mom resolution etc.)
  Kitty_pp.SetEpsilonConv(1e-8);
  Kitty_pp.SetEpsilonProp(1e-8);
  //Kitty_pp.KillTheCat();

  TH2F* h2D_pp_spd = new TH2F("h2D_pp_spd","h2D_pp_spd", NumMomBins, kMin, kMax, NumRadBins, rMin, rMax);
  TGraph2D g2D_pp_spd;
  g2D_pp_spd.SetName("g2D_pp_spd");

  for(unsigned uRad=0; uRad<NumRadBins; uRad++){
    double rad = h2D_pp_spd->GetYaxis()->GetBinCenter(uRad+1);
    printf("rad = %.3f\n",rad);
    if(uRad%10==0)
      Kitty_pp.SetAnaSource(0,rad);
    else
      Kitty_pp.SetAnaSource(0,rad,true);
    Kitty_pp.KillTheCat();

    for(unsigned uMom=0; uMom<NumMomBins; uMom++){
      double CkVal = Kitty_pp.GetCorrFun(uMom);
      double mom = Kitty_pp.GetMomentum(uMom);
      h2D_pp_spd->SetBinContent(uMom+1,uRad+1,CkVal);
      g2D_pp_spd.SetPoint(uRad*NumRadBins+uMom,mom,rad,CkVal);
    }
  }



  TH2F* h2D_pp_s = new TH2F("h2D_pp_s","h2D_pp_s", NumMomBins, kMin, kMax, NumRadBins, rMin, rMax);
  TGraph2D g2D_pp_s;
  g2D_pp_s.SetName("g2D_pp_s");


    Kitty_pp.RemoveShortRangePotential(0,2);
    Kitty_pp.RemoveShortRangePotential(1,1);
    Kitty_pp.RemoveShortRangePotential(2,1);
    Kitty_pp.RemoveShortRangePotential(3,1);

  for(unsigned uRad=0; uRad<NumRadBins; uRad++){
    double rad = h2D_pp_s->GetYaxis()->GetBinCenter(uRad+1);
    printf("rad = %.3f\n",rad);
    if(uRad%10==0)
      Kitty_pp.SetAnaSource(0,rad);
    else
      Kitty_pp.SetAnaSource(0,rad,true);
    Kitty_pp.KillTheCat();

    for(unsigned uMom=0; uMom<NumMomBins; uMom++){
      double CkVal = Kitty_pp.GetCorrFun(uMom);
      double mom = Kitty_pp.GetMomentum(uMom);
      h2D_pp_s->SetBinContent(uMom+1,uRad+1,CkVal);
      g2D_pp_s.SetPoint(uRad*NumRadBins+uMom,mom,rad,CkVal);
    }
  }



  TFile fOutput(TString::Format("%s/OtherTasks//generate_pp_Neelima.root",GetFemtoOutputFolder()), "recreate");
  h2D_pp_spd->Write();
  g2D_pp_spd.Write();
  h2D_pp_s->Write();
  g2D_pp_s.Write();

  //g2D_pp.SetPoint(0,-1,-1,1);
  //g2D_pp.SetPoint(1,1,-1,1);
  //g2D_pp.SetPoint(2,-1,1,2);
  //g2D_pp.SetPoint(3,1,1,2);

  //printf("0, -1: %.2f\n", g2D_pp.Interpolate(0,-1));
  //printf(" 1, 0: %.2f\n", g2D_pp.Interpolate(1,0));
  //printf(" 0, 1: %.2f\n", g2D_pp.Interpolate(0,1));
  //printf("-1, 0: %.2f\n", g2D_pp.Interpolate(-1,0));
  //printf(" 0, 0: %.2f\n", g2D_pp.Interpolate(0,0));
  //printf(" 0, 0.5: %.2f\n", g2D_pp.Interpolate(0,0.5));

  delete h2D_pp_spd;
  delete h2D_pp_s;
}

//run using the old code, to set a baseline to reproduce
void TestLambdaKstar_Baseline(TString OutputFileName){
  unsigned NumMomBins_pp = 75;
  unsigned NumMomBins_pL = 25;
  double kMin = 0;
  double kMax = 300;

  CATS Kitty_pp;
  Kitty_pp.SetMomBins(NumMomBins_pp, kMin, kMax);
  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());
  AnalysisObject.SetUpCats_pp(Kitty_pp,"AV18","Gauss",0,0);
  Kitty_pp.KillTheCat();

  CATS Kitty_pL;
  Kitty_pL.SetMomBins(NumMomBins_pL, kMin, kMax);
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());
  AnalysisObject.SetUpCats_pL(Kitty_pL,"Usmani","Gauss",0,0);
  Kitty_pL.KillTheCat();  

  DLM_Ck Ck_Kitty_pp(Kitty_pp.GetNumSourcePars(),0,Kitty_pp,NumMomBins_pp*3,kMin,kMax*3);
  DLM_Ck Ck_Kitty_pL(Kitty_pL.GetNumSourcePars(),0,Kitty_pL,NumMomBins_pL*3,kMin,kMax*3);

  Ck_Kitty_pp.SetCutOff(300,700);
  Ck_Kitty_pL.SetCutOff(300,700);

  Ck_Kitty_pp.Update();
  Ck_Kitty_pL.Update();

  TH2F* hReso_pp = AnalysisObject.GetResolutionMatrix("pp13TeV_HM_DimiJun20","pp");
  TH2F* hFeed_pp_pL = AnalysisObject.GetResidualMatrix("pp","pLambda");

  DLM_CkDecomposition CkDecomp_pp("pp",3,Ck_Kitty_pp,hReso_pp);
  DLM_CkDecomposition CkDecomp_pL("pLambda",2,Ck_Kitty_pL,NULL);


  CkDecomp_pL.AddContribution(0,0.4,DLM_CkDecomposition::cFeedDown);
  CkDecomp_pL.AddContribution(1,0.1,DLM_CkDecomposition::cFake);

  CkDecomp_pp.AddContribution(0,0.2,DLM_CkDecomposition::cFeedDown,&CkDecomp_pL,hFeed_pp_pL);
  CkDecomp_pp.AddContribution(1,0.1,DLM_CkDecomposition::cFeedDown);
  CkDecomp_pp.AddContribution(2,0.01,DLM_CkDecomposition::cFake);

  CkDecomp_pp.Update();
  CkDecomp_pL.Update();

  TFile fOutput(TString::Format("%s/OtherTasks/TestLambdaKstar/%s.root",GetFemtoOutputFolder(),OutputFileName.Data()), "recreate");
  TGraph g_cats_pp;
  g_cats_pp.SetName("g_cats_pp");
  g_cats_pp.SetLineWidth(6);
  g_cats_pp.SetLineColor(kGreen+1);
  TGraph g_ck_pp;
  g_ck_pp.SetName("g_ck_pp");
  g_ck_pp.SetLineWidth(5);
  g_ck_pp.SetLineColor(kBlue);
  g_ck_pp.SetLineStyle(2);
  TGraph g_dec_pp;
  g_dec_pp.SetName("g_dec_pp");
  g_dec_pp.SetLineWidth(4);
  g_dec_pp.SetLineColor(kRed+1);
  for(unsigned uMom=0; uMom<NumMomBins_pp*3; uMom++){
    double kstar = Ck_Kitty_pp.GetBinCenter(0, uMom);
    if(kstar<kMax) g_cats_pp.SetPoint(uMom, kstar, Kitty_pp.GetCorrFun(uMom));
    g_ck_pp.SetPoint(uMom, kstar, Ck_Kitty_pp.GetBinContent(uMom));
    g_dec_pp.SetPoint(uMom, kstar, CkDecomp_pp.EvalCk(kstar));
  }
  g_cats_pp.Write();
  g_ck_pp.Write();
  g_dec_pp.Write();
}

//use the new code with fixed lambda pars, to see if you reproduce your baseline
void TestLambdaKstar_NewCode(){

}

//introduce some kstar dependance, to see if it makes sence
void TestLambdaKstar_KstarDep(TString OutputFileName){
  unsigned NumMomBins_pp = 75;
  unsigned NumMomBins_pL = 25;
  double kMin = 0;
  double kMax = 300;

  CATS Kitty_pp;
  Kitty_pp.SetMomBins(NumMomBins_pp, kMin, kMax);
  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());
  AnalysisObject.SetUpCats_pp(Kitty_pp,"AV18","Gauss",0,0);
  Kitty_pp.KillTheCat();

  CATS Kitty_pL;
  Kitty_pL.SetMomBins(NumMomBins_pL, kMin, kMax);
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());
  AnalysisObject.SetUpCats_pL(Kitty_pL,"Usmani","Gauss",0,0);
  Kitty_pL.KillTheCat();  

  DLM_Ck Ck_Kitty_pp(Kitty_pp.GetNumSourcePars(),0,Kitty_pp,NumMomBins_pp*3,kMin,kMax*3);
  DLM_Ck Ck_Kitty_pL(Kitty_pL.GetNumSourcePars(),0,Kitty_pL,NumMomBins_pL*3,kMin,kMax*3);

  Ck_Kitty_pp.SetCutOff(300,700);
  Ck_Kitty_pL.SetCutOff(300,700);

  Ck_Kitty_pp.Update();
  Ck_Kitty_pL.Update();

  TH2F* hReso_pp = AnalysisObject.GetResolutionMatrix("pp13TeV_HM_DimiJun20","pp");
  TH2F* hFeed_pp_pL = AnalysisObject.GetResidualMatrix("pp","pLambda");

  DLM_CkDecomposition CkDecomp_pp("pp",3,Ck_Kitty_pp,hReso_pp);
  DLM_CkDecomposition CkDecomp_pL("pLambda",2,Ck_Kitty_pL,NULL);

  DLM_Histo<double> lambda_pL_to_pp;
  lambda_pL_to_pp.SetUp(1);
  lambda_pL_to_pp.SetUp(0,100,0,1000);
  lambda_pL_to_pp.Initialize();
  lambda_pL_to_pp.SetBinContentAll(0.2);
  if(OutputFileName=="Dynamic"){
    for(unsigned uBin=0; uBin<lambda_pL_to_pp.GetNbins(); uBin++){
      if(lambda_pL_to_pp.GetBinCenter(0,uBin)<40){
        lambda_pL_to_pp.SetBinContent(uBin,0.05);
      }
    }
  }

  CkDecomp_pL.AddContribution(0,0.4,DLM_CkDecomposition::cFeedDown);
  CkDecomp_pL.AddContribution(1,0.1,DLM_CkDecomposition::cFake);

  CkDecomp_pp.AddContribution(0,lambda_pL_to_pp,DLM_CkDecomposition::cFeedDown,&CkDecomp_pL,hFeed_pp_pL);
  CkDecomp_pp.AddContribution(1,0.1,DLM_CkDecomposition::cFeedDown);
  CkDecomp_pp.AddContribution(2,0.01,DLM_CkDecomposition::cFake);

  CkDecomp_pp.Update();
  CkDecomp_pL.Update();

  TFile fOutput(TString::Format("%s/OtherTasks/TestLambdaKstar/%s.root",GetFemtoOutputFolder(),OutputFileName.Data()), "recreate");
  TGraph g_cats_pp;
  g_cats_pp.SetName("g_cats_pp");
  g_cats_pp.SetLineWidth(6);
  g_cats_pp.SetLineColor(kGreen+1);
  TGraph g_ck_pp;
  g_ck_pp.SetName("g_ck_pp");
  g_ck_pp.SetLineWidth(5);
  g_ck_pp.SetLineColor(kBlue);
  g_ck_pp.SetLineStyle(2);
  TGraph g_dec_pp;
  g_dec_pp.SetName("g_dec_pp");
  g_dec_pp.SetLineWidth(4);
  g_dec_pp.SetLineColor(kRed+1);
  for(unsigned uMom=0; uMom<NumMomBins_pp*3; uMom++){
    double kstar = Ck_Kitty_pp.GetBinCenter(0, uMom);
    if(kstar<kMax) g_cats_pp.SetPoint(uMom, kstar, Kitty_pp.GetCorrFun(uMom));
    g_ck_pp.SetPoint(uMom, kstar, Ck_Kitty_pp.GetBinContent(uMom));
    g_dec_pp.SetPoint(uMom, kstar, CkDecomp_pp.EvalCk(kstar));
  }
  g_cats_pp.Write();
  g_ck_pp.Write();
  g_dec_pp.Write();
}


void Raffa_pLambda_Gauss(){
    const double Pi = TMath::Pi();
    const double fermi = 1./197.3269631;
    const double m1 = 938.2720813;
    const double m2 = 938.2720813;
    const double m3 = 1115.683;
    const double mu23 = m2*m3/(m2+m3);

    const int nBins = 100;
    const double kmin = 0.;
    const double kmax = 250.;


    const double SourceRadValue_pL = 1.19;//1.23;//1.15; // Effective Source after ERRATUM

    double V0NLO13[] = {-35.13,-30.18,-30.574,-31.851,-34.831,-37.198};
    double r0NLO13[] = {1.375,1.467,1.459,1.434,1.38,1.342};
    double V1NLO13[] = {-23.239,-29.205,-33.839,-36.258,-38.455,-39.143};
    double r1NLO13[] = {1.482,1.338,1.247,1.216,1.183,1.17};
    TString NLO13label[] = {"pLCF_NLO13_450","pLCF_NLO13_500","pLCF_NLO13_550","pLCF_NLO13_600","pLCF_NLO13_650","pLCF_NLO13_700"}; 

    double V0NLO19[] = {-25.954,-28.817,-31.851,-34.831};
    double r0NLO19[] = {1.563,1.495,1.434,1.38};
    double V1NLO19[] = {-38.984,-39.470,-42.055,-40.373};
    double r1NLO19[] = {1.178,1.163,1.126,1.143};
    TString NLO19label[] = {"pLCF_NLO19_500","pLCF_NLO19_550","pLCF_NLO19_600","pLCF_NLO19_650"}; 

    double V0SMS[] = {-31.14,-29.753,-34.273};
    double r0SMS[] = {1.439,1.466,1.382};
    double V1SMS[] = {-27.544,-28.609,-27.392};
    double r1SMS[] = {1.361,1.344,1.364};
    TString SMSlabel[] = {"pLCF_SMS_500","pLCF_SMS_550","pLCF_SMS_600"}; 

    double V0 = -29.9;
    double r0 = 1.47;
    double V1 = -20.6;
    double r1 = 1.547;

    double V00 = -25.954;
    double r00 = 1.563;
    double V10 = -38.984;
    double r10 = 1.178;



    //CATS NLO13[6];
    TGraph *hCFNLO13[6];
    TGraph *hCFNLO19[4];
    TGraph *hCFSMS[3];
    TGraph *hCFUsmani;

    for (int i = 0; i < 6; i++)
    {
        CATS NLO13;
        NLO13.SetMomBins(nBins, kmin, kmax);
        CATSparameters* cParSource = new CATSparameters(CATSparameters::tSource, 1, true);
        cParSource->SetParameter(0, SourceRadValue_pL);
        NLO13.SetAnaSource(GaussSource, *cParSource);
        NLO13.SetQ1Q2(0); // 1 same charge; 0 neutral charge; -1
        NLO13.SetQuantumStatistics(false); // set it to true if you want quantum statistics
        NLO13.SetRedMass(mu23);
        NLO13.SetNumChannels(2);
        NLO13.SetNumPW(0,1);
        NLO13.SetNumPW(1,1);
        NLO13.SetSpin(0,0);
        NLO13.SetSpin(1,1);
        NLO13.SetChannelWeight(0, 1./4.);
        NLO13.SetChannelWeight(1, 3./4.);
        CATSparameters cParPotS0Gauss(CATSparameters::tPotential, 2, true);
        cParPotS0Gauss.SetParameter(0, V0NLO13[i]);
        cParPotS0Gauss.SetParameter(1, 1./pow(r0NLO13[i],2));
        NLO13.SetShortRangePotential(0, 0, Gaussian, cParPotS0Gauss);
        CATSparameters cParPotS1Gauss(CATSparameters::tPotential, 2, true);
        cParPotS1Gauss.SetParameter(0, V1NLO13[i]);
        cParPotS1Gauss.SetParameter(1, 1./pow(r1NLO13[i],2));
        NLO13.SetShortRangePotential(1, 0, Gaussian, cParPotS1Gauss);
        NLO13.SetNotifications(CATS::nWarning);
        NLO13.SetMaxNumThreads(1);
        NLO13.KillTheCat();
        
        hCFNLO13[i] = new TGraph();
        for (int j = 0; j < nBins; j++)
        {
          hCFNLO13[i]->SetPoint(j,NLO13.GetMomentum(j),NLO13.GetCorrFun(j));
        }
    }

    for (int i = 0; i < 4; i++)
    {
        CATS NLO19;
        NLO19.SetMomBins(nBins, kmin, kmax);
        CATSparameters* cParSource = new CATSparameters(CATSparameters::tSource, 1, true);
        cParSource->SetParameter(0, SourceRadValue_pL);
        NLO19.SetAnaSource(GaussSource, *cParSource);
        NLO19.SetQ1Q2(0); // 1 same charge; 0 neutral charge; -1
        NLO19.SetQuantumStatistics(false); // set it to true if you want quantum statistics
        NLO19.SetRedMass(mu23);
        NLO19.SetNumChannels(2);
        NLO19.SetNumPW(0,1);
        NLO19.SetNumPW(1,1);
        NLO19.SetSpin(0,0);
        NLO19.SetSpin(1,1);
        NLO19.SetChannelWeight(0, 1./4.);
        NLO19.SetChannelWeight(1, 3./4.);
        CATSparameters cParPotS0Gauss(CATSparameters::tPotential, 2, true);
        cParPotS0Gauss.SetParameter(0, V0NLO19[i]);
        cParPotS0Gauss.SetParameter(1, 1./pow(r0NLO19[i],2));
        NLO19.SetShortRangePotential(0, 0, Gaussian, cParPotS0Gauss);
        CATSparameters cParPotS1Gauss(CATSparameters::tPotential, 2, true);
        cParPotS1Gauss.SetParameter(0, V1NLO19[i]);
        cParPotS1Gauss.SetParameter(1, 1./pow(r1NLO19[i],2));
        NLO19.SetShortRangePotential(1, 0, Gaussian, cParPotS1Gauss);
        NLO19.SetNotifications(CATS::nWarning);
        NLO19.SetMaxNumThreads(1);
        NLO19.KillTheCat();
        
        hCFNLO19[i] = new TGraph();
        for (int j = 0; j < nBins; j++)
        {
          hCFNLO19[i]->SetPoint(j,NLO19.GetMomentum(j),NLO19.GetCorrFun(j));
        }
    }

    for (int i = 0; i < 3; i++)
    {
        CATS SMS;
        SMS.SetMomBins(nBins, kmin, kmax);
        CATSparameters* cParSource = new CATSparameters(CATSparameters::tSource, 1, true);
        cParSource->SetParameter(0, SourceRadValue_pL);
        SMS.SetAnaSource(GaussSource, *cParSource);
        SMS.SetQ1Q2(0); // 1 same charge; 0 neutral charge; -1
        SMS.SetQuantumStatistics(false); // set it to true if you want quantum statistics
        SMS.SetRedMass(mu23);
        SMS.SetNumChannels(2);
        SMS.SetNumPW(0,1);
        SMS.SetNumPW(1,1);
        SMS.SetSpin(0,0);
        SMS.SetSpin(1,1);
        SMS.SetChannelWeight(0, 1./4.);
        SMS.SetChannelWeight(1, 3./4.);
        CATSparameters cParPotS0Gauss(CATSparameters::tPotential, 2, true);
        cParPotS0Gauss.SetParameter(0, V0SMS[i]);
        cParPotS0Gauss.SetParameter(1, 1./pow(r0SMS[i],2));
        SMS.SetShortRangePotential(0, 0, Gaussian, cParPotS0Gauss);
        CATSparameters cParPotS1Gauss(CATSparameters::tPotential, 2, true);
        cParPotS1Gauss.SetParameter(0, V1SMS[i]);
        cParPotS1Gauss.SetParameter(1, 1./pow(r1SMS[i],2));
        SMS.SetShortRangePotential(1, 0, Gaussian, cParPotS1Gauss);
        SMS.SetNotifications(CATS::nWarning);
        SMS.SetMaxNumThreads(1);
        SMS.KillTheCat();
    
        double ScatLen, EffRan;
        TH1F *hFit;
        TF1 *fitSP;
        GetScattParameters(SMS, ScatLen, EffRan, hFit, fitSP, 3, false, false, 0);
        cout<< i <<" Spin 0 "<<ScatLen<<"    "<<EffRan<<endl;
        delete hFit;
        delete fitSP;

        GetScattParameters(SMS, ScatLen, EffRan, hFit, fitSP, 3, false, false, 1);
        cout<< i<<" Spin 1 "<<ScatLen<<"    "<<EffRan<<endl;
        delete hFit;
        delete fitSP;
     
        hCFSMS[i] = new TGraph();
        for (int j = 0; j < nBins; j++)
        {
          hCFSMS[i]->SetPoint(j,SMS.GetMomentum(j),SMS.GetCorrFun(j));
        }
    }
    CATS Usmani;
    Usmani.SetMomBins(nBins, kmin, kmax);
    CATSparameters* cParSource = new CATSparameters(CATSparameters::tSource, 1, true);
    cParSource->SetParameter(0, SourceRadValue_pL);
    Usmani.SetAnaSource(GaussSource, *cParSource);
    Usmani.SetQ1Q2(0); // 1 same charge; 0 neutral charge; -1
    Usmani.SetQuantumStatistics(false); // set it to true if you want quantum statistics
    Usmani.SetRedMass(mu23);
    Usmani.SetNumChannels(2);
    Usmani.SetNumPW(0,1);
    Usmani.SetNumPW(1,1);
    Usmani.SetSpin(0,0);
    Usmani.SetSpin(1,1);
    Usmani.SetChannelWeight(0, 1./4.);
    Usmani.SetChannelWeight(1, 3./4.);
    CATSparameters cParPotS0Gauss(CATSparameters::tPotential, 2, true);
    cParPotS0Gauss.SetParameter(0, V0);
    cParPotS0Gauss.SetParameter(1, 1./pow(r0,2));
    Usmani.SetShortRangePotential(0, 0, Gaussian, cParPotS0Gauss);
    CATSparameters cParPotS1Gauss(CATSparameters::tPotential, 2, true);
    cParPotS1Gauss.SetParameter(0, V1);
    cParPotS1Gauss.SetParameter(1, 1./pow(r1,2));
    Usmani.SetShortRangePotential(1, 0, Gaussian, cParPotS1Gauss);
    Usmani.SetNotifications(CATS::nWarning);
    Usmani.SetMaxNumThreads(1);
    Usmani.KillTheCat();
    
    
    hCFUsmani = new TGraph();
    for (int j = 0; j < nBins; j++)
    {
        hCFUsmani->SetPoint(j,Usmani.GetMomentum(j),Usmani.GetCorrFun(j));
    }


    CATS NLO190;
    NLO190.SetMomBins(nBins, kmin, kmax);
    NLO190.SetAnaSource(GaussSource, *cParSource);
    NLO190.SetQ1Q2(0); // 1 same charge; 0 neutral charge; -1
    NLO190.SetQuantumStatistics(false); // set it to true if you want quantum statistics
    NLO190.SetRedMass(mu23);
    NLO190.SetNumChannels(2);
    NLO190.SetNumPW(0,1);
    NLO190.SetNumPW(1,1);
    NLO190.SetSpin(0,0);
    NLO190.SetSpin(1,1);
    NLO190.SetChannelWeight(0, 1./4.);
    NLO190.SetChannelWeight(1, 3./4.);
    cParPotS0Gauss.SetParameter(0, V00);
    cParPotS0Gauss.SetParameter(1, 1./pow(r00,2));
    NLO190.SetShortRangePotential(0, 0, Gaussian, cParPotS0Gauss);
    cParPotS1Gauss.SetParameter(0, V10);
    cParPotS1Gauss.SetParameter(1, 1./pow(r10,2));
    NLO190.SetShortRangePotential(1, 0, Gaussian, cParPotS1Gauss);
    NLO190.SetNotifications(CATS::nWarning);
    NLO190.SetMaxNumThreads(1);
    NLO190.KillTheCat();
    
    for (int j = 0; j < nBins; j++)
    {
        hCFNLO19[0]->SetPoint(j,NLO190.GetMomentum(j),NLO190.GetCorrFun(j));
    }


    TFile *OutFile = new TFile("OutRaf.root","recreate");
    OutFile->cd();
        for (int i = 0; i < 6; i++)
        {
            hCFNLO13[i]->Write(Form("%s",NLO13label[i].Data()));
        }
        for (int i = 0; i < 4; i++)
        {
            hCFNLO19[i]->Write(Form("%s",NLO19label[i].Data()));
        }
        for (int i = 0; i < 3; i++)
        {
            hCFSMS[i]->Write(Form("%s",SMSlabel[i].Data()));
        }
        hCFUsmani->Write("pLCF_Usmani");
        
    OutFile->Close();

}


void Quick_pXi(){
    
    double kMin = 0;
    double kMax = 10;
    unsigned kSteps = 1;

    const unsigned NumPotentials = 71;


    DLM_CommonAnaFunctions AnalysisObject;
    AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles/",GetCernBoxDimi()));

    CATS Kitty;
    Kitty.SetMomBins(kSteps, kMin, kMax);
    AnalysisObject.SetUpCats_pXim(Kitty,"pXim_HALQCDPaper2020","Gauss",0, 0);
    Kitty.SetAnaSource(0,1.2);
    Kitty.SetNotifications(CATS::nWarning);
    Kitty.SetEpsilonProp(1e-10);
    Kitty.SetEpsilonConv(1e-10);
    //Kitty.SetQ1Q2(0);

    double min_val = 1000;
    double max_val = 0;
    unsigned min_sol = 1000;
    unsigned max_sol = 1000;

    for(unsigned uPot=0; uPot<NumPotentials; uPot++){
        Kitty.SetShortRangePotential(0, 0, 1, uPot);
        Kitty.SetShortRangePotential(1, 0, 1, uPot);
        Kitty.SetShortRangePotential(2, 0, 1, uPot);
        Kitty.SetShortRangePotential(3, 0, 1, uPot);
        Kitty.KillTheCat();
        if(Kitty.GetCorrFun(0) < min_val){
            min_val = Kitty.GetCorrFun(0);
            min_sol = uPot;
        }
        if(Kitty.GetCorrFun(0) > max_val){
            max_val = Kitty.GetCorrFun(0);
            max_sol = uPot;
        }        
        printf("uPot = %u; C(5) = %.2f\n", uPot, Kitty.GetCorrFun(0));
    }

    printf("MIN: %.2f @ %u\n", min_val, min_sol);
    printf("MAX: %.2f @ %u\n", max_val, max_sol);

}


void pp_reff(double r_core){


    DLM_CleverMcLevyResoTM MagicSource;

    //0 = correct sign
    SetUp_RSM_pp(MagicSource, TString::Format("%s",GetCernBoxDimi()).Data(), 0);

    double r_eff = GetReff(MagicSource, r_core);

    printf("r_core = %.3f --> r_eff = %.3f\n", r_core, r_eff);

}


void pSp_test1(){

    const unsigned NumMomBins = 100;
    const double kMin = 0;
    const double kMax = 400;


    CATSparameters cParSource = CATSparameters(CATSparameters::tSource, 1, true);
    cParSource.SetParameter(0, 1.2);

    CATS Cat_Old;
    Cat_Old.SetMomBins(NumMomBins,kMin,kMax);
    Cat_Old.SetAnaSource(GaussSource, cParSource);


    CATS Cat_New;
    Cat_New.SetMomBins(NumMomBins,kMin,kMax);
    Cat_New.SetAnaSource(GaussSource, cParSource);

    DLM_Histo<complex<double>> ***ExternalWF_Old = Init_pSigmaPlus_Haidenbauer(TString::Format("%s/CatsFiles/Interaction/Haidenbauer/pSigmaPlus/",GetCernBoxDimi()).Data(), Cat_Old, 0);
    DLM_Histo<complex<double>> ***ExternalWF_New = Init_pSigmaPlus_Haidenbauer(TString::Format("%s/CatsFiles/Interaction/Haidenbauer/pSigmaPlus/",GetCernBoxDimi()).Data(), Cat_New, 1);

    Cat_Old.SetExternalWaveFunction(0, 0, ExternalWF_Old[0][0][0], ExternalWF_Old[1][0][0]);
    Cat_Old.SetExternalWaveFunction(1, 0, ExternalWF_Old[0][1][0], ExternalWF_Old[1][1][0]);

    Cat_New.SetExternalWaveFunction(0, 0, ExternalWF_New[0][0][0], ExternalWF_New[1][0][0]);
    Cat_New.SetExternalWaveFunction(1, 0, ExternalWF_New[0][1][0], ExternalWF_New[1][1][0]);

    Cat_Old.KillTheCat();
    Cat_New.KillTheCat();

    TGraph gCat_Old;
    gCat_Old.SetName("gCat_Old");
    TGraph gCat_New;
    gCat_New.SetName("gCat_New");

    for(unsigned uMom=0; uMom<NumMomBins; uMom++){
        double kstar = Cat_Old.GetMomentum(uMom);
        gCat_Old.SetPoint(uMom, kstar, Cat_Old.GetCorrFun(uMom));
        gCat_New.SetPoint(uMom, kstar, Cat_New.GetCorrFun(uMom));
    }

    TFile fOutput(TString::Format("%s/OtherTasks/pSp_test1.root",GetFemtoOutputFolder()), "recreate");
    gCat_Old.Write();
    gCat_New.Write();

}


void pSp_test2(){

    const unsigned NumMomBins = 100;
    const double kMin = 0;
    const double kMax = 400;


    CATSparameters cParSource = CATSparameters(CATSparameters::tSource, 1, true);
    cParSource.SetParameter(0, 1.2);//source size, assuming we use a Gauss

    CATS Cat;
    Cat.SetMomBins(NumMomBins,kMin,kMax);
    Cat.SetAnaSource(GaussSource, cParSource);

    DLM_Histo<complex<double>> ***ExternalWF = Init_pSigmaPlus_Haidenbauer(TString::Format("%s/CatsFiles/Interaction/Haidenbauer/pSigmaPlus/",GetCernBoxDimi()).Data(), Cat, 1);

    Cat.SetExternalWaveFunction(0, 0, ExternalWF[0][0][0], ExternalWF[1][0][0]);
    Cat.SetExternalWaveFunction(1, 0, ExternalWF[0][1][0], ExternalWF[1][1][0]);

    Cat.KillTheCat();

    TGraph gCat;
    gCat.SetName("gCat");

    for(unsigned uMom=0; uMom<NumMomBins; uMom++){
        double kstar = Cat.GetMomentum(uMom);
        gCat.SetPoint(uMom, kstar, Cat.GetCorrFun(uMom));
    }

    TFile fOutput(TString::Format("%s/OtherTasks/pSp_test2.root",GetFemtoOutputFolder()), "recreate");
    gCat.Write();

}

void pp_large_RSM_source(){
    double kMin = 0;
    double kMax = 300;
    unsigned NumMomBins = 75;

    double rMin = 0;
    double rMax = 64;
    unsigned NumRadBins = 4096;

    std::vector<float> CoreSize = {1,2,2.5,3.0,3.5,4};

    TH1F** hSource = new TH1F*[CoreSize.size()];
    for(unsigned uSrc=0; uSrc<CoreSize.size(); uSrc++){
        hSource[uSrc] = new TH1F(
            TString::Format("hSrc_%.2f",CoreSize.at(uSrc)), 
            TString::Format("hSrc_%.2f",CoreSize.at(uSrc)), 
            NumRadBins, rMin, rMax);
    }

    CATS Kitty;
    Kitty.SetMomBins(NumMomBins, kMin, kMax);

    DLM_CommonAnaFunctions AnalysisObject;
    AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles/",GetCernBoxDimi()));
    AnalysisObject.SetUpCats_pp(Kitty,"AV18","McLevy_ResoTM",0, 202);
    for(unsigned uSrc=0; uSrc<CoreSize.size(); uSrc++){
        //AnalysisObject->GetCleverMcLevyResoTM_pp()
        Kitty.SetAnaSource(0, CoreSize.at(uSrc));
        Kitty.KillTheCat();

        for(unsigned uRad=0; uRad<NumRadBins; uRad++){
            double rad = hSource[uSrc]->GetBinCenter(uRad+1);
            double src = Kitty.EvaluateTheSource(1, rad, 0);
            hSource[uSrc]->SetBinContent(uRad+1, src);
        }
    }
    
    TFile fOutput(TString::Format("%s/OtherTasks/pp_large_RSM_source.root",GetFemtoOutputFolder()), "recreate");
    for(unsigned uSrc=0; uSrc<CoreSize.size(); uSrc++){
        hSource[uSrc]->Write();
    }
}


void pd_withQS(){
    unsigned NumMomBins = 40;
    double kMin = 0;
    double kMax = 400;

    std::vector<float> source_size = {1.0/sqrt(2), 1.08/sqrt(2), 1.43/sqrt(2), 2.0/sqrt(2), 3.0/sqrt(2), 5.0/sqrt(2), 5, 7};
    std::vector<int> colors = {kBlack, kCyan+1, kGreen+1, kBlue, kRed, kYellow+1, kRed+2, kBlue-2};



    CATS Kitty;
    Kitty.SetMomBins(NumMomBins, kMin, kMax);
    CATSparameters cParSource = CATSparameters(CATSparameters::tSource, 1, true);
    //cParSource->SetParameter(0, SourceRadValue_pL);


    Kitty.SetAnaSource(GaussSource, cParSource);
    Kitty.SetQ1Q2(1); // 1 same charge; 0 neutral charge; -1
    Kitty.SetQuantumStatistics(true); // set it to true if you want quantum statistics
    Kitty.SetRedMass((Mass_p * Mass_d)/(Mass_p + Mass_d));
    Kitty.SetNumChannels(2);
    Kitty.SetNumPW(0,1);
    Kitty.SetNumPW(1,1);
    Kitty.SetSpin(0,0);
    Kitty.SetSpin(1,1);
    Kitty.SetChannelWeight(0, 1./3.);
    Kitty.SetChannelWeight(1, 2./3.);
    //Kitty.SetNotifications(CATS::nWarning);
    //Usmani.SetMaxNumThreads(1);



    CATS KittyC;
    KittyC.SetMomBins(NumMomBins, kMin, kMax);
    CATSparameters cParSourceC = CATSparameters(CATSparameters::tSource, 1, true);
    //cParSource->SetParameter(0, SourceRadValue_pL);


    KittyC.SetAnaSource(GaussSource, cParSourceC);
    KittyC.SetQ1Q2(1); // 1 same charge; 0 neutral charge; -1
    KittyC.SetQuantumStatistics(false); // set it to true if you want quantum statistics
    KittyC.SetRedMass((Mass_p * Mass_d)/(Mass_p + Mass_d));
    KittyC.SetNumChannels(2);
    KittyC.SetNumPW(0,1);
    KittyC.SetNumPW(1,1);
    KittyC.SetSpin(0,0);
    KittyC.SetSpin(1,1);
    KittyC.SetChannelWeight(0, 1./4.);
    KittyC.SetChannelWeight(1, 3./4.);


    CATS KittyQS;
    KittyQS.SetMomBins(NumMomBins, kMin, kMax);
    CATSparameters cParSourceQS = CATSparameters(CATSparameters::tSource, 1, true);
    //cParSource->SetParameter(0, SourceRadValue_pL);


    KittyQS.SetAnaSource(GaussSource, cParSourceQS);
    KittyQS.SetQ1Q2(0); // 1 same charge; 0 neutral charge; -1
    KittyQS.SetQuantumStatistics(true); // set it to true if you want quantum statistics
    KittyQS.SetRedMass((Mass_p * Mass_d)/(Mass_p + Mass_d));
    KittyQS.SetNumChannels(2);
    KittyQS.SetNumPW(0,1);
    KittyQS.SetNumPW(1,1);
    KittyQS.SetSpin(0,0);
    KittyQS.SetSpin(1,1);
    KittyQS.SetChannelWeight(0, 1./3.);
    KittyQS.SetChannelWeight(1, 2./3.);


    TGraph* gCk_CQS = new TGraph[source_size.size()];
    TGraph* gCk_C = new TGraph[source_size.size()];
    TGraph* gCk_QS = new TGraph[source_size.size()];

    for(unsigned uSrc=0; uSrc<source_size.size(); uSrc++){
        double SrcSize = source_size.at(uSrc);
        gCk_CQS[uSrc].SetName(TString::Format("gCk_CQS_Src%.2fm", SrcSize));
        gCk_CQS[uSrc].SetLineWidth(4);
        gCk_CQS[uSrc].SetLineColor(colors.at(uSrc));
        gCk_CQS[uSrc].SetMarkerStyle(8);
        gCk_CQS[uSrc].SetMarkerColor(colors.at(uSrc));
        gCk_CQS[uSrc].SetMarkerSize(1);

        gCk_QS[uSrc].SetName(TString::Format("gCk_QS_Src%.2fm", SrcSize));
        gCk_QS[uSrc].SetLineWidth(4);
        gCk_QS[uSrc].SetLineColor(colors.at(uSrc));
        gCk_QS[uSrc].SetLineStyle(3);
        gCk_QS[uSrc].SetMarkerStyle(8);
        gCk_QS[uSrc].SetMarkerColor(colors.at(uSrc));
        gCk_QS[uSrc].SetMarkerSize(1);

        gCk_C[uSrc].SetName(TString::Format("gCk_C_Src%.2fm", SrcSize));
        gCk_C[uSrc].SetLineWidth(4);
        gCk_C[uSrc].SetLineColor(colors.at(uSrc));
        gCk_C[uSrc].SetLineStyle(2);
        gCk_C[uSrc].SetMarkerStyle(4);
        gCk_C[uSrc].SetMarkerColor(colors.at(uSrc));
        gCk_C[uSrc].SetMarkerSize(1);

        Kitty.SetAnaSource(0, SrcSize);
        Kitty.KillTheCat();
        for(unsigned uMom=0; uMom<NumMomBins; uMom++){
            double MOM = Kitty.GetMomentum(uMom);
            gCk_CQS[uSrc].SetPoint(uMom, MOM, Kitty.GetCorrFun(uMom));
        }

        KittyC.SetAnaSource(0, SrcSize);
        KittyC.KillTheCat(CATS::kAllChanged);
        for(unsigned uMom=0; uMom<NumMomBins; uMom++){
        double MOM = KittyC.GetMomentum(uMom);
        gCk_C[uSrc].SetPoint(uMom, MOM, KittyC.GetCorrFun(uMom));
        }

        KittyQS.SetAnaSource(0, SrcSize);
        KittyQS.KillTheCat();
        for(unsigned uMom=0; uMom<NumMomBins; uMom++){
            double MOM = Kitty.GetMomentum(uMom);
            gCk_QS[uSrc].SetPoint(uMom, MOM, KittyQS.GetCorrFun(uMom));
        }
    }

    TFile fOutput(TString::Format("%s/OtherTasks/pd_withQS.root",GetFemtoOutputFolder()), "recreate");
    for(unsigned uSrc=0; uSrc<source_size.size(); uSrc++){
        gCk_C[uSrc].Write();
    }
    for(unsigned uSrc=0; uSrc<source_size.size(); uSrc++){
        gCk_QS[uSrc].Write();
    }   
    for(unsigned uSrc=0; uSrc<source_size.size(); uSrc++){
        gCk_CQS[uSrc].Write();
    }

    delete [] gCk_QS;
    delete [] gCk_CQS;
    delete [] gCk_C;


}



//DummyBootTest_pSigma_inspired, singlet
//two gaussians, one dummy to bring it all to zero (fixed pars) and 
//one extra. Amplitude width and mean a linked to mimic femto
double DBST_PS_INSP_sin(double* x, double* pars){
    double& rad = *x;
    double zero_gauss = 1.-exp(-pow(-(rad)/(10.),2.));
    double femto_gauss = pars[0]/15.*exp(-pow(-(rad-pars[0]*1.3)/(pars[0]),2.))+1;
    return zero_gauss*femto_gauss;
}

//similar
double DBST_PS_INSP_tri(double* x, double* pars){
    double& rad = *x;
    double zero_gauss = 1.-exp(-pow(-(rad)/(10.),2.));
    double femto_exp = -exp(-(rad+50.)/(pars[0]))+1;
    //double femto_gauss = -pars[0]/360.*exp(-pow(-(rad)/(pars[0]),2.))+1;
    return zero_gauss*femto_exp;
}

double DBST_PS_INSP_tot(double* x, double* pars){
    return 0.25*DBST_PS_INSP_sin(x, pars) + 0.75*DBST_PS_INSP_tri(x, &pars[1]);
}

//generate Ck as two gaussians, which look a bit like the two channels.
//generate dummy data
//with the bootstrap, see which delta chi2 you get for 1 sigma
void DummyBootTest_pSigma_inspired(){

    double kMin = 20;
    double kMax = 420;
    unsigned NumMomBins = 10;

    //40,80,120,160,200,240,280,320,360,400
    std::vector<double> bin_center = {55,90,125,160,200,240,280,320,360,400};
    double MOM;
    double Ck_s;
    double Ck_t;
    double Ck;
    double Ck_s_rnd;
    double Ck_t_rnd;   
    double Ck_rnd;

    TF1* fSin = new TF1("fSin",DBST_PS_INSP_sin,kMin,kMax,1);
    fSin->SetParameter(0, 40);

    TF1* fTri = new TF1("fTri",DBST_PS_INSP_tri,kMin,kMax,1);
    fTri->SetParameter(0, 120);    

    TF1* fTot = new TF1("fTot",DBST_PS_INSP_tot,kMin,kMax,2);
    fTot->SetParameter(0, 40);
    fTot->SetParameter(1, 120);

    TH2F* h_fit_result = new TH2F("h_fit_result", "h_fit_result", 256, 0, 80, 256, 0, 240);

    TH2F* h_fit_result_v1p0 = (TH2F*)h_fit_result->Clone("h_fit_result_v1p0");
    TH2F* h_fit_result_v2p3 = (TH2F*)h_fit_result->Clone("h_fit_result_v2p3");


    TH2F* h_fit_result_dlt1p0 = (TH2F*)h_fit_result->Clone("h_fit_result_dlt1p0");
    TH2F* h_fit_result_dlt2p3 = (TH2F*)h_fit_result->Clone("h_fit_result_dlt2p3");

    TH2F* h_fit_result_hc1 = (TH2F*)h_fit_result->Clone("h_fit_result_hc1");
    TH2F* h_fit_result_hc2 = (TH2F*)h_fit_result->Clone("h_fit_result_hc2");   


    TRandom3 rangen(23);

    TGraphErrors* gData = new TGraphErrors();
    gData->SetName("gData");
    for(unsigned uMom=0; uMom<NumMomBins; uMom++){
        MOM = bin_center.at(uMom);
        Ck_s = fSin->Eval(MOM);
        Ck_t = fTri->Eval(MOM);
        gData->SetPoint(uMom, MOM, 0.25*Ck_s + 0.75*Ck_t);
    }
    gData->SetPointError(0,0,0.1);
    gData->SetPointError(1,0,0.06);
    gData->SetPointError(2,0,0.05);
    gData->SetPointError(3,0,0.04);
    gData->SetPointError(4,0,0.035);
    gData->SetPointError(5,0,0.03);
    gData->SetPointError(6,0,0.025);
    gData->SetPointError(7,0,0.02);
    gData->SetPointError(8,0,0.015);
    gData->SetPointError(9,0,0.01);

    for(unsigned uMom=0; uMom<NumMomBins; uMom++){
        MOM = bin_center.at(uMom);
        Ck_s = fSin->Eval(MOM);
        Ck_t = fTri->Eval(MOM);
        Ck = 0.25*Ck_s+0.75*Ck_t;
        Ck_rnd = -1; while(Ck_rnd<=0) Ck_rnd = rangen.Gaus(Ck, gData->GetErrorY(uMom));
        gData->SetPoint(uMom, MOM, Ck_rnd);
    }


    //fit the generated data to get the baseline parameters
    gData->Fit(fTot,"S, N, R, M");
    double Base_Chi2 = fTot->GetChisquare();
//Base_Chi2 = NumMomBins-2;
    fSin->SetParameter(0, fTot->GetParameter(0));
    fTri->SetParameter(0, fTot->GetParameter(1));

    TGraphErrors* gBoot = (TGraphErrors*)gData->Clone("gBoot");

    unsigned NumBootIter = 100*1000;
    double BestChi2 = 1e6;

    for(unsigned uBoot=0; uBoot<NumBootIter; uBoot++){
        for(unsigned uMom=0; uMom<NumMomBins; uMom++){
            MOM = bin_center.at(uMom);
            Ck_s = fSin->Eval(MOM);
            Ck_t = fTri->Eval(MOM);
            Ck = 0.25*Ck_s+0.75*Ck_t;
            Ck_rnd = -1; while(Ck_rnd<=0) Ck_rnd = rangen.Gaus(Ck, gData->GetErrorY(uMom));
            gBoot->SetPoint(uMom, MOM, Ck_rnd);
        }
        gBoot->Fit(fTot,"Q, S, N, R, M");
        double Chi2_fit = fTot->GetChisquare();
        double Chi2_dlm=0;
        for(unsigned uMom=0; uMom<NumMomBins; uMom++){
            MOM = bin_center.at(uMom);
            Chi2_dlm += pow((fTot->Eval(MOM)-gBoot->Eval(MOM))/gBoot->GetErrorY(uMom),2.);
        }
        //printf("chi2 = %f %f\n", Chi2_fit, Chi2_dlm);
        h_fit_result->Fill(fTot->GetParameter(0), fTot->GetParameter(1));
        if(Chi2_dlm<1) h_fit_result_v1p0->Fill(fTot->GetParameter(0), fTot->GetParameter(1));
        if(Chi2_dlm<2.3) h_fit_result_v2p3->Fill(fTot->GetParameter(0), fTot->GetParameter(1));

        if(Chi2_dlm<1+Base_Chi2) h_fit_result_dlt1p0->Fill(fTot->GetParameter(0), fTot->GetParameter(1));
        if(Chi2_dlm<2.3+Base_Chi2) h_fit_result_dlt2p3->Fill(fTot->GetParameter(0), fTot->GetParameter(1));


        if(Chi2_dlm<BestChi2) BestChi2 = Chi2_dlm;

    }
    printf("BestChi2 = %f\n",BestChi2);


    double MaxAllowedChi2 = 0;
    for(unsigned uBinX=0; uBinX<h_fit_result->GetXaxis()->GetNbins(); uBinX++){
        for(unsigned uBinY=0; uBinY<h_fit_result->GetYaxis()->GetNbins(); uBinY++){
            if(float(h_fit_result->GetBinContent(uBinX+1,uBinY+1))>float(NumBootIter)*0.00107){
                //printf("%f\n",h_fit_result->GetBinContent(uBinX+1,uBinY+1));
                for(int iEntr=0; iEntr<h_fit_result->GetBinContent(uBinX+1,uBinY+1); iEntr++){
                    h_fit_result_hc1->Fill(h_fit_result->GetXaxis()->GetBinCenter(uBinX+1),h_fit_result->GetYaxis()->GetBinCenter(uBinY+1));
                }
                fTot->FixParameter(0, h_fit_result->GetXaxis()->GetBinCenter(uBinX+1));
                fTot->FixParameter(1, h_fit_result->GetYaxis()->GetBinCenter(uBinY+1));
                gData->Fit(fTot,"Q, S, N, R, M");
                if(MaxAllowedChi2<fTot->GetChisquare()) MaxAllowedChi2 = fTot->GetChisquare();
            }
            if(h_fit_result->GetBinContent(uBinX+1,uBinY+1)>float(NumBootIter)*0.0005)
                for(int iEntr=0; iEntr<h_fit_result->GetBinContent(uBinX+1,uBinY+1); iEntr++)
                    h_fit_result_hc2->Fill(h_fit_result->GetXaxis()->GetBinCenter(uBinX+1),h_fit_result->GetYaxis()->GetBinCenter(uBinY+1));
        }
    }
    printf("chi2: %f %f, dchi2 = %f\n", MaxAllowedChi2, Base_Chi2, MaxAllowedChi2 - Base_Chi2);


    TFile fOutput(TString::Format("%s/OtherTasks/DummyBootTest_pSigma_inspired.root",GetFemtoOutputFolder()), "recreate");
    gData->Write();
    h_fit_result->Write();
    h_fit_result_v1p0->Write();
    h_fit_result_v2p3->Write();
    h_fit_result_dlt1p0->Write();
    h_fit_result_dlt2p3->Write();   
    h_fit_result_hc1->Write();   
    h_fit_result_hc2->Write();   
}


void xCheck_Benedict_pSigma_v1(){
    TString InputFileName = TString::Format("%s/p_Sigma/Benedict/pSigma_Variation_Distribution.root",GetCernBoxDimi());
    TFile InFile(InputFileName, "read");
    //InFile.ls();
    TH2F* h2D_iter = (TH2F*)InFile.Get("var_dist");
    printf("Total integral = %f\n", h2D_iter->Integral());
    double missing_integral = 100. - h2D_iter->Integral();
    double my_integral = 0;

    unsigned NumElements = h2D_iter->GetXaxis()->GetNbins()*h2D_iter->GetYaxis()->GetNbins();
    float* Elements = new float[NumElements];
    for(unsigned uBinX=0; uBinX<h2D_iter->GetXaxis()->GetNbins(); uBinX++){
        for(unsigned uBinY=0; uBinY<h2D_iter->GetYaxis()->GetNbins(); uBinY++){
            Elements[uBinX*h2D_iter->GetYaxis()->GetNbins() + uBinY] = h2D_iter->GetBinContent(uBinX+1, uBinY+1);
            my_integral += h2D_iter->GetBinContent(uBinX+1, uBinY+1);
        }
    }
    printf("My integral = %f\n", my_integral);
    DLM_Sort < float, unsigned > SortTool;
    SortTool.SetData(Elements,NumElements);
    SortTool.MergeSort(true);
    SortTool.GetSortedData(Elements,Elements);
    float cum_integral = 0;
    float cut_off = 0;
    for(unsigned uEl=0; uEl<NumElements; uEl++){
        //printf("%u: %f\n",uEl,Elements[uEl]);
        //if(uEl>200) break;
        cum_integral += Elements[uEl];
        if(cum_integral>68.2){
            cut_off = Elements[uEl];
            break;
        }
    }

    TH2F* h1Sig = (TH2F*)h2D_iter->Clone("h1Sig");
    for(unsigned uBinX=0; uBinX<h2D_iter->GetXaxis()->GetNbins(); uBinX++){
        for(unsigned uBinY=0; uBinY<h2D_iter->GetYaxis()->GetNbins(); uBinY++){
            Elements[uBinX*h2D_iter->GetYaxis()->GetNbins() + uBinY] = h2D_iter->GetBinContent(uBinX+1, uBinY+1);
            my_integral += h2D_iter->GetBinContent(uBinX+1, uBinY+1);
            if(h2D_iter->GetBinContent(uBinX+1, uBinY+1)>cut_off){
                h1Sig->SetBinContent(uBinX+1, uBinY+1, h2D_iter->GetBinContent(uBinX+1, uBinY+1));
            }
            else{
                h1Sig->SetBinContent(uBinX+1, uBinY+1, 0);
            }
        }
    }    

    TFile fOutput(TString::Format("%s/OtherTasks/xCheck_Benedict_pSigma_v1.root",GetFemtoOutputFolder()), "recreate");
    h2D_iter->Write();
    h1Sig->Write();
}


void extrap_me_for_3b_v1(){

  //ppp
  TString InputFileName = TString::Format("%s/3body/LauraPrelim/Before_QM_ME/ppp_syst_defaultLarge.root",GetCernBoxDimi());
  TString InputHistName = "ME_partPlusAntiPart";
  //which bins to is (ROOT numbering) for the fit
  int first_bin_id = 0;
  int last_bin_id = 12;
  TFile fOutput(TString::Format("%s/3body/LauraPrelim/Before_QM_ME/fit_result_ppp.root",GetCernBoxDimi()), "recreate");
  
  //ppL
  //TString InputFileName = TString::Format("%s/3body/LauraPrelim/Before_QM_ME/ppL_syst_defaultLarge_MakeBigger.root",GetCernBoxDimi());
  //TString InputHistName = "ME_partPlusAntiPart";
  //which bins to is (ROOT numbering) for the fit
  //int first_bin_id = 2;
  //int last_bin_id = 20;  
  //TFile fOutput(TString::Format("%s/3body/LauraPrelim/Before_QM_ME/fit_result_ppL.root",GetCernBoxDimi()), "recreate");

  TFile fInput(InputFileName, "read");
  TH1F* hME = (TH1F*)fInput.Get(InputHistName);


  unsigned NumBootIter = 1000;
  float q_min = hME->GetXaxis()->GetBinLowEdge(first_bin_id);
  float q_max = hME->GetXaxis()->GetBinUpEdge(last_bin_id);

  //TF1* fBGHetto = new TF1("fBGHetto", "[0]*x*x*exp(-pow(x/[1],[2]))",q_min,q_max*10);
  TF1* fBGHetto = new TF1("fBGHetto", "[0]*x*x/([1]+[2]*exp(pow(x-[3],[4])))",q_min,q_max);

  fBGHetto->SetParameter(0,1e3);
  fBGHetto->FixParameter(1,1);
  //fBGHetto->SetParLimits(1,0.05,40);
  fBGHetto->SetParameter(2,1);
  //fBGHetto->SetParLimits(2,0.5,40);
  fBGHetto->SetParameter(3,3);
  fBGHetto->FixParameter(4,1);
  //fBGHetto->SetParLimits(4,0.5,2);
  
  hME->Fit(fBGHetto,"S, N, R, M");

  fOutput.cd();
  hME->Write();
  fBGHetto->Write();

  delete fBGHetto;
}



double Emma_potential(double* ppars){
// ppars[0] - radius in fm
// ppars[1] - k* in MeV
// ppars[2] - scaling constant g
// ppars[3] - scaling constant a
// m - mass of mediating particle in MeV: neutral pion

double& r = ppars[0];
double& beta = ppars[2];
double& a1 = ppars[3];
double& b1 = ppars[4];//fm
double& a2 = ppars[5];
double& b2 = ppars[6];//fm
double& a3mpi4 = ppars[7];
double& b3 = ppars[8];//fm
double& mpi = ppars[9];//MeV

double fATFF = pow(1-exp(-pow(r/b3,2.)),2.);

double rad = r;
double result;
//if(r<0.01){
//  rad = 0.01;
//}

//if(r>16){
//  result = 0;
//}
//else{
  result = beta*(a1*exp(-pow(rad/b1,2.))+a2*exp(-pow(rad/b2,2.)))+a3mpi4*fATFF*exp(-2.*mpi*rad/hbarc)/pow(rad,2.);
//}


return result;
}

void Ghetto_pPhi_check_20250426(){
  
  const double source_size = 1.08;

  CATS Kitty;
	CATS smelly_cat;
  double kMin = 0;
  double kMax = 320;
  unsigned NumMomBins = 80;
	smelly_cat.SetMomBins(NumMomBins, kMin, kMax);
	smelly_cat.SetThetaDependentSource(false);
	CATSparameters ppars(CATSparameters::tSource, 1, true);

	smelly_cat.SetAnaSource(GaussSource, ppars);
	smelly_cat.SetAnaSource(0, source_size);
	smelly_cat.SetAutoNormSource(true);
	smelly_cat.SetUseAnalyticSource(true);
	smelly_cat.SetMomentumDependentSource(false);
	smelly_cat.SetExcludeFailedBins(false);
	smelly_cat.SetQ1Q2(0);
	smelly_cat.SetQuantumStatistics(false);
	smelly_cat.SetRedMass((938*1019.461)/(938+1019.461));
// set parameters for the interaction potential
//double& r = ppars[0];
//0 double& beta = ppars[2];
//1 double& a1 = ppars[3];
//2 double& b1 = ppars[4];//fm
//3 double& a2 = ppars[5];
//4 double& b2 = ppars[6];//fm
//5 double& a3mpi4 = ppars[7];
//6 double& b3 = ppars[8];//fm
//7 double& mpi = ppars[9];//MeV
	CATSparameters p01(CATSparameters::tPotential, 8, true);
	p01.SetParameter(0, 1);//beta
	p01.SetParameter(1, -392);//a1
  p01.SetParameter(2, 0.128);//b1
	p01.SetParameter(3, -145);//a2
  p01.SetParameter(4, 0.284);//b2
  p01.SetParameter(5, -83);//a3mpi4
  p01.SetParameter(6, 0.582);//b3
  p01.SetParameter(7, 146);//mpi

	smelly_cat.SetNumChannels(1);
	smelly_cat.SetNumPW(0, 1);
	smelly_cat.SetShortRangePotential(0,0,Emma_potential, p01);
	smelly_cat.SetSpin(0, 0);
	smelly_cat.SetChannelWeight(0, 1);
  smelly_cat.SetMaxRad(64);

	smelly_cat.KillTheCat();

  TFile fOutput(TString::Format("%s/OtherTasks/Ghetto_pPhi_check_20250426.root",GetFemtoOutputFolder()), "recreate");
  TH1F* hCk_pphi = new TH1F("hCk_pphi", "hCk_pphi", NumMomBins, kMin, kMax);

  for(unsigned uMom=0; uMom<NumMomBins; uMom++){
    hCk_pphi->SetBinContent(uMom+1, smelly_cat.GetCorrFun(uMom));

  }

  fOutput.cd();
  hCk_pphi->Write();
  delete hCk_pphi;
}


//save 2D histos, kstar vs rstar for the radial part of the WF
//run for both SI only and SI+coulomb
void pp_asymptotic(){
  TString OutputFileName = TString::Format("%s/OtherTasks/pp_asymptotic_av18.root",GetFemtoOutputFolder());

  unsigned NumMomBins = 60;
  double kMin = 2.5;
  double kMax = 302.5;

  unsigned NumRadBins = 512;
  double rMin = 0;
  double rMax = 24;


	CATS smelly_cat;
// set how many points, k_min, k-max
	smelly_cat.SetMomBins(NumMomBins, kMin, kMax);

  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles/",GetCernBoxDimi()));
  AnalysisObject.SetUpCats_pp(smelly_cat,"AV18","Gauss",0, 0);


  smelly_cat.KillTheCat();

  TFile fOutput(OutputFileName, "recreate");

  TH2F* h_pp_av18_1S0_SI = new TH2F("h_pp_av18_1S0_SI", "h_pp_av18_1S0_SI", NumMomBins, kMin, kMax, NumRadBins, rMin, rMax);
  TH2F* h_pp_av18_1S0_FULL = new TH2F("h_pp_av18_1S0_FULL", "h_pp_av18_1S0_FULL", NumMomBins, kMin, kMax, NumRadBins, rMin, rMax);

  TH2F* h_pp_av18_3P0_SI = new TH2F("h_pp_av18_3P0_SI", "h_pp_av18_3P0_SI", NumMomBins, kMin, kMax, NumRadBins, rMin, rMax);
  TH2F* h_pp_av18_3P0_FULL = new TH2F("h_pp_av18_3P0_FULL", "h_pp_av18_3P0_FULL", NumMomBins, kMin, kMax, NumRadBins, rMin, rMax);

  TH2F* h_pp_av18_3P1_SI = new TH2F("h_pp_av18_3P1_SI", "h_pp_av18_3P1_SI", NumMomBins, kMin, kMax, NumRadBins, rMin, rMax);
  TH2F* h_pp_av18_3P1_FULL = new TH2F("h_pp_av18_3P1_FULL", "h_pp_av18_3P1_FULL", NumMomBins, kMin, kMax, NumRadBins, rMin, rMax);

  TH2F* h_pp_av18_3P2_SI = new TH2F("h_pp_av18_3P2_SI", "h_pp_av18_3P2_SI", NumMomBins, kMin, kMax, NumRadBins, rMin, rMax);
  TH2F* h_pp_av18_3P2_FULL = new TH2F("h_pp_av18_3P2_FULL", "h_pp_av18_3P2_FULL", NumMomBins, kMin, kMax, NumRadBins, rMin, rMax);

  TH2F* h_pp_av18_1D2_SI = new TH2F("h_pp_av18_1D2_SI", "h_pp_av18_1D2_SI", NumMomBins, kMin, kMax, NumRadBins, rMin, rMax);
  TH2F* h_pp_av18_1D2_FULL = new TH2F("h_pp_av18_1D2_FULL", "h_pp_av18_1D2_FULL", NumMomBins, kMin, kMax, NumRadBins, rMin, rMax);

  
  for(unsigned uMom=0; uMom<NumMomBins; uMom++){
    for(unsigned uRad=0; uRad<NumRadBins; uRad++){
      double RAD = h_pp_av18_1S0_SI->GetYaxis()->GetBinCenter(uRad+1);
      h_pp_av18_1S0_FULL->SetBinContent(uMom+1, uRad+1, smelly_cat.EvalRadialWaveFunction(uMom, 0, 0, RAD, true).real());
      h_pp_av18_3P0_FULL->SetBinContent(uMom+1, uRad+1, smelly_cat.EvalRadialWaveFunction(uMom, 1, 1, RAD, true).real());
      h_pp_av18_3P1_FULL->SetBinContent(uMom+1, uRad+1, smelly_cat.EvalRadialWaveFunction(uMom, 2, 1, RAD, true).real());
      h_pp_av18_3P2_FULL->SetBinContent(uMom+1, uRad+1, smelly_cat.EvalRadialWaveFunction(uMom, 3, 1, RAD, true).real());
      h_pp_av18_1D2_FULL->SetBinContent(uMom+1, uRad+1, smelly_cat.EvalRadialWaveFunction(uMom, 0, 2, RAD, true).real());
    }
  }

  smelly_cat.SetQ1Q2(0);
  smelly_cat.KillTheCat();

  for(unsigned uMom=0; uMom<NumMomBins; uMom++){
    for(unsigned uRad=0; uRad<NumRadBins; uRad++){
      double RAD = h_pp_av18_1S0_SI->GetYaxis()->GetBinCenter(uRad+1);
      h_pp_av18_1S0_SI->SetBinContent(uMom+1, uRad+1, smelly_cat.EvalRadialWaveFunction(uMom, 0, 0, RAD, true).real());
      h_pp_av18_3P0_SI->SetBinContent(uMom+1, uRad+1, smelly_cat.EvalRadialWaveFunction(uMom, 1, 1, RAD, true).real());
      h_pp_av18_3P1_SI->SetBinContent(uMom+1, uRad+1, smelly_cat.EvalRadialWaveFunction(uMom, 2, 1, RAD, true).real());
      h_pp_av18_3P2_SI->SetBinContent(uMom+1, uRad+1, smelly_cat.EvalRadialWaveFunction(uMom, 3, 1, RAD, true).real());
      h_pp_av18_1D2_SI->SetBinContent(uMom+1, uRad+1, smelly_cat.EvalRadialWaveFunction(uMom, 0, 2, RAD, true).real());
    }
  }

  h_pp_av18_1S0_SI->GetXaxis()->SetTitle("k* (MeV)");
  h_pp_av18_1S0_FULL->GetXaxis()->SetTitle("k* (MeV)");
  h_pp_av18_1S0_SI->GetYaxis()->SetTitle("r* (fm)");
  h_pp_av18_1S0_FULL->GetYaxis()->SetTitle("r* (fm)");

  h_pp_av18_3P0_SI->GetXaxis()->SetTitle("k* (MeV)");
  h_pp_av18_3P0_FULL->GetXaxis()->SetTitle("k* (MeV)");
  h_pp_av18_3P0_SI->GetYaxis()->SetTitle("r* (fm)");
  h_pp_av18_3P0_FULL->GetYaxis()->SetTitle("r* (fm)");

  h_pp_av18_3P1_SI->GetXaxis()->SetTitle("k* (MeV)");
  h_pp_av18_3P1_FULL->GetXaxis()->SetTitle("k* (MeV)");
  h_pp_av18_3P1_SI->GetYaxis()->SetTitle("r* (fm)");
  h_pp_av18_3P1_FULL->GetYaxis()->SetTitle("r* (fm)");

  h_pp_av18_3P2_SI->GetXaxis()->SetTitle("k* (MeV)");
  h_pp_av18_3P2_FULL->GetXaxis()->SetTitle("k* (MeV)");
  h_pp_av18_3P2_SI->GetYaxis()->SetTitle("r* (fm)");
  h_pp_av18_3P2_FULL->GetYaxis()->SetTitle("r* (fm)");

  h_pp_av18_1D2_SI->GetXaxis()->SetTitle("k* (MeV)");
  h_pp_av18_1D2_FULL->GetXaxis()->SetTitle("k* (MeV)");
  h_pp_av18_1D2_SI->GetYaxis()->SetTitle("r* (fm)");
  h_pp_av18_1D2_FULL->GetYaxis()->SetTitle("r* (fm)");



  h_pp_av18_1S0_SI->Write();
  h_pp_av18_1S0_FULL->Write();

  h_pp_av18_3P0_SI->Write();
  h_pp_av18_3P0_FULL->Write();

  h_pp_av18_3P1_SI->Write();
  h_pp_av18_3P1_FULL->Write();

  h_pp_av18_3P2_SI->Write();
  h_pp_av18_3P2_FULL->Write();

  h_pp_av18_1D2_SI->Write();
  h_pp_av18_1D2_FULL->Write();



}



void pp_asymptotic_check_Epelbaum_v1(){
  TString OutputFileName = TString::Format("%s/OtherTasks/pp_asymptotic_av18_epel.root",GetFemtoOutputFolder());

  //8 start from zero
  //5 start from 2.5
  unsigned NumMomBins = 20;
  double kMin = 0;
  double kMax = 160;

  unsigned NumRadBins = 512;
  double rMin = 0;
  double rMax = 24;

  double source_size = 0.8;

	CATS smelly_cat_av18;
// set how many points, k_min, k-max
	smelly_cat_av18.SetMomBins(NumMomBins, kMin, kMax);

  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles/",GetCernBoxDimi()));
  AnalysisObject.SetUpCats_pp(smelly_cat_av18,"AV18","Gauss",0, 0);
  smelly_cat_av18.SetAnaSource(0, source_size);
  smelly_cat_av18.KillTheCat();

	CATS smelly_cat_epel;
// set how many points, k_min, k-max
	smelly_cat_epel.SetMomBins(NumMomBins, kMin, kMax);

  AnalysisObject.SetUpCats_pp(smelly_cat_epel,"Epelbaum","Gauss",0, 0);
  smelly_cat_epel.SetAnaSource(0, source_size);
  smelly_cat_epel.KillTheCat();

  TFile fOutput(OutputFileName, "recreate");

  TH2F* h_pp_av18_1S0_SI = new TH2F("h_pp_av18_1S0_SI", "h_pp_av18_1S0_SI", NumMomBins, kMin, kMax, NumRadBins, rMin, rMax);
  TH2F* h_pp_av18_1S0_FULL = new TH2F("h_pp_av18_1S0_FULL", "h_pp_av18_1S0_FULL", NumMomBins, kMin, kMax, NumRadBins, rMin, rMax);

  TH2F* h_pp_av18_3P0_SI = new TH2F("h_pp_av18_3P0_SI", "h_pp_av18_3P0_SI", NumMomBins, kMin, kMax, NumRadBins, rMin, rMax);
  TH2F* h_pp_av18_3P0_FULL = new TH2F("h_pp_av18_3P0_FULL", "h_pp_av18_3P0_FULL", NumMomBins, kMin, kMax, NumRadBins, rMin, rMax);

  TH2F* h_pp_av18_3P1_SI = new TH2F("h_pp_av18_3P1_SI", "h_pp_av18_3P1_SI", NumMomBins, kMin, kMax, NumRadBins, rMin, rMax);
  TH2F* h_pp_av18_3P1_FULL = new TH2F("h_pp_av18_3P1_FULL", "h_pp_av18_3P1_FULL", NumMomBins, kMin, kMax, NumRadBins, rMin, rMax);

  TH2F* h_pp_av18_3P2_SI = new TH2F("h_pp_av18_3P2_SI", "h_pp_av18_3P2_SI", NumMomBins, kMin, kMax, NumRadBins, rMin, rMax);
  TH2F* h_pp_av18_3P2_FULL = new TH2F("h_pp_av18_3P2_FULL", "h_pp_av18_3P2_FULL", NumMomBins, kMin, kMax, NumRadBins, rMin, rMax);

  TH2F* h_pp_av18_1D2_SI = new TH2F("h_pp_av18_1D2_SI", "h_pp_av18_1D2_SI", NumMomBins, kMin, kMax, NumRadBins, rMin, rMax);
  TH2F* h_pp_av18_1D2_FULL = new TH2F("h_pp_av18_1D2_FULL", "h_pp_av18_1D2_FULL", NumMomBins, kMin, kMax, NumRadBins, rMin, rMax);

  TH2F* h_pp_epel_1S0_SI = new TH2F("h_pp_epel_1S0_SI", "h_pp_epel_1S0_SI", NumMomBins, kMin, kMax, NumRadBins, rMin, rMax);
  TH2F* h_pp_epel_1S0_FULL = new TH2F("h_pp_epel_1S0_FULL", "h_pp_epel_1S0_FULL", NumMomBins, kMin, kMax, NumRadBins, rMin, rMax);

  TH2F* h_pp_epel_3P0_SI = new TH2F("h_pp_epel_3P0_SI", "h_pp_epel_3P0_SI", NumMomBins, kMin, kMax, NumRadBins, rMin, rMax);
  TH2F* h_pp_epel_3P0_FULL = new TH2F("h_pp_epel_3P0_FULL", "h_pp_epel_3P0_FULL", NumMomBins, kMin, kMax, NumRadBins, rMin, rMax);

  TH2F* h_pp_epel_3P1_SI = new TH2F("h_pp_epel_3P1_SI", "h_pp_epel_3P1_SI", NumMomBins, kMin, kMax, NumRadBins, rMin, rMax);
  TH2F* h_pp_epel_3P1_FULL = new TH2F("h_pp_epel_3P1_FULL", "h_pp_epel_3P1_FULL", NumMomBins, kMin, kMax, NumRadBins, rMin, rMax);

  TH2F* h_pp_epel_3P2_SI = new TH2F("h_pp_epel_3P2_SI", "h_pp_epel_3P2_SI", NumMomBins, kMin, kMax, NumRadBins, rMin, rMax);
  TH2F* h_pp_epel_3P2_FULL = new TH2F("h_pp_epel_3P2_FULL", "h_pp_epel_3P2_FULL", NumMomBins, kMin, kMax, NumRadBins, rMin, rMax);

  TH2F* h_pp_epel_1D2_SI = new TH2F("h_pp_epel_1D2_SI", "h_pp_epel_1D2_SI", NumMomBins, kMin, kMax, NumRadBins, rMin, rMax);
  TH2F* h_pp_epel_1D2_FULL = new TH2F("h_pp_epel_1D2_FULL", "h_pp_epel_1D2_FULL", NumMomBins, kMin, kMax, NumRadBins, rMin, rMax);

  TH1F** h_pp_av18_1S0_FULL_k = new TH1F* [NumMomBins];
  TH1F** h_pp_epel_1S0_FULL_k = new TH1F* [NumMomBins];

  TH1F** h_pp_av18_1D2_FULL_k = new TH1F* [NumMomBins];
  TH1F** h_pp_epel_1D2_FULL_k = new TH1F* [NumMomBins];

  TH1F** h_pp_av18_3P0_FULL_k = new TH1F* [NumMomBins];
  TH1F** h_pp_epel_3P0_FULL_k = new TH1F* [NumMomBins];

  TH1F** h_pp_av18_3P1_FULL_k = new TH1F* [NumMomBins];
  TH1F** h_pp_epel_3P1_FULL_k = new TH1F* [NumMomBins];

  TH1F** h_pp_av18_3P2_FULL_k = new TH1F* [NumMomBins];
  TH1F** h_pp_epel_3P2_FULL_k = new TH1F* [NumMomBins];

  for(unsigned uMom=0; uMom<NumMomBins; uMom++){
    double kstar = h_pp_av18_1S0_FULL->GetXaxis()->GetBinCenter(uMom+1);
    h_pp_av18_1S0_FULL_k[uMom] = new TH1F(TString::Format("h_pp_av18_1S0_FULL_k%.1f",kstar), TString::Format("h_pp_av18_1S0_FULL_k%.1f",kstar), NumRadBins, rMin, rMax);
    h_pp_epel_1S0_FULL_k[uMom] = new TH1F(TString::Format("h_pp_epel_1S0_FULL_k%.1f",kstar), TString::Format("h_pp_epel_1S0_FULL_k%.1f",kstar), NumRadBins, rMin, rMax);
    
    h_pp_av18_1D2_FULL_k[uMom] = new TH1F(TString::Format("h_pp_av18_1D2_FULL_k%.1f",kstar), TString::Format("h_pp_av18_1D2_FULL_k%.1f",kstar), NumRadBins, rMin, rMax);
    h_pp_epel_1D2_FULL_k[uMom] = new TH1F(TString::Format("h_pp_epel_1D2_FULL_k%.1f",kstar), TString::Format("h_pp_epel_1D2_FULL_k%.1f",kstar), NumRadBins, rMin, rMax);

    h_pp_av18_3P0_FULL_k[uMom] = new TH1F(TString::Format("h_pp_av18_3P0_FULL_k%.1f",kstar), TString::Format("h_pp_av18_3P0_FULL_k%.1f",kstar), NumRadBins, rMin, rMax);
    h_pp_epel_3P0_FULL_k[uMom] = new TH1F(TString::Format("h_pp_epel_3P0_FULL_k%.1f",kstar), TString::Format("h_pp_epel_3P0_FULL_k%.1f",kstar), NumRadBins, rMin, rMax);

    h_pp_av18_3P1_FULL_k[uMom] = new TH1F(TString::Format("h_pp_av18_3P1_FULL_k%.1f",kstar), TString::Format("h_pp_av18_3P1_FULL_k%.1f",kstar), NumRadBins, rMin, rMax);
    h_pp_epel_3P1_FULL_k[uMom] = new TH1F(TString::Format("h_pp_epel_3P1_FULL_k%.1f",kstar), TString::Format("h_pp_epel_3P1_FULL_k%.1f",kstar), NumRadBins, rMin, rMax);

    h_pp_av18_3P2_FULL_k[uMom] = new TH1F(TString::Format("h_pp_av18_3P2_FULL_k%.1f",kstar), TString::Format("h_pp_av18_3P2_FULL_k%.1f",kstar), NumRadBins, rMin, rMax);
    h_pp_epel_3P2_FULL_k[uMom] = new TH1F(TString::Format("h_pp_epel_3P2_FULL_k%.1f",kstar), TString::Format("h_pp_epel_3P2_FULL_k%.1f",kstar), NumRadBins, rMin, rMax);

    for(unsigned uRad=0; uRad<NumRadBins; uRad++){
      double RAD = h_pp_av18_1S0_SI->GetYaxis()->GetBinCenter(uRad+1);
      h_pp_av18_1S0_FULL->SetBinContent(uMom+1, uRad+1, smelly_cat_av18.EvalRadialWaveFunction(uMom, 0, 0, RAD, true).real());
      h_pp_av18_1S0_FULL_k[uMom]->SetBinContent(uRad+1, smelly_cat_av18.EvalRadialWaveFunction(uMom, 0, 0, RAD, true).real());
      h_pp_av18_3P0_FULL_k[uMom]->SetBinContent(uRad+1, smelly_cat_av18.EvalRadialWaveFunction(uMom, 1, 1, RAD, true).real());
      h_pp_av18_3P1_FULL_k[uMom]->SetBinContent(uRad+1, smelly_cat_av18.EvalRadialWaveFunction(uMom, 2, 1, RAD, true).real());
      h_pp_av18_3P2_FULL_k[uMom]->SetBinContent(uRad+1, smelly_cat_av18.EvalRadialWaveFunction(uMom, 3, 1, RAD, true).real());
      h_pp_av18_1D2_FULL_k[uMom]->SetBinContent(uRad+1, smelly_cat_av18.EvalRadialWaveFunction(uMom, 0, 2, RAD, true).real());
      h_pp_av18_3P0_FULL->SetBinContent(uMom+1, uRad+1, smelly_cat_av18.EvalRadialWaveFunction(uMom, 1, 1, RAD, true).real());
      h_pp_av18_3P1_FULL->SetBinContent(uMom+1, uRad+1, smelly_cat_av18.EvalRadialWaveFunction(uMom, 2, 1, RAD, true).real());
      h_pp_av18_3P2_FULL->SetBinContent(uMom+1, uRad+1, smelly_cat_av18.EvalRadialWaveFunction(uMom, 3, 1, RAD, true).real());
      h_pp_av18_1D2_FULL->SetBinContent(uMom+1, uRad+1, smelly_cat_av18.EvalRadialWaveFunction(uMom, 0, 2, RAD, true).real());

      h_pp_epel_1S0_FULL->SetBinContent(uMom+1, uRad+1, smelly_cat_epel.EvalRadialWaveFunction(uMom, 0, 0, RAD, true).real());
      h_pp_epel_1S0_FULL_k[uMom]->SetBinContent(uRad+1, smelly_cat_epel.EvalRadialWaveFunction(uMom, 0, 0, RAD, true).real());
      h_pp_epel_3P0_FULL_k[uMom]->SetBinContent(uRad+1, smelly_cat_epel.EvalRadialWaveFunction(uMom, 1, 1, RAD, true).real());
      h_pp_epel_3P1_FULL_k[uMom]->SetBinContent(uRad+1, smelly_cat_epel.EvalRadialWaveFunction(uMom, 4, 1, RAD, true).real());
      h_pp_epel_3P2_FULL_k[uMom]->SetBinContent(uRad+1, smelly_cat_epel.EvalRadialWaveFunction(uMom, 7, 1, RAD, true).real());
      h_pp_epel_1D2_FULL_k[uMom]->SetBinContent(uRad+1, smelly_cat_epel.EvalRadialWaveFunction(uMom, 0, 2, RAD, true).real());

      h_pp_epel_3P0_FULL->SetBinContent(uMom+1, uRad+1, smelly_cat_epel.EvalRadialWaveFunction(uMom, 1, 1, RAD, true).real());
      h_pp_epel_3P1_FULL->SetBinContent(uMom+1, uRad+1, smelly_cat_epel.EvalRadialWaveFunction(uMom, 2, 1, RAD, true).real());
      h_pp_epel_3P2_FULL->SetBinContent(uMom+1, uRad+1, smelly_cat_epel.EvalRadialWaveFunction(uMom, 3, 1, RAD, true).real());
      h_pp_epel_1D2_FULL->SetBinContent(uMom+1, uRad+1, smelly_cat_epel.EvalRadialWaveFunction(uMom, 0, 2, RAD, true).real());

    }
  }

  TGraph gCk_av18;
  gCk_av18.SetName("gCk_av18");
  TGraph gCk_epel;
  gCk_epel.SetName("gCk_epel");

  TGraph gCk_ratio;
  gCk_ratio.SetName("gCk_ratio");

  for(unsigned uMom=0; uMom<NumMomBins; uMom++){
    gCk_av18.SetPoint(uMom,smelly_cat_av18.GetMomentum(uMom),smelly_cat_av18.GetCorrFun(uMom));
    gCk_epel.SetPoint(uMom,smelly_cat_epel.GetMomentum(uMom),smelly_cat_epel.GetCorrFun(uMom));
    gCk_ratio.SetPoint(uMom,smelly_cat_epel.GetMomentum(uMom),smelly_cat_epel.GetCorrFun(uMom)/smelly_cat_av18.GetCorrFun(uMom));
  }

  smelly_cat_av18.SetQ1Q2(0);
  smelly_cat_av18.KillTheCat();

  smelly_cat_epel.SetQ1Q2(0);
  smelly_cat_epel.KillTheCat();

  for(unsigned uMom=0; uMom<NumMomBins; uMom++){
    for(unsigned uRad=0; uRad<NumRadBins; uRad++){
      double RAD = h_pp_av18_1S0_SI->GetYaxis()->GetBinCenter(uRad+1);
      h_pp_av18_1S0_SI->SetBinContent(uMom+1, uRad+1, smelly_cat_av18.EvalRadialWaveFunction(uMom, 0, 0, RAD, true).real());
      h_pp_av18_3P0_SI->SetBinContent(uMom+1, uRad+1, smelly_cat_av18.EvalRadialWaveFunction(uMom, 1, 1, RAD, true).real());
      h_pp_av18_3P1_SI->SetBinContent(uMom+1, uRad+1, smelly_cat_av18.EvalRadialWaveFunction(uMom, 2, 1, RAD, true).real());
      h_pp_av18_3P2_SI->SetBinContent(uMom+1, uRad+1, smelly_cat_av18.EvalRadialWaveFunction(uMom, 3, 1, RAD, true).real());
      h_pp_av18_1D2_SI->SetBinContent(uMom+1, uRad+1, smelly_cat_av18.EvalRadialWaveFunction(uMom, 0, 2, RAD, true).real());

      h_pp_epel_1S0_SI->SetBinContent(uMom+1, uRad+1, smelly_cat_epel.EvalRadialWaveFunction(uMom, 0, 0, RAD, true).real());
      h_pp_epel_3P0_SI->SetBinContent(uMom+1, uRad+1, smelly_cat_epel.EvalRadialWaveFunction(uMom, 1, 1, RAD, true).real());
      h_pp_epel_3P1_SI->SetBinContent(uMom+1, uRad+1, smelly_cat_epel.EvalRadialWaveFunction(uMom, 2, 1, RAD, true).real());
      h_pp_epel_3P2_SI->SetBinContent(uMom+1, uRad+1, smelly_cat_epel.EvalRadialWaveFunction(uMom, 3, 1, RAD, true).real());
      h_pp_epel_1D2_SI->SetBinContent(uMom+1, uRad+1, smelly_cat_epel.EvalRadialWaveFunction(uMom, 0, 2, RAD, true).real());
    }
  }

  //without F-waves
  smelly_cat_epel.SetQ1Q2(1);
  for(unsigned uCh=1; uCh<12; uCh++){
    smelly_cat_epel.SetChannelWeight(uCh,0);
  }
  smelly_cat_epel.SetChannelWeight(3,1./12.);
  smelly_cat_epel.SetChannelWeight(6,3./12.);
  smelly_cat_epel.SetChannelWeight(9,5./12.);
  smelly_cat_epel.KillTheCat();  
  TGraph gCk_epel_spd;
  TGraph gCkNet_epel_spd;
  gCk_epel_spd.SetName("gCk_epel_spd");
  gCkNet_epel_spd.SetName("gCkNet_epel_spd");
  for(unsigned uMom=0; uMom<NumMomBins; uMom++){
    gCk_epel_spd.SetPoint(uMom,smelly_cat_epel.GetMomentum(uMom),smelly_cat_epel.GetCorrFun(uMom));
    double original = gCk_epel.GetPointY(uMom);
    gCkNet_epel_spd.SetPoint(uMom,smelly_cat_epel.GetMomentum(uMom),original-smelly_cat_epel.GetCorrFun(uMom));
  }

  h_pp_av18_1S0_SI->GetXaxis()->SetTitle("k* (MeV)");
  h_pp_av18_1S0_FULL->GetXaxis()->SetTitle("r* (fm)");

  h_pp_av18_3P0_SI->GetXaxis()->SetTitle("k* (MeV)");
  h_pp_av18_3P0_FULL->GetXaxis()->SetTitle("r* (fm)");

  h_pp_av18_3P1_SI->GetXaxis()->SetTitle("k* (MeV)");
  h_pp_av18_3P1_FULL->GetXaxis()->SetTitle("r* (fm)");

  h_pp_av18_3P2_SI->GetXaxis()->SetTitle("k* (MeV)");
  h_pp_av18_3P2_FULL->GetXaxis()->SetTitle("r* (fm)");

  h_pp_av18_1D2_SI->GetXaxis()->SetTitle("k* (MeV)");
  h_pp_av18_1D2_FULL->GetXaxis()->SetTitle("r* (fm)");


  h_pp_epel_1S0_SI->GetXaxis()->SetTitle("k* (MeV)");
  h_pp_epel_1S0_FULL->GetXaxis()->SetTitle("r* (fm)");

  h_pp_epel_3P0_SI->GetXaxis()->SetTitle("k* (MeV)");
  h_pp_epel_3P0_FULL->GetXaxis()->SetTitle("r* (fm)");

  h_pp_epel_3P1_SI->GetXaxis()->SetTitle("k* (MeV)");
  h_pp_epel_3P1_FULL->GetXaxis()->SetTitle("r* (fm)");

  h_pp_epel_3P2_SI->GetXaxis()->SetTitle("k* (MeV)");
  h_pp_epel_3P2_FULL->GetXaxis()->SetTitle("r* (fm)");

  h_pp_epel_1D2_SI->GetXaxis()->SetTitle("k* (MeV)");
  h_pp_epel_1D2_FULL->GetXaxis()->SetTitle("r* (fm)");

  gCk_av18.Write();
  gCk_epel.Write();
  gCk_epel_spd.Write();
  gCkNet_epel_spd.Write();
  gCk_ratio.Write();

  h_pp_av18_1S0_SI->Write();
  h_pp_epel_1S0_SI->Write();
  h_pp_av18_1S0_FULL->Write();
  h_pp_epel_1S0_FULL->Write();

  h_pp_av18_3P0_SI->Write();
  h_pp_epel_3P0_SI->Write();
  h_pp_av18_3P0_FULL->Write();
  h_pp_epel_3P0_FULL->Write();

  h_pp_av18_3P1_SI->Write();
  h_pp_epel_3P1_SI->Write();
  h_pp_av18_3P1_FULL->Write();
  h_pp_epel_3P1_FULL->Write();

  h_pp_av18_3P2_SI->Write();
  h_pp_epel_3P2_SI->Write();
  h_pp_av18_3P2_FULL->Write();
  h_pp_epel_3P2_FULL->Write();

  h_pp_av18_1D2_SI->Write();
  h_pp_epel_1D2_SI->Write();
  h_pp_av18_1D2_FULL->Write();
  h_pp_epel_1D2_FULL->Write();

  for(unsigned uMom=0; uMom<NumMomBins; uMom++){
    h_pp_av18_1S0_FULL_k[uMom]->Write();
    h_pp_epel_1S0_FULL_k[uMom]->Write();

    h_pp_av18_1D2_FULL_k[uMom]->Write();
    h_pp_epel_1D2_FULL_k[uMom]->Write();  
    
    h_pp_av18_3P0_FULL_k[uMom]->Write();
    h_pp_epel_3P0_FULL_k[uMom]->Write();  
    
    h_pp_av18_3P1_FULL_k[uMom]->Write();
    h_pp_epel_3P1_FULL_k[uMom]->Write();  

    h_pp_av18_3P2_FULL_k[uMom]->Write();
    h_pp_epel_3P2_FULL_k[uMom]->Write();  
  }

}


//it sucks to be honest
void test_screen_coulomb(){

  TString OutputFolder = TString::Format("%s/OtherTasks/test_screen_coulomb/",GetFemtoOutputFolder());

  const unsigned NumMomBins = 75;
  const double kMin = 0;
  const double kMax = 300;
  const double source_size = 1.2;

  CATSparameters cSorPars(CATSparameters::tSource,1,true);
  cSorPars.SetParameter(0,source_size);

  CATS cat_pp_c;
  cat_pp_c.SetMomBins(NumMomBins,kMin,kMax);
  cat_pp_c.SetAnaSource(GaussSource, cSorPars);
  cat_pp_c.SetUseAnalyticSource(true);
  cat_pp_c.SetMomentumDependentSource(false);
  cat_pp_c.SetThetaDependentSource(false);
  cat_pp_c.SetExcludeFailedBins(false);
  cat_pp_c.SetQ1Q2(1);
  cat_pp_c.SetQuantumStatistics(false);
  cat_pp_c.SetRedMass( 0.5*Mass_p );
  cat_pp_c.SetNumChannels(1);
  cat_pp_c.SetSpin(0,0);
  cat_pp_c.SetChannelWeight(0, 1.);
  cat_pp_c.KillTheCat();

  CATS cat_pp_sc;
  cat_pp_sc.SetMomBins(NumMomBins,kMin,kMax);
  cat_pp_sc.SetAnaSource(GaussSource, cSorPars);
  cat_pp_sc.SetUseAnalyticSource(true);
  cat_pp_sc.SetMomentumDependentSource(false);
  cat_pp_sc.SetThetaDependentSource(false);
  cat_pp_sc.SetExcludeFailedBins(false);
  cat_pp_sc.SetQ1Q2(0);
  cat_pp_sc.SetQuantumStatistics(false);
  cat_pp_sc.SetRedMass( 0.5*Mass_p );
  cat_pp_sc.SetNumChannels(1);
  cat_pp_sc.SetNumPW(0,1);
  cat_pp_sc.SetSpin(0,0);
  cat_pp_sc.SetChannelWeight(0, 1.);
  CATSparameters cPotPars(CATSparameters::tPotential,3,true);
  cPotPars.SetParameter(0,1);
  cPotPars.SetParameter(1,64);
  cPotPars.SetParameter(2,0.1);
  cat_pp_sc.SetShortRangePotential(0,0,ScreenedCoulomb,cPotPars);
  cat_pp_sc.SetEpsilonConv(1e-12);
  cat_pp_sc.SetEpsilonProp(1e-8);
  cat_pp_sc.SetMaxRad(128);
  cat_pp_sc.SetMaxRho(8*cat_pp_sc.GetMaxRho());
  cat_pp_sc.KillTheCat();

  TGraph gCat_c;
  gCat_c.SetName("gCat_c");

  TGraph gCat_sc;
  gCat_sc.SetName("gCat_sc");

  for(unsigned uMom=0; uMom<NumMomBins; uMom++){
    double kstar = cat_pp_c.GetMomentum(uMom);
    gCat_c.SetPoint(uMom, kstar, cat_pp_c.GetCorrFun(uMom));
    gCat_sc.SetPoint(uMom, kstar, cat_pp_sc.GetCorrFun(uMom));
  }

  TFile fOutput(OutputFolder+"fOutput.root", "recreate");
  gCat_c.Write();
  gCat_sc.Write();
}

void big_core(){
  TString OutputFolder = TString::Format("%s/OtherTasks/",GetFemtoOutputFolder());

  const unsigned NumMomBins = 75;
  const double kMin = 0;
  const double kMax = 300;
  const double source_size = 1.2;

  CATSparameters cSorPars(CATSparameters::tSource,1,true);
  cSorPars.SetParameter(0,source_size);

  CATS cat_pp_c;
  cat_pp_c.SetMomBins(NumMomBins,kMin,kMax);
  cat_pp_c.SetAnaSource(GaussSource, cSorPars);
  cat_pp_c.SetUseAnalyticSource(true);
  cat_pp_c.SetMomentumDependentSource(false);
  cat_pp_c.SetThetaDependentSource(false);
  cat_pp_c.SetExcludeFailedBins(false);
  cat_pp_c.SetQ1Q2(0);
  cat_pp_c.SetQuantumStatistics(false);
  cat_pp_c.SetRedMass( 0.5*Mass_p );
  cat_pp_c.SetNumChannels(1);
  cat_pp_c.SetNumPW(0,1);
  cat_pp_c.SetSpin(0,0);
  cat_pp_c.SetChannelWeight(0, 1.);
  CATSparameters cPotPars(CATSparameters::tPotential,3,true);
  cPotPars.SetParameter(0,1e4);
  cPotPars.SetParameter(1,0.8);
  cPotPars.SetParameter(2,0.01);
  cat_pp_c.SetShortRangePotential(0,0,RepulsiveCore,cPotPars);
  cat_pp_c.KillTheCat();


  TGraph gCat_c1;
  gCat_c1.SetName("gCat_c1");
  for(unsigned uMom=0; uMom<NumMomBins; uMom++){
    double kstar = cat_pp_c.GetMomentum(uMom);
    gCat_c1.SetPoint(uMom, kstar, cat_pp_c.GetCorrFun(uMom));
  }


  TGraph gCat_c2;
  gCat_c2.SetName("gCat_c2");
  cat_pp_c.SetShortRangePotential(0,0,0,1e5);
  cat_pp_c.KillTheCat();  
  for(unsigned uMom=0; uMom<NumMomBins; uMom++){
    double kstar = cat_pp_c.GetMomentum(uMom);
    gCat_c2.SetPoint(uMom, kstar, cat_pp_c.GetCorrFun(uMom));
  }

  TGraph gCat_c3;
  gCat_c3.SetName("gCat_c3");
  cat_pp_c.SetShortRangePotential(0,0,0,1e6);
  cat_pp_c.KillTheCat();  
  for(unsigned uMom=0; uMom<NumMomBins; uMom++){
    double kstar = cat_pp_c.GetMomentum(uMom);
    gCat_c3.SetPoint(uMom, kstar, cat_pp_c.GetCorrFun(uMom));
  }



  TFile fOutput(OutputFolder+"big_core.root", "recreate");
  gCat_c1.Write();
  gCat_c2.Write();
  gCat_c3.Write();
}

CATS* cat_ps=NULL;
double fit_ps_2p(double* x, double* pars){
  if(!cat_ps) return 0;
  cat_ps->SetShortRangePotential(0,0,0,pars[0]);
  cat_ps->SetShortRangePotential(0,0,1,pars[1]);
  cat_ps->KillTheCat();
  return cat_ps->EvalPhaseShift(*x,0,0);
}
double fit_er_2p(double* x, double* pars){
  if(!cat_ps) return 0;
  cat_ps->SetShortRangePotential(0,0,0,pars[0]);
  cat_ps->SetShortRangePotential(0,0,1,pars[1]);
  cat_ps->KillTheCat();
  return *x/tan(cat_ps->EvalPhaseShift(*x,0,0));
}
double fit_er_4p(double* x, double* pars){
  if(!cat_ps) return 0;
  cat_ps->SetShortRangePotential(0,0,0,pars[0]);
  cat_ps->SetShortRangePotential(0,0,1,pars[1]);
  cat_ps->SetShortRangePotential(0,0,2,pars[2]);
  cat_ps->SetShortRangePotential(0,0,3,pars[3]);  
  cat_ps->KillTheCat();
  return *x/tan(cat_ps->EvalPhaseShift(*x,0,0));
}

//
//linear potential, with pars V0,R0,V1,R1
//the potential is equal to V0 at R0, and V1 at R1.
//for r<R0 it becomes constant equal to V0
//for r>R1 it becomes constant equal to 0
double USR_pot_1(double* pars){
  double& r = pars[0];
  double& V0 = pars[2];
  double& R0 = pars[3];
  double& V1 = pars[4];
  double& R1 = pars[5];  
  //the linear parameters
  double p1 = (V1-V0)/(R1-R0);
  double p0 = V0-p1*R0;
  if(r<R0) return V0;
  if(r>R1) return 0;
  return p0+p1*r;
  
  //if(r<pars[2]) return 0;
  //return p0*(1-r*p1);s
}

//super silly potential: a polynomial castrated at both ends
double USR_pot_2(double* pars){
  double& r = pars[0];
  double& max_r = pars[2];
  double& p0 = pars[3];
  double& p1 = pars[4];
  double& p2 = pars[5];
  //the linear parameters
  if(r>max_r) return 0;
  return p0+p1*r+p2*r*r;
  
  //if(r<pars[2]) return 0;
  //return p0*(1-r*p1);
}


//same as 1, but the cut off on the far end now depends on rho = r*kstar rather than r
double USR_pot_1k(double* pars){
  /*
  double& r = pars[0];
  double& k = pars[1];
  double& V0 = pars[2];
  double& rho0 = pars[3];
  double& V1 = pars[4];
  double& rho1 = pars[5];  
  //the linear parameters
  double p1 = (V1-V0)/(rho1-rho0);
  double p0 = V0-p1*R0;
  if(r*k<rho0) return V0;
  if(r*k>rho1) return 0;
  return p0+p1*r;
  */
  return 0;
  
  //if(r<pars[2]) return 0;
  //return p0*(1-r*p1);
}


//and ultra-short-ranged potential with a source which does not start from zero
void USR_Potential_ShiftedSource(){
  
  TString OutputFolder = TString::Format("%s/OtherTasks/USR_Potential_ShiftedSource/",GetFemtoOutputFolder());
  const double potential_depth = -10;
  const double shift = 0.8;
  const double source_size = 0.8;
  const double source_size_def = 1.2;

  unsigned kSteps = 70;
  double kMin = 0;
  double kMax = 280;

  unsigned kFit = 10;
  double kMaxFit = 40;

  unsigned rSteps = 1024;
  double rMin = 0;
  double rMax = 4;
  TH1F* hPot_Toy = new TH1F("hPot_Toy","hPot_Toy",rSteps,rMin,rMax);

  CATS Kitty_AV18_Def;
  Kitty_AV18_Def.SetMomBins(kSteps,kMin,kMax);
  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles/",GetCernBoxDimi()));
  AnalysisObject.SetUpCats_pp(Kitty_AV18_Def,"AV18_s","",0,0);
  CATSparameters cSorPars_Def(CATSparameters::tSource,1,true);
  cSorPars_Def.SetParameter(0,source_size_def);
  Kitty_AV18_Def.SetAnaSource(GaussSource, cSorPars_Def);
  Kitty_AV18_Def.SetUseAnalyticSource(true);
  Kitty_AV18_Def.SetMomentumDependentSource(false);
  Kitty_AV18_Def.SetThetaDependentSource(false);
  Kitty_AV18_Def.SetNormalizedSource(true);
  Kitty_AV18_Def.SetQ1Q2(0);
  Kitty_AV18_Def.KillTheCat();
  TH1F* hPS_AV18_Def = new TH1F("hPS_AV18_Def","hPS_AV18_Def",kSteps,kMin,kMax);
  TH1F* hER_AV18_Def = new TH1F("hER_AV18_Def","hER_AV18_Def",kSteps,kMin,kMax);
  for(unsigned uMom=0; uMom<kSteps; uMom++){
    double kstar = Kitty_AV18_Def.GetMomentum(uMom);
    hPS_AV18_Def->SetBinContent(uMom+1,Kitty_AV18_Def.GetPhaseShift(uMom,0,0));
    //hPS_AV18_Def->SetBinError(uMom+1, 0.01 * sqrt(kstar*0.5));
    hPS_AV18_Def->SetBinError(uMom+1, 0.001 + 0.0001 * kstar);

    double ps = Kitty_AV18_Def.GetPhaseShift(uMom,0,0);
    hER_AV18_Def->SetBinContent(uMom+1,kstar/tan(ps));
    //printf("%f %f %f\n",kstar,ps,kstar/tan(ps));
    double ps_err = kstar>=20?0.001*kstar*kstar/tan(ps):0.001* 20*kstar/tan(ps);
    hER_AV18_Def->SetBinError(uMom+1, ps_err);    
  }
  Kitty_AV18_Def.SetQ1Q2(1);
  Kitty_AV18_Def.KillTheCat();  
  TF1* fPS = new TF1("fPS",fit_ps_2p,kMin,kMaxFit,2);
  TF1* fER = new TF1("fER",fit_er_2p,kMin,kMaxFit,2);

  CATS Kitty_AV18;
  Kitty_AV18.SetMomBins(kSteps,kMin,kMax);
  AnalysisObject.SetUpCats_pp(Kitty_AV18,"AV18_s","",0,0);
  CATSparameters cSorPars(CATSparameters::tSource,2,true);
  cSorPars.SetParameter(0,source_size);
  cSorPars.SetParameter(1,shift);
  Kitty_AV18.SetAnaSource(GaussSourceShifted, cSorPars);
  Kitty_AV18.SetUseAnalyticSource(true);
  Kitty_AV18.SetMomentumDependentSource(false);
  Kitty_AV18.SetThetaDependentSource(false);
  Kitty_AV18.SetNormalizedSource(true);
  Kitty_AV18.KillTheCat();

  CATS Kitty_SW_Def_fit;
  Kitty_SW_Def_fit.SetMomBins(kFit,kMin,kMaxFit);
  AnalysisObject.SetUpCats_pp(Kitty_SW_Def_fit,"SW_s","",0,0);
  Kitty_SW_Def_fit.SetAnaSource(GaussSource, cSorPars_Def);
  Kitty_SW_Def_fit.SetUseAnalyticSource(true);
  Kitty_SW_Def_fit.SetMomentumDependentSource(false);
  Kitty_SW_Def_fit.SetThetaDependentSource(false);
  Kitty_SW_Def_fit.SetNormalizedSource(true);
  Kitty_SW_Def_fit.SetQ1Q2(0);
  Kitty_SW_Def_fit.KillTheCat();

  cat_ps = &Kitty_SW_Def_fit;
  cat_ps->SetNotifications(CATS::nWarning);
  fPS->SetParameter(0,-20);
  fPS->SetParLimits(0,-1000,0);
  fPS->SetParameter(1,1.0);
  fPS->SetParLimits(1,0.01,2);
  hPS_AV18_Def->Fit(fPS,"S, N, R, M");
  CATS Kitty_SW;
  Kitty_SW.SetMomBins(kSteps,kMin,kMax);
  AnalysisObject.SetUpCats_pp(Kitty_SW,"SW_s","",0,0);
  Kitty_SW.SetAnaSource(GaussSourceShifted, cSorPars);
  Kitty_SW.SetUseAnalyticSource(true);
  Kitty_SW.SetMomentumDependentSource(false);
  Kitty_SW.SetThetaDependentSource(false);
  Kitty_SW.SetNormalizedSource(true);
  Kitty_SW.SetShortRangePotential(0,0,0,fPS->GetParameter(0));
  Kitty_SW.SetShortRangePotential(0,0,1,fPS->GetParameter(1));
  Kitty_SW.KillTheCat();

  fER->SetParameter(0,-20);
  fER->SetParLimits(0,-1000,0);
  //fER->SetParameter(1,1.0);
  //fER->SetParLimits(1,0.01,2);
  fER->FixParameter(1,0.5);
  hER_AV18_Def->Fit(fER,"S, N, R, M");
  CATS Kitty_SW_ER;
  Kitty_SW_ER.SetMomBins(kSteps,kMin,kMax);
  AnalysisObject.SetUpCats_pp(Kitty_SW_ER,"SW_s","",0,0);
  Kitty_SW_ER.SetAnaSource(GaussSourceShifted, cSorPars);
  Kitty_SW_ER.SetUseAnalyticSource(true);
  Kitty_SW_ER.SetMomentumDependentSource(false);
  Kitty_SW_ER.SetThetaDependentSource(false);
  Kitty_SW_ER.SetNormalizedSource(true);
  Kitty_SW_ER.SetShortRangePotential(0,0,0,fER->GetParameter(0));
  Kitty_SW_ER.SetShortRangePotential(0,0,1,fER->GetParameter(1));
  Kitty_SW_ER.KillTheCat();

  //const double V_1S0 = -16200./Mass_p*1.;//c.a. 17.27 MeV
  //const double d_1S0 = 2.3039*1.;
  CATS Kitty_SW_Def;
  Kitty_SW_Def.SetMomBins(kSteps,kMin,kMax);
  AnalysisObject.SetUpCats_pp(Kitty_SW_Def,"SW_s","",0,0);
  Kitty_SW_Def.SetAnaSource(GaussSource, cSorPars_Def);
  Kitty_SW_Def.SetUseAnalyticSource(true);
  Kitty_SW_Def.SetMomentumDependentSource(false);
  Kitty_SW_Def.SetThetaDependentSource(false);
  Kitty_SW_Def.SetNormalizedSource(true);
  //Kitty_SW_Def.SetShortRangePotential(0,0,0,-17.27);
  //Kitty_SW_Def.SetShortRangePotential(0,0,1,2.3039);
  Kitty_SW_Def.SetShortRangePotential(0,0,0,-800);
  Kitty_SW_Def.SetShortRangePotential(0,0,1,0.5); 
  Kitty_SW_Def.SetEpsilonConv(1e-9);
  Kitty_SW_Def.SetEpsilonProp(1e-9);
  Kitty_SW_Def.KillTheCat();

  TH1F* hPS_SW_Def = new TH1F("hPS_SW_Def","hPS_SW_Def",kSteps,kMin,kMax);
  TH1F* hER_SW_Def = new TH1F("hER_SW_Def","hER_SW_Def",kSteps,kMin,kMax);
  for(unsigned uMom=0; uMom<kSteps; uMom++){
    double kstar = Kitty_SW_Def.GetMomentum(uMom);
    hPS_SW_Def->SetBinContent(uMom+1,Kitty_SW_Def.GetPhaseShift(uMom,0,0));
    //hPS_SW_Def->SetBinError(uMom+1, 0.01 * sqrt(kstar*0.5));

    double ps = Kitty_SW_Def.GetPhaseShift(uMom,0,0);
    hER_SW_Def->SetBinContent(uMom+1,kstar/tan(ps));
    //hER_SW_Def->SetBinError(uMom+1, 0.001 * sqrt(kstar*0.5) * kstar/tan(ps));    

  }

  TH1F* hPS_SW = new TH1F("hPS_SW","hPS_SW",kSteps,kMin,kMax);
  TH1F* hER_SW = new TH1F("hER_SW","hER_SW",kSteps,kMin,kMax);
  for(unsigned uMom=0; uMom<kSteps; uMom++){
    double kstar = Kitty_SW.GetMomentum(uMom);
    hPS_SW->SetBinContent(uMom+1,Kitty_SW.GetPhaseShift(uMom,0,0));
    double ps = Kitty_SW.GetPhaseShift(uMom,0,0);
    hER_SW->SetBinContent(uMom+1,kstar/tan(ps));
  }

  TH1F* hPS_SW_ER = new TH1F("hPS_SW_ER","hPS_SW_ER",kSteps,kMin,kMax);
  TH1F* hER_SW_ER = new TH1F("hER_SW_ER","hER_SW_ER",kSteps,kMin,kMax);
  for(unsigned uMom=0; uMom<kSteps; uMom++){
    double kstar = Kitty_SW_ER.GetMomentum(uMom);
    hPS_SW_ER->SetBinContent(uMom+1,Kitty_SW_ER.GetPhaseShift(uMom,0,0));
    double ps = Kitty_SW_ER.GetPhaseShift(uMom,0,0);
    hER_SW_ER->SetBinContent(uMom+1,kstar/tan(ps));
  }


  CATS Kitty_ToyFit;
  Kitty_ToyFit.SetMomBins(kFit,kMin,kMaxFit);
  //AnalysisObject.SetUpCats_pp(Kitty_Toy,"","",0,0);
  Kitty_ToyFit.SetAnaSource(GaussSourceShifted, cSorPars);
  Kitty_ToyFit.SetUseAnalyticSource(true);
  Kitty_ToyFit.SetMomentumDependentSource(false);
  Kitty_ToyFit.SetThetaDependentSource(false);
  Kitty_ToyFit.SetNormalizedSource(true);

  Kitty_ToyFit.SetExcludeFailedBins(false);
  Kitty_ToyFit.SetQ1Q2(1);
  Kitty_ToyFit.SetQuantumStatistics(true);
  Kitty_ToyFit.SetRedMass( 0.5*Mass_p );
  Kitty_ToyFit.SetNumChannels(2);
  Kitty_ToyFit.SetNumPW(0,1);
  Kitty_ToyFit.SetSpin(0,0);
  Kitty_ToyFit.SetSpin(1,1);
  Kitty_ToyFit.SetChannelWeight(0, 1./4.);
  Kitty_ToyFit.SetChannelWeight(1, 3./4.);
  CATSparameters cToyPars(CATSparameters::tPotential,4,true);
  Kitty_ToyFit.SetShortRangePotential(0,0,USR_pot_1,cToyPars);

  Kitty_ToyFit.SetShortRangePotential(0,0,0,1000);
  Kitty_ToyFit.SetShortRangePotential(0,0,1,0);  
  Kitty_ToyFit.SetShortRangePotential(0,0,2,-540);
  Kitty_ToyFit.SetShortRangePotential(0,0,3,0.8);

  Kitty_ToyFit.KillTheCat();

  TF1* fER_Toy = new TF1("fER_Toy",fit_er_4p,kMin,kMaxFit,4);

  cat_ps = &Kitty_ToyFit;
  cat_ps->SetNotifications(CATS::nWarning);
  fER_Toy->SetParameter(0,1000);
  fER_Toy->SetParLimits(0,0,10000);
  fER_Toy->SetParameter(1,0);
  fER_Toy->SetParLimits(1,0,0.3);
  fER_Toy->SetParameter(2,-400);
  fER_Toy->SetParLimits(2,-1000,0);
  fER_Toy->SetParameter(3,0.6);
  fER_Toy->SetParLimits(3,0.4,0.8);  

  fER_Toy->FixParameter(0,1000);
  fER_Toy->FixParameter(1,0);
  fER_Toy->FixParameter(2,-520);
  fER_Toy->FixParameter(3,0.8);

  hPS_AV18_Def->Fit(fER_Toy,"S, N, R, M");

  CATS Kitty_Toy;
  Kitty_Toy.SetMomBins(kSteps,kMin,kMax);
  //AnalysisObject.SetUpCats_pp(Kitty_Toy,"","",0,0);
  Kitty_Toy.SetAnaSource(GaussSourceShifted, cSorPars);
  Kitty_Toy.SetUseAnalyticSource(true);
  Kitty_Toy.SetMomentumDependentSource(false);
  Kitty_Toy.SetThetaDependentSource(false);
  Kitty_Toy.SetNormalizedSource(true);

  Kitty_Toy.SetExcludeFailedBins(false);
  Kitty_Toy.SetQ1Q2(1);
  Kitty_Toy.SetQuantumStatistics(true);
  Kitty_Toy.SetRedMass( 0.5*Mass_p );
  Kitty_Toy.SetNumChannels(2);
  Kitty_Toy.SetNumPW(0,1);
  Kitty_Toy.SetSpin(0,0);
  Kitty_Toy.SetSpin(1,1);
  Kitty_Toy.SetChannelWeight(0, 1./4.);
  Kitty_Toy.SetChannelWeight(1, 3./4.);
  Kitty_Toy.SetShortRangePotential(0,0,USR_pot_1,cToyPars);

  //Kitty_Toy.SetShortRangePotential(0,0,0,1000);
  //Kitty_Toy.SetShortRangePotential(0,0,1,0);  
  //Kitty_Toy.SetShortRangePotential(0,0,2,-540);
  //Kitty_Toy.SetShortRangePotential(0,0,3,0.8);
  Kitty_Toy.SetShortRangePotential(0,0,0,fER_Toy->GetParameter(0));
  Kitty_Toy.SetShortRangePotential(0,0,1,fER_Toy->GetParameter(1));  
  Kitty_Toy.SetShortRangePotential(0,0,2,fER_Toy->GetParameter(2));
  Kitty_Toy.SetShortRangePotential(0,0,3,fER_Toy->GetParameter(3));


  

  Kitty_Toy.KillTheCat();

  TH1F* hPS_Toy = new TH1F("hPS_Toy","hPS_Toy",kSteps,kMin,kMax);
  TH1F* hER_Toy = new TH1F("hER_Toy","hER_Toy",kSteps,kMin,kMax);
  for(unsigned uMom=0; uMom<kSteps; uMom++){
    double kstar = Kitty_SW_Def.GetMomentum(uMom);
    hPS_Toy->SetBinContent(uMom+1,Kitty_Toy.GetPhaseShift(uMom,0,0));
    double ps = Kitty_Toy.GetPhaseShift(uMom,0,0);
    hER_Toy->SetBinContent(uMom+1,kstar/tan(ps));
  }

  for(unsigned uRad=0; uRad<rSteps; uRad++){
    double rad = hPot_Toy->GetBinCenter(uRad+1);
    double pot = Kitty_Toy.EvaluateThePotential(0,0,10,rad);
    hPot_Toy->SetBinContent(uRad+1, pot);
  }
  
  TGraph g_AV18_Def;
  g_AV18_Def.SetName("g_AV18_Def");
  g_AV18_Def.SetLineWidth(4);
  g_AV18_Def.SetLineColor(kBlue);

  TGraph g_AV18;
  g_AV18.SetName("g_AV18");
  g_AV18.SetLineWidth(4);
  g_AV18.SetLineColor(kCyan+1);

  TGraph g_SW_Def;
  g_SW_Def.SetName("g_SW_Def");
  g_SW_Def.SetLineWidth(4);
  g_SW_Def.SetLineColor(kRed);

  TGraph g_SW;
  g_SW.SetName("g_SW");
  g_SW.SetLineWidth(4);
  g_SW.SetLineColor(kMagenta);

  TGraph g_Toy;
  g_Toy.SetName("g_Toy");
  g_Toy.SetLineWidth(4);
  g_Toy.SetLineColor(kYellow+2);

  for(unsigned uMom=0; uMom<kSteps; uMom++){
    double kstar = Kitty_AV18_Def.GetMomentum(uMom);
    g_AV18_Def.SetPoint(uMom, kstar, Kitty_AV18_Def.GetCorrFun(uMom));
    g_AV18.SetPoint(uMom, kstar, Kitty_AV18.GetCorrFun(uMom));
    g_SW_Def.SetPoint(uMom, kstar, Kitty_SW_Def.GetCorrFun(uMom));
    g_SW.SetPoint(uMom, kstar, Kitty_SW.GetCorrFun(uMom));
    g_Toy.SetPoint(uMom, kstar, Kitty_Toy.GetCorrFun(uMom));
  }

  TFile fOutput(OutputFolder+"f1.root", "recreate");
  g_AV18_Def.Write();
  g_AV18.Write();
  g_SW_Def.Write();
  g_SW.Write();
  g_Toy.Write();

  hPS_AV18_Def->Write();
  hPS_SW_Def->Write();
  hPS_SW->Write();
  hPS_SW_ER->Write();
  hPS_Toy->Write();

  hER_AV18_Def->Write();
  hER_SW_Def->Write();
  hER_SW->Write();
  hER_SW_ER->Write();
  hER_Toy->Write();

  hPot_Toy->Write();
}


//DG super close to AV18:
//V0 = 324.264
//m0 = 0.16704
//V1 = -33.0925
//m1 = 1.74996
void Toy_Potential_Scan(int SEED, unsigned NumIter){


  TString OutputFile = TString::Format("%s/OtherTasks/Toy_Potential_Scan/dg_seed_%u.root",GetFemtoOutputFolder(),SEED);
  unsigned kSteps = 16;
  double kMin = 0;
  double kMax = 160;
  double kMin_ER = 10;
  double kMax_ER = 80;

  double source_size = 1.2;
  double shift = 0.8;

  CATSparameters cSorPars(CATSparameters::tSource,2,true);
  cSorPars.SetParameter(0,source_size);
  cSorPars.SetParameter(1,shift);


  CATS Kitty_AV18_Def;
  Kitty_AV18_Def.SetMomBins(kSteps,kMin,kMax);
  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles/",GetCernBoxDimi()));
  AnalysisObject.SetUpCats_pp(Kitty_AV18_Def,"AV18_s","",0,0);
  Kitty_AV18_Def.SetAnaSource(GaussSource, cSorPars);
  Kitty_AV18_Def.SetUseAnalyticSource(true);
  Kitty_AV18_Def.SetMomentumDependentSource(false);
  Kitty_AV18_Def.SetThetaDependentSource(false);
  Kitty_AV18_Def.SetNormalizedSource(true);
  Kitty_AV18_Def.SetQ1Q2(0);
  Kitty_AV18_Def.KillTheCat();

  TH1F* hPS_AV18 = new TH1F(TString::Format("hPS_AV18"),TString::Format("hPS_AV18"),kSteps,kMin,kMax);
  TH1F* hER_AV18 = new TH1F(TString::Format("hER_AV18"),TString::Format("hER_AV18"),kSteps,kMin,kMax);
  //TH1F* hCk_AV18 = new TH1F(TString::Format("hCk_AV18"),TString::Format("hCk_AV18"),kSteps,kMin,kMax);;
  for(unsigned uMom=0; uMom<kSteps; uMom++){
    double ps = Kitty_AV18_Def.GetPhaseShift(uMom,0,0);
    double kstar = Kitty_AV18_Def.GetMomentum(uMom);

    hPS_AV18->SetBinContent(uMom+1, ps);
    hPS_AV18->SetBinError(uMom+1, 0.01);

    hER_AV18->SetBinContent(uMom+1, kstar/tan(ps));
    hER_AV18->SetBinError(uMom+1, 1);

  }


  CATS Kitty_Toy;
  Kitty_Toy.SetMomBins(kSteps,kMin,kMax);
  //AnalysisObject.SetUpCats_pp(Kitty_Toy,"","",0,0);
  //Kitty_Toy.SetAnaSource(GaussSourceShifted, cSorPars);
  Kitty_Toy.SetAnaSource(GaussSource, cSorPars);
  Kitty_Toy.SetUseAnalyticSource(true);
  Kitty_Toy.SetMomentumDependentSource(false);
  Kitty_Toy.SetThetaDependentSource(false);
  Kitty_Toy.SetNormalizedSource(true);

  Kitty_Toy.SetExcludeFailedBins(false);
  Kitty_Toy.SetQ1Q2(0);
  Kitty_Toy.SetQuantumStatistics(true);
  Kitty_Toy.SetRedMass( 0.5*Mass_p );
  Kitty_Toy.SetNumChannels(2);
  Kitty_Toy.SetNumPW(0,1);
  Kitty_Toy.SetSpin(0,0);
  Kitty_Toy.SetSpin(1,1);
  Kitty_Toy.SetChannelWeight(0, 1./4.);
  Kitty_Toy.SetChannelWeight(1, 3./4.);
  CATSparameters cToyPars(CATSparameters::tPotential,4,true);
  //Kitty_Toy.SetShortRangePotential(0,0,USR_pot_1,cToyPars);
  Kitty_Toy.SetShortRangePotential(0,0,DoubleGaussSum,cToyPars);

  TFile fOutput(OutputFile, "recreate");

  TNtuple* ntToyPot = new TNtuple("ntToyPot", "ntToyPot","SEED:ID:V0:R0:V1:R1:shift:f0:d0:chi2:diff_av18");
  Float_t ntBuffer[11];

  TRandom3 rangen(SEED);
  double best_diff_av18 = 1e6;
  for(unsigned uIter=0; uIter<NumIter; uIter++){


    TH1F* hPS = new TH1F(TString::Format("hPS_%u_%u",SEED,uIter),TString::Format("hPS_%u_%u",SEED,uIter),kSteps,kMin,kMax);
    TH1F* hER = new TH1F(TString::Format("hER_%u_%u",SEED,uIter),TString::Format("hER_%u_%u",SEED,uIter),kSteps,kMin,kMax);
    TF1* fER = new TF1(TString::Format("fER_%u_%u",SEED,uIter),"197.327*[0]+0.5*[1]/197.327*x*x",kMin_ER,kMax_ER);
    TH1F* hPot = new TH1F(TString::Format("hPot_%u_%u",SEED,uIter),TString::Format("hPot_%u_%u",SEED,uIter),2048,0,2);



    Kitty_Toy.SetShortRangePotential(0,0,0,324.264);
    Kitty_Toy.SetShortRangePotential(0,0,1,0.16704);  
    Kitty_Toy.SetShortRangePotential(0,0,2,-33.0925);
    Kitty_Toy.SetShortRangePotential(0,0,3,1.74996);
    //DG (up is for good descr. of av18)
    //Kitty_Toy.SetShortRangePotential(0,0,0,rangen.Uniform(2000,4500));
    //Kitty_Toy.SetShortRangePotential(0,0,1,rangen.Uniform(0.3,0.7));  
    //Kitty_Toy.SetShortRangePotential(0,0,2,rangen.Uniform(-300,-50));
    //Kitty_Toy.SetShortRangePotential(0,0,3,rangen.Uniform(0.7,2.0));

    //toy1
    //Kitty_Toy.SetShortRangePotential(0,0,0,rangen.Uniform(500,1500));
    //Kitty_Toy.SetShortRangePotential(0,0,1,rangen.Uniform(0,shift*0.8));  
    //Kitty_Toy.SetShortRangePotential(0,0,2,rangen.Uniform(-2000,-800));
    //Kitty_Toy.SetShortRangePotential(0,0,3,rangen.Uniform(shift*0.8,shift));

    //toy1k
    //Kitty_Toy.SetShortRangePotential(0,0,0,rangen.Uniform(0,1500));
    //Kitty_Toy.SetShortRangePotential(0,0,1,rangen.Uniform(0,shift*0.5));  
    //Kitty_Toy.SetShortRangePotential(0,0,2,rangen.Uniform(-1500,0));
    //Kitty_Toy.SetShortRangePotential(0,0,3,rangen.Uniform(shift*0.5,shift));

    //toy2
    //Kitty_Toy.SetShortRangePotential(0,0,0,rangen.Uniform(shift));
    //Kitty_Toy.SetShortRangePotential(0,0,1,rangen.Uniform(-1000,1000));  
    //Kitty_Toy.SetShortRangePotential(0,0,2,rangen.Uniform(-10, 10));
    //Kitty_Toy.SetShortRangePotential(0,0,3,rangen.Uniform(-1, 1));


    //Kitty_Toy.SetShortRangePotential(0,0,0,100);
    //Kitty_Toy.SetShortRangePotential(0,0,1,0);  
    //Kitty_Toy.SetShortRangePotential(0,0,2,-150);
    //Kitty_Toy.SetShortRangePotential(0,0,3,0.8);    
    Kitty_Toy.SetNotifications(CATS::nWarning);
    Kitty_Toy.KillTheCat();

    double diff_av18 = 0;
    for(unsigned uMom=0; uMom<kSteps; uMom++){
      double ps = Kitty_Toy.GetPhaseShift(uMom,0,0);
      double kstar = Kitty_Toy.GetMomentum(uMom);
      //printf("%f %f %f\n",kstar,ps,Kitty_Toy.GetCorrFun(uMom));
      hER->SetBinContent(uMom+1, kstar/tan(ps));
      hER->SetBinError(uMom+1, 1);

      hPS->SetBinContent(uMom+1, ps);
      hPS->SetBinError(uMom+1, 0.01);

      if(kstar>=kMin_ER && kstar<=kMax_ER){
        diff_av18 += pow(ps-hPS_AV18->GetBinContent(uMom+1),2.);
      }
    }
    fER->SetParameter(0, 0);
    fER->SetParLimits(0, -10,10);
    fER->SetParLimits(1, -50,150);

    hER->Fit(fER,"Q, S, N, R, M");\
    //if(1/fER->GetParameter(0)>10)

    for(unsigned uRad=0; uRad<hPot->GetNbinsX(); uRad++){
      float rad = hPot->GetBinCenter(uRad+1);
      hPot->SetBinContent(uRad+1, Kitty_Toy.EvaluateThePotential(0,0,10,rad));
    }

    if(fER->GetChisquare()<100){
      if(diff_av18<best_diff_av18){

        printf("%u: f; d = %.2f/%.2f; chi2=%.0f; diff=%.3f\n", uIter, 1/fER->GetParameter(0), fER->GetParameter(1), fER->GetChisquare(), diff_av18);
        printf("    V0,R0,V1,R1 = %.3f, %.3f, %.3f, %.3f\n", Kitty_Toy.GetPotPar(0,0,0),Kitty_Toy.GetPotPar(0,0,1),Kitty_Toy.GetPotPar(0,0,2),Kitty_Toy.GetPotPar(0,0,3));

        hPS->Write();
        hER->Write();
        fER->Write();
        hPot->Write();

        best_diff_av18 = diff_av18;
      }
      if(diff_av18<0.1){
        ntBuffer[0] = SEED;
        ntBuffer[1] = uIter;
        ntBuffer[2] = Kitty_Toy.GetPotPar(0,0,0);
        ntBuffer[3] = Kitty_Toy.GetPotPar(0,0,1);
        ntBuffer[4] = Kitty_Toy.GetPotPar(0,0,2);
        ntBuffer[5] = Kitty_Toy.GetPotPar(0,0,3);
        ntBuffer[6] = shift;
        ntBuffer[7] = 1/fER->GetParameter(0);
        ntBuffer[8] = fER->GetParameter(1);
        ntBuffer[9] = fER->GetChisquare();
        ntBuffer[10] = diff_av18;
        ntToyPot->Fill(ntBuffer);

      }
    }



    

    
    delete hPS;
    delete hER;
    delete fER;
    delete hPot;
  }
  //Kitty_Toy.SetShortRangePotential(0,0,0,1000);
  //Kitty_Toy.SetShortRangePotential(0,0,1,0);  
  //Kitty_Toy.SetShortRangePotential(0,0,2,-540);
  //Kitty_Toy.SetShortRangePotential(0,0,3,0.8);
  hPS_AV18->Write();
  hER_AV18->Write();
  ntToyPot->Write();

  delete hPS_AV18;
  delete ntToyPot;
}


void Toy_Potential_Compare(){


  TString OutputFile = TString::Format("%s/OtherTasks/Toy_Potential_Scan/Toy_Potential_Compare.root",GetFemtoOutputFolder());
  unsigned kSteps = 120;
  double kMin = 0;
  double kMax = 240;

  double source_size = 0.7;

  CATSparameters cSorPars(CATSparameters::tSource,1,true);
  cSorPars.SetParameter(0,source_size);


  CATS Kitty_AV18_Def;
  Kitty_AV18_Def.SetMomBins(kSteps,kMin,kMax);
  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles/",GetCernBoxDimi()));
  AnalysisObject.SetUpCats_pp(Kitty_AV18_Def,"AV18_s","",0,0);
  Kitty_AV18_Def.SetAnaSource(GaussSource, cSorPars);
  Kitty_AV18_Def.SetUseAnalyticSource(true);
  Kitty_AV18_Def.SetMomentumDependentSource(false);
  Kitty_AV18_Def.SetThetaDependentSource(false);
  Kitty_AV18_Def.SetNormalizedSource(true);
  Kitty_AV18_Def.SetQ1Q2(0);
  Kitty_AV18_Def.SetEpsilonConv(1e-9);
  Kitty_AV18_Def.SetEpsilonProp(1e-9);
  Kitty_AV18_Def.KillTheCat();


  CATS Kitty_Toy;
  Kitty_Toy.SetMomBins(kSteps,kMin,kMax);
  Kitty_Toy.SetAnaSource(GaussSource, cSorPars);
  Kitty_Toy.SetUseAnalyticSource(true);
  Kitty_Toy.SetMomentumDependentSource(false);
  Kitty_Toy.SetThetaDependentSource(false);
  Kitty_Toy.SetNormalizedSource(true);

  Kitty_Toy.SetExcludeFailedBins(false);
  Kitty_Toy.SetQ1Q2(0);
  Kitty_Toy.SetQuantumStatistics(true);
  Kitty_Toy.SetRedMass( 0.5*Mass_p );
  Kitty_Toy.SetNumChannels(2);
  Kitty_Toy.SetNumPW(0,1);
  Kitty_Toy.SetSpin(0,0);
  Kitty_Toy.SetSpin(1,1);
  Kitty_Toy.SetChannelWeight(0, 1./4.);
  Kitty_Toy.SetChannelWeight(1, 3./4.);
  CATSparameters cToyPars(CATSparameters::tPotential,4,true);
  Kitty_Toy.SetShortRangePotential(0,0,DoubleGaussSum,cToyPars);
  Kitty_Toy.SetShortRangePotential(0,0,0,324.264);
  Kitty_Toy.SetShortRangePotential(0,0,1,0.16704);  
  Kitty_Toy.SetShortRangePotential(0,0,2,-33.0925);
  Kitty_Toy.SetShortRangePotential(0,0,3,1.74996);
  Kitty_Toy.SetEpsilonConv(1e-9);
  Kitty_Toy.SetEpsilonProp(1e-9);
  Kitty_Toy.KillTheCat();

  TFile fOutput(OutputFile, "recreate");


  TH1F* hPS_AV18 = new TH1F(TString::Format("hPS_AV18"),TString::Format("hPS_AV18"),kSteps,kMin,kMax);
  TH1F* hER_AV18 = new TH1F(TString::Format("hER_AV18"),TString::Format("hER_AV18"),kSteps,kMin,kMax);
  TH1F* hCk_AV18 = new TH1F(TString::Format("hCk_AV18"),TString::Format("hCk_AV18"),kSteps,kMin,kMax);
  TH1F* hdCdk_AV18 = new TH1F(TString::Format("hdCdk_AV18"),TString::Format("hdCdk_AV18"),kSteps,kMin,kMax);
  TH1F* hPot_AV18 = new TH1F(TString::Format("hPot_AV18"),TString::Format("hPot_AV18"),2048,0,2);

  TH1F* hPS = new TH1F(TString::Format("hPS_DG"),TString::Format("hPS_DG"),kSteps,kMin,kMax);
  TH1F* hER = new TH1F(TString::Format("hER_DG"),TString::Format("hER_DG"),kSteps,kMin,kMax);
  TH1F* hCk = new TH1F(TString::Format("hCk_DG"),TString::Format("hCk_DG"),kSteps,kMin,kMax);
  TH1F* hdCdk = new TH1F(TString::Format("hdCdk"),TString::Format("hdCdk"),kSteps,kMin,kMax);
  TH1F* hPot = new TH1F(TString::Format("hPot_DG"),TString::Format("hPot_DG"),2048,0,2);


  for(unsigned uMom=0; uMom<kSteps; uMom++){
    double ps = Kitty_AV18_Def.GetPhaseShift(uMom,0,0);
    double kstar = Kitty_AV18_Def.GetMomentum(uMom);

    hPS_AV18->SetBinContent(uMom+1, ps);
    hPS_AV18->SetBinError(uMom+1, 0.01);

    hER_AV18->SetBinContent(uMom+1, kstar/tan(ps));
    hER_AV18->SetBinError(uMom+1, 1);

    hCk_AV18->SetBinContent(uMom+1, Kitty_AV18_Def.GetCorrFun(uMom));
    if(uMom){
      double dC = Kitty_AV18_Def.GetCorrFun(uMom) - Kitty_AV18_Def.GetCorrFun(uMom-1);
      double dk = Kitty_AV18_Def.GetMomentum(uMom) - Kitty_AV18_Def.GetMomentum(uMom-1);
      hdCdk_AV18->SetBinContent(uMom, dC/dk );
    }

    ps = Kitty_Toy.GetPhaseShift(uMom,0,0);
    hPS->SetBinContent(uMom+1, ps);
    hPS->SetBinError(uMom+1, 0.01);

    hER->SetBinContent(uMom+1, kstar/tan(ps));
    hER->SetBinError(uMom+1, 1);

    hCk->SetBinContent(uMom+1, Kitty_Toy.GetCorrFun(uMom));
    if(uMom){
      double dC = Kitty_Toy.GetCorrFun(uMom) - Kitty_Toy.GetCorrFun(uMom-1);
      double dk = Kitty_Toy.GetMomentum(uMom) - Kitty_Toy.GetMomentum(uMom-1);
      hdCdk->SetBinContent(uMom, dC/dk );
    }
  }

  for(unsigned uRad=0; uRad<hPot->GetNbinsX(); uRad++){
    float rad = hPot->GetBinCenter(uRad+1);
    hPot_AV18->SetBinContent(uRad+1, Kitty_AV18_Def.EvaluateThePotential(0,0,10,rad));
    hPot->SetBinContent(uRad+1, Kitty_Toy.EvaluateThePotential(0,0,10,rad));
  }


  hPS_AV18->Write();
  hPS->Write();
  hER->Write();
  hER_AV18->Write();
  hCk_AV18->Write();
  hCk->Write();
  hdCdk_AV18->Write();
  hdCdk->Write();
  hPot_AV18->Write();
  hPot->Write();

  delete hPS;
  delete hER;
  delete hCk;
  delete hdCdk;
  delete hPot;

  delete hPS_AV18;
  delete hER_AV18;
  delete hCk_AV18;
  delete hdCdk_AV18;
  delete hPot_AV18;
}


void Toy_Potential_Scan_FloatGauss(int SEED, unsigned NumIter){


  TString OutputFile = TString::Format("%s/OtherTasks/Toy_Potential_Scan_FloatGauss/seed_%u.root",GetFemtoOutputFolder(),SEED);
  unsigned kSteps = 16;
  double kMin = 0;
  double kMax = 80;
  double kMin_ER = 0;
  double kMax_ER = 60;

  double source_size = 1.2;

  CATSparameters cSorPars(CATSparameters::tSource,1,true);
  cSorPars.SetParameter(0,source_size);


  CATS Kitty_AV18_Def;
  Kitty_AV18_Def.SetMomBins(kSteps,kMin,kMax);
  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles/",GetCernBoxDimi()));
  AnalysisObject.SetUpCats_pp(Kitty_AV18_Def,"AV18_s","",0,0);
  Kitty_AV18_Def.SetAnaSource(GaussSource, cSorPars);
  Kitty_AV18_Def.SetUseAnalyticSource(true);
  Kitty_AV18_Def.SetMomentumDependentSource(false);
  Kitty_AV18_Def.SetThetaDependentSource(false);
  Kitty_AV18_Def.SetNormalizedSource(true);
  Kitty_AV18_Def.SetQ1Q2(0);
  Kitty_AV18_Def.KillTheCat();
  TH1F* hPS_AV18_Def = new TH1F("hPS_AV18_Def","hPS_AV18_Def",kSteps,kMin,kMax);
  TH1F* hER_AV18_Def = new TH1F("hER_AV18_Def","hER_AV18_Def",kSteps,kMin,kMax);
  for(unsigned uMom=0; uMom<kSteps; uMom++){
    double kstar = Kitty_AV18_Def.GetMomentum(uMom);
    hPS_AV18_Def->SetBinContent(uMom+1,Kitty_AV18_Def.GetPhaseShift(uMom,0,0));
    //hPS_AV18_Def->SetBinError(uMom+1, 0.01 * sqrt(kstar*0.5));
    hPS_AV18_Def->SetBinError(uMom+1, 0.001 + 0.0001 * kstar);

    double ps = Kitty_AV18_Def.GetPhaseShift(uMom,0,0);
    hER_AV18_Def->SetBinContent(uMom+1,kstar/tan(ps));
    //printf("%f %f %f\n",kstar,ps,kstar/tan(ps));
    double ps_err = kstar>=20?0.001*kstar*kstar/tan(ps):0.001* 20*kstar/tan(ps);
    hER_AV18_Def->SetBinError(uMom+1, ps_err);    
  }
  //Kitty_AV18_Def.SetQ1Q2(1);
  //Kitty_AV18_Def.KillTheCat();  



  CATS Kitty_Toy;
  Kitty_Toy.SetMomBins(kSteps,kMin,kMax);
  Kitty_Toy.SetAnaSource(GaussSource, cSorPars);
  Kitty_Toy.SetUseAnalyticSource(true);
  Kitty_Toy.SetMomentumDependentSource(false);
  Kitty_Toy.SetThetaDependentSource(false);
  Kitty_Toy.SetNormalizedSource(true);

  Kitty_Toy.SetExcludeFailedBins(false);
  Kitty_Toy.SetQ1Q2(0);
  Kitty_Toy.SetQuantumStatistics(true);
  Kitty_Toy.SetRedMass( 0.5*Mass_p );
  Kitty_Toy.SetNumChannels(2);
  Kitty_Toy.SetNumPW(0,1);
  Kitty_Toy.SetSpin(0,0);
  Kitty_Toy.SetSpin(1,1);
  Kitty_Toy.SetChannelWeight(0, 1./4.);
  Kitty_Toy.SetChannelWeight(1, 3./4.);
  CATSparameters cToyPars(CATSparameters::tPotential,10,true);
  cToyPars.SetParameter(0,3);
  Kitty_Toy.SetShortRangePotential(0,0,FloatingGaussians,cToyPars);


  TFile fOutput(OutputFile, "recreate");

  TNtuple* ntToyPot = new TNtuple("ntToyPot", "ntToyPot","SEED:ID:f0:d0:diff_av18:NumG:V0:mu0:r0:V1:mu1:r1:V2:mu2:r2");
  Float_t ntBuffer[15];

  TH1F* hPS_best = NULL;// = new TH1F(TString::Format("hPS_best_%u_%u",SEED,uIter),TString::Format("hPS_best_%u_%u",SEED,uIter),kSteps,kMin,kMax);
  TH1F* hER_best = NULL;// = new TH1F(TString::Format("hER_best_%u_%u",SEED,uIter),TString::Format("hER_best_%u_%u",SEED,uIter),kSteps,kMin,kMax);
  TF1* fER_best = NULL;// = new TF1(TString::Format("fER_best_%u_%u",SEED,uIter),"197.327*[0]+0.5*[1]/197.327*x*x",kMin_ER,kMax_ER);
  TH1F* hPot_best = NULL;// = new TH1F(TString::Format("hPot_best_%u_%u",SEED,uIter),TString::Format("hPot_best_%u_%u",SEED,uIter),2048,0,2);

  TRandom3 rangen(SEED);
  double best_diff_av18 = 1e6;
  for(unsigned uIter=0; uIter<NumIter; uIter++){

    TH1F* hPS = new TH1F(TString::Format("hPS_%u_%u",SEED,uIter),TString::Format("hPS_%u_%u",SEED,uIter),kSteps,kMin,kMax);
    TH1F* hER = new TH1F(TString::Format("hER_%u_%u",SEED,uIter),TString::Format("hER_%u_%u",SEED,uIter),kSteps,kMin,kMax);
    TF1* fER = new TF1(TString::Format("fER_%u_%u",SEED,uIter),"197.327*[0]+0.5*[1]/197.327*x*x",kMin_ER,kMax_ER);
    TH1F* hPot = new TH1F(TString::Format("hPot_%u_%u",SEED,uIter),TString::Format("hPot_%u_%u",SEED,uIter),2048,0,2);

    //toy1
    Kitty_Toy.SetShortRangePotential(0,0,1,rangen.Uniform(-10000,10000));
    Kitty_Toy.SetShortRangePotential(0,0,2,0);  
    Kitty_Toy.SetShortRangePotential(0,0,3,rangen.Uniform(0.01,0.1));
    //Kitty_Toy.SetShortRangePotential(0,0,1,rangen.Uniform(0));
    //Kitty_Toy.SetShortRangePotential(0,0,2,rangen.Uniform(0));  
    //Kitty_Toy.SetShortRangePotential(0,0,3,rangen.Uniform(1));

    Kitty_Toy.SetShortRangePotential(0,0,4,rangen.Uniform(-1000,1000));
    Kitty_Toy.SetShortRangePotential(0,0,5,rangen.Uniform(0.1,0.4));  
    Kitty_Toy.SetShortRangePotential(0,0,6,rangen.Uniform(0.1,0.4));
    //Kitty_Toy.SetShortRangePotential(0,0,4,527.05401);
    //Kitty_Toy.SetShortRangePotential(0,0,5,0.18620);  
    //Kitty_Toy.SetShortRangePotential(0,0,6,0.00734);


    Kitty_Toy.SetShortRangePotential(0,0,7,rangen.Uniform(-1000,1000));
    Kitty_Toy.SetShortRangePotential(0,0,8,rangen.Uniform(0.1,0.4));  
    Kitty_Toy.SetShortRangePotential(0,0,9,rangen.Uniform(0.1,0.4));
    //Kitty_Toy.SetShortRangePotential(0,0,7,rangen.Uniform(-686.79650));
    //Kitty_Toy.SetShortRangePotential(0,0,8,rangen.Uniform(0,0.35244));  
    //Kitty_Toy.SetShortRangePotential(0,0,9,rangen.Uniform(0,0.09530));    

    Kitty_Toy.SetNotifications(CATS::nWarning);
    Kitty_Toy.KillTheCat();


    double diff_av18 = 0;
    for(unsigned uMom=0; uMom<kSteps; uMom++){
      double ps = Kitty_Toy.GetPhaseShift(uMom,0,0);
      double kstar = Kitty_Toy.GetMomentum(uMom);
      //printf("%f %f %f\n",kstar,ps,Kitty_Toy.GetCorrFun(uMom));
      hER->SetBinContent(uMom+1, kstar/tan(ps));
      hER->SetBinError(uMom+1, 1+kstar*0.2);

      hPS->SetBinContent(uMom+1, ps);
      hPS->SetBinError(uMom+1, 0.01);

      if(kstar>=kMin_ER && kstar<=kMax_ER){
        diff_av18 += pow(ps-hPS_AV18_Def->GetBinContent(uMom+1),2.);
      }
    }
    fER->SetParameter(0, 0);
    fER->SetParLimits(0, -10,10);
    fER->SetParLimits(1, -50,150);

    hER->Fit(fER,"Q, S, N, R, M");\
    //if(1/fER->GetParameter(0)>10)

    for(unsigned uRad=0; uRad<hPot->GetNbinsX(); uRad++){
      float rad = hPot->GetBinCenter(uRad+1);
      hPot->SetBinContent(uRad+1, Kitty_Toy.EvaluateThePotential(0,0,10,rad));
    }


    //printf("%u: f; d = %.2f/%.2f; chi2=%.0f; diff=%.3f\n", uIter, 1/fER->GetParameter(0), fER->GetParameter(1), fER->GetChisquare(), diff_av18);
    //printf("    V0,R0,V1,R1 = %.3f, %.3f, %.3f, %.3f\n", Kitty_Toy.GetPotPar(0,0,0),Kitty_Toy.GetPotPar(0,0,1),Kitty_Toy.GetPotPar(0,0,2),Kitty_Toy.GetPotPar(0,0,3));

    //hPS->Write();
    //hER->Write();
    //fER->Write();
    //hPot->Write();

    if(fER->GetChisquare()<100){
      ntBuffer[0] = SEED;
      ntBuffer[1] = uIter;
      ntBuffer[2] = 1/fER->GetParameter(0);
      ntBuffer[3] = fER->GetParameter(1);
      ntBuffer[4] = diff_av18;
      ntBuffer[5] = Kitty_Toy.GetPotPar(0,0,0);
      ntBuffer[6] = Kitty_Toy.GetPotPar(0,0,1);
      ntBuffer[7] = Kitty_Toy.GetPotPar(0,0,2);
      ntBuffer[8] = Kitty_Toy.GetPotPar(0,0,3);
      ntBuffer[9] = Kitty_Toy.GetPotPar(0,0,4);
      ntBuffer[10] = Kitty_Toy.GetPotPar(0,0,5);
      ntBuffer[11] = Kitty_Toy.GetPotPar(0,0,6);
      ntBuffer[12] = Kitty_Toy.GetPotPar(0,0,7);
      ntBuffer[13] = Kitty_Toy.GetPotPar(0,0,8);
      ntBuffer[14] = Kitty_Toy.GetPotPar(0,0,9);


      ntToyPot->Fill(ntBuffer);

      if(diff_av18<best_diff_av18 && 1./fER->GetParameter(0)>2 && fER->GetParameter(1)>2.0 && fER->GetParameter(1)<3.6){

        printf("%u: f; d = %.2f/%.2f; chi2=%.0f; diff=%.3f\n", uIter, 1/fER->GetParameter(0), fER->GetParameter(1), fER->GetChisquare(), diff_av18);
        printf("    V0,M0,R0 = %.5f, %.5f, %.5f\n", Kitty_Toy.GetPotPar(0,0,1),Kitty_Toy.GetPotPar(0,0,2),Kitty_Toy.GetPotPar(0,0,3));
        printf("    V1,M1,R1 = %.5f, %.5f, %.5f\n", Kitty_Toy.GetPotPar(0,0,4),Kitty_Toy.GetPotPar(0,0,5),Kitty_Toy.GetPotPar(0,0,6));
        printf("    V2,M2,R2 = %.5f, %.5f, %.5f\n", Kitty_Toy.GetPotPar(0,0,7),Kitty_Toy.GetPotPar(0,0,8),Kitty_Toy.GetPotPar(0,0,9));

        //hPS->Write();
        //hER->Write();
        //fER->Write();
        //hPot->Write();
        if(hPS_best) delete hPS_best;
        if(hER_best) delete hER_best;
        if(fER_best) delete fER_best;
        if(hPot_best) delete hPot_best;

        hPS_best = (TH1F*)hPS->Clone("hPS_best");
        hER_best = (TH1F*)hER->Clone("hER_best");
        fER_best = (TF1*)fER->Clone("fER_best");
        hPot_best = (TH1F*)hPot->Clone("hPot_best");

        best_diff_av18 = diff_av18;
      }
    }

    delete hPS;
    delete hER;
    delete fER;
    delete hPot;
  }
  //Kitty_Toy.SetShortRangePotential(0,0,0,1000);
  //Kitty_Toy.SetShortRangePotential(0,0,1,0);  
  //Kitty_Toy.SetShortRangePotential(0,0,2,-540);
  //Kitty_Toy.SetShortRangePotential(0,0,3,0.8);
  hPS_AV18_Def->Write();
  if(hPS_best) hPS_best->Write();
  hER_AV18_Def->Write();
  if(hER_best) hER_best->Write();
  if(hPot_best) hPot_best->Write();
  ntToyPot->Write();




  delete ntToyPot;
}


void effect_of_square_barrier(){

  TString OutputFolder = TString::Format("%s/OtherTasks/effect_of_square_barrier/",GetFemtoOutputFolder());

  std::vector<float> positions = {0,0.1,0.2,0.4,0.8,1.6};
  std::vector<float> amplitudes = {-400,-200,-100,-50,-25};
  std::vector<float> half_widths = {0.001,0.1,0.2,0.4,0.8};


  unsigned kSteps = 16;
  double kMin = 0;
  double kMax = 80;
  double kMin_ER = 0;
  double kMax_ER = 60;

  double source_size = 1.2;

  CATSparameters cSorPars(CATSparameters::tSource,1,true);
  cSorPars.SetParameter(0,source_size);


  CATS Kitty_AV18_Def;
  Kitty_AV18_Def.SetMomBins(kSteps,kMin,kMax);
  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles/",GetCernBoxDimi()));
  AnalysisObject.SetUpCats_pp(Kitty_AV18_Def,"AV18_s","",0,0);
  Kitty_AV18_Def.SetAnaSource(GaussSource, cSorPars);
  Kitty_AV18_Def.SetUseAnalyticSource(true);
  Kitty_AV18_Def.SetMomentumDependentSource(false);
  Kitty_AV18_Def.SetThetaDependentSource(false);
  Kitty_AV18_Def.SetNormalizedSource(true);
  Kitty_AV18_Def.SetQ1Q2(0);
  Kitty_AV18_Def.KillTheCat();
  TH1F* hPS_AV18_Def = new TH1F("hPS_AV18_Def","hPS_AV18_Def",kSteps,kMin,kMax);
  TH1F* hER_AV18_Def = new TH1F("hER_AV18_Def","hER_AV18_Def",kSteps,kMin,kMax);
  for(unsigned uMom=0; uMom<kSteps; uMom++){
    double kstar = Kitty_AV18_Def.GetMomentum(uMom);
    hPS_AV18_Def->SetBinContent(uMom+1,Kitty_AV18_Def.GetPhaseShift(uMom,0,0));
    //hPS_AV18_Def->SetBinError(uMom+1, 0.01 * sqrt(kstar*0.5));
    hPS_AV18_Def->SetBinError(uMom+1, 0.001 + 0.0001 * kstar);

    double ps = Kitty_AV18_Def.GetPhaseShift(uMom,0,0);
    hER_AV18_Def->SetBinContent(uMom+1,kstar/tan(ps));
    //printf("%f %f %f\n",kstar,ps,kstar/tan(ps));
    double ps_err = kstar>=20?0.001*kstar*kstar/tan(ps):0.001* 20*kstar/tan(ps);
    hER_AV18_Def->SetBinError(uMom+1, ps_err);    
  }



  CATS Kitty_Toy;
  Kitty_Toy.SetMomBins(kSteps,kMin,kMax);
  Kitty_Toy.SetAnaSource(GaussSource, cSorPars);
  Kitty_Toy.SetUseAnalyticSource(true);
  Kitty_Toy.SetMomentumDependentSource(false);
  Kitty_Toy.SetThetaDependentSource(false);
  Kitty_Toy.SetNormalizedSource(true);

  Kitty_Toy.SetExcludeFailedBins(false);
  Kitty_Toy.SetQ1Q2(0);
  Kitty_Toy.SetQuantumStatistics(false);
  Kitty_Toy.SetRedMass( 0.5*Mass_p );
  Kitty_Toy.SetNumChannels(2);
  Kitty_Toy.SetNumPW(0,1);
  Kitty_Toy.SetSpin(0,0);
  Kitty_Toy.SetSpin(1,1);
  Kitty_Toy.SetChannelWeight(0, 1);
  Kitty_Toy.SetChannelWeight(1, 0);
  CATSparameters cToyPars(CATSparameters::tPotential,3,true);
  Kitty_Toy.SetShortRangePotential(0,0,PotentialBarrier,cToyPars);
  Kitty_Toy.SetNotifications(CATS::nWarning);

  std::vector<TH1F*> hPS;
  std::vector<TH1F*> hER;

  TFile fOutput(OutputFolder+"fOutput.root", "recreate");
  for(int iAmp=0; iAmp<amplitudes.size(); iAmp++){
    //printf("amp = %.0f\n",amplitudes.at(iAmp));
    Kitty_Toy.SetShortRangePotential(0,0,0,amplitudes.at(iAmp));
    for(int iPos=0; iPos<positions.size(); iPos++){
      //printf(" pos = %.2f\n",positions.at(iPos));
      Kitty_Toy.SetShortRangePotential(0,0,1,positions.at(iPos));
      for(int iHw=0; iHw<half_widths.size(); iHw++){
        //printf("  hw = %.3f\n",half_widths.at(iHw));
        Kitty_Toy.SetShortRangePotential(0,0,2,half_widths.at(iHw));
        Kitty_Toy.KillTheCat();
        TH1F* hPS_PB = new TH1F(TString::Format("hPS_PB_A%.0f_R%.3f_HW%.3f",amplitudes.at(iAmp),positions.at(iPos),half_widths.at(iHw)),TString::Format("hPS_PB_A%.0f_R%.3f_HW%.3f",amplitudes.at(iAmp),positions.at(iPos),half_widths.at(iHw)),kSteps,kMin,kMax);
        TH1F* hER_PB = new TH1F(TString::Format("hER_PB_A%.0f_R%.3f_HW%.3f",amplitudes.at(iAmp),positions.at(iPos),half_widths.at(iHw)),TString::Format("hER_PB_A%.0f_R%.3f_HW%.3f",amplitudes.at(iAmp),positions.at(iPos),half_widths.at(iHw)),kSteps,kMin,kMax);
        for(unsigned uMom=0; uMom<kSteps; uMom++){
          double kstar = Kitty_AV18_Def.GetMomentum(uMom);
          hPS_PB->SetBinContent(uMom+1,Kitty_Toy.GetPhaseShift(uMom,0,0));
          //hPS_AV18_Def->SetBinError(uMom+1, 0.01 * sqrt(kstar*0.5));
          //hPS_PB->SetBinError(uMom+1, 0.001 + 0.0001 * kstar);
          

          double ps = Kitty_Toy.GetPhaseShift(uMom,0,0);
          hER_PB->SetBinContent(uMom+1,kstar/tan(ps));
          if(uMom==0 && tan(ps)/kstar*hbarc>2){
            printf("amp = %.0f\n",amplitudes.at(iAmp));
            printf(" pos = %.2f\n",positions.at(iPos));
            printf("  hw = %.3f\n",half_widths.at(iHw));
            printf("   f=%.2f\n",tan(ps)/kstar*hbarc);
          }
            
          //printf("%f %f %f\n",kstar,ps,kstar/tan(ps));
          double ps_err = kstar>=20?0.001*kstar*kstar/tan(ps):0.001* 20*kstar/tan(ps);
          //hER_PB->SetBinError(uMom+1, ps_err);
         
        }
        hPS.push_back(hPS_PB);
        hER.push_back(hER_PB);
      }
    }
  }

  //TFile fOutput(OutputFolder+"effect_of_square_barrier.root","recreate");
  fOutput.cd();
  hPS_AV18_Def->Write();
  for (TH1F* histPS : hPS){
    histPS->Write();
  }
  hER_AV18_Def->Write();
  for (TH1F* histER : hER){
    histER->Write();
  }  

}


std::vector<int> get_available_slots(int* airport_status, int* gate_status, int slot_len){
  std::vector<int> list_of_available_slots;
  for(int iSlot=0; iSlot<144; iSlot++){
    bool status=true;
    //end of the day
    if(iSlot>144-slot_len) {status=false;}
    //runway blocked
    if(airport_status[iSlot]/1000==1) {status=false;}
    //taxi blocked at landing
    if( (airport_status[iSlot+1]/10)%10==1 ) {status=false;}
    //gate blocked at landing
    if( airport_status[iSlot+1]%10==1 || airport_status[iSlot+2]%10==1 ) {status=false;}
    int last_slot = iSlot+slot_len-1;
    int departure_slot = last_slot-3;
    int departure_len = 2;
    if(slot_len<=18) {departure_slot = last_slot-2; departure_len = 1;} //short hauls
    //gate blocked at departure
    if( airport_status[departure_slot]%10==1 ) {status=false;}
    //taxi blocked at departure
    if( (airport_status[departure_slot]/10)%10==1 ) {status=false;}
    if( departure_len==2 && (airport_status[departure_slot+1]/10)%10==1 ) {status=false;}
    //runway blocked at departure
    if( (airport_status[last_slot-1]/100)%10==1) {status=false;}
    if( airport_status[last_slot]/1000==1) {status=false;}
    for(int iSlot2=iSlot; iSlot2<iSlot+slot_len; iSlot2++){
      if( gate_status[iSlot2]!=0 ) {status=false; break;}
    }
    

    //if(iSlot==3&&status==false){
    //  printf("airport_status[iSlot] = %i\n",airport_status[iSlot]);
    //  printf("airport_status[iSlot+1] = %i\n",airport_status[iSlot+1]);
    //  printf("departure_slot = %i\n",departure_slot);
    //  printf("departure_len = %i\n",departure_len);
    //  printf("airport_status[departure_slot] = %i\n",airport_status[departure_slot]);
    //  printf("airport_status[departure_slot+1] = %i\n",airport_status[departure_slot+1]);
    //  printf("last_slot = %i\n",last_slot);
    //  printf("airport_status[last_slot-1] = %i\n",airport_status[last_slot-1]);
    //  printf("airport_status[last_slot] = %i\n",airport_status[last_slot]);     
    //}

    if(status){
      list_of_available_slots.push_back(iSlot);
    }
  }
  return list_of_available_slots;
}

void AirSim_Schedule(int SEED, int nIter, int nSH, int nMH, int nLH){

  const double rating_SH = 1;
  const double rating_MH = 2;
  const double rating_LH = 3;
  const int slot_len_SH = 18;
  const int slot_len_MH = 24;
  const int slot_len_LH = 36;


  double best_rating = 0;

  TRandom3 rangen(SEED);
  for(int iIter=0; iIter<nIter; iIter++){

    int nGate = nSH+nMH+nLH;
    //slot values:
    //0 = available (for the gate -1 is booked, currently at parked)
    //!0 = booked
    //ABCD with 
    //  A runway status: 0 free, 1 blocled
    //  B taxi to runway: 0 free, 1 blocled 
    //  C taxi status: 0 free, 1 blocked (at the cross road)
    //  D terminal status: 0 free, 1 blocked
    //The rule is: 
    int** gate = new int* [nGate];
    int* airport_status = new int [144];
    for(int iSlot=0; iSlot<144; iSlot++){
      airport_status[iSlot] = 0;
    }
    for(int iGate=0; iGate<nGate; iGate++){
      gate[iGate] = new int[144];//24*6
      for(int iSlot=0; iSlot<144; iSlot++){
        gate[iGate][iSlot] = 0;
      }
    }

    const int max_dev = 4;
    double total_rating = 0;
    for(int iLH=0; iLH<nLH; iLH++){
      int counter=0;
      while(true){
        std::vector<int> available_slots = get_available_slots(airport_status, gate[iLH], slot_len_LH);
        //for(int iav=0; iav<available_slots.size(); iav++){
        //  printf("%i\n",available_slots.at(iav));
        //}
        if(available_slots.size()==0) break;
        int new_slot = available_slots.at(rangen.Integer(available_slots.size()<max_dev?available_slots.size():max_dev));
        //int new_slot = available_slots.at(0);
        total_rating += rating_LH;
        gate[iLH][new_slot] = 1000;
        gate[iLH][new_slot+1] = 11;
        for(int iSlot=new_slot+2; iSlot<new_slot+32; iSlot++){
          gate[iLH][iSlot] = -1;
        }
        gate[iLH][new_slot+32] = 11;
        gate[iLH][new_slot+33] = 10;
        gate[iLH][new_slot+34] = 100;
        gate[iLH][new_slot+35] = 1000;

        airport_status[new_slot] = airport_status[new_slot]%1000 + 1000;
        airport_status[new_slot+1] = (airport_status[new_slot+1]/100)*100 + 11;
        airport_status[new_slot+32] = (airport_status[new_slot+32]/100)*100 + 11;
        airport_status[new_slot+33] = airport_status[new_slot+33] - ((airport_status[new_slot+33]/10)%10)*10 + 10;//change the second digit
        airport_status[new_slot+34] = airport_status[new_slot+34] - ((airport_status[new_slot+34]/100)%10)*100 + 100;//change the third digit
        airport_status[new_slot+35] = airport_status[new_slot+35]%1000 + 1000;
        //printf("new_slot = %i\n",new_slot);
        counter++;
//if(counter==2) break;
      }
    }

    for(int iMH=0; iMH<nMH; iMH++){
      int counter=0;
      while(true){
        std::vector<int> available_slots = get_available_slots(airport_status, gate[iMH], slot_len_MH);
        if(available_slots.size()==0) break;
        int new_slot = available_slots.at(rangen.Integer(available_slots.size()<max_dev?available_slots.size():max_dev));
        total_rating += rating_MH;
        gate[iMH][new_slot] = 1000;
        gate[iMH][new_slot+1] = 11;
        for(int iSlot=new_slot+2; iSlot<new_slot+20; iSlot++){
          gate[iMH][iSlot] = -1;
        }
        gate[iMH][new_slot+20] = 11;
        gate[iMH][new_slot+21] = 10;
        gate[iMH][new_slot+22] = 100;
        gate[iMH][new_slot+23] = 1000;

        airport_status[new_slot] = airport_status[new_slot]%1000 + 1000;
        airport_status[new_slot+1] = (airport_status[new_slot+1]/100)*100 + 11;
        airport_status[new_slot+20] = (airport_status[new_slot+20]/100)*100 + 11;
        airport_status[new_slot+21] = airport_status[new_slot+21] - ((airport_status[new_slot+21]/10)%10)*10 + 10;//change the second digit
        airport_status[new_slot+22] = airport_status[new_slot+22] - ((airport_status[new_slot+22]/100)%10)*100 + 100;//change the third digit
        airport_status[new_slot+23] = airport_status[new_slot+23]%1000 + 1000;
        counter++;
      }
    }

    for(int iSH=0; iSH<nSH; iSH++){
      int counter=0;
      while(true){
        std::vector<int> available_slots = get_available_slots(airport_status, gate[iSH], slot_len_SH);
        if(available_slots.size()==0) break;
        int new_slot = available_slots.at(rangen.Integer(available_slots.size()<max_dev?available_slots.size():max_dev));
        total_rating += rating_SH;
        gate[iSH][new_slot] = 1000;
        gate[iSH][new_slot+1] = 11;
        for(int iSlot=new_slot+2; iSlot<new_slot+15; iSlot++){
          gate[iSH][iSlot] = -1;
        }
        gate[iSH][new_slot+15] = 11;
        gate[iSH][new_slot+16] = 100;
        gate[iSH][new_slot+17] = 1000;

        airport_status[new_slot] = airport_status[new_slot]%1000 + 1000;
        airport_status[new_slot+1] = (airport_status[new_slot+1]/100)*100 + 11;
        airport_status[new_slot+15] = (airport_status[new_slot+15]/100)*100 + 11;
        airport_status[new_slot+16] = airport_status[new_slot+16] - ((airport_status[new_slot+16]/100)%10)*100 + 100;//change the third digit
        airport_status[new_slot+17] = airport_status[new_slot+17]%1000 + 1000;
        counter++;
      }
    }

    
    if(total_rating>best_rating){
      best_rating = total_rating;
      printf("best_rating = %f\n",best_rating);
    }
    for(int iGate=0; iGate<nGate; iGate++){
      delete [] gate[iGate];
    }
    delete [] gate;
    delete [] airport_status;
  }
  
}


//what I leaned
//for pp the position of the maximum goes up (in k*, amplitude actually goes down for same f0) as the eff range goes down
//in any case, it will not work by simply shifting the WF, as the division by 1/r kind of makes any inperfections at low kstar to explode
//so probably, I will only get e "perfect" fit in the effective range zero approximation
void Sofia_test1(){
  TString OutputFolder = TString::Format("%s/OtherTasks/Sofia_test1/",GetFemtoOutputFolder());

  unsigned NumMomBins = 60;
  double kMin = 0;
  double kMax = 240;
  double source_size = 1.8;
  double cut_off = 0.8;


  CATSparameters cSorPars_Castrate(CATSparameters::tSource,3,true);
  cSorPars_Castrate.SetParameter(0,source_size);
  cSorPars_Castrate.SetParameter(1,cut_off);
  cSorPars_Castrate.SetParameter(2,0.01);

  CATS Nala;
  Nala.SetMomBins(NumMomBins,kMin,kMax);
  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles/",GetCernBoxDimi()));
  AnalysisObject.SetUpCats_pp(Nala,"AV18_s","",0,0);
  Nala.SetAnaSource(GaussSourceGentleCutOff, cSorPars_Castrate);
  Nala.SetUseAnalyticSource(true);
  Nala.SetMomentumDependentSource(false);
  Nala.SetThetaDependentSource(false);
  Nala.SetNormalizedSource(false);
  Nala.SetAutoNormSource(true);
  //Nala.SetAnaSource(0,source_size);

  //Nala.SetThetaDependentSource(true);
  Nala.ComputeTheSource();
  printf("AV18\n");
  Nala.KillTheCat(); 

  TFile fOutput(OutputFolder+"SOF.root", "recreate");
  TGraph* gComputedSource_Normal = GetSourceFromCATS_TGraph(Nala, 0, 0, "gComputedSource_Normal");

  TF1* fSource = new TF1("fSource",GaussSourceTF1,0,16,1);
  fSource->SetParameter(0,source_size);

  fOutput.cd();
  TGraph g_AV18;
  g_AV18.SetName("g_AV18");
  g_AV18.SetLineWidth(6);
  g_AV18.SetLineColor(kBlue+1);
  for(unsigned uMom=0; uMom<NumMomBins; uMom++){
    double kstar = Nala.GetMomentum(uMom);
    g_AV18.SetPoint(uMom, kstar, Nala.GetCorrFun(uMom));
  }

  Nala.SetShortRangeSofia(cut_off);
  printf("AV18 SOF\n");
  Nala.KillTheCat(); 

  fOutput.cd();
  TGraph* gComputedSource_CutOff = GetSourceFromCATS_TGraph(Nala, 0, 0, "gComputedSource_CutOff");
  fOutput.cd();
  TGraph g_SOF;
  g_SOF.SetName("g_SOF");
  g_SOF.SetLineWidth(4);
  g_SOF.SetLineColor(kGray);
  for(unsigned uMom=0; uMom<NumMomBins; uMom++){
    double kstar = Nala.GetMomentum(uMom);
    g_SOF.SetPoint(uMom, kstar, Nala.GetCorrFun(uMom));
  }  



  CATS NalaG;
  NalaG.SetMomBins(NumMomBins,kMin,kMax);
  CATSparameters sPars(CATSparameters::tSource,1,true);
  //NalaG.SetAnaSource(GaussSource, sPars);
  //NalaG.SetUseAnalyticSource(true);
  //NalaG.SetAnaSource(0,source_size);
  NalaG.SetAnaSource(GaussSourceGentleCutOff, cSorPars_Castrate);
  NalaG.SetUseAnalyticSource(true);
  NalaG.SetMomentumDependentSource(false);
  NalaG.SetThetaDependentSource(false);
  NalaG.SetNormalizedSource(false);
  NalaG.SetAutoNormSource(true);

  NalaG.SetRedMass(0.5*Mass_p);
  NalaG.SetQ1Q2(1);
  NalaG.SetQuantumStatistics(true);
  NalaG.SetNumChannels(2);
  NalaG.SetSpin(0,0);
  NalaG.SetSpin(1,1);
  NalaG.SetNumPW(0,1);
  NalaG.SetNumPW(1,0);
  NalaG.SetChannelWeight(0, 0.25);
  NalaG.SetChannelWeight(1, 0.75);
  double PotPars1S0[2]={-29.32577,1.841348};
  CATSparameters cPotPars1S0(CATSparameters::tPotential,2,true);
  cPotPars1S0.SetParameters(PotPars1S0);
  NalaG.SetShortRangePotential(0,0,SingleGauss,cPotPars1S0);

  printf("Gauss to AV18\n");
  NalaG.KillTheCat(); 

  NalaG.SetShortRangePotential(0,0,0,-29.32577);
  NalaG.SetShortRangePotential(0,0,1,1.841348);
  NalaG.KillTheCat(); 

  fOutput.cd();
  TGraph g_Gauss;
  g_Gauss.SetName("g_Gauss");
  g_Gauss.SetLineWidth(6);
  g_Gauss.SetLineColor(kRed);
  for(unsigned uMom=0; uMom<NumMomBins; uMom++){
    double kstar = Nala.GetMomentum(uMom);
    g_Gauss.SetPoint(uMom, kstar, NalaG.GetCorrFun(uMom));
  }


  //small ~0.8fm effective range
  NalaG.SetShortRangePotential(0,0,0,-360.376357);
  NalaG.SetShortRangePotential(0,0,1,0.5463289867);
  //NalaG.SetAnaSource(0,source_size+0.25);
  printf("0.8 ER\n");
  NalaG.KillTheCat(); 
  

  fOutput.cd();
  TGraph g_Gauss_smallER;
  g_Gauss_smallER.SetName("g_Gauss_smallER");
  g_Gauss_smallER.SetLineWidth(6);
  g_Gauss_smallER.SetLineColor(kRed-2);
  for(unsigned uMom=0; uMom<NumMomBins; uMom++){
    double kstar = Nala.GetMomentum(uMom);
    g_Gauss_smallER.SetPoint(uMom, kstar, NalaG.GetCorrFun(uMom));
  }

  //zero ~0.14fm effective range
  NalaG.SetShortRangePotential(0,0,0,-11094.11455);
  NalaG.SetShortRangePotential(0,0,1,0.09987347);
  //NalaG.SetAnaSource(0,source_size+0.25);
  printf("0.15 ER\n");
  NalaG.KillTheCat(); 

  fOutput.cd();
  TGraph g_Gauss_zeroER;
  g_Gauss_zeroER.SetName("g_Gauss_zeroER");
  g_Gauss_zeroER.SetLineWidth(6);
  g_Gauss_zeroER.SetLineColor(kCyan-1);
  for(unsigned uMom=0; uMom<NumMomBins; uMom++){
    double kstar = Nala.GetMomentum(uMom);
    g_Gauss_zeroER.SetPoint(uMom, kstar, NalaG.GetCorrFun(uMom));
  }


  //NalaG.SetAnaSource(0,source_size);


  printf("SOF\n");
  NalaG.SetShortRangeSofia(cut_off);
  NalaG.KillTheCat(); 


  fOutput.cd();
  TGraph g_Gauss_SOF;
  g_Gauss_SOF.SetName("g_Gauss_SOF");
  g_Gauss_SOF.SetLineWidth(4);
  g_Gauss_SOF.SetLineColor(kRed+3);
  for(unsigned uMom=0; uMom<NumMomBins; uMom++){
    double kstar = Nala.GetMomentum(uMom);
    g_Gauss_SOF.SetPoint(uMom, kstar, NalaG.GetCorrFun(uMom));
  }  


  fOutput.cd();
  g_AV18.Write();
  g_Gauss_smallER.Write();
  g_Gauss_zeroER.Write();
  g_Gauss.Write();
  g_SOF.Write();
  g_Gauss_SOF.Write();
  fSource->Write();
  gComputedSource_Normal->Write();
  gComputedSource_CutOff->Write();

  delete fSource;
  gROOT->cd();
  delete gComputedSource_Normal;
  delete gComputedSource_CutOff;
  
}

//two gaussians of exactly the same scattering pars, but rather different range parameters
//BOTH are SHORT RANGED though!!
void compare_pot_same_shift_diff_range(){
  /*
  CATS NalaG;
  NalaG.SetMomBins(NumMomBins,kMin,kMax);
  CATSparameters sPars(CATSparameters::tSource,1,true);
  //NalaG.SetAnaSource(GaussSource, sPars);
  //NalaG.SetUseAnalyticSource(true);
  //NalaG.SetAnaSource(0,source_size);
  NalaG.SetAnaSource(GaussSourceGentleCutOff, cSorPars_Castrate);
  NalaG.SetUseAnalyticSource(true);
  NalaG.SetMomentumDependentSource(false);
  NalaG.SetThetaDependentSource(false);
  NalaG.SetNormalizedSource(false);
  NalaG.SetAutoNormSource(true);

  NalaG.SetRedMass(0.5*Mass_p);
  NalaG.SetQ1Q2(1);
  NalaG.SetQuantumStatistics(true);
  NalaG.SetNumChannels(2);
  NalaG.SetSpin(0,0);
  NalaG.SetSpin(1,1);
  NalaG.SetNumPW(0,1);
  NalaG.SetNumPW(1,0);
  NalaG.SetChannelWeight(0, 0.25);
  NalaG.SetChannelWeight(1, 0.75);
  double PotPars1S0[2]={-29.32577,1.841348};
  CATSparameters cPotPars1S0(CATSparameters::tPotential,2,true);
  cPotPars1S0.SetParameters(PotPars1S0);
  NalaG.SetShortRangePotential(0,0,SingleGauss,cPotPars1S0);

  printf("Gauss to AV18\n");
  NalaG.KillTheCat(); 
  */
}


void get_eta_fractions(){
  TString InputFileName = TString::Format("%s/CatsFiles/Source/CECA/eta_reso/eta_reso_list.txt",GetCernBoxDimi());
  //stuff with larger ctau will be counted as a lambda par of flat contrib
  const float upper_tau = 10;
  //stuff with smaller fraction will be added into one single effective resonance
  const float upper_fraction = 0.1 / 100.;

  

}






void ThreeBodyDecay(unsigned NumIter, const int SEED){

  //https://link.aps.org/accepted/10.1103/PhysRevLett.114.013003
  //we have 3H -> 3He(+) + e + antineutr
  //so the He is missing 1 e !!!
  double MASS = 2809.4321;//mother
  double mass[3] = {2809.4135-0.51099895,0.51099895,0.0};
  std::vector<double> vec_mass;
  //vec_mass.push_back(1);
  //vec_mass.push_back(1);
  vec_mass.push_back(mass[0]);
  vec_mass.push_back(mass[1]);
  vec_mass.push_back(mass[2]);
  NumIter *= 1000;
  double available_energy = MASS - mass[0] - mass[1] - mass[2];
  const double epsilon = 1e-6;
  printf("available_energy = %.6f\n",available_energy);
  printf(" for the e = %.6f\n",available_energy-mass[1]);

  TRandom3 rangen(SEED);
  double cos_th, sin_th, eta, phi, ptot, px, py, pz, ekin;
  TLorentzVector particle[3];
  TLorentzVector total;

  TFile fOutput(TString::Format("%s/OtherTasks/ThreeBodyDecay/fOutput_%iK_S%i.root",GetFemtoOutputFolder(),NumIter/1000,SEED), "recreate");
  TH1F* electon_spectrum = new TH1F("electon_spectrum","electon_spectrum",1024,0,250);
  TH1F* electon_ekin = new TH1F("electon_ekin","electon_ekin",1024,0,50);
  TH1F* electon_spectrum_cos_th = new TH1F("electon_spectrum_cos_th","electon_spectrum_cos_th",1024,-1,1);

  for(unsigned uIter=0; uIter<NumIter/100; uIter++){
    do{
      //tot_engy = 0;
      total.SetPxPyPzE(0,0,0,0);
      for(unsigned uDaughter=0; uDaughter<3; uDaughter++){
          cos_th = rangen.Uniform(-1,1);
          phi = rangen.Uniform(0,TMath::Pi()*2);
          //ekin = rangen.Uniform(0,available_energy);
          //ptot = rangen.Uniform(0,available_energy);
          //ptot = sqrt(ekin*ekin+2*ekin*mass[uDaughter]);
          ptot = rangen.Uniform(0,sqrt(available_energy*available_energy+2*available_energy*mass[uDaughter]));

          //ptot = fabs(rangen.Gaus(0,available_energy/1.5));
          eta = -0.5*log((1.-cos_th)/(1.+cos_th));
          //sin_th = 2.*exp(-eta)/(1.+exp(-2.*eta));
          sin_th = sqrt(1.-cos_th*cos_th);//okay, as theta is 0 - 180
          px = ptot*cos(phi)*sin_th;
          py = ptot*sin(phi)*sin_th;
          pz = ptot*cos_th;
          //printf(" d%i ptot = %f\n",uDaughter, ptot);
          particle[uDaughter].SetPxPyPzE(px,py,pz,sqrt(ptot*ptot+mass[uDaughter]*mass[uDaughter]));
          //tot_engy += particle[uDaughter].M2();
          total += particle[uDaughter];
        }
       //tot_engy = sqrt(tot_engy);
       //printf("tot_engy = %f; MASS = %f\n",total.M(),MASS);
       //usleep(40e3);
    }
    while(fabs(MASS-total.M())>epsilon);
    if(uIter%1000==0){
      printf("uIter = %u K \n",uIter/1000);
    }
    
    electon_spectrum->Fill( particle[1].P()*1000. );
    electon_spectrum_cos_th->Fill( particle[1].CosTheta() );
    electon_ekin->Fill( (particle[1].E()-particle[1].M())*1000 );
  }



  TH1F* electon_spectrum_dlm = new TH1F("electon_spectrum_dlm","electon_spectrum_dlm",1024,0,250);
  TH1F* electon_ekin_dlm = new TH1F("electon_ekin_dlm","electon_ekin_dlm",1024,0,50);
  DLM_Random dlm_rangen(11);
  for(unsigned uIter=0; uIter<NumIter; uIter++){
    
    CatsParticle cp_triton;
    cp_triton.SetMass(MASS);
    cp_triton.SetWidth(1e-30);
    cp_triton.SetDecayRanGen(dlm_rangen);
    CatsParticle* daughters = cp_triton.DecayN(vec_mass,true);

    double toteng=daughters[0].GetE()+daughters[1].GetE()+daughters[2].GetE();
    //printf("M=%f; toteng=%f\n",MASS,toteng);
    //printf(" daughters[0].GetMass() = %f; p = %f; e = %f\n",daughters[0].GetMass(), daughters[0].GetP(), daughters[0].GetE());
    //printf(" daughters[1].GetMass() = %f; p = %f; e = %f\n",daughters[1].GetMass(), daughters[1].GetP(), daughters[1].GetE());
    //printf(" daughters[2].GetMass() = %f; p = %f; e = %f\n",daughters[2].GetMass(), daughters[2].GetP(), daughters[2].GetE());
    electon_spectrum_dlm->Fill(daughters[1].GetP()*1000);
    electon_ekin_dlm->Fill((daughters[1].GetE()-daughters[1].GetMass())*1000);

    //printf("%u %e\n",uIter,daughters[1].GetP());
    delete [] daughters;
  }


  fOutput.cd();
  electon_spectrum->Write();
  electon_spectrum_cos_th->Write();
  electon_ekin->Write();

  electon_spectrum_dlm->Write();
  electon_ekin_dlm->Write();
  
  delete electon_spectrum;
  delete electon_spectrum_dlm;
}



/*
void ThreeBodyDecay(unsigned NumIter, const int SEED){

  double MASS = 2809.4321;//mother
  double mass[3] = {2809.4135-0.51099895,0.51099895,0.0};
  std::vector<double> vec_mass;
  //vec_mass.push_back(1);
  //vec_mass.push_back(1);
  vec_mass.push_back(mass[0]);
  vec_mass.push_back(mass[1]);
  vec_mass.push_back(mass[2]);
  NumIter *= 1000;
  double available_energy = MASS - mass[0] - mass[1] - mass[2];
  const double epsilon = 1e-6;
  printf("available_energy = %.6f\n",available_energy);
  printf(" for the e = %.6f\n",available_energy-mass[1]);

  TRandom3 rangen(SEED);
  double cos_th, sin_th, eta, phi, ptot, px, py, pz, ekin;
  TLorentzVector particle[3];
  TLorentzVector total;

  TFile fOutput(TString::Format("%s/OtherTasks/ThreeBodyDecay/fOutput_%iK_S%i.root",GetFemtoOutputFolder(),NumIter/1000,SEED), "recreate");
  TH1F* electon_spectrum = new TH1F("electon_spectrum","electon_spectrum",1024,0,250);
  TH1F* electon_ekin = new TH1F("electon_ekin","electon_ekin",1024,0,50);
  TH1F* electon_spectrum_cos_th = new TH1F("electon_spectrum_cos_th","electon_spectrum_cos_th",1024,-1,1);

  for(unsigned uIter=0; uIter<NumIter/100; uIter++){
    do{
      //tot_engy = 0;
      total.SetPxPyPzE(0,0,0,0);
      for(unsigned uDaughter=0; uDaughter<3; uDaughter++){
          cos_th = rangen.Uniform(-1,1);
          phi = rangen.Uniform(0,TMath::Pi()*2);
          //ekin = rangen.Uniform(0,available_energy);
          //ptot = rangen.Uniform(0,available_energy);
          //ptot = sqrt(ekin*ekin+2*ekin*mass[uDaughter]);
          ptot = rangen.Uniform(0,sqrt(available_energy*available_energy+2*available_energy*mass[uDaughter]));

          //ptot = fabs(rangen.Gaus(0,available_energy/1.5));
          eta = -0.5*log((1.-cos_th)/(1.+cos_th));
          //sin_th = 2.*exp(-eta)/(1.+exp(-2.*eta));
          sin_th = sqrt(1.-cos_th*cos_th);//okay, as theta is 0 - 180
          px = ptot*cos(phi)*sin_th;
          py = ptot*sin(phi)*sin_th;
          pz = ptot*cos_th;
          //printf(" d%i ptot = %f\n",uDaughter, ptot);
          particle[uDaughter].SetPxPyPzE(px,py,pz,sqrt(ptot*ptot+mass[uDaughter]*mass[uDaughter]));
          //tot_engy += particle[uDaughter].M2();
          total += particle[uDaughter];
        }
       //tot_engy = sqrt(tot_engy);
       //printf("tot_engy = %f; MASS = %f\n",total.M(),MASS);
       //usleep(40e3);
    }
    while(fabs(MASS-total.M())>epsilon);
    if(uIter%1000==0){
      printf("uIter = %u K \n",uIter/1000);
    }
    
    electon_spectrum->Fill( particle[1].P()*1000. );
    electon_spectrum_cos_th->Fill( particle[1].CosTheta() );
    electon_ekin->Fill( (particle[1].E()-particle[1].M())*1000 );
  }



  TH1F* electon_spectrum_dlm = new TH1F("electon_spectrum_dlm","electon_spectrum_dlm",1024,0,250);
  TH1F* electon_ekin_dlm = new TH1F("electon_ekin_dlm","electon_ekin_dlm",1024,0,50);
  DLM_Random dlm_rangen(11);
  for(unsigned uIter=0; uIter<NumIter; uIter++){
    
    CatsParticle cp_triton;
    cp_triton.SetMass(MASS);
    cp_triton.SetWidth(1e-30);
    cp_triton.SetDecayRanGen(dlm_rangen);
    CatsParticle* daughters = cp_triton.DecayN(vec_mass,true);

    double toteng=daughters[0].GetE()+daughters[1].GetE()+daughters[2].GetE();
    //printf("M=%f; toteng=%f\n",MASS,toteng);
    //printf(" daughters[0].GetMass() = %f; p = %f; e = %f\n",daughters[0].GetMass(), daughters[0].GetP(), daughters[0].GetE());
    //printf(" daughters[1].GetMass() = %f; p = %f; e = %f\n",daughters[1].GetMass(), daughters[1].GetP(), daughters[1].GetE());
    //printf(" daughters[2].GetMass() = %f; p = %f; e = %f\n",daughters[2].GetMass(), daughters[2].GetP(), daughters[2].GetE());
    electon_spectrum_dlm->Fill(daughters[1].GetP()*1000);
    electon_ekin_dlm->Fill((daughters[1].GetE()-daughters[1].GetMass())*1000);

    //printf("%u %e\n",uIter,daughters[1].GetP());
    delete [] daughters;
  }


  fOutput.cd();
  electon_spectrum->Write();
  electon_spectrum_cos_th->Write();
  electon_ekin->Write();

  electon_spectrum_dlm->Write();
  electon_ekin_dlm->Write();
  
  delete electon_spectrum;
  delete electon_spectrum_dlm;
}
*/

void pL_DanielPolish_test(){
  double source_size = 1.0;

  unsigned NumMomBins = 100;
  double kMin = 0;
  double kMax = 125;


  CATSparameters cSorPars(CATSparameters::tSource,1,true);
  cSorPars.SetParameter(0,source_size);

  CATS Nala;
  Nala.SetMomBins(NumMomBins,kMin,kMax);
  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles/",GetCernBoxDimi()));
  AnalysisObject.SetUpCats_pL(Nala,"Usmani","",0,0);
  Nala.SetAnaSource(GaussSource, cSorPars);
  Nala.SetUseAnalyticSource(true);
  Nala.SetMomentumDependentSource(false);
  Nala.SetThetaDependentSource(false);
  Nala.SetNormalizedSource(false);
  Nala.SetAutoNormSource(true);
  Nala.SetChannelWeight(0,1);
  Nala.SetChannelWeight(1,0);
  Nala.KillTheCat(); 

  for(unsigned uMom=0; uMom<NumMomBins; uMom++){
    printf("%.3f   %.3f\n",Nala.GetMomentum(uMom),Nala.GetCorrFun(uMom));
  }

}



void PlaneTimeSlots(){
  

  const int NumIter = 1;

  TRandom3 rangen(23);

  const int NumDest = 15;
  const int NumDays = 7;
  const int NumSlots = NumDays*24*4;
  TString* Destination = new TString [NumDest];
  //if 0 at the hub, else id of the plane counting from 1
  //the departure slot is marked with minus plane id
  int** RouteSlots = new int* [NumDest];
  //if true we depart at that slot
  bool** DeparSlots = new bool* [NumDest];
  //travel time in num slots
  int* TravelSlots = new int [NumDest];
  for(int iDest=0; iDest<NumDest; iDest++){
    RouteSlots[iDest] = new int [NumSlots];
    DeparSlots[iDest] = new bool [NumSlots];
    for(int iSlot=0; iSlot<NumSlots; iSlot++){
      RouteSlots[iDest][iSlot] = 0;
      DeparSlots[iDest][iSlot] = 0;
    }
  }

  const int NumAircraft = 19;
  int* AircraftID = new int [NumAircraft];
  for(int iPlane = 0; iPlane<NumAircraft; iPlane++){
    AircraftID[iPlane] = iPlane + 1;
  }


  
  int di = 0;
  int dep_hour;
  int dep_mins;
  int trv_hour;
  int trv_mins;


  for(int iDay = 0; iDay<NumDays; iDay++){
    dep_hour = 17;
    dep_mins = 0;
    DeparSlots[di][iDay*24*4 + dep_hour*4 + dep_mins/15] = true;
  }
  trv_hour = 31;
  trv_mins = 0;
  TravelSlots[di] = trv_hour*4 + trv_mins%15;
  Destination[di++] = "Honolulu";//0



  for(int iDay = 0; iDay<NumDays; iDay++){
    dep_hour = 10;
    dep_mins = 0;
    DeparSlots[di][iDay*24*4 + dep_hour*4 + dep_mins/15] = true;
  }
  trv_hour = 26;
  trv_mins = 0;
  TravelSlots[di] = trv_hour*4 + trv_mins%15;
  Destination[di++] = "Los Angeles";//1




  for(int iDay = 0; iDay<NumDays; iDay++){
    dep_hour = 15;
    dep_mins = 15;
    DeparSlots[di][iDay*24*4 + dep_hour*4 + dep_mins/15] = true;
  }
  trv_hour = 25;
  trv_mins = 30;  
  TravelSlots[di] = trv_hour*4 + trv_mins%15;
  Destination[di++] = "San Francisco";//2


  for(int iDay = 0; iDay<NumDays; iDay++){
    dep_hour = 13;
    dep_mins = 00;
    DeparSlots[di][iDay*24*4 + dep_hour*4 + dep_mins/15] = true;
  }
  trv_hour = 26;
  trv_mins = 45;    
  TravelSlots[di] = trv_hour*4 + trv_mins%15;
  Destination[di++] = "Mexico City";//3


  for(int iDay = 0; iDay<NumDays; iDay++){
    if(iDay%2==0){
      dep_hour = 9;
      dep_mins = 45;
    }
    else{
      dep_hour = 4;
      dep_mins = 30;      
    }
    DeparSlots[di][iDay*24*4 + dep_hour*4 + dep_mins/15] = true;
  }
  trv_hour = 25;
  trv_mins = 15;    
  TravelSlots[di] = trv_hour*4 + trv_mins%15;
  Destination[di++] = "Bogotá";//4



  for(int iDay = 0; iDay<NumDays; iDay++){
    if(iDay%2==0){
      dep_hour = 4;
      dep_mins = 45;
    }
    else{
      dep_hour = 12;
      dep_mins = 00;      
    }
    DeparSlots[di][iDay*24*4 + dep_hour*4 + dep_mins/15] = true;
  }
  trv_hour = 28;
  trv_mins = 00;    
  TravelSlots[di] = trv_hour*4 + trv_mins%15;  
  Destination[di++] = "Lima";//5

  for(int iDay = 0; iDay<NumDays; iDay++){
    dep_hour = 9;
    dep_mins = 30;
    DeparSlots[di][iDay*24*4 + dep_hour*4 + dep_mins/15] = true;
  }
  trv_hour = 25;
  trv_mins = 00;    
  TravelSlots[di] = trv_hour*4 + trv_mins%15;
  Destination[di++] = "São Paulo";//6

  for(int iDay = 0; iDay<NumDays; iDay++){
    dep_hour = 13;
    dep_mins = 15;
    DeparSlots[di][iDay*24*4 + dep_hour*4 + dep_mins/15] = true;
  }
  trv_hour = 36;
  trv_mins = 45;    
  TravelSlots[di] = trv_hour*4 + trv_mins%15;
  Destination[di++] = "Sydney";//7


  for(int iDay = 0; iDay<NumDays; iDay++){
    dep_hour = 13;
    dep_mins = 00;
    DeparSlots[di][iDay*24*4 + dep_hour*4 + dep_mins/15] = true;
  }
  trv_hour = 36;
  trv_mins = 00;    
  TravelSlots[di] = trv_hour*4 + trv_mins%15; 
  Destination[di++] = "Melbourne";//8


  for(int iDay = 0; iDay<NumDays; iDay++){
    dep_hour = 15;
    dep_mins = 00;
    DeparSlots[di][iDay*24*4 + dep_hour*4 + dep_mins/15] = true;
  }
  trv_hour = 41;
  trv_mins = 30;    
  TravelSlots[di] = trv_hour*4 + trv_mins%15; 
  Destination[di++] = "Aukland";//9

  for(int iDay = 0; iDay<NumDays; iDay++){
    dep_hour = 19;
    dep_mins = 30;
    DeparSlots[di][iDay*24*4 + dep_hour*4 + dep_mins/15] = true;
  }
  trv_hour = 14;
  trv_mins = 00; 
  TravelSlots[di] = trv_hour*4 + trv_mins%15;  
  Destination[di++] = "Mumbai";//10

  for(int iDay = 0; iDay<NumDays; iDay++){
    dep_hour = 20;
    dep_mins = 15;
    DeparSlots[di][iDay*24*4 + dep_hour*4 + dep_mins/15] = true;
  }
  trv_hour = 13;
  trv_mins = 15;   
  TravelSlots[di] = trv_hour*4 + trv_mins%15;  
  Destination[di++] = "New Delhi";//11


  for(int iDay = 0; iDay<NumDays; iDay++){
    dep_hour = 0;
    dep_mins = 30;
    DeparSlots[di][iDay*24*4 + dep_hour*4 + dep_mins/15] = true;
  }
  trv_hour = 11;
  trv_mins = 15;  
  TravelSlots[di] = trv_hour*4 + trv_mins%15;  
  Destination[di++] = "Addis Ababa";//12

  for(int iDay = 0; iDay<NumDays; iDay++){
    dep_hour = 13;
    dep_mins = 00;
    DeparSlots[di][iDay*24*4 + dep_hour*4 + dep_mins/15] = true;
  }
  trv_hour = 7;
  trv_mins = 15;  
  TravelSlots[di] = trv_hour*4 + trv_mins%15;    
  Destination[di++] = "Madrid";//13


  for(int iDay = 0; iDay<NumDays; iDay++){
    dep_hour = 13;
    dep_mins = 30;
    DeparSlots[di][iDay*24*4 + dep_hour*4 + dep_mins/15] = true;
  }
  trv_hour = 8;
  trv_mins = 15;  
  TravelSlots[di] = trv_hour*4 + trv_mins%15;   
  Destination[di++] = "Jeddah";//14


  for(int iIter=0; iIter<NumIter; iIter++){
    bool NoSolution = false;
    for(int iDest=0; iDest<NumDest; iDest++){
      break;
      for(int iSlot=0; iSlot<NumSlots; iSlot++){
        if(DeparSlots[iDest][iSlot]==true)
          printf("%s departs on slot nr: %i\n", Destination[iDest].Data(), iSlot);
      }
      //if(iDest>2) break;
    }
    
    bool** PlaneBusy = new bool* [NumAircraft];
    for(int iPlane=0; iPlane<NumAircraft; iPlane++){
      PlaneBusy[iPlane] = new bool [NumSlots];
      for(int iSlot=0; iSlot<NumSlots; iSlot++){
        PlaneBusy[iPlane][iSlot] = false;
      }
    }

    for(int iSlot=0; iSlot<NumSlots; iSlot++){
      std::vector<int> available_aircraft;
      for(int iPlane=0; iPlane<NumAircraft; iPlane++){
        if(PlaneBusy[iPlane][iSlot]== false)
        available_aircraft.push_back(AircraftID[iPlane]);
      }
      for(int iDest=0; iDest<NumDest; iDest++){
        if(DeparSlots[iDest][iSlot]==true){
          if(available_aircraft.size() == 0){
            NoSolution = true;
            break;
          }
          int rnd_plane_id = available_aircraft.at(rangen.Integer(available_aircraft.size()));
          
        }
      }
      if(NoSolution) break;
    }
    if(NoSolution) continue;
  }




  
  




}


void pSigmap_reso(TString InputFileName, const int REBIN, bool SmoothSmear){

  std::vector<std::string> list_of_graphs;
  list_of_graphs.push_back("Red_total");
  list_of_graphs.push_back("Green_total");
  list_of_graphs.push_back("Purple_total");
  list_of_graphs.push_back("Orange_total");
  list_of_graphs.push_back("Blue_total");
  list_of_graphs.push_back("Julich_total");
  list_of_graphs.push_back("Coulomb");
  list_of_graphs.push_back("Red_dashed");
  list_of_graphs.push_back("Green_dashed");
  list_of_graphs.push_back("Purple_dashed");
  list_of_graphs.push_back("Red_dashdotted");
  list_of_graphs.push_back("Green_dashdotted");
  list_of_graphs.push_back("Purple_dashdotted");

  TString err_descr[3];
  err_descr[0] = "def";
  err_descr[1] = "up";
  err_descr[2] = "low";


  TFile fInputReso(TString::Format("%s/p_Sigma/Benedict/pSigmaResolutionFiles%s.root",GetCernBoxDimi(),SmoothSmear?"_Smooth":""), "read");
  fInputReso.ls();
  TH2F* inverted_Resolution_Matrix = (TH2F*)fInputReso.Get(TString::Format("Resolution_Matrix%s",SmoothSmear?"_smooth":""));
  TH1D* Mixed_Event_Distribution = (TH1D*)fInputReso.Get(TString::Format("Mixed_Event_Distribution%s",SmoothSmear?"_smooth":""));

  TFile fOutput(TString::Format("%s/p_Sigma/Benedict/smeared_result_%iMeV%s.root",GetCernBoxDimi(), REBIN, SmoothSmear?"_Smooth":""), "recreate");
  TH2F* hReso_pSigmap = new TH2F("hReso_pSigmap","hReso_pSigmap", 
    inverted_Resolution_Matrix->GetYaxis()->GetNbins(), inverted_Resolution_Matrix->GetYaxis()->GetBinLowEdge(1), inverted_Resolution_Matrix->GetYaxis()->GetBinUpEdge(inverted_Resolution_Matrix->GetYaxis()->GetNbins()),
    inverted_Resolution_Matrix->GetXaxis()->GetNbins(), inverted_Resolution_Matrix->GetXaxis()->GetBinLowEdge(1), inverted_Resolution_Matrix->GetXaxis()->GetBinUpEdge(inverted_Resolution_Matrix->GetXaxis()->GetNbins()));

  for(unsigned uMomX=0; uMomX<hReso_pSigmap->GetXaxis()->GetNbins(); uMomX++){
    double kstarX = hReso_pSigmap->GetXaxis()->GetBinCenter(uMomX+1);
    for(unsigned uMomY=0; uMomY<hReso_pSigmap->GetYaxis()->GetNbins(); uMomY++){
      double kstarY = hReso_pSigmap->GetYaxis()->GetBinCenter(uMomY+1);
      //if(fabs(kstarX-kstarY)>40) hReso_pSigmap->SetBinContent(uMomX+1,uMomY+1,0);
      //else 
      hReso_pSigmap->SetBinContent(uMomX+1,uMomY+1,inverted_Resolution_Matrix->GetBinContent(uMomY+1,uMomX+1));
    }
  }
  hReso_pSigmap->Rebin2D(REBIN,REBIN);
  hReso_pSigmap->Write();

  TH1F* hME = new TH1F("hME","hME",Mixed_Event_Distribution->GetNbinsX(), Mixed_Event_Distribution->GetBinLowEdge(1), Mixed_Event_Distribution->GetXaxis()->GetBinUpEdge(Mixed_Event_Distribution->GetNbinsX()));
  for(unsigned uMom=0; uMom<hME->GetNbinsX(); uMom++){
    double kstar = hME->GetBinCenter(uMom+1);
    if(fabs(kstar-Mixed_Event_Distribution->GetBinCenter(uMom+1))>1e-3){
      printf("BIG WTF\n");
    }
    hME->SetBinContent(uMom+1, Mixed_Event_Distribution->GetBinContent(uMom+1));
    hME->SetBinError(uMom+1, Mixed_Event_Distribution->GetBinError(uMom+1));
  }
  hME->Rebin(REBIN);
  hME->Write();

  unsigned NumMomBins = 50;
  double kMin = 0;
  double kMax = 400;

  for(unsigned uLine=0; uLine<list_of_graphs.size(); uLine++){
    //
    for(unsigned uErr=0; uErr<3; uErr++){
      DLM_Histo<double> dlm_Theory;
      dlm_Theory.SetUp(1);
      dlm_Theory.SetUp(0,NumMomBins,kMin,kMax);
      dlm_Theory.Initialize();

      TFile fInputLines(InputFileName, "read");
      TGraphAsymmErrors* gTheory = (TGraphAsymmErrors*)fInputLines.Get(list_of_graphs.at(uLine).c_str());

      TH1F* h_Theory = new TH1F(TString::Format("h_%s_%s", list_of_graphs.at(uLine).c_str(),err_descr[uErr].Data()),TString::Format("h_%s_%s", list_of_graphs.at(uLine).c_str(),err_descr[uErr].Data()),NumMomBins,kMin,kMax);
      TGraph* gTheory_BAND = new TGraph();
      gTheory_BAND->SetName("gTheory_BAND");

      for(unsigned uPt=0; uPt<gTheory->GetN(); uPt++){
        double kstar, ck;
        gTheory->GetPoint(uPt,kstar, ck);
        if(uErr==0) gTheory_BAND->SetPoint(uPt,kstar, ck);
        else if(uErr==1) gTheory_BAND->SetPoint(uPt,kstar, ck+gTheory->GetErrorYhigh(uPt));
        else gTheory_BAND->SetPoint(uPt,kstar, ck-gTheory->GetErrorYlow(uPt));
        //printf("%u %f %f\n", uPt, kstar, ck);
        //printf("%u %f %f | %f %f\n", uPt, kstar, ck, gTheory->GetErrorYhigh(uPt), gTheory->GetErrorYlow(uPt));
        
      }

      for(unsigned uMom=0; uMom<NumMomBins; uMom++){
        double kstar, ck;
        kstar = dlm_Theory.GetBinCenter(0, uMom);
        ck = gTheory_BAND->Eval(kstar);
        dlm_Theory.SetBinContent(uMom,ck);
        h_Theory->SetBinContent(uMom+1,ck);
        //printf("%u %f %f %f\n", uMom, kstar, ck, dlm_Theory.Eval(kstar));
      }
      fOutput.cd();
      h_Theory->Write();

    
      DLM_Ck dlmCk_Theory(dlm_Theory);
      dlmCk_Theory.SetCutOff(400,800);
      dlmCk_Theory.Update(true);
      
      DLM_CkDecomposition dlmDec_pSigma("pSigmap",0,dlmCk_Theory,hReso_pSigmap);
      dlmDec_pSigma.AddPhaseSpace(hME);
      dlmDec_pSigma.Update(true, true);


      
      TH1F* h_Theory_Smeared = new TH1F(TString::Format("hsmeared_%s_%s", list_of_graphs.at(uLine).c_str(),err_descr[uErr].Data()),TString::Format("hsmeared_%s_%s", list_of_graphs.at(uLine).c_str(),err_descr[uErr].Data()),NumMomBins,kMin,kMax);
      TH1F* h_Theory_Smearing_Ratio = new TH1F(TString::Format("hratio_%s_%s", list_of_graphs.at(uLine).c_str(),err_descr[uErr].Data()),TString::Format("hratio_%s_%s", list_of_graphs.at(uLine).c_str(),err_descr[uErr].Data()),NumMomBins,kMin,kMax);
      h_Theory_Smearing_Ratio->GetYaxis()->SetTitle("smeared / original");
      for(unsigned uMom=0; uMom<NumMomBins; uMom++){
        double kstar = h_Theory_Smeared->GetBinCenter(uMom+1);
        h_Theory_Smeared->SetBinContent(uMom+1,dlmDec_pSigma.EvalCk(kstar));
        h_Theory_Smearing_Ratio->SetBinContent(uMom+1,dlmDec_pSigma.EvalCk(kstar)/dlmCk_Theory.Eval(kstar));
        //printf("%u %f %f %f %f\n", uMom, kstar, dlm_Theory.Eval(kstar), dlmCk_Theory.Eval(kstar), dlmDec_pSigma.EvalCk(kstar));
      }
      fOutput.cd();
      h_Theory_Smeared->Write();
      h_Theory_Smearing_Ratio->Write();

      delete h_Theory_Smeared;
      delete h_Theory_Smearing_Ratio;
      delete gTheory_BAND;

    }


  }

}

void pSigmap_reso_SYST(TString WorkFolder, TString BaseFileName, bool SmoothSmear, std::string specific_model = ""){

  const int REBIN = 1;

  TRandom3 rangen(23);

  TString InputFileName = WorkFolder+BaseFileName+".root";
  TString OutputFileName; 

  std::vector<std::string> list_of_theories;
  if(specific_model!=""){
    list_of_theories.push_back(specific_model);
    OutputFileName = WorkFolder+"/Output/"+BaseFileName+TString::Format("_smeared%s_%s_REB%i.root",SmoothSmear?"_Smooth":"",specific_model.c_str(), REBIN);
  }
  else{
    list_of_theories.push_back("Sing_ReidA");
    list_of_theories.push_back("Trip_ReidA");
    list_of_theories.push_back("Sing_ReidB");
    list_of_theories.push_back("Trip_ReidB");
    list_of_theories.push_back("Gauss14");
    list_of_theories.push_back("Gauss18");
    OutputFileName = WorkFolder+"/Output/"+BaseFileName+TString::Format("_smeared%s.root",SmoothSmear?"_Smooth":"");
  }



  //TFile fInputReso(TString::Format("%s/p_Sigma/Benedict/pSigmaResolutionFiles.root",GetCernBoxDimi()), "read");
  //TH2F* inverted_Resolution_Matrix = (TH2F*)fInputReso.Get("Resolution_Matrix");
  //TH1D* Mixed_Event_Distribution = (TH1D*)fInputReso.Get("Mixed_Event_Distribution");

  TFile fInputReso(TString::Format("%s/p_Sigma/Benedict/pSigmaResolutionFiles%s.root",GetCernBoxDimi(),SmoothSmear?"_Smooth":""), "read");
  fInputReso.ls();
  TH2F* inverted_Resolution_Matrix = (TH2F*)fInputReso.Get(TString::Format("Resolution_Matrix%s",SmoothSmear?"_smooth":""));
  TH1D* Mixed_Event_Distribution = (TH1D*)fInputReso.Get(TString::Format("Mixed_Event_Distribution%s",SmoothSmear?"_smooth":""));

  TFile fOutput(OutputFileName, "recreate");
  TH2F* hReso_pSigmap = new TH2F("hReso_pSigmap","hReso_pSigmap", 
    inverted_Resolution_Matrix->GetYaxis()->GetNbins(), inverted_Resolution_Matrix->GetYaxis()->GetBinLowEdge(1), inverted_Resolution_Matrix->GetYaxis()->GetBinUpEdge(inverted_Resolution_Matrix->GetYaxis()->GetNbins()),
    inverted_Resolution_Matrix->GetXaxis()->GetNbins(), inverted_Resolution_Matrix->GetXaxis()->GetBinLowEdge(1), inverted_Resolution_Matrix->GetXaxis()->GetBinUpEdge(inverted_Resolution_Matrix->GetXaxis()->GetNbins()));

  for(unsigned uMomX=0; uMomX<hReso_pSigmap->GetXaxis()->GetNbins(); uMomX++){
    double kstarX = hReso_pSigmap->GetXaxis()->GetBinCenter(uMomX+1);
    for(unsigned uMomY=0; uMomY<hReso_pSigmap->GetYaxis()->GetNbins(); uMomY++){
      double kstarY = hReso_pSigmap->GetYaxis()->GetBinCenter(uMomY+1);
      hReso_pSigmap->SetBinContent(uMomX+1,uMomY+1,inverted_Resolution_Matrix->GetBinContent(uMomY+1,uMomX+1));
    }
  }
  hReso_pSigmap->Rebin2D(REBIN,REBIN);
  hReso_pSigmap->Write();


  TH1F* hME = new TH1F("hME","hME",Mixed_Event_Distribution->GetNbinsX(), Mixed_Event_Distribution->GetBinLowEdge(1), Mixed_Event_Distribution->GetXaxis()->GetBinUpEdge(Mixed_Event_Distribution->GetNbinsX()));
  for(unsigned uMom=0; uMom<hME->GetNbinsX(); uMom++){
    double kstar = hME->GetBinCenter(uMom+1);
    if(fabs(kstar-Mixed_Event_Distribution->GetBinCenter(uMom+1))>1e-3){
      printf("BIG WTF\n");
    }
    hME->SetBinContent(uMom+1, Mixed_Event_Distribution->GetBinContent(uMom+1));
    hME->SetBinError(uMom+1, Mixed_Event_Distribution->GetBinError(uMom+1));
  }
  hME->Rebin(REBIN);
  hME->Write();

  TFile fBigInput(InputFileName, "read");
  for(unsigned uTh=0; uTh<list_of_theories.size(); uTh++){
    printf("%s\n",list_of_theories.at(uTh).c_str());
    fBigInput.cd();
    TH3D* big_histo = (TH3D*)fBigInput.Get(list_of_theories.at(uTh).c_str());
    fOutput.cd();
    //TH3D* big_histo_smeared = (TH3D*)big_histo->Clone(TString::Format("%s_smeared", list_of_theories.at(uTh).c_str()));
    
      
    unsigned NumMomBins = big_histo->GetXaxis()->GetNbins()/REBIN;
    unsigned LostBins = big_histo->GetXaxis()->GetNbins()%REBIN;
    double kMin = big_histo->GetXaxis()->GetBinLowEdge(1);
    double kMax = big_histo->GetXaxis()->GetBinUpEdge(big_histo->GetXaxis()->GetNbins()-LostBins);
    printf("kMin/kMax = %.f %.f\n",kMin,kMax);
    DLM_Histo<double> dlm_Theory;
    dlm_Theory.SetUp(1);
    dlm_Theory.SetUp(0,NumMomBins,kMin,kMax);
    dlm_Theory.Initialize();


    TH3D* big_histo_smeared = new TH3D(TString::Format("%s_smeared", list_of_theories.at(uTh).c_str()),TString::Format("%s_smeared", list_of_theories.at(uTh).c_str()),
          NumMomBins,kMin,kMax,
          big_histo->GetYaxis()->GetNbins(),big_histo->GetYaxis()->GetBinLowEdge(1),big_histo->GetYaxis()->GetBinUpEdge(big_histo->GetYaxis()->GetNbins()),
          big_histo->GetZaxis()->GetNbins(),big_histo->GetZaxis()->GetBinLowEdge(1),big_histo->GetZaxis()->GetBinUpEdge(big_histo->GetZaxis()->GetNbins())
        );


    unsigned rnd_var = rangen.Integer(big_histo->GetYaxis()->GetNbins());
    unsigned rnd_rad = rangen.Integer(big_histo->GetZaxis()->GetNbins());
    


    TH1F* one_example = new TH1F(TString::Format("example_%s_smeared_%u_%u", list_of_theories.at(uTh).c_str(),rnd_var,rnd_rad),TString::Format("example_%s_smeared_%u_%u", list_of_theories.at(uTh).c_str(),rnd_var,rnd_rad),NumMomBins,kMin,kMax);
    TH1F* one_exampleO = new TH1F(TString::Format("exampleO_%s_smeared_%u_%u", list_of_theories.at(uTh).c_str(),rnd_var,rnd_rad),TString::Format("exampleO_%s_smeared_%u_%u", list_of_theories.at(uTh).c_str(),rnd_var,rnd_rad),NumMomBins,kMin,kMax);
      

    for(unsigned uVar=0; uVar<big_histo->GetYaxis()->GetNbins(); uVar++){
      //TH1F* h_Theory = new TH1F(TString::Format("h_%s_v%u", list_of_theories.at(uTh).c_str(),TString::Format("h_%s_v%u", list_of_theories.at(uTh).c_str(),uVar)),NumMomBins,kMin,kMax);
      //printf("  var %u\n", uVar);
      for(unsigned uRad=0; uRad<big_histo->GetZaxis()->GetNbins(); uRad++){
        for(unsigned uMom=0; uMom<NumMomBins; uMom++){
          double kstar, ck;
          kstar = dlm_Theory.GetBinCenter(0, uMom);
          unsigned origBin = big_histo->GetXaxis()->FindBin(kstar);
          ck = big_histo->GetBinContent(origBin, uVar+1, uRad+1);
          dlm_Theory.SetBinContent(uMom,ck);
          if(uVar==rnd_var && uRad==rnd_rad){
            one_exampleO->SetBinContent(uMom+1, ck);
          }
          //h_Theory->SetBinContent(uMom+1,ck);
        }
        fOutput.cd();
      
        DLM_Ck dlmCk_Theory(dlm_Theory);
        dlmCk_Theory.SetCutOff(550,800);
        dlmCk_Theory.Update(true);
        
        DLM_CkDecomposition dlmDec_pSigma("pSigmap",0,dlmCk_Theory,hReso_pSigmap);
        dlmDec_pSigma.AddPhaseSpace(hME);
        dlmDec_pSigma.Update(true, true);

        //TH1F* h_Theory_Smeared = new TH1F(TString::Format("hsmeared_%s_%s", list_of_graphs.at(uLine).c_str(),err_descr[uErr].Data()),TString::Format("hsmeared_%s_%s", list_of_graphs.at(uLine).c_str(),err_descr[uErr].Data()),NumMomBins,kMin,kMax);
        //TH1F* h_Theory_Smearing_Ratio = new TH1F(TString::Format("hratio_%s_%s", list_of_graphs.at(uLine).c_str(),err_descr[uErr].Data()),TString::Format("hratio_%s_%s", list_of_graphs.at(uLine).c_str(),err_descr[uErr].Data()),NumMomBins,kMin,kMax);
        //h_Theory_Smearing_Ratio->GetYaxis()->SetTitle("smeared / original");
        for(unsigned uMom=0; uMom<NumMomBins; uMom++){
          double kstar = dlm_Theory.GetBinCenter(0, uMom);
          big_histo_smeared->SetBinContent(uMom+1, uVar+1, uRad+1, dlmDec_pSigma.EvalCk(kstar));
          if(uVar==rnd_var && uRad==rnd_rad){
            one_example->SetBinContent(uMom+1, dlmDec_pSigma.EvalCk(kstar));
          }
          //double kstar = h_Theory_Smeared->GetBinCenter(uMom+1);
          //h_Theory_Smeared->SetBinContent(uMom+1,dlmDec_pSigma.EvalCk(kstar));
          //h_Theory_Smearing_Ratio->SetBinContent(uMom+1,dlmDec_pSigma.EvalCk(kstar)/dlmCk_Theory.Eval(kstar));
          //printf("%u %f %f %f %f\n", uMom, kstar, dlm_Theory.Eval(kstar), dlmCk_Theory.Eval(kstar), dlmDec_pSigma.EvalCk(kstar));
        }
        //fOutput.cd();
      }//uRad
    }//uVar

    TH1F* qa_example = new TH1F(TString::Format("exampleQA_%s_smeared_%u_%u", list_of_theories.at(uTh).c_str(),rnd_var,rnd_rad),TString::Format("exampleQA_%s_smeared_%u_%u", list_of_theories.at(uTh).c_str(),rnd_var,rnd_rad),NumMomBins,kMin,kMax);
    for(unsigned uMom=0; uMom<NumMomBins; uMom++){
      qa_example->SetBinContent(uMom+1, big_histo_smeared->GetBinContent(uMom+1,rnd_var+1,rnd_rad+1));
    }

    fOutput.cd();
    big_histo_smeared->Write();
    one_example->Write();
    one_exampleO->Write();
    qa_example->Write();
    delete big_histo_smeared;
    delete one_example;
    delete one_exampleO;
  }
}




/*
void pSigmap_reso_SYST_TGraph(TString WorkFolder, TString BaseFileName, bool SmoothSmear, std::string specific_model = ""){

  const int REBIN = 8;

  TRandom3 rangen(23);

  TString InputFileName = WorkFolder+BaseFileName+".root";
  TString OutputFileName; 

  std::vector<std::string> list_of_theories;
  if(specific_model!=""){
    list_of_theories.push_back(specific_model);
    OutputFileName = WorkFolder+"/Output/"+BaseFileName+TString::Format("_smeared%s_%s_REB%i.root",SmoothSmear?"_Smooth":"",specific_model.c_str(), REBIN);
  }
  else{
    list_of_theories.push_back("Red_total");
    list_of_theories.push_back("Green_total");
    list_of_theories.push_back("Purple_total");
    list_of_theories.push_back("Orange_total");
    list_of_theories.push_back("Blue_total");
    list_of_theories.push_back("Julich_total");

    list_of_theories.push_back("Coulomb");
    list_of_theories.push_back("Red_dashed");
    list_of_theories.push_back("Green_dashed");
    list_of_theories.push_back("Purple_dashed");
    list_of_theories.push_back("Red_dashdotted");
    list_of_theories.push_back("Green_dashdotted");
    list_of_theories.push_back("Purple_dashdotted");
    OutputFileName = WorkFolder+"/Output/"+BaseFileName+TString::Format("_smeared%s.root",SmoothSmear?"_Smooth":"");
  }



  //TFile fInputReso(TString::Format("%s/p_Sigma/Benedict/pSigmaResolutionFiles.root",GetCernBoxDimi()), "read");
  //TH2F* inverted_Resolution_Matrix = (TH2F*)fInputReso.Get("Resolution_Matrix");
  //TH1D* Mixed_Event_Distribution = (TH1D*)fInputReso.Get("Mixed_Event_Distribution");

  TFile fInputReso(TString::Format("%s/p_Sigma/Benedict/pSigmaResolutionFiles%s.root",GetCernBoxDimi(),SmoothSmear?"_Smooth":""), "read");
  fInputReso.ls();
  TH2F* inverted_Resolution_Matrix = (TH2F*)fInputReso.Get(TString::Format("Resolution_Matrix%s",SmoothSmear?"_smooth":""));
  TH1D* Mixed_Event_Distribution = (TH1D*)fInputReso.Get(TString::Format("Mixed_Event_Distribution%s",SmoothSmear?"_smooth":""));

  TFile fOutput(OutputFileName, "recreate");
  TH2F* hReso_pSigmap = new TH2F("hReso_pSigmap","hReso_pSigmap", 
    inverted_Resolution_Matrix->GetYaxis()->GetNbins(), inverted_Resolution_Matrix->GetYaxis()->GetBinLowEdge(1), inverted_Resolution_Matrix->GetYaxis()->GetBinUpEdge(inverted_Resolution_Matrix->GetYaxis()->GetNbins()),
    inverted_Resolution_Matrix->GetXaxis()->GetNbins(), inverted_Resolution_Matrix->GetXaxis()->GetBinLowEdge(1), inverted_Resolution_Matrix->GetXaxis()->GetBinUpEdge(inverted_Resolution_Matrix->GetXaxis()->GetNbins()));

  for(unsigned uMomX=0; uMomX<hReso_pSigmap->GetXaxis()->GetNbins(); uMomX++){
    double kstarX = hReso_pSigmap->GetXaxis()->GetBinCenter(uMomX+1);
    for(unsigned uMomY=0; uMomY<hReso_pSigmap->GetYaxis()->GetNbins(); uMomY++){
      double kstarY = hReso_pSigmap->GetYaxis()->GetBinCenter(uMomY+1);
      hReso_pSigmap->SetBinContent(uMomX+1,uMomY+1,inverted_Resolution_Matrix->GetBinContent(uMomY+1,uMomX+1));
    }
  }
  hReso_pSigmap->Rebin2D(REBIN,REBIN);
  hReso_pSigmap->Write();


  TH1F* hME = new TH1F("hME","hME",Mixed_Event_Distribution->GetNbinsX(), Mixed_Event_Distribution->GetBinLowEdge(1), Mixed_Event_Distribution->GetXaxis()->GetBinUpEdge(Mixed_Event_Distribution->GetNbinsX()));
  for(unsigned uMom=0; uMom<hME->GetNbinsX(); uMom++){
    double kstar = hME->GetBinCenter(uMom+1);
    if(fabs(kstar-Mixed_Event_Distribution->GetBinCenter(uMom+1))>1e-3){
      printf("BIG WTF\n");
    }
    hME->SetBinContent(uMom+1, Mixed_Event_Distribution->GetBinContent(uMom+1));
    hME->SetBinError(uMom+1, Mixed_Event_Distribution->GetBinError(uMom+1));
  }
  hME->Rebin(REBIN);
  hME->Write();






  TFile fBigInput(InputFileName, "read");
  for(unsigned uTh=0; uTh<list_of_theories.size(); uTh++){
    printf("%s\n",list_of_theories.at(uTh).c_str());
    fBigInput.cd();
    TGraphAsymmErrors* data_histo = (TGraphAsymmErrors*)fBigInput.Get(list_of_theories.at(uTh).c_str());
    fOutput.cd();
    //TH3D* big_histo_smeared = (TH3D*)big_histo->Clone(TString::Format("%s_smeared", list_of_theories.at(uTh).c_str()));
    

    double dummy_x, dummy_y;
    double kMin = 0;
    data_histo->GetPoint(data_histo->GetN(),dummy_x, dummy_y);
    double kMax = dummy_x;
    unsigned NumMomBins = int(kMax)/REBIN + 1;
    kMax = kMin + NumMomBins*REBIN;
    //unsigned LostBins = big_histo->GetXaxis()->GetNbins()%REBIN;
    
    
    printf("kMin/kMax nb%u = %.f %.f\n",kMin,kMax,NumMomBins);
    DLM_Histo<double> dlm_Theory;
    dlm_Theory.SetUp(1);
    dlm_Theory.SetUp(0,NumMomBins,kMin,kMax);
    dlm_Theory.Initialize();


    TH3D* big_histo_smeared = new TH3D(TString::Format("%s_smeared", list_of_theories.at(uTh).c_str()),TString::Format("%s_smeared", list_of_theories.at(uTh).c_str()),
          NumMomBins,kMin,kMax,
          big_histo->GetYaxis()->GetNbins(),big_histo->GetYaxis()->GetBinLowEdge(1),big_histo->GetYaxis()->GetBinUpEdge(big_histo->GetYaxis()->GetNbins()),
          big_histo->GetZaxis()->GetNbins(),big_histo->GetZaxis()->GetBinLowEdge(1),big_histo->GetZaxis()->GetBinUpEdge(big_histo->GetZaxis()->GetNbins())
        );


    unsigned rnd_var = rangen.Integer(big_histo->GetYaxis()->GetNbins());
    unsigned rnd_rad = rangen.Integer(big_histo->GetZaxis()->GetNbins());
    


    TH1F* one_example = new TH1F(TString::Format("example_%s_smeared_%u_%u", list_of_theories.at(uTh).c_str(),rnd_var,rnd_rad),TString::Format("example_%s_smeared_%u_%u", list_of_theories.at(uTh).c_str(),rnd_var,rnd_rad),NumMomBins,kMin,kMax);
    TH1F* one_exampleO = new TH1F(TString::Format("exampleO_%s_smeared_%u_%u", list_of_theories.at(uTh).c_str(),rnd_var,rnd_rad),TString::Format("exampleO_%s_smeared_%u_%u", list_of_theories.at(uTh).c_str(),rnd_var,rnd_rad),NumMomBins,kMin,kMax);
      

    for(unsigned uVar=0; uVar<big_histo->GetYaxis()->GetNbins(); uVar++){
      //TH1F* h_Theory = new TH1F(TString::Format("h_%s_v%u", list_of_theories.at(uTh).c_str(),TString::Format("h_%s_v%u", list_of_theories.at(uTh).c_str(),uVar)),NumMomBins,kMin,kMax);
      //printf("  var %u\n", uVar);
      for(unsigned uRad=0; uRad<big_histo->GetZaxis()->GetNbins(); uRad++){
        for(unsigned uMom=0; uMom<NumMomBins; uMom++){
          double kstar, ck;
          kstar = dlm_Theory.GetBinCenter(0, uMom);
          unsigned origBin = big_histo->GetXaxis()->FindBin(kstar);
          ck = big_histo->GetBinContent(origBin, uVar+1, uRad+1);
          dlm_Theory.SetBinContent(uMom,ck);
          if(uVar==rnd_var && uRad==rnd_rad){
            one_exampleO->SetBinContent(uMom+1, ck);
          }
          //h_Theory->SetBinContent(uMom+1,ck);
        }
        fOutput.cd();
      
        DLM_Ck dlmCk_Theory(dlm_Theory);
        dlmCk_Theory.SetCutOff(550,800);
        dlmCk_Theory.Update(true);
        
        DLM_CkDecomposition dlmDec_pSigma("pSigmap",0,dlmCk_Theory,hReso_pSigmap);
        dlmDec_pSigma.AddPhaseSpace(hME);
        dlmDec_pSigma.Update(true, true);

        //TH1F* h_Theory_Smeared = new TH1F(TString::Format("hsmeared_%s_%s", list_of_graphs.at(uLine).c_str(),err_descr[uErr].Data()),TString::Format("hsmeared_%s_%s", list_of_graphs.at(uLine).c_str(),err_descr[uErr].Data()),NumMomBins,kMin,kMax);
        //TH1F* h_Theory_Smearing_Ratio = new TH1F(TString::Format("hratio_%s_%s", list_of_graphs.at(uLine).c_str(),err_descr[uErr].Data()),TString::Format("hratio_%s_%s", list_of_graphs.at(uLine).c_str(),err_descr[uErr].Data()),NumMomBins,kMin,kMax);
        //h_Theory_Smearing_Ratio->GetYaxis()->SetTitle("smeared / original");
        for(unsigned uMom=0; uMom<NumMomBins; uMom++){
          double kstar = dlm_Theory.GetBinCenter(0, uMom);
          big_histo_smeared->SetBinContent(uMom+1, uVar+1, uRad+1, dlmDec_pSigma.EvalCk(kstar));
          if(uVar==rnd_var && uRad==rnd_rad){
            one_example->SetBinContent(uMom+1, dlmDec_pSigma.EvalCk(kstar));
          }
          //double kstar = h_Theory_Smeared->GetBinCenter(uMom+1);
          //h_Theory_Smeared->SetBinContent(uMom+1,dlmDec_pSigma.EvalCk(kstar));
          //h_Theory_Smearing_Ratio->SetBinContent(uMom+1,dlmDec_pSigma.EvalCk(kstar)/dlmCk_Theory.Eval(kstar));
          //printf("%u %f %f %f %f\n", uMom, kstar, dlm_Theory.Eval(kstar), dlmCk_Theory.Eval(kstar), dlmDec_pSigma.EvalCk(kstar));
        }
        //fOutput.cd();
      }//uRad
    }//uVar

    TH1F* qa_example = new TH1F(TString::Format("exampleQA_%s_smeared_%u_%u", list_of_theories.at(uTh).c_str(),rnd_var,rnd_rad),TString::Format("exampleQA_%s_smeared_%u_%u", list_of_theories.at(uTh).c_str(),rnd_var,rnd_rad),NumMomBins,kMin,kMax);
    for(unsigned uMom=0; uMom<NumMomBins; uMom++){
      qa_example->SetBinContent(uMom+1, big_histo_smeared->GetBinContent(uMom+1,rnd_var+1,rnd_rad+1));
    }

    fOutput.cd();
    big_histo_smeared->Write();
    one_example->Write();
    one_exampleO->Write();
    qa_example->Write();
    delete big_histo_smeared;
    delete one_example;
    delete one_exampleO;
  }
}
*/

DLM_Histo<float>* dlmSrc = NULL;
double IdioticSource(double* pars){
  return dlmSrc->Eval(pars[1]);
}

void Quick_Max_Carla_check(){
  TString InputFileName = TString::Format("%s/Max/CECA_QA/Seed1234.0.root",GetCernBoxDimi());
  TFile fInput(InputFileName, "read");
  //fInput.ls();
  const unsigned NumMt = 7;
  
  const unsigned kSteps = 45;
  const double kMin = 0;
  const double kMax = 180;

  std::vector<float> src_eff_g;
  src_eff_g.push_back(1.37);
  src_eff_g.push_back(1.30);
  src_eff_g.push_back(1.27);
  src_eff_g.push_back(1.21);
  src_eff_g.push_back(1.14);
  src_eff_g.push_back(1.05);
  src_eff_g.push_back(0.96);


  TFile fOutput(TString::Format("%s/Max/CECA_QA/xCheck.root",GetCernBoxDimi()), "recreate");
  
  
  for(unsigned uMt=0; uMt<NumMt; uMt++){
    TH1F* hSrc = (TH1F*)fInput.Get(TString::Format("Source_%u",uMt+1));
    printf("hSrc = %p\n",hSrc);
    hSrc->Scale(1./hSrc->Integral(),"width");
    if(dlmSrc) {delete dlmSrc; dlmSrc=NULL;}
    dlmSrc = Convert_TH1F_DlmHisto(hSrc);
    TGraph* gCk_AI = (TGraph*)fInput.Get(TString::Format("CFpp_check_%u",uMt+1));
    CATS Kitty_AV18;
    Kitty_AV18.SetMomBins(kSteps, kMin, kMax);
    Kitty_AV18.SetNotifications(CATS::nWarning);
    DLM_CommonAnaFunctions AnalysisObject; 
    AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());

    AnalysisObject.SetUpCats_pp(Kitty_AV18,"AV18","",0,0);
    CATSparameters cPars (CATSparameters::tSource, 0, true);
    Kitty_AV18.SetAnaSource(IdioticSource, cPars);
    Kitty_AV18.SetUseAnalyticSource(true);
    Kitty_AV18.KillTheCat(CATS::kAllChanged);

    CATS Kitty_EffG;
    Kitty_EffG.SetMomBins(kSteps, kMin, kMax);
    Kitty_EffG.SetNotifications(CATS::nWarning);

    AnalysisObject.SetUpCats_pp(Kitty_EffG,"AV18","Gauss",0,0);
    Kitty_EffG.SetAnaSource(0, src_eff_g.at(uMt));
    Kitty_EffG.KillTheCat(CATS::kAllChanged);

    fOutput.cd();
    TH1F* hCkXcheck = new TH1F(TString::Format("hCkXcheck_%u",uMt+1),TString::Format("hCkXcheck_%u",uMt+1),kSteps,kMin,kMax);
    TH1F* hCkEffG = new TH1F(TString::Format("hCkEffG_%u",uMt+1),TString::Format("hCkEffG_%u",uMt+1),kSteps,kMin,kMax);
    
    for(unsigned uMom=0; uMom<kSteps; uMom++){
      hCkXcheck->SetBinContent(uMom+1, Kitty_AV18.GetCorrFun(uMom));
      hCkEffG->SetBinContent(uMom+1, Kitty_EffG.GetCorrFun(uMom));
    }

    fOutput.cd();
    gCk_AI->Write();
    hCkXcheck->Write();
    hCkEffG->Write();

    delete hCkXcheck;
    delete hCkEffG;
  }


}




//
int OTHERTASKS(int argc, char *argv[]){

  Quick_Max_Carla_check(); return 0;

  //pSigmap_reso_SYST("/home/dimihayl/CernBox/Sync/p_Sigma/Benedict/", "pSigma_FitFunctions", atoi(argv[2]), argv[1]);  return 0;
  //pSigmap_reso("/home/dimihayl/CernBox/Sync/p_Sigma/Benedict/pSigma_Allmodels.root",8,true); return 0;


  //PlaneTimeSlots(); return 0;

  //pL_DanielPolish_test(); return 0;

  //Toy_Potential_Compare(); return 0;
  //pp_asymptotic_check_Epelbaum_v1(); return 0;
  //Sofia_test1(); return 0;
  //Toy_Potential_Scan(atoi(argv[2]), atoi(argv[1])); return 0;
  //Toy_Potential_Scan_FloatGauss(atoi(argv[2]), atoi(argv[1])); return 0;
  //AirSim_Schedule(23,4000000,4,2,0); return 0;
  //effect_of_square_barrier(); return 0;

  ThreeBodyDecay(atoi(argv[1]), atoi(argv[2])); return 0;

  //USR_Potential_ShiftedSource(); return 0;
  //big_core(); return 0;

  //test_screen_coulomb(); return 0;

  //get_eta_fractions(); return 0;
  //pp_asymptotic(); return 0;
  //Ghetto_pPhi_check_20250426();

  //L_pi_CommonSource_based_on_p_pic(); return 0;
  //extrap_me_for_3b_v1(); return 0;
  //pL_LL_Decay(atoi(argv[1]), atoi(argv[2])); return 0;

  //pp_LL_Decay(atoi(argv[1]), atoi(argv[2])); return 0;



    //pp_large_RSM_source();
    
    //DummyBootTest_pSigma_inspired();
    //xCheck_Benedict_pSigma_v1();

    //pd_withQS();

    //return 0;

  //TestLambdaKstar_Baseline("Baseline");
  //TestLambdaKstar_Baseline("NewCode");
  //TestLambdaKstar_KstarDep("Static");
  //TestLambdaKstar_KstarDep("Dynamic");

//Quick_pXi(); return 0;
    //pp_reff(1.08); //return 0;
    //pp_reff(1.11); //return 0;
    //pp_reff(1.15); return 0;


    //pSp_test2(); return 0;



  
  //Raffa_pLambda_Gauss();
  //return 0;

  //StableDisto_DlmRan_Test(); return 0;
  //pp_Gauss_Cauchy_Mix(); return 0;

  //Unfold_test1();
  //generate_pp_Neelima();

  //pp_pSp_Decay(atoi(argv[1]), atoi(argv[2]));
  //piXi1530_piXi_Decay(atoi(argv[1]), atoi(argv[2]));
  //KplusXi1530_KplusXi_Decay(atoi(argv[1]), atoi(argv[2]));

  //Test_chiral_2023("NLO19_600_291_sd","NLO19_600_141_sd"); 
  //Test_chiral_2023("NLO19_600_253_sd","NLO19_600_141_sd"); 
  //Test_chiral_2023("NLO19_600_291_sd","NLO19_600_132_sd");
  //Test_chiral_2023("NLO19_600_253_sd","NLO19_600_132_sd");
  //Test_chiral_2023("NLO19_600_291_sd","NLO19_600_115_sd");
  //Test_chiral_2023("NLO19_600_253_sd","NLO19_600_115_sd");

  //Test_chiral_2023("NLO19_600_291_s_elastic","NLO19_600_141_s_elastic");
  //Test_chiral_2023("NLO19_600_291_s","NLO19_600_141_s");
  //Test_chiral_2023("NLO19_600_291_sp","NLO19_600_141_sp");
  //Test_chiral_2023("NLO19_600_291_sd","NLO19_600_141_sd");
  //Test_chiral_2023("NLO19_600_291_spd","NLO19_600_141_spd");
  
  return 0;

  //PotentialDesignerEngine(argv[1]); return 0;
  //Test_CoulombEffRangeExp(); return 0;

  //pSigmaPlus_pp(argv[1], atof(argv[2]));
  //pSigmaPlus_NLO19_data();

  //pn_Xchecks_1();

  //
  //TestNewCutOff();

  //DongFang_Example1(); return 0;
  //coal_test_1(); return 0;
  //test_sqwell(); return 0;

  //TestMergeSort();return 0;
  //pp_at_different_radii("Gauss");return 0;
//PlotAv18();return 0;

  //ShowEffectOfCum(); return 0;
  //Jaime_test1();
  //Test_div_any_th1f();
  //Test_align_th1f();
  //TestDLM_Sample();
  //AV18peak();
//return 0;
  //DDstar_TetsuoTest1();
  //return 0;
  //Bhawani_pot("Toy1");
  //Bhawani_pot("ND46");
  //Bhawani_pot("NF48");
    //rootmathboost_test1();
    //Test_unfold_yield();
    //Rafa_2body_expCk_CATS();
    //return 0;

    //pp_CompareToNorfolk();
    //pp_pL_CorrectedMC_EXP();
    //ALL_CorrectedMC_EXP();
    //ReweightME();
    //ParametrizeTemplates("p","p");

    //Fit_pL_MickeyMouse();
    //Fit_pp_MickeyMouse();
    //Compare_RotPhi_ME();
    //Fit_pL_MC_SimpleTemplate();
    //Fit_pp_MC_SimpleTemplate();

    //

    //Fit_pp_CommonAncestorTemplate_dPhi(0);
    //Fit_pp_CommonAncestorTemplate_dPhi(1);
    //Fit_pL_CommonAncestorTemplate_dPhi(0);
    //Fit_pL_CommonAncestorTemplate_dPhi(1);

    //Fit_pp_CommonAncestorTemplate_Ck();
    //Fit_pL_CommonAncestorTemplate_Ck();

    //Fit_pL_CommonAncestorTemplate_CrazyGauss();

    //ReadDeuteronWF("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/Coalesce/WaveFunctions/Machleidt/deuwaves_n4lo500_rspace.d",huWF,hwWF);

    //pp_DariuszSource();
    //Test_flat_pol();
    //TestOmegaOmega();
    //pXi_Bug_vs_True_Potential(0,"I0S0");
    //pXi_Bug_vs_True_Potential(1,"I0S1");
    //pXi_Bug_vs_True_Potential(2,"I1S0");
    //pXi_Bug_vs_True_Potential(3,"I1S1");
    //MemBugHunt();
    //DlmHistoMemBugHunt();
    //SmearNonFlatCorrelation();
    //Georgios_LXi_ResoTest(0);
    //Georgios_LXi_ResoTest(1);
    //Georgios_LXi_ResoTest(2);
    //Ledni_SmallRad("ND46");
    //Ledni_SmallRad("NF48");
    //Ledni_SmallRad("NSC97b");
    //Ledni_SmallRad("GaussEmma");
    //Ledni_SmallRad("custom");
    //Ledni_SmallRad("Toy1");
    //Ledni_SmallRad("Yukawa1");
    //Ledni_SmallRad("pKplusI0");
    //Ledni_SmallRad("pKplusI1");
    //Ledni_SmallRad("pKplusYuki");
    //pXi_Pot();
    //Raffa_Errors();
    //ppp_errors(atoi(argv[1]),atoi(argv[2]));
    //ppL_errors(0,0);
    //ppL_smear_test();

    //Silly2body(atoi(argv[1]),atoi(argv[2]));

    //KilledByFeedDown();
    //SmearEffect_Resolution();
    //SmearEffect_pp();
    //Test_new_fold();
    //TestUnfoldSEME();
    //TestUnfoldSEME_pL();
    //TestUnfoldSEME_pXi();
    //Rafa_2body_expCk();
//return 0;
/*
ppp_errors(2,1);
(0) nsig (pval) = 1.10 (0.27080)
(1) nsig (pval) = 5.19 (0.00000)
(2) nsig (pval) = 1.81 (0.06980)
(3) nsig (pval) = 0.29 (0.76880)
(4) nsig (pval) = 0.50 (0.61760)
(5) nsig (pval) = 0.31 (0.75580)
(6) nsig (pval) = 0.38 (0.70180)
(7) nsig (pval) = 1.58 (0.11440)
(8) nsig (pval) = 0.51 (0.60880)
nsig 3 bins = 4.97
nsig 6 bins = 4.30

ppp_errors(2,0);
(0) nsig (pval) = 0.98 (0.32780)
(1) nsig (pval) = 4.57 (0.00000)
(2) nsig (pval) = 1.65 (0.09860)
(3) nsig (pval) = 0.75 (0.45360)
(4) nsig (pval) = 1.16 (0.24600)
(5) nsig (pval) = 0.07 (0.94340)
(6) nsig (pval) = 0.53 (0.59880)
(7) nsig (pval) = 1.65 (0.09980)
(8) nsig (pval) = 0.65 (0.51260)
nsig 3 bins = 4.28
nsig 6 bins = 3.75
*/

    //pp_for_rock(1.2-0.12);
    //pp_for_rock(1.2);
    //pp_for_rock(1.2+0.12);

    //Ledni_SmallRad_Random(atoi(argv[1]),atoi(argv[2]));
    //MakePotentials(atoi(argv[1]));return 0;
  //RoughPiPiPotScan(atoi(argv[1]),atoi(argv[2]));
  //pi_proton();
    //SelectEmmaPotential();
    //StableDisto_Test();
    //Andi_pDminus_1();
    //Fast_Bootstrap_Example();
    //MakeSmoothAngularSourceDisto(NULL);
    //MaxRamona_piSig_ResoTest();
    //Ramona_pK_ResoTest();

    //EmmaDaniel_piD(0.82);
    //EmmaDaniel_piD(0.89);
    //EmmaDaniel_piD(0.75);
    //EmmaDaniel_piD(1.00);

    //EmmaDaniel_piD(0,true);
    //printf("\n----------------------\n");
    //EmmaDaniel_KD(0.81);
    //EmmaDaniel_KD(0.89);
    //EmmaDaniel_KD(0.74);
    //EmmaDaniel_KD(1.00);

//for Oton
    //EmmaDaniel_KD(1.04);
    //EmmaDaniel_KD(1.10);
    //EmmaDaniel_KD(1.13);
    //EmmaDaniel_KD(1.04);
    //EmmaDaniel_KD(1.04);
    //EmmaDaniel_KD(1.04);
    //pp_QS_Tests();
    //pn_potential();
    //pn_Ck(0.8);
    //pn_Ck(1.0);
    //pn_Ck(1.2);
    //pn_Ck(1.5);
    //pn_Ck(2.0);

    //LamLamStudy1();
    //pXi_BUG_TEST();

    //PlugInWaveFunction();
    //ppSource_bugHunting(true);
    //ppSource_bugHunting(false);

    //pd_delayed(0,false);
    //pd_delayed(100,true);

    //pp_C_vs_SI();
    //p_pi_CommonTest1();




    //EmmaDaniel_KD(0,true);
    //const double RAD = 5.11;
    //printf("u(%.2f) = %.4f\n",RAD,Evaluate_d_u(RAD));
    //printf(" w(%.2f) = %.4f\n",RAD,Evaluate_d_w(RAD));

    return 0;
}

