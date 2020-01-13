
#include "OtherTasks.h"

#include "CATS.h"
#include "CommonAnaFunctions.h"
#include "DLM_Potentials.h"
#include "DLM_Random.h"
#include "DLM_CkDecomposition.h"
#include "DLM_Source.h"
#include "DLM_Potentials.h"
#include "DLM_WfModel.h"
#include "DLM_CkModels.h"

#include "TGraph.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TNtuple.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TNtuple.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TFractionFitter.h"
#include "TGenPhaseSpace.h"
#include "TString.h"

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

    DLM_CommonAnaFunctions AnalysisObject;

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
    RescaleFactor[1] = 1;//pL
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

        fitTemplate[uSpec] = new TF1("fitTemplate_"+SpeciesName[uSpec],dimi_fraction_fitter,350,3500,5);

        fitTemplate[uSpec]->SetParameter(0,0.7);
        fitTemplate[uSpec]->SetParLimits(0,0.05,1.4);

        fitTemplate[uSpec]->SetParameter(1,0.15);
        fitTemplate[uSpec]->SetParLimits(1,0.01,1.00);

        fitTemplate[uSpec]->SetParameter(2,0.0);
        fitTemplate[uSpec]->SetParLimits(2,-0.15,0.15);

        fitTemplate[uSpec]->SetParameter(3,0.0);
        fitTemplate[uSpec]->SetParLimits(3,-1e-4*10,1e-4*10);

        fitTemplate[uSpec]->SetParameter(4,0.0);
        fitTemplate[uSpec]->SetParLimits(4,-1e-7*10,1e-7*10);

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
        fitMickeyMouse[uSpec] = new TF1("fitMickeyMouse_"+SpeciesName[uSpec],dimi_fractionPolQS_fitter,0,2500,10);

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
    MM_CatPL->SetChannelWeight(2,0.25*par[11]);
    MM_CatPL->SetChannelWeight(3,0.75*par[11]);
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
    DLM_CommonAnaFunctions AnalysisObject;
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
    TString OutFileName = "Out.root";
    TFile* OutputFile = new TFile(OutputFolder+OutFileName,"recreate");
    hData_pL->Write();

    CATS AB_pL;
    DLM_Ck* Ck_pL;
    AB_pL.SetMomBins(NumMomBins_pL,MomBins_pL);
    AnalysisObject.SetUpCats_pL(AB_pL,"NLO_Coupled_S",SourceDescription,0,202);//NLO_Coupled_S
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
    fit_pL->SetParLimits(4,20.,400.);

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

    TF1* fit_pL_NonFemto = new TF1("fit_pL_NonFemto","(1.-[0]*exp(-pow(x*[1]/197.327,[2])))*(1.+[3]*TMath::Gaus(x,[4],[5],0))*[6]*(1.+[7]*x+[8]*x*x+[9]*x*x*x)",0,4500);
    fit_pL_NonFemto->SetParameter(0,fit_pL->GetParameter(0));
    fit_pL_NonFemto->SetParameter(1,fit_pL->GetParameter(1));
    fit_pL_NonFemto->SetParameter(2,fit_pL->GetParameter(2));
    fit_pL_NonFemto->SetParameter(3,fit_pL->GetParameter(3));
    fit_pL_NonFemto->SetParameter(4,fit_pL->GetParameter(4));
    fit_pL_NonFemto->SetParameter(5,fit_pL->GetParameter(5));
    fit_pL_NonFemto->SetParameter(6,1);
    fit_pL_NonFemto->SetParameter(7,fit_pL->GetParameter(7));
    fit_pL_NonFemto->SetParameter(8,fit_pL->GetParameter(8));
    fit_pL_NonFemto->SetParameter(9,fit_pL->GetParameter(9));

    TGraph fit_nsigma;
    fit_nsigma.SetName("fit_nsigma");

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
    fit_pL_NonFemto->Write();
    fit_nsigma.Write();
    hfit_Ratio->Write();

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
    //TString DataSample = "pp13TeV_HM_Dec19";
    TString DataSample = "pp13TeV_HM_RotPhiDec19";
    DLM_CommonAnaFunctions AnalysisObject;
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
    //TString SourceDescription = "Gauss";
    TString SourceDescription = "McGauss_ResoTM";
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

    double FitMin = 16;
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
fit_pp->FixParameter(8,0.0);
fit_pp->FixParameter(9,0.0);
fit_pp->FixParameter(10,1000);
fit_pp->FixParameter(11,Ck_pp->GetSourcePar(0));

//fit_pp->FixParameter(12,2.0);

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
    DLM_CommonAnaFunctions AnalysisObject;

    //pp
    //TH1F* hData_pp_ME = AnalysisObject.GetAliceExpCorrFun("pp13TeV_HM_Dec19","pp","_0",0,false,-1);
    //TH1F* hData_pp_RE = AnalysisObject.GetAliceExpCorrFun("pp13TeV_HM_RotPhiDec19","pp","_0",0,false,-1);
    //pLambda
    TH1F* hData_pp_ME = AnalysisObject.GetAliceExpCorrFun("pp13TeV_HM_Dec19","pLambda","_0",3,false,-1);
    TH1F* hData_pp_RE = AnalysisObject.GetAliceExpCorrFun("pp13TeV_HM_RotPhiDec19","pLambda","_0",0,false,-1);
/*
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
*/
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

    const unsigned NumBins = hData_pp_ME->GetNbinsX()>hData_pp_RE->GetNbinsX()?hData_pp_RE->GetNbinsX():hData_pp_ME->GetNbinsX();
    const double kMin = hData_pp_RE->GetBinLowEdge(1);
    const double kMax = hData_pp_RE->GetXaxis()->GetBinUpEdge(NumBins);

    printf("NumBins=%u; kMin kMax = %f %f\n",NumBins,kMin,kMax);

    dlmData_pp_ME.SetUp(1);
    dlmData_pp_ME.SetUp(0,NumBins,kMin,kMax);
    dlmData_pp_ME.Initialize();

    dlmData_pp_RE.SetUp(1);
    dlmData_pp_RE.SetUp(0,NumBins,kMin,kMax);
    dlmData_pp_RE.Initialize();

    for(unsigned uBin=0; uBin<NumBins; uBin++){
        dlmData_pp_ME.SetBinContent(uBin,hData_pp_ME->GetBinContent(uBin+1));
        dlmData_pp_ME.SetBinError(uBin,hData_pp_ME->GetBinError(uBin+1));
        dlmData_pp_RE.SetBinContent(uBin,hData_pp_RE->GetBinContent(uBin+1)+0.001);
        dlmData_pp_RE.SetBinError(uBin,hData_pp_RE->GetBinError(uBin+1));
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

int OTHERTASKS(int narg, char** ARGS){
    //pp_CompareToNorfolk();
    //pp_pL_CorrectedMC_EXP();
    //ALL_CorrectedMC_EXP();
    //ReweightME();
    //ParametrizeTemplates("p","p");

    //Fit_pL_MickeyMouse();
    //Fit_pp_MickeyMouse();
    Compare_RotPhi_ME();
}
