
#include "OtherTasks.h"

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
#include "DLM_RootWrapper.h"
#include "DLM_HistoAnalysis.h"
#include "EnvVars.h"
#include "DLM_Fitters.h"

#include "TGraph.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH1D.h"
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
#include "TList.h"
#include "TROOT.h"

#include "gsl_sf_dawson.h"

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

    const double CoreSize = 0.8934;
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
    TFile* fOutput = new TFile(TString::Format("%s/OtherTasks/Georgios_LXi_ResoTest/fOutput_%i.root",GetFemtoOutputFolder(),int(SmoothSampling)),"recreate");
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

void Eval_ScattParameters(CATS& Kitty, double& ScatLen, double& EffRan){
  Kitty.KillTheCat();
  TF1 fitSP("fitSP","[2]*x*x*x*x+0.5*[1]/197.327*x*x+197.327/[0]", 10, 60);
  fitSP.SetParameter(0,ScatLen);
  fitSP.SetParameter(1,EffRan);
  fitSP.FixParameter(2,0);
  TH1F* hFit = new TH1F("hFit","hFit",Kitty.GetNumMomBins(),Kitty.GetMomBinLowEdge(0),Kitty.GetMomBinUpEdge(Kitty.GetNumMomBins()-1));
  for(unsigned uMom=0; uMom<Kitty.GetNumMomBins(); uMom++){
    hFit->SetBinContent(uMom+1,Kitty.GetMomentum(uMom)/tan(Kitty.GetPhaseShift(uMom,0,0)));
  }
  hFit->Fit(&fitSP, "Q, S, N, R, M");
  ScatLen = fitSP.GetParameter(0);
  EffRan = fitSP.GetParameter(1);
  delete hFit;
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
  const double Radius = 1.08;

  CATSparameters sPars(CATSparameters::tSource,1,true);
  sPars.SetParameter(0,Radius);
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
  CATSparameters pPars(CATSparameters::tPotential,4,true);
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
  else{
    pPars.SetParameter(0,-78.42*0.39*4.5);//0.39,0.4
    pPars.SetParameter(1,1.0*1.35);
    pPars.SetParameter(2,741.76*4.5);
    pPars.SetParameter(3,0.45*1.4);
    //NF46 as a stariting point
    //pPars.SetParameter(0,-1327.26*1.0);
    //pPars.SetParameter(1,0.6);
    //pPars.SetParameter(2,2561.56);
    //pPars.SetParameter(3,0.45);
    c_f0 = 0.350;
    c_d0 = 16.330;

  }
  Kitty_SE.SetEpsilonConv(5e-9);
  Kitty_SE.SetEpsilonProp(5e-9);
  Kitty_SE.SetShortRangePotential(0,0,DoubleGaussSum,pPars);
  Kitty_SE.KillTheCat();

  Eval_ScattParameters(Kitty_SE,c_f0,c_d0);
  printf("c_f0 = %.2f fm\n",c_f0);
  printf("c_d0 = %.2f fm\n",c_d0);

  const double Re_f0 = PotentialName=="NSC97b"?0.397:PotentialName=="NF48"?1.511:c_f0;//NSC97b
  const double Im_f0 = 0;
  double Val_d0 = PotentialName=="NSC97b"?10.36:PotentialName=="NF48"?2.549:c_d0;
  const TString OutputFolder = TString::Format("%s/OtherTasks/Ledni_SmallRad/",GetFemtoOutputFolder());
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

  unsigned NumPts = 0;
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
    g_XXX.SetPoint(NumPts,MOM,(Kitty_SE.GetCorrFun(NumPts)-1.-F1F2)/v_f2);
    //the correction concidered as a general addition factor to C(k)
    g_XXXa.SetPoint(NumPts,MOM,(Kitty_SE.GetCorrFun(NumPts)-1.-v_f2-F1F2));
    //the correction concidered as a general multiplication factor to C(k)
    g_XXXm.SetPoint(NumPts,MOM,Kitty_SE.GetCorrFun(NumPts)/(1.+v_f2+F1F2));

    g_delta.SetPoint(NumPts,MOM,Kitty_SE.GetPhaseShift(NumPts,0,0));
    g_kcotdelta.SetPoint(NumPts,MOM,MOM/tan(Kitty_SE.GetPhaseShift(NumPts,0,0)));
    g_cotdelta.SetPoint(NumPts,MOM,1./tan(Kitty_SE.GetPhaseShift(NumPts,0,0)));
    g_tgdelta.SetPoint(NumPts,MOM,tan(Kitty_SE.GetPhaseShift(NumPts,0,0)));
    g_sindelta.SetPoint(NumPts,MOM,sin(Kitty_SE.GetPhaseShift(NumPts,0,0)));

    NumPts++;
  }

  TFile fOutput(OutputFolder+"fOutput_"+PotentialName+".root","recreate");
  gCk_SE.Write();
  gCk_default.Write();
  gCk_csmall.Write();
  g_csmall.Write();
  g_f2.Write();
  g_f2c.Write();
  g_F1F2.Write();
  g_F1F2_f2c.Write();
  g_XXX.Write();
  g_XXXa.Write();
  g_XXXm.Write();
  g_delta.Write();
  g_kcotdelta.Write();
  g_cotdelta.Write();
  g_tgdelta.Write();
  g_sindelta.Write();
  g_f0d0.Write();
  g_Vmu.Write();

}

int OTHERTASKS(int narg, char** ARGS){
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
    Ledni_SmallRad("NF48");
    Ledni_SmallRad("NSC97b");
    Ledni_SmallRad("emma");
    Ledni_SmallRad("custom");
    Ledni_SmallRad("Toy1");
    Ledni_SmallRad("ND46");
    //StableDisto_Test();
    //Andi_pDminus_1();
    //Fast_Bootstrap_Example();
    //MakeSmoothAngularSourceDisto(NULL);

    //const double RAD = 5.11;
    //printf("u(%.2f) = %.4f\n",RAD,Evaluate_d_u(RAD));
    //printf(" w(%.2f) = %.4f\n",RAD,Evaluate_d_w(RAD));

    return 0;
}
