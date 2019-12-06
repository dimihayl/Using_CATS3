
#include "OtherTasks.h"

#include "CATS.h"
#include "CommonAnaFunctions.h"
#include "DLM_Potentials.h"

#include "TGraph.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TNtuple.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TNtuple.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TFractionFitter.h"

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
    TString DataFolder = "/home/dmihaylov/Downloads/Temp/VALE_DEC19/MyResults/Data/CF/Norm024034/";
    TString McFolder = "/home/dmihaylov/Downloads/Temp/VALE_DEC19/MyResults/MC/CF/AOD_Trains/Norm024034/";
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
double dimi_fraction_fitter(double* x, double* par){
    return  (par[0]*dlm_ff_t1->Eval(x)+par[1]*dlm_ff_t2->Eval(x))+
            (par[2]+par[3]*0.5*(*x)+par[4]*(*x)*(*x));
}

void ALL_CorrectedMC_EXP(){
    TString OutputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/OtherTasks/ALL_CorrectedMC_EXP/";
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
    DataFolder[0] = "/home/dmihaylov/Downloads/Temp/VALE_DEC19/MyResults/Fast_BBar/Data/";
    DataFolder[1] = "/home/dmihaylov/Downloads/Temp/VALE_DEC19/MyResults/Fast_BBar/Data/";
    DataFolder[2] = "/home/dmihaylov/Downloads/Temp/VALE_DEC19/MyResults/Fast_BBar/Data/";
    DataFolder[3] = "/home/dmihaylov/Downloads/Temp/VALE_DEC19/MyResults/Fast_BBar/Data/";
    DataFolder[4] = "/home/dmihaylov/Downloads/Temp/VALE_DEC19/MyResults/Fast_BBar/Data/";
    DataFolder[5] = "/home/dmihaylov/Downloads/Temp/VALE_DEC19/MyResults/Fast_BBar/Data/";
    DataFolder[6] = "/home/dmihaylov/Downloads/Temp/VALE_DEC19/MyResults/Fast_BBar/Data/";
    DataFolder[7] = "/home/dmihaylov/Downloads/Temp/VALE_DEC19/MyResults/Fast_BBar/Data/";

    TString* McFolder = new TString[NumSpecies];
    McFolder[0] = "/home/dmihaylov/Downloads/Temp/VALE_DEC19/MyResults/Fast_BBar/Trains_MCAOD/Norm018028/";
    McFolder[1] = "/home/dmihaylov/Downloads/Temp/VALE_DEC19/MyResults/Fast_BBar/Trains_MCAOD/Norm018028/";
    McFolder[2] = "/home/dmihaylov/Downloads/Temp/VALE_DEC19/MyResults/Fast_BBar/Trains_MCAOD/Norm018028/";
    McFolder[3] = "/home/dmihaylov/Downloads/Temp/VALE_DEC19/MyResults/Fast_BBar/Trains_MCAOD/Norm018028/";
    McFolder[4] = "/home/dmihaylov/Downloads/Temp/VALE_DEC19/MyResults/Fast_BBar/Trains_MCAOD/Norm018028/";
    McFolder[5] = "/home/dmihaylov/Downloads/Temp/VALE_DEC19/MyResults/Fast_BBar/Trains_MCAOD/Norm018028/";
    McFolder[6] = "/home/dmihaylov/Downloads/Temp/VALE_DEC19/MyResults/Fast_BBar/Trains_MCAOD/Norm018028/";
    McFolder[7] = "/home/dmihaylov/Downloads/Temp/VALE_DEC19/MyResults/Fast_BBar/Trains_MCAOD/Norm018028/";

    TString* DataFile = new TString[NumSpecies];
    DataFile[0] = "CFOutput_pp_8.root";
    DataFile[1] = "CFOutput_pL_8.root";
    DataFile[2] = "CFOutput_pXi_8.root";
    DataFile[3] = "CFOutput_LL_8.root";
    DataFile[4] = "CFOutput_pAp_8.root";
    DataFile[5] = "CFOutput_pAL_8.root";
    DataFile[6] = "CFOutput_pAXi_8.root";
    DataFile[7] = "CFOutput_LAL_8.root";

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
    RescaleFactor[0] = 4;
    RescaleFactor[1] = 1;
    RescaleFactor[2] = 1;
    RescaleFactor[3] = 1;
    RescaleFactor[4] = 4;
    RescaleFactor[5] = 1;
    RescaleFactor[6] = 1;
    RescaleFactor[7] = 1;

    TString* DataSeList1 = new TString[NumSpecies];
    for(unsigned uSpec=0; uSpec<NumSpecies; uSpec++) DataSeList1[uSpec] = "PairDist";

    TString* DataSeList2 = new TString[NumSpecies];
    for(unsigned uSpec=0; uSpec<NumSpecies; uSpec++) DataSeList2[uSpec] = "PairReweighted";

    TString* McSeList1 = new TString[NumSpecies];
    for(unsigned uSpec=0; uSpec<NumSpecies; uSpec++) McSeList1[uSpec] = "PairDist";

    TString* McSeList2 = new TString[NumSpecies];
    for(unsigned uSpec=0; uSpec<NumSpecies; uSpec++) McSeList2[uSpec] = "PairReweighted";

    TString* DataSe = new TString[NumSpecies];
    DataSe[0] = "SEDist_Particle0_Particle0_clone_Shifted_FixShifted_Reweighted";
    DataSe[1] = "SEDist_Particle0_Particle2_clone_Shifted_FixShifted_Rebinned_4_Reweighted";
    DataSe[2] = "SEDist_Particle0_Particle4_clone_Shifted_FixShifted_Rebinned_4_Reweighted";
    DataSe[3] = "SEDist_Particle2_Particle2_clone_Shifted_FixShifted_Rebinned_4_Reweighted";
    DataSe[4] = "SEDist_Particle0_Particle1_clone_Shifted_FixShifted_Reweighted";
    DataSe[5] = "SEDist_Particle0_Particle3_clone_Shifted_FixShifted_Rebinned_4_Reweighted";
    DataSe[6] = "SEDist_Particle0_Particle5_clone_Shifted_FixShifted_Rebinned_4_Reweighted";
    DataSe[7] = "SEDist_Particle2_Particle3_clone_Shifted_FixShifted_Rebinned_4_Reweighted";

    TString* McSe = new TString[NumSpecies];
    McSe[0] = "SEDist_Particle0_Particle0_clone_Shifted_FixShifted_Reweighted";
    McSe[1] = "SEDist_Particle0_Particle2_clone_Shifted_FixShifted_Rebinned_4_Reweighted";
    McSe[2] = "SEDist_Particle0_Particle4_clone_Shifted_FixShifted_Rebinned_4_Reweighted";
    McSe[3] = "SEDist_Particle2_Particle2_clone_Shifted_FixShifted_Rebinned_4_Reweighted";
    McSe[4] = "SEDist_Particle0_Particle1_clone_Shifted_FixShifted_Reweighted";
    McSe[5] = "SEDist_Particle0_Particle3_clone_Shifted_FixShifted_Rebinned_4_Reweighted";
    McSe[6] = "SEDist_Particle0_Particle5_clone_Shifted_FixShifted_Rebinned_4_Reweighted";
    McSe[7] = "SEDist_Particle2_Particle3_clone_Shifted_FixShifted_Rebinned_4_Reweighted";

    TString* DataMe = new TString[NumSpecies];
    DataMe[0] = "MEDist_Particle0_Particle0_clone_Shifted_FixShifted_Reweighted";
    DataMe[1] = "MEDist_Particle0_Particle2_clone_Shifted_FixShifted_Rebinned_4_Reweighted";
    DataMe[2] = "MEDist_Particle0_Particle4_clone_Shifted_FixShifted_Rebinned_4_Reweighted";
    DataMe[3] = "MEDist_Particle2_Particle2_clone_Shifted_FixShifted_Rebinned_4_Reweighted";
    DataMe[4] = "MEDist_Particle0_Particle1_clone_Shifted_FixShifted_Reweighted";
    DataMe[5] = "MEDist_Particle0_Particle3_clone_Shifted_FixShifted_Rebinned_4_Reweighted";
    DataMe[6] = "MEDist_Particle0_Particle5_clone_Shifted_FixShifted_Rebinned_4_Reweighted";
    DataMe[7] = "MEDist_Particle2_Particle3_clone_Shifted_FixShifted_Rebinned_4_Reweighted";

    TString* McMe = new TString[NumSpecies];
    McMe[0] = "MEDist_Particle0_Particle0_clone_Shifted_FixShifted_Reweighted";
    McMe[1] = "MEDist_Particle0_Particle2_clone_Shifted_FixShifted_Rebinned_4_Reweighted";
    McMe[2] = "MEDist_Particle0_Particle4_clone_Shifted_FixShifted_Rebinned_4_Reweighted";
    McMe[3] = "MEDist_Particle2_Particle2_clone_Shifted_FixShifted_Rebinned_4_Reweighted";
    McMe[4] = "MEDist_Particle0_Particle1_clone_Shifted_FixShifted_Reweighted";
    McMe[5] = "MEDist_Particle0_Particle3_clone_Shifted_FixShifted_Rebinned_4_Reweighted";
    McMe[6] = "MEDist_Particle0_Particle5_clone_Shifted_FixShifted_Rebinned_4_Reweighted";
    McMe[7] = "MEDist_Particle2_Particle3_clone_Shifted_FixShifted_Rebinned_4_Reweighted";

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
        TList* list1_tmp;
        TList* list2_tmp;
        TH1F* hSe_tmp;
        TH1F* hMe_tmp;
printf("Opening file %s%s\n",DataFolder[uSpec].Data(),DataFile[uSpec].Data());
        TFile datafile(DataFolder[uSpec]+DataFile[uSpec]);
        hCk_tmp = (TH1F*)datafile.Get(DataHisto[uSpec]);
//printf("hCk_tmp->GetNbinsX()=%u\n",hCk_tmp->GetNbinsX());
        hCk_tmp->Rebin(RescaleFactor[uSpec]); hCk_tmp->Scale(1./double(RescaleFactor[uSpec]));
//printf("hCk_tmp->GetNbinsX()=%u\n",hCk_tmp->GetNbinsX());
        list1_tmp = (TList*)datafile.Get(DataSeList1[uSpec]);
        list2_tmp = (TList*)list1_tmp->FindObject(DataSeList2[uSpec]);
        hSe_tmp = (TH1F*)list2_tmp->FindObject(DataSe[uSpec]);
        hSe_tmp->Rebin(RescaleFactor[uSpec]); hSe_tmp->Scale(1./hSe_tmp->Integral(),"width");
        hMe_tmp = (TH1F*)list2_tmp->FindObject(DataMe[uSpec]);
        hMe_tmp->Rebin(RescaleFactor[uSpec]); hMe_tmp->Scale(1./hMe_tmp->Integral(),"width");
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

        TFile mcfile(McFolder[uSpec]+McFile[uSpec]);
        hCk_tmp = (TH1F*)mcfile.Get(McHisto[uSpec]);
        hCk_tmp->Rebin(RescaleFactor[uSpec]); hCk_tmp->Scale(1./double(RescaleFactor[uSpec]));
        list1_tmp = (TList*)mcfile.Get(McSeList1[uSpec]);
        list2_tmp = (TList*)list1_tmp->FindObject(McSeList2[uSpec]);
        hSe_tmp = (TH1F*)list2_tmp->FindObject(McSe[uSpec]);
        hSe_tmp->Rebin(RescaleFactor[uSpec]); hSe_tmp->Scale(1./hSe_tmp->Integral(),"width");
        hMe_tmp = (TH1F*)list2_tmp->FindObject(McMe[uSpec]);
        hMe_tmp->Rebin(RescaleFactor[uSpec]); hMe_tmp->Scale(1./hMe_tmp->Integral(),"width");
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
                                       dlm_RatioSe_Data_Mc[uSpec][uSpec].GetNbins(),dlm_RatioSe_Data_Mc[uSpec][uSpec].GetLowEdge(0),dlm_RatioSe_Data_Mc[uSpec][uSpec].GetUpEdge(0));
        histo_RatioMe_Data_Mc[uSpec][uSpec2] = new TH1F("histo_RatioMe_Data:"+SpeciesName[uSpec]+"_Mc:"+SpeciesName[uSpec2],
                                                        "histo_RatioMe_Data:"+SpeciesName[uSpec]+"_Mc:"+SpeciesName[uSpec2],
                                       dlm_RatioMe_Data_Mc[uSpec][uSpec].GetNbins(),dlm_RatioMe_Data_Mc[uSpec][uSpec].GetLowEdge(0),dlm_RatioMe_Data_Mc[uSpec][uSpec].GetUpEdge(0));
        histo_RatioCk_Data_Mc[uSpec][uSpec2] = new TH1F("histo_RatioCk_Data:"+SpeciesName[uSpec]+"_Mc:"+SpeciesName[uSpec2],
                                                        "histo_RatioCk_Data:"+SpeciesName[uSpec]+"_Mc:"+SpeciesName[uSpec2],
                                       dlm_RatioCk_Data_Mc[uSpec][uSpec].GetNbins(),dlm_RatioCk_Data_Mc[uSpec][uSpec].GetLowEdge(0),dlm_RatioCk_Data_Mc[uSpec][uSpec].GetUpEdge(0));
        histo_RatioSe_Data_Data[uSpec][uSpec2] = new TH1F("histo_RatioSe_Data:"+SpeciesName[uSpec]+"_Data:"+SpeciesName[uSpec2],
                                                        "histo_RatioSe_Data:"+SpeciesName[uSpec]+"_Data:"+SpeciesName[uSpec2],
                                       dlm_RatioSe_Data_Data[uSpec][uSpec].GetNbins(),dlm_RatioSe_Data_Data[uSpec][uSpec].GetLowEdge(0),dlm_RatioSe_Data_Data[uSpec][uSpec].GetUpEdge(0));
        histo_RatioMe_Data_Data[uSpec][uSpec2] = new TH1F("histo_RatioMe_Data:"+SpeciesName[uSpec]+"_Data:"+SpeciesName[uSpec2],
                                                        "histo_RatioMe_Data:"+SpeciesName[uSpec]+"_Data:"+SpeciesName[uSpec2],
                                       dlm_RatioMe_Data_Data[uSpec][uSpec].GetNbins(),dlm_RatioMe_Data_Data[uSpec][uSpec].GetLowEdge(0),dlm_RatioMe_Data_Data[uSpec][uSpec].GetUpEdge(0));
        histo_RatioCk_Data_Data[uSpec][uSpec2] = new TH1F("histo_RatioCk_Data:"+SpeciesName[uSpec]+"_Data:"+SpeciesName[uSpec2],
                                                        "histo_RatioCk_Data:"+SpeciesName[uSpec]+"_Data:"+SpeciesName[uSpec2],
                                       dlm_RatioCk_Data_Data[uSpec][uSpec].GetNbins(),dlm_RatioCk_Data_Data[uSpec][uSpec].GetLowEdge(0),dlm_RatioCk_Data_Data[uSpec][uSpec].GetUpEdge(0));
        histo_RatioSe_Mc_Mc[uSpec][uSpec2] = new TH1F("histo_RatioSe_Mc:"+SpeciesName[uSpec]+"_Mc:"+SpeciesName[uSpec2],
                                                        "histo_RatioSe_Mc:"+SpeciesName[uSpec]+"_Mc:"+SpeciesName[uSpec2],
                                       dlm_RatioSe_Mc_Mc[uSpec][uSpec].GetNbins(),dlm_RatioSe_Mc_Mc[uSpec][uSpec].GetLowEdge(0),dlm_RatioSe_Mc_Mc[uSpec][uSpec].GetUpEdge(0));
        histo_RatioMe_Mc_Mc[uSpec][uSpec2] = new TH1F("histo_RatioMe_Mc:"+SpeciesName[uSpec]+"_Mc:"+SpeciesName[uSpec2],
                                                        "histo_RatioMe_Mc:"+SpeciesName[uSpec]+"_Mc:"+SpeciesName[uSpec2],
                                       dlm_RatioMe_Mc_Mc[uSpec][uSpec].GetNbins(),dlm_RatioMe_Mc_Mc[uSpec][uSpec].GetLowEdge(0),dlm_RatioMe_Mc_Mc[uSpec][uSpec].GetUpEdge(0));
        histo_RatioCk_Mc_Mc[uSpec][uSpec2] = new TH1F("histo_RatioCk_Mc:"+SpeciesName[uSpec]+"_Mc:"+SpeciesName[uSpec2],
                                                        "histo_RatioCk_Mc:"+SpeciesName[uSpec]+"_Mc:"+SpeciesName[uSpec2],
                                       dlm_RatioCk_Mc_Mc[uSpec][uSpec].GetNbins(),dlm_RatioCk_Mc_Mc[uSpec][uSpec].GetLowEdge(0),dlm_RatioCk_Mc_Mc[uSpec][uSpec].GetUpEdge(0));

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
            histo_DataCk[uSpec]->SetBinError(uBin+1,dlm_DataCk[uSpec].GetBinError(uBin));
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

    //FIT pAp
    dlm_ff_t1 = &dlm_McCk[0];
    dlm_ff_t2 = &dlm_McCk[4];
    TF1* fit_pAp = new TF1("fit_pAp",dimi_fraction_fitter,350,3500,5);

    fit_pAp->SetParameter(0,0.7);
    fit_pAp->SetParLimits(0,0.3,1.1);

    fit_pAp->SetParameter(1,0.15);
    fit_pAp->SetParLimits(1,0.05,0.35);

    fit_pAp->SetParameter(2,0.0);
    fit_pAp->SetParLimits(2,-0.1,0.1);

    fit_pAp->SetParameter(3,0.0);
    fit_pAp->SetParLimits(3,-1e-4,1e-4);

    fit_pAp->SetParameter(4,0.0);
    fit_pAp->SetParLimits(4,-1e-7,1e-7);

    //fit_pAp->FixParameter(0,0);
    //fit_pAp->FixParameter(1,1);
    //fit_pAp->FixParameter(2,0);
    //fit_pAp->FixParameter(3,0);
    //fit_pAp->FixParameter(4,0);

    histo_DataCk[4]->Fit(fit_pAp,"S, N, R, M");
    TF1* fitPert_pAp = new TF1("fitPert_pAp","[0]+[1]*0.5*x+[2]*x*x",350,3200);
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
    TH1F* hCorrected_pAp = new TH1F("hCorrected_pAp","hCorrected_pAp",histo_DataCk[4]->GetNbinsX(),
                                    histo_DataCk[4]->GetBinLowEdge(1),histo_DataCk[4]->GetXaxis()->GetBinUpEdge(histo_DataCk[4]->GetNbinsX()));
    for(unsigned uBin=1; uBin<=histo_DataCk[4]->GetNbinsX(); uBin++){
        double MOM = hCorrected_pAp->GetBinCenter(uBin);
        hCorrected_pAp->SetBinContent(uBin,histo_DataCk[4]->GetBinContent(uBin)/fitFull_pAp->Eval(MOM));
        hCorrected_pAp->SetBinError(uBin,histo_DataCk[4]->GetBinError(uBin)/fitFull_pAp->Eval(MOM));
    }
    hCorrected_pAp->Write();


    //FIT pAp


    /////////////////////////////////////////////////////////

}


int OTHERTASKS(int narg, char** ARGS){
    //pp_CompareToNorfolk();
    //pp_pL_CorrectedMC_EXP();
    ALL_CorrectedMC_EXP();
}
