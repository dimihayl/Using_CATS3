
#include "SourceStudies.h"
#include "CommonAnaFunctions.h"
#include "CATS.h"
#include "CATStools.h"
#include "CATSconstants.h"
#include "DLM_Source.h"
#include "DLM_Potentials.h"
#include "DLM_Histo.h"
#include "DLM_Random.h"
#include "FemtoBoyzScripts.h"
#include "DLM_SubPads.h"
#include "DLM_RootWrapper.h"
#include "EnvVars.h"
#include "DLM_Ck.h"
#include "DLM_CkDecomposition.h"
#include "DLM_MathFunctions.h"
#include "DLM_HistoAnalysis.h"
#include "DLM_CkModels.h"
#include "CECA.h"
#include "TREPNI.h"

#include "TRandom3.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TStyle.h"
#include "TLatex.h"

//for the source paper, we compare C(k) and S(r) of pp and pL based on
//Gaussian and Gaussian+Reso sources
void CompareCkAndSr(){

    const TString OutputFolder = "/mnt/Ubuntu_Data/CernBox/Sync/Plots/SourcePaper/PLB_Review2/";

    //const double R_core = 1.00;
    //const double R_pp = 1.32;
    //const double R_pL = 1.48;

    const double R_core = 1.20;
    const double R_pp = 1.280;
    const double R_pL = 1.414;
    //works only for p-p at the moment
    const bool Use_3body = true;
    //the error with the minus sign for pp
    const bool CorrectSign = false;

    const unsigned kNumBins=100;
    const double kMin = 0;
    const double kMax = 200;

    const unsigned rNumBins=100;
    const double rMin = 0;
    const double rMax = 10;
    DLM_Histo<float> rDummy;
    rDummy.SetUp(1);
    rDummy.SetUp(0,rNumBins,rMin,rMax);
    rDummy.Initialize();
    TH1F* hSr_pp_CoreReso = Convert_DlmHisto_TH1F(&rDummy,"hSr_pp_CoreReso");

    DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder("/mnt/Ubuntu_Data/CernBox/Sync/CatsFiles");

    CATS AB_pp;
    AB_pp.SetMomBins(kNumBins,kMin,kMax);
    AnalysisObject.SetUpCats_pp(AB_pp,"AV18","Gauss");
    AB_pp.SetAnaSource(0,R_pp);
    AB_pp.KillTheCat();
    TGraph gCk_pp_Gauss; gCk_pp_Gauss.SetName("gCk_pp_Gauss"); gCk_pp_Gauss.Set(kNumBins);
    for(unsigned uBin=0; uBin<kNumBins; uBin++){
        gCk_pp_Gauss.SetPoint(uBin,AB_pp.GetMomentum(uBin),AB_pp.GetCorrFun(uBin));
    }
    TGraph gSr_pp_Gauss; gSr_pp_Gauss.SetName("gSr_pp_Gauss"); gSr_pp_Gauss.Set(kNumBins);
    for(unsigned uBin=0; uBin<rNumBins; uBin++){
        double RAD = rDummy.GetBinCenter(0,uBin);
        //gSr_pp_Gauss.SetPoint(uBin,rDummy.GetBinCenter(0,uBin),AB_pp.EvaluateTheSource(0,rDummy.GetBinCenter(0,uBin),0)/(4.*TMath::Pi()*RAD*RAD));
        gSr_pp_Gauss.SetPoint(uBin,rDummy.GetBinCenter(0,uBin),AB_pp.EvaluateTheSource(0,rDummy.GetBinCenter(0,uBin),0));
    }
    AB_pp.SetAnaSource(0,R_core);
    AB_pp.KillTheCat();
    TGraph gCk_pp_Core; gCk_pp_Core.SetName("gCk_pp_Core"); gCk_pp_Gauss.Set(kNumBins);
    for(unsigned uBin=0; uBin<kNumBins; uBin++){
        gCk_pp_Core.SetPoint(uBin,AB_pp.GetMomentum(uBin),AB_pp.GetCorrFun(uBin));
    }
    TGraph gSr_pp_Core; gSr_pp_Core.SetName("gSr_pp_Core"); gSr_pp_Core.Set(kNumBins);
    for(unsigned uBin=0; uBin<rNumBins; uBin++){
        double RAD = rDummy.GetBinCenter(0,uBin);
        //gSr_pp_Core.SetPoint(uBin,rDummy.GetBinCenter(0,uBin),AB_pp.EvaluateTheSource(0,rDummy.GetBinCenter(0,uBin),0)/(4.*TMath::Pi()*RAD*RAD));
        gSr_pp_Core.SetPoint(uBin,rDummy.GetBinCenter(0,uBin),AB_pp.EvaluateTheSource(0,rDummy.GetBinCenter(0,uBin),0));
    }

    //AnalysisObject.SetUpCats_pp(AB_pp,"AV18","McGauss_Reso");
    CATS AB_pp_Epos;
    AB_pp_Epos.SetMomBins(kNumBins,kMin,kMax);
printf("Set up p-p EPOS source...\n");
    AnalysisObject.SetUpCats_pp(AB_pp_Epos,"AV18","McGauss_ResoTM",0,Use_3body?204:CorrectSign?-202:202);
printf("Done!\n");
    AB_pp_Epos.SetAnaSource(0,R_core);
    AB_pp_Epos.KillTheCat();
    TGraph gCk_pp_CoreReso; gCk_pp_CoreReso.SetName("gCk_pp_CoreReso"); gCk_pp_CoreReso.Set(kNumBins);
    for(unsigned uBin=0; uBin<kNumBins; uBin++){
        gCk_pp_CoreReso.SetPoint(uBin,AB_pp_Epos.GetMomentum(uBin),AB_pp_Epos.GetCorrFun(uBin));
    }
    TGraph gSr_pp_CoreReso; gSr_pp_CoreReso.SetName("gSr_pp_CoreReso"); gSr_pp_CoreReso.Set(kNumBins);
    for(unsigned uBin=0; uBin<rNumBins; uBin++){
        double RAD = rDummy.GetBinCenter(0,uBin);
        double SVAL = AB_pp_Epos.EvaluateTheSource(0,rDummy.GetBinCenter(0,uBin),0);
        //gSr_pp_CoreReso.SetPoint(uBin,rDummy.GetBinCenter(0,uBin),AB_pp.EvaluateTheSource(0,rDummy.GetBinCenter(0,uBin),0)/(4.*TMath::Pi()*RAD*RAD));
        gSr_pp_CoreReso.SetPoint(uBin,rDummy.GetBinCenter(0,uBin),SVAL);
        hSr_pp_CoreReso->SetBinContent(uBin+1,SVAL);
        hSr_pp_CoreReso->SetBinError(uBin+1,0.1*(0.1-0.283*SVAL));
    }
    TF1* fPP_Gauss = new TF1(TString::Format("fPP_Gauss"),GaussSourceTF1,rMin,rMax,1);
    fPP_Gauss->SetParameter(0,1.0);
    fPP_Gauss->SetParLimits(0,0.5,2.0);
    hSr_pp_CoreReso->Fit(fPP_Gauss,"S, N, R, M");
    printf("pp effective Gaussian source size: %.3f fm\n",fPP_Gauss->GetParameter(0));

    CATS AB_pL;
    AB_pL.SetMomBins(kNumBins,kMin,kMax);
    AnalysisObject.SetUpCats_pL(AB_pL,"NLO_Coupled_S","Gauss");
    AB_pL.SetAnaSource(0,R_pL);
    AB_pL.KillTheCat();
    TGraph gCk_pL_Gauss; gCk_pL_Gauss.SetName("gCk_pL_Gauss"); gCk_pL_Gauss.Set(kNumBins);
    for(unsigned uBin=0; uBin<kNumBins; uBin++){
        gCk_pL_Gauss.SetPoint(uBin,AB_pL.GetMomentum(uBin),AB_pL.GetCorrFun(uBin));
    }
    TGraph gSr_pL_Gauss; gSr_pL_Gauss.SetName("gSr_pL_Gauss"); gSr_pL_Gauss.Set(kNumBins);
    for(unsigned uBin=0; uBin<rNumBins; uBin++){
        double RAD = rDummy.GetBinCenter(0,uBin);
        //gSr_pL_Gauss.SetPoint(uBin,rDummy.GetBinCenter(0,uBin),AB_pL.EvaluateTheSource(0,rDummy.GetBinCenter(0,uBin),0)/(4.*TMath::Pi()*RAD*RAD));
        gSr_pL_Gauss.SetPoint(uBin,rDummy.GetBinCenter(0,uBin),AB_pL.EvaluateTheSource(0,rDummy.GetBinCenter(0,uBin),0));
    }
    AB_pL.SetAnaSource(0,R_core);
    AB_pL.KillTheCat();
    TGraph gCk_pL_Core; gCk_pL_Core.SetName("gCk_pL_Core"); gCk_pL_Gauss.Set(kNumBins);
    for(unsigned uBin=0; uBin<kNumBins; uBin++){
        gCk_pL_Core.SetPoint(uBin,AB_pL.GetMomentum(uBin),AB_pL.GetCorrFun(uBin));
    }
    TGraph gSr_pL_Core; gSr_pL_Core.SetName("gSr_pL_Core"); gSr_pL_Core.Set(kNumBins);
    for(unsigned uBin=0; uBin<rNumBins; uBin++){
        double RAD = rDummy.GetBinCenter(0,uBin);
        //gSr_pL_Core.SetPoint(uBin,RAD,AB_pL.EvaluateTheSource(0,rDummy.GetBinCenter(0,uBin),0)/(4.*TMath::Pi()*RAD*RAD));
        gSr_pL_Core.SetPoint(uBin,RAD,AB_pL.EvaluateTheSource(0,rDummy.GetBinCenter(0,uBin),0));
    }

    CATS AB_pL_Epos;
    AB_pL_Epos.SetMomBins(kNumBins,kMin,kMax);
    //AnalysisObject.SetUpCats_pL(AB_pL_Epos,"NLO_Coupled_S","McGauss_Reso");
    AnalysisObject.SetUpCats_pL(AB_pL_Epos,"NLO_Coupled_S","McGauss_ResoTM",0,202);
    AB_pL_Epos.SetAnaSource(0,R_core);
    AB_pL_Epos.KillTheCat();
    TGraph gCk_pL_CoreReso; gCk_pL_CoreReso.SetName("gCk_pL_CoreReso"); gCk_pL_CoreReso.Set(kNumBins);
    for(unsigned uBin=0; uBin<kNumBins; uBin++){
        gCk_pL_CoreReso.SetPoint(uBin,AB_pL_Epos.GetMomentum(uBin),AB_pL_Epos.GetCorrFun(uBin));
    }
    TGraph gSr_pL_CoreReso; gSr_pL_CoreReso.SetName("gSr_pL_CoreReso"); gSr_pL_CoreReso.Set(kNumBins);
    for(unsigned uBin=0; uBin<rNumBins; uBin++){
        double RAD = rDummy.GetBinCenter(0,uBin);
        //gSr_pL_CoreReso.SetPoint(uBin,rDummy.GetBinCenter(0,uBin),AB_pL_Epos.EvaluateTheSource(0,rDummy.GetBinCenter(0,uBin),0)/(4.*TMath::Pi()*RAD*RAD));
        gSr_pL_CoreReso.SetPoint(uBin,rDummy.GetBinCenter(0,uBin),AB_pL_Epos.EvaluateTheSource(0,rDummy.GetBinCenter(0,uBin),0));
    }

    gSr_pp_Core.SetLineWidth(3);
    gSr_pp_Core.SetLineStyle(2);
    gSr_pp_Core.SetLineColor(kBlack);
    gSr_pp_Core.SetMarkerSize(0);
    gSr_pp_Core.SetMarkerStyle(0);
    gSr_pp_Core.SetMarkerColor(kBlack);

    gSr_pp_Gauss.SetLineWidth(3);
    gSr_pp_Gauss.SetLineStyle(3);
    gSr_pp_Gauss.SetLineColor(kBlue+2);
    gSr_pp_Gauss.SetMarkerSize(0);
    gSr_pp_Gauss.SetMarkerStyle(0);
    gSr_pp_Gauss.SetMarkerColor(kBlue+2);

    gSr_pp_CoreReso.SetLineWidth(0);
    gSr_pp_CoreReso.SetLineStyle(3);
    gSr_pp_CoreReso.SetLineColor(kBlue+2);
    gSr_pp_CoreReso.SetMarkerSize(1);
    gSr_pp_CoreReso.SetMarkerStyle(20);
    gSr_pp_CoreReso.SetMarkerColor(kBlue+2);

    gSr_pL_Gauss.SetLineWidth(2);
    gSr_pL_Gauss.SetLineStyle(3);
    gSr_pL_Gauss.SetLineColor(kRed+1);
    gSr_pL_Gauss.SetMarkerSize(0);
    gSr_pL_Gauss.SetMarkerStyle(0);
    gSr_pL_Gauss.SetMarkerColor(kRed+1);

    gSr_pL_CoreReso.SetLineWidth(0);
    gSr_pL_CoreReso.SetLineStyle(3);
    gSr_pL_CoreReso.SetLineColor(kRed+1);
    gSr_pL_CoreReso.SetMarkerSize(1);
    gSr_pL_CoreReso.SetMarkerStyle(24);
    gSr_pL_CoreReso.SetMarkerColor(kRed+1);

    TFile* fOut = new TFile(OutputFolder+"fOut.root","recreate");
    gCk_pp_Core.Write();
    gCk_pp_CoreReso.Write();
    gCk_pp_Gauss.Write();
    gCk_pL_Core.Write();
    gCk_pL_CoreReso.Write();
    gCk_pL_Gauss.Write();

    gSr_pp_Core.Write();
    gSr_pp_CoreReso.Write();
    hSr_pp_CoreReso->Write();
    fPP_Gauss->Write();
    gSr_pp_Gauss.Write();
    gSr_pL_Core.Write();
    gSr_pL_CoreReso.Write();
    gSr_pL_Gauss.Write();

    SetStyle();

    TCanvas* cSr = new TCanvas("cSr", "cSr", 1);
    cSr->cd(0); cSr->SetCanvasSize(1920/2, 1080/2); cSr->SetMargin(0.15,0.04,0.2,0.05);//lrbt

    TLegend* legSr = new TLegend(0.63,0.625,0.93,0.925);//lbrt
    legSr->SetName(TString::Format("legSr"));
    legSr->SetBorderSize(0);
    legSr->SetTextSize(0.04);
    legSr->AddEntry(&gSr_pp_Core,TString::Format("Gauss r(core)=%.2f fm",R_core));
    legSr->AddEntry(&gSr_pp_CoreReso,"Core + p#minusp resonances");
    legSr->AddEntry(&gSr_pp_Gauss,TString::Format("Gauss r(p#minusp)=%.2f fm",R_pp));
    legSr->AddEntry(&gSr_pL_CoreReso,"Core + p#minus#Lambda resonances");
    legSr->AddEntry(&gSr_pL_Gauss,TString::Format("Gauss r(p#minus#Lambda)=%.2f fm",R_pL));


    TH1F* hSrAxis = new TH1F("hSrAxis", "hSrAxis", rNumBins, rMin, rMax);
    hSrAxis->SetStats(false);
    hSrAxis->SetTitle("");
    hSrAxis->GetXaxis()->SetLabelSize(0.065);
    hSrAxis->GetXaxis()->SetTitle("r* (fm)");
    //hSrAxis->GetXaxis()->CenterTitle();
    hSrAxis->GetXaxis()->SetTitleOffset(1.15);
    hSrAxis->GetXaxis()->SetLabelOffset(0.02);
    hSrAxis->GetXaxis()->SetTitleSize(0.075);
    hSrAxis->GetYaxis()->SetLabelSize(0.065);
    //hAxis->GetYaxis()->SetTitle("(N per bin)x10^{-3}");
    hSrAxis->GetYaxis()->SetTitle("4#pir*^{2}S(r*) (1/fm)");
    //hSrAxis->GetYaxis()->CenterTitle();
    hSrAxis->GetYaxis()->SetTitleOffset(0.90);
    hSrAxis->GetYaxis()->SetTitleSize(0.075);
    //hSrAxis->GetYaxis()->SetLimits(0, hNkSe->GetBinContent(hNkSe->GetMaximumBin())*1.2);
    hSrAxis->GetYaxis()->SetRangeUser(1e-3,1);
    //hSrAxis->GetYaxis()->SetRangeUser(0,0.5);

    hSrAxis->Draw("axis");
    gSr_pp_Core.Draw("same,C");
    gSr_pp_Gauss.Draw("same,C");
    gSr_pL_Gauss.Draw("same,C");
    gSr_pp_CoreReso.Draw("same,P");
    gSr_pL_CoreReso.Draw("same,P");
    legSr->Draw("same");

    cSr->SetLogy(true);
    cSr->SetGridx(true);
    cSr->SetGridy(true);

    cSr->Write();
    cSr->SaveAs(OutputFolder+"cSr.png");
    cSr->SaveAs(OutputFolder+"cSr.pdf");

    delete legSr;
    delete hSrAxis;
    delete cSr;
    delete fOut;
}

void AverageResoApprox_pp(){

    const unsigned NumMomBins = 150;
    const double kMin = 0;
    const double kMax = 300;

    const double Mass_pi0 = 134.9766;
    const double Mass_pic = 139.57018;
    const double Mass_p = 938.272;
    const double Mass_L = 1115.683;

    double PotPars1S0[8]={NN_AV18,v18_Coupled3P2,1,1,1,0,0,0};
    double PotPars3P0[8]={NN_AV18,v18_Coupled3P2,1,1,1,1,1,0};
    double PotPars3P1[8]={NN_AV18,v18_Coupled3P2,1,1,1,1,1,1};
    double PotPars3P2[8]={NN_AV18,v18_Coupled3P2,1,1,1,1,1,2};
    CATSparameters cSorPars_Def(CATSparameters::tSource,2,true);cSorPars_Def.SetParameter(0,1);cSorPars_Def.SetParameter(1,2);
    CATSparameters cPotPars1S0_Def(CATSparameters::tPotential,8,true);cPotPars1S0_Def.SetParameters(PotPars1S0);
    CATSparameters cPotPars3P0_Def(CATSparameters::tPotential,8,true);cPotPars3P0_Def.SetParameters(PotPars3P0);
    CATSparameters cPotPars3P1_Def(CATSparameters::tPotential,8,true);cPotPars3P1_Def.SetParameters(PotPars3P1);
    CATSparameters cPotPars3P2_Def(CATSparameters::tPotential,8,true);cPotPars3P2_Def.SetParameters(PotPars3P2);

    CATSparameters cSorPars_Imp(CATSparameters::tSource,2,true);cSorPars_Imp.SetParameter(0,1);cSorPars_Imp.SetParameter(1,2);
    CATSparameters cPotPars1S0_Imp(CATSparameters::tPotential,8,true);cPotPars1S0_Imp.SetParameters(PotPars1S0);
    CATSparameters cPotPars3P0_Imp(CATSparameters::tPotential,8,true);cPotPars3P0_Imp.SetParameters(PotPars3P0);
    CATSparameters cPotPars3P1_Imp(CATSparameters::tPotential,8,true);cPotPars3P1_Imp.SetParameters(PotPars3P1);
    CATSparameters cPotPars3P2_Imp(CATSparameters::tPotential,8,true);cPotPars3P2_Imp.SetParameters(PotPars3P2);

    DLM_CleverMcLevyReso CleverSource_Def;
    CleverSource_Def.InitStability(21,1,2);
    CleverSource_Def.InitScale(38,0.15,2.0);
    CleverSource_Def.InitRad(257,0,64);
    CleverSource_Def.InitType(2);
    CleverSource_Def.InitReso(0,2);
    CleverSource_Def.InitReso(1,2);
    CleverSource_Def.SetUpReso(0,0,0.6422,1361.52,1.65,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);
    CleverSource_Def.SetUpReso(1,0,0.6422,1361.52,1.65,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);
    DLM_Histo<double>* HISTO_AVG = ConvertThetaAngleHisto(
    "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/EPOS_HM_4PI/pp/DimiPhi_pp.root","h_rkAngle_Mom2",400,580);
    CleverSource_Def.SetUpResoEmission(0,0,HISTO_AVG);
    CleverSource_Def.SetUpResoEmission(1,0,HISTO_AVG);
    CleverSource_Def.InitNumMcIter(100000);

    DLM_CleverMcLevyReso CleverSource_Imp;
    CleverSource_Imp.InitStability(21,1,2);
    CleverSource_Imp.InitScale(38,0.15,2.0);
    CleverSource_Imp.InitRad(257,0,64);
    CleverSource_Imp.InitType(2);
    CleverSource_Imp.InitReso(0,2);
    CleverSource_Imp.InitReso(1,2);
    CleverSource_Imp.SetUpReso(0,0,0.4368,1231.98,1.67,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);
    CleverSource_Imp.SetUpReso(0,1,0.2054,1636.98,1.62,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);
    CleverSource_Imp.SetUpReso(1,0,0.4368,1231.98,1.67,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);
    CleverSource_Imp.SetUpReso(1,1,0.2054,1636.98,1.62,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);
    DLM_Histo<double>* HISTO_DELTA = ConvertThetaAngleHisto(
    "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/EPOS_HM_4PI/pp/DimiPhi_pp.root","h_rkAngle_Mom2",200,400);
    DLM_Histo<double>* HISTO_REST = ConvertThetaAngleHisto(
    "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/EPOS_HM_4PI/pp/DimiPhi_pp.root","h_rkAngle_Mom2",880,1000);
    CleverSource_Imp.SetUpResoEmission(0,0,HISTO_DELTA);
    CleverSource_Imp.SetUpResoEmission(0,1,HISTO_REST);
    CleverSource_Imp.SetUpResoEmission(1,0,HISTO_DELTA);
    CleverSource_Imp.SetUpResoEmission(1,1,HISTO_REST);
    CleverSource_Imp.InitNumMcIter(100000);

    CATS DefaultCat;
    DefaultCat.SetMomBins(NumMomBins,kMin,kMax);
    DefaultCat.SetAnaSource(CatsSourceForwarder, &CleverSource_Def, 2);
    DefaultCat.SetAnaSource(0,1.0);
    DefaultCat.SetAnaSource(1,2.0);
    DefaultCat.SetUseAnalyticSource(true);
    DefaultCat.SetMomentumDependentSource(false);
    DefaultCat.SetThetaDependentSource(false);
    DefaultCat.SetExcludeFailedBins(false);
    DefaultCat.SetQ1Q2(1);
    DefaultCat.SetPdgId(2212, 2212);
    DefaultCat.SetRedMass( 0.5*Mass_p );
    DefaultCat.SetNumChannels(4);
    DefaultCat.SetNumPW(0,1);
    DefaultCat.SetNumPW(1,2);
    DefaultCat.SetNumPW(2,2);
    DefaultCat.SetNumPW(3,2);
    DefaultCat.SetSpin(0,0);
    DefaultCat.SetSpin(1,1);
    DefaultCat.SetSpin(2,1);
    DefaultCat.SetSpin(3,1);
    DefaultCat.SetChannelWeight(0, 3./12.);
    DefaultCat.SetChannelWeight(1, 1./12.);
    DefaultCat.SetChannelWeight(2, 3./12.);
    DefaultCat.SetChannelWeight(3, 5./12.);
    DefaultCat.SetShortRangePotential(0,0,fDlmPot,cPotPars1S0_Def);
    DefaultCat.SetShortRangePotential(1,1,fDlmPot,cPotPars3P0_Def);
    DefaultCat.SetShortRangePotential(2,1,fDlmPot,cPotPars3P1_Def);
    DefaultCat.SetShortRangePotential(3,1,fDlmPot,cPotPars3P2_Def);
    DefaultCat.SetEpsilonConv(1e-8);
    DefaultCat.SetEpsilonProp(1e-8);
    DefaultCat.KillTheCat();

    CATS ImprovedCat;
    ImprovedCat.SetMomBins(NumMomBins,kMin,kMax);
    ImprovedCat.SetAnaSource(CatsSourceForwarder, &CleverSource_Imp, 2);
    ImprovedCat.SetAnaSource(0,1.0);
    ImprovedCat.SetAnaSource(1,2.0);
    ImprovedCat.SetUseAnalyticSource(true);
    ImprovedCat.SetMomentumDependentSource(false);
    ImprovedCat.SetThetaDependentSource(false);
    ImprovedCat.SetExcludeFailedBins(false);
    ImprovedCat.SetQ1Q2(1);
    ImprovedCat.SetPdgId(2212, 2212);
    ImprovedCat.SetRedMass( 0.5*Mass_p );
    ImprovedCat.SetNumChannels(4);
    ImprovedCat.SetNumPW(0,1);
    ImprovedCat.SetNumPW(1,2);
    ImprovedCat.SetNumPW(2,2);
    ImprovedCat.SetNumPW(3,2);
    ImprovedCat.SetSpin(0,0);
    ImprovedCat.SetSpin(1,1);
    ImprovedCat.SetSpin(2,1);
    ImprovedCat.SetSpin(3,1);
    ImprovedCat.SetChannelWeight(0, 3./12.);
    ImprovedCat.SetChannelWeight(1, 1./12.);
    ImprovedCat.SetChannelWeight(2, 3./12.);
    ImprovedCat.SetChannelWeight(3, 5./12.);
    ImprovedCat.SetShortRangePotential(0,0,fDlmPot,cPotPars1S0_Imp);
    ImprovedCat.SetShortRangePotential(1,1,fDlmPot,cPotPars3P0_Imp);
    ImprovedCat.SetShortRangePotential(2,1,fDlmPot,cPotPars3P1_Imp);
    ImprovedCat.SetShortRangePotential(3,1,fDlmPot,cPotPars3P2_Imp);
    ImprovedCat.SetEpsilonConv(1e-8);
    ImprovedCat.SetEpsilonProp(1e-8);
    ImprovedCat.KillTheCat();

    TFile OutputFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/SourceStudies/AverageResoApprox/OutputFile.root","recreate");

    TGraph gDefault;
    gDefault.SetName("gDefault");
    gDefault.Set(NumMomBins);
    TGraph gImproved;
    gImproved.SetName("gImproved");
    gImproved.Set(NumMomBins);

    TGraph gRatio;
    gRatio.SetName("gRatio");
    gRatio.Set(NumMomBins);

    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        gDefault.SetPoint(uBin,DefaultCat.GetMomentum(uBin),DefaultCat.GetCorrFun(uBin));
        gImproved.SetPoint(uBin,ImprovedCat.GetMomentum(uBin),ImprovedCat.GetCorrFun(uBin));
        gRatio.SetPoint(uBin,DefaultCat.GetMomentum(uBin),DefaultCat.GetCorrFun(uBin)/ImprovedCat.GetCorrFun(uBin));
    }

    gDefault.Write();
    gImproved.Write();
    gRatio.Write();
}


void AverageResoApprox_pL(){

    const unsigned NumMomBins = 150;
    const double kMin = 0;
    const double kMax = 300;

    const double Mass_pi0 = 134.9766;
    const double Mass_pic = 139.57018;
    const double Mass_p = 938.272;
    const double Mass_L = 1115.683;

    double PotPars1S0[8]={pL_UsmaniOli,0,0,0,0,0,0,0};
    double PotPars3S1[8]={pL_UsmaniOli,0,0,0,0,1,0,1};
    CATSparameters cSorPars_Def(CATSparameters::tSource,2,true);cSorPars_Def.SetParameter(0,1);cSorPars_Def.SetParameter(1,2);
    CATSparameters cPotPars1S0_Def(CATSparameters::tPotential,8,true);cPotPars1S0_Def.SetParameters(PotPars1S0);
    CATSparameters cPotPars3S1_Def(CATSparameters::tPotential,8,true);cPotPars3S1_Def.SetParameters(PotPars3S1);

    CATSparameters cSorPars_Imp(CATSparameters::tSource,2,true);cSorPars_Imp.SetParameter(0,1);cSorPars_Imp.SetParameter(1,2);
    CATSparameters cPotPars1S0_Imp(CATSparameters::tPotential,8,true);cPotPars1S0_Imp.SetParameters(PotPars1S0);
    CATSparameters cPotPars3S1_Imp(CATSparameters::tPotential,8,true);cPotPars3S1_Imp.SetParameters(PotPars3S1);

    DLM_CleverMcLevyReso CleverSource_Def;
    CleverSource_Def.InitStability(21,1,2);
    CleverSource_Def.InitScale(38,0.15,2.0);
    CleverSource_Def.InitRad(257,0,64);
    CleverSource_Def.InitType(2);
    CleverSource_Def.InitReso(0,1);
    CleverSource_Def.InitReso(1,1);
    CleverSource_Def.SetUpReso(0,0,0.6422,1361.52,1.65,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);
    CleverSource_Def.SetUpReso(1,0,0.6438,1462.93,4.69,Mass_L,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);
    DLM_Histo<double>* HISTO_P_AVG = ConvertThetaAngleHisto(
    "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/EPOS_HM_4PI/pp/DimiPhi_pp.root","h_rkAngle_Mom2",400,580);
    DLM_Histo<double>* HISTO_L_AVG = ConvertThetaAngleHisto(
    "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/EPOS_HM_4PI/pp/DimiPhi_pp.root","h_rkAngle_Mom2",270,470);
    CleverSource_Def.SetUpResoEmission(0,0,HISTO_P_AVG);
    CleverSource_Def.SetUpResoEmission(1,0,HISTO_L_AVG);
    CleverSource_Def.InitNumMcIter(100000);

    DLM_CleverMcLevyReso CleverSource_Imp;
    CleverSource_Imp.InitStability(21,1,2);
    CleverSource_Imp.InitScale(38,0.15,2.0);
    CleverSource_Imp.InitRad(257,0,64);
    CleverSource_Imp.InitType(2);
    CleverSource_Imp.InitReso(0,1);
    //CleverSource_Imp.InitReso(0,2);
    CleverSource_Imp.InitReso(1,2);
    CleverSource_Imp.SetUpReso(0,0,0.6422,1361.52,1.65,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);
    //CleverSource_Imp.SetUpReso(0,0,0.4368,1231.98,1.67,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);
    //CleverSource_Imp.SetUpReso(0,1,0.2054,1636.98,1.62,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);
    CleverSource_Imp.SetUpReso(1,0,0.4864,1384.54,5.34,Mass_L,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);
    CleverSource_Imp.SetUpReso(1,1,0.1573,1705.26,2.70,Mass_L,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);
    DLM_Histo<double>* HISTO_P_DELTA = ConvertThetaAngleHisto(
    "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/EPOS_HM_4PI/pp/DimiPhi_pp.root","h_rkAngle_Mom2",200,400);
    DLM_Histo<double>* HISTO_P_REST = ConvertThetaAngleHisto(
    "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/EPOS_HM_4PI/pp/DimiPhi_pp.root","h_rkAngle_Mom2",880,1000);
    DLM_Histo<double>* HISTO_L_SIGMA = ConvertThetaAngleHisto(
    "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/EPOS_HM_4PI/pp/DimiPhi_pp.root","h_rkAngle_Mom2",160,360);//257
    DLM_Histo<double>* HISTO_L_REST = ConvertThetaAngleHisto(
    "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/EPOS_HM_4PI/pp/DimiPhi_pp.root","h_rkAngle_Mom2",620,820);//722
    CleverSource_Imp.SetUpResoEmission(0,0,HISTO_P_AVG);
    //CleverSource_Imp.SetUpResoEmission(0,0,HISTO_P_DELTA);
    //CleverSource_Imp.SetUpResoEmission(0,1,HISTO_P_REST);
    CleverSource_Imp.SetUpResoEmission(1,0,HISTO_L_SIGMA);
    CleverSource_Imp.SetUpResoEmission(1,1,HISTO_L_REST);
    CleverSource_Imp.InitNumMcIter(100000);

    CATS DefaultCat;
    DefaultCat.SetMomBins(NumMomBins,kMin,kMax);
    DefaultCat.SetAnaSource(CatsSourceForwarder, &CleverSource_Def, 2);
    DefaultCat.SetAnaSource(0,1.0);
    DefaultCat.SetAnaSource(1,2.0);
    DefaultCat.SetUseAnalyticSource(true);
    DefaultCat.SetMomentumDependentSource(false);
    DefaultCat.SetThetaDependentSource(false);
    DefaultCat.SetExcludeFailedBins(false);
    DefaultCat.SetQ1Q2(0);
    DefaultCat.SetPdgId(2212, 3122);
    DefaultCat.SetRedMass( (Mass_p*Mass_L)/(Mass_p+Mass_L) );
    DefaultCat.SetNumChannels(2);
    DefaultCat.SetNumPW(0,1);
    DefaultCat.SetNumPW(1,1);
    DefaultCat.SetSpin(0,0);
    DefaultCat.SetSpin(1,1);
    DefaultCat.SetChannelWeight(0, 1./4.);
    DefaultCat.SetChannelWeight(1, 3./4.);
    DefaultCat.SetShortRangePotential(0,0,fDlmPot,cPotPars1S0_Def);
    DefaultCat.SetShortRangePotential(1,0,fDlmPot,cPotPars3S1_Def);
    DefaultCat.SetEpsilonConv(5e-8);
    DefaultCat.SetEpsilonProp(5e-8);
    DefaultCat.KillTheCat();

    CATS ImprovedCat;
    ImprovedCat.SetMomBins(NumMomBins,kMin,kMax);
    ImprovedCat.SetAnaSource(CatsSourceForwarder, &CleverSource_Imp, 2);
    ImprovedCat.SetAnaSource(0,1.0);
    ImprovedCat.SetAnaSource(1,2.0);
    ImprovedCat.SetUseAnalyticSource(true);
    ImprovedCat.SetMomentumDependentSource(false);
    ImprovedCat.SetThetaDependentSource(false);
    ImprovedCat.SetExcludeFailedBins(false);
    ImprovedCat.SetQ1Q2(0);
    ImprovedCat.SetPdgId(2212, 3122);
    ImprovedCat.SetRedMass( (Mass_p*Mass_L)/(Mass_p+Mass_L) );
    ImprovedCat.SetNumChannels(2);
    ImprovedCat.SetNumPW(0,1);
    ImprovedCat.SetNumPW(1,1);
    ImprovedCat.SetSpin(0,0);
    ImprovedCat.SetSpin(1,1);
    ImprovedCat.SetChannelWeight(0, 1./4.);
    ImprovedCat.SetChannelWeight(1, 3./4.);
    ImprovedCat.SetShortRangePotential(0,0,fDlmPot,cPotPars1S0_Imp);
    ImprovedCat.SetShortRangePotential(1,0,fDlmPot,cPotPars3S1_Imp);
    ImprovedCat.SetEpsilonConv(5e-8);
    ImprovedCat.SetEpsilonProp(5e-8);
    ImprovedCat.KillTheCat();

    TFile OutputFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/SourceStudies/AverageResoApprox/OutputFile_pL.root","recreate");

    TGraph gDefault;
    gDefault.SetName("gDefault");
    gDefault.Set(NumMomBins);
    TGraph gImproved;
    gImproved.SetName("gImproved");
    gImproved.Set(NumMomBins);

    TGraph gRatio;
    gRatio.SetName("gRatio");
    gRatio.Set(NumMomBins);

    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        gDefault.SetPoint(uBin,DefaultCat.GetMomentum(uBin),DefaultCat.GetCorrFun(uBin));
        gImproved.SetPoint(uBin,ImprovedCat.GetMomentum(uBin),ImprovedCat.GetCorrFun(uBin));
        gRatio.SetPoint(uBin,DefaultCat.GetMomentum(uBin),DefaultCat.GetCorrFun(uBin)/ImprovedCat.GetCorrFun(uBin));
    }

    gDefault.Write();
    gImproved.Write();
    gRatio.Write();
}

//TYPE
//0 = back-to-back
//1 = uniform
//2 = uniform back
//3 = EPOS pp and pL window from pp 400-600
//4 = EPOS pL window from pp 270-470
//5 = EPOS pL window from pXim 270-470
void TestDifferentAngularDistributions(const double& SourceSize){

    TFile* fCommonOut = new TFile(
    TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/SourceStudies/TestDifferentAngularDistributions/fCommonOut_%.2f.root",SourceSize),
    "recreate");
    const unsigned NumTypes = 6;
    TH1F* hPP_RAD = new TH1F("hPP_RAD","hPP_RAD",NumTypes,-0.5,5.5);
    TH1F* hPL_RAD = new TH1F("hPL_RAD","hPL_RAD",NumTypes,-0.5,5.5);

    for(int TYPE=0; TYPE<NumTypes; TYPE++){
        const unsigned NumMomBins = 75;
        const double kMin = 0;
        const double kMax = 300;

        const double Mass_pi0 = 134.9766;
        const double Mass_pic = 139.57018;
        const double Mass_p = 938.272;
        const double Mass_L = 1115.683;

        double ppPotPars1S0[8]={NN_AV18,v18_Coupled3P2,1,1,1,0,0,0};
        double ppPotPars3P0[8]={NN_AV18,v18_Coupled3P2,1,1,1,1,1,0};
        double ppPotPars3P1[8]={NN_AV18,v18_Coupled3P2,1,1,1,1,1,1};
        double ppPotPars3P2[8]={NN_AV18,v18_Coupled3P2,1,1,1,1,1,2};
        CATSparameters cPotPars1S0_pp(CATSparameters::tPotential,8,true);cPotPars1S0_pp.SetParameters(ppPotPars1S0);
        CATSparameters cPotPars3P0_pp(CATSparameters::tPotential,8,true);cPotPars3P0_pp.SetParameters(ppPotPars3P0);
        CATSparameters cPotPars3P1_pp(CATSparameters::tPotential,8,true);cPotPars3P1_pp.SetParameters(ppPotPars3P1);
        CATSparameters cPotPars3P2_pp(CATSparameters::tPotential,8,true);cPotPars3P2_pp.SetParameters(ppPotPars3P2);
        //CATSparameters cSorPars_pp(CATSparameters::tSource,2,true);cSorPars_pp.SetParameter(0,1);cSorPars_pp.SetParameter(1,2);

        double pL_PotPars1S0[8]={pL_UsmaniOli,0,0,0,0,0,0,0};
        double pL_PotPars3S1[8]={pL_UsmaniOli,0,0,0,0,1,0,1};
        //CATSparameters cSorPars_pL(CATSparameters::tSource,2,true);cSorPars_pL.SetParameter(0,1);cSorPars_pL.SetParameter(1,2);
        CATSparameters cPotPars1S0_pL(CATSparameters::tPotential,8,true);cPotPars1S0_pL.SetParameters(pL_PotPars1S0);
        CATSparameters cPotPars3S1_pL(CATSparameters::tPotential,8,true);cPotPars3S1_pL.SetParameters(pL_PotPars3S1);


        DLM_CleverMcLevyReso CleverSource_pp;
        CleverSource_pp.InitStability(1,2-1e-6,2+1e-6);
        CleverSource_pp.InitScale(38,0.15,2.0);
        CleverSource_pp.InitRad(257,0,64);
        CleverSource_pp.InitType(2);
        CleverSource_pp.InitReso(0,2);
        CleverSource_pp.InitReso(1,2);

        DLM_CleverMcLevyReso CleverSource_pL;
        CleverSource_pL.InitStability(1,2-1e-6,2+1e-6);
        CleverSource_pL.InitScale(38,0.15,2.0);
        CleverSource_pL.InitRad(257,0,64);
        CleverSource_pL.InitType(2);
        CleverSource_pL.InitReso(0,1);
        //CleverSource_pL_Imp.InitReso(0,2);
        CleverSource_pL.InitReso(1,2);

        DLM_Histo<double>* HISTO_P_AVG;
        DLM_Histo<double>* HISTO_L_AVG;
        if(TYPE==0){
            CleverSource_pp.SetUpReso(0,0,0.6422,1361.52,1.65,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtBackwards);
            CleverSource_pp.SetUpReso(1,0,0.6422,1361.52,1.65,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtBackwards);
            CleverSource_pL.SetUpReso(0,0,0.6422,1361.52,1.65,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtBackwards);
            CleverSource_pL.SetUpReso(1,0,0.6438,1462.93,4.69,Mass_L,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtBackwards);
        }
        else if(TYPE==1){
            CleverSource_pp.SetUpReso(0,0,0.6422,1361.52,1.65,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);
            CleverSource_pp.SetUpReso(1,0,0.6422,1361.52,1.65,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);
            CleverSource_pL.SetUpReso(0,0,0.6422,1361.52,1.65,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);
            CleverSource_pL.SetUpReso(1,0,0.6438,1462.93,4.69,Mass_L,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);
        }
        else if(TYPE==2){
            CleverSource_pp.SetUpReso(0,0,0.6422,1361.52,1.65,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandomBackwards);
            CleverSource_pp.SetUpReso(1,0,0.6422,1361.52,1.65,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandomBackwards);
            CleverSource_pL.SetUpReso(0,0,0.6422,1361.52,1.65,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandomBackwards);
            CleverSource_pL.SetUpReso(1,0,0.6438,1462.93,4.69,Mass_L,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandomBackwards);
        }
        else if(TYPE==3){
            HISTO_P_AVG = ConvertThetaAngleHisto(
            "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/EPOS_HM_4PI/pp/DimiPhi_pp.root","h_rkAngle_Mom2",400,580);
            HISTO_L_AVG = ConvertThetaAngleHisto(
            "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/EPOS_HM_4PI/pp/DimiPhi_pp.root","h_rkAngle_Mom2",400,580);
            CleverSource_pp.SetUpReso(0,0,0.6422,1361.52,1.65,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);
            CleverSource_pp.SetUpReso(1,0,0.6422,1361.52,1.65,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);
            CleverSource_pp.SetUpResoEmission(0,0,HISTO_P_AVG);
            CleverSource_pp.SetUpResoEmission(1,0,HISTO_P_AVG);
            CleverSource_pL.SetUpReso(0,0,0.6422,1361.52,1.65,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);
            CleverSource_pL.SetUpReso(1,0,0.6438,1462.93,4.69,Mass_L,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);
            CleverSource_pL.SetUpResoEmission(0,0,HISTO_P_AVG);
            CleverSource_pL.SetUpResoEmission(1,0,HISTO_L_AVG);
        }
        else if(TYPE==4){
            HISTO_P_AVG = ConvertThetaAngleHisto(
            "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/EPOS_HM_4PI/pp/DimiPhi_pp.root","h_rkAngle_Mom2",400,580);
            HISTO_L_AVG = ConvertThetaAngleHisto(
            "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/EPOS_HM_4PI/pp/DimiPhi_pp.root","h_rkAngle_Mom2",270,470);
            CleverSource_pp.SetUpReso(0,0,0.6422,1361.52,1.65,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);
            CleverSource_pp.SetUpReso(1,0,0.6422,1361.52,1.65,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);
            CleverSource_pp.SetUpResoEmission(0,0,HISTO_P_AVG);
            CleverSource_pp.SetUpResoEmission(1,0,HISTO_P_AVG);
            CleverSource_pL.SetUpReso(0,0,0.6422,1361.52,1.65,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);
            CleverSource_pL.SetUpReso(1,0,0.6438,1462.93,4.69,Mass_L,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);
            CleverSource_pL.SetUpResoEmission(0,0,HISTO_P_AVG);
            CleverSource_pL.SetUpResoEmission(1,0,HISTO_L_AVG);
        }
        else{
            HISTO_P_AVG = ConvertThetaAngleHisto(
            "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/EPOS_HM_4PI/pXi/DimiPhi_pXim.root","h_rkAngle_Mom2",400,580);
            HISTO_L_AVG = ConvertThetaAngleHisto(
            "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/EPOS_HM_4PI/pXi/DimiPhi_pXim.root","h_rkAngle_Mom2",270,470);
            CleverSource_pp.SetUpReso(0,0,0.6422,1361.52,1.65,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);
            CleverSource_pp.SetUpReso(1,0,0.6422,1361.52,1.65,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);
            CleverSource_pp.SetUpResoEmission(0,0,HISTO_P_AVG);
            CleverSource_pp.SetUpResoEmission(1,0,HISTO_P_AVG);
            CleverSource_pL.SetUpReso(0,0,0.6422,1361.52,1.65,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);
            CleverSource_pL.SetUpReso(1,0,0.6438,1462.93,4.69,Mass_L,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);
            CleverSource_pL.SetUpResoEmission(0,0,HISTO_P_AVG);
            CleverSource_pL.SetUpResoEmission(1,0,HISTO_L_AVG);
        }

        CleverSource_pL.InitNumMcIter(262144);

        CATS Kitty_pp;
        Kitty_pp.SetMomBins(NumMomBins,kMin,kMax);
        Kitty_pp.SetAnaSource(CatsSourceForwarder, &CleverSource_pp, 2);
        Kitty_pp.SetAnaSource(0,SourceSize);
        Kitty_pp.SetAnaSource(1,2.0);
        Kitty_pp.SetUseAnalyticSource(true);
        Kitty_pp.SetMomentumDependentSource(false);
        Kitty_pp.SetThetaDependentSource(false);
        Kitty_pp.SetExcludeFailedBins(false);
        Kitty_pp.SetQ1Q2(1);
        Kitty_pp.SetPdgId(2212, 2212);
        Kitty_pp.SetRedMass( 0.5*Mass_p );
        Kitty_pp.SetNumChannels(4);
        Kitty_pp.SetNumPW(0,1);
        Kitty_pp.SetNumPW(1,2);
        Kitty_pp.SetNumPW(2,2);
        Kitty_pp.SetNumPW(3,2);
        Kitty_pp.SetSpin(0,0);
        Kitty_pp.SetSpin(1,1);
        Kitty_pp.SetSpin(2,1);
        Kitty_pp.SetSpin(3,1);
        Kitty_pp.SetChannelWeight(0, 3./12.);
        Kitty_pp.SetChannelWeight(1, 1./12.);
        Kitty_pp.SetChannelWeight(2, 3./12.);
        Kitty_pp.SetChannelWeight(3, 5./12.);
        Kitty_pp.SetShortRangePotential(0,0,fDlmPot,cPotPars1S0_pp);
        Kitty_pp.SetShortRangePotential(1,1,fDlmPot,cPotPars3P0_pp);
        Kitty_pp.SetShortRangePotential(2,1,fDlmPot,cPotPars3P1_pp);
        Kitty_pp.SetShortRangePotential(3,1,fDlmPot,cPotPars3P2_pp);
        //Kitty_pp.SetEpsilonConv(5e-8);
        //Kitty_pp.SetEpsilonProp(5e-8);
        Kitty_pp.SetMaxNumThreads(4);
        Kitty_pp.KillTheCat();

        CATS Kitty_pL;
        Kitty_pL.SetMomBins(NumMomBins,kMin,kMax);
        Kitty_pL.SetAnaSource(CatsSourceForwarder, &CleverSource_pL, 2);
        Kitty_pL.SetAnaSource(0,SourceSize);
        Kitty_pL.SetAnaSource(1,2.0);
        Kitty_pL.SetUseAnalyticSource(true);
        Kitty_pL.SetMomentumDependentSource(false);
        Kitty_pL.SetThetaDependentSource(false);
        Kitty_pL.SetExcludeFailedBins(false);
        Kitty_pL.SetQ1Q2(0);
        Kitty_pL.SetPdgId(2212, 3122);
        Kitty_pL.SetRedMass( (Mass_p*Mass_L)/(Mass_p+Mass_L) );
        Kitty_pL.SetNumChannels(2);
        Kitty_pL.SetNumPW(0,1);
        Kitty_pL.SetNumPW(1,1);
        Kitty_pL.SetSpin(0,0);
        Kitty_pL.SetSpin(1,1);
        Kitty_pL.SetChannelWeight(0, 1./4.);
        Kitty_pL.SetChannelWeight(1, 3./4.);
        Kitty_pL.SetShortRangePotential(0,0,fDlmPot,cPotPars1S0_pL);
        Kitty_pL.SetShortRangePotential(1,0,fDlmPot,cPotPars3S1_pL);
        //Kitty_pL.SetEpsilonConv(5e-8);
        //Kitty_pL.SetEpsilonProp(5e-8);
        Kitty_pL.SetMaxNumThreads(4);
        Kitty_pL.KillTheCat();

        TFile* OutputFile = new TFile(
        TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/SourceStudies/TestDifferentAngularDistributions/TDAD_%.2f_%i.root",SourceSize,TYPE)
                         ,"recreate");

        TH1F* hPP_Source;
        TF1* fPP_Gauss;
        TH1F* hPL_Source;
        TF1* fPL_Gauss;

        const double rMin = 0;
        const double rMax = 8;
        const unsigned rNumBins = 128;
        const double RELERR = 0.002;
        const double ABSERR = 0.002;
        double RAD;
        double SVAL;
        hPP_Source = new TH1F("hPP_Source","hPP_Source",rNumBins,rMin,rMax);
        hPL_Source = new TH1F("hPL_Source","hPL_Source",rNumBins,rMin,rMax);
        for(unsigned uRad=0; uRad<rNumBins; uRad++){
            RAD = hPP_Source->GetBinCenter(uRad+1);
            SVAL = Kitty_pp.EvaluateTheSource(10,RAD,0);
            hPP_Source->SetBinContent(uRad+1,SVAL);
            hPP_Source->SetBinError(uRad+1,RELERR*SVAL+ABSERR);

            RAD = hPL_Source->GetBinCenter(uRad+1);
            SVAL = Kitty_pL.EvaluateTheSource(10,RAD,0);
            hPL_Source->SetBinContent(uRad+1,SVAL);
            hPL_Source->SetBinError(uRad+1,RELERR*SVAL+ABSERR);
        }

        fPP_Gauss = new TF1(TString::Format("fPP_Gauss"),GaussSourceTF1,rMin,rMax,1);
        fPP_Gauss->SetParameter(0,1.0);
        fPP_Gauss->SetParLimits(0,0.5,2.0);

        fPL_Gauss = new TF1(TString::Format("fPL_Gauss"),GaussSourceTF1,rMin,rMax,1);
        fPL_Gauss->SetParameter(0,1.0);
        fPL_Gauss->SetParLimits(0,0.5,2.0);

        hPP_Source->Fit(fPP_Gauss,"S, N, R, M");
        hPL_Source->Fit(fPL_Gauss,"S, N, R, M");

        hPP_RAD->SetBinContent(TYPE+1,fPP_Gauss->GetParameter(0));
        hPP_RAD->SetBinError(TYPE+1,fPP_Gauss->GetParError(0));

        hPL_RAD->SetBinContent(TYPE+1,fPL_Gauss->GetParameter(0));
        hPL_RAD->SetBinError(TYPE+1,fPL_Gauss->GetParError(0));



        TGraph gPP;
        gPP.SetName("gPP");
        gPP.Set(NumMomBins);

        TGraph gPL;
        gPL.SetName("gPL");
        gPL.Set(NumMomBins);

        for(unsigned uBin=0; uBin<NumMomBins; uBin++){
            gPP.SetPoint(uBin,Kitty_pp.GetMomentum(uBin),Kitty_pp.GetCorrFun(uBin));
            gPL.SetPoint(uBin,Kitty_pL.GetMomentum(uBin),Kitty_pL.GetCorrFun(uBin));
        }
        gPP.Write();
        gPL.Write();

        hPP_Source->Write();
        fPP_Gauss->Write();

        hPL_Source->Write();
        fPL_Gauss->Write();

        delete hPP_Source;
        delete fPP_Gauss;
        delete hPL_Source;
        delete fPL_Gauss;
        delete OutputFile;
    }

    fCommonOut->cd();
    hPP_RAD->Write();
    hPL_RAD->Write();

    delete hPP_RAD;
    delete hPL_RAD;
    delete fCommonOut;
}

//TYPE 0 only Gauss
//TYPE 1 is back-to-back
//TYPE 2 is with the proper angular disto
void TestEposDistos(const double& SourceSize){

    const double k_CutOff = 200;
    const TString OutputDir = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/SourceStudies/TestEposDistos/";

    TFile* fCommonOut = new TFile(
    TString::Format("%sfCommonOut_%.2f.root",OutputDir.Data(),SourceSize),
    "recreate");
    const unsigned NumTypes = 3;
    TH1F* hPP_RAD = new TH1F("hPP_RAD","hPP_RAD",NumTypes,-0.5,5.5);
    TH1F* hPL_RAD = new TH1F("hPL_RAD","hPL_RAD",NumTypes,-0.5,5.5);

    for(int TYPE=0; TYPE<NumTypes; TYPE++){
        const unsigned NumMomBins = 75;
        const double kMin = 0;
        const double kMax = 300;

        const double Mass_pi0 = 134.9766;
        const double Mass_pic = 139.57018;
        const double Mass_p = 938.272;
        const double Mass_L = 1115.683;

        double ppPotPars1S0[8]={NN_AV18,v18_Coupled3P2,1,1,1,0,0,0};
        double ppPotPars3P0[8]={NN_AV18,v18_Coupled3P2,1,1,1,1,1,0};
        double ppPotPars3P1[8]={NN_AV18,v18_Coupled3P2,1,1,1,1,1,1};
        double ppPotPars3P2[8]={NN_AV18,v18_Coupled3P2,1,1,1,1,1,2};
        CATSparameters cPotPars1S0_pp(CATSparameters::tPotential,8,true);cPotPars1S0_pp.SetParameters(ppPotPars1S0);
        CATSparameters cPotPars3P0_pp(CATSparameters::tPotential,8,true);cPotPars3P0_pp.SetParameters(ppPotPars3P0);
        CATSparameters cPotPars3P1_pp(CATSparameters::tPotential,8,true);cPotPars3P1_pp.SetParameters(ppPotPars3P1);
        CATSparameters cPotPars3P2_pp(CATSparameters::tPotential,8,true);cPotPars3P2_pp.SetParameters(ppPotPars3P2);
        //CATSparameters cSorPars_pp(CATSparameters::tSource,2,true);cSorPars_pp.SetParameter(0,1);cSorPars_pp.SetParameter(1,2);

        double pL_PotPars1S0[8]={pL_UsmaniOli,0,0,0,0,0,0,0};
        double pL_PotPars3S1[8]={pL_UsmaniOli,0,0,0,0,1,0,1};
        //CATSparameters cSorPars_pL(CATSparameters::tSource,2,true);cSorPars_pL.SetParameter(0,1);cSorPars_pL.SetParameter(1,2);
        CATSparameters cPotPars1S0_pL(CATSparameters::tPotential,8,true);cPotPars1S0_pL.SetParameters(pL_PotPars1S0);
        CATSparameters cPotPars3S1_pL(CATSparameters::tPotential,8,true);cPotPars3S1_pL.SetParameters(pL_PotPars3S1);


        DLM_CleverMcLevyResoTM CleverSource_pp;
        CleverSource_pp.InitStability(1,2-1e-6,2+1e-6);
        CleverSource_pp.InitScale(38,0.15,2.0);
        CleverSource_pp.InitRad(257,0,64);
        CleverSource_pp.InitType(2);

        DLM_CleverMcLevyResoTM CleverSource_pL;
        CleverSource_pL.InitStability(1,2-1e-6,2+1e-6);
        CleverSource_pL.InitScale(38,0.15,2.0);
        CleverSource_pL.InitRad(257,0,64);
        CleverSource_pL.InitType(2);

        DLM_Histo<double>* HISTO_P_AVG;
        DLM_Histo<double>* HISTO_L_AVG;

        DLM_Random RanGen(11);
        double RanVal1;
        double RanVal2;
        double RanVal3;

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

        if(TYPE==0){
            CleverSource_pp.SetUpReso(0,0.6422);
            CleverSource_pp.SetUpReso(1,0.6422);
            CleverSource_pL.SetUpReso(0,0.6438);
            CleverSource_pL.SetUpReso(1,0.6438);
        }
        else if(TYPE==1){
            CleverSource_pp.SetUpReso(0,0.6422);
            CleverSource_pp.SetUpReso(1,0.6422);
            CleverSource_pL.SetUpReso(0,0.6438);
            CleverSource_pL.SetUpReso(1,0.6438);
            const unsigned NumIter = 100;
            for(unsigned uEntry=0; uEntry<NumIter; uEntry++){
                CleverSource_pp.AddBGT_PR(490./1362.*1.65,1.);
                CleverSource_pp.AddBGT_RP(490./1362.*1.65,-1.);
                CleverSource_pp.AddBGT_RR(490./1362.*1.65,-1.,490./1362.*1.65,1.,-1.);

                CleverSource_pL.AddBGT_PR(360./1462.*4.69,1.);
                CleverSource_pL.AddBGT_RP(490./1362.*1.65,-1.);
                CleverSource_pL.AddBGT_RR(490./1362.*1.65,-1.,360./1462.*4.69,1.,-1.);
                //CleverSource_pp.AddBGT_RP(490,1);
            }
        }
        else if(TYPE==2){
            CleverSource_pp.SetUpReso(0,0.6422);
            CleverSource_pp.SetUpReso(1,0.6422);

            TFile* F_EposDisto_p_pReso = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_3/EposDisto_p_pReso.root");
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
                if(k_D>k_CutOff) continue;
                RanVal1 = RanGen.Exponential(fM2/(fP2*Tau2));
                CleverSource_pp.AddBGT_PR(RanVal1,-cos(AngleRcP2));
                CleverSource_pp.AddBGT_RP(RanVal1,cos(AngleRcP2));
            }
            delete F_EposDisto_p_pReso;

            TFile* F_EposDisto_pReso_pReso = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_3/EposDisto_pReso_pReso.root");
            TNtuple* T_EposDisto_pReso_pReso = (TNtuple*)F_EposDisto_pReso_pReso->Get("InfoTuple_ClosePairs");
            unsigned N_EposDisto_pReso_pReso = T_EposDisto_pReso_pReso->GetEntries();
            T_EposDisto_pReso_pReso->SetBranchAddress("k_D",&k_D);
            T_EposDisto_pReso_pReso->SetBranchAddress("P1",&fP1);
            T_EposDisto_pReso_pReso->SetBranchAddress("P2",&fP2);
            T_EposDisto_pReso_pReso->SetBranchAddress("M1",&fM1);
            T_EposDisto_pReso_pReso->SetBranchAddress("M2",&fM2);
            T_EposDisto_pReso_pReso->SetBranchAddress("Tau1",&Tau1);
            T_EposDisto_pReso_pReso->SetBranchAddress("Tau2",&Tau2);
            T_EposDisto_pReso_pReso->SetBranchAddress("AngleRcP1",&AngleRcP1);
            T_EposDisto_pReso_pReso->SetBranchAddress("AngleRcP2",&AngleRcP2);
            T_EposDisto_pReso_pReso->SetBranchAddress("AngleP1P2",&AngleP1P2);
            for(unsigned uEntry=0; uEntry<N_EposDisto_pReso_pReso; uEntry++){
                T_EposDisto_pReso_pReso->GetEntry(uEntry);
Tau1 = 1.65;
Tau2 = 1.65;
                if(k_D>k_CutOff) continue;
                RanVal1 = RanGen.Exponential(fM1/(fP1*Tau1));
                RanVal2 = RanGen.Exponential(fM2/(fP2*Tau2));
                CleverSource_pp.AddBGT_RR(RanVal1,cos(AngleRcP1),RanVal2,cos(AngleRcP2),cos(AngleP1P2));
            }
            delete F_EposDisto_pReso_pReso;

            CleverSource_pL.SetUpReso(0,0.6438);
            CleverSource_pL.SetUpReso(1,0.6438);

            TFile* F_EposDisto_p_LamReso = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_3/EposDisto_p_LamReso.root");
            TNtuple* T_EposDisto_p_LamReso = (TNtuple*)F_EposDisto_p_LamReso->Get("InfoTuple_ClosePairs");
            unsigned N_EposDisto_p_LamReso = T_EposDisto_p_LamReso->GetEntries();
            T_EposDisto_p_LamReso->SetBranchAddress("k_D",&k_D);
            T_EposDisto_p_LamReso->SetBranchAddress("P1",&fP1);
            T_EposDisto_p_LamReso->SetBranchAddress("P2",&fP2);
            T_EposDisto_p_LamReso->SetBranchAddress("M1",&fM1);
            T_EposDisto_p_LamReso->SetBranchAddress("M2",&fM2);
            T_EposDisto_p_LamReso->SetBranchAddress("Tau1",&Tau1);
            T_EposDisto_p_LamReso->SetBranchAddress("Tau2",&Tau2);
            T_EposDisto_p_LamReso->SetBranchAddress("AngleRcP1",&AngleRcP1);
            T_EposDisto_p_LamReso->SetBranchAddress("AngleRcP2",&AngleRcP2);
            T_EposDisto_p_LamReso->SetBranchAddress("AngleP1P2",&AngleP1P2);
            for(unsigned uEntry=0; uEntry<N_EposDisto_p_LamReso; uEntry++){
                T_EposDisto_p_LamReso->GetEntry(uEntry);
Tau1 = 0;
Tau2 = 4.69;
                if(k_D>k_CutOff) continue;
                RanVal2 = RanGen.Exponential(fM2/(fP2*Tau2));
                CleverSource_pL.AddBGT_PR(RanVal2,cos(AngleRcP2));
            }
            delete F_EposDisto_p_LamReso;

            TFile* F_EposDisto_pReso_Lam = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_3/EposDisto_pReso_Lam.root");
            TNtuple* T_EposDisto_pReso_Lam = (TNtuple*)F_EposDisto_pReso_Lam->Get("InfoTuple_ClosePairs");
            unsigned N_EposDisto_pReso_Lam = T_EposDisto_pReso_Lam->GetEntries();
            T_EposDisto_pReso_Lam->SetBranchAddress("k_D",&k_D);
            T_EposDisto_pReso_Lam->SetBranchAddress("P1",&fP1);
            T_EposDisto_pReso_Lam->SetBranchAddress("P2",&fP2);
            T_EposDisto_pReso_Lam->SetBranchAddress("M1",&fM1);
            T_EposDisto_pReso_Lam->SetBranchAddress("M2",&fM2);
            T_EposDisto_pReso_Lam->SetBranchAddress("Tau1",&Tau1);
            T_EposDisto_pReso_Lam->SetBranchAddress("Tau2",&Tau2);
            T_EposDisto_pReso_Lam->SetBranchAddress("AngleRcP1",&AngleRcP1);
            T_EposDisto_pReso_Lam->SetBranchAddress("AngleRcP2",&AngleRcP2);
            T_EposDisto_pReso_Lam->SetBranchAddress("AngleP1P2",&AngleP1P2);
            for(unsigned uEntry=0; uEntry<N_EposDisto_pReso_Lam; uEntry++){
                T_EposDisto_pReso_Lam->GetEntry(uEntry);
Tau1 = 1.65;
Tau2 = 0;
                if(k_D>k_CutOff) continue;
                RanVal1 = RanGen.Exponential(fM1/(fP1*Tau1));
                CleverSource_pL.AddBGT_RP(RanVal1,cos(AngleRcP1));
            }
            delete F_EposDisto_pReso_Lam;

            TFile* F_EposDisto_pReso_LamReso = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_3/EposDisto_pReso_LamReso.root");
            TNtuple* T_EposDisto_pReso_LamReso = (TNtuple*)F_EposDisto_pReso_LamReso->Get("InfoTuple_ClosePairs");
            unsigned N_EposDisto_pReso_LamReso = T_EposDisto_pReso_LamReso->GetEntries();
            T_EposDisto_pReso_LamReso->SetBranchAddress("k_D",&k_D);
            T_EposDisto_pReso_LamReso->SetBranchAddress("P1",&fP1);
            T_EposDisto_pReso_LamReso->SetBranchAddress("P2",&fP2);
            T_EposDisto_pReso_LamReso->SetBranchAddress("M1",&fM1);
            T_EposDisto_pReso_LamReso->SetBranchAddress("M2",&fM2);
            T_EposDisto_pReso_LamReso->SetBranchAddress("Tau1",&Tau1);
            T_EposDisto_pReso_LamReso->SetBranchAddress("Tau2",&Tau2);
            T_EposDisto_pReso_LamReso->SetBranchAddress("AngleRcP1",&AngleRcP1);
            T_EposDisto_pReso_LamReso->SetBranchAddress("AngleRcP2",&AngleRcP2);
            T_EposDisto_pReso_LamReso->SetBranchAddress("AngleP1P2",&AngleP1P2);
            for(unsigned uEntry=0; uEntry<N_EposDisto_pReso_LamReso; uEntry++){
                T_EposDisto_pReso_LamReso->GetEntry(uEntry);
Tau1 = 1.65;
Tau2 = 4.69;
                if(k_D>k_CutOff) continue;
                RanVal1 = RanGen.Exponential(fM1/(fP1*Tau1));
                RanVal2 = RanGen.Exponential(fM2/(fP2*Tau2));
                CleverSource_pL.AddBGT_RR(RanVal1,cos(AngleRcP1),RanVal2,cos(AngleRcP2),cos(AngleP1P2));
            }
            delete F_EposDisto_pReso_LamReso;
        }
        else{
            abort();
        }

        CleverSource_pp.InitNumMcIter(262144);
        CleverSource_pL.InitNumMcIter(262144);

        CATS Kitty_pp;
        Kitty_pp.SetMomBins(NumMomBins,kMin,kMax);
        Kitty_pp.SetAnaSource(CatsSourceForwarder, &CleverSource_pp, 2);
        Kitty_pp.SetAnaSource(0,SourceSize);
        Kitty_pp.SetAnaSource(1,2.0);
        Kitty_pp.SetUseAnalyticSource(true);
        Kitty_pp.SetMomentumDependentSource(false);
        Kitty_pp.SetThetaDependentSource(false);
        Kitty_pp.SetExcludeFailedBins(false);
        Kitty_pp.SetQ1Q2(1);
        Kitty_pp.SetPdgId(2212, 2212);
        Kitty_pp.SetRedMass( 0.5*Mass_p );
        Kitty_pp.SetNumChannels(4);
        Kitty_pp.SetNumPW(0,1);
        Kitty_pp.SetNumPW(1,2);
        Kitty_pp.SetNumPW(2,2);
        Kitty_pp.SetNumPW(3,2);
        Kitty_pp.SetSpin(0,0);
        Kitty_pp.SetSpin(1,1);
        Kitty_pp.SetSpin(2,1);
        Kitty_pp.SetSpin(3,1);
        Kitty_pp.SetChannelWeight(0, 3./12.);
        Kitty_pp.SetChannelWeight(1, 1./12.);
        Kitty_pp.SetChannelWeight(2, 3./12.);
        Kitty_pp.SetChannelWeight(3, 5./12.);
        Kitty_pp.SetShortRangePotential(0,0,fDlmPot,cPotPars1S0_pp);
        Kitty_pp.SetShortRangePotential(1,1,fDlmPot,cPotPars3P0_pp);
        Kitty_pp.SetShortRangePotential(2,1,fDlmPot,cPotPars3P1_pp);
        Kitty_pp.SetShortRangePotential(3,1,fDlmPot,cPotPars3P2_pp);
        //Kitty_pp.SetEpsilonConv(5e-8);
        //Kitty_pp.SetEpsilonProp(5e-8);
        Kitty_pp.SetMaxNumThreads(4);
        Kitty_pp.KillTheCat();

        CATS Kitty_pL;
        Kitty_pL.SetMomBins(NumMomBins,kMin,kMax);
        Kitty_pL.SetAnaSource(CatsSourceForwarder, &CleverSource_pL, 2);
        Kitty_pL.SetAnaSource(0,SourceSize);
        Kitty_pL.SetAnaSource(1,2.0);
        Kitty_pL.SetUseAnalyticSource(true);
        Kitty_pL.SetMomentumDependentSource(false);
        Kitty_pL.SetThetaDependentSource(false);
        Kitty_pL.SetExcludeFailedBins(false);
        Kitty_pL.SetQ1Q2(0);
        Kitty_pL.SetPdgId(2212, 3122);
        Kitty_pL.SetRedMass( (Mass_p*Mass_L)/(Mass_p+Mass_L) );
        Kitty_pL.SetNumChannels(2);
        Kitty_pL.SetNumPW(0,1);
        Kitty_pL.SetNumPW(1,1);
        Kitty_pL.SetSpin(0,0);
        Kitty_pL.SetSpin(1,1);
        Kitty_pL.SetChannelWeight(0, 1./4.);
        Kitty_pL.SetChannelWeight(1, 3./4.);
        Kitty_pL.SetShortRangePotential(0,0,fDlmPot,cPotPars1S0_pL);
        Kitty_pL.SetShortRangePotential(1,0,fDlmPot,cPotPars3S1_pL);
        //Kitty_pL.SetEpsilonConv(5e-8);
        //Kitty_pL.SetEpsilonProp(5e-8);
        Kitty_pL.SetMaxNumThreads(4);
        Kitty_pL.KillTheCat();

        TFile* OutputFile = new TFile(
        TString::Format("%sTDAD_%.2f_%i.root",OutputDir.Data(),SourceSize,TYPE)
                         ,"recreate");

        TH1F* hPP_Source;
        TF1* fPP_Gauss;
        TH1F* hPL_Source;
        TF1* fPL_Gauss;

        const double rMin = 0;
        const double rMax = 8;
        const unsigned rNumBins = 128;
        const double RELERR = 0.002;
        const double ABSERR = 0.002;
        double RAD;
        double SVAL;
        hPP_Source = new TH1F("hPP_Source","hPP_Source",rNumBins,rMin,rMax);
        hPL_Source = new TH1F("hPL_Source","hPL_Source",rNumBins,rMin,rMax);
        for(unsigned uRad=0; uRad<rNumBins; uRad++){
            RAD = hPP_Source->GetBinCenter(uRad+1);
            SVAL = Kitty_pp.EvaluateTheSource(10,RAD,0);
            hPP_Source->SetBinContent(uRad+1,SVAL);
            hPP_Source->SetBinError(uRad+1,RELERR*SVAL+ABSERR);

            RAD = hPL_Source->GetBinCenter(uRad+1);
            SVAL = Kitty_pL.EvaluateTheSource(10,RAD,0);
            hPL_Source->SetBinContent(uRad+1,SVAL);
            hPL_Source->SetBinError(uRad+1,RELERR*SVAL+ABSERR);
        }

        fPP_Gauss = new TF1(TString::Format("fPP_Gauss"),GaussSourceTF1,rMin,rMax,1);
        fPP_Gauss->SetParameter(0,1.0);
        fPP_Gauss->SetParLimits(0,0.5,2.0);

        fPL_Gauss = new TF1(TString::Format("fPL_Gauss"),GaussSourceTF1,rMin,rMax,1);
        fPL_Gauss->SetParameter(0,1.0);
        fPL_Gauss->SetParLimits(0,0.5,2.0);

        hPP_Source->Fit(fPP_Gauss,"S, N, R, M");
        hPL_Source->Fit(fPL_Gauss,"S, N, R, M");

        hPP_RAD->SetBinContent(TYPE+1,fPP_Gauss->GetParameter(0));
        hPP_RAD->SetBinError(TYPE+1,fPP_Gauss->GetParError(0));

        hPL_RAD->SetBinContent(TYPE+1,fPL_Gauss->GetParameter(0));
        hPL_RAD->SetBinError(TYPE+1,fPL_Gauss->GetParError(0));



        TGraph gPP;
        gPP.SetName("gPP");
        gPP.Set(NumMomBins);

        TGraph gPL;
        gPL.SetName("gPL");
        gPL.Set(NumMomBins);

        for(unsigned uBin=0; uBin<NumMomBins; uBin++){
            gPP.SetPoint(uBin,Kitty_pp.GetMomentum(uBin),Kitty_pp.GetCorrFun(uBin));
            gPL.SetPoint(uBin,Kitty_pL.GetMomentum(uBin),Kitty_pL.GetCorrFun(uBin));
        }
        gPP.Write();
        gPL.Write();

        hPP_Source->Write();
        fPP_Gauss->Write();

        hPL_Source->Write();
        fPL_Gauss->Write();

        delete hPP_Source;
        delete fPP_Gauss;
        delete hPL_Source;
        delete fPL_Gauss;
        delete OutputFile;
    }

    fCommonOut->cd();
    hPP_RAD->Write();
    hPL_RAD->Write();

    delete hPP_RAD;
    delete hPL_RAD;
    delete fCommonOut;

}

//for the McLevyNolan_Reso
//void CreateSourceLookUp(const TString& System, const ){

//}
void Compare_2_vs_3_body(){
  TString InputFolder_2body = "/mnt/Ubuntu_Data/CernBox/Sync/Plots/SourcePaper/PLB_Review2/2body/";
  TString InputFolder_3body = "/mnt/Ubuntu_Data/CernBox/Sync/Plots/SourcePaper/PLB_Review2/3body/";
  TString OutputFolder = "/mnt/Ubuntu_Data/CernBox/Sync/Plots/SourcePaper/PLB_Review2a/";

  TGraph* gS2;
  TGraph* gS3;
  TGraph* gC2;
  TGraph* gC3;
  TGraph* gSG2;
  TGraph* gCORE;

  TGraph* gSG2spike = new TGraph();
  gSG2spike->SetName("gSG2spike");

  TGraph* gCOREspike = new TGraph();
  gCOREspike->SetName("gCOREspike");

  TGraph* gS2spike = new TGraph();
  gS2spike->SetName("gS2spike");
  TGraph* gS3spike = new TGraph();
  gS3spike->SetName("gS3spike");

  TGraph* gSratio = new TGraph();
  gSratio->SetName("gSratio");
  TGraph* gCratio = new TGraph();
  gCratio->SetName("gCratio");



  TFile* file2 = new TFile(InputFolder_2body+"fOut.root","read");
  gS2 = (TGraph*)file2->Get("gSr_pp_CoreReso");
  gSG2 = (TGraph*)file2->Get("gSr_pp_Gauss");
  gCORE = (TGraph*)file2->Get("gSr_pp_Core");
  gC2 = (TGraph*)file2->Get("gCk_pp_CoreReso");

  TFile* file3 = new TFile(InputFolder_3body+"fOut.root","read");
  gS3 = (TGraph*)file3->Get("gSr_pp_CoreReso");
  gC3 = (TGraph*)file3->Get("gCk_pp_CoreReso");

  double xval3, yval3;
  double xval2, yval2;
  double xval2G, yval2G;
  double xvalC, yvalC;
  for(unsigned uBin=0; uBin<gS2->GetN(); uBin++){
    gS2->GetPoint(uBin,xval2,yval2);
    gS3->GetPoint(uBin,xval3,yval3);
    gSG2->GetPoint(uBin,xval2G,yval2G);
    gCORE->GetPoint(uBin,xvalC,yvalC);
    gS2spike->SetPoint(uBin,xval2,yval2/(4.*TMath::Pi()*xval2*xval2));
    gSG2spike->SetPoint(uBin,xval2G,yval2G/(4.*TMath::Pi()*xval2G*xval2G));
    gCOREspike->SetPoint(uBin,xvalC,yvalC/(4.*TMath::Pi()*xvalC*xvalC));
    gS3spike->SetPoint(uBin,xval3,yval3/(4.*TMath::Pi()*xval3*xval3));
    gSratio->SetPoint(uBin,xval2,yval3/yval2);
  }
  for(unsigned uBin=0; uBin<gC2->GetN(); uBin++){
    gC2->GetPoint(uBin,xval2,yval2);
    gC3->GetPoint(uBin,xval3,yval3);
    gCratio->SetPoint(uBin,xval2,yval3/yval2);
  }

  TH1F* hC_axis = new TH1F("hC_axis","hC_axis",128,0,200);
  hC_axis->SetStats(false);
  hC_axis->SetTitle("; #it{k*} (MeV/#it{c}); #it{C}(#it{k*})");
  hC_axis->GetXaxis()->SetRangeUser(0, 200);
  hC_axis->GetXaxis()->SetNdivisions(505);
  hC_axis->GetYaxis()->SetRangeUser(-0.1, 4.5);
  for(unsigned uBin=0; uBin<128; uBin++) hC_axis->SetBinContent(uBin+1,1);
  SetStyleHisto2a(hC_axis,2,0);

  DLM_SubPads* DlmPad_Ck = new DLM_SubPads(720,1080);
  DlmPad_Ck->AddSubPad(0,1,0.5,1);
  DlmPad_Ck->AddSubPad(0,1,0.25,0.5);
  DlmPad_Ck->AddSubPad(0,1,0.,0.25);
  DlmPad_Ck->SetMargin(0,0.14,0.02,0.0,0.02);
  DlmPad_Ck->SetMargin(1,0.14,0.02,0.0,0.0);
  DlmPad_Ck->SetMargin(2,0.14,0.02,0.07,0.0);
  DlmPad_Ck->cd(0);
  hC_axis->Draw();
  gC2->SetLineWidth(6);
  gC2->SetLineStyle(1);
  gC2->SetLineColor(kBlue);
  gC2->Draw("same,C");
  gC3->SetLineWidth(4);
  gC3->SetLineStyle(1);
  gC3->SetLineColor(kRed);
  gC3->Draw("same,C");

  //INLET -------------------------------------------------------------------------------------------------------------------

  TH1F* hC_Inlet = (TH1F*)hC_axis->Clone("hC_Inlet");
  hC_Inlet->SetStats(false);
  hC_Inlet->SetMarkerSize(hC_axis->GetMarkerSize()*0.67);
  hC_Inlet->SetLineWidth(hC_axis->GetLineWidth()*0.67);
  hC_Inlet->GetXaxis()->SetRangeUser(0, 200);
  SetStyleHisto2a(hC_Inlet,2,0,2);
  hC_Inlet->GetXaxis()->SetNdivisions(505);
  hC_Inlet->GetYaxis()->SetNdivisions(205);
  hC_Inlet->GetYaxis()->SetRangeUser(0.91, 1.00);

  DlmPad_Ck->cd(1);
  hC_Inlet->Draw();
  gC2->Draw("same,C");
  gC3->Draw("same,C");


  TH1F* hC_Rat = new TH1F("hC_Rat", "hC_Rat", 128, 0, 200);
  hC_Rat->SetStats(false);
  hC_Rat->GetXaxis()->SetRangeUser(0, 200);
  hC_Rat->SetTitle("; #it{k*} (MeV/#it{c}); Ratio");
  hC_Rat->GetYaxis()->SetRangeUser(0.985,1.015);
  hC_Rat->GetXaxis()->SetNdivisions(505);
  for(unsigned uBin=0; uBin<128; uBin++) hC_Rat->SetBinContent(uBin+1,1);
  SetStyleHisto2a(hC_Rat,2,0,2);
  hC_Rat->GetYaxis()->SetNdivisions(504);

  DlmPad_Ck->cd(2);
  hC_Rat->Draw("");
  gCratio->SetLineColor(kBlue+2);
  gCratio->SetLineWidth(6);
  gCratio->Draw("same,C");

  DlmPad_Ck->GetCanvas()->SaveAs(OutputFolder+"DlmPad_Ck.pdf");




  TH1F* hS_axis = new TH1F("hS_axis","hS_axis",128,0,10);
  hS_axis->SetStats(false);
  hS_axis->SetTitle("; #it{r*} (fm); 4#pir*^{2}#it{S(r*)} (1/fm)");
  hS_axis->GetXaxis()->SetRangeUser(0, 10);
  hS_axis->GetXaxis()->SetNdivisions(505);
  hS_axis->GetYaxis()->SetRangeUser(0.001, 1);
  for(unsigned uBin=0; uBin<128; uBin++) hS_axis->SetBinContent(uBin+1,1);
  SetStyleHisto2a(hS_axis,2,0);

  DLM_SubPads* DlmPad_Sr = new DLM_SubPads(720,1080);
  DlmPad_Sr->AddSubPad(0,1,0.5,1);
  DlmPad_Sr->AddSubPad(0,1,0.25,0.5);
  DlmPad_Sr->AddSubPad(0,1,0.,0.25);
  DlmPad_Sr->SetMargin(0,0.14,0.02,0.0,0.02);
  DlmPad_Sr->SetMargin(1,0.14,0.02,0.0,0.0);
  DlmPad_Sr->SetMargin(2,0.14,0.02,0.07,0.0);
  DlmPad_Sr->cd(0);
  DlmPad_Sr->SetLogy(0);
  hS_axis->Draw();

  gSG2->SetLineWidth(3);
  gSG2->SetLineColor(kBlack);
  gSG2->SetLineStyle(1);
  gSG2->Draw("same,C");
  gCORE->SetLineWidth(6);
  gCORE->SetLineStyle(2);
  gCORE->SetLineColor(kBlack);
  gCORE->Draw("same,C");
  gS2->SetLineWidth(6);
  gS2->SetLineStyle(1);
  gS2->SetLineColor(kBlue);
  gS2->Draw("same,C");
  gS3->SetLineWidth(4);
  gS3->SetLineStyle(1);
  gS3->SetLineColor(kRed);
  gS3->Draw("same,C");

  DLM_SubPads* DlmPad_Sr2 = new DLM_SubPads(720,720);
  DlmPad_Sr2->AddSubPad(0,1,0.35,1);
  DlmPad_Sr2->AddSubPad(0,1,0.,0.35);
  DlmPad_Sr2->SetMargin(0,0.14,0.02,0.0,0.02);
  DlmPad_Sr2->SetMargin(1,0.14,0.02,0.09,0.0);
  DlmPad_Sr2->cd(0);
  DlmPad_Sr2->SetLogy(0);
  hS_axis->Draw();
  gSG2->Draw("same,C");
  gCORE->Draw("same,C");
  gS2->Draw("same,C");
  gS3->Draw("same,C");

  //INLET -------------------------------------------------------------------------------------------------------------------

  TH1F* hS_Inlet = (TH1F*)hS_axis->Clone("hS_Inlet");
  hS_Inlet->SetStats(false);
  hS_Inlet->SetTitle("; #it{r*} (fm); #it{S(r*)} (1/fm)^{3}");
  hS_Inlet->SetMarkerSize(hS_axis->GetMarkerSize()*0.67);
  hS_Inlet->SetLineWidth(hS_axis->GetLineWidth()*0.67);
  hS_Inlet->GetXaxis()->SetRangeUser(0, 10);
  SetStyleHisto2a(hS_Inlet,2,0,2);
  hS_Inlet->GetXaxis()->SetNdivisions(505);
  hS_Inlet->GetYaxis()->SetNdivisions(205);
  hS_Inlet->GetYaxis()->SetRangeUser(5e-7, 0.2);

  DlmPad_Sr->cd(1);
  DlmPad_Sr->SetLogy(1);
  hS_Inlet->Draw();
  gSG2spike->SetLineWidth(6);
  gSG2spike->SetLineStyle(1);
  gSG2spike->SetLineColor(kBlack);
  gSG2spike->Draw("same,C");
  gCOREspike->SetLineWidth(6);
  gCOREspike->SetLineStyle(2);
  gCOREspike->SetLineColor(kBlack);
  gCOREspike->Draw("same,C");
  gS2spike->SetLineWidth(6);
  gS2spike->SetLineColor(kBlue);
  gS2spike->Draw("same,C");
  gS3spike->SetLineWidth(4);
  gS3spike->SetLineColor(kRed);
  gS3spike->Draw("same,C");


  TH1F* hS_Rat = new TH1F("hS_Rat", "hS_Rat", 128, 0, 10);
  hS_Rat->SetStats(false);
  hS_Rat->GetXaxis()->SetRangeUser(0, 10);
  hS_Rat->SetTitle("; #it{r*} (fm); Ratio");
  hS_Rat->GetYaxis()->SetRangeUser(0.945,1.055);
  hS_Rat->GetXaxis()->SetNdivisions(505);
  for(unsigned uBin=0; uBin<128; uBin++) hS_Rat->SetBinContent(uBin+1,0);
  SetStyleHisto2a(hS_Rat,2,0,2);
  hS_Rat->GetYaxis()->SetNdivisions(504);

  DlmPad_Sr->cd(2);
  hS_Rat->Draw("");
  gSratio->SetLineColor(kBlue+2);
  gSratio->SetLineWidth(6);
  gSratio->Draw("same,C");
  DlmPad_Sr->GetCanvas()->SaveAs(OutputFolder+"DlmPad_Sr.pdf");

  DlmPad_Sr2->cd(1);
  hS_Rat->Draw("");
  gSratio->SetLineColor(kBlue+2);
  gSratio->SetLineWidth(6);
  gSratio->Draw("same,C");
  DlmPad_Sr2->GetCanvas()->SaveAs(OutputFolder+"DlmPad_Sr2.pdf");



  TFile* fOut = new TFile(OutputFolder+"Compare_2_vs_3_body.root","recreate");
  gS2->Write();
  gS3->Write();
  gSratio->Write();
  gC2->Write();
  gC3->Write();
  gCratio->Write();


  delete hC_axis;
  delete gCratio;
  delete gSratio;
  delete file2;
  delete file3;
  delete fOut;
}


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
double Baseline_BernieSource(double* x, double* par){
    double& k = *x;
    double& p0 = par[0];
    //constrained polynomials
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
DLM_CkDecomposition* SOURCE_FIT=NULL;
//par[0] is the radius
//we keep par[1] free just in case we end up doing a Levy
//from par[2]... the baseline
double Femto_BernieSource(double* x, double* par){
  if(!SOURCE_FIT) return 0;
  SOURCE_FIT->GetCk()->SetSourcePar(0,par[0]);
  SOURCE_FIT->Update(false,false);
  //SOURCE_FIT->Update(true,true);
  //printf(" DECOMP: pL: %.3f; pS0: %.3f; pXim: %.3f; pXi0: %.3f\n",
  //  SOURCE_FIT->GetContribution("pLambda")->EvalCk(*x),
  //  SOURCE_FIT->GetContribution("pSigma0")->EvalCk(*x),
  //  SOURCE_FIT->GetContribution("pXim")->EvalCk(*x),
  //  SOURCE_FIT->GetContribution("pXi0")->EvalCk(*x));
  //printf(" DLMCK:  pL: %.3f; pS0: %.3f; pXim: %.3f; pXi0: %.3f\n",
  //  SOURCE_FIT->GetContribution("pLambda")->GetCk()->Eval(*x),
  //  SOURCE_FIT->GetContribution("pSigma0")->GetCk()->Eval(*x),
  //  SOURCE_FIT->GetContribution("pXim")->GetCk()->Eval(*x),
  //  SOURCE_FIT->GetContribution("pXi0")->GetCk()->Eval(*x));
  return SOURCE_FIT->EvalCk(*x);
  //return SOURCE_FIT->GetCk()->Eval(*x);
}
double Fit_BernieSource(double* x, double* par){
  double Femto = Femto_BernieSource(x,par);
  double Baseline = Baseline_BernieSource(x,&par[2]);
  //printf("k* = %.0f\n",*x);
  //printf(" par = %.3e %.3e %.3e %.3e %.3e %.3e %.3e\n",par[0],par[1],par[2],par[3],par[4],par[5],par[6]);
  //printf(" Femto = %.3f; Baseline = %.3f\n",Femto,Baseline);
  //usleep(100e3);
  return Femto*Baseline;
}

//SmearStrategy == 1 is the new one
void SourcePaper_pp(const TString SourceVar, const int& SmearStrategy, const unsigned DataVar, const int imTbin, const TString OutputFolder){
//printf("0\n");
  gROOT->ProcessLine("gErrorIgnoreLevel = 2001;");
  const bool Silent = true;
  //for the pp source, the error on the sign
  const bool CorrectSign = true;
  //std::vector<TString> SourceVar = {"Gauss","McGauss_ResoTM"};//done
//std::vector<TString> SourceVar = {"McLevy_ResoTM"};
  std::vector<int> pp_lam_var = {0,1,2};//done
//std::vector<int> pp_lam_var = {0};
  std::vector<int> pL_lam_var = {0};
  //0 is the old one, 1 is the new one with the folded ME
  //std::vector<int> SmearStrategy = {0,1};//done
  //std::vector<float> FemtoRegion = {376,352,400};//done
std::vector<float> FemtoRegion = {376};
  const unsigned NumMomBins_pp = TMath::Nint(FemtoRegion.back()/4.);
  const unsigned NumMomBins_feed = TMath::Nint(FemtoRegion.back()/10.);
  std::vector<float> BaselineRegion = {500};
  std::vector<float> CkCutOff = {700};//done
  std::vector<int> pL_pot_var = {11600,-11600};//done
//std::vector<int> pL_pot_var = {11600};
  //pp13TeV_HM_DimiJun20 -> the ME reso matrix
  //pp13TeV_HM_BernieSource -> Bernie's old smearing matrix
  //pp13TeV_HM_BernieSource -> contains the data
  std::vector<TString> MomSmearVar = {"pp13TeV_HM_DimiJun20"};//done
  //and an additional variation will be done on if pS0 is included as a feed or not
  enum BLTYPE { pol0s,pol1s,pol2s,pol3s,dpol2s,dpol3s,dpol4s,pol2e,pol3e,dpol2e,dpol3e,dpol4e,spl1 };
  std::vector<int> BaselineVar = {pol0s,pol1s,dpol3e};
//std::vector<int> BaselineVar = {pol0s,pol1s};
  const bool pS0_Var = true;

  //the region for which the DLM_Ck objects will be defined
  //we put some extra to get away from the edge effects of the smearing
  const double kMin = 0;
  const double kMax = 600;
  const unsigned TotMomBins = TMath::Nint(kMax);
  const unsigned NumMtBins = 7;
  if(imTbin>=NumMtBins){
    printf("\033[1;31mERROR:\033[0m Only %u mT bins are available!\n",NumMtBins);
    return;
  }
  //approximate, based on the output we kind of know we will get
  //std::vector<float> ExpectedRadii = {1.35,1.3,1.25,1.2,1.1,1.05,0.95};
  //const double Scale_pL = 1.1;
  //const double Scale_pS0 = 1.15;
  //const double Scale_pXi = 0.97;
  //const double Scale_core = 0.94;
  std::vector<float> ExpectedRadii = { 1.55, 1.473, 1.421, 1.368, 1.295, 1.220,1.124 };
  //std::vector<float> ExpectedRadii = {1.35,1.3,1.25,1.2,1.1,1.05,0.95};
  const double Scale_pL = 1.0;
  const double Scale_pS0 = 1.0;
  const double Scale_pXi = 1.0;
  const double Scale_core = 1.0;
  //the difference in the effectiv Gaussian compered to pp
//printf("1\n");
  const unsigned NumSourcePars = 1;
  //std::vector<float> pSigma0Radii = {1.55,1.473,1.421,1.368,1.295,1.220,1.124};

  DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());

  //gROOT->cd();
  TFile* fOutputFile = new TFile(OutputFolder+TString::Format("/fOut_%s_SS%i_D%u_mT%i_ppSign%i.root",SourceVar.Data(),SmearStrategy,DataVar,imTbin,CorrectSign),"recreate");

  TH1F* HDATA = NULL;
  //GDATA->Set(TotMomBins);
  TGraph* GFIT = new TGraph();
  //GFIT->Set(TotMomBins);
  TGraph* GBL = new TGraph();
  //GBL->Set(TotMomBins);
  TGraph* GFEMTO = new TGraph();
  //GFEMTO->Set(TotMomBins);
  Int_t BASELINEVAR;
  Float_t FEMTOREGION;
  Int_t PS0;
  Int_t SMEARSTRATEGY;
  Int_t PP_LAM_VAR;
  TString MOMSMEARVAR;
  Float_t CKCUTOFF;
  Int_t PL_POT_VAR;
  TString SOURCEVAR;
  Float_t LAM_PP;
  Float_t LAM_PPL;
  Float_t NSIG;
  Float_t RADIUS;
  Float_t RADERR;

  TTree* ppTree = new TTree("ppTree","ppTree");
  ppTree->Branch("hData","TH1F",&HDATA,32000,0);//
  ppTree->Branch("gFit","TGraph",&GFIT,32000,0);//
  ppTree->Branch("gBl","TGraph",&GBL,32000,0);//
  ppTree->Branch("gFemto","TGraph",&GFEMTO,32000,0);//
  ppTree->Branch("BaselineVar",&BASELINEVAR,"BaselineVar/I");//
  ppTree->Branch("FemtoRegion",&FEMTOREGION,"FemtoRegion/F");//
  ppTree->Branch("pS0",&PS0,"pS0/I");//
  ppTree->Branch("SmearStrategy",&SMEARSTRATEGY,"SmearStrategy/I");//
  ppTree->Branch("pp_lam_var",&PP_LAM_VAR,"pp_lam_var/I");//
  ppTree->Branch("MomSmearVar","TString",&MOMSMEARVAR,8000,0);//
  ppTree->Branch("CkCutOff",&CKCUTOFF,"CkCutOff/F");//
  ppTree->Branch("pL_pot_var",&PL_POT_VAR,"pL_pot_var/I");//
  ppTree->Branch("SourceVar","TString",&SOURCEVAR,8000,0);//
  ppTree->Branch("lam_pp",&LAM_PP,"lam_pp/F");//
  ppTree->Branch("lam_ppL",&LAM_PPL,"lam_ppL/F");//
  ppTree->Branch("nsig",&NSIG,"nsig/F");//
  ppTree->Branch("rad",&RADIUS,"rad/F");//
  ppTree->Branch("raderr",&RADERR,"raderr/F");//

//printf("2\n");
  //for(TString varSource : SourceVar){
    CATS AB_pp;
    AB_pp.SetMomBins(NumMomBins_pp,0,FemtoRegion.back());
    AnalysisObject.SetUpCats_pp(AB_pp,"AV18",SourceVar,0,CorrectSign?-202:202);
    AB_pp.SetAnaSource(0,ExpectedRadii.at(imTbin));
    if(SourceVar.Contains("Levy")) AB_pp.SetAnaSource(1,1.7);
    else AB_pp.SetAnaSource(1,2.0);
    AB_pp.SetNotifications(CATS::nWarning);
    AB_pp.KillTheCat();
//printf("3\n");
    CATS AB_pXim;
    //same binning as pL, as we only use pXim as feed-down
    AB_pXim.SetMomBins(NumMomBins_feed,0,FemtoRegion.back());
    AnalysisObject.SetUpCats_pXim(AB_pXim,"pXim_HALQCDPaper2020","Gauss");
    AB_pXim.SetAnaSource(0,ExpectedRadii.at(imTbin)*Scale_pXi);
    AB_pXim.SetNotifications(CATS::nWarning);
    AB_pXim.KillTheCat();
//printf("4\n");
    CATS AB_pXi1530;
    AB_pXi1530.SetMomBins(NumMomBins_feed,0,FemtoRegion.back());
    AB_pXi1530.SetNotifications(CATS::nWarning);
    AnalysisObject.SetUpCats_pXim(AB_pXi1530,"pXim1530","Gauss");//McLevyNolan_Reso
    AB_pXi1530.SetAnaSource(0,ExpectedRadii.at(imTbin)*Scale_pXi);
    AB_pXi1530.KillTheCat();
//printf("5\n");
    CATS AB_pS0_Chiral;
    //the minus one is to to avoid going above 350 MeV, since we do not have the WF there
    AB_pS0_Chiral.SetMomBins(NumMomBins_feed,0,FemtoRegion.back());
    AnalysisObject.SetUpCats_pS0(AB_pS0_Chiral,"Chiral","Gauss");
    AB_pS0_Chiral.SetAnaSource(0,ExpectedRadii.at(imTbin)*Scale_pS0);
    AB_pS0_Chiral.SetNotifications(CATS::nWarning);
    AB_pS0_Chiral.KillTheCat();
//printf("6\n");
    for(int varPL : pL_pot_var){
      CATS AB_pL;
      AB_pL.SetMomBins(NumMomBins_feed,0,FemtoRegion.back());
      AnalysisObject.SetUpCats_pL(AB_pL,"Chiral_Coupled_SPD","Gauss",varPL,202);//NLO_Coupled_S
      const double CuspWeight = 0.33;//0.54
      if(abs(varPL)>1000){
          AB_pL.SetChannelWeight(7,1./4.*CuspWeight);//1S0 SN(s) -> LN(s)
          AB_pL.SetChannelWeight(8,3./4.*CuspWeight);//3S1 SN(s) -> LN(s)
          AB_pL.SetChannelWeight(10,3./4.*CuspWeight);//3S1 SN(d) -> LN(s)
          AB_pL.SetChannelWeight(13,3./20.*CuspWeight);//3D1 SN(d) -> LN(d)
          AB_pL.SetChannelWeight(15,3./20.*CuspWeight);//3D1 SN(s) -> LN(d)
      }
      AB_pL.SetAnaSource(0,ExpectedRadii.at(imTbin)*Scale_pL);
      AB_pL.SetNotifications(CATS::nError);
      AB_pL.KillTheCat();

      for(float varCutOff : CkCutOff){
        DLM_Ck* Ck_pp = new DLM_Ck(NumSourcePars, 0, AB_pp, TotMomBins, kMin, kMax);
        Ck_pp->SetCutOff(FemtoRegion.back(),varCutOff);
        DLM_Ck* Ck_pL = new DLM_Ck(NumSourcePars, 0, AB_pL, TotMomBins, kMin, kMax);
        Ck_pL->SetCutOff(FemtoRegion.at(0),varCutOff);
        DLM_Ck* Ck_pS0 = new DLM_Ck(NumSourcePars, 0, AB_pS0_Chiral, TotMomBins, kMin, kMax);
        Ck_pS0->SetCutOff(FemtoRegion.at(0),varCutOff);
        DLM_Ck* Ck_pXim = new DLM_Ck(NumSourcePars, 0, AB_pXim, TotMomBins, kMin, kMax);
        Ck_pXim->SetCutOff(FemtoRegion.at(0),varCutOff);
        DLM_Ck* Ck_pXim1530 = new DLM_Ck(NumSourcePars, 0, AB_pXi1530, TotMomBins, kMin, kMax);
        Ck_pXim1530->SetCutOff(FemtoRegion.at(0),varCutOff);

        for(TString varSmear : MomSmearVar){
          TH2F* hResolution_pp = AnalysisObject.GetResolutionMatrix(varSmear,"pp");
          TH2F* hResidual_pp_pL = AnalysisObject.GetResidualMatrix("pp","pLambda");
          TH2F* hResidual_pL_pSigma0 = AnalysisObject.GetResidualMatrix("pLambda","pSigma0");
          TH2F* hResidual_pL_pXim = AnalysisObject.GetResidualMatrix("pLambda","pXim");
          TH2F* hResidual_pXi_pXi1530 = AnalysisObject.GetResidualMatrix("pXim","pXim1530");

          for(int varLam : pp_lam_var){
            double lambda_pp[4];
            AnalysisObject.SetUpLambdaPars_pp("pp13TeV_HM_BernieSource",varLam+imTbin*10,lambda_pp);
            double lambda_pL[5];
            AnalysisObject.SetUpLambdaPars_pL("pp13TeV_HM_Dec19",varLam+imTbin*10,0,lambda_pL);
            double lambda_pXim[5];
            AnalysisObject.SetUpLambdaPars_pXim("pp13TeV_HM_Dec19",varLam+imTbin*10,0,lambda_pXim);
//printf("SetUp\n");
            //for(int varSS : SmearStrategy){
              TH1F* hPhaseSpace_pp=NULL;
              if(SmearStrategy==1){
                TList* list1_tmp;
                TList* list2_tmp;
                TString FileName = TString::Format("%s/CatsFiles/ExpData/Bernie_Source/ppData/mTBin_%i/CFOutput_mT_ppVar%u_HM_%i.root",
                                                    GetCernBoxDimi(),imTbin+1,DataVar,imTbin);
                TFile* inFile = new TFile(FileName,"read");
                //PARTICLES
                list1_tmp = (TList*)inFile->Get("PairDist");
                list2_tmp = (TList*)list1_tmp->FindObject("PairFixShifted");//there is also PairShifted ?? which one ??
                TH1F* hME_PP = (TH1F*)list2_tmp->FindObject(TString::Format("MEPart_mT_%i_FixShifted",imTbin));
                list1_tmp = (TList*)inFile->Get("AntiPairDist");
                list2_tmp = (TList*)list1_tmp->FindObject("PairFixShifted");
                TH1F* hME_APAP = (TH1F*)list2_tmp->FindObject(TString::Format("MEAntiPart_mT_%i_FixShifted",imTbin));
                gROOT->cd();
                hPhaseSpace_pp = (TH1F*)hME_PP->Clone("hPhaseSpace_pp");
                hPhaseSpace_pp->Add(hME_APAP);
                delete inFile;
              }
//printf("GotData\n");
              for(unsigned ipS0=!(pS0_Var); ipS0<2; ipS0++){
              //for(unsigned ipS0=0; ipS0<2; ipS0++){
                DLM_CkDecomposition CkDec_pp("pp",3,*Ck_pp,hResolution_pp);
                DLM_CkDecomposition CkDec_pL("pLambda", 4,*Ck_pL,NULL);
                DLM_CkDecomposition CkDec_pSigma0("pSigma0",0,*Ck_pS0,NULL);
                DLM_CkDecomposition CkDec_pXim("pXim",3,*Ck_pXim,NULL);
                DLM_CkDecomposition CkDec_pXim1530("pXim1530",0,*Ck_pXim1530,NULL);

                CkDec_pp.AddContribution(0,lambda_pp[1],DLM_CkDecomposition::cFeedDown,&CkDec_pL,hResidual_pp_pL);
                CkDec_pp.AddContribution(1,lambda_pp[2],DLM_CkDecomposition::cFeedDown);
                CkDec_pp.AddContribution(2,lambda_pp[3],DLM_CkDecomposition::cFake);
                if(hPhaseSpace_pp){
                  CkDec_pp.AddPhaseSpace(hPhaseSpace_pp);
                  CkDec_pp.AddPhaseSpace(0, hPhaseSpace_pp);
                }

                if(ipS0==0) CkDec_pL.AddContribution(0,lambda_pL[1],DLM_CkDecomposition::cFeedDown);
                else CkDec_pL.AddContribution(0,lambda_pL[1],DLM_CkDecomposition::cFeedDown,&CkDec_pSigma0,hResidual_pL_pSigma0);
                CkDec_pL.AddContribution(1, lambda_pL[2],DLM_CkDecomposition::cFeedDown,&CkDec_pXim,hResidual_pL_pXim);
                CkDec_pL.AddContribution(2, lambda_pL[3],DLM_CkDecomposition::cFeedDown);
                CkDec_pL.AddContribution(3, lambda_pL[4],DLM_CkDecomposition::cFake);
                if(hPhaseSpace_pp){
                  if(ipS0==1) CkDec_pL.AddPhaseSpace(0,hPhaseSpace_pp);
                  CkDec_pL.AddPhaseSpace(1,hPhaseSpace_pp);
                }

                CkDec_pXim.AddContribution(0, lambda_pXim[1],DLM_CkDecomposition::cFeedDown,&CkDec_pXim1530,hResidual_pXi_pXi1530);  //from Xi-(1530)
                CkDec_pXim.AddContribution(1, lambda_pXim[2]+lambda_pXim[3],DLM_CkDecomposition::cFeedDown);  //other feed-down (flat)
                CkDec_pXim.AddContribution(2, lambda_pXim[4], DLM_CkDecomposition::cFake);
                if(hPhaseSpace_pp){
                  CkDec_pXim.AddPhaseSpace(0,hPhaseSpace_pp);
                }

                CkDec_pp.Update();
                CkDec_pL.Update();
                CkDec_pSigma0.Update();
                CkDec_pXim.Update();
                CkDec_pXim1530.Update();
//printf("Updated\n");
                gROOT->cd();
                TH1F* hData = AnalysisObject.GetAliceExpCorrFun("pp13TeV_HM_BernieSource","pp",TString::Format("%u",DataVar),0,0,imTbin);
//printf("GetAliceExpCorrFun %p\n",hData);
                for(float varFit : FemtoRegion){
                  for(int varBL : BaselineVar){
                    double FitRange = varFit;
                    if(varBL==dpol3e) FitRange=BaselineRegion.at(0);
                    //gROOT->cd();
                    fOutputFile->cd();
                    //2 femto and 5 BL fit pars (the BL are norm, pol1,2,3,4)
                    TF1* fData = new TF1("fData",Fit_BernieSource,0,FitRange,7);
                    TF1* fBl = new TF1("fBl",Baseline_BernieSource,0,FitRange,5);
                    fData->SetParameter(0,SourceVar=="Gauss"?ExpectedRadii.at(imTbin):ExpectedRadii.at(imTbin)*Scale_core);
                    fData->SetParLimits(0,fData->GetParameter(0)*0.5,fData->GetParameter(0)*2.0);
                    //fData->SetParLimits(0,fData->GetParameter(0)*0.75,fData->GetParameter(0)*1.25);
//fData->FixParameter(0,1.3);
//fData->SetParameter(0,1.366);
//fData->SetParLimits(0,1.1,1.5);
                    fData->FixParameter(1,0);

                    //BL
                    fData->SetParameter(2,1);

                    //pol0s,pol1s,dpol3e
                    if(varBL==pol1s){
                      fData->SetParameter(3,0);
                      fData->SetParLimits(3,-1e-2,1e-2);
                      fData->FixParameter(4,0);
                      fData->FixParameter(5,0);
                      fData->FixParameter(6,0);
                    }
                    else if(varBL==dpol3e){
                      fData->SetParameter(3,0);
                      fData->SetParLimits(3,-100000,100);
                      fData->SetParameter(4,100);
                      fData->SetParLimits(4,0,400);
                      fData->SetParameter(5,0);
                      fData->SetParLimits(5,-1e-6,1e-6);
                      fData->FixParameter(6,0);
                    }
                    else{
                      fData->FixParameter(3,0);
                      fData->FixParameter(4,0);
                      fData->FixParameter(5,0);
                      fData->FixParameter(6,0);
                    }
                    SOURCE_FIT = &CkDec_pp;
                    printf("BL=%i FIT=%.0f PS0=%i SS=%i LAM=%i SMR=%s PL=%i SRC=%s lam_pp=%.1f lam_ppl=%.1f\n",
                            varBL,varFit,ipS0,SmearStrategy,varLam,varSmear.Data(),varPL,SourceVar.Data(),lambda_pp[0]*100.,lambda_pp[1]*100.);
                    hData->Fit(fData,"Q, S, N, R, M");
//fData->FixParameter(0,1.4);
//fData->FixParameter(2,1);
//TFile TempF(OutputFolder+"/TempF.root","recreate");
//hData->Write();
//fData->Write();
//CkDec_pp.Update(true,true);
//printf("C(20) = %f; cats %f; Ck_pp %f; CkDec_pp %f\n",fData->Eval(20),AB_pp.EvalCorrFun(20),Ck_pp->Eval(20),CkDec_pp.EvalCk(20));
//return;
//printf("fitted\n");
                    fData->SetNpx(TotMomBins);

                    fBl->FixParameter(0,fData->GetParameter(2));
                    fBl->FixParameter(1,fData->GetParameter(3));
                    fBl->FixParameter(2,fData->GetParameter(4));
                    fBl->FixParameter(3,fData->GetParameter(5));
                    fBl->FixParameter(4,fData->GetParameter(6));

printf(" r = %.3f +/- %.3f\n",fData->GetParameter(0),fData->GetParError(0));
                    double Chi2 = fData->GetChisquare();
                    double NDF = fData->GetNDF();
                    double pval = TMath::Prob(Chi2,NDF);
                    double nsig = sqrt(2)*TMath::ErfcInverse(pval);
                    printf(" chi2/ndf = %.2f\n",Chi2/NDF);
                    printf(" nsig = %.2f\n",nsig);

                    //gROOT->cd();
                    fOutputFile->cd();

                    BASELINEVAR = varBL;
                    FEMTOREGION = varFit;
                    PS0 = ipS0;
                    SMEARSTRATEGY = SmearStrategy;
                    PP_LAM_VAR = varLam;
                    MOMSMEARVAR = varSmear;
                    CKCUTOFF = varCutOff;
                    PL_POT_VAR = varPL;
                    SOURCEVAR = SourceVar;
                    LAM_PP = lambda_pp[0];
                    LAM_PPL = lambda_pp[1];
                    NSIG = nsig;
                    RADIUS = fData->GetParameter(0);
                    RADERR = fData->GetParError(0);

                    if(HDATA)delete HDATA;
                    HDATA = (TH1F*)hData->Clone("HDATA");
                    delete GFIT;
                    GFIT = new TGraph();
                    GFIT->SetName("GFIT");
                    delete GBL;
                    GBL = new TGraph();
                    GBL->SetName("GBL");
                    delete GFEMTO;
                    GFEMTO = new TGraph();
                    GFEMTO->SetName("GFEMTO");
                    for(unsigned uBin=0; uBin<NumMomBins_pp; uBin++){
                      double MOM = AB_pp.GetMomentum(uBin);
                      if(MOM>FitRange) break;
                      //GDATA->SetPoint(uBin,MOM,hData->GetBinContent(hData->FindBin(MOM)));
                      GFIT->SetPoint(uBin,MOM,fData->Eval(MOM));
                      GBL->SetPoint(uBin,MOM,fBl->Eval(MOM));
                      GFEMTO->SetPoint(uBin,MOM,fData->Eval(MOM)/fBl->Eval(MOM));
                    }

                    ppTree->Fill();
                    //ppTree->Branch("gData","TH1F",&GDATA,32000,0);//
                    //ppTree->Branch("gFit","TGraph",&GFIT,32000,0);//
                    //ppTree->Branch("gBl","TGraph",&GBL,32000,0);//
                    //ppTree->Branch("gFemto","TGraph",&GFEMTO,32000,0);//

//save output here: OutputFolder

                    TFile outTest(OutputFolder+"/outTest_OLD.root","recreate");

                    TGraph* GFEMTO_TH = new TGraph();
                    TGraph* GFEMTO_SF = new TGraph();
                    GFEMTO_TH->SetName("GFEMTO_TH");
                    GFEMTO_SF->SetName("GFEMTO_SF");
                    CkDec_pp.GetCk()->SetSourcePar(0,RADIUS);
                    CkDec_pp.Update(true,true);
                    for(unsigned uBin=0; uBin<NumMomBins_pp; uBin++){
                      double MOM = AB_pp.GetMomentum(uBin);
                      if(MOM>FitRange) break;
                      GFEMTO_TH->SetPoint(uBin,MOM,CkDec_pp.EvalSignal(MOM)+1);
                      GFEMTO_SF->SetPoint(uBin,MOM,SOURCE_FIT->EvalCk(MOM)+1);
                    }

                    //printf("File opened\n");
                    hData->Write();
                    //printf("Histo written\n");
                    GFIT->Write();
                    GBL->Write();
                    GFEMTO->Write();
                    GFEMTO_TH->Write();
                    GFEMTO_SF->Write();
                    delete GFEMTO_TH;


                    //return;
                    delete fData;
                    delete fBl;
//break;
                  }//varBL (3x)
                }//varFit (3x)
                delete hData;
//break;
              }//ipS0 (2x)
            //}//varSS (2x)
//break;
          }//varLam (3x)
//break;
        }//varSmear (2x)
        delete Ck_pp;
        delete Ck_pL;
        delete Ck_pS0;
        delete Ck_pXim;
        delete Ck_pXim1530;
//break;
      }//varCutOff (1x)
//break;
    }//varPL (2x)
  //}//varSource (2x)

//const unsigned DataVar, const int imTbin
  fOutputFile->cd();
  ppTree->Write();

  delete ppTree;
  delete fOutputFile;
}





void SourcePaper_pL(const TString SourceVar, const int& SmearStrategy, const unsigned DataVar, const int imTbin, const TString OutputFolder){
  gROOT->ProcessLine("gErrorIgnoreLevel = 2001;");
  const bool Silent = true;
  //std::vector<TString> SourceVar = {"Gauss","McGauss_ResoTM"};//done
//std::vector<TString> SourceVar = {"McLevy_ResoTM"};
  //std::vector<int> pp_lam_var = {0,1,2};//done
std::vector<int> pp_lam_var = {0};
//std::vector<int> pp_lam_var = {0};
  std::vector<int> pL_lam_var = {00,01,02,10,11,12,20,21,22};
  //0 is the old one, 1 is the new one with the folded ME
  //std::vector<int> SmearStrategy = {0,1};//done
  //std::vector<float> FemtoRegion = {376,352,400};//done
std::vector<float> FemtoRegion = {204,228,240};
  const unsigned NumMomBins_pL = TMath::Nint(FemtoRegion.back()/12.);
  const unsigned NumMomBins_feed = TMath::Nint(FemtoRegion.back()/12.);
  //std::vector<float> BaselineRegion = {500};
  std::vector<float> CkCutOff = {700};//done
  std::vector<int> pL_pot_var = {11600,-11600};//done
//std::vector<int> pL_pot_var = {11600};
  //pp13TeV_HM_DimiJun20 -> the ME reso matrix
  //pp13TeV_HM_BernieSource -> Bernie's old smearing matrix
  //pp13TeV_HM_BernieSource -> contains the data
  std::vector<TString> MomSmearVar = {"pp13TeV_HM_DimiJun20"};//done
  //and an additional variation will be done on if pS0 is included as a feed or not
  enum BLTYPE { pol0s,pol1s,pol2s,pol3s,dpol2s,dpol3s,dpol4s,pol2e,pol3e,dpol2e,dpol3e,dpol4e,spl1 };
  std::vector<int> BaselineVar = {pol0s,pol1s};
//std::vector<int> BaselineVar = {pol0s,pol1s};
  const bool pS0_Var = true;

  //the region for which the DLM_Ck objects will be defined
  //we put some extra to get away from the edge effects of the smearing
  const double kMin = 0;
  const double kMax = 360;
  const unsigned TotMomBins = TMath::Nint(kMax);
  const unsigned NumMtBins = 6;
  if(imTbin>=NumMtBins){
    printf("\033[1;31mERROR:\033[0m Only %u mT bins are available!\n",NumMtBins);
    return;
  }
  //approximate, based on the output we kind of know we will get
  //std::vector<float> ExpectedRadii = {1.35,1.3,1.25,1.2,1.1,1.05,0.95};
  //const double Scale_pL = 1.1;
  //const double Scale_pS0 = 1.15;
  //const double Scale_pXi = 0.97;
  //const double Scale_core = 0.94;
  std::vector<float> ExpectedRadii = { 1.55, 1.473, 1.421, 1.368, 1.295, 1.220, 1.124 };
  const double Scale_pL = 1.0;
  const double Scale_pS0 = 1.0;
  const double Scale_pXi = 1.0;
  const double Scale_core = 0.85;
  //the difference in the effectiv Gaussian compered to pp

  const unsigned NumSourcePars = 1;
  //std::vector<float> pSigma0Radii = {1.55,1.473,1.421,1.368,1.295,1.220,1.124};

  DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());

  //gROOT->cd();
  TFile* fOutputFile = new TFile(OutputFolder+TString::Format("/fOut_%s_SS%i_D%u_mT%i.root",SourceVar.Data(),SmearStrategy,DataVar,imTbin),"recreate");

  TH1F* HDATA = NULL;
  //GDATA->Set(TotMomBins);
  TGraph* GFIT = new TGraph();
  //GFIT->Set(TotMomBins);
  TGraph* GBL = new TGraph();
  //GBL->Set(TotMomBins);
  TGraph* GFEMTO = new TGraph();
  //GFEMTO->Set(TotMomBins);
  Int_t BASELINEVAR;
  Float_t FEMTOREGION;
  Int_t PS0;
  Int_t SMEARSTRATEGY;
  Int_t PL_LAM_VAR;
  TString MOMSMEARVAR;
  Float_t CKCUTOFF;
  Int_t PL_POT_VAR;
  TString SOURCEVAR;
  Float_t LAM_PL;
  //Float_t LAM_PPL;
  Float_t NSIG;
  Float_t RADIUS;
  Float_t RADERR;

  TTree* pLTree = new TTree("pLTree","pLTree");
  pLTree->Branch("hData","TH1F",&HDATA,32000,0);//
  pLTree->Branch("gFit","TGraph",&GFIT,32000,0);//
  pLTree->Branch("gBl","TGraph",&GBL,32000,0);//
  pLTree->Branch("gFemto","TGraph",&GFEMTO,32000,0);//
  pLTree->Branch("BaselineVar",&BASELINEVAR,"BaselineVar/I");//
  pLTree->Branch("FemtoRegion",&FEMTOREGION,"FemtoRegion/F");//
  pLTree->Branch("pS0",&PS0,"pS0/I");//
  pLTree->Branch("SmearStrategy",&SMEARSTRATEGY,"SmearStrategy/I");//
  pLTree->Branch("pL_lam_var",&PL_LAM_VAR,"pL_lam_var/I");//
  pLTree->Branch("MomSmearVar","TString",&MOMSMEARVAR,8000,0);//
  pLTree->Branch("CkCutOff",&CKCUTOFF,"CkCutOff/F");//
  pLTree->Branch("pL_pot_var",&PL_POT_VAR,"pL_pot_var/I");//
  pLTree->Branch("SourceVar","TString",&SOURCEVAR,8000,0);//
  pLTree->Branch("lam_pL",&LAM_PL,"lam_pL/F");//
  //ppTree->Branch("lam_ppL",&LAM_PPL,"lam_ppL/F");//
  pLTree->Branch("nsig",&NSIG,"nsig/F");//
  pLTree->Branch("rad",&RADIUS,"rad/F");//
  pLTree->Branch("raderr",&RADERR,"raderr/F");//

    CATS AB_pXim;
    //same binning as pL, as we only use pXim as feed-down
    AB_pXim.SetMomBins(NumMomBins_feed,0,FemtoRegion.back());
    AnalysisObject.SetUpCats_pXim(AB_pXim,"pXim_HALQCDPaper2020","Gauss");
    AB_pXim.SetAnaSource(0,ExpectedRadii.at(imTbin)*Scale_pXi);
    AB_pXim.SetNotifications(CATS::nWarning);
    AB_pXim.KillTheCat();

    CATS AB_pXi1530;
    AB_pXi1530.SetMomBins(NumMomBins_feed,0,FemtoRegion.back());
    AB_pXi1530.SetNotifications(CATS::nWarning);
    AnalysisObject.SetUpCats_pXim(AB_pXi1530,"pXim1530","Gauss");//McLevyNolan_Reso
    AB_pXi1530.SetAnaSource(0,ExpectedRadii.at(imTbin)*Scale_pXi);
    AB_pXi1530.KillTheCat();

    CATS AB_pS0_Chiral;
    //the minus one is to to avoid going above 350 MeV, since we do not have the WF there
    AB_pS0_Chiral.SetMomBins(NumMomBins_feed,0,FemtoRegion.back());
    AnalysisObject.SetUpCats_pS0(AB_pS0_Chiral,"Chiral","Gauss");
    AB_pS0_Chiral.SetAnaSource(0,ExpectedRadii.at(imTbin)*Scale_pS0);
    AB_pS0_Chiral.SetNotifications(CATS::nWarning);
    AB_pS0_Chiral.KillTheCat();

    for(int varPL : pL_pot_var){
      CATS AB_pL;
      AB_pL.SetMomBins(NumMomBins_feed,0,FemtoRegion.back());
      AnalysisObject.SetUpCats_pL(AB_pL,"Chiral_Coupled_SPD",SourceVar,varPL,202);//NLO_Coupled_S
      const double CuspWeight = 0.33;//0.54
      if(abs(varPL)>1000){
          AB_pL.SetChannelWeight(7,1./4.*CuspWeight);//1S0 SN(s) -> LN(s)
          AB_pL.SetChannelWeight(8,3./4.*CuspWeight);//3S1 SN(s) -> LN(s)
          AB_pL.SetChannelWeight(10,3./4.*CuspWeight);//3S1 SN(d) -> LN(s)
          AB_pL.SetChannelWeight(13,3./20.*CuspWeight);//3D1 SN(d) -> LN(d)
          AB_pL.SetChannelWeight(15,3./20.*CuspWeight);//3D1 SN(s) -> LN(d)
      }
      if(SourceVar=="Gauss") AB_pL.SetAnaSource(0,ExpectedRadii.at(imTbin));
      else AB_pL.SetAnaSource(0,ExpectedRadii.at(imTbin)*Scale_core);
      AB_pL.SetNotifications(CATS::nError);
      AB_pL.KillTheCat();

      for(float varCutOff : CkCutOff){
        DLM_Ck* Ck_pL = new DLM_Ck(NumSourcePars, 0, AB_pL, TotMomBins, kMin, kMax);
        Ck_pL->SetCutOff(FemtoRegion.at(0),varCutOff);
        DLM_Ck* Ck_pS0 = new DLM_Ck(NumSourcePars, 0, AB_pS0_Chiral, TotMomBins, kMin, kMax);
        Ck_pS0->SetCutOff(FemtoRegion.at(0),varCutOff);
        DLM_Ck* Ck_pXim = new DLM_Ck(NumSourcePars, 0, AB_pXim, TotMomBins, kMin, kMax);
        Ck_pXim->SetCutOff(FemtoRegion.at(0),varCutOff);
        DLM_Ck* Ck_pXim1530 = new DLM_Ck(NumSourcePars, 0, AB_pXi1530, TotMomBins, kMin, kMax);
        Ck_pXim1530->SetCutOff(FemtoRegion.at(0),varCutOff);

        for(TString varSmear : MomSmearVar){
          TH2F* hResolution_pL = AnalysisObject.GetResolutionMatrix(varSmear,"pLambda");
          TH2F* hResidual_pp_pL = AnalysisObject.GetResidualMatrix("pp","pLambda");
          TH2F* hResidual_pL_pSigma0 = AnalysisObject.GetResidualMatrix("pLambda","pSigma0");
          TH2F* hResidual_pL_pXim = AnalysisObject.GetResidualMatrix("pLambda","pXim");
          TH2F* hResidual_pXi_pXi1530 = AnalysisObject.GetResidualMatrix("pXim","pXim1530");

          for(int varLam : pL_lam_var){
            double lambda_pL[5];
            AnalysisObject.SetUpLambdaPars_pL("pp13TeV_HM_Dec19",0,varLam,lambda_pL);
            double lambda_pXim[5];
            AnalysisObject.SetUpLambdaPars_pXim("pp13TeV_HM_Dec19",0,0,lambda_pXim);

            //for(int varSS : SmearStrategy){
              TH1F* hPhaseSpace_pL=NULL;
              if(SmearStrategy==1){
                TList* list1_tmp;
                TList* list2_tmp;
                TString FileName = TString::Format("%s/CatsFiles/ExpData/Bernie_Source/pLData/mTBin_%i/CFOutput_mT_pLVar%u_HM_%i.root",
                                                    GetCernBoxDimi(),imTbin+1,DataVar,imTbin);
                TFile* inFile = new TFile(FileName,"read");
                //PARTICLES
                list1_tmp = (TList*)inFile->Get("PairDist");
                list2_tmp = (TList*)list1_tmp->FindObject("PairFixShifted");//there is also PairShifted ?? which one ??
                TH1F* hME_PP = (TH1F*)list2_tmp->FindObject(TString::Format("MEPart_mT_%i_FixShifted",imTbin));
                list1_tmp = (TList*)inFile->Get("AntiPairDist");
                list2_tmp = (TList*)list1_tmp->FindObject("PairFixShifted");
                TH1F* hME_APAP = (TH1F*)list2_tmp->FindObject(TString::Format("MEAntiPart_mT_%i_FixShifted",imTbin));
                gROOT->cd();
                hPhaseSpace_pL = (TH1F*)hME_PP->Clone("hPhaseSpace_pL");
                hPhaseSpace_pL->Add(hME_APAP);
                delete inFile;
              }
              for(unsigned ipS0=!(pS0_Var); ipS0<2; ipS0++){
              //for(unsigned ipS0=0; ipS0<2; ipS0++){
                DLM_CkDecomposition CkDec_pL("pLambda", 4,*Ck_pL,hResolution_pL);
                DLM_CkDecomposition CkDec_pSigma0("pSigma0",0,*Ck_pS0,NULL);
                DLM_CkDecomposition CkDec_pXim("pXim",3,*Ck_pXim,NULL);
                DLM_CkDecomposition CkDec_pXim1530("pXim1530",0,*Ck_pXim1530,NULL);

                if(ipS0==0) CkDec_pL.AddContribution(0,lambda_pL[1],DLM_CkDecomposition::cFeedDown);
                else CkDec_pL.AddContribution(0,lambda_pL[1],DLM_CkDecomposition::cFeedDown,&CkDec_pSigma0,hResidual_pL_pSigma0);
                CkDec_pL.AddContribution(1, lambda_pL[2],DLM_CkDecomposition::cFeedDown,&CkDec_pXim,hResidual_pL_pXim);
                CkDec_pL.AddContribution(2, lambda_pL[3],DLM_CkDecomposition::cFeedDown);
                CkDec_pL.AddContribution(3, lambda_pL[4],DLM_CkDecomposition::cFake);
                if(hPhaseSpace_pL){
                  if(ipS0==1) CkDec_pL.AddPhaseSpace(0,hPhaseSpace_pL);
                  CkDec_pL.AddPhaseSpace(1,hPhaseSpace_pL);
                }

                CkDec_pXim.AddContribution(0, lambda_pXim[1],DLM_CkDecomposition::cFeedDown,&CkDec_pXim1530,hResidual_pXi_pXi1530);  //from Xi-(1530)
                CkDec_pXim.AddContribution(1, lambda_pXim[2]+lambda_pXim[3],DLM_CkDecomposition::cFeedDown);  //other feed-down (flat)
                CkDec_pXim.AddContribution(2, lambda_pXim[4], DLM_CkDecomposition::cFake);
                if(hPhaseSpace_pL){
                  CkDec_pXim.AddPhaseSpace(0,hPhaseSpace_pL);
                }

                CkDec_pL.Update();
                CkDec_pSigma0.Update();
                CkDec_pXim.Update();
                CkDec_pXim1530.Update();

                gROOT->cd();
                TH1F* hData = AnalysisObject.GetAliceExpCorrFun("pp13TeV_HM_BernieSource","pLambda",TString::Format("%u",DataVar),0,0,imTbin);
                for(float varFit : FemtoRegion){
                  for(int varBL : BaselineVar){
                    double FitRange = varFit;
                    //if(varBL==dpol3e) FitRange=BaselineRegion.at(0);
                    //gROOT->cd();
                    fOutputFile->cd();
                    //2 femto and 5 BL fit pars (the BL are norm, pol1,2,3,4)
                    TF1* fData = new TF1("fData",Fit_BernieSource,0,FitRange,7);
                    TF1* fBl = new TF1("fBl",Baseline_BernieSource,0,FitRange,5);
                    fData->SetParameter(0,SourceVar=="Gauss"?ExpectedRadii.at(imTbin):ExpectedRadii.at(imTbin)*Scale_core);
                    fData->SetParLimits(0,fData->GetParameter(0)*0.5,fData->GetParameter(0)*2.0);
//fData->FixParameter(0,1.366);
                    fData->FixParameter(1,0);

                    //BL
                    fData->SetParameter(2,1);

                    //pol0s,pol1s,dpol3e
                    if(varBL==pol1s){
                      fData->SetParameter(3,0);
                      fData->SetParLimits(3,-1e-2,1e-2);
                      fData->FixParameter(4,0);
                      fData->FixParameter(5,0);
                      fData->FixParameter(6,0);
                    }
                    else if(varBL==dpol3e){
                      fData->SetParameter(3,0);
                      fData->SetParLimits(3,-100000,100);
                      fData->SetParameter(4,100);
                      fData->SetParLimits(4,0,400);
                      fData->SetParameter(5,0);
                      fData->SetParLimits(5,-1e-6,1e-6);
                      fData->FixParameter(6,0);
                    }
                    else{
                      fData->FixParameter(3,0);
                      fData->FixParameter(4,0);
                      fData->FixParameter(5,0);
                      fData->FixParameter(6,0);
                    }
                    SOURCE_FIT = &CkDec_pL;
                    //printf("BL=%i FIT=%.0f PS0=%i SS=%i LAM=%i SMR=%s PL=%i SRC=%s lam_pL=%.1f XXX=%.1f\n",
                    //        varBL,varFit,ipS0,SmearStrategy,varLam,varSmear.Data(),varPL,SourceVar.Data(),lambda_pL[0]*100.,lambda_pL[1]*100.);
                    hData->Fit(fData,"Q, S, N, R, M");
                    fData->SetNpx(TotMomBins);

                    fBl->FixParameter(0,fData->GetParameter(2));
                    fBl->FixParameter(1,fData->GetParameter(3));
                    fBl->FixParameter(2,fData->GetParameter(4));
                    fBl->FixParameter(3,fData->GetParameter(5));
                    fBl->FixParameter(4,fData->GetParameter(6));

                    printf(" r = %.3f +/- %.3f\n",fData->GetParameter(0),fData->GetParError(0));
                    double Chi2 = fData->GetChisquare();
                    double NDF = fData->GetNDF();
                    double pval = TMath::Prob(Chi2,NDF);
                    double nsig = sqrt(2)*TMath::ErfcInverse(pval);
                    //printf(" chi2/ndf = %.2f\n",Chi2/NDF);
                    //printf(" nsig = %.2f\n",nsig);

                    //gROOT->cd();
                    fOutputFile->cd();

                    BASELINEVAR = varBL;
                    FEMTOREGION = varFit;
                    PS0 = ipS0;
                    SMEARSTRATEGY = SmearStrategy;
                    PL_LAM_VAR = varLam;
                    MOMSMEARVAR = varSmear;
                    CKCUTOFF = varCutOff;
                    PL_POT_VAR = varPL;
                    SOURCEVAR = SourceVar;
                    LAM_PL = lambda_pL[0];
                    //LAM_PPL = lambda_pp[1];
                    NSIG = nsig;
                    RADIUS = fData->GetParameter(0);
                    RADERR = fData->GetParError(0);

                    if(HDATA)delete HDATA;
                    HDATA = (TH1F*)hData->Clone("HDATA");
                    delete GFIT;
                    GFIT = new TGraph();
                    GFIT->SetName("GFIT");
                    delete GBL;
                    GBL = new TGraph();
                    GBL->SetName("GBL");
                    delete GFEMTO;
                    GFEMTO = new TGraph();
                    GFEMTO->SetName("GFEMTO");
                    for(unsigned uBin=0; uBin<NumMomBins_pL; uBin++){
                      double MOM = AB_pL.GetMomentum(uBin);
                      if(MOM>FitRange) break;
                      //GDATA->SetPoint(uBin,MOM,hData->GetBinContent(hData->FindBin(MOM)));
                      GFIT->SetPoint(uBin,MOM,fData->Eval(MOM));
                      GBL->SetPoint(uBin,MOM,fBl->Eval(MOM));
                      GFEMTO->SetPoint(uBin,MOM,fData->Eval(MOM)/fBl->Eval(MOM));
                    }

                    pLTree->Fill();
                    //ppTree->Branch("gData","TH1F",&GDATA,32000,0);//
                    //ppTree->Branch("gFit","TGraph",&GFIT,32000,0);//
                    //ppTree->Branch("gBl","TGraph",&GBL,32000,0);//
                    //ppTree->Branch("gFemto","TGraph",&GFEMTO,32000,0);//

//save output here: OutputFolder
/*
                    TFile outTest(OutputFolder+"/outTest_OLD.root","recreate");

                    TGraph* GFEMTO_TH = new TGraph();
                    GFEMTO_TH->SetName("GFEMTO_TH");
                    CkDec_pp.GetCk()->SetSourcePar(0,RADIUS);
                    CkDec_pp.Update(true,true);
                    for(unsigned uBin=0; uBin<NumMomBins_pp; uBin++){
                      double MOM = AB_pp.GetMomentum(uBin);
                      if(MOM>FitRange) break;
                      GFEMTO_TH->SetPoint(uBin,MOM,CkDec_pp.EvalSignal(MOM)+1);
                    }

                    //printf("File opened\n");
                    hData->Write();
                    //printf("Histo written\n");
                    GFIT->Write();
                    GBL->Write();
                    GFEMTO->Write();
                    GFEMTO_TH->Write();
                    delete GFEMTO_TH;
*/
                    //return;
                    delete fData;
                    delete fBl;
//break;
                  }//varBL (3x)
                }//varFit (3x)
                delete hData;
//break;
              }//ipS0 (2x)
            //}//varSS (2x)
//break;
          }//varLam (3x)
//break;
        }//varSmear (2x)
        delete Ck_pL;
        delete Ck_pS0;
        delete Ck_pXim;
        delete Ck_pXim1530;
//break;
      }//varCutOff (1x)
//break;
    }//varPL (2x)
  //}//varSource (2x)

//const unsigned DataVar, const int imTbin
  fOutputFile->cd();
  pLTree->Write();

  delete pLTree;
  delete fOutputFile;
}



//input: mT bin, which data variation, pL intercation, and AnalysisType
//if the mT bin is -1, we fit the total pL correlation (as in paper, only the folded version)
//the pL intercation is (pL_inter):
//      -1: LO 600 only s-wave (as orginally in the paper)
//      0: NLO-600 only s-wave (as orginally in the paper)
//      potential ID as in the pL paper (s and d waves). In practice we will try out:
//      11600: NLO-600
//      131600: NLO-600 with reduced scattering length in the 3S1
//      231600: NLO-600 with reduced scattering length in the 3S1, and charge symmetry breaking (CSB) in 1S0
//AnalysisType consisting of 4 digits (ABCD):
//      A: the smearing strategy, 0 is orginal, 1 is improved based on limited phase space (more accurate)
//      B:  the feed-down, 0 (fss2) and 1 (esc16) are the two original variations for pS0 + the old pXi,
//          while 2 and 3 are the two types used in the pL paper: 2 (0) is pS0 with chEFT , 3 (1) is pS0 flat, in both cases latest pXi is used (or the old one)
//      C: source type. 0 is Gauss, 1 is RSM with EPOS disto, 2 is RSM with flat disto (the latter not done yet)
//      D: related to the fit range, baseline and lambda pars. 0 is the default (up to c.a. 228) while 1 is as in the pL paper (up to c.a. 480).
//         in the original version we use for a baseline either pol0 or pol1, in the new version we ALWAYS use the pol3 without linear term
//         for the lambda pars, 0 is the default in the source paper (not implemented yet), 1 is as used for the pL paper
//      examples: ABCD = 0 is the exact settings as the original, fitting with Gauss
//      ABCD = 1201 are the settings for a Gaussian fit using the set up as in the pL paper
//      ABCD = 1211 are the exact settings of the pL paper (i.e. where RSM was used)
//      N.B. in any case, the one major difference to the pL paper is that we do NOT unfold the data, but apply the momentum smearing on top
//NewAgeFlag: placeholder for anything new we might thing of. So far:
//      last digit: Variations of the feed-down radius.
//                  0: fixed as in the source paper, 1 is equal to the value of reff from the mT(pp) result,
//                  2 is a random value sampled between <r_core(pp)> (old new) and r_eff (pL) as in the source paper. This is a conservative thing to do,
//                  as the eff pS0 should be similar, or a bit larger, than pp, while pXi should be just slightly above the pp_core
//N.B. RndSeed == 1 will make the default variation first
//ErrorType = XY, where X is syst vars, Y is bootstrap
//Timeout in minutes, if negative => DEBUG MODE
void Jaime_pL(int RndSeed, int imTbin, int pL_inter, int AnalysisType, int ErrorType, double Timeout,
  int NewAgeFlag, const char* InputFolder, const char* OutputFolder){

  bool DEBUG = Timeout<0;
  Timeout = fabs(Timeout);

  if(DEBUG){
    printf("Starting Jaime_pL (%i %i %i %i %i %.2f %i %s %s)\n",RndSeed,imTbin,pL_inter,AnalysisType,ErrorType,Timeout,NewAgeFlag,InputFolder,OutputFolder);
  }

  Timeout *= 60;

  DLM_Timer JaimeTimer;

  gROOT->ProcessLine("gErrorIgnoreLevel = 2001;");
  const bool Silent = true;
  TRandom3 RanGen(RndSeed);

  int SystVar = (ErrorType/10)%10;
  int Bootstrap = (ErrorType)%10;

  TString SourceType;
  unsigned NumSourcePars;
  int SourceVar;
  //C: source type. 0 is Gauss, 1 is RSM with EPOS disto, 2 is RSM with flat disto (the latter not done yet)
  if((AnalysisType/10)%10 == 0) {SourceType = "Gauss"; SourceVar = 0; NumSourcePars = 1;}
  else if((AnalysisType/10)%10 == 1) {SourceType = "McGauss_ResoTM"; SourceVar = 202; NumSourcePars = 2;}
  //else if((AnalysisType/10)%10 == 2) {SourceType = "McGauss_ResoTM"; SourceVar = 202;}
  else {printf("Unknown SourceType, aborting!\n"); abort();}

  if(NewAgeFlag%10<0 || NewAgeFlag%10>2){
    printf("\033[1;31mERROR:\033[0m NewAgeFlag, aborting!\n");
    abort();
  }


  int SmearStrategy = (AnalysisType/1000)%10;
  int FeedDownType = (AnalysisType/100)%10;
  int FitType = AnalysisType%10;

  if(SmearStrategy<0 || SmearStrategy>1){
    printf("\033[1;31mERROR:\033[0m SmearStrategy, aborting!\n");
    abort();
  }
  if(DEBUG){
    printf(" SmearStrategy = %i\n",SmearStrategy);
    printf(" FeedDownType = %i\n",FeedDownType);
    printf(" FitType = %i\n",FitType);
    printf(" SourceType = %s\n",SourceType.Data());
  }

  std::vector<float> fitregion_var;
  if(FitType==0){
    fitregion_var.push_back(228);
    fitregion_var.push_back(204);
    fitregion_var.push_back(240);
  }
  else if(FitType==1){
    fitregion_var.push_back(456);
    fitregion_var.push_back(432);
    fitregion_var.push_back(480);
  }
  else{printf("Unknown FitType, aborting!\n"); abort();}

  if(DEBUG){
    printf(" fitregion_var = { ");
    for( float& frv : fitregion_var ){
      printf("%.0f ",frv);
    }
    printf("}\n");
  }

  //DCBA:
  //A:
  //the purity/fraction variations of protons
  //purity is fixed, last digit is the modulation of Lambda:Sigma feeddown ratio
  //second to last digit is the mt bin
  std::vector<int> lampar_p_var;
  if(imTbin<0){
    lampar_p_var.push_back(0);
    lampar_p_var.push_back(1);
    lampar_p_var.push_back(2);
  }
  else{
    //in the old analysis the primary proton fraction was evaluated
    //differentially according to the mtbin. Effect only minor.
    lampar_p_var.push_back(0+10*imTbin);
    lampar_p_var.push_back(1+10*imTbin);
    lampar_p_var.push_back(2+10*imTbin);
  }

  if(DEBUG){
    printf(" lampar_p_var = { ");
    for( int& lpv : lampar_p_var ){
      printf("%i ",lpv);
    }
    printf("}\n");
  }

  std::vector<int> lampar_l_var;
  if(FitType==0){
    //fixed purity, variation of BA:
    //  A is the S0:L, 0 is 1:3, 2 is 40% (kind of our default looking at best fits) 1 is 27%
    //  B: the Xi feed-down, we actually assume that it is fixed, however we change the Xi0 to Xim ratio
    //      this is stupid as it should be 50:50, however in the old analysis we assumed Xi0 was flat and varied this ratio as systematics...
    std::vector<int> var_A = {2,1,0};
    std::vector<int> var_B = {0,1,2};
    for(int& iA : var_A){
      for(int& iB : var_B){
        lampar_l_var.push_back( iA+10*iB );
      }
    }
  }
  else{
    //CBA:  A is the S0:L, 0 is 1:3, 2 is 40% (kind of our default looking at best fits) 1 is 27%
    //      B: the Xi feed-down, 1 and 2 are some old settings, now working with either 0 (fractions averaged over pT)
    //          or 3, which is our new default, where the fraction of lambdas is lowered based on the avg pT of the femto pairs
    //      C: the purity, 4 is 95.3%, 5 is 96.3%, these values are based on the PLB result
    std::vector<int> var_A = {2,1,0};
    std::vector<int> var_B = {3,0};
    std::vector<int> var_C = {4,5};
    for(int& iA : var_A){
      for(int& iB : var_B){
        for(int& iC : var_C){
          lampar_l_var.push_back( iA+10*iB+100*iC );
        }
      }
    }
  }

  if(DEBUG){
    printf(" lampar_l_var = { ");
    for( int& llv : lampar_l_var ){
      printf("%i ",llv);
    }
    printf("}\n");
  }

  //this is related to the feed-down from Xi1530
  //it should be negligible, it was assumed so in the PLB (flag 1)
  //I believe in the source paper it was taken into account (flag 0)
  //here we take both options as a variation, to prove that point
  std::vector<int> lampar_xi_var = {0,1};
  //taken from pS0 analysis, honestly the variations should be negligible
  std::vector<float> lampar_sig_var = {0.18,0.14,0.21};



  std::vector<float> cutoff_var = {700};//done
  std::vector<TString> MomSmearVar = {"pp13TeV_HM_DimiJun20"};//done
  //and an additional variation will be done on if pS0 is included as a feed or not
  enum BLTYPE { pol0s,pol1s,pol2s,pol3s,dpol2s,dpol3s,dpol4s,pol2e,pol3e,dpol2e,dpol3e,dpol4e,spl1 };
  std::vector<int> baseline_var;
  if(FitType==0){
    baseline_var.push_back(pol1s);
    baseline_var.push_back(pol0s);
  }
  else if(FitType==1){
    baseline_var.push_back(dpol3e);
  }

  //the region for which the DLM_Ck objects will be defined
  //we put some extra to get away from the edge effects of the smearing
  const double kMinTot = 0;
  const double kMaxTot = 528;
  const unsigned MomBinsTot = TMath::Nint(kMaxTot/12.);
  const unsigned NumMtBins = 6;
  if(imTbin>=NumMtBins || imTbin<-1){
    printf("\033[1;31mERROR:\033[0m Only %u mT bins are available!\n",NumMtBins);
    abort();
  }

  const double kMinCats = 0;
  const double kMaxCats = 336;
  const unsigned MomBinsCats = TMath::Nint(kMaxCats/12.);
  //const unsigned NumMomBins_feed = TMath::Nint(fitregion_var.back()/12.);

  const double kMinPlot = 0;
  const double kMaxPlot = fitregion_var.at(1);
  const unsigned MomBinsPlot = TMath::Nint(fitregion_var.at(1)/12.);

  int PotVar = pL_inter;
  TString PotType;
  //      -1: LO 600 only s-wave (as orginally in the paper)
  //      0: NLO-600 only s-wave (as orginally in the paper)
  //      potential ID as in the pL paper (s and d waves). In practice we will try out:
  //      11600: NLO-600
  //      131600: NLO-600 with reduced scattering length in the 3S1
  //      231600: NLO-600 with reduced scattering length in the 3S1, and charge symmetry breaking (CSB) in 1S0
  if(pL_inter==-1){
    PotType = "LO";
  }
  else if(pL_inter==0){
    PotType = "NLO";
  }
  else if(pL_inter==11600||pL_inter==131600||pL_inter==231600){
    PotType = "Chiral_Coupled_SPD";
  }
  else{
    printf("\033[1;31mERROR:\033[0m Unknown pL interaction (%i). Aborting!\n",pL_inter);
    abort();
  }

  if(DEBUG){
    printf(" PotType = %s\n",PotType.Data());
    printf("----------------------------------------\n");
  }

  //approximate, based on the output we kind of know we will get
  float ExpectedRadius;
  float FeedDownRadius;

  //values of the radii (as used in the source paper)
  std::vector<float> pS0_rad_srcPaper = { 1.473, 1.421, 1.368, 1.295, 1.220, 1.124 };
  //from the source paper
  std::vector<float> pp_reff = { 1.279, 1.231, 1.184, 1.110, 1.036, 0.948 };
  //averaged (done)
  std::vector<float> pL_reff = { 1.458, 1.398, 1.387, 1.302, 1.194, 1.018 };
  //with the correct sign (done)
  std::vector<float> pp_rcore_new = { 1.114, 1.066, 1.020, 0.946, 0.873, 0.786 };
  //from the source paper published in PLB
  std::vector<float> pp_rcore_srcPaper = { 1.209, 1.158, 1.108, 1.028, 0.950, 0.856 };

  //here we set out FeedDownRadius (will be used both for pS0 and pXi)
  if(imTbin>=0){
    if(NewAgeFlag%10==0){
      FeedDownRadius = pS0_rad_srcPaper.at(imTbin);
    }
    else if(NewAgeFlag%10==1){
      FeedDownRadius = pp_reff.at(imTbin);
    }
    else if(NewAgeFlag%10==2){
      FeedDownRadius = RanGen.Uniform((pp_rcore_new.at(imTbin)+pp_rcore_srcPaper.at(imTbin))*0.5,pL_reff.at(imTbin));
    }
  }
  else{
    //the 3rd mT bin corresponds roughly to the <mT> of the full data sample
    if(NewAgeFlag%10==0){
      FeedDownRadius = pS0_rad_srcPaper.at(3);
    }
    else if(NewAgeFlag%10==1){
      FeedDownRadius = pp_reff.at(3);
    }
    else if(NewAgeFlag%10==2){
      FeedDownRadius = RanGen.Uniform((pp_rcore_new.at(3)+pp_rcore_srcPaper.at(3))*0.5,pL_reff.at(3));
    }
  }

  DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder(InputFolder);

  //gROOT->cd();
//const int RndSeed, const int imTbin, const int pL_inter, const int AnalysisType, const int NewAgeFlag, const char* InputFolder, const char* OutputFolder
  TFile* fOutputFile = new TFile(OutputFolder+TString::Format("/fOut_rnd%i_mT%i_Pot%i_at%i_et%i_naf%i.root",
  RndSeed,imTbin,pL_inter,AnalysisType,ErrorType,NewAgeFlag),"recreate");

//WE NEED TO DELETE!!!!
  TH1F* hData = NULL;
  TH1F* hDataToFit = NULL;
  //GDATA->Set(TotMomBins);
  TGraph* gFit = new TGraph();
  gFit->SetName("gFit");
  //GFIT->Set(TotMomBins);
  TGraph* gBl = new TGraph();
  gBl->SetName("gBl");
  //GBL->Set(TotMomBins);
  TGraph* gFemto = new TGraph();
  gFemto->SetName("gFemto");

  //GFEMTO->Set(TotMomBins);
  int Baseline;
  float FitRegion;
  float CkCutOff;
  float CuspWeight;


  double* lam_pars_pl;
  lam_pars_pl = new double [5];
  double* lam_pars_pxi;
  lam_pars_pxi = new double [5];
  float lam_gen;
  float lam_pS0;
  float lam_pXim;
  float lam_pXi0;
  float lam_flt;
  float lam_mid;
  float lam_pXi_gen;
  float lam_pXi_pXi1530;
  float lam_pXi_flt;
  float lam_pS0_flt;
  float chi2_full;
  int ndf_full;
  int ndp_full;
  float nsig_full;
  float chi2_300;
  float chi2_204;
  float chi2_108;
  float radius;
  float raderr;

  fOutputFile->cd();
  TTree* pLTree = new TTree("pLTree","pLTree");
  pLTree->Branch("hData","TH1F",&hData,32000,0);//
  pLTree->Branch("hDataToFit","TH1F",&hDataToFit,32000,0);//
  pLTree->Branch("gFit","TGraph",&gFit,32000,0);//
  pLTree->Branch("gBl","TGraph",&gBl,32000,0);//
  pLTree->Branch("gFemto","TGraph",&gFemto,32000,0);//
  pLTree->Branch("Baseline",&Baseline,"Baseline/I");//
  pLTree->Branch("FitRegion",&FitRegion,"FitRegion/F");//
  pLTree->Branch("Bootstrap",&Bootstrap,"Bootstrap/I");//
  pLTree->Branch("SystVar",&SystVar,"SystVar/I");//
  pLTree->Branch("AnalysisType",&AnalysisType,"AnalysisType/I");//
  pLTree->Branch("SourceType","TString",&SourceType,8000,0);//
  pLTree->Branch("SourceVar",&SourceVar,"SourceVar/I");//
  pLTree->Branch("PotType","TString",&PotType,8000,0);//
  pLTree->Branch("PotVar",&PotVar,"PotVar/I");//
  pLTree->Branch("FeedDownType",&FeedDownType,"FeedDownType/I");//
  pLTree->Branch("FeedDownRadius",&FeedDownRadius,"FeedDownRadius/F");//
  pLTree->Branch("FitType",&FitType,"FitType/I");//
  pLTree->Branch("SmearStrategy",&SmearStrategy,"SmearStrategy/I");//
  pLTree->Branch("lam_gen",&lam_gen,"lam_gen/F");
  pLTree->Branch("lam_pS0",&lam_pS0,"lam_pS0/F");
  pLTree->Branch("lam_pXim",&lam_pXim,"lam_pXim/F");
  pLTree->Branch("lam_pXi0",&lam_pXi0,"lam_pXi0/F");
  pLTree->Branch("lam_flt",&lam_flt,"lam_flt/F");
  pLTree->Branch("lam_mid",&lam_mid,"lam_mid/F");
  pLTree->Branch("lam_pXi_gen",&lam_pXi_gen,"lam_pXi_gen/F");
  pLTree->Branch("lam_pXi_pXi1530",&lam_pXi_pXi1530,"lam_pXi_pXi1530/F");
  pLTree->Branch("lam_pXi_flt",&lam_pXi_flt,"lam_pXi_flt/F");
  pLTree->Branch("lam_pS0_flt",&lam_pS0_flt,"lam_pS0_flt/F");
  pLTree->Branch("CkCutOff",&CkCutOff,"CkCutOff/F");
  pLTree->Branch("CuspWeight",&CuspWeight,"CuspWeight/F");
  pLTree->Branch("chi2_full",&chi2_full,"chi2_full/F");
  pLTree->Branch("ndf_full",&ndf_full,"ndf_full/I");
  pLTree->Branch("ndp_full",&ndp_full,"ndp_full/I");
  pLTree->Branch("nsig_full",&nsig_full,"nsig_full/F");
  pLTree->Branch("chi2_300",&chi2_300,"chi2_300/F");
  pLTree->Branch("chi2_204",&chi2_204,"chi2_204/F");
  pLTree->Branch("chi2_108",&chi2_108,"chi2_108/F");
  pLTree->Branch("radius",&radius,"radius/F");//
  pLTree->Branch("raderr",&raderr,"raderr/F");//


  CATS AB_pXim;
  //same binning as pL, as we only use pXim as feed-down
  AB_pXim.SetMomBins(MomBinsCats,kMinCats,kMaxCats);
  if(FeedDownType==0||FeedDownType==1){
    //old, historic, used in the source paper
    AnalysisObject.SetUpCats_pXim(AB_pXim,"pXim_HALQCD1","Gauss");
  }
  else{
    //the one used in the Nature paper
    AnalysisObject.SetUpCats_pXim(AB_pXim,"pXim_HALQCDPaper2020","Gauss");
  }
  AB_pXim.SetAnaSource(0,FeedDownRadius);
  AB_pXim.SetNotifications(CATS::nWarning);
  AB_pXim.KillTheCat();
  DLM_Ck* Ck_pXim = new DLM_Ck(1,0,AB_pXim,MomBinsTot,kMinTot,kMaxTot);

  CATS AB_pXi0;
  //same binning as pL, as we only use pXim as feed-down
  AB_pXi0.SetMomBins(MomBinsCats,kMinCats,kMaxCats);
  AnalysisObject.SetUpCats_pXi0(AB_pXi0,"pXim_HALQCDPaper2020","Gauss");
  AB_pXi0.SetAnaSource(0,FeedDownRadius);
  AB_pXi0.SetNotifications(CATS::nWarning);
  AB_pXi0.KillTheCat();
  DLM_Ck* Ck_pXi0 = NULL;
  if(FitType==1) Ck_pXi0 = new DLM_Ck(1,0,AB_pXi0,MomBinsTot,kMinTot,kMaxTot);


  CATS AB_pXi1530;
  AB_pXi1530.SetMomBins(MomBinsCats,kMinCats,kMaxCats);
  AB_pXi1530.SetNotifications(CATS::nWarning);
  AnalysisObject.SetUpCats_pXim(AB_pXi1530,"pXim1530","Gauss");//McLevyNolan_Reso
  AB_pXi1530.SetAnaSource(0,FeedDownRadius);
  AB_pXi1530.SetNotifications(CATS::nWarning);
  AB_pXi1530.KillTheCat();
  DLM_Ck* Ck_pXim1530 = Ck_pXim1530 = new DLM_Ck(1,0,AB_pXi1530,MomBinsTot,kMinTot,kMaxTot);;




  CATS AB_pS0;
  DLM_Ck* Ck_pS0 = NULL;
  if(FeedDownType==1||FeedDownType==2){
    //the minus one is to to avoid going above 350 MeV, since we do not have the WF there
    AB_pS0.SetMomBins(MomBinsCats,kMinCats,kMaxCats);
    if(FeedDownType==1){
      AnalysisObject.SetUpCats_pS0(AB_pS0,"ESC16","Gauss");
    }
    else if(FeedDownType==2){
      AnalysisObject.SetUpCats_pS0(AB_pS0,"Chiral","Gauss");
    }
    AB_pS0.SetAnaSource(0,FeedDownRadius);
    AB_pS0.SetNotifications(CATS::nWarning);
    AB_pS0.KillTheCat();
    Ck_pS0 = new DLM_Ck(1,0,AB_pS0,MomBinsTot,kMinTot,kMaxTot);
  }
  else{
    if(FeedDownType==0){
      Ck_pS0 = new DLM_Ck(1,0,MomBinsTot,kMinTot,kMaxTot,Lednicky_gauss_Sigma0);
      Ck_pS0->SetSourcePar(0,FeedDownRadius);
    }
  }


  CATS AB_pL;
  AB_pL.SetMomBins(MomBinsCats,kMinCats,kMaxCats);
  AnalysisObject.SetUpCats_pL(AB_pL,PotType,SourceType,PotVar,SourceVar);
  if(SourceType=="Gauss"){
    AB_pL.SetAnaSource(0,pL_reff.at(imTbin));
  }
  else{
    AB_pL.SetAnaSource(0,pp_rcore_new.at(imTbin));
    AB_pL.SetAnaSource(1,2.0);
  }
  AB_pL.SetNotifications(CATS::nError);
  AB_pL.KillTheCat();
  DLM_Ck* Ck_pL = new DLM_Ck(NumSourcePars,0,AB_pL,MomBinsTot,kMinTot,kMaxTot);

  TH2F* hResolution_pL = AnalysisObject.GetResolutionMatrix("pp13TeV_HM_DimiJun20","pLambda");
  TH2F* hResidual_pp_pL = AnalysisObject.GetResidualMatrix("pp","pLambda");
  TH2F* hResidual_pL_pSigma0 = AnalysisObject.GetResidualMatrix("pLambda","pSigma0");
  TH2F* hResidual_pL_pXim = AnalysisObject.GetResidualMatrix("pLambda","pXim");
  TH2F* hResidual_pXi_pXi1530 = AnalysisObject.GetResidualMatrix("pXim","pXim1530");

  std::vector<float> cusp_var;
  //the case where we use the new chiral calculations
  if(PotType=="Chiral_Coupled_SPD"){
    cusp_var.push_back(0.40);
    cusp_var.push_back(0.27);
    cusp_var.push_back(0.33);
  }
  //for the old calcs we dont even have this peak
  else{
    cusp_var.push_back(0);
  }

  unsigned CompletedIters = 0;

  long long TotalTime;
  while(double(TotalTime)<Timeout){
    if(DEBUG){
      printf("Iter #%u\n",CompletedIters);
    }
    //printf("t %lld\n",TotalTime);
    int WhichCuspVar;
    int WhichCutOffVar;
    int WhichLambdaVar;
    int WhichProtonVar;
    int WhichXiVar;
    int WhichSigmaVar;
    int WhichBaselineVar;
    int WhichFitRegVar;
    int WhichData;
    if( (CompletedIters==0 && RndSeed==1) || !SystVar){
      WhichCuspVar = 0;
      WhichCutOffVar = 0;
      WhichProtonVar = 0;
      WhichLambdaVar = 0;
      WhichXiVar = 0;
      WhichSigmaVar = 0;
      WhichBaselineVar = 0;
      WhichFitRegVar = 0;
      WhichData = 0;
    }
    else{
      WhichCuspVar = RanGen.Integer(cusp_var.size());
      WhichCutOffVar = RanGen.Integer(cutoff_var.size());
      WhichProtonVar = RanGen.Integer(lampar_p_var.size());
      WhichLambdaVar = RanGen.Integer(lampar_l_var.size());
      WhichSigmaVar = RanGen.Integer(lampar_sig_var.size());
      WhichXiVar = RanGen.Integer(lampar_xi_var.size());
      WhichBaselineVar = RanGen.Integer(baseline_var.size());
      WhichFitRegVar = RanGen.Integer(fitregion_var.size());
      WhichData = RanGen.Integer(43);
    }

    if(DEBUG){
      printf("VAR: %i %i %i %i %i %i %i %i %i\n",
      WhichCuspVar,
      WhichCutOffVar,
      WhichProtonVar,
      WhichLambdaVar,
      WhichSigmaVar,
      WhichXiVar,
      WhichBaselineVar,
      WhichFitRegVar,
      WhichData);
    }

    gROOT->cd();
//what about mt integrated !!!!
    hData = AnalysisObject.GetAliceExpCorrFun("pp13TeV_HM_BernieSource","pLambda",TString::Format("%u",WhichData),0,0,imTbin);
    hDataToFit = (TH1F*)hData->Clone("hDataToFit");
    if(Bootstrap){
      for(unsigned uBin=0; uBin<MomBinsTot; uBin++){
        double BootVal = hData->GetBinContent(uBin+1);
        double BootErr = hData->GetBinError(uBin+1);
        //assuming E/N = sqrt(N)/N = 1/sqrt(N) = relErr; N = 1/relErr^2
        //the effective counts N represent the number of counts in the same event, assuming that
        //the error of the correlation stems ONLY from the same event
        double EffectiveCounts = pow(BootErr/BootVal,-2.);
        //we sample from poisson to avoid the chance of getting negative values
        BootVal = RanGen.Poisson(EffectiveCounts)*BootVal/EffectiveCounts;
        //BootVal = RanGen.Gaus(BootVal,BootErr);
        hDataToFit->SetBinContent(uBin+1,BootVal);
        hDataToFit->SetBinError(uBin+1,BootErr);
      }
    }

    //WhichProtonVar and WhichLambdaVar
    //if(imTbin>=0){
    //  AnalysisObject.SetUpLambdaPars_pL("pp13TeV_HM_BernieSource",lampar_p_var.at(WhichProtonVar),lampar_l_var.at(WhichLambdaVar),lam_pars_pl);
    //}
    //else{
      AnalysisObject.SetUpLambdaPars_pL("pp13TeV_HM_Dec19",lampar_p_var.at(WhichProtonVar),lampar_l_var.at(WhichLambdaVar),lam_pars_pl);
    //}
    AnalysisObject.SetUpLambdaPars_pXim("pp13TeV_HM_Dec19",lampar_p_var.at(WhichProtonVar),0,lam_pars_pxi);

    if(FitType==0){
      lam_gen = lam_pars_pl[0];
      lam_pS0 = lam_pars_pl[1];
      lam_pXim = lam_pars_pl[2];
      lam_pXi0 = 0;
      lam_flt = lam_pars_pl[3];
      lam_mid = lam_pars_pl[4];
    }
    else{
      lam_gen = lam_pars_pl[0];
      lam_pS0 = lam_pars_pl[1];
      lam_pXim = lam_pars_pl[2];
      lam_pXi0 = lam_pars_pl[2];
      lam_flt = lam_pars_pl[3]-lam_pars_pl[2];
      lam_mid = lam_pars_pl[4];
    }
    if(fabs(lam_gen+lam_pS0+lam_pXim+lam_pXi0+lam_flt+lam_mid-1)>1e-6){
      printf("The lambda pars dont sum up to 1 (%f)\n",lam_gen+lam_pS0+lam_pXim+lam_pXi0+lam_flt+lam_mid);
      abort();
    }
    lam_gen = 1.-lam_pS0-lam_pXim-lam_pXi0-lam_flt-lam_mid;
    lam_pXi_gen = lam_pars_pxi[0]/(1.-lam_pars_pxi[4]);
    //the Xi1530 can go to Lambda either -> Xi0 -> L or Xim -> L
    //the fractions to each are a bit different, but the effect on Lambda is always the same,
    //i.e. we assume that 50% of Xi1530 goes through each channel (even if it is 1:2 in truth)
    if(lampar_xi_var.at(WhichXiVar)==0){
      lam_pXi_pXi1530 = 0.5*(lam_pars_pxi[1]+lam_pars_pxi[2])/(1.-lam_pars_pxi[4]);
    }
    else{
      lam_pXi_pXi1530 = 0;
    }
    lam_pXi_flt = 1.-lam_pXi_gen-lam_pXi_pXi1530;

    lam_pS0_flt = lampar_sig_var.at(WhichSigmaVar);

    if(lam_pXi_pXi1530==0){
      delete Ck_pXim1530;
      Ck_pXim1530 = NULL;
    }

    if(DEBUG){
      printf(" lam_gen = %.2f%%\n",lam_gen*100.);
      printf(" lam_pS0 = %.2f%%\n",lam_pS0*100.);
      printf(" lam_pXim= %.2f%%\n",lam_pXim*100.);
      printf(" lam_pXi0= %.2f%%\n",lam_pXi0*100.);
      printf(" lam_flt = %.2f%%\n",lam_flt*100.);
      printf(" lam_mid = %.2f%%\n",lam_mid*100.);
      printf(" lam_pXi_gen = %.2f%%\n",lam_pXi_gen*100.);
      printf(" lam_pXi_pXi1530 = %.2f%%\n",lam_pXi_pXi1530*100.);
      printf(" lam_pXi_flt = %.2f%%\n",lam_pXi_flt*100.);
      printf(" lam_pS0_flt = %.2f%%\n",lam_pS0_flt*100.);
    }

    //WhichCuspVar
    CuspWeight = cusp_var.at(WhichCuspVar);
    if(CuspWeight){
        AB_pL.SetChannelWeight(7,1./4.*CuspWeight);//1S0 SN(s) -> LN(s)
        AB_pL.SetChannelWeight(8,3./4.*CuspWeight);//3S1 SN(s) -> LN(s)
        AB_pL.SetChannelWeight(10,3./4.*CuspWeight);//3S1 SN(d) -> LN(s)
        AB_pL.SetChannelWeight(13,3./20.*CuspWeight);//3D1 SN(d) -> LN(d)
        AB_pL.SetChannelWeight(15,3./20.*CuspWeight);//3D1 SN(s) -> LN(d)
    }
    AB_pL.KillTheCat();
    if(DEBUG){
      printf(" CuspWeight = %.0f%%\n",CuspWeight*100.);
    }

    //WhichCutOffVar and WhichFemtoRegVar
    CkCutOff = cutoff_var.at(WhichCutOffVar);
    FitRegion = fitregion_var.at(WhichFitRegVar);
    Ck_pL->SetCutOff(kMaxCats,CkCutOff);
    if(Ck_pS0) Ck_pS0->SetCutOff(kMaxCats,CkCutOff);
    if(Ck_pXim) Ck_pXim->SetCutOff(kMaxCats,CkCutOff);
    if(Ck_pXi0) Ck_pXi0->SetCutOff(kMaxCats,CkCutOff);
    if(Ck_pXim1530) Ck_pXim1530->SetCutOff(kMaxCats,CkCutOff);

    if(DEBUG){
      printf(" CkCutOff = %.0f MeV\n",CkCutOff);
      printf(" FitRegion = %.0f MeV\n",FitRegion);
      usleep(200e3);
    }

    TH1F* hPhaseSpace_pL=NULL;
    if(SmearStrategy==1){
      TList* list1_tmp;
      TList* list2_tmp;
      TString FileName = TString::Format("%s/ExpData/Bernie_Source/pLData/mTBin_%i/CFOutput_mT_pLVar%u_HM_%i.root",
                                          InputFolder,imTbin+1,WhichData,imTbin);
      TFile* inFile = new TFile(FileName,"read");
      //PARTICLES
      list1_tmp = (TList*)inFile->Get("PairDist");
      list2_tmp = (TList*)list1_tmp->FindObject("PairFixShifted");//there is also PairShifted ?? which one ??
      TH1F* hME_PP = (TH1F*)list2_tmp->FindObject(TString::Format("MEPart_mT_%i_FixShifted",imTbin));
      list1_tmp = (TList*)inFile->Get("AntiPairDist");
      list2_tmp = (TList*)list1_tmp->FindObject("PairFixShifted");
      TH1F* hME_APAP = (TH1F*)list2_tmp->FindObject(TString::Format("MEAntiPart_mT_%i_FixShifted",imTbin));
      gROOT->cd();
      hPhaseSpace_pL = (TH1F*)hME_PP->Clone("hPhaseSpace_pL");
      hPhaseSpace_pL->Add(hME_APAP);
      hPhaseSpace_pL->GetXaxis()->SetLimits(hPhaseSpace_pL->GetXaxis()->GetXmin()*1000.,hPhaseSpace_pL->GetXaxis()->GetXmax()*1000.);
      //printf("hPhaseSpace_pL %p\n",hPhaseSpace_pL);
      delete inFile;
    }

    //fOutputFile->cd();
    //if(hPhaseSpace_pL) hPhaseSpace_pL->Write();

    DLM_CkDecomposition* CkDec_pL = NULL;
    DLM_CkDecomposition* CkDec_pS0 = NULL;
    DLM_CkDecomposition* CkDec_pXim = NULL;
    DLM_CkDecomposition* CkDec_pXi0 = NULL;
    DLM_CkDecomposition* CkDec_pXim1530 = NULL;

    //pS0,pXi0,pXim,flat,misid
    CkDec_pL = new DLM_CkDecomposition("pLambda",5,*Ck_pL,hResolution_pL);
    if(Ck_pS0) CkDec_pS0 = new DLM_CkDecomposition("pSigma0",0,*Ck_pS0,NULL);
    //xi1530,flat
    if(Ck_pXim) CkDec_pXim = new DLM_CkDecomposition("pXim",2,*Ck_pXim,NULL);
    if(Ck_pXi0) CkDec_pXi0 = new DLM_CkDecomposition("pXi0",2,*Ck_pXi0,NULL);
    if(Ck_pXim1530) CkDec_pXim1530 = new DLM_CkDecomposition("pXim1530",0,*Ck_pXim1530,NULL);

    if(DEBUG){
      printf(" lam_pXi_gen = %.2f%%\n",lam_pXi_gen*100.);
      printf(" lam_pXi_pXi1530 = %.2f%%\n",lam_pXi_pXi1530*100.);
      printf(" lam_pXi_flt = %.2f%%\n",lam_pXi_flt*100.);
    }


    if(SmearStrategy==1) {CkDec_pL->AddPhaseSpace(hPhaseSpace_pL);}//momentum smear correction
    CkDec_pL->AddContribution(0,lam_pS0,DLM_CkDecomposition::cFeedDown,CkDec_pS0?CkDec_pS0:NULL,CkDec_pS0?hResidual_pL_pSigma0:NULL);
    if(SmearStrategy==1&&CkDec_pS0) {CkDec_pL->AddPhaseSpace(0,hPhaseSpace_pL);}//feed-down smear correction
    CkDec_pL->AddContribution(1,lam_pXim,DLM_CkDecomposition::cFeedDown,CkDec_pXim?CkDec_pXim:NULL,CkDec_pXim?hResidual_pL_pXim:NULL);
    if(SmearStrategy==1&&CkDec_pXim) {CkDec_pL->AddPhaseSpace(1,hPhaseSpace_pL);}//feed-down smear correction
    CkDec_pL->AddContribution(2,lam_pXi0,DLM_CkDecomposition::cFeedDown,CkDec_pXi0?CkDec_pXi0:NULL,CkDec_pXi0?hResidual_pL_pXim:NULL);
    if(SmearStrategy==1&&CkDec_pXi0) {CkDec_pL->AddPhaseSpace(2,hPhaseSpace_pL);}//feed-down smear correction

    CkDec_pL->AddContribution(3,lam_flt,DLM_CkDecomposition::cFeedDown);
    CkDec_pL->AddContribution(4,lam_mid,DLM_CkDecomposition::cFake);

    if(CkDec_pS0){
      CkDec_pS0->AddContribution(1,lam_pS0_flt,DLM_CkDecomposition::cFeedDown);
    }

    CkDec_pXim->AddContribution(0,lam_pXi_pXi1530,DLM_CkDecomposition::cFeedDown,CkDec_pXim1530?CkDec_pXim1530:NULL,CkDec_pXim1530?hResidual_pXi_pXi1530:NULL);
    if(SmearStrategy==1&&CkDec_pXim1530) {CkDec_pXim->AddPhaseSpace(0,hPhaseSpace_pL);}
    CkDec_pXim->AddContribution(1,lam_pXi_flt,DLM_CkDecomposition::cFeedDown);

    if(CkDec_pXi0){
      CkDec_pXi0->AddContribution(0,lam_pXi_pXi1530,DLM_CkDecomposition::cFeedDown,CkDec_pXim1530?CkDec_pXim1530:NULL,CkDec_pXim1530?hResidual_pXi_pXi1530:NULL);
      if(SmearStrategy==1&&CkDec_pXim1530) {CkDec_pXi0->AddPhaseSpace(0,hPhaseSpace_pL);}
      CkDec_pXi0->AddContribution(1,lam_pXi_flt,DLM_CkDecomposition::cFeedDown);
    }

    if(CkDec_pL) CkDec_pL->Update();
    if(CkDec_pS0) CkDec_pS0->Update();
    if(CkDec_pXim) CkDec_pXim->Update();
    if(CkDec_pXi0) CkDec_pXi0->Update();
    if(CkDec_pXim1530) CkDec_pXim1530->Update();

    SOURCE_FIT = CkDec_pL;
    TF1* fData = new TF1("fData",Fit_BernieSource,0,FitRegion,7);
    TF1* fBl = new TF1("fBl",Baseline_BernieSource,0,FitRegion,5);
    fData->SetParameter(0,SourceType=="Gauss"?pL_reff.at(imTbin):pp_rcore_srcPaper.at(imTbin));
    fData->SetParLimits(0,fData->GetParameter(0)*0.5,fData->GetParameter(0)*2.0);
    fData->FixParameter(1,0);

    //BL
    fData->SetParameter(2,1);

    Baseline = baseline_var.at(WhichBaselineVar);
    //pol0s,pol1s,dpol3e
    if(Baseline==pol1s){
      fData->SetParameter(3,0);
      fData->SetParLimits(3,-1e-2,1e-2);
      fData->FixParameter(4,0);
      fData->FixParameter(5,0);
      fData->FixParameter(6,0);
    }
    else if(Baseline==dpol3e){
      fData->FixParameter(3,0);
      fData->SetParameter(4,100);
      fData->SetParLimits(4,0,400);
      fData->SetParameter(5,0);
      fData->SetParLimits(5,-1e-6,1e-6);
      fData->FixParameter(6,0);
    }
    else{
      fData->FixParameter(3,0);
      fData->FixParameter(4,0);
      fData->FixParameter(5,0);
      fData->FixParameter(6,0);
    }

    //printf("BL=%i FIT=%.0f PS0=%i SS=%i LAM=%i SMR=%s PL=%i SRC=%s lam_pL=%.1f XXX=%.1f\n",
    //        varBL,varFit,ipS0,SmearStrategy,varLam,varSmear.Data(),varPL,SourceVar.Data(),lambda_pL[0]*100.,lambda_pL[1]*100.);
    if(DEBUG) hDataToFit->Fit(fData,"S, N, R, M");
    else hDataToFit->Fit(fData,"Q, S, N, R, M");
    fData->SetNpx(3*MomBinsTot);

    for(unsigned uPar=0; uPar<5; uPar++){
      fBl->FixParameter(uPar,fData->GetParameter(uPar+2));
    }

    for(unsigned uBin=0; uBin<MomBinsPlot; uBin++){
      double Momentum = Ck_pL->GetBinCenter(0,uBin);
      bool MomQA = fabs(hData->GetBinCenter(uBin+1)-Momentum)<1e-6;
      if(!MomQA){
        printf("WARNING: Something is fishy with the binning of the data (%.2f vs %.2f)!!! Thank you GentleFemto!\n",Momentum,hData->GetBinCenter(uBin+1));
      }
      gFit->SetPoint(uBin,Momentum,fData->Eval(Momentum));
      gBl->SetPoint(uBin,Momentum,fBl->Eval(Momentum));
      gFemto->SetPoint(uBin,Momentum,fData->Eval(Momentum)/fBl->Eval(Momentum));
    }


    radius = fData->GetParameter(0);
    raderr = fData->GetParError(0);
    chi2_full = fData->GetChisquare();
    ndf_full = fData->GetNDF();
    ndp_full = fData->GetNumberFitPoints();
    nsig_full = GetNsigma(chi2_full, ndf_full);

    if(DEBUG){
      fOutputFile->cd();
      hData->Write();
      hDataToFit->Write();
      fData->Write();
      fBl->Write();
      gFit->Write();
      gBl->Write();
      gFemto->Write();
    }

    fOutputFile->cd();
    pLTree->Fill();

    TotalTime = JaimeTimer.Stop()/1000000.;
    CompletedIters++;
    if(hData) delete hData;
    if(hDataToFit) delete hDataToFit;
    if(fData) delete fData;
    if(fBl) delete fBl;
    if(CkDec_pL) delete CkDec_pL;
    if(CkDec_pS0) delete CkDec_pS0;
    if(CkDec_pXim) delete CkDec_pXim;
    if(CkDec_pXi0) delete CkDec_pXi0;
    if(CkDec_pXim1530) delete CkDec_pXim1530;
  }//while

  fOutputFile->cd();
  pLTree->Write();

  delete hResolution_pL;
  delete hResidual_pp_pL;
  delete hResidual_pL_pSigma0;
  delete hResidual_pL_pXim;
  delete hResidual_pXi_pXi1530;
  delete [] lam_pars_pl;
  delete [] lam_pars_pxi;
  delete pLTree;
  delete gFit;
  delete gBl;
  delete gFemto;
  if(Ck_pL) delete Ck_pL;
  if(Ck_pS0) delete Ck_pS0;
  if(Ck_pXim) delete Ck_pXim;
  if(Ck_pXi0) delete Ck_pXi0;
  if(Ck_pXim1530) delete Ck_pXim1530;
  delete fOutputFile;
}

void SourcePaper_Published(TString OutputFolder){

  const unsigned NumMtBins_pp = 7;
  double* BinRange_pp = new double[NumMtBins_pp+1];
  double* BinCenter_pp = new double[NumMtBins_pp];
  BinCenter_pp[0] = 1.1077;
  BinCenter_pp[1] = 1.1683;
  BinCenter_pp[2] = 1.2284;
  BinCenter_pp[3] = 1.3156;
  BinCenter_pp[4] = 1.4628;
  BinCenter_pp[5] = 1.6872;
  BinCenter_pp[6] = 2.2116;
  BinRange_pp[0] = BinCenter_pp[0]-BinCenter_pp[1]+BinCenter_pp[0];
  for(int i=1; i<=6; i++) BinRange_pp[i] = (BinCenter_pp[i-1]+BinCenter_pp[i])*0.5;
  BinRange_pp[7] = BinCenter_pp[6]+BinCenter_pp[6]-BinCenter_pp[5];

  const unsigned NumMtBins_pL = 6;
  double* BinRange_pL = new double[NumMtBins_pL+1];
  double* BinCenter_pL = new double[NumMtBins_pL];
  BinCenter_pL[0] = 1.2124;
  BinCenter_pL[1] = 1.2896;
  BinCenter_pL[2] = 1.376;
  BinCenter_pL[3] = 1.5407;
  BinCenter_pL[4] = 1.756;
  BinCenter_pL[5] = 2.2594;
  BinRange_pL[0] = BinCenter_pL[0]-BinCenter_pL[1]+BinCenter_pL[0];
  for(int i=1; i<=5; i++) BinRange_pL[i] = (BinCenter_pL[i-1]+BinCenter_pL[i])*0.5;
  BinRange_pL[6] = BinCenter_pL[5]+BinCenter_pL[5]-BinCenter_pL[4];


  DLM_Histo<float> dlm_rcore_pp;
  dlm_rcore_pp.SetUp(1);
  dlm_rcore_pp.SetUp(0,NumMtBins_pp,BinRange_pp,BinCenter_pp);
  dlm_rcore_pp.Initialize();
  for(unsigned uBin=0; uBin<7; uBin++){
    double val,stat,syst;
    switch (uBin) {
      case 0:   val = 1.3064;
                syst = 0.027805;
                stat = 0.0085539;
                break;
      case 1:   val = 1.2316;
                syst = 0.022773;
                stat = 0.010501;
                break;
      case 2:   val = 1.2006;
                syst = 0.022552;
                stat = 0.014732;
                break;
      case 3:   val = 1.1402;
                syst = 0.025014;
                stat = 0.011086;
                break;
      case 4:   val = 1.0628;
                syst = 0.025221;
                stat = 0.010027;
                break;
      case 5:   val = 0.96238;
                syst = 0.025234;
                stat = 0.0098228;
                break;
      case 6:   val = 0.86503;
                syst = 0.020369;
                stat = 0.010841;
                break;
      default:  break;
    }
    dlm_rcore_pp.SetBinContent(uBin,val);
    dlm_rcore_pp.SetBinError(uBin,sqrt(syst*syst+stat*stat));
  }

  DLM_Histo<float> dlm_rcore_pL;
  dlm_rcore_pL.SetUp(1);
  dlm_rcore_pL.SetUp(0,NumMtBins_pL,BinRange_pL,BinCenter_pL);
  dlm_rcore_pL.Initialize();
  DLM_Histo<float> dlm_rcore_pL_LO;
  dlm_rcore_pL_LO.SetUp(1);
  dlm_rcore_pL_LO.SetUp(0,NumMtBins_pL,BinRange_pL,BinCenter_pL);
  dlm_rcore_pL_LO.Initialize();
  DLM_Histo<float> dlm_rcore_pL_NLO;
  dlm_rcore_pL_NLO.SetUp(1);
  dlm_rcore_pL_NLO.SetUp(0,NumMtBins_pL,BinRange_pL,BinCenter_pL);
  dlm_rcore_pL_NLO.Initialize();
  for(unsigned uBin=0; uBin<6; uBin++){
    double nlo_val,nlo_stat,nlo_syst;
    double lo_val,lo_stat,lo_syst;
    switch (uBin) {
      case 0:   nlo_val = 1.2634;
                nlo_syst = 0.0876;
                nlo_stat = 0.041731;
                lo_val = 1.1816;
                lo_syst = 0.060069;
                lo_stat = 0.042102;
                break;
      case 1:   nlo_val = 1.1926;
                nlo_syst = 0.10035;
                nlo_stat = 0.050217;
                lo_val = 1.1381;
                lo_syst = 0.057701;
                lo_stat = 0.048758;
                break;
      case 2:   nlo_val = 1.1804;
                nlo_syst = 0.065382;
                nlo_stat = 0.042162;
                lo_val = 1.1279;
                lo_syst = 0.070412;
                lo_stat = 0.041525;
                break;
      case 3:   nlo_val = 1.0902;
                nlo_syst = 0.092647;
                nlo_stat = 0.035456;
                lo_val = 1.0516;
                lo_syst = 0.055087;
                lo_stat = 0.033123;
                break;
      case 4:   nlo_val = 0.9717;
                nlo_syst = 0.0825;
                nlo_stat = 0.041914;
                lo_val = 0.94642;
                lo_syst = 0.055;
                lo_stat = 0.038272;
                break;
      case 5:   nlo_val = 0.73832;
                nlo_syst = 0.18077;
                nlo_stat = 0.03668;
                lo_val = 0.82714;
                lo_syst = 0.050213;
                lo_stat = 0.033311;
                break;
      default:  break;
    }
    double avg_val = 0.5*(nlo_val+lo_val);
    double avg_stat = 0.5*(nlo_stat+lo_stat);
    double avg_syst = 0.5*(nlo_syst+lo_syst);
    double diff = nlo_val-lo_val;
    double tot_err = sqrt(avg_stat*avg_stat+avg_syst*avg_syst+diff*diff);
    dlm_rcore_pL.SetBinContent(uBin,avg_val);
    dlm_rcore_pL.SetBinError(uBin,tot_err);

    dlm_rcore_pL_LO.SetBinContent(uBin,lo_val);
    dlm_rcore_pL_LO.SetBinError(uBin,sqrt(lo_syst*lo_syst+lo_stat*lo_stat));

    dlm_rcore_pL_NLO.SetBinContent(uBin,nlo_val);
    dlm_rcore_pL_NLO.SetBinError(uBin,sqrt(nlo_syst*nlo_syst+nlo_stat*nlo_stat));
  }
  TGraphErrors g_rcore_pp;
  g_rcore_pp.SetName("g_rcore_pp");
  g_rcore_pp.SetFillColorAlpha(kBlue+1,0.3);
  g_rcore_pp.SetLineColor(kBlue+1);
  g_rcore_pp.SetLineWidth(3);
  for(unsigned uBin=0; uBin<dlm_rcore_pp.GetNbins(); uBin++){
    g_rcore_pp.SetPoint(uBin,dlm_rcore_pp.GetBinCenter(0,uBin),dlm_rcore_pp.GetBinContent(uBin));
    g_rcore_pp.SetPointError(uBin,0,dlm_rcore_pp.GetBinError(uBin));
  }
  TGraphErrors g_rcore_pL;
  g_rcore_pL.SetName("g_rcore_pL");
  g_rcore_pL.SetFillColorAlpha(kRed+1,0.3);
  g_rcore_pL.SetLineColor(kRed+1);
  g_rcore_pL.SetLineWidth(3);
  for(unsigned uBin=0; uBin<dlm_rcore_pL.GetNbins(); uBin++){
    g_rcore_pL.SetPoint(uBin,dlm_rcore_pL.GetBinCenter(0,uBin),dlm_rcore_pL.GetBinContent(uBin));
    g_rcore_pL.SetPointError(uBin,0,dlm_rcore_pL.GetBinError(uBin));
  }
  TGraphErrors g_rcore_pL_LO;
  g_rcore_pL_LO.SetName("g_rcore_pL_LO");
  g_rcore_pL_LO.SetFillColorAlpha(kRed+1,0.3);
  g_rcore_pL_LO.SetLineColor(kRed+1);
  g_rcore_pL_LO.SetLineWidth(3);
  for(unsigned uBin=0; uBin<dlm_rcore_pL_LO.GetNbins(); uBin++){
    g_rcore_pL_LO.SetPoint(uBin,dlm_rcore_pL_LO.GetBinCenter(0,uBin),dlm_rcore_pL_LO.GetBinContent(uBin));
    g_rcore_pL_LO.SetPointError(uBin,0,dlm_rcore_pL_LO.GetBinError(uBin));
  }
  TGraphErrors g_rcore_pL_NLO;
  g_rcore_pL_NLO.SetName("g_rcore_pL_NLO");
  g_rcore_pL_NLO.SetFillColorAlpha(kRed+1,0.3);
  g_rcore_pL_NLO.SetLineColor(kRed+1);
  g_rcore_pL_NLO.SetLineWidth(3);
  for(unsigned uBin=0; uBin<dlm_rcore_pL_NLO.GetNbins(); uBin++){
    g_rcore_pL_NLO.SetPoint(uBin,dlm_rcore_pL_NLO.GetBinCenter(0,uBin),dlm_rcore_pL_NLO.GetBinContent(uBin));
    g_rcore_pL_NLO.SetPointError(uBin,0,dlm_rcore_pL_NLO.GetBinError(uBin));
  }





  DLM_Histo<float> dlm_reff_pp;
  dlm_reff_pp.SetUp(1);
  dlm_reff_pp.SetUp(0,NumMtBins_pp,BinRange_pp,BinCenter_pp);
  dlm_reff_pp.Initialize();

  for(unsigned uBin=0; uBin<7; uBin++){
    double val,stat,syst;
    switch (uBin) {
      case 0:   val = 1.3717;
                syst = 0.025241;
                stat = 0.0079898;
                break;
      case 1:   val = 1.3005;
                syst = 0.022749;
                stat = 0.0098079;
                break;
      case 2:   val = 1.2709;
                syst = 0.022641;
                stat = 0.013852;
                break;
      case 3:   val = 1.2138;
                syst = 0.022693;
                stat = 0.010282;
                break;
      case 4:   val = 1.1423;
                syst = 0.027803;
                stat = 0.0094262;
                break;
      case 5:   val = 1.0482;
                syst = 0.022569;
                stat = 0.0090911;
                break;
      case 6:   val = 0.95678;
                syst = 0.017544;
                stat = 0.0099002;
                break;
      default:  break;
    }
    dlm_reff_pp.SetBinContent(uBin,val);
    dlm_reff_pp.SetBinError(uBin,sqrt(syst*syst+stat*stat));
  }

  DLM_Histo<float> dlm_reff_pL;
  dlm_reff_pL.SetUp(1);
  dlm_reff_pL.SetUp(0,NumMtBins_pL,BinRange_pL,BinCenter_pL);
  dlm_reff_pL.Initialize();
  DLM_Histo<float> dlm_reff_pL_LO;
  dlm_reff_pL_LO.SetUp(1);
  dlm_reff_pL_LO.SetUp(0,NumMtBins_pL,BinRange_pL,BinCenter_pL);
  dlm_reff_pL_LO.Initialize();
  DLM_Histo<float> dlm_reff_pL_NLO;
  dlm_reff_pL_NLO.SetUp(1);
  dlm_reff_pL_NLO.SetUp(0,NumMtBins_pL,BinRange_pL,BinCenter_pL);
  dlm_reff_pL_NLO.Initialize();
  for(unsigned uBin=0; uBin<6; uBin++){
    double nlo_val,nlo_stat,nlo_syst;
    double lo_val,lo_stat,lo_syst;
    switch (uBin) {
      case 0:   nlo_val = 1.502;
                nlo_syst = 0.085193;
                nlo_stat = 0.043523;
                lo_val = 1.4147;
                lo_syst = 0.060136;
                lo_stat = 0.042372;
                break;
      case 1:   nlo_val = 1.428;
                nlo_syst = 0.098118;
                nlo_stat = 0.051097;
                lo_val = 1.3687;
                lo_syst = 0.052833;
                lo_stat = 0.049813;
                break;
      case 2:   nlo_val = 1.4163;
                nlo_syst = 0.065325;
                nlo_stat = 0.043093;
                lo_val = 1.3582;
                lo_syst = 0.077887;
                lo_stat = 0.041928;
                break;
      case 3:   nlo_val = 1.3256;
                nlo_syst = 0.095299;
                nlo_stat = 0.035395;
                lo_val = 1.2788;
                lo_syst = 0.052665;
                lo_stat = 0.034165;
                break;
      case 4:   nlo_val = 1.2145;
                nlo_syst = 0.090394;
                nlo_stat = 0.040756;
                lo_val = 1.1734;
                lo_syst = 0.057752;
                lo_stat = 0.038961;
                break;
      case 5:   nlo_val = 0.98815;
                nlo_syst = 0.0825;
                nlo_stat = 0.051054;
                lo_val = 1.0473;
                lo_syst = 0.06;
                lo_stat = 0.03544;
                break;
      default:  break;
    }
    double avg_val = 0.5*(nlo_val+lo_val);
    double avg_stat = 0.5*(nlo_stat+lo_stat);
    double avg_syst = 0.5*(nlo_syst+lo_syst);
    double diff = nlo_val-lo_val;
    double tot_err = sqrt(avg_stat*avg_stat+avg_syst*avg_syst+diff*diff);
    dlm_reff_pL.SetBinContent(uBin,avg_val);
    dlm_reff_pL.SetBinError(uBin,tot_err);

    dlm_reff_pL_LO.SetBinContent(uBin,lo_val);
    dlm_reff_pL_LO.SetBinError(uBin,sqrt(lo_syst*lo_syst+lo_stat*lo_stat));

    dlm_reff_pL_NLO.SetBinContent(uBin,nlo_val);
    dlm_reff_pL_NLO.SetBinError(uBin,sqrt(nlo_syst*nlo_syst+nlo_stat*nlo_stat));
  }
  TGraphErrors g_reff_pp;
  g_reff_pp.SetName("g_reff_pp");
  g_reff_pp.SetFillColorAlpha(kBlue+1,0.3);
  g_reff_pp.SetLineColor(kBlue+1);
  g_reff_pp.SetLineWidth(3);
  for(unsigned uBin=0; uBin<dlm_reff_pp.GetNbins(); uBin++){
    g_reff_pp.SetPoint(uBin,dlm_reff_pp.GetBinCenter(0,uBin),dlm_reff_pp.GetBinContent(uBin));
    g_reff_pp.SetPointError(uBin,0,dlm_reff_pp.GetBinError(uBin));
  }
  TGraphErrors g_reff_pL;
  g_reff_pL.SetName("g_reff_pL");
  g_reff_pL.SetFillColorAlpha(kRed+1,0.3);
  g_reff_pL.SetLineColor(kRed+1);
  g_reff_pL.SetLineWidth(3);
  for(unsigned uBin=0; uBin<dlm_reff_pL.GetNbins(); uBin++){
    g_reff_pL.SetPoint(uBin,dlm_reff_pL.GetBinCenter(0,uBin),dlm_reff_pL.GetBinContent(uBin));
    g_reff_pL.SetPointError(uBin,0,dlm_reff_pL.GetBinError(uBin));
  }
  TGraphErrors g_reff_pL_LO;
  g_reff_pL_LO.SetName("g_reff_pL_LO");
  g_reff_pL_LO.SetFillColorAlpha(kRed+1,0.3);
  g_reff_pL_LO.SetLineColor(kRed+1);
  g_reff_pL_LO.SetLineWidth(3);
  for(unsigned uBin=0; uBin<dlm_reff_pL_LO.GetNbins(); uBin++){
    g_reff_pL_LO.SetPoint(uBin,dlm_reff_pL_LO.GetBinCenter(0,uBin),dlm_reff_pL_LO.GetBinContent(uBin));
    g_reff_pL_LO.SetPointError(uBin,0,dlm_reff_pL_LO.GetBinError(uBin));
  }
  TGraphErrors g_reff_pL_NLO;
  g_reff_pL_NLO.SetName("g_reff_pL_NLO");
  g_reff_pL_NLO.SetFillColorAlpha(kRed+1,0.3);
  g_reff_pL_NLO.SetLineColor(kRed+1);
  g_reff_pL_NLO.SetLineWidth(3);
  for(unsigned uBin=0; uBin<dlm_reff_pL_NLO.GetNbins(); uBin++){
    g_reff_pL_NLO.SetPoint(uBin,dlm_reff_pL_NLO.GetBinCenter(0,uBin),dlm_reff_pL_NLO.GetBinContent(uBin));
    g_reff_pL_NLO.SetPointError(uBin,0,dlm_reff_pL_NLO.GetBinError(uBin));
  }

  TFile fOutputFile(OutputFolder+TString::Format("/SourcePaper_Published.root"),"recreate");
  g_rcore_pp.Write();
  g_rcore_pL.Write();
  g_rcore_pL_LO.Write();
  g_rcore_pL_NLO.Write();
  g_reff_pp.Write();
  g_reff_pL.Write();
  g_reff_pL_LO.Write();
  g_reff_pL_NLO.Write();


  delete [] BinRange_pp;
  delete [] BinCenter_pp;
  delete [] BinRange_pL;
  delete [] BinCenter_pL;

}




void TestReadTTree(){
  TString InputFolder = "/home/dimihayl/Software/LocalFemto/Output/SourceStudies/SourcePaper_pp/";
  TString InputFile = "TEST2_Gauss_SS1_D0_mT0.root";

  //gROOT->cd();
  TFile* inFile = new TFile(InputFolder+InputFile,"read");
  //TTree* ppTree = new TTree("ppTree","ppTree");
  //delete ppTree;
  TTree* ppTree = (TTree*)inFile->Get("ppTree");
  unsigned NumEntries = ppTree->GetEntries();

  TGraph* gData;
  TGraph* gFit;
  TGraph* gBl;
  TGraph* gFemto;
  Float_t nsig;
  Float_t rad;
  Float_t FemtoRegion;
  Int_t BaselineVar;

  //ppTree->SetBranchAddress("gData",&gData);
  //ppTree->SetBranchAddress("gFit",&gFit);
  ppTree->SetBranchAddress("gBl",&gBl);
  ppTree->SetBranchAddress("gFemto",&gFemto);
  ppTree->SetBranchAddress("nsig",&nsig);
  ppTree->SetBranchAddress("rad",&rad);
  ppTree->SetBranchAddress("FemtoRegion",&FemtoRegion);
  ppTree->SetBranchAddress("BaselineVar",&BaselineVar);

//  TGraphErrors* geData = NULL;
//  TGraphErrors* geFit = NULL;
  std::vector<double> meanData;
  std::vector<double> stdvData;
  std::vector<unsigned> nentData;

  printf("NumEntries = %u\n\n",NumEntries);

  for(unsigned uEntry=0; uEntry<NumEntries; uEntry++){
    printf("uEntry = %u\n",uEntry);
    ppTree->GetEntry(uEntry);
    printf(" -->\n");
    printf(" FemtoRegion = %f\n",FemtoRegion);
    printf(" N = %i, %i, %i, %i\n",gData->GetN(),gFit->GetN(),gBl->GetN(),gFemto->GetN());
    printf(" r = %f; nsig=%f\n",rad,nsig);
  }

}


//random emission of single particles in the LAB frame
//we combine them in pairs and look at all the angular distribtutions
//relevant for the source.
//the emission follows boltzmann for the pT, uniform eta, uniform phi
//to get the spacial corrdinate, we assume some fixed mass of the particles
//and propagate them assuming they were created at 0,0. The propagation is for a fixed amount of time
void SinglePart_RandomLab(){
  const double Mass = 1000;//MeV
  //if its minus sign, the time is not fixed, rather proportional to abs(Time)/p
  const double Time = -5.0;//fm/c
  const double Temp = 700;//MeV, used for the boltzmann
  const double SP_r0 = 0.4;
  //probability to fix beta
  const double fixed_beta_T = 0;

  const unsigned NumEvents = 1024*32;
  const unsigned Multiplicity = 16;

  const TString OutputFolder = TString::Format("%s/SourceStudies/SinglePart_RandomLab/",GetFemtoOutputFolder());

  DLM_Random rangen(11);

  const unsigned NumMomBins = 1024;
  const double MomMin = 0;
  const double MomMax = 1024*8;
  const unsigned NumEtaBins = 1024;
  const double EtaMin = -2;
  const double EtaMax = 2;
  const unsigned NumPhiBins = 1024;
  const double PhiMin = -Pi;
  const double PhiMax = Pi;
  const unsigned NumThetaBins = 1024;
  const double ThetaMin = 0;
  const double ThetaMax = 1.*Pi;
  const unsigned NumRadBins = 4096;
  const double RadMin = 0;
  const double RadMax = 16;
  const unsigned NumMtBins = 1024;
  const double MtMin = 0;
  const double MtMax = 1024*4;
  const unsigned NumBetaBins = 1024;
  const double betaMin = 0;
  const double betaMax = 2;

  const double FemtoRegion = 200;
  const double EtaRegion = 0.8;
  const double MinMom = 500;

  //SP = SingleParticle
  TH1F* SP_p_LAB = new TH1F("SP_p_LAB","SP_p_LAB",NumMomBins,MomMin,MomMax);
  TH1F* SP_pT_LAB = new TH1F("SP_pT_LAB","SP_pT_LAB",NumMomBins,MomMin,MomMax);
  TH1F* SP_eta_LAB = new TH1F("SP_eta_LAB","SP_eta_LAB",NumEtaBins,EtaMin,EtaMax);
  TH1F* SP_phi_LAB = new TH1F("SP_phi_LAB","SP_phi_LAB",NumPhiBins,PhiMin,PhiMax);
  TH1F* SP_theta_LAB = new TH1F("SP_theta_LAB","SP_theta_LAB",NumThetaBins,ThetaMin,ThetaMax);
  TH1F* SP_costheta_LAB = new TH1F("SP_costheta_LAB","SP_costheta_LAB",NumThetaBins,-1,1);

  TH1F* SP_Sr_LAB = new TH1F("SP_Sr_LAB","SP_Sr_LAB",NumRadBins,RadMin,RadMax);
  TH2F* Pair_mT_Sr_LAB = new TH2F("Pair_mT_Sr_LAB","Pair_mT_Sr_LAB",NumMomBins,MomMin,MomMax,NumRadBins,RadMin,RadMax);
  TH2F* Pair_mT_Sr_CM = new TH2F("Pair_mT_Sr_CM","Pair_mT_Sr_CM",NumMomBins,MomMin,MomMax,NumRadBins,RadMin,RadMax);

  //only for femto pairs
  TH2F* FemPair_mT_Sr_LAB = new TH2F("FemPair_mT_Sr_LAB","FemPair_mT_Sr_LAB",NumMtBins,MtMin,MtMax,NumRadBins,RadMin,RadMax);
  TH2F* FemPair_mT_Sr_CM = new TH2F("FemPair_mT_Sr_CM","FemPair_mT_Sr_CM",NumMtBins,MtMin,MtMax,NumRadBins,RadMin,RadMax);
  TH1F* FemPair_Sr_LAB = new TH1F("FemPair_Sr_LAB","FemPair_Sr_LAB",NumRadBins,RadMin,RadMax);
  TH1F* FemPair_Sr_CM = new TH1F("FemPair_Sr_CM","FemPair_Sr_CM",NumRadBins,RadMin,RadMax);
  TH1F* FemPair_betaT_LAB = new TH1F("FemPair_betaT_LAB","FemPair_betaT_LAB",NumBetaBins,betaMin,betaMax);

  //TH1F* FemPair_Angle_kr_LAB = new TH1F("FemPair_Angle_kr_LAB","FemPair_Angle_kr_LAB",NumThetaBins,ThetaMin,ThetaMax);
  TH1F* FemPair_Angle_kr_CM = new TH1F("FemPair_Angle_kr_CM","FemPair_Angle_kr_CM",NumThetaBins,ThetaMin,ThetaMax);


  double p_tot;
  double p_T;
  double eta;
  double phi;
  double theta;
  double beta;
  double gamma;
  for(unsigned uEvent=0; uEvent<NumEvents; uEvent++){
    CatsEvent* KittyEvent_LAB = new CatsEvent(1,1);
    CatsEvent* KittyEvent_CM = new CatsEvent(1,1);
    double beta_T=0;
    bool FIX_BETA = rangen.Uniform(0,1)<fixed_beta_T;
    for(unsigned uMult=0; uMult<Multiplicity; uMult++){
      p_T = rangen.GaussR(3,0,Temp);
      //common tranverse flow velocity
      if(FIX_BETA){
        //printf("Hello\n");
        //random first particle
        if(uMult==0){
          //theta = acos(rangen.Uniform(-1,1));
          //this will be our leading particle
          theta = acos(rangen.Uniform(-0.7,0.7));
          eta = -log(tan(0.5*theta));
          phi = rangen.Uniform(0.,2.*Pi);
          //fix the beta_T from now on
          p_tot = p_T/sin(theta);
          gamma = sqrt(Mass*Mass+p_tot*p_tot)/Mass;
          beta_T = p_T/sqrt(Mass*Mass+p_tot*p_tot);
        }
        else{
          gamma = p_T/(Mass*beta_T);
          gamma = p_T/(Mass*beta_T);
          beta = sqrt(1.-1./gamma/gamma);
          p_tot = Mass*beta*gamma;
          theta = asin(p_T/p_tot);

          //while(gamma<1 || fabs(cos(theta))>0.7){
          //while(gamma<1){
          while(gamma<1 || fabs(cos(theta))>0.7 || p_T/p_tot<0.7){
            p_T = rangen.GaussR(3,0,Temp);
            gamma = p_T/(Mass*beta_T);
            beta = sqrt(1.-1./gamma/gamma);
            p_tot = Mass*beta*gamma;
            //printf("p_tot = %.2f\n",p_tot);
            theta = asin(p_T/p_tot);
          }

          if(rangen.Uniform(0.,1.)<0.5) theta = Pi-theta;
          eta = -log(tan(0.5*theta));
          phi = rangen.Uniform(0.,2.*Pi);
        }
      }
      else{
        //eta = rangen.Uniform(-0.8,0.8);
        //if we want it flat in all spacial directions
        theta = acos(rangen.Uniform(-1,1));
        eta = -log(tan(0.5*theta));
        phi = rangen.Uniform(0.,2.*Pi);
        p_tot = p_T/sin(theta);
      }

      double time;
      if(Time>0) time=Time;
      else{
        time = fabs(Time*FmToNu/p_tot);
      }
      printf("Time = %.2f\n",Time);
      printf("FmToNu = %.3f\n",FmToNu);
      printf("p_tot = %.2f\n",p_tot);
      printf("time = %.2f\n",time);
      printf("calc = %.2f\n\n",Time*FmToNu/p_tot);
      usleep(100e3);
      TLorentzVector TLM_Part;
      TLM_Part.SetPtEtaPhiM(p_T,eta,phi,Mass);
      //printf(" TLM_Part.P()=%.2f\n",TLM_Part.P());
      //this will come out in fm, it is the distance traveled
      double BGT = TLM_Part.Beta()*TLM_Part.Gamma()*time;
      FemPair_betaT_LAB->Fill(TLM_Part.Pt()/(TLM_Part.Gamma()*TLM_Part.M()));

      //we simply set the radius vector pointing in the direction of travel
      //at the estimated distance of travel assuming constant time for the formation
      TLorentzVector TLS_Part;
      TLS_Part = TLM_Part;
      TLS_Part.SetT(time);
      TLS_Part.SetRho(TLM_Part.Beta()*time+rangen.GaussR(3,0,SP_r0));
      //TLS_Part.SetRho(BGT);
      //TLS_Part.SetRho(BGT+rangen.GaussR(3,0,SP_r0));
      //TLS_Part.SetRho(rangen.GaussR(3,0,BGT)+0*rangen.GaussR(3,0,SP_r0));
      //TLS_Part.SetRho(rangen.GaussR(3,0,Time)+SP_r0);
      //TLS_Part.SetRho(rangen.GaussR(3,0,SP_r0));

      //those two lines are identical, as expected. Good.
      //next to check: does the angular correlation between k and r changes?
      //if so, check if this implies that the angular correlations are preserved compared to e.g. BGT,
      //or are they now different and depending on rho?
      //TLS_Part.SetRho(rangen.GaussR(3,0,SP_r0));
      //TLS_Part.SetXYZT(rangen.Gauss(0,SP_r0),rangen.Gauss(0,SP_r0),rangen.Gauss(0,SP_r0),Time);

      SP_p_LAB->Fill(TLM_Part.P());
      SP_pT_LAB->Fill(TLM_Part.Pt());
      SP_eta_LAB->Fill(TLM_Part.Eta());
      SP_phi_LAB->Fill(TLM_Part.Phi());
      SP_theta_LAB->Fill(TLM_Part.Theta());
      SP_costheta_LAB->Fill(cos(TLM_Part.Theta()));
      SP_Sr_LAB->Fill(TLS_Part.Rho());

      CatsParticle KittyParticle;
      KittyParticle.Set(TLS_Part.T(),TLS_Part.X(),TLS_Part.Y(),TLS_Part.Z(),TLM_Part.E(),TLM_Part.X(),TLM_Part.Y(),TLM_Part.Z());
      KittyParticle.SetPid(1);
      KittyParticle.SetMass(1000);

      //printf("p = %.0f\n",KittyParticle.GetP());
      KittyEvent_LAB->AddParticle(KittyParticle);
      KittyEvent_CM->AddParticle(KittyParticle);
    }
    KittyEvent_LAB->ComputeParticlePairs(false,false);
    KittyEvent_CM->ComputeParticlePairs(false,true);
    //printf("NP = %u\n", KittyEvent->GetNumPairs());
    for(unsigned uPair=0; uPair<KittyEvent_LAB->GetNumPairs(); uPair++){
      //double rstar = KittyEvent_LAB->GetParticlePair(uPair).GetR();
      //double kstar = KittyEvent_LAB->GetParticlePair(uPair).GetP()*0.5;
      //printf("p0 = %.0f\n",KittyEvent->GetParticlePair(uPair).GetParticle(0).GetP());
      //printf("p1 = %.0f\n",KittyEvent->GetParticlePair(uPair).GetParticle(1).GetP());
      //printf("rstar = %.2f\n",rstar);
      //printf("kstar = %.2f\n",kstar);
      //printf("\n");
      //usleep(500e3);

      double kstar = KittyEvent_CM->GetParticlePair(uPair).GetP()*0.5;
      double Eta0 = KittyEvent_LAB->GetParticlePair(uPair).GetParticle(0).GetPseudoRap();
      double Eta1 = KittyEvent_LAB->GetParticlePair(uPair).GetParticle(1).GetPseudoRap();
      double Mom0 = KittyEvent_LAB->GetParticlePair(uPair).GetParticle(0).GetP();
      double Mom1 = KittyEvent_LAB->GetParticlePair(uPair).GetParticle(1).GetP();
      bool FemtoPair = kstar<FemtoRegion && fabs(Eta0)<EtaRegion && fabs(Eta1)<EtaRegion && Mom0>MinMom && Mom1>MinMom;
      //printf("kstar = %.2f\n",kstar);
      //printf("Eta0 = %.2f\n",Eta0);
      //printf("Eta1 = %.2f\n",Eta1);
      //printf("Mom0 = %.2f\n",Mom0);
      //printf("Mom1 = %.2f\n",Mom1);
      //printf("\n");
      //usleep(200e3);

      Pair_mT_Sr_LAB->Fill(KittyEvent_LAB->GetParticlePair(uPair).GetP()*0.5, KittyEvent_LAB->GetParticlePair(uPair).GetR());
      Pair_mT_Sr_CM->Fill(KittyEvent_CM->GetParticlePair(uPair).GetP()*0.5, KittyEvent_CM->GetParticlePair(uPair).GetR());

      if(FemtoPair){
        double avg_mass;
        double avg_pt;
        double mtee;

        avg_mass = KittyEvent_LAB->GetParticlePair(uPair).GetParticle(0).GetMass();
        avg_mass += KittyEvent_LAB->GetParticlePair(uPair).GetParticle(1).GetMass();
        avg_mass *= 0.5;
        avg_pt = KittyEvent_LAB->GetParticlePair(uPair).GetSum().GetPt()*0.5;
        mtee = sqrt(avg_mass*avg_mass+avg_pt*avg_pt);
        FemPair_mT_Sr_LAB->Fill(mtee, KittyEvent_LAB->GetParticlePair(uPair).GetR());

        avg_mass = KittyEvent_CM->GetParticlePair(uPair).GetParticle(0).GetMass();
        avg_mass += KittyEvent_CM->GetParticlePair(uPair).GetParticle(1).GetMass();
        avg_mass *= 0.5;
        avg_pt = KittyEvent_CM->GetParticlePair(uPair).GetSum().GetPt()*0.5;
        mtee = sqrt(avg_mass*avg_mass+avg_pt*avg_pt);
        FemPair_mT_Sr_CM->Fill(mtee, KittyEvent_CM->GetParticlePair(uPair).GetR());

        FemPair_Sr_LAB->Fill(KittyEvent_LAB->GetParticlePair(uPair).GetR());
        FemPair_Sr_CM->Fill(KittyEvent_CM->GetParticlePair(uPair).GetR());

        FemPair_Angle_kr_CM->Fill(KittyEvent_LAB->GetParticlePair(uPair).GetScatAngle());
      }
    }

    delete KittyEvent_LAB;
    delete KittyEvent_CM;
  }

  SP_Sr_LAB->Scale(1./SP_Sr_LAB->Integral(),"width");
  FemPair_Sr_LAB->Scale(1./FemPair_Sr_LAB->Integral(),"width");
  FemPair_Sr_CM->Scale(1./FemPair_Sr_CM->Integral(),"width");
  TF1* f_SP_Sr = new TF1(TString::Format("f_SP_Sr"),GaussSourceTF1,RadMin,RadMax,1);
  f_SP_Sr->SetParameter(0,1.0);
  f_SP_Sr->SetParLimits(0,0.5,2.0);
  SP_Sr_LAB->Fit(f_SP_Sr,"Q, S, N, R, M");
  printf("Single particle source size (LAB): %.3f fm\n",f_SP_Sr->GetParameter(0));

  TFile fOutput(OutputFolder+"fOutput.root","recreate");
  SP_p_LAB->Write();
  SP_pT_LAB->Write();
  SP_eta_LAB->Write();
  SP_phi_LAB->Write();
  SP_theta_LAB->Write();
  SP_costheta_LAB->Write();
  SP_Sr_LAB->Write();
  f_SP_Sr->Write();
  Pair_mT_Sr_LAB->Write();
  Pair_mT_Sr_CM->Write();
  FemPair_betaT_LAB->Write();
  FemPair_mT_Sr_LAB->Write();
  FemPair_mT_Sr_CM->Write();
  FemPair_Sr_LAB->Write();
  FemPair_Sr_CM->Write();
  //FemPair_Angle_kr_LAB->Write();
  FemPair_Angle_kr_CM->Write();

  delete SP_p_LAB;
  delete SP_pT_LAB;
  delete SP_eta_LAB;
  delete SP_phi_LAB;
  delete SP_theta_LAB;
  delete SP_costheta_LAB;
  delete SP_Sr_LAB;
  delete Pair_mT_Sr_LAB;
  delete Pair_mT_Sr_CM;
  delete FemPair_mT_Sr_LAB;
  delete FemPair_mT_Sr_CM;
  delete FemPair_Sr_LAB;
  delete FemPair_Sr_CM;
  //delete FemPair_Angle_kr_LAB;
  delete FemPair_Angle_kr_CM;
  delete FemPair_betaT_LAB;
}


//returns a tgraph errors with the core values
//we estimate that based on the fit with a power law
//oton did not concider the stat err, now we do
//MODE == X0 : power law
//MODE == X1: oton style, i.e. power law + constant
//MODE == 0Y: fit r_core (pp)
//MODE == 1Y: fit r_eff (pp)
DLM_Histo<float>* Fit_mT(const double nsigma, const int mode=0){

  //we get 16x more than the minimum needed
  //e.g. for 1 sig we have c.a. 1:3 pval, or we need c.a. 3 entries to get one in the side region. We will take 3x16 = 48 entries
  //for 3 sigma is c.a. 1:300 => c.a. 300 x 16 = 4800 entires etc.
  //printf("nsigma = %.2f; pval = %.5f = 1:%.0f\n",nsigma,GetPvalFromNsig(nsigma),1./GetPvalFromNsig(nsigma));
  const unsigned NumIter = TMath::Nint(16./GetPvalFromNsig(nsigma));
  const unsigned NumMtBins_pp = 7;
  double* BinRange_pp = new double[NumMtBins_pp+1];
  double* BinCenter_pp = new double[NumMtBins_pp];
  BinCenter_pp[0] = 1.1077;
  BinCenter_pp[1] = 1.1683;
  BinCenter_pp[2] = 1.2284;
  BinCenter_pp[3] = 1.3156;
  BinCenter_pp[4] = 1.4628;
  BinCenter_pp[5] = 1.6872;
  BinCenter_pp[6] = 2.2116;
  BinRange_pp[0] = BinCenter_pp[0]-BinCenter_pp[1]+BinCenter_pp[0];
  for(int i=1; i<=6; i++) BinRange_pp[i] = (BinCenter_pp[i-1]+BinCenter_pp[i])*0.5;
  BinRange_pp[7] = BinCenter_pp[6]+BinCenter_pp[6]-BinCenter_pp[5];

  DLM_Histo<float> dlm_rcore_pp;
  dlm_rcore_pp.SetUp(1);
  dlm_rcore_pp.SetUp(0,NumMtBins_pp,BinRange_pp,BinCenter_pp);
  dlm_rcore_pp.Initialize();

  if(mode/10==0){
    dlm_rcore_pp.SetBinContent(unsigned(0),1.3064);
    dlm_rcore_pp.SetBinError(unsigned(0),sqrt(pow(0.027805,2.)+pow(0.0085539,2.)));
    dlm_rcore_pp.SetBinContent(1,1.2316);
    dlm_rcore_pp.SetBinError(1,sqrt(pow(0.022773,2.)+pow(0.010501,2.)));
    dlm_rcore_pp.SetBinContent(2,1.2006);
    dlm_rcore_pp.SetBinError(2,sqrt(pow(0.022552,2.)+pow(0.014732,2.)));
    dlm_rcore_pp.SetBinContent(3,1.1402);
    dlm_rcore_pp.SetBinError(3,sqrt(pow(0.025014,2.)+pow(0.011086,2.)));
    dlm_rcore_pp.SetBinContent(4,1.0628);
    dlm_rcore_pp.SetBinError(4,sqrt(pow(0.025221,2.)+pow(0.010027,2.)));
    dlm_rcore_pp.SetBinContent(5,0.96238);
    dlm_rcore_pp.SetBinError(5,sqrt(pow(0.025234,2.)+pow(0.0098228,2.)));
    dlm_rcore_pp.SetBinContent(6,0.86503);
    dlm_rcore_pp.SetBinError(6,sqrt(pow(0.020369,2.)+pow(0.010841,2.)));
  }
  else{
    dlm_rcore_pp.SetBinContent(unsigned(0),1.3717);
    dlm_rcore_pp.SetBinError(unsigned(0),sqrt(pow(0.025241,2.)+pow(0.0079898,2.)));
    dlm_rcore_pp.SetBinContent(1,1.3005);
    dlm_rcore_pp.SetBinError(1,sqrt(pow(0.022749,2.)+pow(0.0098079,2.)));
    dlm_rcore_pp.SetBinContent(2,1.2709);
    dlm_rcore_pp.SetBinError(2,sqrt(pow(0.022641,2.)+pow(0.013852,2.)));
    dlm_rcore_pp.SetBinContent(3,1.2138);
    dlm_rcore_pp.SetBinError(3,sqrt(pow(0.022693,2.)+pow(0.010282,2.)));
    dlm_rcore_pp.SetBinContent(4,1.1423);
    dlm_rcore_pp.SetBinError(4,sqrt(pow(0.027803,2.)+pow(0.0094262,2.)));
    dlm_rcore_pp.SetBinContent(5,1.0482);
    dlm_rcore_pp.SetBinError(5,sqrt(pow(0.022569,2.)+pow(0.0090911,2.)));
    dlm_rcore_pp.SetBinContent(6,0.95678);
    dlm_rcore_pp.SetBinError(6,sqrt(pow(0.017544,2.)+pow(0.0099002,2.)));
  }

  TRandom3 RanGen(11);
  double ran_mean;
  double ran_stdv;

  TF1* fPowerLaw = new TF1("fPowerLaw","[0]*pow(x,[1])+[2]",BinRange_pp[0],BinRange_pp[NumMtBins_pp]);
  fPowerLaw->SetParLimits(0,0,100);
  fPowerLaw->SetParLimits(1,-10,0);
  fPowerLaw->SetParameter(2,0);

  TH2F* h_pp_mt = new TH2F("h_pp_mt","h_pp_mt",NumMtBins_pp,BinRange_pp,2048,0.0,2.5);

  //TFile fOutput(TString::Format("%s/SourceStudies/Estimate_Reff/TEMP.root",GetFemtoOutputFolder()),"recreate");

  for(unsigned uIter=0; uIter<NumIter; uIter++){
    TGraphErrors g_bootstrap;
    g_bootstrap.SetName(TString::Format("g_bootstrap_%u",uIter));
    for(unsigned uBin=0; uBin<dlm_rcore_pp.GetNbins(); uBin++){
      ran_mean = dlm_rcore_pp.GetBinContent(uBin);
      ran_stdv = dlm_rcore_pp.GetBinError(uBin);
      g_bootstrap.SetPoint(uBin,dlm_rcore_pp.GetBinCenter(0,uBin),RanGen.Gaus(ran_mean,ran_stdv));
      g_bootstrap.SetPointError(uBin,0,dlm_rcore_pp.GetBinError(uBin));
    }
    g_bootstrap.Fit(fPowerLaw,"Q, S, N, R, M");
    //g_bootstrap.Write();

    for(unsigned uBin=0; uBin<dlm_rcore_pp.GetNbins(); uBin++){
      h_pp_mt->Fill(g_bootstrap.GetPointX(uBin),fPowerLaw->Eval(g_bootstrap.GetPointX(uBin)));
    }
  }

  DLM_Histo<float>* result = new DLM_Histo<float>();
  result->SetUp(1);
  result->SetUp(0,NumMtBins_pp,BinRange_pp,BinCenter_pp);
  result->Initialize();

  for(unsigned uBin=0; uBin<dlm_rcore_pp.GetNbins(); uBin++){
    TH1F* hProj = (TH1F*)h_pp_mt->ProjectionY(TString::Format("hProj"),uBin+1,uBin+1);
    double lowerlimit;
    double upperlimit;
    double median = GetCentralInterval(*hProj, 1.-GetPvalFromNsig(nsigma), lowerlimit, upperlimit, true);
    result->SetBinContent(uBin,(lowerlimit+upperlimit)*0.5);
    result->SetBinError(uBin,(upperlimit-lowerlimit)*0.5);
  }

  delete [] BinRange_pp;
  delete [] BinCenter_pp;
  delete fPowerLaw;
  delete h_pp_mt;

  return result;
}

//type = oton => power law + const
//type = dimi => power law
//model = {RSM_PLB, RSM_GHETTO, RSM_GHETTO_FLAT}
void Estimate_Reff(const TString system, const TString type, const TString model, const double nsigma=3){

  double reff_min;
  double reff_avg;
  double reff_max;
  double reff_err;

  double rcore_min;
  double rcore_avg;
  double rcore_max;
  double rcore_err;

  const unsigned NumRadBins = 256;
  const double rMin = 0;
  const double rMax = 32;

  //there seems to be a consistent systematic shift of c.a. 0.005 between the pp true fit and the
  //GHETTO result. To be safe, we add 3sigma of that as error of the method (in quadriture)
  const double GHETTO_ERROR = 0.005*3;

  double avg_mT;

  if(type!="oton"&&type!="dimi"){
    printf("type == %s ???\n", type.Data());
    return;
  }
  if(model!="RSM_PLB"&&model!="RSM_GHETTO"&&model!="RSM_FLAT"){
    printf("model == %s ???\n", model.Data());
    return;
  }

  DLM_Histo<float>* dlm_fit_rcore_pp= NULL;
  DLM_Histo<float>* dlm_fit_reff_pp= NULL;
  //used to estimate the difference between GHETTO and DEFAULT and adds this as an error
  //DLM_Histo<float>* dlm_fit_rcoreDEF_pp= NULL;
  //the original published result
  if(model=="RSM_PLB"){
    if(type=="oton") dlm_fit_rcore_pp = Fit_mT(nsigma,1);
    else if(type=="dimi") dlm_fit_rcore_pp = Fit_mT(nsigma,0);
    else return;
  }
  //we read off the pp reff, convert it into a rcore by a simple iterative procedure,
  //and use that r_core instead of the published one
  else if(model=="RSM_GHETTO"){
    if(type=="oton") {
      //dlm_fit_rcoreDEF_pp = Fit_mT(nsigma,1);
      dlm_fit_reff_pp = Fit_mT(nsigma,11);
    }
    else if(type=="dimi") {
      //dlm_fit_rcoreDEF_pp = Fit_mT(nsigma,0);
      dlm_fit_reff_pp = Fit_mT(nsigma,10);
    }
    else return;

    DLM_CleverMcLevyResoTM MagicSource_Epos_pp;
    SetUp_RSM_pp(MagicSource_Epos_pp,GetCernBoxDimi(),0);
    dlm_fit_rcore_pp = new DLM_Histo<float> (*dlm_fit_reff_pp);
    for(unsigned umT=0; umT<dlm_fit_reff_pp->GetNbins(); umT++){
      double rmin,rmax;
      rmin = dlm_fit_reff_pp->GetBinContent(umT) - dlm_fit_reff_pp->GetBinError(umT);
      rmax = dlm_fit_reff_pp->GetBinContent(umT) + dlm_fit_reff_pp->GetBinError(umT);
      rmin = GetRcore(MagicSource_Epos_pp,rmin);
      rmax = GetRcore(MagicSource_Epos_pp,rmax);
      dlm_fit_rcore_pp->SetBinContent(umT,(rmax+rmin)*0.5);
      dlm_fit_rcore_pp->SetBinError(umT,(rmax-rmin)*0.5);
    }
  }
  else if(model=="RSM_FLAT"){
    if(type=="oton") {
      //dlm_fit_rcoreDEF_pp = Fit_mT(nsigma,1);
      dlm_fit_reff_pp = Fit_mT(nsigma,11);
    }
    else if(type=="dimi") {
      //dlm_fit_rcoreDEF_pp = Fit_mT(nsigma,0);
      dlm_fit_reff_pp = Fit_mT(nsigma,10);
    }
    else return;

    DLM_CleverMcLevyResoTM MagicSource_Flat_pp;
    SetUp_RSMflat_pp(MagicSource_Flat_pp);
    dlm_fit_rcore_pp = new DLM_Histo<float> (*dlm_fit_reff_pp);
    for(unsigned umT=0; umT<dlm_fit_reff_pp->GetNbins(); umT++){
      double rmin,rmax;
      rmin = dlm_fit_reff_pp->GetBinContent(umT) - dlm_fit_reff_pp->GetBinError(umT);
      rmax = dlm_fit_reff_pp->GetBinContent(umT) + dlm_fit_reff_pp->GetBinError(umT);
      rmin = GetRcore(MagicSource_Flat_pp,rmin);
      rmax = GetRcore(MagicSource_Flat_pp,rmax);
      dlm_fit_rcore_pp->SetBinContent(umT,(rmax+rmin)*0.5);
      dlm_fit_rcore_pp->SetBinError(umT,(rmax-rmin)*0.5);
    }
  }
  else return;

  DLM_CleverMcLevyResoTM MagicSource_system;

  if(system=="pOmega"){
    avg_mT = 2.20;
    if(model=="RSM_FLAT"){
      SetUp_RSMflat_pOmega(MagicSource_system);
    }
    else SetUp_RSM_pOmega(MagicSource_system,GetCernBoxDimi(),0);
  }
  else if(system=="pd"||system=="pd_min"||system=="pd_max"){
    avg_mT = 1.64;//1.65
    if(model=="RSM_FLAT"){
      SetUp_RSMflat_pOmega(MagicSource_system);
    }
    else SetUp_RSM_pOmega(MagicSource_system,GetCernBoxDimi(),0);
    if(system=="pd_min") MagicSource_system.SetUpReso(0,0.6422-0.06422);
    if(system=="pd_max") MagicSource_system.SetUpReso(0,0.6422+0.06422);
  }
  else if(system=="pLambda"){
    avg_mT = 1.55;
    if(model=="RSM_FLAT"){
      SetUp_RSMflat_pL(MagicSource_system);
    }
    else SetUp_RSM_pL(MagicSource_system,GetCernBoxDimi(),0);

  }
  else if(system=="pXi"){
    avg_mT = 1.90;
    if(model=="RSM_FLAT"){
      SetUp_RSMflat_pXi(MagicSource_system);
    }
    else SetUp_RSM_pXi(MagicSource_system,GetCernBoxDimi(),0);
  }

  rcore_avg = dlm_fit_rcore_pp->Eval(&avg_mT);
  rcore_err = dlm_fit_rcore_pp->Eval(&avg_mT,true);
  if(model=="RSM_GHETTO"){
    rcore_err = sqrt(GHETTO_ERROR*GHETTO_ERROR+rcore_err*rcore_err);
  }
  //if(dlm_fit_rcoreDEF_pp){
  //  double rcoredef = dlm_fit_rcoreDEF_pp->Eval(&avg_mT);
  //  double diff = rcore_avg-rcoredef;
  //  rcore_err = sqrt(rcore_err*rcore_err+diff*diff);
  //}
  rcore_min = rcore_avg - rcore_err;
  rcore_max = rcore_avg + rcore_err;


  TH1F* hSource = new TH1F("hSource","hSource",NumRadBins,rMin,rMax);

  for(unsigned uRad=0; uRad<NumRadBins; uRad++){
    double parameters[2];
    double rstar = hSource->GetBinCenter(uRad+1);
    parameters[0] = rcore_min;
    parameters[1] = 2.0;
    double val = MagicSource_system.RootEval(&rstar,parameters);
    hSource->SetBinContent(uRad+1,val);
  }
  reff_min = Get_reff(hSource,1,0.84);

  for(unsigned uRad=0; uRad<NumRadBins; uRad++){
    double parameters[2];
    double rstar = hSource->GetBinCenter(uRad+1);
    parameters[0] = rcore_max;
    parameters[1] = 2.0;
    double val = MagicSource_system.RootEval(&rstar,parameters);
    hSource->SetBinContent(uRad+1,val);
  }
  reff_max = Get_reff(hSource,1,0.84);

  TF1* freff = new TF1();
  delete freff;
  for(unsigned uRad=0; uRad<NumRadBins; uRad++){
    double parameters[2];
    double rstar = hSource->GetBinCenter(uRad+1);
    parameters[0] = rcore_avg;
    parameters[1] = 2.0;
    double val = MagicSource_system.RootEval(&rstar,parameters);
    hSource->SetBinContent(uRad+1,val);
  }


  reff_avg = (reff_max+reff_min)*0.5;
  reff_err = (reff_max-reff_min)*0.5;

  printf("%s source (%s):\n",system.Data(),model.Data());
  printf(" <mT> =  %.2f GeV\n",avg_mT);
  printf(" rcore = %.3f +/- %.3f fm\n", rcore_avg, rcore_err);
  printf("  reff = %.3f +/- %.3f fm\n", reff_avg, reff_err);
  printf("    QA = %.3f\n",Get_reff_TF1(hSource,freff,1,0.84));
  TF1 *freff_full = (TF1*)freff->Clone("freff_full");
  freff_full->SetRange(rMin,rMax);

  TFile fOutput(TString::Format("%s/SourceStudies/Estimate_Reff/Reff_%s_%s_%s_%.1f.root",GetFemtoOutputFolder(),system.Data(),type.Data(),model.Data(),nsigma),"recreate");
  TGraphErrors g_fit_rcore_pp;
  g_fit_rcore_pp.SetName("g_fit_rcore_pp");
  for(unsigned uBin=0; uBin<dlm_fit_rcore_pp->GetNbins(); uBin++){
    g_fit_rcore_pp.SetPoint(uBin,dlm_fit_rcore_pp->GetBinCenter(0,uBin),dlm_fit_rcore_pp->GetBinContent(uBin));
    g_fit_rcore_pp.SetPointError(uBin,0,dlm_fit_rcore_pp->GetBinError(uBin));
  }
  g_fit_rcore_pp.Write();
  hSource->Write();
  if(freff) freff->Write();
  freff_full->Write();

  if(dlm_fit_rcore_pp) delete dlm_fit_rcore_pp;
  if(dlm_fit_reff_pp) delete dlm_fit_reff_pp;
  delete hSource;
  if(freff) delete freff;
}

//a sphere, on which we randomly populate particles
//we build up all their distances (i.e. building pairs) but than we concider
//to have a reconstruction efficiency, where we only keep a fraction of random particles.
//we want to see if the final two-particle source changes.

//RESULT: the efficiency and acceptence do NOT matter for the average distance we see between the pairs.
//Also, what surprised me, the density does not play a direct role as well, i.e. for higher density we still get
//the same average distace between pairs (more pairs = more combinatorics, so averages do not change)
//The ONLY think on which the pair distance depends is the SIZE of the system itself. SIZE (Volume if you will), not density!!!
void SourceDensity(){
  const unsigned NumEvents = 1;
  const double SystemSize = 2.0;
  const unsigned Multiplicity = 10000;
  const double Efficiecy = 0.5;

  const double EtaAcc = 0.5;
  //if the diff in Theta/Phi is no larger than that number, we take the pair
  //simulates the region of homogenity
  const double dAngle = 0.2;

  double Theta,Phi,CosTheta,CosPhi,SinPhi,SinTheta,xcrd,ycrd,zcrd;
  TRandom3 rangen(11);

  TH1F* h_X = new TH1F("h_X","h_X",512,0,32);
  TH1F* h_Y = new TH1F("h_Y","h_Y",512,0,32);
  TH1F* h_Z = new TH1F("h_Z","h_Z",512,0,32);
  TH1F* h_CosTh = new TH1F("h_CosTh","h_CosTh",512,-1,1);
  TH1F* h_Phi = new TH1F("h_Phi","h_Phi",512,0,2.*TMath::Pi());

  TH1F* h_Rad = new TH1F("h_Rad","h_Rad",512,0,32);
  TH1F* h_RadDet = new TH1F("h_RadDet","h_RadDet",512,0,32);


  TH1F* h_Rad_Acc = new TH1F("h_Rad_Acc","h_Rad_Acc",512,0,32);
  TH1F* h_RadDet_Acc = new TH1F("h_RadDet_Acc","h_RadDet_Acc",512,0,32);

  TH1F* h_Rad_AccHom = new TH1F("h_Rad_AccHom","h_Rad_AccHom",512,0,32);
  TH1F* h_RadDet_AccHom = new TH1F("h_RadDet_AccHom","h_RadDet_AccHom",512,0,32);

  CatsLorentzVector* Particle = new CatsLorentzVector[Multiplicity];
  CatsLorentzVector* ParticleDet = new CatsLorentzVector[Multiplicity];

  for(unsigned uEvent=0; uEvent<NumEvents; uEvent++){
    unsigned DetectedPart = 0;
    for(unsigned uPart=0; uPart<Multiplicity; uPart++){
      CosTheta = rangen.Uniform(-1,1);
      Theta = acos(CosTheta);
      SinTheta = sin(Theta);
      Phi = rangen.Uniform(0,2.*TMath::Pi());
      CosPhi = cos(Phi);
      SinPhi = sin(Phi);
      xcrd = SystemSize*CosPhi*SinTheta;
      ycrd = SystemSize*SinPhi*SinTheta;
      zcrd = SystemSize*CosTheta;
      Particle[uPart].SetTXYZ(0,xcrd,ycrd,zcrd);

      //detected
      if(rangen.Uniform()<Efficiecy){
        ParticleDet[DetectedPart].SetTXYZ(0,xcrd,ycrd,zcrd);
        DetectedPart++;
      }

      h_X->Fill(xcrd);
      h_Y->Fill(ycrd);
      h_Z->Fill(zcrd);

      h_CosTh->Fill(CosTheta);
      h_Phi->Fill(Phi);
    }

    CatsLorentzVector Diff;
    //particle pairs
    for(unsigned uPart0=0; uPart0<Multiplicity; uPart0++){
      for(unsigned uPart1=uPart0+1; uPart1<Multiplicity; uPart1++){
        Diff = Particle[uPart0] - Particle[uPart1];
        h_Rad->Fill(Diff.GetR());
        if(fabs(Particle[uPart0].GetPseudoRapR())<EtaAcc && fabs(Particle[uPart1].GetPseudoRapR())<EtaAcc){
          h_Rad_Acc->Fill(Diff.GetR());
          if(fabs(Particle[uPart0].GetPhi()-Particle[uPart1].GetPhi())<dAngle && fabs(Particle[uPart0].GetTheta()-Particle[uPart1].GetTheta())<dAngle){
            h_Rad_AccHom->Fill(Diff.GetR());
          }
        }
      }
    }
    for(unsigned uPart0=0; uPart0<DetectedPart; uPart0++){
      for(unsigned uPart1=uPart0+1; uPart1<DetectedPart; uPart1++){
        Diff = ParticleDet[uPart0] - ParticleDet[uPart1];
        h_RadDet->Fill(Diff.GetR());
        if(fabs(ParticleDet[uPart0].GetPseudoRapR())<EtaAcc && fabs(ParticleDet[uPart1].GetPseudoRapR())<EtaAcc){
          h_RadDet_Acc->Fill(Diff.GetR());
          if(fabs(ParticleDet[uPart0].GetPhi()-ParticleDet[uPart1].GetPhi())<dAngle && fabs(ParticleDet[uPart0].GetTheta()-ParticleDet[uPart1].GetTheta())<dAngle){
            h_RadDet_AccHom->Fill(Diff.GetR());
          }
        }
      }
    }

  }


  TFile fOutput(TString::Format("%s/SourceStudies/SourceDensity.root",GetFemtoOutputFolder()),"recreate");
  h_X->Write();
  h_Y->Write();
  h_Z->Write();
  h_CosTh->Write();
  h_Phi->Write();
  h_Rad->Write();
  h_RadDet->Write();
  h_Rad_Acc->Write();
  h_RadDet_Acc->Write();
  h_Rad_AccHom->Write();
  h_RadDet_AccHom->Write();

  delete [] Particle;
  delete [] ParticleDet;
  delete h_X;
  delete h_Y;
  delete h_Z;
  delete h_CosTh;
  delete h_Phi;
  delete h_Rad;
  delete h_RadDet;
  delete h_Rad_Acc;
  delete h_RadDet_Acc;
  delete h_Rad_AccHom;
  delete h_RadDet_AccHom;
}


void source_pp_syst(){

  std::vector<float> bin_range = {1.02,1.14, 1.20, 1.26, 1.38, 1.56, 1.86, 4.50};
  std::vector<float> bin_center = {1.11,1.17, 1.23, 1.32, 1.46, 1.69, 2.21};

  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());

  bool InfoFor_ppForce = false;
  const unsigned NumDataVar = 27;
  const unsigned RebinFactor = InfoFor_ppForce==true?1:5;

  TFile fOutput(TString::Format("%s/SourceStudies/source_pp_syst/fOutput.root",GetFemtoOutputFolder()),"recreate");

  for(unsigned uMt=0; uMt<bin_center.size(); uMt++){

    TH2F* hSystOnly = NULL;

    double* CkMin = NULL;
    double* CkMax = NULL;

    for(unsigned uDV=0; uDV<NumDataVar; uDV++){
      TH1F* hData = AnalysisObject.GetAliceExpCorrFun("pp13TeV_HM_BernieSource","pp",TString::Format("%u",uDV),0,0,uMt);
      hData->Rebin(RebinFactor);
      if(!hSystOnly){
        unsigned NumBins = hData->GetNbinsX();
        double LowEdge = hData->GetXaxis()->GetBinLowEdge(1);
        double UpEdge = hData->GetXaxis()->GetBinUpEdge(NumBins);
        hSystOnly = new TH2F(TString::Format("hSystOnly_%u",uMt),TString::Format("hSystOnly_%u",uMt),NumBins,LowEdge,UpEdge,65536,0,6);
        CkMin = new double[NumBins];
        CkMax = new double[NumBins];
        for(unsigned uMom=0; uMom<hData->GetNbinsX(); uMom++){
          CkMin[uMom] = 1000;
          CkMax[uMom] = 0;
        }
      }

      for(unsigned uMom=0; uMom<hData->GetNbinsX(); uMom++){
        double MOM = hData->GetBinCenter(uMom+1);
        double CkVal = hData->GetBinContent(uMom+1);
        if(CkVal<CkMin[uMom]){CkMin[uMom] = CkVal;}
        if(CkVal>CkMax[uMom]){CkMax[uMom] = CkVal;}
        hSystOnly->Fill(MOM,CkVal);
      }
      delete hData;
    }//uDV

    TH1F* hDataReb = AnalysisObject.GetAliceExpCorrFun("pp13TeV_HM_BernieSource","pp",TString::Format("%u",0),0,0,uMt);
    hDataReb->SetName(TString::Format("hReb_mT%u",uMt+1));
    hDataReb->Rebin(RebinFactor);
    //TH1F* hSyst = (TH1F*)hDataReb->Clone(TString::Format("hSyst_mT%u",uMt+1));
    TH1F* hRelSystErr = (TH1F*)hDataReb->Clone(TString::Format("hRelSystErr_mT%u",uMt+1));

    TH1F* hDataStat = AnalysisObject.GetAliceExpCorrFun("pp13TeV_HM_BernieSource","pp",TString::Format("%u",0),0,0,uMt);
    hDataStat->SetName(TString::Format("hStat_mT%u",uMt+1));

    TH1F* hDataSyst = AnalysisObject.GetAliceExpCorrFun("pp13TeV_HM_BernieSource","pp",TString::Format("%u",0),0,0,uMt);
    hDataSyst->SetName(TString::Format("hSyst_mT%u",uMt+1));


    TList* list1_tmp;
    TList* list2_tmp;

    TString SystFileName = TString::Format("%s/CatsFiles/ExpData/Bernie_Source/ppData/mTBin_%i/CFOutput_mT_ppVar0_HM_%i.root",GetCernBoxDimi(), uMt+1, uMt);
    TFile* FileROOT = new TFile(SystFileName, "read");
    //FileROOT->ls();

    list1_tmp = (TList*)FileROOT->Get("PairDist");
    list2_tmp = (TList*)list1_tmp->FindObject("PairRebinned");
    TH1F* hMixedEvent_PP = (TH1F*)list2_tmp->FindObject(TString::Format("MEPart_mT_%u_FixShifted_Rebinned_2",uMt));

    list1_tmp = (TList*)FileROOT->Get("AntiPairDist");
    list2_tmp = (TList*)list1_tmp->FindObject("PairRebinned");

    TH1F* hMixedEvent_APAP = (TH1F*)list2_tmp->FindObject(TString::Format("MEAntiPart_mT_%u_FixShifted_Rebinned_2",uMt));
    gROOT->cd();
    TH1F* hMixedEvent = (TH1F*)hMixedEvent_PP->Clone(TString::Format("hMixedEvent_mT%u",uMt+1));
    hMixedEvent->Add(hMixedEvent_APAP);

    list1_tmp = (TList*)FileROOT->Get("PairDist");
    list2_tmp = (TList*)list1_tmp->FindObject("PairRebinned");
    TH1F* hSameEvent_PP = (TH1F*)list2_tmp->FindObject(TString::Format("SEPart_mT_%u_FixShifted_Rebinned_2",uMt));
    list1_tmp = (TList*)FileROOT->Get("AntiPairDist");
    list2_tmp = (TList*)list1_tmp->FindObject("PairRebinned");
    TH1F* hSameEvent_APAP = (TH1F*)list2_tmp->FindObject(TString::Format("SEAntiPart_mT_%u_FixShifted_Rebinned_2",uMt));
    gROOT->cd();
    TH1F* hSameEvent = (TH1F*)hSameEvent_PP->Clone(TString::Format("hSameEvent_mT%u",uMt+1));
    hSameEvent->Add(hSameEvent_APAP);
    //printf("HI\n");
    delete FileROOT;

    printf("mT bin in [%.2f, %.2f)\n",bin_range.at(uMt),bin_range.at(uMt+1));
    printf(" # SE pairs < 100 MeV: %.0fk\n", hSameEvent->Integral(1,hSameEvent->FindBin(0.1))*0.001);
    printf(" # SE pairs < 200 MeV: %.0fk\n", hSameEvent->Integral(1,hSameEvent->FindBin(0.2))*0.001);

    TGraphErrors gDataSyst;
    TGraphErrors gDataStat;

    for(unsigned uMom=0; uMom<hDataReb->GetNbinsX(); uMom++){
      double ERROR = ( (CkMax[uMom]-CkMin[uMom])/sqrt(12.) );
      double CKVAL = hDataReb->GetBinContent(uMom+1);

      hRelSystErr->SetBinContent(uMom+1,ERROR/CKVAL);
      hRelSystErr->SetBinError(uMom+1,ERROR/CKVAL*0.1);
    }//uMom

    TF1* fRelSystErr = new TF1(TString::Format("fRelSystErr_mT%u",uMt+1),"[0]*exp(-[1]*x)+[2]",0,240);
    fRelSystErr->SetParameter(0,0.035);
    fRelSystErr->SetParLimits(0,0.02,0.6);

    fRelSystErr->SetParameter(1,0.05);
    fRelSystErr->SetParLimits(1,0,0.2);

    fRelSystErr->SetParameter(2,hRelSystErr->GetBinContent(hRelSystErr->FindBin(200)));
    fRelSystErr->SetParLimits(2,fRelSystErr->GetParameter(2)*0.5,fRelSystErr->GetParameter(2)*2);


    for(unsigned uMom=0; uMom<hDataStat->GetNbinsX(); uMom++){
      double MOM = hDataStat->GetBinCenter(uMom+1);
      double EVAL_MOM = MOM>240?240:MOM;
      double CKVAL = hDataStat->GetBinContent(uMom+1);
      //hDataStat->SetErrorX(uMom+1);
      double REL_ERR = fRelSystErr->Eval(EVAL_MOM);
      double ERROR = CKVAL*REL_ERR;




      hDataSyst->SetBinError(uMom+1,ERROR);

      gDataStat.SetPoint(uMom,MOM,CKVAL);
      gDataStat.SetPointError(uMom,0,hDataStat->GetBinError(uMom+1));

      gDataSyst.SetPoint(uMom,MOM,CKVAL);
      gDataSyst.SetPointError(uMom,2,ERROR);
      if(uMt==3 && MOM<10){
        //printf("%f +/- %f\n",CKVAL,ERROR);
        //usleep(1500e3);
      }
    }

    SetStyle();


    hDataStat->SetTitle("; #it{k*} (MeV/#it{c}); #it{C}(#it{k*})");
    hDataStat->GetXaxis()->SetRangeUser(0, 400);
    hDataStat->GetXaxis()->SetNdivisions(505);
    if(uMt<=3) hDataStat->GetYaxis()->SetRangeUser(0.85, 3.25);
    else hDataStat->GetYaxis()->SetRangeUser(0.85, 4.85);
    hDataStat->SetFillColor(kGray+1);
    hDataStat->SetMarkerSize(1.5);
    hDataStat->SetLineWidth(2);
    hDataStat->SetMarkerStyle(kOpenCircle);
    hDataStat->SetMarkerColor(kBlack);
    hDataStat->SetLineColor(kBlack);

    TH1F* hAxis = new TH1F("hAxis", "hAxis", 405, 0, 405);
    hAxis->SetTitle("; #it{k*} (MeV/#it{c}); #it{C}(#it{k*})");
    hAxis->GetXaxis()->SetRangeUser(0, 405);
    hAxis->GetXaxis()->SetNdivisions(505);
    if(uMt<=3) hAxis->GetYaxis()->SetRangeUser(0.85, 3.25);
    else hAxis->GetYaxis()->SetRangeUser(0.85, 4.85);
    hAxis->SetFillColor(kGray+1);
    hAxis->SetMarkerSize(1.6);
    hAxis->SetLineWidth(2);
    hAxis->SetMarkerStyle(kOpenCircle);
    hAxis->SetMarkerColor(kBlack);
    hAxis->SetLineColor(kGray+1);
    hAxis->SetLineWidth(0);



    gDataSyst.SetName("gSystError");
    gDataSyst.SetLineColor(kWhite);
    gDataSyst.SetFillColorAlpha(kBlack, 0.4);

    TLegend *legend = new TLegend(0.47,0.785-0.055-0.060*1,0.94,0.785-0.055);//lbrt
    legend->SetBorderSize(0);
    legend->SetTextFont(42);
    legend->SetTextSize(gStyle->GetTextSize()*0.8);
    legend->AddEntry(hAxis,TString::Format("p#minusp #oplus #bar{p}#minus#bar{p}"),"fpe");



//CORR FUN
    TCanvas* c_PP = new TCanvas("CFpp", "CFpp", 0, 0, 650, 650);
    c_PP->SetMargin(0.13,0.025,0.12,0.025);//lrbt

    hAxis->Draw("axis");
    hDataStat->Draw("E0 X0 same");
    gDataSyst.Draw("2 same");
    legend->Draw("same");

    TLatex BeamText;
    TLatex text;
    BeamText.SetTextSize(gStyle->GetTextSize() * 0.75);
    BeamText.SetNDC(kTRUE);
    //BeamText.DrawLatex(0.32, 0.91, Form("#bf{ALICE}"));
    BeamText.DrawLatex(0.49, 0.9,
                       Form("ALICE %s #sqrt{#it{s}} = %i TeV", "pp", (int) 13));
    BeamText.DrawLatex(0.49, 0.85, "High-mult. (0#minus0.17% INEL > 0)");
    text.SetNDC();
    text.SetTextColor(1);
    text.SetTextSize(gStyle->GetTextSize() * 0.75);
    text.DrawLatex(0.49, 0.80, TString::Format("#it{m}_{T}#in [%.2f, %.2f) GeV",bin_range.at(uMt),bin_range.at(uMt+1)));
    text.DrawLatex(0.49, 0.75, TString::Format("<#it{m}_{T}> = %.2f GeV",bin_center.at(uMt)));

    //text.DrawLatex(0.32, 0.73, "Gaussian source");

    if(!InfoFor_ppForce)
      c_PP->SaveAs(TString::Format("%s/SourceStudies/source_pp_syst/c_PP_%u.pdf",GetFemtoOutputFolder(),uMt));


    TString HepFileName = TString::Format("%s/SourceStudies/source_pp_syst/HEP_PP_%u.yaml",GetFemtoOutputFolder(),uMt);
    ofstream hepfile (HepFileName.Data(),ios::out);
    hepfile << "dependent_variables:" << endl;
    hepfile << "- header: {name: proton-proton C(k*) @ <mT>="<< TMath::Nint(bin_center.at(uMt)*1000.) <<" MeV}" << endl;
    hepfile << "  qualifiers:" << endl;
    hepfile << "  - {name: SQRT(S), units: GeV, value: '13000.0'}" << endl;
    hepfile << "  values:" << endl;
    double hep_syst;
    double hep_stat;
    double hep_x;
    double hep_y;
    for(unsigned uBin=1; uBin<hDataStat->GetNbinsX(); uBin++){
      hep_x = hDataStat->GetBinCenter(uBin);
      hep_y = hDataStat->GetBinContent(uBin);
      if(hep_x>400) break;
      hep_syst = gDataSyst.GetErrorY(uBin-1);
      hep_stat = hDataStat->GetBinError(uBin);
      hepfile << "  - errors:" << endl;
      hepfile << "    - {label: stat, symerror: "<<hep_stat<<"}" << endl;
      hepfile << "    - {label: sys, symerror: "<<hep_syst<<"}" << endl;
      hepfile << "    value: "<<hep_y << endl;
    }
    hepfile << "independent_variables:" << endl;
    hepfile << "- header: {name: k* (MeV/c)}" << endl;
    hepfile << "  values:" << endl;
    for(unsigned uBin=1; uBin<hDataStat->GetNbinsX(); uBin++){
        hep_x = hDataStat->GetBinCenter(uBin);
        if(hep_x>400) break;
        hepfile << "  - {value: "<<hep_x<<"}" << endl;
    }
    hepfile.close();


//MIXED EVENTS

hMixedEvent->SetTitle("; #it{k*} (MeV/#it{c}); #it{M}(#it{k*})");
hMixedEvent->GetXaxis()->SetRangeUser(0, 400);
hMixedEvent->GetXaxis()->SetNdivisions(505);
hMixedEvent->GetYaxis()->SetRangeUser(0.9*1000, 20000*1000);
//hMixedEvent->SetFillColor(kGray+1);
hMixedEvent->SetMarkerSize(1.5);
hMixedEvent->SetLineWidth(2);
hMixedEvent->SetMarkerStyle(kOpenCircle);
hMixedEvent->SetMarkerColor(kBlack);
hMixedEvent->SetLineColor(kBlack);
hMixedEvent->GetXaxis()->SetLimits(hMixedEvent->GetXaxis()->GetXmin()*1000.,hMixedEvent->GetXaxis()->GetXmax()*1000.);

//hMixedEvent->Scale(1e-3);

TLegend *legendME = new TLegend(0.47,0.385-0.055-0.060*1,0.94,0.385-0.055);//lbrt
legendME->SetBorderSize(0);
legendME->SetTextFont(42);
legendME->SetTextSize(gStyle->GetTextSize()*0.8);
legendME->AddEntry(hMixedEvent,TString::Format("p#minusp #oplus #bar{p}#minus#bar{p}"),"p");

TCanvas* cME_PP = new TCanvas("cME_PP", "cME_PP", 0, 0, 650, 650);
cME_PP->SetMargin(0.13,0.025,0.12,0.025);//lrbt
cME_PP->SetLogy();

hAxis->SetTitle("; #it{k*} (MeV/#it{c}); #it{M}(#it{k*})");
hAxis->GetXaxis()->SetRangeUser(0, 405);
hAxis->GetXaxis()->SetNdivisions(505);
hAxis->GetYaxis()->SetRangeUser(0.9*1000, 20000*1000);

hAxis->Draw("axis");
hMixedEvent->Draw("E0 X0 same");
//gDataSyst.Draw("2 same");
legendME->Draw("same");

BeamText.DrawLatex(0.49, 0.5,
                   Form("ALICE %s #sqrt{#it{s}} = %i TeV", "pp", (int) 13));
BeamText.DrawLatex(0.49, 0.45, "High-mult. (0#minus0.17% INEL > 0)");
text.DrawLatex(0.49, 0.40, TString::Format("#it{m}_{T}#in [%.2f, %.2f) GeV",bin_range.at(uMt),bin_range.at(uMt+1)));
text.DrawLatex(0.49, 0.35, TString::Format("<#it{m}_{T}> = %.2f GeV",bin_center.at(uMt)));

if(!InfoFor_ppForce)
  cME_PP->SaveAs(TString::Format("%s/SourceStudies/source_pp_syst/cME_PP_%u.pdf",GetFemtoOutputFolder(),uMt));


TString HepFileName_ME = TString::Format("%s/SourceStudies/source_pp_syst/HEP_ME_PP_%u.yaml",GetFemtoOutputFolder(),uMt);
ofstream hepfile_ME (HepFileName_ME.Data(),ios::out);
hepfile_ME << "dependent_variables:" << endl;
hepfile_ME << "- header: {name: proton-proton M(k*) @ <mT>="<< TMath::Nint(bin_center.at(uMt)*1000.) <<" MeV}" << endl;
hepfile_ME << "  qualifiers:" << endl;
hepfile_ME << "  - {name: SQRT(S), units: GeV, value: '13000.0'}" << endl;
hepfile_ME << "  values:" << endl;
for(unsigned uBin=1; uBin<hMixedEvent->GetNbinsX(); uBin++){
  hep_x = hMixedEvent->GetBinCenter(uBin);
  hep_y = hMixedEvent->GetBinContent(uBin);
  if(hep_x>400) break;
  hep_syst = 0;
  hep_stat = sqrt(hep_y);
  hepfile_ME << "  - errors:" << endl;
  hepfile_ME << "    - {label: stat, symerror: "<<hep_stat<<"}" << endl;
  hepfile_ME << "    value: "<<hep_y << endl;
}
hepfile_ME << "independent_variables:" << endl;
hepfile_ME << "- header: {name: k* (MeV/c)}" << endl;
hepfile_ME << "  values:" << endl;
for(unsigned uBin=1; uBin<hMixedEvent->GetNbinsX(); uBin++){
    hep_x = hMixedEvent->GetBinCenter(uBin);
    if(hep_x>400) break;
    hepfile_ME << "  - {value: "<<hep_x<<"}" << endl;
}
hepfile_ME.close();




    fOutput.cd();
    hSystOnly->Write();
    hDataReb->Write();
    //hSyst->Write();
    hRelSystErr->Write();
    fRelSystErr->Write();
    hDataStat->Write();
    hDataSyst->Write();
    hMixedEvent->Write();
    hSameEvent->Write();


    delete [] CkMin;
    delete [] CkMax;
    delete hSystOnly;
    delete hDataReb;
    delete hDataStat;
    delete hDataSyst;
    //delete hSyst;
    delete hRelSystErr;
    delete fRelSystErr;

    delete hAxis;
    delete legend;
    delete legendME;
    delete c_PP;
    delete cME_PP;

    delete hMixedEvent;
    delete hSameEvent;

  }//uMt




  //MOMENTUM SMEAR
  TH2F* hReso = AnalysisObject.GetResolutionMatrix("pp13TeV_HM_DimiJun20","pp");
  std::vector<float> kstar_proj = {10, 20, 50, 100, 200, 400};
  TH1F** hSmear = new TH1F* [kstar_proj.size()];

  fOutput.cd();
  for(unsigned uProj=0; uProj<kstar_proj.size(); uProj++){
    double kstar = kstar_proj.at(uProj);
    unsigned WhichBin = hReso->GetYaxis()->FindBin(kstar);
    hSmear[uProj] = (TH1F*)hReso->ProjectionX(TString::Format("hTransform_k%.0f",kstar),WhichBin,WhichBin);
    hSmear[uProj]->Rebin(4);
    double kstar_min = 0;
    double kstar_max = 0;
    double val_min = 1e6;
    double val_max = 0;
    for(unsigned uMom=0; uMom<hSmear[uProj]->GetNbinsX(); uMom++){
      if(hSmear[uProj]->GetBinContent(uMom+1)<10){
        hSmear[uProj]->SetBinContent(uMom+1,0);
        hSmear[uProj]->SetBinError(uMom+1,0);
        continue;
      }
      hSmear[uProj]->SetBinError(uMom+1,sqrt(hSmear[uProj]->GetBinContent(uMom+1)));
      if(kstar_min==0) kstar_min = hSmear[uProj]->GetXaxis()->GetBinLowEdge(uMom+1);
      kstar_max = hSmear[uProj]->GetXaxis()->GetBinUpEdge(uMom+1);
    }

    double INT_BEFORE = hSmear[uProj]->Integral();
    hSmear[uProj]->Scale(1./hSmear[uProj]->Integral(),"width");
    double INT_AFTER = hSmear[uProj]->Integral();
    double SCALE = INT_AFTER/INT_BEFORE;

    for(unsigned uMom=0; uMom<hSmear[uProj]->GetNbinsX(); uMom++){
      if(hSmear[uProj]->GetBinContent(uMom+1)==0){
        continue;
      }
      if(hSmear[uProj]->GetBinContent(uMom+1)<val_min){
        val_min = hSmear[uProj]->GetBinContent(uMom+1);
      }
      if(hSmear[uProj]->GetBinContent(uMom+1)>val_max){
        val_max = hSmear[uProj]->GetBinContent(uMom+1);
      }
    }



    TH1F* hAxis = new TH1F("hAxis", "hAxis", TMath::Nint((kstar_max-kstar_min)/4.), kstar_min, kstar_max*1.3);
    hAxis->SetTitle("; #it{k*}_{true} (MeV/#it{c}); d#it{N}/d#it{k*}_{true} (#it{c}/MeV)");
    hAxis->GetXaxis()->SetRangeUser(kstar_min, kstar_max*1.3);
    hAxis->GetXaxis()->SetNdivisions(505);
    hAxis->GetYaxis()->SetRangeUser(val_min*0.5, val_max*2.0);
    hAxis->GetYaxis()->SetTitleOffset(1.4);
    hAxis->SetMarkerSize(1.6);
    hAxis->SetLineWidth(2);
    hAxis->SetMarkerStyle(kOpenCircle);
    hAxis->SetMarkerColor(kBlack);
    hAxis->SetLineColor(kGray+1);
    hAxis->SetLineWidth(0);

    hSmear[uProj]->GetXaxis()->SetRangeUser(kstar_min, kstar_max*1.3);
    hSmear[uProj]->GetXaxis()->SetNdivisions(505);
    hSmear[uProj]->GetYaxis()->SetRangeUser(val_min*0.8, val_max*1.2);
    hSmear[uProj]->SetMarkerSize(1.6);
    hSmear[uProj]->SetLineWidth(2);
    hSmear[uProj]->SetMarkerStyle(kOpenCircle);
    hSmear[uProj]->SetMarkerColor(kBlack);
    hSmear[uProj]->SetLineColor(kBlack);
    hSmear[uProj]->SetLineWidth(2);




    TCanvas* cT_PP = new TCanvas("cT_PP", "cT_PP", 0, 0, 650, 650);
    cT_PP->SetMargin(0.16,0.025,0.13,0.025);//lrbt
    cT_PP->SetLogy();

    hAxis->Draw("axis");
    hSmear[uProj]->Draw("E0 X0 same");

    TLatex BeamText;
    BeamText.SetTextSize(gStyle->GetTextSize() * 0.65);
    BeamText.SetNDC(kTRUE);
    BeamText.DrawLatex(0.55, 0.9, Form("ALICE %s #sqrt{#it{s}} = %i TeV", "pp", (int) 13));
    BeamText.DrawLatex(0.55, 0.86, "High-mult. (0#minus0.17% INEL > 0)");
    BeamText.DrawLatex(0.55, 0.82, "Momentum resolution");
    BeamText.DrawLatex(0.55, 0.78, "Simulation for pp pairs");
    TLegend *legendRESO = new TLegend(0.54,0.825-0.055-0.060*1,0.94,0.825-0.055);//lbrt
    legendRESO->SetBorderSize(0);
    legendRESO->SetTextFont(42);
    legendRESO->SetTextSize(gStyle->GetTextSize()*0.65);
    legendRESO->AddEntry(hSmear[uProj],TString::Format("#it{P}(#it{k*}=%.0f MeV/#it{c})",kstar),"fpe");
    legendRESO->Draw("same");

    cT_PP->SaveAs(TString::Format("%s/SourceStudies/source_pp_syst/cT_PP_%.0f.pdf",GetFemtoOutputFolder(),kstar));
    delete hAxis;
    delete legendRESO;
    delete cT_PP;



    TString HepFileName_RESO = TString::Format("%s/SourceStudies/source_pp_syst/HEP_RESO_%.0f_PP.yaml",GetFemtoOutputFolder(),kstar_proj.at(uProj));
    ofstream hepfile_RESO (HepFileName_RESO.Data(),ios::out);
    hepfile_RESO << "dependent_variables:" << endl;
    hepfile_RESO << "- header: {name: proton-proton Momentum resolution @ k*="<<kstar_proj.at(uProj)<<" MeV}" << endl;
    hepfile_RESO << "  qualifiers:" << endl;
    hepfile_RESO << "  - {name: SQRT(S), units: GeV, value: '13000.0'}" << endl;
    hepfile_RESO << "  values:" << endl;
    double hep_x,hep_y,hep_stat;
    for(unsigned uBin=1; uBin<=hSmear[0]->GetNbinsX(); uBin++){
      hep_x = hSmear[uProj]->GetBinCenter(uBin);
      hep_y = hSmear[uProj]->GetBinContent(uBin);
      if(hep_y<10.*SCALE) continue;
      hep_stat = hSmear[uProj]->GetBinError(uBin);
      hepfile_RESO << "  - errors:" << endl;
      hepfile_RESO << "    - {label: stat, symerror: "<<hep_stat<<"}" << endl;
      hepfile_RESO << "    value: "<<hep_y << endl;
    }
    hepfile_RESO << "independent_variables:" << endl;
    hepfile_RESO << "- header: {name: k*_true (MeV/c)}" << endl;
    hepfile_RESO << "  values:" << endl;
    for(unsigned uBin=1; uBin<hSmear[uProj]->GetNbinsX(); uBin++){
        hep_x = hSmear[uProj]->GetBinCenter(uBin);
        hep_y = hSmear[uProj]->GetBinContent(uBin);
        if(hep_y<10.*SCALE) continue;
        hepfile_RESO << "  - {value: "<<hep_x<<"}" << endl;
    }
    hepfile_RESO.close();
  }



  fOutput.cd();
  hReso->Write();

  for(unsigned uProj=0; uProj<kstar_proj.size(); uProj++){
    hSmear[uProj]->Write();
  }

  for(unsigned uProj=0; uProj<kstar_proj.size(); uProj++){
    delete hSmear[uProj];
  }
  delete [] hSmear;

  delete hReso;
}




void source_pL_syst(){

  std::vector<float> bin_range = {1.08, 1.26, 1.32, 1.44, 1.65, 1.9, 4.5 };
  std::vector<float> bin_center = {1.21, 1.29, 1.38, 1.54, 1.76, 2.26};

  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data());

  const unsigned NumDataVar = 43;
  const unsigned RebinFactor = 1;

  TFile fOutput(TString::Format("%s/SourceStudies/source_pL_syst/fOutput.root",GetFemtoOutputFolder()),"recreate");

  for(unsigned uMt=0; uMt<bin_center.size(); uMt++){

    TH2F* hSystOnly = NULL;

    double* CkMin = NULL;
    double* CkMax = NULL;

    for(unsigned uDV=0; uDV<NumDataVar; uDV++){
      TH1F* hData = AnalysisObject.GetAliceExpCorrFun("pp13TeV_HM_BernieSource","pLambda",TString::Format("%u",uDV),0,0,uMt);
      hData->Rebin(RebinFactor);
      if(!hSystOnly){
        unsigned NumBins = hData->GetNbinsX();
        double LowEdge = hData->GetXaxis()->GetBinLowEdge(1);
        double UpEdge = hData->GetXaxis()->GetBinUpEdge(NumBins);
        hSystOnly = new TH2F(TString::Format("hSystOnly_%u",uMt),TString::Format("hSystOnly_%u",uMt),NumBins,LowEdge,UpEdge,65536,0,6);
        CkMin = new double[NumBins];
        CkMax = new double[NumBins];
        for(unsigned uMom=0; uMom<hData->GetNbinsX(); uMom++){
          CkMin[uMom] = 1000;
          CkMax[uMom] = 0;
        }
      }

      for(unsigned uMom=0; uMom<hData->GetNbinsX(); uMom++){
        double MOM = hData->GetBinCenter(uMom+1);
        double CkVal = hData->GetBinContent(uMom+1);
        if(CkVal<CkMin[uMom]){CkMin[uMom] = CkVal;}
        if(CkVal>CkMax[uMom]){CkMax[uMom] = CkVal;}
        hSystOnly->Fill(MOM,CkVal);
      }
      delete hData;
    }//uDV

    TH1F* hDataReb = AnalysisObject.GetAliceExpCorrFun("pp13TeV_HM_BernieSource","pLambda",TString::Format("%u",0),0,0,uMt);
    hDataReb->SetName(TString::Format("hReb_mT%u",uMt+1));
    hDataReb->Rebin(RebinFactor);
    //TH1F* hSyst = (TH1F*)hDataReb->Clone(TString::Format("hSyst_mT%u",uMt+1));
    TH1F* hRelSystErr = (TH1F*)hDataReb->Clone(TString::Format("hRelSystErr_mT%u",uMt+1));

    TH1F* hDataStat = AnalysisObject.GetAliceExpCorrFun("pp13TeV_HM_BernieSource","pLambda",TString::Format("%u",0),0,0,uMt);
    hDataStat->SetName(TString::Format("hStat_mT%u",uMt+1));

    TH1F* hDataSyst = AnalysisObject.GetAliceExpCorrFun("pp13TeV_HM_BernieSource","pLambda",TString::Format("%u",0),0,0,uMt);
    hDataSyst->SetName(TString::Format("hSyst_mT%u",uMt+1));

    TList* list1_tmp;
    TList* list2_tmp;

    TString SystFileName = TString::Format("%s/CatsFiles/ExpData/Bernie_Source/pLData/mTBin_%i/CFOutput_mT_pLVar0_HM_%i.root",GetCernBoxDimi(), uMt+1, uMt);
    TFile* FileROOT = new TFile(SystFileName, "read");
    //FileROOT->ls();

    list1_tmp = (TList*)FileROOT->Get("PairDist");
    list2_tmp = (TList*)list1_tmp->FindObject("PairRebinned");
    //printf("%p\n",list2_tmp);
    TH1F* hMixedEvent_PP = (TH1F*)list2_tmp->FindObject(TString::Format("MEPart_mT_%u_FixShifted_Rebinned_3",uMt));

    list1_tmp = (TList*)FileROOT->Get("AntiPairDist");
    list2_tmp = (TList*)list1_tmp->FindObject("PairRebinned");

    TH1F* hMixedEvent_APAP = (TH1F*)list2_tmp->FindObject(TString::Format("MEAntiPart_mT_%u_FixShifted_Rebinned_3",uMt));
    gROOT->cd();
    TH1F* hMixedEvent = (TH1F*)hMixedEvent_PP->Clone(TString::Format("hMixedEvent_mT%u",uMt+1));
    hMixedEvent->Add(hMixedEvent_APAP);
    //printf("HI\n");
    delete FileROOT;


    TGraphErrors gDataSyst;
    TGraphErrors gDataStat;

    for(unsigned uMom=0; uMom<hDataReb->GetNbinsX(); uMom++){
      double ERROR = ( (CkMax[uMom]-CkMin[uMom])/sqrt(12.) );
      double CKVAL = hDataReb->GetBinContent(uMom+1);

      hRelSystErr->SetBinContent(uMom+1,ERROR/CKVAL);
      hRelSystErr->SetBinError(uMom+1,ERROR/CKVAL*0.1);
    }//uMom

    TF1* fRelSystErr = new TF1(TString::Format("fRelSystErr_mT%u",uMt+1),"[0]*exp(-[1]*x)+[2]",0,240);

    //1  p0           1.14917e-01   1.88225e-03   1.96102e-05   2.68707e-02
    //2  p1           1.21703e-01   1.14106e-03   9.31921e-06  -1.30492e-01
    //3  p2           2.90231e-03   2.98724e-05   4.47792e-08   2.32283e+01

    //1  p0           1.25545e-01   2.15431e-03   1.59325e-05   9.43029e-04
    //2  p1           1.35059e-01   1.32829e-03   8.23132e-06  -2.81160e-03
    //3  p2           4.07229e-03   3.44227e-05   3.76190e-08   7.74786e-01

    //1  p0           1.38030e-01   2.39855e-03   1.99977e-05   1.60832e-02
    //2  p1           1.40216e-01   1.35267e-03   9.44591e-06  -4.48411e-02
    //3  p2           4.11466e-03   3.36705e-05   4.13243e-08   3.17725e+01

    //1  p0           1.52250e-01   2.70805e-03   6.58671e-05  -1.16010e-03
    //2  p1           1.55261e-01   1.54008e-03   3.25274e-05   1.54477e-03
    //3  p2           5.17949e-03   3.81745e-05   1.25282e-07  -2.19199e-01

    fRelSystErr->SetParameter(0,1.52250e-01);
    fRelSystErr->SetParLimits(0,0.02,0.2);

    fRelSystErr->SetParameter(1,1.55261e-01);
    fRelSystErr->SetParLimits(1,0,0.2);

    fRelSystErr->SetParameter(2,hRelSystErr->GetBinContent(hRelSystErr->FindBin(200)));
    fRelSystErr->SetParLimits(2,fRelSystErr->GetParameter(2)*0.5,fRelSystErr->GetParameter(2)*2);


    for(unsigned uMom=0; uMom<hDataStat->GetNbinsX(); uMom++){
      double MOM = hDataStat->GetBinCenter(uMom+1);
      double CKVAL = hDataStat->GetBinContent(uMom+1);
      //hDataStat->SetErrorX(uMom+1);
      double EVAL_MOM = MOM>240?240:MOM;
      double REL_ERR = fRelSystErr->Eval(EVAL_MOM);
      double ERROR = CKVAL*REL_ERR;
      //if(uMt==1 && uMom>=1){
      if(uMom>=1){
        ERROR *= 1.2;
      }
      hDataSyst->SetBinError(uMom+1,ERROR);


      gDataStat.SetPoint(uMom,MOM,CKVAL);
      gDataStat.SetPointError(uMom,0,hDataStat->GetBinError(uMom+1));
      gDataSyst.SetPoint(uMom,MOM,CKVAL);
      gDataSyst.SetPointError(uMom,2,ERROR);


      if(uMt==3 && MOM<10){
        //printf("%f +/- %f\n",CKVAL,ERROR);
        //usleep(1500e3);
      }
    }

    SetStyle();

    hDataStat->SetTitle("; #it{k*} (MeV/#it{c}); #it{C}(#it{k*})");
    hDataStat->GetXaxis()->SetRangeUser(0, 400);
    hDataStat->GetXaxis()->SetNdivisions(505);
    hDataStat->GetYaxis()->SetRangeUser(0.9, 2.2);
    //else hDataStat->GetYaxis()->SetRangeUser(0.9, 2.2);
    hDataStat->SetFillColor(kGray+1);
    hDataStat->SetMarkerSize(1.5);
    hDataStat->SetLineWidth(2);
    hDataStat->SetMarkerStyle(kOpenCircle);
    hDataStat->SetMarkerColor(kBlack);
    hDataStat->SetLineColor(kBlack);

    TH1F* hAxis = new TH1F("hAxis", "hAxis", 405, 0, 405);
    hAxis->SetTitle("; #it{k*} (MeV/#it{c}); #it{C}(#it{k*})");
    hAxis->GetXaxis()->SetRangeUser(0, 405);
    hAxis->GetXaxis()->SetNdivisions(505);
    hAxis->GetYaxis()->SetRangeUser(0.9, 2.2);
    hAxis->SetFillColor(kGray+1);
    hAxis->SetMarkerSize(1.6);
    hAxis->SetLineWidth(2);
    hAxis->SetMarkerStyle(kOpenCircle);
    hAxis->SetMarkerColor(kBlack);
    hAxis->SetLineColor(kGray+1);
    hAxis->SetLineWidth(0);




    gDataSyst.SetName("gSystError");
    gDataSyst.SetLineColor(kWhite);
    gDataSyst.SetFillColorAlpha(kBlack, 0.4);

    TLegend *legend = new TLegend(0.47,0.785-0.055-0.060*1,0.94,0.785-0.055);//lbrt
    legend->SetBorderSize(0);
    legend->SetTextFont(42);
    legend->SetTextSize(gStyle->GetTextSize()*0.8);
    legend->AddEntry(hAxis,TString::Format("p#minus#Lambda #oplus #bar{p}#minus#bar{#Lambda}"),"fpe");

    TCanvas* c_PL = new TCanvas("CFpL", "CFpL", 0, 0, 650, 650);
    c_PL->SetMargin(0.13,0.025,0.12,0.025);//lrbt

    hAxis->Draw("axis");
    hDataStat->Draw("E0 X0 same");
    gDataSyst.Draw("2 same");
    legend->Draw("same");

    TLatex BeamText;
    TLatex text;
    BeamText.SetTextSize(gStyle->GetTextSize() * 0.75);
    BeamText.SetNDC(kTRUE);
    //BeamText.DrawLatex(0.32, 0.91, Form("#bf{ALICE}"));
    BeamText.DrawLatex(0.49, 0.9,
                       Form("ALICE %s #sqrt{#it{s}} = %i TeV", "pp", (int) 13));
    BeamText.DrawLatex(0.49, 0.85, "High-mult. (0#minus0.17% INEL > 0)");
    text.SetNDC();
    text.SetTextColor(1);
    text.SetTextSize(gStyle->GetTextSize() * 0.75);
    text.DrawLatex(0.49, 0.80, TString::Format("#it{m}_{T}#in [%.2f, %.2f) GeV",bin_range.at(uMt),bin_range.at(uMt+1)));
    text.DrawLatex(0.49, 0.75, TString::Format("<#it{m}_{T}> = %.2f GeV",bin_center.at(uMt)));

    c_PL->SaveAs(TString::Format("%s/SourceStudies/source_pL_syst/c_PL_%u.pdf",GetFemtoOutputFolder(),uMt));






//pL_fig1.yaml

    TString HepFileName = TString::Format("%s/SourceStudies/source_pL_syst/HEP_PL_%u.yaml",GetFemtoOutputFolder(),uMt);
    ofstream hepfile (HepFileName.Data(),ios::out);
    hepfile << "dependent_variables:" << endl;
    hepfile << "- header: {name: proton-Lambda C(k*) @ <mT>="<< TMath::Nint(bin_center.at(uMt)*1000.) <<" MeV}" << endl;
    hepfile << "  qualifiers:" << endl;
    hepfile << "  - {name: SQRT(S), units: GeV, value: '13000.0'}" << endl;
    hepfile << "  values:" << endl;
    double hep_syst;
    double hep_stat;
    double hep_x;
    double hep_y;
    for(unsigned uBin=1; uBin<hDataStat->GetNbinsX(); uBin++){
      hep_x = hDataStat->GetBinCenter(uBin);
      hep_y = hDataStat->GetBinContent(uBin);
      if(hep_x>400) break;
      hep_syst = gDataSyst.GetErrorY(uBin-1);
      hep_stat = hDataStat->GetBinError(uBin);
      hepfile << "  - errors:" << endl;
      hepfile << "    - {label: stat, symerror: "<<hep_stat<<"}" << endl;
      hepfile << "    - {label: sys, symerror: "<<hep_syst<<"}" << endl;
      hepfile << "    value: "<<hep_y << endl;
    }
    hepfile << "independent_variables:" << endl;
    hepfile << "- header: {name: k* (MeV/c)}" << endl;
    hepfile << "  values:" << endl;
    for(unsigned uBin=1; uBin<hDataStat->GetNbinsX(); uBin++){
        hep_x = hDataStat->GetBinCenter(uBin);
        if(hep_x>400) break;
        hepfile << "  - {value: "<<hep_x<<"}" << endl;
    }
    hepfile.close();


//MIXED EVENTS

    hMixedEvent->SetTitle("; #it{k*} (MeV/#it{c}); #it{M}(#it{k*})");
    hMixedEvent->GetXaxis()->SetRangeUser(0, 400);
    hMixedEvent->GetXaxis()->SetNdivisions(505);
    hMixedEvent->GetYaxis()->SetRangeUser(0.9*1000, 20000*1000);
    //hMixedEvent->SetFillColor(kGray+1);
    hMixedEvent->SetMarkerSize(1.5);
    hMixedEvent->SetLineWidth(2);
    hMixedEvent->SetMarkerStyle(kOpenCircle);
    hMixedEvent->SetMarkerColor(kBlack);
    hMixedEvent->SetLineColor(kBlack);
    hMixedEvent->GetXaxis()->SetLimits(hMixedEvent->GetXaxis()->GetXmin()*1000.,hMixedEvent->GetXaxis()->GetXmax()*1000.);

    //hMixedEvent->Scale(1e-3);

    TLegend *legendME = new TLegend(0.47,0.385-0.055-0.060*1,0.94,0.385-0.055);//lbrt
    legendME->SetBorderSize(0);
    legendME->SetTextFont(42);
    legendME->SetTextSize(gStyle->GetTextSize()*0.8);
    legendME->AddEntry(hMixedEvent,TString::Format("p#minus#Lambda #oplus #bar{p}#minus#bar{#Lambda}"),"p");

    TCanvas* cME_PP = new TCanvas("cME_PP", "cME_PP", 0, 0, 650, 650);
    cME_PP->SetMargin(0.13,0.025,0.12,0.025);//lrbt
    cME_PP->SetLogy();

    hAxis->SetTitle("; #it{k*} (MeV/#it{c}); #it{M}(#it{k*})");
    hAxis->GetXaxis()->SetRangeUser(0, 405);
    hAxis->GetXaxis()->SetNdivisions(505);
    hAxis->GetYaxis()->SetRangeUser(0.9*1000, 20000*1000);

    hAxis->Draw("axis");
    hMixedEvent->Draw("E0 X0 same");
    //gDataSyst.Draw("2 same");
    legendME->Draw("same");

    BeamText.DrawLatex(0.49, 0.5,
                       Form("ALICE %s #sqrt{#it{s}} = %i TeV", "pp", (int) 13));
    BeamText.DrawLatex(0.49, 0.45, "High-mult. (0#minus0.17% INEL > 0)");
    text.DrawLatex(0.49, 0.40, TString::Format("#it{m}_{T}#in [%.2f, %.2f) GeV",bin_range.at(uMt),bin_range.at(uMt+1)));
    text.DrawLatex(0.49, 0.35, TString::Format("<#it{m}_{T}> = %.2f GeV",bin_center.at(uMt)));

    cME_PP->SaveAs(TString::Format("%s/SourceStudies/source_pL_syst/cME_PL_%u.pdf",GetFemtoOutputFolder(),uMt));


    TString HepFileName_ME = TString::Format("%s/SourceStudies/source_pL_syst/HEP_ME_PL_%u.yaml",GetFemtoOutputFolder(),uMt);
    ofstream hepfile_ME (HepFileName_ME.Data(),ios::out);
    hepfile_ME << "dependent_variables:" << endl;
    hepfile_ME << "- header: {name: proton-Lambda M(k*) @ <mT>="<< TMath::Nint(bin_center.at(uMt)*1000.) <<" MeV}" << endl;
    hepfile_ME << "  qualifiers:" << endl;
    hepfile_ME << "  - {name: SQRT(S), units: GeV, value: '13000.0'}" << endl;
    hepfile_ME << "  values:" << endl;
    for(unsigned uBin=1; uBin<hMixedEvent->GetNbinsX(); uBin++){
      hep_x = hMixedEvent->GetBinCenter(uBin);
      hep_y = hMixedEvent->GetBinContent(uBin);
      if(hep_x>400) break;
      hep_syst = 0;
      hep_stat = sqrt(hep_y);
      hepfile_ME << "  - errors:" << endl;
      hepfile_ME << "    - {label: stat, symerror: "<<hep_stat<<"}" << endl;
      hepfile_ME << "    value: "<<hep_y << endl;
    }
    hepfile_ME << "independent_variables:" << endl;
    hepfile_ME << "- header: {name: k* (MeV/c)}" << endl;
    hepfile_ME << "  values:" << endl;
    for(unsigned uBin=1; uBin<hMixedEvent->GetNbinsX(); uBin++){
        hep_x = hMixedEvent->GetBinCenter(uBin);
        if(hep_x>400) break;
        hepfile_ME << "  - {value: "<<hep_x<<"}" << endl;
    }
    hepfile_ME.close();




    fOutput.cd();
    hSystOnly->Write();
    hDataReb->Write();
    //hSyst->Write();
    hRelSystErr->Write();
    fRelSystErr->Write();
    hDataStat->Write();
    hDataSyst->Write();
    hMixedEvent->Write();
    delete [] CkMin;
    delete [] CkMax;
    delete hSystOnly;
    delete hDataReb;
    delete hDataStat;
    delete hDataSyst;
    //delete hSyst;
    delete hRelSystErr;
    delete fRelSystErr;

    delete hAxis;
    delete legend;
    delete legendME;
    delete c_PL;

    delete hMixedEvent;
  }//uMt

  //MOMENTUM SMEAR
  TH2F* hReso = AnalysisObject.GetResolutionMatrix("pp13TeV_HM_DimiJun20","pLambda");
  std::vector<float> kstar_proj = {10, 20, 50, 100, 200, 400};
  TH1F** hSmear = new TH1F* [kstar_proj.size()];

  fOutput.cd();
  for(unsigned uProj=0; uProj<kstar_proj.size(); uProj++){
    double kstar = kstar_proj.at(uProj);
    unsigned WhichBin = hReso->GetYaxis()->FindBin(kstar);
    hSmear[uProj] = (TH1F*)hReso->ProjectionX(TString::Format("hTransform_k%.0f",kstar),WhichBin,WhichBin);
    hSmear[uProj]->Rebin(4);
    double kstar_min = 0;
    double kstar_max = 0;
    double val_min = 1e6;
    double val_max = 0;
    for(unsigned uMom=0; uMom<hSmear[uProj]->GetNbinsX(); uMom++){
      if(hSmear[uProj]->GetBinContent(uMom+1)<10){
        hSmear[uProj]->SetBinContent(uMom+1,0);
        hSmear[uProj]->SetBinError(uMom+1,0);
        continue;
      }
      hSmear[uProj]->SetBinError(uMom+1,sqrt(hSmear[uProj]->GetBinContent(uMom+1)));
      if(kstar_min==0) kstar_min = hSmear[uProj]->GetXaxis()->GetBinLowEdge(uMom+1);
      kstar_max = hSmear[uProj]->GetXaxis()->GetBinUpEdge(uMom+1);
    }

    double INT_BEFORE = hSmear[uProj]->Integral();
    hSmear[uProj]->Scale(1./hSmear[uProj]->Integral(),"width");
    double INT_AFTER = hSmear[uProj]->Integral();
    double SCALE = INT_AFTER/INT_BEFORE;

    for(unsigned uMom=0; uMom<hSmear[uProj]->GetNbinsX(); uMom++){
      if(hSmear[uProj]->GetBinContent(uMom+1)==0){
        continue;
      }
      if(hSmear[uProj]->GetBinContent(uMom+1)<val_min){
        val_min = hSmear[uProj]->GetBinContent(uMom+1);
      }
      if(hSmear[uProj]->GetBinContent(uMom+1)>val_max){
        val_max = hSmear[uProj]->GetBinContent(uMom+1);
      }
    }



    TH1F* hAxis = new TH1F("hAxis", "hAxis", TMath::Nint((kstar_max-kstar_min)/4.), kstar_min, kstar_max*1.3);
    hAxis->SetTitle("; #it{k*}_{true} (MeV/#it{c}); d#it{N}/d#it{k*}_{true} (#it{c}/MeV)");
    hAxis->GetXaxis()->SetRangeUser(kstar_min, kstar_max*1.3);
    hAxis->GetXaxis()->SetNdivisions(505);
    hAxis->GetYaxis()->SetRangeUser(val_min*0.5, val_max*2.0);
    hAxis->GetYaxis()->SetTitleOffset(1.4);
    hAxis->SetMarkerSize(1.6);
    hAxis->SetLineWidth(2);
    hAxis->SetMarkerStyle(kOpenCircle);
    hAxis->SetMarkerColor(kBlack);
    hAxis->SetLineColor(kGray+1);
    hAxis->SetLineWidth(0);

    hSmear[uProj]->GetXaxis()->SetRangeUser(kstar_min, kstar_max*1.3);
    hSmear[uProj]->GetXaxis()->SetNdivisions(505);
    hSmear[uProj]->GetYaxis()->SetRangeUser(val_min*0.8, val_max*1.2);
    hSmear[uProj]->SetMarkerSize(1.6);
    hSmear[uProj]->SetLineWidth(2);
    hSmear[uProj]->SetMarkerStyle(kOpenCircle);
    hSmear[uProj]->SetMarkerColor(kBlack);
    hSmear[uProj]->SetLineColor(kBlack);
    hSmear[uProj]->SetLineWidth(2);




    TCanvas* cT_PL = new TCanvas("cT_PL", "cT_PL", 0, 0, 650, 650);
    cT_PL->SetMargin(0.16,0.025,0.13,0.025);//lrbt
    cT_PL->SetLogy();

    hAxis->Draw("axis");
    hSmear[uProj]->Draw("E0 X0 same");

    TLatex BeamText;
    BeamText.SetTextSize(gStyle->GetTextSize() * 0.65);
    BeamText.SetNDC(kTRUE);
    BeamText.DrawLatex(0.55, 0.9, Form("ALICE %s #sqrt{#it{s}} = %i TeV", "pp", (int) 13));
    BeamText.DrawLatex(0.55, 0.86, "High-mult. (0#minus0.17% INEL > 0)");
    BeamText.DrawLatex(0.55, 0.82, "Momentum resolution");
    BeamText.DrawLatex(0.55, 0.78, "Simulation for p#Lambda pairs");
    TLegend *legendRESO = new TLegend(0.54,0.825-0.055-0.060*1,0.94,0.825-0.055);//lbrt
    legendRESO->SetBorderSize(0);
    legendRESO->SetTextFont(42);
    legendRESO->SetTextSize(gStyle->GetTextSize()*0.65);
    legendRESO->AddEntry(hSmear[uProj],TString::Format("#it{P}(#it{k*}=%.0f MeV/#it{c})",kstar),"fpe");
    legendRESO->Draw("same");

    cT_PL->SaveAs(TString::Format("%s/SourceStudies/source_pL_syst/cT_PL_%.0f.pdf",GetFemtoOutputFolder(),kstar));
    delete hAxis;
    delete legendRESO;
    delete cT_PL;



    TString HepFileName_RESO = TString::Format("%s/SourceStudies/source_pL_syst/HEP_RESO_%.0f_PL.yaml",GetFemtoOutputFolder(),kstar_proj.at(uProj));
    ofstream hepfile_RESO (HepFileName_RESO.Data(),ios::out);
    hepfile_RESO << "dependent_variables:" << endl;
    hepfile_RESO << "- header: {name: proton-Lambda Momentum resolution @ k*="<<kstar_proj.at(uProj)<<" MeV}" << endl;
    hepfile_RESO << "  qualifiers:" << endl;
    hepfile_RESO << "  - {name: SQRT(S), units: GeV, value: '13000.0'}" << endl;
    hepfile_RESO << "  values:" << endl;
    double hep_x,hep_y,hep_stat;
    for(unsigned uBin=1; uBin<=hSmear[0]->GetNbinsX(); uBin++){
      hep_x = hSmear[uProj]->GetBinCenter(uBin);
      hep_y = hSmear[uProj]->GetBinContent(uBin);
      if(hep_y<10.*SCALE) continue;
      hep_stat = hSmear[uProj]->GetBinError(uBin);
      hepfile_RESO << "  - errors:" << endl;
      hepfile_RESO << "    - {label: stat, symerror: "<<hep_stat<<"}" << endl;
      hepfile_RESO << "    value: "<<hep_y << endl;
    }
    hepfile_RESO << "independent_variables:" << endl;
    hepfile_RESO << "- header: {name: k*_true (MeV/c)}" << endl;
    hepfile_RESO << "  values:" << endl;
    for(unsigned uBin=1; uBin<hSmear[uProj]->GetNbinsX(); uBin++){
        hep_x = hSmear[uProj]->GetBinCenter(uBin);
        hep_y = hSmear[uProj]->GetBinContent(uBin);
        if(hep_y<10.*SCALE) continue;
        hepfile_RESO << "  - {value: "<<hep_x<<"}" << endl;
    }
    hepfile_RESO.close();
  }


  fOutput.cd();
  hReso->Write();

  for(unsigned uProj=0; uProj<kstar_proj.size(); uProj++){
    hSmear[uProj]->Write();
  }

  for(unsigned uProj=0; uProj<kstar_proj.size(); uProj++){
    delete hSmear[uProj];
  }
  delete [] hSmear;

  delete hReso;


}


void estimate_syst_fit(){

//mom0: 1.94 +/- 0.12 --> RE 0.062
//mom1: 1.61 +/- 0.022 --> RE 0.014
//mom2:                 --> RE 0.007
//    3                         0.0045
//    4                         0.0040

  std::vector<float> rel_error_pL = {0.066,0.014,0.007,0.0055,0.0050};
  TH1F* hRE = new TH1F("hRE","hRE",rel_error_pL.size(),0,12.*rel_error_pL.size());
  for(unsigned uBin=0; uBin<rel_error_pL.size(); uBin++){
    hRE->SetBinContent(uBin+1,rel_error_pL.at(uBin));
    hRE->SetBinError(uBin+1,rel_error_pL.at(uBin)*0.01);
  }
  TF1* fRE = new TF1("fRE", "[0]*exp(-[1]*x)+[2]", 0,12.*rel_error_pL.size());

  fRE->SetParameter(0,0.08);
  fRE->SetParLimits(0,0.02,0.2);

  fRE->SetParameter(1,0.1);
  fRE->SetParLimits(1,0,0.2);


  hRE->Fit(fRE,"S, N, R, M");

  TFile fOutput(TString::Format("%s/SourceStudies/source_pL_syst/estimate_syst_fit.root",GetFemtoOutputFolder()),"recreate");
  hRE->Write();
  fRE->Write();
}


void SlopeOfMt_1(){

  //should represent more or less the distances between two
  //neighbouring particles on the hadronization surface.
  const double RH = 0.5;

  //should represent more or less the distances between two
  //neighbouring particles at freeze-out
  const double RF = 1.2;

  //beta*gamma
  const double BG = 1;

  const double h_base = 2.7;
  const double MN = h_base / RH;
  const double M_base = 48.36;

  std::vector<double> mult_val;
  mult_val.push_back(13.75);
  mult_val.push_back(25.48);
  mult_val.push_back(33.16);
  mult_val.push_back(48.36);

  for(unsigned uMult=0; uMult<mult_val.size(); uMult++){
    double MR = sqrt(mult_val.at(uMult)/M_base);
    double h_par = RH*MN*MR;
    double t_par = MN*MR*RF-h_par;
    printf("h = %.2f\n",h_par);
    printf("t = %.2f\n",t_par);
    printf("----------------\n");
  }

}


int SOURCESTUDIES(int argc, char *argv[]){
    //CompareCkAndSr();
    //Compare_2_vs_3_body();

    SlopeOfMt_1(); return 0;

    //ConvertThetaAngleHisto("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/DimiPhi_pp.root","h_rkAngle_Mom2",400,600);
    //AverageResoApprox_pp();
    //AverageResoApprox_pL();
    //TestDifferentAngularDistributions(0.8);
    //TestDifferentAngularDistributions(1.0);
    //TestDifferentAngularDistributions(1.2);

    //TestEposDistos(1.0);

    //void SourcePaper_pp(const TString SourceVar, const int& SmearStrategy, const unsigned DataVar, const int imTbin, const TString OutputFolder)
    //SourcePaper_pp(argv[1],atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),argv[5]);
    //SourcePaper_pp("McGauss_ResoTM",1,0,0,"/home/dimihayl/Software/LocalFemto/Output/SourceStudies/SourcePaper_pp/");
    //SourcePaper_pp("McGauss_ResoTM",0,0,0,"/home/dimihayl/Software/LocalFemto/Output/SourceStudies/SourcePaper_pp/");
    //SourcePaper_pL(argv[1],atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),argv[5]);
    //SourcePaper_pL("McGauss_ResoTM",1,0,0,"/home/dimihayl/Software/LocalFemto/Output/SourceStudies/SourcePaper_pL/");
    //SourcePaper_pp("Gauss",1,0,0,"/home/dimihayl/Software/LocalFemto/Output/SourceStudies/SourcePaper_pp/");
    //SourcePaper_pp("McGauss_ResoTM",1,0,0,"/home/dimihayl/Software/LocalFemto/Output/SourceStudies/SourcePaper_pp/");

    //TestReadTTree();
    //Estimate_Reff("pLambda","oton","RSM_PLB",3.0);
    //Estimate_Reff("pLambda","oton","RSM_GHETTO",3.0);
    //Estimate_Reff("pLambda","oton","RSM_FLAT",3.0);
    //printf("------------------\n");

    //Estimate_Reff("pXi","oton","RSM_PLB",3.0);
    //Estimate_Reff("pXi","oton","RSM_GHETTO",3.0);
    //Estimate_Reff("pXi","oton","RSM_FLAT",3.0);
    //printf("------------------\n");

    //Estimate_Reff("pOmega","oton","RSM_PLB",3.0);
    //Estimate_Reff("pOmega","oton","RSM_GHETTO",3.0);
    //Estimate_Reff("pOmega","oton","RSM_FLAT",3.0);

    //Estimate_Reff("pd_min","oton","RSM_PLB",3.0);
    //Estimate_Reff("pd","oton","RSM_PLB",3.0);
    //Estimate_Reff("pd_max","oton","RSM_PLB",3.0);

    //source_pp_syst();
    //source_pL_syst();
    //estimate_syst_fit();

    //SourceDensity();

    //printf("------------------\n");

    //SinglePart_RandomLab();
    //int RndSeed, int imTbin, int pL_inter, int AnalysisType, int ErrorType, double Timeout, int NewAgeFlag
    ////ErrorType = XY, where X is syst vars, Y is bootstrap
    //Jaime_pL(1, 0, 11600, 1101, 0, -0.1, 0, TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data(), "/home/dimihayl/Software/LocalFemto/Output/SourceStudies/Jaime_pL/");
    //Jaime_pL(1, 0, 11600, 1201, 10, 3, 0, TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data(), "/home/dimihayl/Software/LocalFemto/Output/SourceStudies/Jaime_pL/tmp/");
    //Jaime_pL(1, 0, 0, 0, 0, -0.1, 0, TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data(), "/home/dimihayl/Software/LocalFemto/Output/SourceStudies/Jaime_pL/");
////r,0,11600/131600,12(3)01,11,10,2
    //Jaime_pL(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7]), TString::Format("%s/CatsFiles",GetCernBoxDimi()).Data(), "/home/dimihayl/Software/LocalFemto/Output/SourceStudies/Jaime_pL/310522/");
    //SourcePaper_Published(TString::Format("%s/Jaime",GetCernBoxDimi()).Data());
    //1, 0-6, 0, (0)(0/1)(0)(0) 11 15 0 INPUT OUTPUT
    printf("Terminating\n");
    return 0;
}
