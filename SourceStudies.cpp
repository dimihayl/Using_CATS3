
#include "SourceStudies.h"
#include "CommonAnaFunctions.h"
#include "CATS.h"
#include "CATStools.h"
#include "DLM_Source.h"
#include "DLM_Potentials.h"
#include "DLM_Histo.h"
#include "DLM_Random.h"

#include "TFile.h"
#include "TGraph.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TMath.h"
#include "TNtuple.h"

//for the source paper, we compare C(k) and S(r) of pp and pL based on
//Gaussian and Gaussian+Reso sources
void CompareCkAndSr(){

    const TString OutputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/SourceStudies/CompareCkAndSr/Feb2020/";

    //const double R_core = 1.00;
    //const double R_pp = 1.32;
    //const double R_pL = 1.48;

    const double R_core = 1.20;
    const double R_pp = 1.280;
    const double R_pL = 1.414;

    const unsigned kNumBins=100;
    const double kMin = 0;
    const double kMax = 200;

    const unsigned rNumBins=100;
    const double rMin = 0;
    const double rMax = 10;
    DLM_Histo<double> rDummy;
    rDummy.SetUp(1);
    rDummy.SetUp(0,rNumBins,rMin,rMax);
    rDummy.Initialize();

    DLM_CommonAnaFunctions AnalysisObject;

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
    AnalysisObject.SetUpCats_pp(AB_pp,"AV18","McGauss_ResoTM",0,202);
    AB_pp.SetAnaSource(0,R_core);
    AB_pp.KillTheCat();
    TGraph gCk_pp_CoreReso; gCk_pp_CoreReso.SetName("gCk_pp_CoreReso"); gCk_pp_CoreReso.Set(kNumBins);
    for(unsigned uBin=0; uBin<kNumBins; uBin++){
        gCk_pp_CoreReso.SetPoint(uBin,AB_pp.GetMomentum(uBin),AB_pp.GetCorrFun(uBin));
    }
    TGraph gSr_pp_CoreReso; gSr_pp_CoreReso.SetName("gSr_pp_CoreReso"); gSr_pp_CoreReso.Set(kNumBins);
    for(unsigned uBin=0; uBin<rNumBins; uBin++){
        double RAD = rDummy.GetBinCenter(0,uBin);
        //gSr_pp_CoreReso.SetPoint(uBin,rDummy.GetBinCenter(0,uBin),AB_pp.EvaluateTheSource(0,rDummy.GetBinCenter(0,uBin),0)/(4.*TMath::Pi()*RAD*RAD));
        gSr_pp_CoreReso.SetPoint(uBin,rDummy.GetBinCenter(0,uBin),AB_pp.EvaluateTheSource(0,rDummy.GetBinCenter(0,uBin),0));
    }

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

    //AnalysisObject.SetUpCats_pL(AB_pL,"NLO_Coupled_S","McGauss_Reso");
    AnalysisObject.SetUpCats_pL(AB_pL,"NLO_Coupled_S","McGauss_ResoTM",0,202);
    AB_pL.SetAnaSource(0,R_core);
    AB_pL.KillTheCat();
    TGraph gCk_pL_CoreReso; gCk_pL_CoreReso.SetName("gCk_pL_CoreReso"); gCk_pL_CoreReso.Set(kNumBins);
    for(unsigned uBin=0; uBin<kNumBins; uBin++){
        gCk_pL_CoreReso.SetPoint(uBin,AB_pL.GetMomentum(uBin),AB_pL.GetCorrFun(uBin));
    }
    TGraph gSr_pL_CoreReso; gSr_pL_CoreReso.SetName("gSr_pL_CoreReso"); gSr_pL_CoreReso.Set(kNumBins);
    for(unsigned uBin=0; uBin<rNumBins; uBin++){
        double RAD = rDummy.GetBinCenter(0,uBin);
        //gSr_pL_CoreReso.SetPoint(uBin,rDummy.GetBinCenter(0,uBin),AB_pL.EvaluateTheSource(0,rDummy.GetBinCenter(0,uBin),0)/(4.*TMath::Pi()*RAD*RAD));
        gSr_pL_CoreReso.SetPoint(uBin,rDummy.GetBinCenter(0,uBin),AB_pL.EvaluateTheSource(0,rDummy.GetBinCenter(0,uBin),0));
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
    gSr_pp_Gauss.Write();
    gSr_pL_Core.Write();
    gSr_pL_CoreReso.Write();
    gSr_pL_Gauss.Write();

    TCanvas* cSr = new TCanvas("cSr", "cSr", 1);
    cSr->cd(0); cSr->SetCanvasSize(1920/2, 1080/2); cSr->SetMargin(0.15,0.05,0.2,0.05);//lrbt

    TLegend* legSr = new TLegend(0.65,0.65,0.975,0.975);//lbrt
    legSr->SetName(TString::Format("legSr"));
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

int SOURCESTUDIES(int narg, char** ARGS){
    CompareCkAndSr();
    //ConvertThetaAngleHisto("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/DimiPhi_pp.root","h_rkAngle_Mom2",400,600);
    //AverageResoApprox_pp();
    //AverageResoApprox_pL();
    //TestDifferentAngularDistributions(0.8);
    //TestDifferentAngularDistributions(1.0);
    //TestDifferentAngularDistributions(1.2);

    //TestEposDistos(1.0);

    return 0;
}
