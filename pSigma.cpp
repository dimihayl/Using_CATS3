#include "pSigma.h"

#include "CATS.h"
#include "DLM_Ck.h"


#include "CommonAnaFunctions.h"

#include "TGraph.h"
#include "TFile.h"

void TestExternalWF(){

    const unsigned NumMomBins = 35;
    const double kMin = 0;
    const double kMax = 350;

    CATS Cat_ESC16;
    CATS Cat_CHIRAL;

    Cat_ESC16.SetMomBins(NumMomBins,kMin,kMax);
    Cat_CHIRAL.SetMomBins(NumMomBins,kMin,kMax);

    DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder("/home/dmihaylov/CernBox/CatsFiles");
    AnalysisObject.SetUpCats_pS0(Cat_ESC16,"ESC16","Gauss");
    AnalysisObject.SetUpCats_pS0(Cat_CHIRAL,"Chiral","Gauss");

    Cat_ESC16.KillTheCat();
    Cat_CHIRAL.KillTheCat();

    TGraph gESC16;
    gESC16.SetName("gESC16");
    gESC16.SetLineColor(kTeal+2);
    gESC16.SetLineWidth(3);

    TGraph gCHIRAL;
    gCHIRAL.SetName("gCHIRAL");
    gCHIRAL.SetLineColor(kAzure-8);
    gCHIRAL.SetLineWidth(3);
    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        gESC16.SetPoint(uBin,Cat_ESC16.GetMomentum(uBin),Cat_ESC16.GetCorrFun(uBin));
        gCHIRAL.SetPoint(uBin,Cat_CHIRAL.GetMomentum(uBin),Cat_CHIRAL.GetCorrFun(uBin));
    }

    TFile fOutput("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pSigma/TestExternalWF.root","recreate");
    gESC16.Write();
    gCHIRAL.Write();

}

void BugIn_pL_Fit_Test(){

    const unsigned NumBinsCk = 45; // for dlm_ck
    const double BinWidthCk = 12;
    const double MaxBinValCk = double(NumBinsCk)*BinWidthCk;
    unsigned NumMomBins_pL = 30;// for CATS
    const double MaxBinValCats = double(NumMomBins_pL)*BinWidthCk;

    const double ResidualSourceSizeSigma0 = 1.25;

    DLM_CommonAnaFunctions AnalysisObject; AnalysisObject.SetCatsFilesFolder("/home/dmihaylov/CernBox/CatsFiles");

    CATS AB_pS0_Chiral;
    //the minus one is to to avoid going above 350 MeV, since we do not have the WF there
    AB_pS0_Chiral.SetMomBins(NumBinsCk-1,0,MaxBinValCats-BinWidthCk);
    AnalysisObject.SetUpCats_pS0(AB_pS0_Chiral,"Chiral","Gauss");
    AB_pS0_Chiral.SetAnaSource(0,ResidualSourceSizeSigma0);
    AB_pS0_Chiral.SetNotifications(CATS::nError);
    AB_pS0_Chiral.KillTheCat();

    CATS AB_pS0_ESC16;
    AB_pS0_ESC16.SetMomBins(NumBinsCk-1,0,MaxBinValCats-BinWidthCk);
    AnalysisObject.SetUpCats_pS0(AB_pS0_ESC16,"ESC16","Gauss");
    AB_pS0_ESC16.SetAnaSource(0,ResidualSourceSizeSigma0);
    AB_pS0_ESC16.SetNotifications(CATS::nError);
    AB_pS0_ESC16.KillTheCat();

    DLM_Ck* Ck_pS0_Chiral = new DLM_Ck(AB_pS0_Chiral.GetNumSourcePars(),0,AB_pS0_Chiral,NumBinsCk,0,MaxBinValCk);
    Ck_pS0_Chiral->SetSourcePar(0,ResidualSourceSizeSigma0);
    Ck_pS0_Chiral->SetCutOff(MaxBinValCats-BinWidthCk,700);
    Ck_pS0_Chiral->Update();

    DLM_Ck* Ck_pS0_ESC16 = new DLM_Ck(AB_pS0_ESC16.GetNumSourcePars(),0,AB_pS0_ESC16,NumBinsCk,0,MaxBinValCk);
    Ck_pS0_ESC16->SetSourcePar(0,ResidualSourceSizeSigma0);
    Ck_pS0_ESC16->SetCutOff(MaxBinValCats-BinWidthCk,700);
    Ck_pS0_ESC16->Update();

    printf("CutOff = %f\n",MaxBinValCats-BinWidthCk);

    TGraph gESC16;
    gESC16.SetName("gESC16");
    gESC16.SetLineColor(kTeal+2);
    gESC16.SetLineWidth(3);

    TGraph gCHIRAL;
    gCHIRAL.SetName("gCHIRAL");
    gCHIRAL.SetLineColor(kAzure-8);
    gCHIRAL.SetLineWidth(3);

    TGraph gcESC16;
    gcESC16.SetName("gcESC16");
    gcESC16.SetLineColor(kTeal+2);
    gcESC16.SetLineWidth(3);

    TGraph gcCHIRAL;
    gcCHIRAL.SetName("gcCHIRAL");
    gcCHIRAL.SetLineColor(kAzure-8);
    gcCHIRAL.SetLineWidth(3);

    for(unsigned uBin=0; uBin<NumBinsCk; uBin++){
        gESC16.SetPoint(uBin,Ck_pS0_ESC16->GetBinCenter(0,uBin),Ck_pS0_ESC16->GetBinContent(uBin));
        gCHIRAL.SetPoint(uBin,Ck_pS0_Chiral->GetBinCenter(0,uBin),Ck_pS0_Chiral->GetBinContent(uBin));
    }
    for(unsigned uBin=0; uBin<NumBinsCk-1; uBin++){
        gcESC16.SetPoint(uBin,AB_pS0_ESC16.GetMomentum(uBin),AB_pS0_ESC16.GetCorrFun(uBin));
        gcCHIRAL.SetPoint(uBin,AB_pS0_Chiral.GetMomentum(uBin),AB_pS0_Chiral.GetCorrFun(uBin));
    }

    TFile fOutput("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pSigma/BugIn_pL_Fit_Test.root","recreate");
    gESC16.Write();
    gCHIRAL.Write();
    gcESC16.Write();
    gcCHIRAL.Write();
}



void Main_pSigma(){

    //TestExternalWF();
    BugIn_pL_Fit_Test();

}
