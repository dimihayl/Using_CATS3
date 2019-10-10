
#include "TRandom3.h"
#include "TString.h"
#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "TCanvas.h"
#include "TStyle.h"

double relKcalc(const TLorentzVector& track1, const TLorentzVector& track2)
{

  TLorentzVector trackSum, track1_cms, track2_cms;
  trackSum = track1 + track2;

  double beta = trackSum.Beta();
  double beta_x = beta*cos(trackSum.Phi())*sin(trackSum.Theta());
  double beta_y = beta*sin(trackSum.Phi())*sin(trackSum.Theta());
  double beta_z = beta*cos(trackSum.Theta());

  track1_cms = track1;
  track2_cms = track2;

  track1_cms.Boost(-beta_x,-beta_y,-beta_z);
  track2_cms.Boost(-beta_x,-beta_y,-beta_z);

  TLorentzVector track_relK;

  track_relK = track1_cms - track2_cms;
  double relK = 0.5*track_relK.P();


  return relK;
}

//used to generate the decay matrix
void DecayMatrix(int SEED, int NumIter)
{
    // example of use of TGenPhaseSpace
    //Author: Valerio Filippini
    //further modified my O.W. Arnold and D.L. Mihaylov

    //if (!gROOT->GetClass("TGenPhaseSpace")) gSystem->Load("libPhysics");

    TRandom3 rangen(SEED);
/*
    TRandom3 *t3_ran_p_proton = new TRandom3();
    TRandom3 *t3_ran_p_Xim = new TRandom3();
    TRandom3 *t3_ran_p_Lambda = new TRandom3();
    TRandom3 *t3_ran_p_Sigma0 = new TRandom3();
    TRandom3 *t3_ran_p_Xim1530 = new TRandom3();
    TRandom3 *t3_ran_p_Omegam = new TRandom3();

    t3_ran_p_proton->SetSeed(SEED*10);
    t3_ran_p_Xim->SetSeed(SEED*10+1);
    t3_ran_p_Lambda->SetSeed(SEED*10+2);
    t3_ran_p_Sigma0->SetSeed(SEED*10+3);
    t3_ran_p_Xim1530->SetSeed(SEED*10+4);
    t3_ran_p_Omegam->SetSeed(SEED*10+5);
*/
    TString OutputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/GenBod/ForVale/";
    TString output_string = TString::Format("%sOUTFILE_%i.root",OutputFolder.Data(),SEED);
    TFile *outputfile = new TFile(output_string.Data(),"RECREATE");

    unsigned int NumBins = 1000;
    double kMax = 1000;

    TH2F* hRes_pp_pL = new TH2F("hRes_pp_pL","x: relk pV0, y: relk pp",NumBins,0,kMax,NumBins,0,kMax);
    TH2F* hRes_pL_LL = new TH2F("hRes_pL_LL","x: relk LL, y: relk pL",NumBins,0,kMax,NumBins,0,kMax);
    TH2F* hRes_pp_LL = new TH2F("hRes_pp_LL","x: relk LL, y: relk pp",NumBins,0,kMax,NumBins,0,kMax);
    TH2F* hRes_pL_pSigma0 = new TH2F("hRes_pL_pSigma0","x: relk pSigma0, y: relk pV0",NumBins,0,kMax,NumBins,0,kMax);
    TH2F* hRes_pL_pXim = new TH2F("hRes_pL_pXim","x: relk pXim, y: relk pV0",NumBins,0,kMax,NumBins,0,kMax);
    TH2F* hRes_pXim_pXim1530 = new TH2F("hRes_pXim_pXim1530","x: relk pXim1530, y: relk pXim",NumBins,0,kMax,NumBins,0,kMax);
    TH2F* hRes_pXim_pOmegam = new TH2F("hRes_pXim_pOmegam","x: relk pOmegam, y: relk pXim",NumBins,0,kMax,NumBins,0,kMax);
    TH2F* hRes_pp_ppRESO = new TH2F("hRes_pp_ppRESO","hRes_pp_ppRESO",NumBins,0,kMax,NumBins,0,kMax);
    TH2F* hRes_pp_pRESOpRESO = new TH2F("hRes_pp_pRESOpRESO","hRes_pp_pRESOpRESO",NumBins,0,kMax,NumBins,0,kMax);

    double Xim_decay_daughters[2] = {1.116,0.138};//Lambda-pion
    double Lambda_decay_daughters[2] = {0.938,0.138};//proton-pion
    double Sigma0_decay_daughters[2] = {1.116,0.};//Lambda-photon
    double Xim1530_decay_daughters[2] = {1.322,0.135};//Xim-pion0
    double Omegam_decay_daughters[2] = {1.322,0.135};//Xim-pion0
    double pRESO_decay_daughters[2] = {0.938,0.138};//pp feeddown

    const double momMean = 0.0;
    const double momSpread = 0.35;

    TLorentzVector tvec_proton,tvec_Xim,tvec_Lambda,tvec_Lambda_2,tvec_Sigma0,tvec_Xim1530,tvec_Omegam,tvec_pRESO,tvec_pRESO2;
    TLorentzVector *tvec_pp_pL,*tvec_pL_LL,*tvec_pL_pSigma0,*tvec_pL_pXim,*tvec_pXim_pXim1530,*tvec_pXim_pOmegam,*tvec_pp_ppRESO,*tvec_pp_pRESOpRESO;
    TGenPhaseSpace eventLambda,eventLambda_2,eventSigma0,eventXim,eventXim1530,eventOmegam,event_pRESO,event_pRESO2;

    for(int i=0;i<NumIter;i++){
        //if(i%1000==0) printf("i = %i(%i)\n",i,NumIter);

        double px_proton = rangen.Gaus(momMean,momSpread);
        double py_proton = rangen.Gaus(momMean,momSpread);
        double pz_proton = rangen.Gaus(momMean,momSpread);

        double px_L = rangen.Gaus(momMean,momSpread);
        double py_L = rangen.Gaus(momMean,momSpread);
        double pz_L = rangen.Gaus(momMean,momSpread);

        double px_L_2 = rangen.Gaus(momMean,momSpread);
        double py_L_2 = rangen.Gaus(momMean,momSpread);
        double pz_L_2 = rangen.Gaus(momMean,momSpread);

        double px_Sigma0 = rangen.Gaus(momMean,momSpread);
        double py_Sigma0 = rangen.Gaus(momMean,momSpread);
        double pz_Sigma0 = rangen.Gaus(momMean,momSpread);

        double px_Xim = rangen.Gaus(momMean,momSpread);
        double py_Xim = rangen.Gaus(momMean,momSpread);
        double pz_Xim = rangen.Gaus(momMean,momSpread);

        double px_Xim1530 = rangen.Gaus(momMean,momSpread);
        double py_Xim1530 = rangen.Gaus(momMean,momSpread);
        double pz_Xim1530 = rangen.Gaus(momMean,momSpread);

        double px_Omegam = rangen.Gaus(momMean,momSpread);
        double py_Omegam = rangen.Gaus(momMean,momSpread);
        double pz_Omegam = rangen.Gaus(momMean,momSpread);

        double px_ppRESO = rangen.Gaus(momMean,momSpread);
        double py_ppRESO = rangen.Gaus(momMean,momSpread);
        double pz_ppRESO = rangen.Gaus(momMean,momSpread);

        double px_ppRESO2 = rangen.Gaus(momMean,momSpread);
        double py_ppRESO2 = rangen.Gaus(momMean,momSpread);
        double pz_ppRESO2 = rangen.Gaus(momMean,momSpread);

        tvec_proton.SetXYZM(px_proton,py_proton,pz_proton,0.938);//units are in GeV
        tvec_Lambda.SetXYZM(px_L,py_L,pz_L,1.116);
        tvec_Lambda_2.SetXYZM(px_L_2,py_L_2,pz_L_2,1.116);
        tvec_Sigma0.SetXYZM(px_Sigma0,py_Sigma0,pz_Sigma0,1.193);
        tvec_Xim.SetXYZM(px_Xim,py_Xim,pz_Xim,1.322);
        tvec_Xim1530.SetXYZM(px_Xim1530,py_Xim1530,pz_Xim1530,1.535);
        tvec_Omegam.SetXYZM(px_Omegam,py_Omegam,pz_Omegam,1.672);
        tvec_pRESO.SetXYZM(px_ppRESO,py_ppRESO,pz_ppRESO,1.36152);
        tvec_pRESO2.SetXYZM(px_ppRESO2,py_ppRESO2,pz_ppRESO2,1.36152);

        //Calculate initial momentum of the pair:
        double relK_pL = relKcalc(tvec_proton,tvec_Lambda);
        double relK_LL = relKcalc(tvec_Lambda,tvec_Lambda_2);
        double relK_pSigma0 = relKcalc(tvec_proton,tvec_Sigma0);
        double relK_pXim = relKcalc(tvec_proton,tvec_Xim);
        double relK_pXim1530 = relKcalc(tvec_proton,tvec_Xim1530);
        double relK_pOmegam = relKcalc(tvec_proton,tvec_Omegam);
        double relK_ppRESO = relKcalc(tvec_proton,tvec_pRESO);
        double relK_pRESOpRESO = relKcalc(tvec_pRESO,tvec_pRESO2);

        //p-V0 -> p-p
        eventLambda.SetDecay(tvec_Lambda, 2, Lambda_decay_daughters);
        eventLambda.Generate();
        //momentum of the decay proton
        tvec_pp_pL = eventLambda.GetDecay(0);
        //momentum between decay proton and original proton
        double relK_pp_pL = relKcalc(tvec_proton,*tvec_pp_pL);
        hRes_pp_pL->Fill(relK_pL*1000,relK_pp_pL*1000);

        //V0-V0 -> p-V0
        eventLambda_2.SetDecay(tvec_Lambda_2, 2, Lambda_decay_daughters);
        eventLambda_2.Generate();
        //momentum of the decay proton
        tvec_pL_LL = eventLambda_2.GetDecay(0);
        double relK_pL_LL = relKcalc(tvec_Lambda,*tvec_pL_LL);
        hRes_pL_LL->Fill(relK_LL*1000,relK_pL_LL*1000);

        //V0-V0 -> p-p
        //momentum between the two decay protons
        double relK_pp_LL = relKcalc(*tvec_pp_pL,*tvec_pL_LL);
        hRes_pp_LL->Fill(relK_LL*1000,relK_pp_LL*1000);

        //p-Sigma0 -> p-Lambda
        eventSigma0.SetDecay(tvec_Sigma0, 2, Sigma0_decay_daughters);
        eventSigma0.Generate();
        tvec_pL_pSigma0 = eventSigma0.GetDecay(0);
        double relK_pL_pSigma0 = relKcalc(tvec_proton,*tvec_pL_pSigma0);
        hRes_pL_pSigma0->Fill(relK_pSigma0*1000,relK_pL_pSigma0*1000);

        //Simulate the decay and obtain the final momentum after the decay:
        //p-Xi -> p-V0
        eventXim.SetDecay(tvec_Xim, 2, Xim_decay_daughters);//Let the Xi decay in Lambda pion
        eventXim.Generate();
        tvec_pL_pXim = eventXim.GetDecay(0);
        double relK_pL_pXim = relKcalc(tvec_proton,*tvec_pL_pXim);
        hRes_pL_pXim->Fill(relK_pXim*1000,relK_pL_pXim*1000);

        //p-Xim1530 -> p-Xim
        eventXim1530.SetDecay(tvec_Xim1530, 2, Xim1530_decay_daughters);
        eventXim1530.Generate();
        tvec_pXim_pXim1530 = eventXim1530.GetDecay(0);
        double relK_pXim_pXim1530 = relKcalc(tvec_proton,*tvec_pXim_pXim1530);
        hRes_pXim_pXim1530->Fill(relK_pXim1530*1000,relK_pXim_pXim1530*1000);

        //p-Omegam -> p-Xim
        eventOmegam.SetDecay(tvec_Omegam, 2, Omegam_decay_daughters);
        eventOmegam.Generate();
        tvec_pXim_pOmegam = eventOmegam.GetDecay(0);
        double relK_pXim_pOmegam = relKcalc(tvec_proton,*tvec_pXim_pOmegam);
        hRes_pXim_pOmegam->Fill(relK_pOmegam*1000,relK_pXim_pOmegam*1000);

        //Reso feeddown to p
        event_pRESO.SetDecay(tvec_pRESO, 2, pRESO_decay_daughters);
        event_pRESO.Generate();
        tvec_pp_ppRESO = event_pRESO.GetDecay(0);
        double relK_pp_ppRESO = relKcalc(tvec_proton,*tvec_pp_ppRESO);
        hRes_pp_ppRESO->Fill(relK_ppRESO*1000,relK_pp_ppRESO*1000);

        //ResoReso feeddown to pp
        event_pRESO2.SetDecay(tvec_pRESO2, 2, pRESO_decay_daughters);
        event_pRESO2.Generate();
        tvec_pp_pRESOpRESO = event_pRESO2.GetDecay(0);
        double relK_pp_pRESOpRESO = relKcalc(*tvec_pp_ppRESO,*tvec_pp_pRESOpRESO);
        hRes_pp_pRESOpRESO->Fill(relK_pRESOpRESO*1000,relK_pp_pRESOpRESO*1000);
    }

    outputfile->cd();

    //for(unsigned uBinX=1; uBinX<=hRes_pp_RESO->GetNbinsX(); uBinX++)for(unsigned uBinY=1; uBinY<=hRes_pp_RESO->GetNbinsY(); uBinY++)
    //    if(hRes_pp_RESO->GetBinContent(uBinX,uBinY))hRes_pp_RESO->SetBinContent(uBinX,uBinY,1);

    hRes_pp_ppRESO->SetStats(false);
    hRes_pp_ppRESO->SetTitle("");
    hRes_pp_ppRESO->GetXaxis()->SetTitle("k* (p-Reso) (MeV)");
    hRes_pp_ppRESO->GetXaxis()->SetTitleSize(0.06);
    hRes_pp_ppRESO->GetXaxis()->SetLabelSize(0.06);
    hRes_pp_ppRESO->GetXaxis()->CenterTitle();
    hRes_pp_ppRESO->GetXaxis()->SetTitleOffset(1.3);
    hRes_pp_ppRESO->GetXaxis()->SetLabelOffset(0.02);
    hRes_pp_ppRESO->GetYaxis()->SetTitle("k* (p-p) (MeV)");
    hRes_pp_ppRESO->GetYaxis()->SetTitleSize(0.06);
    hRes_pp_ppRESO->GetYaxis()->SetLabelSize(0.06);
    hRes_pp_ppRESO->GetYaxis()->CenterTitle();
    hRes_pp_ppRESO->GetYaxis()->SetTitleOffset(1.00);
    hRes_pp_ppRESO->GetYaxis()->SetRangeUser(0,kMax);
    //gStyle->SetPalette(51);
    TCanvas* cRes_pp_ppRESO = new TCanvas("cRes_pp_ppRESO", "cRes_pp_ppRESO", 1);
    cRes_pp_ppRESO->cd(0); cRes_pp_ppRESO->SetCanvasSize(1920/2, 1080/2); cRes_pp_ppRESO->SetMargin(0.15,0.05,0.2,0.05);//lrbt
    hRes_pp_ppRESO->Draw("COL");
    cRes_pp_ppRESO->SaveAs(OutputFolder+"cRes_pp_ppRESO.png");


    hRes_pp_pRESOpRESO->SetStats(false);
    hRes_pp_pRESOpRESO->SetTitle("");
    hRes_pp_pRESOpRESO->GetXaxis()->SetTitle("k* (Reso-Reso) (MeV)");
    hRes_pp_pRESOpRESO->GetXaxis()->SetTitleSize(0.06);
    hRes_pp_pRESOpRESO->GetXaxis()->SetLabelSize(0.06);
    hRes_pp_pRESOpRESO->GetXaxis()->CenterTitle();
    hRes_pp_pRESOpRESO->GetXaxis()->SetTitleOffset(1.3);
    hRes_pp_pRESOpRESO->GetXaxis()->SetLabelOffset(0.02);
    hRes_pp_pRESOpRESO->GetYaxis()->SetTitle("k* (p-p) (MeV)");
    hRes_pp_pRESOpRESO->GetYaxis()->SetTitleSize(0.06);
    hRes_pp_pRESOpRESO->GetYaxis()->SetLabelSize(0.06);
    hRes_pp_pRESOpRESO->GetYaxis()->CenterTitle();
    hRes_pp_pRESOpRESO->GetYaxis()->SetTitleOffset(1.00);
    hRes_pp_pRESOpRESO->GetYaxis()->SetRangeUser(0,kMax);
    //gStyle->SetPalette(51);
    TCanvas* cRes_pp_pRESOpRESO = new TCanvas("cRes_pp_pRESOpRESO", "cRes_pp_pRESOpRESO", 1);
    cRes_pp_pRESOpRESO->cd(0); cRes_pp_pRESOpRESO->SetCanvasSize(1920/2, 1080/2); cRes_pp_pRESOpRESO->SetMargin(0.15,0.05,0.2,0.05);//lrbt
    hRes_pp_pRESOpRESO->Draw("COL");
    cRes_pp_pRESOpRESO->SaveAs(OutputFolder+"cRes_pp_pRESOpRESO.png");

    hRes_pp_pL->Write();
    hRes_pL_LL->Write();
    hRes_pp_LL->Write();
    hRes_pL_pSigma0->Write();
    hRes_pL_pXim->Write();
    hRes_pXim_pXim1530->Write();
    hRes_pXim_pOmegam->Write();
    hRes_pp_ppRESO->Write();
    hRes_pp_pRESOpRESO->Write();

    //outputfile->Save();
    //outputfile->Close();

    //hist_pV0_pXim_relK_relK->Draw("colz");

    //delete t3_ran_p_proton;
    //delete t3_ran_p_Xim;
    //delete t3_ran_p_Lambda;
    //delete t3_ran_p_Sigma0;
    //delete t3_ran_p_Xim1530;
    delete hRes_pp_pL;
    delete hRes_pL_LL;
    delete hRes_pp_LL;
    delete hRes_pL_pSigma0;
    delete hRes_pL_pXim;
    delete hRes_pXim_pXim1530;
    delete hRes_pXim_pOmegam;
    delete hRes_pp_ppRESO;
    delete hRes_pp_pRESOpRESO;
    delete outputfile;

}

//used to generate the decay matrix
void DecayMatrixDimi(int SEED, const unsigned int NumBins, const double kMin, const double kMax, const unsigned int NumIterPerBin)
{
    // example of use of TGenPhaseSpace
    //Author: Valerio Filippini
    //further modified my O.W. Arnold and D.L. Mihaylov

    //if (!gROOT->GetClass("TGenPhaseSpace")) gSystem->Load("libPhysics");

    gRandom->SetSeed(SEED);

    TString output_string = TString::Format("../CorrelationFiles_2018/DecayMatrices/run2_DecayMatrixDimi_%i.root",SEED);
    TFile *outputfile = new TFile(output_string.Data(),"RECREATE");

    TH2F* hRes_pp_pL = new TH2F("hRes_pp_pL","x: relk pV0, y: relk pp",NumBins,kMin,kMax,NumBins,kMin,kMax);
    TH2F* hRes_pL_pSigma0 = new TH2F("hRes_pL_pSigma0","x: relk pSigma0, y: relk pV0",NumBins,kMin,kMax,NumBins,kMin,kMax);
    TH2F* hRes_pL_pXim = new TH2F("hRes_pL_pXim","x: relk pXim, y: relk pV0",NumBins,kMin,kMax,NumBins,kMin,kMax);
    TH2F* hRes_pXim_pXim1530 = new TH2F("hRes_pXim_pXim1530","x: relk pXim1530, y: relk pXim",NumBins,kMin,kMax,NumBins,kMin,kMax);
    TH2F* hRes_pXim_pOmegam = new TH2F("hRes_pXim_pOmegam","x: relk pOmegam, y: relk pXim",NumBins,kMin,kMax,NumBins,0,kMax);

    double Xim_decay_daughters[2] = {1116.,138.};//Lambda-pion
    double Lambda_decay_daughters[2] = {938.,138.};//proton-pion
    double Sigma0_decay_daughters[2] = {1116.,0.};//Lambda-photon
    double Xim1530_decay_daughters[2] = {1322.,135.};//Xim-pion0
    double Omegam_decay_daughters[2] = {1322.,135.};//Xim-pion0

    TLorentzVector tvec_proton,tvec_Xim,tvec_Lambda,tvec_Sigma0,tvec_Xim1530,tvec_Omegam;
    TLorentzVector *tvec_pp_pL,*tvec_pL_pSigma0,*tvec_pL_pXim,*tvec_pXim_pXim1530,*tvec_pXim_pOmegam;
    TGenPhaseSpace eventLambda,eventSigma0,eventXim,eventXim1530,eventOmegam;

    double relK_pX;
/*
    double px_proton;
    double py_proton;
    double pz_proton;
    double px_L;
    double py_L;
    double pz_L;
    double px_Sigma0;
    double py_Sigma0;
    double pz_Sigma0;
    double px_Xim;
    double py_Xim;
    double pz_Xim;
    double px_Xim1530;
    double py_Xim1530;
    double pz_Xim1530;
    double px_Omegam;
    double py_Omegam;
    double pz_Omegam;
*/
    for(unsigned int uBin=0; uBin<NumBins; uBin++){
        relK_pX = hRes_pp_pL->GetXaxis()->GetBinCenter(uBin+1);
        for(unsigned int uIter=0; uIter<NumIterPerBin; uIter++){
            //put momenta such that the proton and the X(whatever other particle)
            //are in their CM already here, and the relK corresponds to the uBin-th bin
            tvec_proton.SetXYZM(0,0,relK_pX,938);//units are in MeV
            tvec_Lambda.SetXYZM(0,0,-relK_pX,1116);
            tvec_Sigma0.SetXYZM(0,0,-relK_pX,1193);
            tvec_Xim.SetXYZM(0,0,-relK_pX,1322);
            tvec_Xim1530.SetXYZM(0,0,-relK_pX,1535);
            tvec_Omegam.SetXYZM(0,0,-relK_pX,1672);

//double relK_pL = relKcalc(tvec_proton,tvec_Lambda);
//double relK_pSigma0 = relKcalc(tvec_proton,tvec_Sigma0);
//double relK_pXim = relKcalc(tvec_proton,tvec_Xim);
//double relK_pXim1530 = relKcalc(tvec_proton,tvec_Xim1530);
//double relK_pOmegam = relKcalc(tvec_proton,tvec_Omegam);

//just for cross check
//printf("relK_pX=%.2f\n",relK_pX);
//printf("  relK_pL=%.2f\n", relK_pL);
//printf("  relK_pSigma0=%.2f\n", relK_pSigma0);
//printf("  relK_pXim=%.2f\n", relK_pXim);
//printf("  relK_pXim1530=%.2f\n", relK_pXim1530);
//printf("  relK_pOmegam=%.2f\n", relK_pOmegam);


            //p-V0 -> p-p
            eventLambda.SetDecay(tvec_Lambda, 2, Lambda_decay_daughters);
            eventLambda.Generate();
            tvec_pp_pL = eventLambda.GetDecay(0);
            double relK_pp_pL = relKcalc(tvec_proton,*tvec_pp_pL);
            hRes_pp_pL->Fill(relK_pX,relK_pp_pL);
//printf("    relK_pp_pL=%.2f\n",relK_pp_pL);
            //p-Sigma0 -> p-Lambda
            eventSigma0.SetDecay(tvec_Sigma0, 2, Sigma0_decay_daughters);
            eventSigma0.Generate();
            tvec_pL_pSigma0 = eventSigma0.GetDecay(0);
            double relK_pL_pSigma0 = relKcalc(tvec_proton,*tvec_pL_pSigma0);
            hRes_pL_pSigma0->Fill(relK_pX,relK_pL_pSigma0);

            //Simulate the decay and obtain the final momentum after the decay:
            //p-Xi -> p-V0
            eventXim.SetDecay(tvec_Xim, 2, Xim_decay_daughters);//Let the Xi decay in Lambda pion
            eventXim.Generate();
            tvec_pL_pXim = eventXim.GetDecay(0);
            double relK_pL_pXim = relKcalc(tvec_proton,*tvec_pL_pXim);
            hRes_pL_pXim->Fill(relK_pX,relK_pL_pXim);

            //p-Xim1530 -> p-Xim
            eventXim1530.SetDecay(tvec_Xim1530, 2, Xim1530_decay_daughters);
            eventXim1530.Generate();
            tvec_pXim_pXim1530 = eventXim1530.GetDecay(0);
            double relK_pXim_pXim1530 = relKcalc(tvec_proton,*tvec_pXim_pXim1530);
            hRes_pXim_pXim1530->Fill(relK_pX,relK_pXim_pXim1530);

            //p-Omegam -> p-Xim
            eventOmegam.SetDecay(tvec_Omegam, 2, Omegam_decay_daughters);
            eventOmegam.Generate();
            tvec_pXim_pOmegam = eventOmegam.GetDecay(0);
            double relK_pXim_pOmegam = relKcalc(tvec_proton,*tvec_pXim_pOmegam);
            hRes_pXim_pOmegam->Fill(relK_pX,relK_pXim_pOmegam);

        }
    }

    outputfile->cd();

    hRes_pp_pL->Write();
    hRes_pL_pSigma0->Write();
    hRes_pL_pXim->Write();
    hRes_pXim_pXim1530->Write();
    hRes_pXim_pOmegam->Write();

    //outputfile->Save();
    //outputfile->Close();

    //hist_pV0_pXim_relK_relK->Draw("colz");

    //delete t3_ran_p_proton;
    //delete t3_ran_p_Xim;
    //delete t3_ran_p_Lambda;
    //delete t3_ran_p_Sigma0;
    //delete t3_ran_p_Xim1530;
    delete hRes_pp_pL;
    delete hRes_pL_pSigma0;
    delete hRes_pL_pXim;
    delete hRes_pXim_pXim1530;
    delete hRes_pXim_pOmegam;
    delete outputfile;

}

void Resonances(int SEED, int NumIter)
{

    TRandom3 rangen(SEED);

    TString OutputFolder = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/GenBod/SourceWithReso/";
    TString output_string = TString::Format("%sResonances_%i.root",OutputFolder.Data(),SEED);
    TFile *outputfile = new TFile(output_string.Data(),"RECREATE");

    unsigned int NumBinsX = 250;
    unsigned int NumBinsY = 500;
    double kMaxX = 400;
    double kMaxY = 800;

    TH2F* hRes_pp_pR = new TH2F("hRes_pp_pR","hRes_pp_pR",NumBinsX,0,kMaxX,NumBinsY,0,kMaxY);
    TH2F* hRes_pp_RR = new TH2F("hRes_pp_RR","hRes_pp_RR",NumBinsX,0,kMaxX,NumBinsY,0,kMaxY);
    TH2F* hResBoost_pp_pR = new TH2F("hResBoost_pp_pR","hResBoost_pp_pR",NumBinsX,0,kMaxX,NumBinsY,0,kMaxY);
    TH2F* hResBoost_pp_RR = new TH2F("hResBoost_pp_RR","hResBoost_pp_RR",NumBinsX,0,kMaxX,NumBinsY,0,kMaxY);

    TH2F* hResBoost_pp_pR_R0 = new TH2F("hResBoost_pp_pR_R0","hResBoost_pp_pR_R0",NumBinsX,0,kMaxX,NumBinsY,0,kMaxY);
    TH2F* hResBoost_pp_RR_R0 = new TH2F("hResBoost_pp_RR_R0","hResBoost_pp_RR_R0",NumBinsX,0,kMaxX,NumBinsY,0,kMaxY);
    TH2F* hResBoost_pp_RR_R1 = new TH2F("hResBoost_pp_RR_R1","hResBoost_pp_RR_R1",NumBinsX,0,kMaxX,NumBinsY,0,kMaxY);

    //double R_decay_daughters[2] = {0.938,0.138};//pp feeddown
    //const double RESO_MASS = 1.36152;
    double R_decay_daughters[2] = {1.116,0.138};//pp feeddown
    const double RESO_MASS = 1.463;

    const double momMean = 0.0;
    const double momSpread = 0.4;

    TLorentzVector tvec_p,tvec_R0,tvec_R1;
    //relative momenta of the two (FINAL) protons
    TLorentzVector tvec_RelMom_pp_pR0;
    TLorentzVector tvec_RelMom_pp_R0R1;
    //total momenta of the two (FINAL) protons, to boost with
    TLorentzVector tvec_TotMom_pp_pR0;
    TLorentzVector tvec_TotMom_pp_R0R1;
    TLorentzVector *tvec_p_R0,*tvec_p_R1;

    TGenPhaseSpace eventR0,eventR1;

    for(int i=0;i<NumIter;i++){
        //if(i%1000==0) printf("i = %i(%i)\n",i,NumIter);

        double px_p = rangen.Gaus(momMean,momSpread);
        double py_p = rangen.Gaus(momMean,momSpread);
        double pz_p = rangen.Gaus(momMean,momSpread);

        double px_R0 = rangen.Gaus(momMean,momSpread);
        double py_R0 = rangen.Gaus(momMean,momSpread);
        double pz_R0 = rangen.Gaus(momMean,momSpread);

        double px_R1 = rangen.Gaus(momMean,momSpread);
        double py_R1 = rangen.Gaus(momMean,momSpread);
        double pz_R1 = rangen.Gaus(momMean,momSpread);

        //the momentum vectors of the primary particles
        tvec_p.SetXYZM(px_p,py_p,pz_p,R_decay_daughters[0]);//units are in GeV
        tvec_R0.SetXYZM(px_R0,py_R0,pz_R0,RESO_MASS);
        tvec_R1.SetXYZM(px_R1,py_R1,pz_R1,RESO_MASS);

        //Calculate initial momentum of the (prim-second and second-second) pair:
        double relK_pR0 = relKcalc(tvec_p,tvec_R0);
        double relK_R0R1 = relKcalc(tvec_R0,tvec_R1);

        eventR0.SetDecay(tvec_R0, 2, R_decay_daughters);
        eventR0.Generate();
        tvec_p_R0 = eventR0.GetDecay(0);

        eventR1.SetDecay(tvec_R1, 2, R_decay_daughters);
        eventR1.Generate();
        tvec_p_R1 = eventR1.GetDecay(0);

        tvec_RelMom_pp_pR0 = tvec_p-*tvec_p_R0;
/*
        if(tvec_RelMom_pp_pR0.P()<20e-3){

TLorentzVector p1=tvec_p;
TLorentzVector p2=*tvec_p_R0;
TLorentzVector r2=tvec_R0;
printf("tvec_RelMom_pp_pR0.P()=%.1f | tvec_R0.P()=%.1f\n",tvec_RelMom_pp_pR0.P()*1000,tvec_R0.P()*1000);
printf(" p1.P()=%.1f | p2.P()=%.1f | r2.P()=%.1f\n",p1.P()*1000,p2.P()*1000,r2.P()*1000);
printf(" BOOST into p1+p2:\n");
TLorentzVector BV = p1+p2;
p1.Boost(-BV.BoostVector());
p2.Boost(-BV.BoostVector());
r2.Boost(-BV.BoostVector());
printf("  p1.P()=%.1f | p2.P()=%.1f | r2.P()=%.1f\n",p1.P()*1000,p2.P()*1000,r2.P()*1000);
        }
*/
        tvec_RelMom_pp_R0R1 = *tvec_p_R0-*tvec_p_R1;
        tvec_TotMom_pp_pR0 = tvec_p+*tvec_p_R0;
        tvec_TotMom_pp_R0R1 = *tvec_p_R0+*tvec_p_R1;

        //relative momenta between the FINAL particles we study
        double relK_p_R0 = relKcalc(tvec_p,*tvec_p_R0);
        double relK_R0_R1 = relKcalc(*tvec_p_R0,*tvec_p_R1);

        hRes_pp_pR->Fill(relK_p_R0*1000,relK_pR0*1000);
        hRes_pp_RR->Fill(relK_R0_R1*1000,relK_R0R1*1000);

        //no effect

        tvec_p.Boost(-tvec_TotMom_pp_pR0.BoostVector());
        tvec_R0.Boost(-tvec_TotMom_pp_pR0.BoostVector());
        tvec_R1.Boost(-tvec_TotMom_pp_pR0.BoostVector());
        tvec_p_R0->Boost(-tvec_TotMom_pp_pR0.BoostVector());
        tvec_p_R1->Boost(-tvec_TotMom_pp_pR0.BoostVector());

        relK_pR0 = relKcalc(tvec_p,tvec_R0);
        relK_p_R0 = relKcalc(tvec_p,*tvec_p_R0);
        hResBoost_pp_pR->Fill(relK_p_R0*1000,relK_pR0*1000);
        hResBoost_pp_pR_R0->Fill(relK_p_R0*1000,tvec_R0.P()*1000);

        tvec_p.Boost(tvec_TotMom_pp_pR0.BoostVector());
        tvec_R0.Boost(tvec_TotMom_pp_pR0.BoostVector());
        tvec_R1.Boost(tvec_TotMom_pp_pR0.BoostVector());
        tvec_p_R0->Boost(tvec_TotMom_pp_pR0.BoostVector());
        tvec_p_R1->Boost(tvec_TotMom_pp_pR0.BoostVector());

        tvec_p.Boost(-tvec_TotMom_pp_R0R1.BoostVector());
        tvec_R0.Boost(-tvec_TotMom_pp_R0R1.BoostVector());
        tvec_R1.Boost(-tvec_TotMom_pp_R0R1.BoostVector());
        tvec_p_R0->Boost(-tvec_TotMom_pp_R0R1.BoostVector());
        tvec_p_R1->Boost(-tvec_TotMom_pp_R0R1.BoostVector());

        relK_R0R1 = relKcalc(tvec_R0,tvec_R1);
        relK_R0_R1 = relKcalc(*tvec_p_R0,*tvec_p_R1);

        hResBoost_pp_RR->Fill(relK_R0_R1*1000,relK_R0R1*1000);
        hResBoost_pp_RR_R0->Fill(relK_R0_R1*1000,tvec_R0.P()*1000);
        hResBoost_pp_RR_R1->Fill(relK_R0_R1*1000,tvec_R1.P()*1000);
    }

    outputfile->cd();

    //for(unsigned uBinX=1; uBinX<=hRes_pp_RESO->GetNbinsX(); uBinX++)for(unsigned uBinY=1; uBinY<=hRes_pp_RESO->GetNbinsY(); uBinY++)
    //    if(hRes_pp_RESO->GetBinContent(uBinX,uBinY))hRes_pp_RESO->SetBinContent(uBinX,uBinY,1);

    hRes_pp_pR->SetStats(false);
    hRes_pp_pR->SetTitle("");
    hRes_pp_pR->GetXaxis()->SetTitle("k* (p-p) (MeV)");
    hRes_pp_pR->GetXaxis()->SetTitleSize(0.06);
    hRes_pp_pR->GetXaxis()->SetLabelSize(0.06);
    hRes_pp_pR->GetXaxis()->CenterTitle();
    hRes_pp_pR->GetXaxis()->SetTitleOffset(1.3);
    hRes_pp_pR->GetXaxis()->SetLabelOffset(0.02);
    hRes_pp_pR->GetYaxis()->SetTitle("k* (p-Reso) (MeV)");
    hRes_pp_pR->GetYaxis()->SetTitleSize(0.06);
    hRes_pp_pR->GetYaxis()->SetLabelSize(0.06);
    hRes_pp_pR->GetYaxis()->CenterTitle();
    hRes_pp_pR->GetYaxis()->SetTitleOffset(1.00);
    hRes_pp_pR->GetYaxis()->SetRangeUser(0,kMaxY);
    //gStyle->SetPalette(51);
    TCanvas* cRes_pp_pR = new TCanvas("cRes_pp_pR", "cRes_pp_pR", 1);
    cRes_pp_pR->cd(0); cRes_pp_pR->SetCanvasSize(1920/2, 1080/2); cRes_pp_pR->SetMargin(0.15,0.05,0.2,0.05);//lrbt
    hRes_pp_pR->Draw("COL");
    cRes_pp_pR->SaveAs(OutputFolder+"cRes_pp_pR.png");

    hRes_pp_RR->SetStats(false);
    hRes_pp_RR->SetTitle("");
    hRes_pp_RR->GetXaxis()->SetTitle("k* (p-p) (MeV)");
    hRes_pp_RR->GetXaxis()->SetTitleSize(0.06);
    hRes_pp_RR->GetXaxis()->SetLabelSize(0.06);
    hRes_pp_RR->GetXaxis()->CenterTitle();
    hRes_pp_RR->GetXaxis()->SetTitleOffset(1.3);
    hRes_pp_RR->GetXaxis()->SetLabelOffset(0.02);
    hRes_pp_RR->GetYaxis()->SetTitle("k* (Reso-Reso) (MeV)");
    hRes_pp_RR->GetYaxis()->SetTitleSize(0.06);
    hRes_pp_RR->GetYaxis()->SetLabelSize(0.06);
    hRes_pp_RR->GetYaxis()->CenterTitle();
    hRes_pp_RR->GetYaxis()->SetTitleOffset(1.00);
    hRes_pp_RR->GetYaxis()->SetRangeUser(0,kMaxY);
    //gStyle->SetPalette(51);
    TCanvas* cRes_pp_RR = new TCanvas("cRes_pp_RR", "cRes_pp_RR", 1);
    cRes_pp_RR->cd(0); cRes_pp_RR->SetCanvasSize(1920/2, 1080/2); cRes_pp_RR->SetMargin(0.15,0.05,0.2,0.05);//lrbt
    hRes_pp_RR->Draw("COL");
    cRes_pp_RR->SaveAs(OutputFolder+"cRes_pp_RR.png");

    hResBoost_pp_RR_R0->SetStats(false);
    hResBoost_pp_RR_R0->SetTitle("");
    hResBoost_pp_RR_R0->GetXaxis()->SetTitle("k* (MeV)");
    hResBoost_pp_RR_R0->GetXaxis()->SetTitleSize(0.06);
    hResBoost_pp_RR_R0->GetXaxis()->SetLabelSize(0.06);
    hResBoost_pp_RR_R0->GetXaxis()->CenterTitle();
    hResBoost_pp_RR_R0->GetXaxis()->SetTitleOffset(1.3);
    hResBoost_pp_RR_R0->GetXaxis()->SetLabelOffset(0.02);
    hResBoost_pp_RR_R0->GetYaxis()->SetTitle("p*_{res} (MeV)");
    hResBoost_pp_RR_R0->GetYaxis()->SetTitleSize(0.06);
    hResBoost_pp_RR_R0->GetYaxis()->SetLabelSize(0.06);
    hResBoost_pp_RR_R0->GetYaxis()->CenterTitle();
    hResBoost_pp_RR_R0->GetYaxis()->SetTitleOffset(1.00);
    hResBoost_pp_RR_R0->GetYaxis()->SetRangeUser(0,kMaxY);
    //gStyle->SetPalette(51);
    TCanvas* cResBoost_pp_RR_R0 = new TCanvas("cResBoost_pp_RR_R0", "cResBoost_pp_RR_R0", 1);
    cResBoost_pp_RR_R0->cd(0); cResBoost_pp_RR_R0->SetCanvasSize(1920/2, 1080/2); cResBoost_pp_RR_R0->SetMargin(0.15,0.05,0.2,0.05);//lrbt
    hResBoost_pp_RR_R0->Draw("COL");
    cResBoost_pp_RR_R0->SaveAs(OutputFolder+"cResBoost_pp_RR_R0.png");

    hRes_pp_pR->Write();
    hRes_pp_RR->Write();
    hResBoost_pp_pR->Write();
    hResBoost_pp_RR->Write();
    hResBoost_pp_pR_R0->Write();
    hResBoost_pp_RR_R0->Write();
    hResBoost_pp_RR_R1->Write();

    //outputfile->Save();
    //outputfile->Close();

    //hist_pV0_pXim_relK_relK->Draw("colz");

    //delete t3_ran_p_proton;
    //delete t3_ran_p_Xim;
    //delete t3_ran_p_Lambda;
    //delete t3_ran_p_Sigma0;
    //delete t3_ran_p_Xim1530;
    delete hRes_pp_pR;
    delete hRes_pp_RR;
    delete hResBoost_pp_pR;
    delete hResBoost_pp_RR;
    delete hResBoost_pp_pR_R0;
    delete hResBoost_pp_RR_R0;
    delete hResBoost_pp_RR_R1;
    delete cRes_pp_pR;
    delete cRes_pp_RR;
    delete outputfile;
}




int GENBOD(int narg, char** ARGS){
    //DecayMatrix(atoi(ARGS[1]),atoi(ARGS[2]));
    Resonances(atoi(ARGS[1]),atoi(ARGS[2]));
}
