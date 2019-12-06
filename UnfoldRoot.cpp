
#include <iostream>
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TRandom3.h"
#include "TString.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TLine.h"
#include "TSVDUnfold.h"
#include "TFile.h"

using namespace std;

Double_t Reconstruct( Double_t xt, TRandom3& R )
{
   // apply some Gaussian smearing + bias and efficiency corrections to fake reconstruction
   const Double_t cutdummy = -99999.0;
   Double_t xeff = 0.3 + (1.0 - 0.3)/20.0*(xt + 10.0);  // efficiency
   Double_t x    = R.Rndm();
   if (x > xeff) return cutdummy;
   else {
     Double_t xsmear= R.Gaus(-2.5,0.2); // bias and smear
     return xt+xsmear;
   }
}
void TSVDUnfoldExample()
{
   gROOT->SetStyle("Plain");
   gStyle->SetOptStat(0);
   TRandom3 R;
   const Double_t cutdummy= -99999.0;
   // --------------------------------------
   // Data/MC toy generation
   //
   // The MC input
   Int_t nbins = 40;
   TH1D *xini = new TH1D("xini", "MC truth", nbins, -10.0, 10.0);
   TH1D *bini = new TH1D("bini", "MC reco", nbins, -10.0, 10.0);
   TH2D *Adet = new TH2D("Adet", "detector response", nbins, -10.0, 10.0, nbins, -10.0, 10.0);
   // Data
   TH1D *data = new TH1D("data", "data", nbins, -10.0, 10.0);
   // Data "truth" distribution to test the unfolding
   TH1D *datatrue = new TH1D("datatrue", "data truth", nbins, -10.0, 10.0);
   // Statistical covariance matrix
   TH2D *statcov = new TH2D("statcov", "covariance matrix", nbins, -10.0, 10.0, nbins, -10.0, 10.0);
   // Fill the MC using a Breit-Wigner, mean 0.3 and width 2.5.
   for (Int_t i= 0; i<100000; i++) {
      Double_t xt = R.BreitWigner(0.3, 2.5);
      xini->Fill(xt);
      Double_t x = Reconstruct( xt, R );
      if (x != cutdummy) {
         Adet->Fill(x, xt);
         bini->Fill(x);
      }
   }
   // Fill the "data" with a Gaussian, mean 0 and width 2.
   for (Int_t i=0; i<10000; i++) {
      Double_t xt = R.Gaus(0.0, 2.0);
      datatrue->Fill(xt);
      Double_t x = Reconstruct( xt, R );
      if (x != cutdummy)
      data->Fill(x);
   }
   cout << "Created toy distributions and errors for: " << endl;
   cout << "... \"true MC\"   and \"reconstructed (smeared) MC\"" << endl;
   cout << "... \"true data\" and \"reconstructed (smeared) data\"" << endl;
   cout << "... the \"detector response matrix\"" << endl;
   // Fill the data covariance matrix
   for (int i=1; i<=data->GetNbinsX(); i++) {
       statcov->SetBinContent(i,i,data->GetBinError(i)*data->GetBinError(i));
   }
   // ----------------------------
   // Here starts the actual unfolding
   //
   // Create TSVDUnfold object and initialise
   TSVDUnfold *tsvdunf = new TSVDUnfold( data, statcov, bini, xini, Adet );
   // It is possible to normalise unfolded spectrum to unit area
   tsvdunf->SetNormalize( kFALSE ); // no normalisation here
   // Perform the unfolding with regularisation parameter kreg = 13
   // - the larger kreg, the finer grained the unfolding, but the more fluctuations occur
   // - the smaller kreg, the stronger is the regularisation and the bias
   TH1D* unfres = tsvdunf->Unfold( 13 );
   // Get the distribution of the d to cross check the regularization
   // - choose kreg to be the point where |d_i| stop being statistically significantly >>1
   TH1D* ddist = tsvdunf->GetD();
   // Get the distribution of the singular values
   TH1D* svdist = tsvdunf->GetSV();
   // Compute the error matrix for the unfolded spectrum using toy MC
   // using the measured covariance matrix as input to generate the toys
   // 100 toys should usually be enough
   // The same method can be used for different covariance matrices separately.
   TH2D* ustatcov = tsvdunf->GetUnfoldCovMatrix( statcov, 100 );
   // Now compute the error matrix on the unfolded distribution originating
   // from the finite detector matrix statistics
   TH2D* uadetcov = tsvdunf->GetAdetCovMatrix( 100 );
   // Sum up the two (they are uncorrelated)
   ustatcov->Add( uadetcov );
   //Get the computed regularized covariance matrix (always corresponding to total uncertainty passed in constructor) and add uncertainties from finite MC statistics.
   TH2D* utaucov = tsvdunf->GetXtau();
   utaucov->Add( uadetcov );
   //Get the computed inverse of the covariance matrix
   TH2D* uinvcov = tsvdunf->GetXinv();
   // ---------------------------------
   // Only plotting stuff below

   for (int i=1; i<=unfres->GetNbinsX(); i++) {
      unfres->SetBinError(i, TMath::Sqrt(utaucov->GetBinContent(i,i)));
   }
   // Renormalize just to be able to plot on the same scale
   //xini->Scale(0.7*datatrue->Integral()/xini->Integral());
   //bini->Scale(0.7*datatrue->Integral()/bini->Integral());
   TLegend *leg = new TLegend(0.58,0.60,0.99,0.88);
   leg->SetBorderSize(0);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   leg->AddEntry(unfres,"Unfolded Data","p");
   leg->AddEntry(datatrue,"True Data","l");
   leg->AddEntry(data,"Reconstructed Data","l");
   leg->AddEntry(xini,"True MC","l");
   TCanvas *c1 = new TCanvas( "c1", "Unfolding toy example with TSVDUnfold", 1000, 900 );
   c1->Divide(1,2);
   TVirtualPad * c11 = c1->cd(1);
   TH1D* frame = new TH1D( *unfres );
   frame->SetTitle( "Unfolding toy example with TSVDUnfold" );
   frame->GetXaxis()->SetTitle( "x variable" );
   frame->GetYaxis()->SetTitle( "Events" );
   frame->GetXaxis()->SetTitleOffset( 1.25 );
   frame->GetYaxis()->SetTitleOffset( 1.29 );
   frame->Draw();
   data->SetLineStyle(2);
   data->SetLineColor(4);
   data->SetLineWidth(2);
   unfres->SetMarkerStyle(20);
   datatrue->SetLineColor(2);
   datatrue->SetLineWidth(2);
   xini->SetLineStyle(2);
   xini->SetLineColor(8);
   xini->SetLineWidth(2);
   // ------------------------------------------------------------
   // add histograms
   unfres->Draw("same");
   datatrue->Draw("same");
   data->Draw("same");
   xini->Draw("same");
   leg->Draw();
   // covariance matrix
   TVirtualPad * c12 = c1->cd(2);
   c12->Divide(2,1);
   TVirtualPad * c2 = c12->cd(1);
   c2->SetRightMargin   ( 0.15         );
   TH2D* covframe = new TH2D( *ustatcov );
   covframe->SetTitle( "TSVDUnfold covariance matrix" );
   covframe->GetXaxis()->SetTitle( "x variable" );
   covframe->GetYaxis()->SetTitle( "x variable" );
   covframe->GetXaxis()->SetTitleOffset( 1.25 );
   covframe->GetYaxis()->SetTitleOffset( 1.29 );
   covframe->Draw();
   ustatcov->SetLineWidth( 2 );
   ustatcov->Draw( "colzsame" );
   // distribution of the d quantity
   TVirtualPad * c3 = c12->cd(2);
   c3->SetLogy();
   TLine *line = new TLine( 0.,1.,40.,1. );
   line->SetLineStyle(2);
   TH1D* dframe = new TH1D( *ddist );
   dframe->SetTitle( "TSVDUnfold |d_{i}|" );
   dframe->GetXaxis()->SetTitle( "i" );
   dframe->GetYaxis()->SetTitle( "|d_{i}|" );
   dframe->GetXaxis()->SetTitleOffset( 1.25 );
   dframe->GetYaxis()->SetTitleOffset( 1.29 );
   dframe->SetMinimum( 0.001 );
   dframe->Draw();
   ddist->SetLineWidth( 2 );
   ddist->Draw( "same" );
   line->Draw();

   TFile* fOutput = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/Unfold/Example1.root","recreate");

    c1->Write();
    data->Write();
    statcov->Write();
    bini->Write();
    xini->Write();
    Adet->Write();

    ustatcov->Write();
    utaucov->Write();
    datatrue->Write();
    unfres->Write();

    printf("Integral_McReco = %e\n",bini->Integral());
    printf("Integral_McTrue = %e\n",xini->Integral());

    printf("Integral_Reco = %e\n",data->Integral());
    printf("Integral_True = %e\n",datatrue->Integral());
    printf("Integral_Unfo = %e\n",unfres->Integral());



   delete fOutput;
}

void TSVDUnfoldExampleDimi()
{
   gROOT->SetStyle("Plain");
   gStyle->SetOptStat(0);
   TRandom3 R;
   const Double_t cutdummy= -99999.0;
   const double xmin = 0;
   const double xmax = 10;
   const double RelSmear= 0.1;
   const double AbsSmear= 0.1;
   const double AbsOffset = 0.1;


TFile* fOutput = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/Unfold/ExampleDimi.root","recreate");

   // --------------------------------------
   // Data/MC toy generation
   //
   // The MC input
   Int_t nbins = 40;
   TH1D *GOD = new TH1D("GOD", "GOD", nbins, xmin, xmax);
   TH1D *mctrue = new TH1D("mctrue", "MC truth", nbins, xmin, xmax);
   TH1D *mcreco = new TH1D("mcreco", "MC reco", nbins, xmin, xmax);
   TH2D *Adet = new TH2D("Adet", "detector response", nbins, xmin, xmax, nbins, xmin, xmax);
   // Data
   TH1D *data = new TH1D("data", "data", nbins, xmin, xmax);
   // Data "truth" distribution to test the unfolding
   TH1D *datatrue = new TH1D("datatrue", "data truth", nbins, xmin, xmax);
   // Statistical covariance matrix
   TH2D *statcov = new TH2D("statcov", "covariance matrix", nbins, xmin, xmax, nbins, xmin, xmax);


    // Fill the MC using a Breit-Wigner, mean 0.3 and width 2.5.
    for(int iBin=1; iBin<=nbins; iBin++){
        //double xtrue = TMath::Gaus(mctrue->GetBinCenter(iBin),0.,100.)*(1.+0.0001*mctrue->GetBinCenter(iBin));
        //double xtrue = TMath::BreitWigner(GOD->GetBinCenter(iBin),0.3,2.5);
        double xValue = mctrue->GetBinCenter(iBin);
        double yValue = exp(-pow(xValue/3.,2))+1;
        GOD->SetBinContent(iBin,yValue+1);
        GOD->SetBinError(iBin,0);
    }


    const unsigned NumEntries = 100000;
    for(unsigned uIter=0; uIter<NumEntries; uIter++){
        double RandVal = GOD->GetRandom();
        double SmearedVal;
        SmearedVal = R.Gaus(RandVal+AbsOffset,RandVal*RelSmear+AbsSmear);
        //SmearedVal = Reconstruct( RandVal, R );
        mctrue->Fill(RandVal);
        mcreco->Fill(SmearedVal);
        for(unsigned idet=0; idet<100; idet++){
            SmearedVal = R.Gaus(RandVal+AbsOffset,RandVal*RelSmear+AbsSmear);
            //SmearedVal = Reconstruct( RandVal, R );
            Adet->Fill(SmearedVal, RandVal);
        }
        SmearedVal = R.Gaus(RandVal+AbsOffset,RandVal*RelSmear+AbsSmear);
        //SmearedVal = Reconstruct( RandVal, R );
        datatrue->Fill(RandVal);
        data->Fill(SmearedVal);
    }
    mctrue->Sumw2();
    mcreco->Sumw2();
    data->Sumw2();

    GOD->Scale(1./GOD->Integral());
    datatrue->Scale(1./mctrue->Integral());
    mctrue->Scale(1./mctrue->Integral());
    mcreco->Scale(1./mcreco->Integral());
    data->Scale(1./data->Integral());

/*
   for (Int_t i= 0; i<100000; i++) {
      Double_t xt = R.BreitWigner(0.3, 2.5);
      mctrue->Fill(xt);
      Double_t x = Reconstruct( xt, R );
      if (x != cutdummy) {
         Adet->Fill(x, xt);
         mcreco->Fill(x);
      }
   }
   // Fill the "data" with a Gaussian, mean 0 and width 2.
   for (Int_t i=0; i<10000; i++) {
      Double_t xt = R.Gaus(0.0, 2.0);
      datatrue->Fill(xt);
      Double_t x = Reconstruct( xt, R );
      if (x != cutdummy)
      data->Fill(x);
   }
*/


    //mctrue->Scale(double(NumEntries)/mctrue->Integral());
    //datatrue->Scale(double(NumEntries)/datatrue->Integral());

    for(int iBin=1; iBin<=nbins; iBin++){
        //data->SetBinError(iBin,sqrt(data->GetBinContent(iBin)));
        printf("x=%f; y=%f; e=%f\n",data->GetBinCenter(iBin),data->GetBinContent(iBin),data->GetBinError(iBin));
    }
usleep(1000e3);

    data->Write();
    GOD->Write();
    mctrue->Write();
    Adet->Write();

    cout << "Created toy distributions and errors for: " << endl;
    cout << "... \"true MC\"   and \"reconstructed (smeared) MC\"" << endl;
    cout << "... \"true data\" and \"reconstructed (smeared) data\"" << endl;
    cout << "... the \"detector response matrix\"" << endl;
   // Fill the data covariance matrix
    for (int i=1; i<=data->GetNbinsX(); i++) {
        statcov->SetBinContent(i,i,data->GetBinError(i)*data->GetBinError(i));
    }
    // ----------------------------
    // Here starts the actual unfolding
    //
    // Create TSVDUnfold object and initialise
    TSVDUnfold *tsvdunf = new TSVDUnfold( data, statcov, mcreco, mctrue, Adet );
    // It is possible to normalise unfolded spectrum to unit area
    tsvdunf->SetNormalize( false ); // no normalisation here
    // Perform the unfolding with regularisation parameter kreg = 13
    // - the larger kreg, the finer grained the unfolding, but the more fluctuations occur
    // - the smaller kreg, the stronger is the regularisation and the bias
    TH1D* unfres = tsvdunf->Unfold( 13 );
    // Get the distribution of the d to cross check the regularization
    // - choose kreg to be the point where |d_i| stop being statistically significantly >>1
    TH1D* ddist = tsvdunf->GetD();
    // Get the distribution of the singular values
    TH1D* svdist = tsvdunf->GetSV();
    // Compute the error matrix for the unfolded spectrum using toy MC
    // using the measured covariance matrix as input to generate the toys
    // 100 toys should usually be enough
    // The same method can be used for different covariance matrices separately.
    TH2D* ustatcov = tsvdunf->GetUnfoldCovMatrix( statcov, 100 );
    // Now compute the error matrix on the unfolded distribution originating
    // from the finite detector matrix statistics
    TH2D* uadetcov = tsvdunf->GetAdetCovMatrix( 100 );
    // Sum up the two (they are uncorrelated)
    ustatcov->Add( uadetcov );
    //Get the computed regularized covariance matrix (always corresponding to total uncertainty passed in constructor) and add uncertainties from finite MC statistics.
    TH2D* utaucov = tsvdunf->GetXtau();
    utaucov->Add( uadetcov );
    //Get the computed inverse of the covariance matrix
    TH2D* uinvcov = tsvdunf->GetXinv();
    // ---------------------------------
    // Only plotting stuff below

   for (int i=1; i<=unfres->GetNbinsX(); i++) {
      unfres->SetBinError(i, TMath::Sqrt(utaucov->GetBinContent(i,i)));
   }
   // Renormalize just to be able to plot on the same scale
   //mctrue->Scale(0.7*datatrue->Integral()/mctrue->Integral());
   //mcreco->Scale(0.7*datatrue->Integral()/mcreco->Integral());
   unfres->Scale(1./unfres->Integral());

   TLegend *leg = new TLegend(0.58,0.60,0.99,0.88);
   leg->SetBorderSize(0);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   leg->AddEntry(unfres,"Unfolded Data","p");
   leg->AddEntry(datatrue,"True Data","l");
   leg->AddEntry(data,"Reconstructed Data","l");
   leg->AddEntry(mctrue,"True MC","l");
   TCanvas *c1 = new TCanvas( "c1", "Unfolding toy example with TSVDUnfold", 1000, 900 );
   c1->Divide(1,2);
   TVirtualPad * c11 = c1->cd(1);
   TH1D* frame = new TH1D( *unfres );
   frame->SetTitle( "Unfolding toy example with TSVDUnfold" );
   frame->GetXaxis()->SetTitle( "x variable" );
   frame->GetYaxis()->SetTitle( "Events" );
   frame->GetXaxis()->SetTitleOffset( 1.25 );
   frame->GetYaxis()->SetTitleOffset( 1.29 );
   frame->Draw();
   data->SetLineStyle(2);
   data->SetLineColor(4);
   data->SetLineWidth(2);
   unfres->SetMarkerStyle(20);
   datatrue->SetLineColor(2);
   datatrue->SetLineWidth(2);
   mctrue->SetLineStyle(2);
   mctrue->SetLineColor(8);
   mctrue->SetLineWidth(2);
   // ------------------------------------------------------------
   // add histograms
   unfres->Draw("same");
   datatrue->Draw("same");
   data->Draw("same");
   mctrue->Draw("same");
   leg->Draw();
   // covariance matrix
   TVirtualPad * c12 = c1->cd(2);
   c12->Divide(2,1);
   TVirtualPad * c2 = c12->cd(1);
   c2->SetRightMargin   ( 0.15         );
   TH2D* covframe = new TH2D( *ustatcov );
   covframe->SetTitle( "TSVDUnfold covariance matrix" );
   covframe->GetXaxis()->SetTitle( "x variable" );
   covframe->GetYaxis()->SetTitle( "x variable" );
   covframe->GetXaxis()->SetTitleOffset( 1.25 );
   covframe->GetYaxis()->SetTitleOffset( 1.29 );
   covframe->Draw();
   ustatcov->SetLineWidth( 2 );
   ustatcov->Draw( "colzsame" );
   // distribution of the d quantity
   TVirtualPad * c3 = c12->cd(2);
   c3->SetLogy();
   TLine *line = new TLine( 0.,1.,40.,1. );
   line->SetLineStyle(2);
   TH1D* dframe = new TH1D( *ddist );
   dframe->SetTitle( "TSVDUnfold |d_{i}|" );
   dframe->GetXaxis()->SetTitle( "i" );
   dframe->GetYaxis()->SetTitle( "|d_{i}|" );
   dframe->GetXaxis()->SetTitleOffset( 1.25 );
   dframe->GetYaxis()->SetTitleOffset( 1.29 );
   dframe->SetMinimum( 0.001 );
   dframe->Draw();
   ddist->SetLineWidth( 2 );
   ddist->Draw( "same" );
   line->Draw();



    c1->Write();
    //data->Write();
    statcov->Write();
    mcreco->Write();
    //mctrue->Write();
    //Adet->Write();

    ustatcov->Write();
    utaucov->Write();
    datatrue->Write();
    unfres->Write();

    printf("Integral_McReco = %e\n",mcreco->Integral());
    printf("Integral_McTrue = %e\n",mctrue->Integral());

    printf("Integral_data = %e\n",data->Integral());
    printf("Integral_datatrue = %e\n",datatrue->Integral());
    printf("Integral_unfres = %e\n",unfres->Integral());



   delete fOutput;
}

int UNFOLD_MAIN(int argc, char *argv[]){
    //TSVDUnfoldExample();
    TSVDUnfoldExampleDimi();
}
